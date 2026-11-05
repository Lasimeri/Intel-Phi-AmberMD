/* amview.c: 3D view of molecules on the host's GPU (OpenGL 4.6 core), 1920 x 1080 at 60 fps.
 *
 * Two sources:
 *   card mode (default)  the live MD run on a Xeon Phi: the frame ring (src/card/mdframe.h) and
 *                        the stats block (src/card/mdstat.h) in the shared host-memory region
 *   file mode (--file)   structures and trajectories on disk through molio.c: PDB, mmCIF, SDF,
 *                        GROMACS .gro, mdCATH HDF5 (every temperature and replica)
 * Every atom is a ray-cast sphere impostor (exact surface and depth per pixel). In file mode the
 * backbone is a tube of spheres through the C-alpha atoms, drug molecules (HETATM residues of 6 or
 * more atoms) are drawn larger with green carbons, and playback interpolates between frames.
 *
 * Usage: amview [--card N] [--seconds S] [--radius A] [--snap FILE]
 *        amview --file F [F ...] [--temp K] [--replica R] [--fps F] [--seconds S] [--snap FILE]
 *        amview --file SOLUTE.pdb --follow TRAJ.mdcrd --natoms N   (a sander run as it writes)
 * Keys: mouse drag rotate, wheel zoom, o auto-rotation, q/Esc quit;
 *   file mode: space pause, [ ] speed, , . step, Home restart, c colour (structure / element /
 *   chain rainbow), h hydrogens, a atoms, b backbone tube, n p next/previous file,
 *   t y next temperature / replica (mdCATH).
 * Every second it prints frames drawn, loop interval p50 and p99, GPU time per frame.
 */
#define _GNU_SOURCE
#include <epoxy/gl.h>
#include <GLFW/glfw3.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <time.h>
#include <unistd.h>
#include "../card/mdframe.h"
#include "../card/mdstat.h"
#include "molio.h"

#define WIDTH 1920
#define HEIGHT 1080
#define PERIOD_NS 16666667ll
#define NQ 4
#define TUBE_STEPS 8
#define STICK_RADIUS 0.17f   /* Angstrom */

typedef struct { float m[16]; } mat4;
typedef struct { float x, y, z, r, cr, cg, cb, ca; } inst_t; /* one sphere: centre, radius, colour */

static mat4 perspective(float fovy, float aspect, float zn, float zf) {
    mat4 r = {{0}};
    float f = 1.0f / tanf(fovy / 2);
    r.m[0] = f / aspect;
    r.m[5] = f;
    r.m[10] = (zf + zn) / (zn - zf);
    r.m[11] = -1;
    r.m[14] = 2 * zf * zn / (zn - zf);
    return r;
}

/* camera on a sphere around the origin, looking at it */
static mat4 orbit(float yaw, float pitch, float dist) {
    float cy = cosf(yaw), sy = sinf(yaw), cp = cosf(pitch), sp = sinf(pitch);
    float ex = dist * cp * sy, ey = dist * sp, ez = dist * cp * cy;
    float fx = -ex / dist, fy = -ey / dist, fz = -ez / dist;
    float rx = fz, ry = 0, rz = -fx;
    float rl = sqrtf(rx * rx + rz * rz);
    rx /= rl; rz /= rl;
    float ux = ry * fz - rz * fy, uy = rz * fx - rx * fz, uz = rx * fy - ry * fx;
    mat4 v = {{rx, ux, -fx, 0, ry, uy, -fy, 0, rz, uz, -fz, 0,
               -(rx * ex + ry * ey + rz * ez), -(ux * ex + uy * ey + uz * ez), (fx * ex + fy * ey + fz * ez), 1}};
    return v;
}

static const char *atom_vs =
    "#version 460 core\n"
    "struct Inst { vec4 pr; vec4 c; };\n"
    "layout(std430, binding = 0) readonly buffer Insts { Inst inst[]; };\n"
    "uniform mat4 view; uniform mat4 proj;\n"
    "out vec3 v_center; out vec3 v_pos; out vec3 v_col; out float v_r;\n"
    "const vec2 corner[6] = vec2[](vec2(-1,-1), vec2(1,-1), vec2(1,1), vec2(-1,-1), vec2(1,1), vec2(-1,1));\n"
    "void main() {\n"
    "  Inst a = inst[gl_InstanceID];\n"
    "  vec3 c = (view * vec4(a.pr.xyz, 1.0)).xyz;\n"
    "  vec3 p = c + vec3(corner[gl_VertexID] * a.pr.w * 1.5, a.pr.w);\n"
    "  v_center = c; v_pos = p; v_col = a.c.rgb; v_r = a.pr.w;\n"
    "  gl_Position = proj * vec4(p, 1.0);\n"
    "}\n";

static const char *atom_fs =
    "#version 460 core\n"
    "layout(depth_greater) out float gl_FragDepth;\n"
    "in vec3 v_center; in vec3 v_pos; in vec3 v_col; in float v_r;\n"
    "uniform mat4 proj; uniform vec2 cue; /* depth cueing: view distance where fading starts, where it is full */\n"
    "out vec4 color;\n"
    "void main() {\n"
    "  vec3 d = normalize(v_pos);\n"
    "  float b = dot(d, v_center);\n"
    "  float h = b*b - dot(v_center, v_center) + v_r*v_r;\n"
    "  if (h < 0.0) discard;\n"
    "  vec3 hit = (b - sqrt(h)) * d;\n"
    "  vec3 n = normalize(hit - v_center);\n"
    "  vec3 l = normalize(vec3(0.4, 0.6, 0.7));\n"
    "  float diff = max(dot(n, l), 0.0);\n"
    "  float spec = pow(max(dot(n, normalize(l - d)), 0.0), 48.0);\n"
    "  float rim = pow(1.0 - max(dot(n, -d), 0.0), 3.0);\n"
    "  float fade = clamp((-hit.z - cue.x) / max(cue.y - cue.x, 1e-3), 0.0, 1.0);\n"
    "  vec3 lit = v_col * (0.22 + 0.78 * diff) + vec3(0.35) * spec + v_col * 0.25 * rim;\n"
    "  color = vec4(mix(lit, lit * 0.45, fade), 1.0);\n"
    "  vec4 clip = proj * vec4(hit, 1.0);\n"
    "  gl_FragDepth = clip.z / clip.w * 0.5 + 0.5;\n"
    "}\n";

static const char *box_vs =
    "#version 460 core\n"
    "uniform mat4 view; uniform mat4 proj; uniform vec3 half_box;\n"
    "const ivec2 e[12] = ivec2[](ivec2(0,1), ivec2(2,3), ivec2(4,5), ivec2(6,7), ivec2(0,2), ivec2(1,3),\n"
    "  ivec2(4,6), ivec2(5,7), ivec2(0,4), ivec2(1,5), ivec2(2,6), ivec2(3,7));\n"
    "void main() {\n"
    "  int v = gl_VertexID % 2 == 0 ? e[gl_VertexID / 2].x : e[gl_VertexID / 2].y;\n"
    "  vec3 p = vec3((v & 1) != 0 ? 1.0 : -1.0, (v & 2) != 0 ? 1.0 : -1.0, (v & 4) != 0 ? 1.0 : -1.0) * half_box;\n"
    "  gl_Position = proj * view * vec4(p, 1.0);\n"
    "}\n";

static const char *box_fs =
    "#version 460 core\n"
    "out vec4 color;\n"
    "void main() { color = vec4(0.55, 0.6, 0.7, 1.0); }\n";

/* cylinders: bonds. Endpoints a and b, radius, two colours split at the middle, dashed flag.
 * A view-space quad in front of the cylinder covers it; the fragment shader intersects the eye ray
 * with the cylinder exactly and writes the true depth, like the sphere impostors. */
typedef struct { float ax, ay, az, r, bx, by, bz, dash, car, cag, cab, pad0, cbr, cbg, cbb, pad1; } cyl_t;

static const char *cyl_vs =
    "#version 460 core\n"
    "struct Cyl { vec4 ar; vec4 bd; vec4 ca; vec4 cb; };\n"
    "layout(std430, binding = 1) readonly buffer Cyls { Cyl cyl[]; };\n"
    "uniform mat4 view; uniform mat4 proj;\n"
    "out vec3 v_a; out vec3 v_b; out vec3 v_pos; out float v_r; out float v_dash; out vec3 v_ca; out vec3 v_cb;\n"
    "const vec2 corner[6] = vec2[](vec2(0,-1), vec2(1,-1), vec2(1,1), vec2(0,-1), vec2(1,1), vec2(0,1));\n"
    "void main() {\n"
    "  Cyl c = cyl[gl_InstanceID];\n"
    "  vec3 a = (view * vec4(c.ar.xyz, 1.0)).xyz, b = (view * vec4(c.bd.xyz, 1.0)).xyz;\n"
    "  float r = c.ar.w;\n"
    "  vec2 e = b.xy - a.xy; float el = length(e);\n"
    "  e = el > 1e-5 ? e / el : vec2(1.0, 0.0);\n"
    "  vec2 n = vec2(-e.y, e.x);\n"
    "  vec2 k = corner[gl_VertexID];\n"
    "  vec2 xy = mix(a.xy - e * r, b.xy + e * r, k.x) + n * (k.y * r * 1.2);\n"
    "  vec3 p = vec3(xy, max(a.z, b.z) + r);\n"
    "  v_a = a; v_b = b; v_pos = p; v_r = r; v_dash = c.bd.w; v_ca = c.ca.rgb; v_cb = c.cb.rgb;\n"
    "  gl_Position = proj * vec4(p, 1.0);\n"
    "}\n";

static const char *cyl_fs =
    "#version 460 core\n"
    "layout(depth_greater) out float gl_FragDepth;\n"
    "in vec3 v_a; in vec3 v_b; in vec3 v_pos; in float v_r; in float v_dash; in vec3 v_ca; in vec3 v_cb;\n"
    "uniform mat4 proj; uniform vec2 cue; /* depth cueing: view distance where fading starts, where it is full */\n"
    "out vec4 color;\n"
    "void main() {\n"
    "  vec3 d = normalize(v_pos);\n"
    "  vec3 w = v_b - v_a; float L = length(w); w /= L;\n"
    "  vec3 X = d - dot(d, w) * w;\n"
    "  vec3 Y = -v_a - dot(-v_a, w) * w;\n"
    "  float A = dot(X, X), B = 2.0 * dot(X, Y), C = dot(Y, Y) - v_r * v_r;\n"
    "  float disc = B * B - 4.0 * A * C;\n"
    "  if (disc < 0.0 || A < 1e-9) discard;\n"
    "  float t = (-B - sqrt(disc)) / (2.0 * A);\n"
    "  vec3 hit = t * d;\n"
    "  float s = dot(hit - v_a, w);\n"
    "  if (s < 0.0 || s > L) discard;\n"
    "  if (v_dash > 0.5 && mod(floor(s / L * 7.0), 2.0) > 0.5) discard;\n"
    "  vec3 n = normalize(hit - v_a - s * w);\n"
    "  vec3 col = s < 0.5 * L ? v_ca : v_cb;\n"
    "  vec3 l = normalize(vec3(0.4, 0.6, 0.7));\n"
    "  float diff = max(dot(n, l), 0.0);\n"
    "  float spec = pow(max(dot(n, normalize(l - d)), 0.0), 48.0);\n"
    "  float fade = clamp((-hit.z - cue.x) / max(cue.y - cue.x, 1e-3), 0.0, 1.0);\n"
    "  vec3 lit = col * (0.22 + 0.78 * diff) + vec3(0.30) * spec;\n"
    "  color = vec4(mix(lit, lit * 0.45, fade), 1.0);\n"
    "  vec4 clip = proj * vec4(hit, 1.0);\n"
    "  gl_FragDepth = clip.z / clip.w * 0.5 + 0.5;\n"
    "}\n";

static GLuint shader(GLenum kind, const char *src) {
    GLuint s = glCreateShader(kind);
    glShaderSource(s, 1, &src, NULL);
    glCompileShader(s);
    GLint ok;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[4096];
        glGetShaderInfoLog(s, sizeof log, NULL, log);
        fprintf(stderr, "amview: shader: %s\n", log);
        exit(1);
    }
    return s;
}

static GLuint program(const char *vs, const char *fs) {
    GLuint p = glCreateProgram();
    glAttachShader(p, shader(GL_VERTEX_SHADER, vs));
    glAttachShader(p, shader(GL_FRAGMENT_SHADER, fs));
    glLinkProgram(p);
    GLint ok;
    glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[4096];
        glGetProgramInfoLog(p, sizeof log, NULL, log);
        fprintf(stderr, "amview: link: %s\n", log);
        exit(1);
    }
    return p;
}

static int64_t now_ns(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (int64_t)t.tv_sec * 1000000000ll + t.tv_nsec;
}

static int cmp64(const void *a, const void *b) {
    int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
    return (x > y) - (x < y);
}

static void ramp(float t, float *c) {
    t = t < 0 ? 0 : t > 1 ? 1 : t;
    for (int k = 0; k < 3; k++) {
        float v = 1.5f - fabsf(4.0f * t - (float)(3 - k));
        c[k] = v < 0 ? 0 : v > 1 ? 1 : v;
    }
}

/* ---------------------------------------------------------------- input state */
static float yaw = 0.6f, pitch = 0.35f, dist_scale = 2.4f;
static int autorot = 1, dragging;
static double lastx, lasty;
/* file mode */
static int paused, cmode, show_h, show_atoms = 1, show_tube = 1, file_step, file_switch, temp_switch, rep_switch;
static double speed = 1.0, playhead;

static void on_key(GLFWwindow *w, int key, int sc, int action, int mods) {
    (void)sc; (void)mods;
    if (action != GLFW_PRESS && action != GLFW_REPEAT) return;
    switch (key) {
    case GLFW_KEY_Q: case GLFW_KEY_ESCAPE: glfwSetWindowShouldClose(w, 1); break;
    case GLFW_KEY_O: autorot = !autorot; break;
    case GLFW_KEY_SPACE: paused = !paused; break;
    case GLFW_KEY_RIGHT_BRACKET: speed *= 2; break;
    case GLFW_KEY_LEFT_BRACKET: speed /= 2; break;
    case GLFW_KEY_PERIOD: file_step = 1; break;
    case GLFW_KEY_COMMA: file_step = -1; break;
    case GLFW_KEY_HOME: playhead = 0; break;
    case GLFW_KEY_C: cmode = (cmode + 1) % 3; break;
    case GLFW_KEY_H: show_h = !show_h; break;
    case GLFW_KEY_A: show_atoms = !show_atoms; break;
    case GLFW_KEY_B: show_tube = !show_tube; break;
    case GLFW_KEY_N: file_switch = 1; break;
    case GLFW_KEY_P: file_switch = -1; break;
    case GLFW_KEY_T: temp_switch = 1; break;
    case GLFW_KEY_Y: rep_switch = 1; break;
    default: break;
    }
}
static void on_button(GLFWwindow *w, int button, int action, int mods) {
    (void)mods;
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        dragging = action == GLFW_PRESS;
        glfwGetCursorPos(w, &lastx, &lasty);
    }
}
static void on_cursor(GLFWwindow *w, double x, double y) {
    (void)w;
    if (!dragging) return;
    yaw -= (float)(x - lastx) * 0.005f;
    pitch += (float)(y - lasty) * 0.005f;
    if (pitch > 1.5f) pitch = 1.5f;
    if (pitch < -1.5f) pitch = -1.5f;
    lastx = x; lasty = y;
}
static void on_scroll(GLFWwindow *w, double dx, double dy) {
    (void)w; (void)dx;
    dist_scale *= dy > 0 ? 0.9f : 1.0f / 0.9f;
    if (dist_scale < 0.2f) dist_scale = 0.2f;
    if (dist_scale > 8.0f) dist_scale = 8.0f;
}

/* ---------------------------------------------------------------- file mode: molecules */
typedef struct {
    mol_t m;
    int *ca;            /* C-alpha atom of each residue, -1 when none */
    int *bond, nbond;   /* ligand bonds: atom pairs, inferred from covalent radii or read from a mol2 */
    char *order;        /* per bond: 1, 2, 3, or 4 for aromatic (mol2 "ar"); 1 when inferred */
    int ligand_only;    /* no protein backbone: a drug on its own */
    float radius;       /* bounding radius of frame 0 about its centroid */
    int temps[8], ntemps, ti, rep;
} scene_t;

static float vdw(const char *e) {
    if (!e[1]) switch (e[0]) {
        case 'H': return 1.10f; case 'C': return 1.70f; case 'N': return 1.55f; case 'O': return 1.52f;
        case 'S': return 1.80f; case 'P': return 1.80f; case 'F': return 1.47f; default: break;
        }
    if (!strcmp(e, "CL")) return 1.75f;
    if (!strcmp(e, "FE")) return 1.40f;
    return 1.80f;
}

static void elem_color(const char *e, int ligand, float *c) {
    float g[3] = {0.55f, 0.55f, 0.58f};
    if (!e[1]) switch (e[0]) {
        case 'C': if (ligand) { g[0] = 0.25f; g[1] = 0.95f; g[2] = 0.35f; } break;
        case 'N': g[0] = 0.2f; g[1] = 0.35f; g[2] = 0.95f; break;
        case 'O': g[0] = 0.95f; g[1] = 0.15f; g[2] = 0.1f; break;
        case 'S': g[0] = 0.95f; g[1] = 0.85f; g[2] = 0.2f; break;
        case 'H': g[0] = 0.9f; g[1] = 0.9f; g[2] = 0.9f; break;
        case 'P': g[0] = 1.0f; g[1] = 0.55f; g[2] = 0.1f; break;
        case 'F': g[0] = 0.55f; g[1] = 0.95f; g[2] = 0.95f; break;
        default: g[0] = 0.9f; g[1] = 0.4f; g[2] = 0.75f; break;
        }
    else if (!strcmp(e, "CL")) { g[0] = 0.2f; g[1] = 0.95f; g[2] = 0.2f; }
    else if (!strcmp(e, "FE")) { g[0] = 0.85f; g[1] = 0.45f; g[2] = 0.2f; }
    else { g[0] = 0.9f; g[1] = 0.4f; g[2] = 0.75f; }
    memcpy(c, g, sizeof g);
}

static void ss_color(char s, float *c) {
    switch (s) {
    case 'H': c[0] = 0.92f; c[1] = 0.22f; c[2] = 0.45f; break;
    case 'E': c[0] = 0.97f; c[1] = 0.80f; c[2] = 0.12f; break;
    case 'T': c[0] = 0.40f; c[1] = 0.62f; c[2] = 0.95f; break;
    default: c[0] = 0.78f; c[1] = 0.78f; c[2] = 0.80f; break;
    }
}

/* covalent radii (Cordero 2008), Angstrom */
static float cov(const char *e) {
    if (!e[1]) switch (e[0]) {
        case 'H': return 0.31f; case 'C': return 0.76f; case 'N': return 0.71f; case 'O': return 0.66f;
        case 'S': return 1.05f; case 'P': return 1.07f; case 'F': return 0.57f; default: break;
        }
    if (!strcmp(e, "CL")) return 1.02f;
    return 1.2f;
}

/* bond orders from a mol2 of the same molecule, same atom order (antechamber writes the atoms in
 * the order tleap keeps): each bond already found gets the mol2 order; bonds the mol2 lists that
 * distance inference missed are added. */
static const char *bonds_mol2;
static void mol2_orders(scene_t *s, const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); return; }
    char line[256];
    int in = 0, matched = 0, added = 0;
    while (fgets(line, sizeof line, f)) {
        if (!strncmp(line, "@<TRIPOS>BOND", 13)) { in = 1; continue; }
        if (line[0] == 0x40) { in = 0; continue; }
        if (!in) continue;
        int id, a, b;
        char ty[8];
        if (sscanf(line, "%d %d %d %7s", &id, &a, &b, ty) != 4) continue;
        a--; b--;
        if (a < 0 || b < 0 || a >= s->m.natoms || b >= s->m.natoms) continue;
        char o = !strcmp(ty, "ar") ? 4 : (ty[0] >= 0x31 && ty[0] <= 0x33) ? (char)(ty[0] - 0x30) : 1;
        int found = -1;
        for (int k = 0; k < s->nbond; k++)
            if ((s->bond[2 * k] == a && s->bond[2 * k + 1] == b) || (s->bond[2 * k] == b && s->bond[2 * k + 1] == a)) { found = k; break; }
        if (found < 0) {
            s->bond = realloc(s->bond, sizeof(int) * 2 * (size_t)(s->nbond + 1));
            s->order = realloc(s->order, (size_t)(s->nbond + 1));
            s->bond[2 * s->nbond] = a; s->bond[2 * s->nbond + 1] = b; found = s->nbond++; added++;
        } else matched++;
        s->order[found] = o;
    }
    fclose(f);
    printf("amview: bond orders from %s: %d bonds matched, %d added\n", path, matched, added);
}

static int scene_load(scene_t *s, const char *path) {
    free(s->ca);
    s->ca = NULL;
    mol_free(&s->m);
    int temp = s->ntemps ? s->temps[s->ti] : 0;
    if (mol_load(path, temp, s->rep, &s->m)) return -1;
    s->ca = malloc(sizeof(int) * (size_t)(s->m.nres + 1));
    for (int r = 0; r < s->m.nres; r++) s->ca[r] = -1;
    for (int i = 0; i < s->m.natoms; i++)
        if (!s->m.het[i] && !strcmp(s->m.name[i], "CA") && s->ca[s->m.res[i]] < 0) s->ca[s->m.res[i]] = i;
    /* no protein backbone at all (a drug on its own, e.g. a tleap PDB with ATOM records): draw it as a ligand */
    int nca = 0;
    for (int r = 0; r < s->m.nres; r++) nca += s->ca[r] >= 0;
    s->ligand_only = nca == 0;
    if (nca == 0) for (int i = 0; i < s->m.natoms; i++) s->m.het[i] = 1;
    double cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < s->m.natoms; i++) { cx += s->m.xyz[3 * i]; cy += s->m.xyz[3 * i + 1]; cz += s->m.xyz[3 * i + 2]; }
    cx /= s->m.natoms; cy /= s->m.natoms; cz /= s->m.natoms;
    float r2 = 0;
    for (int i = 0; i < s->m.natoms; i++) {
        float dx = s->m.xyz[3 * i] - (float)cx, dy = s->m.xyz[3 * i + 1] - (float)cy, dz = s->m.xyz[3 * i + 2] - (float)cz;
        float d = dx * dx + dy * dy + dz * dz;
        if (d > r2) r2 = d;
    }
    s->radius = sqrtf(r2) + 2.0f;
    /* ligand bonds: pairs of ligand atoms in the same residue closer than 1.15 x the covalent radii sum */
    free(s->bond);
    s->bond = NULL; s->nbond = 0;
    int cap = 0;
    for (int i = 0; i < s->m.natoms; i++) {
        if (!s->m.het[i]) continue;
        for (int j = i + 1; j < s->m.natoms; j++) {
            if (!s->m.het[j] || s->m.res[j] != s->m.res[i]) continue;
            float dx = s->m.xyz[3 * i] - s->m.xyz[3 * j], dy = s->m.xyz[3 * i + 1] - s->m.xyz[3 * j + 1], dz = s->m.xyz[3 * i + 2] - s->m.xyz[3 * j + 2];
            float lim = 1.15f * (cov(s->m.elem[i]) + cov(s->m.elem[j]));
            if (dx * dx + dy * dy + dz * dz > lim * lim) continue;
            if (s->nbond == cap) { cap = cap ? cap * 2 : 256; s->bond = realloc(s->bond, sizeof(int) * 2 * (size_t)cap); }
            s->bond[2 * s->nbond] = i; s->bond[2 * s->nbond + 1] = j; s->nbond++;
        }
    }
    free(s->order);
    s->order = calloc((size_t)(s->nbond ? s->nbond : 1), 1);
    for (int b = 0; b < s->nbond; b++) s->order[b] = 1;
    if (bonds_mol2) mol2_orders(s, bonds_mol2);
    if (s->ligand_only) { show_h = 1; dist_scale = 2.3f; } /* a lone drug: hydrogens on; with its 25 degree view this is 4.73 x its radius, the vertical fit needs 4.51 */
    printf("amview: loaded %s: %s, %d atoms, %d residues, %d frames\n", path, s->m.title, s->m.natoms, s->m.nres, s->m.nframes);
    return 0;
}

/* the spheres of the scene at a fractional frame */
static float frame_radius; /* the drawn frame's extent about its centre, for the camera */
static int scene_build(const scene_t *s, double fpos, inst_t *out, int cap, cyl_t *cout, int ccap, int *ncyl, int *ss_counts) {
    const mol_t *m = &s->m;
    int n = m->natoms, nf = m->nframes;
    int a = (int)fpos;
    if (a >= nf - 1) a = nf - 1;
    int b = a + 1 < nf ? a + 1 : a;
    float t = (float)(fpos - a);
    if (b == a) t = 0;
    const float *pa = &m->xyz[(size_t)a * n * 3], *pb = &m->xyz[(size_t)b * n * 3];
    double cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < n; i++) {
        cx += pa[3 * i] + t * (pb[3 * i] - pa[3 * i]);
        cy += pa[3 * i + 1] + t * (pb[3 * i + 1] - pa[3 * i + 1]);
        cz += pa[3 * i + 2] + t * (pb[3 * i + 2] - pa[3 * i + 2]);
    }
    cx /= n; cy /= n; cz /= n;
    const char *ss = m->ss ? (m->ss_rows > 1 ? &m->ss[(size_t)(t < 0.5f ? a : b) * m->nres] : m->ss) : NULL;
    int mode = cmode;
    if (mode == 0 && !ss) mode = 2;
    ss_counts[0] = ss_counts[1] = 0;
    if (ss) for (int r = 0; r < m->nres; r++) { ss_counts[0] += ss[r] == 'H'; ss_counts[1] += ss[r] == 'E'; }
    int k = 0;
    *ncyl = 0;
    float col[3];
#define POS(i, d) (pa[3 * (i) + (d)] + t * (pb[3 * (i) + (d)] - pa[3 * (i) + (d)]) - (float)((d) == 0 ? cx : (d) == 1 ? cy : cz))
    for (int i = 0; i < n && k < cap; i++) {
        int lig = m->het[i];
        int hyd = m->elem[i][0] == 'H' && !m->elem[i][1];
        if (!lig && (!show_atoms || (hyd && !show_h))) continue;
        if (lig && hyd && !show_h) continue;
        if (lig || mode == 1) elem_color(m->elem[i], lig && !s->ligand_only, col);
        else if (mode == 0) ss_color(ss[m->res[i]], col);
        else ramp(m->nres > 1 ? (float)m->res[i] / (float)(m->nres - 1) : 0.5f, col);
        float r = vdw(m->elem[i]) * (lig ? 0.24f : (show_tube ? 0.28f : 0.5f)); /* ligands: ball and stick */
        out[k++] = (inst_t){POS(i, 0), POS(i, 1), POS(i, 2), r, col[0], col[1], col[2], 1};
    }
    /* ligand bonds as sticks: dense chains of small spheres, each half in its atom's colour.
     * Single: one stick. Double: two thinner sticks side by side, triple: three. Aromatic: a stick
     * and a thin dashed one beside it. "Beside" is in the plane of the bond and a neighbouring atom. */
    for (int bnd = 0; bnd < s->nbond; bnd++) {
        int i = s->bond[2 * bnd], j = s->bond[2 * bnd + 1];
        int hi = m->elem[i][0] == 'H' && !m->elem[i][1], hj = m->elem[j][0] == 'H' && !m->elem[j][1];
        if ((hi || hj) && !show_h) continue;
        float a[3] = {POS(i, 0), POS(i, 1), POS(i, 2)}, b[3] = {POS(j, 0), POS(j, 1), POS(j, 2)};
        float d[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        float len = sqrtf(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (len < 1e-3f) continue;
        float ci[3], cj[3];
        elem_color(m->elem[i], !s->ligand_only, ci);
        elem_color(m->elem[j], !s->ligand_only, cj);
        int o = s->order ? s->order[bnd] : 1;
        float p[3] = {0, 0, 0};
        if (o > 1) {
            /* a heavy neighbour of either end, not the other end */
            int nb = -1;
            for (int q = 0; q < s->nbond && nb < 0; q++) {
                int x = s->bond[2 * q], y = s->bond[2 * q + 1], other = -1;
                if (x == i && y != j) other = y; else if (y == i && x != j) other = x;
                else if (x == j && y != i) other = y; else if (y == j && x != i) other = x;
                if (other >= 0 && !(m->elem[other][0] == 'H' && !m->elem[other][1])) nb = other;
            }
            float v[3] = {1, 0, 0};
            if (nb >= 0) { v[0] = POS(nb, 0) - a[0]; v[1] = POS(nb, 1) - a[1]; v[2] = POS(nb, 2) - a[2]; }
            float t = (v[0] * d[0] + v[1] * d[1] + v[2] * d[2]) / (len * len);
            for (int c = 0; c < 3; c++) p[c] = v[c] - t * d[c];
            float pl = sqrtf(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            if (pl < 1e-4f) { p[0] = -d[1]; p[1] = d[0]; p[2] = 0; pl = sqrtf(p[0] * p[0] + p[1] * p[1]); }
            for (int c = 0; c < 3; c++) p[c] /= pl;
        }
        /* up to three sticks: offset along p, radius, dashed, shortened at both ends */
        float offs[3], rad[3];
        int dash[3], ns = 0;
        float shrink[3];
        if (o == 2) { offs[0] = -0.13f; offs[1] = 0.13f; rad[0] = rad[1] = 0.10f; dash[0] = dash[1] = 0; shrink[0] = shrink[1] = 0; ns = 2; }
        else if (o == 3) { offs[0] = 0; offs[1] = -0.17f; offs[2] = 0.17f; rad[0] = rad[1] = rad[2] = 0.08f; dash[0] = dash[1] = dash[2] = 0; shrink[0] = shrink[1] = shrink[2] = 0; ns = 3; }
        else if (o == 4) { offs[0] = 0; offs[1] = 0.27f; rad[0] = STICK_RADIUS; rad[1] = 0.065f; dash[0] = 0; dash[1] = 1; shrink[0] = 0; shrink[1] = 0.18f; ns = 2; }
        else { offs[0] = 0; rad[0] = STICK_RADIUS; dash[0] = 0; shrink[0] = 0; ns = 1; }
        for (int st = 0; st < ns && *ncyl < ccap; st++) {
            float u0 = shrink[st], u1 = 1 - shrink[st];
            cyl_t *cy = &cout[(*ncyl)++];
            cy->ax = a[0] + u0 * d[0] + offs[st] * p[0]; cy->ay = a[1] + u0 * d[1] + offs[st] * p[1]; cy->az = a[2] + u0 * d[2] + offs[st] * p[2];
            cy->bx = a[0] + u1 * d[0] + offs[st] * p[0]; cy->by = a[1] + u1 * d[1] + offs[st] * p[1]; cy->bz = a[2] + u1 * d[2] + offs[st] * p[2];
            cy->r = rad[st];
            cy->dash = (float)dash[st];
            cy->car = ci[0]; cy->cag = ci[1]; cy->cab = ci[2];
            cy->cbr = cj[0]; cy->cbg = cj[1]; cy->cbb = cj[2];
            cy->pad0 = cy->pad1 = 0;
        }
    }
    if (show_tube)
        for (int r = 0; r + 1 < m->nres; r++) {
            int i = s->ca[r], j = s->ca[r + 1];
            if (i < 0 || j < 0 || m->chain[r] != m->chain[r + 1]) continue;
            float xi = POS(i, 0), yi = POS(i, 1), zi = POS(i, 2), xj = POS(j, 0), yj = POS(j, 1), zj = POS(j, 2);
            float dd = (xj - xi) * (xj - xi) + (yj - yi) * (yj - yi) + (zj - zi) * (zj - zi);
            if (dd > 4.3f * 4.3f) continue;
            float ci[3], cj[3];
            if (mode == 0) { ss_color(ss[r], ci); ss_color(ss[r + 1], cj); }
            else if (mode == 1) { ci[0] = cj[0] = 0.75f; ci[1] = cj[1] = 0.75f; ci[2] = cj[2] = 0.78f; }
            else { ramp((float)r / (float)(m->nres - 1), ci); ramp((float)(r + 1) / (float)(m->nres - 1), cj); }
            float rad = (mode == 0 && (ss[r] == 'H' || ss[r] == 'E')) ? 0.95f : 0.65f;
            for (int q = 0; q < TUBE_STEPS && k < cap; q++) {
                float u = (float)q / TUBE_STEPS;
                out[k++] = (inst_t){xi + u * (xj - xi), yi + u * (yj - yi), zi + u * (zj - zi), rad,
                                    ci[0] + u * (cj[0] - ci[0]), ci[1] + u * (cj[1] - ci[1]), ci[2] + u * (cj[2] - ci[2]), 1};
            }
        }
    {   /* extent of what was drawn */
        float r2 = 0;
        for (int q = 0; q < k; q++) { float d2 = out[q].x * out[q].x + out[q].y * out[q].y + out[q].z * out[q].z; if (d2 > r2) r2 = d2; }
        frame_radius = sqrtf(r2) + 1.5f;
    }
#undef POS
    return k;
}

/* ---------------------------------------------------------------- file mode: following a running sander */
/* --follow TRAJ --natoms N: the topology is the first --file (a PDB of the solute, as tleap's
 * savepdb writes it); TRAJ is the ASCII mdcrd sander is writing (10F8.3, a box line after each
 * frame). New complete frames are appended to the scene; the solute is the first atoms of each. */
typedef struct {
    const char *path;
    FILE *f;
    long natoms_total;
    int lines_per_frame, got_title;
    char *pend;
    size_t plen, pcap;
    float *frame;   /* one frame being parsed: the solute's coordinates */
    long vals;      /* values read into the current frame */
    int lines;      /* lines read of the current frame */
} follow_t;

static int follow_poll(follow_t *fw, scene_t *sc) {
    if (!fw->f) {
        fw->f = fopen(fw->path, "r");
        if (!fw->f) return 0;
        long v = 3 * fw->natoms_total;
        fw->lines_per_frame = (int)((v + 9) / 10) + 1;          /* coordinates, then the box line */
        fw->frame = malloc(sizeof(float) * 3 * (size_t)sc->m.natoms);
    }
    char buf[1 << 16];
    size_t n;
    int added = 0;
    while ((n = fread(buf, 1, sizeof buf, fw->f)) > 0) {
        if (fw->plen + n + 1 > fw->pcap) { fw->pcap = (fw->plen + n + 1) * 2; fw->pend = realloc(fw->pend, fw->pcap); }
        memcpy(fw->pend + fw->plen, buf, n);
        fw->plen += n;
    }
    clearerr(fw->f);
    /* complete lines only */
    size_t start = 0;
    for (size_t i = 0; i < fw->plen; i++) {
        if (fw->pend[i] != '\n') continue;
        const char *l = fw->pend + start;
        size_t len = i - start;
        start = i + 1;
        if (!fw->got_title) { fw->got_title = 1; continue; }
        fw->lines++;
        if (fw->lines < fw->lines_per_frame) {
            for (size_t p = 0; p + 8 <= len; p += 8) {
                if (fw->vals < 3L * sc->m.natoms) {
                    char num[9];
                    memcpy(num, l + p, 8); num[8] = 0;
                    fw->frame[fw->vals] = (float)atof(num);
                }
                fw->vals++;
            }
            continue;
        }
        /* the box line closes the frame */
        mol_t *m = &sc->m;
        size_t per = (size_t)m->natoms * 3;
        m->xyz = realloc(m->xyz, sizeof(float) * per * (size_t)(m->nframes + 1));
        memcpy(&m->xyz[per * (size_t)m->nframes], fw->frame, sizeof(float) * per);
        m->nframes++;
        fw->vals = 0;
        fw->lines = 0;
        added++;
    }
    memmove(fw->pend, fw->pend + start, fw->plen - start);
    fw->plen -= start;
    return added;
}

/* ---------------------------------------------------------------- main */
int main(int argc, char **argv) {
    int card = 0, nfiles = 0, temp = 0, rep = 0;
    const char *follow_path = NULL;
    long follow_natoms = 0;
    double seconds = 0, fps_traj = 0;
    const char *snap = NULL;
    const char *files[256];
    int snapped = 0;
    float radius = 1.7f;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--card") && i + 1 < argc) card = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seconds") && i + 1 < argc) seconds = atof(argv[++i]);
        else if (!strcmp(argv[i], "--radius") && i + 1 < argc) radius = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap") && i + 1 < argc) snap = argv[++i];
        else if (!strcmp(argv[i], "--temp") && i + 1 < argc) temp = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--replica") && i + 1 < argc) rep = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--fps") && i + 1 < argc) fps_traj = atof(argv[++i]);
        else if (!strcmp(argv[i], "--follow") && i + 1 < argc) follow_path = argv[++i];
        else if (!strcmp(argv[i], "--natoms") && i + 1 < argc) follow_natoms = atol(argv[++i]);
        else if (!strcmp(argv[i], "--bonds") && i + 1 < argc) bonds_mol2 = argv[++i];
        else if (!strcmp(argv[i], "--file")) { while (i + 1 < argc && argv[i + 1][0] != '-' && nfiles < 256) files[nfiles++] = argv[++i]; }
        else {
            fprintf(stderr, "usage: amview [--card N] [--seconds S] [--radius A] [--snap FILE]\n"
                            "       amview --file F [F ...] [--temp K] [--replica R] [--fps F] [--seconds S] [--snap FILE]\n");
            return 2;
        }
    }
    const int filemode = nfiles > 0;
    const volatile struct mdstat *st = NULL;
    const volatile struct mdframe_head *fh = NULL;
    char path[64] = "";
    scene_t sc;
    follow_t fw;
    memset(&fw, 0, sizeof fw);
    int64_t follow_t0 = 0;
    memset(&sc, 0, sizeof sc);
    int cur = 0;
    if (filemode) {
        sc.rep = rep;
        if (mol_mdcath_list(files[0], sc.temps, &sc.ntemps) == 0)
            for (int k = 0; k < sc.ntemps; k++) if (sc.temps[k] == temp) sc.ti = k;
        if (scene_load(&sc, files[0])) return 1;
        if (follow_path) {
            if (follow_natoms < sc.m.natoms) { fprintf(stderr, "amview: --follow needs --natoms (the whole system, at least %d)\n", sc.m.natoms); return 2; }
            fw.path = follow_path;
            fw.natoms_total = follow_natoms;
            sc.m.dt_ns = 0.0002; /* frames every 0.2 ps (scripts/lurasidone-md.sh) */
        }
    } else {
        if (card == 0) snprintf(path, sizeof path, "/dev/shm/phi-hostmem");
        else snprintf(path, sizeof path, "/dev/shm/phi-hostmem-%d", card);
        int fd = open(path, O_RDONLY);
        if (fd < 0) { fprintf(stderr, "amview: %s: %s\n", path, strerror(errno)); return 1; }
        st = mmap(NULL, 4096, PROT_READ, MAP_SHARED, fd, (off_t)MDSTAT_PART);
        fh = mmap(NULL, MDFRAME_LIMIT, PROT_READ, MAP_SHARED, fd, (off_t)MDFRAME_PART);
        close(fd);
        if (st == MAP_FAILED || fh == MAP_FAILED) { fprintf(stderr, "amview: mmap failed\n"); return 1; }
    }

    if (!glfwInit()) { fprintf(stderr, "amview: glfwInit failed\n"); return 1; }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);
    glfwWindowHint(GLFW_SAMPLES, 4);
    GLFWwindow *win = glfwCreateWindow(WIDTH, HEIGHT, filemode ? "amview" : "amview: MD on the Xeon Phi 7120A", NULL, NULL);
    if (!win) { fprintf(stderr, "amview: no GL 4.6 window\n"); return 1; }
    glfwMakeContextCurrent(win);
    glfwSwapInterval(0);
    glfwSetKeyCallback(win, on_key);
    glfwSetMouseButtonCallback(win, on_button);
    glfwSetCursorPosCallback(win, on_cursor);
    glfwSetScrollCallback(win, on_scroll);
    printf("amview: %s, %s, %s\n", glGetString(GL_RENDERER), glGetString(GL_VERSION), filemode ? "file mode" : path);

    GLuint patom = program(atom_vs, atom_fs), pbox = program(box_vs, box_fs), pcyl = program(cyl_vs, cyl_fs);
    GLuint vao, ssbo, ssbo_cyl, q[NQ];
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenBuffers(1, &ssbo);
    glGenBuffers(1, &ssbo_cyl);
    glGenQueries(NQ, q);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glClearColor(0.03f, 0.035f, 0.05f, 1.0f);

    size_t icap = 0;
    inst_t *inst = NULL;
    size_t ccap = 0;
    cyl_t *cyls = NULL;
    int ncyl = 0;
    int ninst = 0;
    uint64_t have = 0, step = 0;
    float box[3] = {0, 0, 0};
    float *stage = NULL;
    size_t stage_cap = 0;
    int64_t t0 = now_ns(), next = t0 + PERIOD_NS, sec0 = t0, last = t0;
    int64_t iv[256];
    int niv = 0, frames = 0, newframes = 0, qn = 0, ssc[2] = {0, 0};
    double gpu_ms_sum = 0;
    int gpu_n = 0;
    uint64_t step0 = 0;
    playhead = 0;

    while (!glfwWindowShouldClose(win)) {
        glfwPollEvents();
        int64_t t = now_ns();
        double dt = (double)(t - last) * 1e-9;
        if (autorot && !dragging) yaw += (float)dt * 0.10f;
        last = t;
        float view_r;

        if (filemode) {
            /* file and dataset switches */
            int reload = 0;
            if (file_switch) { cur = (cur + file_switch + nfiles) % nfiles; file_switch = 0; sc.ntemps = 0; sc.ti = 0;
                               mol_mdcath_list(files[cur], sc.temps, &sc.ntemps); reload = 1; }
            if (temp_switch && sc.ntemps) { sc.ti = (sc.ti + 1) % sc.ntemps; reload = 1; }
            if (rep_switch && sc.ntemps) { sc.rep = (sc.rep + 1) % 5; reload = 1; }
            temp_switch = rep_switch = 0;
            if (reload) { if (scene_load(&sc, files[cur])) glfwSetWindowShouldClose(win, 1); playhead = 0; }
            if (fw.path && now_ns() - follow_t0 > 500000000ll) {
                follow_t0 = now_ns();
                int got = follow_poll(&fw, &sc);
                if (got && sc.m.nframes - 1 - playhead > 25) playhead = sc.m.nframes - 6; /* fell behind: jump close to live */
            }
            int nf = sc.m.nframes;
            double rate = fps_traj > 0 ? fps_traj : (sc.m.dt_ns > 0 ? 10.0 : 2.0);
            if (nf > 1 && !paused) playhead += dt * rate * speed;
            if (file_step) { playhead = floor(playhead) + file_step; file_step = 0; }
            if (fw.path) { if (playhead > nf - 1) playhead = nf - 1; if (playhead < 0) playhead = 0; } /* live: never wrap */
            else if (nf > 1) { while (playhead >= nf - 1) playhead -= nf - 1; while (playhead < 0) playhead += nf - 1; }
            else playhead = 0;
            size_t need = (size_t)sc.m.natoms + (size_t)sc.m.nres * TUBE_STEPS + (size_t)sc.nbond * 64 + 16;
            if (need > icap) { icap = need; inst = realloc(inst, icap * sizeof *inst); }
            if ((size_t)sc.nbond * 3 + 8 > ccap) { ccap = (size_t)sc.nbond * 3 + 8; cyls = realloc(cyls, ccap * sizeof *cyls); }
            ninst = scene_build(&sc, playhead, inst, (int)icap, cyls, (int)ccap, &ncyl, ssc);
            /* the camera follows the drawn extent, smoothed (a drug in water changes shape) */
            static float view_r_s = 0;
            if (view_r_s <= 0 || reload) view_r_s = frame_radius;
            view_r_s += (frame_radius - view_r_s) * (float)(dt * 1.5 < 1 ? dt * 1.5 : 1);
            view_r = view_r_s;
        } else {
            /* the newest whole frame from the card, if it is new */
            if (fh->magic == MDFRAME_MAGIC && fh->latest != have && fh->latest != 0) {
                uint64_t k = fh->latest, n = fh->natoms;
                const volatile struct mdframe_slot *s = &fh->slot[k % MDFRAME_SLOTS];
                if (s->seq == k && n > 0) {
                    if (n * 16 > stage_cap) { stage_cap = n * 16; free(stage); stage = malloc(stage_cap); }
                    if (n > icap) { icap = n; inst = realloc(inst, icap * sizeof *inst); }
                    const char *src = (const char *)fh + fh->data_off + (k % MDFRAME_SLOTS) * fh->slot_bytes;
                    memcpy(stage, src, n * 16);
                    uint64_t st2 = s->step;
                    float b0 = s->box[0], b1 = s->box[1], b2 = s->box[2], wm = s->wmax;
                    if (s->seq == k && fh->latest - k < MDFRAME_SLOTS - 1) {
                        float ws = wm > 0 ? 1.0f / wm : 1.0f, c[3];
                        for (uint64_t i = 0; i < n; i++) {
                            ramp(stage[4 * i + 3] * ws, c);
                            inst[i] = (inst_t){stage[4 * i] - b0 / 2, stage[4 * i + 1] - b1 / 2, stage[4 * i + 2] - b2 / 2, radius,
                                               c[0], c[1], c[2], 1};
                        }
                        ninst = (int)n; have = k; step = st2;
                        box[0] = b0; box[1] = b1; box[2] = b2;
                        newframes++;
                    }
                }
            }
            float bmax = fmaxf(box[0], fmaxf(box[1], box[2]));
            view_r = bmax > 0 ? bmax * 0.87f : 100;
        }

        /* a lone drug: narrow field of view and the camera further back, so near atoms do not balloon */
        float fov = filemode && sc.ligand_only ? 0.4363f : 0.7854f;
        float dist = dist_scale * view_r * 1.1f * (0.41421f / tanf(fov / 2));
        float aim[2] = {0, 0}; /* view-space offset of the point the camera centres */
        if (filemode && sc.ligand_only && ninst) {
            /* exact perspective fit at this rotation: an atom at offset (x, y) and depth z toward the
             * camera stays inside 92 % of the frame when D >= (|y| + r) / tan + z, and likewise for x
             * against the width; D is the largest over all atoms, smoothed */
            mat4 rot = orbit(yaw, pitch, 1.0f);
            float th = tanf(fov / 2) * 0.92f, want = 0;
            float x0 = 1e9f, x1 = -1e9f, y0 = 1e9f, y1 = -1e9f;
            for (int q = 0; q < ninst; q++) {   /* the projected box: aim at its centre, not the centroid */
                const inst_t *a = &inst[q];
                float vx = rot.m[0] * a->x + rot.m[4] * a->y + rot.m[8] * a->z;
                float vy = rot.m[1] * a->x + rot.m[5] * a->y + rot.m[9] * a->z;
                x0 = fminf(x0, vx - a->r); x1 = fmaxf(x1, vx + a->r); y0 = fminf(y0, vy - a->r); y1 = fmaxf(y1, vy + a->r);
            }
            float mx = (x0 + x1) / 2, my = (y0 + y1) / 2;
            for (int q = 0; q < ninst; q++) {
                const inst_t *a = &inst[q];
                float vx = rot.m[0] * a->x + rot.m[4] * a->y + rot.m[8] * a->z - mx;
                float vy = rot.m[1] * a->x + rot.m[5] * a->y + rot.m[9] * a->z - my;
                float vz = rot.m[2] * a->x + rot.m[6] * a->y + rot.m[10] * a->z; /* toward the camera */
                float need = fmaxf((fabsf(vy) + a->r) / th, (fabsf(vx) + a->r) / (th * (float)WIDTH / HEIGHT)) + vz + a->r;
                if (need > want) want = need;
            }
            static float aim_x, aim_y;
            float g = (float)(dt * 2.0 < 1 ? dt * 2.0 : 1);
            aim_x += (mx - aim_x) * g; aim_y += (my - aim_y) * g;
            aim[0] = aim_x; aim[1] = aim_y;
            static float fit = 0;
            fit = fit <= 0 ? want : fit + (want - fit) * (float)(dt * 2.0 < 1 ? dt * 2.0 : 1);
            dist = fit * dist_scale / 2.3f; /* the wheel still zooms: 2.3 is the ligand default */
        }
        mat4 view = orbit(yaw, pitch, dist);
        view.m[12] -= aim[0];
        view.m[13] -= aim[1];
        mat4 proj = perspective(fov, (float)WIDTH / HEIGHT, 0.5f, dist + 40.0f * view_r);
        glBeginQuery(GL_TIME_ELAPSED, q[qn % NQ]);
        glViewport(0, 0, WIDTH, HEIGHT);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        if (ninst) {
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
            glBufferData(GL_SHADER_STORAGE_BUFFER, (GLsizeiptr)((size_t)ninst * sizeof *inst), inst, GL_STREAM_DRAW);
            glUseProgram(patom);
            glUniformMatrix4fv(glGetUniformLocation(patom, "view"), 1, GL_FALSE, view.m);
            glUniformMatrix4fv(glGetUniformLocation(patom, "proj"), 1, GL_FALSE, proj.m);
            glUniform2f(glGetUniformLocation(patom, "cue"), dist - view_r * 0.6f, dist + view_r * 1.4f);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo);
            glDrawArraysInstanced(GL_TRIANGLES, 0, 6, (GLsizei)ninst);
            if (filemode && ncyl) {
                glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_cyl);
                glBufferData(GL_SHADER_STORAGE_BUFFER, (GLsizeiptr)((size_t)ncyl * sizeof *cyls), cyls, GL_STREAM_DRAW);
                glUseProgram(pcyl);
                glUniformMatrix4fv(glGetUniformLocation(pcyl, "view"), 1, GL_FALSE, view.m);
                glUniformMatrix4fv(glGetUniformLocation(pcyl, "proj"), 1, GL_FALSE, proj.m);
                glUniform2f(glGetUniformLocation(pcyl, "cue"), dist - view_r * 0.6f, dist + view_r * 1.4f);
                glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo_cyl);
                glDrawArraysInstanced(GL_TRIANGLES, 0, 6, (GLsizei)ncyl);
            }
            if (!filemode) {
                glUseProgram(pbox);
                glUniformMatrix4fv(glGetUniformLocation(pbox, "view"), 1, GL_FALSE, view.m);
                glUniformMatrix4fv(glGetUniformLocation(pbox, "proj"), 1, GL_FALSE, proj.m);
                glUniform3f(glGetUniformLocation(pbox, "half_box"), box[0] / 2, box[1] / 2, box[2] / 2);
                glDrawArrays(GL_LINES, 0, 24);
            }
        }
        glEndQuery(GL_TIME_ELAPSED);
        qn++;
        if (qn >= NQ) {
            GLuint64 ns = 0;
            GLint ready = 0;
            glGetQueryObjectiv(q[qn % NQ], GL_QUERY_RESULT_AVAILABLE, &ready);
            if (ready) {
                glGetQueryObjectui64v(q[qn % NQ], GL_QUERY_RESULT, &ns);
                gpu_ms_sum += (double)ns * 1e-6;
                gpu_n++;
            }
        }

        /* self-check, once: read back this window's frame, count pixels that are not background */
        if (snap && !snapped && ninst && now_ns() - t0 > 5000000000ll) {
            unsigned char *px = malloc((size_t)WIDTH * HEIGHT * 3);
            glPixelStorei(GL_PACK_ALIGNMENT, 1);
            glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, px);
            long lit = 0;
            for (long i = 0; i < (long)WIDTH * HEIGHT; i++)
                lit += px[3 * i] > 20 || px[3 * i + 1] > 20 || px[3 * i + 2] > 25;
            FILE *o = fopen(snap, "wb");
            if (o) {
                fprintf(o, "P6 %d %d 255\n", WIDTH, HEIGHT);
                for (int y = HEIGHT - 1; y >= 0; y--) fwrite(px + (size_t)y * WIDTH * 3, 1, (size_t)WIDTH * 3, o);
                fclose(o);
            }
            printf("snap: %s, %.1f %% of pixels covered\n", snap, 100.0 * (double)lit / ((double)WIDTH * HEIGHT));
            free(px);
            snapped = 1;
        }

        /* pace to absolute 1/60 s deadlines, then present */
        struct timespec dl = {next / 1000000000ll, next % 1000000000ll};
        clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &dl, NULL);
        glfwSwapBuffers(win);
        int64_t done = now_ns();
        if (niv < 256) iv[niv++] = done - t;
        frames++;
        next += PERIOD_NS;
        if (next < done) next = done + PERIOD_NS - (done - next) % PERIOD_NS;

        if (done - sec0 >= 1000000000ll) {
            double el = (double)(done - sec0) * 1e-9;
            qsort(iv, (size_t)niv, sizeof iv[0], cmp64);
            int i50 = niv / 2, i99 = niv * 99 / 100; /* sample indices */
            char title[320];
            if (filemode) {
                double tns = sc.m.dt_ns > 0 ? playhead * sc.m.dt_ns : 0;
                static const char *cm[3] = {"structure", "element", "rainbow"};
                snprintf(title, sizeof title, "amview: %s  frame %d of %d%s%.1f ns  helix %d strand %d of %d residues  colour %s  speed x%g%s",
                         sc.m.title, (int)playhead + 1, sc.m.nframes, sc.m.dt_ns > 0 ? "  t " : "  ", tns, ssc[0], ssc[1], sc.m.nres,
                         cm[cmode], speed, paused ? "  paused" : "");
                printf("fps %.2f  loop interval p50 %.3f p99 %.3f ms  gpu %.3f ms/frame  spheres %d  %s frame %d\n",
                       frames / el, niv ? iv[i50] * 1e-6 : 0.0, niv ? iv[i99] * 1e-6 : 0.0, gpu_n ? gpu_ms_sum / gpu_n : 0.0,
                       ninst, sc.m.title, (int)playhead + 1);
            } else {
                uint64_t sstep = st->magic == MDSTAT_MAGIC ? st->nstep : 0, slim = st->magic == MDSTAT_MAGIC ? st->nstlim : 0;
                snprintf(title, sizeof title, "amview: step %llu of %llu, %d atoms, %.1f fps, frame %llu from the card",
                         (unsigned long long)sstep, (unsigned long long)slim, ninst, frames / el, (unsigned long long)step);
                printf("fps %.2f  loop interval p50 %.3f p99 %.3f ms  gpu %.3f ms/frame  md frames %.2f/s  md step %llu (+%llu)\n",
                       frames / el, niv ? iv[i50] * 1e-6 : 0.0, niv ? iv[i99] * 1e-6 : 0.0, gpu_n ? gpu_ms_sum / gpu_n : 0.0,
                       newframes / el, (unsigned long long)sstep, (unsigned long long)(sstep - step0));
                step0 = sstep;
            }
            glfwSetWindowTitle(win, title);
            fflush(stdout);
            sec0 = done; frames = 0; newframes = 0; niv = 0; gpu_ms_sum = 0; gpu_n = 0;
        }
        if (seconds > 0 && (done - t0) * 1e-9 >= seconds) break;
    }
    glfwDestroyWindow(win);
    glfwTerminate();
    free(stage);
    free(inst);
    free(cyls);
    return 0;
}
