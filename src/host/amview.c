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
    "uniform mat4 proj;\n"
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
    "  color = vec4(v_col * (0.22 + 0.78 * diff) + vec3(0.35) * spec + v_col * 0.25 * rim, 1.0);\n"
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
    printf("amview: loaded %s: %s, %d atoms, %d residues, %d frames\n", path, s->m.title, s->m.natoms, s->m.nres, s->m.nframes);
    return 0;
}

/* the spheres of the scene at a fractional frame */
static int scene_build(const scene_t *s, double fpos, inst_t *out, int cap, int *ss_counts) {
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
    float col[3];
#define POS(i, d) (pa[3 * (i) + (d)] + t * (pb[3 * (i) + (d)] - pa[3 * (i) + (d)]) - (float)((d) == 0 ? cx : (d) == 1 ? cy : cz))
    for (int i = 0; i < n && k < cap; i++) {
        int lig = m->het[i];
        int hyd = m->elem[i][0] == 'H' && !m->elem[i][1];
        if (!lig && (!show_atoms || (hyd && !show_h))) continue;
        if (lig && hyd && !show_h) continue;
        if (lig || mode == 1) elem_color(m->elem[i], lig, col);
        else if (mode == 0) ss_color(ss[m->res[i]], col);
        else ramp(m->nres > 1 ? (float)m->res[i] / (float)(m->nres - 1) : 0.5f, col);
        float r = vdw(m->elem[i]) * (lig ? 0.42f : (show_tube ? 0.28f : 0.5f));
        out[k++] = (inst_t){POS(i, 0), POS(i, 1), POS(i, 2), r, col[0], col[1], col[2], 1};
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
#undef POS
    return k;
}

/* ---------------------------------------------------------------- main */
int main(int argc, char **argv) {
    int card = 0, nfiles = 0, temp = 0, rep = 0;
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
    memset(&sc, 0, sizeof sc);
    int cur = 0;
    if (filemode) {
        sc.rep = rep;
        if (mol_mdcath_list(files[0], sc.temps, &sc.ntemps) == 0)
            for (int k = 0; k < sc.ntemps; k++) if (sc.temps[k] == temp) sc.ti = k;
        if (scene_load(&sc, files[0])) return 1;
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

    GLuint patom = program(atom_vs, atom_fs), pbox = program(box_vs, box_fs);
    GLuint vao, ssbo, q[NQ];
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenBuffers(1, &ssbo);
    glGenQueries(NQ, q);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glClearColor(0.03f, 0.035f, 0.05f, 1.0f);

    size_t icap = 0;
    inst_t *inst = NULL;
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
            int nf = sc.m.nframes;
            double rate = fps_traj > 0 ? fps_traj : (sc.m.dt_ns > 0 ? 10.0 : 2.0);
            if (nf > 1 && !paused) playhead += dt * rate * speed;
            if (file_step) { playhead = floor(playhead) + file_step; file_step = 0; }
            if (nf > 1) { while (playhead >= nf - 1) playhead -= nf - 1; while (playhead < 0) playhead += nf - 1; }
            else playhead = 0;
            size_t need = (size_t)sc.m.natoms + (size_t)sc.m.nres * TUBE_STEPS + 16;
            if (need > icap) { icap = need; inst = realloc(inst, icap * sizeof *inst); }
            ninst = scene_build(&sc, playhead, inst, (int)icap, ssc);
            view_r = sc.radius;
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

        mat4 view = orbit(yaw, pitch, dist_scale * view_r * 1.1f);
        mat4 proj = perspective(0.7854f, (float)WIDTH / HEIGHT, 0.5f, 40.0f * view_r);
        glBeginQuery(GL_TIME_ELAPSED, q[qn % NQ]);
        glViewport(0, 0, WIDTH, HEIGHT);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        if (ninst) {
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
            glBufferData(GL_SHADER_STORAGE_BUFFER, (GLsizeiptr)((size_t)ninst * sizeof *inst), inst, GL_STREAM_DRAW);
            glUseProgram(patom);
            glUniformMatrix4fv(glGetUniformLocation(patom, "view"), 1, GL_FALSE, view.m);
            glUniformMatrix4fv(glGetUniformLocation(patom, "proj"), 1, GL_FALSE, proj.m);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo);
            glDrawArraysInstanced(GL_TRIANGLES, 0, 6, (GLsizei)ninst);
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
    return 0;
}
