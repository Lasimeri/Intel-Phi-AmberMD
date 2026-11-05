/* amview.c: the host's 3D view of an MD run on the card, drawn by the GPU (OpenGL 4.6 core).
 *
 * Reads the frame ring (src/card/mdframe.h) and the stats block (src/card/mdstat.h) from the
 * shared host-memory region the card writes into, uploads the newest whole frame to a shader
 * storage buffer, and draws every atom as a ray-cast sphere impostor (exact sphere surface and
 * depth per pixel), coloured by the frame's w value (the run's choice: speed for ljmd), inside
 * the periodic box. Window 1920 x 1080, paced to absolute 1/60 s deadlines.
 *
 * Usage: amview [--card N] [--seconds S] [--radius A]
 *   mouse drag: rotate, wheel: zoom, space: auto-rotation on/off, q or Esc: quit
 * Every second it prints: frames drawn, frame interval p50/p99, GPU time per frame (GL timer
 * queries), new MD frames received, and the run's step and rate from the stats block.
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

#define WIDTH 1920
#define HEIGHT 1080
#define PERIOD_NS 16666667ll
#define NQ 4

typedef struct { float m[16]; } mat4;

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
    float fx = -ex / dist, fy = -ey / dist, fz = -ez / dist;          /* forward */
    float rx = fz, ry = 0, rz = -fx;                                   /* right = forward x up */
    float rl = sqrtf(rx * rx + rz * rz);
    rx /= rl; rz /= rl;
    float ux = ry * fz - rz * fy, uy = rz * fx - rx * fz, uz = rx * fy - ry * fx;
    mat4 v = {{rx, ux, -fx, 0, ry, uy, -fy, 0, rz, uz, -fz, 0,
               -(rx * ex + ry * ey + rz * ez), -(ux * ex + uy * ey + uz * ez), (fx * ex + fy * ey + fz * ez), 1}};
    return v;
}

static const char *atom_vs =
    "#version 460 core\n"
    "layout(std430, binding = 0) readonly buffer Atoms { vec4 atom[]; };\n"
    "uniform mat4 view; uniform mat4 proj; uniform float radius; uniform float wscale; uniform vec3 shift;\n"
    "out vec3 v_center; out vec3 v_pos; out float v_w;\n"
    "const vec2 corner[6] = vec2[](vec2(-1,-1), vec2(1,-1), vec2(1,1), vec2(-1,-1), vec2(1,1), vec2(-1,1));\n"
    "void main() {\n"
    "  vec4 a = atom[gl_InstanceID];\n"
    "  vec3 c = (view * vec4(a.xyz - shift, 1.0)).xyz;\n"
    "  vec3 p = c + vec3(corner[gl_VertexID] * radius * 1.5, radius);\n"
    "  v_center = c; v_pos = p; v_w = a.w * wscale;\n"
    "  gl_Position = proj * vec4(p, 1.0);\n"
    "}\n";

static const char *atom_fs =
    "#version 460 core\n"
    "layout(depth_greater) out float gl_FragDepth;\n"
    "in vec3 v_center; in vec3 v_pos; in float v_w;\n"
    "uniform mat4 proj; uniform float radius;\n"
    "out vec4 color;\n"
    "vec3 ramp(float t) { t = clamp(t, 0.0, 1.0);\n"
    "  return clamp(vec3(1.5 - abs(4.0*t - 3.0), 1.5 - abs(4.0*t - 2.0), 1.5 - abs(4.0*t - 1.0)), 0.0, 1.0); }\n"
    "void main() {\n"
    "  vec3 d = normalize(v_pos);\n"
    "  float b = dot(d, v_center);\n"
    "  float h = b*b - dot(v_center, v_center) + radius*radius;\n"
    "  if (h < 0.0) discard;\n"
    "  vec3 hit = (b - sqrt(h)) * d;\n"
    "  vec3 n = normalize(hit - v_center);\n"
    "  vec3 l = normalize(vec3(0.4, 0.6, 0.7));\n"
    "  float diff = max(dot(n, l), 0.0);\n"
    "  float spec = pow(max(dot(n, normalize(l - d)), 0.0), 48.0);\n"
    "  color = vec4(ramp(v_w) * (0.22 + 0.78 * diff) + vec3(0.35) * spec, 1.0);\n"
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

static float yaw = 0.6f, pitch = 0.35f, dist_scale = 2.4f;
static int autorot = 1, dragging;
static double lastx, lasty;

static void on_key(GLFWwindow *w, int key, int sc, int action, int mods) {
    (void)sc; (void)mods;
    if (action != GLFW_PRESS) return;
    if (key == GLFW_KEY_Q || key == GLFW_KEY_ESCAPE) glfwSetWindowShouldClose(w, 1);
    if (key == GLFW_KEY_SPACE) autorot = !autorot;
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
    if (dist_scale < 0.3f) dist_scale = 0.3f;
    if (dist_scale > 6.0f) dist_scale = 6.0f;
}

int main(int argc, char **argv) {
    int card = 0;
    double seconds = 0;
    const char *snap = NULL; /* --snap FILE: after 5 s, the frame as PPM and the covered fraction */
    int snapped = 0;
    float radius = 1.7f;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--card") && i + 1 < argc) card = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seconds") && i + 1 < argc) seconds = atof(argv[++i]);
        else if (!strcmp(argv[i], "--radius") && i + 1 < argc) radius = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap") && i + 1 < argc) snap = argv[++i];
        else { fprintf(stderr, "usage: amview [--card N] [--seconds S] [--radius A] [--snap FILE]\n"); return 2; }
    }
    char path[64];
    if (card == 0) snprintf(path, sizeof path, "/dev/shm/phi-hostmem");
    else snprintf(path, sizeof path, "/dev/shm/phi-hostmem-%d", card);
    int fd = open(path, O_RDONLY);
    if (fd < 0) { fprintf(stderr, "amview: %s: %s\n", path, strerror(errno)); return 1; }
    const volatile struct mdstat *st = mmap(NULL, 4096, PROT_READ, MAP_SHARED, fd, (off_t)MDSTAT_PART);
    const volatile struct mdframe_head *fh = mmap(NULL, MDFRAME_LIMIT, PROT_READ, MAP_SHARED, fd, (off_t)MDFRAME_PART);
    close(fd);
    if (st == MAP_FAILED || fh == MAP_FAILED) { fprintf(stderr, "amview: mmap failed\n"); return 1; }

    if (!glfwInit()) { fprintf(stderr, "amview: glfwInit failed\n"); return 1; }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);
    glfwWindowHint(GLFW_SAMPLES, 4);
    GLFWwindow *win = glfwCreateWindow(WIDTH, HEIGHT, "amview: MD on the Xeon Phi 7120A", NULL, NULL);
    if (!win) { fprintf(stderr, "amview: no GL 4.6 window\n"); return 1; }
    glfwMakeContextCurrent(win);
    glfwSwapInterval(0);
    glfwSetKeyCallback(win, on_key);
    glfwSetMouseButtonCallback(win, on_button);
    glfwSetCursorPosCallback(win, on_cursor);
    glfwSetScrollCallback(win, on_scroll);
    printf("amview: %s, %s, region %s\n", glGetString(GL_RENDERER), glGetString(GL_VERSION), path);

    GLuint patom = program(atom_vs, atom_fs), pbox = program(box_vs, box_fs);
    GLuint vao, ssbo, q[NQ];
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenBuffers(1, &ssbo);
    glGenQueries(NQ, q);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glClearColor(0.03f, 0.035f, 0.05f, 1.0f);

    uint64_t have = 0, natoms = 0, step = 0, cap = 0;
    float box[3] = {0, 0, 0}, wmax = 1;
    float *stage = NULL;
    int64_t t0 = now_ns(), next = t0 + PERIOD_NS, sec0 = t0, last = t0;
    int64_t iv[256];
    int niv = 0, frames = 0, newframes = 0, qn = 0;
    double gpu_ms_sum = 0;
    int gpu_n = 0;
    uint64_t step0 = 0;

    while (!glfwWindowShouldClose(win)) {
        glfwPollEvents();
        /* the newest whole frame, if it is new */
        if (fh->magic == MDFRAME_MAGIC && fh->latest != have && fh->latest != 0) {
            uint64_t k = fh->latest, n = fh->natoms;
            const volatile struct mdframe_slot *s = &fh->slot[k % MDFRAME_SLOTS];
            if (s->seq == k && n > 0) {
                if (n * 16 > cap) {
                    cap = n * 16;
                    free(stage);
                    stage = malloc(cap);
                    glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
                    glBufferData(GL_SHADER_STORAGE_BUFFER, (GLsizeiptr)cap, NULL, GL_DYNAMIC_DRAW);
                }
                const char *src = (const char *)fh + fh->data_off + (k % MDFRAME_SLOTS) * fh->slot_bytes;
                memcpy(stage, src, n * 16);
                uint64_t st2 = s->step;
                float b0 = s->box[0], b1 = s->box[1], b2 = s->box[2], wm = s->wmax;
                if (s->seq == k && fh->latest - k < MDFRAME_SLOTS - 1) {
                    glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
                    glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)(n * 16), stage);
                    have = k; natoms = n; step = st2;
                    box[0] = b0; box[1] = b1; box[2] = b2;
                    wmax = wm > 0 ? wm : 1;
                    newframes++;
                }
            }
        }

        int64_t t = now_ns();
        if (autorot && !dragging) yaw += (float)((t - last) * 1e-9) * 0.10f;
        last = t;
        float bmax = fmaxf(box[0], fmaxf(box[1], box[2]));
        if (bmax <= 0) bmax = 100;
        mat4 view = orbit(yaw, pitch, dist_scale * bmax);
        mat4 proj = perspective(0.7854f, (float)WIDTH / HEIGHT, 1.0f, 20.0f * bmax);

        glBeginQuery(GL_TIME_ELAPSED, q[qn % NQ]);
        glViewport(0, 0, WIDTH, HEIGHT);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        if (natoms) {
            glUseProgram(patom);
            glUniformMatrix4fv(glGetUniformLocation(patom, "view"), 1, GL_FALSE, view.m);
            glUniformMatrix4fv(glGetUniformLocation(patom, "proj"), 1, GL_FALSE, proj.m);
            glUniform1f(glGetUniformLocation(patom, "radius"), radius);
            glUniform1f(glGetUniformLocation(patom, "wscale"), 1.0f / wmax);
            glUniform3f(glGetUniformLocation(patom, "shift"), box[0] / 2, box[1] / 2, box[2] / 2);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo);
            glDrawArraysInstanced(GL_TRIANGLES, 0, 6, (GLsizei)natoms);
            glUseProgram(pbox);
            glUniformMatrix4fv(glGetUniformLocation(pbox, "view"), 1, GL_FALSE, view.m);
            glUniformMatrix4fv(glGetUniformLocation(pbox, "proj"), 1, GL_FALSE, proj.m);
            glUniform3f(glGetUniformLocation(pbox, "half_box"), box[0] / 2, box[1] / 2, box[2] / 2);
            glDrawArrays(GL_LINES, 0, 24);
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
        if (snap && !snapped && natoms && now_ns() - t0 > 5000000000ll) {
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
        if (next < done) next = done + PERIOD_NS - (done - next) % PERIOD_NS; /* fell behind: skip, keep the phase */

        if (done - sec0 >= 1000000000ll) {
            double el = (double)(done - sec0) * 1e-9;
            qsort(iv, (size_t)niv, sizeof iv[0], cmp64);
            int i50 = niv / 2, i99 = niv * 99 / 100; /* sample indices */
            uint64_t sstep = st->magic == MDSTAT_MAGIC ? st->nstep : 0, slim = st->magic == MDSTAT_MAGIC ? st->nstlim : 0;
            char title[256];
            snprintf(title, sizeof title, "amview: step %llu of %llu, %llu atoms, %.1f fps, frame %llu from the card",
                     (unsigned long long)sstep, (unsigned long long)slim, (unsigned long long)natoms, frames / el,
                     (unsigned long long)step);
            glfwSetWindowTitle(win, title);
            printf("fps %.2f  loop interval p50 %.3f p99 %.3f ms  gpu %.3f ms/frame  md frames %.2f/s  md step %llu (+%llu)\n",
                   frames / el, niv ? iv[i50] * 1e-6 : 0.0, niv ? iv[i99] * 1e-6 : 0.0,
                   gpu_n ? gpu_ms_sum / gpu_n : 0.0, newframes / el, (unsigned long long)sstep,
                   (unsigned long long)(sstep - step0));
            fflush(stdout);
            step0 = sstep;
            sec0 = done; frames = 0; newframes = 0; niv = 0; gpu_ms_sum = 0; gpu_n = 0;
        }
        if (seconds > 0 && (done - t0) * 1e-9 >= seconds) break;
    }
    glfwDestroyWindow(win);
    glfwTerminate();
    free(stage);
    return 0;
}
