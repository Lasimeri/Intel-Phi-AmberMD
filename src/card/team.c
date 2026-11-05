/* team.c: the thread team chosen at G3 (results/2026-11-05-g3-team.md): release by one epoch
 * word, finish by an 8-ary gather tree of per-thread answer lines; p50 10,340 cycles per
 * fork-join for 120 threads on the 7120A. Threads are pinned two per core on cores 0..59
 * (CPU 1+4c and 2+4c), core 60 stays with the OS. Each thread counts the TSC ticks it spends
 * inside jobs, so callers can report how busy the team was.
 *
 * Callable from Fortran (bind(C)):
 *   team_init(n)             start n-1 workers; the caller is thread 0. Returns n, or -1.
 *   team_run(fn, ctx)        fn(ctx, tid, n) on every thread, returns when all are done
 *   team_busy_ns()           sum of the threads' time inside jobs, in ns
 *   team_ns()                monotonic clock in ns
 *   team_end()               stop and join the workers
 */
#define _GNU_SOURCE
#include <pthread.h>
#include <sched.h>
#include <stdatomic.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>

#define LINE 64
#define MAXT 240
#define ARITY 8

typedef void (*team_fn)(void *ctx, int tid, int nthr);
typedef struct { _Alignas(LINE) atomic_uint v; uint64_t busy; } line_t;

static line_t epoch, answer[MAXT];
static int nthr = 1;
static team_fn job;
static void *jctx;
static atomic_int quit;
static unsigned master_seen;
static pthread_t th[MAXT];
static double ns_per_tick = 1.0;
/* CPUs online, read once in team_init before anything is pinned: musl answers
 * _SC_NPROCESSORS_ONLN from the calling thread's affinity mask, so a worker created after
 * the master is pinned would see 1 CPU and pin itself to CPU 0 (seen on the 7120A, 2026-11-05:
 * 119 workers on one hardware thread). */
static long ncpu = 1;

static inline uint64_t tsc(void) {
    uint32_t lo, hi;
    __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
    return ((uint64_t)hi << 32) | lo;
}

uint64_t team_ns(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (uint64_t)t.tv_sec * 1000000000u + (uint64_t)t.tv_nsec;
}

static void pin(int t) {
    int cpu = 1 + 4 * (t / 2) + (t % 2);
    if (cpu >= ncpu) cpu = t % (int)ncpu;
    cpu_set_t s;
    CPU_ZERO(&s);
    CPU_SET(cpu, &s);
    sched_setaffinity(0, sizeof s, &s);
}

static inline void gather(int t, unsigned seen) {
    int lo = t * ARITY + 1, hi = t * ARITY + ARITY;
    if (hi > nthr - 1) hi = nthr - 1;
    for (int c = lo; c <= hi; c++)
        while (atomic_load_explicit(&answer[c].v, memory_order_acquire) != seen) { }
}

static void *worker(void *arg) {
    int id = (int)(intptr_t)arg;
    pin(id);
    unsigned seen = 0;
    for (;;) {
        while (atomic_load_explicit(&epoch.v, memory_order_acquire) == seen) { }
        seen++;
        if (atomic_load_explicit(&quit, memory_order_relaxed)) return NULL;
        uint64_t t0 = tsc();
        job(jctx, id, nthr);
        answer[id].busy += tsc() - t0;
        gather(id, seen);
        atomic_store_explicit(&answer[id].v, seen, memory_order_release);
    }
}

int team_init(int n) {
    if (n < 1 || n > MAXT) return -1;
    nthr = n;
    ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    pin(0);
    uint64_t a = team_ns(), t0 = tsc();
    struct timespec d = {0, 50000000};
    nanosleep(&d, NULL);
    uint64_t b = team_ns(), t1 = tsc();
    ns_per_tick = (double)(b - a) / (double)(t1 - t0);
    for (int t = 1; t < n; t++)
        if (pthread_create(&th[t], NULL, worker, (void *)(intptr_t)t)) return -1;
    return n;
}

void team_run(team_fn fn, void *ctx) {
    job = fn;
    jctx = ctx;
    master_seen++;
    atomic_store_explicit(&epoch.v, master_seen, memory_order_release);
    uint64_t t0 = tsc();
    fn(ctx, 0, nthr);
    answer[0].busy += tsc() - t0;
    gather(0, master_seen);
}

uint64_t team_busy_ns(void) {
    uint64_t s = 0;
    for (int t = 0; t < nthr; t++) s += answer[t].busy;
    return (uint64_t)((double)s * ns_per_tick);
}

void team_end(void) {
    atomic_store(&quit, 1);
    master_seen++;
    atomic_store_explicit(&epoch.v, master_seen, memory_order_release);
    for (int t = 1; t < nthr; t++) pthread_join(th[t], NULL);
}
