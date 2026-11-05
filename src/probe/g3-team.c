/* g3-team.c: fork-join cost of a 120-thread team on the Xeon Phi 7120A.
 *
 * One "fork-join" = the master releases the team, every thread does a token of work
 * (adds its id into its own cache line), and the master knows all have finished.
 * Three forms are timed, interleaved, in rounds:
 *   barrier  pthread_barrier_wait twice (the plain form, musl futex barrier)
 *   counter  sense-reversing barrier on one shared atomic counter
 *   flags    master writes an epoch to one line; each worker answers in its own
 *            64-byte line; the master scans the 119 answer lines
 *   treeK    same release; answers gathered up a K-ary tree (K = 4, 8, 16): node t
 *            waits for children tK+1 .. tK+K, then answers; the master waits for its own K
 * Cycles come from rdtsc; ns from CLOCK_MONOTONIC over the whole round, so the
 * TSC rate is measured, not assumed.
 *
 * Placement (7120A CPU map): core c has CPUs 1+4c .. 4+4c for c = 0..59; core 60 is
 * CPU 0 and 241..243 and is left to the OS. Thread t runs on CPU 1 + 4*(t/2) + (t%2).
 * On a host with fewer CPUs the placement wraps (logic check only).
 *
 * Output: one line per form and round: form round n min p50 p99 max (cycles), ns/tick.
 * Exit 1 if any token sum is wrong.
 */
#define _GNU_SOURCE
#include <pthread.h>
#include <sched.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define LINE 64
#define MAXT 240

static int nthr = 120;
static int iters = 10000;

typedef struct { _Alignas(LINE) atomic_uint v; char pad[LINE - sizeof(atomic_uint)]; } line_t;

static line_t epoch;               /* flags form: release word */
static line_t answer[MAXT];        /* flags form: one line per thread */
static line_t work[MAXT];          /* token work, one line per thread */
static line_t ctr, sense;          /* counter form */
static pthread_barrier_t bar;
static atomic_int mode;            /* 0 barrier, 1 counter, 2 flags, 3..5 tree k=4,8,16, 6 quit */
static const int arity[3] = {4, 8, 16};

/* tree forms: node t gathers the answers of children t*k+1 .. t*k+k, then answers itself */
static inline void gather(int t, int k, unsigned seen) {
    int lo = t * k + 1, hi = t * k + k;
    if (hi > nthr - 1) hi = nthr - 1;
    for (int c = lo; c <= hi; c++)
        while (atomic_load_explicit(&answer[c].v, memory_order_acquire) != seen) { }
}

static inline uint64_t tsc(void) {
    uint32_t lo, hi;
    __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
    return ((uint64_t)hi << 32) | lo;
}

static void pin(int t) {
    long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    int cpu = 1 + 4 * (t / 2) + (t % 2);
    if (cpu >= ncpu) cpu = t % ncpu;
    cpu_set_t s;
    CPU_ZERO(&s);
    CPU_SET(cpu, &s);
    sched_setaffinity(0, sizeof s, &s);
}

static void counter_barrier(unsigned *local_sense) {
    unsigned ls = *local_sense ^ 1u;
    *local_sense = ls;
    if (atomic_fetch_add_explicit(&ctr.v, 1, memory_order_acq_rel) == (unsigned)nthr - 1) {
        atomic_store_explicit(&ctr.v, 0, memory_order_relaxed);
        atomic_store_explicit(&sense.v, ls, memory_order_release);
    } else {
        while (atomic_load_explicit(&sense.v, memory_order_acquire) != ls) { }
    }
}

static void *worker(void *arg) {
    int id = (int)(intptr_t)arg;
    pin(id);
    unsigned ls = 0, seen = 0;
    for (;;) {
        pthread_barrier_wait(&bar);              /* start of a measured block */
        int m = atomic_load(&mode);
        if (m == 6) return NULL;
        for (int i = 0; i < iters; i++) {
            if (m == 0) {
                pthread_barrier_wait(&bar);
                atomic_fetch_add_explicit(&work[id].v, (unsigned)id, memory_order_relaxed);
                pthread_barrier_wait(&bar);
            } else if (m == 1) {
                counter_barrier(&ls);
                atomic_fetch_add_explicit(&work[id].v, (unsigned)id, memory_order_relaxed);
                counter_barrier(&ls);
            } else if (m == 2) {
                while (atomic_load_explicit(&epoch.v, memory_order_acquire) == seen) { }
                seen++;
                atomic_fetch_add_explicit(&work[id].v, (unsigned)id, memory_order_relaxed);
                atomic_store_explicit(&answer[id].v, seen, memory_order_release);
            } else {
                while (atomic_load_explicit(&epoch.v, memory_order_acquire) == seen) { }
                seen++;
                atomic_fetch_add_explicit(&work[id].v, (unsigned)id, memory_order_relaxed);
                gather(id, arity[m - 3], seen);
                atomic_store_explicit(&answer[id].v, seen, memory_order_release);
            }
        }
        pthread_barrier_wait(&bar);              /* end of a measured block */
    }
}

static int cmp64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

int main(int argc, char **argv) {
    int rounds = 3;
    if (argc > 1) nthr = atoi(argv[1]);
    if (argc > 2) iters = atoi(argv[2]);
    if (argc > 3) rounds = atoi(argv[3]);
    if (nthr < 2 || nthr > MAXT || iters < 1) { fprintf(stderr, "usage: g3-team [threads<=240] [iters] [rounds]\n"); return 2; }

    uint64_t *dt = malloc(sizeof(uint64_t) * (size_t)iters);
    pthread_barrier_init(&bar, NULL, (unsigned)nthr);
    pthread_t th[MAXT];
    for (int t = 1; t < nthr; t++) pthread_create(&th[t], NULL, worker, (void *)(intptr_t)t);
    pin(0);

    static const char *name[6] = {"barrier", "counter", "flags", "tree4", "tree8", "tree16"};
    unsigned ls = 0, seen = 0;
    int bad = 0;
    for (int r = 1; r <= rounds; r++) {
        for (int m = 0; m < 6; m++) {
            for (int t = 0; t < nthr; t++) atomic_store(&work[t].v, 0);
            atomic_store(&mode, m);
            pthread_barrier_wait(&bar);
            struct timespec a, b;
            clock_gettime(CLOCK_MONOTONIC, &a);
            uint64_t t0 = tsc(), start = t0;
            for (int i = 0; i < iters; i++) {
                if (m == 0) {
                    pthread_barrier_wait(&bar);
                    atomic_fetch_add_explicit(&work[0].v, 0, memory_order_relaxed);
                    pthread_barrier_wait(&bar);
                } else if (m == 1) {
                    counter_barrier(&ls);
                    atomic_fetch_add_explicit(&work[0].v, 0, memory_order_relaxed);
                    counter_barrier(&ls);
                } else if (m == 2) {
                    seen++;
                    atomic_store_explicit(&epoch.v, seen, memory_order_release);
                    for (int t = 1; t < nthr; t++)
                        while (atomic_load_explicit(&answer[t].v, memory_order_acquire) != seen) { }
                } else {
                    seen++;
                    atomic_store_explicit(&epoch.v, seen, memory_order_release);
                    gather(0, arity[m - 3], seen);
                }
                uint64_t t1 = tsc();
                dt[i] = t1 - t0;
                t0 = t1;
            }
            uint64_t end = t0;
            clock_gettime(CLOCK_MONOTONIC, &b);
            pthread_barrier_wait(&bar);
            for (int t = 1; t < nthr; t++)
                if (atomic_load(&work[t].v) != (unsigned)t * (unsigned)iters) bad++;
            double ns = (double)(b.tv_sec - a.tv_sec) * 1e9 + (double)(b.tv_nsec - a.tv_nsec);
            qsort(dt, (size_t)iters, sizeof dt[0], cmp64);
            printf("%s\t%d\t%d\t%llu\t%llu\t%llu\t%llu\t%.6f\n", name[m], r, nthr,
                   (unsigned long long)dt[0], (unsigned long long)dt[iters / 2],
                   (unsigned long long)dt[(size_t)iters * 99 / 100], (unsigned long long)dt[iters - 1],
                   ns / (double)(end - start));
            fflush(stdout);
        }
    }
    atomic_store(&mode, 6);
    pthread_barrier_wait(&bar);
    for (int t = 1; t < nthr; t++) pthread_join(th[t], NULL);
    if (bad) { fprintf(stderr, "g3-team: %d wrong token sums\n", bad); return 1; }
    return 0;
}
