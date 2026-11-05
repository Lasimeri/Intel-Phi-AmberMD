/* mdstat.c: writer of the MD stats block (mdstat.h, docs/md-stats.md). */
#define _GNU_SOURCE
#include "mdstat.h"
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>

static volatile struct mdstat *blk;
static uint64_t seq;

int mdstat_open(const char *name, uint64_t natoms, uint64_t nstlim, uint64_t dt_mfs, uint64_t threads) {
    const char *e = getenv("AMBER_MDSTAT");
    if (e && e[0] == '0') return -1;
    const char *file = getenv("AMBER_MDSTAT_FILE");
    const char *off = getenv("AMBER_MDSTAT_OFFSET");
    int fd = open(file ? file : "/dev/phihost", O_RDWR | O_SYNC);
    if (fd < 0) return -1;
    void *p = mmap(NULL, 4096, PROT_READ | PROT_WRITE, MAP_SHARED, fd,
                   (off_t)(off ? strtoull(off, NULL, 0) : MDSTAT_PART));
    close(fd);
    if (p == MAP_FAILED) return -1;
    blk = p;
    *(volatile uint32_t *)((volatile char *)p + MDSTAT_STOP_OFF) = 0;
    seq = blk->seq | 1u; /* continue from whatever is there, odd = being written */
    blk->seq = seq;
    blk->version = MDSTAT_VERSION;
    blk->nstep = 0;
    blk->nstlim = nstlim;
    blk->dt_mfs = dt_mfs;
    blk->natoms = natoms;
    blk->threads = threads;
    blk->etot_mkcal = 0;
    blk->temp_mk = 0;
    blk->busy_ns = 0;
    blk->step_ns = 0;
    blk->state = MDSTAT_SETUP;
    blk->pid = (uint32_t)getpid();
    char n[32] = {0};
    if (name) strncpy(n, name, sizeof n - 1);
    for (int i = 0; i < 32; i++) blk->name[i] = n[i];
    blk->magic = MDSTAT_MAGIC;
    blk->seq = ++seq; /* even: whole */
    return 0;
}

void mdstat_step(uint64_t nstep, int64_t etot_mkcal, uint64_t temp_mk, uint64_t busy_ns, uint64_t step_ns) {
    if (!blk) return;
    blk->seq = ++seq;
    blk->nstep = nstep;
    blk->etot_mkcal = etot_mkcal;
    blk->temp_mk = temp_mk;
    blk->busy_ns = busy_ns;
    blk->step_ns = step_ns;
    blk->state = MDSTAT_RUNNING;
    blk->seq = ++seq;
}

void mdstat_close(uint32_t state) {
    if (!blk) return;
    blk->seq = ++seq;
    blk->state = state;
    blk->seq = ++seq;
    munmap((void *)blk, 4096);
    blk = NULL;
}

int mdstat_stop_requested(void) {
    if (!blk) return 0;
    return *(volatile uint32_t *)((volatile char *)blk + MDSTAT_STOP_OFF) == MDSTAT_STOP;
}
