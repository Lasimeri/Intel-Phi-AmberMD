/* mdstat.h: the MD stats block an MD program on the card writes for phitop.
 *
 * Where: the shared host-memory region (/dev/phihost on the card, /dev/shm/phi-hostmem on the
 * host), part offset MDSTAT_PART (1320 MiB, unclaimed by the other projects), first 128 bytes.
 * The card maps it uncached (O_SYNC), so the stores reach host memory in program order.
 * Protocol (docs/md-stats.md): seq is odd while the writer changes the block and even when it is
 * whole; a reader copies the block between two equal even reads of seq. All values little-endian.
 */
#ifndef MDSTAT_H
#define MDSTAT_H
#include <stdint.h>

#define MDSTAT_PART (1320ull << 20)
#define MDSTAT_MAGIC 0x53444D41u /* "AMDS" */
#define MDSTAT_VERSION 1u

enum { MDSTAT_SETUP = 0, MDSTAT_RUNNING = 1, MDSTAT_DONE = 2, MDSTAT_FAILED = 3, MDSTAT_STOPPED = 4 };

struct mdstat {
    uint32_t magic;      /*   0 "AMDS" */
    uint32_t version;    /*   4 1 */
    uint64_t seq;        /*   8 odd while being written */
    uint64_t nstep;      /*  16 steps done */
    uint64_t nstlim;     /*  24 steps requested */
    uint64_t dt_mfs;     /*  32 time step in 1/1000 fs (2 fs = 2000) */
    uint64_t natoms;     /*  40 */
    uint64_t threads;    /*  48 team size */
    int64_t etot_mkcal;  /*  56 total energy, 1/1000 kcal/mol */
    uint64_t temp_mk;    /*  64 temperature, 1/1000 K */
    uint64_t busy_ns;    /*  72 sum over team threads of time inside parallel work */
    uint64_t step_ns;    /*  80 wall time inside MD steps */
    uint32_t state;      /*  88 MDSTAT_* */
    uint32_t pid;        /*  92 */
    char name[32];       /*  96 run name, NUL padded */
};
_Static_assert(sizeof(struct mdstat) == 128, "mdstat block is 128 bytes");

/* The host asks a run to stop by writing MDSTAT_STOP at MDSTAT_STOP_OFF of the same page
 * (scripts/md-stop.sh): the one way to end a card run while `phi run` holds the control socket.
 * mdstat_open clears it; mdstat_stop_requested() reads it, once a step. */
#define MDSTAT_STOP_OFF 128
#define MDSTAT_STOP 0x504F5453u /* "STOP" */

/* Fortran calls these through bind(C). mdstat_open returns 0, or -1 when the region cannot be
 * mapped (then the other calls do nothing). AMBER_MDSTAT=0 disables it; AMBER_MDSTAT_FILE and
 * AMBER_MDSTAT_OFFSET name another file and offset (host tests). */
int mdstat_open(const char *name, uint64_t natoms, uint64_t nstlim, uint64_t dt_mfs, uint64_t threads);
void mdstat_step(uint64_t nstep, int64_t etot_mkcal, uint64_t temp_mk, uint64_t busy_ns, uint64_t step_ns);
void mdstat_close(uint32_t state);
int mdstat_stop_requested(void);

#endif
