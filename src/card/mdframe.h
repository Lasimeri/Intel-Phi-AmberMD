/* mdframe.h: atom positions from an MD run on the card to the host viewer (amview).
 *
 * Where: the shared host-memory region, part offset MDFRAME_PART (1320 MiB + 64 KiB, after the
 * stats page of mdstat.h, inside the unclaimed [1292, 1344) MiB). A header page, then
 * MDFRAME_SLOTS slots of natoms x 16 bytes (float x, y, z in Angstrom, float w: a value the
 * viewer colours by), each slot rounded up to 4 KiB.
 *
 * Protocol (docs/md-frames.md): frame k goes to slot k % 3. The writer sets the slot's seq to
 * MDFRAME_BUSY, writes the atoms, the step and the box, sets seq = k, then latest = k. A reader
 * takes latest, copies slot latest % 3, and keeps the copy only if the slot's seq was k before
 * and after. With three slots a slot is rewritten two frames after it was published.
 */
#ifndef MDFRAME_H
#define MDFRAME_H
#include <stdint.h>

#define MDFRAME_PART ((1320ull << 20) + (64ull << 10))
#define MDFRAME_MAGIC 0x46444D41u /* "AMDF" */
#define MDFRAME_VERSION 1u
#define MDFRAME_SLOTS 3u
#define MDFRAME_BUSY UINT64_MAX
#define MDFRAME_LIMIT (52ull << 20) /* the part may not reach 1344 MiB */

struct mdframe_slot {
    uint64_t seq;   /* frame number, or MDFRAME_BUSY while being written */
    uint64_t step;  /* MD step of the frame */
    float box[3];   /* box edges, Angstrom */
    float wmax;     /* largest w in the frame (colour scale) */
    uint64_t pad[4];
};

struct mdframe_head {
    uint32_t magic;       /*  0 "AMDF" */
    uint32_t version;     /*  4 1 */
    uint64_t latest;      /*  8 last complete frame number, 0 = none yet */
    uint64_t natoms;      /* 16 */
    uint32_t nslots;      /* 24 3 */
    uint32_t stride;      /* 28 bytes per atom, 16 */
    uint64_t slot_bytes;  /* 32 bytes per slot, a multiple of 4096 */
    uint64_t data_off;    /* 40 offset of slot 0 from the part, 4096 */
    uint64_t pad[2];      /* 48 */
    struct mdframe_slot slot[MDFRAME_SLOTS]; /* 64, 64 bytes each */
};

/* Card side, callable from Fortran (bind(C)). mdframe_open returns 0 or -1 (no region, too
 * many atoms, or AMBER_MDFRAME=0). mdframe_buffer gives the slot to fill next (natoms x 4
 * floats, written directly into the shared region); mdframe_publish completes it. */
int mdframe_open(uint64_t natoms);
float *mdframe_buffer(void);
void mdframe_publish(uint64_t step, float bx, float by, float bz, float wmax);
void mdframe_close(void);

_Static_assert(sizeof(struct mdframe_slot) == 64, "slot header is 64 bytes");
_Static_assert(sizeof(struct mdframe_head) == 64 + 64 * MDFRAME_SLOTS, "header layout");

#endif
