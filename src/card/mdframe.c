/* mdframe.c: writer of the position frames (mdframe.h, docs/md-frames.md). */
#define _GNU_SOURCE
#include "mdframe.h"
#include <fcntl.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <unistd.h>

static volatile struct mdframe_head *head;
static char *base;
static size_t map_bytes;
static uint64_t frame;

int mdframe_open(uint64_t natoms) {
    const char *e = getenv("AMBER_MDFRAME");
    if (e && e[0] == '0') return -1;
    uint64_t slot_bytes = (natoms * 16 + 4095) & ~4095ull;
    map_bytes = 4096 + MDFRAME_SLOTS * slot_bytes;
    if (natoms == 0 || map_bytes > MDFRAME_LIMIT) return -1;
    const char *file = getenv("AMBER_MDSTAT_FILE");
    const char *off = getenv("AMBER_MDFRAME_OFFSET");
    int fd = open(file ? file : "/dev/phihost", O_RDWR | O_SYNC);
    if (fd < 0) return -1;
    void *p = mmap(NULL, map_bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd,
                   (off_t)(off ? strtoull(off, NULL, 0) : MDFRAME_PART));
    close(fd);
    if (p == MAP_FAILED) return -1;
    base = p;
    head = p;
    head->magic = 0;
    head->latest = 0;
    for (unsigned s = 0; s < MDFRAME_SLOTS; s++) head->slot[s].seq = MDFRAME_BUSY;
    head->version = MDFRAME_VERSION;
    head->natoms = natoms;
    head->nslots = MDFRAME_SLOTS;
    head->stride = 16;
    head->slot_bytes = slot_bytes;
    head->data_off = 4096;
    head->magic = MDFRAME_MAGIC;
    frame = 0;
    return 0;
}

float *mdframe_buffer(void) {
    if (!head) return NULL;
    unsigned s = (unsigned)((frame + 1) % MDFRAME_SLOTS);
    head->slot[s].seq = MDFRAME_BUSY;
    return (float *)(base + head->data_off + s * head->slot_bytes);
}

void mdframe_publish(uint64_t step, float bx, float by, float bz, float wmax) {
    if (!head) return;
    uint64_t k = frame + 1;
    unsigned s = (unsigned)(k % MDFRAME_SLOTS);
    head->slot[s].step = step;
    head->slot[s].box[0] = bx;
    head->slot[s].box[1] = by;
    head->slot[s].box[2] = bz;
    head->slot[s].wmax = wmax;
    head->slot[s].seq = k;
    head->latest = k;
    frame = k;
}

void mdframe_close(void) {
    if (!head) return;
    munmap((void *)head, map_bytes);
    head = NULL;
}
