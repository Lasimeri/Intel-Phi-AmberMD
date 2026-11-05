/* zipget.c: take named entries out of a remote zip without downloading the zip.
 *
 * Reads the end of central directory and the central directory with HTTP range requests, then
 * fetches each matching entry's local header and data by range, inflates it (method 8) or copies
 * it (method 0), checks its CRC-32 and writes it under OUTDIR with its base name.
 *
 * Usage: zipget URL SUBSTRING OUTDIR      (entries whose name contains SUBSTRING; "__MACOSX/"
 *                                          entries are skipped) Zip64 is not handled.
 */
#include <curl/curl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

typedef struct { unsigned char *p; size_t n, cap; } buf_t;

static size_t sink(void *d, size_t sz, size_t nm, void *u) {
    buf_t *b = u;
    size_t k = sz * nm;
    if (b->n + k > b->cap) {
        b->cap = (b->n + k) * 2;
        b->p = realloc(b->p, b->cap);
        if (!b->p) return 0;
    }
    memcpy(b->p + b->n, d, k);
    b->n += k;
    return k;
}

/* bytes [from, to] of the URL (to < 0: the last -from bytes) */
static int fetch(const char *url, long long from, long long to, buf_t *out) {
    char range[64];
    if (to < 0) snprintf(range, sizeof range, "-%lld", from);
    else snprintf(range, sizeof range, "%lld-%lld", from, to);
    CURL *c = curl_easy_init();
    out->n = 0;
    curl_easy_setopt(c, CURLOPT_URL, url);
    curl_easy_setopt(c, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(c, CURLOPT_RANGE, range);
    curl_easy_setopt(c, CURLOPT_USERAGENT, "zipget/1 (libcurl)"); /* Zenodo answers 403 without one */
    curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, sink);
    curl_easy_setopt(c, CURLOPT_WRITEDATA, out);
    CURLcode r = curl_easy_perform(c);
    long code = 0;
    curl_easy_getinfo(c, CURLINFO_RESPONSE_CODE, &code);
    curl_easy_cleanup(c);
    if (r != CURLE_OK || code != 206) {
        fprintf(stderr, "zipget: range %s: curl %d, HTTP %ld\n", range, (int)r, code);
        return -1;
    }
    return 0;
}

static unsigned u16(const unsigned char *p) { return p[0] | p[1] << 8; }
static uint32_t u32(const unsigned char *p) { return (uint32_t)p[0] | (uint32_t)p[1] << 8 | (uint32_t)p[2] << 16 | (uint32_t)p[3] << 24; }

int main(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr, "usage: zipget URL SUBSTRING OUTDIR\n"); return 2; }
    const char *url = argv[1], *want = argv[2], *dir = argv[3];
    curl_global_init(CURL_GLOBAL_DEFAULT);
    buf_t tail = {0}, cd = {0}, ent = {0};
    if (fetch(url, 65558, -1, &tail)) return 1;                 /* EOCD is within the last 64 KiB + 22 */
    long eocd = -1;
    for (long i = (long)tail.n - 22; i >= 0; i--)
        if (u32(tail.p + i) == 0x06054b50) { eocd = i; break; }
    if (eocd < 0) { fprintf(stderr, "zipget: no end of central directory\n"); return 1; }
    unsigned nent = u16(tail.p + eocd + 10);
    uint32_t cdsize = u32(tail.p + eocd + 12), cdoff = u32(tail.p + eocd + 16);
    if (cdoff == 0xffffffffu) { fprintf(stderr, "zipget: zip64 not handled\n"); return 1; }
    if (fetch(url, cdoff, (long long)cdoff + cdsize - 1, &cd)) return 1;
    printf("zipget: %u entries, central directory %u B at %u\n", nent, cdsize, cdoff);
    size_t p = 0;
    int got = 0, bad = 0;
    long long total = 0;
    for (unsigned e = 0; e < nent && p + 46 <= cd.n; e++) {
        if (u32(cd.p + p) != 0x02014b50) { fprintf(stderr, "zipget: bad central entry %u\n", e); return 1; }
        unsigned method = u16(cd.p + p + 10), nl = u16(cd.p + p + 28), xl = u16(cd.p + p + 30), cl = u16(cd.p + p + 32);
        uint32_t crc = u32(cd.p + p + 16), csize = u32(cd.p + p + 20), usize = u32(cd.p + p + 24), loff = u32(cd.p + p + 42);
        char name[512];
        snprintf(name, sizeof name, "%.*s", (int)nl, (const char *)cd.p + p + 46);
        p += 46 + nl + xl + cl;
        if (!strstr(name, want) || strstr(name, "__MACOSX/") || usize == 0) continue;
        /* local header: 30 bytes + its own name and extra lengths, then the data */
        if (fetch(url, loff, (long long)loff + 30 + nl + 1024 + csize - 1, &ent)) return 1;
        if (u32(ent.p) != 0x04034b50) { fprintf(stderr, "zipget: bad local header for %s\n", name); return 1; }
        size_t data = 30 + u16(ent.p + 26) + u16(ent.p + 28);
        unsigned char *outp = malloc(usize ? usize : 1);
        if (method == 0) memcpy(outp, ent.p + data, usize);
        else if (method == 8) {
            z_stream z = {0};
            inflateInit2(&z, -15);
            z.next_in = ent.p + data; z.avail_in = csize;
            z.next_out = outp; z.avail_out = usize;
            int zr = inflate(&z, Z_FINISH);
            inflateEnd(&z);
            if (zr != Z_STREAM_END) { fprintf(stderr, "zipget: inflate %s: %d\n", name, zr); bad++; free(outp); continue; }
        } else { fprintf(stderr, "zipget: %s: method %u not handled\n", name, method); bad++; free(outp); continue; }
        uint32_t c = (uint32_t)crc32(0L, outp, usize);
        const char *base = strrchr(name, '/');
        base = base ? base + 1 : name;
        char path[1024];
        snprintf(path, sizeof path, "%s/%s", dir, base);
        FILE *o = fopen(path, "wb");
        if (!o || fwrite(outp, 1, usize, o) != usize) { fprintf(stderr, "zipget: write %s\n", path); return 1; }
        fclose(o);
        printf("%s  %u B (%u compressed)  crc %08x %s\n", path, usize, csize, c, c == crc ? "ok" : "MISMATCH");
        bad += c != crc;
        got++;
        total += csize;
        free(outp);
    }
    printf("zipget: %d entries written, %lld compressed bytes fetched, %d bad\n", got, total, bad);
    curl_global_cleanup();
    return bad ? 1 : 0;
}
