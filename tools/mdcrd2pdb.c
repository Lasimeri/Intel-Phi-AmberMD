/* mdcrd2pdb.c: the first NSOLUTE atoms of every frame of an ASCII Amber trajectory (mdcrd,
 * 10F8.3, a box line after each frame when the run is periodic) as a multi-model PDB, using the
 * atom records of TEMPLATE.pdb (tleap savepdb of the solute) in order.
 *
 * Usage: mdcrd2pdb TEMPLATE.pdb NATOMS NSOLUTE TRAJ.mdcrd OUT.pdb
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc != 6) { fprintf(stderr, "usage: mdcrd2pdb TEMPLATE.pdb NATOMS NSOLUTE TRAJ.mdcrd OUT.pdb\n"); return 2; }
    int natoms = atoi(argv[2]), nsol = atoi(argv[3]);
    if (natoms <= 0 || nsol <= 0 || nsol > natoms) { fprintf(stderr, "mdcrd2pdb: bad counts\n"); return 2; }
    FILE *tp = fopen(argv[1], "r");
    if (!tp) { perror(argv[1]); return 1; }
    char (*rec)[81] = calloc((size_t)nsol, 81), line[256];
    int k = 0;
    while (k < nsol && fgets(line, sizeof line, tp))
        if (!strncmp(line, "ATOM  ", 6) || !strncmp(line, "HETATM", 6)) { line[strcspn(line, "\n")] = 0; snprintf(rec[k++], 81, "%-80.80s", line); }
    fclose(tp);
    if (k != nsol) { fprintf(stderr, "mdcrd2pdb: template has %d atoms, need %d\n", k, nsol); return 1; }
    FILE *in = fopen(argv[4], "r"), *out = fopen(argv[5], "w");
    if (!in || !out) { perror("mdcrd2pdb"); return 1; }
    if (!fgets(line, sizeof line, in)) return 1; /* title */
    float *c = malloc(sizeof(float) * 3 * (size_t)natoms);
    int frames = 0;
    for (;;) {
        int got = 0;
        while (got < 3 * natoms) {
            if (!fgets(line, sizeof line, in)) break;
            for (int p = 0; p + 8 <= (int)strcspn(line, "\n") && got < 3 * natoms; p += 8) {
                char num[9];
                memcpy(num, line + p, 8); num[8] = 0;
                c[got++] = (float)atof(num);
            }
        }
        if (got < 3 * natoms) break;
        /* the box line, when present, is 3 numbers on its own line: peek and skip it */
        long pos = ftell(in);
        if (fgets(line, sizeof line, in) && strcspn(line, "\n") != 24) fseek(in, pos, SEEK_SET);
        frames++;
        fprintf(out, "MODEL     %4d\n", frames);
        for (int i = 0; i < nsol; i++)
            fprintf(out, "%.30s%8.3f%8.3f%8.3f%s\n", rec[i], c[3 * i], c[3 * i + 1], c[3 * i + 2], rec[i] + 54);
        fprintf(out, "ENDMDL\n");
    }
    fprintf(out, "END\n");
    fclose(out);
    fclose(in);
    printf("mdcrd2pdb: %d frames of %d atoms\n", frames, nsol);
    return frames ? 0 : 1;
}
