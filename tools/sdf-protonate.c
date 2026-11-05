/* sdf-protonate.c: put a proton on the basic amine of a V2000 SDF molecule.
 *
 * The basic amine is a nitrogen with exactly three heavy neighbours, all carbon, each of them
 * sp3 (four neighbours) and carrying at least two hydrogens (an N(CH2)3 centre, as in an
 * alkyl piperazine). Aryl and amide nitrogens fail the test. The new hydrogen goes 1.03 A from
 * the nitrogen, opposite the mean direction to its neighbours; the nitrogen gets charge +1
 * ("M  CHG"). Prints which atom was chosen. Exit 1 unless exactly one nitrogen qualifies.
 *
 * Usage: sdf-protonate IN.sdf OUT.sdf
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXA 512
#define MAXB 1024

int main(int argc, char **argv) {
    if (argc != 3) { fprintf(stderr, "usage: sdf-protonate IN.sdf OUT.sdf\n"); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 1; }
    char head[3][256], cnt[256], atom[MAXA][128], bond[MAXB][128];
    for (int i = 0; i < 3; i++) if (!fgets(head[i], sizeof head[i], f)) return 1;
    if (!fgets(cnt, sizeof cnt, f)) return 1;
    int na = atoi((char[4]){cnt[0], cnt[1], cnt[2], 0}), nb = atoi((char[4]){cnt[3], cnt[4], cnt[5], 0});
    if (na <= 0 || na >= MAXA || nb <= 0 || nb >= MAXB) { fprintf(stderr, "sdf-protonate: bad counts\n"); return 1; }
    double x[MAXA], y[MAXA], z[MAXA];
    char el[MAXA][4];
    for (int i = 0; i < na; i++) {
        if (!fgets(atom[i], sizeof atom[i], f)) return 1;
        x[i] = strtod(atom[i], NULL);
        y[i] = strtod(atom[i] + 10, NULL);
        z[i] = strtod(atom[i] + 20, NULL);
        sscanf(atom[i] + 31, "%3s", el[i]);
    }
    int ba[MAXB], bb[MAXB], bo[MAXB];
    for (int i = 0; i < nb; i++) {
        if (!fgets(bond[i], sizeof bond[i], f)) return 1;
        ba[i] = atoi((char[4]){bond[i][0], bond[i][1], bond[i][2], 0}) - 1;
        bb[i] = atoi((char[4]){bond[i][3], bond[i][4], bond[i][5], 0}) - 1;
        bo[i] = atoi((char[4]){bond[i][6], bond[i][7], bond[i][8], 0});
    }
    /* the rest of the file (properties, data items) is copied as it is, minus M  END */
    static char tail[1 << 20];
    size_t tn = fread(tail, 1, sizeof tail - 1, f);
    tail[tn] = 0;
    fclose(f);
    int deg[MAXA] = {0}, hcnt[MAXA] = {0}, nbr[MAXA][8];
    for (int i = 0; i < nb; i++) {
        int a = ba[i], b = bb[i];
        if (deg[a] < 8) nbr[a][deg[a]] = b;
        if (deg[b] < 8) nbr[b][deg[b]] = a;
        deg[a]++; deg[b]++;
        if (!strcmp(el[a], "H")) hcnt[b]++;
        if (!strcmp(el[b], "H")) hcnt[a]++;
    }
    int pick = -1, found = 0;
    for (int i = 0; i < na; i++) {
        if (strcmp(el[i], "N") || deg[i] != 3) continue;
        int ok = 1;
        for (int k = 0; k < 3; k++) {
            int c = nbr[i][k];
            if (strcmp(el[c], "C") || deg[c] != 4 || hcnt[c] < 2) ok = 0;
        }
        for (int j = 0; j < nb; j++) if ((ba[j] == i || bb[j] == i) && bo[j] != 1) ok = 0;
        if (ok) { pick = i; found++; }
    }
    if (found != 1) { fprintf(stderr, "sdf-protonate: %d candidate nitrogens, need exactly 1\n", found); return 1; }
    double ux = 0, uy = 0, uz = 0;
    for (int k = 0; k < 3; k++) {
        int c = nbr[pick][k];
        double dx = x[c] - x[pick], dy = y[c] - y[pick], dz = z[c] - z[pick], d = sqrt(dx * dx + dy * dy + dz * dz);
        ux += dx / d; uy += dy / d; uz += dz / d;
    }
    double ul = sqrt(ux * ux + uy * uy + uz * uz);
    double hx = x[pick] - 1.03 * ux / ul, hy = y[pick] - 1.03 * uy / ul, hz = z[pick] - 1.03 * uz / ul;
    printf("sdf-protonate: atom %d (N, neighbours %d %d %d, all CH2 or CH3) protonated, H at %.4f %.4f %.4f\n",
           pick + 1, nbr[pick][0] + 1, nbr[pick][1] + 1, nbr[pick][2] + 1, hx, hy, hz);
    FILE *o = fopen(argv[2], "w");
    if (!o) { perror(argv[2]); return 1; }
    for (int i = 0; i < 3; i++) fputs(head[i], o);
    fprintf(o, "%3d%3d%s", na + 1, nb + 1, cnt + 6);
    for (int i = 0; i < na; i++) fputs(atom[i], o);
    fprintf(o, "%10.4f%10.4f%10.4f H   0  0  0  0  0  0  0  0  0  0  0  0\n", hx, hy, hz);
    for (int i = 0; i < nb; i++) fputs(bond[i], o);
    fprintf(o, "%3d%3d  1  0  0  0  0\n", pick + 1, na + 1);
    fprintf(o, "M  CHG  1 %3d   1\n", pick + 1);
    /* copy what followed the bond block, dropping any earlier charge lines (the one above replaces them) */
    char *l = tail;
    while (*l) {
        char *e = strchr(l, '\n');
        size_t len = e ? (size_t)(e - l + 1) : strlen(l);
        if (strncmp(l, "M  CHG", 6)) fwrite(l, 1, len, o);
        l += len;
    }
    fclose(o);
    return 0;
}
