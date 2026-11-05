/* molio.c: readers for molio.h. */
#define _GNU_SOURCE
#include "molio.h"
#include <ctype.h>
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void *xcalloc(size_t n, size_t s) {
    void *p = calloc(n ? n : 1, s);
    if (!p) { fprintf(stderr, "molio: out of memory\n"); exit(1); }
    return p;
}

static void trim(char *dst, const char *src, int n, int cap) {
    while (n > 0 && *src == ' ') { src++; n--; }
    while (n > 0 && (src[n - 1] == ' ' || src[n - 1] == '\n' || src[n - 1] == '\r')) n--;
    if (n > cap - 1) n = cap - 1;
    memcpy(dst, src, (size_t)n);
    dst[n] = 0;
}

/* element from a PDB/GRO atom name when no element column is given */
static void guess_elem(const char *name, char *el) {
    const char *p = name;
    while (*p && isdigit((unsigned char)*p)) p++;
    el[0] = (char)toupper((unsigned char)p[0]);
    el[1] = 0;
    if (!el[0]) { el[0] = 'C'; }
}

void mol_free(mol_t *m) {
    free(m->elem); free(m->name); free(m->resn); free(m->res); free(m->chain); free(m->het); free(m->xyz); free(m->ss);
    memset(m, 0, sizeof *m);
}

/* ------------------------------------------------------------------ PDB */
typedef struct { char chain; int seq; char icode; } rkey_t;

/* het: keep the flag only for residues with 6 or more HETATM atoms (drugs, cofactors), not ions */
static void mark_ligands(mol_t *m) {
    int *cnt = xcalloc((size_t)m->nres + 1, sizeof(int));
    for (int i = 0; i < m->natoms; i++) if (m->het[i]) cnt[m->res[i]]++;
    for (int i = 0; i < m->natoms; i++) m->het[i] = m->het[i] && cnt[m->res[i]] >= 6;
    free(cnt);
}

static int load_pdb_text(const char *text, mol_t *m, int want_frames) {
    /* pass 1: atoms of the first model, models count */
    int nat = 0, models = 0, in_first = 1;
    const char *l = text;
    while (*l) {
        const char *e = strchr(l, '\n');
        size_t len = e ? (size_t)(e - l) : strlen(l);
        if (len >= 6 && !strncmp(l, "MODEL ", 6)) models++;
        else if (len >= 6 && !strncmp(l, "ENDMDL", 6)) in_first = 0;
        else if (in_first && len >= 54 && (!strncmp(l, "ATOM  ", 6) || !strncmp(l, "HETATM", 6))) {
            char rn[4]; trim(rn, l + 17, 3, 4);
            /* the same atoms pass 2 keeps: not water, alternate location blank or A */
            if ((!strncmp(l, "ATOM  ", 6) || strcmp(rn, "HOH")) && (l[16] == ' ' || l[16] == 'A')) nat++;
        }
        if (!e) break;
        l = e + 1;
    }
    if (nat == 0) { fprintf(stderr, "molio: no atoms\n"); return -1; }
    if (models == 0) models = 1;
    if (!want_frames) models = 1;
    m->natoms = nat;
    m->nframes = models;
    m->elem = xcalloc((size_t)nat, sizeof *m->elem);
    m->name = xcalloc((size_t)nat, sizeof *m->name);
    m->resn = xcalloc((size_t)nat, sizeof *m->resn);
    m->res = xcalloc((size_t)nat, sizeof *m->res);
    m->chain = xcalloc((size_t)nat, 1);
    m->het = xcalloc((size_t)nat, 1);
    m->xyz = xcalloc((size_t)nat * models * 3, sizeof(float));
    rkey_t *keys = xcalloc((size_t)nat, sizeof *keys);
    /* pass 2: atoms (first model) and coordinates (every model) */
    int frame = 0, i = 0, nres = 0;
    rkey_t last = {0, -999999, 0};
    l = text;
    while (*l) {
        const char *e = strchr(l, '\n');
        size_t len = e ? (size_t)(e - l) : strlen(l);
        if (len >= 6 && !strncmp(l, "ENDMDL", 6)) {
            if (i > 0) frame++;
            i = 0;
            if (frame >= models) break;
        } else if (len >= 54 && (!strncmp(l, "ATOM  ", 6) || !strncmp(l, "HETATM", 6))) {
            char rn[4]; trim(rn, l + 17, 3, 4);
            if ((!strncmp(l, "ATOM  ", 6) || strcmp(rn, "HOH")) && i < nat) {
                if (l[16] != ' ' && l[16] != 'A') { if (!e) break; l = e + 1; continue; } /* altloc: keep A */
                char num[9];
                float c[3];
                for (int k = 0; k < 3; k++) { memcpy(num, l + 30 + 8 * k, 8); num[8] = 0; c[k] = (float)atof(num); }
                memcpy(&m->xyz[((size_t)frame * nat + i) * 3], c, sizeof c);
                if (frame == 0) {
                    trim(m->name[i], l + 12, 4, 5);
                    strcpy(m->resn[i], rn);
                    if (len >= 78 && l[76] != ' ' && isalpha((unsigned char)l[77])) {
                        m->elem[i][0] = l[76] == ' ' ? l[77] : l[76];
                        m->elem[i][1] = l[76] == ' ' ? 0 : (char)toupper((unsigned char)l[77]);
                    } else if (len >= 78 && isalpha((unsigned char)l[77])) {
                        m->elem[i][0] = l[77];
                    } else guess_elem(m->name[i], m->elem[i]);
                    char seq[5]; memcpy(seq, l + 22, 4); seq[4] = 0;
                    rkey_t k = {l[21], atoi(seq), l[26]};
                    if (k.chain != last.chain || k.seq != last.seq || k.icode != last.icode) { keys[nres] = k; m->chain[nres] = k.chain; nres++; last = k; }
                    m->res[i] = nres - 1;
                    m->het[i] = !strncmp(l, "HETATM", 6);
                }
                i++;
            }
        }
        if (!e) break;
        l = e + 1;
    }
    if (i > 0 && frame < models) frame++;
    m->nframes = frame > 0 ? frame : 1;
    m->nres = nres;
    if (m->natoms != nat) m->natoms = nat;
    /* secondary structure from HELIX / SHEET records */
    m->ss = xcalloc((size_t)nres, 1);
    memset(m->ss, ' ', (size_t)nres);
    m->ss_rows = 1;
    int any = 0;
    l = text;
    while (*l) {
        const char *e = strchr(l, '\n');
        size_t len = e ? (size_t)(e - l) : strlen(l);
        char ch = 0, ch2 = 0, code = 0;
        int a = 0, b = 0;
        char t[6];
        if (len >= 37 && !strncmp(l, "HELIX ", 6)) { ch = l[19]; memcpy(t, l + 21, 4); t[4] = 0; a = atoi(t); ch2 = l[31]; memcpy(t, l + 33, 4); b = atoi(t); code = 'H'; }
        if (len >= 37 && !strncmp(l, "SHEET ", 6)) { ch = l[21]; memcpy(t, l + 22, 4); t[4] = 0; a = atoi(t); ch2 = l[32]; memcpy(t, l + 33, 4); b = atoi(t); code = 'E'; }
        if (code && ch == ch2)
            for (int r = 0; r < nres; r++)
                if (keys[r].chain == ch && keys[r].seq >= a && keys[r].seq <= b) { m->ss[r] = code; any = 1; }
        if (!e) break;
        l = e + 1;
    }
    if (!any) { free(m->ss); m->ss = NULL; m->ss_rows = 0; }
    free(keys);
    mark_ligands(m);
    return 0;
}

static char *slurp(const char *path, size_t *n) {
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); return NULL; }
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *b = xcalloc((size_t)sz + 1, 1);
    if (fread(b, 1, (size_t)sz, f) != (size_t)sz) { fclose(f); free(b); return NULL; }
    fclose(f);
    if (n) *n = (size_t)sz;
    return b;
}

/* ------------------------------------------------------------------ GRO */
static int load_gro(const char *path, mol_t *m) {
    size_t n;
    char *t = slurp(path, &n);
    if (!t) return -1;
    /* count frames: a title line, an atom count, atoms, a box line */
    char *p = t;
    int nat = 0, frames = 0;
    double t0 = -1, t1 = -1;
    while (*p) {
        char *title = p;
        char *e = strchr(p, '\n'); if (!e) break;
        char *cnt = e + 1;
        int k = atoi(cnt);
        if (k <= 0) break;
        if (!nat) nat = k;
        if (k != nat) { fprintf(stderr, "molio: %s: frames differ in size\n", path); free(t); return -1; }
        char *ts = strstr(title, "t=");
        if (ts && ts < e) { if (frames == 0) t0 = atof(ts + 2); else if (frames == 1) t1 = atof(ts + 2); }
        p = cnt;
        for (int i = 0; i < k + 2 && *p; i++) { char *q = strchr(p, '\n'); p = q ? q + 1 : p + strlen(p); }
        frames++;
    }
    if (!nat) { fprintf(stderr, "molio: %s: no frames\n", path); free(t); return -1; }
    m->natoms = nat; m->nframes = frames;
    m->elem = xcalloc((size_t)nat, sizeof *m->elem);
    m->name = xcalloc((size_t)nat, sizeof *m->name);
    m->resn = xcalloc((size_t)nat, sizeof *m->resn);
    m->res = xcalloc((size_t)nat, sizeof *m->res);
    m->chain = xcalloc((size_t)nat, 1);
    m->het = xcalloc((size_t)nat, 1);
    m->xyz = xcalloc((size_t)nat * frames * 3, sizeof(float));
    p = t;
    int nres = 0, lastres = -1;
    for (int f = 0; f < frames; f++) {
        char *q = strchr(p, '\n'); p = q + 1;   /* title */
        q = strchr(p, '\n'); p = q + 1;         /* count */
        for (int i = 0; i < nat; i++) {
            q = strchr(p, '\n');
            char num[9];
            for (int k = 0; k < 3; k++) { memcpy(num, p + 20 + 8 * k, 8); num[8] = 0; m->xyz[((size_t)f * nat + i) * 3 + k] = 10.0f * (float)atof(num); }
            if (f == 0) {
                char rs[6]; memcpy(rs, p, 5); rs[5] = 0;
                int r = atoi(rs);
                trim(m->resn[i], p + 5, 5, 4);
                trim(m->name[i], p + 10, 5, 5);
                guess_elem(m->name[i], m->elem[i]);
                if (r != lastres) { m->chain[nres] = 'A'; nres++; lastres = r; }
                m->res[i] = nres - 1;
            }
            p = q ? q + 1 : p + strlen(p);
        }
        q = strchr(p, '\n'); p = q ? q + 1 : p + strlen(p); /* box */
    }
    m->nres = nres;
    m->dt_ns = (t0 >= 0 && t1 > t0) ? (t1 - t0) * 1e-3 : 0;
    const char *b = strrchr(path, '/');
    snprintf(m->title, sizeof m->title, "%s", b ? b + 1 : path);
    free(t);
    return 0;
}

/* ------------------------------------------------------------------ mmCIF */
/* next whitespace-separated token of a line, quotes honoured; returns its length, 0 at the end */
static int cif_tok(const char **p, const char **start) {
    const char *s = *p;
    while (*s == ' ' || *s == '\t') s++;
    if (!*s || *s == '\n' || *s == '\r') { *p = s; return 0; }
    if (*s == '\'' || *s == '"') {
        char q = *s++;
        const char *e = s;
        while (*e && *e != '\n' && !(*e == q && (e[1] == ' ' || e[1] == '\t' || e[1] == '\n' || e[1] == '\r' || !e[1]))) e++;
        *start = s;
        *p = *e == q ? e + 1 : e;
        return (int)(e - s);
    }
    const char *e = s;
    while (*e && *e != ' ' && *e != '\t' && *e != '\n' && *e != '\r') e++;
    *start = s;
    *p = e;
    return (int)(e - s);
}

static int load_cif(const char *path, mol_t *m) {
    size_t n;
    char *t = slurp(path, &n);
    if (!t) return -1;
    enum { GRP, SYM, ATM, ALT, COMP, ASYM, SEQ, INS, X, Y, Z, MODEL, NCOL };
    static const char *want[NCOL] = {"group_PDB", "type_symbol", "auth_atom_id", "label_alt_id", "auth_comp_id",
                                     "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z",
                                     "pdbx_PDB_model_num"};
    int col[NCOL], ncols = 0;
    for (int k = 0; k < NCOL; k++) col[k] = -1;
    char *p = strstr(t, "\n_atom_site.");
    if (!p) { fprintf(stderr, "molio: %s: no atom_site loop\n", path); free(t); return -1; }
    p++;
    while (!strncmp(p, "_atom_site.", 11)) {
        char *e = strchr(p, '\n');
        int len = (int)(e - p) - 11;
        while (len > 0 && (p[11 + len - 1] == ' ' || p[11 + len - 1] == '\r')) len--;
        for (int k = 0; k < NCOL; k++)
            if ((int)strlen(want[k]) == len && !strncmp(p + 11, want[k], (size_t)len)) col[k] = ncols;
        ncols++;
        p = e + 1;
    }
    if (col[X] < 0 || col[COMP] < 0 || col[SYM] < 0) { fprintf(stderr, "molio: %s: atom_site lacks columns\n", path); free(t); return -1; }
    /* two passes over the rows: count first-model atoms and models, then fill */
    char *rows = p;
    int nat = 0, models = 0, first_model = -1, last_model = -1;
    for (int pass = 0; pass < 2; pass++) {
        if (pass == 1) {
            m->natoms = nat; m->nframes = models ? models : 1;
            m->elem = xcalloc((size_t)nat, sizeof *m->elem);
            m->name = xcalloc((size_t)nat, sizeof *m->name);
            m->resn = xcalloc((size_t)nat, sizeof *m->resn);
            m->res = xcalloc((size_t)nat, sizeof *m->res);
            m->chain = xcalloc((size_t)nat, 1);
            m->het = xcalloc((size_t)nat, 1);
            m->xyz = xcalloc((size_t)nat * m->nframes * 3, sizeof(float));
        }
        int i = 0, frame = -1, nres = 0;
        char lastkey[32] = "";
        last_model = -999;
        for (const char *l = rows; *l && *l != '#' && *l != '_' && strncmp(l, "loop_", 5);) {
            const char *f[64] = {0};
            int fl[64] = {0}, k = 0;
            const char *q = l, *s;
            int len;
            while (k < 64 && (len = cif_tok(&q, &s)) > 0) { f[k] = s; fl[k] = len; k++; }
            const char *e = strchr(q, '\n');
            l = e ? e + 1 : q + strlen(q);
            if (k < ncols) continue;
            char comp[8];
            snprintf(comp, sizeof comp, "%.*s", fl[col[COMP]], f[col[COMP]]);
            if (!strcmp(comp, "HOH")) continue;
            if (col[ALT] >= 0 && !(fl[col[ALT]] == 1 && (f[col[ALT]][0] == '.' || f[col[ALT]][0] == 'A'))) continue;
            int model = col[MODEL] >= 0 ? atoi(f[col[MODEL]]) : 1;
            if (model != last_model) { frame++; i = 0; last_model = model; if (pass == 0) { models++; if (first_model < 0) first_model = model; } }
            if (pass == 0) { if (model == first_model) nat++; continue; }
            if (frame >= m->nframes || i >= nat) continue;
            float *c = &m->xyz[((size_t)frame * nat + i) * 3];
            for (int d = 0; d < 3; d++) { char num[24]; snprintf(num, sizeof num, "%.*s", fl[col[X + d]], f[col[X + d]]); c[d] = (float)atof(num); }
            if (frame == 0) {
                snprintf(m->resn[i], 6, "%.5s", comp);
                snprintf(m->name[i], 5, "%.*s", fl[col[ATM]] < 4 ? fl[col[ATM]] : 4, f[col[ATM]]);
                snprintf(m->elem[i], 3, "%.*s", fl[col[SYM]] < 2 ? fl[col[SYM]] : 2, f[col[SYM]]);
                for (char *x = m->elem[i]; *x; x++) *x = (char)toupper((unsigned char)*x);
                char key[32];
                snprintf(key, sizeof key, "%.*s/%.*s/%.*s", fl[col[ASYM]], f[col[ASYM]], fl[col[SEQ]], f[col[SEQ]],
                         col[INS] >= 0 ? fl[col[INS]] : 0, col[INS] >= 0 ? f[col[INS]] : "");
                if (strcmp(key, lastkey)) { m->chain[nres] = f[col[ASYM]][0]; nres++; snprintf(lastkey, sizeof lastkey, "%s", key); }
                m->res[i] = nres - 1;
                m->het[i] = fl[col[GRP]] == 6 && !strncmp(f[col[GRP]], "HETATM", 6);
            }
            i++;
        }
        if (pass == 1) m->nres = nres;
    }
    free(t);
    if (!m->natoms) { fprintf(stderr, "molio: %s: no atoms\n", path); return -1; }
    mark_ligands(m);
    const char *b = strrchr(path, '/');
    snprintf(m->title, sizeof m->title, "%s%s", b ? b + 1 : path, m->nframes > 1 ? " (models)" : "");
    return 0;
}

/* ------------------------------------------------------------------ SDF */
static int load_sdf(const char *path, mol_t *m) {
    char *t = slurp(path, NULL);
    if (!t) return -1;
    char *l[4] = {t, 0, 0, 0};
    for (int k = 1; k < 4; k++) { char *e = strchr(l[k - 1], '\n'); if (!e) { free(t); return -1; } l[k] = e + 1; }
    char cnt[4]; memcpy(cnt, l[3], 3); cnt[3] = 0;
    int nat = atoi(cnt);
    if (nat <= 0) { fprintf(stderr, "molio: %s: no atoms\n", path); free(t); return -1; }
    m->natoms = nat; m->nframes = 1; m->nres = 1;
    m->elem = xcalloc((size_t)nat, sizeof *m->elem);
    m->name = xcalloc((size_t)nat, sizeof *m->name);
    m->resn = xcalloc((size_t)nat, sizeof *m->resn);
    m->res = xcalloc((size_t)nat, sizeof *m->res);
    m->chain = xcalloc(1, 1);
    m->het = xcalloc((size_t)nat, 1);
    m->xyz = xcalloc((size_t)nat * 3, sizeof(float));
    char *p = strchr(l[3], '\n') + 1;
    for (int i = 0; i < nat && p; i++) {
        char num[11];
        for (int d = 0; d < 3; d++) { memcpy(num, p + 10 * d, 10); num[10] = 0; m->xyz[i * 3 + d] = (float)atof(num); }
        trim(m->elem[i], p + 31, 3, 3);
        for (char *x = m->elem[i]; *x; x++) *x = (char)toupper((unsigned char)*x);
        snprintf(m->name[i], 5, "%s", m->elem[i]);
        strcpy(m->resn[i], "LIG");
        m->het[i] = 1;
        char *e = strchr(p, '\n');
        p = e ? e + 1 : NULL;
    }
    m->chain[0] = 'L';
    const char *b = strrchr(path, '/');
    snprintf(m->title, sizeof m->title, "%s (one molecule)", b ? b + 1 : path);
    free(t);
    return 0;
}

/* ------------------------------------------------------------------ mdCATH */
static int read_vlen_strings(hid_t file, const char *ds_path, char ***out, hsize_t *count) {
    hid_t d = H5Dopen2(file, ds_path, H5P_DEFAULT);
    if (d < 0) return -1;
    hid_t s = H5Dget_space(d);
    hsize_t dims[2] = {1, 1};
    int rank = H5Sget_simple_extent_dims(s, dims, NULL);
    hsize_t n = rank == 2 ? dims[0] * dims[1] : (rank == 1 ? dims[0] : 1);
    hid_t ty = H5Tcopy(H5T_C_S1);
    H5Tset_size(ty, H5T_VARIABLE);
    char **v = xcalloc((size_t)n, sizeof(char *));
    herr_t r = H5Dread(d, ty, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);
    H5Tclose(ty); H5Sclose(s); H5Dclose(d);
    if (r < 0) { free(v); return -1; }
    *out = v; *count = n;
    return 0;
}

static void free_vlen(char **v, hsize_t n) {
    for (hsize_t i = 0; i < n; i++) H5free_memory(v[i]);
    free(v);
}

static int mdcath_domain(hid_t file, char *dom, size_t cap) {
    H5G_info_t gi;
    if (H5Gget_info(file, &gi) < 0 || gi.nlinks < 1) return -1;
    ssize_t k = H5Lget_name_by_idx(file, ".", H5_INDEX_NAME, H5_ITER_INC, 0, dom, cap, H5P_DEFAULT);
    return k > 0 ? 0 : -1;
}

int mol_mdcath_list(const char *path, int *temps, int *ntemps) {
    *ntemps = 0;
    hid_t f = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f < 0) return -1;
    char dom[64];
    if (mdcath_domain(f, dom, sizeof dom)) { H5Fclose(f); return -1; }
    hid_t g = H5Gopen2(f, dom, H5P_DEFAULT);
    H5G_info_t gi;
    H5Gget_info(g, &gi);
    for (hsize_t i = 0; i < gi.nlinks && *ntemps < 8; i++) {
        char nm[32];
        H5Lget_name_by_idx(g, ".", H5_INDEX_NAME, H5_ITER_INC, i, nm, sizeof nm, H5P_DEFAULT);
        if (isdigit((unsigned char)nm[0])) temps[(*ntemps)++] = atoi(nm);
    }
    H5Gclose(g); H5Fclose(f);
    return 0;
}

static int load_mdcath(const char *path, int temp, int rep, mol_t *m) {
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
    hid_t f = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f < 0) { fprintf(stderr, "molio: %s: not HDF5\n", path); return -1; }
    char dom[64], ds[160];
    if (mdcath_domain(f, dom, sizeof dom)) { H5Fclose(f); return -1; }
    /* atom names, residues, elements from the embedded PDB of the protein atoms */
    char **pdb; hsize_t one;
    snprintf(ds, sizeof ds, "/%s/pdbProteinAtoms", dom);
    if (read_vlen_strings(f, ds, &pdb, &one)) { fprintf(stderr, "molio: no pdbProteinAtoms\n"); H5Fclose(f); return -1; }
    int r = load_pdb_text(pdb[0], m, 0);
    free_vlen(pdb, one);
    if (r) { H5Fclose(f); return -1; }
    free(m->ss); m->ss = NULL; m->ss_rows = 0;
    /* coordinates of the temperature and replica */
    snprintf(ds, sizeof ds, "/%s/%d/%d/coords", dom, temp, rep);
    hid_t d = H5Dopen2(f, ds, H5P_DEFAULT);
    if (d < 0) { fprintf(stderr, "molio: %s: no %s\n", path, ds); H5Fclose(f); return -1; }
    hid_t s = H5Dget_space(d);
    hsize_t dims[3];
    H5Sget_simple_extent_dims(s, dims, NULL);
    if ((int)dims[1] != m->natoms) { fprintf(stderr, "molio: %d atoms in the PDB, %llu in coords\n", m->natoms, (unsigned long long)dims[1]); H5Sclose(s); H5Dclose(d); H5Fclose(f); return -1; }
    m->nframes = (int)dims[0];
    free(m->xyz);
    m->xyz = xcalloc((size_t)dims[0] * dims[1] * 3, sizeof(float));
    H5Dread(d, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m->xyz);
    H5Sclose(s); H5Dclose(d);
    /* DSSP per frame and residue */
    char **dssp; hsize_t nd;
    snprintf(ds, sizeof ds, "/%s/%d/%d/dssp", dom, temp, rep);
    if (!read_vlen_strings(f, ds, &dssp, &nd) && nd == (hsize_t)m->nframes * (hsize_t)m->nres) {
        m->ss = xcalloc((size_t)nd, 1);
        m->ss_rows = m->nframes;
        for (hsize_t i = 0; i < nd; i++) {
            char c = dssp[i] && dssp[i][0] ? dssp[i][0] : ' ';
            m->ss[i] = (c == 'H' || c == 'G' || c == 'I') ? 'H' : (c == 'E' || c == 'B') ? 'E' : (c == 'T' || c == 'S') ? 'T' : ' ';
        }
        free_vlen(dssp, nd);
    }
    H5Fclose(f);
    m->dt_ns = 1.0;
    snprintf(m->title, sizeof m->title, "mdCATH %s, %d K, replica %d", dom, temp, rep);
    return 0;
}

int mol_load(const char *path, int temp, int rep, mol_t *m) {
    memset(m, 0, sizeof *m);
    size_t n = strlen(path);
    if (n > 3 && (!strcmp(path + n - 3, ".h5") || !strcmp(path + n - 5, ".hdf5"))) return load_mdcath(path, temp, rep, m);
    if (n > 4 && !strcmp(path + n - 4, ".gro")) return load_gro(path, m);
    if (n > 4 && !strcmp(path + n - 4, ".cif")) return load_cif(path, m);
    if (n > 4 && (!strcmp(path + n - 4, ".sdf") || !strcmp(path + n - 4, ".mol"))) return load_sdf(path, m);
    char *t = slurp(path, NULL);
    if (!t) return -1;
    int r = load_pdb_text(t, m, 1);
    free(t);
    if (!r) {
        const char *b = strrchr(path, '/');
        snprintf(m->title, sizeof m->title, "%s%s", b ? b + 1 : path, m->nframes > 1 ? " (NMR models)" : "");
    }
    return r;
}
