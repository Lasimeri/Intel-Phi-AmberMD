/* molio.h: protein structures and trajectories for amview's file mode.
 * Formats: PDB (every MODEL a frame; HELIX and SHEET records give secondary structure),
 * GROMACS .gro (every frame; nm to Angstrom), mmCIF (atom_site; every model a frame), SDF (one
 * molecule, a ligand), mdCATH HDF5 (one temperature and replica; DSSP per
 * frame). Coordinates in Angstrom. Secondary structure codes: 'H' helix (DSSP H, G, I),
 * 'E' strand (E, B), 'T' turn or bend (T, S), ' ' none or unknown. */
#ifndef MOLIO_H
#define MOLIO_H

typedef struct {
    int natoms, nres, nframes;
    char (*elem)[3];   /* element symbol, upper case */
    char (*name)[5];   /* atom name */
    char (*resn)[6];   /* residue name of the atom (mmCIF ids reach 5 characters) */
    int *res;          /* residue index of the atom, 0 .. nres-1 */
    char *chain;       /* chain of each residue */
    char *het;         /* per atom: 1 for a ligand (HETATM, not water, 6 or more atoms in its residue), else 0 */
    float *xyz;        /* nframes x natoms x 3 */
    char *ss;          /* nframes x nres, or one row shared by all frames (ss_rows 1), or NULL */
    int ss_rows;
    double dt_ns;      /* time between frames, 0 when unknown */
    char title[160];
} mol_t;

int mol_load(const char *path, int temp, int replica, mol_t *m); /* 0, or -1 with a message on stderr */
void mol_free(mol_t *m);
int mol_mdcath_list(const char *path, int *temps, int *ntemps); /* temperatures present (up to 8) */

#endif
