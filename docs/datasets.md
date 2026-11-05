# Protein data for the 3D view

Fetched by `scripts/fetch-proteins.sh` into `input/dl/proteins/` (not in git), each file checked
against `data/*.sha256` or `data/mdcath-pick.tsv`. About 2.1 GB. View with
`out/md/amview --file FILE [FILE ...]` (see `src/host/amview.c` for the keys).

## mdCATH: proteins unfolding (CC BY 4.0)
Mirarchi, Giorgino, De Fabritiis, Sci. Data 2024, https://huggingface.co/datasets/compsciencelab/mdCATH.
All-atom MD of 5,398 CATH domains (CHARMM22*, TIP3P), 5 temperatures (320, 348, 379, 413, 450 K)
x 5 replicas, frames every 1 ns, DSSP per frame and residue. The full set is 3.61 TB (Hugging Face
listing, 2026-11-05); seven domains were taken, 2,046,901,388 B:

| Domain | Protein | Fold | Size |
|---|---|---|---|
| 1e1gA00 | human prion protein (PDB 1E1G, M166V) | 3 helices, 2-strand sheet | 505 MB |
| 2dk1A00 | WW domain of WBP-4 | 3-strand sheet | 240 MB |
| 2a28A00 | SH3 domain of yeast Bzz1 | beta barrel | 234 MB |
| 1em7A00 | protein G B1, helix variant | alpha/beta | 264 MB |
| 1y66A00 | engrailed homeodomain, designed variant | 3 helices | 262 MB |
| 1gjzA00 | ubiquitin N-terminal fragment | hairpin and helix | 255 MB |
| 2gjhA00 | Top7 C-terminal fragment (designed) | alpha/beta | 287 MB |

At 450 K the domains lose their secondary structure over the 400 to 500 ns of a replica; at 320 K
they stay folded. `amview --file F --temp 450 --replica 0`, keys `t` and `y` switch temperature and
replica.

## Chignolin folding (CC BY 4.0)
S. Iida, Zenodo 10.5281/zenodo.6349893: five NVT runs of chignolin (10 residues, Amber
ff99SB-ILDN, TIP3P, 300 K), protein only, 1,501 frames each, 0 to 300 ns, one frame per 200 ps,
plus the native structure 5AWL. Only these files were taken out of the 2,253,101,600-byte zip, with
HTTP range requests (`tools/zipget.c`): 7,360,278 bytes transferred, CRC-32 checked.

## Folded structures (RCSB PDB)
1QM0, 1QLX (human prion protein), 2PRP (hamster prion protein, NMR models), 1E1G (human PrP M166V,
NMR models), 1L2Y (Trp-cage, NMR models), 1YRF (villin headpiece HP-35), 1UBQ (ubiquitin), 1EMA
(GFP), 4HHB (deoxyhaemoglobin), 2LZM (T4 lysozyme), 1PGA (protein G B1), 1ENH (engrailed
homeodomain), 2F21 (Pin1 WW domain), 5AWL (chignolin CLN025). NMR models play as frames.

## Drugs (RCSB PDB, PubChem)
Searched by exact InChIKey through the RCSB search API (2026-11-05):

| Drug | PDB entries | Structure |
|---|---|---|
| selegiline | 2BYB only, as the covalent FAD adduct DPK; free selegiline is in no entry | human MAO-B, X-ray 2.20 A |
| dexmethylphenidate | 8Y2G only (ligand A1D5U) | human dopamine transporter, cryo-EM 2.83 A (mmCIF only) |
| lurasidone | none (exact InChIKey PQXKDMSYBGKCJA-CVTJIBDQSA-N, its first block, and name) | no experimental pose exists |
| ziprasidone (lurasidone's closest relative, same benzisothiazolyl-piperazine) | 8WTY (ligand XEF) | human NET, cryo-EM |

Also: 6CM4 (D2 receptor with risperidone), 6LUQ (D2 with haloperidol), 7XTC (5-HT7 receptor,
active, with Gs), 8WTX (NET with bupropion), and PubChem 3D conformers of selegiline (CID 26757),
dexmethylphenidate (154101), lurasidone (213046), ziprasidone (60854), risperidone (5073).

Drug molecules (HETATM residues with 6 or more atoms) are drawn larger with green carbons.
