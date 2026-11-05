# Lurasidone (+1) in water, first run, 2026-11-05

`scripts/lurasidone-md.sh` (STATE plus1-pubchem), host (Ryzen 5800X, gfortran 16.2.1 sander.MPI, 8 ranks).
Not on the Xeon Phi: sander for the card is stage G2.

| Step | Result |
|---|---|
| start | PubChem CID 213046 3D conformer, aliphatic piperazine N protonated (tools/sdf-protonate.c, atom 5) |
| antechamber | GAFF2 types, AM1-BCC charges via sqm: 72 atoms, charge sum 1.0000 |
| parmchk2 | 0 ATTN lines (no guessed terms) |
| tleap | truncated octahedron, 12 A buffer, 1,838 TIP3P, 1 Cl-: 5,584 atoms, net charge 0 |
| minimise, heat | 2,000 cycles; 20 ps to 310 K, drug restrained |
| NpT 100 ps | mean density 0.9715 g/cm3, mean T 310.29 K |
| production 200 ps | mean T 310.27 K, mean density 0.9759 g/cm3, mean Etot -13,842.01 kcal/mol |
| speed | heating 5.62 ms/step = 30.7 ns/day; production 200 ps in 674.1 s = 25.63 ns/day (amview was rendering beside it) |
| trajectory | 1,000 frames of the drug, one per 0.2 ps (multi-model PDB from tools/mdcrd2pdb.c) |

Raw outputs: logs/lurasidone/ (min, heat, npt, prod .out, antechamber.log, leap.log, lur.frcmod).
The kit's lowest-energy neutral and +1 conformers (input/dl/lurasidone-kit) run next as states
neutral-kit and plus1-kit. This is the free drug in water, not a receptor complex: no
experimental lurasidone pose exists in the PDB (docs/datasets.md).

## Neutral state from the kit (stopped by request, 2026-11-05 14:44)
`STATE=neutral-kit SDF=<kit>/inputs/lurasidone_neutral.sdf CHARGE=0 PROTONATE=0`: antechamber GAFF2/AM1-BCC
on the kit's lowest-energy neutral conformer, 71 atoms; tleap 5,894 atoms. Minimisation and heating
finished; the run was stopped during the 100 ps NpT stage, before production. The +1 run from the
kit's conformer did not start.
