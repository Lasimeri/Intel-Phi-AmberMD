#!/bin/bash
# fetch-proteins.sh: the protein data for amview's file mode, into input/dl/proteins (git-ignored),
# each file checked against data/*.sha256 or data/mdcath-pick.tsv. About 2.1 GB. See docs/datasets.md.
set -euo pipefail
R=$(cd "$(dirname "$0")/.." && pwd)
D=$R/input/dl/proteins
mkdir -p "$D"/{mdcath,pdb,drugs,chignolin} "$R/out"
# mdCATH (CC BY 4.0): one HDF5 file per domain, sha256 as published by Hugging Face
while IFS=$'\t' read -r f size sha; do
    [ -f "$D/mdcath/$f" ] || curl -sS -L -o "$D/mdcath/$f.part" "https://huggingface.co/datasets/compsciencelab/mdCATH/resolve/main/data/$f" \
        && mv -f "$D/mdcath/$f.part" "$D/mdcath/$f" 2>/dev/null || true
    echo "$sha  $D/mdcath/$f" | sha256sum -c - | tail -1
done < "$R/data/mdcath-pick.tsv"
# RCSB PDB entries and PubChem 3D conformers
for f in $(awk '{print $2}' "$R/data/pdb.sha256"); do [ -f "$D/pdb/$f" ] || curl -sS -f -o "$D/pdb/$f" "https://files.rcsb.org/download/$f"; done
declare -A cid=([selegiline.sdf]=26757 [dexmethylphenidate.sdf]=154101 [lurasidone.sdf]=213046 [ziprasidone.sdf]=60854 [risperidone.sdf]=5073)
for f in $(awk '{print $2}' "$R/data/drugs.sha256"); do
    [ -f "$D/drugs/$f" ] && continue
    if [ -n "${cid[$f]:-}" ]; then curl -sS -f -o "$D/drugs/$f" "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid[$f]}/record/SDF?record_type=3d"
    else curl -sS -f -o "$D/drugs/$f" "https://files.rcsb.org/download/$f"; fi
done
(cd "$D/pdb" && sha256sum -c "$R/data/pdb.sha256" | grep -v ': OK$' || true)
(cd "$D/drugs" && sha256sum -c "$R/data/drugs.sha256" | grep -v ': OK$' || true)
# chignolin folding runs (Zenodo 6349893, CC BY 4.0): only the .gro trajectories and 5awl.pdb, by range requests
gcc -std=gnu11 -O2 -Wall -Wextra -Werror "$R/tools/zipget.c" -o "$R/out/zipget" -lcurl -lz
U='https://zenodo.org/records/6349893/files/Chignolin.zip?download=1'
[ -f "$D/chignolin/5.gro" ] || "$R/out/zipget" "$U" .gro "$D/chignolin"
[ -f "$D/chignolin/5awl.pdb" ] || "$R/out/zipget" "$U" 5awl.pdb "$D/chignolin"
(cd "$D/chignolin" && sha256sum -c "$R/data/chignolin.sha256" | grep -v ': OK$' || true)
echo "fetch-proteins: done; NVMe free $(df -B1 "$D" | tail -1 | awk '{print $4}') B"
