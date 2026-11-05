#!/bin/bash
# lurasidone-md.sh: lurasidone in water with AmberTools, on the host (G0 build).
#   1. an SDF (default: PubChem CID 213046 3D conformer, basic piperazine N protonated by
#      tools/sdf-protonate.c, +1); STATE, SDF, CHARGE, PROTONATE select another state or start
#   2. antechamber: GAFF2 atom types, AM1-BCC charges (sqm), residue LUR; parmchk2: missing terms
#   3. tleap: TIP3P truncated octahedron, 12 A buffer, Cl- to neutralise
#   4. sander.MPI (8 ranks): minimise, heat to 310 K with the drug restrained, 100 ps NpT, then
#      production with frames every 0.2 ps (ASCII mdcrd)
#   5. tools/mdcrd2pdb.c: the drug's frames as a multi-model PDB for amview
# Output: $W (default ~/.cache/intel-phi-ambermd/lurasidone) and logs/lurasidone/.
set -euo pipefail
R=$(cd "$(dirname "$0")/.." && pwd)
S=$R/input/dl/ambertools26/src/ambertools26_src
B=${B:-$HOME/.cache/intel-phi-ambermd/g0}
STATE=${STATE:-plus1-pubchem}
SDF=${SDF:-$R/input/dl/proteins/drugs/lurasidone.sdf}
CHARGE=${CHARGE:-1}
PROTONATE=${PROTONATE:-1}   # 1: add the proton (tools/sdf-protonate.c); 0: the SDF already has the intended state
W=${W:-$HOME/.cache/intel-phi-ambermd/lurasidone}
[ "$STATE" = plus1-pubchem ] || W=$HOME/.cache/intel-phi-ambermd/lurasidone-$STATE
PROD_PS=${PROD_PS:-200}
RANKS=${RANKS:-8}
mkdir -p "$W" "$R/logs/lurasidone" "$R/out"
cd "$W"

# an AMBERHOME for antechamber: its helper programs and data files by name
H=$W/amberhome
mkdir -p "$H/bin"
for p in antechamber parmchk2 am1bcc atomtype bondtype prepgen espgen residuegen; do ln -sf "$B/build/AmberTools/src/antechamber/$p" "$H/bin/$p"; done
ln -sf "$B/build/AmberTools/src/sqm/sqm" "$H/bin/sqm"
ln -sfn "$S/dat" "$H/dat"
export AMBERHOME=$H

gcc -std=gnu11 -O2 -Wall -Wextra -Werror "$R/tools/sdf-protonate.c" -o "$R/out/sdf-protonate" -lm
gcc -std=gnu11 -O2 -Wall -Wextra -Werror "$R/tools/mdcrd2pdb.c" -o "$R/out/mdcrd2pdb"
if [ "$PROTONATE" = 1 ]; then "$R/out/sdf-protonate" "$SDF" lurH.sdf; else cp "$SDF" lurH.sdf; fi

if [ ! -s lur.frcmod ]; then
    /usr/bin/time -f 'antechamber %e s' "$H/bin/antechamber" -i lurH.sdf -fi mdl -o lur.mol2 -fo mol2 \
        -c bcc -nc "$CHARGE" -at gaff2 -rn LUR -pf y > antechamber.log 2>&1
    "$H/bin/parmchk2" -i lur.mol2 -f mol2 -o lur.frcmod -s gaff2
fi
grep -c '^ATTN' lur.frcmod | sed 's/^/parmchk2 ATTN lines (guessed parameters needing review): /' || true
awk '/@<TRIPOS>ATOM/{f=1; next} /@<TRIPOS>/{f=0} f{q+=$9; n++} END{printf "mol2: %d atoms, total charge %.4f\n", n, q}' lur.mol2

cat > leap.in <<EOF
source leaprc.gaff2
source leaprc.water.tip3p
LUR = loadmol2 lur.mol2
loadamberparams lur.frcmod
check LUR
saveamberparm LUR lur_vac.prmtop lur_vac.inpcrd
savepdb LUR lur.pdb
sys = copy LUR
solvateoct sys TIP3PBOX 12.0
addions sys Cl- 0
charge sys
saveamberparm sys lur_wat.prmtop lur_wat.inpcrd
savepdb sys lur_wat.pdb
quit
EOF
"$B/build/AmberTools/src/leap/src/leap/teLeap" -I "$S/dat/leap/prep" -I "$S/dat/leap/lib" -I "$S/dat/leap/parm" \
    -I "$S/dat/leap/cmd" -f leap.in > leap.log 2>&1
grep -E 'Added [0-9]+ residues|Total unperturbed charge|Exiting LEaP' leap.log
sed -n 7p lur_wat.prmtop

cat > min.in <<'EOF'
minimise, drug restrained
 &cntrl imin=1, maxcyc=2000, ncyc=500, cut=10.0, ntb=1, ntr=1, restraint_wt=10.0, restraintmask=':LUR', ntpr=500, /
EOF
cat > heat.in <<'EOF'
heat to 310 K, 20 ps, drug restrained
 &cntrl imin=0, irest=0, ntx=1, nstlim=10000, dt=0.002, cut=10.0, ntb=1, ntc=2, ntf=2, ntt=3, gamma_ln=2.0,
  tempi=0.0, temp0=310.0, ig=-1, ntr=1, restraint_wt=5.0, restraintmask=':LUR', ntpr=1000, ntwx=0, /
EOF
cat > npt.in <<'EOF'
100 ps NpT
 &cntrl imin=0, irest=1, ntx=5, nstlim=50000, dt=0.002, cut=10.0, ntb=2, ntp=1, barostat=2, pres0=1.0,
  ntc=2, ntf=2, ntt=3, gamma_ln=2.0, temp0=310.0, ig=-1, ntpr=5000, ntwx=0, /
EOF
steps=$((PROD_PS * 500))
cat > prod.in <<EOF
production, frames every 0.2 ps, unwrapped (iwrap=0) so the drug stays whole
 &cntrl imin=0, irest=1, ntx=5, nstlim=$steps, dt=0.002, cut=10.0, ntb=2, ntp=1, barostat=2, pres0=1.0,
  ntc=2, ntf=2, ntt=3, gamma_ln=2.0, temp0=310.0, ig=-1, ntpr=5000, ntwx=100, ioutfm=0, iwrap=0, /
EOF
M() { mpirun -np "$RANKS" --map-by core --bind-to core "$B/build-mpi/AmberTools/src/sander/sander.MPI" "$@"; }
export AMBERHOME=$S
M -O -i min.in  -p lur_wat.prmtop -c lur_wat.inpcrd -ref lur_wat.inpcrd -r min.rst  -o min.out
M -O -i heat.in -p lur_wat.prmtop -c min.rst -ref min.rst -r heat.rst -o heat.out
M -O -i npt.in  -p lur_wat.prmtop -c heat.rst -r npt.rst -o npt.out
t0=$(date +%s.%N)
M -O -i prod.in -p lur_wat.prmtop -c npt.rst -r prod.rst -o prod.out -x prod.mdcrd
awk -v a="$t0" -v b="$(date +%s.%N)" -v p="$PROD_PS" 'BEGIN{printf "production: %d ps in %.1f s = %.2f ns/day\n", p, b-a, p*1e-3/((b-a)/86400)}'
nsolute=$(grep -c '^HETATM\|^ATOM' lur.pdb)
"$R/out/mdcrd2pdb" lur.pdb "$(sed -n 7p lur_wat.prmtop | awk '{print $1}')" "$nsolute" prod.mdcrd lurasidone-md.pdb
mkdir -p "$R/logs/lurasidone/$STATE"; cp min.out heat.out npt.out prod.out antechamber.log leap.log lur.frcmod "$R/logs/lurasidone/$STATE/"
echo "trajectory: $W/lurasidone-md.pdb"
