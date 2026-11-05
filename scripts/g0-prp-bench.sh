#!/bin/bash
# G0: build the PrP benchmark box on the host, equilibrate it, then time sander on the host.
# Output goes to logs/g0/prp-bench.tsv. Interleaved rounds: serial, 8 ranks, 14 ranks.
set -euo pipefail
R=$(cd "$(dirname "$0")/.." && pwd)
S=$R/input/dl/ambertools26/src/ambertools26_src
B=${B:-$HOME/.cache/intel-phi-ambermd/g0}
W=$B/prp
SANDER=$B/build/AmberTools/src/sander/sander
SMPI=$B/build-mpi/AmberTools/src/sander/sander.MPI
LEAP=$B/build/AmberTools/src/leap/src/leap/teLeap
export AMBERHOME=$S
mkdir -p "$W" "$R/logs/g0"; cd "$W"

# 1QM0 holds residues 125-228. Hydrogens are stripped, teLeap rebuilds them.
awk '/^ATOM/ { e=substr($0,77,2); gsub(/ /,"",e); n=substr($0,13,4); gsub(/ /,"",n);
               if (e=="H" || n ~ /^[0-9]*H/) next;
               r=substr($0,23,4)+0; a=substr($0,18,3); if (r==179||r==214) a="CYX";
               print substr($0,1,17) a substr($0,21) } END{print "END"}' \
    "$R/input/dl/pdb/1QM0.pdb" > prp.pdb
sed "s|@TIP3PBOX@|$S/dat/leap/lib/tip3pbox.off|" "$R/systems/prp-1qm0/build.leap" > build.in
"$LEAP" -I "$S/dat/leap/prep" -I "$S/dat/leap/lib" -I "$S/dat/leap/parm" -I "$S/dat/leap/cmd" -f build.in > leap.log 2>&1
grep -E 'Added [0-9]+ residues|Total unperturbed charge|Exiting' leap.log
sed -n 7p prp.prmtop | cut -c1-16

M() { mpirun -np 14 --use-hwthread-cpus --bind-to hwthread "$SMPI" "$@"; }
if [ ! -f npt.rst ]; then
    M -O -i "$R/systems/prp-1qm0/min.in"  -p prp.prmtop -c prp.inpcrd -ref prp.inpcrd -r min.rst  -o min.out  -inf min.inf
    M -O -i "$R/systems/prp-1qm0/heat.in" -p prp.prmtop -c min.rst -ref min.rst -r heat.rst -o heat.out -inf heat.inf
    M -O -i "$R/systems/prp-1qm0/npt.in"  -p prp.prmtop -c heat.rst   -r npt.rst  -o npt.out  -inf npt.inf
fi

out=$R/logs/g0/prp-bench.tsv
printf 'round\tranks\twall_s\tsteps_per_s\tns_per_day\tload1\n' > "$out"
stamp() { date +%s.%N; }
one() { # round ranks
    local r=$1 n=$2 s e
    read -r l1 _ < /proc/loadavg
    s=$(stamp)
    if [ "$n" = 1 ]; then
        "$SANDER" -O -i "$R/systems/prp-1qm0/bench.in" -p prp.prmtop -c npt.rst -r b.rst -o b$n.$r.out -inf b.inf
    elif [ "$n" = 8 ]; then
        mpirun -np 8 --map-by core --bind-to core "$SMPI" -O -i "$R/systems/prp-1qm0/bench.in" -p prp.prmtop -c npt.rst -r b.rst -o b$n.$r.out -inf b.inf
    else
        mpirun -np "$n" --use-hwthread-cpus --bind-to hwthread "$SMPI" -O -i "$R/systems/prp-1qm0/bench.in" -p prp.prmtop -c npt.rst -r b.rst -o b$n.$r.out -inf b.inf
    fi
    e=$(stamp)
    awk -v r="$r" -v n="$n" -v s="$s" -v e="$e" -v l="$l1" 'BEGIN{w=e-s; st=1000/w; printf "%d\t%d\t%.2f\t%.3f\t%.3f\t%s\n", r,n,w,st,st/5.787037,l}' >> "$out"
}
for r in 1 2 3; do one "$r" 1; one "$r" 8; one "$r" 14; done
cat "$out"
