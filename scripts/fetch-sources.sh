#!/bin/bash
# Fetch the AmberTools26 source and check it. input/dl/ is ignored by git.
# The file is byte-identical to the AmberTools26 release by the md5 on ambermd.org/GetAmber.php.
set -euo pipefail

url=https://ambermd.org/downloads/ambertools26_rc7.tar.bz2
sha=5d46eef3c2bb7d5bf9e8c0c38add34406ea67e3f0e4097ac9d11d8a544538c9c
md5=370620f5ddb646e552a182ce8625b78e

root=$(cd "$(dirname "$0")/.." && pwd)
dir=$root/input/dl/ambertools26
file=$dir/ambertools26_rc7.tar.bz2

mkdir -p "$dir"
[ -f "$file" ] || curl -sS -L -o "$file" "$url"
echo "$sha  $file" | sha256sum -c -
echo "$md5  $file" | md5sum -c -
if [ ! -d "$dir/src/ambertools26_src" ]; then
    mkdir -p "$dir/src"
    tar xjf "$file" -C "$dir/src"
fi
echo "source: $dir/src/ambertools26_src"
