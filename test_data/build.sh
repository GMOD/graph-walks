#!/usr/bin/env bash
# Rebuilds the walk-file fixtures from their GFAs with gfa-to-tabix 4.0, in
# walk format 7. Usage: test_data/build.sh <gfa-to-tabix binary>
set -euo pipefail
cd "$(dirname "$0")"
B=${1:?the gfa-to-tabix binary}
build() {
  local name=$1 gfa=$2; shift 2
  "$B" "$gfa" --walks -o "$name.tmp" "$@" > /dev/null 2>&1
  local out
  for out in "$name".tmp.*.walks.gz; do
    local sample=${out#"$name".tmp.}; sample=${sample%.walks.gz}
    mv "$out" "$name.$sample.walks.gz"
    mv "$out.idx" "$name.$sample.walks.gz.idx"
  done
}
build header header.GRCh38.gfa --refs GRCh38 --tile 1000
build maxnode maxnode.GRCh38.gfa --refs GRCh38 --tile 1000
build revisit revisit.GRCh38.gfa --refs GRCh38 --tile 1000
build farpass farpass.GRCh38.gfa --refs GRCh38 --tile 1000
build follow follow.GRCh38.gfa --refs GRCh38 --tile 1000
build rejoin follow.GRCh38.gfa --refs GRCh38 --tile 10000
build records records.GRCh38.gfa --refs GRCh38 --tile 1000
build skip skip.GRCh38.gfa --refs GRCh38 --tile 10
build span span.gfa --refs GRCh38 --tile 20
build attach attach.GRCh38.gfa --refs GRCh38 --tile 20
build unabridged unabridged.R.gfa --refs R --tile 10
build abridged unabridged.R.gfa --refs R --tile 10 --abridge-bytes 1
build interleaved interleaved.R.gfa --refs R --tile 10 --stay 0 --abridge-steps 3
build ownends ownends.R.gfa --refs R --tile 10 --abridge-bytes 1
build contigends contigends.R.gfa --refs R --tile 10 --abridge-bytes 1
build hprc-chrM hprc-chrM.GRCh38.gfa --refs GRCh38 --sequences
build jumpback jumpback.GRCh38.gfa --refs GRCh38 --tile 20
build loop loop.R.gfa --refs R --tile 10
build ownzero ownzero.R.gfa --refs R --tile 10
build excursion excursion.R.gfa --refs R --tile 64
build hole hole.R.gfa --refs R --tile 16
build holetile hole.R.gfa --refs R --tile 64
build copies copies.R.gfa --refs R --tile 1000
build copyends copyends.R.gfa --refs R --tile 1000 --abridge-bytes 1
build zerostep zerostep.R.gfa --refs R --tile 8
