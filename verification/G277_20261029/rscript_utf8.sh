#!/bin/sh
set -eu
: "${ROUTE_B2_RSCRIPT_REAL:?Set ROUTE_B2_RSCRIPT_REAL to the pinned Rscript executable}"
export LANG=C.UTF-8 LC_ALL=C.UTF-8 LC_CTYPE=C.UTF-8
exec "$ROUTE_B2_RSCRIPT_REAL" "$@"
