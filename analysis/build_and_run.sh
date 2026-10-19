#!/usr/bin/env bash
# Build and run the offline single-fault validation against the OSIDH-LD
# reference eval_action().  Pure C, portable field arithmetic (256-bit),
# runs on any x86-64 PC — no ChipWhisperer needed.
set -e
cd "$(dirname "$0")"
REF=osidh-ld-ref

INC="-I$REF/gf_ref/ld256/include \
     -I$REF/ec_ref/include \
     -I$REF/protocols_ref/include \
     -I$REF/precomp_ref/ld256/include \
     -I$REF/common_generic/include"

SRC="$REF/gf_ref/ld256/fp.c $REF/gf_ref/ld256/fp2.c $REF/gf_ref/ld256/fp_p269.c \
     $REF/ec_ref/lvlx/ec.c $REF/ec_ref/lvlx/kps.c $REF/ec_ref/lvlx/poly-mul.c \
     $REF/ec_ref/lvlx/poly-redc.c $REF/ec_ref/lvlx/xeval.c $REF/ec_ref/lvlx/xisog.c \
     $REF/protocols_ref/lvlx/protocols.c \
     $REF/precomp_ref/ld256/protocol_constants.c $REF/precomp_ref/ld256/torsion_constants.c \
     $REF/common_generic/randombytes_system.c \
     fault_oracle_check.c"

echo "[build] gcc -O2 -DRADIX_64 ..."
gcc -O2 -DRADIX_64 $INC $SRC -o fault_oracle_check 2>/dev/null
echo "[run]"
./fault_oracle_check
