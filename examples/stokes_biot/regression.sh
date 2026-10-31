#!/bin/bash
#
# Regression test for stokes_br.ex: runs a fixed set of cases on levels
# 1, ..., 4 and compares the output with the references in reference/.
#
# Usage: ./regression.sh            compare against reference/
#        ./regression.sh --update   regenerate reference/
#
# Timing columns (t_asm, t_gmres, t_umf) and |xg-xu| are dropped.  ndof must match
# exactly, GMRES iterations to ITS_TOL (UA-AMG varies by up to two iterations
# from run to run), errors to RTOL, rates to RATE_TOL, and max|divT| must not
# exceed DIV_FAC times the reference (it is at the GMRES tolerance level).
# Text lines must match exactly.
#

EXDIR=$(cd "$(dirname "$0")" && pwd)
EXE=$EXDIR/stokes_br.ex
REFDIR=$EXDIR/reference
LEVELS=4
RTOL=1e-3
RATE_TOL=0.02
ITS_TOL=2
DIV_FAC=10

# name:dim:arguments
CASES=(
  "2d_default:2:"
  "2d_test1:2:test=1"
  "2d_test2_mu1e-6:2:test=2 mu=1e-6"
  "2d_full_umfpack:2:full solver=umfpack"
  "2d_ua_amli:2:amg=ua"
  "3d_default:3:"
  "3d_test2_mu1e-6:3:test=2 mu=1e-6"
  "3d_ua_amli:3:amg=ua"
)

update=0
[ "$1" = "--update" ] && update=1

if [ ! -x "$EXE" ]; then
  echo "$EXE not found; run make first"
  exit 1
fi

WORK=$(mktemp -d)
trap "rm -rf $WORK" EXIT
mkdir -p "$WORK/output" "$REFDIR"

# input.dat and stokes_br.dat (defaults) for one run; the program reads both
# from the current directory and the case arguments override stokes_br.dat
write_input() {
  cat > "$WORK/input.dat" <<EOF
print_level	= 0
output_dir	= ./output/
spatial_dim       = $1
refinement_type   = 11
refinement_levels = $LEVELS
boundary_codes    = 1
nquad	= 5
AMG_type       = C
AMG_cycle_type = A
EOF
  cat > "$WORK/stokes_br.dat" <<EOF
viscosity    = 1
test_problem = 0
rhs          = pr
system       = condensed
solver       = gmres
EOF
}

# keep text lines; data lines without timings and |xg-xu|
strip_timings() {
  awk 'NF == 14 && $1 ~ /^[0-9]+$/ {
         print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10; next }
       { print }'
}

# compare two stripped outputs; prints the differing lines, exit 1 if any
compare() {
  awk -v rtol=$RTOL -v ratetol=$RATE_TOL -v itstol=$ITS_TOL -v divfac=$DIV_FAC '
    function abs(x) { return x < 0 ? -x : x }
    function bad(msg) { printf "    line %d: %s\n      ref: %s\n      new: %s\n", FNR, msg, r, $0; nbad++ }
    NR == FNR { ref[FNR] = $0; nref = FNR; next }
    {
      r = ref[FNR]
      if (!(FNR in ref)) { bad("extra line"); next }
      n = split(r, a, " ")
      if (NF == 10 && n == 10 && $1 ~ /^[0-9]+$/) {
        if ($1 != a[1] || $2 != a[2]) { bad("level/ndof"); next }
        for (k = 3; k <= 7; k += 2)
          if (abs($k - a[k]) > rtol * abs(a[k])) { bad("error, column " k); next }
        for (k = 4; k <= 8; k += 2)
          if (abs($k - a[k]) > ratetol) { bad("rate, column " k); next }
        if (abs($9) > divfac * abs(a[9])) { bad("max|divT|"); next }
        if (abs($10 - a[10]) > itstol) { bad("iterations"); next }
      } else if ($0 != r) bad("text")
    }
    END {
      if (FNR < nref) { printf "    missing lines after %d\n", FNR; nbad++ }
      exit nbad > 0
    }' "$1" "$2"
}

PASS=0
FAIL=0
for c in "${CASES[@]}"; do
  IFS=: read -r name dim args <<< "$c"
  write_input "$dim"
  (cd "$WORK" && /usr/bin/time -f "%e %M" -o "$WORK/time" $EXE $args 2>/dev/null) \
    | strip_timings > "$WORK/$name.out"
  # wall time and peak memory of the run
  tm=$(awk 'END { printf "%7.2f s %8.1f MB", $1, $2 / 1024 }' "$WORK/time")
  if [ $update -eq 1 ]; then
    mv "$WORK/$name.out" "$REFDIR/$name.ref"
    printf "updated  %-18s %s\n" "$name" "$tm"
  elif [ ! -f "$REFDIR/$name.ref" ]; then
    printf "FAIL     %-18s (no reference)\n" "$name"
    FAIL=$((FAIL + 1))
  elif compare "$REFDIR/$name.ref" "$WORK/$name.out"; then
    printf "PASS     %-18s %s\n" "$name" "$tm"
    PASS=$((PASS + 1))
  else
    printf "FAIL     %-18s %s\n" "$name" "$tm"
    FAIL=$((FAIL + 1))
  fi
done

[ $update -eq 1 ] && exit 0
echo "passed $PASS, failed $FAIL"
[ $FAIL -eq 0 ]
