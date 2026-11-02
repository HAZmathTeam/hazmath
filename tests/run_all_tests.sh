#!/bin/bash
#
# Rebuilds libhazmath and runs the tests:
#   1. rebuild:   make distclean; make config $HAZ_CONFIG; make install
#                 (in the source tree, as for normal use)
#   2. sanitizer: libhazmath built with AddressSanitizer and
#                 UndefinedBehaviorSanitizer (in the output directory), and the
#                 C tests in tests/sanitizer/ (memory errors, leaks, the HX
#                 apply contract, failed setups)
#   3. examples:  all examples except haznics, built and run against the
#                 rebuilt library, on a copy of examples/ (so the example
#                 output/ directories are not touched), with stokes_biot's
#                 regression.sh
#   4. python:    the haznicsx package in a virtual environment, and
#                 python/tests (needs a network connection for pip)
#
# Usage: tests/run_all_tests.sh [--no-rebuild] [--no-sanitizer]
#                               [--no-examples] [--no-python] [--long]
#   --long also runs "stokes_br.ex solver=compare amg=ua" (UMFPACK on the
#   3D level-5 problem: minutes, and it may run out of memory, which the
#   example reports and survives).
#
# Environment:
#   HAZ_CONFIG        options for "make config"
#                     (default: "suitesparse=yes lapack=yes shared=yes")
#   HAZ_TEST_DIR      output directory (default: a new directory in $TMPDIR
#                     or /tmp); the logs and summary.txt are there
#   HAZ_TEST_TIMEOUT  seconds per run (default 900)
#
# The exit status is 0 when everything passed, 1 otherwise.

set -u
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${HAZ_TEST_DIR:-$(mktemp -d "${TMPDIR:-/tmp}/hazmath-tests.XXXXXX")}
CONFIG=${HAZ_CONFIG:-"suitesparse=yes lapack=yes shared=yes"}
TMAX=${HAZ_TEST_TIMEOUT:-900}
mkdir -p "$OUT/logs"
SUMMARY="$OUT/summary.txt"
: > "$SUMMARY"

do_rebuild=1 do_sanitizer=1 do_examples=1 do_python=1 do_long=0
for a in "$@"; do
  case $a in
    --no-rebuild) do_rebuild=0 ;;
    --no-sanitizer) do_sanitizer=0 ;;
    --no-examples) do_examples=0 ;;
    --no-python) do_python=0 ;;
    --long) do_long=1 ;;
    *) echo "unknown option $a (see the header of $0)"; exit 2 ;;
  esac
done

TIMEOUT=""
if command -v timeout > /dev/null; then TIMEOUT="timeout $TMAX"
elif command -v gtimeout > /dev/null; then TIMEOUT="gtimeout $TMAX"; fi

nfail=0
record() { # name status seconds
  local verdict=PASS
  [ "$2" -ne 0 ] && verdict=FAIL && nfail=$((nfail+1))
  printf "%-4s %-34s exit=%-4s %5ss\n" "$verdict" "$1" "$2" "$3" | tee -a "$SUMMARY"
}
run() { # name dir command... (log in logs/name.log)
  local name=$1 dir=$2; shift 2
  local t0=$(date +%s)
  ( cd "$dir" && $TIMEOUT "$@" ) > "$OUT/logs/$name.log" 2>&1
  local st=$?
  record "$name" $st $(( $(date +%s) - t0 ))
}

echo "HAZmath tests: source $ROOT, output $OUT"

#------------------------------------------------------------------
# 1. rebuild libhazmath in the source tree
#------------------------------------------------------------------
if [ $do_rebuild = 1 ]; then
  run rebuild_distclean "$ROOT" make distclean
  run rebuild_config "$ROOT" make config $CONFIG
  run rebuild_install "$ROOT" make install
  if [ ! -e "$ROOT/lib/libhazmath.so" ] && [ ! -e "$ROOT/lib/libhazmath.a" ] \
       && [ ! -e "$ROOT/lib/libhazmath.dylib" ]; then
    echo "libhazmath was not built; see $OUT/logs/rebuild_*.log"
    exit 1
  fi
fi

#------------------------------------------------------------------
# 2. sanitizer build and tests/sanitizer/*.c
#------------------------------------------------------------------
if [ $do_sanitizer = 1 ]; then
  SAN="-fsanitize=address,undefined -fno-omit-frame-pointer"
  run sanitizer_build "$ROOT" sh -c "cmake -S '$ROOT' -B '$OUT/build-san' \
      -DCMAKE_BUILD_TYPE=Debug -DSHARED=OFF -DUSE_LAPACK=ON \
      -DHAZMATH_RELOCATABLE_INSTALL=ON -DCMAKE_INSTALL_PREFIX='$OUT/prefix-san' \
      -DCMAKE_C_FLAGS='$SAN' && cmake --build '$OUT/build-san' --parallel \
      && cmake --install '$OUT/build-san'"
  if [ -e "$OUT/prefix-san/lib/libhazmath.a" ]; then
    for src in "$ROOT"/tests/sanitizer/*.c; do
      t=$(basename "$src" .c)
      # the HX and RA factories are in src/haznics/helper.c
      run "sanitizer_build_$t" "$OUT" gcc -g -O0 $SAN \
          -I"$OUT/prefix-san/include/hazmath" -I"$ROOT/tests/sanitizer" \
          -o "$OUT/$t" "$src" "$ROOT/src/haznics/helper.c" \
          "$OUT/prefix-san/lib/libhazmath.a" -llapack -lblas -lgfortran -lm
      [ -x "$OUT/$t" ] && ASAN_OPTIONS=detect_leaks=1 \
          UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
          run "sanitizer_$t" "$OUT" "./$t"
    done
  fi
fi

#------------------------------------------------------------------
# 3. examples (no haznics), on a copy linked to the rebuilt library
#------------------------------------------------------------------
if [ $do_examples = 1 ]; then
  EXC="$OUT/ex"
  mkdir -p "$EXC"
  ln -sfn "$ROOT/lib" "$EXC/lib"
  ln -sfn "$ROOT/include" "$EXC/include"
  if command -v rsync > /dev/null; then
    rsync -a --exclude haznics --exclude '*.o' --exclude '*.ex' --exclude '*.so' \
      --exclude 'output/*' --include 'output/.gitignore' "$ROOT/examples/" "$EXC/examples/"
  else
    cp -R "$ROOT/examples" "$EXC/" && rm -rf "$EXC/examples/haznics"
  fi
  EX="$EXC/examples"

  for d in amg_classical amr_grids approximation basic_elliptic elliptic_p1 \
           heat_equation mg_geometric navierstokes solvers stokes stokes_biot \
           xd_1d_models liquidcrystals/PureElastic; do
    run "build_$(echo $d | tr / _)" "$EX/$d" make -B
  done
  run build_xd_1d_meshing "$EX/xd_1d_models" make -B SRCFILE=main_meshing_xd_1d

  run amg_classical_test_amg      "$EX/amg_classical" ./test_amg.ex ../solvers/A.dat
  run amg_classical_test_amg_UA   "$EX/amg_classical" ./test_amg.ex --cycle=UA ../solvers/A.dat
  run amg_classical_ilu_ichol     "$EX/amg_classical" ./test_ilu_ichol.ex ../solvers/A.dat
  for f in "$EX"/amr_grids/input/*.input; do
    run "amr_grids_$(basename "$f" .input)" "$EX/amr_grids" ./amr_grids.ex "input/$(basename "$f")"
  done
  run amr_grids_dgs_adaptive_2d   "$EX/amr_grids" ./test_dgs_adaptive.ex 2 3 1
  run amr_grids_dgs_adaptive_3d   "$EX/amr_grids" ./test_dgs_adaptive.ex 3 2 1
  run amr_grids_msh_roundtrip     "$EX/amr_grids" ./test_msh_roundtrip.ex
  run approximation_test_aaa      "$EX/approximation" ./test_aaa.ex
  run basic_elliptic              "$EX/basic_elliptic" ./basic_elliptic.ex
  run elliptic_p1                 "$EX/elliptic_p1" ./elliptic_p1.ex
  run heat_equation               "$EX/heat_equation" ./heat_equation.ex
  run mg_geometric                "$EX/mg_geometric" ./gmg.ex
  run navierstokes                "$EX/navierstokes" ./navierstokes.ex
  run solvers                     "$EX/solvers" ./solvers.ex A.dat b.dat
  run stokes                      "$EX/stokes" ./stokes.ex
  run liquidcrystals_PureElastic  "$EX/liquidcrystals/PureElastic" ./lcelastic.ex
  run xd_1d_meshing_2d            "$EX/xd_1d_models" python3 -u run_meshing.py -d 2
  run xd_1d_meshing_3d            "$EX/xd_1d_models" python3 -u run_meshing.py -d 3
  run xd_1d_solver                "$EX/xd_1d_models" python3 -u run_solver.py
  run stokes_biot_default         "$EX/stokes_biot" ./stokes_br.ex
  run stokes_biot_ua_amli         "$EX/stokes_biot" ./stokes_br.ex amg=ua
  [ $do_long = 1 ] && run stokes_biot_compare "$EX/stokes_biot" ./stokes_br.ex solver=compare amg=ua
  run stokes_biot_regression      "$EX/stokes_biot" ./regression.sh
fi

#------------------------------------------------------------------
# 4. python bindings
#------------------------------------------------------------------
if [ $do_python = 1 ]; then
  if command -v python3 > /dev/null; then
    # --system-site-packages: DOLFINx and petsc4py, when installed, are found
    # and their tests run; otherwise those tests are skipped
    run python_venv "$OUT" python3 -m venv --system-site-packages "$OUT/venv"
    PY="$OUT/venv/bin/python"
    if [ -x "$PY" ]; then
      # the tests need scipy >= 1.12 (gmres(rtol=...)); numpy is kept at
      # the version already there, which a system DOLFINx was built with
      run python_deps "$OUT" sh -c "'$PY' -c 'import numpy' 2>/dev/null || '$PY' -m pip install -q numpy; \
          '$PY' -c 'import pytest' 2>/dev/null || '$PY' -m pip install -q pytest; \
          '$PY' -c 'import scipy, sys; sys.exit(tuple(int(x) for x in scipy.__version__.split(\".\")[:2]) < (1, 12))' 2>/dev/null \
          || '$PY' -m pip install -q 'scipy>=1.12' \"numpy==\$('$PY' -c 'import numpy; print(numpy.__version__)')\""
      run python_install "$ROOT" "$PY" -m pip install --no-deps "$ROOT" \
          -C build-dir="$OUT/pybuild"
      run python_tests "$ROOT" "$PY" -m pytest -q -p no:cacheprovider python/tests
    fi
  else
    echo "python3 not found: python tests skipped" | tee -a "$SUMMARY"
  fi
fi

echo
echo "$(grep -c '^PASS' "$SUMMARY") passed, $nfail failed; summary $SUMMARY, logs $OUT/logs/"
[ $nfail -eq 0 ]
