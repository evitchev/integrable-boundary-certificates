#!/bin/bash
# Attested wrapper for lab/modp_scan.py (Codex rounds 20-21): validates and prints the
# execution provenance BEFORE computing, and ABORTS on any disagreement, so that the chunk
# log carries an attestation block the certificate can enforce:
#   ATTEST scanner_sha256=<64hex>
#   ATTEST matrix=<file> matrix_sha256=<64hex>
#   ATTEST sidecar_matrix_sha256=<64hex>         (must equal matrix_sha256)
#   ATTEST host=<name> python=<ver> flint=<ver>
#   ATTEST invocation=modp_scan.py W SOL P T0 T1 started=<UTC>
#   ATTEST git=<40hex> dirty=<n>                  (a 40-hex commit and a tracked-clean tree are required)
# The scanner runs isolated and bytecode-free (python -I -B, empty environment) and the
# wrapper refuses untracked/ignored entries or bytecode under lab/ (Codex round 22).
# Usage: modp_scan_attested.sh W SOL P [T0 T1]
set -euo pipefail
W=$1; SOL=$2; P=$3; T0=${4:-}; T1=${5:-}
HERE=$(cd "$(dirname "$0")/.." && pwd)
SCANNER=$HERE/lab/modp_scan.py
MATRIX=$HERE/wolfram/lab/cyl_strata_w${W}_sol${SOL}.m
SIDECAR=$HERE/wolfram/lab/cyl_strata_w${W}_sol${SOL}.json
PY=${MODP_PYTHON:-$HOME/python312-lab/bin/python3}
fail() { echo "ATTEST-ABORT $*" >&2; exit 3; }
[ -f "$SCANNER" ] || fail "scanner missing: $SCANNER"
[ -f "$MATRIX" ] || fail "matrix missing: $MATRIX"
[ -f "$SIDECAR" ] || fail "sidecar missing: $SIDECAR"
ssha=$(sha256sum "$SCANNER" | cut -d' ' -f1); [ ${#ssha} -eq 64 ] || fail "scanner hash"
msha=$(sha256sum "$MATRIX" | cut -d' ' -f1); [ ${#msha} -eq 64 ] || fail "matrix hash"
side=$(grep -o '"matrix_sha256": "[0-9a-f]*' "$SIDECAR" | cut -d'"' -f4); [ ${#side} -eq 64 ] || fail "sidecar hash field"
[ "$msha" = "$side" ] || fail "matrix hash $msha != sidecar $side"
pyv=$("$PY" -I -B -c 'import sys; print(sys.version.split()[0])') || fail "interpreter"
flv=$("$PY" -I -B -c 'import flint; print(flint.__version__)') || fail "python-flint missing"
# Codex round 22: an untracked or ignored lab/fractions.py (or bytecode) would shadow a
# dependency of the scanner; refuse any such entry and run Python isolated (-I: script
# directory off sys.path, PYTHON* ignored) and bytecode-free (-B).
stat_out=$(git -C "$HERE" status --porcelain --untracked-files=all -- lab) || fail "git status under lab/"
untr=$(printf '%s\n' "$stat_out" | grep -c '^??' || true)      # grep -c returns 1 on zero matches; git status itself is checked above
ign_out=$(git -C "$HERE" ls-files --others --ignored --exclude-standard -- lab) || fail "git ls-files under lab/"
ign=$(printf '%s\n' "$ign_out" | grep -c . || true)
byc=$(find "$HERE/lab" \( -name '__pycache__' -o -name '*.pyc' -o -name '*.pyo' \) 2>/dev/null | wc -l)
[ "$untr" -eq 0 ] || fail "untracked entries under lab/ ($untr)"
[ "$ign" -eq 0 ] || fail "ignored entries under lab/ ($ign)"
[ "$byc" -eq 0 ] || fail "bytecode under lab/ ($byc)"
commit=$(git -C "$HERE" rev-parse HEAD 2>/dev/null) || fail "not a git repository: $HERE"
[[ "$commit" =~ ^[0-9a-f]{40}$ ]] || fail "commit not 40-hex: $commit"
dirty=$(git -C "$HERE" status --porcelain --untracked-files=no 2>/dev/null | wc -l) || fail "git status"
[ "$dirty" -eq 0 ] || fail "tracked tree dirty ($dirty entries)"
echo "ATTEST scanner_sha256=$ssha"
echo "ATTEST matrix=$(basename "$MATRIX") matrix_sha256=$msha"
echo "ATTEST sidecar_matrix_sha256=$side"
echo "ATTEST host=$(hostname) python=$pyv flint=$flv"
echo "ATTEST invocation=modp_scan.py $W $SOL $P $T0 $T1 started=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "ATTEST git=$commit dirty=$dirty"
exec env -i HOME="$HOME" PATH="$PATH" "$PY" -I -B -u "$SCANNER" "$W" "$SOL" "$P" $T0 $T1
