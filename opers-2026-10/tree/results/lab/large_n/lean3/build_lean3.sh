#!/usr/bin/env bash
# LEAN3 replay: build both certificates, print the axioms of every theorem,
# build four mutation copies (each must FAIL), and check that the checker
# rejects a false proof.  Writes BUILD_LOG.txt next to this script.
# Exit 0 only if every check came out as required.
set -u
cd "$(dirname "$0")"
W="${LEAN_MATHLIB:-$HOME/projects/lean-benchmarks/lean-mathlib}"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
fail=0
name_at() {  # the theorem/lemma enclosing line $2 of file $1
  awk -v L="$2" 'NR<=L && /^(theorem|lemma)/{t=$2} END{print t}' "$1"
}
mutate() {  # $1 file, $2 old, $3 new, $4 output name
  python3 - "$1" "$2" "$3" "$T/$4.lean" <<'EOF'
import sys
f, a, b, out = sys.argv[1:]
s = open(f).read()
assert s.count(a) == 1, (out, s.count(a))
open(out, "w").write(s.replace(a, b))
EOF
}
{
echo "LEAN3 build log (cognitive, $(date -u +%Y-%m-%dT%H:%MZ))"
echo "toolchain: $("$W" --version)"
echo "Mathlib: v4.24.0, prebuilt oleans via the lean-mathlib wrapper (baked LEAN_PATH; no Lake project)"
echo
for f in Theorem3Intertwining Theorem3FrozenLattice; do
  start=$(date +%s); out=$("$W" $f.lean 2>&1); rc=$?; end=$(date +%s)
  echo "== lean $f.lean  sha256 $(sha256sum $f.lean | cut -c1-64)"
  echo "exit=$rc  seconds=$((end-start))  output_lines=$(printf '%s' "$out" | grep -c .)"
  [ -n "$out" ] && printf '%s\n' "$out"
  [ $rc -eq 0 ] && [ -z "$out" ] || fail=1
  if grep -nE '\bsorry\b|\badmit\b|^axiom |native_decide' $f.lean | grep -v 'no sorry'; then fail=1; fi
done
echo
echo "== #print axioms, every theorem (a copy with #print axioms appended)"
for f in Theorem3Intertwining Theorem3FrozenLattice; do
  cp $f.lean "$T/${f}_axioms.lean"; ns=$(grep -oP '^namespace \K.*' $f.lean)
  for t in $(grep -oP '^theorem \K\S+' $f.lean); do echo "#print axioms $ns.$t" >> "$T/${f}_axioms.lean"; done
  out=$("$W" "$T/${f}_axioms.lean" 2>&1); rc=$?
  echo "$f: exit=$rc, theorems=$(grep -cP '^theorem ' $f.lean)"
  flat=$(printf '%s\n' "$out" | tr '\n' ' ' | sed 's/ \+/ /g; s/'"'"'IntegrableBoundary/\n'"'"'IntegrableBoundary/g' | sed 's/^ *//' | grep .)
  printf '%s\n' "$flat"
  n_ok=$(printf '%s\n' "$flat" | grep -c 'depends on axioms: \[propext, Classical.choice, Quot.sound\] *$')
  [ $rc -eq 0 ] && [ "$n_ok" -eq "$(grep -cP '^theorem ' $f.lean)" ] || fail=1
done
echo
echo "== mutation copies (each must FAIL to build)"
mutate Theorem3Intertwining.lean "theorem m_sub_six (μ : ℂ) : m μ = (μ - 1) * (μ - 4) * m (μ - 6)" \
  "theorem m_sub_six (μ : ℂ) : m μ = (μ - 1) * (μ - 5) * m (μ - 6)" M1_recurrence_mu5
mutate Theorem3Intertwining.lean "PT P c ν + xpow 9 (fun μ => theta c μ + 2 * c μ) ν - E * xpow 6 c ν" \
  "PT P c ν + xpow 9 (fun μ => theta c μ + 3 * c μ) ν - E * xpow 6 c ν" M2_L4_theta_plus_3
mutate Theorem3FrozenLattice.lean "noncomputable def Ind (l0 l1 μ : ℂ) : ℂ := P l0 l1 (s - μ) * (μ - 1) * (μ - 4)" \
  "noncomputable def Ind (l0 l1 μ : ℂ) : ℂ := P l0 l1 (s - μ) * (μ - 1) * (μ - 5)" M3_indicial_mu5
mutate Theorem3FrozenLattice.lean "theorem m_eq_zero_iff (μ : ℂ) : m μ = 0 ↔ ∃ n : ℕ, μ = 1 + 3 * n" \
  "theorem m_eq_zero_iff (μ : ℂ) : m μ = 0 ↔ ∃ n : ℕ, μ = 2 + 3 * n" M4_lattice_2_plus_3Z
for m in "M1_recurrence_mu5:(mu-1)(mu-5) in m_sub_six" "M2_L4_theta_plus_3:x^9(theta+3) in L4" \
         "M3_indicial_mu5:(mu-1)(mu-5) in Ind" "M4_lattice_2_plus_3Z:2+3n in m_eq_zero_iff"; do
  n=${m%%:*}; d=${m#*:}
  "$W" "$T/$n.lean" > "$T/$n.log" 2>&1; rc=$?
  broken=$(grep -oP '^\S+\.lean:\d+' "$T/$n.log" | sort -u | while IFS=: read -r f l; do name_at "$f" "$l"; done | sort -u | tr '\n' ' ')
  echo "$n [$d]: exit=$rc; failing: $broken"
  [ $rc -ne 0 ] || fail=1
done
echo
echo "== checker negative control (a false proof must be rejected)"
printf 'import Mathlib.Tactic\ntheorem bogus : (2:ℕ) + 2 = 5 := by norm_num\n' > "$T/NEG_control.lean"
"$W" "$T/NEG_control.lean" > /dev/null 2>&1; rc=$?
echo "theorem bogus : (2:ℕ) + 2 = 5 := by norm_num  ->  exit=$rc"
[ $rc -ne 0 ] || fail=1
echo
[ $fail -eq 0 ] && echo "LEAN3: ALL CHECKS AS REQUIRED" || echo "LEAN3: A CHECK FAILED"
} > BUILD_LOG.txt 2>&1
exit $fail
