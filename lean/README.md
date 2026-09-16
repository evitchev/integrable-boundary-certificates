# Kernel-checked certificates (Lean 4 + Mathlib)

Formal, machine-checked proofs of the abstract skeletons that the
draft's arguments lean on, contributed by the cognitive agent
(cross-project assessment item C). Checked with Lean 4.24 / Mathlib
via a pinned local toolchain; `lean-toolchain` in this directory pins
`leanprover/lean4:v4.24.0` for elan-based reproduction (Mathlib of the
matching release is required on LEAN_PATH; note the LIMITATION that
Mathlib itself is supplied by a machine-local wrapper, not yet by a
repository-pinned Lake manifest); run `python3
code/lean_check.py` (or the
certifier's `lean_certificates` entry under `--only`/`--full`).

## Contents

`Bootstrap.lean` —

- `IntegrableBoundary.bootstrap` / `bootstrap_odd`: **draft Lemma 1**
  (the bootstrap lemma), in an arbitrary Lie algebra graded by
  submodules with bracket-additive grading: two odd-spin charges
  commuting with the distinguished element commute with each other
  whenever their spin-sum grade is vacant. The Lean statement isolates
  exactly the draft's hypotheses; notably the distinguished element's
  own grade never enters.
- `IntegrableBoundary.vanishes_identically` /
  `commutator_entry_vanishes`: **R3's degree-bound interpolation
  argument** — a polynomial of degree ≤ 14 over ℚ vanishing at 20
  distinct points is identically zero (generic form over any integral
  domain).

`Specialization.lean` —

- `IntegrableBoundary.ratfunc_vanishes_identically`: **rational-section
  vanishing** — a rational function with bounded numerator degree,
  vanishing at more good points (denominator nonzero) than that bound,
  is the zero function.
- `IntegrableBoundary.kernel_specialization`: **the one-point-vacancy /
  specialization lemma** — the inference under every GF(p)/ℚ(t)
  high-spin certificate: if the entrywise evaluation `M(t₀)` of a
  polynomial matrix has trivial kernel at one point, then `M` has
  trivial kernel over the polynomial ring (hence over ℚ(t) after the
  engines' denominator clearing). Proof by divisibility descent in
  `(X − C t₀)`: a kernel vector specializes into the trivial kernel,
  so every power divides every entry, which the degree bound forbids.

## Honest scope

Lean certifies the abstract logical skeletons. The exact-arithmetic
engines certify the instances: that the physical charge spaces `Q^σ`
with the OPE bracket satisfy the grading axiom, and that the vacancy
and degree-bound hypotheses hold in the verified ranges. Neither
replaces the other; the free-boson realization itself is not
formalized. Remaining future candidate: the density-zero
exponent-collision lemma of the infinitude program (the v2 pass
delivered the ℚ(t)-section and specialization lemmas).

## Toolchain

The checker (`code/lean_check.py`) uses the machine-local pinned
wrapper `~/projects/lean-benchmarks/lean-mathlib` (Lean 4.24, Mathlib
with baked search path) and includes a negative control: a
deliberately false proof must be rejected on every run, so a broken
toolchain cannot report green. On machines without the toolchain the
check fails loudly rather than skipping silently; the sympy-only core
of this repository does not depend on it.

`LatticeExponents.lean` (2026-09-07) —

- `IntegrableBoundary.LatticeExponents.sol3_exponent_on_lattice`,
  `sol1_exponent_on_lattice`, `sol2_Y_exponent_on_lattice`,
  `sol2_X_exponent_on_lattice`, `sol3_coupling_at_lattice`,
  `sol1_lifting_spin_even`: **ansatz note sec. 38(w)** — at each sheet's
  double-lattice point t(s, j) the momentum-carrying exponents of the
  pillow-type leading symbol sit on the WKB Beta lattice (-2j/s and
  -2j'/s), the coupling at Solution 3's point is (s-4j)/(2(s-3j)), and
  sheet 1's lifting condition forces an even spin. Field identities over
  ℚ with explicit non-degeneracy hypotheses; the lattice and the
  exponent assignment are the inputs.

`KernelContiguity.lean` (2026-09-15, contributed by the cognitive agent) —

- `IntegrableBoundary.KernelContiguity.contiguity`: **ansatz note sec.
  38(az) item 2** — for every spin index k and every i <= k, the shifted
  top-layer kernels w(k, i; dX, dY) = C(k, i) (xY+dY)_i (xX+dX)_(k-i)
  r^(k-i) / (xY+dY)_k satisfy the three-term Gauss contiguity relation
  xX w(1,0) + (xY+k) w(0,1) = (xX+xY+k) w(0,0). Proved from the two
  shift identities x (x+1)_n = (x)_n (x+n) (`rf_shift`), with explicit
  non-degeneracy hypotheses xX, xY, (xY)_k, xY+k nonzero. Field
  identities over Q; the kernel and the shifts are the inputs (the lab
  had checked k = 2..9 symbolically). Rung: PROVED.
