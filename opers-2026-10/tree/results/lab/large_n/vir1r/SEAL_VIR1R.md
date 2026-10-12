# SEAL VIR1-R -- independent explicit-component replication of Fable's VIR1 (lead's commission, 2026-10-01)

cc (Opus 5.5 seat), sealed BEFORE any computation; the hash goes to the lead before running.
- Fable's VIR1 (~/fable-work/vir1, NOTE 49cc436d, not opened by me): I_3..I_11 of Sols 1-3 lie in Heis_X ⊗ Vir_N modulo total derivatives, identically in t.
- No code is shared with Fable or with mixed_engine for the TEST.
- mixed_engine (w8 clone, read-only) is used ONLY to produce the record's charge representatives (inputs), exactly as the brief asks.

## Object
- **Inputs:** the record's charge densities at a curve point (N, s), N = (t^2-25)/(t^2-1), s = -24t/(t^2-1):
  - I_3 = cyl_P4(sol, N, s);
  - I_5 and I_7 = the 1-dimensional genuine kernels at weights 6, 8 (mixed_engine.genuine_kernel_sparse, as in cyl_membership).
  - Format: (xs, yp), with X derivative orders xs and Y pairs (a, b) = d^aY . d^bY.
- **Translation to EXPLICIT components:** letter (x, p) = d^(p+1) X, and letter (i, p) = d^(p+1) Y_i, i = 1..Ncomp. A pair (a, b) becomes sum_i (i, a-1)(i, b-1). Free-field normal ordering is commutative.
- **Convention check (read, not computed):** mixed_engine's contraction (-1)^(m-1)(m+n-1)! (z-w)^-(m+n) for d^m phi d^n phi equals my winf1 Wick code's (-1)^p (p+q+1)! with m = p + 1. The same OPE, so T_Y = 1/2 sum_i (i,0)(i,0) has c = N.
- **Test space at weight W:** span{ (X monomial of weight W - w2) × (Vir_N field of weight w2) } + d( V_(W-1) ), where:
  - the Vir_N fields are built by MY normal-ordered products (the r = 0 term of my Wick OPE, winf1's `ope`):
    - w2 = 0: 1; 2: T; 3: dT; 4: d^2T, :TT:; 5: d^3T, :(dT)T:;
    - w2 = 6: d^4T, :(d^2T)T:, :(dT)(dT):, :T:TT::;
    - w2 = 7 and 8 analogously (all L_-n words, i.e. normally ordered products of derivatives of T, built recursively);
  - V_(W-1) = ALL weight-(W-1) O(N)-invariant mixed fields (X monomials × products of Y pairs), in explicit components;
  - d is the total derivative.
- **Membership test:** target ∈ test space?
  - Method: evaluation at M random points (letters -> random residues) MODULO two independent primes p1, p2 ~ 2^61, with M = 3 × (number of spanning fields).
  - The column rank of [spanning fields] is compared with that of [spanning fields | target].
- **Certification (disclosed):**
  - OUT (rank increases mod p) is certified unless p divides a denominator of the true coefficients (it must hold for BOTH primes).
  - IN is probabilistic (Schwartz-Zippel; error <= deg/p per point, negligible), and must hold for BOTH primes.

## Points (rational s)

| point | t | s | role |
|---|---|---|---|
| N = 25 | t = 0 | s = 0 | |
| N = 26 | t = 1/5 | s = 5 | |
| N = 26 | t = -1/5 | s = -5 | |
| N = 28 | t = -1/3 | s = -9 | positive control (exact point) |

- Sols 1, 2, 3 at every point.
- I_3 and I_5 at all points. I_7 at N = 25, t = 0, if feasible (stated as feasible / not in the NOTE).
- If N = 25/26 at weight 6 is not feasible: say which, and fall back to the feasible one.

## SEALED PREDICTION (the lead's): IN
- I_3 and I_5 (and I_7 where computed) are IN Heis_X ⊗ Vir_N mod total derivatives, for ALL three solutions at every point.
- Any OUT (for both primes) is a disagreement with VIR1, reported with its point and solution.

## Controls
- **(a) PLANT:** I_5 + (1/1000)(33)_Y, i.e. + (1/1000) sum_i (i,2)(i,2), must be flagged OUT at every point and solution.
- **(b) I_3** must pass (IN) everywhere. This is also the lead's spin-3 confirmation.
- **(c) TAMPER** (after the guards): translate the record's density with Ncomp = N + 1 components and build Vir/V with N + 1 components, labelled N. The result must BREAK (I_5 OUT) for at least the generic points.
- **(d) Sanity:** a known element, d(:TT:) × X-monomial plus an X-only term, must be IN. The Vir field count per weight must equal the generic Virasoro vacuum character (1, 0, 1, 1, 2, 2, 4, 4, 7 for weights 0..8).

## Workspace
~/opus5-work/tq/stageVIR1R. Repo read-only (mixed_engine via the w8 clone). Exact residues mod p. Kill by saved PID only.
