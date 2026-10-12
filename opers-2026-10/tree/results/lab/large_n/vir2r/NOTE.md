# VIR2-R -- independent replication of Fable's VIR2: a two-boson screening system for the cylindrical charges

cc (Opus 5.5 seat), 2026-10-01. Lead's commission.
- Seal: `SEAL_VIR2R.md`, sha256 efa0c7f58f3f0f76821db65ca9ba075a04fc0371b416e557afd318d2ee21f512, sent before computing.
- Fable's stage was NOT opened.
- Exit codes: `ctrl_e.py` 0; `vir2r.py` 0; `posthoc_viralone.py` 0.

## Result: VIR2 REPLICATED IN FULL, as sealed
- **(P1) Conservation.** Every listed screening pair V_+- = exp(+-a X + beta phi) is conserved, exactly, against the record's I_3, I_5 and I_7, at all 3 fibres of each solution.
  - Sol 1 (beta = -p/2, a^2 = 1 - p^2/4): 9 of 9.
  - Sol 2 (beta = -p, a^2 = 1 - p^2): 9 of 9.
  - Sol 3 pair A (beta = 1/p, a^2 = 1/p^2): 9 of 9.
  - Sol 3 pair B (beta = -p, a^2 = p^2): 9 of 9.
  - Both Virasoro screenings exp(p phi) and exp(-(2/p) phi) are conserved for every translated charge (control c).
- **(P2) Determination.** The commutant of {exp(p phi), V_+, V_-} among ALL two-boson densities, modulo d, is **exactly 1-dimensional at weights 4, 6, 8** for each solution (first fibre: Sol 1 p = 8/5, Sol 2 p = 4/5, Sol 3 p = 1/2 with pair A). It is spanned by the record's I_3, I_5, I_7: each is in the commutant and is not a total derivative.
  - The commutant projection ranks are 11, 37, 111, against total-derivative images of 10, 36, 110.
- **VIR1, as a by-product:** at every one of the 9 fibres (non-integer N, e.g. N = 253/400), I_3, I_5, I_7 decompose EXACTLY into Heis_X ⊗ Vir_N modulo d (residual 0). This extends VIR1-R from integer N to generic fibres.
- **Post-hoc (labelled; `posthoc_viralone.py`), Sol 1 p = 8/5:** the commutant mod d at weights 4/6/8 is
  - 5 / 13 / 30 with exp(p phi) alone;
  - 2 / 5 / 11 with exp(p phi) + V_+;
  - 1 / 1 / 1 with the full pair.
  So the X-screening pair is what pins the charges down.

## Method (own code; no code shared with Fable or mixed_engine)
- **Inputs:** the record's I_3 = cyl_P4 and the genuine kernels I_5, I_7 at curve_point(t), from mixed_engine (inputs only). Fibres have rational p, t = (2 + p^2)/(2 - p^2), and N = 1 - 12 rho^2 asserted at each fibre.
- **`onwick.py`:** my own O(N) pair-language Wick engine with formal (rational) N. Dummy labels per pair; Y-Y contractions unify labels; a closed cycle gives N.
  - **Control (e)** (`ctrl_e.py`): EXACTLY equal to my explicit-component engine at N = 3, 4 on :TT:, :T:TT::, :dT T:, the T T 4-pole, and a mixed X-Y product.
- **Decomposition:** I_k = sum x_j (X monomial × Vir word in T_Y) + d(V), solved exactly in pair coordinates.
  - The nullspace (5 / 16 / 45 at W = 4 / 6 / 8) consists of basis combinations that are total derivatives. **EVERY ambiguous direction was checked to map to a total derivative in (X, phi).**
- **Translation:** each Vir word is rebuilt from T_phi = 1/2 (dphi)^2 + rho d^2phi (rho = p/2 - 1/p) with my explicit Wick code.
- **Conservation:** the residue of G(z) exp(alpha.phi)(w), and an exact solve of R = dQ + (alpha.a) Q, in Q(sqrt d) (own field class).
- **Commutant:** an exact Fraction nullspace over all two-boson monomials, with the joint conditions of the three exponentials; dim mod d = projection rank - #monomials(W-1).

## Controls (all pass)
- (a) a^2 + 1/10: FAILS on I_5, all pairs, first fibre of each solution.
- (b) rho -> -rho in the translation: FAILS.
- (c) Virasoro screenings: PASS everywhere.
- (d) generic (a^2, beta) = (2/7, 3/11): FAILS.
- (e) formal-N vs explicit engine: exact.
- (f) a total-derivative sanity element passes.

## Disclosures
- The commutant was run at the first fibre of each solution only, as sealed. I_7 conservation was run at all fibres.
- `posthoc_viralone.py` is post-hoc and informational.
- No deviation from the seal otherwise.

## Files
- Scripts: `wick.py` (= WINF1's Wick code), `onwick.py`, `ctrl_e.py`, `vir2r.py`, `posthoc_viralone.py`.
- Logs: `run_*.log` and `run_*.stdout`.
- JSON: `vir2r_conservation.json`, and `vir2r_translated.json` (the full (X, phi) densities of I_3, I_5, I_7 at all fibres).
