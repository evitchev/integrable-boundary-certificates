# EXC1 (Fable seat, 2026-10-02) -- excited states: the level-1 and level-2 opers of Solution 1

Seal `SEAL_EXC1.md` 7e3297d0... (03:17:57Z); hold-out prediction `PREDICTION_EXC1.json` de8366b1... (03:23:58Z),
registered by the lead before the I_5 level-1 data were computed.

## Verdict

**Solution 1 has BLZ-type excited-state opers, with apparent singularities in the xi = x^(M+1) plane, one per level.**

    -psi'' + [ x^(2M) + alpha x^(M-1) + (lam^2 - 1/4)/x^2 - 2 d_x^2 sum_{k=1}^{L} log(x^(M+1) - z_k) ] psi = E psi,
    M = (t+3)/(t-1),   lam^2 = (M+1) pi^2,   alpha^2 = 4M(M+1) PX^2,

with trivial monodromy at every root of x^(M+1) = z_k.  No fitted number anywhere.

- LEVEL 1 (one point).  Trivial monodromy <=> the quadratic
        2M z^2 + alpha (M-1) z - 2 lam^2 + (M-1)^2/2 = 0
  (derived by hand; re-derived by series at M = 5).  I_1 = Delta + 1 exactly.  Its two roots are the two states of the
  level-1 block:  Res_z(quadratic, e - Lambda_3(z)) == the characteristic polynomial of the I_3 data block, as
  polynomials in e with coefficients in (PX^2, pi^2), at t = 2 AND t = 9/4 (`o2_compare.py` exit 0; tamper fires).
  HOLD-OUT (O3, sealed 80% conditional): the I_5 characteristic polynomial predicted from the oper and hash-registered
  before the I_5 level-1 matrix existed: EQUAL at both fibres (`o3_compare.py` exit 0; tamper fires).
  Branch locus = discriminant of the quadratic: (M-1)^2 alpha^2 + 16 M lam^2 - 4M(M-1)^2 (t = 2: 160(12 PX^2 + 3 pi^2 - 2);
  the data's discriminant is 768 PX^2 (3 pi^2 + 12 PX^2 - 2)/49).  This is the one-root Bethe equation that 38(l) asked for.
- LEVEL 2 (two points; extension beyond the seal's tests, labelled).  The coupled system sealed as (v) (re-derived by
  series at M = 5) has I_1 = Delta + 2 and, at t = 2, bare (P, Q) = (1/2, 1), exact Groebner elimination gives ONE
  irreducible polynomial of degree 5 for e = Lambda_3 -- and it EQUALS the quintic factor of the 6x6 level-2 data block
  (`o5_level2_compare.py` exit 0).  5 = the number of level-2 states of two bosons; the sixth state of the N-component
  block is the level-shift state.  One momentum point, one fibre: a check, not a survey.

## What was wrong in the seal (reported as such)
- O1 ("one point is level 1/2, a vacuum at shifted momenta", 60%): REFUTED.  One point has I_1 = Delta + 1.
- O2 ("two points are level 1", 45%): REFUTED as stated.  Two points have I_1 = Delta + 2; the sealed two-point system
  is correct but it is the LEVEL-2 system.  The reason for my mistake: in Suzuki's normalisation lam^2 = (M+1) pi^2 while
  in BLZ's (l + 1/2)^2 = 4(M+1) Delta, so one xi-point carries the level that 2M+2 x-points carry in BLZ.
- Consequently the level-1 comparison for I_3 is post-hoc with respect to the seal's point count (zero parameters,
  symbolic identity at two fibres); the I_5 comparison is a registered hold-out; level 2 is an unsealed extension.
- The sealed BLZ control (alpha = 0, z2 = -z1) was mis-posed for the same reason and was NOT run.  At alpha = 0 the two
  level-1 roots are +-z and the two eigenvalues coincide, in agreement with the factor PX^2 of the data's discriminant.

## Data side (D1; lab/excited_states.py imported read-only)
Level-1 O(N-1)-singlet block (two states: X oscillator, Y oscillator along the momentum; in two-boson language the
level-1 subspace of the Fock module, i.e. descendants of the primary) of I_3 for Sol 1 and Sol 3 at t = 2, 9/4: irreducible
quadratic in all four cells; the transverse state (a new Vir_N primary of weight h_Y + 1) has the vacuum eigenvalue at
pi^2 -> pi^2 + 2 (level-shift law) in all four.  I_5 (from my VIR1 kernel) for Sol 1 at both fibres: [I_3, I_5] = 0 on the
block.  Level 2, Sol 1, t = 2, (P, Q) = (1/2, 1): characteristic polynomial = (linear, the level-shift value) x (quintic).
Sol 3: the I_3 level-1 blocks are stored (`d1_sol3_t*.json`); its I_5 level-1 block has NOT been computed (kept as a
hold-out for EXC2).

## Engine
`mellin_pot.py` = the Mellin-WKB engine with (a) potential terms Lambda = y^n p (1 + sum_w V_w), (b) the charge
extraction generalised to y-powers a = -i + sigma m.  Vacuum control (`o0_vacuum.py`, exit 0): the Schrodinger equation
with the alpha term reproduces VIR7's Gamma-form S1 charges at spins 1..7 (t = 9/4, 7/3) -- a direct WKB check of the
Suzuki <-> Gamma-form equivalence at all layers.  The apparent singularities enter through the large-xi series
x^2 dV = 2(M+1) sum_k [1 + sum_j (j(M+1)+1)(z_k/xi)^j], order by order (z has weight 1).

## Errors and disclosures
- Two bugs in the first version of the charge extraction (alpha-odd terms have half-integer Beta shifts; they are a total
  derivative at even orders, now checked to integrate to zero); found by the vacuum control, before any excited-state run.
- `sympify('Q')` returns sympy's assumptions object; all parsing uses explicit locals.
- Solution 3: NO oper construction attempted (sealed).  QLIT6 / Masoero-Raimondo 2312.01955 (lead's note) says FFH-type
  apparent singularities have no non-trivial solutions for twisted algebras; the three-term ODE is not an FFH connection,
  so the question is open.  Proposal for EXC2: operators sum_j xi^j A_j(T) + E xi^b (xi - z)^m with deg A_j <= 4, regular
  singular at xi = z with integer exponents and no logarithms, exponents +-l0, +-l1 at 0; target = the stored Sol-3 I_3
  level-1 blocks, hold-out = I_5.

## Exit codes
d1_data (six runs) 0; o0_vacuum 0; o2_compare real 0, tamper 0 (fires); o3_compare real 0, tamper 0 (fires);
o4_level2 0; o5_level2_compare 0.

## Files
`SEAL_EXC1.md`, `PREDICTION_EXC1.json` (+ .sha256, .time), `mellin_pot.py`, `d1_data.py`, `d2_level2.py`, `o0_vacuum.py`,
`o1_engine.py`, `o2_compare.py`, `o3_predict.py`, `o3_compare.py`, `o4_level2.py`, `o5_level2_compare.py`, JSON outputs,
logs, `vir_lib.py`, `mellin_gen.py`, `mellin_wkb.py` (copies), `SHA256SUMS`.
