# JOINT1 (Fable secondary seat, 2026-10-10) -- state-by-state check of the excited-state operators: the (I_3, I_5) values pair up on the same states

Seal: `SEAL_JOINT1.md` 9c58c138... (2026-10-10T19:00:53Z, sent to the lead before any computation).  Prompted by ChatGPT's review (the paper
compares the characteristic polynomials of I_3 and I_5 separately).  Inputs: the archived data blocks and closed forms of EXC1/EXC2/EXC3,
copied unchanged into `inputs/` (`INPUTS.SHA256`; every hash equals the stage manifest's).  Nothing archived was recomputed.
RESULT: all six level-one cells pair correctly, exactly; the crossed pairing and the +1/10 momentum shift fail on every cell.

## Table
| cell | stage, solution, t | block | [M_3, M_5] = 0 (symbolic) | degenerate at the points? | joint check (symbolic / points) | pairs equal (1/2,1/3), (2/5,1/4) | crossed fires | shift fires |
|---|---|---|---|---|---|---|---|---|
| C1 | EXC1, Sol 1, 2   | 2 x 2 | yes | no (disc = 2304 P^2 (4P^2 + Q^2)/49) | SYMBOLIC in (P, Q): tr(M_3 M_5), tr(M_3^2 M_5), tr(M_3 M_5^2) identical | yes, yes | yes | yes |
| C2 | EXC1, Sol 1, 9/4 | 2 x 2 | yes | no (disc = 256 P^2 (64P^2 + 25Q^2)/81) | SYMBOLIC, identical | yes, yes | yes | yes |
| C3 | EXC2, Sol 3, 2   | 2 x 2 | yes | no | SYMBOLIC, identical | yes, yes | yes | yes |
| C4 | EXC2, Sol 3, 9/4 | 2 x 2 | yes | no | SYMBOLIC, identical | yes, yes | yes | yes |
| C5 | EXC3, Sol 2, 2   | 2 x 2 | yes | no | at the two points (exact rationals), identical | yes, yes | yes | yes |
| C6 | EXC3, Sol 2, 9/4 | 2 x 2 | yes | no | at the two points, identical | yes, yes | yes | yes |
Level two: no archived I_5 block exists (EXC1's and EXC4's level-two comparisons are I_3 only), so there is no level-two cell; stated in the seal.

## Method (as sealed)
For a 2 x 2 block with commuting M_3, M_5 the multiset of joint pairs {(a1, b1), (a2, b2)} is fixed by the marginals (tr, det of each, already
compared in the archives) and the MIXED TRACE tr(M_3 M_5) = a1 b1 + a2 b2; the crossed pairing has a1 b2 + a2 b1 = tr M_3 tr M_5 - tr(M_3 M_5),
differing by (a1 - a2)(b1 - b2), which is non-zero off the degeneracy locus.  Operator side: the sum over the operator's two selected solutions of
I_3 I_5 (and of I_3^2 I_5, I_3 I_5^2, I_3^2), computed exactly as traces of multiplication operators in the quotient ring of the solution
polynomial -- EXC1: Q(alpha, lam)[z]/(2M z^2 + alpha(M-1) z - 2 lam^2 + (M-1)^2/2) with R_4, R_6 from `o1_t<t>_n1.json` and EXC1's dictionary;
EXC2: Q(A, B)[U]/(U^2 + c1 U + c0) with I_3, I_5 polynomial in U from `s_closed.json`/`s_sym.json`/`u3_engine_t<t>.json`; EXC3: Q[z]/(quartic),
the four roots counting each unordered solution twice, with the rational-in-z charges from `s_closed_sol2.json`/`g_engine_t<t>.json` (the
inverse of the denominator taken in the quotient ring).  C1-C4 symbolic in bare (P, Q) (polynomial identities); C5-C6 at the two rational points.
Explicit pairs: data -- the I_3 eigenvectors over Q(sqrt(disc)) and the I_5 eigenvalue on each; operator -- the roots of e^2 - s3 e + p3 (its own
s3 and p3 = (s3^2 - tr I_3^2)/2) with b_i = (T11 - a_j s5)/(a_i - a_j); equality as sets of pairs in the quadratic field (sympy radsimp/simplify).
Tampers after the hash guard: 'crossed' (the operator's pairs crossed: would need a1 b2 + a2 b1 = tr(M_3 M_5)) and 'shift' (P -> P + 1/10 on the
operator side).  Every cell prints (a1 - a2)(b1 - b2) so that the crossed tamper is non-vacuous (e.g. C5 at (1/2,1/3): -811400/243).

## Numbers (all exact)
- C1 (1/2,1/3): pairs (11803/45360 -+ 4 sqrt10/7, (17812609 -+ 1632960 sqrt10)/8164800) on both sides, the sign of sqrt10 correlated.
- C2 (1/2,1/3): the discriminant is a rational square there (pairs rational): (-54234361/29568240, 2991211705361/1199458653120),
  (59658119/29568240, 3377876674961/1199458653120) on both sides.
- C3 (1/2,1/3): (4417/1296 -+ 2 sqrt205/9, I_5 correlated); C4, C5, C6: see the logs (`C*_real.log`), quadratic irrationalities in every case.
- Degeneracy loci: C1/C2 P^2 (4P^2 + Q^2) resp. P^2 (64P^2 + 25Q^2) -- the two singlet states are degenerate in I_3 only at P = 0 (where they
  cannot mix; the sealed points avoid it); C3/C4 quartic curves; C5/C6 quartic curves (printed in the logs).

## Priors vs outcome
All six cells pair correctly: 85% -> true.  Both tampers fire on every cell: 95% -> true (12 of 12).

## Scope and disclosures
- This is a check of the archived level-one comparisons at the fibres and momentum points of the record; it adds the joint invariant and the
  explicit pairing that the archives lacked.  It does not extend the operators to new states.
- One script, `joint1.py` (cell, mode); its first draft had two dead lines (an aborted operator-pairs construction) removed before the first run
  and the operator's own determinant p3 = (s3^2 - tr I_3^2)/2 added in place of the data's p3; the smoke test (C1) and the full run used the final
  version.  Logs: `C<n>_<mode>.log`, 18 runs, all exit 0 (real = PASS; crossed/shift = the tamper fires).
- Runs are seconds each; capped (ulimit -v 16 GB) all the same.

## Files
`SEAL_JOINT1.md` (+ .sha256, .time), `joint1.py`, `inputs/` (15 archived files) + `INPUTS.SHA256`, `C1..C6_{real,crossed,shift}.log`, `NOTE.md`, `SHA256SUMS`.
