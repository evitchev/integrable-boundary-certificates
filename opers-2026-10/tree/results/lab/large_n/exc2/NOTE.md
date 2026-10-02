# EXC2 (Fable seat, 2026-10-02) -- Solution 3: the level-1 excited-state oper (fourth-order three-term ODE)

Seals: `SEAL_EXC2.md` cfba205c... (03:32:50Z, one-point class); `SEAL_EXC2_addendum.md` 94b3c010... (03:50:09Z,
two-point class, tests U1/U2/U3); prediction `PREDICTION_EXC2.json` bc0aafcd... (04:48:57Z).  Each was sent to the lead
before the computation it governs.  The Sol-3 I_5 level-1 block was computed at this seat at 04:50:53Z, after the lead's
"registered" (the lead compared first).
STATUS: U1, U2, U3 PASS.  The one-point class of the original seal FAILS its T2 (it is half a level).

## Verdict

**1. Solution 3 has trivial-monodromy excited-state opers at level 1, and they are explicit.**  In the class-U normal
form (xi = x^c up to scale, th = xi d/dxi, b = p^2 = k/(k+1), xi-exponents a0 = p PX, a1 = p pi; A = a0^2, B = a1^2):

    L = (th^2 - A)(th^2 - B) + xi (th + 1/2) - E xi^b + xi sum_{+-} [ R_2 th^2 + R_1 th + R_0 ],

with TWO apparent singularities, at xi = +z and xi = -z, local exponents {-1, 1, 2, 4} at each, no logarithms identically
in E:

    at xi = +z:  R_2 = -b/(xi-z) - 4z/(xi-z)^2,   R_1 = r11/(xi-z) + r12/(xi-z)^2 + 8z^2/(xi-z)^3,
                 R_0 = r01/(xi-z) + r02/(xi-z)^2 + r03/(xi-z)^3 - 8z^3/(xi-z)^4;
    at xi = -z:  the same with z -> -z and (r11, r12, r01, r02, r03) -> (s11, s12, s01, s02, s03);

    (3b-4) u^4 - 8 (b-2)(2A + 2B + b^2 - 2) u^2 - 16 b ((b-2)^2 - 4A)((b-2)^2 - 4B) = 0,          [F]
    z = (u/16)(u^2 - 8(A+B) + 4b^2 - 8b + 8),
    r11 = -s11 = b u/2,   r12 = (u + b + 4) z,   s12 = (u - b - 4) z,   r03 = -z (r12 + 4z),   s03 = z (s12 - 4z),
    r01 = s01 = -(b/16)(3u^2 - 16(A+B) + 4(b-2)^2),
    r02 = -(z/4)(u^2 + (b+2)u - 8(A+B) + 2b^2 - 6b + 16),   s02 = +(z/4)(u^2 - (b+2)u - 8(A+B) + 2b^2 - 6b + 16).

(u, z) and (-u, -z) are the same operator; the TWO roots U = u^2 of [F] are the TWO level-1 singlet states.  No fitted
number.  Charges: the xi -> infinity effective quartic is th^4 - (A + B + 2b) th^2 + AB - (b/8)(3U - 16(A+B) + 4(b-2)^2)
(I_1 = Delta + 1; I_3 linear in U), and the first-order 1/xi term is (b+2) U (U - 8(A+B) + 4b^2 - 8b + 8)/16 * th
(I_5 quadratic in U).

**2. Evidence (all exact).**
- U2, target: trace and determinant of the level-1 I_3 block, as polynomials in (P^2, Q^2), are identical to the stored
  blocks at t = 2 and t = 9/4 (interpolated from 36 momentum points per fibre with held-back points, `u3_interp.py`; and
  again from the closed form with no interpolation, `s_closed2.py`).
- U3, BLIND hold-out: I_5 trace (10 coefficients) and determinant (28 coefficients) at each of the two fibres, plus the
  characteristic polynomial at 72 momentum points: 0 mismatches (`u3_compare.py real`, exit 0).  The lead compared first
  from its own computation of the block: exact at both fibres.  Controls fire: data entry + 1/1000 (1/10 trace, 10/28 det
  coefficients, 36/36 points per fibre); sign of the 1/xi term flipped on the oper side (6/10, 21/28); each fibre against
  the other fibre's prediction (8/10, 26/28).
- Sealed tamper (middle term xi(th + 1) instead of xi(th + 1/2), run through the two-point pipeline in `tamper_mid/`):
  the I_3 eliminant becomes irreducible of degree 14 and the data polynomial no longer divides it.

**3. The selection rule: self-adjointness.**  The sealed two-point system (20 unknowns, 6 fixed by the indicial
conditions, 14 no-logarithm equations for 14 unknowns) is zero-dimensional with 14 unordered solutions at generic
momenta; I_1 = Delta + 1 on all of them; the eliminant of I_3 is (quadratic) x (sextic) and the quadratic is the data
polynomial (U1, U2 as sealed; six generic points, two fibres).  Found after the seal and labelled as such: the adjoint
composed with xi -> -xi maps the vacuum operator to itself and the solution set to itself; the two level-1 states are
EXACTLY the self-adjoint solutions z2 = -z1.  On them r11 + s11 = 0 and z1 + z2 = 0, so the odd-order charges
J_3 = -A3 and J_5 = A2 A3/3 + 15 A3/2 - B2/7 (the t = 2 form; even spins 2 and 4, which Sol 3 does not have) vanish.  The other twelve
come in adjoint pairs (hence the sextic) and carry even-spin charges; they are not states of Sol 3.  I have not
identified them.

**4. One apparent singularity is half a level** (sealed T2 of the original seal FAILS).  Exponent set {-1,1,2,4}: four
solutions (so trivial monodromy is possible with one point), but r21 = -b gives I_1 = Delta + 1/2, the solutions are not
self-adjoint, and the I_3 values are not the data roots.  The other sealed exponent sets have NO solution at t = 2,
(a0, a1) = (1/2, 1/3): {-1,0,3,4} by an explicit obstruction (the coefficient of E at the resonance (-1,3) is -2592 z^6);
{-1,0,2,5} and {-2,1,3,4} by a unit Groebner basis (the last one with Singular; the sympy run had timed out).

**5. The exceptional locus is a = 1 - b/2, and it is a level-1 singular vector** (post-hoc; replaces my interim
statement "a0 = b or a1 = b", which was a coincidence of the fibre t = 2 where 1 - b/2 = b).  The constant term of [F]
vanishes iff a0 = +-(1 - b/2) or a1 = +-(1 - b/2).  There one root is U = 0: the two apparent singularities merge into
xi = 0, and the closed form degenerates to the VACUUM operator with that exponent replaced by 1 + b/2 (the effective
quartic factorises as (th^2 - (1+b/2)^2)(th^2 - B)).  Checked on the record's blocks: at a0 = 1 - b/2 (and likewise a1),
the vacuum eigenvalue at exponent 1 + b/2 is an eigenvalue of the level-1 block, for I_3 and I_5, both fibres, identically
in the other momentum (`x_exceptional.py`, exit 0; control exponent + 1/10 fails).  At t = 9/4 the two-point pipeline
confirms a0 = 8/13 exceptional and a0 = 10/13 = b not.  In momenta: PX = 1/p - p/2 -> 1/p + p/2.

**6. Engine for I_5** (`wkb3x.py`): the three-term WKB of VIR8 extended by a first-order term x^(-c)(B2 T^2 + B1 T + B0).
Validated without data: reproduces EXC1's Schrodinger-form result on the Sol-1 template (`x0_validate.py`); translation
covariance T -> T + d of all charges of orders 2..6, which exercises the T^3 coefficient and the j = 0, 1, 2 paths
(`v2_translation.py`, t = 2 and 9/4; its tamper fails); the 1/xi term enters I_5 only through q11 + q12, which is the
adjoint-invariant combination at both fibres.  At t = 2:
cB6 = -7 A2^3/216 - 15 A2^2/8 + A2 A4/6 - 207 A2/8 + A3^2/12 + 15 A4/2 - B1/16 + 9 B2/16 - 837/8.

## Parameter count against conditions
One point: 10 unknowns, 3 indicial + 7 no-log equations (6 pairs of exponents, split in E and the free Frobenius
constants; exactly one carries E and it is linear): zero-dimensional, 4 solutions.  Two points: 20 unknowns, 6 + 14
equations: zero-dimensional, 14 unordered solutions.  Self-adjoint two points: after the linear eliminations three
unknowns (r11, rho = r12/z, z) and five equations of which three are independent: r11 = b(rho-b-4)/2, z = cubic in rho, and
the quartic [F] (the two remaining equations both reduce to +-[F]).  Nothing is fitted to data at any step.

## Scope: what is and is not shown
- OPERATOR level: the no-logarithm conditions are derived symbolically in A, B and b (`s_sym.py`) and the closed form
  reproduces the numeric two-point computation at all 72 batch points (z-eliminant, `s_closed.py`).
- CHARGE level: I_1, I_3, I_5 on the level-1 singlet block, at t = 2 and t = 9/4 only (the WKB engine runs per fibre).
  Not checked: other fibres, I_7 and above (needs the second order in 1/xi), level 2, the transverse state (its
  eigenvalue is the level-shifted vacuum on the data side; no apparent-singularity oper is involved).
- Masoero-Raimondo: I make no claim against their statement as the lead summarised it (no non-trivial FFH-type
  trivial-monodromy solutions for twisted algebras).  The apparent singularities here sit in the scalar three-term
  (class-U) form with local exponents {-1,1,2,4} in the variable xi; I have not derived the dictionary to an FFH
  connection and name no algebra.  What is established is only: this fourth-order operator has trivial-monodromy
  deformations, and two of them reproduce the level-1 data including a blind hold-out.
- The twelve non-self-adjoint two-point solutions and the four one-point solutions are unidentified.

## Errors, deviations, disclosures
- Seal formula (i) dropped a factor p: the xi-exponents are p PX, p pi (corrected in the addendum before the two-point runs).
- Sealed predictions: T1 "some set has solutions" 35% -> true; T2 given T1 60% -> FALSE for the sealed one-point class;
  the level-1 oper needs two points (addendum U1 55%, U2 40% -> both true).  T3/U3 85% -> true.
- The addendum's remark "105/64 equals the vacuum I_3 at half-shifted exponents" is a coincidence of the momentum point;
  the control `v1_onepoint.py` built on it FAILS (exit 2) and was mis-posed, not an engine fault.
- Interim message to the lead said the exceptional locus is a0 = b or a1 = b; correct statement in item 5.
- `f_template.py` v0 printed FAIL through an overdetermined solve() (log kept); `u3_interp.py` v0 died on sympify('ff')
  without locals (logs kept); `s_sym.py` v0 had a wrong linear-elimination order (log kept).
- `u3_interp.py flip` exits 0: it only checks I_3 and polynomiality, and the flipped I_5 is also a polynomial.  So
  "I_5 interpolates to a polynomial" is NOT a discriminating control; the discriminating ones are in `u3_compare.py`.
- The selection of the level-1 solutions used for the prediction is self-adjointness, not the data I_3 factor named in
  the addendum; at the point (3/7, 2/9), t = 2, both selections were run and give the same I_5 polynomial.
- Two-point system v1 (degree 35) and v2 (14 unknowns) were too slow in Singular/msolve and were killed by saved PID
  (`pids.txt`); the reduced system (4 unknowns) solves in seconds.  `num_two.py` (mpmath Newton) is unused.
- `h_check.py` prints "I_1 = Delta + 1 on all" as fixed text (it follows from r21 = s21 = -b, asserted in `u3_point.py`).
- Exit codes: f_template 0; f_sol3 S_a 0, S_b 3, S_d 3, S_c 124 (timeout) then f_sing 3; g_compare 2 (one-point T2 fails);
  h_check 0 at six generic points, 2 at the exceptional ones; x0_validate 0; v1_onepoint 2 (mis-posed); v2_translation 0
  (tamper 0 = fires); w_selfadj 0 at seven generic points, 2 at five exceptional ones; u3_interp 0 (both fibres);
  u3_compare real 0, tamper 0, flip 0, cross 0 (all three fire); s_closed 0; s_closed2 real 0, tamper 0 (fires);
  x_exceptional 0; tamper_mid h_check 2 (fires).

## Files
`SEAL_EXC2.md`, `SEAL_EXC2_addendum.md`, `PREDICTION_EXC2.json` (+ .sha256, .time), `HOLDOUT_OPENED.time`, `frob.py`,
`f_template.py`, `f_sol3.py`, `f_diag.py`, `f_sym.py`, `f_sing.py`, `g_charges.py`, `wkb3g.py`, `g_compare.py`,
`h_two.py`, `h_two_v2.py`, `h_reduce.py`, `h_reduce2.py`, `h_run.sh`, `h_check.py`, `w_selfadj.py`, `wkb3x.py`,
`x0_validate.py`, `v1_onepoint.py`, `v2_translation.py`, `u3_engine.py`, `u3_predict.py` (first single-point version),
`u3_point.py`, `u3_one.sh`, `u3_points_t*.txt`, `u3_interp.py`, `u3_compare.py`, `d1_data_exc2.py`,
`d1_sol3_t*_with5.json`, `s_sym.py`, `s_closed.py`, `s_closed2.py`, `x_exceptional.py`, `tamper_mid/`, logs, `SHA256SUMS`.
