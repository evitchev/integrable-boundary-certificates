# EXC3 (Fable secondary seat, 2026-10-05) -- Solution 2: the level-one excited-state oper (third-order three-term ODE)

Seal: `SEAL_EXC3.md` 58146a18... (2026-10-05T03:24:08Z, snapshotted by the lead at 03:25:38Z, before any computation);
prediction `PREDICTION_EXC3.json` 7991f7af... (03:49:57Z), registered by the lead at 03:51:27Z; the Solution-2 I_5
level-one blocks were computed at this seat at 03:56:03Z, after the lead's "registered" and after the lead's own comparison.
STATUS: T0, T1, T2, T3 PASS.  The sealed one-point class is half a level (as for Sol 3); the two-point sigma-invariant
class is the level-one oper, with a blind I_5 hold-out passing at both fibres.

## Verdict

**1. Solution 2 has trivial-monodromy level-one opers, and they are explicit.**  Normal form (seal sec. 1, machine-checked
at T0): xi = x^c/c, vartheta = xi d/dxi, b = p^2 = k/(k+1), A = a0^2 = (1-b) P_X^2, B = a1^2 = b pi^2,

    L_0 = vartheta (vartheta^2 - B) + xi ((vartheta + 1/2)^2 - A) - E xi^b .

The level-one oper adds TWO apparent singularities at xi = +z and xi = -z, local exponents {-1, 1, 3} at each, no
logarithms identically in E:

    L = L_0 + xi sum_{+-} [ R_1 vartheta + R_0 ],
    at +z:  R_1 = -b/(xi-z) - 3z/(xi-z)^2,   R_0 = r01/(xi-z) + r02/(xi-z)^2 + 3z^2/(xi-z)^3;
    at -z:  R_1 = -b/(xi+z) + 3z/(xi+z)^2,   R_0 = s01/(xi+z) + s02/(xi+z)^2 + 3z^2/(xi+z)^3,

where U = z^2 is a root of the QUADRATIC (`s_closed_sol2.py`, symbolic in A, B, b; verified against the numeric
z-eliminants at 70 momentum points on two fibres)

    c7 U^2 + c5 U + c3 = 0,
    c7 = -16 (b-3)(b-1) ((1+b)^2 - 4A),
    c3 = -b^2 ((2-b)^2 - 4B) ((2b-1)^2 - 4B)^2,
    c5 = -432 A^2 (b-1)^2 + 288 A B (b-1)^2 + 576 A b^3 - 1224 A b^2 + 720 A b - 72 A + 16 B^2 (b^2 - 2b - 3)
         - 96 B b^4 + 96 B b^3 + 88 B b^2 - 80 B b + 24 B + 32 b^6 - 128 b^5 + 88 b^4 + 104 b^3 - 119 b^2 + 22 b - 3,

and r02, s02 are explicit rational functions of (z, A, B, b) (in `s_closed_sol2.json`; r01 = 2b r02/(3z) - (3-b)z/3 - b(b+3)/3,
s01 likewise with s02 and z -> -z).  Relabelling z -> -z exchanges the two points (r <-> s), so the two roots U are the two
level-one singlet states.  No fitted number anywhere.

**2. Evidence (all exact).**
- T2 (target): trace and determinant of the level-one I_3 block, as polynomials in the bare (P^2, Q^2), are IDENTICAL to the
  stored block (d1_sol2_t*.json, computed with lab/excited_states.py) at t = 2 and t = 9/4: interpolated from 34 and 36
  momentum points (trace fitted on 12, checked on 22/24; det fitted on 24, checked on 10/12), `u_interp.py` exit 0 both.
- T3 (BLIND hold-out): I_5 trace (10 coefficients) and determinant (28) at each fibre, plus the characteristic polynomial
  at 34 + 36 points: 0 mismatches (`u_compare.py real`, exit 0).  The lead computed the blocks and compared FIRST: exact
  at both fibres.  Controls fire: data entry + 1/1000 (1/10, 10/28, all points); each fibre against the other fibre's
  prediction (9/10, 27/28).
- Sealed tamper (middle term xi((vartheta+1)^2 - A) instead of xi((vartheta+1/2)^2 - A)): the symmetry sigma is destroyed
  (T0b), the one-point I_3 changes, and the two-point class has NO solutions (vdim 0).

**3. The symmetry carries over -- the seal's prediction S (90%) holds.**  For the pairing int f g dxi/xi, the formal adjoint
sends the odd P(vartheta) = vartheta(vartheta^2 - B) to -P and the complete-square middle term xi((vartheta+1/2)^2 - A) to
itself; composing with xi -> -xi gives sigma: L -> -(L^+)(-xi) with sigma(L_0) = L_0 (E -> -(-1)^b E), machine-checked
symbolically (`t0b_sigma.py`).  The brief expected the EXC2 rule not to carry over because P is odd; the odd P only
produces the overall sign.  Consequences, all confirmed: sigma reflects exponents about 1 (so {-1,1,3} is invariant and the
one-point solutions come in pairs +-z, eliminant 105 z^2 - 4 at (1/2,1/3), t = 2); the level-one states are the
sigma-INVARIANT two-point solutions (kk = s11 - r11 = 0 on all of them), and on them the even-spin charge J_3 vanishes.

**4. One apparent singularity is HALF a level** (sealed H1: 45% for 1/2).  The one-point solutions of type {-1,1,3} have
I_1 = Delta + 1/2 (the E-coefficient condition forces r11 = -b, exactly as r21 = -b did for Sol 3), I_3 = 2 at (1/2,1/3)
(the same on both members of the +-z pair, so the pair cannot give the data's two distinct roots), and non-zero even-spin
charges.  The other sealed exponent sets at t = 2, (1/2,1/3): {-1,0,4} and {-2,2,3} have no solution (unit Groebner basis);
{-2,1,4}: see sec. 8.

**5. Exceptional locus = three lines, each a level-one singular vector** (post-hoc, `x_exceptional_sol2.py`, exit 0).  The
quadratic degenerates where c7 = 0 or c3 = 0:
    a0 = (1+b)/2  [c7 = 0: one root U -> infinity];   a1 = 1 - b/2  or  a1 = b - 1/2  [c3 = 0: one root U = 0].
These are the batch's two one-solution points at t = 2 (b = 2/3): (5/6, 8/9) has a0 = 5/6 = (1+b)/2 and (5/8, 1/6) has
a1 = 1/6 = b - 1/2.  On each line the lost state is the VACUUM oper with the exponent shifted by a level-one amount
(A -> A + 2(1-b) or B -> B + 2b, i.e. a0 -> (3-b)/2, a1 -> 1 + b/2, a1 -> b + 1/2, each raising the monic I_1 by exactly 2):
checked on the record's blocks for I_3 and I_5, both fibres, identically in the other momentum (12/12; the +1/10 controls
fail).  The EXC2 line a = 1 - b/2 reappears for a1; the other two are new.

**6. The weight count of the seal was wrong, and that mattered.**  Seal sec. 4 said I_5 needs the 1/xi expansion to second
order.  For Sol 2 xi has weight one, so x^(-mc) T^j enters the order-i charge at i = 3 - j + m (j <= 1): I_3 needs
m <= 2 (the second-order marker C1 appears in cB4), I_5 needs m <= 4.  With the engine truncated at second order the I_5
"invariants" over the solutions were NOT polynomial in the momenta (the interpolation failed, as it should); after extending
`wkb3x2.py` to four orders in U = x^(-c) they are.  Validation BEFORE the prediction (`t0e_validate.py`, exit 0): EXC1's
complete Schrodinger-form R_2, R_4, R_6 of the level-one Sol 1 oper at t = 9/4 (M_S = 21/5), including the z^3 and z^4
terms of R_6, reproduced exactly; translation covariance T -> T + d of all charges i = 2..6 with all four orders present;
regression to the two-order engine at order 4.  (At t = 2, M_S = 5, the three-term form of the Sol 1 template has its
nu = 2 gap at spin 5 and returns 0 at order 6, so that fibre cannot test order 6; v1 log kept.)

**7. Dictionary** (seal sec. 1, checked at T0a against the certified vacuum I_1, I_3, I_5 at both fibres, exit 0; tamper
a0^2 x 11/10 fails):  P_xi = vartheta^3 + q2 vartheta^2 + q1 vartheta + q0 -> P_T = T^3 - c q2 T^2 + c^2 q1 T - c^3 q0;
(1/xi^m)(q_m1 vartheta + q_m0) -> x^(-mc) (c^(2+m) q_m1 T - c^(3+m) q_m0);  l0 = c a0 in L(T - c/2) = (T - c/2)^2 - l0^2;
q_mj = r_j1 z^m + m r_j2 z^(m-1) + m(m-1)/2 r_j3 z^(m-2) summed over the two points.  The vacuum's odd-order charges J_3, J_5
vanish (sealed H3, 70%: true).

## Parameter count against conditions
One point: 6 unknowns (r11, r12, r01, r02, r03, z); 2 indicial conditions (exponent sum 3 forces the indicial e^2
coefficient; r12 = -3z, r03 = 3z^2); 4 no-log equations in the 4 remaining unknowns (the E-linear one gives r11 = -b):
zero-dimensional, 2 solutions (+-z).  Two points at +-z: 7 unknowns after the indicial conditions, 8 equations:
zero-dimensional, 4 ordered = 2 unordered solutions.  Closed form: the quartic in z reduces to the quadratic in U = z^2
above; r02 and s02 are then determined linearly (Groebner basis over Q(A, B, b), `s_closed_sol2.out`).  Nothing fitted.

## Scope
- Operator level: the no-log conditions and the closed form are symbolic in A, B, b.  Charge level: I_1, I_3, I_5 at
  t = 2 and 9/4 only.  Not done: other fibres, I_7 (needs the sixth order in 1/xi), level 2, two points at general
  (z1, z2) beyond the count of sec. 8.
- Masoero-Raimondo: as in EXC2, no claim either way; this is the scalar three-term form in xi.

## 8. The remaining counts (completed)
- Exponent set {-2,1,4}, one point, t = 2, (1/2,1/3): NO solution (unit Groebner basis, Singular, `f_sing_ff_sol2.py`, exit 3).
  The sympy Groebner run and the first Singular route both spent their time in the Frobenius recursion (sympy `cancel` on
  nested rational functions, > 50 min); a fraction-free rewrite of the same recursion (`frob_ff.py`, common denominator
  = product of the indicial values, a polynomial in z) generates the conditions in a second and reproduces the {-1,1,3}
  result (2 solutions) as its control.  So of the four sealed sets only {-1,1,3} admits apparent singularities.
- Two points at free positions (z1, z2), both of type {-1,1,3}, t = 2, (1/2,1/3) (`h_two_sol2.py free2`; r11 = s11 = -b
  substituted, which are the local E-linear conditions): zero-dimensional, 12 ordered = 6 unordered solutions.  The I_3
  eliminant is (9e^2 - 120e + 2092)(3e - 38)(3e + 58): the data quadratic times two further values, each shared by an
  adjoint pair; z1 + z2 has the eliminant kk (kk - 8)(kk + 8)(175 kk^2 + 64), so exactly the two sigma-invariant solutions
  (z1 + z2 = 0) are the level-one states and the other four are two sigma-pairs (EXC2: 2 of 14).  The first, unreduced free
  run (9 unknowns) did not finish in 50 min (log kept).
- Closed form at charge level (`s_closed_pts.py`, exit 0 both fibres): at 5 rational momentum points per fibre the resultant
  in z of the quadratic-in-U quartic with e*den - num, for I_3 and for I_5 computed on the closed-form residues, equals
  (characteristic polynomial of the data block)^2 exactly.  The symbolic version of this check (`s_closed_sol2.py` v2) was
  killed by its saved PID after 22 CPU-min as superseded (pids_killed.txt); its v1 used a wrong evenness assumption (log kept).
  `s_closed_num.py`: at all four roots z, I_1 = 11/4 and I_3 = the data roots, and r(-z) = s(z).

## Errors, deviations, disclosures
- Seal sec. 4: "I_5 needs the second order in 1/xi" was wrong (sec. 6); the engine was extended to fourth order and
  re-validated before the prediction; the interim message to the lead stated this.
- `t0e_validate.py` v1 ran at t = 2 and reported FAIL at order 6 for the reason in sec. 6 (a gap of the three-term form,
  not an engine fault); v2 uses t = 9/4.  Both logs kept.
- `s_closed_sol2.py` v1: wrong test (sec. 8); `s_closed_debug.py` was killed by the 600 s harness limit (sympy simplify);
  replaced by `s_closed_num.py` and `s_closed_pts.py`.  `f_sing_sol2.py` (non-fraction-free) timed out on {-2,1,4}; `f_sol2.py`
  (sympy Groebner) too; both logs kept.  The first free two-point run (`h_two_sol2.py free`, 9 unknowns) hit its 50-min cap.
- `u_interp.py` was run once with the two-order engine (I_5 fits failed by huge coefficients; those logs were overwritten
  by the re-run with the same file names -- the failure is recorded in the interim message and in `g_engine_t*_v1.json`,
  the two-order charges, which are kept).
- `d1_data_exc3.py` is EXC1's data script with the seal guard renamed; `vir_lib.py` is this seat's VIR1 engine (data
  side only, as in EXC1/EXC2).
- Exit codes: t0b_sigma 0; t0a_vacuum 0, 0 (tamper 0 = fires); t0d_validate 0; t0e_validate 0 (v1 2 at t = 2, explained);
  f_sol2 {-1,1,3} 0, {-1,0,4} 3, {-2,2,3} 3, {-2,1,4} 124 (timeout) then f_sing_ff_sol2 3; f_sing_ff_sol2 {-1,1,3} 0 (control);
  {-1,1,3} tamper 0 (solutions exist but I_3 differs); g_compare_sol2 2
  (one point = half level, as reported) and tamper 0 (fires); h_two_sol2 sym: vdim 4 at 70 generic points, 2 at the two
  exceptional ones, 0 under the tamper; u_interp 0, 0; u_compare real 0, tamper 0, cross 0 (fire); x_exceptional_sol2 0; s_closed_pts 0, 0; h_two_sol2 free2 vdim 12.

## Files
`SEAL_EXC3.md` (+ .sha256, .time), `PREDICTION_EXC3.json` (+ .sha256, .time), `HOLDOUT_OPENED.time`, `frob.py`, `wkb3x.py`,
`wkb3x2.py` (v2, four orders; `wkb3x2_v1_maxu2.py`), `vir_lib.py`, `d1_data_exc3.py`, `d0_vacuum_exc3.py`, `t0a_vacuum.py`,
`t0b_sigma.py`, `t0d_validate.py`, `t0e_validate.py` (+ v1), `f_sol2.py`, `f_sing_sol2.py`, `g_engine.py` (+ v1 json),
`g_compare_sol2.py` (+ v1), `h_two_sol2.py`, `u_one.sh`, `u_points_t*.txt`, `u_interp.py`, `u_compare.py`, `s_sym_sol2.py`,
`s_closed_sol2.py` (+ .sing/.out/.json), `s_closed_num.py`, `s_closed_pts.py`, `frob_ff.py`, `f_sing_ff_sol2.py`, `x_exceptional_sol2.py`, `pids_killed.txt`, `d0_sol2_t*.json`, `d1_sol2_t*.json`,
`d1_sol2_t*_with5.json`, logs, `SHA256SUMS`.
