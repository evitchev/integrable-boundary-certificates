# SEAL_JOINT1 -- Fable (secondary) seat, 2026-10-10.  Joint (state-by-state) check of the excited-state operators.
Commissioned by the lead (JOINT1_brief_2026-10-10.md, sha256 4f485d93...), prompted by ChatGPT's review: equal characteristic polynomials of
I_3 and of I_5 separately do not show that the operator's (I_3, I_5) values pair up on the same states.  Written before any computation of this
stage.  Archives read-only; nothing archived is recomputed.

## Cells (all level ONE; level two has no archived I_5 block -- EXC1's and EXC4's level-two comparisons are I_3 only -- so it is out of scope)
  C1  EXC1  Sol 1  t = 2      C2  EXC1  Sol 1  t = 9/4
  C3  EXC2  Sol 3  t = 2      C4  EXC2  Sol 3  t = 9/4
  C5  EXC3  Sol 2  t = 2      C6  EXC3  Sol 2  t = 9/4
Data: the archived 2 x 2 singlet blocks of I_3 and I_5 in bare (P, Q) (d1_sol<s>_t<t>_with5.json of EXC1/2/3, symbolic in P, Q; the same files
the comparisons used), each normalised by the P^(2K) coefficient of its vacuum eigenvalue, as in the comparisons.
Operator: the archived closed forms -- EXC1: the level-1 condition 2M z^2 + alpha(M-1) z - 2 lam^2 + (M-1)^2/2 = 0 with R_4(alpha, lam, z), R_6
(o1_t<t>_n1.json) and EXC1's dictionary; EXC2: U = u^2 root of F (s_closed.json), I_3 linear and I_5 quadratic in U (the s_closed2 maps);
EXC3: z root of c7 z^4 + c5 z^2 + c3 (s_closed_sol2.json) with I_3(z), I_5(z) rational in z (the s_closed_pts maps).

## Method
1. Commutation: [M_3, M_5] = 0 exactly (symbolic in P, Q) -- re-checked, not taken from the stored flag.
2. Degeneracy: the discriminant of the I_3 block's characteristic polynomial as a polynomial in (P, Q) (it is not identically zero; the momentum
   points used below avoid its zeros).
3. JOINT INVARIANT, symbolic in (P, Q) where the closed form allows (C1-C4), at rational points for C5-C6 (the EXC3 residues are rational in z):
   for a 2 x 2 block the multiset of pairs {(a1, b1), (a2, b2)} is determined by tr M_3, det M_3, tr M_5, det M_5 (the marginals, already
   compared in the archives) together with the MIXED TRACE  tr(M_3 M_5) = a1 b1 + a2 b2;  the crossed pairing gives a1 b2 + a2 b1, which
   differs by (a1 - a2)(b1 - b2) != 0 off the degeneracy locus.  Operator side: sum over the two selected solutions of I_3 I_5, a symmetric
   function, computed exactly (companion matrix of the quadratic in U for EXC2; trace in Q[z]/(quartic) for EXC3 and in Q[z]/(quadratic) for
   EXC1).  Comparison exact; also tr(M_3^2 M_5) and tr(M_3 M_5^2) for completeness.
4. EXPLICIT PAIRS at one rational momentum point per cell: the I_3 eigenvectors of the data block over the quadratic splitting field, the
   I_5 eigenvalue on each eigenvector; the operator's (I_3(z_i), I_5(z_i)) at the two solutions; equality as sets of pairs in that field.
Tampers (after the hash guard; must FAIL): (i) the crossed pairing {(a1, b2), (a2, b1)} against the data's mixed trace / pairs;
(ii) the operator pairs at a momentum shifted by +1/10 in P (resp. in P_X) against the data at the unshifted point.
Control: on each cell the data's own crossed invariant a1 b2 + a2 b1 is reported with (a1 - a2)(b1 - b2) so that the tamper is seen to be
non-vacuous.
Momentum points for 4: bare (P, Q) = (1/2, 1/3) and (2/5, 1/4) (converted to each stage's variables).

## Priors
All six cells pair correctly (mixed traces identical, explicit pairs equal): 85%.  Each operator was found by matching a one-parameter
family (the quadratic / the quartic) whose two solutions carry both charges jointly, so a crossed match would need a coincidence; still,
no joint check exists and the prior is not 99%.  Both tampers fire on every cell: 95%.

## Conventions
Exact arithmetic (sympy; Gamma atoms do not arise: all data are rational functions of P, Q); runs capped (ulimit -v 16 GB); full logs; SHA256SUMS;
NOTE.md with the table cell / block size / commute / degenerate / pairs match / tamper fires.  Workspace ~/fable-work/joint1.
