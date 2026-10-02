# SEAL_EXC2 addendum -- Fable seat, 2026-10-02.  Written after T0, T1 and the one-point T2; before any two-point computation.

RESULTS SO FAR (logs f_template.log, f_sol3_t2_*.log, f_diag.log, g_compare_t2.log; t = 2, xi-exponents (1/2, 1/3))
- CORRECTION to the seal's formula (i): the exponents of the class-U normal form at xi = 0 are +-p PX, +-p pi (I dropped a
  factor p).  The class is unchanged; "PX, pi" in f_sol3.py are these xi-exponents a0, a1.
- T0 PASSES: the code returns r02 = -2z, r01 = -b and EXC1's quadratic on the Sol-1 template (v0 of the script printed FAIL
  because I asked solve() for an overdetermined pair; the conditions themselves were right; both logs kept).
- T1, set S_a = {-1,1,2,4}: indicial conditions give r22 = -4z, r13 = 8z^2, r04 = -8z^3; the six no-log conditions split into
  SEVEN polynomial equations in the seven remaining unknowns (exactly one of them carries E, and it is linear); the system
  is zero-dimensional with FOUR solutions (z = +-1/18, z^2 = 49/5832 at this point), all with r21 = -b.
  So trivial-monodromy apparent singularities DO exist for the fourth-order three-term operator.
- T2 for S_a FAILS as sealed: I_1 = Delta + 1/2 on every solution (r21 = -b shifts PX^2 + pi^2 by 1, i.e. half a level), and
  the two I_3 values (105/64 and 41/64 in monic normalisation) are not the roots of the level-1 data polynomial.
  The value 105/64 equals the VACUUM I_3 at the shifted exponents (a0 + 1/2, a1 - 1/2).
- T1, set S_b = {-1,0,3,4}: NO solution; the obstruction is explicit: the coefficient of E at the resonance (-1, 3) is
  -2592 z^6, which cannot vanish.  Set S_d = {-1,0,2,5}: NO solution (Groebner basis {1}).  Set S_c: still running.

SEALED NOW (new class, not in the original seal): TWO apparent singularities of type S_a,
    L = L_0 + xi sum_{k=1,2} [ R_2^(k) vartheta^2 + R_1^(k) vartheta + R_0^(k) ],   poles at xi = z_1 and xi = z_2, exponents {-1,1,2,4} at each.
COUNT: 20 unknowns; 6 fixed by the indicial conditions; 14 remain against 2 x 7 = 14 equations.
TESTS (t = 2, same momentum point; exact, Groebner / Singular):
  U1  the system has solutions with z_1 != z_2, z_k != 0, and I_1 = Delta + 1 on them.                                55%
  U2  the set of I_3 values over the solutions contains the two roots of the level-1 data polynomial
      e^2 - 89 e/32 + 20209/4096 (i.e. that quadratic divides the eliminant of I_3).                                   40%
  U3  (hold-out, only if U2 holds) I_5 from the oper for the solutions selected by U2 -> PREDICTION_EXC2.json, hash to
      the lead before the Sol-3 I_5 level-1 block is computed.  Needs the first order in z/xi in the WKB.
If U1 fails: no two-point level-1 oper in this class; reported with the inconsistent subsystem.
If U1 holds and U2 fails: reported with the I_3 values found.
Controls: the one-point code path must reproduce the four one-point solutions when the second point's residues are set to 0.
