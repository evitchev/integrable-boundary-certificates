# SEAL_VIR2 -- Fable seat, written BEFORE any machine computation

Stage VIR2 (lead's go of 2026-10-01, after VIR1): exponential screenings of the
cylindrical charges in the two-boson (X, phi) realisation, phi the Coulomb-gas boson of
the Y-sector's Vir_N.  Workspace ~/fable-work/vir2, head only, repo read-only.

## 0. Disclosure: what I know before sealing

Read (no code run for this stage): memory algebra1-wd21alpha, algebra1-wd21-null-screenings,
paperclip-screenings, screening-charges-are-exponents; ansatz note 38(ew) items 1-5 (BOP
step 2/3) and 38(ex) item 51 head; cyl_hidden_kdv_scan docstring (Part B/C); everything of VIR1.

PENCIL WORK DONE BEFORE SEALING (this is not a blind seal for the items marked [hand]):
I worked out by hand the spin-3 conservation conditions for V = e^{aX} Phi_h, Phi_h a Vir_N
primary, against I_3 = A4 x1^4 + A2 x2^2 + C x1^2 T + M (TT)  (M = 4, C = 2 x cyl_P4's
(11)_X(11)_Y coefficient):
 (i)  generic (Verma) Phi_h: a != 0 forces C = M, i.e. the cross coefficient 2 -- false on all
      three solutions.  So a mixed screening needs a DEGENERATE Phi_h.
 (ii) level-2 degenerate Phi_h (L_-2 = kappa L_-1^2 on it), four conditions E1..E4:
        E1: 2M - C a^2 - 2 M h = 0
        E4: C - 2 a^2 kappa (C - M) = 4 A4
        E2: 2 y1 = 4 A4 a^2 - 4 A2 + 2 C h      (y1 = 2C - 2 kappa (C - M), per unit a)
        E3: a^2 y1 + 2 (C - 2 a^2 kappa (C - M)) = 12 A4 a^2 + 2 C h
      and found, identically in t:
        Sol 1: all four hold with h = (3u^2 - 2)/4, a^2 = 1 - u^2/2;
        Sol 3: all four hold with h = (3/u^2 - 2)/4, a^2 = 1/(2u^2);
        Sol 2: E1 and E4 are incompatible for both level-2 fields.
I have NOT examined I_5 or I_7, level-3 degenerate fields, completeness of the solution set,
or anything on the machine.

## 1. What is new in VIR2 (lead's point 1) -- my reading, confirming the lead's

- BOP (38(ew)): V = e^{aX + i b (u.Y)} and its reflections -- exponentials of the LITERAL Y
  bosons.  For N > 1 such a V generates the full (Verma-type) Vir_N module: the level-2
  singular vector is a nonzero Fock state (it contains transverse oscillators).  By (i)
  above that class cannot contain a mixed screening; BOP's negative is consistent and stands.
  BOP's pure-Y orbit e^{sqrt2 u.Y} is "a = 0, h = 1" here: a weight-1 Vir_N primary, which
  commutes with everything in Heis_X (x) Vir_N -- that is why its weight-4 condition was
  "automatic" and its commutant "characterised nothing".
- cyl_hidden_kdv_scan Part B / BOP step 3: pure-X e^{aX} (beta = 0).  Prior art; a control here.
- ALGEBRA1 (item 51, 122-132): LZ's three screenings act in the pillow's own three-boson
  field content; the comparison with the cylinder was of vacuum eigenvalues and commutant
  characters, never an action on the cylindrical densities.
- NEW here: V = e^{aX + beta phi} = e^{aX} x (Vir_N primary in its Coulomb-gas, i.e.
  degenerate-capable, realisation), acting on the cylindrical densities themselves, which
  VIR1 placed in Heis_X (x) Vir_N.  For generic N this V is not a local functional of Y.
  At N = 1 (rho = 0, phi = Y) the class reduces to BOP's, so BOP's N = 1 solutions are a
  control.

## 2. Conventions (sealed)

- Propagators: dX(z) dX(w) ~ 1/(z-w)^2, d phi(z) d phi(w) ~ 1/(z-w)^2 (the engine's; X and phi
  commute).  X exactly as in cyl_P4: no background charge, densities even in X.
- T_phi = (1/2)(d phi)^2 + rho d^2 phi,  c = 1 - 12 rho^2 = N.
- Parameter u with u^2 = (t-1)/(t+1) (the record's beta^2), t = (1+u^2)/(1-u^2),
  N = 13 - 6(u^2 + u^-2), s = 6(u^2 - u^-2).
  SIGN: rho = (u - 1/u)/sqrt2.  ((u, rho) -> (-u, -rho) is phi -> -phi.)
- h(e^{beta phi}) = beta^2/2 - rho beta.  Virasoro screenings (h = 1): beta_+ = sqrt2 u,
  beta_- = -sqrt2/u.
- Y-sector -> phi: the HV decomposition of VIR1 (density = sum of X-monomial (x) Vir state, mod d),
  each Vir state L_{-n1}..L_{-nk}|0> rebuilt with T_phi.  FAIL-CLOSED GATE: the image of
  (11)(11) - 2(22) is 4 (T_phi T_phi) mod d, and the two-boson images of I_3 and I_5 commute.
- Screening test: R = Res_z P(z) V(w) / V(w) must lie in the image of the twisted derivative
  d + a dX + beta d phi; decided by a unique normal form (no bare d phi letter; bare dX for the
  beta = 0 control).
- A screening COUNTS only if conserved against I_3 AND I_5; I_7 is the hold-out.

## 3. PREDICTIONS

S1 [hand for I_3; I_5, I_7 blind].  Solution 1, every t:
     V_+- = exp(+-a X + beta phi),  beta = -beta_+/2 = -u/sqrt2,  a^2 = 1 - u^2/2 = (t+3)/(2(t+1)).
   Conserved against I_3 (hand), and I PREDICT against I_5 and I_7 (confidence 80%).
   a^2 + beta^2 = 1 for all t.
S3 [hand for I_3; I_5, I_7 blind].  Solution 3, every t:
     V_+- = exp(+-a X + beta phi),  beta = -beta_-/2 = 1/(sqrt2 u),  a^2 = 1/(2u^2) = (t+1)/(2(t-1)).
   Conserved against I_3 (hand); I PREDICT against I_5 and I_7 (confidence 80%).
S2 [blind].  Solution 2: no level-2 (Phi_(1,2)-type) solution [hand].  PREDICTION, confidence
   50%: a mixed pair exists with a LEVEL-3 degenerate phi-exponent, beta = -beta_+ or -beta_-,
   with a^2 -> -1 at N = 1 (BOP's N = 1 control for Solution 2), conserved against I_3 and I_5.
S-complete [blind, confidence 60%].  At generic t the mixed (a != 0, beta != 0) solution set of
   {I_3, I_5} is exactly the pair above on Sols 1 and 3 (plus whatever S2 gives); every other
   I_3-only solution (e.g. the continuation of BOP's N = 1 branch (2,2)) is killed by I_5.
S-Gram [follows from S1, S3 if they hold; printed now so it cannot be adjusted].
   v = (X-component, phi-component).
   Sol 1: {v_phi = (0, sqrt2 u), v_+, v_-}, v_+- = (+-a, -u/sqrt2):  v_phi + v_+ + v_- = 0,
          Gram [[2u^2, -u^2, -u^2], [-u^2, 1, u^2-1], [-u^2, u^2-1, 1]].
   Sol 3: {v_phi' = (0, -sqrt2/u), v_+, v_-}, v_+- = (+-1, 1)/(sqrt2 u):  v_phi' + v_+ + v_- = 0,
          Gram (1/u^2) [[2, -1, -1], [-1, 1, 0], [-1, 0, 1]].
   No algebra is named in this seal.
S-exponents [hand, Sol 1 only].  With the record's (A, B, C) = (eps, -a_X, -a_Y):
   Sol 1: B = -a^2/a^2 = -1,  C = -beta^2/a^2,  A = (v_+ . v_-)/a^2.  For Sol 3 I predict
   B = C corresponds to a^2 = beta^2 (true of S3) and leave the normalisation to the computation.

FAIL FAST (as agreed): no mixed solution conserved against BOTH I_3 and I_5 at generic t on
any solution.  If S1 or S3 fails at I_5, that solution is a FAIL regardless of the I_3 result,
and is reported as "I_3-only".

## 4. Controls

C0  gates of sec. 2 (Virasoro algebra of T_phi with c = 1 - 12 rho^2; the (T_phi T_phi) gate;
    [I_3, I_5] = 0 in the two-boson image; two-boson I_3 has no odd-in-X term).
C1  a = 0: beta = beta_+ and beta_- pass against I_3, I_5, I_7 on all three solutions at every
    tested t; a generic beta does not.
C2  prior art: Sol 1 at (28,-9): (a^2, beta) = (2, 0) passes I_3 and I_5; at a generic fibre
    it does not.  N = 1: (a^2, beta^2) = (1/2, 1/2) on Sols 1, 3 (BOP's control).
C3  tamper A: rho -> -rho with the densities unchanged in (X, T): the sealed (a^2, beta) must
    FAIL and its mirror (a^2, -beta) pass.
    tamper B: rho^2 from N + 1 instead of N: the two-boson images of I_3, I_5 must stop
    commuting and the sealed pairs must fail.
C4  plant: a density that the sealed V does not conserve (I_5 + (1/1000) x a generic HV element)
    is flagged; scoring tamper after the hash guard.
C5  vacuity: a generic (a, beta) fails every condition set; count of conditions printed.
"Identically in t": symbolic in u where feasible, otherwise exact at >= 20 rational u with
both signs of t, none special.

Exit codes: 0 = S1 and S3 hold against I_3 and I_5 (and the I_7 hold-out) with all controls
firing; 2 = fail-fast fired or S1/S3 refuted at I_5 or I_7; 1 = a gate or control failed.
S2 and S-complete are reported either way and do not set the exit code.
