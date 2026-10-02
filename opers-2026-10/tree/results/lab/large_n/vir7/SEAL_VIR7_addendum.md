# SEAL_VIR7 addendum -- Fable seat, 2026-10-01.  Written after h1_scan.py (T0/T1) and before any engine run or any
# comparison beyond the loss-1 layer.

WHAT HAPPENED (logs: h1_scan_real.log, h1_scan_real_v0_Scontrol_scale.log, h2_inverse.log, h2_inverse_v0.log)
1. All seven members of class U pass T0 (top) and FAIL T1 (loss-1 layer).  Class U as sealed is excluded at first quantum
   order.  cc's symbol (control Gcc) fails T1, as expected.  The Sol-2 plant passes and its tamper fails.
2. The control S -- the single-string symbol  H(T) = (T^2 - l1^2) G_a(T - l0),  a = (t-1)/2,  n = 2 + a = (t+3)/2,
   M = 1/a,  sigma = nM = (t+3)/(t-1), which by the class-U derivation is item 293's ODE (Suzuki + centrifugal) --
   PASSES T1: 0 mismatches at K = 2..20 at all five fibres (166..208 coefficients each), with ONE fitted scale that comes
   out C = n^2(M+1)/2 at every fibre, i.e.  l1 = n pi/p,  l0 = n PX/(p sqrt(a_Y))  (the same rule l1 = n pi/p as Sols 2, 3,
   where M + 1 = 1/p^2; here M + 1 = 2/p^2).  My sealed prediction "S fails T1 (90%)" is REFUTED.
   DISCLOSURE: in the first run (v0) the S control returned before testing anything, because I had restricted its X scale
   to a finite set that did not contain the crossed-curve value l0^2 : l1^2 = PX^2/a_Y : pi^2; v1 uses that value.  This
   repair was made after seeing only that v0 had not tested S; member verdicts are unchanged (diffed).
3. Post-hoc inverse problem (even base, n = 2): the unique b(z) reproducing the loss-1 layer is
   b/c2 = (1 + 2 alpha)/24 + ((1 - alpha)/12)(1+u)/(1-u)^2, u = l0^2 z^2, alpha = a_Y: the pair term is exactly a string pair
   of length a_Y with M = Suzuki's M; the constant matches none of the sealed frozen blocks.  (v0 of h2_inverse.py
   mis-indexed the series by one; v1 starts at the spin-1 charge.)
I do not understand why item 293 reports that the lower layers fail for Suzuki + centrifugal; a wrong overall scale of the
momenta (which the top cannot see) would produce exactly that, but I have not checked SUPER3's dictionary.

SEALED NOW (hold-out T2 applied to S, zero parameters)
Operator S1(t):  (T^2 - l1^2) sigma^a R_a((T - l0)/sigma) psi = x^n (x^(n/a) - E) psi  [signs as in the engine],
   a = (t-1)/2, n = (t+3)/2, sigma = n/a, l1 = n pi/p, l0 = n PX/(p sqrt(a_Y)), p^2 = 2(t-1)/(t+1), a_Y = (t-1)/(t+3);
   record variables X = -PX^2, Y = rho^2 - pi^2.
Engine: a copy of mellin_wkb.py generalised to symbols with odd powers of 1/T (the symbol is not even in T), validated
   BEFORE any Sol-1 comparison on (i) the Sol-3 operator against the old engine and (ii) the BLZ point of the class-U
   equivalence: (T^2 - l^2) G_a(T), n = 2 + a, M = 1/a against the Schrodinger equation with M_S = 1 + 2/a.
T2a (not blind): losses 1..4 at K = 2..7 against Codex's Sol-1 Ward amplitudes (items 180/185), and the full polynomials
   at spins 5, 7, 9 against vev_sol1_w6/8/10 (tables this seat has opened before), fibres t = 9/4, 3/2, 4, -2 and 7/3.
T2b (BLIND): spins 11 and 13 at the same fibres -> PREDICTION_VIR7.json, sha256 to the lead; the anchor11/anchor13
   Sol-1 profiles (never opened by this seat) are compared only after the lead's "registered" -- or by the lead.
Also reported: the even-spin (odd-order) charges of S1, which must be nonzero and odd in PX (no tables exist).
PREDICTION: T2a passes at every loss: 55%.  If T2a passes, T2b passes: 85%.  If T2a fails: first failing (K, loss) reported.
Controls: tamper a -> a + 1/10 and tamper l1 scale x 11/10 must fail T2a.
