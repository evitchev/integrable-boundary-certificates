# SEAL_VIR7 -- Fable seat, 2026-10-01.  Solution 1: a sealed family of Gamma-symbol opers, scanned at first quantum order.
Commissioned by the lead (option (c)): "Seal ONE natural class consistent with your Sol 1 top law ... test it zero- or
minimal-parameter at one fibre, with hash-registered hold-outs"; do not duplicate SCAL2 (Suzuki Lambda + kappa sqrt(X) U +
rational V) nor cc's SOL12 class (item 301).  Written before any comparison with Solution-1 data in this stage.

Notation: k = 2(t-1)/(3-t), p^2 = k/(k+1); T = s - theta; a "string" G_a(T - e) = sigma^a R_a((T - e)/sigma),
R_a(u) = Gamma(u + (a+1)/2)/Gamma(u - (a-1)/2), sigma = nM (centred at e; ~ (T-e)^a (1 + c_1(a) sigma^2/(T-e)^2 + ...),
c_1(a) = -a(a^2-1)/24).  K = spin index.  Record exponents of Sol 1: a_X = 1, a_Y = (t-1)/(t+3) = k/(3k+4).

## 0. Disclosure (what was done before this seal)
(a) Read: items 290, 293, 301 of the ansatz note and the shared-memory notes on SOL12/SOL2F and S1-scope.  I know from them:
    Suzuki + centrifugal with M_S = (t+3)/(t-1) reproduces Sol 1's top and FAILS below it; cc's class G with beta = 1
    (symbol (T^2 - l1^2) R_alpha((T -+ l0)/delta) R_(-2 alpha)(T/delta), alpha = a_Y) has no (M, C) fit at spin 3, t = 2.
(b) Built and validated on KNOWN operators only (`g0_validate.py`, exit 0, 0/2013): formula G1, the loss-1 layer of any
    operator H(T) psi = x^n (x^(nM) - E) psi with H_0 = T^n phi(1/T), phi = prod (1 - e_i z)^(alpha_i), z = 1/T, nu = (2K-1)/n:
        Q_2K = [z^2K] phi^nu + (nu c2/24)[z^(2K-2)] phi^nu {(n-1) - nu (n L1 - L2 - L1^2)} + nu [z^(2K-2)] phi^nu b(z),
        L1 = D log phi, L2 = D L1, D = z d/dz, c2 = n^2 (M+1), b(z) = sum over strings c_1(a) sigma^2/(1 - e z)^2.
    It reproduces Codex's loss-1 Ward amplitudes of Sol 3 (item 298 operator) AND of Sol 2 (cc's S2(k): strings of length k
    at +-l0, n = 2k+3, M = 1/k), K = 2..20, five fibres each, with l0^2 = c2 PX^2/alpha, l1^2 = c2 pi^2/beta (alpha, beta =
    total powers of the X and Y factors) and record variables X = -PX^2, Y = rho^2 - pi^2.
(c) Hand derivation (not yet machine-checked), "class U": the three-term equation
        [ P(vartheta) + xi^(1/2) L(vartheta) xi^(1/2) + E xi^b ] phi = 0,   vartheta = xi d/dxi,
    is equivalent, by phi-hat = g eta-hat with g(nu) = prod_(roots r of L) (1-b)^(nu/(1-b)) Gamma((nu - r - 1/2)/(1-b)), to the
    two-term Gamma-symbol equation with symbol P(T) prod_r G_a(T - r), a = b/(1-b), M = 1/a, n = deg P + a deg L.
    In this form: Sol 3 = (P, L) = ((T^2-l0^2)(T^2-l1^2), T), Sol 2 = (T(T^2-l1^2), T^2-l0^2), both with b = p^2;
    Suzuki's ODE with a centrifugal term (item 293) = (T^2 - l1^2, T - x0) with b = 2/(M_S+1), a = 2/(M_S-1) -- so the class
    of item 293 is already a Gamma-symbol operator with ONE string, centred at an exponent linear in PX.
(d) One hand check with Sol-1's TOP LAW only (K = 2, t = 2): the single-string classical symbol (T^2 - Y)(T - x0)^(1/2)
    reproduces the three top ratios.  No Sol-1 lower-layer number has been looked at in this stage; the Sol-1 tables at
    spins 11, 13 (anchor11, anchor13) have never been opened by this seat (W = 6, 8, 10 were, in SHEET1/VIR1).

## 1. The class (ONE class, a finite list of members)
CLASS U for Sol 1: symbol P(T) x prod over roots r of L of G_a(T - r)^(+-1) (numerator roots: strings; denominator roots:
inverse strings), ALL strings of the one length a = 1/M, potential x^n (x^(nM) - E), n = deg P + a deg L, with
    P = (T^2 - l1^2)^(beta_P),      L = (T^2 - l0^2)(T^2 - l1^2)^(beta_L)/T^2
(the X pair only in L, a double inverse string at 0), the classical symbol P L^a being REQUIRED to reproduce the Sol-1
top law identically in t.  With a = k/(uk + v) this leaves exactly (alpha_L = 1):
    U1_u, u = 0, 1, 2, 3:  beta_P = 1, beta_L = 3 - u, a = k/(uk+4)       [u = 3: n = 2, M = (t+3)/(t-1) = Suzuki's M]
    U2_u, u = 0, 1:        beta_P = 2, beta_L = 3 - 2u, a = k/(uk+2)      [u = 1: a = (t-1)/2 = Suzuki's string length]
    U4:                    beta_P = 4, beta_L = 3, a = k.
Repeated roots are taken literally (a unit root of P times a centred string from L at the same exponent; P's repeated
roots unshifted).  Dictionary: l0^2 = C PX^2/alpha, l1^2 = C pi^2/beta, alpha = a, beta = beta_P + a beta_L; ONE scale C
per member (the rule validated on Sols 2, 3 is C = c2 = n^2(M+1); reported separately as the zero-parameter version).
All members are X-even, so they can only give the odd-spin charges; Sol 1's X-odd even-spin charges are outside this
class by construction (stated limitation).
Controls (not members): S = the single-string form of item 293's ODE (P = T^2 - l1^2, L = T - l0, a = (t-1)/2; X scale
fixed by the top) -- must FAIL below the top, as item 293 found; Gcc = cc's beta = 1 symbol with R_(-2 alpha) and
M = 1/alpha -- must FAIL (item 301).

## 2. Tests
T0 (gate, each member): G1's top == the Sol-1 top law K_K(p, q; r), K = 2..20, t = 2, 9/4, 3/2, 4, -2.
T1 (training, first quantum order): G1's loss-1 layer, in record variables, against Codex's Sol-1 loss-1 shift
   amplitudes (item 185 archive, AMPLITUDES_sol1_loss1.json), K = 2..20 at the five fibres.  1/C is fitted from ONE
   coefficient at K = 2 and every other coefficient must then agree exactly.  PASS(member) = 0 mismatches at all five
   fibres.  Also reported: whether the fitted C equals c2.
T2 (hold-out, only for members that pass T1): full polynomials at spins 3..13 from the Mellin-WKB engine (even symbols),
   -> PREDICTION_VIR7.json, sha256 to the lead BEFORE any comparison; then against vev_sol1_w6/8/10 and, after the
   lead's "registered", against the unopened anchor11/anchor13 Sol-1 profiles.  Losses 2..4 also against the item-180
   Sol-1 amplitudes.
If NO member passes T1: clean negative for class U at first quantum order (with the witness: first failing K and
coefficient per member), plus, labelled post-hoc, the residual of the best member expressed as a function b(z).

## 3. Predictions (sealed)
S fails T1: 90%.  Gcc fails T1: 85%.  At least one member of class U passes T1: 35%
(most likely U1_3 or U2_1, 12% each; U1_2 8%; the rest 3% together).  If a member passes T1, it passes T2: 60%.

## 4. Controls
G1 validation on Sols 2, 3 (done, sec. 0b).  Machine check of the class-U equivalence at the BLZ point (alpha-term off):
(T^2 - l^2) G_a(T), n = 2 + a, M = 1/a against the Schrodinger equation with M_S = 1 + 2/a, both by my engine, spins 1..7.
Plant: T1 run on Sol 2 with its true operator must PASS and with string length k + 1/10 must FAIL.  Scoring tamper
(after the seal-hash guard): Codex A_(1,1) + 1/1000 must turn a PASS into a FAIL (run on the Sol-2 plant if no member passes).
Exit codes: h1_scan.py 0 = a member passes T1 and controls behave; 3 = no member passes, controls behave; 1 = a control
misbehaves.  No nsimplify, no lambdify, exact rationals.
