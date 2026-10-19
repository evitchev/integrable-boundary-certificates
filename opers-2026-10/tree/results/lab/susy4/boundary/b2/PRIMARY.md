# B2 PRIMARY (Fable secondary seat, 2026-10-13): is Solution 1 the U(1)-coset image of an N = 2 KdV tower?

PREREG: results/lab/susy4/boundary/B2_PREREG.md, commit c24a7af, sha256 9ece3fda2966938a219b7e1bbd9af5b653a84130e8ab06da97f98aed8cba59fd
(verified on the tree and on the blob).  Erratum (lead, item 348): Zeitlin's quantum vertices are his (22).  Rulings on the
seat's reading points R1-R4 (lead, 2026-10-13): all accepted.  Workspace <workspace>/b2/ (NOTES.md = timeline; logs/ = full
stdout of every run; every run under ulimit -v 8000000 and timeout).  Repository, engines and sources read-only; nothing from
sources/ executed (PDF pages rendered to read KL (21)-(22), Mao-Liu (1), Tian-Liu (1)-(2)).

## VERDICT (PREREG sec. 4)

**FAIL (both towers) for the decisive statement P-II-Z / P-II-K as ruled in R3 (density level, no products), at spins 3, 4, 5,
at all three primary fibres t = 2, 7/2, 11/3.**  Certified negative in the PREREG's words: *Solution 1's tower is not the
U(1)-coset image (H~-vacuum projection P_H I P_H) of Z or K at spins 3 and 5 at the fibres tested.*  What does hold, exactly
and at every fibre: Level I (the screening data of Z coincide with Solution 1's: P-I-Z holds; P-I-K holds); the spin-2
projection of Z lands exactly on Solution 1's X-odd weight-3 charge with c_2(t) = (t-1)/(t+3); every piece of every projection
that lands at weight <= 6 (the T and I_2 levels) and every constant lands on the corresponding Solution 1 charge; the pieces of
weight >= 8 never do.  The coset image is characterised in sec. 8.  Validations V1-V3 pass; V2 identifies Tower Z as the
a = 4 N = 2 KdV hierarchy and Tower K as a = -2 (classical limits, unique fits).  Nothing was refitted; the controls behave as
required.  No PASS branch is claimed at any spin >= 3 under any reading (sec. 7 records a post-hoc "products" reading that
rescues spin 3 only, automatically).

## 1. Dictionary (PREREG sec. 1, last bullet) -- MATCH, no STOP

p^2 = 2(t-1)/(t+1) and rho = p/2 - 1/p (opers paper, Setting); n_FL = -(t-1)/(t+1) = -p^2/2 (opers paper, "Relation to quantum
AKNS", Assis (2.2)-(2.4); item 337).  FL hep-th/0510271 after (83): parafermions at non-integer k = -n-2, so c_W = 2(k-1)/(k+2)
= 2 + 6/n (= KL19 1910.05947 (3.16)).  FNP hep-th/0406017 (2.10)-(2.12): T = T_coset - 1/2 (d phi)^2, c_{N=2} = 3k/(k-2) =
c_coset + 1.  Hence c_{N=2} = 3 + 6/n = 3 - 12/p^2.  Sources' parameters: KL (15) c = 3(1 + 4/b^2) gives b^2 = -p^2; Zeitlin's
beta^2 = b^2 (his (31) exchange exponent -beta^2/2 = KL's S1 S2 exponent -b^2/2).  b, beta imaginary at every fibre.
Fibres: t = 2, 7/2, 11/3: s := p^2 = 2/3, 10/9, 8/7, c_{N=2} = -15, -39/5, -15/2; hold-out t = 13/5: s = 8/9, c = -21/2.

## 2. Frame and V1 (n2frame.py; logs n2frame_v1.log, n2frame_weights.log, n2frame_vertices.log)

Two N = 1 superfields (X_j, psi_j), diagonal metric eta = (s, -s) (engine +log convention), null momenta m_1 = (1/2, 1/2),
m_2 = (1/2, -1/2) (m_i.m_i = 0, m_1.m_2 = s/2 = -b^2/2), m_3 = -m_1, m_4 = -m_2 (KL S1..S4 = (m.psi) e^{m.X}).  Background
charge forced by h(S1) = h(S2) = 1: T = (1/s)[1/2 (dX_0)^2 - 1/2 (dX_1)^2 - 1/2 psi_0 dpsi_0 + 1/2 psi_1 dpsi_1 - d^2 X_0],
G = (1/s)[psi_0 dX_0 - psi_1 dX_1 - 2 dpsi_0].  The second supercurrent, derived by an ansatz solved against the N = 2 OPEs,
is G' = i g' with g' = (1/s)[-psi_0 dX_1 + psi_1 dX_0 - 2 dpsi_1] real (indefinite frame); G+- = (G -+ g')/sqrt2 are real;
J = (G+ G-)_2/2 = (psi_0 psi_1 - 2 dX_1)/s.  **V1: every N = 2 OPE (T, G+, G-, J; KL (9) form; poles 1-5) closes symbolically
in s with c = 3 - 12/s (36/36 checks)**; J G+ ~ +G+, J J ~ (c/3)/z^2.
Weights (symbolic s): S1, S2: h = 1 (primary); S3, S4: h = 0; e^{-X_0}: s/2 - 1; Zeitlin's third vertex W0 = G_{-1/2}[psi_1
e^{-X_0}] = (dX_1 - psi_0 psi_1) e^{-X_0} (item 348 construction; the by-hand ordering gave the same): h = s/2 (R2: p^2/2, not
1); J-charge 0 (J commutes with W0).  Commutation at s = 2/3, 10/9, 8/7, 8/9, 3, 5/7: T, G, g' commute with all five vertices;
J with S1, S2, W0 only (S3, S4 carry J-charge -+2: R1).  W0 built from g' instead of G also commutes with everything.

## 3. Towers (towers.py; logs towers_*.log; json towers_*.json)

Commutants mod d in the full two-superfield monomial space (both Z2 x Z2 parities), s = 2/3, 10/9, 8/7, 3:
| w2 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|---|---|---|
| Z = {S1, S2, W0} | 1 (J) | 2 (G, g') | 1 (T) | 0 | 1 (spin 2) | 0 | 1 (I_3) | 1 (spin 4) | 1 (I_5) |
| K = {S1, S2, S3, S4} | 0 | 2 | 1 | 0 | 0 | 0 | 1 | 0 | 1 |
Identical at all four s.  Every charge commutes with G_0, g'_0, J_0, T (diag3_susy_invariance.log; imposing G-invariance
changes no dimension).  Z has a charge at every integer spin (the a = 4 hierarchy), K at odd spins only.  Z's I_3, I_5 differ
from K's mod d.  Cost: w2 = 8: 0.1 s; w2 = 12: Z 3.7 s, K 11 s; peak RSS 68 MB.

## 4. V3 (v3.py; logs v3_s*.log) -- PASS, strict

K's spin-3 density; Delta(P) and q(P) from OUR T and J by the same evaluator (no use of KL (32)); hw value of KL (22):
-(3/4 + c/8) q^2 + (12 - c) Delta^2 + 3 Delta q^2 [+ 6 Delta q^2 if the sum :L_k J_n J_m: includes k = n = m = 0] + A Delta + B,
A = [576(5+2c) - 48 b^2 (39+20c)]/(1152 b^2), B = [-12c(13+2c) + b^2(123c + 22c^2)]/(1152 b^2).  Result: **exact
proportionality, including the L_0-linear and constant terms, with the reading that the sum EXCLUDES k = n = m = 0** (alpha =
-32, -12000, -3584/3, -288 at s = 2/3, 10/9, 8/7, 3); the other reading fails at every fibre and cannot be rescued by
beta Delta + gamma.  This settles the reading needed by B1.  KL (21) has hw value 0 (nothing to compare).

## 5. V2 (classical.py, diag_flow.py; logs classical_Z.log, classical_K.log, diag_flow_Z.log) -- Z is a = 4, K is a = -2

Classical commutants (one-contraction brackets at eta-hat = (1, -1); the quantum metric is s eta-hat): Z dims 1, 2, 1, 0, 1,
0, 1 and K 0, 2, 1, 0, 0, 0, 1 at w2 = 2..8 (= the quantum ones); J^, G^, g'^, T^ in Z's; G^, g'^, T^ in K's.  The quantum
weight-4 density (normalised (dX_0)^4 = 1) at s = 1/j, j = 2..31, each coefficient reconstructed as a rational function of s
(25 points, verified on 5) and evaluated at s = 0, EQUALS the classical commutant density (both towers).
N = 2 superfield Phi = J^ + theta_1(b_1 G+ + b_2 G-) + theta_2(g_1 G+ + g_2 G-) + theta_1 theta_2(d T^ + e dJ^), fixed by
requiring {oint G^, .} and {oint g'^, .} to act as k_i Q_i (Q_i = d/dtheta_i - theta_i d_x; the graded bracket with an odd
charge carries (-1)^|theta| on the components): four solutions, e = 0, d = +-2i (the covariant derivatives D_i in place of Q_i
give structures that fit no flow, as they should).  Flow: tau {oint D_cl, Phi} versus Mao-Liu 1705.03997 (1)
phi_t = -phi_xxx + 3(phi D_1D_2 phi)_x + (a-1)/2 (D_1D_2 phi^2)_x + 3a phi^2 phi_x with phi = alpha Phi (unknowns a, tau,
alpha; 174 component equations): **Z: a = 4, alpha = -+i/4, tau = 1/8 (unique, all four structures); K: a = -2, alpha = -+i/4,
tau = 1/8 (unique).**  Z's spin-2 charge generates Tian-Liu 1108.2110 (2), Phi_tau = 1/2 (D_1D_2 Phi_y) + 2 Phi Phi_y, the first
flow of the a = 4 hierarchy (alpha = -i/4, tau = i/16; 54 equations).  So Tower Z = Zeitlin's sl^(1)(2|1) tower = the a = 4
("N = 2 NLS") hierarchy, as Krivonos-Sorin's classical pointer says; Tower K = KL's sl(2|2)^(1) tower = a = -2 (Inami-Kanno),
VERIFYING the PRIOR_ART item marked UNVERIFIED.  Zeitlin prints no a; Delduc-Gallot call the sl(2|1) case the N = 2 NLS
hierarchy; the a-value is this computation's.

## 6. Level I (level1.py; log level1.log) -- P-I-Z holds, P-I-K holds; background charges differ (R4)

Bosonisation psi^+- = psi_0 +- psi_1 (null, <psi^+ psi^-> = 2s) = sqrt(2s) e^{+-X_3}, <X_3 X_3> = +log.  Momenta in
Y = (X_0, X_1, X_3), metric (s, -s, 1): S1 (1/2, 1/2, 1), S2 (1/2, -1/2, -1), S3 (-1/2, -1/2, 1), S4 (-1/2, 1/2, -1), W0 (-1, 0, 0).
J -> a.dY, a = (0, -2/s, -1), a.a = (s-4)/s = c/3 (the J J double pole); H~ = a.Y; v = -X_1 - 2X_3 (v.v = 4 - s = 4 kappa^2);
u = X_0.  (u, v, H~) coordinates: S1 (1/2, -1/2, 0), S2 (1/2, 1/2, 0), W0 (-1, 0, 0) [all J-neutral]; S3 (-1/2, (s+4)/(2(s-4)),
-2s/(s-4)), S4 (-1/2, -(s+4)/(2(s-4)), 2s/(s-4)), J-charges -2, +2 (information, no verdict weight).
Gram(S1, S2, W0) = [[1, s/2 - 1, -s/2], [s/2 - 1, 1, -s/2], [-s/2, -s/2, s]] = Gram(V+, V-, e^{p phi}) of Solution 1
(beta = -p/2, kappa^2 = 1 - p^2/4, metric (1, 1)) IDENTICALLY in s.  Isometry u = -p phi, v = -2 kappa X maps S1 -> V+,
S2 -> V-, W0 -> e^{p phi} exactly (not to the dual e^{-2phi/p}: Gram [1, 1, -2, 4/s]).  Background: the coset inherits rho' = 1/p
(c_uv = 2 - 12/p^2 = FL's 2 + 6/n; V+- weight 1, e^{p phi} weight s/2 - 1), the opers paper's Vir_N has rho = p/2 - 1/p (V+- weight
s/4, e^{p phi} weight 1); the difference (p/2 - 2/p) d^2 phi is a total derivative.  K's J-orthogonal system is the pair alone
(two nodes) versus Solution 1's three: P-I-K.  Disclosure: the Gram identity and the isometry were a hand computation (R4) done
before the Level I run; the run reproduced them symbolically.

## 7. Level II (coset.py, sol1.py, level2.py, extras.py, products.py, validate2.py, coset_image.py; logs level2_*.log, ...)

**Objects.**  Projection = to the psi^+- basis -> J-charge decomposition (the charge -+2 parts are single monomials and total
derivatives; the neutral part is kept) -> bosonisation (bilinears by the Schur-polynomial formula; longer neutral Wick words by
the vertex-algebra recursion :B M': = NP(B, M') - corrections) -> rotation X_1 = s(v - 2H~)/(4-s), X_3 = (sH~ - 2v)/(4-s) ->
H~-vacuum expectation at H~-momentum q (cylinder Wick recipe of super_cyl_vev, metric g_H = (s-4)/s per contraction, unpaired
dH~ -> q).  The result is a sum of (u, v) densities labelled by (q-power k, weight w): the H~ contractions produce cylinder
constants, so the q = 0 part has pieces at weights w2, w2 - 4, ... (the PREREG's "lower charges").  Validation (validate2.py, all
three fibres): (i) every neutral fermion word of w2 <= 8 has the same NS vacuum eigenvalue as its bosonic image at X_3-momentum 0
(registered fermionic evaluator versus bosonic one), and the recursion is split-independent; (ii) END TO END: for every charge of
Z (w2 = 4, 6, 8, 10, 12) and K (4, 8, 12), the eigenvalue on |P_0, P_1> x (NS vacuum) by the registered evaluator on the ORIGINAL
two-superfield density equals sum_k q^k vev(D_k)(P_u = P_0, P_v = -P_1) at q = -2P_1/s, exactly.  (iii) J -> dH~ exactly;
the projection of T is T_uv = (du)^2/(2s) + (dv)^2/(2(4-s)) - d^2u/s + constant.
Solution 1's side: certified I_3 (cyl_P4) -> Vir_N ((11)(11) - 2(22) = 4(T_Y T_Y) mod d re-verified with the VIR1 engine) -> (X, phi)
with rho -> (u, v) by the isometry: RATIONAL at every fibre and in the (u, v) commutant of {(1/2, -1/2), (1/2, 1/2), (-1, 0)}.
Momentum dictionary pinned by the exact match of the certified I_3 eigenvalue: our P2 = P_X^2, Q2 = pi^2 - rho^2 with pi the
phi evaluator momentum (the opers paper's pi = P_phi - rho).  Certified I_5 = kernel of ad I_3 at weight 6 (dim 1); its
eigenvalue equals the kit table vev_sol1_w6.json at t = 2, 7/2 (article P = i P_ours); **at t = 11/3 the kit normalisation
(dX)^6 = 1 is singular** (the kit's coefficients carry (3t - 11) denominators; the (dX)^6 and P^4 Q^2 coefficients vanish there),
so I_5 is normalised by its Q^6 eigenvalue coefficient = 1 ("ALT") at every fibre.  (u, v) commutant dimensions at w2 = 4, 6, 8,
10, 12: 1, 1, 1, 1, 1 at t = 7/2, 11/3 and 1, 1, 1, 1, **2** at t = 2 (an extra X-odd element at weight 6; the X-even one is
proportional to the certified I_5).

**Results (identical pattern at t = 2, 7/2, 11/3; Tower Z unless stated).**
| spin (w2) | (k, w) pieces of P_H I P_H and where they land | top piece (0, w2) in Solution 1's commutant? |
|---|---|---|
| 1 (4) | (0,4) = T_uv; (0,0) const | yes (automatic) |
| 2 (6) | (0,6) = c_2 I_2 (X-odd generator); (1,4) = c T_uv; (1,0), (3,0) const | **YES** |
| 3 (8) | (0,8) NO; (0,4) c T_uv; (0,0); (1,6) c I_2; (2,4) c T_uv; (2,0); (4,0) | **NO** |
| 4 (10) | (0,10) NO; (0,6) yes; (1,8) NO; (1,4), (1,0) yes; (2,6) yes; (3,4), (3,0), (5,0) yes | **NO** |
| 5 (12) | (0,12) NO; (0,8) NO; (0,4), (0,0) yes; (1,10) NO; (1,6) yes; (2,8) NO; (2,4), (2,0) yes; (3,6), (4,4), (4,0), (6,0) yes | **NO** |
Every piece of weight <= 6 (T_uv, I_2) and every constant is in Solution 1's commutant at its weight; every piece of weight >= 8
is not, at every q-order.  K (t = 2): spin 3 and 5 as Z's rows (top pieces NO, lower yes): P-II-K as predicted.
Spin-2 constants (canonical normalisation: Z's spin-2 density with coeff[(dX_0)^2 dX_1] = 1, Solution 1's X-odd generator with
coeff[(du)^2 dv] = 1, T_uv with coeff[(du)^2] = 1): c_2 = 1/5, 5/13, 2/5 = s/(4-s) = (t-1)/(t+3); (1,4) piece = -2 c_2 T_uv.
Commutation of the projection, order by order in q (S1uv, S2uv = pure exponentials (1/2, -+1/2); e^{-u} = Solution 1's third
screening; W0proj = (-(s/(4-s)) dv + q s(s-2)/(4-s)) e^{-u}, the projected W0): S1uv, S2uv: every order, every spin (as it
must be: they do not touch H~); e^{-u}: spin 2 every order; spin 3 fails at order 0 only; spin 4 at orders 0, 1; spin 5 at
orders 0, 1, 2; W0proj: fails at the lowest orders at every spin (the projection does not commute with P_H).
Information (R1): Z's unprojected charges do not commute with S3, S4 (K's do); the (u, v)-components of S3, S4 acting on the
top pieces: not twisted-exact.  Z's I_3, I_5 are not K's mod d.
(X, phi)-frame cross-check (sympy over Q(sqrt s, sqrt(4-s))): the spin-3 top piece is not in span(I_3^cert) + Im d (consistent
with the rational test).  Eigenvalue level (R3, no verdict weight): the top piece's eigenvalue is NOT proportional to the
certified E_3; with products allowed it is c E_3 + a E_1^2 + b E_1 + d exactly (E_1 = P_X^2 + pi^2), a != 0, at every fibre
(Z and K).
**Post-hoc products reading (reported, not a verdict; tested after the literal reading failed):** a "polynomial in lower charges"
at the density level = normal-ordered polynomial in the lower densities, which is local.  Spin 3: T^(8) = c_3 I_3^cert +
a :T_uv T_uv: mod d EXACTLY (Z, canonical normalisation (dX_0)^4 = 1 and the kit normalisation of cyl_P4): (c_3, a) = (2/45, 8/5),
(-50/1053, 200/39), (-16/245, 192/35) at t = 2, 7/2, 11/3; K at t = 2: (22/135, 152/135).  This is AUTOMATIC: the top piece lies in
the commutant of the pair {S1uv, S2uv} (dimensions 1, 1, 2, 2, 6 at w2 = 4, 6, 8, 10, 12), which at weight 4 is spanned by I_3 and
:T_uv T_uv:.  Spin 4: T^(10) not in span{I_4, :T_uv I_2:}; spin 5: T^(12) not in span{I_5, :T_uv I_3:, :T_uv :T_uv T_uv::,
::T_uv T_uv: T_uv:} (nor the sub-bases).  The products reading rescues nothing beyond spin 3.
**Controls (PREREG sec. 5), every fibre:** wrong dictionary (tower at s' = s + 1/7 versus Solution 1 at s, compared in the
(X, phi) frame at spin 3 and at the eigenvalue level at spin 5): fails.  Perturbed screening vector (W0 momentum (-1, 0) ->
(-8/7, 0)): the tower commutant is empty at w2 = 8 and 12 (nothing to project): fails.  Perturbed Solution 1 vector (third
vector (-8/7, 0)): the (u, v) commutant is empty at w2 = 8, 12: the top pieces fail (the weight <= 4 pieces still land on T).

## 8. What the coset image is (coset_image.py; logs coset_image_*_t2.log)

The q = 0 top pieces T^(w) of Z's projections at w2 = 4, 6, 8, 10, 12: all lie in the commutant of the N = 2-Liouville pair
{S1uv, S2uv} (the coset W-algebra of (u, v)); T^(4) = T_uv and T^(6) = c_2 I_2 commute with everything of Solution 1's; T^(8),
T^(10), T^(12) commute with S1uv, S2uv (and with the trivial (du)e^{-u} = -d e^{-u}) but with none of e^{-u}, e^{(2/s)u}, e^{-2u},
(dv)e^{-u}, e^{(1/2, 0)}, e^{(0, +-1/2)}; and they do NOT commute with each other ([T^(6), T^(8)], [T^(8), T^(10)], [T^(8), T^(12)],
[T^(10), T^(12)] all non-zero mod d).  So the H~-vacuum projection of the N = 2 tower is not a commuting family on (u, v) beyond
spin 2: the coset image, as the PREREG defines it, is Solution 1's I_1, I_2 plus non-integrable (u, v) densities of the
Liouville-pair commutant (I_3 + :T^2:-type at spin 3).  Same for K at spins 3, 5.

## 9. Hold-out (t = 13/5, p^2 = 8/9): registered predictions, written before any hold-out run

- P-II-Z as ruled (density level): the q = 0 top pieces of the spin-3 and spin-5 projections are NOT in Solution 1's commutant
  at t = 13/5 (and the pieces of weight <= 6 and the constants are); P-II-K the same.
- Spin 2 (functional-form test; a (1,1)-rational function is fixed by three fibres, so the fourth is a genuine test): c_2(13/5) =
  (t-1)/(t+3) = 2/7; the (1,4) piece = -2 c_2 T_uv = -4/7 T_uv (normalisations of sec. 7).
- Products reading, spin 3 (same caveat; low prior, the interpolants are not natural): (1,1) interpolants c_3(t) =
  2(21t - 61)/(3(73t - 431)) -> 32/1809 and a_3(t) = 24(101t - 141)/(53t + 809) -> 2432/789 at t = 13/5.
- The kit normalisation of I_5 is regular at t = 13/5 ((3t - 11) != 0), so the eigenvalue comparison at spin 5 can use it there.

## 10. Disclosures, deviations, errors

- R1-R4 as sent and ruled; R3's "no products" reading is the verdict basis; the products reading is post hoc (sec. 7).
- V2: the first flow fits failed; the cause was the omitted sign (-1)^|theta| of the graded bracket of an ODD charge with a
  component behind theta's, found by observing that the bottom and middle components fitted Tian-Liu (2) while the top did not;
  corrected in classical.sf_bracket; all fits then unique.  The logs of the failed attempts were overwritten by the reruns
  (same file names); the diagnostic logs diag_flow_Z.log, diag_basis_*.log, diag2.log, diag3_susy_invariance.log are kept; the
  sequence is in NOTES.md.  No parameter was refitted: alpha, tau, a are the equation's own unknowns.
- Level II: two corrections of the comparison script during the run, before any verdict was formed: (a) the q^k coefficients
  were first tested against the weight-w2 commutant (KeyError); the pieces must be tested at their own weight (the H~
  contractions lower it); (b) the proportionality constant was first extracted from reduced vectors (wrong scaling, assertion)
  and replaced by an exact linear solve.  The level2_*.log files are the final version's complete stdout.
- The script's line "real projection ... (must be True)" was written under the PREREG's hypothesis; its result False is the
  finding, consistent with the rational test.
- t = 11/3: the kit normalisation of I_5 is singular (sec. 7); the Q^6-normalised I_5 is used; the PREREG's fibre choice hits
  this zero of the (dX)^6 coefficient.  At t = 2 the (u, v) commutant at w2 = 12 is two-dimensional (extra X-odd element).
- The other seat's file code/super2_n2_tower.py (untracked) was read for its frame and snapshotted (inputs/, sha256
  1d2c4b66...); nothing of it was executed or used in the computations; all engines used are the registered ones plus this
  seat's new files.
- Sources read only: Zeitlin hep-th/0606129 (text + (22) by hand), KL 2312.09105 (text + PDF pages 4-5 rendered for (14)-(22)),
  FL, KL19, FNP (texts), Mao-Liu 1705.03997 and Tian-Liu 1108.2110 (texts + PDF pages rendered for (1)-(2)); nothing executed.
- Symbolic-in-p runs: V1 and the weights symbolic in s; Level I symbolic in s; towers and Level II at the fibres (rational).
- Not done: the hold-out (awaits the lead's commit); Z's spin 4/5 and K's spin 5 "products" bases beyond the ones listed.

## 11. Files

NOTES.md, PRIMARY.md, SHA256SUMS; inputs/ (B2_PREREG.md 9ece3fda..., super2_n2_tower snapshot, INPUTS.SHA256); gm.py (metric-aware
layer over the registered engines; classical one-contraction brackets), n2frame.py (frame, V1, weights, vertices), towers.py,
level1.py, v3.py, classical.py (V2), diag_flow.py, diag_basis.py, diag2.py, diag3.py (V2 diagnostics), coset.py (bosonisation,
rotation, projection), sol1.py (certified Solution 1 in the coset frame), level2.py, extras.py, products.py, validate2.py,
coset_image.py; json outputs n2frame_*.json, towers_*.json, v3_s*.json, classical_*.json, sol1_t*.json, level2_*.json,
extras_t*.json, products_*.json; logs/ (every run's full stdout: n2frame_*, towers_*, level1, v3_*, classical_*, diag*,
coset_validate_*, sol1_*, level2_*, validate2_*, coset_image_*, extras_*, products_*, holdout_predictions.log).
