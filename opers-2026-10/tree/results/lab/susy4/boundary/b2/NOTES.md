# B2 running notes (Fable secondary seat) -- started 2026-10-13

Commission: lead (<session>), 2026-10-13.  PREREG = results/lab/susy4/boundary/B2_PREREG.md, commit c24a7af,
sha256 9ece3fda2966938a219b7e1bbd9af5b653a84130e8ab06da97f98aed8cba59fd (VERIFIED on the working tree and on the blob at c24a7af).
Erratum from the lead (register item 348): Zeitlin's quantum vertices are his eq. (22) of the local v1 text (classical ones (11)-(12);
(28) is [t,t] = 0).  Confirmed by reading sources/hep-th_0606129/paper.txt.  Sources are untrusted: read only, nothing executed.
Repository read-only.  Workspace <workspace>/b2/.  Python <python3.12> with code/ on sys.path (read-only).

## (1) Dictionary re-derivation (PREREG sec. 1, last bullet) -- 2026-10-13
Claim: c_{N=2} = 3 - 12/p^2, p^2 = 2(t-1)/(t+1), n_FL = -p^2/2.
- p^2 = 2(t-1)/(t+1): opers paper (papers/opers-2026-10/opers.tex, Setting): "We use p^2=2(t-1)/(t+1), rho=p/2-1/p, N=1-12rho^2".
- n_FL = -(t-1)/(t+1): opers paper sec. "The opers", Relation to quantum AKNS: n = -2/(M+1) = -(t-1)/(t+1) with M = (t+3)/(t-1)
  (Assis (2.2)-(2.4)); register item 337.  Hence n_FL = -p^2/2.  CONSISTENT.
- c_W = 2 + 6/n: FL hep-th/0510271 after (83): the parafermions "extend the notion of the Z_k parafermions ... to non-integer
  k = -n - 2"; Z_k parafermions have c = 2(k-1)/(k+2) [Fateev-Zamolodchikov], so c = 2(-n-3)/(-n) = 2 + 6/n.  KL19 1910.05947
  (3.16): c = 2 + 6/k for the cigar dual (their k = FL's n).  CONSISTENT.
- c_pf = c_{N=2} - 1: FNP hep-th/0406017 (2.10): T = T_{SL(2,R)/U(1)} - 1/2 (d phi)^2, J = i sqrt(c-hat) d phi, c-hat = k/(k-2)
  (2.12); so c_{N=2} = 3k/(k-2) = c_coset + 1.  Equating c_coset = 2 + 6/n gives n = k - 2 (the cigar level), CONSISTENT.
- Therefore c_{N=2} = 3 + 6/n = 3 - 12/p^2.  MATCH with the PREREG dictionary.  No STOP.
- Parameter maps into the sources' conventions (derived from exchange exponents, convention-independent):
  KL 2312.09105 (15) c = 3(1 + 4/b^2)  =>  b^2 = -p^2  (b imaginary at every primary fibre and at the hold-out).
  Zeitlin (23) beta: his fermionic vertices e^{-Phi_1}, e^{-Phi_2} have mutual exponent -beta^2/2 (from (31): W_i W_j = -q^{-1} W_j W_i,
  q = e^{i pi beta^2/2}, with the fermion simple pole), KL's S1 S2 have -b^2/2 ((13),(16)): beta_Z^2 = b^2 = -p^2.
- Primary fibres: t = 2, 7/2, 11/3 -> p^2 = 2/3, 10/9, 8/7 -> c_{N=2} = -15, -39/5, -15/2; hold-out t = 13/5 -> p^2 = 8/9, c = -21/2.

## Frame (exact rational arithmetic at every fibre)
Two N = 1 superfields (X_j, psi_j), j = 0, 1, diagonal metric eta = (p^2, -p^2): <dX_j dX_j> = eta_j/(z-w)^2, <psi_j psi_j> = eta_j/(z-w)
(engine "+log" convention).  Derived from the two null screening momenta m_1, m_2 (KL S1, S2; Zeitlin e^{-Phi_1}, e^{-Phi_2}):
m_i.m_i = 0, m_1.m_2 = -b^2/2 = p^2/2; basis e_+ = m_1 + m_2 (norm p^2), e_- = m_1 - m_2 (norm -p^2); m_1 = (1/2, 1/2),
m_2 = (1/2, -1/2), m_3 = -m_1, m_4 = -m_2 (KL S3, S4); Zeitlin's third vertex: momentum -(m_1 + m_2) = (-1, 0).
Background charge from h(SUSY null vertex) = 1: T = (1/p^2)[1/2 (dX_0)^2 - 1/2 (dX_1)^2 - 1/2 psi_0 dpsi_0 + 1/2 psi_1 dpsi_1 - d^2 X_0],
c = 3 - 12 eta_0 lambda_0^2 = 3 - 12/p^2.  (To be verified by the engine in V1 with p^2 symbolic.)

## Reading points flagged to the lead before Level II (not adapted silently)
R1. KL's affine screenings S3, S4 carry J-charge +-2 (B1 log v2: J does not commute with S3, S4): they are not J-orthogonal, so
    "K's J-orthogonal screening vectors" are S1, S2 only; S3, S4 have no (u,v) image.  P-I-K is tested as: the J-neutral part of
    K's system = the N = 2 Liouville pair alone (two nodes), not Solution 1's three.
R2. Zeitlin's third vertex W_0 = int dtheta (D Phi_1 - D Phi_2) e^{Phi_1+Phi_2} has conformal weight p^2/2 = -b^2/2 (not 1) under the
    N = 2 stress tensor fixed by the null pair; the commutant mod total derivatives of oint W_0 is nevertheless well defined (the
    engine's twisted-exactness condition is weight-blind; the opers paper's V_+- likewise have weight p^2/4 under T_X + T_Y).
R3. Level II "plus a polynomial in Solution 1's lower charges": P_H I P_H is the zero mode of a local (u,v) density, so only
    LINEAR combinations of commutant densities of the same weight can appear at the density level; products of charges are
    non-local.  Tested as: D_proj(q) in the span of Solution 1's commutant at that weight with coefficients polynomial in q
    (q = 0: proportional to I_s with a non-zero constant).  The operator-level statement (products of eigenvalues) is reported
    separately at the eigenvalue level.
R4. Background charge: the coset (u,v) stress tensor inherited from the N = 2 T has background along the e_+ (= phi) direction
    with c_uv = 2 - 12/p^2 (= FL's 2 + 6/n), i.e. rho' = 1/p in the opers paper's phi; the opers paper's Vir_N uses rho = p/2 - 1/p
    (c_{X,phi} = 14 - 3p^2 - 12/p^2).  They differ by a total derivative (rho - rho') d^2 phi: immaterial for the commutant mod d
    and for Level II; reported at Level I as a difference.

## Timeline (UTC, 2026-10-13) and running record
- 13:51 workspace, dictionary re-derivation written; message to the lead (dictionary + R1-R4).  Lead: all four accepted.
- V1: first ansatz run found G' = i g' (imaginary in the indefinite frame); V1 rewritten in the complex basis G+- = (G -+ g')/sqrt2
  with J computed on the fly; ALL OK symbolic in s.  Weights and vertex commutation tables run.
- Towers w2 = 2..8 at s = 2/3, 10/9, 8/7, 3 (Z and K); V3 strict PASS; Level I run.  Message to the lead (cost at spin 3).  Lead:
  no objection; proceed; also report the spin-2 class.  Towers w2 = 10, 12 run (background).
- V2: classical commutants and the quantum -> classical limit (rational reconstruction) OK at once.  The flow fit against
  Mao-Liu (1) failed in three successive forms: (1) alpha pinned to 1 (not allowed: the equation is inhomogeneous) -> freed;
  (2) D_i used for the SUSY action on components -> Q_i added (both kept); (3) STILL no fit; diagnosis: the spin-2 flow fitted
  Tian-Liu (2) in the bottom and both middle components but not the top; all charges verified to commute with G, g', J, T
  (diag3) so the flow IS a superfield; hence the superfield construction was wrong: the graded bracket of an odd charge with a
  component behind theta's carries (-1)^|theta|, omitted.  Fixed in sf_bracket -> a = 4 (Z), a = -2 (K), unique; Tian-Liu (2)
  fits.  The failed runs' logs were overwritten by the reruns (same names); diag logs kept.
- coset.py written and validated (27 random word pairs, VA homomorphism); sol1.py: certified I_3 -> (u,v) rational; momentum
  dictionary found after adding the right candidate (Q2 = pi^2 - rho^2 with pi = evaluator momentum); kit I_5 normalisation
  singular at t = 11/3 ((3t-11)): ALT normalisation by the Q^6 coefficient.
- level2.py: KeyErrors from comparing pieces at the wrong weight (the H~ contractions lower the weight; and q carries weight):
  corrected to (q-power, weight) pieces; proportionality routine replaced by an exact solve.  Then: spin 2 PASS, spins 3-5 top
  pieces FAIL, lower pieces PASS, at all three fibres, Z and K; controls behave.
- validate2.py: end-to-end eigenvalue equality (registered fermionic evaluator vs projected density) for every charge at all
  three fibres: the projection machinery is validated independently of the comparison.
- coset_image.py, extras.py, products.py: characterisation; c_2 = (t-1)/(t+3); products reading rescues spin 3 only (automatic).
- holdout_predictions.log: (1,1) interpolants and the values at t = 13/5 written BEFORE any hold-out run.
- PRIMARY.md written; SHA256SUMS; report to the lead.  Hold-out not run (awaits the lead's commit).
