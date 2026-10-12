# VIR2 (Fable seat, 2026-10-01) -- exponential screenings of the cylindrical charges in two bosons

Seal: `SEAL_VIR2.md`, sha256 `9ac0c5560861245c9e33ee228a96c717fa50555b2bf63d9fdde66633b483a499`,
2026-10-01T17:14:07Z, sent to the lead before any machine computation.  The seal was NOT blind
for spin 3 (pencil work disclosed in it); I_5, I_7, Solution 2 and completeness were blind.

## Verdict

**On every solution and at every t the cylindrical charges commute with a pair of mixed
vertex operators V_+- = exp(+-a X + beta phi), phi the Coulomb-gas boson of the Y-sector's
Vir_N.  The pairs are conserved against I_3 and I_5 identically in t (symbolic), against the
I_7 hold-out at 6 fibres and I_9 at 3 fibres.  The fail-fast did not fire.**

**Conversely the screenings determine the charges: for each solution the space of ALL
two-boson densities (no parity, no Virasoro assumption) commuting with {pair(s), the Virasoro
screening} is one-dimensional modulo total derivatives at weight 4 and at weight 6, and it is
the cylindrical I_3, I_5.**

With p = beta_+ (the Virasoro screening exponent; p^2 = 2u^2 = 2(t-1)/(t+1), rho = p/2 - 1/p,
N = 1 - 12 rho^2, the other Virasoro exponent is beta_- = -2/p):

| | pair | beta | a^2 | a^2 in t | h_phi | norm a^2 + beta^2 | v_+ . v_- | sealed? |
|---|---|---|---|---|---|---|---|---|
| Sol 1 | S1 | -p/2 = -beta_+/2 | 1 - p^2/4 | (t+3)/(2(t+1)) | (3p^2-4)/8 | 1 | (p^2-2)/2 | yes |
| Sol 2 | S2 | -p = -beta_+ | 1 - p^2 | (3-t)/(t+1) | p^2 - 1 | 1 | 2p^2 - 1 | form only (level-3 degenerate, a^2 -> -1 at N = 1) |
| Sol 3 | S3 | 1/p = -beta_-/2 | 1/p^2 | (t+1)/(2(t-1)) | (3/p^2 - 1)/2 | 2/p^2 | 0 | yes |
| Sol 3 | S3' | -p = -beta_+ | p^2 | 2(t-1)/(t+1) | p^2 - 1 | 2p^2 | 0 | NO -- predicted dead at I_5, wrong |

Sealed predictions: S1 HOLDS, S3 HOLDS, S2 HOLDS (50% call), S-Gram HOLDS as printed,
S-exponents HOLDS.  **S-complete is REFUTED on Solution 3**: a second pair S3' exists and
survives I_5, I_7, I_9.  It is the dual of S3: S3'_+- = -p^2 S3_-+ (same two lines in the
(X, phi) plane, exponents with product -2).  Apart from that the solution set is complete:
at 6 fibres (4 real p, 2 imaginary p) Groebner + solve gives exactly these pairs, for I_3 alone
and for {I_3, I_5}; there is no I_3-only solution.

NOT shown: any oper, any ODE.  No algebra is named beyond the printed statement of sec. 4.

## 1. What is new (lead's point 1)

BOP's class e^{aX + i b u.Y} uses the literal Y bosons; for N > 1 it generates the full
Verma-type Vir_N module and (pencil result (i) of the seal, confirmed: a generic primary
would need cross coefficient 2) cannot hold a mixed screening.  The pairs above use DEGENERATE
Vir_N primaries -- level 2 for S1 and S3, level 3 for S2 and S3' -- which exist only in the
Coulomb-gas realisation.  At N = 1 (rho = 0, phi = Y) they reduce to BOP's N = 1 solutions
(control C2 below).  ALGEBRA1's screenings act in the pillow's own three bosons; pure-X
screenings (scan Part B) are prior art and appear here only as a control.

## 2. Method

- Engine: `vir2_lib.py` on top of VIR1's `vir_lib.py`.  Density -> (X, phi): the HV
  decomposition of VIR1 (total derivatives absorbed first), each Vir state rebuilt with
  T_phi = (1/2)(dphi)^2 + rho d2phi.  Gates (`w0_gates.py`, exit 0): Virasoro algebra of T_phi
  with c = 1 - 12 rho^2 on 294 triples; h(e^{beta phi}) = beta^2/2 - rho beta;
  (11)(11) - 2(22) -> exactly 4 (T_phi T_phi); I_3 images even in X.
- Screening test: R = Res P(z)V(w)/V(w) reduced to the unique normal form with no bare dphi
  letter modulo the twisted derivative d + a dX + beta dphi; all normal-form coefficients must
  vanish (5 conditions at spin 3, 16 at spin 5, per pair, generic fibre).
- Machine parametrisation: p instead of the sealed u (p = sqrt2 u), because rho, N, s, the
  densities and all pair data are then rational in p.  Imaginary p covers N > 25 (|t| < 1).
- "Identically in t": `w3_symbolic.py` -- symbolic in p, I_3 and I_5, all four pairs and both
  Virasoro screenings, with four negatives per solution that fail.  Fibre run: `w4_main.py`
  -- 24 real-p + 6 imaginary-p fibres.

## 3. Results in detail

Fibre run (`w4_main.log`): S1, S2, S3, S3' conserved against I_3 and I_5 at 30 of 30 fibres
each; hold-out I_7 6 of 6 each; I_9 3 of 3 each (promoted charges of the repo, read-only).
The two-boson images of I_3 and I_5 commute (gate, all three solutions).

Do the screenings determine the charges?  (`w6_structure.log`, two fibres, post-seal.)
Commutant dimension modulo total derivatives, all two-boson densities:

| screening set | weight 4 | weight 6 |
|---|---|---|
| Virasoro beta_+ alone (= Heis (x) Vir) | 5 | 13 |
| Sol 1 pair alone / Sol 2 pair alone | 2 | 5 |
| Sol 3: S3 alone, or S3 + S3' | 3 | 7 |
| pair(s) + Virasoro beta_+, each solution | **1** | **1** |

and the one-dimensional space is the cylindrical I_3 (I_5) in every case.

Weight 8 (`w7_weight8.log`, p = 3/2, 185 monomials): the commutant of {pair(s), Virasoro} is
again one-dimensional mod d on each solution and equals the repo's promoted I_7.  So I_7 is
REPRODUCED from the screenings alone, with no input from the Ward engine or the ad I_3 kernel.

Spin spectrum (`w8_spectrum.log`, two fibres, weights 2..7), commutant dimensions:
    Sols 2 and 3:  1, 0, 1, 0, 1, 0   (spins 1, 3, 5 only; all X-even)
    Sol 1:         1, 1, 1, 1, 1, 1   -- one charge at EVERY weight; those at odd weight
                                         (spins 2, 4, 6) are X-odd.
My first expectation (1,0,1,0,1,0 everywhere) was wrong for Solution 1 (v0 log kept).  It is
consistent with the record: on B = -1, i.e. all of Solution 1, items 122-131 found a weight-3
current whose zero mode is odd under the X reflection.  The cylindrical family is the X-even
half of Solution 1's commuting family.

Solution 3 in rotated bosons (symbolic in p): with U = (X + phi)/sqrt2, W = (phi - X)/sqrt2,
    I_3 = 2(3p^2 - 2) [(T_U T_U) + (T_W T_W)] - 12 (p^2 - 2) T_U T_W   mod d,
T_U, T_W two commuting Virasoro algebras with rho_U = rho_W = (1/p - p)/sqrt2,
c_U = c_W = -(2p^2 - 3)(3p^2 - 2)/p^2.  (c_U = 1 at the anchor, -2 at the paperclip where the
cross term vanishes, -2 again at t = 5/3.)  So Solution 3 lies in the intersection of
Heis_X (x) Vir_N and Vir_U (x) Vir_W.

Exponents of the record.  For all four pairs, identically in p:
    B : C = a^2 : beta^2        ((A, B, C) = (eps, -a_X, -a_Y), A + B + C = -2),
i.e. the record's X- and Y-exponents are proportional to the squared X- and phi-components
of the screening vector.  Normalisation (B, C) = -k (a^2, beta^2): k = 1/a^2 (S1),
2/(3 - p^2) (S2), 2p^2(p^2-1)/(3p^2-4) (S3), 2(p^2-1)/(p^2(3p^2-4)) (S3').  For Solution 1 the
sealed identity holds in full: (A, B, C) = (v_+ . v_-, -a^2, -beta^2)/a^2.  I have no uniform
formula for k; pattern not derived.

## 4. Gram matrices (lead's point 5; sealed ones reproduced exactly)

v = (X-component, phi-component), G = v_i . v_j in the engine's convention.

- Sol 1 {v_phi+ = (0,p), S1+, S1-}: v_phi+ + S1+ + S1- = 0;
  G = [[p^2, -p^2/2, -p^2/2], [-p^2/2, 1, (p^2-2)/2], [-p^2/2, (p^2-2)/2, 1]]  (= sealed, in u).
- Sol 2 {v_phi+, S2+, S2-}: 2 v_phi+ + S2+ + S2- = 0;
  G = [[p^2, -p^2, -p^2], [-p^2, 1, 2p^2-1], [-p^2, 2p^2-1, 1]].
- Sol 3 {v_phi- = (0,-2/p), S3+, S3-}: v_phi- + S3+ + S3- = 0;
  G = (1/p^2) [[4, -2, -2], [-2, 2, 0], [-2, 0, 2]]  (= sealed (1/u^2)[[2,-1,-1],[-1,1,0],[-1,0,1]]).
  Dual set {v_phi+, S3'+, S3'-}: 2 v_phi+ + S3'+ + S3'- = 0; G = p^2 [[1,-1,-1],[-1,2,0],[-1,0,2]].
  Integer matrices A_ij = 2 G_ij / G_ii, order (S+, v_phi, S-):
  [[2,-2,0],[-1,2,-1],[0,-2,2]] for the first set and its transpose [[2,-1,0],[-2,2,-2],[0,-1,2]]
  for the dual set; both have determinant 0.  I have NOT located a printed table for these in
  papers/, so I name nothing.
- All full Gram matrices have rank 2 (`w5_points.log`).

One statement located in print (Litvinov-Spodyneiko arXiv:1609.06271, page 5, eqs. (2.24),
(3.1), (3.2); text layer and rendered page both read).  LS normalise phi phi ~ -log, so
Gamma = -G; their fermionic root has Theta = -1, i.e. |v|^2 = 1 here.
  Solution 1: (v_phi+, S1+-) is LS's block 3 [[2k, -k], [-k, -1]] with k = -p^2/2;
  (S1+, S1-) is block 5 [[-1, -k'], [-k', -1]] with k' = (p^2-2)/2; and LS's gluing rule (3.2)
  k_1 + k_2 = -1 holds at both fermionic roots.  So Solution 1's three screenings are built
  from LS's REGULAR blocks -- but closed into a 3-cycle (Gram rank 2), where LS treat chains.
  Solution 2: the boson-fermion entries ([[-p^2, p^2], [p^2, -1]], or [[-4/p^2, -2], [-2, -1]]
  with the dual root) are none of LS's six blocks.  Solution 3: nor are its entries.
  (LS's six blocks are derived from the existence of a spin-3 current.)
This matches item 131's finding that Solution 1 (B = -1) is the LS-regular locus.

## 5. The known points (lead's request b) -- values, with conservation recomputed there

(`w5_points.log`; I_5 at kernel-jump points is the limit class of each solution.)

| point | Sol 1: S1 (a^2, beta) | Sol 2: S2 | Sol 3: S3 ; S3' | what happens |
|---|---|---|---|---|
| t = 3, N = -2 (p = 1, rho = -1/2) | (3/4, -1/2) | (0, -1) | (1, 1) ; (1, -1) | Sol 3: beta = beta_+ (h_phi = 1) and beta = 2 rho (h_phi = 0); a = +-1 is in the scan's pure-X quartet {+-1, +-2}: each mixed operator is a product of single-sector c = -2 objects, the pairs are orthogonal -- two decoupled c = -2 sectors.  Sol 2: a -> 0, the pair collapses onto e^{-phi}. |
| t -> inf, N = 1 (p^2 = 2, rho = 0) | (1/2, -1/sqrt2) | (-1, -sqrt2) | (1/2, +1/sqrt2) ; (2, -sqrt2) | Sols 1, 3: norm-1 orthogonal pair and its phi-mirror = BOP's four vectors (+-1/sqrt2, +-1/sqrt2); Sol 2 = BOP's (-1, 2).  T_U T_W coefficient of Sol 3 vanishes: two c = -2 sectors in U, W. |
| t = 5 / -5, N = 0 | (2/3, -1/sqrt3) / (1/4, -sqrt3/2) | (-1/3, -2/sqrt3) / (-2, -sqrt3) | (3/4, sqrt3/2);(4/3, -2/sqrt3) / (1/3, 1/sqrt3);(3, -sqrt3) | For Sol 3 at t = 5 and Sol 2 at t = -5, I_3 = 4 (T_X + T_phi)^2 and the WHOLE curve a^2/2 + h_phi = 1 is conserved against I_3 and I_5 (checked at generic points on it; off it fails): every weight-1 primary of the total c = 1 Virasoro.  The pairs are points of that curve. |
| t = -1/3, N = 28 (p = 2i, rho = 3i/2) | (2, -i) | (5, -2i) | (-1/4, -i/2) ; (-4, -2i) | Sol 1: X-exponent sqrt2 = the prior-art pure-X screening's; (v_+.v_-, -a^2, -beta^2)/a^2 = (-3/2, -1, 1/2), the record's exact-point exponents. |
| t = +1/3, N = 28 (p = i, rho = 3i/2) | (5/4, -i/2) | (2, -i) | (-1, -i) ; (-1, -i) | Sol 2's data = Sol 1's at t = -1/3 (same a^2, beta, rho): the merger.  Sol 3: S3 = S3'. |
| t = 5/3, N = -25/2 (p^2 = 1/2) | (7/8, -sqrt2/4) | (1/2, -1/sqrt2) | (2, sqrt2) ; (1/2, -1/sqrt2) | S2 = S3': the Sol 2 = Sol 3 merger; S2 becomes orthogonal (v_+.v_- = 0). |
| t = -3, N = -2 other sheet (p = 2) | (0, -1) | (-3, -2) | (1/4, 1/2) ; (4, -2) | Sol 1: a -> 0, the pair collapses onto the Virasoro screening beta_- (h_phi = 1); eps has its pole there. |

All entries are conserved against I_3 and I_5 at the point.  I state what the formulas give; I
do not claim to have derived the known opers from them.

## 6. Controls

| | control | result |
|---|---|---|
| C0 | gates (sec. 2); two-boson images commute | PASS |
| C1 | a = 0: both Virasoro screenings pass I_3, I_5 (180 tests), I_7 (18), I_9 (9); beta = p/3 fails everywhere (90) | PASS |
| C2 | Sol 1 at (28,-9): pure-X (a^2, beta) = (2, 0) passes I_3, I_5; a^2 = 3 fails; fails at a generic fibre.  N = 1: (1/2, +-1/sqrt2) on Sols 1, 3 and (-1, 2) on Sol 2 pass | PASS |
| C3 | tamper A (rho -> -rho): each first pair fails, its mirror passes.  Tamper B (rho^2 from N+1, rho = sqrt(11) i/12 at the test fibre): each pair fails I_3 (4 conditions) and I_5 (15).  Tamper B' (rho -> rho + 1/5): the two-boson images of I_3, I_5 stop commuting | PASS |
| C4 | plant I_5 + (1/1000) x1^2 f2^2: flagged (11 conditions) | PASS |
| C5 | generic (a^2, beta) = (5/3, 2/7): 5 of 5 and 16 of 16 conditions nonzero; wrong a^2 fails I_7, I_9 everywhere | PASS |
| scoring tamper | `w4_main.py tamper` (pairs mirrored after the hash guard) | exit 2: S1, S3 REFUTED at I_3/I_5 and at the I_7 hold-out (S2 too); every control still passes |

## 7. Deviations, errors, caveats

- NOT BLIND at spin 3 (disclosed in the seal).  The blind tests are I_5, I_7, I_9, Solution 2,
  completeness, and "the screenings determine the charges".
- S-complete was wrong for Solution 3 (S3').
- Parametrisation changed from the sealed u to p = sqrt2 u for the machine runs; the sealed
  formulas and Gram matrices are reproduced exactly in u (`w5_points.log`).
- Tamper B: the commutation part was run with rho + 1/5 (a rational mismatch) instead of the
  sealed N + 1, because N + 1 gives an irrational rho; the pair-failure part uses exact N + 1.
- First gate run failed on wording: with HV generators ordered first the decomposition
  returned 4 L_-2^2 - 4 L_-4, the same class mod d; fixed by absorbing total derivatives
  first (`w0_gates_v0_FAILED_wording.log` kept).
- `nsimplify` was used in three helper lines of the first main run and in the first structure
  run, where it mangled an exact rational (crash, `v0_nsimplify/w6_structure_v0_FAILED.log`).
  House rule is never to use it.  Removed everywhere; main, tamper and structure re-run with
  exact conversions; the first main log (exit 0, same results) is kept in `v0_nsimplify/`.
- `w5_points.py` first version had a Gram-printing bug (a^2 not substituted in factored
  entries; one false "== sealed: False"); fixed, v0 log kept.
- BOP's remark that the N = 1 branch (2, 2) is "removed with [Q, I_5] = 0" is NOT reproduced:
  here (a^2, beta^2) = (2, 2) is conserved against I_3 and I_5 at N = 1 on Sols 1 and 3 (it is
  S3' at p^2 = 2, a Virasoro-dual screening of a c = -2 sector).  I did not re-run BOP's
  object; possibly a different statement.  Flagged, not resolved.
- I_7, I_9 rest on the repo's promoted charges.  I_9 is post-seal.  Fibre counts for the
  hold-outs are small (6 and 3) because each costs minutes.
- Completeness is established at 6 generic fibres, for a != 0 and beta != 0; at special
  fibres the set is larger (N = 0: a whole curve; N = 1; kernel-jump points).
- Formal N, as in VIR1.

## 8. Files and exit codes

`vir_lib.py` (VIR1 engine, copied), `vir2_lib.py`, `w0_gates.py` (0), `w1_lib.py`,
`w1a_spin3.py`, `w1b_spin5.py` (exploration at three u-fibres, 0), `w2_i5sym.py` (0),
`w3_symbolic.py` (0), `w4_main.py` (0, 13m50s; `tamper` exit 2), `w5_points.py` (0),
`w6_structure.py` (0), `w7_weight8.py` (0), `w8_spectrum.py` (0; v0 exit 1 on a wrong
expectation), logs alongside, `i5_sym.pkl`, `v0_nsimplify/`, `SEAL_VIR2.md`
(+ .sha256, .time), `SHA256SUMS`.
