# Literature memo, 2026-09-11 (23:15): what the recent ODE/IM literature offers the cylindrical-oper search

Read from the papers' own text (arXiv abstract pages and, for two, the PDF via pdftotext); everything else from abstracts only.

## Directly usable: Masoero-Ruzza, arXiv:2605.24563 (May 2026, v3 Sep 8 2026), "ODE/IM correspondence at the free-fermion point"
Proves the BLZ completeness conjecture at c = -2 (alpha = 1): every solution of the quantum-KdV Bethe/QQ system is the spectral
determinant of a Crum-Darboux transformation of the harmonic oscillator -d^2/dx^2 + x^2 + (beta^2 - 1/4)/x^2 - 2 d^2/dx^2 log P(x^4)
(monster potentials = Laguerre Wronskians, parametrized by partitions); diagonalizes the first three KdV Hamiltonians at c = -2
with eigenvalues given by shifted symmetric functions.  RELEVANCE: on our curve t = 3 is (N, s) = (-2, -9), c_{1,2} = -2, the O(N)
free-fermion point, where Solution 3's hierarchy decouples into X (+) Y (sec. 24; cyl_hidden_kdv_scan).  CHECKED TONIGHT (lab
t3_kdv_anchor.py, staged for batch 29): Solution 3's certified e_3, e_5, e_7, e_9 at t = 3 have no mixed terms, the X- and Y-parts
are ONE function of Delta with Delta_X = -X/2 - 1/8 and Delta_Y = -Y/2, and at spins 3 and 5 both parts are EXACTLY BLZ's quantum-KdV
vacuum eigenvalues at c = -2, constants included.  So t = 3 is a SECOND EXACT ANCHOR on the curve: the cylindrical oper there is two
decoupled copies of the c = -2 harmonic-oscillator operator (NOT of pillow form: the pillow oper at t = 3 has c_X = 0 and loses both
momenta), and -- unlike (28, -+9), where only the vacuum is known -- the c = -2 anchor comes with ALL EXCITED STATES, proven
complete.  Uses: (a) test the excited-state data of sec. 38(ab) (level-1 and level-2 singlets) at t = 3 against the shifted-symmetric-
function eigenvalues; (b) the deformation problem now has TWO exact points with different operator forms (a scalar pillow oper at
(28, -+9), a direct sum of two scalar operators at (-2, -9)): any global candidate must interpolate between them, which excludes
'one scalar operator with the pillow symbol everywhere' outright and points at a 2 x 2 (block) structure that is diagonal at t = 3;
(c) Solution 2 at t = 3 has poles in its tables (the X-sector drops out there), so the anchor is Solution 3's.

## Confirmatory, no new tool: Gheisarieha-Yazdi-Arabi Ardehali, arXiv:2509.23869 (read in full, batch 26)
Large-N saddle point confirms LZ's (71) to NLO; states the generic-N ODE is not available; suggests the pillow's spherical limit at
N = 3.  Route 15 (v) filed.

## Method papers (W-algebra ODE/IM), abstracts only
- Kudrna-Prochazka arXiv:2508.20793: systematic WKB algorithm for local IM eigenvalues of Virasoro, W_3, W_4 opers in terms of Bethe
  roots; the 'mirror curve' (three-punctured sphere covered by the WKB curve).  Relevant if the cylindrical algebra is a W(2,4,6)- or
  W(2,4)-type algebra (sec. 38(ah), conjecture): their algorithm gives the eigenvalue/oper dictionary for higher-rank scalar opers,
  i.e. HIGHER-ORDER scalar ODEs -- a class the first-order programme has not touched (all eight identities are for second-order
  scalar operators).  A third-order oper (W_3-type) with the pillow symbol as its leading term is a concrete next class.
- Ito-Zhu arXiv:2408.12917 (WA_r, WD_r) and Ide-Ito-Kono arXiv:2604.07829 (WE_6): period integrals of higher-order ODEs along
  Pochhammer contours reproduce W-algebra IMs on the cylinder to sixth order -- the machinery for testing higher-order opers
  against certified charges, and the D_r case is the natural home of a W(2,4)-type (Solution 3's two-generator ring).
- Tanabe (sole author) arXiv:2604.14899, JHEP 08 (2026) (C(2)^(2), N = 1 SCFT; the memo first misattributed it to "Ito et al." -- corrected 2026-09-12 after the cognitive agent's INSPIRE check): twisted affine superalgebra linear problems with NS/R sectors -- a template
  for two-sector operators.
- Gaiotto-Lee-Wu arXiv:2003.06694, Gaiotto-Lee-Vicedo-Wu arXiv:2010.07325 (Kondo line defects, affine Gaudin, ODE/IM): boundary/
  defect integrable problems in chiral CFT with first-order (matrix) linear problems -- the closest existing framework for a
  NON-scalar oper of a boundary problem.
- Kotousov-Lukyanov-Shabetnik arXiv:2511.20530 (read in batch 13): lattice Z_r models near the free-fermion point; no boundary content.
- Masoero-Raimondo, Feigin-Frenkel-Hernandez opers arXiv:2312.01955, Commun. Math. Phys. 405, 193 (2024) (corrected from "CMP 2025"): the general theory of opers with trivial-monodromy singular terms.

## Not relevant despite the keywords: Komargodski-Popov-Rayhaun arXiv:2508.14963 (2+1d Boson-Kondo defects); Demjaha-Zarembo
arXiv:2506.17955 (AdS string on the Coulomb branch).

## Bottom line
One development changes the search: the proven completeness of ODE/IM at c = -2, together with tonight's exact check that
t = 3 on Solution 3 IS that point in both sectors, gives a second exact anchor with excited states.  The two anchors have
different operator forms, which settles that the global object is not one scalar pillow-symbol operator and makes the block
(2 x 2 or higher-order) structure the next class to build, with the c = -2 excited-state formulas as the first test.

## Citation check (2026-09-12, the cognitive agent, via INSPIRE): all fourteen identifiers resolve; two attribution errors in this memo
corrected above (2604.14899 is Tanabe alone; 2312.01955 is Masoero-Raimondo, CMP 405 (2024)); 2509.23869 is JHEP 04 (2026), 2408.12917
is Nucl. Phys. B 1010 (2025), 2604.07829 is JHEP 07 (2026), 2511.20530 is Nucl. Phys. B 1024 (2026).  New pointer from the rescreen:
Di Francesco-Mathieu-Senechal, "Integrability of the quantum KdV equation at c = -2", Mod. Phys. Lett. A 7 (1992), hep-th/9112063 --
the original free-fermion-point integrability paper behind the t = 3 anchor; to be read before it is cited.  Operational: arXiv's
export API blocks this host for hours after a burst; citation verification goes through INSPIRE.
