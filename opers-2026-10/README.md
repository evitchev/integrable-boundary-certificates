# Verification kit: "Opers for the cylindrical integrals of motion"

This directory lets a reader re-run the computations behind the paper by E. Vitchev, "Opers for the cylindrical
integrals of motion" (2026), and check that the archived results reproduce. Every check is exact rational arithmetic
(sympy), except the numerical spectral check (mpmath), the two vir1r checks (numpy) and the Lean proofs.

## Quick start

    python3 run_all.py --list          # the checks, with the paper section each supports
    python3 run_all.py --quick         # all quick checks (about 2 h on one core)
    python3 run_all.py --drill         # every check re-run with a deliberate tamper; each tampered run must FAIL
    python3 run_all.py --full          # quick checks plus the slower ones (about 12 h on one core: the five spec1_s3_* full
                                       # checks take 40 min to 3 h each, num3_t1_zeros and num3_wkb_pred about 1 h each)
    python3 run_all.py --only gate,thm2_ward --log out.txt

`run_all.py` resolves every path from its own location, so the working directory does not matter. It first checks every
file of the kit against `SHA256SUMS`, and the certified input tables against `INPUTS.SHA256`. It exits 2 if anything is
missing or altered. Each check then runs in a fresh temporary copy of `tree/`, so the archived outputs are never
overwritten. A check PASSES only if:

1. the script exits with the expected code; and
2. its output, with run times removed, equals line by line the archived output of the original run (or, for
   `num1_exact_M1`, every printed difference from the exact solution is below 1e-30).

Exit status: 0 if all requested checks pass, 1 otherwise.

Release v1.2.6 was run in full (all 235 checks, and all drills) under Python 3.8.10 / sympy 1.13.3, and its 58 new checks and their drills also under Python 3.12.8 / sympy 1.14.0, from a fresh copy with HOME pointing to an empty directory (`logs/VERSIONS.txt`, `logs/v126_*`). Release v1.2.5 was run in the same way (`logs/v125_*`). Earlier: tested with Python 3.12.8, sympy 1.14.0 and mpmath 1.3.0. The v1.2.4 quick run (all 55 checks) used Python 3.8.10 / sympy 1.13.3, and the 27 checks new in v1.2.4 were also run under Python 3.12.8 / sympy 1.14.0 (`logs/VERSIONS.txt`). The two `second2_h2b_*` checks also pass, with byte-identical output, under Python 3.8.10 / sympy 1.13.3. Under sympy
1.14.0 the floating-point magnitudes they print as diagnostics (`min ... max ...` of |residual| over the numerical roots)
differ slightly. These numbers are therefore masked in the comparison for those two checks only. The exact counts of
non-zero residuals, computed in Q[z]/(m), are compared in full. Under a different sympy version, printed expressions can differ
in form (term order, factorisation) while being equal. A check that fails only in that way prints the first differing
line, which can then be compared by hand.

## Layout

- `tree/`: the archived stage directories and certified input tables. Each sits at the same relative path as in the
  author's private research repository (`results/lab/...`), so every script runs exactly as archived. Each stage
  directory keeps its own `NOTE.md` (the stage report: method, results, controls, disclosures), its seal and prediction
  files, its archived logs, and, where present, its own `SHA256SUMS` from the time of archiving.
- `INPUTS.SHA256`: the certified vacuum-eigenvalue tables the checks compare against. They are the cylinder vacuum
  eigenvalues of the charges, certified as described in the earlier papers of this series:
  - `vev/vev_sol{1,2,3}_w{6,8,10}.json` (spins 5, 7, 9);
  - `first_order/vev_w4_points.json` (spin 3);
  - `anchor11`/`anchor13`/`anchor15` profiles (spins 11, 13, 15).
- `predictions/INDEX.md`: the registered predictions, with their hashes and registration times.
- `lean/`: the Lean 4 proofs and a Lake project to build them (see below).
- `PATCHES.json`: every file of `tree/` that differs from the archived original, with both hashes (see "Changes to
  archived files").
- `drills.json`: the source edits used by `--drill`.
- `SHA256SUMS` covers every kit file except itself and `logs/` (which has its own); `INPUTS.SHA256` covers the certified tables.
- `logs/`: the outputs of the runs made when the kit was assembled.
- `MANIFEST.json`: the provenance of every input table (its path in the private repository, which `tree/` mirrors).

## What each check establishes

The checks follow the paper's section titles.

| check | paper section | establishes | kind |
|---|---|---|---|
| `gate`, `gate_missing_charge` | Verification (table) | Re-runs the five comparisons of registered predictions with the certified tables (Sol 1 spins 5-13 at 6 fibres; Sol 3 spins 11/13 and 5-13; Sol 3 spin 15 at 9 fibres; Sol 2 spin 11). It fails on any mismatch, on any negative control (the table at t + 1/100) that does not fire, and on any drift from the archived outputs. It also rejects a missing or extra cell, spin or control, overrides included: for every prediction the output's shape (every number replaced by '#') must equal the archived one (gate v2; v1 is kept as `strict_gate_v1.py`). `gate_missing_charge` re-runs the same gate; its drill nulls the VIR7 t = -2 spin-11 prediction, which must be rejected. | registered predictions |
| `sol3_operator_t2` | The opers (Sol 3); Verification | Recomputes, from the Gamma-symbol operator at k = 2 (t = 2), the spin-11 and spin-13 charges. The resulting prediction file must be byte-identical to the archived one. | recomputation |
| `sol3_vir4d_compare` | Verification | Compares the registered Sol 3 spin-11/13 predictions with the certified tables. | registered predictions |
| `sol3_ode_equiv`, `sol2_ode_equiv` | The opers (Theorem 1) | The three-term ODEs of Sols 3 and 2 give, at WKB orders 2-8, the same charges as their Gamma-symbol forms (k = 2 and k = 10/3). | formal WKB identity at two values |
| `sol2_holdout` | Verification | The registered Sol 2 spin-11 prediction against the hold-out table. | registered prediction |
| `sol1_suzuki` | The opers (Sol 1); Verification | Suzuki's equation with a centrifugal term and the stated dictionary reproduces 882 coefficients of the certified tables (spins 3-9 at 19 fibres) with no free parameter. The dictionary was found from these data. | post-hoc identification |
| `thm1_mellin` | The opers (Theorem 1) | The formal Mellin identities, symbolically in a, for deg L = 1, 2, 3, with negative controls, for sigma = n/a > 0 (the sympy proof uses positive symbols; the Lean proof assumes 0 < sigma). | symbolic proof (also in Lean) |
| `thm2_ward` | Verification (Theorem 2) | The closed-form first-order WKB formula G1 equals the Ward-derived loss-one ratios identically in (a, K, t) for all three solutions. G1's derivation from the operator and the Ward derivation are inputs, not proved here (machine-checked in Lean given those inputs). | symbolic proof (also in Lean) |
| `num1_exact_M1` | Numerical spectral check | Control: the numerical determinant equals the exact M = 1 determinant to better than 1e-30. | numerical control |
| `num1_doublet_average` | Numerical spectral check | From the archived high-precision determinant data (M = 5, three momentum points), the ±α average of every extracted grade equals the α-even WKB coefficient. | numerical re-analysis |
| `exc1_level1_I3`, `exc1_level1_I5` | Excited states (Sol 1) | The level-one trivial-monodromy condition reproduces the I_3 (post hoc) and I_5 (registered) level-one characteristic polynomials at t = 2, 9/4. | oper vs data |
| `exc2_closed_form`, `exc2_level1_I5` | Excited states (Sol 3) | The closed-form level-one oper reproduces the I_3 and I_5 blocks; the I_5 prediction was blind. | oper vs data |
| `br2_commutant` | Screening systems; Boundary interaction | At weights 4, 6 and 8 the commutant of the screening triple, modulo total derivatives, is one-dimensional (two fibres per solution); a tampered a gives 0. | exact linear algebra |
| `br2_dims` | Boundary interaction | The boundary dimensions of the screening pairs and their relevance ranges. | symbolic |
| `second1_leading` | The opers (second-order descriptions) | The leading coefficient of the cocycle numerator is ∓256 k^7 (k+2)(k+4) (the constant depends on normalisation; the paper states only the factor k^7 (k+2)(k+4)). | symbolic |
| `second2_gate1` | The opers (second-order descriptions) | Engine control: the Sol 1 dictionary is recovered from spin 3 and a tampered M fails. | control |
| `second2_h2b_t2`, `second2_h2b_plant` | The opers (second-order descriptions) | Sol 3 with a momentum-rotated Suzuki dictionary: at each of the 12 solutions of the fit equations (t = 2) the unfitted coefficients fail. A planted rotated dictionary is recovered exactly by the same pipeline, so the exclusion is not vacuous. | exact exclusion |
| `sol2_ext_exact_point` | The opers (Sol 2); Verification (table) | The Sol 2 operator at the exact point t = 1/3 (k = -1/2, n = 2) reproduces the certified spins 1-9 exactly; the built-in tamper (delta + 1/10) differs. Compared with the archived run (register item 322). | operator vs certified tables |
| `sol2_ext_all_fibres` (full) | The opers (Sol 2); Verification (table) | The same at all 14 fibres of the archived run (t = 1/3, 2, 3/2, 4, -2, 1/2, -1/3, 13/5, 11/2, 21/11, 7, 401/100, 399/100, 701/100). The check reproduces the archived output, which is not "all match": at t = 2 spin 7 (a multiple of n) the engine gives 0, and at the n < 0 cells t = 4 spin 9 and t = 7 spins 3, 9 the engine evaluates to zero while the certified charge is non-zero. Those cells are reported as engine-zero cells, NOT passes. | operator vs certified tables |
| `num2_fit_2_A_60`, `num2_fit_t10o3_A_60` (full: `num2_fit_2_B_60`, `num2_fit_2_C_60`, `num2_fit_t10o3_B_60`, `num2_fit_2_A_120`) | Numerical spectral check (the fourth-order paragraph) | Re-fits the archived high-precision Q-function data of the fourth-order (Solution 3) operator against the WKB coefficients. At t = 2 the sealed basis alone is incoherent; with the REGISTERED E^(-3/2) log E terms the grades agree (P1). At t = 9/4 the sealed fits fall short; the extra terms E^(-2.6), E^(-3.9) of the same family are POST HOC. The registered prediction P2 (non-WKB terms odd under l -> -l) is REFUTED at t = 9/4. | fit of archived data |
| `num2_exact_E0` (full: `num2_x0_control`) | Numerical spectral check | Controls: the solver against the exact E = 0 Q-functions (Meijer G), fourth and sixth order, 55-57 digits; in full mode, the generic solver against the NUM1 determinant. | numerical controls |
| `num2_wrong_ordering` (full) | Numerical spectral check | Negative control: the wrong operator ordering; the check reproduces its archived FAILING output. | control |
| `num2_t1_zeros` (full) | Numerical spectral check (the spectral Theorem 1 test) | The REGISTERED test T1: the Q-functions of the three-term fourth-order ODE and of the sixth-order Gamma-form ODE at t = 2 agree up to E-independent constants (to about 56 digits, two points, E = -20..500); zeros of the projected Q agree with the outward-Wronskian determinant. | numerical (registered) |
| `exc3_level1_I5` | Excited states (Solution 2) | The level-one oper of Solution 2 (two apparent singularities, exponents {-1,1,3}) reproduces the I_5 level-one block, a BLIND registered prediction, at t = 2, 9/4. The I_3 comparison was post hoc. | registered prediction vs data |
| `exc3_lead_blind` | Excited states (Solution 2) | The same prediction compared with independently computed data blocks (the coordinating session compared first). | registered prediction vs data |
| `proof2_p1_factor`, `proof2_p1_proof`, `proof2_p1_proof_b` | Spins that are multiples of n | On the Solution 3 rank-drop locus nu = (2K-1)/n = i, every coefficient of the record's (Ward-derived) loss-1, -2 and -3 layers vanishes: the poles of the Ward laws are enumerated with K or i symbolic, and the finitely many exceptional parameter values are checked cell by cell (i = 1..6). Off the locus the layers do not vanish (control). A theorem about the record's charges; on the operator side, losses >= 2 agree with the record as VERIFIED (spins <= 13), not proved. | symbolic proof |
| `proof2_general_case1` ... `case5`, `proof2_general_case1_near` | The opers (Theorem 1) | The formal WKB equivalence of the three-term ODE and its Gamma-symbol form for general P = (T^2 - l1^2) prod (T - r) and L = (T - l0) prod (T - q), (deg P, deg L) = (3, 1), (3, 2), (4, 2), non-integer k, at every non-degenerate WKB order. At k = 3/2 the Beta-pole cells differ at the pole and agree at k = 3/2 + 1e-6 (`case1_near`). The drill uses the wrong ordering x^c L(T) and must DISAGREE. | formal WKB identity at sample values |
| `proof2_aequiv_*` (Sol 3 at k = 3/2, 5/3, -6/5, 7/2; Sol 2 at k = 3/2, 5/2, -6/5) | The opers (Theorem 1) | The same equivalence for Solutions 3 and 2 at further k, WKB orders 2-8. | formal WKB identity at sample values |
| `thm3_constant` | Theorem 3 (Appendix A) | From the archived fourth- and sixth-order Q-functions at t = 2 (two momentum points, E = -20..500), (Q4_j/Q6_j) / m(theta_j) is the same constant for all four j and all E, to about 56 digits, and equals 1/omega computed from the exact E = 0 solutions. A gamma + 1/7 tamper spreads the values. | numerical check of the proved statement |
| `thm3_kernel` | Theorem 3 (Appendix A) | The Hankel realisation K x^mu = m(mu) x^mu (quadrature, complex mu), the integration-by-parts relation of Lemma 2 (symbolic on monomials; numerical on a non-monomial function), the intertwining of Proposition 1 on monomials (symbolic), the frozen zeros (a = 1..4) and the Laplace value of omega; each with a failing tamper. The docstring says L3 for the operator the paper calls L4. | symbolic + numerical checks |
| `joint1_C1` ... `joint1_C6` | Excited states | For each of the six level-one comparisons (three solutions, two values of t), the operator's (I_3, I_5) eigenvalue PAIRS equal the data's on the same states: the mixed trace tr(I_3 I_5) and the explicit pairs over the quadratic field. The drill crosses the operator's pairing and must fire. | oper vs data (exact) |
| `akns_spin3` | Relation to quantum AKNS | Fateev-Lukyanov's I_3 vacuum eigenvalue of quantum AKNS (hep-th/0510271, eq. (95), with a factor 1/16 restored on the (5n+4) term) equals the certified Solution 1 spin-3 eigenvalue in all six coefficients at all 29 fibres, at n = -(t-1)/(t+1) with one fixed dictionary. The literal printed (95) and a tamper n + 1/10 fail at every fibre. | published formula vs certified table |
| `akns_top` | Relation to quantum AKNS | FL's top layer against Solution 1's top layer, spins 3-9, at t = 7, 11/3, -2: under the assignment FL P^2 ~ Y, Q^2 ~ X the ratios of consecutive coefficients agree (the only scale factor is 1; nan marks a Gamma pole), and under the swapped assignment they do not. | published formula vs top layer |
| `numtab` | Numerical spectral check (evidence table) | Rebuilds the evidence table (167 rows: NUM1-3, fitted vs predicted coefficients) from the archived high-precision data, with its error estimate e (the larger of the window spread, the K2 vs K2 - 4 spread and the 60 vs 120 digit spread); 0 rows in tension. `make_table.py` formats the same table as CSV and LaTeX. | numerical re-analysis |
| `sol3_vir5a_compare` (full) | Verification | All VIR5a cells against the tables, with its own controls. | registered predictions |
| `second2_h2b_t9o4` (full) | The opers (second-order descriptions) | The same exclusion at t = 9/4. | exact exclusion |
| `num1_shoot` (full) | Numerical spectral check | Zeros of D(E) equal the eigenvalues found by outward shooting (three levels, 1e-30 tolerance). | numerical |

### New in v1.2.5

v1.2.5 adds 111 checks (each with a drill), so that every computation the paper reports from the archive can be re-run from the kit
(exceptions: "What the kit does NOT establish"). `python3 run_all.py --list` gives each check's paper section.

| checks | paper section | establishes | kind |
|---|---|---|---|
| `vir1_*`, `vir1r_*`, `vir2r_engine_control` | The Y-sector is Vir_N (Proposition 1) | I_3 identically in t; I_5 symbolically in t; I_7..I_11 symbolically in N on the certified promoted charges (`results/qt_promotion.json`); replications with separate code (a second formal-N engine at three values of t per solution; explicit components with a randomised rank test mod p at N = 25, 26, control 28) | exact; the randomised rank test is probabilistic |
| `vir2_*`, `vir2r_replication`, `vir2r_posthoc_vir_alone`, `br2_direct_paperclip` | Screening systems; Boundary interaction | the screening pairs are conserved by I_3, I_5 (symbolic in p) and I_7, I_9 (sample fibres); Gram matrices and exponents; commutant dimensions at weights 4, 6, 8 at 3-5 fibres per solution (not symbolically in t); LVZ's printed W_4 is conserved by the four hairpin exponentials | exact |
| `nz_*` | Spins that are multiples of n | the operator's coefficient vanishes at the drop points (ten points, i = 1, 2, 3; checked, not proved); the Beta prefactor; the limits at the cells (exact evaluation at two step sizes); the two-dimensional kernel on the jump tower for i = 1, 2, 3 and one-dimensional kernels at the controls; the eigenvalue law (i = 3 a sealed hold-out); the fusion-product law evaluated on the archived kernel records i = 3..6 (difference 0); W_U W_W at i = 1, 2 | exact |
| `red1_*`, `red2_*`, `pillow1_*` | The known points | anchor t = 3: the limiting fourth-order operator and its terminating energies (the two oscillator ladders, N = 0..3); N = 0: the log-radial Kondo limit; paperclip: Solution 1 and Solution 3 relations; exact point: the pillow's recurrence forms and invariants (no relation found) | exact |
| `nzero1_fibres` (full) | Verification (exceptional cells) | Solution 2 cells at n < 0: at t0 + 1/N, N = 25, 50, 100, 200, the operator coefficient equals the certified charge exactly; the differences from t0 shrink like 1/N | exact |
| `num3_*` | Numerical spectral check (third order) | re-fits of the archived Q data on the registered and the extended grids (the figures quoted), the wrong-ordering control, the controls X0/X1/X1', the post-hoc window and E^-3 diagnostics, the T1 ratio (re-verified from the archived Q values and recomputed at the worst cell; `num3_t1_zeros` re-runs the sealed script), one ODE point recomputed (`num3_spot_point`), the data assembly, and the WKB targets (`num3_wkb_pred`) | numerical |
| `exc4_*` | Excited states (Solution 3, level two; the half-level rule) | the level-two data blocks recomputed; the polynomial systems regenerated byte-identically from the archived scripts; the Singular/msolve outputs verified in pure Python (Groebner basis property, membership of the inputs, I_1 = vacuum + 4, the I_3 quintic, 40 = 5 x 8 ordered solutions; at order six: two solutions, half a level); the registered prediction against the blind data | exact |

New kit helpers (not archived files): `tree/results/lab/large_n/kit_num3/`, `tree/results/lab/large_n/kit_exc4/` (verifiers; their reference logs were
made when the kit was assembled) and `tree/lab/kit_gunzip_run.py` (the i = 6 kernel record ships gzipped, `INPUTS.SHA256`; the helper checks the hash
of the decompressed file before the archived script reads it).

### New in v1.2.6

v1.2.6 adds 58 checks (each with a drill): the spectral theorem for positive integer a (the generalisation of Theorem 3, Appendix A of the
revised text; 15 checks, `spec1_*`) and the Solution 1 entry "U(1) coset of N = 2" of the screening table (43 checks, `b2_*`, `b2b_*`).

| checks | paper section | establishes | kind |
|---|---|---|---|
| `spec1_scope`, `spec1_scope_sol1_a3` | Spectral theorem (scope) | Solution 3 at t = (3a+2)/(a+2) and Solution 1 at t = 2a+1 satisfy deg L = 1, sigma > 0, d >= 2 for a = 1..6 (Solution 2 does not: deg L = 2); the no-resonance hypothesis (H) holds at both momentum points at t = 2, 5/3, 11/5, 3, 5, 7 | exact |
| `spec1_kernel` | Spectral theorem; Appendix A (Lemma 2, Proposition A.1, Step 4) | the shift relations m(mu+sigma)/m(mu) = s - mu - r - c/2 and m(mu)/m(mu-n) = G_a(s-mu-r), symbolic in a, sigma, s, r, mu; the frozen exponents are zeros of m and the frozen lattice misses the other frozen exponents (a = 1..8); the intertwining on monomials at six points, symbolic in mu, l0, l1, E; K x^mu = m(mu) x^mu and Lemma 2(i) by 50-digit quadrature at three points | symbolic + numerical |
| `spec1_constant`, `spec1_saddle` | Spectral theorem; Appendix A (Step 5, Proposition A.2) | beta_O, beta_Gamma by exact conjugation (a = 1..6, both solutions, symbolic in l0, l1); the power and saddle identities (symbolic in a); omega(t = 2) = 3^(-1/6) pi^(-1/2); the E-correction exponents; the descent function (lower bound on a grid only) and the Laplace main term by quadrature; omega Q^O(0) = m Q^Gamma(0) from the Meijer-G values at a = 1..6 | symbolic + numerical |
| `spec1_e0_gauss` | Spectral theorem (the exact E = 0 identity) | the a frozen Gamma factors and one ODE factor form a complete Gauss product, and log(Q^O_j(0)/Q^Gamma_j(0)) - log m(theta_j) + log omega = 0 identically in s, r, theta, sigma, d, a | symbolic proof (given Gauss multiplication and the Meijer-G facts) |
| `spec1_ddt` | Solution 1 (relation to DDT at t = 3) | the a = 1 operator maps to the DDT sextic under x = mu y^(2/3), mu^6 = 3, with E_S = -kappa E, kappa = 4/(3 sqrt 3); the gauge x^(-l0/3) maps the Gamma form to DDT's third-order equation | symbolic |
| `spec1_s3_control`, `spec1_s3_archive`, `spec1_s3_spot`, `spec1_s3_*_a*` (5, full) | Spectral theorem (numerics at five points) | the solver reproduces the archived t = 2 Q data; the archived Q values at five points x 8 energies x 2 momenta satisfy the relation with the closed-form omega to the quoted accuracy, and the omega x 8/7 and multiplier-interchange tampers fail at every energy (`spec1_s3_archive`, no ODE solved); a re-solve at two (E, momentum) pairs per point agrees with the archive to >= 60 digits (`spec1_s3_spot`); the full runs, with the dps-140 repeat at E = 37 (`spec1_s3_sol3_a3`, `spec1_s3_sol1_a3`, `spec1_s3_sol3_a1`, `spec1_s3_sol1_a1`, `spec1_s3_sol1_a2`, full mode, 40 min to 3 h each) | numerical |
| `b2_n2frame_v1`, `b2_n2frame_vertices`, `b2_level1` | Screening systems (table entry for Solution 1, and the sentence after it) | the two-superfield frame closes into the N = 2 superconformal algebra with c = 3 - 12/p^2, symbolic in p^2 (36 OPE checks); T, G, G' commute with the five screenings, J only with the three J-neutral ones; after bosonisation and rotation, the Gram matrix of those three equals that of Solution 1's (V+, V-, e^{p phi}), identically in p^2 | exact |
| `b2_classical_Z`, `b2_classical_K` | Screening systems | the classical limits of the two towers fit the N = 2 KdV flow with a = 4 (Zeitlin's sl(2\|1)^(1) screenings) and a = -2 (Kolyaskin-Litvinov's sl(2\|2)^(1)), each fit unique | exact |
| `b2_towers_*`, `b2_sol1_*`, `b2_level2_*`, `b2_validate2_*` | Screening systems | the towers' charges (commutants mod d, w2 = 2..12) and Solution 1's certified I_3, I_5 in the coset frame, recomputed at t = 2, 7/2, 11/3, 13/5; the U(1)-vacuum projections of the spin-3, 4, 5 charges are not in Solution 1's commutant (spin 2 is, with c_2 = (t-1)/(t+3)), for a = 4 and for a = -2 at all four values of t (13/5 = the registered hold-out); the archived controls; an end-to-end eigenvalue validation of the projection | exact |
| `b2b_*` | Screening systems | Solution 1's I_2, I_3, I_5 do not commute with the a = 4 and a = -2 charges (spins 2, 3, 5) at t = 2, 7/2, 13/5: no residue is a total derivative; the consistency rows (each family commutes with itself, also after bosonisation) hold | exact |

New kit helpers (not archived files):
- `tree/results/lab/large_n/kit_spec1/`: `verify_s3.py` (verifier of the archived S3 outputs) and `spot_s3.py` (re-solves the operators at
  two (E, momentum) pairs per point through the archived `spec1_s3.qpair`), with their reference logs made when the kit was assembled, and
  `ORIGIN_SHA256` (the stage's sealed hashes of the archived S3 outputs). `kit_spec1/ref/` holds the references of the archived scripts:
  the archived SPEC1 logs with two changes, the first line (`ARGV: run.sh ...`, written by the author's wrapper) removed and the author's
  absolute workspace prefix replaced by `<stage>/`. The checks mask that path, the interpreter-version line and the solver timings; the
  sha256 of each script printed in the banners is still compared. SPEC1's `deps_snapshot/ode.py`, `ops.py` and `inputs/proof3/t1.json`
  are byte-identical to `proof3/ode.py`, `ops.py` and `proof3/t1.json`; the SPEC1 scripts import their own copies, so the copies ship in `spec1/`.
- `tree/results/lab/susy4/kit_b2/fresh.py`: deletes the archived copy of a script's JSON output before the run and then runs the archived
  script unchanged. With `--same-json`, it also requires that the script writes the file again and that the new file equals the archived
  one as a JSON object (timing keys ignored where stated). Without this, a script that did not write its output could pass, because the
  archived copy would still be in place. The helper writes only to stderr, so it needs no reference log.

Notes on the B2 material:
- `results/lab/susy4/boundary/b2b/inputs/` was not committed with the stage. The kit rebuilds it from byte copies of the `b2/` files;
  their hashes are the entries `inputs/...` of `b2b/SHA256SUMS`.
- The stage notes (`b2/PRIMARY.md`, `b2/NOTES.md`, `b2/HOLDOUT_RESULT.md`, `b2b/NOTES.md`, `b2b/RESULT.md`) are included with private paths
  and session names redacted (`PATCHES.json`). `b2/NOTES.md` (line 69) states the Level II pattern "at all three fibres, Z and K", but the
  stage kept tower K runs only at t = 2 and 13/5. The tower K runs at t = 7/2 and 11/3 (`level2_K_t7-2.*`, `level2_K_t11-3.*`) were made
  on 2026-10-18 with the stage's own unchanged `level2.py` and the stage's inputs, by a replication and its re-run; they show the same
  pattern. `b2/PRIMARY.md` section 7 ("identical pattern at t = 2, 7/2, 11/3") is headed "Tower Z unless stated".
- At t = 11/3 the (dX)^6 coefficient of the kit's normalisation of Solution 1's I_5 vanishes (3t - 11 = 0); `level2.py` then prints the
  kit-normalised c_5 as "undefined at this fibre", and the spin-5 eigenvalue-level ratio is not defined there (c = None). The membership
  tests mod d do not use this normalisation and are unaffected.

## What the kit does NOT establish

- Theorem 3 (Appendix A of the paper) is a hand proof. `thm3_*` check its constant and its algebraic identities numerically and
  symbolically, and `lean/Theorem3*.lean` certify the algebra (see "Lean proofs"); the analytic steps (convergence, growth
  bounds, asymptotics, uniqueness, the quoted existence theorem) are not machine-checked.
- The spectral theorem for general positive integer a is a hand proof, like its t = 2 case. `spec1_*` check its algebraic identities by
  exact symbolic computation (not in Lean; the Lean files cover t = 2 only), its E = 0 case exactly (Gauss multiplication, given the quoted
  Meijer-G facts) and the full statement numerically at five points. The analytic steps (growth bounds, termwise action, the quoted
  existence theorem (T) at Poincare rank Nq in xi = x^(1/N), the uniqueness of Lemma 4 and the Laplace estimate) are not machine-checked.
  **Superseded by paper Appendix A (v1.2.7 note):** the chart is now ξ = x^(σ/M) with M = lcm(d − 1, D) and the contour is the vertical line
  v = 1 + iy; the previous description xi = x^(1/N) and the tilted steepest-descent contour are kept here for provenance.
  The descent constant of Proposition A.2 is shown only on a grid (`spec1_saddle`). The dps-140 repeat at E = 37 is re-run only by the
  full-mode `spec1_s3_*_a*` checks; the quick checks verify the dps-70 data.
- B2 runs at t = 2, 7/2, 11/3, 13/5 only, and B2b at t = 2, 7/2, 13/5 only. At these points the checks are exact, but nothing is proved
  symbolically in t. The dictionary c_{N=2} = 3 - 12/p^2, p^2 = 2(t-1)/(t+1), n_FL = -p^2/2 comes from reading the sources (stage notes);
  the scripts take it as given. The N = 2 KdV flow equation is entered from Mao-Liu (1) as transcribed in `classical.py`; the kit does not
  compare that transcription with Laberge-Mathieu. The post-hoc readings of B2 (products, coset image, extras) and the KL (22) validation
  V3 are not re-run.
- `akns_*` compare published formulas with the certified tables. They do not prove that Solution 1 and quantum AKNS coincide beyond
  spin 3 (lower layers at spins >= 5 are not printed by FL), nor the spectral identification, which Bazhanov-Tsuboi state as a conjecture.

- Agreement of local charges (WKB coefficients) is not equality of spectra. The identification of the WKB coefficients with the
  asymptotic coefficients of the spectral determinants is checked numerically for all three operators at the values of t of the
  paper's Section 6 (`num1_*`, `num2_*`, `num3_*`, `numtab`). Theorem 1 is a formal (Mellin-level) equivalence; its spectral
  version is proved for deg L = 1 and positive integer a (Theorem 3 at t = 2; its generalisation, new in v1.2.6, covers
  Solution 3 at t = (3a+2)/(a+2) and Solution 1 at t = 2a+1, see `spec1_*`), and checked numerically for Sol 3 and Sol 2 at t = 2
  (`num2_*`, `thm3_constant`, `num3_*`) and at five further points (`spec1_s3_*`). Solution 2 (deg L = 2) and non-integer a are not
  covered by the proof.
- The Sol 2 cells at n < 0 (t = 4 spin 9; t = 7 spins 3, 9), where the engine evaluates to zero, are explained in the paper as an
  evaluation artefact at coincident Gamma poles (checked numerically at the fibres t0 + 1/N, `nzero1_*`); a symbolic proof that
  they are removable limits is open.
- The high-precision Q-function integrations of NUM2 (`num2/num2_data.py`, `num2_point.py`) are not re-run; their outputs (`data_*.json`) are inputs to the fits, like the NUM1 determinant data. The extra non-WKB terms at t = 9/4 were found post hoc, and their coefficients have no closed form here.
- Not covered in v1.2.5 (stated in the paper as such, or argued in the stage notes only): the kernel dimension 2 at i = 4, 5, 6 of the
  jump tower and the kernel membership of the fusion-law witnesses at i = 3..6 (archived exact fast-path records; i = 3, 4 also certified by
  separate computations of the research repository not included here); completeness of the Singular/msolve solution sets of EXC4 (the kit
  verifies the bases and solutions, not that no solution is missing) and its ordering tamper; the ODE integrations of NUM3 except one point;
  the anchor statement that the limiting operator is not the least common left multiple of the two oscillator operators; decoupling of the
  other branches at N = 0; the exact-point gauge search beyond its finite grid; the W^ev(c = N, lambda) statement; the doubling argument and
  the Kac labels of the boundary fields.
- The certified tables are inputs. Their certification belongs to the earlier papers and certificates of this
  repository and is not repeated here.
- `num1_doublet_average` re-analyses archived determinant data. Recomputing those data takes many core-hours
  (`num1/num1.py` at dps 60 and 120; one 40-energy point took about 50 minutes at dps 60). The kit does not do this, except for the M = 1 control and, in full mode, the
  shooting check.
- The excited-state data blocks (`exc1/d1_*.json`, `exc2/d1_*.json`) were produced by the author's excited-state engine,
  which is not part of this kit. The checks compare the operators against those archived blocks.
- `second2_h2b_*` re-check the exclusion at the solutions computed by Singular/msolve (`comp_*.txt`, `sing_*.log`,
  `msolve_*.out` are included). The Gröbner computation itself is not re-run, since Singular is not required.
- A "registered" prediction is one whose hash was recorded before the comparison data were opened. The registration
  times are self-recorded (`predictions/INDEX.md`). The first public commitment is the v1.2.1 hash file.
- Several sealed predictions failed, and some results were found post hoc. Each stage's `NOTE.md` says which; the
  paper's verification table repeats this.

## Lean proofs

`lean/OperMellin.lean` (Theorem 1's identities), `lean/WardLossOne.lean` (Theorem 2), `lean/OperTop.lean` (the
rank-drop law of the top layer), and, new in v1.2.4, `lean/Theorem3Intertwining.lean` and `lean/Theorem3FrozenLattice.lean`
(Theorem 3, Appendix A): the multiplier recurrences m(mu+3) = -(mu+2) m(mu) and m(mu) = (mu-1)(mu-4) m(mu-6) for all complex mu,
the three identities of Lemma 2 and Proposition 1 (L4 K = x^6 K x^-6 L6, for any P and E) on coefficient functions; m(mu) = 0
exactly on 1 + 3Z>=0; and, under the no-resonance hypothesis (H), the non-vanishing of the indicial polynomial at the shifted
exponents. They do not cover the analytic steps of the proof. Their stage directory (`tree/results/lab/large_n/lean3/`) has the
author-side replay (four mutation copies that must fail to build, a checker negative control) and its log.

| file | sha256 |
|---|---|
| OperMellin.lean | 15d8ca57d0a2... |
| WardLossOne.lean | 058a244f7f6f... |
| OperTop.lean | 93a895a4c670... |
| Theorem3Intertwining.lean | 98d70c75384d... |
| Theorem3FrozenLattice.lean | 3ee2dd28637a... |

The full hashes are in `SHA256SUMS`. `lean/` is a minimal Lake project:
- `lean-toolchain`: `leanprover/lean4:v4.24.0`;
- `lakefile.toml`: requires Mathlib at tag `v4.24.0`;
- `lake-manifest.json`: pins Mathlib to commit `f897ebcf72cd16f89ab4577d0c826cd14afaafc7` and its dependencies to the
  revisions of that Mathlib release.

To build:

    cd lean
    lake exe cache get    # downloads Mathlib's prebuilt files (several GB)
    lake build

When the kit was assembled, the three files were built with the same Lean and the same Mathlib commit. That build used a
local, already compiled Mathlib (referenced by path instead of git). It completed with `Build completed successfully
(3176 jobs)`; the log is `lean/BUILD_LOG_local_path_mathlib.txt`, and the files contain no `sorry`. For v1.2.4 the five files
were built the same way: `Build completed successfully (3191 jobs)` (`lean/BUILD_LOG_v124_local_path_mathlib.txt`). The git-pinned
project in this directory was NOT itself built from a fresh download here. The Lean statements formalise the abstract
identities; the paper says which hypotheses enter.

## Superseded statements in archived notes

The stage notes are kept as written. Three statements in them were later superseded:
- `proof1/NOTE.md` says the formal-WKB equivalence (both forms give the same coefficient times the same Beta function) is
  "NOT PROVED in general". It was proved at the formal level for all P, L with deg P > deg L and real a in `proof2/NOTE.md`,
  part (2); the paper states this after Theorem 1 ("At the level of formal WKB series the local charges of the two forms coincide ...").
- `proof3/NOTE.md` says that for sigma < 0 the Hankel integral diverges. That is withdrawn: it is false as stated (a counterexample
  at t = 7 has convergent tails); sigma < 0 is simply not covered.
- `proof3/NOTE.md` is the first version of the proof of Theorem 3. The paper's Appendix A supersedes it: there the Laplace step
  uses a vertical contour v = 1 + iy (replacing the NOTE's tilted steepest-descent contour; the NOTE's claim about the expansion
  for |arg y| < pi is not used), the chart is ξ = x^(σ/M) with M = lcm(d − 1, D) (replacing xi = x^(1/N)), the existence of the
  normalised fourth-order solution is derived, the uniqueness lemma divides out the algebraic prefactor, and the quoted
  existence theorem is pinned to Balser (2000), Theorems 42(d) and 44. The old text is kept visible above for provenance.

## Changes to archived files

`PATCHES.json` lists every file of `tree/` that differs from the archived original:
- `lead_exc3_check/lead_compare_exc3.py` read the registered prediction from a path outside the archive. It now reads `PREDICTION_EXC3_snapshot.json` next to it, which has the same sha256 (7991f7af8dc89cc9...); the script asserts that hash.
- Two scripts opened the certified tables through an absolute path on the author's machine:
  `vir4d/d3_compare.py` and `vir5/e3_compare.py`. They now resolve the same files relative to their own location. One
  marked line was added to each; nothing else changed.
- Files that the kit does not execute (logs of superseded runs, helper scripts, process-id files, notes) had private
  absolute paths in their text (46 redaction entries in `PATCHES.json`). These were replaced by `<repo>/` or `<home>/`.
- v1.2.5: 21 archived scripts opened files through an absolute path or imported code from a private read-only clone of the repository's
  `code/` (among them `sol2f/common.py` and `sol2_ext_fable/common.py`, used by checks shipped since v1.2.2: those checks could not run
  outside the author's machine). They now resolve the same files relative to their own location; the code modules are included under
  `tree/code/` and `tree/lab/`, byte-identical to the repository's. Further files with private paths in their text were redacted as above
  (default paths that the checks override at run time, or files the kit does not execute; 41 further redaction entries). Every entry is in `PATCHES.json`.
- v1.2.6: three B2 scripts (`b2/gm.py`, `b2/sol1.py`, `b2b/b2b.py`) opened files through an absolute path or imported from a private copy of
  the stage directory; they now resolve the same files relative to their own location (the imported code is byte-identical). Four B2 stage
  notes had private paths and session names in their text, redacted as above. The SPEC1 files are not changed (their reference logs are
  described under "New in v1.2.6").

Because of these edits, the per-stage `SHA256SUMS` (and `SHA256SUMS_subset`) from the time of archiving no longer match for those files.
`PATCHES.json` gives the archived hash of each one.

New in v1.2.7: `logs/SHA256SUMS` is regenerated to cover every file in `logs/` except itself, and the kit
integrity gate (`run_all.py --quick` starts by checking `SHA256SUMS`, `INPUTS.SHA256`, and now also `logs/SHA256SUMS`) verifies it.
Erratum: v1.2.6 shipped v1.2.5's `logs/SHA256SUMS` byte for byte (sha256 1e9bd0b3...); it listed v1.2.5's `VERSIONS.txt` (1e88cc62...)
instead of v1.2.6's (64c2a36a...) and omitted all 54 `v126_*` log files. Nothing caught this because the opers root `SHA256SUMS`
excludes `logs/` by design and the integrity gate checked only the two root manifests. The B2 row labelled `(must be True)` in the
archived level2 runs (the real-dictionary projection at w2 = 8, printed for consistency with the rational membership test; False in
all eight runs, agreeing with that test's w2 = 8 non-membership) is covered by the new check `b2_level2_w8_agreement`; the archived
output keeps its label. `vir1r.py`'s `import numpy` was moved above the archived-log
open (recorded in `PATCHES.json`), so a missing numpy no longer empties `run_VIR1R*.log` before failing.
`kit_exc4/regen.py` now seeds Python's and sympy's random generators with 0 before it re-runs an archived EXC4 script. sympy's
multivariate factorisation draws random evaluation points, and some draws stall it: `exc4_regen_six` timed out at 3600 s in a v1.2.7
test run and in one v1.2.6 drill run, and one py3.8 run of a four-point check stalled for over 12 minutes (the `exc4_regen_t*` notes).
On `f_six_sol3.py`, 5 of 40 seeds ran past 120 s and the other 35 took about 5 s; every run that finished wrote a byte-identical file.

## AI assistance

As stated in the paper, this work used AI tools (Claude, Anthropic; and other assistants) extensively for computation,
independent replication, literature checks and drafting. The stage directories are the archived work of these
AI-assisted sessions, each with its own report (`NOTE.md`), including errors and corrections. The author takes full
responsibility for all contents.
