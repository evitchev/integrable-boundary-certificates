# Verification kit: "Opers for the cylindrical integrals of motion"

This directory lets a reader re-run the computations behind the paper by E. Vitchev, "Opers for the cylindrical
integrals of motion" (2026), and check that the archived results reproduce. Every check is exact rational arithmetic
(sympy), except the numerical spectral check (mpmath) and the Lean proofs.

## Quick start

    python3 run_all.py --list          # the checks, with the paper section each supports
    python3 run_all.py --quick         # all quick checks (about 18 min on one core)
    python3 run_all.py --drill         # every check re-run with a deliberate tamper; each tampered run must FAIL
    python3 run_all.py --full          # quick checks plus the slower ones (about 1.5 h; num2_t1_zeros alone about 35 min)
    python3 run_all.py --only gate,thm2_ward --log out.txt

`run_all.py` resolves every path from its own location, so the working directory does not matter. It first checks every
file of the kit against `SHA256SUMS`, and the certified input tables against `INPUTS.SHA256`. It exits 2 if anything is
missing or altered. Each check then runs in a fresh temporary copy of `tree/`, so the archived outputs are never
overwritten. A check PASSES only if:

1. the script exits with the expected code; and
2. its output, with run times removed, equals line by line the archived output of the original run (or, for
   `num1_exact_M1`, every printed difference from the exact solution is below 1e-30).

Exit status: 0 if all requested checks pass, 1 otherwise.

Tested with Python 3.12.8, sympy 1.14.0 and mpmath 1.3.0. The two `second2_h2b_*` checks also pass, with byte-identical output, under Python 3.8.10 / sympy 1.13.3. Under sympy
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
| `sol3_vir5a_compare` (full) | Verification | All VIR5a cells against the tables, with its own controls. | registered predictions |
| `second2_h2b_t9o4` (full) | The opers (second-order descriptions) | The same exclusion at t = 9/4. | exact exclusion |
| `num1_shoot` (full) | Numerical spectral check | Zeros of D(E) equal the eigenvalues found by outward shooting (three levels, 1e-30 tolerance). | numerical |

## What the kit does NOT establish

- Agreement of local charges (WKB coefficients) is not equality of spectra. Equality of spectral determinants is checked
  numerically only for Sol 1 at one value of t (`num1_*`). Theorem 1 is a formal (Mellin-level) equivalence; its
  spectral version needs hypotheses the paper does not prove.
- The Sol 2 cells at n < 0 (t = 4 spin 9; t = 7 spins 3, 9) where the engine evaluates to zero are not understood.
- The high-precision Q-function integrations of NUM2 (`num2/num2_data.py`, `num2_point.py`) are not re-run; their outputs (`data_*.json`) are inputs to the fits, like the NUM1 determinant data. The extra non-WKB terms at t = 9/4 were found post hoc, and their coefficients have no closed form here.
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

`lean/OperMellin.lean` (Theorem 1's identities), `lean/WardLossOne.lean` (Theorem 2) and `lean/OperTop.lean` (the
rank-drop law of the top layer).

| file | sha256 |
|---|---|
| OperMellin.lean | 15d8ca57d0a2... |
| WardLossOne.lean | 058a244f7f6f... |
| OperTop.lean | 93a895a4c670... |

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
(3176 jobs)`; the log is `lean/BUILD_LOG_local_path_mathlib.txt`, and the files contain no `sorry`. The git-pinned
project in this directory was NOT itself built from a fresh download here. The Lean statements formalise the abstract
identities; the paper says which hypotheses enter.

## Changes to archived files

`PATCHES.json` lists every file of `tree/` that differs from the archived original:
- `lead_exc3_check/lead_compare_exc3.py` read the registered prediction from a path outside the archive. It now reads `PREDICTION_EXC3_snapshot.json` next to it, which has the same sha256 (7991f7af8dc89cc9...); the script asserts that hash.
- Two scripts opened the certified tables through an absolute path on the author's machine:
  `vir4d/d3_compare.py` and `vir5/e3_compare.py`. They now resolve the same files relative to their own location. One
  marked line was added to each; nothing else changed.
- Files that the kit does not execute (logs of superseded runs, helper scripts, process-id files, notes) had private
  absolute paths in their text (46 redaction entries in `PATCHES.json`). These were replaced by `<repo>/` or `<home>/`.

Because of these edits, the per-stage `SHA256SUMS` (and `SHA256SUMS_subset`) from the time of archiving no longer match for those files.
`PATCHES.json` gives the archived hash of each one.

## AI assistance

As stated in the paper, this work used AI tools (Claude, Anthropic; and other assistants) extensively for computation,
independent replication, literature checks and drafting. The stage directories are the archived work of these
AI-assisted sessions, each with its own report (`NOTE.md`), including errors and corrections. The author takes full
responsibility for all contents.
