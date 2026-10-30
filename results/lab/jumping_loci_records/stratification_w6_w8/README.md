# Multiplicity stratification at weights 6 and 8 (spins 5 and 7)

## What the computation shows

For each sheet sigma = 1, 2, 3 of the cylinder I_3 density and each density weight W = 6, 8 (the charges of spin 5 and 7),
it finds, over all of Q-bar except the excluded poles t = +-1, every point t where the multiplicity m_sigma(t) is larger
than its generic value. Here m_sigma(t) is the dimension of the kernel of ad I_3 on the weight-W charge classes (densities
modulo total derivatives), at N = (t^2-25)/(t^2-1), s = -24t/(t^2-1), beta^2 = (t-1)/(t+1). Result: the generic
multiplicity is 1 on every sheet at both weights. Sheet 1 never jumps. Sheets 2 and 3 jump at t = 5/3 and t = 3 at W = 6,
and only at t = 3 at W = 8, with m = 2 at every jump. No deck image is a jump point: at t = -5/3 and t = -3
(beta^2 = 4 and 2) m = 1 on every sheet. No algebraic point of degree >= 2 jumps. Method, in exact arithmetic
(Python fractions and sympy) with the repository's engine (`code/mixed_engine.py`, `code/invariant_engine.py`) for the
I_3 densities, the Wick residues and the reduction modulo total derivatives:
- S1: the matrix M(t) of ad I_3 from weight-W classes to weight-(W+3) classes, exact over Q(t), saved with a scalar
  denominator cleared (M_poly = (t^2-1)^k c0 M); checked against the engine's numeric kernel and column pipeline at 42
  rational points.
- S2: the rank r over Q(t), bounded below by the rank at seeded random rational points and above by null vectors
  (Cramer vectors) verified identically in Q[t]; the gcd G of the r x r minors met over four row sets has every rank drop
  among its roots.
- S3: every irreducible factor of G other than t -+ 1 (all are linear) decided by the exact rank at its root.
- Checks: direct nullities at seven points; every factor of one minor D(t) of degree <= 8 that the gcd removed has full
  rank over Q[t]/(f); the engine's own kernel at t = +-5/3, +-11/5, +-2.

The scripts and logs call the density weight W "spin W".

Statements of the paper (papers/jumping-loci) that rest on it:
- Section "Record through spin 27", Table 1, footnote a: "That no other point jumps at weight 8 rests on an exact
  stratification over $\overline{\mathbb{Q}}$ \status{lab result}".
- Section "Methods", (ii): algebraic points the residue scans cannot see "are excluded only where the locus is classified
  exactly over $\overline{\mathbb{Q}}$: at weights 10--14 \status{certified}, and at weights 6 and 8 \status{lab result}".
- Section "The points $\beta^2=1/q$ and the open mechanism": "The deck images $\beta^2=q$, which have the same $N$, are
  not jump points wherever tested: [...] by an exact stratification at spins 5 and 7 \status{lab result}".

## Files

- `README.md`: this file.
- `run.sh`: re-runs every script below into a new output directory.
- `compare.sh`: compares that output directory with the records shipped here.
- `code/family.py`: S1 construction of M(t) over Q(t) from the engine; clearing, evaluation, save and load of matrices.
- `code/linalg.py`: exact rank, determinant, solve and interpolation over the rationals.
- `code/numfield.py`: exact rank over a number field Q[t]/(f).
- `code/sizes.py`: sizes of the monomial bases and class spaces at weights 1-11 (engine objects).
- `code/s1_build.py`: S1, builds and saves the six matrices.
- `code/s1_validate.py`: S1 validation at 42 rational points against the engine.
- `code/s2s3.py`: S2 and S3 for one matrix.
- `code/s3_extra.py`: direct nullities at t = -3, -5/3, 0, 9, -9, 5/3, 3 and the factorization of one minor D(t).
- `code/s3_crosscheck.py`: exact rank over Q[t]/(f) for the factors of D(t) of degree <= 8.
- `code/s3_points.py`: nullity of M and the engine's kernel at t = +-5/3, +-11/5, +-2.
- `code/s4_tamper.py`: negative control, builds a tampered matrix (sheet 3, W = 8, one I_3 coefficient times 11/10).
- `data/M_sheet{1,2,3}_W{6,8}.json`: the six matrices M_poly (integer coefficient lists, low to high degree, with row and
  column labels and the clearing scalar).
- `data/M_sheet3_W8_tamper.json`: the tampered matrix of the negative control.
- `results/S23_sheet{1,2,3}_W{6,8}.json`: S2 and S3 per matrix (random points and ranks, row sets, det A, null vectors,
  G, the decision for every factor, the jump set, deck images, run time).
- `results/S23_sheet3_W8_tamper.json`: S2 and S3 on the tampered matrix.
- `logs/03_sizes.txt`, `logs/10_s1_build.txt`, `logs/11_s1_validate.txt`, `logs/20_s2s3_W6.txt`, `logs/21_s2s3_W8.txt`,
  `logs/30_s4_tamper.txt`, `logs/31_s3_extra.txt`, `logs/33_s3_crosscheck.txt`, `logs/34_s3_points.txt`: the standard
  output of the scripts, in the order run.sh runs them.

## How to re-run

From the repository root:

    bash results/lab/jumping_loci_records/stratification_w6_w8/run.sh OUTDIR
    bash results/lab/jumping_loci_records/stratification_w6_w8/compare.sh OUTDIR

OUTDIR must not exist yet (its parent must). Set `PYTHON` to choose the interpreter (default `python3`).
Run the scripts through run.sh: each script writes into its current directory, so a script run by hand from this folder
would overwrite the shipped data/ and results/.
Environment as recorded: Python 3.12.8, sympy 1.14.0 with its pure-Python ground types (gmpy2 not installed); mpmath
(a dependency of sympy, not imported directly) was not recorded, and the same interpreter has 1.3.0. python-flint and
numpy are not used. Python 3.9 or later is needed (`math.lcm`). run.sh sets `SYMPY_GROUND_TYPES=python` unless it is
already set; a re-run with the flint ground types (python-flint 0.9.0) gave the same records. The engine modules import
only the standard library.

Expected runtime: about 70 seconds on one core; peak memory about 75 MB resident. run.sh limits each step to 8 GB of
virtual memory (`ulimit -v`, override with `MEM_KB`) and 900 s (`timeout`, override with `TIMEOUT_S`).

## What "reproduced" means

run.sh completes (it stops if a validation or check does not end with "ALL OK", or if the negative control is not
detected), and compare.sh prints IDENTICAL for every shipped record and the final line "ALL IDENTICAL":
- `data/*.json`: byte for byte;
- `results/*.json`: after the mask `"seconds": <number>` -> `"seconds": MASKED` (the run time of each run);
- `logs/*.txt`: after the mask `(<number> s)` -> `(MASKED s)` (run times printed by the scripts).

The negative control (`code/s4_tamper.py`) multiplies the X(1,1,1,1) coefficient of the sheet-3 I_3 density by 11/10
at W = 8; 106 of the 2059 matrix entries change. run.sh checks two ways in which it is detected, prints "TAMPER DETECTED"
for each, and stops otherwise:
- (a) S2 and S3 on the tampered matrix give another result: the generic rank goes from 28 to 29 (generic multiplicity 0)
  and G gains the factor t - 11/5. These records are shipped and compared.
- (b) The S1 validation against the engine, run in OUTDIR/control/ with the tampered matrix in place of the sheet-3,
  W = 8 matrix, fails at all 7 points of that matrix (at each, 106 entries differ from the engine's) and nowhere else. This output
  is not shipped.

The jump set {t = 3} is unchanged by this tamper.

## Provenance

An independent check made for the paper's revision (October 2026); re-run for release v1.2.8.
