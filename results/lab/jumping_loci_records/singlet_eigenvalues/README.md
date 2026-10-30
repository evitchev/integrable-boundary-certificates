# singlet_eigenvalues

Exact identities behind the free-field (singlet-algebra) reading of the fusion-product law in the paper
"Jumping loci of the cylindrical commuting families: arithmetic ladders on the curve s^2=(N-25)(N-1)" (October 2026).

## What the computation shows

The paper's fusion-product law gives the vacuum eigenvalue of the extra sheet-3 class at beta^2 = 1/q, q = 2m
(confirmed there at q = 4, ..., 14), in the form E = prod_{j=-(m-1)}^{m-1} (P_X + pi - j p)(P_X - pi - j p), with
p^2 = 1/m. `code/b5_fusion.py` proves in exact
symbolic arithmetic, for m = 2, ..., 8 (q = 4, ..., 16), that E = F_m(P_X + pi) F_m(P_X - pi) with
F_m(z) = prod_{j=-(m-1)}^{m-1} (z - j p) = p^(2m-1) (2m-1)! binom(z/p + m - 1, 2m - 1), and that
binom(2P + m - 1, 2m - 1) = 2^(2m-1)/(2m-1)! prod_{j=-(m-1)}^{m-1} (P - j/2). It then takes the two factor eigenvalues
recorded at q = 4 and q = 6, (4/3) P (P^2 - 1/4) and (1/15) P (P^2 - 1)(4 P^2 - 1), and shows that they equal
binom(2P + 1, 3) and binom(2P + 2, 5) exactly, with no free constant: the eigenvalue binom(t, 2m - 1) of the generator H
of the singlet algebra W(2, 2m - 1) on the top level of a Fock module (D. Adamović, J. Algebra 270 (2003) 115-132,
arXiv:math/0207155, Prop. 5.2), here with t = 2P + (m - 1); the products of the two factors reproduce the recorded
eigenvalues -(2/9) E and -(81/4800) E of the extra class at q = 4 and q = 6. `code/b5_vacuity.py` measures how much such
an eigenvalue match can say: with exact fractions it computes the rank of the cylinder vacuum-eigenvalue map on one-boson
fields of weight w = 3, 5, 7, 9 (pairings C(a,b) = (-1)^b R^(a+b)(0), R(eps) = 1/(4 sinh^2(eps/2)) - 1/eps^2; every
unpaired letter must be dV and gives P; kappa = 1). The ranks are 1, 2, 3, 4, inside spaces of odd polynomials of
dimension 2, 3, 4, 5, on 1, 2, 4, 8 field classes modulo total derivatives. The same script recomputes the q = 4 factor
eigenvalue (4/3) P^3 - (1/3) P from the recorded weight-3 density.

## Statements of the paper that rest on it

Section 4 ("The fusion-product law"), paragraph "A free-field reading":

- "F_m(z) = prod_{j=-(m-1)}^{m-1}(z - jp) = p^{2m-1}(2m-1)! binom(z/p+m-1, 2m-1), p^2 = 1/m [lab result: exact at
  q=4,6,...,16]", with the splitting "F_m(P_X+pi) F_m(P_X-pi)": part (i) of `b5_fusion.py`.
- "The recorded factor eigenvalues are binom(2P_+ + 1, 3) at q=4 and binom(2P_+ + 2, 5) at q=6, exactly, so the factors
  are Adamović's H at m=2,3 [lab result]": part (ii) of `b5_fusion.py`, which proves the equality with the binomials; the
  q = 4 eigenvalue is also recomputed from its density at the end of `b5_vacuity.py`.
- "The eigenvalue match alone is weak evidence: at each weight >= 5 it is one linear condition on the density": the
  ranks of `b5_vacuity.py`.

Not resting on this folder: "the span of 49 tested products of a field of U with a field of W meets the kernel modulo
total derivatives in a one-dimensional subspace, spanned by H_U H_W [lab result]" is test P2 of the lab record
`opers-2026-10/tree/results/lab/large_n/wuw2/` (`wuw2.json`, `run_wuw2.log`); the q = 8 statements of Section 4 rest on
a separate computation.

Inputs. The scripts read no file. The factor eigenvalues and the product constants compared in part (ii) are typed in
as printed by the lab records `opers-2026-10/tree/results/lab/large_n/wuw1/run_wuw1.log` and
`opers-2026-10/tree/results/lab/large_n/wuw2/run_wuw2.log` (tests P1 and P3; see `NOTE.md` there). The weight-3 density
(4/3)(dU)^3 + 2 dU d^2U + (1/3) d^3U of the cross-check is `W_U` of `opers-2026-10/tree/results/lab/large_n/wuw1/wuw1.json`.

Notation. The scripts' P_U (q = 4) and P_V (q = 6) are the paper's P_+ = (sqrt(m)/2)(P_X + pi); their P_W and P_Z are
(sqrt(m)/2)(pi - P_X) = -P_-. In part (i) p = 1/sqrt(m); in part (ii) "p_Ad" is the index m of W(2, 2m - 1), so the
"shift p-1" printed there is m - 1. E_i, i = m - 1, is the product E above.

## Files

- `code/b5_fusion.py`: parts (i) and (ii) and the root-string identity; exact sympy.
- `code/b5_vacuity.py`: part (iii), the ranks of the eigenvalue map at weights 3, 5, 7, 9, and the q = 4 cross-check;
  exact fractions (sympy only for the Taylor coefficients of R).
- `logs/b5_fusion.log`: the archived output of `b5_fusion.py` (every check PASS, "ALL PASS", exit 0).
- `logs/b5_vacuity.log`: the archived output of `b5_vacuity.py` (ranks 1, 2, 3, 4; cross-check PASS; exit 0).
- `run.sh`: re-runs both scripts into a new directory.
- `compare.sh`: compares the new logs with `logs/`.
- `README.md`: this file.

## How to re-run

From the repository root, with a new output directory outside the repository:

    PYTHON=python3 bash results/lab/jumping_loci_records/singlet_eigenvalues/run.sh ../singlet_eigenvalues_rerun
    bash results/lab/jumping_loci_records/singlet_eigenvalues/compare.sh ../singlet_eigenvalues_rerun

`run.sh` refuses an output directory that already exists or lies inside this folder or the repository. Each script runs
as `$PYTHON -I -B code/<name>.py` under `ulimit -v 8000000` (KiB) and `timeout 900` (seconds); its log is its stdout and
stderr followed by the line `exit <status>`.

Environment: bash with `timeout` (GNU coreutils) and `sha256sum` or `shasum`; Python with sympy and mpmath, as in the
repository's `requirements.txt`. Recorded environment (the interpreter of the archived run and of the re-run): Python
3.12.8 with sympy 1.14.0, mpmath 1.3.0, python-flint 0.9.0 and numpy 2.5.2 installed. The scripts import only sympy
(which uses mpmath) and the standard library; sympy uses python-flint for its integer arithmetic when it is installed,
and the logs are identical with sympy's pure-Python arithmetic (`SYMPY_GROUND_TYPES=python`). numpy is not used.

Expected cost: about 10 seconds in all (about 7 s and 2 s), peak resident memory about 90 MB.

## What "reproduced" means

`run.sh` exits 0, and `compare.sh` prints `IDENTICAL` for both logs and then `ALL IDENTICAL`.

Mask: the first line of each log is `script sha256 <hash>`, the sha256 of the script file that wrote it. The shipped
scripts differ from the ones that wrote the shipped logs in their first docstring line only (an internal label was
removed for publication), so the shipped logs carry the hashes of the scripts as run (d7e70fc7... and 0cfec3f5...).
`compare.sh` replaces this line by `script sha256 <masked>` in both logs before comparing them, and requires instead that
the first line of each new log equal the sha256 of the shipped script. Every other line, including the final `exit 0`,
must match byte for byte. The logs contain no times or paths.

## Provenance

An independent check made for the paper's revision (October 2026); re-run for release v1.2.8.
