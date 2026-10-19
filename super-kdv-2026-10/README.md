# N = 1 super-KdV integrals of motion: cylinder vacuum eigenvalues, NS and Ramond sectors

This directory holds code and results by E. Vitchev (2026). The code computes the local integrals of motion I_1, I_3, I_5, I_7 and
I_9 of quantum N = 1 super-KdV (the one-parameter "lambda" family of the N = 1 superconformal free-field realisation), and their
eigenvalues on the Fock highest-weight states on the cylinder, in both the Neveu-Schwarz (NS) and the Ramond (R) sectors.
Everything is in exact rational arithmetic, symbolic in the momentum and in the background charge.

The computation was made on 2026-10-12. It was an independent benchmark of the NS- and R-sector eigenvalues given in N. Tanabe,
arXiv:2604.14899 (JHEP 08 (2026) 150). No printed value from that or any other paper is used to produce the results here.

## Results

`results/lab/susy4/tanabe/computed_eigenvalues.txt` (human-readable) and `computed_eigenvalues.json` give I_3, I_5, I_7 and I_9 for
NS and R:
- in the variables (P, l) of the evaluator;
- in (Delta, c), with Delta the conformal weight of the state in its own sector.

I_1 is L_0 - c/24.

Hashes of the two files at the time of computation:

| file | sha256 |
|---|---|
| computed_eigenvalues.json | 88005e3b0618b871436aceff932d65f4204629a7694f38268659e6645c6887f3 |
| computed_eigenvalues.txt | 334fec4a838c5afeb74740efa1d9231550d39ea91d9f67ea56f0699822c7c111 |

## Conventions

- Free fields: a boson X and a Majorana fermion psi (`code/super_lambda.py`).
  - T_l = (1/2) (dX)^2 + l d^2X - (1/2) psi dpsi and G_l = psi dX + 2 l dpsi, with Res(G_l, G_l) = 2 T_l.
  - Central charge: c = 3/2 - 12 l^2.
- P is the zero mode of dX - l on the cylinder: P = p - l for the Fock momentum p.
  - In the notation of arXiv:2604.14899: l = -alpha_0/2 and P = Lambda + alpha_0/2, so P^2 = sigma_1.
- NS sector: L_0 - c/24 = P^2/2 - 1/16.
- R sector: L_0 - c/24 = P^2/2, so Delta_R = Delta_NS(p) + 1/16.
  - The ground space is two-dimensional; the fermion zero mode psi_0 acts as a Clifford generator with psi_0^2 = 1/2, and G_0 = P psi_0.
    In the vacuum eigenvalues of the even densities I_{2k-1}, no product of two unpaired zero modes survives, so the value of
    psi_0^2 does not enter them (a tamper psi_0^2 = 1/3 leaves every result here unchanged); the Ramond kernel does.
- Normalisation: I_3 = :T_l T_l: - (1/4) :G_l dG_l:. For k >= 3, I_{2k-1} is normalised so that the coefficient of (dX)^{2k} is 1/2^k;
  its eigenvalue then starts (P^2/2)^k.

## Method

- `code/super_cyl_vev.py`: the cylinder evaluator. It maps a normal-ordered density to its zero-mode eigenvalue. The kernels are
  derived in its docstring from the mode expansions:
  - NS: 1/(2 sinh(eps/2)) - 1/eps;
  - R: 1/(e^eps - 1) - 1/eps, with psi_0 kept as an operator.
- `code/super_cyl_vev_validate.py`: the evaluator's validation.
  - T_l gives L_0 - c/24 in both sectors, and G_0^2 = L_0 - c/24 in R.
  - Bosonic monomials agree with the independent bosonic evaluator `code/cyl_vev.py`.
  - Total derivatives evaluate to zero.
  - Five negative controls must fail.
- `code/super_cyl_tower.py`: at a rational l, I_{2k-1} is the joint kernel of ad G_l and ad I_3(l) at weight 2k, modulo total
  derivatives. The kernel is asserted to be one-dimensional.
- `code/super_cyl_tower_fit.py`: the symbolic computation.
  - The tower is computed at 61 rational values of l and brought to a canonical representative modulo total derivatives.
  - Every coefficient is reconstructed as a rational function of l: fitted on 40 samples, and all 61 are reproduced.
  - The commutation conditions are then re-checked symbolically in l: [G_l, I] and [I_3, I] are total derivatives for I_5, I_7, I_9,
    and so are [I_5, I_7] and [I_5, I_9].
  - The eigenvalues are computed from the symbolic densities.
- `code/super_engine.py`, `super_lambda.py`, `super_screening_fit.py`, `ope_engine.py`, `sparse_linalg.py`, `cyl_vev.py`: the OPE and
  linear-algebra engines. They are byte-identical to the copies in `code/` at the root of this repository.

Limitations:
- The one-dimensionality of the kernels is established at the 61 sample values of l, not symbolically.
- [I_7, I_9] is a total derivative at three exact held-out values of l; the symbolic check was not completed.
- The eigenvalues are vacuum (Fock highest-weight) eigenvalues only.

## Re-running

    python3 run_checks.py --quick     # check 1: evaluator validation (about 15 s)
    python3 run_checks.py             # checks 1 and 2: also recomputes the tower and the eigenvalues (about 40 min)
    python3 run_checks.py --drill     # the same with a deliberate tamper (Ramond kernel + eps/100); every check must then FAIL

`run_checks.py` first verifies every file against `SHA256SUMS`, then runs each check in a temporary copy, so nothing here is
overwritten.
- Check 1 compares the validation output with `logs/super_cyl_vev_validate.log`, with run times removed.
- Check 2 requires the regenerated `computed_eigenvalues.json` and `.txt` to be byte-identical to the archived files.
- `logs/super_cyl_tower_fit_w20.log` is the output of the original run (2321 s).

Requirements: Python 3 with sympy (tested: Python 3.8.10 with sympy 1.13.3, and Python 3.12.8 with sympy 1.14.0; see `logs/RUNS.md`).

## License and citation

MIT, as the rest of the repository. If you use these results, please cite the repository at its release tag (v1.2.6 or later).
