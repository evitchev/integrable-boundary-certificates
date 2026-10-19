# Runs made when this directory was assembled (2026-10-18)

All runs: `run_checks.py` from this directory, HOME set to an empty directory, `ulimit -v 8000000`.

| run | interpreter | result | file |
|---|---|---|---|
| full (checks 1 and 2) | Python 3.8.10, sympy 1.13.3 | ALL PASS (check 2: 2907 s) | runs/full_py3.8.10.out |
| full (checks 1 and 2) | Python 3.12.8, sympy 1.14.0 | ALL PASS (check 2: 2005 s) | runs/full_py3.12.8.out |
| drill (Ramond kernel + eps/100) | Python 3.12.8, sympy 1.14.0 | both tampered runs FAIL, as required | runs/drill_py3.12.8.out |
| drill (psi_0^2 = 1/3), superseded | Python 3.8.10, sympy 1.13.3 | NOT detected by either check: the value of psi_0^2 does not enter the vacuum eigenvalues (README, Conventions). Replaced by the kernel tamper; kept for the record | runs/drill_psi0sq_py3.8.10_VACUOUS.out |

The two full runs used the first version of `run_checks.py`, which differs from the shipped one only in the drill tamper
(the TAMPER line and the docstring); the checks themselves are unchanged.
