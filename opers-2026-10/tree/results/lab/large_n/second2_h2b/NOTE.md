# SECOND2 sub-case H2b (Sol 3, rotated dictionary): exact decision. Lead batch, 2026-10-02

The seat directory <home>/qwen-work/second2/ was not modified. Inputs were copied here and checked against the seat's
SHA256SUMS (all 16 OK, see `SHA256SUMS.seat`). Seat NOTE kept as `NOTE.seat.md`.

## Verdict: H2b FAIL at t = 2 and at t = 9/4 (exact)

Dictionary (as coded in run_second2.py, `rotated=True`): A = alpha^2 = sX (X + pi_rot)^2 + cX, lam = sY Y + cY, M = (t+3)/(t-1).
Fit (as coded): spin-5 monomials X^2Y, X^2, XY, X plus the spin-7 monomial X^3Y. That is 5 equations in 5 unknowns.
Note that the task text says "from spin 5", but the code also takes one equation from spin 7.

| fibre | raw ideal | saturated ideal (by sX * all fit denominators) | solutions | real | result |
|---|---|---|---|---|---|
| t=2, M=5 | dim 1: 3 minimal primes; the two 1-dim ones lie in sX=0 (excluded) | dim 0, vdim 12 = radical vdim, ONE prime: irreducible deg-12 in pi_rot, shape position | 12 (one Galois orbit) | 6 | all 12 FAIL |
| t=9/4, M=21/5 | dim 1 | dim 0, vdim 12, one prime, deg 12 irreducible | 12 | 6 | all 12 FAIL |

- **First failure.** It comes at spin 5 itself, on monomials the fit did not use:
  - In the certified support (total degree <= 3), X Y^2 is EXACTLY nonzero on the prime. Numerically |monic - cert| is 1.51 to 117 over the 12 roots at t=2, and 1.41 to 23.5 at t=9/4.
  - Also nonzero: Y^3, Y^2, Y, 1.
  - Spin 7 (first: X^2Y^2) and spin 9 (first: X^4Y) also fail exactly.
  - The fitted monomials (spin-5 X^2Y, X^2, XY, X and spin-7 X^3Y) all come out zero, which is a sanity check that the solutions are right.
- **Structural reason (independent of the solve).** Because A is quadratic in X, the mapped spin-(2k-1) charge has X-degree 2k.
  - Its X^(2k) coefficient after the monic normalisation is c_k sX^k / lead.
  - c_k is the A^k coefficient of the WKB charge. At M=5 it is 1/5120, 299/20480000 and 2387/4096000000 for k=3,4,5, all nonzero (also nonzero at M=21/5).
  - The certified tables have total degree k, so the rotated dictionary can never match unless sX = 0, which is excluded.
  - So H2b fails on monomials X^4, X^5, X^6 (spin 5) even before the in-support failures. The in-support failures above show that it fails even if one ignored that.
- **Exactness.** Each saturated prime is in shape position over Q with an irreducible univariate polynomial. Every residual was reduced exactly in Q[z]/(m(z)), so "nonzero" means nonzero at every one of the 12 conjugate solutions (Galois argument). The 60-digit mpmath values at all 12 roots are only an extra numeric display.
- **msolve cross-check** (saturation via 1 - zs*h). At t=2, t=9/4 and t=2 with the tamper: zero-dimensional, parametrisation degree 12, 6 real solutions. This agrees with Singular.

## Controls (same pipeline: gen_eqs.py, then Singular, then check.py)
- **Gate-1 Sol 1 plant** (gate1.py's own cert/wkb_mapped, spin-3 fit, t=2, M=5): PASS. There are 2 rational solutions:
  - {sX=-120, sY=-6, cX=0, cY=15/4} matches spins 5, 7, 9 exactly.
  - {120, 6, 45, -5} has 6/10/15 mismatches, the same as the seat's log.
- **Gate-1 tamper** M=51/10: FAIL (2 solutions, both fail). This changes the outcome.
- **H2a control** (standard dictionary): FAIL at t=2 with 2 solutions (quadratic field; agrees with the seat's sqrt(1281)).
  - First failures: spin 5 X Y^2 with |.|=240/49=4.898, spin 7 X^3Y with 16/253, spin 9 X^4Y with 24/217, as in the seat's log.
  - The run at t=9/4 also FAILs.
- **Rotated-pipeline plant** (plantcert.py: synthetic table from the seat's own wkb_poly + map_and_monic(rotated=True) at sX=-120, sY=-6, cX=7, cY=15/4, pi=2, M=5):
  - The saturated ideal splits into a degree-1 prime equal to the planted point, which matches ALL of spins 5, 7, 9 exactly, plus a degree-11 prime that fails. Result: PASS.
  - With the tamper M -> 51/10 it becomes one degree-12 prime with no match: FAIL. This shows the H2b code path can detect a true positive, so the FAIL above is not vacuous.
- **H2b tamper** M=51/10 at t=2: 12 solutions, still FAIL. This is expected, because a FAIL cannot be "changed" by a tamper; the plant pair above is the discriminating control.

## Disclosures
- Equations were rebuilt by importing the seat's run_second2.py (load_cert, wkb_poly, map_and_monic). For gate1, gate1.py was exec'd up to its `plant_ok =` line, to get its functions without running it. Nothing in the WKB was reimplemented.
- All 5 t=2 equations match the (truncated) printouts in run_second2_h2b.log as prefixes. eq[2] matches in full.
- Bug found and fixed in my first Singular script: in Singular 4.4.1, `sat()` returns an ideal, not a list, so `sat(I,h)[1]` silently took only the first generator. The fixed runs print typeof(sat) = ideal.
- msolve 0.10.1 misparses parenthesised products. It reported a spurious positive dimension until the input was fully expanded, and all msolve inputs here are expanded. This is worth remembering for other batches.
- The raw (unsaturated) prime decomposition was moved to the end of the Singular scripts as informational only. The 4 runs of it (t=2 tamper, t=9/4, both plants) were killed by saved PID (kills.log) after the saturated results had been written. The t=2 raw decomposition did complete: 3 primes, as above.
- sing_h2b_t2/h2a_*/gate1_* .sing files use the earlier layout (raw decomposition first). The math is the same.
- The seat's SEAL lets M be free per fibre. H2b as coded fixes M = (t+3)/(t-1), and that is what was decided here.

## Files
gen_eqs.py, check.py, plantcert.py; per tag: gen_eqs_*.log, eqs_*.json, sing_*.sing/.log, comp_*_j.txt (lex GB of each saturated prime), msolve_*.ms/.out/.log, check_*.log; pids.txt, kills.log, SHA256SUMS.
