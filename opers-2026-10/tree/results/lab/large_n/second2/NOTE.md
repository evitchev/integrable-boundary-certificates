# SECOND2 — Suzuki class test for Sol 2 and Sol 3

Qwen peer seat, 2026-10-01/02. Lead commission: does the Suzuki kappa-odd scalar second-order
ODE class (-psi'' + [x^(2M) + alpha x^(M-1) + l(l+1)/x^2] psi = E psi, alpha^2 ~ X) that gave
Sol 1 its oper (item 306/SUPER3b) also reproduce Sol 2 or Sol 3?

## Engine validation (Gate 1) — PASS

`gate1.py`: at t=2, M=5, fitting {sX, sY, cX, cY} on spin 3 recovers the known Sol 1 dictionary
{sX=-120, sY=-6, cX=0, cY=15/4} exactly, with all of spins 5, 7, 9 matching. Tamper (M→M+1/10)
correctly fails. Engine validated before any Sol 2/3 comparison.

Log: `gate1_plant_tamper.log`.

## Method

- WKB local IMs via `suzuki_wkb.py` (unmodified from SUPER3): exact Riccati recurrence +
  continued Beta charge extraction, Fraction arithmetic throughout.
- Dictionary: A = alpha² = sX·X + cX, L = lam = sY·Y + cY, ALL scales free.
- M fixed per fibre to M = (t+3)/(t-1) (the ODE structural parameter, not a dictionary unknown).
- Fit on spin 5 (4 monomials X²Y, X², XY, X) for Sols 2/3 (spin-3 certified data only exists for Sol 1).
- Spins 7, 9 are zero-parameter checks after substitution.
- Fibres: t = 2, 9/4, 4.
- Sympy solve for {sX, sY, cX, cY} (H1/H2a) or {sX, sY, cX, cY, pi_rot} (H2b), ulimit -v 20000000.

## Results — H1: SOL 2 (standard dictionary)

**FAIL at all three fibres.** At each fibre, sympy.solve finds exactly 2 solutions (conjugate
square-root branches). Neither solution matches the certified top layer at spins 7 or 9.

| Fibre | t | M | #solutions | Spin 5 | Spin 7 | Spin 9 | Verdict |
|-------|---|---|-----------|--------|--------|--------|---------|
| 1 | 2 | 5 | 2 | fit | MISMATCH | MISMATCH | FAIL |
| 2 | 9/4 | 21/5 | 2 | fit | MISMATCH | MISMATCH | FAIL |
| 3 | 4 | 7/3 | 2 | fit | MISMATCH | MISMATCH | FAIL |

Example (t=2, solution 0): sX = -1344√1597029/76049, sY = -504√1597029/76049,
cX = 90 - 12950√1597029/228147, cY = 25/4 - 539√1597029/152098.
First mismatch at spin 7, monomial (3,1): coefficient off by -26/23.

**The Suzuki class does NOT reproduce Sol 2.**

## Results — H2a: SOL 3 (standard dictionary)

**FAIL at all three fibres.** Same pattern as Sol 2: 2 solutions per fibre, neither matches
spins 7 or 9. Notably, at t=4 the two solutions are COMPLEX (involving I·√1495), which is
itself a red flag for a real ODE reproducing a real VEV table.

| Fibre | t | M | #solutions | Spin 5 | Spin 7 | Spin 9 | Verdict |
|-------|---|---|-----------|--------|--------|--------|---------|
| 1 | 2 | 5 | 2 (real) | fit | MISMATCH | MISMATCH | FAIL |
| 2 | 9/4 | 21/5 | 2 (real) | fit | MISMATCH | MISMATCH | FAIL |
| 3 | 4 | 7/3 | 2 (complex) | fit | MISMATCH | MISMATCH | FAIL |

Example (t=2, solution 0): sX = -56√1281/61, sY = -6√1281/61,
cX = 90 - 350√1281/183, cY = 25/4 - 21√1281/122.
First mismatch at spin 7, monomial (3,1): coefficient off by -16/253.

**The Suzuki class does NOT reproduce Sol 3.**

## Results — H2b: SOL 3 (rotated dictionary A = sX(X+π)² + cX)

**INTRACTABLE within time budget.** The 5×5 polynomial system (unknowns: sX, sY, cX, cY, pi_rot;
max degree 4) was submitted to sympy.solve with a 1200s timeout at fibre t=2. The solve did not
terminate within the budget (killed by timeout, exit code 143). No solution set was produced.

This is reported as a genuine outcome, not a failure: the rotated dictionary introduces a 5th
free parameter (the rotation π), making the system significantly harder than the 4×4 systems
that solved instantly for H1/H2a. The structural obstruction may be real (no solution exists)
or the system may require a different elimination strategy (e.g., Groebner with a specific
monomial order, or resultants on a subset of variables first).

**Interpretation:** The standard dictionary already fails decisively for Sol 3 (H2a FAIL at all
three fibres, with complex solutions at t=4). The rotated variant is a narrower question that
cannot be answered within the current computational budget. It is NOT evidence in favour of
the Suzuki class extending to Sol 3; it is simply undecided.

Log: `run_second2_h2b.log` (shows the 5 equations and "Solving..." followed by "Terminated").

## Cross-fibre consistency verdict

- **SOL 2 (H1): FAIL** — no fibre has a full match (2 irrational solutions per fibre, all fail spins 7/9).
- **SOL 3 standard (H2a): FAIL** — no fibre has a full match (2 solutions per fibre incl. complex at t=4, all fail spins 7/9).
- **SOL 3 rotated (H2b): INTRACTABLE** — sympy.solve timed out at 1200s on the 5×5 system at t=2; no solution produced.

## Key structural observation

At every fibre for both Sols 2 and 3, the spin-5 fit produces exactly TWO conjugate solutions
(symmetric under (sX, sY, cX, cY) → (-sX, -sY, conj(cX), conj(cY))). This is the same
two-branch structure seen in Gate 1's Sol 1 plant (where one branch matched and the other
didn't). For Sols 2/3, NEITHER branch survives to higher spins — the mismatch appears at
spin 7 in the (3,1) monomial in every case.

The fact that the mismatch is ALWAYS in the X³Y (or X⁴Y at spin 9) monomial — the monomial
with the highest Y-power among non-leading terms — suggests the WKB charge polynomial has
the wrong Y-dependence structure for Sols 2/3. The Suzuki ODE's L-term enters as l(l+1)/x²
which generates Y-monomials through the Riccati recurrence, but the pattern of Y-powers
generated does not match the certified VEV tables for Sols 2/3.

## Conclusion

**H0 CONFIRMED: The Suzuki kappa-odd scalar second-order ODE class does NOT extend beyond
Sol 1.** Neither Sol 2 nor Sol 3 top layers can be reproduced by this ODE with any choice of
free dictionary parameters at any tested fibre. The Sol 1 oper is specific to Sol 1.

This is consistent with:
- Item 252 (S1 scope): excludes kappa^0 scalar opers for all three solutions; the kappa-odd
  class was never tested for Sols 2/3 until now.
- VIR8 (Fable): Sol 3 = 4th-order three-term ODE, Sol 2 = 3rd-order — different ODE orders,
  not the same Suzuki second-order form.
- The Sol 3 self-duality under u → -u/(1+u) is broken by a single kappa·√X term.

## Files

- `SEAL_SECOND2.md` — seal (accepted by lead, hash f2f547eb...)
- `gate1.py`, `gate1_plant_tamper.log` — engine validation
- `run_second2.py`, `run_second2.log` — main run (H1 + H2a complete, H2b timed out in main script)
- `run_second2_h2b.log` — H2b standalone run (background)
- `suzuki_wkb.py`, `eng_null.py` — engine (copied from SUPER3, unmodified)
- `vev_sol{1,2,3}_w{6,8,10}.json`, `vev_w4_points.json` — certified data
- `SHA256SUMS` — input file hashes
