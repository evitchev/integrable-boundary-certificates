"""Certify the three published I_7 formulas of hep-th/0404195
(eqs. (30) KdV, (33) O(N), (35) A_2^(2)) by solving the kernel of
ad_{I_3} (resp. ad_{I_5}) at spin 7 from scratch in the symbolic-N
invariant algebra, expressing the unique solution in the paper's Q^7
basis (eq. (26)), and reconstructing each coefficient as a function of
N (and s for A_2^(2)).

Consistency lever: KdV and O(N) coincide at N = 1, which already shows
eq. (30)'s printed 2/657 must be 2/675.
"""

from fractions import Fraction as F

from invariant_engine import (canon, gen_inv_basis, genuine_kernel_inv,
                              InvSector, kdv_P4, oN_P4, a2_P6, curve_point,
                              _rref)
from solve_charges import rat_reconstruct, rat_eval

# eq. (26) representatives of Q^7
P8B = [canon(p) for p in [
    ((1, 2), (2, 3)), ((1, 3), (1, 3)), ((1, 1), (3, 3)), ((2, 2), (2, 2)),
    ((4, 4),), ((1, 1), (1, 1), (2, 2)), ((1, 1), (1, 2), (1, 2)),
    ((1, 1), (1, 1), (1, 1), (1, 1)),
]]
NAMES = ["(23)(12)", "(13)(13)", "(33)(11)", "(22)(22)", "(44)",
         "(22)(11)(11)", "(12)(12)(11)", "(11)^4"]

sec8 = InvSector(8)
assert sec8.q_dim() == 8
cols = [sec8.reduce_mod_d(sec8.vec({b: F(1)})) for b in P8B]
rows = [[cols[c][r] for c in range(8)] for r in range(len(cols[0]))]
_, piv = _rref(rows, 8)
assert len(piv) == 8, "eq. (26) transcription is not a basis"


def express(cur):
    target = sec8.reduce_mod_d(sec8.vec(cur))
    aug = [[cols[c][r] for c in range(8)] + [target[r]]
           for r in range(len(target))]
    aug, piv = _rref(aug, 9)
    assert 8 not in piv
    sol = [F(0)] * 8
    for i, p in enumerate(piv):
        sol[p] = aug[i][8]
    return sol


def solve_I7(P_low, Nval):
    gen = genuine_kernel_inv(P_low, gen_inv_basis(8), Nval)
    assert len(gen) == 1, f"kernel dim {len(gen)} != 1 at N={Nval}"
    c = express(gen[0])
    assert c[7] != 0
    return [x / c[7] for x in c]   # normalize coefficient of (11)^4 to 1


def fit_poly(samples, maxdeg=4):
    """Reconstruct a polynomial (as rational function with denom 1)."""
    for d in range(1, maxdeg + 1):
        pq = rat_reconstruct(samples[:d + 2], d, 0)
        if pq and all(rat_eval(pq, n) == f for n, f in samples):
            p, q = pq
            return [x / q[0] for x in p]
    return None


def poly_str(p):
    terms = []
    for i, c in enumerate(p):
        if c == 0:
            continue
        t = f"{c}" if i == 0 else (f"{c}*N" if i == 1 else f"{c}*N^{i}")
        terms.append(t)
    return " + ".join(terms) or "0"


print("== KdV I_7 (certifying eq. 30) ==")
Ns = [F(k) for k in range(1, 9)]
sols = {n: solve_I7(kdv_P4(), n) for n in Ns}
for i in range(8):
    p = fit_poly([(n, sols[n][i]) for n in Ns])
    print(f"  {NAMES[i]:>14}: {poly_str(p)}")

print("\n== O(N) I_7 (certifying eq. 33) ==")
sols = {n: solve_I7(oN_P4(n), n) for n in Ns}
for i in range(8):
    p = fit_poly([(n, sols[n][i]) for n in Ns])
    print(f"  {NAMES[i]:>14}: {poly_str(p)}")

print("\n== A_2^(2) I_7 (certifying eq. 35) ==")
# sample the double cover at +-t; split c = A(N) + B(N) s
ts = [F(x) for x in (2, 3, 4, 6, 7, 8, 9, F(9, 2), F(11, 3), F(13, 4),
                     F(15, 2), F(17, 3))]
solsA, solsB = {}, {}
for t in ts:
    Np, sp = curve_point(t)
    Nm, sm = curve_point(-t)
    assert Np == Nm and sp == -sm
    cp = solve_I7(a2_P6(Np, sp), Np)
    cm = solve_I7(a2_P6(Nm, sm), Nm)
    solsA[Np] = [(a + b) / 2 for a, b in zip(cp, cm)]
    solsB[Np] = [(a - b) / (2 * sp) for a, b in zip(cp, cm)]
for i in range(8):
    sampA = [(n, v[i]) for n, v in solsA.items()]
    sampB = [(n, v[i]) for n, v in solsB.items()]
    pA = fit_poly(sampA)
    pB = fit_poly(sampB)
    strA = poly_str(pA) if pA else "NOT POLYNOMIAL (rational?)"
    strB = poly_str(pB) if pB else "NOT POLYNOMIAL (rational?)"
    print(f"  {NAMES[i]:>14}: [{strA}]  +  s * [{strB}]")
