"""Exact (number-field) + 60-digit numeric check of every solution of the saturated fit ideal against the certified tables,
using the seat's own wkb_poly/load_cert (run_second2.py) or cert (gate1.py), same dictionary map and same monic normalisation
(divide by the coefficient of X^k).  Each saturated minimal prime is zero-dimensional in shape position (lex, last variable
univariate & irreducible), so a rational-coefficient residual vanishes at one of its points iff it vanishes at all (Galois).
Usage: python3 check.py MODE t [Mshift]"""
import sys, json; sys.dont_write_bytecode = True
from fractions import Fraction as F
import sympy as sp, mpmath as mp
import run_second2 as R
from run_second2 import X, Y, A, L, sX, sY, cX, cY, pi_rot, load_cert, wkb_poly
mp.mp.dps = 60
import os
if os.environ.get('PLANT'):
    import plantcert; load_cert = plantcert.load_cert
mode = sys.argv[1]; tv = F(sys.argv[2]); shift = F(sys.argv[3]) if len(sys.argv) > 3 else F(0)
M = (tv + 3) / (tv - 1) + shift
tag = f"{mode}_t{str(tv).replace('/','o')}" + (f"_shift{str(shift).replace('/','o')}" if shift else "") + ("_plant" if os.environ.get('PLANT') else "")
rot = (mode == 'h2b')
unk = [sX, sY, cX, cY] + ([pi_rot] if rot else [])
if mode == 'gate1':
    G = {}; exec(open('gate1.py').read().split('plant_ok =')[0], G); cert = lambda k: G['cert'](k, tv)
else:
    cert = lambda k: load_cert(3, k, tv)
ncomp = int(open(f'ncomp_{tag}.txt').read().strip())
loc = {str(u): u for u in unk}
print(f"=== CHECK {tag}: M = {M}, rotated = {rot}, saturated minimal primes = {ncomp}")
WK = {k: sp.Poly(wkb_poly(M, k)[0], A, L) for k in (3, 4, 5)}
CK = {k: cert(k) for k in (3, 4, 5)}
overall = []
for j in range(1, ncomp + 1):
    gens = [sp.sympify(g.replace('^', '**'), locals=loc) for g in open(f'comp_{tag}_{j}.txt').read().replace('\n', '').split(',') if g.strip()]
    z = unk[-1]
    uni = [g for g in gens if g.free_symbols == {z}]; assert len(uni) == 1, gens
    m = sp.Poly(uni[0], z, domain='QQ'); d = m.degree()
    irr = len(sp.factor_list(m.as_expr())[1]) == 1 and sp.factor_list(m.as_expr())[1][0][1] == 1
    val = {z: sp.Poly(z, z, domain='QQ')}
    for u in unk[:-1]:
        g = [g for g in gens if u in g.free_symbols]; assert len(g) == 1, (u, gens)
        P = sp.Poly(g[0], u, z, domain='QQ'); assert P.degree(u) == 1 and set(P.as_expr().free_symbols) <= {u, z}
        c1 = P.coeff_monomial(u); assert c1.is_number
        val[u] = sp.Poly(-(g[0] - c1 * u) / c1, z, domain='QQ').rem(m)
    print(f"\n--- prime {j}: degree {d} in {z} (irreducible over Q: {irr}); real roots: {sp.Poly(m, z).count_roots()}")
    roots = [r for r in sp.Poly(m.as_expr(), z).nroots(n=60, maxsteps=500)]
    if d <= 2:
        for u in unk: print(f"    {u} = {sp.Poly(val[u], z).as_expr()}   (mod {m.as_expr()})")
    red = lambda p: p.rem(m)
    # A, L as polynomials in X (resp Y) with coefficients in K = Q[z]/(m)
    one = sp.Poly(1, z, domain='QQ'); zero = sp.Poly(0, z, domain='QQ')
    def kmul(P, Q):
        o = {}
        for (a, b), u in P.items():
            for (c, e), v in Q.items():
                o[(a + c, b + e)] = red(o.get((a + c, b + e), zero) + u * v)
        return {k: v for k, v in o.items() if not v.is_zero}
    if rot:
        p = val[pi_rot]
        Apol = {(2, 0): val[sX], (1, 0): red(2 * p * val[sX]), (0, 0): red(p * p * val[sX] + val[cX])}
    else:
        Apol = {(1, 0): val[sX], (0, 0): val[cX]}
    Lpol = {(0, 1): val[sY], (0, 0): val[cY]}
    Apow = [{(0, 0): one}]; Lpow = [{(0, 0): one}]
    for _ in range(6): Apow.append(kmul(Apow[-1], Apol)); Lpow.append(kmul(Lpow[-1], Lpol))
    prime_ok = True
    for k in (3, 4, 5):
        mapped = {}
        for (a, l), c in WK[k].terms():
            for key, v in kmul(Apow[a], Lpow[l]).items():
                mapped[key] = red(mapped.get(key, zero) + v * sp.Rational(c))
        lead = mapped.get((k, 0), zero)
        if lead.is_zero: print(f"  spin {2*k-1}: lead = 0 identically on this prime -> dictionary invalid"); prime_ok = False; continue
        Cp = sp.Poly(CK[k], X, Y); Cd = {mm: sp.Rational(c) for mm, c in zip(Cp.monoms(), Cp.coeffs())}
        keys = sorted(set(mapped) | set(Cd), key=lambda t: (-(t[0] + t[1]), -t[0]))
        res = {key: red(mapped.get(key, zero) - lead * Cd.get(key, 0)) for key in keys}   # = lead*(monic - cert)
        bad = [key for key in keys if not res[key].is_zero]
        bad_in = [key for key in bad if key[0] + key[1] <= k]  # inside the certified table's support (total degree <= k, as for H1/H2a)
        bad_out = [key for key in bad if key[0] + key[1] > k]  # outside it: must vanish for a match
        # numerics at every root: (monic - cert) value
        def num(key, r): return complex(sp.Poly(res[key], z).eval(r) / sp.Poly(lead, z).eval(r))
        print(f"  spin {2*k-1}: {len(keys)} monomials compared, EXACT nonzero residuals: {len(bad)} "
              f"(outside cert support, tot.deg>k: {len(bad_out)} {bad_out[:4]}; inside, tot.deg<=k: {len(bad_in)} {bad_in[:6]})")
        if bad_in:
            key = bad_in[0]
            mags = [abs(num(key, r)) for r in roots]
            print(f"     first in-support failing monomial X^{key[0]}Y^{key[1]}: |monic-cert| over the {d} roots: min {min(mags):.3e} max {max(mags):.3e}")
        if bad_out:
            key = bad_out[0]
            mags = [abs(num(key, r)) for r in roots]
            print(f"     first out-of-support monomial X^{key[0]}Y^{key[1]}: |monic coeff| over roots: min {min(mags):.3e} max {max(mags):.3e}")
        if bad: prime_ok = False
    print(f"  => prime {j} ({d} solutions): {'ALL MATCH' if prime_ok else 'FAILS'}")
    overall.append(prime_ok)
print(f"\nRESULT {tag}: {'PASS (some solution matches all spins)' if any(overall) else 'FAIL (no solution matches)'}")
