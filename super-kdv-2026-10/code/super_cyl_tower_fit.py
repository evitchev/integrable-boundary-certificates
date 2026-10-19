"""Track T2: super-KdV tower of the lambda family SYMBOLIC in l, and its NS / R
cylinder eigenvalues on the Fock highest-weight state, symbolic in P and l.

Route (exact throughout):
 1. at sampled rational l: joint (ad_{G_l}, ad_{I3(l)}) kernel at w2 = 12, 16,
    20 (asserted one-dimensional), normalised to (dX)^{2k} coefficient 1/2^k
    (super_cyl_tower.tower), then brought to a CANONICAL representative mod d:
    full reduction against the l-independent integer d-image pivots of the
    sector (max-index pivots, every pivot index eliminated);
 2. each canonical coefficient is reconstructed as a rational function of l
    (exact linear solve, degrees increased until every held-out sample
    agrees; the fit uses the first N_FIT samples, the rest are held out);
 3. PROOF STEP: the reconstructed symbolic density is checked symbolically in
    l: Res(G_l, I) and Res(I3(l), I) reduce to zero modulo the d-image
    (so the commuting property holds identically in l, not only at samples);
 4. eigenvalues: super_cyl_vev.vev_density on the symbolic density, NS and R
    (P = zero mode of J_rec = p - l).  I3 is built directly symbolically.

Output: results/lab/susy4/tanabe/computed_eigenvalues.json (+ .txt), the file
registered in PREREG.md BEFORE any comparison with printed values.  No printed
value from any paper is used here.
Usage: super_cyl_tower_fit.py [MAXW2=20]"""
import json
import os
import sys
import time
from fractions import Fraction as F

import sympy as sp

from super_engine import SuperSector, cur_add, d_cur, normal_prod, residue_cur, mono_weight2
from super_cyl_tower import tower, DX
import super_cyl_vev as scv

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, os.pardir, 'results', 'lab', 'susy4', 'tanabe')
l, P = sp.symbols('l P')
D2X, PSI, DPSI = (0, 2, 0), (1, 0, 0), (1, 1, 0)

MAXW2 = int(sys.argv[1]) if len(sys.argv) > 1 else 20
SAMPLES = [F(k, 13) for k in range(1, 61) if k % 13] + [F(k, 7) for k in (3, 5, 9, 11, 17)]
N_FIT = 40


def full_reduce(vec, rows, zero, div):
    """complete reduction of {index: coeff} modulo max-index pivot rows
    (integer rows); returns canonical remainder."""
    v = dict(vec)
    for k in sorted(set(v) | set(), reverse=True):
        pass
    keys = sorted(v, reverse=True)
    while keys:
        k = keys.pop(0)
        if k not in v or v[k] == zero:
            v.pop(k, None)
            continue
        if k in rows:
            row = rows[k]
            f = div(v[k], row[k])
            for kk, cc in row.items():
                nv = v.get(kk, zero) - f * cc
                v[kk] = nv
            v.pop(k, None)
            keys = sorted((x for x in v if x < k), reverse=True)
    return {k: c for k, c in v.items() if c != zero}


def canonical(cur, w2):
    sec = SuperSector(w2, 1)
    vec = {sec.index[m]: F(c) for m, c in cur.items() if c}
    red = full_reduce(vec, sec.piv.rows, F(0), lambda a, b: a / b)
    return {sec.basis[i]: c for i, c in red.items()}


def canonical_sym(cur, w2):
    sec = SuperSector(w2, 1)
    vec = {sec.index[m]: sp.sympify(c) for m, c in cur.items()}
    red = full_reduce(vec, sec.piv.rows, sp.Integer(0),
                      lambda a, b: sp.cancel(a / b))
    red = {k: sp.cancel(c) for k, c in red.items()}
    return {sec.basis[i]: c for i, c in red.items() if c != 0}


def rat_fit(points):
    """points: [(x, y)] Fractions; returns sympy rational function of l fitting
    the first N_FIT points with minimal total degree whose prediction matches
    all remaining points; raises if none up to degree 30."""
    fit, hold = points[:N_FIT], points[N_FIT:]
    for tot in range(0, 31):
        for dq in range(0, tot + 1):
            dp = tot - dq
            n = dp + 1 + dq   # q monic: q = l^dq + ...
            if n > len(fit):
                continue
            rows, rhs = [], []
            for x, y in fit:
                rows.append([x ** i for i in range(dp + 1)] + [-y * x ** j for j in range(dq)])
                rhs.append(y * x ** dq)
            M = sp.Matrix(rows)
            b = sp.Matrix(rhs)
            try:
                sol, params = M.gauss_jordan_solve(b)
            except ValueError:
                continue
            if params.shape[0]:
                sol = sol.subs({p_: 0 for p_ in params})
            pn = sum(sp.Rational(sol[i]) * l ** i for i in range(dp + 1))
            qd = l ** dq + sum(sp.Rational(sol[dp + 1 + j]) * l ** j for j in range(dq))
            f = sp.cancel(pn / qd)
            if all(f.subs(l, sp.Rational(x.numerator, x.denominator)) == sp.Rational(y.numerator, y.denominator)
                   for x, y in points):
                return f
    raise RuntimeError("no rational fit up to total degree 30")


def to_sym(cur):
    return {m: (sp.Rational(c.numerator, c.denominator) if isinstance(c, F) else c) for m, c in cur.items()}


def T_ls():
    return {(DX, DX): sp.Rational(1, 2), (D2X,): l, (PSI, DPSI): sp.Rational(-1, 2)}


def G_ls():
    return {(DX, PSI): sp.Integer(1), (DPSI,): 2 * l}


def clean(cur):
    return {m: sp.cancel(c) for m, c in cur.items() if sp.cancel(c) != 0}


def is_total_derivative_sym(cur):
    cur = clean(cur)
    if not cur:
        return True
    w2 = mono_weight2(next(iter(cur)))
    return not canonical_sym(cur, w2)


def main():
    t0 = time.time()
    out = {"maxw2": MAXW2, "convention": "P = zero mode of J_rec = p - l (Tanabe: Lambda + alpha_0/2, l = -alpha_0/2); "
           "densities normalised to (dX)^{2k} coefficient 1/2^k; I3 = :T_l T_l: - 1/4 :G_l dG_l:",
           "samples": [str(x) for x in SAMPLES], "n_fit": N_FIT, "charges": {}}
    # I3 symbolic
    T, G = T_ls(), G_ls()
    I3 = clean(normal_prod(T, T))
    for m, c in normal_prod(G, d_cur(G)).items():
        I3[m] = I3.get(m, 0) + sp.Rational(-1, 4) * c
    I3 = clean(I3)
    dens_sym = {8: canonical_sym(I3, 8)}
    # sampled tower
    samples = {}
    for x in SAMPLES:
        tw = tower(x, MAXW2)
        samples[x] = {w2: canonical(d, w2) for w2, d in tw.items() if w2 >= 12}
    print(f"sampled {len(SAMPLES)} values of l in {time.time()-t0:.1f}s", flush=True)
    for w2 in range(12, MAXW2 + 1, 4):
        monos = sorted(set().union(*(samples[x][w2].keys() for x in SAMPLES)))
        dens = {}
        for m in monos:
            pts = [(x, samples[x][w2].get(m, F(0))) for x in SAMPLES]
            dens[m] = rat_fit(pts)
        dens_sym[w2] = {m: c for m, c in dens.items() if c != 0}
        print(f"w2={w2}: {len(dens_sym[w2])} canonical coefficients reconstructed as rational functions of l "
              f"(fit on {N_FIT}, all {len(SAMPLES)} samples reproduced)", flush=True)
    # proof step: symbolic commutation
    for w2, D in dens_sym.items():
        t1 = time.time()
        okG = is_total_derivative_sym(residue_cur(G, D))
        ok3 = True if w2 == 8 else is_total_derivative_sym(residue_cur(dens_sym[8], D))
        print(f"w2={w2}: symbolic in l: [G_l, I] total derivative: {okG}; [I3(l), I] total derivative: {ok3} "
              f"({time.time()-t1:.1f}s)", flush=True)
        if not (okG and ok3):
            print("SYMBOLIC COMMUTATION FAILED", flush=True)
            sys.exit(1)
    # pairwise among higher charges
    ws = sorted(dens_sym)
    for i, a in enumerate(ws):
        for b in ws[i + 1:]:
            if a == 8:
                continue
            t1 = time.time()
            if a + b <= 32:
                ok = is_total_derivative_sym(residue_cur(dens_sym[a], dens_sym[b]))
                print(f"symbolic in l: [I(w2={a}), I(w2={b})] total derivative: {ok} ({time.time()-t1:.1f}s)", flush=True)
            else:
                # symbolic sympy bracket too slow (>30 min, killed in the first run):
                # exact check at held-out rational l on the reconstructed density
                ok = True
                for x in (F(5, 17), F(-7, 3), F(11, 9)):
                    xs = sp.Rational(x.numerator, x.denominator)
                    A = {m: F(str(sp.cancel(v.subs(l, xs)))) for m, v in dens_sym[a].items()}
                    B = {m: F(str(sp.cancel(v.subs(l, xs)))) for m, v in dens_sym[b].items()}
                    S = residue_cur(A, B)
                    okx = (not S) or SuperSector(mono_weight2(next(iter(S))), 1).is_total_derivative(S)
                    ok = ok and okx
                print(f"[I(w2={a}), I(w2={b})] total derivative at l = 5/17, -7/3, 11/9 (exact, reconstructed "
                      f"densities; not symbolic): {ok} ({time.time()-t1:.1f}s)", flush=True)
            if not ok:
                sys.exit(1)
    # eigenvalues
    c = sp.Rational(3, 2) - 12 * l ** 2
    Dl, cc = sp.symbols('Delta c')
    for w2, D in dens_sym.items():
        name = f"I{w2 // 2 - 1}"
        rec = {"w2": w2, "density_canonical": {repr(m): str(v) for m, v in sorted(D.items())}}
        for sec in ('NS', 'R'):
            t1 = time.time()
            e = sp.factor_terms(sp.expand(scv.scalar_vev(D, sec, [P])))
            e = sp.expand(e)
            # in Delta, c: P^2 = 2 Delta + l^2 - (1/8 if R), l^2 = (3/2 - c)/12
            shift = sp.Rational(1, 8) if sec == 'R' else 0
            e_dc = sp.expand(e.subs(P, sp.sqrt(2 * Dl + l ** 2 - shift)).subs(l, sp.sqrt((sp.Rational(3, 2) - cc) / 12)))
            rec[sec] = {"in_P_l": str(e), "in_Delta_c": str(sp.collect(e_dc, Dl)),
                        "in_sigma1_l": str(sp.expand(e.subs(P, sp.sqrt(sp.Symbol('sigma1')))))}
            print(f"{name} {sec} ({time.time()-t1:.1f}s):\n   in (P, l):      {e}\n   in (Delta, c):  {sp.collect(e_dc, Dl)}",
                  flush=True)
        out["charges"][name] = rec
    json.dump(out, open(os.path.join(OUT, 'computed_eigenvalues.json'), 'w'), indent=1, sort_keys=True)
    with open(os.path.join(OUT, 'computed_eigenvalues.txt'), 'w') as fh:
        fh.write("Computed cylinder eigenvalues (super_cyl_tower_fit.py), Fock highest weight, "
                 "P = p - l, Delta_R = Delta_NS(p) + 1/16; no printed values used.\n")
        for name, rec in out["charges"].items():
            for sec in ('NS', 'R'):
                fh.write(f"{name} {sec} (P,l):     {rec[sec]['in_P_l']}\n")
                fh.write(f"{name} {sec} (Delta,c): {rec[sec]['in_Delta_c']}\n")
    print(f"TOWER FIT COMPLETE ({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
