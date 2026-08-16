"""Q(t)-promotion of sampled existence results.

For the A2^(2) charges I_11, I_13 and the cylindrical charges
I_7, I_9, I_11 (all three solutions): solve the ad-kernel at many
rational points +-t of the double cover s^2 = (N-25)(N-1)
(N = (t^2-25)/(t^2-1), s = -24t/(t^2-1); deck map t -> -t), split each
canonical density coordinate as A(N) + s B(N) using the +-t pairs,
reconstruct A, B as rational functions of N with held-out verification,
then VERIFY SYMBOLICALLY over Q(t) that the reconstructed density
commutes with the low charge --- proving existence at the generic point
of the parameter curve.  Combined with sampled kernel dimension 1
(semicontinuity: generic dim <= sampled dim), this proves the generic
kernel is exactly one-dimensional at these spins.

Output: results/qt_promotion.json (exact A, B strings) and a summary.
Run:  python312 promote_qt.py   (uses 11 parallel workers)
"""

import json
from fractions import Fraction as F
from multiprocessing import Pool

import sympy as sp

from invariant_engine import (InvSector, gen_inv_basis, genuine_kernel_inv,
                              residue_inv_mono, a2_P6, curve_point,
                              spin_inv)
from mixed_engine import (MixedSector, gen_mixed_basis,
                          genuine_kernel_mixed, residue_mixed_mono,
                          cyl_P4, spin_mixed)
from solve_charges import rat_reconstruct, rat_eval

TS = [F(x) for x in (2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16,
                     18, 21, 23,
                     F(5, 2), F(7, 2), F(9, 2), F(11, 2), F(13, 2),
                     F(15, 2), F(21, 2),
                     F(7, 3), F(8, 3), F(10, 3), F(11, 3), F(13, 3),
                     F(14, 3), F(16, 3),
                     F(17, 2), F(19, 2))]

t_sym = sp.Symbol('t')
N_OF_T = (t_sym**2 - 25) / (t_sym**2 - 1)
S_OF_T = -24 * t_sym / (t_sym**2 - 1)


def fit(points, label):
    """Reconstruct a rational function of N from (N, value) samples."""
    for d in [(2, 2), (4, 4), (6, 6), (8, 8), (10, 10), (12, 12)]:
        need = d[0] + d[1] + 2
        if need + 3 > len(points):
            break
        train, held = points[:need + 1], points[need + 1:]
        pq = rat_reconstruct(train, *d)
        if pq is None:
            continue
        if all(rat_eval(pq, x) == y for x, y in held):
            return pq
    raise RuntimeError(f"reconstruction failed: {label}")


def pq_to_expr(pq, var):
    p, q = pq
    num = sum(sp.Rational(c) * var**k for k, c in enumerate(p))
    den = sum(sp.Rational(c) * var**k for k, c in enumerate(q))
    return sp.cancel(num / den)


def task(args):
    try:
        return _task(args)
    except Exception as e:
        fam, sol, sigma = args
        return f"{fam}{sol or ''}_sigma{sigma}", False, 0, {"error": repr(e)}, -1


def _task(args):
    fam, sol, sigma = args
    if fam == "a2":
        low = lambda t: a2_P6(*curve_point(t))
        basis = gen_inv_basis(sigma + 1)
        Sec, res_mono = InvSector, residue_inv_mono
        kern = lambda P, t: genuine_kernel_inv(P, basis, curve_point(t)[0])
        s_low = 6
    else:
        low = lambda t: cyl_P4(sol, *curve_point(t))
        basis = gen_mixed_basis(sigma + 1)
        Sec, res_mono = MixedSector, residue_mixed_mono
        kern = lambda P, t: genuine_kernel_mixed(P, basis, curve_point(t)[0])
        s_low = 4
    sec = Sec(sigma + 1)
    tgt = Sec(s_low + sigma)

    # -------- sampling at +-t pairs
    vecs = {}
    for t0 in TS:
        for t in (t0, -t0):
            g = kern(low(t), t)
            if len(g) != 1:
                break
            vecs[t] = sec.reduce_mod_d(sec.vec(g[0]))
        else:
            continue
        vecs.pop(t0, None)
    pairs = [t0 for t0 in TS if t0 in vecs and -t0 in vecs]
    allts = [t for t0 in pairs for t in (t0, -t0)]
    nz = [i for i in range(len(basis))
          if all(vecs[t][i] for t in allts)]
    if not nz:  # fall back: best coordinate, drop offending pairs
        cnt = {i: sum(1 for t in allts if vecs[t][i])
               for i in range(len(basis))}
        j = max(cnt, key=cnt.get)
        pairs = [t0 for t0 in pairs
                 if vecs[t0][j] and vecs[-t0][j]]
        allts = [t for t0 in pairs for t in (t0, -t0)]
    else:
        j = nz[0]
    for t in allts:
        v = vecs[t]
        vecs[t] = [x / v[j] for x in v]

    # -------- A(N) + s B(N) split and reconstruction
    live = [i for i in range(len(basis))
            if any(vecs[t][i] for t in allts)]
    A, B = {}, {}
    for i in live:
        ptsA, ptsB = [], []
        for t0 in pairs:
            Nv, sv = curve_point(t0)
            vp, vm = vecs[t0][i], vecs[-t0][i]
            ptsA.append((Nv, (vp + vm) / 2))
            ptsB.append((Nv, (vp - vm) / (2 * sv)))
        A[i] = fit(ptsA, f"{fam}{sol or ''} s{sigma} A[{i}]")
        B[i] = fit(ptsB, f"{fam}{sol or ''} s{sigma} B[{i}]")

    # -------- symbolic density over Q(t) and symbolic verification
    coeff = {}
    for i in live:
        e = (pq_to_expr(A[i], sp.Symbol('N'))
             + sp.Symbol('s') * pq_to_expr(B[i], sp.Symbol('N')))
        e = e.subs({sp.Symbol('N'): N_OF_T, sp.Symbol('s'): S_OF_T})
        coeff[i] = sp.cancel(e)
    den = sp.lcm([sp.denom(c) for c in coeff.values()])
    dens = {basis[i]: sp.expand(sp.cancel(c * den))
            for i, c in coeff.items()}
    # low charge symbolically: coefficients rational in t
    if fam == "a2":
        Nsym, ssym = N_OF_T, S_OF_T
        lows = {m: sp.cancel(e) for m, e in _a2_sym(Nsym, ssym).items()}
    else:
        lows = {m: sp.cancel(e) for m, e in _cyl_sym(sol, N_OF_T, S_OF_T).items()}
    dlow = sp.lcm([sp.denom(c) for c in lows.values()])
    lows = {m: sp.expand(sp.cancel(c * dlow)) for m, c in lows.items()}

    out = {}
    Nt = N_OF_T
    for mp, cp in lows.items():
        for mq, cq in dens.items():
            for mono, powd in res_mono(mp, mq).items():
                val = cp * cq * sum(sp.Rational(c.numerator, c.denominator)
                                    * Nt**p for p, c in powd.items())
                out[mono] = out.get(mono, 0) + val
    v = [sp.Integer(0)] * len(tgt.basis)
    for mono, c in out.items():
        v[tgt.index[mono]] = sp.cancel(c)
    for row, p in zip(tgt._img_rows, tgt._img_pivots):
        f0 = v[p]
        if f0 == 0:
            continue
        v = [sp.cancel(a - f0 * sp.Rational(b.numerator, b.denominator))
             for a, b in zip(v, row)]
    ok = all(sp.cancel(x) == 0 for x in v)
    key = f"{fam}{sol or ''}_sigma{sigma}"
    return key, ok, len(pairs), {str(i): [str(pq_to_expr(A[i], sp.Symbol('N'))),
                                          str(pq_to_expr(B[i], sp.Symbol('N')))]
                                 for i in live}, j


def _a2_sym(N, s):
    from invariant_engine import canon, canon_pair
    def M(*pairs):
        return canon(tuple(canon_pair(*p) for p in pairs))
    return {M((3, 3)): (107 + 17 * N + 15 * s) / 96,
            M((2, 2), (1, 1)): sp.Integer(-6),
            M((1, 2), (1, 2)): -(-85 + 17 * N + 15 * s) / 8,
            M((1, 1), (1, 1), (1, 1)): sp.Integer(1)}


def _cyl_sym(sol, N, s):
    from mixed_engine import _M
    base = {_M((), [(1, 1), (1, 1)]): sp.Integer(1),
            _M((), [(2, 2)]): sp.Integer(-2)}
    if sol == 1:
        base.update({_M((2, 2)): (-9 - 3 * N - s) / 6,
                     _M((1, 1), [(1, 1)]): sp.Integer(6),
                     _M((1, 1, 1, 1)): (2 + N + s) / 3})
    elif sol == 2:
        base.update({_M((2, 2)): -(409 + 46 * N + N**2 + (N - 5) * s) / 192,
                     _M((1, 1), [(1, 1)]): (11 + N + s) / 8,
                     _M((1, 1, 1, 1)): (117 + 2 * N + N**2 + (15 + N) * s) / 192})
    else:
        base.update({_M((2, 2)): (13 - 27 * N + 2 * N**2 + (5 - 2 * N) * s) / 6,
                     _M((1, 1), [(1, 1)]): (7 - N + s),
                     _M((1, 1, 1, 1)): sp.Integer(1)})
    return base


if __name__ == "__main__":
    tasks = ([("a2", None, 11), ("a2", None, 13)] +
             [("cyl", sol, sig) for sol in (1, 2, 3) for sig in (7, 9, 11)])
    with Pool(len(tasks)) as pool:
        results = pool.map(task, tasks)
    summary = {}
    for key, ok, npairs, data, j in results:
        print(f"{key}: symbolic Q(t) verification "
              f"{'PASSED' if ok else 'FAILED'}  ({npairs} sample pairs)")
        summary[key] = {"symbolic_verified": ok, "sample_pairs": npairs,
                        "normalized_coordinate": j, "coeffs_A_B_of_N": data}
    ok_all = all(r[1] for r in results)
    dest = "../results/qt_promotion.json" if ok_all \
        else "../results/qt_promotion.partial.json"
    with open(dest, "w") as f:
        json.dump(summary, f, indent=1)
    print(f"written to {dest}")
    if ok_all:
        print("\nALL EXISTENCE RESULTS PROMOTED TO Q(t): generic existence "
              "proven;\ncombined with sampled dim <= 1 (semicontinuity), "
              "generic kernels are exactly one-dimensional.")
    else:
        import sys
        print("\nFAILURES present; successful artifact NOT overwritten.")
        sys.exit(1)
