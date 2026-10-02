"""EXC2 addendum, numerical route (mpmath, 40 digits): residuals of the no-logarithm conditions for N apparent singularities
of type {-1,1,2,4} of the Sol-3 three-term operator, and Newton iteration.  Arithmetic written out directly (no lambdify).
Unknowns per point: (r21, r11, r12, r01, r02, r03, z); r22 = -4z, r13 = 8z^2, r04 = -8z^3."""
import sys
from mpmath import mp, mpf, mpc, binomial, matrix, findroot, norm
sys.dont_write_bytecode = True
mp.dps = 40
EXPS = [-1, 1, 2, 4]


def padd(A, B, cb=1):
    out = dict(A)
    for k, v in B.items():
        out[k] = out.get(k, 0) + cb * v
    return out


def pmul(A, B):
    out = {}
    for k1, v1 in A.items():
        for k2, v2 in B.items():
            k = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2])
            if sum(k) > 2:
                continue
            out[k] = out.get(k, 0) + v1 * v2
    return out


def residuals(points, a0, a1, b):
    """points: list of 7-tuples.  Returns list of 7 residuals per point."""
    res = []
    P4 = [a0 * a0 * a1 * a1, 0, -(a0 * a0 + a1 * a1), 0, 1]
    for idx, pt in enumerate(points):
        r21, r11, r12, r01, r02, r03, z = pt
        own = {(2, 1): r21, (2, 2): -4 * z, (1, 1): r11, (1, 2): r12, (1, 3): 8 * z * z, (0, 1): r01, (0, 2): r02, (0, 3): r03, (0, 4): -8 * z ** 3}
        others = []
        for jdx, q in enumerate(points):
            if jdx == idx:
                continue
            s21, s11, s12, s01, s02, s03, zo = q
            others.append(({(2, 1): s21, (2, 2): -4 * zo, (1, 1): s11, (1, 2): s12, (1, 3): 8 * zo * zo, (0, 1): s01, (0, 2): s02, (0, 3): s03, (0, 4): -8 * zo ** 3}, zo))
        K = 1

        def lk_at(ev):
            """dict k -> (l_k^(0)(ev), E-coefficient) for k = -4..K"""
            powers = [{0: mpf(1)}]
            for j in range(4):
                nxt = {}
                for kk, g in powers[-1].items():
                    f = g * (ev + kk)
                    nxt[kk - 1] = nxt.get(kk - 1, 0) + z * f
                    nxt[kk] = nxt.get(kk, 0) + f
                powers.append({kk: v for kk, v in nxt.items() if kk <= K})
            out = {}
            def addser(ser, coef=1):
                for kk, v in ser.items():
                    if kk <= K:
                        out[kk] = out.get(kk, 0) + coef * v
            def times_xi(ser):
                o = {}
                for kk, v in ser.items():
                    o[kk] = o.get(kk, 0) + z * v
                    o[kk + 1] = o.get(kk + 1, 0) + v
                return o
            def shift(ser, coeffs):
                o = {}
                for kk, v in ser.items():
                    for p_, c in coeffs.items():
                        if kk + p_ <= K:
                            o[kk + p_] = o.get(kk + p_, 0) + v * c
                return o
            for j, c in enumerate(P4):
                if c != 0:
                    addser(powers[j], c)
            addser(times_xi(powers[0]), mpf(1) / 2)
            addser(times_xi(powers[1]), 1)
            for (j, m), c in own.items():
                addser(shift(times_xi(powers[j]), {-m: c}))
            for oth, zo in others:
                Dd = z - zo
                for (j, m), c in oth.items():
                    ser = {i: c * binomial(-m, i) / Dd ** (m + i) for i in range(0, K + 5)}
                    addser(shift(times_xi(powers[j]), ser))
            Epart = {i: -binomial(b, i) / z ** i for i in range(0, K + 1)}      # E-term acts multiplicatively: starts at k = 0
            return out, Epart

        e1 = EXPS[0]
        span = EXPS[-1] - e1
        cache = {}
        c = {0: {(0, 0, 0): mpf(1)}}
        kapidx = {EXPS[1] - e1: 1, EXPS[2] - e1: 2}
        for j in range(1, span + 1):
            rhs = {}
            for kk in range(1, j + 1):
                if (j - kk) not in c:
                    continue
                ev = e1 + j - kk
                if ev not in cache:
                    cache[ev] = lk_at(ev)
                l0, Ep = cache[ev]
                k = kk - 4
                term = {}
                if k in l0 and l0[k] != 0:
                    term = padd(term, {(0, 0, 0): l0[k]})
                if k in Ep:
                    term = padd(term, {(1, 0, 0): Ep[k]})
                if term:
                    rhs = padd(rhs, pmul(term, c[j - kk]))
            if (e1 + j) in EXPS:
                pos = e1 + j
                if pos == EXPS[1]:
                    res.append(rhs.get((0, 0, 0), 0))
                    c[j] = {(0, 1, 0): mpf(1)}
                elif pos == EXPS[2]:
                    res.append(rhs.get((0, 1, 0), 0)); res.append(rhs.get((0, 0, 0), 0))
                    c[j] = {(0, 0, 1): mpf(1)}
                else:
                    res.append(rhs.get((1, 0, 0), 0)); res.append(rhs.get((0, 1, 0), 0)); res.append(rhs.get((0, 0, 1), 0)); res.append(rhs.get((0, 0, 0), 0))
            else:
                ev = e1 + j
                if ev not in cache:
                    cache[ev] = lk_at(ev)
                lead = cache[ev][0][-4]
                c[j] = {kq: -v / lead for kq, v in rhs.items()}
    return res


def newton(points0, a0, a1, b, tol=mpf(10) ** (-30), maxit=60):
    npts = len(points0)
    x0 = [v for pt in points0 for v in pt]
    def F(*x):
        pts = [tuple(x[7 * i:7 * i + 7]) for i in range(npts)]
        return residuals(pts, a0, a1, b)
    try:
        sol = findroot(F, x0, tol=tol, maxsteps=maxit)
    except Exception as ex:
        return None
    x = [sol[i] for i in range(7 * npts)]
    if norm(matrix(F(*x))) > mpf(10) ** (-25):
        return None
    return [tuple(x[7 * i:7 * i + 7]) for i in range(npts)]
