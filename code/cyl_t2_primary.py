"""M2, outcome form (results/cyl_t2_primary_preregistration.md): the
T-primariness split of the spin-9 kernel at t_2 = 7/5 (Solution 3).

Preregistered expectation (one primary class direction at t_2 -- the
Phi_{1,3}-shadow reading) was REFUTED in the lab run; this certificate
pins the observed counts:
  P1  kernel dims: 2 at t_2 = 7/5; 1 at the controls t = -7/5, 7/3, 11/5;
  P2  T_total-primary CLASS count 0 at every point (kernel classes
      modulo Im d; poles k = 3..12 of T(z)v(w) stacked; exact over Q);
  P3  under T_X alone and T_Y alone at t_2: 0 primary classes as well.
Uses code/mixed_pole.py (M1).  Fail-closed.  Marker: T2 PRIMARY SPLIT CERTIFIED."""
import sys
from fractions import Fraction as F

from invariant_engine import _rref
from mixed_engine import cyl_P4, d_mixed_mono, gen_mixed_basis, genuine_kernel_mixed
from mixed_pole import ope_pole_mixed_cur

fails = []


def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)


def point(t):
    t = F(t)
    return (t * t - 25) / (t * t - 1), F(-24) * t / (t * t - 1)


W = 10
BASIS = gen_mixed_basis(W)
IDX = {m: i for i, m in enumerate(BASIS)}
_dvecs = []
for b in gen_mixed_basis(W - 1):
    v = [F(0)] * len(BASIS)
    for m, c in d_mixed_mono(b).items():
        v[IDX[m]] += F(c)
    if any(v):
        _dvecs.append(v)
_rows, _piv = _rref([list(v) for v in _dvecs], len(BASIS))
DBASIS = [_rows[i] for i in range(len(_piv))]

T_TOT = {((1, 1), ()): F(1, 2), ((), ((1, 1),)): F(1, 2)}
T_X = {((1, 1), ()): F(1, 2)}
T_Y = {((), ((1, 1),)): F(1, 2)}


def pole_stack(cur, Nval, T):
    out = {}
    for k in range(3, 13):
        for m, c in ope_pole_mixed_cur(T, cur, k, Nval).items():
            out[(k, m)] = out.get((k, m), 0) + c
    return {key: c for key, c in out.items() if c}


def nullity(vecs, Nval, T, basis=None):
    basis = BASIS if basis is None else basis
    imgs, keys = [], set()
    for v in vecs:
        im = pole_stack({basis[i]: c for i, c in enumerate(v) if c}, Nval, T)
        imgs.append(im)
        keys |= set(im)
    keys = sorted(keys, key=str)
    if not keys:
        return len(vecs)
    _, pv = _rref([[im.get(kk, F(0)) for kk in keys] for im in imgs], len(keys))
    return len(vecs) - len(pv)


def primary_classes(Nval, s, T):
    ker = genuine_kernel_mixed(cyl_P4(3, Nval, s), BASIS, Nval)
    K = []
    for Kv in ker:
        v = [F(0)] * len(BASIS)
        for m, c in Kv.items():
            v[IDX[m]] += F(c)
        K.append(v)
    return len(K), nullity(K + DBASIS, Nval, T) - nullity(DBASIS, Nval, T)


def positive_control(Nval):
    """Codex round 17 (2026-09-02): every asserted count above is zero, so
    nothing showed the detector CAN return nonzero.  At weight 2 the
    invariant space is span{T_X, T_Y} (weight-1 invariants: none, so the
    derivative subspace is empty) and V = T_X - T_Y/N is T_total-primary:
    pole 4 of T_total(z)V(w) is c_X/2 - c_Y/(2N) = 0 with c_X = 1, c_Y = N,
    pole 3 vanishes, pole 2 is 2V.  T_total itself is not primary (pole 4
    = c/2).  Same pole_stack / nullity machinery, weight-2 basis."""
    basis2 = gen_mixed_basis(2)
    require(gen_mixed_basis(1) == [], "P4: weight-1 invariant basis must be empty")
    V = {((1, 1), ()): F(1, 2), ((), ((1, 1),)): F(-1, 2) / Nval}
    stack = pole_stack(V, Nval, T_TOT)
    pole2 = ope_pole_mixed_cur(T_TOT, V, 2, Nval)
    twoV = {m: 2 * c for m, c in V.items()}
    vec = lambda cur: [cur.get(m, F(0)) for m in basis2]
    pair = nullity([vec(T_X), vec(T_Y)], Nval, T_TOT, basis2)
    single = nullity([vec(T_TOT)], Nval, T_TOT, basis2)
    return stack, pole2 == twoV, pair, single


if __name__ == "__main__":
    EXPECT = {"7/5": 2, "-7/5": 1, "7/3": 1, "11/5": 1}
    for tt in EXPECT:
        N, s = point(tt)
        stack, pole2_ok, pair, single = positive_control(N)
        require(stack == {}, f"P4: t = {tt}: T_total(z)V(w) has no poles of order >= 3 (found: {stack})")
        require(pole2_ok, f"P4: t = {tt}: pole 2 of T_total(z)V(w) equals 2V")
        require(pair == 1, f"P4: t = {tt}: primary class count on span{{T_X, T_Y}} is {pair} "
                "(expected 1: the detector must see V)")
        require(single == 0, f"P4: t = {tt}: primary class count on span{{T_total}} is {single} "
                "(expected 0: T_total has pole 4 = c/2)")
    for tt, kd in EXPECT.items():
        N, s = point(tt)
        n, pc = primary_classes(N, s, T_TOT)
        require(n == kd, f"P1: t = {tt}: spin-9 kernel dim {n} (expected {kd})")
        require(pc == 0, f"P2: t = {tt}: T_total-primary class count {pc} (expected 0 -- the "
                "preregistered 'one primary shadow direction at t_2' is refuted)")
    N, s = point("7/5")
    for lab, T in (("T_X", T_X), ("T_Y", T_Y)):
        n, pc = primary_classes(N, s, T)
        require(pc == 0, f"P3: t_2 under {lab} alone: primary class count {pc} (expected 0)")
    print()
    if fails:
        print(f"T2 PRIMARY SPLIT: {len(fails)} PIN(S) NOT MET -- see FINDING lines")
        sys.exit(1)
    print("T2 PRIMARY SPLIT CERTIFIED (spin-9 kernel dims [2,1,1,1] at t_2 and three controls; "
          "ZERO T-primary class directions everywhere, under T_total, T_X and T_Y -- the extra "
          "class at t_2 is not a primary composite field; preregistered shadow reading refuted; "
          "positive control: the same detector returns 1 on span{T_X, T_Y} at weight 2, where "
          "V = T_X - T_Y/N is T_total-primary, and 0 on span{T_total})")
