"""Utilities: kernel of ad_{I_low} on charges, change of representatives,
rational-function reconstruction from samples."""

from fractions import Fraction as F

from ope_engine import Sector, residue_cur, nullspace, rref


def genuine_kernel(P_low, cand_basis, target_sec, cand_sec):
    """Charges oint J, J in span(cand_basis), with [oint P_low, oint J] = 0,
    modulo total derivatives.  Returns list of currents (dicts), one per
    independent genuine commuting charge."""
    cols = []
    for b in cand_basis:
        r = residue_cur(P_low, {b: F(1)})
        cols.append(target_sec.reduce_mod_d(target_sec.vec(r)))
    ker = nullspace(cols, len(target_sec.basis))
    reduced, genuine = [], []
    for v in ker:
        cur = {b: c for b, c in zip(cand_basis, v) if c}
        w = cand_sec.reduce_mod_d(cand_sec.vec(cur))
        if not any(w):
            continue
        rows = [list(x) for x in reduced + [w]]
        _, piv = rref(rows, len(w))
        if len(piv) == len(reduced) + 1:
            reduced.append(w)
            genuine.append(cur)
    return genuine


def express_in(cur, representatives, sector):
    """Coefficients of `cur` in the given representative monomials,
    modulo total derivatives.  Raises if not in their span."""
    cols = [sector.reduce_mod_d(sector.vec({b: F(1)}))
            for b in representatives]
    target = sector.reduce_mod_d(sector.vec(cur))
    nrows = len(target)
    aug = [[cols[c][r] for c in range(len(cols))] + [target[r]]
           for r in range(nrows)]
    aug, piv = rref(aug, len(cols) + 1)
    if len(cols) in piv:
        raise ValueError("current not in span of representatives")
    sol = [F(0)] * len(cols)
    for i, p in enumerate(piv):
        sol[p] = aug[i][len(cols)]
    return sol


def rat_reconstruct(samples, dp, dq):
    """Find p/q, deg p <= dp, deg q <= dq, with p(n_k)/q(n_k) = f_k for all
    (n_k, f_k) in samples (Fractions).  Returns (pcoeffs, qcoeffs) lowest
    degree first, or None.  Caller should verify on held-out samples."""
    rows = []
    for nk, fk in samples:
        row = [nk**i for i in range(dp + 1)]
        row += [-fk * nk**i for i in range(dq + 1)]
        rows.append(row)
    ncols = dp + 1 + dq + 1
    cols = [[rows[r][c] for r in range(len(rows))] for c in range(ncols)]
    ker = nullspace(cols, len(rows))
    if not ker:
        return None
    v = ker[0]
    p, q = v[:dp + 1], v[dp + 1:]
    if not any(q):
        return None
    return p, q


def rat_eval(pq, n):
    p, q = pq
    num = sum(c * n**i for i, c in enumerate(p))
    den = sum(c * n**i for i, c in enumerate(q))
    if den == 0:
        return None
    return num / den
