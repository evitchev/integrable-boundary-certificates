"""BR2 lead's points (a), (b): direct checks independent of VIR2.
(a) LVZ (71): W4sym = (4n^2+9n+4)(d2X)^2 + (4n^2+7n+2)(d2Y)^2 + n(3n+4)(dX)^4 + (n+2)(3n+2)(dY)^4 + 6n(n+2)(dX)^2(dY)^2 (LVZ normalisation, chiral <XX> = -1/2 log)
    must be conserved mod d by the FOUR exponentials e^(+-sqrt(n) X +- i sqrt(n+2) Y) (LVZ (57), (58) and the mirror).  Map to the engine's real +log bosons:
    X = (i/sqrt2) Xt, Y = (i/sqrt2) Yt; at n = -2q^2 with q^2 + r^2 = 1 (rational) the four vectors are (+-q, +-r) in (Xt, Yt), all of norm 1.
    Test at n = -18/25 (q, r = 3/5, 4/5), -32/25 (4/5, 3/5), -50/169 (5/13, 12/13); tamper: coefficient 6n(n+2) -> + 1/7.
(b) N = 1 (p^2 = 2, rho = 0) in rotated coordinates U = (X + phi)/sqrt2, W = (phi - X)/sqrt2: the Sol 1 pair (+-1/sqrt2, -1/sqrt2) and its mirror
    (+-1/sqrt2, +1/sqrt2) become (-1,0), (0,-1) and (0,1), (1,0); the Virasoro screening e^(sqrt2 phi) becomes (1, 1).
    Compare the weight-4 and weight-6 commutants mod d of {four} and of {pair + Virasoro}: same 1-dim class?"""
import sys; sys.dont_write_bytecode = True
from fractions import Fraction as F
from eng import basis, residue, nullspace, deriv, pmul, padd
from br2_commutant import tw, rank, dim_mod_d
def conserved(P, screens, w):
    """is the density P (dict mono -> F) conserved mod d by each screening?  exists Y_s with residue(P, s) = tw(Y_s, s)"""
    Bl = basis(w - 2)
    for s in screens:
        R = {}
        for m, c in P.items():
            for om, cc in residue(m, s).items(): R[om] = R.get(om, 0) + c * cc
        R = {k: v for k, v in R.items() if v}
        cols = [tw({m: F(1)}, s) for m in Bl]
        keys = sorted({k for col in cols for k in col} | set(R))
        rows = [{j: col.get(k, F(0)) for j, col in enumerate(cols) if col.get(k, F(0))} for k in keys]
        # solve sum_j y_j col_j = R: augment
        A = [[col.get(k, F(0)) for col in cols] for k in keys]; b = [R.get(k, F(0)) for k in keys]
        from fractions import Fraction
        M = [row + [bb] for row, bb in zip(A, b)]; n = len(cols); r = 0
        for c in range(n):
            piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
            if piv is None: continue
            M[r], M[piv] = M[piv], M[r]; pv = M[r][c]; M[r] = [x / pv for x in M[r]]
            for i in range(len(M)):
                if i != r and M[i][c] != 0:
                    f = M[i][c]; M[i] = [x - f * y for x, y in zip(M[i], M[r])]
            r += 1
        if any(M[i][n] != 0 for i in range(r, len(M))): return False
    return True
bad = 0
for n_, q, r in ((F(-18, 25), F(3, 5), F(4, 5)), (F(-32, 25), F(4, 5), F(3, 5)), (F(-50, 169), F(5, 13), F(12, 13))):
    assert n_ == -2 * q * q and q * q + r * r == 1
    # (dX)^2 = -1/2 (dXt)^2 etc.; (d2X)^2 = -1/2 (d2Xt)^2
    def W4(tam=0):
        return {((0, 2), (0, 2)): -F(1, 2) * (4 * n_**2 + 9 * n_ + 4), ((1, 2), (1, 2)): -F(1, 2) * (4 * n_**2 + 7 * n_ + 2),
                ((0, 1),) * 4: F(1, 4) * n_ * (3 * n_ + 4), ((1, 1),) * 4: F(1, 4) * (n_ + 2) * (3 * n_ + 2),
                ((0, 1), (0, 1), (1, 1), (1, 1)): F(1, 4) * (6 * n_ * (n_ + 2) + tam)}
    four = [(q, r), (q, -r), (-q, r), (-q, -r)]
    ok = conserved(W4(), four, 4); tamp = conserved(W4(F(1, 7)), four, 4)
    print(f'(a) LVZ W4sym at n = {n_}: conserved mod d by all four hairpin exponentials: {ok}; TAMPER (+1/7): {tamp} (must be False)', flush=True)
    bad += (not ok) or tamp
# (b) N = 1 bookkeeping
pair = [(F(-1), F(0)), (F(0), F(-1))]; mirror = [(F(0), F(1)), (F(1), F(0))]; vir = [(F(1), F(1))]
for w in (4, 6):
    d4 = dim_mod_d(w, pair + mirror); dpv = dim_mod_d(w, pair + vir); dall = dim_mod_d(w, pair + mirror + vir); dp = dim_mod_d(w, pair)
    print(f'(b) N = 1, weight {w}: commutant mod d of pair: {dp}; pair + mirror (four): {d4}; pair + Virasoro: {dpv}; all five: {dall}', flush=True)
print('DIRECT CHECKS', 'PASS' if not bad else 'FAIL')
sys.exit(1 if bad else 0)
