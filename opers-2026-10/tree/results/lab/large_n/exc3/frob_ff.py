"""EXC3: fraction-free variant of frob.nolog_conditions (same mathematics; the Frobenius coefficients are kept as
numerator polynomials over a running common denominator D_j = prod_(i<=j, non-resonant) lead(e1+i), a polynomial in z only).
Returns (lead, conds, kaps) with conds = [(exponent, numerator expression)] -- the numerators are what split_conditions needs."""
import sys
import sympy as sp
import frob as Fb
sys.dont_write_bytecode = True
e = Fb.e


def nolog_conditions_ff(lk, exps, order):
    e1 = exps[0]; span = exps[-1] - e1
    lead = lk.get(-order, 0)
    N = {0: sp.Integer(1)}          # numerators: c_j = N_j / D_j
    D = {0: sp.Integer(1)}
    kap = {}; conds = []
    for j in range(1, span + 1):
        # rhs = sum_k l_k(e1 + j - k) c_(j-k);  common denominator D_(j-1) (all D_i, i < j, divide it)
        rhs = sp.Integer(0)
        for kk in range(1, j + 1):
            lcoef = lk.get(kk - order, 0)
            if lcoef == 0 or (j - kk) not in N:
                continue
            rhs += lcoef.subs(e, e1 + j - kk) * N[j - kk] * sp.cancel(D[j - 1] / D[j - kk])
        rhs = sp.expand(rhs)
        if (e1 + j) in exps:
            conds.append((e1 + j, rhs))
            kap[j] = sp.Symbol('kap%d' % j)
            N[j] = kap[j] * D[j - 1]; D[j] = D[j - 1]        # c_j = kap_j
        else:
            den = sp.expand(lead.subs(e, e1 + j))
            N[j] = sp.expand(-rhs); D[j] = sp.expand(D[j - 1] * den)
    return lead, conds, list(kap.values())
