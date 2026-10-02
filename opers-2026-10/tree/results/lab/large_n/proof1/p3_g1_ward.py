"""PROOF1 (3): G1 (Fable, item 306) == the record's Ward-derived loss-1 layer (Codex, items 180/185), IDENTICALLY in (a, K, t), all three solutions.
Object compared: R(a, K, t) = [X^a Y^(K-1-a)] / [X^a Y^(K-a)] of the record-normalised charge (Codex: 'projection_over_c_a', an explicit rational function).
G1 side, exactly: every term of G1 is [z^N] of prod_e (1 - e z)^(nu alpha_e - d_e) times a monomial.  Paired families (1 -+ l z) are rewritten as
(1 - l^2 z^2)^(beta - D) times a FINITE polynomial; a single family (Sol 1's l0) is one binomial.  So the coefficient of l0^(2i) l1^(2b) is a finite sum of
binomials whose upper/lower arguments differ from the top's by INTEGERS, and the ratio to the top is a rational function of (i, K, t) (helper bratio).
Record conversion (Fable's to_record): l0^2 = -CX X, l1^2 = CY (rho^2 - Y), rho^2 = 2/(t^2 - 1), giving R = -c_L1/(CY c_top) - (K - i) rho^2.
Then sympy.cancel(R_G1 - R_Codex) == 0 is a proof of the identity in Q(i, K, t)."""
import sys, json; sys.dont_write_bytecode = True
import sympy as sp
from itertools import product
i, K, t = sp.symbols('i K t')
def rfi(y, n):        # rising factorial with integer n of either sign
    n = int(n)
    if n >= 0: return sp.Mul(*[y + j for j in range(n)])
    return 1 / sp.Mul(*[y - 1 - j for j in range(-n)])
def bratio(x, dx, m, dm):   # binom(x + dx, m + dm) / binom(x, m), dx, dm integers
    return rfi(x + 1, dx) / rfi(m + 1, dm) / rfi(x - m + 1, dx - dm)
def solution_data(sol):
    rho2 = 2 / (t**2 - 1)
    if sol == 1:
        a_s = (t - 1) / 2; n = (t + 3) / 2; M = 1 / a_s
        CY = n**2 * (M + 1) / 2                    # l1 = n pi / p
        elems = [('l0', 1, a_s), ('l1', 1, 1), ('l1', -1, 1)]; paired = {'l0': False, 'l1': True}
        strings = [('l0', 1, a_s)]
    else:
        k = 2 * (t - 1) / (3 - t); M = 1 / k
        if sol == 3:
            n = k + 4; aX = (t - 3) / (t - 5)
            elems = [('l0', 1, 1), ('l0', -1, 1), ('l1', 1, 1), ('l1', -1, 1)]; strings = [(None, 0, k)]
        else:
            n = 2 * k + 3; aX = -2 * (t - 3) / (t + 5)
            elems = [('l0', 1, k), ('l0', -1, k), ('l1', 1, 1), ('l1', -1, 1)]; strings = [('l0', 1, k), ('l0', -1, k)]
        CY = 2 * n * (1 + M) / aX
        paired = {'l0': True, 'l1': True}
    return dict(n=sp.simplify(n), M=sp.simplify(M), CY=sp.simplify(CY), rho2=rho2, elems=elems, strings=strings, paired=paired)
def terms(D):
    """G1 terms: list of (coef, order N-label ('top' or 'l1'), shifts {(fam, sgn): d}, ell-powers {fam: u})."""
    n, M = D['n'], D['M']; nu = (2 * K - 1) / n; c2 = n**2 * (M + 1); sig = n * M
    out = [(sp.Integer(1), 'top', {}, {})]
    pre = nu * c2 / 24
    out.append((pre * (n - 1), 'l1', {}, {}))
    L1 = [(-al * sg, {(f, sg): 1}, {f: 1}) for f, sg, al in D['elems']]            # -alpha e z/(1 - e z), e = sg * ell_f
    L2 = [(-al * sg, {(f, sg): 2}, {f: 1}) for f, sg, al in D['elems']]
    for c, sh, u in L1: out.append((-pre * nu * n * c, 'l1', sh, u))
    for c, sh, u in L2: out.append((pre * nu * c, 'l1', sh, u))
    for (c1_, s1, u1), (c2_, s2, u2) in product(L1, L1):
        sh = dict(s1)
        for kk, v in s2.items(): sh[kk] = sh.get(kk, 0) + v
        uu = dict(u1)
        for kk, v in u2.items(): uu[kk] = uu.get(kk, 0) + v
        out.append((pre * nu * c1_ * c2_, 'l1', sh, uu))
    for f, sg, a in D['strings']:
        c1 = -a * (a**2 - 1) / 24
        sh = {} if f is None else {(f, sg): 2}
        out.append((nu * c1 * sig**2, 'l1', sh, {}))
    return out
def family_coeff(D, fam, shifts, u, P, Ptop, nu):
    """coefficient of ell^P in the family's factor (with prefactor ell^u), divided by the top coefficient of ell^Ptop. Symbolic P, Ptop (in i, K)."""
    al = [a for f, s, a in D['elems'] if f == fam][0]
    beta = nu * al
    if not D['paired'][fam]:
        d = shifts.get((fam, 1), 0)
        # (1 - ell z)^(beta - d): coefficient of ell^(P - u) = binom(beta - d, P - u)(-1)^(P - u); top: binom(beta, Ptop)(-1)^Ptop, Ptop even
        dP = sp.simplify((P - u) - Ptop); assert dP.is_integer
        return bratio(beta, -d, Ptop, int(dP)) * (1 if int(dP) % 2 == 0 else -1)
    dp, dm = shifts.get((fam, 1), 0), shifts.get((fam, -1), 0); Dm = max(dp, dm)
    poly = sp.Poly(sp.expand((1 - sp.Symbol('w'))**(Dm - dp) * (1 + sp.Symbol('w'))**(Dm - dm)), sp.Symbol('w'))
    tot = 0
    for (m,), c in poly.terms():
        rem = sp.simplify(P - u - m)
        if not (rem / 2).is_integer and not sp.simplify(rem - 2 * sp.floor(rem / 2)) == 0:
            # parity: rem = (even symbolic) + integer offset; odd offset -> odd ell power -> no contribution to the even target
            off = sp.simplify(rem - (P - Ptop) - Ptop)
        jj = sp.simplify(rem / 2); jtop = sp.simplify(Ptop / 2); dj = sp.simplify(jj - jtop)
        if not dj.is_integer: continue      # odd ell power: does not reach the even target
        # (1 - ell^2 z^2)^(beta - Dm): coeff binom(beta - Dm, jj)(-1)^jj; top: binom(beta, jtop)(-1)^jtop
        tot += c * bratio(beta, -Dm, jtop, int(dj)) * (1 if int(dj) % 2 == 0 else -1)
    return tot
def g1_ratio(sol):
    D = solution_data(sol); n = D['n']; nu = (2 * K - 1) / n
    Ptop = {'l0': 2 * i, 'l1': 2 * (K - i)}
    PL1 = {'l0': 2 * i, 'l1': 2 * (K - 1 - i)}
    cL1 = 0
    for coef, kind, sh, u in terms(D):
        if kind == 'top': continue
        val = coef
        for fam in ('l0', 'l1'):
            val *= family_coeff(D, fam, sh, u.get(fam, 0), PL1[fam], Ptop[fam], nu)
        cL1 += val
    # c_L1 here is already divided by c_top (family_coeff returns ratios to the top's binomials)
    return sp.cancel(sp.together(-cL1 / D['CY'] - (K - i) * D['rho2']))
def codex_ratio(sol):
    d = json.load(open(f'codex_sol{sol}_projection_loss1.json'))['projection_over_c_a']
    loc = {'a': i, 'k': K, 't': t}
    return sp.sympify(d['numerator'], locals=loc, rational=True) / sp.sympify(d['denominator'], locals=loc, rational=True)

def main():
  ok = True
  for sol in (1, 2, 3):
    G = g1_ratio(sol); C = codex_ratio(sol)
    diff = sp.cancel(sp.together(G - C))
    print(f'Sol {sol}: G1 ratio - Codex ratio = {diff}   -> IDENTITY IN (a, K, t): {diff == 0}', flush=True)
    ok &= diff == 0
    if diff != 0: print('   G1:', sp.factor(G)); print('   Codex:', sp.factor(C))
  sys.exit(0 if ok else 1)
if __name__ == '__main__': main()
