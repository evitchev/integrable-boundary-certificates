"""VIR2: the screening pairs verified SYMBOLICALLY in p = beta_+ (so u^2 = p^2/2,
t = (2 + p^2)/(2 - p^2), rho = p/2 - 1/p, N = 1 - 12 rho^2, s = 3 p^2 - 12/p^2),
against I_3 and I_5, all three solutions.  Exit 0 iff every listed pair is conserved
identically in p and every negative control fails."""
import sys, time, pickle
import sympy as sp
sys.dont_write_bytecode = True
import vir_lib as V
import vir2_lib as W
from vir2_lib import a, rho

p, t = sp.symbols('p t')
tp = (2 + p**2) / (2 - p**2)
rp = p / 2 - 1 / p
Np = sp.cancel(1 - 12 * rp**2)
sp_ = sp.cancel(3 * p**2 - 12 / p**2)
assert sp.cancel(Np - (tp**2 - 25) / (tp**2 - 1)) == 0 and sp.cancel(sp_ + 24 * tp / (tp**2 - 1)) == 0
print('curve in p verified: N = (t^2-25)/(t^2-1), s = -24t/(t^2-1) with t = (2+p^2)/(2-p^2)')
I5 = {k: [sp.sympify(x) for x in v] for k, v in pickle.load(open('i5_sym.pkl', 'rb')).items()}
B4 = V.basis(4)
i4 = {m: i for i, m in enumerate(B4)}
t0 = time.time()


def P4p(sol):
    N, s = Np, sp_
    base = {((), ((1, 1), (1, 1))): sp.Integer(1), ((), ((2, 2),)): sp.Integer(-2)}
    if sol == 1:
        ex = {((2, 2), ()): (-9 - 3*N - s)/6, ((1, 1), ((1, 1),)): sp.Integer(6), ((1, 1, 1, 1), ()): (2 + N + s)/3}
    elif sol == 2:
        ex = {((2, 2), ()): -(409 + 46*N + N*N + (N - 5)*s)/192, ((1, 1), ((1, 1),)): (11 + N + s)/8,
              ((1, 1, 1, 1), ()): (117 + 2*N + N*N + (15 + N)*s)/192}
    else:
        ex = {((2, 2), ()): (13 - 27*N + 2*N*N + (5 - 2*N)*s)/6, ((1, 1), ((1, 1),)): 7 - N + s, ((1, 1, 1, 1), ()): sp.Integer(1)}
    base.update(ex)
    vec = [sp.Integer(0)] * len(B4)
    for m, c in base.items():
        vec[i4[m]] = sp.cancel(c)
    return vec


def densities(sol):
    d3, _ = W.to_two_boson(P4p(sol), 4)
    v5 = [sp.cancel(c.subs(t, tp) * (2 - p**2)**5) for c in I5[sol]]
    d5, _ = W.to_two_boson(v5, 6)
    sub = lambda d: {m: sp.cancel(c.subs(rho, rp)) for m, c in d.items()}
    return sub(d3), sub(d5)


def conserved(density, bexpr, qexpr):
    """All normal-form coefficients vanish for V = exp(a X + bexpr phi) with a^2 = qexpr."""
    if qexpr == 0:
        nf = W.normal_form(W.screening_residue(density, sp.Integer(0), bexpr), 'f', sp.Integer(0), bexpr)
        return len(nf) == 0, len(nf)
    nf = W.normal_form(W.screening_residue(density, a, bexpr), 'f', a, bexpr)
    bad = 0
    for c in nf.values():
        num = sp.Poly(sp.numer(sp.together(c)), a)
        e = sum(cf * qexpr ** (k // 2) * (a if k % 2 else 1) for (k,), cf in num.terms())
        e = sp.Poly(sp.numer(sp.together(e)), a)
        if any(sp.cancel(cf) != 0 for cf in e.coeffs()):
            bad += 1
    return bad == 0, bad


PAIRS = {1: [('S1  beta = -p/2, a^2 = 1 - p^2/4', -p / 2, 1 - p**2 / 4)],
         2: [('S2  beta = -p,   a^2 = 1 - p^2', -p, 1 - p**2)],
         3: [('S3  beta = 1/p,  a^2 = 1/p^2', 1 / p, 1 / p**2),
             ('S3\' beta = -p,   a^2 = p^2', -p, p**2)]}
VIRA = [('Virasoro beta_+ = p, a = 0', p, sp.Integer(0)), ('Virasoro beta_- = -2/p, a = 0', -2 / p, sp.Integer(0))]
NEG = [('negative: beta = p/3, a = 0', p / 3, sp.Integer(0)),
       ('negative: beta = -p/2, a^2 = 1 - p^2/3', -p / 2, 1 - p**2 / 3),
       ('negative: mirror beta = +p/2, a^2 = 1 - p^2/4', p / 2, 1 - p**2 / 4)]
fails = 0
for sol in (1, 2, 3):
    d3, d5 = densities(sol)
    print('\nSol %d: two-boson images symbolic in p: I_3 %d monomials, I_5 %d monomials  (%.0fs)'
          % (sol, len(d3), len(d5), time.time() - t0), flush=True)
    for name, be, qe in PAIRS[sol] + VIRA:
        r3 = conserved(d3, be, qe)
        r5 = conserved(d5, be, qe)
        print('   %-42s I_3 conserved identically: %s   I_5 conserved identically: %s  (%.0fs)'
              % (name, r3[0], r5[0], time.time() - t0), flush=True)
        fails += (not r3[0]) + (not r5[0])
    for name, be, qe in NEG:
        r3 = conserved(d3, be, qe)
        r5 = conserved(d5, be, qe)
        print('   %-42s I_3: %s (%d conditions nonzero)   I_5: %s (%d nonzero)'
              % (name, r3[0], r3[1], r5[0], r5[1]), flush=True)
        fails += r3[0] or r5[0]
    # cross-solution negative: another solution's pair
    other = PAIRS[1 if sol != 1 else 2][0]
    r3 = conserved(d3, other[1], other[2])
    r5 = conserved(d5, other[1], other[2])
    print('   %-42s I_3: %s   I_5: %s' % ('negative: the pair of Sol %d' % (1 if sol != 1 else 2), r3[0], r5[0]), flush=True)
    fails += r3[0] and r5[0]
print('\nRESULT:', 'FAILED (%d)' % fails if fails else 'all listed pairs conserved against I_3 and I_5 identically in p; negatives fail')
sys.exit(1 if fails else 0)
