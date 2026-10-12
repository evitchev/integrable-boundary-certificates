"""VIR2 main (Fable seat).  Seal SEAL_VIR2.md (hash guard below).

Parametrisation used on the machine: p = beta_+ (so p^2 = 2u^2), t = (2+p^2)/(2-p^2),
rho = p/2 - 1/p, N = 1 - 12 rho^2.  Pairs under test, V = exp(+-a X + beta phi):
   S1  (Sol 1) beta = -p/2, a^2 = 1 - p^2/4        [sealed]
   S3  (Sol 3) beta = 1/p,  a^2 = 1/p^2            [sealed]
   S2  (Sol 2) beta = -p,   a^2 = 1 - p^2          [sealed form: level-3 degenerate, a^2 -> -1 at N = 1]
   S3' (Sol 3) beta = -p,   a^2 = p^2              [found at spin 3; NOT sealed -- S-complete said I_5 kills it]
Exit 0 = S1 and S3 hold against I_3, I_5 and the I_7 hold-out with all controls firing;
     2 = fail-fast or S1/S3 refuted;  1 = a gate or control failed.
"""
import hashlib, json, pickle, sys, time
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
sys.path.insert(0, '<repo>/code')
import vir_lib as V
import vir2_lib as W
from vir2_lib import a, b, rho

SEAL = '9ac0c5560861245c9e33ee228a96c717fa50555b2bf63d9fdde66633b483a499'
if hashlib.sha256(open('SEAL_VIR2.md', 'rb').read()).hexdigest() != SEAL:
    print('SEAL HASH MISMATCH'); sys.exit(1)
print('seal hash ok', SEAL)
TAMPER = len(sys.argv) > 1 and sys.argv[1] == 'tamper'      # scoring tamper: test the MIRROR pairs
t0 = time.time()
ctl_fail, pred_fail = [], []
def ctl(c, msg):
    print(('PASS  ' if c else 'CONTROL-FAIL  ') + msg, flush=True)
    if not c: ctl_fail.append(msg)
def pred(c, msg):
    print(('HOLDS  ' if c else 'REFUTED  ') + msg, flush=True)
    if not c: pred_fail.append(msg)

t = sp.Symbol('t')
q = sp.Symbol('q')
I5SYM = {k: [sp.sympify(x) for x in v] for k, v in pickle.load(open('i5_sym.pkl', 'rb')).items()}
QT = json.load(open('<repo>/results/qt_promotion.json'))
Ns = sp.Symbol('N')
B4 = V.basis(4); i4 = {m: i for i, m in enumerate(B4)}
from mixed_engine import gen_mixed_basis
for w in (8, 10):
    assert V.basis(w) == gen_mixed_basis(w)


def fr(x):
    x = sp.nsimplify(sp.simplify(x))
    assert x.is_Rational, x
    return F(int(x.p), int(x.q))


class Fibre:
    def __init__(self, pval=None, N1=False):
        if N1:                      # the paperclip point (N, s) = (1, 0), t = infinity, rho = 0
            self.p, self.t, self.N, self.s, self.rho = sp.sqrt(2), sp.oo, F(1), F(0), sp.Integer(0)
        else:
            self.p = sp.sympify(pval)
            p2 = sp.simplify(self.p**2)
            self.t = sp.simplify((2 + p2) / (2 - p2))
            self.rho = sp.simplify(self.p / 2 - 1 / self.p)
            self.N = fr(1 - 12 * self.rho**2)
            self.s = fr(3 * p2 - 12 / p2)
            tt = fr(self.t)
            assert V.curve_point(tt) == (self.N, self.s)

    def sub(self, d, rho_val=None):
        rv = self.rho if rho_val is None else rho_val
        return W.clean({m: sp.expand(c.subs(rho, rv)) for m, c in d.items()})

    def I3(self, sol, rho_val=None):
        return self.sub(W.to_two_boson(V.to_row(V.P4(sol, self.N, self.s), i4), 4)[0], rho_val)

    def I5mixed(self, sol):
        if self.t == sp.oo:
            deg = max(sp.Poly(c, t).degree() for c in I5SYM[sol] if c != 0)
            return [sp.Poly(c, t).coeff_monomial(t**deg) if c != 0 else sp.Integer(0) for c in I5SYM[sol]]
        return [c.subs(t, self.t) for c in I5SYM[sol]]

    def I5(self, sol, rho_val=None):
        return self.sub(W.to_two_boson(self.I5mixed(sol), 6)[0], rho_val)

    def Ipromoted(self, sol, sigma, rho_val=None):
        data = QT['cyl%d_sigma%d' % (sol, sigma)]
        Bw = V.basis(sigma + 1)
        vec = [sp.Integer(0)] * len(Bw)
        Nq, sq = sp.Rational(self.N.numerator, self.N.denominator), sp.Rational(self.s.numerator, self.s.denominator)
        for i, (A_, B_) in data['coeffs_A_B_of_N'].items():
            val = sp.sympify(A_, locals={'N': Ns}).subs(Ns, Nq) + sq * sp.sympify(B_, locals={'N': Ns}).subs(Ns, Nq)
            assert val.is_Rational
            vec[int(i)] = val
        return self.sub(W.to_two_boson(vec, sigma + 1)[0], rho_val)


def zero(e):
    return sp.simplify(sp.expand(e)) == 0


def conserved(density, bval, qval, elim='f'):
    """V = exp(a X + bval phi), a^2 = qval.  Returns (bool, #nonzero normal-form coefficients)."""
    if elim == 'x':
        nf = W.normal_form(W.screening_residue(density, a, bval), 'x', a, bval)
    elif qval == 0:
        nf = W.normal_form(W.screening_residue(density, sp.Integer(0), bval), 'f', sp.Integer(0), bval)
        return len([1 for c in nf.values() if not zero(c)]) == 0, len(nf)
    else:
        nf = W.normal_form(W.screening_residue(density, a, bval), 'f', a, bval)
    bad = 0
    for c in nf.values():
        num = sp.Poly(sp.numer(sp.together(c)), a)
        e = sp.Poly(sum(cf * qval ** (k // 2) * (a if k % 2 else 1) for (k,), cf in num.terms()), a)
        if any(not zero(cf) for cf in e.coeffs()):
            bad += 1
    return bad == 0, bad


def PAIRS(sol, p):
    if sol == 1:
        return [('S1', -p / 2, 1 - p**2 / 4)]
    if sol == 2:
        return [('S2', -p, 1 - p**2)]
    return [('S3', 1 / p, 1 / p**2), ("S3'", -p, p**2)]


def joint_solutions(fib, sol, dens):
    """All (q = a^2, b) with a != 0, b != 0 solving every condition of the given densities."""
    polys = []
    for d in dens:
        for c in W.conditions(d):
            pa = sp.Poly(sp.expand(c), a)
            degs = {k % 2 for (k,), _ in pa.terms()}
            assert len(degs) == 1
            sh = 1 if degs == {1} else 0
            polys.append(sp.expand(sum(cf * q ** ((k - sh) // 2) for (k,), cf in pa.terms())))
    gaussian = any(pl.has(sp.I) for pl in polys)
    G = sp.groebner(polys, q, b, order='lex', gaussian=True) if gaussian else sp.groebner(polys, q, b, order='lex')
    sols = sp.solve(list(G.exprs), [q, b], dict=True)
    return sorted({(sp.nsimplify(s_[q]), sp.nsimplify(s_[b])) for s_ in sols}, key=str), len(polys)


# ---------------------------------------------------------------- fibres
REAL = [sp.Rational(x) for x in ('1/3', '-1/3', '2/3', '-2/3', '3/2', '-3/2', '2', '-2', '3', '-3', '4/3', '-4/3',
                                 '5/2', '-5/2', '1/4', '-1/4', '5/3', '-5/3', '7/2', '-7/2', '3/4', '-3/4', '6/5', '-6/5')]
IMAG = [sp.I * sp.Rational(x) for x in ('1/2', '3/2', '3', '2/3', '5/2', '1/3')]
print('%d real-p fibres (N < 1) + %d imaginary-p fibres (N > 25); t < 0 on %d of them'
      % (len(REAL), len(IMAG), sum(1 for p_ in REAL + IMAG if Fibre(p_).t < 0)))

# ---------------------------------------------------------------- C0 gate: [I_3, I_5] = 0 in the two-boson image
print('\n== C0: the two-boson images commute ==')
def commute(d3, d5):
    res = W.ope2(d3, d5, -1)
    return W.is_total_derivative(res, 9), len(res)
for sol in (1, 2, 3):
    fib = Fibre(sp.Rational(3, 2))
    ok, n = commute(fib.I3(sol), fib.I5(sol))
    ctl(ok and n > 0, 'Sol %d at p = 3/2: Res I_3(z) I_5(w) (%d monomials) is a total derivative  (%.0fs)' % (sol, n, time.time() - t0))

# ---------------------------------------------------------------- main: I_3 and I_5 on all fibres
print('\n== sealed pairs against I_3 and I_5; a = 0 control; joint solution sets ==')
tally = {}
for p_ in REAL + IMAG:
    fib = Fibre(p_)
    for sol in (1, 2, 3):
        d3, d5 = fib.I3(sol), fib.I5(sol)
        for name, bv, qv in PAIRS(sol, fib.p):
            if TAMPER:
                bv = -bv
            r3, r5 = conserved(d3, bv, qv)[0], conserved(d5, bv, qv)[0]
            k = tally.setdefault(name, [0, 0, 0])
            k[0] += 1; k[1] += r3; k[2] += r5
        for bv in (fib.p, -2 / fib.p):
            r = conserved(d3, bv, 0)[0] and conserved(d5, bv, 0)[0]
            k = tally.setdefault('Virasoro a=0', [0, 0, 0]); k[0] += 1; k[1] += r; k[2] += r
        r = conserved(d3, fib.p / 3, 0)[0] or conserved(d5, fib.p / 3, 0)[0]
        k = tally.setdefault('generic beta = p/3, a=0 (must fail)', [0, 0, 0]); k[0] += 1; k[1] += r; k[2] += r
    print('   p = %-6s t = %-8s N = %-10s done  (%.0fs)' % (fib.p, fib.t, fib.N, time.time() - t0), flush=True)
for name, (n, n3, n5) in tally.items():
    print('   %-38s fibres x sols tested %3d: conserved vs I_3 %3d, vs I_5 %3d' % (name, n, n3, n5))
n1 = tally['S1']; n3_ = tally['S3']; n2 = tally['S2']; n3p = tally["S3'"]
pred(n1[1] == n1[0] and n1[2] == n1[0], 'S1: Solution 1 pair conserved against I_3 and I_5 at all %d fibres' % n1[0])
pred(n3_[1] == n3_[0] and n3_[2] == n3_[0], 'S3: Solution 3 pair conserved against I_3 and I_5 at all %d fibres' % n3_[0])
print(('HOLDS  ' if n2[1] == n2[0] == n2[2] else 'REFUTED (not exit-setting)  ')
      + 'S2: Solution 2 has a level-3-degenerate pair (beta = -beta_+, a^2 = 1 - p^2 -> -1 at N = 1), I_3 and I_5, all %d fibres' % n2[0])
print(('S-complete REFUTED on Sol 3 (not exit-setting): ' if n3p[2] == n3p[0] else 'S-complete survives: ')
      + "the second Sol-3 pair S3' (beta = -p, a^2 = p^2) is conserved against I_5 at %d of %d fibres" % (n3p[2], n3p[0]))
v = tally['Virasoro a=0']; gneg = tally['generic beta = p/3, a=0 (must fail)']
ctl(v[1] == v[0], 'C1: both Virasoro screenings pass I_3 and I_5 at every fibre, every solution (%d tests)' % v[0])
ctl(gneg[1] == 0, 'C1: a generic pure-phi exponent fails at every fibre')

print('\n== completeness: all (a^2, beta), a != 0, beta != 0 (Groebner + solve) ==')
complete = True
for p_ in [sp.Rational(3, 2), sp.Rational(-2, 3), sp.Rational(3), sp.Rational(1, 4), sp.I * sp.Rational(3, 2), sp.I * 3]:
    fib = Fibre(p_)
    for sol in (1, 2, 3):
        s3, n3c = joint_solutions(fib, sol, [fib.I3(sol)])
        s35, n35 = joint_solutions(fib, sol, [fib.I3(sol), fib.I5(sol)])
        expect = sorted({(sp.nsimplify(qv), sp.nsimplify(bv)) for _, bv, qv in PAIRS(sol, fib.p)}, key=str)
        same = (s35 == expect) and (s3 == expect)
        complete &= same
        print('   p = %-5s Sol %d: I_3 alone (%d conditions): %s ; I_3 and I_5 (%d conditions): %s ; == listed pairs: %s  (%.0fs)'
              % (fib.p, sol, n3c, s3, n35, s35, same, time.time() - t0), flush=True)
print(('HOLDS  ' if complete else 'DIFFERS  ') + 'the mixed solution set of {I_3} and of {I_3, I_5} is exactly S1 | S2 | {S3, S3\'} at 6 fibres')

# ---------------------------------------------------------------- hold-outs I_7 (sealed) and I_9 (bonus)
print('\n== hold-outs: I_7 (sealed), I_9 (post-seal bonus), promoted charges ==')
hold = {}
for sigma, plist in ((7, [sp.Rational(3, 2), sp.Rational(-2, 3), sp.Rational(3), sp.Rational(-5, 2), sp.I * sp.Rational(3, 2), sp.I / 2]),
                     (9, [sp.Rational(3, 2), sp.Rational(-2, 3), sp.I * sp.Rational(3, 2)])):
    for p_ in plist:
        fib = Fibre(p_)
        for sol in (1, 2, 3):
            d = fib.Ipromoted(sol, sigma)
            for name, bv, qv in PAIRS(sol, fib.p):
                if TAMPER:
                    bv = -bv
                r = conserved(d, bv, qv)
                k = hold.setdefault((sigma, name), [0, 0]); k[0] += 1; k[1] += r[0]
            rv = conserved(d, fib.p, 0)[0] and conserved(d, -2 / fib.p, 0)[0]
            k = hold.setdefault((sigma, 'Virasoro'), [0, 0]); k[0] += 1; k[1] += rv
            neg = conserved(d, -fib.p / 2, 1 - fib.p**2 / 3)
            k = hold.setdefault((sigma, 'negative (wrong a^2)'), [0, 0]); k[0] += 1; k[1] += neg[0]
        print('   I_%d at p = %-5s (%d-monomial images) done  (%.0fs)' % (sigma, fib.p, len(d), time.time() - t0), flush=True)
for (sigma, name), (n, k) in sorted(hold.items(), key=str):
    print('   I_%d  %-22s conserved at %d of %d' % (sigma, name, k, n))
pred(hold[(7, 'S1')][0] == hold[(7, 'S1')][1], 'S1 hold-out: conserved against I_7 at all tested fibres')
pred(hold[(7, 'S3')][0] == hold[(7, 'S3')][1], 'S3 hold-out: conserved against I_7 at all tested fibres')
ctl(hold[(7, 'Virasoro')][0] == hold[(7, 'Virasoro')][1], 'C1: Virasoro screenings pass I_7')
ctl(hold[(7, 'negative (wrong a^2)')][1] == 0, 'C5: a wrong a^2 fails I_7 everywhere')

# ---------------------------------------------------------------- C2 prior art
print('\n== C2: prior-art controls ==')
fib = Fibre(2 * sp.I)                      # exact point, Sol 1: t = -1/3
ctl(fib.t == sp.Rational(-1, 3) and fib.N == 28, 'p = 2i is the exact point of Sol 1: t = %s, N = %s, rho = %s' % (fib.t, fib.N, fib.rho))
r3 = conserved(fib.I3(1), sp.Integer(0), sp.Integer(2), elim='x'); r5 = conserved(fib.I5(1), sp.Integer(0), sp.Integer(2), elim='x')
ctl(r3[0] and r5[0], 'Sol 1 at (28,-9): pure-X (a^2, beta) = (2, 0) conserved against I_3 and I_5')
r3n = conserved(fib.I3(1), sp.Integer(0), sp.Integer(3), elim='x')
ctl(not r3n[0], 'Sol 1 at (28,-9): a^2 = 3 is NOT (non-root control)')
fg = Fibre(sp.Rational(3, 2))
r3 = conserved(fg.I3(1), sp.Integer(0), sp.Integer(2), elim='x')
ctl(not r3[0], 'Sol 1 at a generic fibre: pure-X a^2 = 2 is NOT conserved (%d conditions nonzero)' % r3[1])
f1 = Fibre(N1=True)
for sol in (1, 3):
    for sg in (1, -1):
        r3 = conserved(f1.I3(sol), sg / sp.sqrt(2), sp.Rational(1, 2)); r5 = conserved(f1.I5(sol), sg / sp.sqrt(2), sp.Rational(1, 2))
        ctl(r3[0] and r5[0], 'N = 1, Sol %d: BOP\'s (a^2, beta^2) = (1/2, 1/2), beta = %s/sqrt2: I_3 and I_5' % (sol, '+1' if sg > 0 else '-1'))
r3 = conserved(f1.I3(2), -sp.sqrt(2), sp.Integer(-1)); r5 = conserved(f1.I5(2), -sp.sqrt(2), sp.Integer(-1))
ctl(r3[0] and r5[0], 'N = 1, Sol 2: (a^2, beta^2) = (-1, 2) [BOP\'s N = 1 value for Solution 2]: I_3 and I_5')
for sol in (1, 3):
    r3 = conserved(f1.I3(sol), -sp.sqrt(2), sp.Integer(2)); r5 = conserved(f1.I5(sol), -sp.sqrt(2), sp.Integer(2))
    print('   N = 1, Sol %d: the branch (a^2, beta^2) = (2, 2): I_3 %s, I_5 %s' % (sol, r3[0], r5[0]))

# ---------------------------------------------------------------- C3 tampers
print('\n== C3: tampers ==')
fib = Fibre(sp.Rational(3, 2))
for sol in (1, 2, 3):
    name, bv, qv = PAIRS(sol, fib.p)[0]
    d3m, d5m = fib.I3(sol, -fib.rho), fib.I5(sol, -fib.rho)
    fail_sealed = not (conserved(d3m, bv, qv)[0] or conserved(d5m, bv, qv)[0])
    mirror = conserved(d3m, -bv, qv)[0] and conserved(d5m, -bv, qv)[0]
    ctl(fail_sealed and mirror, 'tamper A (rho -> -rho), Sol %d: %s fails, its mirror beta -> -beta passes' % (sol, name))
    rB = sp.sqrt(fib.rho**2 - sp.Rational(1, 12))          # rho^2 from N + 1
    d3b, d5b = fib.I3(sol, rB), fib.I5(sol, rB)
    f3, f5 = conserved(d3b, bv, qv), conserved(d5b, bv, qv)
    ctl(not f3[0] and not f5[0], 'tamper B (rho^2 from N+1, rho = %s), Sol %d: %s fails I_3 (%d) and I_5 (%d conditions nonzero)'
        % (rB, sol, name, f3[1], f5[1]))
    # with the wrong central charge the images stop commuting (rational mismatch rho + 1/5 for the rank test)
    rC = fib.rho + sp.Rational(1, 5)
    okc, _ = commute(fib.I3(sol, rC), fib.I5(sol, rC))
    ctl(not okc, 'tamper B\' (rho -> rho + 1/5, wrong central charge), Sol %d: the two-boson images of I_3, I_5 no longer commute' % sol)

# ---------------------------------------------------------------- C4 plant, C5 vacuity
print('\n== C4 plant, C5 vacuity ==')
for sol in (1, 2, 3):
    name, bv, qv = PAIRS(sol, fib.p)[0]
    d5 = dict(fib.I5(sol))
    W.add_to(d5, ((1, 1), (2, 2)), sp.Rational(1, 1000))
    r = conserved(d5, bv, qv)
    ctl(not r[0], 'plant, Sol %d: I_5 + (1/1000) x1^2 f2^2 is flagged for %s (%d conditions nonzero)' % (sol, name, r[1]))
    r3, r5 = conserved(fib.I3(sol), sp.Rational(2, 7), sp.Rational(5, 3)), conserved(fib.I5(sol), sp.Rational(2, 7), sp.Rational(5, 3))
    ctl(not r3[0] and not r5[0], 'vacuity, Sol %d: a generic (a^2, beta) = (5/3, 2/7) fails: %d of I_3\'s and %d of I_5\'s conditions nonzero'
        % (sol, r3[1], r5[1]))

print('\ntotal %.0fs' % (time.time() - t0))
if ctl_fail:
    print('RESULT: CONTROL FAILURE (%d)' % len(ctl_fail)); sys.exit(1)
if pred_fail:
    print('RESULT: SEALED PREDICTION REFUTED (%d)' % len(pred_fail)); sys.exit(2)
print('RESULT: S1 and S3 CONFIRMED against I_3, I_5 and the I_7 hold-out; all controls fired')
sys.exit(0)
