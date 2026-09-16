"""FIRST-ORDER DEFORMATION OF THE EXACT OPER: five scalar classes excluded by PROVED integrand identities (ansatz note sec. 38(ah)
addendum 1, batches 25-26; Codex reviews of 694a78a and f87dd4f).  At the exact point (N, s) = (28, -9) -- Solution 1 at t0 = -1/3,
where the sec.-38(ad) operator (the pillow-brane vacuum oper at P_1^2 = 1/6) reproduces the certified vacuum eigenvalues
(cyl_exact_point.py) -- consider the operator's deformation along the O(N) curve, t = t0 + h, and a first-order change of its
potential in one of the classes
    (W)   h W(w),   (XY)  h [X f(w) + Y g(w)],   (WXY) h [W(w) + X f(w) + Y g(w)],
    (U)   h U(w)/(kappa^2 Lambda),   (Q) (h/kappa^2 Lambda) [X^2 g_XX(w) + XY g_XY(w) + Y^2 g_YY(w)],
    (L)   the SYMBOL itself, kappa^2 Lambda -> kappa^2 Lambda (1 + h sigma(w))  (any change of the exponents' shape beyond the pillow's
          n(t), eps(t), whose own derivatives are in the targets),
    (Q0)  h [X^2 g_XX(w) + XY g_XY(w) + Y^2 g_YY(w)] at order kappa^0 (a first-order change of the Kac dictionary; it creates a
          degree-(k+1) layer, and its spin-3 responses enter the top degree only),
    (KODD) the kappa-ODD term h kappa Q(w), analytic in h: by parity ((kappa, h) -> (-kappa, -h) exchanges the two WKB branches) its
          first-order odd-spin response is a total derivative, so r(3,0,0) = 0 identically for the class (Codex) -- the relation is
          the single term r(3,0,0) = 0, proved for symbolic exponents with Q = sqrt(Lambda) F absorbed into F's exponents,
with each function a finite combination of monomials w^i (1-w)^j of ANY exponents (Codex f87dd4f: powers of Lambda are absorbed
into the exponents, Lambda = w^(3/2)(1-w) at the point).  The X^k-normalized first-order response r(s, a, b) of the operator's
continued-Beta coefficient I_s (coefficient of X^a Y^b of d/dh [I_s/I_s[X^k]] at h = 0) must VANISH at momentum degrees k and k-1
(the top and sub-leading layers are exact identically in t) and must equal the first-order target M_s'(t0) at degrees <= k-2,
where M_s = e_s/e_s[X^k] - I_s/I_s[X^k] and e_s is the certified vacuum eigenvalue (rational in t).
This certificate establishes, with exact arithmetic throughout (seven potential classes, the kappa-odd class and the symbol class):
 (1) TARGETS: M_s(t0) = 0 and M_s'(t0) for s = 3, 5, 7, 9 from the certified weight-4..10 tables (d/dt by sympy) and the
     operator's coefficients by dual-number arithmetic (a + b h, h^2 = 0) through the Riccati/continued-Beta machinery;
 (2) IDENTITIES: for each class, a linear relation among the r(s, a, b) that holds for the shape w^i (1-w)^j with SYMBOLIC
     exponents (i, j) -- every continued Beta integral at spin s is B_s(i, j) times a rational function of (i, j), and
     B_(s+4)/B_s is rational, so a relation among spins congruent mod 4 is B_s times a rational function, which is shown to
     vanish identically -- for every kind of the class (W; X, Y; W, X, Y; U; XX, XY, YY): hence for every member of the class;
 (3) VIOLATION: with the degree-k and degree-(k-1) responses set to zero, each relation reduces to a condition on the
     degree-(k-2) responses that the targets of (1) violate (the discrepancy is printed exactly).
So no member of any of the five classes is the first-order deformation of the exact oper along the curve.  Fail-closed: any
identity that does not reduce to zero, any target that satisfies a relation, or any nonzero M_s(t0) exits 1.
--control: (a) a wrong coefficient in a relation is rejected by the prover; (b) the W relation applied to the X-shapes is NOT an
identity (the prover distinguishes kinds); (c) ZERO targets satisfy every relation (the certificate can pass); (d) the actual targets
violate the W relation; (e) in the symbol class, sigma = 1 (a rescaling of kappa) has identically zero normalized response; (f) the kappa-odd
source misplaced into rho_1 breaks the parity identity (Kimi: a doubled amplitude cannot, the identity being homogeneous); (g) the local
Riccati balance fixes the source amplitude to F/2 (Codex).  ~2 min.
Usage: cyl_exact_point_deformation.py [--control]"""
import sys, os, time, itertools
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from cyl_first_order import X, Y, cyl, t
CONTROL = '--control' in sys.argv
failures = []
def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok: failures.append(name)
R = lambda v: sp.Rational(v.numerator, v.denominator)
T0 = F(-1, 3); SMAX = 9
# ---------------- (1) targets: dual numbers over Q through the machinery at t = t0 + h ----------------
class D:
    __slots__ = ('a', 'b')
    def __init__(s, a, b=F(0)): s.a = F(a); s.b = F(b)
    def __add__(s, o): o = o if isinstance(o, D) else D(o); return D(s.a + o.a, s.b + o.b)
    __radd__ = __add__
    def __sub__(s, o): o = o if isinstance(o, D) else D(o); return D(s.a - o.a, s.b - o.b)
    def __rsub__(s, o): return D(o) - s
    def __mul__(s, o): o = o if isinstance(o, D) else D(o); return D(s.a * o.a, s.a * o.b + s.b * o.a)
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = o if isinstance(o, D) else D(o); return D(s.a / o.a, s.b / o.a - s.a * o.b / (o.a * o.a))
    def __rtruediv__(s, o): return D(o) / s
    def __neg__(s): return D(-s.a, -s.b)
    def iszero(s): return s.a == 0 and s.b == 0
def dadd(p, q, s=D(1)):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, D(0)) + s * v
        if r[k].iszero(): del r[k]
    return r
def dmul(p, q):
    r = {}
    for (i1, a1, b1), v1 in p.items():
        for (i2, a2, b2), v2 in q.items():
            k = (i1 + i2, a1 + a2, b1 + b2); r[k] = r.get(k, D(0)) + v1 * v2
    return {k: v for k, v in r.items() if not v.iszero()}
def dscale(p, s): return {k: s * v for k, v in p.items()}
def ddw(p):
    r = {}
    for (i, a, b), v in p.items():
        if i == 0: continue
        for k, s in (((i, a, b), i * v), ((i + 1, a, b), -i * v)): r[k] = r.get(k, D(0)) + s
    return {k: v for k, v in r.items() if not v.iszero()}
def dpoch(a, i):
    r = D(1)
    for j in range(i): r = r * (a + j)
    return r
def operator_dual(t0, smax):
    """Solution 1 parameters as dual numbers at t = t0 + h; returns {s: {(a, b): D}} = the X^k-normalized I_s (base, h-part)"""
    tv = D(t0, F(1))
    n, eps, cX, d2 = -(tv + 7) / (tv + 3), D(-4) / (tv + 3), D(1), D(-1) / ((tv + 1) * (tv + 3))
    kap = -(n + 2) / (n - eps); c0 = -eps / 24
    onemw = {(0, 0, 0): D(1), (1, 0, 0): D(-1)}
    V = {(1, 1, 0): -cX / 2}; V = dadd(V, dmul({(1, 0, 1): kap * cX / 2, (1, 0, 0): F(1, 4) + d2}, onemw)); V = dadd(V, dscale(onemw, c0))
    ll = {(0, 0, 0): -eps / 2, (1, 0, 0): n / 2}
    def dz(k, r): return dadd(ddw(r), dmul(r, ll), D(-k))
    rho = {0: dscale(ll, D(F(-1, 2)))}
    rho[1] = dscale(dadd(dadd(V, dmul(rho[0], rho[0]), D(-1)), dz(0, rho[0]), D(-1)), D(F(1, 2)))
    for k in range(1, smax):
        acc = dz(k, rho[k])
        for i in range(0, k + 1): acc = dadd(acc, dmul(rho[i], rho[k - i]))
        rho[k + 1] = dscale(acc, D(F(-1, 2)))
    out = {}
    for s in range(3, smax + 1, 2):
        k = (s + 1) // 2; a0, ab = s * eps / 2, s * n / 2; I = {}
        for (i, a, b), v in rho[s].items():
            den = dpoch(ab, i)
            if den.a == 0: raise ZeroDivisionError("continuation lattice")
            I[(a, b)] = I.get((a, b), D(0)) + v * dpoch(a0, i) / den
        T = I[(k, 0)]; out[s] = {ab_: v / T for ab_, v in I.items()}
    return out
t1 = time.time()
op = operator_dual(T0, SMAX)
tq = R(T0); targets = {}
for s in range(3, SMAX + 1, 2):
    k = (s + 1) // 2
    c = sp.Poly(sp.expand(cyl(1, s + 1)), X, Y, domain=sp.QQ.frac_field(t)); top = c.coeff_monomial(X**k)
    En = sp.cancel(sp.expand(cyl(1, s + 1)) / sp.QQ.frac_field(t).to_sympy(top))
    E0 = sp.Poly(sp.expand(En.subs(t, tq)), X, Y); E1 = sp.Poly(sp.expand(sp.diff(En, t).subs(t, tq)), X, Y)
    I0 = sum(R(v.a) * X**a * Y**b for (a, b), v in op[s].items()); I1 = sum(R(v.b) * X**a * Y**b for (a, b), v in op[s].items())
    M0 = sp.expand(E0.as_expr() - I0); M1 = sp.Poly(sp.expand(E1.as_expr() - I1), X, Y)
    check(f"(1) spin {s}: the mismatch vanishes at the exact point, M_{s}(t0) = 0", M0 == 0)
    check(f"(1) spin {s}: the first-order target M_{s}'(t0) has degree <= k-2 = {k - 2} (top and sub-leading layers exact in t)", M1.total_degree() <= k - 2 or M1.is_zero)
    targets[s] = M1
print(f"   targets: M_3' = {targets[3].as_expr()}; M_5' = {targets[5].as_expr()} ({time.time() - t1:.0f}s)", flush=True)
if CONTROL:
    pass
# ---------------- (2) the identities, proved for symbolic exponents ----------------
ii, jj = sp.symbols('i j')
gens = (ii, jj)
def P(c): return sp.Poly(c, *gens, domain='QQ')
def padd(p, q, s=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, P(0)) + s * v
        if r[k].is_zero: del r[k]
    return r
def pmul(p, q):
    r = {}
    for (i1, j1, a1, b1, f1), v1 in p.items():
        for (i2, j2, a2, b2, f2), v2 in q.items():
            if f1 + f2 > 1: continue
            k = (i1 + i2, j1 + j2, a1 + a2, b1 + b2, f1 + f2); r[k] = r.get(k, P(0)) + v1 * v2
    return {k: v for k, v in r.items() if not v.is_zero}
def pscale(p, s): return {k: v * s for k, v in p.items()}
Pi, Pj = P(ii), P(jj)
def pDw(p):
    r = {}
    for (di, dj, a, b, f), v in p.items():
        I = P(di) + (Pi if f else P(0)); J = P(dj) + (Pj if f else P(0))
        for key, coef in (((di, dj + 1, a, b, f), I * v), ((di + 1, dj, a, b, f), -J * v)):
            r[key] = r.get(key, P(0)) + coef
    return {k: v for k, v in r.items() if not v.is_zero}
def rf(x, m):
    if m >= 0: return sp.Mul(*[x + q for q in range(m)]) if m else sp.Integer(1)
    return 1 / sp.Mul(*[x - q for q in range(1, -m + 1)])
def poch_rat(a, m):
    r = F(1)
    if m >= 0:
        for q in range(m): r *= a + q
    else:
        for q in range(1, -m + 1): r /= (a - q)
    return R(r)
n0, eps0, cX0, d20 = F(-5, 2), F(-3, 2), F(1), F(-9, 16); kap0 = -(n0 + 2) / (n0 - eps0); c00 = -eps0 / 24
MOM = {'W': (0, 0), 'X': (1, 0), 'Y': (0, 1), 'U': (0, 0), 'XX': (2, 0), 'XY': (1, 1), 'YY': (0, 2), 'L': (0, 0), 'XX0': (2, 0), 'XY0': (1, 1), 'YY0': (0, 2), 'KODD': (0, 0)}   # L: the symbol, kappa^2 Lambda -> kappa^2 Lambda (1 + h sigma(w)); XX0, XY0, YY0: kappa^0 momentum-quadratic (the Kac dictionary at first order); KODD: the kappa-odd term h kappa Q, Q = sqrt(Lambda) F (source F/2 in rho_0)
SRCK = ('U', 'XX', 'XY', 'YY')      # kappa^-2 Lambda^-1 sources: enter rho_3 as shape X^a Y^b / 2
_rho_cache = {}
def symbolic_rho(kind, smax):
    if (kind, smax) in _rho_cache: return _rho_cache[(kind, smax)]
    a_, b_ = MOM[kind]; src = kind in SRCK; sym = kind == 'L'; kodd = kind == 'KODD'
    ll = {(0, 0, 0, 0, 0): P(R(-eps0) / 2), (1, 0, 0, 0, 0): P(R(n0) / 2)}
    if sym:   # lamlog -> lamlog + h D sigma / 2, sigma = w^i (1-w)^j with symbolic exponents: D sigma = i w^i (1-w)^(j+1) - j w^(i+1) (1-w)^j
        ll = padd(ll, {(0, 1, 0, 0, 1): Pi * sp.Rational(1, 2), (1, 0, 0, 0, 1): Pj * sp.Rational(-1, 2)})
    V = {(1, 0, 1, 0, 0): P(-R(cX0) / 2), (1, 1, 0, 1, 0): P(R(kap0 * cX0) / 2), (1, 1, 0, 0, 0): P(sp.Rational(1, 4) + R(d20)), (0, 1, 0, 0, 0): P(R(c00))}
    if not src and not sym and not kodd: V[(0, 0, a_, b_, 1)] = P(1)
    def dz(k, r): return padd(pDw(r), pmul(r, ll), -k)
    rho = {0: pscale(ll, sp.Rational(-1, 2))}
    if kodd: rho[0] = padd(rho[0], {(0, 0, 0, 0, 1): P(sp.Rational(1, 2))})     # h kappa Q, Q = sqrt(Lambda) F: rho_0 += F/2 (Codex, review of the queued batch 27)
    rho[1] = pscale(padd(padd(V, pmul(rho[0], rho[0]), -1), dz(0, rho[0]), -1), sp.Rational(1, 2))
    for k in range(1, smax):
        acc = dz(k, rho[k])
        for m in range(0, k + 1): acc = padd(acc, pmul(rho[m], rho[k - m]))
        rho[k + 1] = pscale(acc, sp.Rational(-1, 2))
        if k + 1 == 3 and src: rho[3] = padd(rho[3], {(0, 0, a_, b_, 1): P(sp.Rational(1, 2))})
    _rho_cache[(kind, smax)] = rho; return rho
def prove(kind, terms):
    """terms: [(coef, s, a, b)]; returns (is_identity, numerator polynomial)"""
    spins = sorted({s for _, s, _, _ in terms}); smax = max(spins)
    assert len({s % 4 for s in spins}) == 1
    rho = symbolic_rho(kind, smax)
    s_ref = spins[0]; a_ref, b_ref = -F(3, 4) * s_ref, -F(1, 2) * s_ref
    resp = {}
    for s in spins:
        k = (s + 1) // 2; a0, b0 = -F(3, 4) * s, -F(1, 2) * s
        da, db = int(a0 - a_ref), int(b0 - b_ref); dab = da + db
        Q = rf(R(a_ref) + ii, da) * rf(R(b_ref) + jj, db) / rf(R(a_ref + b_ref) + ii + jj, dab)
        Cs = poch_rat(a_ref + b_ref, dab) / (poch_rat(a_ref, da) * poch_rat(b_ref, db))
        base, lin = {}, {}
        rterms = dict(rho[s])
        if kind == 'L':   # Lambda_h^(-s/2) = Lambda^(-s/2)(1 - h s sigma/2): every base term acquires -(s/2) v w^(di+i) (1-w)^(dj+j)
            for (di, dj, a, b, f), v in rho[s].items():
                if not f:
                    key = (di, dj, a, b, 1); rterms[key] = rterms.get(key, P(0)) + v * sp.Rational(-s, 2)
        for (di, dj, a, b, f), v in rterms.items():
            if f:
                fac = rf(R(a0) + ii, di) * rf(R(b0) + jj, dj) / rf(R(a0 + b0) + ii + jj, di + dj)
                lin[(a, b)] = lin.get((a, b), 0) + v.as_expr() * fac
            else:
                fac = poch_rat(a0, di) * poch_rat(b0, dj) / poch_rat(a0 + b0, di + dj)
                base[(a, b)] = base.get((a, b), 0) + v.as_expr() * fac
        T = base[(k, 0)]; Tlin = lin.get((k, 0), 0)
        for (s_, a, b) in {(s_, a, b) for _, s_, a, b in terms if s_ == s}:
            resp[(s, a, b)] = sp.together((lin.get((a, b), 0) - (base.get((a, b), 0) / T) * Tlin) / T * Q * Cs)
    expr = sp.together(sum(c * resp[(s, a, b)] for c, s, a, b in terms))
    num, den = sp.fraction(expr); numP = sp.Poly(sp.expand(num), ii, jj)
    return numP.is_zero, numP
Rr = sp.Rational
RELATIONS = [
    ("W: h W(w)", ('W',), "spins 3 and 7 (batch 25)",
     [(Rr(1), 3, 0, 0), (Rr(-35, 16), 3, 0, 1), (Rr(-1, 8), 3, 1, 0), (Rr(-627, 224), 7, 0, 3)]),
    ("XY: h [X f(w) + Y g(w)]", ('X', 'Y'), "spins 3 and 7 (batch 25)",
     [(Rr(533113856, 2189105), 7, 0, 3), (Rr(-91706483, 4378210), 7, 0, 4), (Rr(-108288752, 2189105), 7, 1, 2), (Rr(-694893441, 4378210), 7, 1, 3),
      (Rr(-24989712, 2189105), 7, 2, 1), (Rr(97876372, 2189105), 7, 2, 2), (Rr(10412380, 437821), 7, 3, 0), (Rr(21865998, 437821), 7, 3, 1),
      (Rr(-32), 3, 0, 0), (Rr(212514736, 1313463), 3, 0, 1), (Rr(-464875985, 1313463), 3, 0, 2), (Rr(-39575312, 1313463), 3, 1, 0), (Rr(86570995, 1313463), 3, 1, 1)]),
    ("WXY: h [W(w) + X f(w) + Y g(w)]", ('W', 'X', 'Y'), "spins 5 and 9 (batch 25)",
     [(Rr(-16265600, 36243), 9, 0, 5), (Rr(-919360, 108729), 9, 1, 2), (Rr(-4071704, 108729), 9, 1, 3), (Rr(-38240293, 217458), 9, 1, 4),
      (Rr(919360, 108729), 9, 2, 1), (Rr(6950008, 108729), 9, 2, 2), (Rr(12414454, 36243), 9, 2, 3), (Rr(-229840, 36243), 9, 3, 0),
      (Rr(-7814560, 108729), 9, 3, 1), (Rr(-15847910, 36243), 9, 3, 2), (Rr(804440, 12081), 9, 4, 0), (Rr(17898790, 36243), 9, 4, 1),
      (Rr(-153088, 12081), 5, 0, 1), (Rr(593216, 12081), 5, 0, 2), (Rr(-1733602, 12081), 5, 0, 3), (Rr(76544, 12081), 5, 1, 0),
      (Rr(-210496, 36243), 5, 1, 1), (Rr(-31993, 36243), 5, 1, 2), (Rr(-1622296, 181215), 5, 2, 0), (Rr(-53404, 36243), 5, 2, 1)]),
    ("U: h U(w)/(kappa^2 Lambda)", ('U',), "spins 5 and 9 (Codex, review of 694a78a)",
     [(Rr(1), 5, 1, 0), (Rr(-85, 18), 9, 1, 2), (Rr(510, 18), 9, 0, 3)]),
    ("Q: (h/kappa^2 Lambda) [X^2 g_XX + XY g_XY + Y^2 g_YY]", ('XX', 'XY', 'YY'), "spin 3 (batch 26)",
     [(Rr(-32), 3, 0, 0), (Rr(212, 15), 3, 1, 0)]),
    ("L: the symbol, kappa^2 Lambda -> kappa^2 Lambda (1 + h sigma(w))", ('L',), "spin 5 (batch 27)",
     [(Rr(832, 63), 5, 0, 1), (Rr(-21151, 252), 5, 0, 3), (Rr(-416, 63), 5, 1, 0), (Rr(43823, 2520), 5, 1, 2), (Rr(2431, 945), 5, 2, 0), (Rr(-6539, 630), 5, 2, 1)]),
    ("Q0: h [X^2 g_XX(w) + XY g_XY(w) + Y^2 g_YY(w)] at order kappa^0 (the Kac dictionary at first order)", ('XX0', 'XY0', 'YY0'), "spin 3 (batch 28)",
     [(Rr(-32), 3, 0, 0), (Rr(212, 15), 3, 1, 0)]),
    ("KODD: the kappa-odd term h kappa Q(w) (analytic in h; Q = sqrt(Lambda) F)", ('KODD',), "spin 3 (Codex, review of the queued batch 27: the first-order odd-spin response of a kappa-odd source is a total derivative by parity)",
     [(Rr(1), 3, 0, 0)]),
]
def violation(terms, tg):
    """with the degree-k and (k-1) responses zero and the degree-(k-2) responses equal to the targets, the value of the relation"""
    val = 0
    for c, s, a, b in terms:
        k = (s + 1) // 2
        if a + b <= k - 2: val += c * tg[s].coeff_monomial(X**a * Y**b)
    return sp.nsimplify(val)
if CONTROL:
    t2 = time.time()
    bad = list(RELATIONS[0][3]); bad[3] = (Rr(-600, 224), 7, 0, 3)
    ok_a = not prove('W', bad)[0]; check("control (a): a wrong coefficient (-600/224 for -627/224) is REJECTED by the prover", ok_a)
    ok_b = not prove('X', RELATIONS[0][3])[0]; check("control (b): the W relation applied to the X-shapes is NOT an identity (kinds are distinguished)", ok_b)
    zero = {s: sp.Poly(0, X, Y) for s in targets}
    ok_c = all(violation(terms, zero) == 0 for _, _, _, terms in RELATIONS); check("control (c): ZERO targets satisfy every relation (the certificate can pass)", ok_c)
    ok_d = violation(RELATIONS[0][3], targets) != 0; check("control (d): the actual targets violate the W relation", ok_d)
    # (e) sigma = 1 in the symbol class is a rescaling of kappa: every normalized response vanishes -> the relation's terms all vanish at i = j = 0
    rho = symbolic_rho('L', 5); zero_resp = True
    for s_ in (3, 5):
        k = (s_ + 1) // 2; a0, b0 = -F(3, 4) * s_, -F(1, 2) * s_
        base, lin = {}, {}
        rterms = dict(rho[s_])
        for (di, dj, a, b, f), v in rho[s_].items():
            if not f:
                key = (di, dj, a, b, 1); rterms[key] = rterms.get(key, P(0)) + v * sp.Rational(-s_, 2)
        for (di, dj, a, b, f), v in rterms.items():
            val = v.as_expr().subs({ii: 0, jj: 0})
            fac = poch_rat(a0, di) * poch_rat(b0, dj) / poch_rat(a0 + b0, di + dj)
            (lin if f else base)[(a, b)] = (lin if f else base).get((a, b), 0) + val * fac
        T = base[(k, 0)]; Tlin = lin.get((k, 0), 0)
        for ab_ in base:
            if sp.simplify(lin.get(ab_, 0) - base[ab_] / T * Tlin) != 0: zero_resp = False
    check("control (e): sigma = 1 (a rescaling of kappa) has identically zero normalized response at spins 3 and 5", zero_resp)
    # (f) the kappa-odd wiring (Kimi d609d6a F4; Codex d609d6a): a doubled source amplitude cannot break a homogeneous zero identity, but a
    # MISPLACED source (entering rho_1 instead of rho_0) must -- the parity argument needs the source in rho_0 -- and the local Riccati
    # balance fixes the amplitude: with delta rho_0 = a F the order-(h kappa) residual is (2a - 1) sqrt(Lambda) F, zero only at a = 1/2.
    def prove_kodd_misplaced(terms):
        smax = max(s_ for _, s_, _, _ in terms)
        ll = {(0, 0, 0, 0, 0): P(R(-eps0) / 2), (1, 0, 0, 0, 0): P(R(n0) / 2)}
        V = {(1, 0, 1, 0, 0): P(-R(cX0) / 2), (1, 1, 0, 1, 0): P(R(kap0 * cX0) / 2), (1, 1, 0, 0, 0): P(sp.Rational(1, 4) + R(d20)), (0, 1, 0, 0, 0): P(R(c00))}
        def dz(k, r): return padd(pDw(r), pmul(r, ll), -k)
        rho = {0: pscale(ll, sp.Rational(-1, 2))}
        rho[1] = pscale(padd(padd(V, pmul(rho[0], rho[0]), -1), dz(0, rho[0]), -1), sp.Rational(1, 2))
        rho[1] = padd(rho[1], {(0, 0, 0, 0, 1): P(sp.Rational(1, 2))})       # MISPLACED: the source in rho_1
        for k in range(1, smax):
            acc = dz(k, rho[k])
            for m in range(0, k + 1): acc = padd(acc, pmul(rho[m], rho[k - m]))
            rho[k + 1] = pscale(acc, sp.Rational(-1, 2))
        s_ = smax; k = (s_ + 1) // 2; a0, b0 = -F(3, 4) * s_, -F(1, 2) * s_
        base, lin = {}, {}
        for (di, dj, a, b, f), v in rho[s_].items():
            fac = rf(R(a0) + ii, di) * rf(R(b0) + jj, dj) / rf(R(a0 + b0) + ii + jj, di + dj) if f else poch_rat(a0, di) * poch_rat(b0, dj) / poch_rat(a0 + b0, di + dj)
            (lin if f else base)[(a, b)] = (lin if f else base).get((a, b), 0) + v.as_expr() * fac
        T = base[(k, 0)]; Tlin = lin.get((k, 0), 0)
        expr = sp.together(sum(c * (lin.get((a, b), 0) - (base.get((a, b), 0) / T) * Tlin) / T for c, s__, a, b in terms))
        num, den = sp.fraction(expr); return sp.Poly(sp.expand(num), ii, jj).is_zero
    kodd_terms = [terms for name, kinds, where, terms in RELATIONS if kinds == ('KODD',)][0]
    ok_f = not prove_kodd_misplaced(kodd_terms)
    check("control (f): the kappa-odd source MISPLACED into rho_1 breaks the parity identity (the wiring, not the amplitude, is what the identity tests)", ok_f)
    # local balance: rho_0 -> rho_0 + a F for the source h kappa sqrt(Lambda) F; the coefficient of h kappa sqrt(Lambda) F in the Riccati residual is 2a - 1
    # Read the source actually produced by the same adapter used in prove().
    actual_rho0 = symbolic_rho('KODD', 3)[0]
    actual_linear = {key: value for key, value in actual_rho0.items() if key[4] == 1}
    local_residual = padd(pscale(actual_linear, 2), {(0, 0, 0, 0, 1): P(-1)})
    check("control (g): the ACTUAL kappa-odd adapter satisfies 2 delta rho_0 = F in the local Riccati balance", not local_residual)
    if failures: print(f"CONTROL NOT FIRED: {len(failures)}"); sys.exit(1)
    print(f"EXACT POINT DEFORMATION CONTROL FIRED ({time.time() - t2:.0f}s)"); sys.exit(0)
t3 = time.time()
for name, kinds, where, terms in RELATIONS:
    for kind in kinds:
        ok, num = prove(kind, terms)
        check(f"(2) class {name}: the relation over {where} is an IDENTITY in the symbolic exponents for the {kind}-shapes", ok)
    v = violation(terms, targets)
    check(f"(3) class {name}: with the exact layers set to zero the relation is VIOLATED by the targets (value {v} != 0)", v != 0)
print(f"   identities and violations ({time.time() - t3:.0f}s)", flush=True)
if failures:
    print(f"FAILS: {len(failures)}"); sys.exit(1)
print("EXACT POINT DEFORMATION CERTIFIED: at (N, s) = (28, -9) no first-order change of the sec.-38(ad) operator in the classes "
      "h W(w), h [X f + Y g], h [W + X f + Y g], h U/(kappa^2 Lambda), (h/kappa^2 Lambda)[X^2 g + XY g + Y^2 g], nor of its SYMBOL kappa^2 Lambda (1 + h sigma), nor a kappa^0 momentum-quadratic term, nor an analytic kappa-odd term "
      "(any monomial exponents; powers of Lambda absorbed) reproduces the certified first-order eigenvalue targets along the curve: eight proved integrand identities, each violated by the exact targets.")
