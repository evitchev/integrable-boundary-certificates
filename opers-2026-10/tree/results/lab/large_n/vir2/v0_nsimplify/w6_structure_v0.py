"""VIR2 structure (post-seal, labelled): do the screenings DETERMINE the charges?
For each solution, the space of ALL two-boson densities of weight 4 and 6 (no parity or
Virasoro assumption) whose charges commute with the screening set, modulo total
derivatives, at two fibres -- and whether it is spanned by the cylindrical I_3, I_5.
Also: Sol 3's I_3 in the rotated bosons U = (X + phi)/sqrt2, W = (phi - X)/sqrt2.
Exit 0 iff all assertions hold."""
import sys, time
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
sys.path.insert(0, '<repo>/code')
import hashlib, json, pickle
import vir_lib as V
import vir2_lib as W
from vir2_lib import a, b, rho
t = sp.Symbol('t'); q = sp.Symbol('q'); p = sp.Symbol('p')
src = open('w4_main.py').read()
exec(src[src.index('I5SYM = '):src.index('# ---------------------------------------------------------------- fibres')])
t0 = time.time()
fails = 0


def fq(x):
    x = sp.nsimplify(x)
    assert x.is_Rational, x
    return F(int(x.p), int(x.q))


def commutant(w, screens, rho_val):
    """Densities of weight w (two-boson basis, all parities) conserved by every screening
    in `screens` = list of (qval, bval): V = exp(+-a X + bval phi), a^2 = qval (qval = 0: a = 0).
    Returns (dimension modulo total derivatives, kernel basis as Fraction rows, basis)."""
    B = W.two_boson_basis(w)
    rows_per_mono = []
    keys = {}
    for m in B:
        entries = {}
        for si, (qv, bv) in enumerate(screens):
            if qv == 0:
                nf = W.normal_form(W.screening_residue({m: sp.Integer(1)}, sp.Integer(0), bv), 'f', sp.Integer(0), bv)
            else:
                nf = W.normal_form(W.screening_residue({m: sp.Integer(1)}, a, bv), 'f', a, bv)
            for m2, c in nf.items():
                c = sp.expand(c.subs(rho, rho_val))
                pa = sp.Poly(c, a) if c.has(a) else None
                parts = {0: sp.Integer(0), 1: sp.Integer(0)}
                if pa is None:
                    parts[0] = c
                else:
                    for (k,), cf in pa.terms():
                        parts[k % 2] += cf * qv ** (k // 2)
                for par, val in parts.items():
                    if val != 0:
                        entries[(si, m2, par)] = fq(val)
        rows_per_mono.append(entries)
        for k in entries:
            keys.setdefault(k, len(keys))
    cols = []
    for entries in rows_per_mono:
        col = [F(0)] * len(keys)
        for k, v in entries.items():
            col[keys[k]] = v
        cols.append(col)
    ns = V.nullspace(cols, len(keys))
    # modulo total derivatives
    idx = {m: i for i, m in enumerate(B)}
    drows = []
    for m in W.two_boson_basis(w - 1):
        r = [F(0)] * len(B)
        for m2, c in W.d2({m: sp.Integer(1)}).items():
            r[idx[m2]] = fq(c)
        drows.append(r)
    rd = V.rank(drows)
    dim = V.rank(drows + ns) - rd
    return dim, ns, B, drows


def in_span_mod_d(vec, ns, drows):
    return V.rank(drows + ns + [vec]) == V.rank(drows + ns)


for pv in (sp.Rational(3, 2), sp.Rational(-2, 3)):
    fib = Fibre(pv)
    print('\n=== p = %s: t = %s, N = %s, rho = %s' % (pv, fib.t, fib.N, fib.rho), flush=True)
    for sol in (1, 2, 3):
        pairs = [(qv, bv) for _, bv, qv in PAIRS(sol, fib.p)]
        sets = [('pair(s) alone', pairs), ('pair(s) + Virasoro beta_+', pairs + [(0, fib.p)]),
                ('Virasoro beta_+ alone', [(0, fib.p)])]
        if sol == 3:
            sets.insert(1, ('S3 alone', pairs[:1]))
        for w in (4, 6):
            dens = fib.I3(sol) if w == 4 else fib.I5(sol)
            B = W.two_boson_basis(w)
            idx = {m: i for i, m in enumerate(B)}
            dvec = [F(0)] * len(B)
            for m, c in dens.items():
                dvec[idx[m]] = fq(c)
            for name, sc in sets:
                dim, ns, _, drows = commutant(w, sc, fib.rho)
                inside = in_span_mod_d(dvec, ns, drows)
                print('   Sol %d weight %d, %-26s commutant dim mod d = %2d ; cylindrical I_%d inside: %s  (%.0fs)'
                      % (sol, w, name + ':', dim, w - 1, inside, time.time() - t0), flush=True)
                if 'alone' not in name or name.startswith('pair'):
                    fails += not inside
                if name == 'pair(s) + Virasoro beta_+':
                    fails += dim != 1
print()
# ---- Sol 3 in rotated bosons
print('== Sol 3: I_3 in U = (X + phi)/sqrt2, W = (phi - X)/sqrt2, symbolic in p ==')
tp = (2 + p**2) / (2 - p**2); rp = p / 2 - 1 / p
Np = sp.cancel(1 - 12 * rp**2); sP = sp.cancel(3 * p**2 - 12 / p**2)
A2 = sp.cancel((13 - 27 * Np + 2 * Np**2 + (5 - 2 * Np) * sP) / 6); Cx = sp.cancel(7 - Np + sP)
B4 = V.basis(4); i4 = {m: i for i, m in enumerate(B4)}
vec = [sp.Integer(0)] * len(B4)
vec[i4[((), ((1, 1), (1, 1)))]] = 1; vec[i4[((), ((2, 2),))]] = -2
vec[i4[((2, 2), ())]] = A2; vec[i4[((1, 1), ((1, 1),))]] = Cx; vec[i4[((1, 1, 1, 1), ())]] = 1
d3 = {m: sp.cancel(c.subs(rho, rp)) for m, c in W.to_two_boson(vec, 4)[0].items()}
rt2 = sp.sqrt(2)
def rotate(d):
    out = {}
    for (xs, fs), c in d.items():
        terms = [((), (), c)]
        for k in xs:       # X = (U - W)/sqrt2
            terms = [(us + (k,), ws, cc / rt2) for us, ws, cc in terms] + [(us, ws + (k,), -cc / rt2) for us, ws, cc in terms]
        for k in fs:       # phi = (U + W)/sqrt2
            terms = [(us + (k,), ws, cc / rt2) for us, ws, cc in terms] + [(us, ws + (k,), cc / rt2) for us, ws, cc in terms]
        for us, ws, cc in terms:
            W.add_to(out, W.canon2(us, ws), cc)
    return {m: sp.cancel(sp.expand(c)) for m, c in out.items() if sp.cancel(sp.expand(c)) != 0}
d3r = rotate(d3)
rU = (1 / p - p) / rt2
TT = {m: c.subs(rho, rU) for m, c in W.vir_phi((2, 2)).items()}          # (T T) of one boson with rho_U
T1 = {m: c.subs(rho, rU) for m, c in W.vir_phi((2,)).items()}
TU2 = {W.canon2(fs, ()): c for (_, fs), c in TT.items()}
TW2 = {W.canon2((), fs): c for (_, fs), c in TT.items()}
TUW = {}
for (_, f1), c1 in T1.items():
    for (_, f2), c2 in T1.items():
        W.add_to(TUW, W.canon2(f1, f2), c1 * c2)
al, la = sp.symbols('alpha lambda')
target = dict(d3r)
for m, c in TU2.items(): W.add_to(target, m, -al * c)
for m, c in TW2.items(): W.add_to(target, m, -al * c)
for m, c in TUW.items(): W.add_to(target, m, -la * c)
# solve for alpha, lambda and a total derivative
Bw = W.two_boson_basis(4); idx = {m: i for i, m in enumerate(Bw)}
Bm = W.two_boson_basis(3)
ks = sp.symbols('k0:%d' % len(Bm))
eqs = [sp.Integer(0)] * len(Bw)
for m, c in target.items():
    eqs[idx[m]] += c
for kj, m in zip(ks, Bm):
    for m2, c in W.d2({m: sp.Integer(1)}).items():
        eqs[idx[m2]] -= kj * c
solv = sp.solve(eqs, [al, la] + list(ks), dict=True)
ok = len(solv) == 1
print('   I_3 (Sol 3) = alpha [(T_U T_U) + (T_W T_W)] + lambda T_U T_W + total derivative, with rho_U = rho_W = (1/p - p)/sqrt2:', ok)
if ok:
    print('   alpha =', sp.factor(solv[0][al]), '  lambda =', sp.factor(solv[0][la]),
          '  c_U = c_W = 1 - 12 rho_U^2 =', sp.factor(sp.cancel(1 - 12 * rU**2)))
fails += not ok
print('\nSTRUCTURE', 'FAILED (%d)' % fails if fails else 'OK')
sys.exit(1 if fails else 0)
