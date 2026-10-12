"""VIR5d addendum: conjecture P-W.  The PX <-> pi antisymmetric class of the certified kernel at the boundary
drop points t_i = (2i+3)/(2i+1), K = 2i+1, against  prod_{j=-i}^{i} (PX + pi - j p)(PX - pi - j p), p^2 = 1/(i+1).
T-a: i = 2 (t = 7/5), all coefficients.  T-b: i = 3 (t = 9/7) hold-out.  T-c: non-boundary drop points have no
antisymmetric class.  Mode tamper: p^2 -> 1/(i+1) + 1/1000 must fail.
Exit (real): 0 = T-a, T-b (if available), T-c hold; 2 otherwise.  Exit (tamper): 0 = fires."""
import hashlib, json, os, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_VIR5d_addendum.md', 'rb').read()).hexdigest() == open('SEAL_VIR5d_addendum.sha256').read().split()[0]
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
PX, PI = sp.symbols('PX PI')


def load(tstr, K):
    fn = 'z2_kernel_t%s_K%d.json' % (tstr.replace('/', '_'), K)
    if not os.path.exists(fn):
        return None
    kd = json.load(open(fn))
    return [{(int(k.split(',')[0]), int(k.split(',')[1])): sp.Rational(v) for k, v in c['eigenvalue'].items()} for c in kd['classes'] if c['in_HV']], kd['kernel_dimension']


def antisym(ker, K):
    mons = [(2 * a, 2 * b) for a in range(K + 1) for b in range(K + 1 - a)]
    M = sp.Matrix([[e_.get(m, 0) + e_.get((m[1], m[0]), 0) for e_ in ker] for m in mons])
    return [{m: sum(v[i] * ker[i].get(m, 0) for i in range(len(ker))) for m in mons} for v in M.nullspace()]


status = 0
for i, tag in ((1, 'training'), (2, 'T-a (4 of 15 coefficients were training)'), (3, 'T-b HOLD-OUT')):
    t = sp.Rational(2 * i + 3, 2 * i + 1); K = 2 * i + 1
    got = load(str(t), K)
    if got is None:
        print('i = %d (t = %s, K = %d): kernel file not available -- NOT TESTED' % (i, t, K)); continue
    ker, dim = got
    an = antisym(ker, K)
    p2 = sp.Rational(1, i + 1) + (sp.Rational(1, 1000) if mode == 'tamper' else 0)
    u, v = PX + PI, PX - PI
    prod = sp.Integer(1)
    for z in (u, v):
        prod *= z
        for j in range(1, i + 1):
            prod *= (z**2 - j * j * p2)
    Pp = sp.Poly(sp.expand(prod), PX, PI)
    pw = {m: c for m, c in Pp.terms()}
    cU = 13 - 6 * (i + 1) - sp.Rational(6, i + 1)
    pp = sp.Rational(2) * (t - 1) / (t + 1)
    cVir2 = -(2 * pp - 3) * (3 * pp - 2) / pp
    line = 'i = %d (t = %s, K = %d, weight %d) %s: kernel dimension %d; antisymmetric classes %d; c_U = c_W = %s (c_(1,%d) = %s)' % (i, t, K, 2 * K, tag, dim, len(an), cVir2, i + 1, cU)
    if dim != 2 or len(an) != 1:
        print(line + '  -> FAIL'); status = 2; continue
    E = an[0]
    lead = E[(2 * K, 0)]
    mons = set(m for m, c in E.items() if c != 0) | set(pw)
    bad = [m for m in mons if E.get(m, 0) / lead != pw.get(m, 0)]
    print(line + '; coefficients compared %d, mismatches %d' % (len(mons), len(bad)))
    if bad:
        status = 2
for tstr, K in (('11/5', 4), ('17/7', 5), ('13/7', 6), ('23/9', 6)):
    got = load(tstr, K)
    if got is None:
        print('T-c t = %s K = %d: kernel file not available -- NOT TESTED' % (tstr, K)); continue
    ker, dim = got
    an = antisym(ker, K)
    print('T-c t = %s K = %d (non-boundary drop point): kernel dimension %d, antisymmetric classes %d' % (tstr, K, dim, len(an)))
    if dim != 1 or an:
        status = 2
print('P-W (%s):' % mode, 'PASS' if status == 0 else 'FAIL')
sys.exit((0 if status == 0 else 2) if mode == 'real' else (0 if status == 2 else 1))
