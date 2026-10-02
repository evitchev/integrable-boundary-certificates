"""EXC1 oper side, step 1: WKB charges R_2, R_4, R_6 of the Suzuki + centrifugal equation with apparent singularities,
    V = x^(2M) + alpha x^(M-1) + (lam^2 - 1/4)/x^2 - 2 d_x^2 sum_k log(x^(M+1) - z_k),
as polynomials in (alpha, lam, z1, z2) [npts = 2] or (alpha, lam, z1) [npts = 1], via the large-xi series
    x^2 dV = 2(M+1) sum_k [ 1 + sum_{j>=1} (j(M+1)+1) (z_k/xi)^j ].
Usage: o1_engine.py <t> <npts> [tamper]   (tamper: the dV series scaled by 21/20, i.e. double-pole coefficient 2.1)
ODE side only.  Output: o1_t<..>_n<npts>.json with srepr of R_i."""
import hashlib, json, sys, time
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
import mellin_pot as MP
assert hashlib.sha256(open('SEAL_EXC1.md', 'rb').read()).hexdigest() == open('SEAL_EXC1.sha256').read().split()[0]
t = sp.Rational(sys.argv[1]); npts = int(sys.argv[2]); tamper = len(sys.argv) > 3
M = (t + 3) / (t - 1)
Mf = F(int(M.p), int(M.q))
K = 6
t0 = time.time()
h = [MP.ONE, MP.ZERO, MP.P_(-MP.l1**2)] + [MP.ZERO] * (K - 1)
eng = MP.MellinWKB(2, M, h)
scale = sp.Rational(21, 20) if tamper else sp.Integer(1)
V = {1: {(Mf - 1, F(-1)): MP.P_(MP.l0)}}
zs = [MP.z1, MP.z2][:npts]
for j in range(0, K - 1):
    pj = sum(z**j for z in zs) if j > 0 else sp.Integer(npts)
    coef = scale * 2 * (M + 1) * (j * (M + 1) + 1) * pj
    key = (-j * (Mf + 1) - 2, F(-1))
    V.setdefault(2 + j, {})
    V[2 + j] = MP.el_add(V[2 + j], {key: MP.P_(coef)})
eng.solve(K, V)
out = {'t': str(t), 'M': str(M), 'npts': npts, 'tamper': tamper, 'R': {}}
for i in (2, 4, 6):
    R, good = eng.R(i)
    assert good and eng.last_half == 0, (i, good, eng.last_half)
    out['R'][str(i)] = sp.srepr(sp.expand(R))
    print('t = %s M = %s npts = %d order %d: %d monomials  (%.0fs)' % (t, M, npts, i, len(sp.Poly(R, MP.l0, MP.l1, MP.z1, MP.z2).terms()), time.time() - t0), flush=True)
fn = 'o1_t%s_n%d%s.json' % (str(t).replace('/', '_'), npts, '_tamper' if tamper else '')
json.dump(out, open(fn, 'w'), indent=1)
print('written', fn)
