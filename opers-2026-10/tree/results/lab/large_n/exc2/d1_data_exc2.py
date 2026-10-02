"""EXC1 D1: level-1 blocks of I_3 (and optionally I_5) with the record's tool lab/excited_states.py (imported read-only).
Usage: d1_data.py <sol> <t> [with5]
Output JSON: 2x2 singlet matrices (basis: X oscillator, e-oscillator), transverse 1x1, in bare (P, Q); characteristic polynomials.
Variables: bare P = PX, bare Q^2 = pi^2 - rho^2."""
import hashlib, json, os, sys, time
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
os.environ.setdefault('IB_CODE', '<repo>/code')
os.environ.setdefault('IB_LAB', '<repo>/lab')
sys.path.insert(0, os.environ['IB_LAB']); sys.path.insert(0, os.environ['IB_CODE'])
import excited_states as XS
import vir_lib as V
assert hashlib.sha256(open('SEAL_EXC2_addendum.md', 'rb').read()).hexdigest() == open('SEAL_EXC2_addendum.sha256').read().split()[0]
sol = int(sys.argv[1]); t = F(sys.argv[2]); with5 = len(sys.argv) > 3
N, s_ = V.curve_point(t)
t0 = time.time()
charges = {'I3': V.P4(sol, N, s_)}
if with5:
    ker = V.Commutant(6).kernel(charges['I3'], N)
    assert len(ker) == 1, 'kernel dimension %d' % len(ker)
    B6 = V.basis(6)
    charges['I5'] = {B6[i]: c for i, c in enumerate(ker[0]) if c}
out = {'sol': sol, 't': str(t), 'N': str(N), 'blocks': {}}
P, Q = XS.Pm, XS.Qm
for name, dens in charges.items():
    # level 0 check against the vacuum polynomial
    H0, G0 = XS.charge_matrix(dens, N, XS.singlet_states(0))
    vac = sp.expand(H0[0, 0] / G0[0, 0])
    st = XS.singlet_states(1)
    H, Gm = XS.charge_matrix(dens, N, st)
    Mx = sp.simplify(Gm.inv() * H)
    tst = XS.transverse_states(1)
    Ht, Gt = XS.charge_matrix(dens, N, tst)
    tr = sp.expand(Ht[0, 0] / Gt[0, 0])
    shift = sp.expand(vac.subs(Q, sp.sqrt(Q**2 + 2)))        # bare Q^2 -> Q^2 + 2  (article Q^2 -> Q^2 - 2)
    lam = sp.Symbol('lam')
    cp = sp.factor(sp.expand((Mx - lam * sp.eye(2)).det()))
    out['blocks'][name] = {'states': str(st), 'vacuum': str(vac), 'matrix': [[str(sp.expand(Mx[i, j])) for j in range(2)] for i in range(2)],
                           'transverse': str(tr), 'transverse_is_level_shift': bool(sp.expand(tr - shift) == 0), 'charpoly': str(sp.expand((Mx - lam * sp.eye(2)).det()))}
    print('Sol %d t = %s %s: level-1 singlet states %s; transverse eigenvalue == vacuum at Q^2 + 2: %s; charpoly irreducible over Q: %s  (%.0fs)'
          % (sol, t, name, st, out['blocks'][name]['transverse_is_level_shift'], sp.Poly(sp.expand((Mx - lam * sp.eye(2)).det()), lam).is_irreducible, time.time() - t0), flush=True)
    out['blocks'][name]['_M'] = None
    globals()['M_' + name] = Mx
if with5:
    comm = sp.simplify(M_I3 * M_I5 - M_I5 * M_I3)
    out['commute'] = bool(comm == sp.zeros(2, 2))
    print('   [I_3, I_5] = 0 on the singlet block: %s' % out['commute'])
for b in out['blocks'].values():
    b.pop('_M')
tag = 'sol%d_t%s%s' % (sol, str(t).replace('/', '_'), '_with5' if with5 else '')
json.dump(out, open('d1_%s.json' % tag, 'w'), indent=1, sort_keys=True)
print('written d1_%s.json' % tag)
