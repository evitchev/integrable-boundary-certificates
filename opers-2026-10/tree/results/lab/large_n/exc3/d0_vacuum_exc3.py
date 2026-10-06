"""EXC3 T0(a) data side: VACUUM (level-0) eigenvalues of I_1, I_3, I_5 for Solution 2 at a fibre, from the record's tool
lab/excited_states.py (read-only) and the VIR1 kernel for I_5 -- level 0 ONLY (no excited block is touched; the I_5
level-one block is the lead's hold-out).  Output d0_sol2_t<t>.json with the three polynomials in bare (P, Q):
P = P_X, Q^2 = pi^2 - rho^2.   Usage: d0_vacuum_exc3.py <t>"""
import hashlib, json, os, sys, time
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
os.environ.setdefault('IB_CODE', '<repo>/code')
os.environ.setdefault('IB_LAB', '<repo>/lab')
sys.path.insert(0, os.environ['IB_LAB']); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import excited_states as XS
import vir_lib as V
assert hashlib.sha256(open('SEAL_EXC3.md', 'rb').read()).hexdigest() == open('SEAL_EXC3.sha256').read().split()[0]
sol = 2; t = F(sys.argv[1])
N, s_ = V.curve_point(t)
t0 = time.time()
charges = {'I3': V.P4(sol, N, s_)}
ker = V.Commutant(6).kernel(charges['I3'], N)
assert len(ker) == 1, 'kernel dimension %d' % len(ker)
B6 = V.basis(6)
charges['I5'] = {B6[i]: c for i, c in enumerate(ker[0]) if c}
out = {'sol': sol, 't': str(t), 'N': str(N), 'level': 0, 'vacuum': {}}
P, Q = XS.Pm, XS.Qm
st0 = XS.singlet_states(0)
for name, dens in charges.items():
    H0, G0 = XS.charge_matrix(dens, N, st0)
    vac = sp.expand(H0[0, 0] / G0[0, 0])
    out['vacuum'][name] = str(vac)
    print('Sol 2 t = %s %s vacuum eigenvalue: %s  (%.0fs)' % (t, name, vac, time.time() - t0), flush=True)
json.dump(out, open('d0_sol2_t%s.json' % str(t).replace('/', '_'), 'w'), indent=1, sort_keys=True)
print('written d0_sol2_t%s.json' % str(t).replace('/', '_'))
