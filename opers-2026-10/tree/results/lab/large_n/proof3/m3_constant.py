"""PROOF3 M3: the constant check against NUM2's archived T1 data (t1.json, t = 2, k = a = 2, sigma = 3, gamma = 1/3, s = 5/2).
Theorem S predicts  (Q3_j / Q6_j) * Gamma(gamma - theta_j/sigma) * sigma^(-theta_j/sigma) = 1/kappa  for all four dynamical j, at every E (E-independent).
Q3 order: theta = s-l0, s+l0, s-l1, s+l1;  Q6 order: 1, 4, s-l0, s+l0, s-l1, s+l1 (ops.sixth).  Points A (l0, l1) = (1/3, 3/5), B (9/10, 1/5) (x1_control PTS).
Also: the exact E = 0 value from ops.exact_*_E0, and a TAMPER gamma -> gamma + 1/7 that must spread the four values."""
import json, mpmath as mp
mp.mp.dps = 60
from ops import exact_three_term_E0, exact_sixth_E0
d = json.load(open('t1.json'))
PTS = {'A': (mp.mpf(1) / 3, mp.mpf(3) / 5), 'B': (mp.mpf(9) / 10, mp.mpf(1) / 5)}
s, sig, gam = mp.mpf(5) / 2, mp.mpf(3), mp.mpf(1) / 3
def m_inv(theta, g=gam):   # 1/m(theta) = Gamma(g - theta/sigma) * sigma^(-theta/sigma)
    return mp.gamma(g - theta / sig) * sig ** (-theta / sig)
ok = True
for P, (l0, l1) in PTS.items():
    th = [s - l0, s + l0, s - l1, s + l1]
    ex3, _ = exact_three_term_E0(2, l0, l1); ex6, _ = exact_sixth_E0(l0, l1)
    v0 = [ex3[j] / ex6[j + 2] * m_inv(th[j]) for j in range(4)]
    spread0 = max(abs(v - v0[0]) / abs(v0[0]) for v in v0)
    print(f'point {P}: E = 0 exact: 1/kappa = {mp.nstr(v0[0], 30)}; spread over j {mp.nstr(spread0, 3)}')
    ok &= spread0 < mp.mpf(10) ** -50
    worst_spread, worst_E = mp.mpf(0), mp.mpf(0)
    for Es, rec in d[P].items():
        Q3 = [mp.mpf(x) for x in rec['Q3']]; Q6 = [mp.mpf(x) for x in rec['Q6']]
        v = [Q3[j] / Q6[j + 2] * m_inv(th[j]) for j in range(4)]
        sp_ = max(abs(x - v[0]) / abs(v[0]) for x in v); dE = abs(v[0] - v0[0]) / abs(v0[0])
        worst_spread = max(worst_spread, sp_); worst_E = max(worst_E, dE)
        vt = [Q3[j] / Q6[j + 2] * m_inv(th[j], gam + mp.mpf(1) / 7) for j in range(4)]
        spt = max(abs(x - vt[0]) / abs(vt[0]) for x in vt)
        print(f'  E = {Es:>7}: spread over j {mp.nstr(sp_, 3)}; |1/kappa(E) - 1/kappa(0)| rel {mp.nstr(dE, 3)}; TAMPER spread {mp.nstr(spt, 3)}')
        ok &= sp_ < mp.mpf(10) ** -50 and dE < mp.mpf(10) ** -50 and spt > mp.mpf(10) ** -3
    print(f'point {P}: worst j-spread {mp.nstr(worst_spread, 3)}, worst E-drift {mp.nstr(worst_E, 3)}')
print('M3', 'PASS' if ok else 'FAIL')
import sys; sys.exit(0 if ok else 1)
