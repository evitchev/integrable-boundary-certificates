"""Step 1 (v2: exact Rationals in the weight checks; v1 s1_frame.py exit 1 on a float artefact, not shipped): the screening frame at t = 9/7 (q = 8, i = 3, m = 4) RE-DERIVED from the VIR2 record.
Source formulas (results/lab/large_n/vir2/NOTE.md, verdict table, sec. 2, sec. 3 'Solution 3 in rotated bosons', sec. 5;
SEAL_VIR2.md sec. 2), with VIR2's machine parameter p = beta_+ = sqrt2 u, u^2 = (t-1)/(t+1):
  p^2 = 2(t-1)/(t+1), rho = p/2 - 1/p, N = 1 - 12 rho^2, T_phi = (1/2)(dphi)^2 + rho d^2phi, h(e^{b phi}) = b^2/2 - rho b;
  Sol 3 pairs: S3  (beta, a^2) = (1/p, 1/p^2);  S3' (beta, a^2) = (-p, p^2);  Virasoro beta_+ = p, beta_- = -2/p (a = 0);
  rotated bosons U, W = (phi +- X)/sqrt2 with rho_U = rho_W = (1/p - p)/sqrt2, c_U = c_W = -(2p^2-3)(3p^2-2)/p^2.
Conservation of all four Sol-3 screenings against I_3 and I_5 identically in p: vir2/w3_symbolic.log (so at p = 1/2 too).
This script: exact sympy algebra only (no kernel, no OPE). Exit 0 iff every check holds."""
import hashlib, sys
import sympy as sp
print('s1_frame2.py sha256', hashlib.sha256(open(__file__, 'rb').read()).hexdigest())
R = sp.Rational
X, phi, dV, dZ, PX, PI = sp.symbols('X phi dV dZ P_X pi')
checks = []
def chk(name, ok):
    ok = bool(ok); checks.append(ok); print(('OK   ' if ok else 'FAIL ') + name)

def frame(t):
    """VIR2 data at t -> the (V, Z) frame: V = (X + phi)/(2p), Z = (phi - X)/(2p)."""
    p2 = 2 * (t - 1) / (t + 1); p = sp.sqrt(p2)          # VIR2: p = sqrt2 u, u^2 = (t-1)/(t+1), u > 0 branch
    rho = p / 2 - 1 / p; N = sp.expand(1 - 12 * rho**2)
    kappa = 1 / (2 * p2)                                  # <V,V>: dV = (dX + dphi)/(2p), <X,X> = <phi,phi> = 1
    # screening exponents (a, beta) in (X, phi) -> (alpha_V, alpha_Z): aX + beta phi = ((a+beta)/2p) ... see below
    # X = p (V - Z), phi = p (V + Z)  =>  aX + beta phi = p (a + beta) V + p (beta - a) Z
    def vz(a, beta): return sp.simplify(p * (a + beta)), sp.simplify(p * (beta - a))
    S3p, S3m = vz(1 / p, 1 / p), vz(-1 / p, 1 / p)        # S3: beta = 1/p, a = +-1/p
    S3pp, S3pm = vz(p, -p), vz(-p, -p)                    # S3': beta = -p, a = +-p
    Vp, Vm = vz(0, p), vz(0, -2 / p)                      # Virasoro beta_+ = p, beta_- = -2/p
    rhoU = (1 / p - p) / sp.sqrt(2)
    cU = sp.simplify(1 - 12 * rhoU**2)
    # T_U = (1/2)(dU)^2 + rho_U d^2U with U = sqrt2 p V  (U = (phi+X)/sqrt2 = sqrt2 p V):  (dU)^2 = 2 p^2 (dV)^2
    TV_quad = sp.simplify(R(1, 2) * 2 * p2); TV_lin = sp.simplify(rhoU * sp.sqrt(2) * p)   # T_V = TV_quad (dV)^2 + TV_lin d^2V
    # T_total = T_X + T_phi = (1/2)(dX)^2 + (1/2)(dphi)^2 + rho d^2phi, X = p(V - Z), phi = p(V + Z)
    Ttot_quad = sp.simplify(R(1, 2) * (p2 + p2)); Ttot_lin = sp.simplify(rho * p)        # per boson
    PV = (PX + PI) / (2 * p); PZ = (PI - PX) / (2 * p)
    return dict(p=p, p2=p2, rho=rho, N=N, kappa=kappa, S3=(S3p, S3m), S3p=(S3pp, S3pm), Vir=(Vp, Vm), rhoU=rhoU, cU=cU,
                TV=(TV_quad, TV_lin), Ttot=(Ttot_quad, Ttot_lin), PV=sp.simplify(PV), PZ=sp.simplify(PZ))

def h_exp(beta, kappa, quad, lin):
    beta = sp.Rational(beta); kappa = sp.Rational(kappa)
    """weight of e^{beta V} under T = quad (dV)^2 + lin d^2V with <V,V> = kappa: kappa beta^2 quad*... (quad = 1/(2 kappa)): h = kappa beta^2/2 - lin kappa beta"""
    assert sp.simplify(quad - 1 / (2 * kappa)) == 0
    return sp.simplify(kappa * beta**2 / 2 - lin * kappa * beta)

def c_of(kappa, quad, lin):
    kappa = sp.sympify(kappa)
    assert sp.simplify(quad - 1 / (2 * kappa)) == 0
    return sp.simplify(1 - 12 * lin**2 * kappa)

# ---------------- t = 9/7
t = R(9, 7); f = frame(t)
print('\n== t = 9/7 ==')
print('p^2 =', f['p2'], ' p =', f['p'], ' rho =', sp.simplify(f['rho']), ' N =', f['N'], ' kappa = <V,V> =', f['kappa'])
chk('p^2 = 1/4, p = 1/2 (rational: V = X + phi, Z = phi - X exactly)', f['p2'] == R(1, 4) and f['p'] == R(1, 2))
chk('rho = -7/4, N = -143/4 = 13 - 6(beta^2 + beta^-2) with beta^2 = (t-1)/(t+1) = 1/8', sp.simplify(f['rho'] + R(7, 4)) == 0 and f['N'] == R(-143, 4) and 13 - 6 * (R(1, 8) + 8) == R(-143, 4))
u2 = (t - 1) / (t + 1); s_ = 6 * (u2 - 1 / u2)
chk('curve: s^2 = (N - 25)(N - 1) at N = -143/4, s = -189/4', s_ == R(-189, 4) and sp.simplify(s_**2 - (f['N'] - 25) * (f['N'] - 1)) == 0)
chk('kappa = 2', f['kappa'] == 2)
S3p, S3m = f['S3']; S3pp, S3pm = f['S3p']; Vp, Vm = f['Vir']
print('S3  (beta = 1/p = 2, a = +-2):  e^{2X + 2phi} -> (alpha_V, alpha_Z) =', S3p, ';  e^{-2X + 2phi} ->', S3m)
print("S3' (beta = -p = -1/2, a = +-1/2): e^{X/2 - phi/2} ->", S3pp, ';  e^{-X/2 - phi/2} ->', S3pm)
print('Virasoro beta_+ = 1/2 ->', Vp, ';  beta_- = -4 ->', Vm)
chk('pair A = {e^{2V}, e^{2Z}} (S3)', set([S3p, S3m]) == {(2, 0), (0, 2)})
chk("pair B = {e^{-V/2}, e^{-Z/2}} (S3')", set([S3pp, S3pm]) == {(0, R(-1, 2)), (R(-1, 2), 0)})
chk('Virasoro = {e^{(V+Z)/4}, e^{-2(V+Z)}}', Vp == (R(1, 4), R(1, 4)) and Vm == (-2, -2))
q, l = f['TV']; qt, lt = f['Ttot']
print('T_V = %s (dV)^2 + %s d^2V   (from T_U = (1/2)(dU)^2 + rho_U d^2U, rho_U = %s, U = sqrt2 p V)' % (q, l, sp.simplify(f['rhoU'])))
print('T_total per boson = %s (dV)^2 + %s d^2V ;  T_V + T_Z - T_total = %s (d^2V + d^2Z)' % (qt, lt, sp.simplify(l - lt)))
chk('T_V = (1/4)(dV)^2 + (3/4) d^2V', q == R(1, 4) and l == R(3, 4))
chk('c_V = -25/2 = c_(1,4) = VIR2 c_U = -(2p^2-3)(3p^2-2)/p^2 at p^2 = 1/4', c_of(2, q, l) == R(-25, 2) and f['cU'] == R(-25, 2) and 1 - 6 * R(9, 4) == R(-25, 2))
chk('T_total = (1/4)((dV)^2+(dZ)^2) - (7/8)(d^2V + d^2Z); c_total = 1 + N = -139/4', qt == R(1, 4) and lt == R(-7, 8) and 2 * c_of(2, qt, lt) == 1 + f['N'])
chk('T_V + T_Z - T_total = (13/8)(d^2V + d^2Z)', sp.simplify(l - lt) == R(13, 8))
hs = {b: h_exp(b, R(2), q, l) for b in (R(2), R(-1, 2), R(-2), R(1, 4))}
print('weights under T_V: e^{2V}:', hs[2], ' e^{-V/2}:', hs[R(-1, 2)], ' e^{-2V}:', hs[-2], ' e^{V/4}:', hs[R(1, 4)])
chk('e^{2V}, e^{-V/2} weight 1 under T_V (the c_(1,4) Coulomb gas: alpha_+ = 2, alpha_- = -1/2, alpha_+ alpha_- kappa = -2)', hs[2] == 1 and hs[R(-1, 2)] == 1)
chk('e^{-2V} weight 7 = 2p_Ad - 1 (p_Ad = 4): H^(7) = Res e^{2V} e^{-2V} has weight 7', hs[-2] == 7)
ht = {b: h_exp(b, R(2), qt, lt) for b in (R(1, 4), R(-2))}
chk('Virasoro screenings weight 1 under T_total (1/2 per boson): e^{(V+Z)/4}, e^{-2(V+Z)}', ht[R(1, 4)] == R(1, 2) and ht[-2] == R(1, 2))
print('momenta: P_V =', f['PV'], ' P_Z =', f['PZ'])
chk('P_V = P_X + pi, P_Z = pi - P_X', sp.simplify(f['PV'] - (PX + PI)) == 0 and sp.simplify(f['PZ'] - (PI - PX)) == 0)
# the vacuum-eigenvalue predictions in this frame
P = sp.Symbol('P')
binom7 = sp.expand(sp.binomial(2 * P + 3, 7).expand(func=True)) if False else sp.expand(sp.prod([(2 * P + 3 - k) for k in range(7)]) / sp.factorial(7))
prod7 = sp.expand(sp.prod([(P - sp.Rational(j, 2)) for j in range(-3, 4)]))
chk('binom(2P+3, 7) = (8/315) P (P^2-1/4)(P^2-1)(P^2-9/4) = (128/5040) prod_{j=-3}^{3} (P - j/2)', sp.expand(binom7 - R(8, 315) * P * (P**2 - R(1, 4)) * (P**2 - 1) * (P**2 - R(9, 4))) == 0 and sp.expand(binom7 - R(128, 5040) * prod7) == 0)
E3 = sp.expand(sp.prod([(PX + PI - R(j, 2)) * (PX - PI - R(j, 2)) for j in range(-3, 4)]))
bb = sp.expand(binom7.subs(P, PX + PI) * binom7.subs(P, PI - PX))
chk('binom(2P_V+3,7) binom(2P_Z+3,7) = -(64/99225) E_3 (E_3 = prod_j (P_X+pi-j/2)(P_X-pi-j/2), p = 1/2)', sp.expand(bb + R(64, 99225) * E3) == 0)
chk('E_3 antisymmetric under P_X <-> pi', sp.expand(E3.subs({PX: PI, PI: PX}, simultaneous=True) + E3) == 0)
# ---------------- comparison with the design of the test (B5(iv))
design = dict(pairA={(2, 0), (0, 2)}, pairB={(R(-1, 2), 0), (0, R(-1, 2))}, Vir=((R(1, 4), R(1, 4)), (-2, -2)), TV=(R(1, 4), R(3, 4)), c=R(-25, 2), PV=PX + PI, const=R(-64, 99225))
same = (set(f['S3']) == design['pairA'] and set(f['S3p']) == design['pairB'] and f['Vir'] == design['Vir'] and f['TV'] == design['TV'] and c_of(2, q, l) == design['c'] and sp.simplify(f['PV'] - design['PV']) == 0)
chk('DERIVED FRAME == DESIGN FRAME of B5(iv) (no difference to record)', same)
# ---------------- the same pipeline at the records' points (reproduction control)
print('\n== reproduction of the records at t = 5/3 (WUW1) and t = 7/5 (WUW2) by the same pipeline ==')
g = frame(R(5, 3))
chk('t = 5/3: p^2 = 1/2, kappa = 1 (V = U): pair A e^{2U}, pair B e^{-U}, Virasoro e^{(U+W)/2}, e^{-2(U+W)}; T_U = (1/2)(dU)^2 + (1/2)d^2U, c = -2',
    g['p2'] == R(1, 2) and g['kappa'] == 1 and set(g['S3']) == {(2, 0), (0, 2)} and set(g['S3p']) == {(-1, 0), (0, -1)} and g['Vir'] == ((R(1, 2), R(1, 2)), (-2, -2)) and g['TV'] == (R(1, 2), R(1, 2)) and g['cU'] == -2)
g = frame(R(7, 5))
chk('t = 7/5: p^2 = 1/3, kappa = 3/2: pair A e^{2V}, pair B e^{-(2/3)V}, Virasoro e^{(V+Z)/3}, e^{-2(V+Z)}; T_U = (1/3)(dV)^2 + (2/3)d^2V, c = -7; T_total - (T_U+T_W) = -(3/2)(d^2V+d^2Z); P_V = (sqrt3/2)(P_X+pi)',
    g['p2'] == R(1, 3) and g['kappa'] == R(3, 2) and set(g['S3']) == {(2, 0), (0, 2)} and set(g['S3p']) == {(R(-2, 3), 0), (0, R(-2, 3))} and g['Vir'] == ((R(1, 3), R(1, 3)), (-2, -2)) and g['TV'] == (R(1, 3), R(2, 3)) and g['cU'] == -7 and sp.simplify(g['Ttot'][1] - g['TV'][1] + R(3, 2)) == 0 and sp.simplify(g['PV'] - sp.sqrt(3) / 2 * (PX + PI)) == 0)
print('\n== general even q = 2m (p^2 = 1/m, kappa = m/2): pair A e^{2V}, pair B e^{-(2/m)V}, Virasoro e^{(V+Z)/m}, e^{-2(V+Z)}, T_V = (1/m)(dV)^2 + ((m-1)/m) d^2V, c = c_(1,m), P_V = (sqrt m/2)(P_X+pi) ==')
m = sp.Symbol('m', positive=True); gg = frame(sp.Rational(1) * (2 * m + 1) / (2 * m - 1))   # t_i = (2i+3)/(2i+1), m = i+1
chk('general m: p^2 = 1/m, kappa = m/2, pair B exponent -2/m, Virasoro (1/m, 1/m) and (-2,-2), T_V = (1/m)(dV)^2 + ((m-1)/m) d^2V, c_V = 1 - 6(m-1)^2/m',
    sp.simplify(gg['p2'] - 1 / m) == 0 and sp.simplify(gg['kappa'] - m / 2) == 0 and sp.simplify(gg['S3p'][1][0] + 2 / m) == 0 and sp.simplify(gg['Vir'][0][0] - 1 / m) == 0 and gg['Vir'][1] == (-2, -2)
    and sp.simplify(gg['TV'][0] - 1 / m) == 0 and sp.simplify(gg['TV'][1] - (m - 1) / m) == 0 and sp.simplify(c_of(m / 2, gg['TV'][0], gg['TV'][1]) - (1 - 6 * (m - 1)**2 / m)) == 0)
print('\nALL CHECKS PASS' if all(checks) else '\nSOME CHECK FAILED')
sys.exit(0 if all(checks) else 1)
