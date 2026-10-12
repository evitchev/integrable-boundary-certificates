"""RED1 (A), operator level, exact symbolic.
(0) the exact rescaling (R3) of the Sol-3 three-term equation, checked by acting on a monomial x^nu;
(1) limit operator L_inf: (th^2 - PX^2)(th^2 - pi^2) phi = xi (th + 1/2 - En) phi and the oscillator in the same form;
A1: terminating solutions of L_inf exactly at En = N + 1/2 + rho, rho in {+-PX, +-pi} (N = 0..3), and none at a generic En;
A2: the Frobenius recurrence c_m/c_(m-1) = (rho + m - 1/2 - En)/prod_r (rho + m - r): 1F3 with upper parameter rho + 1/2 - En;
    oscillator: 1F1 with the same upper parameter.
Mode tamper: drop the 1/2 (th - En): A1 must then fail.   Exit 0 = as sealed / tamper fires."""
import hashlib, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_RED1.md', 'rb').read()).hexdigest() == open('SEAL_RED1.sha256').read().split()[0]
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
PX, PI, En, nu, k, x, E = sp.symbols('PX PI En nu k x E')
ok = True

# (0) rescaling: act with the three-term operator on x^nu and with (R3) on yt^(nu/n)
n = k + 4; M = 1 / k; c = n * (1 + M); p2 = k / (k + 1); s = sp.Symbol('s')
l0 = n * PX / sp.sqrt(p2); l1 = n * PI / sp.sqrt(p2)
T = s - nu                                   # T on x^nu
lhs = (T**2 - l0**2) * (T**2 - l1**2)        # P(T) x^nu = lhs x^nu
mid = (s - (nu + c / 2))                     # x^(c/2) T x^(c/2) x^nu = x^(nu + c) (s - nu - c/2)
# equation on monomials:  lhs x^nu - mid x^(nu+c) - E x^(nu+n)
tau = T / n
lhsR = n**4 * (tau**2 - PX**2 / p2) * (tau**2 - PI**2 / p2)
midR = n * (tau - (1 + M) / 2)
chk0 = sp.simplify(lhs - lhsR) == 0 and sp.simplify(mid - midR) == 0
# scalings: yt = x^n / lam, lam = n^(3 p2): x^(nu+c) = x^nu yt^(1+M) lam^(1+M); divide by n^4: n lam^(1+M)/n^4 = 1  <=>  lam^(1/p2) = n^3
chk0 &= sp.simplify((1 + M) - 1 / p2) == 0
print('(0) rescaling (R3): symbol identities hold: %s  [P = n^4 (tau^2-PX^2/p^2)(tau^2-pi^2/p^2); middle = n x^c (tau - (1+M)/2); 1 + M = 1/p^2]' % chk0)
ok &= chk0
# coefficientwise limit k -> oo of (R3)
lim_ok = sp.limit(1 / p2, k, sp.oo) == 1 and sp.limit((1 + M) / 2, k, sp.oo) == sp.Rational(1, 2) and sp.limit(1 + M, k, sp.oo) == 1
print('    limit k -> oo: 1/p^2 -> 1, (1+M)/2 -> 1/2, exponent 1+M -> 1: %s' % lim_ok)
ok &= lim_ok

# (1) L_inf on xi^(rho + m): (th^2 - PX^2)(th^2 - PI^2) xi^a = P4(a) xi^a ;  xi (th + 1/2 - En) xi^a = (a + 1/2 - En) xi^(a+1)
half = sp.Rational(1, 2) if mode == 'real' else sp.Integer(0)
def P4(a):
    return (a**2 - PX**2) * (a**2 - PI**2)
def frob(rho, Nmax, EnVal):
    c = [sp.Integer(1)]
    for m in range(1, Nmax + 1):
        c.append(sp.cancel(c[-1] * (rho + m - 1 + half - EnVal) / P4(rho + m)))
    return c
a1 = True
for rho in (PX, -PX, PI, -PI):
    for N in range(0, 4):
        cs = frob(rho, N + 3, N + sp.Rational(1, 2) + rho)
        term = all(sp.simplify(cc) == 0 for cc in cs[N + 1:]) and sp.simplify(cs[N]) != 0
        a1 &= term
gen = frob(PX, 4, En)
a1_generic = all(sp.simplify(cc) != 0 for cc in gen)
print('A1: terminating solutions at En = N + 1/2 + rho for rho in {+-PX, +-pi}, N = 0..3: %s; generic En does not terminate: %s' % (a1, a1_generic))
# A2: hypergeometric identification
m = sp.Symbol('m', positive=True, integer=True)
ratio = (PX + m - 1 + half - En) / P4(PX + m)
# 1F3(a; b1,b2,b3; xi): c_m/c_(m-1) = (a + m - 1)/((b1+m-1)(b2+m-1)(b3+m-1) m), a = rho + 1/2 - En, b = 1 + rho - r' (r' the other three roots)
rho = PX
target = (rho + half - En + m - 1) / (m * (2 * rho + m) * (rho - PI + m) * (rho + PI + m))
a2 = sp.simplify(ratio - target) == 0
# oscillator: (th^2 - rho^2) chi = xi (th + 1/2 - En) chi : c_m/c_(m-1) = (rho + m - 1/2 - En)/(m (2 rho + m)) : 1F1(rho + 1/2 - En; 2 rho + 1; xi)
osc = (rho + m - 1 + half - En) / ((rho + m)**2 - rho**2)
a2 &= sp.simplify(osc - (rho + half - En + m - 1) / (m * (2 * rho + m))) == 0
print('A2: L_inf Frobenius solution = xi^rho 1F3(rho + 1/2 - En; 2rho+1, rho-pi+1, rho+pi+1; xi); oscillator = xi^rho 1F1(rho + 1/2 - En; 2rho+1; xi): %s' % a2)
# oscillator form from the Schrodinger equation: psi = y^(1/2) e^(-y^2/2) chi(xi = y^2)
# v1: done in the variable xi (d/dy = 2 sqrt(xi) d/dxi) with a generic function g(xi); v0 substituted derivatives of
# chi(y^2) by pattern and the pattern did not match sympy's Subs objects (a bug in the CHECK, kept as a_op_v0_osc_check_bug.py).
xiS, lam = sp.symbols('xi lam', positive=True)
g = sp.Function('g')(xiS)
psi = xiS ** sp.Rational(1, 4) * sp.exp(-xiS / 2) * g
ddy = lambda f_: 2 * sp.sqrt(xiS) * sp.diff(f_, xiS)
schr = -ddy(ddy(psi)) + (xiS + (lam**2 - sp.Rational(1, 4)) / xiS) * psi - E * psi
th = lambda f_: xiS * sp.diff(f_, xiS)
form = th(th(g)) - (lam / 2)**2 * g - xiS * (th(g) + g / 2 - E / 4 * g)
a_osc = sp.simplify(sp.expand(schr * xiS / (xiS ** sp.Rational(1, 4) * sp.exp(-xiS / 2)) + 4 * form)) == 0
print('    oscillator: y^2 x [Schrodinger operator on y^(1/2) e^(-y^2/2) chi(y^2)] == -4 x [(th^2 - (lam/2)^2) - xi (th + 1/2 - E/4)] chi: %s' % a_osc)
good = a1 and a1_generic and a2 and a_osc
if mode == 'real':
    print('(A) operator level:', 'AS SEALED' if (ok and good) else 'NOT as sealed')
    sys.exit(0 if (ok and good) else 2)
print('tamper (1/2 dropped): A1 holds: %s  -> %s' % (a1, 'fires' if not a1 else 'DOES NOT FIRE'))
sys.exit(0 if not a1 else 1)
