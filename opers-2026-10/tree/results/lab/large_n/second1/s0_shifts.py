"""SECOND1 pre-seal structural computation (no comparison with Sol 3 data): the tensor (Hadamard) product of two second-order theta-equations
A_pm: theta^2 chi = u_pm theta chi + (lam_pm + x^a - E x^b) chi   (a, b, u, lam, E symbolic).  Which x-monomials x^(i a + j b) occur, with which theta-polynomials?"""
import sys; sys.dont_write_bytecode = True
import sympy as sp
A, B, D, a, b, E = sp.symbols('A B D a b E')
up, um, lp, lm = sp.symbols('u_p u_m lam_p lam_m')
def th(f):  # theta acting on a polynomial in A = x^a, B = x^b
    f = sp.expand(f); P = sp.Poly(f, A, B)
    return sp.expand(sum(c * (i * a + j * b) * A**i * B**j for (i, j), c in P.terms()))
a0p = lp + A - E * B; a0m = lm + A - E * B
def theta_vec(v):
    c1, c2, c3, c4 = v
    out = [th(c1), th(c2), th(c3), th(c4)]
    out[1] += c1; out[2] += c1
    out[0] += c2 * a0p; out[1] += c2 * up; out[3] += c2
    out[0] += c3 * a0m; out[2] += c3 * um; out[3] += c3
    out[2] += c4 * a0p; out[1] += c4 * a0m; out[3] += c4 * (up + um)
    return [sp.expand(o) for o in out]
vecs = [[sp.Integer(1), 0, 0, 0]]
for _ in range(4): vecs.append(theta_vec(vecs[-1]))
ns = sp.Matrix(vecs).T.nullspace(); assert len(ns) == 1
rel = [sp.factor(r / ns[0][4]) for r in ns[0]]
den = sp.lcm([sp.denom(sp.together(r)) for r in rel])
print('common denominator of the relation:', sp.factor(den))
T = sp.expand(sum(sp.cancel(r * den) * D**k for k, r in enumerate(rel)))
P = sp.Poly(T, A, B)
for (i, j), c in sorted(P.terms()):
    print(f'  x^({i} a + {j} b):  {sp.factor(sp.Poly(c, D).as_expr())}')
