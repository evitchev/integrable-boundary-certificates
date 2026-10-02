"""SECOND1 extension (post-seal, labelled): shift structure of the tensor product of a Suzuki-type +-alpha DOUBLET
A_pm: theta^2 chi = u_pm theta chi + (lam_pm + x^a +- al x^ap - E x^b) chi."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
A, Ap, B, D, a, ap, b, E, al = sp.symbols('A Ap B D a ap b E alpha')
up, um, lp, lm = sp.symbols('u_p u_m lam_p lam_m')
def th(f):
    P = sp.Poly(sp.expand(f), A, Ap, B)
    return sp.expand(sum(c * (i * a + j * ap + l * b) * A**i * Ap**j * B**l for (i, j, l), c in P.terms()))
a0p = lp + A + al * Ap - E * B; a0m = lm + A - al * Ap - E * B
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
M = sp.Matrix(vecs).T
ns = M.nullspace(simplify=False)
print('nullspace dimension:', len(ns))
rel = [sp.factor(sp.cancel(r / ns[0][4])) for r in ns[0]]
den = sp.lcm([sp.denom(sp.together(r)) for r in rel])
print('denominator of the relation (non-polynomial coefficients if nontrivial):', sp.factor(den))
T = sp.expand(sum(sp.cancel(r * den) * D**k for k, r in enumerate(rel)))
P = sp.Poly(T, A, Ap, B)
print('x-monomials present (i a + j ap + l b):', sorted(P.monoms()))
