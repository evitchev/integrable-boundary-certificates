"""The fusion-product formula for the extra charges' eigenvalues (ansatz note sec. 31): with q = 1/beta^2 even and
Lambda(P) = P^2 - Q^2 + (q-1)^2/(2q), the extra charge of spin 2q-3 at its first occurrence has eigenvalue
   e = - prod_{m = -(q-2)/2}^{(q-2)/2} Lambda(P + i sqrt2 m beta),
verified exactly against the certified-recipe eigenvalues at q = 4 (weight 6, t = 5/3) and q = 6 (weight 10, t = 7/5)
(results/lab/classical/extra_vev_out.txt); the q = 8 (weight 14, t = 9/7) eigenvalue is a prediction."""
import sympy as sp
P, Q = sp.symbols('P Q')
e6 = -(8*P**2 - 8*Q**2 + 9)*(64*P**4 - 128*P**2*Q**2 + 208*P**2 + 64*Q**4 - 80*Q**2 + 25)/512
e10 = -(12*P**2 - 12*Q**2 + 25)*(48*P**4 - 96*P**2*Q**2 + 232*P**2 + 48*Q**4 - 168*Q**2 + 147)*(48*P**4 - 96*P**2*Q**2 + 328*P**2 + 48*Q**4 - 72*Q**2 + 27)/27648
def fusion(q):
    beta = 1 / sp.sqrt(q); Lam = lambda p: p**2 - Q**2 + sp.Rational((q - 1)**2, 2 * q)
    return sp.expand(sp.prod([Lam(P + sp.I * sp.sqrt(2) * m * beta) for m in range(-(q - 2) // 2, (q - 2) // 2 + 1)]))
if __name__ == "__main__":
    for q, e in ((4, e6), (6, e10)):
        print(f"q = {q}: eigenvalue + fusion product == 0: {sp.expand(e + fusion(q)) == 0}")
    print("q = 8 prediction (weight 14, t = 9/7, sheet 3), e =", sp.expand(-fusion(8)))
