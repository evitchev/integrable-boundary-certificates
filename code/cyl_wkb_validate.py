"""WKB side of the cylindrical oper program (2026-09-02), two certified facts.

P1 -- the WKB machinery (code/wkb.py: Riccati recursion for LVZ-shaped ODEs
Psi'' = [kappa^2 (1+e^z)^m + V] Psi, Beta-function integration) is
validated on hep-th/0404195 sec. 5: for the paperclip ODE (62) in its
corrected form (the ((n+2)Q^2/2 + 1/4) term multiplies e^z/(1+e^z)^2, as
in LVZ (90)), S_1, S_3, S_5, S_7 satisfy S_s = alpha_s I_s + beta_s with
the paper's OWN VEVs (59)-(61) and residuals EXACTLY the Bernoulli
constants of eq. (70): 1/24, -1/2880, 1/40320, -1/215040.  S_1 and S_5
also match the printed (65), (67) verbatim; (66) and (68) each carry one
typo (the P^2Q^2 term of (66) lacks a factor (3n+4); the P^2Q^6 term of
(68) has (2+n)^2 for (2+n)^3), proven by the alpha_s I_s + beta_s
consistency.

P2 -- the MINIMAL LVZ shape is refuted for Solution 3: with V = A e^z/(1+e^z)
+ B e^z/(1+e^z)^2, A and B quadratic in (P, Q), palindromic top-degree
VEVs force C_5(ker l_1) = -32(m+1)/(5(5m+2)(5m+4)(5m+6)(5m+8)) = 0, i.e.
m = -1 on the whole conic; the certified u_4 = 6(t-3)/(t+1) then fixes
the map and PREDICTS u_6 = 3(5t-13)/(t+7), whereas the certified
(cyl_mo2_profiles) u_6 is 15(t-3)/(t+5): they agree only at t -> infinity
(the N = 1 paperclip point) and on the N = 0 fiber t = 5.  Runtime ~1 min."""
import sys
import sympy as sp
from wkb import u, m, wronskian_S

fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg)
    if not cond:
        fails.append(msg)

n, P, Q, A, B, x, y, p, q, t, r = sp.symbols('n P Q A B x y p q t r')
if __name__ == "__main__":
    V = n*P**2/2 * u/(1+u) + ((n+2)*Q**2/2 + sp.Rational(1, 4)) * u/(1+u)**2
    S = wronskian_S(V, 7, n)
    I = {1: -(P**2 + Q**2)/2 - sp.Rational(2, 24),
         3: (4+3*n)/(6*(2+n))*P**4 + (2+3*n)/(6*n)*Q**4 + P**2*Q**2 + (3+2*n)/(6*(2+n))*P**2 + (1+2*n)/(6*n)*Q**2 + (11+36*n+18*n**2)/(360*n*(2+n)),
         5: (-(6+5*n)*(8+5*n)/(120*(2+n)**2)*P**6 - (2+5*n)*(4+5*n)/(120*n**2)*Q**6 - (8+5*n)/(8*(2+n))*P**4*Q**2 - (2+5*n)/(8*n)*P**2*Q**4
             - (4+3*n)*(8+5*n)/(48*(2+n)**2)*P**4 - (2+3*n)*(2+5*n)/(48*n**2)*Q**4 - (2+10*n+5*n**2)/(8*n*(2+n))*P**2*Q**2
             - (30+225*n+250*n**2+76*n**3)/(480*n*(2+n)**2)*P**2 - (28+137*n+206*n**2+76*n**3)/(480*n**2*(2+n))*Q**2
             - (564+3410*n+7385*n**2+5680*n**3+1420*n**4)/(60480*n**2*(2+n)**2)),
         7: (n*(8+7*n)*(10+7*n)/(42*(2+n)**2*(2+7*n))*P**8 + (2+n)*(4+7*n)*(6+7*n)/(42*n**2*(12+7*n))*Q**8
             + 2*n*(10+7*n)/(3*(2+n)*(2+7*n))*P**6*Q**2 + 2*(2+n)*(4+7*n)/(3*n*(12+7*n))*P**2*Q**6 + P**4*Q**4
             + n*(5+4*n)*(10+7*n)/(9*(2+n)**2*(2+7*n))*P**6 + (2+n)*(3+4*n)*(4+7*n)/(9*n**2*(12+7*n))*Q**6
             + (6+49*n+28*n**2)/(3*(2+n)*(2+7*n))*P**4*Q**2 + (20+63*n+28*n**2)/(3*n*(12+7*n))*P**2*Q**4
             + (122+1407*n+1792*n**2+612*n**3)/(180*(2+n)**2*(2+7*n))*P**4 + (420+1583*n+1880*n**2+612*n**3)/(180*n**2*(12+7*n))*Q**4
             + 7*(40+338*n+985*n**2+816*n**3+204*n**4)/(30*n*(2+n)*(2+7*n)*(12+7*n))*P**2*Q**2
             + (2520+26166*n+102459*n**2+136612*n**3+74676*n**4+14552*n**5)/(1260*n*(2+n)**2*(2+7*n)*(12+7*n))*P**2
             + (3720+30202*n+89149*n**2+121284*n**3+70844*n**4+14552*n**5)/(1260*n**2*(2+n)*(2+7*n)*(12+7*n))*Q**2
             + (68760+632142*n+2264647*n**2+4095840*n**3+3708040*n**4+1610448*n**5+268408*n**6)/(151200*n**2*(2+n)**2*(2+7*n)*(12+7*n)))}
    beta = {1: sp.Rational(1, 24), 3: sp.Rational(-1, 2880), 5: sp.Rational(1, 40320), 7: sp.Rational(-1, 215040)}
    for k in (1, 3, 5, 7):
        topm = P**(k + 1)
        alpha = sp.simplify(sp.Poly(sp.expand(S[k]), P, Q).coeff_monomial(topm) / sp.Poly(sp.expand(I[k]), P, Q).coeff_monomial(topm))
        resid = sp.simplify(sp.expand(S[k] - alpha * I[k]))
        require(resid == beta[k], f"P1: S_{k} = alpha_{k} I_{k} + beta_{k} with the paper's I_{k} and beta_{k} = {beta[k]} (alpha_{k} = {sp.factor(alpha)})")
    art65 = P**2/2 + Q**2/2 + sp.Rational(1, 8)
    require(sp.simplify(S[1] - art65) == 0, "P1: S_1 equals the printed eq. (65)")
    d3 = sp.factor(sp.simplify(sp.expand(S[3] - (-n/(24*(2+3*n))*P**4 - (2+n)/(24*(4+3*n))*Q**4 - n*(2+n)/(4*(2+3*n))*P**2*Q**2
            - n*(3+2*n)/(24*(2+3*n)*(4+3*n))*P**2 - (2+n)*(1+2*n)/(24*(2+3*n)*(4+3*n))*Q**2 - (2+6*n+3*n**2)/(192*(2+3*n)*(4+3*n))))))
    require(d3 == sp.factor(3*P**2*Q**2*n*(n+1)*(n+2)/(4*(3*n+2)*(3*n+4))),
            "P1: S_3 differs from the printed eq. (66) ONLY in the P^2Q^2 term, by 3n(n+1)(n+2)P^2Q^2/(4(3n+2)(3n+4)) (the print lacks (3n+4))")
    art67 = (n**2/(40*(2+5*n)*(4+5*n))*P**6 + (2+n)**2/(40*(6+5*n)*(8+5*n))*Q**6
             + 3*n**2*(2+n)/(8*(2+5*n)*(4+5*n)*(6+5*n))*P**4*Q**2 + 3*n*(2+n)**2/(8*(4+5*n)*(6+5*n)*(8+5*n))*P**2*Q**4
             + n**2*(4+3*n)/(16*(2+5*n)*(4+5*n)*(6+5*n))*P**4 + (2+n)**2*(2+3*n)/(16*(4+5*n)*(6+5*n)*(8+5*n))*Q**4
             + 3*n*(2+n)*(2+10*n+5*n**2)/(8*(2+5*n)*(4+5*n)*(6+5*n)*(8+5*n))*P**2*Q**2
             + n*(30+225*n+250*n**2+76*n**3)/(160*(2+5*n)*(4+5*n)*(6+5*n)*(8+5*n))*P**2
             + (2+n)*(28+137*n+206*n**2+76*n**3)/(160*(2+5*n)*(4+5*n)*(6+5*n)*(8+5*n))*Q**2
             + (24+140*n+290*n**2+220*n**3+55*n**4)/(640*(2+5*n)*(4+5*n)*(6+5*n)*(8+5*n)))
    require(sp.simplify(sp.expand(S[5] - art67)) == 0, "P1: S_5 equals the printed eq. (67) verbatim")
    art68 = (-15*n**3/(448*(2+7*n)*(4+7*n)*(6+7*n))*P**8 - 15*(2+n)**3/(448*(8+7*n)*(10+7*n)*(12+7*n))*Q**8
             - 15*n**3*(2+n)/(16*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n))*P**6*Q**2
             - 15*n*(2+n)**2/(16*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*P**2*Q**6
             - 45*n**2*(2+n)**2/(32*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n))*P**4*Q**4
             - 5*n**3*(5+4*n)/(32*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n))*P**6
             - 5*(2+n)**3*(3+4*n)/(32*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*Q**6
             - 15*n**2*(2+n)*(6+49*n+28*n**2)/(32*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n))*P**4*Q**2
             - 15*n*(2+n)**2*(20+63*n+28*n**2)/(32*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*P**2*Q**4
             - n**2*(122+1407*n+1792*n**2+612*n**3)/(128*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n))*P**4
             - (2+n)**2*(420+1583*n+1880*n**2+612*n**3)/(128*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*Q**4
             - 21*n*(2+n)*(40+338*n+985*n**2+816*n**3+204*n**4)/(64*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*P**2*Q**2
             - n*(2520+26166*n+102459*n**2+136612*n**3+74676*n**4+14552*n**5)/(896*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*P**2
             - (2+n)*(3720+30202*n+89149*n**2+121284*n**3+70844*n**4+14552*n**5)/(896*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*Q**2
             - (12240+110628*n+386834*n**2+680568*n**3+606452*n**4+261786*n**5+43631*n**6)/(14336*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n)))
    d7 = sp.factor(sp.simplify(sp.expand(S[7] - art68)))
    require(d7 == sp.factor(-15*P**2*Q**6*n*(n+1)*(n+2)**2/(16*(7*n+6)*(7*n+8)*(7*n+10)*(7*n+12))),
            "P1: S_7 differs from the printed eq. (68) ONLY in the P^2Q^6 term, by -15n(n+1)(n+2)^2 P^2Q^6/(16(7n+6)(7n+8)(7n+10)(7n+12)) (the print has (2+n)^2 for (2+n)^3)")
    # P2: minimal LVZ shape for Solution 3
    Sg = wronskian_S(A*u/(1+u) + B*u/(1+u)**2, 5, m)
    top = lambda k, e: sum(cf*A**i*B**j for (i, j), cf in zip(sp.Poly(sp.expand(e), A, B).monoms(), sp.Poly(sp.expand(e), A, B).coeffs()) if i + j == (k+1)//2)
    Q3, C5 = top(3, Sg[3]), top(5, Sg[5])
    c5m = sp.factor(sp.simplify(C5.subs({A: m, B: -(m + 2)})))
    require(c5m == sp.factor(-32*(m+1)/(5*(5*m+2)*(5*m+4)*(5*m+6)*(5*m+8))),
            f"P2: C_5 on ker(l_1) = {c5m}: palindromic VEVs force m = -1")
    sub = {A: p*(x+y) + q*(x-y), B: -p*(x+y) + q*(x-y)}
    q3 = sp.Poly(sp.expand(Q3.subs(m, -1).subs(sub)), x, y); c5 = sp.Poly(sp.expand(C5.subs(m, -1).subs(sub)), x, y)
    u4 = sp.simplify(q3.coeff_monomial(x*y) / q3.coeff_monomial(x**2)).subs(q, sp.sqrt(r)*p)
    u6 = sp.simplify(c5.coeff_monomial(x**2*y) / c5.coeff_monomial(x**3)).subs(q, sp.sqrt(r)*p)
    rsol = sp.solve(sp.Eq(sp.simplify(u4), 6*(t-3)/(t+1)), r)
    require(len(rsol) == 1 and sp.simplify(rsol[0] - (t-5)/(t-2)) == 0, f"P2: the certified u_4 fixes r = q^2/p^2 = {rsol}")
    pred = sp.factor(sp.simplify(sp.simplify(u6).subs(r, rsol[0])))
    require(sp.simplify(pred - 3*(5*t-13)/(t+7)) == 0, f"P2: the minimal shape predicts u_6 = {pred}")
    diff6 = sp.factor(sp.simplify(pred - 15*(t-3)/(t+5)))
    require(diff6 == sp.factor(-24*(t-5)/((t+5)*(t+7))),
            f"P2: predicted minus certified u_6 = {diff6}: nonzero generically, vanishing only on the N = 0 fiber t = 5 and as t -> infinity (N = 1): minimal LVZ shape REFUTED")
    print()
    if fails:
        print(f"CYL WKB: {len(fails)} FAILURE(S)"); sys.exit(1)
    print("CYL WKB VALIDATED (paperclip S_1..S_7 = alpha_s I_s + Bernoulli beta_s with the paper's VEVs; two printed typos in (66), (68) "
          "identified; minimal LVZ shape forces m = -1 and mispredicts Solution 3's u_6 except on the N = 0 fiber and at N = 1 -- refuted)")
