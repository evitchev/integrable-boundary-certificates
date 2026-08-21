from ope import *
import sympy as sp
Fr1=Fr(1)
Qs = sp.Symbol('Q')
TQ = {mono(1,1): sp.Rational(1,2), mono(2): Qs}
zero_sym = lambda x: sp.simplify(x)==0
# I3(Q) = :T_Q T_Q: mod d
TTQ = nprod(TQ,TQ)
print(":T_Q T_Q: =", {k:sp.factor(v) for k,v in TTQ.items()})
cand = {mono(1,1,1,1):sp.Rational(1,4), mono(2,2):Qs**2-sp.Rational(1,2)}
diff = padd(TTQ,cand,-1)
print("difference is total derivative:", all(len(norm({k:sp.expand(v) for k,v in q.items()}))==0 for q in euler(diff,1).values()))
# I5(Q) = :T(:TT:): + A (dT)^2 ; find A from [I3(Q), I5(Q)] = 0
TTT = nprod(TQ, TTQ)
dT = deriv(TQ)
dTdT = nprod(dT,dT)
I3Q = cand
D0 = euler_vec(bracket_density(I3Q, TTT),1)
D1 = euler_vec(bracket_density(I3Q, dTdT),1)
keys = sorted(set(D0)|set(D1))
ratios=set()
for k in keys:
    a=sp.expand(D0.get(k,0)); b=sp.expand(D1.get(k,0))
    if b==0:
        assert a==0, (k,a)
    else:
        ratios.add(sp.factor(-a/b))
print("A determined by [I3,I5]=0:", ratios)
c = 1-12*Qs**2
print("BLZ (c+2)/12 =", sp.factor((c+2)/12))
A = list(ratios)[0]
I5Q = padd(TTT, dTdT, A)
# expand I5(Q) mod d in the basis J^6, J^2 J'^2, J''^2 : compute Euler image and match
def class_coords(P):
    """coordinates of the class of P (weight 6, even) in the basis J^6, J^2J'^2, J''^2 using Euler images"""
    basis = [{mono(1,1,1,1,1,1):1},{mono(1,1,2,2):1},{mono(3,3):1}]
    ev = [euler_vec(b,1) for b in basis]
    target = euler_vec(P,1)
    idx = sorted(set().union(*[set(e) for e in ev+[target]]))
    Mm = sp.Matrix([[sp.nsimplify(e.get(k,0)) for k in idx] for e in ev]).T
    tv = sp.Matrix([sp.nsimplify(target.get(k,0)) for k in idx])
    sol = sp.linsolve((Mm, tv))
    return sol
print("I5(Q) class coords (J^6, J^2J'^2, J''^2):", class_coords(I5Q))
# specialize Q^2 = 9/8
I5_98 = {k: sp.simplify(v.subs(Qs**2, sp.Rational(9,8)).subs(Qs, sp.sqrt(sp.Rational(9,8)))) for k,v in I5Q.items()}
print("I5 at Q^2=9/8 class coords:", class_coords(I5_98))
# paper-2 P6
P6 = {mono(3,3):Fr(31,24), mono(1,1,2,2):Fr(5,2), mono(1,1,1,1,1,1):Fr(1)}
print("paper2 P6 coords:", class_coords(P6))
# c=1: T0^3 + 41/16 (dT0)^2
T0 = {mono(1,1):Fr(1,2)}
TT0 = nprod(T0,T0); TTT0=nprod(T0,TT0); dT0=deriv(T0); dTdT0=nprod(dT0,dT0)
A2 = padd(TTT0, dTdT0, Fr(41,16))
print("A2(2) c=1 :T^3:+41/16 (dT)^2 coords:", class_coords(A2))
# KdV c=1: T0^3 + 1/4 (dT0)^2
K1 = padd(TTT0, dTdT0, Fr(1,4))
print("KdV c=1 I5 coords:", class_coords(K1))
# commutators
I3p = {mono(1,1,1,1):Fr1, mono(2,2):Fr(5,2)}   # enhanced I'_3
I3k = {mono(1,1,1,1):Fr1, mono(2,2):Fr(-2)}   # KdV c=1
print("[I'_3, P6] total derivative:", is_total_derivative(bracket_density(I3p,P6),1))
print("[I3^KdV(c=1), P6] total derivative:", is_total_derivative(bracket_density(I3k,P6),1))
print("[I'_3, K1] total derivative:", is_total_derivative(bracket_density(I3p,K1),1))
