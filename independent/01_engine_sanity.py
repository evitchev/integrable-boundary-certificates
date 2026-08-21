from ope import *
import sympy as sp
Fr1 = Fr(1)
J = lambda a: ((0,a),)
T0 = {mono(1,1): Fr(1,2)}
# test: [oint T, oint X] is a total derivative for all X up to weight 8
ok = True
for w in range(2,9):
    for m in monomials(w,1):
        d = bracket_density(T0, {m:Fr1})
        if not is_total_derivative(d,1):
            ok=False; print("FAIL", m, d)
print("T commutes with all charges:", ok)
# also T(z)X(w) second-order pole should be weight*X: check coefficient n=-2
for w in range(2,7):
    for m in monomials(w,1):
        c2 = ope_coeff(T0,{m:Fr1},-2)
        assert c2 == {m:Fr(w)}, (m,c2)
print("weights ok")
# Table 1 dims (one boson, Z2-even sector), charge spins 2..10
print("dim Q_sigma sigma=2..10:", [class_dim(s+1,1,(0,)) for s in range(2,11)])
# :TT: at c=1
TT = nprod(T0,T0)
print(":TT: =", TT)
E = euler(TT,1)
print("Euler image of :TT:", E)
# compare with 1/4(J^4 - 2 J'^2)
cand = {mono(1,1,1,1):Fr(1,4), mono(2,2):Fr(-1,2)}
print("TT - cand total derivative?", is_total_derivative(padd(TT,cand,-1),1))
# ker ad_{I3 c=1} dims for sigma=2..10  (Table 1 KdV row)
I3 = {mono(1,1,1,1):Fr(1), mono(2,2):Fr(-2)}
print("KdV c=1 ker ad_I3 dims:", [kernel_of_ad(I3, s+1, 1, (0,))[0] for s in range(2,11)])
