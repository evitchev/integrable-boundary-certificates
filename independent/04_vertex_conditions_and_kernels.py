from ope import *
from qsqrt2 import Q2, SQ2
import sympy as sp, pickle, time
Fr1=Fr(1)
al, kap, a = sp.symbols('alpha kappa a')
# (a) symbolic condition polynomial for I3 = J^4 + kappa J'^2
I3 = {mono(1,1,1,1):sp.Integer(1), mono(2,2):kap}
S = vertex_residue(I3, {0:al})
Us = monomials(2,1)
Ds = [D_alpha({u:sp.Integer(1)}, {0:al}) for u in Us]
idx = sorted(set(S)|set().union(*[set(d) for d in Ds]))
Mat = sp.Matrix([[d.get(k,0) for k in idx] for d in Ds] + [[S.get(k,0) for k in idx]]).T
print("weight-3 monomials:", idx, " matrix shape", Mat.shape)
det = sp.factor(Mat.det())
print("condition polynomial for J^4 + kappa J'^2:", det)
print("  at kappa=5/2:", sp.factor(det.subs(kap, sp.Rational(5,2))))
print("  at kappa=-2 :", sp.factor(det.subs(kap, -2)))
print("  at kappa=a^2-6+4/a^2:", sp.factor(det.subs(kap, a**2-6+4/a**2)))
# (b) Q(sqrt2) numerics: kernels of single momenta and joint kernel, weights 2..12
towerP6, towerI3p, towerI3k = pickle.load(open('towers.pkl','rb'))
one = Q2(1)
alphas = {'1/sqrt2': {0: SQ2*Fr(1,2)}, '2sqrt2': {0: SQ2*2}, 'sqrt2': {0: SQ2}, 'sqrt3ctrl': None}
def is_conserved(P, alpha, w):
    S = vertex_residue({k:one*v for k,v in P.items()}, alpha)
    Us = monomials(w-2,1) if w>=3 else [()]
    Ds = [D_alpha({u:one}, alpha) for u in Us]
    rows, idx = polys_to_matrix(Ds+[S])
    r1 = rank(rows[:-1], len(idx)); r2 = rank(rows, len(idx))
    return r1==r2
print("\nConservation table for tower members (ker ad_I5 reps) and KdV(c=1) I3:")
for w in sorted(towerP6):
    d, reps = towerP6[w]
    if d==0: continue
    P = reps[0]
    res = {name: is_conserved(P, alphas[name], w) for name in ['1/sqrt2','2sqrt2','sqrt2']}
    print(f"  spin {w-1}: {res}")
I3k = {mono(1,1,1,1):Fr1, mono(2,2):Fr(-2)}
print("  KdV c=1 I3:", {name: is_conserved(I3k, alphas[name], 4) for name in ['1/sqrt2','2sqrt2','sqrt2']})
# (c) genuine kernel dims of single momenta and joint kernel, weights 2..12; compare with tower
print("\nkernel dims (density weight w / charge spin w-1):")
t0=time.time()
for w in range(2,13):
    d1, dens1, mons = vertex_kernel(alphas['1/sqrt2'], w, 1, (0,), field_one=one)
    d2, dens2, _ = vertex_kernel(alphas['2sqrt2'], w, 1, (0,), field_one=one)
    d3, dens3, _ = vertex_kernel(alphas['sqrt2'], w, 1, (0,), field_one=one)
    # joint kernel of 1/sqrt2 and 2sqrt2
    if dens1 and dens2:
        inter = intersect_spaces([dens1,dens2], mons)
        dj = len(inter) - len(monomials(w-1,1,(0,)))
    else:
        dj = 0
    # is tower class inside kernel of 1/sqrt2 and inside 2sqrt2?
    dt, reps = towerP6.get(w,(0,[]))
    intower = None
    if dt:
        rows, _ = polys_to_matrix(dens1+[reps[0]], mons)
        in1 = rank(rows,len(mons)) == rank(rows[:-1],len(mons))
        rows, _ = polys_to_matrix(dens2+[reps[0]], mons)
        in2 = rank(rows,len(mons)) == rank(rows[:-1],len(mons))
        intower=(in1,in2)
    print(f"  w={w} spin={w-1}: ker(1/sqrt2)={d1}, ker(2sqrt2)={d2}, joint={dj}, ker(sqrt2)={d3}, tower dim={dt}, tower in ker(1/sqrt2),ker(2sqrt2)={intower}  [{time.time()-t0:.1f}s]")
