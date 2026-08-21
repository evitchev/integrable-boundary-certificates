from ope import *
from qsqrt2 import Q2, SQ2
import pickle, time
one=Q2(1)
towerP6, towerI3p, towerI3k = pickle.load(open('towers.pkl','rb'))
aA = {0: SQ2}            # sqrt2  (c=1 Virasoro screening = SU(2) current)
aB = {0: SQ2*Fr(-1,2)}   # -1/sqrt2 (A2(2) twisted screening at c=1)
t0=time.time()
print("joint kernel of the c=1 A2(2) screening pair {sqrt2, -1/sqrt2} in the Z2-even sector, vs L(1,0) cap KdV(1/sqrt2):")
for w in range(2,15):
    d1, dens1, mons = vertex_kernel(aA, w, 1, (0,), field_one=one)
    d2, dens2, _ = vertex_kernel(aB, w, 1, (0,), field_one=one)
    if dens1 and dens2:
        inter = intersect_spaces([dens1,dens2], mons)
        dj = len(inter) - len(monomials(w-1,1,(0,)))
    else:
        dj = 0
    dt = towerP6.get(w,(0,[]))[0]
    print(f"  spin {w-1}: dim joint = {dj}   (ker sqrt2 = {d1}, ker -1/sqrt2 = {d2}, tower = {dt})  [{time.time()-t0:.1f}s]")
