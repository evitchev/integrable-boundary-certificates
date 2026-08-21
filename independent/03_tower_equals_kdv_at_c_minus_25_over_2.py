from ope import *
import time
Fr1=Fr(1)
P6 = {mono(3,3):Fr(31,24), mono(1,1,2,2):Fr(5,2), mono(1,1,1,1,1,1):Fr(1)}
I3p = {mono(1,1,1,1):Fr1, mono(2,2):Fr(5,2)}   # enhanced I'_3 = KdV I3 at Q^2=9/8
I3k = {mono(1,1,1,1):Fr1, mono(2,2):Fr(-2)}   # KdV c=1
t0=time.time()
towerP6 = {}
towerI3p = {}
towerI3k = {}
for w in range(3,15):
    d1, reps1 = kernel_of_ad(P6, w, 1, (0,))
    d2, reps2 = kernel_of_ad(I3p, w, 1, (0,))
    d3, reps3 = kernel_of_ad(I3k, w, 1, (0,))
    towerP6[w]=(d1,reps1); towerI3p[w]=(d2,reps2); towerI3k[w]=(d3,reps3)
    # compare classes: are the E-image spans equal?
    def span_eq(repsA, repsB):
        if len(repsA)!=len(repsB): return False
        if not repsA: return True
        ev = [euler_vec(r,1) for r in repsA+repsB]
        rows, idx = polys_to_matrix(ev)
        return rank(rows,len(idx)) == len(repsA)
    print(f"spin {w-1}: dim ker ad_I5={d1}, dim ker ad_I'3(KdV c=-25/2)={d2}, dim ker ad_I3(KdV c=1)={d3}, "
          f"ker ad_I5 == ker ad_I'3: {span_eq(reps1,reps2)}, ker ad_I5 == KdV(c=1): {span_eq(reps1,reps3)}   [{time.time()-t0:.1f}s]")
import pickle
pickle.dump((towerP6,towerI3p,towerI3k), open('towers.pkl','wb'))
