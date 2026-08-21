from ope import *
from qsqrt2 import Q2, SQ2
import time, pickle
Fr1=Fr(1); one=Q2(1)
P6 = {mono(3,3):Fr(31,24), mono(1,1,2,2):Fr(5,2), mono(1,1,1,1,1,1):Fr(1)}
I3p = {mono(1,1,1,1):Fr1, mono(2,2):Fr(5,2)}
t0=time.time()
def span_eq(A,B):
    if len(A)!=len(B): return False
    if not A: return True
    ev=[euler_vec(r,1) for r in A+B]; rows,idx=polys_to_matrix(ev)
    return rank(rows,len(idx))==len(A)
def is_conserved(P, alpha, w):
    S = vertex_residue({k:one*v for k,v in P.items()}, alpha)
    Us = monomials(w-2,1)
    Ds = [D_alpha({u:one}, alpha) for u in Us]
    rows, idx = polys_to_matrix(Ds+[S])
    return rank(rows[:-1], len(idx)) == rank(rows, len(idx))
for w in (15,16):
    d1,r1=kernel_of_ad(P6,w,1,(0,))
    print(f"spin {w-1}: dim ker ad_I5 = {d1}  [{time.time()-t0:.1f}s]", flush=True)
    if d1:
        d2,r2=kernel_of_ad(I3p,w,1,(0,))
        print(f"   dim ker ad_I'3 (KdV c=-25/2) = {d2}, classes equal: {span_eq(r1,r2)}  [{time.time()-t0:.1f}s]", flush=True)
        for name,al in [("1/sqrt2",{0:SQ2*Fr(1,2)}),("sqrt2",{0:SQ2})]:
            print(f"   I'_{w-1} conserved under e^{{{name} phi}}: {is_conserved(r1[0],al,w)}  [{time.time()-t0:.1f}s]", flush=True)
