import os
REPO=os.environ.get('REPO_RESULTS', os.path.expanduser('~/integrable-boundary-certificates/results'))+'/'
from ope import *
import json, sympy as sp, time, sys
Fr1=Fr(1)
d=json.load(open(REPO+'paperclip_P8.json'))
nsym=sp.Symbol('n')
def parse_mono(s):
    # "(33|11)" -> X orders (3,3), Y orders (1,1)
    s=s.strip('()'); x,y=s.split('|')
    return M(tuple(int(c) for c in x), tuple(int(c) for c in y))
coeff_exprs=[sp.sympify(c) for c in d['coefficients']]
basis=[parse_mono(b) for b in d['basis']]
def P8_at(nval):
    P={}
    for m,c in zip(basis,coeff_exprs):
        v=sp.Rational(c.subs(nsym,nval))
        P[m]=Fr(int(v.p),int(v.q))
    return P
def screenings(nval):
    nval=Fr(nval)
    a1=-nval/2; a2=(nval+2)/2
    # need rational square roots
    def rsqrt(x):
        from math import isqrt
        p,q=x.numerator,x.denominator
        sp_,sq=isqrt(p),isqrt(q)
        assert sp_*sp_==p and sq*sq==q, x
        return Fr(sp_,sq)
    r1,r2=rsqrt(a1),rsqrt(a2)
    return [{0:s1*r1,1:s2*r2} for s1 in (1,-1) for s2 in (1,-1)]
def joint_vertex_kernel(alphas, w):
    dens_lists=[]; mons=None
    for al in alphas:
        dd, dens, mons = vertex_kernel(al, w, 2, (0,0))
        dens_lists.append(dens)
        if not dens: return 0, [], mons
    inter = intersect_spaces(dens_lists, mons)
    dim = len(inter) - len(monomials(w-1,2,(0,0)))
    # representatives: pick vectors whose E-images are independent
    reps=[]; cur=[]
    for v in inter:
        dd={m:c for c,m in zip(v,mons) if c!=0}
        ev=euler_vec(dd,2)
        trial=cur+[ev]
        rows,idx=polys_to_matrix(trial)
        if rank(rows,len(idx))>len(cur):
            cur=trial; reps.append(dd)
    return dim, reps, mons
for nval in [Fr(-18,25), Fr(-50,169)]:
    print(f"\n===== n = {nval}  screening vectors: {screenings(nval)}")
    t0=time.time()
    P8=P8_at(nval)
    als=screenings(nval)
    # (a) joint kernel dims at weights 2..8 and P8 membership at weight 8
    for w in range(2,9):
        dim, reps, mons = joint_vertex_kernel(als, w)
        extra=""
        if w==4: I3=reps[0]
        if w==6: I5=reps[0]
        if w==8:
            # is P8 in the span of the joint kernel (density level)? test rank
            rows,_=polys_to_matrix(reps+[P8], mons)
            # need full kernel span incl. total derivatives: use inter vectors instead
            # simpler: check P8 conserved under each screening directly
            cons=[]
            for al in als:
                S=vertex_residue(P8, al)
                Us=monomials(6,2,None)
                Ds=[D_alpha({u:Fr1},al) for u in Us]
                rws,idx=polys_to_matrix(Ds+[S])
                cons.append(rank(rws[:-1],len(idx))==rank(rws,len(idx)))
            extra=f"  P8 conserved under the four screenings: {cons}"
        print(f"  spin {w-1}: joint screening kernel dim = {dim}{extra}   [{time.time()-t0:.1f}s]")
    # (b) commutators
    print("  [I3pc, I5pc] total derivative:", is_total_derivative(bracket_density(I3,I5),2))
    print("  [I3pc, P8 ] total derivative:", is_total_derivative(bracket_density(I3,P8),2), f"[{time.time()-t0:.1f}s]")
    print("  [I5pc, P8 ] total derivative:", is_total_derivative(bracket_density(I5,P8),2), f"[{time.time()-t0:.1f}s]")
    # is P8 proportional to the joint-kernel rep at spin 7 (mod d)? compare Euler images
    dim8, reps8, mons8 = joint_vertex_kernel(als, 8)
    ev=[euler_vec(reps8[0],2), euler_vec(P8,2)]
    rows,idx=polys_to_matrix(ev)
    print("  P8 class == joint-kernel class at spin 7:", rank(rows,len(idx))==1, f"[{time.time()-t0:.1f}s]")
    # (c) ker ad_{I3pc} dims spins 2..8
    dims=[]
    for w in range(3,10):
        dd,_=kernel_of_ad(I3, w, 2, (0,0))
        dims.append(dd)
    print("  ker ad_I3pc dims spins 2..8:", dims, f"[{time.time()-t0:.1f}s]")
    print("  dim Q_sigma paperclip spins 2..8:", [class_dim(w,2,(0,0)) for w in range(3,10)])
