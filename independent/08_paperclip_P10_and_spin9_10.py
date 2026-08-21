import os
REPO=os.environ.get('REPO_RESULTS', os.path.expanduser('~/integrable-boundary-certificates/results'))+'/'
from ope import *
import json, sympy as sp, time
Fr1=Fr(1)
nsym=sp.Symbol('n')
def load(path):
    d=json.load(open(path))
    def parse_mono(s):
        s=s.strip('()'); x,y=s.split('|')
        return M(tuple(int(c) for c in x), tuple(int(c) for c in y))
    key='basis' if 'basis' in d else 'basis_pretty'
    return [parse_mono(b) for b in d[key]], [sp.sympify(c) for c in d['coefficients']]
b8,c8=load(REPO+'paperclip_P8.json')
b10,c10=load(REPO+'paperclip_P10.json')
def at(basis,coeffs,nval):
    P={}
    for m,c in zip(basis,coeffs):
        v=sp.Rational(c.subs(nsym,nval)); P[m]=Fr(int(v.p),int(v.q))
    return P
nval=Fr(-18,25)
als=[{0:s1*Fr(3,5),1:s2*Fr(4,5)} for s1 in (1,-1) for s2 in (1,-1)]
P8=at(b8,c8,nval); P10=at(b10,c10,nval)
t0=time.time()
# I3pc from joint kernel at weight 4
def joint_reps(als,w):
    dens_lists=[]; mons=None
    for al in als:
        dd,dens,mons=vertex_kernel(al,w,2,(0,0)); dens_lists.append(dens)
    inter=intersect_spaces(dens_lists,mons)
    dim=len(inter)-len(monomials(w-1,2,(0,0)))
    reps=[];cur=[]
    for v in inter:
        dd={m:c for c,m in zip(v,mons) if c!=0}
        ev=euler_vec(dd,2); trial=cur+[ev]; rows,idx=polys_to_matrix(trial)
        if rank(rows,len(idx))>len(cur): cur=trial; reps.append(dd)
    return dim,reps
_,r4=joint_reps(als,4); I3=r4[0]
_,r6=joint_reps(als,6); I5=r6[0]
cons=[]
for al in als:
    S=vertex_residue(P10,al); Us=monomials(8,2,None); Ds=[D_alpha({u:Fr1},al) for u in Us]
    rws,idx=polys_to_matrix(Ds+[S]); cons.append(rank(rws[:-1],len(idx))==rank(rws,len(idx)))
print("P10 conserved under the four screenings:", cons, f"[{time.time()-t0:.1f}s]")
print("[I3pc,P10] total derivative:", is_total_derivative(bracket_density(I3,P10),2), f"[{time.time()-t0:.1f}s]")
print("[I5pc,P10] total derivative:", is_total_derivative(bracket_density(I5,P10),2), f"[{time.time()-t0:.1f}s]")
print("[P8,P10]  total derivative:", is_total_derivative(bracket_density(P8,P10),2), f"[{time.time()-t0:.1f}s]")
print("[P8,P8-check: [I5,P8]] :", is_total_derivative(bracket_density(I5,P8),2))
for w in (10,11):
    dd,_=kernel_of_ad(I3,w,2,(0,0))
    print(f"ker ad_I3pc at spin {w-1}: {dd}  (dim Q = {class_dim(w,2,(0,0))})  [{time.time()-t0:.1f}s]")
