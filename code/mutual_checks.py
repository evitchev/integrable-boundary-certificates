"""Mutual-commutativity checks among promoted Q(t) charges.

Pairs: cylindrical [I5,I11],[I7,I9],[I7,I11],[I9,I11] x Sol 1,2,3; A2 [I7,I11].
Charges evaluated from STORED Q(t)-symbolic densities (fixed rational
sections; qt_promotion.json), plus the certified A2 I7 closed form and
per-point normalized kernel solves for cyl I5.  Each pair checked at up
to ~70 rational t; entries of the reduced commutator are values of
fixed rational sections of t, so the many-point zeros are strong exact
evidence (NOT a completed interpolation certificate: no
post-substitution degree bound is computed here; the vacancy route in
vacancy_high.py supersedes these checks where completed).
Deferred as engine-infeasible (needs symmetry-collapsed Wick counting):
A2 [I7,I13] and [I11,I13].
"""
import json, sys
from fractions import Fraction as F
from multiprocessing import Pool
import sympy as sp
from invariant_engine import (InvSector, gen_inv_basis, residue_inv_cur,
                              a2_P6, curve_point, canon, canon_pair,
                              genuine_kernel_inv)
from mixed_engine import (MixedSector, gen_mixed_basis, residue_mixed_cur,
                          cyl_P4, genuine_kernel_mixed)

D = json.load(open("../results/qt_promotion.json"))
Nsym = sp.Symbol('N')
CACHE = {}
def stored(key, sigma, mixed):
    if key not in CACHE:
        basis = gen_mixed_basis(sigma+1) if mixed else gen_inv_basis(sigma+1)
        CACHE[key] = [(basis[int(i)], sp.sympify(a, locals={'N': Nsym}),
                       sp.sympify(b, locals={'N': Nsym}))
                      for i,(a,b) in D[key]["coeffs_A_B_of_N"].items()]
    out = CACHE[key]
    def ev(t):
        Nv, sv = curve_point(t)
        NN = sp.Rational(Nv)
        cur = {}
        for m, A, B in out:
            v = sp.Rational(A.subs(Nsym, NN)) + sp.Rational(sv)*sp.Rational(B.subs(Nsym, NN))
            if v != 0:
                cur[m] = F(int(sp.numer(v)), int(sp.denom(v)))
        return cur
    return ev

def a2_I7(t):
    Nv, sv = curve_point(t)
    Bmon = [canon((canon_pair(1,2),canon_pair(2,3))), canon(((1,3),(1,3))),
            canon(((1,1),(3,3))), canon(((2,2),(2,2))), canon(((4,4),)),
            canon(((1,1),(1,1),(2,2))), canon((canon_pair(1,2),canon_pair(1,2),(1,1))),
            canon(((1,1),)*4)]
    A=[F(1207,40)+F(104,3)*Nv-F(65,24)*Nv**2, -F(1057,160)-F(281,48)*Nv+F(65,96)*Nv**2,
       F(131,16)+F(9,16)*Nv, F(2977,160)+F(281,48)*Nv-F(65,96)*Nv**2,
       -F(1881,1600)+F(11,1440)*Nv-F(13,576)*Nv**2, F(-12), -F(9,4)-F(27,4)*Nv, F(1)]
    B=[F(203,40)-F(21,8)*Nv, F(147,160)+F(21,32)*Nv, F(7,16), -F(147,160)-F(21,32)*Nv,
       -F(847,4800)-F(7,320)*Nv, F(0), -F(21,4), F(0)]
    return {m: a+sv*b for m,a,b in zip(Bmon,A,B) if a+sv*b}

def cyl_I5(sol):
    basis = gen_mixed_basis(6); sec = MixedSector(6)
    def ev(t):
        g = genuine_kernel_mixed(cyl_P4(sol,*curve_point(t)), basis, curve_point(t)[0])
        if len(g)!=1: return None
        v = sec.reduce_mod_d(sec.vec(g[0]))
        nz = next(i for i,x in enumerate(v) if x)
        return {basis[i]: x/v[nz] for i,x in enumerate(v) if x}
    return ev

TS = [F(k) for k in range(2,50) if k!=5] + [F(2*k+1,2) for k in range(1,25)]

def task(args):
    fam, sol, s1, s2 = args
    if fam=="a2":
        A = a2_P6 and (a2_I7 if s1==7 else stored(f"a2_sigma{s1}",s1,False))
        Bv = stored(f"a2_sigma{s2}",s2,False)
        res, Sec = residue_inv_cur, InvSector
        getA = a2_I7 if s1==7 else (lambda t: stored(f"a2_sigma{s1}",s1,False)(t))
    else:
        res, Sec = residue_mixed_cur, MixedSector
        getA = cyl_I5(sol) if s1==5 else stored(f"cyl{sol}_sigma{s1}",s1,True)
        Bv = stored(f"cyl{sol}_sigma{s2}",s2,True)
    tgt = Sec(s1+s2+1)
    npass = nskip = 0
    for t in TS:
        a = getA(t)
        if a is None: nskip+=1; continue
        b = Bv(t)
        Nv = curve_point(t)[0]
        r = res(a,b,Nv)
        if not any(tgt.reduce_mod_d(tgt.vec(r))):
            npass+=1
        else:
            print(f"FAIL {fam}{sol or ''} [{s1},{s2}] at t={t}", flush=True)
            return (f"{fam}{sol or ''}[{s1},{s2}]", False, npass, nskip)
    print(f"done {fam}{sol or ''} [I{s1},I{s2}]=0 at {npass} pts ({nskip} skipped)", flush=True)
    return (f"{fam}{sol or ''}[{s1},{s2}]", True, npass, nskip)

if __name__=="__main__":
    tasks=[("cyl",s,a,b) for s in (1,2,3) for a,b in ((5,11),(7,9),(7,11),(9,11))]
    tasks.append(("a2",None,7,11))
    with Pool(len(tasks)) as pool:
        results = pool.map(task, tasks)
    ok = all(r[1] for r in results)
    json.dump({k:{"pass":p,"points":n,"skipped":sk} for k,p,n,sk in results},
              open("../results/mutual_checks.json" if ok else "../results/mutual_checks.partial.json","w"), indent=1)
    print("ALL PASSED" if ok else "FAILURES", flush=True)
    sys.exit(0 if ok else 1)
