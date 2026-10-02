"""Reconstruct fit equations EXACTLY from the seat's own code (run_second2.py / gate1.py functions imported, not reimplemented).
Usage: python3 gen_eqs.py MODE t [Mshift]
  MODE = h2b  : Sol 3, rotated dictionary, 4 spin-5 eqs + spin-7 X^3Y eq (as run_second2.test_fibre rotated=True)
  MODE = h2a  : Sol 3, standard dictionary, 4 spin-5 eqs (control)
  MODE = gate1: Sol 1, standard dictionary, fit on spin 3 monomials XY, X, Y, 1 (as gate1.run_case)
Writes eqs_<tag>.json, sing_<tag>.sing (saturation + minimal primes + lex GB per prime), msolve_<tag>.ms (fully expanded)."""
import sys, json; sys.dont_write_bytecode = True
from fractions import Fraction as F
import sympy as sp
import run_second2 as R
from run_second2 import X, Y, sX, sY, cX, cY, pi_rot, load_cert, wkb_poly, map_and_monic

import os
if os.environ.get('PLANT'):
    import plantcert; load_cert = plantcert.load_cert
mode = sys.argv[1]; tv = F(sys.argv[2]); shift = F(sys.argv[3]) if len(sys.argv) > 3 else F(0)
M = (tv + 3) / (tv - 1) + shift
tag = f"{mode}_t{str(tv).replace('/','o')}" + (f"_shift{str(shift).replace('/','o')}" if shift else "") + ("_plant" if os.environ.get('PLANT') else "")
unk = [sX, sY, cX, cY]
if mode == 'gate1':
    G = {}; exec(open('gate1.py').read().split('plant_ok =')[0], G)   # gate1's own cert() and wkb_mapped(), its run_case not executed
    C3 = G['cert'](2, tv); m3 = G['wkb_mapped'](M, 2)
    d3 = sp.Poly(sp.expand(m3 - C3), X, Y)
    rat = [sp.together(d3.coeff_monomial(m)) for m in (X*Y, X, Y, 1)]
else:
    rot = (mode == 'h2b')
    C5 = load_cert(3, 3, tv); e5, _ = wkb_poly(M, 3); m5 = map_and_monic(e5, 3, rot)
    d5 = sp.Poly(sp.expand(m5 - C5), X, Y)
    rat = [sp.together(d5.coeff_monomial(m)) for m in (X**2*Y, X**2, X*Y, X)]
    if rot:
        C7 = load_cert(3, 4, tv); e7, _ = wkb_poly(M, 4); m7 = map_and_monic(e7, 4, rot)
        d7 = sp.Poly(sp.expand(m7 - C7), X, Y)
        rat.append(sp.together(d7.coeff_monomial(X**3*Y))); unk.append(pi_rot)
eqs = [sp.expand(sp.numer(r)) for r in rat]; dens = [sp.denom(r) for r in rat]
print("tag", tag, "M", M, "unknowns", unk)
for i, (e, d) in enumerate(zip(eqs, dens)):
    print(f"eq[{i}] (total deg {sp.Poly(e, *unk).total_degree()}): {e}")
    print(f"  denom[{i}]: {sp.factor(d)}")
json.dump({'mode': mode, 't': str(tv), 'M': str(M), 'unknowns': [str(u) for u in unk], 'eqs': [str(e) for e in eqs], 'dens': [str(d) for d in dens]},
          open(f'eqs_{tag}.json', 'w'), indent=1)
vs = ','.join(str(u) for u in unk)
satf = [p for p, _ in sp.factor_list(sp.expand(sX * sp.prod(dens)))[1] if p.free_symbols]
hs = '*'.join('(' + str(p).replace('**', '^') + ')' for p in satf)
I = ',\n'.join(str(e).replace('**', '^') for e in eqs)
with open(f'sing_{tag}.sing', 'w') as f:
    f.write(f'ring r = 0, ({vs}), dp;\nideal I = {I};\nLIB "primdec.lib"; LIB "elim.lib";\n')
    f.write('ideal G = std(I); "RAW dim:"; dim(G);\n')
    f.write(f'poly h = {hs};\n')
    f.write('def SS = sat(I, h); "typeof sat:"; typeof(SS); ideal J = SS; ideal GJ = std(J); "SAT dim:"; dim(GJ); "SAT vdim (solutions with mult.):"; vdim(GJ);\n')
    f.write('ideal Jr = radical(J); "radical vdim:"; vdim(std(Jr));\n')
    f.write('list pd = minAssGTZ(J); "SAT minimal primes:"; size(pd);\n')
    f.write(f'ring s = 0, ({vs}), lp; option(redSB); list pd = imap(r, pd); int j;\n')
    f.write(f'for (j=1;j<=size(pd);j++) {{ ideal L = std(pd[j]); "--- SAT prime", j, " vdim", vdim(L); L; write(":w comp_{tag}_" + string(j) + ".txt", L); kill L; }}\n')
    f.write(f'write(":w ncomp_{tag}.txt", size(pd));\n')
    f.write('setring r; "RAW (unsaturated) decomposition, informational:"; list pr = minAssGTZ(I); int i; "RAW minimal primes:"; size(pr);\n')
    f.write('for (i=1;i<=size(pr);i++) { "  raw prime", i, "dim", dim(std(pr[i])), "inside h=0 (excluded):", reduce(h, std(pr[i]))==0, " sX=0 on it:", reduce(sX, std(pr[i]))==0; }\nquit;\n')
with open(f'msolve_{tag}.ms', 'w') as f:   # msolve's parser does not accept parenthesised products: everything expanded, no spaces
    zs = sp.Symbol('zs'); syms = unk + [zs]
    ren = lambda e: str(sp.expand(e)).replace('**', '^').replace(' ', '').replace('pi_rot', 'pr')
    f.write(','.join(str(u) for u in unk).replace('pi_rot', 'pr') + ',zs\n0\n' + ',\n'.join([ren(e) for e in eqs] + [ren(1 - zs * sp.prod(satf))]) + '\n')
