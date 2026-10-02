"""PROOF1 (3) negative controls: the same symbolic pipeline with (i) the momentum scale CY x 11/10, (ii) the string length a + 1/10 (Sol 1) /
k + 1/10 in the strings only (Sols 2, 3).  Each must make the identity FAIL."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
import p3_g1_ward as P
fails = 0
orig = P.solution_data
for name, mod in (('CY x 11/10', lambda D: {**D, 'CY': D['CY'] * sp.Rational(11, 10)}),
                  ('string length + 1/10', lambda D: {**D, 'strings': [(f, s, a + sp.Rational(1, 10)) for f, s, a in D['strings']]})):
    P.solution_data = (lambda m: (lambda sol: m(orig(sol))))(mod)
    for sol in (1, 2, 3):
        d = sp.cancel(sp.together(P.g1_ratio(sol) - P.codex_ratio(sol)))
        print(f'TAMPER {name}: Sol {sol}: identity holds = {d == 0}  (must be False)', flush=True); fails += d == 0
P.solution_data = orig
sys.exit(0 if fails == 0 else 1)
