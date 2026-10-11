"""PROOF2 (1): factor Codex's Sol 3 laws L_m(a, K, t) after t -> t(nu, K) = (10 nu - 6K + 3)/(2 nu - 2K + 1) (nu = (2K-1)/n, n = 2(t-5)/(t-3)).
Print the irreducible factors of numerator and denominator, with their degrees in nu."""
import sys, json; sys.dont_write_bytecode = True
import sympy as sp
a, K, t, nu = sp.symbols('a K t nu')
tnu = (10 * nu - 6 * K + 3) / (2 * nu - 2 * K + 1)
for m in (1, 2, 3):
    d = json.load(open(f'codex_sol3_loss{m}.json'))['projection_over_c_a']
    loc = {'a': a, 'k': K, 't': t}
    L = sp.sympify(d['numerator'], locals=loc, rational=True) / sp.sympify(d['denominator'], locals=loc, rational=True)
    Ln = sp.factor(sp.together(L.subs(t, tnu)))
    num, den = sp.fraction(Ln)
    print(f'== loss {m}:')
    for lab, ex in (('NUM', num), ('DEN', den)):
        for f, e in sp.factor_list(ex)[1]:
            print(f'   {lab} factor (power {e}, deg_nu {sp.degree(f, nu)}): {f if len(str(f)) < 200 else str(f)[:200] + "..."}')
