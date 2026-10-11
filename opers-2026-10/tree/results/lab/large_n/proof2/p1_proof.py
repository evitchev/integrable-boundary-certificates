"""PROOF2 (1): rank-drop vanishing of the RECORD's Sol 3 losses m = 1 (control), 2, 3 on nu = i, 2i <= K - 1, for 0 <= a <= K - m.
coef_m(a) = L_m(a, K, nu) * binom(nu, a) binom(nu, K - a) (up to a-independent nonzero factors).
PART A (symbolic proof): poles of L_m on the locus only come from (nu - j0), j0 < m, and (-K + a + nu + j), j = 1..m (p1_factor.log; (2K-4nu-1), (2K-3nu-1)
never vanish on the locus).  Outside these families L_m is finite and the binomial product has a zero -> coef = 0.  Inside, enumerate the sub-families with
binomial zero order z = 1 (all others have z = 2 against a single pole), and compute the epsilon-order of L_m at nu = i + epsilon SYMBOLICALLY in the free
parameter (K or i).  Coincident-pole parameter values are listed and checked exactly in PART B.
PART B (exact check): every (i, K, a) with i = 1..6, K = 2i+1 .. 2i+16, 0 <= a <= K - m: lim_(eps -> 0) of the exact rational function in eps (binomials are
polynomials in eps at integer a, K)."""
import sys, json; sys.dont_write_bytecode = True
import sympy as sp
a, K, t, nu, eps, i = sp.symbols('a K t nu epsilon i')
tnu = (10 * nu - 6 * K + 3) / (2 * nu - 2 * K + 1)
def law(m):
    d = json.load(open(f'codex_sol3_loss{m}.json'))['projection_over_c_a']; loc = {'a': a, 'k': K, 't': t}
    L = sp.sympify(d['numerator'], locals=loc, rational=True) / sp.sympify(d['denominator'], locals=loc, rational=True)
    return sp.fraction(sp.factor(sp.together(L.subs(t, tnu))))
def eps_order(poly):
    P = sp.Poly(sp.expand(poly), eps)
    for k in range(P.degree() + 1):
        c = sp.factor(P.coeff_monomial(eps**k))
        if c != 0: return k, c
    return None, 0
def binom_eps(av, Kv, iv):   # exact polynomial in eps for integer a, K, i
    f = lambda x, n: sp.Mul(*[(x - j) for j in range(n)]) / sp.factorial(n)
    return sp.expand(f(iv + eps, av) * f(iv + eps, Kv - av))
def main():
  ok = True
  for m in (1, 2, 3):
      N, D = law(m)
      print(f'== loss {m}', flush=True)
      # PART A: z = 1 sub-families
      fams = []
      for j0 in range(1, m):
          # i = j0; a in {0..j0} (fixed integers, K symbolic, K >= 2 j0 + 1 and a <= K - m) and a = K - c, c in [m, j0] (fixed offset)
          for av in range(0, j0 + 1): fams.append((f'(nu - {j0}) pole, i = {j0}, a = {av}, K symbolic', {a: av, nu: j0 + eps}, 'K'))
          for c in range(m, j0 + 1): fams.append((f'(nu - {j0}) pole, i = {j0}, a = K - {c}', {a: K - c, nu: j0 + eps}, 'K'))
      for j in range(1, m + 1):
          for r in range(1, j + 1):      # K = 2i + r, a = K - i - j = i + r - j
              fams.append((f'(-K+a+nu+{j}) pole, K = 2i + {r}, a = i + {r} - {j}, i symbolic', {K: 2 * i + r, a: i + r - j, nu: i + eps}, 'i'))
      for lab, sub, free in fams:
          Ns, Ds = N.subs(sub), D.subs(sub)
          if free == 'K': Ns, Ds = Ns.subs(K, K), Ds
          on, cn = eps_order(Ns); od, cd = eps_order(Ds)
          z = 1
          tot = (on if on is not None else 99) - od + z
          status = 'VANISHES' if tot >= 1 else 'NOT VANISHING'
          roots = sp.solve(cd, sp.Symbol(free)) if cd.free_symbols else []
          print(f'  {lab}: ord_eps(N) = {on}, ord_eps(D) = {od}, z = {z} -> total {tot}: {status}; coincident-pole values of {free} (checked in PART B): {roots}', flush=True)
          ok &= tot >= 1
      # PART B: exact grid
      bad = []
      cnt = 0
      for iv in range(1, 7):
          for Kv in range(2 * iv + 1, 2 * iv + 17):
              for av in range(0, Kv - m + 1):
                  val = sp.cancel(N.subs({a: av, K: Kv, nu: iv + eps}) / D.subs({a: av, K: Kv, nu: iv + eps}) * binom_eps(av, Kv, iv))
                  lim = sp.limit(val, eps, 0) if sp.denom(val).subs(eps, 0) == 0 else val.subs(eps, 0)
                  cnt += 1
                  if lim != 0: bad.append((iv, Kv, av, lim))
      print(f'  PART B: {cnt} exact cells (i = 1..6, K = 2i+1..2i+16): nonzero limits {len(bad)}: {bad[:5]}', flush=True)
      ok &= not bad
      # NEGATIVE CONTROL: off the locus (2i = K, i.e. K even, nu = K/2) the layer must NOT vanish identically
      nz = 0
      for Kv in (4, 6, 8):
          iv = Kv // 2
          for av in range(0, Kv - m + 1):
              val = sp.cancel(N.subs({a: av, K: Kv, nu: iv + eps}) / D.subs({a: av, K: Kv, nu: iv + eps}) * binom_eps(av, Kv, iv))
              lim = sp.limit(val, eps, 0)
              nz += lim != 0
      print(f'  NEGATIVE CONTROL (2i = K, outside the locus): nonzero coefficients found: {nz} (must be > 0)', flush=True)
      ok &= nz > 0
  print('PROOF2 (1)', 'PASS' if ok else 'FAIL'); sys.exit(0 if ok else 1)
if __name__ == '__main__': main()
