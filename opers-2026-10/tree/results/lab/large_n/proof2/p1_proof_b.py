"""PROOF2 (1), supplement: the infinite coincidence families nu = j0 (j0 < m) AND a = K - j0 - j (j = 1..m, j0 + j >= m), K symbolic, where the
binomial zero order is z = 2 (for K >= 2 j0 + j + 1) against two simple poles.  Compute ord_eps(N) - ord_eps(D) + z symbolically in K."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
from p1_proof import law, eps_order, a, K, nu, eps
ok = True
for m in (2, 3):
    N, D = law(m)
    for j0 in range(1, m):
        for j in range(1, m + 1):
            if j0 + j < m: continue
            sub = {a: K - j0 - j, nu: j0 + eps}
            on, cn = eps_order(N.subs(sub)); od, cd = eps_order(D.subs(sub))
            z = 2; tot = (on if on is not None else 99) - od + z
            roots = sp.solve(cd, K) if cd.free_symbols else []
            print(f'loss {m}: nu = {j0}, a = K - {j0 + j} (K >= {2 * j0 + j + 1}): ord N = {on}, ord D = {od}, z = 2 -> total {tot}: {"VANISHES" if tot >= 1 else "NOT VANISHING"}; extra-coincidence K (in PART B range if small): {roots}', flush=True)
            ok &= tot >= 1
print('SUPPLEMENT', 'PASS' if ok else 'FAIL'); sys.exit(0 if ok else 1)
