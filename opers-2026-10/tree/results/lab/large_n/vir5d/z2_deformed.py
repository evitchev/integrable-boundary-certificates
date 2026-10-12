"""VIR5d Z2, operator side: first-order response of the WKB coefficient R_(2K) to a power deformation
    (1 - l0^2/T^2)^(1 + da e) (1 - l1^2/T^2)^(1 + db e)   of the generic-t operator (everything else fixed).
Usage: z2_deformed.py <k> <K> <da> <db> [stripped]   -> z2_def_k<..>_K<..>_<da>_<db>[_stripped].json
ODE side only; opens no data file.  Output: R0 (must vanish at a drop point) and R1 in (PX, pi)."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import mellin_dual as MD
SEAL = open('SEAL_VIR5d.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5d.md', 'rb').read()).hexdigest() == SEAL
k = sp.Rational(sys.argv[1]); K = int(sys.argv[2]); da = sp.Rational(sys.argv[3]); db = sp.Rational(sys.argv[4])
stripped = len(sys.argv) > 5 and sys.argv[5] == 'stripped'
n = k + 4
order = 2 * K
t0 = time.time()
h = MD.deformed_symbol(k, order // 2 + 1, da, db, stripped)
eng = MD.MellinWKB(n, 1 / k, h)
eng.solve(order)
Ri, ok = eng.R(order)
assert ok
PX, PI = sp.symbols('PX PI')
c2 = n**2 * (k + 1) / k
Pl = sp.Poly(Ri, MD.l0, MD.l1, MD.e)
parts = {0: {}, 1: {}}
for (a, b, ee), cf in Pl.terms():
    assert ee < 2 and a % 2 == 0 and b % 2 == 0
    parts[ee]['%d,%d' % (a, b)] = str(cf * c2 ** ((a + b) // 2))
out = {'k': str(k), 'K': K, 'da': str(da), 'db': str(db), 'stripped': stripped, 'R0_identically_zero': len(parts[0]) == 0,
       'R0': parts[0], 'R1': parts[1]}
tag = str(k).replace('/', '_')
fn = 'z2_def_k%s_K%d_%s_%s%s.json' % (tag, K, str(da).replace('-', 'm'), str(db).replace('-', 'm'), '_stripped' if stripped else '')
json.dump(out, open(fn, 'w'), indent=1, sort_keys=True)
print('k = %s K = %d deformation (%s, %s)%s: R0 identically zero: %s; R1 monomials: %d  (%.0fs)'
      % (k, K, da, db, ' stripped' if stripped else '', out['R0_identically_zero'], len(parts[1]), time.time() - t0))
