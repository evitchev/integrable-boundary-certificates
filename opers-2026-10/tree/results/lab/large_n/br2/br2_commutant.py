"""BR2 (2)/(3): independent replication (own engine eng.py, SUPER2) of 'the screening triple determines the family':
commutant MOD TOTAL DERIVATIVES of {V_+, V_-, e^(p phi)} at weights 4, 6 and 8 (8 = hold-out), and of the pair alone at weight 4; tamper a -> a + 1/7.
Densities P (two-boson, weight w) conserved by screening s iff Res P(z) e^{s Phi}(w) = (d + s.dPhi) Y for some weight-(w-2) Y (eng.residue)."""
import sys; sys.dont_write_bytecode = True
from fractions import Fraction as F
from eng import basis, residue, nullspace, deriv, pmul, padd
def tw(Y, s):
    out = deriv(Y)
    for aa in (0, 1):
        if s[aa]: padd(out, pmul({((aa, 1),): s[aa]}, Y))
    return out
def rank(vs):
    if not vs: return 0
    return len(vs) - len(nullspace([{i: x for i, x in enumerate(col) if x} for col in zip(*vs)], len(vs)))
def dim_mod_d(w, screens):
    Bw, Bl = basis(w), basis(w - 2); n = len(Bw); nY = len(Bl); rows = {}
    for si, s in enumerate(screens):
        for j, m in enumerate(Bw):
            for om, c in residue(m, s).items(): rows.setdefault((si, om), {})[j] = c
        for j, m in enumerate(Bl):
            for om, c in tw({m: F(1)}, s).items(): rows.setdefault((si, om), {})[n + si * nY + j] = -c
    ns = nullspace(list(rows.values()), n + len(screens) * nY)
    P = [v[:n] for v in ns]
    Dv = [[deriv({m: F(1)}).get(b, F(0)) for b in Bw] for m in basis(w - 1)]
    return rank(P + Dv) - rank(Dv)
if __name__ == '__main__':
    FIB = {1: [(F(6, 5), F(4, 5)), (F(8, 5), F(3, 5))], 2: [(F(3, 5), F(4, 5)), (F(5, 13), F(12, 13))], 3: [(F(3, 2), None), (F(2, 3), None)]}
    bad = 0
    for sol, fibs in FIB.items():
        for p, a in fibs:
            beta = {1: -p / 2, 2: -p, 3: 1 / p}[sol]; a = a if a is not None else 1 / p
            if sol in (1, 2): assert a * a == {1: 1 - p * p / 4, 2: 1 - p * p}[sol]
            Vp, Vm, Vv = (a, beta), (-a, beta), (F(0), p)
            tri = [dim_mod_d(w, [Vp, Vm, Vv]) for w in (4, 6, 8)]
            pair4 = dim_mod_d(4, [Vp, Vm])
            tam = dim_mod_d(4, [(a + F(1, 7), beta), (-a - F(1, 7), beta), Vv])
            ok = tri == [1, 1, 1] and tam == 0
            bad += not ok
            print(f'Sol {sol}, p = {p}, a = {a}, beta = {beta}: triple mod d at w = 4, 6, 8: {tri} (expect [1,1,1]); pair alone w = 4: {pair4}; TAMPER a + 1/7 (w = 4): {tam} (expect 0)  {"OK" if ok else "FAIL"}', flush=True)
    print('REPLICATION', 'PASS' if not bad else f'FAIL ({bad})'); sys.exit(1 if bad else 0)
