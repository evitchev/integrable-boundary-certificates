"""VIR2-R control (e): the formal-N pair engine (onwick) agrees EXACTLY with the explicit-component engine (wick) at N = 3, 4 on :T T:, :T :T T::,
:dT T:, and on a mixed X-Y product :(dX dX T_Y) T_Y:.  Exit 0 iff all agree."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import onwick as O, wick as Wk
ok = True
def canon(poly): return {tuple(sorted(m, key=str)): v for m, v in poly.items() if v}
for Nc in (3, 4):
    Tp = {((), ((1, 1),)): F(1, 2)}
    Te = {}
    for i in range(Nc): Te = Wk.fadd(Te, {((i, 0), (i, 0)): F(1, 2)})
    tests = []
    TTp = O.ope_pairs(Tp, Tp, Nc, (0,))[0]; TTe = Wk.ope(Te, Te, poles=(0,))[0]; tests.append((':TT:', TTp, TTe))
    TTTp = O.ope_pairs(Tp, TTp, Nc, (0,))[0]; TTTe = Wk.ope(Te, TTe, poles=(0,))[0]; tests.append((':T:TT::', TTTp, TTTe))
    dTp = O.dfield_pairs(Tp); dTe = Wk.dfield(Te)
    tests.append((':dT T:', O.ope_pairs(dTp, Tp, Nc, (0,))[0], Wk.ope(dTe, Te, poles=(0,))[0]))
    tests.append(('T T pole 4', O.ope_pairs(Tp, Tp, Nc, (4,)).get(4, {}), Wk.ope(Te, Te, poles=(4,)).get(4, {})))
    Mp = O.mulx((1, 1), Tp); Me = {}
    for m, v in Te.items(): Me[tuple(sorted(m + ((-1, 0), (-1, 0))))] = v
    tests.append((':(dXdX T) T: pole 2', O.ope_pairs(Mp, Tp, Nc, (2,)).get(2, {}), Wk.ope(Me, Te, poles=(2,)).get(2, {})))
    for name, fp, fe in tests:
        a = canon(O.to_explicit(fp, Nc)); b = canon({tuple(sorted(m, key=str)): v for m, v in fe.items()})
        good = a == b; ok &= good
        print(f'N = {Nc}: {name}: formal-N engine == explicit engine: {good}  ({len(a)} monomials)')
print('control (e):', ok); sys.exit(0 if ok else 1)
