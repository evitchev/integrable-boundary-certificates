"""Direct mutual-commutativity checks for the cylindrical towers
(audit item: the bootstrap lemma covers only [I_3, .] within the
verified vacancy range; [I_5,I_7] and [I_5,I_9] land at spins 12, 14).

For each Solution at curve points t = 2 and t = 7: solve the unique
I_5, I_7, I_9 kernels and verify [I_5, I_7] = [I_5, I_9] = 0 directly.
"""
from fractions import Fraction as F
from mixed_engine import (gen_mixed_basis, genuine_kernel_mixed,
                          MixedSector, residue_mixed_cur, cyl_P4,
                          spin_mixed)
from invariant_engine import curve_point

def is_td(cur, s):
    sec = MixedSector(s)
    return not any(sec.reduce_mod_d(sec.vec(cur)))

for t in (F(2), F(7)):
    N, s = curve_point(t)
    for sol in (1, 2, 3):
        P4 = cyl_P4(sol, N, s)
        ch = {}
        for sig in (5, 7, 9):
            g = genuine_kernel_mixed(P4, gen_mixed_basis(sig + 1), N)
            assert len(g) == 1, (sol, sig, len(g))
            ch[sig] = g[0]
        r57 = residue_mixed_cur(ch[5], ch[7], N)
        r59 = residue_mixed_cur(ch[5], ch[9], N)
        ok57 = is_td(r57, 13)
        ok59 = is_td(r59, 15)
        print(f"t={t} Sol{sol}: [I5,I7]=0: {ok57}  [I5,I9]=0: {ok59}")
        assert ok57 and ok59
print("All direct cylindrical mutual-commutativity checks passed.")
