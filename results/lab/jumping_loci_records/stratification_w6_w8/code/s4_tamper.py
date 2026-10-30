"""stratification_w6_w8 S4(b) tamper: sheet 3, W = 8, the I_3 coefficient of X(1,1,1,1) multiplied by 11/10 (1 -> 11/10).
Builds the tampered M(t), reports exactly what changed (I_3 coefficients, M entries), saves data/M_sheet3_W8_tamper.json.
S2-S3 are then run on it by s2s3.py (same code as the untampered runs)."""
import os
import sys

from sympy import QQ

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family as fam  # noqa: E402

me = fam.me
SHEET, W = 3, 8
key = me._M((1, 1, 1, 1))
P4 = fam.p4_family(SHEET)
P4t = dict(P4)
P4t[key] = P4[key] * fam.K(QQ(11, 10))
changed = [m for m in P4 if P4[m] != P4t[m]]
print(f"sheet {SHEET}: I_3 coefficients changed: {len(changed)}")
for m in changed:
    print(f"    {fam.mono_str(m)}: {P4[m]} -> {P4t[m]}")
assert len(changed) == 1 and set(P4) == set(P4t)
M0, rows, cols = fam.build_M(P4, W)
M1, rows1, cols1 = fam.build_M(P4t, W)
assert rows == rows1 and cols == cols1
diff = [(i, j) for i in range(len(rows)) for j in range(len(cols)) if M0[i][j] != M1[i][j]]
print(f"W = {W}: M entries changed by the tamper: {len(diff)} of {len(rows) * len(cols)}; "
      f"columns touched: {sorted({j for _, j in diff})}")
k, c0, Mpoly = fam.clear(M1)
path = os.path.join(fam.BASE, "data", f"M_sheet{SHEET}_W{W}_tamper.json")
fam.save(path, SHEET, W, k, c0, Mpoly, rows, cols,
         note="TAMPER: I_3 coefficient of X(1,1,1,1) multiplied by 11/10 (sheet 3)")
print(f"saved {os.path.relpath(path, fam.BASE)}; scalar (t^2-1)^{k} * {c0}")
