# The paperclip spin-7 charge density P_8 (new result, 2026-07-27)

The spin-7 local integral of motion of the paperclip family
(hep-th/0404195 §3.2, hep-th/0312168) is $\mathbb{I}_7 = \oint dz\, P_8$,
with $P_8$ below.  It was not published in either paper (only its vacuum
eigenvalue, eq. (61) of hep-th/0404195, appeared).  Computed by exact
kernel solution of $[\mathbb{I}_3, \mathbb{I}_7] = 0$ at 39 rational
values of $n$, rational-function reconstruction with held-out
verification, and final symbolic verification identically in $n$
(`code/compute_p8.py`, `code/verify_symbolic.py`).

Notation: $(n_1\ldots n_a|m_1\ldots m_b) = \partial^{n_1}X\cdots
\partial^{n_a}X\,\partial^{m_1}Y\cdots\partial^{m_b}Y$, as in
hep-th/0404195 eq. (46); Wick pairing normalization of its eq. (17).
Normalization: coefficient of $(11111111|) = (\partial X)^8$ equals 1.
Common denominator: $15\,n^3(7n+8)(7n+10)(7n+12)$.

| monomial | coefficient |
|---|---|
| (33\|11) | $\dfrac{112(n+2)(61n^4+261n^3+357n^2+166n+24)}{n^2(7n+8)(7n+10)(7n+12)}$ |
| (44\|) | $-\dfrac{8(5890n^6+38318n^5+97580n^4+122647n^3+79086n^2+24696n+2880)}{15n^3(7n+8)(7n+10)(7n+12)}$ |
| (13\|13) | $\dfrac{224(n+2)(131n^4+499n^3+586n^2+212n+24)}{n^2(7n+8)(7n+10)(7n+12)}$ |
| (1133\|) | $\dfrac{112(31n^4+113n^3+135n^2+58n+8)}{n^2(7n+8)(7n+10)}$ |
| (1223\|) | $\dfrac{28(128n^4+496n^3+643n^2+314n+48)}{n^2(7n+8)(7n+10)}$ |
| (\|44) | $-\dfrac{8(5890n^6+32362n^5+67800n^4+67673n^3+33284n^2+7884n+720)}{15n^3(7n+8)(7n+10)(7n+12)}$ |
| (14\|12) | $-\dfrac{224(n+1)(n+2)^2(25n^2+50n+12)}{n^2(7n+8)(7n+10)(7n+12)}$ |
| (\|111122) | $-\dfrac{28(n+2)^2(7n+2)(7n+4)(24n^2+59n+30)}{n^3(7n+8)(7n+10)(7n+12)}$ |
| (\|1133) | $\dfrac{112(n+2)(7n+2)(31n^4+135n^3+201n^2+118n+24)}{n^3(7n+8)(7n+10)(7n+12)}$ |
| (11\|1122) | $-\dfrac{168(n+2)^2(7n+2)(32n^2+67n+20)}{n^2(7n+8)(7n+10)(7n+12)}$ |
| (12\|1112) | $-\dfrac{224(n+2)^2(4n+7)(7n+2)(9n+4)}{n^2(7n+8)(7n+10)(7n+12)}$ |
| (11\|33) | $\dfrac{112(n+2)(61n^4+227n^3+255n^2+82n+8)}{n^2(7n+8)(7n+10)(7n+12)}$ |
| (22\|1111) | $-\dfrac{84(n+2)^2(7n+2)(8n^2+17n+4)}{n^2(7n+8)(7n+10)(7n+12)}$ |
| (22\|22) | $-\dfrac{56(n+2)(332n^4+1228n^3+1435n^2+566n+72)}{n^2(7n+8)(7n+10)(7n+12)}$ |
| (1111\|22) | $-\dfrac{84(n+2)(8n^2+15n+2)}{n(7n+8)(7n+10)}$ |
| (1112\|12) | $-\dfrac{224(n+2)(4n+1)(9n+14)}{n(7n+8)(7n+10)}$ |
| (23\|12) | $-\dfrac{224(n+1)(n+2)(3n+2)(25n^2+50n+12)}{n^2(7n+8)(7n+10)(7n+12)}$ |
| (1122\|11) | $-\dfrac{168(n+2)(32n^2+61n+14)}{n(7n+8)(7n+10)}$ |
| (\|1223) | $\dfrac{28(n+2)(7n+2)(128n^4+528n^3+739n^2+402n+72)}{n^3(7n+8)(7n+10)(7n+12)}$ |
| (111122\|) | $-\dfrac{28(24n^2+37n+8)}{n(7n+8)}$ |
| (\|11111111) | $\dfrac{(n+2)^3(7n+2)(7n+4)(7n+6)}{n^3(7n+8)(7n+10)(7n+12)}$ |
| (11\|111111) | $\dfrac{28(n+2)^3(7n+2)(7n+4)}{n^2(7n+8)(7n+10)(7n+12)}$ |
| (1111\|1111) | $\dfrac{42(n+2)^2(7n+2)}{n(7n+8)(7n+10)}$ |
| (111111\|11) | $\dfrac{28(n+2)}{7n+8}$ |
| (11111111\|) | $1$ |

Exact machine-readable form: `paperclip_P8.json`.

## Structural observations

1. **Duality.** Under $n \to -2-n$ the factor sets $\{7n+2, 7n+4, 7n+6\}$
   and $\{7n+12, 7n+10, 7n+8\}$ swap (up to sign) and $n \leftrightarrow
   -(n+2)$: this is the $X \leftrightarrow Y$ (equivalently $P
   \leftrightarrow Q$) symmetry of the paperclip; e.g. the coefficients of
   $(11111111|)$ and $(|11111111)$ map to each other's reciprocals.
2. **Resonance poles.** The zero/pole divisor is built from $7n + 2j = 0$: poles at
   $j = 4,5,6$ (and $n = 0$), zeros at $j = 1,2,3$ --- together
   exactly the factors appearing in the
   quasiclassical Wronskian coefficient $S_7$, eq. (68) of
   hep-th/0404195 — same pattern as $(5n+2j)$ in $P_6$/$S_5$ and
   $(3n+2j)$ in $P_4$/$S_3$.  The density's analytic structure in $n$
   mirrors the spectral data of the paperclip ODE.

## Erratum to hep-th/0404195 found en route

In eq. (51), the coefficient of $P_6^{[10]} = (111111|)$ is printed as
$-\frac{(6+5n)(8+5n)}{120(2+n)}$; solving $[\mathbb{I}_3,\mathbb{I}_5]=0$
from scratch at $n = 2, 3, 7/3$ shows it must be
$-\frac{n(6+5n)(8+5n)}{120(2+n)}$ (the two agree only at $n=1$).
Also, the right-hand side of eq. (50) is printed in terms of
$P_6^{[i]}$ but must read the spin-4 basis $P_4^{[i]}$ of eq. (47).
