# Originality screen of the opers (DDMST (3.18) and Suzuki citation graphs, FRU 2304.04501)

Commissioned 2026-10-02 (the lead's relay of the operator's go). The claims screened
are those of the write-up draft of 2026-10-01:
- C1: the screening systems and their algebras;
- C2: the Gamma-symbol / three-term opers of non-integer effective order;
- C3: Sol 1 as Suzuki's sl(2|1) equation with a centrifugal term.

A screening agent did the discovery. The cognitive agent then re-read every quote
below in the fetched LaTeX source, checked equation numbers against the compiled
.aux labels, and re-read two of the agent's dismissals.

## Verdicts

- **C3: IN PRINT, and Sol 1 is very probably the quantum AKNS structure.** See
  section 1. This is the finding that changes the write-up.
- **C1: PARTIAL.** The three X + O(N) families are in print as first members
  (Vitchev 2004). Screening-pair-plus-third-exponential commutants in two bosons
  are in print for AKNS and Getmanov. The C_2^(1)/D_3^(2) pairing is known in
  general.
- **C2: NOT FOUND.** No ODE/IM operator with a Gamma-ratio symbol or a continuously
  varying order was found. FRU is the nearest framework, and it is not one.

## 1. Sol 1 and the quantum AKNS structure (C3)

**What is in print.**

- **Fateev-Lukyanov, hep-th/0510271v1, "Boundary RG flow associated with the AKNS
  soliton hierarchy".**
  - The local IM commute with the screenings (82): `{\cal V}_{\pm} =
    \re^{\sqrt{n}\,X_R \pm \ri\sqrt{n+2}\,Y_R}`, and with (87): `{\cal V}_{0} =
    \re^{-{2X_R\over \sqrt{n}}}`.
  - The ODE is (138):

        \bigg[\, -{{\rd^2}\over{\rd x^2}} - p^2+2q\, \kappa\, \re^{x}+\kappa^2\, \big(\re^{2x}+\re^{-nx}\big) \, \bigg]\ \Psi(x) = 0

  - Its WKB "exactly reproduces the vacuum eigenvalues" (95) of the local IM: the
    source paragraph after the top-term display.
  - (95) prints the first three vacuum eigenvalues, I_1^vac, I_2^vac (odd in Q)
    and I_3^vac.
- **Assis, 1211.2397v2.** (2.3), `\label{aknsparamap}`, with β² = 1/(M+1) and n = −2β²,
  maps the FL equation to (2.4):

      \left( - \frac{d^2}{d x^2} + { \frac{\ell(\ell+1)}{x^2} } + {x^{2M}} + {\alpha x^{M-1}} - E \right) \psi(x) =0

- **Bazhanov-Tsuboi, 0805.4274v3.** (4.37) is the same equation with ℓ(ℓ+1)/y².
  (4.60) is the correspondence D_i(E, ±r) = vacuum eigenvalues of the U_q(sl(2|1)^)
  Q-operators. It is introduced as "We expect the following exact
  correspondence", i.e. a conjecture.
- **Attribution.** Suzuki quant-ph/0003066 eq. (1) has NO centrifugal term (the
  file was grepped). The term first appears in Dorey-Dunning-Tateo hep-th/0103051
  eq. (2.1), `\label{sh}`.

**The cognitive agent's three checks that Sol 1 = AKNS at n = −(t−1)/(t+1).**

1. **Parameters.** The record's Sol 1 oper has M = (t+3)/(t−1). Assis's map gives
   n = −2/(M+1) = −(t−1)/(t+1).
2. **Top layer** (`sol1_akns_top_check.py`, log beside it). Take Sol 1's F1 top
   layer (item 178) and FL's printed top terms of the normalised vacuum
   eigenvalues, with n = −(t−1)/(t+1), FL's P² along the record's Y and Q² along X.
   Every consecutive coefficient ratio agrees with scale factor exactly 1, for
   K = 2..5 at t = 7, 11/3 and −2. The "nan" entries are Gamma poles at isolated
   (t, K). With P² and Q² swapped the ratios disagree (7-8 distinct factors), so
   the agreement is not automatic. The same n is found as in check 1, with no
   free parameter.
3. **Parity structure.** FL's I_2^vac is odd in Q, as the record's Sol 1 has
   X-odd charges at odd weight (item 288). The J-even AKNS charges I_1, I_3, ...
   are the cylindrical (X-even) tower.

A screening Gram-matrix comparison also matches at n = −p²/2. But the cognitive
agent could not pin the vertex-operator normalisation, and the same comparison
for Sol 2 against Getmanov matches the Gram matrix while the central charges
disagree. It is therefore NOT offered as evidence.

**Decisive test, for the record.** Compare FL's (95) I_3^vac with the record's Sol 1
spin-3 vacuum eigenvalue, and FL's I_2^vac with Sol 1's X-odd weight-3 charge, at
n = −(t−1)/(t+1), P² along Y and Q² along X with one overall scale each. If they
agree:
- Sol 1 is the quantum AKNS integrable structure, continued to n < 0;
- its oper is FL (138) = Assis (2.4) = BT (4.37), printed with the WKB match
  (FL) and the conjectured spectral identification (BT (4.60));
- what remains new for Sol 1 is the realisation of AKNS by X + an O(N) vector
  through Vir_N (item 285), and its place among the three solutions.

## 2. C1

- **Vitchev, hep-th/0404195v1, lines 749-753:** "There are three non-trivial
  solutions for $I_3$ and $I_5$ with $\mathbb Z_2\times O(N)$, and it seems
  plausible that these are the first members of infinite families of commuting
  charges."
  - The coefficients carry `\sqN` = √((N−25)(N−1)) (macro at line 45), which
    N = (t²−25)/(t²−1) rationalises.
  - No screenings, no algebra, no ODE for these families.
- **FL (82) and (87)**, and Kotousov-Lukyanov 1910.05947v2 (3.33)-(3.34)
  (`\label{AKNSint1}`, the AKNS I_1, I_2, I_3 in J and a Coulomb-gas T with
  Q_L = (1+k)/√k), with Getmanov at (3.43)-(3.45): two-boson commutants of a
  screening pair plus a third exponential.
- **Masoero-Raimondo-Valeri 1511.00895v2**, table `\label{table:langlands}`, lists
  C_n^(1) with D_(n+1)^(2): the Sol 3 pairing in general. The X + O(N)
  realisation is not there.

## 3. C2

- **DDMST hep-th/0612298v2 (3.18)** (`\label{sun0}`, compiled 3.18):
  `\Bigl((-1)^{n}D_{n}({\bf g})-P_K(x,E) \Bigr)\psi(x,E,{\bf g} )=0`, with
  `P_K = ( x^{h^{\vee} M/K}-E)^K` (3.17). Integer order.
- **FRU 2304.04501v1:** Gaudin model and Deligne's category. Complex-order
  pseudo-differential operators with regular singularities only: no spectral term
  E x^n, no CFT IM, no Gamma symbol (the one `\Gamma` in the file is a macro
  definition).
- **BKL 1706.09941 (5.69)-(5.72)**, cited in item 312 as a "non-integer-n
  precedent": there n is a non-integer exponent in a SECOND-ORDER equation, so it
  is not a precedent for non-integer order.

## Corrections to the brief and to item 312

- Suzuki's eq. (1) has no centrifugal term; it is due to DDT (2.1).
- hep-th/0105177 is Bazhanov-Hibberd-Khoroshkin's W_3 / quantum Boussinesq paper,
  not an sl(2|1) paper.
- BKL's non-integer n is an exponent, not an order (above).

## Screen

- INSPIRE citer lists: DDMST (69), Suzuki (40), FL (25), Vitchev (2), triaged by
  title and abstract.
- About 40 fulltext queries.
- 25 arXiv sources fetched and grepped.

The full lists are in the cognitive agent's review package. Limits:
- INSPIRE fulltext is token-sensitive, so LaTeX-style queries may give false
  zeros.
- The math-only literature on non-integer-order opers is not indexed by INSPIRE.
- Only four citer lists were screened.

## Addendum 2026-10-10: the decisive Sol 1 test, done

The comparison asked for in section 1 has been run (`sol1_akns_spin3_check.py`, log beside it).

- **Result.** FL's I_3^vac of quantum AKNS equals the record's certified Sol 1
  spin-3 vacuum eigenvalue EXACTLY, in all six coefficients, at all 29 fibres of
  the record's spin-3 table.
- **Parameters.** One fixed, t-independent dictionary: Q_FL = sqrt2 P_X,
  P_FL = sqrt2 pi, with n = −(t−1)/(t+1). This is the n that Assis's map and the
  top-layer check of section 1 gave.
- **Controls.** The literal printed (95) fails at every fibre, and so does the tamper
  n + 1/10.
- **A typo in FL (95).** The printed (5n+4) term lacks a factor 1/16. With it
  restored, (95) agrees with FL's own general top-term formula (checked at three
  values of n); as printed, it contradicts that formula.

**Conclusion.** Sol 1 is the quantum AKNS integrable structure of Fateev-Lukyanov,
continued to n = −(t−1)/(t+1) < 0 for t > 1. That is established at spin 3 in full,
and at the top layer for spins 3-9.
- Its oper is FL (138) = Assis (2.4) = BT (4.37). The WKB match to the local IM is
  printed by FL; the spectral identification is BT's conjecture (4.60).
- **Not checked.** The full lower layers at spins ≥ 5, which FL does not print.
- **New for Sol 1, given this identification:** the realisation of AKNS by X + an
  O(N) vector through Vir_N, and its place among the three solutions.
- **The paper.** Neither version 1 (deposited) nor version 2 (HEAD on 2026-10-10)
  cites FL, BT, Assis, or DDT hep-th/0103051, which first added the centrifugal
  term.
