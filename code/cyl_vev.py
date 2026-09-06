"""Vacuum expectation values of the cylindrical integrals of motion on the
momentum state |P, Q>, X-momentum P, Y-momentum vector Q (O(N) invariant:
only Q^2 enters), by the recipe VALIDATED on hep-th/0404195 eqs. (59)-(60)
(paperclip I_3, I_5 reproduced coefficient by coefficient, 2026-09-02):

  * read the plane density polynomial on the cylinder: d_z^k X -> d_w^{k-1} J;
  * cylinder normal ordering on |P,Q>: Wick sum over partial pairings of the
    letters; an unpaired J gives its zero mode (P or Q_i), an unpaired
    derivative of J gives 0; a contracted pair (d^a J, d^b J) gives
    C(a,b) = (-1)^b R^{(a+b)}(0), R(eps) = 1/(4 sinh^2(eps/2)) - 1/eps^2;
  * O(N): Y letters come in contracted pairs (k,l) = d^k Y . d^l Y; after a
    Wick pairing, each connected component of pairs is a closed cycle
    (factor N) or an open chain whose two end letters must be zero modes
    (factor Q^2);
  * the article's momenta are P_art = i P, Q_art = i Q (P^2 -> -P^2).

The result is returned in the ARTICLE convention.  Exact (Fractions / sympy).
"""
from fractions import Fraction as F
import sympy as sp

_eps = sp.Symbol('eps')
_R = sp.series(1 / (4 * sp.sinh(_eps / 2) ** 2) - 1 / _eps ** 2, _eps, 0, 24).removeO()
_Rd = [F(str(sp.diff(_R, _eps, k).subs(_eps, 0))) for k in range(22)]

def C(a, b):
    return (-1) ** b * _Rd[a + b]

def _pairings(letters):
    """all partial pairings of a list of letter ids (species-compatible pairs only)"""
    def rec(remaining):
        if not remaining:
            yield []
            return
        first, rest = remaining[0], remaining[1:]
        # first unpaired
        for p in rec(rest):
            yield p
        for i, other in enumerate(rest):
            if letters[first][0] == letters[other][0]:
                for p in rec(rest[:i] + rest[i + 1:]):
                    yield [(first, other)] + p
    return rec(list(range(len(letters))))

def vev_monomial(xs, ypairs, N, P2, Q2):
    """eigenvalue of the cylinder-normal-ordered monomial prod d^{k} X prod (d^k Y . d^l Y)
    on |P,Q>, as a polynomial value with P^2 -> P2, Q^2 -> Q2 (our momentum convention)."""
    letters = [('X', k - 1, None) for k in xs]                      # (species, J-derivative order, pair id)
    for pid, (k, l) in enumerate(ypairs):
        letters.append(('Y', k - 1, pid)); letters.append(('Y', l - 1, pid))
    total = F(0)
    for pairing in _pairings(letters):
        paired = {i for pr in pairing for i in pr}
        val = F(1)
        for a, b in pairing:
            val *= C(letters[a][1], letters[b][1])
        # unpaired X letters: zero modes
        nx = 0
        ok = True
        for i, (sp_, j, pid) in enumerate(letters):
            if i in paired or sp_ != 'X':
                continue
            if j != 0:
                ok = False; break
            nx += 1
        if not ok or nx % 2:
            continue
        val *= P2 ** (nx // 2)
        # Y structure: union-find over pair ids via cross-pair contractions
        npairs = len(ypairs)
        parent = list(range(npairs))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for a, b in pairing:
            if letters[a][0] == 'Y':
                pa, pb = find(letters[a][2]), find(letters[b][2])
                if pa != pb:
                    parent[pa] = pb
        comp_unpaired = {}
        for i, (sp_, j, pid) in enumerate(letters):
            if sp_ != 'Y':
                continue
            r = find(pid)
            comp_unpaired.setdefault(r, [])
            if i not in paired:
                comp_unpaired[r].append(j)
        for r, ups in comp_unpaired.items():
            if len(ups) == 0:
                val *= N                         # closed cycle: free index sum
            elif len(ups) == 2:
                if any(j != 0 for j in ups):
                    val = F(0); break
                val *= Q2                        # open chain: Q_i Q_i
            else:
                raise RuntimeError("odd component")
        total += val
    return total

def vev_density(dens, N, P2, Q2):
    """dens: {(xs, ypairs): coeff}; returns the eigenvalue in OUR convention (P2 = P^2, Q2 = Q^2)."""
    return sum(F(c) * vev_monomial(xs, yp, N, P2, Q2) for (xs, yp), c in dens.items())

def vev_coefficients(dens, N, W):
    """Exact coefficients {(a, b): Fraction} of P2^a Q2^b, a + b <= W/2, of the eigenvalue (OUR
    convention): evaluation points are chosen greedily until the exact evaluation matrix has full
    rank, the square system is solved over the rationals, and the result is re-checked on further
    independent points (raises if the eigenvalue is not a polynomial of this degree)."""
    monos = [(a, b) for a in range(W // 2 + 1) for b in range(W // 2 + 1 - a)]
    n = len(monos)
    pool = [(F(i), F(j)) for i in range(1, W // 2 + 4) for j in range(1, W // 2 + 4)]
    pool += [(F(-1), F(2)), (F(3, 2), F(-1, 3)), (F(-2), F(-3)), (F(1, 2), F(5))]
    rows, rhs, echelon, used = [], [], [], []
    for p2, q2 in pool:
        row = [p2 ** a * q2 ** b for a, b in monos]
        # reduce against the echelon basis to test independence
        r = list(row)
        for piv, e in echelon:
            if r[piv] != 0:
                f = r[piv] / e[piv]
                r = [x - f * y for x, y in zip(r, e)]
        piv = next((i for i, x in enumerate(r) if x != 0), None)
        if piv is None:
            continue
        echelon.append((piv, r))
        rows.append(row); rhs.append(F(vev_density(dens, N, p2, q2))); used.append((p2, q2))
        if len(rows) == n:
            break
    if len(rows) < n:
        raise ValueError("could not find a full-rank evaluation set")
    M = [r + [v] for r, v in zip(rows, rhs)]
    for col in range(n):
        piv = next(i for i in range(col, n) if M[i][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        inv = 1 / M[col][col]
        M[col] = [x * inv for x in M[col]]
        for i in range(n):
            if i != col and M[i][col] != 0:
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[col])]
    coeffs = {monos[i]: M[i][n] for i in range(n)}
    check = [pt for pt in pool if pt not in used][:4]
    for p2, q2 in check:
        if sum(c * p2 ** a * q2 ** b for (a, b), c in coeffs.items()) != F(vev_density(dens, N, p2, q2)):
            raise ValueError("eigenvalue is not a polynomial of the assumed degree")
    return coeffs


def vev_polynomial(dens, N, W):
    """symbolic polynomial in the ARTICLE momenta (P^2 -> -P^2, Q^2 -> -Q^2), exact"""
    P, Q = sp.symbols('P Q')
    return sp.expand(sum(sp.Rational(c.numerator, c.denominator) * (-P ** 2) ** a * (-Q ** 2) ** b
                         for (a, b), c in vev_coefficients(dens, N, W).items()))
