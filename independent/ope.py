"""
Minimal exact free-boson OPE engine (multi-species), conventions of the papers:
  <d^a X_i(z) d^b X_j(w)> = delta_ij (-1)^(a-1) (a+b-1)! (z-w)^-(a+b)
Densities are Wick-ordered polynomials in J_{s,a} := d^a X_s (a>=1), represented as
  monomial = sorted tuple of (species, order); poly = dict monomial -> coeff.
Coefficients may be Fraction, int, sympy expressions, or any field-like object.
"""
from fractions import Fraction as Fr
from collections import Counter, defaultdict
from math import factorial
from functools import lru_cache
import itertools

# ---------------- basic poly ops ----------------
def norm(p):
    """drop zero coefficients"""
    out = {}
    for m, c in p.items():
        if c != 0:
            out[m] = c
    return out

def padd(p, q, scale=1):
    out = dict(p)
    for m, c in q.items():
        out[m] = out.get(m, 0) + scale * c
    return norm(out)

def pscale(p, s):
    return norm({m: s * c for m, c in p.items()})

def pmul(p, q):
    """commutative (classical) product of monomials -- used only for J*U in D_alpha"""
    out = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(sorted(m1 + m2))
            out[m] = out.get(m, 0) + c1 * c2
    return norm(out)

def mono(*orders, species=0):
    return tuple(sorted((species, o) for o in orders))

def M(orders_x=(), orders_y=()):
    return tuple(sorted([(0, o) for o in orders_x] + [(1, o) for o in orders_y]))

def weight(m):
    return sum(o for s, o in m)

def deriv_mono(m):
    """d of a monomial: returns poly"""
    out = {}
    cnt = Counter(m)
    for (s, o), c in cnt.items():
        lst = list(m)
        lst.remove((s, o))
        lst.append((s, o + 1))
        mm = tuple(sorted(lst))
        out[mm] = out.get(mm, 0) + c
    return out

def deriv(p):
    out = {}
    for m, c in p.items():
        for mm, k in deriv_mono(m).items():
            out[mm] = out.get(mm, 0) + c * k
    return norm(out)

def partial(p, s, a):
    """partial derivative wrt the variable J_{s,a}"""
    out = {}
    for m, c in p.items():
        k = m.count((s, a))
        if k:
            lst = list(m)
            lst.remove((s, a))
            mm = tuple(lst)
            out[mm] = out.get(mm, 0) + c * k
    return norm(out)

def euler(p, nspecies=2, maxorder=None):
    """Euler operator per species: E_s(P) = sum_a (-d)^a dP/dJ_{s,a}. Returns dict s -> poly.
    P is a total derivative iff all E_s(P) vanish (weight>0)."""
    if maxorder is None:
        maxorder = max([o for m in p for (s, o) in m], default=0)
    res = {}
    for s in range(nspecies):
        tot = {}
        for a in range(1, maxorder + 1):
            q = partial(p, s, a)
            if not q:
                continue
            for _ in range(a):
                q = deriv(q)
            tot = padd(tot, q, scale=(-1) ** a)
        res[s] = tot
    return res

def is_total_derivative(p, nspecies=2):
    e = euler(p, nspecies)
    return all(len(v) == 0 for v in e.values())

# ---------------- Wick / OPE ----------------
@lru_cache(maxsize=None)
def _taylor(Aprime, m):
    """coefficient of t^m in prod_{(s,a) in Aprime} sum_k t^k/k! J_{s,a+k}; returns tuple of (mono, coeff)"""
    cur = [dict() for _ in range(m + 1)]
    cur[0][()] = Fr(1)
    for (s, a) in Aprime:
        new = [dict() for _ in range(m + 1)]
        for d in range(m + 1):
            for k in range(d + 1):
                src = cur[d - k]
                if not src:
                    continue
                f = Fr(1, factorial(k))
                for mm, c in src.items():
                    nm = tuple(sorted(mm + ((s, a + k),)))
                    new[d][nm] = new[d].get(nm, 0) + c * f
        cur = new
    return tuple(cur[m].items())

def _kmatrices(typesA, typesB):
    """enumerate contraction matrices k[i][j] respecting species and multiplicities"""
    nA = [n for (_, n) in typesA]
    mB = [n for (_, n) in typesB]
    I, J = len(typesA), len(typesB)
    allowed = [[typesA[i][0][0] == typesB[j][0][0] for j in range(J)] for i in range(I)]
    result = []
    def rec_row(i, colcap, rows):
        if i == I:
            result.append([r[:] for r in rows])
            return
        # enumerate row i vector
        row = [0] * J
        def rec_col(j, remaining):
            if j == J:
                rows.append(row[:])
                rec_row(i + 1, colcap, rows)
                rows.pop()
                return
            maxk = min(remaining, colcap[j]) if allowed[i][j] else 0
            for k in range(maxk + 1):
                row[j] = k
                colcap[j] -= k
                rec_col(j + 1, remaining - k)
                colcap[j] += k
            row[j] = 0
        rec_col(0, nA[i])
    rec_row(0, mB[:], [])
    return result

def mono_ope_coeff(A, B, n):
    """coefficient of (z-w)^n in the OPE A(z) B(w) (Wick), as poly at w.
    n = -1: residue (zero mode bracket density); n = 0: normal ordered product A_{(-1)}B."""
    cA = Counter(A); cB = Counter(B)
    typesA = list(cA.items()); typesB = list(cB.items())
    out = {}
    for kmat in _kmatrices(typesA, typesB):
        p = 0
        val = 1
        mult = 1
        r = [0] * len(typesA); c = [0] * len(typesB)
        for i, ((sa, a), na) in enumerate(typesA):
            for j, ((sb, b), mb) in enumerate(typesB):
                k = kmat[i][j]
                if k:
                    p += k * (a + b)
                    val *= ((-1) ** (a - 1) * factorial(a + b - 1)) ** k
                    mult = Fr(mult, factorial(k))
                    r[i] += k; c[j] += k
        m = p + n
        if m < 0:
            continue
        for i, ((sa, a), na) in enumerate(typesA):
            mult *= Fr(factorial(na), factorial(na - r[i]))
        for j, ((sb, b), mb) in enumerate(typesB):
            mult *= Fr(factorial(mb), factorial(mb - c[j]))
        Aprime = []
        for i, ((sa, a), na) in enumerate(typesA):
            Aprime += [(sa, a)] * (na - r[i])
        Bprime = []
        for j, ((sb, b), mb) in enumerate(typesB):
            Bprime += [(sb, b)] * (mb - c[j])
        Aprime = tuple(sorted(Aprime))
        coef = mult * val
        for mm, cc in _taylor(Aprime, m):
            nm = tuple(sorted(mm + tuple(Bprime)))
            out[nm] = out.get(nm, 0) + coef * cc
    return norm(out)

def ope_coeff(P, Q, n):
    out = {}
    for A, ca in P.items():
        for B, cb in Q.items():
            for mm, c in mono_ope_coeff(A, B, n).items():
                out[mm] = out.get(mm, 0) + ca * cb * c
    return norm(out)

def bracket_density(P, Q):
    """density of [oint P, oint Q] = oint Res P(z) Q(w)"""
    return ope_coeff(P, Q, -1)

def nprod(P, Q):
    """normal ordered product :PQ: = P_{(-1)}Q"""
    return ope_coeff(P, Q, 0)

# ---------------- vertex operators ----------------
def mono_vertex_residue(A, alpha):
    """S such that Res_{z->w} A(z) e^{alpha.X}(w) = :S e^{alpha.X}:(w).
    alpha: dict species -> value (any ring element)."""
    cA = Counter(A)
    typesA = list(cA.items())
    out = {}
    ranges = [range(na + 1) for (_, na) in typesA]
    for ks in itertools.product(*ranges):
        p = 0
        coef = 1
        Aprime = []
        for ((s, a), na), k in zip(typesA, ks):
            if k:
                p += k * a
                v = alpha[s] * ((-1) ** (a - 1) * factorial(a - 1))
                coef = coef * (v ** k) * (factorial(na) // (factorial(k) * factorial(na - k)))
            Aprime += [(s, a)] * (na - k)
        m = p - 1
        if m < 0:
            continue
        Aprime = tuple(sorted(Aprime))
        for mm, cc in _taylor(Aprime, m):
            out[mm] = out.get(mm, 0) + coef * cc
    return norm(out)

def vertex_residue(P, alpha):
    out = {}
    for A, ca in P.items():
        for mm, c in mono_vertex_residue(A, alpha).items():
            out[mm] = out.get(mm, 0) + ca * c
    return norm(out)

def D_alpha(U, alpha):
    """twisted derivative: d(U) + (alpha . J) U"""
    out = deriv(U)
    for s, a in alpha.items():
        if a != 0:
            out = padd(out, pmul({((s, 1),): 1}, U), scale=a)
    return norm(out)

# ---------------- monomial enumeration ----------------
def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest

def monomials(w, nspecies=1, parity=None):
    """all monomials of total weight w. parity: None (all) or tuple of required parities per species (0=even count)"""
    res = []
    if nspecies == 1:
        for part in partitions(w):
            if parity is None or len(part) % 2 == parity[0]:
                res.append(tuple(sorted((0, o) for o in part)))
        return res
    if nspecies == 2:
        for wx in range(w + 1):
            for px in partitions(wx):
                if parity is not None and len(px) % 2 != parity[0]:
                    continue
                for py in partitions(w - wx):
                    if parity is not None and len(py) % 2 != parity[1]:
                        continue
                    res.append(tuple(sorted([(0, o) for o in px] + [(1, o) for o in py])))
        return res
    raise NotImplementedError

# ---------------- exact linear algebra over a field ----------------
def rref(rows, ncols, zero=lambda x: x == 0):
    """rows: list of lists (field elements). returns (rref_rows, pivot_cols)"""
    A = [r[:] for r in rows]
    pivots = []
    r = 0
    for c in range(ncols):
        pr = None
        for i in range(r, len(A)):
            if not zero(A[i][c]):
                pr = i; break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        inv = 1 / A[r][c] if not isinstance(A[r][c], int) else Fr(1, A[r][c])
        A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and not zero(A[i][c]):
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        pivots.append(c)
        r += 1
        if r == len(A):
            break
    return A[:r], pivots

def rank(rows, ncols, zero=lambda x: x == 0):
    return len(rref(rows, ncols, zero)[1])

def nullspace(rows, ncols, zero=lambda x: x == 0):
    """kernel of the matrix (rows x ncols) acting on column vectors: returns basis vectors (lists)"""
    R, piv = rref(rows, ncols, zero)
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for f in free:
        v = [0] * ncols
        v[f] = 1
        for i, pc in enumerate(piv):
            v[pc] = -R[i][f]
        basis.append(v)
    return basis

def polys_to_matrix(polys, index=None):
    """given list of polys, return (matrix rows = polys as coordinate vectors, index list)"""
    if index is None:
        index = sorted({m for p in polys for m in p})
    pos = {m: i for i, m in enumerate(index)}
    rows = []
    for p in polys:
        r = [0] * len(index)
        for m, c in p.items():
            r[pos[m]] = c
        rows.append(r)
    return rows, index

def euler_vec(p, nspecies, index=None):
    """flatten euler images into one poly with species-tagged monomials (tag by shifting species by 10*s)"""
    e = euler(p, nspecies)
    out = {}
    for s, q in e.items():
        for m, c in q.items():
            out[(s,) + m] = c
    return out

# ---------------- generic charge-space helpers ----------------
def class_dim(w, nspecies=1, parity=None):
    """dim Q_{w-1} = #monomials(w) - #monomials(w-1) in the sector"""
    return len(monomials(w, nspecies, parity)) - len(monomials(w - 1, nspecies, parity))

def kernel_of_ad(P, w, nspecies=1, parity=None, zero=lambda x: x == 0):
    """genuine kernel of X -> [oint P, oint X] on charge classes of density weight w.
    returns (dim, list of E-images (flattened polys) spanning the kernel classes, list of representative densities)"""
    mons = monomials(w, nspecies, parity)
    imgs = [euler_vec(bracket_density(P, {m: 1}), nspecies) for m in mons]
    rows, index = polys_to_matrix(imgs)
    # columns = monomials; we need kernel of map mono-coeff-vector -> image, i.e. nullspace of transpose
    T = [[rows[j][i] for j in range(len(mons))] for i in range(len(index))]
    ns = nullspace(T, len(mons), zero) if index else [[1 if i == j else 0 for i in range(len(mons))] for j in range(len(mons))]
    # densities in the kernel (including total derivatives)
    dens = []
    for v in ns:
        d = {}
        for c, m in zip(v, mons):
            if c != 0:
                d[m] = c
        dens.append(d)
    # genuine classes: E-images of dens, take independent ones
    evecs = [euler_vec(d, nspecies) for d in dens]
    erows, eindex = polys_to_matrix(evecs)
    R, piv = rref(erows, len(eindex), zero) if eindex else ([], [])
    dim = len(piv)
    # pick representatives: rows of erows that are independent
    reps = []
    chosen = []
    cur = []
    for d, ev in zip(dens, erows):
        trial = cur + [ev]
        if rank(trial, len(eindex), zero) > len(cur):
            cur = trial
            reps.append(d)
    return dim, reps

def vertex_kernel(alpha, w, nspecies=1, parity=None, zero=lambda x: x == 0, field_one=1):
    """subspace of densities (weight w, sector) whose charges are first-order conserved under e^{alpha.X}.
    Returns (genuine dim, list of representative densities)."""
    mons = monomials(w, nspecies, parity)
    S = [vertex_residue({m: field_one}, alpha) for m in mons]
    Us = monomials(w - 2, nspecies, None) if w >= 3 else ([()] if w == 2 else [])
    Ds = [D_alpha({u: field_one}, alpha) for u in Us]
    allp = S + Ds
    rows, index = polys_to_matrix(allp)
    ncols = len(mons) + len(Us)
    T = [[rows[j][i] for j in range(ncols)] for i in range(len(index))]
    ns = nullspace(T, ncols, zero) if index else [[1 if i == j else 0 for i in range(ncols)] for j in range(ncols)]
    dens = []
    for v in ns:
        d = {}
        for c, m in zip(v[:len(mons)], mons):
            if not zero(c):
                d[m] = c
        if d:
            dens.append(d)
    # genuine dim: dim of span(dens) - dim Im d
    if not dens:
        return 0, [], mons
    drows, dindex = polys_to_matrix(dens, mons)
    dimspan = rank(drows, len(mons), zero)
    dim_imd = len(monomials(w - 1, nspecies, parity))
    return dimspan - dim_imd, dens, mons

def intersect_spaces(list_of_dens_lists, mons, zero=lambda x: x == 0):
    """intersection of subspaces of span(mons) given by generating sets; returns basis (list of coordinate vectors)"""
    # represent each subspace by its annihilator; intersection = kernel of stacked annihilators
    n = len(mons)
    stacked = []
    for dens in list_of_dens_lists:
        rows, _ = polys_to_matrix(dens, mons)
        # annihilator of row space: nullspace of rows (as matrix acting on column vectors)
        ann = nullspace(rows, n, zero)
        stacked += ann
    if not stacked:
        return [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    return nullspace(stacked, n, zero)
