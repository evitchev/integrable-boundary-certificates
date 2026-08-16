"""N=1 super-OPE engine: free scalar superfields in components.

Letters are triples (p, m, j): parity p (0 boson d^m X^j with m >= 1,
1 fermion d^a psi^j with a >= 0), order, field index.  A monomial is a
tuple of letters in canonical (sorted) order; identical fermion letters
square to zero.  Weights are held doubled (weight2: boson letter 2m,
fermion 2a+1) so all grading stays integer.

Pairings (conventions of the bosonic engine, plus):
  <d^m X^i(z) d^n X^j(w)> = delta (-1)^(m-1) (m+n-1)! (z-w)^-(m+n)
  <d^a psi^i(z) d^b psi^j(w)> = delta (-1)^a (a+b)! (z-w)^-(a+b+1)
Boson-fermion pairings vanish.

Koszul signs enter in exactly two places, kept as small auditable
functions: canonicalization (parity of the odd-letter permutation) and
contraction (word simulation: moving each contracted letter past the
alive odd letters between it and its partner).  The decisive validation
canary is Res(G, G) = 2T exactly (super_validate.py)."""

from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial

from ope_engine import _compositions


def is_odd(letter):
    return letter[0] == 1


def weight2(letter):
    p, m, _ = letter
    return 2 * m + 1 if p else 2 * m


def mono_weight2(mono):
    return sum(weight2(l) for l in mono)


def canonicalize(word):
    """(canonical monomial, sign) or (None, 0) if a fermion repeats.
    Sign = parity of the permutation of the ODD letters from word order
    to sorted order (bosons commute freely)."""
    odds = [l for l in word if is_odd(l)]
    if len(set(odds)) != len(odds):
        return None, 0
    inv = 0
    for i in range(len(odds)):
        for k in range(i + 1, len(odds)):
            if odds[i] > odds[k]:
                inv += 1
    return tuple(sorted(word)), -1 if inv % 2 else 1


def cur_add(a, b, coeff=1):
    for mono, c in b.items():
        v = a.get(mono, 0) + coeff * c
        if v:
            a[mono] = v
        else:
            a.pop(mono, None)
    return a


def d_mono(mono):
    out = {}
    for i, (p, m, j) in enumerate(mono):
        word = mono[:i] + ((p, m + 1, j),) + mono[i + 1:]
        cm, sg = canonicalize(word)
        if cm is not None:
            out[cm] = out.get(cm, 0) + sg
    return {m_: c for m_, c in out.items() if c}


def d_cur(cur):
    out = {}
    for mono, c in cur.items():
        cur_add(out, d_mono(mono), c)
    return out


def _pairing(lp, lq):
    """(numerator, pole) for contracting letters lp(z), lq(w); None if
    they do not pair."""
    p1, m, i = lp
    p2, n, j = lq
    if p1 != p2 or i != j:
        return None
    if p1 == 0:
        return ((-1) ** (m - 1) * factorial(m + n - 1), m + n)
    return ((-1) ** m * factorial(m + n), m + n + 1)


def _koszul(P, Q, psel, qsel):
    """Sign from contracting P[i] with Q[qsel[k]] pairs: simulate the
    word (P letters)(Q letters), moving each contracted P letter past
    the alive odd letters strictly after it in P and strictly before
    its partner in Q.  Order of processing is immaterial; we fix
    increasing P position."""
    aliveP = [True] * len(P)
    aliveQ = [True] * len(Q)
    sign = 1
    for i, jq in sorted(zip(psel, qsel)):
        if is_odd(P[i]):
            cnt = sum(1 for k in range(i + 1, len(P))
                      if aliveP[k] and is_odd(P[k]))
            cnt += sum(1 for k in range(jq)
                       if aliveQ[k] and is_odd(Q[k]))
            if cnt % 2:
                sign = -sign
        aliveP[i] = False
        aliveQ[jq] = False
    return sign


@lru_cache(maxsize=100_000)
def residue_mono(P, Q):
    """Res_{z->w} P(z) Q(w) for canonical super monomials."""
    out = {}
    np_, nq = len(P), len(Q)
    for k in range(1, min(np_, nq) + 1):
        for psel in combinations(range(np_), k):
            for qsel in permutations(range(nq), k):
                num, D = 1, 0
                for a, b in zip(psel, qsel):
                    pr = _pairing(P[a], Q[b])
                    if pr is None:
                        num = 0
                        break
                    num *= pr[0]
                    D += pr[1]
                if not num:
                    continue
                num *= _koszul(P, Q, psel, qsel)
                restP = [P[a] for a in range(np_) if a not in psel]
                restQ = [Q[b] for b in range(nq) if b not in set(qsel)]
                r = len(restP)
                if r == 0:
                    # c-number (z-w)^-D times restQ(w): residue survives
                    # only at the fermionic simple pole D = 1 (bosonic
                    # poles are >= 2), e.g. the zero-mode bracket
                    # Res(psi, psi) = 1.
                    if D == 1:
                        cm, sg = canonicalize(tuple(restQ))
                        if cm is not None:
                            c = out.get(cm, 0) + F(num * sg)
                            if c:
                                out[cm] = c
                            else:
                                out.pop(cm, None)
                    continue
                for ts in _compositions(D - 1, r):
                    den = 1
                    shifted = []
                    for (p, m, j), t in zip(restP, ts):
                        den *= factorial(t)
                        shifted.append((p, m + t, j))
                    word = tuple(shifted) + tuple(restQ)
                    cm, sg = canonicalize(word)
                    if cm is None:
                        continue
                    c = out.get(cm, 0) + F(num * sg, den)
                    if c:
                        out[cm] = c
                    else:
                        out.pop(cm, None)
    return out


def residue_cur(Pc, Qc):
    out = {}
    for mp, cp in Pc.items():
        for mq, cq in Qc.items():
            r = residue_mono(mp, mq)
            if r:
                cur_add(out, r, cp * cq)
    return out


@lru_cache(maxsize=100_000)
def normal_prod_mono(P, Q):
    """Point-split normal product :PQ:(w) = oint dz (z-w)^-1 P(z)Q(w):
    the (z-w)^0 part of the OPE.  Same Wick enumerator as residue_mono
    with Taylor order D instead of D-1, plus the contraction-free
    concatenation term (k = 0).  NOTE: quantum composites like :TT:
    differ from naive concatenation by non-d-exact correction terms —
    conflating the two mislabels charge rays (audit finding, see
    results/super_kdv_calibration.md)."""
    out = {}
    np_, nq = len(P), len(Q)
    for k in range(0, min(np_, nq) + 1):
        for psel in combinations(range(np_), k):
            for qsel in permutations(range(nq), k):
                num, D = 1, 0
                for a, b in zip(psel, qsel):
                    pr = _pairing(P[a], Q[b])
                    if pr is None:
                        num = 0
                        break
                    num *= pr[0]
                    D += pr[1]
                if not num:
                    continue
                num *= _koszul(P, Q, psel, qsel)
                restP = [P[a] for a in range(np_) if a not in psel]
                restQ = [Q[b] for b in range(nq) if b not in set(qsel)]
                r = len(restP)
                if r == 0:
                    if D == 0:
                        cm, sg = canonicalize(tuple(restQ))
                        if cm is not None:
                            c = out.get(cm, 0) + F(num * sg)
                            if c:
                                out[cm] = c
                            else:
                                out.pop(cm, None)
                    continue
                for ts in _compositions(D, r):
                    den = 1
                    shifted = []
                    for (p, m, j), t in zip(restP, ts):
                        den *= factorial(t)
                        shifted.append((p, m + t, j))
                    cm, sg = canonicalize(tuple(shifted) + tuple(restQ))
                    if cm is None:
                        continue
                    c = out.get(cm, 0) + F(num * sg, den)
                    if c:
                        out[cm] = c
                    else:
                        out.pop(cm, None)
    return out


def normal_prod(A, B):
    """Point-split :AB: for currents."""
    out = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            cur_add(out, normal_prod_mono(ma, mb), ca * cb)
    return {m: c for m, c in out.items() if c}


# ------------------------------------------------------------- sectors

def gen_basis_super(w2, nfields):
    """Canonical monomials of doubled weight w2: letters chosen in
    non-decreasing canonical order (so multisets appear once), equal
    consecutive fermion letters forbidden (square zero)."""
    out = []

    def rec2(rem, min_l, acc):
        if rem == 0:
            if acc:
                out.append(tuple(acc))
            return
        for p in (0, 1):
            for j in range(nfields):
                for m in ((range(1, rem // 2 + 1)) if p == 0
                          else range(0, (rem + 1) // 2)):
                    l = (p, m, j)
                    if l < min_l:
                        continue
                    lw = weight2(l)
                    if lw > rem:
                        continue
                    if p == 1 and acc and acc[-1] == l:
                        continue
                    rec2(rem - lw, l, acc + [l])

    rec2(w2, (0, 0, -1), [])
    return sorted(set(out))


def genuine_kernel_super(P_cur, cands, nfields):
    """Genuine charges in span(cands) whose graded bracket with
    oint P_cur is a total derivative.  `cands` is a list of currents
    (dicts), each weight-homogeneous of the same doubled weight, so
    joint kernels can be computed by feeding a previous kernel back in.

    Same single augmented sparse echelon as
    sparse_linalg.genuine_kernel_sparse: target-sector d-columns are
    pre-inserted tagless (value indices shifted above the candidate
    tags), each candidate residue carries an identity tag, and a
    tag-only reduction means the combination's residue is d-exact.
    Genuineness = adds a pivot beyond the candidate-sector d-span.

    Note parity comes for free: a doubled weight w2 is odd iff the
    monomial has an odd number of fermion letters, so every sector is
    parity-homogeneous (odd w2 = fermionic charges of half-integer
    spin)."""
    from sparse_linalg import SparsePivots, _to_int_vec
    w2p = mono_weight2(next(iter(P_cur)))
    w2c = mono_weight2(next(iter(cands[0])))
    nc = len(cands)
    w2t = w2p + w2c - 2
    tgt_index = {m: i for i, m in enumerate(gen_basis_super(w2t, nfields))}
    apiv = SparsePivots()
    for b in gen_basis_super(w2t - 2, nfields):
        col = {}
        for mono, c in d_mono(b).items():
            k = tgt_index[mono] + nc
            col[k] = col.get(k, 0) + c
        col = {k: c for k, c in col.items() if c}
        if col:
            apiv.insert(col)
    kernel_tags, scales = [], []
    for j, cand in enumerate(cands):
        iv, L = _to_int_vec(residue_cur(P_cur, cand), tgt_index)
        scales.append(L)
        aug = {k + nc: c for k, c in iv.items()}
        aug[j] = 1
        r = apiv.reduce(aug)
        if r and max(r) < nc:
            kernel_tags.append(r)
        elif r:
            apiv.rows[max(r)] = r
    cs_index = {m: i for i, m in enumerate(gen_basis_super(w2c, nfields))}
    span = SparsePivots()
    for b in gen_basis_super(w2c - 2, nfields):
        col = {}
        for mono, c in d_mono(b).items():
            col[cs_index[mono]] = col.get(cs_index[mono], 0) + c
        col = {k: c for k, c in col.items() if c}
        if col:
            span.insert(col)
    genuine = []
    for t in kernel_tags:
        cur = {}
        for j, c in t.items():
            cur_add(cur, cands[j], F(c) * scales[j])
        cur = {m: c for m, c in cur.items() if c}
        iv, _ = _to_int_vec(cur, cs_index)
        if span.insert(iv):
            genuine.append(cur)
    return genuine


def q_dim_super(w2, nfields):
    """dim of the charge space at doubled density weight w2: densities
    mod total derivatives."""
    from sparse_linalg import SparsePivots
    basis = gen_basis_super(w2, nfields)
    index = {m: i for i, m in enumerate(basis)}
    piv = SparsePivots()
    for b in gen_basis_super(w2 - 2, nfields):
        col = {}
        for mono, c in d_mono(b).items():
            col[index[mono]] = col.get(index[mono], 0) + c
        col = {k: c for k, c in col.items() if c}
        if col:
            piv.insert(col)
    return len(basis) - len(piv.rows)


class SuperSector:
    """Densities of doubled weight w2 modulo total derivatives."""
    _cache = {}

    def __new__(cls, w2, nfields):
        key = (w2, nfields)
        if key in cls._cache:
            return cls._cache[key]
        self = super().__new__(cls)
        self.basis = gen_basis_super(w2, nfields)
        self.index = {m: i for i, m in enumerate(self.basis)}
        from sparse_linalg import SparsePivots
        self.piv = SparsePivots()
        for b in gen_basis_super(w2 - 2, nfields):
            col = {}
            for mono, c in d_mono(b).items():
                col[self.index[mono]] = col.get(self.index[mono], 0) + c
            col = {k: c for k, c in col.items() if c}
            if col:
                self.piv.insert(col)
        cls._cache[key] = self
        return self

    def is_total_derivative(self, cur):
        from sparse_linalg import _to_int_vec
        if not cur:
            return True
        iv, _ = _to_int_vec(cur, self.index)
        return not self.piv.reduce(iv)
