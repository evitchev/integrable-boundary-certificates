"""Phase-2: pair-class collapse for the invariant (pure-pair) engine.

Exact drop-in for invariant_engine.residue_inv_mono, efficient when the
FIRST argument (the z-side density, e.g. the low charge of a kernel
computation) has few pairs and the second has many identical ones.

Idea: enumerate the small P side raw (subsets of its <= 2*kP slots) and
collapse the large Q side.  A raw matching touches Q only through
"instances": distinct Q pairs receiving one or two contraction
endpoints.  So enumerate
  - subsets S of P slots (|S| = k >= 1);
  - set partitions of S into blocks of size <= 2 (one block = one Q
    instance; <= 76 partitions for k = 6);
  - a Q pair type per instance and a slot arrangement within it,
with multiplicity
  - prod_type falling_factorial(m_type, #instances of that type)
    (instances are distinguishable through their P blocks, so ordered
    assignment of distinct identical-type pairs is exactly a falling
    factorial -- no symmetry division needed), and
  - a factor 2 per instance of a symmetric (u,u) type whose two raw
    slot assignments give identical outcomes.
Then build the index graph over P pairs + instances exactly as the raw
engine does (paths -> template pairs, closed components -> powers of
N), Taylor-distribute the pole budget over free P slots, and append the
untouched Q pairs by type count.  Output format and values identical to
residue_inv_mono; validated by collapsed_pairs_validate.py.
"""

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from math import factorial, comb, perm

import invariant_engine as _ie
from invariant_engine import canon, canon_pair, _compositions

_ORIG = _ie.residue_inv_mono


def _blocks(items):
    """Set partitions of `items` into blocks of size 1 or 2."""
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for sub in _blocks(rest):
        yield [(first,)] + sub
    for i in range(len(rest)):
        pair_block = (first, rest[i])
        remaining = rest[:i] + rest[i + 1:]
        for sub in _blocks(remaining):
            yield [pair_block] + sub


def residue_inv_mono_fast(P, Q):
    kP = len(P)
    PS = [(i, s, P[i][s]) for i in range(kP) for s in (0, 1)]
    qtypes = sorted(Counter(Q).items())          # [(type, count)]
    out = {}

    def emit(num, D, ncyc, templates, freeP, restQ_counts):
        r = len(freeP)
        if r == 0:
            return
        for ts in _compositions(D - 1, r):
            den = 1
            for t in ts:
                den *= factorial(t)
            newpairs = []
            it = iter(range(r))
            shifts = dict(zip(freeP, ts))
            for (o1, id1), (o2, id2) in templates:
                a = o1 + shifts.get(id1, 0)
                b = o2 + shifts.get(id2, 0)
                newpairs.append(canon_pair(a, b))
            mono = canon(tuple(newpairs) +
                         tuple(t for t, c in restQ_counts for _ in range(c)))
            slot = out.setdefault(mono, {})
            slot[ncyc] = slot.get(ncyc, 0) + F(num, den)

    from itertools import combinations, product
    nslots = 2 * kP
    for k in range(1, nslots + 1):
        for S in combinations(range(nslots), k):
            for blocks in _blocks(list(S)):
                nb = len(blocks)
                # assign a Q type to each block, with slot arrangements
                def assign(bi, used, wnum, edges):
                    if bi == nb:
                        mult = 1
                        for t, u in used.items():
                            mult *= perm(dict(qtypes)[t], u)
                        if mult == 0:
                            return
                        finish(wnum * mult, edges, used)
                        return
                    blk = blocks[bi]
                    for (u, v), avail in qtypes:
                        cnt = used.get((u, v), 0)
                        if cnt >= avail:
                            continue
                        inst = (bi, (u, v))
                        if len(blk) == 1:
                            p = blk[0]
                            m = PS[p][2]
                            opts = ([(u, v, 2)] if u == v
                                    else [(u, v, 1), (v, u, 1)])
                            for cslot, fslot, fac in opts:
                                e = [(p, cslot, inst, fslot)]
                                used2 = dict(used)
                                used2[(u, v)] = cnt + 1
                                w = fac * (-1) ** (m - 1) \
                                    * factorial(m + cslot - 1)
                                assign(bi + 1, used2, wnum * w, edges + e)
                        else:
                            p1, p2 = blk
                            m1, m2 = PS[p1][2], PS[p2][2]
                            arr = ([((u, v), 2)] if u == v
                                   else [((u, v), 1), ((v, u), 1)])
                            for (c1, c2), fac in arr:
                                e = [(p1, c1, inst, None),
                                     (p2, c2, inst, None)]
                                used2 = dict(used)
                                used2[(u, v)] = cnt + 1
                                w = fac \
                                    * (-1) ** (m1 - 1) * factorial(m1 + c1 - 1) \
                                    * (-1) ** (m2 - 1) * factorial(m2 + c2 - 1)
                                assign(bi + 1, used2, wnum * w, edges + e)
                    return

                def finish(wnum, edges, used):
                    D = sum(PS[p][2] + c for p, c, _, _ in edges)
                    # graph over P pairs and instances
                    adjP = {i: [] for i in range(kP)}
                    adjI = {}
                    inst_free = {}
                    for p, c, inst, fslot in edges:
                        pi = PS[p][0]
                        adjP[pi].append(inst)
                        adjI.setdefault(inst, []).append(pi)
                        if fslot is not None:
                            inst_free[inst] = fslot
                    Sset = set(S)
                    freeP_slots = {}
                    for idx, (pi, s, o) in enumerate(PS):
                        if idx not in Sset:
                            freeP_slots.setdefault(pi, []).append((o, idx))
                    ncyc = 0
                    templates = []
                    freeP_ids = []
                    seen = set()
                    for start in list(adjP) + list(adjI):
                        if start in seen:
                            continue
                        stack, comp = [start], []
                        while stack:
                            v = stack.pop()
                            if v in seen:
                                continue
                            seen.add(v)
                            comp.append(v)
                            nbrs = adjP[v] if isinstance(v, int) else adjI[v]
                            stack.extend(nbrs)
                        ends = []
                        for v in comp:
                            if isinstance(v, int):
                                for o, idx in freeP_slots.get(v, []):
                                    ends.append((o, ('P', idx)))
                            else:
                                if v in inst_free:
                                    ends.append((inst_free[v], ('Q', v)))
                        if not ends:
                            ncyc += 1
                        else:
                            (o1, id1), (o2, id2) = ends
                            templates.append(((o1, id1), (o2, id2)))
                            for o, ident in ends:
                                if ident[0] == 'P':
                                    freeP_ids.append(ident)
                    restQ = [(t, c - used.get(t, 0)) for t, c in qtypes]
                    restQ = [(t, c) for t, c in restQ if c > 0]
                    emit(wnum, D, ncyc, templates, freeP_ids, restQ)

                assign(0, {}, 1, [])
    return {m: {p: c for p, c in d.items() if c}
            for m, d in out.items() if any(d.values())}


RAW_THRESHOLD = 200_000


def _raw_count(kP, kQ):
    return sum(comb(2 * kP, k) * perm(2 * kQ, k)
               for k in range(min(2 * kP, 2 * kQ) + 1))


@lru_cache(maxsize=4096)
def residue_inv_mono_smart(P, Q):
    if len(P) <= 4 and _raw_count(len(P), len(Q)) > RAW_THRESHOLD:
        return residue_inv_mono_fast(P, Q)
    return _ORIG(P, Q)
