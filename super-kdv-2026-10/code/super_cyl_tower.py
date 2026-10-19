"""Track T2 sampler: the super-KdV tower of the lambda family at a RATIONAL
lambda = l, and its cylinder eigenvalues on the Fock highest-weight state
(NS and R), via super_cyl_vev.py.  Exact (Fractions).

  I3(l) = :T_l T_l: - 1/4 :G_l dG_l:   (point-split; super_lambda.py, kappa = -1/4)
  I_{2k-1}(l), k >= 3: the joint (ad_{G_l}, ad_{I3(l)}) kernel at doubled
       weight w2 = 4k (asserted one-dimensional mod d), normalised so that the
       coefficient of (dX)^{2k} is 1/2^k (the :T^k: normalisation; the
       eigenvalue then starts (P^2/2)^k).

Eigenvalues are independent of the total-derivative representative
(super_cyl_vev_validate.py check (v)).  Output: one JSON per l with
{w2: {"NS": {exp: "num/den"}, "R": {...}}} under
results/lab/susy4/tanabe/samples/.  Usage: super_cyl_tower.py MAXW2 l1 [l2 ...]
(l as 'a/b')."""
import json
import os
import sys
import time
from fractions import Fraction as F

from super_engine import (cur_add, d_cur, gen_basis_super,
                          genuine_kernel_super, normal_prod)
import super_cyl_vev as scv

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, os.pardir, 'results', 'lab', 'susy4', 'tanabe', 'samples')
DX, D2X, PSI, DPSI = (0, 1, 0), (0, 2, 0), (1, 0, 0), (1, 1, 0)


def T_lf(l):
    return {(DX, DX): F(1, 2), (D2X,): l, (PSI, DPSI): F(-1, 2)}


def G_lf(l):
    return {(DX, PSI): F(1), (DPSI,): 2 * l}


def I3_of(l):
    T, G = T_lf(l), G_lf(l)
    cur = dict(normal_prod(T, T))
    cur_add(cur, normal_prod(G, d_cur(G)), F(-1, 4))
    return {m: c for m, c in cur.items() if c}


def normalise(cur, k):
    key = tuple([DX] * (2 * k))
    c = cur.get(key, 0)
    if not c:
        raise RuntimeError(f"(dX)^{2*k} coefficient vanishes")
    s = F(1, 2 ** k) / c
    return {m: v * s for m, v in cur.items()}


def tower(l, maxw2):
    Gl = G_lf(l)
    I3 = I3_of(l)
    out = {8: I3}
    for w2 in range(12, maxw2 + 1, 4):
        t0 = time.time()
        greps = genuine_kernel_super(Gl, [{b: F(1)} for b in gen_basis_super(w2, 1)], 1)
        joint = genuine_kernel_super(I3, greps, 1)
        if len(joint) != 1:
            raise RuntimeError(f"l={l} w2={w2}: joint kernel dim {len(joint)} != 1")
        out[w2] = normalise(joint[0], w2 // 4)
        print(f"  l={l} w2={w2}: |G-kernel|={len(greps)}, joint dim 1, "
              f"{len(out[w2])} monomials, {time.time()-t0:.1f}s", flush=True)
    return out


def main():
    maxw2 = int(sys.argv[1])
    ls = [F(x) for x in sys.argv[2:]]
    os.makedirs(OUTDIR, exist_ok=True)
    for l in ls:
        t0 = time.time()
        tw = tower(l, maxw2)
        res = {}
        for w2, dens in tw.items():
            res[str(w2)] = {}
            for sec in ('NS', 'R'):
                t1 = time.time()
                poly = scv.vev_fraction_poly(dens, sec)
                res[str(w2)][sec] = {str(e): str(v) for e, v in sorted(poly.items())}
                print(f"  l={l} w2={w2} {sec}: eigenvalue {len(poly)} terms, {time.time()-t1:.1f}s", flush=True)
            res[str(w2)]["density"] = {repr(m): str(c) for m, c in sorted(dens.items())}
        path = os.path.join(OUTDIR, f"l_{l.numerator}_{l.denominator}_w{maxw2}.json")
        json.dump({"l": str(l), "maxw2": maxw2, "eig": res}, open(path, "w"), indent=1)
        print(f"l={l}: done in {time.time()-t0:.1f}s -> {path}", flush=True)


if __name__ == "__main__":
    main()
