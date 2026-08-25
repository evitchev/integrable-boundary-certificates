#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ============================================================================
   triangle_targets.wls -- independent Mathematica replication of the
   "verification triangle" targets T1..T6c of
       wolfram/HANDOFF-SPEC.md

   Written from the hand-off spec and the two papers ONLY (no access to the
   project's Python engines, stored results, or target values).

   Engine design
   -------------
   * Normal-ordered products of composite fields (T1, T2, T6b) are computed
     with K. Thielemans' OPEdefs 3.1 (beta 4), loaded unmodified from
     wolfram/OPEdefs/ (licence in the file header: free redistribution on the
     condition the header is kept and modified code is not distributed).
     OPEdefs' NO[A,B] is exactly the spec's point-split product, and nested
     products are right-to-left as required (NO[A,NO[B,C]]).
   * Densities are then mapped to a commutative "Wick" algebra of monomials in
     fld[name,k] = d^k(name).  This is justified by the identity, verified in
     OPEdefs below, that right-nested point-split products of free-field
     derivatives carry no correction terms (all mutual contractions are
     c-numbers), so NO[a,NO[b,NO[c,...]]] = :abc...: as a Wick monomial.
   * Charges are identified modulo total derivatives exactly: at each weight
     the image of d in the monomial space is computed, and class coordinates
     are obtained from a basis of the annihilator (left null space) of that
     image.  A density is a total derivative iff its class coordinates vanish.
   * Residues against vertex operators e^{alpha.phi} (T3, T4, T5, T6c) use a
     generating-function Wick engine written here from the spec's conventions
     (<d phi d phi> = (z-w)^-2, hence <phi phi> = +log(z-w), and
     d^k phi(z) e^{alpha phi(w)} ~ alpha (-1)^(k-1) (k-1)! (z-w)^-k e^{alpha phi}),
     cross-checked monomial by monomial against OPEdefs (V declared as a
     primary with OPE[J,V] = alpha V/(z-w)) through weight 7.
   * Conservation to first order: Res P(z) V_alpha(w) = (d + alpha d phi) J
     for some differential polynomial J  (twisted total derivative).

   Output: results/wolfram_triangle.json (relative to the repository root,
   i.e. the parent of this script's directory).
   ============================================================================ *)

(* ----------------------------------------------------------------------------
   0. Paths, logging, provenance
   ---------------------------------------------------------------------------- *)

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName, First[$ScriptCommandLine]];
scriptPath = ExpandFileName[scriptPath];
scriptDir = DirectoryName[scriptPath];
repoDir = ParentDirectory[scriptDir];
resultsDir = FileNameJoin[{repoDir, "results"}];
outFile = FileNameJoin[{resultsDir, "wolfram_triangle.json"}];
opedefsDir = FileNameJoin[{scriptDir, "OPEdefs"}];

log[msg__] := Print["[triangle] ", StringJoin[ToString /@ {msg}]];
fail[msg__] := (log["FATAL: ", msg]; Exit[1]);
assert[cond_, msg__] := If[TrueQ[cond], Null, fail["assertion failed: ", msg]];

scriptRel = FileNameJoin[{FileNameTake[scriptDir], FileNameTake[scriptPath]}];
log["script: ", scriptRel];
log["kernel: ", $Version];
scriptSHA = FileHash[scriptPath, "SHA256", "HexString"];
runStamp = DateString[Now, {"ISODateTime"}, TimeZone -> 0] <> "Z";

(* ----------------------------------------------------------------------------
   1. OPEdefs (K. Thielemans), vendored unmodified; fetched if absent
   ---------------------------------------------------------------------------- *)

opedefsURL = "https://people.maths.ox.ac.uk/beem/materials/OPEdefs.zip";
opedefsFiles = FileNameJoin[{opedefsDir, #}] & /@ {"Delta.m", "OPEdefs.m"};
If[! AllTrue[opedefsFiles, FileExistsQ],
  log["OPEdefs not found in ", opedefsDir, "; downloading ", opedefsURL];
  Module[{zip = FileNameJoin[{$TemporaryDirectory, "OPEdefs.zip"}], ex},
    URLDownload[opedefsURL, zip];
    ex = ExtractArchive[zip, $TemporaryDirectory, OverwriteTarget -> True];
    Quiet[CreateDirectory[opedefsDir]];
    Do[CopyFile[FileNameJoin[{$TemporaryDirectory, "OPEdefs", f}],
      FileNameJoin[{opedefsDir, f}], OverwriteTarget -> True],
      {f, {"Delta.m", "OPEdefs.m"}}]]];
assert[AllTrue[opedefsFiles, FileExistsQ], "OPEdefs files unavailable"];
opedefsSHA = FileHash[#, "SHA256", "HexString"] & /@ opedefsFiles;

Get[opedefsFiles[[1]]];   (* Delta` *)
Get[opedefsFiles[[2]]];   (* OPEdefs` *)
assert[NumericQ[OPEdefsVersion], "OPEdefs did not load"];
opedefsVersionString = "OPEdefs " <> ToString[OPEdefsVersion] <> " (beta " <>
  ToString[betaVersion] <> ") by K. Thielemans";
log["loaded ", opedefsVersionString];

(* Basic OPEdefs fields.  J = d phi (one boson); JX, JY1..JY3 = dX, dY_i;
   V = vertex operator e^{alpha phi} (declared LAST so that OPEdefs' standard
   form keeps it innermost/rightmost in nested NO's). *)
Bosonic[J, JX, JY1, JY2, JY3, V];
OPE[J, J] = MakeOPE[{One, 0}];
OPE[JX, JX] = MakeOPE[{One, 0}];
OPE[JY1, JY1] = MakeOPE[{One, 0}];
OPE[JY2, JY2] = MakeOPE[{One, 0}];
OPE[JY3, JY3] = MakeOPE[{One, 0}];
OPE[J, V] = MakeOPE[{alpha V}];      (* d phi(z) e^{alpha phi(w)} ~ alpha V/(z-w) *)

engineNO = "Wolfram Language " <> $Version <> " / " <> opedefsVersionString <>
  " (vendored unmodified, wolfram/OPEdefs/) for point-split normal-ordered products" <>
  "; total-derivative quotient by exact linear algebra (this script)";
engineVertex = "Wolfram Language " <> $Version <>
  " / vertex-operator residues by this script's generating-function Wick engine" <>
  " (conventions of the hand-off spec), cross-checked against " <> opedefsVersionString <>
  " through weight 7; twisted-derivative quotient by exact linear algebra";

(* ----------------------------------------------------------------------------
   2. Commutative density algebra:  fld[name, k]  <->  d^k name  (k >= 1)
   ---------------------------------------------------------------------------- *)

(* exponent vector of a single monomial term (coefficient allowed) *)
monoExps[t_] := Module[{vl = Union[Cases[{t}, fld[_, _], Infinity]]},
  If[vl === {}, {{}, {}}, {vl, First[CoefficientRules[t, vl]][[1]]}]];
weightOf[t_] := Module[{ve = monoExps[t]}, Total[MapThread[#1[[2]] #2 &, ve]]];
degreeOf[t_] := Total[Last[monoExps[t]]];

(* all monomials of weight w in the named fields (w >= 0) *)
basis[0, names_List] := {1};
basis[w_Integer?Positive, names_List] := basis[w, names] =
  Union[Flatten[Table[
    Times @@ MapThread[fld[#2, #1] &, {part, tup}],
    {part, IntegerPartitions[w]}, {tup, Tuples[names, Length[part]]}]]];

(* total derivative *)
Dtot[p_] := Module[{vs = Union[Cases[{p}, fld[_, _], Infinity]]},
  Expand[Sum[D[p, v] fld[v[[1]], v[[2]] + 1], {v, vs}]]];

(* robust exact zero test for rational / algebraic / polynomial-in-Q coefficients *)
zeroQ[e_] := Module[{s = Together[Expand[e]]},
  s === 0 || (NumericQ[s] && RootReduce[s] === 0) || Simplify[s] === 0];
allZeroQ[v_List] := AllTrue[v, zeroQ];

(* coordinates of a homogeneous polynomial in a monomial basis *)
coordVector[p_, B_List] := Module[{pe = Expand[p], vl, rules, monos, assoc, vec},
  If[pe === 0, Return[ConstantArray[0, Length[B]]]];
  vl = Union[Cases[B, fld[_, _], Infinity]];
  assert[Complement[Union[Cases[{pe}, fld[_, _], Infinity]], vl] === {},
    "polynomial contains a field variable outside the basis: ", pe];
  rules = CoefficientRules[pe, vl];
  monos = (Times @@ (vl^#)) & /@ rules[[All, 1]];
  assoc = AssociationThread[monos -> rules[[All, 2]]];
  assert[Complement[Keys[assoc], B] === {},
    "polynomial has monomials outside the basis (weight mismatch?)"];
  Lookup[assoc, B, 0]];

(* derivative matrix  d : V_{w-1} -> V_w  (columns = images of basis monomials) *)
dMatrix[w_, names_] := dMatrix[w, names] = Module[{Bw = basis[w, names], Bw1 = basis[w - 1, names]},
  Transpose[coordVector[Dtot[#], Bw] & /@ Bw1]];

(* annihilator of Im d at weight w: rows n with n . (d m) = 0 for all m *)
classAnn[w_, names_] := classAnn[w, names] = Module[{M = dMatrix[w, names], ann},
  ann = NullSpace[Transpose[M]];
  assert[w == 1 || Length[ann] == Length[basis[w, names]] - Length[basis[w - 1, names]],
    "d is not injective at weight ", w - 1, " -- unexpected"];
  ann];

(* class coordinates of a weight-w density: zero iff total derivative *)
classCoords[p_, w_, names_] := classAnn[w, names] . coordVector[p, basis[w, names]];
totalDerivativeQ[p_, w_, names_] := allZeroQ[classCoords[p, w, names]];

(* ----------------------------------------------------------------------------
   3. OPEdefs operator expressions  ->  Wick polynomials
   ---------------------------------------------------------------------------- *)

opeFieldName = <|J -> phi, JX -> X, JY1 -> Y1, JY2 -> Y2, JY3 -> Y3|>;
basicOpQ[f_Symbol] := KeyExistsQ[opeFieldName, f];
basicOpQ[Derivative[_][f_Symbol]] := KeyExistsQ[opeFieldName, f];
basicOpQ[_] := False;

(* vertex operator inside an NO: V -> 1, d^n V -> (d + alpha d phi)^n 1 *)
vertexPoly[V] := 1;
vertexPoly[Derivative[n_][V]] := Nest[Dtot[#] + vAlpha fld[phi, 1] # &, 1, n];

opFactorToPoly[One] := 1;
opFactorToPoly[V] := vertexPoly[V];
opFactorToPoly[Derivative[n_][V]] := vertexPoly[Derivative[n][V]];
opFactorToPoly[f_Symbol?basicOpQ] := fld[opeFieldName[f], 1];
opFactorToPoly[Derivative[n_][f_Symbol?basicOpQ]] := fld[opeFieldName[f], n + 1];
opFactorToPoly[NO[a_, b_]] := (assert[basicOpQ[a],
    "NO not right-nested on a basic field: ", NO[a, b]];
  opFactorToPoly[a] opFactorToPoly[b]);
opFactorToPoly[x_] := fail["unrecognised operator factor: ", x];

operatorFactorQ[f_] := ! FreeQ[f, NO | One | V | Alternatives @@ Keys[opeFieldName]];

termToPoly[t_] := Module[{facs, ops, cs},
  facs = If[Head[t] === Times, List @@ t, {t}];
  ops = Select[facs, operatorFactorQ];
  cs = Select[facs, ! operatorFactorQ[#] &];
  assert[Length[ops] <= 1, "operator factors multiplied with Times: ", t];
  (Times @@ cs) If[ops === {}, 1, opFactorToPoly[First[ops]]]];

opToPoly[expr_] := Module[{e = Expand[expr]},
  If[e === 0, 0, If[Head[e] === Plus, Total[termToPoly /@ (List @@ e)], termToPoly[e]]]];

(* ----------------------------------------------------------------------------
   4. Vertex-operator residues (own Wick engine) and the twisted derivative
   ---------------------------------------------------------------------------- *)

(* Res_{z->w} P(z) e^{alpha name(w)} =: R . e^{alpha name}(w).
   For a Wick monomial P = prod_i d^{k_i} f_i (z):  each factor either contracts
   with the exponential, giving alpha*g*(-1)^(k-1)(k-1)!/(z-w)^k (only if f_i is
   the field `name`; g is the propagator normalisation, g = 1 in the spec), or is
   Taylor-expanded about w.  R is the coefficient of (z-w)^-1 of the product. *)
vertexResidue[P_, name_, alpha_, g_: 1] := Module[{pe = Expand[P], terms},
  If[pe === 0, Return[0]];
  terms = If[Head[pe] === Plus, List @@ pe, {pe}];
  Expand[Total[monomialResidue[#, name, alpha, g] & /@ terms]]];

(* truncated product of polynomials in eps given as coefficient lists (index n+1 <-> eps^n) *)
truncConv[a_List, b_List, deg_] := Table[Expand[Sum[a[[i + 1]] b[[n - i + 1]],
  {i, Max[0, n - Length[b] + 1], Min[n, Length[a] - 1]}]], {n, 0, deg}];

monomialResidue[t_, name_, alpha_, g_] := Module[{facs, coeff, flist, w, deg, lists, prod},
  facs = If[Head[t] === Times, List @@ t, {t}];
  coeff = Times @@ Select[facs, FreeQ[#, fld] &];
  flist = Flatten[Replace[Select[facs, ! FreeQ[#, fld] &],
    {Power[f_fld, e_] :> ConstantArray[f, e], f_fld :> {f}}, {1}]];
  w = Total[weightOf /@ flist];
  If[w == 0, Return[0]];
  deg = w - 1;
  (* (z-w)^k * [d^k f (z) acting on e^{alpha name}(w)] as a polynomial in eps = z-w:
     eps^0: contraction constant (only for the field `name`);  eps^n (n >= k): d^n f / (n-k)! *)
  lists = Function[f, Table[Which[
      n == 0, If[f[[1]] === name, alpha g (-1)^(f[[2]] - 1) (f[[2]] - 1)!, 0],
      n < f[[2]], 0,
      True, fld[f[[1]], n]/(n - f[[2]])!], {n, 0, deg}]] /@ flist;
  prod = Fold[truncConv[#1, #2, deg] &, lists];
  coeff prod[[deg + 1]]];

twistD[p_, name_, alpha_] := Expand[Dtot[p] + alpha fld[name, 1] p];

(* matrix of the twisted derivative D_alpha : V_{w-1} -> V_w *)
twistMatrix[w_, names_, name_, alpha_] := Module[{Bw = basis[w, names], Bw1 = basis[w - 1, names]},
  Transpose[coordVector[twistD[#, name, alpha], Bw] & /@ Bw1]];

(* annihilator rows of Im D_alpha at weight w (exact, over the field generated by alpha) *)
twistAnn[w_, names_, name_, alpha_] := twistAnn[w, names, name, alpha] =
  NullSpace[Transpose[twistMatrix[w, names, name, alpha]]];

(* exact rank of a matrix with entries in Q(Sqrt[d]) (d = 2 or 3) by restriction of scalars:
   a + b Sqrt[d]  ->  [[a, d b], [b, a]];  rank_Q = 2 rank_{Q(Sqrt[d])}.  Rational matrices: direct. *)
algCoords[x_, d_] := Module[{y = ToNumberField[x, Sqrt[d]]},
  Which[Head[y] === AlgebraicNumber, y[[2]], NumberQ[y] && Element[y, Rationals], {y, 0},
    True, fail["entry not in Q(Sqrt[", d, "]): ", x]]];
rankAlg[M_] := Module[{d, A, B, co},
  If[AllTrue[Flatten[M], Element[#, Rationals] &], Return[MatrixRank[M]]];
  d = Which[! FreeQ[M, Sqrt[2]], 2, ! FreeQ[M, Sqrt[3]], 3, True, fail["unsupported field in rankAlg"]];
  co = Map[algCoords[#, d] &, M, {2}];
  A = co[[All, All, 1]]; B = co[[All, All, 2]];
  assert[AllTrue[Flatten[{A, B}], Element[#, Rationals] &], "restriction of scalars failed"];
  Module[{r = MatrixRank[ArrayFlatten[{{A, d B}, {B, A}}]]},
    assert[EvenQ[r], "odd rank after restriction of scalars"]; r/2]];

(* first-order conservation of oint P (weight w) under e^{alpha name}:
   residue lies in Im(d + alpha d name) *)
conservedQ[P_, w_, names_, name_, alpha_, g_: 1] := Module[{R, ann},
  R = vertexResidue[P, name, alpha, g];
  ann = twistAnn[w - 1, names, name, alpha];
  allZeroQ[ann . coordVector[R, basis[w - 1, names]]]];

(* ----------------------------------------------------------------------------
   5. Validation of the machinery on identities derived independently
   ---------------------------------------------------------------------------- *)

log["--- validation ---"];
one = {phi};
u[k_] := fld[phi, k];

(* (a) Virasoro algebra of T_Q with c = 1 - 12 Q^2 (spec convention) *)
TQ = NO[J, J]/2 + Q J';
ttOPE = OPESimplify[OPE[TQ, TQ]];
assert[zeroQ[OPEPole[4][ttOPE] - (1 - 12 Q^2)/2 One], "TT pole 4 (central charge)"];
assert[zeroQ[OPEPole[3][ttOPE]], "TT pole 3"];
assert[OPESimplify[OPEPole[2][ttOPE] - 2 TQ] === 0, "TT pole 2"];
assert[OPESimplify[OPEPole[1][ttOPE] - TQ'] === 0, "TT pole 1"];
log["(a) T_Q T_Q OPE: c = 1 - 12 Q^2, poles 3..1 as required -- OK"];

(* (b) vertex operator weight h = alpha(alpha - 2Q)/2 under T_Q (spec), and
       pole 1 = dV = alpha :J V: *)
tvOPE = OPESimplify[OPE[TQ, V]];
assert[zeroQ[(OPEPole[2][tvOPE] /. V -> 1) - alpha (alpha - 2 Q)/2], "h(V_alpha)"];
assert[OPESimplify[OPEPole[1][tvOPE] - alpha NO[J, V]] === 0, "T V pole 1"];
log["(b) T_Q V_alpha: h = alpha(alpha-2Q)/2 -- OK"];

(* (c) nested point-split products of free-field derivatives are Wick monomials:
       derived by hand: NO[NO[J,J],J] = NO[J,NO[J,J]] + J''  (left nesting carries
       a correction), while right-nested products reorder without corrections. *)
assert[OPESimplify[NO[NO[J, J], J] - NO[J, NO[J, J]] - J''] === 0, "left-nesting correction"];
assert[OPESimplify[NO[J, NO[J', J]] - NO[J', NO[J, J]]] === 0, "reordering NO[J,NO[J',J]]"];
assert[OPESimplify[NO[J'', NO[J, NO[J', J]]] - NO[J, NO[J', NO[J'', J]]]] === 0, "reordering 3"];
assert[OPESimplify[NO[J, NO[J, NO[J, J]]] - NO[J, NO[J, NO[J, NO[J, One]]]]] === 0, "NO with One"];
(* NO[T_0,T_0] = 1/4 J^4 + 1/2 J'' J exactly (hand derivation: double contraction
   4*(1/(z-w)^2)*:J(z)J(w): contributes 2 :J''J: at order 0, times 1/4) *)
assert[zeroQ[opToPoly[NO[NO[J, J]/2, NO[J, J]/2]] - (u[1]^4/4 + u[3] u[1]/2)], "NO[T0,T0]"];
log["(c) nesting/reordering identities and NO[T0,T0] = J^4/4 + J''J/2 -- OK"];

(* (d) quotient mod d: J''J + J'J' is a total derivative, J'J' is not; and
       dimension of classes = p(w) - p(w-1) *)
assert[totalDerivativeQ[u[3] u[1] + u[2]^2, 4, one], "d(J'J) total derivative"];
assert[! totalDerivativeQ[u[2]^2, 4, one], "J'J' not a total derivative"];
assert[Length[classAnn[6, one]] == PartitionsP[6] - PartitionsP[5], "class count at weight 6"];
log["(d) total-derivative quotient -- OK"];

(* (e) own vertex-residue engine vs OPEdefs on every monomial through weight 7,
       symbolic alpha.  OPEdefs: OPEPole[1][P, V] with P right-nested. *)
toNO[m_] := Module[{facs},
  facs = Flatten[Replace[If[Head[m] === Times, List @@ m, {m}],
    {Power[fld[phi, k_], e_] :> ConstantArray[k, e], fld[phi, k_] :> {k}}, {1}]];
  Fold[NO[Derivative[#2 - 1][J], #1] &, One, Reverse[facs]]];
(* P = NO[J^(k1), NO[J^(k2), ... NO[J^(kn), One]]] (right-nested Wick monomial) *)
vAlpha = alpha;   (* momentum used when OPEdefs returns derivatives of V *)
Module[{cnt = 0},
  Do[With[{m = mono},
    Module[{ope, own},
      ope = OPESimplify[OPEPole[1][toNO[m], V]];
      own = vertexResidue[m, phi, alpha];
      assert[zeroQ[opToPoly[ope] - own], "vertex residue mismatch at ", m];
      cnt++]], {w, 1, 7}, {mono, basis[w, one]}];
  log["(e) vertex residues: own engine == OPEdefs on ", cnt, " monomials (weights 1..7) -- OK"]];

(* (f) residue of a total derivative vanishes identically (Res d_z S(z) V(w) = 0) *)
Do[assert[zeroQ[vertexResidue[Dtot[s], phi, alpha]], "Res of total derivative at ", s],
  {s, basis[5, one]}];
log["(f) Res of total derivatives vanishes -- OK"];

(* (g) screening: oint :T_Q T_Q: is conserved under V_alpha when h_Q(alpha) = 1
       (a screening charge commutes with the whole Virasoro vertex algebra).
       At Q = 0: alpha^2 = 2.  Check with alpha = Sqrt[2] and a control alpha = Sqrt[3]. *)
tt0 = opToPoly[NO[NO[J, J]/2, NO[J, J]/2]];
assert[conservedQ[tt0, 4, one, phi, Sqrt[2]], "KdV I3 conserved under e^{sqrt2 phi}"];
assert[conservedQ[tt0, 4, one, phi, -Sqrt[2]], "KdV I3 conserved under e^{-sqrt2 phi}"];
assert[! conservedQ[tt0, 4, one, phi, Sqrt[3]], "KdV I3 not conserved under e^{sqrt3 phi} (control)"];
log["(g) screening test for oint :T_0 T_0: -- OK"];

(* (h) Feigin-Fuchs screening: at Q = 1/2, alpha with alpha(alpha-1)/2 = 1 -> alpha = 2 or -1:
       oint :T_Q T_Q: conserved; control alpha = 3 not. *)
ttQ = opToPoly[NO[TQ, TQ]] /. Q -> 1/2;
assert[conservedQ[ttQ, 4, one, phi, 2], "FF screening alpha=2 at Q=1/2"];
assert[conservedQ[ttQ, 4, one, phi, -1], "FF screening alpha=-1 at Q=1/2"];
assert[! conservedQ[ttQ, 4, one, phi, 3], "control alpha=3 at Q=1/2"];
log["(h) Feigin-Fuchs screening test at Q = 1/2 -- OK"];

log["--- validation complete ---"];

(* ----------------------------------------------------------------------------
   6. Targets
   ---------------------------------------------------------------------------- *)

results = <||>;
record[id_, value_, notes_, engine_] := (results[id] = <|"value" -> value, "exact" -> True,
  "notes" -> notes, "engine" -> engine, "script_sha256" -> scriptSHA, "run" -> runStamp|>;
  log[id, " = ", ToString[value, InputForm]]);
inp[e_] := ToString[e, InputForm];

(* ---------- T1 ---------- *)
Module[{lhs, rhs, cc, ok, notes},
  lhs = opToPoly[NO[TQ, TQ]];
  rhs = (u[1]^4 + (4 Q^2 - 2) u[2]^2)/4;
  cc = Expand /@ classCoords[lhs - rhs, 4, one];
  ok = allZeroQ[cc];
  notes = "T_Q = NO[J,J]/2 + Q J' with J = d phi (OPEdefs); :T_Q T_Q: = NO[T_Q,T_Q] (point-split). " <>
    "Class coordinates of NO[T_Q,T_Q] - (1/4)[(d phi)^4 + (4Q^2-2)(d2 phi)^2] modulo total derivatives, " <>
    "computed with symbolic Q: " <> inp[cc] <> " (all identically zero in Q). " <>
    "Explicit NO[T_Q,T_Q] as a Wick polynomial: " <> inp[Collect[lhs, Q]] <> ". " <>
    "Monomials (d phi)^4, (d2 phi)^2 are read as Wick-normal-ordered monomials, equal to right-nested " <>
    "point-split products of the basic field (identity verified in validation (c)).";
  record["T1", ok, notes, engineNO]];

(* ---------- T2 ---------- *)
Module[{expr, P6, lhsQ, ccQ, oks, notes, qs, ccFac},
  expr = NO[TQ, NO[TQ, TQ]] + 7/8 NO[TQ', TQ'];
  P6 = u[1]^6 + 5/2 u[1]^2 u[2]^2 + 31/24 u[3]^2;
  lhsQ = opToPoly[expr];                       (* symbolic Q *)
  ccQ = Together /@ classCoords[lhsQ - P6/8, 6, one];
  ccFac = Factor /@ ccQ;
  qs = {3/(2 Sqrt[2]), -3/(2 Sqrt[2])};
  oks = allZeroQ[ccQ /. Q -> #] & /@ qs;
  notes = "Left side NO[T_Q, NO[T_Q,T_Q]] + (7/8) NO[T_Q', T_Q'] (OPEdefs, right-nested), T_Q' = derivative " <>
    "of the composite (Leibniz).  Evaluated at both roots Q = +-3/(2 Sqrt[2]) of Q^2 = 9/8: " <> inp[oks] <>
    ".  Class coordinates of (lhs - P6/8) modulo d as functions of symbolic Q (factored): " <>
    inp[ccFac] <> " -- vanish exactly on 8Q^2 - 9 = 0.  " <>
    "P6 = (d phi)^6 + (5/2)(d phi)^2 (d2 phi)^2 + (31/24)(d3 phi)^2 read as Wick monomials.";
  record["T2", And @@ oks, notes, engineNO]];

(* ---------- T3 ---------- *)
Module[{a = 2 Sqrt[2], eps, res, Jd, ann, cR, cJ, rSol, r, consistent, notes, evenPart, oddPart, parityEvenTD},
  (* e^{a phi(z)} e^{-a phi(w)} = (z-w)^(-a^2) :exp(a Sum_k (z-w)^k/k! d^k phi(w)):  (<phi phi> = +log(z-w),
     which is the sign fixed by <d phi d phi> = +(z-w)^-2).  a^2 = 8: residue = coefficient of (z-w)^7. *)
  res = Expand[SeriesCoefficient[Series[Exp[a Sum[eps^k/k! u[k], {k, 1, 7}]], {eps, 0, 7}], 7]];
  Jd = 64/7 u[1]^7 - 40 u[1]^3 u[2]^2 + u[1] u[3]^2;
  ann = classAnn[7, one];
  cR = ann . coordVector[res, basis[7, one]];
  cJ = ann . coordVector[Jd, basis[7, one]];
  assert[! allZeroQ[cJ], "J_density is a total derivative?!"];
  (* solve cR == r cJ: r from the first nonzero component, then check all *)
  rSol = Module[{i = SelectFirst[Range[Length[cJ]], ! zeroQ[cJ[[#]]] &]}, RootReduce[cR[[i]]/cJ[[i]]]];
  consistent = allZeroQ[cR - rSol cJ];
  (* parity split: even part (even number of factors) must be a total derivative *)
  evenPart = Total[Select[List @@ res, EvenQ[degreeOf[#]] &]];
  oddPart = Expand[res - evenPart];
  parityEvenTD = totalDerivativeQ[evenPart, 7, one];
  r = If[consistent, rSol, Null];
  notes = "Residue computed from the explicit free-field formula e^{a phi(z)} e^{-a phi(w)} = (z-w)^(-a^2) " <>
    ":exp(a Sum_{k>=1} (z-w)^k/k! d^k phi(w)): with a = 2 Sqrt[2] (a^2 = 8, eighth-order pole; " <>
    "residue = (z-w)^7 coefficient of the normal-ordered exponential), the sign <phi(z) phi(w)> = +log(z-w) " <>
    "being forced by <d phi d phi> = +(z-w)^-2 and by the spec's weight formula h = alpha^2/2 at Q = 0 " <>
    "(checked in validation (b) via OPEdefs).  The charge class of the residue is compared with that of " <>
    "J_density modulo total derivatives; proportionality holds for the whole residue: " <> inp[consistent] <>
    ".  Residue (Wick polynomial): " <> inp[res] <> ".  Parity-even part of the residue is a total derivative: " <>
    inp[parityEvenTD] <> ".  Ambiguity: the ordering e^{+a}(z) e^{-a}(w) is taken literally from the spec; " <>
    "the opposite ordering would flip the sign of r.  Class coordinates: residue " <> inp[cR] <> ", J " <> inp[cJ] <> ".";
  record["T3", If[r === Null, Null, inp[r]], notes, engineVertex]];

(* ---------- T4 ---------- *)
Module[{P, R, ann, cond, condFac, roots, alpha2, notes, Pk, condK, condKFac},
  P = u[1]^4 + 5/2 u[2]^2;
  R = vertexResidue[P, phi, alpha];
  ann = twistAnn[3, one, phi, alpha];          (* symbolic alpha; D_alpha injective for all alpha *)
  assert[Length[ann] == Length[basis[3, one]] - Length[basis[2, one]], "twisted annihilator size"];
  cond = Together /@ (ann . coordVector[R, basis[3, one]]);
  cond = DeleteCases[Numerator /@ cond, 0];
  condFac = Factor /@ cond;
  (* common roots of all condition polynomials *)
  roots = Select[Union[Flatten[alpha /. Solve[# == 0, alpha] & /@ cond]],
    Function[rt, AllTrue[cond, zeroQ[# /. alpha -> rt] &]]];
  alpha2 = Union[RootReduce[#^2] & /@ DeleteCases[roots, 0]];
  (* generic spin-3 ray for the notes *)
  Pk = u[1]^4 + kappa u[2]^2;
  condK = Together /@ (ann . coordVector[vertexResidue[Pk, phi, alpha], basis[3, one]]);
  condKFac = Factor /@ DeleteCases[Numerator /@ condK, 0];
  notes = "Residue of [(d phi)^4 + (5/2)(d2 phi)^2](z) e^{alpha phi(w)} (symbolic alpha, Q = 0) reduced modulo " <>
    "the twisted image Im(d + alpha d phi) at weight 3 via a basis of the annihilator (the twisted derivative " <>
    "V_2 -> V_3 is injective for every alpha, so no parameter-dependent pivots).  Condition polynomial(s): " <>
    inp[condFac] <> ".  Full root set in alpha: " <> inp[roots] <> ".  VALUE = sorted set of alpha^2 over the " <>
    "NONZERO roots (the spec names {1/2, 8} and mentions alpha = 0 parenthetically; the trivial root alpha = 0 " <>
    "is recorded here and excluded from the value).  For the generic ray (d phi)^4 + kappa (d2 phi)^2 the " <>
    "condition is " <> inp[condKFac] <> ".";
  record["T4", inp /@ alpha2, notes, engineVertex]];

(* ---------- T5 ---------- *)
jointKernelDim[w_, names_, name_, alphas_List, g_: 1] := Module[{Bw, Bw1, anns, Rmats, comb, rank, kdim, Md},
  Bw = basis[w, names]; Bw1 = basis[w - 1, names];
  anns = twistAnn[w - 1, names, name, #] & /@ alphas;
  Rmats = Function[al, Transpose[coordVector[vertexResidue[#, name, al, g], Bw1] & /@ Bw]] /@ alphas;
  comb = Join @@ MapThread[If[#1 === {}, {}, #1 . #2] &, {anns, Rmats}];
  If[comb === {}, Return[{Length[Bw] - Length[Bw1], Length[Bw], Length[Bw1]}]];
  (* consistency: total derivatives lie in the kernel *)
  Md = dMatrix[w, names];
  assert[allZeroQ[Flatten[comb . Md]], "Im d not in the joint kernel at weight ", w];
  rank = rankAlg[comb];
  kdim = Length[Bw] - rank;
  {kdim - Length[Bw1], Length[Bw], Length[Bw1]}];

Module[{dims, dimsRescaled, notes, ok, ws = Range[2, 14], t1, t2},
  {t1, dims} = AbsoluteTiming[Table[First[jointKernelDim[w, one, phi, {Sqrt[2], -Sqrt[2]}]], {w, ws}]];
  log["T5 over Q(Sqrt[2]) took ", t1, " s: ", dims];
  (* independent arithmetic path: rescale phi -> phi/alpha, i.e. propagator g = alpha^2 = 2 and unit
     momentum (+1, -1); everything rational. Kernel dimension is invariant under the field rescaling. *)
  {t2, dimsRescaled} = AbsoluteTiming[Table[First[jointKernelDim[w, one, phi, {1, -1}, 2]], {w, ws}]];
  log["T5 rescaled rational path took ", t2, " s: ", dimsRescaled];
  ok = dims === dimsRescaled;
  assert[ok, "T5: algebraic and rescaled-rational paths disagree"];
  notes = "For each density weight w = 2..14: the space of classes [P] in V_w / d V_{w-1} (V_w = all Wick monomials " <>
    "of weight w in d^k phi, both parities) such that Res P(z) e^{+-Sqrt[2] phi(w)} lies in Im(d +- Sqrt[2] d phi). " <>
    "Computed as dim ker(combined annihilator . residue matrices on V_w) - dim V_{w-1} (d is injective and Im d lies " <>
    "in the kernel, which is asserted).  Dimensions over Q(Sqrt[2]).  Cross-checked by an all-rational path " <>
    "(boson rescaled so that the propagator is 2/(z-w)^2 and the momenta are +-1; identical dimensions): " <>
    inp[ok] <> ".  dim V_w = p(w) (partitions).";
  record["T5", dims, notes, engineVertex]];

(* ---------- T6a ---------- *)
cylCoeffs[sol_, N_, s_] := Switch[sol,
  1, {(2 + N + s)/3, 6, 1, (-9 - 3 N - s)/6, -2},
  2, {(117 + 2 N + N^2 + (15 + N) s)/192, (11 + N + s)/8, 1, -(409 + 46 N + N^2 + (N - 5) s)/192, -2},
  3, {1, 7 - N + s, 1, (13 - 27 N + 2 N^2 + (5 - 2 N) s)/6, -2}];
(* coefficient order: (dX)^4, (dX)^2 (dY.dY), (dY.dY)^2, (d2X)^2, (d2Y.d2Y) *)
Module[{c1, c2, ok, notes},
  c1 = Simplify /@ cylCoeffs[1, 28, -9];
  c2 = Simplify /@ cylCoeffs[2, 28, 9];
  ok = And @@ MapThread[zeroQ[#1 - #2] &, {c1, c2}];
  notes = "Coefficients in the monomial order [(dX)^4, (dX)^2 (dY.dY), (dY.dY)^2, (d2X)^2, (d2Y.d2Y)] taken from " <>
    "the spec's Sol 1 / Sol 2 formulas: Sol 1 at (N,s) = (28,-9): " <> inp[c1] <> "; Sol 2 at (28,9): " <> inp[c2] <>
    ".  Both points satisfy s^2 = (N-25)(N-1): " <> inp[(-9)^2 == (28 - 25) (28 - 1) && 9^2 == (28 - 25) (28 - 1)] <>
    ".  The Y-part (dY.dY)^2 - 2 (d2Y.d2Y) is as stated in the spec's normalisation rule.";
  record["T6a", ok, notes, "Wolfram Language " <> $Version <> " (exact rational arithmetic on the spec's coefficient formulas)"]];

(* ---------- T6b ---------- *)
cylDensity[N_, s_, ynames_List] := Module[{c = cylCoeffs[1, N, s], yy, y22},
  yy = Total[fld[#, 1]^2 & /@ ynames]; y22 = Total[fld[#, 2]^2 & /@ ynames];
  Expand[c[[1]] fld[X, 1]^4 + c[[2]] fld[X, 1]^2 yy + c[[3]] yy^2 + c[[4]] fld[X, 2]^2 + c[[5]] y22]];
Module[{checks, notes, ok, details = {}},
  checks = Table[
    Module[{ynames, yops, TX, TY, comp, dens, names, diff, td},
      ynames = Take[{Y1, Y2, Y3}, n]; yops = Take[{JY1, JY2, JY3}, n];
      names = Join[{X}, ynames];
      TX = NO[JX, JX]/2; TY = Total[NO[#, #] & /@ yops]/2;
      comp = opToPoly[4 (7 NO[TX, TX] + 6 NO[TX, TY] + NO[TY, TY])];
      dens = cylDensity[28, -9, ynames];
      diff = Expand[comp - dens];
      td = totalDerivativeQ[diff, 4, names];
      AppendTo[details, "N_explicit=" <> ToString[n] <> ": 4[7 NO[TX,TX] + 6 NO[TX,TY] + NO[TY,TY]] = " <>
        inp[comp] <> "; difference from Sol-1 density: " <> inp[diff] <> "; total derivative: " <> inp[td]];
      td], {n, {2, 3}}];
  ok = And @@ checks;
  notes = "T_X = NO[JX,JX]/2, T_Y = Sum_i NO[JY_i,JY_i]/2 (OPEdefs, explicit bosons; NO[TX,TY] computed by OPEdefs, " <>
    "it has no contractions).  Compared with the Sol-1 density at (28,-9) modulo total derivatives in the " <>
    "multi-boson Wick algebra, with N_explicit = 2 and 3 (the spec's permitted simplification; coefficient " <>
    "values at N = 28 kept).  " <> StringRiffle[details, "  "];
  record["T6b", ok, notes, engineNO]];

(* ---------- T6c ---------- *)
Module[{res, notes, pair, details = {}},
  res = Table[
    Module[{ynames = Take[{Y1, Y2, Y3}, n], names, dens},
      names = Join[{X}, ynames];
      dens = cylDensity[28, -9, ynames];
      Table[conservedQ[dens, 4, names, X, al], {al, {Sqrt[2], -Sqrt[2], Sqrt[3], -Sqrt[3]}}]],
    {n, {2, 3}}];
  (* res[[i]] = {+sqrt2, -sqrt2, +sqrt3, -sqrt3} for N_explicit = 2, 3 *)
  assert[res[[1]] === res[[2]], "T6c: N_explicit = 2 and 3 disagree"];
  pair = {res[[1, 1]] && res[[1, 2]], res[[1, 3]] && res[[1, 4]]};
  notes = "Density = Sol 1 at (28,-9) (equal to Sol 2 at (28,9) by T6a), weight 4, fields X, Y_1..Y_n.  " <>
    "Residue against e^{alpha X} contracts only dX-factors (Y's are spectators); conservation iff the residue " <>
    "lies in Im(d + alpha dX) at weight 3 in the multi-boson Wick algebra.  Results per momentum " <>
    "[+Sqrt[2], -Sqrt[2], +Sqrt[3], -Sqrt[3]]: N_explicit=2: " <> inp[res[[1]]] <> "; N_explicit=3: " <> inp[res[[2]]] <>
    ".  VALUE = [conserved under BOTH e^{+Sqrt[2]X} and e^{-Sqrt[2]X}, conserved under BOTH e^{+Sqrt[3]X} and " <>
    "e^{-Sqrt[3]X}] (the density is even in X, so the two signs agree).";
  record["T6c", pair, notes, engineVertex]];

(* ----------------------------------------------------------------------------
   7. Provenance and output
   ---------------------------------------------------------------------------- *)

results["_provenance"] = <|
  "script" -> scriptRel,  (* repository-relative: artifacts are clean at source *)
  "script_sha256" -> scriptSHA,
  "kernel_version" -> $Version,
  "kernel_version_number" -> StringTrim[ToString[NumberForm[$VersionNumber, {6, 1}]]],
  "run_utc" -> runStamp,
  "ope_package" -> opedefsVersionString,
  "ope_package_files" -> AssociationThread[FileNameTake /@ opedefsFiles -> opedefsSHA],
  "ope_package_source" -> opedefsURL,
  "ope_package_license" -> "From the OPEdefs.m header: \"You're free to redistribute this package, on the only condition you keep this header and don't distribute modified code.\" (no warranty; author asks to be contacted by users). Used unmodified.",
  "independence" -> "Independent implementation with specified regression targets: inputs wolfram/HANDOFF-SPEC.md (which names the target quantities and their expected values, pinned from the Python stack before this script ran), draft/enhanced_c1.tex, draft/main.tex. No project code or results were read. This is an independent replication, not a blind prediction.",
  (* hostname and absolute command line deliberately omitted: both are
     privacy-screen tokens, and neither adds provenance beyond the
     script hash, kernel version, and run stamp recorded above. *)
  "ambiguities" -> {
    "Monomials such as (d phi)^4 in the spec/papers are read as Wick-normal-ordered monomials; verified in OPEdefs that these coincide with right-nested point-split products of the basic field derivatives (no correction terms).",
    "Sign of the vertex OPE: d phi(z) e^{alpha phi(w)} ~ +alpha e^{alpha phi}(w)/(z-w) (i.e. <phi phi> = +log(z-w)), fixed by the spec's <d phi d phi> = +(z-w)^-2 and its weight formula h = alpha(alpha-2Q)/2, which OPEdefs reproduces with this sign.",
    "T3: the ordering e^{+a phi}(z) e^{-a phi}(w) of the residue is taken literally; the opposite ordering flips the sign of r.",
    "T4: value = set of alpha^2 over the nonzero roots of the condition polynomial; alpha = 0 is a root too and is recorded in the notes.",
    "T6c: value[i] = conserved under BOTH signs of the momentum; all four individual results are in the notes.",
    "T6b/T6c: N_explicit = 2 (and 3 as a check) Y bosons with the N = 28 coefficient values, as permitted by the spec."}
|>;

Quiet[CreateDirectory[resultsDir]];
jsonString = StringReplace[ExportString[results, "JSON", "Compact" -> False], "\\/" -> "/"];
Export[outFile, jsonString, "Text"];
log["wrote ", outFile, " (", StringLength[jsonString], " bytes)"];
(* verify the file parses as JSON and has every target *)
Module[{back = Import[outFile, "RawJSON"]},
  assert[AssociationQ[back], "output is not valid JSON"];
  assert[SubsetQ[Keys[back], {"T1", "T2", "T3", "T4", "T5", "T6a", "T6b", "T6c", "_provenance"}], "missing target ids"]];
log["done"];
