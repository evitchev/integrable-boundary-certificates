#!/usr/bin/env wolframscript
(* Lab opener, Mathematica side: exact multiplicity stratification of the
   cylindrical curve at one weight.

   Input:  wolfram/lab/cyl_strata_w{W}_sol{k}.m  (integer polynomial matrices
           in t of ad_{I_3} in class coordinates; built and interpolation-
           verified by cyl_strata_build.py).
   Method: generic rank r0 (agreement at several rational t); a full-rank
           r0 x r0 submatrix chosen at a generic t; its determinant D(t) is a
           maximal minor, so EVERY rank-drop point is a root of D (complete
           candidate set); D is factored over Q and, for each irreducible
           factor, the rank of M at one root is computed exactly (Galois
           conjugates share the rank); nullity - r0-nullity = multiplicity
           jump.  Points with no drop are reported as rejected candidates.
   Output: wolfram/lab/cyl_strata_w{W}.json with provenance.
   Status: LAB artifact -- hypotheses until certified in the open stack.  *)

w = ToExpression[$ScriptCommandLine[[2]]];
lab = DirectoryName[$InputFileName];
root = FileNameJoin[{lab, "..", ".."}];
out = <||>;
probes = {7/3, 11/5, -13/7, 29/11};

curveN[t_] := (t^2 - 25)/(t^2 - 1);
curveS[t_] := -24 t/(t^2 - 1);

rankAt[M_, t0_] := MatrixRank[M /. t -> t0];   (* rational points only *)

(* rank of M over the number field Q[t]/(f), f irreducible: Gaussian
   elimination with entries reduced mod f and pivot inverses from the
   extended Euclidean algorithm -- no algebraic numbers, no RootReduce *)
(* finite-field screen: rank over GF(p)[t]/(f).  Rank can only DROP under
   reduction mod p, so rank_p == r0 PROVES no drop over Q(t0); only
   candidates with rank_p < r0 go to the exact rational elimination. *)
(* Every pivot must be a UNIT of Z_p[t]/(f) (gcd with f constant): f may
   split mod p, and with unit pivots the elimination is valid in every
   component field; a unit minor mod (p, f) lifts, f being primitive and
   irreducible over Q (Gauss), to a nonzero minor in Q[t]/(f).  A column
   with nonzero but no unit entries aborts the screen ($Failed -> exact
   fallback).  Also require p not dividing the leading coefficient of f. *)
rankModFp[M_, f_, p_] := Module[{A, nr, nc, r = 0, c, i, piv, g, ab, inv, k, unitQ},
  If[Mod[Coefficient[f, t, Exponent[f, t]], p] == 0, Return[$Failed]];
  A = Map[PolynomialRemainder[#, f, t, Modulus -> p] &, M, {2}];
  unitQ[e_] := e =!= 0 && Exponent[PolynomialGCD[e, f, Modulus -> p], t] == 0;
  {nr, nc} = Dimensions[A];
  For[c = 1, c <= nc && r < nr, c++,
    piv = 0;
    For[i = r + 1, i <= nr, i++, If[unitQ[A[[i, c]]], piv = i; Break[]]];
    If[piv == 0,
      If[AnyTrue[A[[r + 1 ;;, c]], # =!= 0 &], Return[$Failed]];   (* nonzero non-units only: not a field argument *)
      Continue[]];
    If[piv != r + 1, A[[{r + 1, piv}]] = A[[{piv, r + 1}]]];
    {g, ab} = PolynomialExtendedGCD[A[[r + 1, c]], f, t, Modulus -> p];
    inv = PolynomialRemainder[PolynomialMod[ab[[1]] PowerMod[g, -1, p], p], f, t, Modulus -> p];
    A[[r + 1]] = PolynomialRemainder[Expand[inv #], f, t, Modulus -> p] & /@ A[[r + 1]];
    For[k = r + 2, k <= nr, k++,
      If[A[[k, c]] =!= 0, A[[k]] = PolynomialRemainder[Expand[#], f, t, Modulus -> p] & /@ (A[[k]] - A[[k, c]] A[[r + 1]])]];
    r++];
  r];

rankModF[M_, f_] := Module[{A = Map[PolynomialRemainder[#, f, t] &, M, {2}], nr, nc, r = 0, c, i, piv, g, ab, inv, k},
  {nr, nc} = Dimensions[A];
  For[c = 1, c <= nc && r < nr, c++,
    piv = 0;
    For[i = r + 1, i <= nr, i++, If[A[[i, c]] =!= 0, piv = i; Break[]]];
    If[piv == 0, Continue[]];
    If[piv != r + 1, A[[{r + 1, piv}]] = A[[{piv, r + 1}]]];
    {g, ab} = PolynomialExtendedGCD[A[[r + 1, c]], f, t];
    inv = PolynomialRemainder[ab[[1]]/g, f, t];
    A[[r + 1]] = PolynomialRemainder[Expand[inv #], f, t] & /@ A[[r + 1]];
    For[k = r + 2, k <= nr, k++,
      If[A[[k, c]] =!= 0, A[[k]] = PolynomialRemainder[Expand[#], f, t] & /@ (A[[k]] - A[[k, c]] A[[r + 1]])]];
    r++];
  r];

Do[
  file = FileNameJoin[{lab, "cyl_strata_w" <> ToString[w] <> "_sol" <> ToString[sol] <> ".m"}];
  M = Get[file];
  {nr, nc} = Dimensions[M];
  ranks = MatrixRank[M /. t -> #] & /@ probes;
  r0 = Max[ranks];
  If[Min[ranks] =!= r0, Print["WARNING: probe ranks disagree ", ranks]];
  nullity0 = nc - r0;
  (* choose a full-rank r0 x r0 submatrix at a generic rational point:
     pivot columns of the RREF are independent columns; same on the
     transpose for independent rows *)
  pivotCols[A_] := Module[{rr = RowReduce[A]},
     Function[row, First[Select[Range[Length[row]], row[[#]] =!= 0 &]]] /@ Select[rr, AnyTrue[#, Function[x, x =!= 0]] &]];
  Mg = M /. t -> probes[[1]];
  colPiv = pivotCols[Mg];
  rowPiv = pivotCols[Transpose[Mg]];
  sub = M[[rowPiv, colPiv]];
  If[MatrixRank[sub /. t -> probes[[1]]] =!= r0, Print["submatrix selection failed"]; Abort[]];
  tD = AbsoluteTiming[detD = TimeConstrained[Det[sub], 600, $Failed];][[1]];
  If[detD === $Failed,
    (* fallback: evaluation / interpolation with a rigorous degree bound
       (sum over rows of the maximal entry degree), exact rational
       determinants at bound+1 points (avoiding the poles t = +-1),
       verified at 3 held-out points -- fail closed *)
    degBound = Total[Max[Exponent[#, t] & /@ row] & /@ sub];
    pts = Select[Range[2, degBound + 12], # =!= 1 && # =!= -1 &];
    nodes = Take[pts, degBound + 1]; extra = Take[Drop[pts, degBound + 1], 3];
    tD = AbsoluteTiming[
      vals = Det[sub /. t -> #] & /@ nodes;
      detD = Expand[InterpolatingPolynomial[Transpose[{nodes, vals}], t]];][[1]];
    If[AnyTrue[extra, (detD /. t -> #) =!= Det[sub /. t -> #] &],
      Print["weight ", w, " Solution ", sol, ": interpolation verification FAILED"]; Continue[]];
    Print["   (Det by interpolation: degree bound ", degBound, ", ", Length[nodes], " nodes, verified at 3 extra points)"]];
  detD = Expand[detD];
  facs = FactorList[detD];
  found = {}; rejected = {};
  Do[
    f = fac[[1]];
    If[FreeQ[f, t], Continue[]];
    deg = Exponent[f, t];
    t0 = If[deg == 1, t /. First[Solve[f == 0, t]], Root[Function[tt, f /. t -> tt], 1]];
    (* t = +-1 are the poles of the parametrization (N -> infinity): not curve
       points; they appear only because denominators were cleared.  Record,
       do not count.  (N = 1 is t = infinity, checked separately in Python.) *)
    If[deg == 1 && (t0 === 1 || t0 === -1), AppendTo[rejected,
        <|"factor" -> ToString[f, InputForm], "note" -> "parametrization pole t=+-1, not a curve point"|>]; Continue[]];
    rk = If[deg == 1, rankAt[M, t0],
      Module[{rp = rankModFp[M, f, 2147483647]},
        If[rp === r0, r0, rankModF[M, f]]]];   (* GF(p) screen ($Failed or < r0 -> exact fallback) *)
    nul = nc - rk;
    rec = <|"factor" -> ToString[f, InputForm], "degree" -> deg,
            "t0" -> ToString[t0, InputForm],
            "N" -> ToString[RootReduce[curveN[t0]], InputForm],
            "s" -> ToString[RootReduce[curveS[t0]], InputForm],
            "rank" -> rk, "nullity" -> nul, "jump" -> nul - nullity0|>;
    If[rk < r0, AppendTo[found, rec], AppendTo[rejected, rec]],
    {fac, facs}];
  Print["weight ", w, " Solution ", sol, ": ", nr, "x", nc, ", generic rank ", r0,
        " (generic multiplicity ", nullity0, "), Det time ", Round[tD, 0.01], " s, deg D = ", Exponent[detD, t]];
  Do[Print["   JUMP at t0 = ", r["t0"], "  (N,s) = (", r["N"], ", ", r["s"], ")  nullity ", r["nullity"], " (jump +", r["jump"], ")"], {r, found}];
  Print["   rejected candidate factors: ", Length[rejected]];
  out["sol" <> ToString[sol]] = <|"rows" -> nr, "cols" -> nc, "generic_rank" -> r0,
     "generic_multiplicity" -> nullity0, "det_degree" -> Exponent[detD, t],
     "jumps" -> found, "rejected_candidates" -> rejected|>,
  {sol, 1, 3}];

out["provenance"] = <|"weight" -> w, "engine" -> $Version,
  "script_sha256" -> Hash[Import[$InputFileName, "Text"], "SHA256", "HexString"],
  "matrix_sha256" -> Association@Table["sol" <> ToString[sol] ->
      Hash[Import[FileNameJoin[{lab, "cyl_strata_w" <> ToString[w] <> "_sol" <> ToString[sol] <> ".m"}], "Text"], "SHA256", "HexString"], {sol, 1, 3}],
  "status" -> "LAB HYPOTHESIS -- certify in the open stack before claiming"|>;
Export[FileNameJoin[{lab, "cyl_strata_w" <> ToString[w] <> ".json"}], out, "JSON"];
