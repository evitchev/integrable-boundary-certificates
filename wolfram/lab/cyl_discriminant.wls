(* Discriminant factorization at weights 10 and 12 (preregistered:
   results/cyl_discriminant_preregistration.md).

   For each weight w and Solution sol: load the pinned ad_I3 matrix
   M(t) (wolfram/lab/cyl_strata_w<w>_sol<sol>.m), find the generic rank
   r0 by rational probes, choose a full-rank r0 x r0 submatrix at a
   generic point; its determinant D(t) is a maximal minor, so EVERY
   rank-drop point over Q-bar is a root of D.  Every irreducible factor
   of D is classified by exact rank: rational roots directly, irrational
   factors over Q[t]/(f) with a GF(p) unit-pivot screen (rank_p == r0
   proves no drop).  Each factor also gets its beta^2-image polynomial,
   x = (t-1)/(t+1).  Emits cyl_discriminant_w1012_wls.json next to the
   matrices.  Rank helpers copied verbatim from cyl_strata.wls (the
   certified-w14 lab method).  Final marker: DISCRIMINANT-WLS DONE. *)

lab = DirectoryName[$InputFileName];
args = Rest[$ScriptCommandLine];
weights = If[Length[args] > 0, ToExpression /@ args, {10, 12}];
out = <||>;
probes = {7/3, 11/5, -13/7, 29/11};

curveN[t_] := (t^2 - 25)/(t^2 - 1);
curveS[t_] := -24 t/(t^2 - 1);

rankAt[M_, t0_] := MatrixRank[M /. t -> t0];

rankModFp[M_, f_, p_] := Module[{A, nr, nc, r = 0, c, i, piv, g, ab, inv, k, unitQ},
  If[Mod[Coefficient[f, t, Exponent[f, t]], p] == 0, Return[$Failed]];
  A = Map[PolynomialRemainder[#, f, t, Modulus -> p] &, M, {2}];
  unitQ[e_] := e =!= 0 && Exponent[PolynomialGCD[e, f, Modulus -> p], t] == 0;
  {nr, nc} = Dimensions[A];
  For[c = 1, c <= nc && r < nr, c++,
    piv = 0;
    For[i = r + 1, i <= nr, i++, If[unitQ[A[[i, c]]], piv = i; Break[]]];
    If[piv == 0,
      If[AnyTrue[A[[r + 1 ;;, c]], # =!= 0 &], Return[$Failed]];
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

beta2Poly[f_] := Module[{g = Numerator[Together[f /. t -> (1 + x)/(1 - x)]]},
  ToString[Expand[g/PolynomialGCD @@ CoefficientList[g, x]] // Factor, InputForm]];

Do[
 Do[
  file = FileNameJoin[{lab, "cyl_strata_w" <> ToString[w] <> "_sol" <> ToString[sol] <> ".m"}];
  M = Get[file];
  {nr, nc} = Dimensions[M];
  ranks = MatrixRank[M /. t -> #] & /@ probes;
  r0 = Max[ranks];
  If[Min[ranks] =!= r0, Print["WARNING: probe ranks disagree ", ranks]];
  nullity0 = nc - r0;
  pivotCols[A_] := Module[{rr = RowReduce[A]},
     Function[row, First[Select[Range[Length[row]], row[[#]] =!= 0 &]]] /@ Select[rr, AnyTrue[#, Function[xx, xx =!= 0]] &]];
  Mg = M /. t -> probes[[1]];
  colPiv = pivotCols[Mg];
  rowPiv = pivotCols[Transpose[Mg]];
  sub = M[[rowPiv, colPiv]];
  If[MatrixRank[sub /. t -> probes[[1]]] =!= r0, Print["submatrix selection failed"]; Quit[1]];
  tD = AbsoluteTiming[detD = TimeConstrained[Det[sub], 3600, $Failed];][[1]];
  If[detD === $Failed,
    degBound = Total[Function[row, Max[Exponent[#, t] & /@ row]] /@ sub];
    pts = Select[Range[2, degBound + 12], # =!= 1 && # =!= -1 &];
    nodes = Take[pts, degBound + 1]; extra = Take[Drop[pts, degBound + 1], 3];
    tD = AbsoluteTiming[
      vals = Det[sub /. t -> #] & /@ nodes;
      detD = Expand[InterpolatingPolynomial[Transpose[{nodes, vals}], t]];][[1]];
    If[AnyTrue[extra, (detD /. t -> #) =!= Det[sub /. t -> #] &],
      Print["weight ", w, " Solution ", sol, ": interpolation verification FAILED"]; Quit[1]];
    Print["   (Det by interpolation: degree bound ", degBound, ")"]];
  detD = Expand[detD];
  facs = FactorList[detD];
  jumps = {}; nonjump = {}; poles = {}; unclassified = {};
  Do[
    f = fac[[1]];
    If[FreeQ[f, t], Continue[]];
    deg = Exponent[f, t];
    If[deg == 1 && ((t /. First[Solve[f == 0, t]]) === 1 || (t /. First[Solve[f == 0, t]]) === -1),
      AppendTo[poles, ToString[f, InputForm]]; Continue[]];
    rk = If[deg == 1, rankAt[M, t /. First[Solve[f == 0, t]]],
      Module[{rp = rankModFp[M, f, 2147483647]},
        If[rp === r0, r0, rankModF[M, f]]]];
    If[! IntegerQ[rk], AppendTo[unclassified, ToString[f, InputForm]]; Continue[]];
    rec = <|"factor" -> ToString[f, InputForm], "degree" -> deg,
            "beta2_poly" -> beta2Poly[f], "rank" -> rk, "nullity" -> nc - rk|>;
    If[deg == 1,
      rec["t0"] = ToString[t /. First[Solve[f == 0, t]], InputForm];
      rec["x0"] = ToString[Together[((t /. First[Solve[f == 0, t]]) - 1)/((t /. First[Solve[f == 0, t]]) + 1)], InputForm]];
    If[rk < r0, AppendTo[jumps, rec], AppendTo[nonjump, rec]],
    {fac, facs}];
  Print["w", w, " sol", sol, ": ", nr, "x", nc, " rank ", r0, " detdeg ", Exponent[detD, t],
        " Det ", Round[tD, 0.01], "s jumps ", Length[jumps], " nonjump ", Length[nonjump],
        " poles ", Length[poles], " unclassified ", Length[unclassified]];
  out["w" <> ToString[w] <> "_sol" <> ToString[sol]] = <|"rows" -> nr, "cols" -> nc,
     "generic_rank" -> r0, "generic_multiplicity" -> nullity0,
     "det_degree" -> Exponent[detD, t], "det_seconds" -> Round[tD, 0.01],
     "jumps" -> jumps, "nonjump_factors" -> nonjump, "poles" -> poles,
     "unclassified" -> unclassified|>,
  {sol, 1, 3}],
 {w, weights}];

out["provenance"] = <|"engine" -> $Version,
  "script_sha256" -> Hash[Import[$InputFileName, "Text"], "SHA256", "HexString"],
  "matrix_sha256" -> Association @ Flatten @ Table[
     "w" <> ToString[w] <> "_sol" <> ToString[sol] ->
       Hash[Import[FileNameJoin[{lab, "cyl_strata_w" <> ToString[w] <> "_sol" <> ToString[sol] <> ".m"}], "Text"], "SHA256", "HexString"],
     {w, weights}, {sol, 1, 3}]|>;
Export[FileNameJoin[{lab, "cyl_discriminant_w1012_wls.json"}], out, "JSON"];
Print["DISCRIMINANT-WLS DONE"];
