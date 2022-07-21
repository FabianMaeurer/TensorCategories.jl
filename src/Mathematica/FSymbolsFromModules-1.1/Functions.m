(* ::Package:: *)

polarForm=Expand[#/.z_?NumericQ:>Abs[z] Exp[I Arg[z]]]&;


(* ::Section:: *)
(*Fusion Ring*)


fusionproduct[a_,b_]/;MemberQ[obs,a]&&MemberQ[obs,b]:=DeleteCases[Table[If[V[a,b][c]=!=0,c],{c,obs}],Null]

If[!ValueQ[unit],unit:=unit=Do[If[DeleteDuplicates[Flatten[Table[fusionproduct[a,u]=={a}==fusionproduct[u,a],{a,obs}]]]=={True},Return[u]],{u,obs}]];

dim[a_]/;MemberQ[obs,a]:=Module[{d0,dimeqs,dimsol,da},
dimeqs=Join[
Outer[d0[#1]d0[#2]==Sum[V[#1,#2][x]d0[x],{x,fusionproduct[#1,#2]}]&,obs,obs]//Flatten
,d0[#]>=1&/@obs];
da:=d0[a]//.Solve[dimeqs][[1]]//RootReduce//FullSimplify//ToRadicals;
Return[da]];

dtot[]:=Sqrt[Total[dim[#]^2&/@obs]]//RootReduce//FullSimplify//ToRadicals
dtot[X_]:=Sqrt[Total[dim[#]^2&/@X]]//RootReduce//FullSimplify//ToRadicals

dual[a_]/;MemberQ[obs,a]:=Module[{v=Abs[V[a,#][unit]]&/@obs},obs[[Position[v,1][[1,1]]]]]
OverBar[a_]/;MemberQ[obs,a]:=dual[a]

dim[{A__}]:=Sum[dim[a],{a,{A}}]


(* ::Section:: *)
(*Fusion category*)


F[a_,b_,c_][d_][\[Alpha]0_,e_,\[Beta]0_][\[Mu]0_,f_,\[Nu]0_]:=0

F[a_,b_,c_][d_][\[Alpha]0_,e_,\[Beta]0_][\[Mu]0_,f_,\[Nu]0_]/;MemberQ[fusionproduct[a,b],e]&&MemberQ[fusionproduct[b,c],f]&&MemberQ[Intersection[fusionproduct[e,c],fusionproduct[a,f]],d]&&\[Alpha]0<=V[a,b][e]&&\[Beta]0<=V[e,c][d]&&\[Mu]0<=V[b,c][f]&&\[Nu]0<=V[a,f][d]:=F[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0]=
X0[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0]//.FData

NF[a_,b_,c_][d_][\[Alpha]0_,e_,\[Beta]0_][\[Sigma]0_,f_,\[Tau]0_]:=NF[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]=N[F[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0],20]


pentagonQ:={True}==(Flatten[Table[
RootReduce[
Sum[F[a5,a2,a3][a4][\[Alpha]1,a6,\[Alpha]2][\[Gamma]0,c0,\[Delta]0]F[a0,a1,c0][a4][\[Alpha]0,a5,\[Delta]0][\[Gamma]1,c1,\[Gamma]2],{\[Delta]0,V[a5,c0][a4]}]
-
Sum[F[a0,a1,a2][a6][\[Alpha]0,a5,\[Alpha]1][\[Beta]0,b0,\[Beta]1]F[a0,b0,a3][a4][\[Beta]1,a6,\[Alpha]2][\[Beta]2,c1,\[Gamma]2]F[a1,a2,a3][c1][\[Beta]0,b0,\[Beta]2][\[Gamma]0,c0,\[Gamma]1],{b0,Intersection[fusionproduct[a1,a2],fusionproduct[dual[a0],a6],fusionproduct[c1,dual[a3]]]},{\[Beta]0,V[a1,a2][b0]},{\[Beta]1,V[a0,b0][a6]},{\[Beta]2,V[b0,a3][c1]}]
]==0
,
{a0,obs},{a1,obs},{a2,obs},{a3,obs},{a5,fusionproduct[a0,a1]},{a6,fusionproduct[a5,a2]},{a4,fusionproduct[a6,a3]},{c0,fusionproduct[a2,a3]},{c1,Intersection[fusionproduct[a1,c0],fusionproduct[dual[a0],a4]]},
{\[Alpha]0,V[a0,a1][a5]},{\[Alpha]1,V[a5,a2][a6]},{\[Alpha]2,V[a6,a3][a4]},{\[Gamma]0,V[a2,a3][c0]},{\[Gamma]1,V[a1,c0][c1]},{\[Gamma]2,V[a0,c1][a4]}
]]//DeleteDuplicates)

NpentagonQ:={True}==(Flatten[Table[
Chop[
Sum[NF[a5,a2,a3][a4][\[Alpha]1,a6,\[Alpha]2][\[Gamma]0,c0,\[Delta]0]NF[a0,a1,c0][a4][\[Alpha]0,a5,\[Delta]0][\[Gamma]1,c1,\[Gamma]2],{\[Delta]0,V[a5,c0][a4]}]
-
Sum[NF[a0,a1,a2][a6][\[Alpha]0,a5,\[Alpha]1][\[Beta]0,b0,\[Beta]1]NF[a0,b0,a3][a4][\[Beta]1,a6,\[Alpha]2][\[Beta]2,c1,\[Gamma]2]NF[a1,a2,a3][c1][\[Beta]0,b0,\[Beta]2][\[Gamma]0,c0,\[Gamma]1],{b0,Intersection[fusionproduct[a1,a2],fusionproduct[dual[a0],a6],fusionproduct[c1,dual[a3]]]},{\[Beta]0,V[a1,a2][b0]},{\[Beta]1,V[a0,b0][a6]},{\[Beta]2,V[b0,a3][c1]}]
,10^-16]==0
,
{a0,obs},{a1,obs},{a2,obs},{a3,obs},{a5,fusionproduct[a0,a1]},{a6,fusionproduct[a5,a2]},{a4,fusionproduct[a6,a3]},{c0,fusionproduct[a2,a3]},{c1,Intersection[fusionproduct[a1,c0],fusionproduct[dual[a0],a4]]},
{\[Alpha]0,V[a0,a1][a5]},{\[Alpha]1,V[a5,a2][a6]},{\[Alpha]2,V[a6,a3][a4]},{\[Gamma]0,V[a2,a3][c0]},{\[Gamma]1,V[a1,c0][c1]},{\[Gamma]2,V[a0,c1][a4]}
]]//DeleteDuplicates)

Fs:=Table[F[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0],{a,obs},{b,obs},{c,obs},{e,fusionproduct[a,b]},{f,fusionproduct[b,c]},{d,Intersection[fusionproduct[a,f],fusionproduct[e,c]]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,V[e,c][d]},{\[Mu]0,V[b,c][f]},{\[Nu]0,V[a,f][d]}]//Flatten
NFs:=Table[NF[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0],{a,obs},{b,obs},{c,obs},{e,fusionproduct[a,b]},{f,fusionproduct[b,c]},{d,Intersection[fusionproduct[a,f],fusionproduct[e,c]]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,V[e,c][d]},{\[Mu]0,V[b,c][f]},{\[Nu]0,V[a,f][d]}]//Flatten

\[Kappa][x_]:=F[x,dual[x],x][x][1,unit,1][1,unit,1]dim[x]//FullSimplify


uF1[a_,b_,c_][d_]:=
Table[
Sum[F[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[F[a,b,c][d][\[Alpha]0p,ep,\[Beta]0p][\[Sigma]0,f,\[Tau]0]],{f,Intersection[fusionproduct[b,c],fusionproduct[dual[a],d]]},{\[Sigma]0,V[b,c][f]},{\[Tau]0,V[a,f][d]}]
-If[e===ep&&\[Alpha]0===\[Alpha]0p&&\[Beta]0===\[Beta]0p,1,0]
,{e,fusionproduct[a,b]},{ep,fusionproduct[a,b]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,V[e,c][d]},{\[Alpha]0p,V[a,b][ep]},{\[Beta]0p,V[ep,c][d]}]

NuF1[a_,b_,c_][d_]:=
Table[
Sum[NF[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[NF[a,b,c][d][\[Alpha]0p,ep,\[Beta]0p][\[Sigma]0,f,\[Tau]0]],{f,Intersection[fusionproduct[b,c],fusionproduct[dual[a],d]]},{\[Sigma]0,V[b,c][f]},{\[Tau]0,V[a,f][d]}]
-If[e===ep&&\[Alpha]0===\[Alpha]0p&&\[Beta]0===\[Beta]0p,1,0]
,{e,fusionproduct[a,b]},{ep,fusionproduct[a,b]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,V[e,c][d]},{\[Alpha]0p,V[a,b][ep]},{\[Beta]0p,V[ep,c][d]}]

uF2[a_,b_,c_][d_]:=
Table[
Sum[F[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[F[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0p,fp,\[Tau]0p]],{e,fusionproduct[a,b]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,V[e,c][d]}]
-If[f===fp&&\[Sigma]0===\[Sigma]0p&&\[Tau]0===\[Tau]0p,1,0]
,{f,fusionproduct[b,c]},{fp,fusionproduct[b,c]},{\[Sigma]0,V[b,c][f]},{\[Tau]0,V[a,f][d]},{\[Sigma]0p,V[b,c][fp]},{\[Tau]0p,V[fp,c][d]}]

NuF2[a_,b_,c_][d_]:=
Table[
Sum[NF[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[NF[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0p,fp,\[Tau]0p]],{e,fusionproduct[a,b]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,V[e,c][d]}]
-If[f===fp&&\[Sigma]0===\[Sigma]0p&&\[Tau]0===\[Tau]0p,1,0]
,{f,fusionproduct[b,c]},{fp,fusionproduct[b,c]},{\[Sigma]0,V[b,c][f]},{\[Tau]0,V[a,f][d]},{\[Sigma]0p,V[b,c][fp]},{\[Tau]0p,V[fp,c][d]}]


unitaryQ:={True}==Module[{X},
Flatten[Table[
X=Flatten[{uF1[a,b,c][d],uF2[a,b,c][d]}];
RootReduce[X==ConstantArray[0,Dimensions[X]]]
,{a,obs},{b,obs},{c,obs},{d,DeleteDuplicates[Flatten[fusionproduct[#,c]&/@fusionproduct[a,b]]]}]]//DeleteDuplicates
]

NunitaryQ:={True}==Module[{X},
Flatten[Table[
X=Flatten[{Chop[NuF1[a,b,c][d],10^-16],Chop[NuF2[a,b,c][d],10^-16]}];
X==ConstantArray[0,Dimensions[X]]
,{a,obs},{b,obs},{c,obs},{d,DeleteDuplicates[Flatten[fusionproduct[#,c]&/@fusionproduct[a,b]]]}]]//DeleteDuplicates
]


leftBasis[a_,b_,c_][d_]:=Flatten[Table[{\[Alpha]0,e,\[Beta]0},{e,Intersection[fusionproduct[a,b],fusionproduct[d,dual[c]]]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,V[e,c][d]}],2]
rightBasis[a_,b_,c_][d_]:=Flatten[Table[{\[Sigma]0,f,\[Tau]0},{f,Intersection[fusionproduct[b,c],fusionproduct[dual[a],d]]},{\[Sigma]0,V[b,c][f]},{\[Tau]0,V[a,f][d]}],2]

Fblk[a_,b_,c_][d_]:=
Table[
F[a,b,c][d][l[[1]],l[[2]],l[[3]]][r[[1]],r[[2]],r[[3]]],
{l,leftBasis[a,b,c][d]},{r,rightBasis[a,b,c][d]}
]

NFblk[a_,b_,c_][d_]:=
Table[
NF[a,b,c][d][l[[1]],l[[2]],l[[3]]][r[[1]],r[[2]],r[[3]]],
{l,leftBasis[a,b,c][d]},{r,rightBasis[a,b,c][d]}
]


(* ::Section:: *)
(*Module category*)


fusionproductM[a_,b_]/;MemberQ[obs,a]&&MemberQ[obsM,b]:=DeleteCases[Table[If[W[a,b][c]=!=0,c],{c,obsM}],Null]

dimM[a_]/;MemberQ[obsM,a]:=dimM[a]=Module[{d0,dimeqs,dimsol,da},
dimeqs=Join[
Outer[dim[#1]d0[#2]==Sum[W[#1,#2][x]d0[x],{x,fusionproductM[#1,#2]}]&,obs,obsM]//Flatten
,d0[#]>0&/@obsM,{Total[d0[#]^2&/@obsM]==dtot[]^2}];
da:=d0[a]//.Solve[dimeqs][[1]]//RootReduce//FullSimplify//ToRadicals;
Return[da]];


L[a_,b_,c_][d_][\[Alpha]0_,e_,\[Beta]0_][\[Mu]0_,f_,\[Nu]0_]/;MemberQ[fusionproduct[a,b],e]&&MemberQ[fusionproductM[b,c],f]&&MemberQ[Intersection[fusionproductM[e,c],fusionproductM[a,f]],d]&&\[Alpha]0<=V[a,b][e]&&\[Beta]0<=W[e,c][d]&&\[Mu]0<=W[b,c][f]&&\[Nu]0<=W[a,f][d]:=L[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0]=L0[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0]//.MData;
L[a_,b_,c_][d_][\[Alpha]0_,e_,\[Beta]0_][\[Mu]0_,f_,\[Nu]0_]:=0;

NL[a_,b_,c_][d_][\[Alpha]0_,e_,\[Beta]0_][\[Sigma]0_,f_,\[Tau]0_]:=NL[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]=N[L[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0],20]


pentagonMQ:={True}==(Flatten[Table[
RootReduce[
Sum[L[a5,a2,a3][a4][\[Alpha]1,a6,\[Alpha]2][\[Gamma]0,c0,\[Delta]0]L[a0,a1,c0][a4][\[Alpha]0,a5,\[Delta]0][\[Gamma]1,c1,\[Gamma]2],{\[Delta]0,W[a5,c0][a4]}]
-
Sum[F[a0,a1,a2][a6][\[Alpha]0,a5,\[Alpha]1][\[Beta]0,b0,\[Beta]1]L[a0,b0,a3][a4][\[Beta]1,a6,\[Alpha]2][\[Beta]2,c1,\[Gamma]2]L[a1,a2,a3][c1][\[Beta]0,b0,\[Beta]2][\[Gamma]0,c0,\[Gamma]1],
{b0,Intersection[fusionproduct[a1,a2],fusionproduct[dual[a0],a6]]},{\[Beta]0,V[a1,a2][b0]},{\[Beta]1,V[a0,b0][a6]},{\[Beta]2,W[b0,a3][c1]}]
]==0
,
{a0,obs},{a1,obs},{a2,obs},{a3,obsM},{a5,fusionproduct[a0,a1]},{a6,fusionproduct[a5,a2]},{a4,fusionproductM[a6,a3]},{c0,fusionproductM[a2,a3]},{c1,Intersection[fusionproductM[a1,c0],fusionproductM[dual[a0],a4]]},
{\[Alpha]0,V[a0,a1][a5]},{\[Alpha]1,V[a5,a2][a6]},{\[Alpha]2,W[a6,a3][a4]},{\[Gamma]0,W[a2,a3][c0]},{\[Gamma]1,W[a1,c0][c1]},{\[Gamma]2,W[a0,c1][a4]}
]]//DeleteDuplicates)

NpentagonMQ:={True}==(Flatten[Table[
Chop[
Sum[NL[a5,a2,a3][a4][\[Alpha]1,a6,\[Alpha]2][\[Gamma]0,c0,\[Delta]0]NL[a0,a1,c0][a4][\[Alpha]0,a5,\[Delta]0][\[Gamma]1,c1,\[Gamma]2],{\[Delta]0,W[a5,c0][a4]}]
-
Sum[NF[a0,a1,a2][a6][\[Alpha]0,a5,\[Alpha]1][\[Beta]0,b0,\[Beta]1]NL[a0,b0,a3][a4][\[Beta]1,a6,\[Alpha]2][\[Beta]2,c1,\[Gamma]2]NL[a1,a2,a3][c1][\[Beta]0,b0,\[Beta]2][\[Gamma]0,c0,\[Gamma]1],
{b0,Intersection[fusionproduct[a1,a2],fusionproduct[dual[a0],a6]]},{\[Beta]0,V[a1,a2][b0]},{\[Beta]1,V[a0,b0][a6]},{\[Beta]2,W[b0,a3][c1]}]
,10^-16]==0
,
{a0,obs},{a1,obs},{a2,obs},{a3,obsM},{a5,fusionproduct[a0,a1]},{a6,fusionproduct[a5,a2]},{a4,fusionproductM[a6,a3]},{c0,fusionproductM[a2,a3]},{c1,Intersection[fusionproductM[a1,c0],fusionproductM[dual[a0],a4]]},
{\[Alpha]0,V[a0,a1][a5]},{\[Alpha]1,V[a5,a2][a6]},{\[Alpha]2,W[a6,a3][a4]},{\[Gamma]0,W[a2,a3][c0]},{\[Gamma]1,W[a1,c0][c1]},{\[Gamma]2,W[a0,c1][a4]}
]]//DeleteDuplicates)

Ls:=Table[L[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0],{a,obs},{b,obs},{c,obsM},{e,fusionproduct[a,b]},{f,fusionproductM[b,c]},{d,Intersection[fusionproductM[a,f],fusionproductM[e,c]]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,W[e,c][d]},{\[Mu]0,W[b,c][f]},{\[Nu]0,W[a,f][d]}]//Flatten
NLs:=Table[NL[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0],{a,obs},{b,obs},{c,obsM},{e,fusionproduct[a,b]},{f,fusionproductM[b,c]},{d,Intersection[fusionproductM[a,f],fusionproductM[e,c]]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,W[e,c][d]},{\[Mu]0,W[b,c][f]},{\[Nu]0,W[a,f][d]}]//Flatten


uL[a_,b_,c_][d_]:=
Table[
Sum[L[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[L[a,b,c][d][\[Alpha]0p,ep,\[Beta]0p][\[Sigma]0,f,\[Tau]0]],{f,Intersection[fusionproductM[b,c],fusionproductM[dual[a],d]]},{\[Sigma]0,W[b,c][f]},{\[Tau]0,W[a,f][d]}]
-If[e===ep&&\[Alpha]0===\[Alpha]0p&&\[Beta]0===\[Beta]0p,1,0]
,{e,fusionproduct[a,b]},{ep,fusionproduct[a,b]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,W[e,c][d]},{\[Alpha]0p,V[a,b][ep]},{\[Beta]0p,W[ep,c][d]}]

NuL[a_,b_,c_][d_]:=
Table[
Sum[NF[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[NF[a,b,c][d][\[Alpha]0p,ep,\[Beta]0p][\[Sigma]0,f,\[Tau]0]],{f,Intersection[fusionproduct[b,c],fusionproduct[dual[a],d]]},{\[Sigma]0,V[b,c][f]},{\[Tau]0,V[a,f][d]}]
-If[e===ep&&\[Alpha]0===\[Alpha]0p&&\[Beta]0===\[Beta]0p,1,0]
,{e,fusionproduct[a,b]},{ep,fusionproduct[a,b]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,V[e,c][d]},{\[Alpha]0p,V[a,b][ep]},{\[Beta]0p,V[ep,c][d]}]


uL1[a_,b_,c_][d_]:=
Table[
Sum[L[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[L[a,b,c][d][\[Alpha]0p,ep,\[Beta]0p][\[Sigma]0,f,\[Tau]0]],{f,Intersection[fusionproductM[b,c],fusionproductM[dual[a],d]]},{\[Sigma]0,W[b,c][f]},{\[Tau]0,W[a,f][d]}]
-If[e===ep&&\[Alpha]0===\[Alpha]0p&&\[Beta]0===\[Beta]0p,1,0]
,{e,fusionproduct[a,b]},{ep,fusionproduct[a,b]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,W[e,c][d]},{\[Alpha]0p,V[a,b][ep]},{\[Beta]0p,W[ep,c][d]}]

NuL1[a_,b_,c_][d_]:=
Table[
Sum[NL[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[NL[a,b,c][d][\[Alpha]0p,ep,\[Beta]0p][\[Sigma]0,f,\[Tau]0]],{f,Intersection[fusionproductM[b,c],fusionproductM[dual[a],d]]},{\[Sigma]0,W[b,c][f]},{\[Tau]0,W[a,f][d]}]
-If[e===ep&&\[Alpha]0===\[Alpha]0p&&\[Beta]0===\[Beta]0p,1,0]
,{e,fusionproduct[a,b]},{ep,fusionproduct[a,b]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,W[e,c][d]},{\[Alpha]0p,V[a,b][ep]},{\[Beta]0p,W[ep,c][d]}]

uL2[a_,b_,c_][d_]:=
Table[
Sum[L[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[L[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0p,fp,\[Tau]0p]],{e,fusionproduct[a,b]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,W[e,c][d]}]
-If[f===fp&&\[Sigma]0===\[Sigma]0p&&\[Tau]0===\[Tau]0p,1,0]
,{f,fusionproductM[b,c]},{fp,fusionproductM[b,c]},{\[Sigma]0,W[b,c][f]},{\[Tau]0,W[a,f][d]},{\[Sigma]0p,W[b,c][fp]},{\[Tau]0p,W[fp,c][d]}]

NuL2[a_,b_,c_][d_]:=
Table[
Sum[NL[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[NL[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0p,fp,\[Tau]0p]],{e,fusionproduct[a,b]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,W[e,c][d]}]
-If[f===fp&&\[Sigma]0===\[Sigma]0p&&\[Tau]0===\[Tau]0p,1,0]
,{f,fusionproductM[b,c]},{fp,fusionproductM[b,c]},{\[Sigma]0,W[b,c][f]},{\[Tau]0,W[a,f][d]},{\[Sigma]0p,W[b,c][fp]},{\[Tau]0p,W[fp,c][d]}]


unitaryMQ:={True}==Module[{X},
Flatten[Table[
X=Flatten[{uL1[a,b,c][d],uL2[a,b,c][d]}];
RootReduce[X==ConstantArray[0,Dimensions[X]]]
,{a,obs},{b,obs},{c,obsM},{d,DeleteDuplicates[Flatten[fusionproductM[#,c]&/@fusionproduct[a,b]]]}]]//DeleteDuplicates
]

NunitaryMQ:={True}==Module[{X},
Flatten[Table[
X=Flatten[{Chop[NuL1[a,b,c][d],10^-16],Chop[NuL2[a,b,c][d],10^-16]}];
X==ConstantArray[0,Dimensions[X]]
,{a,obs},{b,obs},{c,obsM},{d,DeleteDuplicates[Flatten[fusionproductM[#,c]&/@fusionproduct[a,b]]]}]]//DeleteDuplicates
]


leftBasisM[a_,b_,c_][d_]:=Flatten[Table[ba[\[Alpha]0,e,\[Beta]0],{e,Select[fusionproduct[a,b],W[#,c][d]>0&]},{\[Alpha]0,V[a,b][e]},{\[Beta]0,W[e,c][d]}]]
rightBasisM[a_,b_,c_][d_]:=Flatten[Table[ba[\[Sigma]0,f,\[Tau]0],{f,Intersection[fusionproductM[b,c],fusionproductM[dual[a],d]]},{\[Sigma]0,W[b,c][f]},{\[Tau]0,W[a,f][d]}]]

Lblk[a_,b_,c_][d_]:=
Table[
L[a,b,c][d][l[[1]],l[[2]],l[[3]]][r[[1]],r[[2]],r[[3]]],
{l,leftBasisM[a,b,c][d]},{r,rightBasisM[a,b,c][d]}
]

NLblk[a_,b_,c_][d_]:=
Table[
NL[a,b,c][d][l[[1]],l[[2]],l[[3]]][r[[1]],r[[2]],r[[3]]],
{l,leftBasisM[a,b,c][d]},{r,rightBasisM[a,b,c][d]}
]
