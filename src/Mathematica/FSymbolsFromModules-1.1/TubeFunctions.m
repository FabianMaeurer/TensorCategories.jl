(* ::Package:: *)

mult[tub[mp_,np_][pp_,qp_][\[Alpha]0p_,xp_,\[Beta]0p_],tub[m_,n_][p_,q_][\[Alpha]0_,x_,\[Beta]0_]]:=mult[tub[mp,np][pp,qp][\[Alpha]0p,xp,\[Beta]0p],tub[m,n][p,q][\[Alpha]0,x,\[Beta]0]]=If[(p===mp&&q===np),Sum[Sqrt[(dim[x]dim[xp])/dim[y]]Conjugate[L[xp,x,n][qp][\[Zeta]0,y,\[Tau]0][\[Beta]0,q,\[Beta]0p]]L[xp,x,m][pp][\[Zeta]0,y,\[Sigma]0][\[Alpha]0,p,\[Alpha]0p]tub[m,n][pp,qp][\[Sigma]0,y,\[Tau]0]
,{y,Select[fusionproduct[xp,x],W[#,n][qp]>0&&W[#,m][pp]>0&]},{\[Sigma]0,W[y,m][pp]},{\[Tau]0,W[y,n][qp]},{\[Zeta]0,V[xp,x][y]}],0]

mult[\[Alpha]0_ A_,B_]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0 mult[A,B]
mult[A_,\[Alpha]0_ B_]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0 mult[A,B]
mult[A_+B_,C_]:=mult[A,C]+mult[B,C]
mult[A_,B_+C_]:=mult[A,B]+mult[A,C]

mult[\[Alpha]0_,B_]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0 B
mult[A_,\[Alpha]0_]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0 A


Nmult[tub[mp_,np_][pp_,qp_][\[Alpha]0p_,xp_,\[Beta]0p_],tub[m_,n_][p_,q_][\[Alpha]0_,x_,\[Beta]0_]]:=Nmult[tub[mp,np][pp,qp][\[Alpha]0p,xp,\[Beta]0p],tub[m,n][p,q][\[Alpha]0,x,\[Beta]0]]=If[(p===mp&&q===np),Sum[Sqrt[(dim[x]dim[xp])/dim[y]]Conjugate[NL[xp,x,n][qp][\[Zeta]0,y,\[Tau]0][\[Beta]0,q,\[Beta]0p]]NL[xp,x,m][pp][\[Zeta]0,y,\[Sigma]0][\[Alpha]0,p,\[Alpha]0p]tub[m,n][pp,qp][\[Sigma]0,y,\[Tau]0]
,{y,Select[fusionproduct[xp,x],W[#,n][qp]>0&&W[#,m][pp]>0&]},{\[Sigma]0,W[y,m][pp]},{\[Tau]0,W[y,n][qp]},{\[Zeta]0,V[xp,x][y]}],0]

Nmult[\[Alpha]0_ A_,B_]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0 Nmult[A,B]
Nmult[A_,\[Alpha]0_ B_]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0 Nmult[A,B]
Nmult[A_+B_,C_]:=Nmult[A,C]+Nmult[B,C]
Nmult[A_,B_+C_]:=Nmult[A,B]+Nmult[A,C]

Nmult[\[Alpha]0_,B_]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0 B
Nmult[A_,\[Alpha]0_]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0 A

idtube:=Sum[tube[m,n][m,n][1,unit,1],{m,obsM},{n,obsM}];


fusionMultiplicity[r1_,r2_][r3_]:=fusionMultiplicity[r1,r2][r3]=Module[{basis=tensorRepBasis[r1,r2],d,T=Variables[tensorRepBasis[r1,r2]],M={},vi,x},
d=Length[basis];
While[
vi=basis . RandomComplex[{-(1+I)/Sqrt[2],(1+I)/Sqrt[2]},d]//Simplify;(*Pick a random vector in r1\[CircleTimes]0r2*)
x=(Chop@*Coefficient)[Nmult[e[r3][0,0],vi],T];(*Project onto rep r3, compute the associated column vector in the picture basis*)

If[M=={},M={x};True,MatrixRank[M]<MatrixRank[AppendTo[M,x]]](*If the new vector is linearly indep, continue*)
];
Return[MatrixRank[M]]
]


VEnd[a_,b_][c_]:=VEnd[a,b][c]=fusionMultiplicity[a,b][c]

fusionproductEnd[a_,b_]/;MemberQ[reps,a]&&MemberQ[reps,b]:=fusionproductEnd[a,b]=DeleteCases[Table[If[VEnd[a,b][c]=!=0,c],{c,reps}],Null]

unitEnd:=unitEnd=Do[If[DeleteDuplicates[Flatten[Table[fusionproductEnd[a,u]=={a}==fusionproductEnd[u,a],{a,reps}]]]=={True},Return[u]],{u,reps}]

dimEnd[a_]/;MemberQ[reps,a]:=Module[{d0,dimeqs,dimsol,da},
dimeqs=Join[
Outer[d0[#1]d0[#2]==Sum[VEnd[#1,#2][x]d0[x],{x,fusionproductEnd[#1,#2]}]&,reps,reps]//Flatten
,d0[#]>=1&/@reps];
da:=d0[a]//.Solve[dimeqs][[1]]//RootReduce//FullSimplify//ToRadicals;
Return[da]];

dtotEnd[]:=Sqrt[Total[dimEnd[#]^2&/@reps]]//RootReduce//FullSimplify//ToRadicals

dualEnd[a_]/;MemberQ[reps,a]:=Module[{v=Abs[VEnd[a,#][unitEnd]]&/@reps},reps[[Position[v,1][[1,1]]]]]
OverBar[a_]/;MemberQ[reps,a]:=dualEnd[a]


tube[m_,n_][p_,q_][\[Alpha]0_,x_,\[Beta]0_]/;(\[Alpha]0<=W[x,m][p]&&\[Beta]0<=W[x,n][q]):=tub[m,n][p,q][\[Alpha]0,x,\[Beta]0];
tube[m_,n_][p_,q_][\[Alpha]0_,x_,\[Beta]0_]:=0;

tubes[m_,n_][p_,q_]:=DeleteCases[Table[tube[m,n][p,q][\[Alpha]0,x,\[Beta]0],{x,obs},{\[Alpha]0,W[x,m][p]},{\[Beta]0,W[x,n][q]}],0]

pictureBasis[]:=pictureBasis[]=DeleteDuplicates[Flatten[Table[tubes[m,n][p,q],{m,obsM},{n,obsM},{p,obsM},{q,obsM}]]];


star[tub[m_,n_][p_,q_][\[Alpha]0_,x_,\[Beta]0_]]:=star[tub[m,n][p,q][\[Alpha]0,x,\[Beta]0]]=Sum[dim[x]Sqrt[(dimM[m]dimM[n])/(dimM[p]dimM[q])]Conjugate[L[dual[x],x,m][m][1,unit,1][\[Alpha]0,p,\[Sigma]0]]L[dual[x],x,n][n][1,unit,1][\[Beta]0,q,\[Tau]0]tub[p,q][m,n][\[Sigma]0,dual[x],\[Tau]0],{\[Sigma]0,W[dual[x],p][m]},{\[Tau]0,W[dual[x],q][n]}]//RootReduce
star[\[Alpha]0_ A_]/;FreeQ[\[Alpha]0,tub]:=Conjugate[\[Alpha]0]star[A]
star[A_+B_]:=star[A]+star[B]
star[\[Alpha]0_]/;FreeQ[\[Alpha]0,tub]:=Conjugate[\[Alpha]0]


\[Omega]0[tub[m_,n_][p_,q_][\[Alpha]0_,x_,\[Beta]0_]]:=If[x===unit,dimM[m]dimM[n],0]

\[Omega]0[\[Alpha]0_ A_]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0 \[Omega]0[A]
\[Omega]0[A_+B_]:=\[Omega]0[A]+\[Omega]0[B]
\[Omega]0[\[Alpha]0_ ]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0

dot[tub[mp_,np_][pp_,qp_][\[Alpha]0p_,xp_,\[Beta]0p_],tub[m_,n_][p_,q_][\[Alpha]0_,x_,\[Beta]0_]]:=dot[tub[mp,np][pp,qp][\[Alpha]0p,xp,\[Beta]0p],tub[m,n][p,q][\[Alpha]0,x,\[Beta]0]]=\[Omega]0[mult[star[tub[mp,np][pp,qp][\[Alpha]0p,xp,\[Beta]0p]],tub[m,n][p,q][\[Alpha]0,x,\[Beta]0]]]//RootReduce
dot[\[Alpha]0_ A_,B_]/;FreeQ[\[Alpha]0,tub]:=Conjugate[\[Alpha]0]dot[A,B]
dot[A_,\[Alpha]0_ B_]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0 dot[A,B]

dot[\[Alpha]0_,\[Beta]0_]/;FreeQ[\[Alpha]0,tub]&&FreeQ[\[Beta]0,tub]:=Conjugate[\[Alpha]0]\[Beta]0

dot[A_+B_,C_]:=dot[A,C]+dot[B,C]
dot[A_,B_+C_]:=dot[A,B]+dot[A,C]

Ndot[tub[mp_,np_][pp_,qp_][\[Alpha]0p_,xp_,\[Beta]0p_],tub[m_,n_][p_,q_][\[Alpha]0_,x_,\[Beta]0_]]:=dot[tub[mp,np][pp,qp][\[Alpha]0p,xp,\[Beta]0p],tub[m,n][p,q][\[Alpha]0,x,\[Beta]0]]=\[Omega]0[Nmult[star[tub[mp,np][pp,qp][\[Alpha]0p,xp,\[Beta]0p]],tub[m,n][p,q][\[Alpha]0,x,\[Beta]0]]]//RootReduce
Ndot[\[Alpha]0_ A_,B_]/;FreeQ[\[Alpha]0,tub]:=Conjugate[\[Alpha]0]Ndot[A,B]
Ndot[A_,\[Alpha]0_ B_]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0 Ndot[A,B]

Ndot[A_+B_,C_]:=Ndot[A,C]+Ndot[B,C]
Ndot[A_,B_+C_]:=Ndot[A,B]+Ndot[A,C]

norm[a_]:=Sqrt[dot[a,a]]


compose[tub[m_,n_][p_,q_][\[Alpha]0_,x_,\[Beta]0_],tub[mp_,np_][pp_,qp_][\[Alpha]0p_,xp_,\[Beta]0p_]]:=If[q===pp,TensorProduct[tub[m,n][p,q][\[Alpha]0,x,\[Beta]0],tub[mp,np][pp,qp][\[Alpha]0p,xp,\[Beta]0p]],0]

compose[\[Alpha]0_ A_,B_]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0 compose[A,B]
compose[A_,\[Alpha]0_ B_]/;FreeQ[\[Alpha]0,tub]:=\[Alpha]0 compose[A,B]
compose[A_+B_,C_]:=compose[A,C]+compose[B,C]
compose[A_,B_+C_]:=compose[A,B]+compose[A,C]


mult[tub[mpp_,npp_][ppp_,qpp_][\[Alpha]0pp_,xpp_,\[Beta]0pp_],TensorProduct[tub[m_,n_][p_,q_][\[Alpha]0_,x_,\[Beta]0_],tub[mp_,np_][pp_,qp_][\[Alpha]0p_,xp_,\[Beta]0p_]]]:=
mult[tub[mpp,npp][ppp,qpp][\[Alpha]0pp,xpp,\[Beta]0pp],TensorProduct[tub[m,n][p,q][\[Alpha]0,x,\[Beta]0],tub[mp,np][pp,qp][\[Alpha]0p,xp,\[Beta]0p]]]=(If[p===mpp&&qp===npp,Sum[
Sqrt[(dimM[r]dim[x]dim[xpp]dim[xp]dim[xpp])/(dim[xpp]dimM[q]dim[y]dim[yp])]
Conjugate[L[xpp,x,n][r][\[Zeta]1,y,\[Tau]0][\[Beta]0,q,\[Zeta]0]]
L[xpp,x,m][ppp][\[Zeta]1,y,\[Sigma]0][\[Alpha]0,p,\[Alpha]0pp]
Conjugate[L[xpp,xp,np][qpp][\[Zeta]2,yp,\[Tau]0p][\[Beta]0p,qp,\[Beta]0pp]]
L[xpp,xp,mp][r][\[Zeta]2,yp,\[Sigma]0p][\[Alpha]0p,q,\[Zeta]0]
TensorProduct[tub[m,n][ppp,r][\[Sigma]0,y,\[Tau]0],tub[mp,np][r,qpp][\[Sigma]0p,yp,\[Tau]0p]],
{r,fusionproductM[xpp,q]},{y,Select[fusionproduct[xpp,x],W[#,n][r]>0&&W[#,m][ppp]>0&]},{yp,Select[fusionproduct[xpp,xp],W[#,np][qpp]>0&&W[#,mp][r]>0&]},
{\[Zeta]0,W[xpp,q][r]},{\[Zeta]1,V[xpp,x][y]},{\[Zeta]2,V[xpp,xp][yp]},{\[Sigma]0,W[y,m][ppp]},{\[Sigma]0p,W[yp,mp][r]},{\[Tau]0,W[y,n][r]},{\[Tau]0p,W[yp,np][qpp]}]
,0]
)

Nmult[tub[mpp_,npp_][ppp_,qpp_][\[Alpha]0pp_,xpp_,\[Beta]0pp_],TensorProduct[tub[m_,n_][p_,q_][\[Alpha]0_,x_,\[Beta]0_],tub[mp_,np_][pp_,qp_][\[Alpha]0p_,xp_,\[Beta]0p_]]]:=
Nmult[tub[mpp,npp][ppp,qpp][\[Alpha]0pp,xpp,\[Beta]0pp],TensorProduct[tub[m,n][p,q][\[Alpha]0,x,\[Beta]0],tub[mp,np][pp,qp][\[Alpha]0p,xp,\[Beta]0p]]]=(If[p===mpp&&qp===npp,Sum[
Sqrt[(dimM[r]dim[x]dim[xpp]dim[xp]dim[xpp])/(dim[xpp]dimM[q]dim[y]dim[yp])]
Conjugate[NL[xpp,x,n][r][\[Zeta]1,y,\[Tau]0][\[Beta]0,q,\[Zeta]0]]
NL[xpp,x,m][ppp][\[Zeta]1,y,\[Sigma]0][\[Alpha]0,p,\[Alpha]0pp]
Conjugate[NL[xpp,xp,np][qpp][\[Zeta]2,yp,\[Tau]0p][\[Beta]0p,qp,\[Beta]0pp]]
NL[xpp,xp,mp][r][\[Zeta]2,yp,\[Sigma]0p][\[Alpha]0p,q,\[Zeta]0]
TensorProduct[tub[m,n][ppp,r][\[Sigma]0,y,\[Tau]0],tub[mp,np][r,qpp][\[Sigma]0p,yp,\[Tau]0p]],
{r,fusionproductM[xpp,q]},{y,Select[fusionproduct[xpp,x],W[#,n][r]>0&&W[#,m][ppp]>0&]},{yp,Select[fusionproduct[xpp,xp],W[#,np][qpp]>0&&W[#,mp][r]>0&]},
{\[Zeta]0,W[xpp,q][r]},{\[Zeta]1,V[xpp,x][y]},{\[Zeta]2,V[xpp,xp][yp]},{\[Sigma]0,W[y,m][ppp]},{\[Sigma]0p,W[yp,mp][r]},{\[Tau]0,W[y,n][r]},{\[Tau]0p,W[yp,np][qpp]}]
,0]//Simplify
)


dot[TensorProduct[tub[r_,s_][t_,u_][\[Sigma]0_,y_,\[Tau]0_],tub[rp_,sp_][tp_,up_][\[Sigma]0p_,yp_,\[Tau]0p_]],TensorProduct[tub[m_,n_][p_,q_][\[Alpha]0_,x_,\[Beta]0_],tub[mp_,np_][pp_,qp_][\[Alpha]0p_,xp_,\[Beta]0p_]]]:=(dot[tub[r,s][t,u][\[Sigma]0,y,\[Tau]0],tub[m,n][p,q][\[Alpha]0,x,\[Beta]0]]dot[tub[rp,sp][tp,up][\[Sigma]0p,yp,\[Tau]0p],tub[mp,np][pp,qp][\[Alpha]0p,xp,\[Beta]0p]])/Sqrt[dimM[q]dimM[u]]
Ndot[TensorProduct[tub[r_,s_][t_,u_][\[Sigma]0_,y_,\[Tau]0_],tub[rp_,sp_][tp_,up_][\[Sigma]0p_,yp_,\[Tau]0p_]],TensorProduct[tub[m_,n_][p_,q_][\[Alpha]0_,x_,\[Beta]0_],tub[mp_,np_][pp_,qp_][\[Alpha]0p_,xp_,\[Beta]0p_]]]:=(Ndot[tub[r,s][t,u][\[Sigma]0,y,\[Tau]0],tub[m,n][p,q][\[Alpha]0,x,\[Beta]0]]Ndot[tub[rp,sp][tp,up][\[Sigma]0p,yp,\[Tau]0p],tub[mp,np][pp,qp][\[Alpha]0p,xp,\[Beta]0p]])/Sqrt[dimM[q]dimM[u]]


tubVars[X_]:=Cases[Variables[X],_?(!FreeQ[#,tub]&)];
toVec[X_]:=Coefficient[X,tubVars[X]]//RootReduce;

toNum[X_]:=Module[{basis=tubVars[X]},N[toVec[X]] . basis]

clean[X_]:=Block[{Y=X},
Y=Collect[Y,tubVars[Y],RootReduce@*Simplify@*RootReduce];
Return[Y]
];

multclean[i_,j_]:=multclean[i,j]=(clean@*mult)[i,j]


e[r_][i_,j_]/;(j>i):=e[r][i,j]=star[e[r][j,i]]//RootReduce//Simplify
e[r_][i_,j_]/;(0<=j<=i<dimVS[r]):=e[r][i,j]=multclean[e[r][i,0],e[r][0,j]]


v[r_][i_]/;MemberQ[reps,r]&&i<dimVS[r]:=v[r][i]=e[r][i,0]//clean
v[x_,y_][i_,j_]:=v[x,y][i,j]=(clean@*compose)[v[x][i],v[y][j]]

repBasis[r_]/;MemberQ[reps,r]:=(v[r][#]&/@Range[0,dimVS[r]-1])
tensorRepBasis[x_,y_]:=tensorRepBasis[x,y]=DeleteCases[DeleteDuplicates[Flatten[Table[(clean@*compose)[v0,v1],{v0,repBasis[x]},{v1,repBasis[y]}]]],0](*Note, this basis may be overcomplete. That won't matter for our purposes*)


isometricQ:=isometricQ=DeleteDuplicates[Flatten[Table[
RootReduce[KroneckerDelta[x,y](dot)[v[\[Gamma]0][i],v[\[Gamma]0][j]]
==
RootReduce[(dot)[VEmbedding[\[Gamma]0->\[Alpha]0\[CircleTimes]\[Beta]0,x][[;;,i+1]] . tensorRepBasis[\[Alpha]0,\[Beta]0],VEmbedding[\[Gamma]0->\[Alpha]0\[CircleTimes]\[Beta]0,y][[;;,j+1]] . tensorRepBasis[\[Alpha]0,\[Beta]0]]]]
,
{\[Alpha]0,reps},{\[Beta]0,reps},{\[Gamma]0,fusionproductEnd[\[Alpha]0,\[Beta]0]},{i,0,dimVS[\[Gamma]0]-1},{j,0,dimVS[\[Gamma]0]-1},{x,VEnd[\[Alpha]0,\[Beta]0][\[Gamma]0]},{y,VEnd[\[Alpha]0,\[Beta]0][\[Gamma]0]}
]]]

NisometricQ:=NisometricQ=DeleteDuplicates[Flatten[Table[
Chop[KroneckerDelta[x,y](Ndot)[v[\[Gamma]0][i],v[\[Gamma]0][j]]
-
(Ndot)[VEmbedding[\[Gamma]0->\[Alpha]0\[CircleTimes]\[Beta]0,x][[;;,i+1]] . tensorRepBasis[\[Alpha]0,\[Beta]0],VEmbedding[\[Gamma]0->\[Alpha]0\[CircleTimes]\[Beta]0,y][[;;,j+1]] . tensorRepBasis[\[Alpha]0,\[Beta]0]]]==0
,
{\[Alpha]0,reps},{\[Beta]0,reps},{\[Gamma]0,fusionproductEnd[\[Alpha]0,\[Beta]0]},{i,0,dimVS[\[Gamma]0]-1},{j,0,dimVS[\[Gamma]0]-1},{x,VEnd[\[Alpha]0,\[Beta]0][\[Gamma]0]},{y,VEnd[\[Alpha]0,\[Beta]0][\[Gamma]0]}
]]]


(*Since we can have A\[CircleTimes]0B==0 for nonzero A and B, we need to embed the basis of nonzero tensor products into the space of all Subscript[v, i]\[CircleTimes]0Subscript[v, j] to produce the 3-tensors called Subsuperscript[V, \[Alpha]0,\[Beta]0, \[Gamma]0;x] in the main text*)
embed[\[Alpha]0_,\[Beta]0_]:=Module[{a,b,r=1,c=1,x,M=ConstantArray[0,{dimVS[\[Alpha]0]dimVS[\[Beta]0],Length[tensorRepBasis[\[Alpha]0,\[Beta]0]]}]},
Do[
If[(clean@*compose)[v[\[Alpha]0][a],v[\[Beta]0][b]]===0,
r++
,
M[[r,c]]=1;
r++;c++
];
,{a,0,dimVS[\[Alpha]0]-1},{b,0,dimVS[\[Beta]0]-1}];
Return[ArrayReshape[M,{dimVS[\[Alpha]0],dimVS[\[Beta]0],Length[tensorRepBasis[\[Alpha]0,\[Beta]0]]}]]
]

VEmbeddingFull[\[Gamma]0_->\[Alpha]0_\[CircleTimes]\[Beta]0_,x_]/;(MemberQ[fusionproductEnd[\[Alpha]0,\[Beta]0],\[Gamma]0]&&x<=VEnd[\[Alpha]0,\[Beta]0][\[Gamma]0]):=VEmbeddingFull[\[Gamma]0->\[Alpha]0\[CircleTimes]\[Beta]0,x]=embed[\[Alpha]0,\[Beta]0] . VEmbedding[\[Gamma]0->\[Alpha]0\[CircleTimes]\[Beta]0,x]

leftBasisEnd[a_,b_,c_][d_]:=Flatten[Table[{\[Alpha]0,e,\[Beta]0},{e,Intersection[fusionproductEnd[a,b],fusionproductEnd[d,dualEnd[c]]]},{\[Alpha]0,VEnd[a,b][e]},{\[Beta]0,VEnd[e,c][d]}],2]
rightBasisEnd[a_,b_,c_][d_]:=Flatten[Table[{\[Sigma]0,f,\[Tau]0},{f,Intersection[fusionproductEnd[b,c],fusionproductEnd[dualEnd[a],d]]},{\[Sigma]0,VEnd[b,c][f]},{\[Tau]0,VEnd[a,f][d]}],2]


leftTree[\[Alpha]0_,\[Beta]0_,\[Gamma]0_][\[Delta]0_][i_,\[Sigma]0_,j_]:=VEmbeddingFull[\[Sigma]0->\[Alpha]0\[CircleTimes]\[Beta]0,i] . VEmbeddingFull[\[Delta]0->\[Sigma]0\[CircleTimes]\[Gamma]0,j]
rightTree[\[Alpha]0_,\[Beta]0_,\[Gamma]0_][\[Delta]0_][k_,\[Tau]0_,l_]:=TensorTranspose[VEmbeddingFull[\[Tau]0->\[Beta]0\[CircleTimes]\[Gamma]0,k] . TensorTranspose[VEmbeddingFull[\[Delta]0->\[Alpha]0\[CircleTimes]\[Tau]0,l],{2,1,3}],{2,3,1,4}]


solveForF[\[Alpha]0_,\[Beta]0_,\[Gamma]0_][\[Delta]0_]:=solveForF[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0]=Module[{A,B,sol,i,j,lB=leftBasisEnd[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0],rB=rightBasisEnd[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0],l,r},
A=Table[RootReduce[Flatten[leftTree[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0][l[[1]],l[[2]],l[[3]]]]],{l,lB}];
B=Table[RootReduce[Flatten[rightTree[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0][r[[1]],r[[2]],r[[3]]]]],{r,rB}];
sol=Transpose[LinearSolve[Transpose[B],Transpose[A]]]//RootReduce;
Flatten[Table[X0End[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0][lB[[i,1]],lB[[i,2]],lB[[i,3]]][rB[[j,1]],rB[[j,2]],rB[[j,3]]]->sol[[i,j]],{i,Length[lB]},{j,Length[rB]}]]
]
solveForF[]:=solveForF[]=Table[solveForF[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0],{\[Alpha]0,reps},{\[Beta]0,reps},{\[Gamma]0,reps},{\[Delta]0,Intersection[Flatten[(fusionproductEnd[#,\[Gamma]0]&/@fusionproductEnd[\[Alpha]0,\[Beta]0])],Flatten[(fusionproductEnd[\[Alpha]0,#]&/@fusionproductEnd[\[Beta]0,\[Gamma]0])]]}]//Flatten

solveForFN[\[Alpha]0_,\[Beta]0_,\[Gamma]0_][\[Delta]0_]:=solveForFN[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0]=Module[{A,B,sol,i,j,lB=leftBasisEnd[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0],rB=rightBasisEnd[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0],l,r},
A=Table[Flatten[N[leftTree[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0][l[[1]],l[[2]],l[[3]]],20]],{l,lB}];
B=Table[Flatten[N[rightTree[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0][r[[1]],r[[2]],r[[3]]],20]],{r,rB}];
sol=Transpose[LinearSolve[Transpose[B],Transpose[A]]]//RootReduce;
Flatten[Table[NX0End[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0][lB[[i,1]],lB[[i,2]],lB[[i,3]]][rB[[j,1]],rB[[j,2]],rB[[j,3]]]->sol[[i,j]],{i,Length[lB]},{j,Length[rB]}]]
]
solveForFN[]:=solveForFN[]=Table[solveForFN[\[Alpha]0,\[Beta]0,\[Gamma]0][\[Delta]0],{\[Alpha]0,reps},{\[Beta]0,reps},{\[Gamma]0,reps},{\[Delta]0,Intersection[Flatten[(fusionproductEnd[#,\[Gamma]0]&/@fusionproductEnd[\[Alpha]0,\[Beta]0])],Flatten[(fusionproductEnd[\[Alpha]0,#]&/@fusionproductEnd[\[Beta]0,\[Gamma]0])]]}]//Flatten


FEnd[a_,b_,c_][d_][\[Alpha]0_,e_,\[Beta]0_][\[Mu]0_,f_,\[Nu]0_]:=0

FEnd[a_,b_,c_][d_][\[Alpha]0_,e_,\[Beta]0_][\[Mu]0_,f_,\[Nu]0_]/;MemberQ[fusionproductEnd[a,b],e]&&MemberQ[fusionproductEnd[b,c],f]&&MemberQ[Intersection[fusionproductEnd[e,c],fusionproductEnd[a,f]],d]&&\[Alpha]0<=VEnd[a,b][e]&&\[Beta]0<=VEnd[e,c][d]&&\[Mu]0<=VEnd[b,c][f]&&\[Nu]0<=VEnd[a,f][d]:=FEnd[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0]=
X0End[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0]//.solveForF[]

NFEnd[a_,b_,c_][d_][\[Alpha]0_,e_,\[Beta]0_][\[Mu]0_,f_,\[Nu]0_]:=0

NFEnd[a_,b_,c_][d_][\[Alpha]0_,e_,\[Beta]0_][\[Mu]0_,f_,\[Nu]0_]/;MemberQ[fusionproductEnd[a,b],e]&&MemberQ[fusionproductEnd[b,c],f]&&MemberQ[Intersection[fusionproductEnd[e,c],fusionproductEnd[a,f]],d]&&\[Alpha]0<=VEnd[a,b][e]&&\[Beta]0<=VEnd[e,c][d]&&\[Mu]0<=VEnd[b,c][f]&&\[Nu]0<=VEnd[a,f][d]:=NFEnd[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0]=
NX0End[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0]//.solveForFN[]


pentagonEndQ:={True}==(Flatten[Table[
RootReduce[
Sum[FEnd[a5,a2,a3][a4][\[Alpha]1,a6,\[Alpha]2][\[Gamma]0,c0,\[Delta]0]FEnd[a0,a1,c0][a4][\[Alpha]0,a5,\[Delta]0][\[Gamma]1,c1,\[Gamma]2],{\[Delta]0,VEnd[a5,c0][a4]}]
-
Sum[FEnd[a0,a1,a2][a6][\[Alpha]0,a5,\[Alpha]1][\[Beta]0,b0,\[Beta]1]FEnd[a0,b0,a3][a4][\[Beta]1,a6,\[Alpha]2][\[Beta]2,c1,\[Gamma]2]FEnd[a1,a2,a3][c1][\[Beta]0,b0,\[Beta]2][\[Gamma]0,c0,\[Gamma]1],{b0,Intersection[fusionproductEnd[a1,a2],fusionproductEnd[dualEnd[a0],a6],fusionproductEnd[c1,dualEnd[a3]]]},{\[Beta]0,VEnd[a1,a2][b0]},{\[Beta]1,VEnd[a0,b0][a6]},{\[Beta]2,VEnd[b0,a3][c1]}]
]==0
,
{a0,reps},{a1,reps},{a2,reps},{a3,reps},{a5,fusionproductEnd[a0,a1]},{a6,fusionproductEnd[a5,a2]},{a4,fusionproductEnd[a6,a3]},{c0,fusionproductEnd[a2,a3]},{c1,Intersection[fusionproductEnd[a1,c0],fusionproductEnd[dualEnd[a0],a4]]},
{\[Alpha]0,VEnd[a0,a1][a5]},{\[Alpha]1,VEnd[a5,a2][a6]},{\[Alpha]2,VEnd[a6,a3][a4]},{\[Gamma]0,VEnd[a2,a3][c0]},{\[Gamma]1,VEnd[a1,c0][c1]},{\[Gamma]2,VEnd[a0,c1][a4]}
]]//DeleteDuplicates)

NpentagonEndQ:={True}==(Flatten[Table[
Chop[
Sum[NFEnd[a5,a2,a3][a4][\[Alpha]1,a6,\[Alpha]2][\[Gamma]0,c0,\[Delta]0]NFEnd[a0,a1,c0][a4][\[Alpha]0,a5,\[Delta]0][\[Gamma]1,c1,\[Gamma]2],{\[Delta]0,VEnd[a5,c0][a4]}]
-
Sum[NFEnd[a0,a1,a2][a6][\[Alpha]0,a5,\[Alpha]1][\[Beta]0,b0,\[Beta]1]NFEnd[a0,b0,a3][a4][\[Beta]1,a6,\[Alpha]2][\[Beta]2,c1,\[Gamma]2]NFEnd[a1,a2,a3][c1][\[Beta]0,b0,\[Beta]2][\[Gamma]0,c0,\[Gamma]1],{b0,Intersection[fusionproductEnd[a1,a2],fusionproductEnd[dualEnd[a0],a6],fusionproductEnd[c1,dualEnd[a3]]]},{\[Beta]0,VEnd[a1,a2][b0]},{\[Beta]1,VEnd[a0,b0][a6]},{\[Beta]2,VEnd[b0,a3][c1]}]
,10^-16]==0
,
{a0,reps},{a1,reps},{a2,reps},{a3,reps},{a5,fusionproductEnd[a0,a1]},{a6,fusionproductEnd[a5,a2]},{a4,fusionproductEnd[a6,a3]},{c0,fusionproductEnd[a2,a3]},{c1,Intersection[fusionproductEnd[a1,c0],fusionproductEnd[dualEnd[a0],a4]]},
{\[Alpha]0,VEnd[a0,a1][a5]},{\[Alpha]1,VEnd[a5,a2][a6]},{\[Alpha]2,VEnd[a6,a3][a4]},{\[Gamma]0,VEnd[a2,a3][c0]},{\[Gamma]1,VEnd[a1,c0][c1]},{\[Gamma]2,VEnd[a0,c1][a4]}
]]//DeleteDuplicates)

FEnds:=Table[FEnd[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0],{a,reps},{b,reps},{c,reps},{e,fusionproductEnd[a,b]},{f,fusionproductEnd[b,c]},{d,Intersection[fusionproductEnd[a,f],fusionproductEnd[e,c]]},{\[Alpha]0,VEnd[a,b][e]},{\[Beta]0,VEnd[e,c][d]},{\[Mu]0,VEnd[b,c][f]},{\[Nu]0,VEnd[a,f][d]}]//Flatten
NFEnds:=Table[NFEnd[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Mu]0,f,\[Nu]0],{a,reps},{b,reps},{c,reps},{e,fusionproductEnd[a,b]},{f,fusionproductEnd[b,c]},{d,Intersection[fusionproductEnd[a,f],fusionproductEnd[e,c]]},{\[Alpha]0,VEnd[a,b][e]},{\[Beta]0,VEnd[e,c][d]},{\[Mu]0,VEnd[b,c][f]},{\[Nu]0,VEnd[a,f][d]}]//Flatten

\[Kappa]End[x_]:=FEnd[x,dualEnd[x],x][x][1,unitEnd,1][1,unitEnd,1]dimEnd[x]//FullSimplify
N\[Kappa]End[x_]:=NFEnd[x,dualEnd[x],x][x][1,unitEnd,1][1,unitEnd,1]dimEnd[x]//FullSimplify


uFEnd1[a_,b_,c_][d_]:=
Table[
Sum[FEnd[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[FEnd[a,b,c][d][\[Alpha]0p,ep,\[Beta]0p][\[Sigma]0,f,\[Tau]0]],{f,Intersection[fusionproductEnd[b,c],fusionproductEnd[dualEnd[a],d]]},{\[Sigma]0,VEnd[b,c][f]},{\[Tau]0,VEnd[a,f][d]}]
-If[e===ep&&\[Alpha]0===\[Alpha]0p&&\[Beta]0===\[Beta]0p,1,0]
,{e,fusionproductEnd[a,b]},{ep,fusionproductEnd[a,b]},{\[Alpha]0,VEnd[a,b][e]},{\[Beta]0,VEnd[e,c][d]},{\[Alpha]0p,VEnd[a,b][ep]},{\[Beta]0p,VEnd[ep,c][d]}]

NuFEnd1[a_,b_,c_][d_]:=
Table[
Sum[NFEnd[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[NFEnd[a,b,c][d][\[Alpha]0p,ep,\[Beta]0p][\[Sigma]0,f,\[Tau]0]],{f,Intersection[fusionproductEnd[b,c],fusionproductEnd[dualEnd[a],d]]},{\[Sigma]0,VEnd[b,c][f]},{\[Tau]0,VEnd[a,f][d]}]
-If[e===ep&&\[Alpha]0===\[Alpha]0p&&\[Beta]0===\[Beta]0p,1,0]
,{e,fusionproductEnd[a,b]},{ep,fusionproductEnd[a,b]},{\[Alpha]0,VEnd[a,b][e]},{\[Beta]0,VEnd[e,c][d]},{\[Alpha]0p,VEnd[a,b][ep]},{\[Beta]0p,VEnd[ep,c][d]}]

uFEnd2[a_,b_,c_][d_]:=
Table[
Sum[FEnd[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[FEnd[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0p,fp,\[Tau]0p]],{e,fusionproductEnd[a,b]},{\[Alpha]0,VEnd[a,b][e]},{\[Beta]0,VEnd[e,c][d]}]
-If[f===fp&&\[Sigma]0===\[Sigma]0p&&\[Tau]0===\[Tau]0p,1,0]
,{f,fusionproductEnd[b,c]},{fp,fusionproductEnd[b,c]},{\[Sigma]0,VEnd[b,c][f]},{\[Tau]0,VEnd[a,f][d]},{\[Sigma]0p,VEnd[b,c][fp]},{\[Tau]0p,VEnd[fp,c][d]}]

NuFEnd2[a_,b_,c_][d_]:=
Table[
Sum[NFEnd[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0,f,\[Tau]0]Conjugate[NFEnd[a,b,c][d][\[Alpha]0,e,\[Beta]0][\[Sigma]0p,fp,\[Tau]0p]],{e,fusionproductEnd[a,b]},{\[Alpha]0,VEnd[a,b][e]},{\[Beta]0,VEnd[e,c][d]}]
-If[f===fp&&\[Sigma]0===\[Sigma]0p&&\[Tau]0===\[Tau]0p,1,0]
,{f,fusionproductEnd[b,c]},{fp,fusionproductEnd[b,c]},{\[Sigma]0,VEnd[b,c][f]},{\[Tau]0,VEnd[a,f][d]},{\[Sigma]0p,VEnd[b,c][fp]},{\[Tau]0p,VEnd[fp,c][d]}]


unitaryEndQ:={True}==Module[{X},
Flatten[Table[
X=Flatten[{uFEnd1[a,b,c][d],uFEnd2[a,b,c][d]}];
RootReduce[X==ConstantArray[0,Dimensions[X]]]
,{a,reps},{b,reps},{c,reps},{d,DeleteDuplicates[Flatten[fusionproductEnd[#,c]&/@fusionproductEnd[a,b]]]}]]//DeleteDuplicates
]

NunitaryEndQ:={True}==Module[{X},
Flatten[Table[
X=Flatten[{Chop[NuFEnd1[a,b,c][d],10^-16],Chop[NuFEnd2[a,b,c][d],10^-16]}];
X==ConstantArray[0,Dimensions[X]]
,{a,reps},{b,reps},{c,reps},{d,DeleteDuplicates[Flatten[fusionproductEnd[#,c]&/@fusionproductEnd[a,b]]]}]]//DeleteDuplicates
]
