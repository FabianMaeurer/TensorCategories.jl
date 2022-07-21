(* ::Package:: *)

(* ::Input::Initialization:: *)
obs=Range[0,1];unit=0;

V[a_,b_][c_]:=If[c==Mod[a+b,2],1,0];

objNames=<|0->0,1->1|>;

FData:={X0[a_,b_,c_][d_][1,e_,1][1,f_,1]:>1};
