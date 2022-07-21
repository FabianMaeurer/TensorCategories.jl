(* ::Package:: *)

(* ::Subsubsection:: *)
(*General data*)


reps={1,\[Mu],\[Nu]};
dimVS[1]:=3;dimVS[\[Mu]]:=3;dimVS[\[Nu]]:=5;(*Vector space dimensions*)

repNames[x_]/;MemberQ[reps,x]:=x


(* ::Subsubsection:: *)
(*Representations*)


(* ::Input::Initialization:: *)
e[1][0,0]:=tub[0,0][0,0][1,0,1]
e[1][1,0]:=tub[0,0][1,1][1,1,1]
e[1][2,0]:=tub[0,0][2,2][1,2,1]/Sqrt[2]

e[\[Mu]][0,0]:=tub[0,1][0,1][1,0,1]

e[\[Mu]][1,0]:=tub[0,1][1,0][1,1,1]

e[\[Mu]][2,0]:=tub[0,1][2,2][1,2,1]/Sqrt[2]

e[\[Nu]][0,0]:=tub[0,2][0,2][1,0,1]
e[\[Nu]][1,0]:=tub[0,2][1,2][1,1,1]
e[\[Nu]][2,0]:=tub[0,2][2,0][1,2,1]
e[\[Nu]][3,0]:=tub[0,2][2,1][1,2,1]
e[\[Nu]][4,0]:=tub[0,2][2,2][1,2,1]/2^(1/4)


(* ::Subsubsection:: *)
(*Embedding maps*)


VEmbedding[1->1\[CircleTimes]1,1]={{1,0,0},{0,1,0},{0,0,Sqrt[2]}};
VEmbedding[\[Mu]->1\[CircleTimes]\[Mu],1]={{1,0,0},{0,1,0},{0,0,Sqrt[2]}};
VEmbedding[\[Nu]->1\[CircleTimes]\[Nu],1]={{1,0,0,0,0},{0,1,0,0,0},{0,0,Sqrt[2],0,0},{0,0,0,Sqrt[2],0},{0,0,0,0,Sqrt[2]}};

VEmbedding[\[Mu]->\[Mu]\[CircleTimes]1,1]={{1,0,0},{0,1,0},{0,0,Sqrt[2]}};
VEmbedding[1->\[Mu]\[CircleTimes]\[Mu],1]={{1,0,0},{0,1,0},{0,0,Sqrt[2]}};
VEmbedding[\[Nu]->\[Mu]\[CircleTimes]\[Nu],1]={{1,0,0,0,0},{0,1,0,0,0},{0,0,Sqrt[2],0,0},{0,0,0,Sqrt[2],0},{0,0,0,0,-Sqrt[2]}};

VEmbedding[\[Nu]->\[Nu]\[CircleTimes]1,1]={{Sqrt[2],0,0,0,0},{0,Sqrt[2],0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,Sqrt[2]}};
VEmbedding[\[Nu]->\[Nu]\[CircleTimes]\[Mu],1]={{Sqrt[2],0,0,0,0},{0,Sqrt[2],0,0,0},{0,0,0,1,0},{0,0,1,0,0},{0,0,0,0,-Sqrt[2]}};

VEmbedding[1->\[Nu]\[CircleTimes]\[Nu],1]={{1/Sqrt[2],0,0},{0,0,0},{0,0,0},{0,0,0},{0,1/Sqrt[2],0},{0,0,0},{0,0,1/4},{0,0,1/4},{0,0,0},{0,0,0},{0,0,1/2}};
VEmbedding[\[Mu]->\[Nu]\[CircleTimes]\[Nu],1]={{0,0,0},{1/Sqrt[2],0,0},{0,0,0},{0,1/Sqrt[2],0},{0,0,0},{0,0,0},{0,0,1/4},{0,0,1/4},{0,0,0},{0,0,0},{0,0,-(1/2)}};
VEmbedding[\[Nu]->\[Nu]\[CircleTimes]\[Nu],1]={{0,0,0,0,0},{0,0,0,0,0},{1,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,-1,0,0,0},{0,0,0,0,1/2},{0,0,0,0,-(1/2)},{0,0,1,0,0},{0,0,0,-1,0},{0,0,0,0,0}};
