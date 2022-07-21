(* ::Package:: *)

(* ::Subsubsection:: *)
(*General data*)


reps={1,\[Psi]};
dimVS[1]:=1;dimVS[\[Psi]]:=1;

repNames[x_]/;MemberQ[reps,x]:=x


(* ::Subsubsection:: *)
(*Representations*)


e[1][0,0]:=1/2 (tub[0,0][0,0][1,0,1]+tub[0,0][0,0][1,1,1]);

e[\[Psi]][0,0]:=1/2 (tub[0,0][0,0][1,0,1]-tub[0,0][0,0][1,1,1]);


(* ::Subsubsection:: *)
(*Embedding maps*)


(* ::Input::Initialization:: *)
VEmbedding[1->1\[CircleTimes]1,1]={{2^(1/4)}};

VEmbedding[\[Psi]->1\[CircleTimes]\[Psi],1]={{2^(1/4)}};

VEmbedding[\[Psi]->\[Psi]\[CircleTimes]1,1]={{2^(1/4)}};

VEmbedding[1->\[Psi]\[CircleTimes]\[Psi],1]={{2^(1/4)}};
