(* ::Package:: *)

(* ::Subsubsection:: *)
(*General data*)


reps={1,\[Psi],\[DoubledPi]};
dimVS[1]:=1;dimVS[\[Psi]]:=1;dimVS[\[DoubledPi]]:=2;(*Vector space dimensions*)

repNames[x_]/;MemberQ[reps,x]:=x


(* ::Subsubsection:: *)
(*Representations*)


e[1][0,0]:={1/6,1/6,1/6,1/6,1/6,1/6} . pictureBasis[]
e[\[Psi]][0,0]:={1/6,1/6,1/6,-(1/6),-(1/6),-(1/6)} . pictureBasis[]


e[\[DoubledPi]][0,0]:={1/3,-(1/6),-(1/6),-(1/3),1/6,1/6} . pictureBasis[]
e[\[DoubledPi]][1,0]:={0,1/(2 Sqrt[3]),-(1/(2 Sqrt[3])),0,-(1/(2 Sqrt[3])),1/(2 Sqrt[3])} . pictureBasis[]


(* ::Subsubsection:: *)
(*Embedding maps*)


VEmbedding[1->1\[CircleTimes]1,1]={{6^(1/4) }};
VEmbedding[\[Psi]->1\[CircleTimes]\[Psi],1]={{6^(1/4) }};
VEmbedding[\[DoubledPi]->1\[CircleTimes]\[DoubledPi],1]={{6^(1/4) ,0},{0,6^(1/4) }};

VEmbedding[\[Psi]->\[Psi]\[CircleTimes]1,1]={{6^(1/4) }};
VEmbedding[1->\[Psi]\[CircleTimes]\[Psi],1]={{6^(1/4) }};
VEmbedding[\[DoubledPi]->\[Psi]\[CircleTimes]\[DoubledPi],1]={{0,-6^(1/4) },{6^(1/4),0}};

VEmbedding[\[DoubledPi]->\[DoubledPi]\[CircleTimes]1,1]={{6^(1/4) ,0},{0,6^(1/4) }};
VEmbedding[\[DoubledPi]->\[DoubledPi]\[CircleTimes]\[Psi],1]={{0,-6^(1/4) },{6^(1/4) ,0}};

VEmbedding[1->\[DoubledPi]\[CircleTimes]\[DoubledPi],1]={{1/2 (3/2)^(1/4) },{0},{0},{1/2 (3/2)^(1/4)}};
VEmbedding[\[Psi]->\[DoubledPi]\[CircleTimes]\[DoubledPi],1]={{0},{1/2 (3/2)^(1/4) },{-(1/2) (3/2)^(1/4) },{0}};
VEmbedding[\[DoubledPi]->\[DoubledPi]\[CircleTimes]\[DoubledPi],1]={{0,3^(1/4)/2^(3/4)},{3^(1/4)/2^(3/4),0},{(3^(1/4)) /2^(3/4),0},{0,Root[-3 + 8 #^4& , 1, 0] }};
