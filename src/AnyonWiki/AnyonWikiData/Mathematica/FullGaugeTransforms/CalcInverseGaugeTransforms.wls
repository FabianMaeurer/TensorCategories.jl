PacletDirectoryLoad["~/Projects"];
<<Anyonica`

gtdir = "/home/gvercley/Projects/MultFreeCenters/Mathematica/FullGaugeTransforms/GaugeSplitTransforms";
expdir = "/home/gvercley/Projects/MultFreeCenters/Mathematica/FullGaugeTransforms/InverseGaugeSplitTransforms";

(* gsts = Import[gtdir<>"/AllTransforms.wdx"]; *) 

memInv[mat_] := memInv[mat] = SparseArray @ Inverse[mat,Method->"OneStepRowReduction" ]; 

Do[ 
  fn = expdir <> "/inverse_transform_FCL_"<>ToString[i]<>".wdx";
  
  { mat, n } = Import[gtdir<>"transform_FCL_"<>ToString[i]<>".wds"]; 

  Export[ fn, Inverse[mat,Method->"OneStepRowReduction"] ];
  Print["Exported for i = "<>ToString[i] ];

  , { i, {902} }
]
