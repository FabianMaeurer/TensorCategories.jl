PacletDirectoryLoad["~/Projects"];

<<Anyonica`

Print["Loaded Anyonica"];

dir = "~/Projects/MultFreeCenters/Mathematica/FullGaugeTransforms/GaugeSplitTransforms";

fileNames = FileNames[ All, dir ];

memGST[ x__ ] := 
  memGST[ x ] = 
  MapAt[ SparseArray, GaugeSplitTransform[x], {1} ];

Do[
  fn = FileNameJoin[ { dir, "transform_FCL_"<>ToString[i]<>".wdx" } ];
  If[ 
    (* FreeQ[fn] @ fileNames *) True
    ,
    Print["Exporting data for FCL[["<>ToString[i]<>"]]"];
    cat       = FCL[[i]];
    ring      = FusionRing @ cat; 
    zeros     = Keys @ Select[ FSymbols @ cat, #[[2]] === 0& ];
    include   = 
      If[ 
        BraidedQ @ cat, 
        { "FSymbols", "RSymbols", "PSymbols" }, 
        { "FSymbols", "PSymbols" }  
      ];
  
    gst = memGST[ ring, "Zeros" -> zeros, "IncludeOnly" -> include ]

    Export[ fn, gst ]
    ,
    Print["Data already exists for FCL[["<>ToString[i]<>"]]"];
  ]
  ,
  { i, {902} }
]
