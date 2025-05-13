PacletDirectoryLoad @ "~/Projects";
<<Anyonica`

Print["Loaded Anyonica"];

(* Import gauge split transforms *)
gtdir  = "/home/gvercley/Projects/MultFreeCenters/Mathematica/FullGaugeTransforms/GaugeSplitTransforms";
igtdir = "/home/gvercley/Projects/MultFreeCenters/Mathematica/FullGaugeTransforms/InverseGaugeSplitTransforms";
gsts   = Import[gtdir<>"/AllTransforms.wdx"];
igsts  = Import[ igtdir <> "/AllInvTransforms.wdx" ];
expdir = "/home/gvercley/Projects/MultFreeCenters/Mathematica/FullFields/Mathematica/NewSymbols/";

Print["Imported Data"];

(* Helper functions *)

ClearAll[applyTransform];
applyTransform[ { V_, n_ }, symbols_ ] := 
  applyTransform[ { V, n }, symbols ] =
  With[ { nV = Normal @ V }, 
    Map[ PowerDot[ symbols, Transpose[#] ]&,
      { nV[[;;,;;n]], nV[[;;,n+1;;]] }
    ]
  ];
  
getSymbols[cat_] :=
  If[ 
    BraidedQ @ cat,
    Comap[ { FSymbols, RSymbols, PSymbols } ]
    ,
    Comap[ { FSymbols, PSymbols } ]
  ];    

PSymbols[r_FusionRing] := 
	Array[ \[ScriptP], Rank @ ring ];
PSymbols[c_FusionCategory]:=
	PivotalStructure@c; 
	
(* Code to export values using field extensions  *)

ClearAll[ListToNumberField]
(* Divide an conquer method for converting list of numbers to 
	 algebraic numbers over the same field. This is much faster
	 than applying ToNumberField to the whole list *)
ListToNumberField[l_] := 
	If[ 
		Length[l] < 5
		,
		ToNumberField[ l, All ]
		,
		ToNumberField @ 
		Join[
			ListToNumberField @ l[[ ;; Floor[ Length[l]/2 ] ]],
			ListToNumberField @ l[[ Floor[ Length[l]/2 ] + 1 ;; ]] 
		]
	];

ClearAll[memoizedToNumberField];	
memoizedToNumberField[x_] := 
	memoizedToNumberField[x] = 
	ToNumberField[x,All];

Options[TryCyclotomics] = 
	{ "StartingValue" -> 2 };
	
TryCyclotomics[ vals_, opts:OptionsPattern[] ] :=
	Module[{ i = 1 + Floor[OptionValue["StartingValue"]/2], an, memToNumberField, result },
		result = $Failed;
		While[ 
			i < 528, 
			an = ToNumberField[(-1)^(1/i)]; 
			ClearAll[memToNumberField];
			memToNumberField[ x_ ] := 
				memToNumberField[x] = 
				ToNumberField[ x, an ];
			result = FailOnMessage[ { Reverse[ memToNumberField /@ Reverse[vals] ], i }, $Failed ];
			If[ result =!= $Failed, i = 528 ];
			i++;
		];
		result	
	];
	
power = Function[ { x, y }, If[ x === 0 && y === 0, 1, x^y ] ];

SafePowerDot[ a_, b_ ] :=
	Inner[ power, a, Transpose[b], Times ];

cyclodegree = 1;

fns = FileNames[ All, expdir ];

Do[
	fn = expdir<>"new_symbols_FCL_"<>ToString[i]<>".wdx";
	If[ 
    FreeQ[ fns, fn ] && FreeQ[ Range[390,461], i ] && FreeQ[ Range[548,572], i ] && FreeQ[ Range[683,686], i ]
  ,		
		prevCat = FCL[[i-1]];
    cat     = FCL[[i]];
		Print["Finding field for FCL[[" <>ToString[i]<>"]]"];
    newRingQ = FC[prevCat][[;;4]] =!= FC[cat][[;;4]];
    
    ring    = FusionRing @ cat;
    symbols = Join @@ getSymbols[cat] @ ring;
    values  = Join @@ getSymbols[cat] @ cat;

    (* set up gauge split basis *)
    symbbasis = applyTransform[ gsts[[i]], symbols ]; 
    gsb       = symbbasis/.Dispatch[values];
    
    (* Find the field extension *)
		If[ newRingQ, ClearAll[memRootReduce] ];
		memRootReduce[ x__ ] := memRootReduce[x] = RootReduce[x];    
    
    reducedgsb = memRootReduce /@ gsb[[1]];
        
		symbdegree = 
			If[ 
				newRingQ,
				TryCyclotomics[ reducedgsb ],
				TryCyclotomics[ reducedgsb, "StartingValue" -> cyclodegree ]
			];
		
		If[ 
			symbdegree === $Failed, 
			Print["Gauge-split basis of FCL[["<>ToString[i]<>" can't be expressed in cyclotomic field of small size."];
			symbdegree = ListToNumberField[ reducedgsb ],
			{ fieldSymbols, cyclodegree } = symbdegree; 
		]; 
	
		newsymb =
		Thread[
			symbols ->
			SafePowerDot[ 
				fieldSymbols ~ Join ~ Array[ u, Length @ gsb[[2]]], 
				Transpose @ igsts[[i]] 
			]
		];
		
		(* Demand vacuum F-s and R-s to be 1 *)
		vacUVals = 
			Quiet[ 
				Solve[ 
					Thread[ (Cases[ symbols, $VacuumFPattern | $VacuumRPattern ]/.Dispatch[newsymb]) == 1 ],
					Array[ u, Length @ gsb[[2]] ] 
				],
				Solve::svars (* Often not all values of u will be fixed this way *)
			];
			
		If[ Length @ vacUVals > 0, vacUVals = First @ vacUVals ];
		
		Export[fn,newsymb/.vacUVals/.u[_]->1 ]
	]                                                                                                             
  , { i, 2, Length @ FCL } 
];

