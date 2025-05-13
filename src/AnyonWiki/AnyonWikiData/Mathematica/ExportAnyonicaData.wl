(* ::Package:: *)

(*
	The first export function will export numbers to the formal 
	used by Calcium. 
	This is saddly the least optimal implementation of roots of
	polynomials but at the moment the only one for which I can 
	guarantee that we can export all data.
	
	Later we should add functions that export to cyclotomics. 
*)


(* Load Anyonica *)
PacletDirectoryLoad["~/Projects"];
<<Anyonica`


(* ::Section::Closed:: *)
(*Converting lists of numbers to algebraic numbers*)


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




SetAttributes[ safeSimplify, Listable ];
safeSimplify[ f_ -> x_ ] :=
	With[ { newx = FullSimplify[ x ] }, 
		If[
			InfN[x,100] === InfN[newx,100],
			f -> newx,
			f -> x
		]	
	];


(*
indices = {392,393,396,397,684,686};	

Monitor[
	Table[
		Module[
			{fs, nums,invariants,zeros,gis,ring,cat,fgis,rgis,algfgis,algrgis},
			cat = FCL[[i]];
			ring = FusionRing @ cat;
			fs = safeSimplify @ FSymbols @ cat;
			zeros = Keys @ Select[ fs, #[[2]] == 0 & ];
			fgis = 
				Union @ 
				ReplaceAll[ Dispatch @ Table[ frule[[1]] -> memoizedToNumberField[frule[[2]]], { frule, fs } ] ] @ 
				GaugeInvariants[ ring, "IncludeOnly" -> "FSymbols", "Zeros" -> zeros ]; 
			algfgis = ParallelTable[ ToNumberField[ x, All ], { x, fgis } ];
			
			If[ 
				BraidedQ @ cat, 
				rgis = 
					Union @ ReplaceAll[ Dispatch @ Table[ rrule[[1]] -> memoizedToNumberField[rrule[[2]]], { rrule, safeSimplify @ RSymbols @ cat } ] ] @ GaugeInvariants[ ring, "IncludeOnly" -> "RSymbols" ];
				algrgis = ParallelTable[ ToNumberField[ x, All ], { x, rgis } ];
			]; 			
			
			Export[
				"/home/gert/Projects/MultFreeCenters/Mathematica/NumberField"<>ToString[i]<>".wdx"
				,
				{ 
					i,
					FormalCode @ cat,
					AbsoluteTiming @ ListToNumberField[ algfgis ],
					If[
						BraidedQ @ cat, 
						AbsoluteTiming @ ListToNumberField[ algrgis ],
						Missing["NonBraidedCategory"]
					]
				}
			]
		],
		{i, indices }
	],
	i 
];
*)


(* ::Section::Closed:: *)
(*Write the F-symbols in our new fields (Breaks Unitarity)*)


numberFieldPols = 
	Quiet[
		Prepend[1]@ (* I appended 1 at some point rather than prepending it *)
		Table[ 
			Import[ "/home/gert/Projects/MultFreeCenters/Mathematica/simplerfields"<>ToString[i]<>".txt" ], 
			{ i, 1, Length @ FCL } 
		] /. s_String :> ToExpression[ StringReplace[ s, "\n" -> "" ] ]/.$Failed -> 1
	][[;;-2]];


PowerDot[ a_, b_ ] :=
  If[ 
    MatchQ[ a, { 1 .. } ],
    ConstantArray[1,Length[b]],
    Inner[ Power, a, Transpose @ b, Times ]
  ];
ClearAll[TransposeInverse]
TransposeInverse[mat_] := TransposeInverse[mat] = Transpose[Inverse[mat]];


indices = {392,393,396,397,684,686};

fSymbols = 
Monitor[
	Table[ 
		FailOnMessage[ 
		Block[ { transfMat, cat, ring, zeroFs, n, gis, m, fs, newFs, rt, gvs,algebraicFs, $IterationLimit = 2^20,validSolution },
			cat = FCL[[i]];
			ring = Echo @ FusionRing @ FCL[[i]];
			fs = FSymbols @ cat;
			zeroFs = Keys @ Cases[ FSymbols @ cat, rule_/; rule[[2]] == 0 ];
			{ transfMat, n } = GaugeSplitTransform[ ring, "IncludeOnly" -> "FSymbols", "Zeros" -> zeroFs ];
			{ gis, gvs } = GaugeSplitBasis[ ring, "IncludeOnly" -> "FSymbols", "Zeros" -> zeroFs ];

			newFs =
				PowerDot[ gis ~ Join ~ ConstantArray[ 1, Length @ gvs], TransposeInverse @ transfMat ];

			validSolution = 
				CheckPentagonEquations[ 
					ring,
					Thread[ Keys @ fs -> (newFs/.Dispatch[fs]) ], 
					"PreEqualCheck" -> ( InfN[ #, 100 ]& ) 
				];
			
			If[ 
				!TrueQ[validSolution], 
				Return[ 
					gvs/.\[ScriptCapitalF][i__]:>GaugeTransform[u][\[ScriptCapitalF][i]]
				]
			];
			
			algebraicFs = 
			If[ 
				numberFieldPols[[i]] === 1
				,
				Thread[ FSymbols[ring] -> (newFs/.Dispatch[fs]) ]
				,
				rt =
					(ReverseSortBy[ { RealValuedNumberQ @ #, If[ RealValuedNumberQ[#], N[#], {Abs[N[#]], Arg[N[#]] } ] }& ] @ 
					Values @ 
					Flatten @ 
					Solve[ numberFieldPols[[i]] == 0, x ]
					)[[2]];
				Thread[ FSymbols[ring] -> Table[ ToNumberField[ f, rt ], { f, newFs/.Dispatch[fs] } ] ]
			];
			
			Export["~/FSymbols_"<>ToString[i]<>".wdx", algebraicFs ]
		],
		Echo[i]
		]
		,
		{ i, probPos } 
	], 
	Column[ { ProgressIndicator[i/Length[FCL]], i } ]
]


fSymb =
	Table[
		Import["/home/gert/Projects/MultFreeCenters/Mathematica/FSymbols/FSymbols_"<>ToString[i]<>".wdx"],
		{ i, Length @ FCL }
	];


truthTab = Monitor[ 
	Table[ 
		{ i, CheckPentagonEquations[ FusionRing @ FCL[[i]], fSymb[[i]], "PreEqualCheck" -> (InfN[#,64]&) ] },
		{ i, probPos }
	], 
	ProgressIndicator[ i/Length[FCL] ]
]


(* ::Section::Closed:: *)
(*Finding the gauge transforms between old F-symbols and new ones*)


gaugeTransforms1 =
	QuietEcho @
	Monitor[
		Table[ 
			FailOnMessage[ 
				gt[ 
					FusionRing @ FCL[[i]], FSymbols @ FCL[[i]], fSymb[[i]], u, "PreEqualCheck" -> (InfN[ #, 8 ]&)
				],
				Echo[i]
			], 
			{ i, Length @ FCL } 
		],
		ProgressIndicator[i/Length[FCL]]
	];



Export["/home/gert/Projects/MultFreeCenters/Mathematica/NecessaryGaugeTransf_1.wdx", gaugeTransforms1]


failedGTQ[ l_ ] := !MatchQ[l,{__}];


gaugeTransforms2 =
	QuietEcho @
	Monitor[
		Do[ 
			If[
				failedGTQ[ gaugeTransforms2[[i]] ], 
				FailOnMessage[ 
					gaugeTransforms2[[i]] = 
						gt[ 
							FusionRing @ FCL[[i]], FSymbols @ FCL[[i]], fSymb[[i]], u
						];
				Print["Fixed ",FCL[[i]]]
			]
			], 
			{ i, Length @ gaugeTransforms1 } 
		],
		Column[ { ProgressIndicator[i/Length[gaugeTransforms1]], i, FCL[[i]] } ]
	];


gaugeTransforms2 = Import["/home/gert/Projects/MultFreeCenters/Mathematica/numericGaugeTransforms.wdx"];


gaugeTransforms2 = Table[ Missing[], { i, Length[gaugeTransforms1] } ];


	QuietEcho @ 
	Monitor[
		Do[ 
			Block[{$MaxExtraPrecision = 2048},
			If[
				MissingQ[gaugeTransforms2[[i]]] || 
				failedGTQ[ gaugeTransforms2[[i]] ], 
				FailOnMessage[ 
					gaugeTransforms2[[i]] = 
						WhichGaugeTransform[ 
							FusionRing @ FCL[[i]], 
							FSymbols @ FCL[[i]], 
							fSymb[[i]], 
							u, 
							"Numeric" -> True, 
							"Accuracy" -> 1028, 
							"PreEqualCheck" -> (InfN[#,64]&) 
						];
				Print["Numerically Fixed ",FCL[[i]]]
				,
				gaugeTransforms2[[i]] = gaugeTransforms1[[i]]
			]
			]
			], 
			{ i, Length @ gaugeTransforms1 } 
		],
		Column[ { ProgressIndicator[i/Length[gaugeTransforms1]], i, FCL[[i]] } ]
	];


Export["/home/gert/Projects/MultFreeCenters/Mathematica/numericGaugeTransforms.wdx", gaugeTransforms2 ]


gaugeTransforms2 = Import["/home/gert/Projects/MultFreeCenters/Mathematica/numericGaugeTransforms.wdx"];


numPos = Flatten @ Position[ gaugeTransforms2,gTrans_/; Or@@(InexactNumberQ/@Values[gTrans]),1,Heads->False];


rtApproximants = 
ParallelTable[ 
	MapAt[ RootApproximant, gaugeTransforms2[[i]], { All, 2 } ]
	,
	{ i, numPos }
];


properApproximants = Table[ True , { Length @ FCL }];
properApproximants[[numPos]] = ConstantArray[ Missing[], Length @ numPos ];


	Monitor[
		Do[
			If[ 902 > numPos[[i]] > 897 && Mod[ i, 3 ] == 0, CloseKernels[]; LaunchKernels[] ];
			If[ 
			MissingQ[ properApproximants[[numPos[[i]]]] ],
			reduced = 
				(GaugeTransform[u]/@FSymbols[FusionRing @ FCL[[numPos[[i]]]]])/.
				Dispatch[rtApproximants[[i]]]/.
				Dispatch[FSymbols[FCL[[numPos[[i]]]]]];
			newFVals = Values @ fSymb[[numPos[[i]]]];
			diff = Numerator @ Together[ reduced - newFVals ];
			properApproximants[[numPos[[i]]]] = (ParallelDo[ If[ RootReduce[d] =!= 0, Throw["Problem"] ], { d, diff} ]/.{$Aborted->False,Null->True});
			Export["/home/gert/Projects/MultFreeCenters/Mathematica/ChecksCorrectGT/"<>ToString[numPos[[i]]]<>":"<>ToString[properApproximants[[numPos[[i]]]]]<>".txt",""]
			]
			,
			{ i, Length[numPos]/2+.5, Length[numPos] }
		],
		Column[{ ProgressIndicator[i/Length[numPos]],i,FCL[[numPos[[i]]]]}]
	];


Export["/home/gert/Projects/MultFreeCenters/Mathematica/gaugeChecks.wdx",properApproximants]


(* Some of these were not recognized as proper gauge transforms because their F-symbols 
contained Abs[Root[...]] and RootReduce does not reduce such expressions. They are 
actually proper (use SafeFullSimplify) *)
checks = Import["/home/gert/Projects/MultFreeCenters/Mathematica/gaugeChecks.wdx"];
Flatten @ Position[checks,False]


gaugeTransforms3 = gaugeTransforms2;
gaugeTransforms3[[numPos]] = rtApproximants;


Export["/home/gert/Projects/MultFreeCenters/Mathematica/GaugeTransforms.wdx", gaugeTransforms3 ];


(* ::Section:: *)
(*Setting up full data of category*)


Quit[]


PacletDirectoryLoad @ "~/Projects";
<<Anyonica`

(* Import gauge split transforms *)
gtdir = "/home/gert/Projects/MultFreeCenters/Mathematica/FullGaugeTransforms/";
igtdir = "/home/gert/Projects/MultFreeCenters/Mathematica/FullGaugeTransforms/";

gsts = Import[gtdir<>"/AllTransforms.wdx"];
igsts = Import[ igtdir <> "/AllInvTransforms.wdx" ];


expdir = "/home/gert/Projects/MultFreeCenters/Mathematica/FullFields/Mathematica/NewSymbols/";


ClearAll[FixSparseArray]
FixSparseArray[ array_ ] :=
	With[ { dims = Dimensions @ array}, 
		SparseArray[ ArrayRules[array]/.s_String->1, dims ]
	]


(* Helper functions *)

ClearAll[applyTransform];
applyTransform[ { V_, n_ }, symbols_ ] := 
  applyTransform[ { V, n }, symbols ] =
  With[ { nV = Normal @ V }, 
    Map[ PowerDot[ symbols, Transpose[#] ]&,
      { nV[[;;,;;n]], nV[[;;,n+1;;]] }
    ]
  ];

ClearAll[getSymbols];
getSymbols[cat_] :=
  If[ 
    BraidedQ @ cat,
    Comap[ { FSymbols, PSymbols, RSymbols } ]
    ,
    Comap[ { FSymbols, PSymbols } ]
  ];    

PSymbols[r_FusionRing] := 
	Array[ \[ScriptP], Rank @ r ];
PSymbols[c_FusionCategory]:=
	PivotalStructure@c; 
	
(* Code to export values using field extensions  *)

ClearAll[ListToNumberField]
(* Divide an conquer method for converting list of numbers to 
	 algebraic numbers over the same field. This is much faster
	 than applying ToNumberField to the whole list *)
ListToNumberField[l_] := 
	If[ 
		Length[l] == 1
		,
		memoizedToNumberField[ l ]
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

\[Zeta][n_] := ToNumberField[Exp[ 2 Pi I / n ]];	
	
TryCyclotomics[ vals_, opts:OptionsPattern[] ] :=
	Module[{ i = 3, memToNumberField, result, uniqueNontrivialVals, an },
		result = $Failed;
		uniqueNontrivialVals = 
			SortBy[ 
				Cases[ 
					DeleteDuplicates @ vals, 
					x_/; 
					!(Head[x]===Integer) && 
					!(Head[x]===Rational)
				], 
				{ Head[#]===Root, LeafCount[#] }&
			];
		
		While[ 
			i < 128, 
			PrintTemporary[i];
			an = \[Zeta][i];
			result = FailOnMessage[ {  ToNumberField[#,an]& /@ uniqueNontrivialVals , i }, $Failed ];
			If[ result =!= $Failed, Break[] ];
			i++;
		];
		
		If[ result === $Failed, Return @ $Failed ];
	
		memToNumberField[ x_ ] := 
			memToNumberField[x] = 
			ToNumberField[ x, an ];
		
		{ memToNumberField /@ vals, i }	
	];
	
power = Function[ { x, y }, If[ x === 0 && y === 0, 1, x^y ] ];

SafePowerDot[ a_, b_ ] :=
	Inner[ power, a, Transpose[b], Times ];
	
SetAttributes[ safeSimplify, Listable ];
safeSimplify[ f_ -> x_ ] :=
	With[ { newx = FullSimplify[ x ] }, 
		If[
			InfN[x,100] === InfN[newx,100],
			f -> newx,
			f -> x
		]	
	];
	
fixRoots[expr_] := ReplaceAll[ expr, Root[a_,b_,c_]:> Root[a,b] ];

checkSymbols[i_,symb_]:=
Module[ { check1, check2, check3, bq, fsymb,rsymb, ring, psymb,vs },
	If[symb === $Failed, Return @ $Failed];
	ring  = FusionRing @ FCL[[i]];
	fsymb = FilterFRules[symb];
	rsymb = FilterRRules[symb];
	psymb = FilterRules[ symb, \[ScriptP][_] ];
	If[ !TrueQ[ CheckPivotalEquations[ ring, fsymb, psymb ] ], Return @ {False,"Pivotal"} ];
	If[ !TrueQ[ vs = CheckPentagonEquations[ ring, fsymb, "PreEqualCheck" -> (InfN[#,64]&) ] ], Return @ {False,"Pentagon", vs[[2]]} ];
	If[ rsymb=!= {} && !TrueQ[ vs = CheckHexagonEquations[ ring, fsymb, rsymb, "PreEqualCheck" -> (InfN[#,64]&) ] ], Return @ {False,"Hexagon",vs[[2]]} ];
	True
];


fns = FileNames[ All, expdir ];

MemoryConstrained[
Do[
	fn = expdir<>"new_symbols_FCL_"<>ToString[i]<>".wdx";
	If[ 
		FreeQ[ fns, fn ],		
		prevCat = FCL[[i-1]];
    cat     = FCL[[i]];
		Print["Finding field for FCL[["<>ToString[i]<>"]]: ", FCL[[i]]];
    newRingQ = FC[prevCat][[;;4]] =!= FC[cat][[;;4]];
    
    ring    = FusionRing @ cat;
    symbols = Join @@ getSymbols[cat] @ ring;
    
    Print["Calculating values of symbols"];
    values  = DynamicMap[ #[[1]]->ToNumberField[#[[2]],All]&, fixRoots[ Union[ Join @@ getSymbols[cat] @ cat ] ] ];

    (* set up gauge split basis *)
    symbbasis = applyTransform[ gsts[[i]], symbols ]; 
    gsb       = symbbasis/.Dispatch[values];
    
    Print["Simplifying values of symbols"];
    (* Find the field extension *)
		If[ newRingQ, ClearAll[memRootReduce] ];
		memRootReduce[ x__ ] := memRootReduce[x] = RootReduce[x];    
    
    reducedgsb = fixRoots[ DynamicMap[ RootReduce, gsb[[1]] ] ];
    
    If[ 
    	(*FreeQ[ Range[398,461], i ] && (* TriCritIsing *) 
			FreeQ[ Range[548,572], i ] &&*) (* Fib x PSU : no cyclotomics found with rank < 500! *) 
			FreeQ[ Range[683,686], i ] (*&&
			FreeQ[ Range[860,863], i ] && 
			FreeQ[ Range[880,887], i ]*),
			symbdegree = 
				If[ 
					newRingQ,
					TryCyclotomics[ reducedgsb ],
					TryCyclotomics[ reducedgsb, "StartingValue" -> cyclodegree ]
				]
			,
			symbdegree = $Failed
		];
		
		If[ 
			symbdegree === $Failed, 
			Print["Gauge-split basis of FCL[["<>ToString[i]<>"]] can't be expressed in cyclotomic field of small size."];
			fieldSymbols = ToNumberField[#,All]& /@ reducedgsb,
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
		
		newSymbols = newsymb/.vacUVals/.u[_]->1;
		If[ 
			TrueQ[checkSymbols[i,newSymbols]],
			Export[fn,newsymb/.vacUVals/.u[_]->1],
			Print[Style["Equations not satisfied for i = "<>ToString[i],Red]]
		]
	]                                                                                                             
  , { i, 895, 901 } 
], 100 * 10^9
]



(*
ISSUES:
683::686 : Haagerup
{37,38,39,40,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,880,881,882,883,884,885,886,887}:
	symbols are not in proper form due to bug with rootreduce and tonumberfield
*)


(* ::Section::Closed:: *)
(*Putting the R Symbols in the new gauges *)


gaugeTransforms = Import["/home/gert/Projects/MultFreeCenters/Mathematica/GaugeTransforms.wdx"];


newRSymbols =
	Monitor[
	Table[ 
		cat = FCL[[i]];
		rs = RSymbols @ cat;
		Thread[
			Keys @ rs ->
			ToNumberField @ 
			RootReduce[
				GaugeTransform[ u ] /@ Keys @ rs/.
				Dispatch[gaugeTransforms[[i]]]/.
				Dispatch[rs]
			]
		]
		,
		{ i, 20 }
		],
		ProgressIndicator[i/Length[FCL]]
	];


Export["/home/gert/Projects/MultFreeCenters/Mathematica/NewPivots.wdx",newPivots];


newPivots = Import["/home/gert/Projects/MultFreeCenters/Mathematica/NewPivots.wdx"];


(* ::Section:: *)
(*Testing whether the symbols are still valid*)


checkSymbols[i_,symb_] :=
Module[ { check1, check2, check3, bq, fsymb,rsymb, ring, psymb,vs },
	If[symb === $Failed, Return @ $Failed];
	ring  = FusionRing @ FCL[[i]];
	fsymb = FilterFRules[symb];
	rsymb = FilterRRules[symb];
	psymb = FilterRules[ symb, \[ScriptP][_] ];
	If[ !TrueQ[ CheckPivotalEquations[ ring, fsymb, psymb ] ], Return @ {False,"Pivotal"} ];
	If[ !TrueQ[ vs = CheckPentagonEquations[ ring, fsymb, "PreEqualCheck" -> (InfN[#,64]&) ] ], Return @ {False,"Pentagon", vs[[2]]} ];
	If[ rsymb=!= {} && !TrueQ[ vs = CheckHexagonEquations[ ring, fsymb, rsymb, "PreEqualCheck" -> (InfN[#,64]&) ] ], Return @ {False,"Hexagon",vs[[2]]} ];
	True
];
checkSymbols[i_] := checkSymbols[ i, symbols[[i]] ];


checks = QuietEcho@ DynamicMap[ checkSymbols, Range[1,300] ];
probpos = Position[checks,{False,__},1]//Flatten;


FCL[[451]]


(* ::Section:: *)
(*Exporting data to Julia files*)


(* ::Text:: *)
(*We will export the F-symbols, R-symbols and P-symbols  to a format Julia can read. *)
(*There are 2 scenarios per set of symbols *)
(*1. All are expressed in a cyclotomic field. In this case we export the numeric values and the degree (n) of the cyclotomic numbers, i.e. the n'th root of unity that allows us to express the solutions*)
(*2. Not all numbers are expressed in a cyclotomic field. In this case we export the numeric values and the highest degree of the polynomials required to express the solution*)


symbols = (* DynamicMap[Import["/home/gert/Projects/MultFreeCenters/Mathematica/FullFields/Mathematica/NewSymbols/new_symbols_FCL_"<>ToString[#]<>".wdx"]&,Range @ Length @ FCL];*)
	Import["/home/gert/Projects/MultFreeCenters/Mathematica/FullFields/Mathematica/AllNewSymbols.mx"];


Export["/home/gert/Projects/MultFreeCenters/Mathematica/FullFields/Mathematica/AllNewSymbols.wdx",symbols];


Table[Export["/home/gert/Projects/MultFreeCenters/Mathematica/FullFields/Mathematica/NewSymbols/new_symbols_FCL_"<>ToString[i]<>".wdx",symbols],{i,572,573}]


(* ::Subsubsection:: *)
(*Determining which symbols are expressed using cyclotomics*)


numberField[symb_] := DeleteDuplicatesBy[ Cases[ symb, _AlgebraicNumber, Infinity ], First ];

cycloPols = Association @@ Table[ Cyclotomic[i,x]->i, { i, 3, 127 } ]; 

cycloDegree[ a_AlgebraicNumber ] :=
	With[ { pol = a[[1]] }, 
		If[
			NumericQ[pol], 
			cycloPols[ MinimalPolynomial[pol][x] ],
			cycloPols[pol]
		]
	];

CyclotomicQ[ symb_ ] := 
	With[{ nf = numberField[symb] }, 
		!MemberQ[Values@symb,_Root]&&
		properSymbolsQ[symb] && (
			Length[nf] === 0 ||
			Length[nf] === 1 && IntegerQ[ cycloDegree[ First @ nf ] ]
		)
	];

properSymbolsQ[ symb_ ] :=
	With[ {vls = Values @ symb},
		Catch[
			Do[ If[ MatchQ[vl,HoldPattern[Times[___,_Root,___]|Times[___,_AlgebraicNumber,___]]], Throw @ False ] , {vl, vls} ];
			True
		] 
	];


cycloPos = Position[ symbols, x_/; CyclotomicQ[x], 1, Heads->False ]//Flatten;	
nonCycloPos = Complement[Range[977],cycloPos]
(* nonCycloPos = Flatten @ Position[symbols,s_/;s=!=$Failed && MemberQ[Values @ s,HoldPattern[ Times[ AlgebraicNumber[x__],AlgebraicNumber[y__],___] ] ],1,Heads->False ]*)


Export["/home/gert/Projects/MultFreeCenters/Mathematica/FullFields/Mathematica/cycloPos.txt",StringRiffle[ToString/@cycloPos,"\n"] ]


(* ::Subsection:: *)
(*Export the cyclotomics*)


symbols[[4]]


degree[ symb_ ] := 
	Module[ { nf },
		nf = numberField @ symb;
		If[ Length @ nf === 0, Return @ 0 ];
		If[ MatchQ[ First @ nf, _AlgebraicNumber], cycloDegree[ First @ nf ] ]
	]


cycloData = Transpose[ { cycloPos, symbols[[cycloPos]], degree/@symbols[[cycloPos]] }];


(* Sometimes the generator of the field changed to another unit root >.<  *)
findPower[ a_AlgebraicNumber ] := Numerator[If[ !MatchQ[a,_Complex|_Integer|_Rational], Rationalize[InfN[ Arg[a[[1]]]/(2Pi), 16]] ]];

cycloString[n_][ a_Rational ] := numString[a] ;
cycloString[n_][ a_Integer ]  := numString[a] ;

cycloString[n_][ a_AlgebraicNumber ] :=
	If[
		MatchQ[ a[[1]], _Complex ],
		complexCycloString[n][a],
		With[{pow = findPower @ a },
			StringRiffle[ DeleteCases[ listToCyclo[n,pow][a[[2]]], "0" ], " + " ]
		]
	];
	
listToCyclo[n_,pow_][ l_List ] := 
	combine[n] /@ Transpose[ { pow *( Range @ Length[l] - 1 ), l } ];

(* Construct string for i'th term *)
combine[n_][{i_,val_}] :=
	If[ 
		i == 0, 
		numString[val],
		StringReplace[ If[ val == 0, "0", numString[val]<>" * (z"<>ToString[n]<>")^"<>ToString[i] ], "1 * " -> "" ]
	]; 

complexCycloString[ n_ ][ z_ ] := 
	Module[ { a, b , l, list },
		{ a, b } = Last @ z; 
		l = n/4+1; 
		If[ !IntegerQ[l], Throw["Problem"]];
		list = ConstantArray[0,n];
		list[[1]] = a;
		list[[l]] = b;
		list;
		StringRiffle[ DeleteCases[ listToCyclo[n][list], "0" ], " + " ]
	]
	

numString[a_Rational] := ToString[Numerator[a]] <> "//" <> ToString[Denominator[a]]; 
numString[a_Integer]  := ToString[a];


ToDict[str_][ data_ ] :=
	Module[
		{ symbols, degr, cStrings, indStrings, contents }, 
		symbols = 
			Switch[ str,
				"F", FilterFRules @ data[[2]],
				"R", FilterRRules @ data[[2]],
				"P", Cases[ data[[2]], HoldPattern[ \[ScriptP][a_] -> b_ ] ]
			]; 
		degr = data[[3]];
		cStrings = cycloString[degr] /@ Echo @Values[symbols];
		indStrings = (StringReplace[#,{ "{"->"[", "}"->"]" }]&)@*ToString /@ ( List @@@ Keys[symbols] );
		contents = StringRiffle[ MapThread[ StringJoin[ #1, " => ", #2 ]&, { indStrings, cStrings } ], ", " ];
		"# zn = exp( 2 pi i / n )\n\n"<>  
		"Dict(\n\t"<>contents <>"\n)"
	];


z16 = Exp[2 Pi I / 16 ];
RootReduce[ z16^4 - AlgebraicNumber[Root[1 + #^8& , 8, 0], {0, 0, 0, 0, 1, 0, 0, 0}]]


QuietEcho @ 
Monitor[Do[ 
	Export["/home/gert/Projects/MultFreeCenters/cyclic_6j_Symbols/cat_"<>ToString[ cycloData[[i,1]] ]<> ".jl", ToDict["F"][cycloData[[i]]] , "Text" ], 
	{ i, 246, Length @ cycloData } 
],i]


QuietEcho @ 
Monitor[
	Do[
		If[ BraidedQ[ FCL[[cycloData[[i,1]] ]]],
	  Export[
	  "/home/gert/Projects/MultFreeCenters/cyclic_R_Symbols/cat_"<>ToString[ cycloData[[i,1]] ]<> ".jl", ToDict["R"][cycloData[[i]]] , 
	  "Text" 
	  ]
	  ], 
	{ i, 1, Length @ cycloData } 
],i
]


QuietEcho @ 
Monitor[
Do[
	  Export[
	  "/home/gert/Projects/MultFreeCenters/cyclic_P_Symbols/cat_"<>ToString[ cycloData[[i,1]] ]<> ".jl", ToDict["P"][cycloData[[i]]] , 
	  "Text" 
	  ],
	{ i, 1, Length @ cycloData } 
],i
]


(* ::Subsection:: *)
(*Non-cyclotomic data*)


(* First we get the highest degrees of the polynomials and bundle all info together *)


maxDegree[ symb_ ] := 
	Max @ 
	Cases[ 
		Comap[ First /@ Union[ Cases[ symbols[[nonCycloPos[[3]]]], _Root, Infinity ]] ][x], 
		Power[x,a_]:>a,
		Infinity 
	];
	
nonCycloData = 
	Transpose[{ nonCycloPos, symbols[[nonCycloPos]], maxDegree /@ symbols[[nonCycloPos]] }];


formatNumber[x_] :=
	"\""<>StringRiffle[( ToString @* DecimalForm /@ ReIm @ InfN[x,256] ), " + " ] <> "im\"";

nonCycloToDict[str_][ data_ ] :=
	Module[
		{ symbols, degr, numStrings, indStrings, contents }, 
		symbols = 
			Switch[ str,
				"F", FilterFRules @ data[[2]],
				"R", FilterRRules @ data[[2]],
				"P", Cases[ data[[2]], HoldPattern[ \[ScriptP][a_] -> b_ ] ]
			]; 
			
		indStrings = (StringReplace[#,{ "{"->"[", "}"->"]" }]&)@*ToString /@ ( List @@@ Keys[symbols] );
		numStrings = formatNumber /@ Values[symbols];
		contents = StringRiffle[ MapThread[ StringJoin[ #1, " => ", #2 ]&, { indStrings, numStrings } ], ", " ];
		"# "<>ToString[data[[3]] ]<> "\n" <>
		"Dict(\n\t"<>contents <>"\n)"
	];


Export[ "/home/gert/Projects/MultFreeCenters/nonCycloPos.txt", StringRiffle[ First /@ nonCycloData, "\n" ] ]


Monitor[Do[ 
	Export[
		"/home/gert/Projects/MultFreeCenters/non_cyclic_6j_Symbols/cat_"<>ToString[ nonCycloData[[i,1]] ]<> ".jl", 
		nonCycloToDict["F"][nonCycloData[[i]]] , 
		"Text" 
	], 
	{ i, Length @ nonCycloData } 
],i]


Monitor[
	Do[
		If[ BraidedQ[ FCL[[nonCycloData[[i,1]] ]]],
	  Export[
	  "/home/gert/Projects/MultFreeCenters/non_cyclic_R_Symbols/cat_"<>ToString[ nonCycloData[[i,1]] ]<> ".jl", nonCycloToDict["R"][nonCycloData[[i]]] , 
	  "Text" 
	  ]
	  ], 
	{ i, 1, Length @ nonCycloData } 
],i
]


Monitor[
Do[
	  Export[
	  "/home/gert/Projects/MultFreeCenters/non_cyclic_P_Symbols/cat_"<>ToString[ nonCycloData[[i,1]] ]<> ".jl", nonCycloToDict["P"][nonCycloData[[i]]] , 
	  "Text" 
	  ],
	{ i, 1, Length @ nonCycloData } 
],i
]


(* Export the number fields *)
(* Convert a polynomial with variables of the form s[i] to C string *)
PolToString[ i_Integer, s_Symbol ]:=
	ToString[i];
PolToString[ pol_, s_Symbol ] :=
	StringJoin @@
		ToString /@ ({
			pol //
			Expand//
			ReplacePlus//
			ReplaceTimes//
			ReplacePower//
			ReplaceBrackets[s]
			}// Flatten
		);
		
ReplacePlus[expr_]:=
	expr/.Plus[a_,b__] :> Riffle[ { a, b }, " + " ];
ReplaceTimes[expr_]:=
	expr/.Times[a_,b__] :> Riffle[{a,b},"*"];
ReplacePower[expr_]:=
	expr/.Power[a_,b__] :> ToString[a]<>"^"<>ToString[b];
ReplaceBrackets[s_Symbol][expr_]:=
	expr/.s[i_] :> SymbolName[s]<>"_"<>ToString[i];

getRtPol[ symbols_ ] := 
	FirstCase[ 
		symbols, 
		AlgebraicNumber[ rt_, _ ] :> { formatNumber[rt], MinimalPolynomial[rt][x] }, 
		{formatNumber[1],1},
		Infinity 
	];

{ rootNums, polStrings } = Transpose @  
	Table[
		{ root, pol } = getRtPol[ fSymb[[i]] ];
		{ root, PolToString[ pol, x ] },
		{ i, Length @ FCL }
	];


Export["/home/gert/Projects/MultFreeCenters/6jSymbols/Roots.dat", rootNums ]
Export["/home/gert/Projects/MultFreeCenters/6jSymbols/Polynomials", StringRiffle[ polStrings, "\n" ], "Text" ]


Export["/home/gert/Projects/MultFreeCenters/cat_properties.dat",	
	"# Braided, Unitary, Spherical, Ribbon, Modular\n" <>
	StringReplace[
		StringRiffle[
			Map[ ToString, Comap[{BraidedQ,UnitaryQ,SphericalQ,RibbonQ,ModularQ}] /@ FCL ],
			"\n"
		],
		{"T"->"t","F"->"f","{"->"","}"->""}
	]
]


(* ::Section::Closed:: *)
(*Legacy code*)


(* F Symbols for PSU(2)_{13} *)
ring = FusionRing @ FCL[[900]];
(* There are no 0 F-symbols for these cats *)
{ sTransMat, n } = 
	EchoPerformance @ 
	SparseArray @ Developer`ToPackedArray @
	GaugeSplitTransform[ ring, "IncludeOnly" -> "FSymbols" ];


invTransMat = Inverse[ sTransMat ];


Dimensions[sTransMat]


{ gis, gvs } = 
	EchoPerformance @ 
	GaugeSplitBasis[ ring, "IncludeOnly" -> "FSymbols" ];


replaceByKnownAlgebraicNumbers[ vals_, algs_ ] := 
	With[ 
		{ algReps = Dispatch @ Thread[ InfN[100] /@ algs -> algs ],
		  valReps = Dispatch @ Thread[ InfN[100] /@ vals -> vals ] },
		ReplaceAll[ InfN[100] /@ vals, algReps ]/.valReps (* If no alg number found: revert to exact value *)
	];


replaceByKnownAlgebraicNumbers[ gis[[;;100]]/.FSymbols[FCL[[900]]], algNums ]


algNums = %[[3,2]];


Monitor[
	Table[ 
		Block[ {  cat, ring, zeroFs, fs, newFs, rt,algebraicFs, $IterationLimit = 2^20, algNums, newFVals },
			cat = FCL[[i]];
			fs = FSymbols @ cat;
			
			algNums = Union @ 
				Import["/home/gert/Projects/MultFreeCenters/Mathematica/NumberField"<>ToString[i]<>".wdx"][[3,2]];
			
			newFVals = EchoPerformance @ 
				replaceByKnownAlgebraicNumbers[
					PowerDot[ gis ~ Join ~ ConstantArray[ 1, Length @ gvs], Transpose @ invTransMat ] // 
					ReplaceAll[ Dispatch @ fs ]
					,
					algNums
				];
			
			algebraicFs = 
			If[ 
				numberFieldPols[[i]] === 1
				,
				Thread[ FSymbols[ring] -> newFVals ]
				,
				rt =
					First @ 
					ReverseSortBy[ { RealValuedNumberQ@#, If[ RealValuedNumberQ[#], N[#], {Abs[N[#]], Arg[N[#]] } ] }& ] @ 
					Values @ 
					Flatten @ 
					Solve[ numberFieldPols[[i]] == 0, x ];
				Thread[ FSymbols[ring] -> ParallelTable[ ToNumberField[ f, rt ], { f, newFVals } ] ]
			];
			Export["/home/gert/Projects/MultFreeCenters/Mathematica/FSymbols_"<>ToString[i]<>".wdx", algebraicFs ]
		],
		{ i, 896, 901 } 
	], 
	Column[ { ProgressIndicator[(i-896)/6], FCL[[i]] } ]
]


algNumbers = (FirstCase[#,_AlgebraicNumber]&/@ allFields[[;;,3,2]])/.Missing[_]->1;


minPol[x_]:= minPol[x] = MinimalPolynomial[x];


minPols = Monitor[ Table[ If[ MatchQ[y,_AlgebraicNumber], minPol[First @ y]@x, y ], { y, algNumbers } ], x ];


minPols


Do[ 
	Module[{polstr,code}, 
	If[ 
		minPols[[i]] =!= 1,
		polstr = PolToString @ minPols[[i]];
		code = 
			StringJoin[ 
				"K,a = number_field("<>polstr<>");\n",
				"K2 = simplify(K)[1];\n",
				"writedlm(\"/home/gert/Projects/MultFreeCenters/Mathematica/simplerfields\" * string("<>ToString[i]<>") * \".txt\"",
				", string(defining_polynomial(K2)));"
			];
		
		Export["/home/gert/Projects/MultFreeCenters/Mathematica/Polynomial"<>ToString[i]<>".txt", 
		code ]
	]
	],	{ i, Length @ minPols }
	
];


(* Indices of nonconstant polynomials *)
indices = Complement[ Range @ Length @ minPols, Flatten @ Position[ minPols, 1 ] ];


(*
Evaluate in julia:
> using Oscar
> using DelimitedFiles
> indices = [9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,32,33,34,35,36,37,38,39,46,47,48,49,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,170,171,172,173,174,175,176,177,178,179,180,181,222,223,224,225,226,227,228,229,233,234,235,236,237,238,241,242,243,244,245,246,247,248,249,250,251,252,253,254,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,526,527,528,529,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,589,590,591,592,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,734,735,736,737,738,739,740,741,742,743,744,745,746,747,748,749,750,751,752,753,754,755,756,757,758,759,760,761,762,763,764,765,766,767,768,769,770,771,772,773,774,775,776,777,778,779,780,781,782,783,784,785,786,787,788,789,790,791,792,793,794,795,796,797,798,799,800,801,861,862,863,864,865,866,867,868,869,870,871,872,873,874,875,876,877,878,879,880,881,882,883,884,885,886,887,888,889,890,891,892,893,894,895,896,897,898,899,900,901,902,903,904,905,906,907,908,909,910,911,912,913,914,915,916,935,936,937,938,939,940,941,942,943,944,945,946,947,948,949,950,951,952,953,954,955,956,957,958,969,970,971,972,973,974,975,976]
> for i in indices
  include(dir*"Polynomial"*string(i)*".txt");
  end
*)


simplerFieldsStrings = 
	Table[ 
		Import["/home/gert/Projects/MultFreeCenters/Mathematica/simplerfields"<>ToString[i]<>".txt"],
		{ i, Complement[ Range @ Length @ minPols, Flatten @ Position[ minPols, 1 ] ] }
	];
simplerFields = 
	Association @ Thread[ indices -> ToExpression@*StringReplace["\n"->""]/@simplerFieldsStrings ];


simplerFields


rt = First @ Values @ Solve[1-4 x+3 x^2+x^3==0,x][[3]]


newGaugeInvariants[i_] := 
	Module[{ rt, cat, ring, fs },
		rt = First @ Values @ Last @ Solve[ simplerFields[i] == 0, x ];
		cat = FCL[[i+1]]; (* +1 because we didn't include trivial cat *)
		ring = FusionRing @ cat;
		fs = FSymbols @ cat;
		ToNumberField[ 
			GaugeInvariants[ ring, "IncludeOnly" -> "FSymbols" ]/.Dispatch[fs],
			rt 
		]
	]


newGaugeInvariants[34]//Length


{M,n} = GaugeSplitTransform[ FusionRing @ FCL[[35]], "IncludeOnly" -> "FSymbols" ]


git = GaugeInvariants[ FusionRing @ FCL[[35]], "IncludeOnly" -> "FSymbols" ]


fst = FSymbols @ FusionRing @ FCL[[35]]


PowerDot[ FSymbols @ FusionRing @ FCL[[35]], Echo @ GaugeSplitTransform[ FusionRing @ FCL[[35]], "IncludeOnly" -> "FSymbols" ][[1]] ]


PowerDot[Join[ newGaugeInvariants[34], Array[u,11] ], Transpose @ Inverse @ M ]
