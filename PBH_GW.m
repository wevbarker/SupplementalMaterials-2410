(*==========*)
(*  PBH_GW  *)
(*==========*)

<<xAct`xPlain`;

Title@"Gravitational wave integration";
$Compute=False;
$TheProcessorCount=100;
If[$Compute,
	Unprotect@$ProcessorCount;$ProcessorCount=$TheProcessorCount;Protect@$ProcessorCount;
	LaunchKernels[$ProcessorCount];
];

AllNames={"1e6g","1e6g_lower","1e6g_upper","1e7g","1e7g_lower","1e7g_upper","1e8g","1e8g_lower","1e8g_upper","1e9g","1e9g_lower","1e9g_upper"};

ProcessSpectrum[InputString_]:=Module[{},
	Code@"First we import the data from the power spectrum file.";
	Code[
		filePath=InputString<>"_PowerSpectrum.txt";
		data=Import[filePath,"Data"];
	];
	Comment@"If you want to extract the individual columns for further processing."
	Code[
		frequencies=data[[All,1]];
		MSPowerSolutions1=data[[All,2]];
	];
	Comment@"For example,you can create a plot."
	Code[
		Expr=ListLogLogPlot[Transpose[{frequencies,MSPowerSolutions1}],Joined->True,AxesLabel->{"f/Hz","h2OmegaGW2(k)"},PlotLabel->"Log-Log Plot of OmegaGW2(k)",PlotRange->{Automatic,Automatic}];
	];
	DisplayExpression@Expr;
	Comment@"Constants and Functions."
	Code[
		MSPowerSolutions=Transpose[{frequencies,MSPowerSolutions1}];
		cg=0.4;
		OmegaR0=2.47*10^-5;
		Ic[d_,s_]:=-36*Pi*((s^2+d^2-2)^2/(s^2-d^2)^3)*UnitStep[s-1];
		Is[d_,s_]:=-36*((s^2+d^2-2)/(s^2-d^2)^2)*((s^2+d^2-2)/(s^2-d^2)*Log[(1-d^2)/Abs[s^2-1]]+2);
	];
	Comment@"Interpolate the power spectrum so we have P as a function of k.";
	Code[
		LogPInterpolation=Interpolation[Log10[MSPowerSolutions],InterpolationOrder->1];
		P[k_]:=10^LogPInterpolation[Log10[k]];
		kMin=Min[MSPowerSolutions[[All,1]]];
		kMax=Max[MSPowerSolutions[[All,1]]];
	];
	Comment@"Plot P[k] over the range of k values.";
	Code[
		Expr=LogLogPlot[P[k],{k,kMin,kMax},PlotRange->All,AxesLabel->{"k","P(k)"},PlotLabel->"PowerSpectrumP(k)",GridLines->Automatic];
	];
	DisplayExpression@Expr;
	Comment@"Calculate GW signature.";
	Code[
		Integrand2[k_,d_,s_]:=((d^2-1/3)*(s^2-1/3)/(s^2-d^2))^2*P[k*Sqrt[3]*(s+d)/2]*P[k*Sqrt[3]*(s-d)/2]*(Is[d,s]^2+Ic[d,s]^2);
		h2OmegaGW2[k_]:=cg*(OmegaR0/36)*(NIntegrate[Integrand2[k,d,s],{d,0,1/Sqrt[3]},{s,1/Sqrt[3],Infinity},WorkingPrecision->40]);
	];
	Comment@"Local Adaptive can often taken a while to run, in which case Adaptive Monte Carlo is recommended for quick but reliable estimates.";
	Comment@"Frequency region to plot GW signature in.";
	Code[
		lowF=1;
		highF=10^6;
		conv=1.546*10^-15;
		k1=10^Range[Log10[lowF/(conv)],Log10[highF/(conv)],(Log10[highF/(conv)]-Log10[lowF/(conv)])/1000];
	];
	DisplayExpression@k1;
	Comment@"Calculate the GW signature for the given frequency range.";
	If[$Compute,
		Code[
			Comment@"The model is:";
			Print@InputString;
			DistributeDefinitions[h2OmegaGW2];
			h2OmegaGW2Sol=Map[
				(ParallelSubmit@(Quiet@h2OmegaGW2[#]))&,
				k1];
			h2OmegaGW2Sol//=WaitAll;
			DumpSave[FileNameJoin@{Directory[],InputString<>"_GW.mx"},{h2OmegaGW2Sol}];
		];
	,
		Code[
			Get[FileNameJoin@{Directory[],InputString<>"_GW.mx"}];
		];
	];
	Comment@"Plot.";
	Code[
		Expr=ListLogLogPlot[Transpose[{conv*k1,h2OmegaGW2Sol}],Joined->True,AxesLabel->{"f/Hz","h2OmegaGW2(k)"},PlotLabel->"Log-LogPlotofOmegaGW2(k)",PlotRange->{Automatic,Automatic}];
		Export[InputString<>"_GW.csv",Transpose[{conv*k1,h2OmegaGW2Sol}],"CSV"];
	];
	DisplayExpression@Expr;
];

ProcessSpectrum/@AllNames;

Supercomment@"This is the end of the GW integration script.";
Quit[];

(*===========================================================================================================================*)
(*  Benjamin, I've not examined any of the content below here, so if you want it to be displayed then feel free to clean it  *)
(*===========================================================================================================================*)

(*Set the working directory to the specified path*)
SetDirectory["C:\\Users\\Bubba\\Documents\\VSCode\\PBH\\PBH\\"];

(*Confirm the current working directory*)
Directory[]
dataToSave=Transpose[{conv*k1,h2OmegaGW2Sol}];
Export["GW_signature_data_1e8_lower.txt",dataToSave,"Table"]

SetDirectory[NotebookDirectory[]];
data=Import["plis_HLVK.dat","Table","HeaderLines"->14];

lowF=10;
highF=10^6;
frequency=10^data[[All,1]];
energydensity=10^data[[All,2]];

(*Plot*)
ListLogLogPlot[{Transpose[{conv*k1,h2OmegaGW2Sol}],Transpose[{frequency,energydensity}]},Joined->{True,True},PlotStyle->{Blue,Red},AxesLabel->{"f/Hz","h2OmegaGW2(k)"},FrameLabel->{"Frequency(Hz)","EnergyDensity"},PlotLabel->"Log-Log Plot of OmegaGW2(k)",PlotRange->{{lowF,highF},{10^-10,10^-5}},Filling->{2->1},FillingStyle->Orange,PlotLegends->{"h2OmegaGW2(k)","aLIGO Sensitivity"}];
