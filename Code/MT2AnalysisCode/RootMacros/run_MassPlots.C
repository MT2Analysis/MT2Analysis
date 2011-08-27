{
  TString outputdir = "../MassPlots/";
  //TString samples   = "samples/samples_V02-04-01_highMT2_oldQCD_fastSim_allLM.dat";
  TString samples   = "../MT2Code/samples/samples_1479_highMT2_oldQCD_fastsim.dat";
//  TString samples   = "samples/samples_1479_lowMT2_oldQCD_fastsim.dat";
//	TString samples   = "samples/samples_skimmed.dat";
  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "k");
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  int gNMT2bins                   = 19;
  double  gMT2bins[gNMT2bins+1]   = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 135, 180, 260, 360, 650}; 	

  int gNMT2bins_l                   = 14;
  double  gMT2bins_l[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500}; 	

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
  tA->setVerbose(verbose);
  tA->init(samples);
  
	std::ostringstream cutStream;
	cutStream << " " 
          << "misc.MET>=30"                            << "&&"
	  << "misc.HT > 300 "                          << "&&"
	  << "misc.caloHT50_ID >600 "                  << "&&"
	  << "misc.Jet0Pass ==1"                       << "&&"
	  << "misc.Jet1Pass ==1"                       << "&&"
//	  << "misc.LeadingJPt >150"                    << "&&"
	  << "misc.SecondJPt  >100"                    << "&&"
//	  << "NBJets >0"                               << "&&"
	  << "misc.PassJetID ==1"                      << "&&"
	  << "misc.Vectorsumpt < 70"                   << "&&"
	  << "misc.MinMetJetDPhi >0.3"                 << "&&"
	  << "misc.HBHENoiseFlag == 1"                 << "&&"
	  << "misc.MT2 > 400"                           << "&&"
//	  << "misc.MT2 > 200 && misc.MT2<400"          << "&&"
//	  << "misc.MT2 > 100 && misc.MT2<150"          << "&&"
	  << "misc.CrazyHCAL==0";

	std::ostringstream triggerStream;
	triggerStream << "( "
	<< "(trigger.HLT_HT440_v2==1 && misc.Run<=161176)" << "||"
	<< "(trigger.HLT_HT450_v2==1 && (misc.Run>=161217 && misc.Run<=163261))" << "||"
	<< "(trigger.HLT_HT500_v3==1 && (misc.Run>=163270 && misc.Run<=163869))" << "||"
	<< "(trigger.HLT_HT550_v4==1 && (misc.Run>=165088 && misc.Run<=165633))" << "||"
	<< "(trigger.HLT_HT550_v5==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" << "||"
	<< "(trigger.HLT_HT550_v6==1 && (misc.Run==166346))"                     << "||"
	<< "(trigger.HLT_HT550_v7==1 && (misc.Run>=167078 && misc.Run<=167913))" << "||"
	<< "(trigger.HLT_HT550_v8==1 && (misc.Run>=169561 && misc.Run<=172619))" << " )";
	TString trigger = triggerStream.str().c_str();

  TString cuts = cutStream.str().c_str();
  TString ele_cuts = cuts+"&&ele[0].lv.Pt()>10";
    
//       samples , variable,                     cuts,njet, nlep,   HLT,    xtitle             nbins  bins          flip_order,    log  , comp ,  ratio, stack, overlay
//         tA->makePlot("jet.bTagProbSSVHP",    cuts,     -4,  2,    trigger , "bTagProbSSVHP"        ,50,  -2, 7,          false,     true,  true,   true,  true,  true , 1);
//         tA->makePlot("jet.bTagProbTCHP",     cuts,    -4,  2,    trigger , "bTagProbTCHP"         ,25,  -10, 15,          false,     true,  true,   true,  true,  true , 1);
//          tA->makePlot("misc.LeadingJPt",       cuts,    -3, 0,    trigger , "leading JPt" , 100,  0, 1000,    false,     false,  true,   true,  true,  true , 1);
//           tA->makePlot("ele[0].lv.Pt()",       cuts,   -4, -11,    trigger , "ele Pt"      , 15,  0, 300,    false,     false,  true,   true,  true,  true , 1);
//           tA->makePlot("muo[0].lv.Pt()",      cuts,    -4, -13,    trigger , "muo Pt"      , 15,  0, 300,    false,     false,  true,   true,  true,  true , 1);
//         tA->makePlot("muo[0].MT",           cuts,    -4, -13,    trigger , "muo MT"      , 15,  0, 300,    false,     false,  true,   true,  true,  true , 1);
//         tA->makePlot("ele[0].MT",           cuts,    -4, -11,    trigger , "ele MT"      , 15,  0, 300,    false,     false,  true,   true,  true,  true , 1);
//           tA->makePlot("misc.MET",             cuts,    -3,   0,    trigger , "MET"         ,15, 0, 500,      false,     true ,  true,   true,  true,  true, 1);
//         tA->makePlot("misc.caloHT50",        cuts,    -3,   0,    trigger , "M_{T2}"      , 50, 400, 2000,      false,     true ,  true,   true,  true,  true, 1);
//	tA->makePlot("GetSqrtS(0,false,1,20,2.4,1)",  cuts,    -3,  0 , trigger , "#sqrt{s_{min}}"  , 80, 0, 4000,      false,     false ,  true,   true,  true,  true, 1);
//         tA->makePlot("misc.MT2",             cuts,    -3,  0 ,    trigger , "M_{T2}"      , 8, 0, 500,      false,     true ,  true,   true,  true,  true, 1);
//          tA->makePlot("misc.MT2",             cuts,    -2,  -13 ,    trigger , "M_{T2}"      , 10, 0, 400,      false,     true ,  true,   true,  true,  true, 1);
           tA->makePlot("misc.MT2",             cuts,    -3,  0 ,    trigger , "M_{T2}"      , 70, 0, 700,      false,     true ,  true,   true,  true,  true, 1);
//          tA->makePlot("misc.MT2",             cuts,    -4,  0 ,    trigger , "M_{T2}"      , 32, 0, 480,      false,     true ,  true,   true,  true,  true, 1);
//          tA->makePlot("misc.HT",             cuts,    -3,  -1 ,    trigger , "HT"     ,       80, 0, 1000,      false,     true ,  true,   true,  true,  true, 1);
//         tA->makePlot("misc.MT2",               ele_cuts,    -3,  -11,   trigger , "M_{T2}"      , 20, 0, 500,     false,     false ,  true,   true,  true,  true, 1);
//       tA->makePlot("misc.MinMetJetDPhi",   cuts,    -2,    0,    trigger , "minDPhi"      , 50, 0, 3.2,      false,     true ,  true,   true,  true,  true, 1);
//         tA->makePlot("Znunu.RecoOSee_mll",   cuts,    -2,   2,    trigger , "m_{ll}"    ,  10, 71, 111,               true,     true ,  true,   true,  true,  true, 1);
//         tA->makePlot("pileUp.NVertices"     , cuts,    -3,   0,    trigger , "NVertices"    , 20, 0, 20,                  true,     true ,  true,   true,  true,  true, 1);
//         tA->makePlot("NJets",                cuts,    -2,   -11,    trigger , "NJets"    , 10, 0, 10,                  false,     true ,  true,   true,  true,  true, 1);
//           tA->MakePlot("misc.MT2",              cuts,      -3,   0,    trigger , "M_{T2}"      ,gNMT2bins, gMT2bins,false,     true ,  true,   true,  true,  true, 1);
//         tA->makePlot("HemiMassTop()",          cuts,    -4,  0 ,    trigger , "Mass"       , 60, 0, 600,      false,     true ,  true,   true,  true,  true, 1);
//         tA->makePlot("LeptJetDR(13,1,1)",        cuts,    -4,   -13,  trigger , "minDR(muo, B-jet)" , 10, 0, 5,            false,     false ,  true,   true,  true,  true, 1);
//       tA->PrintCutFlow(-3, 0, "HT");
//       tA->PrintCutFlow(-3, 0, "MHT_HT");

}
