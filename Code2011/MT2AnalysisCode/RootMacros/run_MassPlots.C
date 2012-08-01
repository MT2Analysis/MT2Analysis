{

  TString outputdir = "../MassPlots/";
  //TString samples   = "samples/samples_V02-04-01_highMT2_oldQCD_fastSim_allLM.dat";//
  //  TString samples   = "../samples/samples_highMT2_4600_pnef.dat";
//  TString samples   = "samples/samples_1479_lowMT2_oldQCD_fastsim.dat";
//	TString samples   = "samples/samples_skimmed.dat";
  TString samples   = "../samples/samples_Had.dat";
  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "k");
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  int gNMT2bins                   = 20;
  double  gMT2bins[gNMT2bins+1]   = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 135, 180, 260, 360, 500, 650}; 	
  int gNMT2Bbins                   = 15;
  double  gMT2Bbins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500}; 	

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
  tA->SetSave(false);
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
  
  std::ostringstream cutStream;
  cutStream << "  misc.Run < 179959 && " 
    //    << "(misc.ProcessID!=6 || (misc.Event!=814918 && misc.Event!=5500089 && misc.Event!=1934425))" << "&&"
	    << "NJetsIDLoose40 > 2 && "
	    << " ( ((jet.ChHadFrac+jet.NeuHadFrac)/(jet.ChEmFrac+jet.NeuEmFrac))>=1.1  && ((jet.ChHadFrac+jet.NeuHadFrac)/(jet.ChEmFrac+jet.NeuEmFrac))<=1.3 ) !=0 && "
    //<< "  !( GetMHT(1, 20, 2.4,false)<240 || GetMHT(1, 20, 2.4,false)>260) &&"
	    << "misc.MET>=30"                            << "&&"
	    << "misc.HT > 750 "                          << "&&"
    //  	  << "misc.HT > 750 && misc.HT<=950"                          << "&&"
	    << "misc.Jet0Pass ==1"                       << "&&"
	    << "misc.Jet1Pass ==1"                       << "&&"
    //<< "misc.LeadingJPt >150"                    << "&&"
    
	    << "misc.SecondJPt  >100"                    << "&&"
	  //<<"NJetsIDLoose>3 " << " && "
	  //  << "NBJets > 0"                               << "&&"
	    << "misc.PassJetID ==1"                      << "&&"
	    << "misc.Vectorsumpt < 70"                   << "&&"
    //  << " (misc.MinMetJetDPhi > 0.3 ||  misc.MinMetJetDPhiIndex>3)"                 << "&&"
	    << "misc.MinMetJetDPhi >0.3 "                 << "&&"
	    << "misc.HBHENoiseFlagIso == 0"                 << "&&"
	  //<< "misc.MT2 >= 180 && misc.MT2<=220 "                           << "&&"
	    <<" (NEles==0 || ele[0].lv.Pt()<10) && " 
	    << "(NMuons==0 || muo[0].lv.Pt()<10) && "
    //<< "misc.MT2 >=180 && misc.MT2<220"          << "&&"
	    << "misc.MT2 > 150 "          << "&&"
	    << " misc.CSCTightHaloID==0 && " 
	    << "misc.CrazyHCAL==0";
	
  std::ostringstream triggerStream;
  
  triggerStream << "( "
		<< "(trigger.HLT_HT440_v2 ==1 && misc.Run<161216)" << "||"
		<< "(trigger.HLT_HT450_v2 ==1 && (misc.Run>=161216 && misc.Run< 163269))" << "||"
		<< "(trigger.HLT_HT500_v3 ==1 && (misc.Run>=163269 && misc.Run<=163869))" << "||"
		<< "(trigger.HLT_HT500_v4 ==1 && (misc.Run>=165088 && misc.Run< 165970))" << "||"
		<< "(trigger.HLT_HT550_v5 ==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" << "||"
		<< "(trigger.HLT_HT550_v6 ==1 && (misc.Run==166346))" << "||"
		<< "(trigger.HLT_HT550_v7 ==1 && (misc.Run>=167078 && misc.Run< 170249))" << "||"
		<< "(trigger.HLT_HT550_v8 ==1 && (misc.Run>=170249 && misc.Run< 173236))" << "||"
		<< "(trigger.HLT_HT600_v1 ==1 && (misc.Run>=173236 && misc.Run< 178420))" << ""
    //		<< " || misc.Run< 178420 "
		<< "|| (trigger.HLT_HT650_v4 ==1 && (misc.Run>=178420 && misc.Run< 179959))" << ""
    //    		<< "|| (trigger.HLT_PFHT650_v1==1 && misc.Run>=179959)" 
		<< " || (misc.Run<179959) )"; 
  //			<< " )";
  TString cuts = cutStream.str().c_str();
  TString ele_cuts = cuts+"&&ele[0].lv.Pt()>10";

  //cuts += " && jet[2].lv.Pt()>20";


	std::ostringstream triggerStream;
	triggerStream << "( "
		<< "(trigger.HLT_HT440_v2 ==1 && misc.Run<161216)" << "||"
		<< "(trigger.HLT_HT450_v2 ==1 && (misc.Run>=161216 && misc.Run< 163269))" << "||"
		<< "(trigger.HLT_HT500_v3 ==1 && (misc.Run>=163269 && misc.Run<=163869))" << "||"
		<< "(trigger.HLT_HT500_v4 ==1 && (misc.Run>=165088 && misc.Run< 165970))" << "||"
		<< "(trigger.HLT_HT550_v5 ==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" << "||"
		      << "(trigger.HLT_HT550_v6 ==1 && (misc.Run==166346))" << "||"
		      << "(trigger.HLT_HT550_v7 ==1 && (misc.Run>=167078 && misc.Run< 170249))" << "||"
		      << "(trigger.HLT_HT550_v8 ==1 && (misc.Run>=170249 && misc.Run< 173236))" << "||"
		      << "(trigger.HLT_HT600_v1 ==1 && (misc.Run>=173236 && misc.Run< 178420))" << "||"
		      << "(trigger.HLT_HT650_v4 ==1 && (misc.Run>=178420 && misc.Run< 179959))" << ")";
	TString trigger = triggerStream.str().c_str();
	  

	  //       samples , variable,                     cuts,njet, nlep,   HLT,    xtitle             nbins  bins          flip_order,    log  , comp ,  ratio, stack, overlay
	  //         tA->makePlot("jet.bTagProbSSVHP",    cuts,     -4,  2,    trigger , "bTagProbSSVHP"        ,50,  -2, 7,          false,     true,  true,   true,  true,  true , 1);
	  //         tA->makePlot("jet.bTagProbTCHP",     cuts,    -4,  2,    trigger , "bTagProbTCHP"         ,25,  -10, 15,          false,     true,  true,   true,  true,  true , 1);
	  //          tA->makePlot("misc.LeadingJPt",       cuts,    -3, 0,    trigger , "leading JPt" , 100,  0, 1000,    false,     false,  true,   true,  true,  true , 1);
	  //           tA->makePlot("ele[0].lv.Pt()",       cuts,   -4, -11,    trigger , "ele Pt"      , 15,  0, 300,    false,     false,  true,   true,  true,  true , 1);
	  //           tA->makePlot("muo[0].lv.Pt()",      cuts,    -4, -13,    trigger , "muo Pt"      , 15,  0, 300,    false,     false,  true,   true,  true,  true , 1);
	  //         tA->makePlot("muo[0].MT",           cuts,    -4, -13,    trigger , "muo MT"      , 15,  0, 300,    false,     false,  true,   true,  true,  true , 1);
	  //         tA->makePlot("ele[0].MT",           cuts,    -4, -11,    trigger , "ele MT"      , 15,  0, 300,    false,     false,  true,   true,  true,  true , 1);
	  //           tA->makePlot("misc.MET",             cuts,    -3,   -13,    trigger , "MET"         ,15, 0, 500,      false,     true ,  true,   true,  true,  true, 1);
	  //         tA->makePlot("misc.caloHT50",        cuts,    -3,   0,    trigger , "M_{T2}"      , 50, 400, 2000,      false,     true ,  true,   true,  true,  true, 1);
	  //	tA->makePlot("GetSqrtS(0,false,1,20,2.4,1)",  cuts,    -3,  0 , trigger , "#sqrt{s_{min}}"  , 80, 0, 4000,      false,     false ,  true,   true,  true,  true, 1);
	  //         tA->makePlot("misc.MT2",             cuts,    -3,  0 ,    trigger , "M_{T2}"      , 8, 0, 500,      false,     true ,  true,   true,  true,  true, 1);
	  //          tA->makePlot("misc.MT2",             cuts,    -2,  -13 ,    trigger , "M_{T2}"      , 10, 0, 400,      false,     true ,  true,   true,  true,  true, 1);
	  
	//tA->makePlot("jet.lv.Eta():jet.lv.Phi()",       cuts,   -4, -11,    trigger , "ele Pt"      , 15,  0, 300,    false,     false,  true,   true,  true,  true , 1);
	
	//	tA->makePlot("(jet.ChHadFrac+jet.NeuHadFrac)/(jet.ChEmFrac+jet.NeuEmFrac):misc.MT2",       "jet.lv.Pt()>200 && "+cuts,    -3,   0,  trigger , "# jets Pt>50" , 20,0,400,            false,     false ,  true,   true,  true,  true, 1);


	//tA->MakePlot("misc.MT2",              cuts,      -3,   -10,    trigger , "M_{T2}"      ,gNMT2bins, gMT2bins,false,     true ,  true,   true,  true,  true, 1);  
	//tA->makePlot("misc.MT2",             cuts,    -3,  -10 ,    trigger , "M_{T2}"      , 70, 0, 700,      false,     true ,  true,   true,  true,  true, 1);
	tA->makePlot("jet.lv.Pt()",             cuts,    -3,  -10 ,    trigger , "jet pt"      , 25, 0, 500,      false,     true ,  true,   true,  true,  true, 1);


	//       tA->makePlot("misc.MinMetJetDPhi",   cuts,    -2,    0,    trigger , "minDPhi"      , 50, 0, 3.2,      false,     true ,  true,   true,  true,  true, 1);
//         tA->makePlot("Znunu.RecoOSee_mll",   cuts,    -2,   2,    trigger , "m_{ll}"    ,  10, 71, 111,               true,     true ,  true,   true,  true,  true, 1);
//  tA->makePlot("pileUp.NVertices"     , cuts,    -3,   -13,    trigger , "NVertices"    , 50, 0, 50,                  true,     true ,  true,   true,  true,  true, 1);
//         tA->makePlot("NJets",                cuts,    -2,   -11,    trigger , "NJets"    , 10, 0, 10,                  false,     true ,  true,   true,  true,  true, 1);
//         tA->makePlot("HemiMassTop()",          cuts,    -4,  0 ,    trigger , "Mass"       , 60, 0, 600,      false,     true ,  true,   true,  true,  true, 1);
//         tA->makePlot("LeptJetDR(13,1,1)",        cuts,    -4,   -13,  trigger , "minDR(muo, B-jet)" , 10, 0, 5,            false,     false ,  true,   true,  true,  true, 1);

    //tA->makePlot("Sum$(jet[].lv.Pt()>50)",       cuts,    -3,   0,  trigger , "# jets Pt>50" , 15,0,15,            false,     false ,  true,   true,  true,  true, 1);
  //tA->makePlot("Sum$(jet[].lv.Pt()>40)",       cuts,    -3,   0,  trigger , "# jets Pt>40" , 15,0,15,            false,     false ,  true,   true,  true,  true, 1);
  //tA->makePlot("Sum$(jet[].lv.Pt()>30)",       cuts,    -3,   0,  trigger , "# jets Pt>30" , 15,0,15,            false,     false ,  true,   true,  true,  true, 1);
  //tA->makePlot("Sum$(jet[].lv.Pt()>20)",       cuts,    -3,   0,  trigger , "# jets Pt>20" , 15,0,15,            false,     false ,  true,   true,  true,  true, 1);


	   //      tA->PrintCutFlow(-3, 0, trigger, cuts);

  //  tA->PrintCutFlowMT2vsHT(trigger, cuts);


//       tA->PrintCutFlow(-3, 0, "MHT_HT");
//  gSystem->Exit(0);


}
