{
  gROOT->ProcessLine(".x rootlogon.C");
  

  TString outputdir = "../Histos/";

  //TString samples   = "../samples/samples_highMT2_4600_pnef.dat";
  TString samples   = "../samples/samples_Had.dat";
  //TString samples   = "../samples/samples_scale.dat";

  // TString samples   = "../samples/samples_highMT2_lept_4600.dat";

  TString MT2_REGIME = "MT2b";
  TString HT_LABEL = "950";
  TString LABEL = "test"; 
  TString FNAME = MT2_REGIME+"-HT_"+HT_LABEL+"-"+LABEL;


  int verbose =0;

  gSystem->CompileMacro("../MT2Code/src/AnalysisAndHistos.cc", "k");

  AnalysisAndHistos *tA = new AnalysisAndHistos(outputdir, "histos-"+FNAME+".root");

  tA->setVerbose(verbose);
  tA->init(samples);
  
  string HTCut = "";
  if(HT_LABEL=="750to950") HTCut = "misc.HT > 750 && misc.HT<=950 ";
  else if (HT_LABEL=="950") HTCut = "misc.HT>950 ";
  else if (HT_LABEL=="750") HTCut = "misc.HT>750 ";
  else {cout << "HT_LABEL not supported, exiting"; return;}

  std::ostringstream cutStream;
  cutStream <<      "  misc.Run < 179959 && " 
	    << " misc.MT2>125 && misc.MT2<9999 && "
	    << "misc.MET>=30"                            << "&&"
    	    <<  HTCut                          << "&&"
	    << "misc.Jet0Pass ==1"                       << "&&"
	    << "misc.Jet1Pass ==1"                       << "&&"
	    << "misc.SecondJPt  >100"                    << "&&"
	    << "misc.PassJetID ==1"                      << "&&"
	    << "misc.Vectorsumpt < 70 "                   << "&&"
	    <<" (NEles==0 || ele[0].lv.Pt()<10) && " 
	    << "(NMuons==0 || muo[0].lv.Pt()<10) && "
	    << "misc.HBHENoiseFlagIso == 0"                 << "&&"
	    << " misc.CSCTightHaloID==0 && " 
	    << "misc.CrazyHCAL==0";
  
  if(MT2_REGIME=="MT2"){
    cutStream  << " && misc.MinMetJetDPhi >= 0.3"
	       << " && NJetsIDLoose40 > 2 ";
  }
  else if(MT2_REGIME=="MT2b"){
    cutStream  << " && misc.MinMetJetDPhi4 >= 0.3"
	       << " && NBJets != 0"                        
	       << " && NJetsIDLoose40 > 3 ";
  }

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
		<< "|| (trigger.HLT_PFHT650_v1==1 && misc.Run>=179959)" 
		<< " )";
  TString trigger = triggerStream.str().c_str();
  
  TString cuts = cutStream.str().c_str();
  
  //cuts += " && jet[2].lv.Pt()>20";
	
  tA->Analysis(MT2_REGIME, HT_LABEL,3, 0, 0, trigger, cuts, false);


}
