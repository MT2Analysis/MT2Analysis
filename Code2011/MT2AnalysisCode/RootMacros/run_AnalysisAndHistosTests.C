{
  gROOT->ProcessLine(".x rootlogon.C");
  

  TString outputdir = "../Histos/";

  //TString samples   = "../samples/samples_highMT2_4600_pnef.dat";
  TString samples   = "../samples/samples_Had.dat";
  //TString samples   = "../samples/samples_scale.dat";

  // TString samples   = "../samples/samples_highMT2_lept_4600.dat";

  TString MT2_REGIME = "MT2";
  TString HT_LABEL = "750to950";
  TString LABEL = "test"; 
  TString FNAME = MT2_REGIME+"-HT_"+HT_LABEL+"-"+LABEL;


  int verbose =0;

  gSystem->CompileMacro("../MT2Code/src/AnalysisAndHistos.cc", "k");

  //  AnalysisAndHistos *tA = new AnalysisAndHistos(outputdir, "histos-MT2_150_9999-Jet40-HT_750-4400pb-Jet380HoE1.0to1.4-Jet100Eta0.2.root");//-NoMHT_240to260.root");//-NoSSVHEMBJets.root");//JetPt100_Phi02or29");//JetVeto_Pt30Eta3_ChfGT085_NhfGT04_NefGT07.root");//NOJetPt180to240_Phi02or29.root");

  AnalysisAndHistos *tA = new AnalysisAndHistos(outputdir, "histos-"+FNAME+".root");

  tA->setVerbose(verbose);
  tA->init(samples);
  
  std::ostringstream cutStream;
  cutStream <<  " " 
    "  misc.Run < 179959 && " 
    //    "  misc.Run <= 167043 && " 
    //" ((misc.isData && misc.Run>=160405 && misc.Run <= 161176 ) || (!misc.isData && 1==1) )&& " 
    //" ((misc.isData && misc.Run>=161216 && misc.Run <=163261 ) || (!misc.isData && 1==1) )&& " 
    //" ((misc.isData && misc.Run>=163270 && misc.Run <=163869 ) || (!misc.isData && 1==1) )&& " 
    // " ((misc.isData && misc.Run>=165088 && misc.Run <=165633 ) || (!misc.isData && 1==1) )&& " 
    //" ((misc.isData && misc.Run>=167078 && misc.Run <=173198) || (!misc.isData && 1==1) )&& " 
    //" ((misc.isData && misc.Run>=173236 && misc.Run <=178380) || (!misc.isData && 1==1) )&& " 
    //						<< " Sum$(jet.lv.Pt()>100. && (abs(jet.lv.Phi())>2.9 || abs(jet.lv.Phi())<0.2  ))!=0 &&"
    //" !( GetMHT(1, 20, 2.4,false)<240 || GetMHT(1, 20, 2.4,false)>260) &&"
    
    //<< " MHT[0].Phi() < 2.0 &&" 
    << " misc.MT2>125 && misc.MT2<9999 && "
    
    //<< " Sum$( jet.lv.Pt()>100 && abs(jet.lv.Eta())<3.0   && (jet.ChHadFrac>=0.85 || jet.NeuHadFrac>=0.4 || jet.NeuEmFrac>=0.7) )==0 &&" 
    //<< "( (misc.MT2>150 && misc.MT2<180 ) || ( misc.MT2>220 && misc.MT2<9999 ) )&& "
    //	    << "(misc.ProcessID!=6 || (misc.Event!=814918 && misc.Event!=5500089 && misc.Event!=1934425)) && "
	    << "NJetsIDLoose40 > 3 && "
    //<< " Sum$( jet.lv.Pt()>100 && abs( jet.lv.Eta()) < 0.2  )!=0  &&" 
    //<< "NJetsIDLoose40 > 2 && "
	    << "misc.MET>=30"                            << "&&"
    //<< "misc.caloHT50_ID > 750 "                          << "&&"
	    << "misc.HT>950  "                          << "&&"
    //	    << "!(misc.MET>220 && misc.MET<260)  &&" 
    //<< "( (misc.HT - misc.caloHT50_ID)>0  && (misc.HT - misc.caloHT50_ID)<=40 ) "                          << "&&"
    //      << " Sum$(jet.lv.Pt()>380 && ((jet.ChHadFrac+jet.NeuHadFrac)/(jet.ChEmFrac+jet.NeuEmFrac))>=1.  && ((jet.ChHadFrac+jet.NeuHadFrac)/(jet.ChEmFrac+jet.NeuEmFrac))<=1.4 ) !=0 && "

	    << "misc.Jet0Pass ==1"                       << "&&"
	    << "misc.Jet1Pass ==1"                       << "&&"
    //<< " Sum$(jet.lv.Pt() >=180 && jet.lv.Pt() <=240 && jet.isPFIDLoose==1 ) ==0 "                    << "&&"
    //	    << " Sum$(jet.lv.Pt() >160 && jet.lv.Pt() <300 && jet.isPFIDLoose==1 && jet.lv.Phi()>-0.127 &&  jet.lv.Phi()<0.384 &&   jet.lv.Eta()>-0.363 &&  jet.lv.Eta()<-0.128 ) ==0 "                    << "&&"
    //	    << " Sum$(jet.lv.Pt() >100 && jet.lv.Pt() <300 && jet.isPFIDLoose==1 && ( (jet.ChHadFrac+jet.NeuHadFrac)/(jet.NeuEmFrac+jet.ChEmFrac)>5.2 )) ==0 "                    << "&&"
	    << "misc.SecondJPt  >100"                    << "&&"
    //<< "NBJets == 0"                               << "&&"
	    << "misc.PassJetID ==1"                      << "&&"
	    << "misc.Vectorsumpt < 70 "                   << "&&"
    //<< "misc.Vectorsumpt > 30 "                   << "&&"
	    << "misc.MinMetJetDPhi4 >= 0.3"                 << "&&"
    //<< " (misc.MinMetJetDPhi >0.3||misc.MinMetJetDPhiIndex>3) "                 << "&&"
	    <<" (NEles==0 || ele[0].lv.Pt()<10) && " 
	    << "(NMuons==0 || muo[0].lv.Pt()<10) && "
	    << "misc.HBHENoiseFlagIso == 0"                 << "&&"
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
		      << "|| (trigger.HLT_PFHT650_v1==1 && misc.Run>=179959)" 
		      << " )";
	TString trigger = triggerStream.str().c_str();
	
	TString cuts = cutStream.str().c_str();
		
	//cuts += " && jet[2].lv.Pt()>20";
	
	tA->Analysis(3, 0, 0, trigger, cuts, false);


}
