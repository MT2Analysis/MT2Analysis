run_ZInvEstFromW(float MT2_low, float MT2_high, TString LABEL, TString HTString){
  cout << LABEL << endl;

  gROOT->ProcessLine(".x rootlogon.C");
  gSystem->CompileMacro("../MT2Code/src/ZInvEstFromW.cc", "k");
  //gROOT->ProcessLine(".x SetStyle_PRD.C");


  TString outputdir = "../Histos/";
  
  string MT2_REGIME = "MT2"; //MT2 or MT2b
  int NJets = 3;
  bool EFF_CORR=false;
  //  float MT2_low = 275;
  //float MT2_high = 375;
  int NEles = 0;
  int NMuons = 1;
  
  bool IS_MC = false; //MC closure computation
  
  ostringstream MT2_low_s, MT2_high_s;
  MT2_low_s << MT2_low;
  MT2_high_s << MT2_high;
  //TString outname = "ZInv-"+MT2_REGIME+"_"+MT2_low_s.str()+"_"+MT2_high_s.str()+"-Jet40-HT_750-1Mu_pt10_4400pb_NoMETCorr_Zk1_Wk1_noBTag_NoPU_NoDPhi";//_NoVSPtNoDPhi";
  //TString outname = "ZInv-"+MT2_REGIME+"_"+MT2_low_s.str()+"_"+MT2_high_s.str()+"-Jet40-HT_750-1Mu_pt20_BJet40_4400pb_Zk1_Wk1_NoMT_NoDPhi_NoCorr_3D";//_NoVSPtNoDPhi";
  TString outname = "ZInv-"+MT2_REGIME+"_"+MT2_low_s.str()+"_"+MT2_high_s.str()+"-"+LABEL.Data();
//Jet40-HT_750to950-1Mu_pt10_4400pb_Zk1_Wk1_3D";//_NoVSPtNoDPhi";
  
  //TString outname = "test";
  if(IS_MC) outname += "-MC_Closure";
  
  TString samples   = "../samples/samples_highMT2_lept_3D_4600.dat";
  int verbose = 0;
  
  //gSystem->CompileMacro("../MT2Code/src/ZInvEstFromW.cc", "k");
  //gROOT->ProcessLine(".x SetStyle_PRD.C");
  
  ZInvEstFromW *tA = new ZInvEstFromW(outputdir);
  
  
  
  tA->setVerbose(verbose);
  tA->init(samples);
  
  string HTCut = "";
  if(HTString=="750to950") HTCut = "misc.HT > 750 && misc.HT<=950 ";
  else if (HTString=="950") HTCut = "misc.HT>950 ";
  else if (HTString=="750") HTCut = "misc.HT>750 ";
  else {cout << "HTString not supported, exiting"; return;}

  std::ostringstream cutStream;
  if(MT2_REGIME=="MT2"){
    cutStream << "   misc.Run < 179959 && " 
      //<< "Sum$(jet[].lv.Pt()>40)>2"                            << "&&"
      //<< "misc.MET>=240 && misc.MET<=320"                            << "&&"
      //<< " misc.MT2 > 100"                           << "&&"
      	      << "NJetsIDLoose40 > 2 && "
      //      << "misc.HT > 750 "                          << "&&"
      //<< "misc.HT > 950 "                          << "&&"
      //	      << "misc.HT > 750 && misc.HT<=950 "                          << "&&"
      //<< "misc.caloHT50_ID >750 "                  << "&&"
	      << HTCut << " && " 
	      << "misc.Jet0Pass ==1"                       << "&&"
	      << "misc.Jet1Pass ==1"                       << "&&"
	      << "misc.PassJetID ==1"                      << "&&"
      //	    << "misc.Vectorsumpt < 70"                   << "&&"
      //<< "misc.MinMetJetDPhi <=0.3"                 << "&&"
	      << "misc.HBHENoiseFlagIso == 0"                 << "&&"
      
	      << "misc.CrazyHCAL==0 && "
      //    << "misc.LeadingJPt >150"                    << "&&"
	      << " misc.CSCTightHaloID==0 && " 
	      << "misc.SecondJPt  >100";//                    << "&&"
    //<< "NBJets >0"                               << "&&"
    //	  << "misc.MT2 > 200 && misc.MT2<400"          << "&&"
    //	  << "misc.MT2 > 100 && misc.MT2<150"          << "&&"
    //	  << "misc.CrazyHCAL==0";
  }
  else{
    cutStream << "  misc.Run < 179959 && " 
      //<< "jet[2].lv.Pt()>40 && "
      //	      << "jet[3].lv.Pt()>40 && "
      	      << "NJetsIDLoose40 > 3 && "
	      << "misc.HT > 750 && misc.HT<=950"                          << "&&"
      //<< "misc.HT > 950 "                          << "&&"
	      << "misc.Jet0Pass ==1"                       << "&&"
	      << "misc.Jet1Pass ==1"                       << "&&"
      //<< "misc.PassJetID ==1"                      << "&&"
      //	    << "misc.Vectorsumpt < 70"                   << "&&"
      //    << "misc.MinMetJetDPhi >0.3"                 << "&&"
	      << "misc.HBHENoiseFlag == 0"                 << "&&"
	      << "misc.CrazyHCAL==0 && "
	      << "misc.LeadingJPt >150"                    << "&&"
	      << " misc.CSCTightHaloID==0 && " 
	      << "misc.SecondJPt  >100";// 
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
    		<< "|| (trigger.HLT_HT650_v4 ==1 && (misc.Run>=178420 && misc.Run< 179959))" << ""
    //<< "|| (trigger.HLT_PFHT650_v1==1 && misc.Run>=179959)" 
    //		<< " || (misc.Run<179959) )" ;
    		<< " )";
  TString trigger = triggerStream.str().c_str();
  
  TString cuts = cutStream.str().c_str();
  
  //actual call
  cout << "\n\n\n%TWISTY{mode=\"div\" showlink=\"MT2 " << MT2_low<<"-"<< MT2_high <<"\" hidelink=\"Hide MT2 " << MT2_low<<"-"<< MT2_high <<" \"}%" << endl;
  
  float *ttbarRes = new float[10];
  ttbarRes = tA->Analysis(MT2_REGIME, IS_MC, "histos_"+outname+"_1b.root", NJets, NEles, NMuons, MT2_low, MT2_high, trigger, cuts, true, EFF_CORR);
  tA->Analysis(MT2_REGIME, IS_MC, "histos_"+outname+"_0b.root", NJets, NEles, NMuons, MT2_low, MT2_high, trigger, cuts, false, EFF_CORR, ttbarRes);
  cout << "\n%ENDTWISTY%\n\n\n\n" << endl;
  
  
}
