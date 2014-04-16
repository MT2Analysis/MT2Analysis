{

  //this macro is basically the same as run_MassPlots_HT.C
  //this particular function is used to run over data triggered by the MET/HTMHT and(!) HT trigger stream

  //run via root -l run_MassPlots_METHT.C


  TString outputdir = "Filtered/data_20fb-1_METandHT/";

  bool highHT   = false;
  bool MT2b     = false;
  bool leptonic = false;
  TString htcut = "";
  htcut = "misc.HT>450";
  outputdir += htcut!="" ? "" : highHT ? "hHT/" : "lHT/";
  if(htcut=="misc.HT<750&&misc.HT>450") outputdir += "lHT/";
  if(htcut=="misc.HT>450") outputdir += "aHT/";

  TString samples  = "samples/samples_HTandMET_filter.dat";

  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "k");
  gROOT->ProcessLine(".x SetStyle_PRD.C");


    const int gNMT2bins_2j0b_lHT                      = 8;
    double  gMT2bins_2j0b_lHT[gNMT2bins_2j0b_lHT+1]   = {200, 240, 290, 350, 420, 490, 570, 650, 750};
    const int gNMT2bins_2j1b_lHT                      = 6;
    double  gMT2bins_2j1b_lHT[gNMT2bins_2j1b_lHT+1]   = {200, 250, 310, 380, 450, 550, 700};
    const int gNMT2bins_3j0b_lHT                      = 8;
    double  gMT2bins_3j0b_lHT[gNMT2bins_3j0b_lHT+1]   = {200, 240, 290, 350, 420, 490, 570, 650, 750};
    const int gNMT2bins_3j1b_lHT                      = 6;
    double  gMT2bins_3j1b_lHT[gNMT2bins_3j1b_lHT+1]   = {200, 250, 310, 380, 460, 550, 700};
    const int gNMT2bins_3j2b_lHT                      = 4;
    double  gMT2bins_3j2b_lHT[gNMT2bins_3j2b_lHT+1]   = {200, 250, 325, 425, 550};
    const int gNMT2bins_6j0b_lHT                      = 3;
    double  gMT2bins_6j0b_lHT[gNMT2bins_6j0b_lHT+1]   = {200, 280, 380, 520};
    const int gNMT2bins_6j1b_lHT                      = 3;
    double  gMT2bins_6j1b_lHT[gNMT2bins_6j1b_lHT+1]   = {200, 250, 325, 450};
    const int gNMT2bins_6j2b_lHT                      = 3;
    double  gMT2bins_6j2b_lHT[gNMT2bins_6j2b_lHT+1]   = {200, 250, 300, 400};
    const int gNMT2bins_3b_lHT                        = 2;
    double  gMT2bins_3b_lHT  [gNMT2bins_3b_lHT+1]     = {200, 280, 400};
    // HT > 1200
    const int gNMT2bins_2j0b_hHT                      = 6;
    double  gMT2bins_2j0b_hHT[gNMT2bins_2j0b_hHT+1]   = {120, 150, 200, 260, 350, 550, 900};
    const int gNMT2bins_2j1b_hHT                      = 2;
    double  gMT2bins_2j1b_hHT[gNMT2bins_2j1b_hHT+1]   = {100, 180, 350};
    const int gNMT2bins_3j0b_hHT                      = 7;
    double  gMT2bins_3j0b_hHT[gNMT2bins_3j0b_hHT+1]   = {160, 185, 220, 270, 350, 450, 650, 1000};
    const int gNMT2bins_3j1b_hHT                      = 4;
    double  gMT2bins_3j1b_hHT[gNMT2bins_3j1b_hHT+1]   = {150, 180, 230, 350, 550};
    const int gNMT2bins_3j2b_hHT                      = 2;
    double  gMT2bins_3j2b_hHT[gNMT2bins_3j2b_hHT+1]   = {130, 200, 350};
    const int gNMT2bins_6j0b_hHT                      = 3;
    double  gMT2bins_6j0b_hHT[gNMT2bins_6j0b_hHT+1]   = {160, 200, 300, 500};
    const int gNMT2bins_6j1b_hHT                      = 3;
    double  gMT2bins_6j1b_hHT[gNMT2bins_6j1b_hHT+1]   = {150, 200, 300, 500};
    const int gNMT2bins_6j2b_hHT                      = 2;
    double  gMT2bins_6j2b_hHT[gNMT2bins_6j2b_hHT+1]   = {130, 200, 350};
    const int gNMT2bins_3b_hHT                        = 1;
    double  gMT2bins_3b_hHT[gNMT2bins_3b_hHT+1]       = {125, 300};
    // HT > 750 && HT < 1200
    const int gNMT2bins_2j0b_mHT                      = 9;
    double  gMT2bins_2j0b_mHT[gNMT2bins_2j0b_mHT+1]   = {125, 150, 180, 220, 270, 325, 425, 580, 780, 1000};
    const int gNMT2bins_2j1b_mHT                      = 5;
    double  gMT2bins_2j1b_mHT[gNMT2bins_2j1b_mHT+1]   = {100, 135, 170, 260, 450, 700};
    const int gNMT2bins_3j0b_mHT                      = 9;
    double  gMT2bins_3j0b_mHT[gNMT2bins_3j0b_mHT+1]   = {160, 185, 215, 250, 300, 370, 480, 640, 800, 1000};
    const int gNMT2bins_3j1b_mHT                      = 6;
    double  gMT2bins_3j1b_mHT[gNMT2bins_3j1b_mHT+1]   = {150, 175, 210, 270, 380, 600, 900};
    const int gNMT2bins_3j2b_mHT                      = 5;
    double  gMT2bins_3j2b_mHT[gNMT2bins_3j2b_mHT+1]   = {130, 160, 200, 270, 370, 500};
    const int gNMT2bins_6j0b_mHT                      = 5;
    double  gMT2bins_6j0b_mHT[gNMT2bins_6j0b_mHT+1]   = {160, 200, 250, 325, 425, 600};
    const int gNMT2bins_6j1b_mHT                      = 4;
    double  gMT2bins_6j1b_mHT[gNMT2bins_6j1b_mHT+1]   = {150, 190, 250, 350, 500};
    const int gNMT2bins_6j2b_mHT                      = 4;
    double  gMT2bins_6j2b_mHT[gNMT2bins_6j2b_mHT+1]   = {130, 170, 220, 300, 450};
    const int gNMT2bins_3b_mHT                        = 3;
    double  gMT2bins_3b_mHT[gNMT2bins_3b_mHT+1]       = {125, 175, 275, 450};

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
  tA->SetEventsPerGeV(true);
  tA->SetPileUpReweight(true);
  tA->SetbSFWeights(true);
  //tA->SetSave(false);
  tA->SetSave(true);


  std::ostringstream cutStream;
  cutStream 
      << " " 
   //   << "misc.MT2>200"                            << "&&"
   //   << "misc.MET>200"                            << "&&"
      << "misc.HT>450"                             << "&&"
      << "((misc.MET>200&&misc.HT<750)||(misc.MET>30&&misc.HT>=750))"                            << "&&"

      << "misc.Jet0Pass ==1"                       << "&&"
      << "misc.Jet1Pass ==1"                       << "&&"
      << "misc.SecondJPt  > 100"                   << "&&"
      << "misc.MinMetJetDPhi4Pt40 >0.3"            << "&&"
      << "misc.Vectorsumpt < 70";
  std::ostringstream basecutStream;
  basecutStream 
      << " " 
      << "misc.PassJet40ID ==1"                    << "&&"
      << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"  << "&&" // for rare SM samples
      << "misc.CSCTightHaloIDFlag == 0"             << "&&"
      << "misc.trackingFailureFlag==0"              << "&&"
      << "misc.eeBadScFlag==0"                      << "&&"
      << "misc.EcalDeadCellTriggerPrimitiveFlag==0" << "&&"
      << "misc.TrackingManyStripClusFlag==0"             << "&&"
      << "misc.TrackingTooManyStripClusFlag==0"          << "&&"
      << "misc.TrackingLogErrorTooManyClustersFlag==0"   << "&&"
      << "misc.CrazyHCAL==0";
basecutStream << "&& misc.MET/misc.CaloMETRaw<=2.";
  
  if (leptonic){
  //  cutStream << "&& NMuons==0 && NEles==1 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NMuons==1 && NEles==0 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NMuons==0 && NEles==0 && NTausIDLoose3Hits==1";
  }
  else{
  //  cutStream << "&& NMuons==0 && NEles==0";
  //  cutStream << "&& NMuons==0 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NTausIDLoose3Hits==0 && NEles==0";
    cutStream << "&& NMuons==0 && NEles==0 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NMuons==0 && NEles==1 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NMuons==1 && NEles==0 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NMuons==0 && NEles==0 && NTausIDLoose3Hits==1";//&&NTaus==1";
  //  cutStream << "&& ((NMuons==1 && NEles==0 && NTausIDLoose3Hits==0)||(NMuons==0 && NEles==1 && NTausIDLoose3Hits==0))";
 //   cutStream << "&& (NMuons+NEles+NTausIDLoose3Hits)==1";
  //  cutStream << "&& ((NMuons==1&&muo[0].MT<100)||(NEles==1&&ele[0].MT<100)||(NTausIDLoose3Hits==1&&tau[0].MT<100))";
  }

  if(htcut!="")
    cutStream << "&&" << htcut;
  else if (highHT)
    cutStream << "&& misc.HT>950";
  else
    cutStream << "&& misc.HT>750 && misc.HT<=950";

  std::ostringstream triggerStream;
  triggerStream << "( ( ( "
		<< "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
		<< "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1)"
		<< "||( ( "
		<< "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
  		<< "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )"
		<< "||("
		<< "trigger.HLT_PFHT350_PFMET100_v3==1 || trigger.HLT_PFHT350_PFMET100_v4==1 || trigger.HLT_PFHT350_PFMET100_v5==1 || "
		<< "trigger.HLT_PFHT350_PFMET100_v6==1 || trigger.HLT_PFHT350_PFMET100_v7==1 || trigger.HLT_PFNoPUHT350_PFMET100_v1==1 || "
		<< "trigger.HLT_PFNoPUHT350_PFMET100_v3==1 || trigger.HLT_PFNoPUHT350_PFMET100_v4==1 ) ) ) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";// && misc.HT>9999";

  TString trigger = triggerStream.str().c_str();
  
  TString cuts = cutStream.str().c_str();
  TString basecuts = basecutStream.str().c_str();
  TString ele_cuts = "NEles==1&&ele[0].MT<100&& NMuons==0 && NTausIDLoose3Hits==0&&"+cuts;
  TString muo_cuts = "NMuons==1&&muo[0].MT<100&& NTausIDLoose3Hits==0 && NEles==0&&"+cuts;
  TString tau_cuts = "NTausIDLoose3Hits==1&&tau[0].MT<100&& NMuons==0 && NEles==0&&"+cuts;
  TString lep_cuts = "((NTausIDLoose3Hits==1&&tau[0].MT<100&& NMuons==0 && NEles==0)||(NMuons==1&&muo[0].MT<100&& NTausIDLoose3Hits==0 && NEles==0)||(NEles==1&&ele[0].MT<100&& NMuons==0 && NTausIDLoose3Hits==0))&&"+basecuts;
        
  TString higgscut = "NJetsIDLoose40>=4&&NBJetsCSVM>=2&&"+cuts;
    
//   samples , variable , cuts,njet, nlep, HLT     ,  xtitle        nbins   bins   flip_order,  log , comp, ratio,  stack,overlay,times,underflow
//tA->makeplot("misc.MT2", cuts, basecuts,   -2,-10, -10, trigger , "M_{T2} [GeV]" , 40,  200,  800,     false,  true,  true,  false,  true,   true,  1  ,   true,  true);

}
