{

  //this macro is basically the same as run_MassPlots_HT.C
  //this particular function is used to run over data triggered by the MET/HTMHT trigger stream

  //run via root -l run_MassPlots_MET.C


  TString outputdir = "Filtered/data_20fb-1_MET/";

  bool highHT   = false;
  bool leptonic = false;
  TString htcut = "";
  htcut = "misc.HT<750&&misc.HT>450";

  if(htcut=="misc.HT<750&&misc.HT>450") outputdir += "lHT/";
  if(htcut=="misc.HT>450") outputdir += "aHT/";

  TString samples  = "samples/samples_MET_filter.dat";

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
      << "misc.MT2>200"                            << "&&"
      << "misc.MET>200"                            << "&&"
      << "misc.HT>450"                             << "&&"
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
  else {
  //  cutStream << "&& NMuons==0 && NEles==0";
  //  cutStream << "&& NMuons==0 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NTausIDLoose3Hits==0 && NEles==0";
    cutStream << "&& NMuons==0 && NEles==0 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NMuons==0 && NEles==1 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NMuons==1 && NEles==0 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NMuons==0 && NEles==0 && NTausIDLoose3Hits==1";//&&NTaus==1";
  //  cutStream << "&& ((NMuons==1 && NEles==0 && NTausIDLoose3Hits==0)||(NMuons==0 && NEles==1 && NTausIDLoose3Hits==0))";
  //  cutStream << "&& (NMuons+NEles+NTausIDLoose3Hits)==1";
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
		<< "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
  		<< "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )"
		<< "||("
		<< "trigger.HLT_PFHT350_PFMET100_v3==1 || trigger.HLT_PFHT350_PFMET100_v4==1 || trigger.HLT_PFHT350_PFMET100_v5==1 || "
		<< "trigger.HLT_PFHT350_PFMET100_v6==1 || trigger.HLT_PFHT350_PFMET100_v7==1 || trigger.HLT_PFNoPUHT350_PFMET100_v1==1 || "
		<< "trigger.HLT_PFNoPUHT350_PFMET100_v3==1 || trigger.HLT_PFNoPUHT350_PFMET100_v4==1 ) ) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";// && misc.HT>9999";

  TString trigger = triggerStream.str().c_str();
  
  TString cuts = cutStream.str().c_str();
  TString basecuts = basecutStream.str().c_str();
  TString ele_cuts = "NEles==1&&ele[0].MT<100&& NMuons==0 && NTausIDLoose3Hits==0&&"+cuts;
  TString muo_cuts = "NMuons==1&&muo[0].MT<100&& NEles==0 && NTausIDLoose3Hits==0&&"+cuts;
  TString tau_cuts = "NTausIDLoose3Hits==1&&tau[0].MT<100&& NMuons==0 && NEles==0&&"+cuts;
  TString nol_cuts = "NTausIDLoose3Hits==0&& NMuons==0 && NEles==0&&"+cuts;
        
  TString higgscut = "NJetsIDLoose40>=4&&NBJetsCSVM>=2&&"+cuts;
    

//   samples , variable , cuts,njet, nlep, HLT     ,  xtitle        nbins   bins   flip_order,  log , comp, ratio,  stack,overlay,times,underflow
//tA->makeplot("misc.MT2", cuts, basecuts,   2,  -10, -10, trigger , "M_{T2} [GeV]" , 40,  200,  800,     false,  true,  true,  true,  true,   true,  1  ,   true);
//tA->makeplot("misc.MT2", cuts, basecuts,   35,  -10, -10, trigger , "M_{T2} [GeV]" , 40,  200,  800,     false,  true,  true,  true,  true,   true,  1  ,   true);
//tA->makeplot("misc.MT2", cuts, basecuts,   -6,  -10, -10, trigger , "M_{T2} [GeV]" , 40,  200,  800,     false,  true,  true,  true,  true,   true,  1  ,   true);
//tA->makeplot("misc.MT2", cuts, basecuts,   -2,  0, -10, trigger , "M_{T2} [GeV]" , 40,  200,  800,     false,  true,  true,  true,  true,   true,  1  ,   true);
//tA->makeplot("misc.MT2", cuts, basecuts,   -2,  1, -10, trigger , "M_{T2} [GeV]" , 40,  200,  800,     false,  true,  true,  true,  true,   true,  1  ,   true);
//tA->makeplot("misc.MT2", cuts, basecuts,   -2,  2, -10, trigger , "M_{T2} [GeV]" , 40,  200,  800,     false,  true,  true,  true,  true,   true,  1  ,   true);
//tA->makeplot("misc.MT2", cuts, basecuts,   -2,  -3, -10, trigger , "M_{T2} [GeV]" , 40,  200,  800,     false,  true,  true,  true,  true,   true,  1  ,   true);
}
