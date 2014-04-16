{

  //this macro is basically the same as run_MassPlots_HT.C
  //this particular function is used to run over data triggered by the dileptonic trigger stream

  //run via root -l run_MassPlots_2l.C


  TString outputdir = "data_20fb-1_2l/";

  bool highHT   = false;
  bool leptonic = false;
  TString htcut = "";
//  htcut = "misc.HT<750";
  htcut = "misc.HT>450";
//  htcut = "misc.HT>450";

//  TString samples  = "samples/samples_MET.dat";
//   TString samples  = "samples/samples_2l_highHT.dat";
//  TString samples  = "samples/samples_emu.dat";
  TString samples  = "samples/samples_mumu.dat";
//  TString samples  = "samples/samples_mumu_allVSPT.dat";
//  TString samples  = "samples/samples_ee.dat";
//  TString samples  = "samples/samples_HT_filter_allVSPT.dat";


  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "k");
  gROOT->ProcessLine(".x SetStyle_PRD.C");

    const int gNMT2bins_2j0b_lHT                      = 9;
    double  gMT2bins_2j0b_lHT[gNMT2bins_2j0b_lHT+1]   = {190,200, 240, 290, 350, 420, 490, 570, 650, 750};
    const int gNMT2bins_2j1b_lHT                      = 7;
    double  gMT2bins_2j1b_lHT[gNMT2bins_2j1b_lHT+1]   = {190,200, 250, 310, 380, 450, 550, 700};
    const int gNMT2bins_3j0b_lHT                      = 9;
    double  gMT2bins_3j0b_lHT[gNMT2bins_3j0b_lHT+1]   = {190,200, 240, 290, 350, 420, 490, 570, 650, 750};
    const int gNMT2bins_3j1b_lHT                      = 7;
    double  gMT2bins_3j1b_lHT[gNMT2bins_3j1b_lHT+1]   = {190,200, 250, 310, 380, 460, 550, 700};
    const int gNMT2bins_3j2b_lHT                      = 5;
    double  gMT2bins_3j2b_lHT[gNMT2bins_3j2b_lHT+1]   = {190,200, 250, 325, 425, 550};
    const int gNMT2bins_6j0b_lHT                      = 4;
    double  gMT2bins_6j0b_lHT[gNMT2bins_6j0b_lHT+1]   = {190,200, 280, 380, 520};
    const int gNMT2bins_6j1b_lHT                      = 4;
    double  gMT2bins_6j1b_lHT[gNMT2bins_6j1b_lHT+1]   = {190,200, 250, 325, 450};
    const int gNMT2bins_6j2b_lHT                      = 4;
    double  gMT2bins_6j2b_lHT[gNMT2bins_6j2b_lHT+1]   = {190,200, 250, 300, 400};
    const int gNMT2bins_3b_lHT                        = 3;
    double  gMT2bins_3b_lHT  [gNMT2bins_3b_lHT+1]     = {190,200, 280, 400};
    // HT > 750 && HT < 1200
    const int gNMT2bins_2j0b_mHT                      = 10;
    double  gMT2bins_2j0b_mHT[gNMT2bins_2j0b_mHT+1]   = {95,125, 150, 180, 220, 270, 325, 425, 580, 780, 1000};
    const int gNMT2bins_2j1b_mHT                      = 6;
    double  gMT2bins_2j1b_mHT[gNMT2bins_2j1b_mHT+1]   = {95,100, 135, 170, 260, 450, 700};
    const int gNMT2bins_3j0b_mHT                      = 10;
    double  gMT2bins_3j0b_mHT[gNMT2bins_3j0b_mHT+1]   = {95,160, 185, 215, 250, 300, 370, 480, 640, 800, 1000};
    const int gNMT2bins_3j1b_mHT                      = 7;
    double  gMT2bins_3j1b_mHT[gNMT2bins_3j1b_mHT+1]   = {95,150, 175, 210, 270, 380, 600, 900};
    const int gNMT2bins_3j2b_mHT                      = 6;
    double  gMT2bins_3j2b_mHT[gNMT2bins_3j2b_mHT+1]   = {95,130, 160, 200, 270, 370, 500};
    const int gNMT2bins_6j0b_mHT                      = 6;
    double  gMT2bins_6j0b_mHT[gNMT2bins_6j0b_mHT+1]   = {95,160, 200, 250, 325, 425, 600};
    const int gNMT2bins_6j1b_mHT                      = 5;
    double  gMT2bins_6j1b_mHT[gNMT2bins_6j1b_mHT+1]   = {95,150, 190, 250, 350, 500};
    const int gNMT2bins_6j2b_mHT                      = 5;
    double  gMT2bins_6j2b_mHT[gNMT2bins_6j2b_mHT+1]   = {95,130, 170, 220, 300, 450};
    const int gNMT2bins_3b_mHT                        = 4;
    double  gMT2bins_3b_mHT[gNMT2bins_3b_mHT+1]       = {95,125, 175, 275, 450};

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
  tA->SetEventsPerGeV(true);
  tA->SetPileUpReweight(true);
  tA->SetbSFWeights(true);
  tA->SetSave(false);
  //tA->SetSave(true);

  std::ostringstream cutStream;
  cutStream 
      << " " 
      << "misc.HT>450"                             << "&&"
      << "misc.Jet0Pass ==1"                       << "&&"
      << "misc.Jet1Pass ==1"                       << "&&"
      << "misc.SecondJPt  > 100"                   << "&&"
      << "misc.MinMetJetDPhi4Pt40 >0.3"            << "&&"
      << "misc.PassJet40ID ==1"                    << "&&"
      << "misc.Vectorsumpt < 70";
  std::ostringstream basecutStream;
  basecutStream 
      << " " 
      << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"  << "&&" // for rare SM samples
      << "misc.CSCTightHaloIDFlag == 0"             << "&&"
      << "misc.trackingFailureFlag==0"              << "&&"
      << "misc.eeBadScFlag==0"                      << "&&"
      << "misc.EcalDeadCellTriggerPrimitiveFlag==0" << "&&"
      << "misc.TrackingManyStripClusFlag==0"             << "&&"
      << "misc.TrackingTooManyStripClusFlag==0"          << "&&"
      << "misc.TrackingLogErrorTooManyClustersFlag==0"   << "&&"
      << "misc.CrazyHCAL==0";
//      << "&&(type1pfmet[0].Pt()>30||type1pfmet[0].Pt()/misc.CaloMETRaw<=2.)";
  
  if (leptonic){
  //  cutStream << "&& NMuons==0 && NEles==1 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NMuons==1 && NEles==0 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NMuons==0 && NEles==0 && NTausIDLoose3Hits==1";
  }
  else
    cutStream << "&& NMuons==2 && NEles==0 && NTausIDLoose3Hits==0&&muo[0].lv.Pt()>20&&muo[1].lv.Pt()>10";
//    cutStream << "&& NMuons==0 && NEles==2 && NTausIDLoose3Hits==0&&ele[0].lv.Pt()>20&&ele[1].lv.Pt()>20&&ele[0].IDLoose&&ele[1].IDLoose";
//    cutStream << "&& (NMuons+NEles+NTausIDLoose3Hits)==2";
//      cutStream << " ";

  if(htcut!="")
    cutStream << "&&" << htcut;
  else if (highHT)
    cutStream << "&& misc.HT>950";
  else
    cutStream << "&& misc.HT>750 && misc.HT<=950";

//    cutStream << "&&ZllRecalculate()";

  std::ostringstream triggerStream;

  triggerStream << "( ( "
/*		<< "(trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
  		<< "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )"
		<< "||("
		<< "trigger.HLT_PFHT350_PFMET100_v3==1 || trigger.HLT_PFHT350_PFMET100_v4==1 || trigger.HLT_PFHT350_PFMET100_v5==1 || "
		<< "trigger.HLT_PFHT350_PFMET100_v6==1 || trigger.HLT_PFHT350_PFMET100_v7==1 || trigger.HLT_PFNoPUHT350_PFMET100_v1==1 || "
		<< "trigger.HLT_PFNoPUHT350_PFMET100_v3==1 || trigger.HLT_PFNoPUHT350_PFMET100_v4==1 ) ) )";//&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
*/
/*  triggerStream << "( ( "
		<< "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
		<< "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
*/
//  << " trigger.HLT_EMu&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0 ))";
  << " trigger.HLT_DiMuons))";//&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0 ))";
//  << " trigger.HLT_DiElectrons&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0 ))";

  TString trigger = triggerStream.str().c_str();
  
  TString cuts = cutStream.str().c_str();
  TString basecuts = basecutStream.str().c_str();

TString rcuts = "ZllRecalculate()&&GetDiLeptonPt(0,1,0,10,76,106)>20&&" + cuts;//ll
//TString rcuts = cuts + "&&ZllRecalculate()&&GetDiLeptonPt(0,0,0,10,76,106)>180";//emu
TString rrcuts = rcuts + "&&GetDiLeptonInvMass(0,1,0,10,1)>76&&GetDiLeptonInvMass(0,1,0,10,1)<106";

TString triggerE = "trigger.HLT_DiElectrons";
TString triggerM = "trigger.HLT_DiMuons";
TString triggerL = "trigger.HLT_EMu";
TString ecuts = "NEles==2&&ele[0].lv.Pt()>20&&ele[1].lv.Pt()>20&&ele[0].IDLoose&&ele[1].IDLoose&&"+rrcuts;
TString mcuts = "NMuons==2&&muo[0].lv.Pt()>20&&muo[1].lv.Pt()>20&&"+rrcuts;
TString lcuts = "NEles==1&&NMuons==1&&ele[0].lv.Pt()>20&&muo[0].lv.Pt()>20&&ele[0].IDLoose&&"+rrcuts;
    
//   tA->makePlot("GetDiLeptonPt(0,1,0,10,76,106)", ecuts, -2, -10, -10, triggerE , "Z p_{T} [GeV]" , 35, 125, 1000, false, true, true, true, true, true,  1 ,  false);
//   tA->makePlot("GetDiLeptonPt(0,1,0,10,76,106)", mcuts, -2, -10, -10, triggerM , "Z p_{T} [GeV]" , 35, 125, 1000, false, true, true, true, true, true,  1 ,  false);

}
