{

  //this macro is basically the same as run_MassPlots_HT.C
  //this particular function is used to run over data triggered by the MET trigger stream and it runs over MT2 trees of tag V02-02-00(!!!)

  //run via root -l run_MassPlots_MET_V02-02-00.C

  TString outputdir = "data_11fb-1_MET/";

  bool highHT   = false;
  bool leptonic = false;
  TString htcut = "";
  htcut = "misc.HT>=450&&misc.HT<750";

  outputdir += leptonic ? "lep/" : "had/";
  outputdir += htcut!="" ? "" : highHT ? "hHT/" : "lHT/";

  TString samples  = "samples/samples_type1MET_CHS_53X_MET_Run2012ABC.dat";
  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "k");
  gROOT->ProcessLine(".x SetStyle_PRD.C");

   //note that this is not the current MT2 binning of 8 TeV analysis
   int gNMT2bins_2j0b                  = 12;
   double  gMT2bins_2j0b[gNMT2bins_2j0b+1]   = {100, 160, 200, 220, 240, 260, 280, 300, 330, 370, 425, 600, 900};
   int gNMT2bins_2j1b                  = 9;
   double  gMT2bins_2j1b[gNMT2bins_2j1b+1]   = {100, 160, 200, 220, 240, 260, 290, 340, 440, 600};
   int gNMT2bins_3j0b                  = 9;
   double  gMT2bins_3j0b[gNMT2bins_3j0b+1]   = {100, 160, 200, 225, 260, 300, 360, 450, 560, 700};
   int gNMT2bins_3j1b                  = 7;
   double  gMT2bins_3j1b[gNMT2bins_3j1b+1]   = {100, 160, 200, 230, 275, 340, 450, 650};
   int gNMT2bins_3j2b                  = 5;
   double  gMT2bins_3j2b[gNMT2bins_3j2b+1]   = {100, 160, 200, 240, 320, 450};
   int gNMT2bins_6j0b                  = 6;
   double  gMT2bins_6j0b[gNMT2bins_6j0b+1]   = {100, 160, 200, 240, 280, 350, 500};
   int gNMT2bins_6j1b                  = 5;
   double  gMT2bins_6j1b[gNMT2bins_6j1b+1]   = {100, 160, 200, 230, 280, 500};
   int gNMT2bins_6j2b                  = 3;
   double  gMT2bins_6j2b[gNMT2bins_6j2b+1]   = {100, 200, 260, 400};
   int gNMT2bins_3b                  = 3;
   double  gMT2bins_3b[gNMT2bins_3b+1]   = {100, 200, 260, 400};


  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
  tA->SetEventsPerGeV(true);
  tA->SetPileUpReweight(true);
  tA->SetbSFWeights(true);
  tA->SetSave(true);
  //tA->SetSave(true);



  std::ostringstream cutStream;
  std::ostringstream cutStreamBase;
  cutStream << " " 
    //<< "misc.MT2>125"                           << "&&"
    << "misc.Jet0Pass ==1"                      << "&&"
    << "misc.Jet1Pass ==1"                      << "&&"
    << "misc.Vectorsumpt < 70"                  << "&&"
    << "misc.MET>200";
    cutStream << "&& misc.MinMetJetDPhi4 >0.3";
  }
  if (leptonic){
    //cutStream << "&&(NMuons>0||NEles >0)";
    //cutStream << "&& NMuons==1 && NEles == 0";
    //cutStream << "&& NMuons==0 && NEles == 1";
    //cutStream << "&& NMuons >0";
    //cutStream << "&& NEles >0";
    cutStream << "&& NMuons+NEles == 1";
  }
  else
    cutStream << "&& NMuons==0 && NEles==0";

  if(htcut!="")
    cutStream << "&&" << htcut;
  else if (highHT)
    cutStream << "&& misc.HT>950";
  else
    cutStream << "&& misc.HT>750 && misc.HT<=950";

  cutStreamBase << " " 
    << "misc.PassJetID ==1"                      << "&&"
    << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"  << "&&" // for rare SM samples
    << "misc.CSCTightHaloIDFlag == 0"                     << "&&"
    << "(misc.hcalLaserEventFlag==0|| misc.ProcessID==10)"<< "&&"
    << "misc.trackingFailureFlag==0"                      << "&&"
    << "misc.eeBadScFlag==0"                              << "&&"
    << "misc.EcalDeadCellTriggerPrimitiveFlag==0"         << "&&"
    << "misc.CrazyHCAL==0";

  std::ostringstream triggerStream;
  triggerStream << "( "
		<< "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
  		<< "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )";
  //triggerStream << "( HLT_QuadJet60_DiJet20_v1||trigger.HLT_QuadJet60_DiJet20_v2 )";
  TString trigger = triggerStream.str().c_str();
  
  TString cuts = cutStream.str().c_str();
  TString basecuts = cutStreamBase.str().c_str();

    
//    samples , variable , cuts, basecuts, njet, nlep, HLT     ,  xtitle        nbins   bins   flip_order,  log , comp, ratio,  stack,overlay,times,underflow
  //tA->Makeplot("misc.MT2", cuts, basecuts,   -3,  -10, trigger , "M_{T2} [GeV]" ,gNMT2bins,gMT2bins,false, false,  true,  true,  true,   true,  1  ,   false);
  //tA->makeplot("misc.MT2", cuts, basecuts,    -2, -10, -10, trigger , "M_{T2} [GeV]" , 50,  0,  900,     false,  true,  true,  false,  true,   true,  1  ,   false);
 
}
