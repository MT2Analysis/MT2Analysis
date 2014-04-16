{
  //this macro calls one of the make plot functions of the MassPlotter.cc
  //thus it makes one plot of chosen variable and selection per line
  //this is our standard way to do data/MC plots (in one dimension)
  //note that the way plot is produced is not fully 'TDR style'
  //but that has to be worked on in MassPlotter.cc

  //this particular function is used to run over data triggered by the HT trigger stream

  //run via root -l run_MassPlots_HT.C


  TString outputdir = "Filtered/data_20fb-1_HT/";

  bool highHT   = false;//I don't use that flag, but it would discriminate between  medium HT and high HT
  bool leptonic = false;//I don't use that flag, but it would discriminate between hadronic and leptonic region
  TString htcut = "";
  htcut = "misc.HT>750";//HT cut you want to apply
//  htcut = "misc.HT>750&&misc.HT<1200";
//  htcut = "misc.HT>1200";
//  outputdir += leptonic ? "lep/" : "had/";
  if(htcut=="misc.HT>750")               outputdir += "mhHT/";
  if(htcut=="misc.HT>750&&misc.HT<1200") outputdir += "mHT/";
  if(htcut=="misc.HT>1200")              outputdir += "hHT/";

//  outputdir += "jets/;//add what you need to discriminate the event

  TString samples  = "samples/samples_HT_filter.dat";//sample you want to read in

  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "k");//call the MassPlotter.cc
  gROOT->ProcessLine(".x SetStyle_PRD.C");

     //MT2 binning
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

  //define the MassPlotter
  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
  tA->setVerbose(verbose);
  tA->init(samples);//load the samples
  tA->SetIsPhoton(false);//if run over photon selection set this to true, usually it is false
  tA->SetEventsPerGeV(true);//for variables that are measured in GeV (like MT2/HT/MET/...) set this to true, for others (like angles/NJets/NBJets/...) set it to false
  tA->SetPileUpReweight(true);//do PU reweighting, default is true
  tA->SetbSFWeights(true);//do BTV SF reweighting, default is true
  //tA->SetQCDfromData(true);//this flag is used only by Bruno for some special case to show how good the QCD estimation method works compared to pure simulation
  //tA->SetSave(true);//if you want to save the plot
  tA->SetSave(false);//if you don't want to save the plot


  std::ostringstream cutStream;
  cutStream 
      << " " 
    //<< "NJetsIDLoose40>2"                       << "&&"
  //    << "misc.MT2>125"                            << "&&"
  //    << "misc.MT2>200"                            << "&&"
      << "misc.MET>30"                             << "&&"
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
  else
  //  cutStream << "&& NMuons==0 && NEles==0";
  //  cutStream << "&& NMuons==0 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NTausIDLoose3Hits==0 && NEles==0";
    cutStream << "&& NMuons==0 && NEles==0 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NMuons==0 && NEles==1 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NMuons==1 && NEles==0 && NTausIDLoose3Hits==0";
  //  cutStream << "&& NMuons==0 && NEles==0 && NTausIDLoose3Hits==1";//&&tau[0].isLooseID3Hits";
   // cutStream << "&& ((NMuons==1 && NEles==0 && NTausIDLoose3Hits==0)||(NMuons==0 && NEles==1 && NTausIDLoose3Hits==0))";
  //  cutStream << "&& (NMuons+NEles+NTausIDLoose3Hits)==1";

  if(htcut!="")
    cutStream << "&&" << htcut;
  else if (highHT)
    cutStream << "&& misc.HT>950";
  else
    cutStream << "&& misc.HT>750 && misc.HT<=950";

//photons
//   cutStream << "&&NPhotons==1&&photon[0].isLooseIso";

//cutStream << "&&NGenLepts>0&&Sum$(abs(genlept.ID)==16&&abs(genlept.MID)==15&&(abs(genlept.GMID)==23||abs(genlept.GMID)==22))";

  std::ostringstream triggerStream;
  triggerStream << "( ( "
		<< "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
		<< "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1))";// &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";

  TString trigger = triggerStream.str().c_str();
  
  TString cuts = cutStream.str().c_str();
  TString basecuts = basecutStream.str().c_str();
  TString ele_cuts = "NEles==1&&"+cuts;
  TString muo_cuts = "NMuons==1&&"+cuts;
  TString tau_cuts = "NTausIDLoose3Hits==1&&"+cuts;
  TString elemt_cuts = "NEles==1&&ele[0].MT<100&&"+cuts;
  TString muomt_cuts = "NMuons==1&&muo[0].MT<100&&"+cuts;
  TString taumt_cuts = "NTausIDLoose3Hits==1&&tau[0].MT<100&&"+cuts;
  TString phot_cuts = cuts+"&&NPhotons==1&&photon[0].isLooseID";
  TString med_cuts = cuts+"&&misc.HT<1200";
  TString hig_cuts = cuts+"&&misc.HT>=1200";
  TString mt2_cuts = "misc.MT2>100&&"+cuts;

  TString higgscut = "NJetsIDLoose40>=4&&NBJetsCSVM>=2&&"+cuts;

    
//the line below defines all the flags
//   samples , variable , cuts,njet,nbjet, nlep, HLT     ,  xtitle        nbins   bins   flip_order,  log , comp, ratio,  stack,overlay,times,underflow, save
//variable is the variable to plot
//cuts is the event selection - in some functions cuts is separated into cuts and basecuts (that separation is only important for saving issues as too long file names cannot be handled)
//njet is the selection on NJets: positive means ==NJets, negative means >=NJets, positive and 2 digits means for XY: X <=NJets<= Y
//nbjet means selection on NBJets, like NJets but special: -10 means no selection on NBJets (i.e. NBJets>=0)
//nlep means selection on leptons(e+mu), special next to NBJets: -11: 1 e, -13: 1 mu, -15: 1 tau and no e+mu (only place of tau in nlep)
//HLT means trigger selection
//xtitle means the title on x-axis of plot
//nbins bins - or - nbins binlow binup - means the binning of the plot
//flip_order means stack the mc backgrounds in opposite order to the standard way, default = false
//log means that y-axis of plot is log-style
//comp means composite (i.e. plot samples of same type like QCD as one histogram and not each as a separate histogram) - default = true
//ratio means make also plot with a MCsum/data ratio, for log plot default is true
//stack means stacking MC samples on top of each other - default is always true
//overlay means that susy signal is overlayed and not stacked on top of simulation - default is true
//times is the scale factor of susy signal cross section to nominal cross section - default is 1
//underflow means that the underflow is added to lowest bin - default depends on variable you plot
//save means on top of making the plot storing the most important histograms into a root file
tA->makeplot("misc.Vectorsumpt", mt2_cuts, basecuts,   -2,0, -1, trigger , "|MHT-MET| [GeV]" , 20,  0,  200,     false,  true,  true,  true,  true,   true,  1  ,   true,  false);
tA->makeplot("misc.Vectorsumpt", mt2_cuts, basecuts,   -2,-1, -1, trigger , "|MHT-MET| [GeV]" , 20,  0,  200,     false,  true,  true,  true,  true,   true,  1  ,   true,  false);

/*
     tA->Makeplot("misc.MT2", med_cuts, basecuts,   2, 0, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_2j0b_mHT, gMT2bins_2j0b_mHT, false,  false,  true,  false,  true,   true,  1  ,   false);
     tA->Makeplot("misc.MT2", med_cuts, basecuts,   2,-1, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_2j1b_mHT, gMT2bins_2j1b_mHT, false,  false,  true,  false,  true,   true,  1  ,   false);
     tA->Makeplot("misc.MT2", med_cuts, basecuts,  35, 0, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_3j0b_mHT, gMT2bins_3j0b_mHT, false,  false,  true,  false,  true,   true,  1  ,   false);
     tA->Makeplot("misc.MT2", med_cuts, basecuts,  35, 1, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_3j1b_mHT, gMT2bins_3j1b_mHT, false,  false,  true,  false,  true,   true,  1  ,   false);
     tA->Makeplot("misc.MT2", med_cuts, basecuts,  35, 2, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_3j2b_mHT, gMT2bins_3j2b_mHT, false,  false,  true,  false,  true,   true,  1  ,   false);
     tA->Makeplot("misc.MT2", med_cuts, basecuts,  -6, 0, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_6j0b_mHT, gMT2bins_6j0b_mHT, false,  false,  true,  false,  true,   true,  1  ,   false);
     tA->Makeplot("misc.MT2", med_cuts, basecuts,  -6, 1, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_6j1b_mHT, gMT2bins_6j1b_mHT, false,  false,  true,  false,  true,   true,  1  ,   false);
     tA->Makeplot("misc.MT2", med_cuts, basecuts,  -6, 2, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_6j2b_mHT, gMT2bins_6j2b_mHT, false,  false,  true,  false,  true,   true,  1  ,   false);
     tA->Makeplot("misc.MT2", med_cuts, basecuts,  -3,-3, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_3b_mHT,   gMT2bins_3b_mHT,   false,  false,  true,  false,  true,   true,  1  ,   false);


    tA->Makeplot("misc.MT2", hig_cuts,  basecuts,  2, 0, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_2j0b_hHT, gMT2bins_2j0b_hHT, false,  false,  true,  false,  true,   true,  1  ,   false);
    tA->Makeplot("misc.MT2", hig_cuts, basecuts,   2,-1, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_2j1b_hHT, gMT2bins_2j1b_hHT, false,  false,  true,  false,  true,   true,  1  ,   false);
    tA->Makeplot("misc.MT2", hig_cuts, basecuts,  35, 0, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_3j0b_hHT, gMT2bins_3j0b_hHT, false,  false,  true,  false,  true,   true,  1  ,   false);
    tA->Makeplot("misc.MT2", hig_cuts, basecuts,  35, 1, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_3j1b_hHT, gMT2bins_3j1b_hHT, false,  false,  true,  false,  true,   true,  1  ,   false);
    tA->Makeplot("misc.MT2", hig_cuts, basecuts,  35, 2, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_3j2b_hHT, gMT2bins_3j2b_hHT, false,  false,  true,  false,  true,   true,  1  ,   false);
    tA->Makeplot("misc.MT2", hig_cuts, basecuts,  -6, 0, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_6j0b_hHT, gMT2bins_6j0b_hHT, false,  false,  true,  false,  true,   true,  1  ,   false);
    tA->Makeplot("misc.MT2", hig_cuts, basecuts,  -6, 1, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_6j1b_hHT, gMT2bins_6j1b_hHT, false,  false,  true,  false,  true,   true,  1  ,   false);
    tA->Makeplot("misc.MT2", hig_cuts, basecuts,  -6, 2, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_6j2b_hHT, gMT2bins_6j2b_hHT, false,  false,  true,  false,  true,   true,  1  ,   false);
    tA->Makeplot("misc.MT2", hig_cuts, basecuts,  -3,-3, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_3b_hHT,   gMT2bins_3b_hHT,   false,  false,  true,  false,  true,   true,  1  ,   false);
*/

//   tA->makePlot("misc.MT2", cuts,   2, 0, -10, trigger , "M_{T2} [GeV]" , 50,  0,  900,     false,  true,  true,  true,  true,   true,  1  ,   false);
//   tA->makePlot("misc.MT2", cuts,   2,-1, -10, trigger , "M_{T2} [GeV]" , 25,  0,  500,     false,  true,  true,  true,  true,   true,  1  ,   false);
//   tA->makePlot("misc.MT2", cuts,  35, 0, -10, trigger , "M_{T2} [GeV]" , 60,  0,  900,     false,  true,  true,  true,  true,   true,  1  ,   false);
//   tA->makePlot("misc.MT2", cuts,  35, 1, -10, trigger , "M_{T2} [GeV]" , 25,  0,  500,     false,  true,  true,  true,  true,   true,  1  ,   false);
//   tA->makePlot("misc.MT2", cuts,  35, 2, -10, trigger , "M_{T2} [GeV]" , 18,  0,  360,     false,  true,  true,  true,  true,   true,  1  ,   false);
//   tA->makePlot("misc.MT2", cuts,  -6, 0, -10, trigger , "M_{T2} [GeV]" , 40,  0,  600,     false,  true,  true,  true,  true,   true,  1  ,   false);
//   tA->makePlot("misc.MT2", cuts,  -6, 1, -10, trigger , "M_{T2} [GeV]" , 25,  0,  500,     false,  true,  true,  true,  true,   true,  1  ,   false);
//   tA->makePlot("misc.MT2", cuts,  -6, 2, -10, trigger , "M_{T2} [GeV]" , 14,  0,  350,     false,  true,  true,  true,  true,   true,  1  ,   false);
//   tA->makePlot("misc.MT2", cuts,  -3,-3, -10, trigger , "M_{T2} [GeV]" , 12,  0,  300,     false,  true,  true,  true,  true,   true,  1  ,   false);

//     tA->MakePlot("misc.MT2", cuts,   2, 0, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_2j0b, gMT2bins_2j0b, false,  false,  true,  false,  true,   true,  1  ,   false);
//     tA->MakePlot("misc.MT2", cuts,   2,-1, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_2j1b, gMT2bins_2j1b, false,  false,  true,  false,  true,   true,  1  ,   false);
//     tA->MakePlot("misc.MT2", cuts,  35, 0, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_3j0b, gMT2bins_3j0b, false,  false,  true,  false,  true,   true,  1  ,   false);
//     tA->MakePlot("misc.MT2", cuts,  35, 1, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_3j1b, gMT2bins_3j1b, false,  false,  true,  false,  true,   true,  1  ,   false);
//     tA->MakePlot("misc.MT2", cuts,  35, 2, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_3j2b, gMT2bins_3j2b, false,  false,  true,  false,  true,   true,  1  ,   false);
//     tA->MakePlot("misc.MT2", cuts,  -6, 0, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_6j0b, gMT2bins_6j0b, false,  false,  true,  false,  true,   true,  1  ,   false);
//     tA->MakePlot("misc.MT2", cuts,  -6, 1, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_6j1b, gMT2bins_6j1b, false,  false,  true,  false,  true,   true,  1  ,   false);
//     tA->MakePlot("misc.MT2", cuts,  -6, 2, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_6j2b, gMT2bins_6j2b, false,  false,  true,  false,  true,   true,  1  ,   false);
//    tA->MakePlot("misc.MT2", cuts,  -3,-3, -10, trigger , "M_{T2} [GeV]" , gNMT2bins_3b,   gMT2bins_3b,   false,  false,  true,  false,  true,   true,  1  ,   false);

//as you can see their are different ways of make plot: 
//Makeplot(...) has cuts and base cuts as well as nbins bins
//makeplot(...) has cuts and base cuts as well as nbins binlow binup
//makePlot(...) has only cuts as well as nbins binlow binup
//MakePlot(...) has only cuts as well as nbins bins
}
