run_ScanAnalysis(TString TYPE, bool DOBRew=true, TString SCALE="", TString OTHER=""){

  // TYPE: mSugraTanB10, T1, T1bbbb ...
  // DoBRew: do btag reweighting
  // SCALE: jes/th/b_up/down

  cout << "Starting... " << endl;

  //gROOT->ProcessLine(".x rootlogon.C");
  //gSystem->CompileMacro("../MT2Code/src/ScanAnalysis.cc","kf");
  gSystem->Load("/shome/leo/Analysis/MT2_V01-05-01_dev/MT2Analysis/Code/MT2AnalysisCode/MT2Code/../MT2Code/src/ScanAnalysis_cc.so");

  TString outputdir = ".";

  //TString TYPE="T1bbbb"; //mSugra or T1, T2bb, etc....

  TString basename = "histoScan_"+TYPE+"_2012StudyAt7TeV_JES42V23_NLL-MinDPhi4-2";
  TString outname = basename;
  
  //TString OTHER = "ISR"; //Generic label...

  TString beffFile="";
  //TString SCALE=""; //th_up/down, b_up/down
  //bool DOBRew = true;

  TString samples   = "../samples/samples_scan_jes.dat";
  //samples   = "../samples/samples_scan_jes_ISR.dat";

  TString NLOfile, LOfile;

  cout << SCALE << endl;
  if(SCALE=="-") SCALE="";

  cout << SCALE << endl;


  if(TYPE=="mSugraTanB10"){
    NLOfile = "combined_cross_section_Errmsugra_m0_m12_10_0_1.txt";
    LOfile = "msugra_m0_m12_10_0_1_LOPythia.txt";
  }
  else{
    NLOfile = "smsXSec.root";
    LOfile = "";
  }

  if(OTHER=="ISR"){
    outname+="_"+OTHER;
    samples   = "../samples/samples_scan_jes_ISR.dat";
    if(SCALE=="jes_up") samples   = "../samples/samples_scan_jesup_ISR.dat";
    if(SCALE=="jes_down") samples   = "../samples/samples_scan_jesdown_ISR.dat";
  }
  else{
    if(SCALE=="jes_up") samples   = "../samples/samples_scan_jesup.dat";
    if(SCALE=="jes_down") samples   = "../samples/samples_scan_jesdown.dat";
  }
  if(SCALE!="") outname+="_"+SCALE;


  cout << "Samples file: " <<samples << endl;

  if(DOBRew){
    outname += "_wBWeight";
    if(OTHER=="ISR")     beffFile = "/shome/leo/Analysis/MT2_V01-05-01_dev/MT2Analysis/Code/MT2AnalysisCode/MT2Code/"+basename+"_ISR.root";
    else
      beffFile = "/shome/leo/Analysis/MT2_V01-05-01_dev/MT2Analysis/Code/MT2AnalysisCode/MT2Code/"+basename+".root";
  }
  
  int verbose = 3;

  //  return;

  
  ScanAnalysis *tA = new ScanAnalysis(outputdir);
  
  tA->setVerbose(verbose);
  //tA->PDF_UNC_FILE = "";
  
  //if(TYPE=="mSugra")

  //return;
  tA->init(TYPE, samples, NLOfile, LOfile,SCALE,DOBRew,beffFile, OTHER );
  
  std::ostringstream cutStream;
  cutStream 
    //<< "   (Susy.M0>180 && Susy.M0<220) && (Susy.M12>580 && Susy.M12<620) && " 
    //<< "Sum$(jet[].lv.Pt()>40)>2"                            << "&&"
    << "misc.MET>=30 "                            << "&&"
    << " misc.MT2 > 80"                           << "&&"
    << "NJetsIDLoose40 > 1 && "
    << "misc.HT > 750 "                          << "&&"
    //<< "((misc.MinMetJetDPhi >0.3&&NBJets==0)||NBJets>=1)" << "&&"
    << "misc.MinMetJetDPhi4 >0.3 " << "&&"
 //<< "misc.HT > 950 "                          << "&&"
    //	      << "misc.HT > 750 && misc.HT<=950 "                          << "&&"
    //<< "misc.caloHT50_ID >750 "                  << "&&"
    << "misc.Jet0Pass ==1"                       << "&&"
    << "misc.Jet1Pass ==1"                       << "&&"
    << "misc.PassJetID ==1"                      << "&&"
    << "misc.Vectorsumpt < 70"                   << "&&"
    //<< "misc.MinMetJetDPhi <=0.3"                 << "&&"
    //<< "misc.HBHENoiseFlagIso == 0"                 << "&&"
    //      << "misc.CrazyHCAL==0 && "
    //    << "misc.LeadingJPt >150"                    << "&&"
    //      << " misc.CSCTightHaloID==0 && " 
    << "misc.SecondJPt  >100";//                    << "&&"
  //<< "NBJets >0"                               << "&&"
  //	  << "misc.MT2 > 200 && misc.MT2<400"          << "&&"
  //	  << "misc.MT2 > 100 && misc.MT2<150"          << "&&"
  //	  << "misc.CrazyHCAL==0";
  
  TString cuts = cutStream.str().c_str();
  

  tA->Analysis(outname+".root",cuts);
  
  
  
}
