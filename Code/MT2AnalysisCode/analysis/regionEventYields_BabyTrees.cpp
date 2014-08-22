#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

//#include "MT2tree.hh"
#include "helper/Utilities.hh"
#include "Utilities.hh"

#include <TTreeFormula.h>

#include "interface/MT2Common.h"
#include "interface/MT2Region.h"
#include "interface/MT2RegionAnalysisUtilities.h"

#define treeProducerSusyFullHad_cxx
#include "treeProducerSusyFullHad.h"


//bool fFast = true;
//bool fISRreweight = true;
//bool fbTagReweight = true;


MT2YieldAnalysis* computeYield( const MT2SampleBaby_basic& sample, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions );
std::vector<MT2YieldPreAnalysis*> getEventYield( std::vector<std::string> fakeID, TTree* tree, MT2SampleBaby_basic sample, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions );
//float getISRCorrection( MT2tree* fMT2tree, const MT2SampleBaby_basic& sample );
//void getBTagScaleFactor( MT2tree* fMT2tree, MT2Region* region, float& btagSF, float& btagSFerr );
MT2YieldAnalysis* mergeYields( std::vector<MT2YieldAnalysis*> EventYield, const std::string& n1, const std::string& n2="", const std::string& n3="", const std::string& n4="", const std::string& n5="" );
std::vector<TH1D*> getYieldHistos( const std::string& prefix, const std::string& fakeID, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions, MT2YieldAnalysis* EventYield_tot, MT2YieldAnalysis* EventYield_bg, std::ofstream& logfile );
std::vector<TH1D*> getSimTruthYieldHistos( const std::string& prefix, const std::string& fakeID, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions, MT2YieldAnalysis* EventYield_tot );




int main( int argc, char* argv[] ) {


  if( argc!=2 ) {
    std::cout << "USAGE: ./regionEventYields_BabyTrees [samplesFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }

  std::string sampleName(argv[1]);

  std::string samplesFileName = "samples/samples_" + sampleName + ".dat";
  std::cout << std::endl << std::endl;
  std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;

  std::vector<MT2SampleBaby_basic> fSamples = MT2Common::loadSamplesBaby_basic(samplesFileName);



  std::ostringstream HLT_metHT;
  HLT_metHT << "";
  //HLT_metHT << "( ( ("
  //    << "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
  //    << "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )"
  //    << "||("
  //    << "trigger.HLT_PFHT350_PFMET100_v3==1 || trigger.HLT_PFHT350_PFMET100_v4==1 || trigger.HLT_PFHT350_PFMET100_v5==1 || "
  //    << "trigger.HLT_PFHT350_PFMET100_v6==1 || trigger.HLT_PFHT350_PFMET100_v7==1 || trigger.HLT_PFNoPUHT350_PFMET100_v1==1 || "
  //    << "trigger.HLT_PFNoPUHT350_PFMET100_v3==1 || trigger.HLT_PFNoPUHT350_PFMET100_v4==1 ) ) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";

  std::ostringstream HLT_HT;
  HLT_HT << "";
  //HLT_HT << "( ("
  //    << "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || "
  //    << "trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
  //    << "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";



  std::vector<MT2HTRegion> HTRegions;
  HTRegions.push_back(MT2HTRegion("lowHT",    450.,    750., 200., HLT_metHT.str()));
  HTRegions.push_back(MT2HTRegion("mediumHT", 750.,   1200.,  30., HLT_HT.str()));
  HTRegions.push_back(MT2HTRegion("highHT",  1200., 100000.,  30., HLT_HT.str()));

  std::vector<MT2SignalRegion> signalRegions;
  signalRegions.push_back(MT2SignalRegion(2, 2, 0, 0));  // 2j0b
  signalRegions.push_back(MT2SignalRegion(2, 2, 1, 2));  // 2j1to2b
  signalRegions.push_back(MT2SignalRegion(3, 5, 0, 0));  // 3to5j0b
  signalRegions.push_back(MT2SignalRegion(3, 5, 1, 1));  // 3to5j1b
  signalRegions.push_back(MT2SignalRegion(3, 5, 2, 2));  // 3to5j2b
  signalRegions.push_back(MT2SignalRegion(6, -1, 0, 0));  // 6j0b
  signalRegions.push_back(MT2SignalRegion(6, -1, 1, 1));  // 6j1b
  signalRegions.push_back(MT2SignalRegion(6, -1, 2, 2));  // 6j2b
  signalRegions.push_back(MT2SignalRegion(-1, -1, 3, -1));  // 3b


  TH1::AddDirectory(kFALSE); // stupid ROOT memory allocation needs this

  
  std::vector<MT2YieldAnalysis*> EventYield;
  for( unsigned i=0; i<fSamples.size(); ++i )
    EventYield.push_back( computeYield( fSamples[i], HTRegions, signalRegions ) );
  

  system( "rm tmp.root" );



  std::cout << "-> Done looping on samples. Start merging." << std::endl;

  MT2YieldAnalysis* EventYield_data  = mergeYields( EventYield, "HT-Data", "MET-Data" );
  //MT2YieldAnalysis* EventYield_data  = mergeYields( EventYield, "HT-Data");
  MT2YieldAnalysis* EventYield_top   = mergeYields( EventYield, "Top" );
  MT2YieldAnalysis* EventYield_qcd   = mergeYields( EventYield, "QCD" );
  MT2YieldAnalysis* EventYield_wjets = mergeYields( EventYield, "WJets" );
  MT2YieldAnalysis* EventYield_other = mergeYields( EventYield, "DY", "VV" );

  MT2YieldAnalysis* EventYield_signal = mergeYields( EventYield, "Signal" );
 	       
  MT2YieldAnalysis* EventYield_topW   = new MT2YieldAnalysis(*EventYield_top  + *EventYield_wjets);
  MT2YieldAnalysis* EventYield_bg     = new MT2YieldAnalysis(*EventYield_qcd + *EventYield_other);
  MT2YieldAnalysis* EventYield_allMC  = new MT2YieldAnalysis(*EventYield_topW + *EventYield_qcd + *EventYield_other);



  std::cout << "-> Done merging. Start computing event yields." << std::endl;


  std::string outputdir = "EventYields_WJets-HT400to600";
  system(Form("mkdir -p %s", outputdir.c_str()));

  TFile* outfile = TFile::Open(Form("%s/EventYields_%s.root", outputdir.c_str(), sampleName.c_str()), "recreate");
  outfile->cd();
  
  std::ofstream logfile;
  logfile.open(Form("%s/EventYields_%s.txt", outputdir.c_str(), sampleName.c_str()));
  

  // For data yield:
  //std::vector<TH1D*> vh1_data = getYieldHistos( "EventYield_data", "", HTRegions, signalRegions, EventYield_data, EventYield_allMC, logfile );
  
  // For topW yield:
  std::vector<TH1D*> vh1_data = getYieldHistos( "EventYield_topW", "", HTRegions, signalRegions, EventYield_topW, EventYield_bg, logfile );
  
  // For signal yield:
  //std::vector<TH1D*> vh1_data = getYieldHistos( "EventYield_signal", "", HTRegions, signalRegions, EventYield_signal, EventYield_allMC, logfile ); 

  // For all SM yield:
  std::vector<TH1D*> vh1_mc   = getYieldHistos( "EventYield_MC", "", HTRegions, signalRegions, EventYield_allMC, EventYield_bg, logfile );
  //std::vector<TH1D*> vh1_mc   = getYieldHistos( "EventYield_WJets", "", HTRegions, signalRegions, EventYield_wjets, EventYield_bg, logfile );

  std::vector<TH1D*> vh1_sim  = getSimTruthYieldHistos( "SimulationTruthEventYield", "", HTRegions, signalRegions, EventYield_allMC );
  //std::vector<TH1D*> vh1_sim  = getSimTruthYieldHistos( "SimulationTruthEventYield", "", HTRegions, signalRegions, EventYield_signal );
 

  for( unsigned i=0; i<HTRegions.size(); ++i ) {

    vh1_data[i]->Write();
    
    vh1_mc  [i]->Write();
    
    vh1_sim [i]->Write();
 
  }


  outfile->Close();

  std::cout << "-> Done computing preditions. Written to file: " << outfile->GetName() << std::endl;

  return 0;

}




MT2YieldAnalysis* computeYield( const MT2SampleBaby_basic& sample, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions ) {


  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("treeProducerSusyFullHad");
  
  treeProducerSusyFullHad myTree;
  myTree.Init(tree);


  std::ostringstream preselectionStream;
  preselectionStream << " " 
    << "(nTaus20==0 && nMuons10==0 && nElectrons10==0)" << " && "
    << "(nVert > 0)"                      << " && "
    << "(nJet40 > 1)"                     << " && "
    << "(jet_pt[1] > 100)"                << " && "
    << "(ht > 450)"                       << " && "
    << "((ht < 750 && met_pt > 200) || (ht > 750 && met_pt > 30))" << " && "
    << "(deltaPhiMin > 0.3)"              << " && " 
    << "(diffMetMht < 70)";
  //if( fFast ) preselectionStream << " && (misc.MT2>=100.) ";//lowest border in MT2

  /*
  // add "base" cuts:
  preselectionStream << " && " 
    << "(misc.PassJet40ID ==1)"                             << " && "
    << "((misc.HBHENoiseFlag == 0 || misc.ProcessID==10))"  << " && " // signal samples (fastsim) do not have this filter
    << "(misc.CSCTightHaloIDFlag == 0)"                     << " && "
    << "(misc.trackingFailureFlag==0)"                      << " && "
    << "(misc.eeBadScFlag==0)"                              << " && "
    << "(misc.EcalDeadCellTriggerPrimitiveFlag==0)"         << " && "
    << "(misc.TrackingManyStripClusFlag==0)"                << " && "
    << "(misc.TrackingTooManyStripClusFlag==0)"             << " && "
    << "(misc.TrackingLogErrorTooManyClustersFlag==0)"      << " && "
    << "(misc.CrazyHCAL==0)";
  preselectionStream << " && (misc.MET/misc.CaloMETRaw<=2.)"; //HO cut
  */
  

  TString preselection = preselectionStream.str().c_str();
  TString cuts = preselection;

  TFile* tmpFile = TFile::Open("tmp.root", "recreate");
  tmpFile->cd();
  TTree* tree_reduced = tree->CopyTree(cuts);

  std::vector<std::string> fakeID;
  fakeID.push_back("");
    
  std::vector<MT2YieldPreAnalysis*> v_EventYield = getEventYield( fakeID, tree_reduced, sample, HTRegions, signalRegions );

  MT2YieldAnalysis* EventYield = new MT2YieldAnalysis(sample.sname);
  for( unsigned i=0; i<fakeID.size(); ++i )
    EventYield->f[ fakeID[i].c_str() ] = v_EventYield[i];
  
  delete tree;
  //delete tree_reduced;

  tmpFile->Close();
  delete tmpFile;

  file->Close();
  delete file;
  
  return EventYield;

}




std::vector<MT2YieldPreAnalysis*> getEventYield( std::vector<std::string> fakeID, TTree* tree, MT2SampleBaby_basic sample, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions ) {


  bool isData = false;

  // global sample weight:
  //Double_t weight = sample.xsection * sample.kfact * sample.lumi / (sample.nevents*sample.PU_avg_weight);
  Double_t weight = sample.xsection * sample.kfact * sample.lumi / (sample.nevents);
  cout << endl << "Weight " << weight <<endl; 

  // return event yield:
  std::vector<MT2YieldPreAnalysis*> EventYield;
  
  for( unsigned i=0; i<fakeID.size(); ++i ) {
    MT2YieldPreAnalysis* this_EventYield = new MT2YieldPreAnalysis( sample.sname, HTRegions, signalRegions );
    EventYield.push_back( this_EventYield );
  }

  //MT2YieldAnalysis* this_yield = new MT2YieldAnalysis( sample.sname, HTRegions, signalRegions );
  ////yield.push_back( this_yield );
  

  // instantiate HLT TTreeFormulas
  std::vector<TTreeFormula*> ttf_hlt;
  if( isData ) {
    for( unsigned iHT=0; iHT<HTRegions.size(); ++iHT ) {
      TTreeFormula* this_ttf = new TTreeFormula(Form("ttf_hlt_%s", HTRegions[iHT].getName().c_str()), HTRegions[iHT].HLT_selection.c_str(), tree);
      ttf_hlt.push_back(this_ttf);
    }
  }


  //MT2tree* fMT2tree = new MT2tree();
  //tree->SetBranchAddress("MT2tree", &fMT2tree);

  treeProducerSusyFullHad myTree;
  myTree.Init(tree);


  int nentries = tree->GetEntries();

  float fullweight = weight;
  float fullweight_btagUp = weight;
  float fullweight_btagDown = weight;



  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    myTree.GetEntry(iEntry);


    //float ht   = fMT2tree->misc.HT;
    //float met  = fMT2tree->misc.MET;
    //float mt2  = fMT2tree->misc.MT2;
    //int njets  = fMT2tree->NJetsIDLoose40;
    //int nbjets = fMT2tree->NBJets40CSVM;

    float ht   = myTree.ht;
    float met  = myTree.met_pt;
    float mt2  = myTree.mt2;
    int njets  = myTree.nJet40;
    int nbjets = myTree.nBJet40;

    fullweight = weight;
    fullweight_btagUp = weight;
    fullweight_btagDown = weight;


    bool foundRegion = false;

    for( unsigned iHT=0; iHT<HTRegions.size() && !(foundRegion); ++iHT ) {

      if( ht<HTRegions[iHT].htMin || ht>=HTRegions[iHT].htMax || met<HTRegions[iHT].metMin ) continue;

      if( isData ) { // apply HLT
        if( ttf_hlt[iHT]->EvalInstance()==0. ) continue;
      }


      for( unsigned iSR=0; iSR<signalRegions.size() && !(foundRegion); ++iSR ) {

        if( signalRegions[iSR].nJetsMin  >= 0  &&  njets  < signalRegions[iSR].nJetsMin ) continue;
        if( signalRegions[iSR].nJetsMax  >= 0  &&  njets  > signalRegions[iSR].nJetsMax ) continue;
        if( signalRegions[iSR].nBJetsMin >= 0  &&  nbjets < signalRegions[iSR].nBJetsMin ) continue;
        if( signalRegions[iSR].nBJetsMax >= 0  &&  nbjets > signalRegions[iSR].nBJetsMax ) continue;


        // found region

        MT2Region thisRegion( &HTRegions[iHT], &signalRegions[iSR] );


        //// compute weight
        //if( !isData ) {
	//
        //  //float puweight =  fMT2tree->pileUp.Weight;
        //  float puweight  = 1.0;
	//  //float isrweight = getISRCorrection( fMT2tree, sample );
	//  float isrweight = 1.0;
        //  float btagSF = 1.0, btagSFerr = 0.0;
        //  //getBTagScaleFactor( fMT2tree, &thisRegion, btagSF, btagSFerr );
	//
        //  fullweight          = weight * puweight * isrweight * btagSF;
        //  fullweight_btagUp   = weight * puweight * isrweight *(btagSF + btagSFerr);
        //  fullweight_btagDown = weight * puweight * isrweight *(btagSF - btagSFerr);
	//
        //}
	
	fullweight = weight;

	
	for( unsigned iF=0; iF<fakeID.size(); ++iF ) {

	  MT2RegionAnalysis*     thisYield         = EventYield[iF]->getRegion   ( thisRegion.getName() );
	  MT2RegionAnalysisBase* thisSimTruthYield = EventYield[iF]->getRegionGen( thisRegion.getName() );
	  
	  
	  //thisYield->yield->Fill( mt2, fullweight );
	  //thisYield->yield_btagUp->Fill( mt2, fullweight_btagUp );
	  //thisYield->yield_btagDown->Fill( mt2, fullweight_btagDown );
	  //
	  //thisSimTruthYield->yield->Fill( mt2, fullweight );

	  thisYield->yield->Fill( 1, weight );
          thisYield->yield_btagUp->Fill( 1, fullweight_btagUp );
          thisYield->yield_btagDown->Fill( 1, fullweight_btagDown );
	  
          thisSimTruthYield->yield->Fill( 1, weight );

	}
	
	foundRegion = true; // found region, move on
	
      }  // for iSR
      
    }  // for iHT
    
  } // for entries

  for( unsigned i=0; i<EventYield.size(); ++i )
    EventYield[i]->addOverflow();
  
  delete tree;
  
  return EventYield;

}




MT2YieldAnalysis* mergeYields( std::vector<MT2YieldAnalysis*> EventYield, const std::string& n1, const std::string& n2, const std::string& n3, const std::string& n4, const std::string& n5 ) {

  std::vector<std::string> snames;
  snames.push_back(n1);
  snames.push_back(n2);
  snames.push_back(n3);
  snames.push_back(n4);
  snames.push_back(n5);

  std::string newname = (snames.size()>0) ? (snames[0]) : "";
  for( unsigned i=1; i<snames.size(); ++i ) newname += "_" + snames[i];
  newname = "merge_" + newname;

  MT2YieldAnalysis* return_EventYield = new MT2YieldAnalysis(newname);

  for( unsigned i=0; i<EventYield.size(); ++i ) {

    for( unsigned iname=0; iname<snames.size(); ++iname ) {

      if( EventYield[i]->SName == snames[iname] ) {

	return_EventYield->add(*(EventYield[i]));
	
      }

    } // for snames

  } // for EventYield


  return return_EventYield;

}





std::vector<TH1D*> getYieldHistos( const std::string& prefix, const std::string& fakeID, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions, MT2YieldAnalysis* EventYield_tot, MT2YieldAnalysis* EventYield_bg, std::ofstream& logfile ) {


  int nHistos = HTRegions.size();
  int nBins = signalRegions.size();


  std::vector<TH1D*> histos;
  
  std::cout << std::endl << std::endl;
  std::cout << "Event yield for sample " << prefix.c_str() << std::endl;
  
  logfile << std::endl << std::endl <<  "Event yield for sample " << prefix.c_str() << std::endl;
  
  for( int i=0; i<nHistos; ++i ) {

    TH1D* h1 = new TH1D(Form("%s_%s", prefix.c_str(), HTRegions[i].name.c_str()), "", nBins, 0., nBins);
    h1->Sumw2();

    //float EventYield_sum = 0;
    //float statErr_sum = 0;

    std::cout << std::endl << std::endl;
    std::cout << "HT region: " << HTRegions[i].name.c_str() << std::endl;
    logfile << std::endl << "HT region: " << HTRegions[i].name.c_str() << std::endl;

    for( int j=0; j<nBins; ++j ) {
  
      int iBin = j+1;

      MT2Region thisRegion( &HTRegions[i], &signalRegions[j] );

      h1->GetXaxis()->SetBinLabel(iBin, signalRegions[j].getName().c_str());

      //int tot = 0;
      //
      //if(EventYield_tot!=0){
      //	MT2RegionAnalysis* thisYield = EventYield_tot->getRegion( thisRegion.getName() );
      //	
      //	cout<< "DEBUG1" << endl;
      //	
      //  (TH1D*) (thisYield->yield)->Clone();
      //	
      //	cout<< "DEBUG2"<< endl;
      //	
      //}
      //else ;

      double err_histo = 0;
      float  tot = (EventYield_tot!=0 && EventYield_tot->f[fakeID.c_str()]!=0) ? EventYield_tot->f[fakeID.c_str()]->getRegion(thisRegion.getName())->yield->IntegralAndError(0, 2, err_histo) : 0.;
      //float  bg = (EventYield_bg !=0 && EventYield_bg ->f[fakeID.c_str()]!=0) ? EventYield_bg ->f[fakeID.c_str()]->getRegion(thisRegion.getName())->yield->Integral() : 0.;
      //float  bg_btagUp = (EventYield_bg !=0 && EventYield_bg ->f[fakeID.c_str()]!=0) ? EventYield_bg ->f[fakeID.c_str()]->getRegion(thisRegion.getName())->yield_btagUp->Integral() : 0.;
      //float  bg_btagDown = (EventYield_bg !=0 && EventYield_bg ->f[fakeID.c_str()]!=0) ? EventYield_bg ->f[fakeID.c_str()]->getRegion(thisRegion.getName())->yield_btagDown->Integral() : 0.;

      //float EventYield = tot-bg;
      
      float EventYield = tot;
      h1->SetBinContent( iBin, EventYield );
      h1->SetBinError( iBin, err_histo );

      //cout << "tot " << tot << endl;
      
      ////float  tot = (EventYield_tot!=0) ? EventYield_tot->getRegion(thisRegion.getName())->yield->Integral() : 0.;      
      //float  bg  = (EventYield_bg !=0) ? EventYield_bg->getRegion(thisRegion.getName())->yield->Integral() : 0.;
      //float  bg_btagUp   = (EventYield_bg !=0) ? EventYield_bg->getRegion(thisRegion.getName())->yield_btagUp->Integral() : 0.;
      //float  bg_btagDown = (EventYield_bg !=0) ? EventYield_bg->getRegion(thisRegion.getName())->yield_btagDown->Integral() : 0.;
      //
      //float EventYield = (tot-bg);
      //
      //cout<< "Credo il problema sia qui sopra.";
      //
      //h1->SetBinContent( iBin, EventYield );

      std::cout << std::endl;
      std::cout << "Signal region: " << signalRegions[j].getName().c_str() << std::endl;  
     
      //std::cout << "tot: " << tot << std::endl;
      //std::cout << "bg: " << bg << std::endl;
      
      std::cout << "yield: " << EventYield << std::endl;
      std::cout << "hist_err: " << err_histo << std::endl;
           
      if(EventYield >= 10)
	logfile << fixed << setprecision(1) << EventYield << " \\pm " << err_histo << " & ";
      else if(EventYield < 10)
	logfile << fixed << setprecision(2) << EventYield << " \\pm " << err_histo << " & ";

      //float statErr = sqrt(tot);
      //float sysErr_dbg = 0.5*bg; // 50% uncert on BG estimate
      //float bg_SFerrUp = fabs(bg_btagUp-bg);
      //float bg_SFerrDown = fabs(bg_btagDown-bg);
      //float bg_SFerr = 0.5*(bg_SFerrUp+bg_SFerrDown);
      //float sysErr_btag = bg_SFerr;

      //float sysErr_dLL = nWT_goodrecoevt_dLL*(1.*effLept_err/(effLept*effLept*effMT));//additional error due to double lost with 100% uncert

      //float sysErr = sqrt( sysErr_dbg*sysErr_dbg + sysErr_btag*sysErr_btag );
      //std::cout << "sysErr_dbg: " << sysErr_dbg << std::endl;
      //std::cout << "sysErr_btag: " << sysErr_btag << std::endl;
      //std::cout << "sysErr: " << sysErr<< std::endl;
      //std::cout << "statErr: " << statErr<< std::endl; 

      //h1->SetBinError( iBin, sqrt(statErr*statErr + sysErr*sysErr) );
      h1->SetBinError( iBin, err_histo );

    }

    histos.push_back(h1);

  } // for

  return histos;

}





std::vector<TH1D*> getSimTruthYieldHistos( const std::string& prefix, const std::string& fakeID, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions, MT2YieldAnalysis* EventYield_tot ) {


  int nHistos = HTRegions.size();
  int nBins = signalRegions.size();


  
  std::vector<TH1D*> histos;

  for( int i=0; i<nHistos; ++i ) {

    TH1D* h1 = new TH1D(Form("%s_%s", prefix.c_str(), HTRegions[i].name.c_str()), "", nBins, 0., nBins);
    h1->Sumw2();

    for( int j=0; j<nBins; ++j ) {
  
      int iBin = j+1;

      MT2Region thisRegion( &HTRegions[i], &signalRegions[j] );

      h1->GetXaxis()->SetBinLabel(iBin, signalRegions[j].getName().c_str());

      float totYield = (EventYield_tot!=0) ? EventYield_tot->f[fakeID.c_str()]->getRegionGen(thisRegion.getName())->yield->Integral() : 0.;

      h1->SetBinContent( iBin, totYield );

    }

    histos.push_back(h1);

  } // for

  return histos;

}





//float getISRCorrection( MT2tree* fMT2tree, const MT2SampleBaby_basic& sample ) {
//
//  //assign all samples to its type
//  string sampletype = (string)sample.type;
//  if(sampletype==(string)"mc"){
//    if(sample.sname=="QCD")         sampletype = (string)"QCD";
//    else if(sample.sname=="Wtolnu") sampletype = (string)"WJets";
//    else if(sample.sname=="DY")     sampletype = (string)"ZJets";
//    else if(sample.name=="TTbar")   sampletype = (string)"TTbar";
//    else if(sample.name=="TTbar_Madgraph0l")   sampletype = (string)"TTbar";
//    else if(sample.name=="TTbar_Madgraph1l")   sampletype = (string)"TTbar";
//    else if(sample.name=="TTbar_Madgraph2l")   sampletype = (string)"TTbar";
//    else if(sample.sname=="Top")    sampletype = (string)"SingleTop";//no ttbar, includes TTZ, TTW
//    else sampletype = (string)"Other";
//  }
//
//  TLorentzVector hardgenlv; hardgenlv.SetPtEtaPhiM(0,0,0,0);
//  if(sampletype=="WJets"){
//    bool foundW(false);
//    for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){//lepton from W
//      int ID  =abs(fMT2tree->genlept[ngl].ID);
//      int MID =abs(fMT2tree->genlept[ngl].MID);
//      if((ID == 11 || ID == 12 || ID == 13 || ID == 14 || ID == 15 || ID == 16) && MID==24){
//        hardgenlv = fMT2tree->genlept[ngl].Mlv; foundW = true;
//      }
//      if(foundW) break;
//    }
//    if(!foundW){
//      for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){//lepton from tau, tau from W
//        int ID  =abs(fMT2tree->genlept[ngl].ID);
//        int MID =abs(fMT2tree->genlept[ngl].MID);
//        int GMID=abs(fMT2tree->genlept[ngl].GMID);
//        if((ID == 11 || ID == 12 || ID == 13 || ID == 14 || ID == 16) && MID==15 && GMID==24){
//          hardgenlv = fMT2tree->genlept[ngl].Mlv; foundW = true;
//        }
//        if(foundW) break;
//      }
//    }
//  } else if(sampletype=="ZJets"){
//    hardgenlv = fMT2tree->GenZ[0];
//  } else if(sampletype=="TTbar"){
//    TLorentzVector top1(0.,0.,0.,0.), top2(0.,0.,0.,0.);
//    bool top1f(false), top2f(false);
//    for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
//      int id   = abs(fMT2tree->genlept[ngl].ID);
//      if(id!=5) continue;
//      int mid  = fMT2tree->genlept[ngl].MID;//from b
//      if(mid==6&&top1f) continue;
//      else if(mid==6) { top1 = fMT2tree->genlept[ngl].Mlv; top1f = true; }
//      if(mid==-6&&top2f) continue;
//      else if(mid==-6){ top2 = fMT2tree->genlept[ngl].Mlv; top2f = true; }
//      if(top1f&&top2f) {
//        hardgenlv = top1+top2;
//        break;
//      }
//    }
//  } else if(sampletype=="SingleTop"){
//    TString samplename_tstr(sample.name);
//    if(samplename_tstr.Contains("tW")){//t + W
//      TLorentzVector top(0.,0.,0.,0.), W(0.,0.,0.,0.);
//      bool topf(false), Wf(false);
//      for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
//        int id    = abs(fMT2tree->genlept[ngl].ID);
//        int mid   = abs(fMT2tree->genlept[ngl].MID);//from b
//        int gmid  = abs(fMT2tree->genlept[ngl].GMID);//from b
//        if(mid==6&&topf) continue;
//        else if(mid==6) { top = fMT2tree->genlept[ngl].Mlv; topf = true; }
//        if(mid==24&&gmid!=6&&Wf) continue;
//        if((id == 11 || id == 12 || id == 13 || id == 14 || id == 15 || id == 16) && mid==24 && gmid!=6){
//          W = fMT2tree->genlept[ngl].Mlv; Wf = true;
//        }
//        if(topf&&Wf){
//          hardgenlv = top+W;
//          break;
//        }
//      }
//      if(!Wf){//this might be wrong - but influence negligible anyway
//        for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
//          int id    = abs(fMT2tree->genlept[ngl].ID);
//          int mid   = abs(fMT2tree->genlept[ngl].MID);//from b
//          int gmid  = abs(fMT2tree->genlept[ngl].GMID);//from b
//          if(mid==6||gmid==6) continue;
//          if(gmid==24&&gmid==15&&Wf) continue;
//          if((id == 11 || id == 12 || id == 13 || id == 14 || id == 15 || id == 16) && mid==15 && gmid==24){
//            W = fMT2tree->genlept[ngl].Mlv; Wf = true;
//          }
//          if(topf&&Wf){
//            hardgenlv = top+W;
//            break;
//          }
//        }
//      }
//    
//    } else {
//      TLorentzVector top(0.,0.,0.,0.);
//      bool topf(false), Wf(false);
//      for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
//        int id    = abs(fMT2tree->genlept[ngl].ID);
//        int mid   = abs(fMT2tree->genlept[ngl].MID);//from b
//        int gmid  = abs(fMT2tree->genlept[ngl].GMID);//from b
//        if(mid==6&&topf) continue;
//        else if(mid==6) { top = fMT2tree->genlept[ngl].Mlv; topf = true; }
//        if(topf){
//          hardgenlv = top;
//          break;
//        }
//      }
//    }
//  }
//
//  float ISRweight=1.;
//  if(hardgenlv.Pt()>250.) ISRweight = 0.8;
//  else if(hardgenlv.Pt()>150.) ISRweight = 0.9;
//  else if(hardgenlv.Pt()>120.) ISRweight = 0.95;
//  else                         ISRweight = 1.;
//      
//  return ISRweight;
//
//}
//
//
//
//
//void getBTagScaleFactor( MT2tree* fMT2tree, MT2Region* region, float& btagSF, float& btagSFerr ) {
//
//  int nbjetsmin = region->nBJetsMin();
//  int nbjetsmax = region->nBJetsMax();
//
//  if( nbjetsmin != nbjetsmax ) {  // use "ge"
//
//    if( nbjetsmin >= 3 ) {
//
//      btagSF = fMT2tree->SFWeight.BTagCSV40ge3;
//      btagSFerr = fMT2tree->SFWeight.BTagCSV40ge3Error;
//
//    } else if( nbjetsmin == 2 ) {
//
//      btagSF = fMT2tree->SFWeight.BTagCSV40ge2;
//      btagSFerr = fMT2tree->SFWeight.BTagCSV40ge2Error;
//
//    } else if( nbjetsmin == 1 ) {
//
//      btagSF = fMT2tree->SFWeight.BTagCSV40ge1;
//      btagSFerr = fMT2tree->SFWeight.BTagCSV40ge1Error;
//
//    } else { // use "eq"
//
//      btagSF = 1.;
//      btagSFerr = 0.;
//
//    }
//
//  } else { // if nbjetsmin==nbjetsmax
//
//    if( nbjetsmin >= 3 ) {
//
//      btagSF = fMT2tree->SFWeight.BTagCSV40eq3;
//      btagSFerr = fMT2tree->SFWeight.BTagCSV40eq3Error;
//
//    } else if( nbjetsmin == 2 ) {
//
//      btagSF = fMT2tree->SFWeight.BTagCSV40eq2;
//      btagSFerr = fMT2tree->SFWeight.BTagCSV40eq2Error;
//
//    } else if( nbjetsmin == 1 ) {
//
//      btagSF = fMT2tree->SFWeight.BTagCSV40eq1;
//      btagSFerr = fMT2tree->SFWeight.BTagCSV40eq1Error;
//
//    } else {
//
//      btagSF = 1.;
//      btagSFerr = 0.;
//
//    }
//
//  }
//
//}
