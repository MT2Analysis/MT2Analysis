#include <iostream>
#include <sstream>

#include "MT2tree.hh"
#include "helper/Utilities.hh"
#include "Utilities.hh"

#include "interface/MT2Common.h"
#include "interface/MT2Region.h"
#include "interface/MT2LostLeptonUtilities.h"



bool fFast = false;
bool fISRreweight = true;
bool fbTagReweight = true;
bool fIncludeTaus = true;





MT2LostLeptonEstimate* computeLostLepton( const MT2Sample& sample, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions );
//void computeSignalRegionEstimate( TFile* outfile, const MT2Sample& sample, const std::string& leptType, int njet_min, int njet_max, int nbjet_min, int nbjet_max, float ht_min, float ht_max );
MT2LeptonTypeLLEstimate* getLeptonTypeLLEstimate( const std::string& leptType, TTree* tree, MT2Sample sample, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions );
//float getISRCorrection( MT2tree* fMT2tree, const MT2Sample& sample, const std::string& sampletype );
//void getBTagSF( int njets, int nbjets, double& btagSF, double& btagSFerr );
MT2LostLeptonEstimate* mergeEstimates( std::vector<MT2LostLeptonEstimate*> llest, const std::string& n1, const std::string& n2="", const std::string& n3="", const std::string& n4="", const std::string& n5="" );



int main() {



  std::string samplesFileName = "samples_HT_filter_prova_fast.dat";
  std::vector<MT2Sample> fSamples = MT2Common::loadSamples(samplesFileName);

 
  std::vector<MT2HTRegion> HTRegions;
  HTRegions.push_back(MT2HTRegion("lowHT",    450.,    750., 200.));
  HTRegions.push_back(MT2HTRegion("mediumHT", 750.,   1200.,  30.));
  HTRegions.push_back(MT2HTRegion("highHT",  1200., 100000.,  30.));

  std::vector<MT2SignalRegion> signalRegions;
  signalRegions.push_back(MT2SignalRegion(3, 3, 1, 1));  // 3j1b
  signalRegions.push_back(MT2SignalRegion(2, 2, 0, 0));  // 2j0b
  signalRegions.push_back(MT2SignalRegion(3, 3, 0, 0));  // 3j0b


  TFile* outfile = TFile::Open("outfile.root", "recreate");
  
  std::vector<MT2LostLeptonEstimate*> llest;
  for( unsigned i=0; i<fSamples.size(); ++i )
    llest.push_back( computeLostLepton( fSamples[i], HTRegions, signalRegions ) );
  

  MT2LostLeptonEstimate* ll_data  = mergeEstimates( llest, "HT-Data" );
  MT2LostLeptonEstimate* ll_top   = mergeEstimates( llest, "Top" );
  //MT2LostLeptonEstimate* ll_qcd   = mergeEstimates( llest, "QCD" );
  //MT2LostLeptonEstimate* ll_wjets = mergeEstimates( llest, "Wtolnu" );
  //MT2LostLeptonEstimate* ll_other = mergeEstimates( llest, "DY", "VV" );
  
  //MT2LostLeptonEstimate* ll_topW   = new MT2LostLeptonEstimate(*ll_top  + *ll_wjets);
  //MT2LostLeptonEstimate* ll_bg     = new MT2LostLeptonEstimate(*ll_qcd  + *ll_other);
  //MT2LostLeptonEstimate* ll_allMC  = new MT2LostLeptonEstimate(*ll_topW + *ll_bg);

  MT2Region r( &HTRegions[0], &signalRegions[0] );

  std::cout << "yields:" << std::endl;
  std::cout << llest[1]->l["Ele"] << std::endl;
  std::cout << llest[1]->l["Ele"]->getRecoRegion(r.getName()) << std::endl;
  std::cout << llest[1]->l["Ele"]->getRecoRegion(r.getName())->yield << std::endl;
  std::cout << "all: " << llest[1]->l["Ele"]->getRecoRegion(r.getName())->yield->Integral() << std::endl;
  std::cout << ll_top->l["Ele"] << std::endl;
  std::cout << ll_top->l["Ele"]->getRecoRegion(r.getName()) << std::endl;
  std::cout << ll_top->l["Ele"]->getRecoRegion(r.getName())->yield << std::endl;
  std::cout << "top: " << ll_top->l["Ele"]->getRecoRegion(r.getName())->yield->Integral() << std::endl;
  //std::cout << "data: " << ll_data->l["Ele"]->getRecoRegion(r.getName())->yield->Integral() << std::endl;

  // then create summary histos and write to outfile
//  std::vector<TH1D*> histos_ele_reco = createPredictionHistos("Prediction"     , "Ele", llest );
//  std::vector<TH1D*> histos_gen  = createPredictionHistos("SimulationTruth", leptType, HTRegions, signalRegions);


  return 0;

}




MT2LostLeptonEstimate* computeLostLepton( const MT2Sample& sample, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions ) {


  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  //definition of MT2tree
  MT2tree* fMT2tree = new MT2tree();
  sample.tree->SetBranchAddress("MT2tree", &fMT2tree);


  std::ostringstream preselectionStream;
  preselectionStream << " " 
    << "NTausIDLoose3Hits==0"                   << "&&"
    << "misc.Jet0Pass ==1"                      << "&&"
    << "misc.Jet1Pass ==1"                      << "&&"
    << "misc.Vectorsumpt < 70"                  << "&&" 
    << "misc.MinMetJetDPhi4Pt40 >0.3"           << "&&"
    << "misc.MET>30.";
  if( fFast ) preselectionStream << "&&misc.MT2>=100.";//lowest border in MT2

  // add "base" cuts:
  preselectionStream << " && " 
    << "misc.PassJet40ID ==1"                             << "&&"
    << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"  << "&&" // signal samples (fastsim) do not have this filter
    << "misc.CSCTightHaloIDFlag == 0"                     << "&&"
    << "misc.trackingFailureFlag==0"                      << "&&"
    << "misc.eeBadScFlag==0"                              << "&&"
    << "misc.EcalDeadCellTriggerPrimitiveFlag==0"         << "&&"
    << "misc.TrackingManyStripClusFlag==0"                << "&&"
    << "misc.TrackingTooManyStripClusFlag==0"             << "&&"
    << "misc.TrackingLogErrorTooManyClustersFlag==0"      << "&&"
    << "misc.CrazyHCAL==0";
  preselectionStream << "&&misc.MET/misc.CaloMETRaw<=2.";//HO cut

  

  std::ostringstream triggerStream;
  //if( isMET ) {

  //  triggerStream << "( ( ("
  //    << "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
  //    << "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )"
  //    << "||("
  //    << "trigger.HLT_PFHT350_PFMET100_v3==1 || trigger.HLT_PFHT350_PFMET100_v4==1 || trigger.HLT_PFHT350_PFMET100_v5==1 || "
  //    << "trigger.HLT_PFHT350_PFMET100_v6==1 || trigger.HLT_PFHT350_PFMET100_v7==1 || trigger.HLT_PFNoPUHT350_PFMET100_v1==1 || "
  //    << "trigger.HLT_PFNoPUHT350_PFMET100_v3==1 || trigger.HLT_PFNoPUHT350_PFMET100_v4==1 ) ) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";

  //} else  {

    triggerStream << "( ("
      << "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || "
      << "trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
      << "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";

  //}


  TString preselection = preselectionStream.str().c_str();
  TString trigger = triggerStream.str().c_str();
  TString cuts = ( sample.type=="data") ? (preselection + " && " + trigger) : preselection;
  
  TString cuts_ele = cuts + " && NEles==1 && NMuons==0";
  TString cuts_muo = cuts + " && NEles==0 && NMuons==1";

  TTree* tree_ele = sample.tree->CopyTree(cuts_ele);
  TTree* tree_muo = sample.tree->CopyTree(cuts_muo);


  MT2LeptonTypeLLEstimate* llest_ele = getLeptonTypeLLEstimate( "Ele", tree_ele, sample, HTRegions, signalRegions );
  MT2LeptonTypeLLEstimate* llest_muo = getLeptonTypeLLEstimate( "Muo", tree_muo, sample, HTRegions, signalRegions );

  MT2LostLeptonEstimate* llest = new MT2LostLeptonEstimate(sample.name, sample.sname);
  llest->l["Ele"] = llest_ele;
  llest->l["Muo"] = llest_muo;

  return llest;

}




MT2LeptonTypeLLEstimate* getLeptonTypeLLEstimate( const std::string& leptType, TTree* tree, MT2Sample sample, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions ) {


  int pdgIdLept = -999;
  if(leptType=="Ele") {
    pdgIdLept = 11;
  } else if(leptType=="Muo") {
    pdgIdLept = 13;
  } else {
    std::cout << "Unkown lepton type: " << leptType << "! Exiting." << std::endl;
    exit(111);
  }

  std::cout << "  -> Starting lepton type: " << leptType << std::endl;


  // global sample weight:
  Double_t weight = sample.xsection * sample.kfact * sample.lumi / (sample.nevents*sample.PU_avg_weight);


  // return estimate:
  std::string fullName = leptType + "_" + sample.name;
  MT2LeptonTypeLLEstimate* llest = new MT2LeptonTypeLLEstimate( fullName, sample.sname, HTRegions, signalRegions );


  MT2tree* fMT2tree = new MT2tree();
  sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

  int nentries = tree->GetEntries();


  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 10000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

    tree->GetEntry(iEntry);

    float ht   = fMT2tree->misc.HT;
    float met  = fMT2tree->misc.MET;
    float mt2  = fMT2tree->misc.MT2;
    int njets  = fMT2tree->NJetsIDLoose40;
    int nbjets = fMT2tree->NBJets40CSVM;

    float mt = (leptType=="Ele") ? fMT2tree->ele[0].MT : fMT2tree->muo[0].MT;

    int ngenlept = fMT2tree->GenNumLeptFromW(pdgIdLept, 0, 1000, fIncludeTaus);


    for( unsigned iHT=0; iHT<HTRegions.size(); ++iHT ) {


      if( ht<HTRegions[iHT].htMin || ht>=HTRegions[iHT].htMax || met<HTRegions[iHT].metMin ) continue;

      for( unsigned iSR=0; iSR<signalRegions.size(); ++iSR ) {

        if( njets <signalRegions[iSR].nJetsMin   || njets >signalRegions[iSR].nJetsMax ) continue;
        if( nbjets<signalRegions[iSR].nBJetsMin  || nbjets>signalRegions[iSR].nBJetsMax ) continue;

        MT2Region thisRegion( &HTRegions[iHT], &signalRegions[iSR] );

        MT2SingleLLEstimate* thisRecoEst = llest->getRecoRegion( thisRegion.getName() );

        thisRecoEst->yield->Fill( mt2, weight );

        thisRecoEst->effMT_tot->Fill( mt2, weight );
        if( mt<100. ) thisRecoEst->effMT_pass->Fill( mt2, weight );

        thisRecoEst->effLept_tot->Fill( mt2, weight );
        if( ngenlept>0 ) thisRecoEst->effLept_pass->Fill( mt2, weight );

      }  // for iSR

    }  // for iHT

  } // for entries


  std::cout << "  -> Finished lepton type: " << leptType << std::endl;
  std::cout << std::endl << std::endl;
  
  return llest;

}





MT2LostLeptonEstimate* mergeEstimates( std::vector<MT2LostLeptonEstimate*> llest, const std::string& n1, const std::string& n2, const std::string& n3, const std::string& n4, const std::string& n5 ) {

  std::vector<std::string> snames;
  snames.push_back(n1);
  snames.push_back(n2);
  snames.push_back(n3);
  snames.push_back(n4);
  snames.push_back(n5);

  std::string newname = (snames.size()>0) ? (snames[0]) : "";
  for( unsigned i=1; i<snames.size(); ++i ) newname += "_" + snames[i];
  newname = "merge_" + newname;

  MT2LostLeptonEstimate* return_llest = new MT2LostLeptonEstimate(newname);

  for( unsigned i=0; i<llest.size(); ++i ) {

    for( unsigned iname=0; iname<snames.size(); ++iname ) {

      if( llest[i]->SName == snames[iname] ) {
        return_llest->add(*(llest[i]));
      }

    } // for snames

  } // for llest


  return return_llest;

}




std::vector<TH1D*> createPredictionHistos( const std::string& prefix, const std::string& leptType, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions) {

  int nHistos = HTRegions.size();
  int nBins = signalRegions.size();

  
  std::vector<TH1D*> histos;

  for( unsigned i=0; i<nHistos; ++i ) {

    TH1D* h1 = new TH1D(Form("%s_%s_%s", prefix.c_str(), leptType.c_str(), HTRegions[i].name.c_str()), "", nBins, 0., nBins);
    h1->Sumw2();

    for( unsigned iBin=1; iBin<(nBins+1); ++iBin ) 
      h1->GetXaxis()->SetBinLabel(iBin, signalRegions[iBin].getName().c_str());

    histos.push_back(h1);

  } // for

  return histos;

}





    
/*

void computeSignalRegionEstimate( TFile* outfile, MT2Sample sample, int njet_min, int njet_max, int nbjet_min, int nbjet_max, float ht_min, float ht_max ) {

  bool isMET = (ht_min==450.);

  std::string signal_region = MT2Common::getSignalRegion(ht_min, njet_min, njet_max, nbjet_min, nbjet_max);

  int nBins;
  Double_t *bins;
  MT2Common::getBins(signal_region, nBins, bins);


  std::string hs = "_" + leptType + "_" + signal_region + "_" + sample.type;

  std::string histoName = "RecoLepEvents" + hs;
  TH1D* histo = new TH1D( histoName.c_str(), "", nBins, bins );
  histo->Sumw2();




  std::ostringstream cutStream;
  cutStream << " " 
    << "NTausIDLoose3Hits==0"                   << "&&"
    << "misc.Jet0Pass ==1"                      << "&&"
    << "misc.Jet1Pass ==1"                      << "&&"
    << "misc.Vectorsumpt < 70"                  << "&&" 
    << "misc.MinMetJetDPhi4Pt40 >0.3"           << "&&"
    << "misc.HT>" << ht_min << "&& misc.HT<=" << ht_max;
  if( isMET ) cutStream << "&&misc.MET>200.";
  else        cutStream << "&&misc.MET>30.";
  if( fFast ) cutStream << "&&misc.MT2>=100.";//lowest border in MT2
  
  std::ostringstream cutStreamBase;
  cutStreamBase << " " 
    << "misc.PassJet40ID ==1"                             << "&&"
    << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"  << "&&" // signal samples (fastsim) do not have this filter
    << "misc.CSCTightHaloIDFlag == 0"                     << "&&"
    << "misc.trackingFailureFlag==0"                      << "&&"
    << "misc.eeBadScFlag==0"                              << "&&"
    << "misc.EcalDeadCellTriggerPrimitiveFlag==0"         << "&&"
    << "misc.TrackingManyStripClusFlag==0"                << "&&"
    << "misc.TrackingTooManyStripClusFlag==0"             << "&&"
    << "misc.TrackingLogErrorTooManyClustersFlag==0"      << "&&"
    << "misc.CrazyHCAL==0";
  cutStreamBase << "&&misc.MET/misc.CaloMETRaw<=2.";//HO cut
  

  std::ostringstream triggerStream;
  if( isMET ) {

    triggerStream << "( ( ("
      << "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
      << "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )"
      << "||("
      << "trigger.HLT_PFHT350_PFMET100_v3==1 || trigger.HLT_PFHT350_PFMET100_v4==1 || trigger.HLT_PFHT350_PFMET100_v5==1 || "
      << "trigger.HLT_PFHT350_PFMET100_v6==1 || trigger.HLT_PFHT350_PFMET100_v7==1 || trigger.HLT_PFNoPUHT350_PFMET100_v1==1 || "
      << "trigger.HLT_PFNoPUHT350_PFMET100_v3==1 || trigger.HLT_PFNoPUHT350_PFMET100_v4==1 ) ) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";

  } else  {

    triggerStream << "( ("
      << "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || "
      << "trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
      << "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";

  }

  TString cuts = cutStream.str().c_str();
  TString basecuts = cutStreamBase.str().c_str();
  TString trigger = triggerStream.str().c_str();


  Long64_t nbytes = 0, nb = 0;
  int nev =0;

  TString myCuts = cuts + "&&" + basecuts;
  if( sample.type=="data") myCuts += " && " + trigger; //cuts to be aplied only on data
 
  std::cout << "Cuts for Flow: " << myCuts << std::endl;
  TTree* tree_sel = sample.tree->CopyTree(myCuts);
 
  Long64_t nentries = tree_sel->GetEntries();

  std::cout << "Filtering done, size=" << nentries  << std::endl;

  //run only over events passing event selection
  for( unsigned jentry=0; jentry<nentries; ++jentry ) {

    tree_sel->GetEntry(jentry); 

    if ( jentry % 100000 == 0 )  std::cout << "+++ Proccessing event " << jentry << " / " << nentries << std::endl;

    Double_t weight = sample_weight;
    bool isData = fMT2tree->misc.isData;

    if ( !(isData) ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 

    float ht = fMT2tree->misc.HT;
    float mt2 = fMT2tree->misc.MT2;

    std::string sampletype = MT2Common::getSampleType(sample);
    if(ht>=750. && ht<1200. && fMT2tree->NJetsIDLoose40==2 && fMT2tree->NBJets40CSVM>=1 && mt2 >=135.&& mt2 <=170. && sampletype=="SingleTop" && fMT2tree->misc.LumiSection==149 && fMT2tree->misc.Event==44699) continue;//sample T_t_highHT - fix for one 8 TeV MC sample with crazy weight

    //the if clause below computes the weight due to the 'ISR recipe'
    Double_t ISRweight(1.); 
    Double_t weightnoISR = weight;
    if(fISRreweight && !isData) ISRweight = getISRCorrection(fMT2tree, sample, sampletype);
  
  
    int nele = fMT2tree->NEles;
    int nmu  = fMT2tree->NMuons;
    Bool_t recoedele       = (nele ==1 && nmu ==0 && fMT2tree->ele[0].MT<100. );
    Bool_t recoedmuo       = (nmu  ==1 && nele==0 && fMT2tree->muo[0].MT<100. );
    Bool_t recoedelenomt   = (nele ==1 && nmu ==0                             );
    Bool_t recoedmuonomt   = (nmu  ==1 && nele==0                             );
    Bool_t norecolep       = (nmu  ==0 && nele==0                             );

    std::string slep;
    if(recoedele) slep = "_Ele";
    if(recoedmuo) slep = "_Muo";

    std::string sHT;
    if( isMET ){
      if( ht<=450. )      sHT = "_HTge0";
      else if( ht<=750. ) sHT = "_HTge450";
    } else {
      if( ht>1200. )      sHT = "_HTge1200";
      else if( ht>750. )  sHT = "_HTge750";
    }


    //get BTV weight
    double btagSF(1.), btagSFerr(0.);
    if(!isData && fbTagReweight)
      getBTagSF( fMT2tree->NJetsIDLoose40, fMT2tree->NBJets40CSVM, btagSF, btagSFerr );
     
    weight      *= btagSF;
    weightnoISR *= btagSF;

    std::string hh = slep + sHT + "_" + signal_region + "_" + sampletype;
    std::string hhnolep  =  sHT + "_" + signal_region + "_" + sampletype;//needed for genlept filling

    if(recoedele||recoedmuo) histos[(string)"RecoLepEvents"        + hh]->Fill(fMT2tree->misc.MT2, weight);//no bg, lepton
    if(recoedele||recoedmuo) histos[(string)"noISR_RecoLepEvents"        + hh]->Fill(fMT2tree->misc.MT2, weightnoISR);//no bg, lepton
    if(fbTagError && !fMT2tree->misc.isData && btagSF!=0){//data has no BTV SF weights
    if(recoedele||recoedmuo) histos[(string)"RecoLepEventsUp"         + hh]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF+btagSFerr));
    if(recoedele||recoedmuo) histos[(string)"RecoLepEventsDown"       + hh]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF-btagSFerr));
    if(recoedele||recoedmuo) histos[(string)"noISR_RecoLepEventsUp"   + hh]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF+btagSFerr));
    if(recoedele||recoedmuo) histos[(string)"noISR_RecoLepEventsDown" + hh]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF-btagSFerr));
    }
    
    if(recoedele||recoedmuo) histos[(string)"LeptonEvents"               + hh     ]->Fill(fMT2tree->misc.MT2,weight);//no bg, lepton
    if(recoedelenomt)        histos[(string)"LeptonEventsNoMT_Ele"       + hhnolep]->Fill(fMT2tree->misc.MT2,weight);//no bg, lepton <-- here I need hhnolep
    if(recoedmuonomt)        histos[(string)"LeptonEventsNoMT_Muo"       + hhnolep]->Fill(fMT2tree->misc.MT2,weight);//no bg, lepton <-- here I need hhnolep
    if(recoedele||recoedmuo) histos[(string)"noISR_LeptonEvents"         + hh     ]->Fill(fMT2tree->misc.MT2,weightnoISR);//no bg, lepton
    if(recoedelenomt)        histos[(string)"noISR_LeptonEventsNoMT_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2,weightnoISR);//no bg, lepton <-- here I need hhnolep
    if(recoedmuonomt)        histos[(string)"noISR_LeptonEventsNoMT_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2,weightnoISR);//no bg, lepton <-- here I need hhnolep

  } // for entries

} 







float getISRCorrection( MT2tree* fMT2tree, const MT2Sample& sample, const std::string& sampletype ) {

  TLorentzVector hardgenlv; hardgenlv.SetPtEtaPhiM(0,0,0,0);
  if(sampletype=="WJets"){
    bool foundW(false);
    for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){//lepton from W
      int ID  =abs(fMT2tree->genlept[ngl].ID);
      int MID =abs(fMT2tree->genlept[ngl].MID);
      if((ID == 11 || ID == 12 || ID == 13 || ID == 14 || ID == 15 || ID == 16) && MID==24){
        hardgenlv = fMT2tree->genlept[ngl].Mlv; foundW = true;
      }
      if(foundW) break;
    }
    if(!foundW){
      for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){//lepton from tau, tau from W
        int ID  =abs(fMT2tree->genlept[ngl].ID);
        int MID =abs(fMT2tree->genlept[ngl].MID);
        int GMID=abs(fMT2tree->genlept[ngl].GMID);
        if((ID == 11 || ID == 12 || ID == 13 || ID == 14 || ID == 16) && MID==15 && GMID==24){
          hardgenlv = fMT2tree->genlept[ngl].Mlv; foundW = true;
        }
        if(foundW) break;
      }
    }
  } else if(sampletype=="ZJets"){
    hardgenlv = fMT2tree->GenZ[0];
  } else if(sampletype=="TTbar"){
    TLorentzVector top1(0.,0.,0.,0.), top2(0.,0.,0.,0.);
    bool top1f(false), top2f(false);
    for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
      int id   = abs(fMT2tree->genlept[ngl].ID);
      if(id!=5) continue;
      int mid  = fMT2tree->genlept[ngl].MID;//from b
      if(mid==6&&top1f) continue;
      else if(mid==6) { top1 = fMT2tree->genlept[ngl].Mlv; top1f = true; }
      if(mid==-6&&top2f) continue;
      else if(mid==-6){ top2 = fMT2tree->genlept[ngl].Mlv; top2f = true; }
      if(top1f&&top2f) {
        hardgenlv = top1+top2;
        break;
      }
    }
  } else if(sampletype=="SingleTop"){
    TString samplename_tstr(sample.name);
    if(samplename_tstr.Contains("tW")){//t + W
      TLorentzVector top(0.,0.,0.,0.), W(0.,0.,0.,0.);
      bool topf(false), Wf(false);
      for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
        int id    = abs(fMT2tree->genlept[ngl].ID);
        int mid   = abs(fMT2tree->genlept[ngl].MID);//from b
        int gmid  = abs(fMT2tree->genlept[ngl].GMID);//from b
        if(mid==6&&topf) continue;
        else if(mid==6) { top = fMT2tree->genlept[ngl].Mlv; topf = true; }
        if(mid==24&&gmid!=6&&Wf) continue;
        if((id == 11 || id == 12 || id == 13 || id == 14 || id == 15 || id == 16) && mid==24 && gmid!=6){
          W = fMT2tree->genlept[ngl].Mlv; Wf = true;
        }
        if(topf&&Wf){
          hardgenlv = top+W;
          break;
        }
      }
      if(!Wf){//this might be wrong - but influence negligible anyway
        for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
          int id    = abs(fMT2tree->genlept[ngl].ID);
          int mid   = abs(fMT2tree->genlept[ngl].MID);//from b
          int gmid  = abs(fMT2tree->genlept[ngl].GMID);//from b
          if(mid==6||gmid==6) continue;
          if(gmid==24&&gmid==15&&Wf) continue;
          if((id == 11 || id == 12 || id == 13 || id == 14 || id == 15 || id == 16) && mid==15 && gmid==24){
            W = fMT2tree->genlept[ngl].Mlv; Wf = true;
          }
          if(topf&&Wf){
            hardgenlv = top+W;
            break;
          }
        }
      }
    
    } else {
      TLorentzVector top(0.,0.,0.,0.);
      bool topf(false), Wf(false);
      for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
        int id    = abs(fMT2tree->genlept[ngl].ID);
        int mid   = abs(fMT2tree->genlept[ngl].MID);//from b
        int gmid  = abs(fMT2tree->genlept[ngl].GMID);//from b
        if(mid==6&&topf) continue;
        else if(mid==6) { top = fMT2tree->genlept[ngl].Mlv; topf = true; }
        if(topf){
          hardgenlv = top;
          break;
        }
      }
    }
  }

  float ISRweight=1.;
  if(hardgenlv.Pt()>250.) ISRweight = 0.8;
  else if(hardgenlv.Pt()>150.) ISRweight = 0.9;
  else if(hardgenlv.Pt()>120.) ISRweight = 0.95;
  else                         ISRweight = 1.;
      
  return ISRweight;

}



void getBTagSF( int njets, int nbjets, double& btagSF, double& btagSFerr ) {


  if( njets < 2 ) { // shouldnt be possible
    btagSF=0.;
    btagSFerr=0.;
    return;
  }



  if( nbjets >= 3 ) {

    btagSF = fMT2tree->SFWeight.BTagCSV40ge3; 
    btagSFerr = fMT2tree->SFWeight.BTagCSV40ge3Error; 

  } else {

    if( njets == 2 ) {  
      if( nbjets == 0) { 
        btagSF = fMT2tree->SFWeight.BTagCSV40eq0; 
        btagSFerr = fMT2tree->SFWeight.BTagCSV40eq0Error; 
      } else {
        btagSF = fMT2tree->SFWeight.BTagCSV40ge1; 
        btagSFerr = fMT2tree->SFWeight.BTagCSV40ge1Error; 
      }
    } else if( njets >= 3 && njets <= 5 ) {
      if( nbjets == 0) { 
        btagSF = fMT2tree->SFWeight.BTagCSV40eq0; 
        btagSFerr = fMT2tree->SFWeight.BTagCSV40eq0Error;
      } else if( nbjets == 1) { 
        btagSF = fMT2tree->SFWeight.BTagCSV40eq1; 
        btagSFerr = fMT2tree->SFWeight.BTagCSV40eq1Error;
      } else { // this is only nbjets==2 (see above for nbjets>=3)
        btagSF = fMT2tree->SFWeight.BTagCSV40eq2; 
        btagSFerr = fMT2tree->SFWeight.BTagCSV40eq2Error;
      }
    } else { // this is njets >=6 (but still nbjets < 3 )
      if( nbjets == 0) { 
        btagSF = fMT2tree->SFWeight.BTagCSV40eq0; 
        btagSFerr = fMT2tree->SFWeight.BTagCSV40eq0Error;
      } else if( nbjets == 1) { 
        btagSF = fMT2tree->SFWeight.BTagCSV40eq1; 
        btagSFerr = fMT2tree->SFWeight.BTagCSV40eq1Error; 
      } else { // this is only nbjets==2 (see above for nbjets>=3)
        btagSF = fMT2tree->SFWeight.BTagCSV40eq2; 
        btagSFerr = fMT2tree->SFWeight.BTagCSV40eq2Error; 
      }
    }

  }

}

*/
