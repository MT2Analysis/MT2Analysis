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
MT2LeptonTypeLLEstimate* getLeptonTypeLLEstimate( const std::string& leptType, TTree* tree, MT2Sample sample, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions );
//float getISRCorrection( MT2tree* fMT2tree, const MT2Sample& sample, const std::string& sampletype );
//void getBTagSF( int njets, int nbjets, double& btagSF, double& btagSFerr );
MT2LostLeptonEstimate* mergeEstimates( std::vector<MT2LostLeptonEstimate*> llest, const std::string& n1, const std::string& n2="", const std::string& n3="", const std::string& n4="", const std::string& n5="" );
std::vector<TH1D*> getPredictionHistos( const std::string& prefix, const std::string& leptType, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions, MT2LostLeptonEstimate* ll_tot, MT2LostLeptonEstimate* ll_bg, MT2LostLeptonEstimate* ll_eff );




int main( int argc, char* argv[] ) {


  if( argc!=2 ) {
    std::cout << "USAGE: ./computeLostLepton [samplesFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }


  std::string sampleName(argv[1]);

  std::string samplesFileName = "samples/samples_" + sampleName + ".dat";
  std::cout << std::endl << std::endl;
  std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;

  std::vector<MT2Sample> fSamples = MT2Common::loadSamples(samplesFileName);

 
  std::vector<MT2HTRegion> HTRegions;
  HTRegions.push_back(MT2HTRegion("lowHT",    450.,    750., 200.));
  HTRegions.push_back(MT2HTRegion("mediumHT", 750.,   1200.,  30.));
  HTRegions.push_back(MT2HTRegion("highHT",  1200., 100000.,  30.));

  std::vector<MT2SignalRegion> signalRegions;
  signalRegions.push_back(MT2SignalRegion(2, 2, 0, 0));  // 2j0b
  signalRegions.push_back(MT2SignalRegion(2, 2, 1, 2));  // 2j1to2b
  signalRegions.push_back(MT2SignalRegion(-1, -1, 3, -1));  // 3b
  signalRegions.push_back(MT2SignalRegion(3, 5, 0, 0));  // 3to5j0b
  signalRegions.push_back(MT2SignalRegion(3, 5, 1, 1));  // 3to5j1b
  signalRegions.push_back(MT2SignalRegion(3, 5, 2, 2));  // 3to5j2b
  signalRegions.push_back(MT2SignalRegion(6, 6, 0, 0));  // 6j0b
  signalRegions.push_back(MT2SignalRegion(6, 6, 1, 1));  // 6j1b
  signalRegions.push_back(MT2SignalRegion(6, 6, 2, 2));  // 6j2b



  
  std::vector<MT2LostLeptonEstimate*> llest;
  for( unsigned i=0; i<fSamples.size(); ++i )
    llest.push_back( computeLostLepton( fSamples[i], HTRegions, signalRegions ) );
  

  std::cout << "-> Done looping on samples. Start merging." << std::endl;

  MT2LostLeptonEstimate* ll_data  = mergeEstimates( llest, "HT-Data" );
  MT2LostLeptonEstimate* ll_top   = mergeEstimates( llest, "Top" );
  MT2LostLeptonEstimate* ll_qcd   = mergeEstimates( llest, "QCD" );
  MT2LostLeptonEstimate* ll_wjets = mergeEstimates( llest, "Wtolnu" );
  MT2LostLeptonEstimate* ll_other = mergeEstimates( llest, "DY", "VV" );
  
  MT2LostLeptonEstimate* ll_topW   = new MT2LostLeptonEstimate(*ll_top  + *ll_wjets);
  MT2LostLeptonEstimate* ll_bg     = new MT2LostLeptonEstimate(*ll_qcd  + *ll_other);
  MT2LostLeptonEstimate* ll_allMC  = new MT2LostLeptonEstimate(*ll_topW + *ll_bg);

  std::cout << "-> Done merging. Start computing preditions." << std::endl;


  std::vector<TH1D*> vh1_data_ele = getPredictionHistos( "Prediction", "Ele", HTRegions, signalRegions, ll_data, ll_bg, ll_topW );
  std::vector<TH1D*> vh1_data_mu  = getPredictionHistos( "Prediction", "Muo", HTRegions, signalRegions, ll_data, ll_bg, ll_topW );

  std::vector<TH1D*> vh1_mc_ele   = getPredictionHistos( "SimulationTruth", "Ele", HTRegions, signalRegions, ll_allMC, ll_bg, ll_topW );
  std::vector<TH1D*> vh1_mc_mu    = getPredictionHistos( "SimulationTruth", "Muo", HTRegions, signalRegions, ll_allMC, ll_bg, ll_topW );



  std::string outputdir = "LostLeptonFiles";
  system(Form("mkdir -p %s", outputdir.c_str()));

  
  TFile* outfile = TFile::Open(Form("%s/ll_%s.root", outputdir.c_str(), sampleName.c_str()), "recreate");
  outfile->cd();

  for( unsigned i=0; i<HTRegions.size(); ++i ) {

    vh1_data_ele[i]->Write();
    vh1_data_mu [i]->Write();

    vh1_mc_ele  [i]->Write();
    vh1_mc_mu   [i]->Write();

  }


  outfile->Close();

  std::cout << "-> Done computing preditions. Written to file: " << outfile->GetName() << std::endl;

  return 0;

}




MT2LostLeptonEstimate* computeLostLepton( const MT2Sample& sample, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions ) {


  std::cout << std::endl << std::endl;
  std::cout << "-> Starting computation for sample: " << sample.name << std::endl;

  TFile* file = TFile::Open(sample.file.c_str());
  TTree* tree = (TTree*)file->Get("MassTree");



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

  gROOT->cd();

  TTree* tree_ele = tree->CopyTree(cuts_ele);
  TTree* tree_muo = tree->CopyTree(cuts_muo);


  MT2LeptonTypeLLEstimate* llest_ele = getLeptonTypeLLEstimate( "Ele", tree_ele, sample, HTRegions, signalRegions );
  MT2LeptonTypeLLEstimate* llest_muo = getLeptonTypeLLEstimate( "Muo", tree_muo, sample, HTRegions, signalRegions );

  MT2LostLeptonEstimate* llest = new MT2LostLeptonEstimate(sample.name, sample.sname);
  llest->l["Ele"] = llest_ele;
  llest->l["Muo"] = llest_muo;

  delete tree;
  delete tree_ele;
  delete tree_muo;
  file->Close();

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
  tree->SetBranchAddress("MT2tree", &fMT2tree);

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

        if( signalRegions[iSR].nJetsMin  >= 0  &&  njets  < signalRegions[iSR].nJetsMin ) continue;
        if( signalRegions[iSR].nJetsMax  >= 0  &&  njets  > signalRegions[iSR].nJetsMax ) continue;
        if( signalRegions[iSR].nBJetsMin >= 0  &&  nbjets < signalRegions[iSR].nBJetsMin ) continue;
        if( signalRegions[iSR].nBJetsMax >= 0  &&  nbjets > signalRegions[iSR].nBJetsMax ) continue;

        MT2Region thisRegion( &HTRegions[iHT], &signalRegions[iSR] );

        MT2SingleLLEstimate* thisRecoEst = llest->getRegion( thisRegion.getName() );

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








std::vector<TH1D*> getPredictionHistos( const std::string& prefix, const std::string& leptType, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions, MT2LostLeptonEstimate* ll_tot, MT2LostLeptonEstimate* ll_bg, MT2LostLeptonEstimate* ll_eff ) {


  int nHistos = HTRegions.size();
  int nBins = signalRegions.size();

  
  std::vector<TH1D*> histos;

  for( unsigned i=0; i<nHistos; ++i ) {

    TH1D* h1 = new TH1D(Form("%s_%s_%s", prefix.c_str(), leptType.c_str(), HTRegions[i].name.c_str()), "", nBins, 0., nBins);
    h1->Sumw2();

    for( unsigned j=0; j<nBins; ++j ) {
  
      int iBin = j+1;

      MT2Region thisRegion( &HTRegions[i], &signalRegions[j] );

      h1->GetXaxis()->SetBinLabel(iBin, signalRegions[j].getName().c_str());

      float totYield = (ll_tot->l[leptType.c_str()]!=0) ? ll_tot->l[leptType.c_str()]->getRegion(thisRegion.getName())->yield->Integral() : 0.;
      float  bgYield = (ll_bg ->l[leptType.c_str()]!=0) ? ll_bg ->l[leptType.c_str()]->getRegion(thisRegion.getName())->yield->Integral() : 0.;

      float  effMT   = (ll_eff->l[leptType.c_str()]!=0) ? ll_eff->l[leptType.c_str()]->getRegion(thisRegion.getName())->effMT()  ->GetEfficiency(1) : 0.;
      float  effLept = (ll_eff->l[leptType.c_str()]!=0) ? ll_eff->l[leptType.c_str()]->getRegion(thisRegion.getName())->effLept()->GetEfficiency(1) : 0.;

      float pred = (effLept>0. && effMT>0.) ? (totYield-bgYield)*(1.-effLept)/(effLept*effMT) : 0.;

      h1->SetBinContent( iBin, pred );

    }

    histos.push_back(h1);

  } // for

  return histos;

}








/*

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
