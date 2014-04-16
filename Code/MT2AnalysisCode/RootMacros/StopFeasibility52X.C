#include "TEventList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TCut.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TMap.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLine.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string>
#include <cmath>
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path of MT2tree.hh
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/MT2Code/include/MassPlotter.hh"//use your own path of MassPlotter.hh
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/MT2Code/include/BTagWeight.hh"
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/Stops/SolveTTbar.h"
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/Stops/SolveTTbar.c"
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/Stops/SolveTTbarNew.h"
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/Stops/SolveTTbarNew.C"

//use via root -l -b -q StopFeasibility52X.C++
//this code is superseeded by StopFeasibility.C (uses 53X)
//therefore no comments are done here
//only kept if I need it for my thesis

using namespace std;

void load(const char* filename = "/shome/haweber/CMSSW_4_2_3/src/DiLeptonAnalysis/NTupleProducer/MT2Analysis/Code/MT2AnalysisCode/MT2Code/samples/datasamples/samples_2141_dataonly.dat");
void StopFeasibility52X();
void MakeCorrelationPlots(map<string, TH2D*> histos2, const int leptontypesize, string *lepton_type, vector<string> histonames, TString outputdirectory);
void MakePlots(MassPlotter *tA, map<string, TH1D*> histos, const int leptontypesize, string *lepton_type, map<string, THStack*> stacks, TLegend *Legend1, vector<string> histonames);
void Make3DPlots(map<string, TH3D*> histos3, const int leptontypesize, string *lepton_type, vector<string> histonames, TString outputdirectory);

struct sample{
	TString name;
	TString sname;
	TString shapename;
	TString type;
	TFile *file;
	TTree *tree;
	float xsection;
	float nevents;
	float kfact;
	float PU_avg_weight;
	float lumi;
	int color;
};

vector < sample >  fSamples;
const int fVerbose = 3; 
TString fPath;

TString btagging_file = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/Efficiencies/Stop_BEfficiencies_PTbinned_nocuts_SSVHEM.root";
TFile* btagefffile;
TString tagger = "SSVHEM";// SSVHPT or SSVHEM
TH1D*  hbeff;
TH1D*  hceff;
TH1D*  hleff;
int njets = -2;
int ntags = 2;//- means >=, while + means ==; -99 means no tagging requirement
int Tagger = 2;//TCHE:0 TCHP:1 SSVHE:2 SSVHP:2
float discr = 1.74;
float bpt = 40.;

TString samples                   = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_Stops_presel.dat";//dummy
Bool_t  runData                   = true; // if you do only MCClosure you can set this to false
Bool_t  calcsusy                  = true; // not used - can set false
TString outputdir                 = "Stops/testeq2bge2jStop400LSP75_20120614";//default
TString outputname                = "DiLeptonicStops.root";//default
Bool_t  dofastbtagSFreweighting   = false;
Bool_t  logflag                   = true;
Bool_t  fSave                     = true;
Bool_t  debug                     = false;//don't run over mc to save time
Bool_t  calcD2          = false;
Bool_t  debugD2     = false;
Bool_t  calcD2Luc	  = false;
Bool_t  calcD2Minuit      = false;
Bool_t  plotonlywithratio         = true;
Bool_t  plotonlyLL                = true;
Bool_t  stupidchecks              = false;
//tail plots
Bool_t  tailMT2bb                 = false;
Bool_t  tailMT2bb_linMET          = false;
Bool_t  tailMT2ll                 = false;
Bool_t  tailMT2lb_massive         = false;
Bool_t  tailMT2bb_massless        = false;


void StopFeasibility52X(){

	if(stupidchecks) outputdir = "Stops/stupidchecks_ge2btags";

gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");//Needed
  tA->SetSave(true);
  tA->setVerbose(fVerbose);
  tA->init(samples);
  tA->SetIsPhoton(false);

  btagefffile	= TFile::Open(btagging_file);
  hbeff		= (TH1D*)btagefffile->Get("h_beff_NJets40ge2");//also NJets40eq2
  hceff		= (TH1D*)btagefffile->Get("h_ceff_NJets40ge2");
  hleff		= (TH1D*)btagefffile->Get("h_leff_NJets40ge2");

    Util::MakeOutputDir(outputdir);
    map<string, TH1D*> histos;
    map<string, TH2D*> histos2;
    map<string, TH3D*> histos3;
    map<string, THStack*> hstacks;
    map<string, THStack*> stacks;
    TLegend *Legend1 = new TLegend(.71,.54,.91,.92);
    vector<string> histonames; histonames.clear();
    TLegend* legend = new TLegend(.71,.54,.91,.92);
    legend -> SetFillColor(0);
    legend -> SetBorderSize(0);

    std::ostringstream baseStream;
    std::ostringstream cutStream;
    TString basecuts;
if(!stupidchecks){
	baseStream << " " 
    << "misc.MET>=30"                                                  << "&&"//minimal MET for MT2
    << "misc.PassJetID ==1"                                            << "&&"
    << "NJetsIDLoose40>=2"                                             << "&&"
    << "(Sum$(ele.lv.Pt()>10) + Sum$(muo.lv.Pt()>10))>=2"              << "&&"//veto secondary leptons
    // Noise
    << "misc.CSCTightHaloID==0"                                        << "&&"
    << "misc.CrazyHCAL==0";
	basecuts = baseStream.str().c_str();
	cutStream << " " 
    << "NJetsIDLoose40>=2"                                             << "&&"
    << "GetNBtags(2,1.74,40.,2.4,1)==2"                                << "&&"//change this also below for NumBJets
    << "(Sum$(ele.lv.Pt()>10) + Sum$(muo.lv.Pt()>10))==2";
}
else{
	baseStream << " " 
    << "misc.MET>=10"                                                  << "&&"//minimal MET for MT2
    << "misc.PassJetID ==1"                                            << "&&"
    << "NJetsIDLoose>=2"                                               << "&&"
    << "(Sum$(ele.lv.Pt()>8) + Sum$(muo.lv.Pt()>8))>=2"              << "&&"//veto secondary leptons
    // Noise
    << "misc.CSCTightHaloID==0"                                        << "&&"
    << "misc.CrazyHCAL==0";
	basecuts = baseStream.str().c_str();
	//Analysis cuts
	cutStream << " " 
    << "GetNBtags(2,1.74,20.,2.4,1)>=2";
}
	TString cuts = cutStream.str().c_str();
	cuts = basecuts + (string)"&&" + cuts;

	std::ostringstream triggerStreamEE;
	triggerStreamEE << "( "
	<< "(trigger.HLT_DiElectrons==1)" << " )";
	TString triggerEE = triggerStreamEE.str().c_str();
	std::ostringstream triggerStreamEMu;
	triggerStreamEMu << "( "
	<< "(trigger.HLT_EMu==1)" << " )";
	TString triggerEMu = triggerStreamEMu.str().c_str();
	std::ostringstream triggerStreamMuMu;
	triggerStreamMuMu << "( "
	<< "(trigger.HLT_DiMuons==1)" << " )";
	TString triggerMuMu = triggerStreamMuMu.str().c_str();


    const int sampletypesize = 10;
    string sample_type[sampletypesize] = {"WJets", "ZJets", "TTbar", "SingleTop", "TTbarV", "VV/VVV", "Other", "mc", "Stop", "data"};

    const int leptontypesize = 4;
    string lepton_type[leptontypesize] = {"MuMu", "EMu", "EE", "LL"};

	cout << " initializing histograms ...";

	vector<string> vs; vs.clear();
	bool av = true;//append vector
	for(int is = 0; is<sampletypesize; ++is){
           for(int il = 0; il<leptontypesize; ++il){
		string hs = string("_") + lepton_type[il] + string("_") + sample_type[is];
		string mapname;
		string title, xtitle, ytitle, ztitle;
		bool firstentry = true;
		if(il<leptontypesize-1) firstentry = false;
		mapname = "LepPt"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Lepton Pt"; xtitle = "Lepton p_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "LepMT"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Lepton MT"; xtitle = "1st Lepton m_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "B-Jet Pt - the ones which make MT2lbmin"; xtitle = "B p_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "AllBPt"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "B-Jet Pt - all B jets"; xtitle = "B p_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NJetsID20"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of Jets(Pt>20 GeV,passing Loose JID)"; xtitle = "#jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NJetsID30"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of Jets(Pt>30 GeV,passing Loose JID)"; xtitle = "#jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NJetsID40"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of Jets(Pt>40 GeV,passing Loose JID)"; xtitle = "#jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NJetsID50"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of Jets(Pt>50 GeV,passing Loose JID)"; xtitle = "#jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NBJets20"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of B-Jets(Pt>20 GeV)"; xtitle = "#b jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 5, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NBJets30"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of B-Jets(Pt>30 GeV)"; xtitle = "#b jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 5, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NBJets40"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of B-Jets(Pt>40 GeV)"; xtitle = "#b jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 5, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NBJets50"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of B-Jets(Pt>50 GeV)"; xtitle = "#b jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 5, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BtaggingDiscriminant"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "BTagging Discriminant"; xtitle = "b-tag discriminant";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BtaggingDiscriminant_MT2minbjets"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "BTagging Discriminant - jets fot MT2min"; xtitle = "b-tag discriminant";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "HT"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "HT"; xtitle = "H_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 85, 0, 1275);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MET"; xtitle = "E_{T}^{miss} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 1000);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2min"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2(min)"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2min_massless"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2(min)"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2min_withMlbcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2(min) after M(lb) cut"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2min_massless_withMlbcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2(min) after M(lb) cut"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_allcombi"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant Mass of lepton and B-Jets"; xtitle = "M_{lb} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_combiforMT2min"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant Mass of lepton B-Jets making minimal MT2"; xtitle = "M_{lb} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant DiLepton Mass"; xtitle = "M_{ll} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 105, 0, 1050);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mbb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant B-Jet Mass"; xtitle = "M_{bb} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 1050);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MbbAll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant B-Jet Mass - all combinations"; xtitle = "M_{bb} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 1050);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2(ll)"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_afterMlbcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2(ll) after M(lb) cut"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2(bb)"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb_lInMET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2(bb)"; xtitle = "M_{T2} (GeV), leps in MET";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb_lInMET_testmass80"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2(bb), l's in MET, testmass 80.4"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb_lInMET_testmass80_afterMlbcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2(bb), l's in MET, testmass 80.4 after M(lb) cut"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "JetsPt"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Jets Pt of subleading jets"; xtitle = " jet p_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "UTM"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Upstream Transverse Momentum"; xtitle = "UTM (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 500);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DPhi_lb_allcombi"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "DeltaPhi(Lepton, B-Jet)"; xtitle = "#Delta#phi(l,b)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DPhi_lb_combiforMT2min"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "DeltaPhi(Lepton, B-Jet) making minimal MT2"; xtitle = "#Delta#phi(l,b)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DPhi_ll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "DeltaPhi(leptons)"; xtitle = "#Delta#phi(l,l)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DPhi_bb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "DeltaPhi(B-Jets) for MT2min"; xtitle = "#Delta#phi(b,b)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DPhi_bbAll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "DeltaPhi(B-Jets), all combinations"; xtitle = "#Delta#phi(b,b)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DR_lb_allcombi"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "DeltaR(Lepton, B-Jet)"; xtitle = "#Delta R(l,b)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DR_lb_combiforMT2min"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "DeltaR(Lepton, B-Jet) making minimal MT2"; xtitle = "#Delta R(l,b)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DR_ll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "DeltaR(leptons)"; xtitle = "#Delta R(l,l)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DR_bb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "DeltaR(B-Jets) for MT2min"; xtitle = "#Delta R(b,b)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DR_bbAll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "DeltaR(B-Jets) all combinations"; xtitle = "#Delta R(b,b)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_genmatchingTop"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M(lb) genflavour-matching Top"; xtitle = "M_{lb} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_genmatchingTops"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2 - l1bl2b genflavour-matching tops"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2min_withMlb_NJets20ge2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2min(with Mlb) after NJets20>=2"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMlb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Mlb - genleveltruth"; xtitle = "M_{lb}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMT2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2 - genleveltruth"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMT2_withMlbcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2 after Mlb cut - genleveltruth"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMlb_noPtcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Mlb - genleveltruth_noPtcut"; xtitle = "M_{lb}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMT2_noPtcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2 - genleveltruth_noPtcut"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMT2_withMlbcut_noPtcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2 after Mlb cut - genleveltruth_noPtcut"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMtop_noPtcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "topMass - genleveltruth_noPtcut"; xtitle = "M_{bl#nu} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMtop"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "topMass - genleveltruth"; xtitle = "M_{bl#nu} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2min_withMlbcut_PUle5"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2min for #PV<=5 after Mlb cut"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2min_withMlbcut_PUge9"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2min for #PV>=9 after Mlb cut"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2min_withMlbcut_woZveto"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2min including Z Mll window after Mlbcut"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2min_woZveto"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2min including Z Mll window"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2min_massless_woZveto"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2min including Z Mll window"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb_woZveto"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2bb including Z Mll window"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_woZveto"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2ll including Z Mll window"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "No_PV"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of primary vertices"; xtitle = "#vertices";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 25);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2min_afterMlbcut_afterMT2llge90"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2min after Mlbcut and MT2(ll)>=90"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2min_afterMT2llge90"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2min after MT2(ll)>=90"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2min_massless_afterMT2llge90"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2minMassless after MT2(ll)>=90"; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "new Mass Fit Discriminant"; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_withMlb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "new Mass Fit Discriminant after Mlb cut"; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());

		//new 14/06/2012
		mapname = "MT2lb_bothbtags_truebjets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) (GeV)";        
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_onebtag_truebjet"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_nonebtags_truebjets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_bothbtags_matchedtobfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_onebtag_matchedtobfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_nonebtags_matchedtobfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "RecoPt_over_GenPartonPt_truebfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "Reco-p_{T} / genparton-p_{T}";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 40, 0, 4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "RecoPt_over_GenPartonPt_truebfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "Reco-p_{T} / genparton-p_{T}";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 40, 0, 4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_bothbtags_truebjets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_onebtag_truebjet"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_nonebtags_truebjets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_bothbtags_matchedtobfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_onebtag_matchedtobfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_nonebtags_matchedtobfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "RecoPt_over_GenPartonPt_truebfromtop_Mlb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "Reco-p_{T} / genparton-p_{T}";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 40, 0, 4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "RecoPt_over_GenPartonPt_truebfromtop_Mlb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "Reco-p_{T} / genparton-p_{T}";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 40, 0, 4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "p_{T} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "p_{T} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "TrueBfromtopEta_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "TrueBfromtopEta_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "TrueBfromtopPt_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "p_{T} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "TrueBfromtopPt_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "p_{T} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname  = "PFMEToverGenMET_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "genMET / PFMET";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 40, 0, 4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname  = "PFMEToverGenMET_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "genMET / PFMET";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 40, 0, 4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname  = "DeltaGenMET_PFMET_Phi_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta#phi(genMET,PFMET)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname  = "DeltaGenMET_PFMET_Phi_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta#phi(genMET,PFMET)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "AltMT2lb_MT2lb_gt_200_usingtruebpartons"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "AltMT2lb_MT2lb_gt_200_usingjetsmatchtotruebpartons"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaPt_TrueBparton_selectedjet_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta p_{T} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 350);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaPhi_TrueBparton_selectedjet_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta#phi";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaPt_bjet_matchedtrueparton_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta p_{T} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 350);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaEta_TrueBparton_selectedjet_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaPt_TrueBparton_selectedjet_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta p_{T}  (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 350);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaPt_bjet_matchedtrueparton_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta#phi";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 350);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaPhi_TrueBparton_selectedjet_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta#phi";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaEta_TrueBparton_selectedjet_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_bothleptons_true"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_onelepton_true"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_nolepton_true"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discriminant";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discriminant";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		//new in week 18/06/2012
		mapname = "DeltaPt_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta p_{T} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 350);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaPhi_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta#phi";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaEta_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaPt_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta p_{T} (GeV)";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 350);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaPhi_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta#phi";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaEta_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discriminant";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "B p_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discriminant";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "B p_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discriminant";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "B p_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discriminant";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "B p_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BJetID_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "jet ID";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 6, 0, 6);
		histos[mapname]->GetXaxis()->SetBinLabel(1, "none");   histos[mapname]->GetXaxis()->SetBinLabel(2, "loose"); histos[mapname]->GetXaxis()->SetBinLabel(3, "medium"); histos[mapname]->GetXaxis()->SetBinLabel(4, "tight");
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos[mapname]->GetXaxis()->SetLabelSize(0.03);
		mapname = "BJetID_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "jet ID";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 6, 0, 6);
		histos[mapname]->GetXaxis()->SetBinLabel(1, "none");   histos[mapname]->GetXaxis()->SetBinLabel(2, "loose"); histos[mapname]->GetXaxis()->SetBinLabel(3, "medium"); histos[mapname]->GetXaxis()->SetBinLabel(4, "tight");
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos[mapname]->GetXaxis()->SetLabelSize(0.03);
		//8 special histos: matching jet which was not selected as it is also matched to other bfromtop
		mapname = "BDiscr_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discriminant";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "B p_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discriminant";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "B p_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BJetID_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "jet ID";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 6, 0, 6);
		histos[mapname]->GetXaxis()->SetBinLabel(1, "none");   histos[mapname]->GetXaxis()->SetBinLabel(2, "loose"); histos[mapname]->GetXaxis()->SetBinLabel(3, "medium"); histos[mapname]->GetXaxis()->SetBinLabel(4, "tight");
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos[mapname]->GetXaxis()->SetLabelSize(0.03);
		mapname = "BJetID_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "jet ID";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 6, 0, 6);
		histos[mapname]->GetXaxis()->SetBinLabel(1, "none");   histos[mapname]->GetXaxis()->SetBinLabel(2, "loose"); histos[mapname]->GetXaxis()->SetBinLabel(3, "medium"); histos[mapname]->GetXaxis()->SetBinLabel(4, "tight");
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos[mapname]->GetXaxis()->SetLabelSize(0.03);
		//CHECKS1
		mapname = "X1_D2min_MT2llge100"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_D2min_MT2llge100_Mlble180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_D2min_MT2lbminge225"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_D2min_MT2lbminge225_Mlble180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_D2min_Mlble180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_D2min_MT2lbminge225_MT2llge100"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_D2min_MT2lbminge225_MT2llge100_Mlble180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_D2min"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2ll_D2minge50"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2ll_D2minge50_Mlble180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2ll_MT2lbminge225"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2ll_MT2lbminge225_Mlble180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2ll_Mlble180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2ll_D2minge50_MT2lbminge225"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2ll_D2minge50_MT2lbminge225_Mlble180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(ll) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2ll"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(ll) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2lbmin_MT2llge100"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2lbmin_MT2llge100_Mlble180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2lbmin_D2minge50"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2lbmin_D2minge50_Mlble180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2lbmin_Mlble180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2lbmin_D2minge50_MT2llge100"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2lbmin_D2minge50_MT2llge100_Mlble180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_Mlb_MT2lbminge225"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;; xtitle = "M_{lb} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_Mlb_D2minge50"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;; xtitle = "M_{lb} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_Mlb_MT2llge100"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;; xtitle = "M_{lb} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_Mlb_D2minge50_MT2lbminge225"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;; xtitle = "M_{lb} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_Mlb_MT2llge100_MT2lbminge225"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;; xtitle = "M_{lb} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_Mlb_D2minge50_MT2llge100"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;; xtitle = "M_{lb} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_Mlb_D2minge50_MT2llge100_MT2lbminge225"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;; xtitle = "M_{lb} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "X1_Mlb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;; xtitle = "M_{lb} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		//Check2
		mapname = "Y1_MT2lbcombi_Mlble180ForOneCombiOnly"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Y1_MT2lbmin_Mlble180ForOneCombiOnly"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Y1_MT2lbmin_Mlble180ForBothCombi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Y1_MT2lbmin_Mlbge180ForBothCombi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Y1_MT2lbmin_Mlbge180ForEverything"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		//Check3
		mapname = "Z1_LepPt_looser2ndLep"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "Lepton p_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_LepPt_looser1stLep"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "Lepton p_{T} (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_afterMlbcut_looser1stLep"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_looser1stLep"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_D2min_looser1stLep"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2ll_looser1stLep"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_afterMlbcut_looser2ndLep"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_looser2ndLep"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_D2min_looser2ndLep"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2ll_looser2ndLep"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_NBJets_tighterNBJets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#b jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 5, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_afterMlbcut_tighterNBJets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_tighterNBJets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_D2min_tighterNBJets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2ll_tighterNBJets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_BPt_looserBjet_30"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;; xtitle = "b jet p_{T}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_BPt_looserBjet_20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;; xtitle = "b jet p_{T}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 750);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_afterMlbcut_looserBJet_30"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_looserBJet_30"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_D2min_looserBJet_30"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2ll_looserBJet_30"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_afterMlbcut_looserBJet_20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_looserBJet_20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_D2min_looserBJet_20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2ll_looserBJet_20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MET_MET20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "E_T^{miss}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 1000);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MET_MET10"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "E_{T}^{miss}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 1000);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_afterMlbcut_MET10"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_MET10"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_D2min_MET10"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2ll_MET10"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_afterMlbcut_MET20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_MET20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_D2min_MET20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2ll_MET20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_afterMlbcut"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_D2min"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 100, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2ll"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_massless"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_massless_MET10"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_massless_MET20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_massless_looserBJet_20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_massless_looserBJet_30"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_massless_tighterNBJets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_massless_looser1stLep"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Z1_MT2lbmin_massless_looser2ndLep"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		//2d correlation checks1
		mapname = "X2_D2min_vs_MT2ll"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}"; ytitle = "M_{T2}(l,l) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 100, 0, 300, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "X2_D2min_vs_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}"; ytitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 100, 0, 300, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "X2_D2min_vs_MT2lbmin_afterMlb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D2"; ytitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 30, 0, 300, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "X2_D2min_vs_Mlb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}"; ytitle = "M_{lb} (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 100, 0, 300, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "X2_MT2ll_vs_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)"; ytitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "X2_MT2ll_vs_MT2lbmin_afterMlb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)"; ytitle = "M_{T2}(lb,lb) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "X2_MT2ll_vs_Mlb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l,l) (GeV)"; ytitle = "M_{lb} (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "X2_MT2lbmin_vs_Mlb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{lb} (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "X2_MT2lbmin_afterMlb_vs_Mlb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{lb} (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());

		mapname = "MT2min_vs_UTM"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massive MT2(min) vs UTM"; xtitle = "M_{T2} (GeV)"; ytitle = "UTM";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 50, 0, 500);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_withMlb_vs_UTM"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massive MT2(min) after M(l,b) cut vs UTM"; xtitle = "M_{T2} (GeV)"; ytitle = "UTM";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 50, 0, 500);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "DR_lb_configuration"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Closest DeltaR configuration";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 2, 0, 2, 2, 0, 2);
		histos2[mapname]->GetXaxis()->SetBinLabel(1, "l1 closer to b1"); histos2[mapname]->GetXaxis()->SetBinLabel(2, "l1 closer to b2"); 
		histos2[mapname]->GetYaxis()->SetBinLabel(1, "l2 closer to b1"); histos2[mapname]->GetYaxis()->SetBinLabel(2, "l2 closer to b2"); 
		histos2[mapname]->GetXaxis()->SetLabelSize(0.03); histos2[mapname]->GetYaxis()->SetLabelSize(0.03);
		mapname = "trueb1b2_configuration"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "match b1 b2 to true b";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 2, 0, 2, 2, 0, 2);
		histos2[mapname]->GetXaxis()->SetBinLabel(1, "b1"); histos2[mapname]->GetXaxis()->SetBinLabel(2, "b2"); 
		histos2[mapname]->GetYaxis()->SetBinLabel(1, "true"); histos2[mapname]->GetYaxis()->SetBinLabel(2, "fake"); 
		histos2[mapname]->GetXaxis()->SetLabelSize(0.03); histos2[mapname]->GetYaxis()->SetLabelSize(0.03);
		mapname = "trueTop_configuration"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Top matching via charge, genflavour";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 3, 0, 3, 3, 0, 3);
		histos2[mapname]->GetXaxis()->SetBinLabel(1, "l1b1"); histos2[mapname]->GetXaxis()->SetBinLabel(2, "l1b but not b1"); histos2[mapname]->GetXaxis()->SetBinLabel(3, "l1 no b match");
		histos2[mapname]->GetYaxis()->SetBinLabel(1, "l2b2"); histos2[mapname]->GetYaxis()->SetBinLabel(2, "l2b but not b2"); histos2[mapname]->GetYaxis()->SetBinLabel(3, "l2 no b match"); 
		histos2[mapname]->GetXaxis()->SetLabelSize(0.03); histos2[mapname]->GetYaxis()->SetLabelSize(0.03);
		//new in week 18/06/2012
		mapname = "Matchingconfiguration_MT2lb_lt_200"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "jet matching to genID,MID";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 3, 0, 3, 3, 0, 3);
		histos2[mapname]->GetXaxis()->SetBinLabel(1, "trueb1 matched"); histos2[mapname]->GetXaxis()->SetBinLabel(2, "trueb1 wrong matched"); histos2[mapname]->GetXaxis()->SetBinLabel(3, "trueb1 not matched");
		histos2[mapname]->GetYaxis()->SetBinLabel(1, "trueb2 matched"); histos2[mapname]->GetYaxis()->SetBinLabel(2, "trueb2 wrong matched"); histos2[mapname]->GetYaxis()->SetBinLabel(3, "trueb2 not matched"); 
		histos2[mapname]->GetXaxis()->SetLabelSize(0.03); histos2[mapname]->GetYaxis()->SetLabelSize(0.03);
		mapname = "Matchingconfiguration_MT2lb_gt_200"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "jet matching to genID,MID";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 3, 0, 3, 3, 0, 3);
		histos2[mapname]->GetXaxis()->SetBinLabel(1, "trueb1 matched"); histos2[mapname]->GetXaxis()->SetBinLabel(2, "trueb1 wrong matched"); histos2[mapname]->GetXaxis()->SetBinLabel(3, "trueb1 not matched");
		histos2[mapname]->GetYaxis()->SetBinLabel(1, "trueb2 matched"); histos2[mapname]->GetYaxis()->SetBinLabel(2, "trueb2 wrong matched"); histos2[mapname]->GetYaxis()->SetBinLabel(3, "trueb2 not matched"); 
		histos2[mapname]->GetXaxis()->SetLabelSize(0.03); histos2[mapname]->GetYaxis()->SetLabelSize(0.03);
		mapname = "MT2min_massless_vs_btaggingdiscriminant"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2lb vs btagging discriminant"; xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "b-tag discr";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 7);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_massless_vs_MT2ll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2lb vs MT2l"; xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(l,l) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_massless_vs_MT2bb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2lb vs MT2b"; xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(b,b) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_massless_vs_MT2bb_lInMET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2lb vs MT2b"; xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(b,b) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_massless_vs_MT2bb_lInMET_testmass80"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2lb vs MT2b(l in MET)"; xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(b,b) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_massless_vs_MET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2lb vs MET"; xtitle = "M_{T2} (GeV)"; ytitle = "E_{T}^{miss} (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 50, 0, 1000);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_massless_vs_HT"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2lb vs HT"; xtitle = "M_{T2} (GeV)"; ytitle = "H_{T} (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 85, 0, 1275);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_massless_vs_PV"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2lb vs #PV"; xtitle = "M_{T2} (GeV)"; ytitle = "#vertices";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 25, 0, 25);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_massless_vs_D2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2lb vs new Massfit discriminant"; xtitle = "M_{T2} (GeV)"; ytitle = "D_{2}";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 100, 0, 300);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_massless_vs_Mlbsamecombi"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2lb vs Mlb same combi"; xtitle = "M_{T2} (GeV)"; ytitle = "M_{lb} (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 100, 0, 300);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_vs_btaggingdiscriminant"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2lbmin vs btagging discriminant"; xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "b-tag discr";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 7);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_vs_MT2ll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2lb vs MT2l"; xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(l,l) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_vs_MT2bb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2lb vs MT2b"; xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(b,b) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_vs_MT2bb_lInMET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2lb vs MT2b"; xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(b,b) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_vs_MT2bb_lInMET_testmass80"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2lb vs MT2b(l in MET)"; xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(b,b) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_vs_MET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2lb vs MET"; xtitle = "M_{T2} (GeV)"; ytitle = "E_{T}^{miss} (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 50, 0, 1000);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_vs_HT"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2lb vs HT"; xtitle = "M_{T2} (GeV)"; ytitle = "H_{T} (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 85, 0, 1275);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_vs_PV"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2lb vs #PV"; xtitle = "M_{T2} (GeV)"; ytitle = "#vertices";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 25, 0, 25);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_vs_D2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2lb vs new Massfit discriminant"; xtitle = "M_{T2} (GeV)"; ytitle = "D_{2}";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 100, 0, 300);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_vs_Mlbsamecombi"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2lb vs Mlb same combi"; xtitle = "M_{T2} (GeV)"; ytitle = "M_{lb} (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 100, 0, 300);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_withMlbcut_vs_MT2ll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2min after Mlb cut vs MT2(ll)"; xtitle = "M_{T2} (GeV)"; ytitle = "M_{T2}(l,l) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_withMlbcut_vs_MT2bb_lInMET_testmass80"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2min after Mlb cut vs MT2(bb)"; xtitle = "M_{T2} (GeV)"; ytitle = "M_{T2}(b,b) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_withMlbcut_vs_MET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2min after Mlb cut vs MET"; xtitle = "M_{T2} (GeV)"; ytitle = "E_{T}^{miss} (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 50, 0, 1000);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_withMlbcut_vs_HT"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2min after Mlb cut vs HT"; xtitle = "M_{T2} (GeV)"; ytitle = "H_{T} (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 85, 0, 1275);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_withMlbcut_vs_PV"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2min after Mlb cut vs #PV"; xtitle = "M_{T2} (GeV)"; ytitle = "#vertices";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 25, 0, 25);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2min_withMlbcut_vs_D2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2min after Mlb cut vs new Massfit discriminant"; xtitle = "M_{T2} (GeV)"; ytitle = "D_{2}";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 100, 0, 300);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2ll_vs_MT2bb_lInMET_testmass80"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MTll vs MT2bb"; xtitle = "M_{T2}(l,l) (GeV)"; ytitle = "M_{T2}(b,b) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2ll_vs_MT2bb_lInMET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2ll vs MT2bb"; xtitle = "M_{T2}(l,l) (GeV)"; ytitle = "M_{T2}(b,b) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2ll_vs_MT2bb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2ll vs MT2bb"; xtitle = "M_{T2}(l,l) (GeV)"; ytitle = "M_{T2}(b,b) (GeV)";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2ll_vs_D2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2ll vs D2"; xtitle = "M_{T2}(l,l) (GeV)"; ytitle = "D_{2}";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 100, 0, 300);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2bb_vs_D2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2bb vs D2"; xtitle = "M_{T2}(b,b) (GeV)"; ytitle = "D_{2}";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 100, 0, 300);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2bb_vs_btaggingdiscriminant"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "MT2bb vs btagging discriminant"; xtitle = "M_{T2}(b,b) (GeV)"; ytitle = "b-tag discr";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 7);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "Mlb_vs_btaggingdiscriminant"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Mlb vs btagging discriminant"; xtitle = "M_{lb} (GeV)"; ytitle = "b-tag discr";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 7);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "D2_vs_btaggingdiscriminant"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "D2 vs btagging discriminant"; xtitle = "D_{2}"; ytitle = "b-tag discr";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 100, 0, 300, 70, 0, 7);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "PUWeight_vs_PV"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "PU Weight vs #PV"; xtitle = "PU Weight"; ytitle = "#vertices";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 40, 0, 2, 25, 0, 25);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		//tail plots
		if(tailMT2bb){
		mapname = "Tail_MT2bbgt125_B1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_B2Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_J1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_JPt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_L1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_L2Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_B1Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_B2Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_L1Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_L2Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_B1Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_B2Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_L1Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_L2Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_JbDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_B1bDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_B2bDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_MuoPtErrOverPt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 20, 0, 0.4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_DPhiBB"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_DPhiLL"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_DPhiLB_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_DPhiLB_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_DRBB"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_DRLL"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_DRLB_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_DRLB_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_HT"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 1250);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_MET"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 500);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_PV"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 25);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_Mll"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_Mbb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_Mlb_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_Mlb_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_NJets20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_NJets40"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_masslessMT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_massiveMT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_MT2l"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbgt125_D2"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 80, 0, 400);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		} if(tailMT2bb_linMET){
		mapname = "Tail_MT2bbLinMETgt150_B1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_B2Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_J1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_JPt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_L1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_L2Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_B1Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_B2Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_L1Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_L2Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_B1Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_B2Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_L1Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_L2Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_JbDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_B1bDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_B2bDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_MuoPtErrOverPt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 20, 0, 0.4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_DPhiBB"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_DPhiLL"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_DPhiLB_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_DPhiLB_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_DRBB"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_DRLL"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_DRLB_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_DRLB_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_HT"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 1250);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_MET"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 500);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_PV"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 25);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_Mll"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_Mbb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_Mlb_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_Mlb_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_NJets20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_NJets40"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_masslessMT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_massiveMT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_MT2l"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2bbLinMETgt150_D2"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 80, 0, 400);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		} if(tailMT2ll){
		mapname = "Tail_MT2llgt85_B1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_B2Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_J1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_JPt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_L1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_L2Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_B1Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_B2Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_L1Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_L2Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_B1Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_B2Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_L1Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_L2Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_JbDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_B1bDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_B2bDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_MuoPtErrOverPt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 20, 0, 0.4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_DPhiBB"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_DPhiLL"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_DPhiLB_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_DPhiLB_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_DRBB"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_DRLL"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_DRLB_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_DRLB_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_HT"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 1250);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_MET"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 500);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_PV"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 25);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_Mll"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_Mbb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_Mlb_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_Mlb_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_NJets20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_NJets40"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_masslessMT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_massiveMT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_MT2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_MT2b_linMET"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_MT2llgt85_D2"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 80, 0, 400);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		} if(tailMT2lb_massive){
		mapname = "Tail_massiveMT2lbgt210_B1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_B2Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_J1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_JPt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_L1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_L2Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_B1Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_B2Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_L1Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_L2Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_B1Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_B2Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_L1Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_L2Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_JbDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_B1bDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_B2bDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_MuoPtErrOverPt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 20, 0, 0.4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_DPhiBB"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_DPhiLL"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_DPhiLB_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_DPhiLB_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_DRBB"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_DRLL"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_DRLB_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_DRLB_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_HT"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 1250);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_MET"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 500);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_PV"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 25);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_Mll"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_Mbb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_Mlb_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_Mlb_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_NJets20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_NJets40"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_masslessMT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_MT2l"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_MT2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_MT2b_linMET"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_massiveMT2lbgt210_D2"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 80, 0, 400);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		} if(tailMT2bb_massless){
		mapname = "Tail_masslessMT2lbgt125_B1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_B2Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_J1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_JPt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_L1Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_L2Pt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_B1Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_B2Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_L1Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_L2Phi"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 41, -3.2, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_B1Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_B2Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_L1Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_L2Eta"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, -2.5, 4.1);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_JbDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_B1bDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_B2bDiscr"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_MuoPtErrOverPt"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 20, 0, 0.4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_DPhiBB"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_DPhiLL"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_DPhiLB_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_DPhiLB_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_DRBB"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_DRLL"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_DRLB_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_DRLB_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 33, 0, 6.6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_HT"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 1250);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_MET"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 500);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_PV"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 25);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_Mll"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_Mbb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_Mlb_all"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_Mlb_MT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_NJets20"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_NJets40"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_massiveMT2lbmin"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_MT2l"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_MT2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_MT2b_linMET"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Tail_masslessMT2lbgt125_D2"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = " ";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 80, 0, 400);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		}
		mapname = "Tail_DR_lb_configuration"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Closest DeltaR configuration";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 2, 0, 2, 2, 0, 2);
		histos2[mapname]->GetXaxis()->SetBinLabel(1, "l1 closer to b1"); histos2[mapname]->GetXaxis()->SetBinLabel(2, "l1 closer to b2"); 
		histos2[mapname]->GetYaxis()->SetBinLabel(1, "l2 closer to b1"); histos2[mapname]->GetYaxis()->SetBinLabel(2, "l2 closer to b2"); 
		histos2[mapname]->GetXaxis()->SetLabelSize(0.03); histos2[mapname]->GetYaxis()->SetLabelSize(0.03);
		mapname = "Tail_trueb1b2_configuration"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "match b1 b2 to true b";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 2, 0, 2, 2, 0, 2);
		histos2[mapname]->GetXaxis()->SetBinLabel(1, "b1"); histos2[mapname]->GetXaxis()->SetBinLabel(2, "b2"); 
		histos2[mapname]->GetYaxis()->SetBinLabel(1, "true"); histos2[mapname]->GetYaxis()->SetBinLabel(2, "fake"); 
		histos2[mapname]->GetXaxis()->SetLabelSize(0.03); histos2[mapname]->GetYaxis()->SetLabelSize(0.03);
		mapname = "Tail_trueTop_configuration"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Top matching via charge, genflavour";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 3, 0, 3, 3, 0, 3);
		histos2[mapname]->GetXaxis()->SetBinLabel(1, "l1b1"); histos2[mapname]->GetXaxis()->SetBinLabel(2, "l1b but not b1"); histos2[mapname]->GetXaxis()->SetBinLabel(3, "l1 no b match");
		histos2[mapname]->GetYaxis()->SetBinLabel(1, "l2b2"); histos2[mapname]->GetYaxis()->SetBinLabel(2, "l2b but not b2"); histos2[mapname]->GetYaxis()->SetBinLabel(3, "l2 no b match"); 
		histos2[mapname]->GetXaxis()->SetLabelSize(0.03); histos2[mapname]->GetYaxis()->SetLabelSize(0.03);
		mapname = "MT2lb_MT2ll_MT2bb"; if(av) vs.push_back(mapname); mapname = mapname+hs; title = "MT2lb vs. MT2ll vs. MT2bb";
		xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(l,l) (GeV)"; ztitle = "M_{T2}(b,b) (GeV)";
        	if(histos3.count(mapname) == 0 ) histos3[mapname] = new TH3D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700,70,0,700);
		histos3[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos3[mapname]->GetYaxis()->SetTitle(ytitle.c_str()); histos3[mapname]->GetZaxis()->SetTitle(ztitle.c_str());
		mapname = "MT2lb_MT2ll_MT2bb_lInMET"; if(av) vs.push_back(mapname); mapname = mapname+hs; title = "MT2lb vs. MT2ll vs. MT2bb with leps in MET";
		xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(l,l) (GeV)"; ztitle = "M_{T2}(b,b) (GeV)";
        	if(histos3.count(mapname) == 0 ) histos3[mapname] = new TH3D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700,70,0,700);
		histos3[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos3[mapname]->GetYaxis()->SetTitle(ytitle.c_str()); histos3[mapname]->GetZaxis()->SetTitle(ztitle.c_str());
		mapname = "MT2lb_MT2ll_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname); mapname = mapname+hs; title = "MT2lb vs. MT2ll vs. MT2bb with leps in MET, testmass mW";
		xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(l,l) (GeV)"; ztitle = "M_{T2}(b,b) (GeV)";
        	if(histos3.count(mapname) == 0 ) histos3[mapname] = new TH3D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700,70,0,700);
		histos3[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos3[mapname]->GetYaxis()->SetTitle(ytitle.c_str()); histos3[mapname]->GetZaxis()->SetTitle(ztitle.c_str());
		mapname = "MT2lbMassless_MT2ll_MT2bb"; if(av) vs.push_back(mapname); mapname = mapname+hs; title = "massless MT2lb vs. MT2ll vs. MT2bb";
		xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(l,l) (GeV)"; ztitle = "M_{T2}(b,b) (GeV)";
        	if(histos3.count(mapname) == 0 ) histos3[mapname] = new TH3D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700,70,0,700);
		histos3[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos3[mapname]->GetYaxis()->SetTitle(ytitle.c_str()); histos3[mapname]->GetZaxis()->SetTitle(ztitle.c_str());
		mapname = "MT2lbMassless_MT2ll_MT2bb_lInMET"; if(av) vs.push_back(mapname); mapname = mapname+hs; title = "massless MT2lb vs. MT2ll vs. MT2bb with leps in MET";
		xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(l,l) (GeV)"; ztitle = "M_{T2}(b,b) (GeV)";
        	if(histos3.count(mapname) == 0 ) histos3[mapname] = new TH3D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700,70,0,700);
		histos3[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos3[mapname]->GetYaxis()->SetTitle(ytitle.c_str()); histos3[mapname]->GetZaxis()->SetTitle(ztitle.c_str());
		mapname = "MT2lbMassless_MT2ll_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname); mapname = mapname+hs; title = "massless MT2lb vs. MT2ll vs. MT2bb with leps in MET, testmass mW";
		xtitle = "M_{T2}(lb,lb) (GeV)"; ytitle = "M_{T2}(l,l) (GeV)"; ztitle = "M_{T2}(b,b) (GeV)";
        	if(histos3.count(mapname) == 0 ) histos3[mapname] = new TH3D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700,70,0,700);
		histos3[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos3[mapname]->GetYaxis()->SetTitle(ytitle.c_str()); histos3[mapname]->GetZaxis()->SetTitle(ztitle.c_str());


		av = false;

	}
	}

	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Sumw2();
		h->second->GetYaxis()->SetTitle("Events"); 
		h->second-> SetStats(false);
	}
	for(map<string,TH2D*>::iterator h=histos2.begin(); h!=histos2.end();++h){
		h->second->Sumw2();
		h->second-> SetStats(false);
	}
	for(map<string,TH3D*>::iterator h=histos3.begin(); h!=histos3.end();++h){
		h->second->Sumw2();
		h->second-> SetStats(false);
	}

	cout << " done." << endl;

	histonames = vs;

	cout << "have " << histonames.size() << " different histgram types. " << endl;
	cout << "check: total numbers of 1d histograms " << (int)histos.size() << endl;
	cout << "check: total numbers of 2d histograms " << (int)histos2.size() << endl;
	cout << "check total " << (int)histos.size()+(int)histos2.size() << " = " << histonames.size() * sampletypesize * leptontypesize << endl;

	load(samples.Data());

    for(size_t i = 0; i < fSamples.size(); ++i){
        
	    int ssb = 0; int osb = 0;
\	    int vetoeddiscrcount = 0; int discrcount = 0; int globsol = 0;
	    int statvec[3]; statvec[0]=0; statvec[1]=0; statvec[2]=0;
	    int nstptot = 0; int ndertot = 0; int nrottot = 0; int nscntot = 0; int globalsolution = 0;
   	    if(runData==false && fSamples[i].type=="data") continue;
	    if(calcsusy==false && fSamples[i].type=="susy") continue;
	    if(debug && fSamples[i].type=="mc") continue;
        
	    string sampletype = (string)fSamples[i].type;
	    string leptype = "LL";
	    if(sampletype==(string)"mc"){
		if(fSamples[i].sname=="Wtolnu") sampletype = (string)"WJets";
		else if(fSamples[i].sname=="DY") sampletype = (string)"ZJets";
		else if(fSamples[i].name=="TTbar") sampletype = (string)"TTbar";
		else if(fSamples[i].sname=="TTbarV") sampletype = (string)"TTbarV";
		else if(fSamples[i].sname=="Top") sampletype = (string)"SingleTop";//no ttbar
		else if(fSamples[i].sname=="VV" || fSamples[i].sname=="VVV") sampletype = (string)"VV/VVV";
		else sampletype = (string)"Other";
	    }
	    if(sampletype==(string)"susy") sampletype=(string)"Stop";
        
	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "StopSearch: looping over " << fSamples[i].name << " added in " << sampletype << endl;
	    if(fVerbose>2) cout << "            sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

	    double stopreweightfactor = 1.;//need this since simplified model is not taken as whole but only one point
	    if(sampletype=="Stop")                                                                          stopreweightfactor = fSamples[i].nevents/50000.;
	    sample_weight = sample_weight * stopreweightfactor;

	    SolveTTbarNew *sttb = new SolveTTbarNew();//discriminant solver
	    SolveTTbar *sttbLuc = new SolveTTbar();//discriminant solver

	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts;
        
	    if( fSamples[i].type=="data" && fSamples[i].sname=="EE-Data") { myCuts += " && " + triggerEE; }//cuts to be aplied only on data
	    else if( fSamples[i].type=="data" && fSamples[i].sname=="EMu-Data") { myCuts += " && " + triggerEMu; }//cuts to be aplied only on data
	    else if( fSamples[i].type=="data" && fSamples[i].sname=="MuMu-Data") { myCuts += " && " + triggerMuMu; } //cuts to be aplied only on data
	    else if(fSamples[i].type=="data") {cout << "data not usuable" << " type " << fSamples[i].type << " Sname " << fSamples[i].sname << endl; continue; }//not usuable
        
            if(fSamples[i].type=="susy") cout << "no doing " << sampletype << endl;

   	    cout << "Cuts for Flow: " << myCuts << endl;
   	    fSamples[i].tree->Draw(">>selList", myCuts);

	    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

	    fSamples[i].tree->SetEventList(myEvtList);

	    int counter=0;
	    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
	    if(myEvtList->GetSize()==0) continue;
        while(myEvtList->GetEntry(counter++) !=-1){	
      		int jentry = myEvtList->GetEntry(counter-1);
            
            nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
            fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
            
            if ( fVerbose>2 && counter % 50000 == 0 )  cout << "+++ Proccessing event " << counter << endl;

	    if((fSamples[i].sname!="VVV"&&fSamples[i].sname!="TTbarV" && fSamples[i].type!="susy"&&sampletype!="Stop")  && fMT2tree->misc.HBHENoiseFlagIso!=0) continue;
            Double_t weight = sample_weight;
            if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 
            Bool_t recoedee   = false;// exact 2 ele, 0 muo
            Bool_t recoedemu  = false;// exact 1 ele, 1 muo
            Bool_t recoedmumu = false;// exact 0 ele, 2 muo
            if(fMT2tree->NEles>=2 &&fMT2tree->ele[0].lv.Pt()>20&&fMT2tree->ele[1].lv.Pt()>10&&fMT2tree->ele[2].lv.Pt()<10&&(fMT2tree->NMuons==0||fMT2tree->muo[0].lv.Pt()<10) && fMT2tree->ele[0].ID90==7 && fMT2tree->ele[1].ID90==7) recoedee   = true;
            if(fMT2tree->NMuons>=2&&fMT2tree->muo[0].lv.Pt()>20&&fMT2tree->muo[1].lv.Pt()>10&&fMT2tree->muo[2].lv.Pt()<10&&(fMT2tree->NEles ==0||fMT2tree->ele[0].lv.Pt()<10)) recoedmumu = true;
            if( (fMT2tree->NMuons>=1&&fMT2tree->muo[0].lv.Pt()>20&&fMT2tree->muo[1].lv.Pt()<10&&fMT2tree->NEles>=1 &&fMT2tree->ele[0].lv.Pt()>10&&fMT2tree->ele[1].lv.Pt()<10 && fMT2tree->ele[0].ID90==7) || (fMT2tree->NMuons>=1&&fMT2tree->muo[0].lv.Pt()>10&&fMT2tree->muo[1].lv.Pt()<10&&fMT2tree->NEles>=1 &&fMT2tree->ele[0].lv.Pt()>20&&fMT2tree->ele[1].lv.Pt()<10 && fMT2tree->ele[0].ID90==7) ) recoedemu  = true;
	    if(stupidchecks){
		if(fMT2tree->NEles>=2 &&fMT2tree->ele[0].lv.Pt()>17&&fMT2tree->ele[1].lv.Pt()>8&&fMT2tree->ele[2].lv.Pt()<8&&(fMT2tree->NMuons==0||fMT2tree->muo[0].lv.Pt()<8) && fMT2tree->ele[0].ID90==7 && fMT2tree->ele[1].ID90==7) recoedee   = true;
		if(fMT2tree->NMuons>=2&&fMT2tree->muo[0].lv.Pt()>17&&fMT2tree->muo[1].lv.Pt()>8&&fMT2tree->muo[2].lv.Pt()<8&&(fMT2tree->NEles ==0||fMT2tree->ele[0].lv.Pt()<8)) recoedmumu = true;
		if( (fMT2tree->NMuons>=1&&fMT2tree->muo[0].lv.Pt()>17&&fMT2tree->muo[1].lv.Pt()<8&&fMT2tree->NEles>=1 &&fMT2tree->ele[0].lv.Pt()>8&&fMT2tree->ele[1].lv.Pt()<8 && fMT2tree->ele[0].ID90==7) || (fMT2tree->NMuons>=1&&fMT2tree->muo[0].lv.Pt()>8&&fMT2tree->muo[1].lv.Pt()<8&&fMT2tree->NEles>=1 &&fMT2tree->ele[0].lv.Pt()>17&&fMT2tree->ele[1].lv.Pt()<8 && fMT2tree->ele[0].ID90==7) ) recoedemu  = true;
	    }

            Bool_t recoedosee   = false;// opposite sign
            Bool_t recoedosemu  = false;// opposite sign
            Bool_t recoedosmumu = false;// opposite sign
	    if(recoedee   && (fMT2tree->ele[0].Charge)*(fMT2tree->ele[1].Charge)==(-1)) recoedosee   = true;
	    if(recoedemu  && (fMT2tree->ele[0].Charge)*(fMT2tree->muo[0].Charge)==(-1)) recoedosemu  = true;
	    if(recoedmumu && (fMT2tree->muo[0].Charge)*(fMT2tree->muo[1].Charge)==(-1)) recoedosmumu = true;
            Bool_t recoedeenZ   = false;// off Z
            Bool_t recoedemunZ  = false;// off Z -> is not requirement, so this variable should be useless
            Bool_t recoedmumunZ = false;// off Z
	    if(recoedee   && ((fMT2tree->ele[0].lv + fMT2tree->ele[1].lv).M()<81 || (fMT2tree->ele[0].lv + fMT2tree->ele[1].lv).M()>101) ) recoedeenZ   = true;
	    if(recoedemu  && ((fMT2tree->ele[0].lv + fMT2tree->muo[0].lv).M()<81 || (fMT2tree->ele[0].lv + fMT2tree->muo[0].lv).M()>101) ) recoedemunZ  = true;
	    if(recoedmumu && ((fMT2tree->muo[0].lv + fMT2tree->muo[1].lv).M()<81 || (fMT2tree->muo[0].lv + fMT2tree->muo[1].lv).M()>101) ) recoedmumunZ = true;

	    if(!(recoedee)   && !(recoedemu)   && !(recoedmumu)  ) continue;//require dilepton
	    if(!(recoedosee) && !(recoedosemu) && !(recoedosmumu) ) continue;//require os-dilepton, maybe comment this for background est.
	    //if(!(recoedeenZ) ||                   !(recoedmumunZ) ) continue;//require of Zpeak

	    if(recoedee)   leptype = "EE";
	    if(recoedemu)  leptype = "EMu";
	    if(recoedmumu) leptype = "MuMu";

		//trigger weight; EMuweight; MUMUweight, if different lumis;
		//taken from AN2011_466_v4
	    if(recoedee   && !(fMT2tree->misc.isData) ) weight = 1.00 * weight;
	    if(recoedemu  && !(fMT2tree->misc.isData) ) weight = 0.95 * weight;
	    if(recoedmumu && !(fMT2tree->misc.isData) ) weight = 0.92 * weight;
	    if(fMT2tree->misc.isData && recoedee   && fSamples[i].sname!="EE-Data"  ) continue;
	    if(fMT2tree->misc.isData && recoedemu  && fSamples[i].sname!="EMu-Data" ) continue;
	    if(fMT2tree->misc.isData && recoedmumu && fSamples[i].sname!="MuMu-Data") continue;

	    unsigned int NumBJets = fMT2tree->GetNBtags(Tagger,discr,bpt,2.4,1);

	    float SFweightErr = 0;
	    float SFweight = 1;//outside, since need it there later
	    if(!(dofastbtagSFreweighting)){
	    if((!fMT2tree->misc.isData)){
		vector<float> jetEff;
		vector<float> jetEffErr;
		vector<float> jetSF;
		vector<float> jetSFErr;
		int njetsusuable = 0;
		for(int n = 0; n<fMT2tree->NJets; ++n){
			if(fMT2tree->jet[n].isPFIDLoose==false) continue;
			if(fMT2tree->jet[n].Flavour<=-7777) continue;///XXXX: might need change
			float effPt  = fMT2tree->jet[n].lv.Pt();
			float effEta = fabs(fMT2tree->jet[n].lv.Eta());
			++njetsusuable;
			if(abs(fMT2tree->jet[n].Flavour)==5){///XXXX: might need change
					jetEff.push_back( float(hbeff->GetBinContent(hbeff->FindBin(effPt))) );
					jetEffErr.push_back( float(hbeff->GetBinError(hbeff->FindBin(effPt))) );
				float SFErr;
				float SF = getBTagSF(SFErr, tagger, effPt, effEta);
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr );//first implementation
			}
			else if(abs(fMT2tree->jet[n].Flavour)==4){///XXXX: might need change
					jetEff.push_back( float(hceff->GetBinContent(hceff->FindBin(effPt))) );
					jetEffErr.push_back( float(hceff->GetBinError(hceff->FindBin(effPt))) );
				float SFErr;
				float SF = getBTagSF(SFErr, tagger, effPt, effEta );
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr*2.);
			}
			else {
					jetEff.push_back( float(hleff->GetBinContent(hleff->FindBin(effPt))) );
					jetEffErr.push_back( float(hleff->GetBinError(hleff->FindBin(effPt))) );
				float SFErr;
				float SF = getMistagSF(SFErr, tagger, effPt, effEta, 0 );
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr );
			}
		}
		SFweight = getBTagEventWeightError(SFweightErr, jetEff, jetEffErr, jetSF, jetSFErr, ntags);
		if(SFweight==0){
			if(njetsusuable!=0){
				cout << "Event has zero weight, do not use it" << endl;
				continue;
			}
			else { //event has no flavour information, use average event weight
				SFweight = 0.945572*abs(ntags);
				if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
			}
		}
		weight  = weight * SFweight;
	    }
            }
            else{
	     //from ttbar payload
             if(tagger=="SSVHPT"){
		SFweight = pow(0.95*abs(ntags),2);
		if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
		weight = weight * SFweight;
	     }
             if(tagger=="SSVHEM"){
		SFweight = pow(0.96,abs(ntags));
		if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
		weight = weight * SFweight;
	     }
	    }

	   TLorentzVector l1, l2, met;
	   vector<TLorentzVector> bjets; bjets.clear();
	   vector<int> bjetsflavour; bjetsflavour.clear();
	   vector<float> bjetsdiscriminant; bjetsdiscriminant.clear();
	   int l1c(-999), l2c(-999)
	   Bool_t l1e(true), l2e(true);
	   vector<int> jind; jind.clear();
	   vector<int> bind; bind.clear();
	   if(recoedee){ 
		l1 = (fMT2tree->ele[0]).lv; l2 = (fMT2tree->ele[1]).lv;
		l1c = (fMT2tree->ele[0]).Charge; l2c = fMT2tree->ele[1].Charge;
		l1e=true; l2e=true;}
	   else if(recoedemu){
		if(fMT2tree->ele[0].lv.Pt()>fMT2tree->muo[0].lv.Pt() ) { 
			l1 = (fMT2tree->ele[0]).lv; l2 = (fMT2tree->muo[0]).lv; 
			l1c = (fMT2tree->ele[0]).Charge; l2c = fMT2tree->muo[0].Charge;
			l1e=true; l2e=false;}
		else  { 
			l1 = (fMT2tree->muo[0]).lv; l2 = (fMT2tree->ele[0]).lv; 
			l1c = fMT2tree->muo[0].Charge; l2c = fMT2tree->ele[0].Charge;
			l1e=false; l2e=true;}
	   }
	   else if(recoedmumu) { 
		l1 = (fMT2tree->muo[0]).lv; l2 = (fMT2tree->muo[1]).lv; 
		l1c = fMT2tree->muo[0].Charge; l2c = fMT2tree->muo[1].Charge;
		l1e=false; l2e=false;}
	   for(int n = 0; n<fMT2tree->NJets; ++n){//loop only over bjets
		if(fMT2tree->jet[n].isPFIDLoose==false) continue;
		if(fMT2tree->jet[n].lv.Pt()<bpt) continue;
		if(fabs(fMT2tree->jet[n].lv.Eta())>2.4) continue;
		float btempdiscr = (tagger=="SSVHPT" ? fMT2tree->jet[n].bTagProbSSVHP : tagger=="SSVHEM" ? fMT2tree->jet[n].bTagProbSSVHE :fMT2tree->jet[n].bTagProbSSVHE);//default SSVHE
		jind.push_back(n);//NOTE: add here an additional pt constraint?????
		if(btempdiscr<discr) continue;
		bjets.push_back(fMT2tree->jet[n].lv); bjetsflavour.push_back(fMT2tree->jet[n].Flavour); bjetsdiscriminant.push_back(btempdiscr); bind.push_back(n);
	   }
	   if(NumBJets!=bjets.size()) cout << "NumBJets " << NumBJets << " bjets.size() " << bjets.size() << endl;
	   if(NumBJets<2) cout << "NumBJets " << NumBJets << endl;
	   //sort 
	   for(size_t n =0; n<bjets.size();++n){
		for(size_t m =n+1; m<bjets.size();++m){
			if(bjetsdiscriminant[m]<=bjetsdiscriminant[n]) continue;
			swap(bjetsdiscriminant[n],bjetsdiscriminant[m]);
			swap(bjets[n],bjets[m]);
			swap(bjetsflavour[n],bjetsflavour[m]);
			swap(bind[n],bind[m]);
		}
	  }

	   met = fMT2tree->pfmet[0];
	   //genlevel plots only
	   TLorentzVector genl1, genl2, genb1, genb2, genmet;
	   int gfl1(-99), gfl2(-99), gfb1(-99), gfb2(-99);
	   int gil1(-99), gil2(-99), gib1(-99), gib2(-99);
	   int gin1(-99), gin2(-99);
	    vector<TLorentzVector> genlv; genlv.clear();
	    vector<int> genflav; genflav.clear();
	    vector<int> genindex; genindex.clear();
	    vector<int> genneutrinoindex; genneutrinoindex.clear();
	    vector<int> genneutrinoflavour; genneutrinoflavour.clear();
	   for(int n = 0; n<25; ++n){
		bool add = true;
		int glid = fMT2tree->genlept[n].ID; int glmid = fMT2tree->genlept[n].MID; int glgmid = fMT2tree->genlept[n].GMID;
		//bool addneutrino = true;
		     if( (fabs(glid==11||glid==13) && ((glid>0&&glmid==-24&&glgmid==-6)||(glid<0&&glmid==24&&glmid==6)) && fMT2tree->genlept[n].lv.Pt()>1) ||
			 (((glid==5&&glmid==6)||(glid==-5&&glmid==-6)) && fMT2tree->genlept[n].lv.Pt()>1) ){
				for(unsigned int m = 0; m<genindex.size(); ++m){
					if(fMT2tree->genlept[n].ID  != fMT2tree->genlept[genindex[m] ].ID  ) continue;
					if(fMT2tree->genlept[n].MID != fMT2tree->genlept[genindex[m] ].MID ) continue;
					if(fMT2tree->genlept[n].GMID!= fMT2tree->genlept[genindex[m] ].GMID) continue;
					if(fMT2tree->genlept[n].lv.DeltaR(  fMT2tree->genlept[genindex[m] ].lv)<0.1 ||  
					   fMT2tree->genlept[n].lv.DeltaPhi(fMT2tree->genlept[genindex[m] ].lv)<0.1 ||
					   fabs(fMT2tree->genlept[n].lv.Eta()-fMT2tree->genlept[genindex[m] ].lv.Eta() )<0.07 ) {
						if(fabs(fMT2tree->genlept[n].lv.Pt()-fMT2tree->genlept[genindex[m] ].lv.Pt() )<0.1*fMT2tree->genlept[genindex[m] ].lv.Pt() ){
						//don't add duplicates
						add = false;
						break;}}
				}
				if(add){ genflav.push_back(fMT2tree->genlept[n].ID); genlv.push_back(fMT2tree->genlept[n].lv); genindex.push_back(n); }
		     }
		     if( (fabs(glid)==12||fabs(glid)==14) && ((glid>0&&glmid==24&&glmid==6) || (glid<0&&glmid==-24&&glgmid==-6)) && fMT2tree->genlept[n].lv.Pt()>1) {
				for(unsigned int m = 0; m<genindex.size(); ++m){
					if(fMT2tree->genlept[n].ID  != fMT2tree->genlept[genindex[m] ].ID  ) continue;
					if(fMT2tree->genlept[n].MID != fMT2tree->genlept[genindex[m] ].MID ) continue;
					if(fMT2tree->genlept[n].GMID!= fMT2tree->genlept[genindex[m] ].GMID) continue;
					if(fMT2tree->genlept[n].lv.DeltaR(  fMT2tree->genlept[genindex[m] ].lv)<0.1 || 
					   fMT2tree->genlept[n].lv.DeltaPhi(fMT2tree->genlept[genindex[m] ].lv)<0.1 ||
					   fabs(fMT2tree->genlept[n].lv.Eta()-fMT2tree->genlept[genindex[m] ].lv.Eta() )<0.07 ) {
						if(fabs(fMT2tree->genlept[n].lv.Pt()-fMT2tree->genlept[genindex[m] ].lv.Pt() )<0.1*fMT2tree->genlept[genindex[m] ].lv.Pt() ){
						//don't add duplicates
						add = false;
						break;}}
				}
				if(add){ genneutrinoflavour.push_back(fMT2tree->genlept[n].ID); genneutrinoindex.push_back(n); }
		     }
	   }
	if(genflav.size()==4){
		for(int n =0; n<4; ++n){
			if(abs(genflav[n])==11 || abs(genflav[n])==13){
				if(gfl1<-50){ gfl1 = genflav[n]; genl1 = genlv[n];  gil1 = genindex[n]; }
				else if(gfl2<-50){ genl2 = genlv[n]; gfl2 = genflav[n]; gil2 = genindex[n];}
			}
		}
		for(int n =0; n<4; ++n){
			if(gfl1<-50 || gfl2<-50) continue;
			if(abs(genflav[n])==5){
				if(genflav[n]<0){ if(gfl1>0){ gfb1 = genflav[n]; genb1 = genlv[n]; gib1 = genindex[n]; }
						  else if(gfl2>0){  gfb2 = genflav[n]; genb2 = genlv[n]; gib2 = genindex[n]; } 
				}
				else if(genflav[n]>0){ if(gfl1<0){  gfb1 = genflav[n]; genb1 = genlv[n]; gib1 = genindex[n]; } 
						  	else if(gfl2<0){  gfb2 = genflav[n]; genb2 = genlv[n]; gib2 = genindex[n]; } 
				}
			}
		}
		if(genneutrinoindex.size()==2){
			for(int n = 0; n<2;++n){
				if(gfl1*genneutrinoflavour[n]==-182||gfl1*genneutrinoflavour[n]==-132) {  gin1=genneutrinoindex[n]; }
				if(gfl2*genneutrinoflavour[n]==-182||gfl2*genneutrinoflavour[n]==-132) { gin2=genneutrinoindex[n]; }
			}
		}
	}

	bool truettbarEE(false), truettbarEMu(false), truettbarMuMu(false);
	if(gfl1>-50 && gfl2>-50 && gfb1 >-50 && gfb2>-50){
		if((abs(gfl1)==11 && gfl1*gfb1<0) && (gfl1*gfl2==-121 && gfl2*gfb2<0))                   truettbarEE   = true;
		if((abs(gfl1)==13 && gfl1*gfb1<0) && (gfl1*gfl2==-169 && gfl2*gfb2<0))                   truettbarMuMu = true;
		if( ((abs(gfl1)==11||abs(gfl1)==13) && gfl1*gfb1<0) && (gfl2*gfl1==-143 && gfl2*gfb2<0)) truettbarEMu  = true;
	}

	float genMl1b1(-99.), genMl2b2(-99.), genMT2(-99.), genMT2massless(-99.), genM1(-99.), genM2(-99.);
	float genMl1b1_noPtcut(-99.), genMl2b2_noPtcut(-99.), genMT2_noPtcut(-99.), genMT2massless_noPtcut(-99.), genM1_noPtcut(-99.), genM2_noPtcut(-99.);
	if(truettbarEE || truettbarEMu || truettbarMuMu){
	   genMl1b1_noPtcut      = (genl1+genb1).M();
	   genMl2b2_noPtcut      = (genl2+genb2).M();
	   genMT2_noPtcut        = fMT2tree->CalcMT2(0., true,  genl1+genb1, genl2+genb2, fMT2tree->genmet[0]);
	   genMT2massless_noPtcut= fMT2tree->CalcMT2(0., false, genl1+genb1, genl2+genb2, fMT2tree->genmet[0]);
	   if(gin2>-50 && gin1>-50){
		genM1_noPtcut = (genl1+genb1+fMT2tree->genlept[gin1].lv).M();
		genM2_noPtcut = (genl2+genb2+fMT2tree->genlept[gin2].lv).M();
	   }
	   if( ((genl1.Pt()>20&&genl2.Pt()>10)||(genl2.Pt()>20&&genl1.Pt()>10)) && genb1.Pt()>bpt && genb2.Pt()>bpt){
	   	genMl1b1      = (genl1+genb1).M();
	   	genMl2b2      = (genl2+genb2).M();
	   	genMT2        = fMT2tree->CalcMT2(0., true,  genl1+genb1, genl2+genb2, fMT2tree->genmet[0]);
	   	genMT2massless= fMT2tree->CalcMT2(0., false, genl1+genb1, genl2+genb2, fMT2tree->genmet[0]);
	   	if(gin2>-50 && gin1>-50){
			genM1 = (genl1+genb1+fMT2tree->genlept[gin1].lv).M();
			genM2 = (genl2+genb2+fMT2tree->genlept[gin2].lv).M();
	   	}
	   }
	}
	//finished genlevel plots
	   int njetsID20(0), njetsID30(0), njetsID40(0), njetsID50(0), nbjets20(0), nbjets30(0), nbjets40(0), nbjets50(0);
	   for(int n = 0; n<fMT2tree->NJets; ++n){
		if(fabs(fMT2tree->jet[n].lv.Eta())>=2.4) continue; 
		float jpt = fMT2tree->jet[n].lv.Pt();
		if(fMT2tree->jet[n].isPFIDLoose==false) continue;
		++njetsID20;
		if(fMT2tree->jet[n].bTagProbSSVHE>1.74)            ++nbjets20;
		if(fMT2tree->jet[n].bTagProbSSVHE>1.74 && jpt>30.) ++nbjets30;
		if(fMT2tree->jet[n].bTagProbSSVHE>1.74 && jpt>40.) ++nbjets40;
		if(fMT2tree->jet[n].bTagProbSSVHE>1.74 && jpt>50.) ++nbjets50;
		if(jpt>30.) ++njetsID30;
		if(jpt>40.) ++njetsID40;
		if(jpt>50.) ++njetsID50;
	   }
	   if(njetsID20 != fMT2tree->NJetsIDLoose)   cout << "njetsID20 " << njetsID20 << " fMT2tree->NJetsIDLoose "   << fMT2tree->NJetsIDLoose   << endl;
	   if(njetsID40 != fMT2tree->NJetsIDLoose40) cout << "njetsID40 " << njetsID40 << " fMT2tree->NJetsIDLoose40 " << fMT2tree->NJetsIDLoose40 << endl;
	   if(njetsID50 != fMT2tree->NJetsIDLoose50) cout << "njetsID50 " << njetsID50 << " fMT2tree->NJetsIDLoose50 " << fMT2tree->NJetsIDLoose50 << endl;
	   if(njetsID40<njetsID50)                   cout << "njetsID40 " << njetsID40 << " njetsID50 " << njetsID50 << endl;
	//lb with b vector
	   float MT2minV(-999.99), MT2minV_massless(-999.99);
	   float MT2minV_withMlb(-999.99), MT2minV_massless_withMlb(-999.99);//as this might be difficult to get
	   float MT2_l1b[NumBJets][NumBJets];//l1 with first b[], l2 with second b[];
	   float MT2_massless_l1b[NumBJets][NumBJets];//l1 with first b[], l2 with second b[];
	   float Ml1b[NumBJets];
	   float Ml2b[NumBJets];
	   float Mbbarr[NumBJets][NumBJets];
	   float MT2bbarr[NumBJets][NumBJets];
	   float MT2bb_linMET[NumBJets][NumBJets];
	   float MT2bb_mW_linMET[NumBJets][NumBJets];
	   float MT2bbmin(-999.99), MT2bbmin_mW_linMET(-999.99), MT2bbmin_linMET(-999.99);
	   float DRbbarr[NumBJets][NumBJets], DPhibbarr[NumBJets][NumBJets];
	   float DRl1b[NumBJets], DRl2b[NumBJets], DPhil1b[NumBJets], DPhil2b[NumBJets];
	   int ind_l1bxMT2min(-1), ind_l2bxMT2min(-1), ind_l1bxMT2minMlb(-1), ind_l2bxMT2minMlb(-1);
	   int ind_b1x_MT2bbmin(-1), ind_b2x_MT2bbmin(-1), ind_b1x_MT2bbmin_mW_linMET(-1), ind_b2x_MT2bbmin_mW_linMET(-1), ind_b1x_MT2bbmin_linMET(-1), ind_b2x_MT2bbmin_linMET(-1), ind_l1bxMT2minMassless(-1), ind_l2bxMT2minMassless(-1);
	   float MT2ll(-999.99), Mll(-999.99), DPhill(-999.99), DRll(-999.99);//check that this is not defined above anymore

	   float UTMvec[NumBJets][NumBJets];
	   bool  truel1bxtop[NumBJets], truel2bxtop[NumBJets], truebtop[NumBJets];//use truebtop as truebfromtop --> do genfitting
	   int  truel1flav(-99), truel2flav(-99);
	   float dR1 = 99, dR2 = 99; float dRb[NumBJets];
	   int match1(-1),match2(-1); int matchb[NumBJets];

	   //reset all arrays
	   for(unsigned int i1 = 0; i1<NumBJets; ++i1){
		DRl1b[i1] = -999.99; DRl2b[i1] = -999.99; DPhil1b[i1] = -999.99; DPhil2b[i1] = -999.99; 
		Ml1b[i1] = -999.99; Ml2b[i1] = -999.99;
		truel1bxtop[i1] = false; truel2bxtop[i1] = false; truebtop[i1] = false;
		dRb[i1] = 99.; matchb[i1] = -1;
		for(unsigned int i2 = 0; i2<NumBJets; ++i2){
			MT2_l1b[i1][i2] = -999.99; MT2_massless_l1b[i1][i2] = -999.99;
			MT2bbarr[i1][i2] = -999.99; MT2bb_mW_linMET[i1][i2] = -999.99; MT2bb_linMET[i1][i2] = -999.99; 
			DRbbarr[i1][i2] = -999.99; DPhibbarr[i1][i2] = -999.99; UTMvec[i1][i2]=-999.99; Mbbarr[i1][i2] = -999.;
		}
	   }

	   //now do genlevelmatching first - dR matching only for leptons dR<0.3, dR<0.5 for jets
	   for(int n = 0; n<25; ++n){
		//leptons
		if(fMT2tree->genlept[n].lv.Pt()<1.) continue;
		if(fMT2tree->genlept[n].lv.DeltaR(l1)<0.3 && fMT2tree->genlept[n].lv.DeltaR(l1)<dR1 &&
		   (abs(fMT2tree->genlept[n].ID)==11||abs(fMT2tree->genlept[n].ID)==13||abs(fMT2tree->genlept[n].ID)==15) ){
			match1=n; dR1 = fMT2tree->genlept[n].lv.DeltaR(l1); }
		if(fMT2tree->genlept[n].lv.DeltaR(l2)<0.3 && fMT2tree->genlept[n].lv.DeltaR(l2)<dR2 &&
		   (abs(fMT2tree->genlept[n].ID)==11||abs(fMT2tree->genlept[n].ID)==13||abs(fMT2tree->genlept[n].ID)==15) ){
			match2=n; dR2 = fMT2tree->genlept[n].lv.DeltaR(l2); }
	   	for(unsigned int i1 = 0; i1<NumBJets; ++i1){
			//bjets
			if((abs(fMT2tree->genlept[n].ID)==5) &&fMT2tree->genlept[n].lv.DeltaR(bjets[i1])<0.3 && fMT2tree->genlept[n].lv.DeltaR(bjets[i1])<dRb[i1]){
				matchb[i1]=n; dRb[i1] = fMT2tree->genlept[n].lv.DeltaR(bjets[i1]); }
		}
	   }
	   if(match1>=0) truel1flav = fMT2tree->genlept[match1].ID;
	   if(match2>=0) truel2flav = fMT2tree->genlept[match2].ID;

	   MT2ll  = fMT2tree->CalcMT2(0., true,  l1, l2, met);
	   Mll    = (l1+l2).M();
	   DPhill = fabs(l1.DeltaPhi(l2));
	   DRll   = l1.DeltaR(l2);

	   int ntruetops(0), ntruel1tops(0),ntruel2tops(0);
	   for(unsigned int i1 = 0; i1<NumBJets; ++i1){
		//first genlevelthings , 30 = 5*6 = -5*-6, -55 = 5*-11=-5*11, -65 = 5*-13 = 13*-5
		if(matchb[i1]>=0){
			int n = matchb[i1];//gen parton b is matched to this jet
			if(abs(fMT2tree->genlept[n].ID)==5 && fMT2tree->genlept[n].ID*fMT2tree->genlept[n].MID==30) { truebtop[i1] = true;
				if((fMT2tree->genlept[n].ID*truel1flav==-55 || fMT2tree->genlept[n].ID*truel1flav==-65)&& truebtop[i1]){
					//b from top, now test also lepton, if W->tau->l say also true, might fail for TTV
					if(fMT2tree->genlept[match1].GMID==fMT2tree->genlept[n].MID&&abs(fMT2tree->genlept[match1].MID)==24) truel1bxtop[i1]=true;
					else if(fabs(fMT2tree->genlept[match1].GMID)==24&&fabs(fMT2tree->genlept[match1].MID)==15) truel1bxtop[i1]=true;//believe
				}
				if((fMT2tree->genlept[n].ID*truel2flav==-55 || fMT2tree->genlept[n].ID*truel2flav==-65)&& truebtop[i1]){
					//b from top, now test also lepton, if W->tau->l say also true, might fail for TTV
					if(fMT2tree->genlept[match2].GMID==fMT2tree->genlept[n].MID&&abs(fMT2tree->genlept[match2].MID)==24) truel2bxtop[i1]=true;
					else if(fabs(fMT2tree->genlept[match2].GMID)==24&&fabs(fMT2tree->genlept[match2].MID)==15) truel2bxtop[i1]=true;//believe
				}
			}
		}
		//now real things
		DRl1b[i1]   = l1.DeltaR(bjets[i1]);
		DRl2b[i1]   = l2.DeltaR(bjets[i1]);
		DPhil1b[i1] = fabs(l1.DeltaPhi(bjets[i1]));
		DPhil2b[i1] = fabs(l2.DeltaPhi(bjets[i1]));
		Ml1b[i1]    = (l1+bjets[i1]).M();
		Ml2b[i1]    = (l2+bjets[i1]).M();
	   }

	   for(unsigned int i1 = 0; i1<NumBJets; ++i1){//need to do this like this since we need Mlb
		for(unsigned int i2 = 0; i2<NumBJets; ++i2){
			if(i2==i1) continue;
			DRbbarr[i1][i2] = (bjets[i1]).DeltaR(bjets[i2]);
			DPhibbarr[i1][i2] = fabs((bjets[i1]).DeltaPhi(bjets[i2]));
			MT2_l1b[i1][i2]              = fMT2tree->CalcMT2(0., true,  l1+bjets[i1], l2+bjets[i2], met);
	   		MT2_massless_l1b[i1][i2]     = fMT2tree->CalcMT2(0., false, l1+bjets[i1], l2+bjets[i2], met);
			MT2bbarr[i1][i2]                   = fMT2tree->CalcMT2(0.,   true,  bjets[i1], bjets[i2], met);
			MT2bb_mW_linMET[i1][i2] = fMT2tree->CalcMT2(80.4, true,  bjets[i1], bjets[i2], met+l1+l2);
			MT2bb_linMET[i1][i2] = fMT2tree->CalcMT2(0., true,  bjets[i1], bjets[i2], met+l1+l2);
			UTMvec[i1][i2] = sqrt(pow((l1+l2+bjets[i1]+bjets[i2]+met).Px(),2) + pow((l1+l2+bjets[i1]+bjets[i2]+met).Py(),2));
			Mbbarr[i1][i2] = (bjets[i1]+bjets[i2]).M();

		}
	   }
	   for(unsigned int i1 = 0; i1<NumBJets; ++i1){
		for(unsigned int i2 = 0; i2<NumBJets; ++i2){
			if(i2==i1) continue;
			if(MT2minV<0. && MT2_l1b[i1][i2]>=0.) { MT2minV = MT2_l1b[i1][i2];ind_l1bxMT2min = i1; ind_l2bxMT2min = i2; }
			else if (MT2minV>=0. && MT2_l1b[i1][i2]>=0. && MT2_l1b[i1][i2]<MT2minV)   {  
				MT2minV = MT2_l1b[i1][i2]; ind_l1bxMT2min = i1; ind_l2bxMT2min = i2; }
			if(MT2minV_massless<0. && MT2_massless_l1b[i1][i2]>=0.) { MT2minV_massless = MT2_massless_l1b[i1][i2]; ind_l1bxMT2minMassless = i1; ind_l2bxMT2minMassless = i2; }
			else if (MT2minV_massless>=0. && MT2_massless_l1b[i1][i2]>=0. && MT2_massless_l1b[i1][i2]<MT2minV) {   
				 MT2minV_massless = MT2_massless_l1b[i1][i2]; ind_l1bxMT2minMassless = i1; ind_l2bxMT2minMassless = i2; }
			if(MT2minV_withMlb<0. && MT2_l1b[i1][i2]>=0. && Ml1b[i1]<180. && Ml2b[i2]<180. ) {
				 MT2minV_withMlb = MT2_l1b[i1][i2]; ind_l1bxMT2minMlb = i1; ind_l2bxMT2minMlb = i2; }
			else if(MT2minV_withMlb>=0. && MT2_l1b[i1][i2]>=0. && Ml1b[i1]<180. && Ml2b[i2]<180. && MT2minV_withMlb>MT2_l1b[i1][i2]) {
				 MT2minV_withMlb = MT2_l1b[i1][i2]; ind_l1bxMT2minMlb = i1; ind_l2bxMT2minMlb = i2; }
			if(MT2minV_massless_withMlb<0. && MT2_massless_l1b[i1][i2]>=0. && Ml1b[i1]<180. && Ml2b[i2]<180. ) 
				MT2minV_massless_withMlb = MT2_massless_l1b[i1][i2];
			else if(MT2minV_massless_withMlb>=0. && MT2_massless_l1b[i1][i2]>=0. && Ml1b[i1]<180. && Ml2b[i2]<180. && MT2minV_massless_withMlb>MT2_massless_l1b[i1][i2]) MT2minV_massless_withMlb = MT2_massless_l1b[i1][i2];

			if(MT2bbmin<0. && MT2bbarr[i1][i2]>=0.) { MT2bbmin = MT2bbarr[i1][i2]; ind_b1x_MT2bbmin = i1; ind_b2x_MT2bbmin = i2; }
			else if(MT2bbmin>=0. && MT2bbarr[i1][i2]>=0. && MT2bbarr[i1][i2]<MT2bbmin) {  
				MT2bbmin = MT2bbarr[i1][i2]; ind_b1x_MT2bbmin = i1; ind_b2x_MT2bbmin = i2; }
			if(MT2bbmin_mW_linMET<0. && MT2bb_mW_linMET[i1][i2]>=0.){ 
				MT2bbmin_mW_linMET = MT2bb_mW_linMET[i1][i2]; ind_b1x_MT2bbmin_mW_linMET = i1; ind_b2x_MT2bbmin_mW_linMET = i2; }
			else if(MT2bbmin_mW_linMET>=0. && MT2bb_mW_linMET[i1][i2]>=0. && MT2bb_mW_linMET[i1][i2]<MT2bbmin_mW_linMET) { 
				MT2bbmin_mW_linMET = MT2bb_mW_linMET[i1][i2]; ind_b1x_MT2bbmin_mW_linMET = i1; ind_b2x_MT2bbmin_mW_linMET = i2; }
			if(MT2bbmin_linMET<0. && MT2bb_linMET[i1][i2]>=0.){ 
				MT2bbmin_linMET = MT2bb_linMET[i1][i2]; ind_b1x_MT2bbmin_linMET = i1; ind_b2x_MT2bbmin_linMET = i2; }
			else if(MT2bbmin_linMET>=0. && MT2bb_linMET[i1][i2]>=0. && MT2bb_linMET[i1][i2]<MT2bbmin_linMET) { 
				MT2bbmin_linMET = MT2bb_linMET[i1][i2]; ind_b1x_MT2bbmin_linMET = i1; ind_b2x_MT2bbmin_linMET = i2; }
		}
	   }

	//new checks 14/06/2012
	//first gen reco matching if possible -- using produced vectors genlv, genflav, genindex
	//these vectors contain all b's and leptons coming from top(!!!!), via flav can also match to same top//no tau decays in here
	//match indices via dR only, no pT
	int matchedgenl_tol1_regardlessFlavourAndCharge(-1), matchedgenl_tol2_regardlessFlavourAndCharge(-1);
	int matchedgenl_tol1_regardlessCharge(-1), matchedgenl_tol2_regardlessCharge(-1);
	int matchedgenl_tol1(-1), matchedgenl_tol2(-1);
	int matchedgenp_tob1(-1), matchedgenp_tob2(-1);//b1 = bind[ind_l1bxMT2min], b2 = bind[ind_l2bxMT2min];
	float dR_ml1fc(99.), dR_ml2fc(99.), dR_ml1c(99.), dR_ml2c(99.), dR_ml1(99.), dR_ml2(99.), dR_mb1(99.), dR_mb2(99.);
	for(int n = 0; n<25; ++n){
		if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
		int genlepID = fMT2tree->genlept[n].ID;
		int genlepMID = fMT2tree->genlept[n].MID;
		int genlepGMID = fMT2tree->genlept[n].GMID;
		if(((abs(genlepID)==11)||(abs(genlepID)==13))){
			if(abs(genlepMID)!=24 && abs(genlepMID)!=15 )     continue;
			if(abs(genlepGMID)!=6 && abs(genlepMID)==24 )     continue;
			if(abs(genlepGMID)!=24&& abs(genlepMID)==15 )     continue;
			if(abs(genlepMID)==24 && genlepMID *genlepID>=0 ) continue;
			if(abs(genlepMID)==15 && genlepMID *genlepID<=0 ) continue;
			if(genlepGMID*genlepID>=0) continue;
			float dRl1 = fMT2tree->genlept[n].lv.DeltaR(l1); 
			if(dRl1<0.3) {
				if(dRl1>=dR_ml1fc) continue;
				matchedgenl_tol1_regardlessFlavourAndCharge = n; dR_ml1fc = dRl1; }
			if(dRl1<0.3) {
				if(dRl1>=dR_ml1c) continue;
				if(l1e==true &&abs(genlepID)!=11) continue;
				if(l1e==false&&abs(genlepID)!=13) continue;
				matchedgenl_tol1_regardlessCharge = n; dR_ml1c= dRl1; }
			if(dRl1<0.3) {
				if(dRl1>=dR_ml1) continue;
				if(l1e==true &&abs(genlepID)!=11) continue;
				if(l1e==false&&abs(genlepID)!=13) continue;
				if(l1c*genlepID>0) continue;
				matchedgenl_tol1 = n; dR_ml1 = dRl1; }//11*-1 or -11*1 --> true
		}
	}
	for(int n = 0; n<25; ++n){
		if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
		int genlepID = fMT2tree->genlept[n].ID;
		int genlepMID = fMT2tree->genlept[n].MID;
		int genlepGMID = fMT2tree->genlept[n].GMID;
		if(((abs(genlepID)==11)||(abs(genlepID)==13))){
			//include tau decays with good hopes
			if(abs(genlepMID)!=24 && abs(genlepMID)!=15 )     continue;
			if(abs(genlepGMID)!=6 && abs(genlepMID)==24 )     continue;
			if(abs(genlepGMID)!=24&& abs(genlepMID)==15 )     continue;
			if(abs(genlepMID)==24 && genlepMID *genlepID>=0 ) continue;
			if(abs(genlepMID)==15 && genlepMID *genlepID<=0 ) continue;
			if(genlepGMID*genlepID>=0) continue;
			float dRl2 = fMT2tree->genlept[n].lv.DeltaR(l2);
			if(dRl2<0.3 && dRl2<dR_ml2fc) {
				if(dRl2>=dR_ml2fc) continue;
				matchedgenl_tol2_regardlessFlavourAndCharge = n; dR_ml2fc = dRl2; }
			if(dRl2<0.3 && dRl2<dR_ml2c) {
				if(dRl2>=dR_ml2c) continue;
				if(l2e==true &&abs(genlepID)!=11 && abs(genlepMID)!=15) continue;
				if(l2e==false&&abs(genlepID)!=13 && abs(genlepMID)!=15) continue;
				matchedgenl_tol2_regardlessCharge = n; dR_ml2c = dRl2; }
			if(dRl2<0.3 && dRl2<dR_ml2) {
				if(dRl2>=dR_ml2) continue;
				if(l2e==true &&abs(genlepID)!=11 && abs(genlepMID)!=15) continue;
				if(l2e==false&&abs(genlepID)!=13 && abs(genlepMID)!=15) continue;
				if(l2c*genlepID>0) continue;
				matchedgenl_tol2 = n; dR_ml2 = dRl2; }//11*-1 or -11*1 --> true
		}
	}
	for(int n = 0; n<25; ++n){
		if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
		int genlepID = fMT2tree->genlept[n].ID;
		int genlepMID = fMT2tree->genlept[n].MID;
		//int genlepGMID = fMT2tree->genlept[n].GMID;
		if(abs(genlepID)==5){
			if(abs(genlepMID)!=6)     continue;
			if(genlepMID*genlepID<=0) continue;
			float dRb1 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l1bxMT2min]);
			if(dRb1<0.5 && dRb1<dR_mb1){
				matchedgenp_tob1 = n; dR_mb1 = dRb1; }
			float dRb2 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l2bxMT2min]);
			if(dRb2<0.5 && dRb2<dR_mb2){
				matchedgenp_tob2 = n; dR_mb2 = dRb1; }
		}
	}

	int gl1(-1), gl2(-1), gb1(-1), gb2(-1);
	vector<int> gvv; gvv.clear();
	   for(int n = 0; n<25; ++n){
		bool add = true;
		if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
		int genlepID = fMT2tree->genlept[n].ID;
		int genlepMID = fMT2tree->genlept[n].MID;
		int genlepGMID = fMT2tree->genlept[n].GMID;
		if(((abs(genlepID)==11)||(abs(genlepID)==13))){
			//include tau decays with good hopes
			if(abs(genlepMID)!=24 && abs(genlepMID)!=15 )     continue;
			if(abs(genlepGMID)!=6 && abs(genlepMID)==24 )     continue;
			if(abs(genlepGMID)!=24&& abs(genlepMID)==15 )     continue;
			if(abs(genlepMID)==24 && genlepMID *genlepID>=0 ) continue;
			if(abs(genlepMID)==15 && genlepMID *genlepID<=0 ) continue;
			if(genlepGMID*genlepID>=0) continue;
			for(unsigned int nn=0; nn<gvv.size(); ++nn){
				if(genlepID  != fMT2tree->genlept[gvv[nn] ].ID  ) continue;
				if(genlepMID != fMT2tree->genlept[gvv[nn] ].MID ) continue;
				if(genlepGMID!= fMT2tree->genlept[gvv[nn] ].GMID) continue;
				if(fMT2tree->genlept[n].lv.DeltaR(fMT2tree->genlept[gvv[nn] ].lv)<0.15) continue;
				add = false;
			}
			if(add) gvv.push_back(n);
		}
	   }
	   for(int n = 0; n<25; ++n){
		bool add = true;
		if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
		int genlepID = fMT2tree->genlept[n].ID;
		int genlepMID = fMT2tree->genlept[n].MID;
		//int genlepGMID = fMT2tree->genlept[n].GMID;
		if(abs(genlepID)==5){
			if(abs(genlepMID)!=6)     continue;
			if(genlepMID*genlepID<=0) continue;
			for(unsigned int nn=0; nn<gvv.size(); ++nn){
				if(genlepID  != fMT2tree->genlept[gvv[nn] ].ID  ) continue;
				if(genlepMID != fMT2tree->genlept[gvv[nn] ].MID ) continue;
				//if(genlepGMID!= fMT2tree->genlept[gvv[nn] ].GMID) continue;
				if(fMT2tree->genlept[n].lv.DeltaR(fMT2tree->genlept[gvv[nn] ].lv)<0.15) continue;
				add = false;
			}
			if(add) gvv.push_back(n);
		}
	   }
	//match l's and b's to OS
	if(gvv.size()>=4){
		for(unsigned int n = 0; n<gvv.size(); ++n){
			int genlepID = fMT2tree->genlept[gvv[n] ].ID;
			if(((abs(genlepID)==11)||(abs(genlepID)==13))&&gl1==-1) gl1 = gvv[n];
			if(((abs(genlepID)==11)||(abs(genlepID)==13))&&gl2==-1){
				if(gl1>=0&&gl1!=gvv[n]){
					if(genlepID*fMT2tree->genlept[gl1 ].ID < 0) gl2 = gvv[n];
				}
			}
		}
		for(unsigned int n = 0; n<gvv.size(); ++n){
			int genlepID = fMT2tree->genlept[gvv[n] ].ID;
			if((abs(genlepID)==5)&&gb1==-1) gb1 = gvv[n];
			if((abs(genlepID)==5)&&gb2==-1){
				if(gb1>=0&&gb1!=gvv[n]){
					if(genlepID*fMT2tree->genlept[gb1 ].ID < 0) gb2 = gvv[n];
				}
			}
		}
	}
	//match l,b together
	if(gl1>=0&&gl2>=0&&gb1>=0&&gb2>=0){
		if(fMT2tree->genlept[gl1].ID*fMT2tree->genlept[gb1].ID>0){//have to swap gb1 and gb2
			int gb2t=gb2; int gb1t = gb1;
			gb1 = gb2t;
			gb2 = gb1t;
		}
	}
	float MT2lbgen(-999.), MT2lbmasslessgen(-999.), MT2bbgen(-999.), MT2llgen(-999.), Mlb1gen(-999.), Mlb2gen(-999.);
	//get MT2 etc
	if(gl1>=0&&gl2>=0&&gb1>=0&&gb2>=0){
		MT2lbgen = fMT2tree->CalcMT2(0., true, fMT2tree->genlept[gl1].lv + fMT2tree->genlept[gb1].lv, fMT2tree->genlept[gl2].lv + fMT2tree->genlept[gb2].lv, fMT2tree->genmet[0]);
		MT2lbmasslessgen = fMT2tree->CalcMT2(0., false, fMT2tree->genlept[gl1].lv + fMT2tree->genlept[gb1].lv, fMT2tree->genlept[gl2].lv + fMT2tree->genlept[gb2].lv, fMT2tree->genmet[0]);
		MT2bbgen = fMT2tree->CalcMT2(0., true,  fMT2tree->genlept[gb1].lv, fMT2tree->genlept[gb2].lv, fMT2tree->genmet[0]);
		MT2llgen = fMT2tree->CalcMT2(0., true,  fMT2tree->genlept[gl1].lv, fMT2tree->genlept[gl2].lv, fMT2tree->genmet[0]);
		Mlb1gen  = (fMT2tree->genlept[gl1].lv+fMT2tree->genlept[gb1].lv).M();
		Mlb2gen  = (fMT2tree->genlept[gl2].lv+fMT2tree->genlept[gb1].lv).M();
	//	cout << "id gen b1/b2/l1/l2 " << gb1<<"/"<<gb2<<"/"<<gl1<<"/"<<gl2<< " and pt " << fMT2tree->genlept[gb1].lv.Pt() <<"/"<< fMT2tree->genlept[gb2].lv.Pt()<<"/"<< fMT2tree->genlept[gl1].lv.Pt() <<"/"<< fMT2tree->genlept[gl2].lv.Pt()<<"/"<<fMT2tree->genmet[0].Pt()<<endl;
	//	cout << "MT2lbgen "<<MT2lbgen<<"  MT2lbmasslessgen " << MT2lbmasslessgen << "  MT2bbgen " << MT2bbgen <<"  MT2llgen " << MT2llgen << "  Mlb1gen " << Mlb1gen << "  Mlb2gen " << endl;
	}
	int trueb1match(-1), trueb2match(-1);//jet matched to b-parton
	float minDRtruthmatching = 9999.;
	for(int n = 0; n<25; ++n){
		if(gb1<0) break;
		if(fMT2tree->jet[n].lv.Pt()<10.) continue;//veto 0
		float dRb1 = fMT2tree->genlept[gb1].lv.DeltaR(fMT2tree->jet[n].lv);
		if(dRb1<0.5 && dRb1<minDRtruthmatching){
			if(trueb1match>=0 && fMT2tree->jet[trueb1match].isPFIDLoose && !(fMT2tree->jet[n].isPFIDLoose)) continue;//matched already a good jet
			if(trueb1match>=0 && fMT2tree->jet[trueb1match].IsBJet(2) && !(fMT2tree->jet[n].IsBJet(2))) continue;//matched already a good b-tagged jet
			trueb1match = n; minDRtruthmatching = dRb1; }
	}
	minDRtruthmatching = 9999.;
	for(int n = 0; n<25; ++n){
		if(gb2<0) break;
		if(fMT2tree->jet[n].lv.Pt()<10.) continue;//veto 0
		float dRb1 = fMT2tree->genlept[gb2].lv.DeltaR(fMT2tree->jet[n].lv);
		if(dRb1<0.5 && dRb1<minDRtruthmatching){
			if(trueb2match>=0 && fMT2tree->jet[trueb2match].isPFIDLoose && !(fMT2tree->jet[n].isPFIDLoose)) continue;//matched already a good jet
			if(trueb2match>=0 && fMT2tree->jet[trueb2match].IsBJet(2) && !(fMT2tree->jet[n].IsBJet(2))) continue;//matched already a good b-tagged jet
			trueb2match = n; minDRtruthmatching = dRb1; }
	}


	//new discriminant implementation
	double b1arr[4], b2arr[4], l1arr[4], l2arr[4], metarr[4];
	b1arr[0] = bjets[ind_l1bxMT2min].Px(); b1arr[1] = bjets[ind_l1bxMT2min].Py();
	b1arr[2] = bjets[ind_l1bxMT2min].Pz(); b1arr[3] = bjets[ind_l1bxMT2min].P(); 
	b2arr[0] = bjets[ind_l2bxMT2min].Px(); b2arr[1] = bjets[ind_l2bxMT2min].Py(); 
	b2arr[2] = bjets[ind_l2bxMT2min].Pz(); b2arr[3] = bjets[ind_l2bxMT2min].P(); 
	l1arr[0]  = l1.Px();  l1arr[1]  = l1.Py();  l1arr[2]  = l1.Pz();  l1arr[3]  = l1.P();
	l2arr[0]  = l2.Px();  l2arr[1]  = l2.Py();  l2arr[2]  = l2.Pz();  l2arr[3]  = l2.P();
	metarr[0] = met.Px(); metarr[1] = met.Py(); metarr[2] = met.Pz(); metarr[3] = met.P();


	double discriminant(-999.99), edm(9999.),q1z(0.),q2z(0.); int status(-99);
	double discriminant1(-999.99);//, discriminant2(-999.99);
	int status1(-99);//, status2(-99);
	double edm1(9999.), q1z1(0.),q2z1(0.);
	int statvec1[3]; 
	bool failedroot(false);
	if(calcD2&&calcD2Minuit){//NOTE: here also Luc's implementation needed, for 4th order equation solving
		sttbLuc->SetNscan(1000.);
		sttbLuc->SetDiscrConv(0.05);
	    	int globalsolution1 = 0;// int globalsolution2 = 0;
		int nsol1 = sttbLuc->SolveEqn4(b1arr, b2arr, l1arr, l2arr, metarr);
		if(nsol1>0) {	discriminant1 = 0.; double qq11[]={0.,0.,0.}, qq22[]={0.,0.,0.}; 
				//get only first solution
				sttbLuc->GetMomSol(0,qq11,qq22); q1z1 = qq11[2]; q2z1 = qq22[2]; ++globalsolution1;
		}
		else{
			discriminant1 = sttb->NumericalMinimization(status1, edm1, q1z1, q2z1, b1arr, b2arr, l1arr, l2arr, metarr);
			if(fabs(q1z1)>=8000.||fabs(q2z1)>=8000.||q1z1==0.||q2z1==0.||edm1>0.01||discriminant1>800.||status1!=0) { 
				q1z=q1z1;q2z=q2z1;
				discriminant = sttb->NumericalMinimization(status, edm, q1z, q2z, b1arr, b2arr, l1arr, l2arr, metarr, "Minuit2", "Migrad");
				if(fabs(q1z)<8000.&&fabs(q2z)<8000.&&q1z!=0.&&q2z!=0.&&edm<=0.01&&status==0){
					discriminant1 = discriminant; status1=status;edm1=edm;q1z1=q1z;q2z1=q2z;
				}
				else if(status1==0&&status!=0){ discriminant1 = discriminant; status1=status;edm1=edm;q1z1=q1z;q2z1=q2z; }
			}//killed Combined and Fumili2 --> Put back in??
			//NOTE: in loop might change sth. here
			if((fabs(q1z1)>=8000.||fabs(q2z1)>=8000.||q1z1==0.||q2z1==0.||edm1>0.01||status1!=0)) ++vetoeddiscrcount; 
			else ++discrcount;
			globsol += globalsolution1;
		}
	}

	if(calcD2&&calcD2Luc){
	    	statvec1[0]=0; statvec1[1]=0; statvec1[2]=0; int nstptot1 = 0; int ndertot1 = 0; int nrottot1 = 0; int nscntot1 = 0;
	    	int globalsolution1 = 0;// int globalsolution2 = 0;
		sttbLuc->SetNscan(1000.);
		sttbLuc->SetDiscrConv(0.05);
		int nsol1 = sttbLuc->SolveEqn4(b1arr, b2arr, l1arr, l2arr, metarr);
		if(nsol1>0) {	discriminant1 = 0.; double qq11[]={0.,0.,0.}, qq22[]={0.,0.,0.}; 
				//get only first solution
				sttbLuc->GetMomSol(0,qq11,qq22); q1z1 = qq11[2]; q2z1 = qq22[2]; ++globalsolution1;
		}
		else {  discriminant1 = sttbLuc->IterDeriv(b1arr, b2arr, l1arr, l2arr, metarr);//luc's code
			q1z1 = sttbLuc->GetXMin(); q2z1 = sttbLuc->GetYMin();
			int stat1 = sttbLuc->GetStatus ();    if (stat1 >= 0) statvec1[stat1]++;
			nstptot1 += sttbLuc->GetNSteps ();    ndertot1 += sttbLuc->GetNDerSteps ();
			nrottot1 += sttbLuc->GetNRotSteps (); nscntot1 += sttbLuc->GetNScanSteps ();
		}
		if(failedroot ) cout << "Luc: D2 " << discriminant1 << " q1z/q2z " <<q1z1<<"/"<<q2z1<<" solution/conv/loop/fail " << globalsolution1<<"/"<<statvec1[0]<<"/"<<statvec1[1]<<"/"<<statvec1[2]<<endl<<endl;
			statvec[0] += statvec1[0]; statvec[1] += statvec1[1]; statvec[2] += statvec1[2];
			nstptot += nstptot1; ndertot += ndertot1; nrottot += nrottot1; nscntot += nscntot1;
			globalsolution += globalsolution1;

	double D2 = discriminant1;//stupid, but for simpleneww use this

	if(debugD2){
		bool plotsomething = false;
		double dd1(-999.9), ddtempA(-99.99), ddtempB(-99.99);
		int stattempB(-99); double edmtempB(9999.);
		double z11(0.), z1tempA(0.), z1tempB(0.);
		double z21(0.), z2tempA(0.), z2tempB(0.);
		double z11_2(0.), z11_3(0.), z11_4(0.);
		double z21_2(0.), z21_3(0.), z21_4(0.);
		int staticvector1[3]; staticvector1[0]=0;staticvector1[1]=0;staticvector1[2]=0;//0:conv,1:loop;3.iter(fail)
		int nsolution1 = sttbLuc->SolveEqn4(b1arr, b2arr, l1arr, l2arr, metarr);
		if(nsolution1>0) { dd1 = 0.; double qq11[]={0.,0.,0.}, qq22[]={0.,0.,0.};
			sttbLuc->GetMomSol(0,qq11,qq22); z11 = qq11[2]; z21 = qq22[2];
			sttbLuc->GetMomSol(1,qq11,qq22); z11_2 = qq11[2]; z21_2 = qq22[2];
			sttbLuc->GetMomSol(2,qq11,qq22); z11_3 = qq11[2]; z21_3 = qq22[2];
			sttbLuc->GetMomSol(3,qq11,qq22); z11_4 = qq11[2]; z21_4 = qq22[2]; }
		ddtempA = sttbLuc->IterDeriv(b1arr, b2arr, l1arr, l2arr, metarr);//luc's code
		z1tempA = sttbLuc->GetXMin(); z2tempA = sttbLuc->GetYMin();
		int staticc = sttbLuc->GetStatus ();
		if (staticc >= 0) staticvector1[staticc]++;
		if(nsolution1>0){
			if(ddtempA>0.3){
				cout << __LINE__ << " 4th order equation solved: D = " << dd1 << " at q1z/q2z = " << z11 << " / " << z21 << " and " << z11_2 << " / " << z21_2 <<  " and " << z11_3 << " / " << z21_3 <<  " and " << z11_4 << " / " << z21_4 <<  endl;
				cout << " But Luc's minimizer found D = " << ddtempA << " at q1z/q2z = " << z1tempA << " / " << z2tempA << " conv/loop/fail = " << staticvector1[0]<<"/"<<staticvector1[1]<<"/"<<staticvector1[2]<< endl;
				plotsomething = true; } }
		ddtempB = sttb->NumericalMinimization(stattempB, edmtempB, z1tempB, z2tempB, b1arr, b2arr, l1arr, l2arr, metarr);
		if(nsolution1>0){
			if(ddtempB>0.3){
				if(ddtempA<=0.3) cout << __LINE__ << " 4th order equation solved: D = " << dd1 << " at q1z/q2z = " << z11 << " / " << z21 << " and " << z11_2 << " / " << z21_2 <<  " and " << z11_3 << " / " << z21_3 <<  " and " << z11_4 << " / " << z21_4 <<  endl;
				cout << " But Minuit minimizer found D = " << ddtempB << " at q1z/q2z = " << z1tempB << " / " << z2tempB << " status/edm = " << stattempB<<"/"<<edmtempB<< endl;
				plotsomething = true;  } }
		if((fabs(ddtempA-ddtempB)>0.03*ddtempA&&ddtempA>0.3&&ddtempB>0.3)){
				if(nsolution1<=0) cout << __LINE__ << " Luc's minimizer found D = " << ddtempA << " at q1z/q2z = " << z1tempA << " / " << z2tempA << " conv/loop/fail = " << staticvector1[0]<<"/"<<staticvector1[1]<<"/"<<staticvector1[2]<< endl;
				if(nsolution1<=0) cout << " But Minuit minimizer found D = " << ddtempB << " at q1z/q2z = " << z1tempB << " / " << z2tempB << " status/edm = " << stattempB<<"/"<<edmtempB<< endl;
				cout << "difference/Luc D " << fabs(ddtempA-ddtempB)/ddtempA << " q1z " << fabs(z1tempA-z1tempB)/z1tempA << " q2z " << fabs(z2tempA-z2tempB)/z2tempA << endl;
				plotsomething = true; }
		if(fabs(z1tempB)>=8000.||fabs(z2tempB)>=8000.||z1tempB==0.||z2tempB==0.||(edmtempB>0.01)||ddtempB>500.||stattempB!=0  || ddtempA>500. || staticvector1[2]>0 ){//non-convergence
			if((!(nsolution1>0&&(ddtempA>0.3||ddtempB>0.3))) || (!(nsolution1<=0&&fabs(ddtempA-ddtempB)>0.03*ddtempA&&ddtempA>0.3&&ddtempB>0.3)) ){ 
			cout << __LINE__ << " interesting event: " << endl;
				cout << " Luc's minimizer found D = " << ddtempA << " at q1z/q2z = " << z1tempA << " / " << z2tempA << " conv/loop/fail = " << staticvector1[0]<<"/"<<staticvector1[1]<<"/"<<staticvector1[2]<< endl;
				cout << " Minuit minimizer found D = " << ddtempB << " at q1z/q2z = " << z1tempB << " / " << z2tempB << " status/edm = " << stattempB<<"/"<<edmtempB<< endl;
				cout << "difference/Luc D " << fabs(ddtempA-ddtempB)/ddtempA << " q1z " << fabs(z1tempA-z1tempB)/z1tempA << " q2z " << fabs(z2tempA-z2tempB)/z2tempA << endl;
				plotsomething = true;
			}
		}
		if(fabs(z1tempB)>=8000.||fabs(z2tempB)>=8000.||z1tempB==0.||z2tempB==0.||edmtempB>0.01||ddtempB>300.||stattempB!=0){
			double z1 = z1tempB; double z2 = z2tempB; double dnew(-9999.); int statnew(99.); double edmnew(9999.);
			dnew = sttb->NumericalMinimization(statnew, edmnew, z1, z2, b1arr, b2arr, l1arr, l2arr, metarr, "Minuit2", "Migrad");
			if((fabs(z1)<8000.&&fabs(z2)<8000.&&z1!=0.&&z2!=0.&&edmnew<=0.01&&statnew==0)||(fabs(dnew-ddtempB)>0.03*ddtempB&&ddtempB>0.3)){
				cout << __LINE__ << " Minuit minimizer found D = " << ddtempB << " at q1z/q2z = " << z1tempB << " / " << z2tempB << " status/edm = " << stattempB<<"/"<<edmtempB<< endl;
				cout << " But Minuit/Migrad new minimizer found D = " << dnew << " at q1z/q2z = " << z1 << " / " << z2 << " status/edm = " << statnew<<"/"<<edmnew<< endl;
				z1tempB = z1; z2tempB = z2; stattempB = statnew; edmtempB = edmnew; ddtempB = dnew;
				plotsomething = true;
			}
		}
		if(fabs(z1tempB)>=8000.||fabs(z2tempB)>=8000.||z1tempB==0.||z2tempB==0.||edmtempB>0.01||ddtempB>300.||stattempB!=0){
			double z1 = z1tempB; double z2 = z2tempB; double dnew(-9999.); int statnew(99.); double edmnew(9999.);
			dnew = sttb->NumericalMinimization(statnew, edmnew, z1, z2, b1arr, b2arr, l1arr, l2arr, metarr, "Minuit2", "Combined");
			if((fabs(z1)<8000.&&fabs(z2)<8000.&&z1!=0.&&z2!=0.&&edmnew<=0.01&&statnew==0)||(fabs(dnew-ddtempB)>0.03*ddtempB&&ddtempB>0.3)){
				cout << " But Minuit/Combined new minimizer found D = " << dnew << " at q1z/q2z = " << z1 << " / " << z2 << " status/edm = " << statnew<<"/"<<edmnew<< endl;
				z1tempB = z1; z2tempB = z2; stattempB = statnew; edmtempB = edmnew; ddtempB = dnew;
				plotsomething = true;
			}
		}
		if(fabs(z1tempB)>=8000.||fabs(z2tempB)>=8000.||z1tempB==0.||z2tempB==0.||edmtempB>0.01||ddtempB>300.||stattempB!=0){
			double z1 = z1tempB; double z2 = z2tempB; double dnew(-9999.); int statnew(99.); double edmnew(9999.);
			dnew = sttb->NumericalMinimization(statnew, edmnew, z1, z2, b1arr, b2arr, l1arr, l2arr, metarr, "Minuit2", "Fumili2");
			if((fabs(z1)<8000.&&fabs(z2)<8000.&&z1!=0.&&z2!=0.&&edmnew<=0.01&&statnew==0)||(fabs(dnew-ddtempB)>0.03*ddtempB&&ddtempB>0.3)){
				cout << " But Minuit/Fumili2 minimizer found D = " << dnew << " at q1z/q2z = " << z1 << " / " << z2 << " status/edm = " << statnew<<"/"<<edmnew<< endl;
				z1tempB = z1; z2tempB = z2; stattempB = statnew; edmtempB = edmnew; ddtempB = dnew;
				plotsomething = true;
			}
		}
		if(plotsomething){
			cout << "Event content:" << endl;
			cout << "b1 PxPyPzP " << b1arr[0] << ", " << b1arr[1] << ", " << b1arr[2] << ", " << b1arr[3] << " Discr " << bjetsdiscriminant[ind_l1bxMT2min] << " index " << ind_l1bxMT2min << endl;
			cout << "b2 PxPyPzP " << b2arr[0] << ", " << b2arr[1] << ", " << b2arr[2] << ", " << b2arr[3] << " Discr " << bjetsdiscriminant[ind_l2bxMT2min] << " index " << ind_l2bxMT2min << endl;
			cout << "l1 PxPyPzP " << l1arr[0] << ", " << l1arr[1] << ", " << l1arr[2] << ", " << l1arr[3] << endl;
			cout << "l2 PxPyPzP " << l2arr[0] << ", " << l2arr[1] << ", " << l2arr[2] << ", " << l2arr[3] << endl;
			cout << "met PxPy   " << metarr[0]<< ", " << metarr[1]<< endl;
			cout <<  endl;
		}
	}

	//DO STUPID CHECKS FLAGS --> kill at a later stage --> need to be fixed for multiple bjets
	//NOTE: PUT THIS LATER WHEN HAVING CORRECT B ETC...
	bool bjet1Pt40(false),bjet1Pt30(false),bjet1Pt20(false), bjet2Pt40(false),bjet2Pt30(false),bjet2Pt20(false), nbjetsge2(false), nbjetseq2(false);
	bool lep1Pt20(false), lep1Pt17(false), lep2Pt10(false), lep2Pt8(false);
	bool MET30(false), MET20(false), MET10(false);
	bool nbjetsge2_30(false), nbjetseq2_30(false), nbjetsge2_20(false), nbjetseq2_20(false);
	if(met.Pt()>30){ MET30 = true; MET20 = true; MET10 = true;}
	if(met.Pt()>20){ MET20 = true; MET10 = true;}
	if(met.Pt()>10){ MET10 = true;}
	if(bjets[ind_l1bxMT2min].Pt()>40) { bjet1Pt40=true;bjet1Pt30=true;bjet1Pt20=true;}
	if(bjets[ind_l1bxMT2min].Pt()>30) { bjet1Pt30=true;bjet1Pt20=true;}
	if(bjets[ind_l1bxMT2min].Pt()>20) bjet1Pt20=true;
	if(bjets[ind_l2bxMT2min].Pt()>40) { bjet2Pt40=true;bjet2Pt30=true;bjet2Pt20=true;}
	if(bjets[ind_l2bxMT2min].Pt()>30) { bjet2Pt30=true;bjet2Pt20=true;}
	if(bjets[ind_l2bxMT2min].Pt()>20) bjet2Pt20=true;
	if(fMT2tree->GetNBtags(2,1.74,40.,2.4,1)==2) nbjetseq2 = true;
	if(fMT2tree->GetNBtags(2,1.74,40.,2.4,1)>=2) nbjetsge2 = true;
	if(fMT2tree->GetNBtags(2,1.74,30.,2.4,1)==2) nbjetseq2_30 = true;
	if(fMT2tree->GetNBtags(2,1.74,30.,2.4,1)>=2) nbjetsge2_30 = true;
	if(fMT2tree->GetNBtags(2,1.74,20.,2.4,1)==2) nbjetseq2_20 = true;
	if(fMT2tree->GetNBtags(2,1.74,20.,2.4,1)>=2) nbjetsge2_20 = true;
	if(l1.Pt()>20){ lep1Pt20 = true; lep1Pt17 = true;}
	if(l1.Pt()>17) lep1Pt17 = true;
	if(l2.Pt()>10){ lep2Pt10 = true; lep2Pt8 = true;}
	if(l2.Pt()>8) lep2Pt8 = true;

	bool all(false),allbutMET30(false), allbutMET20(false), allbutlep1Pt(false),allbutlep2Pt(false),allbutNBJets(false),allbutbjet1Pt40(false), allbutbjet2Pt40(false), allbutbjetPt40(false),allbutbjet1Pt30(false),allbutbjet2Pt30(false), allbutbjetPt30(false);// allbutlepPt(false);
	if(bjet1Pt40&&bjet2Pt40&&nbjetsge2   &&lep1Pt20&&lep2Pt10&&MET30) all = true;
	if(bjet1Pt40&&bjet2Pt40&&nbjetsge2   &&lep1Pt20&&lep2Pt10&&MET20) allbutMET30 = true;
	if(bjet1Pt40&&bjet2Pt40&&nbjetsge2   &&lep1Pt20&&lep2Pt10&&MET10) allbutMET20 = true;
	if(bjet1Pt40&&bjet2Pt40&&nbjetsge2   &&lep1Pt17&&lep2Pt10&&MET30) allbutlep1Pt = true;
	if(bjet1Pt40&&bjet2Pt40&&nbjetsge2   &&lep1Pt20&&lep2Pt8 &&MET30) allbutlep2Pt = true;
	if(bjet1Pt40&&bjet2Pt40&&nbjetseq2   &&lep1Pt20&&lep2Pt10&&MET30) allbutNBJets = true;
	if(bjet1Pt30&&bjet2Pt40&&nbjetsge2_30&&lep1Pt20&&lep2Pt10&&MET30) allbutbjet1Pt40 = true;
	if(bjet1Pt40&&bjet2Pt30&&nbjetsge2_30&&lep1Pt20&&lep2Pt10&&MET30) allbutbjet2Pt40 = true;
	if(bjet1Pt30&&bjet2Pt30&&nbjetsge2_30&&lep1Pt20&&lep2Pt10&&MET30) allbutbjetPt40 = true;
	if(bjet1Pt20&&bjet2Pt40&&nbjetsge2_20&&lep1Pt20&&lep2Pt10&&MET30) allbutbjet1Pt30 = true;
	if(bjet1Pt40&&bjet2Pt20&&nbjetsge2_20&&lep1Pt20&&lep2Pt10&&MET30) allbutbjet2Pt30 = true;
	if(bjet1Pt20&&bjet2Pt20&&nbjetsge2_20&&lep1Pt20&&lep2Pt10&&MET30) allbutbjetPt30 = true;

	if(!(all) && (allbutMET30)){
		if(met.Pt()<20 || met.Pt()>30) cout << "met " << met.Pt() << endl << "b1/b2_Pt40,b1/b2_Pt30,b1/b2_Pt20 " << bjet1Pt40<<"/"<<bjet2Pt40<<","<<bjet1Pt30<<","<<bjet2Pt30<<"/"<<bjet1Pt20<<","<<bjet2Pt20 << " nbjetseq2,ge2/_30/_20 " << nbjetseq2 << ","<<nbjetsge2<< "/"<<nbjetsge2_30<<"/"<<nbjetsge2_20<< " l1/2_Pt20/10,17/8 " << lep1Pt20<<"/"<<lep2Pt10<<","<<lep1Pt17<<"/"<<lep2Pt8<<" met30/20/10 " << MET30<<"/"<<MET20<<"/"<<MET10<<" " << int(true) << endl;
	}

	string hs = string("_") + leptype + string("_") + sampletype;

	int nnbbjjeettss20 = fMT2tree->GetNBtags(2,1.74,20.,2.4,1);
	int nnbbjjeettss30 = fMT2tree->GetNBtags(2,1.74,30.,2.4,1);
	int nnbbjjeettss40 = fMT2tree->GetNBtags(2,1.74,40.,2.4,1);

	if(njets>=0 && fMT2tree->NJetsIDLoose40==    njets ) histos["Mll" + hs]->Fill(Mll, weight);//NOTE: might change this depending on jetthreshold
	if(njets< 0 && fMT2tree->NJetsIDLoose40>=abs(njets)) histos["Mll" + hs]->Fill(Mll, weight);//NOTE: might change this depending on jetthreshold
	if(njets>=0 && fMT2tree->NJetsIDLoose40==    njets  && MT2minV_withMlb>=0 ) histos["MT2min_withMlbcut_woZveto" + hs]->Fill(MT2minV_withMlb, weight);
	if(njets< 0 && fMT2tree->NJetsIDLoose40>=abs(njets) && MT2minV_withMlb>=0 ) histos["MT2min_withMlbcut_woZveto" + hs]->Fill(MT2minV_withMlb, weight);
	if(njets>=0 && fMT2tree->NJetsIDLoose40==   njets   ) histos["MT2min_woZveto" + hs]->Fill(MT2minV, weight);
	if(njets< 0 && fMT2tree->NJetsIDLoose40>=abs(njets) ) histos["MT2min_woZveto" + hs]->Fill(MT2minV, weight);
	if(njets>=0 && fMT2tree->NJetsIDLoose40==    njets  && MT2minV_massless>=0 ) histos["MT2min_massless_woZveto" + hs]->Fill(MT2minV_massless, weight);
	if(njets< 0 && fMT2tree->NJetsIDLoose40>=abs(njets) && MT2minV_massless>=0 ) histos["MT2min_massless_woZveto" + hs]->Fill(MT2minV_massless, weight);
	if(njets>=0 && fMT2tree->NJetsIDLoose40==    njets  && MT2bbmin>=0 ) histos["MT2bb_woZveto" + hs]->Fill(MT2bbmin, weight);
	if(njets< 0 && fMT2tree->NJetsIDLoose40>=abs(njets) && MT2bbmin>=0 ) histos["MT2bb_woZveto" + hs]->Fill(MT2bbmin, weight);
	if(njets>=0 && fMT2tree->NJetsIDLoose40==    njets  && MT2ll>=0 ) histos["MT2ll_woZveto" + hs]->Fill(MT2ll, weight);
	if(njets< 0 && fMT2tree->NJetsIDLoose40>=abs(njets) && MT2ll>=0 ) histos["MT2ll_woZveto" + hs]->Fill(MT2ll, weight);

	if(!(recoedeenZ) &&                   !(recoedmumunZ) && !(recoedemu) ) continue;//require of Zpeakveto

	//STUPID CHECKS
	if(stupidchecks){
		if(Ml1b[ind_l1bxMT2minMlb]<180 && Ml2b[ind_l2bxMT2minMlb]<180){
			bool onecombi(true), allcombi(true);
			for(unsigned int i1 =0; i1<NumBJets; ++i1){
				bool onecombi1(true);
				if(i1==unsigned int(ind_l1bxMT2minMlb)) continue;
				if(Ml1b[i1]<180) onecombi1 = false;
				else             allcombi = false;
				for(unsigned int i2 =0; i2<NumBJets; ++i2){
					if(i2==unsigned int(ind_l1bxMT2minMlb)) continue;
					if(i1==i2) continue;
					if(Ml1b[i2]<180 && onecombi1==false) onecombi = false;//there is another combi possible;
					else if(Ml1b[i2]>=180) allcombi = false;
				}
			}
			if(onecombi) histos["Y1_MT2lbcombi_Mlble180ForOneCombiOnly"+hs]->Fill(MT2minV_withMlb, weight);
			if(allcombi) histos["Y1_MT2lbmin_Mlble180ForBothCombi"+hs]->Fill(MT2minV, weight);//RENAME
		}
		bool bothcombi(true),allcombi(true);
		for(unsigned int i1 =0; i1<NumBJets; ++i1){
			bool bothcombi1(true);
			if(Ml1b[i1]<180) { allcombi = false; bothcombi1 = false; }
			for(unsigned int i2 =0; i2<NumBJets; ++i2){
				if(i2==unsigned int(ind_l1bxMT2minMlb)) continue;
				if(i1==i2) continue;
				if(Ml1b[i2]<180 && bothcombi1==false) bothcombi = false;//there is another combi possible;
				if(Ml1b[i2]<180) allcombi = false;
			}
		}
		if(bothcombi) histos["Y1_MT2lbmin_Mlbge180ForBothCombi"+hs]->Fill(MT2minV, weight);// RENAME
		if(allcombi)  histos["Y1_MT2lbmin_Mlbge180ForEverything"+hs]->Fill(MT2minV, weight);

	if(all){

		histos["Z1_MT2lbmin" + hs]->Fill(MT2minV, weight);
		histos["Z1_MT2lbmin_afterMlbcut" + hs]->Fill(MT2minV_withMlb, weight);
		histos["Z1_MT2lbmin_massless" + hs]->Fill(MT2minV_massless, weight);
		histos["Z1_D2min" + hs]->Fill(D2, weight);
		histos["Z1_MT2ll" + hs]->Fill(MT2ll, weight);
	}
	if(allbutMET30){
		histos["Z1_MT2lbmin_MET20" + hs]->Fill(MT2minV, weight);
		histos["Z1_MT2lbmin_massless_MET20" + hs]->Fill(MT2minV_massless, weight);
		histos["Z1_MT2lbmin_afterMlbcut_MET20" + hs]->Fill(MT2minV_withMlb, weight);
		histos["Z1_D2min_MET20" + hs]->Fill(D2, weight);
		histos["Z1_MT2ll_MET20" + hs]->Fill(MT2ll, weight);
		histos["Z1_MET_MET20" + hs]->Fill(met.Pt(), weight);
	} if(allbutMET20){
		if(met.Pt()<10) cout << "Err met " << met.Pt() << endl;
		if(bjets[ind_l1bxMT2min].Pt()<40 || bjets[ind_l2bxMT2min].Pt()<40 || l1.Pt()<20 || l2.Pt()<10) cout << __LINE__ << " and nbjets20/30/40 " << nnbbjjeettss20 <<"/"<<nnbbjjeettss30<<"/"<<nnbbjjeettss40<< " b1/b2 " << bjets[ind_l1bxMT2min].Pt()<<"/"<<bjets[ind_l2bxMT2min].Pt() << " l1/l2 " << l1.Pt() <<"/"<< l2.Pt() << " met " << met.Pt() << " should be all but met>10 " << endl;
		histos["Z1_MT2lbmin_afterMlbcut_MET10" + hs]->Fill(MT2minV_withMlb, weight);
		histos["Z1_MT2lbmin_MET10" + hs]->Fill(MT2minV, weight);
		histos["Z1_MT2lbmin_massless_MET10" + hs]->Fill(MT2minV_massless, weight);
		histos["Z1_D2min_MET10" + hs]->Fill(D2, weight);
		histos["Z1_MT2ll_MET10" + hs]->Fill(MT2ll, weight);
		histos["Z1_MET_MET10" + hs]->Fill(met.Pt(), weight);
	} if(allbutbjetPt30){
		if(bjets[ind_l1bxMT2min].Pt()<20 ||bjets[ind_l2bxMT2min].Pt()<20) cout << "Err b1 " << bjets[ind_l1bxMT2min].Pt() << " b2 " << bjets[ind_l2bxMT2min].Pt() << endl;
		if(met.Pt()<30 || l1.Pt()<20 || l2.Pt()<10) cout << __LINE__ << " and nbjets20/30/40 " << nnbbjjeettss20 <<"/"<<nnbbjjeettss30<<"/"<<nnbbjjeettss40<< " b1/b2 " << bjets[ind_l1bxMT2min].Pt()<<"/"<<bjets[ind_l2bxMT2min].Pt() << " l1/l2 " << l1.Pt() <<"/"<< l2.Pt() << " met " << met.Pt() << " should be all but b>20 " << endl;
		histos["Z1_MT2lbmin_afterMlbcut_looserBJet_20" + hs]->Fill(MT2minV_withMlb, weight);
		histos["Z1_MT2lbmin_looserBJet_20" + hs]->Fill(MT2minV, weight);
		histos["Z1_MT2lbmin_massless_looserBJet_20" + hs]->Fill(MT2minV_massless, weight);
		histos["Z1_D2min_looserBJet_20" + hs]->Fill(D2, weight);
		histos["Z1_MT2ll_looserBJet_20" + hs]->Fill(MT2ll, weight);
		histos["Z1_BPt_looserBjet_20" + hs]->Fill(bjets[ind_l1bxMT2min].Pt(), weight);
		histos["Z1_BPt_looserBjet_20" + hs]->Fill(bjets[ind_l2bxMT2min].Pt(), weight);
	} if(allbutbjetPt40){
		if(bjets[ind_l1bxMT2min].Pt()<30 ||bjets[ind_l2bxMT2min].Pt()<30) cout << "Err b1 " << bjets[ind_l1bxMT2min].Pt() << " b2 " << bjets[ind_l2bxMT2min].Pt() << endl;
		if(met.Pt()<30 || l1.Pt()<20 || l2.Pt()<10) cout << __LINE__ << " and nbjets20/30/40 " << nnbbjjeettss20 <<"/"<<nnbbjjeettss30<<"/"<<nnbbjjeettss40<< " b1/b2 " << bjets[ind_l1bxMT2min].Pt()<<"/"<<bjets[ind_l2bxMT2min].Pt() << " l1/l2 " << l1.Pt() <<"/"<< l2.Pt() << " met " << met.Pt() << " should be all but b>30 " << endl;
		histos["Z1_MT2lbmin_afterMlbcut_looserBJet_30" + hs]->Fill(MT2minV_withMlb, weight);
		histos["Z1_MT2lbmin_looserBJet_30" + hs]->Fill(MT2minV, weight);
		histos["Z1_MT2lbmin_massless_looserBJet_30" + hs]->Fill(MT2minV_massless, weight);
		histos["Z1_D2min_looserBJet_30" + hs]->Fill(D2, weight);
		histos["Z1_MT2ll_looserBJet_30" + hs]->Fill(MT2ll, weight);
		histos["Z1_BPt_looserBjet_30" + hs]->Fill(bjets[ind_l1bxMT2min].Pt(), weight);
		histos["Z1_BPt_looserBjet_30" + hs]->Fill(bjets[ind_l2bxMT2min].Pt(), weight);
	} if(allbutNBJets){//this should be soon absolete anyway
		if(bjets[ind_l1bxMT2min].Pt()<40 || bjets[ind_l2bxMT2min].Pt()<40 || l1.Pt()<20 || l2.Pt()<10 || met.Pt()<30) cout << __LINE__ << " and nbjets20/30/40 " << nnbbjjeettss20 <<"/"<<nnbbjjeettss30<<"/"<<nnbbjjeettss40<< " b1/b2 " << bjets[ind_l1bxMT2min].Pt()<<"/"<<bjets[ind_l2bxMT2min].Pt() << " l1/l2 " << l1.Pt() <<"/"<< l2.Pt() << " met " << met.Pt() << " should be all but nbjets>=2 " << endl;
		histos["Z1_NBJets_tighterNBJets" + hs]->Fill(NumBJets, weight);
		histos["Z1_MT2lbmin_afterMlbcut_tighterNBJets" + hs]->Fill(MT2minV_withMlb, weight);
		histos["Z1_MT2lbmin_tighterNBJets" + hs]->Fill(MT2minV, weight);
		histos["Z1_MT2lbmin_massless_tighterNBJets" + hs]->Fill(MT2minV_massless, weight);
		histos["Z1_D2min_tighterNBJets" + hs]->Fill(D2, weight);
		histos["Z1_MT2ll_tighterNBJets" + hs]->Fill(MT2ll, weight);
	} if(allbutlep1Pt){
		if(l1.Pt()<17) cout << "Err l1 " << l1.Pt() << endl;
		if(bjets[ind_l2bxMT2min].Pt()<40 || met.Pt()<30 || bjets[ind_l1bxMT2min].Pt()<40 || l2.Pt()<10) cout << __LINE__ << " and nbjets20/30/40 " << nnbbjjeettss20 <<"/"<<nnbbjjeettss30<<"/"<<nnbbjjeettss40<< " b1/b2 " << bjets[ind_l1bxMT2min].Pt()<<"/"<<bjets[ind_l2bxMT2min].Pt() << " l1/l2 " << l1.Pt() <<"/"<< l2.Pt() << " met " << met.Pt() << " should be all but l1>17 " << endl;
		histos["Z1_LepPt_looser1stLep" + hs]->Fill(l1.Pt(), weight);
		histos["Z1_MT2lbmin_afterMlbcut_looser1stLep" + hs]->Fill(MT2minV_withMlb, weight);
		histos["Z1_MT2lbmin_looser1stLep" + hs]->Fill(MT2minV, weight);
		histos["Z1_MT2lbmin_massless_looser1stLep" + hs]->Fill(MT2minV_massless, weight);
		histos["Z1_D2min_looser1stLep" + hs]->Fill(D2, weight);
		histos["Z1_MT2ll_looser1stLep" + hs]->Fill(MT2ll, weight);
	} if(allbutlep2Pt){
		if(l2.Pt()<8) cout << "Err l2 " << l2.Pt() << endl;
		if(bjets[ind_l2bxMT2min].Pt()<40 || met.Pt()<30 || bjets[ind_l1bxMT2min].Pt()<40 || l1.Pt()<20) cout << __LINE__ << " and nbjets20/30/40 " << nnbbjjeettss20 <<"/"<<nnbbjjeettss30<<"/"<<nnbbjjeettss40<< " b1/b2 " << bjets[ind_l1bxMT2min].Pt()<<"/"<<bjets[ind_l2bxMT2min].Pt() << " l1/l2 " << l1.Pt() <<"/"<< l2.Pt() << " met " << met.Pt() << " should be all but l2>8 " << endl;
		histos["Z1_LepPt_looser2ndLep" + hs]->Fill(l2.Pt(), weight);
		histos["Z1_MT2lbmin_afterMlbcut_looser2ndLep" + hs]->Fill(MT2minV_withMlb, weight);
		histos["Z1_MT2lbmin_massless_looser2ndLep" + hs]->Fill(MT2minV_massless, weight);
		histos["Z1_MT2lbmin_looser2ndLep" + hs]->Fill(MT2minV, weight);
		histos["Z1_D2min_looser2ndLep" + hs]->Fill(D2, weight);
		histos["Z1_MT2ll_looser2ndLep" + hs]->Fill(MT2ll, weight);
	}
	}//stupidchecks


	//CHECKS BORDERS -- want to lower it - do this later to get not to manny bugs in the code
	if(MT2ll>=100.) histos["X1_D2min_MT2llge100" + hs]->Fill(D2, weight);
	if(MT2ll>=100. && Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180.) histos["X1_D2min_MT2llge100_Mlble180" + hs]->Fill(D2, weight);
	if(MT2minV>=225.) histos["X1_D2min_MT2lbminge225" + hs]->Fill(D2, weight);
	if(MT2minV>=225. && Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180.) histos["X1_D2min_MT2lbminge225_Mlble180" + hs]->Fill(D2, weight);
	if(Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180.) histos["X1_D2min_Mlble180" + hs]->Fill(D2, weight);
	if(MT2minV>=225. && MT2ll>=100.) histos["X1_D2min_MT2lbminge225_MT2llge100" + hs]->Fill(D2, weight);
	if(MT2minV>=225. && MT2ll>=100. && Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180.) histos["X1_D2min_MT2lbminge225_MT2llge100_Mlble180" + hs]->Fill(D2, weight);
	histos["X1_D2min" + hs]->Fill(D2, weight);
	//here it is a bit stupid with Mlb, but use the ones for MT2minV
	if(D2>50.) histos["X1_MT2ll_D2minge50" + hs]->Fill(MT2ll, weight);
	if(D2>50. && Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180.) histos["X1_MT2ll_D2minge50_Mlble180" + hs]->Fill(MT2ll, weight);
	if(MT2minV>=225.) histos["X1_MT2ll_MT2lbminge225" + hs]->Fill(MT2ll, weight);
	if(MT2minV>=225. && Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180.) histos["X1_MT2ll_MT2lbminge225_Mlble180" + hs]->Fill(MT2ll, weight);
	if(Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180.) histos["X1_MT2ll_Mlble180" + hs]->Fill(MT2ll, weight); 
	if(D2>50. && MT2minV>=225.) histos["X1_MT2ll_D2minge50_MT2lbminge225" + hs]->Fill(MT2ll, weight);
	if(D2>50. && MT2minV>=225. && Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180.) histos["X1_MT2ll_D2minge50_MT2lbminge225_Mlble180" + hs]->Fill(MT2ll, weight);
	histos["X1_MT2ll" + hs]->Fill(MT2ll, weight);
	if(MT2ll>=100.) histos["X1_MT2lbmin_MT2llge100" + hs]->Fill(MT2minV, weight);
	if(MT2ll>=100. && Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180.) histos["X1_MT2lbmin_MT2llge100_Mlble180" + hs]->Fill(MT2minV, weight);
	if(D2>50.) histos["X1_MT2lbmin_D2minge50" + hs]->Fill(MT2minV, weight);
	if(D2>50. && Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180.)  histos["X1_MT2lbmin_D2minge50_Mlble180" + hs]->Fill(MT2minV, weight);
	if(Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180.) histos["X1_MT2lbmin_Mlble180" + hs]->Fill(MT2minV, weight);
	if(D2>50.&&MT2ll>=100.) histos["X1_MT2lbmin_D2minge50_MT2llge100" + hs]->Fill(MT2minV, weight);
	if(D2>50.&&MT2ll>=100. && Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180.) histos["X1_MT2lbmin_D2minge50_MT2llge100_Mlble180" + hs]->Fill(MT2minV, weight);
	histos["X1_MT2lbmin" + hs]->Fill(MT2minV, weight);
	//this is ill-defined somehow - use the two Mlb which make MT2minV
	if(MT2minV>=225.) histos["X1_Mlb_MT2lbminge225" + hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
	if(MT2minV>=225.) histos["X1_Mlb_MT2lbminge225" + hs]->Fill(Ml2b[ind_l2bxMT2min], weight);
	if(D2>50.) histos["X1_Mlb_D2minge50" + hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
	if(D2>50.) histos["X1_Mlb_D2minge50" + hs]->Fill(Ml2b[ind_l2bxMT2min], weight);
	if(MT2ll>=100.) histos["X1_Mlb_MT2llge100" + hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
	if(MT2ll>=100.) histos["X1_Mlb_MT2llge100" + hs]->Fill(Ml2b[ind_l2bxMT2min], weight);
	if(MT2minV>=225.&&D2>50.) histos["X1_Mlb_D2minge50_MT2lbminge225" + hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
	if(MT2minV>=225.&&D2>50.) histos["X1_Mlb_D2minge50_MT2lbminge225" + hs]->Fill(Ml2b[ind_l2bxMT2min], weight);
	if(MT2minV>=225.&&MT2ll>=100.) histos["X1_Mlb_MT2llge100_MT2lbminge225" + hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
	if(MT2minV>=225.&&MT2ll>=100.) histos["X1_Mlb_MT2llge100_MT2lbminge225" + hs]->Fill(Ml2b[ind_l2bxMT2min], weight);
	if(MT2ll>=100.&&D2>50.) histos["X1_Mlb_D2minge50_MT2llge100" + hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
	if(MT2ll>=100.&&D2>50.) histos["X1_Mlb_D2minge50_MT2llge100" + hs]->Fill(Ml2b[ind_l2bxMT2min], weight);
	if(MT2minV>=225.&&MT2ll>=100.&&D2>50.) histos["X1_Mlb_D2minge50_MT2llge100_MT2lbminge225" + hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
	if(MT2minV>=225.&&MT2ll>=100.&&D2>50.) histos["X1_Mlb_D2minge50_MT2llge100_MT2lbminge225" + hs]->Fill(Ml2b[ind_l2bxMT2min], weight);
	histos["X1_Mlb"+hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
	histos["X1_Mlb"+hs]->Fill(Ml2b[ind_l2bxMT2min], weight);
	histos2["X2_D2min_vs_MT2ll" + hs]->Fill(D2, MT2ll, weight);
	histos2["X2_D2min_vs_MT2lbmin" + hs]->Fill(D2, MT2minV, weight);
	histos2["X2_D2min_vs_MT2lbmin_afterMlb" + hs]->Fill(D2, MT2minV_withMlb, weight);//new
	histos2["X2_D2min_vs_Mlb" + hs]->Fill(D2, Ml1b[ind_l1bxMT2min], weight);
	histos2["X2_D2min_vs_Mlb" + hs]->Fill(D2, Ml2b[ind_l2bxMT2min], weight);
	histos2["X2_MT2ll_vs_MT2lbmin" + hs]->Fill(MT2ll, MT2minV, weight);
	histos2["X2_MT2ll_vs_MT2lbmin_afterMlb" + hs]->Fill(MT2ll, MT2minV_withMlb, weight);//new
	histos2["X2_MT2ll_vs_Mlb" + hs]->Fill(D2, Ml1b[ind_l1bxMT2min], weight);
	histos2["X2_MT2ll_vs_Mlb" + hs]->Fill(D2, Ml2b[ind_l2bxMT2min], weight);
	histos2["X2_MT2lbmin_vs_Mlb" + hs]->Fill(MT2minV, Ml1b[ind_l1bxMT2min], weight);
	histos2["X2_MT2lbmin_vs_Mlb" + hs]->Fill(MT2minV, Ml2b[ind_l2bxMT2min], weight);
	histos2["X2_MT2lbmin_afterMlb_vs_Mlb" + hs]->Fill(MT2minV_withMlb, Ml1b[ind_l1bxMT2min], weight);
	histos2["X2_MT2lbmin_afterMlb_vs_Mlb" + hs]->Fill(MT2minV_withMlb, Ml2b[ind_l2bxMT2min], weight);

	bool truel1notb1(false), truel2notb2(false);//need to be defined outside
	//now real plots
	histos["NJetsID20" + hs]->Fill(njetsID20, weight);
	histos["NJetsID30" + hs]->Fill(njetsID30, weight);
	histos["NJetsID40" + hs]->Fill(njetsID40, weight);
	histos["NJetsID50" + hs]->Fill(njetsID50, weight);
	//NOTE: THESE TWO LINES NEEDS TO BE CROSS CHECKED: jetthreshold
	if(truettbarEE){
		if(genMl1b1>=0.) histos[(string)"GenMlb_EE_"+sampletype]->Fill(genMl1b1, weight);
		if(genMl2b2>=0.) histos[(string)"GenMlb_EE_"+sampletype]->Fill(genMl2b2, weight);
		if(genMT2>=0.) histos[(string)"GenMT2_EE_"+sampletype]->Fill(genMT2, weight);
		if(genM1>=0.) histos[(string)"GenMtop_EE_"+sampletype]->Fill(genM1, weight);
		if(genM2>=0.) histos[(string)"GenMtop_EE_"+sampletype]->Fill(genM1, weight);
	} if(truettbarEMu){
		if(genMl1b1>=0.) histos[(string)"GenMlb_EMu_"+sampletype]->Fill(genMl1b1, weight);
		if(genMl2b2>=0.) histos[(string)"GenMlb_EMu_"+sampletype]->Fill(genMl2b2, weight);
		if(genMT2>=0.) histos[(string)"GenMT2_EMu_"+sampletype]->Fill(genMT2, weight);
		if(genM1>=0.) histos[(string)"GenMtop_EMu_"+sampletype]->Fill(genM1, weight);
		if(genM2>=0.) histos[(string)"GenMtop_EMu_"+sampletype]->Fill(genM1, weight);
	} if(truettbarMuMu){
		if(genMl1b1>=0.) histos[(string)"GenMlb_MuMu_"+sampletype]->Fill(genMl1b1, weight);
		if(genMl2b2>=0.) histos[(string)"GenMlb_MuMu_"+sampletype]->Fill(genMl2b2, weight);
		if(genMT2>=0.) histos[(string)"GenMT2_MuMu_"+sampletype]->Fill(genMT2, weight);
		if(genM1>=0.) histos[(string)"GenMtop_MuMu_"+sampletype]->Fill(genM1, weight);
		if(genM2>=0.) histos[(string)"GenMtop_MuMu_"+sampletype]->Fill(genM1, weight);
	} if(truettbarEE){
		if(genMl1b1_noPtcut>=0.) histos[(string)"GenMlb_noPtcut_EE_"+sampletype]->Fill(genMl1b1_noPtcut, weight);
		if(genMl2b2_noPtcut>=0.) histos[(string)"GenMlb_noPtcut_EE_"+sampletype]->Fill(genMl2b2_noPtcut, weight);
		if(genMT2_noPtcut>=0.) histos[(string)"GenMT2_noPtcut_EE_"+sampletype]->Fill(genMT2_noPtcut, weight);
		if(genM1_noPtcut>=0.) histos[(string)"GenMtop_noPtcut_EE_"+sampletype]->Fill(genM1_noPtcut, weight);
		if(genM2_noPtcut>=0.) histos[(string)"GenMtop_noPtcut_EE_"+sampletype]->Fill(genM1_noPtcut, weight);
	} if(truettbarEMu){
		if(genMl1b1_noPtcut>=0.) histos[(string)"GenMlb_noPtcut_EMu_"+sampletype]->Fill(genMl1b1_noPtcut, weight);
		if(genMl2b2_noPtcut>=0.) histos[(string)"GenMlb_noPtcut_EMu_"+sampletype]->Fill(genMl2b2_noPtcut, weight);
		if(genMT2_noPtcut>=0.) histos[(string)"GenMT2_noPtcut_EMu_"+sampletype]->Fill(genMT2_noPtcut, weight);
		if(genM1_noPtcut>=0.) histos[(string)"GenMtop_noPtcut_EMu_"+sampletype]->Fill(genM1_noPtcut, weight);
		if(genM2_noPtcut>=0.) histos[(string)"GenMtop_noPtcut_EMu_"+sampletype]->Fill(genM1_noPtcut, weight);
	} if(truettbarMuMu){
		if(genMl1b1_noPtcut>=0.) histos[(string)"GenMlb_noPtcut_MuMu_"+sampletype]->Fill(genMl1b1_noPtcut, weight);
		if(genMl2b2_noPtcut>=0.) histos[(string)"GenMlb_noPtcut_MuMu_"+sampletype]->Fill(genMl2b2_noPtcut, weight);
		if(genMT2_noPtcut>=0.) histos[(string)"GenMT2_noPtcut_MuMu_"+sampletype]->Fill(genMT2_noPtcut, weight);
		if(genM1_noPtcut>=0.) histos[(string)"GenMtop_noPtcut_MuMu_"+sampletype]->Fill(genM1_noPtcut, weight);
		if(genM2_noPtcut>=0.) histos[(string)"GenMtop_noPtcut_MuMu_"+sampletype]->Fill(genM1_noPtcut, weight);
	}
	//if(njets>=0 && fMT2tree->NJetsIDLoose40!=njets) continue;//NOTE: UNDO THIS
	//if(njets< 0 && fMT2tree->NJetsIDLoose40< abs(njets)) continue;
	if(met.Pt()<30) continue;
	if(l1.Pt()<20) continue;
	if(l2.Pt()<10) continue;
	if(bjets[ind_l1bxMT2min].Pt()<40) continue;
	if(bjets[ind_l2bxMT2min].Pt()<40) continue;
//new 14/06/2012
	if(abs(bjetsflavour[ind_l1bxMT2min])==5&&abs(bjetsflavour[ind_l2bxMT2min])==5)
		histos["MT2lb_bothbtags_truebjets"+hs]->Fill(MT2minV, weight);
	else if((abs(bjetsflavour[ind_l1bxMT2min])!=5&&abs(bjetsflavour[ind_l2bxMT2min])==5) || (abs(bjetsflavour[ind_l1bxMT2min])==5&&abs(bjetsflavour[ind_l2bxMT2min])!=5))
		histos["MT2lb_onebtag_truebjet"+hs]->Fill(MT2minV, weight);
	else if(abs(bjetsflavour[ind_l1bxMT2min])!=5&&abs(bjetsflavour[ind_l2bxMT2min])!=5)
		histos["MT2lb_nonebtags_truebjets"+hs]->Fill(MT2minV, weight);
	if(matchedgenp_tob1>=0 && matchedgenp_tob2>=0 && matchedgenp_tob1!=matchedgenp_tob2)//not match to same b
		histos["MT2lb_bothbtags_matchedtobfromtop"+hs]->Fill(MT2minV, weight);
	else if((matchedgenp_tob1==-1 && matchedgenp_tob2>=0)||(matchedgenp_tob1>=0 && matchedgenp_tob2==10))
		histos["MT2lb_onebtag_matchedtobfromtop"+hs]->Fill(MT2minV, weight);
	else if(matchedgenp_tob1==-1 && matchedgenp_tob2==-1)
		histos["MT2lb_nonebtags_matchedtobfromtop"+hs]->Fill(MT2minV, weight);
	if(MT2minV>200.){
		if(matchedgenp_tob1>=0)
			histos["RecoPt_over_GenPartonPt_truebfromtop_MT2lb_gt_200"+hs]->Fill(bjets[ind_l1bxMT2min].Pt()/fMT2tree->genlept[matchedgenp_tob1].lv.Pt(), weight);
		if(matchedgenp_tob2>=0)
			histos["RecoPt_over_GenPartonPt_truebfromtop_MT2lb_gt_200"+hs]->Fill(bjets[ind_l2bxMT2min].Pt()/fMT2tree->genlept[matchedgenp_tob2].lv.Pt(), weight);
		if(matchedgenp_tob1>=0)
			histos["DeltaPt_bjet_matchedtrueparton_MT2lb_gt_200"+hs]->Fill(fabs(bjets[ind_l1bxMT2min].Pt()-fMT2tree->genlept[matchedgenp_tob1].lv.Pt()), weight);
		if(matchedgenp_tob2>=0)
			histos["DeltaPt_bjet_matchedtrueparton_MT2lb_gt_200"+hs]->Fill(fabs(bjets[ind_l2bxMT2min].Pt()-fMT2tree->genlept[matchedgenp_tob2].lv.Pt()), weight);
	}
	else {
		if(matchedgenp_tob1>=0)
			histos["RecoPt_over_GenPartonPt_truebfromtop_MT2lb_lt_200"+hs]->Fill(bjets[ind_l1bxMT2min].Pt()/fMT2tree->genlept[matchedgenp_tob1].lv.Pt(), weight);
		if(matchedgenp_tob2>=0)
			histos["RecoPt_over_GenPartonPt_truebfromtop_MT2lb_lt_200"+hs]->Fill(bjets[ind_l2bxMT2min].Pt()/fMT2tree->genlept[matchedgenp_tob2].lv.Pt(), weight);
		if(matchedgenp_tob1>=0)
			histos["DeltaPt_bjet_matchedtrueparton_MT2lb_lt_200"+hs]->Fill(fabs(bjets[ind_l1bxMT2min].Pt()-fMT2tree->genlept[matchedgenp_tob1].lv.Pt()), weight);
		if(matchedgenp_tob2>=0)
			histos["DeltaPt_bjet_matchedtrueparton_MT2lb_lt_200"+hs]->Fill(fabs(bjets[ind_l2bxMT2min].Pt()-fMT2tree->genlept[matchedgenp_tob2].lv.Pt()), weight);
	}
	if(abs(bjetsflavour[ind_l1bxMT2min])==5&&abs(bjetsflavour[ind_l2bxMT2min])==5){
		histos["Mlb_bothbtags_truebjets"+hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
		histos["Mlb_bothbtags_truebjets"+hs]->Fill(Ml2b[ind_l2bxMT2min], weight);}
	else if((abs(bjetsflavour[ind_l1bxMT2min])!=5&&abs(bjetsflavour[ind_l2bxMT2min])==5) || (abs(bjetsflavour[ind_l1bxMT2min])==5&&abs(bjetsflavour[ind_l2bxMT2min])!=5)){
		histos["Mlb_onebtag_truebjet"+hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
		histos["Mlb_onebtag_truebjet"+hs]->Fill(Ml2b[ind_l2bxMT2min], weight);}
	else if(abs(bjetsflavour[ind_l1bxMT2min])!=5&&abs(bjetsflavour[ind_l2bxMT2min])!=5){
		histos["Mlb_nonebtags_truebjets"+hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
		histos["Mlb_nonebtags_truebjets"+hs]->Fill(Ml2b[ind_l2bxMT2min], weight);}
	if(matchedgenp_tob1>=0 && matchedgenp_tob2>=0 && matchedgenp_tob1!=matchedgenp_tob2){//not match to same b
		histos["Mlb_bothbtags_matchedtobfromtop"+hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
		histos["Mlb_bothbtags_matchedtobfromtop"+hs]->Fill(Ml2b[ind_l2bxMT2min], weight);}
	else if((matchedgenp_tob1==-1 && matchedgenp_tob2>=0)||(matchedgenp_tob1>=0 && matchedgenp_tob2==-1)){
		histos["Mlb_onebtag_matchedtobfromtop"+hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
		histos["Mlb_onebtag_matchedtobfromtop"+hs]->Fill(Ml2b[ind_l2bxMT2min], weight);}
	else if(matchedgenp_tob1==-1 && matchedgenp_tob2==-1){
		histos["Mlb_nonebtags_matchedtobfromtop"+hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
		histos["Mlb_nonebtags_matchedtobfromtop"+hs]->Fill(Ml2b[ind_l2bxMT2min], weight);}
	if(Ml1b[ind_l1bxMT2min]>200.){
		if(matchedgenp_tob1>=0)
			histos["RecoPt_over_GenPartonPt_truebfromtop_Mlb_gt_200"+hs]->Fill(bjets[ind_l1bxMT2min].Pt()/fMT2tree->genlept[matchedgenp_tob1].lv.Pt(), weight);
	}
	else {
		if(matchedgenp_tob1>=0)
			histos["RecoPt_over_GenPartonPt_truebfromtop_Mlb_lt_200"+hs]->Fill(bjets[ind_l1bxMT2min].Pt()/fMT2tree->genlept[matchedgenp_tob1].lv.Pt(), weight);
	}
	if(Ml2b[ind_l2bxMT2min]>200.){
		if(matchedgenp_tob2>=0)
			histos["RecoPt_over_GenPartonPt_truebfromtop_Mlb_gt_200"+hs]->Fill(bjets[ind_l2bxMT2min].Pt()/fMT2tree->genlept[matchedgenp_tob2].lv.Pt(), weight);
	}
	else {
		if(matchedgenp_tob2>=0)
			histos["RecoPt_over_GenPartonPt_truebfromtop_Mlb_lt_200"+hs]->Fill(bjets[ind_l2bxMT2min].Pt()/fMT2tree->genlept[matchedgenp_tob2].lv.Pt(), weight);
	}
	if(MT2minV>200.){
		histos["BEta_MT2lb_gt_200"+hs]->Fill(bjets[ind_l1bxMT2min].Eta(),weight);
		histos["BEta_MT2lb_gt_200"+hs]->Fill(bjets[ind_l2bxMT2min].Eta(),weight);
		histos["BPt_MT2lb_gt_200"+hs]->Fill(bjets[ind_l1bxMT2min].Pt(),weight);
		histos["BPt_MT2lb_gt_200"+hs]->Fill(bjets[ind_l2bxMT2min].Pt(),weight);
		histos["BDiscr_MT2lb_gt_200"+hs]->Fill(bjetsdiscriminant[ind_l1bxMT2min],weight);
		histos["BDiscr_MT2lb_gt_200"+hs]->Fill(bjetsdiscriminant[ind_l2bxMT2min],weight);
		if(gb1>=0){
			histos["TrueBfromtopEta_MT2lb_gt_200"+hs]->Fill(fMT2tree->genlept[gb1].lv.Eta(),weight);
			histos["TrueBfromtopPt_MT2lb_gt_200"+hs]->Fill(fMT2tree->genlept[gb1].lv.Pt(),weight);
		} if(gb2>=0){
			if(gb1>=0&&gl2>=0&&gl1>=0){
				histos["AltMT2lb_MT2lb_gt_200_usingtruebpartons"+hs]->Fill(MT2lbgen, weight);
			}
			histos["TrueBfromtopEta_MT2lb_gt_200"+hs]->Fill(fMT2tree->genlept[gb2].lv.Eta(),weight);
			histos["TrueBfromtopPt_MT2lb_gt_200"+hs]->Fill(fMT2tree->genlept[gb2].lv.Pt(),weight);
		}
		if(fMT2tree->genmet[0].Pt()>0.001){
			histos["PFMEToverGenMET_MT2lb_gt_200"+hs]->Fill(met.Pt()/fMT2tree->genmet[0].Pt(),weight);
			histos["DeltaGenMET_PFMET_Phi_MT2lb_gt_200"+hs]->Fill(fabs(met.DeltaPhi(fMT2tree->genmet[0])),weight);
		}
		if(gl1>=0&&gl2>=0&&gb1>=0&&gb2>=0){
			if(sqrt(pow(l1.DeltaR(fMT2tree->genlept[gl1].lv),2)+pow(l2.DeltaR(fMT2tree->genlept[gl2].lv),2))<sqrt(pow(l2.DeltaR(fMT2tree->genlept[gl1].lv),2)+pow(l1.DeltaR(fMT2tree->genlept[gl2].lv),2))){
				//l1 matched to genl1, l2 matched to genl2 better than other way round
				histos["DeltaPt_TrueBparton_selectedjet_MT2lb_gt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb2].lv.Pt()-bjets[ind_l2bxMT2min].Pt()),weight);
				histos["DeltaPt_TrueBparton_selectedjet_MT2lb_gt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb1].lv.Pt()-bjets[ind_l1bxMT2min].Pt()),weight);
				histos["DeltaEta_TrueBparton_selectedjet_MT2lb_gt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb2].lv.Eta()-bjets[ind_l2bxMT2min].Eta()),weight);
				histos["DeltaEta_TrueBparton_selectedjet_MT2lb_gt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb1].lv.Eta()-bjets[ind_l1bxMT2min].Eta()),weight);
				histos["DeltaPhi_TrueBparton_selectedjet_MT2lb_gt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb2].lv.DeltaPhi(bjets[ind_l2bxMT2min])),weight);
				histos["DeltaPhi_TrueBparton_selectedjet_MT2lb_gt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb1].lv.DeltaPhi(bjets[ind_l1bxMT2min])),weight);
			} else{
				//l1 matched to genl2, l2 matched to genl1 better than other way round
				histos["DeltaPt_TrueBparton_selectedjet_MT2lb_gt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb1].lv.Pt()-bjets[ind_l2bxMT2min].Pt()),weight);
				histos["DeltaPt_TrueBparton_selectedjet_MT2lb_gt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb2].lv.Pt()-bjets[ind_l1bxMT2min].Pt()),weight);
				histos["DeltaEta_TrueBparton_selectedjet_MT2lb_gt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb1].lv.Eta()-bjets[ind_l2bxMT2min].Eta()),weight);
				histos["DeltaEta_TrueBparton_selectedjet_MT2lb_gt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb2].lv.Eta()-bjets[ind_l1bxMT2min].Eta()),weight);
				histos["DeltaPhi_TrueBparton_selectedjet_MT2lb_gt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb1].lv.DeltaPhi(bjets[ind_l2bxMT2min])),weight);
				histos["DeltaPhi_TrueBparton_selectedjet_MT2lb_gt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb2].lv.DeltaPhi(bjets[ind_l1bxMT2min])),weight);
			}
		}
		if((matchedgenp_tob1==-1 && matchedgenp_tob2>=0)||(matchedgenp_tob1>=0 && matchedgenp_tob2==-1)){
			int bshortind, blongind, truematchedb,truenotmatchedb, notselectedbutmatchedjet, otherb;
			if(matchedgenp_tob1>=0){//b1 is matched, probe b2
				bshortind = ind_l2bxMT2min; blongind = bind[ind_l2bxMT2min]; otherb = bind[ind_l1bxMT2min];
				truematchedb = matchedgenp_tob1;
				truenotmatchedb = gb1; notselectedbutmatchedjet = trueb1match;
				if(gb1 == matchedgenp_tob1) { truenotmatchedb = gb2; notselectedbutmatchedjet = trueb2match; }
			}
			else if(matchedgenp_tob2>=0){//b2 is matched, probe b1
				bshortind = ind_l1bxMT2min; blongind = bind[ind_l1bxMT2min]; otherb = bind[ind_l2bxMT2min];
				truematchedb = matchedgenp_tob2;
				truenotmatchedb = gb1; notselectedbutmatchedjet = trueb1match;
				if(gb1 == matchedgenp_tob2) { truenotmatchedb = gb2; notselectedbutmatchedjet = trueb2match; }
			}
			else cout << __LINE__ << "  this should not happen: matchedgenp_tob1 " << matchedgenp_tob1 << " matchedgenp_tob2 " << matchedgenp_tob2 << endl;
			//now fill
			if(truenotmatchedb>=0) histos["DeltaPt_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(fabs(fMT2tree->genlept[truenotmatchedb].lv.Pt() - fMT2tree->jet[blongind].lv.Pt()),weight);
			if(truenotmatchedb>=0) histos["DeltaPhi_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(fabs(fMT2tree->genlept[truenotmatchedb].lv.DeltaPhi(fMT2tree->jet[blongind].lv)),weight);
			if(truenotmatchedb>=0) histos["DeltaEta_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(fabs(fMT2tree->genlept[truenotmatchedb].lv.Eta() - fMT2tree->jet[blongind].lv.Eta()),weight);
			if(truenotmatchedb>=0) histos["BDiscr_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(bjetsdiscriminant[bshortind],weight);
			if(truenotmatchedb>=0) histos["BEta_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(fMT2tree->jet[blongind].lv.Eta(),weight);
			if(truenotmatchedb>=0) histos["BPt_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(fMT2tree->jet[blongind].lv.Pt(),weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb) histos["BDiscr_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].bTagProbSSVHE,weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb) histos["BEta_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].lv.Eta(),weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb) histos["BPt_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].lv.Pt(),weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb&&fMT2tree->jet[notselectedbutmatchedjet].isPFIDTight) histos["BJetID_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(3.5,weight);
			else if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb&&fMT2tree->jet[notselectedbutmatchedjet].isPFIDMedium) histos["BJetID_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(2.5,weight);
			else if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb&&fMT2tree->jet[notselectedbutmatchedjet].isPFIDLoose) histos["BJetID_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(1.5,weight);
			else if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb&&!(fMT2tree->jet[notselectedbutmatchedjet].isPFIDLoose)) histos["BJetID_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(0.5,weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb) histos["BDiscr_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].bTagProbSSVHE,weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb) histos["BEta_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].lv.Eta(),weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb) histos["BPt_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].lv.Pt(),weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb&&fMT2tree->jet[notselectedbutmatchedjet].isPFIDTight) histos["BJetID_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(3.5,weight);
			else if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb&&fMT2tree->jet[notselectedbutmatchedjet].isPFIDMedium) histos["BJetID_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(2.5,weight);
			else if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb&&fMT2tree->jet[notselectedbutmatchedjet].isPFIDLoose) histos["BJetID_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(1.5,weight);
			else if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb&&!(fMT2tree->jet[notselectedbutmatchedjet].isPFIDLoose)) histos["BJetID_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(0.5,weight);
		}
		if(trueb1match<0 && trueb2match<0) histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(2.5,2.5,weight);//no matched jet
		else if(trueb1match<0){
			if(trueb2match==bind[ind_l1bxMT2min] || trueb2match==bind[ind_l2bxMT2min])//matched jet is also selected
				histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(2.5,0.5,weight);
			else histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(2.5,1.5,weight);//matched jet is not selected
		}
		else if(trueb2match<0){
			if(trueb1match==bind[ind_l1bxMT2min] || trueb1match==bind[ind_l2bxMT2min])//matched jet is also selected
				histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(0.5,2.5,weight);
			else histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(1.5,2.5,weight);//matched jet is not selected
		}
		else{//both genb have a matched jet
			if((trueb1match==bind[ind_l1bxMT2min]&&trueb2match==bind[ind_l2bxMT2min]) || (trueb2match==bind[ind_l1bxMT2min]&&trueb1match==bind[ind_l2bxMT2min]))
				histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(0.5,0.5,weight);//both matched jets are selected
			else if(trueb1match==bind[ind_l1bxMT2min] || trueb1match==bind[ind_l2bxMT2min])
				histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(0.5,1.5,weight);//one matched jet is selected, other not
			else if(trueb2match==bind[ind_l1bxMT2min] || trueb2match==bind[ind_l2bxMT2min])
				histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(1.5,0.5,weight);//one matched jet is selected, other not
			else 
				histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(1.5,1.5,weight);//both matched jets are not selected
		}
	} //MT2lb>200
	else {
		histos["BEta_MT2lb_lt_200"+hs]->Fill(bjets[ind_l1bxMT2min].Eta(),weight);
		histos["BEta_MT2lb_lt_200"+hs]->Fill(bjets[ind_l2bxMT2min].Eta(),weight);
		histos["BPt_MT2lb_lt_200"+hs]->Fill(bjets[ind_l1bxMT2min].Pt(),weight);
		histos["BPt_MT2lb_lt_200"+hs]->Fill(bjets[ind_l2bxMT2min].Pt(),weight);
		histos["BDiscr_MT2lb_lt_200"+hs]->Fill(bjetsdiscriminant[ind_l1bxMT2min],weight);
		histos["BDiscr_MT2lb_lt_200"+hs]->Fill(bjetsdiscriminant[ind_l2bxMT2min],weight);
		if(gb1>=0){
			histos["TrueBfromtopEta_MT2lb_lt_200"+hs]->Fill(fMT2tree->genlept[gb1].lv.Eta(),weight);
			histos["TrueBfromtopPt_MT2lb_lt_200"+hs]->Fill(fMT2tree->genlept[gb1].lv.Pt(),weight);
		} if(gb2>=0){
			histos["TrueBfromtopEta_MT2lb_lt_200"+hs]->Fill(fMT2tree->genlept[gb2].lv.Eta(),weight);
			histos["TrueBfromtopPt_MT2lb_lt_200"+hs]->Fill(fMT2tree->genlept[gb2].lv.Pt(),weight);
		}
		if(fMT2tree->genmet[0].Pt()>0.001){
			histos["PFMEToverGenMET_MT2lb_lt_200"+hs]->Fill(met.Pt()/fMT2tree->genmet[0].Pt(),weight);
			histos["DeltaGenMET_PFMET_Phi_MT2lb_lt_200"+hs]->Fill(fabs(met.DeltaPhi(fMT2tree->genmet[0])),weight);
		}
		if(gl1>=0&&gl2>=0&&gb1>=0&&gb2>=0){
			if(sqrt(pow(l1.DeltaR(fMT2tree->genlept[gl1].lv),2)+pow(l2.DeltaR(fMT2tree->genlept[gl2].lv),2))<sqrt(pow(l2.DeltaR(fMT2tree->genlept[gl1].lv),2)+pow(l1.DeltaR(fMT2tree->genlept[gl2].lv),2))){
				//l1 matched to genl1, l2 matched to genl2 better than other way round
				histos["DeltaPt_TrueBparton_selectedjet_MT2lb_lt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb2].lv.Pt()-bjets[ind_l2bxMT2min].Pt()),weight);
				histos["DeltaPt_TrueBparton_selectedjet_MT2lb_lt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb1].lv.Pt()-bjets[ind_l1bxMT2min].Pt()),weight);
				histos["DeltaEta_TrueBparton_selectedjet_MT2lb_lt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb2].lv.Eta()-bjets[ind_l2bxMT2min].Eta()),weight);
				histos["DeltaEta_TrueBparton_selectedjet_MT2lb_lt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb1].lv.Eta()-bjets[ind_l1bxMT2min].Eta()),weight);
				histos["DeltaPhi_TrueBparton_selectedjet_MT2lb_lt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb2].lv.DeltaPhi(bjets[ind_l2bxMT2min])),weight);
				histos["DeltaPhi_TrueBparton_selectedjet_MT2lb_lt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb1].lv.DeltaPhi(bjets[ind_l1bxMT2min])),weight);
			} else{
				//l1 matched to genl2, l2 matched to genl1 better than other way round
				histos["DeltaPt_TrueBparton_selectedjet_MT2lb_lt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb1].lv.Pt()-bjets[ind_l2bxMT2min].Pt()),weight);
				histos["DeltaPt_TrueBparton_selectedjet_MT2lb_lt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb2].lv.Pt()-bjets[ind_l1bxMT2min].Pt()),weight);
				histos["DeltaEta_TrueBparton_selectedjet_MT2lb_lt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb1].lv.Eta()-bjets[ind_l2bxMT2min].Eta()),weight);
				histos["DeltaEta_TrueBparton_selectedjet_MT2lb_lt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb2].lv.Eta()-bjets[ind_l1bxMT2min].Eta()),weight);
				histos["DeltaPhi_TrueBparton_selectedjet_MT2lb_lt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb1].lv.DeltaPhi(bjets[ind_l2bxMT2min])),weight);
				histos["DeltaPhi_TrueBparton_selectedjet_MT2lb_lt_200"+hs]->Fill(fabs(fMT2tree->genlept[gb2].lv.DeltaPhi(bjets[ind_l1bxMT2min])),weight);
			}
		}
		if((matchedgenp_tob1==-1 && matchedgenp_tob2>=0)||(matchedgenp_tob1>=0 && matchedgenp_tob2==-1)){
			int bshortind, blongind, truematchedb,truenotmatchedb, notselectedbutmatchedjet, otherb;
			if(matchedgenp_tob1>=0){//b1 is matched, probe b2
				bshortind = ind_l2bxMT2min; blongind = bind[ind_l2bxMT2min]; otherb = bind[ind_l1bxMT2min];
				truematchedb = matchedgenp_tob1;
				truenotmatchedb = gb1; notselectedbutmatchedjet = trueb1match;
				if(gb1 == matchedgenp_tob1) { truenotmatchedb = gb2; notselectedbutmatchedjet = trueb2match; }
			}
			else if(matchedgenp_tob2>=0){//b2 is matched, probe b1
				bshortind = ind_l1bxMT2min; blongind = bind[ind_l1bxMT2min]; otherb = bind[ind_l2bxMT2min];
				truematchedb = matchedgenp_tob2;
				truenotmatchedb = gb1; notselectedbutmatchedjet = trueb1match;
				if(gb1 == matchedgenp_tob2) { truenotmatchedb = gb2; notselectedbutmatchedjet = trueb2match; }
			}
			else cout << __LINE__ << "  this should not happen: matchedgenp_tob1 " << matchedgenp_tob1 << " matchedgenp_tob2 " << matchedgenp_tob2 << endl;
			//now fill
			if(truenotmatchedb>=0) histos["DeltaPt_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(fabs(fMT2tree->genlept[truenotmatchedb].lv.Pt() - fMT2tree->jet[blongind].lv.Pt()),weight);
			if(truenotmatchedb>=0) histos["DeltaPhi_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(fabs(fMT2tree->genlept[truenotmatchedb].lv.DeltaPhi(fMT2tree->jet[blongind].lv)),weight);
			if(truenotmatchedb>=0) histos["DeltaEta_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(fabs(fMT2tree->genlept[truenotmatchedb].lv.Eta() - fMT2tree->jet[blongind].lv.Eta()),weight);
			if(truenotmatchedb>=0) histos["BDiscr_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(bjetsdiscriminant[bshortind],weight);
			if(truenotmatchedb>=0) histos["BEta_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(fMT2tree->jet[blongind].lv.Eta(),weight);
			if(truenotmatchedb>=0) histos["BPt_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(fMT2tree->jet[blongind].lv.Pt(),weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb) histos["BDiscr_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].bTagProbSSVHE,weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb) histos["BEta_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].lv.Eta(),weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb) histos["BPt_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].lv.Pt(),weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb&&fMT2tree->jet[notselectedbutmatchedjet].isPFIDTight) histos["BJetID_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(3.5,weight);
			else if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb&&fMT2tree->jet[notselectedbutmatchedjet].isPFIDMedium) histos["BJetID_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(2.5,weight);
			else if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb&&fMT2tree->jet[notselectedbutmatchedjet].isPFIDLoose) histos["BJetID_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(1.5,weight);
			else if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb&&!(fMT2tree->jet[notselectedbutmatchedjet].isPFIDLoose)) histos["BJetID_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(0.5,weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb) histos["BDiscr_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].bTagProbSSVHE,weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb) histos["BEta_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].lv.Eta(),weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb) histos["BPt_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].lv.Pt(),weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb&&fMT2tree->jet[notselectedbutmatchedjet].isPFIDTight) histos["BJetID_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(3.5,weight);
			else if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb&&fMT2tree->jet[notselectedbutmatchedjet].isPFIDMedium) histos["BJetID_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(2.5,weight);
			else if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb&&fMT2tree->jet[notselectedbutmatchedjet].isPFIDLoose) histos["BJetID_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(1.5,weight);
			else if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet==otherb&&!(fMT2tree->jet[notselectedbutmatchedjet].isPFIDLoose)) histos["BJetID_notselectedbutdoublematchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(0.5,weight);
		}
		if(trueb1match<0 && trueb2match<0) histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(2.5,2.5,weight);//no matched jet
		else if(trueb1match<0){
			if(trueb2match==bind[ind_l1bxMT2min] || trueb2match==bind[ind_l2bxMT2min])//matched jet is also selected
				histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(2.5,0.5,weight);
			else histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(2.5,1.5,weight);//matched jet is not selected
		}
		else if(trueb2match<0){
			if(trueb1match==bind[ind_l1bxMT2min] || trueb1match==bind[ind_l2bxMT2min])//matched jet is also selected
				histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(0.5,2.5,weight);
			else histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(1.5,2.5,weight);//matched jet is not selected
		}
		else{//both genb have a matched jet
			if((trueb1match==bind[ind_l1bxMT2min]&&trueb2match==bind[ind_l2bxMT2min]) || (trueb2match==bind[ind_l1bxMT2min]&&trueb1match==bind[ind_l2bxMT2min]))
				histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(0.5,0.5,weight);//both matched jets are selected
			else if(trueb1match==bind[ind_l1bxMT2min] || trueb1match==bind[ind_l2bxMT2min])
				histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(0.5,1.5,weight);//one matched jet is selected, other not
			else if(trueb2match==bind[ind_l1bxMT2min] || trueb2match==bind[ind_l2bxMT2min])
				histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(1.5,0.5,weight);//one matched jet is selected, other not
			else 
				histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(1.5,1.5,weight);//both matched jets are not selected
		}
	}//MT2lb <=200
	if(matchedgenl_tol1>=0 && matchedgenl_tol2>=0 && matchedgenl_tol1!=matchedgenl_tol2)//not match to same b
		histos["MT2lb_bothleptons_true"+hs]->Fill(MT2minV, weight);
	else if((matchedgenl_tol1==-1 && matchedgenl_tol2>=0)||(matchedgenl_tol1>=0 && matchedgenl_tol2==-1))
		histos["MT2lb_onelepton_true"+hs]->Fill(MT2minV, weight);
	else if(matchedgenl_tol1==-1 && matchedgenl_tol2==-1)
		histos["MT2lb_nolepton_true"+hs]->Fill(MT2minV, weight);

//rest
	if(l1e) histos["LepMT"+hs]->Fill(fMT2tree->ele[0].MT, weight);//only leading lepton
	else    histos["LepMT"+hs]->Fill(fMT2tree->muo[0].MT, weight);
	histos["LepPt" + hs]->Fill(l1.Pt(), weight);
	histos["LepPt" + hs]->Fill(l2.Pt(), weight);
	histos["BPt" + hs]->Fill(bjets[ind_l1bxMT2min].Pt(), weight);
	histos["BPt" + hs]->Fill(bjets[ind_l2bxMT2min].Pt(), weight);
	for(unsigned int n = 0; n<bjets.size(); ++n) histos["AllBPt" + hs]->Fill(bjets[n].Pt(), weight);
	histos["NBJets20" + hs]->Fill(nbjets20, weight);
	histos["NBJets30" + hs]->Fill(nbjets30, weight);
	histos["NBJets40" + hs]->Fill(nbjets40, weight);
	histos["NBJets50" + hs]->Fill(nbjets50, weight);
	for(unsigned int n = 0; n<jind.size(); ++n) histos["BtaggingDiscriminant" + hs]->Fill(fMT2tree->jet[jind[n] ].bTagProbSSVHE, weight);
	histos["BtaggingDiscriminant_MT2minbjets"+hs]->Fill(bjetsdiscriminant[ind_l1bxMT2min], weight);
	histos["BtaggingDiscriminant_MT2minbjets"+hs]->Fill(bjetsdiscriminant[ind_l2bxMT2min], weight);
	histos["HT" + hs]->Fill(fMT2tree->misc.HT, weight);
	histos["MET" + hs]->Fill(fMT2tree->misc.MET, weight);
	histos["MT2min" + hs]->Fill(MT2minV, weight);
	histos["MT2min_massless" + hs]->Fill(MT2minV_massless, weight);
	if(MT2minV_withMlb>=0 ) histos["MT2min_withMlbcut" + hs]->Fill(MT2minV_withMlb, weight);
	if(MT2minV_massless_withMlb>=0 ) histos["MT2min_massless_withMlbcut" + hs]->Fill(MT2minV_massless_withMlb, weight);
	if(MT2minV_withMlb>=0 && fMT2tree->pileUp.NVertices<=5) histos["MT2min_withMlbcut_PUle5" + hs]->Fill(MT2minV_withMlb, weight);
	if(MT2minV_withMlb>=0 && fMT2tree->pileUp.NVertices>=9) histos["MT2min_withMlbcut_PUge9" + hs]->Fill(MT2minV_withMlb, weight);
	histos["No_PV" + hs]->Fill(fMT2tree->pileUp.NVertices, weight);
	for(unsigned int n = 0; n<NumBJets; ++n){
		histos["Mlb_allcombi" + hs]->Fill(Ml1b[n], weight);
		histos["Mlb_allcombi" + hs]->Fill(Ml2b[n], weight);
	}
	histos["Mlb_combiforMT2min" + hs]->Fill(Ml1b[ind_l1bxMT2min], weight);
	histos["Mlb_combiforMT2min" + hs]->Fill(Ml2b[ind_l2bxMT2min], weight);
	histos["Mbb" + hs]->Fill(Mbbarr[ind_l1bxMT2min][ind_l2bxMT2min], weight);
	for(int i1 =0; i1<int(NumBJets-1.); ++i1){ for(unsigned int i2 = i1+1; i2<NumBJets; ++i2){ histos["MbbAll" + hs]->Fill(Mbbarr[i1][i2], weight); } }//new
	histos["MT2ll" + hs]->Fill(MT2ll, weight);
	if(sampletype=="data" && MT2ll>90.) cout <<endl<<endl << "Event/Lumi/Run " << fMT2tree->misc.Event<<"/"<<fMT2tree->misc.LumiSection<<"/"<<fMT2tree->misc.Run<<" MT2ll " << MT2ll << endl << endl << endl;
	if(Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180. )histos["MT2ll_afterMlbcut" + hs]->Fill(MT2ll, weight);
	if(MT2minV_withMlb>=0 && MT2ll>=90. )histos["MT2min_afterMlbcut_afterMT2llge90" + hs]->Fill(MT2minV_withMlb, weight);
	if(MT2minV>=0 && MT2ll>=90. )histos["MT2min_afterMT2llge90" + hs]->Fill(MT2minV, weight);
	if(MT2minV_massless>=0 && MT2ll>=90. )histos["MT2min_massless_afterMT2llge90" + hs]->Fill(MT2minV_massless, weight);
	if(MT2bbmin>=0.) histos["MT2bb" + hs]->Fill(MT2bbmin, weight);;
	if(MT2bbmin_mW_linMET>=0.) histos["MT2bb_lInMET_testmass80" + hs]->Fill(MT2bbmin_mW_linMET, weight);
	if(MT2bbmin_linMET>=0.) histos["MT2bb_lInMET" + hs]->Fill(MT2bbmin_linMET, weight);
	if(MT2bbmin_mW_linMET>=0. && Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180.)histos["MT2bb_lInMET_testmass80_afterMlbcut" + hs]->Fill(MT2bbmin_mW_linMET, weight);
	for(unsigned int n = 0; n<jind.size(); ++n) histos["JetsPt" + hs]->Fill(fMT2tree->jet[jind[n] ].lv.Pt(), weight);
	histos["UTM" + hs]->Fill(UTMvec[ind_l1bxMT2min][ind_l2bxMT2min], weight);
	if(D2!=9999999.) histos["D2" + hs]->Fill(D2, weight);
	if(Ml1b[ind_l1bxMT2min]<180. && Ml2b[ind_l2bxMT2min]<180. &&D2!=9999999.) histos["D2_withMlb"+hs]->Fill(D2, weight);
	for(unsigned int n=0; n<NumBJets; ++n){
		histos["DPhi_lb_allcombi" + hs]->Fill(DPhil1b[n], weight);
		histos["DPhi_lb_allcombi" + hs]->Fill(DPhil2b[n], weight);
		histos["DR_lb_allcombi" + hs]->Fill(DRl1b[n], weight);
		histos["DR_lb_allcombi" + hs]->Fill(DRl1b[n], weight);
		for(unsigned int n2 = n+1; n2<NumBJets; ++n2){
			histos["DPhi_bbAll" + hs]->Fill(fabs(bjets[n].DeltaPhi(bjets[n2])), weight);
			histos["DR_bbAll" + hs]->Fill(DRbbarr[n][n2], weight);
		}
	}
	histos["DPhi_lb_combiforMT2min" + hs]->Fill(DPhil1b[ind_l1bxMT2min], weight);
	histos["DPhi_lb_combiforMT2min" + hs]->Fill(DPhil2b[ind_l2bxMT2min], weight);
	histos["DR_lb_combiforMT2min" + hs]->Fill(DRl1b[ind_l1bxMT2min], weight);
	histos["DR_lb_combiforMT2min" + hs]->Fill(DRl2b[ind_l2bxMT2min], weight);
	histos["DPhi_ll" + hs]->Fill(DPhill, weight);
	histos["DPhi_bb" + hs]->Fill(fabs(bjets[ind_l1bxMT2min].DeltaPhi(bjets[ind_l2bxMT2min])), weight);
	histos["DR_ll" + hs]->Fill(DRll, weight);
	histos["DR_bb" + hs]->Fill(DRbbarr[ind_l1bxMT2min][ind_l2bxMT2min], weight);
	for(unsigned int i1 = 0; i1<NumBJets; ++ i1){
	for(unsigned int i2 = 0; i2<NumBJets; ++ i2){
		if(i1==i2) continue;
		if(truel1bxtop[i1]) histos["Mlb_genmatchingTop"+hs]->Fill(Ml1b[i1],weight);
		if(truel2bxtop[i2]) histos["Mlb_genmatchingTop"+hs]->Fill(Ml2b[i2],weight);
		if(truel1bxtop[i1] && truel2bxtop[i2]) histos["MT2_genmatchingTops"+hs]->Fill(MT2_l1b[i1][i2],weight);
		if(i1!=unsigned int(ind_l1bxMT2min) && truel1bxtop[i1]) truel1notb1 = true;
		if(i2!=unsigned int(ind_l2bxMT2min) && truel2bxtop[i2]) truel2notb2 = true;
	}}
	histos2["MT2min_vs_UTM"+hs]->Fill(MT2minV, UTMvec[ind_l1bxMT2min][ind_l2bxMT2min], weight);
	if(MT2minV_withMlb>=0) histos2["MT2min_withMlb_vs_UTM"+hs]->Fill(MT2minV_withMlb, UTMvec[ind_l1bxMT2minMlb][ind_l1bxMT2minMlb], weight);
	//x is l1b1,l1b2, y is l2b1,l2b2
	if(DRl1b[ind_l1bxMT2min]<DRl1b[ind_l2bxMT2min] && DRl2b[ind_l1bxMT2min]<DRl2b[ind_l2bxMT2min]) histos2["DR_lb_configuration" + hs]->Fill(0.5,0.5, weight);
	else if (DRl1b[ind_l1bxMT2min]<DRl1b[ind_l2bxMT2min] && DRl2b[ind_l1bxMT2min]>=DRl2b[ind_l2bxMT2min]) histos2["DR_lb_configuration" + hs]->Fill(0.5,1.5, weight);
	else if(DRl1b[ind_l1bxMT2min]>=DRl1b[ind_l2bxMT2min] && DRl2b[ind_l1bxMT2min]<DRl2b[ind_l2bxMT2min]) histos2["DR_lb_configuration" + hs]->Fill(1.5,0.5, weight);
	else if(DRl1b[ind_l1bxMT2min]>=DRl1b[ind_l2bxMT2min] && DRl2b[ind_l1bxMT2min]>=DRl2b[ind_l2bxMT2min]) histos2["DR_lb_configuration" + hs]->Fill(1.5,1.5, weight);
	if(abs(bjetsflavour[ind_l1bxMT2min])==5) histos2["trueb1b2_configuration"+hs]->Fill(0.5,0.5,weight);
	else                                     histos2["trueb1b2_configuration"+hs]->Fill(0.5,1.5,weight);
	if(abs(bjetsflavour[ind_l2bxMT2min])==5) histos2["trueb1b2_configuration"+hs]->Fill(1.5,0.5,weight);
	else                                     histos2["trueb1b2_configuration"+hs]->Fill(1.5,1.5,weight);
	if(truel1bxtop[ind_l1bxMT2min] && truel2bxtop[ind_l2bxMT2min]) histos2["trueTop_configuration"+hs]->Fill(0.5, 0.5, weight);//changed definition
	if(truel1bxtop[ind_l1bxMT2min] && truel2notb2) histos2["trueTop_configuration"+hs]->Fill(0.5, 1.5, weight);
	if(truel1bxtop[ind_l1bxMT2min] && !(truel2bxtop[ind_l2bxMT2min]) && !(truel2notb2)) histos2["trueTop_configuration"+hs]->Fill(0.5, 2.5, weight);
	if(truel1notb1 && truel2bxtop[ind_l2bxMT2min]) histos2["trueTop_configuration"+hs]->Fill(1.5, 0.5, weight);
	if(truel1notb1 && truel2notb2) histos2["trueTop_configuration"+hs]->Fill(1.5, 1.5, weight);
	if(truel1notb1 && !(truel2bxtop[ind_l2bxMT2min]) && !(truel2notb2)) histos2["trueTop_configuration"+hs]->Fill(1.5, 2.5, weight);
	if(!(truel1bxtop[ind_l1bxMT2min]) && !(truel1notb1) && truel2bxtop[ind_l2bxMT2min]) histos2["trueTop_configuration"+hs]->Fill(2.5, 0.5, weight);
	if(!(truel1bxtop[ind_l1bxMT2min]) && !(truel1notb1) && truel2notb2) histos2["trueTop_configuration"+hs]->Fill(2.5, 1.5, weight);
	if(!(truel1bxtop[ind_l1bxMT2min]) && !(truel1notb1) && !(truel2bxtop[ind_l2bxMT2min]) && !(truel2notb2)) histos2["trueTop_configuration"+hs]->Fill(2.5, 2.5, weight);
	if(MT2minV_withMlb>=0) histos2["MT2min_withMlbcut_vs_MT2ll"+hs]->Fill(MT2minV_withMlb, MT2ll, weight);
	if(MT2minV_withMlb>=0) histos2["MT2min_withMlbcut_vs_MT2bb_lInMET_testmass80"+hs]->Fill(MT2minV_withMlb, MT2bb_mW_linMET[ind_l1bxMT2minMlb][ind_l2bxMT2minMlb], weight);
	if(MT2minV_withMlb>=0) histos2["MT2min_withMlbcut_vs_MET"+hs]->Fill(MT2minV_withMlb, fMT2tree->misc.MET, weight);
	if(MT2minV_withMlb>=0) histos2["MT2min_withMlbcut_vs_HT"+hs]->Fill(MT2minV_withMlb, fMT2tree->misc.HT, weight);
	if(MT2minV_withMlb>=0) histos2["MT2min_withMlbcut_vs_PV"+hs]->Fill(MT2minV_withMlb, fMT2tree->pileUp.NVertices, weight);
	if(MT2minV_withMlb>=0) histos2["MT2min_withMlbcut_vs_D2"+hs]->Fill(MT2minV_withMlb ,D2, weight);
	if(MT2minV>=0) histos2["MT2min_vs_MT2ll"+hs]->Fill(MT2minV, MT2ll, weight);
	if(MT2minV>=0) histos2["MT2min_vs_MT2bb_lInMET_testmass80"+hs]->Fill(MT2minV, MT2bb_mW_linMET[ind_l1bxMT2min][ind_l2bxMT2min], weight);
	if(MT2minV>=0) histos2["MT2min_vs_MT2bb_lInMET"+hs]->Fill(MT2minV, MT2bb_linMET[ind_l1bxMT2min][ind_l2bxMT2min], weight);
	if(MT2minV>=0) histos2["MT2min_vs_MT2bb"+hs]->Fill(MT2minV, MT2bbarr[ind_l1bxMT2min][ind_l2bxMT2min], weight);//new
	if(MT2minV>=0) histos2["MT2min_vs_MET"+hs]->Fill(MT2minV, fMT2tree->misc.MET, weight);
	if(MT2minV>=0) histos2["MT2min_vs_HT"+hs]->Fill(MT2minV, fMT2tree->misc.HT, weight);
	if(MT2minV>=0) histos2["MT2min_vs_PV"+hs]->Fill(MT2minV, fMT2tree->pileUp.NVertices, weight);
	if(MT2minV>=0) histos2["MT2min_vs_D2"+hs]->Fill(MT2minV, D2, weight);
	if(MT2minV>=0) histos2["MT2min_vs_Mlbsamecombi"+hs]->Fill(MT2minV, Ml1b[ind_l1bxMT2min], weight);
	if(MT2minV>=0) histos2["MT2min_vs_Mlbsamecombi"+hs]->Fill(MT2minV, Ml2b[ind_l2bxMT2min], weight);
	if(MT2minV_massless>=0) histos2["MT2min_massless_vs_MT2ll"+hs]->Fill(MT2minV_massless, MT2ll, weight);
	if(MT2minV_massless>=0) histos2["MT2min_massless_vs_MT2bb_lInMET_testmass80"+hs]->Fill(MT2minV_massless, MT2bb_mW_linMET[ind_l1bxMT2minMassless][ind_l2bxMT2minMassless], weight);
	if(MT2minV_massless>=0) histos2["MT2min_massless_vs_MT2bb_lInMET"+hs]->Fill(MT2minV_massless, MT2bb_linMET[ind_l1bxMT2minMassless][ind_l2bxMT2minMassless], weight);
	if(MT2minV_massless>=0) histos2["MT2min_massless_vs_MT2bb"+hs]->Fill(MT2minV_massless, MT2bbarr[ind_l1bxMT2minMassless][ind_l2bxMT2minMassless], weight);//new
	if(MT2minV_massless>=0) histos2["MT2min_massless_vs_MET"+hs]->Fill(MT2minV_massless, fMT2tree->misc.MET, weight);
	if(MT2minV_massless>=0) histos2["MT2min_massless_vs_HT"+hs]->Fill(MT2minV_massless, fMT2tree->misc.HT, weight);
	if(MT2minV_massless>=0) histos2["MT2min_massless_vs_PV"+hs]->Fill(MT2minV_massless, fMT2tree->pileUp.NVertices, weight);
	if(MT2minV_massless>=0) histos2["MT2min_massless_vs_D2"+hs]->Fill(MT2minV_massless, D2, weight);
	if(MT2minV_massless>=0) histos2["MT2min_massless_vs_Mlbsamecombi"+hs]->Fill(MT2minV_massless, Ml1b[ind_l1bxMT2minMassless], weight);
	if(MT2minV_massless>=0) histos2["MT2min_massless_vs_Mlbsamecombi"+hs]->Fill(MT2minV_massless, Ml2b[ind_l2bxMT2minMassless], weight);
	if(MT2minV_massless>=0) histos2["MT2min_massless_vs_btaggingdiscriminant"+hs]->Fill(MT2minV_massless, bjetsdiscriminant[ind_l1bxMT2minMassless], weight);
	if(MT2minV_massless>=0) histos2["MT2min_massless_vs_btaggingdiscriminant"+hs]->Fill(MT2minV_massless, bjetsdiscriminant[ind_l2bxMT2minMassless], weight);
	histos2["MT2ll_vs_MT2bb_lInMET_testmass80"+hs]->Fill(MT2ll, MT2bbmin_mW_linMET, weight);
	histos2["MT2ll_vs_MT2bb_lInMET"+hs]->Fill(MT2ll, MT2bbmin_linMET, weight);
	histos2["MT2ll_vs_MT2bb"+hs]->Fill(MT2ll, MT2bbmin, weight);
	histos2["MT2ll_vs_D2"+hs]->Fill(MT2ll, D2, weight);
	histos2["MT2bb_vs_D2"+hs]->Fill(MT2bbmin, D2, weight);
	histos2["MT2bb_vs_btaggingdiscriminant"+hs]->Fill(MT2bbmin, bjetsdiscriminant[ind_b1x_MT2bbmin], weight);
	histos2["MT2bb_vs_btaggingdiscriminant"+hs]->Fill(MT2bbmin, bjetsdiscriminant[ind_b2x_MT2bbmin], weight);
	histos2["MT2min_vs_btaggingdiscriminant"+hs]->Fill(MT2minV, bjetsdiscriminant[ind_l1bxMT2min], weight);
	histos2["MT2min_vs_btaggingdiscriminant"+hs]->Fill(MT2minV, bjetsdiscriminant[ind_l2bxMT2min], weight);
	histos2["Mlb_vs_btaggingdiscriminant"+hs]->Fill(Ml1b[ind_l1bxMT2min], bjetsdiscriminant[ind_l1bxMT2min], weight);
	histos2["Mlb_vs_btaggingdiscriminant"+hs]->Fill(Ml2b[ind_l2bxMT2min], bjetsdiscriminant[ind_l2bxMT2min], weight);
	histos2["D2_vs_btaggingdiscriminant"+hs]->Fill(D2, bjetsdiscriminant[ind_l1bxMT2min], weight);
	histos2["D2_vs_btaggingdiscriminant"+hs]->Fill(D2, bjetsdiscriminant[ind_l2bxMT2min], weight);
	histos2["PUWeight_vs_PV"+hs]->Fill(fMT2tree->pileUp.Weight ,fMT2tree->pileUp.NVertices, weight);

	histos3["MT2lb_MT2ll_MT2bb"+hs]->Fill(MT2minV, MT2ll, MT2bbarr[ind_l1bxMT2min][ind_l2bxMT2min], weight);
	histos3["MT2lb_MT2ll_MT2bb_lInMET"+hs]->Fill(MT2minV, MT2ll, MT2bb_linMET[ind_l1bxMT2min][ind_l2bxMT2min], weight);
	histos3["MT2lb_MT2ll_MT2bb_lInMET_mW"+hs]->Fill(MT2minV, MT2ll, MT2bb_mW_linMET[ind_l1bxMT2min][ind_l2bxMT2min], weight);

	histos3["MT2lbMassless_MT2ll_MT2bb"+hs]->Fill(MT2minV_massless, MT2ll, MT2bbarr[ind_l1bxMT2minMassless][ind_l2bxMT2minMassless], weight);
	histos3["MT2lbMassless_MT2ll_MT2bb_lInMET"+hs]->Fill(MT2minV_massless, MT2ll, MT2bb_linMET[ind_l1bxMT2minMassless][ind_l2bxMT2minMassless], weight);
	histos3["MT2lbMassless_MT2ll_MT2bb_lInMET_mW"+hs]->Fill(MT2minV_massless, MT2ll, MT2bb_mW_linMET[ind_l1bxMT2minMassless][ind_l2bxMT2minMassless], weight);
	//Tail plots
	if(tailMT2bb && MT2bbmin>125.){
		histos["Tail_MT2bbgt125_B1Pt"+hs]->Fill(bjets[ind_b1x_MT2bbmin].Pt(),weight);
		histos["Tail_MT2bbgt125_B2Pt"+hs]->Fill(bjets[ind_b2x_MT2bbmin].Pt(),weight);
		histos["Tail_MT2bbgt125_J1Pt"+hs]->Fill(fMT2tree->jet[jind[0] ].lv.Pt(),weight);
		for(unsigned int n=0; n<jind.size(); ++n) histos["Tail_MT2bbgt125_JPt"+hs]->Fill(fMT2tree->jet[jind[n] ].lv.Pt(),weight);
		histos["Tail_MT2bbgt125_L1Pt"+hs]->Fill(l1.Pt(),weight);
		histos["Tail_MT2bbgt125_L2Pt"+hs]->Fill(l2.Pt(),weight);
		histos["Tail_MT2bbgt125_B1Phi"+hs]->Fill(bjets[ind_b1x_MT2bbmin].Phi(), weight);
		histos["Tail_MT2bbgt125_B2Phi"+hs]->Fill(bjets[ind_b2x_MT2bbmin].Phi(), weight);
		histos["Tail_MT2bbgt125_L1Phi"+hs]->Fill(l1.Phi(), weight);
		histos["Tail_MT2bbgt125_L2Phi"+hs]->Fill(l2.Phi(), weight);
		histos["Tail_MT2bbgt125_B1Eta"+hs]->Fill(bjets[ind_b1x_MT2bbmin].Eta(), weight);
		histos["Tail_MT2bbgt125_B2Eta"+hs]->Fill(bjets[ind_b2x_MT2bbmin].Eta(), weight);
		histos["Tail_MT2bbgt125_L1Eta"+hs]->Fill(l1.Eta(), weight);
		histos["Tail_MT2bbgt125_L2Eta"+hs]->Fill(l2.Eta(), weight);
		for(unsigned int n=0; n<jind.size(); ++n) histos["Tail_MT2bbgt125_JbDiscr"+hs]->Fill(fMT2tree->jet[jind[n] ].bTagProbSSVHE,weight);
		histos["Tail_MT2bbgt125_B1bDiscr"+hs]->Fill(fMT2tree->jet[bind[ind_b1x_MT2bbmin] ].bTagProbSSVHE,weight);
		histos["Tail_MT2bbgt125_B2bDiscr"+hs]->Fill(fMT2tree->jet[bind[ind_b2x_MT2bbmin] ].bTagProbSSVHE,weight);
		if((!(l1e))||(!(l2e))) histos["Tail_MT2bbgt125_MuoPtErrOverPt"+hs]->Fill(fMT2tree->muo[0].PtErr/fMT2tree->muo[0].lv.Pt(),weight);
		if((!(l1e))&&(!(l2e))) histos["Tail_MT2bbgt125_MuoPtErrOverPt"+hs]->Fill(fMT2tree->muo[1].PtErr/fMT2tree->muo[1].lv.Pt(),weight);
		histos["Tail_MT2bbgt125_DPhiBB"+hs]->Fill(DPhibbarr[ind_b1x_MT2bbmin][ind_b2x_MT2bbmin],weight);
		histos["Tail_MT2bbgt125_DPhiLL"+hs]->Fill(DPhill,weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2bbgt125_DPhiLB_all"+hs]->Fill(DPhil1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2bbgt125_DPhiLB_all"+hs]->Fill(DPhil2b[n],weight);
		histos["Tail_MT2bbgt125_DPhiLB_MT2lbmin"+hs]->Fill(DPhil1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_MT2bbgt125_DPhiLB_MT2lbmin"+hs]->Fill(DPhil2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_MT2bbgt125_DRBB"+hs]->Fill(DRbbarr[ind_b1x_MT2bbmin][ind_b2x_MT2bbmin],weight);
		histos["Tail_MT2bbgt125_DRLL"+hs]->Fill(DRll,weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2bbgt125_DRLB_all"+hs]->Fill(DRl1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2bbgt125_DRLB_all"+hs]->Fill(DRl2b[n],weight);
		histos["Tail_MT2bbgt125_DRLB_MT2lbmin"+hs]->Fill(DRl1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_MT2bbgt125_DRLB_MT2lbmin"+hs]->Fill(DRl2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_MT2bbgt125_HT"+hs]->Fill(fMT2tree->misc.HT,weight);
		histos["Tail_MT2bbgt125_MET"+hs]->Fill(met.Pt(),weight);
		histos["Tail_MT2bbgt125_PV"+hs]->Fill(fMT2tree->pileUp.NVertices,weight);
		histos["Tail_MT2bbgt125_Mll"+hs]->Fill(Mll,weight);
		histos["Tail_MT2bbgt125_Mbb"+hs]->Fill(Mbbarr[ind_b1x_MT2bbmin][ind_b2x_MT2bbmin],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2bbgt125_Mlb_all"+hs]->Fill(Ml1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2bbgt125_Mlb_all"+hs]->Fill(Ml1b[n],weight);
		histos["Tail_MT2bbgt125_Mlb_MT2lbmin"+hs]->Fill(Ml1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_MT2bbgt125_Mlb_MT2lbmin"+hs]->Fill(Ml2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_MT2bbgt125_NJets20"+hs]->Fill(fMT2tree->NJetsIDLoose,weight);
		histos["Tail_MT2bbgt125_NJets40"+hs]->Fill(fMT2tree->NJetsIDLoose40,weight);
		histos["Tail_MT2bbgt125_masslessMT2lbmin"+hs]->Fill(MT2minV_massless,weight);
		histos["Tail_MT2bbgt125_massiveMT2lbmin"+hs]->Fill(MT2minV,weight);
		histos["Tail_MT2bbgt125_MT2l"+hs]->Fill(MT2ll,weight);
		histos["Tail_MT2bbgt125_D2"+hs]->Fill(D2,weight);
	}
	if(tailMT2bb_linMET && MT2bbmin_linMET>150.){
		histos["Tail_MT2bbLinMETgt150_B1Pt"+hs]->Fill(bjets[ind_b1x_MT2bbmin_linMET].Pt(),weight);
		histos["Tail_MT2bbLinMETgt150_B2Pt"+hs]->Fill(bjets[ind_b2x_MT2bbmin_linMET].Pt(),weight);
		histos["Tail_MT2bbLinMETgt150_J1Pt"+hs]->Fill(fMT2tree->jet[jind[0] ].lv.Pt(),weight);
		for(unsigned int n=0; n<jind.size(); ++n) histos["Tail_MT2bbLinMETgt150_JPt"+hs]->Fill(fMT2tree->jet[jind[n] ].lv.Pt(),weight);
		histos["Tail_MT2bbLinMETgt150_L1Pt"+hs]->Fill(l1.Pt(),weight);
		histos["Tail_MT2bbLinMETgt150_L2Pt"+hs]->Fill(l2.Pt(),weight);
		histos["Tail_MT2bbLinMETgt150_B1Phi"+hs]->Fill(bjets[ind_b1x_MT2bbmin_linMET].Phi(), weight);
		histos["Tail_MT2bbLinMETgt150_B2Phi"+hs]->Fill(bjets[ind_b2x_MT2bbmin_linMET].Phi(), weight);
		histos["Tail_MT2bbLinMETgt150_L1Phi"+hs]->Fill(l1.Phi(), weight);
		histos["Tail_MT2bbLinMETgt150_L2Phi"+hs]->Fill(l2.Phi(), weight);
		histos["Tail_MT2bbLinMETgt150_B1Eta"+hs]->Fill(bjets[ind_b1x_MT2bbmin_linMET].Eta(), weight);
		histos["Tail_MT2bbLinMETgt150_B2Eta"+hs]->Fill(bjets[ind_b2x_MT2bbmin_linMET].Eta(), weight);
		histos["Tail_MT2bbLinMETgt150_L1Eta"+hs]->Fill(l1.Eta(), weight);
		histos["Tail_MT2bbLinMETgt150_L2Eta"+hs]->Fill(l2.Eta(), weight);
		for(unsigned int n=0; n<jind.size(); ++n) histos["Tail_MT2bbLinMETgt150_JbDiscr"+hs]->Fill(fMT2tree->jet[jind[n] ].bTagProbSSVHE,weight);
		histos["Tail_MT2bbLinMETgt150_B1bDiscr"+hs]->Fill(fMT2tree->jet[bind[ind_b1x_MT2bbmin_linMET] ].bTagProbSSVHE,weight);
		histos["Tail_MT2bbLinMETgt150_B2bDiscr"+hs]->Fill(fMT2tree->jet[bind[ind_b2x_MT2bbmin_linMET] ].bTagProbSSVHE,weight);
		if((!(l1e))||(!(l2e))) histos["Tail_MT2bbLinMETgt150_MuoPtErrOverPt"+hs]->Fill(fMT2tree->muo[0].PtErr/fMT2tree->muo[0].lv.Pt(),weight);
		if((!(l1e))&&(!(l2e))) histos["Tail_MT2bbLinMETgt150_MuoPtErrOverPt"+hs]->Fill(fMT2tree->muo[1].PtErr/fMT2tree->muo[1].lv.Pt(),weight);
		histos["Tail_MT2bbLinMETgt150_DPhiBB"+hs]->Fill(DPhibbarr[ind_b1x_MT2bbmin_linMET][ind_b2x_MT2bbmin_linMET],weight);
		histos["Tail_MT2bbLinMETgt150_DPhiLL"+hs]->Fill(DPhill,weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2bbLinMETgt150_DPhiLB_all"+hs]->Fill(DPhil1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2bbLinMETgt150_DPhiLB_all"+hs]->Fill(DPhil2b[n],weight);
		histos["Tail_MT2bbLinMETgt150_DPhiLB_MT2lbmin"+hs]->Fill(DPhil1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_MT2bbLinMETgt150_DPhiLB_MT2lbmin"+hs]->Fill(DPhil2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_MT2bbLinMETgt150_DRBB"+hs]->Fill(DRbbarr[ind_b1x_MT2bbmin_linMET][ind_b2x_MT2bbmin_linMET],weight);
		histos["Tail_MT2bbLinMETgt150_DRLL"+hs]->Fill(DRll,weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2bbLinMETgt150_DRLB_all"+hs]->Fill(DRl1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2bbLinMETgt150_DRLB_all"+hs]->Fill(DRl2b[n],weight);
		histos["Tail_MT2bbLinMETgt150_DRLB_MT2lbmin"+hs]->Fill(DRl1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_MT2bbLinMETgt150_DRLB_MT2lbmin"+hs]->Fill(DRl2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_MT2bbLinMETgt150_HT"+hs]->Fill(fMT2tree->misc.HT,weight);
		histos["Tail_MT2bbLinMETgt150_MET"+hs]->Fill(met.Pt(),weight);
		histos["Tail_MT2bbLinMETgt150_PV"+hs]->Fill(fMT2tree->pileUp.NVertices,weight);
		histos["Tail_MT2bbLinMETgt150_Mll"+hs]->Fill(Mll,weight);
		histos["Tail_MT2bbLinMETgt150_Mbb"+hs]->Fill(Mbbarr[ind_b1x_MT2bbmin_linMET][ind_b2x_MT2bbmin_linMET],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2bbLinMETgt150_Mlb_all"+hs]->Fill(Ml1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2bbLinMETgt150_Mlb_all"+hs]->Fill(Ml1b[n],weight);
		histos["Tail_MT2bbLinMETgt150_Mlb_MT2lbmin"+hs]->Fill(Ml1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_MT2bbLinMETgt150_Mlb_MT2lbmin"+hs]->Fill(Ml2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_MT2bbLinMETgt150_NJets20"+hs]->Fill(fMT2tree->NJetsIDLoose,weight);
		histos["Tail_MT2bbLinMETgt150_NJets40"+hs]->Fill(fMT2tree->NJetsIDLoose40,weight);
		histos["Tail_MT2bbLinMETgt150_masslessMT2lbmin"+hs]->Fill(MT2minV_massless,weight);
		histos["Tail_MT2bbLinMETgt150_massiveMT2lbmin"+hs]->Fill(MT2minV,weight);
		histos["Tail_MT2bbLinMETgt150_MT2l"+hs]->Fill(MT2ll,weight);
		histos["Tail_MT2bbLinMETgt150_D2"+hs]->Fill(D2,weight);
	}
	if(tailMT2ll && MT2ll>85.){
		histos["Tail_MT2llgt85_B1Pt"+hs]->Fill(bjets[ind_l1bxMT2minMassless].Pt(),weight);
		histos["Tail_MT2llgt85_B2Pt"+hs]->Fill(bjets[ind_l2bxMT2minMassless].Pt(),weight);
		histos["Tail_MT2llgt85_J1Pt"+hs]->Fill(fMT2tree->jet[jind[0] ].lv.Pt(),weight);
		for(unsigned int n=0; n<jind.size(); ++n) histos["Tail_MT2llgt85_JPt"+hs]->Fill(fMT2tree->jet[jind[n] ].lv.Pt(),weight);
		histos["Tail_MT2llgt85_L1Pt"+hs]->Fill(l1.Pt(),weight);
		histos["Tail_MT2llgt85_L2Pt"+hs]->Fill(l2.Pt(),weight);
		histos["Tail_MT2llgt85_B1Phi"+hs]->Fill(bjets[ind_l1bxMT2minMassless].Phi(), weight);
		histos["Tail_MT2llgt85_B2Phi"+hs]->Fill(bjets[ind_l2bxMT2minMassless].Phi(), weight);
		histos["Tail_MT2llgt85_L1Phi"+hs]->Fill(l1.Phi(), weight);
		histos["Tail_MT2llgt85_L2Phi"+hs]->Fill(l2.Phi(), weight);
		histos["Tail_MT2llgt85_B1Eta"+hs]->Fill(bjets[ind_l1bxMT2minMassless].Eta(), weight);
		histos["Tail_MT2llgt85_B2Eta"+hs]->Fill(bjets[ind_l2bxMT2minMassless].Eta(), weight);
		histos["Tail_MT2llgt85_L1Eta"+hs]->Fill(l1.Eta(), weight);
		histos["Tail_MT2llgt85_L2Eta"+hs]->Fill(l2.Eta(), weight);
		for(unsigned int n=0; n<jind.size(); ++n) histos["Tail_MT2llgt85_JbDiscr"+hs]->Fill(fMT2tree->jet[jind[n] ].bTagProbSSVHE,weight);
		histos["Tail_MT2llgt85_B1bDiscr"+hs]->Fill(fMT2tree->jet[bind[ind_l1bxMT2minMassless] ].bTagProbSSVHE,weight);
		histos["Tail_MT2llgt85_B2bDiscr"+hs]->Fill(fMT2tree->jet[bind[ind_l2bxMT2minMassless] ].bTagProbSSVHE,weight);
		if((!(l1e))||(!(l2e))) histos["Tail_MT2llgt85_MuoPtErrOverPt"+hs]->Fill(fMT2tree->muo[0].PtErr/fMT2tree->muo[0].lv.Pt(),weight);
		if((!(l1e))&&(!(l2e))) histos["Tail_MT2llgt85_MuoPtErrOverPt"+hs]->Fill(fMT2tree->muo[1].PtErr/fMT2tree->muo[1].lv.Pt(),weight);
		histos["Tail_MT2llgt85_DPhiBB"+hs]->Fill(DPhibbarr[ind_l1bxMT2minMassless][ind_l2bxMT2minMassless],weight);
		histos["Tail_MT2llgt85_DPhiLL"+hs]->Fill(DPhill,weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2llgt85_DPhiLB_all"+hs]->Fill(DPhil1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2llgt85_DPhiLB_all"+hs]->Fill(DPhil2b[n],weight);
		histos["Tail_MT2llgt85_DPhiLB_MT2lbmin"+hs]->Fill(DPhil1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_MT2llgt85_DPhiLB_MT2lbmin"+hs]->Fill(DPhil2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_MT2llgt85_DRBB"+hs]->Fill(DRbbarr[ind_l1bxMT2minMassless][ind_l2bxMT2minMassless],weight);
		histos["Tail_MT2llgt85_DRLL"+hs]->Fill(DRll,weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2llgt85_DRLB_all"+hs]->Fill(DRl1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2llgt85_DRLB_all"+hs]->Fill(DRl2b[n],weight);
		histos["Tail_MT2llgt85_DRLB_MT2lbmin"+hs]->Fill(DRl1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_MT2llgt85_DRLB_MT2lbmin"+hs]->Fill(DRl2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_MT2llgt85_HT"+hs]->Fill(fMT2tree->misc.HT,weight);
		histos["Tail_MT2llgt85_MET"+hs]->Fill(met.Pt(),weight);
		histos["Tail_MT2llgt85_PV"+hs]->Fill(fMT2tree->pileUp.NVertices,weight);
		histos["Tail_MT2llgt85_Mll"+hs]->Fill(Mll,weight);
		histos["Tail_MT2llgt85_Mbb"+hs]->Fill(Mbbarr[ind_l1bxMT2minMassless][ind_l2bxMT2minMassless],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2llgt85_Mlb_all"+hs]->Fill(Ml1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_MT2llgt85_Mlb_all"+hs]->Fill(Ml1b[n],weight);
		histos["Tail_MT2llgt85_Mlb_MT2lbmin"+hs]->Fill(Ml1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_MT2llgt85_Mlb_MT2lbmin"+hs]->Fill(Ml2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_MT2llgt85_NJets20"+hs]->Fill(fMT2tree->NJetsIDLoose,weight);
		histos["Tail_MT2llgt85_NJets40"+hs]->Fill(fMT2tree->NJetsIDLoose40,weight);
		histos["Tail_MT2llgt85_masslessMT2lbmin"+hs]->Fill(MT2minV_massless,weight);
		histos["Tail_MT2llgt85_massiveMT2lbmin"+hs]->Fill(MT2minV,weight);
		histos["Tail_MT2llgt85_MT2b"+hs]->Fill(MT2bbmin,weight);
		histos["Tail_MT2llgt85_MT2b_linMET"+hs]->Fill(MT2bbmin_linMET,weight);
		histos["Tail_MT2llgt85_D2"+hs]->Fill(D2,weight);
	}
	if(tailMT2lb_massive && MT2minV>210.){
		histos["Tail_massiveMT2lbgt210_B1Pt"+hs]->Fill(bjets[ind_l1bxMT2min].Pt(),weight);
		histos["Tail_massiveMT2lbgt210_B2Pt"+hs]->Fill(bjets[ind_l2bxMT2min].Pt(),weight);
		histos["Tail_massiveMT2lbgt210_J1Pt"+hs]->Fill(fMT2tree->jet[jind[0] ].lv.Pt(),weight);
		for(unsigned int n=0; n<jind.size(); ++n) histos["Tail_massiveMT2lbgt210_JPt"+hs]->Fill(fMT2tree->jet[jind[n] ].lv.Pt(),weight);
		histos["Tail_massiveMT2lbgt210_L1Pt"+hs]->Fill(l1.Pt(),weight);
		histos["Tail_massiveMT2lbgt210_L2Pt"+hs]->Fill(l2.Pt(),weight);
		histos["Tail_massiveMT2lbgt210_B1Phi"+hs]->Fill(bjets[ind_l1bxMT2min].Phi(), weight);
		histos["Tail_massiveMT2lbgt210_B2Phi"+hs]->Fill(bjets[ind_l2bxMT2min].Phi(), weight);
		histos["Tail_massiveMT2lbgt210_L1Phi"+hs]->Fill(l1.Phi(), weight);
		histos["Tail_massiveMT2lbgt210_L2Phi"+hs]->Fill(l2.Phi(), weight);
		histos["Tail_massiveMT2lbgt210_B1Eta"+hs]->Fill(bjets[ind_l1bxMT2min].Eta(), weight);
		histos["Tail_massiveMT2lbgt210_B2Eta"+hs]->Fill(bjets[ind_l2bxMT2min].Eta(), weight);
		histos["Tail_massiveMT2lbgt210_L1Eta"+hs]->Fill(l1.Eta(), weight);
		histos["Tail_massiveMT2lbgt210_L2Eta"+hs]->Fill(l2.Eta(), weight);
		for(unsigned int n=0; n<jind.size(); ++n) histos["Tail_massiveMT2lbgt210_JbDiscr"+hs]->Fill(fMT2tree->jet[jind[n] ].bTagProbSSVHE,weight);
		histos["Tail_massiveMT2lbgt210_B1bDiscr"+hs]->Fill(fMT2tree->jet[bind[ind_l1bxMT2min] ].bTagProbSSVHE,weight);
		histos["Tail_massiveMT2lbgt210_B2bDiscr"+hs]->Fill(fMT2tree->jet[bind[ind_l2bxMT2min] ].bTagProbSSVHE,weight);
		if((!(l1e))||(!(l2e))) histos["Tail_massiveMT2lbgt210_MuoPtErrOverPt"+hs]->Fill(fMT2tree->muo[0].PtErr/fMT2tree->muo[0].lv.Pt(),weight);
		if((!(l1e))&&(!(l2e))) histos["Tail_massiveMT2lbgt210_MuoPtErrOverPt"+hs]->Fill(fMT2tree->muo[1].PtErr/fMT2tree->muo[1].lv.Pt(),weight);
		histos["Tail_massiveMT2lbgt210_DPhiBB"+hs]->Fill(DPhibbarr[ind_l1bxMT2min][ind_l2bxMT2min],weight);
		histos["Tail_massiveMT2lbgt210_DPhiLL"+hs]->Fill(DPhill,weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_massiveMT2lbgt210_DPhiLB_all"+hs]->Fill(DPhil1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_massiveMT2lbgt210_DPhiLB_all"+hs]->Fill(DPhil2b[n],weight);
		histos["Tail_massiveMT2lbgt210_DPhiLB_MT2lbmin"+hs]->Fill(DPhil1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_massiveMT2lbgt210_DPhiLB_MT2lbmin"+hs]->Fill(DPhil2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_massiveMT2lbgt210_DRBB"+hs]->Fill(DRbbarr[ind_l1bxMT2min][ind_l2bxMT2min],weight);
		histos["Tail_massiveMT2lbgt210_DRLL"+hs]->Fill(DRll,weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_massiveMT2lbgt210_DRLB_all"+hs]->Fill(DRl1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_massiveMT2lbgt210_DRLB_all"+hs]->Fill(DRl2b[n],weight);
		histos["Tail_massiveMT2lbgt210_DRLB_MT2lbmin"+hs]->Fill(DRl1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_massiveMT2lbgt210_DRLB_MT2lbmin"+hs]->Fill(DRl2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_massiveMT2lbgt210_HT"+hs]->Fill(fMT2tree->misc.HT,weight);
		histos["Tail_massiveMT2lbgt210_MET"+hs]->Fill(met.Pt(),weight);
		histos["Tail_massiveMT2lbgt210_PV"+hs]->Fill(fMT2tree->pileUp.NVertices,weight);
		histos["Tail_massiveMT2lbgt210_Mll"+hs]->Fill(Mll,weight);
		histos["Tail_massiveMT2lbgt210_Mbb"+hs]->Fill(Mbbarr[ind_l1bxMT2min][ind_l2bxMT2min],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_massiveMT2lbgt210_Mlb_all"+hs]->Fill(Ml1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_massiveMT2lbgt210_Mlb_all"+hs]->Fill(Ml1b[n],weight);
		histos["Tail_massiveMT2lbgt210_Mlb_MT2lbmin"+hs]->Fill(Ml1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_massiveMT2lbgt210_Mlb_MT2lbmin"+hs]->Fill(Ml2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_massiveMT2lbgt210_NJets20"+hs]->Fill(fMT2tree->NJetsIDLoose,weight);
		histos["Tail_massiveMT2lbgt210_NJets40"+hs]->Fill(fMT2tree->NJetsIDLoose40,weight);
		histos["Tail_massiveMT2lbgt210_masslessMT2lbmin"+hs]->Fill(MT2minV_massless,weight);
		histos["Tail_massiveMT2lbgt210_MT2l"+hs]->Fill(MT2ll,weight);
		histos["Tail_massiveMT2lbgt210_MT2b"+hs]->Fill(MT2bbmin,weight);
		histos["Tail_massiveMT2lbgt210_MT2b_linMET"+hs]->Fill(MT2bbmin_linMET,weight);
		histos["Tail_massiveMT2lbgt210_D2"+hs]->Fill(D2,weight);
	}
	if(tailMT2bb_massless && MT2minV_massless>125.){
		histos["Tail_masslessMT2lbgt125_B1Pt"+hs]->Fill(bjets[ind_l1bxMT2minMassless].Pt(),weight);
		histos["Tail_masslessMT2lbgt125_B2Pt"+hs]->Fill(bjets[ind_l2bxMT2minMassless].Pt(),weight);
		histos["Tail_masslessMT2lbgt125_J1Pt"+hs]->Fill(fMT2tree->jet[jind[0] ].lv.Pt(),weight);
		for(unsigned int n=0; n<jind.size(); ++n) histos["Tail_masslessMT2lbgt125_JPt"+hs]->Fill(fMT2tree->jet[jind[n] ].lv.Pt(),weight);
		histos["Tail_masslessMT2lbgt125_L1Pt"+hs]->Fill(l1.Pt(),weight);
		histos["Tail_masslessMT2lbgt125_L2Pt"+hs]->Fill(l2.Pt(),weight);
		histos["Tail_masslessMT2lbgt125_B1Phi"+hs]->Fill(bjets[ind_l1bxMT2minMassless].Phi(), weight);
		histos["Tail_masslessMT2lbgt125_B2Phi"+hs]->Fill(bjets[ind_l2bxMT2minMassless].Phi(), weight);
		histos["Tail_masslessMT2lbgt125_L1Phi"+hs]->Fill(l1.Phi(), weight);
		histos["Tail_masslessMT2lbgt125_L2Phi"+hs]->Fill(l2.Phi(), weight);
		histos["Tail_masslessMT2lbgt125_B1Eta"+hs]->Fill(bjets[ind_l1bxMT2minMassless].Eta(), weight);
		histos["Tail_masslessMT2lbgt125_B2Eta"+hs]->Fill(bjets[ind_l2bxMT2minMassless].Eta(), weight);
		histos["Tail_masslessMT2lbgt125_L1Eta"+hs]->Fill(l1.Eta(), weight);
		histos["Tail_masslessMT2lbgt125_L2Eta"+hs]->Fill(l2.Eta(), weight);
		for(unsigned int n=0; n<jind.size(); ++n) histos["Tail_masslessMT2lbgt125_JbDiscr"+hs]->Fill(fMT2tree->jet[jind[n] ].bTagProbSSVHE,weight);
		histos["Tail_masslessMT2lbgt125_B1bDiscr"+hs]->Fill(fMT2tree->jet[bind[ind_l1bxMT2minMassless] ].bTagProbSSVHE,weight);
		histos["Tail_masslessMT2lbgt125_B2bDiscr"+hs]->Fill(fMT2tree->jet[bind[ind_l2bxMT2minMassless] ].bTagProbSSVHE,weight);
		if((!(l1e))||(!(l2e))) histos["Tail_masslessMT2lbgt125_MuoPtErrOverPt"+hs]->Fill(fMT2tree->muo[0].PtErr/fMT2tree->muo[0].lv.Pt(),weight);
		if((!(l1e))&&(!(l2e))) histos["Tail_masslessMT2lbgt125_MuoPtErrOverPt"+hs]->Fill(fMT2tree->muo[1].PtErr/fMT2tree->muo[1].lv.Pt(),weight);
		histos["Tail_masslessMT2lbgt125_DPhiBB"+hs]->Fill(DPhibbarr[ind_l1bxMT2minMassless][ind_l2bxMT2minMassless],weight);
		histos["Tail_masslessMT2lbgt125_DPhiLL"+hs]->Fill(DPhill,weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_masslessMT2lbgt125_DPhiLB_all"+hs]->Fill(DPhil1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_masslessMT2lbgt125_DPhiLB_all"+hs]->Fill(DPhil2b[n],weight);
		histos["Tail_masslessMT2lbgt125_DPhiLB_MT2lbmin"+hs]->Fill(DPhil1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_masslessMT2lbgt125_DPhiLB_MT2lbmin"+hs]->Fill(DPhil2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_masslessMT2lbgt125_DRBB"+hs]->Fill(DRbbarr[ind_l1bxMT2minMassless][ind_l2bxMT2minMassless],weight);
		histos["Tail_masslessMT2lbgt125_DRLL"+hs]->Fill(DRll,weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_masslessMT2lbgt125_DRLB_all"+hs]->Fill(DRl1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_masslessMT2lbgt125_DRLB_all"+hs]->Fill(DRl2b[n],weight);
		histos["Tail_masslessMT2lbgt125_DRLB_MT2lbmin"+hs]->Fill(DRl1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_masslessMT2lbgt125_DRLB_MT2lbmin"+hs]->Fill(DRl2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_masslessMT2lbgt125_HT"+hs]->Fill(fMT2tree->misc.HT,weight);
		histos["Tail_masslessMT2lbgt125_MET"+hs]->Fill(met.Pt(),weight);
		histos["Tail_masslessMT2lbgt125_PV"+hs]->Fill(fMT2tree->pileUp.NVertices,weight);
		histos["Tail_masslessMT2lbgt125_Mll"+hs]->Fill(Mll,weight);
		histos["Tail_masslessMT2lbgt125_Mbb"+hs]->Fill(Mbbarr[ind_l1bxMT2minMassless][ind_l2bxMT2minMassless],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_masslessMT2lbgt125_Mlb_all"+hs]->Fill(Ml1b[n],weight);
		for(unsigned int n=0; n<NumBJets; ++n) histos["Tail_masslessMT2lbgt125_Mlb_all"+hs]->Fill(Ml1b[n],weight);
		histos["Tail_masslessMT2lbgt125_Mlb_MT2lbmin"+hs]->Fill(Ml1b[ind_l1bxMT2minMassless],weight);
		histos["Tail_masslessMT2lbgt125_Mlb_MT2lbmin"+hs]->Fill(Ml2b[ind_l2bxMT2minMassless],weight);
		histos["Tail_masslessMT2lbgt125_NJets20"+hs]->Fill(fMT2tree->NJetsIDLoose,weight);
		histos["Tail_masslessMT2lbgt125_NJets40"+hs]->Fill(fMT2tree->NJetsIDLoose40,weight);
		histos["Tail_masslessMT2lbgt125_massiveMT2lbmin"+hs]->Fill(MT2minV,weight);
		histos["Tail_masslessMT2lbgt125_MT2l"+hs]->Fill(MT2ll,weight);
		histos["Tail_masslessMT2lbgt125_MT2b"+hs]->Fill(MT2bbmin,weight);
		histos["Tail_masslessMT2lbgt125_MT2b_linMET"+hs]->Fill(MT2bbmin_linMET,weight);
		histos["Tail_masslessMT2lbgt125_D2"+hs]->Fill(D2,weight);
	}

	}//while

	cout << "two osflavour b " << osb << " - two sameflavour b " << ssb << endl;
	cout << "events with non-converging discriminant " << vetoeddiscrcount << " - converging events " << discrcount << " 4thorder solution " << globsol << endl;
	cout << " Conv = " << statvec[0] << ", loop = " << statvec[1] << ", iter = " << statvec[2] << " global_4th_order_equation_solved = " << globalsolution <<  endl;
	delete fMT2tree;
	delete sttb;//discriminant solver
	delete sttbLuc;//discriminant solver


	}//for sample

	cout << "adding overflow/underflow to 1st, last bin" << endl;

	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){

		h->second->SetBinContent(1,
				    h->second->GetBinContent(0) + h->second->GetBinContent(1));
		h->second->SetBinError(1,
				  sqrt(h->second->GetBinError(0)*h->second->GetBinError(0)+
				       h->second->GetBinError(1)*h->second->GetBinError(1) ));
		h->second->SetBinContent(h->second->GetNbinsX(),
					    h->second->GetBinContent(h->second->GetNbinsX()  )+ 
					    h->second->GetBinContent(h->second->GetNbinsX()+1) );
		h->second->SetBinError(h->second->GetNbinsX(),
					  sqrt(h->second->GetBinError(h->second->GetNbinsX()  )*
					       h->second->GetBinError(h->second->GetNbinsX()  )+
					       h->second->GetBinError(h->second->GetNbinsX()+1)*
					       h->second->GetBinError(h->second->GetNbinsX()+1)  ));
	}

	cout << "Adding MuMu, EE, EMu to LL" << endl;
	for(unsigned int n =0; n<histonames.size(); ++n){
	   for(int is = 0; is<sampletypesize; ++is){
		string hs1 = string("_") + string("EE") + string("_") + sample_type[is];
		string hs2 = string("_") + string("EMu") + string("_") + sample_type[is];
		string hs3 = string("_") + string("MuMu") + string("_") + sample_type[is];
		string mapname = histonames[n];
		string h = string("_") + string("LL") + string("_") + sample_type[is];
		if(histos.count(mapname+hs1)>0) histos[mapname+h]->Add(histos[mapname+hs1],1);
		if(histos.count(mapname+hs2)>0) histos[mapname+h]->Add(histos[mapname+hs2],1);
		if(histos.count(mapname+hs3)>0) histos[mapname+h]->Add(histos[mapname+hs3],1);
		if(histos2.count(mapname+hs1)>0) histos2[mapname+h]->Add(histos2[mapname+hs1],1);
		if(histos2.count(mapname+hs2)>0) histos2[mapname+h]->Add(histos2[mapname+hs2],1);
		if(histos2.count(mapname+hs3)>0) histos2[mapname+h]->Add(histos2[mapname+hs3],1);
		if(histos3.count(mapname+hs1)>0) histos3[mapname+h]->Add(histos3[mapname+hs1],1);
		if(histos3.count(mapname+hs2)>0) histos3[mapname+h]->Add(histos3[mapname+hs2],1);
		if(histos3.count(mapname+hs3)>0) histos3[mapname+h]->Add(histos3[mapname+hs3],1);
	}
	}
	cout << "Adding all mc samples to mc histo" << endl;
	for(unsigned int n =0; n<histonames.size(); ++n){
	   for(int is = 0; is<sampletypesize; ++is){
	   if(sample_type[is]=="mc") continue;
	   if(sample_type[is]=="Stop") continue;
	   if(sample_type[is]=="data") continue;
	   for(int il = 0; il<leptontypesize; ++il){
		string hs = string("_") + lepton_type[il] + string("_") + sample_type[is];
		string mapname = histonames[n];
		string h = string("_") + lepton_type[il] + string("_") + (string)"mc";
		if(histos.count(mapname+hs)>0) histos[mapname+h]->Add(histos[mapname+hs],1);
		if(histos2.count(mapname+hs)>0) histos2[mapname+h]->Add(histos2[mapname+hs],1);
		if(histos3.count(mapname+hs)>0) histos3[mapname+h]->Add(histos3[mapname+hs],1);

	   }}
	}
	cout << "setting stack color" << endl;
	Legend1 -> SetFillColor(0);
   	Legend1 -> SetBorderSize(0);
	for(unsigned int n = 0; n<histonames.size(); ++n){
	   	for(int il = 0; il<leptontypesize; ++il){
		for(int is = 0; is<sampletypesize; ++is){
			string hs = string("_") + lepton_type[il] + string("_") + sample_type[is];
			if(sample_type[is]=="QCD"){
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs]->SetFillColor(kViolet);//changed to kViolet from 401 -> now TopColor
					histos[histonames[n]+hs]->SetLineColor(kViolet);
				} if(histos2.count(histonames[n]+hs)>0){
					histos2[histonames[n]+hs]->SetFillColor(kViolet);
					histos2[histonames[n]+hs]->SetLineColor(kViolet);
				} if(histos3.count(histonames[n]+hs)>0){
					histos3[histonames[n]+hs]->SetFillColor(kViolet);
					histos3[histonames[n]+hs]->SetLineColor(kViolet);
				}
			} else if(sample_type[is]=="WJets"){
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs]->SetFillColor(417);
					histos[histonames[n]+hs]->SetLineColor(417);
				} if(histos2.count(histonames[n]+hs)>0){
					histos2[histonames[n]+hs]->SetFillColor(417);
					histos2[histonames[n]+hs]->SetLineColor(417);
				} if(histos3.count(histonames[n]+hs)>0){
					histos3[histonames[n]+hs]->SetFillColor(417);
					histos3[histonames[n]+hs]->SetLineColor(417);
				}
			} else if(sample_type[is]=="ZJets"){
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs]->SetFillColor(419);
					histos[histonames[n]+hs]->SetLineColor(419);
				} if(histos2.count(histonames[n]+hs)>0){
					histos2[histonames[n]+hs]->SetFillColor(419);
					histos2[histonames[n]+hs]->SetLineColor(419);
				} if(histos3.count(histonames[n]+hs)>0){
					histos3[histonames[n]+hs]->SetFillColor(419);
					histos3[histonames[n]+hs]->SetLineColor(419);
				}
			} else if(sample_type[is]=="TTbar"){
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs]->SetFillColor(401);//changed to 401 from 600
					histos[histonames[n]+hs]->SetLineColor(401);
				} if(histos2.count(histonames[n]+hs)>0){
					histos2[histonames[n]+hs]->SetFillColor(401);
					histos2[histonames[n]+hs]->SetLineColor(401);
				} if(histos3.count(histonames[n]+hs)>0){
					histos3[histonames[n]+hs]->SetFillColor(401);
					histos3[histonames[n]+hs]->SetLineColor(401);
				}
			} else if(sample_type[is]=="SingleTop"){
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs]->SetFillColor(595);
					histos[histonames[n]+hs]->SetLineColor(595);
				} if(histos2.count(histonames[n]+hs)>0){
					histos2[histonames[n]+hs]->SetFillColor(595);
					histos2[histonames[n]+hs]->SetLineColor(595);
				} if(histos3.count(histonames[n]+hs)>0){
					histos3[histonames[n]+hs]->SetFillColor(595);
					histos3[histonames[n]+hs]->SetLineColor(595);
				}
			} else if(sample_type[is]=="TTbarV"){
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs]->SetFillColor(65);
					histos[histonames[n]+hs]->SetLineColor(65);
				} if(histos2.count(histonames[n]+hs)>0){
					histos2[histonames[n]+hs]->SetFillColor(65);
					histos2[histonames[n]+hs]->SetLineColor(65);
				} if(histos3.count(histonames[n]+hs)>0){
					histos3[histonames[n]+hs]->SetFillColor(65);
					histos3[histonames[n]+hs]->SetLineColor(65);
				}
			} else if(sample_type[is]=="VV/VVV"){
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs]->SetFillColor(635);
					histos[histonames[n]+hs]->SetLineColor(635);
				} if(histos2.count(histonames[n]+hs)>0){
					histos2[histonames[n]+hs]->SetFillColor(635);
					histos2[histonames[n]+hs]->SetLineColor(635);
				} if(histos3.count(histonames[n]+hs)>0){
					histos3[histonames[n]+hs]->SetFillColor(635);
					histos3[histonames[n]+hs]->SetLineColor(635);
				}
			} else if(sample_type[is]=="Other"){
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs]->SetFillColor(603);
					histos[histonames[n]+hs]->SetLineColor(603);
				} if(histos2.count(histonames[n]+hs)>0){
					histos2[histonames[n]+hs]->SetFillColor(603);
					histos2[histonames[n]+hs]->SetLineColor(603);
				} if(histos3.count(histonames[n]+hs)>0){
					histos3[histonames[n]+hs]->SetFillColor(603);
					histos3[histonames[n]+hs]->SetLineColor(603);
				}
			} else if(sample_type[is]=="mc"){
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs]->SetFillColor(603);
					histos[histonames[n]+hs]->SetLineColor(603);
				} if(histos2.count(histonames[n]+hs)>0){
					histos2[histonames[n]+hs]->SetFillColor(603);
					histos2[histonames[n]+hs]->SetLineColor(603);
				} if(histos3.count(histonames[n]+hs)>0){
					histos3[histonames[n]+hs]->SetFillColor(603);
					histos3[histonames[n]+hs]->SetLineColor(603);
				}
			} else if(sample_type[is]=="Stop"){
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs]->SetLineStyle(kDotted);
					histos[histonames[n]+hs]->SetLineColor(kBlack);
					histos[histonames[n]+hs]->SetLineWidth(4);
				} if(histos2.count(histonames[n]+hs)>0){
					histos2[histonames[n]+hs]->SetLineStyle(kDotted);
					histos2[histonames[n]+hs]->SetLineColor(kBlack);
					histos2[histonames[n]+hs]->SetLineWidth(4);
				} if(histos3.count(histonames[n]+hs)>0){
					histos3[histonames[n]+hs]->SetLineStyle(kDotted);
					histos3[histonames[n]+hs]->SetLineColor(kBlack);
					histos3[histonames[n]+hs]->SetLineWidth(4);
				} 
			} else if(sample_type[is]=="data"){
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs]->SetLineColor(kBlack);
					histos[histonames[n]+hs]->SetMarkerStyle(20);
					histos[histonames[n]+hs]->SetMarkerColor(kBlack);
				} if(histos2.count(histonames[n]+hs)>0){
					histos2[histonames[n]+hs]->SetLineColor(kBlack);
					histos2[histonames[n]+hs]->SetMarkerStyle(20);
					histos2[histonames[n]+hs]->SetMarkerColor(kBlack);
				} if(histos3.count(histonames[n]+hs)>0){
					histos3[histonames[n]+hs]->SetLineColor(kBlack);
					histos3[histonames[n]+hs]->SetMarkerStyle(20);
					histos3[histonames[n]+hs]->SetMarkerColor(kBlack);
				}
			} 
		}}
	}
	cout << "setting stacks and legend" << endl;
	bool leggy = true;
	for(unsigned int n = 0; n<histonames.size(); ++n){
	   	for(int il = 0; il<leptontypesize; ++il){
			string h = string("_") + lepton_type[il];
			if(stacks.count(histonames[n]+h)==0) stacks[(histonames[n])+h ] = new THStack((histonames[n] +h).c_str(), (histonames[n] +h).c_str());
			//stack filling
		   	for(int is = 0; is<sampletypesize; ++is){
				if(sample_type[is]=="Stop") continue;
				if(sample_type[is]=="data") continue;
				if(sample_type[is]=="mc") continue;
				string hs = string("_") + lepton_type[il] + string("_") + sample_type[is];
				if(histos.count(histonames[n]+hs)==0) continue;//use only 1d for HStack
				stacks[histonames[n]+h] ->Add(histos[histonames[n]+hs]);
			}
		}
		if(leggy){
	   	for(int is = 0; is<sampletypesize; ++is){
			string hs = string("_") + string("LL") + string("_") + sample_type[is];
			if(sample_type[is]=="mc") continue;
			if(histos.count(histonames[n]+hs)==0) continue;
			//legend filling
			if(sample_type[is]!="data") Legend1->AddEntry(histos[histonames[n]+hs], (sample_type[is]).c_str(), "f");
			else Legend1->AddEntry(histos[histonames[n]+hs], (sample_type[is]).c_str(), "l");
			leggy = false;
		}
		}
	}

	cout << "saving " << endl;
    TFile *fsavefile = new TFile(outputdir + outputname,"RECREATE");//do not delete this one
	fsavefile->cd();
	for(map<string,TH3D*>::iterator h=histos3.begin(); h!=histos3.end();++h){  h->second->Write(); }
	for(map<string,TH2D*>::iterator h=histos2.begin(); h!=histos2.end();++h){  h->second->Write(); }
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){  h->second->Write(); }
	for(map<string,THStack*>::iterator h=stacks.begin(); h!=stacks.end();++h){  h->second->Write(); }
	fsavefile->Close();
	cout << "Saved histograms in " << outputdir << outputname << endl;

	cout << "plotting histograms ..." << endl;
	MakePlots(tA, histos, leptontypesize, lepton_type, stacks, Legend1, histonames);
        MakeCorrelationPlots(histos2, leptontypesize, lepton_type, histonames, outputdir);
	Make3DPlots(histos3, leptontypesize, lepton_type, histonames, outputdir);

	delete tA;
	delete Legend1;
	delete legend;
	delete fsavefile;

}//function

void MakePlots(MassPlotter *tA, map<string, TH1D*> histos, const int leptontypesize, string *lepton_type, map<string, THStack*> stacks, TLegend *Legend1, vector<string> histonames){

	cout << "Plotting 1D histograms" << endl;

	for(unsigned int n = 0; n<histonames.size(); ++n){//plot 1d
	   	for(int il = 0; il<leptontypesize; ++il){
			if(plotonlyLL && lepton_type[il]!="LL") continue;
			string name = histonames[n];
			string h = string("_") + lepton_type[il];
			string hs1 = string("_") + lepton_type[il] + string("_") + string("mc");
			string hs2 = string("_") + lepton_type[il] + string("_") + string("data");
			string hs3 = string("_") + lepton_type[il] + string("_") + string("Stop");
			if(histos.count(name+hs1)==0) continue;
			TString ytitle = "Events";//histos[name+hs1]->GetYaxis()->GetTitle();
			TString xtitle = histos[name+hs1]->GetXaxis()->GetTitle();

			TString outname = name + hs3 + (logflag ? "_log" : "") + "_overlay";

		if(histos[name+hs2]->Integral()>0 || histos[name+hs1]->Integral()>0 || histos[name+hs3]->Integral()>0){
			if(!plotonlywithratio) tA->printHisto(stacks[name+h], histos[name+hs2],histos[name+hs1],histos[name+hs3], Legend1, outname, "hist", logflag, xtitle, ytitle, 2, 2, 1.);
			tA->plotRatioStack(stacks[name+h], histos[name+hs1], histos[name+hs2], histos[name+hs3], logflag, false, outname, Legend1, xtitle, ytitle, 2, 2, 1.);
		}

	}}

}

void MakeCorrelationPlots(map<string, TH2D*> histos2, const int leptontypesize, string *lepton_type, vector<string> histonames, TString outputdirectory){

	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	cout << "Do 2D Correlation plots" << endl;
	for(unsigned int n = 0; n<histonames.size(); ++n){//plot 1d
	   	for(int il = 0; il<leptontypesize; ++il){
			if(plotonlyLL && lepton_type[il]!="LL") continue;
			string name = histonames[n];
			string hs0 = name + string("_") + lepton_type[il] + string("_") + string("TTbar");
			string hs1 = name + string("_") + lepton_type[il] + string("_") + string("mc");
			string hs2 = name + string("_") + lepton_type[il] + string("_") + string("data");
			string hs3 = name + string("_") + lepton_type[il] + string("_") + string("Stop");
			col->cd();
			if(logflag) gPad->SetLogz(1);
			if(histos2.count(hs0)>0) {
				histos2[hs0]->SetMinimum(0.005);
				col->SetName(hs0.c_str() );
				col->SetTitle(hs0.c_str() );
				histos2[hs0]->Draw("COLZ");
				col->Update();
				if(histos2[hs0]->GetEntries()>0) Util::PrintNoEPS(col, hs0, outputdirectory, false);
				if(histos2[hs0]->GetEntries()>0) Util::PrintEPS(col, hs0, outputdirectory);
				col->Clear();
			}
			if(histos2.count(hs1)>0) {
				histos2[hs1]->SetMinimum(0.005);
				col->SetName(hs1.c_str() );
				col->SetTitle(hs1.c_str() );
				histos2[hs1]->Draw("COLZ");
				col->Update();
				if(histos2[hs1]->GetEntries()>0) Util::PrintNoEPS(col, hs1, outputdirectory, false);
				if(histos2[hs1]->GetEntries()>0) Util::PrintEPS(col, hs1, outputdirectory);
				col->Clear();
			}
			if(histos2.count(hs2)>0) {
				histos2[hs2]->SetMinimum(0.5);
				col->SetName(hs2.c_str() );
				col->SetTitle(hs2.c_str() );
				histos2[hs2]->Draw("COLZ");
				col->Update();
				if(histos2[hs2]->GetEntries()>0) Util::PrintNoEPS(col, hs2, outputdirectory, false);
				if(histos2[hs2]->GetEntries()>0) Util::PrintEPS(col, hs2, outputdirectory);
				col->Clear();
			}
			if(histos2.count(hs3)>0) {
				histos2[hs3]->SetMinimum(5.*10e-5);
				col->SetName(hs3.c_str() );
				col->SetTitle(hs3.c_str() );
				histos2[hs3]->Draw("COLZ");
				col->Update();
				if(histos2[hs3]->GetEntries()>0) Util::PrintNoEPS(col, hs3, outputdirectory, false);
				if(histos2[hs3]->GetEntries()>0) Util::PrintEPS(col, hs3, outputdirectory);
				col->Clear();
			}
		}
	}
	delete col;
}

void Make3DPlots(map<string, TH3D*> histos3, const int leptontypesize, string *lepton_type, vector<string> histonames, TString outputdirectory){

	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	cout << "Do 3D Correlation plots" << endl;
	for(unsigned int n = 0; n<histonames.size(); ++n){//plot 1d
		//cout << histonames[n] << endl;
	   	for(int il = 0; il<leptontypesize; ++il){
			if(plotonlyLL && lepton_type[il]!="LL") continue;
			string name = histonames[n];
			string hs0 = name + string("_") + lepton_type[il] + string("_") + string("TTbar");
			string hs1 = name + string("_") + lepton_type[il] + string("_") + string("mc");
			string hs2 = name + string("_") + lepton_type[il] + string("_") + string("data");
			string hs3 = name + string("_") + lepton_type[il] + string("_") + string("Stop");
			col->cd();
			if(histos3.count(hs0)>0) {
				col->SetName(hs0.c_str() );
				col->SetTitle(hs0.c_str() );
				histos3[hs0]->Draw("BOX");
				col->Update();
				if(histos3[hs0]->GetEntries()>0) Util::PrintNoEPS(col, hs0, outputdirectory, false);
				if(histos3[hs0]->GetEntries()>0) Util::PrintEPS(col, hs0, outputdirectory);
				col->Clear();
			}
			if(histos3.count(hs1)>0) {
				col->SetName(hs1.c_str() );
				col->SetTitle(hs1.c_str() );
				histos3[hs1]->Draw("BOX");
				col->Update();
				if(histos3[hs1]->GetEntries()>0) Util::PrintNoEPS(col, hs1, outputdirectory, false);
				if(histos3[hs1]->GetEntries()>0) Util::PrintEPS(col, hs1, outputdirectory);
				col->Clear();
			}
			if(histos3.count(hs2)>0) {
				col->SetName(hs2.c_str() );
				col->SetTitle(hs2.c_str() );
				histos3[hs2]->Draw("BOX");
				col->Update();
				if(histos3[hs2]->GetEntries()>0) Util::PrintNoEPS(col, hs2, outputdirectory, false);
				if(histos3[hs2]->GetEntries()>0) Util::PrintEPS(col, hs2, outputdirectory);
				col->Clear();
			}
			if(histos3.count(hs3)>0) {
				col->SetName(hs3.c_str() );
				col->SetTitle(hs3.c_str() );
				histos3[hs3]->Draw("BOX");
				col->Update();
				if(histos3[hs3]->GetEntries()>0) Util::PrintNoEPS(col, hs3, outputdirectory, false);
				if(histos3[hs3]->GetEntries()>0) Util::PrintEPS(col, hs3, outputdirectory);
				col->Clear();
			}
		}
	}
	delete col;
}

void load(const char* filename){
	
	fSamples.clear();
	char buffer[200];
	ifstream IN(filename);

	char ParName[100], StringValue[1000];
	float ParValue;

	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Sample File  " << filename << endl;
	int counter(0);
	
	while( IN.getline(buffer, 200, '\n') ){
		// ok = false;
		if (buffer[0] == '#') {
			continue; // Skip lines commented with '#'
		}
		if( !strcmp(buffer, "GENERAL") ){
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Path\t%s", StringValue);
			fPath = StringValue;	
			cout << fPath << endl;
			
			if(fVerbose >0){
				cout << " ----  " << endl;
				cout << "  Path " << fPath << endl;
			}

		}
		if( !strcmp(buffer, "SAMPLE")){

			sample s;
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Name\t%s", StringValue);
			s.name = TString(StringValue);

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "SName\t%s", StringValue);
			s.sname = TString(StringValue);

			IN.getline(buffer, 200, '\n');			//new
			sscanf(buffer, "ShapeName\t%s", StringValue);	//new
			s.shapename = TString(StringValue);		//new

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "File\t%s", StringValue);
			TString file =fPath+StringValue;
			TFile *f = TFile::Open(file);
			s.file = f;
			s.tree = (TTree*)f->Get("MassTree");
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Xsection\t%f", &ParValue);
			s.xsection = ParValue;
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Kfact\t%f", &ParValue);
			s.kfact = ParValue;
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Lumi\t%f", &ParValue);
			s.lumi = ParValue;

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Type\t%s", StringValue);
			s.type = StringValue;
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Color\t%f", &ParValue);
			s.color = ParValue;

			TH1F *h_PUWeights = (TH1F*) s.file->Get("h_PUWeights");
			TH1F *h_Events    = (TH1F*) s.file->Get("h_Events");
			if(h_PUWeights==0 || h_Events==0){
				cout << "ERROR: sample " << (s.file)->GetName() << " does not have PU and NEvents histos! " << endl;
				exit(1);
			}
			s.type!="data" ? s.PU_avg_weight = h_PUWeights->GetMean()    : s.PU_avg_weight =1;
			s.type!="data" ? s.nevents       = h_Events   ->GetEntries() : s.nevents       =1;
			delete h_PUWeights;
			delete h_Events;
			if(fVerbose > 0){
				cout << " ---- " << endl;
				cout << "  New sample added: " << s.name << endl;
				cout << "   Sample no.      " << counter << endl;
				cout << "   Short name:     " << s.sname << endl;
				cout << "   File:           " << (s.file)->GetName() << endl;
				cout << "   Events:         " << s.nevents  << endl;
				cout << "   Events in tree: " << s.tree->GetEntries() << endl; 
				cout << "   Xsection:       " << s.xsection << endl;
				cout << "   Lumi:           " << s.lumi << endl;
				cout << "   kfactor:        " << s.kfact << endl;
				cout << "   avg PU weight:  " << s.PU_avg_weight << endl;
				cout << "   type:           " << s.type << endl;
				cout << "   Color:          " << s.color << endl;
			}
			fSamples.push_back(s);
			counter++;
		}
	}
	if(fVerbose > 0) cout << "------------------------------------" << endl;
}
