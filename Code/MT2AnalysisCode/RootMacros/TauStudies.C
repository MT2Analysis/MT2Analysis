#include "TEventList.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"
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
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string>
#include <cmath>
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path for MT2tree.hh
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MassPlotter.hh"//use printHisto macro

//run via root -l -b -q TauStudies.C++

using namespace std;

void load(const char* filename = "samples/datasamples/samples_2141_dataonly.dat");
void Make1DPlotsRatio(THStack *hstack, TH1D *histdata, TH1D* histmc, vector<TH1D*> histsusy, int nsig,  TLegend *leg,  TString outname, TString outputdir, Option_t *drawopt,  bool logFlag, TString xtitle, TString ytitle, int njets, int nbjets, int nleps, float overlayScale);
void printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, vector<TH1D*> h_sig, int nsig, TLegend* legend,  TString canvname, TString outputdir, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale, double maxi);
void TauStudies();

//struct combining MT2trees with important information like cross section
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
const int fVerbose = 3; // already defined below
TString fPath;

//do some studies on the reconstructed hadronic tau vs. generator information
//used during implementation of hadronic taus
void TauStudies(){

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

  	gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "k");
  	gROOT->ProcessLine(".x SetStyle_PRD.C");
	

	bool MET  = true; //use MET trigger streams (low HT)
	bool HT   = false;//use HT trigger stream (medium+high HT)
	bool T2bb = true; //overlay T2bb signal models
	bool T2tt = false;//overlay T2tt signal models
	bool correcttau = true;//try to correct tau pT for invisible decay (tau neutrino in hadronic tau), from gen-info, use at own risk
	
	//samples.dat and output directory
	TString samples          = "samples/samples_type1MET_CHS_53X_HT_T2tt.dat";//dummy
	if(HT){
		if(T2tt) samples = "samples/samples_type1MET_CHS_53X_HT_T2tt.dat";
		if(T2bb) samples = "samples/samples_type1MET_CHS_53X_HT_T2bb.dat";
	}
	if(MET){
		if(T2tt) samples = "samples/samples_type1MET_CHS_53X_MET_T2tt.dat";
		if(T2bb) samples = "samples/samples_type1MET_CHS_53X_MET_T2bb.dat";
	}

	TString outputdir = "dummy";
	if(HT){
		if(T2tt) outputdir = "data_11fb-1_HT/T2tt/TauStudies/";
		if(T2bb) outputdir = "data_11fb-1_HT/T2bb/TauStudies/";
	}
	if(MET){
		if(T2tt) outputdir = "data_11fb-1_MET/T2tt/TauStudies/";
		if(T2bb) outputdir = "data_11fb-1_MET/T2bb/TauStudies/";
	}
	if(correcttau) outputdir = outputdir + "CorrectTauPt/";
    	Util::MakeOutputDir(outputdir);

	MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
	tA->setVerbose(fVerbose);
	tA->init(samples);
	tA->SetIsPhoton(false);
	tA->SetEventsPerGeV(true);
	tA->SetPileUpReweight(true);
	tA->SetSave(true);

	map<string, TH1D*>    histos;
	map<string, THStack*> stacks;
	TLegend *Legend1 = new TLegend(.71,.54,.91,.92);
	vector<string> histonames; histonames.clear();
	TLegend* legend = new TLegend(.71,.54,.91,.92);
	legend -> SetFillColor(0);
	legend -> SetBorderSize(0);


	//event selection
  std::ostringstream cutStream;
  std::ostringstream cutStreamBase;
  cutStream << " " 
    << "misc.Jet0Pass ==1"                      << "&&"
    << "misc.Jet1Pass ==1"                      << "&&"
    << "misc.Vectorsumpt < 70";
    cutStream << "&& misc.MinMetJetDPhi4 >0.3";
    cutStream << "&& NMuons==0 && NEles==0";
    if(MET) cutStream << "&&misc.MET>200&&misc.HT>450&&misc.HT<=750";
    if(HT ) cutStream << "&&misc.HT>750&&misc.MET>30";

  cutStreamBase << " " 
    << "misc.PassJetID ==1"                      << "&&"
    << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"  << "&&" // for rare SM samples
    << "misc.CSCTightHaloIDFlag == 0"                     << "&&"
    << "(misc.hcalLaserEventFlag==0|| misc.ProcessID==10)"<< "&&"
    << "misc.trackingFailureFlag==0"                      << "&&"
    << "misc.eeBadScFlag==0"                              << "&&"
    << "misc.EcalDeadCellTriggerPrimitiveFlag==0"         << "&&"
    << "misc.CrazyHCAL==0"
    << "&& misc.MET/misc.CaloMETRaw<=2.";

  std::ostringstream triggerStream;
  if(MET){
  triggerStream << "( ("
		<< "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
  		<< "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )  &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
  }
  if(HT){
  triggerStream << "( ("
		<< "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || "
		<< "trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
		<< "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
  }
  TString trigger = triggerStream.str().c_str();
  
  TString cuts = cutStream.str().c_str();
  TString basecuts = cutStreamBase.str().c_str();

    const int sampletypesize = 10;
    string sample_type[sampletypesize] = {"QCD", "WJets", "ZJets", "Top", "mc", "Stop1", "Stop2", "Stop3", "Stop4", "data"};
	//have some signal points
    if(T2tt){ sample_type[5] = "3xT2tt_250_50" ; sample_type[6] = "3xT2tt_400_125"; sample_type[7] = "3xT2tt_500_150"; sample_type[8] = "3xT2tt_650_0"; }
    if(T2bb){ sample_type[5] = "3xT2bb_400_200"; sample_type[6] = "3xT2bb_600_300"; sample_type[7] = "3xT2bb_750_0";   sample_type[8] = "empty";        }
    const int signalregionsize = 9;
    string signalregions[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b"};


	//get all histogams
	bool av = true;//append vector
	for(int is = 0; is<sampletypesize; ++is){
	for(int js = 0; js<signalregionsize; ++js){
		string hs = string("_") + signalregions[js] + string("_") + sample_type[is];
		string mapname;
		string title, xtitle, ytitle, ztitle;
		mapname = "MT2_RecoTausge1"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (reco-taus #ge 1)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_RecoTauseq0"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (reco-taus = 0)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenTausge1"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-taus #ge 1)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenTauseq0"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-taus = 0)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenHadTausge1"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-had-taus #ge 1)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenHadTauseq0"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-had-taus = 0)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenTausge1_acc"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (acc gen-taus #ge 1)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenTauseq0_acc"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (acc gen-taus = 0)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenHadTausge1_acc"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (acc gen-had-taus #ge 1)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenHadTauseq0_acc"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (acc gen-had-taus = 0)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());

		mapname = "MT2_RecoTau_matchGenTau"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (reco-tau match gen-tau)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_RecoTau_matchGenEle"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (reco-tau match gen-ele)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_RecoTau_matchGenMuo"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (reco-tau match gen-muo)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_RecoTau_matchGenJet"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (reco-tau match gen-jet)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_RecoTau_matchGenHadTau"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (reco-tau match gen-had.tau)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_RecoTau_matchGenNonHadTau"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (reco-tau match not gen-had-tau)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_RecoTau_NomatchGenTau"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (reco-tau match not gen-tau)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_RecoTau_NomatchGenHadTau"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (reco-tau no match gen-had.tau)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_RecoTau_Nomatch"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (reco-tau no match)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());

		mapname = "MT2_GenTau_matchRecoTau"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-tau match reco-tau)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenHadTau_matchRecoTau"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-had-tau match reco-tau)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenTau_NomatchRecoTau"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-tau no match reco-tau)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenHadTau_NomatchRecoTau"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-had-tau no match reco-tau)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenTau_matchRecoJet"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-tau match reco-jet)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenHadTau_matchRecoJet"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-had-tau match reco-jet)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenTau_matchRecoEle"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-tau match reco-ele)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenHadTau_matchRecoEle"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-had-tau match reco-ele)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenTau_matchRecoMuo"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-tau match reco-muo)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenHadTau_matchRecoMuo"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-had-tau match reco-muo)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenTau_Nomatch"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-tau no match)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_GenHadTau_Nomatch"; if(av) histonames.push_back(mapname); title = ""; mapname = mapname + hs; xtitle = "M_{T2} [GeV] (gen-had-tau no match)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		av = false;
	}}
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Sumw2();
		h->second->GetYaxis()->SetTitle("Events"); 
		h->second-> SetStats(false);
	}
	cout << "have " << histonames.size() << " different histgram types. " << endl;
	load(samples.Data());
   	for(size_t i = 0; i < fSamples.size(); ++i){
	    string sampletype = (string)fSamples[i].type;
	    string leptype = "LL";
	    if(sampletype==(string)"mc"){
		if(fSamples[i].sname=="QCD") sampletype = (string)fSamples[i].sname;
		else if(fSamples[i].sname=="Wtolnu") sampletype = (string)"WJets";
		else if(fSamples[i].sname=="DY"    ) sampletype = (string)"ZJets";
		else if(fSamples[i].sname =="Top"  ) sampletype = (string)"Top";
		else if(fSamples[i].sname=="TTbarV") sampletype = (string)"TTbarV";
		else sampletype = (string)"Other";
	    }
	    if(sampletype==(string)"susy") sampletype=fSamples[i].sname;
		if(sampletype == (string)"Other") cout << "EEEEEERRRRRRRROOOOOOORRRRRRRR" << endl;
	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "StopSearch: looping over " << fSamples[i].name << " added in " << sampletype << endl;
	    if(fVerbose>2) cout << "            sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts;
	    if(fSamples[i].type=="data") myCuts += " && "+trigger;
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
            Double_t weight = sample_weight;
            if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 

		string signalregiont;
		if(fMT2tree->NJetsIDLoose40 == 2 &&                                  fMT2tree->NBJets40CSVM == 0) signalregiont = "2j0b";
		if(fMT2tree->NJetsIDLoose40 == 2 &&                                  fMT2tree->NBJets40CSVM >= 1) signalregiont = "2j1to2b";
		if(fMT2tree->NJetsIDLoose40 >= 3 && fMT2tree->NJetsIDLoose40 <= 5 && fMT2tree->NBJets40CSVM == 0) signalregiont = "3to5j0b";
		if(fMT2tree->NJetsIDLoose40 >= 3 && fMT2tree->NJetsIDLoose40 <= 5 && fMT2tree->NBJets40CSVM == 1) signalregiont = "3to5j1b";
		if(fMT2tree->NJetsIDLoose40 >= 3 && fMT2tree->NJetsIDLoose40 <= 5 && fMT2tree->NBJets40CSVM == 2) signalregiont = "3to5j2b";
		if(                                                                  fMT2tree->NBJets40CSVM >= 3) signalregiont = "3b";
		if(fMT2tree->NJetsIDLoose40 >= 6 &&                                  fMT2tree->NBJets40CSVM == 0) signalregiont = "6j0b";
		if(fMT2tree->NJetsIDLoose40 >= 6 &&                                  fMT2tree->NBJets40CSVM == 1) signalregiont = "6j1b";
		if(fMT2tree->NJetsIDLoose40 >= 6 &&                                  fMT2tree->NBJets40CSVM == 2) signalregiont = "6j2b";


		string hh = "_" + signalregiont + "_" + sampletype;

		//get gen- and reco-tau vectors
	   vector<int> gentauind, geneleind, genelefromtauind, genmuoind, genmuofromtauind, realgentauind, looserecotauind;
	   gentauind.clear(); geneleind.clear(); genelefromtauind.clear(); genmuoind.clear(); genmuofromtauind.clear(); realgentauind.clear(); looserecotauind.clear();
	   vector<int> geneleneufromtauind, genmuoneufromtauind; geneleneufromtauind.clear(); genmuoneufromtauind.clear();
	   for(int j = 0; j<20; ++j){
		if(fMT2tree->genlept[j].lv.Pt()<0.01) continue;//veto zero genlepts
		if(fabs(fMT2tree->genlept[j].lv.Eta())>99999.) continue;//veto zero genlepts
		if(abs(fMT2tree->genlept[j].ID)==16&&abs(fMT2tree->genlept[j].MID)==15) gentauind.push_back(j);  //tau neutrino coming from tau decay
		if(abs(fMT2tree->genlept[j].ID)==11){
			geneleind.push_back(j);//genele
			if(abs(fMT2tree->genlept[j].MID)==15) genelefromtauind.push_back(j);
		}
		if(abs(fMT2tree->genlept[j].ID)==13){
			genmuoind.push_back(j);//genele
			if(abs(fMT2tree->genlept[j].MID)==15) genmuofromtauind.push_back(j);
		}
		if(abs(fMT2tree->genlept[j].ID)==12){
			if(abs(fMT2tree->genlept[j].MID)==15) geneleneufromtauind.push_back(j);
		}
		if(abs(fMT2tree->genlept[j].ID)==14){
			if(abs(fMT2tree->genlept[j].MID)==15) genmuoneufromtauind.push_back(j);
		}
		if(abs(fMT2tree->genlept[j].ID)==15) realgentauind.push_back(j);
	  }
	  for(int j = 0; j<8; ++j){
		if(fMT2tree->tau[j].lv.Pt()<20 )        continue;
		if(fabs(fMT2tree->tau[j].lv.Eta())>2.4) continue;
		if(fMT2tree->tau[j].Isolation<2)        continue;
		if(fMT2tree->tau[j].ElectronRej < 1)    continue;
		if(fMT2tree->tau[j].MuonRej < 3)        continue;
		looserecotauind.push_back(j);
          }

	int NRecoTaus   = looserecotauind.size();
	int NGenTaus    = gentauind.size();
	int NGenTauEles = genelefromtauind.size();
	int NGenTauMuos = genmuofromtauind.size();
	int NGenHadTaus = NGenTaus - NGenTauEles - NGenTauMuos;
	int NGenEles    = geneleind.size();
	int NGenMuos    = genmuoind.size();
	int NGenJets    = fMT2tree->NGenJets;
	int NRecoJets   = fMT2tree->NJetsIDLoose;
	if(NGenHadTaus<0) cout << "ERROR: GenHadTaus " << NGenHadTaus << " gentaus " << NGenTaus << " gentaueles " << NGenTauEles << " gentaumuos " << NGenTauMuos << endl;

	//initialize tau arrays
	bool  recomatchgen[NRecoTaus];
	bool  recomatchgenhad[NRecoTaus];
	bool  recomatchgenele[NRecoTaus];
	bool  recomatchgenmuo[NRecoTaus];
	bool  recomatchgenjet[NRecoTaus];
	int   recobestfit[NRecoTaus];
	float recobestfitdR[NRecoTaus];

	bool  genmatchreco[NGenTaus];
	bool  genmatchjet[NGenTaus];
	bool  genmatchmuo[NGenTaus];
	bool  genmatchele[NGenTaus];
	bool  genacc[NGenTaus];
	bool  genalmostacc[NGenTaus];
	bool  genhad[NGenTaus];
	bool  geneletau[NGenTaus];
	bool  genmuotau[NGenTaus];
	int   genbestfit[NGenTaus];
	float genbestfitdR[NGenTaus];
	for(int j =0; j<NRecoTaus; ++j){
		 recomatchgen[j]    = false;
		 recomatchgenhad[j] = false;
		 recomatchgenele[j] = false;
		 recomatchgenmuo[j] = false;
		 recomatchgenjet[j] = false;
		 recobestfit[j]   = -1;//0: GenHadTau, 1: GenJet, 2: GenEleTau, 3: GenMuoTau
		 recobestfitdR[j] = 999.;
	}
	for(int j = 0; j< NGenTaus; ++j){
		 genmatchreco[j] = false;
		 genmatchjet[j]  = false;
		 genmatchmuo[j]  = false;
		 genmatchele[j]  = false;
		 genacc[j]       = false;//d
		 genalmostacc[j] = false;//d
		 genhad[j]       = false;//d
		 geneletau[j]    = false;//d
		 genmuotau[j]    = false;//d
		 genbestfit[j]   = -1;// 0: RecoTau, 1: RecoJet, 2: RecoEle, 3: RecoMuo
		 genbestfitdR[j] = 999.;
	}

	int  nnrecomatchgen = 0;
	int  nnrecomatchgenhad = 0;
	int  nnrecomatchgenele = 0;
	int  nnrecomatchgenmuo = 0;
	int  nnrecomatchgenjet = 0;

	int  nngenmatchreco = 0;
	int  nngenmatchjet = 0;
	int  nngenmatchmuo = 0;
	int  nngenmatchele = 0;
	int  nngenacc = 0;
	int  nngenalmostacc = 0;
	int  nngenhad = 0;
	int  nngeneletau = 0;
	int  nngenmuotau = 0;
	int  nngenhadacc = 0;
	//check
	int grt = -1; int gt = -1;
	for(int j = 0; j<realgentauind.size(); ++j){
		TLorentzVector rgtau = fMT2tree->genlept[ realgentauind[j] ].lv;
		if(rgtau.Pt()<0.001) continue;
		for(int k = 0; k< NGenTaus; ++k){
			TLorentzVector gtau = fMT2tree->genlept[ gentauind[k] ].Mlv;
			if(gtau.Pt()<0.001) continue;
			if(fabs(rgtau.Pt()  - gtau.Pt()) >0.001) continue;
			if(fabs(rgtau.Eta() - gtau.Eta())>0.001) continue;
			if(fabs(rgtau.Phi() - gtau.Phi())>0.001) continue;
			if(fMT2tree->genlept[ gentauind[k] ].MID!=fMT2tree->genlept[ realgentauind[j] ].ID) continue;
			if(fMT2tree->genlept[ gentauind[k] ].GMID!=fMT2tree->genlept[ realgentauind[j] ].MID) continue;
			grt = j; gt = k;
		}
	}
	if(realgentauind.size()>0 && grt==-1){
		cout << "Have generated stable tau which does not match a gentau through neutrino decay " << endl;
	}

	//now start analysis:
	vector<TLorentzVector> gtaulv; gtaulv.clear();
	vector<TLorentzVector> ghadtauunlv; ghadtauunlv.clear();
	vector<TLorentzVector> ghadtaulv; ghadtaulv.clear();
	//get gen tau decau mode
	for(int j = 0; j< NGenTaus; ++j){
		TLorentzVector gtau = fMT2tree->genlept[ gentauind[j] ].Mlv;
		if(gtau.Pt()<0.001) continue;
		if(fabs(gtau.Eta())<2.4 && gtau.Pt()>20) { genacc[j] = true; ++nngenacc;}
		if(fabs(gtau.Eta())<2.7 && gtau.Pt()>15) { genalmostacc[j] = true; ++nngenalmostacc;}
		for(int k = 0; k < NGenTauEles; ++k){
			//check if mother has same Pt, eta, phi, ID within rounding errors
			if(fabs(fMT2tree->genlept[ genelefromtauind[k] ].Mlv.Pt()  - gtau.Pt()) >0.001) continue;
			if(fabs(fMT2tree->genlept[ genelefromtauind[k] ].Mlv.Eta() - gtau.Eta())>0.001) continue;
			if(fabs(fMT2tree->genlept[ genelefromtauind[k] ].Mlv.Phi() - gtau.Phi())>0.001) continue;
			if(fMT2tree->genlept[ genelefromtauind[k] ].MID!= fMT2tree->genlept[ gentauind[j] ].MID) continue;
			genhad[j] = false;
			geneletau[j] = true; ++nngeneletau;
			break;
		}
		for(int k = 0; k < NGenTauMuos; ++k){
			//check if mother has same Pt, eta, phi, ID within rounding errors
			if(fabs(fMT2tree->genlept[ genmuofromtauind[k] ].Mlv.Pt()  - gtau.Pt()) >0.001) continue;
			if(fabs(fMT2tree->genlept[ genmuofromtauind[k] ].Mlv.Eta() - gtau.Eta())>0.001) continue;
			if(fabs(fMT2tree->genlept[ genmuofromtauind[k] ].Mlv.Phi() - gtau.Phi())>0.001) continue;
			if(fMT2tree->genlept[ genmuofromtauind[k] ].MID!= fMT2tree->genlept[ gentauind[j] ].MID) continue;
			genhad[j] = false;
			genmuotau[j] = true; ++nngenmuotau;
			break;
		}
		if(genhad[j]==true){
			ghadtauunlv.push_back(gtau);
			ghadtaulv.push_back(gtau - fMT2tree->genlept[ gentauind[j] ].lv);
		}
		if(correcttau){
		if(geneletau[j]==true){
			for(int k = 0; k < geneleneufromtauind.size(); ++k){
				//check if mother has same Pt, eta, phi, ID within rounding errors
				if(fabs(fMT2tree->genlept[ geneleneufromtauind[k] ].Mlv.Pt()  - gtau.Pt()) >0.001) continue;
				if(fabs(fMT2tree->genlept[ geneleneufromtauind[k] ].Mlv.Eta() - gtau.Eta())>0.001) continue;
				if(fabs(fMT2tree->genlept[ geneleneufromtauind[k] ].Mlv.Phi() - gtau.Phi())>0.001) continue;
				if(fMT2tree->genlept[ geneleneufromtauind[k] ].MID!= fMT2tree->genlept[ gentauind[j] ].MID) continue;
				if(correcttau) gtau = gtau - fMT2tree->genlept[ geneleneufromtauind[k] ].lv;
				break;
			}
		}
		if(genmuotau[j]==true){
			for(int k = 0; k < genmuoneufromtauind.size(); ++k){
				//check if mother has same Pt, eta, phi, ID within rounding errors
				if(fabs(fMT2tree->genlept[ genmuoneufromtauind[k] ].Mlv.Pt()  - gtau.Pt()) >0.001) continue;
				if(fabs(fMT2tree->genlept[ genmuoneufromtauind[k] ].Mlv.Eta() - gtau.Eta())>0.001) continue;
				if(fabs(fMT2tree->genlept[ genmuoneufromtauind[k] ].Mlv.Phi() - gtau.Phi())>0.001) continue;
				if(fMT2tree->genlept[ genmuoneufromtauind[k] ].MID!= fMT2tree->genlept[ gentauind[j] ].MID) continue;
				if(correcttau) gtau = gtau - fMT2tree->genlept[ genmuoneufromtauind[k] ].lv;
				break;
			}
		}
		if(correcttau) gtau = gtau - fMT2tree->genlept[ gentauind[j] ].lv;
		}
		gtaulv.push_back(gtau);
	}
	//get matching of gentau to reotau
	for(int j = 0; j< NGenTaus; ++j){
		TLorentzVector gtau = gtaulv[j];
		if(gtau.Pt()<0.001) continue;
		if(genmuotau[j] == false && geneletau[j] == false){ genhad[j] = true; ++nngenhad; if(fabs(gtau.Eta())<2.4 && gtau.Pt()>20) {++nngenhadacc;} }
		for(int k = 0; k < fMT2tree->NEles; ++k){
			if(fMT2tree->ele[k].IDVeto==false) continue;
			if(fMT2tree->ele[k].lv.Pt()<10) continue;
			if(fabs(fMT2tree->ele[k].lv.Eta())>2.4) continue;
			float dR = fMT2tree->ele[k].lv.DeltaR(gtau);
			if(dR<0.3){
				genmatchele[j] = true; ++nngenmatchele;
				if(dR<=genbestfitdR[j]){
					genbestfitdR[j] = dR; genbestfit[j] = 2;
		}}}
		for(int k = 0; k < fMT2tree->NMuons; ++k){
			if(fMT2tree->muo[k].lv.Pt()<10) continue;
			if(fabs(fMT2tree->muo[k].lv.Eta())>2.4) continue;
			float dR = fMT2tree->muo[k].lv.DeltaR(gtau);
			if(dR<0.3){
				genmatchmuo[j] = true; ++nngenmatchmuo;
				if(dR<=genbestfitdR[j]){
					genbestfitdR[j] = dR; genbestfit[j] = 3;
		}}}
		for(int k = 0; k < fMT2tree->NJets; ++k){
			if(fMT2tree->jet[k].isPFIDLoose==false) continue;
			if(fMT2tree->jet[k].lv.Pt()<20) continue;
			if(fabs(fMT2tree->jet[k].lv.Eta())>2.4) continue;
			float dR = fMT2tree->jet[k].lv.DeltaR(gtau);
			if(dR<0.5){
				genmatchjet[j] = true; ++nngenmatchjet;
				if(dR<=genbestfitdR[j]){
					genbestfitdR[j] = dR; genbestfit[j] = 1;
		}}}
		for(int k = 0; k<NRecoTaus; ++k){
			if(fMT2tree->tau[ looserecotauind[k] ].lv.Pt()<20) continue;
			float dR = fMT2tree->tau[ looserecotauind[k] ].lv.DeltaR(gtau);
			if(dR<0.3){
				genmatchreco[j] = true; ++nngenmatchreco;
				if(dR<=genbestfitdR[j]){
					genbestfitdR[j] = dR; genbestfit[j] = 0;
		}}}
	}
	for(int j =0; j<NRecoTaus; ++j){
		TLorentzVector rtau = fMT2tree->tau[ looserecotauind[j] ].lv;
		if(rtau.Pt()<0.001) continue;
		for(int k = 0; k< NGenEles; ++k){
			float dR = fMT2tree->genlept[ geneleind[k] ].lv.DeltaR(rtau);
			if(dR<0.3){
				recomatchgenele[j] = true; ++nnrecomatchgenele;
				if(dR<=recobestfitdR[j]){
					recobestfitdR[j] = dR; recobestfit[j] = 22;
		}}}
		for(int k = 0; k< NGenMuos; ++k){
			float dR = fMT2tree->genlept[ genmuoind[k] ].lv.DeltaR(rtau);
			if(dR<0.3){
				recomatchgenmuo[j] = true; ++nnrecomatchgenmuo;
				if(dR<=recobestfitdR[j]){
					recobestfitdR[j] = dR; recobestfit[j] = 33;
		}}}
		for(int k = 0; k<fMT2tree->NGenJets; ++k){
			if(fMT2tree->genjet[k].lv.Pt()<0.001) continue;//save against 0
			float dR = fMT2tree->genjet[ k ].lv.DeltaR(rtau);
			if(dR<0.3){
				recomatchgenjet[j] = true; ++nnrecomatchgenjet;
				if(dR<=recobestfitdR[j]){
					recobestfitdR[j] = dR; recobestfit[j] = 1;
		}}}
		for(int k = 0; k< NGenTaus; ++k){
			TLorentzVector gtau = gtaulv[k];
			if(gtau.Pt()<0.001) continue;
			float dR = gtau.DeltaR(rtau);
			if(dR<0.3){
				recomatchgen[j] = true; ++nnrecomatchgen;
				if(genhad[k]) { recomatchgenhad[j] = true; ++nnrecomatchgenhad;}
				if(dR<=recobestfitdR[j]){
					recobestfitdR[j] = dR; if(genhad[k]) recobestfit[j] = 0; else recobestfit[j] = 4;
		}}}
	}
	//checks
	if(nngeneletau!=NGenTauEles) {cout << " nngeneletau " << nngeneletau << " NGenTauEles " << NGenTauEles << endl; cout << " NGenTaus " << NGenTaus << endl;}
	if(nngenmuotau!=NGenTauMuos) {cout << " nngenmuotau " << nngenmuotau << " NGenTauMuos " << NGenTauMuos << endl; cout << " NGenTaus " << NGenTaus << endl;}
	if(nngenhad   !=NGenHadTaus) {cout << " nngenhadtau " << nngenhad    << " NGenHadTaus " << NGenHadTaus << endl; cout << " NGenTaus " << NGenTaus << endl;}

	if(NRecoTaus>0) histos["MT2_RecoTausge1" + hh]->Fill(fMT2tree->misc.MT2, weight);
	else            histos["MT2_RecoTauseq0" + hh]->Fill(fMT2tree->misc.MT2, weight);
	if(NGenTaus>0)  histos["MT2_GenTausge1" + hh]->Fill(fMT2tree->misc.MT2, weight);
	else            histos["MT2_GenTauseq0" + hh]->Fill(fMT2tree->misc.MT2, weight);
	if(nngenhad>0)  histos["MT2_GenHadTausge1" + hh]->Fill(fMT2tree->misc.MT2, weight);
	else            histos["MT2_GenHadTauseq0" + hh]->Fill(fMT2tree->misc.MT2, weight);
	if(nngenacc>0)  histos["MT2_GenTausge1_acc" + hh]->Fill(fMT2tree->misc.MT2, weight);
	else            histos["MT2_GenTauseq0_acc" + hh]->Fill(fMT2tree->misc.MT2, weight);
	if(nngenhadacc>0)  histos["MT2_GenHadTausge1_acc" + hh]->Fill(fMT2tree->misc.MT2, weight);
	else            histos["MT2_GenHadTauseq0_acc" + hh]->Fill(fMT2tree->misc.MT2, weight);

	//order reco
	//reco match to genhadtau, gennonhadtau, genjet, ele, muo, none; gentaus separate
	if(NRecoTaus>0){
		if(nnrecomatchgen>0) histos["MT2_RecoTau_matchGenTau" + hh]->Fill(fMT2tree->misc.MT2, weight);
		else histos["MT2_RecoTau_NomatchGenTau" + hh]->Fill(fMT2tree->misc.MT2, weight);
		if(nnrecomatchgenhad>0)  histos["MT2_RecoTau_matchGenHadTau" + hh]->Fill(fMT2tree->misc.MT2, weight);
		else{
			histos["MT2_RecoTau_NomatchGenHadTau" + hh]->Fill(fMT2tree->misc.MT2, weight);
			if(nnrecomatchgen>0) histos["MT2_RecoTau_matchGenNonHadTau" + hh]->Fill(fMT2tree->misc.MT2, weight);
			else if(nnrecomatchgenjet>0) histos["MT2_RecoTau_matchGenJet" + hh]->Fill(fMT2tree->misc.MT2, weight);
			else if(nnrecomatchgenele>0) histos["MT2_RecoTau_matchGenEle" + hh]->Fill(fMT2tree->misc.MT2, weight);
			else if(nnrecomatchgenmuo>0) histos["MT2_RecoTau_matchGenMuo" + hh]->Fill(fMT2tree->misc.MT2, weight);
			else histos["MT2_RecoTau_Nomatch" + hh]->Fill(fMT2tree->misc.MT2, weight);
		}
	}
	//order gen
	// gen matched to tau,jet,ele,muo,none
	if(NGenTaus>0){
		if(nnrecomatchgenhad>0)  histos["MT2_GenTau_matchRecoTau" + hh]->Fill(fMT2tree->misc.MT2, weight);
		else{
			histos["MT2_GenTau_NomatchRecoTau" + hh]->Fill(fMT2tree->misc.MT2, weight);
			if(nnrecomatchgenjet>0) histos["MT2_GenTau_matchRecoJet" + hh]->Fill(fMT2tree->misc.MT2, weight);
			else if(nnrecomatchgenele>0) histos["MT2_GenTau_matchRecoEle" + hh]->Fill(fMT2tree->misc.MT2, weight);
			else if(nnrecomatchgenmuo>0) histos["MT2_GenTau_matchRecoMuo" + hh]->Fill(fMT2tree->misc.MT2, weight);
			else histos["MT2_GenTau_Nomatch" + hh]->Fill(fMT2tree->misc.MT2, weight);
		}
	}
	//order genhad
	// gen had match to tau, jet, ele,muo,none
	if(nngenhad>0){
		if(nnrecomatchgenhad>0)  histos["MT2_GenHadTau_matchRecoTau" + hh]->Fill(fMT2tree->misc.MT2, weight);
		else{
			histos["MT2_GenHadTau_NomatchRecoTau" + hh]->Fill(fMT2tree->misc.MT2, weight);
			if(nnrecomatchgenjet>0) histos["MT2_GenHadTau_matchRecoJet" + hh]->Fill(fMT2tree->misc.MT2, weight);
			else if(nnrecomatchgenele>0) histos["MT2_GenHadTau_matchRecoEle" + hh]->Fill(fMT2tree->misc.MT2, weight);
			else if(nnrecomatchgenmuo>0) histos["MT2_GenHadTau_matchRecoMuo" + hh]->Fill(fMT2tree->misc.MT2, weight);
			else histos["MT2_GenHadTau_Nomatch" + hh]->Fill(fMT2tree->misc.MT2, weight);
		}
	}

	}} //while and for

	cout << "adding overflow/underflow to 1st, last bin" << endl;

	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){

		h->second->SetBinContent(h->second->GetNbinsX(),
					    h->second->GetBinContent(h->second->GetNbinsX()  )+ 
					    h->second->GetBinContent(h->second->GetNbinsX()+1) );
		h->second->SetBinError(h->second->GetNbinsX(),
					  sqrt(h->second->GetBinError(h->second->GetNbinsX()  )*
					       h->second->GetBinError(h->second->GetNbinsX()  )+
					       h->second->GetBinError(h->second->GetNbinsX()+1)*
					       h->second->GetBinError(h->second->GetNbinsX()+1)  ));
		//ADD FUNCTIONS TO HTSTack AND TLEGEND
	}

	cout << "Adding all mc samples to mc histo" << endl;
	for(unsigned int n =0; n<histonames.size(); ++n){
	   for(int il = 0; il<signalregionsize; ++il){
		string hs1 = string("_") + signalregions[il] + string("_") + "QCD";
		string hs2 = string("_") + signalregions[il] + string("_") + "WJets";
		string hs3 = string("_") + signalregions[il] + string("_") + "ZJets";
		string hs4 = string("_") + signalregions[il] + string("_") + "Top";
		string mapname = histonames[n];
		string h  = string("_") + signalregions[il] + string("_") + (string)"mc";
		if(histos.count(mapname+hs1)>0)  histos[mapname+h] ->Add(histos[mapname+hs1], 1);
		if(histos.count(mapname+hs2)>0)  histos[mapname+h] ->Add(histos[mapname+hs2], 1);
		if(histos.count(mapname+hs3)>0)  histos[mapname+h] ->Add(histos[mapname+hs3], 1);
		if(histos.count(mapname+hs4)>0)  histos[mapname+h] ->Add(histos[mapname+hs4], 1);
	   }
	}

	//get tau efficiencies
	TH1D *passgentaus[signalregionsize];
	TH1D  *totgentaus[signalregionsize];
	TH1D *passgentaus_acc[signalregionsize];
	TH1D  *totgentaus_acc[signalregionsize];
	TH1D *passgenhadtaus[signalregionsize];
	TH1D  *totgenhadtaus[signalregionsize];
	TH1D *passgenhadtaus_acc[signalregionsize];
	TH1D  *totgenhadtaus_acc[signalregionsize];

	TH1D *matchedpassgentaus[signalregionsize];
	TH1D  *matchedtotgentaus[signalregionsize];
	TH1D *matchedpassgentaus_acc[signalregionsize];
	TH1D  *matchedtotgentaus_acc[signalregionsize];
	TH1D *matchedpassgenhadtaus[signalregionsize];
	TH1D  *matchedtotgenhadtaus[signalregionsize];
	TH1D *matchedpassgenhadtaus_acc[signalregionsize];
	TH1D  *matchedtotgenhadtaus_acc[signalregionsize];

	TEfficiency *taueff[signalregionsize];
	TEfficiency *hadtaueff[signalregionsize];
	TEfficiency *tauwithinacceff[signalregionsize];
	TEfficiency *hadtauwithinacceff[signalregionsize];
	TEfficiency *matchedtaueff[signalregionsize];
	TEfficiency *matchedhadtaueff[signalregionsize];
	TEfficiency *matchedtauwithinacceff[signalregionsize];
	TEfficiency *matchedhadtauwithinacceff[signalregionsize];

	//MT2 distributions
	for(int il = 0; il<signalregionsize; ++il){

		string pgt    = "MT2_RecoTausge1_"           + signalregions[il] + "_mc"; passgentaus[il] = (TH1D*)histos[pgt]->Clone((pgt+(string)"_copy").c_str());
		string tgt    = "MT2_GenTausge1"             + signalregions[il] + "_mc"; totgentaus[il] = (TH1D*)histos[tgt]->Clone((tgt+(string)"_copy").c_str());
		string pgta   = "MT2_RecoTausge1_"           + signalregions[il] + "_mc"; passgentaus_acc[il] = (TH1D*)histos[pgta]->Clone((pgta+(string)"_copy").c_str());
		string tgta   = "MT2_GenTausge1_acc"         + signalregions[il] + "_mc"; totgentaus_acc[il] = (TH1D*)histos[tgta]->Clone((tgta+(string)"_copy").c_str());
		string pght   = "MT2_RecoTausge1_"           + signalregions[il] + "_mc"; passgenhadtaus[il] = (TH1D*)histos[pght]->Clone((pght+(string)"_copy").c_str());
		string tght   = "MT2_GenHadTausge1"          + signalregions[il] + "_mc"; totgenhadtaus[il] = (TH1D*)histos[tght]->Clone((tght+(string)"_copy").c_str());
		string pghta  = "MT2_RecoTausge1_"           + signalregions[il] + "_mc"; passgenhadtaus_acc[il] = (TH1D*)histos[pghta]->Clone((pghta+(string)"_copy").c_str());
		string tghta  = "MT2_GenHadTausge1_acc"      + signalregions[il] + "_mc"; totgenhadtaus_acc[il] = (TH1D*)histos[tghta]->Clone((tghta+(string)"_copy").c_str());

		string mpgt   = "MT2_RecoTau_matchGenTau"    + signalregions[il] + "_mc"; matchedpassgentaus[il] = (TH1D*)histos[mpgt]->Clone((mpgt+(string)"_copy").c_str());
		string mtgt   = "MT2_GenTausge1"             + signalregions[il] + "_mc"; matchedtotgentaus[il] = (TH1D*)histos[mtgt]->Clone((mtgt+(string)"_copy").c_str());
		string mpgta  = "MT2_RecoTau_matchGenTau"    + signalregions[il] + "_mc"; matchedpassgentaus_acc[il] = (TH1D*)histos[mpgta]->Clone((mpgta+(string)"_copy").c_str());
		string mtgta  = "MT2_GenTausge1_acc"         + signalregions[il] + "_mc"; matchedtotgentaus_acc[il] = (TH1D*)histos[mtgta]->Clone((mtgta+(string)"_copy").c_str());
		string mpght  = "MT2_RecoTau_matchGenHadTau" + signalregions[il] + "_mc"; matchedpassgenhadtaus[il] = (TH1D*)histos[mpght]->Clone((mpght+(string)"_copy").c_str());
		string mtght  = "MT2_GenHadTausge1"          + signalregions[il] + "_mc"; matchedtotgenhadtaus[il] = (TH1D*)histos[mtght]->Clone((mtght+(string)"_copy").c_str());
		string mpghta = "MT2_RecoTau_matchGenHadTau" + signalregions[il] + "_mc"; matchedpassgenhadtaus_acc[il] = (TH1D*)histos[mpghta]->Clone((mpghta+(string)"_copy").c_str());
		string mtghta = "MT2_GenHadTausge1_acc"      + signalregions[il] + "_mc"; matchedtotgenhadtaus_acc[il] = (TH1D*)histos[mtghta]->Clone((mtghta+(string)"_copy").c_str());

		taueff[il] = new TEfficiency((*(passgentaus[il])), (*(totgentaus[il])));
		taueff[il]->SetNameTitle(("taueff"+signalregions[il]).c_str(), ("taueff"+signalregions[il]).c_str());
		tauwithinacceff[il] = new TEfficiency((*(passgentaus_acc[il])), (*(totgentaus_acc[il])));
		tauwithinacceff[il]->SetNameTitle(("tauwithinacceff"+signalregions[il]).c_str(), ("tauwithinacceff"+signalregions[il]).c_str());
		hadtaueff[il] = new TEfficiency((*(passgenhadtaus[il])), (*(totgenhadtaus[il])));
		hadtaueff[il]->SetNameTitle(("hadtaueff"+signalregions[il]).c_str(), ("hadtaueff"+signalregions[il]).c_str());
		hadtauwithinacceff[il] = new TEfficiency((*(passgenhadtaus_acc[il])), (*(totgenhadtaus_acc[il])));
		hadtauwithinacceff[il]->SetNameTitle(("hadtauwithinacceff"+signalregions[il]).c_str(), ("hadtauwithinacceff"+signalregions[il]).c_str());
		matchedtaueff[il] = new TEfficiency((*(matchedpassgentaus[il])), (*(matchedtotgentaus[il])));
		matchedtaueff[il]->SetNameTitle(("matchedtaueff"+signalregions[il]).c_str(), ("matchedtaueff"+signalregions[il]).c_str());
		matchedtauwithinacceff[il] = new TEfficiency((*(matchedpassgentaus_acc[il])), (*(matchedtotgentaus_acc[il])));
		matchedhadtaueff[il]->SetNameTitle(("matchedhadtaueff"+signalregions[il]).c_str(), ("matchedhadtaueff"+signalregions[il]).c_str());
		matchedtauwithinacceff[il] = new TEfficiency((*(matchedpassgenhadtaus[il])), (*(matchedtotgenhadtaus[il])));
		matchedtauwithinacceff[il]->SetNameTitle(("matchedtauwithinacceff"+signalregions[il]).c_str(), ("matchedtauwithinacceff"+signalregions[il]).c_str());
		matchedhadtauwithinacceff[il] = new TEfficiency((*(matchedpassgenhadtaus_acc[il])), (*(matchedtotgenhadtaus_acc[il])));
		matchedhadtauwithinacceff[il]->SetNameTitle(("matchedhadtauwithinacceff"+signalregions[il]).c_str(), ("matchedhadtauwithinacceff"+signalregions[il]).c_str());
	}


	cout << "setting stack color" << endl;
	Legend1 -> SetFillColor(0);
   	Legend1 -> SetBorderSize(0);
	for(unsigned int n = 0; n<histonames.size(); ++n){
	   	for(int il = 0; il<signalregionsize; ++il){
		for(int is = 0; is<sampletypesize; ++is){
			string hs = string("_") + signalregions[il] + string("_") + sample_type[is];
			Color_t colour = 603;//set default others if we have some failure
			     if(is==0)     colour = 401;
			else if(is==1)     colour = 417;
			else if(is==2)     colour = 419;
			else if(is==3)     colour = 855;
			else if(is==4)     colour = 0;
			else if(is==5)     colour = kRed+3;
			else if(is==6)     colour = kOrange-3;
			else if(is==7)     colour = kViolet-6;
			else if(is==8)     colour = kBlack;
			else if(is==9)     colour = 1;
			if(is<4){//mc
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs] ->SetFillColor(colour);
					histos[histonames[n]+hs] ->SetLineColor(colour);
				}
			} else if(is>4&&is<8){//signal
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs]->SetLineStyle(kDotted);
					histos[histonames[n]+hs]->SetLineColor(colour);
					histos[histonames[n]+hs]->SetLineWidth(4);
				}
			} else if(is==9){//data
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs]->SetLineColor(kBlack);
					histos[histonames[n]+hs]->SetMarkerStyle(20);
					histos[histonames[n]+hs]->SetMarkerColor(kBlack);
				}
			} 
		}}
	}
	cout << "setting stacks and legend" << endl;
	bool leggy = true;
	for(unsigned int n = 0; n<histonames.size(); ++n){
		   	for(int js = 0; js<signalregionsize; ++js){
				string hss = histonames[n] + string("_") + signalregions[js];
				if(stacks.count(hss)==0) stacks[hss ] = new THStack((hss).c_str(), "");
				stacks[hss ] ->Add(histos[(hss + string("_") + "QCD")]);
				stacks[hss ] ->Add(histos[(hss + string("_") + "WJets")]);
				stacks[hss ] ->Add(histos[(hss + string("_") + "ZJets")]);
				stacks[hss ] ->Add(histos[(hss + string("_") + "Top")]);
			}
		if(leggy){
	   	for(int is = 0; is<sampletypesize; ++is){
			string hs = string("_") + string("3b") + string("_") + sample_type[is];
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
    	TFile *fsavefile = new TFile(outputdir + "Histograms.root","RECREATE");//do not delete this one
	fsavefile->cd();
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Write();
	}
	for(map<string,THStack*>::iterator h=stacks.begin(); h!=stacks.end();++h){
		h->second->Write();
	}
	for(int il = 0; il<signalregionsize; ++il){
		taueff[il]->Write();
		hadtaueff[il]->Write();
		tauwithinacceff[il]->Write();
		hadtauwithinacceff[il]->Write();
		matchedtaueff[il]->Write();
		matchedhadtaueff[il]->Write();
		matchedtauwithinacceff[il]->Write();
		matchedhadtauwithinacceff[il]->Write();
	}
	fsavefile->Close();
	cout << "Saved histograms in " << fsavefile->GetName() << endl;

	//now print them
	for(int n = 0; n<histonames.size(); ++n){
	for(int js = 0; js<signalregionsize; ++js){
		string hsmc    = string("_") + signalregions[js] + string("_") + "mc";
		string hsdata  = string("_") + signalregions[js] + string("_") + "data";
		string hssusy1 = string("_") + signalregions[js] + string("_") + sample_type[5];
		string hssusy2 = string("_") + signalregions[js] + string("_") + sample_type[6];
		string hssusy3 = string("_") + signalregions[js] + string("_") + sample_type[7];
		string hssusy4 = string("_") + signalregions[js] + string("_") + sample_type[8];
		string hsstack = string("_") + signalregions[js];
		string mapname = histonames[n];
		TH1D* h_data = (TH1D*)histos[(histonames[n]+hsdata)]->Clone("h_data");
		TH1D* h_mc_sum = (TH1D*)histos[(histonames[n]+hsmc)]->Clone("h_mc");
		vector<TH1D*> h_signals; h_signals.clear();
		if(T2tt||T2bb) h_signals.push_back((TH1D*)histos[(histonames[n]+hssusy1)]->Clone("h_sig1"));
		if(T2tt||T2bb) h_signals.push_back((TH1D*)histos[(histonames[n]+hssusy2)]->Clone("h_sig2"));
		if(T2tt||T2bb) h_signals.push_back((TH1D*)histos[(histonames[n]+hssusy3)]->Clone("h_sig3"));
	 	TString xtitle = h_mc_sum->GetXaxis()->GetTitle();
		TString ytitle = h_mc_sum->GetYaxis()->GetTitle();
		THStack *h_stack = (THStack*)stacks[histonames[n]+hsstack ]->Clone("h_stack");;
		int njets, nbjets;
		if(js<2) njets = 2;
		else if(js<5) njets = 35;
		else if(js<8) njets = -6;
		else njets = -3;
		if(js==0||js==2||js==5) nbjets = 0;
		else if(js==1) nbjets = -1;
		else if(js==3||js==6) nbjets = 1;
		else if(js==4||js==7) nbjets = 2;
		else nbjets = -3;
		int nleps = -10;
		int n_sig = h_signals.size();
		TString outname = histonames[n]+hsstack;
		bool logflag = true;
		double maxi = 0;
		for(int nbin = 1; nbin<= h_mc_sum->GetNbinsX(); ++nbin){
			if(h_mc_sum->GetBinContent(nbin) > maxi) maxi = h_mc_sum->GetBinContent(nbin);
		}
		maxi = maxi * 2.;
		h_stack->SetMaximum(maxi);
		h_mc_sum->SetMaximum(maxi);
		h_data->SetMaximum(maxi);
		if(T2tt||T2bb) h_signals[0]->SetMaximum(maxi);
		if(T2tt||T2bb) h_signals[1]->SetMaximum(maxi);
		if(T2tt||T2bb) h_signals[2]->SetMaximum(maxi);
		if(h_mc_sum->Integral()>0||h_data->Integral()>0){
	 		tA->printHisto(h_stack, h_data, h_mc_sum, h_signals, n_sig, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nbjets, nleps, 1);
		}

	}}

	//make other plots
	TCanvas* c1 = new TCanvas("dummy","",0,0,600,600);
	for(int il = 0; il<signalregionsize; ++il){
		TString canvname;
		c1->cd();
		taueff[il]->Draw();
		c1->Update();
		canvname = taueff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		hadtaueff[il]->Draw();
		c1->Update();
		canvname = hadtaueff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		tauwithinacceff[il]->Draw();
		c1->Update();
		canvname = tauwithinacceff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		hadtauwithinacceff[il]->Draw();
		c1->Update();
		canvname = hadtauwithinacceff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		matchedtaueff[il]->Draw();
		c1->Update();
		canvname = matchedtaueff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		matchedhadtaueff[il]->Draw();
		c1->Update();
		canvname = matchedhadtaueff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		matchedtauwithinacceff[il]->Draw();
		c1->Update();
		canvname = matchedtauwithinacceff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		matchedhadtauwithinacceff[il]->Draw();
		c1->Update();
		canvname = matchedhadtauwithinacceff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
	}

}

//this is a copy from the MassPlotter.cc function - no detailed comments
void Make1DPlotsRatio(THStack *hstack, TH1D *histdata, TH1D* histmc, vector<TH1D*> histsusy, int nsig,  TLegend *leg,  TString outname, TString outputdir, Option_t *drawopt,  bool logFlag, TString xtitle, TString ytitle, int njets, int nbjets, int nleps, float overlayScale){


	// define canvas and pads 
	TH1D *h1 = (TH1D*)histmc->Clone("h1_copy");
	TH1D *h2 = (TH1D*)histdata->Clone("h2_copy");

	h1->SetTitle("");
	h2->SetTitle("");
	hstack->SetTitle("");

	h1->SetStats(0);	
	h2->SetStats(0);	
	h1->SetMarkerStyle(20);
	h2->SetMarkerStyle(20);

	TCanvas* c1 = new TCanvas(outname+"_ratio","",0,0,600,600);
	c1->SetFrameLineWidth(1);
	c1 -> cd();

 	TPad *p_plot  = new TPad(outname+"_plotpad",  "Pad containing the overlay plot", 0,0.211838,1,1);
	p_plot->SetLeftMargin(0.131579);
	p_plot->SetRightMargin(0.08);
	p_plot->SetTopMargin(0.06895515);
	p_plot->SetBottomMargin(0.07206074);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad(outname+"_ratiopad", "Pad containing the ratio",   0,0.01863354,0.9967105,0.2189441);
	p_ratio->SetLeftMargin(0.1336634);	
	p_ratio->SetRightMargin(0.075);
	p_ratio->SetTopMargin(0.06976745);
	p_ratio->SetBottomMargin(0.2790698);

	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	if(logFlag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logFlag) max = 2.5*max;
	else max = 1.5*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	hstack->SetMaximum(max);
	hstack->SetMinimum(0.02);
	hstack->Draw("hist");
	hstack->GetYaxis()->SetTitle(ytitle);
	hstack->GetYaxis()->SetLabelSize(0.05);
	hstack->GetYaxis()->SetTitleSize(0.05);
	hstack->GetYaxis()->SetTitleOffset(1.3);

	h2    ->Draw("sameE");
	for(int n=0; n<nsig;++n){
		histsusy[n]->Scale(overlayScale ? overlayScale : h2->Integral() / histsusy[n]->Integral());
		histsusy[n]->SetFillColor(0);
		histsusy[n]->SetLineStyle(kDotted);
		histsusy[n]->Draw("samehist");
	}

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.03);

	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d jets",abs(njets)) : TString::Format("%d jets",abs(njets));
	text += nbjets==-10 ? "" : nbjets < 0 ? TString::Format(", #geq %d b-tag",abs(nbjets)) : TString::Format(", %d b-tag",abs(nbjets));
	text += nleps == 1 ? ", 1 lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());//standard
	TLatex LumiBox;
	LumiBox.SetNDC();
	LumiBox.SetTextSize(0.0305);
	TString lumi = outname;
	LumiBox.DrawLatex(0.45,0.943,outname.Data());//standard
	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
	//gPad->SetLogy();
 
	TH1D *h_ratio = (TH1D*)histdata->Clone("h2_copy_2");	
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);

 	h_ratio ->Divide(h2, h1);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	//MC with errors
	TH1D*h_ratio_mc = (TH1D*)histmc->Clone("h1_copy_2");
	h_ratio_mc->Divide(h1);
	h_ratio_mc->GetYaxis()->SetRangeUser(0,2);
	h_ratio_mc->GetXaxis()->SetLabelSize( 0.);
	h_ratio_mc->GetYaxis()->SetTitle("Data / MC");
	h_ratio_mc->GetXaxis()->SetTitle(xtitle);
	h_ratio_mc->GetXaxis()->SetTitleSize(0.2);
	h_ratio_mc->GetXaxis()->SetTitleOffset(0.5);
	h_ratio_mc->GetYaxis()->SetLabelSize(0.19);
	h_ratio_mc->GetXaxis()->SetTickLength(0.09);
	h_ratio_mc->GetYaxis()->SetTitleSize(0.18);
	h_ratio_mc->GetYaxis()->SetTitleOffset(0.36);
	h_ratio_mc->GetYaxis()->SetNdivisions(509);//xxxnew

	h_ratio   ->SetTitle("");
	h_ratio_mc->SetTitle("");
	
	h_ratio_mc->SetFillStyle(3001);
	h_ratio_mc->Draw("E2");
	h_ratio ->DrawCopy("Esame");//LEO MOD
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=outname+"_ratio";
	Util::Print(c1, save, outputdir);

}

//copy from a MassPlotter.cc function - no detailed comments
void printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, vector<TH1D*> h_sig, int nsig, TLegend* legend,  TString canvname, TString outputdir, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale, double maxi){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 600, 600);
	col->SetRightMargin(0.08);
	col->SetLeftMargin(0.18);
	col->SetTopMargin(0.07);
	col->SetBottomMargin(0.17);

	TLegend *leg = (TLegend*) legend->Clone("leg");

	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h        -> SetMinimum(0.05);
		h_mc_sum -> SetMinimum(0.05);
		h_data   -> SetMinimum(0.05);
		for (int i=0; i<nsig; i++)  h_sig[i]-> SetMinimum(0.05);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h_mc_sum->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h_mc_sum->SetMaximum(max);
	h       ->SetMaximum(max);
	h_data  ->SetMaximum(maxi);
	h_mc_sum->SetMaximum(maxi);
	h       ->SetMaximum(maxi);

	h->Draw(drawopt);
	if(h_data->Integral()>0) {
		h_data       ->Draw("sameE");
	}
	for (int i=0; i<nsig; i++) {
	  h_sig[i]->Scale(overlayScale ? overlayScale : h_data->Integral() / h_sig[i]->Integral());
	  h_sig[i]->SetLineStyle(2);
	  h_sig[i]->SetFillColor(0);
	  h_sig[i]->Draw("samehist");
	}
	if(leg != NULL ){
		leg -> SetY1NDC(0.68);
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.0305);
	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d jets",abs(njets)) : TString::Format("%d jets",abs(njets));
	text += nbjets==-10 ? "" : nbjets < 0 ? TString::Format(", #geq %d b-tag",abs(nbjets)) : TString::Format(", %d b-tag",abs(nbjets));
	text += nleps == 1 ? ", 1 lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());//standard
	TLatex LumiBox;
	LumiBox.SetNDC();
	LumiBox.SetTextSize(0.0305);
	TString lumi = canvname;
	LumiBox.DrawLatex(0.45,0.943,lumi.Data());//standard
	h->GetXaxis()->SetTitle(xtitle);
	h->GetXaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetTitleOffset(1.1);
	
	stringstream yTitle;
	yTitle << ytitle.Data();
	h->GetYaxis()->SetTitle(yTitle.str().c_str());
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleOffset(1.47);

	gPad->RedrawAxis();
	Util::PrintNoEPS(col, canvname, outputdir);
	Util::PrintEPS(col, canvname, outputdir);

}

//reads the samples.dat
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

			if(s.type!="data"){
			TH1F *h_PUWeights = (TH1F*) s.file->Get("h_PUWeights");
			TH1F *h_Events    = (TH1F*) s.file->Get("h_Events");
			if(h_PUWeights==0 || h_Events==0){
				cout << "ERROR: sample " << (s.file)->GetName() << " does not have PU and NEvents histos! " << endl;
				exit(1);
			}
			s.type!="data" ? s.PU_avg_weight = h_PUWeights->GetMean()    : s.PU_avg_weight =1;
			s.type!="data" ? s.nevents       = h_Events   ->GetEntries() : s.nevents       =1;

			// DON'T DO THIS AT HOME !!!!
			if ( s.name == "T1tttt_mGlu-1000_mLSP-400" ) s.nevents = 20000;
			if ( s.name.Contains("T2bb") ) s.nevents = 10000;
			if ( s.name.Contains("T2tt") ) s.nevents = 50000;
			// DON'T DO THAT AT HOME !!!!

			delete h_PUWeights;
			delete h_Events;
			} else {
			s.PU_avg_weight =1;
			s.nevents       =1;
			}
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

