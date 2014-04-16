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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path of MT2tree.hh
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/BTagWeight.hh"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/Stops/SolveTTbar.h"//use this for D_2 calculation
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/Stops/SolveTTbar.c"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/Stops/SolveTTbarNew.h"//alternative calculation of D_2 calculation, however less stable
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/Stops/SolveTTbarNew.C"

using namespace std;

void load(const char* filename = "samples_2141_dataonly.dat");
void StopFeasibility();
void MakeCorrelationPlots(map<string, TH2D*> histos2, const int leptontypesize, string *lepton_type, vector<string> histonames, TString outputdirectory);
void MakePlots(map<string, TH1D*> histos, const int leptontypesize, string *lepton_type, map<string, THStack*> stacks, TLegend *Legend1, vector<string> histonames);
void Make3DPlots(map<string, TH3D*> histos3, const int leptontypesize, string *lepton_type, vector<string> histonames, TString outputdirectory);
void Make1DPlotsRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, TH1D* histsusy, bool logFlag, bool normalize, TString outname, TLegend *leg, TString xtitle, TString ytitle, float overlayScale);
void Make1DPlotsNoRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, TH1D* histsusy, bool logFlag, bool normalize, TString outname, TLegend *leg, TString xtitle, TString ytitle, float overlayScale);


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

//in this code the calculation of the BTV weights is done interactively - not needed for 'new' (by now old) defintion of MT2trees
TString btagging_file = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/Efficiencies/Stop_BEfficiencies_PTbinned_nocuts_SSVHEM.root";
TFile*  btagefffile;
TH1D*   hbeff;
TH1D*   hceff;
TH1D*   hleff;
TString tagger = "CSVM";// SSVHPT or SSVHEM
//adds some flexibility for jet and b-jet selection
int     njets  = -2;
int     ntags  = -2;//- means >=, while + means ==; -99 means no tagging requirement
int     Tagger = 4;//TCHE:0 TCHP:1 SSVHE:2 SSVHP:3, CSV: 4
float   discr  = 0.679;
float   bpt    = 40.;//add in cutStream
float   beta   = 2.4;
float   jpt    = 40.;
float   jeta   = 2.4;

TString samples                   = "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/samples/samples_Stop_53X_0bMET30_corr.dat";//samples.dat to load
Bool_t  runData                   = true; // if you do only MCClosure you can set this to false
Bool_t  calcsusy                  = true; // run over susy sample - default is true
TString outputdir                 = "Stops/test53X_corr_Stop400bChi300_SlepSnu200";//directory where histogram plots are stored, defines also root file of all histograms
TString outputname                = "DiLeptonicStops.root";//default
Bool_t  dofastbtagSFreweighting   = false;// do not calculate correct BTV SF weights, but some average weights - makes code faster - not needed for 'new' (by now old) defintion of MT2trees
Bool_t  logflag                   = true; // plots have log style one y-axis
Bool_t  fSave                     = true; // save plots
Bool_t  debug                     = false;//don't run over mc to save time
Bool_t  D2calc                    = false;//do the D2 calculation - default = true, note however that this takes a lot of time
Bool_t  D2calcLuc                 = true; //use Luc's implementation (robust, but slower) of D2 calculation, needs D2calc=true, default = true
Bool_t  D2calcMinuit              = false;//use my Minuit implementation(fast, but less robust) of D2 calculation, needs D2calc=true, default = false
Bool_t  saveD2                    = false;//make this false if D2calc = false, otherwise get nan histograms
Bool_t  debugD2                   = false;//print out some addtional information on D2 calculation to debug the two calculations by comparing them to each other
Bool_t  plotonlywithratio         = true; //make only plots complemented by a ratio plot, otherwise also make the plots without ratio, default = true
Bool_t  plotonlyLL                = true; //make only plots with LL(==ee+mumu+emu), default = true
Bool_t  plotsingletop             = false;//for 2d correlation plots, do a separate version for single top, default = false
//the next for flags depend on chosen stop samples in the samples.dat: As I do genlevel plots, these plots depend on decay chain of the stop
Bool_t  stopTotopLSP              = false;//stop --> top + LSP
Bool_t  stopTobChargino           = true; //stop --> b + chargino
Bool_t  stopTobCharginoWLSP       = true; //chargino --> W + LSP, needs stopTobChargino = true
Bool_t  stopTobCharginoslepsnu    = false;//chargino --> l + snu, or slep + nu (implementation not complete though), needs stopTobChargino = true
Bool_t  doplotting                = true; // default = true, make plots (otherwise the histograms are only stored in the root file)
Bool_t  ftype1met                 = false;//set met explicitly to type1 met, default = false (note that the default MET is anyway type1met - so usually this flag does not do anything)
Bool_t  frawmet                   = false;//set met explicitly to raw pf met, default = false, used only if type1met == false

//this function produces various tons of plots used for performing
//my dileptonic stop feasibility study
//note that we decided to stop this due to manpower issues ;-(
void StopFeasibility(){

gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

   if(!saveD2) outputname = (TString)"WithoutD2_"+outputname;

  //histograms for BTV SF calculation
  btagefffile	= TFile::Open(btagging_file);
  hbeff		= (TH1D*)btagefffile->Get("h_beff_NJets40ge2");//also NJets40eq2
  hceff		= (TH1D*)btagefffile->Get("h_ceff_NJets40ge2");
  hleff		= (TH1D*)btagefffile->Get("h_leff_NJets40ge2");

    Util::MakeOutputDir(outputdir);
    map<string, TH1D*>    histos;
    map<string, TH2D*>    histos2;
    map<string, TH3D*>    histos3;
    map<string, THStack*> hstacks;
    map<string, THStack*> stacks;
    TLegend *Legend1 = new TLegend(.71,.54,.91,.92);
    vector<string> histonames; histonames.clear();
    TLegend* legend = new TLegend(.71,.54,.91,.92);
    legend -> SetFillColor(0);
    legend -> SetBorderSize(0);

	//event selection
    std::ostringstream baseStream;
    std::ostringstream cutStream;
    TString basecuts;
	baseStream << " " 
    << "misc.MET>=30"                                                      << "&&"//minimal MET for MT2
    << "misc.PassJetID ==1"                                                << "&&"
    << "NJetsIDLoose>=2"                                                   << "&&"
    << "(NEles+NMuons)>=2"                  << "&&"//veto secondary leptons
    // Noise
    << "(misc.isFastSim || misc.HBHENoiseFlag == 0)"                       << "&&"
    << "(misc.isFastSim || misc.CSCTightHaloIDFlag == 0)"                  << "&&"
    << "(misc.isFastSim || misc.hcalLaserEventFlag == 0)"                  << "&&"
    << "misc.NegativeJEC==0"                                               << "&&"
    << "misc.trackingFailureFlag==0"                                       << "&&"
    << "misc.eeBadScFlag==0"                                               << "&&"
    << "misc.EcalDeadCellTriggerPrimitiveFlag==0"                          << "&&"
    << "misc.CrazyHCAL==0";
	basecuts = baseStream.str().c_str();
	//Analysis cuts
	cutStream << " " 
    << "(NEles+NMuons)==2"                 << "&&"
    << "NBJets40CSVM>=2";
    if(jpt<40)      cutStream << "&& NJetsIDLoose  >=2";
    else if(jpt<50) cutStream << "&& NJetsIDLoose40>=2";
    else            cutStream << "&& NJetsIDLoose50>=2";


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

    const int leptontypesize = 5;
    string lepton_type[leptontypesize] = {"MuMu", "EMu", "EE", "SF", "LL"};

	cout << " initializing histograms ...";//including axis title - again: really tons of plots, I guess by now one could kill a lot of them, but I leave that open
	vector<string> vs; vs.clear();
	bool av = true;//append vector
	for(int is = 0; is<sampletypesize; ++is){
           for(int il = 0; il<leptontypesize; ++il){
		string hs = string("_") + lepton_type[il] + string("_") + sample_type[is];
		string mapname;
		string title, xtitle, ytitle, ztitle;
		mapname = "LepPt"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Lepton p_{T}"; xtitle = "Lepton p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, 0, 600);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "LepMT"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Lepton M_{T}"; xtitle = "Lepton M_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "N jet p_{T}, the ones which make M_{T2}(lb)"; xtitle = "b-p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "AllBPt"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "B jet p_{T}, all b jets"; xtitle = "b-p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NJetsID20"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of jets(p_{T}>20 GeV, ;oose JID)"; xtitle = "#jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NJetsID30"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of jets(p_{T}>30 GeV, loose JID)"; xtitle = "#jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NJetsID40"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of jets(p_{T}>40 GeV, loose JID)"; xtitle = "#jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NJetsID50"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of jets(p_{T}>50 GeV, loose JID)"; xtitle = "#jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 15, 0, 15);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NBJets20"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of b jets(p_{T}>20 GeV)"; xtitle = "#b jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 5, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NBJets30"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of b jets(p_{T}>30 GeV)"; xtitle = "#b jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 5, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NBJets40"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of b jets(p_{T}>40 GeV)"; xtitle = "#b jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 5, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "NBJets50"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of b jets(p_{T}>50 GeV)"; xtitle = "#b jets";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 5, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "AllBTagDiscr"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "B-tagging discriminant"; xtitle = "b-tag discr";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 1.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BTagDiscr"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "B-tagging discriminant - jets fot M_{T2}(lb)"; xtitle = "b-tag discr";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 1.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "HT"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "H_{T}"; xtitle = "H_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 85, 0, 1275);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "E_{T}^{miss}"; xtitle = "E_{T}^{miss} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 35, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb)"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb)"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlbcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_withMlbcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2(lb) after M(lb) cut"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MlbAll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant mass of leptons and b jets"; xtitle = "M_{lb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant mass of leptons and b jets making M_{T2}(lb)"; xtitle = "M_{lb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant dilepton mass"; xtitle = "M_{ll} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 105, 0, 1050);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mbb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant di-b jet mass"; xtitle = "M_{bb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 1050);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MbbAll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant di-b jet mass - all combinations"; xtitle = "M_{bb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 1050);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l)"; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_withMlbcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l) after M(lb) cut"; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(b)"; xtitle = "M_{T2}(b) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb_lInMET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(b)"; xtitle = "M_{T2}(b) [GeV], l's in E_{T}^{miss}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb_lInMET_mW"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(b), l's in E_{T}^{miss}, testmass 80.4"; xtitle = "M_{T2}(b) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb_lInMET_mW_withMlbcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(b), l's in E_{T}^{miss}, testmass 80.4 after M(lb) cut"; xtitle = "M_{T2}(b) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "JetPt"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "jet p_{T} of all jets"; xtitle = " jet p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "UTM"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Upstream transverse momentum"; xtitle = "UTM [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 500);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DPhi_lbAll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "#Delta#phi_{lb}"; xtitle = "#Delta#phi_{lb}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DPhi_lb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "#Delta#phi_{lb} making M_{T2}(lb)"; xtitle = "#Delta#phi_{lb}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DPhi_ll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "#Delta#phi_{ll}"; xtitle = "#Delta#phi_{ll}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DPhi_bb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "#Delta#phi_{bb} for M_{T2}(lb)"; xtitle = "#Delta#phi_{bb}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DPhi_bbAll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "#Delta#phi_{ll}, all combinations"; xtitle = "#Delta#phi(b,b)";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DR_lbAll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "#Delta R_{lb}"; xtitle = "#Delta R_{lb}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DR_lb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "#Delta R_{lb} making M_{T2}(lb)"; xtitle = "#Delta R_{lb}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DR_ll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "#Delta R_{ll}"; xtitle = "#Delta R_{ll}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DR_bb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "#Delta R_{bb} for M_{T2}(lb)"; xtitle = "#Delta R_{bb}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DR_bbAll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "#Delta R_{bb}, all combinations"; xtitle = "#Delta R_{bb}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_genmatchingTop"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{lb}, genflavour-matching Top"; xtitle = "M_{lb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2_genmatchingTops"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb), genflavour-matching Tops"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMlb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{lb} - genleveltruth"; xtitle = "M_{lb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMT2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) - genleveltruth"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMT2_withMlbcut"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut - genleveltruth"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMlb_allAcceptance"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{lb} - genleveltruth, full acceptance"; xtitle = "M_{lb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMT2_allAcceptance"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) - genleveltruth, full acceptance"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMT2_withMlbcut_allAcceptance"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut - genleveltruth, full acceptance"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMtop_allAcceptance"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{top} - genleveltruth, full acceptance"; xtitle = "M_{bl#nu} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "GenMtop"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{top} - genleveltruth"; xtitle = "M_{bl#nu} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlbcut_PUle5"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) for #PV<=5 after Mlb cut"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlbcut_PUge9"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) for #PV>=9 after Mlb cut"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlbcut_woZveto"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) including Z M_{ll} window after M_{lb} cut"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_woZveto"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) including Z M_{ll} window"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_woZveto"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) including Z M_{ll} window"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb_woZveto"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(b) including Z M_{ll} window"; xtitle = "M_{T2}(b) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb_linMET_woZveto"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(b), l in E_{T}^{miss} including Z M_{ll} window"; xtitle = "M_{T2}(b) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_woZveto"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l) including Z M_{ll} window"; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "No_PV"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Number of primary vertices"; xtitle = "#vertices";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 40, 0, 40);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Mass fit discriminant D_{2}"; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_withMlb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Mass fit discriminant after M_{lb} cut"; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());

		mapname = "MT2lb_withMlbcut_afterMT2llge85"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_afterMT2llge85"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_afterMT2llge85"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlbcut_afterMT2bb_linMETge150"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_afterMT2bb_linMETge150"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_afterMT2bb_linMETge150"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_afterMT2bb_linMETge150"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_afterMT2lbge200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_afterMT2lb_withMlbcutge180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_afterMT2lb_masslessge140"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb_linMET_afterMT2llge85"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(b) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb_linMET_afterMT2lbge200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(b) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb_linMET_afterMT2lb_withMlbcutge180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(b) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2bb_linMET_afterMT2lb_masslessge140"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(b) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_afterMT2llge85"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 20, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_afterMT2bb_linMETge150"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 20, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_afterMT2lbge200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 20, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_afterMT2lb_withMlbcutge180"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 20, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_afterMT2lb_masslessge140"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 20, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());

		//new 14/06/2012
		mapname = "MT2lb_bothbtags_truebjets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_onebtag_truebjet"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_nonebtags_truebjets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_bothbtags_matchedtobfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_onebtag_matchedtobfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_nonebtags_matchedtobfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "RecoPt_over_GenPartonPt_truebfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "Reco-p_{T} / genparton-p_{T}";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 40, 0, 4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "RecoPt_over_GenPartonPt_truebfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "Reco-p_{T} / genparton-p_{T}";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 40, 0, 4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_bothbtags_truebjets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_onebtag_truebjet"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_nonebtags_truebjets"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_bothbtags_matchedtobfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_onebtag_matchedtobfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_nonebtags_matchedtobfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-p_{T} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-p_{T} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discr";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 1.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discr";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 1.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "TrueBfromtopEta_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "TrueBfromtopEta_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "TrueBfromtopPt_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-p_{T} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "TrueBfromtopPt_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-p_{T} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname  = "PFMEToverGenMET_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "gen-E_T^{miss} / PF-E_T^{miss}";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 40, 0, 4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname  = "PFMEToverGenMET_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "gen-E_T^{miss} / PF-E_T^{miss}";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 40, 0, 4);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname  = "DeltaGenMET_PFMET_Phi_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta#phi(gen-E_{T}^{miss},PF-E_{T}^{miss})";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname  = "DeltaGenMET_PFMET_Phi_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta#phi(gen-E_{T}^{miss},PF-E_{T}^{miss})";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 50, 0, 5.);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "AltMT2lb_MT2lb_gt_200_usingtruebpartons"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 25, 0, 500);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaR_bjet_matchedtrueparton_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta R";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 13);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaR_bjet_matchedtrueparton_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta R";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 13);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaPt_bjet_matchedtrueparton_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta p_{T} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 350);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaPt_bjet_matchedtrueparton_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta p_{T} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 350);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_bothleptons_true_trueb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_onelepton_true_trueb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_nolepton_true_trueb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_bothb_true_truel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_oneb_true_truel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_nob_true_truel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		//new week 24/07/2012
		mapname = "masslessMT2lb_bothleptons_true_trueb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "masslessMT2lb_onelepton_true_trueb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "masslessMT2lb_nolepton_true_trueb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "masslessMT2lb_bothb_true_truel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "masslessMT2lb_oneb_true_truel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "masslessMT2lb_nob_true_truel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2l_bothleptons_true_trueb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2l_onelepton_true_trueb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2l_nolepton_true_trueb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2l_bothb_true_truel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2l_oneb_true_truel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2l_nob_true_truel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_bothleptons_true_trueb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_onelepton_true_trueb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_nolepton_true_trueb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_bothb_true_truel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_oneb_true_truel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_nob_true_truel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{lb} [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());

		//new in week 18/06/2012
		mapname = "DeltaR_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta R";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "DeltaR_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#Delta R";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 65, 0, 6.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discr";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 1.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discr";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 1.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discriminant";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 1.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discriminant";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 1.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_bothbtags_matchedtosamebfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_matchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discriminant";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 75, 0, 1.5);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_matchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_matchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BDiscr_matchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs;  xtitle = "b-tag discr";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 7);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BEta_matchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "BPt_matchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Eta_truebfromtop_matchedwithreco_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Pt_truebfromtop_matchedwithreco_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Eta_truebfromtop_notmatchedwithreco_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Pt_truebfromtop_notmatchedwithreco_MT2lb_lt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Eta_truebfromtop_matchedwithreco_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Pt_truebfromtop_matchedwithreco_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Eta_truebfromtop_notmatchedwithreco_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Pt_truebfromtop_notmatchedwithreco_MT2lb_gt_200"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());

	//new single top checks 05/07/2012
		mapname = "SingleTop_BPt_truebfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_BPt_otherb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_LPt_truelfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "Lepton p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_LPt_otherl"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "Lepton p_{T} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_BEta_truebfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_BEta_otherb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "b-#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_LEta_truelfromtop"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_LEta_otherl"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "#eta";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 60, -6, 6);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_LID_otherl"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "ID";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 30);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_LMID_otherl"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "mother-ID";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 30);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_LGMID_otherl"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "grandmother-ID";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 30);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_BID_otherb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "ID";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 30);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_BMID_otherb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "mother-ID";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 30);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_BGMID_otherb"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "grandmother-ID";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 30);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_MT2lb_tchannel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_MT2lb_tWchannel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "SingleTop_MT2lb_schannel"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());

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
		//new single top 05/07/2012
		mapname = "LOriginID"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Matched lepton ID (unmatched==0)"; xtitle = "L1 ID"; ytitle = "L2 ID";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 60, -30, 30, 60, -30, 30);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "LOriginMID"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Matched lepton mother-ID (unmatched==0)"; xtitle = "L1 MID"; ytitle = "L2 MID";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 60, -30, 30, 60, -30, 30);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "LOriginGMID"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Matched lepton grandmother-ID (unmatched==0)"; xtitle = "L1 GMID"; ytitle = "L2 GMID";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 60, -30, 30, 60, -30, 30);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "BOriginID"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Matched b-tag ID (unmatched==0)"; xtitle = "B1 ID"; ytitle = "B2 ID";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 60, -30, 30, 60, -30, 30);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "BOriginMID"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Matched b-tag mother-ID (unmatched==0)"; xtitle = "B1 MID"; ytitle = "B2 MID";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 60, -30, 30, 60, -30, 30);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "BOriginGMID"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Matched b-tag grandmother-ID (unmatched==0)"; xtitle = "B1 GMID"; ytitle = "B2 GMID";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 60, -30, 30, 60, -30, 30);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "LOrigin"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Matched lepton identification (unmatched==0)"; xtitle = "Lepton identification"; ytitle = "PID";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 4, 0, 4, 60, -30, 30);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		histos2[mapname]->GetXaxis()->SetBinLabel(1, "ID");   histos2[mapname]->GetXaxis()->SetBinLabel(2, "MID"); histos2[mapname]->GetXaxis()->SetBinLabel(3, "GMID"); histos2[mapname]->GetXaxis()->SetBinLabel(4, "recoID");
		mapname = "BOrigin"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Matched b-tag identification (unmatched==0)"; xtitle = "b-tag identification"; ytitle = "PID";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 4, 0, 4, 60, -30, 30);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		histos2[mapname]->GetXaxis()->SetBinLabel(1, "ID");   histos2[mapname]->GetXaxis()->SetBinLabel(2, "MID"); histos2[mapname]->GetXaxis()->SetBinLabel(3, "GMID"); histos2[mapname]->GetXaxis()->SetBinLabel(4, "PAT ID");
		mapname = "MT2lb_vs_UTM"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) vs UTM"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "UTM [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 50, 0, 500);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlb_vs_UTM"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut vs UTM"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "UTM [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 50, 0, 500);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_vs_BTagDiscr"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs b-tagging discriminant"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "b-tag discr";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 75, 0, 1.5);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_vs_MT2ll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_vs_MT2bb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs M_{T2}(b)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(b) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_vs_MT2bb_lInMET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs M_{T2}(b)(l in E_{T}^{miss})"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(b) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_vs_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs M_{T2}(b)(l in E_{T}^{miss}), M_{W}"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(b) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_vs_MET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs E_T^{miss}"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "E_{T}^{miss} [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 50, 0, 1000);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_vs_HT"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs H_T"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "H_{T} [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 85, 0, 1275);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_vs_PV"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs #PV"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "#vertices";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 25, 0, 25);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_withMlb_vs_MT2ll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_vs_D2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs mass fit discriminant D_{2}"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "D_{2}";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 30, 0, 300);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_vs_Mlb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}lb vs M_{lb}"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{lb} [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 90, 0, 900);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_vs_BTagDiscr"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) vs b-tagging discriminant"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "b-tag discr";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 75, 0, 1.5);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_vs_MT2ll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_vs_MT2bb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) vs M_{T2}(b)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(b) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_vs_MT2bb_lInMET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) vs M_{T2}(b)(l in E_{T}^{miss})"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(b) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_vs_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) vs M_{T2}(b)(l in E_{T}^{miss}), M_{W}"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(b) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_vs_MET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) vs E_{T}^{miss}"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "E_{T}^{miss} [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 50, 0, 1000);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_vs_HT"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) vs H_{T}"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "H_{T} [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 85, 0, 1275);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_vs_PV"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) vs #PV"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "#vertices";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 25, 0, 25);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_vs_D2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) vs mass fit discriminant D_{2}"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "D_{2}";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 30, 0, 300);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_vs_Mlb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) vs M_{lb}"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{lb} [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlbcut_vs_MT2ll"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlbcut_vs_MT2bb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut  vs M_{T2}(b)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(b) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlbcut_vs_MT2bb_lInMET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut  vs M_{T2}(b) (l in E_{T}^{miss})"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(b) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlbcut_vs_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_lb} cut vs M_{T2}(b) (l in E_{T}^{miss}), M_{W}"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(b) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlbcut_vs_MET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut vs E_{T}^{miss}"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "E_{T}^{miss} [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 50, 0, 1000);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlbcut_vs_HT"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut vs H_{T}"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "H_{T} [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 85, 0, 1275);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlbcut_vs_PV"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut vs #PV"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "#vertices";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 25, 0, 25);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlbcut_vs_D2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut vs mass fit discriminant D_{2}"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "D_{2}";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 30, 0, 300);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlbcut_vs_MT2lb_massless"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut vs massless M_{T2}(lb)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "massless M_{T2}(lb) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_vs_MT2lb_massless"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) vs massless M_{T2}(lb)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "massless M_{T2}(lb) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2ll_vs_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l) vs M_{T2}(b), l's in E_{T}^{miss}, M_{W}"; xtitle = "M_{T2}(l) [GeV]"; ytitle = "M_{T2}(b) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 140, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2ll_vs_MT2bb_lInMET"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l) vs M_{T2}(b), l's in E_{T}^{miss}"; xtitle = "M_{T2}(l) [GeV]"; ytitle = "M_{T2}(b) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 140, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2ll_vs_MT2bb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l) vs M_{T2}(b)"; xtitle = "M_{T2}(l) [GeV]"; ytitle = "M_{T2}(b) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 140, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2ll_vs_D2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l) vs D_{2}"; xtitle = "M_{T2}(l) [GeV]"; ytitle = "D_{2}";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 140, 0, 700, 30, 0, 300);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2ll_vs_Mlb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l) vs M_{lb}"; xtitle = "M_{T2}(l) [GeV]"; ytitle = "M_{lb} [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 140, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2bb_vs_D2"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(b) vs D_{2}"; xtitle = "M_{T2}(b) [GeV]"; ytitle = "D_{2}";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 30, 0, 300);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2bb_vs_Mlb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(b) vs M_{lb}"; xtitle = "M_{T2}(b) [GeV]"; ytitle = "M_{lb} [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2bb_lInMET_vs_Mlb"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(b) (l in E_{T}^{miss}) vs Mlb"; xtitle = "M_{T2}(b) [GeV]"; ytitle = "M_{lb} [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2bb_vs_BTagDiscr"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(b) vs b-tagging discriminant"; xtitle = "M_{T2}(b) [GeV]"; ytitle = "b-tag discr";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 75, 0, 1.5);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "Mlb_vs_BTagDiscr"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{lb} vs b-tagging discriminant"; xtitle = "M_{lb} [GeV]"; ytitle = "b-tag discr";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 75, 0, 1.5);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "D2_vs_BTagDiscr"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "D_{2} vs b-tagging discriminant"; xtitle = "D_{2}"; ytitle = "b-tag discr";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 30, 0, 300, 75, 0, 1.5);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "PUWeight_vs_PV"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "PU Weight vs #PV"; xtitle = "PU Weight"; ytitle = "#vertices";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 40, 0, 2, 25, 0, 25);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());

		mapname = "MT2lb_MT2ll_MT2bb"; if(av) vs.push_back(mapname); mapname = mapname+hs; title = "M_{T2}(lb) vs. M_{T2}(l) vs. M_{T2}(b)";
		xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]"; ztitle = "M_{T2}(b) [GeV]";
        	if(histos3.count(mapname) == 0 ) histos3[mapname] = new TH3D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700,70,0,700);
		histos3[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos3[mapname]->GetYaxis()->SetTitle(ytitle.c_str()); histos3[mapname]->GetZaxis()->SetTitle(ztitle.c_str());
		mapname = "MT2lb_MT2ll_MT2bb_lInMET"; if(av) vs.push_back(mapname); mapname = mapname+hs; title = "M_{T2}(lb) vs. M_{T2}(l) vs. M_{T2}(b) with l's in E_{T}^{miss}";
		xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]"; ztitle = "M_{T2}(b) [GeV]";
        	if(histos3.count(mapname) == 0 ) histos3[mapname] = new TH3D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700,70,0,700);
		histos3[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos3[mapname]->GetYaxis()->SetTitle(ytitle.c_str()); histos3[mapname]->GetZaxis()->SetTitle(ztitle.c_str());
		mapname = "MT2lb_MT2ll_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname); mapname = mapname+hs; title = "M_{T2}(lb) vs. M_{T2}(l) vs. M_{T2}(b) with l's in E_{T}^{miss}, M_{W}";
		xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]"; ztitle = "M_{T2}(b) [GeV]";
        	if(histos3.count(mapname) == 0 ) histos3[mapname] = new TH3D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700,70,0,700);
		histos3[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos3[mapname]->GetYaxis()->SetTitle(ytitle.c_str()); histos3[mapname]->GetZaxis()->SetTitle(ztitle.c_str());
		mapname = "MT2lbMassless_MT2ll_MT2bb"; if(av) vs.push_back(mapname); mapname = mapname+hs; title = "massless M_{T2}(lb) vs. M_{T2}(l) vs. M_{T2}(b)";
		xtitle = "M_{T2}(lb,lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]"; ztitle = "M_{T2}(b) [GeV]";
        	if(histos3.count(mapname) == 0 ) histos3[mapname] = new TH3D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700,70,0,700);
		histos3[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos3[mapname]->GetYaxis()->SetTitle(ytitle.c_str()); histos3[mapname]->GetZaxis()->SetTitle(ztitle.c_str());
		mapname = "MT2lbMassless_MT2ll_MT2bb_lInMET"; if(av) vs.push_back(mapname); mapname = mapname+hs; title = "massless M_{T2}(lb) vs. M_{T2}(l) vs. M_{T2}(b) with l's in E_{T}^{miss}";
		xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]"; ztitle = "M_{T2}(b) [GeV]";
        	if(histos3.count(mapname) == 0 ) histos3[mapname] = new TH3D(mapname.c_str(), title.c_str(), 70, 0, 700, 70, 0, 700,70,0,700);
		histos3[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos3[mapname]->GetYaxis()->SetTitle(ytitle.c_str()); histos3[mapname]->GetZaxis()->SetTitle(ztitle.c_str());
		mapname = "MT2lbMassless_MT2ll_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname); mapname = mapname+hs; title = "massless M_{T2}(lb) vs. M_{T2}(l) vs. M_{T2}(b) with l's in E_{T}^{miss}, M_{W}";
		xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]"; ztitle = "M_{T2}(b) [GeV]";
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
        
	    //some variables needed later for D2 calculation
	    int ssb = 0; int osb = 0;
	    int vetoeddiscrcount = 0; int discrcount = 0; int globsol = 0;
	    int statvec[3]; statvec[0]=0; statvec[1]=0; statvec[2]=0;
	    int nstptot = 0; int ndertot = 0; int nrottot = 0; int nscntot = 0; int globalsolution = 0;
   	    if(runData ==false && fSamples[i].type=="data") continue;
	    if(calcsusy==false && fSamples[i].type=="susy") continue;
	    if(debug           && fSamples[i].type=="mc"  ) continue;
        
		//get sample type
	    string sampletype = (string)fSamples[i].type;
	    string leptype = "LL";
	    if(sampletype==(string)"mc"){
		     if(fSamples[i].sname=="Wtolnu") sampletype = (string)"WJets";
		else if(fSamples[i].sname=="DY"    ) sampletype = (string)"ZJets";
		else if(fSamples[i].name =="TTbar" ) sampletype = (string)"TTbar";
		else if(fSamples[i].sname=="TTbarV") sampletype = (string)"TTbarV";
		else if(fSamples[i].sname=="Top"   ) sampletype = (string)"SingleTop";//no ttbar
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
        
	    if(      fSamples[i].type=="data" && fSamples[i].sname=="EE-Data"  ) { myCuts += " && " + triggerEE;   } //cuts to be aplied only on data
	    else if( fSamples[i].type=="data" && fSamples[i].sname=="EMu-Data" ) { myCuts += " && " + triggerEMu;  } //cuts to be aplied only on data
	    else if( fSamples[i].type=="data" && fSamples[i].sname=="MuMu-Data") { myCuts += " && " + triggerMuMu; } //cuts to be aplied only on data
	    else if( fSamples[i].type=="data") {cout << "data not usuable" << " type " << fSamples[i].type << " Sname " << fSamples[i].sname << endl; continue; }//not usuable
        
            if(fSamples[i].type=="susy") cout << "no doing " << sampletype << endl;

   	    cout << "Cuts for Flow: " << myCuts << endl;
   	    fSamples[i].tree->Draw(">>selList", myCuts);

	    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

	    fSamples[i].tree->SetEventList(myEvtList);

	    int counter=0;
	    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
	    if(myEvtList->GetSize()==0) continue;
	//run over selected events
        while(myEvtList->GetEntry(counter++) !=-1){	
      	    int jentry = myEvtList->GetEntry(counter-1);
            nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
            fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
            if ( fVerbose>2 && counter % 50000 == 0 )  cout << "+++ Proccessing event " << counter << endl;
	    //if((fSamples[i].sname!="VVV"&&fSamples[i].sname!="TTbarV" && fSamples[i].type!="susy"&&sampletype!="Stop")  && fMT2tree->misc.HBHENoiseFlag!=0     ) continue;//test
	    //if((fSamples[i].sname!="VVV"&&fSamples[i].sname!="TTbarV" && fSamples[i].type!="susy"&&sampletype!="Stop")  && fMT2tree->misc.hcalLaserEventFlag!=0) continue;//test
            Double_t weight = sample_weight;
            if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 

            Bool_t recoedee   = false;// exact 2 ele, 0 muo
            Bool_t recoedemu  = false;// exact 1 ele, 1 muo
            Bool_t recoedmumu = false;// exact 0 ele, 2 muo
		//could change this - might have some inefficiencies - three leptons with one not passing IDMedium==1
            if(fMT2tree->NEles>=2 &&fMT2tree->ele[0].lv.Pt()>20&&fMT2tree->ele[1].lv.Pt()>20&&fMT2tree->ele[2].lv.Pt()<10&&(fMT2tree->NMuons==0||fMT2tree->muo[0].lv.Pt()<10) && fMT2tree->ele[0].IDMedium==1 && fMT2tree->ele[1].IDMedium==1) recoedee   = true;
            if(fMT2tree->NMuons>=2&&fMT2tree->muo[0].lv.Pt()>20&&fMT2tree->muo[1].lv.Pt()>20&&fMT2tree->muo[2].lv.Pt()<10&&(fMT2tree->NEles ==0||fMT2tree->ele[0].lv.Pt()<10)) recoedmumu = true;
            if(fMT2tree->NMuons>=1&&fMT2tree->muo[0].lv.Pt()>20&&fMT2tree->muo[1].lv.Pt()<20&&fMT2tree->NEles>=1 &&fMT2tree->ele[0].lv.Pt()>20&&fMT2tree->ele[1].lv.Pt()<20 && fMT2tree->ele[0].IDMedium==1) recoedemu  = true;

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


	    if(!(recoedee)   && !(recoedemu)   && !(recoedmumu)   ) continue;//require dilepton
	    if(!(recoedosee) && !(recoedosemu) && !(recoedosmumu) ) continue;//require os-dilepton, maybe comment this for background est.
	    //if(!(recoedeenZ) ||                   !(recoedmumunZ) ) continue;//require of Zpeak --> done after Mll histo filling

	    // Mll < 10 GeV  cut - safety cut against low mass resonances
	    if(recoedee   && ((fMT2tree->ele[0].lv + fMT2tree->ele[1].lv).M() )<10) continue;
	    if(recoedemu  && ((fMT2tree->ele[0].lv + fMT2tree->muo[0].lv).M() )<10) continue;
	    if(recoedmumu && ((fMT2tree->muo[0].lv + fMT2tree->muo[1].lv).M() )<10) continue;

	    if(recoedee)   leptype = "EE";
	    if(recoedemu)  leptype = "EMu";
	    if(recoedmumu) leptype = "MuMu";

		//'trigger weights'
		//EEweight; EMuweight; MUMUweight, if different lumis;
		//taken from AN2011_466_v4
		// weight = triggerweight * weight * lumiweight
		// check triggerweight = 0.96 | 0.934 | 0.875
		// trigger efficiency = 0.948 | 0.919 | 0.941
	    if(recoedee   && !(fMT2tree->misc.isData) ) weight = 0.960/*0.948*/ * weight * 0.953366;//from SSb Florida
	    if(recoedemu  && !(fMT2tree->misc.isData) ) weight = 0.934/*0.919*/ * weight;//from SSb Florida
	    if(recoedmumu && !(fMT2tree->misc.isData) ) weight = 0.875/*0.941*/ * weight * 1.031213;//from SSb Florida
	    if(fMT2tree->misc.isData && recoedee   && fSamples[i].sname!="EE-Data"  ) continue;
	    if(fMT2tree->misc.isData && recoedemu  && fSamples[i].sname!="EMu-Data" ) continue;
	    if(fMT2tree->misc.isData && recoedmumu && fSamples[i].sname!="MuMu-Data") continue;

	    unsigned int NumBJets = fMT2tree->NBJets40CSVM;
		//do BTV SF weights - now done interactively - by now this is obsolete due to 'new' (actually old) MT2trees
	    float SFweightErr = 0;
	    float SFweight = 1;//outside, since need it there later
	    if((!fMT2tree->misc.isData)){
	    if(!(dofastbtagSFreweighting)){
		vector<float> jetEff;
		vector<float> jetEffErr;
		vector<float> jetSF;
		vector<float> jetSFErr;
		vector<float> jetEffup;
		vector<float> jetEffdown;
		vector<float> jetSFup;
		vector<float> jetSFdown;
		int njetsusuable = 0;
		for(int n = 0; n<fMT2tree->NJets; ++n){
			if(fMT2tree->jet[n].isPFIDLoose==false) continue;
			if(fMT2tree->jet[n].Flavour<=-7777)     continue;
			float effPt  = fMT2tree->jet[n].lv.Pt();
			float effEta = fabs(fMT2tree->jet[n].lv.Eta());
			if(effPt<20.)  continue;
			if(effEta>2.4) continue;
			++njetsusuable;
			if(abs(fMT2tree->jet[n].Flavour)==5){
				jetEff.push_back(    float(hbeff->GetBinContent(hbeff->FindBin(effPt))) );
				jetEffErr.push_back( float(hbeff->GetBinError(  hbeff->FindBin(effPt))) );
				float SFErr;
				float SF = getBTagSF(SFErr, tagger, effPt, effEta);
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr );//first implementation
			}
			else if(abs(fMT2tree->jet[n].Flavour)==4){
				jetEff.push_back(    float(hceff->GetBinContent(hceff->FindBin(effPt))) );
				jetEffErr.push_back( float(hceff->GetBinError(  hceff->FindBin(effPt))) );
				float SFErr;
				float SF = getBTagSF(SFErr, tagger, effPt, effEta );
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr*2.);
			}
			else {
				jetEff.push_back(    float(hleff->GetBinContent(hleff->FindBin(effPt))) );
				jetEffErr.push_back( float(hleff->GetBinError(  hleff->FindBin(effPt))) );
				float SFErr;
				float SF = getMistagSF(SFErr, tagger, effPt, effEta, 0 );
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr );
			}
		}
		SFweight      = getBTagEventWeightError(SFweightErr, jetEff, jetEffErr, jetSF,     jetSFErr, ntags);
		if(SFweight==0){
			if(njetsusuable!=0){
				cout << "Event has zero weight, do not use it" << endl;
				continue;
			}
			else { //event has no flavour information, use average event weight
				SFweight = pow(0.95,abs(ntags));
				if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
				else if(fabs(ntags)>=2) SFweightErr = sqrt((0.0257166*0.0257166 + 0.0370919+0.0370919)*((abs(ntags)-1)*SFweight)*abs(ntags));
				else                    SFweightErr = sqrt( 0.0257166*0.0257166 + 0.0370919+0.0370919);
			}
		}
		weight  = weight * SFweight;
	   }
           else{
	     //from ttbar payload
             if(tagger=="SSVHPT"){
		SFweight = pow(0.95,abs(ntags));
		if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
		weight = weight * SFweight;
	     }
             if(tagger=="SSVHEM"){
		SFweight = pow(0.96,abs(ntags));
		if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
		weight = weight * SFweight;
	     }
             if(tagger=="CSVM"){
		SFweight = pow(0.972,abs(ntags));
		if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
		weight = weight * SFweight;
	     }
	  }
          }
 
	//help vectors, lepton1/2 (including flavour and charge), met, bjets (and there flavour,discriminant), (b)jet indices
	   TLorentzVector l1, l2, met;
	   vector<TLorentzVector> bjets; bjets.clear();
	   vector<int> bjetsflavour; bjetsflavour.clear();
	   vector<float> bjetsdiscriminant; bjetsdiscriminant.clear();
	   int l1c(-999), l2c(-999);
	   Bool_t l1e(true), l2e(true);
	   vector<int> jind; jind.clear();
	   vector<int> bind; bind.clear();
	   if(recoedee){ 
		l1  = (fMT2tree->ele[0]).lv;     l2  = (fMT2tree->ele[1]).lv;
		l1c = (fMT2tree->ele[0]).Charge; l2c = fMT2tree->ele[1].Charge;
		l1e = true;                      l2e = true; }
	   else if(recoedemu){//NOTE: Consider changing this to ele first, muon second (or other way round) --> would get rid of l1e,l2e
		if(fMT2tree->ele[0].lv.Pt()>fMT2tree->muo[0].lv.Pt() ) { 
			l1  = (fMT2tree->ele[0]).lv;     l2  = (fMT2tree->muo[0]).lv; 
			l1c = (fMT2tree->ele[0]).Charge; l2c = fMT2tree->muo[0].Charge;
			l1e = true;                      l2e = false;}
		else  { 
			l1  = (fMT2tree->muo[0]).lv;   l2  = (fMT2tree->ele[0]).lv; 
			l1c = fMT2tree->muo[0].Charge; l2c = fMT2tree->ele[0].Charge;
			l1e = false;                   l2e = true;}
	   }
	   else if(recoedmumu) { 
		l1  = (fMT2tree->muo[0]).lv;   l2  = (fMT2tree->muo[1]).lv; 
		l1c = fMT2tree->muo[0].Charge; l2c = fMT2tree->muo[1].Charge;
		l1e = false;                   l2e = false;}
	   //bjets vector
	   for(int n = 0; n<fMT2tree->NJets; ++n){//loop only over bjets
		if(fMT2tree->jet[n].isPFIDLoose==false)  continue;
		if(fMT2tree->jet[n].lv.Pt()<bpt)         continue;
		if(fabs(fMT2tree->jet[n].lv.Eta())>beta) continue;
		jind.push_back(n);//NOTE: add here an additional pt constraint?????
		float btempdiscr = (Tagger==3 ? fMT2tree->jet[n].bTagProbSSVHP : Tagger==2 ? fMT2tree->jet[n].bTagProbSSVHE : Tagger==4 ? fMT2tree->jet[n].bTagProbCSV :fMT2tree->jet[n].bTagProbCSV);//default CSV
		if(btempdiscr<discr)                     continue;
		bjets.push_back(fMT2tree->jet[n].lv); bjetsflavour.push_back(fMT2tree->jet[n].Flavour); bjetsdiscriminant.push_back(btempdiscr); bind.push_back(n);
	   }
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
	  for(size_t n =0; n<bjets.size()-1;++n) {if(bjetsdiscriminant[n+1]>bjetsdiscriminant[n]) cout << __LINE__ << "discr(ind) " << bjetsdiscriminant[n] << "("<< n<<") " << bjetsdiscriminant[n+1] << "("<<n+1<<")"<<endl;
	  }

	   met = fMT2tree->pfmet[0];
	   if(ftype1met)    met = fMT2tree->type1pfmet[0];
	   else if(frawmet) met = fMT2tree->rawpfmet[0];

	//help vectors for gen-level information
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
				if(gfl1<-50){      genl1 = genlv[n]; gfl1 = genflav[n]; gil1 = genindex[n]; }
				else if(gfl2<-50){ genl2 = genlv[n]; gfl2 = genflav[n]; gil2 = genindex[n]; }
				//else { cout << "this shouldn't happen: flav1 " << gfl1 << " flav2 " << gfl2 << " this flav " << genflav[n] << endl; }
			}
		}
		for(int n =0; n<4; ++n){
			if(gfl1<-50 || gfl2<-50) continue;
			if(abs(genflav[n])==5){
				//if(gfb1>-50 && gfb2>-50) { cout << "this shouldn't happen: bflav1 " << gfb1 << " bflav2 " << gfb2 << " this bflav " << genflav[n] << endl; }
				if(genflav[n]<0){ if(gfl1>0){ //is bbar, has to go 
						 	 gfb1 = genflav[n]; genb1 = genlv[n]; gib1 = genindex[n]; }
						  else if(gfl2>0){ 
							 gfb2 = genflav[n]; genb2 = genlv[n]; gib2 = genindex[n]; } 
				}
				else if(genflav[n]>0){ if(gfl1<0){ 
								 gfb1 = genflav[n]; genb1 = genlv[n]; gib1 = genindex[n]; } 
						  	else if(gfl2<0){ 
								 gfb2 = genflav[n]; genb2 = genlv[n]; gib2 = genindex[n]; } 
				}
			}
		}
		if(genneutrinoindex.size()==2){
			for(int n = 0; n<2;++n){
				if(gfl1*genneutrinoflavour[n]==-182||gfl1*genneutrinoflavour[n]==-132) { //13*-14/-13*14 12*-11/-12*11
					gin1=genneutrinoindex[n]; }
				if(gfl2*genneutrinoflavour[n]==-182||gfl2*genneutrinoflavour[n]==-132) {//13*-14/-13*14 12*-11/-12*11
					gin2=genneutrinoindex[n]; }
			}
		}
	}

	bool truettbarEE(false), truettbarEMu(false), truettbarMuMu(false);
	if(gfl1>-50 && gfl2>-50 && gfb1 >-50 && gfb2>-50){
		//if(gfl1*gfl2>0 || gfb1*gfb2>0) cout << "why flavl1 " << gfl1 << " flavl2 " << gfl2 << " gfb1 " << gfb1 << " gfb2 " << gfb2 << endl;
		if((abs(gfl1)==11 && gfl1*gfb1<0) && (gfl1*gfl2==-121 && gfl2*gfb2<0))                   truettbarEE   = true;
		if((abs(gfl1)==13 && gfl1*gfb1<0) && (gfl1*gfl2==-169 && gfl2*gfb2<0))                   truettbarMuMu = true;
		if( ((abs(gfl1)==11||abs(gfl1)==13) && gfl1*gfb1<0) && (gfl2*gfl1==-143 && gfl2*gfb2<0)) truettbarEMu  = true;
	}

	float genMl1b1(-99.), genMl2b2(-99.), genMT2(-99.), genM1(-99.), genM2(-99.);
	float genMl1b1_allAcceptance(-99.), genMl2b2_allAcceptance(-99.), genMT2_allAcceptance(-99.), genM1_allAcceptance(-99.), genM2_allAcceptance(-99.);
	if(truettbarEE || truettbarEMu || truettbarMuMu){
	   genMl1b1_allAcceptance      = (genl1+genb1).M();
	   genMl2b2_allAcceptance      = (genl2+genb2).M();
	   genMT2_allAcceptance        = fMT2tree->CalcMT2(0., true,  genl1+genb1, genl2+genb2, fMT2tree->genmet[0]);
	   if(gin2>-50 && gin1>-50){
		genM1_allAcceptance    = (genl1+genb1+fMT2tree->genlept[gin1].lv).M();
		genM2_allAcceptance    = (genl2+genb2+fMT2tree->genlept[gin2].lv).M();
	   }
	   if( ((genl1.Pt()>20&&genl2.Pt()>10)||(genl2.Pt()>20&&genl1.Pt()>10)) && genb1.Pt()>bpt && genb2.Pt()>bpt){
	   	genMl1b1      = (genl1+genb1).M();
	   	genMl2b2      = (genl2+genb2).M();
	   	genMT2        = fMT2tree->CalcMT2(0., true,  genl1+genb1, genl2+genb2, fMT2tree->genmet[0]);
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
		if(fMT2tree->jet[n].isPFIDLoose==false)  continue;
		float jetpt = fMT2tree->jet[n].lv.Pt();
		++njetsID20;
		float btempdiscr = (Tagger==3 ? fMT2tree->jet[n].bTagProbSSVHP : Tagger==2 ? fMT2tree->jet[n].bTagProbSSVHE : Tagger==4 ? fMT2tree->jet[n].bTagProbCSV :fMT2tree->jet[n].bTagProbCSV);//default CSV
		if(btempdiscr>discr          ) ++nbjets20;
		if(btempdiscr>discr && jetpt>30.) ++nbjets30;
		if(btempdiscr>discr && jetpt>40.) ++nbjets40;
		if(btempdiscr>discr && jetpt>50.) ++nbjets50;
		if(jetpt>30.) ++njetsID30;
		if(jetpt>40.) ++njetsID40;
		if(jetpt>50.) ++njetsID50;
	   }

		//these are all the needed high level objects (except D2)
	   float MT2lbV(-999.99), MT2lbV_massless(-999.99);
	   float MT2lbV_withMlb(-999.99), MT2lbV_massless_withMlb(-999.99);//as this might be difficult to get
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
	   int ind_l1bxMT2lb(-1), ind_l2bxMT2lb(-1), ind_l1bxMT2lbMlb(-1), ind_l2bxMT2lbMlb(-1);
	   int ind_b1x_MT2bbmin(-1), ind_b2x_MT2bbmin(-1), ind_b1x_MT2bbmin_mW_linMET(-1), ind_b2x_MT2bbmin_mW_linMET(-1), ind_b1x_MT2bbmin_linMET(-1), ind_b2x_MT2bbmin_linMET(-1), ind_l1bxMT2lbMassless(-1), ind_l2bxMT2lbMassless(-1);
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
			MT2_l1b[i1][i2]  = -999.99; MT2_massless_l1b[i1][i2] = -999.99;
			MT2bbarr[i1][i2] = -999.99; MT2bb_mW_linMET[i1][i2]  = -999.99; MT2bb_linMET[i1][i2] = -999.99; 
			DRbbarr[i1][i2]  = -999.99; DPhibbarr[i1][i2]        = -999.99; UTMvec[i1][i2]       = -999.99; Mbbarr[i1][i2] = -999.;
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
		MT2ll  = fMT2tree->CalcMT2(0., true,  l1, l2, met);
		Mll    = (l1+l2).M();
		DPhill = fabs(l1.DeltaPhi(l2));
		DRll   = l1.DeltaR(l2);
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
			DRbbarr[i1][i2]              = (bjets[i1]).DeltaR(bjets[i2]);
			DPhibbarr[i1][i2]            = fabs((bjets[i1]).DeltaPhi(bjets[i2]));
			MT2_l1b[i1][i2]              = fMT2tree->CalcMT2(0.,   true,  l1+bjets[i1], l2+bjets[i2], met);
	   		MT2_massless_l1b[i1][i2]     = fMT2tree->CalcMT2(0.,  false,  l1+bjets[i1], l2+bjets[i2], met);
			MT2bbarr[i1][i2]             = fMT2tree->CalcMT2(0.,   true,  bjets[i1],    bjets[i2],    met);
			MT2bb_mW_linMET[i1][i2]      = fMT2tree->CalcMT2(80.4, true,  bjets[i1],    bjets[i2],    met+l1+l2);
			MT2bb_linMET[i1][i2]         = fMT2tree->CalcMT2(0.,   true,  bjets[i1],    bjets[i2],    met+l1+l2);
			UTMvec[i1][i2]               = sqrt(pow((l1+l2+bjets[i1]+bjets[i2]+met).Px(),2) + pow((l1+l2+bjets[i1]+bjets[i2]+met).Py(),2));
			Mbbarr[i1][i2]               = (bjets[i1]+bjets[i2]).M();

		}
	   }
	   //get minimal values
	   for(unsigned int i1 = 0; i1<NumBJets; ++i1){
		for(unsigned int i2 = 0; i2<NumBJets; ++i2){
			if(i2==i1) continue;
			if(MT2lbV<0. && MT2_l1b[i1][i2]>=0.) { MT2lbV = MT2_l1b[i1][i2];ind_l1bxMT2lb = i1; ind_l2bxMT2lb = i2; }
			else if (MT2lbV>=0. && MT2_l1b[i1][i2]>=0. && MT2_l1b[i1][i2]<MT2lbV)   {  
				MT2lbV = MT2_l1b[i1][i2]; ind_l1bxMT2lb = i1; ind_l2bxMT2lb = i2; }
			if(MT2lbV_massless<0. && MT2_massless_l1b[i1][i2]>=0.) { MT2lbV_massless = MT2_massless_l1b[i1][i2]; ind_l1bxMT2lbMassless = i1; ind_l2bxMT2lbMassless = i2; }
			else if (MT2lbV_massless>=0. && MT2_massless_l1b[i1][i2]>=0. && MT2_massless_l1b[i1][i2]<MT2lbV) {   
				 MT2lbV_massless = MT2_massless_l1b[i1][i2]; ind_l1bxMT2lbMassless = i1; ind_l2bxMT2lbMassless = i2; }
			if(MT2lbV_withMlb<0. && MT2_l1b[i1][i2]>=0. && Ml1b[i1]<180. && Ml2b[i2]<180. ) {
				 MT2lbV_withMlb = MT2_l1b[i1][i2]; ind_l1bxMT2lbMlb = i1; ind_l2bxMT2lbMlb = i2; }
			else if(MT2lbV_withMlb>=0. && MT2_l1b[i1][i2]>=0. && Ml1b[i1]<180. && Ml2b[i2]<180. && MT2lbV_withMlb>MT2_l1b[i1][i2]) {
				 MT2lbV_withMlb = MT2_l1b[i1][i2]; ind_l1bxMT2lbMlb = i1; ind_l2bxMT2lbMlb = i2; }
			if(MT2lbV_massless_withMlb<0. && MT2_massless_l1b[i1][i2]>=0. && Ml1b[i1]<180. && Ml2b[i2]<180. ) 
				MT2lbV_massless_withMlb = MT2_massless_l1b[i1][i2];
			else if(MT2lbV_massless_withMlb>=0. && MT2_massless_l1b[i1][i2]>=0. && Ml1b[i1]<180. && Ml2b[i2]<180. && MT2lbV_massless_withMlb>MT2_massless_l1b[i1][i2]) MT2lbV_massless_withMlb = MT2_massless_l1b[i1][i2];

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

	//new checks 14/06/2012 w.r.t. to gen-level plots - this should be deleted
	//first gen reco matching if possible -- using produced vectors genlv, genflav, genindex
	//these vectors contain all b's and leptons coming from top(!!!!), via flav can also match to same top//no tau decays in here
	//match indices via dR only, no pT
	int matchedgenl_tol1_regardlessFlavourAndCharge(-1), matchedgenl_tol2_regardlessFlavourAndCharge(-1);
	int matchedgenl_tol1_regardlessCharge(-1), matchedgenl_tol2_regardlessCharge(-1);
	int matchedgenl_tol1_regardlessFlavour(-1), matchedgenl_tol2_regardlessFlavour(-1);//still require e or mu
	int matchedgenl_tol1(-1), matchedgenl_tol2(-1);
	int matchedgenp_tob1(-1), matchedgenp_tob2(-1);//b1 = bind[ind_l1bxMT2lb], b2 = bind[ind_l2bxMT2lb];
	float dR_ml1fc(99.), dR_ml2fc(99.), dR_ml1c(99.), dR_ml2c(99.), dR_ml1f(99.), dR_ml2f(99.), dR_ml1(99.), dR_ml2(99.), dR_mb1(99.), dR_mb2(99.);
	//start l,b true top
	if(sampletype !="Stop" || stopTotopLSP){
		for(int n = 0; n<25; ++n){
			if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
			int genlepID   = fMT2tree->genlept[n].ID;
			int genlepMID  = fMT2tree->genlept[n].MID;
			int genlepGMID = fMT2tree->genlept[n].GMID;
			if(((abs(genlepID)==11)||(abs(genlepID)==13))){
				if(abs(genlepMID)!=24 && abs(genlepMID)!=15 )     continue;
				if(abs(genlepGMID)!=6 && abs(genlepMID)==24 )     continue;
				if(abs(genlepGMID)!=24&& abs(genlepMID)==15 )     continue;
				if(abs(genlepMID)==24 && genlepMID *genlepID>=0 ) continue;
				if(abs(genlepMID)==15 && genlepMID *genlepID<=0 ) continue;
				if(genlepGMID*genlepID>=0)                        continue;
				float dRl1 = fMT2tree->genlept[n].lv.DeltaR(l1); 
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1fc) continue;
					matchedgenl_tol1_regardlessFlavourAndCharge = n; dR_ml1fc = dRl1; }
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1c)                 continue;
					if(l1e==true &&abs(genlepID)!=11) continue;
					if(l1e==false&&abs(genlepID)!=13) continue;
					matchedgenl_tol1_regardlessCharge = n; dR_ml1c= dRl1; }
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1f)  continue;
					if(l1c*genlepID>0) continue;
					matchedgenl_tol1_regardlessFlavour = n; dR_ml1f= dRl1; }
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1)                  continue;
					if(l1e==true &&abs(genlepID)!=11) continue;
					if(l1e==false&&abs(genlepID)!=13) continue;
					if(l1c*genlepID>0)                continue;
					matchedgenl_tol1 = n; dR_ml1 = dRl1; }//11*-1 or -11*1 --> true
			}
		}
		for(int n = 0; n<25; ++n){
			if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
			int genlepID   = fMT2tree->genlept[n].ID;
			int genlepMID  = fMT2tree->genlept[n].MID;
			int genlepGMID = fMT2tree->genlept[n].GMID;
			if(((abs(genlepID)==11)||(abs(genlepID)==13))){
				//include tau decays with good hopes
				if(abs(genlepMID)!=24 && abs(genlepMID)!=15 )     continue;
				if(abs(genlepGMID)!=6 && abs(genlepMID)==24 )     continue;
				if(abs(genlepGMID)!=24&& abs(genlepMID)==15 )     continue;
				if(abs(genlepMID)==24 && genlepMID *genlepID>=0 ) continue;
				if(abs(genlepMID)==15 && genlepMID *genlepID<=0 ) continue;
				if(genlepGMID*genlepID>=0)                        continue;
				float dRl2 = fMT2tree->genlept[n].lv.DeltaR(l2);
				if(dRl2<0.3 && dRl2<dR_ml2fc) {
					if(dRl2>=dR_ml2fc) continue;
					matchedgenl_tol2_regardlessFlavourAndCharge = n; dR_ml2fc = dRl2; }
				if(dRl2<0.3 && dRl2<dR_ml2c) {
					if(dRl2>=dR_ml2c)                                       continue;
					if(l2e==true &&abs(genlepID)!=11 && abs(genlepMID)!=15) continue;
					if(l2e==false&&abs(genlepID)!=13 && abs(genlepMID)!=15) continue;
					matchedgenl_tol2_regardlessCharge = n; dR_ml2c = dRl2; }
				if(dRl2<0.3) {
					if(dRl2>=dR_ml2f)  continue;
					if(l2c*genlepID>0) continue;
					matchedgenl_tol2_regardlessFlavour = n; dR_ml2f= dRl2; }
				if(dRl2<0.3 && dRl2<dR_ml2) {
					if(dRl2>=dR_ml2)                                        continue;
					if(l2e==true &&abs(genlepID)!=11 && abs(genlepMID)!=15) continue;
					if(l2e==false&&abs(genlepID)!=13 && abs(genlepMID)!=15) continue;
					if(l2c*genlepID>0)                                      continue;
					matchedgenl_tol2 = n; dR_ml2 = dRl2; }//11*-1 or -11*1 --> true
			}
		}
		for(int n = 0; n<25; ++n){
			if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
			int genlepID  = fMT2tree->genlept[n].ID;
			int genlepMID = fMT2tree->genlept[n].MID;
			//int genlepGMID = fMT2tree->genlept[n].GMID;
			if(abs(genlepID)==5){
				if(abs(genlepMID)!=6)     continue;
				if(genlepMID*genlepID<=0) continue;
				float dRb1 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l1bxMT2lb]);
				if(dRb1<0.5 && dRb1<dR_mb1){ matchedgenp_tob1 = n; dR_mb1 = dRb1; }
				float dRb2 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l2bxMT2lb]);
				if(dRb2<0.5 && dRb2<dR_mb2){ matchedgenp_tob2 = n; dR_mb2 = dRb1; }
			}
		}
	//end true b,l top
	} else if(sampletype =="Stop"){
	   if(stopTobCharginoWLSP){
		for(int n = 0; n<25; ++n){
			if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
			int genlepID   = fMT2tree->genlept[n].ID;
			int genlepMID  = fMT2tree->genlept[n].MID;
			int genlepGMID = fMT2tree->genlept[n].GMID;
			if(((abs(genlepID)==11)||(abs(genlepID)==13))){
				if(abs(genlepMID)!=24       && abs(genlepMID)!=15 )     continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepMID)==24 )     continue;
				if(abs(genlepGMID)!=24      && abs(genlepMID)==15 )     continue;
				if(abs(genlepMID)==24       && genlepMID *genlepID>=0 ) continue;
				if(abs(genlepMID)==15       && genlepMID *genlepID<=0 ) continue;
				if(genlepGMID*genlepID>=0)                              continue;
				float dRl1 = fMT2tree->genlept[n].lv.DeltaR(l1); 
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1fc) continue;
					matchedgenl_tol1_regardlessFlavourAndCharge = n; dR_ml1fc = dRl1; }
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1c)                 continue;
					if(l1e==true &&abs(genlepID)!=11) continue;
					if(l1e==false&&abs(genlepID)!=13) continue;
					matchedgenl_tol1_regardlessCharge = n; dR_ml1c= dRl1; }
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1f)  continue;
					if(l1c*genlepID>0) continue;
					matchedgenl_tol1_regardlessFlavour = n; dR_ml1f= dRl1; }
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1)                  continue;
					if(l1e==true &&abs(genlepID)!=11) continue;
					if(l1e==false&&abs(genlepID)!=13) continue;
					if(l1c*genlepID>0)                continue;
					matchedgenl_tol1 = n; dR_ml1 = dRl1; }//11*-1 or -11*1 --> true
			}
		}
		for(int n = 0; n<25; ++n){
			if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
			int genlepID   = fMT2tree->genlept[n].ID;
			int genlepMID  = fMT2tree->genlept[n].MID;
			int genlepGMID = fMT2tree->genlept[n].GMID;
			if(((abs(genlepID)==11)||(abs(genlepID)==13))){
				//include tau decays with good hopes
				if(abs(genlepMID)!=24       && abs(genlepMID)!=15 )     continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepMID)==24 )     continue;
				if(abs(genlepGMID)!=24      && abs(genlepMID)==15 )     continue;
				if(abs(genlepMID)==24       && genlepMID *genlepID>=0 ) continue;
				if(abs(genlepMID)==15       && genlepMID *genlepID<=0 ) continue;
				if(genlepGMID*genlepID>=0)                              continue;
				float dRl2 = fMT2tree->genlept[n].lv.DeltaR(l2);
				if(dRl2<0.3 && dRl2<dR_ml2fc) {
					if(dRl2>=dR_ml2fc) continue;
					matchedgenl_tol2_regardlessFlavourAndCharge = n; dR_ml2fc = dRl2; }
				if(dRl2<0.3 && dRl2<dR_ml2c) {
					if(dRl2>=dR_ml2c)                                       continue;
					if(l2e==true &&abs(genlepID)!=11 && abs(genlepMID)!=15) continue;
					if(l2e==false&&abs(genlepID)!=13 && abs(genlepMID)!=15) continue;
					matchedgenl_tol2_regardlessCharge = n; dR_ml2c = dRl2; }
				if(dRl2<0.3) {
					if(dRl2>=dR_ml2f)  continue;
					if(l2c*genlepID>0) continue;
					matchedgenl_tol2_regardlessFlavour = n; dR_ml2f= dRl2; }
				if(dRl2<0.3 && dRl2<dR_ml2) {
					if(dRl2>=dR_ml2)                                        continue;
					if(l2e==true &&abs(genlepID)!=11 && abs(genlepMID)!=15) continue;
					if(l2e==false&&abs(genlepID)!=13 && abs(genlepMID)!=15) continue;
					if(l2c*genlepID>0)                                      continue;
					matchedgenl_tol2 = n; dR_ml2 = dRl2; }//11*-1 or -11*1 --> true
			}
		}
		for(int n = 0; n<25; ++n){
			if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
			int genlepID  = fMT2tree->genlept[n].ID;
			int genlepMID = fMT2tree->genlept[n].MID;
			//int genlepGMID = fMT2tree->genlept[n].GMID;
			if(abs(genlepID)==5){
				if(abs(genlepMID)==2000006) cout << "Cool generated t_2" << endl;;
				if(abs(genlepMID)!=1000006) continue;
				if(genlepMID*genlepID<=0)   continue;
				float dRb1 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l1bxMT2lb]);
				if(dRb1<0.5 && dRb1<dR_mb1){
					matchedgenp_tob1 = n; dR_mb1 = dRb1; }
				float dRb2 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l2bxMT2lb]);
				if(dRb2<0.5 && dRb2<dR_mb2){
					matchedgenp_tob2 = n; dR_mb2 = dRb1; }
			}
		}
	   }//stopTobCharginoWLSP
	   else if(stopTobCharginoslepsnu){
		for(int n = 0; n<25; ++n){
			if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
			int genlepID   = fMT2tree->genlept[n].ID;
			int genlepMID  = fMT2tree->genlept[n].MID;
			int genlepGMID = fMT2tree->genlept[n].GMID;
			if(((abs(genlepID)==11)||(abs(genlepID)==13))){
				if(abs(genlepMID )!=1000024 && abs(genlepMID)!=15 && 
                                   abs(genlepMID )!=1000011 && abs(genlepMID)!=1000013 )                        continue;
				if(abs(genlepGMID)!=1000006 && abs(genlepMID)==1000024 )                        continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepMID)==1000013 )                        continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepMID)==1000011 )                        continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepGMID)!=1000015 && abs(genlepMID)==15 ) continue;
				if(abs(genlepMID )==1000024 && genlepMID *genlepID>=0  )                        continue;
				if(abs(genlepMID )==1000013 && genlepMID *genlepID<=0  )                        continue;
				if(abs(genlepMID )==1000011 && genlepMID *genlepID<=0  )                        continue;
				if(abs(genlepMID )==15      && genlepMID *genlepID<=0  )                        continue;
				if(genlepGMID*genlepID>=0) continue;
				float dRl1 = fMT2tree->genlept[n].lv.DeltaR(l1); 
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1fc) continue;
					matchedgenl_tol1_regardlessFlavourAndCharge = n; dR_ml1fc = dRl1; }
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1c)                 continue;
					if(l1e==true &&abs(genlepID)!=11) continue;
					if(l1e==false&&abs(genlepID)!=13) continue;
					matchedgenl_tol1_regardlessCharge = n; dR_ml1c= dRl1; }
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1f)  continue;
					if(l1c*genlepID>0) continue;
					matchedgenl_tol1_regardlessFlavour = n; dR_ml1f= dRl1; }
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1)                  continue;
					if(l1e==true &&abs(genlepID)!=11) continue;
					if(l1e==false&&abs(genlepID)!=13) continue;
					if(l1c*genlepID>0)                continue;
					matchedgenl_tol1 = n; dR_ml1 = dRl1; }//11*-1 or -11*1 --> true
			}
		}
		for(int n = 0; n<25; ++n){
			if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
			int genlepID   = fMT2tree->genlept[n].ID;
			int genlepMID  = fMT2tree->genlept[n].MID;
			int genlepGMID = fMT2tree->genlept[n].GMID;
			if(((abs(genlepID)==11)||(abs(genlepID)==13))){
				//include tau decays with good hopes
				if(abs(genlepMID )!=1000024 && abs(genlepMID)!=15 && 
				   abs(genlepMID )!=1000011 && abs(genlepMID)!=1000013 )                        continue;
				if(abs(genlepGMID)!=1000006 && abs(genlepMID)==1000024 )                        continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepMID)==1000013 )                        continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepMID)==1000011 )                        continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepGMID)!=1000015 && abs(genlepMID)==15 ) continue;
				if(abs(genlepMID )==1000024 && genlepMID *genlepID>=0  )                        continue;
				if(abs(genlepMID )==1000013 && genlepMID *genlepID<=0  )                        continue;
				if(abs(genlepMID )==1000011 && genlepMID *genlepID<=0  )                        continue;
				if(abs(genlepMID )==15      && genlepMID *genlepID<=0  )                        continue;
				if(genlepGMID*genlepID>=0) continue;
				float dRl2 = fMT2tree->genlept[n].lv.DeltaR(l2);
				if(dRl2<0.3 && dRl2<dR_ml2fc) {
					if(dRl2>=dR_ml2fc) continue;
					matchedgenl_tol2_regardlessFlavourAndCharge = n; dR_ml2fc = dRl2; }
				if(dRl2<0.3 && dRl2<dR_ml2c) {
					if(dRl2>=dR_ml2c)                                       continue;
					if(l2e==true &&abs(genlepID)!=11 && abs(genlepMID)!=15) continue;
					if(l2e==false&&abs(genlepID)!=13 && abs(genlepMID)!=15) continue;
					matchedgenl_tol2_regardlessCharge = n; dR_ml2c = dRl2; }
				if(dRl2<0.3) {
					if(dRl2>=dR_ml2f)  continue;
					if(l2c*genlepID>0) continue;
					matchedgenl_tol2_regardlessFlavour = n; dR_ml2f= dRl2; }
				if(dRl2<0.3 && dRl2<dR_ml2) {
					if(dRl2>=dR_ml2)                                        continue;
					if(l2e==true &&abs(genlepID)!=11 && abs(genlepMID)!=15) continue;
					if(l2e==false&&abs(genlepID)!=13 && abs(genlepMID)!=15) continue;
					if(l2c*genlepID>0)                                      continue;
					matchedgenl_tol2 = n; dR_ml2 = dRl2; }//11*-1 or -11*1 --> true
			}
		}
		for(int n = 0; n<25; ++n){
			if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
			int genlepID  = fMT2tree->genlept[n].ID;
			int genlepMID = fMT2tree->genlept[n].MID;
			//int genlepGMID = fMT2tree->genlept[n].GMID;
			if(abs(genlepID)==5){
				if(abs(genlepMID)==2000006) cout << "Cool generated t_2" << endl;;
				if(abs(genlepMID)!=1000006) continue;
				if(genlepMID*genlepID<=0)   continue;
				float dRb1 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l1bxMT2lb]);
				if(dRb1<0.5 && dRb1<dR_mb1){
					matchedgenp_tob1 = n; dR_mb1 = dRb1; }
				float dRb2 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l2bxMT2lb]);
				if(dRb2<0.5 && dRb2<dR_mb2){
					matchedgenp_tob2 = n; dR_mb2 = dRb1; }
			}
		}
	   }//stopTobCharginoslepsnu
	   else if(stopTobChargino){//so generic, that signal above might be killed
		for(int n = 0; n<25; ++n){
			if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
			int genlepID   = fMT2tree->genlept[n].ID;
			int genlepMID  = fMT2tree->genlept[n].MID;
			int genlepGMID = fMT2tree->genlept[n].GMID;
			if(((abs(genlepID)==11)||(abs(genlepID)==13))){
				if(abs(genlepGMID)!=1000006 &&abs(genlepGMID)!=1000024 &&abs(genlepGMID)!=1000015) continue;
				if(abs(genlepGMID)==1000006 &&abs(genlepMID) !=1000024 )                           continue;
				if(abs(genlepGMID)==1000015 &&abs(genlepMID) !=15      )                           continue;
				if(genlepGMID*genlepID>=0)                                                         continue;
				float dRl1 = fMT2tree->genlept[n].lv.DeltaR(l1); 
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1fc) continue;
					matchedgenl_tol1_regardlessFlavourAndCharge = n; dR_ml1fc = dRl1; }
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1c)                 continue;
					if(l1e==true &&abs(genlepID)!=11) continue;
					if(l1e==false&&abs(genlepID)!=13) continue;
					matchedgenl_tol1_regardlessCharge = n; dR_ml1c= dRl1; }
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1f)  continue;
					if(l1c*genlepID>0) continue;
					matchedgenl_tol1_regardlessFlavour = n; dR_ml1f= dRl1; }
				if(dRl1<0.3) {
					if(dRl1>=dR_ml1)                  continue;
					if(l1e==true &&abs(genlepID)!=11) continue;
					if(l1e==false&&abs(genlepID)!=13) continue;
					if(l1c*genlepID>0)                continue;
					matchedgenl_tol1 = n; dR_ml1 = dRl1; }//11*-1 or -11*1 --> true
			}
		}
		for(int n = 0; n<25; ++n){
			if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
			int genlepID   = fMT2tree->genlept[n].ID;
			int genlepMID  = fMT2tree->genlept[n].MID;
			int genlepGMID = fMT2tree->genlept[n].GMID;
			if(((abs(genlepID)==11)||(abs(genlepID)==13))){
				//include tau decays with good hopes
				if(abs(genlepGMID)!=1000006 &&abs(genlepGMID)!=1000024 &&abs(genlepGMID)!=1000015) continue;
				if(abs(genlepGMID)==1000006 &&abs(genlepMID) !=1000024) continue;
				if(abs(genlepGMID)==1000015 &&abs(genlepMID) !=15     ) continue;
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
				if(dRl2<0.3) {
					if(dRl2>=dR_ml2f) continue;
					if(l2c*genlepID>0) continue;
					matchedgenl_tol2_regardlessFlavour = n; dR_ml2f= dRl2; }
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
			int genlepID  = fMT2tree->genlept[n].ID;
			int genlepMID = fMT2tree->genlept[n].MID;
			//int genlepGMID = fMT2tree->genlept[n].GMID;
			if(abs(genlepID)==5){
				if(abs(genlepMID)==2000006) cout << "Cool generated t_2" << endl;;
				if(abs(genlepMID)!=1000006) continue;
				if(genlepMID*genlepID<=0)   continue;
				float dRb1 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l1bxMT2lb]);
				if(dRb1<0.5 && dRb1<dR_mb1){
					matchedgenp_tob1 = n; dR_mb1 = dRb1; }
				float dRb2 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l2bxMT2lb]);
				if(dRb2<0.5 && dRb2<dR_mb2){
					matchedgenp_tob2 = n; dR_mb2 = dRb1; }
			}
		}
	   }//stopTobChargino, too generic, might not pick up everything in order not to keep fakes
	}

	int gl1(-1), gl2(-1), gb1(-1), gb2(-1);
	vector<int> gvv; gvv.clear();
	   for(int n = 0; n<25; ++n){
		bool add = true;
		if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
		int genlepID   = fMT2tree->genlept[n].ID;
		int genlepMID  = fMT2tree->genlept[n].MID;
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
		int genlepID  = fMT2tree->genlept[n].ID;
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
				if(gl1>=0&&gl1!=gvv[n]){ if(genlepID*fMT2tree->genlept[gl1 ].ID < 0) gl2 = gvv[n]; }
			}
		}
		for(unsigned int n = 0; n<gvv.size(); ++n){
			int genlepID = fMT2tree->genlept[gvv[n] ].ID;
			if((abs(genlepID)==5)&&gb1==-1) gb1 = gvv[n];
			if((abs(genlepID)==5)&&gb2==-1){
				if(gb1>=0&&gb1!=gvv[n]){ if(genlepID*fMT2tree->genlept[gb1 ].ID < 0) gb2 = gvv[n]; }
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

	//05/07/2012 new matching - matching to b, l(e,mu) -> store matched genlep ID,MID,GMID, matched yes or no - 3rd gen-level selection
	int b1_matchedID(-99), b1_matchedMID(-99), b1_matchedGMID(-99), b1_PF2PATmatch(-99); bool b1_matched(false);//allow only id 1-9,21
	int b2_matchedID(-99), b2_matchedMID(-99), b2_matchedGMID(-99), b2_PF2PATmatch(-99); bool b2_matched(false);//allow only id 1-9,21
	int l1_matchedID(-99), l1_matchedMID(-99), l1_matchedGMID(-99), l1_flavour(-99);     bool l1_matched(false);//allow only id 11,13
	int l2_matchedID(-99), l2_matchedMID(-99), l2_matchedGMID(-99), l2_flavour(-99);     bool l2_matched(false);//allow only id 11,13
	double dR_match05072012(99.);
	int ind_match05072012(-1);
	for(int n = 0; n<25; ++n){
		if(fMT2tree->genlept[n].lv.Pt()<=0.01                                  ) continue;
		if(fabs(fMT2tree->genlept[n].lv.Eta())>=1000.                          ) continue;
		if(abs(fMT2tree->genlept[n].ID)!=11 && abs(fMT2tree->genlept[n].ID)!=13) continue;
		double dRtemp = fMT2tree->genlept[n].lv.DeltaR(l1);
		if(dRtemp<=0.3 && dR_match05072012>dRtemp){ dR_match05072012 = dRtemp; ind_match05072012 = n; }
	}
	if(l1e){
		if(l1c>0) l1_flavour = -11;
		else      l1_flavour =  11;
	}
	else{
		if(l1c>0) l1_flavour = -13;
		else      l1_flavour =  13;
	}
	if(ind_match05072012>=0){
		l1_matched      = true;
		l1_matchedID    = fMT2tree->genlept[ind_match05072012].ID;
		l1_matchedMID   = fMT2tree->genlept[ind_match05072012].MID;
		l1_matchedGMID  = fMT2tree->genlept[ind_match05072012].GMID;
	}
	dR_match05072012 = 99.; ind_match05072012 = -1;
	for(int n = 0; n<25; ++n){
		if(fMT2tree->genlept[n].lv.Pt()<=0.01                                  ) continue;
		if(fabs(fMT2tree->genlept[n].lv.Eta())>=1000.                          ) continue;
		if(abs(fMT2tree->genlept[n].ID)!=11 && abs(fMT2tree->genlept[n].ID)!=13) continue;
		double dRtemp = fMT2tree->genlept[n].lv.DeltaR(l2);
		if(dRtemp<=0.3 && dR_match05072012>dRtemp){ dR_match05072012 = dRtemp; ind_match05072012 = n; }
	}
	if(l2e){
		if(l2c>0) l2_flavour = -11;
		else      l2_flavour =  11;
	}
	else{
		if(l2c>0) l2_flavour = -13;
		else      l2_flavour =  13;
	}
	if(ind_match05072012>=0){
		l2_matched      = true;
		l2_matchedID    = fMT2tree->genlept[ind_match05072012].ID;
		l2_matchedMID   = fMT2tree->genlept[ind_match05072012].MID;
		l2_matchedGMID  = fMT2tree->genlept[ind_match05072012].GMID;
	}
	dR_match05072012 = 99.; ind_match05072012 = -1;
	for(int n = 0; n<25; ++n){
		if(fMT2tree->genlept[n].lv.Pt()<=0.01                                  ) continue;
		if(fabs(fMT2tree->genlept[n].lv.Eta())>=1000.                          ) continue;
		if(abs(fMT2tree->genlept[n].ID)>9 && abs(fMT2tree->genlept[n].ID)!=21  ) continue;
		if(fMT2tree->genlept[n].ID==0                                          ) continue;
		double dRtemp = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l1bxMT2lb]);
		if(dRtemp<=0.5 && dR_match05072012>dRtemp){ dR_match05072012 = dRtemp; ind_match05072012 = n; }
	}
	b1_PF2PATmatch = fMT2tree->jet[bind[ind_l1bxMT2lb] ].Flavour;
	if(ind_match05072012>=0){
		b1_matched      = true;
		b1_matchedID    = fMT2tree->genlept[ind_match05072012].ID;
		b1_matchedMID   = fMT2tree->genlept[ind_match05072012].MID;
		b1_matchedGMID  = fMT2tree->genlept[ind_match05072012].GMID;
	}
	dR_match05072012 = 99.; ind_match05072012 = -1;
	for(int n = 0; n<25; ++n){
		if(fMT2tree->genlept[n].lv.Pt()<=0.01                                  ) continue;
		if(fabs(fMT2tree->genlept[n].lv.Eta())>=1000.                          ) continue;
		if(abs(fMT2tree->genlept[n].ID)>9 && abs(fMT2tree->genlept[n].ID)!=21  ) continue;
		if(fMT2tree->genlept[n].ID==0                                          ) continue;
		double dRtemp = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l2bxMT2lb]);
		if(dRtemp<=0.5 && dR_match05072012>dRtemp){ dR_match05072012 = dRtemp; ind_match05072012 = n; }
	}
	b2_PF2PATmatch = fMT2tree->jet[bind[ind_l2bxMT2lb] ].Flavour;
	if(ind_match05072012>=0){
		b2_matched      = true;
		b2_matchedID    = fMT2tree->genlept[ind_match05072012].ID;
		b2_matchedMID   = fMT2tree->genlept[ind_match05072012].MID;
		b2_matchedGMID  = fMT2tree->genlept[ind_match05072012].GMID;
	}

	//new D2 implementation
	double b1arr[4], b2arr[4], l1arr[4], l2arr[4], metarr[4];
	b1arr[0] = bjets[ind_l1bxMT2lb].Px(); b1arr[1] = bjets[ind_l1bxMT2lb].Py();
	b1arr[2] = bjets[ind_l1bxMT2lb].Pz(); b1arr[3] = bjets[ind_l1bxMT2lb].P(); 
	b2arr[0] = bjets[ind_l2bxMT2lb].Px(); b2arr[1] = bjets[ind_l2bxMT2lb].Py(); 
	b2arr[2] = bjets[ind_l2bxMT2lb].Pz(); b2arr[3] = bjets[ind_l2bxMT2lb].P(); 
	l1arr[0]  = l1.Px();  l1arr[1]  = l1.Py();  l1arr[2]  = l1.Pz();  l1arr[3]  = l1.P();
	l2arr[0]  = l2.Px();  l2arr[1]  = l2.Py();  l2arr[2]  = l2.Pz();  l2arr[3]  = l2.P();
	metarr[0] = met.Px(); metarr[1] = met.Py(); metarr[2] = met.Pz(); metarr[3] = met.P();
	double discriminant(-999.99), edm(9999.),q1z(0.),q2z(0.); int status(-99);
	double discriminant1(-999.99);//, discriminant2(-999.99);
	int status1(-99);//, status2(-99);
	double edm1(9999.), q1z1(0.),q2z1(0.);
	int statvec1[3]; 
	if(D2calc&&D2calcMinuit){//NOTE: here also Luc's implementation needed, for 4th order equation solving
		sttbLuc->SetNscan(1000.);
		sttbLuc->SetDiscrConv(0.05);
	    	int globalsolution1 = 0;// int globalsolution2 = 0;
		int nsol1 = sttbLuc->SolveEqn4(b1arr, b2arr, l1arr, l2arr, metarr);
		if(nsol1>0) {	discriminant1 = 0.; double qq11[]={0.,0.,0.}, qq22[]={0.,0.,0.}; 
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
			if((fabs(q1z1)>=8000.||fabs(q2z1)>=8000.||q1z1==0.||q2z1==0.||edm1>0.01||status1!=0)) ++vetoeddiscrcount; 
			else ++discrcount;
			globsol += globalsolution1;
		}
	}

	if(D2calc&&D2calcLuc){
	    	statvec1[0]=0; statvec1[1]=0; statvec1[2]=0; int nstptot1 = 0; int ndertot1 = 0; int nrottot1 = 0; int nscntot1 = 0;
	    	//int statvec2[3]; statvec2[0]=0; statvec2[1]=0; statvec2[2]=0; int nstptot2 = 0; int ndertot2 = 0; int nrottot2 = 0; int nscntot2 = 0;
	    	int globalsolution1 = 0;// int globalsolution2 = 0;
		sttbLuc->SetNscan(1000.);
		sttbLuc->SetDiscrConv(0.05);
		int nsol1 = sttbLuc->SolveEqn4(b1arr, b2arr, l1arr, l2arr, metarr);
		if(nsol1>0) {	discriminant1 = 0.; double qq11[]={0.,0.,0.}, qq22[]={0.,0.,0.}; 
				sttbLuc->GetMomSol(0,qq11,qq22); q1z1 = qq11[2]; q2z1 = qq22[2]; ++globalsolution1;
		}
		else {  discriminant1 = sttbLuc->IterDeriv(b1arr, b2arr, l1arr, l2arr, metarr);//luc's code
			q1z1 = sttbLuc->GetXMin(); q2z1 = sttbLuc->GetYMin();
			int stat1 = sttbLuc->GetStatus ();    if (stat1 >= 0) statvec1[stat1]++;
			nstptot1 += sttbLuc->GetNSteps ();    ndertot1 += sttbLuc->GetNDerSteps ();
			nrottot1 += sttbLuc->GetNRotSteps (); nscntot1 += sttbLuc->GetNScanSteps ();
		}
		statvec[0] += statvec1[0]; statvec[1] += statvec1[1]; statvec[2] += statvec1[2];
		nstptot += nstptot1; ndertot += ndertot1; nrottot += nrottot1; nscntot += nscntot1;
		globalsolution += globalsolution1;
	}
	double D2 = discriminant1;//stupid, but for simpleneww use this

	if(debugD2){
		//this is really for debugging purposes
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
			cout << "b1 PxPyPzP " << b1arr[0] << ", " << b1arr[1] << ", " << b1arr[2] << ", " << b1arr[3] << " Discr " << bjetsdiscriminant[ind_l1bxMT2lb] << " index " << ind_l1bxMT2lb << endl;
			cout << "b2 PxPyPzP " << b2arr[0] << ", " << b2arr[1] << ", " << b2arr[2] << ", " << b2arr[3] << " Discr " << bjetsdiscriminant[ind_l2bxMT2lb] << " index " << ind_l2bxMT2lb << endl;
			cout << "l1 PxPyPzP " << l1arr[0] << ", " << l1arr[1] << ", " << l1arr[2] << ", " << l1arr[3] << endl;
			cout << "l2 PxPyPzP " << l2arr[0] << ", " << l2arr[1] << ", " << l2arr[2] << ", " << l2arr[3] << endl;
			cout << "met PxPy   " << metarr[0]<< ", " << metarr[1]<< endl;
			cout <<  endl;
		}
	}


	//Finally - start filling the histograms
	//as you can see several (basically equal) histograms are filled with some small cuts changed to see (really see in plots) their influence
	string hs = string("_") + leptype + string("_") + sampletype;

	//NOTE: might change this depending on jetthreshold
	if(njets>=0 && fMT2tree->NJetsIDLoose40==   njets   )                       histos["Mll"                      + hs]->Fill(Mll,             weight);
	if(njets<=0 && fMT2tree->NJetsIDLoose40>=abs(njets) )                       histos["Mll"                      + hs]->Fill(Mll,             weight);
	if(njets>=0 && fMT2tree->NJetsIDLoose40==    njets  && MT2lbV_withMlb>=0 )  histos["MT2lb_withMlbcut_woZveto" + hs]->Fill(MT2lbV_withMlb,  weight);
	if(njets< 0 && fMT2tree->NJetsIDLoose40>=abs(njets) && MT2lbV_withMlb>=0 )  histos["MT2lb_withMlbcut_woZveto" + hs]->Fill(MT2lbV_withMlb,  weight);
	if(njets>=0 && fMT2tree->NJetsIDLoose40==   njets   )                       histos["MT2lb_woZveto"            + hs]->Fill(MT2lbV,          weight);
	if(njets< 0 && fMT2tree->NJetsIDLoose40>=abs(njets) )                       histos["MT2lb_woZveto"            + hs]->Fill(MT2lbV,          weight);
	if(njets>=0 && fMT2tree->NJetsIDLoose40==    njets  && MT2lbV_massless>=0 ) histos["MT2lb_massless_woZveto"   + hs]->Fill(MT2lbV_massless, weight);
	if(njets< 0 && fMT2tree->NJetsIDLoose40>=abs(njets) && MT2lbV_massless>=0 ) histos["MT2lb_massless_woZveto"   + hs]->Fill(MT2lbV_massless, weight);
	if(njets>=0 && fMT2tree->NJetsIDLoose40==    njets  && MT2bbmin>=0 )        histos["MT2bb_woZveto"            + hs]->Fill(MT2bbmin,        weight);
	if(njets< 0 && fMT2tree->NJetsIDLoose40>=abs(njets) && MT2bbmin>=0 )        histos["MT2bb_woZveto"            + hs]->Fill(MT2bbmin,        weight);
	if(njets>=0 && fMT2tree->NJetsIDLoose40==    njets  && MT2bbmin_linMET>=0 ) histos["MT2bb_linMET_woZveto"     + hs]->Fill(MT2bbmin_linMET, weight);
	if(njets< 0 && fMT2tree->NJetsIDLoose40>=abs(njets) && MT2bbmin_linMET>=0 ) histos["MT2bb_linMET_woZveto"     + hs]->Fill(MT2bbmin_linMET, weight);
	if(njets>=0 && fMT2tree->NJetsIDLoose40==    njets  && MT2ll>=0 )           histos["MT2ll_woZveto"            + hs]->Fill(MT2ll,           weight);
	if(njets< 0 && fMT2tree->NJetsIDLoose40>=abs(njets) && MT2ll>=0 )           histos["MT2ll_woZveto"            + hs]->Fill(MT2ll,           weight);

	if(!(recoedeenZ) &&                   !(recoedmumunZ) && !(recoedemu) ) continue;//require of Zpeakveto

	bool truel1notb1(false), truel2notb2(false);//need to be defined outside
	//now real plots
	histos["NJetsID20" + hs]->Fill(njetsID20, weight);
	histos["NJetsID30" + hs]->Fill(njetsID30, weight);
	histos["NJetsID40" + hs]->Fill(njetsID40, weight);
	histos["NJetsID50" + hs]->Fill(njetsID50, weight);
	if(truettbarEE){
		if(genMl1b1>=0.) histos[ (string)"GenMlb_EE_"  +sampletype]->Fill(genMl1b1, weight);
		if(genMl2b2>=0.) histos[ (string)"GenMlb_EE_"  +sampletype]->Fill(genMl2b2, weight);
		if(genMT2>=0.)   histos[ (string)"GenMT2_EE_"  +sampletype]->Fill(genMT2,   weight);
		if(genM1>=0.)    histos[(string)"GenMtop_EE_"  +sampletype]->Fill(genM1,    weight);
		if(genM2>=0.)    histos[(string)"GenMtop_EE_"  +sampletype]->Fill(genM1,    weight);
	} if(truettbarEMu){
		if(genMl1b1>=0.) histos[ (string)"GenMlb_EMu_" +sampletype]->Fill(genMl1b1, weight);
		if(genMl2b2>=0.) histos[ (string)"GenMlb_EMu_" +sampletype]->Fill(genMl2b2, weight);
		if(genMT2>=0.)   histos[ (string)"GenMT2_EMu_" +sampletype]->Fill(genMT2,   weight);
		if(genM1>=0.)    histos[(string)"GenMtop_EMu_" +sampletype]->Fill(genM1,    weight);
		if(genM2>=0.)    histos[(string)"GenMtop_EMu_" +sampletype]->Fill(genM1,    weight);
	} if(truettbarMuMu){
		if(genMl1b1>=0.) histos[ (string)"GenMlb_MuMu_"+sampletype]->Fill(genMl1b1, weight);
		if(genMl2b2>=0.) histos[ (string)"GenMlb_MuMu_"+sampletype]->Fill(genMl2b2, weight);
		if(genMT2>=0.)   histos[ (string)"GenMT2_MuMu_"+sampletype]->Fill(genMT2,   weight);
		if(genM1>=0.)    histos[(string)"GenMtop_MuMu_"+sampletype]->Fill(genM1,    weight);
		if(genM2>=0.)    histos[(string)"GenMtop_MuMu_"+sampletype]->Fill(genM1,    weight);
	} if(truettbarEE){
		if(genMl1b1_allAcceptance>=0.) histos[ (string)"GenMlb_allAcceptance_EE_"  +sampletype]->Fill(genMl1b1_allAcceptance, weight);
		if(genMl2b2_allAcceptance>=0.) histos[ (string)"GenMlb_allAcceptance_EE_"  +sampletype]->Fill(genMl2b2_allAcceptance, weight);
		if(genMT2_allAcceptance>=0.)   histos[ (string)"GenMT2_allAcceptance_EE_"  +sampletype]->Fill(genMT2_allAcceptance,   weight);
		if(genM1_allAcceptance>=0.)    histos[(string)"GenMtop_allAcceptance_EE_"  +sampletype]->Fill(genM1_allAcceptance,    weight);
		if(genM2_allAcceptance>=0.)    histos[(string)"GenMtop_allAcceptance_EE_"  +sampletype]->Fill(genM1_allAcceptance,    weight);
	} if(truettbarEMu){
		if(genMl1b1_allAcceptance>=0.) histos[ (string)"GenMlb_allAcceptance_EMu_" +sampletype]->Fill(genMl1b1_allAcceptance, weight);
		if(genMl2b2_allAcceptance>=0.) histos[ (string)"GenMlb_allAcceptance_EMu_" +sampletype]->Fill(genMl2b2_allAcceptance, weight);
		if(genMT2_allAcceptance>=0.)   histos[ (string)"GenMT2_allAcceptance_EMu_" +sampletype]->Fill(genMT2_allAcceptance,   weight);
		if(genM1_allAcceptance>=0.)    histos[(string)"GenMtop_allAcceptance_EMu_" +sampletype]->Fill(genM1_allAcceptance,    weight);
		if(genM2_allAcceptance>=0.)    histos[(string)"GenMtop_allAcceptance_EMu_" +sampletype]->Fill(genM1_allAcceptance,    weight);
	} if(truettbarMuMu){
		if(genMl1b1_allAcceptance>=0.) histos[ (string)"GenMlb_allAcceptance_MuMu_"+sampletype]->Fill(genMl1b1_allAcceptance, weight);
		if(genMl2b2_allAcceptance>=0.) histos[ (string)"GenMlb_allAcceptance_MuMu_"+sampletype]->Fill(genMl2b2_allAcceptance, weight);
		if(genMT2_allAcceptance>=0.)   histos[ (string)"GenMT2_allAcceptance_MuMu_"+sampletype]->Fill(genMT2_allAcceptance,   weight);
		if(genM1_allAcceptance>=0.)    histos[(string)"GenMtop_allAcceptance_MuMu_"+sampletype]->Fill(genM1_allAcceptance,    weight);
		if(genM2_allAcceptance>=0.)    histos[(string)"GenMtop_allAcceptance_MuMu_"+sampletype]->Fill(genM1_allAcceptance,    weight);
	}

//new 14/06/2012
	if(abs(bjetsflavour[ind_l1bxMT2lb])==5&&abs(bjetsflavour[ind_l2bxMT2lb])==5)
		histos["MT2lb_bothbtags_truebjets"                   +hs]->Fill(MT2lbV,              weight);
	else if((abs(bjetsflavour[ind_l1bxMT2lb])!=5&&abs(bjetsflavour[ind_l2bxMT2lb])==5) || (abs(bjetsflavour[ind_l1bxMT2lb])==5&&abs(bjetsflavour[ind_l2bxMT2lb])!=5))
		histos["MT2lb_onebtag_truebjet"                      +hs]->Fill(MT2lbV,              weight);
	else if(abs(bjetsflavour[ind_l1bxMT2lb])!=5&&abs(bjetsflavour[ind_l2bxMT2lb])!=5)
		histos["MT2lb_nonebtags_truebjets"                   +hs]->Fill(MT2lbV,              weight);
	if(matchedgenp_tob1>=0 && matchedgenp_tob2>=0 && matchedgenp_tob1!=matchedgenp_tob2){//not match to same b
		histos["MT2lb_bothbtags_matchedtobfromtop"           +hs]->Fill(MT2lbV,              weight);
		if(matchedgenl_tol1>=0 && matchedgenl_tol2>=0 && matchedgenl_tol1!=matchedgenl_tol2) {//not match to same b
			histos["MT2lb_bothleptons_true_trueb"        +hs]->Fill(MT2lbV,              weight);
			histos["masslessMT2lb_bothleptons_true_trueb"+hs]->Fill(MT2lbV_massless,     weight);
			histos["MT2l_bothleptons_true_trueb"         +hs]->Fill(MT2ll,               weight);
			histos["Mlb_bothleptons_true_trueb"          +hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
			histos["Mlb_bothleptons_true_trueb"          +hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);
		} else if((matchedgenl_tol1==-1 && matchedgenl_tol2>=0)||(matchedgenl_tol1>=0 && matchedgenl_tol2==-1)) {
			histos["MT2lb_onelepton_true_trueb"          +hs]->Fill(MT2lbV,              weight);
			histos["masslessMT2lb_onelepton_true_trueb"  +hs]->Fill(MT2lbV_massless,     weight);
			histos["MT2l_onelepton_true_trueb"           +hs]->Fill(MT2ll,               weight);
			histos["Mlb_onelepton_true_trueb"            +hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
			histos["Mlb_onelepton_true_trueb"            +hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);
		} else if(matchedgenl_tol1==-1 && matchedgenl_tol2==-1) {
			histos["MT2lb_nolepton_true_trueb"           +hs]->Fill(MT2lbV,              weight);
			histos["masslessMT2lb_nolepton_true_trueb"   +hs]->Fill(MT2lbV_massless,     weight);
			histos["MT2l_nolepton_true_trueb"            +hs]->Fill(MT2ll,               weight);
			histos["Mlb_nolepton_true_trueb"             +hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
			histos["Mlb_nolepton_true_trueb"             +hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);
		}
	}
	else if((matchedgenp_tob1==-1 && matchedgenp_tob2>=0)||(matchedgenp_tob1>=0 && matchedgenp_tob2==-1))
		histos["MT2lb_onebtag_matchedtobfromtop"      +hs]->Fill(MT2lbV, weight);
	else if(matchedgenp_tob1==-1 && matchedgenp_tob2==-1)
		histos["MT2lb_nonebtags_matchedtobfromtop"    +hs]->Fill(MT2lbV, weight);
	else if(matchedgenp_tob1==matchedgenp_tob2 && matchedgenp_tob1!=-1) 
		histos["MT2lb_bothbtags_matchedtosamebfromtop"+hs]->Fill(MT2lbV, weight);
	else cout << __LINE__<< " somthing went wrong MT2 lb = " << MT2lbV << endl;

	if(matchedgenl_tol1>=0 && matchedgenl_tol2>=0 && matchedgenl_tol1!=matchedgenl_tol2){//not match to same b
		if(matchedgenp_tob1>=0 && matchedgenp_tob2>=0 && matchedgenp_tob1!=matchedgenp_tob2) {
			histos["MT2lb_bothb_true_truel"        +hs]->Fill(MT2lbV,              weight);
			histos["masslessMT2lb_bothb_true_truel"+hs]->Fill(MT2lbV_massless,     weight);
			histos["MT2l_bothb_true_truel"         +hs]->Fill(MT2ll,               weight);
			histos["Mlb_bothb_true_truel"          +hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
			histos["Mlb_bothb_true_truel"          +hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);
		} else if((matchedgenp_tob1==-1 && matchedgenp_tob2>=0)||(matchedgenp_tob1>=0 && matchedgenp_tob2==-1)) {
			histos["MT2lb_oneb_true_truel"         +hs]->Fill(MT2lbV,              weight);
			histos["masslessMT2lb_oneb_true_truel" +hs]->Fill(MT2lbV_massless,     weight);
			histos["MT2l_oneb_true_truel"          +hs]->Fill(MT2ll,               weight);
			histos["Mlb_oneb_true_truel"           +hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
			histos["Mlb_oneb_true_truel"           +hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);
		} else if(matchedgenp_tob1==-1 && matchedgenp_tob2==-1) {
			histos["MT2lb_nob_true_truel"          +hs]->Fill(MT2lbV,              weight);
			histos["masslessMT2lb_nob_true_truel"  +hs]->Fill(MT2lbV_massless,     weight);
			histos["MT2l_nob_true_truel"           +hs]->Fill(MT2ll,               weight);
			histos["Mlb_nob_true_truel"            +hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
			histos["Mlb_nob_true_truel"            +hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);
		}
	}
	if(abs(bjetsflavour[ind_l1bxMT2lb])==5&&abs(bjetsflavour[ind_l2bxMT2lb])==5){
		histos["Mlb_bothbtags_truebjets"        +hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
		histos["Mlb_bothbtags_truebjets"        +hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);}
	else if((abs(bjetsflavour[ind_l1bxMT2lb])!=5&&abs(bjetsflavour[ind_l2bxMT2lb])==5) || (abs(bjetsflavour[ind_l1bxMT2lb])==5&&abs(bjetsflavour[ind_l2bxMT2lb])!=5)){
		histos["Mlb_onebtag_truebjet"           +hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
		histos["Mlb_onebtag_truebjet"           +hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);}
	else if(abs(bjetsflavour[ind_l1bxMT2lb])!=5&&abs(bjetsflavour[ind_l2bxMT2lb])!=5){
		histos["Mlb_nonebtags_truebjets"        +hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
		histos["Mlb_nonebtags_truebjets"        +hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);}
	if(matchedgenp_tob1>=0 && matchedgenp_tob2>=0 && matchedgenp_tob1!=matchedgenp_tob2){//not match to same b
		histos["Mlb_bothbtags_matchedtobfromtop"+hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
		histos["Mlb_bothbtags_matchedtobfromtop"+hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);}
	else if((matchedgenp_tob1==-1 && matchedgenp_tob2>=0)||(matchedgenp_tob1>=0 && matchedgenp_tob2==-1)){
		histos["Mlb_onebtag_matchedtobfromtop"  +hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
		histos["Mlb_onebtag_matchedtobfromtop"  +hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);}
	else if(matchedgenp_tob1==-1 && matchedgenp_tob2==-1){
		histos["Mlb_nonebtags_matchedtobfromtop"+hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
		histos["Mlb_nonebtags_matchedtobfromtop"+hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);}
	else{
		if(!(matchedgenp_tob1==matchedgenp_tob2 && matchedgenp_tob1!=-1)) 
			cout << __LINE__<< " somthing went wrong Mlb = " << Ml1b[ind_l1bxMT2lb] << " & " << Ml2b[ind_l2bxMT2lb] << endl;
	}

	if(MT2lbV>200.){
		if(matchedgenp_tob1>=0){
			histos["RecoPt_over_GenPartonPt_truebfromtop_MT2lb_gt_200"+hs]->Fill(bjets[ind_l1bxMT2lb].Pt()/fMT2tree->genlept[matchedgenp_tob1].lv.Pt(), weight);
			histos["DeltaPt_bjet_matchedtrueparton_MT2lb_gt_200"      +hs]->Fill(fabs(bjets[ind_l1bxMT2lb].Pt()-fMT2tree->genlept[matchedgenp_tob1].lv.Pt()), weight);
			histos["DeltaR_bjet_matchedtrueparton_MT2lb_gt_200"       +hs]->Fill(bjets[ind_l1bxMT2lb].DeltaR(fMT2tree->genlept[matchedgenp_tob1].lv), weight);
		}
		if(matchedgenp_tob2>=0){
			histos["RecoPt_over_GenPartonPt_truebfromtop_MT2lb_gt_200"+hs]->Fill(bjets[ind_l2bxMT2lb].Pt()/fMT2tree->genlept[matchedgenp_tob2].lv.Pt(), weight);
			histos["DeltaPt_bjet_matchedtrueparton_MT2lb_gt_200"      +hs]->Fill(fabs(bjets[ind_l2bxMT2lb].Pt()-fMT2tree->genlept[matchedgenp_tob2].lv.Pt()), weight);
			histos["DeltaR_bjet_matchedtrueparton_MT2lb_gt_200"       +hs]->Fill(bjets[ind_l2bxMT2lb].DeltaR(fMT2tree->genlept[matchedgenp_tob2].lv), weight);
		}
		histos["BEta_MT2lb_gt_200"  +hs]->Fill(bjets[ind_l1bxMT2lb].Eta(),       weight);
		histos["BEta_MT2lb_gt_200"  +hs]->Fill(bjets[ind_l2bxMT2lb].Eta(),       weight);
		histos["BPt_MT2lb_gt_200"   +hs]->Fill(bjets[ind_l1bxMT2lb].Pt(),        weight);
		histos["BPt_MT2lb_gt_200"   +hs]->Fill(bjets[ind_l2bxMT2lb].Pt(),        weight);
		histos["BDiscr_MT2lb_gt_200"+hs]->Fill(bjetsdiscriminant[ind_l1bxMT2lb], weight);
		histos["BDiscr_MT2lb_gt_200"+hs]->Fill(bjetsdiscriminant[ind_l2bxMT2lb], weight);
		if(gb1>=0){
			histos["TrueBfromtopEta_MT2lb_gt_200"+hs]->Fill(fMT2tree->genlept[gb1].lv.Eta(), weight);
			histos["TrueBfromtopPt_MT2lb_gt_200" +hs]->Fill(fMT2tree->genlept[gb1].lv.Pt(),  weight);
		} if(gb2>=0){
			if(gb1>=0&&gl2>=0&&gl1>=0){
				if(MT2lbgen<0) cout << "MT2lbgen " << MT2lbgen << endl;
				histos["AltMT2lb_MT2lb_gt_200_usingtruebpartons"+hs]->Fill(MT2lbgen, weight);
			}
			histos["TrueBfromtopEta_MT2lb_gt_200"+hs]->Fill(fMT2tree->genlept[gb2].lv.Eta(), weight);
			histos["TrueBfromtopPt_MT2lb_gt_200" +hs]->Fill(fMT2tree->genlept[gb2].lv.Pt(),  weight);
		}
		if(fMT2tree->genmet[0].Pt()>0.001){
			histos["PFMEToverGenMET_MT2lb_gt_200"      +hs]->Fill(met.Pt()/fMT2tree->genmet[0].Pt(),       weight);
			histos["DeltaGenMET_PFMET_Phi_MT2lb_gt_200"+hs]->Fill(fabs(met.DeltaPhi(fMT2tree->genmet[0])), weight);
		}
		if((matchedgenp_tob1==-1 && matchedgenp_tob2>=0)||(matchedgenp_tob1>=0 && matchedgenp_tob2==-1)){
			int bshortind, blongind, truematchedb,truenotmatchedb, notselectedbutmatchedjet, otherb, otherbshort;
			if(matchedgenp_tob1>=0){//b1 is matched, probe b2
				bshortind = ind_l2bxMT2lb; blongind = bind[ind_l2bxMT2lb]; otherb = bind[ind_l1bxMT2lb]; otherbshort = ind_l1bxMT2lb;
				truematchedb = matchedgenp_tob1;
				truenotmatchedb = gb1; notselectedbutmatchedjet = trueb1match;
				if(gb1 == matchedgenp_tob1) { truenotmatchedb = gb2; notselectedbutmatchedjet = trueb2match; }
			}
			else if(matchedgenp_tob2>=0){//b2 is matched, probe b1
				bshortind = ind_l1bxMT2lb; blongind = bind[ind_l1bxMT2lb]; otherb = bind[ind_l2bxMT2lb]; otherbshort = ind_l2bxMT2lb;
				truematchedb = matchedgenp_tob2;
				truenotmatchedb = gb1; notselectedbutmatchedjet = trueb1match;
				if(gb1 == matchedgenp_tob2) { truenotmatchedb = gb2; notselectedbutmatchedjet = trueb2match; }
			}
			else cout << __LINE__ << "  this should not happen: matchedgenp_tob1 " << matchedgenp_tob1 << " matchedgenp_tob2 " << matchedgenp_tob2 << endl;
			//now fill
			if(truenotmatchedb>=0) histos["DeltaR_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(fabs(fMT2tree->genlept[truenotmatchedb].lv.DeltaR(fMT2tree->jet[blongind].lv)), weight);
			if(truenotmatchedb>=0) histos["BDiscr_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"            +hs] ->Fill(bjetsdiscriminant[bshortind],weight);
			if(truenotmatchedb>=0) histos["BEta_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"              +hs] ->Fill(fMT2tree->jet[blongind].lv.Eta(),weight);
			if(truenotmatchedb>=0) histos["BPt_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"               +hs] ->Fill(fMT2tree->jet[blongind].lv.Pt(),weight);
			if(truenotmatchedb>=0) histos["Eta_truebfromtop_notmatchedwithreco_MT2lb_gt_200"                       +hs] ->Fill(fMT2tree->genlept[truenotmatchedb].lv.Eta(),weight);
			if(truenotmatchedb>=0) histos["Pt_truebfromtop_notmatchedwithreco_MT2lb_gt_200"                        +hs] ->Fill(fMT2tree->genlept[truenotmatchedb].lv.Pt(),weight);
			if(truematchedb>=0) histos["BDiscr_matchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs]->Fill(bjetsdiscriminant[otherbshort],          weight);
			if(truematchedb>=0) histos["BEta_matchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"  +hs]->Fill(fMT2tree->jet[otherb].lv.Eta(),          weight);
			if(truematchedb>=0) histos["BPt_matchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"   +hs]->Fill(fMT2tree->jet[otherb].lv.Pt(),           weight);
			if(truematchedb>=0) histos["Eta_truebfromtop_matchedwithreco_MT2lb_gt_200"           +hs]->Fill(fMT2tree->genlept[truematchedb].lv.Eta(),weight);
			if(truematchedb>=0) histos["Pt_truebfromtop_matchedwithreco_MT2lb_gt_200"            +hs]->Fill(fMT2tree->genlept[truematchedb].lv.Pt() ,weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb) histos["BDiscr_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].bTagProbCSV,weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb) histos["BEta_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"  +hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].lv.Eta(),weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb) histos["BPt_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"  +hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].lv.Pt(),weight);
		}
		if(trueb1match<0 && trueb2match<0) histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(2.5,2.5,weight);//no matched jet
		else if(trueb1match<0){
			if(trueb2match==bind[ind_l1bxMT2lb] || trueb2match==bind[ind_l2bxMT2lb])//matched jet is also selected
				histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(2.5,0.5,weight);
			else    histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(2.5,1.5,weight);//matched jet is not selected
		}
		else if(trueb2match<0){
			if(trueb1match==bind[ind_l1bxMT2lb] || trueb1match==bind[ind_l2bxMT2lb])//matched jet is also selected
				histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(0.5,2.5,weight);
			else    histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(1.5,2.5,weight);//matched jet is not selected
		}
		else{//both genb have a matched jet
			if((trueb1match==bind[ind_l1bxMT2lb]&&trueb2match==bind[ind_l2bxMT2lb]) || (trueb2match==bind[ind_l1bxMT2lb]&&trueb1match==bind[ind_l2bxMT2lb]))
				histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(0.5,0.5,weight);//both matched jets are selected
			else if(trueb1match==bind[ind_l1bxMT2lb] || trueb1match==bind[ind_l2bxMT2lb])
				histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(0.5,1.5,weight);//one matched jet is selected, other not
			else if(trueb2match==bind[ind_l1bxMT2lb] || trueb2match==bind[ind_l2bxMT2lb])
				histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(1.5,0.5,weight);//one matched jet is selected, other not
			else 
				histos2["Matchingconfiguration_MT2lb_gt_200"+hs]->Fill(1.5,1.5,weight);//both matched jets are not selected
			if(trueb1match<0 || trueb2match<0) cout << "error line " << __LINE__ << endl;
		}
	} //MT2lb>200
	else {
		if(matchedgenp_tob1>=0){
			histos["RecoPt_over_GenPartonPt_truebfromtop_MT2lb_lt_200"+hs]->Fill(bjets[ind_l1bxMT2lb].Pt()/fMT2tree->genlept[matchedgenp_tob1].lv.Pt(), weight);
			histos["DeltaPt_bjet_matchedtrueparton_MT2lb_lt_200"      +hs]->Fill(fabs(bjets[ind_l1bxMT2lb].Pt()-fMT2tree->genlept[matchedgenp_tob1].lv.Pt()), weight);
			histos["DeltaR_bjet_matchedtrueparton_MT2lb_lt_200"       +hs]->Fill(bjets[ind_l1bxMT2lb].DeltaR(fMT2tree->genlept[matchedgenp_tob1].lv), weight);
		}
		if(matchedgenp_tob2>=0){
			histos["RecoPt_over_GenPartonPt_truebfromtop_MT2lb_lt_200"+hs]->Fill(bjets[ind_l2bxMT2lb].Pt()/fMT2tree->genlept[matchedgenp_tob2].lv.Pt(), weight);
			histos["DeltaPt_bjet_matchedtrueparton_MT2lb_lt_200"      +hs]->Fill(fabs(bjets[ind_l2bxMT2lb].Pt()-fMT2tree->genlept[matchedgenp_tob2].lv.Pt()), weight);
			histos["DeltaR_bjet_matchedtrueparton_MT2lb_lt_200"       +hs]->Fill(bjets[ind_l2bxMT2lb].DeltaR(fMT2tree->genlept[matchedgenp_tob2].lv), weight);
		}
		histos["BEta_MT2lb_lt_200"  +hs]->Fill(bjets[ind_l1bxMT2lb].Eta(),       weight);
		histos["BEta_MT2lb_lt_200"  +hs]->Fill(bjets[ind_l2bxMT2lb].Eta(),       weight);
		histos["BPt_MT2lb_lt_200"   +hs]->Fill(bjets[ind_l1bxMT2lb].Pt(),        weight);
		histos["BPt_MT2lb_lt_200"   +hs]->Fill(bjets[ind_l2bxMT2lb].Pt(),        weight);
		histos["BDiscr_MT2lb_lt_200"+hs]->Fill(bjetsdiscriminant[ind_l1bxMT2lb], weight);
		histos["BDiscr_MT2lb_lt_200"+hs]->Fill(bjetsdiscriminant[ind_l2bxMT2lb], weight);
		if(gb1>=0){
			histos["TrueBfromtopEta_MT2lb_lt_200"      +hs]->Fill(fMT2tree->genlept[gb1].lv.Eta(),         weight);
			histos["TrueBfromtopPt_MT2lb_lt_200"       +hs]->Fill(fMT2tree->genlept[gb1].lv.Pt(),          weight);
		} if(gb2>=0){
			histos["TrueBfromtopEta_MT2lb_lt_200"      +hs]->Fill(fMT2tree->genlept[gb2].lv.Eta(),         weight);
			histos["TrueBfromtopPt_MT2lb_lt_200"       +hs]->Fill(fMT2tree->genlept[gb2].lv.Pt(),          weight);
		}
		if(fMT2tree->genmet[0].Pt()>0.001){
			histos["PFMEToverGenMET_MT2lb_lt_200"      +hs]->Fill(met.Pt()/fMT2tree->genmet[0].Pt(),       weight);
			histos["DeltaGenMET_PFMET_Phi_MT2lb_lt_200"+hs]->Fill(fabs(met.DeltaPhi(fMT2tree->genmet[0])), weight);
		}
		if((matchedgenp_tob1==-1 && matchedgenp_tob2>=0)||(matchedgenp_tob1>=0 && matchedgenp_tob2==-1)){
			int bshortind, blongind, truematchedb,truenotmatchedb, notselectedbutmatchedjet, otherb, otherbshort;
			if(matchedgenp_tob1>=0){//b1 is matched, probe b2
				bshortind = ind_l2bxMT2lb; blongind = bind[ind_l2bxMT2lb]; otherb = bind[ind_l1bxMT2lb]; otherbshort = ind_l1bxMT2lb;
				truematchedb = matchedgenp_tob1;
				truenotmatchedb = gb1; notselectedbutmatchedjet = trueb1match;
				if(gb1 == matchedgenp_tob1) { truenotmatchedb = gb2; notselectedbutmatchedjet = trueb2match; }
			}
			else if(matchedgenp_tob2>=0){//b2 is matched, probe b1
				bshortind = ind_l1bxMT2lb; blongind = bind[ind_l1bxMT2lb]; otherb = bind[ind_l2bxMT2lb]; otherbshort = ind_l2bxMT2lb;
				truematchedb = matchedgenp_tob2;
				truenotmatchedb = gb1; notselectedbutmatchedjet = trueb1match;
				if(gb1 == matchedgenp_tob2) { truenotmatchedb = gb2; notselectedbutmatchedjet = trueb2match; }
			}
			else cout << __LINE__ << "  this should not happen: matchedgenp_tob1 " << matchedgenp_tob1 << " matchedgenp_tob2 " << matchedgenp_tob2 << endl;
			//now fill
			if(truenotmatchedb>=0) histos["DeltaR_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(fabs(fMT2tree->genlept[truenotmatchedb].lv.DeltaR(fMT2tree->jet[blongind].lv)),weight);
			if(truenotmatchedb>=0) histos["BDiscr_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(bjetsdiscriminant[bshortind],weight);
			if(truenotmatchedb>=0) histos["BEta_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"  +hs] ->Fill(fMT2tree->jet[blongind].lv.Eta(),weight);
			if(truenotmatchedb>=0) histos["BPt_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"   +hs] ->Fill(fMT2tree->jet[blongind].lv.Pt(),weight);
			if(truenotmatchedb>=0) histos["Eta_truebfromtop_notmatchedwithreco_MT2lb_lt_200"           +hs] ->Fill(fMT2tree->genlept[truenotmatchedb].lv.Eta(),weight);
			if(truenotmatchedb>=0) histos["Pt_truebfromtop_notmatchedwithreco_MT2lb_lt_200"            +hs] ->Fill(fMT2tree->genlept[truenotmatchedb].lv.Pt(),weight);
			if(truematchedb>=0) histos["BDiscr_matchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs]->Fill(bjetsdiscriminant[otherbshort],          weight);
			if(truematchedb>=0) histos["BEta_matchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"  +hs]->Fill(fMT2tree->jet[otherb].lv.Eta(),          weight);
			if(truematchedb>=0) histos["BPt_matchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"   +hs]->Fill(fMT2tree->jet[otherb].lv.Pt(),           weight);
			if(truematchedb>=0) histos["Eta_truebfromtop_matchedwithreco_MT2lb_lt_200"           +hs]->Fill(fMT2tree->genlept[truematchedb].lv.Eta(),weight);
			if(truematchedb>=0) histos["Pt_truebfromtop_matchedwithreco_MT2lb_lt_200"            +hs]->Fill(fMT2tree->genlept[truematchedb].lv.Pt(), weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb) histos["BDiscr_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"+hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].bTagProbCSV,weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb) histos["BEta_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"  +hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].lv.Eta(),weight);
			if(notselectedbutmatchedjet>=0 && notselectedbutmatchedjet!=otherb) histos["BPt_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"   +hs] ->Fill(fMT2tree->jet[notselectedbutmatchedjet].lv.Pt(),weight);
		}
		if(trueb1match<0 && trueb2match<0) histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(2.5,2.5,weight);//no matched jet
		else if(trueb1match<0){
			if(trueb2match==bind[ind_l1bxMT2lb] || trueb2match==bind[ind_l2bxMT2lb])//matched jet is also selected
				histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(2.5,0.5,weight);
			else    histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(2.5,1.5,weight);//matched jet is not selected
		}
		else if(trueb2match<0){
			if(trueb1match==bind[ind_l1bxMT2lb] || trueb1match==bind[ind_l2bxMT2lb])//matched jet is also selected
				histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(0.5,2.5,weight);
			else    histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(1.5,2.5,weight);//matched jet is not selected
		}
		else{//both genb have a matched jet
			if((trueb1match==bind[ind_l1bxMT2lb]&&trueb2match==bind[ind_l2bxMT2lb]) || (trueb2match==bind[ind_l1bxMT2lb]&&trueb1match==bind[ind_l2bxMT2lb]))
				histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(0.5,0.5,weight);//both matched jets are selected
			else if(trueb1match==bind[ind_l1bxMT2lb] || trueb1match==bind[ind_l2bxMT2lb])
				histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(0.5,1.5,weight);//one matched jet is selected, other not
			else if(trueb2match==bind[ind_l1bxMT2lb] || trueb2match==bind[ind_l2bxMT2lb])
				histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(1.5,0.5,weight);//one matched jet is selected, other not
			else 
				histos2["Matchingconfiguration_MT2lb_lt_200"+hs]->Fill(1.5,1.5,weight);//both matched jets are not selected
		}
	}//MT2lb <=200

	//new single top 05/07/2012
	if(fSamples[i].name=="Tbar_s"  || fSamples[i].name=="T_s" ) histos["SingleTop_MT2lb_schannel" +hs]->Fill(MT2lbV, weight);
	if(fSamples[i].name=="Tbar_t"  || fSamples[i].name=="T_t" ) histos["SingleTop_MT2lb_tchannel" +hs]->Fill(MT2lbV, weight);
	if(fSamples[i].name=="Tbar_tW" || fSamples[i].name=="T_tW") histos["SingleTop_MT2lb_tWchannel"+hs]->Fill(MT2lbV, weight);
	if(sampletype == "SingleTop"){
		if(matchedgenp_tob1>=0 && matchedgenp_tob2==-1){
			histos["SingleTop_BPt_truebfromtop" +hs]->Fill(bjets[ind_l1bxMT2lb].Pt(),  weight);
			histos["SingleTop_BEta_truebfromtop"+hs]->Fill(bjets[ind_l1bxMT2lb].Eta(), weight);
			histos["SingleTop_BPt_otherb"       +hs]->Fill(bjets[ind_l2bxMT2lb].Pt(),  weight);
			histos["SingleTop_BEta_otherb"      +hs]->Fill(bjets[ind_l2bxMT2lb].Eta(), weight);
		}
		else if(matchedgenp_tob2>=0 && matchedgenp_tob1==-1){
			histos["SingleTop_BPt_truebfromtop" +hs]->Fill(bjets[ind_l2bxMT2lb].Pt(),  weight);
			histos["SingleTop_BEta_truebfromtop"+hs]->Fill(bjets[ind_l2bxMT2lb].Eta(), weight);
			histos["SingleTop_BPt_otherb"       +hs]->Fill(bjets[ind_l1bxMT2lb].Pt(),  weight);
			histos["SingleTop_BEta_otherb"      +hs]->Fill(bjets[ind_l1bxMT2lb].Eta(), weight);
		}
		if(matchedgenl_tol1_regardlessFlavour>=0 && matchedgenl_tol2_regardlessFlavour==-1){
			histos["SingleTop_LPt_truelfromtop" +hs]->Fill(l1.Pt(),  weight);
			histos["SingleTop_LEta_truelfromtop"+hs]->Fill(l1.Eta(), weight);
			histos["SingleTop_LPt_otherl"       +hs]->Fill(l2.Pt(),  weight);
			histos["SingleTop_LEta_otherl"      +hs]->Fill(l2.Eta(), weight);
		}
		else if(matchedgenl_tol2_regardlessFlavour>=0 && matchedgenl_tol1_regardlessFlavour==-1){
			histos["SingleTop_LPt_truelfromtop" +hs]->Fill(l2.Pt(),  weight);
			histos["SingleTop_LEta_truelfromtop"+hs]->Fill(l2.Eta(), weight);
			histos["SingleTop_LPt_otherl"       +hs]->Fill(l1.Pt(),  weight);
			histos["SingleTop_LEta_otherl"      +hs]->Fill(l1.Eta(), weight);
		}
		if(abs(b1_matchedID)==5 && abs(b1_matchedMID)==6 && b1_matchedID*b1_matchedMID>0){
			if(abs(b2_matchedID)!=5 || abs(b2_matchedMID)!=6){
				if(b2_matched==false){
					histos["SingleTop_BID_otherb"  +hs] ->Fill(0.,  weight);
					histos["SingleTop_BMID_otherb" +hs] ->Fill(0.,  weight);
					histos["SingleTop_BGMID_otherb"+hs] ->Fill(0.,  weight);
				}
				else{
					histos["SingleTop_BID_otherb"  +hs] ->Fill(b2_matchedID,   weight);
					histos["SingleTop_BMID_otherb" +hs] ->Fill(b2_matchedMID,  weight);
					histos["SingleTop_BGMID_otherb"+hs] ->Fill(b2_matchedGMID, weight);
				}
			}
		}
		if(abs(b2_matchedID)==5 && abs(b2_matchedMID)==6 && b2_matchedID*b2_matchedMID>0){
			if(abs(b1_matchedID)!=5 || abs(b1_matchedMID)!=6){
				if(b1_matched==false){
					histos["SingleTop_BID_otherb"  +hs] ->Fill(0.,  weight);
					histos["SingleTop_BMID_otherb" +hs] ->Fill(0.,  weight);
					histos["SingleTop_BGMID_otherb"+hs] ->Fill(0.,  weight);
				}
				else{
					histos["SingleTop_BID_otherb"  +hs] ->Fill(b1_matchedID,   weight);
					histos["SingleTop_BMID_otherb" +hs] ->Fill(b1_matchedMID,  weight);
					histos["SingleTop_BGMID_otherb"+hs] ->Fill(b1_matchedGMID, weight);
				}
			}
		}
		if(abs(l1_matchedGMID)==6 && abs(l1_matchedMID)==24 && l1_matchedGMID*l1_matchedMID<0){
			if(abs(l2_matchedGMID)!=6){
				if(l2_matched==false){
					histos["SingleTop_LID_otherl"  +hs] ->Fill(0.,  weight);
					histos["SingleTop_LMID_otherl" +hs] ->Fill(0.,  weight);
					histos["SingleTop_LGMID_otherl"+hs] ->Fill(0.,  weight);
				}
				else{
					histos["SingleTop_LID_otherl"  +hs] ->Fill(l2_matchedID,   weight);
					histos["SingleTop_LMID_otherl" +hs] ->Fill(l2_matchedMID,  weight);
					histos["SingleTop_LGMID_otherl"+hs] ->Fill(l2_matchedGMID, weight);
				}
			}
		}
		if(abs(l2_matchedGMID)==6 && abs(l2_matchedMID)==24 && l2_matchedGMID*l1_matchedMID<0){
			if(abs(l1_matchedGMID)!=6){
				if(l1_matched==false){
					histos["SingleTop_LID_otherl"  +hs] ->Fill(0.,  weight);
					histos["SingleTop_LMID_otherl" +hs] ->Fill(0.,  weight);
					histos["SingleTop_LGMID_otherl"+hs] ->Fill(0.,  weight);
				}
				else{
					histos["SingleTop_LID_otherl"  +hs] ->Fill(l1_matchedID,   weight);
					histos["SingleTop_LMID_otherl" +hs] ->Fill(l1_matchedMID,  weight);
					histos["SingleTop_LGMID_otherl"+hs] ->Fill(l1_matchedGMID, weight);
				}
			}
		}
	}

	if(!l1_matched && !l2_matched){
		histos2["LOriginID"  +hs] ->Fill(0., 0., weight);
		histos2["LOriginMID" +hs] ->Fill(0., 0., weight);
		histos2["LOriginGMID"+hs] ->Fill(0., 0., weight);
	}
	else if(l1_matched && !l2_matched){
		histos2["LOriginID"  +hs] ->Fill((double)l1_matchedID,   0., weight);
		histos2["LOriginMID" +hs] ->Fill((double)l1_matchedMID,  0., weight);
		histos2["LOriginGMID"+hs] ->Fill((double)l1_matchedGMID, 0., weight);
	}
	else if(!l1_matched && l2_matched){
		histos2["LOriginID"  +hs] ->Fill(0., (double)l2_matchedID,   weight);
		histos2["LOriginMID" +hs] ->Fill(0., (double)l2_matchedMID,  weight);
		histos2["LOriginGMID"+hs] ->Fill(0., (double)l2_matchedGMID, weight);
	}
	else{
		histos2["LOriginID"  +hs] ->Fill((double)l1_matchedID,   (double)l2_matchedID,   weight);
		histos2["LOriginMID" +hs] ->Fill((double)l1_matchedMID,  (double)l2_matchedMID,  weight);
		histos2["LOriginGMID"+hs] ->Fill((double)l1_matchedGMID, (double)l2_matchedGMID, weight);
	}
	if(!b1_matched && !b2_matched){
		histos2["BOriginID"  +hs] ->Fill(0., 0., weight);
		histos2["BOriginMID" +hs] ->Fill(0., 0., weight);
		histos2["BOriginGMID"+hs] ->Fill(0., 0., weight);
	}
	else if(b1_matched && !b2_matched){
		histos2["BOriginID"  +hs] ->Fill((double)b1_matchedID,   0., weight);
		histos2["BOriginMID" +hs] ->Fill((double)b1_matchedMID,  0., weight);
		histos2["BOriginGMID"+hs] ->Fill((double)b1_matchedGMID, 0., weight);
	}
	else if(!b1_matched && b2_matched){
		histos2["BOriginID"  +hs] ->Fill(0., (double)b2_matchedID,   weight);
		histos2["BOriginMID" +hs] ->Fill(0., (double)b2_matchedMID,  weight);
		histos2["BOriginGMID"+hs] ->Fill(0., (double)b2_matchedGMID, weight);
	}
	else{
		histos2["BOriginID"  +hs] ->Fill((double)b1_matchedID,   (double)b2_matchedID,   weight);
		histos2["BOriginMID" +hs] ->Fill((double)b1_matchedMID,  (double)b2_matchedMID,  weight);
		histos2["BOriginGMID"+hs] ->Fill((double)b1_matchedGMID, (double)b2_matchedGMID, weight);
	}
	if(b1_matched){
		histos2["BOrigin"+hs] ->Fill(0.5, (double)b1_matchedID,   weight);
		histos2["BOrigin"+hs] ->Fill(1.5, (double)b1_matchedMID,  weight);
		histos2["BOrigin"+hs] ->Fill(2.5, (double)b1_matchedGMID, weight);
		histos2["BOrigin"+hs] ->Fill(3.5, (double)fMT2tree->jet[bind[ind_l1bxMT2lb] ].Flavour, weight);
	}
	else{
		histos2["BOrigin"+hs] ->Fill(0.5, 0., weight);
		histos2["BOrigin"+hs] ->Fill(1.5, 0., weight);
		histos2["BOrigin"+hs] ->Fill(2.5, 0., weight);
		histos2["BOrigin"+hs] ->Fill(3.5, (double)fMT2tree->jet[bind[ind_l1bxMT2lb] ].Flavour, weight);
	}
	if(b2_matched){
		histos2["BOrigin"+hs] ->Fill(0.5, (double)b2_matchedID,   weight);
		histos2["BOrigin"+hs] ->Fill(1.5, (double)b2_matchedMID,  weight);
		histos2["BOrigin"+hs] ->Fill(2.5, (double)b2_matchedGMID, weight);
		histos2["BOrigin"+hs] ->Fill(3.5, (double)fMT2tree->jet[bind[ind_l2bxMT2lb] ].Flavour, weight);
	}
	else{
		histos2["BOrigin"+hs] ->Fill(0.5, 0., weight);
		histos2["BOrigin"+hs] ->Fill(1.5, 0., weight);
		histos2["BOrigin"+hs] ->Fill(2.5, 0., weight);
		histos2["BOrigin"+hs] ->Fill(3.5, (double)fMT2tree->jet[bind[ind_l2bxMT2lb] ].Flavour, weight);
	}
	if(l1_matched){
		histos2["LOrigin"+hs] ->Fill(0.5, (double)l1_matchedID,   weight);
		histos2["LOrigin"+hs] ->Fill(1.5, (double)l1_matchedMID,  weight);
		histos2["LOrigin"+hs] ->Fill(2.5, (double)l1_matchedGMID, weight);
		histos2["LOrigin"+hs] ->Fill(3.5, (double)l1_flavour,     weight);
	}
	else{
		histos2["LOrigin"+hs] ->Fill(0.5, 0., weight);
		histos2["LOrigin"+hs] ->Fill(1.5, 0., weight);
		histos2["LOrigin"+hs] ->Fill(2.5, 0., weight);
		histos2["LOrigin"+hs] ->Fill(3.5, (double)l1_flavour, weight);
	}
	if(l2_matched){
		histos2["LOrigin"+hs] ->Fill(0.5, (double)l2_matchedID,   weight);
		histos2["LOrigin"+hs] ->Fill(1.5, (double)l2_matchedMID,  weight);
		histos2["LOrigin"+hs] ->Fill(2.5, (double)l2_matchedGMID, weight);
		histos2["LOrigin"+hs] ->Fill(3.5, (double)l2_flavour,     weight);
	}
	else{
		histos2["LOrigin"+hs] ->Fill(0.5, 0., weight);
		histos2["LOrigin"+hs] ->Fill(1.5, 0., weight);
		histos2["LOrigin"+hs] ->Fill(2.5, 0., weight);
		histos2["LOrigin"+hs] ->Fill(3.5, (double)l2_flavour, weight);
	}
// 
//rest
	for(unsigned int n = 0; n<jind.size(); ++n)               histos["AllBTagDiscr"              + hs]->Fill(fMT2tree->jet[jind[n] ].bTagProbCSV, weight);
	for(unsigned int n = 0; n<bjets.size(); ++n)              histos["AllBPt"                    + hs]->Fill(bjets[n].Pt(),                       weight);
	for(unsigned int n = 0; n<jind.size(); ++n)               histos["JetPt"                     + hs]->Fill(fMT2tree->jet[jind[n] ].lv.Pt(),     weight);
	if(l1e) histos["LepMT"  + hs]->Fill(fMT2tree->ele[0].MT,                  weight);//only leading lepton
	else    histos["LepMT"  + hs]->Fill(fMT2tree->muo[0].MT,                  weight);
	histos["LepPt"          + hs]->Fill(l1.Pt(),                              weight);
	histos["LepPt"          + hs]->Fill(l2.Pt(),                              weight);
	histos["BPt"            + hs]->Fill(bjets[ind_l1bxMT2lb].Pt(),            weight);
	histos["BPt"            + hs]->Fill(bjets[ind_l2bxMT2lb].Pt(),            weight);
	histos["NBJets20"       + hs]->Fill(nbjets20,                             weight);
	histos["NBJets30"       + hs]->Fill(nbjets30,                             weight);
	histos["NBJets40"       + hs]->Fill(nbjets40,                             weight);
	histos["NBJets50"       + hs]->Fill(nbjets50,                             weight);
	histos["BTagDiscr"      + hs]->Fill(bjetsdiscriminant[ind_l1bxMT2lb],     weight);
	histos["BTagDiscr"      + hs]->Fill(bjetsdiscriminant[ind_l2bxMT2lb],     weight);
	histos["HT"             + hs]->Fill(fMT2tree->misc.HT,                    weight);
	histos["MET"            + hs]->Fill(fMT2tree->misc.MET,                   weight);
	histos["MT2ll"          + hs]->Fill(MT2ll,                                weight);
	histos["MT2lb"          + hs]->Fill(MT2lbV,                               weight);
	histos["MT2lb_massless" + hs]->Fill(MT2lbV_massless,                      weight);
	if(Ml1b[ind_l1bxMT2lb]<180. && Ml2b[ind_l2bxMT2lb]<180. ) histos["MT2ll_withMlbcut"          + hs]->Fill(MT2ll, weight);
	if(MT2lbV_withMlb>=0 )                                    histos["MT2lb_withMlbcut"          + hs]->Fill(MT2lbV_withMlb,          weight);
	if(MT2lbV_massless_withMlb>=0 )                           histos["MT2lb_massless_withMlbcut" + hs]->Fill(MT2lbV_massless_withMlb, weight);
	if(MT2lbV_withMlb>=0 && fMT2tree->pileUp.NVertices<=5)    histos["MT2lb_withMlbcut_PUle5"    + hs]->Fill(MT2lbV_withMlb,          weight);
	if(MT2lbV_withMlb>=0 && fMT2tree->pileUp.NVertices>=9)    histos["MT2lb_withMlbcut_PUge9"    + hs]->Fill(MT2lbV_withMlb,          weight);
	if(MT2bbmin>=0.)                                          histos["MT2bb"                     + hs]->Fill(MT2bbmin,                weight);
	if(MT2bbmin_mW_linMET>=0.)                                histos["MT2bb_lInMET_mW"           + hs]->Fill(MT2bbmin_mW_linMET,      weight);
	if(MT2bbmin_linMET>=0.)                                   histos["MT2bb_lInMET"              + hs]->Fill(MT2bbmin_linMET,         weight);
	if(MT2bbmin_mW_linMET>=0. && Ml1b[ind_l1bxMT2lb]<180. && Ml2b[ind_l2bxMT2lb]<180.) histos["MT2bb_lInMET_mW_withMlbcut" + hs]->Fill(MT2bbmin_mW_linMET, weight);
	for(unsigned int n = 0; n<NumBJets; ++n){
		histos["MlbAll" + hs]->Fill(Ml1b[n],                              weight);
		histos["MlbAll" + hs]->Fill(Ml2b[n],                              weight);
	}
	histos["Mlb"            + hs]->Fill(Ml1b[ind_l1bxMT2lb],                  weight);
	histos["Mlb"            + hs]->Fill(Ml2b[ind_l2bxMT2lb],                  weight);
	histos["Mbb"            + hs]->Fill(Mbbarr[ind_l1bxMT2lb][ind_l2bxMT2lb], weight);
	for(int i1 =0; i1<int(NumBJets-1.); ++i1){ for(unsigned int i2 = i1+1; i2<NumBJets; ++i2){ histos["MbbAll" + hs]->Fill(Mbbarr[i1][i2], weight); } }//new
	histos["No_PV"          + hs]->Fill(fMT2tree->pileUp.NVertices,           weight);
	histos["UTM"            + hs]->Fill(UTMvec[ind_l1bxMT2lb][ind_l2bxMT2lb], weight);
	if(saveD2&&D2!=9999999.)                                                           histos["D2"                         + hs]->Fill(D2,                 weight);
	if(saveD2&&Ml1b[ind_l1bxMT2lb]<180. && Ml2b[ind_l2bxMT2lb]<180. &&D2!=9999999.)    histos["D2_withMlb"                 + hs]->Fill(D2,                 weight);
	for(unsigned int n=0; n<NumBJets; ++n){
		histos["DPhi_lbAll"         + hs]->Fill(DPhil1b[n],                              weight);
		histos["DPhi_lbAll"         + hs]->Fill(DPhil2b[n],                              weight);
		histos["DR_lbAll"           + hs]->Fill(DRl1b[n],                                weight);
		histos["DR_lbAll"           + hs]->Fill(DRl1b[n],                                weight);
		for(unsigned int n2 = n+1; n2<NumBJets; ++n2){
			histos["DPhi_bbAll" + hs]->Fill(DPhibbarr[n][n2],                        weight);
			histos["DR_bbAll"   + hs]->Fill(DRbbarr[n][n2],                          weight);
		}
	}
	histos["DPhi_lb"                    + hs]->Fill(DPhil1b[ind_l1bxMT2lb],                  weight);
	histos["DPhi_lb"                    + hs]->Fill(DPhil2b[ind_l2bxMT2lb],                  weight);
	histos["DR_lb"                      + hs]->Fill(DRl1b[ind_l1bxMT2lb],                    weight);
	histos["DR_lb"                      + hs]->Fill(DRl2b[ind_l2bxMT2lb],                    weight);
	histos["DPhi_ll"                    + hs]->Fill(DPhill,                                  weight);
	histos["DPhi_bb"                    + hs]->Fill(DPhibbarr[ind_l1bxMT2lb][ind_l2bxMT2lb], weight);
	histos["DR_ll"                      + hs]->Fill(DRll,                                    weight);
	histos["DR_bb"                      + hs]->Fill(DRbbarr[ind_l1bxMT2lb][ind_l2bxMT2lb],   weight);


	//distributions if other 'dicriminating' cut is applied
	if(MT2lbV_withMlb>=0   && MT2ll>=85.)            histos["MT2lb_withMlbcut_afterMT2llge85"         + hs]->Fill(MT2lbV_withMlb,          weight);
	if(MT2lbV>=0           && MT2ll>=85.)            histos["MT2lb_afterMT2llge85"                    + hs]->Fill(MT2lbV,                  weight);
	if(MT2lbV_massless>=0  && MT2ll>=85.)            histos["MT2lb_massless_afterMT2llge85"           + hs]->Fill(MT2lbV_massless,         weight);
	if(MT2lbV_massless_withMlb>=0 && MT2ll>=85.)     histos["MT2lb_massless_withMlbcut_afterMT2llge85"+ hs]->Fill(MT2lbV_massless_withMlb, weight);
	if(MT2lbV_withMlb>=0   && MT2bbmin_linMET>=150.) histos["MT2lb_withMlbcut_afterMT2bb_linMETge150" + hs]->Fill(MT2lbV_withMlb,          weight);
	if(MT2lbV>=0           && MT2bbmin_linMET>=150.) histos["MT2lb_afterMT2bb_linMETge150"            + hs]->Fill(MT2lbV,                  weight);
	if(MT2lbV_massless>=0  && MT2bbmin_linMET>=150.) histos["MT2lb_massless_afterMT2bb_linMETge150"   + hs]->Fill(MT2lbV_massless,         weight);
	if(MT2ll>=0            && MT2bbmin_linMET>=150.) histos["MT2ll_afterMT2bb_linMETge150"            + hs]->Fill(MT2ll,                   weight);
	if(MT2ll>=0            && MT2lbV>=200.)          histos["MT2ll_afterMT2lbge200"                   + hs]->Fill(MT2ll,                   weight);
	if(MT2ll>=0            && MT2lbV_withMlb>=180.)  histos["MT2ll_afterMT2lb_withMlbcutge180"        + hs]->Fill(MT2ll,                   weight);
	if(MT2ll>=0            && MT2lbV_massless>=140.) histos["MT2ll_afterMT2lb_masslessge140"          + hs]->Fill(MT2ll,                   weight);
	if(MT2bbmin_linMET>=0  && MT2ll>=85.)            histos["MT2bb_linMET_afterMT2llge85"             + hs]->Fill(MT2bbmin_linMET,         weight);
	if(MT2bbmin_linMET>=0  && MT2lbV>=200.)          histos["MT2bb_linMET_afterMT2lbge200"            + hs]->Fill(MT2bbmin_linMET,         weight);
	if(MT2bbmin_linMET>=0  && MT2lbV_withMlb>=180.)  histos["MT2bb_linMET_afterMT2lb_withMlbcutge180" + hs]->Fill(MT2bbmin_linMET,         weight);
	if(MT2bbmin_linMET>=0  && MT2lbV_massless>=140.) histos["MT2bb_linMET_afterMT2lb_masslessge140"   + hs]->Fill(MT2bbmin_linMET,         weight);
	if(saveD2&&D2!=9999999.&& MT2ll>=85.)            histos["D2_afterMT2llge85"                       + hs]->Fill(D2,                      weight);
	if(saveD2&&D2!=9999999.&& MT2bbmin_linMET>=150.) histos["D2_afterMT2bb_linMETge150"               + hs]->Fill(D2,                      weight);
	if(saveD2&&D2!=9999999.&& MT2lbV>=200.)          histos["D2_afterMT2lbge200"                      + hs]->Fill(D2,                      weight);
	if(saveD2&&D2!=9999999.&& MT2lbV_withMlb>=180.)  histos["D2_afterMT2lb_withMlbcutge180"           + hs]->Fill(D2,                      weight);
	if(saveD2&&D2!=9999999.&& MT2lbV_massless>=140.) histos["D2_afterMT2lb_masslessge140"             + hs]->Fill(D2,                      weight);

	for(unsigned int i1 = 0; i1<NumBJets; ++ i1){
	for(unsigned int i2 = 0; i2<NumBJets; ++ i2){
		if(i1==i2) continue;
		if(truel1bxtop[i1])                      histos["Mlb_genmatchingTop"  + hs]->Fill(Ml1b[i1],        weight);
		if(truel2bxtop[i2])                      histos["Mlb_genmatchingTop"  + hs]->Fill(Ml2b[i2],        weight);
		if(truel1bxtop[i1] && truel2bxtop[i2])   histos["MT2_genmatchingTops" + hs]->Fill(MT2_l1b[i1][i2], weight);
		if((int)i1!=ind_l1bxMT2lb && truel1bxtop[i1]) truel1notb1 = true;
		if((int)i2!=ind_l2bxMT2lb && truel2bxtop[i2]) truel2notb2 = true;
	}}
	histos2["MT2lb_vs_UTM"                               + hs]->Fill(MT2lbV,         UTMvec[ind_l1bxMT2lb][ind_l2bxMT2lb],       weight);
	if(MT2lbV_withMlb>=0) histos2["MT2lb_withMlb_vs_UTM" + hs]->Fill(MT2lbV_withMlb, UTMvec[ind_l1bxMT2lbMlb][ind_l1bxMT2lbMlb], weight);
	//x is l1b1,l1b2, y is l2b1,l2b2
	if(DRl1b[ind_l1bxMT2lb]<DRl1b[ind_l2bxMT2lb] && DRl2b[ind_l1bxMT2lb]<DRl2b[ind_l2bxMT2lb])        histos2["DR_lb_configuration"    + hs]->Fill(0.5,0.5, weight);
	else if(DRl1b[ind_l1bxMT2lb]< DRl1b[ind_l2bxMT2lb] && DRl2b[ind_l1bxMT2lb]>=DRl2b[ind_l2bxMT2lb]) histos2["DR_lb_configuration"    + hs]->Fill(0.5,1.5, weight);
	else if(DRl1b[ind_l1bxMT2lb]>=DRl1b[ind_l2bxMT2lb] && DRl2b[ind_l1bxMT2lb]< DRl2b[ind_l2bxMT2lb]) histos2["DR_lb_configuration"    + hs]->Fill(1.5,0.5, weight);
	else if(DRl1b[ind_l1bxMT2lb]>=DRl1b[ind_l2bxMT2lb] && DRl2b[ind_l1bxMT2lb]>=DRl2b[ind_l2bxMT2lb]) histos2["DR_lb_configuration"    + hs]->Fill(1.5,1.5, weight);
	if(abs(bjetsflavour[ind_l1bxMT2lb])==5)                                                           histos2["trueb1b2_configuration" + hs]->Fill(0.5,0.5, weight);
	else                                                                                              histos2["trueb1b2_configuration" + hs]->Fill(0.5,1.5, weight);
	if(abs(bjetsflavour[ind_l2bxMT2lb])==5)                                                           histos2["trueb1b2_configuration" + hs]->Fill(1.5,0.5, weight);
	else                                                                                              histos2["trueb1b2_configuration" + hs]->Fill(1.5,1.5, weight);
	if(truel1bxtop[ind_l1bxMT2lb] && truel2bxtop[ind_l2bxMT2lb])                                      histos2["trueTop_configuration"  + hs]->Fill(0.5,0.5, weight);
	if(truel1bxtop[ind_l1bxMT2lb] && truel2notb2)                                                     histos2["trueTop_configuration"  + hs]->Fill(0.5,1.5, weight);
	if(truel1bxtop[ind_l1bxMT2lb] && !(truel2bxtop[ind_l2bxMT2lb]) && !(truel2notb2))                 histos2["trueTop_configuration"  + hs]->Fill(0.5,2.5, weight);
	if(truel1notb1 && truel2bxtop[ind_l2bxMT2lb])                                                     histos2["trueTop_configuration"  + hs]->Fill(1.5,0.5, weight);
	if(truel1notb1 && truel2notb2)                                                                    histos2["trueTop_configuration"  + hs]->Fill(1.5,1.5, weight);
	if(truel1notb1 && !(truel2bxtop[ind_l2bxMT2lb]) && !(truel2notb2))                                histos2["trueTop_configuration"  + hs]->Fill(1.5,2.5, weight);
	if(!(truel1bxtop[ind_l1bxMT2lb]) && !(truel1notb1) && truel2bxtop[ind_l2bxMT2lb])                 histos2["trueTop_configuration"  + hs]->Fill(2.5,0.5, weight);
	if(!(truel1bxtop[ind_l1bxMT2lb]) && !(truel1notb1) && truel2notb2)                                histos2["trueTop_configuration"  + hs]->Fill(2.5,1.5, weight);
	if(!(truel1bxtop[ind_l1bxMT2lb] || truel1notb1 || truel2bxtop[ind_l2bxMT2lb] || truel2notb2) )    histos2["trueTop_configuration"  + hs]->Fill(2.5,2.5, weight);
	if(MT2lbV_massless_withMlb>=0 && MT2ll>=0) 
		               histos2["MT2lb_massless_withMlb_vs_MT2ll"     + hs]->Fill(MT2lbV_massless_withMlb, MT2ll,                                        weight);
	if(MT2lbV_withMlb>=0)  histos2["MT2lb_withMlbcut_vs_MT2ll"           + hs]->Fill(MT2lbV_withMlb,  MT2ll,                                                weight);
	if(MT2lbV_withMlb>=0)  histos2["MT2lb_withMlbcut_vs_MT2bb"           + hs]->Fill(MT2lbV_withMlb,  MT2bbarr[ind_l1bxMT2lbMlb][ind_l2bxMT2lbMlb],         weight);
	if(MT2lbV_withMlb>=0)  histos2["MT2lb_withMlbcut_vs_MT2bb_lInMET"    + hs]->Fill(MT2lbV_withMlb,  MT2bb_linMET[ind_l1bxMT2lbMlb][ind_l2bxMT2lbMlb],     weight);
	if(MT2lbV_withMlb>=0)  histos2["MT2lb_withMlbcut_vs_MT2bb_lInMET_mW" + hs]->Fill(MT2lbV_withMlb,  MT2bb_mW_linMET[ind_l1bxMT2lbMlb][ind_l2bxMT2lbMlb],  weight);
	if(MT2lbV_withMlb>=0)  histos2["MT2lb_withMlbcut_vs_MET"             + hs]->Fill(MT2lbV_withMlb,  fMT2tree->misc.MET,                                   weight);
	if(MT2lbV_withMlb>=0)  histos2["MT2lb_withMlbcut_vs_HT"              + hs]->Fill(MT2lbV_withMlb,  fMT2tree->misc.HT,                                    weight);
	if(MT2lbV_withMlb>=0)  histos2["MT2lb_withMlbcut_vs_PV"              + hs]->Fill(MT2lbV_withMlb,  fMT2tree->pileUp.NVertices,                           weight);
	if(saveD2&&MT2lbV_withMlb>=0) histos2["MT2lb_withMlbcut_vs_D2"       + hs]->Fill(MT2lbV_withMlb,  D2,                                                   weight);
	if(MT2lbV_withMlb>=0)  histos2["MT2lb_withMlbcut_vs_MT2lb_massless"  + hs]->Fill(MT2lbV_withMlb,  MT2_massless_l1b[ind_l1bxMT2lbMlb][ind_l2bxMT2lbMlb], weight);
	if(MT2lbV>=0)          histos2["MT2lb_vs_MT2lb_massless"             + hs]->Fill(MT2lbV,          MT2_massless_l1b[ind_l1bxMT2lb][ind_l2bxMT2lb],       weight);
	if(MT2lbV>=0)          histos2["MT2lb_vs_MT2ll"                      + hs]->Fill(MT2lbV,          MT2ll,                                                weight);
	if(MT2lbV>=0)          histos2["MT2lb_vs_MT2bb_lInMET_mW"            + hs]->Fill(MT2lbV,          MT2bb_mW_linMET[ind_l1bxMT2lb][ind_l2bxMT2lb],        weight);
	if(MT2lbV>=0)          histos2["MT2lb_vs_MT2bb_lInMET"               + hs]->Fill(MT2lbV,          MT2bb_linMET[ind_l1bxMT2lb][ind_l2bxMT2lb],           weight);
	if(MT2lbV>=0)          histos2["MT2lb_vs_MT2bb"                      + hs]->Fill(MT2lbV,          MT2bbarr[ind_l1bxMT2lb][ind_l2bxMT2lb],               weight);
	if(MT2lbV>=0)          histos2["MT2lb_vs_MET"                        + hs]->Fill(MT2lbV,          fMT2tree->misc.MET,                                   weight);
	if(MT2lbV>=0)          histos2["MT2lb_vs_HT"                         + hs]->Fill(MT2lbV,          fMT2tree->misc.HT,                                    weight);
	if(MT2lbV>=0)          histos2["MT2lb_vs_PV"                         + hs]->Fill(MT2lbV,          fMT2tree->pileUp.NVertices,                           weight);
	if(saveD2&&MT2lbV>=0)  histos2["MT2lb_vs_D2"                         + hs]->Fill(MT2lbV,          D2,                                                   weight);
	if(MT2lbV>=0)          histos2["MT2lb_vs_Mlb"                        + hs]->Fill(MT2lbV,          Ml1b[ind_l1bxMT2lb],                                  weight);
	if(MT2lbV>=0)          histos2["MT2lb_vs_Mlb"                        + hs]->Fill(MT2lbV,          Ml2b[ind_l2bxMT2lb],                                  weight);
	if(MT2lbV_massless>=0) histos2["MT2lb_massless_vs_MT2ll"             + hs]->Fill(MT2lbV_massless, MT2ll,                                                weight);
	if(MT2lbV_massless>=0) histos2["MT2lb_massless_vs_MT2bb_lInMET_mW"   + hs]->Fill(MT2lbV_massless, MT2bb_mW_linMET[ind_l1bxMT2lbMassless][ind_l2bxMT2lbMassless], weight);
	if(MT2lbV_massless>=0) histos2["MT2lb_massless_vs_MT2bb_lInMET"      + hs]->Fill(MT2lbV_massless, MT2bb_linMET[ind_l1bxMT2lbMassless][ind_l2bxMT2lbMassless],    weight);
	if(MT2lbV_massless>=0) histos2["MT2lb_massless_vs_MT2bb"             + hs]->Fill(MT2lbV_massless, MT2bbarr[ind_l1bxMT2lbMassless][ind_l2bxMT2lbMassless],        weight);//new
	if(MT2lbV_massless>=0) histos2["MT2lb_massless_vs_MET"               + hs]->Fill(MT2lbV_massless, fMT2tree->misc.MET,                                   weight);
	if(MT2lbV_massless>=0) histos2["MT2lb_massless_vs_HT"                + hs]->Fill(MT2lbV_massless, fMT2tree->misc.HT,                                    weight);
	if(MT2lbV_massless>=0) histos2["MT2lb_massless_vs_PV"                + hs]->Fill(MT2lbV_massless, fMT2tree->pileUp.NVertices,                           weight);
	if(saveD2&&MT2lbV_massless>=0) histos2["MT2lb_massless_vs_D2"        + hs]->Fill(MT2lbV_massless, D2,                                                   weight);
	if(MT2lbV_massless>=0) histos2["MT2lb_massless_vs_Mlb"               + hs]->Fill(MT2lbV_massless, Ml1b[ind_l1bxMT2lbMassless],                          weight);
	if(MT2lbV_massless>=0) histos2["MT2lb_massless_vs_Mlb"               + hs]->Fill(MT2lbV_massless, Ml2b[ind_l2bxMT2lbMassless],                          weight);
	if(MT2lbV_massless>=0) histos2["MT2lb_massless_vs_BTagDiscr"         + hs]->Fill(MT2lbV_massless, bjetsdiscriminant[ind_l1bxMT2lbMassless],             weight);
	if(MT2lbV_massless>=0) histos2["MT2lb_massless_vs_BTagDiscr"         + hs]->Fill(MT2lbV_massless, bjetsdiscriminant[ind_l2bxMT2lbMassless],             weight);
	                       histos2["MT2ll_vs_MT2bb_lInMET_mW"            + hs]->Fill(MT2ll,           MT2bbmin_mW_linMET,                                   weight);
	                       histos2["MT2ll_vs_MT2bb_lInMET"               + hs]->Fill(MT2ll,           MT2bbmin_linMET,                                      weight);
	                       histos2["MT2ll_vs_MT2bb"                      + hs]->Fill(MT2ll,           MT2bbmin,                                             weight);
	if(saveD2)             histos2["MT2ll_vs_D2"                         + hs]->Fill(MT2ll,           D2,                                                   weight);
	if(saveD2)             histos2["MT2bb_vs_D2"                         + hs]->Fill(MT2bbmin,        D2,                                                   weight);
	if(MT2lbV>=0)          histos2["MT2ll_vs_Mlb"                        + hs]->Fill(MT2ll,           Ml1b[ind_l1bxMT2lb],                                  weight);
	if(MT2lbV>=0)          histos2["MT2ll_vs_Mlb"                        + hs]->Fill(MT2ll,           Ml2b[ind_l2bxMT2lb],                                  weight);
	if(MT2lbV>=0) histos2["MT2bb_vs_Mlb"        + hs]->Fill(MT2bbarr[ind_l1bxMT2lb][ind_l2bxMT2lb],       Ml1b[ind_l1bxMT2lb],                              weight);
	if(MT2lbV>=0) histos2["MT2bb_vs_Mlb"        + hs]->Fill(MT2bbarr[ind_l1bxMT2lb][ind_l2bxMT2lb],       Ml2b[ind_l2bxMT2lb],                              weight);
	if(MT2lbV>=0) histos2["MT2bb_lInMET_vs_Mlb" + hs]->Fill(MT2bb_linMET[ind_l1bxMT2lb][ind_l2bxMT2lb],   Ml1b[ind_l1bxMT2lb],                              weight);
	if(MT2lbV>=0) histos2["MT2bb_lInMET_vs_Mlb" + hs]->Fill(MT2bb_linMET[ind_l1bxMT2lb][ind_l2bxMT2lb],   Ml2b[ind_l2bxMT2lb],                              weight);
	                       histos2["MT2bb_vs_BTagDiscr"                  + hs]->Fill(MT2bbmin,            bjetsdiscriminant[ind_b1x_MT2bbmin],              weight);
	                       histos2["MT2bb_vs_BTagDiscr"                  + hs]->Fill(MT2bbmin,            bjetsdiscriminant[ind_b2x_MT2bbmin],              weight);
	                       histos2["MT2lb_vs_BTagDiscr"                  + hs]->Fill(MT2lbV,              bjetsdiscriminant[ind_l1bxMT2lb],                 weight);
	                       histos2["MT2lb_vs_BTagDiscr"                  + hs]->Fill(MT2lbV,              bjetsdiscriminant[ind_l2bxMT2lb],                 weight);
	                       histos2["Mlb_vs_BTagDiscr"                    + hs]->Fill(Ml1b[ind_l1bxMT2lb], bjetsdiscriminant[ind_l1bxMT2lb],                 weight);
	                       histos2["Mlb_vs_BTagDiscr"                    + hs]->Fill(Ml2b[ind_l2bxMT2lb], bjetsdiscriminant[ind_l2bxMT2lb],                 weight);
	if(saveD2)             histos2["D2_vs_BTagDiscr"                     + hs]->Fill(D2,                  bjetsdiscriminant[ind_l1bxMT2lb],                 weight);
	if(saveD2)             histos2["D2_vs_BTagDiscr"                     + hs]->Fill(D2,                  bjetsdiscriminant[ind_l2bxMT2lb],                 weight);
	                       histos2["PUWeight_vs_PV"                      + hs]->Fill(fMT2tree->pileUp.Weight,   fMT2tree->pileUp.NVertices,                 weight);

	histos3["MT2lb_MT2ll_MT2bb"                   + hs]->Fill(MT2lbV,          MT2ll, MT2bbarr[ind_l1bxMT2lb][ind_l2bxMT2lb],                        weight);
	histos3["MT2lb_MT2ll_MT2bb_lInMET"            + hs]->Fill(MT2lbV,          MT2ll, MT2bb_linMET[ind_l1bxMT2lb][ind_l2bxMT2lb],                    weight);
	histos3["MT2lb_MT2ll_MT2bb_lInMET_mW"         + hs]->Fill(MT2lbV,          MT2ll, MT2bb_mW_linMET[ind_l1bxMT2lb][ind_l2bxMT2lb],                 weight);
	histos3["MT2lbMassless_MT2ll_MT2bb"           + hs]->Fill(MT2lbV_massless, MT2ll, MT2bbarr[ind_l1bxMT2lbMassless][ind_l2bxMT2lbMassless],        weight);
	histos3["MT2lbMassless_MT2ll_MT2bb_lInMET"    + hs]->Fill(MT2lbV_massless, MT2ll, MT2bb_linMET[ind_l1bxMT2lbMassless][ind_l2bxMT2lbMassless],    weight);
	histos3["MT2lbMassless_MT2ll_MT2bb_lInMET_mW" + hs]->Fill(MT2lbV_massless, MT2ll, MT2bb_mW_linMET[ind_l1bxMT2lbMassless][ind_l2bxMT2lbMassless], weight);

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
		//ADD FUNCTIONS TO HTSTack AND TLEGEND
	}

	cout << "Adding MuMu, EE, EMu to LL" << endl;
	for(unsigned int n =0; n<histonames.size(); ++n){
	   for(int is = 0; is<sampletypesize; ++is){
		string hs1 = string("_") + string("EE")   + string("_") + sample_type[is];
		string hs2 = string("_") + string("EMu")  + string("_") + sample_type[is];
		string hs3 = string("_") + string("MuMu") + string("_") + sample_type[is];
		string mapname = histonames[n];
		string h   = string("_") + string("LL") + string("_") + sample_type[is];
		if(histos.count(mapname+hs1)>0)  histos[mapname+h] ->Add(histos[mapname+hs1], 1);
		if(histos.count(mapname+hs2)>0)  histos[mapname+h] ->Add(histos[mapname+hs2], 1);
		if(histos.count(mapname+hs3)>0)  histos[mapname+h] ->Add(histos[mapname+hs3], 1);
		if(histos2.count(mapname+hs1)>0) histos2[mapname+h]->Add(histos2[mapname+hs1],1);
		if(histos2.count(mapname+hs2)>0) histos2[mapname+h]->Add(histos2[mapname+hs2],1);
		if(histos2.count(mapname+hs3)>0) histos2[mapname+h]->Add(histos2[mapname+hs3],1);
		if(histos3.count(mapname+hs1)>0) histos3[mapname+h]->Add(histos3[mapname+hs1],1);
		if(histos3.count(mapname+hs2)>0) histos3[mapname+h]->Add(histos3[mapname+hs2],1);
		if(histos3.count(mapname+hs3)>0) histos3[mapname+h]->Add(histos3[mapname+hs3],1);
		string hsf = string("_") + string("SF") + string("_") + sample_type[is];
		if(histos.count(mapname+hs1)>0)  histos[mapname+hsf] ->Add(histos[mapname+hs1], 1);
		if(histos.count(mapname+hs3)>0)  histos[mapname+hsf] ->Add(histos[mapname+hs3], 1);
		if(histos2.count(mapname+hs1)>0) histos2[mapname+hsf]->Add(histos2[mapname+hs1],1);
		if(histos2.count(mapname+hs3)>0) histos2[mapname+hsf]->Add(histos2[mapname+hs3],1);
		if(histos3.count(mapname+hs1)>0) histos3[mapname+hsf]->Add(histos3[mapname+hs1],1);
		if(histos3.count(mapname+hs3)>0) histos3[mapname+hsf]->Add(histos3[mapname+hs3],1);
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
		string h  = string("_") + lepton_type[il] + string("_") + (string)"mc";

		if(histos.count(mapname+hs)>0)  histos[mapname+h] ->Add(histos[mapname+hs], 1);
		if(histos2.count(mapname+hs)>0) histos2[mapname+h]->Add(histos2[mapname+hs],1);
		if(histos3.count(mapname+hs)>0) histos3[mapname+h]->Add(histos3[mapname+hs],1);

	   }}
	}
	cout << "setting stack color" << endl;//several options implemented as nobody liked the default ones
	Legend1 -> SetFillColor(0);
   	Legend1 -> SetBorderSize(0);
	for(unsigned int n = 0; n<histonames.size(); ++n){
	   	for(int il = 0; il<leptontypesize; ++il){
		for(int is = 0; is<sampletypesize; ++is){
			string hs = string("_") + lepton_type[il] + string("_") + sample_type[is];
			Color_t colour = 603;//set default others if we have some failure
			//v0 0th option
/* */			     if(sample_type[is]=="QCD")       colour = 880;
			else if(sample_type[is]=="WJets")     colour = 417;
			else if(sample_type[is]=="ZJets")     colour = 419;
			else if(sample_type[is]=="TTbar")     colour = 401;
			else if(sample_type[is]=="SingleTop") colour = 595;
			else if(sample_type[is]=="TTbarV")    colour = 65;
			else if(sample_type[is]=="VV/VVV")    colour = 635;
			else if(sample_type[is]=="Other")     colour = 603;
			else if(sample_type[is]=="mc")        colour = 603;
/*			//v0b Z <--> W 5th option
			     if(sample_type[is]=="QCD")       colour = 880;
			else if(sample_type[is]=="WJets")     colour = 419;
			else if(sample_type[is]=="ZJets")     colour = 417;
			else if(sample_type[is]=="TTbar")     colour = 401;
			else if(sample_type[is]=="SingleTop") colour = 595;
			else if(sample_type[is]=="TTbarV")    colour = 65;
			else if(sample_type[is]=="VV/VVV")    colour = 635;
			else if(sample_type[is]=="Other")     colour = 603;
			else if(sample_type[is]=="mc")        colour = 603;
			//v1 top blau 1st option, close with 2nd option
			     if(sample_type[is]=="QCD")       colour = 880;
			else if(sample_type[is]=="WJets")     colour = 401;
			else if(sample_type[is]=="ZJets")     colour = 800;
			else if(sample_type[is]=="TTbar")     colour = 65;
			else if(sample_type[is]=="SingleTop") colour = 855;
			else if(sample_type[is]=="TTbarV")    colour = 603;
			else if(sample_type[is]=="VV/VVV")    colour = 635;
			else if(sample_type[is]=="Other")     colour = 419;
			else if(sample_type[is]=="mc")        colour = 603;
			//v2 top orange  2nd option, close with 1st option
			     if(sample_type[is]=="QCD")       colour = 880;
			else if(sample_type[is]=="WJets")     colour = 419;
			else if(sample_type[is]=="ZJets")     colour = 65;
			else if(sample_type[is]=="TTbar")     colour = 800;
			else if(sample_type[is]=="SingleTop") colour = 417;
			else if(sample_type[is]=="TTbarV")    colour = 808;
			else if(sample_type[is]=="VV/VVV")    colour = 635;
			else if(sample_type[is]=="Other")     colour = 223;
			else if(sample_type[is]=="mc")        colour = 603;
			//v3 top rot 4th option
			     if(sample_type[is]=="QCD")       colour = 880;
			else if(sample_type[is]=="WJets")     colour = 215;
			else if(sample_type[is]=="ZJets")     colour = 65;
			else if(sample_type[is]=="TTbar")     colour = 625;
			else if(sample_type[is]=="SingleTop") colour = 613;
			else if(sample_type[is]=="TTbarV")    colour = 882;
			else if(sample_type[is]=="VV/VVV")    colour = 417;
			else if(sample_type[is]=="Other")     colour = 401;
			else if(sample_type[is]=="mc")        colour = 603;
			//v4 top gruen 3rd option
			     if(sample_type[is]=="QCD")       colour = 603;
			else if(sample_type[is]=="WJets")     colour = 797;
			else if(sample_type[is]=="ZJets")     colour = 401;
			else if(sample_type[is]=="TTbar")     colour = 417;
			else if(sample_type[is]=="SingleTop") colour = 419;
			else if(sample_type[is]=="TTbarV")    colour = 65;
			else if(sample_type[is]=="VV/VVV")    colour = 635;
			else if(sample_type[is]=="Other")     colour = 600;
			else if(sample_type[is]=="mc")        colour = 603;
			//v5 top violett 6th option
			     if(sample_type[is]=="QCD")       colour = 600;
			else if(sample_type[is]=="WJets")     colour = 419;
			else if(sample_type[is]=="ZJets")     colour = 65;
			else if(sample_type[is]=="TTbar")     colour = 223;
			else if(sample_type[is]=="SingleTop") colour = 871;
			else if(sample_type[is]=="TTbarV")    colour = 600;
			else if(sample_type[is]=="VV/VVV")    colour = 417;
			else if(sample_type[is]=="Other")     colour = 401;
			else if(sample_type[is]=="mc")        colour = 603;
 */
			if(sample_type[is]=="QCD" ||sample_type[is]=="WJets" ||sample_type[is]=="ZJets" ||sample_type[is]=="TTbar" ||sample_type[is]=="SingleTop" ||sample_type[is]=="TTbarV" ||sample_type[is]=="VV/VVV" ||sample_type[is]=="Other" ||sample_type[is]=="mc"){
				if(histos.count(histonames[n]+hs)>0){
					histos[histonames[n]+hs] ->SetFillColor(colour);
					histos[histonames[n]+hs] ->SetLineColor(colour);
				} if(histos2.count(histonames[n]+hs)>0){
					histos2[histonames[n]+hs]->SetFillColor(colour);
					histos2[histonames[n]+hs]->SetLineColor(colour);
				} if(histos3.count(histonames[n]+hs)>0){
					histos3[histonames[n]+hs]->SetFillColor(colour);
					histos3[histonames[n]+hs]->SetLineColor(colour);
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
				if(sample_type[is]=="mc")   continue;
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
	for(map<string,TH3D*>::iterator h=histos3.begin(); h!=histos3.end();++h){
		h->second->Write();
	}
	for(map<string,TH2D*>::iterator h=histos2.begin(); h!=histos2.end();++h){
		h->second->Write();
	}
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Write();
	}
	for(map<string,THStack*>::iterator h=stacks.begin(); h!=stacks.end();++h){
		h->second->Write();
	}
	fsavefile->Close();
	cout << "Saved histograms in " << outputdir << outputname << endl;

	if(doplotting){
		cout << "plotting histograms ..." << endl;
		MakePlots(histos, leptontypesize, lepton_type, stacks, Legend1, histonames);
        	MakeCorrelationPlots(histos2, leptontypesize, lepton_type, histonames, outputdir);
	}

	delete Legend1;
	delete legend;
	delete fsavefile;

}//function

//this function calls Make1DPlotsRatio(...)
void MakePlots(map<string, TH1D*> histos, const int leptontypesize, string *lepton_type, map<string, THStack*> stacks, TLegend *Legend1, vector<string> histonames){

	cout << "Plotting 1D histograms" << endl;

	for(unsigned int n = 0; n<histonames.size(); ++n){//plot 1d
	   	for(int il = 0; il<leptontypesize; ++il){
			if(plotonlyLL && lepton_type[il]!="LL") continue;
			string name = histonames[n];
			if(!(saveD2) && ( name=="D2" || name=="D2_withMlb" || name=="D2_afterMT2llge85" || name=="D2_afterMT2bb_linMETge150" || name=="D2_afterMT2lbge200" || name=="D2_afterMT2lb_withMlbcutge180" || name=="D2_afterMT2lb_masslessge140" ) ) continue;

			string h = string("_") + lepton_type[il];
			string hs1 = string("_") + lepton_type[il] + string("_") + string("mc");
			string hs2 = string("_") + lepton_type[il] + string("_") + string("data");
			string hs3 = string("_") + lepton_type[il] + string("_") + string("Stop");
			if(histos.count(name+hs1)==0) continue;
			TString ytitle = "Events";//histos[name+hs1]->GetYaxis()->GetTitle();
			TString xtitle = histos[name+hs1]->GetXaxis()->GetTitle();
			TString outname = name + hs3 + (logflag ? "_log" : "") + "_overlay";

		if(histos[name+hs2]->Integral()>0 || histos[name+hs1]->Integral()>0 || histos[name+hs3]->Integral()>0){
			if(!plotonlywithratio) Make1DPlotsNoRatio(stacks[name+h], histos[name+hs1], histos[name+hs2], histos[name+hs3], logflag, false, outname, Legend1, xtitle, ytitle, 1.);
			Make1DPlotsRatio(stacks[name+h], histos[name+hs1], histos[name+hs2], histos[name+hs3], logflag, false, outname, Legend1, xtitle, ytitle, 1.);
		}

	}}

}

//makes the 1D plots including ratio - as copy from MassPlotter.cc function - no detailed comments
void Make1DPlotsRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, TH1D* histsusy, bool logFlag, bool normalize, TString outname, TLegend *leg, TString xtitle, TString ytitle, float overlayScale){

	// define canvas and pads 
	TH1D *h1 = (TH1D*)histmc->Clone("h1_copy");
	TH1D *h2 = (TH1D*)histdata->Clone("h2_copy");
	TH1D *h3 = (TH1D*)histsusy->Clone("h3_copy");

	h1->SetTitle("");
	h2->SetTitle("");
	h3->SetTitle("");
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
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2.5*max;
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
	h3->Scale(overlayScale ? overlayScale : h2->Integral() / h3->Integral());
	h3->SetFillColor(0);
	h3->SetLineStyle(kDotted);
	h3->Draw("samehist");

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.03);
	TString text ="";
	text = outname;

	TitleBox.DrawLatex(0.13,0.943,text.Data());

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
	if(fSave)Util::Print(c1, save, outputdir);

	

}

//make 1D plot without ratio - as copy from a MassPlotter.cc function, no detailed comments
void Make1DPlotsNoRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, TH1D* histsusy, bool logFlag, bool normalize, TString outname, TLegend *leg, TString xtitle, TString ytitle, float overlayScale){

	TCanvas *col = new TCanvas(outname, "", 0, 0, 600, 600);
	col->SetRightMargin(0.08);
	col->SetLeftMargin(0.18);
	col->SetTopMargin(0.07);
	col->SetBottomMargin(0.17);

	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logFlag) {
		gPad->SetLogy(1);
		hstack     -> SetMinimum(0.05);
		histmc     -> SetMinimum(0.05);
		histdata   -> SetMinimum(0.05);
		histsusy   -> SetMinimum(0.05);
	}else{
		hstack->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = histdata->GetMaximum();
	double max2 = histmc->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	histdata  ->SetMaximum(max);
	histmc    ->SetMaximum(max);
	hstack    ->SetMaximum(max);

	hstack->Draw("hist");
	if(histdata->Integral()>0) {
		histdata       ->Draw("sameE");
	}
	histsusy->Scale(overlayScale ? overlayScale : histdata->Integral() / histsusy->Integral());
	histsusy->SetLineStyle(kDotted);
	histsusy->SetFillColor(0);
	histsusy->Draw("samehist");
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
	TString text = outname;

	TitleBox.DrawLatex(0.18,0.943,text.Data());


	hstack->GetXaxis()->SetTitle(xtitle);
	hstack->GetXaxis()->SetLabelSize(0.05);
	hstack->GetXaxis()->SetTitleSize(0.05);
	hstack->GetXaxis()->SetTitleOffset(1.1);

	hstack->GetYaxis()->SetTitle(ytitle);
	hstack->GetYaxis()->SetLabelSize(0.05);
	hstack->GetYaxis()->SetTitleSize(0.05);
	hstack->GetYaxis()->SetTitleOffset(1.47);

	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, outname, outputdir);
	if(fSave)Util::PrintEPS(col, outname, outputdir);
}

//does 2D correlation plots
void MakeCorrelationPlots(map<string, TH2D*> histos2, const int leptontypesize, string *lepton_type, vector<string> histonames, TString outputdirectory){

	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	cout << "Do 2D Correlation plots" << endl;
	for(unsigned int n = 0; n<histonames.size(); ++n){//plot 1d
	   	for(int il = 0; il<leptontypesize; ++il){
			if(plotonlyLL && lepton_type[il]!="LL") continue;
			string name = histonames[n];
			if(!(saveD2) && ( name=="MT2lb_massless_vs_D2" || name=="MT2lb_vs_D2" || name=="MT2lb_withMlbcut_vs_D2" || name=="MT2ll_vs_D2" || name=="MT2bb_vs_D2" || name=="D2_vs_BTagDiscr" ) ) continue;

			string hs0 = name + string("_") + lepton_type[il] + string("_") + string("TTbar");
			string hs1 = name + string("_") + lepton_type[il] + string("_") + string("mc");
			string hs2 = name + string("_") + lepton_type[il] + string("_") + string("data");
			string hs3 = name + string("_") + lepton_type[il] + string("_") + string("Stop");
			string hs4 = name + string("_") + lepton_type[il] + string("_") + string("SingleTop");
			col->cd();
			if(logflag) gPad->SetLogz(1);
			if(histos2.count(hs0)>0) {
				histos2[hs0]->SetMinimum(0.005);
				col->SetName(hs0.c_str() );
				col->SetTitle(hs0.c_str() );
				histos2[hs0]->Draw("COLZ");
				col->Update();
				if(histos2[hs0]->GetEntries()>0) Util::PrintNoEPS(col, hs0, outputdirectory/*, false*/);
				if(histos2[hs0]->GetEntries()>0) Util::PrintEPS(col, hs0, outputdirectory);
				col->Clear();
			}
			if(histos2.count(hs1)>0) {
				histos2[hs1]->SetMinimum(0.005);
				col->SetName(hs1.c_str() );
				col->SetTitle(hs1.c_str() );
				histos2[hs1]->Draw("COLZ");
				col->Update();
				if(histos2[hs1]->GetEntries()>0) Util::PrintNoEPS(col, hs1, outputdirectory/*, false*/);
				if(histos2[hs1]->GetEntries()>0) Util::PrintEPS(col, hs1, outputdirectory);
				col->Clear();
			}
			if(histos2.count(hs2)>0) {
				histos2[hs2]->SetMinimum(0.5);
				col->SetName(hs2.c_str() );
				col->SetTitle(hs2.c_str() );
				histos2[hs2]->Draw("COLZ");
				col->Update();
				if(histos2[hs2]->GetEntries()>0) Util::PrintNoEPS(col, hs2, outputdirectory/*, false*/);
				if(histos2[hs2]->GetEntries()>0) Util::PrintEPS(col, hs2, outputdirectory);
				col->Clear();
			}
			if(histos2.count(hs3)>0) {
				histos2[hs3]->SetMinimum(5.*10e-5);
				col->SetName(hs3.c_str() );
				col->SetTitle(hs3.c_str() );
				histos2[hs3]->Draw("COLZ");
				col->Update();
				if(histos2[hs3]->GetEntries()>0) Util::PrintNoEPS(col, hs3, outputdirectory/*, false*/);
				if(histos2[hs3]->GetEntries()>0) Util::PrintEPS(col, hs3, outputdirectory);
				col->Clear();
			}
			if(plotsingletop && histos2.count(hs4)>0) {
				histos2[hs4]->SetMinimum(5.*10e-5);
				col->SetName(hs4.c_str() );
				col->SetTitle(hs4.c_str() );
				histos2[hs4]->Draw("COLZ");
				col->Update();
				if(histos2[hs4]->GetEntries()>0) Util::PrintNoEPS(col, hs4, outputdirectory/*, false*/);
				if(histos2[hs4]->GetEntries()>0) Util::PrintEPS(col, hs4, outputdirectory);
				col->Clear();
			}
		}
	}
	delete col;
}

//does 3D correlation plots
void Make3DPlots(map<string, TH3D*> histos3, const int leptontypesize, string *lepton_type, vector<string> histonames, TString outputdirectory){

	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	cout << "Do 3D Correlation plots" << endl;
	for(unsigned int n = 0; n<histonames.size(); ++n){//plot 1d
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
				if(histos3[hs0]->GetEntries()>0) Util::PrintNoEPS(col, hs0, outputdirectory/*, false*/);
				if(histos3[hs0]->GetEntries()>0) Util::PrintEPS(col, hs0, outputdirectory);
				col->Clear();
			}
			if(histos3.count(hs1)>0) {
				col->SetName(hs1.c_str() );
				col->SetTitle(hs1.c_str() );
				histos3[hs1]->Draw("BOX");
				col->Update();
				if(histos3[hs1]->GetEntries()>0) Util::PrintNoEPS(col, hs1, outputdirectory/*, false*/);
				if(histos3[hs1]->GetEntries()>0) Util::PrintEPS(col, hs1, outputdirectory);
				col->Clear();
			}
			if(histos3.count(hs2)>0) {
				col->SetName(hs2.c_str() );
				col->SetTitle(hs2.c_str() );
				histos3[hs2]->Draw("BOX");
				col->Update();
				if(histos3[hs2]->GetEntries()>0) Util::PrintNoEPS(col, hs2, outputdirectory/*, false*/);
				if(histos3[hs2]->GetEntries()>0) Util::PrintEPS(col, hs2, outputdirectory);
				col->Clear();
			}
			if(histos3.count(hs3)>0) {
				col->SetName(hs3.c_str() );
				col->SetTitle(hs3.c_str() );
				histos3[hs3]->Draw("BOX");
				col->Update();
				if(histos3[hs3]->GetEntries()>0) Util::PrintNoEPS(col, hs3, outputdirectory/*, false*/);
				if(histos3[hs3]->GetEntries()>0) Util::PrintEPS(col, hs3, outputdirectory);
				col->Clear();
			}
		}
	}
	delete col;
}

//standard function to read samples.dat
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
