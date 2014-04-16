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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/BEfficiency/BTagSFWeight.hh"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/Stops/SolveTTbar.h"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/Stops/SolveTTbar.c"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/Stops/SolveTTbarNew.h"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/Stops/SolveTTbarNew.C"

//run via root -l -b -q StopFeasibilityBTagVariable.C++
//as this code is very similar to StopFeasibility.C (only here you can have variable btag requirement (like define MT2(lb) without b-jet requirement))
//no detailed comments are made
//note that not all plots of StopFeasibility.C are done
//and some plots are done for cases ==1b, >=1b, ==2b, >=2b simulatneously

using namespace std;

void load(const char* filename = "samples_2141_dataonly.dat");
void StopFeasibilityBTagVariable();
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
const int fVerbose = 3; // already defined below
TString fPath;

TString btagging_file = "/shome/haweber/MT2Analysis_8TeV/Code/Efficiencies/BEfficiencies_DileptonicStop_CSVM.root";
TFile*  btagefffile;
TH1D*   hbeff;
TH1D*   hceff;
TH1D*   hleff;
TString taggerName("CSVM");// SSVHPT or SSVHEM-->Define this later via Tagger and discr!
//NOTE: In cutstream this has also to be corrected
int     njets  = -2;
int     ntags  = -2;//- means >=, while + means ==; -99 means no tagging requirement
int     Tagger = 4;//TCHE:0 TCHP:1 SSVHE:2 SSVHP:3, CSV: 4
//medium: 0.629, try also: tight: 0.898, loose: 0.244
float   discr  = 0.629;
float   bpt    = 40.;//add in cutStream, not used in cutstream
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

//use (for MT2(lb),etc) any jet if true, if false use only bjets, except for ==1b case: here use 1b and an additional jet with same pT
Bool_t  anyjet                    = false;

//as this code is very similar to StopFeasibility.C (only here you can have variable btag requirement (like define MT2(lb) without b-jet requirement))
//no detailed comments are made
//note that not all plots of StopFeasibility.C are done
//and some plots are done for cases ==1b, >=1b, ==2b, >=2b simulatneously
void StopFeasibilityBTagVariable(){

gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

  btagefffile	= TFile::Open(btagging_file);
  hbeff		= (TH1D*)btagefffile->Get("BEfficiency_mc");//also NJets40eq2
  hceff		= (TH1D*)btagefffile->Get("CEfficiency_mc");
  hleff		= (TH1D*)btagefffile->Get("LEfficiency_mc");
  float hptupperedge = hbeff->GetBinLowEdge(hbeff->GetNbinsX()) + hbeff->GetBinWidth(hbeff->GetNbinsX());//upper pt edge

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
    << "(NEles+NMuons)==2"                 << "&&"//veto secondary leptons //
    << "NBJets40CSVM>=1";

    if(jpt<40)      cutStream << "&& NJetsIDLoose  >=2";
    else if(jpt<50) cutStream << "&& NJetsIDLoose40>=2";
    else            cutStream << "&& NJetsIDLoose50>=2";

	TString cuts = cutStream.str().c_str();
	cuts = basecuts + (string)"&&" + cuts;

	std::ostringstream triggerStreamEE;
	triggerStreamEE << "( "
	<< "(trigger.HLT_DiElectrons==1)" << " )";//updateXXX
	TString triggerEE = triggerStreamEE.str().c_str();
	std::ostringstream triggerStreamEMu;
	triggerStreamEMu << "( "
	<< "(trigger.HLT_EMu==1)" << " )";//updateXXX
	TString triggerEMu = triggerStreamEMu.str().c_str();
	std::ostringstream triggerStreamMuMu;
	triggerStreamMuMu << "( "
	<< "(trigger.HLT_DiMuons==1)" << " )";//updateXXX
	TString triggerMuMu = triggerStreamMuMu.str().c_str();


    const int sampletypesize = 10;
    string sample_type[sampletypesize] = {"WJets", "ZJets", "TTbar", "SingleTop", "TTbarV", "VV/VVV", "Other", "mc", "Stop", "data"};

    const int leptontypesize = 5;
    string lepton_type[leptontypesize] = {"MuMu", "EMu", "EE", "SF", "LL"};

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
		//>=0b
		mapname = "Mll_0b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant dilepton mass, >=0b"; xtitle = "M_{ll} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 105, 0, 1050);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		// >= 1b
		mapname = "MT2lb_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb), >=1b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb), >=1b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlbcut_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut, >=1b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_withMlbcut_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2(lb) after M(lb) cut, >=1b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MlbAll_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant mass of leptons and b jets, >=1b"; xtitle = "M_{lb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant mass of leptons and b jets making M_{T2}(lb), >=1b"; xtitle = "M_{lb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mll_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant dilepton mass, >=1b"; xtitle = "M_{ll} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 105, 0, 1050);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l), >=1b"; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_withMlbcut_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l) after M(lb) cut, >=1b"; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Mass fit discriminant D_{2}, >=1b"; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_withMlb_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Mass fit discriminant after M_{lb} cut, >=1b"; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		// == 1b
		mapname = "MT2lb_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb), ==1b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb), ==1b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlbcut_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut, ==1b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_withMlbcut_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2(lb) after M(lb) cut, ==1b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MlbAll_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant mass of leptons and b jets, ==1b"; xtitle = "M_{lb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant mass of leptons and b jets making M_{T2}(lb), ==1b"; xtitle = "M_{lb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mll_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant dilepton mass, ==1b"; xtitle = "M_{ll} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 105, 0, 1050);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l), ==1b"; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_withMlbcut_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l) after M(lb) cut, ==1b"; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Mass fit discriminant D_{2}, ==1b"; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_withMlb_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Mass fit discriminant after M_{lb} cut, ==1b"; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		// >=2b
		mapname = "MT2lb_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb), >=2b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb), >=2b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlbcut_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut, >=2b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_withMlbcut_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2(lb) after M(lb) cut, >=2b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MlbAll_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant mass of leptons and b jets, >=2b"; xtitle = "M_{lb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant mass of leptons and b jets making M_{T2}(lb), >=2b"; xtitle = "M_{lb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mll_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant dilepton mass, >=2b"; xtitle = "M_{ll} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 105, 0, 1050);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l), >=2b"; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_withMlbcut_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l) after M(lb) cut, >=2b"; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Mass fit discriminant D_{2}, >=2b"; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_withMlb_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Mass fit discriminant after M_{lb} cut, >=2b"; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		// ==2b
		mapname = "MT2lb_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb), ==2b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb), ==2b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlbcut_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) after M_{lb} cut, ==2b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_withMlbcut_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless MT2(lb) after M(lb) cut, ==2b"; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MlbAll_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant mass of leptons and b jets, ==2b"; xtitle = "M_{lb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 90, 0, 900);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mlb_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant mass of leptons and b jets making M_{T2}(lb), ==2b"; xtitle = "M_{lb} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "Mll_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Invariant dilepton mass, ==2b"; xtitle = "M_{ll} [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 105, 0, 1050);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l), ==2b"; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2ll_withMlbcut_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(l) after M(lb) cut, ==2b"; xtitle = "M_{T2}(l) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Mass fit discriminant D_{2}, ==2b"; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "D2_withMlb_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "Mass fit discriminant after M_{lb} cut, ==2b"; xtitle = "D_{2}";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 30, 0, 300);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		//>=1b
		mapname = "MT2lb_massless_vs_MT2ll_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_withMlb_vs_MT2ll_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) with M_{lb} cut vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlb_vs_MT2ll_1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) with M_{lb} cut vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlbcut_afterMT2llge85_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_afterMT2llge85_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		//>=2b
		mapname = "MT2lb_massless_vs_MT2ll_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_withMlb_vs_MT2ll_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) with M_{lb} cut vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlb_vs_MT2ll_2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) with M_{lb} cut vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlbcut_afterMT2llge85_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_afterMT2llge85_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		//==1b
		mapname = "MT2lb_massless_vs_MT2ll_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_withMlb_vs_MT2ll_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) with M_{lb} cut vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlb_vs_MT2ll_eq1b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) with M_{lb} cut vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlbcut_afterMT2llge85_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_afterMT2llge85_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		// ==2b
		mapname = "MT2lb_massless_vs_MT2ll_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_massless_withMlb_vs_MT2ll_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "massless M_{T2}(lb) with M_{lb} cut vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlb_vs_MT2ll_eq2b"; if(av) vs.push_back(mapname); mapname = mapname + hs; title = "M_{T2}(lb) with M_{lb} cut vs M_{T2}(l)"; xtitle = "M_{T2}(lb) [GeV]"; ytitle = "M_{T2}(l) [GeV]";
        	if(histos2.count(mapname) == 0 ) histos2[mapname] = new TH2D(mapname.c_str(), title.c_str(), 70, 0, 700, 140, 0, 700);
		histos2[mapname]->GetXaxis()->SetTitle(xtitle.c_str()); histos2[mapname]->GetYaxis()->SetTitle(ytitle.c_str());
		mapname = "MT2lb_withMlbcut_afterMT2llge85_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_afterMT2llge85_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 70, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		//NEW1 03/08/2012
		// >=1b
		mapname = "MT2lb_withMlb_bothl_true_trueb_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_onel_true_trueb_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_nol_true_trueb_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_bothb_true_truel_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_oneb_true_truel_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_nob_true_truel_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_bothl_true_trueb_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_onel_true_trueb_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_nol_true_trueb_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_bothb_true_truel_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_oneb_true_truel_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_nob_true_truel_1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		// ==1b
		mapname = "MT2lb_withMlb_bothl_true_trueb_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_onel_true_trueb_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_nol_true_trueb_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_bothb_true_truel_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_oneb_true_truel_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_nob_true_truel_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_bothl_true_trueb_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_onel_true_trueb_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_nol_true_trueb_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_bothb_true_truel_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_oneb_true_truel_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_nob_true_truel_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		//>=2b
		mapname = "MT2lb_withMlb_bothl_true_trueb_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_onel_true_trueb_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_nol_true_trueb_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_bothb_true_truel_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_oneb_true_truel_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_nob_true_truel_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_bothl_true_trueb_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_onel_true_trueb_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_nol_true_trueb_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_bothb_true_truel_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_oneb_true_truel_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_nob_true_truel_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		//==2b
		mapname = "MT2lb_withMlb_bothl_true_trueb_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_onel_true_trueb_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_nol_true_trueb_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_bothb_true_truel_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_oneb_true_truel_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_nob_true_truel_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_bothl_true_trueb_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_onel_true_trueb_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_nol_true_trueb_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_bothb_true_truel_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_oneb_true_truel_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT2lb_withMlb_afterMT2llge85_nob_true_truel_eq2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T2}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		//new 08/10/2012
		mapname = "MT_lb_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT_lb_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT_lMET_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T}(l) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT_lMET_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T}(l) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT_bMET_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T}(b) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT_bMET_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T}(b) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT_lb_withMlb_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT_lb_withMlb_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T}(lb) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT_lMET_withMlb_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T}(l) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT_lMET_withMlb_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T}(l) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT_bMET_withMlb_2b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T}(b) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());
		mapname = "MT_bMET_withMlb_eq1b"; if(av) vs.push_back(mapname); title = mapname; mapname = mapname + hs; xtitle = "M_{T}(b) [GeV]";
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), title.c_str(), 140, 0, 700);
		histos[mapname]->GetXaxis()->SetTitle(xtitle.c_str());

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



    if(!(runAnalysis)) return;

    for(size_t i = 0; i < fSamples.size(); ++i){
        
		int samplecounter = 0;
	    int ssb = 0; int osb = 0;
	    int vetoeddiscrcount = 0; int discrcount = 0; int globsol = 0;
	    int statvec[3]; statvec[0]=0; statvec[1]=0; statvec[2]=0;
	    int nstptot = 0; int ndertot = 0; int nrottot = 0; int nscntot = 0; int globalsolution = 0;
   	    if(runData==false && fSamples[i].type=="data")  continue;
	    if(calcsusy==false && fSamples[i].type=="susy") continue;
	    if(debug && fSamples[i].type=="mc")             continue;
        
	    string sampletype = (string)fSamples[i].type;
	    string leptype = "LL";
	    if(sampletype==(string)"mc"){
		     if(fSamples[i].sname=="Wtolnu") sampletype = (string)"WJets";
		else if(fSamples[i].sname=="DY")     sampletype = (string)"ZJets";
		else if(fSamples[i].name=="TTbar")   sampletype = (string)"TTbar";
		else if(fSamples[i].sname=="TTbarV") sampletype = (string)"TTbarV";
		else if(fSamples[i].sname=="Top")    sampletype = (string)"SingleTop";//no ttbar
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
            ++samplecounter;
            if ( fVerbose>2 && counter % 50000 == 0 )  cout << "+++ Proccessing event " << counter << endl;

            Double_t weight = sample_weight;
            if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 
            Bool_t recoedee   = false;// exact 2 ele, 0 muo
            Bool_t recoedemu  = false;// exact 1 ele, 1 muo
            Bool_t recoedmumu = false;// exact 0 ele, 2 muo
		//change this - might have some inefficiencies - three leptons with one not passing IDMedium==1
		//20-20 selection --> third lepton is 10 GeV 
            if(fMT2tree->NEles>=2 &&fMT2tree->ele[0].lv.Pt()>20&&fMT2tree->ele[1].lv.Pt()>20&&fMT2tree->ele[2].lv.Pt()<20&&(fMT2tree->NMuons==0||fMT2tree->muo[0].lv.Pt()<10) && fMT2tree->ele[0].IDMedium==1 && fMT2tree->ele[1].IDMedium==1) recoedee   = true;
            if(fMT2tree->NMuons>=2&&fMT2tree->muo[0].lv.Pt()>20&&fMT2tree->muo[1].lv.Pt()>20&&fMT2tree->muo[2].lv.Pt()<10&&(fMT2tree->NEles ==0||fMT2tree->ele[0].lv.Pt()<10)) recoedmumu = true;
            if(fMT2tree->NMuons>=1&&fMT2tree->muo[0].lv.Pt()>20&&fMT2tree->muo[1].lv.Pt()<10&&fMT2tree->NEles>=1 &&fMT2tree->ele[0].lv.Pt()>20&&fMT2tree->ele[1].lv.Pt()<10 && fMT2tree->ele[0].IDMedium==1) recoedemu  = true;//XXX change emu to 20/20 selection, left ee,mumu at 20/10 for now
	    //id cuts?

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


	    if(recoedee   && ((fMT2tree->ele[0].lv + fMT2tree->ele[1].lv).M() )<10) continue;
	    if(recoedemu  && ((fMT2tree->ele[0].lv + fMT2tree->muo[0].lv).M() )<10) continue;
	    if(recoedmumu && ((fMT2tree->muo[0].lv + fMT2tree->muo[1].lv).M() )<10) continue;

	    if(!(recoedee)   && !(recoedemu)   && !(recoedmumu)  ) continue;//require dilepton
	    if(!(recoedosee) && !(recoedosemu) && !(recoedosmumu) ) continue;//require os-dilepton, maybe comment this for background est.
	    //if(!(recoedeenZ) ||                   !(recoedmumunZ) ) continue;//require of Zpeak

	    if(recoedee)   leptype = "EE";
	    if(recoedemu)  leptype = "EMu";
	    if(recoedmumu) leptype = "MuMu";

		//trigger weights: EEweight; EMuweight; MUMUweight, if different lumis;
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


	    unsigned int NumBJetsReal = fMT2tree->NBJets40CSVM;
	    if(NumBJetsReal<1) continue;//require at least one b-jet 

		//three different weights, depending on b-jet (1,2,>=2)
	    TString tagger = taggerName;//default
	    float SFweightErr = 0;
	    float SFweight = 1;//outside, since need it there later
	    float SFweight1 = 1;//outside, since need it there later
	    float SFweight2 = 1;//outside, since need it there later
	    if(!(dofastbtagSFreweighting)){
	    if((!fMT2tree->misc.isData)){
		vector<float> jetEff; jetEff.clear();
		vector<float> jetEffErr; jetEffErr.clear();
		vector<float> jetSF; jetSF.clear();
		vector<float> jetSFErr; jetSFErr.clear();
		int njetsusuable = 0;
		for(int n = 0; n<fMT2tree->NJets; ++n){
			if(fMT2tree->jet[n].isPFIDLoose==false) continue;
			if(fMT2tree->jet[n].Flavour<=-7777)     continue;
			float effPt  = fMT2tree->jet[n].lv.Pt();
			float effEta = fabs(fMT2tree->jet[n].lv.Eta());
			if(effPt<20.)  continue;
			if(effEta>2.4) continue;
			++njetsusuable;
			if(effPt>=hptupperedge) effPt = hptupperedge - 1.;//so now it lies in last bin
			if(abs(fMT2tree->jet[n].Flavour)==5){
				jetEff.push_back( float(hbeff->GetBinContent(hbeff->FindBin(effPt))) );
				jetEffErr.push_back( float(hbeff->GetBinError(hbeff->FindBin(effPt))) );
				float SFErr;
				float SF = getBtagSF(SFErr, tagger, effPt, effEta);
				float SFFSErr = 0;
				if(fMT2tree->misc.isFastSim){
					float SFFS = FastSimCorrectionFactor(SFFSErr, tagger, 5, effPt, effEta);
					SF = SFFS*SF;
					SFErr = SFFS*SFFS;//is that correct? FSerr is separate
				}
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr );//first implementation



			}
			else if(abs(fMT2tree->jet[n].Flavour)==4){
				jetEff.push_back( float(hceff->GetBinContent(hceff->FindBin(effPt))) );
				jetEffErr.push_back( float(hceff->GetBinError(hceff->FindBin(effPt))) );
				float SFErr;
				float SF = getBtagSF(SFErr, tagger, effPt, effEta);
				if(fMT2tree->misc.isFastSim){
					float SFFS = FastSimCorrectionFactor(SFFSErr, tagger, 4, effPt, effEta);
					SF = SFFS*SF;
					SFErr = SFErr*SFFS;//is that correct? FSerr is separate
				}
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr*2.);
			}
			else {
					jetEff.push_back( float(hleff->GetBinContent(hleff->FindBin(effPt))) );
					jetEffErr.push_back( float(hleff->GetBinError(hleff->FindBin(effPt))) );
				float SFErr;
				float SFmin = getMistagSF(tagger, effPt, effEta, -1 );//second implementation
				float SFmax = getMistagSF(tagger, effPt, effEta, 1 );//second implementation
				float SF = getMistagSF(tagger, effPt, effEta, 0 );

				float SFFSErr = 0;
				if(fMT2tree->misc.isFastSim){
					float SFFS = FastSimCorrectionFactor(SFFSErr, tagger, 3, effPt, effEta);
					SFErr= SF-SFmin;
					SFErr = sqrt(SFErr*SFErr*SFFS*SFFS + SFFSErr*SFFSErr*SF*SF);//is SFFSerr only relative
					SFErr = SFmax-SF;
					SFErr = sqrt(SFErr*SFErr*SFFS*SFFS + SFFSErr*SFFSErr*SF*SF);//is SFFSerr only relative
					SF = SFFS*SF;
				}
				if(SF-SFmin > SFmax-SF) SFErr = SF-SFmin;
				else SFErr = SFmax-SF;
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr );

			}
		}
		SFweight     = getBTagEventWeightError(SFweightErr, jetEff, jetEffErr, jetSF/*, jetSFErr*/, ntags);
		SFweight1    = getBTagEventWeightError(SFweightErr, jetEff, jetEffErr, jetSF/*, jetSFErr*/, -1);
		SFweight2    = getBTagEventWeightError(SFweightErr, jetEff, jetEffErr, jetSF/*, jetSFErr*/, -2);
		if(SFweight==0){
			if(njetsusuable!=0){
				cout << "Event has zero weight, do not use it" << endl;
				continue;
			}
			else { //event has no flavour information, use average event weight
				SFweight = pow(0.95,abs(ntags));
				if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
			}
		}
	    }
            }
            else{
	    if((!fMT2tree->misc.isData)){
	     //from ttbar payload
             if(tagger=="SSVHPT"){
		SFweight = pow(0.95,abs(ntags));
		if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
		SFweight1    = 0.95;
		SFweight2    = 0.95*0.95;
	     }
             if(tagger=="SSVHEM"){
		SFweight = pow(0.96,abs(ntags));
		if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
		SFweight1    = 0.96;
		SFweight2    = 0.96*0.96;
	     }
             if(tagger=="CSVM"){
		SFweight = pow(0.972,abs(ntags));
		if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
		SFweight1    = 0.972;
		SFweight2    = 0.972*0.972;
	     }
             if(tagger=="CSVT"){
		SFweight = pow(0.947,abs(ntags));
		if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
		SFweight1    = 0.947;
		SFweight2    = 0.947*0.947;
	     }
             if(tagger=="CSVL"){
		SFweight = pow(1.02,abs(ntags));
		if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
		SFweight1    = 1.02;
		SFweight2    = 1.02*1.02;
	     }
             if(tagger=="JPM"){
		SFweight = pow(0.95,abs(ntags));
		if(ntags==(-99)) {SFweight = 1; SFweightErr = 0; }
		SFweight1    = 0.95;
		SFweight2    = 0.95*0.95;
	     }
	    }
	    }

	   TLorentzVector l1, l2, met;
	   vector<TLorentzVector> bjets; bjets.clear();
	   vector<int> bjetsflavour; bjetsflavour.clear();
	   vector<float> bjetsdiscriminant; bjetsdiscriminant.clear();
	   vector<int> bjetistrueb; bjetistrueb.clear();//use bool to int
	   int l1c(-999), l2c(-999);
	   Bool_t l1e(true), l2e(true);
	   vector<int> jind; jind.clear();
	   vector<int> bind; bind.clear();
	   if(recoedee){ 
		l1  = (fMT2tree->ele[0]).lv;     l2  = (fMT2tree->ele[1]).lv;
		l1c = (fMT2tree->ele[0]).Charge; l2c = fMT2tree->ele[1].Charge;
		l1e = true;                      l2e = true;}
	   else if(recoedemu){
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
	   for(int n = 0; n<fMT2tree->NJets; ++n){//loop only over bjets
		if(fMT2tree->jet[n].isPFIDLoose==false) continue;
		if(fMT2tree->jet[n].lv.Pt()<jpt) continue;
		if(fabs(fMT2tree->jet[n].lv.Eta())>jeta) continue;
		float btempdiscr = (Tagger==3 ? fMT2tree->jet[n].bTagProbSSVHP : Tagger==2 ? fMT2tree->jet[n].bTagProbSSVHE : Tagger==4 ? fMT2tree->jet[n].bTagProbCSV : Tagger==5 ? fMT2tree->jet[n].bTagProbJProb :fMT2tree->jet[n].bTagProbCSV);//default CSV
		jind.push_back(n);//NOTE: add here an additional pt constraint?????
		bjets.push_back(fMT2tree->jet[n].lv); bjetsflavour.push_back(fMT2tree->jet[n].Flavour); bjetsdiscriminant.push_back(btempdiscr); bind.push_back(n);
		if(btempdiscr>=discr && fMT2tree->jet[n].lv.Pt()>=bpt && fabs(fMT2tree->jet[n].lv.Eta())<=beta) bjetistrueb.push_back(true);
		else bjetistrueb.push_back(false);
	   }
	   //sort 
	   for(size_t n =0; n<bjets.size();++n){
		for(size_t m =n+1; m<bjets.size();++m){
			if(bjetsdiscriminant[m]<=bjetsdiscriminant[n]) continue;
			swap(bjetsdiscriminant[n],bjetsdiscriminant[m]);
			swap(bjets[n],bjets[m]);
			swap(bjetsflavour[n],bjetsflavour[m]);
			swap(bind[n],bind[m]);
			swap(bjetistrueb[n],bjetistrueb[m]);
		}
	  }
	   for(size_t n =0; n<bjets.size();++n){
		for(size_t m =n+1; m<bjets.size();++m){
			if(bjetsdiscriminant[m]==bjetsdiscriminant[n]){
				if(bjets[m].Pt()<=bjets[n].Pt()) continue;
				swap(bjetsdiscriminant[n],bjetsdiscriminant[m]);
				swap(bjets[n],bjets[m]);
				swap(bjetsflavour[n],bjetsflavour[m]);
				swap(bind[n],bind[m]);
				swap(bjetistrueb[n],bjetistrueb[m]);
			}
		}
	  }
	   unsigned int NumBJets = bjets.size();

	   met = fMT2tree->pfmet[0];

	//lb with b vector
	   float MT2lbV(-999.99), MT2lbV_massless(-999.99);
	   float MT2lbV_withMlb(-999.99), MT2lbV_massless_withMlb(-999.99);//as this might be difficult to get
	   float MT2_l1b[NumBJets][NumBJets];//l1 with first b[], l2 with second b[];
	   float MT2_massless_l1b[NumBJets][NumBJets];//l1 with first b[], l2 with second b[];
	   float Ml1b[NumBJets];
	   float Ml2b[NumBJets];
	   float MT_l1b[NumBJets];
	   float MT_l2b[NumBJets];
	   float Mbbarr[NumBJets][NumBJets];
	   float MT2bbarr[NumBJets][NumBJets];
	   float MT2bb_linMET[NumBJets][NumBJets];
	   float MT2bb_mW_linMET[NumBJets][NumBJets];
	   float MT2bbmin(-999.99), MT2bbmin_mW_linMET(-999.99), MT2bbmin_linMET(-999.99);
	   int ind_l1bxMT2lb(-1), ind_l2bxMT2lb(-1), ind_l1bxMT2lbMlb(-1), ind_l2bxMT2lbMlb(-1);
	   int ind_b1x_MT2bbmin(-1), ind_b2x_MT2bbmin(-1), ind_b1x_MT2bbmin_mW_linMET(-1), ind_b2x_MT2bbmin_mW_linMET(-1), ind_b1x_MT2bbmin_linMET(-1), ind_b2x_MT2bbmin_linMET(-1), ind_l1bxMT2lbMassless(-1), ind_l2bxMT2lbMassless(-1);
	   float MT2ll(-999.99), Mll(-999.99);//check that this is not defined above anymore
	   float MT_l1(-999.99), MT_l2(-999.99); float MT_b[NumBJets]; 
	   bool  truel1bxtop[NumBJets], truel2bxtop[NumBJets], truebtop[NumBJets];//use truebtop as truebfromtop --> do genfitting
	   int  truel1flav(-99), truel2flav(-99);
	   int match1(-1),match2(-1); int matchb[NumBJets];
	   double dR1(999.), dR2(999.); double dRb[NumBJets];
	   //reset all arrays
	   for(unsigned int i1 = 0; i1<NumBJets; ++i1){
		Ml1b[i1] = -999.99; Ml2b[i1] = -999.99;
		MT_l1b[i1] = -999.99; MT_l2b[i1] = -999.99; MT_b[i1] = -999.99;
		dRb[i1] = 999.; matchb[i1] = -999.99;
		truel1bxtop[i1] = false; truel2bxtop[i1] = false; truebtop[i1] = false;
		for(unsigned int i2 = 0; i2<NumBJets; ++i2){
			MT2_l1b[i1][i2]  = -999.99; MT2_massless_l1b[i1][i2] = -999.99;
			MT2bbarr[i1][i2] = -999.99; MT2bb_mW_linMET[i1][i2]  = -999.99; MT2bb_linMET[i1][i2] = -999.99; 
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
	//	MT_l1  = (l1+met).Mt();
	//	MT_l2  = (l2+met).Mt();
		MT_l1  = fMT2tree->GetMT(l1,l1.M(),met,0.);
		MT_l2  = fMT2tree->GetMT(l2,l2.M(),met,0.);
		Ml1b[i1]    = (l1+bjets[i1]).M();
		Ml2b[i1]    = (l2+bjets[i1]).M();
		MT_l1b[i1]  = fMT2tree->GetMT(l1,l1.M(),bjets[i1],bjets[i1].M());
		MT_l2b[i1]  = fMT2tree->GetMT(l2,l2.M(),bjets[i1],bjets[i1].M());
		MT_b[i1]    = fMT2tree->GetMT(met,0.,bjets[i1],bjets[i1].M());
	   }
	   for(unsigned int i1 = 0; i1<NumBJets; ++i1){//need to do this like this since we need Mlb
		for(unsigned int i2 = 0; i2<NumBJets; ++i2){
			if(i2==i1) continue;
			MT2_l1b[i1][i2]              = fMT2tree->CalcMT2(0.,   true, l1+bjets[i1], l2+bjets[i2], met);
	   		MT2_massless_l1b[i1][i2]     = fMT2tree->CalcMT2(0.,  false, l1+bjets[i1], l2+bjets[i2], met);
			MT2bbarr[i1][i2]             = fMT2tree->CalcMT2(0.,   true, bjets[i1],    bjets[i2],    met);
			MT2bb_mW_linMET[i1][i2]      = fMT2tree->CalcMT2(80.4, true, bjets[i1],    bjets[i2],    met+l1+l2);
			MT2bb_linMET[i1][i2]         = fMT2tree->CalcMT2(0.,   true, bjets[i1],    bjets[i2],    met+l1+l2);
			Mbbarr[i1][i2]               = (bjets[i1]+bjets[i2]).M();

		}
	   }
	   if(anyjet){
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
           }//if(anyjet)
	   else {//use only bjets
		if(NumBJetsReal>=2){//can use only bjets
			for(unsigned int i1 = 0; i1<NumBJets; ++i1){
				if(bjetsdiscriminant[i1]<discr) continue;
				if(bjets[i1].Pt()       <bpt  ) continue;
				if(fabs(bjets[i1].Eta())>beta ) continue;
				if(!(bjetistrueb[i1]) ) cout << "Error line " << __LINE__ << ": This should be a true b" << endl;
				for(unsigned int i2 = 0; i2<NumBJets; ++i2){
					if(i2==i1) continue;
					if(bjetsdiscriminant[i2]<discr) continue;
					if(bjets[i2].Pt()       <bpt  ) continue;
					if(fabs(bjets[i2].Eta())>beta ) continue;
					//additional checks
					if(bjetsdiscriminant[i1]==bjetsdiscriminant[i2] && bjets[i1].Pt()==bjets[i2].Pt() && bjets[i1].Eta()==bjets[i2].Eta()) continue;
					if(!(bjetistrueb[i2]) ) cout << "Error line " << __LINE__ << ": This should be a true b" << endl;
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
		}//if(NumBJetsReal>=2)
		else if(NumBJetsReal==1){
			//first jet (i.e. combined with l1) is b
			for(unsigned int i1 = 0; i1<NumBJets; ++i1){//only one jet should pass this
				if(bjetsdiscriminant[i1]<discr) continue;
				if(bjets[i1].Pt()       <bpt  ) continue;
				if(fabs(bjets[i1].Eta())>beta ) continue;
				if(!(bjetistrueb[i1]) ) cout << "Error line " << __LINE__ << ": This should be a true b" << endl;
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
			//second jet (i.e. combined with l2) is b
			for(unsigned int i1 = 0; i1<NumBJets; ++i1){
				for(unsigned int i2 = 0; i2<NumBJets; ++i2){//only one jet should pass this
					if(i2==i1) continue;
					if(bjetsdiscriminant[i2]<discr) continue;
					if(bjets[i2].Pt()       <bpt  ) continue;
					if(fabs(bjets[i2].Eta())>beta ) continue;
					if(!(bjetistrueb[i2]) ) cout << "Error line " << __LINE__ << ": This should be a true b" << endl;
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
		}//if(NumBJetsReal==1)
	   }//if not anyjet, i.e. only bjets

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
		if(MT2lbV_withMlb>=0){
			for(int n = 0; n<25; ++n){
				if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
				int genlepID  = fMT2tree->genlept[n].ID;
				int genlepMID = fMT2tree->genlept[n].MID;
				//int genlepGMID = fMT2tree->genlept[n].GMID;
				if(abs(genlepID)==5){
					if(abs(genlepMID)!=6)     continue;
					if(genlepMID*genlepID<=0) continue;
					float dRb1 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l1bxMT2lbMlb]);
					if(dRb1<0.5 && dRb1<dR_mb1){ matchedgenp_tob1 = n; dR_mb1 = dRb1; }
					float dRb2 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l2bxMT2lbMlb]);
					if(dRb2<0.5 && dRb2<dR_mb2){ matchedgenp_tob2 = n; dR_mb2 = dRb1; }
				}
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
		if(MT2lbV_withMlb>=0){
			for(int n = 0; n<25; ++n){
				if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
				int genlepID = fMT2tree->genlept[n].ID;
				int genlepMID = fMT2tree->genlept[n].MID;
				//int genlepGMID = fMT2tree->genlept[n].GMID;
				if(abs(genlepID)==5){
					if(abs(genlepMID)==2000006) cout << "Cool generated t_2" << endl;;
					if(abs(genlepMID)!=1000006) continue;
					if(genlepMID*genlepID<=0)   continue;
					float dRb1 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l1bxMT2lbMlb]);
					if(dRb1<0.5 && dRb1<dR_mb1){ matchedgenp_tob1 = n; dR_mb1 = dRb1; }
					float dRb2 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l2bxMT2lbMlb]);
					if(dRb2<0.5 && dRb2<dR_mb2){ matchedgenp_tob2 = n; dR_mb2 = dRb1; }
				}
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
				if(abs(genlepMID) !=1000024 && abs(genlepMID)!=15 && abs(genlepMID)!=1000011 && abs(genlepMID)!=1000013)      continue;
				if(abs(genlepGMID)!=1000006 && abs(genlepMID)==1000024) continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepMID)==1000013) continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepMID)==1000011) continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepGMID)!=1000015 && abs(genlepMID)==15 )          continue;
				if(abs(genlepMID) ==1000024 && genlepMID *genlepID>=0 ) continue;
				if(abs(genlepMID) ==1000013 && genlepMID *genlepID<=0 ) continue;
				if(abs(genlepMID) ==1000011 && genlepMID *genlepID<=0 ) continue;
				if(abs(genlepMID) ==15      && genlepMID *genlepID<=0 ) continue;
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
				if(abs(genlepMID) !=1000024 && abs(genlepMID)!=15 && abs(genlepMID)!=1000011 && abs(genlepMID)!=1000013)      continue;
				if(abs(genlepGMID)!=1000006 && abs(genlepMID)==1000024) continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepMID)==1000013) continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepMID)==1000011) continue;
				if(abs(genlepGMID)!=1000024 && abs(genlepGMID)!=1000015 && abs(genlepMID)==15 )          continue;
				if(abs(genlepMID) ==1000024 && genlepMID *genlepID>=0 ) continue;
				if(abs(genlepMID) ==1000013 && genlepMID *genlepID<=0 ) continue;
				if(abs(genlepMID) ==1000011 && genlepMID *genlepID<=0 ) continue;
				if(abs(genlepMID) ==15      && genlepMID *genlepID<=0 ) continue;
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
		if(MT2lbV_withMlb>=0){
			for(int n = 0; n<25; ++n){
				if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
				int genlepID  = fMT2tree->genlept[n].ID;
				int genlepMID = fMT2tree->genlept[n].MID;
				//int genlepGMID = fMT2tree->genlept[n].GMID;
				if(abs(genlepID)==5){
					if(abs(genlepMID)==2000006) cout << "Cool generated t_2" << endl;;
					if(abs(genlepMID)!=1000006) continue;
					if(genlepMID*genlepID<=0)   continue;
					float dRb1 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l1bxMT2lbMlb]);
					if(dRb1<0.5 && dRb1<dR_mb1){ matchedgenp_tob1 = n; dR_mb1 = dRb1; }
					float dRb2 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l2bxMT2lbMlb]);
					if(dRb2<0.5 && dRb2<dR_mb2){ matchedgenp_tob2 = n; dR_mb2 = dRb1; }
				}
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
				if(abs(genlepGMID)!=1000006 && abs(genlepGMID)!=1000024 && abs(genlepGMID)!=1000015) continue;
				if(abs(genlepGMID)==1000006 && abs(genlepMID) !=1000024) continue;
				if(abs(genlepGMID)==1000015 && abs(genlepMID) !=15     ) continue;
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
				if(abs(genlepGMID)!=1000006 && abs(genlepGMID)!=1000024 && abs(genlepGMID)!=1000015) continue;
				if(abs(genlepGMID)==1000006 && abs(genlepMID) !=1000024) continue;
				if(abs(genlepGMID)==1000015 && abs(genlepMID) !=15     ) continue;
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
		if(MT2lbV_withMlb>=0){
			for(int n = 0; n<25; ++n){
				if(fMT2tree->genlept[n].lv.Pt()<=0.01) continue;//veto 0
				int genlepID = fMT2tree->genlept[n].ID;
				int genlepMID = fMT2tree->genlept[n].MID;
				//int genlepGMID = fMT2tree->genlept[n].GMID;
				if(abs(genlepID)==5){
					if(abs(genlepMID)==2000006) cout << "Cool generated t_2" << endl;;
					if(abs(genlepMID)!=1000006) continue;
					if(genlepMID*genlepID<=0)   continue;
					float dRb1 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l1bxMT2lbMlb]);
					if(dRb1<0.5 && dRb1<dR_mb1){ matchedgenp_tob1 = n; dR_mb1 = dRb1; }
					float dRb2 = fMT2tree->genlept[n].lv.DeltaR(bjets[ind_l2bxMT2lbMlb]);
					if(dRb2<0.5 && dRb2<dR_mb2){ matchedgenp_tob2 = n; dR_mb2 = dRb1; }
				}
			}
		}
	   }//stopTobChargino, too generic, might not pick up everything in order not to keep fakes
	}

	//new discriminant implementation
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
	bool failedroot(false);
	if(calcD2&&calcD2Minuit){//NOTE: here also Luc's implementation needed, for 4th order equation solving
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

	if(calcD2&&calcD2Luc){
	    	statvec1[0]=0; statvec1[1]=0; statvec1[2]=0; int nstptot1 = 0; int ndertot1 = 0; int nrottot1 = 0; int nscntot1 = 0;
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
		if(failedroot ) cout << "Luc: D2 " << discriminant1 << " q1z/q2z " <<q1z1<<"/"<<q2z1<<" solution/conv/loop/fail " << globalsolution1<<"/"<<statvec1[0]<<"/"<<statvec1[1]<<"/"<<statvec1[2]<<endl<<endl;
			statvec[0] += statvec1[0]; statvec[1] += statvec1[1]; statvec[2] += statvec1[2];
			nstptot += nstptot1; ndertot += ndertot1; nrottot += nrottot1; nscntot += nscntot1;
			globalsolution += globalsolution1;

	}

	double D2 = discriminant1;//stupid, but for simpleneww use this

	string hs = string("_") + leptype + string("_") + sampletype;
	histos["Mll_0b" + hs]->Fill(Mll, weight);//NOTE: might change this depending on jetthreshold
	if(NumBJetsReal>=1){
		if(SFweight1!=0) weight = weight * SFweight1;
		histos["Mll_1b" + hs]->Fill(Mll, weight);//NOTE: might change this depending on jetthreshold
		if(SFweight1!=0) weight = weight / SFweight1;
	}
	if(NumBJetsReal>=2){
		if(SFweight2!=0) weight = weight * SFweight2;
		histos["Mll_2b" + hs]->Fill(Mll, weight);//NOTE: might change this depending on jetthreshold
		if(SFweight2!=0) weight = weight / SFweight2;
	}
	if(NumBJetsReal==1){
		if(SFweight1!=0) weight = weight * SFweight1;
		histos["Mll_eq1b" + hs]->Fill(Mll, weight);//NOTE: might change this depending on jetthreshold
		if(SFweight1!=0) weight = weight / SFweight1;
	}
	if(NumBJetsReal==2){
		if(SFweight2!=0) weight = weight * SFweight2;
		histos["Mll_eq2b" + hs]->Fill(Mll, weight);//NOTE: might change this depending on jetthreshold
		if(SFweight2!=0) weight = weight / SFweight2;
	}
	if(!(recoedeenZ) &&                   !(recoedmumunZ) && !(recoedemu) ) continue;//require of Zpeakveto

	bool truel1notb1(false), truel2notb2(false);//need to be defined outside
	//now real plots
	for(unsigned int i1 = 0; i1<NumBJets; ++ i1){
	for(unsigned int i2 = 0; i2<NumBJets; ++ i2){
		if(i1==i2) continue;
		if((int)i1!=ind_l1bxMT2lb && truel1bxtop[i1]) truel1notb1 = true;
		if((int)i2!=ind_l2bxMT2lb && truel2bxtop[i2]) truel2notb2 = true;
	}}

	if(NumBJetsReal>=1){
		if(SFweight1!=0) weight = weight * SFweight1;
		histos["MT2lb_1b"          + hs]->Fill(MT2lbV,              weight);
		histos["MT2lb_massless_1b" + hs]->Fill(MT2lbV_massless,     weight);
		if(MT2lbV_withMlb>=0 )          histos["MT2lb_withMlbcut_1b"          + hs]->Fill(MT2lbV_withMlb,          weight);
		if(MT2lbV_massless_withMlb>=0 ) histos["MT2lb_massless_withMlbcut_1b" + hs]->Fill(MT2lbV_massless_withMlb, weight);
		for(unsigned int n = 0; n<NumBJets; ++n){
			histos["MlbAll_1b" + hs]->Fill(Ml1b[n],             weight);
			histos["MlbAll_1b" + hs]->Fill(Ml2b[n],             weight);
		}
		histos["Mlb_1b"            + hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
		histos["Mlb_1b"            + hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);
		histos["MT2ll_1b"          + hs]->Fill(MT2ll,               weight);
		if(Ml1b[ind_l1bxMT2lb]<180. && Ml2b[ind_l2bxMT2lb]<180.)                histos["MT2ll_withMlbcut_1b" + hs]->Fill(MT2ll, weight);
		if(D2!=9999999.)                                                        histos["D2_1b"               + hs]->Fill(D2,    weight);
		if(Ml1b[ind_l1bxMT2lb]<180. && Ml2b[ind_l2bxMT2lb]<180. &&D2!=9999999.) histos["D2_withMlb_1b"       + hs]->Fill(D2,    weight);

		if(MT2lbV_withMlb>=0 && MT2ll>=85.)          histos["MT2lb_withMlbcut_afterMT2llge85_1b"          + hs]->Fill(MT2lbV_withMlb,          weight);
		if(MT2lbV_massless_withMlb>=0 && MT2ll>=85.) histos["MT2lb_massless_withMlbcut_afterMT2llge85_1b" + hs]->Fill(MT2lbV_massless_withMlb, weight);
		if(MT2lbV_massless>=0 && MT2ll>=85.)         histos["MT2lb_massless_afterMT2llge85_1b"            + hs]->Fill(MT2lbV_massless,         weight);

		if(MT2lbV_withMlb>=0 && MT2ll>=0)           histos2["MT2lb_withMlb_vs_MT2ll_1b"                   + hs]->Fill(MT2lbV_withMlb,          MT2ll, weight);
		if(MT2lbV_massless_withMlb>=0 && MT2ll>=0)  histos2["MT2lb_massless_withMlb_vs_MT2ll_1b"          + hs]->Fill(MT2lbV_massless_withMlb, MT2ll, weight);
		if(MT2lbV_massless>=0 && MT2ll>=0)          histos2["MT2lb_massless_vs_MT2ll_1b"                  + hs]->Fill(MT2lbV_massless,         MT2ll, weight);

		if(matchedgenp_tob1>=0 && matchedgenp_tob2>=0 && matchedgenp_tob1!=matchedgenp_tob2){//not match to same b
			if(matchedgenl_tol1>=0 && matchedgenl_tol2>=0 && matchedgenl_tol1!=matchedgenl_tol2) {//not match to same b
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_bothl_true_trueb_1b"                + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_bothl_true_trueb_1b" + hs]->Fill(MT2lbV_withMlb, weight);
			} else if((matchedgenl_tol1==-1 && matchedgenl_tol2>=0)||(matchedgenl_tol1>=0 && matchedgenl_tol2==-1)) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_onel_true_trueb_1b"                 + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_onel_true_trueb_1b"  + hs]->Fill(MT2lbV_withMlb, weight);
			} else if(matchedgenl_tol1==-1 && matchedgenl_tol2==-1) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_nol_true_trueb_1b"                  + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_nol_true_trueb_1b"   + hs]->Fill(MT2lbV_withMlb, weight);
			}
		}
		if(matchedgenl_tol1>=0 && matchedgenl_tol2>=0 && matchedgenl_tol1!=matchedgenl_tol2){//not match to same b
			if(matchedgenp_tob1>=0 && matchedgenp_tob2>=0 && matchedgenp_tob1!=matchedgenp_tob2) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_bothb_true_truel_1b"                + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_bothb_true_truel_1b" + hs]->Fill(MT2lbV_withMlb, weight);
			} else if((matchedgenp_tob1==-1 && matchedgenp_tob2>=0)||(matchedgenp_tob1>=0 && matchedgenp_tob2==-1)) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_oneb_true_truel_1b"                 + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_oneb_true_truel_1b"  + hs]->Fill(MT2lbV_withMlb, weight);
			} else if(matchedgenp_tob1==-1 && matchedgenp_tob2==-1) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_nob_true_truel_1b"                  + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_nob_true_truel_1b"   + hs]->Fill(MT2lbV_withMlb, weight);
			}
		}

		if(SFweight1!=0) weight = weight / SFweight1;
	}
	if(NumBJetsReal>=2){
		if(SFweight2!=0) weight = weight * SFweight2;
		histos["MT2lb_2b"          + hs]->Fill(MT2lbV,              weight);
		histos["MT2lb_massless_2b" + hs]->Fill(MT2lbV_massless,     weight);
		if(MT2lbV_withMlb>=0 )          histos["MT2lb_withMlbcut_2b"          + hs]->Fill(MT2lbV_withMlb,          weight);
		if(MT2lbV_massless_withMlb>=0 ) histos["MT2lb_massless_withMlbcut_2b" + hs]->Fill(MT2lbV_massless_withMlb, weight);
		for(unsigned int n = 0; n<NumBJets; ++n){
			histos["MlbAll_2b" + hs]->Fill(Ml1b[n],             weight);
			histos["MlbAll_2b" + hs]->Fill(Ml2b[n],             weight);
		}
		histos["Mlb_2b"            + hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
		histos["Mlb_2b"            + hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);
		histos["MT2ll_2b"          + hs]->Fill(MT2ll,               weight);
		if(Ml1b[ind_l1bxMT2lb]<180. && Ml2b[ind_l2bxMT2lb]<180.)                histos["MT2ll_withMlbcut_2b" + hs]->Fill(MT2ll, weight);
		if(D2!=9999999.)                                                        histos["D2_2b"               + hs]->Fill(D2,    weight);
		if(Ml1b[ind_l1bxMT2lb]<180. && Ml2b[ind_l2bxMT2lb]<180. &&D2!=9999999.) histos["D2_withMlb_2b"       + hs]->Fill(D2,    weight);

		if(MT2lbV_withMlb>=0 && MT2ll>=85.)          histos["MT2lb_withMlbcut_afterMT2llge85_2b"          + hs]->Fill(MT2lbV_withMlb,          weight);
		if(MT2lbV_massless_withMlb>=0 && MT2ll>=85.) histos["MT2lb_massless_withMlbcut_afterMT2llge85_2b" + hs]->Fill(MT2lbV_massless_withMlb, weight);
		if(MT2lbV_massless>=0 && MT2ll>=85.)         histos["MT2lb_massless_afterMT2llge85_2b"            + hs]->Fill(MT2lbV_massless,         weight);

		if(MT2lbV_withMlb>=0 && MT2ll>=0)           histos2["MT2lb_withMlb_vs_MT2ll_2b"                   + hs]->Fill(MT2lbV_withMlb,          MT2ll, weight);
		if(MT2lbV_massless_withMlb>=0 && MT2ll>=0)  histos2["MT2lb_massless_withMlb_vs_MT2ll_2b"          + hs]->Fill(MT2lbV_massless_withMlb, MT2ll, weight);
		if(MT2lbV_massless>=0 && MT2ll>=0)          histos2["MT2lb_massless_vs_MT2ll_2b"                  + hs]->Fill(MT2lbV_massless,         MT2ll, weight);

		if(matchedgenp_tob1>=0 && matchedgenp_tob2>=0 && matchedgenp_tob1!=matchedgenp_tob2){//not match to same b
			if(matchedgenl_tol1>=0 && matchedgenl_tol2>=0 && matchedgenl_tol1!=matchedgenl_tol2) {//not match to same b
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_bothl_true_trueb_2b"                + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_bothl_true_trueb_2b" + hs]->Fill(MT2lbV_withMlb, weight);
			} else if((matchedgenl_tol1==-1 && matchedgenl_tol2>=0)||(matchedgenl_tol1>=0 && matchedgenl_tol2==-1)) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_onel_true_trueb_2b"                 + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_onel_true_trueb_2b"  + hs]->Fill(MT2lbV_withMlb, weight);
			} else if(matchedgenl_tol1==-1 && matchedgenl_tol2==-1) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_nol_true_trueb_2b"                  + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_nol_true_trueb_2b"   + hs]->Fill(MT2lbV_withMlb, weight);
			}
		}
		if(matchedgenl_tol1>=0 && matchedgenl_tol2>=0 && matchedgenl_tol1!=matchedgenl_tol2){//not match to same b
			if(matchedgenp_tob1>=0 && matchedgenp_tob2>=0 && matchedgenp_tob1!=matchedgenp_tob2) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_bothb_true_truel_2b"                + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_bothb_true_truel_2b" + hs]->Fill(MT2lbV_withMlb, weight);
			} else if((matchedgenp_tob1==-1 && matchedgenp_tob2>=0)||(matchedgenp_tob1>=0 && matchedgenp_tob2==-1)) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_oneb_true_truel_2b"                 + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_oneb_true_truel_2b"  + hs]->Fill(MT2lbV_withMlb, weight);
			} else if(matchedgenp_tob1==-1 && matchedgenp_tob2==-1) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_nob_true_truel_2b"                  + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_nob_true_truel_2b"   + hs]->Fill(MT2lbV_withMlb, weight);
			}
		}
		//new 08/10/2012
		if(MT2lbV>=0){
			histos["MT_lb_2b"           + hs]->Fill(MT_l1b[ind_l1bxMT2lb],    weight);
			histos["MT_lb_2b"           + hs]->Fill(MT_l2b[ind_l2bxMT2lb],    weight);
			histos["MT_bMET_2b"         + hs]->Fill(MT_b[ind_l1bxMT2lb],      weight);
			histos["MT_bMET_2b"         + hs]->Fill(MT_b[ind_l2bxMT2lb],      weight);
			histos["MT_lMET_2b"         + hs]->Fill(MT_l1,                    weight);
			histos["MT_lMET_2b"         + hs]->Fill(MT_l2,                    weight);
		}
		if(MT2lbV_withMlb>=0){
			histos["MT_lb_withMlb_2b"   + hs]->Fill(MT_l1b[ind_l1bxMT2lbMlb], weight);
			histos["MT_lb_withMlb_2b"   + hs]->Fill(MT_l2b[ind_l2bxMT2lbMlb], weight);
			histos["MT_bMET_withMlb_2b" + hs]->Fill(MT_b[ind_l1bxMT2lbMlb],   weight);
			histos["MT_bMET_withMlb_2b" + hs]->Fill(MT_b[ind_l2bxMT2lbMlb],   weight);
			histos["MT_lMET_withMlb_2b" + hs]->Fill(MT_l1,                    weight);
			histos["MT_lMET_withMlb_2b" + hs]->Fill(MT_l2,                    weight);
		}
		if(SFweight2!=0) weight = weight / SFweight2;

	}

	if(NumBJetsReal==1){
		if(SFweight1!=0) weight = weight * SFweight1;
		histos["MT2lb_eq1b"          + hs]->Fill(MT2lbV,              weight);
		histos["MT2lb_massless_eq1b" + hs]->Fill(MT2lbV_massless,     weight);
		if(MT2lbV_withMlb>=0 )          histos["MT2lb_withMlbcut_eq1b"          + hs]->Fill(MT2lbV_withMlb,          weight);
		if(MT2lbV_massless_withMlb>=0 ) histos["MT2lb_massless_withMlbcut_eq1b" + hs]->Fill(MT2lbV_massless_withMlb, weight);
		for(unsigned int n = 0; n<NumBJets; ++n){
			histos["MlbAll_eq1b" + hs]->Fill(Ml1b[n],             weight);
			histos["MlbAll_eq1b" + hs]->Fill(Ml2b[n],             weight);
		}
		histos["Mlb_eq1b"            + hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
		histos["Mlb_eq1b"            + hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);
		histos["MT2ll_eq1b"          + hs]->Fill(MT2ll,               weight);
		if(Ml1b[ind_l1bxMT2lb]<180. && Ml2b[ind_l2bxMT2lb]<180.)                histos["MT2ll_withMlbcut_eq1b" + hs]->Fill(MT2ll, weight);
		if(D2!=9999999.)                                                        histos["D2_eq1b"               + hs]->Fill(D2,    weight);
		if(Ml1b[ind_l1bxMT2lb]<180. && Ml2b[ind_l2bxMT2lb]<180. &&D2!=9999999.) histos["D2_withMlb_eq1b"       + hs]->Fill(D2,    weight);

		if(MT2lbV_withMlb>=0 && MT2ll>=85.)          histos["MT2lb_withMlbcut_afterMT2llge85_eq1b"          + hs]->Fill(MT2lbV_withMlb,          weight);
		if(MT2lbV_massless_withMlb>=0 && MT2ll>=85.) histos["MT2lb_massless_withMlbcut_afterMT2llge85_eq1b" + hs]->Fill(MT2lbV_massless_withMlb, weight);
		if(MT2lbV_massless>=0 && MT2ll>=85.)         histos["MT2lb_massless_afterMT2llge85_eq1b"            + hs]->Fill(MT2lbV_massless,         weight);

		if(MT2lbV_withMlb>=0 && MT2ll>=0)           histos2["MT2lb_withMlb_vs_MT2ll_eq1b"                   + hs]->Fill(MT2lbV_withMlb,          MT2ll, weight);
		if(MT2lbV_massless_withMlb>=0 && MT2ll>=0)  histos2["MT2lb_massless_withMlb_vs_MT2ll_eq1b"          + hs]->Fill(MT2lbV_massless_withMlb, MT2ll, weight);
		if(MT2lbV_massless>=0 && MT2ll>=0)          histos2["MT2lb_massless_vs_MT2ll_eq1b"                  + hs]->Fill(MT2lbV_massless,         MT2ll, weight);

		if(matchedgenp_tob1>=0 && matchedgenp_tob2>=0 && matchedgenp_tob1!=matchedgenp_tob2){//not match to same b
			if(matchedgenl_tol1>=0 && matchedgenl_tol2>=0 && matchedgenl_tol1!=matchedgenl_tol2) {//not match to same b
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_bothl_true_trueb_eq1b"                + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_bothl_true_trueb_eq1b" + hs]->Fill(MT2lbV_withMlb, weight);
			} else if((matchedgenl_tol1==-1 && matchedgenl_tol2>=0)||(matchedgenl_tol1>=0 && matchedgenl_tol2==-1)) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_onel_true_trueb_eq1b"                 + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_onel_true_trueb_eq1b"  + hs]->Fill(MT2lbV_withMlb, weight);
			} else if(matchedgenl_tol1==-1 && matchedgenl_tol2==-1) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_nol_true_trueb_eq1b"                  + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_nol_true_trueb_eq1b"   + hs]->Fill(MT2lbV_withMlb, weight);
			}
		}
		if(matchedgenl_tol1>=0 && matchedgenl_tol2>=0 && matchedgenl_tol1!=matchedgenl_tol2){//not match to same b
			if(matchedgenp_tob1>=0 && matchedgenp_tob2>=0 && matchedgenp_tob1!=matchedgenp_tob2) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_bothb_true_truel_eq1b"                + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_bothb_true_truel_eq1b" + hs]->Fill(MT2lbV_withMlb, weight);
			} else if((matchedgenp_tob1==-1 && matchedgenp_tob2>=0)||(matchedgenp_tob1>=0 && matchedgenp_tob2==-1)) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_oneb_true_truel_eq1b"                 + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_oneb_true_truel_eq1b"  + hs]->Fill(MT2lbV_withMlb, weight);
			} else if(matchedgenp_tob1==-1 && matchedgenp_tob2==-1) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_nob_true_truel_eq1b"                  + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_nob_true_truel_eq1b"   + hs]->Fill(MT2lbV_withMlb, weight);
			}
		}
		//new 08/10/2012
		if(MT2lbV>=0){
			histos["MT_lb_eq1b"           + hs]->Fill(MT_l1b[ind_l1bxMT2lb],    weight);
			histos["MT_lb_eq1b"           + hs]->Fill(MT_l2b[ind_l2bxMT2lb],    weight);
			histos["MT_bMET_eq1b"         + hs]->Fill(MT_b[ind_l1bxMT2lb],      weight);
			histos["MT_bMET_eq1b"         + hs]->Fill(MT_b[ind_l2bxMT2lb],      weight);
			histos["MT_lMET_eq1b"         + hs]->Fill(MT_l1,                    weight);
			histos["MT_lMET_eq1b"         + hs]->Fill(MT_l2,                    weight);
		}
		if(MT2lbV_withMlb>=0){
			histos["MT_lb_withMlb_eq1b"   + hs]->Fill(MT_l1b[ind_l1bxMT2lbMlb], weight);
			histos["MT_lb_withMlb_eq1b"   + hs]->Fill(MT_l2b[ind_l2bxMT2lbMlb], weight);
			histos["MT_bMET_withMlb_eq1b" + hs]->Fill(MT_b[ind_l1bxMT2lbMlb],   weight);
			histos["MT_bMET_withMlb_eq1b" + hs]->Fill(MT_b[ind_l2bxMT2lbMlb],   weight);
			histos["MT_lMET_withMlb_eq1b" + hs]->Fill(MT_l1,                    weight);
			histos["MT_lMET_withMlb_eq1b" + hs]->Fill(MT_l2,                    weight);
		}
		if(SFweight1!=0) weight = weight / SFweight1;
	}
	if(NumBJetsReal==2){
		if(SFweight2!=0) weight = weight * SFweight2;
		histos["MT2lb_eq2b"          + hs]->Fill(MT2lbV,              weight);
		histos["MT2lb_massless_eq2b" + hs]->Fill(MT2lbV_massless,     weight);
		if(MT2lbV_withMlb>=0 )          histos["MT2lb_withMlbcut_eq2b"          + hs]->Fill(MT2lbV_withMlb,          weight);
		if(MT2lbV_massless_withMlb>=0 ) histos["MT2lb_massless_withMlbcut_eq2b" + hs]->Fill(MT2lbV_massless_withMlb, weight);
		for(unsigned int n = 0; n<NumBJets; ++n){
			histos["MlbAll_eq2b" + hs]->Fill(Ml1b[n],             weight);
			histos["MlbAll_eq2b" + hs]->Fill(Ml2b[n],             weight);
		}
		histos["Mlb_eq2b"            + hs]->Fill(Ml1b[ind_l1bxMT2lb], weight);
		histos["Mlb_eq2b"            + hs]->Fill(Ml2b[ind_l2bxMT2lb], weight);
		histos["MT2ll_eq2b"          + hs]->Fill(MT2ll,               weight);
		if(Ml1b[ind_l1bxMT2lb]<180. && Ml2b[ind_l2bxMT2lb]<180.)                histos["MT2ll_withMlbcut_eq2b" + hs]->Fill(MT2ll, weight);
		if(D2!=9999999.)                                                        histos["D2_eq2b"               + hs]->Fill(D2,    weight);
		if(Ml1b[ind_l1bxMT2lb]<180. && Ml2b[ind_l2bxMT2lb]<180. &&D2!=9999999.) histos["D2_withMlb_eq2b"       + hs]->Fill(D2,    weight);

		if(MT2lbV_withMlb>=0 && MT2ll>=85.)          histos["MT2lb_withMlbcut_afterMT2llge85_eq2b"          + hs]->Fill(MT2lbV_withMlb,          weight);
		if(MT2lbV_massless_withMlb>=0 && MT2ll>=85.) histos["MT2lb_massless_withMlbcut_afterMT2llge85_eq2b" + hs]->Fill(MT2lbV_massless_withMlb, weight);
		if(MT2lbV_massless>=0 && MT2ll>=85.)         histos["MT2lb_massless_afterMT2llge85_eq2b"            + hs]->Fill(MT2lbV_massless,         weight);

		if(MT2lbV_withMlb>=0 && MT2ll>=0)           histos2["MT2lb_withMlb_vs_MT2ll_eq2b"                   + hs]->Fill(MT2lbV_withMlb,          MT2ll, weight);
		if(MT2lbV_massless_withMlb>=0 && MT2ll>=0)  histos2["MT2lb_massless_withMlb_vs_MT2ll_eq2b"          + hs]->Fill(MT2lbV_massless_withMlb, MT2ll, weight);
		if(MT2lbV_massless>=0 && MT2ll>=0)          histos2["MT2lb_massless_vs_MT2ll_eq2b"                  + hs]->Fill(MT2lbV_massless,         MT2ll, weight);

		if(matchedgenp_tob1>=0 && matchedgenp_tob2>=0 && matchedgenp_tob1!=matchedgenp_tob2){//not match to same b
			if(matchedgenl_tol1>=0 && matchedgenl_tol2>=0 && matchedgenl_tol1!=matchedgenl_tol2) {//not match to same b
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_bothl_true_trueb_eq2b"                + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_bothl_true_trueb_eq2b" + hs]->Fill(MT2lbV_withMlb, weight);
			} else if((matchedgenl_tol1==-1 && matchedgenl_tol2>=0)||(matchedgenl_tol1>=0 && matchedgenl_tol2==-1)) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_onel_true_trueb_eq2b"                 + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_onel_true_trueb_eq2b"  + hs]->Fill(MT2lbV_withMlb, weight);
			} else if(matchedgenl_tol1==-1 && matchedgenl_tol2==-1) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_nol_true_trueb_eq2b"                  + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_nol_true_trueb_eq2b"   + hs]->Fill(MT2lbV_withMlb, weight);
			}
		}
		if(matchedgenl_tol1>=0 && matchedgenl_tol2>=0 && matchedgenl_tol1!=matchedgenl_tol2){//not match to same b
			if(matchedgenp_tob1>=0 && matchedgenp_tob2>=0 && matchedgenp_tob1!=matchedgenp_tob2) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_bothb_true_truel_eq2b"                + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_bothb_true_truel_eq2b" + hs]->Fill(MT2lbV_withMlb, weight);
			} else if((matchedgenp_tob1==-1 && matchedgenp_tob2>=0)||(matchedgenp_tob1>=0 && matchedgenp_tob2==-1)) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_oneb_true_truel_eq2b"                 + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_oneb_true_truel_eq2b"  + hs]->Fill(MT2lbV_withMlb, weight);
			} else if(matchedgenp_tob1==-1 && matchedgenp_tob2==-1) {
				if(MT2lbV_withMlb>=0              ) histos["MT2lb_withMlb_nob_true_truel_eq2b"                  + hs]->Fill(MT2lbV_withMlb, weight);
				if(MT2lbV_withMlb>=0 && MT2ll>=85.) histos["MT2lb_withMlb_afterMT2llge85_nob_true_truel_eq2b"   + hs]->Fill(MT2lbV_withMlb, weight);
			}
		}
		if(SFweight2!=0) weight = weight / SFweight2;

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
		if(histos.count(mapname+hs1)>0)  histos[mapname+h] ->Add(histos[mapname+hs1],1);
		if(histos.count(mapname+hs2)>0)  histos[mapname+h] ->Add(histos[mapname+hs2],1);
		if(histos.count(mapname+hs3)>0)  histos[mapname+h] ->Add(histos[mapname+hs3],1);
		if(histos2.count(mapname+hs1)>0) histos2[mapname+h]->Add(histos2[mapname+hs1],1);
		if(histos2.count(mapname+hs2)>0) histos2[mapname+h]->Add(histos2[mapname+hs2],1);
		if(histos2.count(mapname+hs3)>0) histos2[mapname+h]->Add(histos2[mapname+hs3],1);
		if(histos3.count(mapname+hs1)>0) histos3[mapname+h]->Add(histos3[mapname+hs1],1);
		if(histos3.count(mapname+hs2)>0) histos3[mapname+h]->Add(histos3[mapname+hs2],1);
		if(histos3.count(mapname+hs3)>0) histos3[mapname+h]->Add(histos3[mapname+hs3],1);
		string hsf = string("_") + string("SF") + string("_") + sample_type[is];
		if(histos.count(mapname+hs1)>0)  histos[mapname+hsf] ->Add(histos[mapname+hs1],1);
		if(histos.count(mapname+hs3)>0)  histos[mapname+hsf] ->Add(histos[mapname+hs3],1);
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
			Color_t colour = 603;//set default others if we have some failure
			     if(sample_type[is]=="QCD")       colour = 880;
			else if(sample_type[is]=="WJets")     colour = 417;
			else if(sample_type[is]=="ZJets")     colour = 419;
			else if(sample_type[is]=="TTbar")     colour = 401;
			else if(sample_type[is]=="SingleTop") colour = 595;
			else if(sample_type[is]=="TTbarV")    colour = 65;
			else if(sample_type[is]=="VV/VVV")    colour = 635;
			else if(sample_type[is]=="Other")     colour = 603;
			else if(sample_type[is]=="mc")        colour = 603;
			if(sample_type[is]!="Stop" && sample_type[is]!="data"){
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
				if(sample_type[is]=="mc")  continue;
				string hs = string("_") + lepton_type[il] + string("_") + sample_type[is];
				if(histos.count(histonames[n]+hs)==0) continue;//use only 1d for HStack
				stacks[histonames[n]+h] ->Add(histos[histonames[n]+hs]);
			}
		}
		if(leggy){
	   	for(int is = 0; is<sampletypesize; ++is){
			string hs = string("_") + string("LL") + string("_") + sample_type[is];
			if(sample_type[is]=="mc")             continue;
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
		Make3DPlots(histos3, leptontypesize, lepton_type, histonames, outputdir);
	}

	delete Legend1;
	delete legend;
	delete fsavefile;

}//function

void MakePlots(map<string, TH1D*> histos, const int leptontypesize, string *lepton_type, map<string, THStack*> stacks, TLegend *Legend1, vector<string> histonames){

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
			TString outname = name + hs3 + (flip_order ? "_flipped" : "") + (logflag ? "_log" : "") + (composited ? "_comp" : "") + "_overlay";

		if(histos[name+hs2]->Integral()>0 || histos[name+hs1]->Integral()>0 || histos[name+hs3]->Integral()>0){
			if(!plotonlywithratio) Make1DPlotsNoRatio(stacks[name+h], histos[name+hs1], histos[name+hs2], histos[name+hs3], logflag, false, outname, Legend1, xtitle, ytitle, 1.);
			Make1DPlotsRatio(stacks[name+h], histos[name+hs1], histos[name+hs2], histos[name+hs3], logflag, false, outname, Legend1, xtitle, ytitle, 1.);
		}

	}}

}

void Make1DPlotsRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, TH1D* histsusy, bool logFlag, bool normalize, TString outname, TLegend *leg, TString xtitle, TString ytitle, float overlayScale){

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
	TString text = outname;	TitleBox.DrawLatex(0.18,0.943,text.Data());


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
