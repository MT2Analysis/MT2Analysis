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
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string>
#include <cmath>
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"
//#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MassPlotter.hh"//only used if fMakeEfficienciesPresel==true;
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"
//#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/Utilities.hh"

using namespace std;

//void load(const char* filename = "/shome/haweber/CMSSW_4_2_3/src/DiLeptonAnalysis/NTupleProducer/MT2Analysis/Code/MT2AnalysisCode/MT2Code/samples/datasamples/samples_2141_dataonly.dat");
//inline vector<string> initializeHistos(map<string, TH1D*> histos, map<string, TH2D*> histos2, map<string, THStack*> hstacks, TLegend *legend, const int sampletypesize, string *sample_type, const int leptontypesize, string *lepton_type);
void PlotStopFeasibility();
//void MakeControlPlots(map<string, TH1D*> histos, map<string, TH2D*> histos2, map<string, THStack*> hstacks, TString basecuts, TString triggerEE, TString triggerEMu, TString triggerMuMu);
//void MakePlots(MassPlotter *tA, map<string, TH1D*> histos, const int sampletypesize, string *sample_type, const int leptontypesize, string *lepton_type, map<string, THStack*> stacks, TLegend *Legend1, vector<string> histonames);
void MakePlots(map<string, TH1D*> histos, const int sampletypesize, string *sample_type, const int leptontypesize, string *lepton_type, map<string, THStack*> stacks, TLegend *Legend1, vector<string> histonames);
void MakeCorrelationPlots(map<string, TH1D*> histos, map<string, TH2D*> histos2, const int sampletypesize, string *sample_type, const int leptontypesize, string *lepton_type, vector<string> histonames, TString outputdir);
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

//TString btagging_file          = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/Efficiencies/BEffHistos_PTbinned_allHT_SSVHPT.root";
//TFile* btagefffile;
//TString tagger = "SSVHEM";//or SSVHEM
//TH1D*  hbeff;
//TH1D*  hceff;
//TH1D*  hleff;
//int njets = -2;
//int ntags = -1;//- means >=, while + means ==; -99 means no tagging requirement
//int Tagger = 2;
//float discr = 1.74;
//float pt = 40.;

TString samples                   = "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/samples/samples_Stops_presel.dat";//dummy//updateXXX
//Bool_t  runData                   = true; // if you do only MCClosure you can set this to false
//Bool_t  calcsusy                  = true; // not used - can set false
//TString outputdir                 = "Stops/testeq2bge2jStop400bChi_WLSP_20120720_2_withFilters/noBRequirement40/addendum";// /nolog";//default
//TString outputdirold              = "Stops/testeq2bge2jStop400bChi_WLSP_20120720_2_withFilters/noBRequirement40/addendum";//default
TString outputdir                 = "Stops/test53X_corr_Stop400bChi300_WLSP/D2";// /nolog";//default
TString outputdirold              = "Stops/test53X_corr_Stop400bChi300_WLSP";//default
//TString outputdir                 = "Stops/testeq2bge2jStop400bChi_WLSP_20120816/Mll/nolog";// /nolog";//default
//TString outputdirold              = "Stops/testeq2bge2jStop400bChi_WLSP_20120816/Mll";//default
TString outputname                = "DiLeptonicStops.root";//default
//TString outputname                = "WithoutD2_DiLeptonicStops.root";//need to set plotD2 = false;
Bool_t  setYlog                   = true;
Bool_t  makectrlplots             = false;//set to true if want to do control plots
//Bool_t  runAnalysis               = true;//set to false if only want to do control plots
//Bool_t  dofastbtagSFreweighting   = true;
Bool_t  logflag                   = true;
Bool_t  flip_order                = false;//dummy
Bool_t  composited                = false;//dummy
Bool_t  fSave                     = true;
Bool_t  debug                     = false;//don't run over mc to save time
Bool_t  plotD2                    = true;
Bool_t  plotonlywithratio         = true;
Bool_t  plotonlyLL                = false;
Bool_t  plotonlyOF                = false;
Bool_t  plotonlySF                = true;
Bool_t  plotsingletop             = false;//for 2d plots only


void PlotStopFeasibility(){

gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
/*
  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
  tA->SetSave(true);
  tA->setVerbose(fVerbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
*/
//  btagefffile	= TFile::Open(btagging_file);
//  hbeff		= (TH1D*)btagefffile->Get("h_beff");
//  hceff		= (TH1D*)btagefffile->Get("h_ceff");
//  hleff		= (TH1D*)btagefffile->Get("h_leff");

    Util::MakeOutputDir(outputdir);
    map<string, TH1D*> histos;
    map<string, TH2D*> histos2;
    map<string, TH3D*> histos3;
    map<string, THStack*> stacks;
    map<string, THStack*> hstacks;
    TLegend *Legend1 = new TLegend(.71,.54,.91,.92);
    vector<string> histonames; histonames.clear();
    //Hemisphere *hemi_l1b1_l2b2;
    //Hemisphere *hemi_l1b2_l2b1;
    Legend1 -> SetFillColor(0);
    Legend1 -> SetBorderSize(0);

    TFile *oldfile = TFile::Open(outputdirold + outputname);


    const int sampletypesize = 10;//11;
    string sample_type[sampletypesize] = {/*"QCD", */"WJets", "ZJets", "TTbar", "SingleTop", "TTbarV", "VV/VVV", "Other", "mc", "Stop", "data"};

    const int leptontypesize = 5;//4;
    string lepton_type[leptontypesize] = {"MuMu", "EMu", "EE", "SF", "LL"};

	vector<string> vs; vs.clear();
	vector<string> vs1; vs1.clear();
	vector<string> vs2; vs2.clear();
	vector<string> vs3; vs3.clear();
	bool av = true;//append vector
		string mapname;

/*
		//from StopFeasibilityMll.C
		mapname = "Mll_0b_all"; if(av) vs.push_back(mapname);
		mapname = "Mll_0b_zoom"; if(av) vs.push_back(mapname);
		mapname = "Mll_0b_all_NJetsle2"; if(av) vs.push_back(mapname);
		mapname = "Mll_0b_zoom_NJetsle2"; if(av) vs.push_back(mapname);
		mapname = "Mll_0b_all_NJetsgt2"; if(av) vs.push_back(mapname);
		mapname = "Mll_0b_zoom_NJetsgt2"; if(av) vs.push_back(mapname);
		mapname = "Mll_0b_all_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_0b_zoom_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_0b_all_NJetsle2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_0b_zoom_NJetsle2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_0b_all_NJetsgt2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_0b_zoom_NJetsgt2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_1b_all"; if(av) vs.push_back(mapname);
		mapname = "Mll_1b_zoom"; if(av) vs.push_back(mapname);
		mapname = "Mll_1b_all_NJetsle2"; if(av) vs.push_back(mapname);
		mapname = "Mll_1b_zoom_NJetsle2"; if(av) vs.push_back(mapname);
		mapname = "Mll_1b_all_NJetsgt2"; if(av) vs.push_back(mapname);
		mapname = "Mll_1b_zoom_NJetsgt2"; if(av) vs.push_back(mapname);
		mapname = "Mll_1b_all_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_1b_zoom_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_1b_all_NJetsle2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_1b_zoom_NJetsle2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_1b_all_NJetsgt2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_1b_zoom_NJetsgt2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_2b_all"; if(av) vs.push_back(mapname);
		mapname = "Mll_2b_zoom"; if(av) vs.push_back(mapname);
		mapname = "Mll_2b_all_NJetsle2"; if(av) vs.push_back(mapname);
		mapname = "Mll_2b_zoom_NJetsle2"; if(av) vs.push_back(mapname);
		mapname = "Mll_2b_all_NJetsgt2"; if(av) vs.push_back(mapname);
		mapname = "Mll_2b_zoom_NJetsgt2"; if(av) vs.push_back(mapname);
		mapname = "Mll_2b_all_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_2b_zoom_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_2b_all_NJetsle2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_2b_zoom_NJetsle2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_2b_all_NJetsgt2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_2b_zoom_NJetsgt2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq0b_all"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq0b_zoom"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq0b_all_NJetsle2"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq0b_zoom_NJetsle2"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq0b_all_NJetsgt2"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq0b_zoom_NJetsgt2"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq0b_all_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq0b_zoom_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq0b_all_NJetsle2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq0b_zoom_NJetsle2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq0b_all_NJetsgt2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq0b_zoom_NJetsgt2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq1b_all"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq1b_zoom"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq1b_all_NJetsle2"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq1b_zoom_NJetsle2"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq1b_all_NJetsgt2"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq1b_zoom_NJetsgt2"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq1b_all_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq1b_zoom_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq1b_all_NJetsle2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq1b_zoom_NJetsle2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq1b_all_NJetsgt2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq1b_zoom_NJetsgt2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq2b_all"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq2b_zoom"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq2b_all_NJetsle2"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq2b_zoom_NJetsle2"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq2b_all_NJetsgt2"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq2b_zoom_NJetsgt2"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq2b_all_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq2b_zoom_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq2b_all_NJetsle2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq2b_zoom_NJetsle2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq2b_all_NJetsgt2_MET100"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq2b_zoom_NJetsgt2_MET100"; if(av) vs.push_back(mapname);
		mapname = "MT2l_0b_NJetsgt2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_0b_NJetsgt2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_0b_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_0b_NJetsle2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_0b_NJetsle2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_0b_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_0b_NJetsgt2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_0b_NJetsgt2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_0b_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_0b_NJetsle2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_0b_NJetsle2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_0b_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_1b_NJetsgt2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_1b_NJetsgt2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_1b_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_1b_NJetsle2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_1b_NJetsle2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_1b_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_1b_NJetsgt2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_1b_NJetsgt2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_1b_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_1b_NJetsle2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_1b_NJetsle2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_1b_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_2b_NJetsgt2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_2b_NJetsgt2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_2b_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_2b_NJetsle2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_2b_NJetsle2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_2b_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_2b_NJetsgt2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_2b_NJetsgt2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_2b_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_2b_NJetsle2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_2b_NJetsle2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_2b_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq0b_NJetsgt2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq0b_NJetsgt2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq0b_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq0b_NJetsle2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq0b_NJetsle2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq0b_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq0b_NJetsgt2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq0b_NJetsgt2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq0b_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq0b_NJetsle2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq0b_NJetsle2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq0b_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq1b_NJetsgt2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq1b_NJetsgt2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq1b_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq1b_NJetsle2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq1b_NJetsle2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq1b_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq1b_NJetsgt2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq1b_NJetsgt2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq1b_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq1b_NJetsle2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq1b_NJetsle2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq1b_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq2b_NJetsgt2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq2b_NJetsgt2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq2b_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq2b_NJetsle2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq2b_NJetsle2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2l_eq2b_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq2b_NJetsgt2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq2b_NJetsgt2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq2b_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq2b_NJetsle2_MET100_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq2b_NJetsle2_Mll20_70"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq2b_Mll20_70"; if(av) vs.push_back(mapname);
*/
/*
		//from StopFeasibilityBTagVariable.C --> add everything from untitled + eq1, eq2
		mapname = "MT2lb_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlbcut_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_withMlbcut_1b"; if(av) vs.push_back(mapname);
		mapname = "MlbAll_1b"; if(av) vs.push_back(mapname);
		mapname = "Mlb_1b"; if(av) vs.push_back(mapname);
		mapname = "Mll_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2ll_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2ll_withMlbcut_1b"; if(av) vs.push_back(mapname);
		if(plotD2){ mapname = "D2_1b"; if(av) vs.push_back(mapname); }
		if(plotD2){ mapname = "D2_withMlb_1b"; if(av) vs.push_back(mapname); }
		mapname = "MT2lb_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlbcut_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_withMlbcut_2b"; if(av) vs.push_back(mapname);
		mapname = "MlbAll_2b"; if(av) vs.push_back(mapname);
		mapname = "Mlb_2b"; if(av) vs.push_back(mapname);
		mapname = "Mll_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2ll_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2ll_withMlbcut_2b"; if(av) vs.push_back(mapname);
		if(plotD2){ mapname = "D2_2b"; if(av) vs.push_back(mapname); }
		if(plotD2){ mapname = "D2_withMlb_2b"; if(av) vs.push_back(mapname); }
		mapname = "MT2lb_withMlbcut_afterMT2llge85_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_afterMT2llge85_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlbcut_afterMT2llge85_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_afterMT2llge85_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlbcut_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_withMlbcut_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MlbAll_eq1b"; if(av) vs.push_back(mapname);
		mapname = "Mlb_eq1b"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2ll_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2ll_withMlbcut_eq1b"; if(av) vs.push_back(mapname);
		if(plotD2){ mapname = "D2_eq1b"; if(av) vs.push_back(mapname); }
		if(plotD2){ mapname = "D2_withMlb_eq1b"; if(av) vs.push_back(mapname); }
		mapname = "MT2lb_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlbcut_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_withMlbcut_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MlbAll_eq2b"; if(av) vs.push_back(mapname);
		mapname = "Mlb_eq2b"; if(av) vs.push_back(mapname);
		mapname = "Mll_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2ll_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2ll_withMlbcut_eq2b"; if(av) vs.push_back(mapname);
		if(plotD2){ mapname = "D2_eq2b"; if(av) vs.push_back(mapname); }
		if(plotD2){ mapname = "D2_withMlb_eq2b"; if(av) vs.push_back(mapname); }
		mapname = "MT2lb_withMlbcut_afterMT2llge85_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_afterMT2llge85_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlbcut_afterMT2llge85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_afterMT2llge85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85_eq2b"; if(av) vs.push_back(mapname);
		//>=1b
		mapname = "MT2lb_withMlb_bothl_true_trueb_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_onel_true_trueb_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_nol_true_trueb_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_bothb_true_truel_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_oneb_true_truel_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_nob_true_truel_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_bothl_true_trueb_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_onel_true_trueb_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_nol_true_trueb_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_bothb_true_truel_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_oneb_true_truel_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_nob_true_truel_1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_bothl_true_trueb_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_onel_true_trueb_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_nol_true_trueb_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_bothb_true_truel_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_oneb_true_truel_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_nob_true_truel_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_bothl_true_trueb_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_onel_true_trueb_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_nol_true_trueb_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_bothb_true_truel_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_oneb_true_truel_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_nob_true_truel_eq1b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_bothl_true_trueb_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_onel_true_trueb_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_nol_true_trueb_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_bothb_true_truel_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_oneb_true_truel_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_nob_true_truel_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_bothl_true_trueb_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_onel_true_trueb_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_nol_true_trueb_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_bothb_true_truel_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_oneb_true_truel_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_nob_true_truel_2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_bothl_true_trueb_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_onel_true_trueb_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_nol_true_trueb_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_bothb_true_truel_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_oneb_true_truel_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_nob_true_truel_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_bothl_true_trueb_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_onel_true_trueb_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_nol_true_trueb_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_bothb_true_truel_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_oneb_true_truel_eq2b"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlb_afterMT2llge85_nob_true_truel_eq2b"; if(av) vs.push_back(mapname);
		//>=1b
		mapname = "LepPt_1b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_1b"; if(av) vs.push_back(mapname);
		mapname = "BPt_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_1b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_1b"; if(av) vs.push_back(mapname);
		mapname = "BEta_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_1b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_gt200_1b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_gt200_1b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_gt200_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_gt200_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_gt200_1b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_gt200_1b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_gt200_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_gt200_1b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_lt200_1b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_lt200_1b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_lt200_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_lt200_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_lt200_1b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_lt200_1b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_lt200_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_lt200_1b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_gt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_gt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_gt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_gt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_gt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_gt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_gt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_gt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_lt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_lt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_lt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_lt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_lt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_lt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_lt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_lt180_MT2l_gt85_1b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_eq1b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BPt_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BEta_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_eq1b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_gt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_gt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_gt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_gt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_gt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_gt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_gt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_gt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_lt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_lt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_lt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_lt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_lt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_lt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_lt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_lt200_eq1b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_gt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_gt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_gt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_gt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_gt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_gt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_gt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_gt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_lt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_lt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_lt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_lt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_lt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_lt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_lt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_lt180_MT2l_gt85_eq1b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_2b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_2b"; if(av) vs.push_back(mapname);
		mapname = "BPt_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_2b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_2b"; if(av) vs.push_back(mapname);
		mapname = "BEta_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_2b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_gt200_2b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_gt200_2b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_gt200_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_gt200_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_gt200_2b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_gt200_2b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_gt200_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_gt200_2b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_lt200_2b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_lt200_2b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_lt200_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_lt200_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_lt200_2b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_lt200_2b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_lt200_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_lt200_2b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_gt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_gt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_gt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_gt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_gt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_gt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_gt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_gt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_lt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_lt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_lt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_lt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_lt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_lt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_lt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_lt180_MT2l_gt85_2b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_eq2b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BPt_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BEta_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_eq2b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_gt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_gt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_gt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_gt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_gt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_gt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_gt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_gt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_lt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_lt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_lt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_lt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_lt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_lt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_lt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_lt200_eq2b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_gt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_gt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_gt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_gt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_gt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_gt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_gt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_gt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "LepPt_MT2lbwMlb_lt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "LepEta_MT2lbwMlb_lt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lbwMlb_lt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelPt_MT2lbwMlb_lt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelTagDiscr_MT2lbwMlb_lt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr_MT2lbwMlb_lt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lbwMlb_lt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);
		mapname = "JSelEta_MT2lbwMlb_lt180_MT2l_gt85_eq2b"; if(av) vs.push_back(mapname);

*/


		//from StopFeasibilityNew.C
		//1d
/*		mapname = "LepPt"; if(av) vs.push_back(mapname);
		mapname = "LepMT"; if(av) vs.push_back(mapname);
		mapname = "BPt"; if(av) vs.push_back(mapname);
		mapname = "AllBPt"; if(av) vs.push_back(mapname);
		mapname = "NJetsID20"; if(av) vs.push_back(mapname);
		mapname = "NJetsID30"; if(av) vs.push_back(mapname);
		mapname = "NJetsID40"; if(av) vs.push_back(mapname);
		mapname = "NJetsID50"; if(av) vs.push_back(mapname);
		mapname = "NBJets20"; if(av) vs.push_back(mapname);
		mapname = "NBJets30"; if(av) vs.push_back(mapname);
		mapname = "NBJets40"; if(av) vs.push_back(mapname);
		mapname = "NBJets50"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr"; if(av) vs.push_back(mapname);
		mapname = "AllBTagDiscr"; if(av) vs.push_back(mapname);
		mapname = "BTagDiscr"; if(av) vs.push_back(mapname);
		mapname = "HT"; if(av) vs.push_back(mapname);
		mapname = "MET"; if(av) vs.push_back(mapname);
*/	//	mapname = "MT2lb"; if(av) vs.push_back(mapname);
	//	mapname = "MT2lb_massless"; if(av) vs.push_back(mapname);
	//	mapname = "MT2lb_withMlbcut";  if(av) vs.push_back(mapname);
//		mapname = "MT2lb_massless_withMlbcut"; if(av) vs.push_back(mapname);
//		mapname = "MlbAll"; if(av) vs.push_back(mapname);
	//	mapname = "Mlb"; if(av) vs.push_back(mapname);
	//	mapname = "Mll"; if(av) vs.push_back(mapname);
//		mapname = "Mbb"; if(av) vs.push_back(mapname);
//		mapname = "MbbAll"; if(av) vs.push_back(mapname);
	//	mapname = "MT2ll"; if(av) vs.push_back(mapname);
/*		mapname = "MT2ll_withMlbcut"; if(av) vs.push_back(mapname);
		mapname = "MT2bb"; if(av) vs.push_back(mapname);
		mapname = "MT2bb_lInMET"; if(av) vs.push_back(mapname);
		mapname = "MT2bb_lInMET_mW"; if(av) vs.push_back(mapname);
		mapname = "MT2bb_lInMET_mW_withMlbcut"; if(av) vs.push_back(mapname);
		mapname = "JetPt"; if(av) vs.push_back(mapname);
		mapname = "UTM"; if(av) vs.push_back(mapname);
		mapname = "DPhi_lbAll"; if(av) vs.push_back(mapname);
		mapname = "DPhi_lb"; if(av) vs.push_back(mapname);
		mapname = "DPhi_ll"; if(av) vs.push_back(mapname);
		mapname = "DPhi_bb"; if(av) vs.push_back(mapname);
		mapname = "DPhi_bbAll"; if(av) vs.push_back(mapname);
		mapname = "DR_lbAll"; if(av) vs.push_back(mapname);
		mapname = "DR_lb"; if(av) vs.push_back(mapname);
		mapname = "DR_ll"; if(av) vs.push_back(mapname);
		mapname = "DR_bb"; if(av) vs.push_back(mapname);
		mapname = "DR_bbAll"; if(av) vs.push_back(mapname);
		mapname = "Mlb_genmatchingTop"; if(av) vs.push_back(mapname);
		mapname = "MT2_genmatchingTops"; if(av) vs.push_back(mapname);
		mapname = "GenMlb"; if(av) vs.push_back(mapname);
		mapname = "GenMT2"; if(av) vs.push_back(mapname);
		mapname = "GenMT2_withMlbcut"; if(av) vs.push_back(mapname);
		mapname = "GenMlb_allAcceptance"; if(av) vs.push_back(mapname);
		mapname = "GenMT2_allAcceptance"; if(av) vs.push_back(mapname);
		mapname = "GenMT2_withMlbcut_allAcceptance"; if(av) vs.push_back(mapname);
		mapname = "GenMtop_allAcceptance"; if(av) vs.push_back(mapname);
		mapname = "GenMtop"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlbcut_PUle5"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlbcut_PUge9"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlbcut_woZveto"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_woZveto"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_woZveto"; if(av) vs.push_back(mapname);
		mapname = "MT2bb_woZveto"; if(av) vs.push_back(mapname);
		mapname = "MT2bb_linMET_woZveto"; if(av) vs.push_back(mapname);
		mapname = "MT2ll_woZveto"; if(av) vs.push_back(mapname);
		mapname = "No_PV"; if(av) vs.push_back(mapname);
*/		if(plotD2){ mapname = "D2"; if(av) vs.push_back(mapname); }
/*		if(plotD2){ mapname = "D2_withMlb"; if(av) vs.push_back(mapname); }

		mapname = "MT2lb_withMlbcut_afterMT2llge85"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_afterMT2llge85"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_afterMT2llge85"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_withMlbcut_afterMT2bb_linMETge150"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_afterMT2bb_linMETge150"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_massless_afterMT2bb_linMETge150"; if(av) vs.push_back(mapname);
		mapname = "MT2ll_afterMT2bb_linMETge150"; if(av) vs.push_back(mapname);
		mapname = "MT2ll_afterMT2lbge200"; if(av) vs.push_back(mapname);
		mapname = "MT2ll_afterMT2lb_withMlbcutge180"; if(av) vs.push_back(mapname);
		mapname = "MT2ll_afterMT2lb_masslessge140"; if(av) vs.push_back(mapname);
		mapname = "MT2bb_linMET_afterMT2llge85"; if(av) vs.push_back(mapname);
		mapname = "MT2bb_linMET_afterMT2lbge200"; if(av) vs.push_back(mapname);
		mapname = "MT2bb_linMET_afterMT2lb_withMlbcutge180"; if(av) vs.push_back(mapname);
		mapname = "MT2bb_linMET_afterMT2lb_masslessge140"; if(av) vs.push_back(mapname);
*/		if(plotD2){ mapname = "D2_afterMT2llge85"; if(av) vs.push_back(mapname); }
/*		if(plotD2){ mapname = "D2_afterMT2bb_linMETge150"; if(av) vs.push_back(mapname); }
		if(plotD2){ mapname = "D2_afterMT2lbge200"; if(av) vs.push_back(mapname); }
		if(plotD2){ mapname = "D2_afterMT2lb_withMlbcutge180"; if(av) vs.push_back(mapname); }
		if(plotD2){ mapname = "D2_afterMT2lb_masslessge140"; if(av) vs.push_back(mapname); }
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85"; if(av) vs.push_back(mapname);

		mapname = "MT2lb_bothbtags_truebjets"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_onebtag_truebjet"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_nonebtags_truebjets"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_bothbtags_matchedtobfromtop"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_onebtag_matchedtobfromtop"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_nonebtags_matchedtobfromtop"; if(av) vs.push_back(mapname);
		mapname = "RecoPt_over_GenPartonPt_truebfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "RecoPt_over_GenPartonPt_truebfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "Mlb_bothbtags_truebjets"; if(av) vs.push_back(mapname);
		mapname = "Mlb_onebtag_truebjet"; if(av) vs.push_back(mapname);
		mapname = "Mlb_nonebtags_truebjets"; if(av) vs.push_back(mapname);
		mapname = "Mlb_bothbtags_matchedtobfromtop"; if(av) vs.push_back(mapname);
		mapname = "Mlb_onebtag_matchedtobfromtop"; if(av) vs.push_back(mapname);
		mapname = "Mlb_nonebtags_matchedtobfromtop"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "BEta_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "BPt_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "BDiscr_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "BDiscr_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "TrueBfromtopEta_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "TrueBfromtopEta_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "TrueBfromtopPt_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "TrueBfromtopPt_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname  = "DeltaGenMET_PFMET_Phi_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname  = "DeltaGenMET_PFMET_Phi_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "AltMT2lb_MT2lb_gt_200_usingtruebpartons"; if(av) vs.push_back(mapname);
		mapname = "DeltaR_bjet_matchedtrueparton_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "DeltaR_bjet_matchedtrueparton_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "DeltaPt_bjet_matchedtrueparton_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "DeltaPt_bjet_matchedtrueparton_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_bothleptons_true_trueb"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_onelepton_true_trueb"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_nolepton_true_trueb"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_bothb_true_truel"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_oneb_true_truel"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_nob_true_truel"; if(av) vs.push_back(mapname);
		mapname = "DeltaR_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "DeltaR_TrueBparton_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "BDiscr_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "BEta_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "BPt_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "BDiscr_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "BEta_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "BPt_notmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "BDiscr_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "BEta_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "BPt_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "BDiscr_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "BEta_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "BPt_notselectedbutmatchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "MT2lb_bothbtags_matchedtosamebfromtop"; if(av) vs.push_back(mapname);
		mapname = "BDiscr_matchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "BPt_matchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "BEta_matchedjet_onebtag_matchedtobfromtop_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "BDiscr_matchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "BPt_matchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "BEta_matchedjet_onebtag_matchedtobfromtop_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "Eta_truebfromtop_matchedwithreco_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "Pt_truebfromtop_matchedwithreco_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "Eta_truebfromtop_notmatchedwithreco_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "Pt_truebfromtop_notmatchedwithreco_MT2lb_lt_200"; if(av) vs.push_back(mapname);
		mapname = "Eta_truebfromtop_matchedwithreco_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "Pt_truebfromtop_matchedwithreco_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "Eta_truebfromtop_notmatchedwithreco_MT2lb_gt_200"; if(av) vs.push_back(mapname);
		mapname = "Pt_truebfromtop_notmatchedwithreco_MT2lb_gt_200"; if(av) vs.push_back(mapname);

		mapname = "SingleTop_BPt_truebfromtop"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_BPt_otherb"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_LPt_truelfromtop"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_LPt_otherl"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_BEta_truebfromtop"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_BEta_otherb"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_LEta_truelfromtop"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_LEta_otherl"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_LID_otherl"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_LMID_otherl"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_LGMID_otherl"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_BID_otherb"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_BMID_otherb"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_BGMID_otherb"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_MT2lb_tchannel"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_MT2lb_tWchannel"; if(av) vs.push_back(mapname);
		mapname = "SingleTop_MT2lb_schannel"; if(av) vs.push_back(mapname);
*/
		vs1 = vs;
/*		//2d
		mapname = "DR_lb_configuration"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "trueb1b2_configuration"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "trueTop_configuration"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "Matchingconfiguration_MT2lb_lt_200"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "Matchingconfiguration_MT2lb_gt_200"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "Tail_DR_lb_configuration"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "Tail_trueb1b2_configuration"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "Tail_trueTop_configuration"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "LOriginID"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "LOriginMID"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "LOriginGMID"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "BOriginID"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "BOriginMID"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "BOriginGMID"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "LOrigin"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "BOrigin"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_vs_UTM"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_withMlb_vs_UTM"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_vs_BTagDiscr"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_vs_MT2ll"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_vs_MT2bb"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_vs_MT2bb_lInMET"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_vs_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_vs_MET"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_vs_HT"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_vs_PV"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		if(plotD2){ mapname = "MT2lb_massless_vs_D2"; if(av) vs.push_back(mapname);  vs2.push_back(mapname); }
		mapname = "MT2lb_massless_vs_Mlb"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_vs_BTagDiscr"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_vs_MT2ll"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_vs_MT2bb"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_vs_MT2bb_lInMET"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_vs_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_vs_MET"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_vs_HT"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_vs_PV"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		if(plotD2){ mapname = "MT2lb_vs_D2"; if(av) vs.push_back(mapname);  vs2.push_back(mapname); }
		mapname = "MT2lb_vs_Mlb"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_withMlbcut_vs_MT2ll"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_withMlbcut_vs_MT2bb"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_withMlbcut_vs_MT2bb_lInMET"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_withMlbcut_vs_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_withMlbcut_vs_MET"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_withMlbcut_vs_HT"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_withMlbcut_vs_PV"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		if(plotD2){ mapname = "MT2lb_withMlbcut_vs_D2"; if(av) vs.push_back(mapname);  vs2.push_back(mapname); }
		mapname = "MT2lb_withMlbcut_vs_MT2lb_massless"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_vs_MT2lb_massless"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2ll_vs_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2ll_vs_MT2bb_lInMET"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2ll_vs_MT2bb"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		if(plotD2){ mapname = "MT2ll_vs_D2"; if(av) vs.push_back(mapname);  vs2.push_back(mapname); }
		mapname = "MT2ll_vs_Mlb"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		if(plotD2){ mapname = "MT2bb_vs_D2"; if(av) vs.push_back(mapname);  vs2.push_back(mapname); }
		mapname = "MT2bb_vs_Mlb"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2bb_lInMET_vs_Mlb"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2bb_vs_BTagDiscr"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "Mlb_vs_BTagDiscr"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		if(plotD2){ mapname = "D2_vs_BTagDiscr"; if(av) vs.push_back(mapname);  vs2.push_back(mapname); }
		mapname = "PUWeight_vs_PV"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_withMlb_vs_MT2ll"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);

		mapname = "MT2lb_massless_vs_MT2ll_1b"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_withMlb_vs_MT2ll_1b"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_withMlb_vs_MT2ll_1b"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_vs_MT2ll_2b"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_withMlb_vs_MT2ll_2b"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_withMlb_vs_MT2ll_2b"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_vs_MT2ll_eq1b"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_withMlb_vs_MT2ll_eq1b"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_withMlb_vs_MT2ll_eq1b"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_vs_MT2ll_eq2b"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_massless_withMlb_vs_MT2ll_eq2b"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);
		mapname = "MT2lb_withMlb_vs_MT2ll_eq2b"; if(av) vs.push_back(mapname);  vs2.push_back(mapname);

		//3d
		mapname = "MT2lb_MT2ll_MT2bb"; if(av) vs.push_back(mapname);  vs3.push_back(mapname);
		mapname = "MT2lb_MT2ll_MT2bb_lInMET"; if(av) vs.push_back(mapname);  vs3.push_back(mapname);
		mapname = "MT2lb_MT2ll_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname);  vs3.push_back(mapname);
		mapname = "MT2lbMassless_MT2ll_MT2bb"; if(av) vs.push_back(mapname);  vs3.push_back(mapname);
		mapname = "MT2lbMassless_MT2ll_MT2bb_lInMET"; if(av) vs.push_back(mapname);  vs3.push_back(mapname);
		mapname = "MT2lbMassless_MT2ll_MT2bb_lInMET_mW"; if(av) vs.push_back(mapname);  vs3.push_back(mapname);
*/
		av = false;
		histonames = vs;

	for(int is = 0; is<sampletypesize; ++is){
           for(int il = 0; il<leptontypesize; ++il){
		string hs = string("_") + lepton_type[il] + string("_") + sample_type[is];
		for(unsigned int in = 0; in<vs1.size(); ++in){
			string name = vs1[in] + hs;
			if(histos.count(name)==0) histos[name ] = (TH1D*)oldfile->Get((name).c_str() );
		}
		for(unsigned int in = 0; in<vs2.size(); ++in){
			string name = vs2[in] + hs;
			if(histos2.count(name)==0){ histos2[name ] = (TH2D*)oldfile->Get((name).c_str() );
			}
		}
		for(unsigned int in = 0; in<vs3.size(); ++in){
			string name = vs3[in] + hs;
			if(histos3.count(name)==0){ histos3[name ] = (TH3D*)oldfile->Get((name).c_str() );
			}
		}
	   }
	}

	for(unsigned int n = 0; n<histonames.size(); ++n){
	   	for(int il = 0; il<leptontypesize; ++il){
			string h = string("_") + lepton_type[il];
			if(stacks.count(histonames[n]+h)==0) stacks[(histonames[n])+h ] = (THStack*)oldfile->Get((histonames[n]+h).c_str() );
		}
	}
	bool leggy = true;
	for(unsigned int n = 0; n<histonames.size(); ++n){
		if(leggy==false) break;
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

	cout << "plotting histograms ..." << endl;
	//MakePlots(tA, histos, histos2, sampletypesize, sample_type, leptontypesize, lepton_type, stacks, Legend1, histonames);
	MakePlots(histos, sampletypesize, sample_type, leptontypesize, lepton_type, stacks, Legend1, histonames);
	MakeCorrelationPlots(histos, histos2, sampletypesize, sample_type, leptontypesize, lepton_type, histonames, outputdir);
	Make3DPlots(histos3, leptontypesize, lepton_type, histonames, outputdir);

}

//void MakePlots(MassPlotter *tA, map<string, TH1D*> histos, const int sampletypesize, string *sample_type, const int leptontypesize, string *lepton_type, map<string, THStack*> stacks, TLegend *Legend1, vector<string> histonames){
void MakePlots(map<string, TH1D*> histos, const int sampletypesize, string *sample_type, const int leptontypesize, string *lepton_type, map<string, THStack*> stacks, TLegend *Legend1, vector<string> histonames){



	for(int n = 0; n<histonames.size(); ++n){//plot 1d
	   	for(int il = 0; il<leptontypesize; ++il){
			if(plotonlyLL && lepton_type[il]!="LL")  continue;
			if(plotonlyOF && lepton_type[il]!="EMu") continue;
			if(plotonlySF && lepton_type[il]!="SF")  continue;
			//if(lepton_type[il]=="EE" || lepton_type[il]=="MuMu") continue;
			string name = histonames[n];
			string h = string("_") + lepton_type[il];
			string hs1 = string("_") + lepton_type[il] + string("_") + string("mc");
			string hs2 = string("_") + lepton_type[il] + string("_") + string("data");
			string hs3 = string("_") + lepton_type[il] + string("_") + string("Stop");
			if(histos.count(name+hs1)==0) continue;
			TString ytitle = histos[name+hs1]->GetYaxis()->GetTitle();
			TString xtitle = histos[name+hs1]->GetXaxis()->GetTitle();

			TString outname = name + hs3 + (flip_order ? "_flipped" : "") + (logflag ? "_log" : "") + (composited ? "_comp" : "") + "_overlay";

//void MassPlotter::plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets, int nleps)
		if(histos[name+hs2]->Integral()>0 || histos[name+hs1]->Integral()>0){
		//	cout << name+hs3 << ": mc = " << histos[name+hs1]->Integral() << ", data = " << histos[name+hs2]->Integral() << ", stop = " << histos[name+hs3]->Integral() << endl;
			//if(!plotonlywithratio) tA->printHisto(stacks[name+h], histos[name+hs2],histos[name+hs1],histos[name+hs3], Legend1, outname, "hist", logflag, xtitle, ytitle, 2, 2, 1.);
			if(!plotonlywithratio) Make1DPlotsNoRatio(stacks[name+h], histos[name+hs1], histos[name+hs2], histos[name+hs3], logflag, false, outname, Legend1, xtitle, ytitle, 1.);

			//tA->plotRatioStack(stacks[name+h], histos[name+hs1], histos[name+hs2], histos[name+hs3], logflag, false, outname, Legend1, xtitle, ytitle, 2, 2, 1.);
			Make1DPlotsRatio(stacks[name+h], histos[name+hs1], histos[name+hs2], histos[name+hs3], logflag, false, outname, Legend1, xtitle, ytitle, 1.);

		}

	}}

}

void Make1DPlotsRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, TH1D* histsusy, bool logFlag, bool normalize, TString outname, TLegend *leg, TString xtitle, TString ytitle, float overlayScale){

	// define canvas and pads 
	TH1D *h1 = (TH1D*)histmc->Clone("h1_copy");
	TH1D *h2 = (TH1D*)histdata->Clone("h2_copy");
	TH1D *h3 = (TH1D*)histsusy->Clone("h3_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	h1->SetMarkerStyle(20);
	h2->SetMarkerStyle(20);
	
	//TCanvas* c1 = new TCanvas(name,"", 20,100,1000,700);
	TCanvas* c1 = new TCanvas(outname+"_ratio","",0,0,600,600 /*37, 60,636,670*/);
	c1->SetFrameLineWidth(1);
	c1 -> cd();
	
	float border = 0.2;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad(outname+"_plotpad",  "Pad containing the overlay plot", 0,0.211838,1,1 /*0.00, border, 1.00, 1.00, 0, 0*/);
 	//p_plot->SetBottomMargin(0.05);
	//p_plot->SetTopMargin(0.09);
	//p_plot->SetLeftMargin(0.1669107);
      	//p_plot->SetRightMargin(0.02);
	p_plot->SetLeftMargin(0.131579);
	p_plot->SetRightMargin(0.08);
	p_plot->SetTopMargin(0.06895515);
	p_plot->SetBottomMargin(0.07206074);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad(outname+"_ratiopad", "Pad containing the ratio",   0,0.01863354,0.9967105,0.2189441/*     0.00, 0.05, 1.00, border, 0, 0*/);
 	//p_ratio->SetTopMargin(0.03);
 	//p_ratio->SetBottomMargin(0.05/*5*/);
	//p_ratio->SetRightMargin(0.02);
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

//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);


	
	//MT2_bSel[0]->SetTitleSize(0.03);
	///MT2_bSel[0]->SetTitleOffset(1.);
	hstack->SetTitle("");
	h1    ->SetTitle("");
	h2    ->SetTitle("");
	h3    ->SetTitle("");
	hstack->SetMinimum(0.02);
	hstack->Draw("hist");
	hstack->GetYaxis()->SetTitle(ytitle.Data());
	hstack->GetYaxis()->SetLabelSize(0.05);
	hstack->GetYaxis()->SetTitleSize(0.05);
	hstack->GetYaxis()->SetTitleOffset(1.3);

/*	stringstream yTitle;
	if(fEventsPerGeV){
		if(fabs(h1_orig->GetBinWidth(1) -h1_orig->GetBinWidth(h1_orig->GetNbinsX()-1))<0.01){
		double binwidth = h1_orig->GetBinWidth(1);
		yTitle.precision(3);
		yTitle << ytitle.Data();
		yTitle << " / ";
		yTitle << binwidth;
		yTitle << " GeV";
		} else{
		cout << h1_orig->GetBinWidth(1) << " " << h1_orig->GetBinWidth(h1_orig->GetNbinsX()-1) << endl;
		}
	}else{
		yTitle << ytitle.Data();
	}

	hstack->GetYaxis()->SetTitle(yTitle.str().c_str());
*/

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
/*	text = fMT2Analysis?  "M_{T2} Analysis                                          ":"";
	text +=fMT2bAnalysis? "M_{T2}b Analysis                                         ":"";
	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
	text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";
*/
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
	h_ratio ->SetTitle("");	
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	//MC with errors
	TH1D*h_ratio_mc = (TH1D*)histmc->Clone("h1_copy_2");
	h_ratio_mc->Divide(h1);
	h_ratio_mc->SetTitle("");	
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

	//TLegend *leg = (TLegend*) legend->Clone("leg");

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
	//histmc -> Draw("same, E2");
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
/*
	TString text = fMT2Analysis? "M_{T2} Analysis        ":"";
	text +=fMT2bAnalysis? "M_{T2}b Analysis     ":"";
	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
	text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";
*/
	TitleBox.DrawLatex(0.18,0.943,text.Data());


	hstack->GetXaxis()->SetTitle(xtitle);
	hstack->GetXaxis()->SetLabelSize(0.05);
	hstack->GetXaxis()->SetTitleSize(0.05);
	hstack->GetXaxis()->SetTitleOffset(1.1);
/*	
	stringstream yTitle;
	if(fEventsPerGeV){
		if(fabs(histdata->GetBinWidth(1) -histdata->GetBinWidth(histdata->GetNbinsX()-1))<0.01){
		double binwidth = histdata->GetBinWidth(1);
		yTitle.precision(3);
		yTitle << ytitle.Data();
		yTitle << " / ";
		yTitle << binwidth;
		yTitle << " GeV";
		} else{
		cout << histdata->GetBinWidth(1) << " " << histdata->GetBinWidth(histdata->GetNbinsX()-1) << endl;
		}
	}else{
		yTitle << ytitle.Data();
	}
	hstack->GetYaxis()->SetTitle(yTitle.str().c_str());
*/
	hstack->GetYaxis()->SetTitle(ytitle);
	hstack->GetYaxis()->SetLabelSize(0.05);
	hstack->GetYaxis()->SetTitleSize(0.05);
	hstack->GetYaxis()->SetTitleOffset(1.47);

	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, outname, outputdir);
	if(fSave)Util::PrintEPS(col, outname, outputdir);
// 	delete col;

}


void MakeCorrelationPlots(map<string, TH1D*> histos, map<string, TH2D*> histos2, const int sampletypesize, string *sample_type, const int leptontypesize, string *lepton_type, vector<string> histonames, TString outputDir){

	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	cout << "Do Correlation plots" << endl;
	for(int n = 0; n<histonames.size(); ++n){//plot 1d
		//cout << histonames[n] << endl;
	   	for(int il = 0; il<leptontypesize; ++il){
			if(plotonlyLL && lepton_type[il]!="LL") continue;
			if(plotonlyOF && lepton_type[il]!="EMu") continue;
			if(plotonlySF && lepton_type[il]!="SF")  continue;
			string name = histonames[n];
			if(!(plotD2) && ( name=="MT2lb_massless_vs_D2" || name=="MT2lb_vs_D2" || name=="MT2lb_withMlbcut_vs_D2" || name=="MT2ll_vs_D2" || name=="MT2bb_vs_D2" || name=="D2_vs_BTagDiscr" ) ) continue;

			string hs0 = name + string("_") + lepton_type[il] + string("_") + string("TTbar");
			string hs1 = name + string("_") + lepton_type[il] + string("_") + string("mc");
			string hs2 = name + string("_") + lepton_type[il] + string("_") + string("data");
			string hs3 = name + string("_") + lepton_type[il] + string("_") + string("Stop");
			string hs4 = name + string("_") + lepton_type[il] + string("_") + string("SingleTop");
			col->cd();
			//cout << hs0 << " ";
			//cout << histos2.count(hs0) << " ";
			//cout << histos2[hs0]->GetEntries() << endl;
			if(logflag) gPad->SetLogz(1);
			if(histos2.count(hs0)>0) {
				col->SetName(hs0.c_str() );
				col->SetTitle(hs0.c_str() );
				histos2[hs0]->Draw("COLZ");
				col->Update();
				if(histos2[hs0]->GetEntries()>3) Util::PrintNoEPS(col, hs0, outputDir);
				if(histos2[hs0]->GetEntries()>3) Util::PrintEPS(col, hs0, outputDir);
				col->Clear();
			}
			if(histos2.count(hs1)>0) {
				col->SetName(hs1.c_str() );
				col->SetTitle(hs1.c_str() );
				histos2[hs1]->Draw("COLZ");
				col->Update();
				if(histos2[hs1]->GetEntries()>3) Util::PrintNoEPS(col, hs1, outputDir);
				if(histos2[hs1]->GetEntries()>3) Util::PrintEPS(col, hs1, outputDir);
				col->Clear();
			}
			if(histos2.count(hs2)>0) {
				col->SetName(hs2.c_str() );
				col->SetTitle(hs2.c_str() );
				histos2[hs2]->Draw("COLZ");
				col->Update();
				if(histos2[hs2]->GetEntries()>3) Util::PrintNoEPS(col, hs2, outputDir);
				if(histos2[hs2]->GetEntries()>3) Util::PrintEPS(col, hs2, outputDir);
				col->Clear();
			}
			if(histos2.count(hs3)>0) {
				col->SetName(hs3.c_str() );
				col->SetTitle(hs3.c_str() );
				histos2[hs3]->Draw("COLZ");
				col->Update();
				if(histos2[hs3]->GetEntries()>3) Util::PrintNoEPS(col, hs3, outputDir);
				if(histos2[hs3]->GetEntries()>3) Util::PrintEPS(col, hs3, outputDir);
				col->Clear();
			}
			if(plotsingletop && histos2.count(hs4)>0) {
				histos2[hs4]->SetMinimum(5.*10e-5);
				col->SetName(hs4.c_str() );
				col->SetTitle(hs4.c_str() );
				histos2[hs4]->Draw("COLZ");
				col->Update();
				if(histos2[hs4]->GetEntries()>0) Util::PrintNoEPS(col, hs4, outputDir);
				if(histos2[hs4]->GetEntries()>0) Util::PrintEPS(col, hs4, outputDir);
				col->Clear();
			}
		}
	}
	delete col;
}

void Make3DPlots(map<string, TH3D*> histos3, const int leptontypesize, string *lepton_type, vector<string> histonames, TString outputdirectory){

	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	cout << "Do Correlation plots" << endl;
	for(unsigned int n = 0; n<histonames.size(); ++n){//plot 1d
		//cout << histonames[n] << endl;
	   	for(int il = 0; il<leptontypesize; ++il){
			if(plotonlyLL && lepton_type[il]!="LL") continue;
			if(plotonlyOF && lepton_type[il]!="EMu") continue;
			if(plotonlySF && lepton_type[il]!="SF")  continue;
			string name = histonames[n];
			string hs0 = name + string("_") + lepton_type[il] + string("_") + string("TTbar");
			string hs1 = name + string("_") + lepton_type[il] + string("_") + string("mc");
			string hs2 = name + string("_") + lepton_type[il] + string("_") + string("data");
			string hs3 = name + string("_") + lepton_type[il] + string("_") + string("Stop");
			string hs4 = name + string("_") + lepton_type[il] + string("_") + string("SingleTop");
			col->cd();
			//cout << hs0 << " ";
			//cout << histos2.count(hs0) << " ";
			//cout << histos2[hs0]->GetEntries() << endl;
			//if(logflag) gPad->SetLogz(1);
			if(histos3.count(hs0)>0) {
			//	histos3[hs0]->SetMinimum(0.005);
				col->SetName(hs0.c_str() );
				col->SetTitle(hs0.c_str() );
				histos3[hs0]->Draw("BOX");
				col->Update();
				if(histos3[hs0]->GetEntries()>0) Util::PrintNoEPS(col, hs0, outputdirectory);
				if(histos3[hs0]->GetEntries()>0) Util::PrintEPS(col, hs0, outputdirectory);
				col->Clear();
			}
			if(histos3.count(hs1)>0) {
			//	histos3[hs1]->SetMinimum(0.005);
				col->SetName(hs1.c_str() );
				col->SetTitle(hs1.c_str() );
				histos3[hs1]->Draw("BOX");
				col->Update();
				if(histos3[hs1]->GetEntries()>0) Util::PrintNoEPS(col, hs1, outputdirectory);
				if(histos3[hs1]->GetEntries()>0) Util::PrintEPS(col, hs1, outputdirectory);
				col->Clear();
			}
			if(histos3.count(hs2)>0) {
			//	histos3[hs2]->SetMinimum(0.5);
				col->SetName(hs2.c_str() );
				col->SetTitle(hs2.c_str() );
				histos3[hs2]->Draw("BOX");
				col->Update();
				if(histos3[hs2]->GetEntries()>0) Util::PrintNoEPS(col, hs2, outputdirectory);
				if(histos3[hs2]->GetEntries()>0) Util::PrintEPS(col, hs2, outputdirectory);
				col->Clear();
			}
			if(histos3.count(hs3)>0) {
			//	histos3[hs3]->SetMinimum(5.*10e-5);
				col->SetName(hs3.c_str() );
				col->SetTitle(hs3.c_str() );
				histos3[hs3]->Draw("BOX");
				col->Update();
				if(histos3[hs3]->GetEntries()>0) Util::PrintNoEPS(col, hs3, outputdirectory);
				if(histos3[hs3]->GetEntries()>0) Util::PrintEPS(col, hs3, outputdirectory);
				col->Clear();
			}
			if(plotsingletop && histos3.count(hs4)>0) {
				histos3[hs4]->SetMinimum(5.*10e-5);
				col->SetName(hs4.c_str() );
				col->SetTitle(hs4.c_str() );
				histos3[hs4]->Draw("BOX");
				col->Update();
				if(histos3[hs4]->GetEntries()>0) Util::PrintNoEPS(col, hs4, outputdirectory);
				if(histos3[hs4]->GetEntries()>0) Util::PrintEPS(col, hs4, outputdirectory);
				col->Clear();
			}
		}
	}
	delete col;
}
