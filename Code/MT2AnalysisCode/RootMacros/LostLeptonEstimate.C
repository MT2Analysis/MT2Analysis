#include "TEfficiency.h"
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
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TMap.h"
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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/Utilities.hh"


using namespace std;

//use via root -l -b -q LostLeptonEstimate.C++

void load(const char* filename = "samples/datasamples/samples_2141_dataonly.dat");
void LostLeptonEstimate();
TH1D* RebinThisHistogram(TH1D* histogram, int nrebinnedbins, double* rebinnedbins);
TEfficiency* RebinThisEfficiency(TEfficiency* tefficiency, int nrebinnedbins, double* rebinnedbins);
void GetLostLeptonEstimate(map<string, TH1D*> histograms, map<string, TEfficiency*> tefficiencies, double rel_sys_uncert,  double rel_sys_uncert_bg, int version, Bool_t PrintSummaryTable, Bool_t PrintYieldTable, Bool_t PrintPredictionCard, Bool_t makeFullPrintout=false, Bool_t SavePrediction=true, Bool_t PlotPrediction=false);
void PredictionCard(vector<int> sr, vector<int> htr, vector<int> lepr, vector<int> MT2bin, vector<double> mctruth, vector<double> mctrutherr, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, vector<double> LLeff, vector<double> LLefferr, vector<double> MTeff, vector<double> MTefferr);
void PredictionCardSplitted(vector<int> sr, vector<int> htr, vector<int> lepr, vector<int> MT2bin, vector<double> mctruth, vector<double> mctrutherr, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, vector<double> LLeff, vector<double> LLefferr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruthfractionW, vector<double> mctruthfractionTT, vector<double> numW, vector<double> numTT, vector<double> numMC);
void YieldTable(vector<int> sr, vector<int> htr, vector<int> lepr, vector<double> MT2low, vector<double> MT2up, vector<double> numQCD, vector<double> numZ, vector<double> numW, vector<double> numT, vector<double> numOther, vector<double> numMC, vector<double> BGerr, vector<double> numData);
void SummaryTable(int version, vector<int> sr, vector<int> htr, vector<int> lepr, vector<double> MT2low, vector<double> MT2up, vector<double> numtrueW, vector<double> numtrueW_bg, vector<double> numData, vector<double> numBG, vector<double> LLeff, vector<double> LLefferr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruth, vector<double> mctrutherr, vector<double> MCpred, vector<double> MCpred_stat_err, vector<double> MCpred_syst_err, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err);
void PredictionFile(int version, vector<int> sr, vector<int> htr, vector<int> lepr, vector<double> MT2low, vector<double> MT2up, vector<double> numtrueW, vector<double> numtrueW_bg, vector<double> numData, vector<double> numBG, vector<double> LLeff, vector<double> LLefferr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruth, vector<double> mctrutherr, vector<double> MCpred, vector<double> MCpred_stat_err, vector<double> MCpred_syst_err, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, Bool_t PlotPrediction);

//standard sample struct that combines the MT2trees with necessary information like sample type, cross section, etc.
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

vector < sample >  fSamples;//note that in this code the samples.dat is hard-coded
const int fVerbose = 3;
TString fPath;

Bool_t fMET        = false;//use MET triggers (i.e. used for low HT region)
Bool_t fHT         = true; //use HT triggers (i.e. used for medium+high HT region), will be set false automatically if fMET==true
Bool_t fFast       = true;//includes a lower MT2 cut (MT2>100 GeV), default = true, as low MT2 region not part of any signal region
Bool_t fRebin      = true;//do not bin-by-bin estimation but one estimate per HT/topological region (this means not along MT2)

Bool_t fISRreweight=true; //reweight MC according to SUSY's 'ISR recipe' - influence on estimate is minimal, but MC truth changes
                          //keep that flag true - then you will save both the ISR reweighted and non-reweighted style, later use LostLeptonEstimateFromRootfile to get the number for your chosen case

Bool_t fData       = true;//set false if don't want to run over data - only if doing checks
Bool_t fQCD        = false;//set false if don't want to run over QCD samples - only if doing checks and want to save time
Bool_t fSusy       = false;//set false if don't want to run over Susy samples - default is false

Bool_t fbTagReweight= true;//reweight MC according to BTV SF weights - default is true
Bool_t fbTagError   = true;//compute additionally error due to BTV SF uncertainty - default is true

Double_t frel_sys_uncert    = 0.05;//uncertainty on lepton efficiency (reconstruction and acceptance), default is 5 percent
Double_t frel_sys_uncert_bg =  0.5;//uncertainty due to background subtraction, default is 50 percent
Double_t fdoubeLLerr        =  1.0;//uncertainty due to double lost leptons, default is 100 percent - use at own risk (there is no good recipe for that)
Double_t frelMTerr          = 0.05;//uncertainty due to MT cut on lepton-MET system, default is 5 percent

std::ostringstream* fLogStream     = 0;
Bool_t  fWriteToFile               = false; // writes couts to a file
Bool_t  fAppend                    = true;  // append at end of file (if existing), otherwise delete previous content - needs fWriteToFile = true
TString  outputdir = "Filtered/LostLepton/tryout/";//create output directory

Bool_t fUsedoubleLL               = false;	// event called good if there are two generated leptons (i.e. not only single lepton events) - this is optional and has no big influence
Bool_t fDoubleLL_BG               = false;	// events that are double lost leptons (two lost leptons) are background, if false (default) they are part of yield you want to predict
Bool_t fIncludeTaus               = true;	// include leptonic tau decays to LostLeptonEstimate // default = true;
Bool_t fTopEfficencies            = true;	//lepton efficiency from top sample only, if false only from W sample; works only if fWeightedProb=false; used for debugging default = true;
Bool_t  fWeightedProb             = true;	// this should be true at all times. Take LostLepton efficiency from both W+Top sample
Bool_t  fIncludeTop               = true;	// this should be true at all times. LostLeption yield estimated for Top+W sample; if false take only W sample
Bool_t  fTopOnly                  = false;	// this should be false at all times. LostLeption yield estimated for Top sample only
Bool_t  fIncludeSingleTop         = true;	// this should be true at all times. Top sample contains also single top, not only ttbar.

    //definition of all signal bins
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

    // HT > 450 && HT < 750
    const int gNMT2bins_2j0b_lHT                      = 8;
    double  gMT2bins_2j0b_lHT[gNMT2bins_2j0b_lHT+1]   = {200, 240, 290, 350, 420, 490, 570, 650, 750};
    const int gNMT2bins_2j1b_lHT                      = 6;
    double  gMT2bins_2j1b_lHT[gNMT2bins_2j1b_lHT+1]   = {200, 250, 310, 380, 450, 550, 700};
    const int gNMT2bins_3j0b_lHT                      = 8;
    double  gMT2bins_3j0b_lHT[gNMT2bins_3j0b_lHT+1]   = {200, 240, 290, 350, 420, 490, 570, 650, 750};
    const int gNMT2bins_3j1b_lHT                      = 6;
    double  gMT2bins_3j1b_lHT[gNMT2bins_3j1b_lHT+1]   = {200, 250, 310, 380, 460, 550, 700};
    const int gNMT2bins_3j2b_lHT                      = 4;
    double  gMT2bins_3j2b_lHT[gNMT2bins_3j2b_lHT+1]   = {200, 250, 325, 425, 550};
    const int gNMT2bins_6j0b_lHT                      = 3;
    double  gMT2bins_6j0b_lHT[gNMT2bins_6j0b_lHT+1]   = {200, 280, 380, 520};
    const int gNMT2bins_6j1b_lHT                      = 3;
    double  gMT2bins_6j1b_lHT[gNMT2bins_6j1b_lHT+1]   = {200, 250, 325, 450};
    const int gNMT2bins_6j2b_lHT                      = 3;
    double  gMT2bins_6j2b_lHT[gNMT2bins_6j2b_lHT+1]   = {200, 250, 300, 400};
    const int gNMT2bins_3b_lHT                        = 2;
    double  gMT2bins_3b_lHT  [gNMT2bins_3b_lHT+1]     = {200, 280, 400};

//definition of all samples (and helpsamples like WandTop)
const int sampletypesize = 12;
string sample_type[sampletypesize] = {"QCD", "WJets", "ZJets", "TTbar", "SingleTop", "Top", "WandTop", "noWandTop", "Other", "mc", "susy", "data"};
//definition of topological regions
const int signalregionsize = 9;
string signal_region[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b"};
//definition of HT regions, the names are set by the fMET/fHT flags above
const int HTbinsize = 2;
string HT_bin[HTbinsize] = {"HTge450", "HTge750"};//dummy
//definition of lepton types, taus are separate, as they have other uncertainties since their tagging nature (i.e. fake tau rate) is much higher
const int leptontypesize = 2;
string lepton_type[leptontypesize] = {"Muo", "Ele"};


//This is the main code!!!
//This code does the LostLeptonEstimate for electrons and muons, taus are done with TauEstimation.C (due to different type of calling the gen-info and different uncertainties --> an integration here would for sure be possible if wanted)
//Basically only some histograms are filled and later (in another function of this macro) all relevant quantities are extracted out of all these histograms, depending on the flags you set before.
void LostLeptonEstimate(){

	// logStream - store it to file or write it out to terminal - depends on fWriteToFile
	fLogStream = new std::ostringstream();

	if(fMET){ HT_bin[0] = "HTge450"; }
	if(fHT) { HT_bin[0] = "HTge750"; HT_bin[1] = "HTge1200"; }
	if(fMET==true) fHT = false;				//either fHT or fMET has to be false
	if(fMET==false && fHT==false) fHT = true;		//either fHT or fMET has to be true
	if(fIncludeTop==false) fIncludeSingleTop = false;	//cannot have single top if rejecting top
	if(fIncludeTop==false) fTopOnly = false;		//cannot compute purely top if rejecting top
	if(fbTagReweight==false) fbTagError = false;		//there are no btagging errors without btagging weights
  	gROOT->ProcessLine(".x SetStyle_PRD.C");

	//create output directory
	if(fMET) outputdir = outputdir + "MET/";
	if(fHT)  outputdir = outputdir + "HT/";
	if(fUsedoubleLL) outputdir += "used_dLL_forEff/";
	if(fDoubleLL_BG) outputdir += "veto_dLL/";
    	Util::MakeOutputDir(outputdir);

	//file where all histograms are stored - you can use this file later to do all calculations - i.e. you do not need to rerun every time you reset a flag
	TString outputname = "LostLeptonHistograms.root";

	//the samples.dat files containing all samples you want to run over
	TString  samples = "dummy.dat";//only dummy
	if(fMET) samples = "samples/samples_MET_filter.dat";
	if(fHT)  samples = "samples/samples_HT_filter.dat";

	map<string, TH1D*>    histos;
	map<string, TEfficiency*>    teff;
	vector<string> histonames; histonames.clear();
	//get all histogram names you need
	histonames.push_back("RecoLepEvents");			//all reco leptons (i.e. 1muo, 0ele or 0muo, 1 ele, incl. MT cut)
	histonames.push_back("LeptonEvents");			//all 1 lepton events (i.e. == RecoLepEvents) - histograms used for different purposes, therefore kept both
	histonames.push_back("LeptonEventsNoMT");		//all 1 lepton events (i.e. == RecoLepEvents) without cutting on MT
	histonames.push_back("LeptonEventsBG");			//all 1 lepton event, but where the lepton is fake (not from W), like non-prompt leptons from heavy flavour decays
	histonames.push_back("LeptonEvents_doubleLL");		//all 1 lepton events from W/Top events with one lost lepton (double leptonic)
	histonames.push_back("NoLeptonEvents");			//all 0 lepton events but where there is one genlepton - i.e. these are the lost leptons
	histonames.push_back("NoLeptonEvents_doubleLL");	//all 0 lepton events but where there is two genlepton from W/Top decay (double[two] lost leptons)
	histonames.push_back("WEvents");			//all events with one W(lnu) (also from  top decays)
	histonames.push_back("WEventsAcc");			//all events with one W(lnu), where the charged lepton is within detector acceptance
	histonames.push_back("WEventsReco");			//all events with one W(lnu), where the charge lepton is reconstructed (within detector acceptance), and MT<100 GeV
	histonames.push_back("WEventsRecoNoMT");		//as above, but use all leptons, i.e. also MT>100 GeV are accepted
	histonames.push_back("WEvents_dLL");			//as above, but select events with more than one W(lnu)
	histonames.push_back("WEventsReco_dLL");		//as above, but select events with more than one W(lnu)
	//the following histograms are needed to assign the uncertainty due to BTV SF
	histonames.push_back("LeptonEventsBGUp");		//B-Tag SF + SFerrup   <-- Getting SF error
	histonames.push_back("RecoLepEventsUp");		//B-Tag SF + SFerrup   <-- Getting SF error
	histonames.push_back("LeptonEventsBGDown");		//B-Tag SF - SFerrdown <-- Getting SF error
	histonames.push_back("RecoLepEventsDown");		//B-Tag SF - SFerrdown <-- Getting SF error
	histonames.push_back("NoLeptonEvents_SFup");		//B-Tag SF + SFerrup   <-- Getting SF error
	histonames.push_back("NoLeptonEvents_SFdown");		//B-Tag SF - SFerrdown <-- Getting SF error
	histonames.push_back("LeptonEvents_doubleLL_SFup");	//B-Tag SF + SFerrup   <-- Getting SF error
	histonames.push_back("LeptonEvents_doubleLL_SFdown");	//B-Tag SF - SFerrdown <-- Getting SF error
	histonames.push_back("NoLeptonEvents_doubleLL_SFup");	//B-Tag SF + SFerrup   <-- Getting SF error
	histonames.push_back("NoLeptonEvents_doubleLL_SFdown");	//B-Tag SF - SFerrdown <-- Getting SF error
	vector<string> teffnames; teffnames.clear();//define them later

	//now all histograms are created systematically
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTbinsize;        ++i3){
	for(int i4 = 0; i4<leptontypesize;   ++i4){
		if(fMET && i3==1) continue;
		int NMT2bins;
		if(fMET){
			if(signal_region[i2]=="2j0b")    NMT2bins = gNMT2bins_2j0b_lHT;
			if(signal_region[i2]=="2j1to2b") NMT2bins = gNMT2bins_2j1b_lHT;
			if(signal_region[i2]=="3to5j0b") NMT2bins = gNMT2bins_3j0b_lHT;
			if(signal_region[i2]=="3to5j1b") NMT2bins = gNMT2bins_3j1b_lHT;
			if(signal_region[i2]=="3to5j2b") NMT2bins = gNMT2bins_3j2b_lHT;
			if(signal_region[i2]=="6j0b")    NMT2bins = gNMT2bins_6j0b_lHT;
			if(signal_region[i2]=="6j1b")    NMT2bins = gNMT2bins_6j1b_lHT;
			if(signal_region[i2]=="6j2b")    NMT2bins = gNMT2bins_6j2b_lHT;
			if(signal_region[i2]=="3b")      NMT2bins = gNMT2bins_3b_lHT;
		} if(fHT){
		   if(i3==0){
			if(signal_region[i2]=="2j0b")    NMT2bins = gNMT2bins_2j0b_mHT;
			if(signal_region[i2]=="2j1to2b") NMT2bins = gNMT2bins_2j1b_mHT;
			if(signal_region[i2]=="3to5j0b") NMT2bins = gNMT2bins_3j0b_mHT;
			if(signal_region[i2]=="3to5j1b") NMT2bins = gNMT2bins_3j1b_mHT;
			if(signal_region[i2]=="3to5j2b") NMT2bins = gNMT2bins_3j2b_mHT;
			if(signal_region[i2]=="6j0b")    NMT2bins = gNMT2bins_6j0b_mHT;
			if(signal_region[i2]=="6j1b")    NMT2bins = gNMT2bins_6j1b_mHT;
			if(signal_region[i2]=="6j2b")    NMT2bins = gNMT2bins_6j2b_mHT;
			if(signal_region[i2]=="3b")      NMT2bins = gNMT2bins_3b_mHT;
		   } if(i3==1){
			if(signal_region[i2]=="2j0b")    NMT2bins = gNMT2bins_2j0b_hHT;
			if(signal_region[i2]=="2j1to2b") NMT2bins = gNMT2bins_2j1b_hHT;
			if(signal_region[i2]=="3to5j0b") NMT2bins = gNMT2bins_3j0b_hHT;
			if(signal_region[i2]=="3to5j1b") NMT2bins = gNMT2bins_3j1b_hHT;
			if(signal_region[i2]=="3to5j2b") NMT2bins = gNMT2bins_3j2b_hHT;
			if(signal_region[i2]=="6j0b")    NMT2bins = gNMT2bins_6j0b_hHT;
			if(signal_region[i2]=="6j1b")    NMT2bins = gNMT2bins_6j1b_hHT;
			if(signal_region[i2]=="6j2b")    NMT2bins = gNMT2bins_6j2b_hHT;
			if(signal_region[i2]=="3b")      NMT2bins = gNMT2bins_3b_hHT;
		   }
		}
  		double MT2bins[NMT2bins+1];
		if(fMET){
			if(signal_region[i2]=="2j0b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_2j0b_lHT[i0]; }
			if(signal_region[i2]=="2j1to2b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_2j1b_lHT[i0]; }
			if(signal_region[i2]=="3to5j0b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j0b_lHT[i0]; }
			if(signal_region[i2]=="3to5j1b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j1b_lHT[i0]; }
			if(signal_region[i2]=="3to5j2b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j2b_lHT[i0]; }
			if(signal_region[i2]=="6j0b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j0b_lHT[i0]; }
			if(signal_region[i2]=="6j1b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j1b_lHT[i0]; }
			if(signal_region[i2]=="6j2b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j2b_lHT[i0]; }
			if(signal_region[i2]=="3b")      { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3b_lHT[i0];   }
		} if(fHT){
		   if(i3==0){
			if(signal_region[i2]=="2j0b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_2j0b_mHT[i0];  }
			if(signal_region[i2]=="2j1to2b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_2j1b_mHT[i0];  }
			if(signal_region[i2]=="3to5j0b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j0b_mHT[i0];  }
			if(signal_region[i2]=="3to5j1b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j1b_mHT[i0];  }
			if(signal_region[i2]=="3to5j2b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j2b_mHT[i0];  }
			if(signal_region[i2]=="6j0b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j0b_mHT[i0];  }
			if(signal_region[i2]=="6j1b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j1b_mHT[i0];  }
			if(signal_region[i2]=="6j2b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j2b_mHT[i0];  }
			if(signal_region[i2]=="3b")      { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3b_mHT[i0];    }
		   } if(i3==1){
			if(signal_region[i2]=="2j0b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_2j0b_hHT[i0];  }
			if(signal_region[i2]=="2j1to2b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_2j1b_hHT[i0];  }
			if(signal_region[i2]=="3to5j0b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j0b_hHT[i0];  }
			if(signal_region[i2]=="3to5j1b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j1b_hHT[i0];  }
			if(signal_region[i2]=="3to5j2b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j2b_hHT[i0];  }
			if(signal_region[i2]=="6j0b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j0b_hHT[i0];  }
			if(signal_region[i2]=="6j1b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j1b_hHT[i0];  }
			if(signal_region[i2]=="6j2b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j2b_hHT[i0];  }
			if(signal_region[i2]=="3b")      { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3b_hHT[i0];    }
		   }
		}
		string hs = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_") + sample_type[i1];
		string mapname;
		for(unsigned int i0 = 0; i0<histonames.size(); ++i0){
			mapname = histonames[i0] + hs;
			if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", NMT2bins, MT2bins);
			mapname = "noISR_"+histonames[i0] + hs;
			if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", NMT2bins, MT2bins);
		}
	}}}}

	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Sumw2();}

	//selection cuts, you only run over the events passing this selection
	//this selection contains both the "lost leptons" and the "one lepton control sample"
	std::ostringstream cutStream;
	std::ostringstream cutStreamBase;
	cutStream << " " 
		<< "NTausIDLoose3Hits==0"                   << "&&"
		<< "misc.Jet0Pass ==1"                      << "&&"
		<< "misc.Jet1Pass ==1"                      << "&&"
		<< "misc.Vectorsumpt < 70";
		cutStream  << "&& misc.MinMetJetDPhi4Pt40 >0.3";
//	if(fMET) cutStream << "&&misc.MET>200&&misc.HT<=750";
	if(fMET) cutStream << "&&misc.MET>200&&misc.HT<=750&&misc.HT>=450&&misc.MT2>200";
	if(fHT ) cutStream << "&&misc.HT>750&&misc.MET>30";
	if(fFast)cutStream << "&&misc.MT2>=100";//lowest border in MT2
	
	cutStreamBase << " " 
      << "misc.PassJet40ID ==1"                      << "&&"
      << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"  << "&&" // signal samples (fastsim) do not have this filter
      << "misc.CSCTightHaloIDFlag == 0"             << "&&"
      << "misc.trackingFailureFlag==0"              << "&&"
      << "misc.eeBadScFlag==0"                      << "&&"
      << "misc.EcalDeadCellTriggerPrimitiveFlag==0" << "&&"
      << "misc.TrackingManyStripClusFlag==0"             << "&&"
      << "misc.TrackingTooManyStripClusFlag==0"          << "&&"
      << "misc.TrackingLogErrorTooManyClustersFlag==0"   << "&&"
      << "misc.CrazyHCAL==0";
      cutStreamBase << "&&misc.MET/misc.CaloMETRaw<=2.";//HO cut
	
	std::ostringstream triggerStream;
	if(fMET){
	triggerStream << "( ( ("
			<< "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
			<< "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )"
			<< "||("
			<< "trigger.HLT_PFHT350_PFMET100_v3==1 || trigger.HLT_PFHT350_PFMET100_v4==1 || trigger.HLT_PFHT350_PFMET100_v5==1 || "
			<< "trigger.HLT_PFHT350_PFMET100_v6==1 || trigger.HLT_PFHT350_PFMET100_v7==1 || trigger.HLT_PFNoPUHT350_PFMET100_v1==1 || "
			<< "trigger.HLT_PFNoPUHT350_PFMET100_v3==1 || trigger.HLT_PFNoPUHT350_PFMET100_v4==1 ) ) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	}
	if(fHT){
	triggerStream << "( ("
			<< "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || "
			<< "trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
			<< "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	}
	TString trigger = triggerStream.str().c_str();
	
	TString cuts = cutStream.str().c_str();
	TString basecuts = cutStreamBase.str().c_str();

	//load the samples you want to run over
	load(samples.Data());

   	for(size_t i = 0; i < fSamples.size(); ++i){
        
   	    if(fData==false && fSamples[i].type=="data") continue;
	    if(fSusy==false && fSamples[i].type=="susy") continue;
	    if(fQCD==false  && fSamples[i].sname=="QCD") continue;
        
		//assign all samples to its type
	    string sampletype = (string)fSamples[i].type;
	    if(sampletype==(string)"mc"){
		if(fSamples[i].sname=="QCD")         sampletype = (string)"QCD";
		else if(fSamples[i].sname=="Wtolnu") sampletype = (string)"WJets";
		else if(fSamples[i].sname=="DY")     sampletype = (string)"ZJets";
		else if(fSamples[i].name=="TTbar")   sampletype = (string)"TTbar";
		else if(fSamples[i].name=="TTbar_Madgraph0l")   sampletype = (string)"TTbar";
		else if(fSamples[i].name=="TTbar_Madgraph1l")   sampletype = (string)"TTbar";
		else if(fSamples[i].name=="TTbar_Madgraph2l")   sampletype = (string)"TTbar";
		else if(fSamples[i].sname=="Top")    sampletype = (string)"SingleTop";//no ttbar, includes TTZ, TTW
		else sampletype = (string)"Other";
	    }

	    //global sample weight
	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "LostLepton: looping over " << fSamples[i].name << " added in " << sampletype << endl;
	    if(fVerbose>2) cout << "            sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

	    //definition of MT2tree
	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts + "&&" + basecuts;
        
	    if( fSamples[i].type=="data") myCuts += " && " + trigger; //cuts to be aplied only on data
        
   	    cout << "Cuts for Flow: " << myCuts << endl;
   	    fSamples[i].tree->Draw(">>selList", myCuts);
        
	    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
	    fSamples[i].tree->SetEventList(myEvtList);
	    int counter=0;
	    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
	    if(myEvtList->GetSize()==0) continue;
	//run only over events passing event selection
        while(myEvtList->GetEntry(counter++) !=-1){	
		int jentry = myEvtList->GetEntry(counter-1);
            
		nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
		fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
		
		if ( fVerbose>2 && !fFast && counter % 5000 == 0  )  cout << "+++ Proccessing event " << counter << endl;
		
		Double_t weight = sample_weight;
		if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 

		//the if clause below computes the weight due to the 'ISR recipe'
		Double_t ISRweight(1.); Double_t weightnoISR = weight;
		if(fISRreweight && !fMT2tree->misc.isData){
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
			} if(sampletype=="ZJets"){
				hardgenlv = fMT2tree->GenZ[0];
			} if(sampletype=="TTbar"){
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
			} if(sampletype=="SingleTop"){
				if((fSamples[i].name).Contains("tW")){//t + W
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
			if(hardgenlv.Pt()>250.) ISRweight = 0.8;
			else if(hardgenlv.Pt()>150.) ISRweight = 0.9;
			else if(hardgenlv.Pt()>120.) ISRweight = 0.95;
			else                         ISRweight = 1.;
			weight = weight * ISRweight;
		}

		if(fMT2tree->misc.HT>=750.&&fMT2tree->misc.HT<1200.&&fMT2tree->NJetsIDLoose40==2&&fMT2tree->NBJets40CSVM>=1&&fMT2tree->misc.MT2>=135.&&fMT2tree->misc.MT2<=170.&&sampletype=="SingleTop"&&fMT2tree->misc.LumiSection==149&&fMT2tree->misc.Event==44699) continue;//sample T_t_highHT - fix for one 8 TeV MC sample with crazy weight

		//define lepton case (reconstruction)
		Bool_t recoedele       = false;// exact 1 ele, 0 muo including MT cut on selected ele (MT<100)
            	Bool_t recoedmuo       = false;// exact 1 muo, 0 ele including MT cut on selected muo (MT<100)
		Bool_t recoedelenomt   = false;// exact 1 ele, 0 muo without MT cut on selected ele
            	Bool_t recoedmuonomt   = false;// exact 1 muo, 0 ele without MT cut on selected muo
		Bool_t norecolep       = false;// 0 ele and 0 muo
            	if(fMT2tree->NEles ==1 && fMT2tree->NMuons==0 && fMT2tree->ele[0].MT<100. ) recoedele     = true;
            	if(fMT2tree->NMuons==1 && fMT2tree->NEles ==0 && fMT2tree->muo[0].MT<100. ) recoedmuo     = true;
            	if(fMT2tree->NEles ==1 && fMT2tree->NMuons==0                             ) recoedelenomt = true;
            	if(fMT2tree->NMuons==1 && fMT2tree->NEles ==0                             ) recoedmuonomt = true;
		if(fMT2tree->NMuons==0 && fMT2tree->NEles ==0)                              norecolep     = true;

		//define region for histogram naming
		string slep;
		if(recoedele) slep = "_Ele";
		if(recoedmuo) slep = "_Muo";
		string sHT;
		if(fMET){
			if(fMT2tree->misc.HT<=450.)      sHT = "_HTge0";
			else if(fMT2tree->misc.HT<=750.) sHT = "_HTge450";
		} if(fHT){
			if(fMT2tree->misc.HT>1200.)      sHT = "_HTge1200";
			else if(fMT2tree->misc.HT>750.)  sHT = "_HTge750";
		}
		string ssignal;
		if(fMT2tree->NJetsIDLoose40 == 2 &&                                  fMT2tree->NBJets40CSVM == 0) ssignal = "_2j0b";
		if(fMT2tree->NJetsIDLoose40 == 2 &&                                  fMT2tree->NBJets40CSVM >= 1) ssignal = "_2j1to2b";
		if(fMT2tree->NJetsIDLoose40 >= 3 && fMT2tree->NJetsIDLoose40 <= 5 && fMT2tree->NBJets40CSVM == 0) ssignal = "_3to5j0b";
		if(fMT2tree->NJetsIDLoose40 >= 3 && fMT2tree->NJetsIDLoose40 <= 5 && fMT2tree->NBJets40CSVM == 1) ssignal = "_3to5j1b";
		if(fMT2tree->NJetsIDLoose40 >= 3 && fMT2tree->NJetsIDLoose40 <= 5 && fMT2tree->NBJets40CSVM == 2) ssignal = "_3to5j2b";
		if(                                                                  fMT2tree->NBJets40CSVM >= 3) ssignal = "_3b";
		if(fMT2tree->NJetsIDLoose40 >= 6 &&                                  fMT2tree->NBJets40CSVM == 0) ssignal = "_6j0b";
		if(fMT2tree->NJetsIDLoose40 >= 6 &&                                  fMT2tree->NBJets40CSVM == 1) ssignal = "_6j1b";
		if(fMT2tree->NJetsIDLoose40 >= 6 &&                                  fMT2tree->NBJets40CSVM == 2) ssignal = "_6j2b";

		//get BTV weight
		double btagSF(1.), btagSFerr(0.);
		if(!fMT2tree->misc.isData && fbTagReweight){
		if(fMT2tree->NJetsIDLoose40 == 2 &&                                  fMT2tree->NBJets40CSVM == 0) { 
			btagSF = fMT2tree->SFWeight.BTagCSV40eq0; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq0Error; }
		if(fMT2tree->NJetsIDLoose40 == 2 &&                                  fMT2tree->NBJets40CSVM >= 1) { 
			btagSF = fMT2tree->SFWeight.BTagCSV40ge1; btagSFerr = fMT2tree->SFWeight.BTagCSV40ge1Error; }
		if(fMT2tree->NJetsIDLoose40 >= 3 && fMT2tree->NJetsIDLoose40 <= 5 && fMT2tree->NBJets40CSVM == 0) { 
			btagSF = fMT2tree->SFWeight.BTagCSV40eq0; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq0Error; }
		if(fMT2tree->NJetsIDLoose40 >= 3 && fMT2tree->NJetsIDLoose40 <= 5 && fMT2tree->NBJets40CSVM == 1) { 
			btagSF = fMT2tree->SFWeight.BTagCSV40eq1; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq1Error; }
		if(fMT2tree->NJetsIDLoose40 >= 3 && fMT2tree->NJetsIDLoose40 <= 5 && fMT2tree->NBJets40CSVM == 2) { 
			btagSF = fMT2tree->SFWeight.BTagCSV40eq2; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq2Error; }
		if(                                                                  fMT2tree->NBJets40CSVM >= 3) { 
			btagSF = fMT2tree->SFWeight.BTagCSV40ge3; btagSFerr = fMT2tree->SFWeight.BTagCSV40ge3Error; }
		if(fMT2tree->NJetsIDLoose40 >= 6 &&                                  fMT2tree->NBJets40CSVM == 0) { 
			btagSF = fMT2tree->SFWeight.BTagCSV40eq0; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq0Error; }
		if(fMT2tree->NJetsIDLoose40 >= 6 &&                                  fMT2tree->NBJets40CSVM == 1) { 
			btagSF = fMT2tree->SFWeight.BTagCSV40eq1; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq1Error; }
		if(fMT2tree->NJetsIDLoose40 >= 6 &&                                  fMT2tree->NBJets40CSVM == 2) { 
			btagSF = fMT2tree->SFWeight.BTagCSV40eq2; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq2Error; }
			weight = weight * btagSF;
			weightnoISR = weightnoISR * btagSF;
		}

		string hh = slep + sHT + ssignal + string("_") + sampletype;
		string hhnolep  =  sHT + ssignal + string("_") + sampletype;//needed for genlept filling

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

		if(sampletype=="QCD" || sampletype=="Other" || sampletype=="data" || fMT2tree->misc.isData) continue;

		int genele = fMT2tree->GenNumLeptFromW(11, 0, 1000, fIncludeTaus);
		int genmuo = fMT2tree->GenNumLeptFromW(13, 0, 1000, fIncludeTaus);
		int genlepemu = fMT2tree->GenNumLeptFromW(1113,0,1000,  fIncludeTaus);
		int genemuWtau  = fMT2tree->GenNumLeptFromW(1113,0,1000,  true);
		int genemuNotau = fMT2tree->GenNumLeptFromW(1113,0,1000,  false);
		int genemutau = genemuWtau-genemuNotau;
		int gentau = fMT2tree->GenNumLeptFromW(16,  0,1000, false);
		int genhadtau = gentau-genemutau;
		int genlep = fMT2tree->GenNumLeptFromW(1113,0,1000, fIncludeTaus);//genlep is genele+genmuo

               if(recoedele && (genele==0||genlep>2)){
                    	histos[(string)"LeptonEventsBG"   + hh]->Fill(fMT2tree->misc.MT2, weight);//bg, lepton
                    	histos[(string)"noISR_LeptonEventsBG"   + hh]->Fill(fMT2tree->misc.MT2, weightnoISR);//bg, lepton
		} if(recoedmuo && (genmuo==0||genlep>2)){
                    	histos[(string)"LeptonEventsBG"   + hh]->Fill(fMT2tree->misc.MT2, weight);//bg, lepton
                    	histos[(string)"noISR_LeptonEventsBG"   + hh]->Fill(fMT2tree->misc.MT2, weightnoISR);//bg, lepton
		} if(recoedele && genlep==2){
                    	histos[(string)"LeptonEvents_doubleLL"   + hh]->Fill(fMT2tree->misc.MT2, weight);//bg, lepton
                    	histos[(string)"noISR_LeptonEvents_doubleLL"   + hh]->Fill(fMT2tree->misc.MT2, weightnoISR);//bg, lepton
		} if(recoedmuo && genlep==2){
                    	histos[(string)"LeptonEvents_doubleLL"   + hh]->Fill(fMT2tree->misc.MT2, weight);//bg, lepton
                    	histos[(string)"noISR_LeptonEvents_doubleLL"   + hh]->Fill(fMT2tree->misc.MT2, weightnoISR);//bg, lepton
		}
		if(fbTagError && !fMT2tree->misc.isData && btagSF!=0){// btagSFerr=0 for data
			if(recoedele && (genele==0||genlep>2)){
				histos[(string)"LeptonEventsBGUp"   + hh]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF+btagSFerr));
				histos[(string)"noISR_LeptonEventsBGUp"   + hh]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF+btagSFerr));
			} if(recoedmuo && (genmuo==0||genlep>2)){
				histos[(string)"LeptonEventsBGUp"   + hh]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF+btagSFerr));
				histos[(string)"noISR_LeptonEventsBGUp"   + hh]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF+btagSFerr));
			} if(recoedele && (genele==0||genlep>2)){
				histos[(string)"LeptonEventsBGDown"   + hh]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF-btagSFerr));
				histos[(string)"noISR_LeptonEventsBGDown"   + hh]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF-btagSFerr));
			} if(recoedmuo && (genmuo==0||genlep>2)){
				histos[(string)"LeptonEventsBGDown"   + hh]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF-btagSFerr));
				histos[(string)"noISR_LeptonEventsBGDown"   + hh]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF-btagSFerr));
			} if(recoedele && genlep==2){
				histos[(string)"LeptonEvents_doubleLL_SFup"   + hh]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF+btagSFerr));
				histos[(string)"noISR_LeptonEvents_doubleLL_SFup"   + hh]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF+btagSFerr));
			} if(recoedmuo && genlep==2){
				histos[(string)"LeptonEvents_doubleLL_SFup"   + hh]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF+btagSFerr));
				histos[(string)"noISR_LeptonEvents_doubleLL_SFup"   + hh]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF+btagSFerr));
			} if(recoedele && genlep==2){
				histos[(string)"LeptonEvents_doubleLL_SFdown" + hh]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF-btagSFerr));
				histos[(string)"noISR_LeptonEvents_doubleLL_SFdown" + hh]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF-btagSFerr));
			} if(recoedmuo && genlep==2){
				histos[(string)"LeptonEvents_doubleLL_SFdown" + hh]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF-btagSFerr));
				histos[(string)"noISR_LeptonEvents_doubleLL_SFdown" + hh]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF-btagSFerr));
			}
		}

	//shistos for efficiencies
	bool leptfoundnomt_ele = false;
	bool leptfound_ele     = false;
	bool eventgood_ele     = false;
	bool acceptance_ele    = false;
	bool leptfoundnomt_muo = false;
	bool leptfound_muo     = false;
	bool eventgood_muo     = false;
	bool acceptance_muo    = false;

		if(fUsedoubleLL){//double lost lepton counting not harmful due to ratio
		if(fMT2tree->GenLeptFromW(11, 0 , 1000,fIncludeTaus)==true&&(genlep==1||genlep==2)) eventgood_ele     =true;
		if(fMT2tree->GenLeptFromW(13, 0 , 1000,fIncludeTaus)==true&&(genlep==1||genlep==2)) eventgood_muo     =true;
		if(fMT2tree->GenLeptFromW(11,10,  2.4 ,fIncludeTaus)==true&&(genlep==1||genlep==2)) acceptance_ele    =true;
		if(fMT2tree->GenLeptFromW(13,10,  2.4 ,fIncludeTaus)==true&&(genlep==1||genlep==2)) acceptance_muo    =true;
		} else {
		if(fMT2tree->GenLeptFromW(11, 0 , 1000,fIncludeTaus)==true&&genlep==1)              eventgood_ele     =true;
		if(fMT2tree->GenLeptFromW(13, 0 , 1000,fIncludeTaus)==true&&genlep==1)              eventgood_muo     =true;
		if(fMT2tree->GenLeptFromW(11,10,  2.4 ,fIncludeTaus)==true&&genlep==1)              acceptance_ele    =true;
		if(fMT2tree->GenLeptFromW(13,10,  2.4 ,fIncludeTaus)==true&&genlep==1)              acceptance_muo    =true;
		}
		if(recoedele==true && eventgood_ele==true)                                          leptfound_ele     =true;
		if(recoedmuo==true && eventgood_muo==true)                                          leptfound_muo     =true;
		if(recoedelenomt==true && eventgood_ele==true)                                      leptfoundnomt_ele =true;
		if(recoedmuonomt==true && eventgood_muo==true)                                      leptfoundnomt_muo =true;


	if(eventgood_ele)  histos[(string)"WEvents"          + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);
	if(eventgood_muo)  histos[(string)"WEvents"          + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);
	if(acceptance_ele) histos[(string)"WEventsAcc"       + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);
	if(acceptance_muo) histos[(string)"WEventsAcc"       + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);
	if(leptfound_ele)  histos[(string)"WEventsReco"      + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);
	if(leptfound_muo)  histos[(string)"WEventsReco"      + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);
	if(leptfoundnomt_ele)  histos[(string)"WEventsRecoNoMT"      + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);
	if(leptfoundnomt_muo)  histos[(string)"WEventsRecoNoMT"      + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);
	if(eventgood_ele&&genlep>1)  histos[(string)"WEvents_dLL"     + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);
	if(eventgood_muo&&genlep>1)  histos[(string)"WEvents_dLL"     + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);
	if(leptfound_ele&&genlep>1)  histos[(string)"WEventsReco_dLL" + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);
	if(leptfound_muo&&genlep>1)  histos[(string)"WEventsReco_dLL" + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);

	if(eventgood_ele)  histos[(string)"noISR_WEvents"          + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);
	if(eventgood_muo)  histos[(string)"noISR_WEvents"          + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);
	if(acceptance_ele) histos[(string)"noISR_WEventsAcc"       + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);
	if(acceptance_muo) histos[(string)"noISR_WEventsAcc"       + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);
	if(leptfound_ele)  histos[(string)"noISR_WEventsReco"      + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);
	if(leptfound_muo)  histos[(string)"noISR_WEventsReco"      + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);
	if(leptfoundnomt_ele)  histos[(string)"noISR_WEventsRecoNoMT"      + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);
	if(leptfoundnomt_muo)  histos[(string)"noISR_WEventsRecoNoMT"      + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);
	if(eventgood_ele&&genlep>1)  histos[(string)"noISR_WEvents_dLL"     + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);
	if(eventgood_muo&&genlep>1)  histos[(string)"noISR_WEvents_dLL"     + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);
	if(leptfound_ele&&genlep>1)  histos[(string)"noISR_WEventsReco_dLL" + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);
	if(leptfound_muo&&genlep>1)  histos[(string)"noISR_WEventsReco_dLL" + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);

	if(norecolep && (genlep==1&&fMT2tree->GenLeptFromW(11, 0, 1000, fIncludeTaus)==1)){
		histos[(string)"NoLeptonEvents" + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);//no bg (signal), no lepton
		histos[(string)"noISR_NoLeptonEvents" + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);//no bg (signal), no lepton
	} if(norecolep && (genlep==1&&fMT2tree->GenLeptFromW(13, 0, 1000, fIncludeTaus)==1)){
		histos[(string)"NoLeptonEvents" + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);//no bg (signal), no lepton
		histos[(string)"noISR_NoLeptonEvents" + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);//no bg (signal), no lepton
	} if(norecolep && (genlep==2&&fMT2tree->GenLeptFromW(11, 0, 1000, fIncludeTaus)==1) && !(genlep==2&&fMT2tree->GenLeptFromW(13, 0, 1000, fIncludeTaus)==1)){
		histos[(string)"NoLeptonEvents_doubleLL" + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);//no bg (signal), no lepton
		histos[(string)"noISR_NoLeptonEvents_doubleLL" + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);//no bg (signal), no lepton
	} if(norecolep && (genlep==2&&fMT2tree->GenLeptFromW(13, 0, 1000, fIncludeTaus)==1) && !(genlep==2&&fMT2tree->GenLeptFromW(11, 0, 1000, fIncludeTaus)==1)){
		histos[(string)"NoLeptonEvents_doubleLL" + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weight);//no bg (signal), no lepton
		histos[(string)"noISR_NoLeptonEvents_doubleLL" + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR);//no bg (signal), no lepton
	} if(norecolep && (genlep==2&&fMT2tree->GenLeptFromW(13, 0, 1000, fIncludeTaus)==1&&fMT2tree->GenLeptFromW(11, 0, 1000, fIncludeTaus)==1)){
		histos[(string)"NoLeptonEvents_doubleLL" + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, 0.5*weight);//no bg (signal), no lepton
		histos[(string)"NoLeptonEvents_doubleLL" + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, 0.5*weight);//no bg (signal), no lepton
		histos[(string)"noISR_NoLeptonEvents_doubleLL" + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, 0.5*weightnoISR);//no bg (signal), no lepton
		histos[(string)"noISR_NoLeptonEvents_doubleLL" + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, 0.5*weightnoISR);//no bg (signal), no lepton
	}
	if(fbTagError && !fMT2tree->misc.isData && btagSF!=0){
		if(norecolep && (genlep==1&&fMT2tree->GenLeptFromW(11, 0, 1000, fIncludeTaus)==1)){
			histos[(string)"NoLeptonEvents_SFup" + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF+btagSFerr));
			histos[(string)"noISR_NoLeptonEvents_SFup" + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF+btagSFerr));
		} if(norecolep && (genlep==1&&fMT2tree->GenLeptFromW(13, 0, 1000, fIncludeTaus)==1)){
			histos[(string)"NoLeptonEvents_SFup" + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF+btagSFerr));
			histos[(string)"noISR_NoLeptonEvents_SFup" + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF+btagSFerr));
		} if(norecolep && (genlep==2&&fMT2tree->GenLeptFromW(11, 0, 1000, fIncludeTaus)==1) && !(genlep==2&&fMT2tree->GenLeptFromW(13, 0, 1000, fIncludeTaus)==1)){
			histos[(string)"NoLeptonEvents_doubleLL_SFup_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF+btagSFerr));
			histos[(string)"noISR_NoLeptonEvents_doubleLL_SFup_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF+btagSFerr));
		} if(norecolep && (genlep==2&&fMT2tree->GenLeptFromW(13, 0, 1000, fIncludeTaus)==1) && !(genlep==2&&fMT2tree->GenLeptFromW(11, 0, 1000, fIncludeTaus)==1)){
			histos[(string)"NoLeptonEvents_doubleLL_SFup_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF+btagSFerr));
			histos[(string)"noISR_NoLeptonEvents_doubleLL_SFup_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF+btagSFerr));
		} if(norecolep && (genlep==2&&fMT2tree->GenLeptFromW(13, 0, 1000, fIncludeTaus)==1&&fMT2tree->GenLeptFromW(11, 0, 1000, fIncludeTaus)==1)){
			histos[(string)"NoLeptonEvents_doubleLL_SFup_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, 0.5*weight/btagSF*(btagSF+btagSFerr));
			histos[(string)"NoLeptonEvents_doubleLL_SFup_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, 0.5*weight/btagSF*(btagSF+btagSFerr));
			histos[(string)"noISR_NoLeptonEvents_doubleLL_SFup_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, 0.5*weightnoISR/btagSF*(btagSF+btagSFerr));
			histos[(string)"noISR_NoLeptonEvents_doubleLL_SFup_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, 0.5*weightnoISR/btagSF*(btagSF+btagSFerr));
		}
		if(norecolep && (genlep==1&&fMT2tree->GenLeptFromW(11, 0, 1000, fIncludeTaus)==1)){
			histos[(string)"NoLeptonEvents_SFdown" + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF-btagSFerr));
			histos[(string)"noISR_NoLeptonEvents_SFdown" + (string)"_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF-btagSFerr));
		} if(norecolep && (genlep==1&&fMT2tree->GenLeptFromW(13, 0, 1000, fIncludeTaus)==1)){
			histos[(string)"NoLeptonEvents_SFdown" + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF-btagSFerr));
			histos[(string)"noISR_NoLeptonEvents_SFdown" + (string)"_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF-btagSFerr));
		} if(norecolep && (genlep==2&&fMT2tree->GenLeptFromW(11, 0, 1000, fIncludeTaus)==1) && !(genlep==2&&fMT2tree->GenLeptFromW(13, 0, 1000, fIncludeTaus)==1)){
			histos[(string)"NoLeptonEvents_doubleLL_SFdown_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF-btagSFerr));
			histos[(string)"noISR_NoLeptonEvents_doubleLL_SFdown_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF-btagSFerr));
		} if(norecolep && (genlep==2&&fMT2tree->GenLeptFromW(13, 0, 1000, fIncludeTaus)==1) && !(genlep==2&&fMT2tree->GenLeptFromW(11, 0, 1000, fIncludeTaus)==1)){
			histos[(string)"NoLeptonEvents_doubleLL_SFdown_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weight/btagSF*(btagSF-btagSFerr));
			histos[(string)"noISR_NoLeptonEvents_doubleLL_SFdown_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, weightnoISR/btagSF*(btagSF-btagSFerr));
		} if(norecolep && (genlep==2&&fMT2tree->GenLeptFromW(13, 0, 1000, fIncludeTaus)==1&&fMT2tree->GenLeptFromW(11, 0, 1000, fIncludeTaus)==1)){
			histos[(string)"NoLeptonEvents_doubleLL_SFdown_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, 0.5*weight/btagSF*(btagSF-btagSFerr));
			histos[(string)"NoLeptonEvents_doubleLL_SFdown_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, 0.5*weight/btagSF*(btagSF-btagSFerr));
			histos[(string)"noISR_NoLeptonEvents_doubleLL_SFdown_Ele" + hhnolep]->Fill(fMT2tree->misc.MT2, 0.5*weightnoISR/btagSF*(btagSF-btagSFerr));
			histos[(string)"noISR_NoLeptonEvents_doubleLL_SFdown_Muo" + hhnolep]->Fill(fMT2tree->misc.MT2, 0.5*weightnoISR/btagSF*(btagSF-btagSFerr));
		}
	}

	} //while(myEvtList->GetEntry(counter++) !=-1)
	delete fMT2tree;
	delete fSamples[i].tree;

	} //for(size_t i = 0; i < fSamples.size(); ++i)

	//add overflow to last bin
	cout << "add overflow to last bin -> does not work for TEff" << endl;
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->SetBinContent(h->second->GetNbinsX(),
					 h->second->GetBinContent(h->second->GetNbinsX()  )+ 
					 h->second->GetBinContent(h->second->GetNbinsX()+1)  );
		h->second->SetBinError(  h->second->GetNbinsX(),
					 sqrt(h->second->GetBinError(h->second->GetNbinsX()  )*
					      h->second->GetBinError(h->second->GetNbinsX()  )+
					      h->second->GetBinError(h->second->GetNbinsX()+1)*
					      h->second->GetBinError(h->second->GetNbinsX()+1)  ));
	}

	cout << "add all samples to mc, etc." << endl;
	//add TTbar and SingleTop to Top
	//add WJets and Top to WandTop
	//add rest to noWandTop
	//add noWandTop and WandTop to mc
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTbinsize;        ++i3){
	for(int i4 = 0; i4<leptontypesize;   ++i4){
		if(fMET && i3==1) continue;
		string hsq   = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_QCD");
		string hsw   = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_WJets");
		string hsz   = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_ZJets");
		string hstt  = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_TTbar");
		string hsst  = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_SingleTop");
		string hst   = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_Top");
		string hswt  = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_WandTop");
		string hsnwt = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_noWandTop");
		string hso   = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_Other");
		string hsmc  = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_mc");
		for(unsigned int i0 = 0; i0<histonames.size();++i0){
			histos[(histonames[i0]+hsnwt)]->Add(histos[(histonames[i0]+hsq)],  1);//noWandTop
			histos[(histonames[i0]+hsnwt)]->Add(histos[(histonames[i0]+hsz)],  1);
			histos[(histonames[i0]+hsnwt)]->Add(histos[(histonames[i0]+hso)],  1);
			histos[(histonames[i0]+hst)  ]->Add(histos[(histonames[i0]+hsst)], 1);//Top
			histos[(histonames[i0]+hst)  ]->Add(histos[(histonames[i0]+hstt)], 1);
			histos[(histonames[i0]+hswt) ]->Add(histos[(histonames[i0]+hst)],  1);//WandTop
			histos[(histonames[i0]+hswt) ]->Add(histos[(histonames[i0]+hsw)],  1);
			histos[(histonames[i0]+hsmc) ]->Add(histos[(histonames[i0]+hswt)], 1);//mc
			histos[(histonames[i0]+hsmc) ]->Add(histos[(histonames[i0]+hsnwt)],1);
			string noisr = "noISR_";
			histos[noisr + (histonames[i0]+hsnwt)]->Add(histos[noisr + (histonames[i0]+hsq)],  1);//noWandTop
			histos[noisr + (histonames[i0]+hsnwt)]->Add(histos[noisr + (histonames[i0]+hsz)],  1);
			histos[noisr + (histonames[i0]+hsnwt)]->Add(histos[noisr + (histonames[i0]+hso)],  1);
			histos[noisr + (histonames[i0]+hst)  ]->Add(histos[noisr + (histonames[i0]+hsst)], 1);//Top
			histos[noisr + (histonames[i0]+hst)  ]->Add(histos[noisr + (histonames[i0]+hstt)], 1);
			histos[noisr + (histonames[i0]+hswt) ]->Add(histos[noisr + (histonames[i0]+hst)],  1);//WandTop
			histos[noisr + (histonames[i0]+hswt) ]->Add(histos[noisr + (histonames[i0]+hsw)],  1);
			histos[noisr + (histonames[i0]+hsmc) ]->Add(histos[noisr + (histonames[i0]+hswt)], 1);//mc
			histos[noisr + (histonames[i0]+hsmc) ]->Add(histos[noisr + (histonames[i0]+hsnwt)],1);
		}
	}}}

	cout << "finalize efficiencies" << endl;
	//now create efficiencies histograms and tefficiencies directly from histograms;
	teffnames.push_back("TEff_MTEfficiency");
	teffnames.push_back("TEff_MTEfficiency2");
	teffnames.push_back("TEff_WEfficiencyRecoVsAll");
	teffnames.push_back("TEff_WEfficiencyRecoVsAcc");
	teffnames.push_back("TEff_WEfficiencyRecoVsAllNoMT");
	teffnames.push_back("TEff_WEfficiencyRecoVsAccNoMT");
	teffnames.push_back("TEff_WEfficiencyAccVsAll");
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTbinsize;        ++i3){
	for(int i4 = 0; i4<leptontypesize;   ++i4){
		if(fMET && i3==1) continue;
		string hs   = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_") + sample_type[i1];
		string mapname;
		//don't need no isr as flag does already the job
		TH1D *passRecoVsAll     = (TH1D*)histos["WEventsReco"     +hs]->Clone("passRecoVsAll");
		TH1D  *totRecoVsAll     = (TH1D*)histos["WEvents"         +hs]->Clone( "totRecoVsAll");
		TH1D *passRecoVsAcc     = (TH1D*)histos["WEventsReco"     +hs]->Clone("passRecoVsAcc");
		TH1D  *totRecoVsAcc     = (TH1D*)histos["WEventsAcc"      +hs]->Clone( "totRecoVsAcc");
		TH1D *passRecoVsAllNoMT = (TH1D*)histos["WEventsRecoNoMT" +hs]->Clone("passRecoVsAllNoMT");
		TH1D  *totRecoVsAllNoMT = (TH1D*)histos["WEvents"         +hs]->Clone( "totRecoVsAllNoMT");
		TH1D *passRecoVsAccNoMT = (TH1D*)histos["WEventsRecoNoMT" +hs]->Clone("passRecoVsAccNoMT");
		TH1D  *totRecoVsAccNoMT = (TH1D*)histos["WEventsAcc"      +hs]->Clone( "totRecoVsAccNoMT");
		TH1D *passAccVsAll      = (TH1D*)histos["WEventsAcc"      +hs]->Clone("passAccVsAll" );
		TH1D  *totAccVsAll      = (TH1D*)histos["WEvents"         +hs]->Clone( "totAccVsAll" );
		TH1D *passMT            = (TH1D*)histos["LeptonEvents"    +hs]->Clone("passMT");
		TH1D  *totMT            = (TH1D*)histos["LeptonEventsNoMT"+hs]->Clone("totMT" );
		TH1D *passMT2           = (TH1D*)histos["WEventsReco"     +hs]->Clone("passMT2");
		TH1D  *totMT2           =(TH1D*)histos["WEventsRecoNoMT"  +hs]->Clone("totMT2" );
		mapname = "TEff_WEfficiencyRecoVsAllNoMT";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passRecoVsAllNoMT), (*totRecoVsAllNoMT));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_WEfficiencyRecoVsAccNoMT";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passRecoVsAccNoMT), (*totRecoVsAccNoMT));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_WEfficiencyRecoVsAll";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passRecoVsAll), (*totRecoVsAll));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_WEfficiencyRecoVsAcc";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passRecoVsAcc), (*totRecoVsAcc));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_WEfficiencyAccVsAll";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passAccVsAll),  (*totAccVsAll) );
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_MTEfficiency";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passMT),  (*totMT) );
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_MTEfficiency2";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passMT2),  (*totMT2) );
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
	}}}}

	cout << "Saving." << endl;
	if(!fISRreweight) outputname = "NoISR_" + outputname;
    	TFile *fsavefile = new TFile(outputdir + outputname,"RECREATE");
	fsavefile->cd();
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Write();
	}
	for(map<string,TEfficiency*>::iterator h=teff.begin(); h!=teff.end();++h){
		h->second->Write();
	}
	fsavefile->Close();
	cout << "Saved histograms in " << outputdir << outputname << endl;

	//do estimation
	cout << "do estimation" << endl;

	if(fRebin){//do not estimate along MT2
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTbinsize;        ++i3){
	for(int i4 = 0; i4<leptontypesize;   ++i4){
		if(fMET && i3==1) continue;//MET triggers trigger only one HT bin
		string hs   = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_") + sample_type[i1];
		string mapname;
		double startbin = histos["RecoLepEvents"+hs]->GetBinLowEdge(1);
		double endbin = histos["RecoLepEvents"+hs]->GetBinLowEdge(histos["RecoLepEvents"+hs]->GetNbinsX())+histos["RecoLepEvents"+hs]->GetBinWidth(histos["RecoLepEvents"+hs]->GetNbinsX());
		double rebin[2] = {startbin, endbin};
		int nrebin = 1;
		histos["RecoLepEvents"+hs]                  = RebinThisHistogram(histos["RecoLepEvents"+hs], nrebin, rebin);
		histos["LeptonEvents"+hs]                   = RebinThisHistogram(histos["LeptonEvents"+hs], nrebin, rebin);
		histos["LeptonEventsNoMT"+hs]               = RebinThisHistogram(histos["LeptonEventsNoMT"+hs], nrebin, rebin);
		histos["LeptonEventsBG"+hs]                 = RebinThisHistogram(histos["LeptonEventsBG"+hs], nrebin, rebin);
		histos["LeptonEvents_doubleLL"+hs]          = RebinThisHistogram(histos["LeptonEvents_doubleLL"+hs], nrebin, rebin);
		histos["NoLeptonEvents"+hs]                 = RebinThisHistogram(histos["NoLeptonEvents"+hs], nrebin, rebin);
		histos["NoLeptonEvents_SFup"+hs]            = RebinThisHistogram(histos["NoLeptonEvents_SFup"+hs], nrebin, rebin);
		histos["NoLeptonEvents_SFdown"+hs]          = RebinThisHistogram(histos["NoLeptonEvents_SFdown"+hs], nrebin, rebin);
		histos["NoLeptonEvents_doubleLL"+hs]        = RebinThisHistogram(histos["NoLeptonEvents_doubleLL"+hs], nrebin, rebin);
		histos["NoLeptonEvents_doubleLL_SFup"+hs]   = RebinThisHistogram(histos["NoLeptonEvents_doubleLL_SFup"+hs], nrebin, rebin);
		histos["NoLeptonEvents_doubleLL_SFdown"+hs] = RebinThisHistogram(histos["NoLeptonEvents_doubleLL_SFdown"+hs], nrebin, rebin);
		histos["WEvents"+hs]                        = RebinThisHistogram(histos["WEvents"+hs], nrebin, rebin);
		histos["WEventsAcc"+hs]                     = RebinThisHistogram(histos["WEventsAcc"+hs], nrebin, rebin);
		histos["WEventsReco"+hs]                    = RebinThisHistogram(histos["WEventsReco"+hs], nrebin, rebin);
		histos["WEvents_dLL"+hs]                    = RebinThisHistogram(histos["WEvents_dLL"+hs], nrebin, rebin);
		histos["WEventsReco_dLL"+hs]                = RebinThisHistogram(histos["WEventsReco_dLL"+hs], nrebin, rebin);
		histos["LeptonEventsBGUp"+hs]               = RebinThisHistogram(histos["LeptonEventsBGUp"+hs], nrebin, rebin);
		histos["RecoLepEventsUp"+hs]                = RebinThisHistogram(histos["RecoLepEventsUp"+hs], nrebin, rebin);
		histos["LeptonEventsBGDown"+hs]             = RebinThisHistogram(histos["LeptonEventsBGDown"+hs], nrebin, rebin);
		histos["RecoLepEventsDown"+hs]              = RebinThisHistogram(histos["RecoLepEventsDown"+hs], nrebin, rebin);
		histos["LeptonEvents_doubleLL_SFup"+hs]     = RebinThisHistogram(histos["LeptonEvents_doubleLL_SFup"+hs], nrebin, rebin);
		histos["LeptonEvents_doubleLL_SFdown"+hs]   = RebinThisHistogram(histos["LeptonEvents_doubleLL_SFdown"+hs], nrebin, rebin);
		histos["WEventsRecoNoMT"+hs]                = RebinThisHistogram(histos["WEventsRecoNoMT"+hs], nrebin, rebin);
		teff["TEff_MTEfficiency"+hs]                = RebinThisEfficiency(teff["TEff_MTEfficiency"+hs], nrebin, rebin);
		teff["TEff_WEfficiencyRecoVsAll"+hs]        = RebinThisEfficiency(teff["TEff_WEfficiencyRecoVsAll"+hs], nrebin, rebin);
		teff["TEff_WEfficiencyRecoVsAcc"+hs]        = RebinThisEfficiency(teff["TEff_WEfficiencyRecoVsAcc"+hs], nrebin, rebin);
		teff["TEff_WEfficiencyAccVsAll"+hs]         = RebinThisEfficiency(teff["TEff_WEfficiencyAccVsAll"+hs], nrebin, rebin);
		teff["TEff_WEfficiencyRecoVsAllNoMT"+hs]    = RebinThisEfficiency(teff["TEff_WEfficiencyRecoVsAllNoMT"+hs], nrebin, rebin);
		teff["TEff_WEfficiencyRecoVsAccNoMT"+hs]    = RebinThisEfficiency(teff["TEff_WEfficiencyRecoVsAccNoMT"+hs], nrebin, rebin);
		teff["TEff_MTEfficiency2"+hs]               = RebinThisEfficiency(teff["TEff_MTEfficiency2"+hs], nrebin, rebin);

	}}}}
	}

	//void GetLostLeptonEstimate(histograms,tefficiencies,rel_sys_uncert,rel_sys_uncert_bg,version,PrintSummaryTable,PrintYieldTable,PrintPredictionCard,makeFullPrintout=false);
	//versions for SummaryTable
	//version == 0,1: Summary table contains also MCPred (i.e. method applied on simulation)
	//version == 2,3: Summary table contains only data prediction and MC truth, but not MCPred
	//version == 0,2: Summary table contains the factor R_LL (i.e. 1-e / e) instead of e along (e=efficiency)
	//version == 1,3: Summary table contains the efficiency e
	GetLostLeptonEstimate(histos, teff, frel_sys_uncert,  frel_sys_uncert_bg, 0, true, true, true, true, true, true);

	if(fWriteToFile && fAppend){
		TString logname =outputdir + "results.log"; 
		ofstream f_log (logname.Data(), ios::app);
		f_log << fLogStream->str();
		cout << "wrote results into  " << logname << " (appended at the end of old file)" << endl;
	}else if(fWriteToFile){
		TString logname =outputdir + "results.log"; 
		ofstream f_log (logname.Data(), ios::trunc);
		f_log << fLogStream->str();
		cout << "wrote results into  " << logname <<  " (old file replaced)" << endl;
	} else{
		cout << fLogStream->str();
	}
	delete fLogStream;

}//void LostLeptonEstimate()

TH1D* RebinThisHistogram(TH1D* histogram, int nrebinnedbins, double* rebinnedbins){

	string name = histogram->GetName();
	histogram->SetName((name + "_original").c_str());
	TH1D* temphists = (TH1D*)histogram->Rebin(nrebinnedbins, (name + "_rebinned").c_str(), rebinnedbins);
	temphists->SetName((name).c_str());
	return temphists;

}

TEfficiency* RebinThisEfficiency(TEfficiency* tefficiency, int nrebinnedbins, double* rebinnedbins){

	string name = tefficiency->GetName();
	TH1D *pass = (TH1D*)tefficiency->GetCopyPassedHisto(); TH1D *total = (TH1D*)tefficiency->GetCopyTotalHisto();
	name = pass->GetName();
	TH1D *passrebinned = (TH1D*)pass->Rebin(nrebinnedbins, (name + "_rebinned").c_str(), rebinnedbins);
	name = total->GetName();
	TH1D *totalrebinned = (TH1D*)total->Rebin(nrebinnedbins, (name + "_rebinned").c_str(), rebinnedbins);
	name = tefficiency->GetName();
	tefficiency->SetNameTitle((name+"_original").c_str(), (name+"_original").c_str());
	TEfficiency *temp = new TEfficiency((*passrebinned), (*totalrebinned));
	temp->SetNameTitle((name+"_temp").c_str(), (name).c_str());
	TEfficiency *clone = (TEfficiency*)temp->Clone((name).c_str());
	delete temp;
	return clone;
}


// this function uses all histograms and efficiencies to compute the amount of lost lepton
//NOTE: You should implement a way to store the results into a histogram using the vectors just defined in the beginning - to make the hardcoded stuff of LLVisualization.C obsolete.
void GetLostLeptonEstimate(map<string, TH1D*> histograms, map<string, TEfficiency*> tefficiencies, double rel_sys_uncert,  double rel_sys_uncert_bg, int version, Bool_t PrintSummaryTable, Bool_t PrintYieldTable, Bool_t PrintPredictionCard, Bool_t makeFullPrintout, Bool_t SavePrediction, Bool_t PlotPrediction){

	//copy the imput variables
	map<string, TH1D*> hists = histograms;
	map<string, TEfficiency*> teffs = tefficiencies;

	//these are all the numbers we want to have in a table later
	vector<int> sr; sr.clear();//signal region
	vector<int> htr; htr.clear();//ht region
	vector<int> lepr; lepr.clear();//lepton type
	vector<double> MT2low; MT2low.clear();//lower bound of bin
	vector<double> MT2up; MT2up.clear();//upper bound of bin // for final bin store 10000.
	vector<int> MT2bin; MT2bin.clear();//store only binnumber --> used for datacard
	//printoutnumbers - store numbers which might be used in prinouts
	vector<double> mctruth; mctruth.clear();
	vector<double> mctrutherr; mctrutherr.clear();
	vector<double> mctruthfractionW;  mctruthfractionW.clear();//new
	vector<double> mctruthfractionTT; mctruthfractionTT.clear();//new
	vector<double> datapred; datapred.clear();
	vector<double> datapred_stat_err; datapred_stat_err.clear();
	vector<double> datapred_syst_err; datapred_syst_err.clear();
	vector<double> MCpred; MCpred.clear();
	vector<double> MCpred_stat_err; MCpred_stat_err.clear();
	vector<double> MCpred_syst_err; MCpred_syst_err.clear();
	vector<double> LLeff; LLeff.clear();
	vector<double> LLefferr; LLefferr.clear();
	vector<double> MTeff; MTeff.clear();
	vector<double> MTefferr; MTefferr.clear();
	vector<double> numtrueW; numtrueW.clear();
	vector<double> numtrueW_bg; numtrueW_bg.clear();
	vector<double> numData; numData.clear();
	vector<double> numBG; numBG.clear();
	vector<double> numQCD; numQCD.clear();
	vector<double> numZ; numZ.clear();
	vector<double> numW; numW.clear();
	vector<double> numT; numT.clear();
	vector<double> numTT; numTT.clear();//new
	vector<double> numOther; numOther.clear();
	vector<double> numMC; numMC.clear();
	vector<double> BGerr; BGerr.clear();

	double dummy;//this variable is temporary, but needed
	//start printoutloop;

	//maybe you want to change the order of saving, i.e. the order of the three for loops, play with it if you like
	for(int i4 = 0; i4<leptontypesize;   ++i4){
	for(int i3 = 0; i3<HTbinsize;        ++i3){
	for(int i2 = 0; i2<signalregionsize; ++i2){
		if(fMET && i3==1) continue;
	//as all histograms have the same binning just take first histogram 
	string hshl = string("_") + lepton_type[i4] + string("_") + HT_bin[i3] + string("_") + signal_region[i2];
	for(int nx = 1; nx<=hists[string("RecoLepEvents"+hshl+"_mc")]->GetNbinsX(); ++nx){
		//first push_back the signal region information
		sr.push_back(i2);
		if(fMET) htr.push_back(i3);
		else     htr.push_back(i3+1);
		lepr.push_back(i4);
		MT2low.push_back(hists[string("RecoLepEvents"+hshl+"_mc")]->GetBinLowEdge(nx) );
		if(nx!=hists[string("RecoLepEvents"+hshl+"_mc")]->GetNbinsX())
		   MT2up.push_back( hists[string("RecoLepEvents"+hshl+"_mc")]->GetBinLowEdge(nx) + hists[string("RecoLepEvents"+hshl+"_mc")]->GetBinWidth(nx) );
		else
		   MT2up.push_back(10000.);
		MT2bin.push_back(nx);

		//efficiencies
		double W_prob(0.),   W_prob_err(0.);
		double W_acc(0.),    W_acc_err(0.);
		double W_rec(0.),    W_rec_err(0.);
		double Top_prob(0.), Top_prob_err(0.);
		double Top_acc(0.),  Top_acc_err(0.);
		double Top_rec(0.),  Top_rec_err(0.);
		double WT_prob(0.),  WT_prob_err(0.);
		double WT_acc(0.),   WT_acc_err(0.);
		double WT_rec(0.),   WT_rec_err(0.);
		//event yields
		double  nW(0.),  nW_bg(0.),  nW_leptveto(0.),  nW_leptveto_SFup(0.),  nW_leptveto_SFdown(0.),  nW_leptveto_err(0.);
		double nTT(0.), nTT_bg(0.), nTT_leptveto(0.), nTT_leptveto_SFup(0.), nTT_leptveto_SFdown(0.), nTT_leptveto_err(0.);//TTbar
		double nST(0.), nST_bg(0.), nST_leptveto(0.), nST_leptveto_SFup(0.), nST_leptveto_SFdown(0.), nST_leptveto_err(0.);//SingleTop
		double nW_goodevt_dLL(0.), nT_goodevt_dLL(0.), nST_goodevt_dLL(0.),nTT_goodevt_dLL(0.), nWT_goodevt_dLL(0.);
		double nW_goodrecoevt_dLL(0.), nT_goodrecoevt_dLL(0.), nST_goodrecoevt_dLL(0.),nTT_goodrecoevt_dLL(0.), nWT_goodrecoevt_dLL(0.);
		double  nT(0.),  nT_bg(0.),  nT_leptveto(0.),  nT_leptveto_SFup(0.),  nT_leptveto_SFdown(0.),  nT_leptveto_err(0.);//Top
		double nWT(0.), nWT_bg(0.), nWT_leptveto(0.), nWT_leptveto_SFup(0.), nWT_leptveto_SFdown(0.), nWT_leptveto_err(0.);//WandTop
		//doubleLL numbers
		double  nW_dLL(0.),  nW_dLL_SFup(0.),  nW_dLL_SFdown(0.),  nW_dLL_lv(0.),  nW_dLL_lv_SFup(0.),  nW_dLL_lv_SFdown(0.),  nW_dLL_lv_err(0.);
		double nTT_dLL(0.), nTT_dLL_SFup(0.), nTT_dLL_SFdown(0.), nTT_dLL_lv(0.), nTT_dLL_lv_SFup(0.), nTT_dLL_lv_SFdown(0.), nTT_dLL_lv_err(0.);
		double nST_dLL(0.), nST_dLL_SFup(0.), nST_dLL_SFdown(0.), nST_dLL_lv(0.), nST_dLL_lv_SFup(0.), nST_dLL_lv_SFdown(0.), nST_dLL_lv_err(0.);
		double  nT_dLL(0.),  nT_dLL_SFup(0.),  nT_dLL_SFdown(0.),  nT_dLL_lv(0.),  nT_dLL_lv_SFup(0.),  nT_dLL_lv_SFdown(0.),  nT_dLL_lv_err(0.);
		double nWT_dLL(0.), nWT_dLL_SFup(0.), nWT_dLL_SFdown(0.), nWT_dLL_lv(0.), nWT_dLL_lv_SFup(0.), nWT_dLL_lv_SFdown(0.), nWT_dLL_lv_err(0.);
		//event yields - SF weights
		double nW_bg_SFup(0.),   nTT_bg_SFup(0.),   nST_bg_SFup(0.),   nT_bg_SFup(0.),   nWT_bg_SFup(0.);
		double nW_bg_SFdown(0.), nTT_bg_SFdown(0.), nST_bg_SFdown(0.), nT_bg_SFdown(0.), nWT_bg_SFdown(0.);
		//backgrounds
		double QCD_bg(0.), Z_bg(0.), Other_bg(0.), TT_bg(0.), ST_bg(0.), T_bg(0.), W_bg(0.), WT_bg(0.), nonWT_bg(0.);
		//backgrounds - SF weights
		double QCD_bg_SFup(0.),   Z_bg_SFup(0.),   Other_bg_SFup(0.),   T_bg_SFup(0.),   TT_bg_SFup(0.),   ST_bg_SFup(0.),   W_bg_SFup(0.),   WT_bg_SFup(0.),   nonWT_bg_SFup(0.);
		double QCD_bg_SFdown(0.), Z_bg_SFdown(0.), Other_bg_SFdown(0.), T_bg_SFdown(0.), TT_bg_SFdown(0.), ST_bg_SFdown(0.), W_bg_SFdown(0.), WT_bg_SFdown(0.), nonWT_bg_SFdown(0.);
		//data yields, mc yields and tau_leptveto
		double nData(0.);
		double MC_bg(0.), MC_bg_err(0.);
		double mteff(0.), mtefferr(0.);
		double mteffdata(0.), mtefferrdata(0.);
		double  W_evts(0.),  W_evts_acc(0.),  W_evts_reco(0.),  W_evts_nomt_reco(0.);
		double TT_evts(0.), TT_evts_acc(0.), TT_evts_reco(0.), TT_evts_nomt_reco(0.);
		double WT_evts(0.), WT_evts_acc(0.), WT_evts_reco(0.), WT_evts_nomt_reco(0.);
		double mtpass(0.), mttot(0.);
		double W_prob_nomt(0.),   W_prob_nomt_err(0.);
		double W_rec_nomt(0.),    W_rec_nomt_err(0.);
		double Top_prob_nomt(0.), Top_prob_nomt_err(0.);
		double Top_rec_nomt(0.),  Top_rec_nomt_err(0.);
		double WT_prob_nomt(0.),  WT_prob_nomt_err(0.);
		double WT_rec_nomt(0.),   WT_rec_nomt_err(0.);
		double mteff2(0.), mtefferr2(0.);//W+T
		double mteff3(0.), mtefferr3(0.);//MC
		double mtpass2(0.), mttot2(0.);
		double mtpass3(0.), mttot3(0.);

		//the lines below initializes the variables above to the correct value
		W_evts       = hists["WEvents"    +hshl+"_WJets"   ]->GetBinContent(nx);
		W_evts_acc   = hists["WEventsAcc" +hshl+"_WJets"   ]->GetBinContent(nx);
		W_evts_reco  = hists["WEventsReco"+hshl+"_WJets"   ]->GetBinContent(nx);
		TT_evts      = hists["WEvents"    +hshl+"_TTbar"   ]->GetBinContent(nx);
		TT_evts_acc  = hists["WEventsAcc" +hshl+"_TTbar"   ]->GetBinContent(nx);
		TT_evts_reco = hists["WEventsReco"+hshl+"_TTbar"   ]->GetBinContent(nx);
		WT_evts      = hists["WEvents"    +hshl+"_WandTop" ]->GetBinContent(nx);
		WT_evts_acc  = hists["WEventsAcc" +hshl+"_WandTop" ]->GetBinContent(nx);
		WT_evts_reco = hists["WEventsReco"+hshl+"_WandTop" ]->GetBinContent(nx);
		mtpass       = hists["LeptonEvents"    +hshl+"_mc"      ]->GetBinContent(nx);
		mttot        = hists["LeptonEventsNoMT"+hshl+"_mc"      ]->GetBinContent(nx);
		W_evts_nomt_reco  = hists["WEventsRecoNoMT"+hshl+"_WJets"   ]->GetBinContent(nx);
		TT_evts_nomt_reco = hists["WEventsRecoNoMT"+hshl+"_TTbar"   ]->GetBinContent(nx);
		WT_evts_nomt_reco = hists["WEventsRecoNoMT"+hshl+"_WandTop" ]->GetBinContent(nx);
		mtpass3       = hists["WEventsReco"    +hshl+"_mc"      ]->GetBinContent(nx);
		mttot3        = hists["WEventsRecoNoMT"+hshl+"_mc"      ]->GetBinContent(nx);
		mtpass2       = hists["WEventsReco"    +hshl+"_WandTop"      ]->GetBinContent(nx);
		mttot2        = hists["WEventsRecoNoMT"+hshl+"_WandTop"      ]->GetBinContent(nx);

		mteff2        = teffs["TEff_MTEfficiency2"+hshl+"_WandTop"  ]->GetEfficiency(nx);
		mtefferr2     = teffs["TEff_MTEfficiency2"+hshl+"_WandTop"  ]->GetEfficiencyErrorLow(nx);
		dummy         = teffs["TEff_MTEfficiency2"+hshl+"_WandTop"  ]->GetEfficiencyErrorUp(nx);
		if(dummy>mtefferr2  ) mtefferr2   = dummy;
		mteff3        = teffs["TEff_MTEfficiency2"+hshl+"_mc"  ]->GetEfficiency(nx);
		mtefferr3     = teffs["TEff_MTEfficiency2"+hshl+"_mc"  ]->GetEfficiencyErrorLow(nx);
		dummy         = teffs["TEff_MTEfficiency2"+hshl+"_mc"  ]->GetEfficiencyErrorUp(nx);
		if(dummy>mtefferr3  ) mtefferr3   = dummy;

		mteff        = teffs["TEff_MTEfficiency"+hshl+"_mc"  ]->GetEfficiency(nx);
		mtefferr     = teffs["TEff_MTEfficiency"+hshl+"_mc"  ]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_MTEfficiency"+hshl+"_mc"  ]->GetEfficiencyErrorUp(nx);
		if(dummy>mtefferr  ) mtefferr   = dummy;
		mteffdata    = teffs["TEff_MTEfficiency"+hshl+"_data"]->GetEfficiency(nx);
		mtefferrdata = teffs["TEff_MTEfficiency"+hshl+"_data"]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_MTEfficiency"+hshl+"_data"]->GetEfficiencyErrorUp(nx);
		if(dummy>mtefferrdata) mtefferrdata = dummy;

		W_prob       = teffs["TEff_WEfficiencyRecoVsAll"+hshl+"_WJets"]->GetEfficiency(nx);
		W_prob_err   = teffs["TEff_WEfficiencyRecoVsAll"+hshl+"_WJets"]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_WEfficiencyRecoVsAll"+hshl+"_WJets"]->GetEfficiencyErrorUp(nx);
		if(dummy>W_prob_err  ) W_prob_err   = dummy;
		Top_prob     = teffs["TEff_WEfficiencyRecoVsAll"+hshl+"_Top"  ]->GetEfficiency(nx);		//contains SingleTop
		Top_prob_err = teffs["TEff_WEfficiencyRecoVsAll"+hshl+"_Top"  ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy        = teffs["TEff_WEfficiencyRecoVsAll"+hshl+"_Top"  ]->GetEfficiencyErrorUp(nx);	//contains SingleTop
		if(dummy>Top_prob_err) Top_prob_err = dummy;
		W_acc        = teffs["TEff_WEfficiencyAccVsAll" +hshl+"_WJets"]->GetEfficiency(nx);
		W_acc_err    = teffs["TEff_WEfficiencyAccVsAll" +hshl+"_WJets"]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_WEfficiencyAccVsAll" +hshl+"_WJets"]->GetEfficiencyErrorUp(nx);
		if(dummy>W_acc_err   ) W_acc_err    = dummy;
		Top_acc      = teffs["TEff_WEfficiencyAccVsAll" +hshl+"_Top"  ]->GetEfficiency(nx);		//contains SingleTop
		Top_acc_err  = teffs["TEff_WEfficiencyAccVsAll" +hshl+"_Top"  ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy        = teffs["TEff_WEfficiencyAccVsAll" +hshl+"_Top"  ]->GetEfficiencyErrorUp(nx);	//contains SingleTop
		if(dummy>Top_acc_err ) Top_acc_err  = dummy;
		W_rec        = teffs["TEff_WEfficiencyRecoVsAcc"+hshl+"_WJets"]->GetEfficiency(nx);
		W_rec_err    = teffs["TEff_WEfficiencyRecoVsAcc"+hshl+"_WJets"]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_WEfficiencyRecoVsAcc"+hshl+"_WJets"]->GetEfficiencyErrorUp(nx);
		if(dummy>W_rec_err   ) W_rec_err    = dummy;
		Top_rec      = teffs["TEff_WEfficiencyRecoVsAcc"+hshl+"_Top"  ]->GetEfficiency(nx);		//contains SingleTop
		Top_rec_err  = teffs["TEff_WEfficiencyRecoVsAcc"+hshl+"_Top"  ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy        = teffs["TEff_WEfficiencyRecoVsAcc"+hshl+"_Top"  ]->GetEfficiencyErrorUp(nx);	//contains SingleTop
		if(dummy>Top_rec_err ) Top_rec_err  = dummy;
		WT_prob      = teffs["TEff_WEfficiencyRecoVsAll"+hshl+"_WandTop"]->GetEfficiency(nx);
		WT_prob_err  = teffs["TEff_WEfficiencyRecoVsAll"+hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_WEfficiencyRecoVsAll"+hshl+"_WandTop"]->GetEfficiencyErrorUp(nx);
		if(dummy>WT_prob_err ) WT_prob_err  = dummy;
		WT_acc       = teffs["TEff_WEfficiencyAccVsAll" +hshl+"_WandTop"]->GetEfficiency(nx);
		WT_acc_err   = teffs["TEff_WEfficiencyAccVsAll" +hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_WEfficiencyAccVsAll" +hshl+"_WandTop"]->GetEfficiencyErrorUp(nx);
		if(dummy>WT_acc_err  ) WT_acc_err   = dummy;
		WT_rec       = teffs["TEff_WEfficiencyRecoVsAcc"+hshl+"_WandTop"]->GetEfficiency(nx);
		WT_rec_err   = teffs["TEff_WEfficiencyRecoVsAcc"+hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_WEfficiencyRecoVsAcc"+hshl+"_WandTop"]->GetEfficiencyErrorUp(nx);
		if(dummy>WT_rec_err  ) WT_rec_err   = dummy;
		W_prob_nomt       = teffs["TEff_WEfficiencyRecoVsAllNoMT"+hshl+"_WJets"]->GetEfficiency(nx);
		W_prob_nomt_err   = teffs["TEff_WEfficiencyRecoVsAllNoMT"+hshl+"_WJets"]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEfficiencyRecoVsAllNoMT"+hshl+"_WJets"]->GetEfficiencyErrorUp(nx);
		if(dummy>W_prob_nomt_err  ) W_prob_nomt_err   = dummy;
		Top_prob_nomt     = teffs["TEff_WEfficiencyRecoVsAllNoMT"+hshl+"_Top"  ]->GetEfficiency(nx);		//contains SingleTop
		Top_prob_nomt_err = teffs["TEff_WEfficiencyRecoVsAllNoMT"+hshl+"_Top"  ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy             = teffs["TEff_WEfficiencyRecoVsAllNoMT"+hshl+"_Top"  ]->GetEfficiencyErrorUp(nx);	//contains SingleTop
		if(dummy>Top_prob_nomt_err) Top_prob_nomt_err = dummy;
		W_rec_nomt        = teffs["TEff_WEfficiencyRecoVsAccNoMT"+hshl+"_WJets"]->GetEfficiency(nx);
		W_rec_nomt_err    = teffs["TEff_WEfficiencyRecoVsAccNoMT"+hshl+"_WJets"]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEfficiencyRecoVsAccNoMT"+hshl+"_WJets"]->GetEfficiencyErrorUp(nx);
		if(dummy>W_rec_nomt_err   ) W_rec_nomt_err    = dummy;
		Top_rec_nomt      = teffs["TEff_WEfficiencyRecoVsAccNoMT"+hshl+"_Top"  ]->GetEfficiency(nx);		//contains SingleTop
		Top_rec_nomt_err  = teffs["TEff_WEfficiencyRecoVsAccNoMT"+hshl+"_Top"  ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy             = teffs["TEff_WEfficiencyRecoVsAccNoMT"+hshl+"_Top"  ]->GetEfficiencyErrorUp(nx);	//contains SingleTop
		if(dummy>Top_rec_nomt_err ) Top_rec_nomt_err  = dummy;
		WT_prob_nomt      = teffs["TEff_WEfficiencyRecoVsAllNoMT"+hshl+"_WandTop"]->GetEfficiency(nx);
		WT_prob_nomt_err  = teffs["TEff_WEfficiencyRecoVsAllNoMT"+hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEfficiencyRecoVsAllNoMT"+hshl+"_WandTop"]->GetEfficiencyErrorUp(nx);
		if(dummy>WT_prob_nomt_err ) WT_prob_nomt_err  = dummy;
		WT_rec_nomt       = teffs["TEff_WEfficiencyRecoVsAccNoMT"+hshl+"_WandTop"]->GetEfficiency(nx);
		WT_rec_nomt_err   = teffs["TEff_WEfficiencyRecoVsAccNoMT"+hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEfficiencyRecoVsAccNoMT"+hshl+"_WandTop"]->GetEfficiencyErrorUp(nx);
		if(dummy>WT_rec_nomt_err  ) WT_rec_nomt_err   = dummy;
		nW  = hists["LeptonEvents" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT = hists["LeptonEvents" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST = hists["LeptonEvents" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT  = hists["LeptonEvents" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT = hists["LeptonEvents" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_bg  = hists["LeptonEventsBG" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_bg = hists["LeptonEventsBG" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_bg = hists["LeptonEventsBG" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_bg  = hists["LeptonEventsBG" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_bg = hists["LeptonEventsBG" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_dLL  = hists["LeptonEvents_doubleLL" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_dLL = hists["LeptonEvents_doubleLL" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_dLL = hists["LeptonEvents_doubleLL" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_dLL  = hists["LeptonEvents_doubleLL" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_dLL = hists["LeptonEvents_doubleLL" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_bg_SFup    = hists["LeptonEventsBGUp"   + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_bg_SFup   = hists["LeptonEventsBGUp"   + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_bg_SFup   = hists["LeptonEventsBGUp"   + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_bg_SFup    = hists["LeptonEventsBGUp"   + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_bg_SFup   = hists["LeptonEventsBGUp"   + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_bg_SFdown  = hists["LeptonEventsBGDown" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_bg_SFdown = hists["LeptonEventsBGDown" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_bg_SFdown = hists["LeptonEventsBGDown" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_bg_SFdown  = hists["LeptonEventsBGDown" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_bg_SFdown = hists["LeptonEventsBGDown" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_dLL_SFdown  = hists["LeptonEvents_doubleLL_SFdown" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_dLL_SFdown = hists["LeptonEvents_doubleLL_SFdown" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_dLL_SFdown = hists["LeptonEvents_doubleLL_SFdown" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_dLL_SFdown  = hists["LeptonEvents_doubleLL_SFdown" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_dLL_SFdown = hists["LeptonEvents_doubleLL_SFdown" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_dLL_SFup  = hists["LeptonEvents_doubleLL_SFup" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_dLL_SFup = hists["LeptonEvents_doubleLL_SFup" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_dLL_SFup = hists["LeptonEvents_doubleLL_SFup" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_dLL_SFup  = hists["LeptonEvents_doubleLL_SFup" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_dLL_SFup = hists["LeptonEvents_doubleLL_SFup" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_leptveto  = hists["NoLeptonEvents" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_leptveto = hists["NoLeptonEvents" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_leptveto = hists["NoLeptonEvents" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_leptveto  = hists["NoLeptonEvents" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_leptveto = hists["NoLeptonEvents" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_leptveto_SFup  = hists["NoLeptonEvents_SFup" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_leptveto_SFup = hists["NoLeptonEvents_SFup" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_leptveto_SFup = hists["NoLeptonEvents_SFup" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_leptveto_SFup  = hists["NoLeptonEvents_SFup" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_leptveto_SFup = hists["NoLeptonEvents_SFup" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_leptveto_SFdown  = hists["NoLeptonEvents_SFdown" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_leptveto_SFdown = hists["NoLeptonEvents_SFdown" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_leptveto_SFdown = hists["NoLeptonEvents_SFdown" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_leptveto_SFdown  = hists["NoLeptonEvents_SFdown" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_leptveto_SFdown = hists["NoLeptonEvents_SFdown" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_leptveto_err  = hists["NoLeptonEvents" + hshl+"_WJets"    ]->GetBinError(nx);
		nTT_leptveto_err = hists["NoLeptonEvents" + hshl+"_TTbar"    ]->GetBinError(nx);
		nST_leptveto_err = hists["NoLeptonEvents" + hshl+"_SingleTop"]->GetBinError(nx);
		nT_leptveto_err  = hists["NoLeptonEvents" + hshl+"_Top"      ]->GetBinError(nx);
		nWT_leptveto_err = hists["NoLeptonEvents" + hshl+"_WandTop"  ]->GetBinError(nx);
		nW_dLL_lv  = hists["NoLeptonEvents_doubleLL" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_dLL_lv = hists["NoLeptonEvents_doubleLL" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_dLL_lv = hists["NoLeptonEvents_doubleLL" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_dLL_lv  = hists["NoLeptonEvents_doubleLL" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_dLL_lv = hists["NoLeptonEvents_doubleLL" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_dLL_lv_SFup  = hists["NoLeptonEvents_doubleLL_SFup" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_dLL_lv_SFup = hists["NoLeptonEvents_doubleLL_SFup" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_dLL_lv_SFup = hists["NoLeptonEvents_doubleLL_SFup" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_dLL_lv_SFup  = hists["NoLeptonEvents_doubleLL_SFup" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_dLL_lv_SFup = hists["NoLeptonEvents_doubleLL_SFup" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_dLL_lv_SFdown  = hists["NoLeptonEvents_doubleLL_SFdown" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_dLL_lv_SFdown = hists["NoLeptonEvents_doubleLL_SFdown" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_dLL_lv_SFdown = hists["NoLeptonEvents_doubleLL_SFdown" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_dLL_lv_SFdown  = hists["NoLeptonEvents_doubleLL_SFdown" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_dLL_lv_SFdown = hists["NoLeptonEvents_doubleLL_SFdown" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_dLL_lv_err  = hists["NoLeptonEvents_doubleLL" + hshl+"_WJets"    ]->GetBinError(nx);
		nTT_dLL_lv_err = hists["NoLeptonEvents_doubleLL" + hshl+"_TTbar"    ]->GetBinError(nx);
		nST_dLL_lv_err = hists["NoLeptonEvents_doubleLL" + hshl+"_SingleTop"]->GetBinError(nx);
		nT_dLL_lv_err  = hists["NoLeptonEvents_doubleLL" + hshl+"_Top"      ]->GetBinError(nx);
		nWT_dLL_lv_err = hists["NoLeptonEvents_doubleLL" + hshl+"_WandTop"  ]->GetBinError(nx);

		nW_goodevt_dLL  = hists["WEvents_dLL" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_goodevt_dLL = hists["WEvents_dLL" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_goodevt_dLL = hists["WEvents_dLL" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_goodevt_dLL  = hists["WEvents_dLL" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_goodevt_dLL = hists["WEvents_dLL" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nW_goodrecoevt_dLL  = hists["WEventsReco_dLL" + hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_goodrecoevt_dLL = hists["WEventsReco_dLL" + hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_goodrecoevt_dLL = hists["WEventsReco_dLL" + hshl+"_SingleTop"]->GetBinContent(nx);
		nT_goodrecoevt_dLL  = hists["WEventsReco_dLL" + hshl+"_Top"      ]->GetBinContent(nx);
		nWT_goodrecoevt_dLL = hists["WEventsReco_dLL" + hshl+"_WandTop"  ]->GetBinContent(nx);

		QCD_bg   = hists["RecoLepEvents" + hshl+"_QCD"      ]->GetBinContent(nx);
		Z_bg     = hists["RecoLepEvents" + hshl+"_ZJets"    ]->GetBinContent(nx);
		Other_bg = hists["RecoLepEvents" + hshl+"_Other"    ]->GetBinContent(nx);
		TT_bg    = hists["RecoLepEvents" + hshl+"_TTbar"    ]->GetBinContent(nx);
		ST_bg    = hists["RecoLepEvents" + hshl+"_SingleTop"]->GetBinContent(nx);
		T_bg     = hists["RecoLepEvents" + hshl+"_Top"      ]->GetBinContent(nx);
		W_bg     = hists["RecoLepEvents" + hshl+"_WJets"    ]->GetBinContent(nx);
		WT_bg    = hists["RecoLepEvents" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nonWT_bg = hists["RecoLepEvents" + hshl+"_noWandTop"]->GetBinContent(nx);
		QCD_bg_SFup   = hists["RecoLepEventsUp" + hshl+"_QCD"      ]->GetBinContent(nx);
		Z_bg_SFup     = hists["RecoLepEventsUp" + hshl+"_ZJets"    ]->GetBinContent(nx);
		Other_bg_SFup = hists["RecoLepEventsUp" + hshl+"_Other"    ]->GetBinContent(nx);
		TT_bg_SFup    = hists["RecoLepEventsUp" + hshl+"_TTbar"    ]->GetBinContent(nx);
		ST_bg_SFup    = hists["RecoLepEventsUp" + hshl+"_SingleTop"]->GetBinContent(nx);
		T_bg_SFup     = hists["RecoLepEventsUp" + hshl+"_Top"      ]->GetBinContent(nx);
		W_bg_SFup     = hists["RecoLepEventsUp" + hshl+"_WJets"    ]->GetBinContent(nx);
		WT_bg_SFup    = hists["RecoLepEventsUp" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nonWT_bg_SFup = hists["RecoLepEventsUp" + hshl+"_noWandTop"]->GetBinContent(nx);
		QCD_bg_SFdown   = hists["RecoLepEventsDown" + hshl+"_QCD"      ]->GetBinContent(nx);
		Z_bg_SFdown     = hists["RecoLepEventsDown" + hshl+"_ZJets"    ]->GetBinContent(nx);
		Other_bg_SFdown = hists["RecoLepEventsDown" + hshl+"_Other"    ]->GetBinContent(nx);
		TT_bg_SFdown    = hists["RecoLepEventsDown" + hshl+"_TTbar"    ]->GetBinContent(nx);
		ST_bg_SFdown    = hists["RecoLepEventsDown" + hshl+"_SingleTop"]->GetBinContent(nx);
		T_bg_SFdown     = hists["RecoLepEventsDown" + hshl+"_Top"      ]->GetBinContent(nx);
		W_bg_SFdown     = hists["RecoLepEventsDown" + hshl+"_WJets"    ]->GetBinContent(nx);
		WT_bg_SFdown    = hists["RecoLepEventsDown" + hshl+"_WandTop"  ]->GetBinContent(nx);
		nonWT_bg_SFdown = hists["RecoLepEventsDown" + hshl+"_noWandTop"]->GetBinContent(nx);

		nData = hists["RecoLepEvents"+hshl+"_data"]->GetBinContent(nx);
		MC_bg = hists["RecoLepEvents"+hshl+"_mc"]->GetBinContent(nx);
		MC_bg_err = hists["RecoLepEvents"+hshl+"_mc"]->GetBinError(nx);


		//the code below is very ugly
		//on the other hand it is very flexible as all flags set at the beginning of this macros are considered.
		Top_acc = Top_acc*mteff;
		W_acc   = W_acc  *mteff;
		WT_acc   = WT_acc *mteff;
		Top_acc_err = sqrt(pow(Top_acc_err*mteff,2) + pow(Top_acc*mtefferr, 2));
		W_acc_err   = sqrt(pow(  W_acc_err*mteff,2) + pow(  W_acc*mtefferr, 2));
		WT_acc_err  = sqrt(pow( WT_acc_err*mteff,2) + pow( WT_acc*mtefferr, 2));

		//get variables
		double nW_scaled(0.);//number of true one lepton events (i.e. with a gen W-->lnu)
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_scaled = nWT;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_scaled = nT;
		else{	if(fIncludeTop) nW_scaled += nTT;	if(fIncludeSingleTop) nW_scaled += nST;	if(!fTopOnly) nW_scaled += nW;	}
		double nW_bg_scaled(0.);//background to the number above
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_scaled = nWT_bg;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_scaled = nT_bg;
		else{	if(fIncludeTop) nW_bg_scaled += nTT_bg;	if(fIncludeSingleTop) nW_bg_scaled += nST_bg;	if(!fTopOnly) nW_bg_scaled += nW_bg;	}
		double nW_leptveto_scaled(0.);//number of lost leptons
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_leptveto_scaled = nWT_leptveto;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_leptveto_scaled = nT_leptveto;
		else{	if(fIncludeTop) nW_leptveto_scaled += nTT_leptveto;	if(fIncludeSingleTop) nW_leptveto_scaled += nST_leptveto;	if(!fTopOnly) nW_leptveto_scaled += nW_leptveto;	}
		double nW_leptveto_scaled_err(0.);//statistical error
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_leptveto_scaled_err = pow(nWT_leptveto_err,2);
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_leptveto_scaled_err = pow(nT_leptveto_err,2);
		else{	if(fIncludeTop) nW_leptveto_scaled_err += pow(nTT_leptveto_err,2);	if(fIncludeSingleTop) nW_leptveto_scaled_err += pow(nST_leptveto_err,2);	if(!fTopOnly) nW_leptveto_scaled_err += pow(nW_leptveto_err,2);	}
		nW_leptveto_scaled_err = sqrt(nW_leptveto_scaled_err);
		double nW_leptveto_scaled_SFup(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_leptveto_scaled_SFup = pow(nWT_leptveto_SFup,2);
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_leptveto_scaled_SFup = pow(nT_leptveto_SFup,2);
		else{	if(fIncludeTop) nW_leptveto_scaled_SFup += pow(nTT_leptveto_SFup,2);	if(fIncludeSingleTop) nW_leptveto_scaled_SFup += pow(nST_leptveto_SFup,2);	if(!fTopOnly) nW_leptveto_scaled_SFup += pow(nW_leptveto_SFup,2);	}
		nW_leptveto_scaled_SFup = sqrt(nW_leptveto_scaled_SFup);
		double nW_leptveto_scaled_SFdown(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_leptveto_scaled_SFdown = pow(nWT_leptveto_SFdown,2);
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_leptveto_scaled_SFdown = pow(nT_leptveto_SFdown,2);
		else{	if(fIncludeTop) nW_leptveto_scaled_SFdown += pow(nTT_leptveto_SFdown,2);	if(fIncludeSingleTop) nW_leptveto_scaled_SFdown += pow(nST_leptveto_SFdown,2);	if(!fTopOnly) nW_leptveto_scaled_SFdown += pow(nW_leptveto_SFdown,2);	}
		nW_leptveto_scaled_SFdown = sqrt(nW_leptveto_scaled_SFdown);
		double nW_leptveto_scaled_SFerr = fabs(nW_leptveto_scaled-nW_leptveto_scaled_SFup) > fabs(nW_leptveto_scaled-nW_leptveto_scaled_SFdown) ? fabs(nW_leptveto_scaled-nW_leptveto_scaled_SFup) : fabs(nW_leptveto_scaled-nW_leptveto_scaled_SFdown);
		double nW_bg_scaled_SFup(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_scaled_SFup = nWT_bg_SFup;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_scaled_SFup = nT_bg_SFup;
		else{	if(fIncludeTop) nW_bg_scaled_SFup += nTT_bg_SFup;	if(fIncludeSingleTop) nW_bg_scaled_SFup += nST_bg_SFup;	if(!fTopOnly) nW_bg_scaled_SFup += nW_bg_SFup;	}
		double nW_bg_scaled_SFdown(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_scaled_SFdown = nWT_bg_SFdown;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_scaled_SFdown = nT_bg_SFdown;
		else{	if(fIncludeTop) nW_bg_scaled_SFdown += nTT_bg_SFdown;	if(fIncludeSingleTop) nW_bg_scaled_SFdown += nST_bg_SFdown;	if(!fTopOnly) nW_bg_scaled_SFdown += nW_bg_SFdown;	}

		//numbers in case of double lost lepton
		double nW_dLL_scaled(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_dLL_scaled = nWT_dLL;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_dLL_scaled = nT_dLL;
		else{	if(fIncludeTop) nW_dLL_scaled += nTT_dLL;	if(fIncludeSingleTop) nW_dLL_scaled += nST_dLL;	if(!fTopOnly) nW_dLL_scaled += nW_dLL;	}
		double nW_dLL_lv_scaled(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_dLL_lv_scaled = nWT_dLL_lv;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_dLL_lv_scaled = nT_dLL_lv;
		else{	if(fIncludeTop) nW_dLL_lv_scaled += nTT_dLL_lv;	if(fIncludeSingleTop) nW_dLL_lv_scaled += nST_dLL_lv;	if(!fTopOnly) nW_dLL_lv_scaled += nW_dLL_lv;	}
		double nW_dLL_lv_err_scaled(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_dLL_lv_err_scaled = nWT_dLL_lv_err;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_dLL_lv_err_scaled = nT_dLL_lv_err;
		else{	if(fIncludeTop) nW_dLL_lv_err_scaled += nTT_dLL_lv_err;	if(fIncludeSingleTop) nW_dLL_lv_err_scaled += nST_dLL_lv_err;	if(!fTopOnly) nW_dLL_lv_err_scaled += nW_dLL_lv_err;	}
		double nW_dLL_lv_scaled_SFup(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_dLL_lv_scaled_SFup = pow(nWT_dLL_lv_SFup,2);
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_dLL_lv_scaled_SFup = pow(nT_dLL_lv_SFup,2);
		else{	if(fIncludeTop) nW_dLL_lv_scaled_SFup += pow(nTT_dLL_lv_SFup,2);	if(fIncludeSingleTop) nW_dLL_lv_scaled_SFup += pow(nST_dLL_lv_SFup,2);	if(!fTopOnly) nW_dLL_lv_scaled_SFup += pow(nW_dLL_lv_SFup,2);	}
		nW_dLL_lv_scaled_SFup = sqrt(nW_dLL_lv_scaled_SFup);
		double nW_dLL_lv_scaled_SFdown(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_dLL_lv_scaled_SFdown = pow(nWT_dLL_lv_SFdown,2);
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_dLL_lv_scaled_SFdown = pow(nT_dLL_lv_SFdown,2);
		else{	if(fIncludeTop) nW_dLL_lv_scaled_SFdown += pow(nTT_dLL_lv_SFdown,2);	if(fIncludeSingleTop) nW_dLL_lv_scaled_SFdown += pow(nST_dLL_lv_SFdown,2);	if(!fTopOnly) nW_dLL_lv_scaled_SFdown += pow(nW_dLL_lv_SFdown,2);	}
		nW_dLL_lv_scaled_SFdown = sqrt(nW_dLL_lv_scaled_SFdown);
		double nW_dLL_scaled_SFup(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_dLL_scaled_SFup = nWT_dLL_SFup;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_dLL_scaled_SFup = nT_dLL_SFup;
		else{	if(fIncludeTop) nW_dLL_scaled_SFup += nTT_dLL_SFup;	if(fIncludeSingleTop) nW_dLL_scaled_SFup += nST_dLL_SFup;	if(!fTopOnly) nW_dLL_scaled_SFup += nW_dLL_SFup;	}
		double nW_dLL_scaled_SFdown(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_dLL_scaled_SFdown = nWT_dLL_SFdown;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_dLL_scaled_SFdown = nT_dLL_SFdown;
		else{	if(fIncludeTop) nW_dLL_scaled_SFdown += nTT_dLL_SFdown;	if(fIncludeSingleTop) nW_dLL_scaled_SFdown += nST_dLL_SFdown;	if(!fTopOnly) nW_dLL_scaled_SFdown += nW_dLL_SFdown;	}

		//depending what you do with double lost lepton, they are background events (to this method), or events that you want to predict
		if(fDoubleLL_BG){
			nW_bg_scaled += nW_dLL_scaled;
			nW_bg_scaled_SFup += nW_dLL_scaled_SFup;
			nW_bg_scaled_SFdown += nW_dLL_scaled_SFdown;
		} else {
			nW_leptveto_scaled += nW_dLL_lv_scaled;
			nW_leptveto_scaled_err = sqrt(pow(nW_leptveto_scaled_err,2)+pow(nW_dLL_lv_err_scaled,2));
			nW_leptveto_scaled_SFup += nW_dLL_scaled_SFup;
			nW_leptveto_scaled_SFdown += nW_dLL_scaled_SFdown;
			nW_leptveto_scaled_SFerr = fabs(nW_leptveto_scaled-nW_leptveto_scaled_SFup) > fabs(nW_leptveto_scaled-nW_leptveto_scaled_SFdown) ? fabs(nW_leptveto_scaled-nW_leptveto_scaled_SFup) : fabs(nW_leptveto_scaled-nW_leptveto_scaled_SFdown);
		}

		//total background is background from W+Top sample, as well other backgrounds
		//note here, that we do not consider double lost lepton from Z (they have no genuine MET)
		double bg = nW_bg_scaled + Z_bg + QCD_bg + Other_bg;
		if( fTopOnly)          bg +=  W_bg;
		if(!fIncludeTop)       bg += TT_bg;
		if(!fIncludeSingleTop) bg += ST_bg;
		double bg_SFup = nW_bg_scaled_SFup + Z_bg_SFup + QCD_bg_SFup + Other_bg_SFup;
		if( fTopOnly)          bg_SFup +=  W_bg_SFup;
		if(!fIncludeTop)       bg_SFup += TT_bg_SFup;
		if(!fIncludeSingleTop) bg_SFup += ST_bg_SFup;
		double bg_SFdown = nW_bg_scaled_SFdown + Z_bg_SFdown + QCD_bg_SFdown + Other_bg_SFdown;
		if( fTopOnly)          bg_SFdown +=  W_bg_SFdown;
		if(!fIncludeTop)       bg_SFdown += TT_bg_SFdown;
		if(!fIncludeSingleTop) bg_SFdown += ST_bg_SFdown;
        	double allMC = WT_bg + nonWT_bg;
        	double allMC_SFup   = WT_bg_SFup   + nonWT_bg_SFup;  
        	double allMC_SFdown = WT_bg_SFdown + nonWT_bg_SFdown;
		double allMC_SFerr = fabs(allMC-allMC_SFup) > fabs(allMC-allMC_SFdown) ? fabs(allMC-allMC_SFup) : fabs(allMC-allMC_SFdown);
		double bg_SFerr = fabs(bg-bg_SFup) > fabs(bg-bg_SFdown) ? fabs(bg-bg_SFup) : fabs(bg-bg_SFdown);

		//This is at your own risk!!!!
		//the background is rescaled by the factor (data yield)/(expected yield from simulation)
		//this was implemented, as in some cases the background was larger than the signal
		//because the total MC yield >> data yield (e.g. strong down fluctuations in data).
		//Thus this rescaling 'improves' estimate in case data<<MC, however 'punishes' if data>>MC
		if(allMC>0){
			bg = bg * nData / allMC;
			bg_SFup = bg_SFup * nData / allMC;
			bg_SFdown = bg_SFdown * nData / allMC;
			bg_SFerr = bg_SFerr * nData / allMC;
		}

		//as this is for 'closure tests', additional uncertainties on efficiencies, etc. are set to 0
		//these are the lepton efficiencies (reconstruction and acceptance)
		double prob_MC(-1.);
		if(fWeightedProb) prob_MC = WT_prob;
		else{
			if(fTopEfficencies)  prob_MC = Top_prob;
			else                prob_MC = W_prob;
		}
		double prob_MC_err_sys;
		if(fWeightedProb) prob_MC_err_sys = WT_prob_err;
		else{
			if(fTopEfficencies) prob_MC_err_sys = Top_prob_err;
			else prob_MC_err_sys = W_prob_err;
		}
		//noMTtest
		double prob_MC_nomt(-1.);
		if(fWeightedProb) prob_MC_nomt = WT_prob_nomt;
		else{
			if(fTopEfficencies)  prob_MC_nomt = Top_prob_nomt;
			else                prob_MC_nomt = W_prob_nomt;
		}
		double prob_MC_nomt_err_sys;
		if(fWeightedProb) prob_MC_nomt_err_sys = WT_prob_nomt_err;
		else{
			if(fTopEfficencies) prob_MC_nomt_err_sys = Top_prob_nomt_err;
			else prob_MC_nomt_err_sys = W_prob_nomt_err;
		}

		//in order to not double count MC statistics, the statistical error is set to 0, also systematical error due to data/MC differences are set to 0 (this is 'for closure tests')
		double pred_MC_nomt1 = (allMC-bg)*(1-prob_MC_nomt)/(prob_MC_nomt*mteff);
		double pred_MC_nomt1_error_stat = 0;
		double pred_MC_nomt1_error_sys = sqrt(pow(pred_MC_nomt1_error_sys,2) + pow(MC_bg_err*(1-prob_MC_nomt)/(prob_MC_nomt*mteff),2));
		double pred_MC_nomt2 = (allMC-bg)*(1-prob_MC_nomt)/(prob_MC_nomt*mteff2);
		double pred_MC_nomt2_error_stat = 0;
		double pred_MC_nomt2_error_sys = sqrt(pow(pred_MC_nomt2_error_sys,2) + pow((MC_bg_err/mteff)*(1-prob_MC_nomt)/prob_MC_nomt,2) + pow((mtefferr/(mteff*mteff))*(allMC-bg)*(1-prob_MC_nomt)/prob_MC_nomt,2));
		double pred_MC = (allMC-bg)*(1-prob_MC_nomt)/prob_MC;
		double pred_MC_error_stat = 0;
		double pred_MC_error_sys = sqrt(pow(pred_MC_error_sys,2) + pow(MC_bg_err*(1-prob_MC_nomt)/prob_MC,2));
		if(fbTagError) pred_MC_error_sys = sqrt(pow(pred_MC_error_sys,2) + pow(allMC_SFerr*(1-prob_MC_nomt)/prob_MC,2) + pow(bg_SFerr*(1-prob_MC_nomt)/prob_MC,2));

		//these are the lepton efficiencies (reconstruction and acceptance) for data (i.e. including extra uncertainties set at the beginning of this macro)
		double prob(-1.);
		if(fWeightedProb) prob = WT_prob;
		else{
			if(fTopEfficencies)  prob = Top_prob;
			else                prob = W_prob;
		}
		double prob_err_sys;
		//factor 2 on rel_sys_uncert, one for lepton efficiency, one for MT cut efficiency
		if(fWeightedProb) prob_err_sys = sqrt(pow(2.*rel_sys_uncert*WT_prob,2)+ pow(WT_prob_err,2));
		else{
			if(fTopEfficencies){
				prob_err_sys = sqrt(pow(2.*rel_sys_uncert*Top_prob,2)+ pow(Top_prob_err,2));
			}else{
				prob_err_sys = sqrt(pow(2.*rel_sys_uncert*W_prob,2)+ pow(W_prob_err,2));
			}
		}
		//no MT means that there was no MT cut - used for cross checks
		double prob_nomt(-1.);
		if(fWeightedProb) prob_nomt = WT_prob_nomt;
		else{
			if(fTopEfficencies) prob_nomt = Top_prob_nomt;
			else                prob_nomt = W_prob_nomt;
		}
		double prob_nomt_err_sys;
		//factor 2 on rel_sys_uncert, one for lepton efficiency, one for MT cut efficiency
		if(fWeightedProb) prob_nomt_err_sys = sqrt(pow(rel_sys_uncert*WT_prob_nomt,2)+ pow(WT_prob_nomt_err,2));
		else{
			if(fTopEfficencies){
				prob_nomt_err_sys = sqrt(pow(rel_sys_uncert*Top_prob_nomt,2)+ pow(Top_prob_nomt_err,2));
			}else{
				prob_nomt_err_sys = sqrt(pow(rel_sys_uncert*W_prob_nomt,2)+ pow(W_prob_nomt_err,2));
			}
		}
		//these are the predicted numbers
		double pred = (nData - bg)*(1.-prob_nomt)/(prob_nomt*mteff);
		double pred_error_stat = fabs(sqrt(nData)*(1.-prob_nomt)/(prob_nomt*mteff));
		double pred_error_sys  = sqrt(  pow((nData-bg) *(prob_nomt_err_sys/(prob_nomt*prob_nomt*mteff)),2) + pow((nData-bg) *(mtefferr*(1.-prob_nomt)/(prob_nomt*mteff*mteff)),2) + pow(rel_sys_uncert_bg*bg*(1.-prob_nomt)/(prob_nomt*mteff),2));
		if(mtefferr==1) pred_error_sys  = sqrt(  pow((nData-bg) *(prob_nomt_err_sys/(prob_nomt*prob_nomt*mteff)),2) + pow(rel_sys_uncert_bg*bg*(1.-prob_nomt)/(prob_nomt*mteff),2));
		if(fbTagError) pred_error_sys = sqrt( pow(pred_error_sys,2) + pow(bg_SFerr *(1-prob_nomt)/(prob_nomt*mteff),2) );
		pred_error_sys = sqrt(pred_error_sys*pred_error_sys +    pow((frelMTerr*(1.-prob_nomt)/(prob_nomt*mteff*mteff)),2) );//addtional MT cut error
		pred_error_sys = sqrt(pred_error_sys*pred_error_sys + pow(nWT_goodrecoevt_dLL*(fdoubeLLerr*prob_nomt_err_sys/(prob_nomt*prob_nomt*mteff)),2) );//additional error due to double lost lepton

	//here we store all variables that go into the tables, naming should be clear
	mctruth.push_back(nW_leptveto_scaled);
	if(!fbTagError) mctrutherr.push_back(nW_leptveto_scaled_err);
	else            mctrutherr.push_back(sqrt(pow(nW_leptveto_scaled_err,2)+pow(nW_leptveto_scaled_SFerr,2) ) );
	mctruthfractionW.push_back(nW_leptveto);
	mctruthfractionTT.push_back(nTT_leptveto);
	datapred.push_back(pred);
	datapred_stat_err.push_back(pred_error_stat);
	datapred_syst_err.push_back(pred_error_sys );
	MCpred.push_back(pred_MC);
	MCpred_stat_err.push_back(pred_MC_error_stat );
	MCpred_syst_err.push_back(pred_MC_error_sys );
	LLeff.push_back(prob_nomt);
	LLefferr.push_back(prob_nomt_err_sys);
	MTeff.push_back(mteff);
	MTefferr.push_back(mtefferr);
	numtrueW.push_back(nW_scaled);
	numtrueW_bg.push_back(nW_bg_scaled);
	numData.push_back(nData);
	numBG.push_back(bg);
	numQCD.push_back(QCD_bg);
	numZ.push_back(Z_bg);
	numW.push_back(W_bg);
	numT.push_back(T_bg);
	numTT.push_back(TT_bg);
	numOther.push_back(Other_bg);
	numMC.push_back(allMC);
	BGerr.push_back(MC_bg_err);

	//the full printouts are useful for debugging
	//they are very chaotic, but the user should self define if he/she wants to use them
	//and if yes, modify them to there likings
	if(makeFullPrintout){
		*fLogStream << "*******************************************" << endl;
		*fLogStream << "************* FULL PRINTOUT ***************" << endl;
		*fLogStream << "*******************************************" << endl;
		*fLogStream << "Signal region " << signal_region[i2] << ", lepton type " << lepton_type[i4] << ", HT bin " << HT_bin[i3] << " and MT2bin " << nx << " = (" << hists[string("RecoLepEvents"+hshl+"_mc")]->GetBinLowEdge(nx) << "," << hists[string("RecoLepEvents"+hshl+"_mc")]->GetBinLowEdge(nx) + hists[string("RecoLepEvents"+hshl+"_mc")]->GetBinWidth(nx) << "):" << endl;
		*fLogStream << "One Lepton yield: QCD = " << QCD_bg << ", Z = " << Z_bg << ", Other = " << Other_bg << " ==> non WandTop = " << nonWT_bg << endl;
		*fLogStream << "                  TTbar = " << TT_bg << ", SingleTop = " << ST_bg << ", W = " << W_bg << " --> Top = " << T_bg << " ==> WandTop = " << WT_bg << endl;
		*fLogStream << "                  Data = " << nData <<  ", total mc = " << MC_bg << "+/-" << MC_bg_err << endl;
		*fLogStream << "1l yield SF down: QCD = " << QCD_bg_SFdown << ", Z = " << Z_bg_SFdown << ", Other = " << Other_bg_SFdown << " ==> non WandTop = " << nonWT_bg_SFdown << endl;
		*fLogStream << "                  TTbar = " << TT_bg_SFdown << ", SingleTop = " << ST_bg_SFdown << ", W = " << W_bg_SFdown << " --> Top = " << T_bg_SFdown << " ==> WandTop = " << WT_bg_SFdown << endl;
		*fLogStream << "1l yield SF up:   QCD = " << QCD_bg_SFup << ", Z = " << Z_bg_SFup << ", Other = " << Other_bg_SFup << " ==> non WandTop = " << nonWT_bg_SFup << endl;
		*fLogStream << "                  TTbar = " << TT_bg_SFup << ", SingleTop = " << ST_bg_SFup << ", W = " << W_bg_SFup << " --> Top = " << T_bg_SFup << " ==> WandTop = " << WT_bg_SFup << endl;
		*fLogStream << "True one lepton events:" << endl;
		*fLogStream << "                  TTbar = " << nTT << ", SingleTop = " << nST << ", W = " << nW << " --> Top = " << nT << " ==> WandTop = " << nWT << endl;
		*fLogStream << "one lepton events without 1 gen lepton from W:" << endl;
		*fLogStream << "                  TTbar = " << nTT_bg << ", SingleTop = " << nST_bg << ", W = " << nW_bg << " --> Top = " << nT_bg << " ==> WandTop = " << nWT_bg << endl;
		*fLogStream << "SFdown            TTbar = " << nTT_bg_SFdown << ", SingleTop = " << nST_bg_SFdown << ", W = " << nW_bg_SFdown << " --> Top = " << nT_bg_SFdown << " ==> WandTop = " << nWT_bg_SFdown << endl;
		*fLogStream << "SFup              TTbar = " << nTT_bg_SFup << ", SingleTop = " << nST_bg_SFup << ", W = " << nW_bg_SFup << " --> Top = " << nT_bg_SFup << " ==> WandTop = " << nWT_bg_SFup << endl;
		*fLogStream << "Zero lepton events, but with one gen lepton from W:" << endl;
		*fLogStream << "                  TTbar = " << nTT_leptveto << "+/-" << nTT_leptveto_err << ", SingleTop = " << nST_leptveto << "+/-" << nST_leptveto_err << ", W = " << nW_leptveto << "+/-" << nW_leptveto_err << " --> Top = " << nT_leptveto << "+/-" << nT_leptveto_err << " ==> WandTop = " << nWT_leptveto << "+/-" << nWT_leptveto_err << endl;
		*fLogStream << endl;
		*fLogStream << "Events with gen lepton from W:                       W = " << W_evts      << ", TTbar = " << TT_evts      << ", WandTop = " << WT_evts      << endl;
		*fLogStream << "Events with gen lepton from W within acceptance:     W = " << W_evts_acc  << ", TTbar = " << TT_evts_acc  << ", WandTop = " << WT_evts_acc  << endl;
		*fLogStream << "Events with gen lepton from W with also reco lepton: W = " << W_evts_reco << ", TTbar = " << TT_evts_reco << ", WandTop = " << WT_evts_reco << endl;
		*fLogStream << "Events as line before, but no MT cut:                W = " << W_evts_nomt_reco << ", TTbar = " << TT_evts_nomt_reco << ", WandTop = " << WT_evts_nomt_reco << endl;
		*fLogStream << "events with double LL (reconstructed):" << endl;
		*fLogStream << "        TTbar = " << nTT_goodevt_dLL << " (" << nTT_goodrecoevt_dLL << "), SingleTop = " << nST_goodevt_dLL << " (" << nST_goodrecoevt_dLL << "), W = " << nW_goodevt_dLL << " (" << nW_goodrecoevt_dLL << ") --> Top = " << nT_goodevt_dLL << " (" << nT_goodrecoevt_dLL << ") ==> WandTop = " << nWT_goodevt_dLL << " (" << nWT_goodrecoevt_dLL << ")" << endl;
		*fLogStream << "events with double LL (one lep reco):" << endl;
		*fLogStream << "        TTbar = " << nTT_dLL << ", SingleTop = " << nST_dLL << ", W = " << nW_dLL << " --> Top = " << nT_dLL << " ==> WandTop = " << nWT_dLL << endl;
		*fLogStream << "events with double LL (no lep reco):" << endl;
		*fLogStream << "        TTbar = " << nTT_dLL_lv << ", SingleTop = " << nST_dLL_lv << ", W = " << nW_dLL_lv << " --> Top = " << nT_dLL_lv << " ==> WandTop = " << nWT_dLL_lv << endl;

		*fLogStream << "Evts passing MT cut       = " << mtpass << ", all leptonic events (no MT cut) = " << mttot << endl;
		*fLogStream << "Evts passing MT cut 2 WT  = " << mtpass2 << ", all leptonic events (no MT cut) = " << mttot2 << endl;
		*fLogStream << "Evts passing MT cut 2 MC  = " << mtpass3 << ", all leptonic events (no MT cut) = " << mttot3 << endl;
		*fLogStream << "MT efficiency: MC         = " << mteff << "+/-" << mtefferr << ", data = " << mteffdata << "+/-" << mtefferrdata << endl;
		*fLogStream << "MT efficiency2: W+Top     = " << mteff2 << "+/-" << mtefferr2 << ", MC   = " << mteff3 << "+/-" << mtefferr3 << endl;
		*fLogStream << "Acceptance:             W = " << W_acc << "+/-" << W_acc_err << ", Top = " << Top_acc << "+/-" << Top_acc_err << ", WandTop = " << WT_acc << "+/-" << WT_acc_err << endl;
		*fLogStream << "Reco Efficiency:        W = " << W_rec << "+/-" << W_rec_err << ", Top = " << Top_rec << "+/-" << Top_rec_err << ", WandTop = " << WT_rec << "+/-" << WT_rec_err << endl;
		*fLogStream << "Total Efficiency:       W = " << W_prob << "+/-" << W_prob_err << ", Top = " << Top_prob << "+/-" << Top_prob_err << ", WandTop = " << WT_prob << "+/-" << WT_prob_err << endl;
		*fLogStream << "Reco Efficiency, noMT:  W = " << W_rec_nomt << "+/-" << W_rec_nomt_err << ", Top = " << Top_rec_nomt << "+/-" << Top_rec_nomt_err << ", WandTop = " << WT_rec_nomt << "+/-" << WT_rec_nomt_err << endl;
		*fLogStream << "Total Efficiency, noMT: W = " << W_prob_nomt << "+/-" << W_prob_nomt_err << ", Top = " << Top_prob_nomt << "+/-" << Top_prob_nomt_err << ", WandTop = " << WT_prob_nomt << "+/-" << WT_prob_nomt_err << endl;

		*fLogStream << endl;
		*fLogStream << "Number of leptonic events from W decay = " << nW_scaled << endl;
		*fLogStream << "Background to that number (i.e leptonic events that are not true genlepts) = " << nW_bg_scaled << " (SFdown/up = " << nW_bg_scaled_SFdown << "/" << nW_bg_scaled_SFup << ")" << endl;
		*fLogStream << "Total background    = " << bg << " (SFerr = " << bg_SFerr << " due to SFdown/up = " << bg_SFdown << "/" << bg_SFup << ")" << endl;
		*fLogStream << "Total MC            = " << allMC << " (SFerr = " << allMC_SFerr << " due to SFdown/up = " << allMC_SFdown << "/" << allMC_SFup << ")" << endl;
		*fLogStream << "prob (epsilon) MC   = " << prob_MC_nomt << "+/-" << prob_MC_nomt_err_sys << " ==> R_LL = " << (1.-prob_MC_nomt)/(prob_MC_nomt*mteff) << "+/-" << sqrt(pow(prob_MC_nomt_err_sys/(mteff*pow(prob_MC_nomt,2)),2) + pow(mtefferr*(1.-prob_MC_nomt)/(prob_MC_nomt*pow(mteff,2)),2)) << endl;
		*fLogStream << "prob (epsilon) data = " << prob_nomt << "+/-" << prob_nomt_err_sys << " ==> R_LL = " << (1.-prob_nomt)/(prob_nomt*mteff) << "+/-" << sqrt(pow(prob_nomt_err_sys/(mteff*pow(prob_nomt,2)),2) + pow(mtefferr*(1.-prob_nomt)/(prob_nomt*pow(mteff,2)),2)) << endl;
		//*fLogStream << "probNoMT (eps) MC 1 = " << prob_MC_nomt << "+/-" << prob_MC_nomt_err_sys << " ==> R_LL = " << (1.-prob_MC_nomt)/prob_MC_nomt << "+/-" << prob_MC_nomt_err_sys/(pow(prob_MC_nomt,2)) << " ==> R_LL*mteff = " << mteff*(1.-prob_MC_nomt)/prob_MC_nomt << "+/-" << sqrt(pow(mteff*prob_MC_nomt_err_sys/(pow(prob_MC_nomt,2)),2)+pow(mtefferr*(1.-prob_MC_nomt)/prob_MC_nomt,2)) << endl;
		//*fLogStream << "probNoMT*MTeff MC   = " << prob_MC_nomt*mteff << "+/-" << sqrt(pow(mteff*prob_MC_nomt_err_sys,2)+pow(mtefferr*prob_MC_nomt,2)) << " ==> R_LL = " << (1.-mteff*prob_MC_nomt)/(mteff*prob_MC_nomt) << "+/-" << sqrt(pow(mteff*prob_MC_nomt_err_sys,2)+pow(mtefferr*prob_MC_nomt,2))/(pow(mteff*prob_MC_nomt,2)) << endl;
		*fLogStream << "prediction     MC   = " << pred_MC << " +/- " << pred_MC_error_stat << "(stat) +/- " << pred_MC_error_sys << "(syst)" << endl;
		*fLogStream << "prediction     data = " << pred << " +/- " << pred_error_stat << "(stat) +/- " << pred_error_sys << "(syst)" << endl;
		*fLogStream << "Xprediction    MC 1 = " << pred_MC_nomt1 << " +/- " << pred_MC_nomt1_error_stat << "(stat) +/- " << pred_MC_nomt1_error_sys << "(syst)" << endl;
		*fLogStream << "Xprediction/mteff MC2=" << pred_MC_nomt2 << " +/- " << pred_MC_nomt2_error_stat << "(stat) +/- " << pred_MC_nomt2_error_sys << "(syst)" << endl;
		*fLogStream << "MC truth of lost leptons = " << nW_leptveto_scaled << " +/- " << nW_leptveto_scaled_err << endl;
		*fLogStream << "MC truth / MC pred  = " << nW_leptveto_scaled/pred_MC << " +/- " << sqrt(pow(nW_leptveto_scaled_err/pred_MC,2) + pow(nW_leptveto_scaled*sqrt(pred_MC_error_stat*pred_MC_error_stat+pred_MC_error_sys*pred_MC_error_sys)/(pred_MC*pred_MC),2)) << endl;
		*fLogStream << "MC pred  / MC truth = " << pred_MC/nW_leptveto_scaled << " +/- " << sqrt( pow(sqrt(pred_MC_error_stat*pred_MC_error_stat+pred_MC_error_sys*pred_MC_error_sys) / nW_leptveto_scaled,2) + pow(pred_MC*nW_leptveto_scaled_err/(nW_leptveto_scaled*nW_leptveto_scaled),2)) << endl;
		*fLogStream << "XMC 1 truth / MC pred = " << nW_leptveto_scaled/pred_MC_nomt1 << " +/- " << sqrt(pow(nW_leptveto_scaled_err/pred_MC_nomt1,2) + pow(nW_leptveto_scaled*sqrt(pred_MC_nomt1_error_stat*pred_MC_nomt1_error_stat+pred_MC_nomt1_error_sys*pred_MC_nomt1_error_sys)/(pred_MC_nomt1*pred_MC_nomt1),2)) << endl;
		*fLogStream << "XMC 2 truth / MC pred = " << nW_leptveto_scaled/pred_MC_nomt2 << " +/- " << sqrt(pow(nW_leptveto_scaled_err/pred_MC_nomt2,2) + pow(nW_leptveto_scaled*sqrt(pred_MC_nomt2_error_stat*pred_MC_nomt2_error_stat+pred_MC_nomt2_error_sys*pred_MC_nomt2_error_sys)/(pred_MC_nomt2*pred_MC_nomt2),2)) << endl;
		*fLogStream << "R_LL truth = MC truth / (nW-nBg) = " << nW_leptveto_scaled/(nW_scaled-nW_bg_scaled) << endl;
		*fLogStream << "R_LL truth * mteff = " << (nW_leptveto_scaled/(nW_scaled-nW_bg_scaled))*mteff << endl;
		*fLogStream << "==> prob_eff truth               = " << 1./((nW_leptveto_scaled/(nW_scaled-nW_bg_scaled))+1.) << endl;
		*fLogStream << "==> prob_eff truth * mteff       = " << (1./((nW_leptveto_scaled/(nW_scaled-nW_bg_scaled))+1.))*mteff << endl;
		*fLogStream << "==> prob_eff truth mteff-corr    = " << 1./((nW_leptveto_scaled/(nW_scaled-nW_bg_scaled))*mteff+1.) << endl;
		*fLogStream << "*******************************************" << endl << endl << endl;

	}
	
	}//for(int nx = 1; nx<=hists[string("RecoLepEvents"+hshl+"_mc")]->GetNbinsX(); ++nx)
	}}}// i2, i3, i4

	//after you just stored all variables for the tables you want, you produce now all tables you want (as set by the flags of function GetLostLeptonEstimate(...)

	//this makes a compact card of signal region vs MCtruth +/- err vs. pred. +/- err
	if(PrintPredictionCard) PredictionCard(sr, htr, lepr, MT2bin, mctruth, mctrutherr, datapred, datapred_stat_err, datapred_syst_err, LLeff, LLefferr, MTeff, MTefferr);
	//this makes a compact card of signal region vs MCtruth +/- err vs. pred. +/- err, and also prints out (from MC truth) the fraction of expected lost lepton from W,Top,rest(dibosons?)
	if(PrintPredictionCard) PredictionCardSplitted(sr, htr, lepr, MT2bin, mctruth, mctrutherr, datapred, datapred_stat_err, datapred_syst_err, LLeff, LLefferr, MTeff, MTefferr, mctruthfractionW, mctruthfractionTT, numW, numTT, numMC);
	//this just makes the one lepton yield table (usually with MT cut applied for one lepton selection, unless flags are set differently)
	if(PrintYieldTable) YieldTable(sr, htr, lepr, MT2low, MT2up, numQCD, numZ, numW, numT, numOther, numMC, BGerr, numData);
	//this line produces the final result tables
	if(PrintSummaryTable) SummaryTable(version, sr, htr, lepr, MT2low, MT2up, numtrueW, numtrueW_bg, numData, numBG, LLeff, LLefferr, MTeff, MTefferr, mctruth, mctrutherr, MCpred, MCpred_stat_err, MCpred_syst_err, datapred, datapred_stat_err, datapred_syst_err);
	//this function stores results into a root file, also can produce MCtruth / prediction plots
	if(SavePrediction) PredictionFile(version, sr, htr, lepr, MT2low, MT2up, numtrueW, numtrueW_bg, numData, numBG, LLeff, LLefferr, MTeff, MTefferr, mctruth, mctrutherr, MCpred, MCpred_stat_err, MCpred_syst_err, datapred, datapred_stat_err, datapred_syst_err, PlotPrediction);

}

//this function produces the final result tables
//version == 0,1: with MCPred
//version == 2,3: without MCPred
//version == 0,2: with R_LL (i.e. 1-e / e)
//version == 1,3: with e
void SummaryTable(int version, vector<int> sr, vector<int> htr, vector<int> lepr, vector<double> MT2low, vector<double> MT2up, vector<double> numtrueW, vector<double> numtrueW_bg, vector<double> numData, vector<double> numBG, vector<double> LLeff, vector<double> LLefferr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruth, vector<double> mctrutherr, vector<double> MCpred, vector<double> MCpred_stat_err, vector<double> MCpred_syst_err, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err){

	*fLogStream << "*********************************************************************" << endl;
	*fLogStream << "\%BEGINLATEX\%"             << endl;
	*fLogStream << "\\begin{table}"             << endl
	     << "\\begin{center}"            << endl
	     << "\\small"                    << endl;
       	*fLogStream << "\\begin{tabular}{lccccccc}" << endl;	     
	*fLogStream << "\\hline\\hline"             << endl;

	if(!fRebin) *fLogStream << "$M_{T2}$ (GeV) & ";
	else        *fLogStream << "signal region       & ";
	if(fIncludeTop && !fTopOnly) *fLogStream << "$N^{MC}(W \\& Top)"<< "$ & $";
	if(fTopOnly)                 *fLogStream << "N^{MC}(Top)"<< "$ & $";
	if(!fIncludeTop)             *fLogStream << "N^{MC}(W)"<< "$ & $";
	                             *fLogStream << "N^{reco}" << "$ &  $" << "N^{bg}"  <<  "$  &     $";  
	if(version==1||version==3)   *fLogStream <<  "\\varepsilon"  << "$ &    $" << "N^{pass}$ MC      " << " & $";
	if(version==0||version==2)   *fLogStream <<  "R_{LL}"   << "$    &      $" << "N^{pass}$ MC      " << " &     $";
	if(version==0||version==1)   *fLogStream << "N^{pass}$ MCPred    " << " &                    $";
	                             *fLogStream << "N^{pass}$ Pred                      " << " \\\\" << endl;
	string oldlep = "dummy";
	string oldsr = "dummy";
	string oldHT = "dummy";
	for(unsigned int n = 0; n<sr.size(); ++n){
		string lep;
		if(lepr[n]==0) lep="Muo";
		else if(lepr[n]==1) lep="Ele";
		else continue;
		string sigreg;
		if(sr[n]==0) sigreg="2j, 0b";
		if(sr[n]==1) sigreg="2j, $\\ge 1$b";
		if(sr[n]==2) sigreg="$3-5$j, 0b";
		if(sr[n]==3) sigreg="$3-5$j, 1b";
		if(sr[n]==4) sigreg="$3-5$j, 2b";
		if(sr[n]==5) sigreg="$\\ge 6$j, 0b";
		if(sr[n]==6) sigreg="$\\ge 6$j, 1b";
		if(sr[n]==7) sigreg="$\\ge 6$j, 2b";
		if(sr[n]==8) sigreg="$\\ge 3$j, $\\ge 3$b";
		string htreg;
		if(htr[n]==0 && fMET) htreg = "450 GeV $\\leq H_{T} < 750$ GeV";
		if(htr[n]==1 && fHT) htreg = "750 GeV $\\leq H_{T} < 1200$ GeV";
		if(htr[n]==2 && fHT) htreg = "$H_{T}\\ge 1200$ GeV";
		if(htreg!=oldHT){//this if makes sure that HT region is only printed out, when the vector jumps to the next HT region
			*fLogStream << " \\hline " << endl << "\\multicolumn{7}{l}{" <<  htreg << "} \\\\ " << endl << "\\hline" << endl;
			oldHT = htreg;
		}
		if(lep!=oldlep){//this if makes sure that leptoninc region is only printed out, when the vector jumps to the next leptonic region
			*fLogStream << " \\hline " << endl << lep << "\\\\ " << endl << "\\hline" << endl;
			oldlep = lep;
		}
		if(!fRebin){
		if(sigreg!=oldsr){//this if makes sure that topological region is only printed out, when the vector jumps topological the next HT region
			*fLogStream << " \\hline  " << endl << sigreg << "\\\\" /* << endl << "\\hline"*/ << endl;
			oldsr = sigreg;
		}
		if(MT2up[n]==10000.) *fLogStream << "$" << int(MT2low[n]) << "-" << "\\infty$" << " " << setw(4) << " & ";
		else                 *fLogStream << "$" << int(MT2low[n]) << "-" << int(MT2up[n]) << "$" << " " << setw(7) << " & ";
		}
		//after running the code ones successfully you will understand the below
		else *fLogStream << " " << setw(18) << sigreg << " & ";
		*fLogStream << fixed << setprecision(2) << " " << setw(17) << numtrueW[n]-numtrueW_bg[n] << " & ";
		*fLogStream << " " << setw(10) << int(numData[n]) << " & ";
		*fLogStream << fixed << setprecision(2) << " " << setw(8) << numBG[n] << " & ";
		if(version==1||version==3) *fLogStream << fixed << setprecision(2) << "$" << LLeff[n] << " \\pm " << LLefferr[n] << "$" << " & ";//MT EFFICIENCY MISSING
		if(version==0||version==2) *fLogStream << fixed << setprecision(2) << "$" << (1.-LLeff[n])/(LLeff[n]*MTeff[n]) << " \\pm " << sqrt(pow(LLefferr[n]/(MTeff[n]*pow(LLeff[n],2)),2) + pow(MTefferr[n]*(1.-LLeff[n])/(LLeff[n]*pow(MTeff[n],2)),2)) << "$" << " & ";
		*fLogStream << "$" << " " << setw(9) << mctruth[n] << " \\pm " << " " << setw(6) <<  mctrutherr[n] << "$ & ";
		if(version==0||version==1) *fLogStream << fixed << setprecision(2) << "$" << " " << setw(9) << MCpred[n] << " \\pm " << " " << setw(7) << MCpred_syst_err[n] << "$" <<  " & ";
        	*fLogStream << fixed << setprecision(2) << " " << " " << setw(9) << datapred[n] << " $\\pm$ " << " " << setw(8) << datapred_stat_err[n] << " (stat) $\\pm$ " << " " << setw(8) << datapred_syst_err[n] << " (sys)"   << " \\\\" << endl;
	}
	*fLogStream << "\\hline\\hline"                                                                                         << endl
	     << "\\end{tabular}"                                                                                                << endl
	     << "\\end{center}"                                                                                                 << endl
	     << "\\end{table}"                                                                                                  << endl
	     << "\%ENDLATEX\%"                                                                                                  << endl
	     << endl;
	*fLogStream << endl << endl;
}

//this stores prediction into File, so that one can use it for prediction purposes
void PredictionFile(int version, vector<int> sr, vector<int> htr, vector<int> lepr, vector<double> MT2low, vector<double> MT2up, vector<double> numtrueW, vector<double> numtrueW_bg, vector<double> numData, vector<double> numBG, vector<double> LLeff, vector<double> LLefferr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruth, vector<double> mctrutherr, vector<double> MCpred, vector<double> MCpred_stat_err, vector<double> MCpred_syst_err, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, Bool_t PlotPrediction){

	TString filename = outputdir;
	TString plotdirectory = outputdir + "plots/";
	if(fRebin) plotdirectory = plotdirectory + "plots/";
	else       plotdirectory = plotdirectory + "plotsMT2binned/";
	if(PlotPrediction) Util::MakeOutputDir(plotdirectory);
	if(fRebin) filename = filename + "LostLeptonPredictionFile.root";
	else       filename = filename + "FineBinnedLostLeptonPredictionFile.root";

	map<string, TH1D*>    hs;
	//for now it is only truth +/- error
	//and prediction +/- error stored
	for(int i4 = 0; i4<leptontypesize;   ++i4){
	for(int i3 = 0; i3<HTbinsize;        ++i3){
		string htreg;
		if(i3==0 && fMET) htreg = "lowHT";
		if(i3==0 && fHT) htreg = "mediumHT";
		if(i3==1 && fHT) htreg = "highHT";
		if(fRebin){
			string mapname;
			mapname = "Prediction_"+lepton_type[i4]+"_"+htreg;
			if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", 9, 0, 9); hs[mapname]->Sumw2();
			hs[mapname]->GetXaxis()->SetTitle("signal region");
			hs[mapname]->GetXaxis()->SetBinLabel(1,"2j,0b");     hs[mapname]->GetXaxis()->SetBinLabel(2,"2j,#geq1b"); hs[mapname]->GetXaxis()->SetBinLabel(9,"#geq3j,#geq3b");
			hs[mapname]->GetXaxis()->SetBinLabel(3,"3-5j,0b");   hs[mapname]->GetXaxis()->SetBinLabel(4,"3-5j,1b");   hs[mapname]->GetXaxis()->SetBinLabel(5,"3-5j,2b");
			hs[mapname]->GetXaxis()->SetBinLabel(6,"#geq6j,0b"); hs[mapname]->GetXaxis()->SetBinLabel(7,"#geq6j,1b"); hs[mapname]->GetXaxis()->SetBinLabel(8,"#geq6j,2b");
			hs[mapname]->SetMarkerStyle(20);
			hs[mapname]->SetMarkerColor(kBlack);
			hs[mapname]->SetLineColor(kBlack);
			hs[mapname]->SetLineWidth(3);
			mapname = "Prediction_AllLep_"+htreg;
			if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", 9, 0, 9); hs[mapname]->Sumw2();
			hs[mapname]->GetXaxis()->SetTitle("signal region");
			hs[mapname]->GetXaxis()->SetBinLabel(1,"2j,0b");     hs[mapname]->GetXaxis()->SetBinLabel(2,"2j,#geq1b"); hs[mapname]->GetXaxis()->SetBinLabel(9,"#geq3j,#geq3b");
			hs[mapname]->GetXaxis()->SetBinLabel(3,"3-5j,0b");   hs[mapname]->GetXaxis()->SetBinLabel(4,"3-5j,1b");   hs[mapname]->GetXaxis()->SetBinLabel(5,"3-5j,2b");
			hs[mapname]->GetXaxis()->SetBinLabel(6,"#geq6j,0b"); hs[mapname]->GetXaxis()->SetBinLabel(7,"#geq6j,1b"); hs[mapname]->GetXaxis()->SetBinLabel(8,"#geq6j,2b");
			hs[mapname]->SetMarkerStyle(20);
			hs[mapname]->SetMarkerColor(kBlack);
			hs[mapname]->SetLineColor(kBlack);
			hs[mapname]->SetLineWidth(3);
			mapname = "SimulationTruth_"+lepton_type[i4]+"_"+htreg;
			if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", 9, 0, 9); hs[mapname]->Sumw2();
			hs[mapname]->GetXaxis()->SetTitle("signal region");
			hs[mapname]->GetXaxis()->SetBinLabel(1,"2j,0b");     hs[mapname]->GetXaxis()->SetBinLabel(2,"2j,#geq1b"); hs[mapname]->GetXaxis()->SetBinLabel(9,"#geq3j,#geq3b");
			hs[mapname]->GetXaxis()->SetBinLabel(3,"3-5j,0b");   hs[mapname]->GetXaxis()->SetBinLabel(4,"3-5j,1b");   hs[mapname]->GetXaxis()->SetBinLabel(5,"3-5j,2b");
			hs[mapname]->GetXaxis()->SetBinLabel(6,"#geq6j,0b"); hs[mapname]->GetXaxis()->SetBinLabel(7,"#geq6j,1b"); hs[mapname]->GetXaxis()->SetBinLabel(8,"#geq6j,2b");
			hs[mapname]->SetMarkerColor(kBlue);
			hs[mapname]->SetLineColor(kBlue);
			hs[mapname]->SetFillColor(kBlue);
			hs[mapname]->SetLineWidth(0);
			hs[mapname]->SetFillStyle(3002); 
			mapname = "SimulationTruth_AllLep_"+htreg;
			if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", 9, 0, 9); hs[mapname]->Sumw2();
			hs[mapname]->GetXaxis()->SetTitle("signal region");
			hs[mapname]->GetXaxis()->SetBinLabel(1,"2j,0b");     hs[mapname]->GetXaxis()->SetBinLabel(2,"2j,#geq1b"); hs[mapname]->GetXaxis()->SetBinLabel(9,"#geq3j,#geq3b");
			hs[mapname]->GetXaxis()->SetBinLabel(3,"3-5j,0b");   hs[mapname]->GetXaxis()->SetBinLabel(4,"3-5j,1b");   hs[mapname]->GetXaxis()->SetBinLabel(5,"3-5j,2b");
			hs[mapname]->GetXaxis()->SetBinLabel(6,"#geq6j,0b"); hs[mapname]->GetXaxis()->SetBinLabel(7,"#geq6j,1b"); hs[mapname]->GetXaxis()->SetBinLabel(8,"#geq6j,2b");
			hs[mapname]->SetMarkerColor(kBlue);
			hs[mapname]->SetLineColor(kBlue);
			hs[mapname]->SetFillColor(kBlue);
			hs[mapname]->SetLineWidth(0);
			hs[mapname]->SetFillStyle(3002); 
		}
	for(int i2 = 0; i2<signalregionsize; ++i2){
		if(fRebin&&i4==0&&i3==0){
			string sigreg;
			if(i2==0) sigreg="2j0b";
			if(i2==1) sigreg="2j1to2b";
			if(i2==2) sigreg="3to5j0b";
			if(i2==3) sigreg="3to5j1b";
			if(i2==4) sigreg="3to5j2b";
			if(i2==5) sigreg="6j0b";
			if(i2==6) sigreg="6j1b";
			if(i2==7) sigreg="6j2b";
			if(i2==8) sigreg="3b";
			string mapname;
			mapname = "Prediction_"+sigreg;
			if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", 9, 0, 9); hs[mapname]->Sumw2();
			hs[mapname]->GetXaxis()->SetTitle("signal region");
			hs[mapname]->GetXaxis()->SetBinLabel(1,"e,lH_{T}");    hs[mapname]->GetXaxis()->SetBinLabel(2,"e,mH_{T}");    hs[mapname]->GetXaxis()->SetBinLabel(3,"e,hH_{T}");
			hs[mapname]->GetXaxis()->SetBinLabel(4,"#mu,lH_{T}");  hs[mapname]->GetXaxis()->SetBinLabel(5,"#mu,mH_{T}");  hs[mapname]->GetXaxis()->SetBinLabel(6,"#mu,hH_{T}");
			hs[mapname]->GetXaxis()->SetBinLabel(7,"#tau,lH_{T}"); hs[mapname]->GetXaxis()->SetBinLabel(8,"#tau,mH_{T}"); hs[mapname]->GetXaxis()->SetBinLabel(9,"#tau,hH_{T}");
			hs[mapname]->SetMarkerStyle(20);
			hs[mapname]->SetMarkerColor(kBlack);
			hs[mapname]->SetLineColor(kBlack);
			hs[mapname]->SetLineWidth(3);
			mapname = "SimulationTruth_"+sigreg;
			if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", 9, 0, 9); hs[mapname]->Sumw2();
			hs[mapname]->GetXaxis()->SetTitle("signal region");
			hs[mapname]->GetXaxis()->SetBinLabel(1,"e,lH_{T}");    hs[mapname]->GetXaxis()->SetBinLabel(2,"e,mH_{T}");    hs[mapname]->GetXaxis()->SetBinLabel(3,"e,hH_{T}");
			hs[mapname]->GetXaxis()->SetBinLabel(4,"#mu,lH_{T}");  hs[mapname]->GetXaxis()->SetBinLabel(5,"#mu,mH_{T}");  hs[mapname]->GetXaxis()->SetBinLabel(6,"#mu,hH_{T}");
			hs[mapname]->GetXaxis()->SetBinLabel(7,"#tau,lH_{T}"); hs[mapname]->GetXaxis()->SetBinLabel(8,"#tau,mH_{T}"); hs[mapname]->GetXaxis()->SetBinLabel(9,"#tau,hH_{T}");
			hs[mapname]->SetMarkerColor(kBlue);
			hs[mapname]->SetLineColor(kBlue);
			hs[mapname]->SetFillColor(kBlue);
			hs[mapname]->SetLineWidth(0);
			hs[mapname]->SetFillStyle(3002); 
		}
	}
	}}
	string oldlep = "dummy";
	string oldsr = "dummy";
	string oldHT = "dummy";
	vector<double> dp, dpe, st, ste, mt2binning; dp.clear(); dpe.clear(); st.clear(); ste.clear(); mt2binning.clear();
	for(unsigned int n = 0; n<sr.size(); ++n){
		string lep;
		if(lepr[n]==0) lep="Muo";
		else if(lepr[n]==1) lep="Ele";
		else continue;
		string sigreg;
		if(sr[n]==0) sigreg="2j0b";
		if(sr[n]==1) sigreg="2j1to2b";
		if(sr[n]==2) sigreg="3to5j0b";
		if(sr[n]==3) sigreg="3to5j1b";
		if(sr[n]==4) sigreg="3to5j2b";
		if(sr[n]==5) sigreg="6j0b";
		if(sr[n]==6) sigreg="6j1b";
		if(sr[n]==7) sigreg="6j2b";
		if(sr[n]==8) sigreg="3b";
		string htreg;
		if(htr[n]==0 && fMET) htreg = "lowHT";
		if(htr[n]==1 && fHT) htreg = "mediumHT";
		if(htr[n]==2 && fHT) htreg = "highHT";
		int mt2lastbin = 9999;
		if(sr[n]==0) { if(htr[n]==0) mt2lastbin = 750; if(htr[n]==1) mt2lastbin = 1000; if(htr[n]==2) mt2lastbin =  900;}
		if(sr[n]==1) { if(htr[n]==0) mt2lastbin = 700; if(htr[n]==1) mt2lastbin =  700; if(htr[n]==2) mt2lastbin =  350;}
		if(sr[n]==2) { if(htr[n]==0) mt2lastbin = 750; if(htr[n]==1) mt2lastbin = 1000; if(htr[n]==2) mt2lastbin = 1000;}
		if(sr[n]==3) { if(htr[n]==0) mt2lastbin = 700; if(htr[n]==1) mt2lastbin =  900; if(htr[n]==2) mt2lastbin =  550;}
		if(sr[n]==4) { if(htr[n]==0) mt2lastbin = 550; if(htr[n]==1) mt2lastbin =  500; if(htr[n]==2) mt2lastbin =  350;}
		if(sr[n]==5) { if(htr[n]==0) mt2lastbin = 520; if(htr[n]==1) mt2lastbin =  600; if(htr[n]==2) mt2lastbin =  500;}
		if(sr[n]==6) { if(htr[n]==0) mt2lastbin = 450; if(htr[n]==1) mt2lastbin =  500; if(htr[n]==2) mt2lastbin =  500;}
		if(sr[n]==7) { if(htr[n]==0) mt2lastbin = 400; if(htr[n]==1) mt2lastbin =  450; if(htr[n]==2) mt2lastbin =  350;}
		if(sr[n]==8) { if(htr[n]==0) mt2lastbin = 400; if(htr[n]==1) mt2lastbin =  450; if(htr[n]==2) mt2lastbin =  300;}
		if(htreg!=oldHT) oldHT = htreg;
		if(lep!=oldlep) oldlep = lep;
		if(fRebin){
			string mapname = "Prediction_"+lep+"_"+htreg;
			hs[mapname]->SetBinContent(sr[n]+1, datapred[n]);
			hs[mapname]->SetBinError(sr[n]+1, sqrt(pow(datapred_stat_err[n],2)+pow(datapred_syst_err[n],2)) );
			mapname = "SimulationTruth_"+lep+"_"+htreg;
			hs[mapname]->SetBinContent(sr[n]+1, mctruth[n]);
			hs[mapname]->SetBinError(sr[n]+1, mctrutherr[n]);
		} else{
			mt2binning.push_back(MT2low[n]);
			dp.push_back(datapred[n]);
			dpe.push_back(sqrt(pow(datapred_stat_err[n],2)+pow(datapred_syst_err[n],2)) );
			st.push_back(mctruth[n]);
			ste.push_back(mctrutherr[n]);
			if(MT2up[n]==10000.) { 
				mt2binning.push_back(mt2lastbin);
				int Numbins = mt2binning.size()-1;
				int numarray = mt2binning.size();
				double binns[numarray];
				for(int ii = 0; ii<=Numbins; ++ii) binns[ii] = mt2binning[ii];
				string mapname = "Prediction_"+lep+"_"+htreg+"_"+sigreg;
				if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", Numbins, binns); hs[mapname]->Sumw2();
				hs[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); hs[mapname]->GetYaxis()->SetTitle("events");
				hs[mapname]->SetMarkerStyle(20);
				hs[mapname]->SetMarkerColor(kBlack);
				hs[mapname]->SetLineColor(kBlack);
				hs[mapname]->SetLineWidth(3);
				for(int ii = 1; ii<=Numbins; ++ii){
					hs[mapname]->SetBinContent(ii, dp[ii-1]);
					hs[mapname]->SetBinError(ii, dpe[ii-1]);
				}
				mapname = "Prediction_AllLep_"+htreg+"_"+sigreg;
				if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", Numbins, binns); hs[mapname]->Sumw2();
				hs[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); hs[mapname]->GetYaxis()->SetTitle("events");
				hs[mapname]->SetMarkerStyle(20);
				hs[mapname]->SetMarkerColor(kBlack);
				hs[mapname]->SetLineColor(kBlack);
				hs[mapname]->SetLineWidth(3);
				mapname = "SimulationTruth_"+lep+"_"+htreg+"_"+sigreg;
				if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", Numbins, binns); hs[mapname]->Sumw2();
				hs[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); hs[mapname]->GetYaxis()->SetTitle("events");
				hs[mapname]->SetMarkerColor(kBlue);
				hs[mapname]->SetLineColor(kBlue);
				hs[mapname]->SetFillColor(kBlue);
				hs[mapname]->SetLineWidth(0);
				hs[mapname]->SetFillStyle(3002); 
				for(int ii = 1; ii<=Numbins; ++ii){
					hs[mapname]->SetBinContent(ii, st[ii-1]);
					hs[mapname]->SetBinError(ii, ste[ii-1]);
				}
				mapname = "SimulationTruth_AllLep_"+htreg+"_"+sigreg;
				if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", Numbins, binns); hs[mapname]->Sumw2();
				hs[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); hs[mapname]->GetYaxis()->SetTitle("events");
				hs[mapname]->SetMarkerColor(kBlue);
				hs[mapname]->SetLineColor(kBlue);
				hs[mapname]->SetFillColor(kBlue);
				hs[mapname]->SetLineWidth(0);
				hs[mapname]->SetFillStyle(3002); 
				mt2binning.clear();
				dp.clear(); dpe.clear();
				st.clear(); ste.clear();
			}
		}
	}
	//adding ele+muo to alllep
	for(int i3 = 0; i3<HTbinsize;        ++i3){
		string htreg;
		if(i3==0 && fMET) htreg = "lowHT";
		if(i3==0 && fHT) htreg = "mediumHT";
		if(i3==1 && fHT) htreg = "highHT";
		if(fRebin){
			string hse = "Prediction_Ele_"+htreg;
			string hsm = "Prediction_Muo_"+htreg;
			string hsa = "Prediction_AllLep_"+htreg;
			hs[hsa]->Add(hs[hse]);
			hs[hsa]->Add(hs[hsm]);
			hse = "SimulationTruth_Ele_"+htreg;
			hsm = "SimulationTruth_Muo_"+htreg;
			hsa = "SimulationTruth_AllLep_"+htreg;
			hs[hsa]->Add(hs[hse]);
			hs[hsa]->Add(hs[hsm]);
		}
		for(int i2 = 0; i2<signalregionsize; ++i2){
			string sigreg;
			if(i2==0) sigreg="2j0b";
			if(i2==1) sigreg="2j1to2b";
			if(i2==2) sigreg="3to5j0b";
			if(i2==3) sigreg="3to5j1b";
			if(i2==4) sigreg="3to5j2b";
			if(i2==5) sigreg="6j0b";
			if(i2==6) sigreg="6j1b";
			if(i2==7) sigreg="6j2b";
			if(i2==8) sigreg="3b";
			if(!fRebin){
				string hse = "Prediction_Ele_"+htreg+"_"+sigreg;
				string hsm = "Prediction_Muo_"+htreg+"_"+sigreg;
				string hsa = "Prediction_AllLep_"+htreg+"_"+sigreg;
				hs[hsa]->Add(hs[hse]);
				hs[hsa]->Add(hs[hsm]);
				hse = "SimulationTruth_Ele_"+htreg+"_"+sigreg;
				hsm = "SimulationTruth_Muo_"+htreg+"_"+sigreg;
				hsa = "SimulationTruth_AllLep_"+htreg+"_"+sigreg;
				hs[hsa]->Add(hs[hse]);
				hs[hsa]->Add(hs[hsm]);
			}
		}
	}

	//now store all the shit
   	TFile *fpredictionfile = new TFile(filename.Data(),"RECREATE");
	fpredictionfile->cd();
	for(map<string,TH1D*>::iterator h=hs.begin(); h!=hs.end();++h){
		h->second->Write();
	}
	fpredictionfile->Close();
	*fLogStream << "Saved prediction/mc truth histograms in " << filename.Data() << endl;


	if(PlotPrediction){
		//format histogram
		for(map<string,TH1D*>::iterator h=hs.begin(); h!=hs.end();++h){
			TString helperstring = h->first;
			h->second->SetStats(0);
			h->second->SetLineStyle(0);
			h->second->GetXaxis()->SetLabelOffset(0.007);
			h->second->GetXaxis()->SetLabelSize(0.05);
			h->second->GetXaxis()->SetTitleSize(0.06);
			h->second->GetXaxis()->SetTitleOffset(0.9);
			h->second->GetXaxis()->SetTitleFont(42);
			h->second->GetYaxis()->SetTitle("Events");
			h->second->GetYaxis()->SetLabelFont(42);
			h->second->GetYaxis()->SetLabelOffset(0.007);
			h->second->GetYaxis()->SetLabelSize(0.05);
			h->second->GetYaxis()->SetTitleSize(0.06);
			h->second->GetYaxis()->SetTitleOffset(1.25);
			h->second->GetYaxis()->SetTitleFont(42);
			h->second->GetZaxis()->SetLabelFont(42);
			h->second->GetZaxis()->SetLabelOffset(0.007);
			h->second->GetZaxis()->SetLabelSize(0.05);
			h->second->GetZaxis()->SetTitleSize(0.06);
			h->second->GetZaxis()->SetTitleFont(42);
			h->second->GetXaxis()->SetLabelFont(42);
//			if(h->first[0]==(string)"S"){//simulation truth gets plotted first
			if(helperstring.Contains("Simulation")){//simulation truth gets plotted first
			h->second->SetMarkerStyle(1);
			h->second->SetMarkerColor(kBlue);
			h->second->SetLineColor(kBlue);
			h->second->SetFillColor(kBlue);
			h->second->SetLineWidth(0);
			h->second->SetFillStyle(3002); 
			if(helperstring.Contains("b") && fRebin){
				h->second->GetXaxis()->SetTitle("signal region");
				h->second->GetXaxis()->SetBinLabel(1,"e,lH_{T}");
				h->second->GetXaxis()->SetBinLabel(2,"e,mH_{T}");
				h->second->GetXaxis()->SetBinLabel(3,"e,hH_{T}");
				h->second->GetXaxis()->SetBinLabel(4,"#mu,lH_{T}");
				h->second->GetXaxis()->SetBinLabel(5,"#mu,mH_{T}");
				h->second->GetXaxis()->SetBinLabel(6,"#mu,hH_{T}");
				h->second->GetXaxis()->SetBinLabel(7,"#tau,lH_{T}");
				h->second->GetXaxis()->SetBinLabel(8,"#tau,mH_{T}");
				h->second->GetXaxis()->SetBinLabel(9,"#tau,hH_{T}");
		//	} else if(h->first[h->first.size()]==(string)"b" && !fRebin){
			} else if(helperstring.Contains("b") && !fRebin){
				h->second->GetXaxis()->SetTitle("M_{T2} [GeV]");
			} else{//fRebin and last letter is a "T"
				h->second->GetXaxis()->SetTitle("signal region");
				h->second->GetXaxis()->SetBinLabel(1,"2j,0b");
				h->second->GetXaxis()->SetBinLabel(2,"2j,#geq1b");
				h->second->GetXaxis()->SetBinLabel(3,"3-5j,0b");
				h->second->GetXaxis()->SetBinLabel(4,"3-5j,1b");
				h->second->GetXaxis()->SetBinLabel(5,"3-5j,2b");
				h->second->GetXaxis()->SetBinLabel(6,"#geq6j,0b");
				h->second->GetXaxis()->SetBinLabel(7,"#geq6j,1b");
				h->second->GetXaxis()->SetBinLabel(8,"#geq6j,2b");
				h->second->GetXaxis()->SetBinLabel(9,"#geq3j,#geq3b");
			}
			} else{
			h->second->SetMarkerStyle(20);
			h->second->SetMarkerColor(kBlack);
			h->second->SetLineColor(kBlack);
			h->second->SetLineWidth(3);
			}
		}
		TCanvas *c1 = new TCanvas("c1", "",60,22,600,600);
		gStyle->SetOptFit(1);
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		c1->SetFillColor(0);
		c1->SetBorderMode(0);
		c1->SetBorderSize(2);
		c1->SetTickx(1);
		c1->SetTicky(1);
		c1->SetLeftMargin(0.18);
		c1->SetRightMargin(0.05);
		c1->SetTopMargin(0.07);
		c1->SetBottomMargin(0.15);
		c1->SetFrameFillStyle(0);
		c1->SetFrameBorderMode(0);
		c1->SetFrameFillStyle(0);
		c1->SetFrameBorderMode(0);
		c1->cd();

		TLatex *   tex = new TLatex(0.328859,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV");
		tex->SetNDC();
		tex->SetTextFont(42);
		tex->SetTextSize(0.04181185);
		tex->SetLineWidth(2);
		TLatex toplep;
		toplep.SetNDC();
		toplep.SetTextAlign(31);
		toplep.SetTextFont(42);
		toplep.SetTextSize(0.04181185);
		toplep.SetLineWidth(2);
		TLatex ht;
		ht.SetNDC();
		ht.SetTextAlign(31);
		ht.SetTextFont(42);
		ht.SetTextSize(0.04181185);
		ht.SetLineWidth(2);

		TLegend *leg = new TLegend(0.6252416,0.7657343,0.824906,0.9003497,NULL,"brNDC");
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04181185);
		leg->SetTextFont(42);
		leg->SetLineColor(1);
		leg->SetLineStyle(1);
		leg->SetLineWidth(2);
		leg->SetFillColor(0);
		leg->SetFillStyle(1001);
   		TLegendEntry *entry=leg->AddEntry("NULL","simulation truth","f");
		entry->SetMarkerColor(kBlue);
		entry->SetLineColor(kBlue);
		entry->SetFillColor(kBlue);
		entry->SetLineWidth(0);
		entry->SetFillStyle(3002); 
   		TLegendEntry *entry2=leg->AddEntry("NULL","data prediction","p");
		entry2->SetLineColor(1);
		entry2->SetLineStyle(1);
		entry2->SetLineWidth(2);
		entry2->SetMarkerColor(1);
		entry2->SetMarkerStyle(20);
		entry2->SetMarkerSize(1);

		string hsp, hss;
		TString texttoplep, textht;
		string outname;
		double max(0.); double min(99999.);
		double maxp(0.); double minp(99999.);
		double maxs(0.); double mins(99999.);
		for(int i3 = 0; i3<HTbinsize;        ++i3){
			string htreg;
			if(i3==0 && fMET) { htreg = "lowHT"; textht = "low H_{T}"; }
			if(i3==0 && fHT)  { htreg = "mediumHT"; textht = "medium H_{T}"; }
			if(i3==1 && fHT)  { htreg = "highHT"; textht = "high H_{T}"; }
			if(fRebin){
				hsp = "Prediction_Ele_"+htreg;
				hss = "SimulationTruth_Ele_"+htreg;
				texttoplep = "1 electron";
				min = 0.;
				maxp = hs[hsp]->GetBinContent(hs[hsp]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hsp]->GetMaximumBin());
				maxs = hs[hss]->GetBinContent(hs[hss]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hss]->GetMaximumBin());
				max  = (maxp>maxs)?maxp:maxs;
				max = 1.5*max;
				hs[hss]->SetMaximum(max);
				c1->Clear();
				c1->cd();
				hs[hss]->SetFillStyle(3013);
				hs[hss]->Draw("E2");
				hs[hsp]->Draw("sameE");
				tex->Draw();
				leg->Draw();
				toplep.DrawLatex(0.6,0.8548951,texttoplep.Data());
				ht.DrawLatex(0.6,0.791958,textht.Data());
				outname = plotdirectory+hsp + ".eps";
				c1->SaveAs(outname.c_str());
				hsp = "Prediction_Muo_"+htreg;
				hss = "SimulationTruth_Muo_"+htreg;
				texttoplep = "1 muon";
				min = 0.;
				maxp = hs[hsp]->GetBinContent(hs[hsp]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hsp]->GetMaximumBin());
				maxs = hs[hss]->GetBinContent(hs[hss]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hss]->GetMaximumBin());
				max  = (maxp>maxs)?maxp:maxs;
				max = 1.5*max;
				hs[hss]->SetMaximum(max);
				c1->Clear();
				c1->cd();
				hs[hss]->SetFillStyle(3013);
				hs[hss]->Draw("E2");
				hs[hsp]->Draw("sameE");
				tex->Draw();
				leg->Draw();
				toplep.DrawLatex(0.6,0.8548951,texttoplep.Data());
				ht.DrawLatex(0.6,0.791958,textht.Data());
				outname = plotdirectory+hsp + ".eps";
				c1->SaveAs(outname.c_str());
				hsp = "Prediction_AllLep_"+htreg;
				hss = "SimulationTruth_AllLep_"+htreg;
				texttoplep = "1 lepton";
				min = 0.;
				maxp = hs[hsp]->GetBinContent(hs[hsp]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hsp]->GetMaximumBin());
				maxs = hs[hss]->GetBinContent(hs[hss]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hss]->GetMaximumBin());
				max  = (maxp>maxs)?maxp:maxs;
				max = 1.5*max;
				hs[hss]->SetMaximum(max);
				c1->Clear();
				c1->cd();
				hs[hss]->SetFillStyle(3013);
				hs[hss]->Draw("E2");
				hs[hsp]->Draw("sameE");
				tex->Draw();
				leg->Draw();
				toplep.DrawLatex(0.6,0.8548951,texttoplep.Data());
				ht.DrawLatex(0.6,0.791958,textht.Data());
				outname = plotdirectory+hsp + ".eps";
				c1->SaveAs(outname.c_str());
			}
			for(int i2 = 0; i2<signalregionsize; ++i2){
				string sigreg; TString regsig;
				if(i2==0) { sigreg="2j0b";    regsig="2j, 0b"; }
				if(i2==1) { sigreg="2j1to2b"; regsig="2j, #geq1b"; }
				if(i2==2) { sigreg="3to5j0b"; regsig="3-5j, 0b"; }
				if(i2==3) { sigreg="3to5j1b"; regsig="3-5j, 1b"; }
				if(i2==4) { sigreg="3to5j2b"; regsig="3-5j, 2b"; }
				if(i2==5) { sigreg="6j0b";    regsig="#geq6j, 0b"; }
				if(i2==6) { sigreg="6j1b";    regsig="#geq6j, 1b"; }
				if(i2==7) { sigreg="6j2b";    regsig="#geq6j, 2b"; }
				if(i2==8) { sigreg="3b";      regsig="#geq3j, #geq3b"; }
				if(!fRebin){
					hsp = "Prediction_Ele_"+htreg+"_"+sigreg;
					hss = "SimulationTruth_Ele_"+htreg+"_"+sigreg;
					texttoplep = "1 e, "+ regsig;
					min = 0.;
					maxp = hs[hsp]->GetBinContent(hs[hsp]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hsp]->GetMaximumBin());
					maxs = hs[hss]->GetBinContent(hs[hss]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hss]->GetMaximumBin());
					max  = (maxp>maxs)?maxp:maxs;
					max = 1.25*max;
					hs[hss]->SetMaximum(max);
					c1->Clear();
					c1->cd();
					hs[hss]->SetFillStyle(3013);
					hs[hss]->Draw("E2");
					hs[hsp]->Draw("sameE");
					tex->Draw();
					leg->Draw();
					toplep.DrawLatex(0.6,0.8548951,texttoplep.Data());
					ht.DrawLatex(0.6,0.791958,textht.Data());
					outname = plotdirectory+hsp + ".eps";
					c1->SaveAs(outname.c_str());
					hsp = "Prediction_Muo_"+htreg+"_"+sigreg;
					hss = "SimulationTruth_Muo_"+htreg+"_"+sigreg;
					texttoplep = "1 #mu, "+ regsig;
					min = 0.;
					maxp = hs[hsp]->GetBinContent(hs[hsp]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hsp]->GetMaximumBin());
					maxs = hs[hss]->GetBinContent(hs[hss]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hss]->GetMaximumBin());
					max  = (maxp>maxs)?maxp:maxs;
					max = 1.25*max;
					hs[hss]->SetMaximum(max);
					c1->Clear();
					c1->cd();
					hs[hss]->SetFillStyle(3013);
					hs[hss]->Draw("E2");
					hs[hsp]->Draw("sameE");
					tex->Draw();
					leg->Draw();
					toplep.DrawLatex(0.6,0.8548951,texttoplep.Data());
					ht.DrawLatex(0.6,0.791958,textht.Data());
					outname = plotdirectory+hsp + ".eps";
					c1->SaveAs(outname.c_str());
					hsp = "Prediction_AllLep_"+htreg+"_"+sigreg;
					hss = "SimulationTruth_AllLep_"+htreg+"_"+sigreg;
					texttoplep = "1 l, "+ regsig;
					min = 0.;
					maxp = hs[hsp]->GetBinContent(hs[hsp]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hsp]->GetMaximumBin());
					maxs = hs[hss]->GetBinContent(hs[hss]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hss]->GetMaximumBin());
					max  = (maxp>maxs)?maxp:maxs;
					max = 1.25*max;
					hs[hss]->SetMaximum(max);
					c1->Clear();
					c1->cd();
					hs[hss]->SetFillStyle(3013);
					hs[hss]->Draw("E2");
					hs[hsp]->Draw("sameE");
					tex->Draw();
					leg->Draw();
					toplep.DrawLatex(0.6,0.8548951,texttoplep.Data());
					ht.DrawLatex(0.6,0.791958,textht.Data());
					outname = plotdirectory+hsp + ".eps";
					c1->SaveAs(outname.c_str());
				}
				if(fRebin && i3==0){
					if(i2==0) { sigreg="2j0b";    regsig="2 jets, 0 b jets"; }
					if(i2==1) { sigreg="2j1to2b"; regsig="2 jets, #geq1 b jets"; }
					if(i2==2) { sigreg="3to5j0b"; regsig="3-5 jets, 0 b jets"; }
					if(i2==3) { sigreg="3to5j1b"; regsig="3-5 jets, 1 b jet"; }
					if(i2==4) { sigreg="3to5j2b"; regsig="3-5 jets, 2 b jets"; }
					if(i2==5) { sigreg="6j0b";    regsig="#geq6 jets, 0 b jets"; }
					if(i2==6) { sigreg="6j1b";    regsig="#geq6 jets, 1 b jet"; }
					if(i2==7) { sigreg="6j2b";    regsig="#geq6 jets, 2 b jets"; }
					if(i2==8) { sigreg="3b";      regsig="#geq3 jets, #geq 3b jets"; }
					hsp = "Prediction_"+sigreg;
					hss = "SimulationTruth_"+sigreg;
					texttoplep = regsig;
					min = 0.;
					maxp = hs[hsp]->GetBinContent(hs[hsp]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hsp]->GetMaximumBin());
					maxs = hs[hss]->GetBinContent(hs[hss]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hss]->GetMaximumBin());
					max  = (maxp>maxs)?maxp:maxs;
					max = 1.5*max;
					hs[hss]->SetMaximum(max);
					c1->Clear();
					c1->cd();
					hs[hss]->SetFillStyle(3013);
					hs[hss]->Draw("E2");
					hs[hsp]->Draw("sameE");
					tex->Draw();
					leg->Draw();
					toplep.DrawLatex(0.6,0.8548951,texttoplep.Data());
					outname = plotdirectory+hsp + ".eps";
					c1->SaveAs(outname.c_str());	
				}
			}
		}
	}
}

//this function just makes the one lepton yield table (usually with MT cut applied for one lepton selection, unless flags are set differently)
//it contains no comments as this function is from style very similar to SummaryTable(...)
void YieldTable(vector<int> sr, vector<int> htr, vector<int> lepr, vector<double> MT2low, vector<double> MT2up, vector<double> numQCD, vector<double> numZ, vector<double> numW, vector<double> numT, vector<double> numOther, vector<double> numMC, vector<double> BGerr, vector<double> numData){

	*fLogStream << "*********************************************************************" << endl;
	*fLogStream << "\%BEGINLATEX\%"              << endl;
	*fLogStream << "\\begin{table}"              << endl
	     << "\\begin{center}"             << endl
	     << "\\small"                     << endl
             << "\\begin{tabular}{lccccccc}"  << endl
	     << "\\hline\\hline"              << endl;

	if(!fRebin) *fLogStream << "$M_{T2}$ (GeV)";
	else        *fLogStream << "signal region       ";
	*fLogStream  << " & $" << "N^{QCD}" << "$ & $" << "N^{Z}" << "$  & $" << "N^{W}" << "$  & $" << "N^{Top}" << "$  & $" << "N^{Other}" << "$ &           $"  << "N^{MC}" << "$      & $" << "N^{data}" << "$ ";
	*fLogStream << "\\\\" << endl << "\\hline\\hline"             << endl;
	string oldlep = "dummy";
	string oldsr = "dummy";
	string oldHT = "dummy";
	for(unsigned int n = 0; n<sr.size(); ++n){
		string lep;
		if(lepr[n]==0) lep="Muo";
		else if(lepr[n]==1) lep="Ele";
		else continue;
		string sigreg;
		if(sr[n]==0) sigreg="2j, 0b";
		if(sr[n]==1) sigreg="2j, $\\ge 1$b";
		if(sr[n]==2) sigreg="$3-5$j, 0b";
		if(sr[n]==3) sigreg="$3-5$j, 1b";
		if(sr[n]==4) sigreg="$3-5$j, 2b";
		if(sr[n]==5) sigreg="$\\ge 6$j, 0b";
		if(sr[n]==6) sigreg="$\\ge 6$j, 1b";
		if(sr[n]==7) sigreg="$\\ge 6$j, 2b";
		if(sr[n]==8) sigreg="$\\ge 3$j, $\\ge 3$b";
		string htreg;
		if(htr[n]==0 && fMET) htreg = "450 GeV $\\leq H_{T} < 750$ GeV";
		if(htr[n]==1 && fHT) htreg = "750 GeV $\\leq H_{T} < 1200$ GeV";
		if(htr[n]==2 && fHT) htreg = "$H_{T}\\ge 1200$ GeV";
		if(htreg!=oldHT){
			*fLogStream << " \\hline " << endl << "\\multicolumn{7}{l}{" <<  htreg << "} \\\\ " << endl << "\\hline" << endl;
			oldHT = htreg;
		}
		if(lep!=oldlep){
			*fLogStream << " \\hline " << endl << lep << "\\\\ " << endl << "\\hline" << endl;
			oldlep = lep;
		}
		if(!fRebin){
		if(sigreg!=oldsr){
			*fLogStream << " \\hline  " << endl << sigreg << "\\\\" << endl;
			oldsr = sigreg;
		}
		if(MT2up[n]==10000.) *fLogStream << "$" << int(MT2low[n]) << "-" << "\\infty$" << " " << setw(4) << " & ";
		else            *fLogStream << "$" << int(MT2low[n]) << "-" << int(MT2up[n]) << "$" << " " << setw(7) << " & ";
		}
		else *fLogStream << " " << setw(18) << sigreg << " & ";
		*fLogStream << fixed << setprecision(2)
		<< " " << setw(7) << numQCD[n] << " & " << " " << setw(7) << numZ[n] << " & " << " " << setw(7) << numW[n] << " &  " << " " << setw(8) << numT[n] << " & " << " " << setw(10) << numOther[n] << " & " << " " << setw(10) 
		<< numMC[n] << "$\\pm" << " " << setw(7) << BGerr[n] << "$ & ";
		*fLogStream << " " << setw(10) << int(numData[n]) << " \\\\" << endl;
	}
	*fLogStream << "\\hline\\hline"                                                                                            << endl
		<< "\\end{tabular}"                                                                                                << endl
		<< "\\end{center}"                                                                                                 << endl
		<< "\\end{table}"                                                                                                  << endl
		<< "\%ENDLATEX\%"                                                                                                  << endl
		<< endl;
	*fLogStream << endl << endl;
}

//this function makes a compact card of signal region vs MCtruth +/- err vs. pred. +/- err
//it contains no comments as this function is from style very similar to SummaryTable(...)
void PredictionCard(vector<int> sr, vector<int> htr, vector<int> lepr, vector<int> MT2bin, vector<double> mctruth, vector<double> mctrutherr, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, vector<double> LLeff, vector<double> LLefferr, vector<double> MTeff, vector<double> MTefferr){
	
	*fLogStream << "*********************************************************************" << endl;
	*fLogStream << " Name     " << " Lep  " << " SR " << " HTbin " << " MT2bin   " << " MCpred   " << " DataDrivenPred " << " DataDrivenPredError " << " ScaleFactor " << " ScaleFactorError " << endl;
	for(unsigned int n = 0; n<sr.size(); ++n){
		string lep;
		if(lepr[n]==0) lep="Muo";
		else if(lepr[n]==1) lep="Ele";
		else continue;
		if(LLeff[n]!=0){
		*fLogStream << "LostLepton " << lep << " " << setw(3) << sr[n] << " " << setw(5) << htr[n] << " " << setw(6) << MT2bin[n]-1 << " " << setw(14) << fixed << setprecision(4) << mctruth[n] <<  " " << setw(12) << datapred[n] << " " << setw(17) << sqrt( pow(datapred_stat_err[n],2) + pow(datapred_syst_err[n],2) );
		*fLogStream << " " << setw(15) << fixed << setprecision(4) << (1.-LLeff[n])/(LLeff[n]*MTeff[n]) << " " << setw(14) << sqrt(pow(LLefferr[n]/(MTeff[n]*pow(LLeff[n],2)),2) + pow(MTefferr[n]*(1.-LLeff[n])/(LLeff[n]*pow(MTeff[n],2)),2)) << endl;
		}
		else{
		*fLogStream << "LostLepton " << lep << " " << setw(3) << sr[n] << " " << setw(5) << htr[n] << " " << setw(6) << MT2bin[n]-1 << " " << setw(14) << fixed << setprecision(4) << mctruth[n] << " " << setw(12) << "-9" << " " << setw(17) << "-9";
            	*fLogStream << " " << setw(15) << "-9" << " " << setw(14) << "-9" << endl;
		}
	}
	*fLogStream << endl << endl;

}

//this function makes a compact card of signal region vs MCtruth +/- err vs. pred. +/- err, and also prints out (from MC truth) the fraction of expected lost lepton from W,Top,rest(dibosons?)
//it contains no comments as this function is from style very similar to SummaryTable(...)
void PredictionCardSplitted(vector<int> sr, vector<int> htr, vector<int> lepr, vector<int> MT2bin, vector<double> mctruth, vector<double> mctrutherr, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, vector<double> LLeff, vector<double> LLefferr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruthfractionW, vector<double> mctruthfractionTT, vector<double> numW, vector<double> numTT, vector<double> numMC){
	
	*fLogStream << "*********************************************************************" << endl;
	*fLogStream << " Name     " << " Lep  " << " SR " << " HTbin " << " MT2bin   " << " MCpred   " << " DataDrivenPred " << " DataDrivenPredError " << " ScaleFactor " << " ScaleFactorError      " << " W:TTbar:Other(truth)          " << " W:TTbar:Other(lepton yield) " << endl;
	for(unsigned int n = 0; n<sr.size(); ++n){
		string lep;
		if(lepr[n]==0) lep="Muo";
		else if(lepr[n]==1) lep="Ele";
		else continue;
		if(LLeff[n]!=0){
		*fLogStream << "LostLepton " << lep << " " << setw(3) << sr[n] << " " << setw(5) << htr[n] << " " << setw(6) << MT2bin[n]-1 << " " << setw(14) << fixed << setprecision(4) << mctruth[n] << " " << setw(12) << datapred[n] << " " << setw(17) << sqrt( pow(datapred_stat_err[n],2) + pow(datapred_syst_err[n],2) );
		*fLogStream << " " << setw(15) << fixed << setprecision(4) << (1.-LLeff[n])/(LLeff[n]*MTeff[n]) << " " << setw(14) << sqrt(pow(LLefferr[n]/(MTeff[n]*pow(LLeff[n],2)),2) + pow(MTefferr[n]*(1.-LLeff[n])/(LLeff[n]*pow(MTeff[n],2)),2));
		*fLogStream << " " << setw(16) << fixed << setprecision(4) << mctruthfractionW[n]/mctruth[n]<< " : " <<mctruthfractionTT[n]/mctruth[n]<< " : " << (mctruth[n]-mctruthfractionW[n]-mctruthfractionTT[n])/mctruth[n] << " " << setw(15) << numW[n]/numMC[n]<< " : " <<numTT[n]/numMC[n]<< " : " <<(numMC[n]-numW[n]-numTT[n])/numMC[n] << endl;
		}
		else{
		*fLogStream << "LostLepton " << lep << " " << setw(3) << sr[n] << " " << setw(5) << htr[n] << " " << setw(6) << MT2bin[n]-1 << " " << setw(14) << fixed << setprecision(4) << mctruth[n] << " " << setw(12) << "-9" << " " << setw(17) << "-9";
            	*fLogStream << " " << setw(15) << "-9" << " " << setw(14) << "-9";
		*fLogStream << " " << setw(16) << fixed << setprecision(4) << mctruthfractionW[n]/mctruth[n]<<" : "<<mctruthfractionTT[n]/mctruth[n]<<" : "<<(mctruth[n]-mctruthfractionW[n]-mctruthfractionTT[n])/mctruth[n] << " " << setw(26) << "- : - : -" << endl;
		}
	}
	*fLogStream << endl << endl;

}

//standard load function that reads out the samples.dat into the fSamples vector
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
