// #ifndef __CINT__
// #include "RooGlobalFunc.h"
// #endif
// #include "RooFit.h"
// #include "RooRealVar.h"
// #include "RooDataHist.h"
// #include "RooHistPdf.h"
// #include "RooAddPdf.h"
// #include "RooFitResult.h"
// #include "RooPlot.h"
// #include "RooArgSet.h"
//how do I do the includes below without explicitly calling the directories?????
#ifndef __CINT__
#include "/swshare/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms18/include/RooGlobalFunc.h"
#endif
#include "/swshare/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms18/include/RooFit.h"
#include "/swshare/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms18/include/RooRealVar.h"
#include "/swshare/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms18/include/RooDataHist.h"
#include "/swshare/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms18/include/RooHistPdf.h"
#include "/swshare/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms18/include/RooAddPdf.h"
#include "/swshare/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms18/include/RooFitResult.h"
#include "/swshare/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms18/include/RooPlot.h"
#include "/swshare/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms18/include/RooArgSet.h"
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
#include "TMap.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
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

using namespace RooFit;
using namespace std;

//use via root -l -b -q ZnunuFromGammaEstimateFromRootfile.C++

void ZnunuFromGammaEstimateFromRootfile();
TH1D* RebinThisHistogram(TH1D* histogram, int nrebinnedbins, double* rebinnedbins);
TH1D* RescaleHisto(TH1D* histMC1, TH1D* histMC2, TH1D* hist_data_EB, TH1D* hist_data_EE);
TH1D* GetScaledHisto(TH1D* histo, float scale_fact, float scale_fact_err);
TH1D* RefillRatioHisto(TH1D* histo, float lastbinlowedge);
TH1D* RescaleUncertaintyHisto(TH1D* histo, float bulk_adduncertainty, float tail_adduncertainty, float separatingbin_lowedge);
float GetZnunuPreditionErrorStat(int i, TH1D* hData,  TH1D* ratio);
float GetZnunuPreditionErrorSysClosure(int i, TH1D* hPhotons, TH1D* ratio);
float GetZnunuPreditionErrorSys(int i, TH1D* hData, TH1D* hOther, TH1D* hQCD, TH1D* ratio);
float GetZnunuPreditionErrorSysRZG(int i, TH1D* hData, TH1D* hOther, TH1D* hQCD, TH1D* ratio);
void DoZnunuEstimate(map<string, TH1D*> histograms, map<string, TH1D*> histogramsZG, Bool_t PrintFullPrintOut, Bool_t PrintPredictionCard, Bool_t SavePrediction, Bool_t PrintSummaryTable, Bool_t PlotPrediction);
void PredictionCard(vector<int> srbin, vector<int> htbin, vector<int> mt2bin, vector<double> znunugen, vector<double> znunugenerr, vector<double> znunupred, vector<double> znunuprederrstat, vector<double> znunuprederrsyst, vector<double> ZGratio, vector<double> ZGratioerr, vector<double> modtabSF, vector<double> modtabSFerr);
void PredictionFile(vector<int> srbin, vector<int> htbin, vector<double> mt2low, vector<double> mt2up, vector<int> mt2last, vector<double> znunugen, vector<double> znunugenerr, vector<double> znunupred, vector<double> znunuprederrstat, vector<double> znunuprederrsyst, vector<double> znunuprederr, vector<double> znunuprederrsystRZG, vector<double> znunuprederrsyst1b0b, Bool_t PlotPrediction);
void SummaryTable(vector<int> srbin, vector<int> htbin, vector<double> mt2low, vector<double> mt2up, vector<int> mt2last, vector<double> znunugen, vector<double> znunugenerr, vector<double> ndata, vector<double> ngamma, vector<double> ngammaerr, vector<double> nqcd, vector<double> nqcderr, vector<double> nother, vector<double> nothererr, vector<double> ZGratio, vector<double> ZGratioerr, vector<double> modtabSF, vector<double> modtabSFerr, vector<double> znunupred, vector<double> znunuprederrstat, vector<double> znunuprederrsyst, vector<double> znunuprederr, vector<double> znunuprederrsystRZG, vector<double> znunuprederrsyst1b0b);

Bool_t  fSeparateEBEE              = true;  // set to false if EE and EB should not be separated in purity measurement of sigmaietaieta fit, default = true
Bool_t  fDoPhotonSigmaIEtaIEta     = true;  // set to "true" to perfom sigmaIetaIeta fit to obtain QCD normalization, default = true
Bool_t  fDoPhotonSignalRegion      = true; // get shapes for photons in signal region <-- needed for background prediction, default = true
Bool_t  fDoHadronicSignalRegion    = true; // get shapes for hadronic signal region <-- needed for background prediction, default = true
Bool_t  fDoPrediction              = true; // calculate the Znunu background from photon sample, needs fDoPhotonSignalRegion,fDoHadronicSignalRegion to be true, default = true
Bool_t  fPrintPredictionCard       = true; // gay printout - not needed anymore
Bool_t  fPrintFullPrintOut         = true;  // full printout during running
Bool_t  fMakeFinalTable            = true;  // make final prediction table
Bool_t  fDoVariableMT2bins         = true;  // variable MT2bins - not used
Bool_t  fDoPileUpWeights           = true;  // apply PU weights - default = true
Bool_t  fbSFReWeight               = true;  // apply BTV SF weights --> need to adjust fNBJets = {-10,[-3,3]}, default = true
Bool_t  fMrennaHack                = true;  // Steve Mrenna Status 2 parton jets - prompt photon dR - this is used to remove overlap of gamma MC and QCD MC
Bool_t  fUseConstantZToGammaR      = false; // Constant Z/gamma ratio: see analysis note!! - this should be always false
Bool_t  fUseConstZToGammaRdynamic  = true; // Constant Z/gamma ratio: see analysis note!! <-- from fMinMT2forConstR on, before Z/gamma is bin-by-bin - default = true
Float_t fMinMT2forConstR           = 350;   // min MT2 value after which a constant Z/gamma ratio is used, needs fUseConstZToGammaRdynamic = true, default = 350 for 8 TeV
// uncertainty
Bool_t  fAddFitIntrinsicUncert     = false; // if you'd like to add an intrinsic uncertainty to the template fit
Int_t   fVerbose                   = 9;
Float_t fFitIntrinsicUncert        = 0.0;
Bool_t  fAddRMCUncertainty         = true;  // systematic uncertainty on Z/gamma ratio
Float_t fRMCUncertainty            = 0.2;   // 0.3 for MT2 > 275, otherwise 0.2 (see AN)
Float_t fRMCUncertaintyTail        = 0.3;   // 0.3 for MT2 > 275, otherwise 0.2 (see AN)
// options ---------------
//TString fOutDir                    = "../Results/Filtered/GammaJetsPrediction/20130617_0b1b_test/";    // directory where shapes will be saved
TString fOutDir                    = "dummy/";    // directory where shapes will be saved
Bool_t  fSaveResults               = true;  // default = true //this one stores everything
Bool_t  fSaveResultsPlots          = true;  // default = true //this one is specific to final numbers of prediction and MC truth plotting
Bool_t  fSaveZnunuToGammaRatio     = true;  // default = true //a special flag
Bool_t  fDraw                      = true;  // draw histos or simply save them
Bool_t  fWriteToFile               = false; // writes couts to file, otherwise just prints it to terminal, default = false
Bool_t  fAppend                    = false; // writes couts to file, append to existing file
//the 4 numbers below are a relict from 7 TeV analysis, but are used if fUseConstantZToGammaR = true, which might be interesting for debugging isues
Float_t fConstantZToGammaR_LowHT   = 0.458; // 750 < HT < 950
Float_t fConstantZToGammaErr_LowHT = 0.057; // abs uncertainty on fConstantZToGammaR the relative uncertainty fRMCUncertainty is added in quadrature if fAddRMCUncertainty==true
Float_t fConstantZToGammaR_HighHT  = 0.628; // HT > 950
Float_t fConstantZToGammaErr_HighHT= 0.120; // abs uncertainty on fConstantZToGammaR 
// ------
Bool_t  fHT                        = true;   //run over HT dataset
Bool_t  fhighHT                    = true;   //medium or high HT.
Bool_t  fMET                       = false;  //run over single photon dataset
Bool_t  fISRreweight               = false;  //apply an effective ISR weight (effective means not reweight gen-photon/Z but observed distribution), default = false
Bool_t  fReScaleAfterMT2           = false;  //rescales MC histograms to fit data normalization - as histograms are with MT2 cut, this means renormalizating after applying MT2 cut, default = false
Bool_t  fdontRescaleRatio          = true;   //keep this true, for Z/gamma calculation use photon MC not normalized to data in sigmaietaieta fit - idea is that both the gamma MC and Z MC should ne on equal footing
Float_t fGammakFactor              = 1.0;//this is the kFaktor to gamma MC from the samples.dat, is needed to have correct Z/G ratio

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


//definition of topological regions
const int signalregionsize = 12;
string signal_region[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b", "2j1to2bmod", "3to5j1bmod", "6j1bmod"};
//definition of HT regions, the names are set by the fMET/fHT flags above
const int HTbinsize = 2;
string HT_bin[HTbinsize] = {"HTge450", "HTge750"};//dummy
const int sampletypesize = 6;
string sample_type[sampletypesize] = {"QCD", "Photon", "Znunu", "Other", "Susy", "Data"};

std::ostringstream* fLogStream     = 0;

//see ZnunuFromGammaEstimate.C for the comments
void ZnunuFromGammaEstimateFromRootfile(){

	gSystem->Load("libPhysics");
	gSystem->Load("libRooFit") ;
	fLogStream = new std::ostringstream();

	if(fMET) fDoPhotonSigmaIEtaIEta = false;
	if(!fDoPhotonSigmaIEtaIEta) cout << "YOUR ARE NOT DOING fDoPhotonSigmaIEtaIEta " << endl;
	//adjust k factor
	if(fMET)      fGammakFactor = 0.9;
	else if (fHT) fGammakFactor = 0.8;
	cout << "photonMC k factor is " << fGammakFactor << endl;

	if(fMET){ HT_bin[0] = "HTge450"; }
	if(fHT) { HT_bin[0] = "HTge750"; HT_bin[1] = "HTge1200"; }
	if(fMET==true) fHT = false;				//either fHT or fMET has to be false
	if(fMET==false && fHT==false) fHT = true;		//either fHT or fMET has to be true

	TString filedirectory;
	if(fISRreweight) fOutDir = fOutDir + "ISR/";
	if(fMET) fOutDir = fOutDir + "MET/";
	if(fHT)  fOutDir = fOutDir + "HT/";
	filedirectory = fOutDir;
//	if(fHT){
//		if(fhighHT) filedirectory = filedirectory + "high/";
//		else        filedirectory = filedirectory + "medium/";
//	}
    	Util::MakeOutputDir(fOutDir);
	TString outputdir = fOutDir;

	TString inputname = "ZnunuGammaHistograms.root";
	if(!fISRreweight) inputname = "NoISR_" + inputname;

	map<string, TH1D*>    histos;
	map<string, TH1D*>    histosZG;
	vector<string> histonames; histonames.clear();
	//get all histogram names you need
	histonames.push_back("SigmaIEtaIEta_EB");			//photon inclusive, sigmaietaieta distribution used for fitting, EB
	histonames.push_back("SigmaIEtaIEta_EE");			//photon inclusive, sigmaietaieta distribution used for fitting, EE
	histonames.push_back("PhotonicSignalRegion_EB");		//photon signal, MT2 distribution used for Estimate, EB
	histonames.push_back("PhotonicSignalRegion_EE");		//photon signal, MT2 distribution used for Estimate, EE
	histonames.push_back("HadronicRegion");				//hadronic signal MC truth, MT2 distribution you want to estimate

	TString fileinputname = filedirectory+inputname;
	TFile *inputfile = TFile::Open(fileinputname.Data());
	cout << fileinputname << endl;

	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTbinsize;        ++i3){
		if(fMET && i3==1) continue;
//		if(fHT &&  fhighHT && i3==0) continue;
//		if(fHT && !fhighHT && i3==1) continue;
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
			if(signal_region[i2]=="2j1to2bmod") NMT2bins = gNMT2bins_2j0b_lHT;
			if(signal_region[i2]=="3to5j1bmod") NMT2bins = gNMT2bins_3j0b_lHT;
			if(signal_region[i2]=="6j1bmod")    NMT2bins = gNMT2bins_6j0b_lHT;
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
			if(signal_region[i2]=="2j1to2bmod") NMT2bins = gNMT2bins_2j0b_mHT;
			if(signal_region[i2]=="3to5j1bmod") NMT2bins = gNMT2bins_3j0b_mHT;
			if(signal_region[i2]=="6j1bmod")    NMT2bins = gNMT2bins_6j0b_mHT;
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
			if(signal_region[i2]=="2j1to2bmod") NMT2bins = gNMT2bins_2j0b_hHT;
			if(signal_region[i2]=="3to5j1bmod") NMT2bins = gNMT2bins_3j0b_hHT;
			if(signal_region[i2]=="6j1bmod")    NMT2bins = gNMT2bins_6j0b_hHT;
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
			if(signal_region[i2]=="2j1to2bmod") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_2j0b_lHT[i0]; }
			if(signal_region[i2]=="3to5j1bmod") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j0b_lHT[i0]; }
			if(signal_region[i2]=="6j1bmod")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j0b_lHT[i0]; }

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
			if(signal_region[i2]=="2j1to2bmod") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_2j0b_mHT[i0]; }
			if(signal_region[i2]=="3to5j1bmod") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j0b_mHT[i0]; }
			if(signal_region[i2]=="6j1bmod")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j0b_mHT[i0]; }
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
			if(signal_region[i2]=="2j1to2bmod") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_2j0b_hHT[i0]; }
			if(signal_region[i2]=="3to5j1bmod") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j0b_hHT[i0]; }
			if(signal_region[i2]=="6j1bmod")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j0b_hHT[i0]; }
		   }
		}
		string hs = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_") + sample_type[i1];
		string mapname;
		for(unsigned int i0 = 0; i0<histonames.size(); ++i0){
			if(histonames[i0]=="SigmaIEtaIEta_EB" || histonames[i0]=="SigmaIEtaIEta_EE"){
				mapname = histonames[i0] + hs;
				if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)inputfile->Get(mapname.c_str());
			} else {
				mapname = histonames[i0] + hs;
				if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)inputfile->Get(mapname.c_str());
			}
		}
		if(i1==0) {//not sampletype dependent
			mapname = "ZGRatio"+string("_") + HT_bin[i3] + string("_") + signal_region[i2];
			if(histosZG.count(mapname) == 0 ){ histosZG[mapname] = new TH1D(mapname.c_str(), "", NMT2bins, MT2bins); //cout << mapname << " " << histos.count(mapname) << endl;
			                                   histosZG[mapname]->GetXaxis()->SetTitle("M_{T2}");
			                                   histosZG[mapname]->GetYaxis()->SetTitle("Z(#nu#bar{#nu})/#gamma"); }
		}
	}}}
//	inputfile->Close();
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h) { cout << h->first << " "; cout << h->second->Integral() << endl; }

	//do estimation
	cout << "do estimation" << endl;
	DoZnunuEstimate(histos, histosZG, fPrintFullPrintOut, fPrintPredictionCard, fSaveResults, fMakeFinalTable, fSaveResultsPlots);

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

}//ZnunuFromGammaEstimate()

//this function is actually doing the estimate
void DoZnunuEstimate(map<string, TH1D*> histograms, map<string, TH1D*> histogramsZG, Bool_t PrintFullPrintOut, Bool_t PrintPredictionCard, Bool_t SavePrediction, Bool_t PrintSummaryTable, Bool_t PlotPrediction){

	map<string, TH1D*> histos   = histograms;
	map<string, TH1D*> histosZG = histogramsZG;//note that they are empty and get stored in this function

	string helpstring = "";
	double helpdouble = -1.;

	//note: keep in mind you want to do a table in the end
	vector<int>     htbin;
	vector<int>     mt2bin;
	vector<int>     srbin;
	vector<double>  znunugen;
	vector<double>  znunugenerr;
	vector<double>  znunupred;
	vector<double>  znunuprederr;
	vector<double>  znunuprederrstat;
	vector<double>  znunuprederrsyst;
	vector<double>  znunuprederrsystRZG;
	vector<double>  znunuprederrsyst1b0b;
	vector<double>  ZGratio;
	vector<double>  ZGratioerr;
	vector<double>  mt2low;
	vector<double>  mt2up;
	vector<int>     mt2last;
	vector<double>  numdata;
	vector<double>  numdataerr;
	vector<double>  numqcd;
	vector<double>  numqcderr;
	vector<double>  numother;
	vector<double>  numothererr;
	vector<double>  numgamma;
	vector<double>  numgammaerr;
	vector<double>  modtabSF;
	vector<double>  modtabSFerr;
	htbin.clear(); mt2bin.clear(); mt2low.clear(); mt2up.clear(); mt2last.clear(); srbin.clear();
	znunugen.clear(); znunugenerr.clear(); 
	znunupred.clear(); znunuprederr.clear(); znunuprederrstat.clear(); znunuprederrsyst.clear(); znunuprederrsyst1b0b.clear(); znunuprederrsystRZG.clear();
	ZGratio.clear(); ZGratioerr.clear(); modtabSF.clear(); modtabSFerr.clear();
	numdata.clear(); numdataerr.clear(); numqcd.clear(); numqcderr.clear(); numother.clear(); numothererr.clear(); numgamma.clear(); numgammaerr.clear();

	for(int i3 = 0; i3<HTbinsize;        ++i3){
	for(int i2 = 0; i2<signalregionsize; ++i2){
//sigmaietaieta code crashes for highHT, topological region >= 3to5j0b - WHY??????
		if(fMET && i3==1) continue;
//		if( fhighHT && i3==0) continue;
//		if(!fhighHT && i3==1) continue;
		//as all histograms have the same binning just take first histogram 
		string hshl = string("_") + HT_bin[i3] + string("_") + signal_region[i2];

		//first do sigmaietaieta - if flag fSeparateEBEE = false --> all stored in EB and EE is copy of EE
		//ML means sigmaIetaIeta Maximum Likelihood fit
		float MLPhotonEB_SF(-1), MLPhotonEB_SFerr(-1), MLQCDEB_SF(-1), MLQCDEB_SFerr(-1), MLOtherEB_SF(-1), MLOtherEB_SFerr(-1);
		float MLPhotonEE_SF(-1), MLPhotonEE_SFerr(-1), MLQCDEE_SF(-1), MLQCDEE_SFerr(-1), MLOtherEE_SF(-1), MLOtherEE_SFerr(-1);
		float MLPhotonEB(-1), MLPhotonEBerr(-1), MLQCDEB(-1), MLQCDEBerr(-1), MLOtherEB(-1), MLOtherEBerr(-1), MLDataEB(-1), MLDataEBerr(-1);
		float MLPhotonEE(-1), MLPhotonEEerr(-1), MLQCDEE(-1), MLQCDEEerr(-1), MLOtherEE(-1), MLOtherEEerr(-1), MLDataEE(-1), MLDataEEerr(-1);
		if(fDoPhotonSigmaIEtaIEta){
			string h = "SigmaIEtaIEta_EB" + hshl + "_QCD";    TH1D *hQCD       = (TH1D*)histos[h]->Clone("QCD");
			       h = "SigmaIEtaIEta_EB" + hshl + "_Photon"; TH1D *hPhotons   = (TH1D*)histos[h]->Clone("Photons");
			       h = "SigmaIEtaIEta_EB" + hshl + "_Data";   TH1D *hData      = (TH1D*)histos[h]->Clone("Data");
			       h = "SigmaIEtaIEta_EB" + hshl + "_Other";  TH1D *hOther     = (TH1D*)histos[h]->Clone("Other");
			       h = "SigmaIEtaIEta_EE" + hshl + "_QCD";    TH1D *hQCDEE     = (TH1D*)histos[h]->Clone("QCDEE");
			       h = "SigmaIEtaIEta_EE" + hshl + "_Photon"; TH1D *hPhotonsEE = (TH1D*)histos[h]->Clone("PhotonsEE");
			       h = "SigmaIEtaIEta_EE" + hshl + "_Data";   TH1D *hDataEE    = (TH1D*)histos[h]->Clone("DataEE");
			       h = "SigmaIEtaIEta_EE" + hshl + "_Other";  TH1D *hOtherEE   = (TH1D*)histos[h]->Clone("OtherEE");

			if(!fSeparateEBEE){
				hQCD    ->Add(hQCDEE);
				hPhotons->Add(hPhotonsEE);
				hData   ->Add(hDataEE);
				hOther  ->Add(hOtherEE);
			}
			// Add Pedestal to QCD MC PDF in order to avoid bins in PDF with zero entries
			// but data in same bin! this causes problems!
			for(int i=0; i<hQCD->GetNbinsX(); ++i){
				if(hData->GetBinContent(i)==0) continue;
				double MCcontent = hPhotons->GetBinContent(i)+hQCD->GetBinContent(i)+hOther->GetBinContent(i);
				if(MCcontent==0) {hQCD->SetBinContent(i, 1E-01);hQCD->SetBinError(i, 1E-01);}
			}
			for(int i=0; i<hQCDEE->GetNbinsX(); ++i){
				if(hDataEE->GetBinContent(i)==0) continue;
				double MCcontent = hPhotonsEE->GetBinContent(i)+hQCDEE->GetBinContent(i)+hOtherEE->GetBinContent(i);
				if(MCcontent==0) {hQCDEE->SetBinContent(i, 1E-01);hQCDEE->SetBinError(i, 1E-01);}
			}
			if(hOther->Integral()<=0){
				hOther->Add(hQCD);
				if(hOther->Integral()>1.) hOther->Scale(0.01/hOther->Integral());//total hOther yield is <= 0.01;
				else                      hOther->Scale(0.01);
			}
			if(hOtherEE->Integral()<=0){
				hOtherEE->Add(hQCDEE);
				if(hOtherEE->Integral()>1.) hOther->Scale(0.01/hOtherEE->Integral());//total hOther yield is <= 0.01;
				else                        hOther->Scale(0.01);
			}
			RooRealVar sigmaietaieta("sigmaietaieta","sigmaietaieta",0.,0.08) ; // contained in histos
		
			RooDataHist Data   ("data"   ,"data"  ,sigmaietaieta,hData) ;    // define RooDataHists
			RooDataHist Photons("photons","photon",sigmaietaieta,hPhotons);
			RooDataHist QCD    ("QCD"    ,"QCD"   ,sigmaietaieta,hQCD);
			RooDataHist Other  ("Other"  ,"Other" ,sigmaietaieta,hOther);//check this!!!
		
			RooHistPdf Photons_pdf("photons_pdf","photons_pdf",sigmaietaieta,Photons); // define PDFs for signal and bkg
			RooHistPdf QCD_pdf    ("qcd_pdf"    ,"qcd_pdf"    ,sigmaietaieta,QCD    ); 
			RooHistPdf Other_pdf  ("Other_pdf"  ,"other_pdf"  ,sigmaietaieta,Other  );
		
			RooRealVar nsig       ("nsig"   ,"number of signal events",     hPhotons->Integral()  ,  hPhotons->Integral()*0.1,hData->Integral());
			RooRealVar nqcd       ("nqcd"   ,"number of QCD events",        hQCD->Integral()      ,  hQCD->Integral()    *0.5,hData->Integral());
			RooRealVar nother     ("nother" ,"number of Other SM events",   hOther->Integral()); nother.setConstant(kTRUE);

			// model(x) = nsig*Photons_pdf(x) + nqcd*QCD_pdf(x) + nother*Other_pdf(x), where nother is fixed to nominal contribution
			RooAddPdf model("model","model", RooArgList(Photons_pdf,QCD_pdf,Other_pdf), RooArgList(nsig, nqcd, nother));
			model.defaultPrintStream(fLogStream);
cout << __LINE__ << " " << hshl << endl;
		
//the code crashes here for high HT, topological region >=3to5j0b
			// perform fit
			RooFitResult* fitres = model.fitTo(Data,SumW2Error(kFALSE),Extended(), Save(kTRUE)); 
			// if I'm not mistaken: SumW2==false is the right option, as mc-histos already contain proper weights. 
			// SumW2Error == true would be needed if input comes from a TTree with (observable, weight) for each entry. 
			// then, data.setWeightVar(y) would also be needed. 
cout << __LINE__ << " " << hshl << endl;

			TH1D* dataclone = (TH1D*)hData->Clone("dataclone");
			dataclone->Rebin(dataclone->GetNbinsX());
			MLPhotonEB       = nsig.getVal();
			MLPhotonEBerr    = nsig.getError();
			MLPhotonEB_SF    = nsig.getVal()/hPhotons->Integral();
			MLPhotonEB_SFerr = nsig.getError()/hPhotons->Integral();
			MLQCDEB          = nqcd.getVal();
			MLQCDEBerr       = nqcd.getError();
			MLQCDEB_SF       = nqcd.getVal()/hQCD->Integral();
			MLQCDEB_SFerr    = nqcd.getError()/hQCD->Integral();
			MLOtherEB        = nother.getVal();
			MLOtherEBerr     = nother.getError();
			MLOtherEB_SF     = nother.getVal()/hOther->Integral();
			MLOtherEB_SFerr  = nother.getError()/hOther->Integral();
			MLDataEB         = dataclone->GetBinContent(1);
			MLDataEBerr      = dataclone->GetBinError(1);
			delete dataclone;
			if(fAddFitIntrinsicUncert){ // adding uncertainty to scale factors due to intrinsic method uncertainty
				*fLogStream << "+++ Adding a uncertainty of " << fFitIntrinsicUncert << "  percent fit Scale Factors" << endl;
				MLPhotonEB_SFerr = sqrt(pow(MLPhotonEB_SFerr,2)+pow(fFitIntrinsicUncert*MLPhotonEB_SF,2));
				MLQCDEB_SFerr    = sqrt(pow(MLQCDEB_SFerr   ,2)+pow(fFitIntrinsicUncert*MLQCDEB_SF   ,2));
			}
			if(fSaveResults){
				TString canvname = "sigmaietaietaEB"+hshl;
				TCanvas* canv = new TCanvas(canvname.Data(),"", 0, 0, 500, 500 );
				RooPlot* frame = sigmaietaieta.frame();
				Data.plotOn(frame, Name("Data")) ;
				model.plotOn(frame,Components(RooArgSet(Photons_pdf,QCD_pdf,Other_pdf)), Name("Model"));
				model.plotOn(frame,Components(QCD_pdf),        LineStyle(kDotted), LineColor(kMagenta));
				model.plotOn(frame,Components(Photons_pdf),    LineStyle(kDotted), LineColor(kGreen));
				frame->Draw();
				Double_t chi2 = frame->chiSquare("Model", "Data", 3);
				if(fPrintFullPrintOut){
					*fLogStream << "-----------------------------------------------------------------" << endl;
					*fLogStream << "Fit result for: " <<canvname                                       << endl; 
					fitres->Print("v");
					*fLogStream << "ChiSquare of fit: " << chi2                                        << endl;
					*fLogStream << "-----------------------------------------------------------------" << endl;
				}
				// save RooFit output:
				TString outdir = fOutDir+"sigmaietaietafits/";
				Util::MakeOutputDir(outdir);
				TString filename=outdir+"FitEBEE_"+signal_region[i2]+"_"+HT_bin[i3]+".root";
				if(fSeparateEBEE) filename=outdir+"FitEB_"+signal_region[i2]+"_"+HT_bin[i3]+".root";
				TFile *file = new TFile(filename.Data(), "UPDATE");
				fitres->Write();
				frame->Write();
				canv->Write();
			//	Data.Write();
			//	Photons.Write();
			//	QCD.Write();
			//	Photons_pdf.Write();
			//	QCD_pdf.Write();
			//	Other_pdf.Write();
			//	nsig.Write();
			//	nqcd.Write();
			//	nother.Write();
			//	model.Write();
				file->Close();
				delete file;
				canv->Close();
				delete canv;
			}
		//	delete fitres;
			if(fSeparateEBEE){//separate for EE
				RooRealVar sigmaietaietaEE("sigmaietaietaEE","sigmaietaieta",0.,0.08) ; // contained in histos
			
				RooDataHist DataEE   ("dataEE"   ,"data"  ,sigmaietaietaEE,hDataEE) ;    // define RooDataHists
				RooDataHist PhotonsEE("photonsEE","photon",sigmaietaietaEE,hPhotonsEE);
				RooDataHist QCDEE    ("QCDEE"    ,"QCD"   ,sigmaietaietaEE,hQCDEE);
				RooDataHist OtherEE  ("OtherEE"  ,"Other" ,sigmaietaietaEE,hOtherEE);//check this!!!
			
				RooHistPdf PhotonsEE_pdf("photonsEE_pdf","photons_pdf",sigmaietaietaEE,PhotonsEE); // define PDFs for signal and bkg
				RooHistPdf QCDEE_pdf    ("qcdEE_pdf"    ,"qcd_pdf"    ,sigmaietaietaEE,QCDEE    ); 
				RooHistPdf OtherEE_pdf  ("OtherEE_pdf"  ,"other_pdf"  ,sigmaietaietaEE,OtherEE  );
			
				RooRealVar nsigEE       ("nsigEE"   ,"number of signal events",     hPhotonsEE->Integral()  ,  hPhotonsEE->Integral()*0.1,hDataEE->Integral());
				RooRealVar nqcdEE       ("nqcdEE"   ,"number of QCD events",        hQCDEE->Integral()      ,  hQCDEE->Integral()    *0.5,hDataEE->Integral());
				RooRealVar notherEE     ("notherEE" ,"number of Other SM events",   hOtherEE->Integral()); notherEE.setConstant(kTRUE);
			
				// model(x) = nsig*Photons_pdf(x) + nqcd*QCD_pdf(x) + nother*Other_pdf(x), where nother is fixed to nominal contribution
				RooAddPdf modelEE("modelEE","model", RooArgList(PhotonsEE_pdf,QCDEE_pdf,OtherEE_pdf), RooArgList(nsigEE, nqcdEE, notherEE));
				modelEE.defaultPrintStream(fLogStream);
cout << __LINE__ << " " << hshl << endl;
				
				// perform fit
				RooFitResult* fitresEE = modelEE.fitTo(DataEE,SumW2Error(kFALSE),Extended(), Save(kTRUE)); 
				// if I'm not mistaken: SumW2==false is the right option, as mc-histos already contain proper weights. 
				// SumW2Error == true would be needed if input comes from a TTree with (observable, weight) for each entry. 
				// then, data.setWeightVar(y) would also be needed. 
cout << __LINE__ << " " << hshl << endl;

				TH1D* datacloneEE = (TH1D*)hDataEE->Clone("datacloneEE");
				datacloneEE->Rebin(datacloneEE->GetNbinsX());
				MLPhotonEE       = nsigEE.getVal();
				MLPhotonEEerr    = nsigEE.getError();
				MLPhotonEE_SF    = nsigEE.getVal()/hPhotonsEE->Integral();
				MLPhotonEE_SFerr = nsigEE.getError()/hPhotonsEE->Integral();
				MLQCDEE          = nqcdEE.getVal();
				MLQCDEEerr       = nqcdEE.getError();
				MLQCDEE_SF       = nqcdEE.getVal()/hQCDEE->Integral();
				MLQCDEE_SFerr    = nqcdEE.getError()/hQCDEE->Integral();
				MLOtherEE        = notherEE.getVal();
				MLOtherEEerr     = notherEE.getError();
				MLOtherEE_SF     = notherEE.getVal()/hOtherEE->Integral();
				MLOtherEE_SFerr  = notherEE.getError()/hOtherEE->Integral();
				MLDataEE         = datacloneEE->GetBinContent(1);
				MLDataEEerr      = datacloneEE->GetBinError(1);
				delete datacloneEE;

				if(fAddFitIntrinsicUncert){ // adding uncertainty to scale factors due to intrinsic method uncertainty
					*fLogStream << "+++ Adding a uncertainty of " << fFitIntrinsicUncert << "  percent fit Scale Factors" << endl;
					MLPhotonEE_SFerr = sqrt(pow(MLPhotonEE_SFerr,2)+pow(fFitIntrinsicUncert*MLPhotonEE_SF,2));
					MLQCDEE_SFerr    = sqrt(pow(MLQCDEE_SFerr   ,2)+pow(fFitIntrinsicUncert*MLQCDEE_SF   ,2));
				}
				if(fSaveResults){
					TString canvname = "sigmaietaietaEE"+hshl;
					TCanvas* canv = new TCanvas(canvname.Data(),"", 0, 0, 500, 500 );
					RooPlot* frame = sigmaietaietaEE.frame();
					DataEE.plotOn(frame, Name("Data")) ;
					modelEE.plotOn(frame,Components(RooArgSet(PhotonsEE_pdf,QCDEE_pdf,OtherEE_pdf)), Name("Model"));
					modelEE.plotOn(frame,Components(QCDEE_pdf),        LineStyle(kDotted), LineColor(kMagenta));
					modelEE.plotOn(frame,Components(PhotonsEE_pdf),    LineStyle(kDotted), LineColor(kGreen));
					frame->Draw();
					Double_t chi2 = frame->chiSquare("Model", "Data", 3);
					if(fPrintFullPrintOut){
						*fLogStream << "-----------------------------------------------------------------" << endl;
						*fLogStream << "Fit result for: " << canvname                                      << endl; 
						fitresEE->Print("v");
						*fLogStream << "ChiSquare of fit: " << chi2                                        << endl;
						*fLogStream << "-----------------------------------------------------------------" << endl;
					}
					// save RooFit output:
					TString outdir = fOutDir+"sigmaietaietafits/";
					Util::MakeOutputDir(outdir);
					TString filename=outdir+"FitEE_"+signal_region[i2]+"_"+HT_bin[i3]+".root";
					TFile *file = new TFile(filename.Data(), "UPDATE");
					fitresEE->Write();
					frame->Write();
					canv->Write();
				//	DataEE.Write();
				//	PhotonsEE.Write();
				//	QCDEE.Write();
				//	PhotonsEE_pdf.Write();
				//	QCDEE_pdf.Write();
				//	OtherEE_pdf.Write();
				//	nsigEE.Write();
				//	nqcdEE.Write();
				//	notherEE.Write();
				//	modelEE.Write();
					file->Close();
					delete file;
					canv->Close();
					delete canv;
				}
				//delete fitresEE;
			}//if(fSeparateEBEE)
			else {
				MLPhotonEE_SF =  MLPhotonEB_SF;
				MLQCDEE_SF    =  MLQCDEB_SF;
				MLOtherEE_SF  =  MLOtherEB_SF;
				MLPhotonEE_SFerr = MLPhotonEB_SFerr;
				MLQCDEE_SFerr    =  MLQCDEB_SFerr;
				MLOtherEE_SFerr  =  MLOtherEB_SFerr;
			}
			//delete fitres;
			delete hQCD;
			delete hPhotons;
			delete hData;
			delete hOther;
			delete hQCDEE;
			delete hPhotonsEE;
			delete hDataEE;
			delete hOtherEE;
		} else {
			MLPhotonEB_SF =  1.2; MLPhotonEE_SF =  1.2;
			MLQCDEB_SF    =  1.3; MLQCDEE_SF    =  1.3;
			MLOtherEB_SF  =  1.0; MLOtherEE_SF  =  1.0;
			MLPhotonEB_SFerr = 0.05; MLPhotonEE_SFerr = 0.05;
			MLQCDEB_SFerr    =  1.3; MLQCDEE_SFerr    =  1.3;
			MLOtherEB_SFerr  =  0.0; MLOtherEE_SFerr  =  0.0;
		}//else w.r.t. if(fDoPhotonSigmaIEtaIEta)
		//sigmaietaieta fit is done.

		//now you need to obtain the ZGamma ratio
		if(fUseConstantZToGammaR){
			for(int i = 1; i<=histosZG["ZGRatio"+hshl]->GetNbinsX(); ++i){
				if(fHT){
					histosZG["ZGRatio"+hshl]->SetBinContent(i,fConstantZToGammaR_HighHT);
					histosZG["ZGRatio"+hshl]->SetBinError(  i,fConstantZToGammaErr_HighHT);
				} if(fMET){
					histosZG["ZGRatio"+hshl]->SetBinContent(i,fConstantZToGammaR_LowHT);
					histosZG["ZGRatio"+hshl]->SetBinError(  i,fConstantZToGammaErr_LowHT);
				}
			}
		} else {
			if(histos["PhotonicSignalRegion_EB"+hshl+"_Photon"]==0||histos["PhotonicSignalRegion_EE"+hshl+"_Photon"]==0||histos["HadronicRegion"+hshl+"_Znunu"]==0) { 
				cout << "GetMCZnunuToPhotonRatio: received 0 pointer!" << endl;
				exit(-1);
			}
			double rescaleGammaEB = MLPhotonEB_SF; double rescaleGammaEBerr = MLPhotonEB_SFerr;
			double rescaleGammaEE = MLPhotonEE_SF; double rescaleGammaEEerr = MLPhotonEE_SFerr;
			double rescaleQCDEB   = MLQCDEB_SF;    double rescaleQCDEBerr   = MLQCDEB_SFerr;
			double rescaleQCDEE   = MLQCDEE_SF;    double rescaleQCDEEerr   = MLQCDEE_SFerr;
			if(fdontRescaleRatio){
				rescaleGammaEB = 1./fGammakFactor; rescaleGammaEBerr = 0.;
				rescaleGammaEE = 1./fGammakFactor; rescaleGammaEEerr = 0.;
				rescaleQCDEB   = 1.; rescaleQCDEBerr   = 0.;
				rescaleQCDEE   = 1.; rescaleQCDEEerr   = 0.;
			}
			//the GetScaledHisto2 keeps the binning and scales the histogram with a factor (which here is obtained from the sigmaietaieta fit)
			TH1D *currPhotons= GetScaledHisto(histos["PhotonicSignalRegion_EB"+hshl+"_Photon"] , rescaleGammaEB, rescaleGammaEBerr); // EB
			TH1D *currPhotonsEE= GetScaledHisto(histos["PhotonicSignalRegion_EE"+hshl+"_Photon"] , rescaleGammaEE, rescaleGammaEEerr); // EE
		//	for(int nnn = 1; nnn<=currPhotons->GetNbinsX();++nnn){
		//		*fLogStream << "bin " << currPhotons->GetBinLowEdge(nnn) << "-" << currPhotons->GetBinLowEdge(nnn)+currPhotons->GetBinWidth(nnn) << ": ";
		//		*fLogStream << "EB " << currPhotons->GetBinContent(nnn) << ", EE " << currPhotonsEE->GetBinContent(nnn) << endl;
		//	}
			currPhotons->Add(  currPhotonsEE); //EE
			if(fReScaleAfterMT2){
				TH1D *currQCD= GetScaledHisto(histos["PhotonicSignalRegion_EB"+hshl+"_QCD"] , rescaleQCDEB, rescaleQCDEBerr); // EB
				currQCD->Add(  GetScaledHisto(histos["PhotonicSignalRegion_EE"+hshl+"_QCD"] , rescaleQCDEE, rescaleQCDEEerr)); //EE
				currPhotons = RescaleHisto(currPhotons, currQCD, histos["PhotonicSignalRegion_EB"+hshl+"_Data"], histos["PhotonicSignalRegion_EE"+hshl+"_Data"]);
				delete currQCD;
			}
			TH1D *currZnunu  = GetScaledHisto(histos["HadronicRegion"+hshl+"_Znunu"], 1.                   , 0.);
		//	for(int nnn = 1; nnn<=currZnunu->GetNbinsX();++nnn){
		//		*fLogStream << "bin " << currZnunu->GetBinLowEdge(nnn) << "-" << currZnunu->GetBinLowEdge(nnn)+currZnunu->GetBinWidth(nnn) << ": ";
		//		*fLogStream << "Znunu " << currZnunu->GetBinContent(nnn) << ", Gamma " << currPhotons->GetBinContent(nnn) << " ==>> ";
		//		*fLogStream << "ratio " << currZnunu->GetBinContent(nnn)/currPhotons->GetBinContent(nnn) << endl;
		//	}
			//constant Z/G ratio for high MT2 is achieved by setting the bins of the input histograms with fMinMT2forConstR to the sum of those bins
			if(fUseConstZToGammaRdynamic) currPhotons = RefillRatioHisto(currPhotons,fMinMT2forConstR);
			if(fUseConstZToGammaRdynamic) currZnunu   = RefillRatioHisto(currZnunu  ,fMinMT2forConstR);
	
			//get the ratio
			helpstring = "ZGRatio"+hshl;
			TH1D *ratio = (TH1D*)currZnunu->Clone("Znunu_To_Photon_ratio");
			ratio->Divide(currPhotons);
			if(fAddRMCUncertainty){
				if(fPrintFullPrintOut){
					*fLogStream << "------------------------------------------------------------------------------" << endl;
					*fLogStream << "+++ adding in quadrature " << fRMCUncertainty-fRMCUncertaintyTail << " percent uncertainty on R   " << endl;
					*fLogStream << "------------------------------------------------------------------------------" << endl;
				}
				ratio = RescaleUncertaintyHisto(ratio, fRMCUncertainty, fRMCUncertaintyTail, fMinMT2forConstR);
			}
			histosZG["ZGRatio"+hshl] = (TH1D*)ratio->Clone(helpstring.c_str());

			if(fSaveZnunuToGammaRatio){
				TString outdir = fOutDir+"ZnunuToGammaRatio/";
				Util::MakeOutputDir(outdir);
				TString filename=outdir+ "Ratio_"+signal_region[i2]+"_"+HT_bin[i3]+".root";
				TFile *file = new TFile(filename.Data(), "RECREATE");
				currZnunu->Write();
				currPhotons->Write();
				ratio->Write();
				histosZG["ZGRatio"+hshl]->Write();
				file->Close();
				delete file;
			}
			delete currPhotons;
			delete currZnunu;
			delete ratio;
		}
		//ZGammaRatio is done.

		//now do the prediction
		if(!fDoPrediction) continue;
		if(fPrintFullPrintOut){
			*fLogStream << "******************************* Prediction ****************************************" << endl;
			*fLogStream << "Photonic Region where EB Normalization was extracted: (ML-fit result) ------------" << endl;
			*fLogStream << "  NData: "    << MLDataEB    << " pm " << MLDataEBerr       << endl;
			*fLogStream << "  NPhotons: " << MLPhotonEB  << " pm " << MLPhotonEBerr     << endl;
			*fLogStream << "  NOther: "   << MLOtherEB   << " pm " << MLOtherEBerr      << endl;
			*fLogStream << "  NQCD: "     << MLQCDEB     << " pm " << MLQCDEBerr        << endl;
			*fLogStream << "Photonic Region where EE Normalization was extracted: (ML-fit result) ------------" << endl;
			*fLogStream << "  NData: "    << MLDataEE    << " pm " << MLDataEEerr       << endl;
			*fLogStream << "  NPhotons: " << MLPhotonEE  << " pm " << MLPhotonEEerr     << endl;
			*fLogStream << "  NOther: "   << MLOtherEE   << " pm " << MLOtherEEerr      << endl;
			*fLogStream << "  NQCD: "     << MLQCDEE     << " pm " << MLQCDEEerr        << endl;
			*fLogStream << "Scale factors are  ----------------------------------------------------------------" << endl;
			*fLogStream << "EB ---------------------------------------------------------------------------------" << endl;
			*fLogStream << "  NPhotons: " << MLPhotonEB_SF  << " pm " << MLPhotonEB_SFerr     << endl;
			*fLogStream << "  NQCD: "     << MLQCDEB_SF     << " pm " << MLQCDEB_SFerr        << endl;
			*fLogStream << "EE ---------------------------------------------------------------------------------" << endl;
			*fLogStream << "  NPhotons: " << MLPhotonEE_SF  << " pm " << MLPhotonEE_SFerr     << endl;
			*fLogStream << "  NQCD: "     << MLQCDEE_SF     << " pm " << MLQCDEE_SFerr        << endl;
		}
		//the GetScaledHisto keeps the binning and scales the histogram with a factor (which here is obtained from the sigmaietaieta fit)
		// EB
		TH1D *hPhotons_EB = GetScaledHisto(histos["PhotonicSignalRegion_EB"+hshl+"_Photon"], MLPhotonEB_SF, MLPhotonEB_SFerr);
		TH1D *hOther_EB   = GetScaledHisto(histos["PhotonicSignalRegion_EB"+hshl+"_Other" ], MLOtherEB_SF,  MLOtherEB_SFerr);
		TH1D *hQCD_EB     = GetScaledHisto(histos["PhotonicSignalRegion_EB"+hshl+"_QCD"   ], MLQCDEB_SF,    MLQCDEB_SFerr);
		TH1D *hData_EB    = GetScaledHisto(histos["PhotonicSignalRegion_EB"+hshl+"_Data"  ], 1.                   , 0.);
		// EE
		TH1D *hPhotons_EE = GetScaledHisto(histos["PhotonicSignalRegion_EE"+hshl+"_Photon"], MLPhotonEE_SF, MLPhotonEE_SFerr);
		TH1D *hOther_EE   = GetScaledHisto(histos["PhotonicSignalRegion_EE"+hshl+"_Other" ], MLOtherEE_SF,  MLOtherEE_SFerr);
		TH1D *hQCD_EE     = GetScaledHisto(histos["PhotonicSignalRegion_EE"+hshl+"_QCD"   ], MLQCDEE_SF,    MLQCDEE_SFerr);
		TH1D *hData_EE    = GetScaledHisto(histos["PhotonicSignalRegion_EE"+hshl+"_Data"  ], 1.                   , 0.);
		
		TH1D* hPhotons = (TH1D*)hPhotons_EB->Clone("hPhotons"); hPhotons->Add(hPhotons_EE);
		TH1D* hOther   = (TH1D*)hOther_EB->Clone("hOther");     hOther  ->Add(hOther_EE);
		TH1D* hQCD     = (TH1D*)hQCD_EB->Clone("hQCD");         hQCD    ->Add(hQCD_EE);
		TH1D* hData    = (TH1D*)hData_EB->Clone("hData");       hData   ->Add(hData_EE);
		TH1D* hRatio   = (TH1D*)histosZG["ZGRatio"+hshl]->Clone("hRatio");
		TH1D* hZnunu   = (TH1D*)histos["HadronicRegion"+hshl+"_Znunu"]->Clone("hZnunu");
	
		if(fReScaleAfterMT2){
			//this is discouraged - meaning  see at the beginning at fReScaleAfterMT2 setting
			TH1D *temp = (TH1D*)hPhotons->Clone("hPhotons_Clone");
			hPhotons = RescaleHisto(hPhotons, hQCD, histos["PhotonicSignalRegion_EB"+hshl+"_Data"  ], histos["PhotonicSignalRegion_EE"+hshl+"_Data"  ]);
			hQCD     = RescaleHisto(hQCD    , temp, histos["PhotonicSignalRegion_EB"+hshl+"_Data"  ], histos["PhotonicSignalRegion_EE"+hshl+"_Data"  ]);
			delete temp;
		}

		float modtableScaleFactor = 1.; float modtableScaleFactorError = -99.;
		if(signal_region[i2]=="2j1to2bmod" && HT_bin[i3]== "HTge450"){ modtableScaleFactor = 0.0937533; modtableScaleFactorError = 0.0224717; }
		if(signal_region[i2]=="3to5j1bmod" && HT_bin[i3]== "HTge450"){ modtableScaleFactor =   0.16062; modtableScaleFactorError = 0.0238675; }
		if(signal_region[i2]==   "6j1bmod" && HT_bin[i3]== "HTge450"){ modtableScaleFactor = 0.2693407; modtableScaleFactorError =  0.167893; }
		if(signal_region[i2]=="2j1to2bmod" && HT_bin[i3]== "HTge750"){ modtableScaleFactor = 0.0937533; modtableScaleFactorError = 0.0503878; }
		if(signal_region[i2]=="3to5j1bmod" && HT_bin[i3]== "HTge750"){ modtableScaleFactor =   0.16062; modtableScaleFactorError = 0.0181368; }
		if(signal_region[i2]==   "6j1bmod" && HT_bin[i3]== "HTge750"){ modtableScaleFactor = 0.2693407; modtableScaleFactorError =  0.142773; }
		if(signal_region[i2]=="2j1to2bmod" && HT_bin[i3]=="HTge1200"){ modtableScaleFactor = 0.0937533; modtableScaleFactorError = 0.0597187; }
		if(signal_region[i2]=="3to5j1bmod" && HT_bin[i3]=="HTge1200"){ modtableScaleFactor =   0.16062; modtableScaleFactorError = 0.0303366; }
		if(signal_region[i2]==   "6j1bmod" && HT_bin[i3]=="HTge1200"){ modtableScaleFactor = 0.2693407; modtableScaleFactorError =  0.201307; }
		for(int i = 1; i<=hData->GetNbinsX(); ++i){
			float PredictedZnunu                 = hRatio->GetBinContent(i)*(hData->GetBinContent(i)-hOther->GetBinContent(i)-hQCD->GetBinContent(i));
			float PredictedZnunu_ErrSys          = GetZnunuPreditionErrorSys (i, hData, hOther, hQCD, hRatio);
			float PredictedZnunu_ErrSysClosure   = GetZnunuPreditionErrorSysClosure (i, hPhotons, hRatio);
			float PredictedZnunu_ErrStat         = GetZnunuPreditionErrorStat(i, hData, hRatio);
			float PredictedZnunu_ErrSysRZG       = GetZnunuPreditionErrorSysRZG (i, hData, hOther, hQCD, hRatio);
			float PredictedZnunu_ErrSys1b0b      = -1.;
			if(modtableScaleFactorError>=0){//the factor is already applied!!
				PredictedZnunu               = PredictedZnunu * modtableScaleFactor;
				PredictedZnunu_ErrSys        = sqrt(pow(PredictedZnunu_ErrSys*modtableScaleFactor,2)+pow(PredictedZnunu/modtableScaleFactor*modtableScaleFactorError,2));
				PredictedZnunu_ErrSysClosure = PredictedZnunu_ErrSysClosure * modtableScaleFactor;
				PredictedZnunu_ErrStat       = PredictedZnunu_ErrStat * modtableScaleFactor;
				PredictedZnunu_ErrSysRZG     = PredictedZnunu_ErrSysRZG * modtableScaleFactor;
				PredictedZnunu_ErrSys1b0b    = PredictedZnunu / modtableScaleFactor * modtableScaleFactorError;
			}

			if(fPrintFullPrintOut){
		//	for(int nnn = 1; nnn<=hPhotons_EB->GetNbinsX();++nnn){
		//		*fLogStream << "bin " << hPhotons_EB->GetBinLowEdge(nnn) << "-" << hPhotons_EB->GetBinLowEdge(nnn)+hPhotons_EB->GetBinWidth(nnn) << ": ";
		//		*fLogStream << "EB " << hPhotons_EB->GetBinContent(nnn) << ", EE " << hPhotons_EE->GetBinContent(nnn) << endl;
		//	}
				*fLogStream << "For MT2=" << hPhotons->GetBinLowEdge(i) << "-" << hPhotons->GetBinLowEdge(i+1) <<  "(" << i << "/" << hData->GetNbinsX() << ")" << endl;
				*fLogStream << "Photon events in Photon Signal region (data-bg): " <<  hData->GetBinContent(i) << " - " << hOther->GetBinContent(i)+hQCD->GetBinContent(i) << " MC: " << hPhotons->GetBinContent(i) << endl;
				*fLogStream << "Predicted N Znunu events in Hadronic Signal region: ratio " << hRatio->GetBinContent(i) << " pm " <<hRatio->GetBinError(i) << endl;
				*fLogStream << "Prediction: "                                                                                            << endl;
				*fLogStream << "  " << PredictedZnunu << " pm " << PredictedZnunu_ErrSys  << " pm " << PredictedZnunu_ErrStat << " stat "        << endl;
				*fLogStream << "True N Znunu events:"                                                                                    << endl;
				*fLogStream << "  " << hZnunu->GetBinContent(i)                                                           << endl;
				*fLogStream << "MC closure:"                                                                                              << endl;
				*fLogStream << "  " << hRatio->GetBinContent(i)*(hPhotons->GetBinContent(i))   << " pm " << PredictedZnunu_ErrSysClosure << " (sys) " << endl; 
			}
			srbin.push_back(i2);
			if(fMET) htbin.push_back(i3);
			else     htbin.push_back(i3+1);
			mt2bin.push_back(i);
			mt2low.push_back(hData->GetBinLowEdge(i));
			mt2up.push_back(hData->GetBinLowEdge(i)+hData->GetBinWidth(i));
			if(i!=hData->GetNbinsX()) mt2last.push_back(0);
			else                      mt2last.push_back(1);
			znunugen.push_back(hZnunu->GetBinContent(i));
			znunugenerr.push_back(hZnunu->GetBinError(i));
			znunupred.push_back(PredictedZnunu);
			znunuprederrstat.push_back(PredictedZnunu_ErrStat);
			znunuprederrsyst.push_back(PredictedZnunu_ErrSys);
			znunuprederr.push_back(sqrt(pow(PredictedZnunu_ErrStat,2)+pow(PredictedZnunu_ErrSys,2) ) );
			znunuprederrsystRZG.push_back(PredictedZnunu_ErrSysRZG);
			znunuprederrsyst1b0b.push_back(PredictedZnunu_ErrSys1b0b);
			ZGratio.push_back(hRatio->GetBinContent(i));
			ZGratioerr.push_back(hRatio->GetBinError(i));
			numdata.push_back(hData->GetBinContent(i));
			numdataerr.push_back(hData->GetBinError(i));
			numqcd.push_back(hQCD->GetBinContent(i));
			numqcderr.push_back(hQCD->GetBinError(i));
			numother.push_back(hOther->GetBinContent(i));
			numothererr.push_back(hOther->GetBinError(i));
			numgamma.push_back(hPhotons->GetBinContent(i));
			numgammaerr.push_back(hPhotons->GetBinError(i));
			modtabSF.push_back(modtableScaleFactor);
			modtabSFerr.push_back(modtableScaleFactorError);
		}//for(int i = 1; i<=hData->GetNbinsX(); ++i)

		delete hPhotons_EB;
		delete hOther_EB;
		delete hQCD_EB;
		delete hData_EB;
		delete hPhotons_EE;
		delete hOther_EE;
		delete hQCD_EE;
		delete hData_EE;
		delete hPhotons;
		delete hOther;
		delete hQCD;
		delete hData;
		delete hRatio;
		delete hZnunu;

	}}//for HT and topological region bins

	//makes prediction card according to Bruno's wishes (7 TeV)
	if(PrintPredictionCard) PredictionCard(srbin, htbin, mt2bin, znunugen, znunugenerr, znunupred, znunuprederrstat, znunuprederrsyst, ZGratio, ZGratioerr, modtabSF, modtabSFerr);
	//saves a root file of prediction and MC truth and possibly makes plots
	if(SavePrediction) PredictionFile(srbin, htbin, mt2low, mt2up, mt2last, znunugen, znunugenerr, znunupred, znunuprederrstat, znunuprederrsyst, znunuprederr, znunuprederrsystRZG, znunuprederrsyst1b0b, PlotPrediction);
	//this one does the final prediction table
	if(PrintSummaryTable) SummaryTable(srbin, htbin, mt2low, mt2up, mt2last, znunugen, znunugenerr, numdata, numgamma, numgammaerr, numqcd, numqcderr, numother, numothererr, ZGratio, ZGratioerr, modtabSF, modtabSFerr, znunupred, znunuprederrstat, znunuprederrsyst, znunuprederr, znunuprederrsystRZG, znunuprederrsyst1b0b);
//NOTE: Add a function to store results into a root file with histograms, with possibility of plotting
//CONTINUE HERE
//one function for "Prediction card"
//one function for summary table
//one function for root file & plots

}

//this function produces the final result tables
void SummaryTable(vector<int> srbin, vector<int> htbin, vector<double> mt2low, vector<double> mt2up, vector<int> mt2last, vector<double> znunugen, vector<double> znunugenerr, vector<double> ndata, vector<double> ngamma, vector<double> ngammaerr, vector<double> nqcd, vector<double> nqcderr, vector<double> nother, vector<double> nothererr, vector<double> ZGratio, vector<double> ZGratioerr, vector<double> modtabSF, vector<double> modtabSFerr, vector<double> znunupred, vector<double> znunuprederrstat, vector<double> znunuprederrsyst, vector<double> znunuprederr, vector<double> znunuprederrsystRZG, vector<double> znunuprederrsyst1b0b){

	*fLogStream << "*********************************************************************" << endl;
	*fLogStream << "\%BEGINLATEX\%"                    << endl
		    << "\\begin{table}[!htb]"              << endl
		    << "\\begin{center}"                   << endl
		    << "\\begin{tabular}{|r|cccccc|}"      << endl
		    << "\\hline\\hline"                    << endl;
	*fLogStream << " $M_\\mathrm{T2}$ [GeV] & $N^{\\gamma}_{\\mathrm{data}}$         & $N^{QCD}_{\\mathrm{bg}}$&           $N^{other}_{\\mathrm{bg}}$      &    $R_{MC}(Z(\\nu\\bar{\\nu})/\\gamma)$ & data prediction               & sim. truth \\\\" << endl;
	//*fLogStream << "\\hline" << endl;// don't need this line
	string oldsr = "dummy";
	string oldHT = "dummy";
	for(unsigned int n = 0; n<srbin.size(); ++n){
		string sigreg;
		if(srbin[n]==0) sigreg="2 jets, 0 b jets";
		if(srbin[n]==1) sigreg="2 jets, $\\geq 1$ b jets";
		if(srbin[n]==2) sigreg="$3-5$ jets, 0 b jets";
		if(srbin[n]==3) sigreg="$3-5$ jets, 1 b jet";
		if(srbin[n]==4) sigreg="$3-5$ jets, 2 b jets";
		if(srbin[n]==5) sigreg="$\\geq 6$ jets, 0 b jets";
		if(srbin[n]==6) sigreg="$\\geq 6$ jets, 1 b jet";
		if(srbin[n]==7) sigreg="$\\geq 6$ jets, 2 b jets";
		if(srbin[n]==8) sigreg="$\\geq 3$ jets, $\\geq 3$ b jets";
		if(srbin[n]>8) continue;//these are the 'modified regions', i.e. 0b scaled by Z(1b)/Z(0b) --> done later!
		string htreg;
		if(htbin[n]==0 && fMET) htreg = "low $H_{\\mathrm{T}}$";
		if(htbin[n]==1 && fHT)  htreg = "medium $H_{\\mathrm{T}}$";
		if(htbin[n]==2 && fHT)  htreg = "high $H_{\\mathrm{T}}$";
		bool alreadyHTchange = false;
		if(htreg!=oldHT){//this if makes sure that HT region is only printed out, when the vector jumps to the next HT region
			*fLogStream << " \\hline\\hline " << endl << " & \\multicolumn{6}{l}{" <<  sigreg << ", " << htreg << "} \\\\ " << endl << "\\hline" << endl;
			oldHT = htreg;
			alreadyHTchange = true;
		}
		if(sigreg!=oldsr){//this if makes sure that topological region is only printed out, when the vector jumps topological the next HT region
			if(alreadyHTchange==false) *fLogStream << " \\hline\\hline " << endl << " & \\multicolumn{6}{l}{" <<  sigreg << ", " << htreg << "} \\\\ " << endl << "\\hline" << endl;
			oldsr = sigreg;
		}
		if(mt2last[n]==1)   *fLogStream  << " $"            <<                            "\\geq" << int(mt2low[n]) << "$" << " " << setw(7) << " & ";
		else                *fLogStream  << " $"            <<          int(mt2low[n]) << "-" << int(mt2up[n])      << "$" << " " << setw(4) << " & ";
		                    *fLogStream  << " " << setw(12) <<          int(ndata[n])                                      << " " << setw(1)  << "& ";
		                    *fLogStream  << fixed << setprecision(2);
		if(nqcd[ n]>0)      *fLogStream  << " " << setw(2)  <<   "$" << nqcd[ n]       << " \\pm "<< nqcderr[ n]    << "$" << " " << setw(1)  << "& ";
		else                *fLogStream  << " " << setw(2)  <<   "$" <<                     "-"                     << "$" << " " << setw(1)  << "& ";
		if(nother[ n]>0)    *fLogStream  << " " << setw(2)  <<   "$" << nother[ n]     << " \\pm "<< nothererr[ n]  << "$" << " " << setw(1)  << "& ";
		else                *fLogStream  << " " << setw(2)  <<   "$" <<                     "-"                     << "$" << " " << setw(1)  << "& ";
		                    *fLogStream  << " " << setw(4)  <<   "$" << ZGratio[ n]    << " \\pm "<< ZGratioerr[ n] << "$" << " " << setw(1)  << "& ";
		if(znunupred[ n]>0) *fLogStream  << " " << setw(7)  <<   "$" << znunupred[ n]  << " \\pm " << znunuprederrstat[ n] << " \\pm " << znunuprederrsyst[ n] << "$ & ";
		else                *fLogStream  << " " << setw(7)  <<   "$" <<                     "-"                     << "$" << " " << setw(1)  << "& ";
		                    *fLogStream  << " " << setw(4)  << znunugen[ n] << " \\\\ " << endl;
	}
	*fLogStream << "\\hline\\hline"  << endl
		    << "\\end{tabular}"  << endl
		    << "\\end{center}"   << endl
		    << "\\end{table}"    << endl
		    << "\%ENDLATEX\%"    << endl
		    << endl;
	if(!fDoPhotonSigmaIEtaIEta) cout << "FYI: YOUR ARE NOT DOING fDoPhotonSigmaIEtaIEta " << endl << endl;//just a reminder


	*fLogStream << "***** now the 1b/0b modified - not that MC truth is wrong, but can be copied from the above ******" << endl;
	*fLogStream << "\%BEGINLATEX\%"                    << endl
		    << "\\begin{table}[!htb]"              << endl
		    << "\\begin{center}"                   << endl
		    << "\\begin{tabular}{|r|cccccc|}"      << endl
		    << "\\hline\\hline"                    << endl;
	*fLogStream << " $M_\\mathrm{T2}$ [GeV] & $N^{\\gamma}_{\\mathrm{data}}$         & $N^{QCD}_{\\mathrm{bg}}$&           $R^{Z}_{\\mathrm{data}}(1b/0b)$      &    $R_{MC}(Z(\\nu\\bar{\\nu})/\\gamma)$ & data prediction               & sim. truth \\\\" << endl;
	//*fLogStream << "\\hline" << endl;// don't need this line
	for(unsigned int n = 0; n<srbin.size(); ++n){
		string sigreg;
		if(srbin[n]== 9) sigreg="2 jets, $\\geq 1$ b jets";
		if(srbin[n]==10) sigreg="$3-5$ jets, 1 b jet";
		if(srbin[n]==11) sigreg="$\\geq 6$ jets, 1 b jet";
		if(srbin[n]<= 8) continue;//these are the 'modified regions', i.e. 0b scaled by Z(1b)/Z(0b) --> done later!
		string htreg;
		if(htbin[n]==0 && fMET) htreg = "low $H_{\\mathrm{T}}$";
		if(htbin[n]==1 && fHT)  htreg = "medium $H_{\\mathrm{T}}$";
		if(htbin[n]==2 && fHT)  htreg = "high $H_{\\mathrm{T}}$";
		bool alreadyHTchange = false;
		if(htreg!=oldHT){//this if makes sure that HT region is only printed out, when the vector jumps to the next HT region
			*fLogStream << " \\hline\\hline " << endl << " & \\multicolumn{6}{l}{" <<  sigreg << ", " << htreg << "} \\\\ " << endl << "\\hline" << endl;
			oldHT = htreg;
			alreadyHTchange = true;
		}
		if(sigreg!=oldsr){//this if makes sure that topological region is only printed out, when the vector jumps topological the next HT region
			if(alreadyHTchange==false) *fLogStream << " \\hline\\hline " << endl << " & \\multicolumn{6}{l}{" <<  sigreg << ", " << htreg << "} \\\\ " << endl << "\\hline" << endl;
			oldsr = sigreg;
		}
		if( mt2last[n] ==1) *fLogStream  << " $"            <<                            "\\geq" << int(mt2low[n]) << "$" << " " << setw(7) << " & ";
		else                *fLogStream  << " $"            <<          int(mt2low[n]) << "-" << int(mt2up[n])      << "$" << " " << setw(4) << " & ";
		                    *fLogStream  << " " << setw(12) <<          int(ndata[ n])                                     << " " << setw(1)  << "& ";
		                    *fLogStream  << fixed << setprecision(2);
		if(nqcd[ n]>0)      *fLogStream  << " " << setw(2)  <<   "$" << nqcd[ n]       << " \\pm "<< nqcderr[ n]    << "$" << " " << setw(1)  << "& ";
		else                *fLogStream  << " " << setw(2)  <<   "$" <<                     "-"                     << "$" << " " << setw(1)  << "& ";
		                    *fLogStream  << fixed << setprecision(3);
		                    *fLogStream  << " " << setw(2)  <<   "$" << modtabSF[ n]   << " \\pm "<< modtabSFerr[ n]<< "$" << " " << setw(1)  << "& ";
		                    *fLogStream  << fixed << setprecision(2);
		                    *fLogStream  << " " << setw(4)  <<   "$" << ZGratio[ n]    << " \\pm "<< ZGratioerr[ n] << "$" << " " << setw(1)  << "& ";
		if(znunupred[ n]>0) *fLogStream  << " " << setw(7)  <<   "$" << znunupred[ n]  << " \\pm " << znunuprederrstat[ n] << " \\pm " << znunuprederrsyst[ n] << "$ & ";
		else                *fLogStream  << " " << setw(7)  <<   "$" <<                     "-"                     << "$" << " " << setw(1)  << "& ";
		                    *fLogStream  << " " << setw(4)  << znunugen[ n] << " \\\\ " << endl;
	}
	*fLogStream << "\\hline\\hline"  << endl
		    << "\\end{tabular}"  << endl
		    << "\\end{center}"   << endl
		    << "\\end{table}"    << endl
		    << "\%ENDLATEX\%"    << endl;
	*fLogStream << "***** this was the 1b/0b modified - not that MC truth is wrong, but can be copied from the above ******" << endl
		    << endl;

	if(!fDoPhotonSigmaIEtaIEta) cout << "FYI: YOUR ARE NOT DOING fDoPhotonSigmaIEtaIEta " << endl << endl;//just a reminder

}

//this stores prediction into File, so that one can use it for prediction purposes - also can do prediction vs. mctruth plots
void PredictionFile(vector<int> srbin, vector<int> htbin, vector<double> mt2low, vector<double> mt2up, vector<int> mt2last, vector<double> znunugen, vector<double> znunugenerr, vector<double> znunupred, vector<double> znunuprederrstat, vector<double> znunuprederrsyst, vector<double> znunuprederr, vector<double> znunuprederrsystRZG, vector<double> znunuprederrsyst1b0b, Bool_t PlotPrediction){

	TString filename = fOutDir;
	TString plotdirectory = fOutDir + "plots/";
	if(PlotPrediction) Util::MakeOutputDir(plotdirectory);
	filename = filename + "ZnunuFromGammaPredictionFile.root";

	map<string, TH1D*>    hs;
	//for now it is only truth +/- error
	//and prediction +/- error stored
	string sigreg; string regsig;
	vector<double> dp, dpe, st, ste, mt2binning; dp.clear(); dpe.clear(); st.clear(); ste.clear(); mt2binning.clear();
	vector <double> dpeStat, dpeSyst, dpeRZG, dpeZ1b0b; dpeStat.clear(); dpeSyst.clear(); dpeRZG.clear(); dpeZ1b0b.clear();
	string oldHT, oldsr;
	for(unsigned int n = 0; n<srbin.size(); ++n){
		if(srbin[n]==0) { regsig="2 jets, 0 b jets";         sigreg="2j0b";        }
		if(srbin[n]==1) { regsig="2 jets, #geq1 b jets";     sigreg="2j1to2b";     }
		if(srbin[n]==2) { regsig="3-5 jets, 0 b jets";       sigreg="3to5j0b";     }
		if(srbin[n]==3) { regsig="3-5 jets, 1 b jet";        sigreg="3to5j1b";     }
		if(srbin[n]==4) { regsig="3-5 jets, 2 b jets";       sigreg="3to5j2b";     }
		if(srbin[n]==5) { regsig="#geq6 jets, 0 b jets";     sigreg="6j0b";        }
		if(srbin[n]==6) { regsig="#geq6 jets, 1 b jet";      sigreg="6j1b";        }
		if(srbin[n]==7) { regsig="#geq6 jets, 2 b jets";     sigreg="6j2b";        }
		if(srbin[n]==8) { regsig="#geq3 jets, #geq3 b jets"; sigreg="3b";          }
		if(srbin[n]== 9){ regsig="2 jets, #geq1 b jets";     sigreg="2j1to2bmod";  }
		if(srbin[n]==10){ regsig="3-5 jets, 1 b jet";        sigreg="3to5j1bmod";  }
		if(srbin[n]==11){ regsig="#geq6 jets, 1 b jet";      sigreg="6j1bmod";     }
		string htreg; string reght;
		if(htbin[n]==0 && fMET) { reght = "low H_{T}";       htreg="lowHT";       }
		if(htbin[n]==1 && fHT)  { reght = "medium H_{T}";    htreg="mediumHT";    }
		if(htbin[n]==2 && fHT)  { reght = "high H_{T}";      htreg="highHT";      }
		bool alreadyHTchange = false;
		bool sigchange = false;
		if(reght!=oldHT){//this if makes sure that HT region is only printed out, when the vector jumps to the next HT region
			oldHT = reght;
		}
		if(regsig!=oldsr){//this if makes sure that topological region is only printed out, when the vector jumps topological the next HT region
			oldsr = regsig;
			sigchange = true;
		}
		mt2binning.push_back(mt2low[n]);
		dp.push_back(znunupred[n]);
		dpe.push_back(sqrt(pow(znunuprederrstat[n],2)+pow(znunuprederrsyst[n],2)) );
		dpeStat.push_back(znunuprederrstat[n]);
		dpeSyst.push_back(znunuprederrsyst[n]);
		dpeRZG.push_back(znunuprederrsystRZG[n]);
		double znunuprederrsyst1b0btemp = znunuprederrsyst1b0b[n];
		if(znunuprederrsyst1b0btemp<0) znunuprederrsyst1b0btemp = 0;
		dpeZ1b0b.push_back(znunuprederrsyst1b0btemp);
		st.push_back(znunugen[n]);
		ste.push_back(znunugenerr[n]);
		if(mt2last[n]==1) { 
			mt2binning.push_back(mt2up[n]);
			int Numbins = mt2binning.size()-1;
			int numarray = mt2binning.size();
			double binns[numarray];
			for(int ii = 0; ii<=Numbins; ++ii) binns[ii] = mt2binning[ii];
			string mapname = "Prediction_"+htreg+"_"+sigreg;
			if(hs.count(mapname) == 0 ) {hs[mapname] = new TH1D(mapname.c_str(), "", Numbins, binns); hs[mapname]->Sumw2();
						hs[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); hs[mapname]->GetYaxis()->SetTitle("events");
						hs[mapname]->SetMarkerStyle(20);
						hs[mapname]->SetMarkerColor(kBlack);
						hs[mapname]->SetLineColor(kBlack);
						hs[mapname]->SetLineWidth(3);
			}
			for(int ii = 1; ii<=Numbins; ++ii){
				hs[mapname]->SetBinContent(ii, dp[ii-1]);
				hs[mapname]->SetBinError(ii, dpe[ii-1]);
			}
			mapname = "Prediction_RZGonly_"+htreg+"_"+sigreg;
			if(hs.count(mapname) == 0 ) {hs[mapname] = new TH1D(mapname.c_str(), "", Numbins, binns); hs[mapname]->Sumw2();
						hs[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); hs[mapname]->GetYaxis()->SetTitle("events");
						hs[mapname]->SetMarkerStyle(20);
						hs[mapname]->SetMarkerColor(kBlack);
						hs[mapname]->SetLineColor(kBlack);
						hs[mapname]->SetLineWidth(3);
			}
			for(int ii = 1; ii<=Numbins; ++ii){
				hs[mapname]->SetBinContent(ii, dp[ii-1]);
				hs[mapname]->SetBinError(ii, dpeRZG[ii-1]);
			}
			mapname = "Prediction_Z1bZ0bonly_"+htreg+"_"+sigreg;
			if(hs.count(mapname) == 0 ) {hs[mapname] = new TH1D(mapname.c_str(), "", Numbins, binns); hs[mapname]->Sumw2();
						hs[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); hs[mapname]->GetYaxis()->SetTitle("events");
						hs[mapname]->SetMarkerStyle(20);
						hs[mapname]->SetMarkerColor(kBlack);
						hs[mapname]->SetLineColor(kBlack);
						hs[mapname]->SetLineWidth(3);
			}
			for(int ii = 1; ii<=Numbins; ++ii){
				hs[mapname]->SetBinContent(ii, dp[ii-1]);
				hs[mapname]->SetBinError(ii, dpeZ1b0b[ii-1]);
			}
			mapname = "Prediction_statonly_"+htreg+"_"+sigreg;
			if(hs.count(mapname) == 0 ) {hs[mapname] = new TH1D(mapname.c_str(), "", Numbins, binns); hs[mapname]->Sumw2();
						hs[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); hs[mapname]->GetYaxis()->SetTitle("events");
						hs[mapname]->SetMarkerStyle(20);
						hs[mapname]->SetMarkerColor(kBlack);
						hs[mapname]->SetLineColor(kBlack);
						hs[mapname]->SetLineWidth(3);
			}
			for(int ii = 1; ii<=Numbins; ++ii){
				hs[mapname]->SetBinContent(ii, dp[ii-1]);
				hs[mapname]->SetBinError(ii, dpeStat[ii-1]);
			}
			mapname = "Prediction_syst_"+htreg+"_"+sigreg;
			if(hs.count(mapname) == 0 ) {hs[mapname] = new TH1D(mapname.c_str(), "", Numbins, binns); hs[mapname]->Sumw2();
						hs[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); hs[mapname]->GetYaxis()->SetTitle("events");
						hs[mapname]->SetMarkerStyle(20);
						hs[mapname]->SetMarkerColor(kBlack);
						hs[mapname]->SetLineColor(kBlack);
						hs[mapname]->SetLineWidth(3);
			}
			for(int ii = 1; ii<=Numbins; ++ii){
				hs[mapname]->SetBinContent(ii, dp[ii-1]);
				hs[mapname]->SetBinError(ii, dpeSyst[ii-1]);
			}
			mapname = "SimulationTruth_"+htreg+"_"+sigreg;
			if(hs.count(mapname) == 0 ) { hs[mapname] = new TH1D(mapname.c_str(), "", Numbins, binns); hs[mapname]->Sumw2();
						hs[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); hs[mapname]->GetYaxis()->SetTitle("events");
						hs[mapname]->SetMarkerColor(kViolet-3);
						hs[mapname]->SetLineColor(kViolet-3);
						hs[mapname]->SetFillColor(kViolet-3);
						hs[mapname]->SetLineWidth(0);
						hs[mapname]->SetFillStyle(3002);
			} 
			for(int ii = 1; ii<=Numbins; ++ii){
				hs[mapname]->SetBinContent(ii, st[ii-1]);
				hs[mapname]->SetBinError(ii, ste[ii-1]);
			}
			mt2binning.clear();
			dp.clear(); dpe.clear();
			st.clear(); ste.clear();
		}
	}//all histograms are defined
	for(map<string,TH1D*>::iterator h=hs.begin(); h!=hs.end();++h){
		//general style
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
	}
	//now save everything
   	TFile *fpredictionfile = new TFile(filename.Data(),"RECREATE");
	fpredictionfile->cd();
	for(map<string,TH1D*>::iterator h=hs.begin(); h!=hs.end();++h){
		h->second->Write();
	}
	fpredictionfile->Close();
	*fLogStream << "Saved prediction/mc truth histograms in " << filename.Data() << endl;

	//now do the plotting
	if(PlotPrediction){
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
		entry->SetMarkerColor(kViolet-3);
		entry->SetLineColor(kViolet-3);
		entry->SetFillColor(kViolet-3);
		entry->SetLineWidth(0);
		entry->SetFillStyle(3002); 
   		TLegendEntry *entry2=leg->AddEntry("NULL","data prediction","p");
		entry2->SetLineColor(1);
		entry2->SetLineStyle(1);
		entry2->SetLineWidth(2);
		entry2->SetMarkerColor(1);
		entry2->SetMarkerStyle(20);
		entry2->SetMarkerSize(1);

		TString texttoplep, textht;
		string outname;
		double max(0.); double min(99999.);
		double maxp(0.); double minp(99999.);
		double maxs(0.); double mins(99999.);
		for(int ns = 0; ns < signalregionsize;++ns){
		for(int nh = 0; nh < HTbinsize; ++nh){
			if(fMET && nh==1) continue;
			if( fhighHT && nh==0) continue;
			if(!fhighHT && nh==1) continue;
			if(ns==0) { regsig="2 jets, 0 b jets";                 sigreg="2j0b";        }
			if(ns==1) { regsig="2 jets, $\\geq 1$ b jets";         sigreg="2j1to2b";     }
			if(ns==2) { regsig="$3-5$ jets, 0 b jets";             sigreg="3to5j0b";     }
			if(ns==3) { regsig="$3-5$ jets, 1 b jet";              sigreg="3to5j1b";     }
			if(ns==4) { regsig="$3-5$ jets, 2 b jets";             sigreg="3to5j2b";     }
			if(ns==5) { regsig="$\\geq 6$ jets, 0 b jets";         sigreg="6j0b";        }
			if(ns==6) { regsig="$\\geq 6$ jets, 1 b jet";          sigreg="6j1b";        }
			if(ns==7) { regsig="$\\geq 6$ jets, 2 b jets";         sigreg="6j2b";        }
			if(ns==8) { regsig="$\\geq 3$ jets, $\\geq 3$ b jets"; sigreg="3b";          }
			if(ns== 9){ regsig="2 jets, $\\geq 1$ b jets";         sigreg="2j1to2bmod";  }
			if(ns==10){ regsig="$3-5$ jets, 1 b jet";              sigreg="3to5j1bmod";  }
			if(ns==11){ regsig="$\\geq 6$ jets, 1 b jet";          sigreg="6j1bmod";     }
			string htreg; string reght;
			if(nh==0 && fMET) { reght = "low $H_{\\mathrm{T}}$";    htreg="lowHT";       }
			if(nh==0 && fHT)  { reght = "medium $H_{\\mathrm{T}}$"; htreg="mediumHT";    }
			if(nh==1 && fHT)  { reght = "high $H_{\\mathrm{T}}$";   htreg="highHT";      }
			string nameprediction = "Prediction_"+htreg+"_"+sigreg;
			string namesimulation = "SimulationTruth_"+htreg+"_"+sigreg;
			if(ns== 9) namesimulation = "SimulationTruth_"+htreg+"_2j1to2b";
			if(ns==10) namesimulation = "SimulationTruth_"+htreg+"_3to5j1b";
			if(ns==11) namesimulation = "SimulationTruth_"+htreg+"_6j1b";
			texttoplep = (TString)regsig;
			textht = (TString)reght;
			min = 0.;
			maxp = hs[nameprediction]->GetBinContent(hs[nameprediction]->GetMaximumBin())+hs[nameprediction]->GetBinError(hs[nameprediction]->GetMaximumBin());
			maxs = hs[namesimulation]->GetBinContent(hs[namesimulation]->GetMaximumBin())+hs[namesimulation]->GetBinError(hs[namesimulation]->GetMaximumBin());
			max  = (maxp>maxs)?maxp:maxs;
			max = 1.5*max;
			hs[nameprediction]->SetMaximum(max);
			hs[namesimulation]->SetMaximum(max);
			c1->Clear();
			c1->cd();
			hs[namesimulation]->SetFillStyle(3013);
			hs[namesimulation]->Draw("E2");
			hs[nameprediction]->Draw("sameE");
			tex->Draw();
			leg->Draw();
			toplep.DrawLatex(0.6,0.8548951,texttoplep.Data());
			ht.DrawLatex(0.6,0.791958,textht.Data());
			outname = plotdirectory+nameprediction + ".eps";
			c1->SaveAs(outname.c_str());
		}}
	}//if(PlotPrediction)
}

//this function makes a compact card of signal region vs MCtruth +/- err vs. pred. +/- err
void PredictionCard(vector<int> srbin, vector<int> htbin, vector<int> mt2bin, vector<double> znunugen, vector<double> znunugenerr, vector<double> znunupred, vector<double> znunuprederrstat, vector<double> znunuprederrsyst, vector<double> ZGratio, vector<double> ZGratioerr, vector<double> modtabSF, vector<double> modtabSFerr){

	*fLogStream << "*********************************************************************" << endl;
	*fLogStream << " Name     " << " SR " << " HTbin " << " MT2bin   " << " MCpred   " << "MCpredErr   " << " DataDrivenPred " << " DataDrivenPredError " << " ScaleFactor " << " ScaleFactorError " << endl;
	for(unsigned int n = 0; n<srbin.size(); ++n){
		if(srbin[n]>8) continue;//these are the 'modified regions', i.e. 0b scaled by Z(1b)/Z(0b) --> done later!
		*fLogStream << "ZnunuG " << setw(3) << srbin[n] << " " << setw(5) << htbin[n] << " " << setw(6) << int(mt2bin[n]-1) << " " << setw(14) << fixed << setprecision(4) << znunugen[n] <<  " " <<  setw(12) << znunugenerr[n] << " " <<  setw(12) << znunupred[n] << " " << setw(17) << sqrt( pow(znunuprederrstat[n],2) + pow(znunuprederrsyst[n],2) ) << " " << setw(15) << fixed << setprecision(4) << ZGratio[n] << " " << setw(14) << ZGratioerr[n] << endl;
	}
	*fLogStream << "****** now the 1b/0b modified - not that MC truth is wrong, but can be copied from the above *****" << endl;
	for(unsigned int n = 0; n<srbin.size(); ++n){
		if(srbin[n]<=8) continue;//these are the 'modified regions', i.e. 0b scaled by Z(1b)/Z(0b) --> done later!
		*fLogStream << "ZnunuG " << setw(3) << srbin[n] << " " << setw(5) << htbin[n] << " " << setw(6) << mt2bin[n]-1 << " " << setw(14) << fixed << setprecision(4) << znunugen[n] <<  " " <<  setw(12) << znunugenerr[n] << " " <<  setw(12) << znunupred[n] << " " << setw(17) << sqrt( pow(znunuprederrstat[n],2) + pow(znunuprederrsyst[n],2) ) << " " << setw(15) << fixed << setprecision(4) << ZGratio[n]*modtabSF[n] << " " << setw(14) << sqrt(pow(ZGratioerr[n]*modtabSF[n],2)+pow(ZGratio[n]*modtabSFerr[n],2)) << endl;
	}
	*fLogStream << endl << endl;

}

// histograms must be scales with prober uncertainties
float GetZnunuPreditionErrorSys(int i, TH1D* hData, TH1D* hOther, TH1D* hQCD, TH1D* ratio){

	if(i>hData->GetNbinsX() || i<=0) return -999.;
	float pred_err_2 = pow(ratio->GetBinError(i)   * (hData->GetBinContent(i) - hOther->GetBinContent(i) - hQCD->GetBinContent(i)),2);
	pred_err_2      += pow(ratio->GetBinContent(i) * hOther->GetBinError(i),2);
	pred_err_2      += pow(ratio->GetBinContent(i) * hQCD->GetBinError(i),2);

	return sqrt(pred_err_2);
}

// computes systematic uncerainty on Znunu prediction for MC closure test:
// takes into account MLfit uncertainty on Photon-MC yield, MC statistics, statistical uncerainty on Znunu/photon ratio
float GetZnunuPreditionErrorSysClosure(int i, TH1D* hPhotons, TH1D* ratio){

	if(i>hPhotons->GetNbinsX() || i<=0) return -999.;
	float pred_err_2 = pow(ratio->GetBinError(i)     * hPhotons->GetBinContent(i),2);
	pred_err_2      += pow(ratio->GetBinContent(i)   * hPhotons->GetBinError(i)  ,2);
	return sqrt(pred_err_2);
}

// statistical uncertainty on prediction due to data statistics
float GetZnunuPreditionErrorStat(int i, TH1D* hData,  TH1D* ratio){

	if(i>hData->GetNbinsX() || i<=0) return -999.;
	float pred_err   = ratio->GetBinContent(i) * hData->GetBinError(i);
	return pred_err;
}

float GetZnunuPreditionErrorSysRZG(int i, TH1D* hData, TH1D* hOther, TH1D* hQCD, TH1D* ratio){
	if(i>hData->GetNbinsX() || i<=0) return -999.;
	float pred_err_2 = pow(ratio->GetBinError(i)  * (hData->GetBinContent(i) - hOther->GetBinContent(i) - hQCD->GetBinContent(i)),2);

	return sqrt(pred_err_2);
}

//rebin histogram
TH1D* RebinThisHistogram(TH1D* histogram, int nrebinnedbins, double* rebinnedbins){

	string name = histogram->GetName();
	histogram->SetName((name + "_original").c_str());
	TH1D* temphists = (TH1D*)histogram->Rebin(nrebinnedbins, (name + "_rebinned").c_str(), rebinnedbins);
	temphists->SetName((name).c_str());
	return temphists;

}

// takes histo, scale factor and uncertainty on scale factor
// and returns scaled histo with uncertainty on scale factor propagated.
TH1D* GetScaledHisto(TH1D* histo, float scale_fact, float scale_fact_err){
	TH1D *h        = (TH1D*) histo->Clone(histo->GetName());
	for(int i = 1; i<=h->GetNbinsX(); ++i){
	h->SetBinError(i, sqrt(h->GetBinError(i)*  h->GetBinError(i)   *scale_fact    *scale_fact + 
			  h->GetBinContent(i)*h->GetBinContent(i) *scale_fact_err*scale_fact_err));
	h->SetBinContent(i, h->GetBinContent(i)*scale_fact);
	}
	return h;
}

//rescale histMC1+histMC2 to event yield of hist_data_EB+hist_data_EE
TH1D* RescaleHisto(TH1D* histMC1, TH1D* histMC2, TH1D* hist_data_EB, TH1D* hist_data_EE){
	// takes histo, scale factor and uncertainty on scale factor
	// and returns scaled and rebinned histo with 1 bin with uncertainty on scale factor propagated.
	TH1D *hData    = (TH1D*) hist_data_EB->Clone("hData");
	TH1D *hDataEE  = (TH1D*) hist_data_EE->Clone("hData_EE");
	TH1D *h        = (TH1D*) histMC1->Clone(histMC1->GetName());
	TH1D *hMC      = (TH1D*) histMC1->Clone("hMC");
	TH1D *htemp    = (TH1D*) histMC2->Clone("htemp");
	hData->Add(hDataEE);
	hMC  ->Add(htemp);
	double scalefactor = 1;
	if(hMC->Integral()>0&&hData->Integral()>0) scalefactor = hData->Integral()/hMC->Integral();
	if(fVerbose>4) cout << "after MT2>" << hData->GetBinLowEdge(1) << ": data " << hData->Integral() << ", MC " << hMC->Integral() << " --> SF = " << scalefactor << endl;
	h->Scale(scalefactor);
	delete hData;
	delete hDataEE;
	delete hMC;
	delete htemp;
	return h;
}

//this function adds up the contents and uncertainties for bins > lastbinlowedge
//and replaces those contents and uncertainties by their (quadratic) sum
TH1D* RefillRatioHisto(TH1D* histo, float lastbinlowedge){
	TH1D *h        = (TH1D*) histo->Clone(histo->GetName());
	if(h->GetNbinsX()<2) return h;
	if(h->GetBinLowEdge(h->GetNbinsX())<lastbinlowedge) return h;
	double content(0.), error2(0.);
	for(int i = 1; i<=h->GetNbinsX(); ++i){
		//note: int should make sure that no rounding error screws this up
		if((int)h->GetBinLowEdge(i)>=(int)lastbinlowedge){
			content += h->GetBinContent(i);
			error2  += pow(h->GetBinError(i),2);
		}
	}
	for(int i = 1; i<=h->GetNbinsX(); ++i){
		//note: int should make sure that no rounding error screws this up
		if((int)h->GetBinLowEdge(i)>=(int)lastbinlowedge){
			h->SetBinContent(i, content);
			h->SetBinError(i, sqrt(error2));
		}
	}
	return h;
}

//this functions add an uncertainty of bulk_adduncertainty/tail_adduncertainty to 
//the histograms uncertainty, separated by separatingbin_lowedge
//if a bin is mostly in 'tail' although lowedge < separatingbin_lowedge it is assessed to the tail
TH1D* RescaleUncertaintyHisto(TH1D* histo, float bulk_adduncertainty, float tail_adduncertainty, float separatingbin_lowedge){

	TH1D *h        = (TH1D*) histo->Clone(histo->GetName());
	//note: int should make sure that no rounding error screws this up
	for(int i = 1; i<=h->GetNbinsX(); ++i){
		//if bin is above threshold, or last bin
		if((int)h->GetBinLowEdge(i)>=(int)separatingbin_lowedge || i==h->GetNbinsX() )
			h->SetBinError(i, sqrt(pow(h->GetBinError(i),2)+pow(h->GetBinContent(i)*tail_adduncertainty,2)));
		//if bin is last bin below threshold, but its majority is above
		else if((int)h->GetBinLowEdge(i+1)>=(int)separatingbin_lowedge && (int)h->GetBinLowEdge(i)<(int)separatingbin_lowedge &&
		(h->GetBinLowEdge(i+1)-separatingbin_lowedge)>=(separatingbin_lowedge-h->GetBinLowEdge(i) ) )
			h->SetBinError(i, sqrt(pow(h->GetBinError(i),2)+pow(h->GetBinContent(i)*tail_adduncertainty,2)));
		else
			h->SetBinError(i, sqrt(pow(h->GetBinError(i),2)+pow(h->GetBinContent(i)*bulk_adduncertainty,2)));
	}
	return h;
}


