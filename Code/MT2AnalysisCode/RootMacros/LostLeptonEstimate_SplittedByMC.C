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
#include "TPad.h"
#include "TLine.h"
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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path for MT2tree.hh
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/Utilities.hh"

//use via root -l LostLeptonEstimate_SplittedbyMC.C++

using namespace std;

//this code takes the lost lepton estimates produced with LostLeptonEstimate.C (actually at the moment hard-coded)
//and the MT2 shapes produced with  TTbarStudies.C (or TTbarStudiesISR.C)
//and reweights the shapes to fit the normalization obtained by the prediction.
//This is then the final prediction of the LostLeptonEstimate (excluding shape uncertainties that come later into the game)
//the ouputfile is used for producing final result tables and figures
//Also a 'final' table for the LostLeptonEstimate is produced
void LostLeptonEstimate_SplittedByMC();

const int fVerbose = 3;

const int HTregionsize = 3;
string HT_region[HTregionsize] = {"lowHT", "mediumHT", "highHT"};
const int signalregionsize = 9;
string signal_region[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b"};

std::ostringstream* fLogStream     = 0;

void LostLeptonEstimate_SplittedByMC(){

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	bool fixemptybin = true; //if true use MC +/- 100% if data prediction bin is empty. Note that this is checked separately for ele,muo, and tau - default = true

	bool onlyttbar = false;//if true: onlytop = false, onlyW = false - use only ttbar shape for splitting LostLeptonEstimate - default = false
	bool onlytop   = false;//if true: onlyW = false - use only top (ttbar+single top) shape for splitting LostLeptonEstimate - default = false
	bool onlyW     = false;// - use only W shape for splitting LostLeptonEstimate - default = false

	//you should remember this from LostLeptonEstimate.C
	fLogStream = new std::ostringstream();
	Bool_t  fWriteToFile               = false; // writes couts to file
	Bool_t  fAppend                    = true; // append at end of file (if existing), otherwise delete previous content


	bool WnoScaleMatching = true;//true if want W matchingup/down scaleup/down, Wincl sample in Nominal1Sample (default == false due to statistics of those W samples) // here only in plotting
				     //actually this flag is irrelevant - unless the TFile with this flag set false does not exist
	bool ISRusage = true;//For the shape use as default shape the MT2 shape after 'ISR recipe' has been applied, default is true.

	//you need the shapes done with TTbarStudies.C
	//you should also implement to get the LostLeptonEstimate from a TFile, as at the moment this is hard-coded
	TString   inputdir = "Filtered/TTbarStudies/";
	TString  inputname = "TTbarStudiesHistograms_all.root";
	if(WnoScaleMatching)  inputname = "TTbarStudiesHistograms_all_noWscaleupdown.root";
	if(ISRusage) inputname = "TTbarStudiesHistograms_all_noWscaleupdown_ISRdefault.root";

	TFile *infile = TFile::Open(inputdir + inputname);

	//save everything
	TString outputname = "LostLeptonEstimate_SplittedbyMC.root";
	if(ISRusage) outputname = "LostLeptonEstimate_SplittedbyMC_ISR.root";
	if(fixemptybin){ 
		outputname = "LostLeptonEstimate_SplittedbyMC_fixemptybin.root";
		if(ISRusage) outputname = "LostLeptonEstimate_SplittedbyMC_fixemptybin_ISR.root";
	}

	map<string, TH1D*>    histos;
	map<string, TH1D*>    histosnormalized;

	string samplekind;
	if(onlyttbar)    samplekind = "TTbar";
	else if(onlytop) samplekind = "allTop";
	else if(onlyW)   samplekind = "WJets";
	else             samplekind = "allMC";

	//load the shapes
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		string hs = string("_") + "Nominal" + string("_") + samplekind + string("_") + HT_region[i3] + string("_") + signal_region[i2];
		string mapname = "MT2" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)infile->Get(mapname.c_str());
		string mapnamenorm = "Norm"+mapname;
		if(histosnormalized.count(mapnamenorm) == 0 ) histosnormalized[mapnamenorm] = (TH1D*)infile->Get(mapnamenorm.c_str());
	}}

   
	//load the estimates - taken from LLVisualization.C - rewrite this section to read in the histograms from the rootfile produced by LLVisualization.C
    TH1D *lHT_ele = new TH1D("lHT_ele","",9,0,9); lHT_ele->SetMarkerStyle(20), lHT_ele->SetMarkerColor(kBlack); lHT_ele->SetLineWidth(3); lHT_ele->SetLineColor(kBlack);
    TH1D *lHT_muo = new TH1D("lHT_muo","",9,0,9); lHT_muo->SetMarkerStyle(20), lHT_muo->SetMarkerColor(kBlack); lHT_muo->SetLineWidth(3); lHT_muo->SetLineColor(kBlack);
    TH1D *lHT_tau = new TH1D("lHT_tau","",9,0,9); lHT_tau->SetMarkerStyle(20), lHT_tau->SetMarkerColor(kBlack); lHT_tau->SetLineWidth(3); lHT_tau->SetLineColor(kBlack);
    TH1D *mHT_ele = new TH1D("mHT_ele","",9,0,9); mHT_ele->SetMarkerStyle(20), mHT_ele->SetMarkerColor(kBlack); mHT_ele->SetLineWidth(3); mHT_ele->SetLineColor(kBlack);
    TH1D *mHT_muo = new TH1D("mHT_muo","",9,0,9); mHT_muo->SetMarkerStyle(20), mHT_muo->SetMarkerColor(kBlack); mHT_muo->SetLineWidth(3); mHT_muo->SetLineColor(kBlack);
    TH1D *mHT_tau = new TH1D("mHT_tau","",9,0,9); mHT_tau->SetMarkerStyle(20), mHT_tau->SetMarkerColor(kBlack); mHT_tau->SetLineWidth(3); mHT_tau->SetLineColor(kBlack);
    TH1D *hHT_ele = new TH1D("hHT_ele","",9,0,9); hHT_ele->SetMarkerStyle(20), hHT_ele->SetMarkerColor(kBlack); hHT_ele->SetLineWidth(3); hHT_ele->SetLineColor(kBlack);
    TH1D *hHT_muo = new TH1D("hHT_muo","",9,0,9); hHT_muo->SetMarkerStyle(20), hHT_muo->SetMarkerColor(kBlack); hHT_muo->SetLineWidth(3); hHT_muo->SetLineColor(kBlack);
    TH1D *hHT_tau = new TH1D("hHT_tau","",9,0,9); hHT_tau->SetMarkerStyle(20), hHT_tau->SetMarkerColor(kBlack); hHT_tau->SetLineWidth(3); hHT_tau->SetLineColor(kBlack);
    TH1D *lHT_eleMC = new TH1D("lHT_eleMC","",9,0,9); lHT_eleMC->SetFillStyle(3001); lHT_eleMC->SetFillColor(kBlue);
    TH1D *lHT_muoMC = new TH1D("lHT_muoMC","",9,0,9); lHT_muoMC->SetFillStyle(3001); lHT_muoMC->SetFillColor(kBlue);
    TH1D *lHT_tauMC = new TH1D("lHT_tauMC","",9,0,9); lHT_tauMC->SetFillStyle(3001); lHT_tauMC->SetFillColor(kBlue);
    TH1D *mHT_eleMC = new TH1D("mHT_eleMC","",9,0,9); mHT_eleMC->SetFillStyle(3001); mHT_eleMC->SetFillColor(kBlue);
    TH1D *mHT_muoMC = new TH1D("mHT_muoMC","",9,0,9); mHT_muoMC->SetFillStyle(3001); mHT_muoMC->SetFillColor(kBlue);
    TH1D *mHT_tauMC = new TH1D("mHT_tauMC","",9,0,9); mHT_tauMC->SetFillStyle(3001); mHT_tauMC->SetFillColor(kBlue);
    TH1D *hHT_eleMC = new TH1D("hHT_eleMC","",9,0,9); hHT_eleMC->SetFillStyle(3001); hHT_eleMC->SetFillColor(kBlue);
    TH1D *hHT_muoMC = new TH1D("hHT_muoMC","",9,0,9); hHT_muoMC->SetFillStyle(3001); hHT_muoMC->SetFillColor(kBlue);
    TH1D *hHT_tauMC = new TH1D("hHT_tauMC","",9,0,9); hHT_tauMC->SetFillStyle(3001); hHT_tauMC->SetFillColor(kBlue);

    // error on pred: stat + syst + ISR + doubleLL
    lHT_ele  ->SetBinContent(1,215.07); lHT_ele  ->SetBinError(1,sqrt(pow(12.46,2)+pow(27.09,2)+pow(4.02,2)+pow(0.1,2)));
    lHT_ele  ->SetBinContent(2, 18.84); lHT_ele  ->SetBinError(2,sqrt(pow( 3.59,2)+pow( 2.99,2)+pow(0.37,2)+pow(0.2,2)));
    lHT_ele  ->SetBinContent(3,399.46); lHT_ele  ->SetBinError(3,sqrt(pow(19.51,2)+pow(43.41,2)+pow(7.28,2)+pow(0.36,2)));
    lHT_ele  ->SetBinContent(4,128.33); lHT_ele  ->SetBinError(4,sqrt(pow(11.26,2)+pow(14.73,2)+pow(2.25,2)+pow(0.56,2)));
    lHT_ele  ->SetBinContent(5, 29.55); lHT_ele  ->SetBinError(5,sqrt(pow( 5.98,2)+pow( 3.74,2)+pow(0.45,2)+pow(0.53,2)));
    lHT_ele  ->SetBinContent(6,  9.70); lHT_ele  ->SetBinError(6,sqrt(pow( 3.50,2)+pow( 1.58,2)+pow(0.01,2)+pow(0.32,2)));
    lHT_ele  ->SetBinContent(7, 15.81); lHT_ele  ->SetBinError(7,sqrt(pow( 4.56,2)+pow( 2.36,2)+pow(0.00,2)+pow(0.43,2)));
    lHT_ele  ->SetBinContent(8,  9.66); lHT_ele  ->SetBinError(8,sqrt(pow( 3.22,2)+pow( 1.75,2)+pow(0.07,2)+pow(0.01,2)));
    lHT_ele  ->SetBinContent(9,  4.79); lHT_ele  ->SetBinError(9,sqrt(pow( 2.14,2)+pow( 1.10,2)+pow(0.02,2)+pow(0.07,2)));
    lHT_eleMC->SetBinContent(1,247.94); lHT_eleMC->SetBinError(1, 4.93);
    lHT_eleMC->SetBinContent(2, 19.37); lHT_eleMC->SetBinError(2, 2.09);
    lHT_eleMC->SetBinContent(3,404.54); lHT_eleMC->SetBinError(3, 9.28);
    lHT_eleMC->SetBinContent(4,113.31); lHT_eleMC->SetBinError(4, 5.93);
    lHT_eleMC->SetBinContent(5, 38.79); lHT_eleMC->SetBinError(5, 2.84);
    lHT_eleMC->SetBinContent(6, 11.04); lHT_eleMC->SetBinError(6, 1.21);
    lHT_eleMC->SetBinContent(7, 14.37); lHT_eleMC->SetBinError(7, 1.29);
    lHT_eleMC->SetBinContent(8,  8.99); lHT_eleMC->SetBinError(8, 0.89);
    lHT_eleMC->SetBinContent(9,  4.66); lHT_eleMC->SetBinError(9, 0.78);
  
    lHT_muo  ->SetBinContent(1,174.00); lHT_muo  ->SetBinError(1,sqrt(pow( 8.55,2)+pow(30.66,2)+pow(4.00,2)+pow(0.60,2)));
    lHT_muo  ->SetBinContent(2, 17.60); lHT_muo  ->SetBinError(2,sqrt(pow( 2.60,2)+pow( 3.72,2)+pow(0.35,2)+pow(0.14,2)));
    lHT_muo  ->SetBinContent(3,351.08); lHT_muo  ->SetBinError(3,sqrt(pow(14.26,2)+pow(50.43,2)+pow(7.40,2)+pow(0.82,2)));
    lHT_muo  ->SetBinContent(4,106.24); lHT_muo  ->SetBinError(4,sqrt(pow( 8.05,2)+pow(20.92,2)+pow(0.26,2)+pow(0.51,2)));
    lHT_muo  ->SetBinContent(5, 28.73); lHT_muo  ->SetBinError(5,sqrt(pow( 4.25,2)+pow( 4.65,2)+pow(0.48,2)+pow(0.59,2)));
    lHT_muo  ->SetBinContent(6,  8.33); lHT_muo  ->SetBinError(6,sqrt(pow( 2.32,2)+pow( 1.51,2)+pow(0.17,2)+pow(0.03,2)));
    lHT_muo  ->SetBinContent(7, 11.24); lHT_muo  ->SetBinError(7,sqrt(pow( 2.81,2)+pow( 1.95,2)+pow(0.03,2)+pow(0.28,2)));
    lHT_muo  ->SetBinContent(8,  1.93); lHT_muo  ->SetBinError(8,sqrt(pow( 1.13,2)+pow( 0.37,2)+pow(0.02,2)+pow(0.18,2)));
    lHT_muo  ->SetBinContent(9,  6.76); lHT_muo  ->SetBinError(9,sqrt(pow( 2.06,2)+pow( 1.57,2)+pow(0.06,2)+pow(0.02,2)));
    lHT_muoMC->SetBinContent(1,203.15); lHT_muoMC->SetBinError(1, 5.25);
    lHT_muoMC->SetBinContent(2, 14.01); lHT_muoMC->SetBinError(2, 1.88);
    lHT_muoMC->SetBinContent(3,340.25); lHT_muoMC->SetBinError(3, 9.84);
    lHT_muoMC->SetBinContent(4,105.22); lHT_muoMC->SetBinError(4,14.10);
    lHT_muoMC->SetBinContent(5, 33.30); lHT_muoMC->SetBinError(5, 3.12);
    lHT_muoMC->SetBinContent(6,  8.99); lHT_muoMC->SetBinError(6, 0.98);
    lHT_muoMC->SetBinContent(7, 11.70); lHT_muoMC->SetBinError(7, 1.28);
    lHT_muoMC->SetBinContent(8,  7.42); lHT_muoMC->SetBinError(8, 0.84);
    lHT_muoMC->SetBinContent(9,  4.02); lHT_muoMC->SetBinError(9, 0.75);

    lHT_tau  ->SetBinContent(1,274.04); lHT_tau  ->SetBinError(1,sqrt(pow(49.46,2)+pow(32.70,2)+pow(4.35,2)+pow(0.22,2)));
    lHT_tau  ->SetBinContent(2, 14.99); lHT_tau  ->SetBinError(2,sqrt(pow(13.74,2)+pow( 4.24,2)+pow(0.09,2)+pow(0.91,2)));
    lHT_tau  ->SetBinContent(3,468.93); lHT_tau  ->SetBinError(3,sqrt(pow(42.55,2)+pow(54.41,2)+pow(7.93,2)+pow(0.24,2)));
    lHT_tau  ->SetBinContent(4,170.09); lHT_tau  ->SetBinError(4,sqrt(pow(27.99,2)+pow(24.89,2)+pow(2.79,2)+pow(0.12,2)));
    lHT_tau  ->SetBinContent(5, 84.82); lHT_tau  ->SetBinError(5,sqrt(pow(29.44,2)+pow(16.38,2)+pow(1.15,2)+pow(1.18,2)));
    lHT_tau  ->SetBinContent(6, 29.44); lHT_tau  ->SetBinError(6,sqrt(pow( 8.12,2)+pow( 5.05,2)+pow(0.02,2)+pow(0.15,2)));
    lHT_tau  ->SetBinContent(7, 18.85); lHT_tau  ->SetBinError(7,sqrt(pow( 6.95,2)+pow( 3.62,2)+pow(0.11,2)+pow(0.11,2)));
    lHT_tau  ->SetBinContent(8,  7.09); lHT_tau  ->SetBinError(8,sqrt(pow( 5.53,2)+pow( 1.40,2)+pow(0.02,2)+pow(0.16,2)));
    lHT_tau  ->SetBinContent(9,  8.12); lHT_tau  ->SetBinError(9,sqrt(pow( 6.59,2)+pow( 2.25,2)+pow(0.23,2)+pow(0.49,2))); 
    lHT_tauMC->SetBinContent(1,262.54); lHT_tauMC->SetBinError(1, 5.22);
    lHT_tauMC->SetBinContent(2, 20.20); lHT_tauMC->SetBinError(2, 2.07);
    lHT_tauMC->SetBinContent(3,452.74); lHT_tauMC->SetBinError(3,10.89);
    lHT_tauMC->SetBinContent(4,157.92); lHT_tauMC->SetBinError(4,16.79);
    lHT_tauMC->SetBinContent(5, 44.68); lHT_tauMC->SetBinError(5, 4.66);
    lHT_tauMC->SetBinContent(6, 15.72); lHT_tauMC->SetBinError(6, 1.42);
    lHT_tauMC->SetBinContent(7, 18.19); lHT_tauMC->SetBinError(7, 1.52);
    lHT_tauMC->SetBinContent(8, 13.98); lHT_tauMC->SetBinError(8, 1.41);
    lHT_tauMC->SetBinContent(9,  7.23); lHT_tauMC->SetBinError(9, 1.16);
    

    mHT_ele  ->SetBinContent(1, 76.94); mHT_ele  ->SetBinError(1,sqrt(pow( 7.14,2)+pow(10.52,2)+pow(3.19,2)+pow(0.20,2)));
    mHT_ele  ->SetBinContent(2,  5.94); mHT_ele  ->SetBinError(2,sqrt(pow( 2.00,2)+pow( 1.04,2)+pow(0.33,2)+pow(0.19,2)));
    mHT_ele  ->SetBinContent(3,126.84); mHT_ele  ->SetBinError(3,sqrt(pow(10.27,2)+pow(15.07,2)+pow(3.26,2)+pow(0.12,2)));
    mHT_ele  ->SetBinContent(4, 46.06); mHT_ele  ->SetBinError(4,sqrt(pow( 5.86,2)+pow( 8.85,2)+pow(0.80,2)+pow(1.54,2)));
    mHT_ele  ->SetBinContent(5, 19.15); mHT_ele  ->SetBinError(5,sqrt(pow( 4.33,2)+pow( 2.70,2)+pow(0.38,2)+pow(0.44,2)));
    mHT_ele  ->SetBinContent(6,  8.60); mHT_ele  ->SetBinError(6,sqrt(pow( 3.25,2)+pow( 1.33,2)+pow(0.08,2)+pow(0.20,2)));
    mHT_ele  ->SetBinContent(7, 19.44); mHT_ele  ->SetBinError(7,sqrt(pow( 4.95,2)+pow( 2.82,2)+pow(0.12,2)+pow(0.16,2)));
    mHT_ele  ->SetBinContent(8, 21.57); mHT_ele  ->SetBinError(8,sqrt(pow( 5.59,2)+pow( 2.74,2)+pow(0.10,2)+pow(0.62,2)));
    mHT_ele  ->SetBinContent(9,  3.02); mHT_ele  ->SetBinError(9,sqrt(pow( 1.78,2)+pow( 0.58,2)+pow(0.01,2)+pow(0.49,2)));  
    mHT_eleMC->SetBinContent(1,106.35); mHT_eleMC->SetBinError(1, 3.19);
    mHT_eleMC->SetBinContent(2, 13.02); mHT_eleMC->SetBinError(2, 1.35);
    mHT_eleMC->SetBinContent(3,164.25); mHT_eleMC->SetBinError(3, 4.76);
    mHT_eleMC->SetBinContent(4, 56.24); mHT_eleMC->SetBinError(4, 3.73);
    mHT_eleMC->SetBinContent(5, 28.95); mHT_eleMC->SetBinError(5, 2.53);
    mHT_eleMC->SetBinContent(6, 15.03); mHT_eleMC->SetBinError(6, 1.42);
    mHT_eleMC->SetBinContent(7, 21.58); mHT_eleMC->SetBinError(7, 1.71);
    mHT_eleMC->SetBinContent(8, 26.12); mHT_eleMC->SetBinError(8, 2.34);
    mHT_eleMC->SetBinContent(9,  8.19); mHT_eleMC->SetBinError(9, 1.41);
 
    mHT_muo  ->SetBinContent(1, 78.12); mHT_muo  ->SetBinError(1,sqrt(pow( 5.92,2)+pow(13.82,2)+pow(3.23,2)+pow(0.03,2)));
    mHT_muo  ->SetBinContent(2, 13.99); mHT_muo  ->SetBinError(2,sqrt(pow( 2.62,2)+pow( 2.80,2)+pow(0.83,2)+pow(0.05,2)));
    mHT_muo  ->SetBinContent(3,122.86); mHT_muo  ->SetBinError(3,sqrt(pow( 8.12,2)+pow(19.05,2)+pow(3.51,2)+pow(0.44,2)));
    mHT_muo  ->SetBinContent(4, 56.04); mHT_muo  ->SetBinError(4,sqrt(pow( 5.70,2)+pow( 8.79,2)+pow(0.16,2)+pow(0.68,2)));
    mHT_muo  ->SetBinContent(5, 20.46); mHT_muo  ->SetBinError(5,sqrt(pow( 3.46,2)+pow( 3.43,2)+pow(0.24,2)+pow(1.19,2)));
    mHT_muo  ->SetBinContent(6, 11.21); mHT_muo  ->SetBinError(6,sqrt(pow( 2.82,2)+pow( 1.91,2)+pow(0.26,2)+pow(0.20,2)));
    mHT_muo  ->SetBinContent(7, 21.37); mHT_muo  ->SetBinError(7,sqrt(pow( 4.00,2)+pow( 3.33,2)+pow(0.21,2)+pow(0.14,2)));
    mHT_muo  ->SetBinContent(8, 14.96); mHT_muo  ->SetBinError(8,sqrt(pow( 3.12,2)+pow( 2.41,2)+pow(0.00,2)+pow(0.05,2)));
    mHT_muo  ->SetBinContent(9, 11.17); mHT_muo  ->SetBinError(9,sqrt(pow( 2.89,2)+pow( 2.25,2)+pow(0.06,2)+pow(0.01,2)));   
    mHT_muoMC->SetBinContent(1, 96.83); mHT_muoMC->SetBinError(1, 2.89);
    mHT_muoMC->SetBinContent(2, 12.83); mHT_muoMC->SetBinError(2, 1.58);
    mHT_muoMC->SetBinContent(3,137.49); mHT_muoMC->SetBinError(3, 4.37);
    mHT_muoMC->SetBinContent(4, 50.56); mHT_muoMC->SetBinError(4, 5.13);
    mHT_muoMC->SetBinContent(5, 24.55); mHT_muoMC->SetBinError(5, 4.39);
    mHT_muoMC->SetBinContent(6, 12.27); mHT_muoMC->SetBinError(6, 1.33);
    mHT_muoMC->SetBinContent(7, 17.53); mHT_muoMC->SetBinError(7, 1.34);
    mHT_muoMC->SetBinContent(8, 16.71); mHT_muoMC->SetBinError(8, 1.60);
    mHT_muoMC->SetBinContent(9,  7.82); mHT_muoMC->SetBinError(9, 1.24);

    mHT_tau  ->SetBinContent(1,121.86); mHT_tau  ->SetBinError(1,sqrt(pow(33.40,2)+pow(15.75,2)+pow(2.21,2)+pow(0.04,2)));
    //mHT_tau  ->SetBinContent(2,); mHT_tau  ->SetBinError(2,sqrt(pow( ,2)+pow(,2)+pow(,2)+pow(,2)));
    mHT_tau  ->SetBinContent(3,155.19); mHT_tau  ->SetBinError(3,sqrt(pow(23.89,2)+pow(18.69,2)+pow(3.71,2)+pow(0.37,2)));
    mHT_tau  ->SetBinContent(4, 78.34); mHT_tau  ->SetBinError(4,sqrt(pow(20.67,2)+pow(11.16,2)+pow(0.37,2)+pow(0.73,2)));
    mHT_tau  ->SetBinContent(5, 28.22); mHT_tau  ->SetBinError(5,sqrt(pow(12.58,2)+pow( 5.16,2)+pow(0.64,2)+pow(0.71,2)));
    mHT_tau  ->SetBinContent(6, 26.23); mHT_tau  ->SetBinError(6,sqrt(pow( 7.60,2)+pow( 4.09,2)+pow(0.41,2)+pow(0.15,2)));
    mHT_tau  ->SetBinContent(7, 18.06); mHT_tau  ->SetBinError(7,sqrt(pow( 7.16,2)+pow( 2.95,2)+pow(0.23,2)+pow(0.23,2)));
    mHT_tau  ->SetBinContent(8, 27.39); mHT_tau  ->SetBinError(8,sqrt(pow(10.16,2)+pow( 4.24,2)+pow(0.26,2)+pow(0.48,2)));
    mHT_tau  ->SetBinContent(9, 36.27); mHT_tau  ->SetBinError(9,sqrt(pow(17.27,2)+pow( 9.11,2)+pow(0.54,2)+pow(0.41,2))); 
    mHT_tauMC->SetBinContent(1,107.15); mHT_tauMC->SetBinError(1, 3.05);
    mHT_tauMC->SetBinContent(2, 15.22); mHT_tauMC->SetBinError(2, 1.71);
    mHT_tauMC->SetBinContent(3,180.91); mHT_tauMC->SetBinError(3, 5.49);
    mHT_tauMC->SetBinContent(4, 72.15); mHT_tauMC->SetBinError(4, 6.12);
    mHT_tauMC->SetBinContent(5, 33.24); mHT_tauMC->SetBinError(5, 4.17);
    mHT_tauMC->SetBinContent(6, 19.87); mHT_tauMC->SetBinError(6, 1.72);
    mHT_tauMC->SetBinContent(7, 32.69); mHT_tauMC->SetBinError(7, 2.57);
    mHT_tauMC->SetBinContent(8, 32.17); mHT_tauMC->SetBinError(8, 2.56);
    mHT_tauMC->SetBinContent(9, 12.50); mHT_tauMC->SetBinError(9, 1.93);


    hHT_ele  ->SetBinContent(1, 12.08); hHT_ele  ->SetBinError(1,sqrt(pow( 2.81,2)+pow( 2.07,2)+pow(0.46,2)+pow(0.06,2)));
    hHT_ele  ->SetBinContent(2,  0.74); hHT_ele  ->SetBinError(2,sqrt(pow( 0.78,2)+pow( 0.24,2)+pow(0.07,2)+pow(0.03,2)));
    hHT_ele  ->SetBinContent(3, 16.34); hHT_ele  ->SetBinError(3,sqrt(pow( 3.54,2)+pow( 2.42,2)+pow(0.33,2)+pow(0.06,2)));
    hHT_ele  ->SetBinContent(4,  2.58); hHT_ele  ->SetBinError(4,sqrt(pow( 1.30,2)+pow( 0.58,2)+pow(0.07,2)+pow(0.12,2)));
    hHT_ele  ->SetBinContent(5,  4.08); hHT_ele  ->SetBinError(5,sqrt(pow( 2.35,2)+pow( 1.31,2)+pow(0.06,2)+pow(0.14,2)));
    hHT_ele  ->SetBinContent(6,  2.88); hHT_ele  ->SetBinError(6,sqrt(pow( 2.88,2)+pow( 0.86,2)+pow(0.14,2)+pow(0.25,2)));
    hHT_ele  ->SetBinContent(7,  2.97); hHT_ele  ->SetBinError(7,sqrt(pow( 2.10,2)+pow( 0.70,2)+pow(0.07,2)+pow(0.05,2)));
    hHT_ele  ->SetBinContent(8,  3.85); hHT_ele  ->SetBinError(8,sqrt(pow( 2.23,2)+pow( 1.01,2)+pow(0.00,2)+pow(0.17,2)));
//    hHT_ele  ->SetBinContent(9,); hHT_ele  ->SetBinError(9,sqrt(pow( ,2)+pow( ,2)+pow(,2)+pow(,2)));
    hHT_eleMC->SetBinContent(1, 13.82); hHT_eleMC->SetBinError(1, 1.06);
    hHT_eleMC->SetBinContent(2,  2.42); hHT_eleMC->SetBinError(2, 0.49);
    hHT_eleMC->SetBinContent(3, 21.06); hHT_eleMC->SetBinError(3, 1.37);
    hHT_eleMC->SetBinContent(4,  5.33); hHT_eleMC->SetBinError(4, 0.68);
    hHT_eleMC->SetBinContent(5,  3.37); hHT_eleMC->SetBinError(5, 0.65);
    hHT_eleMC->SetBinContent(6,  5.04); hHT_eleMC->SetBinError(6, 1.02);
    hHT_eleMC->SetBinContent(7,  4.39); hHT_eleMC->SetBinError(7, 0.64);
    hHT_eleMC->SetBinContent(8,  3.07); hHT_eleMC->SetBinError(8, 0.54);
    hHT_eleMC->SetBinContent(9,  0.83); hHT_eleMC->SetBinError(9, 0.28);

    hHT_muo  ->SetBinContent(1,  7.16); hHT_muo  ->SetBinError(1,sqrt(pow( 1.74,2)+pow( 1.46,2)+pow(0.33,2)+pow(0.08,2)));
    hHT_muo  ->SetBinContent(2,  2.20); hHT_muo  ->SetBinError(2,sqrt(pow( 1.10,2)+pow( 0.71,2)+pow(0.08,2)+pow(0.06,2)));
    hHT_muo  ->SetBinContent(3, 18.36); hHT_muo  ->SetBinError(3,sqrt(pow( 2.99,2)+pow( 3.30,2)+pow(0.42,2)+pow(0.05,2)));
    hHT_muo  ->SetBinContent(4,  7.97); hHT_muo  ->SetBinError(4,sqrt(pow( 2.22,2)+pow( 1.69,2)+pow(0.07,2)+pow(0.41,2)));
    hHT_muo  ->SetBinContent(5,  0.47); hHT_muo  ->SetBinError(5,sqrt(pow( 0.47,2)+pow( 0.17,2)+pow(0.00,2)+pow(0.27,2)));
    hHT_muo  ->SetBinContent(6,  1.95); hHT_muo  ->SetBinError(6,sqrt(pow( 1.41,2)+pow( 0.56,2)+pow(0.10,2)+pow(0.12,2)));
    hHT_muo  ->SetBinContent(7,  2.69); hHT_muo  ->SetBinError(7,sqrt(pow( 1.55,2)+pow( 0.76,2)+pow(0.05,2)+pow(0.12,2)));
    hHT_muo  ->SetBinContent(8,  2.50); hHT_muo  ->SetBinError(8,sqrt(pow( 1.44,2)+pow( 0.73,2)+pow(0.03,2)+pow(0.16,2)));
    //hHT_muo  ->SetBinContent(9,); hHT_muo  ->SetBinError(9,sqrt(pow( ,2)+pow( ,2)+pow(,2)+pow(,2)));
    hHT_muoMC->SetBinContent(1, 12.55); hHT_muoMC->SetBinError(1, 1.01);
    hHT_muoMC->SetBinContent(2,  1.92); hHT_muoMC->SetBinError(2, 0.46);
    hHT_muoMC->SetBinContent(3, 18.47); hHT_muoMC->SetBinError(3, 1.24);
    hHT_muoMC->SetBinContent(4,  6.71); hHT_muoMC->SetBinError(4, 1.30);
    hHT_muoMC->SetBinContent(5,  2.15); hHT_muoMC->SetBinError(5, 0.52);
    hHT_muoMC->SetBinContent(6,  2.42); hHT_muoMC->SetBinError(6, 0.48);
    hHT_muoMC->SetBinContent(7,  3.87); hHT_muoMC->SetBinError(7, 0.76);
    hHT_muoMC->SetBinContent(8,  3.23); hHT_muoMC->SetBinError(8, 0.71);
    hHT_muoMC->SetBinContent(9,  1.30); hHT_muoMC->SetBinError(9, 0.37);

    hHT_tau  ->SetBinContent(1, 10.97); hHT_tau  ->SetBinError(1,sqrt(pow( 9.17,2)+pow( 2.61,2)+pow(0.00,2)+pow(0.00,2)));
    hHT_tau  ->SetBinContent(2,  7.68); hHT_tau  ->SetBinError(2,sqrt(pow( 8.58,2)+pow( 4.78,2)+pow(0.13,2)+pow(0.04,2)));
    hHT_tau  ->SetBinContent(3, 22.64); hHT_tau  ->SetBinError(3,sqrt(pow( 8.42,2)+pow( 3.53,2)+pow(0.31,2)+pow(0.11,2)));
    hHT_tau  ->SetBinContent(4,  3.29); hHT_tau  ->SetBinError(4,sqrt(pow( 4.84,2)+pow( 0.89,2)+pow(0.02,2)+pow(0.00,2)));
  //  hHT_tau  ->SetBinContent(5,); hHT_tau  ->SetBinError(5,sqrt(pow( ,2)+pow( ,2)+pow(,2)+pow(,2)));
    hHT_tau  ->SetBinContent(6,  7.85); hHT_tau  ->SetBinError(6,sqrt(pow( 3.54,2)+pow( 2.16,2)+pow(0.11,2)+pow(0.11,2)));
 //   hHT_tau  ->SetBinContent(7,); hHT_tau  ->SetBinError(7,sqrt(pow( ,2)+pow( ,2)+pow(,2)+pow(,2)));
    hHT_tau  ->SetBinContent(8,  6.12); hHT_tau  ->SetBinError(8,sqrt(pow(7.77,2)+pow( 2.38,2)+pow(0.14,2)+pow(0.98,2)));
 //   hHT_tau  ->SetBinContent(9,); hHT_tau  ->SetBinError(9,sqrt(pow( ,2)+pow( ,2)+pow(,2)+pow(,2)));
    hHT_tauMC->SetBinContent(1, 12.48); hHT_tauMC->SetBinError(1, 0.99);
    hHT_tauMC->SetBinContent(2,  2.04); hHT_tauMC->SetBinError(2, 0.49);
    hHT_tauMC->SetBinContent(3, 23.54); hHT_tauMC->SetBinError(3, 1.49);
    hHT_tauMC->SetBinContent(4,  9.09); hHT_tauMC->SetBinError(4, 1.08);
    hHT_tauMC->SetBinContent(5,  3.03); hHT_tauMC->SetBinError(5, 0.57);
    hHT_tauMC->SetBinContent(6,  3.45); hHT_tauMC->SetBinError(6, 0.61);
    hHT_tauMC->SetBinContent(7,  5.93); hHT_tauMC->SetBinError(7, 0.77);
    hHT_tauMC->SetBinContent(8,  4.80); hHT_tauMC->SetBinError(8, 0.67);
    hHT_tauMC->SetBinContent(9,  2.09); hHT_tauMC->SetBinError(9, 0.49);


	//get normalization from background estimate as performed via LostLeptonEstimate.C and visualized with LLVisualization.C
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		string hs0, hs1, hs2, hs3; int numhistos;
		string hs = string("_") + samplekind + string("_") + HT_region[i3] + string("_") + signal_region[i2];
		hs0 = "_Nominal"+hs;
		int bin(0);
		if(signal_region[i2]=="2j0b")    bin = 1;
		if(signal_region[i2]=="2j1to2b") bin = 2;
		if(signal_region[i2]=="3to5j0b") bin = 3;
		if(signal_region[i2]=="3to5j1b") bin = 4;
		if(signal_region[i2]=="3to5j2b") bin = 5;
		if(signal_region[i2]=="6j0b")    bin = 6;
		if(signal_region[i2]=="6j1b")    bin = 7;
		if(signal_region[i2]=="6j2b")    bin = 8;
		if(signal_region[i2]=="3b")      bin = 9;
		double scale(1.), scaleerr(0.);
		double scaleMC(1.), scaleerrMC(0.);
		//load the corresponding yields and add them together (i.e. ele+muo+tau)
		if(HT_region[i3]=="lowHT"){
			double ele    = lHT_ele->GetBinContent(bin);
			double eleerr = lHT_ele->GetBinError(  bin);
			double muo    = lHT_muo->GetBinContent(bin);
			double muoerr = lHT_muo->GetBinError(  bin);
			double tau    = lHT_tau->GetBinContent(bin);
			double tauerr = lHT_tau->GetBinError(  bin);
			double MCele    = lHT_eleMC->GetBinContent(bin);
			double MCeleerr = lHT_eleMC->GetBinError(  bin);
			double MCmuo    = lHT_muoMC->GetBinContent(bin);
			double MCmuoerr = lHT_muoMC->GetBinError(  bin);
			double MCtau    = lHT_tauMC->GetBinContent(bin);
			double MCtauerr = lHT_tauMC->GetBinError(  bin);
			if(fixemptybin&&ele==0) { ele = MCele; eleerr = MCele;}
			if(fixemptybin&&muo==0) { muo = MCmuo; eleerr = MCmuo;}
			if(fixemptybin&&tau==0) { tau = MCtau; eleerr = MCtau;}
			scaleMC = MCele+MCmuo+MCtau;
			scaleerrMC = sqrt(MCeleerr*MCeleerr+MCmuoerr*MCmuoerr+MCtauerr*MCtauerr);
			scale = ele+muo+tau;
			scaleerr = sqrt(eleerr*eleerr+muoerr*muoerr+tauerr*tauerr);
			cout << "lowHT ele " << bin << ": " << ele << " +/- " << eleerr << endl;
			cout << "lowHT muo " << bin << ": " << muo << " +/- " << muoerr << endl;
			cout << "lowHT tau " << bin << ": " << tau << " +/- " << tauerr << endl;
		}
		if(HT_region[i3]=="mediumHT"){
			double ele    = mHT_ele->GetBinContent(bin);
			double eleerr = mHT_ele->GetBinError(  bin);
			double muo    = mHT_muo->GetBinContent(bin);
			double muoerr = mHT_muo->GetBinError(  bin);
			double tau    = mHT_tau->GetBinContent(bin);
			double tauerr = mHT_tau->GetBinError(  bin);
			double MCele    = mHT_eleMC->GetBinContent(bin);
			double MCeleerr = mHT_eleMC->GetBinError(  bin);
			double MCmuo    = mHT_muoMC->GetBinContent(bin);
			double MCmuoerr = mHT_muoMC->GetBinError(  bin);
			double MCtau    = mHT_tauMC->GetBinContent(bin);
			double MCtauerr = mHT_tauMC->GetBinError(  bin);
			if(fixemptybin&&ele==0) { ele = MCele; eleerr = MCele;}
			if(fixemptybin&&muo==0) { muo = MCmuo; eleerr = MCmuo;}
			if(fixemptybin&&tau==0) { tau = MCtau; eleerr = MCtau;}
			scaleMC = MCele+MCmuo+MCtau;
			scaleerrMC = sqrt(MCeleerr*MCeleerr+MCmuoerr*MCmuoerr+MCtauerr*MCtauerr);
			scale = ele+muo+tau;
			scaleerr = sqrt(eleerr*eleerr+muoerr*muoerr+tauerr*tauerr);
			cout << "mediumHT ele " << bin << ": " << ele << " +/- " << eleerr << endl;
			cout << "mediumHT muo " << bin << ": " << muo << " +/- " << muoerr << endl;
			cout << "mediumHT tau " << bin << ": " << tau << " +/- " << tauerr << endl;
		}
		if(HT_region[i3]=="highHT"){
			double ele    = hHT_ele->GetBinContent(bin);
			double eleerr = hHT_ele->GetBinError(  bin);
			double muo    = hHT_muo->GetBinContent(bin);
			double muoerr = hHT_muo->GetBinError(  bin);
			double tau    = hHT_tau->GetBinContent(bin);
			double tauerr = hHT_tau->GetBinError(  bin);
			double MCele    = hHT_eleMC->GetBinContent(bin);
			double MCeleerr = hHT_eleMC->GetBinError(  bin);
			double MCmuo    = hHT_muoMC->GetBinContent(bin);
			double MCmuoerr = hHT_muoMC->GetBinError(  bin);
			double MCtau    = hHT_tauMC->GetBinContent(bin);
			double MCtauerr = hHT_tauMC->GetBinError(  bin);
			if(fixemptybin&&ele==0) { ele = MCele; eleerr = MCele;}
			if(fixemptybin&&muo==0) { muo = MCmuo; eleerr = MCmuo;}
			if(fixemptybin&&tau==0) { tau = MCtau; eleerr = MCtau;}
			scaleMC = MCele+MCmuo+MCtau;
			scaleerrMC = sqrt(MCeleerr*MCeleerr+MCmuoerr*MCmuoerr+MCtauerr*MCtauerr);
			scale = ele+muo+tau;
			scaleerr = sqrt(eleerr*eleerr+muoerr*muoerr+tauerr*tauerr);
			cout << "highHT ele " << bin << ": " << ele << " +/- " << eleerr << endl;
			cout << "highHT muo " << bin << ": " << muo << " +/- " << muoerr << endl;
			cout << "highHT tau " << bin << ": " << tau << " +/- " << tauerr << endl;
		}

		//rescaled normalized shape you obtained from TTbarStudies.C
		double squaredsum = 0;
		for(int n = 1; n<= histosnormalized["NormMT2"+hs0]->GetNbinsX(); ++n){
			double content = histosnormalized["NormMT2"+hs0]->GetBinContent(n);
			histosnormalized["NormMT2"+hs0]->SetBinContent(n, content*scale);
			squaredsum += pow(content*scale,2);
			histosnormalized["NormMT2"+hs0]->SetBinError(  n, content*scaleerr);
		}
	}}


	//plot out something like a final table
	*fLogStream << "*********************************************************************" << endl;
	*fLogStream << "\%BEGINLATEX\%"             << endl;
	*fLogStream << "\\begin{table}"             << endl
	     << "\\begin{center}"            << endl
	     << "\\footnotesize"                    << endl;
       	*fLogStream << "\\begin{tabular}{|c|cc||c|cc||c|cc|}" << endl;	     
	*fLogStream << "\\hline\\hline"             << endl;
	*fLogStream << " $M_\\mathrm{T2}$ & \\multicolumn{2}{|c||}{low $H_\\mathrm{T}$} & $M_\\mathrm{T2}$ & \\multicolumn{2}{|c||}{medium $H_\\mathrm{T}$} & $M_\\mathrm{T2}$ & \\multicolumn{2}{|c|}{high $H_\\mathrm{T}$} \\\\" << endl;
	*fLogStream << "$\\text{[GeV]}$ & $\\mathrm{N}^\\mathrm{MC}_\\mathrm{truth}$ & $\\mathrm{N}^\\mathrm{data}_\\mathrm{prediction}$ & "
	            << "$\\text{[GeV]}$ & $\\mathrm{N}^\\mathrm{MC}_\\mathrm{truth}$ & $\\mathrm{N}^\\mathrm{data}_\\mathrm{prediction}$ & "
	            << "$\\text{[GeV]}$ & $\\mathrm{N}^\\mathrm{MC}_\\mathrm{truth}$ & $\\mathrm{N}^\\mathrm{data}_\\mathrm{prediction}$ "
	            << " \\\\" << endl;
	*fLogStream << "\\hline\\hline"             << endl;
	for(int i2 = 0; i2<signalregionsize; ++i2){
		string hsl = "_Nominal"+string("_") + samplekind + string("_") + HT_region[0] + string("_") + signal_region[i2];
		string hsm = "_Nominal"+string("_") + samplekind + string("_") + HT_region[1] + string("_") + signal_region[i2];
		string hsh = "_Nominal"+string("_") + samplekind + string("_") + HT_region[2] + string("_") + signal_region[i2];
		string bin;
		if(signal_region[i2]=="2j0b")    bin = "2 jets, 0 b-jets";
		if(signal_region[i2]=="2j1to2b") bin = "2 jets, $\\geq$1 b-jets";
		if(signal_region[i2]=="3to5j0b") bin = "3-5 jets, 0 b-jets";
		if(signal_region[i2]=="3to5j1b") bin = "3-5 jets, 1 b-jets";
		if(signal_region[i2]=="3to5j2b") bin = "3-5 jets, 2 b-jets";
		if(signal_region[i2]=="6j0b")    bin = "$\\geq$6 jets, 0 b-jets";
		if(signal_region[i2]=="6j1b")    bin = "$\\geq$6 jets, 1 b-jets";
		if(signal_region[i2]=="6j2b")    bin = "$\\geq$6 jets, 2 b-jets";
		if(signal_region[i2]=="3b")      bin = "$\\geq$3 jets, $\\geq$3 b-jets";
		//REMOVE THIS WHEN MOVING TO DIFFERENT SAMPLES
		if(signal_region[i2]=="2j1to2b"){//remove high weight single top - remove this is hardcoded bullshit, only valid for 8 TeV analysis
			histos["MT2"+hsm]->SetBinContent(2, histos["MT2"+hsm]->GetBinContent(2)-13.42);
			histos["MT2"+hsm]->SetBinError(2, sqrt(pow(histos["MT2"+hsm]->GetBinError(2),2)-pow(12.42,2) ) );
		}

		*fLogStream << " \\multicolumn{9}{|l|}{" << bin << "} \\\\" << endl;
		   *fLogStream << "\\hline" <<  fixed << setprecision(2) << endl;
		
		int nbinsl = histosnormalized["NormMT2"+hsl]->GetNbinsX();
		int nbinsm = histosnormalized["NormMT2"+hsm]->GetNbinsX();
		int nbinsh = histosnormalized["NormMT2"+hsh]->GetNbinsX();
		int nbins(0);
		nbins  = (nbinsl>nbinsm)?nbinsl:nbinsm;
		nbins  = (nbins>nbinsh)?nbins:nbinsh;
		for(int n = 1; n<=nbins; ++n){
			if(n>nbinsl) *fLogStream << " & & & ";
			else{
				int binlow   = histosnormalized["NormMT2"+hsl]->GetBinLowEdge(n);
				int binwidth = histosnormalized["NormMT2"+hsl]->GetBinWidth(n);
				double mctruth = histos["MT2"+hsl]->GetBinContent(n);
				double datapred = histosnormalized["NormMT2"+hsl]->GetBinContent(n);
				double dataprederr = histosnormalized["NormMT2"+hsl]->GetBinError(n);
				*fLogStream << binlow << "-";
				if(n==nbinsl) *fLogStream << "Inf";
				else *fLogStream << binlow+binwidth;
				*fLogStream << " & " << mctruth << " & " << datapred << "$\\pm$" << dataprederr << " & ";
			}
			if(n>nbinsm) *fLogStream << " & & & ";
			else{
				int binlow   = histosnormalized["NormMT2"+hsm]->GetBinLowEdge(n);
				int binwidth = histosnormalized["NormMT2"+hsm]->GetBinWidth(n);
				double mctruth = histos["MT2"+hsm]->GetBinContent(n);
				double datapred = histosnormalized["NormMT2"+hsm]->GetBinContent(n);
				double dataprederr = histosnormalized["NormMT2"+hsm]->GetBinError(n);
				*fLogStream << binlow << "-";
				if(n==nbinsm) *fLogStream << "Inf";
				else *fLogStream << binlow+binwidth;
				*fLogStream << " & " << mctruth << " & " << datapred << "$\\pm$" << dataprederr << " & ";
			}
			if(n>nbinsh) *fLogStream << " & &  \\\\" << endl;
			else{
				int binlow   = histosnormalized["NormMT2"+hsh]->GetBinLowEdge(n);
				int binwidth = histosnormalized["NormMT2"+hsh]->GetBinWidth(n);
				double mctruth = histos["MT2"+hsh]->GetBinContent(n);
				double datapred = histosnormalized["NormMT2"+hsh]->GetBinContent(n);
				double dataprederr = histosnormalized["NormMT2"+hsh]->GetBinError(n);
				*fLogStream << binlow << "-";
				if(n==nbinsh) *fLogStream << "Inf";
				else *fLogStream << binlow+binwidth;
				*fLogStream << " & " << mctruth << " & " << datapred << "$\\pm$" << dataprederr << " \\\\" << endl;
			}
		}
		*fLogStream << "\\hline" << endl;
	}
	*fLogStream << "\\hline" << endl;
	*fLogStream     << "\\end{tabular}"     << endl
	     << "\\end{center}"                 << endl
	     << "\\end{table}"                  << endl
	     << "\%ENDLATEX\%"                  << endl
	     << endl;

    	TFile *fsavefile = new TFile(inputdir + outputname,"RECREATE");
	fsavefile->cd();
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		string hs0, hs1, hs2, hs3; int numhistos;
		string hs = string("_") + samplekind + string("_") + HT_region[i3] + string("_") + signal_region[i2];
		hs0 = "_Nominal"+hs;
		histosnormalized["NormMT2"+hs0]->Write();
		histos["MT2"+hs0]->Write();
	}}
	lHT_ele->Write();
	lHT_muo->Write();
	lHT_tau->Write();
	mHT_ele->Write();
	mHT_muo->Write();
	mHT_tau->Write();
	hHT_ele->Write();
	hHT_muo->Write();
	hHT_tau->Write();
	lHT_eleMC->Write();
	lHT_muoMC->Write();
	lHT_tauMC->Write();
	mHT_eleMC->Write();
	mHT_muoMC->Write();
	mHT_tauMC->Write();
	hHT_eleMC->Write();
	hHT_muoMC->Write();
	hHT_tauMC->Write();
	fsavefile->Close();
	cout << "Saved histograms in " << inputdir << outputname << endl;

	if(fWriteToFile && fAppend){
		TString logname =inputdir + "LostLeptonEstimateTryoutSplittedbyMC.log"; 
		ofstream f_log (logname.Data(), ios::app);
		f_log << fLogStream->str();
		cout << "wrote results into  " << logname << " (appended at the end of old file)" << endl;
	}else if(fWriteToFile){
		TString logname =inputdir + "LostLeptonEstimateSplittedTryoutbyMC.log"; 
		ofstream f_log (logname.Data(), ios::trunc);
		f_log << fLogStream->str();
		cout << "wrote results into  " << logname <<  " (old file replaced)" << endl;
	} else{
		cout << fLogStream->str();
	}
	delete fLogStream;
}
