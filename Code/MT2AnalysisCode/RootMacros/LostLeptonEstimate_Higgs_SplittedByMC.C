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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path of MT2tree.hh
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/Utilities.hh"

//use via root -l LostLeptonEstimate_Higgs_SplittedByMC.C++

using namespace std;

void LostLeptonEstimate_Higgs_SplittedByMC();

const int HTregionsize = 2;
string HT_region[HTregionsize] = {"lowHT", "highHT"};

std::ostringstream* fLogStream     = 0;

//this code takes the lost lepton estimates produced with LostLeptonEstimate_Higgs.C (actually at the moment hard-coded)
//and the MT2 shapes produced with  TTbarStudiesHiggs.C (or TTbarStudieHiggsISR.C)
//and reweights the shapes to fit the normalization obtained by the prediction.
//This is then the final prediction of the LostLeptonEstimate (excluding shape uncertainties that come later into the game)
void LostLeptonEstimate_Higgs_SplittedByMC(){

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
	TString   inputdir = "Filtered/TTbarStudies/Higgs/";
	TString  inputname = "TTbarStudiesHistograms_all_ISR_limitMbb.root";
	if(WnoScaleMatching)  inputname = "TTbarStudiesHistograms_all_noWscaleupdown_ISR_limitMbb.root";
	TString  outputname = "LostLeptonEstimate_Higgs_SplittedByMC.root";
	if(WnoScaleMatching)  outputname = "LostLeptonEstimate_Higgs_Higgs_SplittedByMC_noWscaleupdown_ISR_limitMbb.root";


	TFile *infile = TFile::Open(inputdir + inputname);

	map<string, TH1D*>    histos;
	map<string, TH1D*>    histosnormalized;
	vector<string> histonamestosave; histonamestosave.clear();

	string samplekind;
	if(onlyttbar)    samplekind = "TTbar";
	else if(onlytop) samplekind = "allTop";
	else if(onlyW)   samplekind = "WJets";
	else             samplekind = "allMC";
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		string hs = string("_") + sample_type[i1] + string("_") + "Nominal" + string("_") + HT_region[i3];// + string("_") + signal_region[i2];
		string mapname = "MT2" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)infile->Get(mapname.c_str());
		string mapnamenorm = "Norm"+mapname;
		if(histosnormalized.count(mapnamenorm) == 0 ) histosnormalized[mapnamenorm] = (TH1D*)infile->Get(mapnamenorm.c_str());
	}

	//get the LostLepton estimate histograms for Higgs from LostLeptonEstimate_Higgs.C (at the moment this is hardcoded)
    TH1D *hLLest    = new TH1D("hLLest"   ,"",6,0,6); hLLest   ->SetMarkerStyle(20), hLLest   ->SetMarkerColor(kBlack); hLLest->SetLineWidth(3); hLLest->SetLineColor(kBlack);
    TH1D *hLLest_MC = new TH1D("hLLest_MC","",6,0,6); hLLest_MC->SetFillStyle(3001); hLLest_MC->SetFillColor(kBlue);

    hLLest_MC->SetBinContent(1,14.30); hLLest_MC->SetBinError(1,1.65);// e,   lHT
    hLLest   ->SetBinContent(1,15.54); hLLest   ->SetBinError(1,sqrt(pow(4.33,2)+pow(2.43,2)+pow(0.25,2)+pow(0.,2)));
    hLLest_MC->SetBinContent(2,11.15); hLLest_MC->SetBinError(2,1.41);// mu,  lHT
    hLLest   ->SetBinContent(2,12.49); hLLest   ->SetBinError(2,sqrt(pow(3.05,2)+pow(2.19,2)+pow(0.12,2)+pow(0.,2)));
    hLLest_MC->SetBinContent(3,18.31); hLLest_MC->SetBinError(3,1.81);// tau, lHT
//    hLLest   ->SetBinContent(3,6.83); hLLest   ->SetBinError(3,sqrt(pow(6.30,2)+pow(1.46,2)+pow(0.72,2)+pow(0.,2)));
    hLLest   ->SetBinContent(3,9.05); hLLest   ->SetBinError(3,sqrt(pow(6.26,2)+pow(1.60,2)+pow(0.11,2)+pow(0.,2)));
    hLLest_MC->SetBinContent(4,18.98); hLLest_MC->SetBinError(4,1.64);// e,   hHT
    hLLest   ->SetBinContent(4,11.86); hLLest   ->SetBinError(4,sqrt(pow(3.79,2)+pow(1.67,2)+pow(0.12,2)+pow(0.,2)));
    hLLest_MC->SetBinContent(5,17.96); hLLest_MC->SetBinError(5,2.44);// mu,  hHT
    hLLest   ->SetBinContent(5,22.08); hLLest   ->SetBinError(5,sqrt(pow(4.04,2)+pow(3.40,2)+pow(0.11,2)+pow(0.,2)));
    hLLest_MC->SetBinContent(6,29.98); hLLest_MC->SetBinError(6,3.78);// tau, hHT
    hLLest   ->SetBinContent(6,30.87); hLLest   ->SetBinError(6,sqrt(pow(14.04,2)+pow(5.24,2)+pow(0.32,2)+pow(0.,2)));

	for(int i3 = 0; i3<HTregionsize;     ++i3){
		string hs0;
		string hs = string("_") + samplekind + string("_") + HT_region[i3];// + string("_") + signal_region[i2];
		hs0 = "_Nominal"+hs;
		double scale(1.), scaleerr(0.);
		double scaleMC(1.), scaleerrMC(0.);
		if(HT_region[i3]=="lowHT"){
			double ele    = hLLest->GetBinContent(1);
			double eleerr = hLLest->GetBinError(  1);
			double muo    = hLLest->GetBinContent(2);
			double muoerr = hLLest->GetBinError(  2);
			double tau    = hLLest->GetBinContent(3);
			double tauerr = hLLest->GetBinError(  3);
			double MCele    = hLLest_MC->GetBinContent(1);
			double MCeleerr = hLLest_MC->GetBinError(  1);
			double MCmuo    = hLLest_MC->GetBinContent(2);
			double MCmuoerr = hLLest_MC->GetBinError(  2);
			double MCtau    = hLLest_MC->GetBinContent(3);
			double MCtauerr = hLLest_MC->GetBinError(  3);
			if(fixemptybin&&ele==0) { ele = MCele; eleerr = MCele;}
			if(fixemptybin&&muo==0) { muo = MCmuo; eleerr = MCmuo;}
			if(fixemptybin&&tau==0) { tau = MCtau; eleerr = MCtau;}
			scaleMC = MCele+MCmuo+MCtau;
			scaleerrMC = sqrt(MCeleerr*MCeleerr+MCmuoerr*MCmuoerr+MCtauerr*MCtauerr);
			scale = ele+muo+tau;
			scaleerr = sqrt(eleerr*eleerr+muoerr*muoerr+tauerr*tauerr);
		}
		if(HT_region[i3]=="highHT"){
			double ele    = hLLest->GetBinContent(4);
			double eleerr = hLLest->GetBinError(  4);
			double muo    = hLLest->GetBinContent(5);
			double muoerr = hLLest->GetBinError(  5);
			double tau    = hLLest->GetBinContent(6);
			double tauerr = hLLest->GetBinError(  6);
			double MCele    = hLLest_MC->GetBinContent(4);
			double MCeleerr = hLLest_MC->GetBinError(  4);
			double MCmuo    = hLLest_MC->GetBinContent(5);
			double MCmuoerr = hLLest_MC->GetBinError(  5);
			double MCtau    = hLLest_MC->GetBinContent(6);
			double MCtauerr = hLLest_MC->GetBinError(  6);
			if(fixemptybin&&ele==0) { ele = MCele; eleerr = MCele;}
			if(fixemptybin&&muo==0) { muo = MCmuo; eleerr = MCmuo;}
			if(fixemptybin&&tau==0) { tau = MCtau; eleerr = MCtau;}
			scaleMC = MCele+MCmuo+MCtau;
			scaleerrMC = sqrt(MCeleerr*MCeleerr+MCmuoerr*MCmuoerr+MCtauerr*MCtauerr);
			scale = ele+muo+tau;
			scaleerr = sqrt(eleerr*eleerr+muoerr*muoerr+tauerr*tauerr);
		}
		double squaredsum = 0;
			for(int n = 1; n<= histosnormalized["NormMT2"+hs0]->GetNbinsX(); ++n){
				double content = histosnormalized["NormMT2"+hs0]->GetBinContent(n);
				histosnormalized["NormMT2"+hs0]->SetBinContent(n, content*scale);
				histosnormalized["NormMT2"+hs0]->SetBinError(  n, content*scaleerr);
				bool existing = false;
				for(int m = 0; m<histonamestosave.size();++m){
					if(histonamestosave[m]=="NormMT2"+hs0) existing = true;
				}
				if(!existing) histonamestosave.push_back("NormMT2"+hs0);
			}
	}}
    	TFile *fsavefile = new TFile(inputdir + outputname,"RECREATE");
	fsavefile->cd();
	for(int i2 = 0; i2<histonamestosave.size(); ++i2){
		histosnormalized[histonamestosave[i2] ]->Write();
	}
	hLLest->Write();
	hLLest_MC->Write();
	cout << "files saved in " << fsavefile->GetName() << endl;
}
