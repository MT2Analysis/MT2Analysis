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
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TSystem.h"
#include "THStack.h"
#include "TMap.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TArrow.h"
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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"//use your owhn path


//run via root -l -b -q MT2Results_PlotsAndTables.C++

using namespace std;

//two functions for statistical poissonian uncertainties
//these should be obsolete in next ROOT version
//where TH1's can have poissonian uncertainties
float statErrorN(float x){return x -
0.5*TMath::ChisquareQuantile(0.3173/2.,2.*x);}
float statErrorP(float x){return
0.5*TMath::ChisquareQuantile(1.-0.3173/2.,2.*(x+1.))-x;}


//this function uses all background estimates (stored in root histograms)
//the data and possible also MC (also stored in root histograms)
//and produces final result tables (three versions) and plots (also three versions)
void MT2Results_PlotsAndTables(){


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


const int HTregionsize = 3;
string HT_bin[HTregionsize] = {"lowHT", "mediumHT", "highHT"};
const int signalregionsize = 9;
string signal_region[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b"};


	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	bool TOBTECreweight = true;//reweight MC due to TOBTEC filter efficiency - should be obsolete for next generation of analyses

	bool fixemptybin = true; //in estimates (Znunu and LostLepton) take MC yield +/-100percent if data prediction yields 0. default is true
	bool useMCyields = false;//this decides if in your result tables you also want to quote the MC truth numbers, default is false
	bool ISRusage    = true; //this decides if the MC shape for the LostLeptonEstimate are obtained after reweighting them due to the 'ISR recipe'. default is true

	bool saveeps = true; //save result plots as eps
	bool saveasC = true; //save result plots as .C root macros
	bool poisson = false;//plot data with poissonian statistics, default is false
	//works better for ROOT 5.34 where poissonian uncertainties work on TH1's, for now work-around with TGraphAsymmErrors implemented
	//note that data predictions (at least Z(nunu),LostLepton are done using gaussian statistics

	bool dummyy = false;           //this was a try to implement correctly (i.e. fully correlated) the 20percent systematics on Z/gamma ratio. Use at own risk, default is false
	bool summaryscale2bins = true; // in three summary plots we scale down the 2j0b and 3to5j0b yield and prediction so that the other bins are better visible, default is true

	//histograms containing data prediction, data, and MC yield
	TFile *ZnunuNumbers;
	if(TOBTECreweight) ZnunuNumbers = TFile::Open("../Results/Filtered/GammaJetsPrediction/20130617_test/ZnunuNumbers_TOBTECreweight.root");
	else               ZnunuNumbers = TFile::Open("../Results/Filtered/GammaJetsPrediction/20130617_test/ZnunuNumbers.root");
	TFile *ZinvPrediction           = TFile::Open("../Results/Filtered/GammaJetsPrediction/20130617_test/ZinvPredictionNumbers.root");
	TFile *LostLeptonPrediction;
	if(fixemptybin) LostLeptonPrediction         = TFile::Open("Filtered/TTbarStudies/LostLeptonEstimate_SplittedByMC_fixemptybin.root");
	else            LostLeptonPrediction         = TFile::Open("Filtered/TTbarStudies/LostLeptonEstimate_SplittedByMC.root");
	if(ISRusage){
		if(fixemptybin) LostLeptonPrediction = TFile::Open("Filtered/TTbarStudies/LostLeptonEstimate_SplittedByMC_fixemptybin_ISR.root");
		else            LostLeptonPrediction = TFile::Open("Filtered/TTbarStudies/LostLeptonEstimate_SplittedByMC_ISR.root");
	}
	TFile *QCDPredictionLow     = TFile::Open("~casal/QCDfromData_final.root");
	TFile *QCDPredictionMedium  = TFile::Open("~casal/qcd_filter/QCDfromData_mediumHT.root");
	TFile *QCDPredictionHigh    = TFile::Open("~casal/qcd_filter/QCDfromData_highHT.root");
	TFile *DataFile             = TFile::Open("../Results/Filtered/NewDataNumbers.root");
	TFile *MCFile;
	if(TOBTECreweight) MCFile   = TFile::Open("../Results/Filtered/MCNumbers_TOBTECreweight.root");
	else               MCFile   = TFile::Open("../Results/Filtered/MCNumbers.root");

	//file name for storing all elements used for final prediction
	TString newfilename                 = "../Results/Filtered/PredictionFile";
	if(fixemptybin)        newfilename += "_fixemptybin";
	if(ISRusage)           newfilename += "_ISR";
	if(TOBTECreweight)     newfilename += "_TOBTECreweight";
	if(poisson)            newfilename += "_PoissonErrors";
        TString allfile = newfilename + "_AllInOne.root";
	                       newfilename += ".root";
	//directory where result plots will be stored
	string outputdir                  = "../Results/Filtered/FinalPrediction/";
	if(fixemptybin)        outputdir += "fixemptybin/";
	if(dummyy)             outputdir = "dummy/";
	if(ISRusage)           outputdir += "ISR/";
	if(TOBTECreweight)     outputdir += "TOBTECreweight/";
	if(poisson)            outputdir += "PoissonErrors/";
    	Util::MakeOutputDir(outputdir);

	//histograms containing all signal bins, good for using in LEE studies
	TH1D *hZinvAll = new TH1D("hZinvAll", "", 123, 1., 124.); hZinvAll->Sumw2();
	TH1D *hLLAll   = new TH1D("hLLAll",   "", 123, 1., 124.); hLLAll  ->Sumw2();
	TH1D *hQCDAll  = new TH1D("hQCDAll",  "", 123, 1., 124.); hQCDAll ->Sumw2();
	TH1D *hPredAll = new TH1D("hPredAll", "", 123, 1., 124.); hPredAll->Sumw2();
	TH1D *hDataAll = new TH1D("hDataAll", "", 123, 1., 124.); hDataAll->Sumw2();
	//pull histograms
	TH1D *hPullData = new TH1D("hPullData", "", 50, -5., 5.); hPullData->Sumw2();
	TH1D *hPullMC   = new TH1D("hPullMC"  , "", 50, -5., 5.); hPullMC  ->Sumw2();
	TH1D *hPullMCDa = new TH1D("hPullMCDa", "", 50, -5., 5.); hPullMCDa->Sumw2();

	//all histograms
	map<string, TH1D*>    histos;
	map<string, TGraphAsymmErrors*>    graphData;//if next ROOT version has poissonian uncertainties for TH1's, graph will be obsolete
	map<string, TH1D*>    histosMC;
	map<string, THStack*> stackPrediction;

	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		int NMT2bins;
		if(i3==0){
			if(signal_region[i2]=="2j0b")    NMT2bins = gNMT2bins_2j0b_lHT;
			if(signal_region[i2]=="2j1to2b") NMT2bins = gNMT2bins_2j1b_lHT;
			if(signal_region[i2]=="3to5j0b") NMT2bins = gNMT2bins_3j0b_lHT;
			if(signal_region[i2]=="3to5j1b") NMT2bins = gNMT2bins_3j1b_lHT;
			if(signal_region[i2]=="3to5j2b") NMT2bins = gNMT2bins_3j2b_lHT;
			if(signal_region[i2]=="6j0b")    NMT2bins = gNMT2bins_6j0b_lHT;
			if(signal_region[i2]=="6j1b")    NMT2bins = gNMT2bins_6j1b_lHT;
			if(signal_region[i2]=="6j2b")    NMT2bins = gNMT2bins_6j2b_lHT;
			if(signal_region[i2]=="3b")      NMT2bins = gNMT2bins_3b_lHT;
		} if(i3==1){
			if(signal_region[i2]=="2j0b")    NMT2bins = gNMT2bins_2j0b_mHT;
			if(signal_region[i2]=="2j1to2b") NMT2bins = gNMT2bins_2j1b_mHT;
			if(signal_region[i2]=="3to5j0b") NMT2bins = gNMT2bins_3j0b_mHT;
			if(signal_region[i2]=="3to5j1b") NMT2bins = gNMT2bins_3j1b_mHT;
			if(signal_region[i2]=="3to5j2b") NMT2bins = gNMT2bins_3j2b_mHT;
			if(signal_region[i2]=="6j0b")    NMT2bins = gNMT2bins_6j0b_mHT;
			if(signal_region[i2]=="6j1b")    NMT2bins = gNMT2bins_6j1b_mHT;
			if(signal_region[i2]=="6j2b")    NMT2bins = gNMT2bins_6j2b_mHT;
			if(signal_region[i2]=="3b")      NMT2bins = gNMT2bins_3b_mHT;
		} if(i3==2){
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
  		double MT2bins[NMT2bins+1];
		if(i3==0){
			if(signal_region[i2]=="2j0b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_2j0b_lHT[i0]; }
			if(signal_region[i2]=="2j1to2b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_2j1b_lHT[i0]; }
			if(signal_region[i2]=="3to5j0b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j0b_lHT[i0]; }
			if(signal_region[i2]=="3to5j1b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j1b_lHT[i0]; }
			if(signal_region[i2]=="3to5j2b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j2b_lHT[i0]; }
			if(signal_region[i2]=="6j0b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j0b_lHT[i0]; }
			if(signal_region[i2]=="6j1b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j1b_lHT[i0]; }
			if(signal_region[i2]=="6j2b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j2b_lHT[i0]; }
			if(signal_region[i2]=="3b")      { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3b_lHT[i0];   }
		} if(i3==1){
			if(signal_region[i2]=="2j0b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_2j0b_mHT[i0];  }
			if(signal_region[i2]=="2j1to2b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_2j1b_mHT[i0];  }
			if(signal_region[i2]=="3to5j0b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j0b_mHT[i0];  }
			if(signal_region[i2]=="3to5j1b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j1b_mHT[i0];  }
			if(signal_region[i2]=="3to5j2b") { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3j2b_mHT[i0];  }
			if(signal_region[i2]=="6j0b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j0b_mHT[i0];  }
			if(signal_region[i2]=="6j1b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j1b_mHT[i0];  }
			if(signal_region[i2]=="6j2b")    { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_6j2b_mHT[i0];  }
			if(signal_region[i2]=="3b")      { for(int i0 = 0; i0<=NMT2bins; ++i0) MT2bins[i0] = gMT2bins_3b_mHT[i0];    }
		} if(i3==2){
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

		string hs = string("_") + HT_bin[i3] + string("_") + signal_region[i2];
		string loadstring = "MT2pred" + hs;
		string mapname = "MT2PredZinv" + hs;
		//Zinv
		if(i2!=4&&i2!=7&&i2!=8){ //ZinvPrediction
			ZinvPrediction->cd();
			TH1D *dummyZ = (TH1D*)ZinvPrediction->Get(loadstring.c_str());
			if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)dummyZ->Clone(mapname.c_str());
			histos[mapname]->Sumw2();
			histos[mapname]->SetFillColor(418);//419 is Z MC
			histos[mapname]->SetMarkerColor(418);
			histos[mapname]->SetLineWidth(1);
			histos[mapname]->SetLineColor(kBlack);
		} else { //ZnunuNumbers - needed for >=2b region
			loadstring = "MT2" + hs;
			ZnunuNumbers->cd();
			TH1D *dummyZ = (TH1D*)ZnunuNumbers->Get(loadstring.c_str());
			if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)dummyZ->Clone(mapname.c_str());
			for(int i = 1; i<=histos[mapname]->GetNbinsX(); ++i) {
				histos[mapname]->SetBinError(i, histos[mapname]->GetBinContent(i));
			}
			histos[mapname]->Sumw2();
			histos[mapname]->SetFillColor(418);//419 is Z MC
			histos[mapname]->SetMarkerColor(418);
			histos[mapname]->SetLineWidth(1);
			histos[mapname]->SetLineColor(kBlack);
		}
		//LostLepton
		loadstring = "NormMT2_Nominal"+string("_") + "allMC" + string("_") + HT_bin[i3] + string("_") + signal_region[i2];
		mapname = "MT2PredLostLepton" + hs;
		LostLeptonPrediction->cd();
		TH1D *dummyLL = (TH1D*)LostLeptonPrediction->Get(loadstring.c_str());
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)dummyLL->Clone(mapname.c_str());
		histos[mapname]->Sumw2();
		histos[mapname]->SetFillColor(430);//417 is W MC, 855 is Top MC
		histos[mapname]->SetMarkerColor(430);//this is kCyan-2
		histos[mapname]->SetLineWidth(1);
		histos[mapname]->SetLineColor(kBlack);
		//QCD
		loadstring = "QCDfromData" + hs;
		mapname = "MT2PredQCD" + hs;
		if(HT_bin[i3]=="lowHT"){
			QCDPredictionLow->cd();
			TH1D *dummyQCD = (TH1D*)QCDPredictionLow->Get(loadstring.c_str());
			if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)dummyQCD->Clone(mapname.c_str());
		}
		if(HT_bin[i3]=="mediumHT"){
			QCDPredictionMedium->cd();
			TH1D *dummyQCD = (TH1D*)QCDPredictionMedium->Get(loadstring.c_str());
			if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)dummyQCD->Clone(mapname.c_str());
		}
		if(HT_bin[i3]=="highHT"){
			QCDPredictionHigh->cd();
			TH1D *dummyQCD = (TH1D*)QCDPredictionHigh->Get(loadstring.c_str());
			if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)dummyQCD->Clone(mapname.c_str());
		}
		histos[mapname]->Sumw2();
		histos[mapname]->SetFillColor(402);//401 is QCD MC
		histos[mapname]->SetMarkerColor(402);
		histos[mapname]->SetLineWidth(1);
		histos[mapname]->SetLineColor(kBlack);
		for(int nb = 1; nb<=histos[mapname]->GetNbinsX(); ++nb){
			if(histos[mapname]->GetBinLowEdge(nb)>=200){//at high MT2>200 GeV, uncertainty is at least 50 percent
				if(histos[mapname]->GetBinError(nb)<0.5*histos[mapname]->GetBinContent(nb)){
					histos[mapname]->SetBinError(nb,0.5*histos[mapname]->GetBinContent(nb));
				}
			} else if((histos[mapname]->GetBinLowEdge(nb)+histos[mapname]->GetBinWidth(nb))>200){//at border of MT2 == 200 GeV, uncertainty is at least 30 percent - this is a dirty hack
				if(histos[mapname]->GetBinError(nb)<0.3*histos[mapname]->GetBinContent(nb)){
					histos[mapname]->SetBinError(nb,0.3*histos[mapname]->GetBinContent(nb));
				}
			}
		}
		//Data
		loadstring = "MT2" + hs;
		mapname = "MT2Data" + hs;
		DataFile->cd();
		TH1D *dummyData = (TH1D*)DataFile->Get(loadstring.c_str());
//		if(poisson) dummyData->SetBinErrorOption((TH1::EBinErrorOpt)1);//This function should work from next ROOT version on
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)dummyData->Clone(mapname.c_str());
//		if(poisson) histos[mapname]->SetBinErrorOption((TH1::EBinErrorOpt)1);//This function should work from next ROOT version on
		histos[mapname]->Sumw2();
		histos[mapname]->SetMarkerStyle(20);
		histos[mapname]->SetLineWidth(2);
		histos[mapname]->SetLineColor(kBlack);
		histos[mapname]->SetMarkerColor(kBlack);
		if(graphData.count(mapname) == 0 ){//if next ROOT version has poissonian uncertainties for TH1's, graph will be obsolete
			graphData[mapname] = new TGraphAsymmErrors();
			string graphname = "Graph"+mapname;
			graphData[mapname]->SetName(graphname.c_str());
			for(int ig = 0; ig<histos[mapname]->GetNbinsX(); ++ig){
				int nnn = graphData[mapname]->GetN();
				graphData[mapname]->SetPoint(nnn, histos[mapname]->GetBinCenter(ig+1), histos[mapname]->GetBinContent(ig+1));
				double pointerrorX = 0.5*histos[mapname]->GetBinWidth(ig+1);
				double pointerrorYlow = statErrorN(histos[mapname]->GetBinContent(ig+1));
				double pointerrorYup = statErrorP(histos[mapname]->GetBinContent(ig+1));
				graphData[mapname]->SetPointError(nnn, pointerrorX, pointerrorX, pointerrorYlow, pointerrorYup);
				graphData[mapname]->SetMarkerStyle(20);
				graphData[mapname]->SetLineWidth(2);
				graphData[mapname]->SetLineColor(kBlack);
				graphData[mapname]->SetMarkerColor(kBlack);
			}
			cout << endl;
		}
		//PredSum - this histogram will contain the prediction sum in order to plot the total uncertainty
		mapname = "MT2PredSum" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", NMT2bins, MT2bins);
		histos[mapname]->Sumw2();
		histos[mapname]->SetFillStyle(3013);
		histos[mapname]->SetFillColor(1);
		//THStack
		mapname = "MT2PredStack" + hs;
		if(stackPrediction.count(mapname)==0) stackPrediction[mapname] = new THStack(mapname.c_str(), "");
		//MC histograms
		if(useMCyields){
			string hsqcd = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + "_QCD";
			string hsZ   = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + "_ZJets";
			string hsW   = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + "_WJets";
			string hsTop = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + "_Top";
			string hsO   = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + "_Other";
			string hsMC  = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + "_mc";
			loadstring = "MT2" + hsqcd;
			mapname = "MT2MCQCD" + hs;
			TH1D *dummyMC1 = (TH1D*)MCFile->Get(loadstring.c_str());
			if(histosMC.count(mapname) == 0 ) histosMC[mapname] = (TH1D*)dummyMC1->Clone(mapname.c_str());
			loadstring = "MT2" + hsW;
			mapname = "MT2MCW" + hs;
			TH1D *dummyMC2 = (TH1D*)MCFile->Get(loadstring.c_str());
			if(histosMC.count(mapname) == 0 ) histosMC[mapname] = (TH1D*)dummyMC2->Clone(mapname.c_str());
			loadstring = "MT2" + hsTop;
			mapname = "MT2MCTop" + hs;
			TH1D *dummyMC3 = (TH1D*)MCFile->Get(loadstring.c_str());
			if(histosMC.count(mapname) == 0 ) histosMC[mapname] = (TH1D*)dummyMC3->Clone(mapname.c_str());
			loadstring = "MT2" + hsZ;
			mapname = "MT2MCZ" + hs;
			TH1D *dummyMC4 = (TH1D*)MCFile->Get(loadstring.c_str());
			if(histosMC.count(mapname) == 0 ) histosMC[mapname] = (TH1D*)dummyMC4->Clone(mapname.c_str());
			loadstring = "MT2" + hsO;
			mapname = "MT2MCOther" + hs;
			TH1D *dummyMC5 = (TH1D*)MCFile->Get(loadstring.c_str());
			if(histosMC.count(mapname) == 0 ) histosMC[mapname] = (TH1D*)dummyMC5->Clone(mapname.c_str());
			loadstring = "MT2" + hsMC;
			mapname = "MT2MCSum" + hs;
			TH1D *dummyMC6 = (TH1D*)MCFile->Get(loadstring.c_str());
			if(histosMC.count(mapname) == 0 ) histosMC[mapname] = (TH1D*)dummyMC6->Clone(mapname.c_str());
		}
	}}
	//fill stack with individual background estimates, add individual estimates to prediction sum
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		string hs = string("_") + HT_bin[i3] + string("_") + signal_region[i2];
		string mapnameZ     = "MT2PredZinv" + hs;
		string mapnameLL    = "MT2PredLostLepton" + hs;
		string mapnameQCD   = "MT2PredQCD" + hs;
		string mapnamePred  = "MT2PredSum" + hs;
		string mapnameData = "MT2Data" + hs;
		if(fixemptybin){
			for(int bn = 1; bn<=histos[mapnameZ]->GetNbinsX(); ++bn){
				if(histos[mapnameZ]->GetBinContent(bn)==0){
					string loadstring = "MT2" + hs;
					ZnunuNumbers->cd();
					TH1D *dummyZ = (TH1D*)ZnunuNumbers->Get(loadstring.c_str());
					histos[mapnameZ]->SetBinContent(bn, dummyZ->GetBinContent(bn));
					histos[mapnameZ]->SetBinError(bn, histos[mapnameZ]->GetBinContent(bn));
				}
			}
		}
		histos[mapnamePred]->Add(histos[mapnameQCD],1);
		histos[mapnamePred]->Add(histos[mapnameLL] ,1);
		histos[mapnamePred]->Add(histos[mapnameZ]  ,1);
		string mapnameStack = "MT2PredStack" + hs;
		stackPrediction[mapnameStack]->Add(histos[mapnameQCD]);
		stackPrediction[mapnameStack]->Add(histos[mapnameLL] );
		stackPrediction[mapnameStack]->Add(histos[mapnameZ]  );
		TString histoTitle = ";M_{T2} [GeV];Events";
		stackPrediction[mapnameStack]->SetTitle(histoTitle);

		//here we store all prediction into one big histogram - can be used if want to do LLE studies
		int binstart(0);
		if(i3==0){//lHT
		if(i2==0) binstart = 1; if(i2==1) binstart = 9; if(i2==2) binstart = 15; if(i2==3) binstart = 23; if(i2==4) binstart = 29; if(i2==5) binstart = 33; if(i2==6) binstart = 36; if(i2==7) binstart = 39; if(i2==8) binstart = 42;
		}
		if(i3==1){//mHT
		if(i2==0) binstart = 44; if(i2==1) binstart = 53; if(i2==2) binstart = 58; if(i2==3) binstart = 67; if(i2==4) binstart = 73; if(i2==5) binstart = 78; if(i2==6) binstart = 83; if(i2==7) binstart = 87; if(i2==8) binstart = 91;
		}
		if(i3==2){//hHT
		if(i2==0) binstart = 94; if(i2==1) binstart = 100; if(i2==2) binstart = 102; if(i2==3) binstart = 109; if(i2==4) binstart = 113; if(i2==5) binstart = 115; if(i2==6) binstart = 118; if(i2==7) binstart = 121; if(i2==8) binstart = 123;
		}
		for(int bn = 1; bn<=histos[mapnameZ]->GetNbinsX(); ++bn){
			hZinvAll->SetBinContent(binstart+bn-1, histos[mapnameZ]->GetBinContent(bn));
			hZinvAll->SetBinError(binstart+bn-1, histos[mapnameZ]->GetBinError(bn));
			hLLAll  ->SetBinContent(binstart+bn-1, histos[mapnameLL]->GetBinContent(bn));
			hLLAll  ->SetBinError(binstart+bn-1, histos[mapnameLL]->GetBinError(bn));
			hQCDAll ->SetBinContent(binstart+bn-1, histos[mapnameQCD]->GetBinContent(bn));
			hQCDAll ->SetBinError(binstart+bn-1, histos[mapnameQCD]->GetBinError(bn));
			hPredAll->SetBinContent(binstart+bn-1, histos[mapnamePred]->GetBinContent(bn));
			hPredAll->SetBinError(binstart+bn-1, histos[mapnamePred]->GetBinError(bn));
			hDataAll->SetBinContent(binstart+bn-1, histos[mapnameData]->GetBinContent(bn));
			hDataAll->SetBinError(binstart+bn-1, histos[mapnameData]->GetBinError(bn));
		}
	}}


	TLegend *leg = new TLegend(0.6652416,0.6311189,0.864906,0.9003497,NULL,"brNDC");
	leg->SetBorderSize(0);	leg->SetTextSize(0.04181185);	leg->SetTextFont(42);	leg->SetLineColor(1);
	leg->SetLineStyle(1);	leg->SetLineWidth(2);		leg->SetFillColor(0);	leg->SetFillStyle(1001);
	TLegend *legs = new TLegend(0.6992416,0.6311189,0.898906,0.9003497,NULL,"brNDC");//summary plot legend
	legs->SetBorderSize(0);	legs->SetTextSize(0.04181185);	legs->SetTextFont(42);	legs->SetLineColor(1);
	legs->SetLineStyle(1);	legs->SetLineWidth(2);		legs->SetFillColor(0);	legs->SetFillStyle(1001);
	string dummy = "MT2PredQCD_mediumHT_3to5j0b";   leg->AddEntry(histos[dummy], "Multijet", "f");             legs->AddEntry(histos[dummy], "Multijet", "f");
	dummy = "MT2PredLostLepton_lowHT_2j0b";         leg->AddEntry(histos[dummy], "Lost lepton", "f");          legs->AddEntry(histos[dummy], "Lost lepton", "f");
	dummy = "MT2PredZinv_lowHT_2j0b";               leg->AddEntry(histos[dummy], "Z(#nu#bar{#nu})+jets", "f"); legs->AddEntry(histos[dummy], "Z(#nu#bar{#nu})+jets", "f");
	dummy = "MT2Data_lowHT_2j0b";                   leg->AddEntry(histos[dummy], "Data", "lp");                legs->AddEntry(histos[dummy], "Data", "lp");
    
	//plots for three summary plots
	TH1D *predQCDlHT  = new TH1D("predQCDlHT", "",9,0,9);predQCDlHT ->SetFillColor(402);predQCDlHT ->SetMarkerColor(402);predQCDlHT ->SetLineWidth(1);predQCDlHT ->SetLineColor(1);
	TH1D *predQCDmHT  = new TH1D("predQCDmHT", "",9,0,9);predQCDmHT ->SetFillColor(402);predQCDmHT ->SetMarkerColor(402);predQCDmHT ->SetLineWidth(1);predQCDmHT ->SetLineColor(1);
	TH1D *predQCDhHT  = new TH1D("predQCDhHT", "",9,0,9);predQCDhHT ->SetFillColor(402);predQCDhHT ->SetMarkerColor(402);predQCDhHT ->SetLineWidth(1);predQCDhHT ->SetLineColor(1);
	TH1D *predLLlHT   = new TH1D("predLLlHT",  "",9,0,9);predLLlHT  ->SetFillColor(430);predLLlHT  ->SetMarkerColor(430);predLLlHT  ->SetLineWidth(1);predLLlHT  ->SetLineColor(1);
	TH1D *predLLmHT   = new TH1D("predLLmHT",  "",9,0,9);predLLmHT  ->SetFillColor(430);predLLmHT  ->SetMarkerColor(430);predLLmHT  ->SetLineWidth(1);predLLMHT  ->SetLineColor(1);
	TH1D *predLLhHT   = new TH1D("predLLhHT",  "",9,0,9);predLLhHT  ->SetFillColor(430);predLLhHT  ->SetMarkerColor(430);predLLhHT  ->SetLineWidth(1);predLLHHT  ->SetLineColor(1);
	TH1D *predZinvlHT = new TH1D("predZinvlHT","",9,0,9);predZinvlHT->SetFillColor(418);predZinvlHT->SetMarkerColor(418);predZinvlHT->SetLineWidth(1);predZinvlHT->SetLineColor(1);
	TH1D *predZinvmHT = new TH1D("predZinvmHT","",9,0,9);predZinvmHT->SetFillColor(418);predZinvmHT->SetMarkerColor(418);predZinvmHT->SetLineWidth(1);predZinvmHT->SetLineColor(1);
	TH1D *predZinvhHT = new TH1D("predZinvhHT","",9,0,9);predZinvhHT->SetFillColor(418);predZinvhHT->SetMarkerColor(418);predZinvhHT->SetLineWidth(1);predZinvhHT->SetLineColor(1);
	TH1D *predSumlHT  = new TH1D("predSumlHT", "",9,0,9);predSumlHT ->SetFillColor(  1);predSumlHT ->SetFillStyle(3013);
	TH1D *predSummHT  = new TH1D("predSummHT", "",9,0,9);predSummHT ->SetFillColor(  1);predSummHT ->SetFillStyle(3013);
	TH1D *predSumhHT  = new TH1D("predSumhHT", "",9,0,9);predSumhHT ->SetFillColor(  1);predSumhHT ->SetFillStyle(3013);
	TH1D *predDatalHT = new TH1D("predDatalHT","",9,0,9);predDatalHT->SetLineWidth(  2);predDatalHT->SetMarkerStyle(20);predDatalHT->SetLineColor(1);
	TH1D *predDatamHT = new TH1D("predDatamHT","",9,0,9);predDatamHT->SetLineWidth(  2);predDatamHT->SetMarkerStyle(20);predDatamHT->SetLineColor(1);
	TH1D *predDatahHT = new TH1D("predDatahHT","",9,0,9);predDatahHT->SetLineWidth(  2);predDatahHT->SetMarkerStyle(20);predDatahHT->SetLineColor(1);
	THStack *SumlHT = new THStack("SumlHT",""); THStack *SummHT = new THStack("SummHT",""); THStack *SumhHT = new THStack("SumhHT","");

	predQCDlHT ->GetXaxis()->SetBinLabel(1, "2j,0b");         predQCDmHT ->GetXaxis()->SetBinLabel(1, "2j,0b");         predQCDhHT ->GetXaxis()->SetBinLabel(1, "2j,0b");
	predQCDlHT ->GetXaxis()->SetBinLabel(2, "2j,#geq1b");     predQCDmHT ->GetXaxis()->SetBinLabel(2, "2j,#geq1b");     predQCDhHT ->GetXaxis()->SetBinLabel(2, "2j,#geq1b");
	predQCDlHT ->GetXaxis()->SetBinLabel(3, "3-5j,0b");       predQCDmHT ->GetXaxis()->SetBinLabel(3, "3-5j,0b");       predQCDhHT ->GetXaxis()->SetBinLabel(3, "3-5j,0b");
	predQCDlHT ->GetXaxis()->SetBinLabel(4, "3-5j,1b");       predQCDmHT ->GetXaxis()->SetBinLabel(4, "3-5j,1b");       predQCDhHT ->GetXaxis()->SetBinLabel(4, "3-5j,1b");
	predQCDlHT ->GetXaxis()->SetBinLabel(5, "3-5j,2b");       predQCDmHT ->GetXaxis()->SetBinLabel(5, "3-5j,2b");       predQCDhHT ->GetXaxis()->SetBinLabel(5, "3-5j,2b");
	predQCDlHT ->GetXaxis()->SetBinLabel(6, "#geq6j,0b");     predQCDmHT ->GetXaxis()->SetBinLabel(6, "#geq6j,0b");     predQCDhHT ->GetXaxis()->SetBinLabel(6, "#geq6j,0b");
	predQCDlHT ->GetXaxis()->SetBinLabel(7, "#geq6j,1b");     predQCDmHT ->GetXaxis()->SetBinLabel(7, "#geq6j,1b");     predQCDhHT ->GetXaxis()->SetBinLabel(7, "#geq6j,1b");
	predQCDlHT ->GetXaxis()->SetBinLabel(8, "#geq6j,2b");     predQCDmHT ->GetXaxis()->SetBinLabel(8, "#geq6j,2b");     predQCDhHT ->GetXaxis()->SetBinLabel(8, "#geq6j,2b");
	predQCDlHT ->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b"); predQCDmHT ->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b"); predQCDhHT ->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b");
	predLLlHT  ->GetXaxis()->SetBinLabel(1, "2j,0b");         predLLmHT  ->GetXaxis()->SetBinLabel(1, "2j,0b");         predLLhHT  ->GetXaxis()->SetBinLabel(1, "2j,0b");
	predLLlHT  ->GetXaxis()->SetBinLabel(2, "2j,#geq1b");     predLLmHT  ->GetXaxis()->SetBinLabel(2, "2j,#geq1b");     predLLhHT  ->GetXaxis()->SetBinLabel(2, "2j,#geq1b");
	predLLlHT  ->GetXaxis()->SetBinLabel(3, "3-5j,0b");       predLLmHT  ->GetXaxis()->SetBinLabel(3, "3-5j,0b");       predLLhHT  ->GetXaxis()->SetBinLabel(3, "3-5j,0b");
	predLLlHT  ->GetXaxis()->SetBinLabel(4, "3-5j,1b");       predLLmHT  ->GetXaxis()->SetBinLabel(4, "3-5j,1b");       predLLhHT  ->GetXaxis()->SetBinLabel(4, "3-5j,1b");
	predLLlHT  ->GetXaxis()->SetBinLabel(5, "3-5j,2b");       predLLmHT  ->GetXaxis()->SetBinLabel(5, "3-5j,2b");       predLLhHT  ->GetXaxis()->SetBinLabel(5, "3-5j,2b");
	predLLlHT  ->GetXaxis()->SetBinLabel(6, "#geq6j,0b");     predLLmHT  ->GetXaxis()->SetBinLabel(6, "#geq6j,0b");     predLLhHT  ->GetXaxis()->SetBinLabel(6, "#geq6j,0b");
	predLLlHT  ->GetXaxis()->SetBinLabel(7, "#geq6j,1b");     predLLmHT  ->GetXaxis()->SetBinLabel(7, "#geq6j,1b");     predLLhHT  ->GetXaxis()->SetBinLabel(7, "#geq6j,1b");
	predLLlHT  ->GetXaxis()->SetBinLabel(8, "#geq6j,2b");     predLLmHT  ->GetXaxis()->SetBinLabel(8, "#geq6j,2b");     predLLhHT  ->GetXaxis()->SetBinLabel(8, "#geq6j,2b");
	predLLlHT  ->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b"); predLLmHT  ->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b"); predLLhHT  ->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b");
	predZinvlHT->GetXaxis()->SetBinLabel(1, "2j,0b");         predZinvmHT->GetXaxis()->SetBinLabel(1, "2j,0b");         predZinvhHT->GetXaxis()->SetBinLabel(1, "2j,0b");
	predZinvlHT->GetXaxis()->SetBinLabel(2, "2j,#geq1b");     predZinvmHT->GetXaxis()->SetBinLabel(2, "2j,#geq1b");     predZinvhHT->GetXaxis()->SetBinLabel(2, "2j,#geq1b");
	predZinvlHT->GetXaxis()->SetBinLabel(3, "3-5j,0b");       predZinvmHT->GetXaxis()->SetBinLabel(3, "3-5j,0b");       predZinvhHT->GetXaxis()->SetBinLabel(3, "3-5j,0b");
	predZinvlHT->GetXaxis()->SetBinLabel(4, "3-5j,1b");       predZinvmHT->GetXaxis()->SetBinLabel(4, "3-5j,1b");       predZinvhHT->GetXaxis()->SetBinLabel(5, "3-5j,1b");
	predZinvlHT->GetXaxis()->SetBinLabel(5, "3-5j,2b");       predZinvmHT->GetXaxis()->SetBinLabel(5, "3-5j,2b");       predZinvhHT->GetXaxis()->SetBinLabel(5, "3-5j,2b");
	predZinvlHT->GetXaxis()->SetBinLabel(6, "#geq6j,0b");     predZinvmHT->GetXaxis()->SetBinLabel(6, "#geq6j,0b");     predZinvhHT->GetXaxis()->SetBinLabel(6, "#geq6j,0b");
	predZinvlHT->GetXaxis()->SetBinLabel(7, "#geq6j,1b");     predZinvmHT->GetXaxis()->SetBinLabel(7, "#geq6j,1b");     predZinvhHT->GetXaxis()->SetBinLabel(7, "#geq6j,1b");
	predZinvlHT->GetXaxis()->SetBinLabel(8, "#geq6j,2b");     predZinvmHT->GetXaxis()->SetBinLabel(8, "#geq6j,2b");     predZinvhHT->GetXaxis()->SetBinLabel(8, "#geq6j,2b");
	predZinvlHT->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b"); predZinvmHT->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b"); predZinvhHT->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b");
	predDatalHT->GetXaxis()->SetBinLabel(1, "2j,0b");         predDatamHT->GetXaxis()->SetBinLabel(1, "2j,0b");         predDatahHT->GetXaxis()->SetBinLabel(1, "2j,0b");
	predDatalHT->GetXaxis()->SetBinLabel(2, "2j,#geq1b");     predDatamHT->GetXaxis()->SetBinLabel(2, "2j,#geq1b");     predDatahHT->GetXaxis()->SetBinLabel(2, "2j,#geq1b");
	predDatalHT->GetXaxis()->SetBinLabel(3, "3-5j,0b");       predDatamHT->GetXaxis()->SetBinLabel(3, "3-5j,0b");       predDatahHT->GetXaxis()->SetBinLabel(3, "3-5j,0b");
	predDatalHT->GetXaxis()->SetBinLabel(4, "3-5j,1b");       predDatamHT->GetXaxis()->SetBinLabel(4, "3-5j,1b");       predDatahHT->GetXaxis()->SetBinLabel(4, "3-5j,1b");
	predDatalHT->GetXaxis()->SetBinLabel(5, "3-5j,2b");       predDatamHT->GetXaxis()->SetBinLabel(5, "3-5j,2b");       predDatahHT->GetXaxis()->SetBinLabel(5, "3-5j,2b");
	predDatalHT->GetXaxis()->SetBinLabel(6, "#geq6j,0b");     predDatamHT->GetXaxis()->SetBinLabel(6, "#geq6j,0b");     predDatahHT->GetXaxis()->SetBinLabel(6, "#geq6j,0b");
	predDatalHT->GetXaxis()->SetBinLabel(7, "#geq6j,1b");     predDatamHT->GetXaxis()->SetBinLabel(7, "#geq6j,1b");     predDatahHT->GetXaxis()->SetBinLabel(7, "#geq6j,1b");
	predDatalHT->GetXaxis()->SetBinLabel(8, "#geq6j,2b");     predDatamHT->GetXaxis()->SetBinLabel(8, "#geq6j,2b");     predDatahHT->GetXaxis()->SetBinLabel(8, "#geq6j,2b");
	predDatalHT->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b"); predDatamHT->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b"); predDatahHT->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b");
	predSumlHT ->GetXaxis()->SetBinLabel(1, "2j,0b");         predSummHT ->GetXaxis()->SetBinLabel(1, "2j,0b");         predSumhHT ->GetXaxis()->SetBinLabel(1, "2j,0b");
	predSumlHT ->GetXaxis()->SetBinLabel(2, "2j,#geq1b");     predSummHT ->GetXaxis()->SetBinLabel(2, "2j,#geq1b");     predSumhHT ->GetXaxis()->SetBinLabel(2, "2j,#geq1b");
	predSumlHT ->GetXaxis()->SetBinLabel(3, "3-5j,0b");       predSummHT ->GetXaxis()->SetBinLabel(3, "3-5j,0b");       predSumhHT ->GetXaxis()->SetBinLabel(3, "3-5j,0b");
	predSumlHT ->GetXaxis()->SetBinLabel(4, "3-5j,1b");       predSummHT ->GetXaxis()->SetBinLabel(4, "3-5j,1b");       predSumhHT ->GetXaxis()->SetBinLabel(4, "3-5j,1b");
	predSumlHT ->GetXaxis()->SetBinLabel(5, "3-5j,2b");       predSummHT ->GetXaxis()->SetBinLabel(5, "3-5j,2b");       predSumhHT ->GetXaxis()->SetBinLabel(5, "3-5j,2b");
	predSumlHT ->GetXaxis()->SetBinLabel(6, "#geq6j,0b");     predSummHT ->GetXaxis()->SetBinLabel(6, "#geq6j,0b");     predSumhHT ->GetXaxis()->SetBinLabel(6, "#geq6j,0b");
	predSumlHT ->GetXaxis()->SetBinLabel(7, "#geq6j,1b");     predSummHT ->GetXaxis()->SetBinLabel(7, "#geq6j,1b");     predSumhHT ->GetXaxis()->SetBinLabel(7, "#geq6j,1b");
	predSumlHT ->GetXaxis()->SetBinLabel(8, "#geq6j,2b");     predSummHT ->GetXaxis()->SetBinLabel(8, "#geq6j,2b");     predSumhHT ->GetXaxis()->SetBinLabel(8, "#geq6j,2b");
	predSumlHT ->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b"); predSummHT ->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b"); predSumhHT ->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b");

    
	//plots for single summary plot
    	TH1D *predQCD  = new TH1D("predQCD", "",27,0,27);predQCD ->SetFillColor(402);predQCD ->SetMarkerColor(402);predQCD ->SetLineWidth(2);predQCD ->SetLineColor(kBlack);
	TH1D *predLL   = new TH1D("predLL",  "",27,0,27);predLL  ->SetFillColor(430);predLL  ->SetMarkerColor(430);predLL  ->SetLineWidth(2);predLL ->SetLineColor(kBlack);
	TH1D *predZinv = new TH1D("predZinv","",27,0,27);predZinv->SetFillColor(418);predZinv->SetMarkerColor(418);predZinv->SetLineWidth(2);predZinv ->SetLineColor(kBlack);
	TH1D *predSum  = new TH1D("predSum", "",27,0,27);predSum ->SetFillColor(  1);predSum ->SetFillStyle(3013); predSum->Sumw2();
	TH1D *predData = new TH1D("predData","",27,0,27);predData->SetLineWidth(  2);predData->SetMarkerStyle(20); predData->Sumw2(); predData->SetLineColor(  kBlack);predData->SetMarkerColor(kBlack);
	//TH1D *predRatio = new TH1D("predRatio","",27,0,27);predRatio->SetLineWidth(  2);predRatio->SetMarkerStyle(20); predRatio->Sumw2();
	//TH1D *predRatioR = new TH1D("predRatioR","",27,0,27);predRatioR->SetFillColor(  1);predRatioR ->SetFillStyle(3013); predRatioR->Sumw2();
	THStack *Sum = new THStack("Sum","");
    	predQCD->SetMinimum(2.5);    predLL->SetMinimum(2.5);    predZinv->SetMinimum(2.5);    predSum->SetMinimum(2.5);    predData->SetMinimum(2.5);

	predQCD ->GetXaxis()->SetBinLabel(1, "2j,0b");      predQCD ->GetXaxis()->SetBinLabel(2, "2j,#geq1b");    	predQCD ->GetXaxis()->SetBinLabel(3, "3-5j,0b");      
	predQCD ->GetXaxis()->SetBinLabel(4, "3-5j,1b");    predQCD ->GetXaxis()->SetBinLabel(5, "3-5j,2b");      	predQCD ->GetXaxis()->SetBinLabel(6, "#geq6j,0b");    
	predQCD ->GetXaxis()->SetBinLabel(7, "#geq6j,1b");  predQCD ->GetXaxis()->SetBinLabel(8, "#geq6j,2b");    	predQCD ->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b");
	predQCD ->GetXaxis()->SetBinLabel(10, "2j,0b");     predQCD ->GetXaxis()->SetBinLabel(11, "2j,#geq1b");    	predQCD ->GetXaxis()->SetBinLabel(12, "3-5j,0b");      
	predQCD ->GetXaxis()->SetBinLabel(13, "3-5j,1b");   predQCD ->GetXaxis()->SetBinLabel(14, "3-5j,2b");      	predQCD ->GetXaxis()->SetBinLabel(15, "#geq6j,0b");    
	predQCD ->GetXaxis()->SetBinLabel(16, "#geq6j,1b"); predQCD ->GetXaxis()->SetBinLabel(17, "#geq6j,2b");    	predQCD ->GetXaxis()->SetBinLabel(18, "#geq3j,#geq3b");
	predQCD ->GetXaxis()->SetBinLabel(19, "2j,0b");	    predQCD ->GetXaxis()->SetBinLabel(20, "2j,#geq1b");		predQCD ->GetXaxis()->SetBinLabel(21, "3-5j,0b");
	predQCD ->GetXaxis()->SetBinLabel(22, "3-5j,1b");   predQCD ->GetXaxis()->SetBinLabel(23, "3-5j,2b");		predQCD ->GetXaxis()->SetBinLabel(24, "#geq6j,0b");
	predQCD ->GetXaxis()->SetBinLabel(25, "#geq6j,1b"); predQCD ->GetXaxis()->SetBinLabel(26, "#geq6j,2b");		predQCD ->GetXaxis()->SetBinLabel(27, "#geq3j,#geq3b");
	predLL  ->GetXaxis()->SetBinLabel(1, "2j,0b");      predLL  ->GetXaxis()->SetBinLabel(2, "2j,#geq1b");    	predLL  ->GetXaxis()->SetBinLabel(3, "3-5j,0b");      
	predLL  ->GetXaxis()->SetBinLabel(4, "3-5j,1b");    predLL  ->GetXaxis()->SetBinLabel(5, "3-5j,2b");      	predLL  ->GetXaxis()->SetBinLabel(6, "#geq6j,0b");    
	predLL  ->GetXaxis()->SetBinLabel(7, "#geq6j,1b");  predLL  ->GetXaxis()->SetBinLabel(8, "#geq6j,2b");    	predLL  ->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b");
	predLL  ->GetXaxis()->SetBinLabel(10, "2j,0b");     predLL  ->GetXaxis()->SetBinLabel(11, "2j,#geq1b");    	predLL  ->GetXaxis()->SetBinLabel(12, "3-5j,0b");      
	predLL  ->GetXaxis()->SetBinLabel(13, "3-5j,1b");   predLL  ->GetXaxis()->SetBinLabel(14, "3-5j,2b");      	predLL  ->GetXaxis()->SetBinLabel(15, "#geq6j,0b");    
	predLL  ->GetXaxis()->SetBinLabel(16, "#geq6j,1b"); predLL  ->GetXaxis()->SetBinLabel(17, "#geq6j,2b");    	predLL  ->GetXaxis()->SetBinLabel(18, "#geq3j,#geq3b");
	predLL  ->GetXaxis()->SetBinLabel(19, "2j,0b");     predLL  ->GetXaxis()->SetBinLabel(20, "2j,#geq1b");		predLL  ->GetXaxis()->SetBinLabel(21, "3-5j,0b");
	predLL  ->GetXaxis()->SetBinLabel(22, "3-5j,1b");   predLL  ->GetXaxis()->SetBinLabel(23, "3-5j,2b");		predLL  ->GetXaxis()->SetBinLabel(24, "#geq6j,0b");
	predLL  ->GetXaxis()->SetBinLabel(25, "#geq6j,1b"); predLL  ->GetXaxis()->SetBinLabel(26, "#geq6j,2b");		predLL  ->GetXaxis()->SetBinLabel(27, "#geq3j,#geq3b");
	predZinv->GetXaxis()->SetBinLabel(1, "2j,0b");      predZinv->GetXaxis()->SetBinLabel(2, "2j,#geq1b");    	predZinv->GetXaxis()->SetBinLabel(3, "3-5j,0b");      
	predZinv->GetXaxis()->SetBinLabel(4, "3-5j,1b");    predZinv->GetXaxis()->SetBinLabel(5, "3-5j,2b");      	predZinv->GetXaxis()->SetBinLabel(6, "#geq6j,0b");    
	predZinv->GetXaxis()->SetBinLabel(7, "#geq6j,1b");  predZinv->GetXaxis()->SetBinLabel(8, "#geq6j,2b");    	predZinv->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b");
	predZinv->GetXaxis()->SetBinLabel(10, "2j,0b");     predZinv->GetXaxis()->SetBinLabel(11, "2j,#geq1b");    	predZinv->GetXaxis()->SetBinLabel(12, "3-5j,0b");      
	predZinv->GetXaxis()->SetBinLabel(13, "3-5j,1b");   predZinv->GetXaxis()->SetBinLabel(14, "3-5j,2b");      	predZinv->GetXaxis()->SetBinLabel(15, "#geq6j,0b");    
	predZinv->GetXaxis()->SetBinLabel(16, "#geq6j,1b"); predZinv->GetXaxis()->SetBinLabel(17, "#geq6j,2b");    	predZinv->GetXaxis()->SetBinLabel(18, "#geq3j,#geq3b");
	predZinv->GetXaxis()->SetBinLabel(19, "2j,0b");	    predZinv->GetXaxis()->SetBinLabel(20, "2j,#geq1b");		predZinv->GetXaxis()->SetBinLabel(21, "3-5j,0b");
	predZinv->GetXaxis()->SetBinLabel(22, "3-5j,1b");   predZinv->GetXaxis()->SetBinLabel(23, "3-5j,2b");		predZinv->GetXaxis()->SetBinLabel(24, "#geq6j,0b");
	predZinv->GetXaxis()->SetBinLabel(25, "#geq6j,1b"); predZinv->GetXaxis()->SetBinLabel(26, "#geq6j,2b");		predZinv->GetXaxis()->SetBinLabel(27, "#geq3j,#geq3b");
	predData->GetXaxis()->SetBinLabel(1, "2j,0b");      predData->GetXaxis()->SetBinLabel(2, "2j,#geq1b");    	predData->GetXaxis()->SetBinLabel(3, "3-5j,0b");      
	predData->GetXaxis()->SetBinLabel(4, "3-5j,1b");    predData->GetXaxis()->SetBinLabel(5, "3-5j,2b");      	predData->GetXaxis()->SetBinLabel(6, "#geq6j,0b");    
	predData->GetXaxis()->SetBinLabel(7, "#geq6j,1b");  predData->GetXaxis()->SetBinLabel(8, "#geq6j,2b");    	predData->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b");
	predData->GetXaxis()->SetBinLabel(10, "2j,0b");     predData->GetXaxis()->SetBinLabel(11, "2j,#geq1b");    	predData->GetXaxis()->SetBinLabel(12, "3-5j,0b");      
	predData->GetXaxis()->SetBinLabel(13, "3-5j,1b");   predData->GetXaxis()->SetBinLabel(14, "3-5j,2b");      	predData->GetXaxis()->SetBinLabel(15, "#geq6j,0b");    
	predData->GetXaxis()->SetBinLabel(16, "#geq6j,1b"); predData->GetXaxis()->SetBinLabel(17, "#geq6j,2b");    	predData->GetXaxis()->SetBinLabel(18, "#geq3j,#geq3b");
	predData->GetXaxis()->SetBinLabel(19, "2j,0b");	    predData->GetXaxis()->SetBinLabel(20, "2j,#geq1b");		predData->GetXaxis()->SetBinLabel(21, "3-5j,0b");
	predData->GetXaxis()->SetBinLabel(22, "3-5j,1b");   predData->GetXaxis()->SetBinLabel(23, "3-5j,2b");		predData->GetXaxis()->SetBinLabel(24, "#geq6j,0b");
	predData->GetXaxis()->SetBinLabel(25, "#geq6j,1b"); predData->GetXaxis()->SetBinLabel(26, "#geq6j,2b");		predData->GetXaxis()->SetBinLabel(27, "#geq3j,#geq3b");
	predSum ->GetXaxis()->SetBinLabel(1, "2j,0b");      predSum ->GetXaxis()->SetBinLabel(2, "2j,#geq1b");    	predSum ->GetXaxis()->SetBinLabel(3, "3-5j,0b");      
	predSum ->GetXaxis()->SetBinLabel(4, "3-5j,1b");    predSum ->GetXaxis()->SetBinLabel(5, "3-5j,2b");      	predSum ->GetXaxis()->SetBinLabel(6, "#geq6j,0b");    
	predSum ->GetXaxis()->SetBinLabel(7, "#geq6j,1b");  predSum ->GetXaxis()->SetBinLabel(8, "#geq6j,2b");    	predSum ->GetXaxis()->SetBinLabel(9, "#geq3j,#geq3b");
	predSum ->GetXaxis()->SetBinLabel(10, "2j,0b");     predSum ->GetXaxis()->SetBinLabel(11, "2j,#geq1b");    	predSum ->GetXaxis()->SetBinLabel(12, "3-5j,0b");      
	predSum ->GetXaxis()->SetBinLabel(13, "3-5j,1b");   predSum ->GetXaxis()->SetBinLabel(14, "3-5j,2b");      	predSum ->GetXaxis()->SetBinLabel(15, "#geq6j,0b");    
	predSum ->GetXaxis()->SetBinLabel(16, "#geq6j,1b"); predSum ->GetXaxis()->SetBinLabel(17, "#geq6j,2b");    	predSum ->GetXaxis()->SetBinLabel(18, "#geq3j,#geq3b");
	predSum ->GetXaxis()->SetBinLabel(19, "2j,0b");     predSum ->GetXaxis()->SetBinLabel(20, "2j,#geq1b");		predSum ->GetXaxis()->SetBinLabel(21, "3-5j,0b");
	predSum ->GetXaxis()->SetBinLabel(22, "3-5j,1b");   predSum ->GetXaxis()->SetBinLabel(23, "3-5j,2b");		predSum ->GetXaxis()->SetBinLabel(24, "#geq6j,0b");
	predSum ->GetXaxis()->SetBinLabel(25, "#geq6j,1b"); predSum ->GetXaxis()->SetBinLabel(26, "#geq6j,2b");		predSum ->GetXaxis()->SetBinLabel(27, "#geq3j,#geq3b");
    
	//add single bin predictions to summary plot histograms
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		string hs = string("_") + HT_bin[i3] + string("_") + signal_region[i2];
		string mapnameZ     = "MT2PredZinv" + hs;
		double Zval(0.); double Zvalerr(0.); double Zval20percenterr(0.); double Zval20percenterr2(0.);
		for(int in = 1; in<=histos[mapnameZ]->GetNbinsX();++in){
			Zval += histos[mapnameZ]->GetBinContent(in);
			if(dummyy){//use at own risk - tried to implement the correlation on Z/gamma ratio
			Zvalerr += pow(histos[mapnameZ]->GetBinError(in),2);//error 
			Zval20percenterr += (0.2)*histos[mapnameZ]->GetBinError(in);//error only 20% syst on R(Z/g)
			Zval20percenterr2 += pow((0.2)*histos[mapnameZ]->GetBinError(in),2);//error only 20% syst on R(Z/g)
			} else{
			Zvalerr += pow(histos[mapnameZ]->GetBinError(in),2);
			}
		}
	//	if(dummyy) Zvalerr = sqrt(Zvalerr+Zval20percenterr*Zval20percenterr);
		if(dummyy) {Zvalerr = sqrt(Zvalerr-Zval20percenterr2+Zval20percenterr*Zval20percenterr);
	//	cout << Zval20percenterr2 << " " << Zval20percenterr*Zval20percenterr << endl;
		}
		else       Zvalerr = sqrt(Zvalerr);
		string mapnameQCD   = "MT2PredQCD" + hs;
		double QCDval(0.); double QCDvalerr(0.);
		for(int in = 1; in<=histos[mapnameQCD]->GetNbinsX();++in){
			QCDval += histos[mapnameQCD]->GetBinContent(in);
			QCDvalerr += pow(histos[mapnameQCD]->GetBinError(in),2);
		}
		QCDvalerr = sqrt(QCDvalerr);
		string mapnameLL    = "MT2PredLostLepton" + hs;
		double LLval(0.); double LLvalerr(0.);//note: LLerr is relative:
		//Delta t /t = Delta bi / bi --> Sum(Delta bi) = sum(bi) * Delta t / t = Delta t
		for(int in = 1; in<=histos[mapnameLL]->GetNbinsX();++in){
			LLval += histos[mapnameLL]->GetBinContent(in);
			LLvalerr += histos[mapnameLL]->GetBinError(in);//no pow(,2) !!
		}
		string mapnameD = "MT2Data" + hs;
		double Dval(0.); double Dvalerr(0.);
		for(int in = 1; in<=histos[mapnameD]->GetNbinsX();++in){
			Dval += histos[mapnameD]->GetBinContent(in);
		}
		Dvalerr = sqrt(Dval);//gaussian error
		if(i3==0){//lHT
			predQCDlHT ->SetBinContent(i2+1, QCDval);            predQCDlHT ->SetBinError(i2+1, QCDvalerr);
			predLLlHT  ->SetBinContent(i2+1, LLval );            predLLlHT  ->SetBinError(i2+1, LLvalerr );
			predZinvlHT->SetBinContent(i2+1, Zval  );            predZinvlHT->SetBinError(i2+1, Zvalerr  );
			predSumlHT ->SetBinContent(i2+1, QCDval+LLval+Zval); predSumlHT ->SetBinError(i2+1, sqrt(pow(QCDvalerr,2)+pow(LLvalerr,2)+pow(Zvalerr,2) ) );
			predDatalHT->SetBinContent(i2+1, Dval  );            predDatalHT->SetBinError(i2+1, Dvalerr  );
            		predQCD ->SetBinContent(i2+1, QCDval);            predQCD ->SetBinError(i2+1, QCDvalerr);
			predLL  ->SetBinContent(i2+1, LLval );            predLL  ->SetBinError(i2+1, LLvalerr );
			predZinv->SetBinContent(i2+1, Zval  );            predZinv->SetBinError(i2+1, Zvalerr  );
			predSum ->SetBinContent(i2+1, QCDval+LLval+Zval); predSum ->SetBinError(i2+1, sqrt(pow(QCDvalerr,2)+pow(LLvalerr,2)+pow(Zvalerr,2) ) );
			predData->SetBinContent(i2+1, Dval  );            predData->SetBinError(i2+1, Dvalerr  );
		} if(i3==1){//mHT
			predQCDmHT ->SetBinContent(i2+1, QCDval);            predQCDmHT ->SetBinError(i2+1, QCDvalerr);
			predLLmHT  ->SetBinContent(i2+1, LLval );            predLLmHT  ->SetBinError(i2+1, LLvalerr );
			predZinvmHT->SetBinContent(i2+1, Zval  );            predZinvmHT->SetBinError(i2+1, Zvalerr  );
			predSummHT ->SetBinContent(i2+1, QCDval+LLval+Zval); predSummHT ->SetBinError(i2+1, sqrt(pow(QCDvalerr,2)+pow(LLvalerr,2)+pow(Zvalerr,2) ) );
			predDatamHT->SetBinContent(i2+1, Dval  );            predDatamHT->SetBinError(i2+1, Dvalerr  );
            		predQCD ->SetBinContent(i2+10, QCDval);            predQCD ->SetBinError(i2+10, QCDvalerr);
			predLL  ->SetBinContent(i2+10, LLval );            predLL  ->SetBinError(i2+10, LLvalerr );
			predZinv->SetBinContent(i2+10, Zval  );            predZinv->SetBinError(i2+10, Zvalerr  );
			predSum ->SetBinContent(i2+10, QCDval+LLval+Zval); predSum ->SetBinError(i2+10, sqrt(pow(QCDvalerr,2)+pow(LLvalerr,2)+pow(Zvalerr,2) ) );
			predData->SetBinContent(i2+10, Dval  );            predData->SetBinError(i2+10, Dvalerr  );
		} if(i3==2){//hHT
			predQCDhHT ->SetBinContent(i2+1, QCDval);            predQCDhHT ->SetBinError(i2+1, QCDvalerr);
			predLLhHT  ->SetBinContent(i2+1, LLval );            predLLhHT  ->SetBinError(i2+1, LLvalerr );
			predZinvhHT->SetBinContent(i2+1, Zval  );            predZinvhHT->SetBinError(i2+1, Zvalerr  );
			predSumhHT ->SetBinContent(i2+1, QCDval+LLval+Zval); predSumhHT ->SetBinError(i2+1, sqrt(pow(QCDvalerr,2)+pow(LLvalerr,2)+pow(Zvalerr,2) ) );
			predDatahHT->SetBinContent(i2+1, Dval  );            predDatahHT->SetBinError(i2+1, Dvalerr  );
            		predQCD ->SetBinContent(i2+19, QCDval);            predQCD ->SetBinError(i2+19, QCDvalerr);
			predLL  ->SetBinContent(i2+19, LLval );            predLL  ->SetBinError(i2+19, LLvalerr );
			predZinv->SetBinContent(i2+19, Zval  );            predZinv->SetBinError(i2+19, Zvalerr  );
			predSum ->SetBinContent(i2+19, QCDval+LLval+Zval); predSum ->SetBinError(i2+19, sqrt(pow(QCDvalerr,2)+pow(LLvalerr,2)+pow(Zvalerr,2) ) );
			predData->SetBinContent(i2+19, Dval  );            predData->SetBinError(i2+19, Dvalerr  );
		}

	}}
    TString histoTitle2 = ";signal region;ratio of data / prediction";
    predSum->SetTitle(histoTitle2);
    predData->SetTitle(histoTitle2);
    TH1D* predRatioR = (TH1D*)predSum ->Clone("predRatioR");
    TH1D* predRatio  = (TH1D*)predData->Clone("predRatio");
    for(int i=1;i<=predRatioR->GetNbinsX();++i){
	if(predSum->GetBinContent(i)!=0){
		predRatioR->SetBinContent(i, 1.);
		predRatioR->SetBinError(  i, predSum ->GetBinError(  i)/predSum->GetBinContent(i));
		predRatio ->SetBinContent(i, predData->GetBinContent(i)/predSum->GetBinContent(i));
		predRatio ->SetBinError(  i, predData->GetBinError(  i)/predSum->GetBinContent(i));
	} else {
		predRatioR->SetBinContent(i, 1.);
		predRatioR->SetBinError(  i, 0.);
		predRatio ->SetBinContent(i, 1.);
		predRatio ->SetBinError(  i, 0.);
	}
    }

	//for three summary plots we scale the 2j0b and 3to5j0b yield down, so that other bins are also visible, scales are hardcoded
	if(summaryscale2bins){
		predQCDlHT ->SetBinContent(1, predQCDlHT ->GetBinContent(1)*(1./6.)); predQCDlHT ->SetBinError(1, predQCDlHT ->GetBinError(1)*(1./6.));
		predQCDmHT ->SetBinContent(1, predQCDmHT ->GetBinContent(1)*(1./4.)); predQCDmHT ->SetBinError(1, predQCDmHT ->GetBinError(1)*(1./4.));
		predQCDhHT ->SetBinContent(1, predQCDhHT ->GetBinContent(1)*(1./2.)); predQCDhHT ->SetBinError(1, predQCDhHT ->GetBinError(1)*(1./2.));
		predLLlHT  ->SetBinContent(1, predLLlHT  ->GetBinContent(1)*(1./6.)); predLLlHT  ->SetBinError(1, predLLlHT  ->GetBinError(1)*(1./6.));
		predLLmHT  ->SetBinContent(1, predLLmHT  ->GetBinContent(1)*(1./4.)); predLLmHT  ->SetBinError(1, predLLmHT  ->GetBinError(1)*(1./4.));
		predLLhHT  ->SetBinContent(1, predLLhHT  ->GetBinContent(1)*(1./2.)); predLLhHT  ->SetBinError(1, predLLhHT  ->GetBinError(1)*(1./2.));
		predZinvlHT->SetBinContent(1, predZinvlHT->GetBinContent(1)*(1./6.)); predZinvlHT->SetBinError(1, predZinvlHT->GetBinError(1)*(1./6.));
		predZinvmHT->SetBinContent(1, predZinvmHT->GetBinContent(1)*(1./4.)); predZinvmHT->SetBinError(1, predZinvmHT->GetBinError(1)*(1./4.));
		predZinvhHT->SetBinContent(1, predZinvhHT->GetBinContent(1)*(1./2.)); predZinvhHT->SetBinError(1, predZinvhHT->GetBinError(1)*(1./2.));
		predSumlHT ->SetBinContent(1, predSumlHT ->GetBinContent(1)*(1./6.)); predSumlHT ->SetBinError(1, predSumlHT ->GetBinError(1)*(1./6.));
		predSummHT ->SetBinContent(1, predSummHT ->GetBinContent(1)*(1./4.)); predSummHT ->SetBinError(1, predSummHT ->GetBinError(1)*(1./4.));
		predSumhHT ->SetBinContent(1, predSumhHT ->GetBinContent(1)*(1./2.)); predSumhHT ->SetBinError(1, predSumhHT ->GetBinError(1)*(1./2.));
		predDatalHT->SetBinContent(1, predDatalHT->GetBinContent(1)*(1./6.)); predDatalHT->SetBinError(1, predDatalHT->GetBinError(1)*(1./6.));
		predDatamHT->SetBinContent(1, predDatamHT->GetBinContent(1)*(1./4.)); predDatamHT->SetBinError(1, predDatamHT->GetBinError(1)*(1./4.));
		predDatahHT->SetBinContent(1, predDatahHT->GetBinContent(1)*(1./2.)); predDatahHT->SetBinError(1, predDatahHT->GetBinError(1)*(1./2.));

		predQCDlHT ->SetBinContent(3, predQCDlHT ->GetBinContent(3)*(1./9.)); predQCDlHT ->SetBinError(3, predQCDlHT ->GetBinError(3)*(1./9.));
		predQCDmHT ->SetBinContent(3, predQCDmHT ->GetBinContent(3)*(1./6.)); predQCDmHT ->SetBinError(3, predQCDmHT ->GetBinError(3)*(1./6.));
		predQCDhHT ->SetBinContent(3, predQCDhHT ->GetBinContent(3)*(1./3.)); predQCDhHT ->SetBinError(3, predQCDhHT ->GetBinError(3)*(1./3.));
		predLLlHT  ->SetBinContent(3, predLLlHT  ->GetBinContent(3)*(1./9.)); predLLlHT  ->SetBinError(3, predLLlHT  ->GetBinError(3)*(1./9.));
		predLLmHT  ->SetBinContent(3, predLLmHT  ->GetBinContent(3)*(1./6.)); predLLmHT  ->SetBinError(3, predLLmHT  ->GetBinError(3)*(1./6.));
		predLLhHT  ->SetBinContent(3, predLLhHT  ->GetBinContent(3)*(1./3.)); predLLhHT  ->SetBinError(3, predLLhHT  ->GetBinError(3)*(1./3.));
		predZinvlHT->SetBinContent(3, predZinvlHT->GetBinContent(3)*(1./9.)); predZinvlHT->SetBinError(3, predZinvlHT->GetBinError(3)*(1./9.));
		predZinvmHT->SetBinContent(3, predZinvmHT->GetBinContent(3)*(1./6.)); predZinvmHT->SetBinError(3, predZinvmHT->GetBinError(3)*(1./6.));
		predZinvhHT->SetBinContent(3, predZinvhHT->GetBinContent(3)*(1./3.)); predZinvhHT->SetBinError(3, predZinvhHT->GetBinError(3)*(1./3.));
		predSumlHT ->SetBinContent(3, predSumlHT ->GetBinContent(3)*(1./9.)); predSumlHT ->SetBinError(3, predSumlHT ->GetBinError(3)*(1./9.));
		predSummHT ->SetBinContent(3, predSummHT ->GetBinContent(3)*(1./6.)); predSummHT ->SetBinError(3, predSummHT ->GetBinError(3)*(1./6.));
		predSumhHT ->SetBinContent(3, predSumhHT ->GetBinContent(3)*(1./3.)); predSumhHT ->SetBinError(3, predSumhHT ->GetBinError(3)*(1./3.));
		predDatalHT->SetBinContent(3, predDatalHT->GetBinContent(3)*(1./9.)); predDatalHT->SetBinError(3, predDatalHT->GetBinError(3)*(1./9.));
		predDatamHT->SetBinContent(3, predDatamHT->GetBinContent(3)*(1./6.)); predDatamHT->SetBinError(3, predDatamHT->GetBinError(3)*(1./6.));
		predDatahHT->SetBinContent(3, predDatahHT->GetBinContent(3)*(1./3.)); predDatahHT->SetBinError(3, predDatahHT->GetBinError(3)*(1./3.));

		predQCDlHT ->SetBinContent(4, predQCDlHT ->GetBinContent(4)*(1./3.)); predQCDlHT ->SetBinError(4, predQCDlHT ->GetBinError(4)*(1./3.));
		predQCDmHT ->SetBinContent(4, predQCDmHT ->GetBinContent(4)*(1./2.)); predQCDmHT ->SetBinError(4, predQCDmHT ->GetBinError(4)*(1./2.));
		predLLlHT  ->SetBinContent(4, predLLlHT  ->GetBinContent(4)*(1./3.)); predLLlHT  ->SetBinError(4, predLLlHT  ->GetBinError(4)*(1./3.));
		predLLmHT  ->SetBinContent(4, predLLmHT  ->GetBinContent(4)*(1./2.)); predLLmHT  ->SetBinError(4, predLLmHT  ->GetBinError(4)*(1./2.));
		predZinvlHT->SetBinContent(4, predZinvlHT->GetBinContent(4)*(1./3.)); predZinvlHT->SetBinError(4, predZinvlHT->GetBinError(4)*(1./3.));
		predZinvmHT->SetBinContent(4, predZinvmHT->GetBinContent(4)*(1./2.)); predZinvmHT->SetBinError(4, predZinvmHT->GetBinError(4)*(1./2.));
		predSumlHT ->SetBinContent(4, predSumlHT ->GetBinContent(4)*(1./3.)); predSumlHT ->SetBinError(4, predSumlHT ->GetBinError(4)*(1./3.));
		predSummHT ->SetBinContent(4, predSummHT ->GetBinContent(4)*(1./2.)); predSummHT ->SetBinError(4, predSummHT ->GetBinError(4)*(1./2.));
		predDatalHT->SetBinContent(4, predDatalHT->GetBinContent(4)*(1./3.)); predDatalHT->SetBinError(4, predDatalHT->GetBinError(4)*(1./3.));
		predDatamHT->SetBinContent(4, predDatamHT->GetBinContent(4)*(1./2.)); predDatamHT->SetBinError(4, predDatamHT->GetBinError(4)*(1./2.));

	}

	SumlHT->Add(predQCDlHT ); SumlHT->Add(predLLlHT  ); SumlHT->Add(predZinvlHT);
	SummHT->Add(predQCDmHT ); SummHT->Add(predLLmHT  ); SummHT->Add(predZinvmHT);
	SumhHT->Add(predQCDhHT ); SumhHT->Add(predLLhHT  ); SumhHT->Add(predZinvhHT);
	Sum->Add(predQCD ); Sum->Add(predLL  ); Sum->Add(predZinv);
	TString histoTitle = ";signal region;Events";
	SumlHT->SetTitle(histoTitle); SummHT->SetTitle(histoTitle); SumhHT->SetTitle(histoTitle); 
	Sum->SetTitle(histoTitle);

	//write to file
	TFile *newfile = new TFile(newfilename.Data(), "RECREATE");
	newfile->cd();
	leg->Write();
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Write();}
	for(map<string,THStack*>::iterator h=stackPrediction.begin(); h!=stackPrediction.end();++h){
		h->second->Write();}
		predQCDlHT ->Write();
		predQCDmHT ->Write();
		predQCDhHT ->Write();
		predLLlHT  ->Write();
		predLLmHT  ->Write();
		predLLhHT  ->Write();
		predZinvlHT->Write();
		predZinvmHT->Write();
		predZinvhHT->Write();
		predSumlHT  ->Write();
		predSummHT  ->Write();
		predSumhHT  ->Write();
		predDatalHT->Write();
		predDatamHT->Write();
		predDatahHT->Write();
		SumlHT      ->Write();
		SummHT      ->Write();
		SumhHT      ->Write();
		predQCD ->Write();
		predLL  ->Write();
		predZinv->Write();
		predSum  ->Write();
		predData->Write();
		Sum      ->Write();
	newfile->Close();

	//use later for LEE studies
	TFile *secondfile = new TFile(allfile.Data(), "RECREATE");
	secondfile->cd();
	hZinvAll->Write();
	hLLAll  ->Write();
	hQCDAll ->Write();
	hPredAll->Write();
	hDataAll->Write();
	secondfile->Close();

	//another option to set poissonian uncertainties also for summary plots - again should become obsolete in next ROOT version
	TGraphAsymmErrors *gpredData = new TGraphAsymmErrors((TH1D*)predData);
	TGraphAsymmErrors *gpredDatalHT = new TGraphAsymmErrors((TH1D*)predDatalHT);
	TGraphAsymmErrors *gpredDatamHT = new TGraphAsymmErrors((TH1D*)predDatamHT);
	TGraphAsymmErrors *gpredDatahHT = new TGraphAsymmErrors((TH1D*)predDatahHT);
	for(int ig = 0; ig<predData->GetNbinsX(); ++ig){
		int nnn = ig;
		gpredData->SetPoint(nnn, predData->GetBinCenter(ig+1), predData->GetBinContent(ig+1));
		double pointerrorX = 0.5*predData->GetBinWidth(ig+1);
		double pointerrorYlow = statErrorN(predData->GetBinContent(ig+1));
		double pointerrorYup = statErrorP(predData->GetBinContent(ig+1));
		gpredData->SetPointError(nnn, pointerrorX, pointerrorX, pointerrorYlow, pointerrorYup);
		gpredData->SetMarkerStyle(20);
		gpredData->SetLineWidth(2);
		gpredData->SetLineColor(kBlack);
		gpredData->SetMarkerColor(kBlack);
	}
	for(int ig = 0; ig<predDatalHT->GetNbinsX(); ++ig){
		int nnn = ig;
		gpredDatalHT->SetPoint(nnn, predDatalHT->GetBinCenter(ig+1), predDatalHT->GetBinContent(ig+1));
		double pointerrorX = 0.5*predDatalHT->GetBinWidth(ig+1);
		double pointerrorYlow = statErrorN(predDatalHT->GetBinContent(ig+1));
		double pointerrorYup = statErrorP(predDatalHT->GetBinContent(ig+1));
		gpredDatalHT->SetPointError(nnn, pointerrorX, pointerrorX, pointerrorYlow, pointerrorYup);
		gpredDatalHT->SetMarkerStyle(20);
		gpredDatalHT->SetLineWidth(2);
		gpredDatalHT->SetLineColor(kBlack);
		gpredDatalHT->SetMarkerColor(kBlack);
	}
	for(int ig = 0; ig<predDatamHT->GetNbinsX(); ++ig){
		int nnn = ig;
		gpredDatamHT->SetPoint(nnn, predDatamHT->GetBinCenter(ig+1), predDatamHT->GetBinContent(ig+1));
		double pointerrorX = 0.5*predDatamHT->GetBinWidth(ig+1);
		double pointerrorYlow = statErrorN(predDatamHT->GetBinContent(ig+1));
		double pointerrorYup = statErrorP(predDatamHT->GetBinContent(ig+1));
		gpredDatamHT->SetPointError(nnn, pointerrorX, pointerrorX, pointerrorYlow, pointerrorYup);
		gpredDatamHT->SetMarkerStyle(20);
		gpredDatamHT->SetLineWidth(2);
		gpredDatamHT->SetLineColor(kBlack);
		gpredDatamHT->SetMarkerColor(kBlack);
	}
	for(int ig = 0; ig<predDatahHT->GetNbinsX(); ++ig){
		int nnn = ig;
		gpredDatahHT->SetPoint(nnn, predDatahHT->GetBinCenter(ig+1), predDatahHT->GetBinContent(ig+1));
		double pointerrorX = 0.5*predDatahHT->GetBinWidth(ig+1);
		double pointerrorYlow = statErrorN(predDatahHT->GetBinContent(ig+1));
		double pointerrorYup = statErrorP(predDatahHT->GetBinContent(ig+1));
		gpredDatahHT->SetPointError(nnn, pointerrorX, pointerrorX, pointerrorYlow, pointerrorYup);
		gpredDatahHT->SetMarkerStyle(20);
		gpredDatahHT->SetLineWidth(2);
		gpredDatahHT->SetLineColor(kBlack);
		gpredDatahHT->SetMarkerColor(kBlack);
	}


	cout << endl << endl << "now the results table:" << endl << endl;

	//this first table is a 'summary table': numbers without MT2 binning - table does not contain MC yields, even if flag is set
	cout << "\%BEGINLATEX\%"             << endl;
	cout << "\\begin{table}"             << endl
	<< "\\begin{center}"            << endl;
	string htreg = "";
	cout << "\\caption{Estimated background event yields in the various $H_\\mathrm{T}$ and topological regions. The uncertainties are the quadratic sum of statistical and systematic uncertainties.}" << endl;
	cout << "\\label{table:datadrivenMT2_Summary}" << endl;
	if(!useMCyields) cout << "\\begin{tabular}{|l|l|c|c|c|c|c|}" << endl;	     
	else             cout << "\\begin{tabular}{|l|l|cc|cc|cc|cc|c|}" << endl;
	cout << "\\hline\\hline"             << endl;
	int ncolumns = 7; if(useMCyields) ncolumns = 11;
	if(useMCyields){
		cout << " signal & lowest & \\multicolumn{2}{|c|}{Multijet} & \\multicolumn{2}{|c|}{Lost lepton} & \\multicolumn{2}{|c|}{$Z(\\nu\\bar{\\nu})$+jets} & \\multicolumn{2}{|c|}{Total bkg.} & Data \\\\" << endl;
		cout << " region & $M_\\mathrm{T2}$ & MC & data pred. & MC & data pred. & MC & data pred. & MC & data pred. & \\\\" << endl;
	} else {
		cout << " signal region & lowest $M_\\mathrm{T2}$ & Multijet & Lost lepton & $Z(\\nu\\bar{\\nu})$+jets & Total bkg. & Data \\\\" << endl;
	}
	cout << "\\hline\\hline"             << endl;
	//low HT
	cout << " \\multicolumn{" << ncolumns << "}{|c|}{low $H_\\mathrm{T}$ region: $450~\\mathrm{GeV} \\leq H_\\mathrm{T} < 750~\\mathrm{GeV}$, $\\MET > 200~\\mathrm{GeV}$} \\\\" << endl;
	cout << "\\hline";
	for(int i2 = 0;  i2<signalregionsize; ++i2){
		if(i2==0) cout << " $2$ jets, $0$ b jets & 200 GeV";
		if(i2==1) cout << " $2$ jets, $\\geq 1$ b jets & 200 GeV";
		if(i2==2) cout << " $3-5$ jets, $0$ b jets & 200 GeV";
		if(i2==3) cout << " $3-5$ jets, $1$ b jets & 200 GeV";
		if(i2==4) cout << " $3-5$ jets, $2$ b jets & 200 GeV";
		if(i2==5) cout << " $\\geq 6$ jets, $0$ b jets & 200 GeV";
		if(i2==6) cout << " $\\geq 6$ jets, $1$ b jets & 200 GeV";
		if(i2==7) cout << " $\\geq 6$ jets, $2$ b jets & 200 GeV";
		if(i2==8) cout << " $\\geq 3$ jets, $\\geq 3$ b jets & 200 GeV";
		cout << fixed << setprecision(2)
		     << " & " << predQCDlHT->GetBinContent(i2+1) << "$\\pm$" << predQCDlHT->GetBinError(i2+1);
		cout << fixed << setprecision(2)
		     << " & " << predLLlHT->GetBinContent(i2+1) << "$\\pm$" << predLLlHT->GetBinError(i2+1)
		     << " & " << predZinvlHT->GetBinContent(i2+1) << "$\\pm$" << predZinvlHT->GetBinError(i2+1)
		     << " & " << predSumlHT->GetBinContent(i2+1) << "$\\pm$" << predSumlHT->GetBinError(i2+1)
		     << " & " << int(predDatalHT->GetBinContent(i2+1)) << " \\\\" << endl;
	}
	cout << "\\hline\\hline" << endl;
	//medium HT
	cout << " \\multicolumn{" << ncolumns << "}{|c|}{medium $H_\\mathrm{T}$ region: $750~\\mathrm{GeV} \\leq H_\\mathrm{T} < 1200~\\mathrm{GeV}$, $\\MET > 30~\\mathrm{GeV}$} \\\\" << endl;
	cout << "\\hline" << endl;
	for(int i2 = 0;  i2<signalregionsize; ++i2){
		if(i2==0) cout << " $2$ jets, $0$ b jets & 125 GeV";
		if(i2==1) cout << " $2$ jets, $\\geq 1$ b jets & 100 GeV";
		if(i2==2) cout << " $3-5$ jets, $0$ b jets & 160 GeV";
		if(i2==3) cout << " $3-5$ jets, $1$ b jets & 150 GeV";
		if(i2==4) cout << " $3-5$ jets, $2$ b jets & 130 GeV";
		if(i2==5) cout << " $\\geq 6$ jets, $0$ b jets & 160 GeV";
		if(i2==6) cout << " $\\geq 6$ jets, $1$ b jets & 150 GeV";
		if(i2==7) cout << " $\\geq 6$ jets, $2$ b jets & 130 GeV";
		if(i2==8) cout << " $\\geq 3$ jets, $\\geq 3$ b jets & 125 GeV";
		cout << fixed << setprecision(2)
		     << " & " << predQCDmHT->GetBinContent(i2+1) << "$\\pm$" << predQCDmHT->GetBinError(i2+1);
		cout << fixed << setprecision(2)
		     << " & " << predLLmHT->GetBinContent(i2+1) << "$\\pm$" << predLLmHT->GetBinError(i2+1)
		     << " & " << predZinvmHT->GetBinContent(i2+1) << "$\\pm$" << predZinvmHT->GetBinError(i2+1)
		     << " & " << predSummHT->GetBinContent(i2+1) << "$\\pm$" << predSummHT->GetBinError(i2+1)
		     << " & " << int(predDatamHT->GetBinContent(i2+1)) << " \\\\" << endl;
	}
	cout << "\\hline\\hline" << endl;
	//high HT
	cout << " \\multicolumn{" << ncolumns << "}{|c|}{high $H_\\mathrm{T}$ region: $H_\\mathrm{T} \\geq 1200~\\mathrm{GeV}$, $\\MET > 30~\\mathrm{GeV}$} \\\\" << endl;
	cout << "\\hline" << endl;
	for(int i2 = 0;  i2<signalregionsize; ++i2){
		if(i2==0) cout << " $2$ jets, $0$ b jets & 120 GeV";
		if(i2==1) cout << " $2$ jets, $\\geq 1$ b jets & 100 GeV";
		if(i2==2) cout << " $3-5$ jets, $0$ b jets & 160 GeV";
		if(i2==3) cout << " $3-5$ jets, $1$ b jets & 150 GeV";
		if(i2==4) cout << " $3-5$ jets, $2$ b jets & 130 GeV";
		if(i2==5) cout << " $\\geq 6$ jets, $0$ b jets & 160 GeV";
		if(i2==6) cout << " $\\geq 6$ jets, $1$ b jets & 150 GeV";
		if(i2==7) cout << " $\\geq 6$ jets, $2$ b jets & 130 GeV";
		if(i2==8) cout << " $\\geq 3$ jets, $\\geq 3$ b jets & 125 GeV";
		cout << fixed << setprecision(2)
		     << " & " << predQCDhHT->GetBinContent(i2+1) << "$\\pm$" << predQCDhHT->GetBinError(i2+1);
		cout << fixed << setprecision(2)
		     << " & " << predLLhHT->GetBinContent(i2+1) << "$\\pm$" << predLLhHT->GetBinError(i2+1)
		     << " & " << predZinvhHT->GetBinContent(i2+1) << "$\\pm$" << predZinvhHT->GetBinError(i2+1)
		     << " & " << predSumhHT->GetBinContent(i2+1) << "$\\pm$" << predSumhHT->GetBinError(i2+1)
		     << " & " << int(predDatahHT->GetBinContent(i2+1)) << " \\\\" << endl;
	}
	cout << "\\hline\\hline"             << endl;
	cout << "\\end{tabular}" << endl
	     << "\\end{center}"  << endl
	     << "\\end{table}"   << endl
	     << "\%ENDLATEX\%"   << endl
	     << endl;




	//the detailed tables - all signal bins splitted in all backgrounds - has also option to plot MC expectation
	//NOTE - within this table code is the filling of the pulls !!!!!!!!! don't delete that
	cout << endl << endl << "now the results table:" << endl << endl;
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		cout << "\%BEGINLATEX\%"             << endl;
		cout << "\\begin{table}"             << endl
		<< "\\begin{center}"            << endl;
	//	cout << "\\small"                    << endl;
		string httreg = "";
		if(i3==0) httreg = "low"; if(i3==1) httreg = "medium"; if(i3==2) httreg = "high";
		cout << "\\caption{Estimated background event yields in the various $M_\\mathrm{T2}$ and topological bins for the " << httreg << " $H_\\mathrm{T}$ region. The uncertainties are the quadratic sum of statistical and systematic uncertainties.}" << endl;
		cout << "\\label{table:datadrivenMT2_" << httreg << "HT}" << endl;
		if(!useMCyields) cout << "\\begin{tabular}{|l|c|c|c|c|c|}" << endl;	     
		else             cout << "\\begin{tabular}{|l|cc|cc|cc|cc|c|}" << endl;
		cout << "\\hline\\hline"             << endl;
		int nncolumns = 6; if(useMCyields) nncolumns = 10;
		if(i3==0) cout << " \\multicolumn{" << nncolumns << "}{|c|}{low $H_\\mathrm{T}$ region: $450~\\mathrm{GeV} \\leq H_\\mathrm{T} < 750~\\mathrm{GeV}$, $\\MET > 200~\\mathrm{GeV}$} \\\\" << endl;
		if(i3==1) cout << " \\multicolumn{" << nncolumns << "}{|c|}{medium $H_\\mathrm{T}$ region: $750~\\mathrm{GeV} \\leq H_\\mathrm{T} < 1200~\\mathrm{GeV}$, $\\MET > 30~\\mathrm{GeV}$} \\\\" << endl;
		if(i3==2) cout << " \\multicolumn{" << nncolumns << "}{|c|}{high $H_\\mathrm{T}$ region: $H_\\mathrm{T} \\geq 1200~\\mathrm{GeV}$, $\\MET > 30~\\mathrm{GeV}$} \\\\" << endl;
		cout << "\\hline\\hline"             << endl;
		if(useMCyields){
			cout << " & \\multicolumn{2}{|c|}{Multijet} & \\multicolumn{2}{|c|}{Lost lepton} & \\multicolumn{2}{|c|}{$Z(\\nu\\bar{\\nu})$+jets} & \\multicolumn{2}{|c|}{Total bkg.} & Data \\\\" << endl;
			cout << " $M_\\mathrm{T2}$ [GeV] & sim. & data pred. & sim. & data pred. & sim. & data pred. & sim. & data pred. & \\\\" << endl;
		} else {
			cout << " $M_\\mathrm{T2}$ [GeV] & Multijet & Lost lepton & $Z(\\nu\\bar{\\nu})$+jets & Total bkg. & Data \\\\" << endl;
		}
		cout << "\\hline\\hline"             << endl;
	for(int i2 = 0; i2<signalregionsize; ++i2){
		cout << " & \\multicolumn{" << nncolumns-1 << "}{|l|}{";
		if(i2==0) cout << "$2$ jets, $0$ b jets";
		if(i2==1) cout << "$2$ jets, $\\geq 1$ b jets";
		if(i2==2) cout << "$3-5$ jets, $0$ b jets";
		if(i2==3) cout << "$3-5$ jets, $1$ b jets";
		if(i2==4) cout << "$3-5$ jets, $2$ b jets";
		if(i2==5) cout << "$\\geq 6$ jets, $0$ b jets";
		if(i2==6) cout << "$\\geq 6$ jets, $1$ b jets";
		if(i2==7) cout << "$\\geq 6$ jets, $2$ b jets";
		if(i2==8) cout << "$\\geq 3$ jets, $\\geq 3$ b jets";
		cout << "} \\\\" << endl;
		cout << "\\hline" << endl;
		string hs = string("_") + HT_bin[i3] + string("_") + signal_region[i2];
		string mapnameZ     = "MT2PredZinv" + hs;
		string mapnameLL    = "MT2PredLostLepton" + hs;
		string mapnameQCD   = "MT2PredQCD" + hs;
		string mapnamePred  = "MT2PredSum" + hs;
		string mapnameData  = "MT2Data" + hs;
		for(int n = 1; n<=histos[mapnameData]->GetNbinsX(); ++n){
			cout << " $" << int(histos[mapnameData]->GetBinLowEdge(n)) << "-";
			if(n!=histos[mapnameData]->GetNbinsX()) cout << int(histos[mapnameData]->GetBinLowEdge(n)+histos[mapnameData]->GetBinWidth(n)) << "$";
			else cout << "$Inf";
			if(useMCyields){
				string mapnameMCQCD = "MT2MCQCD" + hs;
				string mapnameMCW   = "MT2MCW" + hs;
				string mapnameMCTop = "MT2MCTop" + hs;
				string mapnameMCZ   = "MT2MCZ" + hs;
				string mapnameMCO   = "MT2MCOther" + hs;
				string mapnameMCsum = "MT2MCSum" + hs;
				double LLyield(0.);// double LLerr(0.);
				LLyield = histosMC[mapnameMCW]->GetBinContent(n) + histosMC[mapnameMCTop]->GetBinContent(n) + histosMC[mapnameMCO]->GetBinContent(n);
				cout << fixed << setprecision(2)
				     << " & " << histosMC[mapnameMCQCD]->GetBinContent(n);
				cout << fixed << setprecision(2)
				     << " & " << histos[mapnameQCD]->GetBinContent(n) << "$\\pm$" << histos[mapnameQCD]->GetBinError(n);
				cout << fixed << setprecision(2)
				     << " & " << LLyield
				     << " & " << histos[mapnameLL]->GetBinContent(n) << "$\\pm$" << histos[mapnameLL]->GetBinError(n)
				     << " & " << histosMC[mapnameMCZ]->GetBinContent(n)
				     << " & " << histos[mapnameZ]->GetBinContent(n) << "$\\pm$" << histos[mapnameZ]->GetBinError(n)
				     << " & " << histosMC[mapnameMCsum]->GetBinContent(n)
				     << " & " << histos[mapnamePred]->GetBinContent(n) << "$\\pm$" << histos[mapnamePred]->GetBinError(n)
				     << " & " << int(histos[mapnameData]->GetBinContent(n)) << " \\\\" << endl;
			}
			else {
				cout << fixed << setprecision(2)
				     << " & " << histos[mapnameQCD]->GetBinContent(n) << "$\\pm$" << histos[mapnameQCD]->GetBinError(n);
				cout << fixed << setprecision(2)
				     << " & " << histos[mapnameLL]->GetBinContent(n) << "$\\pm$" << histos[mapnameLL]->GetBinError(n)
				     << " & " << histos[mapnameZ]->GetBinContent(n) << "$\\pm$" << histos[mapnameZ]->GetBinError(n)
				     << " & " << histos[mapnamePred]->GetBinContent(n) << "$\\pm$" << histos[mapnamePred]->GetBinError(n)
				     << " & " << int(histos[mapnameData]->GetBinContent(n)) << " \\\\" << endl;
			}
			//for pulls
			if(histos[mapnamePred]->GetBinContent(n)>0&&histos[mapnameData]->GetBinContent(n)>0){
				double pulldata = (histos[mapnameData]->GetBinContent(n) - histos[mapnamePred]->GetBinContent(n)) / sqrt(pow(histos[mapnamePred]->GetBinError(n),2) + histos[mapnameData]->GetBinContent(n) );
				if(poisson){//this is a private implementation of pulls that seem to be correct for poissonian uncertainties - however I am not sure, use at own risk
				//POISSONIAN PULL
				if(histos[mapnameData]->GetBinContent(n) - histos[mapnamePred]->GetBinContent(n)>0)
				 pulldata = (histos[mapnameData]->GetBinContent(n) - histos[mapnamePred]->GetBinContent(n)) / sqrt(pow(histos[mapnamePred]->GetBinError(n),2) + pow(statErrorN(histos[mapnamePred]->GetBinContent(n)),2) );
				else
				pulldata = (histos[mapnameData]->GetBinContent(n) - histos[mapnamePred]->GetBinContent(n)) / sqrt(pow(histos[mapnamePred]->GetBinError(n),2) + pow(statErrorP(histos[mapnamePred]->GetBinContent(n)),2) );
				}
				hPullData->Fill(pulldata);
			}
			string mapnameMCSum = "MT2MCSum" + hs;
			if(useMCyields && histos[mapnamePred]->GetBinContent(n)>0&&histosMC[mapnameMCSum]->GetBinContent(n)>0){
				double pulldata = (histosMC[mapnameMCSum]->GetBinContent(n) - histos[mapnamePred]->GetBinContent(n)) / sqrt(pow(histos[mapnamePred]->GetBinError(n),2) + pow(histosMC[mapnameMCSum]->GetBinError(n),2) );
				hPullMC->Fill(pulldata);
			} if(useMCyields && histos[mapnameData]->GetBinContent(n)>0&&histosMC[mapnameMCSum]->GetBinContent(n)>0){
				double pulldata = (histos[mapnameData]->GetBinContent(n) - histosMC[mapnameMCSum]->GetBinContent(n)) / sqrt(histos[mapnameData]->GetBinContent(n) + pow(histosMC[mapnameMCSum]->GetBinError(n),2) );
				hPullMCDa->Fill(pulldata);
			}
		}//for(int n = 1; n<=histos[mapnameData]->GetNbinsX(); ++n)
		cout << "\\hline\\hline" << endl;
	}
	cout << "\\end{tabular}" << endl
	     << "\\end{center}"  << endl
	     << "\\end{table}"   << endl
	     << "\%ENDLATEX\%"   << endl
	     << endl;
	}


	cout << endl << endl << "now the results table summary:" << endl << endl;
	//grand prediction table - all bins - but not splitted to each background component
	cout << "\%BEGINLATEX\%"             << endl;
	cout << "\\begin{table}" << endl;
	cout << "\\begin{center}" << endl;
	cout << "\\caption{Estimated background event yields in all regions. The uncertainties are the quadratic sum of statistical and systematic uncertainties.}" << endl;
	cout << "\\label{table:datadrivenMT2_short}" << endl;
	cout << "\\setlength{\\tabcolsep}{4pt}" << endl;
	cout << "\\begin{tabular}{|l||l|c|c||l|c|c||l|c|c|}" << endl;
	cout << "\\hline\\hline" << endl;
	cout << "\\multirow{2}{*}{\\begin{minipage}{1.5cm}signal\\\\ region\\end{minipage}}" << endl;
	cout << " & \\multicolumn{3}{|c||}{low $H_\\mathrm{T}$ region} & \\multicolumn{3}{|c||}{medium $H_\\mathrm{T}$ region} & \\multicolumn{3}{|c|}{high $H_\\mathrm{T}$ region} \\\\" << endl;
	cout << "\\cline{2-10}" << endl;
	cout << " & $M_\\mathrm{T2}$ [GeV] & prediction & data & $M_\\mathrm{T2}$ [GeV] & prediction & data & $M_\\mathrm{T2}$ [GeV] & prediction & data \\\\" << endl;
	cout << "\\hline\\hline" << endl;
	for(int i2 = 0; i2<signalregionsize; ++i2){
		string srname;
		if(i2==0) srname = "$2$ jets,\\\\ $0$ b jets";
		if(i2==1) srname = "$2$ jets,\\\\ $\\geq 1$ b jets";
		if(i2==2) srname = "$3-5$ jets,\\\\ $0$ b jets";
		if(i2==3) srname = "$3-5$ jets,\\\\ $1$ b jets";
		if(i2==4) srname = "$3-5$ jets,\\\\ $2$ b jets";
		if(i2==5) srname = "$\\geq 6$ jets,\\\\ $0$ b jets";
		if(i2==6) srname = "$\\geq 6$ jets,\\\\ $1$ b jets";
		if(i2==7) srname = "$\\geq 6$ jets,\\\\ $2$ b jets";
		if(i2==8) srname = "$\\geq 3$ jets,\\\\ $\\geq 3$ b jets";
		string hsl = string("_") + HT_bin[0] + string("_") + signal_region[i2];
		string hsm = string("_") + HT_bin[1] + string("_") + signal_region[i2];
		string hsh = string("_") + HT_bin[2] + string("_") + signal_region[i2];
		string mapnamePredl = "MT2PredSum" + hsl;
		string mapnameDatal = "MT2Data" + hsl;
		string mapnamePredm = "MT2PredSum" + hsm;
		string mapnameDatam = "MT2Data" + hsm;
		string mapnamePredh = "MT2PredSum" + hsh;
		string mapnameDatah = "MT2Data" + hsh;
		//low or medium HT has always higher Nbins than high HT
		int maxbin = (histos[mapnamePredl]->GetNbinsX()>histos[mapnamePredm]->GetNbinsX())?histos[mapnamePredl]->GetNbinsX():histos[mapnamePredm]->GetNbinsX();
		cout << "\\multirow{" << maxbin << "}{*}{\\begin{minipage}{1.5cm}" << srname << "\\end{minipage}}" << endl;
		for(int in = 1; in <=maxbin; ++in){
			//low HT
			if(in>histos[mapnamePredl]->GetNbinsX()){
				cout << " &                    &                                &         ";
			} else {
				if(in==histos[mapnamePredl]->GetNbinsX())
					cout << " & $\\geq"<< (int)histos[mapnamePredl]->GetBinLowEdge(in) << "$ ";
				else 
				cout << " & $"<< (int)histos[mapnamePredl]->GetBinLowEdge(in)<< "-" << (int)histos[mapnamePredl]->GetBinLowEdge(in+1) << "$ ";
				cout << "& " << histos[mapnamePredl]->GetBinContent(in) << "$\\pm$" << histos[mapnamePredl]->GetBinError(in) << " ";
				cout << "& " << (int)histos[mapnameDatal]->GetBinContent(in) << " ";
			}
			//medium HT
			if(in>histos[mapnamePredm]->GetNbinsX()){
				cout << " &                    &                                &         ";
			} else {
				if(in==histos[mapnamePredm]->GetNbinsX())
					cout << " & $\\geq"<< (int)histos[mapnamePredm]->GetBinLowEdge(in) << "$ ";
				else 
				cout << " & $"<< (int)histos[mapnamePredm]->GetBinLowEdge(in)<< "-" << (int)histos[mapnamePredm]->GetBinLowEdge(in+1) << "$ ";
				cout << "& " << histos[mapnamePredm]->GetBinContent(in) << "$\\pm$" << histos[mapnamePredm]->GetBinError(in) << " ";
				cout << "& " << (int)histos[mapnameDatam]->GetBinContent(in) << " ";
			}
			//high HT
			if(in>histos[mapnamePredh]->GetNbinsX()){
				cout << " &                    &                                &         ";
			} else {
				if(in==histos[mapnamePredh]->GetNbinsX())
					cout << " & $\\geq"<< (int)histos[mapnamePredh]->GetBinLowEdge(in) << "$ ";
				else 
				cout << " & $"<< (int)histos[mapnamePredh]->GetBinLowEdge(in)<< "-" << (int)histos[mapnamePredh]->GetBinLowEdge(in+1) << "$ ";
				cout << "& " << histos[mapnamePredh]->GetBinContent(in) << "$\\pm$" << histos[mapnamePredh]->GetBinError(in) << " ";
				cout << "& " << (int)histos[mapnameDatah]->GetBinContent(in) << " ";
			}
		cout << "\\\\" << endl;
		}
		cout << "\\hline" << endl;
	}//signalregionsize
	cout << "\\hline" << endl;	
	cout << "\\end{tabular}" << endl
	     << "\\end{center}"  << endl
	     << "\\end{table}"   << endl
	     << "\%ENDLATEX\%"   << endl
	     << endl;

	TLatex *TitleBox = new TLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = 19.5 fb^{-1}");
	TitleBox->SetNDC();
	TitleBox->SetTextFont(42);
	TitleBox->SetTextSize(0.04181185);
	TitleBox->SetLineWidth(2);
	TLatex *   TitleBoxS = new TLatex(0.328418,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = 19.5 fb^{-1}");
	TitleBoxS->SetNDC();
	TitleBoxS->SetTextFont(42);
	TitleBoxS->SetTextSize(0.04181185);
	TitleBoxS->SetLineWidth(2);
	//jet
	TLatex box1;//0.64,0.8548951
	box1.SetNDC();
   	box1.SetTextAlign(31);
   	box1.SetTextFont(42);
   	box1.SetTextSize(0.04181185);
   	box1.SetLineWidth(2);
	TString text1 = "";
	//b jet
	TLatex box2;//0.64,0.791958
	box2.SetNDC();
   	box2.SetTextAlign(31);
   	box2.SetTextFont(42);
   	box2.SetTextSize(0.04181185);
   	box2.SetLineWidth(2);
	TString text2 = "";
	//HT
	TLatex box3;//0.64,0.7290209
	box3.SetNDC();
   	box3.SetTextAlign(31);
   	box3.SetTextFont(42);
   	box3.SetTextSize(0.04181185);
   	box3.SetLineWidth(2);
	TString text3 = "";

	TLatex scale2j0b;
	scale2j0b.SetNDC();
	scale2j0b.SetTextFont(42);
   	scale2j0b.SetTextSize(0.04181185);
	scale2j0b.SetLineWidth(2);

	TLatex HTbox;
	HTbox.SetTextFont(42);
   	HTbox.SetTextSize(0.04181185);
	HTbox.SetLineWidth(2);


	//all plots but grand summary plot done with this canvas
	TCanvas *c1 = new TCanvas("c1", "c1",485,220,600,600);
   	gStyle->SetOptFit(1);						gStyle->SetOptStat(0);		gStyle->SetOptTitle(0);
	c1->Range(82.71719,-0.4425771,532.9945,2.212885);
	c1->SetFillColor(0);		c1->SetBorderMode(0);		c1->SetBorderSize(2);		c1->SetTickx(1);		c1->SetTicky(1);
	c1->SetLeftMargin(0.18);	c1->SetRightMargin(0.05);	c1->SetTopMargin(0.07);		c1->SetBottomMargin(0.15);
	c1->SetFrameFillStyle(0);	c1->SetFrameBorderMode(0);	c1->SetFrameFillStyle(0);	c1->SetFrameBorderMode(0);

	//prediction plots (the real ones)
	for(int i3 = 0; i3<HTregionsize;     ++i3){
	for(int i2 = 0; i2<signalregionsize; ++i2){
		c1->Clear();
		c1->cd();
		string hs = string("_") + HT_bin[i3] + string("_") + signal_region[i2];
		string mapnameP = "MT2PredStack" + hs;
		string mapnameS = "MT2PredSum" + hs;
		string mapnameD = "MT2Data" + hs;
		double max1 = histos[mapnameS]->GetMaximum();
		double max2 = histos[mapnameD]->GetMaximum();
		double max  = (max1>max2)?max1:max2;
		max = 1.3*max;
		stackPrediction[mapnameP]->SetMaximum(max);
		stackPrediction[mapnameP]->SetMinimum(0.);

		stackPrediction[mapnameP]->Draw("hist");
		TH1D *haxis = (TH1D*)stackPrediction[mapnameP]->GetHistogram();
		haxis->GetXaxis()->SetLabelSize(0.05);
		haxis->GetXaxis()->SetLabelOffset(0.007);
		haxis->GetXaxis()->SetLabelFont(42);
		haxis->GetXaxis()->SetTitleFont(42);
		haxis->GetXaxis()->SetTitleSize(0.06);
		haxis->GetXaxis()->SetTitleOffset(0.9);
		haxis->GetYaxis()->SetLabelSize(0.05);
		haxis->GetYaxis()->SetLabelOffset(0.007);
		haxis->GetYaxis()->SetLabelFont(42);
		haxis->GetYaxis()->SetTitleFont(42);
		haxis->GetYaxis()->SetTitleSize(0.06);
		haxis->GetYaxis()->SetTitleOffset(1.25);
		stackPrediction[mapnameP]->SetHistogram(haxis);
		stackPrediction[mapnameP]->Draw("hist");
		if(!poisson) histos[mapnameD]->Draw("sameE1");//from next ROOT version this should also work for poissonian uncertainties
		else graphData[mapnameD]->Draw("P");
		histos[mapnameS]->Draw("sameE2");
		leg->Draw();
		if(i2==0) {text1 = "2 jets";     text2 = "0 b jets";}
		if(i2==1) {text1 = "2 jets";     text2 = "#geq1 b jets";}
		if(i2==2) {text1 = "3-5 jets";   text2 = "0 b jets";}
		if(i2==3) {text1 = "3-5 jets";   text2 = "1 b jets";}
		if(i2==4) {text1 = "3-5 jets";   text2 = "2 b jets";}
		if(i2==5) {text1 = "#geq6 jets"; text2 = "0 b jets";}
		if(i2==6) {text1 = "#geq6 jets"; text2 = "1 b jets";}
		if(i2==7) {text1 = "#geq6 jets"; text2 = "2 b jets";}
		if(i2==8) {text1 = "#geq3 jets"; text2 = "#geq3 b jets";}
		if(i3==0) text3 = "low H_{T}";
		if(i3==1) text3 = "medium H_{T}";
		if(i3==2) text3 = "high H_{T}";
		box1.DrawLatex(0.64,0.8548951,text1.Data();
		box2.DrawLatex(0.64,0.791958 ,text2.Data();
		box2.DrawLatex(0.64,0.7290209,text2.Data();
		TitleBox->Draw();
		string outname = histos[mapnameS]->GetName();
		string cname = outputdir + outname + ".C";
		outname = outputdir + outname + ".eps";
		if(saveeps) c1->SaveAs(outname.c_str());
		if(saveasC) c1->SaveAs(cname.c_str());
	}}

	//pull plots
	TLatex pullbox2;
	pullbox2.SetNDC();
   	pullbox2.SetTextAlign(12);
   	pullbox2.SetTextFont(42);
   	pullbox2.SetTextSize(0.04181185);
   	pullbox2.SetLineWidth(2);
	TLatex pullbox1;
	pullbox1.SetNDC();
   	pullbox1.SetTextAlign(12);
   	pullbox1.SetTextFont(42);
   	pullbox1.SetTextSize(0.04181185);
   	pullbox1.SetLineWidth(2);
	TString text3 = "";
	if(useMCyields){
		c1->Clear();
		c1->cd();
		hPullMCDa->GetXaxis()->SetTitle("Pull");
		hPullMCDa->GetYaxis()->SetTitle("Entries");
		hPullMCDa->GetXaxis()->SetLabelSize(0.05);
		hPullMCDa->GetXaxis()->SetLabelOffset(0.007);
		hPullMCDa->GetXaxis()->SetLabelFont(42);
		hPullMCDa->GetXaxis()->SetTitleFont(42);
		hPullMCDa->GetXaxis()->SetTitleSize(0.06);
		hPullMCDa->GetXaxis()->SetTitleOffset(0.9);
		hPullMCDa->GetYaxis()->SetLabelSize(0.05);
		hPullMCDa->GetYaxis()->SetLabelOffset(0.007);
		hPullMCDa->GetYaxis()->SetLabelFont(42);
		hPullMCDa->GetYaxis()->SetTitleFont(42);
		hPullMCDa->GetYaxis()->SetTitleSize(0.06);
		hPullMCDa->GetYaxis()->SetTitleOffset(1.25);
		hPullMCDa->Draw("hist");
		stringstream pullstrDa;
		pullstrDa << "Mean: " << fixed << setprecision(2) << hPullMCDa->GetMean();
		text1 = pullstrDa.str().c_str();
		stringstream pullstrDa2;
		pullstrDa2 << "RMS: " << fixed << setprecision(2) << hPullMCDa->GetRMS();
		text2 = pullstrDa2.str().c_str();
		TitleBox->Draw();
		pullbox1.Draw(0.7181208,0.8776224,text1.Data());
		pullbox2.Draw(0.7181208,0.8286713,text2.Data());
		string outname = hPullMCDa->GetName();
		string cname = outputdir + outname + ".C";
		outname = outputdir + outname + ".eps";
		if(saveeps) c1->SaveAs(outname.c_str());
		if(saveasC) c1->SaveAs(cname.c_str());
		cout << hPullMCDa->GetName() << " has mean " << hPullMCDa->GetMean() << " +/- " << hPullMCDa->GetMeanError() << " and rms " << hPullMCDa->GetRMS() << " +/- " << hPullMCDa->GetRMSError() << endl;
	
		c1->Clear();
		c1->cd();
		hPullMC->GetXaxis()->SetTitle("Pull");
		hPullMC->GetYaxis()->SetTitle("Entries");
		hPullMC->GetXaxis()->SetLabelSize(0.05);
		hPullMC->GetXaxis()->SetLabelOffset(0.007);
		hPullMC->GetXaxis()->SetLabelFont(42);
		hPullMC->GetXaxis()->SetTitleFont(42);
		hPullMC->GetXaxis()->SetTitleSize(0.06);
		hPullMC->GetXaxis()->SetTitleOffset(0.9);
		hPullMC->GetYaxis()->SetLabelSize(0.05);
		hPullMC->GetYaxis()->SetLabelOffset(0.007);
		hPullMC->GetYaxis()->SetLabelFont(42);
		hPullMC->GetYaxis()->SetTitleFont(42);
		hPullMC->GetYaxis()->SetTitleSize(0.06);
		hPullMC->GetYaxis()->SetTitleOffset(1.25);
		hPullMC->Draw("hist");
		stringstream pullstrMC;
		pullstrMC << "Mean: " << fixed << setprecision(2) << hPullMCDa->GetMean();
		text1 = pullstrMC.str().c_str();
		stringstream pullstrMC2;
		pullstrMC2 << "RMS: " << fixed << setprecision(2) << hPullMCDa->GetRMS();
		text2 = pullstrMC2.str().c_str();
		TitleBox->Draw();
		pullbox1.Draw(0.7181208,0.8776224,text1.Data());
		pullbox2.Draw(0.7181208,0.8286713,text2.Data());
		outname = hPullMC->GetName();
		cname = outputdir + outname + ".C";
		outname = outputdir + outname + ".eps";
		if(saveeps) c1->SaveAs(outname.c_str());	
		if(saveasC) c1->SaveAs(cname.c_str());	
		cout << hPullMC->GetName() << " has mean " << hPullMC->GetMean() << " +/- " << hPullMC->GetMeanError() << " and rms " << hPullMC->GetRMS() << " +/- " << hPullMC->GetRMSError() << endl;

	}
	c1->Clear();
	c1->cd();
	hPullData->GetXaxis()->SetTitle("Pull");
	hPullData->GetYaxis()->SetTitle("Entries");
	hPullData->GetXaxis()->SetLabelSize(0.05);
	hPullData->GetXaxis()->SetLabelOffset(0.007);
	hPullData->GetXaxis()->SetLabelFont(42);
	hPullData->GetXaxis()->SetTitleFont(42);
	hPullData->GetXaxis()->SetTitleSize(0.06);
	hPullData->GetXaxis()->SetTitleOffset(0.9);
	hPullData->GetYaxis()->SetLabelSize(0.05);
	hPullData->GetYaxis()->SetLabelOffset(0.007);
	hPullData->GetYaxis()->SetLabelFont(42);
	hPullData->GetYaxis()->SetTitleFont(42);
	hPullData->GetYaxis()->SetTitleSize(0.06);
	hPullData->GetYaxis()->SetTitleOffset(1.25);
	hPullData->Draw("hist");
	stringstream pullstr;
	pullstr << "Mean: " << fixed << setprecision(2) << hPullMCDa->GetMean();
	text1 = pullstr.str().c_str();
	stringstream pullstr2;
	pullstr2 << "RMS: " << fixed << setprecision(2) << hPullMCDa->GetRMS();
	text2 = pullstr2.str().c_str();
	TitleBox->Draw();
	pullbox1.Draw(0.7181208,0.8776224,text1.Data());
	pullbox2.Draw(0.7181208,0.8286713,text2.Data());
	string outname = hPullData->GetName();
	string cname = outputdir + outname + ".C";
	outname = outputdir + outname + ".eps";
	if(saveeps) c1->SaveAs(outname.c_str());
	if(saveasC) c1->SaveAs(cname.c_str());
	cout << hPullData->GetName() << " has mean " << hPullData->GetMean() << " +/- " << hPullData->GetMeanError() << " and rms " << hPullData->GetRMS() << " +/- " << hPullData->GetRMSError() << endl;

	double max1,max2,max;
	TH1D *haxis;

	//summary plots - the three ones
	c1->Clear();
	c1->cd();
	max1 = predSumlHT->GetMaximum();
	max2 = predDatalHT->GetMaximum();
	max  = (max1>max2)?max1:max2;
	max = 1.3*max;
	SumlHT->SetMaximum(max);
	SumlHT->SetMinimum(0.);
	SumlHT->Draw("hist");
	haxis->GetXaxis()->SetLabelSize(0.05);
	haxis->GetXaxis()->SetLabelOffset(0.007);
	haxis->GetXaxis()->SetLabelFont(42);
	haxis->GetXaxis()->SetTitleFont(42);
	haxis->GetXaxis()->SetTitleSize(0.06);
	haxis->GetXaxis()->SetTitleOffset(0.9);
	haxis->GetYaxis()->SetLabelSize(0.05);
	haxis->GetYaxis()->SetLabelOffset(0.007);
	haxis->GetYaxis()->SetLabelFont(42);
	haxis->GetYaxis()->SetTitleFont(42);
	haxis->GetYaxis()->SetTitleSize(0.06);
	haxis->GetYaxis()->SetTitleOffset(1.25);
	SumlHT->SetHistogram(haxis);
	SumlHT->Draw("hist");
	if(!poisson) predDatalHT->Draw("sameE1");//from next ROOT version this should also work for poissonian uncertainties
	else gpredDatalHT->Draw("P");
	predSumlHT->Draw("sameE2");
	leg->Draw();
	text1 = "low H_{T}";
	box1.DrawLatex(0.64,0.8548951,text1.Data());
	TitleBox->Draw();
	if(summaryscale2bins){
		scale2j0b.DrawLatex(0.20,.80,"x#frac{1}{6}");
		scale2j0b.DrawLatex(0.37,.80,"x#frac{1}{9}");
		scale2j0b.DrawLatex(0.46,.80,"x#frac{1}{3}");
	}
	outname = predDatalHT->GetName();
	cname = outputdir + outname + ".C";
	outname = outputdir + outname + ".eps";
	if(saveeps) c1->SaveAs(outname.c_str());
	if(saveasC) c1->SaveAs(cname.c_str());

	c1->Clear();
	c1->cd();
	max1 = predSummHT->GetMaximum();
	max2 = predDatamHT->GetMaximum();
	max  = (max1>max2)?max1:max2;
	max = 1.3*max;
	SummHT->SetMaximum(max);
	SummHT->SetMinimum(0.);
	SummHT->Draw("hist");
	haxis = (TH1D*)SummHT->GetHistogram();
	haxis->GetXaxis()->SetLabelSize(0.05);
	haxis->GetXaxis()->SetLabelOffset(0.007);
	haxis->GetXaxis()->SetLabelFont(42);
	haxis->GetXaxis()->SetTitleFont(42);
	haxis->GetXaxis()->SetTitleSize(0.06);
	haxis->GetXaxis()->SetTitleOffset(0.9);
	haxis->GetYaxis()->SetLabelSize(0.05);
	haxis->GetYaxis()->SetLabelOffset(0.007);
	haxis->GetYaxis()->SetLabelFont(42);
	haxis->GetYaxis()->SetTitleFont(42);
	haxis->GetYaxis()->SetTitleSize(0.06);
	haxis->GetYaxis()->SetTitleOffset(1.25);
	SummHT->SetHistogram(haxis);
	SummHT->Draw("hist");
	if(!poisson) predDatamHT->Draw("sameE1");//from next ROOT version this should also work for poissonian uncertainties
	else gpredDatamHT->Draw("P");
	predSummHT->Draw("sameE2");
	leg->Draw();
	text1 = "medium H_{T}";
	box1.DrawLatex(0.64,0.8548951,text1.Data());
	TitleBox->Draw();
	if(summaryscale2bins){
		scale2j0b.DrawLatex(0.20,.80,"x#frac{1}{4}");
		scale2j0b.DrawLatex(0.37,.80,"x#frac{1}{6}");
		scale2j0b.DrawLatex(0.46,.80,"x#frac{1}{2}");
	}
	outname = predDatamHT->GetName();
	cname = outputdir + outname + ".C";
	outname = outputdir + outname + ".eps";
	if(saveeps) c1->SaveAs(outname.c_str());
	if(saveasC) c1->SaveAs(cname.c_str());

	c1->Clear();
	c1->cd();
	max1 = predSumhHT->GetMaximum();
	max2 = predDatahHT->GetMaximum();
	max  = (max1>max2)?max1:max2;
	max = 1.3*max;
	SumhHT->SetMaximum(max);
	SumhHT->SetMinimum(0.);
	SumhHT->Draw("hist");
	haxis = (TH1D*)SumhHT->GetHistogram();
	haxis->GetXaxis()->SetLabelSize(0.05);
	haxis->GetXaxis()->SetLabelOffset(0.007);
	haxis->GetXaxis()->SetLabelFont(42);
	haxis->GetXaxis()->SetTitleFont(42);
	haxis->GetXaxis()->SetTitleSize(0.06);
	haxis->GetXaxis()->SetTitleOffset(0.9);
	haxis->GetYaxis()->SetLabelSize(0.05);
	haxis->GetYaxis()->SetLabelOffset(0.007);
	haxis->GetYaxis()->SetLabelFont(42);
	haxis->GetYaxis()->SetTitleFont(42);
	haxis->GetYaxis()->SetTitleSize(0.06);
	haxis->GetYaxis()->SetTitleOffset(1.25);
	SumhHT->SetHistogram(haxis);
	SumhHT->Draw("hist");
	if(!poisson) predDatahHT->Draw("sameE1");//from next ROOT version this should also work for poissonian uncertainties
	else gpredDatahHT->Draw("P");
	predSumhHT->Draw("sameE2");
	leg->Draw();
	text1 = "medium H_{T}";
	box1.DrawLatex(0.64,0.8548951,text1.Data());
	TitleBox->Draw();
	if(summaryscale2bins){
		scale2j0b.DrawLatex(0.20,.80,"x#frac{1}{2}");
		scale2j0b.DrawLatex(0.37,.80,"x#frac{1}{3}");
	}
	outname = predDatahHT->GetName();
	cname = outputdir + outname + ".C";
	outname = outputdir + outname + ".eps";
	if(saveeps) c1->SaveAs(outname.c_str());
	if(saveasC) c1->SaveAs(cname.c_str());

    //redefine all latex here
    //summary plots - the big one
	TCanvas *c2 = new TCanvas("c2", "c2",56,22,870,580);
	gStyle->SetOptFit(1);	gStyle->SetOptStat(0);	gStyle->SetOptTitle(0);
	c2->Range(-6.311689,-0.7369661,28.75325,4.476528);
	c2->SetFillColor(0);		c2->SetBorderMode(0);		c2->SetBorderSize(2);		c2->SetLogy();	c2->SetTickx(1);	c2->SetTicky(1);
	c2->SetLeftMargin(0.18);	c2->SetRightMargin(0.05);	c2->SetTopMargin(0.07);		c2->SetBottomMargin(0.3);
	c2->SetFrameFillStyle(0);	c2->SetFrameBorderMode(0);	c2->SetFrameFillStyle(0);	c2->SetFrameBorderMode(0);

	c2->Clear();
	c2->cd();
	gPad->SetLogy();
	max1 = predSum->GetMaximum();
	max2 = predData->GetMaximum();
	max  = (max1>max2)?max1:max2;
	max = 2.6*max;
	Sum->SetMaximum(max);
	Sum->SetMinimum(2.5);
	Sum->Draw("hist");
	haxis = (TH1D*)Sum->GetHistogram();
	haxis->GetXaxis()->SetBit(TAxis::kLabelsVert);
	haxis->GetXaxis()->SetLabelFont(42);
	haxis->GetXaxis()->SetLabelOffset(0.007);
	haxis->GetXaxis()->SetLabelSize(0.05);
	haxis->GetXaxis()->SetTitleSize(0.06);
	haxis->GetXaxis()->SetTitleOffset(2.6);
	haxis->GetXaxis()->SetTitleFont(42);
	haxis->GetYaxis()->SetLabelFont(42);
	haxis->GetYaxis()->SetLabelOffset(0.007);
	haxis->GetYaxis()->SetLabelSize(0.05);
	haxis->GetYaxis()->SetTitleSize(0.06);
	haxis->GetYaxis()->SetTitleOffset(1.25);
	haxis->GetYaxis()->SetTitleFont(42);
	haxis->GetXaxis()->LabelsOption("v");
	Sum->SetHistogram(haxis);
    	Sum->SetMinimum(3.);
	Sum->Draw("hist");
	if(!poisson) predData->Draw("sameE1");//from next ROOT version this should also work for poissonian uncertainties
	else gpredData->Draw("P");
	predSum->Draw("sameE2");
	leg->Draw();
	TitleBoxS->Draw();
	TArrow *arrow = new TArrow(0.05,0.0725,8.95,0.0725,0.02,"<>");
	arrow->SetFillColor(1); arrow->SetFillStyle(1001);  arrow->SetLineStyle(3); arrow->SetLineWidth(2);
	arrow->Draw();
	arrow = new TArrow(9.05,0.0725,17.95,0.0725,0.02,"<>");
	arrow->SetFillColor(1); arrow->SetFillStyle(1001);  arrow->SetLineStyle(3); arrow->SetLineWidth(2);
	arrow->Draw();
	arrow = new TArrow(18.05,0.0725,26.95,0.0725,0.02,"<>");
	arrow->SetFillColor(1); arrow->SetFillStyle(1001);  arrow->SetLineStyle(3); arrow->SetLineWidth(2);
	arrow->Draw();
	TLine *line = new TLine(9.0,0.06,9.0,13617.89);
	line->SetLineStyle(2);
	line->Draw();
	line = new TLine(18.,0.06,18.,13617.89);
	line->SetLineStyle(2);
	line->Draw();
	line = new TLine(0.0,2.5,0.0,0.06);
	line->SetLineStyle(2);
	line->Draw();
	line = new TLine(27,2.5,27,0.06);
	line->SetLineStyle(2);
	line->Draw();
	HTBox.DrawLatex(2.9,0.094,"low H_{T}");
	HTBox.DrawLatex(10.5,0.094,"medium H_{T}");
	HTBox.DrawLatex(20.75,0.094,"high H_{T}");
	outname = predData->GetName();outname = outputdir + outname + ".C";
	if(saveasC) c2->SaveAs(outname.c_str());
	outname = predData->GetName();outname = outputdir + outname + ".eps";
	if(saveeps) c2->SaveAs(outname.c_str());
    
    
    //ratio plot - the big one
	TCanvas *c3 = new TCanvas("c3", "c3",67,27,870,580);
	gStyle->SetOptFit(1);	gStyle->SetOptStat(0);	gStyle->SetOptTitle(0);
	c3->Range(-6.311689,-0.4807692,28.75325,2.724359);
	c3->SetFillColor(0);		c3->SetBorderMode(0);		c3->SetBorderSize(2);		c3->SetTickx(1);	c3->SetTicky(1);
	c3->SetLeftMargin(0.18);	c3->SetRightMargin(0.05);	c3->SetTopMargin(0.07);		c3->SetBottomMargin(0.3);
	c3->SetFrameFillStyle(0);	c3->SetFrameBorderMode(0);	c3->SetFrameFillStyle(0);	c3->SetFrameBorderMode(0);
	c3->Clear();
	c3->cd();
	gPad->SetLogy(0);

	predRatioR->SetMaximum(2.);
	predRatioR->SetMinimum(0.);
	predRatioR->Draw("hist");
	predRatioR->GetXaxis()->SetBit(TAxis::kLabelsVert);
	predRatioR->GetXaxis()->SetLabelFont(42);
	predRatioR->GetXaxis()->SetLabelOffset(0.007);
	predRatioR->GetXaxis()->SetLabelSize(0.05);
	predRatioR->GetXaxis()->SetTitleSize(0.06);
	predRatioR->GetXaxis()->SetTitleOffset(2.6);
	predRatioR->GetXaxis()->SetTitleFont(42);
	predRatioR->GetYaxis()->SetTitle("ratio of data / prediction");
	predRatioR->GetYaxis()->SetLabelFont(42);
	predRatioR->GetYaxis()->SetLabelOffset(0.007);
	predRatioR->GetYaxis()->SetLabelSize(0.05);
	predRatioR->GetYaxis()->SetTitleSize(0.06);
	predRatioR->GetYaxis()->SetTitleOffset(1.25);
	predRatioR->GetYaxis()->SetTitleFont(42);
	predRatioR->GetZaxis()->SetLabelFont(42);
	predRatioR->GetZaxis()->SetLabelOffset(0.007);
	predRatioR->GetZaxis()->SetLabelSize(0.05);
	predRatioR->GetZaxis()->SetTitleSize(0.06);
	predRatioR->GetZaxis()->SetTitleFont(42);
	predRatioR->GetXaxis()->LabelsOption("v");
	predRatioR->SetMinimum(0.);
	predRatioR->SetMaximum(2.5);
	predRatioR->Draw("E2");
	line = new TLine(0.0,1.0,27.0,1.0);
	line->SetLineWidth(2);
	line->SetLineColor(kRed);
	line->Draw();
	predRatio->Draw("sameE1");
	TitleBoxS->Draw();
	arrow = new TArrow(0.05,-0.73,8.95,-0.73,0.02,"<>");
	arrow->SetFillColor(1); arrow->SetFillStyle(1001);  arrow->SetLineStyle(3); arrow->SetLineWidth(2);
	arrow->Draw();
	arrow = new TArrow(9.05,-0.73,17.95,-0.73,0.02,"<>");
	arrow->SetFillColor(1); arrow->SetFillStyle(1001);  arrow->SetLineStyle(3); arrow->SetLineWidth(2);
	arrow->Draw();
	arrow = new TArrow(18.05,-0.73,26.95,-0.73,0.02,"<>");
	arrow->SetFillColor(1); arrow->SetFillStyle(1001);  arrow->SetLineStyle(3); arrow->SetLineWidth(2);
	arrow->Draw();
	line = new TLine(9.0,-0.76,9.0,2.5);
	line->SetLineStyle(2);
	line->Draw();
	line = new TLine(18.,-0.76,18.,2.5);
	line->SetLineStyle(2);
	line->Draw();
	line = new TLine(0.0,-0.76,0.0,0.0);
	line->SetLineStyle(2);
	line->Draw();
	line = new TLine(27,-0.76,27,0.0);
	line->SetLineStyle(2);
	line->Draw();
	HTBox.DrawLatex(2.9,-0.66,"low H_{T}");
	HTBox.DrawLatex(11.25,-0.66,"medium H_{T}");
	HTBox.DrawLatex(20.75,-0.66,"high H_{T}");
	outname = predRatio->GetName();outname = outputdir + outname + ".C";
	if(saveasC) c3->SaveAs(outname.c_str());
	outname = predRatio->GetName();outname = outputdir + outname + ".eps";
	if(saveeps) c3->SaveAs(outname.c_str());


}//void MT2Results_PlotsAndTables()