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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"//use your own path

//run via root -l -b -q ZinvVisualization.C++

//this makes Z(nunu) prediction vs. MC truth plots for the AN.
//at the moment all numbers are hard-coded and this code needs major revision, the mc truth numbers are loaded from a root file created by ZnunuNumbers.C
//furthermore the stored root file is used for both producing the final result plots/tables and also interpretation
void ZinvVisualization(){

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

const int signalregionsize = 9;
string signal_region[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b"};
const int HTbinsize = 3;
string HT_bin[HTbinsize] = {"lowHT", "mediumHT", "highHT"};

    	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

    TLegend *leg = new TLegend(0.6551724,0.7299578,0.8706897,0.8987342,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(62);
    leg->SetTextSize(0.04575163);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);

//load MC truth numbers
TFile *fmctruth = TFile::Open("../Results/Filtered/GammaJetsPrediction/20130617_test/ZnunuNumbers.root");

	//define all histograms
	map<string, TH1D*>    histos;
	for(int i3 = 0; i3<HTbinsize;        ++i3){
	for(int i2 = 0; i2<signalregionsize; ++i2){
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

		string hs = string("_") + HT_bin[i3] + string("_") + signal_region[i2];// + string("_") + sample_type[i1];
		string mapname = "MT2" + hs;
		TH1D *h = (TH1D*)fmctruth->Get(mapname.c_str());
		mapname = "MT2truth" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)h->Clone(mapname.c_str());
		histos[mapname]->SetFillColor(kViolet-3); histos[mapname]->SetFillStyle(3001);
		//get MC truth histo style
		if(i3==0&&i2==0) leg->AddEntry(histos[mapname], "MC truth", "f");
			histos[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]");
			histos[mapname]->GetXaxis()->SetLabelFont(42);
			histos[mapname]->GetXaxis()->SetLabelSize(0.06);
			histos[mapname]->GetXaxis()->SetTitleSize(0.06);
			histos[mapname]->GetXaxis()->SetLabelOffset(0.01);
			histos[mapname]->GetXaxis()->SetTitleOffset(1.2);
			histos[mapname]->GetXaxis()->SetTitleFont(42);
			histos[mapname]->GetYaxis()->SetTitle("Z(#nu#nu) yield");
			histos[mapname]->GetYaxis()->SetLabelFont(42);
			histos[mapname]->GetYaxis()->SetLabelSize(0.06);
			histos[mapname]->GetYaxis()->SetTitleSize(0.06);
			histos[mapname]->GetYaxis()->SetLabelOffset(0.01);
			histos[mapname]->GetYaxis()->SetTitleOffset(1.2);
			histos[mapname]->GetYaxis()->SetTitleFont(42);
			histos[mapname]->GetZaxis()->SetLabelFont(42);
			histos[mapname]->GetZaxis()->SetLabelSize(0.035);
			histos[mapname]->GetZaxis()->SetTitleSize(0.035);
			histos[mapname]->GetZaxis()->SetTitleFont(42);
		mapname = "MT2pred" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", NMT2bins, MT2bins);
		histos[mapname]->SetMarkerStyle(20), histos[mapname]->SetMarkerColor(kBlack); histos[mapname]->SetLineWidth(3); histos[mapname]->SetLineColor(kBlack);
		//get prediction histo style
		if(i3==0&&i2==0) leg->AddEntry(histos[mapname], "data prediction", "lp");
			histos[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]");
			histos[mapname]->GetXaxis()->SetLabelFont(42);
			histos[mapname]->GetXaxis()->SetLabelSize(0.06);
			histos[mapname]->GetXaxis()->SetTitleSize(0.06);
			histos[mapname]->GetXaxis()->SetLabelOffset(0.01);
			histos[mapname]->GetXaxis()->SetTitleOffset(1.2);
			histos[mapname]->GetXaxis()->SetTitleFont(42);
			histos[mapname]->GetYaxis()->SetTitle("Z(#nu#nu) yield");
			histos[mapname]->GetYaxis()->SetLabelFont(42);
			histos[mapname]->GetYaxis()->SetLabelSize(0.06);
			histos[mapname]->GetYaxis()->SetTitleSize(0.06);
			histos[mapname]->GetYaxis()->SetLabelOffset(0.01);
			histos[mapname]->GetYaxis()->SetTitleOffset(1.2);
			histos[mapname]->GetYaxis()->SetTitleFont(42);
			histos[mapname]->GetZaxis()->SetLabelFont(42);
			histos[mapname]->GetZaxis()->SetLabelSize(0.035);
			histos[mapname]->GetZaxis()->SetTitleSize(0.035);
			histos[mapname]->GetZaxis()->SetTitleFont(42);
	}}

	//hard-coded numbers
	for(int i3 = 0; i3<HTbinsize;        ++i3){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	   string hs = string("_") + HT_bin[i3] + string("_") + signal_region[i2];
	   string mapname = "MT2pred" + hs;
	   if(i3==0&&i2==0){//lHT,2j,0b
		histos[mapname]->SetBinContent(1,304.27); histos[mapname]->SetBinError(1,sqrt(pow(11.84,2)+pow(62.55,2)));
		histos[mapname]->SetBinContent(2,236.52); histos[mapname]->SetBinError(2,sqrt(pow(10.55,2)+pow(48.41,2)));
		histos[mapname]->SetBinContent(3,182.98); histos[mapname]->SetBinError(3,sqrt(pow( 9.30,2)+pow(37.50,2)));
		histos[mapname]->SetBinContent(4,160.92); histos[mapname]->SetBinError(4,sqrt(pow( 9.09,2)+pow(50.38,2)));
		histos[mapname]->SetBinContent(5,116.69); histos[mapname]->SetBinError(5,sqrt(pow( 7.46,2)+pow(35.07,2)));
		histos[mapname]->SetBinContent(6, 47.95); histos[mapname]->SetBinError(6,sqrt(pow( 4.86,2)+pow(14.57,2)));
		histos[mapname]->SetBinContent(7, 11.74); histos[mapname]->SetBinError(7,sqrt(pow( 2.39,2)+pow( 3.56,2)));
		histos[mapname]->SetBinContent(8,  3.28); histos[mapname]->SetBinError(8,sqrt(pow( 1.24,2)+pow( 0.98,2)));
	   } if(i3==0&&i2==2){//lHT,35j,0b
		histos[mapname]->SetBinContent(1,458.79); histos[mapname]->SetBinError(1,sqrt(pow(16.03,2)+pow(97.64,2)));
		histos[mapname]->SetBinContent(2,366.96); histos[mapname]->SetBinError(2,sqrt(pow(14.70,2)+pow(80.26,2)));
		histos[mapname]->SetBinContent(3,301.04); histos[mapname]->SetBinError(3,sqrt(pow(12.94,2)+pow(62.16,2)));
		histos[mapname]->SetBinContent(4,173.25); histos[mapname]->SetBinError(4,sqrt(pow(10.18,2)+pow(54.86,2)));
		histos[mapname]->SetBinContent(5, 93.60); histos[mapname]->SetBinError(5,sqrt(pow( 7.29,2)+pow(28.28,2)));
		histos[mapname]->SetBinContent(6, 41.36); histos[mapname]->SetBinError(6,sqrt(pow( 4.88,2)+pow(12.64,2)));
		histos[mapname]->SetBinContent(7,  9.80); histos[mapname]->SetBinError(7,sqrt(pow( 2.38,2)+pow( 3.01,2)));
		histos[mapname]->SetBinContent(8,  2.73); histos[mapname]->SetBinError(8,sqrt(pow( 1.22,2)+pow( 0.82,2)));
	   } if(i3==0&&i2==5){//lHT,6j,0b
		histos[mapname]->SetBinContent(1, 12.02); histos[mapname]->SetBinError(1,sqrt(pow( 2.76,2)+pow( 2.58,2)));
		histos[mapname]->SetBinContent(2,  5.13); histos[mapname]->SetBinError(2,sqrt(pow( 1.95,2)+pow( 1.29,2)));
		histos[mapname]->SetBinContent(3,  3.53); histos[mapname]->SetBinError(3,sqrt(pow( 1.77,2)+pow( 1.24,2)));
	   } if(i3==1&&i2==0){//mHT,2j,0b
		histos[mapname]->SetBinContent(1, 76.13); histos[mapname]->SetBinError(1,sqrt(pow( 5.38,2)+pow(15.32,2)));
		histos[mapname]->SetBinContent(2, 62.24); histos[mapname]->SetBinError(2,sqrt(pow( 4.92,2)+pow(12.55,2)));
		histos[mapname]->SetBinContent(3, 38.49); histos[mapname]->SetBinError(3,sqrt(pow( 3.95,2)+pow( 7.78,2)));
		histos[mapname]->SetBinContent(4, 40.84); histos[mapname]->SetBinError(4,sqrt(pow( 4.19,2)+pow( 8.26,2)));
		histos[mapname]->SetBinContent(5, 20.11); histos[mapname]->SetBinError(5,sqrt(pow( 3.00,2)+pow( 4.08,2)));
		histos[mapname]->SetBinContent(6, 30.21); histos[mapname]->SetBinError(6,sqrt(pow( 3.82,2)+pow( 9.14,2)));
		histos[mapname]->SetBinContent(7, 12.73); histos[mapname]->SetBinError(7,sqrt(pow( 2.44,2)+pow( 3.84,2)));
		histos[mapname]->SetBinContent(8, 15.96); histos[mapname]->SetBinError(8,sqrt(pow( 2.79,2)+pow( 4.83,2)));
		histos[mapname]->SetBinContent(9,  2.83); histos[mapname]->SetBinError(9,sqrt(pow( 1.15,2)+pow( 0.85,2)));
	   } if(i3==1&&i2==2){//mHT,35j,0b
		histos[mapname]->SetBinContent(1, 88.17); histos[mapname]->SetBinError(1,sqrt(pow( 6.50,2)+pow(17.75,2)));
		histos[mapname]->SetBinContent(2, 73.11); histos[mapname]->SetBinError(2,sqrt(pow( 6.38,2)+pow(14.78,2)));
		histos[mapname]->SetBinContent(3, 63.34); histos[mapname]->SetBinError(3,sqrt(pow( 5.72,2)+pow(12.78,2)));
		histos[mapname]->SetBinContent(4, 59.03); histos[mapname]->SetBinError(4,sqrt(pow( 5.71,2)+pow(11.95,2)));
		histos[mapname]->SetBinContent(5, 50.85); histos[mapname]->SetBinError(5,sqrt(pow( 5.36,2)+pow(10.32,2)));
		histos[mapname]->SetBinContent(6, 43.54); histos[mapname]->SetBinError(6,sqrt(pow( 4.88,2)+pow(13.11,2)));
		histos[mapname]->SetBinContent(7, 23.40); histos[mapname]->SetBinError(7,sqrt(pow( 3.63,2)+pow( 7.07,2)));
		histos[mapname]->SetBinContent(8,  6.43); histos[mapname]->SetBinError(8,sqrt(pow( 1.85,2)+pow( 1.93,2)));
		histos[mapname]->SetBinContent(9,  2.68); histos[mapname]->SetBinError(9,sqrt(pow( 1.20,2)+pow( 0.81,2)));
	  } if(i3==1&&i2==5){//mHT,6j,0b
		histos[mapname]->SetBinContent(1, 10.51); histos[mapname]->SetBinError(1,sqrt(pow( 2.84,2)+pow( 2.43,2)));
		histos[mapname]->SetBinContent(2,  3.03); histos[mapname]->SetBinError(2,sqrt(pow( 1.35,2)+pow( 0.68,2)));
		histos[mapname]->SetBinContent(3,  4.05); histos[mapname]->SetBinError(3,sqrt(pow( 1.65,2)+pow( 0.94,2)));
		histos[mapname]->SetBinContent(4,  2.61); histos[mapname]->SetBinError(4,sqrt(pow( 1.51,2)+pow( 0.88,2)));
		histos[mapname]->SetBinContent(5,  0.66); histos[mapname]->SetBinError(5,sqrt(pow( 0.66,2)+pow( 0.23,2)));
	   } if(i3==2&&i2==0){//hHT,2j,0b
		histos[mapname]->SetBinContent(1, 10.05); histos[mapname]->SetBinError(1,sqrt(pow( 1.99,2)+pow( 2.09,2)));
		histos[mapname]->SetBinContent(2, 10.34); histos[mapname]->SetBinError(2,sqrt(pow( 2.07,2)+pow( 2.14,2)));
		histos[mapname]->SetBinContent(3,  9.22); histos[mapname]->SetBinError(3,sqrt(pow( 2.06,2)+pow( 1.97,2)));
		histos[mapname]->SetBinContent(4,  3.52); histos[mapname]->SetBinError(4,sqrt(pow( 1.25,2)+pow( 0.76,2)));
		histos[mapname]->SetBinContent(5,  2.94); histos[mapname]->SetBinError(5,sqrt(pow( 1.20,2)+pow( 0.91,2)));
		histos[mapname]->SetBinContent(6,  2.45); histos[mapname]->SetBinError(6,sqrt(pow( 1.10,2)+pow( 0.76,2)));
	   } if(i3==2&&i2==2){//hHT,35j,0b
		histos[mapname]->SetBinContent(1,  8.73); histos[mapname]->SetBinError(1,sqrt(pow( 2.13,2)+pow( 1.86,2)));
		histos[mapname]->SetBinContent(2, 11.85); histos[mapname]->SetBinError(2,sqrt(pow( 2.42,2)+pow( 2.48,2)));
		histos[mapname]->SetBinContent(3, 10.16); histos[mapname]->SetBinError(3,sqrt(pow( 2.33,2)+pow( 2.13,2)));
		histos[mapname]->SetBinContent(4,  8.54); histos[mapname]->SetBinError(4,sqrt(pow( 2.20,2)+pow( 1.82,2)));
		histos[mapname]->SetBinContent(5,  4.85); histos[mapname]->SetBinError(5,sqrt(pow( 1.73,2)+pow( 1.52,2)));
		histos[mapname]->SetBinContent(6,  2.88); histos[mapname]->SetBinError(6,sqrt(pow( 1.29,2)+pow( 0.88,2)));
		histos[mapname]->SetBinContent(7,  2.88); histos[mapname]->SetBinError(7,sqrt(pow( 1.29,2)+pow( 0.88,2)));
	   } if(i3==2&&i2==5){//hHT,6j,0b
		histos[mapname]->SetBinContent(1,  2.80); histos[mapname]->SetBinError(1,sqrt(pow( 1.40,2)+pow( 0.80,2)));
		histos[mapname]->SetBinContent(2,  3.43); histos[mapname]->SetBinError(2,sqrt(pow( 1.98,2)+pow( 0.93,2)));
	//	histos[mapname]->SetBinContent(3,); histos[mapname]->SetBinError(3,sqrt(pow(,2)+pow(,2)));
	  } if(i3==0&&i2==1){//lHT,2j,1b
		histos[mapname]->SetBinContent(1, 33.99); histos[mapname]->SetBinError(1,sqrt(pow( 1.21,2)+pow(10.70,2)));
		histos[mapname]->SetBinContent(2, 22.83); histos[mapname]->SetBinError(2,sqrt(pow( 1.00,2)+pow( 7.19,2)));
		histos[mapname]->SetBinContent(3, 17.90); histos[mapname]->SetBinError(3,sqrt(pow( 0.91,2)+pow( 6.97,2)));
		histos[mapname]->SetBinContent(4, 13.92); histos[mapname]->SetBinError(4,sqrt(pow( 0.82,2)+pow( 5.44,2)));
		histos[mapname]->SetBinContent(5,  9.25); histos[mapname]->SetBinError(5,sqrt(pow( 0.65,2)+pow( 3.56,2)));
		histos[mapname]->SetBinContent(6,  1.91); histos[mapname]->SetBinError(6,sqrt(pow( 0.30,2)+pow( 0.74,2)));
	   } if(i3==0&&i2==3){//lHT,35j,1b
		histos[mapname]->SetBinContent(1, 87.24); histos[mapname]->SetBinError(1,sqrt(pow( 2.82,2)+pow(22.79,2)));
		histos[mapname]->SetBinContent(2, 64.63); histos[mapname]->SetBinError(2,sqrt(pow( 2.45,2)+pow(16.77,2)));
		histos[mapname]->SetBinContent(3, 42.61); histos[mapname]->SetBinError(3,sqrt(pow( 1.94,2)+pow(14.56,2)));
		histos[mapname]->SetBinContent(4, 23.77); histos[mapname]->SetBinError(4,sqrt(pow( 1.52,2)+pow( 8.25,2)));
		histos[mapname]->SetBinContent(5, 11.56); histos[mapname]->SetBinError(5,sqrt(pow( 1.04,2)+pow( 3.90,2)));
		histos[mapname]->SetBinContent(6,  2.66); histos[mapname]->SetBinError(6,sqrt(pow( 0.50,2)+pow( 0.91,2)));
	   } if(i3==0&&i2==6){//lHT,6j,1b
		histos[mapname]->SetBinContent(1,  2.12); histos[mapname]->SetBinError(1,sqrt(pow( 0.61,2)+pow( 1.40,2)));
		histos[mapname]->SetBinContent(2,  1.84); histos[mapname]->SetBinError(2,sqrt(pow( 0.53,2)+pow( 1.22,2)));
		histos[mapname]->SetBinContent(3,  1.55); histos[mapname]->SetBinError(3,sqrt(pow( 0.64,2)+pow( 1.11,2)));
	   } if(i3==1&&i2==1){//mHT,2j,1b
		histos[mapname]->SetBinContent(1, 12.04); histos[mapname]->SetBinError(1,sqrt(pow( 0.61,2)+pow( 6.91,2)));
		histos[mapname]->SetBinContent(2,  8.03); histos[mapname]->SetBinError(2,sqrt(pow( 0.54,2)+pow( 4.61,2)));
		histos[mapname]->SetBinContent(3,  8.67); histos[mapname]->SetBinError(3,sqrt(pow( 0.58,2)+pow( 4.97,2)));
		histos[mapname]->SetBinContent(4,  5.38); histos[mapname]->SetBinError(4,sqrt(pow( 0.49,2)+pow( 3.31,2)));
		histos[mapname]->SetBinContent(5,  2.65); histos[mapname]->SetBinError(5,sqrt(pow( 0.34,2)+pow( 1.63,2)));
	   } if(i3==1&&i2==3){//mHT,35j,1b
		histos[mapname]->SetBinContent(1, 16.86); histos[mapname]->SetBinError(1,sqrt(pow( 1.14,2)+pow( 3.89,2)));
		histos[mapname]->SetBinContent(2, 15.66); histos[mapname]->SetBinError(2,sqrt(pow( 1.16,2)+pow( 3.62,2)));
		histos[mapname]->SetBinContent(3, 16.99); histos[mapname]->SetBinError(3,sqrt(pow( 1.18,2)+pow( 3.92,2)));
		histos[mapname]->SetBinContent(4, 13.70); histos[mapname]->SetBinError(4,sqrt(pow( 1.13,2)+pow( 3.17,2)));
		histos[mapname]->SetBinContent(5,  8.91); histos[mapname]->SetBinError(5,sqrt(pow( 0.89,2)+pow( 2.87,2)));
		histos[mapname]->SetBinContent(6,  2.27); histos[mapname]->SetBinError(6,sqrt(pow( 0.45,2)+pow( 0.73,2)));
	   } if(i3==1&&i2==6){//mHT,6j,1b
		histos[mapname]->SetBinContent(1,  2.41); histos[mapname]->SetBinError(1,sqrt(pow( 0.70,2)+pow( 1.39,2)));
		histos[mapname]->SetBinContent(2,  1.55); histos[mapname]->SetBinError(2,sqrt(pow( 0.52,2)+pow( 0.89,2)));
		histos[mapname]->SetBinContent(3,  1.11); histos[mapname]->SetBinError(3,sqrt(pow( 0.45,2)+pow( 0.64,2)));
		histos[mapname]->SetBinContent(4,  0.85); histos[mapname]->SetBinError(4,sqrt(pow( 0.42,2)+pow( 0.53,2)));
	   } if(i3==2&&i2==1){//hHT,2j,1b
		histos[mapname]->SetBinContent(1,  2.69); histos[mapname]->SetBinError(1,sqrt(pow( 0.31,2)+pow( 1.80,2)));
		histos[mapname]->SetBinContent(2,  2.25); histos[mapname]->SetBinError(2,sqrt(pow( 0.31,2)+pow( 1.51,2)));
	   } if(i3==2&&i2==3){//hHT,35j,1b
		histos[mapname]->SetBinContent(1,  2.16); histos[mapname]->SetBinError(1,sqrt(pow( 0.41,2)+pow( 0.61,2)));
		histos[mapname]->SetBinContent(2,  2.40); histos[mapname]->SetBinError(2,sqrt(pow( 0.43,2)+pow( 0.67,2)));
		histos[mapname]->SetBinContent(3,  2.57); histos[mapname]->SetBinError(3,sqrt(pow( 0.49,2)+pow( 0.72,2)));
		histos[mapname]->SetBinContent(4,  1.70); histos[mapname]->SetBinError(4,sqrt(pow( 0.40,2)+pow( 0.61,2)));
	   } if(i3==2&&i2==6){//hHT,6j,1b
		histos[mapname]->SetBinContent(1,  1.06); histos[mapname]->SetBinError(1,sqrt(pow( 0.44,2)+pow( 0.84,2)));
		histos[mapname]->SetBinContent(2,  0.92); histos[mapname]->SetBinError(2,sqrt(pow( 0.53,2)+pow( 0.73,2)));
	//	histos[mapname]->SetBinContent(3,); histos[mapname]->SetBinError(3,sqrt(pow(,2)+pow(,2)));
	   }
	}}

	cout << "Saving." << endl;
    	TFile *fsavefile = new TFile("../Results/Filtered/GammaJetsPrediction/20130617_test/ZinvPredictionNumbers.root","RECREATE");
	fsavefile->cd();
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Write();
	}
	fsavefile->Close();
	cout << "Saved histograms in " << fsavefile->GetName() << endl;

	//make the plots - not TDR style
    TLatex TitleBox;
	TitleBox.SetNDC();
    TitleBox.SetNDC();
    TitleBox.SetTextAlign(12);
    TitleBox.SetTextFont(42);
    TitleBox.SetTextSize(0.04219409);
    TitleBox.SetLineWidth(2);
	TString text;
    string outname;
    double max = 0.;
    double max1,max2;
    string outputdir = "../Results/Filtered/GammaJetsPrediction/20130617_test/ZinvPredictionPlots/";
    	Util::MakeOutputDir(outputdir);

   TCanvas *c1 = new TCanvas("c1", "c1",485,220,700,504);
   c1->Range(82.71719,-0.4425771,532.9945,2.212885);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
   c1->SetLeftMargin(0.1494253);
   c1->SetRightMargin(0.07327586);
   c1->SetTopMargin(0.08016878);
   c1->SetBottomMargin(0.1666667);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);
	for(int i3 = 0; i3<HTbinsize;        ++i3){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	   string hs = string("_") + HT_bin[i3] + string("_") + signal_region[i2];// + string("_") + sample_type[i1];
	   string mapname   = "MT2pred"  + hs;
	   string mapnameMC = "MT2truth" + hs;
	   if(i3==0&&i2==0){//lHT,2j,0b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "2 jets, 0 b-jets, low H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "2 jets, 0 b-jets, low H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==0&&i2==2){//lHT,35j,0b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "3-5 jets, 0 b-jets, low H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "3-5 jets, 0 b-jets, low H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==0&&i2==5){//lHT,6j,0b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "#geq6 jets, 0 b-jets, low H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "#geq6 jets, 0 b-jets, low H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==1&&i2==0){//mHT,2j,0b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "2 jets, 0 b-jets, medium H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "2 jets, 0 b-jets, medium H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==1&&i2==2){//mHT,35j,0b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "3-5 jets, 0 b-jets, medium H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "3-5 jets, 0 b-jets, medium H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==1&&i2==5){//mHT,6j,0b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "#geq6 jets, 0 b-jets, medium H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.2);
		histos[mapnameMC]->SetMinimum(0.2);
		text = "#geq6 jets, 0 b-jets, medium H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==2&&i2==0){//hHT,2j,0b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "2 jets, 0 b-jets, high H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "2 jets, 0 b-jets, high H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==2&&i2==2){//hHT,35j,0b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "3-5 jets, 0 b-jets, high H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "3-5 jets, 0 b-jets, high H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==2&&i2==5){//hHT,6j,0b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
	//	max = 1.5*max;
		max = 5.5;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "#geq6 jets, 0 b-jets, high H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "#geq6 jets, 0 b-jets, high H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==0&&i2==1){//lHT,2j,1b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "2 jets, #geq1 b-jets, low H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "2 jets, #geq1 b-jets, low H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==0&&i2==3){//lHT,35j,1b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "3-5 jets, 1 b-jet, low H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "3-5 jets, 1 b-jet, low H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==0&&i2==6){//lHT,6j,1b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
	//	max = 1.5*max;
		max = 3.5;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "#geq6 jets, 1 b-jet, low H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.3);
		histos[mapnameMC]->SetMinimum(0.3);
		text = "#geq6 jets, 1 b-jet, low H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==1&&i2==1){//mHT,2j,1b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "2 jets, #geq1 b-jets, medium H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "2 jets, #geq1 b-jets, medium H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==1&&i2==3){//mHT,35j,1b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "3-5 jets, 1 b-jet, medium H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "3-5 jets, 1 b-jet, medium H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==1&&i2==6){//mHT,6j,1b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "#geq6 jets, 1 b-jet, medium H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.3);
		histos[mapnameMC]->SetMinimum(0.3);
		text = "#geq6 jets, 1 b-jet, medium H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==2&&i2==1){//hHT,2j,1b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
	//	max = 1.5*max;
		max = 6.5;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "2 jets, #geq1 b-jets, high H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "2 jets, #geq1 b-jets, high H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==2&&i2==3){//hHT,35j,1b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 1.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "3-5 jets, 1 b-jet, high H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(1.);
		histos[mapnameMC]->SetMinimum(1.);
		text = "3-5 jets, 1 b-jet, high H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
		c1->Clear();
	   } if(i3==2&&i2==6){//hHT,6j,1b
		c1->Clear();
		c1->cd();
    		gPad->SetLogy(0);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
	//	max = 1.5*max;
		max = 2.;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.);
		histos[mapnameMC]->SetMinimum(0.);
		text = "#geq6 jets, 1 b-jet, high H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + ".eps";
		c1->SaveAs(outname.c_str());

		c1->Clear();
		c1->cd();
    		gPad->SetLogy(1);
		max1=histos[mapname]  ->GetMaximum();
		max2=histos[mapnameMC]->GetMaximum();
		max  = (max1>max2)?max1:max2;
		max = 2.5*max;
		histos[mapname  ]->SetMaximum(max);
		histos[mapnameMC]->SetMaximum(max);
		histos[mapname  ]->SetMinimum(0.2);
		histos[mapnameMC]->SetMinimum(0.2);
		text = "#geq6 jets, 1 b-jet, high H_{T}";
		histos[mapnameMC]->Draw("E2");
		histos[mapname  ]->Draw("sameE1");
		leg->Draw();
		TitleBox.DrawLatex(0.1494253,0.9493671,text.Data());
		outname = histos[mapname  ]->GetName();outname = outputdir + outname + "_log.eps";
		c1->SaveAs(outname.c_str());
	//	c1->Clear();
	   }
	}}

}