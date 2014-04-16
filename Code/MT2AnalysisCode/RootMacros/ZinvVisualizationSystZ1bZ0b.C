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

//run via root -l -b -q ZinvVisualizationZ1b0b.C++

//this code is the same as ZinvVisualization.C, but the uncertainty on the data prediction is only due to the Zll(1b)/Zll(0b) scaling
//the produced root file is used for the interpretation
//this code does not produce any plot
void ZinvVisualizationZ1b0b(){

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

TFile *fmctruth = TFile::Open("../Results/GammaJetsPrediction/20130617_test/ZnunuNumbers.root");

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
	//	cout << mapname << endl;
		TH1D *h = (TH1D*)fmctruth->Get(mapname.c_str());
	//	cout << h->Integral() << endl;
		mapname = "MT2truth" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)h->Clone(mapname.c_str());
		histos[mapname]->SetFillColor(kViolet-3); histos[mapname]->SetFillStyle(3001);
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

	for(int i3 = 0; i3<HTbinsize;        ++i3){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	   string hs = string("_") + HT_bin[i3] + string("_") + signal_region[i2];// + string("_") + sample_type[i1];
	   string mapname = "MT2pred" + hs;
	   if(i3==0&&i2==0){//lHT,2j,0b
		histos[mapname]->SetBinContent(1,304.27); histos[mapname]->SetBinError(1,0.*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2,236.52); histos[mapname]->SetBinError(2,0.*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3,182.98); histos[mapname]->SetBinError(3,0.*histos[mapname]->GetBinContent(3));
		histos[mapname]->SetBinContent(4,160.92); histos[mapname]->SetBinError(4,0.*histos[mapname]->GetBinContent(4));
		histos[mapname]->SetBinContent(5,116.69); histos[mapname]->SetBinError(5,0.*histos[mapname]->GetBinContent(5));
		histos[mapname]->SetBinContent(6, 47.95); histos[mapname]->SetBinError(6,0.*histos[mapname]->GetBinContent(6));
		histos[mapname]->SetBinContent(7, 11.74); histos[mapname]->SetBinError(7,0.*histos[mapname]->GetBinContent(7));
		histos[mapname]->SetBinContent(8,  3.28); histos[mapname]->SetBinError(8,0.*histos[mapname]->GetBinContent(8));
	   } if(i3==0&&i2==2){//lHT,35j,0b
		histos[mapname]->SetBinContent(1,458.79); histos[mapname]->SetBinError(1,0.*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2,366.96); histos[mapname]->SetBinError(2,0.*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3,301.04); histos[mapname]->SetBinError(3,0.*histos[mapname]->GetBinContent(3));
		histos[mapname]->SetBinContent(4,173.25); histos[mapname]->SetBinError(4,0.*histos[mapname]->GetBinContent(4));
		histos[mapname]->SetBinContent(5, 93.60); histos[mapname]->SetBinError(5,0.*histos[mapname]->GetBinContent(5));
		histos[mapname]->SetBinContent(6, 41.36); histos[mapname]->SetBinError(6,0.*histos[mapname]->GetBinContent(6));
		histos[mapname]->SetBinContent(7,  9.80); histos[mapname]->SetBinError(7,0.*histos[mapname]->GetBinContent(7));
		histos[mapname]->SetBinContent(8,  2.73); histos[mapname]->SetBinError(8,0.*histos[mapname]->GetBinContent(8));
	   } if(i3==0&&i2==5){//lHT,6j,0b
		histos[mapname]->SetBinContent(1, 12.02); histos[mapname]->SetBinError(1,0.*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2,  5.13); histos[mapname]->SetBinError(2,0.*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3,  3.53); histos[mapname]->SetBinError(3,0.*histos[mapname]->GetBinContent(3));
	   } if(i3==1&&i2==0){//mHT,2j,0b
		histos[mapname]->SetBinContent(1, 76.13); histos[mapname]->SetBinError(1,0.*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2, 62.24); histos[mapname]->SetBinError(2,0.*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3, 38.49); histos[mapname]->SetBinError(3,0.*histos[mapname]->GetBinContent(3));
		histos[mapname]->SetBinContent(4, 40.84); histos[mapname]->SetBinError(4,0.*histos[mapname]->GetBinContent(4));
		histos[mapname]->SetBinContent(5, 20.11); histos[mapname]->SetBinError(5,0.*histos[mapname]->GetBinContent(5));
		histos[mapname]->SetBinContent(6, 30.21); histos[mapname]->SetBinError(6,0.*histos[mapname]->GetBinContent(6));
		histos[mapname]->SetBinContent(7, 12.73); histos[mapname]->SetBinError(7,0.*histos[mapname]->GetBinContent(7));
		histos[mapname]->SetBinContent(8, 15.96); histos[mapname]->SetBinError(8,0.*histos[mapname]->GetBinContent(8));
		histos[mapname]->SetBinContent(9,  2.83); histos[mapname]->SetBinError(9,0.*histos[mapname]->GetBinContent(9));
	   } if(i3==1&&i2==2){//mHT,35j,0b
		histos[mapname]->SetBinContent(1, 88.17); histos[mapname]->SetBinError(1,0.*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2, 73.11); histos[mapname]->SetBinError(2,0.*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3, 63.34); histos[mapname]->SetBinError(3,0.*histos[mapname]->GetBinContent(3));
		histos[mapname]->SetBinContent(4, 59.03); histos[mapname]->SetBinError(4,0.*histos[mapname]->GetBinContent(4));
		histos[mapname]->SetBinContent(5, 50.85); histos[mapname]->SetBinError(5,0.*histos[mapname]->GetBinContent(5));
		histos[mapname]->SetBinContent(6, 43.54); histos[mapname]->SetBinError(6,0.*histos[mapname]->GetBinContent(6));
		histos[mapname]->SetBinContent(7, 23.40); histos[mapname]->SetBinError(7,0.*histos[mapname]->GetBinContent(7));
		histos[mapname]->SetBinContent(8,  6.43); histos[mapname]->SetBinError(8,0.*histos[mapname]->GetBinContent(8));
		histos[mapname]->SetBinContent(9,  2.68); histos[mapname]->SetBinError(9,0.*histos[mapname]->GetBinContent(9));
	  } if(i3==1&&i2==5){//mHT,6j,0b
		histos[mapname]->SetBinContent(1, 10.51); histos[mapname]->SetBinError(1,0.*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2,  3.03); histos[mapname]->SetBinError(2,0.*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3,  4.05); histos[mapname]->SetBinError(3,0.*histos[mapname]->GetBinContent(3));
		histos[mapname]->SetBinContent(4,  2.61); histos[mapname]->SetBinError(4,0.*histos[mapname]->GetBinContent(4));
		histos[mapname]->SetBinContent(5,  0.66); histos[mapname]->SetBinError(5,0.*histos[mapname]->GetBinContent(5));
	   } if(i3==2&&i2==0){//hHT,2j,0b
		histos[mapname]->SetBinContent(1, 10.05); histos[mapname]->SetBinError(1,0.*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2, 10.34); histos[mapname]->SetBinError(2,0.*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3,  9.22); histos[mapname]->SetBinError(3,0.*histos[mapname]->GetBinContent(3));
		histos[mapname]->SetBinContent(4,  3.52); histos[mapname]->SetBinError(4,0.*histos[mapname]->GetBinContent(4));
		histos[mapname]->SetBinContent(5,  2.94); histos[mapname]->SetBinError(5,0.*histos[mapname]->GetBinContent(5));
		histos[mapname]->SetBinContent(6,  2.45); histos[mapname]->SetBinError(6,0.*histos[mapname]->GetBinContent(6));
	   } if(i3==2&&i2==2){//hHT,35j,0b
		histos[mapname]->SetBinContent(1,  8.73); histos[mapname]->SetBinError(1,0.*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2, 11.85); histos[mapname]->SetBinError(2,0.*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3, 10.16); histos[mapname]->SetBinError(3,0.*histos[mapname]->GetBinContent(3));
		histos[mapname]->SetBinContent(4,  8.54); histos[mapname]->SetBinError(4,0.*histos[mapname]->GetBinContent(4));
		histos[mapname]->SetBinContent(5,  4.85); histos[mapname]->SetBinError(5,0.*histos[mapname]->GetBinContent(5));
		histos[mapname]->SetBinContent(6,  2.88); histos[mapname]->SetBinError(6,0.*histos[mapname]->GetBinContent(6));
		histos[mapname]->SetBinContent(7,  2.88); histos[mapname]->SetBinError(7,0.*histos[mapname]->GetBinContent(7));
	   } if(i3==2&&i2==5){//hHT,6j,0b
		histos[mapname]->SetBinContent(1,  2.80); histos[mapname]->SetBinError(1,0.*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2,  3.43); histos[mapname]->SetBinError(2,0.*histos[mapname]->GetBinContent(2));
	//	histos[mapname]->SetBinContent(3,); histos[mapname]->SetBinError(3,sqrt(pow(,2)+pow(,2)));
	  } if(i3==0&&i2==1){//lHT,2j,1b
		histos[mapname]->SetBinContent(1, 33.99); histos[mapname]->SetBinError(1,(0.022/0.093)*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2, 22.83); histos[mapname]->SetBinError(2,(0.022/0.093)*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3, 17.90); histos[mapname]->SetBinError(3,(0.022/0.093)*histos[mapname]->GetBinContent(3));
		histos[mapname]->SetBinContent(4, 13.92); histos[mapname]->SetBinError(4,(0.022/0.093)*histos[mapname]->GetBinContent(4));
		histos[mapname]->SetBinContent(5,  9.25); histos[mapname]->SetBinError(5,(0.022/0.093)*histos[mapname]->GetBinContent(5));
		histos[mapname]->SetBinContent(6,  1.91); histos[mapname]->SetBinError(6,(0.022/0.093)*histos[mapname]->GetBinContent(6));
	   } if(i3==0&&i2==3){//lHT,35j,1b
		histos[mapname]->SetBinContent(1, 87.24); histos[mapname]->SetBinError(1,(0.024/0.162)*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2, 64.63); histos[mapname]->SetBinError(2,(0.024/0.162)*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3, 42.61); histos[mapname]->SetBinError(3,(0.024/0.162)*histos[mapname]->GetBinContent(3));
		histos[mapname]->SetBinContent(4, 23.77); histos[mapname]->SetBinError(4,(0.024/0.162)*histos[mapname]->GetBinContent(4));
		histos[mapname]->SetBinContent(5, 11.56); histos[mapname]->SetBinError(5,(0.024/0.162)*histos[mapname]->GetBinContent(5));
		histos[mapname]->SetBinContent(6,  2.66); histos[mapname]->SetBinError(6,(0.024/0.162)*histos[mapname]->GetBinContent(6));
	   } if(i3==0&&i2==6){//lHT,6j,1b
		histos[mapname]->SetBinContent(1,  2.12); histos[mapname]->SetBinError(1,(0.168/0.269)*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2,  1.84); histos[mapname]->SetBinError(2,(0.168/0.269)*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3,  1.55); histos[mapname]->SetBinError(3,(0.168/0.269)*histos[mapname]->GetBinContent(3));
	   } if(i3==1&&i2==1){//mHT,2j,1b
		histos[mapname]->SetBinContent(1, 12.04); histos[mapname]->SetBinError(1,(0.05/0.093)*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2,  8.03); histos[mapname]->SetBinError(2,(0.05/0.093)*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3,  8.67); histos[mapname]->SetBinError(3,(0.05/0.093)*histos[mapname]->GetBinContent(3));
		histos[mapname]->SetBinContent(4,  5.38); histos[mapname]->SetBinError(4,(0.05/0.093)*histos[mapname]->GetBinContent(4));
		histos[mapname]->SetBinContent(5,  2.65); histos[mapname]->SetBinError(5,(0.05/0.093)*histos[mapname]->GetBinContent(5));
	   } if(i3==1&&i2==3){//mHT,35j,1b
		histos[mapname]->SetBinContent(1, 16.86); histos[mapname]->SetBinError(1,(0.018/0.162)*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2, 15.66); histos[mapname]->SetBinError(2,(0.018/0.162)*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3, 16.99); histos[mapname]->SetBinError(3,(0.018/0.162)*histos[mapname]->GetBinContent(3));
		histos[mapname]->SetBinContent(4, 13.70); histos[mapname]->SetBinError(4,(0.018/0.162)*histos[mapname]->GetBinContent(4));
		histos[mapname]->SetBinContent(5,  8.91); histos[mapname]->SetBinError(5,(0.018/0.162)*histos[mapname]->GetBinContent(5));
		histos[mapname]->SetBinContent(6,  2.27); histos[mapname]->SetBinError(6,(0.018/0.162)*histos[mapname]->GetBinContent(6));
	   } if(i3==1&&i2==6){//mHT,6j,1b
		histos[mapname]->SetBinContent(1,  2.41); histos[mapname]->SetBinError(1,(0.143/0.269)*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2,  1.55); histos[mapname]->SetBinError(2,(0.143/0.269)*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3,  1.11); histos[mapname]->SetBinError(3,(0.143/0.269)*histos[mapname]->GetBinContent(3));
		histos[mapname]->SetBinContent(4,  0.85); histos[mapname]->SetBinError(4,(0.143/0.269)*histos[mapname]->GetBinContent(4));
	   } if(i3==2&&i2==1){//hHT,2j,1b
		histos[mapname]->SetBinContent(1,  2.69); histos[mapname]->SetBinError(1,(0.06/0.093)*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2,  2.25); histos[mapname]->SetBinError(2,(0.06/0.093)*histos[mapname]->GetBinContent(2));
	   } if(i3==2&&i2==3){//hHT,35j,1b
		histos[mapname]->SetBinContent(1,  2.16); histos[mapname]->SetBinError(1,(0.03/0.162)*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2,  2.40); histos[mapname]->SetBinError(2,(0.03/0.162)*histos[mapname]->GetBinContent(2));
		histos[mapname]->SetBinContent(3,  2.57); histos[mapname]->SetBinError(3,(0.03/0.162)*histos[mapname]->GetBinContent(3));
		histos[mapname]->SetBinContent(4,  1.70); histos[mapname]->SetBinError(4,(0.03/0.162)*histos[mapname]->GetBinContent(4));
	   } if(i3==2&&i2==6){//hHT,6j,1b
		histos[mapname]->SetBinContent(1,  1.06); histos[mapname]->SetBinError(1,(0.201/0.269)*histos[mapname]->GetBinContent(1));
		histos[mapname]->SetBinContent(2,  0.92); histos[mapname]->SetBinError(2,(0.201/0.269)*histos[mapname]->GetBinContent(2));
	//	histos[mapname]->SetBinContent(3,); histos[mapname]->SetBinError(3,sqrt(pow(,2)+pow(,2)));
	   }
	}}

	cout << "Saving." << endl;
	//if(!fISRreweight) outputname = "NoISR_" + outputname;
    	TFile *fsavefile = new TFile("../Results/Filtered/GammaJetsPrediction/20130617_test/ZinvPredictionNumbersSystOnlyZllErr.root","RECREATE");
	fsavefile->cd();
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Write();
	}
	fsavefile->Close();
	cout << "Saved histograms in " << fsavefile->GetName() << endl;

}