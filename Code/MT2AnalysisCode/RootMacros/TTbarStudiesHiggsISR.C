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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use own path of MT2tree.hh
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/Utilities.hh"

//run via root -l -b -q TTbarStudiesHiggs.C++

//this function is basically a copy of TTbarStudies.C
//but shapes are done for Mbb instead of MT2, event selection is the MT2 Higgs selection instead of inclusive MT2 analysis
//therefore no additional comments are given

using namespace std;

void load(const char* filename = "samples/datasamples/samples_2141_dataonly.dat");
void TTbarStudiesHiggs();

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

Bool_t fbTagReweight    = true;//for all shapes (except BTV SF systematics) apply BTV SF weights
Bool_t fISRreweight     = true;//for nominal shape (and successive other shapes except ISR systematics) apply the 'ISR reweighting', default = true, but consider also false
Bool_t WnoScaleMatching = true;//do not compute scale/matching up/down (default == true due to statistics of those W samples)

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

double ISRweights[4] = {1.00,0.95,0.90,0.80};
double ISupperpts[4] = {120.,150.,250.,99999.};

const int sampletypesize = 32;
string sample_type[sampletypesize] = {"Nominal", "MatchingUp", "MatchingDown", "ScaleUp", "ScaleDown", "MassUp", "MassDown", "JESUp", "JESDown", "METUp", "METDown", "BSFUp", "BSFDown", "ISRUp", "ISR", "ISRDown", "PUUp", "PUDown", "TopPtReweighted", "nominal0l", "nominal1l", "nominal2l", "nominal1Sample", "nominalPowheg", "WUp", "WDown", "TopUp", "TopDown", "SingleTopUp", "SingleTopDown", "TTVUp", "TTVDown"};//maybe not all filled  - ! include TTbar and WJets
const int samplekindsize = 6;
string sample_kind[samplekindsize] = {"allMC","allTop", "TTbar", "SingleTop", "TTV", "WJets"};
const int HTregionsize = 3;
string HT_region[HTregionsize] = {"lowHT", "mediumHT", "highHT"};
const int signalregionsize = 9;
string signal_region[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b"};

void TTbarStudiesHiggsISR(){


	TFile *pu = TFile::Open("/shome/haweber/MT2Analysis_8TeV/Code/Certification/pileup_weights.root");
	TH1D *pu_up = (TH1D*)pu->Get("pileup_weight_xsecPlusFivePercent");
	TH1D *pu_down = (TH1D*)pu->Get("pileup_weight_xsecMinusFivePercent");

//	gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
//	gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
//	gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
//	gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

//	JetCorrectionUncertainty *totalunc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("/shome/haweber/MT2Analysis_8TeV/Code/JetEnergyCorrection/FALL12_V7_Uncertainties/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt", "Total")));
	

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
//	if(fMET==true) fHT = false;
//	if(fMET==false && fHT==false) fHT = true;
	
	TString  outputdir = "Filtered/TTbarStudies/Higgs/";
    	Util::MakeOutputDir(outputdir);

	TString outputname = "TTbarStudiesHistograms_all_noWscaleupdown.root";
	if(!WnoScaleMatching) outputname = "TTbarStudiesHistograms_all.root";
	if(fISRreweight){
		outputname = "TTbarStudiesHistograms_all_noWscaleupdown_ISRdefault.root";
		if(!WnoScaleMatching) outputname = "TTbarStudiesHistograms_all_ISRdefault.root";
	}

	TString  samples = "samples/samples_TTbar_filter.dat";//only dummy
	map<string, TH1D*>    histos;
	map<string, TH1D*>    histosnormalized;
	for(int i0 = 0; i0<samplekindsize;   ++i0){
	for(int i1 = 0; i1<sampletypesize;   ++i1){
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
  		double MT2bins[NMT2bins];
		if(i3==0){
			if(signal_region[i2]=="2j0b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_2j0b_lHT[in]; }
			if(signal_region[i2]=="2j1to2b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_2j1b_lHT[in]; }
			if(signal_region[i2]=="3to5j0b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3j0b_lHT[in]; }
			if(signal_region[i2]=="3to5j1b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3j1b_lHT[in]; }
			if(signal_region[i2]=="3to5j2b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3j2b_lHT[in]; }
			if(signal_region[i2]=="6j0b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_6j0b_lHT[in]; }
			if(signal_region[i2]=="6j1b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_6j1b_lHT[in]; }
			if(signal_region[i2]=="6j2b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_6j2b_lHT[in]; }
			if(signal_region[i2]=="3b")      { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3b_lHT[in];   }
		} if(i3==1){
			if(signal_region[i2]=="2j0b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_2j0b_mHT[in];  }
			if(signal_region[i2]=="2j1to2b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_2j1b_mHT[in];  }
			if(signal_region[i2]=="3to5j0b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3j0b_mHT[in];  }
			if(signal_region[i2]=="3to5j1b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3j1b_mHT[in];  }
			if(signal_region[i2]=="3to5j2b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3j2b_mHT[in];  }
			if(signal_region[i2]=="6j0b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_6j0b_mHT[in];  }
			if(signal_region[i2]=="6j1b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_6j1b_mHT[in];  }
			if(signal_region[i2]=="6j2b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_6j2b_mHT[in];  }
			if(signal_region[i2]=="3b")      { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3b_mHT[in];    }
		} if(i3==2){
			if(signal_region[i2]=="2j0b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_2j0b_hHT[in];  }
			if(signal_region[i2]=="2j1to2b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_2j1b_hHT[in];  }
			if(signal_region[i2]=="3to5j0b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3j0b_hHT[in];  }
			if(signal_region[i2]=="3to5j1b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3j1b_hHT[in];  }
			if(signal_region[i2]=="3to5j2b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3j2b_hHT[in];  }
			if(signal_region[i2]=="6j0b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_6j0b_hHT[in];  }
			if(signal_region[i2]=="6j1b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_6j1b_hHT[in];  }
			if(signal_region[i2]=="6j2b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_6j2b_hHT[in];  }
			if(signal_region[i2]=="3b")      { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3b_hHT[in];    }
		}
		//for(int i = 0; i < NMT2bins; ++i) cout << signal_region[i2] << ": " << i+1 << "/" << NMT2bins << ": " << MT2bins[i] << endl;

		string hs = string("_") + sample_type[i1] + string("_") + sample_kind[i0] + string("_") + HT_region[i3] + string("_") + signal_region[i2];
		string mapname = "MT2" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", NMT2bins, MT2bins);
		Color_t colour = 1;
	//	if(i1==0) colour = 600; if(i1==1) colour = 417; if(i1==2) colour = 419; if(i1==3 ) colour = 401; if(i1==4 ) colour = 403; if(i1==5 ) colour = 632; if(i1==6 ) colour = 1;
	//	if(i1==7) colour = 419; if(i1==8) colour = 417; if(i1==9) colour = 403; if(i1==10) colour = 401; if(i1==11) colour = 620; if(i1==12) colour = 616; if(i1==13) colour = 632; if(i1==14) colour = 800;
		histos[mapname]->SetLineColor(colour); histos[mapname]->SetMarkerColor(colour); histos[mapname]->SetLineWidth(2);
		histos[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); histos[mapname]->GetYaxis()->SetTitle("Events"); 

	}}}}

	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		//cout << h->first << "    " << h->second->GetName() << endl;
		h->second->Sumw2();}

	std::ostringstream cutStream;
	std::ostringstream cutStreamBase;
	cutStream << " " 
		<< "NJetsIDLoose>=2"                        << "&&"//preselection only
		<< "NTausIDLoose3Hits==0"                   << "&&"
		<< "NEles==0"                               << "&&"
		<< "NMuons==0"                              << "&&"
		<< "misc.HT>400" << "&&"//;//preselection only
		<< "misc.Jet0Pass ==1"                      << "&&"
		<< "misc.Jet1Pass ==1"                      << "&&"
		<< "misc.Vectorsumpt < 95";//precut
	cutStream  << "&& misc.MinMetJetDPhi4Pt40 >0.3";//have nowhere angular uncertainties
	//if(fFast)cutStream << "&&misc.MT2>=100";//lowest border in MT2
	cutStream << "&&((misc.HT>450&&misc.HT<750&&misc.MET>160)||(misc.HT>750&&misc.MET>10))";
	
	cutStreamBase << " " 
      << "misc.PassJet40ID ==1"                      << "&&"
      << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"  << "&&" // for rare SM samples
      << "misc.CSCTightHaloIDFlag == 0"             << "&&"
      << "misc.trackingFailureFlag==0"              << "&&"
      << "misc.eeBadScFlag==0"                      << "&&"
      << "misc.EcalDeadCellTriggerPrimitiveFlag==0" << "&&"
      << "misc.TrackingManyStripClusFlag==0"             << "&&"
      << "misc.TrackingTooManyStripClusFlag==0"          << "&&"
      << "misc.TrackingLogErrorTooManyClustersFlag==0"   << "&&"
      << "misc.CrazyHCAL==0";
/*
	std::ostringstream triggerStream;
	if(fMET){
	triggerStream << "( ( "
			<< "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
			<< "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )"
			<< "||("
			<< "trigger.HLT_PFHT350_PFMET100_v3==1 || trigger.HLT_PFHT350_PFMET100_v4==1 || trigger.HLT_PFHT350_PFMET100_v5==1 || "
			<< "trigger.HLT_PFHT350_PFMET100_v6==1 || trigger.HLT_PFHT350_PFMET100_v7==1 || trigger.HLT_PFNoPUHT350_PFMET100_v1==1 || "
			<< "trigger.HLT_PFNoPUHT350_PFMET100_v3==1 || trigger.HLT_PFNoPUHT350_PFMET100_v4==1 ) )";
	}
	if(fHT){
	triggerStream << "( "
			<< "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || "
			<< "trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
			<< "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1)";
	}
	TString trigger = triggerStream.str().c_str();
*/
	TString cuts = cutStream.str().c_str();
	TString basecuts = cutStreamBase.str().c_str();

	load(samples.Data());

   	for(size_t i = 0; i < fSamples.size(); ++i){
        


		TString kindsample = fSamples[i].sname;
		kindsample.ReplaceAll("_lowHT","");
		kindsample.ReplaceAll("_highHT","");

		TString namesample = fSamples[i].name;
		namesample.ReplaceAll("_lowHT","");
		namesample.ReplaceAll("_highHT","");
		namesample.ReplaceAll("v1","");
		namesample.ReplaceAll("v2","");

	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "SampleName: looping over " << fSamples[i].name << " added in " << kindsample << endl;
	    if(fVerbose>2) cout << "            sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

		string nominal = " ";
		if(namesample=="TTbar_Madgraph0l" || namesample=="TTbar_Madgraph1l" || namesample=="TTbar_Madgraph2l")                          nominal = "TTbar";
		if(namesample=="WJetsHT200"||namesample=="WJetsHT250"||namesample=="WJetsHT300"||namesample=="WJetsHT400")                      nominal = "WJets";
		if(namesample=="T_s"||namesample=="Tbar_s"||namesample=="T_t"||namesample=="Tbar_t"||namesample=="T_tW"||namesample=="Tbar_tW") nominal = "SingleTop";
		if(namesample=="WW"||namesample=="WZ"||namesample=="ZZ"||namesample=="TTW"||namesample=="TTZ"||namesample=="TTG")               nominal = "TTV";
		cout << "this sample is " << nominal << " of kind " << kindsample << " with name " << namesample << endl;

	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts + "&&" + basecuts;
        
	//    if( fSamples[i].type=="data") myCuts += " && " + trigger; //cuts to be aplied only on data
        
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
		
		if ( fVerbose>2 && counter % 50000 == 0 )  cout << "+++ Proccessing event " << counter << endl;
		
		Double_t weight = sample_weight;
		if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 

		//these are the values that can change
		int NumBJets  = fMT2tree->NBJets40CSVM;
		int NumJets   = fMT2tree->NJetsIDLoose40;
		float HTval   = fMT2tree->misc.HT;
		float METval  = fMT2tree->misc.MET;
		float MinDPhi = fMT2tree->misc.MinMetJetDPhi4Pt40;
		float VSPT    = fMT2tree->misc.Vectorsumpt;
		float J1Pt    = fMT2tree->misc.LeadingJPt;
		float J2Pt    = fMT2tree->misc.SecondJPt;
		float MT2val  = fMT2tree->misc.MT2;

//		if ( fVerbose>2 && counter % 500 == 0 )  cout << NumJets << "j,"<<NumBJets<<"b, HT "<<HTval<<", MET " << METval << ", MT2 " << MT2val << ", MinDPhi " << MinDPhi << ", VSPT " << VSPT << ", J1/2Pt " << J1Pt << "/" << J2Pt << endl;


		int signalbin(-1); string ssignal = "";//note: this can change, depending on redefinitions
		if(NumJets == 2 &&                 NumBJets == 0) { signalbin = 0; ssignal = "_2j0b";    }
		if(NumJets == 2 &&                 NumBJets >= 1) { signalbin = 1; ssignal = "_2j1to2b"; }
		if(NumJets >= 3 && NumJets <= 5 && NumBJets == 0) { signalbin = 2; ssignal = "_3to5j0b"; }
		if(NumJets >= 3 && NumJets <= 5 && NumBJets == 1) { signalbin = 3; ssignal = "_3to5j1b"; }
		if(NumJets >= 3 && NumJets <= 5 && NumBJets == 2) { signalbin = 4; ssignal = "_3to5j2b"; }
		if(                                NumBJets >= 3) { signalbin = 5; ssignal = "_3b";      }
		if(NumJets >= 6 &&                 NumBJets == 0) { signalbin = 6; ssignal = "_6j0b";    }
		if(NumJets >= 6 &&                 NumBJets == 1) { signalbin = 7; ssignal = "_6j1b";    }
		if(NumJets >= 6 &&                 NumBJets == 2) { signalbin = 8; ssignal = "_6j2b";    }
		double btagSF(1.), btagSFerr(0.);
		if(!fMT2tree->misc.isData && fbTagReweight){
		if(NumJets == 2 &&                 NumBJets == 0) { btagSF = fMT2tree->SFWeight.BTagCSV40eq0; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq0Error; }
		if(NumJets == 2 &&                 NumBJets >= 1) { btagSF = fMT2tree->SFWeight.BTagCSV40ge1; btagSFerr = fMT2tree->SFWeight.BTagCSV40ge1Error; }
		if(NumJets >= 3 && NumJets <= 5 && NumBJets == 0) { btagSF = fMT2tree->SFWeight.BTagCSV40eq0; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq0Error; }
		if(NumJets >= 3 && NumJets <= 5 && NumBJets == 1) { btagSF = fMT2tree->SFWeight.BTagCSV40eq1; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq1Error; }
		if(NumJets >= 3 && NumJets <= 5 && NumBJets == 2) { btagSF = fMT2tree->SFWeight.BTagCSV40eq2; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq2Error; }
		if(NumJets >= 3 &&                 NumBJets >= 3) { btagSF = fMT2tree->SFWeight.BTagCSV40ge3; btagSFerr = fMT2tree->SFWeight.BTagCSV40ge3Error; }
		if(NumJets >= 6 &&                 NumBJets == 0) { btagSF = fMT2tree->SFWeight.BTagCSV40eq0; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq0Error; }
		if(NumJets >= 6 &&                 NumBJets == 1) { btagSF = fMT2tree->SFWeight.BTagCSV40eq1; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq1Error; }
		if(NumJets >= 6 &&                 NumBJets == 2) { btagSF = fMT2tree->SFWeight.BTagCSV40eq2; btagSFerr = fMT2tree->SFWeight.BTagCSV40eq2Error; }
			weight = weight * btagSF;
			//if(!fbTagError) btagSFerr = 0;//don't store btag error is equivalent to no error at all
		}
		int HTbin(-1); string sHT = "";
		if(HTval>=450.  && HTval<750.  && METval>200.) { HTbin = 0; sHT = "_lowHT";    }
		if(HTval>=750.  && HTval<1200. && METval>30. ) { HTbin = 1; sHT = "_mediumHT"; }
		if(HTval>=1200.                && METval>30. ) { HTbin = 2; sHT = "_highHT";   }

		string hh = string("_") + namesample.Data() + string("_") + kindsample.Data() + sHT + ssignal;
//bool printtt = false;
//if(nominal=="TTbar" && HTval>1200.&&METval>30.&&signalbin<=2&&MT2val>100) {cout << "this is ttbar, 2j,   high HT, weight " << weight << " " << MT2val << " " << hh << endl; printtt=true;}
//if(nominal=="TTbar" && HTval>1200.&&METval>30.&&signalbin==5&&MT2val>125) {cout << "this is ttbar, 3b,   high HT, weight " << weight << " " << MT2val << " " << hh << endl; printtt=true;}
//if(nominal=="TTbar" && HTval>1200.&&METval>30.&&signalbin==8&&MT2val>130) {cout << "this is ttbar, 6j2b, high HT, weight " << weight << " " << MT2val << " " << hh << endl; printtt=true;}

	//	if ( fVerbose>2 && counter % 5000 == 0 ) cout << "hh " << hh << endl;
            
            
            Double_t ISRweight(1.); Double_t weightnoISR = weight;
            if(!fMT2tree->misc.isData){
                TLorentzVector hardgenlv; hardgenlv.SetPtEtaPhiM(0,0,0,0);
                if(namesample.Contains("WJets")||namesample.Contains("W-")){
                    if(counter<2) cout << "Is WJets" << endl;
                    bool foundW(false);
                    for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
                        int ID  =abs(fMT2tree->genlept[ngl].ID);
                        int MID =abs(fMT2tree->genlept[ngl].MID);
                        if((ID == 11 || ID == 12 || ID == 13 || ID == 14 || ID == 15 || ID == 16) && MID==24){
                            hardgenlv = fMT2tree->genlept[ngl].Mlv; foundW = true;
                        }
                        if(foundW) break;
                    }
                    if(!foundW){
                        for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
                            int ID  =abs(fMT2tree->genlept[ngl].ID);
                            int MID =abs(fMT2tree->genlept[ngl].MID);
                            int GMID=abs(fMT2tree->genlept[ngl].GMID);
                            if((ID == 11 || ID == 12 || ID == 13 || ID == 14 || ID == 16) && MID==15 && GMID==24){
                                hardgenlv = fMT2tree->genlept[ngl].Mlv; foundW = true;
                            }
                            if(foundW) break;
                        }
                    }
                } if(namesample.Contains("ZJets")){
                    if(counter<2) cout << "Is ZJets" << endl;
                    hardgenlv = fMT2tree->GenZ[0];
                } if(namesample.Contains("TTbar")||namesample.Contains("TT-")){
                    if(counter<2) cout << "Is TTbar" << endl;
                    TLorentzVector top1(0.,0.,0.,0.), top2(0.,0.,0.,0.);
                    bool top1f(false), top2f(false);
                    for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
                        int id   = abs(fMT2tree->genlept[ngl].ID);
                        if(id!=5) continue;
                        int mid  = fMT2tree->genlept[ngl].MID;//from b
                        if(mid==6&&top1f) continue;
                        else if(mid==6) { top1 = fMT2tree->genlept[ngl].Mlv; top1f = true; }
                        //else if(mid==6) { top1.SetPtEtaPhiM(fMT2tree->genlept[n].Mlv.Pt(), fMT2tree->genlept[n].Mlv.Eta(), fMT2tree->genlept[n].Mlv.Phi(), fMT2tree->genlept[n].Mlv.M()); top1f = true; }
                        if(mid==-6&&top2f) continue;
                        else if(mid==-6){ top2 = fMT2tree->genlept[ngl].Mlv; top2f = true; }
                        //else if(mid==-6) { top2.SetPtEtaPhiM(fMT2tree->genlept[n].Mlv.Pt(), fMT2tree->genlept[n].Mlv.Eta(), fMT2tree->genlept[n].Mlv.Phi(), fMT2tree->genlept[n].Mlv.M()); top2f = true; }
                        if(top1f&&top2f) {
                            hardgenlv = top1+top2;
                            break;
                        }
                    }
                } if(namesample.Contains("T_t")||namesample.Contains("T_s")||namesample.Contains("Tbar_t")||namesample.Contains("Tbar_s")||namesample.Contains("T-t")||namesample.Contains("T-s")||namesample.Contains("Tbar-t")||namesample.Contains("Tbar-s")){
                    if(counter<2) cout << "Is Singletop" << endl;
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
                        if(!Wf){//this might be wrong
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
                if(fISRreweight) weight = weight * ISRweight;
            }
            
            

		string helptype = " ";
		//bool nominal = false;
// 		string nominal = " ";
// 		//cout << "MT2" << hh << endl;
// 		//			    Up/Down/...                     TTbar/WJets/...                 l/m/h HT                      2j0b ... 6j2b
// 		//string hs = string("_") + sample_type[i1] + string("_") + sample_kind[i0] + string("_") + HT_region[i3] + string("_") + signal_region[i2];
// // || namesample = "TTbar_MadgraphInclv1" || namesample = "TTbar_MadgraphInclv1" || namesample = "TTbar_Powhegv1" || namesample = "TTbar_Powhegv2" NOT NOMINAL
// 		if(namesample=="TTbar_Madgraph0l" || namesample=="TTbar_Madgraph1l" || namesample=="TTbar_Madgraph2l")                          nominal = "TTbar";
// 		if(namesample=="WJetsHT200"||namesample=="WJetsHT250"||namesample=="WJetsHT300"||namesample=="WJetsHT400")                      nominal = "WJets";
// 		if(namesample=="T_s"||namesample=="Tbar_s"||namesample=="T_t"||namesample=="Tbar_t"||namesample=="T_tW"||namesample=="Tbar_tW") nominal = "SingleTop";
// 		if(namesample=="WW"||namesample=="WZ"||namesample=="ZZ"||namesample=="TTW"||namesample=="TTZ"||namesample=="TTG")               nominal = "TTV";
// 		
		//nominal TTbar
		if(signalbin>=0 && HTbin>=0 && VSPT<70. && J1Pt>100. && J2Pt>100. && MinDPhi>0.3 && MT2val>100.){
		//can put lower cut on MT2 here, as no event modification is made here; only weights change
		//nominal shape, xsec up/down for each process, ISR, TopMass, Matching, Scale
		if(nominal=="TTbar"){
			if(namesample=="TTbar_Madgraph0l") helptype = "_nominal0l";
			if(namesample=="TTbar_Madgraph1l") helptype = "_nominal1l";
			if(namesample=="TTbar_Madgraph2l") helptype = "_nominal2l";
			hh = "_TTbar"+sHT+ssignal;
			if (histos.count("MT2"+helptype     +hh)==0) cout << "MT2"+helptype     +hh << endl;

			histos["MT2"+helptype     +hh]->Fill(MT2val, weight);
//			if(MT2val>100) cout << "MT2_Nominal"      +hh << " should be filled with " << MT2val << " with weight " << weight << ". Before this event integral " << histos["MT2_Nominal"      +hh]->Integral() << endl;
			histos["MT2_Nominal"      +hh]->Fill(MT2val, weight);
//			if(MT2val>100){
//			for(int mm=0; mm<=histos["MT2_Nominal"      +hh]->GetNbinsX();++mm){
//				cout << histos["MT2_Nominal"      +hh]->GetBinLowEdge(mm) << "-" << histos["MT2_Nominal"      +hh]->GetBinLowEdge(mm)+histos["MT2_Nominal"      +hh]->GetBinWidth(mm) << ": " << histos["MT2_Nominal"      +hh]->GetBinContent(mm) << "; " ;}
//			cout << "After this event integral " << histos["MT2_Nominal"      +hh]->Integral() << endl;
//			}
			histos["MT2_TopUp"        +hh]->Fill(MT2val, weight*1.07);
			histos["MT2_TopDown"      +hh]->Fill(MT2val, weight*0.93);
			histos["MT2_SingleTopUp"  +hh]->Fill(MT2val, weight);
			histos["MT2_SingleTopDown"+hh]->Fill(MT2val, weight);
			histos["MT2_WUp"          +hh]->Fill(MT2val, weight);
			histos["MT2_WDown"        +hh]->Fill(MT2val, weight);
			histos["MT2_TTVUp"        +hh]->Fill(MT2val, weight);
			histos["MT2_TTVDown"      +hh]->Fill(MT2val, weight);
			//TopPT up down
			//so far only possible for leptonic decaying W as for rest we miss genid (so both b's and leptonic WJets)
			//procedure here:
			//looking for leptonic top decay --> find matching b-jet
			//				     if lepton is matched to jet (e.g. taus or also other leptons), find that
			//				     assume all met due to leptonic top decay - use met
			// merge them to a top pT --> reweight if: pT > 200 by { 1.-(pT-200)/800. } using AN2012_493_v8
/*
			int lepindtop(-1), bindtop(-1);
			int leptop(-1), btop(-1);
			double lDR(99.), bDR(99.);
			for(int n = 0; n<20; ++n){
				int id   = fMT2tree->genlept[n].ID;
				int mid  = fMT2tree->genlept[n].ID;
				int gmid = fMT2tree->genlept[n].ID;
				if((abs(id)==11||abs(id)==13)&&((abs(mid)==24&&abs(gmid)==6)||(abs(mid)==15&&abs(gmid)==6))) { lepindtop = n; break; }
			}
			if(lepindtop<0){
			for(int n = 0; n<20; ++n){
				int id   = fMT2tree->genlept[n].ID;
				int mid  = fMT2tree->genlept[n].ID;
				int gmid = fMT2tree->genlept[n].ID;
				if(abs(id)==15&&abs(mid)==24&&abs(gmid)==6) lepindtop = n;
			}
			}
			for(int n = 0; n<20; ++n){
				int id   = fMT2tree->genlept[n].ID;
				int mid  = fMT2tree->genlept[n].ID;
				int lepid = 0;
				if(lepindtop>=0) lepid = fMT2tree->genlept[lepindtop].ID;
				if(lepid>0&&mid!=(-6)) continue;
				if(lepid<0&&mid!=  6 ) continue;
				if(abs(id)==5&&abs(mid)==6) bindtop = n;
			}
			for(int n = 0; n<fMT2tree->NJets; ++n){//no jetid cut here, as lepton might miss that
				if(fMT2tree->jet[n].lv.Pt()<20.) continue;
				if(fabs(fMT2tree->jet[n].lv.Eta())>3.0) continue;
				double min = TMath::Min(0.3, lDR);
				if(lepindtop>=0&&fMT2tree->jet[n].lv.DeltaR(fMT2tree->genlept[lepindtop].lv)<min){
					lDR = fMT2tree->jet[n].lv.DeltaR(fMT2tree->genlept[lepindtop].lv); leptop = n; }
				min = TMath::Min(0.3, bDR);
				if(bindtop>=0&&fMT2tree->jet[n].lv.DeltaR(fMT2tree->genlept[bindtop].lv)<min){
					bDR = fMT2tree->jet[n].lv.DeltaR(fMT2tree->genlept[bindtop].lv); btop = n; }
			}
			TLorentzVector l(0.,0.,0.,0.), b(0.,0.,0.,0.), m(0.,0.,0.,0.);
			m.SetPtEtaPhiM(fMT2tree->misc.MET,0.,fMT2tree->misc.METPhi,0.);
			if(lDR<0.3) l.SetPtEtaPhiM(fMT2tree->jet[leptop].lv.Pt(),0.,fMT2tree->jet[leptop].lv.Phi(),0.);
			if(bDR<0.3) l.SetPtEtaPhiM(fMT2tree->jet[btop].lv.Pt(),0.,fMT2tree->jet[btop].lv.Phi(),0.);
			double tpt = (m+l+b).Pt();
			double tweight = 1;
			if(tpt>200.&&tpt<600.) tweight = 1.-(tpt-200.)/800.
			if(tpt>600.)           tweight = 0.25;
*/
			//second option reweight using gentop lv only
			int t1weight(1), t2weight(1);
			bool t1(false), t2(false);
			for(int n = 0; n<20; ++n){
				int mid  = fMT2tree->genlept[n].MID;
				if(mid==6&&t1) continue;
				else if(mid==6){
					double t1pt = fMT2tree->genlept[n].Mlv.Pt();
					if(t1pt>200.&&t1pt<600.) t1weight = 1.-(t1pt-200.)/800.;
					if(t1pt>600.)            t1weight = 0.25;
					t1 = true;
				}
				if(mid==-6&&t2) continue;
				else if(mid==-6){
					double t2pt = fMT2tree->genlept[n].Mlv.Pt();
					if(t2pt>200.&&t2pt<600.) t2weight = 1.-(t2pt-200.)/800.;
					if(t2pt>600.)            t2weight = 0.25;
					t2 = true;
				}
				if(t1&&t2) break;
			}
			histos["MT2_TopPtReweighted"+hh]->Fill(MT2val, weight*t1weight*t2weight);
			//ISR: defined by ttbar
		/*	TLorentzVector top1(0.,0.,0.,0.), top2(0.,0.,0.,0.), ttbar(0.,0.,0.,0.);
			bool top1f(false), top2f(false);
			for(int n = 0; n<20; ++n){
				int id   = abs(fMT2tree->genlept[n].ID);
				if(id!=5) continue;
				int mid  = fMT2tree->genlept[n].MID;//from b
				if(mid==6&&top1f) continue;
				else if(mid==6) { top1.SetPtEtaPhiM(fMT2tree->genlept[n].Mlv.Pt(), fMT2tree->genlept[n].Mlv.Eta(), fMT2tree->genlept[n].Mlv.Phi(), fMT2tree->genlept[n].Mlv.M()); top1f = true; }
				if(mid==-6&&top2f) continue;
				else if(mid==-6) { top2.SetPtEtaPhiM(fMT2tree->genlept[n].Mlv.Pt(), fMT2tree->genlept[n].Mlv.Eta(), fMT2tree->genlept[n].Mlv.Phi(), fMT2tree->genlept[n].Mlv.M()); top2f = true; }
				if(top1f&&top2f) break;
			}
			ttbar = top1 + top2;
			double ISRweight(1.);
			if(ttbar.Pt()>250.) ISRweight = 0.8;
			else if(ttbar.Pt()>150.) ISRweight = 0.9;
			else if(ttbar.Pt()>120.) ISRweight = 0.95;
			histos["MT2_ISR"+hh]->Fill(MT2val, weight*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weight);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weight*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
        */    histos["MT2_ISR"+hh]->Fill(MT2val, weightnoISR*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weightnoISR);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weightnoISR*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
		}
		//nominal W
		if(nominal=="WJets"){
			hh = "_WJets"+sHT+ssignal;
//			if(MT2val>100) cout << "MT2_Nominal"      +hh << " should be filled with " << MT2val << " with weight " << weight << ". Before this event integral " << histos["MT2_Nominal"      +hh]->Integral() << endl;
			histos["MT2_Nominal"      +hh]->Fill(MT2val, weight);
//			if(MT2val>100){
//			for(int mm=0; mm<=histos["MT2_Nominal"      +hh]->GetNbinsX();++mm){
//				cout << histos["MT2_Nominal"      +hh]->GetBinLowEdge(mm) << "-" << histos["MT2_Nominal"      +hh]->GetBinLowEdge(mm)+histos["MT2_Nominal"      +hh]->GetBinWidth(mm) << ": " << histos["MT2_Nominal"      +hh]->GetBinContent(mm) << "; " ;}
//			cout << "After this event integral " << histos["MT2_Nominal"      +hh]->Integral() << endl;
//			}
			histos["MT2_WUp"          +hh]->Fill(MT2val, weight*1.05);
			histos["MT2_WDown"        +hh]->Fill(MT2val, weight*0.95);
			histos["MT2_TopUp"        +hh]->Fill(MT2val, weight);
			histos["MT2_TopDown"      +hh]->Fill(MT2val, weight);
			histos["MT2_SingleTopUp"  +hh]->Fill(MT2val, weight);
			histos["MT2_SingleTopDown"+hh]->Fill(MT2val, weight);
			histos["MT2_TTVUp"        +hh]->Fill(MT2val, weight);
			histos["MT2_TTVDown"      +hh]->Fill(MT2val, weight);
			histos["MT2_MassUp"       +hh]->Fill(MT2val, weight);
			histos["MT2_MassDown"     +hh]->Fill(MT2val, weight);
			//ISR: defined by W (only important W events have neutrinos)
		/*	TLorentzVector W(0.,0.,0.,0.);
			for(int n = 0; n<20; ++n){
				int mid  = fMT2tree->genlept[n].MID;//from neutrino
				if(abs(mid)==24) {
					W.SetPtEtaPhiM(fMT2tree->genlept[n].Mlv.Pt(), fMT2tree->genlept[n].Mlv.Eta(), fMT2tree->genlept[n].Mlv.Phi(), fMT2tree->genlept[n].Mlv.M());
					break;
				}
			}
			double ISRweight(1.);
			if(W.Pt()>250.) ISRweight = 0.8;
			else if(W.Pt()>150.) ISRweight = 0.9;
			else if(W.Pt()>120.) ISRweight = 0.95;
			histos["MT2_ISR"+hh]->Fill(MT2val, weight*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weight);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weight*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
         */
            histos["MT2_ISR"+hh]->Fill(MT2val, weightnoISR*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weightnoISR);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weightnoISR*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)  
		}
		//nominal singleTop
		if(nominal=="SingleTop"){
			hh = "_SingleTop"+sHT+ssignal;
			histos["MT2_Nominal"      +hh]->Fill(MT2val, weight);
			histos["MT2_WUp"          +hh]->Fill(MT2val, weight);
			histos["MT2_WDown"        +hh]->Fill(MT2val, weight);
			histos["MT2_TopUp"        +hh]->Fill(MT2val, weight);
			histos["MT2_TopDown"      +hh]->Fill(MT2val, weight);
			histos["MT2_SingleTopUp"  +hh]->Fill(MT2val, weight*1.1);
			histos["MT2_SingleTopDown"+hh]->Fill(MT2val, weight*0.9);
			histos["MT2_TTVUp"        +hh]->Fill(MT2val, weight);
			histos["MT2_TTVDown"      +hh]->Fill(MT2val, weight);
			histos["MT2_MassUp"       +hh]->Fill(MT2val, weight);
			histos["MT2_MassDown"     +hh]->Fill(MT2val, weight);
			histos["MT2_MatchingUp"   +hh]->Fill(MT2val, weight);
			histos["MT2_MatchingDown" +hh]->Fill(MT2val, weight);
			histos["MT2_ScaleUp"      +hh]->Fill(MT2val, weight);
			histos["MT2_ScaleDown"    +hh]->Fill(MT2val, weight);
			//ISR: defined by top (maybe add W if possible)
	/*		TLorentzVector top1(0.,0.,0.,0.), W(0.,0.,0.,0.), top(0.,0.,0.,0.);
			bool top1f(false), Wf(false);
			for(int n = 0; n<20; ++n){
				int mid  = fMT2tree->genlept[n].MID;//from b
				int gmid = fMT2tree->genlept[n].GMID;//W not from top
				if(abs(mid)==6&&top1f) continue;
				else if(abs(mid)==6) { top1.SetPtEtaPhiM(fMT2tree->genlept[n].Mlv.Pt(), fMT2tree->genlept[n].Mlv.Eta(), fMT2tree->genlept[n].Mlv.Phi(), fMT2tree->genlept[n].Mlv.M()); top1f = true; }
				if(abs(mid)==24&&abs(gmid)!=6&&Wf) continue;
				else if(abs(mid)==24&&abs(gmid)!=6) { W.SetPtEtaPhiM(fMT2tree->genlept[n].Mlv.Pt(), fMT2tree->genlept[n].Mlv.Eta(), fMT2tree->genlept[n].Mlv.Phi(), fMT2tree->genlept[n].Mlv.M()); Wf = true; }
				if(top1f&&Wf) break;
			}
			top = top1 + W;
			double ISRweight(1.);
			if(top.Pt()>250.) ISRweight = 0.8;
			else if(top.Pt()>150.) ISRweight = 0.9;
			else if(top.Pt()>120.) ISRweight = 0.95;
			histos["MT2_ISR"+hh]->Fill(MT2val, weight*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weight);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weight*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
         */   histos["MT2_ISR"+hh]->Fill(MT2val, weightnoISR*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weightnoISR);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weightnoISR*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)

		}
		//nominal rest
		if(nominal=="TTV"){
			hh = "_TTV"+sHT+ssignal;
			histos["MT2_Nominal"      +hh]->Fill(MT2val, weight);
			histos["MT2_WUp"          +hh]->Fill(MT2val, weight);
			histos["MT2_WDown"        +hh]->Fill(MT2val, weight);
			histos["MT2_TopUp"        +hh]->Fill(MT2val, weight);
			histos["MT2_TopDown"      +hh]->Fill(MT2val, weight);
			histos["MT2_SingleTopUp"  +hh]->Fill(MT2val, weight);
			histos["MT2_SingleTopDown"+hh]->Fill(MT2val, weight);
			histos["MT2_TTVUp"        +hh]->Fill(MT2val, weight*1.3);
			histos["MT2_TTVDown"      +hh]->Fill(MT2val, weight*0.7);
			histos["MT2_MassUp"       +hh]->Fill(MT2val, weight);
			histos["MT2_MassDown"     +hh]->Fill(MT2val, weight);
			histos["MT2_MatchingUp"   +hh]->Fill(MT2val, weight);
			histos["MT2_MatchingDown" +hh]->Fill(MT2val, weight);
			histos["MT2_ScaleUp"      +hh]->Fill(MT2val, weight);
			histos["MT2_ScaleDown"    +hh]->Fill(MT2val, weight);
			//ttbar: Define by ttbar, if possible also add V, can also have dibosons
		/*	TLorentzVector top1(0.,0.,0.,0.), top2(0.,0.,0.,0.), V1(0.,0.,0.,0.), V2(0.,0.,0.,0.), ttbar(0.,0.,0.,0.);
			bool top1f(false), top2f(false), V1f(false), V2f(false);
			for(int n = 0; n<20; ++n){
				int mid  = fMT2tree->genlept[n].MID;//from b
				int gmid  = fMT2tree->genlept[n].MID;//from b
				if(mid==6&&top1f) continue;
				else if(mid==6) { top1.SetPtEtaPhiM(fMT2tree->genlept[n].Mlv.Pt(), fMT2tree->genlept[n].Mlv.Eta(), fMT2tree->genlept[n].Mlv.Phi(), fMT2tree->genlept[n].Mlv.M()); top1f = true; }
				if(mid==-6&&top2f) continue;
				else if(mid==-6) { top2.SetPtEtaPhiM(fMT2tree->genlept[n].Mlv.Pt(), fMT2tree->genlept[n].Mlv.Eta(), fMT2tree->genlept[n].Mlv.Phi(), fMT2tree->genlept[n].Mlv.M()); top2f = true; }
				if(((abs(mid)==24&&abs(gmid)!=6)||abs(mid)==23||abs(mid)==22)&&(V1f&&V2f)) continue;
				else if(((abs(mid)==24&&abs(gmid)!=6)||abs(mid)==23||abs(mid)==22)&&V1f==false){ V1.SetPtEtaPhiM(fMT2tree->genlept[n].Mlv.Pt(), fMT2tree->genlept[n].Mlv.Eta(), fMT2tree->genlept[n].Mlv.Phi(), fMT2tree->genlept[n].Mlv.M()); V1f = true; }
				else if(((abs(mid)==24&&abs(gmid)!=6)||abs(mid)==23||abs(mid)==22)){ V2.SetPtEtaPhiM(fMT2tree->genlept[n].Mlv.Pt(), fMT2tree->genlept[n].Mlv.Eta(), fMT2tree->genlept[n].Mlv.Phi(), fMT2tree->genlept[n].Mlv.M()); V2f = true; }
				if(top1f&&top2f&&V1f&&V2f) break;
			}
			ttbar = top1 + top2 + V1 + V2;
			double ISRweight(1.);
			if(ttbar.Pt()>250.) ISRweight = 0.8;
			else if(ttbar.Pt()>150.) ISRweight = 0.9;
			else if(ttbar.Pt()>120.) ISRweight = 0.95;
			histos["MT2_ISR"+hh]->Fill(MT2val, weight*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weight);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weight*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
         */   histos["MT2_ISR"+hh]->Fill(MT2val, weightnoISR*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weightnoISR);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weightnoISR*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
		}
		//scale up/down
		if(namesample=="TT-scaleup"||namesample=="TT-scaledown"||namesample=="W-scaleup"||namesample=="W-scaledown"){
			hh = sHT+ssignal;
			if(namesample=="TT-scaleup"  ) histos["MT2_ScaleUp_TTbar"  +hh]->Fill(MT2val,weight);
			if(namesample=="TT-scaledown") histos["MT2_ScaleDown_TTbar"+hh]->Fill(MT2val,weight);
			if(!WnoScaleMatching && namesample== "W-scaleup"  ) histos["MT2_ScaleUp_WJets"  +hh]->Fill(MT2val,weight);
			if(!WnoScaleMatching && namesample== "W-scaledown") histos["MT2_ScaleDown_WJets"+hh]->Fill(MT2val,weight);
		}
		//mass up/down
		if(namesample=="TT-massup"||namesample=="TT-massdown"){
			hh = "_TTbar"+sHT+ssignal;
			if(namesample=="TT-massup"  ) histos["MT2_MassUp"  +hh]->Fill(MT2val,weight);
			if(namesample=="TT-massdown") histos["MT2_MassDown"+hh]->Fill(MT2val,weight);
		}
		//matching up/down
		if(namesample=="TT-matchingup"||namesample=="TT-matchingdown"||namesample=="W-matchingup"||namesample=="W-matchingdown"){
			hh = sHT+ssignal;
			if(namesample=="TT-matchingup"  ) histos["MT2_MatchingUp_TTbar"  +hh]->Fill(MT2val,weight);
			if(namesample=="TT-matchingdown") histos["MT2_MatchingDown_TTbar"+hh]->Fill(MT2val,weight);
			if(!WnoScaleMatching && namesample== "W-matchingup"  ) histos["MT2_MatchingUp_WJets"  +hh]->Fill(MT2val,weight);
			if(!WnoScaleMatching && namesample== "W-matchingdown") histos["MT2_MatchingDown_WJets"+hh]->Fill(MT2val,weight);
		}
		if(namesample=="WJets_Incl"){//moved up
			hh = "_WJets"+sHT+ssignal;
			if(!WnoScaleMatching) histos["MT2_nominal1Sample"+hh]->Fill(MT2val, weight);
		}
		if(nominal=="WJets" && WnoScaleMatching){//only if not using scale up/down for WJets above due to too small statistics
			hh = sHT+ssignal;
			histos["MT2_ScaleUp_WJets"  +hh]->Fill(MT2val,weight);
			histos["MT2_ScaleDown_WJets"+hh]->Fill(MT2val,weight);
			histos["MT2_MatchingUp_WJets"  +hh]->Fill(MT2val,weight);
			histos["MT2_MatchingDown_WJets"+hh]->Fill(MT2val,weight);
			histos["MT2_nominal1Sample_WJets"+hh]->Fill(MT2val, weight);
		}
		//JES up/down
		if(namesample=="TTbar_Madgraph0l_JESdown" || namesample=="TTbar_Madgraph1l_JESdown" || namesample=="TTbar_Madgraph2l_JESdown") helptype = "TTbar_JESdown";
		if(namesample=="TTbar_Madgraph0l_JESup" || namesample=="TTbar_Madgraph1l_JESup" || namesample=="TTbar_Madgraph2l_JESup") helptype = "TTbar_JESup";
		if(namesample=="WJetsHT200_JESdown"||namesample=="WJetsHT250_JESdown"||namesample=="WJetsHT300_JESdown"||namesample=="WJetsHT400_JESdown") helptype = "WJets_JESdown";
		if(namesample=="WJetsHT200_JESup"||namesample=="WJetsHT250_JESup"||namesample=="WJetsHT300_JESup"||namesample=="WJetsHT400_JESup") helptype = "WJets_JESup";
		if(namesample=="T-s_JESdown"||namesample=="Tbar-s_JESdown"||namesample=="T-t_JESdown"||namesample=="Tbar-t_JESdown"||namesample=="T-tW_JESdown"||namesample=="Tbar-tW_JESdown") helptype = "SingleTop_JESdown";
		if(namesample=="T-s_JESup"||namesample=="Tbar-s_JESup"||namesample=="T-t_JESup"||namesample=="Tbar-t_JESup"||namesample=="T-tW_JESup"||namesample=="Tbar-tW_JESup") helptype = "SingleTop_JESup";
		if(namesample=="WW_JESdown"||namesample=="WZ_JESdown"||namesample=="ZZ_JESdown"||namesample=="TTW_JESdown"||namesample=="TTZ_JESdown"||namesample=="TTG_JESdown") helptype = "TTV_JESdown";
		if(namesample=="WW_JESup"||namesample=="WZ_JESup"||namesample=="ZZ_JESup"||namesample=="TTW_JESup"||namesample=="TTZ_JESup"||namesample=="TTG_JESup") helptype = "TTV_JESup";

		if(helptype=="TTbar_JESdown"||helptype=="TTbar_JESup"||helptype=="WJets_JESdown"||helptype=="WJets_JESup"||helptype=="SingleTop_JESdown"||helptype=="SingleTop_JESup"||helptype=="TTV_JESdown"||helptype=="TTV_JESup"){
			hh = sHT+ssignal;
			if(helptype=="TTbar_JESdown")     histos["MT2_JESDown_TTbar"    +hh]->Fill(MT2val,weight);
			if(helptype=="TTbar_JESup"  )     histos["MT2_JESUp_TTbar"      +hh]->Fill(MT2val,weight);
			if(helptype=="WJets_JESdown")     histos["MT2_JESDown_WJets"    +hh]->Fill(MT2val,weight);
			if(helptype=="WJets_JESup"  )     histos["MT2_JESUp_WJets"      +hh]->Fill(MT2val,weight);
			if(helptype=="SingleTop_JESdown") histos["MT2_JESDown_SingleTop"+hh]->Fill(MT2val,weight);
			if(helptype=="SingleTop_JESup"  ) histos["MT2_JESUp_SingleTop"  +hh]->Fill(MT2val,weight);
			if(helptype=="TTV_JESdown")       histos["MT2_JESDown_TTV"      +hh]->Fill(MT2val,weight);
			if(helptype=="TTV_JESup"  )       histos["MT2_JESUp_TTV"        +hh]->Fill(MT2val,weight);
		}
		//BTV, PU, different TTbar MC
		if(nominal=="TTbar"||nominal=="WJets"||nominal=="SingleTop"||nominal=="TTV"){
			hh = "_"+nominal+sHT+ssignal;
			//BTV
			if(btagSFerr>0 && btagSF!=0){
			histos["MT2_BSFUp"         +hh]->Fill(MT2val, weight*(btagSF+btagSFerr)/btagSF);
			histos["MT2_BSFDown"       +hh]->Fill(MT2val, weight*(btagSF-btagSFerr)/btagSF);
			} else {
			histos["MT2_BSFUp"         +hh]->Fill(MT2val, weight);
			histos["MT2_BSFDown"       +hh]->Fill(MT2val, weight);
			}
			//PUWeight
			int pubin = pu_up->FindBin(fMT2tree->pileUp.PUtrueNumInt);
			double weightup = sample_weight * ISRweight * pu_up->GetBinContent(pubin);
			double weightdown = sample_weight * ISRweight * pu_down->GetBinContent(pubin);
			histos["MT2_PUUp"          +hh]->Fill(MT2val, weightup);
			histos["MT2_PUDown"        +hh]->Fill(MT2val, weightdown);
			if(nominal!="TTbar"){
			if(nominal!="WJets") histos["MT2_nominal1Sample"+hh]->Fill(MT2val, weight);
			histos["MT2_nominalPowheg" +hh]->Fill(MT2val, weight);
			}
		}
		if(namesample=="TTbar_MadgraphIncl"){
			hh = "_TTbar"+sHT+ssignal;
			histos["MT2_nominal1Sample"+hh]->Fill(MT2val, weight);
		}
		if(namesample=="TTbar_Powheg"){
			hh = "_TTbar"+sHT+ssignal;
			histos["MT2_nominalPowheg" +hh]->Fill(MT2val, weight);
		}
		}//if(signalbin>=0 && HTbin>=0 && VSPT<70. && J1Pt>100. && J2Pt>100. && MinDPhi<0.3)
		//unclustered MET
		if(nominal=="TTbar"||nominal=="WJets"||nominal=="SingleTop"||nominal=="TTV"){
			hh = "_"+nominal+sHT+ssignal;
			//unclustered MET
			if(signalbin>=0 && J1Pt>100. && J2Pt>100.&& MinDPhi>0.3){
				//no angular uncertainty
				TLorentzVector METunclusteredup = fMT2tree->pfmet[0] - fMT2tree->MHT[0];
				METunclusteredup.SetPtEtaPhiM(1.1*METunclusteredup.Pt(), METunclusteredup.Eta(), METunclusteredup.Phi(), METunclusteredup.M());
//				if(fabs(METunclusteredup.Pt()-fMT2tree->Vectorsumpt)>7.)
//					cout << "VSPT " << fMT2tree->Vectorsumpt << " new up VSPT " << METunclusteredup.Pt() << endl;
				if(METunclusteredup.Pt()<70){//VSPT cut
					METunclusteredup += fMT2tree->MHT[0];
//					if(fabs(METunclusteredup.Pt()-fMT2tree->misc.MET)>7.)
//						cout << "pfMET " << fMT2tree->misc.MET << " MHT " << fMT2tree->MHT[0].Pt() << " new up MET " << METunclusteredup.Pt() <<endl;
					int uHTbin(-1); string usHT;
					if(HTval>=450.  && HTval<750.  && METunclusteredup.Pt()>200.) { uHTbin = 0; usHT = "_lowHT";    }
					if(HTval>=750.  && HTval<1200. && METunclusteredup.Pt()>30. ) { uHTbin = 1; usHT = "_mediumHT"; }
					if(HTval>=1200.                && METunclusteredup.Pt()>30. ) { uHTbin = 2; usHT = "_highHT";   }
					if(uHTbin>=0){
						hh = "_"+nominal+usHT+ssignal;
						double MT2up = fMT2tree->CalcMT2(0,false,fMT2tree->hemi[0].lv1, fMT2tree->hemi[0].lv2, METunclusteredup);
						histos["MT2_METUp"+hh]->Fill(MT2up, weight);
					}
				}
				TLorentzVector METunclustereddown = fMT2tree->pfmet[0] - fMT2tree->MHT[0];
				METunclustereddown.SetPtEtaPhiM(0.9*METunclustereddown.Pt(), METunclustereddown.Eta(), METunclustereddown.Phi(), METunclustereddown.M());
//				if(fabs(METunclustereddown.Pt()-fMT2tree->Vectorsumpt)>7.)
//					cout << "VSPT " << fMT2tree->Vectorsumpt << " new down VSPT " << METunclustereddown.Pt() << endl;
				if(METunclustereddown.Pt()<70){//VSPT cut
					METunclustereddown += fMT2tree->MHT[0];
//					if(fabs(METunclustereddown.Pt()-fMT2tree->misc.MET)>7.)
//						cout << "pfMET " << fMT2tree->misc.MET << " MHT " << fMT2tree->MHT[0].Pt() << " new down MET " << METunclustereddown.Pt() <<endl;
					int dHTbin(-1); string dsHT;
					if(HTval>=450.  && HTval<750.  && METunclustereddown.Pt()>200.) { dHTbin = 0; dsHT = "_lowHT";    }
					if(HTval>=750.  && HTval<1200. && METunclustereddown.Pt()>30. ) { dHTbin = 1; dsHT = "_mediumHT"; }
					if(HTval>=1200.                && METunclustereddown.Pt()>30. ) { dHTbin = 2; dsHT = "_highHT";   }
					if(dHTbin>=0){
						hh = "_"+nominal+dsHT+ssignal;
						double MT2down = fMT2tree->CalcMT2(0,false,fMT2tree->hemi[0].lv1, fMT2tree->hemi[0].lv2, METunclustereddown);
						histos["MT2_METDown"+hh]->Fill(MT2down, weight);
					}
				}
			}
		}
		//all direct changes

 	}//while(myEvtList->GetEntry(counter++) !=-1)
	delete fMT2tree;
	}//for(size_t i = 0; i < fSamples.size(); ++i)
	cout << "add overflow to last bin" << endl;
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
// 	cout << "normalized ISR histos" << endl;
// 	for(int i1 = 0; i1<samplekindsize;   ++i1){
// 	for(int i2 = 0; i2<signalregionsize; ++i2){
// 	for(int i3 = 0; i3<HTregionsize;     ++i3){
// 		string hs = string("_") + sample_kind[i1] + string("_") + HT_region[i3] + string("_") + signal_region[i2];
// 		string mapnameISR     = "MT2_ISR"     + hs;
// 		string mapnameISRdown = "MT2_ISRDown" + hs;
// 		string mapnameISRup   = "MT2_ISRUp"   + hs;
// 		if(histos[mapnameISRdown]->Integral()>0) histos[mapnameISRdown]->Scale(histos[mapnameISR]->Integral()/histos[mapnameISRdown]->Integral());
// 		if(histos[mapnameISRup  ]->Integral()>0) histos[mapnameISRup  ]->Scale(histos[mapnameISR]->Integral()/histos[mapnameISRup  ]->Integral());
// 	}}}
	cout << "add samples to allMC" << endl;
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		string ht = string("MT2_") + sample_type[i1];
		string hs = string("_") + HT_region[i3] + string("_") + signal_region[i2];
		string mapnameallMC   = ht + string("_allMC")     + hs;
		string mapnameTop     = ht + string("_allTop")    + hs;
		string mapnameTTbar   = ht + string("_TTbar")     + hs;
		string mapnameSingleT = ht + string("_SingleTop") + hs;
		string mapnameTTV     = ht + string("_TTV")       + hs;
		string mapnameWJets   = ht + string("_WJets")     + hs;
		histos[mapnameTop]  ->Add(histos[mapnameTTbar]);
		histos[mapnameTop]  ->Add(histos[mapnameSingleT]);
		histos[mapnameallMC]->Add(histos[mapnameTTbar]);
		histos[mapnameallMC]->Add(histos[mapnameSingleT]);
		histos[mapnameallMC]->Add(histos[mapnameTTV]);
		histos[mapnameallMC]->Add(histos[mapnameWJets]);
	}}}
	cout << "create normalized distribution" << endl;
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		string mapname    = h->first;
		string mapnamenorm = "Norm"+mapname;
		TH1D *hcopy = (TH1D*)h->second->Clone(mapnamenorm.c_str());
		if(hcopy->Integral()>0) hcopy->Scale(1./hcopy->Integral());
		if(histosnormalized.count(mapnamenorm) == 0 ) histosnormalized[mapnamenorm] = hcopy;
	}
	cout << "Saving." << endl;
    	TFile *fsavefile = new TFile(outputdir + outputname,"RECREATE");
	fsavefile->cd();
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Write();
	}
	for(map<string,TH1D*>::iterator h=histosnormalized.begin(); h!=histosnormalized.end();++h){
		h->second->Write();
	}
	fsavefile->Close();
	cout << "Saved histograms in " << outputdir << outputname << endl;
/*
	//do plotting
	vector<string> names; names.clear();
	names.push_back("ComparingDifferentTT");
	names.push_back("ComparingDifferentTTNorm");
	names.push_back("ComparingInclvsComb");
	names.push_back("ComparingInclvsCombNorm");
	names.push_back("ComparingMatchScale");
	names.push_back("ComparingMatchScaleNorm");
	names.push_back("ComparingMasses");
	names.push_back("ComparingMassesNorm");
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(unsigned int i0 = 0; i0<names.size(); ++i0){

	string name = names[i0] + "_" + signal_region[i2];

	TCanvas* c1 = new TCanvas((name+"_ratio").c_str(),"",0,0,1000,1000);
	c1->SetFrameLineWidth(1);
	c1 -> cd();

//cout<<__LINE__<<endl;
	TH1D *h[5];
	int numhistos = 1;
//cout<<__LINE__<<endl;
//cout << "MT2_MadgraphIncl_"<<signal_region[i2] <<" ->count " << histos.count("MT2_MadgraphIncl_"+signal_region[i2]) << endl;
	h[0] = (TH1D*)histos["MT2_MadgraphIncl_"+signal_region[i2]]->Clone(("MadgraphIncl_"+signal_region[i2]).c_str());
//cout<<__LINE__<<endl;
	if(names[i0]=="ComparingDifferentTT" || names[i0]=="ComparingDifferentTTNorm"){
		h[1] = (TH1D*)histos["MT2_MCatNLO_"+signal_region[i2]]->Clone(("MCatNLO_"+signal_region[i2]).c_str());
		h[2] = (TH1D*)histos["MT2_Powheg_"+signal_region[i2]]->Clone(("Powheg_"+signal_region[i2]).c_str());
		numhistos = 3;
	}
//cout<<__LINE__<<endl;
	if(names[i0]=="ComparingInclvsComb" || names[i0]=="ComparingInclvsCombNorm"){
		h[1] = (TH1D*)histos["MT2_MadGraphComb_"+signal_region[i2]]->Clone(("MadGraphComb_"+signal_region[i2]).c_str());
		h[2] = (TH1D*)histos["MT2_MadgraphHad_"+signal_region[i2]]->Clone(("MadgraphHad_"+signal_region[i2]).c_str());
		h[3] = (TH1D*)histos["MT2_MadgraphSemiLep_"+signal_region[i2]]->Clone(("MadgraphSemiLep_"+signal_region[i2]).c_str());
		h[4] = (TH1D*)histos["MT2_MadGraphFullyLep_"+signal_region[i2]]->Clone(("MadGraphFullyLep_"+signal_region[i2]).c_str());
		numhistos = 5;
	}
//cout<<__LINE__<<endl;
	if(names[i0]=="ComparingMatchScale" || names[i0]=="ComparingMatchScaleNorm"){
		h[1] = (TH1D*)histos["MT2_MatchingUp_"+signal_region[i2]]->Clone(("MatchingUp_"+signal_region[i2]).c_str());
		h[2] = (TH1D*)histos["MT2_MatchingDown_"+signal_region[i2]]->Clone(("MatchingDown_"+signal_region[i2]).c_str());
		h[3] = (TH1D*)histos["MT2_ScaleUp_"+signal_region[i2]]->Clone(("ScaleUp_"+signal_region[i2]).c_str());
		h[4] = (TH1D*)histos["MT2_ScaleDown_"+signal_region[i2]]->Clone(("ScaleDown_"+signal_region[i2]).c_str());
		numhistos = 5;
	}
//cout<<__LINE__<<endl;
	if(names[i0]=="ComparingMasses" || names[i0]=="ComparingMassesNorm"){
		h[1] = (TH1D*)histos["MT2_Mass166_"+signal_region[i2]]->Clone(("Mass166_"+signal_region[i2]).c_str());
		h[2] = (TH1D*)histos["MT2_Mass169_"+signal_region[i2]]->Clone(("Mass169_"+signal_region[i2]).c_str());
		h[3] = (TH1D*)histos["MT2_Mass175_"+signal_region[i2]]->Clone(("Mass175_"+signal_region[i2]).c_str());
		h[4] = (TH1D*)histos["MT2_Mass178_"+signal_region[i2]]->Clone(("Mass178_"+signal_region[i2]).c_str());
		numhistos = 5;
	}
//cout<<__LINE__<<endl;
	if(names[i0]=="ComparingDifferentTTNorm" || names[i0]=="ComparingInclvsCombNorm" || names[i0]=="ComparingMatchScaleNorm" || names[i0]=="ComparingMassesNorm"){
		for(int n = 0; n<numhistos; ++n) { h[n]->Scale(1.0/(h[n]->Integral())); h[n]->GetYaxis()->SetTitle("Events (norm.)"); h[n]->SetMinimum(0.001); }
	}
	else{   for(int n = 0; n<numhistos; ++n) h[n]->SetMinimum(0.02); }
//cout<<__LINE__<<endl;

	TLegend *Legend1 = new TLegend(.75,.65,.88,.88);
	Legend1 -> SetFillColor(0);
	Legend1 -> SetBorderSize(0);
	for(int n = 0; n<numhistos; ++n) Legend1->AddEntry(h[n], h[n]->GetName(), "f");
 	TPad *p_plot  = new TPad((name+"_plotpad").c_str(),  "Pad containing the overlay plot", 0,0.211838,1,1 );
 	//p_plot->SetBottomMargin(0.05);
	//p_plot->SetTopMargin(0.09);
	//p_plot->SetLeftMargin(0.1669107);
      	//p_plot->SetRightMargin(0.02);
	p_plot->SetLeftMargin(0.131579);
	p_plot->SetRightMargin(0.08);
	p_plot->SetTopMargin(0.06895515);
	p_plot->SetBottomMargin(0.07206074);
 	p_plot->Draw();
 	TPad *p_rat = new TPad((name+"_ratiopad").c_str(), "Pad containing the ratio",   0,0.01863354,0.9967105,0.2189441);
 	//p_rat->SetTopMargin(0.03);
 	//p_rat->SetBottomMargin(0.05);
	//p_rat->SetRightMargin(0.02);
	p_rat->SetLeftMargin(0.1336634);	
	p_rat->SetRightMargin(0.075);
	p_rat->SetTopMargin(0.06976745);
	p_rat->SetBottomMargin(0.2790698);
 	p_rat->Draw();

	// draw overlay plot
 	p_plot ->cd();

	gPad->SetLogy(1);
	gPad->SetFillStyle(0);
	double max = 0;
	for(int n = 0; n<numhistos; ++n){
		double max1 = h[n]->GetMaximum();
		double max2 = max;
		max  = (max1>max2)?max1:max2;
	}
	max = 100.*max;
	for(int n = 0; n<numhistos; ++n) h[n]->SetMaximum(max);
	h[0]->Draw();
	for(int n = 1; n<numhistos; ++n) h[n]->Draw("same");

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.0305);
	//TString text = fMT2Analysis? "M_{T2} Analysis        ":"";
	//text +=fMT2bAnalysis? "M_{T2}b Analysis     ":"";
	string text = signal_region[i2];
	TitleBox.DrawLatex(0.18,0.943,text.c_str());

 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(Legend1 != NULL ){
		Legend1 -> SetFillColor(0);
		Legend1 -> SetBorderSize(0);
		Legend1 -> Draw();
	} 
	
	// draw the ratio plot
 	p_rat ->cd();
//cout<<__LINE__<<endl;
	//MC with errors
	TH1D *h_ratio[5];
	h_ratio[0] = (TH1D*)h[0]->Clone("h0_copy");
	h_ratio[0]->Divide(h[0]);
	h_ratio[0]->GetYaxis()->SetRangeUser(0.,5.);
	h_ratio[0]->GetXaxis()->SetLabelSize( 0.);
	h_ratio[0]->GetYaxis()->SetTitle("X / MgIncl");
	h_ratio[0]->GetXaxis()->SetTitle(h[0]->GetXaxis()->GetTitle() );
	h_ratio[0]->GetXaxis()->SetTitleSize(0.2);
	h_ratio[0]->GetXaxis()->SetTitleOffset(0.5);
	h_ratio[0]->GetYaxis()->SetLabelSize(0.19);
	h_ratio[0]->GetXaxis()->SetTickLength(0.09);
	h_ratio[0]->GetYaxis()->SetTitleSize(0.18);
	h_ratio[0]->GetYaxis()->SetTitleOffset(0.36);
	h_ratio[0]->GetYaxis()->SetNdivisions(509);

	h_ratio[0]-> SetFillStyle(3004);
	h_ratio[0]-> SetFillColor(kBlack);

	for(int n = 1; n<numhistos; ++n){
	h_ratio[n] = (TH1D*)h[n]->Clone(((string)h[n]->GetName() + "_copy").c_str());	
	h_ratio[n] ->SetStats(0);
	h_ratio[n] ->SetMarkerStyle(20);

 	h_ratio[n] ->Divide(h[n], h[0]);
 	h_ratio[n] ->SetMinimum(0.4);
 	h_ratio[n] ->SetMaximum(3.0);
	h_ratio[n] ->GetYaxis()->SetTitleOffset(h_ratio[n]->GetYaxis()->GetTitleOffset());
	}
	c1->cd();
 	p_rat ->cd();

	h_ratio[0]->SetFillStyle(3001);
	h_ratio[0]->Draw("E2");
	for(int n = 1; n<numhistos; ++n) h_ratio[n] ->Draw("E same");//LEO MOD

//	for(int n = 0; n<numhistos; ++n) cout << h_ratio[n] ->Integral() << " ";//LEO MOD
//	cout << endl;
 
	TLine *l3 = new TLine(h_ratio[0]->GetXaxis()->GetXmin(), 1.00, h_ratio[0]->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	p_rat->Update();
	gPad->RedrawAxis();
 	c1->Update();
//cout << "now ratio should be drawned now " << h_ratio[0]->Integral() << endl;
 	p_rat ->Draw();
 	p_rat ->Draw();
 	c1->Update();
//cout << "now ratio should have been drawned" << endl;
	TString save=name+"_ratio";
	Util::Print(c1, save, outputdir);


//cout<<__LINE__<<endl;

//	delete c1;
//	delete l3;
//	delete Legend1;
//	delete p_plot;
//	delete p_ratio;
	

	}//for(unsigend int i0 = 0; i0<names.size(); ++i0)
	}//for(int i2 = 0; i2<signalregionsize; ++i2)
*/

}//void TTbarStudies



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
