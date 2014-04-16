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

//run via root -l -b -q ZnunuNumbers.C++

void load(const char* filename = "samples/datasamples/samples_2141_dataonly.dat");
void ZnunuNumbers();//extendable to all MC numbers if wanted

//struct that combines MT2tree.hh and necessary information like cross section
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
const int fVerbose = 3; 
TString fPath;

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

const int sampletypesize = 8;
string sample_type[sampletypesize] = {"QCD", "WJets", "ZJets", "Top", "Other", "mc", "susy", "data"};
const int signalregionsize = 9;
string signal_region[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b"};
const int HTbinsize = 3;
string HT_bin[HTbinsize] = {"lowHT", "mediumHT", "highHT"};

//this is a dummy code to produce Z(nunu) numbers from MC for >=2 bjet prediction
//and stores them into a root file - dummy code in a sense that it has correct histogram naming
//such that it is easy to load these numbers for result plots/tables and interpreation
void ZnunuNumbers(){

	bool fbTagReweight  = true; //default = true, reweight MC due to BTV SF weights
	bool TOBTECreweight = true; //default = true, reweight MC due to tobtec filter efficiency, this will be obsolete from CMSSW 6X series on

  	gROOT->ProcessLine(".x SetStyle_PRD.C");
	TString fOutDir                    = "../Results/Filtered/GammaJetsPrediction/20130617_test/";
	TString outputname = "ZnunuNumbers.root";
	TString samples = "samples/samples_Znunu_HTMET_filter.dat";
	if(TOBTECreweight){
		outputname = "ZnunuNumbers_TOBTECreweight.root";
		samples = "samples/samples_Znunu_HTMET.dat";
	}
	//get tobtec reweighting map
	TFile *_file0 = TFile::Open("/shome/haweber/TOBTECfiles/PassingEventsStudy2_NoMinDPhiNoLeptonVetoIncludingPhotons.root");
	TEfficiency *effmap = (TEfficiency*)_file0->Get("SJetPt100Etaeff"); 

	//load histograms
	map<string, TH1D*>    histos;
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTbinsize;        ++i3){
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
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", NMT2bins, MT2bins);
	}}
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Sumw2();}

	//signal event selection
std::ostringstream  fCutStreamSignal;
	fCutStreamSignal << " " 
	  << "misc.MT2>=100"                                               << "&&"
	  << "misc.MET>=30"                                                << "&&"
	  << "misc.HT >=450"                                               << "&&"
	  << "NEles==0"                                                    << "&&"
	  << "NMuons==0"                                                   << "&&"
	  << "NTausIDLoose3Hits==0"                                        << "&&"
	  << "misc.Jet0Pass ==1"                                           << "&&"
	  << "misc.Jet1Pass ==1"                                           << "&&"
	  << "misc.PassJet40ID ==1"                                        << "&&"
	  << "NJetsIDLoose40 >=2"                                          << "&&"
	  << "misc.MinMetJetDPhi4Pt40 >0.3"                                << "&&"
	  << "misc.Vectorsumpt < 70"                                       << "&&"
	  // Noise
	  << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"             << "&&" // for rare SM samples
	  << "misc.CSCTightHaloIDFlag == 0"                                << "&&"
	  << "misc.trackingFailureFlag==0"                                 << "&&"
	  << "misc.eeBadScFlag==0"                                         << "&&"
	  << "misc.EcalDeadCellTriggerPrimitiveFlag==0"                    << "&&"
	  << "misc.TrackingManyStripClusFlag==0"                           << "&&"
	  << "misc.TrackingTooManyStripClusFlag==0"                        << "&&"
	  << "misc.TrackingLogErrorTooManyClustersFlag==0"                 << "&&"
          << "misc.CrazyHCAL==0";
	fCutStreamSignal << "&&((misc.MET>=200&&misc.MT2>=200&&misc.HT<750)||(misc.HT>750))";
	fCutStreamSignal << "&&NJetsIDLoose40>=2";
	fCutStreamSignal << "&&misc.MET/misc.CaloMETRaw<=2.";
	TString cuts = fCutStreamSignal.str().c_str();
	//no trigger for MC

	load(samples.Data());

  	for(size_t i = 0; i < fSamples.size(); ++i){

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
	if(sampletype!="ZJets") continue;//run only over Z events

	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "ZnunuNumbers: looping over " << fSamples[i].name << " added in " << sampletype << endl;
	    if(fVerbose>2) cout << "              sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts;// + "&&" + basecuts;
        
        
   	    cout << "Cuts for Flow: " << myCuts << endl;
   	    fSamples[i].tree->Draw(">>selList", myCuts);
        
	    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
	    fSamples[i].tree->SetEventList(myEvtList);
	    int counter=0;
	    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
	    if(myEvtList->GetSize()==0) continue;
	//run over all events
        while(myEvtList->GetEntry(counter++) !=-1){	
		int jentry = myEvtList->GetEntry(counter-1);
            
		nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
		fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
		
		if ( fVerbose>2 && counter % 5000 == 0  )  cout << "+++ Proccessing event " << counter << endl;
		
		Double_t weight = sample_weight;
		if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 
		Double_t ISRweight(1.); Double_t weightnoISR = weight;
		//decided to not do ISR reweighting here. It was not clear to us if ISRweights are valid in b-enriched region.
		/*if(fISRreweight && !fMT2tree->misc.isData){
			TLorentzVector hardgenlv; hardgenlv.SetPtEtaPhiM(0,0,0,0);
			if(sampletype=="WJets"){
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
			weight = weight * ISRweight;
		}*/

		if(TOBTECreweight){
			float ptmax(-1), etamin(99); int nind(-1);
			for(int n = 0; n<fMT2tree->NJets; ++n){
				if(!(fMT2tree->jet[n].isPFIDLoose))     continue;
				if(fMT2tree->jet[n].lv.Pt()<100)        continue;
				if(fabs(fMT2tree->jet[n].lv.Eta())>2.4) continue;
				float deltaetamin = fabs(1.5-etamin);
				float deltaeta = fabs(1.5-fabs(fMT2tree->jet[n].lv.Eta()));
				float eta = fabs(fMT2tree->jet[n].lv.Eta());
				float pt = fMT2tree->jet[n].lv.Pt();
				if(     deltaeta< 0.2&&deltaetamin>=0.2) {ptmax = pt; etamin = eta; nind = n;}
				else if(deltaeta< 0.2&&deltaetamin< 0.2) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }
				else if(deltaeta>=0.2&&deltaetamin>=0.2){
				if(     deltaeta< 0.4&&deltaetamin>=0.4) {ptmax = pt; etamin = eta; nind = n;}
				else if(deltaeta< 0.4&&deltaetamin< 0.4) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }
				else if(deltaeta>=0.4&&deltaetamin>=0.4){
				if(     deltaeta< 0.6&&deltaetamin>=0.6) {ptmax = pt; etamin = eta; nind = n;}
				else if(deltaeta< 0.6&&deltaetamin< 0.6) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }
				else if(deltaeta>=0.6&&deltaetamin>=0.6){
				if(     deltaeta< 0.8&&deltaetamin>=0.8) {ptmax = pt; etamin = eta; nind = n;}
				else if(deltaeta< 0.8&&deltaetamin< 0.8) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }
				else if(deltaeta>=0.8&&deltaetamin>=0.8){
				if(     deltaeta< 1.0&&deltaetamin>=1.0) {ptmax = pt; etamin = eta; nind = n;}
				else if(deltaeta< 1.0&&deltaetamin< 1.0) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }
				else if(deltaeta>=1.0&&deltaetamin>=1.0){
				if(     deltaeta< 1.2&&deltaetamin>=1.2) {ptmax = pt; etamin = eta; nind = n;}
				else if(deltaeta< 1.2&&deltaetamin< 1.2) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }
				else if(deltaeta>=1.2&&deltaetamin>=1.2){
				if(     deltaeta< 1.4&&deltaetamin>=1.4) {ptmax = pt; etamin = eta; nind = n;}
				else if(deltaeta< 1.4&&deltaetamin< 1.4) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }
				else if(deltaeta>=1.4&&deltaetamin>=1.4){
				if(     deltaeta< 1.6&&deltaetamin>=1.6) {ptmax = pt; etamin = eta; nind = n;}
				else if(deltaeta< 1.6&&deltaetamin< 1.6) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }  }}}}}}}
			}
			
			int bin = effmap->FindFixBin(ptmax,etamin);
			float eff = effmap->GetEfficiency(bin);
			float efflow = effmap->GetEfficiencyErrorLow(bin);
			float effup = effmap->GetEfficiencyErrorUp(bin); 

			if(eff>0.) weight = weight * eff;
		}

		//get HT x topological region
		string sHT;
		if(fMT2tree->misc.HT<450.)      sHT = "_HTge0";
		else if(fMT2tree->misc.HT<750.) sHT = "_lowHT";
		else if(fMT2tree->misc.HT<1200.)sHT = "_mediumHT";
		else                            sHT = "_highHT";

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
	//		if(!fbTagError) btagSFerr = 0;//don't store btag error is equivalent to no error at all
		}
		string hh = sHT + ssignal;
		//fill histograms
		histos[(string)"MT2"        + hh]->Fill(fMT2tree->misc.MT2, weight);

	}//while
	delete fMT2tree;
	delete fSamples[i].tree;

	}//for samples

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

	cout << "Saving." << endl;
    	TFile *fsavefile = new TFile(fOutDir + outputname,"RECREATE");
	fsavefile->cd();
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Write();
	}
	fsavefile->Close();
	cout << "Saved histograms in " << fsavefile->GetName() << endl;

}

//reads in samples.dat
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
