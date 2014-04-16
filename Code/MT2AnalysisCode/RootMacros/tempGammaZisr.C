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

//run via root -l -b -q tempGammaZisr.C++

using namespace std;

void load(const char* filename = "samples/datasamples/samples_2141_dataonly.dat");
void tempGammaZisr();

//combines MT2trees with necessary information like cross section
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

const int gNMT2bins_2j0b_MET                    = 12;
double  gMT2bins_2j0b_MET[gNMT2bins_2j0b_MET+1] = {100, 160, 200, 220, 240, 260, 280, 300, 330, 370, 425, 600, 900};
const int gNMT2bins_2j1b_MET                    = 9;
double  gMT2bins_2j1b_MET[gNMT2bins_2j1b_MET+1] = {100, 160, 200, 220, 240, 260, 290, 340, 440, 600};
const int gNMT2bins_3j0b_MET                    = 9;
double  gMT2bins_3j0b_MET[gNMT2bins_3j0b_MET+1] = {100, 160, 200, 225, 260, 300, 360, 450, 560, 700};
const int gNMT2bins_3j1b_MET                    = 7;
double  gMT2bins_3j1b_MET[gNMT2bins_3j1b_MET+1] = {100, 160, 200, 230, 275, 340, 450, 650};
const int gNMT2bins_3j2b_MET                    = 5;
double  gMT2bins_3j2b_MET[gNMT2bins_3j2b_MET+1] = {100, 160, 200, 240, 320, 450};
const int gNMT2bins_6j0b_MET                    = 6;
double  gMT2bins_6j0b_MET[gNMT2bins_6j0b_MET+1] = {100, 160, 200, 240, 280, 350, 500};
const int gNMT2bins_6j1b_MET                    = 5;
double  gMT2bins_6j1b_MET[gNMT2bins_6j1b_MET+1] = {100, 160, 200, 230, 280, 500};
const int gNMT2bins_6j2b_MET                    = 3;
double  gMT2bins_6j2b_MET[gNMT2bins_6j2b_MET+1] = {100, 200, 260, 400};
const int gNMT2bins_3b_MET                      = 3;
double  gMT2bins_3b_MET[gNMT2bins_3b_MET+1]     = {100, 200, 260, 400};


double ISRweights[4] = {1.00,0.95,0.90,0.80};
double ISupperpts[4] = {120.,150.,250.,99999.};


const int sampletypesize = 4;
string sample_type[sampletypesize] = {"ISRUp", "ISR", "ISRDown", "ISRdivNoISR"};//maybe not all filled  - ! include TTbar and WJets
const int samplekindsize = 2;
string sample_kind[samplekindsize] = {"G","Z"};
const int HTregionsize = 4;
string HT_region[HTregionsize] = {"lowHT", "mediumHT", "highHT", "mhHT"};
const int signalregionsize = 9;
string signal_region[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b"};

//this function runs the 'ISR recipe' (incl. up/down systematics) on the signal MT2 distributions
//via gen-Z and gen-photon information for Z(nunu) and single photon data sample
//then the ratio MT2(ISR reweighted) / MT2 (nominal) is built
//that histogram will be used for effective 'ISR reweighting' for run_GammaJetsToZnunu.C
void tempGammaZisr(){


	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	TString  outputdir = "";
	TString outputname = "ZGammaISRCorrectionsHistos.root";

	TString  samplesGlowHT  = "samples/samples_1g_GEst_MET.dat";//only dummy
	TString  samplesGhighHT = "samples/samples_1g_GEst_HT.dat";//only dummy
	TString  samplesZlowHT  = "samples/samples_1g_GEst_MET.dat";//only dummy
	TString  samplesZhighHT = "samples/samples_1g_GEst_HT.dat";//only dummy
	map<string, TH1D*>    histos;
	map<string, TH1D*>    histosnormalized;
	for(int i0 = 0; i0<samplekindsize;   ++i0){
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		if(i3==3) continue;
		int NMT2bins;
		if(i3==0){
			if(signal_region[i2]=="2j0b")    NMT2bins = gNMT2bins_2j0b_MET;
			if(signal_region[i2]=="2j1to2b") NMT2bins = gNMT2bins_2j1b_MET;
			if(signal_region[i2]=="3to5j0b") NMT2bins = gNMT2bins_3j0b_MET;
			if(signal_region[i2]=="3to5j1b") NMT2bins = gNMT2bins_3j1b_MET;
			if(signal_region[i2]=="3to5j2b") NMT2bins = gNMT2bins_3j2b_MET;
			if(signal_region[i2]=="6j0b")    NMT2bins = gNMT2bins_6j0b_MET;
			if(signal_region[i2]=="6j1b")    NMT2bins = gNMT2bins_6j1b_MET;
			if(signal_region[i2]=="6j2b")    NMT2bins = gNMT2bins_6j2b_MET;
			if(signal_region[i2]=="3b")      NMT2bins = gNMT2bins_3b_MET;
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
			if(signal_region[i2]=="2j0b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_2j0b_MET[in]; }
			if(signal_region[i2]=="2j1to2b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_2j1b_MET[in]; }
			if(signal_region[i2]=="3to5j0b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3j0b_MET[in]; }
			if(signal_region[i2]=="3to5j1b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3j1b_MET[in]; }
			if(signal_region[i2]=="3to5j2b") { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3j2b_MET[in]; }
			if(signal_region[i2]=="6j0b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_6j0b_MET[in]; }
			if(signal_region[i2]=="6j1b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_6j1b_MET[in]; }
			if(signal_region[i2]=="6j2b")    { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_6j2b_MET[in]; }
			if(signal_region[i2]=="3b")      { for(int in = 0; in<=NMT2bins; ++in) MT2bins[in] = gMT2bins_3b_MET[in];   }
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

		string hs = string("_") + sample_type[i1] + string("_") + sample_kind[i0] + string("_") + HT_region[i3] + string("_") + signal_region[i2];
		string mapname = "MT2" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", NMT2bins, MT2bins);
		Color_t colour = 1;
		histos[mapname]->SetLineColor(colour); histos[mapname]->SetMarkerColor(colour); histos[mapname]->SetLineWidth(2);
		histos[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); histos[mapname]->GetYaxis()->SetTitle("Events"); 

	}}}}

	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		cout << h->first << "    " << h->second->GetName() << endl;
		h->second->Sumw2();}

	std::ostringstream cutStream;
	std::ostringstream cutStreamBase;
	cutStream << " " 
		<< "NJetsIDLoose>=2"                        << "&&"//preselection only
		<< "NTausIDLoose3Hits==0"                   << "&&"
		<< "NEles==0"                               << "&&"
		<< "NMuons==0"                              << "&&"
		<< "misc.Jet0Pass ==1"                      << "&&"
		<< "misc.Jet1Pass ==1"                      << "&&"
		<< "misc.Vectorsumpt < 70";//precut
	cutStream  << "&& misc.MinMetJetDPhi4Pt40 >0.3";
	cutStream << "&&((misc.HT>450&&misc.HT<750&&misc.MET>200)||(misc.HT>750&&misc.MET>30))";
 //       cutStream << "&&NPhotons==1&&photon[0].isLooseID";

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
      << "misc.CrazyHCAL==0"
      << "&&misc.MET/misc.CaloMETRaw<=2.";

	TString cuts = cutStream.str().c_str();
	TString basecuts = cutStreamBase.str().c_str();

	load(samplesGlowHT.Data());
	//this is the photon selection - low HT
   	for(size_t i = 0; i < fSamples.size(); ++i){

		if(fSamples[i].sname!="Photons") continue;

	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "SampleName: looping over " << fSamples[i].name <<  endl;
	    if(fVerbose>2) cout << "            sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts + "&&" + basecuts + "&&NPhotons==1&&photon[0].isLooseID&&photon[0].lv.Pt()>180";
        
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
		if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight*fMT2tree->SFWeight.BTagCSV40eq0; // pile-up reweighting for MC 

		int NumJets   = fMT2tree->NJetsIDLoose40;
		int NumBJets  = fMT2tree->NBJets40CSVM;
		float HTval   = fMT2tree->misc.HT;
		float METval  = fMT2tree->misc.MET;
		float MT2val  = fMT2tree->misc.MT2;
		int signalbin(-1); string ssignal;//note: this can change, depending on redefinitions
		if(NumJets == 2 &&                 NumBJets == 0) { signalbin = 0; ssignal = "_2j0b";    }
		if(NumJets >= 3 && NumJets <= 5 && NumBJets == 0) { signalbin = 2; ssignal = "_3to5j0b"; }
		if(NumJets >= 6 &&                 NumBJets == 0) { signalbin = 6; ssignal = "_6j0b";    }
		int HTbin(-1); string sHT;
		if(HTval>=450.  && HTval<750.  && METval>200.) { HTbin = 0; sHT = "_lowHT";    }
		if(HTval>=750.  && HTval<1200. && METval>30. ) { HTbin = 1; sHT = "_mediumHT"; }
		if(HTval>=1200.                && METval>30. ) { HTbin = 2; sHT = "_highHT";   }
		
		if(NumJets>=2&&NumBJets>=0&&HTval>=450.&&METval>=200.){
		string hh = "_G_lowHT_ge2j";

		if(fMT2tree->GenPhoton[0].Pt()>0){
			double ISRweight(1.);
			if(fMT2tree->GenPhoton[0].Pt()>250.) ISRweight = 0.8;
			else if(fMT2tree->GenPhoton[0].Pt()>150.) ISRweight = 0.9;
			else if(fMT2tree->GenPhoton[0].Pt()>120.) ISRweight = 0.95;
			//store all lowHT histos
			for(int tt =0; tt<signalregionsize; ++tt){
			hh = "_G_lowHT_"+signal_region[tt];
			histos["MT2_ISR"+hh]->Fill(MT2val, weight*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weight);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weight*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
			}
		}
		}
	}
	delete fMT2tree;
	}
	load(samplesGhighHT.Data());
	//this is the photon selection - Medium and high HT
   	for(size_t i = 0; i < fSamples.size(); ++i){

		if(fSamples[i].sname!="Photons") continue;

	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "SampleName: looping over " << fSamples[i].name <<  endl;
	    if(fVerbose>2) cout << "            sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts + "&&" + basecuts + "&&NPhotons==1&&photon[0].isLooseID&&photon[0].lv.Pt()>20";
        
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
		if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight*fMT2tree->SFWeight.BTagCSV40eq0; // pile-up reweighting for MC 

		int NumJets   = fMT2tree->NJetsIDLoose40;
		int NumBJets  = fMT2tree->NBJets40CSVM;
		float HTval   = fMT2tree->misc.HT;
		float METval  = fMT2tree->misc.MET;
		float MT2val  = fMT2tree->misc.MT2;
		int signalbin(-1); string ssignal;//note: this can change, depending on redefinitions
		if(NumJets == 2 &&                 NumBJets == 0) { signalbin = 0; ssignal = "_2j0b";    }
		if(NumJets >= 3 && NumJets <= 5 && NumBJets == 0) { signalbin = 2; ssignal = "_3to5j0b"; }
		if(NumJets >= 6 &&                 NumBJets == 0) { signalbin = 6; ssignal = "_6j0b";    }
		int HTbin(-1); string sHT;
		if(HTval>=450.  && HTval<750.  && METval>200.) { HTbin = 0; sHT = "_lowHT";    }
		if(HTval>=750.  && HTval<1200. && METval>30. ) { HTbin = 1; sHT = "_mediumHT"; }
		if(HTval>=1200.                && METval>30. ) { HTbin = 2; sHT = "_highHT";   }

		if(NumJets>=2&&NumBJets>=0&&HTval>750.&&METval>30.){
		string hh = "_G_mhHT_ge2j";

		if(fMT2tree->GenPhoton[0].Pt()>0){
			double ISRweight(1.);
			if(fMT2tree->GenPhoton[0].Pt()>250.) ISRweight = 0.8;
			else if(fMT2tree->GenPhoton[0].Pt()>150.) ISRweight = 0.9;
			else if(fMT2tree->GenPhoton[0].Pt()>120.) ISRweight = 0.95;
			//store all lowHT histos
			for(int tt =0; tt<signalregionsize; ++tt){
			hh = "_G_mediumHT_"+signal_region[tt];
			histos["MT2_ISR"+hh]->Fill(MT2val, weight*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weight);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weight*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
			hh = "_G_highHT_"+signal_region[tt];
			histos["MT2_ISR"+hh]->Fill(MT2val, weight*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weight);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weight*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
			}
		}
		}
	}
	delete fMT2tree;
	}


//Z

	load(samplesZlowHT.Data());
	//this is the Z selection - low HT
   	for(size_t i = 0; i < fSamples.size(); ++i){

		if(fSamples[i].shapename!="ZJetsToNuNu") continue;

	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "SampleName: looping over " << fSamples[i].name <<  endl;
	    if(fVerbose>2) cout << "            sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts + "&&" + basecuts;

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
		if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight*fMT2tree->SFWeight.BTagCSV40eq0; // pile-up reweighting for MC 

		int NumJets   = fMT2tree->NJetsIDLoose40;
		int NumBJets  = fMT2tree->NBJets40CSVM;
		float HTval   = fMT2tree->misc.HT;
		float METval  = fMT2tree->misc.MET;
		float MT2val  = fMT2tree->misc.MT2;
		int signalbin(-1); string ssignal;//note: this can change, depending on redefinitions
		if(NumJets == 2 &&                 NumBJets == 0) { signalbin = 0; ssignal = "_2j0b";    }
		if(NumJets >= 3 && NumJets <= 5 && NumBJets == 0) { signalbin = 2; ssignal = "_3to5j0b"; }
		if(NumJets >= 6 &&                 NumBJets == 0) { signalbin = 6; ssignal = "_6j0b";    }
		int HTbin(-1); string sHT;
		if(HTval>=450.  && HTval<750.  && METval>200.) { HTbin = 0; sHT = "_lowHT";    }
		if(HTval>=750.  && HTval<1200. && METval>30. ) { HTbin = 1; sHT = "_mediumHT"; }
		if(HTval>=1200.                && METval>30. ) { HTbin = 2; sHT = "_highHT";   }
		
		if(NumJets>=2&&NumBJets>=0&&HTval>=450.&&METval>=200.){
		string hh = "_Z_lowHT_ge2j";
		if(fMT2tree->GenZ[0].Pt()>0){
			double ISRweight(1.);
			if(fMT2tree->GenZ[0].Pt()>250.) ISRweight = 0.8;
			else if(fMT2tree->GenZ[0].Pt()>150.) ISRweight = 0.9;
			else if(fMT2tree->GenZ[0].Pt()>120.) ISRweight = 0.95;
			//store all lowHT histos
			for(int tt =0; tt<signalregionsize; ++tt){
			hh = "_Z_lowHT_"+signal_region[tt];
			histos["MT2_ISR"+hh]->Fill(MT2val, weight*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weight);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weight*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
			}
		}
		}
	}
	delete fMT2tree;
	}
	load(samplesZhighHT.Data());
	//this is the Z selection - medium and high HT
   	for(size_t i = 0; i < fSamples.size(); ++i){

		if(fSamples[i].shapename!="ZJetsToNuNu") continue;

	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "SampleName: looping over " << fSamples[i].name <<  endl;
	    if(fVerbose>2) cout << "            sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts + "&&" + basecuts;
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
		if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight*fMT2tree->SFWeight.BTagCSV40eq0; // pile-up reweighting for MC 

		int NumJets   = fMT2tree->NJetsIDLoose40;
		int NumBJets  = fMT2tree->NBJets40CSVM;
		float HTval   = fMT2tree->misc.HT;
		float METval  = fMT2tree->misc.MET;
		float MT2val  = fMT2tree->misc.MT2;
		int signalbin(-1); string ssignal;//note: this can change, depending on redefinitions
		if(NumJets == 2 &&                 NumBJets == 0) { signalbin = 0; ssignal = "_2j0b";    }
		if(NumJets >= 3 && NumJets <= 5 && NumBJets == 0) { signalbin = 2; ssignal = "_3to5j0b"; }
		if(NumJets >= 6 &&                 NumBJets == 0) { signalbin = 6; ssignal = "_6j0b";    }
		int HTbin(-1); string sHT;
		if(HTval>=450.  && HTval<750.  && METval>200.) { HTbin = 0; sHT = "_lowHT";    }
		if(HTval>=750.  && HTval<1200. && METval>30. ) { HTbin = 1; sHT = "_mediumHT"; }
		if(HTval>=1200.                && METval>30. ) { HTbin = 2; sHT = "_highHT";   }

		if(NumJets>=2&&NumBJets>=0&&HTval>750.&&METval>30.){
		string hh = "_Z_mhHT_ge2j";
		if(fMT2tree->GenZ[0].Pt()>0){
			double ISRweight(1.);
			if(fMT2tree->GenZ[0].Pt()>250.) ISRweight = 0.8;
			else if(fMT2tree->GenZ[0].Pt()>150.) ISRweight = 0.9;
			else if(fMT2tree->GenZ[0].Pt()>120.) ISRweight = 0.95;
			//store all lowHT histos
			for(int tt =0; tt<signalregionsize; ++tt){
			hh = "_Z_mediumHT_"+signal_region[tt];
			histos["MT2_ISR"+hh]->Fill(MT2val, weight*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weight);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weight*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
			hh = "_Z_highHT_"+signal_region[tt];
			histos["MT2_ISR"+hh]->Fill(MT2val, weight*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(MT2val, weight);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(MT2val, weight*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
			}
		}
		}
	}
	delete fMT2tree;
	}

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
	cout << "make new vs. old ratio" << endl;
	//that histogram will be used for effective 'ISR reweighting' for run_GammaJetsToZnunu.C
	for(int i1 = 0; i1<samplekindsize;   ++i1){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		if(i3==3) continue;
		string hs = string("_") + sample_kind[i1] + string("_") + HT_region[i3] + string("_") + signal_region[i2];
		string mapnameISR     = "MT2_ISR"     + hs;
		string mapnameISRnew  = "MT2_ISRdivNoISR" + hs;
		string mapnameISRup   = "MT2_ISRUp"   + hs;
		if(histos[mapnameISRup  ]->Integral()>0) histos[mapnameISRnew  ]->Divide(histos[mapnameISR],histos[mapnameISRup  ]);
	}}}
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
