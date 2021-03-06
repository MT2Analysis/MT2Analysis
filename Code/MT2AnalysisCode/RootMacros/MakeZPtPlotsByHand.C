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
#include "THStack.h"
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

//run via root -l -b -q MakeZPtPlotsByHand.C++

using namespace std;

void load(const char* filename = "samples/datasamples/samples_2141_dataonly.dat");
void MakeZPtPlotsByHand();
void Make1DPlotsRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, bool logFlag, bool normalize, TString outname, TLegend *leg, TString xtitle, TString ytitle, float overlayScale);
void MakePlots(map<string, TH1D*> histos, map<string, THStack*> stacks, TLegend *Legend1);

//struct that combines MT2trees with necessary information
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

    const int sampletypesize = 7;//11;
    string sample_type[sampletypesize] = {"QCD", "WJets", "ZJets", "Top", "Other", "mc", "data"};

    const int leptontypesize = 5;//4;
    string lepton_type[leptontypesize] = {"MuMu", "EE", "EMu", "LL", "had"};

const int HTregionsize = 4;
string HT_region[HTregionsize] = {"lowHT", "mediumHT", "highHT", "allHT"};
const int signalregionsize = 4;
string signal_region[signalregionsize] = {"2j", "3to5j", "ge6j", "ge2j"};

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


Bool_t ISR_W      = true;//apply 'ISR recipe' on W+jets, default = true
Bool_t ISR_Z      = true;//apply 'ISR recipe' on Z+jets, default = true
Bool_t ISR_Top    = true;//apply 'ISR recipe' on ttbar+jets, default = true

TString  outputdir = "stupidplots/Z/ISRstudy/";
bool fSave = true;//save plots, default = true

//code creates ZPt spectrum before and after 'ISR reweighting' to compare the influence
//plots both spectra
//at the moment a boosted selection is chosen (in order to compare with singlephoton yield (gamma pt > 180 GeV)
void MakeZPtPlotsByHand(){


	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

    	Util::MakeOutputDir(outputdir);

	TString outputname = "ZllPtHistos.root";

	TString  samples = "samples/samples_2l_noHTdata_filter.dat";//only dummy

	map<string, TH1D*>    histos;
	for(int il = 0; il<leptontypesize;   ++il){
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;     ++i3){

		string hs = string("_") + lepton_type[il] + "_" + signal_region[i2] + "_" + HT_region[i3] + string("_") + sample_type[i1];
		string mapname;
		mapname = "ZPt" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", 20, 150.,750.); //boosted: 20, 150.,750.
		mapname = "ZPt_ISR" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", 20, 150.,750.); //boosted: 20, 150.,750.

	}}}}

	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Sumw2();}

	std::ostringstream cutStream;
	std::ostringstream cutStreamBase;
	cutStream << " " 
		<< "(NEles+NMuons)==2"                   << "&&"
		<< "(NTausIDLoose3Hits)==0"                   << "&&"
	//	<< "ZllRecalculate()"                    << "&&"
		<< "misc.Jet0Pass ==1"                      << "&&"
		<< "misc.Jet1Pass ==1"                      << "&&"
		<< "misc.Vectorsumpt < 70";
	cutStream << "&&misc.HT>450";
	
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
      cutStreamBase << "&&misc.MET/misc.CaloMETRaw<=2.";
	

	TString triggerEE = "(trigger.HLT_DiElectrons&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	TString triggerMM = "(trigger.HLT_DiMuons&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	TString triggerEM = "(trigger.HLT_EMu&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	TString cuts = cutStream.str().c_str();
	TString basecuts = cutStreamBase.str().c_str();

	load(samples.Data());

   	for(size_t i = 0; i < fSamples.size(); ++i){
        
	    string sampletype = (string)fSamples[i].type;
	    if(sampletype==(string)"mc"){
		if(fSamples[i].sname=="QCD")         sampletype = (string)"QCD";
		else if(fSamples[i].sname=="Wtolnu") sampletype = (string)"Other";
		else if(fSamples[i].sname=="TTbar")   sampletype = (string)"Top";
		else if(fSamples[i].name=="TTbar_Madgraph0l")   sampletype = (string)"Top";
		else if(fSamples[i].name=="TTbar_Madgraph1l")   sampletype = (string)"Top";
		else if(fSamples[i].name=="TTbar_Madgraph2l")   sampletype = (string)"Top";
		else if(fSamples[i].sname=="Top")    sampletype = (string)"Top";
//		else if(fSamples[i].sname=="ZJetsToLL_PtZ100")	continue;//nonboosted
		else if(fSamples[i].sname=="ZJetsToLL_PtZ100")	sampletype = (string)"ZJets";//boosted
		else if(fSamples[i].sname=="ZJetsToLL")         continue;//boosted
		else sampletype = (string)"Other";
	    }

	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "          : looping over " << fSamples[i].name << " added in " << sampletype << endl;
	    if(fVerbose>2) cout << "            sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts + "&&" + basecuts;
        
	    if( fSamples[i].type=="data") {
		//leptonic trigger
		if(fSamples[i].sname=="EE-Data") myCuts += " && " + triggerEE + "&&NEles==2&&NMuons==0"; //cuts to be aplied only on data
		if(fSamples[i].sname=="MuMu-Data") myCuts += " && " + triggerMM + "&&NEles==0&&NMuons==2"; //cuts to be aplied only on data
		if(fSamples[i].sname=="EMu-Data") myCuts += " && " + triggerEM + "&&NEles==1&&NMuons==1"; //cuts to be aplied only on data
            }
   	    cout << "Cuts for Flow: " << myCuts << endl;
   	    fSamples[i].tree->Draw(">>selList", myCuts);
        
	    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
	    fSamples[i].tree->SetEventList(myEvtList);
	    int counter=0;
	    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
	    if(myEvtList->GetSize()==0) continue;
	//run over selected events
        while(myEvtList->GetEntry(counter++) !=-1){	
		int jentry = myEvtList->GetEntry(counter-1);
            
		nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
		fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
		
		if ( fVerbose>2 && counter % 5000 == 0  )  cout << "+++ Proccessing event " << counter << endl;
		
		Double_t weight = sample_weight;
		//get ISR weight
		if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 
		Double_t ISRweight(1.); Double_t weightnoISR = weight; Double_t weightISR = weight;
		if(!fMT2tree->misc.isData){
			TLorentzVector hardgenlv; hardgenlv.SetPtEtaPhiM(0,0,0,0);
			if(ISR_W && sampletype=="WJets"){
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
			} if(ISR_Z && sampletype=="ZJets"){
				hardgenlv = fMT2tree->GenZ[0];
			} if(ISR_Top && sampletype=="Top"){
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
			}
			if(hardgenlv.Pt()>250.) ISRweight = 0.8;
			else if(hardgenlv.Pt()>150.) ISRweight = 0.9;
			else if(hardgenlv.Pt()>120.) ISRweight = 0.95;
			else                         ISRweight = 1.;
			weightISR = weight * ISRweight;
		}

		//select leptons within trigger selection
		string slep = "_had";
		if(fMT2tree->NEles==2)         {
			if(fMT2tree->ele[0].lv.Pt()<20) continue;
			if(fMT2tree->ele[1].lv.Pt()<20) continue;
			if(!(fMT2tree->ele[0].IDLoose)) continue;
			if(!(fMT2tree->ele[1].IDLoose)) continue;
			slep = "_EE";
		} else if(fMT2tree->NMuons==2) {
			if(fMT2tree->muo[0].lv.Pt()<20) continue;
			if(fMT2tree->muo[1].lv.Pt()<20) continue;
			slep = "_MuMu";
		} else                         { 
			if(!(fMT2tree->muo[0].lv.Pt()>20||fMT2tree->ele[0].lv.Pt()>20)) continue;
			if(fMT2tree->muo[0].lv.Pt()<20||fMT2tree->ele[0].lv.Pt()<20) continue;
			if(!(fMT2tree->ele[0].IDLoose)) continue;
			slep = "_EMu";
		}
		string sHT;

		fMT2tree->ZllRecalculate();//add leptons to MET and recalculate event variables like MT2
		if(fMT2tree->misc.HT>=450.&&fMT2tree->misc.HT<750.&&fMT2tree->misc.MET>200.) sHT = "_lowHT";
		else if(fMT2tree->misc.HT>=750.&&fMT2tree->misc.HT<1200.&&fMT2tree->misc.MET>30.) sHT = "_mediumHT";
		else if(fMT2tree->misc.HT>=1200.&&fMT2tree->misc.MET>30.) sHT = "_highHT";
		else continue;//for now
		string ssignal;

		//further signal region definition
		if(fMT2tree->NJetsIDLoose40 == 2                                  ) ssignal = "_2j";
		if(fMT2tree->NJetsIDLoose40 >= 3 && fMT2tree->NJetsIDLoose40 <= 5 ) ssignal = "_3to5j";
		if(fMT2tree->NJetsIDLoose40 >= 6                                  ) ssignal = "_ge6j";

		double btagSF(1.), btagSFerr(0.);
		if(!fMT2tree->misc.isData) btagSF = fMT2tree->SFWeight.BTagCSV40ge2; btagSFerr = fMT2tree->SFWeight.BTagCSV40ge2Error;

		string hh = slep+ssignal+sHT+"_"+sampletype;
		if(slep!=""&&sHT!=""&&ssignal!=""){
			if(slep=="_EMu"){
				continue;//dummy, dont du anything
			} else {
				//double Minv = fMT2tree->GetDiLeptonInvMass(0,1,0,10,true);//already in Z_pt
				double Z_pt = fMT2tree->GetDiLeptonPt(0,1,0,10,76,106);
				//if(Minv>=76.&&Minv<=106.){
					histos["ZPt"+hh]->Fill(Z_pt ,weight);//save histograms
					histos["ZPt_ISR"+hh]->Fill(Z_pt ,weightISR);
				//}
			}
		}
	} //while(myEvtList->GetEntry(counter++) !=-1)
	delete fMT2tree;
	delete fSamples[i].tree;

	} //for(size_t i = 0; i < fSamples.size(); ++i)


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

	cout << "add Ele and Muo to Lep" << endl;
	//add ele and muo to lep
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;        ++i3){
		string hse = string("_EE_") + signal_region[i2] + string("_") + HT_region[i3] + string("_") + sample_type[i1];
		string hsem= string("_EMu_") + signal_region[i2] + string("_") + HT_region[i3] + string("_") + sample_type[i1];
		string hsm = string("_MuMu_") + signal_region[i2] + string("_") + HT_region[i3] + string("_") + sample_type[i1];
		string hsl = string("_LL_") + signal_region[i2] + string("_") + HT_region[i3] + string("_") + sample_type[i1];
		histos["ZPt"+hsl]->Add(histos["ZPt"+hse],1);
		histos["ZPt"+hsl]->Add(histos["ZPt"+hsm],1);
		histos["ZPt"+hsl]->Add(histos["ZPt"+hsem],1);
		histos["ZPt_ISR"+hsl]->Add(histos["ZPt_ISR"+hse],1);
		histos["ZPt_ISR"+hsl]->Add(histos["ZPt_ISR"+hsm],1);
		histos["ZPt_ISR"+hsl]->Add(histos["ZPt_ISR"+hsem],1);
	}}}
	cout << "add HT regions" << endl;
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i4 = 0; i4<leptontypesize;   ++i4){
		string hsl = string("_") + lepton_type[i4] + string("_") + signal_region[i2] + string("_") + "lowHT" + string("_") + sample_type[i1];
		string hsm = string("_") + lepton_type[i4] + string("_") + signal_region[i2] + string("_") + "mediumHT" + string("_") + sample_type[i1];
		string hsh = string("_") + lepton_type[i4] + string("_") + signal_region[i2] + string("_") + "highHT" + string("_") + sample_type[i1];
		string hsa = string("_") + lepton_type[i4] + string("_") + signal_region[i2] + string("_") + "allHT" + string("_") + sample_type[i1];
		histos["ZPt"+hsa]->Add(histos["ZPt"+hsl],1);
		histos["ZPt"+hsa]->Add(histos["ZPt"+hsm],1);
		histos["ZPt"+hsa]->Add(histos["ZPt"+hsh],1);
		histos["ZPt_ISR"+hsa]->Add(histos["ZPt_ISR"+hsl],1);
		histos["ZPt_ISR"+hsa]->Add(histos["ZPt_ISR"+hsm],1);
		histos["ZPt_ISR"+hsa]->Add(histos["ZPt_ISR"+hsh],1);
	}}}
	cout << "add NJet regions" << endl;
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i3 = 0; i3<HTregionsize;        ++i3){
	for(int i4 = 0; i4<leptontypesize;   ++i4){
		string hs2 = string("_") + lepton_type[i4] + string("_") + "2j" + string("_") + HT_region[i3] + string("_") + sample_type[i1];
		string hs3 = string("_") + lepton_type[i4] + string("_") + "3to5j" + string("_") + HT_region[i3] + string("_") + sample_type[i1];
		string hs6 = string("_") + lepton_type[i4] + string("_") + "ge6j" + string("_") + HT_region[i3] + string("_") + sample_type[i1];
		string hsa = string("_") + lepton_type[i4] + string("_") + "ge2j" + string("_") + HT_region[i3] + string("_") + sample_type[i1];
		histos["ZPt"+hsa]->Add(histos["ZPt"+hs2],1);
		histos["ZPt"+hsa]->Add(histos["ZPt"+hs3],1);
		histos["ZPt"+hsa]->Add(histos["ZPt"+hs6],1);
		histos["ZPt_ISR"+hsa]->Add(histos["ZPt_ISR"+hs2],1);
		histos["ZPt_ISR"+hsa]->Add(histos["ZPt_ISR"+hs3],1);
		histos["ZPt_ISR"+hsa]->Add(histos["ZPt_ISR"+hs6],1);
	}}}

	cout << "add all samples to mc, etc." << endl;
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;        ++i3){
	for(int i4 = 0; i4<leptontypesize;   ++i4){
		string hsq   = string("_") + lepton_type[i4] + string("_") + signal_region[i2] + string("_") + HT_region[i3] + string("_QCD");
		string hsw   = string("_") + lepton_type[i4] + string("_") + signal_region[i2] + string("_") + HT_region[i3] + string("_WJets");
		string hsz   = string("_") + lepton_type[i4] + string("_") + signal_region[i2] + string("_") + HT_region[i3] + string("_ZJets");
		string hst   = string("_") + lepton_type[i4] + string("_") + signal_region[i2] + string("_") + HT_region[i3] + string("_Top");
		string hso   = string("_") + lepton_type[i4] + string("_") + signal_region[i2] + string("_") + HT_region[i3] + string("_Other");
		string hsm   = string("_") + lepton_type[i4] + string("_") + signal_region[i2] + string("_") + HT_region[i3] + string("_mc");
			histos[("ZPt"+hsm)]->Add(histos[("ZPt"+hsq)],  1);
			histos[("ZPt"+hsm)]->Add(histos[("ZPt"+hsw)],  1);
			histos[("ZPt"+hsm)]->Add(histos[("ZPt"+hsz)],  1);
			histos[("ZPt"+hsm)]->Add(histos[("ZPt"+hst)],  1);
			histos[("ZPt"+hsm)]->Add(histos[("ZPt"+hso)],  1);
			histos[("ZPt_ISR"+hsm)]->Add(histos[("ZPt_ISR"+hsq)],  1);
			histos[("ZPt_ISR"+hsm)]->Add(histos[("ZPt_ISR"+hsw)],  1);
			histos[("ZPt_ISR"+hsm)]->Add(histos[("ZPt_ISR"+hsz)],  1);
			histos[("ZPt_ISR"+hsm)]->Add(histos[("ZPt_ISR"+hst)],  1);
			histos[("ZPt_ISR"+hsm)]->Add(histos[("ZPt_ISR"+hso)],  1);
	}}}

	cout << "Saving." << endl;
    	TFile *fsavefile = new TFile(outputdir + outputname,"RECREATE");
	fsavefile->cd();
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Write();
	}
	fsavefile->Close();
	cout << "Saved histograms in " << outputdir << outputname << endl;

	//prepare legend for plotting
    TLegend *Legend1 = new TLegend(.71,.71,.91,.92);
	Legend1 -> SetFillColor(0);
   	Legend1 -> SetBorderSize(0);
	for(int il = 0; il<leptontypesize;   ++il){
	for(int is = 0; is<sampletypesize;   ++is){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		string hs = string("_") + lepton_type[il] + string("_") + signal_region[i2] + string("_") + HT_region[i3] + string("_") + sample_type[is];
		Color_t colour = 603;		     if(sample_type[is]=="QCD")       colour = 401;
		else if(sample_type[is]=="WJets")     colour = 417;
		else if(sample_type[is]=="ZJets")     colour = 419;
		else if(sample_type[is]=="Top")       colour = 855;
		else if(sample_type[is]=="Other")     colour = 603;
		else if(sample_type[is]=="mc")        colour = 603;
		if(sample_type[is]!="data"){
			histos["ZPt"+hs] ->SetFillColor(colour);
			histos["ZPt"+hs] ->SetLineColor(colour);
			histos["ZPt_ISR"+hs] ->SetFillColor(colour);
			histos["ZPt_ISR"+hs] ->SetLineColor(colour);
		} else{
			histos["ZPt"+hs]->SetLineColor(kBlack);
			histos["ZPt"+hs]->SetMarkerStyle(20);
			histos["ZPt"+hs]->SetMarkerColor(kBlack);
			histos["ZPt_ISR"+hs]->SetLineColor(kBlack);
			histos["ZPt_ISR"+hs]->SetMarkerStyle(20);
			histos["ZPt_ISR"+hs]->SetMarkerColor(kBlack);
		}
	}}}}
	//prepare stack for plotting
    map<string, THStack*> stacks;
	bool leggy = true;
	for(unsigned int i3 = 0; i3<HTregionsize; ++i3){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	   	for(int il = 0; il<leptontypesize; ++il){
			string h = string("_") + lepton_type[il] + string("_") + signal_region[i2] + string("_") + HT_region[i3];
			if(stacks.count("ZPt"+h)==0) stacks["ZPt"+h ] = new THStack(("ZPt" +h).c_str(), ("ZPt" +h).c_str());
			if(stacks.count("ZPt_ISR"+h)==0) stacks["ZPt_ISR"+h ] = new THStack(("ZPt_ISR" +h).c_str(), ("ZPt_ISR" +h).c_str());
			//stack filling
		   	for(int is = 0; is<sampletypesize; ++is){
				if(sample_type[is]=="data") continue;
				if(sample_type[is]=="mc")   continue;
				string hs = string("_") + lepton_type[il] + string("_") + signal_region[i2] + string("_") + HT_region[i3] + string("_") + sample_type[is];
				if(histos.count("ZPt"+hs)!=0){ stacks["ZPt"+h] ->Add(histos["ZPt"+hs]);}
				if(histos.count("ZPt_ISR"+hs)!=0){ stacks["ZPt_ISR"+h] ->Add(histos["ZPt_ISR"+hs]);}
			}
		}

	}
	}
	if(leggy){
   	for(int is = 0; is<sampletypesize; ++is){
		string hs = string("_") + string("LL_ge2j_allHT") + string("_") + sample_type[is];
		if(sample_type[is]=="mc") continue;
		if(sample_type[is]!="data" && sample_type[is]!="QCD" && sample_type[is]!="WJets") Legend1->AddEntry(histos["ZPt"+hs], (sample_type[is]).c_str(), "f");
		else if(sample_type[is]=="data") Legend1->AddEntry(histos["ZPt"+hs], (sample_type[is]).c_str(), "lp");
		leggy = false;
	}
	}
	MakePlots(histos, stacks, Legend1);//do the plotting

}

//function to invoke Make1DPlotsRatio(...)
void MakePlots(map<string, TH1D*> histos, map<string, THStack*> stacks, TLegend *Legend1){

	cout << "Plotting 1D histograms" << endl;

	for(int i0=0;i0<leptontypesize;++i0){
	for(int i1=0;i1<HTregionsize;++i1){
	for(int i2=0;i2<signalregionsize;++i2){

	string name = "ZPt";
	string h = "_"+lepton_type[i0]+"_"+signal_region[i2]+"_"+HT_region[i1];
	string hd = "_"+lepton_type[i0]+"_"+signal_region[i2]+"_"+HT_region[i1]+"_data";
	string hm = "_"+lepton_type[i0]+"_"+signal_region[i2]+"_"+HT_region[i1]+"_mc";
	TString ytitle = "Events"; TString xtitle = "Z_{ll}-p_{T} [GeV]";
	TString outname = name + h;
	if(histos[name+hm]->Integral()>0||histos[name+hd]->Integral()>0)
		Make1DPlotsRatio(stacks[name+h], histos[name+hm], histos[name+hd], true, false, outname, Legend1, xtitle, ytitle, 1.);
	name = "ZPt_ISR";
	outname = name + h;
	if(histos[name+hm]->Integral()>0||histos[name+hd]->Integral()>0)
		Make1DPlotsRatio(stacks[name+h], histos[name+hm], histos[name+hd], true, false, outname, Legend1, xtitle, ytitle, 1.);
	}}}

}

//copy out of MassPlotter.cc - no detail comments
void Make1DPlotsRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, bool logFlag, bool normalize, TString outname, TLegend *leg, TString xtitle, TString ytitle, float overlayScale){

	// define canvas and pads 
	TH1D *h1 = (TH1D*)histmc->Clone("h1_copy");
	TH1D *h2 = (TH1D*)histdata->Clone("h2_copy");

	h1->SetTitle("");
	h2->SetTitle("");
	hstack->SetTitle("");

	h1->SetStats(0);	
	h2->SetStats(0);	
	h1->SetMarkerStyle(20);
	h2->SetMarkerStyle(20);

	TCanvas* c1 = new TCanvas(outname+"_ratio","",0,0,600,600);
	c1->SetFrameLineWidth(1);
	c1 -> cd();

 	TPad *p_plot  = new TPad(outname+"_plotpad",  "Pad containing the overlay plot", 0,0.211838,1,1 );

	p_plot->SetLeftMargin(0.131579);
	p_plot->SetRightMargin(0.08);
	p_plot->SetTopMargin(0.06895515);
	p_plot->SetBottomMargin(0.07206074);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad(outname+"_ratiopad", "Pad containing the ratio",   0,0.01863354,0.9967105,0.2189441);
	p_ratio->SetLeftMargin(0.1336634);	
	p_ratio->SetRightMargin(0.075);
	p_ratio->SetTopMargin(0.06976745);
	p_ratio->SetBottomMargin(0.2790698);

	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	if(logFlag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logFlag) max = 2.5*max;
	else max = 1.5*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	hstack->SetMaximum(max);

	hstack->SetMinimum(0.02);
	hstack->Draw("hist");
	hstack->GetYaxis()->SetTitle(ytitle);
	hstack->GetYaxis()->SetLabelSize(0.05);
	hstack->GetYaxis()->SetTitleSize(0.05);
	hstack->GetYaxis()->SetTitleOffset(1.3);

	h2    ->Draw("sameE");

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.03);
	TString text ="";
	text = outname;
	TitleBox.DrawLatex(0.13,0.943,text.Data());

 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
	//gPad->SetLogy();
 
	TH1D *h_ratio = (TH1D*)histdata->Clone("h2_copy_2");	
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);

 	h_ratio ->Divide(h2, h1);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	//MC with errors
	TH1D*h_ratio_mc = (TH1D*)histmc->Clone("h1_copy_2");
	h_ratio_mc->Divide(h1);
	h_ratio_mc->GetYaxis()->SetRangeUser(0,2);
	h_ratio_mc->GetXaxis()->SetLabelSize( 0.);
	h_ratio_mc->GetYaxis()->SetTitle("Data / MC");
	h_ratio_mc->GetXaxis()->SetTitle(xtitle);
	h_ratio_mc->GetXaxis()->SetTitleSize(0.2);
	h_ratio_mc->GetXaxis()->SetTitleOffset(0.5);
	h_ratio_mc->GetYaxis()->SetLabelSize(0.19);
	h_ratio_mc->GetXaxis()->SetTickLength(0.09);
	h_ratio_mc->GetYaxis()->SetTitleSize(0.18);
	h_ratio_mc->GetYaxis()->SetTitleOffset(0.36);
	h_ratio_mc->GetYaxis()->SetNdivisions(509);//xxxnew

	h_ratio   ->SetTitle("");
	h_ratio_mc->SetTitle("");
	
	h_ratio_mc->SetFillStyle(3001);
	h_ratio_mc->Draw("E2");
	h_ratio ->DrawCopy("Esame");//LEO MOD
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=outname+"_ratio";
	if(fSave)Util::Print(c1, save, outputdir);
}

//standard function to load samples.dat
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


