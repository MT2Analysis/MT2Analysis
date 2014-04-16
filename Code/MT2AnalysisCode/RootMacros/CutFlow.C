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
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path of MT2tree.hh
//I think only one of the bottom two includes where needed in order to compute the BTV SF weights
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/MT2Code/include/Utilities.hh"


using namespace std;

//call this macro via root -l -b -q CutFlow.C++
void load(const char* filename = "samples_2141_dataonly.dat");
void GetContent2DRange(TH2D* num, double &eff, double &err, double xlow, double xup, double ylow, double yup);
void GetContent2DRangeBin(TH2D* num, double &eff, double &err, int xlow, int xup, int ylow, int yup);
inline void PrintHeadFoot(Bool_t head, Bool_t foot, Bool_t HTline, Bool_t HTlow, Bool_t headerline, Bool_t SUSY);
void GetCutFlowLine(map<string, TH2D*> histograms, Bool_t withuncertainty, Bool_t SUSY, double MT2low, double MT2up, double HTlow, double HTup);

//This code is very old, but may serve as a template.
//Basically it takes simulation samples and data and plot either
//the cutflow for bins chosen
//or just the event count per bin chosen
//where cut flow means the expected event counts for each MC background type and data
void CutFlow();

//the MT2 sample
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

//Flags - these are for the 7 TeV analysis
Bool_t  lowMT2                    = false;//Cutflow for MT2b analysis (7 TeV)
Bool_t  highMT2                   = true; //Cutflow for MT2 analysis (7 TeV)
Bool_t  lowHT                     = false;//Cutflow only for the low HT region (7 TeV, i.e. 750<HT<950 GeV)
Bool_t  highHT                    = false;//Cutflow only for the high HT region (7 TeV, i.e. HT>950 GeV)
Bool_t  allHT                     = true; //Cutflow for both HT regions (7 TeV)
Bool_t  dofastestimate            = false;//do an MT2 cut from the beginning (starting with lowest MT2 that belongs to signal region
Bool_t  withuncert                = false;//Plot also the uncertainties in the cutflow table (default = false), note that this implementation had not been done, i.e. errors are only statistical
Bool_t  susytable                 = false;//Plot also a possible susy signal yield in the cutflow table (default = false)
TString samples; 			  //loads the samples you want to be in the cutflow (this is hardcoded below)
//this is obsolete, but it loads the b-tagging efficiencies for calculating the BTV SF.
TString btagging_hadfile          = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/Efficiencies/BEffHistos_PTbinned_allHT_SSVHPT.root";
TString tagger                    = "SSVHPT";


void CutFlow(){

    //hardcoded samples.dat
    samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/";
    if(lowHT)             samples += "samples_4400_noMinDPhi.dat";
    if(highHT)            samples += "samples_4600_noMinDPhi.dat";
    if(highMT2 && lowHT)  samples += "samples_4400_noMinDPhi_MT2.dat";
    if(highMT2 && allHT)  samples += "samples_4400_noMinDPhi_MT2.dat";
    if(highMT2 && highHT) samples += "samples_4400_noMinDPhi_MT2_highHT.dat";
    if(lowMT2  && lowHT)  samples += "samples_4400_noMinDPhi_MT2b.dat";
    if(lowMT2  && allHT)  samples += "samples_4400_noMinDPhi_MT2b.dat";
    if(lowMT2  && highHT) samples += "samples_4400_noMinDPhi_MT2b_highHT.dat";

    const int sampletypesize = 11;
    const int sampletypesizemc = 5;//first 5 samples only used for cutflow table
    string sample_type[sampletypesize] = {"QCD", "WJets", "ZJetsToNuNu", "Top", "Other", "mc", "LM6", "LM7", "LM8", "LM9", "data"};

    //Use 2D histograms in HT and MT2
    map<string, TH2D*> histos;

    for(int is = 0; is<sampletypesize; ++is){
	string hs = string("_") + sample_type[is];
	string mapname;
	mapname = "MT2" + hs;
        if(histos.count(mapname) == 0 ) histos[mapname] = new TH2D(mapname.c_str(), mapname.c_str(), 700, 0., 3500., 125, 750., 7000.);          
    }

	for(map<string,TH2D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Sumw2();}
	//get b-/c-/light-tagging efficiencies
	TFile* btagefffile;
	TH1D*  hbeff;
	TH1D*  hceff;
	TH1D*  hleff;
	if(lowMT2){
		btagefffile	= TFile::Open(btagging_hadfile);
		hbeff		= (TH1D*)btagefffile->Get("h_beff");
		hceff		= (TH1D*)btagefffile->Get("h_ceff");
		hleff		= (TH1D*)btagefffile->Get("h_leff");
	}
    

    std::ostringstream cutStream;
	cutStream << " " 
    << "misc.MET>=30"                                                  << "&&"
    << "misc.HT >= 750 "                                               << "&&"//trigger cut
    << "NJetsIDLoose40 >= 3"                                           << "&&"
    << "misc.Jet0Pass==1"                                              << "&&"
    << "misc.Jet1Pass==1"                                              << "&&"
    << "misc.SecondJPt  >100"                                          << "&&"
    << "misc.PassJetID ==1"                                            << "&&"
    << "misc.Vectorsumpt < 70"                                         << "&&"
    << "(NMuons==0 || muo[0].lv.Pt()<10)"                              << "&&"
    << "(NEles==0  || ele[0].lv.Pt()<10)"                              << "&&"
    // Noise
    << "misc.HBHENoiseFlagIso==0"                                      << "&&"
    << "misc.CSCTightHaloID==0"                                        << "&&"
    << "misc.CrazyHCAL==0";
    if(dofastestimate){
        if(lowMT2){         cutStream << "&&misc.MT2>=125";//lowest MT2 for MT2b
        } else if(highMT2){ cutStream << "&&misc.MT2>=150"; }//lowest MT2 for MT2
    } if(lowMT2){
        cutStream << "&&"
        << "misc.LeadingJPt  >150"                                     << "&&"
        << "NJetsIDLoose40 >=4"                                        << "&&"
        << "misc.MinMetJetDPhi4 >0.3"                                  << "&&"
        << "NBJets >0";
    } if(highMT2){
        cutStream << "&&"
        << "NJetsIDLoose40 >=3"                                        << "&&"
        << "misc.MinMetJetDPhi >0.3";
    } if(highHT){ cutStream << "&&misc.HT>=950";
    } if(lowHT){  cutStream << "&&misc.HT<950&&misc.HT>=750";  }
	TString cuts = cutStream.str().c_str();
   
	std::ostringstream triggerStream;
	triggerStream << "( "
	<< "(trigger.HLT_HT440_v2 ==1 && misc.Run<161216)" << "||"
	<< "(trigger.HLT_HT450_v2 ==1 && (misc.Run>=161216 && misc.Run< 163269))" << "||"
	<< "(trigger.HLT_HT500_v3 ==1 && (misc.Run>=163269 && misc.Run<=163869))" << "||"
	<< "(trigger.HLT_HT500_v4 ==1 && (misc.Run>=165088 && misc.Run< 165970))" << "||"
	<< "(trigger.HLT_HT550_v5 ==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" << "||"
	<< "(trigger.HLT_HT550_v6 ==1 && (misc.Run==166346))" << "||"
	<< "(trigger.HLT_HT550_v7 ==1 && (misc.Run>=167078 && misc.Run< 170249))" << "||"
	<< "(trigger.HLT_HT550_v8 ==1 && (misc.Run>=170249 && misc.Run< 173236))" << "||"
	<< "(trigger.HLT_HT600_v1 ==1 && (misc.Run>=173236 && misc.Run< 178420))" << "||"
        << "(trigger.HLT_HT650_v4 ==1 && (misc.Run>=178420 && misc.Run< 179959))" << " )";

	TString trigger = triggerStream.str().c_str();

	//load the samples --> go into fSamples
	load(samples.Data());

    for(size_t i = 0; i < fSamples.size(); ++i){

	    string sampletype = (string)fSamples[i].type;
	    //define the background types
	    if(sampletype==(string)"mc"){
		if(fSamples[i].sname=="QCD") sampletype = (string)"QCD";
		else if(fSamples[i].sname=="Wtolnu") sampletype = (string)"WJets";
		else if(fSamples[i].shapename=="ZJetsToNuNu") sampletype = (string)"ZJetsToNuNu";
		else if(fSamples[i].sname=="Top") sampletype = (string)"Top";//no ttbar, includes TTZ, TTW
		else sampletype = (string)"Other";
	    }
	    if(sampletype==(string)"susy") sampletype = (string)fSamples[i].sname;
        
	    //get global event weight
	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "CutFlow:    looping over " << fSamples[i].name << " added in " << sampletype <<  endl;
	    if(fVerbose>2) cout << "            sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

	    //Define the MT2tree
	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts;
        
	    if( fSamples[i].type=="data") myCuts += " && " + trigger; //cuts to be aplied only on data
        
   	    cout << "Cuts for Flow: " << myCuts << endl;
   	    fSamples[i].tree->Draw(">>selList", myCuts);
        
	    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
	    fSamples[i].tree->SetEventList(myEvtList);
	    int counter=0;
	    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
	    if(myEvtList->GetSize()==0) continue;
	//run only over those events passing event selection
        while(myEvtList->GetEntry(counter++) !=-1){	
      		int jentry = myEvtList->GetEntry(counter-1);
            
            nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
            fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
            
            if ( fVerbose>2 && jentry % 50000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;
            
            Double_t weight = sample_weight;
            if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 

	//REWEIGHT DUE TO BTV SF -> ONLY  FOR MT2b
	//This code is obsolete for >=8 TeV analyses
	    //need to be initialized so in case of highMT2 nothing happens
	    float SFweightErr = 0;
	    float SFweight = 1;//outside, since need it there later
	    float SFweightup = 1;
	    float SFweightdown = 1;
	    if(lowMT2 && (!fMT2tree->misc.isData)){
		vector<float> jetEff;
		vector<float> jetEffErr;
		vector<float> jetSF;
		vector<float> jetSFErr;
		vector<float> jetEffup;
		vector<float> jetEffdown;
		vector<float> jetSFup;
		vector<float> jetSFdown;
		int njetsusuable = 0;
		for(int n = 0; n<fMT2tree->NJets; ++n){
			if(fMT2tree->jet[n].isPFIDLoose==false) continue;
			if(fMT2tree->jet[n].Flavour<=-7777) continue;//safety cut for uninitialized jetFlavour
			float effPt  = fMT2tree->jet[n].lv.Pt();
			float effEta = fabs(fMT2tree->jet[n].lv.Eta());
			++njetsusuable;
			if(abs(fMT2tree->jet[n].Flavour)==5){//b-jets
				jetEff.push_back( float(hbeff->GetBinContent(hbeff->FindBin(effPt))) );
				jetEffErr.push_back( float(hbeff->GetBinError(hbeff->FindBin(effPt))) );
				float SFErr;
				float SF = getBTagSF(SFErr, tagger, effPt, effEta);
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr );
			}
			else if(abs(fMT2tree->jet[n].Flavour)==4){//c-jets
				jetEff.push_back( float(hceff->GetBinContent(hceff->FindBin(effPt))) );
				jetEffErr.push_back( float(hceff->GetBinError(hceff->FindBin(effPt))) );
				float SFErr;
				float SF = getBTagSF(SFErr, tagger, effPt, effEta );
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr*2.);
			}
			else {//light-jets
				jetEff.push_back( float(hleff->GetBinContent(hleff->FindBin(effPt))) );
				jetEffErr.push_back( float(hleff->GetBinError(hleff->FindBin(effPt))) );
				float SFErr;
				float SF = getMistagSF(SFErr, tagger, effPt, effEta, 0 );
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr );
			}
		}
		SFweight = getBTagEventWeightError(SFweightErr, jetEff, jetEffErr, jetSF, jetSFErr, -1);
		if(SFweight==0){
			if(njetsusuable!=0){
				cout << "Event has zero weight, do not use it" << endl;
				continue;
			}
			else { //event has no flavour information, use average event weight - this never happened and was just a safety measure in case b-tagging file had empty histograms
				SFweight = 0.945572;
				SFweightErr = sqrt(0.0257166*0.0257166 + 0.0370919+0.0370919);
			}
		}
		weight  = weight * SFweight;
	    }
		//fill histograms
		histos[(string)"MT2_" + sampletype]->Fill(fMT2tree->misc.MT2, fMT2tree->misc.HT, weight);
	}//while
	}//for sample
	//add all MC sampes to MT2_mc
	for(int ists = 0; ists<sampletypesizemc; ++ists){
		histos[string("MT2_mc")]->Add(histos[string("MT2_") + sample_type[ists]], 1.);
	}

	PrintHeadFoot(true, false, false, false, true, susytable);//prints the "head" of Latex table
	PrintHeadFoot(false, false,true, true, false, susytable);//prints line for "low HT"
	GetCutFlowLine(histos, withuncert, susytable, 0., 3500., 750., 950.);//first line was without any MT2 cut... unless flag was set, then it was sum of all MT2 signal bins only
	if(lowMT2){
	GetCutFlowLine(histos, withuncert, susytable, 125., 150., 750., 950.);
	GetCutFlowLine(histos, withuncert, susytable, 150., 200., 750., 950.);
	GetCutFlowLine(histos, withuncert, susytable, 200., 300., 750., 950.);
	GetCutFlowLine(histos, withuncert, susytable, 300., 3500., 750., 950.);
	}
	else if(highMT2){
	GetCutFlowLine(histos, withuncert, susytable, 150., 200., 750., 950.);
	GetCutFlowLine(histos, withuncert, susytable, 200., 275., 750., 950.);
	GetCutFlowLine(histos, withuncert, susytable, 275., 375., 750., 950.);
	GetCutFlowLine(histos, withuncert, susytable, 375., 500., 750., 950.);
	GetCutFlowLine(histos, withuncert, susytable, 500., 3500., 750., 950.);
	}
	PrintHeadFoot(false, false, true, false, false, susytable);//prints line for "high HT"
	GetCutFlowLine(histos, withuncert, susytable, 0., 3500., 950., 7000.);
	if(lowMT2){
	GetCutFlowLine(histos, withuncert, susytable, 125., 150., 950., 7000.);
	GetCutFlowLine(histos, withuncert, susytable, 150., 180., 950., 7000.);
	GetCutFlowLine(histos, withuncert, susytable, 180., 260., 950., 7000.);
	GetCutFlowLine(histos, withuncert, susytable, 260., 3500., 950., 7000.);
	}
	else if(highMT2){
	GetCutFlowLine(histos, withuncert, susytable, 150., 200., 950., 7000.);
	GetCutFlowLine(histos, withuncert, susytable, 200., 275., 950., 7000.);
	GetCutFlowLine(histos, withuncert, susytable, 275., 375., 950., 7000.);
	GetCutFlowLine(histos, withuncert, susytable, 375., 500., 950., 7000.);
	GetCutFlowLine(histos, withuncert, susytable, 500., 3500., 950., 7000.);
	}
	PrintHeadFoot(false, true, false, false, false, susytable);//prints the "foot" of Latex table

	cout << endl;

}

//print "head", "foot" and "HT line" for table
inline void PrintHeadFoot(Bool_t head, Bool_t foot, Bool_t HTline, Bool_t HTlow, Bool_t headerline, Bool_t SUSY){

	if(head){
		cout << "\\begin{table}[htb!]" << endl << "\\begin{center}" << endl;
		cout << "\\caption{Expected background event yields and observed number of events in data ...ADD SOMETHINH HERE... in the low and high $H_T$ regions.}" << endl;
		cout << "\\label{table:ADDCORRECTLABELHERE}" << endl;
		cout << "\\begin{tabular}{lcccccccccccccccc}" << endl;
		cout << "\\hline\\hline" << endl;
	}
	if(headerline && SUSY==false) cout << " & QCD& $W+$jets& Top & $Z(\\nu\\nu)+$jets& Other& MC& Data \\\\ \\hline \\hline" << endl;
	else if(headerline && SUSY)   cout << " & LM6 & LM7 & LM8 & LM9 \\\\ \\hline \\hline" << endl;
	if(HTline && HTlow)           cout << " $ 750 \\leq H_T \\leq 950$ GeV \\\\ \\hline" << endl;
	if(HTline && HTlow==false)    cout << " \\hline $H_T \\geq 950$ GeV \\\\ \\hline" << endl;

	if(foot) cout << " \\hline\\hline " << endl << " \\end{tabular} " << endl << " \\end{center} " << endl << " \\end{table}" << endl;

}

//prints one line in cutflow table using previously filled histograms
void GetCutFlowLine(map<string, TH2D*> histograms, Bool_t withuncertainty, Bool_t SUSY, double MT2low, double MT2up, double HTlow, double HTup){

	map<string, TH2D*> hists = histograms;
	if(MT2low==0. && MT2up==3500.) cout << " All Selections      & ";
	else if(MT2up==3500.)          cout << " $M_{T2}$ (" << int(MT2low) << ", $\\infty$] GeV & ";
	else                           cout << " $M_{T2}$ (" << int(MT2low) << ", " << int(MT2up) << "] GeV & ";

	double QCD, QCDerr, W, Werr, Z, Zerr, Top, Toperr, Other, Othererr, MC, MCerr;
	double data, dataerr, LM6, LM6err, LM7, LM7err, LM8, LM8err, LM9, LM9err;

	if(SUSY==false){
		GetContent2DRange(hists[string("MT2_QCD")],       QCD,   QCDerr, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("MT2_WJets")],       W,     Werr, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("MT2_ZJetsToNuNu")], Z,     Zerr, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("MT2_Top")],       Top,   Toperr, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("MT2_Other")],   Other, Othererr, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("MT2_mc")],         MC,    MCerr, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("MT2_data")],     data,  dataerr, MT2low, MT2up, HTlow, HTup);
	
		if(withuncertainty){
			cout << fixed << setprecision(2) << QCD << " \\pm " << QCDerr << " & " << W << " \\pm " << Werr << " & " << Top << " \\pm " << Toperr << " & " << Z << " \\pm " << Zerr << " & ";
			cout << fixed << setprecision(2) << Other << " \\pm " << Othererr << " & " << MC << " \\pm " << MCerr << " & ";
			cout << (int)data << " \\pm " << fixed << setprecision(2) << dataerr << " \\\\ " << endl;
		}
		else{
			cout << fixed << setprecision(2) << QCD << " & " << W << " & " << Top << " & " << Z << " & " << Other << " & " << MC << " & ";
			cout << (int)data << " \\\\ " << endl;
		}
	}
	else{
		GetContent2DRange(hists[string("MT2_LM6")],       LM6,   LM6err, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("MT2_LM7")],       LM7,   LM7err, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("MT2_LM8")],       LM8,   LM8err, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("MT2_LM9")],       LM9,   LM9err, MT2low, MT2up, HTlow, HTup);

	
		if(withuncertainty){
			cout << fixed << setprecision(2) << LM6 << " \\pm " << LM6err << " & " << LM7 << " \\pm " << LM7err << " & " << LM8 << " \\pm " << LM8err << " & " << LM9 << " \\pm " << LM9err << " & " << " \\\\ " << endl;
		}
		else{
			cout << fixed << setprecision(2) << LM6 << " & " << LM7 << " & " << LM8 << " & " << LM9  << " \\\\ " << endl;
		}
	}
}

//this function actually only finds the bins of the 2D histogram that should be selected for the event count printout
void GetContent2DRange(TH2D* num, double &eff, double &err, double xlow, double xup, double ylow, double yup){
    int xbinlow, xbinup, ybinlow, ybinup;
    //change by 10e-6 to avoid border effects where one bin is included just as it is the lower border of a bin, so it in fact should not be included
    xbinlow = num->GetXaxis()->FindBin(xlow+10e-6);
    xbinup  = num->GetXaxis()->FindBin(xup -10e-6);
    ybinlow = num->GetYaxis()->FindBin(ylow+10e-6);
    ybinup  = num->GetYaxis()->FindBin(yup -10e-6);
    GetContent2DRangeBin(num, eff, err, xbinlow, xbinup, ybinlow, ybinup);
}
//this selects the necessary histogram bins selected for printout
//and stores them in a separate histogram that gets rebinned afterwards 
//--> thus I make sure that the statistical error is calculated correctly
void GetContent2DRangeBin(TH2D* num, double &eff, double &err, int xlow, int xup, int ylow, int yup){
    TH2D *num_c;
    num_c = (TH2D*) num->Clone();
    vector<double> xx, yy;
    for(int i = xlow; i<=xup; ++i){
        xx.push_back(num_c->GetXaxis()->GetBinLowEdge(i));}
    xx.push_back(num_c->GetXaxis()->GetBinLowEdge(xup)+num_c->GetXaxis()->GetBinWidth(xup));//upper border of last bin
    for(int i = ylow; i<=yup; ++i){
        yy.push_back(num_c->GetYaxis()->GetBinLowEdge(i));}
    yy.push_back(num_c->GetYaxis()->GetBinLowEdge(yup)+num_c->GetYaxis()->GetBinWidth(yup));//upper border of last bin
    int sizex = xx.size();
    int sizey = yy.size();
    double x[sizex], y[sizey];
    for(int i = 0; i<sizex; ++i){ x[i] = xx[i]; }
    for(int i = 0; i<sizey; ++i){ y[i] = yy[i]; }
    TH2D *num_sm = new TH2D("num_sm", "num_sm", sizex-1, x, sizey-1, y); num_sm->Sumw2();
    for(int nx=1; nx<=num_sm->GetNbinsX(); ++nx){
        for(int ny=1; ny<=num_sm->GetNbinsY(); ++ny){
            num_sm->SetBinContent(nx, ny, num_c->GetBinContent(num_c->FindBin(num_sm->GetXaxis()->GetBinCenter(nx), num_sm->GetYaxis()->GetBinCenter(ny) )));
            num_sm->SetBinError(nx, ny, num_c->GetBinError(num_c->FindBin(num_sm->GetXaxis()->GetBinCenter(nx), num_sm->GetYaxis()->GetBinCenter(ny) )));
        }
    }
    num_sm->RebinX( num_sm->GetNbinsX() );
    num_sm->RebinY( num_sm->GetNbinsY() );
    (eff) = (double)num_sm->GetBinContent(1,1);
    (err) = (double)num_sm->GetBinError(1,1);
    delete num_sm;
}

//load the MT2trees from the samples.dat
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






