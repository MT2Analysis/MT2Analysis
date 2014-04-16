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
#include "TGraph.h"
#include "TGraphErrors.h"
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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"

//run via root -l -b -q TauEstimationFromRootfile.C++

using namespace std;

void TauEstimationFromRootfile();
TH1D* RebinThisHistogram(TH1D* histogram, int nrebinnedbins, double* rebinnedbins);
TEfficiency* RebinThisEfficiency(TEfficiency* tefficiency, int nrebinnedbins, double* rebinnedbins);
void GetTauEstimate(map<string, TH1D*> histograms, map<string, TEfficiency*> tefficiencies, double rel_sys_uncert,  double rel_sys_uncert_bg, int version, Bool_t PrintSummaryTable, Bool_t PrintYieldTable, Bool_t PrintPredictionCard, Bool_t makeFullPrintout=false, Bool_t SavePrediction=true, Bool_t PlotPrediction=false);
void PredictionCard(vector<int> sr, vector<int> htr, vector<int> MT2bin, vector<double> mctruth, vector<double> mctrutherr, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, vector<double> LLeff, vector<double> LLefferr, vector<double> RLL, vector<double> RLLerr, vector<double> MTeff, vector<double> MTefferr);
void PredictionCardSplitted(vector<int> sr, vector<int> htr, vector<int> MT2bin, vector<double> mctruth, vector<double> mctrutherr, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, vector<double> LLeff, vector<double> LLefferr, vector<double> RLL, vector<double> RLLerr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruthW, vector<double> mctruthTT, vector<double> numW, vector<double> numTT, vector<double> numMC);//neew
void YieldTable(vector<int> sr, vector<int> htr, vector<double> MT2low, vector<double> MT2up, vector<double> numQCD, vector<double> numZ, vector<double> numW, vector<double> numT, vector<double> numOther, vector<double> numMC, vector<double> BGerr, vector<double> numData, vector<double> numBGLL);
void SummaryTable(int version, vector<int> sr, vector<int> htr, vector<double> MT2low, vector<double> MT2up, vector<double> numtrueW, vector<double> numtrueW_bg, vector<double> numData, vector<double> numBG, vector<double> numBGLL, vector<double> LLeff, vector<double> LLefferr, vector<double> RLL, vector<double> RLLerr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruth, vector<double> mctrutherr, vector<double> MCpred, vector<double> MCpred_stat_err, vector<double> MCpred_syst_err, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err);
void PredictionFile(int version, vector<int> sr, vector<int> htr, vector<double> MT2low, vector<double> MT2up, vector<double> numtrueW, vector<double> numtrueW_bg, vector<double> numData, vector<double> numBG, vector<double> numBGLL, vector<double> LLeff, vector<double> LLefferr, vector<double> RLL, vector<double> RLLerr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruth, vector<double> mctrutherr, vector<double> MCpred, vector<double> MCpred_stat_err, vector<double> MCpred_syst_err, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, Bool_t PlotPrediction);


const int fVerbose = 3; // already defined below
TString fPath;


Bool_t fMET        = false;//use MET triggers (i.e. used for low HT region)
Bool_t fHT         = true; //use HT triggers (i.e. used for medium+high HT region), will be set false automatically if fMET==true
Bool_t fRebin      = true;//do not bin-by-bin estimation but one estimate per HT/topological region (this means not along MT2)

Bool_t fISRreweight=true; //reweight MC according to SUSY's 'ISR recipe' - influence on estimate is minimal, but MC truth changes
                          //note, that if this flag is false you need to redo the TEfficiencies (fdoTEff=true)

Bool_t fresave     = false;//keep that false - you don't need it
TString foldfile   = "TauEstimation/HT/NoISR_TauEstimationHistograms.root";
//all the below should be done, unless it was some emergency saved file
Bool_t foverflow   = false;//add overflow or not, default = false (as this should have been done before)
Bool_t faddMC      = false;//add all MC samples to mc, default = false (as this should have been done before)
Bool_t fdoTEff     = false;//compute TEfficiencies
TString  outputdir = "Filtered/TauEstimation/";

Bool_t fbTagReweight= true;//reweight MC according to BTV SF weights - default is true
Bool_t fbTagError   = true;//compute additionally error due to BTV SF uncertainty - default is true

std::ostringstream* fLogStream     = 0;
Bool_t  fWriteToFile               = false; // writes couts to file
Bool_t  fAppend                    = true; // append at end of file (if existing), otherwise delete previous content

Bool_t fData       = true;//set false if don't want to run over data
Bool_t fQCD        = true;//set false if don't want to run over QCD samples
Bool_t fSusy       = false;//set false if don't want to run over Susy samples

Bool_t fbTagReweight = true;//reweight MC according to BTV SF weights - default is true
Bool_t fbTagError    = true;//compute additionally error due to BTV SF uncertainty - default is true

Bool_t fVetoLL       = true;//is true if vetoing lost electron/muon by geninfo (only in final printout), see fLLError, default is true
Bool_t fPlotLLTab    = true;//specify LL components in outcometables - is false if fVetoLL==false
Bool_t fUseMTcut     = true;//require MT cut for Taus to reduce signal contamination, default is true

Bool_t fUseDoubleTau = false;	// event called good if there are two generated taus (i.e. not only single tau events) - this is optional and has no too big influence
Bool_t fDoubleTau_BG = false;	// events that are double lost taus (two lost taus) are background, if false (default) they are part of yield you want to predict

Double_t frel_sys_uncert    = 0.05;//uncertainty on lepton efficiency (reconstruction and acceptance), default is 5 percent
Double_t frel_sys_uncert_bg =  0.2;//uncertainty due to background subtraction, default is 50 percent
Double_t fdoubeLLerr        =  1.0;//uncertainty due to double lost leptons, default is 100 percent - use at own risk (there is no good recipe for that)
Double_t frelMTerr          = 0.05;//uncertainty due to MT cut on lepton-MET system, default is 5 percent
Double_t fLLError           =  0.5;//if you have lost electron/muon and lost tau, by default we count it as lost electron/muon, thus we veto it here using MC truth, make this error depending on LL estimate, default is 50 percent

Bool_t fTopEfficencies      = true; //lepton efficiency from top sample only, if false only from W sample; works only if fWeightedProb=false; used for debugging default = true;
Bool_t  fWeightedProb       = true; // this should be true at all times. Take LostLepton efficiency from both W+Top sample
Bool_t  fIncludeTop         = true; // this should be true at all times. LostLeption yield estimated for Top+W sample; if false take only W sample
Bool_t  fTopOnly            = false;// this should be false at all times. LostLeption yield estimated for Top sample only
Bool_t  fIncludeSingleTop   = true; // this should be true at all times. Top sample contains also single top, not only ttbar.

const int sampletypesize = 12;
string sample_type[sampletypesize] = {"QCD", "WJets", "ZJets", "TTbar", "SingleTop", "Top", "WandTop", "noWandTop", "Other", "mc", "susy", "data"};
const int signalregionsize = 9;
string signal_region[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b"};
const int HTbinsize = 2;
string HT_bin[HTbinsize] = {"HTge450", "HTge750"};//dummy

//This function does exactly the same as TauEstimation.C, but it uses the histograms you created before with TauEstimation.C
//i.e. you do not run over the MT2trees anymore.
//this is very useful if you want to do the estimate with flags set differently (and that is also the reason why I kept all possible histograms)
//This code has no explanations as all relevant information can be found in TauEstimation.C
void TauEstimationFromRootfile(){

	fLogStream = new std::ostringstream();

	if(fMET){ HT_bin[0] = "HTge450"; }
	if(fHT) { HT_bin[0] = "HTge750"; HT_bin[1] = "HTge1200"; }
	if(fMET==true) fHT = false;
	if(fMET==false && fHT==false) fHT = true;
	if(fbTagReweight==false) fbTagError = false;
	if(fIncludeTop==false) fIncludeSingleTop = false;
	if(fIncludeTop==false) fTopOnly = false;
  	gROOT->ProcessLine(".x SetStyle_PRD.C");

	if(fMET) outputdir = outputdir + "MET/";
	if(fHT)  outputdir = outputdir + "HT/";
    	Util::MakeOutputDir(outputdir);

	TString outputname = "TauEstimationHistogramsNew.root";

	map<string, TH1D*>    histos;
	map<string, TEfficiency*>    teff;
	vector<string> histonames; histonames.clear();
	histonames.push_back("TauEvents");			//all   1 reco tau events incl. MT cut
	histonames.push_back("TauEvents_noMT");			//all   1 reco tau events w/o   MT cut
	histonames.push_back("TauEvents_LL");			//all   1 reco tau events incl. MT cut, containing LostLepton (via genlept)
	histonames.push_back("TauEvents_noMT_LL");		//all   1 reco tau events w/o   MT cut, containing LostLepton (via genlept)
	histonames.push_back("TauEventsBG");			//all   1 reco tau events where there is no genhadtau from W/Top decay
	histonames.push_back("TauEventsBG_LL");			//all   1 reco tau events where there is no genhadtau from W/Top decay, but genlep
	histonames.push_back("TauEvents_dTau");			//all   1 reco tau events where there is 2 genhadtau from W/Top decay
	histonames.push_back("TauEventsBG_noMT");		//all   1 reco tau events where there is no genhadtau from W/Top decay
	histonames.push_back("TauEventsBG_noMT_LL");		//all   1 reco tau events where there is no genhadtau from W/Top decay, but genlep
	histonames.push_back("TauEvents_noMT_dTau");		//all   1 reco tau events where there is 2 genhadtau from W/Top decay
	histonames.push_back("WEvents");			//genhadtau from W(Top) Decay - all
	histonames.push_back("WEvents_dTau");			//genhadtau from W(Top) Decay - all, double tau
	histonames.push_back("WEvents_noLL");			//genhadtau from W(Top) Decay - all, vetoing LostLepton (via genlept)
	histonames.push_back("WEventsAcc");			//genhadtau from W(Top) Decay - within acceptance
	histonames.push_back("WEventsAcc_noLL");		//genhadtau from W(Top) Decay - within acceptance, vetoing LostLepton (via genlept)
	histonames.push_back("WEventsReco");			//genhadtau from W(Top) Decay - recoed incl. MT cut
	histonames.push_back("WEventsReco_dTau");		//genhadtau from W(Top) Decay - recoed incl. MT cut
	histonames.push_back("WEventsReco_noLL");		//genhadtau from W(Top) Decay - recoed incl. MT cut, vetoing LostLepton (via genlept)
	histonames.push_back("WEventsReco_noMT");		//genhadtau from W(Top) Decay - recoed w/o   MT cut
	histonames.push_back("WEventsReco_noMT_noLL");		//genhadtau from W(Top) Decay - recoed w/o   MT cut, vetoing LostLepton (via genlept)
	histonames.push_back("NoRecoButGenTauEvents");		//events with no reco taus but gen had taus
	histonames.push_back("NoRecoButGenTauEvents_LL");	//events with no reco taus but gen had taus, vetoing LostLepton (via genlept)
	histonames.push_back("NoRecoButGenTauEvents_dTau");	//events with no reco taus but 2 gen had taus

	histonames.push_back("TrueRecoTauEvents");		//all reco taus w/o MT cut matched to gentaus//not this is not on event basis
	histonames.push_back("FakeRecoTauEvents");		//all reco taus w/o MT cut not matched to gentaus//not this is not on eventbasis
	histonames.push_back("AllRecoTauEvents");		//all reco taus w/o MT cut//not this is not on eventbasis
	histonames.push_back("AllRecoJetEvents");		//all reco jets w/o MT cut//not this is not on eventbasis
	histonames.push_back("AllGenTauAccEvents");		//all gen had taus within acceptance//not this is not on eventbasis

	histonames.push_back("TauEvents_SFup");
	histonames.push_back("TauEvents_noMT_SFup");
	histonames.push_back("TauEventsBG_SFup");
	histonames.push_back("TauEventsBG_LL_SFup");
	histonames.push_back("TauEvents_LL_SFup");
	histonames.push_back("TauEvents_dTau_SFup");
	histonames.push_back("TauEventsBG_noMT_SFup");
	histonames.push_back("TauEventsBG_noMT_LL_SFup");
	histonames.push_back("TauEvents_noMT_dTau_SFup");
	histonames.push_back("TauEvents_noMT_LL_SFup");
	histonames.push_back("NoRecoButGenTauEvents_SFup");
	histonames.push_back("NoRecoButGenTauEvents_LL_SFup");
	histonames.push_back("NoRecoButGenTauEvents_dTau_SFup");
	histonames.push_back("WEvents_SFup");
	histonames.push_back("WEvents_dTau_SFup");
	histonames.push_back("WEvents_noLL_SFup");

	histonames.push_back("TauEvents_SFdown");
	histonames.push_back("TauEvents_noMT_SFdown");
	histonames.push_back("TauEventsBG_SFdown");
	histonames.push_back("TauEventsBG_LL_SFdown");
	histonames.push_back("TauEvents_LL_SFdown");
	histonames.push_back("TauEvents_dTau_SFdown");
	histonames.push_back("TauEventsBG_noMT_SFdown");
	histonames.push_back("TauEventsBG_noMT_LL_SFdown");
	histonames.push_back("TauEvents_noMT_dTau_SFdown");
	histonames.push_back("TauEvents_noMT_LL_SFdown");
	histonames.push_back("NoRecoButGenTauEvents_SFdown");
	histonames.push_back("NoRecoButGenTauEvents_LL_SFdown");
	histonames.push_back("NoRecoButGenTauEvents_dTau_SFdown");
	histonames.push_back("WEvents_SFdown");
	histonames.push_back("WEvents_dTau_SFdown");
	histonames.push_back("WEvents_noLL_SFdown");


	vector<string> teffnames; teffnames.clear();


	TFile *oldfile = TFile::Open(foldfile);
//	bool av = true;//append vector --> maybe store histonames vector
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTbinsize;        ++i3){
		if(fMET && i3==1) continue;
		string hs = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_") + sample_type[i1];
		string mapname;
		for(unsigned int i0 = 0; i0<histonames.size(); ++i0){
			mapname = histonames[i0] + hs;
			string mapnamenosir = "noISR_" + histonames[i0] + hs;
			if(fISRreweight && histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)oldfile->Get(mapname.c_str());
			else if(           histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)oldfile->Get(mapnamenosir.c_str());
		}
	}}}
	//CHECK THAT
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Sumw2();}

	if(foverflow){
	//add overflow to last bin
	cout << "add overflow to last bin -> does not work for TEff" << endl;
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
	}

	if(faddMC){
	cout << "add all samples to mc, etc." << endl;
	//add TTbar and SingleTop to Top
	//add WJets and Top to WandTop
	//add rest to noWandTop
	//add noWandTop and WandTop to mc
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTbinsize;        ++i3){
		if(fMET && i3==1) continue;
		string hsq   = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_QCD");
		string hsw   = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_WJets");
		string hsz   = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_ZJets");
		string hstt  = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_TTbar");
		string hsst  = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_SingleTop");
		string hst   = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_Top");
		string hswt  = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_WandTop");
		string hsnwt = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_noWandTop");
		string hso   = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_Other");
		string hsmc  = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_mc");
		for(unsigned int i0 = 0; i0<histonames.size();++i0){
			histos[(histonames[i0]+hsnwt)]->Add(histos[(histonames[i0]+hsq)],  1);//noWandTop
			histos[(histonames[i0]+hsnwt)]->Add(histos[(histonames[i0]+hsz)],  1);
			histos[(histonames[i0]+hsnwt)]->Add(histos[(histonames[i0]+hso)],  1);
			histos[(histonames[i0]+hst)  ]->Add(histos[(histonames[i0]+hsst)], 1);//Top
			histos[(histonames[i0]+hst)  ]->Add(histos[(histonames[i0]+hstt)], 1);
			histos[(histonames[i0]+hswt) ]->Add(histos[(histonames[i0]+hst)],  1);//WandTop
			histos[(histonames[i0]+hswt) ]->Add(histos[(histonames[i0]+hsw)],  1);
			histos[(histonames[i0]+hsmc) ]->Add(histos[(histonames[i0]+hswt)], 1);//mc
			histos[(histonames[i0]+hsmc) ]->Add(histos[(histonames[i0]+hsnwt)],1);
		}
	}}
	}
	if(fdoTEff){
	cout << "finalize efficiencies" << endl;
	//now create efficiencies histograms and tefficiencies directly from histograms;
	teffnames.push_back("TEff_WEffRecoVsAll");
	teffnames.push_back("TEff_WEffRecoVsAcc");
	teffnames.push_back("TEff_WEffAccVsAll" );
	teffnames.push_back("TEff_WEffRecoVsAll_noMT");
	teffnames.push_back("TEff_WEffRecoVsAcc_noMT");
	teffnames.push_back("TEff_WEffRecoVsAll_noLL");
	teffnames.push_back("TEff_WEffRecoVsAcc_noLL");
	teffnames.push_back("TEff_WEffAccVsAll_noLL" );
	teffnames.push_back("TEff_WEffRecoVsAll_noMT_noLL");
	teffnames.push_back("TEff_WEffRecoVsAcc_noMT_noLL");
	teffnames.push_back("TEff_MTEff_onlygoodevts");
	teffnames.push_back("TEff_MTEff_onlygoodevts_noLL");
	teffnames.push_back("TEff_MTEff");
	teffnames.push_back("TEff_TauPurity");//XXXX
	teffnames.push_back("TEff_TauEfficiency");//XXXX
	teffnames.push_back("TEff_TauJetFakeRate");
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTbinsize;        ++i3){
		if(fMET && i3==1) continue;
		string hs   = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_") + sample_type[i1];
		string mapname;
		TH1D *passRecoVsAll          = (TH1D*)histos["WEventsReco"          +hs]->Clone("passRecoVsAll");
		TH1D  *totRecoVsAll          = (TH1D*)histos["WEvents"              +hs]->Clone( "totRecoVsAll");
		TH1D *passRecoVsAcc          = (TH1D*)histos["WEventsReco"          +hs]->Clone("passRecoVsAcc");
		TH1D  *totRecoVsAcc          = (TH1D*)histos["WEventsAcc"           +hs]->Clone( "totRecoVsAcc");
		TH1D *passRecoVsAll_noMT     = (TH1D*)histos["WEventsReco_noMT"     +hs]->Clone("passRecoVsAll_noMT");
		TH1D  *totRecoVsAll_noMT     = (TH1D*)histos["WEvents"              +hs]->Clone( "totRecoVsAll_noMT");
		TH1D *passRecoVsAcc_noMT     = (TH1D*)histos["WEventsReco_noMT"     +hs]->Clone("passRecoVsAcc_noMT");
		TH1D  *totRecoVsAcc_noMT     = (TH1D*)histos["WEventsAcc"           +hs]->Clone( "totRecoVsAcc_noMT");
		TH1D *passAccVsAll           = (TH1D*)histos["WEventsAcc"           +hs]->Clone("passAccVsAll" );
		TH1D  *totAccVsAll           = (TH1D*)histos["WEvents"              +hs]->Clone( "totAccVsAll" );
		TH1D *passRecoVsAll_noLL     = (TH1D*)histos["WEventsReco_noLL"     +hs]->Clone("passRecoVsAll_noLL");
		TH1D  *totRecoVsAll_noLL     = (TH1D*)histos["WEvents_noLL"         +hs]->Clone( "totRecoVsAll_noLL");
		TH1D *passRecoVsAcc_noLL     = (TH1D*)histos["WEventsReco_noLL"     +hs]->Clone("passRecoVsAcc_noLL");
		TH1D  *totRecoVsAcc_noLL     = (TH1D*)histos["WEventsAcc_noLL"      +hs]->Clone( "totRecoVsAcc_noLL");
		TH1D *passRecoVsAll_noMT_noLL= (TH1D*)histos["WEventsReco_noMT_noLL"+hs]->Clone("passRecoVsAll_noMT_noLL");
		TH1D  *totRecoVsAll_noMT_noLL= (TH1D*)histos["WEvents_noLL"         +hs]->Clone( "totRecoVsAll_noMT_noLL");
		TH1D *passRecoVsAcc_noMT_noLL= (TH1D*)histos["WEventsReco_noMT_noLL"+hs]->Clone("passRecoVsAcc_noMT_noLL");
		TH1D  *totRecoVsAcc_noMT_noLL= (TH1D*)histos["WEventsAcc_noLL"      +hs]->Clone( "totRecoVsAcc_noMT_noLL");
		TH1D *passAccVsAll_noLL      = (TH1D*)histos["WEventsAcc_noLL"      +hs]->Clone("passAccVsAll_noLL" );
		TH1D  *totAccVsAll_noLL      = (TH1D*)histos["WEvents_noLL"         +hs]->Clone( "totAccVsAll_noLL" );
		TH1D *passMT_2               = (TH1D*)histos["WEventsReco"          +hs]->Clone("passMT_2");
		TH1D  *totMT_2               = (TH1D*)histos["WEventsReco_noMT"     +hs]->Clone( "totMT_2");
		TH1D *passMT_2_noLL          = (TH1D*)histos["WEventsReco_noLL"     +hs]->Clone("passMT_2_noLL");
		TH1D  *totMT_2_noLL          = (TH1D*)histos["WEventsReco_noMT_noLL"+hs]->Clone( "totMT_2_noLL");
		TH1D *passMT                 = (TH1D*)histos["TauEvents"            +hs]->Clone("passMT");
		TH1D  *totMT                 = (TH1D*)histos["TauEvents_noMT"       +hs]->Clone( "totMT");
		TH1D *passtrue               = (TH1D*)histos["TrueRecoTauEvents"    +hs]->Clone("passtrue");
		TH1D  *tottrue               = (TH1D*)histos["AllRecoTauEvents"     +hs]->Clone( "tottrue");
		TH1D *passeff                = (TH1D*)histos["TrueRecoTauEvents"    +hs]->Clone("passeff");
		TH1D  *toteff                = (TH1D*)histos["AllGenTauAccEvents"   +hs]->Clone( "toteff");
		TH1D *passfake               = (TH1D*)histos["FakeRecoTauEvents"    +hs]->Clone("passfake");
		TH1D  *totfake               = (TH1D*)histos["AllRecoJetEvents"     +hs]->Clone( "totfake");
		mapname = "TEff_WEffRecoVsAll";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passRecoVsAll), (*totRecoVsAll));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_WEffRecoVsAcc";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passRecoVsAcc), (*totRecoVsAcc));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_WEffAccVsAll" ;
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passAccVsAll), (*totAccVsAll));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_WEffRecoVsAll_noMT";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passRecoVsAll_noMT), (*totRecoVsAll_noMT));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_WEffRecoVsAcc_noMT";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passRecoVsAcc_noMT), (*totRecoVsAcc_noMT));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_WEffRecoVsAll_noLL";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passRecoVsAll_noLL), (*totRecoVsAll_noLL));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_WEffRecoVsAcc_noLL";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passRecoVsAcc_noLL), (*totRecoVsAcc_noLL));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_WEffAccVsAll_noLL" ;
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passAccVsAll_noLL), (*totAccVsAll_noLL));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_WEffRecoVsAll_noMT_noLL";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passRecoVsAll_noMT_noLL), (*totRecoVsAll_noMT_noLL));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_WEffRecoVsAcc_noMT_noLL";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passRecoVsAcc_noMT_noLL), (*totRecoVsAcc_noMT_noLL));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_MTEff_onlygoodevts";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passMT_2), (*totMT_2));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_MTEff_onlygoodevts_noLL";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passMT_2_noLL), (*totMT_2_noLL));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_MTEff";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passMT), (*totMT));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_TauPurity";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passtrue), (*tottrue));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_TauEfficiency";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passeff), (*toteff));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
		mapname = "TEff_TauJetFakeRate";
		if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = new TEfficiency((*passfake), (*totfake));
		if(teff.count(mapname+hs) == 1 ) teff[mapname+hs]->SetNameTitle((mapname+hs).c_str(), (mapname+hs).c_str());
	}}}
	}
	else {
	//load efficiencies instead or redoing them!
		teffnames.push_back("TEff_WEffRecoVsAll");
		teffnames.push_back("TEff_WEffRecoVsAcc");
		teffnames.push_back("TEff_WEffAccVsAll" );
		teffnames.push_back("TEff_WEffRecoVsAll_noMT");
		teffnames.push_back("TEff_WEffRecoVsAcc_noMT");
		teffnames.push_back("TEff_WEffRecoVsAll_noLL");
		teffnames.push_back("TEff_WEffRecoVsAcc_noLL");
		teffnames.push_back("TEff_WEffAccVsAll_noLL" );
		teffnames.push_back("TEff_WEffRecoVsAll_noMT_noLL");
		teffnames.push_back("TEff_WEffRecoVsAcc_noMT_noLL");
		teffnames.push_back("TEff_MTEff_onlygoodevts");
		teffnames.push_back("TEff_MTEff_onlygoodevts_noLL");
		teffnames.push_back("TEff_MTEff");
		teffnames.push_back("TEff_TauPurity");//XXXX
		teffnames.push_back("TEff_TauEfficiency");//XXXX
		teffnames.push_back("TEff_TauJetFakeRate");
		for(int i1 = 0; i1<sampletypesize;   ++i1){
		for(int i2 = 0; i2<signalregionsize; ++i2){
		for(int i3 = 0; i3<HTbinsize;        ++i3){
			if(fMET && i3==1) continue;
			string hs   = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_") + sample_type[i1];
			string mapname;
			for(unsigned int i0 = 0; i0<teffnames.size(); ++i0){
				mapname = teffnames[i0];
				if(teff.count(mapname+hs) == 0 ) teff[mapname+hs] = (TEfficiency*)oldfile->Get(mapname.c_str());
			}
		}}}
	}
	if(fresave){
		cout << "Saving." << endl;
		TFile *fsavefile = new TFile(outputdir + outputname,"RECREATE");
		fsavefile->cd();
		for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
			h->second->Write();
		}
		for(map<string,TEfficiency*>::iterator h=teff.begin(); h!=teff.end();++h){
			h->second->Write();
		}
		fsavefile->Close();
		cout << "Saved histograms in " << outputdir << outputname << endl;
	}


	//do estimation
	cout << "do estimation" << endl;
//versions for SummaryTable
//version == 0,1: with MCPred
//version == 2,3: without MCPred
//version == 0,2: with R_LL (i.e. 1-e / e)
//version == 1,3: with e
	if(fRebin){
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTbinsize;        ++i3){
		if(fMET && i3==1) continue;
		string hs   = string("_") + HT_bin[i3] + string("_") + signal_region[i2] + string("_") + sample_type[i1];
		string mapname;
		double startbin = histos["TauEvents"+hs]->GetBinLowEdge(1);
		double endbin = histos["TauEvents"+hs]->GetBinLowEdge(histos["TauEvents"+hs]->GetNbinsX())+histos["TauEvents"+hs]->GetBinWidth(histos["TauEvents"+hs]->GetNbinsX());
		double rebin[2] = {startbin, endbin};
		int nrebin = 1;
		histos["TauEvents"                        +hs] = RebinThisHistogram(histos["TauEvents"                        +hs], nrebin, rebin);
		histos["TauEvents_noMT"                   +hs] = RebinThisHistogram(histos["TauEvents_noMT"                   +hs], nrebin, rebin);
		histos["TauEvents_LL"                     +hs] = RebinThisHistogram(histos["TauEvents_LL"                     +hs], nrebin, rebin);
		histos["TauEvents_noMT_LL"                +hs] = RebinThisHistogram(histos["TauEvents_noMT_LL"                +hs], nrebin, rebin);
		histos["TauEventsBG"                      +hs] = RebinThisHistogram(histos["TauEventsBG"                      +hs], nrebin, rebin);
		histos["TauEventsBG_LL"                   +hs] = RebinThisHistogram(histos["TauEventsBG_LL"                   +hs], nrebin, rebin);
		histos["TauEvents_dTau"                   +hs] = RebinThisHistogram(histos["TauEvents_dTau"                   +hs], nrebin, rebin);
		histos["TauEventsBG_noMT"                 +hs] = RebinThisHistogram(histos["TauEventsBG_noMT"                 +hs], nrebin, rebin);
		histos["TauEventsBG_noMT_LL"              +hs] = RebinThisHistogram(histos["TauEventsBG_noMT_LL"              +hs], nrebin, rebin);
		histos["TauEvents_noMT_dTau"              +hs] = RebinThisHistogram(histos["TauEvents_noMT_dTau"              +hs], nrebin, rebin);
		histos["WEvents"                          +hs] = RebinThisHistogram(histos["WEvents"                          +hs], nrebin, rebin);
		histos["WEvents_noLL"                     +hs] = RebinThisHistogram(histos["WEvents_noLL"                     +hs], nrebin, rebin);
		histos["WEvents_dTau"                     +hs] = RebinThisHistogram(histos["WEvents_dTau"                     +hs], nrebin, rebin);
		histos["WEventsReco_dTau"                 +hs] = RebinThisHistogram(histos["WEventsReco_dTau"                 +hs], nrebin, rebin);
		histos["WEventsAcc"                       +hs] = RebinThisHistogram(histos["WEventsAcc"                       +hs], nrebin, rebin);
		histos["WEventsAcc_noLL"                  +hs] = RebinThisHistogram(histos["WEventsAcc_noLL"                  +hs], nrebin, rebin);
		histos["WEventsReco"                      +hs] = RebinThisHistogram(histos["WEventsReco"                      +hs], nrebin, rebin);
		histos["WEventsReco_noLL"                 +hs] = RebinThisHistogram(histos["WEventsReco_noLL"                 +hs], nrebin, rebin);
		histos["WEventsReco_noMT"                 +hs] = RebinThisHistogram(histos["WEventsReco_noMT"                 +hs], nrebin, rebin);
		histos["WEventsReco_noMT_noLL"            +hs] = RebinThisHistogram(histos["WEventsReco_noMT_noLL"            +hs], nrebin, rebin);
		histos["NoRecoButGenTauEvents"            +hs] = RebinThisHistogram(histos["NoRecoButGenTauEvents"            +hs], nrebin, rebin);
		histos["NoRecoButGenTauEvents_LL"         +hs] = RebinThisHistogram(histos["NoRecoButGenTauEvents_LL"         +hs], nrebin, rebin);
		histos["NoRecoButGenTauEvents_dTau"       +hs] = RebinThisHistogram(histos["NoRecoButGenTauEvents_dTau"       +hs], nrebin, rebin);
		histos["TrueRecoTauEvents"                +hs] = RebinThisHistogram(histos["TrueRecoTauEvents"                +hs], nrebin, rebin);
		histos["FakeRecoTauEvents"                +hs] = RebinThisHistogram(histos["FakeRecoTauEvents"                +hs], nrebin, rebin);
		histos["AllRecoTauEvents"                 +hs] = RebinThisHistogram(histos["AllRecoTauEvents"                 +hs], nrebin, rebin);
		histos["AllRecoJetEvents"                 +hs] = RebinThisHistogram(histos["AllRecoJetEvents"                 +hs], nrebin, rebin);
		histos["AllGenTauAccEvents"               +hs] = RebinThisHistogram(histos["AllGenTauAccEvents"               +hs], nrebin, rebin);
		histos["TauEvents_SFup"                   +hs] = RebinThisHistogram(histos["TauEvents_SFup"                   +hs], nrebin, rebin);
		histos["TauEvents_noMT_SFup"              +hs] = RebinThisHistogram(histos["TauEvents_noMT_SFup"              +hs], nrebin, rebin);
		histos["TauEventsBG_SFup"                 +hs] = RebinThisHistogram(histos["TauEventsBG_SFup"                 +hs], nrebin, rebin);
		histos["TauEventsBG_LL_SFup"              +hs] = RebinThisHistogram(histos["TauEventsBG_LL_SFup"              +hs], nrebin, rebin);
		histos["TauEvents_LL_SFup"                +hs] = RebinThisHistogram(histos["TauEvents_LL_SFup"                +hs], nrebin, rebin);
		histos["TauEvents_dTau_SFup"              +hs] = RebinThisHistogram(histos["TauEvents_dTau_SFup"              +hs], nrebin, rebin);
		histos["TauEventsBG_noMT_SFup"            +hs] = RebinThisHistogram(histos["TauEventsBG_noMT_SFup"            +hs], nrebin, rebin);
		histos["TauEventsBG_noMT_LL_SFup"         +hs] = RebinThisHistogram(histos["TauEventsBG_noMT_LL_SFup"         +hs], nrebin, rebin);
		histos["TauEvents_noMT_LL_SFup"           +hs] = RebinThisHistogram(histos["TauEvents_noMT_LL_SFup"           +hs], nrebin, rebin);
		histos["TauEvents_noMT_dTau_SFup"         +hs] = RebinThisHistogram(histos["TauEvents_noMT_dTau_SFup"         +hs], nrebin, rebin);
		histos["NoRecoButGenTauEvents_SFup"       +hs] = RebinThisHistogram(histos["NoRecoButGenTauEvents_SFup"       +hs], nrebin, rebin);
		histos["NoRecoButGenTauEvents_LL_SFup"    +hs] = RebinThisHistogram(histos["NoRecoButGenTauEvents_LL_SFup"    +hs], nrebin, rebin);
		histos["NoRecoButGenTauEvents_dTau_SFup"  +hs] = RebinThisHistogram(histos["NoRecoButGenTauEvents_dTau_SFup"  +hs], nrebin, rebin);
		histos["WEvents_SFup"                     +hs] = RebinThisHistogram(histos["WEvents_SFup"                     +hs], nrebin, rebin);
		histos["WEvents_noLL_SFup"                +hs] = RebinThisHistogram(histos["WEvents_noLL_SFup"                +hs], nrebin, rebin);
		histos["WEvents_dTau_SFup"                +hs] = RebinThisHistogram(histos["WEvents_dTau_SFup"                +hs], nrebin, rebin);
		histos["TauEvents_SFdown"                 +hs] = RebinThisHistogram(histos["TauEvents_SFdown"                 +hs], nrebin, rebin);
		histos["TauEvents_noMT_SFdown"            +hs] = RebinThisHistogram(histos["TauEvents_noMT_SFdown"            +hs], nrebin, rebin);
		histos["TauEventsBG_SFdown"               +hs] = RebinThisHistogram(histos["TauEventsBG_SFdown"               +hs], nrebin, rebin);
		histos["TauEventsBG_LL_SFdown"            +hs] = RebinThisHistogram(histos["TauEventsBG_LL_SFdown"            +hs], nrebin, rebin);
		histos["TauEvents_dTau_SFdown"            +hs] = RebinThisHistogram(histos["TauEvents_dTau_SFdown"            +hs], nrebin, rebin);
		histos["TauEvents_LL_SFdown"              +hs] = RebinThisHistogram(histos["TauEvents_LL_SFdown"              +hs], nrebin, rebin);
		histos["TauEventsBG_noMT_SFdown"          +hs] = RebinThisHistogram(histos["TauEventsBG_noMT_SFdown"          +hs], nrebin, rebin);
		histos["TauEventsBG_noMT_LL_SFdown"       +hs] = RebinThisHistogram(histos["TauEventsBG_noMT_LL_SFdown"       +hs], nrebin, rebin);
		histos["TauEvents_noMT_dTau_SFdown"       +hs] = RebinThisHistogram(histos["TauEvents_noMT_dTau_SFdown"       +hs], nrebin, rebin);
		histos["TauEvents_noMT_LL_SFdown"         +hs] = RebinThisHistogram(histos["TauEvents_noMT_LL_SFdown"         +hs], nrebin, rebin);
		histos["NoRecoButGenTauEvents_SFdown"     +hs] = RebinThisHistogram(histos["NoRecoButGenTauEvents_SFdown"     +hs], nrebin, rebin);
		histos["NoRecoButGenTauEvents_LL_SFdown"  +hs] = RebinThisHistogram(histos["NoRecoButGenTauEvents_LL_SFdown"  +hs], nrebin, rebin);
		histos["NoRecoButGenTauEvents_dTau_SFdown"+hs] = RebinThisHistogram(histos["NoRecoButGenTauEvents_dTau_SFdown"+hs], nrebin, rebin);
		histos["WEvents_SFdown"                   +hs] = RebinThisHistogram(histos["WEvents_SFdown"                   +hs], nrebin, rebin);
		histos["WEvents_noLL_SFdown"              +hs] = RebinThisHistogram(histos["WEvents_noLL_SFdown"              +hs], nrebin, rebin);
		histos["WEvents_dTau_SFdown"              +hs] = RebinThisHistogram(histos["WEvents_dTau_SFdown"              +hs], nrebin, rebin);
		teff["TEff_WEffRecoVsAll"          +hs] = RebinThisEfficiency(teff["TEff_WEffRecoVsAll"          +hs], nrebin, rebin);
		teff["TEff_WEffRecoVsAcc"          +hs] = RebinThisEfficiency(teff["TEff_WEffRecoVsAcc"          +hs], nrebin, rebin);
		teff["TEff_WEffAccVsAll"           +hs] = RebinThisEfficiency(teff["TEff_WEffAccVsAll"           +hs], nrebin, rebin);
		teff["TEff_WEffRecoVsAll_noMT"     +hs] = RebinThisEfficiency(teff["TEff_WEffRecoVsAll_noMT"     +hs], nrebin, rebin);
		teff["TEff_WEffRecoVsAcc_noMT"     +hs] = RebinThisEfficiency(teff["TEff_WEffRecoVsAcc_noMT"     +hs], nrebin, rebin);
		teff["TEff_WEffRecoVsAll_noLL"     +hs] = RebinThisEfficiency(teff["TEff_WEffRecoVsAll_noLL"     +hs], nrebin, rebin);
		teff["TEff_WEffRecoVsAcc_noLL"     +hs] = RebinThisEfficiency(teff["TEff_WEffRecoVsAcc_noLL"     +hs], nrebin, rebin);
		teff["TEff_WEffAccVsAll_noLL"      +hs] = RebinThisEfficiency(teff["TEff_WEffAccVsAll_noLL"      +hs], nrebin, rebin);
		teff["TEff_WEffRecoVsAll_noMT_noLL"+hs] = RebinThisEfficiency(teff["TEff_WEffRecoVsAll_noMT_noLL"+hs], nrebin, rebin);
		teff["TEff_WEffRecoVsAcc_noMT_noLL"+hs] = RebinThisEfficiency(teff["TEff_WEffRecoVsAcc_noMT_noLL"+hs], nrebin, rebin);
		teff["TEff_MTEff_onlygoodevts"     +hs] = RebinThisEfficiency(teff["TEff_MTEff_onlygoodevts"     +hs], nrebin, rebin);
		teff["TEff_MTEff_onlygoodevts_noLL"+hs] = RebinThisEfficiency(teff["TEff_MTEff_onlygoodevts_noLL"+hs], nrebin, rebin);
		teff["TEff_MTEff"                  +hs] = RebinThisEfficiency(teff["TEff_MTEff"                  +hs], nrebin, rebin);
		teff["TEff_TauPurity"              +hs] = RebinThisEfficiency(teff["TEff_TauPurity"              +hs], nrebin, rebin);
		teff["TEff_TauEfficiency"          +hs] = RebinThisEfficiency(teff["TEff_TauEfficiency"          +hs], nrebin, rebin);
		teff["TEff_TauJetFakeRate"         +hs] = RebinThisEfficiency(teff["TEff_TauJetFakeRate"         +hs], nrebin, rebin);
	}}}
	}

//	GetTauEstimatehistograms,tefficiencies,rel_sys_err,rel_sys_err_bg,version,SummaryTable,YieldTable,PredictionCard,makeFullPrintout)
	GetTauEstimate(histos, teff, frel_sys_uncert, frel_sys_uncert_bg, 0, true, true, true, true, true, true);
	//do full printout
	//GetTauEstimate(histos, teff, frel_sys_uncert, frel_sys_uncert_bg, 0, true, true, true, true, true, true);


	if(fWriteToFile && fAppend){
		TString logname =outputdir + "results.log"; 
		ofstream f_log (logname.Data(), ios::app);
		f_log << fLogStream->str();
		cout << "wrote results into  " << logname << " (appended at the end of old file)" << endl;
	}else if(fWriteToFile){
		TString logname =outputdir + "results.log"; 
		ofstream f_log (logname.Data(), ios::trunc);
		f_log << fLogStream->str();
		cout << "wrote results into  " << logname <<  " (old file replaced)" << endl;
	} else{
		cout << fLogStream->str();
	}
	delete fLogStream;

}//void TauEstimationFromRootfile()

TH1D* RebinThisHistogram(TH1D* histogram, int nrebinnedbins, double* rebinnedbins){

	string name = histogram->GetName();
	histogram->SetName((name + "_original").c_str());
	TH1D* temphists = (TH1D*)histogram->Rebin(nrebinnedbins, (name + "_rebinned").c_str(), rebinnedbins);
	temphists->SetName((name).c_str());
	return temphists;
}

TEfficiency* RebinThisEfficiency(TEfficiency* tefficiency, int nrebinnedbins, double* rebinnedbins){

	string name = tefficiency->GetName();
	TH1D *pass = (TH1D*)tefficiency->GetCopyPassedHisto(); TH1D *total = (TH1D*)tefficiency->GetCopyTotalHisto();
	name = pass->GetName();
	TH1D *passrebinned = (TH1D*)pass->Rebin(nrebinnedbins, (name + "_rebinned").c_str(), rebinnedbins);
	name = total->GetName();
	TH1D *totalrebinned = (TH1D*)total->Rebin(nrebinnedbins, (name + "_rebinned").c_str(), rebinnedbins);
	name = tefficiency->GetName();
	tefficiency->SetNameTitle((name+"_original").c_str(), (name+"_original").c_str());
	TEfficiency *temp = new TEfficiency((*passrebinned), (*totalrebinned));
	temp->SetNameTitle((name+"_temp").c_str(), (name).c_str());
	TEfficiency *clone = (TEfficiency*)temp->Clone((name).c_str());
	delete temp;
	return clone;
}

//includeTop, includeSingleTop, TopOnly, OnlyMCClosure, MakeEfficienciesPablo, WeightedProb, TopEfficencies are outside
//versions for SummaryTable
//version == 0,1: with MCPred
//version == 2,3: without MCPred
//version == 0,2: with R_LL (i.e. 1-e / e)
//version == 1,3: with e
void GetTauEstimate(map<string, TH1D*> histograms, map<string, TEfficiency*> tefficiencies, double rel_sys_uncert,  double rel_sys_uncert_bg, int version, Bool_t PrintSummaryTable, Bool_t PrintYieldTable, Bool_t PrintPredictionCard, Bool_t makeFullPrintout, Bool_t SavePrediction, Bool_t PlotPrediction){

	map<string, TH1D*> hists = histograms;
	map<string, TEfficiency*> teffs = tefficiencies;

	double dummy;

	//printoutnumbers - store them as vector
	//printoutnumbers - store region
	vector<int> sr; sr.clear();//signal region
	vector<int> htr; htr.clear();//ht region
	vector<double> MT2low; MT2low.clear();//lower bound of bin
	vector<double> MT2up; MT2up.clear();//upper bound of bin // for final bin store 10000.
	vector<int> MT2bin; MT2bin.clear();//store only binnumber --> used for datacard
	//printoutnumbers - store numbers which might be used in prinouts
	vector<double> mctruth; mctruth.clear();
	vector<double> mctrutherr; mctrutherr.clear();
	vector<double> mctruthTT; mctruthTT.clear();//neew
	vector<double> mctruthW;  mctruthW.clear();//neew
	vector<double> datapred; datapred.clear();
	vector<double> datapred_stat_err; datapred_stat_err.clear();
	vector<double> datapred_syst_err; datapred_syst_err.clear();
	vector<double> MCpred; MCpred.clear();
	vector<double> MCpred_stat_err; MCpred_stat_err.clear();
	vector<double> MCpred_syst_err; MCpred_syst_err.clear();
	vector<double> LLeff; LLeff.clear();
	vector<double> LLefferr; LLefferr.clear();
	vector<double> RLL; RLL.clear();
	vector<double> RLLerr; RLLerr.clear();
	vector<double> MTeff; MTeff.clear();
	vector<double> MTefferr; MTefferr.clear();
	vector<double> numtrueW; numtrueW.clear();
	vector<double> numtrueW_bg; numtrueW_bg.clear();
	vector<double> numData; numData.clear();
	vector<double> numBG; numBG.clear();
	vector<double> numBGLL; numBGLL.clear();
	vector<double> numQCD; numQCD.clear();
	vector<double> numZ; numZ.clear();
	vector<double> numW; numW.clear();
	vector<double> numT; numT.clear();
	vector<double> numTT; numTT.clear();//neew
	vector<double> numOther; numOther.clear();
	vector<double> numMC; numMC.clear();
	vector<double> BGerr; BGerr.clear();

	//start printoutloop;
	for(int i3 = 0; i3<HTbinsize;        ++i3){
	for(int i2 = 0; i2<signalregionsize; ++i2){
		if(fMET && i3==1) continue;
	//as all histograms have the same binning just take first histogram 
	string hshl = string("_") + HT_bin[i3] + string("_") + signal_region[i2];
	for(int nx = 1; nx<=hists[string("TauEvents"+hshl+"_mc")]->GetNbinsX(); ++nx){
		//first push_back the region information
		sr.push_back(i2);
		if(fMET) htr.push_back(i3);
		else     htr.push_back(i3+1);
		MT2low.push_back(hists[string("TauEvents"+hshl+"_mc")]->GetBinLowEdge(nx) );
		if(nx!=hists[string("TauEvents"+hshl+"_mc")]->GetNbinsX())
		   MT2up.push_back( hists[string("TauEvents"+hshl+"_mc")]->GetBinLowEdge(nx) + hists[string("TauEvents"+hshl+"_mc")]->GetBinWidth(nx) );
		else
		   MT2up.push_back(10000.);
		MT2bin.push_back(nx);

		//efficiencies
		double W_prob(0.),   W_prob_err(0.);     double W_prob_noLL(0.),   W_prob_err_noLL(0.);
		double W_acc(0.),    W_acc_err(0.);      double W_acc_noLL(0.),    W_acc_err_noLL(0.);
		double W_rec(0.),    W_rec_err(0.);      double W_rec_noLL(0.),    W_rec_err_noLL(0.);
		double Top_prob(0.), Top_prob_err(0.);   double Top_prob_noLL(0.), Top_prob_err_noLL(0.);
		double Top_acc(0.),  Top_acc_err(0.);    double Top_acc_noLL(0.),  Top_acc_err_noLL(0.);
		double Top_rec(0.),  Top_rec_err(0.);    double Top_rec_noLL(0.),  Top_rec_err_noLL(0.);
		double WT_prob(0.),  WT_prob_err(0.);    double WT_prob_noLL(0.),  WT_prob_err_noLL(0.);
		double WT_acc(0.),   WT_acc_err(0.);     double WT_acc_noLL(0.),   WT_acc_err_noLL(0.);
		double WT_rec(0.),   WT_rec_err(0.);     double WT_rec_noLL(0.),   WT_rec_err_noLL(0.);
		//the below should be the standard - no matter if using MT or not, acc is from above as there is no reco
		double W_prob_nomt(0.),   W_prob_err_nomt(0.);      double W_prob_nomt_noLL(0.),   W_prob_err_nomt_noLL(0.);
		double W_rec_nomt(0.),    W_rec_err_nomt(0.);       double W_rec_nomt_noLL(0.),    W_rec_err_nomt_noLL(0.);
		double Top_prob_nomt(0.), Top_prob_err_nomt(0.);    double Top_prob_nomt_noLL(0.), Top_prob_err_nomt_noLL(0.);
		double Top_rec_nomt(0.),  Top_rec_err_nomt(0.);     double Top_rec_nomt_noLL(0.),  Top_rec_err_nomt_noLL(0.);
		double WT_prob_nomt(0.),  WT_prob_err_nomt(0.);     double WT_prob_nomt_noLL(0.),  WT_prob_err_nomt_noLL(0.);
		double WT_rec_nomt(0.),   WT_rec_err_nomt(0.);      double WT_rec_nomt_noLL(0.),   WT_rec_err_nomt_noLL(0.);

		double mteff(0.),       mtefferr(0.);
		double mteffdata(0.),   mtefferrdata(0.);
		double mteff2(0.),      mtefferr2(0.);//only good events
		double mteff2_noLL(0.), mtefferr2_noLL(0.);//only good events
		double taueff(0.),      tauefferr(0.);
		double taupur(0.),      taupurerr(0.);
		double taufake(0.),     taufakeerr(0.);
		//event yields
		double  nW(0.),  nW_bg(0.),  nW_dTau(0.),  nW_LL(0.),  nW_bg_LL(0.);//W
		double nTT(0.), nTT_bg(0.), nTT_dTau(0.), nTT_LL(0.), nTT_bg_LL(0.);//TTbar
		double nST(0.), nST_bg(0.), nST_dTau(0.), nST_LL(0.), nST_bg_LL(0.);//SingleTop
		double  nT(0.),  nT_bg(0.),  nT_dTau(0.),  nT_LL(0.),  nT_bg_LL(0.);//Top=SingleTop+TTbar
		double nWT(0.), nWT_bg(0.), nWT_dTau(0.), nWT_LL(0.), nWT_bg_LL(0.);//W+Top
		double  nW_lv(0.),  nW_lv_err(0.),  nW_lv_LL(0.),  nW_lv_LL_err(0.),  nW_lv_dTau(0.),  nW_lv_dTau_err(0.);
		double nTT_lv(0.), nTT_lv_err(0.), nTT_lv_LL(0.), nTT_lv_LL_err(0.), nTT_lv_dTau(0.), nTT_lv_dTau_err(0.);
		double nST_lv(0.), nST_lv_err(0.), nST_lv_LL(0.), nST_lv_LL_err(0.), nST_lv_dTau(0.), nST_lv_dTau_err(0.);
		double  nT_lv(0.),  nT_lv_err(0.),  nT_lv_LL(0.),  nT_lv_LL_err(0.),  nT_lv_dTau(0.),  nT_lv_dTau_err(0.);
		double nWT_lv(0.), nWT_lv_err(0.), nWT_lv_LL(0.), nWT_lv_LL_err(0.), nWT_lv_dTau(0.), nWT_lv_dTau_err(0.);
		double  nW_nomt(0.),  nW_bg_nomt(0.),  nW_dTau_nomt(0.),  nW_LL_nomt(0.),  nW_bg_LL_nomt(0.);
		double nTT_nomt(0.), nTT_bg_nomt(0.), nTT_dTau_nomt(0.), nTT_LL_nomt(0.), nTT_bg_LL_nomt(0.);
		double nST_nomt(0.), nST_bg_nomt(0.), nST_dTau_nomt(0.), nST_LL_nomt(0.), nST_bg_LL_nomt(0.);
		double  nT_nomt(0.),  nT_bg_nomt(0.),  nT_dTau_nomt(0.),  nT_LL_nomt(0.),  nT_bg_LL_nomt(0.);
		double nWT_nomt(0.), nWT_bg_nomt(0.), nWT_dTau_nomt(0.), nWT_LL_nomt(0.), nWT_bg_LL_nomt(0.);
		//for BTV SF uncertainty
		double  nW_bg_SFup(0.),  nW_dTau_SFup(0.),  nW_LL_SFup(0.),  nW_bg_LL_SFup(0.);
		double nTT_bg_SFup(0.), nTT_dTau_SFup(0.), nTT_LL_SFup(0.), nTT_bg_LL_SFup(0.);
		double nST_bg_SFup(0.), nST_dTau_SFup(0.), nST_LL_SFup(0.), nST_bg_LL_SFup(0.);
		double  nT_bg_SFup(0.),  nT_dTau_SFup(0.),  nT_LL_SFup(0.),  nT_bg_LL_SFup(0.);
		double nWT_bg_SFup(0.), nWT_dTau_SFup(0.), nWT_LL_SFup(0.), nWT_bg_LL_SFup(0.);
		double  nW_lv_SFup(0.),  nW_lv_LL_SFup(0.),  nW_lv_dTau_SFup(0.);
		double nTT_lv_SFup(0.), nTT_lv_LL_SFup(0.), nTT_lv_dTau_SFup(0.);
		double nST_lv_SFup(0.), nST_lv_LL_SFup(0.), nST_lv_dTau_SFup(0.);
		double  nT_lv_SFup(0.),  nT_lv_LL_SFup(0.),  nT_lv_dTau_SFup(0.);
		double nWT_lv_SFup(0.), nWT_lv_LL_SFup(0.), nWT_lv_dTau_SFup(0.);
		double  nW_bg_nomt_SFup(0.),  nW_dTau_nomt_SFup(0.),  nW_LL_nomt_SFup(0.),  nW_bg_LL_nomt_SFup(0.);
		double nTT_bg_nomt_SFup(0.), nTT_dTau_nomt_SFup(0.), nTT_LL_nomt_SFup(0.), nTT_bg_LL_nomt_SFup(0.);
		double nST_bg_nomt_SFup(0.), nST_dTau_nomt_SFup(0.), nST_LL_nomt_SFup(0.), nST_bg_LL_nomt_SFup(0.);
		double  nT_bg_nomt_SFup(0.),  nT_dTau_nomt_SFup(0.),  nT_LL_nomt_SFup(0.),  nT_bg_LL_nomt_SFup(0.);
		double nWT_bg_nomt_SFup(0.), nWT_dTau_nomt_SFup(0.), nWT_LL_nomt_SFup(0.), nWT_bg_LL_nomt_SFup(0.);
		double  nW_bg_SFdown(0.),  nW_dTau_SFdown(0.),  nW_LL_SFdown(0.),  nW_bg_LL_SFdown(0.);
		double nTT_bg_SFdown(0.), nTT_dTau_SFdown(0.), nTT_LL_SFdown(0.), nTT_bg_LL_SFdown(0.);
		double nST_bg_SFdown(0.), nST_dTau_SFdown(0.), nST_LL_SFdown(0.), nST_bg_LL_SFdown(0.);
		double  nT_bg_SFdown(0.),  nT_dTau_SFdown(0.),  nT_LL_SFdown(0.),  nT_bg_LL_SFdown(0.);
		double nWT_bg_SFdown(0.), nWT_dTau_SFdown(0.), nWT_LL_SFdown(0.), nWT_bg_LL_SFdown(0.);
		double  nW_lv_SFdown(0.),  nW_lv_LL_SFdown(0.),  nW_lv_dTau_SFdown(0.);
		double nTT_lv_SFdown(0.), nTT_lv_LL_SFdown(0.), nTT_lv_dTau_SFdown(0.);
		double nST_lv_SFdown(0.), nST_lv_LL_SFdown(0.), nST_lv_dTau_SFdown(0.);
		double  nT_lv_SFdown(0.),  nT_lv_LL_SFdown(0.),  nT_lv_dTau_SFdown(0.);
		double nWT_lv_SFdown(0.), nWT_lv_LL_SFdown(0.), nWT_lv_dTau_SFdown(0.);
		double  nW_bg_nomt_SFdown(0.),  nW_dTau_nomt_SFdown(0.),  nW_LL_nomt_SFdown(0.),  nW_bg_LL_nomt_SFdown(0.);
		double nTT_bg_nomt_SFdown(0.), nTT_dTau_nomt_SFdown(0.), nTT_LL_nomt_SFdown(0.), nTT_bg_LL_nomt_SFdown(0.);
		double nST_bg_nomt_SFdown(0.), nST_dTau_nomt_SFdown(0.), nST_LL_nomt_SFdown(0.), nST_bg_LL_nomt_SFdown(0.);
		double  nT_bg_nomt_SFdown(0.),  nT_dTau_nomt_SFdown(0.),  nT_LL_nomt_SFdown(0.),  nT_bg_LL_nomt_SFdown(0.);
		double nWT_bg_nomt_SFdown(0.), nWT_dTau_nomt_SFdown(0.), nWT_LL_nomt_SFdown(0.), nWT_bg_LL_nomt_SFdown(0.);
		//backgrounds - this does not contain any LL or dTau by definition
		double QCD_bg(0.), Z_bg(0.), Other_bg(0.), TT_bg(0.), ST_bg(0.), T_bg(0.), W_bg(0.), WT_bg(0.), nonWT_bg(0.);
		double QCD_bg_nomt(0.), Z_bg_nomt(0.), Other_bg_nomt(0.), TT_bg_nomt(0.), ST_bg_nomt(0.), T_bg_nomt(0.), W_bg_nomt(0.), WT_bg_nomt(0.), nonWT_bg_nomt(0.);
		//backgrounds - SF weights
		double QCD_bg_SFup(0.),   Z_bg_SFup(0.),   Other_bg_SFup(0.),   T_bg_SFup(0.),   TT_bg_SFup(0.),   ST_bg_SFup(0.),   W_bg_SFup(0.),   WT_bg_SFup(0.),   nonWT_bg_SFup(0.);
		double QCD_bg_SFdown(0.), Z_bg_SFdown(0.), Other_bg_SFdown(0.), T_bg_SFdown(0.), TT_bg_SFdown(0.), ST_bg_SFdown(0.), W_bg_SFdown(0.), WT_bg_SFdown(0.), nonWT_bg_SFdown(0.);
		double QCD_bg_nomt_SFup(0.),   Z_bg_nomt_SFup(0.),   Other_bg_nomt_SFup(0.),   T_bg_nomt_SFup(0.),   TT_bg_nomt_SFup(0.),   ST_bg_nomt_SFup(0.),   W_bg_nomt_SFup(0.),   WT_bg_nomt_SFup(0.),   nonWT_bg_nomt_SFup(0.);
		double QCD_bg_nomt_SFdown(0.), Z_bg_nomt_SFdown(0.), Other_bg_nomt_SFdown(0.), T_bg_nomt_SFdown(0.), TT_bg_nomt_SFdown(0.), ST_bg_nomt_SFdown(0.), W_bg_nomt_SFdown(0.), WT_bg_nomt_SFdown(0.), nonWT_bg_nomt_SFdown(0.);
		//data yields, mc yields and tau_lv
		double nData(0.), nData_nomt(0.);
		double MC_bg(0.), MC_bg_err(0.), MC_bg_nomt(0.), MC_bg_err_nomt(0.);

		taupur        = teffs["TEff_TauPurity"     +hshl+"_mc"  ]->GetEfficiency(nx);
		taupurerr     = teffs["TEff_TauPurity"     +hshl+"_mc"  ]->GetEfficiencyErrorLow(nx);
		dummy         = teffs["TEff_TauPurity"     +hshl+"_mc"  ]->GetEfficiencyErrorUp(nx); if(dummy>taupurerr   ) taupurerr    = dummy;
		taueff        = teffs["TEff_TauEfficiency" +hshl+"_mc"  ]->GetEfficiency(nx);
		tauefferr     = teffs["TEff_TauEfficiency" +hshl+"_mc"  ]->GetEfficiencyErrorLow(nx);
		dummy         = teffs["TEff_TauEfficiency" +hshl+"_mc"  ]->GetEfficiencyErrorUp(nx); if(dummy>tauefferr   ) tauefferr    = dummy;
		taufake       = teffs["TEff_TauJetFakeRate"+hshl+"_mc"  ]->GetEfficiency(nx);
		taufakeerr    = teffs["TEff_TauJetFakeRate"+hshl+"_mc"  ]->GetEfficiencyErrorLow(nx);
		dummy         = teffs["TEff_TauJetFakeRate"+hshl+"_mc"  ]->GetEfficiencyErrorUp(nx); if(dummy>taufakeerr  ) taufakeerr   = dummy;
		mteff         = teffs["TEff_MTEff"         +hshl+"_mc"  ]->GetEfficiency(nx);
		mtefferr      = teffs["TEff_MTEff"         +hshl+"_mc"  ]->GetEfficiencyErrorLow(nx);
		dummy         = teffs["TEff_MTEff"         +hshl+"_mc"  ]->GetEfficiencyErrorUp(nx); if(dummy>mtefferr    ) mtefferr     = dummy;
		mteffdata     = teffs["TEff_MTEff"         +hshl+"_data"]->GetEfficiency(nx);
		mtefferrdata  = teffs["TEff_MTEff"         +hshl+"_data"]->GetEfficiencyErrorLow(nx);
		dummy         = teffs["TEff_MTEff"         +hshl+"_data"]->GetEfficiencyErrorUp(nx); if(dummy>mtefferrdata) mtefferrdata = dummy;
		mteff2        = teffs["TEff_MTEff_onlygoodevts"     +hshl+"_WandTop"]->GetEfficiency(nx);
		mtefferr2     = teffs["TEff_MTEff_onlygoodevts"     +hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy         = teffs["TEff_MTEff_onlygoodevts"     +hshl+"_WandTop"]->GetEfficiencyErrorUp(nx); if(dummy>mtefferr2     ) mtefferr2      = dummy;
		mteff2_noLL   = teffs["TEff_MTEff_onlygoodevts_noLL"+hshl+"_WandTop"]->GetEfficiency(nx);
		mtefferr2_noLL= teffs["TEff_MTEff_onlygoodevts_noLL"+hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy         = teffs["TEff_MTEff_onlygoodevts_noLL"+hshl+"_WandTop"]->GetEfficiencyErrorUp(nx); if(dummy>mtefferr2_noLL) mtefferr2_noLL = dummy;
		W_prob       = teffs["TEff_WEffRecoVsAll"+hshl+"_WJets"  ]->GetEfficiency(nx);
		W_prob_err   = teffs["TEff_WEffRecoVsAll"+hshl+"_WJets"  ]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_WEffRecoVsAll"+hshl+"_WJets"  ]->GetEfficiencyErrorUp(nx); if(dummy>W_prob_err  ) W_prob_err   = dummy;
		Top_prob     = teffs["TEff_WEffRecoVsAll"+hshl+"_Top"    ]->GetEfficiency(nx);		//contains SingleTop
		Top_prob_err = teffs["TEff_WEffRecoVsAll"+hshl+"_Top"    ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy        = teffs["TEff_WEffRecoVsAll"+hshl+"_Top"    ]->GetEfficiencyErrorUp(nx); if(dummy>Top_prob_err) Top_prob_err = dummy;
		W_acc        = teffs["TEff_WEffAccVsAll" +hshl+"_WJets"  ]->GetEfficiency(nx);
		W_acc_err    = teffs["TEff_WEffAccVsAll" +hshl+"_WJets"  ]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_WEffAccVsAll" +hshl+"_WJets"  ]->GetEfficiencyErrorUp(nx); if(dummy>W_acc_err   ) W_acc_err    = dummy;
		Top_acc      = teffs["TEff_WEffAccVsAll" +hshl+"_Top"    ]->GetEfficiency(nx);		//contains SingleTop
		Top_acc_err  = teffs["TEff_WEffAccVsAll" +hshl+"_Top"    ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy        = teffs["TEff_WEffAccVsAll" +hshl+"_Top"    ]->GetEfficiencyErrorUp(nx); if(dummy>Top_acc_err ) Top_acc_err  = dummy;
		W_rec        = teffs["TEff_WEffRecoVsAcc"+hshl+"_WJets"  ]->GetEfficiency(nx);
		W_rec_err    = teffs["TEff_WEffRecoVsAcc"+hshl+"_WJets"  ]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_WEffRecoVsAcc"+hshl+"_WJets"  ]->GetEfficiencyErrorUp(nx); if(dummy>W_rec_err   ) W_rec_err    = dummy;
		Top_rec      = teffs["TEff_WEffRecoVsAcc"+hshl+"_Top"    ]->GetEfficiency(nx);		//contains SingleTop
		Top_rec_err  = teffs["TEff_WEffRecoVsAcc"+hshl+"_Top"    ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy        = teffs["TEff_WEffRecoVsAcc"+hshl+"_Top"    ]->GetEfficiencyErrorUp(nx); if(dummy>Top_rec_err ) Top_rec_err  = dummy;
		WT_prob      = teffs["TEff_WEffRecoVsAll"+hshl+"_WandTop"]->GetEfficiency(nx);
		WT_prob_err  = teffs["TEff_WEffRecoVsAll"+hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_WEffRecoVsAll"+hshl+"_WandTop"]->GetEfficiencyErrorUp(nx); if(dummy>WT_prob_err ) WT_prob_err  = dummy;
		WT_acc       = teffs["TEff_WEffAccVsAll" +hshl+"_WandTop"]->GetEfficiency(nx);
		WT_acc_err   = teffs["TEff_WEffAccVsAll" +hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_WEffAccVsAll" +hshl+"_WandTop"]->GetEfficiencyErrorUp(nx); if(dummy>WT_acc_err  ) WT_acc_err   = dummy;
		WT_rec       = teffs["TEff_WEffRecoVsAcc"+hshl+"_WandTop"]->GetEfficiency(nx);
		WT_rec_err   = teffs["TEff_WEffRecoVsAcc"+hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy        = teffs["TEff_WEffRecoVsAcc"+hshl+"_WandTop"]->GetEfficiencyErrorUp(nx); if(dummy>WT_rec_err  ) WT_rec_err   = dummy;
		W_prob_nomt       = teffs["TEff_WEffRecoVsAll_noMT"+hshl+"_WJets"  ]->GetEfficiency(nx);
		W_prob_err_nomt   = teffs["TEff_WEffRecoVsAll_noMT"+hshl+"_WJets"  ]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEffRecoVsAll_noMT"+hshl+"_WJets"  ]->GetEfficiencyErrorUp(nx); if(dummy>W_prob_err_nomt  ) W_prob_err_nomt   = dummy;
		Top_prob_nomt     = teffs["TEff_WEffRecoVsAll_noMT"+hshl+"_Top"    ]->GetEfficiency(nx);		//contains SingleTop
		Top_prob_err_nomt = teffs["TEff_WEffRecoVsAll_noMT"+hshl+"_Top"    ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy             = teffs["TEff_WEffRecoVsAll_noMT"+hshl+"_Top"    ]->GetEfficiencyErrorUp(nx); if(dummy>Top_prob_err_nomt) Top_prob_err_nomt = dummy;
		W_rec_nomt        = teffs["TEff_WEffRecoVsAcc_noMT"+hshl+"_WJets"  ]->GetEfficiency(nx);
		W_rec_err_nomt    = teffs["TEff_WEffRecoVsAcc_noMT"+hshl+"_WJets"  ]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEffRecoVsAcc_noMT"+hshl+"_WJets"  ]->GetEfficiencyErrorUp(nx); if(dummy>W_rec_err_nomt   ) W_rec_err_nomt    = dummy;
		Top_rec_nomt      = teffs["TEff_WEffRecoVsAcc_noMT"+hshl+"_Top"    ]->GetEfficiency(nx);		//contains SingleTop
		Top_rec_err_nomt  = teffs["TEff_WEffRecoVsAcc_noMT"+hshl+"_Top"    ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy             = teffs["TEff_WEffRecoVsAcc_noMT"+hshl+"_Top"    ]->GetEfficiencyErrorUp(nx); if(dummy>Top_rec_err_nomt ) Top_rec_err_nomt  = dummy;
		WT_prob_nomt      = teffs["TEff_WEffRecoVsAll_noMT"+hshl+"_WandTop"]->GetEfficiency(nx);
		WT_prob_err_nomt  = teffs["TEff_WEffRecoVsAll_noMT"+hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEffRecoVsAll_noMT"+hshl+"_WandTop"]->GetEfficiencyErrorUp(nx); if(dummy>WT_prob_err_nomt ) WT_prob_err_nomt  = dummy;
		WT_rec_nomt       = teffs["TEff_WEffRecoVsAcc_noMT"+hshl+"_WandTop"]->GetEfficiency(nx);
		WT_rec_err_nomt   = teffs["TEff_WEffRecoVsAcc_noMT"+hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEffRecoVsAcc_noMT"+hshl+"_WandTop"]->GetEfficiencyErrorUp(nx); if(dummy>WT_rec_err_nomt  ) WT_rec_err_nomt   = dummy;
		W_prob_noLL       = teffs["TEff_WEffRecoVsAll_noLL"+hshl+"_WJets"  ]->GetEfficiency(nx);
		W_prob_err_noLL   = teffs["TEff_WEffRecoVsAll_noLL"+hshl+"_WJets"  ]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEffRecoVsAll_noLL"+hshl+"_WJets"  ]->GetEfficiencyErrorUp(nx); if(dummy>W_prob_err_noLL  ) W_prob_err_noLL   = dummy;
		Top_prob_noLL     = teffs["TEff_WEffRecoVsAll_noLL"+hshl+"_Top"    ]->GetEfficiency(nx);		//contains SingleTop
		Top_prob_err_noLL = teffs["TEff_WEffRecoVsAll_noLL"+hshl+"_Top"    ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy             = teffs["TEff_WEffRecoVsAll_noLL"+hshl+"_Top"    ]->GetEfficiencyErrorUp(nx); if(dummy>Top_prob_err_noLL) Top_prob_err_noLL = dummy;
		W_acc_noLL        = teffs["TEff_WEffAccVsAll_noLL" +hshl+"_WJets"  ]->GetEfficiency(nx);
		W_acc_err_noLL    = teffs["TEff_WEffAccVsAll_noLL" +hshl+"_WJets"  ]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEffAccVsAll_noLL" +hshl+"_WJets"  ]->GetEfficiencyErrorUp(nx); if(dummy>W_acc_err_noLL   ) W_acc_err_noLL    = dummy;
		Top_acc_noLL      = teffs["TEff_WEffAccVsAll_noLL" +hshl+"_Top"    ]->GetEfficiency(nx);		//contains SingleTop
		Top_acc_err_noLL  = teffs["TEff_WEffAccVsAll_noLL" +hshl+"_Top"    ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy             = teffs["TEff_WEffAccVsAll_noLL" +hshl+"_Top"    ]->GetEfficiencyErrorUp(nx); if(dummy>Top_acc_err_noLL ) Top_acc_err_noLL  = dummy;
		W_rec_noLL        = teffs["TEff_WEffRecoVsAcc_noLL"+hshl+"_WJets"  ]->GetEfficiency(nx);
		W_rec_err_noLL    = teffs["TEff_WEffRecoVsAcc_noLL"+hshl+"_WJets"  ]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEffRecoVsAcc_noLL"+hshl+"_WJets"  ]->GetEfficiencyErrorUp(nx); if(dummy>W_rec_err_noLL   ) W_rec_err_noLL    = dummy;
		Top_rec_noLL      = teffs["TEff_WEffRecoVsAcc_noLL"+hshl+"_Top"    ]->GetEfficiency(nx);		//contains SingleTop
		Top_rec_err_noLL  = teffs["TEff_WEffRecoVsAcc_noLL"+hshl+"_Top"    ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy             = teffs["TEff_WEffRecoVsAcc_noLL"+hshl+"_Top"    ]->GetEfficiencyErrorUp(nx); if(dummy>Top_rec_err_noLL ) Top_rec_err_noLL  = dummy;
		WT_prob_noLL      = teffs["TEff_WEffRecoVsAll_noLL"+hshl+"_WandTop"]->GetEfficiency(nx);
		WT_prob_err_noLL  = teffs["TEff_WEffRecoVsAll_noLL"+hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEffRecoVsAll_noLL"+hshl+"_WandTop"]->GetEfficiencyErrorUp(nx); if(dummy>WT_prob_err_noLL ) WT_prob_err_noLL  = dummy;
		WT_acc_noLL       = teffs["TEff_WEffAccVsAll_noLL" +hshl+"_WandTop"]->GetEfficiency(nx);
		WT_acc_err_noLL   = teffs["TEff_WEffAccVsAll_noLL" +hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEffAccVsAll_noLL" +hshl+"_WandTop"]->GetEfficiencyErrorUp(nx); if(dummy>WT_acc_err_noLL  ) WT_acc_err_noLL   = dummy;
		WT_rec_noLL       = teffs["TEff_WEffRecoVsAcc_noLL"+hshl+"_WandTop"]->GetEfficiency(nx);
		WT_rec_err_noLL   = teffs["TEff_WEffRecoVsAcc_noLL"+hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy             = teffs["TEff_WEffRecoVsAcc_noLL"+hshl+"_WandTop"]->GetEfficiencyErrorUp(nx); if(dummy>WT_rec_err_noLL  ) WT_rec_err_noLL   = dummy;
		W_prob_nomt_noLL       = teffs["TEff_WEffRecoVsAll_noMT_noLL"+hshl+"_WJets"  ]->GetEfficiency(nx);
		W_prob_err_nomt_noLL   = teffs["TEff_WEffRecoVsAll_noMT_noLL"+hshl+"_WJets"  ]->GetEfficiencyErrorLow(nx);
		dummy                  = teffs["TEff_WEffRecoVsAll_noMT_noLL"+hshl+"_WJets"  ]->GetEfficiencyErrorUp(nx); if(dummy>W_prob_err_nomt_noLL  ) W_prob_err_nomt_noLL   = dummy;
		Top_prob_nomt_noLL     = teffs["TEff_WEffRecoVsAll_noMT_noLL"+hshl+"_Top"    ]->GetEfficiency(nx);		//contains SingleTop
		Top_prob_err_nomt_noLL = teffs["TEff_WEffRecoVsAll_noMT_noLL"+hshl+"_Top"    ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy                  = teffs["TEff_WEffRecoVsAll_noMT_noLL"+hshl+"_Top"    ]->GetEfficiencyErrorUp(nx); if(dummy>Top_prob_err_nomt_noLL) Top_prob_err_nomt_noLL = dummy;
		W_rec_nomt_noLL        = teffs["TEff_WEffRecoVsAcc_noMT_noLL"+hshl+"_WJets"  ]->GetEfficiency(nx);
		W_rec_err_nomt_noLL    = teffs["TEff_WEffRecoVsAcc_noMT_noLL"+hshl+"_WJets"  ]->GetEfficiencyErrorLow(nx);
		dummy                  = teffs["TEff_WEffRecoVsAcc_noMT_noLL"+hshl+"_WJets"  ]->GetEfficiencyErrorUp(nx); if(dummy>W_rec_err_nomt_noLL   ) W_rec_err_nomt_noLL    = dummy;
		Top_rec_nomt_noLL      = teffs["TEff_WEffRecoVsAcc_noMT_noLL"+hshl+"_Top"    ]->GetEfficiency(nx);		//contains SingleTop
		Top_rec_err_nomt_noLL  = teffs["TEff_WEffRecoVsAcc_noMT_noLL"+hshl+"_Top"    ]->GetEfficiencyErrorLow(nx);	//contains SingleTop
		dummy                  = teffs["TEff_WEffRecoVsAcc_noMT_noLL"+hshl+"_Top"    ]->GetEfficiencyErrorUp(nx); if(dummy>Top_rec_err_nomt_noLL ) Top_rec_err_nomt_noLL  = dummy;
		WT_prob_nomt_noLL      = teffs["TEff_WEffRecoVsAll_noMT_noLL"+hshl+"_WandTop"]->GetEfficiency(nx);
		WT_prob_err_nomt_noLL  = teffs["TEff_WEffRecoVsAll_noMT_noLL"+hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy                  = teffs["TEff_WEffRecoVsAll_noMT_noLL"+hshl+"_WandTop"]->GetEfficiencyErrorUp(nx); if(dummy>WT_prob_err_nomt_noLL ) WT_prob_err_nomt_noLL= dummy;
		WT_rec_nomt_noLL       = teffs["TEff_WEffRecoVsAcc_noMT_noLL"+hshl+"_WandTop"]->GetEfficiency(nx);
		WT_rec_err_nomt_noLL   = teffs["TEff_WEffRecoVsAcc_noMT_noLL"+hshl+"_WandTop"]->GetEfficiencyErrorLow(nx);
		dummy                  = teffs["TEff_WEffRecoVsAcc_noMT_noLL"+hshl+"_WandTop"]->GetEfficiencyErrorUp(nx); if(dummy>WT_rec_err_nomt_noLL  ) WT_rec_err_nomt_noLL = dummy;
		 nW = hists["TauEvents"+hshl+"_WJets"               ]->GetBinContent(nx);  nW_nomt       = hists["TauEvents_noMT"+hshl+"_WJets"         ]->GetBinContent(nx);
		nTT = hists["TauEvents"+hshl+"_TTbar"               ]->GetBinContent(nx); nTT_nomt       = hists["TauEvents_noMT"+hshl+"_TTbar"         ]->GetBinContent(nx);
		nST = hists["TauEvents"+hshl+"_SingleTop"           ]->GetBinContent(nx); nST_nomt       = hists["TauEvents_noMT"+hshl+"_SingleTop"     ]->GetBinContent(nx);
		 nT = hists["TauEvents"+hshl+"_Top"                 ]->GetBinContent(nx);  nT_nomt       = hists["TauEvents_noMT"+hshl+"_Top"           ]->GetBinContent(nx);
		nWT = hists["TauEvents"+hshl+"_WandTop"             ]->GetBinContent(nx); nWT_nomt       = hists["TauEvents_noMT"+hshl+"_WandTop"       ]->GetBinContent(nx);
		 nW_bg = hists["TauEventsBG"+hshl+"_WJets"          ]->GetBinContent(nx);  nW_bg_nomt    = hists["TauEventsBG_noMT"+hshl+"_WJets"       ]->GetBinContent(nx);
		nTT_bg = hists["TauEventsBG"+hshl+"_TTbar"          ]->GetBinContent(nx); nTT_bg_nomt    = hists["TauEventsBG_noMT"+hshl+"_TTbar"       ]->GetBinContent(nx);
		nST_bg = hists["TauEventsBG"+hshl+"_SingleTop"      ]->GetBinContent(nx); nST_bg_nomt    = hists["TauEventsBG_noMT"+hshl+"_SingleTop"   ]->GetBinContent(nx);
		 nT_bg = hists["TauEventsBG"+hshl+"_Top"            ]->GetBinContent(nx);  nT_bg_nomt    = hists["TauEventsBG_noMT"+hshl+"_Top"         ]->GetBinContent(nx);
		nWT_bg = hists["TauEventsBG"+hshl+"_WandTop"        ]->GetBinContent(nx); nWT_bg_nomt    = hists["TauEventsBG_noMT"+hshl+"_WandTop"     ]->GetBinContent(nx);
		 nW_dTau = hists["TauEvents_dTau"+hshl+"_WJets"     ]->GetBinContent(nx);  nW_dTau_nomt  = hists["TauEvents_noMT_dTau"+hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_dTau = hists["TauEvents_dTau"+hshl+"_TTbar"     ]->GetBinContent(nx); nTT_dTau_nomt  = hists["TauEvents_noMT_dTau"+hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_dTau = hists["TauEvents_dTau"+hshl+"_SingleTop" ]->GetBinContent(nx); nST_dTau_nomt  = hists["TauEvents_noMT_dTau"+hshl+"_SingleTop"]->GetBinContent(nx);
		 nT_dTau = hists["TauEvents_dTau"+hshl+"_Top"       ]->GetBinContent(nx);  nT_dTau_nomt  = hists["TauEvents_noMT_dTau"+hshl+"_Top"      ]->GetBinContent(nx);
		nWT_dTau = hists["TauEvents_dTau"+hshl+"_WandTop"   ]->GetBinContent(nx); nWT_dTau_nomt  = hists["TauEvents_noMT_dTau"+hshl+"_WandTop"  ]->GetBinContent(nx);
		 nW_LL = hists["TauEvents_LL"+hshl+"_WJets"         ]->GetBinContent(nx);  nW_LL_nomt    = hists["TauEvents_noMT_LL"+hshl+"_WJets"      ]->GetBinContent(nx);
		nTT_LL = hists["TauEvents_LL"+hshl+"_TTbar"         ]->GetBinContent(nx); nTT_LL_nomt    = hists["TauEvents_noMT_LL"+hshl+"_TTbar"      ]->GetBinContent(nx);
		nST_LL = hists["TauEvents_LL"+hshl+"_SingleTop"     ]->GetBinContent(nx); nST_LL_nomt    = hists["TauEvents_noMT_LL"+hshl+"_SingleTop"  ]->GetBinContent(nx);
		 nT_LL = hists["TauEvents_LL"+hshl+"_Top"           ]->GetBinContent(nx);  nT_LL_nomt    = hists["TauEvents_noMT_LL"+hshl+"_Top"        ]->GetBinContent(nx);
		nWT_LL = hists["TauEvents_LL"+hshl+"_WandTop"       ]->GetBinContent(nx); nWT_LL_nomt    = hists["TauEvents_noMT_LL"+hshl+"_WandTop"    ]->GetBinContent(nx);
		 nW_bg_LL = hists["TauEventsBG_LL"+hshl+"_WJets"    ]->GetBinContent(nx);  nW_bg_LL_nomt = hists["TauEventsBG_noMT_LL"+hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_bg_LL = hists["TauEventsBG_LL"+hshl+"_TTbar"    ]->GetBinContent(nx); nTT_bg_LL_nomt = hists["TauEventsBG_noMT_LL"+hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_bg_LL = hists["TauEventsBG_LL"+hshl+"_SingleTop"]->GetBinContent(nx); nST_bg_LL_nomt = hists["TauEventsBG_noMT_LL"+hshl+"_SingleTop"]->GetBinContent(nx);
		 nT_bg_LL = hists["TauEventsBG_LL"+hshl+"_Top"      ]->GetBinContent(nx);  nT_bg_LL_nomt = hists["TauEventsBG_noMT_LL"+hshl+"_Top"      ]->GetBinContent(nx);
		nWT_bg_LL = hists["TauEventsBG_LL"+hshl+"_WandTop"  ]->GetBinContent(nx); nWT_bg_LL_nomt = hists["TauEventsBG_noMT_LL"+hshl+"_WandTop"  ]->GetBinContent(nx);
		//+/- SF
		 nW_bg_SFup = hists["TauEventsBG_SFup"+hshl+"_WJets"    ]->GetBinContent(nx);  nW_bg_nomt_SFup = hists["TauEventsBG_noMT_SFup"+hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_bg_SFup = hists["TauEventsBG_SFup"+hshl+"_TTbar"    ]->GetBinContent(nx); nTT_bg_nomt_SFup = hists["TauEventsBG_noMT_SFup"+hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_bg_SFup = hists["TauEventsBG_SFup"+hshl+"_SingleTop"]->GetBinContent(nx); nST_bg_nomt_SFup = hists["TauEventsBG_noMT_SFup"+hshl+"_SingleTop"]->GetBinContent(nx);
		 nT_bg_SFup = hists["TauEventsBG_SFup"+hshl+"_Top"      ]->GetBinContent(nx);  nT_bg_nomt_SFup = hists["TauEventsBG_noMT_SFup"+hshl+"_Top"      ]->GetBinContent(nx);
		nWT_bg_SFup = hists["TauEventsBG_SFup"+hshl+"_WandTop"  ]->GetBinContent(nx); nWT_bg_nomt_SFup = hists["TauEventsBG_noMT_SFup"+hshl+"_WandTop"  ]->GetBinContent(nx);
		 nW_dTau_SFup = hists["TauEvents_dTau_SFup"+hshl+"_WJets"    ]->GetBinContent(nx);  nW_dTau_nomt_SFup = hists["TauEvents_noMT_dTau_SFup"+hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_dTau_SFup = hists["TauEvents_dTau_SFup"+hshl+"_TTbar"    ]->GetBinContent(nx); nTT_dTau_nomt_SFup = hists["TauEvents_noMT_dTau_SFup"+hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_dTau_SFup = hists["TauEvents_dTau_SFup"+hshl+"_SingleTop"]->GetBinContent(nx); nST_dTau_nomt_SFup = hists["TauEvents_noMT_dTau_SFup"+hshl+"_SingleTop"]->GetBinContent(nx);
		 nT_dTau_SFup = hists["TauEvents_dTau_SFup"+hshl+"_Top"      ]->GetBinContent(nx);  nT_dTau_nomt_SFup = hists["TauEvents_noMT_dTau_SFup"+hshl+"_Top"      ]->GetBinContent(nx);
		nWT_dTau_SFup = hists["TauEvents_dTau_SFup"+hshl+"_WandTop"  ]->GetBinContent(nx); nWT_dTau_nomt_SFup = hists["TauEvents_noMT_dTau_SFup"+hshl+"_WandTop"  ]->GetBinContent(nx);
		 nW_LL_SFup = hists["TauEvents_LL_SFup"+hshl+"_WJets"    ]->GetBinContent(nx);  nW_LL_nomt_SFup = hists["TauEvents_noMT_LL_SFup"+hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_LL_SFup = hists["TauEvents_LL_SFup"+hshl+"_TTbar"    ]->GetBinContent(nx); nTT_LL_nomt_SFup = hists["TauEvents_noMT_LL_SFup"+hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_LL_SFup = hists["TauEvents_LL_SFup"+hshl+"_SingleTop"]->GetBinContent(nx); nST_LL_nomt_SFup = hists["TauEvents_noMT_LL_SFup"+hshl+"_SingleTop"]->GetBinContent(nx);
		 nT_LL_SFup = hists["TauEvents_LL_SFup"+hshl+"_Top"      ]->GetBinContent(nx);  nT_LL_nomt_SFup = hists["TauEvents_noMT_LL_SFup"+hshl+"_Top"      ]->GetBinContent(nx);
		nWT_LL_SFup = hists["TauEvents_LL_SFup"+hshl+"_WandTop"  ]->GetBinContent(nx); nWT_LL_nomt_SFup = hists["TauEvents_noMT_LL_SFup"+hshl+"_WandTop"  ]->GetBinContent(nx);
		 nW_bg_LL_SFup = hists["TauEventsBG_LL_SFup"+hshl+"_WJets"    ]->GetBinContent(nx);  nW_bg_LL_nomt_SFup = hists["TauEventsBG_noMT_LL_SFup"+hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_bg_LL_SFup = hists["TauEventsBG_LL_SFup"+hshl+"_TTbar"    ]->GetBinContent(nx); nTT_bg_LL_nomt_SFup = hists["TauEventsBG_noMT_LL_SFup"+hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_bg_LL_SFup = hists["TauEventsBG_LL_SFup"+hshl+"_SingleTop"]->GetBinContent(nx); nST_bg_LL_nomt_SFup = hists["TauEventsBG_noMT_LL_SFup"+hshl+"_SingleTop"]->GetBinContent(nx);
		 nT_bg_LL_SFup = hists["TauEventsBG_LL_SFup"+hshl+"_Top"      ]->GetBinContent(nx);  nT_bg_LL_nomt_SFup = hists["TauEventsBG_noMT_LL_SFup"+hshl+"_Top"      ]->GetBinContent(nx);
		nWT_bg_LL_SFup = hists["TauEventsBG_LL_SFup"+hshl+"_WandTop"  ]->GetBinContent(nx); nWT_bg_LL_nomt_SFup = hists["TauEventsBG_noMT_LL_SFup"+hshl+"_WandTop"  ]->GetBinContent(nx);
		 nW_bg_SFdown = hists["TauEventsBG_SFdown"+hshl+"_WJets"    ]->GetBinContent(nx);  nW_bg_nomt_SFdown = hists["TauEventsBG_noMT_SFdown"+hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_bg_SFdown = hists["TauEventsBG_SFdown"+hshl+"_TTbar"    ]->GetBinContent(nx); nTT_bg_nomt_SFdown = hists["TauEventsBG_noMT_SFdown"+hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_bg_SFdown = hists["TauEventsBG_SFdown"+hshl+"_SingleTop"]->GetBinContent(nx); nST_bg_nomt_SFdown = hists["TauEventsBG_noMT_SFdown"+hshl+"_SingleTop"]->GetBinContent(nx);
		 nT_bg_SFdown = hists["TauEventsBG_SFdown"+hshl+"_Top"      ]->GetBinContent(nx);  nT_bg_nomt_SFdown = hists["TauEventsBG_noMT_SFdown"+hshl+"_Top"      ]->GetBinContent(nx);
		nWT_bg_SFdown = hists["TauEventsBG_SFdown"+hshl+"_WandTop"  ]->GetBinContent(nx); nWT_bg_nomt_SFdown = hists["TauEventsBG_noMT_SFdown"+hshl+"_WandTop"  ]->GetBinContent(nx);
		 nW_dTau_SFdown = hists["TauEvents_dTau_SFdown"+hshl+"_WJets"    ]->GetBinContent(nx);  nW_dTau_nomt_SFdown = hists["TauEvents_noMT_dTau_SFdown"+hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_dTau_SFdown = hists["TauEvents_dTau_SFdown"+hshl+"_TTbar"    ]->GetBinContent(nx); nTT_dTau_nomt_SFdown = hists["TauEvents_noMT_dTau_SFdown"+hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_dTau_SFdown = hists["TauEvents_dTau_SFdown"+hshl+"_SingleTop"]->GetBinContent(nx); nST_dTau_nomt_SFdown = hists["TauEvents_noMT_dTau_SFdown"+hshl+"_SingleTop"]->GetBinContent(nx);
		 nT_dTau_SFdown = hists["TauEvents_dTau_SFdown"+hshl+"_Top"      ]->GetBinContent(nx);  nT_dTau_nomt_SFdown = hists["TauEvents_noMT_dTau_SFdown"+hshl+"_Top"      ]->GetBinContent(nx);
		nWT_dTau_SFdown = hists["TauEvents_dTau_SFdown"+hshl+"_WandTop"  ]->GetBinContent(nx); nWT_dTau_nomt_SFdown = hists["TauEvents_noMT_dTau_SFdown"+hshl+"_WandTop"  ]->GetBinContent(nx);
		 nW_LL_SFdown = hists["TauEvents_LL_SFdown"+hshl+"_WJets"    ]->GetBinContent(nx);  nW_LL_nomt_SFdown = hists["TauEvents_noMT_LL_SFdown"+hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_LL_SFdown = hists["TauEvents_LL_SFdown"+hshl+"_TTbar"    ]->GetBinContent(nx); nTT_LL_nomt_SFdown = hists["TauEvents_noMT_LL_SFdown"+hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_LL_SFdown = hists["TauEvents_LL_SFdown"+hshl+"_SingleTop"]->GetBinContent(nx); nST_LL_nomt_SFdown = hists["TauEvents_noMT_LL_SFdown"+hshl+"_SingleTop"]->GetBinContent(nx);
		 nT_LL_SFdown = hists["TauEvents_LL_SFdown"+hshl+"_Top"      ]->GetBinContent(nx);  nT_LL_nomt_SFdown = hists["TauEvents_noMT_LL_SFdown"+hshl+"_Top"      ]->GetBinContent(nx);
		nWT_LL_SFdown = hists["TauEvents_LL_SFdown"+hshl+"_WandTop"  ]->GetBinContent(nx); nWT_LL_nomt_SFdown = hists["TauEvents_noMT_LL_SFdown"+hshl+"_WandTop"  ]->GetBinContent(nx);
		 nW_bg_LL_SFdown = hists["TauEventsBG_LL_SFdown"+hshl+"_WJets"    ]->GetBinContent(nx);  nW_bg_LL_nomt_SFdown = hists["TauEventsBG_noMT_LL_SFdown"+hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_bg_LL_SFdown = hists["TauEventsBG_LL_SFdown"+hshl+"_TTbar"    ]->GetBinContent(nx); nTT_bg_LL_nomt_SFdown = hists["TauEventsBG_noMT_LL_SFdown"+hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_bg_LL_SFdown = hists["TauEventsBG_LL_SFdown"+hshl+"_SingleTop"]->GetBinContent(nx); nST_bg_LL_nomt_SFdown = hists["TauEventsBG_noMT_LL_SFdown"+hshl+"_SingleTop"]->GetBinContent(nx);
		 nT_bg_LL_SFdown = hists["TauEventsBG_LL_SFdown"+hshl+"_Top"      ]->GetBinContent(nx);  nT_bg_LL_nomt_SFdown = hists["TauEventsBG_noMT_LL_SFdown"+hshl+"_Top"      ]->GetBinContent(nx);
		nWT_bg_LL_SFdown = hists["TauEventsBG_LL_SFdown"+hshl+"_WandTop"  ]->GetBinContent(nx); nWT_bg_LL_nomt_SFdown = hists["TauEventsBG_noMT_LL_SFdown"+hshl+"_WandTop"  ]->GetBinContent(nx);
		 nW_lv          = hists["NoRecoButGenTauEvents"+hshl+"_WJets"         ]->GetBinContent(nx);
		nTT_lv          = hists["NoRecoButGenTauEvents"+hshl+"_TTbar"         ]->GetBinContent(nx);
		nST_lv          = hists["NoRecoButGenTauEvents"+hshl+"_SingleTop"     ]->GetBinContent(nx);
		 nT_lv          = hists["NoRecoButGenTauEvents"+hshl+"_Top"           ]->GetBinContent(nx);
		nWT_lv          = hists["NoRecoButGenTauEvents"+hshl+"_WandTop"       ]->GetBinContent(nx);
		 nW_lv_err      = hists["NoRecoButGenTauEvents"+hshl+"_WJets"         ]->GetBinError(nx);
		nTT_lv_err      = hists["NoRecoButGenTauEvents"+hshl+"_TTbar"         ]->GetBinError(nx);
		nST_lv_err      = hists["NoRecoButGenTauEvents"+hshl+"_SingleTop"     ]->GetBinError(nx);
		 nT_lv_err      = hists["NoRecoButGenTauEvents"+hshl+"_Top"           ]->GetBinError(nx);
		nWT_lv_err      = hists["NoRecoButGenTauEvents"+hshl+"_WandTop"       ]->GetBinError(nx);
		 nW_lv_LL       = hists["NoRecoButGenTauEvents_LL"+hshl+"_WJets"      ]->GetBinContent(nx);
		nTT_lv_LL       = hists["NoRecoButGenTauEvents_LL"+hshl+"_TTbar"      ]->GetBinContent(nx);
		nST_lv_LL       = hists["NoRecoButGenTauEvents_LL"+hshl+"_SingleTop"  ]->GetBinContent(nx);
		 nT_lv_LL       = hists["NoRecoButGenTauEvents_LL"+hshl+"_Top"        ]->GetBinContent(nx);
		nWT_lv_LL       = hists["NoRecoButGenTauEvents_LL"+hshl+"_WandTop"    ]->GetBinContent(nx);
		 nW_lv_LL_err   = hists["NoRecoButGenTauEvents_LL"+hshl+"_WJets"      ]->GetBinError(nx);
		nTT_lv_LL_err   = hists["NoRecoButGenTauEvents_LL"+hshl+"_TTbar"      ]->GetBinError(nx);
		nST_lv_LL_err   = hists["NoRecoButGenTauEvents_LL"+hshl+"_SingleTop"  ]->GetBinError(nx);
		 nT_lv_LL_err   = hists["NoRecoButGenTauEvents_LL"+hshl+"_Top"        ]->GetBinError(nx);
		nWT_lv_LL_err   = hists["NoRecoButGenTauEvents_LL"+hshl+"_WandTop"    ]->GetBinError(nx);
		 nW_lv_dTau     = hists["NoRecoButGenTauEvents_dTau"+hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_lv_dTau     = hists["NoRecoButGenTauEvents_dTau"+hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_lv_dTau     = hists["NoRecoButGenTauEvents_dTau"+hshl+"_SingleTop"]->GetBinContent(nx);
		 nT_lv_dTau     = hists["NoRecoButGenTauEvents_dTau"+hshl+"_Top"      ]->GetBinContent(nx);
		nWT_lv_dTau     = hists["NoRecoButGenTauEvents_dTau"+hshl+"_WandTop"  ]->GetBinContent(nx);
		 nW_lv_dTau_err = hists["NoRecoButGenTauEvents_dTau"+hshl+"_WJets"    ]->GetBinError(nx);
		nTT_lv_dTau_err = hists["NoRecoButGenTauEvents_dTau"+hshl+"_TTbar"    ]->GetBinError(nx);
		nST_lv_dTau_err = hists["NoRecoButGenTauEvents_dTau"+hshl+"_SingleTop"]->GetBinError(nx);
		 nT_lv_dTau_err = hists["NoRecoButGenTauEvents_dTau"+hshl+"_Top"      ]->GetBinError(nx);
		nWT_lv_dTau_err = hists["NoRecoButGenTauEvents_dTau"+hshl+"_WandTop"  ]->GetBinError(nx);
		//+/- SF
		 nW_lv_SFup        = hists["NoRecoButGenTauEvents_SFup"+hshl+"_WJets"           ]->GetBinContent(nx);
		nTT_lv_SFup        = hists["NoRecoButGenTauEvents_SFup"+hshl+"_TTbar"           ]->GetBinContent(nx);
		nST_lv_SFup        = hists["NoRecoButGenTauEvents_SFup"+hshl+"_SingleTop"       ]->GetBinContent(nx);
		 nT_lv_SFup        = hists["NoRecoButGenTauEvents_SFup"+hshl+"_Top"             ]->GetBinContent(nx);
		nWT_lv_SFup        = hists["NoRecoButGenTauEvents_SFup"+hshl+"_WandTop"         ]->GetBinContent(nx);
		 nW_lv_LL_SFup     = hists["NoRecoButGenTauEvents_LL_SFup"+hshl+"_WJets"        ]->GetBinContent(nx);
		nTT_lv_LL_SFup     = hists["NoRecoButGenTauEvents_LL_SFup"+hshl+"_TTbar"        ]->GetBinContent(nx);
		nST_lv_LL_SFup     = hists["NoRecoButGenTauEvents_LL_SFup"+hshl+"_SingleTop"    ]->GetBinContent(nx);
		 nT_lv_LL_SFup     = hists["NoRecoButGenTauEvents_LL_SFup"+hshl+"_Top"          ]->GetBinContent(nx);
		nWT_lv_LL_SFup     = hists["NoRecoButGenTauEvents_LL_SFup"+hshl+"_WandTop"      ]->GetBinContent(nx);
		 nW_lv_dTau_SFup   = hists["NoRecoButGenTauEvents_dTau_SFup"+hshl+"_WJets"      ]->GetBinContent(nx);
		nTT_lv_dTau_SFup   = hists["NoRecoButGenTauEvents_dTau_SFup"+hshl+"_TTbar"      ]->GetBinContent(nx);
		nST_lv_dTau_SFup   = hists["NoRecoButGenTauEvents_dTau_SFup"+hshl+"_SingleTop"  ]->GetBinContent(nx);
		 nT_lv_dTau_SFup   = hists["NoRecoButGenTauEvents_dTau_SFup"+hshl+"_Top"        ]->GetBinContent(nx);
		nWT_lv_dTau_SFup   = hists["NoRecoButGenTauEvents_dTau_SFup"+hshl+"_WandTop"    ]->GetBinContent(nx);
		 nW_lv_SFdown      = hists["NoRecoButGenTauEvents_SFdown"+hshl+"_WJets"         ]->GetBinContent(nx);
		nTT_lv_SFdown      = hists["NoRecoButGenTauEvents_SFdown"+hshl+"_TTbar"         ]->GetBinContent(nx);
		nST_lv_SFdown      = hists["NoRecoButGenTauEvents_SFdown"+hshl+"_SingleTop"     ]->GetBinContent(nx);
		 nT_lv_SFdown      = hists["NoRecoButGenTauEvents_SFdown"+hshl+"_Top"           ]->GetBinContent(nx);
		nWT_lv_SFdown      = hists["NoRecoButGenTauEvents_SFdown"+hshl+"_WandTop"       ]->GetBinContent(nx);
		 nW_lv_LL_SFdown   = hists["NoRecoButGenTauEvents_LL_SFdown"+hshl+"_WJets"      ]->GetBinContent(nx);
		nTT_lv_LL_SFdown   = hists["NoRecoButGenTauEvents_LL_SFdown"+hshl+"_TTbar"      ]->GetBinContent(nx);
		nST_lv_LL_SFdown   = hists["NoRecoButGenTauEvents_LL_SFdown"+hshl+"_SingleTop"  ]->GetBinContent(nx);
		 nT_lv_LL_SFdown   = hists["NoRecoButGenTauEvents_LL_SFdown"+hshl+"_Top"        ]->GetBinContent(nx);
		nWT_lv_LL_SFdown   = hists["NoRecoButGenTauEvents_LL_SFdown"+hshl+"_WandTop"    ]->GetBinContent(nx);
		 nW_lv_dTau_SFdown = hists["NoRecoButGenTauEvents_dTau_SFdown"+hshl+"_WJets"    ]->GetBinContent(nx);
		nTT_lv_dTau_SFdown = hists["NoRecoButGenTauEvents_dTau_SFdown"+hshl+"_TTbar"    ]->GetBinContent(nx);
		nST_lv_dTau_SFdown = hists["NoRecoButGenTauEvents_dTau_SFdown"+hshl+"_SingleTop"]->GetBinContent(nx);
		 nT_lv_dTau_SFdown = hists["NoRecoButGenTauEvents_dTau_SFdown"+hshl+"_Top"      ]->GetBinContent(nx);
		nWT_lv_dTau_SFdown = hists["NoRecoButGenTauEvents_dTau_SFdown"+hshl+"_WandTop"  ]->GetBinContent(nx);
		QCD_bg   = hists["TauEvents"+hshl+"_QCD"      ]->GetBinContent(nx); QCD_bg_nomt   = hists["TauEvents_noMT"+hshl+"_QCD"      ]->GetBinContent(nx);
		Z_bg     = hists["TauEvents"+hshl+"_ZJets"    ]->GetBinContent(nx); Z_bg_nomt     = hists["TauEvents_noMT"+hshl+"_ZJets"    ]->GetBinContent(nx);
		Other_bg = hists["TauEvents"+hshl+"_Other"    ]->GetBinContent(nx); Other_bg_nomt = hists["TauEvents_noMT"+hshl+"_Other"    ]->GetBinContent(nx);
		TT_bg    = hists["TauEvents"+hshl+"_TTbar"    ]->GetBinContent(nx); TT_bg_nomt    = hists["TauEvents_noMT"+hshl+"_TTbar"    ]->GetBinContent(nx);
		ST_bg    = hists["TauEvents"+hshl+"_SingleTop"]->GetBinContent(nx); ST_bg_nomt    = hists["TauEvents_noMT"+hshl+"_SingleTop"]->GetBinContent(nx);
		T_bg     = hists["TauEvents"+hshl+"_Top"      ]->GetBinContent(nx); T_bg_nomt     = hists["TauEvents_noMT"+hshl+"_Top"      ]->GetBinContent(nx);
		W_bg     = hists["TauEvents"+hshl+"_WJets"    ]->GetBinContent(nx); W_bg_nomt     = hists["TauEvents_noMT"+hshl+"_WJets"    ]->GetBinContent(nx);
		WT_bg    = hists["TauEvents"+hshl+"_WandTop"  ]->GetBinContent(nx); WT_bg_nomt    = hists["TauEvents_noMT"+hshl+"_WandTop"  ]->GetBinContent(nx);
		nonWT_bg = hists["TauEvents"+hshl+"_noWandTop"]->GetBinContent(nx); nonWT_bg_nomt = hists["TauEvents_noMT"+hshl+"_noWandTop"]->GetBinContent(nx);
		//+/- SF
		QCD_bg_SFup     = hists["TauEvents_SFup"+hshl+"_QCD"        ]->GetBinContent(nx); QCD_bg_nomt_SFup     = hists["TauEvents_noMT_SFup"+hshl+"_QCD"        ]->GetBinContent(nx);
		Z_bg_SFup       = hists["TauEvents_SFup"+hshl+"_ZJets"      ]->GetBinContent(nx); Z_bg_nomt_SFup       = hists["TauEvents_noMT_SFup"+hshl+"_ZJets"      ]->GetBinContent(nx);
		Other_bg_SFup   = hists["TauEvents_SFup"+hshl+"_Other"      ]->GetBinContent(nx); Other_bg_nomt_SFup   = hists["TauEvents_noMT_SFup"+hshl+"_Other"      ]->GetBinContent(nx);
		TT_bg_SFup      = hists["TauEvents_SFup"+hshl+"_TTbar"      ]->GetBinContent(nx); TT_bg_nomt_SFup      = hists["TauEvents_noMT_SFup"+hshl+"_TTbar"      ]->GetBinContent(nx);
		ST_bg_SFup      = hists["TauEvents_SFup"+hshl+"_SingleTop"  ]->GetBinContent(nx); ST_bg_nomt_SFup      = hists["TauEvents_noMT_SFup"+hshl+"_SingleTop"  ]->GetBinContent(nx);
		T_bg_SFup       = hists["TauEvents_SFup"+hshl+"_Top"        ]->GetBinContent(nx); T_bg_nomt_SFup       = hists["TauEvents_noMT_SFup"+hshl+"_Top"        ]->GetBinContent(nx);
		W_bg_SFup       = hists["TauEvents_SFup"+hshl+"_WJets"      ]->GetBinContent(nx); W_bg_nomt_SFup       = hists["TauEvents_noMT_SFup"+hshl+"_WJets"      ]->GetBinContent(nx);
		WT_bg_SFup      = hists["TauEvents_SFup"+hshl+"_WandTop"    ]->GetBinContent(nx); WT_bg_nomt_SFup      = hists["TauEvents_noMT_SFup"+hshl+"_WandTop"    ]->GetBinContent(nx);
		nonWT_bg_SFup   = hists["TauEvents_SFup"+hshl+"_noWandTop"  ]->GetBinContent(nx); nonWT_bg_nomt_SFup   = hists["TauEvents_noMT_SFup"+hshl+"_noWandTop"  ]->GetBinContent(nx);
		QCD_bg_SFdown   = hists["TauEvents_SFdown"+hshl+"_QCD"      ]->GetBinContent(nx); QCD_bg_nomt_SFdown   = hists["TauEvents_noMT_SFdown"+hshl+"_QCD"      ]->GetBinContent(nx);
		Z_bg_SFdown     = hists["TauEvents_SFdown"+hshl+"_ZJets"    ]->GetBinContent(nx); Z_bg_nomt_SFdown     = hists["TauEvents_noMT_SFdown"+hshl+"_ZJets"    ]->GetBinContent(nx);
		Other_bg_SFdown = hists["TauEvents_SFdown"+hshl+"_Other"    ]->GetBinContent(nx); Other_bg_nomt_SFdown = hists["TauEvents_noMT_SFdown"+hshl+"_Other"    ]->GetBinContent(nx);
		TT_bg_SFdown    = hists["TauEvents_SFdown"+hshl+"_TTbar"    ]->GetBinContent(nx); TT_bg_nomt_SFdown    = hists["TauEvents_noMT_SFdown"+hshl+"_TTbar"    ]->GetBinContent(nx);
		ST_bg_SFdown    = hists["TauEvents_SFdown"+hshl+"_SingleTop"]->GetBinContent(nx); ST_bg_nomt_SFdown    = hists["TauEvents_noMT_SFdown"+hshl+"_SingleTop"]->GetBinContent(nx);
		T_bg_SFdown     = hists["TauEvents_SFdown"+hshl+"_Top"      ]->GetBinContent(nx); T_bg_nomt_SFdown     = hists["TauEvents_noMT_SFdown"+hshl+"_Top"      ]->GetBinContent(nx);
		W_bg_SFdown     = hists["TauEvents_SFdown"+hshl+"_WJets"    ]->GetBinContent(nx); W_bg_nomt_SFdown     = hists["TauEvents_noMT_SFdown"+hshl+"_WJets"    ]->GetBinContent(nx);
		WT_bg_SFdown    = hists["TauEvents_SFdown"+hshl+"_WandTop"  ]->GetBinContent(nx); WT_bg_nomt_SFdown    = hists["TauEvents_noMT_SFdown"+hshl+"_WandTop"  ]->GetBinContent(nx);
		nonWT_bg_SFdown = hists["TauEvents_SFdown"+hshl+"_noWandTop"]->GetBinContent(nx); nonWT_bg_nomt_SFdown = hists["TauEvents_noMT_SFdown"+hshl+"_noWandTop"]->GetBinContent(nx);

		nData     = hists["TauEvents"+hshl+"_data"]->GetBinContent(nx); nData_nomt     = hists["TauEvents_noMT"+hshl+"_data"]->GetBinContent(nx);
		MC_bg     = hists["TauEvents"+hshl+"_mc"  ]->GetBinContent(nx); MC_bg_nomt     = hists["TauEvents_noMT"+hshl+"_mc"  ]->GetBinContent(nx);
		MC_bg_err = hists["TauEvents"+hshl+"_mc"  ]->GetBinError(nx);   MC_bg_err_nomt = hists["TauEvents_noMT"+hshl+"_mc"  ]->GetBinError(nx);


		//the code below is very ugly
		//on the other hand it is very flexible as all flags set at the beginning of this macros are considered.
		if(fUseMTcut){
			//only used for full printouts
			Top_acc = Top_acc*mteff;
			W_acc   = W_acc  *mteff;
			WT_acc   = WT_acc *mteff;
			Top_acc_err = sqrt(pow(Top_acc_err*mteff,2) + pow(Top_acc*mtefferr, 2));
			W_acc_err   = sqrt(pow(  W_acc_err*mteff,2) + pow(  W_acc*mtefferr, 2));
			WT_acc_err  = sqrt(pow( WT_acc_err*mteff,2) + pow( WT_acc*mtefferr, 2));
		}

		//get variables
		//cases of !fUseMTcut or !fVetoLL are dealt with later
		double nW_scaled(0.);//number of true one lepton events (i.e. with a gen W-->lnu)
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_scaled = nWT;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_scaled = nT;
		else{ if(fIncludeTop) nW_scaled += nTT;	if(fIncludeSingleTop) nW_scaled += nST;	if(!fTopOnly) nW_scaled += nW;	}
		double nW_bg_scaled(0.);//background to the number above
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_scaled = nWT_bg;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_scaled = nT_bg;
		else{ if(fIncludeTop) nW_bg_scaled += nTT_bg;	if(fIncludeSingleTop) nW_bg_scaled += nST_bg;	if(!fTopOnly) nW_bg_scaled += nW_bg;	}
		double nW_bg_scaled_SFup(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_scaled_SFup = nWT_bg_SFup;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_scaled_SFup = nT_bg_SFup;
		else{ if(fIncludeTop) nW_bg_scaled_SFup += nTT_bg_SFup;	if(fIncludeSingleTop) nW_bg_scaled_SFup += nST_bg_SFup;	if(!fTopOnly) nW_bg_scaled_SFup += nW_bg_SFup;	}
		double nW_bg_scaled_SFdown(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_scaled_SFdown = nWT_bg_SFdown;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_scaled_SFdown = nT_bg_SFdown;
		else{ if(fIncludeTop) nW_bg_scaled_SFdown += nTT_bg_SFdown;	if(fIncludeSingleTop) nW_bg_scaled_SFdown += nST_bg_SFdown;	if(!fTopOnly) nW_bg_scaled_SFdown += nW_bg_SFdown;	}
		double nW_lv_scaled(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_lv_scaled = nWT_lv;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_lv_scaled = nT_lv;
		else{ if(fIncludeTop) nW_lv_scaled += nTT_lv;	if(fIncludeSingleTop) nW_lv_scaled += nST_lv;	if(!fTopOnly) nW_lv_scaled += nW_lv;	}
		double nW_lv_scaled_err(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_lv_scaled_err = pow(nWT_lv_err,2);
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_lv_scaled_err = pow(nT_lv_err,2);
		else{ if(fIncludeTop) nW_lv_scaled_err += pow(nTT_lv_err,2);	if(fIncludeSingleTop) nW_lv_scaled_err += pow(nST_lv_err,2);	if(!fTopOnly) nW_lv_scaled_err += pow(nW_lv_err,2);	}
		nW_lv_scaled_err = sqrt(nW_lv_scaled_err);
		double nW_lv_scaled_SFup(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_lv_scaled_SFup = nWT_lv_SFup;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_lv_scaled_SFup = nT_lv_SFup;
		else{ if(fIncludeTop) nW_lv_scaled_SFup += nTT_lv_SFup;	if(fIncludeSingleTop) nW_lv_scaled_SFup += nST_lv_SFup;	if(!fTopOnly) nW_lv_scaled_SFup += nW_lv_SFup;	}
		double nW_lv_scaled_SFdown(0.);
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_lv_scaled_SFdown = nWT_lv_SFdown;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_lv_scaled_SFdown = nT_lv_SFdown;
		else{ if(fIncludeTop) nW_lv_scaled_SFdown += nTT_lv_SFdown;	if(fIncludeSingleTop) nW_lv_scaled_SFdown += nST_lv_SFdown;	if(!fTopOnly) nW_lv_scaled_SFdown += nW_lv_SFdown;	}
		double nW_lv_scaled_SFerr = fabs(nW_lv_scaled-nW_lv_scaled_SFup) > fabs(nW_lv_scaled-nW_lv_scaled_SFdown) ? fabs(nW_lv_scaled-nW_lv_scaled_SFup) : fabs(nW_lv_scaled-nW_lv_scaled_SFdown);
		double nW_dTau_scaled(0.);//double genhadtau
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_dTau_scaled = nWT_dTau;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_dTau_scaled = nT_dTau;
		else{ if(fIncludeTop) nW_dTau_scaled += nTT_dTau;	if(fIncludeSingleTop) nW_dTau_scaled += nST_dTau;	if(!fTopOnly) nW_dTau_scaled += nW_dTau;	}
		double nW_dTau_scaled_SFup(0.);//double genhadtau
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_dTau_scaled_SFup = nWT_dTau_SFup;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_dTau_scaled_SFup = nT_dTau_SFup;
		else{ if(fIncludeTop) nW_dTau_scaled_SFup += nTT_dTau_SFup;	if(fIncludeSingleTop) nW_dTau_scaled_SFup += nST_dTau_SFup;	if(!fTopOnly) nW_dTau_scaled_SFup += nW_dTau_SFup;	}
		double nW_dTau_scaled_SFdown(0.);//double genhadtau
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_dTau_scaled_SFdown = nWT_dTau_SFdown;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_dTau_scaled_SFdown = nT_dTau_SFdown;
		else{ if(fIncludeTop) nW_dTau_scaled_SFdown += nTT_dTau_SFdown;	if(fIncludeSingleTop) nW_dTau_scaled_SFdown += nST_dTau_SFdown;	if(!fTopOnly) nW_dTau_scaled_SFdown += nW_dTau_SFdown;	}
		double nW_lv_dTau_scaled(0.);//double genhadtau
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_lv_dTau_scaled = nWT_lv_dTau;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_lv_dTau_scaled = nT_lv_dTau;
		else{ if(fIncludeTop) nW_lv_dTau_scaled += nTT_lv_dTau;	if(fIncludeSingleTop) nW_lv_dTau_scaled += nST_lv_dTau;	if(!fTopOnly) nW_lv_dTau_scaled += nW_lv_dTau;	}
		double nW_lv_dTau_scaled_err(0.);//double genhadtau
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_lv_dTau_scaled_err = pow(nWT_lv_dTau_err,2);
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_lv_dTau_scaled_err = pow(nT_lv_dTau_err,2);
		else{ if(fIncludeTop) nW_lv_dTau_scaled_err += pow(nTT_lv_dTau_err,2);	if(fIncludeSingleTop) nW_lv_dTau_scaled_err += pow(nST_lv_dTau_err,2);	if(!fTopOnly) nW_lv_dTau_scaled_err += pow(nW_lv_dTau_err,2);	}
		nW_lv_dTau_scaled_err = sqrt(nW_lv_dTau_scaled_err);
		double nW_lv_dTau_scaled_SFup(0.);//double genhadtau
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_lv_dTau_scaled_SFup = nWT_lv_dTau_SFup;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_lv_dTau_scaled_SFup = nT_lv_dTau_SFup;
		else{ if(fIncludeTop) nW_lv_dTau_scaled_SFup += nTT_lv_dTau_SFup;	if(fIncludeSingleTop) nW_lv_dTau_scaled_SFup += nST_lv_dTau_SFup;	if(!fTopOnly) nW_lv_dTau_scaled_SFup += nW_lv_dTau_SFup;	}
		double nW_lv_dTau_scaled_SFdown(0.);//double genhadtau
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_lv_dTau_scaled_SFdown = nWT_lv_dTau_SFdown;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_lv_dTau_scaled_SFdown = nT_lv_dTau_SFdown;
		else{ if(fIncludeTop) nW_lv_dTau_scaled_SFdown += nTT_lv_dTau_SFdown;	if(fIncludeSingleTop) nW_lv_dTau_scaled_SFdown += nST_lv_dTau_SFdown;	if(!fTopOnly) nW_lv_dTau_scaled_SFdown += nW_lv_dTau_SFdown;	}
		double nW_lv_dTau_scaled_SFerr = fabs(nW_lv_dTau_scaled-nW_lv_dTau_scaled_SFup) > fabs(nW_lv_dTau_scaled-nW_lv_dTau_scaled_SFdown) ? fabs(nW_lv_dTau_scaled-nW_lv_dTau_scaled_SFup) : fabs(nW_lv_dTau_scaled-nW_lv_dTau_scaled_SFdown);
		double nW_LL_scaled(0.);//events having a LostLepton (e,mu)
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_LL_scaled = nWT_LL;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_LL_scaled = nT_LL;
		else{ if(fIncludeTop) nW_LL_scaled += nTT_LL;	if(fIncludeSingleTop) nW_LL_scaled += nST_LL;	if(!fTopOnly) nW_LL_scaled += nW_LL;	}
		double nW_LL_scaled_SFup(0.);//events having a LostLepton (e,mu)
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_LL_scaled_SFup = nWT_LL_SFup;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_LL_scaled_SFup = nT_LL_SFup;
		else{ if(fIncludeTop) nW_LL_scaled_SFup += nTT_LL_SFup;	if(fIncludeSingleTop) nW_LL_scaled_SFup += nST_LL_SFup;	if(!fTopOnly) nW_LL_scaled_SFup += nW_LL_SFup;	}
		double nW_LL_scaled_SFdown(0.);//events having a LostLepton (e,mu)
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_LL_scaled_SFdown = nWT_LL_SFdown;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_LL_scaled_SFdown = nT_LL_SFdown;
		else{ if(fIncludeTop) nW_LL_scaled_SFdown += nTT_LL_SFdown;	if(fIncludeSingleTop) nW_LL_scaled_SFdown += nST_LL_SFdown;	if(!fTopOnly) nW_LL_scaled_SFdown += nW_LL_SFdown;	}
		double nW_bg_LL_scaled(0.);//events having a LostLepton (e,mu)
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_LL_scaled = nWT_bg_LL;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_LL_scaled = nT_bg_LL;
		else{ if(fIncludeTop) nW_bg_LL_scaled += nTT_bg_LL;	if(fIncludeSingleTop) nW_bg_LL_scaled += nST_bg_LL;	if(!fTopOnly) nW_bg_LL_scaled += nW_bg_LL;	}
		double nW_bg_LL_scaled_SFup(0.);//events having a LostLepton (e,mu)
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_LL_scaled_SFup = nWT_bg_LL_SFup;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_LL_scaled_SFup = nT_bg_LL_SFup;
		else{ if(fIncludeTop) nW_bg_LL_scaled_SFup += nTT_bg_LL_SFup;	if(fIncludeSingleTop) nW_bg_LL_scaled_SFup += nST_bg_LL_SFup;	if(!fTopOnly) nW_bg_LL_scaled_SFup += nW_bg_LL_SFup;	}
		double nW_bg_LL_scaled_SFdown(0.);//events having a LostLepton (e,mu)
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_LL_scaled_SFdown = nWT_bg_LL_SFdown;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_LL_scaled_SFdown = nT_bg_LL_SFdown;
		else{ if(fIncludeTop) nW_bg_LL_scaled_SFdown += nTT_bg_LL_SFdown;	if(fIncludeSingleTop) nW_bg_LL_scaled_SFdown += nST_bg_LL_SFdown;	if(!fTopOnly) nW_bg_LL_scaled_SFdown += nW_bg_LL_SFdown;	}
		double nW_lv_LL_scaled(0.);//events having a LostLepton (e,mu)
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_lv_LL_scaled = nWT_lv_LL;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_lv_LL_scaled = nT_lv_LL;
		else{ if(fIncludeTop) nW_lv_LL_scaled += nTT_lv_LL;	if(fIncludeSingleTop) nW_lv_LL_scaled += nST_lv_LL;	if(!fTopOnly) nW_lv_LL_scaled += nW_lv_LL;	}
		double nW_lv_LL_scaled_err(0.);//events having a LostLepton (e,mu)
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_lv_LL_scaled_err = pow(nWT_lv_LL_err,2);
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_lv_LL_scaled_err = pow(nT_lv_LL_err,2);
		else{ if(fIncludeTop) nW_lv_LL_scaled_err += pow(nTT_lv_LL_err,2);	if(fIncludeSingleTop) nW_lv_LL_scaled_err += pow(nST_lv_LL_err,2);	if(!fTopOnly) nW_lv_LL_scaled_err += pow(nW_lv_LL_err,2);	}
		nW_lv_LL_scaled_err = sqrt(nW_lv_LL_scaled_err);
		double nW_lv_LL_scaled_SFup(0.);//events having a LostLepton (e,mu)
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_lv_LL_scaled_SFup = nWT_lv_LL_SFup;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_lv_LL_scaled_SFup = nT_lv_LL_SFup;
		else{ if(fIncludeTop) nW_lv_LL_scaled_SFup += nTT_lv_LL_SFup;	if(fIncludeSingleTop) nW_lv_LL_scaled_SFup += nST_lv_LL_SFup;	if(!fTopOnly) nW_lv_LL_scaled_SFup += nW_lv_LL_SFup;	}
		double nW_lv_LL_scaled_SFdown(0.);//events having a LostLepton (e,mu)
		if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_lv_LL_scaled_SFdown = nWT_lv_LL_SFdown;
		else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_lv_LL_scaled_SFdown = nT_lv_LL_SFdown;
		else{ if(fIncludeTop) nW_lv_LL_scaled_SFdown += nTT_lv_LL_SFdown;	if(fIncludeSingleTop) nW_lv_LL_scaled_SFdown += nST_lv_LL_SFdown;	if(!fTopOnly) nW_lv_LL_scaled_SFdown += nW_lv_LL_SFdown;	}
		double nW_lv_LL_scaled_SFerr = fabs(nW_lv_LL_scaled-nW_lv_LL_scaled_SFup) > fabs(nW_lv_LL_scaled-nW_lv_LL_scaled_SFdown) ? fabs(nW_lv_LL_scaled-nW_lv_LL_scaled_SFup) : fabs(nW_lv_LL_scaled-nW_lv_LL_scaled_SFdown);
		if(!fUseMTcut){
			if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_scaled = nWT_nomt;
			else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_scaled = nT_nomt;
			else{ if(fIncludeTop) nW_scaled += nTT_nomt;	if(fIncludeSingleTop) nW_scaled += nST_nomt;	if(!fTopOnly) nW_scaled += nW_nomt;	}
			if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_scaled = nWT_bg_nomt;
			else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_scaled = nT_bg_nomt;
			else{ if(fIncludeTop) nW_bg_scaled += nTT_bg_nomt;	if(fIncludeSingleTop) nW_bg_scaled += nST_bg_nomt;	if(!fTopOnly) nW_bg_scaled += nW_bg_nomt;	}
			if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_scaled_SFup = nWT_bg_nomt_SFup;
			else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_scaled_SFup = nT_bg_nomt_SFup;
			else{ if(fIncludeTop) nW_bg_scaled_SFup += nTT_bg_nomt_SFup;	if(fIncludeSingleTop) nW_bg_scaled_SFup += nST_bg_nomt_SFup;	if(!fTopOnly) nW_bg_scaled_SFup += nW_bg_nomt_SFup;	}
			if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_scaled_SFdown = nWT_bg_nomt_SFdown;
			else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_scaled_SFdown = nT_bg_nomt_SFdown;
			else{ if(fIncludeTop) nW_bg_scaled_SFdown += nTT_bg_nomt_SFdown;	if(fIncludeSingleTop) nW_bg_scaled_SFdown += nST_bg_nomt_SFdown;	if(!fTopOnly) nW_bg_scaled_SFdown += nW_bg_nomt_SFdown;	}
			//double genhadtau
			if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_dTau_scaled = nWT_dTau_nomt;
			else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_dTau_scaled = nT_dTau_nomt;
			else{ if(fIncludeTop) nW_dTau_scaled += nTT_dTau_nomt;	if(fIncludeSingleTop) nW_dTau_scaled += nST_dTau_nomt;	if(!fTopOnly) nW_dTau_scaled += nW_dTau_nomt;	}
			if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_dTau_scaled_SFup = nWT_dTau_nomt_SFup;
			else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_dTau_scaled_SFup = nT_dTau_nomt_SFup;
			else{ if(fIncludeTop) nW_dTau_scaled_SFup += nTT_dTau_nomt_SFup;	if(fIncludeSingleTop) nW_dTau_scaled_SFup += nST_dTau_nomt_SFup;	if(!fTopOnly) nW_dTau_scaled_SFup += nW_dTau_nomt_SFup;	}
			if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_dTau_scaled_SFdown = nWT_dTau_nomt_SFdown;
			else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_dTau_scaled_SFdown = nT_dTau_nomt_SFdown;
			else{ if(fIncludeTop) nW_dTau_scaled_SFdown += nTT_dTau_nomt_SFdown;	if(fIncludeSingleTop) nW_dTau_scaled_SFdown += nST_dTau_nomt_SFdown;	if(!fTopOnly) nW_dTau_scaled_SFdown += nW_dTau_nomt_SFdown;	}
			//events having a LostLepton (e,mu)
			if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_LL_scaled = nWT_LL_nomt;
			else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_LL_scaled = nT_LL_nomt;
			else{ if(fIncludeTop) nW_LL_scaled += nTT_LL_nomt;	if(fIncludeSingleTop) nW_LL_scaled += nST_LL_nomt;	if(!fTopOnly) nW_LL_scaled += nW_LL_nomt;	}
			if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_LL_scaled_SFup = nWT_LL_nomt_SFup;
			else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_LL_scaled_SFup = nT_LL_nomt_SFup;
			else{ if(fIncludeTop) nW_LL_scaled_SFup += nTT_LL_nomt_SFup;	if(fIncludeSingleTop) nW_LL_scaled_SFup += nST_LL_nomt_SFup;	if(!fTopOnly) nW_LL_scaled_SFup += nW_LL_nomt_SFup;	}
			if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_LL_scaled_SFdown = nWT_LL_nomt_SFdown;
			else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_LL_scaled_SFdown = nT_LL_nomt_SFdown;
			else{ if(fIncludeTop) nW_LL_scaled_SFdown += nTT_LL_nomt_SFdown;	if(fIncludeSingleTop) nW_LL_scaled_SFdown += nST_LL_nomt_SFdown;	if(!fTopOnly) nW_LL_scaled_SFdown += nW_LL_nomt_SFdown;	}
			if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_LL_scaled = nWT_bg_LL_nomt;
			else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_LL_scaled = nT_bg_LL_nomt;
			else{ if(fIncludeTop) nW_bg_LL_scaled += nTT_bg_LL_nomt;	if(fIncludeSingleTop) nW_bg_LL_scaled += nST_bg_LL_nomt;	if(!fTopOnly) nW_bg_LL_scaled += nW_bg_LL_nomt;	}
			if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_LL_scaled_SFup = nWT_bg_LL_nomt_SFup;
			else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_LL_scaled_SFup = nT_bg_LL_nomt_SFup;
			else{ if(fIncludeTop) nW_bg_LL_scaled_SFup += nTT_bg_LL_nomt_SFup;	if(fIncludeSingleTop) nW_bg_LL_scaled_SFup += nST_bg_LL_nomt_SFup;	if(!fTopOnly) nW_bg_LL_scaled_SFup += nW_bg_LL_nomt_SFup;	}
			if(      fIncludeTop && fIncludeSingleTop && !fTopOnly) nW_bg_LL_scaled_SFdown = nWT_bg_LL_nomt_SFdown;
			else if( fIncludeTop && fIncludeSingleTop &&  fTopOnly) nW_bg_LL_scaled_SFdown = nT_bg_LL_nomt_SFdown;
			else{ if(fIncludeTop) nW_bg_LL_scaled_SFdown += nTT_bg_LL_nomt_SFdown;	if(fIncludeSingleTop) nW_bg_LL_scaled_SFdown += nST_bg_LL_nomt_SFdown;	if(!fTopOnly) nW_bg_LL_scaled_SFdown += nW_bg_LL_nomt_SFdown;	}
		}


		//keep MT and non-MT prob separate, but not noLL/LL
		double prob(-1.);
		if(fWeightedProb)            prob = WT_prob_noLL;
		else{	if(fTopEfficencies)  prob = Top_prob_noLL;
			else                 prob = W_prob_noLL;
		}
		double prob_err_sys;
		//factor 2 on rel_sys_uncert, one for T&P, one for MT cut
		if(fWeightedProb)           prob_err_sys = sqrt(pow(2.*rel_sys_uncert*WT_prob_noLL, 2) + pow(WT_prob_err_noLL, 2));
		else{	if(fTopEfficencies) prob_err_sys = sqrt(pow(2.*rel_sys_uncert*Top_prob_noLL,2) + pow(Top_prob_err_noLL,2));
			else                prob_err_sys = sqrt(pow(2.*rel_sys_uncert*W_prob_noLL,  2) + pow(W_prob_err_noLL,  2));
		}
		double prob_MC_err_sys;
		//no factor on rel_sys_uncert, as MC closure
		if(fWeightedProb)           prob_MC_err_sys = WT_prob_err_noLL;
		else{	if(fTopEfficencies) prob_MC_err_sys = Top_prob_err_noLL;
			else                prob_MC_err_sys = W_prob_err_noLL;
		}
		double prob_nomt(-1.);
		if(fWeightedProb)         prob_nomt = WT_prob_nomt_noLL;
		else{ if(fTopEfficencies) prob_nomt = Top_prob_nomt_noLL;
			else              prob_nomt = W_prob_nomt_noLL;
		}
		double prob_nomt_err_sys;
		//factor 2 on rel_sys_uncert, one for T&P, one for MT cut
		if(fWeightedProb)         prob_nomt_err_sys = sqrt(pow(rel_sys_uncert*WT_prob_nomt_noLL, 2) + pow(WT_prob_err_nomt_noLL, 2));
		else{ if(fTopEfficencies) prob_nomt_err_sys = sqrt(pow(rel_sys_uncert*Top_prob_nomt_noLL,2) + pow(Top_prob_err_nomt_noLL,2));
			else              prob_nomt_err_sys = sqrt(pow(rel_sys_uncert*W_prob_nomt_noLL,  2) + pow(W_prob_err_nomt_noLL,  2));
		}
		double prob_MC_nomt_err_sys;
		//no factor on rel_sys_uncert, as MC closure
		if(fWeightedProb)         prob_MC_nomt_err_sys = WT_prob_err_nomt_noLL;
		else{ if(fTopEfficencies) prob_MC_nomt_err_sys = Top_prob_err_nomt_noLL;
			else              prob_MC_nomt_err_sys = W_prob_err_nomt_noLL;
		}

		if(!fVetoLL){
			if(fWeightedProb)         prob = WT_prob;
			else{ if(fTopEfficencies) prob = Top_prob;
				else              prob = W_prob; }
			if(fWeightedProb)           prob_err_sys = sqrt(pow(2.*rel_sys_uncert*WT_prob, 2) + pow(WT_prob_err, 2));
			else{	if(fTopEfficencies) prob_err_sys = sqrt(pow(2.*rel_sys_uncert*Top_prob,2) + pow(Top_prob_err,2));
				else                prob_err_sys = sqrt(pow(2.*rel_sys_uncert*W_prob,  2) + pow(W_prob_err,  2)); }
			if(fWeightedProb)           prob_MC_err_sys = WT_prob_err;
			else{	if(fTopEfficencies) prob_MC_err_sys = Top_prob_err;
				else                prob_MC_err_sys = W_prob_err; }
			if(fWeightedProb)           prob_nomt = WT_prob_nomt;
			else{	if(fTopEfficencies) prob_nomt = Top_prob_nomt;
				else                prob_nomt = W_prob_nomt; }
			if(fWeightedProb)           prob_nomt_err_sys = sqrt(pow(rel_sys_uncert*WT_prob_nomt, 2) + pow(WT_prob_err_nomt, 2));
			else{	if(fTopEfficencies) prob_nomt_err_sys = sqrt(pow(rel_sys_uncert*Top_prob_nomt,2) + pow(Top_prob_err_nomt,2));
				else                prob_nomt_err_sys = sqrt(pow(rel_sys_uncert*W_prob_nomt,  2) + pow(W_prob_err_nomt,  2)); }
			if(fWeightedProb)           prob_MC_nomt_err_sys = WT_prob_err_nomt;
			else{	if(fTopEfficencies) prob_MC_nomt_err_sys = Top_prob_err_nomt;
				else                prob_MC_nomt_err_sys = W_prob_err_nomt; }
		}
//cases:
//1) veto taus (default) or not //what do I predict
//2) veto LL (default) or not // define rel error on LL veto
//3) use MT cut to reduce signal contamination or not (to be tested)
//4) predict also doubleLostTaus 
//--> 2x2x2 = 8 cases * 2 (data/MC) = 16 implementations into 2 formula (pred and predMC)
		if(fDoubleTau_BG){//does not effect in any way LL, but noLL
			nW_bg_scaled        += nW_dTau_scaled;		//double genhadtau added to background
			nW_bg_scaled_SFup   += nW_dTau_scaled_SFup;
			nW_bg_scaled_SFdown += nW_dTau_scaled_SFdown;
			//if fUseDoubleTau==true: WEvents and NoRecoButGenTaus contains doubleTau, else not
			if(fUseDoubleTau) {
				//if event is doubleTau it cannot have LL (for TTbar like (not for TTV)!
				//nW_lv_scaled contains double gentaus --> need to subtract it here
				nW_lv_scaled       -= nW_lv_dTau_scaled;
				double nW_lv_err2   = TMath::Max(0.,(pow(nW_lv_scaled_err,  2)-pow(nW_lv_dTau_scaled_err,  2) ) );//is this correct???
				nW_lv_scaled_err    = nW_lv_err2;
				nW_lv_scaled_SFup   = TMath::Max(0.,nW_lv_scaled_SFup  -nW_lv_dTau_scaled_SFup);
				nW_lv_scaled_SFdown = TMath::Max(0.,nW_lv_scaled_SFdown-nW_lv_dTau_scaled_SFdown);
				nW_lv_scaled_SFerr  = fabs(nW_lv_scaled-nW_lv_scaled_SFup) > fabs(nW_lv_scaled-nW_lv_scaled_SFdown) ? fabs(nW_lv_scaled-nW_lv_scaled_SFup) : fabs(nW_lv_scaled-nW_lv_scaled_SFdown);
				//nW_lv_noLL_scaled does not exist - is calculated with [nW_lv_scaled-nW_lv_LL_scaled] and nW_lv_LL_scaled does not need it
			}
		} else {
			//double genhad tau is signal - does not effect LL
			//if fUseDoubleTau==true: WEvents and NoRecoButGenTaus contains dTau, else not
			if(!fUseDoubleTau){
				//if event is doubleTau it cannot have LL (for TTbar like (not for TTV)!
				nW_lv_scaled        += nW_lv_dTau_scaled;
				nW_lv_scaled_err     = sqrt(pow(nW_lv_scaled_err,2)+pow(nW_lv_dTau_scaled_err,2));
				nW_lv_scaled_SFup   += nW_dTau_scaled_SFup;
				nW_lv_scaled_SFdown += nW_dTau_scaled_SFdown;
				nW_lv_scaled_SFerr   = fabs(nW_lv_scaled-nW_lv_scaled_SFup) > fabs(nW_lv_scaled-nW_lv_scaled_SFdown) ? fabs(nW_lv_scaled-nW_lv_scaled_SFup) : fabs(nW_lv_scaled-nW_lv_scaled_SFdown);
				//nW_lv_noLL_scaled does not exist - is calculated with [nW_lv_scaled-nW_lv_LL_scaled] and nW_lv_LL_scaled does not need it
			}
		}
		double allMC(0.), allMC_SFup(0.), allMC_SFdown(0.), allMC_err(0.);
		double numberData(0.);
		numberData = nData;
		allMC        = WT_bg + nonWT_bg;//contains all cases, below, sometimes single top missing
		allMC_SFup   = WT_bg_SFup   + nonWT_bg_SFup;
		allMC_SFdown = WT_bg_SFdown + nonWT_bg_SFdown;
		allMC_err    = MC_bg_err;
		double allMC_SFerr = fabs(allMC-allMC_SFup) > fabs(allMC-allMC_SFdown) ? fabs(allMC-allMC_SFup) : fabs(allMC-allMC_SFdown);
		double R_LL(0.), R_LL_err(0.), R_LL_mc(0.), R_LL_mc_err(0.);
		R_LL     = (1.-prob_nomt)/(prob_nomt*mteff);
		R_LL_err = sqrt(pow((prob_nomt_err_sys/(prob_nomt*prob_nomt*mteff)),2)
			   +    pow((rel_sys_uncert*(1.-prob_nomt)/(prob_nomt*mteff)),2)
			   +    pow((frelMTerr*(1.-prob_nomt)/(prob_nomt*mteff*mteff)),2) );
		if(mtefferr<1) R_LL_err = sqrt(pow(R_LL_err,2) + pow((mtefferr*(1.-prob_nomt)/(prob_nomt*mteff*mteff)),2));
		R_LL_mc_err = sqrt(pow((prob_MC_nomt_err_sys/(prob_nomt*prob_nomt*mteff)),2)
			      +    pow((rel_sys_uncert*(1.-prob_nomt)/(prob_nomt*mteff)),2) );
		if(mtefferr<1) R_LL_mc_err = sqrt(pow(R_LL_mc_err,2) + pow((mtefferr*(1.-prob_nomt)/(prob_nomt*mteff*mteff)),2));
		R_LL_mc     = (1-prob_nomt)/prob;
		double bg(0.), bg_SFup(0.), bg_SFdown(0.);//also contains LL for //all 1 reco tau events where there is no genhadtau from W/Top decay
		bg        = nW_bg_scaled        + Z_bg        + QCD_bg        + Other_bg;
		bg_SFup   = nW_bg_scaled_SFup   + Z_bg_SFup   + QCD_bg_SFup   + Other_bg_SFup;
		bg_SFdown = nW_bg_scaled_SFdown + Z_bg_SFdown + QCD_bg_SFdown + Other_bg_SFdown;
		if( fTopOnly)  bg        +=  W_bg;        if(!fIncludeTop) bg        += TT_bg;        if(!fIncludeSingleTop) bg        += ST_bg;
		if( fTopOnly)  bg_SFup   +=  W_bg_SFup;   if(!fIncludeTop) bg_SFup   += TT_bg_SFup;   if(!fIncludeSingleTop) bg_SFup   += ST_bg_SFup;
		if( fTopOnly)  bg_SFdown +=  W_bg_SFdown; if(!fIncludeTop) bg_SFdown += TT_bg_SFdown; if(!fIncludeSingleTop) bg_SFdown += ST_bg_SFdown;
		if(!fUseMTcut) {
			numberData = nData_nomt;
			allMC        = WT_bg_nomt + nonWT_bg_nomt;
			allMC_SFup   = WT_bg_nomt_SFup   + nonWT_bg_nomt_SFup;
			allMC_SFdown = WT_bg_nomt_SFdown + nonWT_bg_nomt_SFdown;
			allMC_err    = MC_bg_err_nomt;
			R_LL     = (1.-prob_nomt)/(prob_nomt);
			R_LL_err = sqrt(pow((prob_nomt_err_sys/(prob_nomt*prob_nomt)),2)
				   +    pow((rel_sys_uncert*(1.-prob_nomt)/(prob_nomt)),2) );
			R_LL_mc_err = sqrt(pow((prob_MC_nomt_err_sys/(prob_nomt*prob_nomt)),2)
				      +    pow((rel_sys_uncert*(1.-prob_nomt)/(prob_nomt)),2) );
			R_LL_mc     = R_LL;
			bg        = nW_bg_scaled        + Z_bg_nomt        + QCD_bg_nomt        + Other_bg_nomt;
			bg_SFup   = nW_bg_scaled_SFup   + Z_bg_nomt_SFup   + QCD_bg_nomt_SFup   + Other_bg_nomt_SFup;
			bg_SFdown = nW_bg_scaled_SFdown + Z_bg_nomt_SFdown + QCD_bg_nomt_SFdown + Other_bg_nomt_SFdown;
			if( fTopOnly) bg        +=  W_bg_nomt;        if(!fIncludeTop) bg        += TT_bg_nomt;        if(!fIncludeSingleTop) bg        += ST_bg_nomt;
			if( fTopOnly) bg_SFup   +=  W_bg_nomt_SFup;   if(!fIncludeTop) bg_SFup   += TT_bg_nomt_SFup;   if(!fIncludeSingleTop) bg_SFup   += ST_bg_nomt_SFup;
			if( fTopOnly) bg_SFdown +=  W_bg_nomt_SFdown; if(!fIncludeTop) bg_SFdown += TT_bg_nomt_SFdown; if(!fIncludeSingleTop) bg_SFdown += ST_bg_nomt_SFdown;
		}
		double bg_SFerr = fabs(bg-bg_SFup) > fabs(bg-bg_SFdown) ? fabs(bg-bg_SFup) : fabs(bg-bg_SFdown);
		//only BG for //all   1 reco tau events containing LostLepton - those with no gentau (no double counting)
		//nW_bg_LL_scaled is already contained in bg --> need to subtract it from bg_KK
		// MT cut already applied before
		double bg_LL = nW_LL_scaled - nW_bg_LL_scaled;// Z, QCD, Others don't contain LL
		double bg_LL_SFup = nW_LL_scaled_SFup - nW_bg_LL_scaled_SFup;
		double bg_LL_SFdown = nW_LL_scaled_SFdown - nW_bg_LL_scaled_SFdown;
		double bg_LL_SFerr = fabs(bg_LL-bg_LL_SFup) > fabs(bg_LL-bg_LL_SFdown) ? fabs(bg_LL-bg_LL_SFup) : fabs(bg_LL-bg_LL_SFdown);
		double bgk(0.), bgkerr(0.);//bg for mt/nomt; dTau/nodTau already done, bg does not change for veto/noveto taus
		if(fVetoLL){//ll added to background
			bgk = bg + bg_LL;//THINK HERE IF REL BG ERROR SHOULD CONTAIN BTAGGING ERROR????
			if(fbTagError) bgkerr = sqrt(pow(bg_SFerr,2)+pow(bg_LL_SFerr,2)+pow(rel_sys_uncert_bg*bg,2)+pow(bg_LL*fLLError,2));
			else           bgkerr = sqrt(pow(rel_sys_uncert_bg*bg,2)+pow(bg_LL*fLLError,2));
		} else {//do nothing
			bgk = bg;
			if(fbTagError) bgkerr = sqrt(pow(bg_SFerr,2) + pow(rel_sys_uncert_bg*bg,2));
			else           bgkerr = sqrt(pow(rel_sys_uncert_bg*bg,2));
		}
		if(allMC>0){
			bgk = bgk * numberData / allMC;
			bgkerr = bgkerr * numberData / allMC;
		}

		double pred    = (numberData - bgk) * R_LL;//all cases already implemented here
		double pred_MC = (allMC   - bgk) * R_LL_mc;
		double pred_error_stat    = fabs(sqrt(numberData)*R_LL);
		double pred_MC_error_stat = 0.;
		double pred_error_sys    = sqrt(pow((numberData - bgk)*R_LL_err,2) + pow(bgkerr*R_LL,2) + pow(nWT_goodrecoevt_dTau*R_LL_err*fdoubeLLerr,2 ));
		double pred_MC_error_sys =  0.;
		if(fUseMTcut)       pred_MC_error_sys = sqrt(pow((allMC   - bgk)*R_LL_mc_err,2) + pow(MC_bg_err     *R_LL_mc,2));
		else                pred_MC_error_sys = sqrt(pow((allMC   - bgk)*R_LL_mc_err,2) + pow(MC_bg_err_nomt*R_LL_mc,2));

		if(fbTagError){
			nW_lv_LL_scaled_err        = sqrt(pow(nW_lv_LL_scaled_err,2)+pow(nW_lv_LL_scaled_SFerr,2));
			nW_lv_dTau_scaled_err      = sqrt(pow(nW_lv_dTau_scaled_err,2)+pow(nW_lv_dTau_scaled_SFerr,2));
			nW_lv_scaled_err           = sqrt(pow(nW_lv_scaled_err,2)+pow(nW_lv_scaled_SFerr,2));
			allMC_err                  = sqrt(pow(allMC_err,2)+pow(allMC_SFerr,2));
		}

	datapred.push_back(pred);
	datapred_stat_err.push_back(pred_error_stat);
	datapred_syst_err.push_back(pred_error_sys );
	MCpred.push_back(pred_MC);
	MCpred_stat_err.push_back(pred_MC_error_stat );
	MCpred_syst_err.push_back(pred_MC_error_sys );
	MTeff.push_back(mteff);
	MTefferr.push_back(mtefferr);
	RLL.push_back(R_LL);
	RLLerr.push_back(R_LL_err);
	if(fVetoLL){
		mctruth.push_back(nW_lv_scaled-nW_lv_LL_scaled);
		double mctrutherrtemp = TMath::Max(0., pow(nW_lv_scaled_err,2) - pow(nW_lv_LL_scaled_err,2));
		mctrutherr.push_back(sqrt(mctrutherrtemp));
		mctruthTT.push_back(nTT_lv-nTT_lv_LL);
		mctruthW.push_back(  nW_lv- nW_lv_LL);
	} else {
		mctruth.push_back(nW_lv_scaled);
		mctrutherr.push_back(nW_lv_scaled_err);
		mctruthTT.push_back(nTT_lv);
		mctruthW.push_back(  nW_lv);
	}
	numtrueW.push_back(nW_scaled);
	numtrueW_bg.push_back(nW_bg_scaled);
	numData.push_back(numberData);
	numBG.push_back(bgk);
	numBGLL.push_back(bg_LL);//maybe here add bgk error due to LL
	if(fUseMTcut){
		LLeff.push_back(prob_nomt);
		LLefferr.push_back(prob_err_sys);
		numQCD.push_back(QCD_bg);
		numZ.push_back(Z_bg);
		numW.push_back(W_bg);
		numT.push_back(T_bg);
		numTT.push_back(TT_bg);
		numOther.push_back(Other_bg);
	} else {
		LLeff.push_back(prob);
		LLefferr.push_back(prob_nomt_err_sys);
		numQCD.push_back(QCD_bg_nomt);
		numZ.push_back(Z_bg_nomt);
		numW.push_back(W_bg_nomt);
		numT.push_back(T_bg_nomt);
		numTT.push_back(TT_bg_nomt);
		numOther.push_back(Other_bg_nomt);
	}
	numMC.push_back(allMC);
	BGerr.push_back(allMC_err);


	//this is there for debugging issues only
	//unfortunately some of the namings are wrong.
	//but you should define these printouts by yourselve, so I did not correct this
	if(makeFullPrintout){
		*fLogStream << "*******************************************" << endl;
		*fLogStream << "************* FULL PRINTOUT ***************" << endl;
		*fLogStream << "*******************************************" << endl;
		*fLogStream << "Signal region " << signal_region[i2] << ", HT bin " << HT_bin[i3] << " and MT2bin " << nx << " = (" << hists[string("TauEvents"+hshl+"_mc")]->GetBinLowEdge(nx) << "," << hists[string("TauEvents"+hshl+"_mc")]->GetBinLowEdge(nx) + hists[string("TauEvents"+hshl+"_mc")]->GetBinWidth(nx) << "):" << endl;

		*fLogStream << "1tau yield:       QCD = " << QCD_bg << ", Z = " << Z_bg << ", Other = " << Other_bg << " ==> non WandTop = " << nonWT_bg << endl;
		*fLogStream << "                  TTbar = " << TT_bg << ", SingleTop = " << ST_bg << ", W = " << W_bg << " --> Top = " << T_bg << " ==> WandTop = " << WT_bg << endl;
		*fLogStream << "                  Data = " << nData <<  ", total mc = " << MC_bg << "+/-" << MC_bg_err << endl;
		*fLogStream << "1tau yield, noMT: QCD = " << QCD_bg_nomt << ", Z = " << Z_bg_nomt << ", Other = " << Other_bg_nomt << " ==> non WandTop = " << nonWT_bg_nomt << endl;
		*fLogStream << "                  TTbar = " << TT_bg_nomt << ", SingleTop = " << ST_bg_nomt << ", W = " << W_bg_nomt << " --> Top = " << T_bg_nomt << " ==> WandTop = " << WT_bg_nomt << endl;
		*fLogStream << "                  Data = " << nData_nomt <<  ", total mc = " << MC_bg_nomt << "+/-" << MC_bg_err_nomt << endl;
		*fLogStream << "1tau yield SFdown:QCD = " << QCD_bg_SFdown << ", Z = " << Z_bg_SFdown << ", Other = " << Other_bg_SFdown << " ==> non WandTop = " << nonWT_bg_SFdown << endl;
		*fLogStream << "                  TTbar = " << TT_bg_SFdown << ", SingleTop = " << ST_bg_SFdown << ", W = " << W_bg_SFdown << " --> Top = " << T_bg_SFdown << " ==> WandTop = " << WT_bg_SFdown << endl;
		*fLogStream << "1tau yield SFup:  QCD = " << QCD_bg_SFup << ", Z = " << Z_bg_SFup << ", Other = " << Other_bg_SFup << " ==> non WandTop = " << nonWT_bg_SFup << endl;
		*fLogStream << "                  TTbar = " << TT_bg_SFup << ", SingleTop = " << ST_bg_SFup << ", W = " << W_bg_SFup << " --> Top = " << T_bg_SFup << " ==> WandTop = " << WT_bg_SFup << endl;
		*fLogStream << "1tau events:         TTbar = " << nTT << ", SingleTop = " << nST << ", W = " << nW << " --> Top = " << nT << " ==> WandTop = " << nWT << endl;
		*fLogStream << "1 tau evts noMT:     TTbar = " << nTT_nomt << ", SingleTop = " << nST_nomt << ", W = " << nW_nomt << " --> Top = " << nT_nomt << " ==> WandTop = " << nWT_nomt << endl;
		*fLogStream << "1 tau evts LL:       TTbar = " << nTT_LL << ", SingleTop = " << nST_LL << ", W = " << nW_LL << " --> Top = " << nT_LL << " ==> WandTop = " << nWT_LL << endl;
		*fLogStream << "1 tau evts noMT LL:  TTbar = " << nTT_LL_nomt << ", SingleTop = " << nST_LL_nomt << ", W = " << nW_LL_nomt << " --> Top = " << nT_LL_nomt << " ==> WandTop = " << nWT_LL_nomt << endl;
		*fLogStream << "1tau events without 1gentau from W:" << endl;
		*fLogStream << "            TTbar = " << nTT_bg << ", SingleTop = " << nST_bg << ", W = " << nW_bg << " --> Top = " << nT_bg << " ==> WandTop = " << nWT_bg << endl;
		*fLogStream << "noMT:       TTbar = " << nTT_bg_nomt << ", SingleTop = " << nST_bg_nomt << ", W = " << nW_bg_nomt << " --> Top = " << nT_bg_nomt << " ==> WandTop = " << nWT_bg_nomt << endl;
		*fLogStream << "LL              TTbar = " << nTT_bg_LL << ", SingleTop = " << nST_bg_LL << ", W = " << nW_bg_LL << " --> Top = " << nT_bg_LL << " ==> WandTop = " << nWT_bg_LL << endl;
		*fLogStream << "noMT: LL:   TTbar = " << nTT_bg_LL_nomt << ", SingleTop = " << nST_bg_LL_nomt << ", W = " << nW_bg_LL_nomt << " --> Top = " << nT_bg_LL_nomt << " ==> WandTop = " << nWT_bg_LL_nomt << endl;
		*fLogStream << "0tau but with 1gentau from W:" << endl;
		*fLogStream << "        TTbar = " << nTT_lv << "+/-" << nTT_lv_err << ", SingleTop = " << nST_lv << "+/-" << nST_lv_err << ", W = " << nW_lv << "+/-" << nW_lv_err << " --> Top = " << nT_lv << "+/-" << nT_lv_err << " ==> WandTop = " << nWT_lv << "+/-" << nWT_lv_err << endl;
		*fLogStream << "LL:     TTbar = " << nTT_lv_LL << "+/-" << nTT_lv_LL_err << ", SingleTop = " << nST_lv_LL << "+/-" << nST_lv_LL_err << ", W = " << nW_lv_LL << "+/-" << nW_lv_LL_err << " --> Top = " << nT_lv_LL << "+/-" << nT_lv_LL_err << " ==> WandTop = " << nWT_lv_LL << "+/-" << nWT_lv_LL_err << endl;
		*fLogStream << "Events with 1gentau from W:" << endl;
		*fLogStream << "events no reco tau but two gen tau: " << nWT_lv_dTau << endl;
		*fLogStream << "events  1 reco tau but two gen tau: " << nWT_dTau << endl;
		*fLogStream << endl;
		*fLogStream << "MT efficiency:               MC = " << mteff << "+/-" << mtefferr << ", data = " << mteffdata << "+/-" << mtefferrdata << endl;
		*fLogStream << "MT efficiency2:           W+Top = " << mteff2 << "+/-" << mtefferr2 << ", noLL   = " << mteff2_noLL << "+/-" << mtefferr2_noLL << endl;
		*fLogStream << "Acceptance:                   W = " << W_acc << "+/-" << W_acc_err << ", Top = " << Top_acc << "+/-" << Top_acc_err << ", WandTop = " << WT_acc << "+/-" << WT_acc_err << endl;
		*fLogStream << "Reco Efficiency:              W = " << W_rec << "+/-" << W_rec_err << ", Top = " << Top_rec << "+/-" << Top_rec_err << ", WandTop = " << WT_rec << "+/-" << WT_rec_err << endl;
		*fLogStream << "Total Efficiency:             W = " << W_prob << "+/-" << W_prob_err << ", Top = " << Top_prob << "+/-" << Top_prob_err << ", WandTop = " << WT_prob << "+/-" << WT_prob_err << endl;
		*fLogStream << "Reco Efficiency, noMT:        W = " << W_rec_nomt << "+/-" << W_rec_err_nomt << ", Top = " << Top_rec_nomt << "+/-" << Top_rec_err_nomt << ", WandTop = " << WT_rec_nomt << "+/-" << WT_rec_err_nomt << endl;
		*fLogStream << "Total Efficiency, noMT:       W = " << W_prob_nomt << "+/-" << W_prob_err_nomt << ", Top = " << Top_prob_nomt << "+/-" << Top_prob_err_nomt << ", WandTop = " << WT_prob_nomt << "+/-" << WT_prob_err_nomt << endl;
		*fLogStream << "Acceptance, noLL:             W = " << W_acc_noLL << "+/-" << W_acc_err_noLL << ", Top = " << Top_acc_noLL << "+/-" << Top_acc_err_noLL << ", WandTop = " << WT_acc_noLL << "+/-" << WT_acc_err_noLL << endl;
		*fLogStream << "Reco Efficiency, noLL:        W = " << W_rec_noLL << "+/-" << W_rec_err_noLL << ", Top = " << Top_rec_noLL << "+/-" << Top_rec_err_noLL << ", WandTop = " << WT_rec_noLL << "+/-" << WT_rec_err_noLL << endl;
		*fLogStream << "Total Efficiency, noLL:       W = " << W_prob_noLL << "+/-" << W_prob_err_noLL << ", Top = " << Top_prob_noLL << "+/-" << Top_prob_err_noLL << ", WandTop = " << WT_prob_noLL << "+/-" << WT_prob_err_noLL << endl;
		*fLogStream << "Reco Efficiency, noMT, noLL:  W = " << W_rec_nomt_noLL << "+/-" << W_rec_err_nomt_noLL << ", Top = " << Top_rec_nomt_noLL << "+/-" << Top_rec_err_nomt_noLL << ", WandTop = " << WT_rec_nomt_noLL << "+/-" << WT_rec_err_nomt_noLL << endl;
		*fLogStream << "Total Efficiency, noMT, noLL: W = " << W_prob_nomt_noLL << "+/-" << W_prob_err_nomt_noLL << ", Top = " << Top_prob_nomt_noLL << "+/-" << Top_prob_err_nomt_noLL << ", WandTop = " << WT_prob_nomt_noLL << "+/-" << WT_prob_err_nomt_noLL << endl;
		*fLogStream << "Tau Purity        = " << taupur << "+/-" << taupurerr << endl;
		*fLogStream << "Tau Efficiency    = " << taueff << "+/-" << tauefferr << endl;
		*fLogStream << "Jet Tau Fake Rate = " << taufake << "+/-" << taufakeerr << endl;
		*fLogStream << endl;
		*fLogStream << "Number of leptonic events from W decay LL = " << nW_LL_scaled << endl;
		*fLogStream << "Number of leptonic events from W decay = " << nW_scaled << "; LL: " << nW_LL_scaled << endl;
		*fLogStream << "Background to that number (i.e leptonic events that are not true genlept) = " << nW_bg_scaled << " (SFdown/up = " << nW_bg_scaled_SFdown << "/" << nW_bg_scaled_SFup << endl;
		*fLogStream << "Total background    = " << bg << " (SFerr = " << bg_SFerr << " due to SFdown/up = " << bg_SFdown << "/" << bg_SFup << "); LL = " << bg_LL << endl;
		*fLogStream << "Total MC            = " << allMC << endl;
		*fLogStream << "prob (epsilon)      = " << prob << "+/-" << prob_err_sys << endl;
		*fLogStream << "prob (epsilon) noMT = " << prob_nomt << "+/-" << prob_nomt_err_sys << endl;
		*fLogStream << "       R_LL         = " << (1.-prob_nomt)/(prob_nomt*mteff) << "+/-" << fabs(prob_nomt_err_sys/(prob_nomt*prob_nomt*mteff)) << endl;
		*fLogStream << "       R_LL (noMT)  = " << (1.-prob_nomt)/(prob_nomt) << "+/-" << fabs(prob_nomt_err_sys/(prob_nomt*prob_nomt)) << endl;
		*fLogStream << "prediction     MC   = " << pred_MC << " +/- " << pred_MC_error_stat << "(stat) +/- " << pred_MC_error_sys << "(syst)" << endl;
		*fLogStream << "prediction     data = " << pred << " +/- " << pred_error_stat << "(stat) +/- " << pred_error_sys << "(syst)" << endl;
		*fLogStream << "MC truth of lost taus       = " << nW_lv_scaled << " +/- " << nW_lv_scaled_err << endl;
		*fLogStream << "MC truth of lost taus, noLL = " << nW_lv_scaled-nW_lv_LL_scaled << " +/- " << sqrt(nW_lv_scaled_err*nW_lv_scaled_err-nW_lv_LL_scaled_err*nW_lv_LL_scaled_err) << endl;
		*fLogStream << "*******************************************" << endl << endl << endl;
	}
	
	}//for(int nx = 1; nx<=hists[string("RecoTauEvents"+hshl+"_mc")]->GetNbinsX(); ++nx)
	}}// i2, i3, i4

	//this makes a compact card of signal region vs MCtruth +/- err vs. pred. +/- err
	if(PrintPredictionCard) PredictionCard(sr, htr, MT2bin, mctruth, mctrutherr, datapred, datapred_stat_err, datapred_syst_err, LLeff, LLefferr, RLL, RLLerr, MTeff, MTefferr);
	//this makes a compact card of signal region vs MCtruth +/- err vs. pred. +/- err, and also prints out (from MC truth) the fraction of expected lost lepton from W,Top,rest(dibosons?)
	if(PrintPredictionCard) PredictionCardSplitted(sr, htr, MT2bin, mctruth, mctrutherr, datapred, datapred_stat_err, datapred_syst_err, LLeff, LLefferr, RLL, RLLerr, MTeff, MTefferr, mctruthW, mctruthTT, numW, numTT, numMC);
	//this just makes the one lepton yield table (usually with MT cut applied for one lepton selection, unless flags are set differently)
	if(PrintYieldTable) YieldTable(sr, htr, MT2low, MT2up, numQCD, numZ, numW, numT, numOther, numMC, BGerr, numData, numBGLL);
	//this line produces the final result tables
	if(PrintSummaryTable) SummaryTable(version, sr, htr, MT2low, MT2up, numtrueW, numtrueW_bg, numData, numBG, numBGLL, LLeff, LLefferr, RLL, RLLerr, MTeff, MTefferr, mctruth, mctrutherr, MCpred, MCpred_stat_err, MCpred_syst_err, datapred, datapred_stat_err, datapred_syst_err);
	//this function stores results into a root file, also can produce MCtruth / prediction plots
	if(SavePrediction) PredictionFile(version, sr, htr, MT2low, MT2up, numtrueW, numtrueW_bg, numData, numBG, numBGLL, LLeff, LLefferr, RLL, RLLerr, MTeff, MTefferr, mctruth, mctrutherr, MCpred, MCpred_stat_err, MCpred_syst_err, datapred, datapred_stat_err, datapred_syst_err, PlotPrediction);

}


//version == 0,1: with MCPred
//version == 2,3: without MCPred
//version == 0,2: with R_LL (i.e. 1-e / e)
//version == 1,3: with e (e = efficiency)
//this function makes the final result tables
void SummaryTable(int version, vector<int> sr, vector<int> htr, vector<double> MT2low, vector<double> MT2up, vector<double> numtrueW, vector<double> numtrueW_bg, vector<double> numData, vector<double> numBG, vector<double> numBGLL, vector<double> LLeff, vector<double> LLefferr, vector<double> RLL, vector<double> RLLerr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruth, vector<double> mctrutherr, vector<double> MCpred, vector<double> MCpred_stat_err, vector<double> MCpred_syst_err, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err){

	*fLogStream << "*********************************************************************" << endl;
	*fLogStream << "\%BEGINLATEX\%"             << endl;
	*fLogStream << "\\begin{table}"             << endl
	     << "\\begin{center}"            << endl
	     << "\\small"                    << endl;
       	*fLogStream << "\\begin{tabular}{lccccccc}" << endl;	     
	*fLogStream << "\\hline\\hline"             << endl;

	if(!fRebin) *fLogStream << "$M_{T2}$ (GeV) & ";
	else        *fLogStream << "signal region       & ";
	if(fIncludeTop && !fTopOnly) *fLogStream << "$N^{MC}(W \\& Top)"<< "$ &  $";
	if(fTopOnly)                 *fLogStream << "N^{MC}(Top)"<< "$ &  $";
	if(!fIncludeTop)             *fLogStream << "N^{MC}(W)"<< "$ &  $";
	                             *fLogStream << "N^{reco}" << "$ &  $" << "N^{bg}"  <<  "$ ";
	if(fPlotLLTab)               *fLogStream << "      ";
	                             *fLogStream << " &     $";  
	if(version==1||version==3)   *fLogStream <<  "\\varepsilon"  << "$ &    $" << "N^{pass}$ MC      " << " & $";
	if(version==0||version==2)   *fLogStream <<  "R_{LL}"   << "$    &      $" << "N^{pass}$ MC      " << " &     $";
	if(version==0||version==1)   *fLogStream << "N^{pass}$ MCPred    " << " &                    $";
	                             *fLogStream << "N^{pass}$ Pred                      " << " \\\\" << endl;
	string oldlep = "dummy";
	string oldsr = "dummy";
	string oldHT = "dummy";
	for(unsigned int n = 0; n<sr.size(); ++n){
		string sigreg;
		if(sr[n]==0) sigreg="2j, 0b";
		if(sr[n]==1) sigreg="2j, $\\ge 1$b";
		if(sr[n]==2) sigreg="$3-5$j, 0b";
		if(sr[n]==3) sigreg="$3-5$j, 1b";
		if(sr[n]==4) sigreg="$3-5$j, 2b";
		if(sr[n]==5) sigreg="$\\ge 6$j, 0b";
		if(sr[n]==6) sigreg="$\\ge 6$j, 1b";
		if(sr[n]==7) sigreg="$\\ge 6$j, 2b";
		if(sr[n]==8) sigreg="$\\ge 3$j, $\\ge 3$b";
		string htreg;
		if(htr[n]==0 && fMET) htreg = "450 GeV $\\leq H_{T} < 750$ GeV";
		if(htr[n]==1 && fHT) htreg = "750 GeV $\\leq H_{T} < 1200$ GeV";
		if(htr[n]==2 && fHT) htreg = "$H_{T}\\ge 1200$ GeV";
		if(htreg!=oldHT){//this if makes sure that HT region is only printed out, when the vector jumps to the next HT region
			*fLogStream << " \\hline " << endl << "\\multicolumn{7}{l}{" <<  htreg << "} \\\\ " << endl << "\\hline" << endl;
			oldHT = htreg;
		}
		if(!fRebin){
		if(sigreg!=oldsr){//this if makes sure that topological region is only printed out, when the vector jumps topological the next HT region
			*fLogStream << " \\hline  " << endl << sigreg << "\\\\" /* << endl << "\\hline"*/ << endl;
			oldsr = sigreg;
		}
		if(MT2up[n]==10000.) *fLogStream << "$" << int(MT2low[n]) << "-" << "\\infty$" << " " << setw(4) << " & ";
		else            *fLogStream << "$" << int(MT2low[n]) << "-" << int(MT2up[n]) << "$" << " " << setw(7) << " & ";
		}
		else *fLogStream << " " << setw(18) << sigreg << " & ";

		if(fVetoLL) *fLogStream << fixed << setprecision(2) << " " << setw(17) << numtrueW[n]-numtrueW_bg[n]-numBGLL[n] << " & ";
		else        *fLogStream << fixed << setprecision(2) << " " << setw(17) << numtrueW[n]-numtrueW_bg[n] << " & ";
		*fLogStream << " " << setw(10) << int(numData[n]) << " & ";
		if(fPlotLLTab) *fLogStream << fixed << setprecision(2) << " " << setw(8) << numBG[n] << " (" << numBGLL[n] << ") & ";
		else *fLogStream << fixed << setprecision(2) << " " << setw(9) << numBG[n] << " & ";
		if(version==1||version==3) *fLogStream << fixed << setprecision(2) << "$" << LLeff[n] << " \\pm " << LLefferr[n] << "$" << " & ";//MT EFFICIENCY MISSING
		if(version==0||version==2) *fLogStream << fixed << setprecision(2) << "$" << RLL[n] << " \\pm " << RLLerr[n] << "$" << " & ";
		*fLogStream << "$" << " " << setw(9) << mctruth[n] << " \\pm " << " " << setw(6) <<  mctrutherr[n] << "$ & ";
		if(version==0||version==1) *fLogStream << fixed << setprecision(2) << "$" << " " << setw(9) << MCpred[n] << " \\pm " << " " << setw(7) << MCpred_syst_err[n] << "$" <<  " & ";
        	*fLogStream << fixed << setprecision(2) << " " << " " << setw(9) << datapred[n] << " $\\pm$ " << " " << setw(8) << datapred_stat_err[n] << " (stat) $\\pm$ " << " " << setw(8) << datapred_syst_err[n] << " (sys)"   << " \\\\" << endl;//think of adding syst + stat
	}
	*fLogStream << "\\hline\\hline"                                                                                                << endl
	     << "\\end{tabular}"                                                                                                << endl
	     << "\\end{center}"                                                                                                 << endl
	     << "\\end{table}"                                                                                                  << endl
	     << "\%ENDLATEX\%"                                                                                                  << endl
	     << endl;
	*fLogStream << endl << endl;
}

//stores prediction and MC truth into a file and possibly makes plots
void PredictionFile(int version, vector<int> sr, vector<int> htr, vector<double> MT2low, vector<double> MT2up, vector<double> numtrueW, vector<double> numtrueW_bg, vector<double> numData, vector<double> numBG, vector<double> numBGLL, vector<double> LLeff, vector<double> LLefferr, vector<double> RLL, vector<double> RLLerr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruth, vector<double> mctrutherr, vector<double> MCpred, vector<double> MCpred_stat_err, vector<double> MCpred_syst_err, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, Bool_t PlotPrediction){

	TString filename = outputdir;
	TString plotdirectory = outputdir + "plots/";
	if(fRebin) plotdirectory = plotdirectory + "plots/";
	else       plotdirectory = plotdirectory + "plotsMT2binned/";
	if(PlotPrediction) Util::MakeOutputDir(plotdirectory);
	if(fRebin) filename = filename + "LostTauPredictionFile.root";
	else       filename = filename + "FineBinnedLostTauPredictionFile.root";


	map<string, TH1D*>    hs;
	//for now it is only truth +/- error
	//and prediction +/- error stored
	for(int i3 = 0; i3<HTbinsize;        ++i3){
		string htreg;
		if(i3==0 && fMET) htreg = "lowHT";
		if(i3==0 && fHT) htreg = "mediumHT";
		if(i3==1 && fHT) htreg = "highHT";
		if(fRebin){
			string mapname;
			mapname = "Prediction_Tau_"+htreg;
			if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", 9, 0, 9); hs[mapname]->Sumw2();
			hs[mapname]->GetXaxis()->SetTitle("signal region");
			hs[mapname]->GetXaxis()->SetBinLabel(1,"2j,0b");     hs[mapname]->GetXaxis()->SetBinLabel(2,"2j,#geq1b"); hs[mapname]->GetXaxis()->SetBinLabel(9,"#geq3j,#geq3b");
			hs[mapname]->GetXaxis()->SetBinLabel(3,"3-5j,0b");   hs[mapname]->GetXaxis()->SetBinLabel(4,"3-5j,1b");   hs[mapname]->GetXaxis()->SetBinLabel(5,"3-5j,2b");
			hs[mapname]->GetXaxis()->SetBinLabel(6,"#geq6j,0b"); hs[mapname]->GetXaxis()->SetBinLabel(7,"#geq6j,1b"); hs[mapname]->GetXaxis()->SetBinLabel(8,"#geq6j,2b");
			hs[mapname]->SetMarkerStyle(20);
			hs[mapname]->SetMarkerColor(kBlack);
			hs[mapname]->SetLineColor(kBlack);
			hs[mapname]->SetLineWidth(3);
			mapname = "SimulationTruth_Tau_"+htreg;
			if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", 9, 0, 9); hs[mapname]->Sumw2();
			hs[mapname]->GetXaxis()->SetTitle("signal region");
			hs[mapname]->GetXaxis()->SetBinLabel(1,"2j,0b");     hs[mapname]->GetXaxis()->SetBinLabel(2,"2j,#geq1b"); hs[mapname]->GetXaxis()->SetBinLabel(9,"#geq3j,#geq3b");
			hs[mapname]->GetXaxis()->SetBinLabel(3,"3-5j,0b");   hs[mapname]->GetXaxis()->SetBinLabel(4,"3-5j,1b");   hs[mapname]->GetXaxis()->SetBinLabel(5,"3-5j,2b");
			hs[mapname]->GetXaxis()->SetBinLabel(6,"#geq6j,0b"); hs[mapname]->GetXaxis()->SetBinLabel(7,"#geq6j,1b"); hs[mapname]->GetXaxis()->SetBinLabel(8,"#geq6j,2b");
			hs[mapname]->SetMarkerColor(kBlue);
			hs[mapname]->SetLineColor(kBlue);
			hs[mapname]->SetFillColor(kBlue);
			hs[mapname]->SetLineWidth(0);
			hs[mapname]->SetFillStyle(3002); 
		}
	for(int i2 = 0; i2<signalregionsize; ++i2){
		if(fRebin&&i4==0&&i3==0){
			string sigreg;
			if(i2==0) sigreg="2j0b";
			if(i2==1) sigreg="2j1to2b";
			if(i2==2) sigreg="3to5j0b";
			if(i2==3) sigreg="3to5j1b";
			if(i2==4) sigreg="3to5j2b";
			if(i2==5) sigreg="6j0b";
			if(i2==6) sigreg="6j1b";
			if(i2==7) sigreg="6j2b";
			if(i2==8) sigreg="3b";
			string mapname;
			mapname = "Prediction_"+sigreg;
			if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", 9, 0, 9); hs[mapname]->Sumw2();
			hs[mapname]->GetXaxis()->SetTitle("signal region");
			hs[mapname]->GetXaxis()->SetBinLabel(1,"e,lH_{T}");    hs[mapname]->GetXaxis()->SetBinLabel(2,"e,mH_{T}");    hs[mapname]->GetXaxis()->SetBinLabel(3,"e,hH_{T}");
			hs[mapname]->GetXaxis()->SetBinLabel(4,"#mu,lH_{T}");  hs[mapname]->GetXaxis()->SetBinLabel(5,"#mu,mH_{T}");  hs[mapname]->GetXaxis()->SetBinLabel(6,"#mu,hH_{T}");
			hs[mapname]->GetXaxis()->SetBinLabel(7,"#tau,lH_{T}"); hs[mapname]->GetXaxis()->SetBinLabel(8,"#tau,mH_{T}"); hs[mapname]->GetXaxis()->SetBinLabel(9,"#tau,hH_{T}");
			hs[mapname]->SetMarkerStyle(20);
			hs[mapname]->SetMarkerColor(kBlack);
			hs[mapname]->SetLineColor(kBlack);
			hs[mapname]->SetLineWidth(3);
			mapname = "SimulationTruth_"+sigreg;
			if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", 9, 0, 9); hs[mapname]->Sumw2();
			hs[mapname]->GetXaxis()->SetTitle("signal region");
			hs[mapname]->GetXaxis()->SetBinLabel(1,"e,lH_{T}");    hs[mapname]->GetXaxis()->SetBinLabel(2,"e,mH_{T}");    hs[mapname]->GetXaxis()->SetBinLabel(3,"e,hH_{T}");
			hs[mapname]->GetXaxis()->SetBinLabel(4,"#mu,lH_{T}");  hs[mapname]->GetXaxis()->SetBinLabel(5,"#mu,mH_{T}");  hs[mapname]->GetXaxis()->SetBinLabel(6,"#mu,hH_{T}");
			hs[mapname]->GetXaxis()->SetBinLabel(7,"#tau,lH_{T}"); hs[mapname]->GetXaxis()->SetBinLabel(8,"#tau,mH_{T}"); hs[mapname]->GetXaxis()->SetBinLabel(9,"#tau,hH_{T}");
			hs[mapname]->SetMarkerColor(kBlue);
			hs[mapname]->SetLineColor(kBlue);
			hs[mapname]->SetFillColor(kBlue);
			hs[mapname]->SetLineWidth(0);
			hs[mapname]->SetFillStyle(3002); 
		}
	}
	}

	string oldlep = "dummy";
	string oldsr = "dummy";
	string oldHT = "dummy";
	for(unsigned int n = 0; n<sr.size(); ++n){
		string lep = "Tau";
		string sigreg;
		if(sr[n]==0) sigreg="2j0b";
		if(sr[n]==1) sigreg="2j1to2b";
		if(sr[n]==2) sigreg="3to5j0b";
		if(sr[n]==3) sigreg="3to5j1b";
		if(sr[n]==4) sigreg="3to5j2b";
		if(sr[n]==5) sigreg="6j0b";
		if(sr[n]==6) sigreg="6j1b";
		if(sr[n]==7) sigreg="6j2b";
		if(sr[n]==8) sigreg="3b";
		string htreg;
		if(htr[n]==0 && fMET) htreg = "lowHT";
		if(htr[n]==1 && fHT) htreg = "mediumHT";
		if(htr[n]==2 && fHT) htreg = "highHT";
		int mt2lastbin = 9999;
		if(sr[n]==0) { if(htr[n]==0) mt2lastbin = 750; if(htr[n]==1) mt2lastbin = 1000; if(htr[n]==2) mt2lastbin =  900;}
		if(sr[n]==1) { if(htr[n]==0) mt2lastbin = 700; if(htr[n]==1) mt2lastbin =  700; if(htr[n]==2) mt2lastbin =  350;}
		if(sr[n]==2) { if(htr[n]==0) mt2lastbin = 750; if(htr[n]==1) mt2lastbin = 1000; if(htr[n]==2) mt2lastbin = 1000;}
		if(sr[n]==3) { if(htr[n]==0) mt2lastbin = 700; if(htr[n]==1) mt2lastbin =  900; if(htr[n]==2) mt2lastbin =  550;}
		if(sr[n]==4) { if(htr[n]==0) mt2lastbin = 550; if(htr[n]==1) mt2lastbin =  500; if(htr[n]==2) mt2lastbin =  350;}
		if(sr[n]==5) { if(htr[n]==0) mt2lastbin = 520; if(htr[n]==1) mt2lastbin =  600; if(htr[n]==2) mt2lastbin =  500;}
		if(sr[n]==6) { if(htr[n]==0) mt2lastbin = 450; if(htr[n]==1) mt2lastbin =  500; if(htr[n]==2) mt2lastbin =  500;}
		if(sr[n]==7) { if(htr[n]==0) mt2lastbin = 400; if(htr[n]==1) mt2lastbin =  450; if(htr[n]==2) mt2lastbin =  350;}
		if(sr[n]==8) { if(htr[n]==0) mt2lastbin = 400; if(htr[n]==1) mt2lastbin =  450; if(htr[n]==2) mt2lastbin =  300;}
		if(htreg!=oldHT) oldHT = htreg;
		if(lep!=oldlep) oldlep = lep;
		if(fRebin){
			string mapname = "Prediction_"+lep+"_"+htreg;
			hs[mapname]->SetBinContent(sr[n]+1, datapred[n]);
			hs[mapname]->SetBinError(sr[n]+1, sqrt(pow(datapred_stat_err[n],2)+pow(datapred_syst_err[n],2)) );
			mapname = "SimulationTruth_"+lep+"_"+htreg;
			hs[mapname]->SetBinContent(sr[n]+1, mctruth[n]);
			hs[mapname]->SetBinError(sr[n]+1, mctrutherr[n]);
		} else{
			mt2binning.push_back(MT2low[n]);
			dp.push_back(datapred[n]);
			dpe.push_back(sqrt(pow(datapred_stat_err[n],2)+pow(datapred_syst_err[n],2)) );
			st.push_back(mctruth[n]);
			ste.push_back(mctrutherr[n]);
			if(MT2up[n]==10000.) { 
				mt2binning.push_back(mt2lastbin);
				int Numbins = mt2binning.size()-1;
				int numarray = mt2binning.size();
				double binns[numarray];
				for(int ii = 0; ii<=Numbins; ++ii) binns[ii] = mt2binning[ii];
				string mapname = "Prediction_"+lep+"_"+htreg+"_"+sigreg;
				if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", Numbins, binns); hs[mapname]->Sumw2();
				hs[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); hs[mapname]->GetYaxis()->SetTitle("events");
				hs[mapname]->SetMarkerStyle(20);
				hs[mapname]->SetMarkerColor(kBlack);
				hs[mapname]->SetLineColor(kBlack);
				hs[mapname]->SetLineWidth(3);
				for(int ii = 1; ii<=Numbins; ++ii){
					hs[mapname]->SetBinContent(ii, dp[ii-1]);
					hs[mapname]->SetBinError(ii, dpe[ii-1]);
				}
				mapname = "SimulationTruth_"+lep+"_"+htreg+"_"+sigreg;
				if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", Numbins, binns); hs[mapname]->Sumw2();
				hs[mapname]->GetXaxis()->SetTitle("M_{T2} [GeV]"); hs[mapname]->GetYaxis()->SetTitle("events");
				hs[mapname]->SetMarkerColor(kBlue);
				hs[mapname]->SetLineColor(kBlue);
				hs[mapname]->SetFillColor(kBlue);
				hs[mapname]->SetLineWidth(0);
				hs[mapname]->SetFillStyle(3002); 
				for(int ii = 1; ii<=Numbins; ++ii){
					hs[mapname]->SetBinContent(ii, st[ii-1]);
					hs[mapname]->SetBinError(ii, ste[ii-1]);
				}
				mt2binning.clear();
				dp.clear(); dpe.clear();
				st.clear(); ste.clear();
			}
		}
	}
	//now store all the shit
   	TFile *fpredictionfile = new TFile(filename.Data(),"RECREATE");
	fpredictionfile->cd();
	for(map<string,TH1D*>::iterator h=hs.begin(); h!=hs.end();++h){
		h->second->Write();
	}
	fpredictionfile->Close();
	*fLogStream << "Saved prediction/mc truth histograms in " << filename.Data() << endl;


	if(PlotPrediction){
		//format histogram
		for(map<string,TH1D*>::iterator h=hs.begin(); h!=hs.end();++h){
			TString helperstring = h->first;
			h->second->SetStats(0);
			h->second->SetLineStyle(0);
			h->second->GetXaxis()->SetLabelOffset(0.007);
			h->second->GetXaxis()->SetLabelSize(0.05);
			h->second->GetXaxis()->SetTitleSize(0.06);
			h->second->GetXaxis()->SetTitleOffset(0.9);
			h->second->GetXaxis()->SetTitleFont(42);
			h->second->GetYaxis()->SetTitle("Events");
			h->second->GetYaxis()->SetLabelFont(42);
			h->second->GetYaxis()->SetLabelOffset(0.007);
			h->second->GetYaxis()->SetLabelSize(0.05);
			h->second->GetYaxis()->SetTitleSize(0.06);
			h->second->GetYaxis()->SetTitleOffset(1.25);
			h->second->GetYaxis()->SetTitleFont(42);
			h->second->GetZaxis()->SetLabelFont(42);
			h->second->GetZaxis()->SetLabelOffset(0.007);
			h->second->GetZaxis()->SetLabelSize(0.05);
			h->second->GetZaxis()->SetTitleSize(0.06);
			h->second->GetZaxis()->SetTitleFont(42);
			h->second->GetXaxis()->SetLabelFont(42);
//			if(h->first[0]==(string)"S"){//simulation truth gets plotted first
			if(helperstring.Contains("Simulation")){//simulation truth gets plotted first
			h->second->SetMarkerStyle(1);
			h->second->SetMarkerColor(kBlue);
			h->second->SetLineColor(kBlue);
			h->second->SetFillColor(kBlue);
			h->second->SetLineWidth(0);
			h->second->SetFillStyle(3002); 
			if(helperstring.Contains("b") && fRebin){
				h->second->GetXaxis()->SetTitle("signal region");
				h->second->GetXaxis()->SetBinLabel(1,"e,lH_{T}");
				h->second->GetXaxis()->SetBinLabel(2,"e,mH_{T}");
				h->second->GetXaxis()->SetBinLabel(3,"e,hH_{T}");
				h->second->GetXaxis()->SetBinLabel(4,"#mu,lH_{T}");
				h->second->GetXaxis()->SetBinLabel(5,"#mu,mH_{T}");
				h->second->GetXaxis()->SetBinLabel(6,"#mu,hH_{T}");
				h->second->GetXaxis()->SetBinLabel(7,"#tau,lH_{T}");
				h->second->GetXaxis()->SetBinLabel(8,"#tau,mH_{T}");
				h->second->GetXaxis()->SetBinLabel(9,"#tau,hH_{T}");
		//	} else if(h->first[h->first.size()]==(string)"b" && !fRebin){
			} else if(helperstring.Contains("b") && !fRebin){
				h->second->GetXaxis()->SetTitle("M_{T2} [GeV]");
			} else{//fRebin and last letter is a "T"
				h->second->GetXaxis()->SetTitle("signal region");
				h->second->GetXaxis()->SetBinLabel(1,"2j,0b");
				h->second->GetXaxis()->SetBinLabel(2,"2j,#geq1b");
				h->second->GetXaxis()->SetBinLabel(3,"3-5j,0b");
				h->second->GetXaxis()->SetBinLabel(4,"3-5j,1b");
				h->second->GetXaxis()->SetBinLabel(5,"3-5j,2b");
				h->second->GetXaxis()->SetBinLabel(6,"#geq6j,0b");
				h->second->GetXaxis()->SetBinLabel(7,"#geq6j,1b");
				h->second->GetXaxis()->SetBinLabel(8,"#geq6j,2b");
				h->second->GetXaxis()->SetBinLabel(9,"#geq3j,#geq3b");
			}
			} else{
			h->second->SetMarkerStyle(20);
			h->second->SetMarkerColor(kBlack);
			h->second->SetLineColor(kBlack);
			h->second->SetLineWidth(3);
			}
		}
		TCanvas *c1 = new TCanvas("c1", "",60,22,600,600);
		gStyle->SetOptFit(1);
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		c1->SetFillColor(0);
		c1->SetBorderMode(0);
		c1->SetBorderSize(2);
		c1->SetTickx(1);
		c1->SetTicky(1);
		c1->SetLeftMargin(0.18);
		c1->SetRightMargin(0.05);
		c1->SetTopMargin(0.07);
		c1->SetBottomMargin(0.15);
		c1->SetFrameFillStyle(0);
		c1->SetFrameBorderMode(0);
		c1->SetFrameFillStyle(0);
		c1->SetFrameBorderMode(0);
		c1->cd();

		TLatex *   tex = new TLatex(0.328859,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV");
		tex->SetNDC();
		tex->SetTextFont(42);
		tex->SetTextSize(0.04181185);
		tex->SetLineWidth(2);
		TLatex toplep;
		toplep.SetNDC();
		toplep.SetTextAlign(31);
		toplep.SetTextFont(42);
		toplep.SetTextSize(0.04181185);
		toplep.SetLineWidth(2);
		TLatex ht;
		ht.SetNDC();
		ht.SetTextAlign(31);
		ht.SetTextFont(42);
		ht.SetTextSize(0.04181185);
		ht.SetLineWidth(2);

		TLegend *leg = new TLegend(0.6252416,0.7657343,0.824906,0.9003497,NULL,"brNDC");
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04181185);
		leg->SetTextFont(42);
		leg->SetLineColor(1);
		leg->SetLineStyle(1);
		leg->SetLineWidth(2);
		leg->SetFillColor(0);
		leg->SetFillStyle(1001);
   		TLegendEntry *entry=leg->AddEntry("NULL","simulation truth","f");
		entry->SetMarkerColor(kBlue);
		entry->SetLineColor(kBlue);
		entry->SetFillColor(kBlue);
		entry->SetLineWidth(0);
		entry->SetFillStyle(3002); 
   		TLegendEntry *entry2=leg->AddEntry("NULL","data prediction","p");
		entry2->SetLineColor(1);
		entry2->SetLineStyle(1);
		entry2->SetLineWidth(2);
		entry2->SetMarkerColor(1);
		entry2->SetMarkerStyle(20);
		entry2->SetMarkerSize(1);

		string hsp, hss;
		TString texttoplep, textht;
		string outname;
		double max(0.); double min(99999.);
		double maxp(0.); double minp(99999.);
		double maxs(0.); double mins(99999.);
		for(int i3 = 0; i3<HTbinsize;        ++i3){
			string htreg;
			if(i3==0 && fMET) { htreg = "lowHT"; textht = "low H_{T}"; }
			if(i3==0 && fHT)  { htreg = "mediumHT"; textht = "medium H_{T}"; }
			if(i3==1 && fHT)  { htreg = "highHT"; textht = "high H_{T}"; }
			if(fRebin){
				hsp = "Prediction_Tau_"+htreg;
				hss = "SimulationTruth_Tau_"+htreg;
				texttoplep = "1 tau";
				min = 0.;
				maxp = hs[hsp]->GetBinContent(hs[hsp]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hsp]->GetMaximumBin());
				maxs = hs[hss]->GetBinContent(hs[hss]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hss]->GetMaximumBin());
				max  = (maxp>maxs)?maxp:maxs;
				max = 1.5*max;
				hs[hss]->SetMaximum(max);
				c1->Clear();
				c1->cd();
				hs[hss]->SetFillStyle(3013);
				hs[hss]->Draw("E2");
				hs[hsp]->Draw("sameE");
				tex->Draw();
				leg->Draw();
				toplep.DrawLatex(0.6,0.8548951,texttoplep.Data());
				ht.DrawLatex(0.6,0.791958,textht.Data());
				outname = plotdirectory+hsp + ".eps";
				c1->SaveAs(outname.c_str());
			}
			for(int i2 = 0; i2<signalregionsize; ++i2){
				string sigreg; TString regsig;
				if(i2==0) { sigreg="2j0b";    regsig="2j, 0b"; }
				if(i2==1) { sigreg="2j1to2b"; regsig="2j, #geq1b"; }
				if(i2==2) { sigreg="3to5j0b"; regsig="3-5j, 0b"; }
				if(i2==3) { sigreg="3to5j1b"; regsig="3-5j, 1b"; }
				if(i2==4) { sigreg="3to5j2b"; regsig="3-5j, 2b"; }
				if(i2==5) { sigreg="6j0b";    regsig="#geq6j, 0b"; }
				if(i2==6) { sigreg="6j1b";    regsig="#geq6j, 1b"; }
				if(i2==7) { sigreg="6j2b";    regsig="#geq6j, 2b"; }
				if(i2==8) { sigreg="3b";      regsig="#geq3j, #geq3b"; }
				if(!fRebin){
					hsp = "Prediction_Tau_"+htreg+"_"+sigreg;
					hss = "SimulationTruth_Tau_"+htreg+"_"+sigreg;
					texttoplep = "1 tau, "+ regsig;
					min = 0.;
					maxp = hs[hsp]->GetBinContent(hs[hsp]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hsp]->GetMaximumBin());
					maxs = hs[hss]->GetBinContent(hs[hss]->GetMaximumBin())+hs[hsp]->GetBinError(hs[hss]->GetMaximumBin());
					max  = (maxp>maxs)?maxp:maxs;
					max = 1.25*max;
					hs[hss]->SetMaximum(max);
					c1->Clear();
					c1->cd();
					hs[hss]->SetFillStyle(3013);
					hs[hss]->Draw("E2");
					hs[hsp]->Draw("sameE");
					tex->Draw();
					leg->Draw();
					toplep.DrawLatex(0.6,0.8548951,texttoplep.Data());
					ht.DrawLatex(0.6,0.791958,textht.Data());
					outname = plotdirectory+hsp + ".eps";
					c1->SaveAs(outname.c_str());
				}
			}
		}
	}
}


//this function just makes the one lepton yield table (usually with MT cut applied for one lepton selection, unless flags are set differently)
//it contains no comments as this function is from style very similar to SummaryTable(...)
void YieldTable(vector<int> sr, vector<int> htr, vector<double> MT2low, vector<double> MT2up, vector<double> numQCD, vector<double> numZ, vector<double> numW, vector<double> numT, vector<double> numOther, vector<double> numMC, vector<double> BGerr, vector<double> numData, vector<double> numBGLL){

	*fLogStream << "*********************************************************************" << endl;
	*fLogStream << "\%BEGINLATEX\%"              << endl;
	*fLogStream << "\\begin{table}"              << endl
	     << "\\begin{center}"             << endl
	     << "\\small"                     << endl;
	if(fPlotLLTab){
              *fLogStream << "\\begin{tabular}{lcccccccc}" << endl
	                  << "\\hline\\hline"              << endl;
	} else {
              *fLogStream << "\\begin{tabular}{lccccccc}"  << endl
	                  << "\\hline\\hline"              << endl;
	}

	if(!fRebin)    *fLogStream << "$M_{T2}$ (GeV)";
	else           *fLogStream << "signal region       ";
	               *fLogStream  << " & $" << "N^{QCD}" << "$ & $" << "N^{Z}" << "$  & $" << "N^{W}" << "$  & $" << "N^{Top}" << "$  & $" << "N^{Other}" << "$ &";
	if(fPlotLLTab) *fLogStream << " $" << "N_{LL}^{MC}" << "$ &";
	               *fLogStream  << "            $"  << "N^{MC}" << "$      & $" << "N^{data}" << "$ ";
	               *fLogStream << "\\\\" << endl << "\\hline\\hline"             << endl;
	string oldlep = "dummy";
	string oldsr = "dummy";
	string oldHT = "dummy";
	for(unsigned int n = 0; n<sr.size(); ++n){
		string sigreg;
		if(sr[n]==0) sigreg="2j, 0b";
		if(sr[n]==1) sigreg="2j, $\\ge 1$b";
		if(sr[n]==2) sigreg="$3-5$j, 0b";
		if(sr[n]==3) sigreg="$3-5$j, 1b";
		if(sr[n]==4) sigreg="$3-5$j, 2b";
		if(sr[n]==5) sigreg="$\\ge 6$j, 0b";
		if(sr[n]==6) sigreg="$\\ge 6$j, 1b";
		if(sr[n]==7) sigreg="$\\ge 6$j, 2b";
		if(sr[n]==8) sigreg="$\\ge 3$j, $\\ge 3$b";
		string htreg;
		if(htr[n]==0 && fMET) htreg = "450 GeV $\\leq H_{T} < 750$ GeV";
		if(htr[n]==1 && fHT) htreg = "750 GeV $\\leq H_{T} < 1200$ GeV";
		if(htr[n]==2 && fHT) htreg = "$H_{T}\\ge 1200$ GeV";
		if(htreg!=oldHT){
			*fLogStream << " \\hline " << endl << "\\multicolumn{7}{l}{" <<  htreg << "} \\\\ " << endl << "\\hline" << endl;
			oldHT = htreg;
		}
		if(!fRebin){
			if(sigreg!=oldsr){
				*fLogStream << " \\hline  " << endl << sigreg << "\\\\" << endl;
				oldsr = sigreg;
			}
			if(MT2up[n]==10000.) *fLogStream << "$" << int(MT2low[n]) << "-" << "\\infty$" << " " << setw(4) << " & ";
			else                 *fLogStream << "$" << int(MT2low[n]) << "-" << int(MT2up[n]) << "$" << " " << setw(7) << " & ";
		}
		else                         *fLogStream << " " << setw(18) << sigreg << " & ";
		               *fLogStream << fixed << setprecision(2) << " " << setw(9) << numQCD[n] << " & " << " " << setw(7) << numZ[n] << " & " << " " << setw(7) << numW[n] << " &  " << " " << setw(8) << numT[n] << " & " << " " << setw(10) << numOther[n] << " & ";
		if(fPlotLLTab) *fLogStream << fixed << setprecision(2) << " " << setw(12) << numBGLL[n] << " & ";
		               *fLogStream << fixed << setprecision(2) << " " << setw(10) << numMC[n] << "$\\pm" << " " << setw(7) << BGerr[n] << "$ & ";
		               *fLogStream << " " << setw(10) << int(numData[n]) << " \\\\" << endl;
	}
	*fLogStream << "\\hline\\hline"                                                                                                << endl
		<< "\\end{tabular}"                                                                                                << endl
		<< "\\end{center}"                                                                                                 << endl
		<< "\\end{table}"                                                                                                  << endl
		<< "\%ENDLATEX\%"                                                                                                  << endl
		<< endl;
	*fLogStream << endl << endl;
}

//this function makes a compact card of signal region vs MCtruth +/- err vs. pred. +/- err
//it contains no comments as this function is from style very similar to SummaryTable(...)
void PredictionCard(vector<int> sr, vector<int> htr, vector<int> MT2bin, vector<double> mctruth, vector<double> mctrutherr, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, vector<double> LLeff, vector<double> LLefferr, vector<double> RLL, vector<double> RLLerr, vector<double> MTeff, vector<double> MTefferr){
	
	*fLogStream << "*********************************************************************" << endl;
	*fLogStream << " Name            " << " " << " SR " << " HTbin " << " MT2bin   " << " MCpred   " << " DataDrivenPred " << " DataDrivenPredError " << " ScaleFactor " << " ScaleFactorError " << endl;
	for(unsigned int n = 0; n<sr.size(); ++n){


		if(LLeff[n]!=0){
		*fLogStream << "LostLepton Tau  " << " " << setw(3) << sr[n] << " " << setw(5) << htr[n] << " " << setw(6) << MT2bin[n]-1 << " " << setw(14) << fixed << setprecision(4) << mctruth[n] <<  " " << setw(12) << datapred[n] << " " << setw(17) << sqrt( pow(datapred_stat_err[n],2) + pow(datapred_syst_err[n],2) );
		*fLogStream <<  " " << setw(22) << fixed << setprecision(4) << RLL[n] << " " << setw(14) << RLLerr[n] << endl;
		}
		else{
		*fLogStream << "LostLepton Tau  " << " " << setw(3) << sr[n] << " " << setw(5) << htr[n] << " " << setw(6) << MT2bin[n]-1 << " " << setw(14) << fixed << setprecision(4) << mctruth[n] << " " << setw(12) << "-9" << " " << setw(17) << "-9";
            	*fLogStream << " " << setw(22) << "-9" << " " << setw(14) << "-9" << endl;
		}
	}
	*fLogStream << endl << endl;

}

//this function makes a compact card of signal region vs MCtruth +/- err vs. pred. +/- err, and also prints out (from MC truth) the fraction of expected lost lepton from W,Top,rest(dibosons?)
//it contains no comments as this function is from style very similar to SummaryTable(...)
void PredictionCardSplitted(vector<int> sr, vector<int> htr, vector<int> MT2bin, vector<double> mctruth, vector<double> mctrutherr, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, vector<double> LLeff, vector<double> LLefferr, vector<double> RLL, vector<double> RLLerr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruthW, vector<double> mctruthTT, vector<double> numW, vector<double> numTT, vector<double> numMC){
	
	*fLogStream << "*********************************************************************" << endl;
	*fLogStream << " Name            " << " " << " SR " << " HTbin " << " MT2bin   " << " MCpred   " << " DataDrivenPred " << " DataDrivenPredError " << " ScaleFactor " << " ScaleFactorError      " << " W:TTbar:Other(truth)          " << " W:TTbar:Other(lepton yield) " << endl;

	for(unsigned int n = 0; n<sr.size(); ++n){

		if(LLeff[n]!=0){
		*fLogStream << "LostLepton Tau  " << " " << setw(3) << sr[n] << " " << setw(5) << htr[n] << " " << setw(6) << MT2bin[n]-1 << " " << setw(14) << fixed << setprecision(4) << mctruth[n] << " " << setw(12) << datapred[n] << " " << setw(17) << sqrt( pow(datapred_stat_err[n],2) + pow(datapred_syst_err[n],2) );
		*fLogStream <<  " " << setw(22) << fixed << setprecision(4) << RLL[n] << " " << setw(14) << RLLerr[n];
		*fLogStream << " " << setw(23) << fixed << setprecision(4) << mctruthW[n]/mctruth[n]<< " : " <<mctruthTT[n]/mctruth[n]<< " : " << (mctruth[n]-mctruthW[n]-mctruthTT[n])/mctruth[n] << " " << setw(15) << numW[n]/numMC[n]<< " : " <<numTT[n]/numMC[n]<< " : " <<(numMC[n]-numW[n]-numTT[n])/numMC[n] << endl;
		}
		else{
		*fLogStream << "LostLepton Tau  " << " " << setw(3) << sr[n] << " " << setw(5) << htr[n] << " " << setw(6) << MT2bin[n]-1 << " " << setw(14) << fixed << setprecision(4) << mctruth[n] << " " << setw(12) << "-9" << " " << setw(17) << "-9";
            	*fLogStream << " " << setw(22) << "-9" << " " << setw(14) << "-9";
		*fLogStream << " " << setw(23) << fixed << setprecision(4) << mctruthW[n]/mctruth[n]<<" : "<<mctruthTT[n]/mctruth[n]<<" : "<<(mctruth[n]-mctruthW[n]-mctruthTT[n])/mctruth[n] << " " << setw(26) << "- : - : -" << endl;
		}
	}
	*fLogStream << endl << endl;

}