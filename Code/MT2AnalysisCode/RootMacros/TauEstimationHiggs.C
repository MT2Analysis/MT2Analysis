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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path of MT2tree.hh
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"

//run via root -l -b -q TauEstimationHiggs.C++

using namespace std;

void load(const char* filename = "samples/datasamples/samples_2141_dataonly.dat");
void TauEstimationHiggs();
TH1D* RebinThisHistogram(TH1D* histogram, int nrebinnedbins, double* rebinnedbins);
TEfficiency* RebinThisEfficiency(TEfficiency* tefficiency, int nrebinnedbins, double* rebinnedbins);
void GetTauEstimate(map<string, TH1D*> histograms, map<string, TEfficiency*> tefficiencies, double rel_sys_uncert,  double rel_sys_uncert_bg, int version, Bool_t PrintSummaryTable, Bool_t PrintYieldTable, Bool_t PrintPredictionCard, Bool_t makeFullPrintout=false, Bool_t SavePrediction=true, Bool_t PlotPrediction=false);
void PredictionCard(vector<int> sr, vector<int> htr, vector<int> MT2bin, vector<double> mctruth, vector<double> mctrutherr, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, vector<double> LLeff, vector<double> LLefferr, vector<double> RLL, vector<double> RLLerr, vector<double> MTeff, vector<double> MTefferr);
void PredictionCardSplitted(vector<int> sr, vector<int> htr, vector<int> MT2bin, vector<double> mctruth, vector<double> mctrutherr, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, vector<double> LLeff, vector<double> LLefferr, vector<double> RLL, vector<double> RLLerr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruthW, vector<double> mctruthTT, vector<double> numW, vector<double> numTT, vector<double> numMC);//neew
void YieldTable(vector<int> sr, vector<int> htr, vector<double> MT2low, vector<double> MT2up, vector<double> numQCD, vector<double> numZ, vector<double> numW, vector<double> numT, vector<double> numOther, vector<double> numMC, vector<double> BGerr, vector<double> numData, vector<double> numBGLL);
void SummaryTable(int version, vector<int> sr, vector<int> htr, vector<double> MT2low, vector<double> MT2up, vector<double> numtrueW, vector<double> numtrueW_bg, vector<double> numData, vector<double> numBG, vector<double> numBGLL, vector<double> LLeff, vector<double> LLefferr, vector<double> RLL, vector<double> RLLerr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruth, vector<double> mctrutherr, vector<double> MCpred, vector<double> MCpred_stat_err, vector<double> MCpred_syst_err, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err);
void PredictionFile(int version, vector<int> sr, vector<int> htr, vector<double> MT2low, vector<double> MT2up, vector<double> numtrueW, vector<double> numtrueW_bg, vector<double> numData, vector<double> numBG, vector<double> numBGLL, vector<double> LLeff, vector<double> LLefferr, vector<double> RLL, vector<double> RLLerr, vector<double> MTeff, vector<double> MTefferr, vector<double> mctruth, vector<double> mctrutherr, vector<double> MCpred, vector<double> MCpred_stat_err, vector<double> MCpred_syst_err, vector<double> datapred, vector<double> datapred_stat_err, vector<double> datapred_syst_err, Bool_t PlotPrediction);


//standard struct combining MT2trees with important information like cross section
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

Bool_t fMET         = false;//use MET triggers (i.e. used for low HT region)
Bool_t fHT          = true; //use HT triggers (i.e. used for medium+high HT region), will be set false automatically if fMET==true
Bool_t fFast        = true;//includes a lower MT2 cut (MT2>100 GeV), default = true, as low MT2 region not part of any signal region
Bool_t fRebin       = true;//do not bin-by-bin estimation but one estimate per HT/topological region (this means not along MT2)

Bool_t fISRreweight =true; //reweight MC according to SUSY's 'ISR recipe' - influence on estimate is minimal, but MC truth changes
                          //keep that flag true - then you will save both the ISR reweighted and non-reweighted style, later use TauEstimationFromRootfile to get the number for your chosen case

std::ostringstream* fLogStream     = 0;
Bool_t  fWriteToFile               = false; // writes couts to a file
Bool_t  fAppend                    = true;  // append at end of file (if existing), otherwise delete previous content - needs fWriteToFile = true

Bool_t fData         = true;//set false if don't want to run over data - only if doing checks
Bool_t fQCD          = true;//set false if don't want to run over QCD samples - only if doing checks and want to save time
Bool_t fSusy         = false;//set false if don't want to run over Susy samples - default is false

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

//LostTau estimation for the MT2 Higgs selection - for detailed comments see TauEstimation.C (which does LostTau for inclusive MT2 analysis)
//this function does the LostLeptonEstimate for (hadronic) taus, electrons and muons are done with LostLeptonEstimate.C, not so detailed commented as LostLeptonEstimate.C
//Basically only some histograms are filled and later (in another function of this macro) all relevant quantities are extracted out of all these histograms, depending on the flags you set before.
void TauEstimationHiggs(){

	// logStream
	fLogStream = new std::ostringstream();

	if(fMET==true) fHT = false;
	if(fMET==false && fHT==false) fHT = true;
	if(fbTagReweight==false) fbTagError = false;
	if(fIncludeTop==false) fIncludeSingleTop = false;
	if(fIncludeTop==false) fTopOnly = false;
  	gROOT->ProcessLine(".x SetStyle_PRD.C");

	TString  outputdir = "Filtered/TauEstimation/Higgs/";
	if(fMET) outputdir = outputdir + "MET/";
	if(fHT)  outputdir = outputdir + "HT/";
	if(fUseDoubleTau) outputdir += "used_dTau/";
	if(fDoubleTau_BG) outputdir += "veto_dTau/";
    	Util::MakeOutputDir(outputdir);

	TString outputname = "TauEstimationHistograms_limitMbb.root";

	TString  samples = "samples/dummy_filter.dat";//only dummy
	if(fMET) samples = "samples/samples_MET_filter.dat";
	if(fHT)  samples = "samples/samples_HT_filter.dat";

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
	//the next histograms are not used in prediction but important for checks, get tau efficiency, (jet-->)tau fake rate, tau purity
	histonames.push_back("TrueRecoTauEvents");		//all reco taus w/o MT cut matched to gentaus     //note, this is not on event basis
	histonames.push_back("FakeRecoTauEvents");		//all reco taus w/o MT cut not matched to gentaus //note, this is not on eventbasis
	histonames.push_back("AllRecoTauEvents");		//all reco taus w/o MT cut                        //note, this is not on eventbasis
	histonames.push_back("AllRecoJetEvents");		//all reco jets w/o MT cut                        //note, this is not on eventbasis
	histonames.push_back("AllGenTauAccEvents");		//all gen had taus within acceptance              //note, this is not on eventbasis
	//below are histograms used for BTV SF uncertainty, first SFup
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
	//now SFdown
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
	vector<string> teffnames; teffnames.clear();//filled later

	for(int i1 = 0; i1<sampletypesize;   ++i1){
		string hs = string("_")+ sample_type[i1];
		string mapname;
		for(unsigned int i0 = 0; i0<histonames.size(); ++i0){
			mapname = histonames[i0] + hs;
			if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", 33,5,500);//histogram has larger range than is defined in event selection
			mapname = "noISR_" + histonames[i0] + hs;
			if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", 33,5,500);
		}
	}

	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Sumw2();}

	std::ostringstream cutStream;
	std::ostringstream cutStreamBase;
	cutStream << " " 
		<< "NEles+NMuons  ==0"                      << "&&"
		<< "misc.Jet0Pass ==1"                      << "&&"
		<< "misc.Jet1Pass ==1"                      << "&&"
		<< "misc.Vectorsumpt < 70";
		cutStream  << "&& misc.MinMetJetDPhi4Pt40 >0.3";
	if(fMET) cutStream << "&&misc.MET>200&&misc.HT<=750&&misc.HT>=450&&misc.MT2>200";
	if(fHT ) cutStream << "&&misc.HT>750&&misc.MET>30&&misc.MT2>125";
	if(fFast)cutStream << "&&misc.MT2>=100";//lowest border in MT2
	cutStream <<"&&NJetsIDLoose40 >= 4"<<"&&"
    		  <<"NBJetsCSVM>=2";
	cutStream << "&&GetSelBBMinv()>=20&&GetSelBBMinv()<=200";

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
	
	std::ostringstream triggerStream;
	if(fMET){
	triggerStream << "( ( ("
			<< "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
			<< "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )"
			<< "||("
			<< "trigger.HLT_PFHT350_PFMET100_v3==1 || trigger.HLT_PFHT350_PFMET100_v4==1 || trigger.HLT_PFHT350_PFMET100_v5==1 || "
			<< "trigger.HLT_PFHT350_PFMET100_v6==1 || trigger.HLT_PFHT350_PFMET100_v7==1 || trigger.HLT_PFNoPUHT350_PFMET100_v1==1 || "
			<< "trigger.HLT_PFNoPUHT350_PFMET100_v3==1 || trigger.HLT_PFNoPUHT350_PFMET100_v4==1 ) ) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	}
	if(fHT){
	triggerStream << "( ("
			<< "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || "
			<< "trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
			<< "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	}
	TString trigger = triggerStream.str().c_str();
	
	TString cuts = cutStream.str().c_str();
	TString basecuts = cutStreamBase.str().c_str();

	load(samples.Data());

   	for(size_t i = 0; i < fSamples.size(); ++i){
        
   	    if(fData==false && fSamples[i].type=="data") continue;
	    if(fSusy==false && fSamples[i].type=="susy") continue;
	    if(fQCD==false  && fSamples[i].sname=="QCD") continue;
        
		//get sample name
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
		//global event weight
	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "TauEstimate: looping over " << fSamples[i].name << " added in " << sampletype << endl;
	    if(fVerbose>2) cout << "             sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts + "&&" + basecuts;
        
	    if( fSamples[i].type=="data") myCuts += " && " + trigger; //cuts to be aplied only on data
        
   	    cout << "Cuts for Flow: " << myCuts << endl;
   	    fSamples[i].tree->Draw(">>selList", myCuts);
        
	    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
	    fSamples[i].tree->SetEventList(myEvtList);
	    int counter=0;
	    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
	    if(myEvtList->GetSize()==0) continue;
	//run over all selected events
        while(myEvtList->GetEntry(counter++) !=-1){	
		int jentry = myEvtList->GetEntry(counter-1);
            
		nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
		fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
		
		if ( fVerbose>2 && counter % 50000 == 0 )  cout << "+++ Proccessing event " << counter << endl;
		
		Double_t weight = sample_weight;
		if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 
		Double_t ISRweight(1.); Double_t weightnoISR = weight;
		//get 'ISR weight'
		if(fISRreweight && !fMT2tree->misc.isData){
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
		}
		double btagSF(1.), btagSFerr(0.);
		if(!fMT2tree->misc.isData && fbTagReweight){
			btagSF = fMT2tree->SFWeight.BTagCSV40ge2; btagSFerr = fMT2tree->SFWeight.BTagCSV40ge2Error;
			weight = weight * btagSF;
			weightnoISR = weightnoISR * btagSF;
		}

		string hh = string("_") + sampletype;

		double Mbb = fMT2tree->GetSelBBMinv();
		if(Mbb<0) continue;

		Bool_t recoedtau     = false;// exact 1 reconstructed tau including MT cut on selected tau
		Bool_t recoedtaunomt = false;// exact 1 reconstructed tau without MT cut on selected tau
		Bool_t norecotau     = false;//       0 reconstructed taus
		int NumRecoTaus = 0;
		int NumRecoTausMT = 0;
			for(int n = 0; n<fMT2tree->NTaus; ++n){
				if(!(fMT2tree->tau[n].isLooseID3Hits)) continue;
	//			++NumRecoTaus;
				if(fMT2tree->tau[n].MT>100) continue;
				++NumRecoTausMT;
			}
		NumRecoTaus = fMT2tree->NTausIDLoose3Hits;
		if(NumRecoTausMT==1) recoedtau     = true;
		if(NumRecoTaus  ==1) recoedtaunomt = true;
		if(NumRecoTaus  ==0) norecotau     = true;

		if(recoedtau     ) histos[(string)"TauEvents"        + hh]->Fill(Mbb, weight);
		if(recoedtaunomt ) histos[(string)"TauEvents_noMT"   + hh]->Fill(Mbb, weight);
		if(recoedtau     ) histos[(string)"noISR_TauEvents"        + hh]->Fill(Mbb, weightnoISR);
		if(recoedtaunomt ) histos[(string)"noISR_TauEvents_noMT"   + hh]->Fill(Mbb, weightnoISR);
		if(fbTagError && !fMT2tree->misc.isData && btagSF!=0){
			if(recoedtau     ) histos[(string)"TauEvents_SFup"       + hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
			if(recoedtaunomt ) histos[(string)"TauEvents_noMT_SFup"  + hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
			if(recoedtau     ) histos[(string)"TauEvents_SFdown"     + hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
			if(recoedtaunomt ) histos[(string)"TauEvents_noMT_SFdown"+ hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
			if(recoedtau     ) histos[(string)"noISR_TauEvents_SFup"       + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
			if(recoedtaunomt ) histos[(string)"noISR_TauEvents_noMT_SFup"  + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
			if(recoedtau     ) histos[(string)"noISR_TauEvents_SFdown"     + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
			if(recoedtaunomt ) histos[(string)"noISR_TauEvents_noMT_SFdown"+ hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
		}
		if(/*sampletype=="QCD" || sampletype=="Other" ||*/ sampletype=="data" || fMT2tree->misc.isData) continue;

		int genele = fMT2tree->GenNumLeptFromW(11,  0, 1000, true);
		int genmuo = fMT2tree->GenNumLeptFromW(13,  0, 1000, true);
		int genlep = fMT2tree->GenNumLeptFromW(1113,0,1000,  true);
		int genelenotau = fMT2tree->GenNumLeptFromW(11,  0, 1000, false);
		int genmuonotau = fMT2tree->GenNumLeptFromW(13,  0, 1000, false);
		int genlepnotau = fMT2tree->GenNumLeptFromW(1113,0,1000,  false);
		int geneletau = genele-genelenotau;
		int genmuotau = genmuo-genmuonotau;
		int genleptau = genlep-genlepnotau;//electron/muons from W-->nutau-->nununulep, i.e. from tau decays
		int gentau = fMT2tree->GenNumLeptFromW(16,  0,1000, false);//using tau neutrinos
		int genhadtau = gentau-genleptau;
	bool leptfoundnomt  = false;
	bool leptfound      = false;
	bool eventgood      = false;
	bool acceptance     = false;
	bool leptfoundnomt_noLL  = false;
	bool leptfound_noLL      = false;
	bool eventgood_noLL      = false;
	bool acceptance_noLL     = false;

	if(fUseDoubleTau){
		if((genhadtau==1||genhadtau==2))                    eventgood      = true;
		if((genhadtau==1||genhadtau==2) &&genele+genmuo==0) eventgood_noLL = true;
	} else {
		if(genhadtau==1)                                    eventgood      = true;
		if(genhadtau==1 &&genele+genmuo==0)                 eventgood_noLL = true;//note: vetoing leptons by NEles+NMuons==0
	}
	if(recoedtau && eventgood)          leptfound          = true;
	if(recoedtaunomt && eventgood)      leptfoundnomt      = true;
	if(recoedtau && eventgood_noLL)     leptfound_noLL     = true;
	if(recoedtaunomt && eventgood_noLL) leptfoundnomt_noLL = true;

	Bool_t LLfound = false;
	if(genlep-fMT2tree->NEles-fMT2tree->NMuons>0) LLfound   = true;//more genele+genmuo than recoele+recomuo --> Lost electron+muon

	//do fancy via genlept loop to get true genlept.Pt() (corrected for momentum from neutrino)
	int gentauwithinacc = 0;
	int recotaumatchedtrue = 0;
	int recotaumatchedfake = 0;
	for(int j = 0; j< fMT2tree->NGenLepts; ++j){
		if(abs(fMT2tree->genlept[j].ID  )!=16) continue;
		if(abs(fMT2tree->genlept[j].MID )!=15) continue;
		if(abs(fMT2tree->genlept[j].GMID)!=24) continue;
		TLorentzVector gtau = fMT2tree->genlept[j].Mlv;//gentau
		if(gtau.Pt()<0.001) continue;
		bool gtauishad = true;//assume all taus are hadronic taus
		for(int jj = 0; jj< fMT2tree->NGenLepts; ++jj){//filter leptonc taus --> e/mu
			if(abs(fMT2tree->genlept[jj].ID  )!=11 && abs(fMT2tree->genlept[jj].ID  )!=13) continue;
			if(abs(fMT2tree->genlept[jj].MID )!=15)  continue;
			if(abs(fMT2tree->genlept[jj].GMID)!=24)  continue;
			if(fMT2tree->genlept[j].MID!= fMT2tree->genlept[jj].MID) continue;//match 'charge'
			TLorentzVector gtau2 = fMT2tree->genlept[jj].Mlv;//match mother-lv
			if(fabs(gtau2.Pt()  - gtau.Pt()) >0.001) continue;
			if(fabs(gtau2.Eta() - gtau.Eta())>0.001) continue;
			if(fabs(gtau2.Phi() - gtau.Phi())>0.001) continue;
			gtauishad = false;//gtau matched to genele/genmuo from tau decay
			break;//it is leptonic
		}
		if(!gtauishad) continue;//continue if tau is not hadronic
		TLorentzVector visiblegtau = fMT2tree->genlept[j].Mlv - fMT2tree->genlept[j].lv; //subtract neutrino from gentau
		//get match to reco
		int matched = false;
		for(int n = 0; n<fMT2tree->NTaus; ++n){
			if(fMT2tree->tau[n].lv.Pt()<20.)       continue;
			if(!(fMT2tree->tau[n].isLooseID3Hits)) continue;
			if(fMT2tree->tau[n].lv.DeltaR(visiblegtau) < 0.5) { matched = true; break;}//found match)
		}
		if(matched) ++recotaumatchedtrue;
		else        ++recotaumatchedfake;
		//get acceptance of visible gen tau
		if(visiblegtau.Pt()<20)         continue;
		if(fabs(visiblegtau.Eta())>2.4) continue;
		acceptance = true;
		++gentauwithinacc;
		if(eventgood_noLL) acceptance_noLL = true;
	}

	if(recoedtau     && genhadtau==0)            histos[(string)"TauEventsBG"           + hh]->Fill(Mbb, weight);//bg
	if(recoedtaunomt && genhadtau==0)            histos[(string)"TauEventsBG_noMT"      + hh]->Fill(Mbb, weight);//bg
	if(recoedtau     && genhadtau==0 && LLfound) histos[(string)"TauEventsBG_LL"        + hh]->Fill(Mbb, weight);//bg
	if(recoedtaunomt && genhadtau==0 && LLfound) histos[(string)"TauEventsBG_noMT_LL"   + hh]->Fill(Mbb, weight);//bg
	if(recoedtau     && LLfound)                 histos[(string)"TauEvents_LL"          + hh]->Fill(Mbb, weight);//signal
	if(recoedtaunomt && LLfound)                 histos[(string)"TauEvents_noMT_LL"     + hh]->Fill(Mbb, weight);//signal
	if(recoedtau     && genhadtau==2)            histos[(string)"TauEvents_dTau"        + hh]->Fill(Mbb, weight);//signal
	if(recoedtaunomt && genhadtau==2)            histos[(string)"TauEvents_noMT_dTau"   + hh]->Fill(Mbb, weight);//signal
	if(recoedtau     && genhadtau==0)            histos[(string)"noISR_TauEventsBG"           + hh]->Fill(Mbb, weightnoISR);
	if(recoedtaunomt && genhadtau==0)            histos[(string)"noISR_TauEventsBG_noMT"      + hh]->Fill(Mbb, weightnoISR);
	if(recoedtau     && genhadtau==0 && LLfound) histos[(string)"noISR_TauEventsBG_LL"        + hh]->Fill(Mbb, weightnoISR);
	if(recoedtaunomt && genhadtau==0 && LLfound) histos[(string)"noISR_TauEventsBG_noMT_LL"   + hh]->Fill(Mbb, weightnoISR);
	if(recoedtau     && LLfound)                 histos[(string)"noISR_TauEvents_LL"          + hh]->Fill(Mbb, weightnoISR);
	if(recoedtaunomt && LLfound)                 histos[(string)"noISR_TauEvents_noMT_LL"     + hh]->Fill(Mbb, weightnoISR);
	if(recoedtau     && genhadtau==2)            histos[(string)"noISR_TauEvents_dTau"        + hh]->Fill(Mbb, weightnoISR);
	if(recoedtaunomt && genhadtau==2)            histos[(string)"noISR_TauEvents_noMT_dTau"   + hh]->Fill(Mbb, weightnoISR);
	if(fbTagError && !fMT2tree->misc.isData && btagSF!=0){
		if(recoedtau     && genhadtau==0)            histos[(string)"TauEventsBG_SFup"         + hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(recoedtaunomt && genhadtau==0)            histos[(string)"TauEventsBG_noMT_SFup"    + hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(recoedtau     && genhadtau==0 && LLfound) histos[(string)"TauEventsBG_LL_SFup"      + hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(recoedtaunomt && genhadtau==0 && LLfound) histos[(string)"TauEventsBG_noMT_LL_SFup" + hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(recoedtau     && LLfound)                 histos[(string)"TauEvents_LL_SFup"        + hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(recoedtaunomt && LLfound)                 histos[(string)"TauEvents_noMT_LL_SFup"   + hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(recoedtau     && genhadtau==2)            histos[(string)"TauEvents_dTau_SFup"      + hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(recoedtaunomt && genhadtau==2)            histos[(string)"TauEvents_noMT_dTau_SFup" + hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(recoedtau     && genhadtau==0)            histos[(string)"TauEventsBG_SFdown"        +hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(recoedtaunomt && genhadtau==0)            histos[(string)"TauEventsBG_noMT_SFdown"   +hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(recoedtau     && genhadtau==0 && LLfound) histos[(string)"TauEventsBG_LL_SFdown"     +hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(recoedtaunomt && genhadtau==0 && LLfound) histos[(string)"TauEventsBG_noMT_LL_SFdown"+hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(recoedtau     && LLfound)                 histos[(string)"TauEvents_LL_SFdown"       +hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(recoedtaunomt && LLfound)                 histos[(string)"TauEvents_noMT_LL_SFdown"  +hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(recoedtau     && genhadtau==2)            histos[(string)"TauEvents_dTau_SFdown"     +hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(recoedtaunomt && genhadtau==2)            histos[(string)"TauEvents_noMT_dTau_SFdown"+hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(recoedtau     && genhadtau==0)            histos[(string)"noISR_TauEventsBG_SFup"         + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(recoedtaunomt && genhadtau==0)            histos[(string)"noISR_TauEventsBG_noMT_SFup"    + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(recoedtau     && genhadtau==0 && LLfound) histos[(string)"noISR_TauEventsBG_LL_SFup"      + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(recoedtaunomt && genhadtau==0 && LLfound) histos[(string)"noISR_TauEventsBG_noMT_LL_SFup" + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(recoedtau     && LLfound)                 histos[(string)"noISR_TauEvents_LL_SFup"        + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(recoedtaunomt && LLfound)                 histos[(string)"noISR_TauEvents_noMT_LL_SFup"   + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(recoedtau     && genhadtau==2)            histos[(string)"noISR_TauEvents_dTau_SFup"      + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(recoedtaunomt && genhadtau==2)            histos[(string)"noISR_TauEvents_noMT_dTau_SFup" + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(recoedtau     && genhadtau==0)            histos[(string)"noISR_TauEventsBG_SFdown"        +hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
		if(recoedtaunomt && genhadtau==0)            histos[(string)"noISR_TauEventsBG_noMT_SFdown"   +hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
		if(recoedtau     && genhadtau==0 && LLfound) histos[(string)"noISR_TauEventsBG_LL_SFdown"     +hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
		if(recoedtaunomt && genhadtau==0 && LLfound) histos[(string)"noISR_TauEventsBG_noMT_LL_SFdown"+hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
		if(recoedtau     && LLfound)                 histos[(string)"noISR_TauEvents_LL_SFdown"       +hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
		if(recoedtaunomt && LLfound)                 histos[(string)"noISR_TauEvents_noMT_LL_SFdown"  +hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
		if(recoedtau     && genhadtau==2)            histos[(string)"noISR_TauEvents_dTau_SFdown"     +hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
		if(recoedtaunomt && genhadtau==2)            histos[(string)"noISR_TauEvents_noMT_dTau_SFdown"+hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
	}

	if(eventgood)               histos[(string)"WEvents"             +hh]->Fill(Mbb, weight);
	if(eventgood_noLL)          histos[(string)"WEvents_noLL"        +hh]->Fill(Mbb, weight);
	if(genhadtau==2)            histos[(string)"WEvents_dTau"        +hh]->Fill(Mbb, weight);//double tau == no LL, as we have at max 2 W decays
	if(acceptance)              histos[(string)"WEventsAcc"          +hh]->Fill(Mbb, weight);
	if(acceptance_noLL)         histos[(string)"WEventsAcc_noLL"     +hh]->Fill(Mbb, weight);
	if(leptfound)               histos[(string)"WEventsReco"         +hh]->Fill(Mbb, weight);
	if(leptfound&&genhadtau==2) histos[(string)"WEventsReco_dTau"    +hh]->Fill(Mbb, weight);
	if(leptfound_noLL)          histos[(string)"WEventsReco_noLL"    +hh]->Fill(Mbb, weight);
	if(leptfoundnomt)           histos[(string)"WEventsReco_noMT"     +hh]->Fill(Mbb, weight);
	if(leptfoundnomt_noLL)      histos[(string)"WEventsReco_noMT_noLL"+hh]->Fill(Mbb, weight);
	if(eventgood)               histos[(string)"noISR_WEvents"             +hh]->Fill(Mbb, weightnoISR);
	if(eventgood_noLL)          histos[(string)"noISR_WEvents_noLL"        +hh]->Fill(Mbb, weightnoISR);
	if(genhadtau==2)            histos[(string)"noISR_WEvents_dTau"        +hh]->Fill(Mbb, weightnoISR);
	if(acceptance)              histos[(string)"noISR_WEventsAcc"          +hh]->Fill(Mbb, weightnoISR);
	if(acceptance_noLL)         histos[(string)"noISR_WEventsAcc_noLL"     +hh]->Fill(Mbb, weightnoISR);
	if(leptfound)               histos[(string)"noISR_WEventsReco"         +hh]->Fill(Mbb, weightnoISR);
	if(leptfound&&genhadtau==2) histos[(string)"noISR_WEventsReco_dTau"    +hh]->Fill(Mbb, weightnoISR);
	if(leptfound_noLL)          histos[(string)"noISR_WEventsReco_noLL"    +hh]->Fill(Mbb, weightnoISR);
	if(leptfoundnomt)           histos[(string)"noISR_WEventsReco_noMT"     +hh]->Fill(Mbb, weightnoISR);
	if(leptfoundnomt_noLL)      histos[(string)"noISR_WEventsReco_noMT_noLL"+hh]->Fill(Mbb, weightnoISR);
	if(fbTagError && !fMT2tree->misc.isData && btagSF!=0){
		if(eventgood)               histos[(string)"WEvents_SFup"       +hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(eventgood_noLL)          histos[(string)"WEvents_noLL_SFup"  +hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(eventgood&&genhadtau==2) histos[(string)"WEvents_dTau_SFup"  +hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(eventgood)               histos[(string)"WEvents_SFdown"     +hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(eventgood_noLL)          histos[(string)"WEvents_noLL_SFdown"+hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(eventgood&&genhadtau==2) histos[(string)"WEvents_dTau_SFdown"+hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(eventgood)               histos[(string)"noISR_WEvents_SFup"       +hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(eventgood_noLL)          histos[(string)"noISR_WEvents_noLL_SFup"  +hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(eventgood&&genhadtau==2) histos[(string)"noISR_WEvents_dTau_SFup"  +hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(eventgood)               histos[(string)"noISR_WEvents_SFdown"     +hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
		if(eventgood_noLL)          histos[(string)"noISR_WEvents_noLL_SFdown"+hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
		if(eventgood&&genhadtau==2) histos[(string)"noISR_WEvents_dTau_SFdown"+hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
	}
	//add LLveto for recotaus
	if(norecotau &&  eventgood           ) histos[(string)"NoRecoButGenTauEvents"      + hh]->Fill(Mbb, weight);
	if(norecotau &&  eventgood && LLfound) histos[(string)"NoRecoButGenTauEvents_LL"   + hh]->Fill(Mbb, weight);
	if(norecotau &&  genhadtau==2        ) histos[(string)"NoRecoButGenTauEvents_dTau" + hh]->Fill(Mbb, weight);
	if(norecotau &&  eventgood           ) histos[(string)"noISR_NoRecoButGenTauEvents"      + hh]->Fill(Mbb, weightnoISR);
	if(norecotau &&  eventgood && LLfound) histos[(string)"noISR_NoRecoButGenTauEvents_LL"   + hh]->Fill(Mbb, weightnoISR);
	if(norecotau &&  genhadtau==2        ) histos[(string)"noISR_NoRecoButGenTauEvents_dTau" + hh]->Fill(Mbb, weightnoISR);
	if(fbTagError && !fMT2tree->misc.isData && btagSF!=0){
		if(norecotau &&  eventgood           ) histos[(string)"NoRecoButGenTauEvents_SFup"        + hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(norecotau &&  eventgood && LLfound) histos[(string)"NoRecoButGenTauEvents_LL_SFup"     + hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(norecotau &&  genhadtau==2        ) histos[(string)"NoRecoButGenTauEvents_dTau_SFup"   + hh]->Fill(Mbb, weight/btagSF*(btagSF+btagSFerr));
		if(norecotau &&  eventgood           ) histos[(string)"NoRecoButGenTauEvents_SFdown"      + hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(norecotau &&  eventgood && LLfound) histos[(string)"NoRecoButGenTauEvents_LL_SFdown"   + hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(norecotau &&  genhadtau==2        ) histos[(string)"NoRecoButGenTauEvents_dTau_SFdown" + hh]->Fill(Mbb, weight/btagSF*(btagSF-btagSFerr));
		if(norecotau &&  eventgood           ) histos[(string)"noISR_NoRecoButGenTauEvents_SFup"        + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(norecotau &&  eventgood && LLfound) histos[(string)"noISR_NoRecoButGenTauEvents_LL_SFup"     + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(norecotau &&  genhadtau==2        ) histos[(string)"noISR_NoRecoButGenTauEvents_dTau_SFup"   + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF+btagSFerr));
		if(norecotau &&  eventgood           ) histos[(string)"noISR_NoRecoButGenTauEvents_SFdown"      + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
		if(norecotau &&  eventgood && LLfound) histos[(string)"noISR_NoRecoButGenTauEvents_LL_SFdown"   + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
		if(norecotau &&  genhadtau==2        ) histos[(string)"noISR_NoRecoButGenTauEvents_dTau_SFdown" + hh]->Fill(Mbb, weightnoISR/btagSF*(btagSF-btagSFerr));
	}
	//these histograms are not used for prediction, but for important checks like tau purity
	histos[(string)"TrueRecoTauEvents"      + hh]->Fill(Mbb, weight*((double)recotaumatchedtrue));
	histos[(string)"FakeRecoTauEvents"      + hh]->Fill(Mbb, weight*((double)recotaumatchedfake));
	histos[(string)"AllRecoTauEvents"       + hh]->Fill(Mbb, weight*((double)(recotaumatchedtrue+recotaumatchedfake)));
	histos[(string)"AllRecoJetEvents"       + hh]->Fill(Mbb, weight*((double)fMT2tree->NJetsIDLoose));
	histos[(string)"AllGenTauAccEvents"     + hh]->Fill(Mbb, weight*((double)gentauwithinacc));
	histos[(string)"noISR_TrueRecoTauEvents"      + hh]->Fill(Mbb, weightnoISR*((double)recotaumatchedtrue));
	histos[(string)"noISR_FakeRecoTauEvents"      + hh]->Fill(Mbb, weightnoISR*((double)recotaumatchedfake));
	histos[(string)"noISR_AllRecoTauEvents"       + hh]->Fill(Mbb, weightnoISR*((double)(recotaumatchedtrue+recotaumatchedfake)));
	histos[(string)"noISR_AllRecoJetEvents"       + hh]->Fill(Mbb, weightnoISR*((double)fMT2tree->NJetsIDLoose));
	histos[(string)"noISR_AllGenTauAccEvents"     + hh]->Fill(Mbb, weightnoISR*((double)gentauwithinacc));
	} //while(myEvtList->GetEntry(counter++) !=-1)
	delete fMT2tree;
	delete fSamples[i].tree;
	} //for(size_t i = 0; i < fSamples.size(); ++i)

	//add overflow to last bin
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
		h->second->SetBinContent(0, h->second->GetBinContent(0 )+  h->second->GetBinContent(1)  );
		h->second->SetBinError(  0, sqrt(pow(h->second->GetBinError(0  ),2)+pow(h->second->GetBinError(1),2)));
	}


	cout << "add all samples to mc, etc." << endl;
	//add TTbar and SingleTop to Top
	//add WJets and Top to WandTop
	//add rest to noWandTop
	//add noWandTop and WandTop to mc
//	for(int i2 = 0; i2<signalregionsize; ++i2){
//	for(int i3 = 0; i3<HTbinsize;        ++i3){
//		if(fMET && i3==1) continue;
		string hsq   = string("_QCD");
		string hsw   = string("_WJets");
		string hsz   = string("_ZJets");
		string hstt  = string("_TTbar");
		string hsst  = string("_SingleTop");
		string hst   = string("_Top");
		string hswt  = string("_WandTop");
		string hsnwt = string("_noWandTop");
		string hso   = string("_Other");
		string hsmc  = string("_mc");
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
			string noisr = "noISR_";
			histos[noisr + (histonames[i0]+hsnwt)]->Add(histos[noisr + (histonames[i0]+hsq)],  1);//noWandTop
			histos[noisr + (histonames[i0]+hsnwt)]->Add(histos[noisr + (histonames[i0]+hsz)],  1);
			histos[noisr + (histonames[i0]+hsnwt)]->Add(histos[noisr + (histonames[i0]+hso)],  1);
			histos[noisr + (histonames[i0]+hst)  ]->Add(histos[noisr + (histonames[i0]+hsst)], 1);//Top
			histos[noisr + (histonames[i0]+hst)  ]->Add(histos[noisr + (histonames[i0]+hstt)], 1);
			histos[noisr + (histonames[i0]+hswt) ]->Add(histos[noisr + (histonames[i0]+hst)],  1);//WandTop
			histos[noisr + (histonames[i0]+hswt) ]->Add(histos[noisr + (histonames[i0]+hsw)],  1);
			histos[noisr + (histonames[i0]+hsmc) ]->Add(histos[noisr + (histonames[i0]+hswt)], 1);//mc
			histos[noisr + (histonames[i0]+hsmc) ]->Add(histos[noisr + 	(histonames[i0]+hsnwt)],1);
		}
	//}}
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
		string hs   = string("_") + sample_type[i1];
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
	}

	cout << "Saving." << endl;
	if(!fISRreweight) outputname = "NoISR_" + outputname;
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

	//do estimation
	cout << "do estimation" << endl;

	if(fRebin){
	for(int i1 = 0; i1<sampletypesize;   ++i1){
		string hs   = string("_") + sample_type[i1];
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
		histos["WEventsAcc"                       +hs] = RebinThisHistogram(histos["WEventsAcc"                       +hs], nrebin, rebin);
		histos["WEventsAcc_noLL"                  +hs] = RebinThisHistogram(histos["WEventsAcc_noLL"                  +hs], nrebin, rebin);
		histos["WEventsReco"                      +hs] = RebinThisHistogram(histos["WEventsReco"                      +hs], nrebin, rebin);
		histos["WEventsReco_dTau"                 +hs] = RebinThisHistogram(histos["WEventsReco_dTau"                 +hs], nrebin, rebin);
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
	}
	}

	//versions for SummaryTable
	//version == 0,1: with MCPred
	//version == 2,3: without MCPred
	//version == 0,2: with R_LL (i.e. 1-e / e)
	//version == 1,3: with e
//	GetTauEstimatehistograms,tefficiencies,rel_sys_err,rel_sys_err_bg,version,SummaryTable,YieldTable,PredictionCard,makeFullPrintout)
	GetTauEstimate(histos, teff, frel_sys_uncert, frel_sys_uncert_bg, 0, true, true, true, true, true, true);


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

}//void TauEstimationHiggs()

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
	string hshl = "";
	for(int nx = 1; nx<=hists[string("TauEvents"+hshl+"_mc")]->GetNbinsX(); ++nx){
		//first push_back the region information
		sr.push_back(0);
		if(fMET) htr.push_back(0);
		else     htr.push_back(1);
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
		*fLogStream << "Signal region " << 0 << ", HT bin " << "X" << " and MT2bin " << nx << " = (" << hists[string("TauEvents"+hshl+"_mc")]->GetBinLowEdge(nx) << "," << hists[string("TauEvents"+hshl+"_mc")]->GetBinLowEdge(nx) + hists[string("TauEvents"+hshl+"_mc")]->GetBinWidth(nx) << "):" << endl;

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
//	}}// i2, i3, i4

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

	if(!fRebin) *fLogStream << "$M_{bb}$ (GeV) & ";
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
		if(htr[n]==1 && fHT) htreg = "$H_{T}\\ge 750 GeV";
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
	if(fRebin) filename = filename + "LostTauPredictionHiggsFile.root";
	else       filename = filename + "FineBinnedLostTauPredictionHiggsFile.root";

	map<string, TH1D*>    hs;
	//for now it is only truth +/- error
	//and prediction +/- error stored
		string htreg;
		if(fMET) htreg = "lowHT";
		if(fHT) htreg = "highHT";
		if(!fRebin){
			string mapname;
			mapname = "Prediction_Tau_"+htreg;
			if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "",12, 20, 200); hs[mapname]->Sumw2();
			hs[mapname]->GetXaxis()->SetTitle("M_{bb} [GeV]");
			hs[mapname]->SetMarkerStyle(20);
			hs[mapname]->SetMarkerColor(kBlack);
			hs[mapname]->SetLineColor(kBlack);
			hs[mapname]->SetLineWidth(3);
			mapname = "SimulationTruth_Tau_"+htreg;
			if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "",12, 20, 200); hs[mapname]->Sumw2();
			hs[mapname]->GetXaxis()->SetTitle("M_{bb} [GeV]");
			hs[mapname]->SetMarkerColor(kBlue);
			hs[mapname]->SetLineColor(kBlue);
			hs[mapname]->SetFillColor(kBlue);
			hs[mapname]->SetLineWidth(0);
			hs[mapname]->SetFillStyle(3002); 
		}
		if(fRebin){
			mapname = "Prediction_everything";
			if(hs.count(mapname) == 0 ) hs[mapname] = new TH1D(mapname.c_str(), "", 6, 0, 6); hs[mapname]->Sumw2();
			hs[mapname]->GetXaxis()->SetTitle("signal region");
			hs[mapname]->GetXaxis()->SetBinLabel(1,"e,lH_{T}");    hs[mapname]->GetXaxis()->SetBinLabel(2,"e,hH_{T}");
			hs[mapname]->GetXaxis()->SetBinLabel(3,"#mu,lH_{T}");  hs[mapname]->GetXaxis()->SetBinLabel(4,"#mu,hH_{T}");
			hs[mapname]->GetXaxis()->SetBinLabel(5,"#tau,lH_{T}"); hs[mapname]->GetXaxis()->SetBinLabel(6,"#tau,hH_{T}");
			hs[mapname]->SetMarkerStyle(20);
			hs[mapname]->SetMarkerColor(kBlack);
			hs[mapname]->SetLineColor(kBlack);
			hs[mapname]->SetLineWidth(3);
			mapname = "SimulationTruth_everything";
			hs[mapname]->GetXaxis()->SetTitle("signal region");
			hs[mapname]->GetXaxis()->SetBinLabel(1,"e,lH_{T}");    hs[mapname]->GetXaxis()->SetBinLabel(2,"e,hH_{T}");
			hs[mapname]->GetXaxis()->SetBinLabel(3,"#mu,lH_{T}");  hs[mapname]->GetXaxis()->SetBinLabel(4,"#mu,hH_{T}");
			hs[mapname]->GetXaxis()->SetBinLabel(5,"#tau,lH_{T}"); hs[mapname]->GetXaxis()->SetBinLabel(6,"#tau,hH_{T}");
			hs[mapname]->SetMarkerColor(kBlue);
			hs[mapname]->SetLineColor(kBlue);
			hs[mapname]->SetFillColor(kBlue);
			hs[mapname]->SetLineWidth(0);
			hs[mapname]->SetFillStyle(3002); 
		}

	vector<double> dp, dpe, st, ste, mt2binning; dp.clear(); dpe.clear(); st.clear(); ste.clear(); mt2binning.clear();
	for(unsigned int n = 0; n<sr.size(); ++n){
		string lep = "Tau";
		string htreg;
		if(fMET) htreg = "lowHT";
		if(fHT) htreg = "highHT";
		int mt2lastbin = 9999;
		if(htreg!=oldHT) oldHT = htreg;
		if(lep!=oldlep) oldlep = lep;
		if(!fRebin){
			string mapname = "Prediction_"+lep+"_"+htreg;
			hs[mapname]->SetBinContent(sr[n]+1, datapred[n]);
			hs[mapname]->SetBinError(sr[n]+1, sqrt(pow(datapred_stat_err[n],2)+pow(datapred_syst_err[n],2)) );
			mapname = "SimulationTruth_"+lep+"_"+htreg;
			hs[mapname]->SetBinContent(sr[n]+1, mctruth[n]);
			hs[mapname]->SetBinError(sr[n]+1, mctrutherr[n]);
		} else{
			if(fMET && lep=="Ele"){
				string mapname = "Prediction_everything";
				hs[mapname]->SetBinContent(1, datapred[n]);
				hs[mapname]->SetBinError(1, sqrt(pow(datapred_stat_err[n],2)+pow(datapred_syst_err[n],2)) );
				mapname = "SimulationTruth_everything";
				hs[mapname]->SetBinContent(1, mctruth[n]);
				hs[mapname]->SetBinError(1, mctrutherr[n]);
			}
			if(fHT && lep=="Ele"){
				string mapname = "Prediction_everything";
				hs[mapname]->SetBinContent(2, datapred[n]);
				hs[mapname]->SetBinError(2, sqrt(pow(datapred_stat_err[n],2)+pow(datapred_syst_err[n],2)) );
				mapname = "SimulationTruth_everything";
				hs[mapname]->SetBinContent(2, mctruth[n]);
				hs[mapname]->SetBinError(2, mctrutherr[n]);
			}
			if(fMET && lep=="Muo"){
				string mapname = "Prediction_everything";
				hs[mapname]->SetBinContent(3, datapred[n]);
				hs[mapname]->SetBinError(3, sqrt(pow(datapred_stat_err[n],2)+pow(datapred_syst_err[n],2)) );
				mapname = "SimulationTruth_everything";
				hs[mapname]->SetBinContent(3, mctruth[n]);
				hs[mapname]->SetBinError(3, mctrutherr[n]);
			}
			if(fHT && lep=="Muo"){
				string mapname = "Prediction_everything";
				hs[mapname]->SetBinContent(4, datapred[n]);
				hs[mapname]->SetBinError(4, sqrt(pow(datapred_stat_err[n],2)+pow(datapred_syst_err[n],2)) );
				mapname = "SimulationTruth_everything";
				hs[mapname]->SetBinContent(4, mctruth[n]);
				hs[mapname]->SetBinError(4, mctrutherr[n]);
			}
			if(fMET && lep=="Tau"){
				string mapname = "Prediction_everything";
				hs[mapname]->SetBinContent(5, datapred[n]);
				hs[mapname]->SetBinError(5, sqrt(pow(datapred_stat_err[n],2)+pow(datapred_syst_err[n],2)) );
				mapname = "SimulationTruth_everything";
				hs[mapname]->SetBinContent(5, mctruth[n]);
				hs[mapname]->SetBinError(5, mctrutherr[n]);
			}
			if(fHT && lep=="Tau"){
				string mapname = "Prediction_everything";
				hs[mapname]->SetBinContent(6, datapred[n]);
				hs[mapname]->SetBinError(6, sqrt(pow(datapred_stat_err[n],2)+pow(datapred_syst_err[n],2)) );
				mapname = "SimulationTruth_everything";
				hs[mapname]->SetBinContent(6, mctruth[n]);
				hs[mapname]->SetBinError(6, mctrutherr[n]);
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
			//axis titles defined before
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
			string htreg;
			if(fMET) { htreg = "lowHT"; textht = "low H_{T}"; }
			if(fHT)  { htreg = "highHT"; textht = "high H_{T}"; }
			if(!fRebin){
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
				hsp = "Prediction_everything";
				hss = "SimulationTruth_everything";
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
				outname = plotdirectory+hsp + ".eps";
				c1->SaveAs(outname.c_str());
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
              *fLogStream << "\\begin{tabular}{lcccccccc}" << endl	     //deleted two cc
	                  << "\\hline\\hline"              << endl;
	} else {
              *fLogStream << "\\begin{tabular}{lccccccc}"  << endl	     //deleted two cc
	                  << "\\hline\\hline"              << endl;
	}

	if(!fRebin)    *fLogStream << "$M_{bb}$ (GeV)";
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
		if(htr[n]==1 && fHT) htreg = "$H_{T}\\ge 750 GeV";
		if(htreg!=oldHT){
			*fLogStream << " \\hline " << endl << "\\multicolumn{7}{l}{" <<  htreg << "} \\\\ " << endl << "\\hline" << endl;
			oldHT = htreg;
		}
		if(!fRebin){
		if(sigreg!=oldsr){
			*fLogStream << " \\hline  " << endl << sigreg << "\\\\"  << endl;
			oldsr = sigreg;
		}
		if(MT2up[n]==10000.) *fLogStream << "$" << int(MT2low[n]) << "-" << "\\infty$" << " " << setw(4) << " & ";
		else            *fLogStream << "$" << int(MT2low[n]) << "-" << int(MT2up[n]) << "$" << " " << setw(7) << " & ";
		}
		else *fLogStream << " " << setw(18) << sigreg << " & ";
		*fLogStream << fixed << setprecision(2)
		<< " " << setw(9) << numQCD[n] << " & " << " " << setw(7) << numZ[n] << " & " << " " << setw(7) << numW[n] << " &  " << " " << setw(8) << numT[n] << " & " << " " << setw(10) << numOther[n] << " & ";
		if(fPlotLLTab) *fLogStream << fixed << setprecision(2) << " " << setw(12) << numBGLL[n] << " & ";
		*fLogStream << fixed << setprecision(2) << " " << setw(10) << numMC[n] << "$\\pm" << " " << setw(7) << BGerr[n] << "$ & ";
		if(!fMCClosure) *fLogStream << " " << setw(10) << int(numData[n]) << " \\\\" << endl;
		else *fLogStream << " " << setw(10) << "-" << " \\\\" << endl;
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

//standard load function that reads out the samples.dat into the fSamples vector
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
			counter++;
			fSamples.push_back(s);
		}
	}
	if(fVerbose > 0) cout << "------------------------------------" << endl;
}
