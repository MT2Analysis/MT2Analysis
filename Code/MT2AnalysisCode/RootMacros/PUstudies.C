#include "TEfficiency.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TProfile.h"
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
#include "TLatex.h"
#include "TLegend.h"
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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path for MT2tree.hh
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/Utilities.hh"

//run via root -l -b -q PUstudies,C++

using namespace std;

void load(const char* filename = "samples/datasamples/samples_2141_dataonly.dat");
void PUstudies();

//struct that combines MT2trees with important information like x section
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


Int_t  fHTbin      = 9; //0: lowHT, 1: mediumHT, 2: highHT, -1: no HT cut, everything else: medium+high HT cut
Int_t fNJets       = -2;//positive: NJets==X, negative: NJets>=X
Int_t fNBJets      = -10;//same as NJets, if ==-10, no cut on NBJets
Bool_t fPUReweight = true; //do PU reweighting
Bool_t fbReweight  = true; //do BTV SF reweighting

Int_t fNBins       = 27; //number of bins for your histograms, default = 27
Double_t fBinMin   = 0;  //lower x-axis border (default = 0)
Double_t fBinMax   = 810;//upper x-axis border (default = 810)

//plots MT2 vs. NVertices or PUnumberOfInteractions
//for Znunu events
//to check stability of MT2 vs. PU
void PUstudies(){

  //plot and cut variable (i.e. MT2 and NVertices/PUnumInt)
  TString variable   = "misc.MT2";
  TString varname    = "MT2";
  TString puvariable = "pileUp.NVertices";
  TString puvarname  = "NVertices";
//  TString puvariable = "pileUp.PUnumInt";
//  TString puvarname  = "NPUInt";

  TString samples    = "samples/samples_Znunu_HT.dat";
  if(fHTbin==0) samples    = "samples/samples_Znunu_MET.dat";

	//event selection
  std::ostringstream cutStream;
  cutStream 
      << " " 
 //     << "misc.MT2>50"                           << "&&"
      << "misc.MET>30"                            << "&&"
      << "misc.Jet0Pass ==1"                       << "&&"
      << "misc.Jet1Pass ==1"                       << "&&"
      << "misc.SecondJPt  > 100"                   << "&&"
      << "misc.PassJet40ID ==1"                      << "&&"
      << "misc.Vectorsumpt < 70";
    cutStream << "&& misc.MinMetJetDPhi4Pt40 >0.3";
    cutStream << "&& NMuons==0 && NEles==0 && NTausIDLoose3Hits==0";
    if(fHTbin==0)       cutStream << "&&misc.MET>200&&misc.HT<750&&misc.HT>=450";
    else if(fHTbin==1)  cutStream << "&&misc.HT>=750&&misc.HT<=1050";
    else if(fHTbin==2)  cutStream << "&&misc.HT>1050";
    else if(fHTbin!=-1) cutStream << "&&misc.HT>750";


  std::ostringstream basecutStream;
  basecutStream 
      << " " 
      << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"  << "&&" // for rare SM samples
      << "misc.CSCTightHaloIDFlag == 0"             << "&&"
      << "misc.trackingFailureFlag==0"              << "&&"
      << "misc.eeBadScFlag==0"                      << "&&"
      << "misc.EcalDeadCellTriggerPrimitiveFlag==0" << "&&"
      << "misc.TrackingManyStripClusFlag==0"             << "&&"
      << "misc.TrackingTooManyStripClusFlag==0"          << "&&"
      << "misc.TrackingLogErrorTooManyClustersFlag==0"   << "&&"
      << "misc.CrazyHCAL==0";

	TString nJets, nJetsVar = "NJetsIDLoose40";
	if (fNJets>=10) {
	  nJets =  "(" + nJetsVar + TString::Format(">=%d",fNJets/10);
	  nJets += "&&"+ nJetsVar + TString::Format("<=%d",fNJets%10)+")";
	}
	else{
	  nJets = nJetsVar + (fNJets < 0 ? ">=" : "==");
	  nJets = nJets + TString::Format("%d",abs(fNJets));
	}

	TString nBJets = "NBJets40CSVM";    // nbjets = -10  --> >=0 b-tags
	nBJets += fNBJets < 0 ? ">=" : "==";
	nBJets += fNBJets==-10 ? "0" : TString::Format("%d",abs(fNBJets));

  std::ostringstream triggerStream;
  if(fHTbin==0){
  triggerStream << "( ( "
		<< "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
  		<< "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )"
		<< "||("
		<< "trigger.HLT_PFHT350_PFMET100_v3==1 || trigger.HLT_PFHT350_PFMET100_v4==1 || trigger.HLT_PFHT350_PFMET100_v5==1 || "
		<< "trigger.HLT_PFHT350_PFMET100_v6==1 || trigger.HLT_PFHT350_PFMET100_v7==1 || trigger.HLT_PFNoPUHT350_PFMET100_v1==1 || "
		<< "trigger.HLT_PFNoPUHT350_PFMET100_v3==1 || trigger.HLT_PFNoPUHT350_PFMET100_v4==1 ) )";
  } else {
  triggerStream << "( "
		<< "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
		<< "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1)";
  }
  TString trigger = triggerStream.str().c_str();
  TString cuts = cutStream.str().c_str();
  TString basecuts = basecutStream.str().c_str();

	//define histograms
  	TH1D*    hMT2_PU0toInf  = new TH1D   ("hMT2_PU0toInf" , "", fNBins, fBinMin, fBinMax );
  	TH1D*    hMT2_PU0to8    = new TH1D   ("hMT2_PU0to8"   , "", fNBins, fBinMin, fBinMax );
  	TH1D*    hMT2_PU9to12   = new TH1D   ("hMT2_PU9to12"  , "", fNBins, fBinMin, fBinMax );
  	TH1D*    hMT2_PU13to16  = new TH1D   ("hMT2_PU13to16" , "", fNBins, fBinMin, fBinMax );
  	TH1D*    hMT2_PU17to20  = new TH1D   ("hMT2_PU17to20" , "", fNBins, fBinMin, fBinMax );
  	TH1D*    hMT2_PU21to28  = new TH1D   ("hMT2_PU21to28" , "", fNBins, fBinMin, fBinMax );
  	TH1D*    hMT2_PU29toInf = new TH1D   ("hMT2_PU29toInf", "", fNBins, fBinMin, fBinMax );//for NVertices
  	TH1D*    hMT2_PU29to40  = new TH1D   ("hMT2_PU29to40" , "", fNBins, fBinMin, fBinMax );//for PUnumInt, PUtrueNumInt
  	TH1D*    hMT2_PU41toInf = new TH1D   ("hMT2_PU41toInf", "", fNBins, fBinMin, fBinMax );//for PUnumInt, PUtrueNumInt

	TH2D*    hMT2vsPU       = new TH2D   ("hMT2vsPU",       "", 160, 0, 800 , 30, 0, 60);

	hMT2_PU0toInf ->Sumw2(); hMT2_PU0toInf ->SetStats(false); hMT2_PU0toInf ->SetMarkerStyle(20), hMT2_PU0toInf ->SetMarkerColor(kBlack);
	hMT2_PU0to8   ->Sumw2(); hMT2_PU0to8   ->SetStats(false); hMT2_PU0to8   ->SetMarkerStyle(20), hMT2_PU0to8   ->SetMarkerColor(kYellow+3);
	hMT2_PU9to12  ->Sumw2(); hMT2_PU9to12  ->SetStats(false); hMT2_PU9to12  ->SetMarkerStyle(20), hMT2_PU9to12  ->SetMarkerColor(kCyan+2);
	hMT2_PU13to16 ->Sumw2(); hMT2_PU13to16 ->SetStats(false); hMT2_PU13to16 ->SetMarkerStyle(20), hMT2_PU13to16 ->SetMarkerColor(kMagenta+2);
	hMT2_PU17to20 ->Sumw2(); hMT2_PU17to20 ->SetStats(false); hMT2_PU17to20 ->SetMarkerStyle(20), hMT2_PU17to20 ->SetMarkerColor(kRed);
	hMT2_PU21to28 ->Sumw2(); hMT2_PU21to28 ->SetStats(false); hMT2_PU21to28 ->SetMarkerStyle(20), hMT2_PU21to28 ->SetMarkerColor(kBlue);
	hMT2_PU29toInf->Sumw2(); hMT2_PU29toInf->SetStats(false); hMT2_PU29toInf->SetMarkerStyle(20), hMT2_PU29toInf->SetMarkerColor(kGreen+2);
	hMT2_PU29to40 ->Sumw2(); hMT2_PU29to40 ->SetStats(false); hMT2_PU29to40 ->SetMarkerStyle(20), hMT2_PU29to40 ->SetMarkerColor(kGreen+2);
	hMT2_PU41toInf->Sumw2(); hMT2_PU41toInf->SetStats(false); hMT2_PU41toInf->SetMarkerStyle(20), hMT2_PU41toInf->SetMarkerColor(kOrange+3);
 	hMT2vsPU->Sumw2();       hMT2vsPU      ->SetStats(false);

	load(samples.Data());

    	hMT2vsPU->GetDirectory()->cd();

	for(size_t i = 0; i < fSamples.size(); ++i){

		Double_t weight=0;
		//get global weight
		if(fPUReweight) weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
		else            weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);

		TString theCuts = nJets + "&&" + nBJets;
		theCuts = theCuts + "&&" + cuts;
		if(basecuts!="") theCuts = theCuts + "&&" + basecuts;
		if(fSamples[i].type=="data" && trigger!="") theCuts += " &&("+trigger+")"; // triggers for data

		TString btagweight = "1.00"; //stored btag weights up to >=3, default is weight==1 to avoid non-existing weights
		if(fNBJets>=0 && fNBJets<=3) btagweight = TString::Format("SFWeight.BTagCSV40eq%d",abs(fNBJets));
		else if(fNBJets>=-3)         btagweight = TString::Format("SFWeight.BTagCSV40ge%d",abs(fNBJets));

		TString selection; TString var; int nev;
		TString Cutstemp;
		//all PU
		//define selection
		Cutstemp = theCuts;
		if(     fSamples[i].type!="data" && fPUReweight && fbReweight) selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), Cutstemp.Data());
		else if(fSamples[i].type!="data" && fPUReweight              ) selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",   weight,                    Cutstemp.Data());
		else if(fSamples[i].type!="data" &&                fbReweight) selection = TString::Format("(%.15f*%s) * (%s)",              weight, btagweight.Data(), Cutstemp.Data());
		else                                                           selection = TString::Format("(%.15f) * (%s)",                 weight,                    Cutstemp.Data()); 
		//define variable
		var  = TString::Format("%s>>+%s",variable.Data(),hMT2_PU0toInf->GetName());
		if(fVerbose>2) cout << "  +++++++ Drawing      " << var  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
		//fill histogram using tree->Draw()
		nev = fSamples[i].tree->Draw(var.Data(),selection.Data(),"goff");
		cout << nev << " events" << endl;
		cout << hMT2_PU0toInf->Integral() << " weighted events" << endl;
		// 0 <= PU <= 8
		Cutstemp = theCuts + "&&" + puvariable + "<=8";
		if(     fSamples[i].type!="data" && fPUReweight && fbReweight) selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), Cutstemp.Data());
		else if(fSamples[i].type!="data" && fPUReweight              ) selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",   weight,                    Cutstemp.Data());
		else if(fSamples[i].type!="data" &&                fbReweight) selection = TString::Format("(%.15f*%s) * (%s)",              weight, btagweight.Data(), Cutstemp.Data());
		else                                                           selection = TString::Format("(%.15f) * (%s)",                 weight,                    Cutstemp.Data()); 
		var  = TString::Format("%s>>+%s",variable.Data(),hMT2_PU0to8->GetName());
		if(fVerbose>2) cout << "  +++++++ Drawing      " << var  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
		nev = fSamples[i].tree->Draw(var.Data(),selection.Data(),"goff");
		cout << nev << " events" << endl;
		// 9 <= PU <= 12
		Cutstemp = theCuts + "&&" + puvariable + ">=9&&" + puvariable + "<=12";
		if(     fSamples[i].type!="data" && fPUReweight && fbReweight) selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), Cutstemp.Data());
		else if(fSamples[i].type!="data" && fPUReweight              ) selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",   weight,                    Cutstemp.Data());
		else if(fSamples[i].type!="data" &&                fbReweight) selection = TString::Format("(%.15f*%s) * (%s)",              weight, btagweight.Data(), Cutstemp.Data());
		else                                                           selection = TString::Format("(%.15f) * (%s)",                 weight,                    Cutstemp.Data()); 
		var  = TString::Format("%s>>+%s",variable.Data(),hMT2_PU9to12->GetName());
		if(fVerbose>2) cout << "  +++++++ Drawing      " << var  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
		nev = fSamples[i].tree->Draw(var.Data(),selection.Data(),"goff");
		// 13 <= PU <= 16
		Cutstemp = theCuts + "&&" + puvariable + ">=13&&" + puvariable + "<=16";
		if(     fSamples[i].type!="data" && fPUReweight && fbReweight) selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), Cutstemp.Data());
		else if(fSamples[i].type!="data" && fPUReweight              ) selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",   weight,                    Cutstemp.Data());
		else if(fSamples[i].type!="data" &&                fbReweight) selection = TString::Format("(%.15f*%s) * (%s)",              weight, btagweight.Data(), Cutstemp.Data());
		else                                                           selection = TString::Format("(%.15f) * (%s)",                 weight,                    Cutstemp.Data()); 
		var  = TString::Format("%s>>+%s",variable.Data(),hMT2_PU13to16->GetName());
		if(fVerbose>2) cout << "  +++++++ Drawing      " << var  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
		nev = fSamples[i].tree->Draw(var.Data(),selection.Data(),"goff");
		cout << nev << " events" << endl;
		// 17 <= PU <= 20
		Cutstemp = theCuts + "&&" + puvariable + ">=17&&" + puvariable + "<=20";
		if(     fSamples[i].type!="data" && fPUReweight && fbReweight) selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), Cutstemp.Data());
		else if(fSamples[i].type!="data" && fPUReweight              ) selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",   weight,                    Cutstemp.Data());
		else if(fSamples[i].type!="data" &&                fbReweight) selection = TString::Format("(%.15f*%s) * (%s)",              weight, btagweight.Data(), Cutstemp.Data());
		else                                                           selection = TString::Format("(%.15f) * (%s)",                 weight,                    Cutstemp.Data()); 
		var  = TString::Format("%s>>+%s",variable.Data(),hMT2_PU17to20->GetName());
		if(fVerbose>2) cout << "  +++++++ Drawing      " << var  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
		nev = fSamples[i].tree->Draw(var.Data(),selection.Data(),"goff");
		cout << nev << " events" << endl;
		// 21 <= PU <= 28
		Cutstemp = theCuts + "&&" + puvariable + ">=21&&" + puvariable + "<=28";
		if(     fSamples[i].type!="data" && fPUReweight && fbReweight) selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), Cutstemp.Data());
		else if(fSamples[i].type!="data" && fPUReweight              ) selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",   weight,                    Cutstemp.Data());
		else if(fSamples[i].type!="data" &&                fbReweight) selection = TString::Format("(%.15f*%s) * (%s)",              weight, btagweight.Data(), Cutstemp.Data());
		else                                                           selection = TString::Format("(%.15f) * (%s)",                 weight,                    Cutstemp.Data()); 
		var  = TString::Format("%s>>+%s",variable.Data(),hMT2_PU21to28->GetName());
		if(fVerbose>2) cout << "  +++++++ Drawing      " << var  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
		nev = fSamples[i].tree->Draw(var.Data(),selection.Data(),"goff");
		cout << nev << " events" << endl;
		// 29 <= PU
		Cutstemp = theCuts + "&&" + puvariable + ">=29";
		if(     fSamples[i].type!="data" && fPUReweight && fbReweight) selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), Cutstemp.Data());
		else if(fSamples[i].type!="data" && fPUReweight              ) selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",   weight,                    Cutstemp.Data());
		else if(fSamples[i].type!="data" &&                fbReweight) selection = TString::Format("(%.15f*%s) * (%s)",              weight, btagweight.Data(), Cutstemp.Data());
		else                                                           selection = TString::Format("(%.15f) * (%s)",                 weight,                    Cutstemp.Data()); 
		var  = TString::Format("%s>>+%s",variable.Data(),hMT2_PU29toInf->GetName());
		if(fVerbose>2) cout << "  +++++++ Drawing      " << var  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
		nev = fSamples[i].tree->Draw(var.Data(),selection.Data(),"goff");
		cout << nev << " events" << endl;
		// 29 <= PU <= 40
		Cutstemp = theCuts + "&&" + puvariable + ">=29&&" + puvariable + "<=40";
		if(     fSamples[i].type!="data" && fPUReweight && fbReweight) selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), Cutstemp.Data());
		else if(fSamples[i].type!="data" && fPUReweight              ) selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",   weight,                    Cutstemp.Data());
		else if(fSamples[i].type!="data" &&                fbReweight) selection = TString::Format("(%.15f*%s) * (%s)",              weight, btagweight.Data(), Cutstemp.Data());
		else                                                           selection = TString::Format("(%.15f) * (%s)",                 weight,                    Cutstemp.Data()); 
		var  = TString::Format("%s>>+%s",variable.Data(),hMT2_PU29to40->GetName());
		if(fVerbose>2) cout << "  +++++++ Drawing      " << var  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
		nev = fSamples[i].tree->Draw(var.Data(),selection.Data(),"goff");
		cout << nev << " events" << endl;
		// 41 <= PU
		Cutstemp = theCuts + "&&" + puvariable + ">=41";
		if(     fSamples[i].type!="data" && fPUReweight && fbReweight) selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), Cutstemp.Data());
		else if(fSamples[i].type!="data" && fPUReweight              ) selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",   weight,                    Cutstemp.Data());
		else if(fSamples[i].type!="data" &&                fbReweight) selection = TString::Format("(%.15f*%s) * (%s)",              weight, btagweight.Data(), Cutstemp.Data());
		else                                                           selection = TString::Format("(%.15f) * (%s)",                 weight,                    Cutstemp.Data()); 
		var  = TString::Format("%s>>+%s",variable.Data(),hMT2_PU41toInf->GetName());
		if(fVerbose>2) cout << "  +++++++ Drawing      " << var  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
		nev = fSamples[i].tree->Draw(var.Data(),selection.Data(),"goff");
		cout << nev << " events" << endl;

		// MT2 vs. PU
		Cutstemp = theCuts;
		if(     fSamples[i].type!="data" && fPUReweight && fbReweight) selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), Cutstemp.Data());
		else if(fSamples[i].type!="data" && fPUReweight              ) selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",   weight,                    Cutstemp.Data());
		else if(fSamples[i].type!="data" &&                fbReweight) selection = TString::Format("(%.15f*%s) * (%s)",              weight, btagweight.Data(), Cutstemp.Data());
		else                                                           selection = TString::Format("(%.15f) * (%s)",                 weight,                    Cutstemp.Data()); 
		TString tempvar = puvariable + ":misc.MT2";
		var  = TString::Format("%s>>+%s",tempvar.Data(),hMT2vsPU->GetName());
		if(fVerbose>2) cout << "  +++++++ Drawing      " << var  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
		nev = fSamples[i].tree->Draw(var.Data(),selection.Data(),"goff");
		cout << nev << " events" << endl;
	}

	TProfile *hMT2vsPUMT2proj = (TProfile*)hMT2vsPU->ProfileX(); hMT2vsPUMT2proj->SetStats(false);
	TProfile *hMT2vsPUPUproj  = (TProfile*)hMT2vsPU->ProfileY(); hMT2vsPUPUproj ->SetStats(false);

	TString puvarsavename = puvarname;

	if(fHTbin==0)       puvarsavename += "_lowHT";
	else if(fHTbin==1)  puvarsavename += "_mediumHT";
	else if(fHTbin==2)  puvarsavename += "_highHT";
	else if(fHTbin==-1) puvarsavename += "_allHT";
	else                puvarsavename += "_mediumhighHT";


	TString filename = "PUstudies/PUstudies_"+puvarsavename+".root";
	TFile* newfile = new TFile (filename.Data(), "RECREATE");
	newfile->cd();
  	hMT2_PU0toInf   ->Write();
  	hMT2_PU0to8     ->Write();
  	hMT2_PU9to12    ->Write();
  	hMT2_PU13to16   ->Write();
  	hMT2_PU17to20   ->Write();
  	hMT2_PU21to28   ->Write();
  	hMT2_PU29toInf  ->Write();
  	hMT2_PU29to40   ->Write();
  	hMT2_PU41toInf  ->Write();
	hMT2vsPU        ->Write();
	hMT2vsPUMT2proj ->Write();
	hMT2vsPUPUproj  ->Write();
	newfile->Close();
	cout << "histograms stored in " << newfile->GetName() << endl;

	//create normalized histograms to compare shapes
  	TH1D *hMT2_PU0toInf_norm   = (TH1D*)hMT2_PU0toInf  ->Clone("hMT2_PU0toInf_norm");
  	TH1D *hMT2_PU0to8_norm     = (TH1D*)hMT2_PU0to8    ->Clone("hMT2_PU0to8_norm");
  	TH1D *hMT2_PU9to12_norm    = (TH1D*)hMT2_PU9to12   ->Clone("hMT2_PU9to12_norm");
  	TH1D *hMT2_PU13to16_norm   = (TH1D*)hMT2_PU13to16  ->Clone("hMT2_PU13to16_norm");
  	TH1D *hMT2_PU17to20_norm   = (TH1D*)hMT2_PU17to20  ->Clone("hMT2_PU17to20_norm");
  	TH1D *hMT2_PU21to28_norm   = (TH1D*)hMT2_PU21to28  ->Clone("hMT2_PU21to28_norm");
  	TH1D *hMT2_PU29toInf_norm  = (TH1D*)hMT2_PU29toInf ->Clone("hMT2_PU29toInf_norm");
  	TH1D *hMT2_PU29to40_norm   = (TH1D*)hMT2_PU29to40  ->Clone("hMT2_PU29to40_norm");
  	TH1D *hMT2_PU41toInf_norm  = (TH1D*)hMT2_PU41toInf ->Clone("hMT2_PU41toInf_norm");
	hMT2_PU0toInf_norm ->Scale(1.0/hMT2_PU0toInf_norm ->Integral());
	hMT2_PU0to8_norm   ->Scale(1.0/hMT2_PU0to8_norm   ->Integral());
	hMT2_PU9to12_norm  ->Scale(1.0/hMT2_PU9to12_norm  ->Integral());
	hMT2_PU13to16_norm ->Scale(1.0/hMT2_PU13to16_norm ->Integral());
	hMT2_PU17to20_norm ->Scale(1.0/hMT2_PU17to20_norm ->Integral());
	hMT2_PU21to28_norm ->Scale(1.0/hMT2_PU21to28_norm ->Integral());
	hMT2_PU29toInf_norm->Scale(1.0/hMT2_PU29toInf_norm->Integral());
	hMT2_PU29to40_norm ->Scale(1.0/hMT2_PU29to40_norm ->Integral());
	hMT2_PU41toInf_norm->Scale(1.0/hMT2_PU41toInf_norm->Integral());
	hMT2_PU0to8_norm   ->Divide(hMT2_PU0toInf_norm);
	hMT2_PU9to12_norm  ->Divide(hMT2_PU0toInf_norm);
	hMT2_PU13to16_norm ->Divide(hMT2_PU0toInf_norm);
	hMT2_PU17to20_norm ->Divide(hMT2_PU0toInf_norm);
	hMT2_PU21to28_norm ->Divide(hMT2_PU0toInf_norm);
	hMT2_PU29toInf_norm->Divide(hMT2_PU0toInf_norm);
	hMT2_PU29to40_norm ->Divide(hMT2_PU0toInf_norm);
	hMT2_PU41toInf_norm->Divide(hMT2_PU0toInf_norm);
	hMT2_PU0toInf_norm ->Divide(hMT2_PU0toInf_norm);


	//the following makes the plot
	//as this is a copy from a MassPlotter.cc function, no detailed comments

	TCanvas* c1 = new TCanvas("c_ratio","",0,0,600,600);
	c1->SetFrameLineWidth(1);
	c1 -> cd();

	TH1D *haxis = (TH1D*)hMT2_PU0toInf->Clone("haxis");
 	TPad *p_plot  = new TPad("plotpad",  "Pad containing the overlay plot", 0,0.211838,1,1);
	p_plot->SetLeftMargin(0.131579);
	p_plot->SetRightMargin(0.08);
	p_plot->SetTopMargin(0.06895515);
	p_plot->SetBottomMargin(0.07206074);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio",   0,0.01863354,0.9967105,0.2189441);
	p_ratio->SetLeftMargin(0.1336634);	
	p_ratio->SetRightMargin(0.075);
	p_ratio->SetTopMargin(0.06976745);
	p_ratio->SetBottomMargin(0.2790698);
	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	gPad->SetLogy(1);
	gPad->SetFillStyle(0);
	gStyle->SetErrorX(0);

	//hardcoded
	haxis->SetMaximum(2.*1000.);//adjusted by hand: medium+high HT
	haxis->SetMinimum(1.*0.01);//adjusted by hand
	stringstream yTitle;
	if(fabs(haxis->GetBinWidth(1) -haxis->GetBinWidth(haxis->GetNbinsX()-1))<0.01){
		double binwidth = haxis->GetBinWidth(1);
		yTitle.precision(3);
		yTitle << "Events";
		yTitle << " / ";
		yTitle << binwidth;
		yTitle << " GeV";
	}
	haxis->GetYaxis()->SetTitle(yTitle.str().c_str());
	haxis->GetYaxis()->SetLabelSize(0.05);
	haxis->GetYaxis()->SetTitleSize(0.05);
	haxis->GetYaxis()->SetTitleOffset(1.3);
	haxis->Draw("axis");
  	hMT2_PU0toInf   ->Draw("sameX0");
  	hMT2_PU0to8     ->Draw("sameX0");
  	hMT2_PU9to12    ->Draw("sameX0");
  	hMT2_PU13to16   ->Draw("sameX0");
  	hMT2_PU17to20   ->Draw("sameX0");
  	hMT2_PU21to28   ->Draw("sameX0");
  	if(puvarname=="NVertices") hMT2_PU29toInf  ->Draw("sameX0");//for NVertices
	else{
  	hMT2_PU29to40   ->Draw("sameX0");//for N true PU Int
  	hMT2_PU41toInf  ->Draw("sameX0");//for N true PU Int
	}
	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.03);

	TString text;
	if (fNJets>=10)
	  text = TString::Format("%d-%d jets",fNJets/10,fNJets%10);
	else
	  text = fNJets < 0 ? TString::Format("#geq %d jets",abs(fNJets)) : TString::Format("%d jets",abs(fNJets));
	text += fNBJets==-10 ? "" : fNBJets < 0 ? TString::Format(", #geq %d b-tag",abs(fNBJets)) : TString::Format(", %d b-tag",abs(fNBJets));
	TitleBox.DrawLatex(0.13,0.943,text.Data());
	TLatex LumiBox;
	LumiBox.SetNDC();
	LumiBox.SetTextSize(0.0305);
	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
	LumiBox.DrawLatex(0.68,0.943,"#sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//standard
 	p_plot ->Draw();
	gPad->RedrawAxis();

	TLegend* leg = new TLegend(.71,.68,.91,.92);
	leg->SetName("legend");
	leg -> SetFillColor(0);
	leg -> SetBorderSize(0);
	// NVertices
	if(puvarname=="NVertices"){
	leg ->AddEntry(hMT2_PU0toInf,  "Z_{#nu#nu} incl",               "p");
	leg ->AddEntry(hMT2_PU0to8,    "Z_{#nu#nu} NV #leq 8",          "p");
	leg ->AddEntry(hMT2_PU9to12,   "Z_{#nu#nu} 9 #leq NV #leq 12",  "p");
	leg ->AddEntry(hMT2_PU13to16,  "Z_{#nu#nu} 13 #leq NV #leq 16", "p");
	leg ->AddEntry(hMT2_PU17to20,  "Z_{#nu#nu} 17 #leq NV #leq 20", "p");
	leg ->AddEntry(hMT2_PU21to28,  "Z_{#nu#nu} 21 #leq NV #leq 28", "p");
	leg ->AddEntry(hMT2_PU29toInf, "Z_{#nu#nu} NV #geq 29",         "p");
	} else {
	//PU int
	leg ->AddEntry(hMT2_PU0toInf,  "Z_{#nu#nu} incl",               "p");
	leg ->AddEntry(hMT2_PU0to8,    "Z_{#nu#nu} PU #leq 8",          "p");
	leg ->AddEntry(hMT2_PU9to12,   "Z_{#nu#nu} 9 #leq PU #leq 12",  "p");
	leg ->AddEntry(hMT2_PU13to16,  "Z_{#nu#nu} 13 #leq PU #leq 16", "p");
	leg ->AddEntry(hMT2_PU17to20,  "Z_{#nu#nu} 17 #leq PU #leq 20", "p");
	leg ->AddEntry(hMT2_PU21to28,  "Z_{#nu#nu} 21 #leq PU #leq 28", "p");
	leg ->AddEntry(hMT2_PU29to40,  "Z_{#nu#nu} 29 #leq PU #leq 40", "p");
	leg ->AddEntry(hMT2_PU41toInf, "Z_{#nu#nu} PU #geq 41",         "p");
	}
	leg -> Draw();

 	p_ratio ->cd();
	gStyle->SetErrorX(0);
	TH1D *haxis_norm = (TH1D*)hMT2_PU0toInf_norm->Clone("haxis_norm");
	haxis_norm->GetYaxis()->SetRangeUser(0.5,1.5);
	haxis_norm->GetXaxis()->SetLabelSize( 0.);
	if(puvarname=="NVertices") haxis_norm->GetYaxis()->SetTitle("NV bin / Incl");
	else haxis_norm->GetYaxis()->SetTitle("PU bin / Incl");
	haxis_norm->GetXaxis()->SetTitle("M_{T2} [GeV]");
	haxis_norm->GetXaxis()->SetTitleSize(0.2);
	haxis_norm->GetXaxis()->SetTitleOffset(0.5);
	haxis_norm->GetYaxis()->SetLabelSize(0.19);
	haxis_norm->GetXaxis()->SetTickLength(0.09);
	haxis_norm->GetYaxis()->SetTitleSize(0.15);
	haxis_norm->GetYaxis()->SetTitleOffset(0.36);
	haxis_norm->GetYaxis()->SetNdivisions(509);
	haxis_norm ->Draw("axis");

	hMT2_PU0toInf_norm->SetFillStyle(3001);
  	hMT2_PU0toInf_norm   ->Draw("E2same");//==1 by definition
  	hMT2_PU0to8_norm     ->Draw("X0same");
  	hMT2_PU9to12_norm    ->Draw("X0same");
  	hMT2_PU13to16_norm   ->Draw("X0same");
  	hMT2_PU17to20_norm   ->Draw("X0same");
  	hMT2_PU21to28_norm   ->Draw("X0same");
  	if(puvarname=="NVertices") hMT2_PU29toInf_norm  ->Draw("X0same");//for NVertices
	else{
  	hMT2_PU29to40_norm   ->Draw("X0same");//for N true PU Int
  	hMT2_PU41toInf_norm  ->Draw("X0same");//for N true PU Int
	}

	TLine *l3 = new TLine(haxis->GetXaxis()->GetXmin(), 1.00, haxis->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	Util::Print(c1, "PUstudies_"+puvarsavename, "PUstudies");

}//void PUstudies()

//read out samples.dat
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
