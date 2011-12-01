#include "TFile.h"
#include "TLegend.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <vector>
#include "TROOT.h"
//#include "../MT2Code/include/MT2Shapes.hh"
using namespace std;
using namespace RooFit ;

// User Input:  ----------------------------------------------
// -----
//TString fSamplesRemovedPhotons   ="./samples/samples_RemovedPhotonsHerwig.dat";
//TString fSamplesRemovedPhotons     ="./samples/samples_RemovedPhotons.dat";
//TString fSamplesHadronic           ="./samples/samples_PhotonsHadronic.dat";
TString fSamplesRemovedPhotons     ="./samples/samples_RemovedPhotonsMGQCD.dat";
TString fSamplesHadronic           ="./samples/samples_PhotonsHadronicMGQCD.dat";
// -----
Bool_t  fDoPhotonSigmaIEtaIEta     = true; 
Bool_t  fDoPhotonSignalRegion      = false; 
Bool_t  fDoHadronicSignalRegion    = false;
Bool_t  fDoPrediction              = false;
Int_t   fVerbose                   = 6;
Bool_t  fSaveResults               = true;
Bool_t  fMT2b                      = true;
Bool_t  fSeparateEBEE              = true;
Bool_t  fSaveZnunuToGammaRatio     = true;
Bool_t  fDraw                      = true;
Bool_t  fAddFitIntrinsicUncert     = true;
Float_t fFitIntrinsicUncert        = 0.1;
TString fOutDir                    = "GammaJetsPrediction/20111130_MG_Photons_MG_QCD_test";
Bool_t  fWriteToFile               = false;
// MT2b--------
Bool_t  fMT2b                      = true;
Bool_t  fDoBRatioCorrection        = true;
Float_t fRB_MCtoData               = 1.23;
Float_t fRB_MCtoDataErr            = 0.27;
// -------
Float_t fHTmin                     = 700;
Float_t fHTmax                     = 1E+8;
Float_t fMT2min                    = 200;
Float_t fMT2max                    = 1E+8;
// ------

// Global Variables ------------------------------------
std::ostringstream  fTriggerStream;
std::ostringstream  fCutStreamPhotons;
std::ostringstream  fCutStreamPhotonsMT2;
std::ostringstream  fCutStreamSignal;
std::ostringstream* fLogStream     = 0;

// Cut Streams
void DefineCutStreams(float HTmin, float HTmax, float MT2min, float MT2max){
	// Trigger Stream ---------------------------------------------------------------
	//
	fTriggerStream << "((trigger.HLT_HT440_v2 ==1 && misc.Run<161216)" << "||"
		<< "(trigger.HLT_HT450_v2 ==1 && (misc.Run>=161216 && misc.Run< 163269))" << "||"
		<< "(trigger.HLT_HT500_v3 ==1 && (misc.Run>=163269 && misc.Run<=163869))" << "||"
		<< "(trigger.HLT_HT500_v4 ==1 && (misc.Run>=165088 && misc.Run< 165970))" << "||"
		<< "(trigger.HLT_HT550_v5 ==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" << "||"
		<< "(trigger.HLT_HT550_v6 ==1 && (misc.Run==166346))" << "||"
		<< "(trigger.HLT_HT550_v7 ==1 && (misc.Run>=167078 && misc.Run< 170249))" << "||"
		<< "(trigger.HLT_HT550_v8 ==1 && (misc.Run>=170249 && misc.Run< 173236))" << "||"
		<< "(trigger.HLT_HT600_v1 ==1 && (misc.Run>=173236 && misc.Run< 178420))" << "||"
		<< "(trigger.HLT_HT650_v4 ==1 && (misc.Run>=178420 && misc.Run< 179959))" << "||"
		<< "(trigger.HLT_PFHT650_v1==1 && misc.Run>=179959)" << " )";

	// CutStream for SigmaIEtaIEta ------------------------------------------
	fCutStreamPhotons << " " 
	  << "misc.MET>=30"                                                << "&&"
	  << "misc.HT >= " << HTmin << " && misc.HT<= " << HTmax           << "&&"
	  << "(NEles==0  || ele[0].lv.Pt()<10)"                            << "&&"
	  << "(NMuons==0 || muo[0].lv.Pt()<10)"                            << "&&"
	  << "misc.Jet0Pass ==1"                                           << "&&"
	  << "misc.Jet1Pass ==1"                                           << "&&"
	  << "misc.SecondJPt >100"                                         << "&&"
	  << "misc.PassJetID ==1"                                          << "&&"
	  << "misc.Vectorsumpt < 70"                                       << "&&"
	  << "NPhotons ==1 "                                               << "&&"
	  << "(misc.ProcessID!=6||photon[0].MCmatchexitcode!=1)"           << "&&"
	  << "rawpfmet[0].Pt()<50"                                         << "&&"
	  // Noise
	  << "misc.HBHENoiseFlag ==0"                                      << "&&"
	  << "misc.CSCTightHaloID==0"                                      << "&&"
	  << "misc.CrazyHCAL==0";

	if(fMT2b){
		fCutStreamPhotons << " &&"
		  << "NJetsIDLoose40 >=4"                                          << "&&"
		  << "(misc.MinMetJetDPhi >0.3||misc.MinMetJetDPhiIndex>3)"        << "&&"
		  << "misc.LeadingJPt >150";
	}else{
		fCutStreamPhotons << " &&"
		  << "NJetsIDLoose40 >=3"                                          << "&&"
		  << "misc.MinMetJetDPhi >0.3";
	}

	// CutStream for Photon Signal Region ------------------------------------------ 
	fCutStreamPhotonsMT2 << " " 
	  << "misc.MET>=30"                                                 << "&&"
	  << "misc.HT >=" << HTmin  << " && misc.HT <=" << HTmax            << "&&"
	  << "misc.MT2>=" << MT2min << " && misc.MT2<=" << MT2max           << "&&"
	  << "(NEles==0  || ele[0].lv.Pt()<10)"                             << "&&"
	  << "(NMuons==0 || muo[0].lv.Pt()<10)"                             << "&&"
	  << "misc.Jet0Pass ==1"                                            << "&&"
	  << "misc.Jet1Pass ==1"                                            << "&&"
	  << "misc.SecondJPt >100"                                          << "&&"
	  << "misc.PassJetID ==1"                                           << "&&"
	  << "misc.Vectorsumpt < 70"                                        << "&&"
	  << "NPhotons ==1 "                                                << "&&"
	  << "photon[0].isEGMloose==1"                                      << "&&"
	  << "(misc.ProcessID!=6||photon[0].MCmatchexitcode!=1)"            << "&&"
	  << "rawpfmet[0].Pt()<50"                                          << "&&"
	  // Noise
	  << "misc.HBHENoiseFlag ==0"                                      << "&&"
	  << "misc.CSCTightHaloID==0"                                      << "&&"
	  << "misc.CrazyHCAL==0";

	if(fMT2b){
		fCutStreamPhotonsMT2 << " &&"
		  << "NJetsIDLoose40 >=4"                                          << "&&"
		  << "(misc.MinMetJetDPhi >0.3||misc.MinMetJetDPhiIndex>3)"        << "&&"
		  << "misc.LeadingJPt >150";
	}else{
		fCutStreamPhotonsMT2 << " &&"
		  << "NJetsIDLoose40 >=3"                                          << "&&"
		  << "misc.MinMetJetDPhi >0.3";
	}

	// CutStream for Hadronic Signal Region ----------------------------------------------
	fCutStreamSignal << " " 
	  << "misc.MET>=30"                                                << "&&"
	  << "misc.HT >=" << HTmin  << " && misc.HT <=" << HTmax           << "&&"
	  << "misc.MT2>=" << MT2min << " && misc.MT2<=" << MT2max          << "&&"
	  << "(NEles==0  || ele[0].lv.Pt()<10)"                            << "&&"
	  << "(NMuons==0 || muo[0].lv.Pt()<10)"                            << "&&"
	  << "misc.Jet0Pass ==1"                                           << "&&"
	  << "misc.Jet1Pass ==1"                                           << "&&"
	  << "misc.SecondJPt  >100"                                        << "&&"
	  << "misc.PassJetID ==1"                                          << "&&"
	  << "misc.Vectorsumpt < 70"                                       << "&&"
	  // Noise
	  << "misc.HBHENoiseFlag ==0"                                      << "&&"
	  << "misc.CSCTightHaloID==0"                                      << "&&"
	  << "misc.CrazyHCAL==0";
	if(fMT2b){
		fCutStreamSignal << " &&"
		  << "NBJets>0"                                                    << "&&"
		  << "NJetsIDLoose40 >=4"                                          << "&&"
		  << "(misc.MinMetJetDPhi >0.3||misc.MinMetJetDPhiIndex>3)"        << "&&"
		  << "misc.LeadingJPt >150";
	}else{
		fCutStreamSignal << " &&"
		  << "NJetsIDLoose40 >=3"                                          << "&&"
		  << "misc.MinMetJetDPhi >0.3";
	}
}

// *********************************** run_GammaJetsToZnunu: specify user input here *************************************************
void run_GammaJetsToZnunu(){
	gSystem->Load("libPhysics");
	gSystem->Load("libRooFit") ;
	gSystem->CompileMacro("../MT2Code/src/MT2Shapes.cc", "k");
	
	// logStream
	fLogStream = new std::ostringstream();
	
	// define cutsteams
	DefineCutStreams(fHTmin, fHTmax, fMT2min, fMT2max);

	// fix output dir
	if(fHTmax <10000) fOutDir= TString::Format("%s_%d_HT_%d",  fOutDir.Data(), abs(fHTmin),  abs(fHTmax));
	else              fOutDir= TString::Format("%s_%d_HT_%s",  fOutDir.Data(), abs(fHTmin),  "Inf");
	if(fMT2max<10000) fOutDir= TString::Format("%s_%d_MT2_%d", fOutDir.Data(), abs(fMT2min), abs(fMT2max));
	else              fOutDir= TString::Format("%s_%d_MT2_%s", fOutDir.Data(), abs(fMT2min), "Inf");
	
	if(!fMT2b || !fDoBRatioCorrection){
		fRB_MCtoData=1; fRB_MCtoDataErr=0;
	}

	// log MT2 and HT cuts
	*fLogStream << "------------------------------------------------------------------------------------------------" << endl;
	*fLogStream << "+++ new Znunu with Gamma+jets prediction                                                     +++" << endl;
	*fLogStream << "+++ outputdir: " << fOutDir <<                                                              "+++" << endl; 
	*fLogStream << "------------------------------------------------------------------------------------------------" << endl;
	if(fMT2b){
	*fLogStream << "+++ using Photon+jets MC to data b-tag correction: " << fRB_MCtoData << " pm " << fRB_MCtoDataErr << "+++" << endl;
	*fLogStream << "------------------------------------------------------------------------------------------------" << endl;
	}


	// new prediction ------------------------------------------------------
	Prediction* prediction = new Prediction();
	prediction->fVerbose=fVerbose;
	prediction->fSave   =fSaveResults;
	prediction->fOutputDir=fOutDir;


	// SigmaIEtaIEta Fit *********************************************************************
	// Get Photon Normalization: EB
	if(fDoPhotonSigmaIEtaIEta){
		std::ostringstream cutStreamPhotons_EB;
		if(fSeparateEBEE){
			cutStreamPhotons_EB << fCutStreamPhotons.str() << "&&abs(photon[0].lv.Eta())<1.4442";
		} else{
			cutStreamPhotons_EB << fCutStreamPhotons.str() ;
		}
		prediction->PhotonSigmaIEtaIEta_EB = new Channel("SigmaIEtaIEta_EB", "photon[0].SigmaIEtaIEta", 
								 cutStreamPhotons_EB.str().c_str(), fTriggerStream.str().c_str(),fSamplesRemovedPhotons);
		prediction->PhotonSigmaIEtaIEta_EB->fVerbose =prediction->fVerbose;
		prediction->PhotonSigmaIEtaIEta_EB->fOutputDir=prediction->fOutputDir;
		prediction->PhotonSigmaIEtaIEta_EB->fAddMCPedestal = true; // set this to avoid
		prediction->PhotonSigmaIEtaIEta_EB->fRootFile="SigmaIEtaIEta_EB_Shapes.root";
		prediction->PhotonSigmaIEtaIEta_EB->GetShapes("SigmaIEtaIEta_EB", "SigmaIEtaIEta", 30, 0, 0.08);
		prediction->GetPhotonNormalization(prediction->PhotonSigmaIEtaIEta_EB, prediction->MLRes_EB);

		// Get Photon Normalization: EE
		if(fSeparateEBEE){
			std::ostringstream cutStreamPhotons_EE;
			cutStreamPhotons_EE << fCutStreamPhotons.str() << "&&abs(photon[0].lv.Eta())>1.566";
			prediction->PhotonSigmaIEtaIEta_EE = new Channel("SigmaIEtaIEta_EE", "photon[0].SigmaIEtaIEta", 
									 cutStreamPhotons_EE.str().c_str(), fTriggerStream.str().c_str(),fSamplesRemovedPhotons);
			prediction->PhotonSigmaIEtaIEta_EE->fVerbose =prediction->fVerbose;
			prediction->PhotonSigmaIEtaIEta_EE->fOutputDir=prediction->fOutputDir;
			prediction->PhotonSigmaIEtaIEta_EE->fAddMCPedestal = true; // set this to avoid
			prediction->PhotonSigmaIEtaIEta_EE->fRootFile="SigmaIEtaIEta_EE_Shapes.root";
			prediction->PhotonSigmaIEtaIEta_EE->GetShapes("SigmaIEtaIEta_EE", "SigmaIEtaIEta", 30, 0, 0.08);
			prediction->GetPhotonNormalization(prediction->PhotonSigmaIEtaIEta_EE, prediction->MLRes_EE);
		}
	}

	// Photon Signal Region ******************************************************************************************
	if(fDoPhotonSignalRegion){	
		// Get Photon Selection Signal Region: EB
		std::ostringstream cutStreamPhotonsMT2_EB;
		cutStreamPhotonsMT2_EB << fCutStreamPhotonsMT2.str() << "&&abs(photon[0].lv.Eta())<1.4442";
		prediction->PhotonicSignalRegion_EB = new Channel("PhotonicSignalRegion_EB", "misc.MT2", cutStreamPhotonsMT2_EB.str().c_str(), 
							       fTriggerStream.str().c_str(), fSamplesRemovedPhotons);
		prediction->PhotonicSignalRegion_EB->fVerbose =prediction->fVerbose;
		prediction->PhotonicSignalRegion_EB->fOutputDir=prediction->fOutputDir;
		prediction->PhotonicSignalRegion_EB->fRootFile="SignalRegionRemovedPhotons_EB_Shapes.root";
		prediction->PhotonicSignalRegion_EB->GetShapes("PhotonicSignalRegion_EB", "MT2 (GeV)", 30, 0, 800);
		
		// Get Photon Selection Signal Region: EE
		std::ostringstream cutStreamPhotonsMT2_EE;
		cutStreamPhotonsMT2_EE << fCutStreamPhotonsMT2.str() << "&&abs(photon[0].lv.Eta())>1.566";
		prediction->PhotonicSignalRegion_EE = new Channel("PhotonicSignalRegion_EE", "misc.MT2", cutStreamPhotonsMT2_EE.str().c_str(), 
							       fTriggerStream.str().c_str(), fSamplesRemovedPhotons);
		prediction->PhotonicSignalRegion_EE->fVerbose =prediction->fVerbose;
		prediction->PhotonicSignalRegion_EE->fOutputDir=prediction->fOutputDir;
		prediction->PhotonicSignalRegion_EE->fRootFile="SignalRegionRemovedPhotons_EE_Shapes.root";
		prediction->PhotonicSignalRegion_EE->GetShapes("PhotonicSignalRegion_EE", "MT2 (GeV)", 30, 0, 800);
	}

	// Hadronic Signal Region ********************************************************************************************* 
	if(fDoHadronicSignalRegion){
		prediction->HadronicSignalRegion = new Channel("HadronicSignalRegion","misc.MT2", fCutStreamSignal.str().c_str(), 
							       fTriggerStream.str().c_str(), fSamplesHadronic);
		prediction->HadronicSignalRegion->fRootFile="HadronicMT2Shapes.root";
		prediction->HadronicSignalRegion->fVerbose =prediction->fVerbose;
		prediction->HadronicSignalRegion->fOutputDir=prediction->fOutputDir;
		prediction->HadronicSignalRegion->GetShapes("HadronicRegion", "MT2 (GeV)", 30, 0, 800);

		if(fDoPrediction)  {
		// compute MC Znunu/Photon ratio --------------------------------------------------
		prediction->GetMCZnunuToPhotonRatio();
		// make Prediction ----------------------------------------------------------------
		prediction->MakePrediction();
		}
	}
	

	if(fWriteToFile){
		TString logname =fOutDir + ".log"; 
		ofstream f_log (logname.Data(), ios::app);
		f_log << fLogStream->str();
	} else{
		cout << fLogStream->str();
	}
}


// ****************************************** Definition of classes and methods *********************************************************
class MLResult {
public:
	MLResult();
	~MLResult();
	float fMLNumPhotons;  
	float fMLNumPhotonsErr;
	float fMLPhotonScaleFactor;
	float fMLPhotonScaleFactorErr;
	float fMLNumQCD;
	float fMLNumQCDErr;
	float fMLQCDScaleFactor;
	float fMLQCDScaleFactorErr;
	float fMLNumOther;
	float fMLNumOtherErr;
	float fMLOtherScaleFactor;
	float fMLOtherScaleFactorErr;
	float fMLNumData;
	float fMLNumDataErr;
};
MLResult::MLResult(){};
MLResult::~MLResult(){};

class Channel {
public:
	Channel(TString name, TString variable, TString cuts, TString trigger, TString samples);
	~Channel();
	TH1D* hQCD; 
	TH1D* hData; 
	TH1D* hPhotons;
	TH1D* hWJets;
	TH1D* hZJetsToLL;
	TH1D* hZJetsToNuNu;
	TH1D* hTop;
	TH1D* hSignal;
	TH1D* hOther;
	TH1D* hTotalSM;
	TString fName;
	TString fCuts;
	TString fTrigger;
	TString fSamples;
	TString fVariable;
	TString fOutputDir;
	TString fRootFile;
	vector<TH1D*> h_shapes;
	void GetShapes(TString SelectionName, TString xtitle, int nbins, float binmin, float binmax);
	int fVerbose;
	bool fGotShapes;
	bool fAddMCPedestal;

};

Channel::Channel(TString name, TString variable, TString cuts, TString trigger, TString samples){
	fName=name;
	fCuts=cuts;
	fTrigger=trigger;
	fSamples=samples;
	fVariable=variable;
	hQCD=0;
	hData=0;
	hPhotons=0;
	hWJets=0;
	hZJetsToLL=0;
	hZJetsToNuNu=0;
	hTop=0;
	hSignal=0;
	hOther=0;
	hTotalSM=0;
	fOutputDir="test";
	fRootFile="test.root";
	fVerbose=4;
	fGotShapes=false;
	fAddMCPedestal=false;
}
Channel::~Channel(){};

void Channel::GetShapes(TString SelectionName, TString xtitle, int nbins, float binmin, float binmax){ 
	MT2Shapes *tA;
	if(fWriteToFile) tA = new MT2Shapes(fOutputDir, fRootFile, fLogStream);
	else             tA = new MT2Shapes(fOutputDir, fRootFile);
	tA->setVerbose(fVerbose);
	tA->init(fSamples);
	tA->SetPrintSummary(true);
	tA->SetDraw(false);
	tA->SetWrite(false);
  
//                    variable,      cuts,  njet, nlep,  selection_name,      HLT,   xtitle   nbins  bins   
        tA->GetShapes(fVariable,  fCuts,    -1,  -10  , SelectionName,    fTrigger , xtitle , nbins, binmin, binmax);

	// retrieve shapes
	for(int i=0; i<tA->GetNShapes(); ++i){
		TString name =tA->fh_shapes[i]->GetName();
		if      (name.Contains("QCD_"))        {hQCD         = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hQCD->SetDirectory(0);}
		else if (name.Contains("PhotonsJets_")){hPhotons     = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hPhotons->SetDirectory(0);}
		else if (name.Contains("Data_"))       {hData        = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hData->SetDirectory(0);}
		else if (name.Contains("WJets_"))      {hWJets       = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hWJets->SetDirectory(0);}
		else if (name.Contains("ZJetsToLL_"))  {hZJetsToLL   = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hZJetsToLL->SetDirectory(0);}
		else if (name.Contains("ZJetsToNuNu_")){hZJetsToNuNu = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hZJetsToNuNu->SetDirectory(0);}
		else if (name.Contains("Top_"))        {hTop         = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hTop->SetDirectory(0);}
		else if (name.Contains("Signal_"))     {hSignal      = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hSignal->SetDirectory(0);}
		else if (name.Contains("Other_"))      {hOther       = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hOther->SetDirectory(0);}
	}
	delete tA;
	fGotShapes=true;

	// fix colors
	if(hQCD!=0)        {hQCD    ->SetLineColor(kYellow+1);    hQCD    ->SetFillColor(kYellow+1);      hQCD    ->SetFillStyle(3001);}
	if(hPhotons!=0)    {hPhotons->SetLineColor(kViolet-3);    hPhotons->SetFillColor(kViolet-3);      hPhotons->SetFillStyle(3001);}
	if(hOther!=0)      {hOther  ->SetLineColor(kCyan+2);      hOther  ->SetFillColor(kCyan+2);        hOther  ->SetFillStyle(3001);}
	if(hZJetsToNuNu!=0){hZJetsToNuNu->SetLineColor(kGreen+1); hZJetsToNuNu->SetFillColor(kGreen+1);   hZJetsToNuNu  ->SetFillStyle(3001);}
	if(hSignal!=0)     {hSignal ->SetLineColor(kBlack);       hSignal  ->SetFillColor(kBlack);        hSignal  ->SetFillStyle(3001);}
	if(hTop!=0)        {hTop    ->SetLineColor(600);          hTop     ->SetFillColor(600);           hTop->SetFillStyle(3001);}

	if(fSaveResults){
		TString filename=fOutputDir+"/"+fRootFile;
		TFile *file = new TFile(filename.Data(), "UPDATE");
		if(hQCD!=0)         hQCD->Write();
		if(hPhotons!=0)     hPhotons->Write();
		if(hOther!=0)       hOther->Write();
		if(hData!=0)        hData->Write();
		if(hTop!=0)         hTop->Write();
		if(hZJetsToNuNu!=0) hZJetsToNuNu->Write();
		if(hSignal!=0)      hSignal->Write();
		file->Close();
	}
	if(fDraw){
		if(hQCD!=0)        DrawHisto(hQCD,           hQCD->GetName(),            "hist", this);
		if(hPhotons!=0)    DrawHisto(hPhotons,       hPhotons->GetName(),        "hist", this);
		if(hOther!=0)      DrawHisto(hOther,         hOther->GetName(),          "hist", this);
		if(hData!=0)       DrawHisto(hData,          hData->GetName(),           "EXO", this);
		if(hTop!=0)        DrawHisto(hTop,           hTop ->GetName(),           "hist", this);
		if(hZJetsToNuNu!=0)DrawHisto(hZJetsToNuNu,   hZJetsToNuNu->GetName(),    "hist", this);
		if(hSignal!=0)     DrawHisto(hSignal,        hSignal->GetName(),         "hist", this);
	}
}

class Prediction {
public:
	Prediction();
	~Prediction();
	
	void GetPhotonNormalization(Channel* channel, MLResult* MLRes);
	void GetMCZnunuToPhotonRatio();
	void MakePrediction();
	float GetZnunuPreditionErrorSys(TH1D* hData, TH1D* hOther, TH1D* hQCD);
	float GetZnunuPreditionErrorSysClosure(TH1D* hPhotons);
	float GetZnunuPreditionErrorStat(TH1D* hData);
	TH1D* GetScaledHisto(TH1D* histo, float scale_fact, float scale_fact_err);
	TH1D* GetScaledHisto(TH1D* histo, float scale_fact, float scale_fact_err, int ngroup);
	Channel* HadronicSignalRegion;
	Channel* PhotonicSignalRegion_EE;
	Channel* PhotonicSignalRegion_EB;
	Channel* PhotonSigmaIEtaIEta_EE;
	Channel* PhotonSigmaIEtaIEta_EB;
	MLResult* MLRes_EB;
	MLResult* MLRes_EE;

	bool fSave;
	int fVerbose;
	float fMCZnunuPhotonRatio;
	float fMCZnunuPhotonRatioErr;
	float fFractionPhotons;
	TString fOutputDir;

	
};

Prediction::Prediction(){
	fVerbose =0;
	fSave = false;
	fOutputDir = "./GammaJetsPrediction/";
	MLRes_EB = new MLResult();
	MLRes_EE = new MLResult();
}
Prediction::~Prediction(){}
void Prediction::GetPhotonNormalization(Channel* channel, MLResult* MLRes){
	if (channel->fGotShapes=false)                        {cerr << "GetPhotonNormalization: need to get Shapes first!"  << endl; exit(-1);}
	if (    channel->hQCD==0  || channel->hPhotons==0 
	     || channel->hData==0 || channel->hOther==0  )    {cerr << "GetPhotonNormalization: ERROR: received 0 pointer!" << endl; exit(-1);}
	*fLogStream <<"\n**************************************************************************************************\n"
	     <<"Starting to extract Photon Normalization.....\n";

	TH1D *hQCD     = (TH1D*)channel->hQCD->Clone("QCD");
	TH1D *hPhotons = (TH1D*)channel->hPhotons->Clone("Photons");
	TH1D *hData    = (TH1D*)channel->hData->Clone("Data");
	TH1D *hOther   = (TH1D*)channel->hOther->Clone("Other");

	if(channel->fAddMCPedestal){ // Add Pedestal to Other MC PDF in order to avoid bins in PDF with zero entries
		                     // but data in same bin! this causes problems!
		for(int i=0; i<hQCD->GetNbinsX(); ++i){
			if(hData->GetBinContent(i)==0) continue;
			double MCcontent = hPhotons->GetBinContent(i)+hQCD->GetBinContent(i)+hOther->GetBinContent(i);
			if(MCcontent==0) {hOther->SetBinContent(i, 1E-01);hOther->SetBinError(i, 1E-01);}
		}
	}
	
	RooRealVar sigmaietaieta("sigmaietaieta","sigmaietaieta",0.,0.08) ; // contained in histos

	RooDataHist Data   ("data"   ,"data"  ,sigmaietaieta,hData) ;    // define RooDataHists
	RooDataHist Photons("photons","photon",sigmaietaieta,hPhotons);
	RooDataHist QCD    ("QCD"    ,"QCD"   ,sigmaietaieta,hQCD);
	RooDataHist Other  ("QCD"    ,"QCD"   ,sigmaietaieta,hOther);

	RooHistPdf Photons_pdf("photons_pdf","photons_pdf",sigmaietaieta,Photons); // define PDFs for signal and bkg
	RooHistPdf QCD_pdf    ("qcd_pdf"    ,"qcd_pdf"    ,sigmaietaieta,QCD    ); 
	RooHistPdf Other_pdf  ("Other_pdf"  ,"other_pdf"  ,sigmaietaieta,Other  );

	RooRealVar nsig       ("nsig"   ,"number of signal events",     hPhotons->Integral()  ,  hPhotons->Integral()*0.1,hData->Integral());
	RooRealVar nqcd       ("nqcd"   ,"number of QCD events",        hQCD->Integral()      ,  hQCD->Integral()    *0.5,hData->Integral());
	RooRealVar nother     ("nother" ,"number of Other SM events",   hOther->Integral()); nother.setConstant(kTRUE);

	// model(x) = nsig*Photons_pdf(x) + nqcd*QCD_pdf(x) + nother*Other_pdf(x), where nother is fixed to nominal contribution
	RooAddPdf model("model","model", RooArgList(Photons_pdf,QCD_pdf,Other_pdf), RooArgList(nsig, nqcd, nother));
	model.defaultPrintStream(fLogStream);
	
	// preform fit
	RooFitResult* fitres = model.fitTo(Data,SumW2Error(kFALSE),Extended(), Save(kTRUE)); 
	// if I'm not mistaken: SumW2==false is the right option, as mc-histos already contain proper weights. 
	// SumW2Error == true would be needed if input comes from a TTree with (observable, weight) for each entry. 
	// then, data.setWeightVar(y) would also be needed. 


	// make plot
	TCanvas* canv = new TCanvas(channel->fName,"", 0, 0, 500, 500 );
	RooPlot* frame = sigmaietaieta.frame();
	Data.plotOn(frame, Name("Data")) ;
	model.plotOn(frame,Components(RooArgSet(Photons_pdf,QCD_pdf,Other_pdf)), Name("Model"));
	model.plotOn(frame,Components(QCD_pdf),        LineStyle(kDotted), LineColor(kMagenta));
	model.plotOn(frame,Components(Photons_pdf),    LineStyle(kDotted), LineColor(kGreen));
	frame->Draw();
	Double_t chi2 = frame->chiSquare("Model", "Data", 3);
	*fLogStream << "-----------------------------------------------------------------" << endl;
	*fLogStream << "Fit result for: " <<channel->fName                                 << endl; 
	fitres->Print("v");
	*fLogStream << "ChiSquare of fit: " << chi2                                        << endl;
	*fLogStream << "-----------------------------------------------------------------" << endl;
	
	// save RooFit output:
	if(fSave){
		TString filename=channel->fOutputDir+"/"+channel->fRootFile;
		TFile *file = new TFile(filename.Data(), "UPDATE");
		fitres->Write();
		frame->Write();
		file->Close();
	}

	// compute photon contribution and relative normaizations
	TH1D* dataclone = (TH1D*) hData->Clone("dataclone");
	dataclone->Rebin(dataclone->GetNbinsX());

	MLRes->fMLNumPhotons           = nsig.getVal();
	MLRes->fMLNumPhotonsErr        = nsig.getError();
	MLRes->fMLPhotonScaleFactor    = nsig.getVal()/channel->hPhotons->Integral();
	MLRes->fMLPhotonScaleFactorErr = nsig.getError()/channel->hPhotons->Integral();
	MLRes->fMLNumQCD               = nqcd.getVal();
	MLRes->fMLNumQCDErr            = nqcd.getError();
	MLRes->fMLQCDScaleFactor       = nqcd.getVal()/channel->hQCD->Integral();
	MLRes->fMLQCDScaleFactorErr    = nqcd.getError()/channel->hQCD->Integral();
	MLRes->fMLNumOther             = nother.getVal();
	MLRes->fMLNumOtherErr          = nother.getError();
	MLRes->fMLOtherScaleFactor     = nother.getVal()/channel->hOther->Integral();
	MLRes->fMLOtherScaleFactorErr  = nother.getError()/channel->hOther->Integral();
	MLRes->fMLNumData              = dataclone->GetBinContent(1);
	MLRes->fMLNumDataErr           = dataclone->GetBinError(1);

	delete dataclone;
	delete fitres;

	if(fVerbose>4){
		*fLogStream << "-------------------------------------------------------" << endl;
		*fLogStream << "fMLNumPhotons " <<  MLRes->fMLNumPhotons  << " pm " << MLRes->fMLNumPhotonsErr << endl; 
		*fLogStream << "fMLNumQCD "     <<  MLRes->fMLNumQCD      << " pm " << MLRes->fMLNumQCDErr << endl; 
		*fLogStream << "fMLNumOther "   <<  MLRes->fMLNumOther    << " pm " << MLRes->fMLNumOtherErr << endl; 
		*fLogStream << "fMLNumData "    <<  MLRes->fMLNumData     << " pm " << MLRes->fMLNumDataErr << endl; 
		*fLogStream << "-------------------------------------------------------" << endl;
	}
	if(fAddFitIntrinsicUncert){ // adding uncertainty to scale factors due to intrinsic method uncertainty
		*fLogStream << "+++ Adding a uncertainty of " << fFitIntrinsicUncert << " \% fit Scale Factors" << endl;
		MLRes->fMLPhotonScaleFactorErr = sqrt(pow(MLRes->fMLPhotonScaleFactorErr,2)+pow(fFitIntrinsicUncert*MLRes->fMLPhotonScaleFactor,2));
		MLRes->fMLQCDScaleFactorErr    = sqrt(pow(MLRes->fMLQCDScaleFactorErr   ,2)+pow(fFitIntrinsicUncert*MLRes->fMLQCDScaleFactor   ,2));
	}
	
	if(fVerbose>4){
		*fLogStream << "-------------------------------------------------------" << endl;
		*fLogStream << "Photon Scale factor: " << MLRes->fMLPhotonScaleFactor << " pm " << MLRes->fMLPhotonScaleFactorErr << endl;
		*fLogStream << "QCD    Scale factor: " << MLRes->fMLQCDScaleFactor    << " pm " << MLRes->fMLQCDScaleFactorErr    << endl;
		*fLogStream << "Other  Scale factor: " << MLRes->fMLOtherScaleFactor  << " pm " << MLRes->fMLOtherScaleFactorErr  << endl;
		*fLogStream << "-------------------------------------------------------" << endl;
	}

}

void Prediction::GetMCZnunuToPhotonRatio(){
	// compute MC Znunu/Photon ratio
	// for this PhotonicSignalRegion->fGotShapes==1 and HadronicSignalRegion->fGotShapes==1 is required!
	if(PhotonicSignalRegion_EE->fGotShapes==false || PhotonicSignalRegion_EB->fGotShapes==false || HadronicSignalRegion->fGotShapes==false){
		cerr << "ERROR in GetMCZnunuToPhotonRatio: fGotShape found to be false" << endl;
		exit(-1);
	}
	*fLogStream << "------------------------------------------------------------------------" << endl;
	*fLogStream << "GetMCZnunuToPhotonRatio: Computing MC Znunu to Photon Ratio             " << endl;           

	// GetScaled histos with only one bin and probagated errors: stat error and error on scale factor
	TH1D *currPhotons= GetScaledHisto(PhotonicSignalRegion_EB->hPhotons , MLRes_EB->fMLPhotonScaleFactor, MLRes_EB->fMLPhotonScaleFactorErr); // EB
	currPhotons->Add(  GetScaledHisto(PhotonicSignalRegion_EE->hPhotons , MLRes_EE->fMLPhotonScaleFactor, MLRes_EE->fMLPhotonScaleFactorErr)); //EE
	TH1D *currZnunu  = GetScaledHisto(HadronicSignalRegion->hZJetsToNuNu, 1                   , 0);

	TH1D *ratio = (TH1D*) currZnunu->Clone("Znunu_To_Photon_ratio");
	ratio->Divide(currPhotons);
	fMCZnunuPhotonRatio    = ratio->GetBinContent(1);
	fMCZnunuPhotonRatioErr = ratio->GetBinError(1);

	*fLogStream << "Ratio: " << fMCZnunuPhotonRatio << " pm " << fMCZnunuPhotonRatioErr << endl;
	*fLogStream << "------------------------------------------------------------------------" << endl;

	if(fSaveZnunuToGammaRatio){
		TH1D *currPhotons2= GetScaledHisto(PhotonicSignalRegion_EB->hPhotons , MLRes_EB->fMLPhotonScaleFactor, MLRes_EB->fMLPhotonScaleFactorErr, 2); // EB
		currPhotons2->Add(  GetScaledHisto(PhotonicSignalRegion_EE->hPhotons , MLRes_EE->fMLPhotonScaleFactor, MLRes_EE->fMLPhotonScaleFactorErr, 2)); //EE
		TH1D *currZnunu2  = GetScaledHisto(HadronicSignalRegion->hZJetsToNuNu, 1                               , 0                              , 2);
		TH1D *ratio2 = (TH1D*) currZnunu2->Clone("Znunu_To_Photon_Ratio");
		ratio2->Divide(currPhotons2);
		TString filename=fOutputDir+"/ZnunuToGammaRatio.root";
		TFile *file = new TFile(filename.Data(), "RECREATE");
		currZnunu2->Write();
		currPhotons2->Write();
		ratio2->Write();
		file->Close();
		delete currPhotons2;
		delete currZnunu2;
		delete ratio2;
	}
	delete currPhotons;
	delete currZnunu;
	delete ratio;
}

void Prediction::MakePrediction(){
	if(fMCZnunuPhotonRatio==0)        {*fLogStream << "ERROR in MakePrediction: fMCZnunuPhotonRatio==0"        << endl; exit(-1); }
	*fLogStream << "******************************* Prediction ****************************************" << endl;

	*fLogStream << "Photonic Region where EB Normalization was extracted: (ML-fit result) ------------" << endl;
	*fLogStream << "  NData: "    << MLRes_EB->fMLNumData    << " pm " << MLRes_EB->fMLNumDataErr       << endl;
	*fLogStream << "  NPhotons: " << MLRes_EB->fMLNumPhotons << " pm " << MLRes_EB->fMLNumPhotonsErr    << endl;
	*fLogStream << "  NOther: "   << MLRes_EB->fMLNumOther   << " pm " << MLRes_EB->fMLNumOtherErr      << endl;
	*fLogStream << "  NQCD: "     << MLRes_EB->fMLNumQCD	  << " pm " << MLRes_EB->fMLNumQCDErr        << endl;
	*fLogStream << "Photonic Region where EE Normalization was extracted: (ML-fit result) ------------" << endl;
	*fLogStream << "  NData: "    << MLRes_EE->fMLNumData    << " pm " << MLRes_EE->fMLNumDataErr       << endl;
	*fLogStream << "  NPhotons: " << MLRes_EE->fMLNumPhotons << " pm " << MLRes_EE->fMLNumPhotonsErr    << endl;
	*fLogStream << "  NOther: "   << MLRes_EE->fMLNumOther   << " pm " << MLRes_EE->fMLNumOtherErr      << endl;
	*fLogStream << "  NQCD: "     << MLRes_EE->fMLNumQCD	  << " pm " << MLRes_EE->fMLNumQCDErr        << endl;
	
	// scale MC event yields in PhotonicSignalRegion to data with MC scale factor as measured in 
	// GetPhotonNormalization()
	
	// Get Scales 1-bin histograms!
	// EB
	TH1D *hPhotons_EB = GetScaledHisto(PhotonicSignalRegion_EB->hPhotons, MLRes_EB->fMLPhotonScaleFactor, MLRes_EB->fMLPhotonScaleFactorErr);
	TH1D *hOther_EB   = GetScaledHisto(PhotonicSignalRegion_EB->hOther  , MLRes_EB->fMLOtherScaleFactor , MLRes_EB->fMLOtherScaleFactorErr);
	TH1D *hQCD_EB     = GetScaledHisto(PhotonicSignalRegion_EB->hQCD    , MLRes_EB->fMLQCDScaleFactor   , MLRes_EB->fMLQCDScaleFactorErr);
	TH1D *hData_EB    = GetScaledHisto(PhotonicSignalRegion_EB->hData   , 1                   , 0);
	// EE
	TH1D *hPhotons_EE = GetScaledHisto(PhotonicSignalRegion_EE->hPhotons, MLRes_EE->fMLPhotonScaleFactor, MLRes_EE->fMLPhotonScaleFactorErr);
	TH1D *hOther_EE   = GetScaledHisto(PhotonicSignalRegion_EE->hOther  , MLRes_EE->fMLOtherScaleFactor , MLRes_EE->fMLOtherScaleFactorErr);
	TH1D *hQCD_EE     = GetScaledHisto(PhotonicSignalRegion_EE->hQCD    , MLRes_EE->fMLQCDScaleFactor   , MLRes_EE->fMLQCDScaleFactorErr);
	TH1D *hData_EE    = GetScaledHisto(PhotonicSignalRegion_EE->hData   , 1                   , 0);
	
	TH1D* hPhotons = hPhotons_EB->Clone("hPhotons"); hPhotons->Add(hPhotons_EE);
	TH1D* hOther   = hOther_EB->Clone("hOther");     hOther  ->Add(hOther_EE);
	TH1D* hQCD     = hQCD_EB->Clone("hQCD");         hQCD    ->Add(hQCD_EE);
	TH1D* hData    = hData_EB->Clone("hData");       hData   ->Add(hData_EE);


	*fLogStream << "Photonic Signal Region: event yields:-------------------" << endl;
	*fLogStream << "EB:                                                     " << endl;
	*fLogStream << "  NData:    "    << hData_EB->GetBinContent(1)    << " pm " << hData_EB   ->GetBinError(1)  << endl;
	*fLogStream << "  NPhotons: "    << hPhotons_EB->GetBinContent(1) << " pm " << hPhotons_EB->GetBinError(1)  << endl;
	*fLogStream << "  NQCD:     "    << hQCD_EB->GetBinContent(1)     << " pm " << hQCD_EB    ->GetBinError(1)  << endl;
	*fLogStream << "  NOther:   "    << hOther_EB->GetBinContent(1)   << " pm " << hOther_EB  ->GetBinError(1)  << endl;
	*fLogStream << "EE:                                                     " << endl;
	*fLogStream << "  NData:    "    << hData_EE->GetBinContent(1)    << " pm " << hData_EE   ->GetBinError(1)  << endl;
	*fLogStream << "  NPhotons: "    << hPhotons_EE->GetBinContent(1) << " pm " << hPhotons_EE->GetBinError(1)  << endl;
	*fLogStream << "  NQCD:     "    << hQCD_EE->GetBinContent(1)     << " pm " << hQCD_EE    ->GetBinError(1)  << endl;
	*fLogStream << "  NOther:   "    << hOther_EE->GetBinContent(1)   << " pm " << hOther_EE  ->GetBinError(1)  << endl;
	*fLogStream << "total:                                                  " << endl;
	*fLogStream << "  NData:    "    << hData->GetBinContent(1)    << " pm " << hData   ->GetBinError(1)  << endl;
	*fLogStream << "  NPhotons: "    << hPhotons->GetBinContent(1) << " pm " << hPhotons->GetBinError(1)  << endl;
	*fLogStream << "  NQCD:     "    << hQCD->GetBinContent(1)     << " pm " << hQCD    ->GetBinError(1)  << endl;
	*fLogStream << "  NOther:   "    << hOther->GetBinContent(1)   << " pm " << hOther  ->GetBinError(1)  << endl;
	*fLogStream << "where the following scale factors were used:            " << endl;
	*fLogStream << "EB:                                                     " << endl;
	*fLogStream << "  Photon Scale Factor: " << MLRes_EB->fMLPhotonScaleFactor << " pm " << MLRes_EB->fMLPhotonScaleFactorErr            << endl;
	*fLogStream << "  QCD    Scale Factor: " << MLRes_EB->fMLQCDScaleFactor    << " pm " << MLRes_EB->fMLQCDScaleFactorErr               << endl;
	*fLogStream << "  Other  Scake Factor: " << MLRes_EB->fMLOtherScaleFactor  << " pm " << MLRes_EB->fMLOtherScaleFactorErr             << endl;
	*fLogStream << "EE:                                                     " << endl;
	*fLogStream << "  Photon Scale Factor: " << MLRes_EE->fMLPhotonScaleFactor << " pm " << MLRes_EE->fMLPhotonScaleFactorErr            << endl;
	*fLogStream << "  QCD    Scale Factor: " << MLRes_EE->fMLQCDScaleFactor    << " pm " << MLRes_EE->fMLQCDScaleFactorErr               << endl;
	*fLogStream << "  Other  Scake Factor: " << MLRes_EE->fMLOtherScaleFactor  << " pm " << MLRes_EE->fMLOtherScaleFactorErr             << endl;
	
	*fLogStream << "MC Znunu To Photon ratio ----------------------------------- " << endl;
	*fLogStream << fMCZnunuPhotonRatio  << " pm " << fMCZnunuPhotonRatioErr        << endl;	
	
	if(fMT2b && fDoBRatioCorrection){
	*fLogStream << "MC to Data Photon+jets B-tag correction: " << fRB_MCtoData << " pm " << fRB_MCtoDataErr << endl;
	}

	float PredictedZnunu                = fRB_MCtoData*fMCZnunuPhotonRatio*(hData->GetBinContent(1)-hOther->GetBinContent(1)-hQCD->GetBinContent(1));


	float PredictedZnunu_ErrSys         = GetZnunuPreditionErrorSys (hData, hOther, hQCD);
	float PredictedZnunu_ErrSysClosure  = GetZnunuPreditionErrorSysClosure (hPhotons);
	float PredictedZnunu_ErrStat        = GetZnunuPreditionErrorStat(hData);
	*fLogStream << "Predicted N Znunu events in Hadronic Signal region           "                                           << endl;
	*fLogStream << "Prediction: "                                                                                            << endl;
	*fLogStream << "  " << PredictedZnunu << " pm " << PredictedZnunu_ErrSys  << " pm " << PredictedZnunu_ErrStat << " stat "        << endl;
	*fLogStream << "True N Znunu events:"                                                                                    << endl;
	*fLogStream << "  " << HadronicSignalRegion->hZJetsToNuNu->Integral()                                                            << endl;
	*fLogStream << "MC cloure:"                                                                                              << endl;
	*fLogStream << "  " << fMCZnunuPhotonRatio*(hPhotons->GetBinContent(1))   << " pm " << PredictedZnunu_ErrSysClosure << " (sys) " << endl; 

}

TH1D* Prediction::GetScaledHisto(TH1D* histo, float scale_fact, float scale_fact_err){
	// takes histo, scale factor and uncertainty on scale factor
	// and returns scaled and rebinned histo with 1 bin with uncertainty on scale factor propagated.
	TH1D *h        = (TH1D*) histo->Clone(histo->GetName());
	h->Rebin(h->GetNbinsX());
	h->SetBinError(1, sqrt(h->GetBinError(1)*  h->GetBinError(1)   *scale_fact    *scale_fact + 
			  h->GetBinContent(1)*h->GetBinContent(1) *scale_fact_err*scale_fact_err));
	h->SetBinContent(1, h->GetBinContent(1)*scale_fact);
	return h;
}

TH1D* Prediction::GetScaledHisto(TH1D* histo, float scale_fact, float scale_fact_err, int ngroup){
	// takes histo, scale factor and uncertainty on scale factor
	// and returns scaled and rebinned histo with ngroup bins merged into 1 with uncertainty on scale factor propagated.
	if(ngroup>=histo->GetNbinsX()) ngroup=histo->GetNbinsX();
	TH1D *h        = (TH1D*) histo->Clone(histo->GetName());
	h->Rebin(ngroup);
	for(int i=0; i<h->GetNbinsX(); ++i){
		h->SetBinError(i, sqrt(h->GetBinError(i)*  h->GetBinError(i)   *scale_fact    *scale_fact + 
				  h->GetBinContent(i)*h->GetBinContent(i) *scale_fact_err*scale_fact_err));
		h->SetBinContent(i, h->GetBinContent(i)*scale_fact);
	}
	return h;
}

float Prediction::GetZnunuPreditionErrorSys(TH1D* hData, TH1D* hOther, TH1D* hQCD){
	// histograms must be scales with prober uncertainties
	float pred_err_2 = pow(fMCZnunuPhotonRatioErr * (hData->GetBinContent(1) - hOther->GetBinContent(1) - hQCD->GetBinContent(1))*fRB_MCtoData,2);
	pred_err_2      += pow(fMCZnunuPhotonRatio*fRB_MCtoData * hOther->GetBinError(1),2);
	pred_err_2      += pow(fMCZnunuPhotonRatio*fRB_MCtoData * hQCD->GetBinError(1),2);
	pred_err_2      += pow(fRB_MCtoDataErr* fMCZnunuPhotonRatio*(hData->GetBinContent(1) - hOther->GetBinContent(1) - hQCD->GetBinContent(1)),2);

	return sqrt(pred_err_2);
}

float Prediction::GetZnunuPreditionErrorSysClosure(TH1D* hPhotons){
	// computes systematic uncerainty on Znunu prediction for MC closure test:
	// takes into account MLfit uncertainty on Photon-MC yield, MC statistics, statistical uncerainty on Znunu/photon ratio
	float pred_err_2 = pow(fMCZnunuPhotonRatioErr* hPhotons->GetBinContent(1)*fRB_MCtoData,2);
	pred_err_2      += pow(fMCZnunuPhotonRatio   * hPhotons->GetBinError(1)  *fRB_MCtoData,2);
	pred_err_2      += pow(fMCZnunuPhotonRatio   * hPhotons->GetBinContent(1)*fRB_MCtoDataErr,2);
	return sqrt(pred_err_2);
}

float Prediction::GetZnunuPreditionErrorStat(TH1D* hData){
	// statistical uncertainty on prediction due to data statistics
	float pred_err   = fRB_MCtoData* fMCZnunuPhotonRatio * hData->GetBinError(1);
	return pred_err;
}

void DrawHisto(TH1* h_orig, TString name,  Option_t *drawopt, Channel* channel){
	TH1D* h = (TH1D*)h_orig->Clone(name);
	TString canvname = "canv_"+name;
	TCanvas *col = new TCanvas(canvname, "", 0, 0, 500, 500);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	gStyle->SetPalette(1);
	gPad->SetRightMargin(0.15);
	h->DrawCopy(drawopt);
	gPad->RedrawAxis();
	if(fSaveResults){
		TString filename=channel->fOutputDir+"/"+channel->fRootFile;
		TFile *file = new TFile(filename.Data(), "UPDATE");
		col->Write();
		file->Close();
	}
}
