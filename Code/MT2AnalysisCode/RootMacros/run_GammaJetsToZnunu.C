
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
//#include "../MT2Code/include/MT2Shapes.hh"
using namespace std;
using namespace RooFit ;

// User Input:  ----------------------------------------------
TString fSamplesRemovedPhotons     ="./samples/samples_RemovedPhotons.dat";
//TString fSamplesRemovedPhotons   ="./samples/samples_RemovedPhotonsHerwig.dat";
TString fSamplesHadronic           ="./samples/samples_highMT2noQCD.dat";
Bool_t  fDoPhotonSignalRegion      = true;
Bool_t  fDoHadronicSignalRegion    = true;
Bool_t  fDoPrediction              = true;
Int_t   fVerbose                   = 1;
Bool_t  fSaveMLFitResult           = true;
Bool_t  fMT2b                      = true;

// Global Variables ------------------------------------
std::ostringstream fTriggerStream;
std::ostringstream fCutStreamPhotons;
std::ostringstream fCutStreamPhotonsMT2;
std::ostringstream fCutStreamSignal;

// Cut Streams
void DefineCutStreams(){
	// Trigger Stream ---------------------------------------------------------------
	fTriggerStream << "( "
		<< "(trigger.HLT_HT440_v2==1 && misc.Run<=161176)" << "||"
		<< "(trigger.HLT_HT450_v2==1 && (misc.Run>=161217 && misc.Run<=163261))" << "||"
		<< "(trigger.HLT_HT500_v3==1 && (misc.Run>=163270 && misc.Run<=163869))" << "||"
		<< "(trigger.HLT_HT550_v4==1 && (misc.Run>=165088 && misc.Run<=165633))" << "||"
		<< "(trigger.HLT_HT550_v5==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" << "||"
		<< "(trigger.HLT_HT550_v6==1 && (misc.Run==166346))"                     << "||"
		<< "(trigger.HLT_HT550_v7==1 && (misc.Run>=167078 && misc.Run<=167913))" << "||"
		<< "(trigger.HLT_HT550_v8==1 && (misc.Run>=169561 && misc.Run<=173198))" << "||"
		<< "(trigger.HLT_HT600_v1==1 && (misc.Run>=173236 && misc.Run<=173692))" << " )";

	// CutStream for SigmaIEtaIEta ------------------------------------------
	fCutStreamPhotons << " " 
	  << "misc.MET>=30"                                                << "&&"
	  << "misc.HT > 700 "                                              << "&&"
	  << "NJetsIDLoose >=3"                                            << "&&"
	  << "(NEles==0  || ele[0].lv.Pt()<10)"                            << "&&"
	  << "(NMuons==0 || muo[0].lv.Pt()<10)"                            << "&&"
	  << "misc.Jet0Pass ==1"                                           << "&&"
	  << "misc.Jet1Pass ==1"                                           << "&&"
	  << "misc.SecondJPt >100"                                         << "&&"
	  << "misc.PassJetID ==1"                                          << "&&"
	  << "misc.Vectorsumpt < 70"                                       << "&&"
	  << "misc.MinMetJetDPhi >0.3"                                     << "&&"
	  << "NPhotons ==1 "                                               << "&&"
	  << "(misc.ProcessID!=6||photon[0].MCmatchexitcode!=1)"           << "&&"
	  << "rawpfmet[0].Pt()<50"                                         << "&&"
	  << "misc.HBHENoiseFlag ==0"                                      << "&&"
	  << "misc.CSCTightHaloID==0"                                      << "&&"
	  << "misc.CrazyHCAL==0";

	// CutStream for Photon Signal Region ------------------------------------------ 
	fCutStreamPhotonsMT2 << " " 
	  << "misc.MET>=30"                                                   << "&&"
	  << "misc.HT > 700 "                                                 << "&&"
	  << "misc.MT2 > 200"                                                 << "&&"
	  << "NJetsIDLoose >=3"                                               << "&&"
	  << "(NEles==0  || ele[0].lv.Pt()<10)"                               << "&&"
	  << "(NMuons==0 || muo[0].lv.Pt()<10)"                               << "&&"
	  << "misc.Jet0Pass ==1"                                              << "&&"
	  << "misc.Jet1Pass ==1"                                              << "&&"
	  << "misc.SecondJPt >100"                                            << "&&"
	  << "misc.PassJetID ==1"                                             << "&&"
	  << "misc.Vectorsumpt < 70"                                          << "&&"
	  << "misc.MinMetJetDPhi >0.3"                                        << "&&"
	  << "NPhotons ==1 "                                                  << "&&"
	  << "photon[0].isEGMloose==1"                                        << "&&"
	  << "(misc.ProcessID!=6||photon[0].MCmatchexitcode!=1)"              << "&&"
	  << "rawpfmet[0].Pt()<50"                                            << "&&"
	  << "misc.HBHENoiseFlag ==0"                                         << "&&"
	  << "misc.CSCTightHaloID==0"                                         << "&&"
	  << "misc.CrazyHCAL==0";

	// CutStream for Hadronic Signal Region ----------------------------------------------
	fCutStreamSignal << " " 
	  << "misc.MET>=30"                            << "&&"
	  << "misc.HT > 700 "                          << "&&"
	  << "misc.MT2 > 200"                          << "&&"
	  << "NJetsIDLoose >=3"                        << "&&"
	  << "(NEles==0  || ele[0].lv.Pt()<10)"        << "&&"
	  << "(NMuons==0 || muo[0].lv.Pt()<10)"        << "&&"
	  << "misc.Jet0Pass ==1"                       << "&&"
	  << "misc.Jet1Pass ==1"                       << "&&"
	  << "misc.SecondJPt  >100"                    << "&&"
	  << "misc.PassJetID ==1"                      << "&&"
	  << "misc.Vectorsumpt < 70"                   << "&&"
	  << "misc.MinMetJetDPhi >0.3"                 << "&&"
	  << "misc.HBHENoiseFlag ==0"                  << "&&"
	  << "misc.CSCTightHaloID==0"                  << "&&"
	  << "misc.CrazyHCAL==0";
}

// *********************************** run_GammaJetsToZnunu: specify user input here *************************************************
void run_GammaJetsToZnunu(){
	gSystem->Load("libPhysics");
	gSystem->Load("libRooFit") ;
	gSystem->CompileMacro("../MT2Code/src/MT2Shapes.cc", "k");
	
	// define cutsteams
	DefineCutStreams();

	// new prediction ------------------------------------------------------
	Prediction* prediction = new Prediction();
	prediction->fVerbose=fVerbose;
	prediction->fSave   =fSaveMLFitResult;

	// Get Photon Normalization: EB
	std::ostringstream cutStreamPhotons_EB;
       	cutStreamPhotons_EB << fCutStreamPhotons.str() << "&&abs(photon[0].lv.Eta())<1.4442";
	prediction->PhotonSigmaIEtaIEta_EB = new Channel("SigmaIEtaIEta_EB", "photon[0].SigmaIEtaIEta", 
			                                 cutStreamPhotons_EB.str().c_str(), fTriggerStream.str().c_str(),fSamplesRemovedPhotons);
	prediction->PhotonSigmaIEtaIEta_EB->fVerbose =prediction->fVerbose;
	prediction->PhotonSigmaIEtaIEta_EB->fRootFile="SigmaIEtaIEta_EB_Shapes.root";
	prediction->PhotonSigmaIEtaIEta_EB->GetShapes("SigmaIEtaIEta_EB", "SigmaIEtaIEta", 30, 0, 0.08);
	prediction->GetPhotonNormalization(prediction->PhotonSigmaIEtaIEta_EB, prediction->MLRes_EB);

	// Get Photon Normalization: EE
	std::ostringstream cutStreamPhotons_EE;
	cutStreamPhotons_EE << fCutStreamPhotons.str() << "&&abs(photon[0].lv.Eta())>1.566";
	prediction->PhotonSigmaIEtaIEta_EE = new Channel("SigmaIEtaIEta_EE", "photon[0].SigmaIEtaIEta", 
			                                 cutStreamPhotons_EE.str().c_str(), fTriggerStream.str().c_str(),fSamplesRemovedPhotons);
	prediction->PhotonSigmaIEtaIEta_EE->fVerbose =prediction->fVerbose;
	prediction->PhotonSigmaIEtaIEta_EE->fRootFile="SigmaIEtaIEta_EE_Shapes.root";
	prediction->PhotonSigmaIEtaIEta_EE->GetShapes("SigmaIEtaIEta_EE", "SigmaIEtaIEta", 30, 0, 0.08);
	prediction->GetPhotonNormalization(prediction->PhotonSigmaIEtaIEta_EE, prediction->MLRes_EE);

	// Photon Signal Region ******************************************************************************************
	if(fDoPhotonSignalRegion){	
		// Get Photon Selection Signal Region: EB
		std::ostringstream cutStreamPhotonsMT2_EB;
		cutStreamPhotonsMT2_EB << fCutStreamPhotonsMT2.str() << "&&abs(photon[0].lv.Eta())<1.4442";
		prediction->PhotonicSignalRegion_EB = new Channel("PhotonicSignalRegion_EB", "misc.MT2", cutStreamPhotonsMT2_EB.str().c_str(), 
							       fTriggerStream.str().c_str(), fSamplesRemovedPhotons);
		prediction->PhotonicSignalRegion_EB->fVerbose =prediction->fVerbose;
		prediction->PhotonicSignalRegion_EB->fRootFile="SignalRegionRemovedPhotons_EB_Shapes.root";
		prediction->PhotonicSignalRegion_EB->GetShapes("PhotonicSignalRegion_EB", "MT2 (GeV)", 30, 0, 800);
		
		// Get Photon Selection Signal Region: EE
		std::ostringstream cutStreamPhotonsMT2_EE;
		cutStreamPhotonsMT2_EE << fCutStreamPhotonsMT2.str() << "&&abs(photon[0].lv.Eta())>1.566";
		prediction->PhotonicSignalRegion_EE = new Channel("PhotonicSignalRegion_EE", "misc.MT2", cutStreamPhotonsMT2_EE.str().c_str(), 
							       fTriggerStream.str().c_str(), fSamplesRemovedPhotons);
		prediction->PhotonicSignalRegion_EE->fVerbose =prediction->fVerbose;
		prediction->PhotonicSignalRegion_EE->fRootFile="SignalRegionRemovedPhotons_EE_Shapes.root";
		prediction->PhotonicSignalRegion_EE->GetShapes("PhotonicSignalRegion_EE", "MT2 (GeV)", 30, 0, 800);
	}

	// Hadronic Signal Region ********************************************************************************************* 
	if(fDoHadronicSignalRegion){
		prediction->HadronicSignalRegion = new Channel("HadronicSignalRegion","misc.MT2", fCutStreamSignal.str().c_str(), 
							       fTriggerStream.str().c_str(), fSamplesHadronic);
		prediction->HadronicSignalRegion->fRootFile="HadronicMT2Shapes.root";
		prediction->HadronicSignalRegion->fVerbose =prediction->fVerbose;
		prediction->HadronicSignalRegion->GetShapes("HadronicRegion", "MT2 (GeV)", 30, 0, 800);

		if(fDoPrediction)  {
		// compute MC Znunu/Photon ratio --------------------------------------------------
		prediction->GetMCZnunuToPhotonRatio();
		// make Prediction ----------------------------------------------------------------
		prediction->MakePrediction();
		}
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
	fOutputDir="GammaJetsPrediction";
	fRootFile="test.root";
	fVerbose=4;
	fGotShapes=false;
}
Channel::~Channel(){};

void Channel::GetShapes(TString SelectionName, TString xtitle, int nbins, float binmin, float binmax){ 
	MT2Shapes *tA = new MT2Shapes(fOutputDir, fRootFile);
	tA->setVerbose(fVerbose);
	tA->init(fSamples);
	tA->SetPrintSummary(true);
  
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
	
};

Prediction::Prediction(){
	fVerbose =0;
	fSave = false;
	MLRes_EB = new MLResult();
	MLRes_EE = new MLResult();
}
Prediction::~Prediction(){}
void Prediction::GetPhotonNormalization(Channel* channel, MLResult* MLRes){
	if (channel->fGotShapes=false)                        {cerr << "GetPhotonNormalization: need to get Shapes first!"  << endl; exit(-1);}
	if (    channel->hQCD==0  || channel->hPhotons==0 
	     || channel->hData==0 || channel->hOther==0  )    {cerr << "GetPhotonNormalization: ERROR: received 0 pointer!" << endl; exit(-1);}
	cout <<"\n**************************************************************************************************\n"
	     <<"Starting to extract Photon Normalization.....\n";

	TH1D *hQCD     = (TH1D*)channel->hQCD->Clone("QCD");
	TH1D *hPhotons = (TH1D*)channel->hPhotons->Clone("Photons");
	TH1D *hData    = (TH1D*)channel->hData->Clone("Data");
	TH1D *hOther   = (TH1D*)channel->hOther->Clone("Other");

	RooRealVar sigmaietaieta("sigmaietaieta","sigmaietaieta",0.,0.08) ; // contained in histos

	RooDataHist Data   ("data"   ,"data"  ,sigmaietaieta,hData) ;    // define RooDataHists
	RooDataHist Photons("photons","photon",sigmaietaieta,hPhotons);
	RooDataHist QCD    ("QCD"    ,"QCD"   ,sigmaietaieta,hQCD);
	RooDataHist Other  ("QCD"    ,"QCD"   ,sigmaietaieta,hOther);

	RooHistPdf Photons_pdf("photons_pdf","photons_pdf",sigmaietaieta,Photons); // define PDFs for signal and bkg
	RooHistPdf QCD_pdf    ("qcd_pdf"    ,"qcd_pdf"    ,sigmaietaieta,QCD    ); 
	RooHistPdf Other_pdf  ("Other_pdf"  ,"other_pdf"  ,sigmaietaieta,Other  );

	RooRealVar nsig       ("nsig"   ,"number of signal events",     hPhotons->Integral()  ,  0.,hData->Integral());
	RooRealVar nqcd       ("nqcd"   ,"number of QCD events",        hQCD->Integral()      ,  0.,hData->Integral());
	RooRealVar nother     ("nother" ,"number of Other SM events",   hOther->Integral()); nother.setConstant(kTRUE);

	// model(x) = nsig*Photons_pdf(x) + nqcd*QCD_pdf(x) + nother*Other_pdf(x), where nother is fixed to nominal contribution
	RooAddPdf model("model","model", RooArgList(Photons_pdf,QCD_pdf,Other_pdf), RooArgList(nsig, nqcd, nother));
	
	// preform fit
	RooFitResult* fitres = model.fitTo(Data,SumW2Error(kFALSE),Extended(), Save(kTRUE)); 
	// if I'm not mistaken: SumW2==false is the right option, as mc-histos already contain proper weights. 
	// SumW2Error == true would be needed if input comes from a TTree with (observable, weight) for each entry. 
	// then, data.setWeightVar(y) would also be needed. 


	// make plot
	RooPlot* frame = sigmaietaieta.frame();
	Data.plotOn(frame, Name("Data")) ;
	model.plotOn(frame,Components(RooArgSet(Photons_pdf,QCD_pdf,Other_pdf)), Name("Model"));
	model.plotOn(frame,Components(QCD_pdf),        LineStyle(kDotted), LineColor(kMagenta));
	model.plotOn(frame,Components(Photons_pdf),    LineStyle(kDotted), LineColor(kGreen));
	frame->Draw();
	Double_t chi2 = frame->chiSquare("Model", "Data", 3);
	cout << "-----------------------------------------------------------------" << endl;
	cout << "Fit result for: " <<channel->fName                                 << endl; 
	fitres->Print("v");
	cout << "ChiSquare of fit: " << chi2                                        << endl;
	cout << "-----------------------------------------------------------------" << endl;
	
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
	
	if(fVerbose>4){
		cout << "fMLNumPhotons " <<  MLRes->fMLNumPhotons  << " pm " << MLRes->fMLNumPhotonsErr << endl; 
		cout << "fMLNumQCD "     <<  MLRes->fMLNumQCD      << " pm " << MLRes->fMLNumQCDErr << endl; 
		cout << "fMLNumOther "   <<  MLRes->fMLNumOther    << " pm " << MLRes->fMLNumOtherErr << endl; 
		cout << "fMLNumData "    <<  MLRes->fMLNumData     << " pm " << MLRes->fMLNumDataErr << endl; 
		cout << "-------------------------------------------------------" << endl;
		cout << "Photon Scale factor: " << MLRes->fMLPhotonScaleFactor << " pm " << MLRes->fMLPhotonScaleFactorErr << endl;
		cout << "QCD    Scale factor: " << MLRes->fMLQCDScaleFactor    << " pm " << MLRes->fMLQCDScaleFactorErr    << endl;
		cout << "Other  Scale factor: " << MLRes->fMLOtherScaleFactor  << " pm " << MLRes->fMLOtherScaleFactorErr  << endl;
		cout << "-------------------------------------------------------" << endl;
	}

}

void Prediction::GetMCZnunuToPhotonRatio(){
	// compute MC Znunu/Photon ratio
	// for this PhotonicSignalRegion->fGotShapes==1 and HadronicSignalRegion->fGotShapes==1 is required!
	if(PhotonicSignalRegion_EE->fGotShapes==false || PhotonicSignalRegion_EB->fGotShapes==false || HadronicSignalRegion->fGotShapes==false){
		cerr << "ERROR in GetMCZnunuToPhotonRatio: fGotShape found to be false" << endl;
		exit(-1);
	}
	cout << "------------------------------------------------------------------------" << endl;
	cout << "GetMCZnunuToPhotonRatio: Computing MC Znunu to Photon Ratio             " << endl;           

	// GetScaled histos with only one bin and probagated errors: stat error and error on scale factor
	TH1D *currPhotons= GetScaledHisto(PhotonicSignalRegion_EB->hPhotons , MLRes_EB->fMLPhotonScaleFactor, MLRes_EB->fMLPhotonScaleFactorErr); // EB
	currPhotons->Add(  GetScaledHisto(PhotonicSignalRegion_EE->hPhotons , MLRes_EE->fMLPhotonScaleFactor, MLRes_EE->fMLPhotonScaleFactorErr)); //EE
	TH1D *currZnunu  = GetScaledHisto(HadronicSignalRegion->hZJetsToNuNu, 1                   , 0);

	TH1D *ratio = (TH1D*) currZnunu->Clone("Znunu_To_Photon_ratio");
	ratio->Divide(currPhotons);
	fMCZnunuPhotonRatio    = ratio->GetBinContent(1);
	fMCZnunuPhotonRatioErr = ratio->GetBinError(1);

	cout << "Ratio: " << fMCZnunuPhotonRatio << " pm " << fMCZnunuPhotonRatioErr << endl;
	cout << "------------------------------------------------------------------------" << endl;
	delete currPhotons;
	delete currZnunu;
	delete ratio;
}

void Prediction::MakePrediction(){
	if(fMCZnunuPhotonRatio==0)        {cout << "ERROR in MakePrediction: fMCZnunuPhotonRatio==0"        << endl; exit(-1); }
	cout << "******************************* Prediction ****************************************" << endl;

	cout << "Photonic Region where EB Normalization was extracted: (ML-fit result) ------------" << endl;
	cout << "  NData: "    << MLRes_EB->fMLNumData    << " pm " << MLRes_EB->fMLNumDataErr       << endl;
	cout << "  NPhotons: " << MLRes_EB->fMLNumPhotons << " pm " << MLRes_EB->fMLNumPhotonsErr    << endl;
	cout << "  NOther: "   << MLRes_EB->fMLNumOther   << " pm " << MLRes_EB->fMLNumOtherErr      << endl;
	cout << "  NQCD: "     << MLRes_EB->fMLNumQCD	  << " pm " << MLRes_EB->fMLNumQCDErr        << endl;
	cout << "Photonic Region where EE Normalization was extracted: (ML-fit result) ------------" << endl;
	cout << "  NData: "    << MLRes_EE->fMLNumData    << " pm " << MLRes_EE->fMLNumDataErr       << endl;
	cout << "  NPhotons: " << MLRes_EE->fMLNumPhotons << " pm " << MLRes_EE->fMLNumPhotonsErr    << endl;
	cout << "  NOther: "   << MLRes_EE->fMLNumOther   << " pm " << MLRes_EE->fMLNumOtherErr      << endl;
	cout << "  NQCD: "     << MLRes_EE->fMLNumQCD	  << " pm " << MLRes_EE->fMLNumQCDErr        << endl;
	
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


	cout << "Photonic Signal Region: event yields:-------------------" << endl;
	cout << "EB:                                                     " << endl;
	cout << "  NData:    "    << hData_EB->GetBinContent(1)    << " pm " << hData_EB   ->GetBinError(1)  << endl;
	cout << "  NPhotons: "    << hPhotons_EB->GetBinContent(1) << " pm " << hPhotons_EB->GetBinError(1)  << endl;
	cout << "  NQCD:     "    << hQCD_EB->GetBinContent(1)     << " pm " << hQCD_EB    ->GetBinError(1)  << endl;
	cout << "  NOther:   "    << hOther_EB->GetBinContent(1)   << " pm " << hOther_EB  ->GetBinError(1)  << endl;
	cout << "EE:                                                     " << endl;
	cout << "  NData:    "    << hData_EE->GetBinContent(1)    << " pm " << hData_EE   ->GetBinError(1)  << endl;
	cout << "  NPhotons: "    << hPhotons_EE->GetBinContent(1) << " pm " << hPhotons_EE->GetBinError(1)  << endl;
	cout << "  NQCD:     "    << hQCD_EE->GetBinContent(1)     << " pm " << hQCD_EE    ->GetBinError(1)  << endl;
	cout << "  NOther:   "    << hOther_EE->GetBinContent(1)   << " pm " << hOther_EE  ->GetBinError(1)  << endl;
	cout << "total:                                                  " << endl;
	cout << "  NData:    "    << hData->GetBinContent(1)    << " pm " << hData   ->GetBinError(1)  << endl;
	cout << "  NPhotons: "    << hPhotons->GetBinContent(1) << " pm " << hPhotons->GetBinError(1)  << endl;
	cout << "  NQCD:     "    << hQCD->GetBinContent(1)     << " pm " << hQCD    ->GetBinError(1)  << endl;
	cout << "  NOther:   "    << hOther->GetBinContent(1)   << " pm " << hOther  ->GetBinError(1)  << endl;
	cout << "where the following scale factors were used:            " << endl;
	cout << "EB:                                                     " << endl;
	cout << "  Photon Scale Factor: " << MLRes_EB->fMLPhotonScaleFactor << " pm " << MLRes_EB->fMLPhotonScaleFactorErr            << endl;
	cout << "  QCD    Scale Factor: " << MLRes_EB->fMLQCDScaleFactor    << " pm " << MLRes_EB->fMLQCDScaleFactorErr               << endl;
	cout << "  Other  Scake Factor: " << MLRes_EB->fMLOtherScaleFactor  << " pm " << MLRes_EB->fMLOtherScaleFactorErr             << endl;
	cout << "EE:                                                     " << endl;
	cout << "  Photon Scale Factor: " << MLRes_EE->fMLPhotonScaleFactor << " pm " << MLRes_EE->fMLPhotonScaleFactorErr            << endl;
	cout << "  QCD    Scale Factor: " << MLRes_EE->fMLQCDScaleFactor    << " pm " << MLRes_EE->fMLQCDScaleFactorErr               << endl;
	cout << "  Other  Scake Factor: " << MLRes_EE->fMLOtherScaleFactor  << " pm " << MLRes_EE->fMLOtherScaleFactorErr             << endl;
	
	cout << "MC Znunu To Photon ratio ----------------------------------- " << endl;
	cout << fMCZnunuPhotonRatio  << " pm " << fMCZnunuPhotonRatioErr        << endl;	

	float PredictedZnunu                = fMCZnunuPhotonRatio*(hData->GetBinContent(1)-hOther->GetBinContent(1)-hQCD->GetBinContent(1));
	float PredictedZnunu_ErrSys         = GetZnunuPreditionErrorSys (hData, hOther, hQCD);
	float PredictedZnunu_ErrSysClosure  = GetZnunuPreditionErrorSysClosure (hPhotons);
	float PredictedZnunu_ErrStat        = GetZnunuPreditionErrorStat(hData);
	cout << "Predicted N Znunu events in Hadronic Signal region           "                                           << endl;
	cout << "Prediction: "                                                                                            << endl;
	cout << PredictedZnunu << " pm " << PredictedZnunu_ErrSys  << " pm " << PredictedZnunu_ErrStat << " stat "        << endl;
	cout << "True N Znunu events:"                                                                                    << endl;
	cout << HadronicSignalRegion->hZJetsToNuNu->Integral()                                                            << endl;
	cout << "MC cloure:"                                                                                              << endl;
	cout << fMCZnunuPhotonRatio*(hPhotons->GetBinContent(1))   << " pm " << PredictedZnunu_ErrSysClosure << " (sys) " << endl; 

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

float Prediction::GetZnunuPreditionErrorSys(TH1D* hData, TH1D* hOther, TH1D* hQCD){
	// histograms must be scales with prober uncertainties
	float pred_err_2 = pow(fMCZnunuPhotonRatioErr,2)*pow(hData->GetBinContent(1) - hOther->GetBinContent(1) - hQCD->GetBinContent(1),2);
	pred_err_2      += pow(fMCZnunuPhotonRatio,2)   *pow(hOther->GetBinError(1),2);
	pred_err_2      += pow(fMCZnunuPhotonRatio,2)   *pow(hQCD->GetBinError(1),2);

	return sqrt(pred_err_2);
}

float Prediction::GetZnunuPreditionErrorSysClosure(TH1D* hPhotons){
	// computes systematic uncerainty on Znunu prediction for MC closure test:
	// takes into account MLfit uncertainty on Photon-MC yield, MC statistics, statistical uncerainty on Znunu/photon ratio
	float pred_err_2 = pow(fMCZnunuPhotonRatioErr,2)*pow(hPhotons->GetBinContent(1),2);
	pred_err_2      += pow(fMCZnunuPhotonRatio,2)   *pow(hPhotons->GetBinError(1),2);
	return sqrt(pred_err_2);
}

float Prediction::GetZnunuPreditionErrorStat(TH1D* hData){
	// statistical uncertainty on prediction due to data statistics
	float pred_err   =  fMCZnunuPhotonRatio * hData->GetBinError(1);
	return pred_err;
}

