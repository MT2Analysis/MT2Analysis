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
using namespace std;
using namespace RooFit ;

//run via root -l -b -q cracked_GammaJetsToZnunu.C

//note that there won't much documentation
//all this is a 'cracked' version of run_GammaJetsToZnunu.C
//to predict not Z(nunu) but Z(ll) data
//this was done to do checks on 'ISR reweighting' and debugging purposes
//this is NOT suited for any kind of prediction

// User Input:  ----------------------------------------------
// -----
//TString fSamplesRemovedPhotons     ="samples/samples_1g_GEst_MET.dat";     
//TString fSamplesHadronic           ="samples/samples_2l_GEst_MET_2.dat";    
TString fSamplesRemovedPhotons     ="samples/samples_1g_GEst_HT.dat";     
TString fSamplesHadronic           ="samples/samples_2l_GEst_HT_2.dat";    
//TString fSamplesRemovedPhotons     ="samples/samples_1g_GEst_extremeHT.dat";
//TString fSamplesHadronic           ="samples/samples_2l_GEst_HT.dat";    
// stearing ------------
Bool_t  fSeparateEBEE              = true; 
Bool_t  fDoPhotonSigmaIEtaIEta     = true; 
Bool_t  fDoPhotonSignalRegion      = true; 
Bool_t  fDoHadronicSignalRegion    = true; 
Bool_t  fDoPrediction              = true; 
Bool_t  fPrintBrunoTable           = true; 
Bool_t  fDoVariableMT2bins         = true; 
Bool_t  fDoPileUpWeights           = true; 
Bool_t  fbSFReWeight               = true; 
Bool_t  fMrennaHack                = true; 
Bool_t  fEnforceAbsIso             = false;
Bool_t  fUseConstantZToGammaR      = false;
Bool_t  fUseConstantZToGammaRdynamic=true; 
Float_t fMinMT2forConstR           = 260;  
// uncertainty
Bool_t  fAddFitIntrinsicUncert     = false; 
Float_t fFitIntrinsicUncert        = 0.0;
Bool_t  fAddRMCUncertainty         = true;  
Float_t fRMCUncertainty            = 0.3;   
// options ---------------
TString fOutDir                    = "../cracked/GammaJetsPrediction/20130514_test/";    
Int_t   fVerbose                   = 6; 
Bool_t  fSaveResults               = true;  
Bool_t  fSaveZnunuToGammaRatio     = true;  
Bool_t  fDraw                      = true;  
Bool_t  fWriteToFile               = false; 
// -- specify constant Z/gamma ratio
// pileup weights 4.7 fb-1: new Lumi        
Float_t fConstantZToGammaR_LowHT   = 0.458; 
Float_t fConstantZToGammaErr_LowHT = 0.057; 
Float_t fConstantZToGammaR_HighHT  = 0.628; 
Float_t fConstantZToGammaErr_HighHT= 0.120; 
// MT2b specific --------
Bool_t  fMT2b                      = false; 
Bool_t  fDoBRatioCorrection        = true;  
Float_t fRB_MCtoData               = 1.23;  
Float_t fRB_MCtoDataErr            = 0.27;  
Bool_t  fBTagforMLFit              = true;  
Bool_t  fUseFlatMCBRatio           = true;  
//Float_t fFlatMCBRatio              = 
// -------
Float_t fHTmin                     = 750;   
Float_t fHTmax                     = 1200;
Float_t fMT2min                    = 0;
Float_t fMT2max                    = 1E+8;
Int_t   fNBJets                    = 0;
Int_t   fNJets                     = 35;
Int_t   fSR                        = -1;
// ------
Bool_t  fHT                        = true;
Bool_t  fMET                       = false;
Bool_t  fISRreweight               = true;
// ------
/// high HT (HT>1200)
//    int gNMT2bins_2j0b                  = 6;
//    double  gMT2bins_2j0b[gNMT2bins_2j0b+1]   = {120, 150, 200, 260, 350, 550, 900};
//    int gNMT2bins_2j1b                  = 2;
//    double  gMT2bins_2j1b[gNMT2bins_2j1b+1]   = {100, 180, 350};
//    int gNMT2bins_3j0b                  = 7;
//    double  gMT2bins_3j0b[gNMT2bins_3j0b+1]   = {160, 185, 220, 270, 350, 450, 650, 1000};
//    int gNMT2bins_3j1b                  = 4;
//    double  gMT2bins_3j1b[gNMT2bins_3j1b+1]   = {150, 180, 230, 350, 550};
//    int gNMT2bins_3j2b                  = 2;
//    double  gMT2bins_3j2b[gNMT2bins_3j2b+1]   = {130, 200, 350};
//    int gNMT2bins_6j0b                  = 3;
//    double  gMT2bins_6j0b[gNMT2bins_6j0b+1]   = {160, 200, 300, 500};
//    int gNMT2bins_6j1b                  = 3;
//    double  gMT2bins_6j1b[gNMT2bins_6j1b+1]   = {150, 200, 300, 500};
//    int gNMT2bins_6j2b                  = 2;
//    double  gMT2bins_6j2b[gNMT2bins_6j2b+1]   = {130, 200, 350};
//    int gNMT2bins_3b                  = 1;
//    double  gMT2bins_3b[gNMT2bins_3b+1]   = {125, 300};
/// medium HT (750<HT<1200)
//    int gNMT2bins_2j0b                  = 9;
//    double  gMT2bins_2j0b[gNMT2bins_2j0b+1]   = {125, 150, 180, 220, 270, 325, 425, 580, 780, 1000};
//    int gNMT2bins_2j1b                  = 5;
//    double  gMT2bins_2j1b[gNMT2bins_2j1b+1]   = {100, 135, 170, 260, 450, 700};
//    int gNMT2bins_3j0b                  = 9;
//    double  gMT2bins_3j0b[gNMT2bins_3j0b+1]   = {160, 185, 215, 250, 300, 370, 480, 640, 800, 1000};
//    int gNMT2bins_3j1b                  = 6;
//    double  gMT2bins_3j1b[gNMT2bins_3j1b+1]   = {150, 175, 210, 270, 380, 600, 900};
//    int gNMT2bins_3j2b                  = 5;
//    double  gMT2bins_3j2b[gNMT2bins_3j2b+1]   = {130, 160, 200, 270, 370, 500};
//    int gNMT2bins_6j0b                  = 5;
//    double  gMT2bins_6j0b[gNMT2bins_6j0b+1]   = {160, 200, 250, 325, 425, 600};
//    int gNMT2bins_6j1b                  = 4;
//    double  gMT2bins_6j1b[gNMT2bins_6j1b+1]   = {150, 190, 250, 350, 500};
//    int gNMT2bins_6j2b                  = 4;
//    double  gMT2bins_6j2b[gNMT2bins_6j2b+1]   = {130, 170, 220, 300, 450};
//    int gNMT2bins_3b                  = 3;
//    double  gMT2bins_3b[gNMT2bins_3b+1]   = {125, 175, 275, 450};
/// low HT (450<HT<750)
//     const int gNMT2bins_2j0b_lHT                      = 8;
//     double  gMT2bins_2j0b_lHT[gNMT2bins_2j0b_lHT+1]   = {200, 240, 290, 350, 420, 490, 570, 650, 750};
//     const int gNMT2bins_2j1b_lHT                      = 6;
//     double  gMT2bins_2j1b_lHT[gNMT2bins_2j1b_lHT+1]   = {200, 250, 310, 380, 450, 550, 700};
//     const int gNMT2bins_3j0b_lHT                      = 8;
//     double  gMT2bins_3j0b_lHT[gNMT2bins_3j0b_lHT+1]   = {200, 240, 290, 350, 420, 490, 570, 650, 750};
//     const int gNMT2bins_3j1b_lHT                      = 6;
//     double  gMT2bins_3j1b_lHT[gNMT2bins_3j1b_lHT+1]   = {200, 250, 310, 380, 460, 550, 700};
//     const int gNMT2bins_3j2b_lHT                      = 4;
//     double  gMT2bins_3j2b_lHT[gNMT2bins_3j2b_lHT+1]   = {200, 250, 325, 425, 550};
//     const int gNMT2bins_6j0b_lHT                      = 3;
//     double  gMT2bins_6j0b_lHT[gNMT2bins_6j0b_lHT+1]   = {200, 280, 380, 520};
//     const int gNMT2bins_6j1b_lHT                      = 3;
//     double  gMT2bins_6j1b_lHT[gNMT2bins_6j1b_lHT+1]   = {200, 250, 325, 450};
//     const int gNMT2bins_6j2b_lHT                      = 3;
//     double  gMT2bins_6j2b_lHT[gNMT2bins_6j2b_lHT+1]   = {200, 250, 300, 400};
//     const int gNMT2bins_3b_lHT                        = 2;
//     double  gMT2bins_3b_lHT  [gNMT2bins_3b_lHT+1]     = {200, 280, 400};

const int gNMT2bins = 9;
    double  gMT2bins[gNMT2bins+1]   = {160, 185, 215, 250, 300, 370, 480, 640, 800, 1000};
fMinMT2forConstR           = 270;
// Global Variables ------------------------------------
std::ostringstream  fTriggerStream;
std::ostringstream  fTriggerStreamPhotons;
std::ostringstream  fCutStreamPhotons;
std::ostringstream  fCutStreamPhotonsMT2;
std::ostringstream  fCutStreamSignal;
std::ostringstream* fLogStream     = 0;

// Cut Streams
void DefineCutStreams(float HTmin, float HTmax, float MT2min, float MT2max){
	// Trigger Stream ---------------------------------------------------------------
	//
	if(fHT){
	fTriggerStream << "( ( "
                       << "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || "
                       << "trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
                       << "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1) &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	fTriggerStreamPhotons << "( ("
                              << "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || "
                              << "trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
                              << "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1)&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	}
	if(fMET){
  	fTriggerStream << "( ("
                       << "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
                       << "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	fTriggerStreamPhotons << "(trigger.HLT_SinglePhotons == 1 &&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	}

	// CutStream for SigmaIEtaIEta ------------------------------------------
	fCutStreamPhotons << " " 
	  << "misc.MET>=30"                                                << "&&"
	  << "misc.HT >= " << HTmin << " && misc.HT<= " << HTmax           << "&&"
	  << "NEles==0"                                                    << "&&"
	  << "NMuons==0"                                                   << "&&"
	  << "NTausIDLoose3Hits==0"                                        << "&&"
	  << "misc.Jet0Pass ==1"                                           << "&&"
	  << "misc.Jet1Pass ==1"                                           << "&&"
	  << "misc.PassJet40ID ==1"                                        << "&&"
	  << "misc.Vectorsumpt < 70"                                       << "&&"
	  << "NJetsIDLoose40 >=2"                                          << "&&"
	  << "misc.MinMetJetDPhi4Pt40 >0.3"                                << "&&"
	  << "NPhotons ==1 "                                               << "&&"
	  << "photon[0].isLooseIso"                                        << "&&"
	  << "photon[0].HoverE2012<=0.05"                                  << "&&"
	  << "type1pfmet[0].Pt()<100"                                      << "&&"
	  // Noise
	  << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"             << "&&" // for rare SM samples
	  << "misc.CSCTightHaloIDFlag == 0"                                << "&&"
	  << "misc.trackingFailureFlag==0"                                 << "&&"
	  << "misc.eeBadScFlag==0"                                         << "&&"
	  << "misc.EcalDeadCellTriggerPrimitiveFlag==0"                    << "&&"
	  << "misc.TrackingManyStripClusFlag==0"                           << "&&"
	  << "misc.TrackingTooManyStripClusFlag==0"                        << "&&"
	  << "misc.TrackingLogErrorTooManyClustersFlag==0"                 << "&&"
          << "misc.CrazyHCAL==0"
	  << "&&(type1pfmet[0].Pt()>30||type1pfmet[0].Pt()/misc.CaloMETRaw<=2.)";
	if(fMET){
		fCutStreamPhotons << "&& photon[0].lv.Pt()>=200 ";//&& misc.MET>200&&misc.MT2>=200";//need misc.MET??
	}
	if(fNBJets==0){//improve
		fCutStreamPhotons << "&&NBJets40CSVM==0";
	} else if(fNBJets==1) fCutStreamPhotons << "&&NBJets40CSVM==1";
	else if(fNBJets==2)   fCutStreamPhotons << "&&NBJets40CSVM==2";
	else if(fNBJets==-3)  fCutStreamPhotons << "&&NBJets40CSVM>=3";
	else if(fNBJets==-1)  fCutStreamPhotons << "&&NBJets40CSVM>=1";
	if(!fMrennaHack){
		fCutStreamPhotons << "&&"	
	  	<< "(misc.ProcessID!=6||photon[0].MCmatchexitcode!=1)";
	}else{
		fCutStreamPhotons << "&&"	
		<< "(misc.ProcessID!=6||(photon[0].MCmatchexitcode!=1||photon[0].GenJetMinDR<0.3))"           << "&&"
		<< "(misc.ProcessID!=5||(photon[0].GenJetMinDR>0.3))";
	}
	if(fEnforceAbsIso){
		fCutStreamPhotons << "&&"	
		<< "photon[0].isEGMlooseIso==1";
	}
	if(fNJets==35){
		fCutStreamPhotons << "&&NJetsIDLoose40>=3&&NJetsIDLoose40<=5";
	} else if(fNJets==2){
		fCutStreamPhotons << "&&NJetsIDLoose40==2";
	} else if(fNJets==-6){
		fCutStreamPhotons << "&&NJetsIDLoose40>=6";
	}

	// CutStream for Photon Signal Region ------------------------------------------ 
	fCutStreamPhotonsMT2 << " " 
	  << "misc.MT2>=" << gMT2bins[0]                                   << "&&"
	  << "misc.MET>=30"                                                << "&&"
	  << "misc.HT >=" << HTmin  << " && misc.HT <=" << HTmax           << "&&"
	  << "misc.MT2>=" << MT2min << " && misc.MT2<=" << MT2max          << "&&"
	  << "NEles==0"                                                    << "&&"
	  << "NMuons==0"                                                   << "&&"
	  << "NTausIDLoose3Hits==0"                                        << "&&"
	  << "misc.Jet0Pass ==1"                                           << "&&"
	  << "misc.Jet1Pass ==1"                                           << "&&"
	  << "misc.PassJet40ID ==1"                                        << "&&"
	  << "misc.Vectorsumpt < 70"                                       << "&&"
	  << "NJetsIDLoose40 >=2"                                          << "&&"
	  << "misc.MinMetJetDPhi4Pt40 >0.3"                                << "&&"
	  << "NPhotons ==1 "                                               << "&&"
	  << "photon[0].isLooseID==1"                                      << "&&"//id and iso
	  << "type1pfmet[0].Pt()<100"                                      << "&&"
	  // Noise
	  << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"             << "&&" // for rare SM samples
	  << "misc.CSCTightHaloIDFlag == 0"                                << "&&"
	  << "misc.trackingFailureFlag==0"                                 << "&&"
	  << "misc.eeBadScFlag==0"                                         << "&&"
	  << "misc.EcalDeadCellTriggerPrimitiveFlag==0"                    << "&&"
	  << "misc.TrackingManyStripClusFlag==0"                           << "&&"
	  << "misc.TrackingTooManyStripClusFlag==0"                        << "&&"
	  << "misc.TrackingLogErrorTooManyClustersFlag==0"                 << "&&"
          << "misc.CrazyHCAL==0"
	  << "&&(type1pfmet[0].Pt()>30||type1pfmet[0].Pt()/misc.CaloMETRaw<=2.)";
	if(fMET){
		fCutStreamPhotonsMT2 << "&& photon[0].lv.Pt()>=200 && misc.MET>200";//&&misc.MT2>=200";
	}
	if(fNBJets==0){//improve
		fCutStreamPhotonsMT2 << "&&NBJets40CSVM==0";
	} else if(fNBJets==1) fCutStreamPhotonsMT2 << "&&NBJets40CSVM==1";
	else if(fNBJets==2)   fCutStreamPhotonsMT2 << "&&NBJets40CSVM==2";
	else if(fNBJets==-3)  fCutStreamPhotonsMT2 << "&&NBJets40CSVM>=3";
	else if(fNBJets==-1)  fCutStreamPhotonsMT2 << "&&NBJets40CSVM>=1";
	if(!fMrennaHack){
		fCutStreamPhotonsMT2<< "&&"	
	  	<< "(misc.ProcessID!=6||photon[0].MCmatchexitcode!=1)";
	}else{
		fCutStreamPhotonsMT2<< "&&"	
		<< "(misc.ProcessID!=6||(photon[0].MCmatchexitcode!=1||photon[0].GenJetMinDR<0.3))"           << "&&"
		<< "(misc.ProcessID!=5||(photon[0].GenJetMinDR>0.3))";
	}
	if(fEnforceAbsIso){
		fCutStreamPhotonsMT2<< "&&"	
		<< "photon[0].isEGMlooseIso==1";
	}
	if(fNJets==35){
		fCutStreamPhotonsMT2 << "&&NJetsIDLoose40>=3&&NJetsIDLoose40<=5";
	} else if(fNJets==2){
		fCutStreamPhotonsMT2 << "&&NJetsIDLoose40==2";
	} else if(fNJets==-6){
		fCutStreamPhotonsMT2 << "&&NJetsIDLoose40>=6";
	}

	// CutStream for Hadronic Signal Region ----------------------------------------------
	fCutStreamSignal << " " 
	  << "(NEles+NMuons)==2"                                           << "&&"//it's not hadronic but dileptonic region
	  << "GetDiLeptonPt(0,1,0,10,76,106)>=20"                          << "&&"
	  << "ZllRecalculate()&&"
	  << "misc.MT2>=" << gMT2bins[0]                                   << "&&"
	  << "misc.MET>=30"                                                << "&&"
	  << "misc.HT >=" << HTmin  << " && misc.HT <=" << HTmax           << "&&"
	  << "NTausIDLoose3Hits==0"                                        << "&&"
	  << "misc.Jet0Pass ==1"                                           << "&&"
	  << "misc.Jet1Pass ==1"                                           << "&&"
	  << "misc.PassJet40ID ==1"                                        << "&&"
	  << "NJetsIDLoose40 >=2"                                          << "&&"
	  // Noise
	  << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"             << "&&" // for rare SM samples
	  << "misc.CSCTightHaloIDFlag == 0"                                << "&&"
	  << "misc.trackingFailureFlag==0"                                 << "&&"
	  << "misc.eeBadScFlag==0"                                         << "&&"
	  << "misc.EcalDeadCellTriggerPrimitiveFlag==0"                    << "&&"
	  << "misc.TrackingManyStripClusFlag==0"                           << "&&"
	  << "misc.TrackingTooManyStripClusFlag==0"                        << "&&"
	  << "misc.TrackingLogErrorTooManyClustersFlag==0"                 << "&&"
          << "misc.CrazyHCAL==0"
	  << "&&misc.MET/misc.CaloMETRaw<=2.";
	if(fMET){
		fCutStreamSignal << "&& misc.MET>200";//&&misc.MT2>=200";
	}
	if(fNBJets==0){//improve
		fCutStreamSignal << "&&NBJets40CSVM==0";
	} else if(fNBJets==1) fCutStreamSignal << "&&NBJets40CSVM==1";
	else if(fNBJets==2)   fCutStreamSignal << "&&NBJets40CSVM==2";
	else if(fNBJets==-3)  fCutStreamSignal << "&&NBJets40CSVM>=3";
	else if(fNBJets==-1)  fCutStreamSignal << "&&NBJets40CSVM>=1";
	if(fNJets==35){
		fCutStreamSignal << "&&NJetsIDLoose40>=3&&NJetsIDLoose40<=5";
	} else if(fNJets==2){
		fCutStreamSignal << "&&NJetsIDLoose40==2";
	} else if(fNJets==-6){
		fCutStreamSignal << "&&NJetsIDLoose40>=6";
	}

}

// *********************************** run_GammaJetsToZnunu: specify user input here *************************************************
void cracked_GammaJetsToZnunu(){
	gSystem->Load("libPhysics");
	gSystem->Load("libRooFit") ;
	gSystem->CompileMacro("../MT2Code/src/MT2Shapes.cc", "k");
	
	// logStream
	fLogStream = new std::ostringstream();
	
	// define cutsteams
	DefineCutStreams(fHTmin, fHTmax, fMT2min, fMT2max);

	// fix output dir
	if(fNJets <0)         fOutDir= TString::Format("%s_ge%dj",     fOutDir.Data(), abs(fNJets));
	else                  fOutDir= TString::Format("%s_%dj",       fOutDir.Data(), fNJets);
	if(fNBJets>=0)        fOutDir= TString::Format("%s_%db",       fOutDir.Data(), fNBJets);
	else if(fNBJets!=-10) fOutDir= TString::Format("%s_ge%db",     fOutDir.Data(), abs(fNBJets));
	else                  fOutDir= TString::Format("%s_%d_HT_%s",  fOutDir.Data(), abs(fHTmin),  "Inf");
	if(fHTmax <10000)     fOutDir= TString::Format("%s_%d_HT_%d",  fOutDir.Data(), abs(fHTmin),  abs(fHTmax));
	else                  fOutDir= TString::Format("%s_%d_HT_%s",  fOutDir.Data(), abs(fHTmin),  "Inf");
	if(fMT2max<10000)     fOutDir= TString::Format("%s_%d_MT2_%d", fOutDir.Data(), abs(fMT2min), abs(fMT2max));
	else                  fOutDir= TString::Format("%s_%d_MT2_%s", fOutDir.Data(), abs(fMT2min), "Inf");
	if(fISRreweight)      fOutDir= TString::Format("%s_%s",       "ISRreweighted", fOutDir.Data());
	
	if(!fMT2b || !fDoBRatioCorrection){
		fRB_MCtoData=1; fRB_MCtoDataErr=0;
	}
	if(fNBJets== 0 && fNJets== 2) fSR = 0;
	if(fNBJets==-1 && fNJets== 2) fSR = 1;
	if(fNBJets== 0 && fNJets==35) fSR = 2;
	if(fNBJets== 1 && fNJets==35) fSR = 3;
	if(fNBJets== 2 && fNJets==35) fSR = 4;
	if(fNBJets== 0 && fNJets==-6) fSR = 5;
	if(fNBJets== 1 && fNJets==-6) fSR = 6;
	if(fNBJets== 2 && fNJets==-6) fSR = 7;
	if(fNBJets==-3 && fNJets==-3) fSR = 8;

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
		std::ostringstream cutStreamPhotons_EE;
		if(fSeparateEBEE){
			//idea: add here MT2 bin with for loop
			//then: run loop; get MLRes_EB/EE --> set a normalization histogram
			//i.e. MLRes class --> add histogram for every float --> initialize histogram (here) --> run code with additional MT2 cuts --> take floats of result and fill histogram bin
			cutStreamPhotons_EB << fCutStreamPhotons.str() << "&&abs(photon[0].lv.Eta())<1.4442";//&&misc.MT2>binloweredge&&misc.MT2<=binupperedge
			cutStreamPhotons_EE << fCutStreamPhotons.str() << "&&abs(photon[0].lv.Eta())>1.566";
		} else{
			cutStreamPhotons_EB << fCutStreamPhotons.str() ;
			cutStreamPhotons_EE << fCutStreamPhotons.str() ;
		}
		prediction->PhotonSigmaIEtaIEta_EB = new Channel("SigmaIEtaIEta_EB", "photon[0].SigmaIEtaIEta", 
								 cutStreamPhotons_EB.str().c_str(), fTriggerStreamPhotons.str().c_str(),fSamplesRemovedPhotons);
		prediction->PhotonSigmaIEtaIEta_EB->fVerbose =prediction->fVerbose;
		prediction->PhotonSigmaIEtaIEta_EB->fOutputDir=prediction->fOutputDir;
		prediction->PhotonSigmaIEtaIEta_EB->fAddMCPedestal = true; // set this to avoid
		prediction->PhotonSigmaIEtaIEta_EB->fRootFile="SigmaIEtaIEta_EB_Shapes.root";
		prediction->PhotonSigmaIEtaIEta_EB->GetShapes("SigmaIEtaIEta_EB", "SigmaIEtaIEta", 30, 0, 0.08);
		prediction->GetPhotonNormalization(prediction->PhotonSigmaIEtaIEta_EB, prediction->MLRes_EB);

		// Get Photon Normalization: EE
			prediction->PhotonSigmaIEtaIEta_EE = new Channel("SigmaIEtaIEta_EE", "photon[0].SigmaIEtaIEta", 
									 cutStreamPhotons_EE.str().c_str(), fTriggerStreamPhotons.str().c_str(),fSamplesRemovedPhotons);
			prediction->PhotonSigmaIEtaIEta_EE->fVerbose =prediction->fVerbose;
			prediction->PhotonSigmaIEtaIEta_EE->fOutputDir=prediction->fOutputDir;
			prediction->PhotonSigmaIEtaIEta_EE->fAddMCPedestal = true; // set this to avoid
			prediction->PhotonSigmaIEtaIEta_EE->fRootFile="SigmaIEtaIEta_EE_Shapes.root";
			prediction->PhotonSigmaIEtaIEta_EE->GetShapes("SigmaIEtaIEta_EE", "SigmaIEtaIEta", 30, 0, 0.08);
			prediction->GetPhotonNormalization(prediction->PhotonSigmaIEtaIEta_EE, prediction->MLRes_EE);
			delete prediction->PhotonSigmaIEtaIEta_EE;
		delete prediction->PhotonSigmaIEtaIEta_EB;
	}
	else {
		prediction->SetPhotonNormalizationStupid(prediction->MLRes_EB);
		prediction->SetPhotonNormalizationStupid(prediction->MLRes_EE);
	}

	// Photon Signal Region ******************************************************************************************
	if(fDoPhotonSignalRegion){	
		// Get Photon Selection Signal Region: EB
		std::ostringstream cutStreamPhotonsMT2_EB;
		cutStreamPhotonsMT2_EB << fCutStreamPhotonsMT2.str() << "&&abs(photon[0].lv.Eta())<1.4442";
		prediction->PhotonicSignalRegion_EB = new Channel("PhotonicSignalRegion_EB", "misc.MT2", cutStreamPhotonsMT2_EB.str().c_str(), 
							       fTriggerStreamPhotons.str().c_str(), fSamplesRemovedPhotons);
		prediction->PhotonicSignalRegion_EB->fVerbose =prediction->fVerbose;
		prediction->PhotonicSignalRegion_EB->fOutputDir=prediction->fOutputDir;
		prediction->PhotonicSignalRegion_EB->fRootFile="SignalRegionRemovedPhotons_EB_Shapes.root";
		if(fDoVariableMT2bins) prediction->PhotonicSignalRegion_EB->GetShapes("PhotonicSignalRegion_EB", "MT2 (GeV)", gNMT2bins, gMT2bins );
		else                   prediction->PhotonicSignalRegion_EB->GetShapes("PhotonicSignalRegion_EB", "MT2 (GeV)", 30, 0, 800);
		
		// Get Photon Selection Signal Region: EE
		std::ostringstream cutStreamPhotonsMT2_EE;
		cutStreamPhotonsMT2_EE << fCutStreamPhotonsMT2.str() << "&&abs(photon[0].lv.Eta())>1.566";
		prediction->PhotonicSignalRegion_EE = new Channel("PhotonicSignalRegion_EE", "misc.MT2", cutStreamPhotonsMT2_EE.str().c_str(), 
							       fTriggerStreamPhotons.str().c_str(), fSamplesRemovedPhotons);
		prediction->PhotonicSignalRegion_EE->fVerbose =prediction->fVerbose;
		prediction->PhotonicSignalRegion_EE->fOutputDir=prediction->fOutputDir;
		prediction->PhotonicSignalRegion_EE->fRootFile="SignalRegionRemovedPhotons_EE_Shapes.root";
		if(fDoVariableMT2bins) prediction->PhotonicSignalRegion_EE->GetShapes("PhotonicSignalRegion_EE", "MT2 (GeV)", gNMT2bins, gMT2bins );
		else                   prediction->PhotonicSignalRegion_EE->GetShapes("PhotonicSignalRegion_EE", "MT2 (GeV)", 30, 0, 800);
		
	}

	// Hadronic Signal Region ********************************************************************************************* 
	if(fDoHadronicSignalRegion){
		prediction->HadronicSignalRegion = new Channel("HadronicSignalRegion","misc.MT2", fCutStreamSignal.str().c_str(), 
							       fTriggerStream.str().c_str(), fSamplesHadronic);
		prediction->HadronicSignalRegion->fRootFile="HadronicMT2Shapes.root";
		prediction->HadronicSignalRegion->fVerbose =prediction->fVerbose;
		prediction->HadronicSignalRegion->fOutputDir=prediction->fOutputDir;
		if(fDoVariableMT2bins) prediction->HadronicSignalRegion->GetShapes("HadronicRegion", "MT2 (GeV)", gNMT2bins, gMT2bins );
		else                   prediction->HadronicSignalRegion->GetShapes("HadronicRegion", "MT2 (GeV)", 30, 0, 800);
		// compute MC Znunu/Photon ratio --------------------------------------------------
		prediction->GetMCZnunuToPhotonRatio();
		if(fDoPrediction)  {
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
	delete prediction;
	delete fLogStream;
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
	void GetShapes(TString SelectionName, TString xtitle, const int nbins, float binmin, float binmax);
	void GetShapes(TString SelectionName, TString xtitle, const int nbins, const double *bins );
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

void Channel::GetShapes(TString SelectionName, TString xtitle, const int nbins, const double *bins ){ 
	if(SelectionName!="HadronicRegion"||fHT){//for fHT no statistics?
	MT2Shapes *tA;
	if(fWriteToFile) tA = new MT2Shapes(fOutputDir, fRootFile, fLogStream);
	else             tA = new MT2Shapes(fOutputDir, fRootFile);
	tA->setVerbose(fVerbose);
	tA->init(fSamples);
	tA->SetPrintSummary(true);
	tA->SetDraw(false);
	tA->SetWrite(false);
	tA->SetPileUpWeights(fDoPileUpWeights);
	tA->SetbSFWeights(fbSFReWeight);
  

//                    variable,    cuts,    njet, nbjets,  nlep, selection_name,      HLT,    xtitle   nbins  bins   
        tA->GetShapes(fVariable,  fCuts,    fNJets,  fNBJets, -10  , SelectionName,    fTrigger , xtitle , nbins, bins);

	// retrieve shapes
	for(int i=0; i<tA->GetNShapes(); ++i){
		TString name =tA->fh_shapes[i]->GetName();
		if      (name.Contains("QCD_"))        {hQCD         = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hQCD->SetDirectory(0);}
		else if (name.Contains("PhotonsJets_")){hPhotons     = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hPhotons->SetDirectory(0);}
		else if (name.Contains("Data_"))       {hData        = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hData->SetDirectory(0);}
		else if (name.Contains("WJets_"))      {hWJets       = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hWJets->SetDirectory(0);}
		else if (name.Contains("ZJetsToNuNu_")){cout << "yeah" << endl; hZJetsToNuNu = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName());cout << hZJetsToNuNu->Integral() << endl; hZJetsToNuNu->SetDirectory(0);  }
		else if (name.Contains("ZJetsToLL_"))  {hZJetsToLL   = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hZJetsToLL->SetDirectory(0);}
		else if (name.Contains("Top_"))        {hTop         = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hTop->SetDirectory(0);}
		else if (name.Contains("Signal_"))     {hSignal      = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hSignal->SetDirectory(0);}
		else if (name.Contains("Other_"))      {hOther       = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hOther->SetDirectory(0);}
	}
	delete tA;
} else {
//if use Zll 
	MT2Shapes *tA;
	if(fWriteToFile) tA = new MT2Shapes(fOutputDir, fRootFile, fLogStream);
	else             tA = new MT2Shapes(fOutputDir, fRootFile);
	tA->setVerbose(fVerbose);
	tA->init(fSamples);
	tA->SetPrintSummary(true);
	tA->SetDraw(false);
	tA->SetWrite(false);
	tA->SetPileUpWeights(fDoPileUpWeights);
	tA->SetbSFWeights(fbSFReWeight);
	TString cutsA = fCuts + "&&((NEles==2&&ele[0].lv.Pt()>20&&ele[1].lv.Pt()>20&&ele[0].IDLoose&&ele[1].IDLoose&&NMuons==2)||(NMuons==2&&NEles==0&&muo[0].lv.Pt()>20&&muo[1].lv.Pt()>20))";
	TString triggerA = "trigger.HLT_DiElectrons==1";
      tA->GetShapes(fVariable,  cutsA,    fNJets,  fNBJets, -10  , SelectionName,    triggerA , xtitle , nbins, bins);
	for(int i=0; i<tA->GetNShapes(); ++i){
		TString name =tA->fh_shapes[i]->GetName();
		if (name.Contains("ZJetsToNuNu_")){cout << "yeah" << endl; hZJetsToNuNu = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName());cout << hZJetsToNuNu->Integral() << endl; hZJetsToNuNu->SetDirectory(0);  }
	}
}

	fGotShapes=true;

	if(fISRreweight){//hadronic region here is MC
		cout << "am here1" << endl;
		//reweight only photon jets and Znunu
		TFile *f = TFile::Open("ZGammaISRCorrectionsHistos.root");
		string sig_reg = "";
		string g_reg, z_reg;
		     if(fNJets== 2 && abs(fNBJets)==0) sig_reg = "2j0b";
		else if(fNJets== 2 && abs(fNBJets)==1) sig_reg = "2j1to2b";
		else if(fNJets==35 &&     fNBJets ==0) sig_reg = "3to5j0b";
		else if(fNJets==35 &&     fNBJets ==1) sig_reg = "3to5j1b";
		else if(fNJets==35 &&     fNBJets ==2) sig_reg = "3to5j2b";
		else if(fNJets==-6 &&     fNBJets ==0) sig_reg = "6j0b";
		else if(fNJets==-6 &&     fNBJets ==1) sig_reg = "6j1b";
		else if(fNJets==-6 &&     fNBJets ==2) sig_reg = "6j2b";
		else if(fNJets==-3 &&     fNBJets==-3) sig_reg = "3b";
		//else sig_reg = "3to5j0b";//dummy
		if(fMET){
			g_reg = "MT2_ISRdivNoISR_G_lowHT_" + sig_reg;
			z_reg = "MT2_ISRdivNoISR_Z_lowHT_" + sig_reg;
		}
		if(fHT){
			if(fHTmin==1200){
			g_reg = "MT2_ISRdivNoISR_G_highHT_" + sig_reg;
			z_reg = "MT2_ISRdivNoISR_Z_highHT_" + sig_reg;
			} else{
			g_reg = "MT2_ISRdivNoISR_G_mediumHT_" + sig_reg;
			z_reg = "MT2_ISRdivNoISR_Z_mediumHT_" + sig_reg;
			}
		}
		TH1D *hZ = (TH1D*)f->Get(z_reg.c_str());
		TH1D *hG = (TH1D*)f->Get(g_reg.c_str());
		for(int nbin = 1; nbin<=hZJetsToNuNu->GetNbinsX(); ++nbin){//do it like that in order to not change the error
			if(hPhotons!=0)     hPhotons->SetBinContent(nbin, hPhotons->GetBinContent(nbin)*hG->GetBinContent(nbin));
			if(hPhotons!=0)     hPhotons->SetBinError(  nbin, hPhotons->GetBinError(  nbin)*hG->GetBinContent(nbin));
			if(hZJetsToNuNu!=0) hZJetsToNuNu->SetBinContent(nbin, hZJetsToNuNu->GetBinContent(nbin)*hZ->GetBinContent(nbin));
			if(hZJetsToNuNu!=0) hZJetsToNuNu->SetBinError(  nbin, hZJetsToNuNu->GetBinError(  nbin)*hZ->GetBinContent(nbin));
		}
		f->Close();
	}

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
		delete file;
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
	delete tA;
}
void Channel::GetShapes(TString SelectionName, TString xtitle, const int nbins, float binmin, float binmax){ 
	double bins[nbins];
	bins[0] = binmin;
	for(int i=1; i<=nbins; i++) bins[i] = binmin+i*(binmax-binmin)/nbins;
	GetShapes(SelectionName, xtitle, nbins, bins);
}

class Prediction {
public:
	Prediction();
	~Prediction();
	
	void GetPhotonNormalization(Channel* channel, MLResult* MLRes);
	void SetPhotonNormalizationStupid(MLResult* MLRes);
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
	TH1D* fMCZnunuPhotonRatioHisto;

	
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

	if(channel->fAddMCPedestal){ // Add Pedestal to QCD MC PDF in order to avoid bins in PDF with zero entries
		                     // but data in same bin! this causes problems!
		for(int i=0; i<hQCD->GetNbinsX(); ++i){
			if(hData->GetBinContent(i)==0) continue;
			double MCcontent = hPhotons->GetBinContent(i)+hQCD->GetBinContent(i)+hOther->GetBinContent(i);
			if(MCcontent==0) {hQCD->SetBinContent(i, 1E-01);hQCD->SetBinError(i, 1E-01);}
		}
	}
	
	RooRealVar sigmaietaieta("sigmaietaieta","sigmaietaieta",0.,0.08) ; // contained in histos

	RooDataHist Data   ("data"   ,"data"  ,sigmaietaieta,hData) ;    // define RooDataHists
	RooDataHist Photons("photons","photon",sigmaietaieta,hPhotons);
	RooDataHist QCD    ("QCD"    ,"QCD"   ,sigmaietaieta,hQCD);
	RooDataHist Other  ("Other"    ,"Other"   ,sigmaietaieta,hOther);//check this!!!

	RooHistPdf Photons_pdf("photons_pdf","photons_pdf",sigmaietaieta,Photons); // define PDFs for signal and bkg
	RooHistPdf QCD_pdf    ("qcd_pdf"    ,"qcd_pdf"    ,sigmaietaieta,QCD    ); 
	RooHistPdf Other_pdf  ("Other_pdf"  ,"other_pdf"  ,sigmaietaieta,Other  );

	RooRealVar nsig       ("nsig"   ,"number of signal events",     hPhotons->Integral()  ,  hPhotons->Integral()*0.1,hData->Integral());
	RooRealVar nqcd       ("nqcd"   ,"number of QCD events",        hQCD->Integral()      ,  hQCD->Integral()    *0.5,hData->Integral());
	RooRealVar nother     ("nother" ,"number of Other SM events",   hOther->Integral()); nother.setConstant(kTRUE);

	// model(x) = nsig*Photons_pdf(x) + nqcd*QCD_pdf(x) + nother*Other_pdf(x), where nother is fixed to nominal contribution
	RooAddPdf model("model","model", RooArgList(Photons_pdf,QCD_pdf,Other_pdf), RooArgList(nsig, nqcd, nother));
	model.defaultPrintStream(fLogStream);
	
	// perform fit
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
		delete file;
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

void Prediction::SetPhotonNormalizationStupid(MLResult* MLRes){

	MLRes->fMLPhotonScaleFactor    = 1;
	MLRes->fMLPhotonScaleFactorErr = 0.05;
	MLRes->fMLQCDScaleFactor       = 1;
	MLRes->fMLQCDScaleFactorErr    = 0.1;
	MLRes->fMLOtherScaleFactor     = 1;
	MLRes->fMLOtherScaleFactorErr  = 0;

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
	if(fUseConstantZToGammaR){
		Float_t fConstantZToGammaR   = 0;
		Float_t fConstantZToGammaErr = 0;
		if(fHTmin == 750 && fHTmax == 950){
			fConstantZToGammaR   = fConstantZToGammaR_LowHT;
			fConstantZToGammaErr = fConstantZToGammaErr_LowHT;
		}else if(fHTmin == 950 && fHTmax > 10000){
			fConstantZToGammaR   = fConstantZToGammaR_HighHT;
			fConstantZToGammaErr = fConstantZToGammaErr_HighHT;
		}else{
			cout << "GetMCZnunuToPhotonRatio: ERROR: cannot use flat ratio for this HT binning" << endl;
			exit(-1);
		}
		*fLogStream << ">>> Using constant Z(nunu) to Gamma Ratio= " << fConstantZToGammaR << " pm " << fConstantZToGammaErr << endl; 
		if(fMinMT2forConstR>fMT2min) {
			cout << "GetMCZnunuToPhotonRatio:ERROR: cannot use constant Znunu to Photon ratio for MT2 > " << fMT2min << endl;
			exit(-1);	
		}
		fMCZnunuPhotonRatio    = fConstantZToGammaR;
		if(fAddRMCUncertainty){
		*fLogStream << ">>> adding an additional uncertainty of " << fRMCUncertainty << "\% on (nunu) to Gamma Ratio" << endl;
		fMCZnunuPhotonRatioErr = sqrt(fConstantZToGammaErr*fConstantZToGammaErr + fRMCUncertainty*fRMCUncertainty*fConstantZToGammaR*fConstantZToGammaR);
		}else{
		fMCZnunuPhotonRatioErr = fConstantZToGammaErr;
		}
		*fLogStream << "Ratio: " << fMCZnunuPhotonRatio << " pm " << fMCZnunuPhotonRatioErr << endl;
		*fLogStream << "------------------------------------------------------------------------" << endl;
	}else{
		if(PhotonicSignalRegion_EB->hPhotons==0 || PhotonicSignalRegion_EE->hPhotons==0 || HadronicSignalRegion->hZJetsToNuNu==0){
			cout << "GetMCZnunuToPhotonRatio: received 0 pointer!" << endl;
			exit(-1);
		}

		// GetScaled histos with only one bin and probagated errors: stat error and error on scale factor
		TH1D *currPhotons= GetScaledHisto2(PhotonicSignalRegion_EB->hPhotons , MLRes_EB->fMLPhotonScaleFactor, MLRes_EB->fMLPhotonScaleFactorErr); // EB
		currPhotons->Add(  GetScaledHisto2(PhotonicSignalRegion_EE->hPhotons , MLRes_EE->fMLPhotonScaleFactor, MLRes_EE->fMLPhotonScaleFactorErr)); //EE
		TH1D *currZnunu  = GetScaledHisto2(HadronicSignalRegion->hZJetsToNuNu, 1                   , 0);
		if(fUseConstantZToGammaRdynamic) currPhotons = RefillRatioHisto(currPhotons,fMinMT2forConstR);
		if(fUseConstantZToGammaRdynamic) currZnunu   = RefillRatioHisto(currZnunu  ,fMinMT2forConstR);

		TH1D *ratio = (TH1D*) currZnunu->Clone("Znunu_To_Photon_ratio");
		ratio->Divide(currPhotons);
		fMCZnunuPhotonRatio    = ratio->GetBinContent(1);
		fMCZnunuPhotonRatioErr = ratio->GetBinError(1);
		fMCZnunuPhotonRatioHisto = (TH1D*)ratio->Clone("Znunu_To_Photon_ratio_histo");


		if(fAddRMCUncertainty){
		*fLogStream << "------------------------------------------------------------------------------" << endl;
		*fLogStream << "+++ adding in quadrature " << fRMCUncertainty << " percent uncertainty on R   " << endl;
		*fLogStream << "------------------------------------------------------------------------------" << endl;
		fMCZnunuPhotonRatioErr = sqrt(pow(fMCZnunuPhotonRatioErr,2)+pow(fMCZnunuPhotonRatio*fRMCUncertainty,2));
		for(int i = 1; i<= fMCZnunuPhotonRatioHisto->GetNbinsX(); ++i) fMCZnunuPhotonRatioHisto->SetBinError(i, sqrt(pow(fMCZnunuPhotonRatioHisto->GetBinError(i),2)+pow(fMCZnunuPhotonRatioHisto->GetBinContent(i)*fRMCUncertainty,2)));
		//TOBEUPDATED: running additional uncertainty
		}

		*fLogStream << "Ratio: " << fMCZnunuPhotonRatio << " pm " << fMCZnunuPhotonRatioErr << endl;
		*fLogStream << "------------------------------------------------------------------------" << endl;

		if(fSaveZnunuToGammaRatio){
			TH1D *currPhotons2= GetScaledHisto(PhotonicSignalRegion_EB->hPhotons , MLRes_EB->fMLPhotonScaleFactor, MLRes_EB->fMLPhotonScaleFactorErr, 1); // EB
			currPhotons2->Add(  GetScaledHisto(PhotonicSignalRegion_EE->hPhotons , MLRes_EE->fMLPhotonScaleFactor, MLRes_EE->fMLPhotonScaleFactorErr, 1)); //EE
			TH1D *currZnunu2  = GetScaledHisto(HadronicSignalRegion->hZJetsToNuNu, 1                               , 0                              , 1);
			TH1D *ratio2 = (TH1D*) currZnunu2->Clone("Znunu_To_Photon_Ratio");
			ratio2->Divide(currPhotons2);
			TString filename=fOutputDir+"/ZnunuToGammaRatio.root";
			TFile *file = new TFile(filename.Data(), "RECREATE");
			currZnunu2->Write();
			currPhotons2->Write();
			ratio2->Write();
			fMCZnunuPhotonRatioHisto->Write();
			file->Close();
			delete currPhotons2;
			delete currZnunu2;
			delete ratio2;
		}
		delete currPhotons;
		delete currZnunu;
		delete ratio;
	}
}

void Prediction::MakePrediction(){
	if(fMCZnunuPhotonRatio==0)        {*fLogStream << "ERROR in MakePrediction: fMCZnunuPhotonRatio==0"        << endl; exit(-1); }
	*fLogStream << "******************************* Prediction ****************************************" << endl;

	*fLogStream << "Photonic Region where EB Normalization was extracted: (ML-fit result) ------------" << endl;
	*fLogStream << "  NData: "    << MLRes_EB->fMLNumData    << " pm " << MLRes_EB->fMLNumDataErr       << endl;
	*fLogStream << "  NPhotons: " << MLRes_EB->fMLNumPhotons << " pm " << MLRes_EB->fMLNumPhotonsErr    << endl;
	*fLogStream << "  NOther: "   << MLRes_EB->fMLNumOther   << " pm " << MLRes_EB->fMLNumOtherErr      << endl;
	*fLogStream << "  NQCD: "     << MLRes_EB->fMLNumQCD	 << " pm " << MLRes_EB->fMLNumQCDErr        << endl;
	*fLogStream << "Photonic Region where EE Normalization was extracted: (ML-fit result) ------------" << endl;
	*fLogStream << "  NData: "    << MLRes_EE->fMLNumData    << " pm " << MLRes_EE->fMLNumDataErr       << endl;
	*fLogStream << "  NPhotons: " << MLRes_EE->fMLNumPhotons << " pm " << MLRes_EE->fMLNumPhotonsErr    << endl;
	*fLogStream << "  NOther: "   << MLRes_EE->fMLNumOther   << " pm " << MLRes_EE->fMLNumOtherErr      << endl;
	*fLogStream << "  NQCD: "     << MLRes_EE->fMLNumQCD	 << " pm " << MLRes_EE->fMLNumQCDErr        << endl;
	
	// EB
	TH1D *hPhotons_EB = GetScaledHisto2(PhotonicSignalRegion_EB->hPhotons, MLRes_EB->fMLPhotonScaleFactor, MLRes_EB->fMLPhotonScaleFactorErr);
	TH1D *hOther_EB   = GetScaledHisto2(PhotonicSignalRegion_EB->hOther  , MLRes_EB->fMLOtherScaleFactor , MLRes_EB->fMLOtherScaleFactorErr);
	TH1D *hQCD_EB     = GetScaledHisto2(PhotonicSignalRegion_EB->hQCD    , MLRes_EB->fMLQCDScaleFactor, MLRes_EB->fMLQCDScaleFactorErr);
	TH1D *hData_EB    = GetScaledHisto2(PhotonicSignalRegion_EB->hData   , 1                   , 0);
	// EE
	TH1D *hPhotons_EE = GetScaledHisto2(PhotonicSignalRegion_EE->hPhotons, MLRes_EE->fMLPhotonScaleFactor, MLRes_EE->fMLPhotonScaleFactorErr);
	TH1D *hOther_EE   = GetScaledHisto2(PhotonicSignalRegion_EE->hOther  , MLRes_EE->fMLOtherScaleFactor , MLRes_EE->fMLOtherScaleFactorErr);
	TH1D *hQCD_EE     = GetScaledHisto2(PhotonicSignalRegion_EE->hQCD    , MLRes_EE->fMLQCDScaleFactor, MLRes_EE->fMLQCDScaleFactorErr);
	TH1D *hData_EE    = GetScaledHisto2(PhotonicSignalRegion_EE->hData   , 1                   , 0);
	
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
	if(fUseConstantZToGammaR){
	*fLogStream << " >> ---------- using constant ratio ---------------- <<      " << endl;
	}
	*fLogStream << fMCZnunuPhotonRatio  << " pm " << fMCZnunuPhotonRatioErr        << endl;	
	
	if(fMT2b && fDoBRatioCorrection){
	*fLogStream << "MC to Data Photon+jets B-tag correction: " << fRB_MCtoData << " pm " << fRB_MCtoDataErr << endl;
	}

	float PredictedZnunu                = fRB_MCtoData*fMCZnunuPhotonRatio*(hData->GetBinContent(1)-hOther->GetBinContent(1)-hQCD->GetBinContent(1));


	float PredictedZnunu_ErrSys         = GetZnunuPreditionErrorSys (hData, hOther, hQCD);
	float PredictedZnunu_ErrSysClosure  = GetZnunuPreditionErrorSysClosure (hPhotons);
	float PredictedZnunu_ErrStat        = GetZnunuPreditionErrorStat(hData);

	vector<int> htbin, mt2bin;
	vector<double> znunugen, znunupred, znunuprederr, znunuprederrstat, znunuprederrsyst, ZGratio, ZGratioerr, mt2low, mt2up, ndata, ndataerr, nqcd, nqcderr, nother, nothererr;
	htbin.clear(), mt2bin.clear(), znunugen.clear(), znunupred.clear(), znunuprederr.clear(), znunuprederrstat.clear(), znunuprederrsyst.clear(), ZGratio.clear(), ZGratioerr.clear(), mt2low.clear(), mt2up.clear(), ndata.clear(), ndataerr.clear(), nqcd.clear(), nqcderr.clear(), nother.clear(), nothererr.clear();
	for(int i = 1; i<=hData->GetNbinsX(); ++i){

	PredictedZnunu                = fRB_MCtoData*fMCZnunuPhotonRatioHisto->GetBinContent(i)*(hData->GetBinContent(i)-hOther->GetBinContent(i)-hQCD->GetBinContent(i));
	PredictedZnunu_ErrSys         = GetZnunuPreditionErrorSys (i, hData, hOther, hQCD);
	PredictedZnunu_ErrSysClosure  = GetZnunuPreditionErrorSysClosure (i, hPhotons);
	PredictedZnunu_ErrStat        = GetZnunuPreditionErrorStat(i, hData);

	*fLogStream << "For " << hPhotons->GetBinLowEdge(i) << " < MT2 < " << hPhotons->GetBinLowEdge(i)+hPhotons->GetBinWidth(i) <<  "(" << i << "/" << hData->GetNbinsX() << ")" << endl;

	*fLogStream << "Photon events in Photon Signal region (data-bg): " <<  hData->GetBinContent(i) << " - " << hOther->GetBinContent(i)+hQCD->GetBinContent(i) << " MC: " << hPhotons->GetBinContent(i) << endl;
	*fLogStream << "Predicted N Znunu events in Hadronic Signal region: ratio " << fMCZnunuPhotonRatioHisto->GetBinContent(i) << " pm " <<fMCZnunuPhotonRatioHisto->GetBinError(i) << endl;
	*fLogStream << "Prediction: "                                                                                            << endl;
	*fLogStream << "  " << PredictedZnunu << " pm " << PredictedZnunu_ErrSys  << " pm " << PredictedZnunu_ErrStat << " stat "        << endl;
	*fLogStream << "True N Znunu events:"                                                                                    << endl;
	*fLogStream << "  " << HadronicSignalRegion->hZJetsToNuNu->GetBinContent(i)                                                           << endl;
	*fLogStream << "MC closure:"                                                                                              << endl;
	*fLogStream << "  " << fMCZnunuPhotonRatioHisto->GetBinContent(i)*(hPhotons->GetBinContent(i))   << " pm " << PredictedZnunu_ErrSysClosure << " (sys) " << endl; 

	if(fPrintBrunoTable){//note that this is differently organized than in run_GammaJetsToZnunu.C - for prediction we need fPrintBrunoTable = true
		int MT2bin;
		int HTbin;
		if      (fMT2min == 150 && fMT2max == 200) MT2bin=0;
		else if (fMT2min == 200 && fMT2max == 275) MT2bin=1;
		else if (fMT2min == 275 && fMT2max == 375) MT2bin=2;
		else if (fMT2min == 375 && fMT2max == 500) MT2bin=3;
		else if (fMT2min == 500)                   MT2bin=4;
		else if (fMT2min ==   0 && fMT2max ==1e+8) MT2bin=-1;
		else    {cout << "MT2bin not valid for printout! " << endl; exit(-1);}
		if      (fHTmin  == 750 && fHTmax == 1200) HTbin =1;
		else if (fHTmin  == 1200)                  HTbin =2;
		else if (fHTmin  == 450 && fHTmax == 750 ) HTbin =0;
		else if (fHTmin  == 750 )                  HTbin =10;
		else    {cout << "HTbin not valid for printout! " << endl; exit(-1);}
		*fLogStream << "************************** Bruno-Gay Printout ******************************" << endl;
	        (*fLogStream).precision(3) ;
		*fLogStream << "ZinvFromG " << HTbin << " " << MT2bin << " " << HadronicSignalRegion->hZJetsToNuNu->GetBinContent(i)
		            << " " << PredictedZnunu << " " << sqrt(PredictedZnunu_ErrSys*PredictedZnunu_ErrSys+PredictedZnunu_ErrStat*PredictedZnunu_ErrStat)
		    	    << " " << fMCZnunuPhotonRatioHisto->GetBinContent(i) << " " << fMCZnunuPhotonRatioHisto->GetBinError(i) << endl;	    
		*fLogStream << "----------------------------------------------------------------------------" << endl;
		*fLogStream << "$" << hPhotons->GetBinLowEdge(i) << "-" << hPhotons->GetBinLowEdge(i)+hPhotons->GetBinWidth(i) << "$ &" << hData->GetBinContent(i)  << " $\\pm$ " << hData->GetBinError(i) << " & " 
		            << hQCD->GetBinContent(i) << " $\\pm$ " << hQCD->GetBinError(i)  << " & "
		            << hOther->GetBinContent(i) << " $\\pm$ " << hOther->GetBinError(i)  << " & "
			    << fMCZnunuPhotonRatioHisto->GetBinContent(i)      << " $\\pm$ " << fMCZnunuPhotonRatioHisto->GetBinError(i)  << " & "
			    << PredictedZnunu           << " $\\pm$ " << PredictedZnunu_ErrSys   << " $\\pm$ " << PredictedZnunu_ErrStat   << " & "
			    << HadronicSignalRegion->hZJetsToNuNu->GetBinContent(i)  << " \\\\ " << endl;
		*fLogStream << "************************** Bruno-Gay Printout ******************************" << endl;
		htbin.push_back(HTbin);
		mt2bin.push_back(MT2bin);
		znunugen.push_back(HadronicSignalRegion->hZJetsToNuNu->GetBinContent(i));
		znunupred.push_back(PredictedZnunu);
		znunuprederr.push_back(sqrt(PredictedZnunu_ErrSys*PredictedZnunu_ErrSys+PredictedZnunu_ErrStat*PredictedZnunu_ErrStat));
		znunuprederrstat.push_back(PredictedZnunu_ErrStat);
		znunuprederrsyst.push_back(PredictedZnunu_ErrSys);
		ZGratio.push_back(fMCZnunuPhotonRatioHisto->GetBinContent(i));
		ZGratioerr.push_back(fMCZnunuPhotonRatioHisto->GetBinError(i));
		mt2low.push_back(hPhotons->GetBinLowEdge(i));
		mt2up.push_back(hPhotons->GetBinLowEdge(i)+hPhotons->GetBinWidth(i));
		ndata.push_back(hData->GetBinContent(i));
		ndataerr.push_back(hData->GetBinError(i));
		nqcd.push_back(hQCD->GetBinContent(i));
		nqcderr.push_back(hQCD->GetBinError(i));
		nother.push_back(hOther->GetBinContent(i));
		nothererr.push_back(hOther->GetBinError(i));
	}
	}
	*fLogStream << endl << "************************** Bruno-Gay Printout ******************************" << endl;
        *fLogStream << fixed << setprecision(2) ;
	*fLogStream << " Name     " << " SR " << " HTbin " << " MT2bin   " << " MCpred   " << " DataDrivenPred " << " DataDrivenPredError " << " ScaleFactor " << " ScaleFactorError " << endl;
	for(int jj = 0; jj<htbin.size(); ++jj){
		*fLogStream << "ZinvFromG   " << fSR << " " << setw(5) << htbin[jj] << " " << setw(4) << mt2bin[jj] << " " << setw(13) << znunugen[jj]
		            << " " << setw(13) << znunupred[jj] << " " << setw(18) << znunuprederr[jj]
		    	    << " " << setw(13) << ZGratio[jj]   << " " << setw(16) << ZGratioerr[jj] << endl;
	}
		*fLogStream << "----------------------------------------------------------------------------" << endl;
	*fLogStream << "\%BEGINLATEX\%"                    << endl
		    << "\\begin{table}[!htb]"              << endl
		    << "\\begin{center}"                   << endl
		    << "\\begin{tabular}{l|cccccc}"        << endl
		    << "\\hline\\hline"                    << endl
		    << "\\multirow{2}{*}{$M_{T2}$\\ bin}";
		    if(fHTmin==750 && fHTmax==950) *fLogStream << " & \\multicolumn{6}{c}{$750 < H_T<950$ GeV}  \\";
	*fLogStream << "  & $N^{\\gamma}$         & $QCD^{bkg}$&           $EWK^{bkg}$      &    $R^{MC}(Z(\\nu\\nu)/\\gamma)$ & data pred                   & MC estimate \\\\" << endl
		    << "\\hline" << endl;
	for(int jj = 0; jj<htbin.size(); ++jj){
		*fLogStream << "$" << int(mt2low[jj]) << "-" << int(mt2up[jj]) << "$" 
		            << " " << setw(19) << "& $" << ndata[jj]      << " \\pm "<< ndataerr[jj]   << "$" << " " << setw(1)  << "& " 
		            << " " << setw(2)  <<   "$" << nqcd[jj]       << " \\pm "<< nqcderr[jj]    << "$" << " " << setw(1)  << "& "
		            << " " << setw(2)  <<   "$" << nother[jj]     << " \\pm "<< nothererr[jj]  << "$" << " " << setw(1)  << "& "
		            << " " << setw(4)  <<   "$" << ZGratio[jj]    << " \\pm "<< ZGratioerr[jj] << "$" << " " << setw(1)  << "& "
		            << " " << setw(7)  <<   "$" << znunupred[jj]  << " \\pm " << znunuprederrstat[jj] << " \\pm " << znunuprederrsyst[jj] << "$ & "
			    << " " << setw(4)  << znunugen[jj] << " \\\\ " << endl;
	}
	*fLogStream << "\\hline\\hline"  << endl
		    << "\\end{tabular}"  << endl
		    << "\\end{center}"   << endl
		    << "\\end{table}"    << endl
		    << "\%ENDLATEX\%"    << endl
		    << endl;
	*fLogStream << "************************** Bruno-Gay Printout ******************************" << endl;

}

TH1D* Prediction::RebinHisto(TH1D* histo, float lastbinlowedge){
	TH1D *h        = (TH1D*) histo->Clone(histo->GetName());
	if(h->GetNbinsX()<2) return h;
	if(h->GetBinLowEdge(h->GetNbinsX())<lastbinlowedge) return h;
	vector<double> bintemp; bintemp.clear();
	bintemp.push_back(h->GetBinLowEdge(1));//this has to be there
	for(int i = 2; i<=h->GetNbinsX(); ++i){
		//note: int should make sure that no rounding error screws this up
		if((int)h->GetBinLowEdge(i)<=(int)lastbinlowedge) bintemp.push_back(h->GetBinLowEdge(i));
	}
	bintemp.push_back(h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(h->GetNbinsX()));
	const int arraysize = (bintemp.size()-1);
	const int arraysize2 = bintemp.size();
	double array[arraysize2];
	for(int i = 0; i<=arraysize; ++i) array[i] = bintemp[i];
	string newname = h->GetName();
	newname += "_withconstratio";
	TH1D *hnew = (TH1D*)h->Rebin(arraysize, newname.c_str(), array);
	return hnew;
}

TH1D* Prediction::RefillRatioHisto(TH1D* histo, float lastbinlowedge){
	TH1D *h        = (TH1D*) histo->Clone(histo->GetName());
	if(h->GetNbinsX()<2) return h;
	if(h->GetBinLowEdge(h->GetNbinsX())<lastbinlowedge) return h;
	double content(0.), error2(0.);
	for(int i = 1; i<=h->GetNbinsX(); ++i){
		//note: int should make sure that no rounding error screws this up
		if((int)h->GetBinLowEdge(i)>=(int)lastbinlowedge){
			content += h->GetBinContent(i);
			error2  += pow(h->GetBinError(i),2);
		}
	}
	for(int i = 1; i<=h->GetNbinsX(); ++i){
		//note: int should make sure that no rounding error screws this up
		if((int)h->GetBinLowEdge(i)>=(int)lastbinlowedge){
			h->SetBinContent(i, content);
			h->SetBinError(i, sqrt(error2));
		}
	}
	return h;
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

TH1D* Prediction::GetScaledHisto2(TH1D* histo, float scale_fact, float scale_fact_err){
	// takes histo, scale factor and uncertainty on scale factor
	// and returns scaled and rebinned histo with 1 bin with uncertainty on scale factor propagated.
	TH1D *h        = (TH1D*) histo->Clone(histo->GetName());
	for(int i = 1; i<=h->GetNbinsX(); ++i){
	h->SetBinError(i, sqrt(h->GetBinError(i)*  h->GetBinError(i)   *scale_fact    *scale_fact + 
			  h->GetBinContent(i)*h->GetBinContent(i) *scale_fact_err*scale_fact_err));
	h->SetBinContent(i, h->GetBinContent(i)*scale_fact);
	}
	return h;
}

TH1D* Prediction::GetScaledHisto(TH1D* histo, float scale_fact, float scale_fact_err, int ngroup){
	// takes histo, scale factor and uncertainty on scale factor
	// and returns scaled and rebinned histo with ngroup bins merged into 1 with uncertainty on scale factor propagated.
	if(ngroup>=histo->GetNbinsX()) ngroup=histo->GetNbinsX();
	TH1D *h        = (TH1D*) histo->Clone(histo->GetName());
	h->Rebin(ngroup);
	for(int i=1; i<=h->GetNbinsX(); ++i){
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

float Prediction::GetZnunuPreditionErrorSys(int i, TH1D* hData, TH1D* hOther, TH1D* hQCD){
	// histograms must be scales with prober uncertainties
	if(i>hData->GetNbinsX() || i<=0) return -999.;
	float pred_err_2 = pow(fMCZnunuPhotonRatioHisto->GetBinError(i)  * (hData->GetBinContent(i) - hOther->GetBinContent(i) - hQCD->GetBinContent(i))*fRB_MCtoData,2);
	pred_err_2      += pow(fMCZnunuPhotonRatioHisto->GetBinContent(i)*fRB_MCtoData * hOther->GetBinError(i),2);
	pred_err_2      += pow(fMCZnunuPhotonRatioHisto->GetBinContent(i)*fRB_MCtoData * hQCD->GetBinError(i),2);
	pred_err_2      += pow(fRB_MCtoDataErr* fMCZnunuPhotonRatioHisto->GetBinContent(i)*(hData->GetBinContent(i) - hOther->GetBinContent(i) - hQCD->GetBinContent(i)),2);

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

float Prediction::GetZnunuPreditionErrorSysClosure(int i, TH1D* hPhotons){
	// computes systematic uncerainty on Znunu prediction for MC closure test:
	// takes into account MLfit uncertainty on Photon-MC yield, MC statistics, statistical uncerainty on Znunu/photon ratio
	if(i>hPhotons->GetNbinsX() || i<=0) return -999.;
	float pred_err_2 = pow(fMCZnunuPhotonRatioHisto->GetBinError(i)     * hPhotons->GetBinContent(i)*fRB_MCtoData,2);
	pred_err_2      += pow(fMCZnunuPhotonRatioHisto->GetBinContent(i)   * hPhotons->GetBinError(i)  *fRB_MCtoData,2);
	pred_err_2      += pow(fMCZnunuPhotonRatioHisto->GetBinContent(i)   * hPhotons->GetBinContent(i)*fRB_MCtoDataErr,2);
	return sqrt(pred_err_2);
}

float Prediction::GetZnunuPreditionErrorStat(TH1D* hData){
	// statistical uncertainty on prediction due to data statistics
	float pred_err   = fRB_MCtoData* fMCZnunuPhotonRatio * hData->GetBinError(1);
	return pred_err;
}

float Prediction::GetZnunuPreditionErrorStat(int i, TH1D* hData){
	// statistical uncertainty on prediction due to data statistics
	if(i>hData->GetNbinsX() || i<=0) return -999.;
	float pred_err   = fRB_MCtoData* fMCZnunuPhotonRatioHisto->GetBinContent(i) * hData->GetBinError(i);
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
		delete file;
	}
}
