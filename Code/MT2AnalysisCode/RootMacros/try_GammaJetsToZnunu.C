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

//run via root -l -b -q try_GammaJetsToZnunu.C

//this is a copy of run_GammaJetsToZnunu.C where I try to be more efficient running over all topological regions.
//use at own risk - easiest is checking with the results of run_GammaJetsToZnunu.C


// User Input:  ----------------------------------------------
// -----
//TString fSamplesRemovedPhotons     ="samples/samples_1g_GEst_MET_filter.dat";        // samples with photons added to MET - low HT region
//TString fSamplesHadronic           ="samples/samples_had_GEst_MET_filter.dat";       // hadronic selection                - low HT region
TString fSamplesRemovedPhotons     ="samples/samples_1g_GEst_HT_filter.dat";           // samples with photons added to MET - medium HT region
TString fSamplesHadronic           ="samples/samples_had_GEst_HT_filter.dat";          // hadronic selection                - medium HT region
//TString fSamplesRemovedPhotons     ="samples/samples_1g_GEst_extremeHT_filter.dat";  // samples with photons added to MET - high HT region (note that all sample types QCD/Gamma/Other need some entries for the fitting)
//TString fSamplesHadronic           ="samples/samples_had_GEst_extremeHT_filter.dat"; // hadronic selection                - high HT selection
// stearing ------------
Bool_t  fSeparateEBEE              = true;  // set to false if EE and EB should not be separated in purity measurement of sigmaietaieta fit, default = true
Bool_t  fDoPhotonSigmaIEtaIEta     = true;  // set to "true" to perfom sigmaIetaIeta fit to obtain QCD normalization, default = true
Bool_t  fDoPhotonSignalRegion      = true; // get shapes for photons in signal region <-- needed for background prediction, default = true
Bool_t  fDoHadronicSignalRegion    = true; // get shapes for hadronic signal region <-- needed for background prediction, default = true
Bool_t  fDoPrediction              = true; // calculate the Znunu background from photon sample, needs fDoPhotonSignalRegion,fDoHadronicSignalRegion to be true, default = true
Bool_t  fPrintBrunoTable           = true; // gay printout - not needed anymore
Bool_t  fPrintFullPrintOut         = true;  // full printout during running
Bool_t  fMakeFinalTable            = true;  // make final prediction table
Bool_t  fDoVariableMT2bins         = true;  // variable MT2bins - needs to be true
Bool_t  fDoPileUpWeights           = true;  // apply PU weights - default = true
Bool_t  fbSFReWeight               = true;  // apply BTV SF weights --> need to adjust fNBJets = {-10,[-3,3]}, default = true
Bool_t  fMrennaHack                = true;  // Steve Mrenna Status 2 parton jets - prompt photon dR - this is used to remove overlap of gamma MC and QCD MC
Bool_t  fUseConstantZToGammaR      = false; // Constant Z/gamma ratio: see analysis note!! - this should be always false
Bool_t  fUseConstZToGammaRdynamic  = true; // Constant Z/gamma ratio: see analysis note!! <-- from fMinMT2forConstR on, before Z/gamma is bin-by-bin - default = true
Float_t fMinMT2forConstR           = 350;   // min MT2 value after which a constant Z/gamma ratio is used, needs fUseConstZToGammaRdynamic = true, default = 350 for 8 TeV
// uncertainty
Bool_t  fAddFitIntrinsicUncert     = false; // if you'd like to add an intrinsic uncertainty to the template fit
Float_t fFitIntrinsicUncert        = 0.0;
Bool_t  fAddRMCUncertainty         = true;  // systematic uncertainty on Z/gamma ratio
Float_t fRMCUncertainty            = 0.2;   // 0.3 for MT2 > 275, otherwise 0.2 (see AN)
Float_t fRMCUncertaintyTail        = 0.3;   // 0.3 for MT2 > 275, otherwise 0.2 (see AN)
// options ---------------
//TString fOutDir                    = "../Results/Filtered/GammaJetsPrediction/20130617_0b1b_test/";    // directory where shapes will be saved
TString fOutDir                    = "../dummy/";    // directory where shapes will be saved
Int_t   fVerbose                   = 9; 
Bool_t  fSaveResults               = true;  // default = true
Bool_t  fSaveZnunuToGammaRatio     = true;  // default = true
Bool_t  fDraw                      = true;  // draw histos or simply save them
Bool_t  fWriteToFile               = false; // writes couts to file, otherwise just prints it to terminal, default = false
//the 4 numbers below are a relict from 7 TeV analysis, but are used if fUseConstantZToGammaR = true, which might be interesting for debugging isues
Float_t fConstantZToGammaR_LowHT   = 0.458; // 750 < HT < 950
Float_t fConstantZToGammaErr_LowHT = 0.057; // abs uncertainty on fConstantZToGammaR the relative uncertainty fRMCUncertainty is added in quadrature if fAddRMCUncertainty==true
Float_t fConstantZToGammaR_HighHT  = 0.628; // HT > 950
Float_t fConstantZToGammaErr_HighHT= 0.120; // abs uncertainty on fConstantZToGammaR 
// ------
Bool_t  fHT                        = true;   //run over HT dataset
Bool_t  fMET                       = false;  //run over single photon dataset
Bool_t  fISRreweight               = false;  //apply an effective ISR weight (effective means not reweight gen-photon/Z but observed distribution), default = false
Bool_t  fReScaleAfterMT2           = false;  //rescales MC histograms to fit data normalization - as histograms are with MT2 cut, this means renormalizating after applying MT2 cut, default = false
Bool_t  fdontRescaleRatio          = true;   //keep this true, for Z/gamma calculation use photon MC not normalized to data in sigmaietaieta fit - idea is that both the gamma MC and Z MC should ne on equal footing

   //definition of all signal bins
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

//definition of topological regions
const int signalregionsize = 12;
string signal_region[signalregionsize] = {"2j0b", "2j1to2b", "2j1to2bmod", "3to5j0b", "3to5j1b", "3to5j1bmod", "3to5j2b", "6j0b", "6j1b", "6j1bmod", "6j2b", "3b"};
//definition of HT regions, the names are set by the fMET/fHT flags above
const int HTbinsize = 2;
string HT_bin[HTbinsize] = {"HTge450", "HTge750"};//dummy

// Some Global Variables ------------------------------------
std::ostringstream  fTriggerStream;
std::ostringstream  fTriggerStreamPhotons;
std::ostringstream  fCutStreamPhotons;
std::ostringstream  fCutStreamPhotonsMT2;
std::ostringstream  fCutStreamSignal;
std::ostringstream* fLogStream     = 0;
// cut imput for HT, NBJets, NJets, fSR - all are dummies now
Float_t fHTmin;
Float_t fHTmax;
Float_t fMT2min;
Float_t fMT2max;
Int_t   fNBJets;
Int_t   fNJets ;
Int_t   fSR;
Double_t fGammakFactor;
Double_t modtablescale;
Double_t modtablescaleerr;
Bool_t   fModTable;
//for final table ---------------------------------------------
vector<int>    fhtbin;
vector<int>    fmt2bin;
vector<int>    fNJ;
vector<int>    fNBJ;
vector<int>    fSigReg;
vector<double> fznunugen;
vector<double> fznunupred;
vector<double> fznunuprederr;
vector<double> fznunuprederrstat;
vector<double> fznunuprederrsyst;
vector<double> fZGratio;
vector<double> fZGratioerr;
vector<double> fmt2low;
vector<double> fmt2up;
vector<double> fndata;
vector<double> fndataerr;
vector<double> fnqcd;
vector<double> fnqcderr;
vector<double> fnother;
vector<double> fnothererr;
vector<double> fngamma;
vector<double> fngammaerr;
vector<double> fModV;
vector<double> fModVE;
vector<bool>   fMod;

void MakeFinalPredictionTable();

// this function calls all Cut Streams - needed as function to separate fHT/fMET cases
void DefineCutStreams(float HTmin, float HTmax, float MT2min, float MT2max, int NJets, int NBJets){
	fTriggerStream = 0;
	fTriggerStreamPhotons = 0;
	fCutStreamPhotons = 0;
	fCutStreamPhotonsMT2 = 0;
	fCutStreamSignal = 0;
	// Trigger Stream ---------------------------------------------------------------
	//
	if(fHT){
	fTriggerStream << "( ("
                       << "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || "
                       << "trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
                       << "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1)"
                       << "&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	fTriggerStreamPhotons << "( ("
                              << "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || "
                              << "trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
                              << "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1)"
                              << "&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	}
	if(fMET){
  	fTriggerStream << "( ("
                       << "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
                       << "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
	fTriggerStreamPhotons << "(trigger.HLT_SinglePhotons == 1&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0)";
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
	  << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"             << "&&" // for fastsim samples
	  << "misc.CSCTightHaloIDFlag == 0"                                << "&&"
	  << "misc.trackingFailureFlag==0"                                 << "&&"
	  << "misc.eeBadScFlag==0"                                         << "&&"
	  << "misc.EcalDeadCellTriggerPrimitiveFlag==0"                    << "&&"
	  << "misc.TrackingManyStripClusFlag==0"                           << "&&"
	  << "misc.TrackingTooManyStripClusFlag==0"                        << "&&"
	  << "misc.TrackingLogErrorTooManyClustersFlag==0"                 << "&&"
          << "misc.CrazyHCAL==0"
	  << "&&(type1pfmet[0].Pt()<30||type1pfmet[0].Pt()/misc.CaloMETRaw<=2.)";
	if(fMET) fCutStreamPhotons << "&& photon[0].lv.Pt()>=180 ";
	if (NJets>=10)   fCutStreamPhotons << "&&NJetsIDLoose40>=" << NJets/10 << "&&NJetsIDLoose40<=" << NJets%10;
	else if(NJets>0) fCutStreamPhotons << "&&NJetsIDLoose40==" << NJets;
	else             fCutStreamPhotons << "&&NJetsIDLoose40>=" << abs(NJets);
	if(NBJets>=0)        fCutStreamPhotons << "&&NBJets40CSVM==" << NBJets;
	else if(NBJets!=-10) fCutStreamPhotons << "&&NBJets40CSVM>=" << NBJets;

	if(!fMrennaHack){
		fCutStreamPhotons << "&&"	
	  	<< "(misc.ProcessID!=6||photon[0].MCmatchexitcode!=1)";
	}else{
		fCutStreamPhotons << "&&"	
		<< "(misc.ProcessID!=6||(photon[0].MCmatchexitcode!=1||photon[0].GenJetMinDR<0.3))"           << "&&"
		<< "(misc.ProcessID!=5||(photon[0].GenJetMinDR>0.3))";
	}

	// CutStream for Photon Signal Region ------------------------------------------ 
	fCutStreamPhotonsMT2 << " " 
	  //<< "misc.MT2>=" << gMT2bins[0]                                   << "&&"
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
	  << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"             << "&&" // for fastsim samples
	  << "misc.CSCTightHaloIDFlag == 0"                                << "&&"
	  << "misc.trackingFailureFlag==0"                                 << "&&"
	  << "misc.eeBadScFlag==0"                                         << "&&"
	  << "misc.EcalDeadCellTriggerPrimitiveFlag==0"                    << "&&"
	  << "misc.TrackingManyStripClusFlag==0"                           << "&&"
	  << "misc.TrackingTooManyStripClusFlag==0"                        << "&&"
	  << "misc.TrackingLogErrorTooManyClustersFlag==0"                 << "&&"
          << "misc.CrazyHCAL==0"
	  << "&&(type1pfmet[0].Pt()<30||type1pfmet[0].Pt()/misc.CaloMETRaw<=2.)";
	if(fMET){
		fCutStreamPhotonsMT2 << "&& photon[0].lv.Pt()>=180 && misc.MET>200";
	}
	if (NJets>=10)   fCutStreamPhotonsMT2 << "&&NJetsIDLoose40>=" << NJets/10 << "&&NJetsIDLoose40<=" << NJets%10;
	else if(NJets>0) fCutStreamPhotonsMT2 << "&&NJetsIDLoose40==" << NJets;
	else             fCutStreamPhotonsMT2 << "&&NJetsIDLoose40>=" << abs(NJets);
	if(NBJets>=0)        fCutStreamPhotonsMT2 << "&&NBJets40CSVM==" << NBJets;
	else if(NBJets!=-10) fCutStreamPhotonsMT2 << "&&NBJets40CSVM>=" << NBJets;
	if(!fMrennaHack){
		fCutStreamPhotonsMT2<< "&&"	
	  	<< "(misc.ProcessID!=6||photon[0].MCmatchexitcode!=1)";
	}else{
		fCutStreamPhotonsMT2<< "&&"	
		<< "(misc.ProcessID!=6||(photon[0].MCmatchexitcode!=1||photon[0].GenJetMinDR<0.3))"           << "&&"
		<< "(misc.ProcessID!=5||(photon[0].GenJetMinDR>0.3))";
	}

	// CutStream for Hadronic Signal Region ----------------------------------------------
	fCutStreamSignal << " " 
	  //<< "misc.MT2>=" << gMT2bins[0]                                   << "&&"
	  << "misc.MET>=30"                                                << "&&"
	  << "misc.HT >=" << HTmin  << " && misc.HT <=" << HTmax           << "&&"
	  << "misc.MT2>=" << MT2min << " && misc.MT2<=" << MT2max          << "&&"
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
	  << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"             << "&&" // for fastsim samples
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
		fCutStreamSignal << "&& misc.MET>200";
	}
	if (NJets>=10)   fCutStreamSignal << "&&NJetsIDLoose40>=" << NJets/10 << "&&NJetsIDLoose40<=" << NJets%10;
	else if(NJets>0) fCutStreamSignal << "&&NJetsIDLoose40==" << NJets;
	else             fCutStreamSignal << "&&NJetsIDLoose40>=" << abs(NJets);
	if(NBJets>=0)        fCutStreamSignal << "&&NBJets40CSVM==" << NBJets;
	else if(NBJets!=-10) fCutStreamSignal << "&&NBJets40CSVM>=" << NBJets;

}

// *********************** run_GammaJetsToZnunu: this is the main function calling also other functions ********************************
//this function calls all functions needed to do the Z(nunu) estimate from photon data sample, note that there is a lot of user input at the beginning
void try_GammaJetsToZnunu(){
	gSystem->Load("libPhysics");
	gSystem->Load("libRooFit") ;
	gSystem->CompileMacro("../MT2Code/src/MT2Shapes.cc", "k");//needed to efficiently get the shapes for MC and data along sigmaietaieta or MT2

	fhtbin.clear(), fmt2bin.clear(), fznunugen.clear(), fznunupred.clear(), fznunuprederr.clear(), fznunuprederrstat.clear(), fznunuprederrsyst.clear(), fZGratio.clear(), fZGratioerr.clear(), fmt2low.clear(), fmt2up.clear(), fndata.clear(), fndataerr.clear(), fnqcd.clear(), fnqcderr.clear(), fnother.clear(), fnothererr.clear(), fngamma.clear(), fngammaerr.clear(), fNJ.clear(); fNBJ.clear(), fMod.clear(), fModV.clear(), fModVE.clear(), fSigReg.clear();

	//at the moment samples.dat are hardcoded
	if(fMET) fSamplesRemovedPhotons     ="samples/samples_1g_GEst_MET_filter.dat";
	if(fMET) fSamplesHadronic           ="samples/samples_had_GEst_MET_filter.dat";
	if(fHT ) fSamplesRemovedPhotons     ="samples/samples_1g_GEst_HT_filter.dat";
	if(fHT ) fSamplesHadronic           ="samples/samples_had_GEst_HT_filter.dat";
                 fSamplesHadronic           ="samples/samples_Znunu_HTMET_filter.dat";//only Z(nunu) needed for hadronic region
	cout << "your samples are " << endl;
	cout << "    " << fSamplesRemovedPhotons << endl;
	cout << "    " << fSamplesHadronic << endl;

	//for fMET we need singlephotons with pT>180 GeV. There sigmaietaieta is not discriminant - don't do the fit - in samples.dat the kFactor are the 'hardcoded' normalization
	if(fMET) fDoPhotonSigmaIEtaIEta = false;
	if(!fDoPhotonSigmaIEtaIEta) cout << "YOUR ARE NOT DOING fDoPhotonSigmaIEtaIEta " << endl;
	//adjust k factor - for fHT the k factor of gamma is different as for fMET
	//this is the kFaktor to gamma MC from the samples.dat, is needed to have correct Z/G ratio
	if(fMET)      fGammakFactor = 0.9;
	else if (fHT) fGammakFactor = 0.8;
	cout << "photonMC k factor is " << fGammakFactor << endl;
	// logStream
	fLogStream = new std::ostringstream();
	TString OutputDirOrig = fOutDir;

	//we make two prediction tables
	//the second table uses the 0 b estimate and scales it with the Zll(1b)/Zll(0b) ratio obtained with functions like GammaVsZllStudies.C GammaVsZllStudiesRatios.C
	//the corresponding ratios and their uncertainties are hard-coded here
	for(int i1 = 0; i1<=HTbinsize; ++i1){
		if(fMET && i1==1) continue;
	for(int i2 = 0; i2<=signalregionsize; ++i2){
		int gNMT2bins;
		if(fMET){
			if(signal_region[i2]=="2j0b")       gNMT2bins = gNMT2bins_2j0b_lHT;
			if(signal_region[i2]=="2j1to2b")    gNMT2bins = gNMT2bins_2j1b_lHT;
			if(signal_region[i2]=="2j1to2bmod") { gNMT2bins = gNMT2bins_2j1b_lHT; modtablescale = 0.0937533; modtablescaleerr = 0.0224717;}
			if(signal_region[i2]=="3to5j0b")    gNMT2bins = gNMT2bins_3j0b_lHT;
			if(signal_region[i2]=="3to5j1b")    gNMT2bins = gNMT2bins_3j1b_lHT;
			if(signal_region[i2]=="3to5j1bmod") { gNMT2bins = gNMT2bins_3j1b_lHT; modtablescale = 0.16062;   modtablescaleerr = 0.0238675;}
			if(signal_region[i2]=="3to5j2b")    gNMT2bins = gNMT2bins_3j2b_lHT;
			if(signal_region[i2]=="6j0b")       gNMT2bins = gNMT2bins_6j0b_lHT;
			if(signal_region[i2]=="6j1b")       gNMT2bins = gNMT2bins_6j1b_lHT;
			if(signal_region[i2]=="6j1bmod")    { gNMT2bins = gNMT2bins_6j1b_lHT; modtablescale = 0.2693407; modtablescaleerr = 0.167893; }
			if(signal_region[i2]=="6j2b")       gNMT2bins = gNMT2bins_6j2b_lHT;
			if(signal_region[i2]=="3b")         gNMT2bins = gNMT2bins_3b_lHT;
		} if(fHT){
		   if(i1==0){
			if(signal_region[i2]=="2j0b")       gNMT2bins = gNMT2bins_2j0b_mHT;
			if(signal_region[i2]=="2j1to2b")    gNMT2bins = gNMT2bins_2j1b_mHT;
			if(signal_region[i2]=="2j1to2bmod") { gNMT2bins = gNMT2bins_2j1b_mHT; modtablescale = 0.0937533; modtablescaleerr = 0.0503878;}
			if(signal_region[i2]=="3to5j0b")    gNMT2bins = gNMT2bins_3j0b_mHT;
			if(signal_region[i2]=="3to5j1b")    gNMT2bins = gNMT2bins_3j1b_mHT;
			if(signal_region[i2]=="3to5j1bmod") { gNMT2bins = gNMT2bins_3j1b_mHT; modtablescale = 0.16062;   modtablescaleerr = 0.0181368;}
			if(signal_region[i2]=="3to5j2b")    gNMT2bins = gNMT2bins_3j2b_mHT;
			if(signal_region[i2]=="6j0b")       gNMT2bins = gNMT2bins_6j0b_mHT;
			if(signal_region[i2]=="6j1b")       gNMT2bins = gNMT2bins_6j1b_mHT;
			if(signal_region[i2]=="6j1bmod")    { gNMT2bins = gNMT2bins_6j1b_mHT; modtablescale = 0.2693407; modtablescaleerr = 0.142773; }
			if(signal_region[i2]=="6j2b")       gNMT2bins = gNMT2bins_6j2b_mHT;
			if(signal_region[i2]=="3b")         gNMT2bins = gNMT2bins_3b_mHT;
		   } if(i1==1){
			if(signal_region[i2]=="2j0b")       gNMT2bins = gNMT2bins_2j0b_hHT;
			if(signal_region[i2]=="2j1to2b")    gNMT2bins = gNMT2bins_2j1b_hHT;
			if(signal_region[i2]=="2j1to2bmod") { gNMT2bins = gNMT2bins_2j1b_hHT; modtablescale = 0.0937533; modtablescaleerr = 0.0597187;}
			if(signal_region[i2]=="3to5j0b")    gNMT2bins = gNMT2bins_3j0b_hHT;
			if(signal_region[i2]=="3to5j1b")    gNMT2bins = gNMT2bins_3j1b_hHT;
			if(signal_region[i2]=="3to5j1bmod") { gNMT2bins = gNMT2bins_3j1b_hHT; modtablescale = 0.16062;   modtablescaleerr = 0.0303366;}
			if(signal_region[i2]=="3to5j2b")    gNMT2bins = gNMT2bins_3j2b_hHT;
			if(signal_region[i2]=="6j0b")       gNMT2bins = gNMT2bins_6j0b_hHT;
			if(signal_region[i2]=="6j1b")       gNMT2bins = gNMT2bins_6j1b_hHT;
			if(signal_region[i2]=="6j1bmod")    { gNMT2bins = gNMT2bins_6j1b_hHT; modtablescale = 0.2693407; modtablescaleerr = 0.201307; }
			if(signal_region[i2]=="6j2b")       gNMT2bins = gNMT2bins_6j2b_hHT;
			if(signal_region[i2]=="3b")         gNMT2bins = gNMT2bins_3b_hHT;
		   }
		}
		const int NMT2bins = gNMT2bins+1;
  		double gMT2bins[NMT2bins];
		if(fMET){
			if(signal_region[i2]=="2j0b")       { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_2j0b_lHT[i0]; }
			if(signal_region[i2]=="2j1to2b")    { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_2j1b_lHT[i0]; }
			if(signal_region[i2]=="2j1to2bmod") { fModTable = true;  for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_2j1b_lHT[i0]; }
			if(signal_region[i2]=="3to5j0b")    { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3j0b_lHT[i0]; }
			if(signal_region[i2]=="3to5j1b")    { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3j1b_lHT[i0]; }
			if(signal_region[i2]=="3to5j1bmod") { fModTable = true;  for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3j1b_lHT[i0]; }
			if(signal_region[i2]=="3to5j2b")    { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3j2b_lHT[i0]; }
			if(signal_region[i2]=="6j0b")       { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_6j0b_lHT[i0]; }
			if(signal_region[i2]=="6j1b")       { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_6j1b_lHT[i0]; }
			if(signal_region[i2]=="6j1bmod")    { fModTable = true;  for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_6j1b_lHT[i0]; }
			if(signal_region[i2]=="6j2b")       { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_6j2b_lHT[i0]; }
			if(signal_region[i2]=="3b")         { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3b_lHT[i0];   }
		} if(fHT){
		   if(i1==0){
			if(signal_region[i2]=="2j0b")       { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_2j0b_mHT[i0];  }
			if(signal_region[i2]=="2j1to2b")    { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_2j1b_mHT[i0];  }
			if(signal_region[i2]=="2j1to2bmod") { fModTable = true;  for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_2j1b_mHT[i0];  }
			if(signal_region[i2]=="3to5j0b")    { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3j0b_mHT[i0];  }
			if(signal_region[i2]=="3to5j1b")    { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3j1b_mHT[i0];  }
			if(signal_region[i2]=="3to5j1bmod") { fModTable = true;  for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3j1b_mHT[i0];  }
			if(signal_region[i2]=="3to5j2b")    { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3j2b_mHT[i0];  }
			if(signal_region[i2]=="6j0b")       { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_6j0b_mHT[i0];  }
			if(signal_region[i2]=="6j1b")       { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_6j1b_mHT[i0];  }
			if(signal_region[i2]=="6j1bmod")    { fModTable = true;  for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_6j1b_mHT[i0];  }
			if(signal_region[i2]=="6j2b")       { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_6j2b_mHT[i0];  }
			if(signal_region[i2]=="3b")         { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3b_mHT[i0];    }
		   } if(i1==1){
			if(signal_region[i2]=="2j0b")       { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_2j0b_hHT[i0];  }
			if(signal_region[i2]=="2j1to2b")    { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_2j1b_hHT[i0];  }
			if(signal_region[i2]=="2j1to2bmod") { fModTable = true;  for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_2j1b_hHT[i0];  }
			if(signal_region[i2]=="3to5j0b")    { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3j0b_hHT[i0];  }
			if(signal_region[i2]=="3to5j1b")    { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3j1b_hHT[i0];  }
			if(signal_region[i2]=="3to5j1bmod") { fModTable = true;  for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3j1b_hHT[i0];  }
			if(signal_region[i2]=="3to5j2b")    { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3j2b_hHT[i0];  }
			if(signal_region[i2]=="6j0b")       { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_6j0b_hHT[i0];  }
			if(signal_region[i2]=="6j1b")       { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_6j1b_hHT[i0];  }
			if(signal_region[i2]=="6j1bmod")    { fModTable = true;  for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_6j1b_hHT[i0];  }
			if(signal_region[i2]=="6j2b")       { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_6j2b_hHT[i0];  }
			if(signal_region[i2]=="3b")         { fModTable = false; for(int i0 = 0; i0<=gNMT2bins; ++i0) gMT2bins[i0] = gMT2bins_3b_hHT[i0];    }
		   }
		}
		if(fModTable==false) { modtablescale = 1.0; modtablescaleerr = 0.0; }
		if(fMET) {             fHTmin =  450.; fHTmax =  750.;}
		if(fHT)  { if(i1==0) { fHTmin =  750.; fHTmax = 1200.;}  
		           if(i1==1) { fHTmin = 1200.; fHTmax =  1E+8;} }
		fMT2min = gMT2bins[0]; fMT2max = 1E+8;
		//set topological region, NJets, NBJets
		if(signal_region[i2]=="2j0b")       { fNJets =  2; fNBJets =  0; fSR = 0; }
		if(signal_region[i2]=="2j1to2b")    { fNJets =  2; fNBJets = -1; fSR = 1; }
		if(signal_region[i2]=="2j1to2bmod") { fNJets =  2; fNBJets =  0; fSR = 1; }//fake NBJets
		if(signal_region[i2]=="3to5j0b")    { fNJets = 35; fNBJets =  0; fSR = 2; }
		if(signal_region[i2]=="3to5j1b")    { fNJets = 35; fNBJets =  1; fSR = 3; }
		if(signal_region[i2]=="3to5j1bmod") { fNJets = 35; fNBJets =  0; fSR = 3; }//fake NBJets
		if(signal_region[i2]=="3to5j2b")    { fNJets = 35; fNBJets =  2; fSR = 4; }
		if(signal_region[i2]=="6j0b")       { fNJets = -6; fNBJets =  0; fSR = 5; }
		if(signal_region[i2]=="6j1b")       { fNJets = -6; fNBJets =  1; fSR = 6; }
		if(signal_region[i2]=="6j1bmod")    { fNJets = -6; fNBJets =  0; fSR = 6; }//fake NBJets
		if(signal_region[i2]=="6j2b")       { fNJets = -6; fNBJets =  2; fSR = 7; }
		if(signal_region[i2]=="3b")         { fNJets = -3; fNBJets = -3; fSR = 8; }

		// define cutsteams
		DefineCutStreams(fHTmin, fHTmax, fMT2min, fMT2max, fNJets, fNBJets);
	
		fOutDir = OutputDirOrig;
		// fix output dir
		if(fReScaleAfterMT2)  fOutDir= TString::Format("%s_%s",       "RescaledMT2",   fOutDir.Data());
		if(fISRreweight)      fOutDir= TString::Format("%s_%s",       "ISRreweighted", fOutDir.Data());
		if(fNJets <0)         fOutDir= TString::Format("%s_ge%dj",     fOutDir.Data(), abs(fNJets));
		else                  fOutDir= TString::Format("%s_%dj",       fOutDir.Data(), fNJets);
		if(fNBJets>=0)        fOutDir= TString::Format("%s_%db",       fOutDir.Data(), fNBJets);
		else if(fNBJets!=-10) fOutDir= TString::Format("%s_ge%db",     fOutDir.Data(), abs(fNBJets));
		else                  fOutDir= TString::Format("%s_%d_HT_%s",  fOutDir.Data(), abs(fHTmin),  "Inf");
		if(fHTmax <10000)     fOutDir= TString::Format("%s_%d_HT_%d",  fOutDir.Data(), abs(fHTmin),  abs(fHTmax));
		else                  fOutDir= TString::Format("%s_%d_HT_%s",  fOutDir.Data(), abs(fHTmin),  "Inf");
		if(fMT2max<10000)     fOutDir= TString::Format("%s_%d_MT2_%d", fOutDir.Data(), abs(fMT2min), abs(fMT2max));
		else                  fOutDir= TString::Format("%s_%d_MT2_%s", fOutDir.Data(), abs(fMT2min), "Inf");
	
		// log MT2 and HT cuts
		*fLogStream << "------------------------------------------------------------------------------------------------" << endl;
		*fLogStream << "+++ new Znunu with Gamma+jets prediction                                                     +++" << endl;
		*fLogStream << "+++ outputdir: " << fOutDir <<                                                              "+++" << endl; 
		*fLogStream << "------------------------------------------------------------------------------------------------" << endl;
	
		// new prediction class ------------------------------------------------------
		Prediction* prediction = new Prediction();
		prediction->fVerbose=fVerbose;
		prediction->fSave   =fSaveResults;
		prediction->fOutputDir=fOutDir;
	
	
		// SigmaIEtaIEta Fit *********************************************************************
		// Get Photon Normalization: EB or EB+EE
		if(fDoPhotonSigmaIEtaIEta){
			std::ostringstream cutStreamPhotons_EB;
			std::ostringstream cutStreamPhotons_EE;
			if(fSeparateEBEE){
				cutStreamPhotons_EB << fCutStreamPhotons.str() << "&&abs(photon[0].lv.Eta())<1.4442";
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
			//use hard-coded normalization
			prediction->SetPhotonNormalizationStupid(prediction->MLRes_EB);
			prediction->SetPhotonNormalizationStupid(prediction->MLRes_EE);
		}
	
		// Get Photon signal yield (in data and MC)
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
	
		// get hadronic Znunu yield - note: you should modify the samples.dat to contain only Znunu MC, as rest is not needed
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
			delete prediction->HadronicSignalRegion;
		}
		if(fDoPhotonSignalRegion){//clean up
			delete prediction->PhotonicSignalRegion_EB;
			delete prediction->PhotonicSignalRegion_EE;
		}
		delete prediction;
	}}//topological,HT

	if(fMakeFinalTable) MakeFinalPredictionTable();//uses only global variables
	
	if(fWriteToFile){
		TString logname =fOutDir + ".log"; 
		ofstream f_log (logname.Data(), ios::app);
		f_log << fLogStream->str();
	} else{
		cout << fLogStream->str();
	}
	delete fLogStream;
}


// ****************************************** Definition of classes and methods *********************************************************
// class to contain results of maximum likelihood fit to sigmaietaieta variable
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

//channels are sigmaietaieta region, photon signal region, hadronic region (x2 for EE and EB)
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

//get shapes is used to obtain the shapes along sigmaietaieta or MT2 for given selection and binning
void Channel::GetShapes(TString SelectionName, TString xtitle, const int nbins, const double *bins ){ 
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

	bool issigmaietaieta = false;
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
		if(name.Contains("SigmaIEtaIEta"))      issigmaietaieta = true;
	}
	delete tA;
	fGotShapes=true;

	if(fISRreweight){
		//reweight only photon jets and Znunu via an effective scaling - this is stupidly coded as not important
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
		if(fPrintFullPrintOut){
			for(int nbin = 1; nbin<=hPhotons->GetNbinsX(); ++nbin){//do it like that in order to not change the error
				cout << "Photon " << hPhotons->GetBinContent(nbin) << " Zinv " << hZJetsToNuNu->GetBinContent(nbin) << endl;
				cout << "SF G   " << hG->GetBinContent(nbin) << " SF Z " << hZ->GetBinContent(nbin) << endl;
				hPhotons->SetBinContent(nbin, hPhotons->GetBinContent(nbin)*hG->GetBinContent(nbin));
				hPhotons->SetBinError(  nbin, hPhotons->GetBinError(  nbin)*hG->GetBinContent(nbin));
				hZJetsToNuNu->SetBinContent(nbin, hZJetsToNuNu->GetBinContent(nbin)*hZ->GetBinContent(nbin));
				hZJetsToNuNu->SetBinError(  nbin, hZJetsToNuNu->GetBinError(  nbin)*hZ->GetBinContent(nbin));
				cout << "Photon " << hPhotons->GetBinContent(nbin) << " Zinv " << hZJetsToNuNu->GetBinContent(nbin) << endl;
			}
		}
		f->Close();
	}
	if(issigmaietaieta){
		if(hOther==0){
			hOther=(TH1D*)hQCD->Clone("Other");
			if(hOther->Integral()>1.) hOther->Scale(0.01/hOther->Integral());//total hOther yield is <= 0.01;
			else                      hOther->Scale(0.01);
		} else if(hOther->Integral()<=0){
			hOther=(TH1D*)hQCD->Clone("Other");
			if(hOther->Integral()>1.) hOther->Scale(0.01/hOther->Integral());//total hOther yield is <= 0.01;
			else                      hOther->Scale(0.01);
		}
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
//help function if binning is given differently
void Channel::GetShapes(TString SelectionName, TString xtitle, const int nbins, float binmin, float binmax){ 
	double bins[nbins];
	bins[0] = binmin;
	for(int i=1; i<=nbins; i++) bins[i] = binmin+i*(binmax-binmin)/nbins;
	GetShapes(SelectionName, xtitle, nbins, bins);
}

//this class is there to really do the predictions - combines channels, fit results and important functions
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
	TH1D* RescaleHisto(TH1D* histMC1, TH1D* histMC2, TH1D* hist_data_EB, TH1D* hist_data_EE);
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

//this function gets the photon normalization (performs sigmaietaieta fit) by using the sigmaietaieta channels and filling results in MLResult
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
	if(hOther->Integral()<=0){
		hOther->Add(hQCD);
		if(hOther->Integral()>1.) hOther->Scale(0.01/hOther->Integral());//total hOther yield is <= 0.01;
		else                      hOther->Scale(0.01);
	}

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

//if sigmaietaieta fit is not performed hardcode all scale factors
//values obtained by an extrapolation - i.e. this might change for 13 TeV
void Prediction::SetPhotonNormalizationStupid(MLResult* MLRes){

	//these are the necessary information for the prediction
	MLRes->fMLPhotonScaleFactor    = 1.2;
	MLRes->fMLPhotonScaleFactorErr = 0.05;
	MLRes->fMLQCDScaleFactor       = 1.3;
	MLRes->fMLQCDScaleFactorErr    = 1.3;
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

// compute MC Znunu/Photon ratio
void Prediction::GetMCZnunuToPhotonRatio(){
	// for this PhotonicSignalRegion->fGotShapes==1 and HadronicSignalRegion->fGotShapes==1 is required!
	if(PhotonicSignalRegion_EE->fGotShapes==false || PhotonicSignalRegion_EB->fGotShapes==false || HadronicSignalRegion->fGotShapes==false){
		cerr << "ERROR in GetMCZnunuToPhotonRatio: fGotShape found to be false" << endl;
		exit(-1);
	}
	*fLogStream << "------------------------------------------------------------------------" << endl;
	*fLogStream << "GetMCZnunuToPhotonRatio: Computing MC Znunu to Photon Ratio             " << endl;           
	if(fUseConstantZToGammaR){//dummy thing to do - but might be useful for debugging reasons
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
	}else{  //this is the important part of this function - get ratio from photon MC and Znunu MC
		if(PhotonicSignalRegion_EB->hPhotons==0 || PhotonicSignalRegion_EE->hPhotons==0 || HadronicSignalRegion->hZJetsToNuNu==0){
			cout << "GetMCZnunuToPhotonRatio: received 0 pointer!" << endl;
			exit(-1);
		}

		// GetScaled histos with only one bin and propagated errors: stat error and error on scale factor
		double rescaleGammaEB = MLRes_EB->fMLPhotonScaleFactor; double rescaleGammaEBerr = MLRes_EB->fMLPhotonScaleFactorErr;
		double rescaleGammaEE = MLRes_EE->fMLPhotonScaleFactor; double rescaleGammaEEerr = MLRes_EE->fMLPhotonScaleFactorErr;
		double rescaleQCDEB   = MLRes_EB->fMLQCDScaleFactor;    double rescaleQCDEBerr   = MLRes_EB->fMLQCDScaleFactorErr;
		double rescaleQCDEE   = MLRes_EE->fMLQCDScaleFactor;    double rescaleQCDEEerr   = MLRes_EE->fMLQCDScaleFactorErr;
		if(fdontRescaleRatio){
			rescaleGammaEB = 1./fGammakFactor; rescaleGammaEBerr = 0.;
			rescaleGammaEE = 1./fGammakFactor; rescaleGammaEEerr = 0.;
			rescaleQCDEB   = 1.; rescaleQCDEBerr   = 0.;
			rescaleQCDEE   = 1.; rescaleQCDEEerr   = 0.;
		}
		//the GetScaledHisto2 keeps the binning and scales the histogram with a factor (which here is obtained from the sigmaietaieta fit)
		TH1D *currPhotons= GetScaledHisto2(PhotonicSignalRegion_EB->hPhotons , rescaleGammaEB, rescaleGammaEBerr); // EB
		currPhotons->Add(  GetScaledHisto2(PhotonicSignalRegion_EE->hPhotons , rescaleGammaEE, rescaleGammaEEerr)); //EE
		if(fReScaleAfterMT2){
			TH1D *currQCD= GetScaledHisto2(PhotonicSignalRegion_EB->hQCD , rescaleQCDEB, rescaleQCDEBerr); // EB
			currQCD->Add(  GetScaledHisto2(PhotonicSignalRegion_EE->hQCD , rescaleQCDEE, rescaleQCDEEerr)); //EE
			currPhotons = RescaleHisto(currPhotons, currQCD, PhotonicSignalRegion_EB->hData, PhotonicSignalRegion_EE->hData);
		}
		TH1D *currZnunu  = GetScaledHisto2(HadronicSignalRegion->hZJetsToNuNu, 1                   , 0);

		//constant Z/G ratio for high MT2 is achieved by setting the bins of the input histograms with fMinMT2forConstR to the sum of those bins
		if(fUseConstZToGammaRdynamic) currPhotons = RefillRatioHisto(currPhotons,fMinMT2forConstR);
		if(fUseConstZToGammaRdynamic) currZnunu   = RefillRatioHisto(currZnunu  ,fMinMT2forConstR);

		//get the ratio
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
		for(int i = 1; i<= fMCZnunuPhotonRatioHisto->GetNbinsX(); ++i) {
			//either bin is above threshold, last bin, or bin is just below threshold, but majority of the bin is in increased uncertainty window
			if(fMCZnunuPhotonRatioHisto->GetBinLowEdge(i)>=(fMinMT2forConstR-0.01) || i==fMCZnunuPhotonRatioHisto->GetNbinsX() || 
			  (fMCZnunuPhotonRatioHisto->GetBinLowEdge(i+1)>=(fMinMT2forConstR-0.01) && fMCZnunuPhotonRatioHisto->GetBinLowEdge(i)<(fMinMT2forConstR-0.01) &&
			  (fMCZnunuPhotonRatioHisto->GetBinLowEdge(i+1)-fMinMT2forConstR)>=(fMinMT2forConstR-fMCZnunuPhotonRatioHisto->GetBinLowEdge(i) ) ) )
				fMCZnunuPhotonRatioHisto->SetBinError(i, sqrt(pow(fMCZnunuPhotonRatioHisto->GetBinError(i),2)+pow(fMCZnunuPhotonRatioHisto->GetBinContent(i)*fRMCUncertaintyTail,2)));
			else
				fMCZnunuPhotonRatioHisto->SetBinError(i, sqrt(pow(fMCZnunuPhotonRatioHisto->GetBinError(i),2)+pow(fMCZnunuPhotonRatioHisto->GetBinContent(i)*fRMCUncertainty,2)));

		}
		}//fAddRMCUncertainty

		*fLogStream << "Ratio: " << fMCZnunuPhotonRatio << " pm " << fMCZnunuPhotonRatioErr << endl;
		*fLogStream << "------------------------------------------------------------------------" << endl;

		if(fSaveZnunuToGammaRatio){
			//the GetScaledHisto2 keeps the binning and scales the histogram with a factor (which here is obtained from the sigmaietaieta fit)
			TH1D *currPhotons2= GetScaledHisto2(PhotonicSignalRegion_EB->hPhotons , rescaleGammaEB, rescaleGammaEBerr); // EB
			currPhotons2->Add(  GetScaledHisto2(PhotonicSignalRegion_EE->hPhotons , rescaleGammaEE, rescaleGammaEEerr)); //EE
			if(fReScaleAfterMT2){
				TH1D *currQCD= GetScaledHisto2(PhotonicSignalRegion_EB->hQCD , rescaleQCDEB, rescaleQCDEBerr); // EB
				currQCD->Add(  GetScaledHisto2(PhotonicSignalRegion_EE->hQCD , rescaleQCDEE, rescaleQCDEEerr)); //EE
				currPhotons2 = RescaleHisto(currPhotons2, currQCD, PhotonicSignalRegion_EB->hData, PhotonicSignalRegion_EE->hData);
			}
			//this GetScaledHisto rebins histogram by nothing(dummy) and scales the histogram with a factor (which here is obtained from the sigmaietaieta fit)
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

//after obtaining all normalizations and the Z/G ratio can perform final prediction
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
	
	//the GetScaledHisto2 keeps the binning and scales the histogram with a factor (which here is obtained from the sigmaietaieta fit)
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

	if(fReScaleAfterMT2){
		//this is discouraged - meaning  see at the beginning at fReScaleAfterMT2 setting
		TH1D *temp = hPhotons->Clone("hPhotons_Clone");
		hPhotons = RescaleHisto(hPhotons, hQCD, PhotonicSignalRegion_EB->hData, PhotonicSignalRegion_EE->hData);
		hQCD     = RescaleHisto(hQCD    , temp, PhotonicSignalRegion_EB->hData, PhotonicSignalRegion_EE->hData);
	}

	*fLogStream << "Photonic Signal Region: event yields (only first bin(:-------------------" << endl;
	*fLogStream << "EB:                                                     " << endl;
	*fLogStream << "  NData:    "    << hData_EB   ->GetBinContent(1) << " pm " << hData_EB   ->GetBinError(1)  << endl;
	*fLogStream << "  NPhotons: "    << hPhotons_EB->GetBinContent(1) << " pm " << hPhotons_EB->GetBinError(1)  << endl;
	*fLogStream << "  NQCD:     "    << hQCD_EB    ->GetBinContent(1) << " pm " << hQCD_EB    ->GetBinError(1)  << endl;
	*fLogStream << "  NOther:   "    << hOther_EB  ->GetBinContent(1) << " pm " << hOther_EB  ->GetBinError(1)  << endl;
	*fLogStream << "EE:                                                     " << endl;
	*fLogStream << "  NData:    "    << hData_EE   ->GetBinContent(1) << " pm " << hData_EE   ->GetBinError(1)  << endl;
	*fLogStream << "  NPhotons: "    << hPhotons_EE->GetBinContent(1) << " pm " << hPhotons_EE->GetBinError(1)  << endl;
	*fLogStream << "  NQCD:     "    << hQCD_EE    ->GetBinContent(1) << " pm " << hQCD_EE    ->GetBinError(1)  << endl;
	*fLogStream << "  NOther:   "    << hOther_EE  ->GetBinContent(1) << " pm " << hOther_EE  ->GetBinError(1)  << endl;
	*fLogStream << "total:                                                  " << endl;
	*fLogStream << "  NData:    "    << hData   ->GetBinContent(1)    << " pm " << hData   ->GetBinError(1)  << endl;
	*fLogStream << "  NPhotons: "    << hPhotons->GetBinContent(1)    << " pm " << hPhotons->GetBinError(1)  << endl;
	*fLogStream << "  NQCD:     "    << hQCD    ->GetBinContent(1)    << " pm " << hQCD    ->GetBinError(1)  << endl;
	*fLogStream << "  NOther:   "    << hOther  ->GetBinContent(1)    << " pm " << hOther  ->GetBinError(1)  << endl;
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
	

	//obtain the predicted yields - Pascal's style (only first bin)
	float PredictedZnunu                = fMCZnunuPhotonRatio*(hData->GetBinContent(1)-hOther->GetBinContent(1)-hQCD->GetBinContent(1));
	float PredictedZnunu_ErrSys         = GetZnunuPreditionErrorSys (hData, hOther, hQCD);
	float PredictedZnunu_ErrSysClosure  = GetZnunuPreditionErrorSysClosure (hPhotons);
	float PredictedZnunu_ErrStat        = GetZnunuPreditionErrorStat(hData);

	//obtained the predicted yields along MT2
	vector<int> htbin, mt2bin;
	vector<double> znunugen, znunupred, znunuprederr, znunuprederrstat, znunuprederrsyst, ZGratio, ZGratioerr, mt2low, mt2up, ndata, ndataerr, nqcd, nqcderr, nother, nothererr, ngamma, ngammaerr;
	htbin.clear(), mt2bin.clear(), znunugen.clear(), znunupred.clear(), znunuprederr.clear(), znunuprederrstat.clear(), znunuprederrsyst.clear(), ZGratio.clear(), ZGratioerr.clear(), mt2low.clear(), mt2up.clear(), ndata.clear(), ndataerr.clear(), nqcd.clear(), nqcderr.clear(), nother.clear(), nothererr.clear(), ngamma.clear(), ngammaerr.clear();
	for(int i = 1; i<=hData->GetNbinsX(); ++i){

	PredictedZnunu                = fMCZnunuPhotonRatioHisto->GetBinContent(i)*(hData->GetBinContent(i)-hOther->GetBinContent(i)-hQCD->GetBinContent(i));
	PredictedZnunu_ErrSys         = GetZnunuPreditionErrorSys (i, hData, hOther, hQCD);
	PredictedZnunu_ErrSysClosure  = GetZnunuPreditionErrorSysClosure (i, hPhotons);
	PredictedZnunu_ErrStat        = GetZnunuPreditionErrorStat(i, hData);

	//make here a if case - this a detailed printout 
	*fLogStream << "For " << hPhotons->GetBinLowEdge(i) << " < MT2 < " << hPhotons->GetBinLowEdge(i)+hPhotons->GetBinWidth(i) <<  "(" << i << "/" << hData->GetNbinsX() << ")" << endl;

	*fLogStream << "Photon events in Photon Signal region (data-bg): " <<  hData->GetBinContent(i) << " - " << hOther->GetBinContent(i)+hQCD->GetBinContent(i) << " MC: " << hPhotons->GetBinContent(i) << endl;
	*fLogStream << "Predicted N Znunu events in Hadronic Signal region: ratio " << fMCZnunuPhotonRatioHisto->GetBinContent(i) << " pm " <<fMCZnunuPhotonRatioHisto->GetBinError(i) << endl;
	*fLogStream << "Prediction: "                                                                                            << endl;
	*fLogStream << "  " << PredictedZnunu << " pm " << PredictedZnunu_ErrSys  << " pm " << PredictedZnunu_ErrStat << " stat "        << endl;
	*fLogStream << "True N Znunu events:"                                                                                    << endl;
	*fLogStream << "  " << HadronicSignalRegion->hZJetsToNuNu->GetBinContent(i)                                                           << endl;
	*fLogStream << "MC closure:"                                                                                              << endl;
	*fLogStream << "  " << fMCZnunuPhotonRatioHisto->GetBinContent(i)*(hPhotons->GetBinContent(i))   << " pm " << PredictedZnunu_ErrSysClosure << " (sys) " << endl; 

	int MT2bin = 99;
	int HTbin  = 99;
	if      (fMT2min == 150 && fMT2max == 200) MT2bin=0;
	else if (fMT2min == 200 && fMT2max == 275) MT2bin=1;
	else if (fMT2min == 275 && fMT2max == 375) MT2bin=2;
	else if (fMT2min == 375 && fMT2max == 500) MT2bin=3;
	else if (fMT2min == 500)                   MT2bin=4;
	else if (fMT2min ==   0 && fMT2max ==1e+8) MT2bin=-1;
	//else    {cout << "MT2bin not valid for printout! " << endl; exit(-1);} //this is stupid
	if      (fHTmin  == 750 && fHTmax == 1200) HTbin =1;
	else if (fHTmin  == 1200)                  HTbin =2;
	else if (fHTmin  == 450 && fHTmax == 750 ) HTbin =0;
	else if (fHTmin  == 750 )                  HTbin =10;
	else if (fHTmin  == 450 )                  HTbin =-10;
	//else    {cout << "HTbin not valid for printout! " << endl; exit(-1);} //this is stupid
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
	ngamma.push_back(hPhotons->GetBinContent(i));
	ngammaerr.push_back(hPhotons->GetBinError(i));

	fhtbin.push_back(HTbin);
	fmt2bin.push_back(i);
	fznunugen.push_back(HadronicSignalRegion->hZJetsToNuNu->GetBinContent(i));
	fznunupred.push_back(PredictedZnunu);
	fznunuprederr.push_back(sqrt(PredictedZnunu_ErrSys*PredictedZnunu_ErrSys+PredictedZnunu_ErrStat*PredictedZnunu_ErrStat));
	fznunuprederrstat.push_back(PredictedZnunu_ErrStat);
	fznunuprederrsyst.push_back(PredictedZnunu_ErrSys);
	fZGratio.push_back(fMCZnunuPhotonRatioHisto->GetBinContent(i));
	fZGratioerr.push_back(fMCZnunuPhotonRatioHisto->GetBinError(i));
	fmt2low.push_back(hPhotons->GetBinLowEdge(i));
	if(i!=hPhotons->GetNbinsX()) fmt2up.push_back(hPhotons->GetBinLowEdge(i)+hPhotons->GetBinWidth(i));
	else                         fmt2up.push_back(9999999.);
	fndata.push_back(hData->GetBinContent(i));
	fndataerr.push_back(hData->GetBinError(i));
	fnqcd.push_back(hQCD->GetBinContent(i));
	fnqcderr.push_back(hQCD->GetBinError(i));
	fnother.push_back(hOther->GetBinContent(i));
	fnothererr.push_back(hOther->GetBinError(i));
	fngamma.push_back(hPhotons->GetBinContent(i));
	fngammaerr.push_back(hPhotons->GetBinError(i));
	fNJ.push_back(fNJets);
	fNBJ.push_back(fNBJets);
	fMod.push_back(modtablescale);
	fModV.push_back(modtablescaleerr);
	fModVE.push_back(fModTable);
	fSigReg.push_back(fSR);

	if(fPrintBrunoTable){//see fPrintBrunoTable at beginning
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
			    << HadronicSignalRegion->hZJetsToNuNu->GetBinContent(i) 
			    <<  "after fit " << hPhotons->GetBinContent(i)
			    << " before fit " << PhotonicSignalRegion_EB->hPhotons->GetBinContent(i)+ PhotonicSignalRegion_EE->hPhotons->GetBinContent(i)
			    << " \\\\ " << endl;
		*fLogStream << "************************** Bruno-Gay Printout ******************************" << endl;
	}
	}
	if(fPrintBrunoTable){//see fPrintBrunoTable at beginning
		*fLogStream << endl << "************************** Bruno-Gay Printout ******************************" << endl;
        	*fLogStream << fixed << setprecision(2) ;
		*fLogStream << " Name     " << " SR " << " HTbin " << " MT2bin   " << " MCpred   " << " DataDrivenPred " << " DataDrivenPredError " << " ScaleFactor " << " ScaleFactorError " /*<< " Ngamma"*/ << endl;
		for(int jj = 0; jj<htbin.size(); ++jj){
			*fLogStream << "ZinvFromG   " << fSR << " " << setw(5) << htbin[jj] << " " << setw(4) << mt2bin[jj] << " " << setw(13) << znunugen[jj]
			            << " " << setw(13) << znunupred[jj] << " " << setw(18) << znunuprederr[jj]
		    		    << " " << setw(13) << ZGratio[jj]   << " " << setw(16) << ZGratioerr[jj] /*<<  "            " << ngamma[jj] << "+/-" << ngammaerr[jj]*/ << endl;
		}
		*fLogStream << "----------------------------------------------------------------------------" << endl;
	}

	//this is the final prediction table for Z(nunu) from photon
	*fLogStream << endl;
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
		            << " " << setw(19) << "& " << int(ndata[jj])   << " " << setw(1)  << "& " 
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


	*fLogStream << endl << "************************** Modified table for 1b est ******************************" << endl;
        *fLogStream << fixed << setprecision(2) ;
	*fLogStream << " Name     " << " SR " << " HTbin " << " MT2bin   " << " MCpred(wrong)   " << " DataDrivenPred " << " DataDrivenPredError(slightwrong) " << " ScaleFactor " << " ScaleFactorError " << endl;
	for(int jj = 0; jj<htbin.size(); ++jj){
		*fLogStream << "ZinvFromG   " << fSR << " " << setw(5) << htbin[jj] << " " << setw(4) << mt2bin[jj] << " " << setw(13) << znunugen[jj]
		            << " " << setw(13) << znunupred[jj]*modtablescale << " " << setw(18) << sqrt(pow(znunuprederr[jj]*modtablescale,2)+pow(znunupred[jj]*modtablescaleerr,2))
		    	    << " " << setw(13) << ZGratio[jj]*modtablescale   << " " << setw(16) << sqrt(pow(ZGratioerr[jj]*modtablescale,2) + pow(ZGratio[jj]*modtablescaleerr,2))  << endl;
	}
		*fLogStream << "----------------------------------------------------------------------------" << endl;
	*fLogStream << "\%BEGINLATEX\%"                    << endl
		    << "\\begin{table}[!htb]"              << endl
		    << "\\begin{center}"                   << endl
		    << "\\begin{tabular}{l|cccccc}"        << endl
		    << "\\hline\\hline"                    << endl
		    << "\\multirow{2}{*}{$M_{T2}$\\ bin}";
		    if(fHTmin==750 && fHTmax==950) *fLogStream << " & \\multicolumn{6}{c}{$750 < H_T<950$ GeV}  \\";
	*fLogStream << "  & $N^{\\gamma}$         & $QCD^{bkg}$&           $R_Z^{data}(1b/0b)$      &    $R^{MC}(Z(\\nu\\nu)/\\gamma)$ & data pred                   & MC estimate(is 0b not 1b) \\\\" << endl
		    << "\\hline" << endl;
	for(int jj = 0; jj<htbin.size(); ++jj){
		*fLogStream << "$" << int(mt2low[jj]) << "-" << int(mt2up[jj]) << "$" 
		            << " " << setw(19) << "& " << int(ndata[jj]) << " " << setw(1)  << "& " 
		            << " " << setw(2)  <<   "$" << nqcd[jj]       << " \\pm "<< nqcderr[jj]    << "$" << " " << setw(1)  << "& "
		            << " " << setw(2) <<  fixed << setprecision(3) <<   "$" << modtablescale     << " \\pm "<< modtablescaleerr  << "$" << " " << setw(1)  << "& "
		            << " " << setw(4) <<  fixed << setprecision(2) <<   "$" << ZGratio[jj]    << " \\pm "<< ZGratioerr[jj] << "$" << " " << setw(1)  << "& "
		            << " " << setw(7)  <<   "$" << znunupred[jj]*modtablescale  << " \\pm " << znunuprederrstat[jj]*modtablescale << " \\pm " << sqrt(pow(znunuprederrsyst[jj]*modtablescale,2)+pow(znunupred[jj]*modtablescaleerr,2)) << "$ & "
			    << " " << setw(4)  << znunugen[jj] << " \\\\ " << endl;
	}
	*fLogStream << "\\hline\\hline"  << endl
		    << "\\end{tabular}"  << endl
		    << "\\end{center}"   << endl
		    << "\\end{table}"    << endl
		    << "\%ENDLATEX\%"    << endl
		    << endl;
	*fLogStream << endl << "************************** Modified table for 1b est ******************************" << endl;

	if(!fDoPhotonSigmaIEtaIEta) cout << "YOUR ARE NOT DOING fDoPhotonSigmaIEtaIEta " << endl;//just a reminder

}

TH1D* Prediction::RebinHisto(TH1D* histo, float lastbinlowedge){
	//rebin a histogram in such a way that all bins < lastbinlowedge are kept as they are
	//but all bins > lastbinlowedge are rebinned into one single bin
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
	//this function adds up the contents and uncertainties for bins > lastbinlowedge
	//and replaces those contents and uncertainties by their (quadratic) sum
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
	// and returns scaled histo with uncertainty on scale factor propagated.
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
	// and returns scaled and rebinned histo with ngroup bins merged into 1 bin with uncertainty on scale factor propagated.
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
	float pred_err_2 = pow(fMCZnunuPhotonRatioErr * (hData->GetBinContent(1) - hOther->GetBinContent(1) - hQCD->GetBinContent(1)),2);
	pred_err_2      += pow(fMCZnunuPhotonRatio * hOther->GetBinError(1),2);
	pred_err_2      += pow(fMCZnunuPhotonRatio * hQCD->GetBinError(1),2);

	return sqrt(pred_err_2);
}

float Prediction::GetZnunuPreditionErrorSys(int i, TH1D* hData, TH1D* hOther, TH1D* hQCD){
	// histograms must be scales with prober uncertainties
	if(i>hData->GetNbinsX() || i<=0) return -999.;
	float pred_err_2 = pow(fMCZnunuPhotonRatioHisto->GetBinError(i)  * (hData->GetBinContent(i) - hOther->GetBinContent(i) - hQCD->GetBinContent(i)),2);
	pred_err_2      += pow(fMCZnunuPhotonRatioHisto->GetBinContent(i) * hOther->GetBinError(i),2);
	pred_err_2      += pow(fMCZnunuPhotonRatioHisto->GetBinContent(i) * hQCD->GetBinError(i),2);

	return sqrt(pred_err_2);
}

float Prediction::GetZnunuPreditionErrorSysClosure(TH1D* hPhotons){
	// computes systematic uncerainty on Znunu prediction for MC closure test:
	// takes into account MLfit uncertainty on Photon-MC yield, MC statistics, statistical uncerainty on Znunu/photon ratio
	float pred_err_2 = pow(fMCZnunuPhotonRatioErr* hPhotons->GetBinContent(1),2);
	pred_err_2      += pow(fMCZnunuPhotonRatio   * hPhotons->GetBinError(1)  ,2);
	return sqrt(pred_err_2);
}

float Prediction::GetZnunuPreditionErrorSysClosure(int i, TH1D* hPhotons){
	// computes systematic uncerainty on Znunu prediction for MC closure test:
	// takes into account MLfit uncertainty on Photon-MC yield, MC statistics, statistical uncerainty on Znunu/photon ratio
	if(i>hPhotons->GetNbinsX() || i<=0) return -999.;
	float pred_err_2 = pow(fMCZnunuPhotonRatioHisto->GetBinError(i)     * hPhotons->GetBinContent(i),2);
	pred_err_2      += pow(fMCZnunuPhotonRatioHisto->GetBinContent(i)   * hPhotons->GetBinError(i)  ,2);
	return sqrt(pred_err_2);
}

float Prediction::GetZnunuPreditionErrorStat(TH1D* hData){
	// statistical uncertainty on prediction due to data statistics
	float pred_err   = fMCZnunuPhotonRatio * hData->GetBinError(1);
	return pred_err;
}

float Prediction::GetZnunuPreditionErrorStat(int i, TH1D* hData){
	// statistical uncertainty on prediction due to data statistics
	if(i>hData->GetNbinsX() || i<=0) return -999.;
	float pred_err   = fMCZnunuPhotonRatioHisto->GetBinContent(i) * hData->GetBinError(i);
	return pred_err;
}

TH1D* Prediction::RescaleHisto(TH1D* histMC1, TH1D* histMC2, TH1D* hist_data_EB, TH1D* hist_data_EE){
	// takes histo, scale factor and uncertainty on scale factor
	// and returns scaled and rebinned histo with 1 bin with uncertainty on scale factor propagated.
	TH1D *hData    = (TH1D*) hist_data_EB->Clone("hData");
	TH1D *hDataEE  = (TH1D*) hist_data_EE->Clone("hData_EE");
	TH1D *h        = (TH1D*) histMC1->Clone(histMC1->GetName());
	TH1D *hMC      = (TH1D*) histMC1->Clone("hMC");
	TH1D *htemp    = (TH1D*) histMC2->Clone("htemp");
	hData->Add(hDataEE);
	hMC  ->Add(htemp);
	double scalefactor = 1;
	if(hMC->Integral()>0&&hData->Integral()>0) scalefactor = hData->Integral()/hMC->Integral();
	if(fVerbose>4) cout << "after MT2>" << hData->GetBinLowEdge(1) << ": data " << hData->Integral() << ", MC " << hMC->Integral() << " --> SF = " << scalefactor << endl;
	h->Scale(scalefactor);
	return h;
}

void DrawHisto(TH1* h_orig, TString name,  Option_t *drawopt, Channel* channel){
	//draws histograms
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


void MakeFinalPredictionTable(){

	if(fPrintBrunoTable){//see fPrintBrunoTable at beginning
		*fLogStream << endl << "************************** Bruno-Gay Printout ******************************" << endl;
        	*fLogStream << fixed << setprecision(2) ;
		*fLogStream << " Name     " << " SR " << " HTbin " << " MT2bin   " << " MCpred   " << " DataDrivenPred " << " DataDrivenPredError " << " ScaleFactor " << " ScaleFactorError " /*<< " Ngamma"*/ << endl;
		for(int jj = 0; jj<fhtbin.size(); ++jj){
			if(fMod[jj]==true) continue;
			*fLogStream << "ZinvFromG   " << fSigReg[jj] << " " << setw(5) << fhtbin[jj] << " " << setw(4) << fmt2bin[jj] << " " << setw(13) << fznunugen[jj]
			            << " " << setw(13) << fznunupred[jj] << " " << setw(18) << fznunuprederr[jj]
		    		    << " " << setw(13) << fZGratio[jj]   << " " << setw(16) << fZGratioerr[jj] /*<<  "            " << ngamma[jj] << "+/-" << ngammaerr[jj]*/ << endl;
		}
		*fLogStream << "----------------------------------------------------------------------------" << endl;
	}

	//this is the final prediction table for Z(nunu) from photon
	*fLogStream << endl;
	*fLogStream << "\%BEGINLATEX\%"                    << endl
		    << "\\begin{table}[!htb]"              << endl
		    << "\\begin{center}"                   << endl
		    << "\\begin{tabular}{r|cccccc}"        << endl
		    << "\\hline\\hline"                    << endl;
	*fLogStream << "$M_\\mathrm{T2}$ [GeV] & $N^{\\gamma}_\\mathrm{data}$ & $N^{QCD}_\\mathrm{bkg}$ &           $N^{EWK}^\\mathrm{bkg}$      &    $R_\\mathrm{sim}(Z(\\nu\\bar{\\nu})/\\gamma)$ & data prediction             & sim. truth \\\\" << endl
		    << "\\hline\\hline" << endl;
	int oldHT = -999; int oldSR = -999;
	for(int jj = 0; jj<fhtbin.size(); ++jj){
		if(fMod[jj]==true) continue;
		if(fhtbin[jj]!=oldHT){
		    oldHT = fhtbin[jj];
	            if(fhtbin[jj]==0) *fLogStream << " & \\multicolumn{6}{c}{low H_\\mathrm{T}} \\\\ \\hline\\hline" << endl;
	            if(fhtbin[jj]==1) *fLogStream << " & \\multicolumn{6}{c}{medium H_\\mathrm{T}} \\\\ \\hline\\hline" << endl;
	            if(fhtbin[jj]==2) *fLogStream << " & \\multicolumn{6}{c}{high H_\\mathrm{T}} \\\\ \\hline\\hline" << endl;
		}
		if(fSigReg[jj]!=oldSR){
		   oldSR = fSigReg[jj];
		   if(fSigReg[jj]==0) *fLogStream << " & \\multicolumn{6}{c}{2 jets, 0 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==1) *fLogStream << " & \\multicolumn{6}{c}{2 jets, $\\geq1$ b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==2) *fLogStream << " & \\multicolumn{6}{c}{3-5 jets, 0 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==3) *fLogStream << " & \\multicolumn{6}{c}{3-5 jets, 1 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==4) *fLogStream << " & \\multicolumn{6}{c}{3-5 jets, 2 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==5) *fLogStream << " & \\multicolumn{6}{c}{$\\geq6$ jets, 0 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==6) *fLogStream << " & \\multicolumn{6}{c}{$\\geq6$ jets, 1 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==7) *fLogStream << " & \\multicolumn{6}{c}{$\\geq6$ jets, 2 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==8) *fLogStream << " & \\multicolumn{6}{c}{$\\geq3$ jets, $\\geq3$ b jets} \\\\ \\hline" << endl;
		};
		if(fmt2up<100000.)  *fLogStream << "$" << int(fmt2low[jj]) << "-" << int(fmt2up[jj]) << "$" 
		else                *fLogStream << "$\\geq " << int(fmt2low[jj]) << "$";
		                    *fLogStream << " " << setw(19) << "& " << int(fndata[jj])   << " " << setw(1)  << "& ";
		if(fnqcd[jj]>0)     *fLogStream << <<  fixed << setprecision(2) << " " << setw(2)  <<   "$" << fnqcd[jj]       << " \\pm "<< fnqcderr[jj]    << "$" << " " << setw(1)  << "& ";
		else                *fLogStream << " " << setw(8) << "$-$" << " " << setw(7)  << "& ";
		if(fnother[jj]>0)   *fLogStream << " " << setw(2)  <<   "$" << fnother[jj]     << " \\pm "<< fnothererr[jj]  << "$" << " " << setw(1)  << "& ";
		else                *fLogStream << " " << setw(8) << "$-$" << " " << setw(7)  << "& ";
		                    *fLogStream << " " << setw(4)  <<   "$" << fZGratio[jj]    << " \\pm "<< fZGratioerr[jj] << "$" << " " << setw(1)  << "& "
		                                << " " << setw(7)  <<   "$" << fznunupred[jj]  << " \\pm " << fznunuprederrstat[jj] << " \\pm " << fznunuprederrsyst[jj] << "$ & "
			                        << " " << setw(4)  << fznunugen[jj] << " \\\\ " << endl;
	}
	*fLogStream << "\\hline\\hline"  << endl
		    << "\\end{tabular}"  << endl
		    << "\\end{center}"   << endl
		    << "\\end{table}"    << endl
		    << "\%ENDLATEX\%"    << endl
		    << endl;
	if(!fDoPhotonSigmaIEtaIEta) cout << "reminder: YOUR ARE NOT DOING fDoPhotonSigmaIEtaIEta " << endl;//just a reminder

	*fLogStream << endl << "************************** Modified table for 1b est ******************************" << endl;
	if(fPrintBrunoTable){//see fPrintBrunoTable at beginning
		*fLogStream << fixed << setprecision(2) ;
		*fLogStream << " Name     " << " SR " << " HTbin " << " MT2bin   " << " MCpred(wrong)   " << " DataDrivenPred " << " DataDrivenPredError(slightwrong) " << " ScaleFactor " << " ScaleFactorError " << endl;
		for(int jj = 0; jj<fhtbin.size(); ++jj){
			if(fMod[jj]==false) continue;
			*fLogStream << "ZinvFromG   " << fSigReg[jj] << " " << setw(5) << fhtbin[jj] << " " << setw(4) << fmt2bin[jj] << " " << setw(13) << fznunugen[jj]
				<< " " << setw(13) << fznunupred[jj]*fModV[jj] << " " << setw(18) << sqrt(pow(fznunuprederr[jj]*fModV[jj],2)+pow(fznunupred[jj]*fModVE[jj],2))
				<< " " << setw(13) << fZGratio[jj]*fModV[jj]   << " " << setw(16) << sqrt(pow(fZGratioerr[jj]*fModV[jj],2) + pow(fZGratio[jj]*fModVE[jj],2))  << endl;
		}
	}
	*fLogStream << "----------------------------------------------------------------------------" << endl << endl;

	*fLogStream << "$M_\\mathrm{T2}$ [GeV] & $N^{\\gamma}_\\mathrm{data}$ & $N^{QCD}_\\mathrm{bkg}$ &    $R^Z_\\mathrm{data}(1b/0b)$   &    $R_\\mathrm{sim}(Z(\\nu\\bar{\\nu})/\\gamma)$ & data prediction             & sim. truth \\\\" << endl
		    << "\\hline\\hline" << endl;
	oldHT = -999; oldSR = -999;
	for(int jj = 0; jj<fhtbin.size(); ++jj){
		if(fMod[jj]==false) continue;
		if(fhtbin[jj]!=oldHT){
		    oldHT = fhtbin[jj];
	            if(fhtbin[jj]==0) *fLogStream << " & \\multicolumn{6}{c}{low H_\\mathrm{T}} \\\\ \\hline\\hline" << endl;
	            if(fhtbin[jj]==1) *fLogStream << " & \\multicolumn{6}{c}{medium H_\\mathrm{T}} \\\\ \\hline\\hline" << endl;
	            if(fhtbin[jj]==2) *fLogStream << " & \\multicolumn{6}{c}{high H_\\mathrm{T}} \\\\ \\hline\\hline" << endl;
		}
		if(fSigReg[jj]!=oldSR){
		   oldSR = fSigReg[jj];
		   if(fSigReg[jj]==0) *fLogStream << " & \\multicolumn{6}{c}{2 jets, 0 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==1) *fLogStream << " & \\multicolumn{6}{c}{2 jets, $\\geq1$ b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==2) *fLogStream << " & \\multicolumn{6}{c}{3-5 jets, 0 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==3) *fLogStream << " & \\multicolumn{6}{c}{3-5 jets, 1 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==4) *fLogStream << " & \\multicolumn{6}{c}{3-5 jets, 2 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==5) *fLogStream << " & \\multicolumn{6}{c}{$\\geq6$ jets, 0 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==6) *fLogStream << " & \\multicolumn{6}{c}{$\\geq6$ jets, 1 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==7) *fLogStream << " & \\multicolumn{6}{c}{$\\geq6$ jets, 2 b jets} \\\\ \\hline" << endl;
		   if(fSigReg[jj]==8) *fLogStream << " & \\multicolumn{6}{c}{$\\geq3$ jets, $\\geq3$ b jets} \\\\ \\hline" << endl;
		};
		if(fmt2up<100000.)  *fLogStream << "$" << int(fmt2low[jj]) << "-" << int(fmt2up[jj]) << "$" 
		else                *fLogStream << "$\\geq " << int(fmt2low[jj]) << "$";
		                    *fLogStream << " " << setw(19) << "& " << int(fndata[jj])   << " " << setw(1)  << "& ";
		if(fnqcd[jj]>0)     *fLogStream << <<  fixed << setprecision(2) << " " << setw(2)  <<   "$" << fnqcd[jj]       << " \\pm "<< fnqcderr[jj]    << "$" << " " << setw(1)  << "& ";
		else                *fLogStream << " " << setw(8) << "$-$" << " " << setw(7)  << "& ";
		                    *fLogStream << " " << setw(2) <<  fixed << setprecision(3) <<   "$" << fModV[jj]     << " \\pm "<< fModVE[jj]  << "$" << " " << setw(1)  << "& ";
		                    *fLogStream << " " << setw(4) <<  fixed << setprecision(2) <<   "$" << fZGratio[jj]    << " \\pm "<< fZGratioerr[jj] << "$" << " " << setw(1)  << "& "
		                    *fLogStream << " " << setw(7)  <<   "$" << fznunupred[jj]*fModV[jj]  << " \\pm " << fznunuprederrstat[jj]*fModV[jj] << " \\pm " << sqrt(pow(fznunuprederrsyst[jj]*fModV[jj],2)+pow(fznunupred[jj]*fModVE[jj],2)) 
                                                << "$ & " << " " << setw(4)  << fznunugen[jj] << " \\\\ " << endl;
	}
	*fLogStream << "\\hline\\hline"  << endl
		    << "\\end{tabular}"  << endl
		    << "\\end{center}"   << endl
		    << "\\end{table}"    << endl
		    << "\%ENDLATEX\%"    << endl
		    << endl;
	*fLogStream << endl << "************************** Modified table for 1b est ******************************" << endl;

	if(!fDoPhotonSigmaIEtaIEta) cout << "reminder: YOUR ARE NOT DOING fDoPhotonSigmaIEtaIEta " << endl;//just a reminder

}//MakeFinalPredictionTable()
