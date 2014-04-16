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
#include "TChain.h"
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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path of MT2tree.cc
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/Utilities.hh"

//use via root -l -b -q TOBTECtest.C++

using namespace std;

//this code is obsolete (we have own TOBTEC noise filter), therefore no comments
//this code tested correlations among variables used for RA2b's TOBTEC filter (and also HO filter)
void TOBTECtest(){

gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

gROOT->cd();
TChain *data = new TChain("MassTree");
data->Add("~mmasciov/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/highHT/HT*.root");
data->Add("~mmasciov/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/highHT/JetHT*.root");
data->Add("~mmasciov/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/MET*.root");

gROOT->cd();

TH1D *hPFMEToverCaloMET = new TH1D("hPFMEToverCaloMET","",50,0,10); hPFMEToverCaloMET->SetTitle(";PFMET/CaloMET; Events");
TH1D *hPFMEToverCaloMETMuCorr = new TH1D("hPFMEToverCaloMETMuCorr","",50,0,10); hPFMEToverCaloMETMuCorr->SetTitle(";PFMET/CaloMET(#mu-corr); Events");
TH2D *hPFMETvsCaloMET = new TH2D("hPFMETvsCaloMET","",100,0,1500,100,0,1500); hPFMETvsCaloMET->SetTitle(";PFMET [GeV];CaloMET [GeV]");
TH2D *hPFMETvsCaloMETMuCorr = new TH2D("hPFMETvsCaloMETMuCorr","",100,0,1500,100,0,1500); hPFMETvsCaloMETMuCorr->SetTitle(";PFMET [GeV];CaloMET(#mu-corr) [GeV]");
TH1D *hMT2PFMEToverCaloMET = new TH1D("hMT2PFMEToverCaloMET","",50,0,10); hMT2PFMEToverCaloMET->SetTitle(";PFMET/CaloMET; Events");
TH1D *hMT2PFMEToverCaloMETMuCorr = new TH1D("hMT2PFMEToverCaloMETMuCorr","",50,0,10); hMT2PFMEToverCaloMETMuCorr->SetTitle(";PFMET/CaloMET(#mu-corr); Events");
TH2D *hMT2PFMETvsCaloMET = new TH2D("hMT2PFMETvsCaloMET","",100,0,1500,100,0,1500); hMT2PFMETvsCaloMET->SetTitle(";PFMET [GeV];CaloMET [GeV]");
TH2D *hMT2PFMETvsCaloMETMuCorr = new TH2D("hMT2PFMETvsCaloMETMuCorr","",100,0,1500,100,0,1500); hMT2PFMETvsCaloMETMuCorr->SetTitle(";PFMET [GeV];CaloMET(#mu-corr) [GeV]");

TH2D *hChvsNeu = new TH2D("hChvsNeu","",40,0,200,20,0,100); hChvsNeu->SetTitle(";N_{charged};N_{neutral}");
TH2D *hChvsAll = new TH2D("hChvsAll","",40,0,200,60,0,300); hChvsAll->SetTitle(";N_{charged};N_{charged}+N_{neutral}");
TH2D *hNeuvsAll = new TH2D("hNeuvsAll","",20,0,100,60,0,300); hNeuvsAll->SetTitle(";N_{neutral};N_{charged}+N_{neutral}");
TH2D *hChMNeuvsAll = new TH2D("hChMNeuvsAll","",40,0,200,60,0,300); hChMNeuvsAll->SetTitle(";N_{charged}-N_{neutral};N_{charged}+N_{neutral}");
TH2D *hChvsMass = new TH2D("hChvsMass","",40,0,200,60,0,300); hChvsMass->SetTitle(";N_{charged};Mass [GeV]");
TH2D *hAllvsMass = new TH2D("hAllvsMass","",60,0,300,60,0,300); hAllvsMass->SetTitle(";N_{charged}+N_{neutral};Mass [GeV]");
TH2D *hChMNeuvsMass = new TH2D("hChMNeuvsMass","",40,0,200,60,0,300); hChMNeuvsMass->SetTitle(";N_{charged}-N_{neutral};Mass [GeV]");
TH2D *hChvsMETDPhi = new TH2D("hChvsMETDPhi","",40,0,200,16,0,3.2); hChvsMETDPhi->SetTitle(";N_{charged};MET-#Delta#phi");
TH2D *hAllvsMETDPhi = new TH2D("hAllvsMETDPhi","",60,0,300,16,0,3.2); hAllvsMETDPhi->SetTitle(";N_{charged}+N_{neutral};MET-#Delta#phi");
TH2D *hChMNeuvsMETDPhi = new TH2D("hChMNeuvsMETDPhi","",40,0,200,16,0,3.2); hChMNeuvsMETDPhi->SetTitle(";N_{charged}-N_{neutral};MET-#Delta#phi");
TH2D *hChvsCHF = new TH2D("hChvsCHF","",40,0,200,20,0,1); hChvsCHF->SetTitle(";N_{charged};CHF");
TH2D *hNeuvsCHF = new TH2D("hNeuvsCHF","",20,0,100,20,0,1); hNeuvsCHF->SetTitle(";N_{neutral};CHF");
TH2D *hAllvsCHF = new TH2D("hAllvsCHF","",60,0,300,20,0,1); hAllvsCHF->SetTitle(";N_{charged}+N_{neutral};CHF");
TH2D *hChMNeuvsCHF = new TH2D("hChMNeuvsCHF","",40,0,200,20,0,1); hChMNeuvsCHF->SetTitle(";N_{charged}-N_{neutral};CHF");
TH2D *hMassvsCHF = new TH2D("hMassvsCHF","",60,0,300,20,0,1); hMassvsCHF->SetTitle(";Mass [GeV];CHF");
TH2D *hMETDPhivsCHF = new TH2D("hMETDPhivsCHF","",16,0,3.2,20,0,1); hMETDPhivsCHF->SetTitle(";MET-#Delta#phi;CHF");
TH2D *hMETDPhivsMass = new TH2D("hMETDPhivsMass","",16,0,3.2,60,0,300); hMETDPhivsMass->SetTitle(";MET-#Delta#phi;Mass [GeV]");

TH2D *hMT2ChvsNeu = new TH2D("hMT2ChvsNeu","",40,0,200,20,0,100); hMT2ChvsNeu->SetTitle(";N_{charged};N_{neutral}");
TH2D *hMT2ChvsAll = new TH2D("hMT2ChvsAll","",40,0,200,60,0,300); hMT2ChvsAll->SetTitle(";N_{charged};N_{charged}+N_{neutral}");
TH2D *hMT2NeuvsAll = new TH2D("hMT2NeuvsAll","",20,0,100,60,0,300); hMT2NeuvsAll->SetTitle(";N_{neutral};N_{charged}+N_{neutral}");
TH2D *hMT2ChMNeuvsAll = new TH2D("hMT2ChMNeuvsAll","",40,0,200,60,0,300); hMT2ChMNeuvsAll->SetTitle(";N_{charged}-N_{neutral};N_{charged}+N_{neutral}");
TH2D *hMT2ChvsMass = new TH2D("hMT2ChvsMass","",40,0,200,60,0,300); hMT2ChvsMass->SetTitle(";N_{charged};Mass [GeV]");
TH2D *hMT2AllvsMass = new TH2D("hMT2AllvsMass","",60,0,300,60,0,300); hMT2AllvsMass->SetTitle(";N_{charged}+N_{neutral};Mass [GeV]");
TH2D *hMT2ChMNeuvsMass = new TH2D("hMT2ChMNeuvsMass","",40,0,200,60,0,300); hMT2ChMNeuvsMass->SetTitle(";N_{charged}-N_{neutral};Mass [GeV]");
TH2D *hMT2ChvsMETDPhi = new TH2D("hMT2ChvsMETDPhi","",40,0,200,16,0,3.2); hMT2ChvsMETDPhi->SetTitle(";N_{charged};MET-#Delta#phi");
TH2D *hMT2AllvsMETDPhi = new TH2D("hMT2AllvsMETDPhi","",60,0,300,16,0,3.2); hMT2AllvsMETDPhi->SetTitle(";N_{charged}+N_{neutral};MET-#Delta#phi");
TH2D *hMT2ChMNeuvsMETDPhi = new TH2D("hMT2ChMNeuvsMETDPhi","",40,0,200,16,0,3.2); hMT2ChMNeuvsMETDPhi->SetTitle(";N_{charged}-N_{neutral};MET-#Delta#phi");
TH2D *hMT2ChvsCHF = new TH2D("hMT2ChvsCHF","",40,0,200,20,0,1); hMT2ChvsCHF->SetTitle(";N_{charged};CHF");
TH2D *hMT2NeuvsCHF = new TH2D("hMT2NeuvsCHF","",20,0,100,20,0,1); hMT2NeuvsCHF->SetTitle(";N_{neutral};CHF");
TH2D *hMT2AllvsCHF = new TH2D("hMT2AllvsCHF","",60,0,300,20,0,1); hMT2AllvsCHF->SetTitle(";N_{charged}+N_{neutral};CHF");
TH2D *hMT2ChMNeuvsCHF = new TH2D("hMT2ChMNeuvsCHF","",40,0,200,20,0,1); hMT2ChMNeuvsCHF->SetTitle(";N_{charged}-N_{neutral};CHF");
TH2D *hMT2MassvsCHF = new TH2D("hMT2MassvsCHF","",60,0,300,20,0,1); hMT2MassvsCHF->SetTitle(";Mass [GeV];CHF");
TH2D *hMT2METDPhivsCHF = new TH2D("hMT2METDPhivsCHF","",16,0,3.2,20,0,1); hMT2METDPhivsCHF->SetTitle(";MET-#Delta#phi;CHF");
TH2D *hMT2METDPhivsMass = new TH2D("hMT2METDPhivsMass","",16,0,3.2,60,0,300); hMT2METDPhivsMass->SetTitle(";MET-#Delta#phi;Mass [GeV]");

TH2D *hChvsEta = new TH2D("hChvsEta","",40,0,200,50,-5,5); hChvsEta->SetTitle(";N_{charged};#eta");
TH2D *hNeuvsEta = new TH2D("hNeuvsEta","",20,0,100,50,-5,5); hNeuvsEta->SetTitle(";N_{neutral};#eta");
TH2D *hAllvsEta = new TH2D("hAllvsEta","",60,0,300,50,-5,5); hAllvsEta->SetTitle(";N_{charged}+N_{neutral};#eta");
TH2D *hChMNeuvsEta = new TH2D("hChMNeuvsEta","",40,0,200,50,-5,5); hChMNeuvsEta->SetTitle(";N_{charged}-N_{neutral};#eta");
TH2D *hMassvsEta = new TH2D("hMassvsEta","",60,0,300,50,-5,5); hMassvsEta->SetTitle(";Mass [GeV];#eta");
TH2D *hMETDPhivsEta = new TH2D("hMETDPhivsEta","",16,0,3.2,50,-5,5); hMETDPhivsEta->SetTitle(";MET-#Delta#phi;#eta");
TH2D *hCHFvsEta = new TH2D("hCHFvsEta","",20,0,1,50,-5,5); hCHFvsEta->SetTitle(";CHF;#eta");

TH2D *hMT2ChvsEta = new TH2D("hMT2ChvsEta","",40,0,200,50,-5,5); hMT2ChvsEta->SetTitle(";N_{charged};#eta");
TH2D *hMT2NeuvsEta = new TH2D("hMT2NeuvsEta","",20,0,100,50,-5,5); hMT2NeuvsEta->SetTitle(";N_{neutral};#eta");
TH2D *hMT2AllvsEta = new TH2D("hMT2AllvsEta","",60,0,300,50,-5,5); hMT2AllvsEta->SetTitle(";N_{charged}+N_{neutral};#eta");
TH2D *hMT2ChMNeuvsEta = new TH2D("hMT2ChMNeuvsEta","",40,0,200,50,-5,5); hMT2ChMNeuvsEta->SetTitle(";N_{charged}-N_{neutral};#eta");
TH2D *hMT2MassvsEta = new TH2D("hMT2MassvsEta","",60,0,300,50,-5,5); hMT2MassvsEta->SetTitle(";Mass [GeV];#eta");
TH2D *hMT2METDPhivsEta = new TH2D("hMT2METDPhivsEta","",16,0,3.2,50,-5,5); hMT2METDPhivsEta->SetTitle(";MET-#Delta#phi;#eta");
TH2D *hMT2CHFvsEta = new TH2D("hMT2CHFvsEta","",20,0,1,50,-5,5); hMT2CHFvsEta->SetTitle(";CHF;#eta");

TH2D *hChvsNeu_etacut = new TH2D("hChvsNeu_etacut","",40,0,200,20,0,100); hChvsNeu_etacut->SetTitle(";N_{charged};N_{neutral}");
TH2D *hChvsAll_etacut = new TH2D("hChvsAll_etacut","",40,0,200,60,0,300); hChvsAll_etacut->SetTitle(";N_{charged};N_{charged}+N_{neutral}");
TH2D *hNeuvsAll_etacut = new TH2D("hNeuvsAll_etacut","",20,0,100,60,0,300); hNeuvsAll_etacut->SetTitle(";N_{neutral};N_{charged}+N_{neutral}");
TH2D *hChMNeuvsAll_etacut = new TH2D("hChMNeuvsAll_etacut","",40,0,200,60,0,300); hChMNeuvsAll_etacut->SetTitle(";N_{charged}-N_{neutral};N_{charged}+N_{neutral}");
TH2D *hChvsMass_etacut = new TH2D("hChvsMass_etacut","",40,0,200,60,0,300); hChvsMass_etacut->SetTitle(";N_{charged};Mass [GeV]");
TH2D *hAllvsMass_etacut = new TH2D("hAllvsMass_etacut","",60,0,300,60,0,300); hAllvsMass_etacut->SetTitle(";N_{charged}+N_{neutral};Mass [GeV]");
TH2D *hChMNeuvsMass_etacut = new TH2D("hChMNeuvsMass_etacut","",40,0,200,60,0,300); hChMNeuvsMass_etacut->SetTitle(";N_{charged}-N_{neutral};Mass [GeV]");
TH2D *hChvsMETDPhi_etacut = new TH2D("hChvsMETDPhi_etacut","",40,0,200,16,0,3.2); hChvsMETDPhi_etacut->SetTitle(";N_{charged};MET-#Delta#phi");
TH2D *hAllvsMETDPhi_etacut = new TH2D("hAllvsMETDPhi_etacut","",60,0,300,16,0,3.2); hAllvsMETDPhi_etacut->SetTitle(";N_{charged}+N_{neutral};MET-#Delta#phi");
TH2D *hChMNeuvsMETDPhi_etacut = new TH2D("hChMNeuvsMETDPhi_etacut","",40,0,200,16,0,3.2); hChMNeuvsMETDPhi_etacut->SetTitle(";N_{charged}-N_{neutral};MET-#Delta#phi");
TH2D *hChvsCHF_etacut = new TH2D("hChvsCHF_etacut","",40,0,200,20,0,1); hChvsCHF_etacut->SetTitle(";N_{charged};CHF");
TH2D *hNeuvsCHF_etacut = new TH2D("hNeuvsCHF_etacut","",20,0,100,20,0,1); hNeuvsCHF_etacut->SetTitle(";N_{neutral};CHF");
TH2D *hAllvsCHF_etacut = new TH2D("hAllvsCHF_etacut","",60,0,300,20,0,1); hAllvsCHF_etacut->SetTitle(";N_{charged}+N_{neutral};CHF");
TH2D *hChMNeuvsCHF_etacut = new TH2D("hChMNeuvsCHF_etacut","",40,0,200,20,0,1); hChMNeuvsCHF_etacut->SetTitle(";N_{charged}-N_{neutral};CHF");
TH2D *hMassvsCHF_etacut = new TH2D("hMassvsCHF_etacut","",60,0,300,20,0,1); hMassvsCHF_etacut->SetTitle(";Mass [GeV];CHF");
TH2D *hMETDPhivsCHF_etacut = new TH2D("hMETDPhivsCHF_etacut","",16,0,3.2,20,0,1); hMETDPhivsCHF_etacut->SetTitle(";MET-#Delta#phi;CHF");
TH2D *hMETDPhivsMass_etacut = new TH2D("hMETDPhivsMass_etacut","",16,0,3.2,60,0,300); hMETDPhivsMass_etacut->SetTitle(";MET-#Delta#phi;Mass [GeV]");

TH2D *hMT2ChvsNeu_etacut = new TH2D("hMT2ChvsNeu_etacut","",40,0,200,20,0,100); hMT2ChvsNeu_etacut->SetTitle(";N_{charged};N_{neutral}");
TH2D *hMT2ChvsAll_etacut = new TH2D("hMT2ChvsAll_etacut","",40,0,200,60,0,300); hMT2ChvsAll_etacut->SetTitle(";N_{charged};N_{charged}+N_{neutral}");
TH2D *hMT2NeuvsAll_etacut = new TH2D("hMT2NeuvsAll_etacut","",20,0,100,60,0,300); hMT2NeuvsAll_etacut->SetTitle(";N_{neutral};N_{charged}+N_{neutral}");
TH2D *hMT2ChMNeuvsAll_etacut = new TH2D("hMT2ChMNeuvsAll_etacut","",40,0,200,60,0,300); hMT2ChMNeuvsAll_etacut->SetTitle(";N_{charged}-N_{neutral};N_{charged}+N_{neutral}");
TH2D *hMT2ChvsMass_etacut = new TH2D("hMT2ChvsMass_etacut","",40,0,200,60,0,300); hMT2ChvsMass_etacut->SetTitle(";N_{charged};Mass [GeV]");
TH2D *hMT2AllvsMass_etacut = new TH2D("hMT2AllvsMass_etacut","",60,0,300,60,0,300); hMT2AllvsMass_etacut->SetTitle(";N_{charged}+N_{neutral};Mass [GeV]");
TH2D *hMT2ChMNeuvsMass_etacut = new TH2D("hMT2ChMNeuvsMass_etacut","",40,0,200,60,0,300); hMT2ChMNeuvsMass_etacut->SetTitle(";N_{charged}-N_{neutral};Mass [GeV]");
TH2D *hMT2ChvsMETDPhi_etacut = new TH2D("hMT2ChvsMETDPhi_etacut","",40,0,200,16,0,3.2); hMT2ChvsMETDPhi_etacut->SetTitle(";N_{charged};MET-#Delta#phi");
TH2D *hMT2AllvsMETDPhi_etacut = new TH2D("hMT2AllvsMETDPhi_etacut","",60,0,300,16,0,3.2); hMT2AllvsMETDPhi_etacut->SetTitle(";N_{charged}+N_{neutral};MET-#Delta#phi");
TH2D *hMT2ChMNeuvsMETDPhi_etacut = new TH2D("hMT2ChMNeuvsMETDPhi_etacut","",40,0,200,16,0,3.2); hMT2ChMNeuvsMETDPhi_etacut->SetTitle(";N_{charged}-N_{neutral};MET-#Delta#phi");
TH2D *hMT2ChvsCHF_etacut = new TH2D("hMT2ChvsCHF_etacut","",40,0,200,20,0,1); hMT2ChvsCHF_etacut->SetTitle(";N_{charged};CHF");
TH2D *hMT2NeuvsCHF_etacut = new TH2D("hMT2NeuvsCHF_etacut","",20,0,100,20,0,1); hMT2NeuvsCHF_etacut->SetTitle(";N_{neutral};CHF");
TH2D *hMT2AllvsCHF_etacut = new TH2D("hMT2AllvsCHF_etacut","",60,0,300,20,0,1); hMT2AllvsCHF_etacut->SetTitle(";N_{charged}+N_{neutral};CHF");
TH2D *hMT2ChMNeuvsCHF_etacut = new TH2D("hMT2ChMNeuvsCHF_etacut","",40,0,200,20,0,1); hMT2ChMNeuvsCHF_etacut->SetTitle(";N_{charged}-N_{neutral};CHF");
TH2D *hMT2MassvsCHF_etacut = new TH2D("hMT2MassvsCHF_etacut","",60,0,300,20,0,1); hMT2MassvsCHF_etacut->SetTitle(";Mass [GeV];CHF");
TH2D *hMT2METDPhivsCHF_etacut = new TH2D("hMT2METDPhivsCHF_etacut","",16,0,3.2,20,0,1); hMT2METDPhivsCHF_etacut->SetTitle(";MET-#Delta#phi;CHF");
TH2D *hMT2METDPhivsMass_etacut = new TH2D("hMT2METDPhivsMass_etacut","",16,0,3.2,60,0,300); hMT2METDPhivsMass_etacut->SetTitle(";MET-#Delta#phi;Mass [GeV]");
gROOT->cd();

TString precuts = "((misc.HT>750&&misc.MET>30)||(misc.HT>450&&misc.HT<750&&misc.MET>200))&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0";
TString preMT2cuts = "misc.MT2>100&&((misc.HT>750&&misc.MET>30)||(misc.HT>450&&misc.HT<750&&misc.MET>200))&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0";

gROOT->cd();

TString variable; TString selection;

	    MT2tree* fMT2tree = new MT2tree();
	    data->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  data->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;

   	    data->Draw(">>selList", precuts);
        
	    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
	    data->SetEventList(myEvtList);
	    int counter=0;
	    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
//	    if(myEvtList->GetSize()==0) continue;
        while(myEvtList->GetEntry(counter++) !=-1){	

		int jentry = myEvtList->GetEntry(counter-1);
            
		nb =  data->GetEntry(jentry);   nbytes += nb;
		data->SetBranchAddress("MT2tree", &fMT2tree);
		
		if (counter  <= 1  )  cout << "Starting Loop " << counter << endl;
		if (counter % 100000 == 0  )  cout << "+++ Proccessing event " << counter << endl;

		if(fMT2tree->misc.CaloMETRaw      >0) hPFMEToverCaloMET->Fill(fMT2tree->misc.MET/fMT2tree->misc.CaloMETRaw);
		if(fMT2tree->misc.CaloMETMuJesCorr>0) hPFMEToverCaloMETMuCorr->Fill(fMT2tree->misc.MET/fMT2tree->misc.CaloMETMuJesCorr);
		if(fMT2tree->misc.CaloMETRaw      >0) hPFMETvsCaloMET->Fill(fMT2tree->misc.MET,fMT2tree->misc.CaloMETRaw);
		if(fMT2tree->misc.CaloMETMuJesCorr>0) hPFMETvsCaloMETMuCorr->Fill(fMT2tree->misc.MET,fMT2tree->misc.CaloMETMuJesCorr);
		if(fMT2tree->misc.MT2>200){
			if(fMT2tree->misc.CaloMETRaw      >0) hMT2PFMEToverCaloMET->Fill(fMT2tree->misc.MET/fMT2tree->misc.CaloMETRaw);
			if(fMT2tree->misc.CaloMETMuJesCorr>0) hMT2PFMEToverCaloMETMuCorr->Fill(fMT2tree->misc.MET/fMT2tree->misc.CaloMETMuJesCorr);
			if(fMT2tree->misc.CaloMETRaw      >0) hMT2PFMETvsCaloMET->Fill(fMT2tree->misc.MET,fMT2tree->misc.CaloMETRaw);
			if(fMT2tree->misc.CaloMETMuJesCorr>0) hMT2PFMETvsCaloMETMuCorr->Fill(fMT2tree->misc.MET,fMT2tree->misc.CaloMETMuJesCorr);
		}


	for(int i =0; i<fMT2tree->NJets; ++i){
		if(fMT2tree->jet[i].isPFIDLoose==false) continue;
		if(fMT2tree->jet[i].lv.Pt()<40.       ) continue;

		hChvsEta->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].lv.Eta());
		hNeuvsEta->Fill(fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].lv.Eta());
		hAllvsEta->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].lv.Eta());
		hChMNeuvsEta->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].lv.Eta());
		hMassvsEta->Fill(fMT2tree->jet[i].lv.M(),fMT2tree->jet[i].lv.Eta());
		hMETDPhivsEta->Fill(fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]),fMT2tree->jet[i].lv.Eta());
		hCHFvsEta->Fill(fMT2tree->jet[i].ChHadFrac,fMT2tree->jet[i].lv.Eta());
		if(fMT2tree->misc.MT2>200.){
			hMT2ChvsEta->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].lv.Eta());
			hMT2NeuvsEta->Fill(fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].lv.Eta());
			hMT2AllvsEta->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].lv.Eta());
			hMT2ChMNeuvsEta->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].lv.Eta());
			hMT2MassvsEta->Fill(fMT2tree->jet[i].lv.M(),fMT2tree->jet[i].lv.Eta());
			hMT2METDPhivsEta->Fill(fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]),fMT2tree->jet[i].lv.Eta());
			hMT2CHFvsEta->Fill(fMT2tree->jet[i].ChHadFrac,fMT2tree->jet[i].lv.Eta());
		}
		if(fabs(fMT2tree->jet[i].lv.Eta())>2.4) continue;

		hChvsNeu->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].NeuMult);
		hChvsAll->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].NConstituents);
		hNeuvsAll->Fill(fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].NConstituents);
		hChMNeuvsAll->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].NConstituents);
		hChvsMass->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].lv.M());
		hAllvsMass->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].lv.M());
		hChMNeuvsMass->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].lv.M());
		hChvsMETDPhi->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]));
		hAllvsMETDPhi->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]));
		hChMNeuvsMETDPhi->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]));
		hChvsCHF->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].ChHadFrac);
		hNeuvsCHF->Fill(fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].ChHadFrac);
		hAllvsCHF->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].ChHadFrac);
		hChMNeuvsCHF->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].ChHadFrac);
		hMassvsCHF->Fill(fMT2tree->jet[i].lv.M(),fMT2tree->jet[i].ChHadFrac);
		hMETDPhivsCHF->Fill(fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]),fMT2tree->jet[i].ChHadFrac);
		hMETDPhivsMass->Fill(fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]),fMT2tree->jet[i].lv.M());
		if(fMT2tree->misc.MT2>200.){
			hMT2ChvsNeu->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].NeuMult);
			hMT2ChvsAll->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].NConstituents);
			hMT2NeuvsAll->Fill(fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].NConstituents);
			hMT2ChMNeuvsAll->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].NConstituents);
			hMT2ChvsMass->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].lv.M());
			hMT2AllvsMass->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].lv.M());
			hMT2ChMNeuvsMass->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].lv.M());
			hMT2ChvsMETDPhi->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]));
			hMT2AllvsMETDPhi->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]));
			hMT2ChMNeuvsMETDPhi->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]));
			hMT2ChvsCHF->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].ChHadFrac);
			hMT2NeuvsCHF->Fill(fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].ChHadFrac);
			hMT2AllvsCHF->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].ChHadFrac);
			hMT2ChMNeuvsCHF->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].ChHadFrac);
			hMT2MassvsCHF->Fill(fMT2tree->jet[i].lv.M(),fMT2tree->jet[i].ChHadFrac);
			hMT2METDPhivsCHF->Fill(fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]),fMT2tree->jet[i].ChHadFrac);
			hMT2METDPhivsMass->Fill(fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]),fMT2tree->jet[i].lv.M());
		}
		if(fabs(fMT2tree->jet[i].lv.Eta())>1.9) continue;
		if(fabs(fMT2tree->jet[i].lv.Eta())<0.9) continue;
		hChvsNeu_etacut->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].NeuMult);
		hChvsAll_etacut->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].NConstituents);
		hNeuvsAll_etacut->Fill(fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].NConstituents);
		hChMNeuvsAll_etacut->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].NConstituents);
		hChvsMass_etacut->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].lv.M());
		hAllvsMass_etacut->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].lv.M());
		hChMNeuvsMass_etacut->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].lv.M());
		hChvsMETDPhi_etacut->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]));
		hAllvsMETDPhi_etacut->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]));
		hChMNeuvsMETDPhi_etacut->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]));
		hChvsCHF_etacut->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].ChHadFrac);
		hNeuvsCHF_etacut->Fill(fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].ChHadFrac);
		hAllvsCHF_etacut->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].ChHadFrac);
		hChMNeuvsCHF_etacut->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].ChHadFrac);
		hMassvsCHF_etacut->Fill(fMT2tree->jet[i].lv.M(),fMT2tree->jet[i].ChHadFrac);
		hMETDPhivsCHF_etacut->Fill(fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]),fMT2tree->jet[i].ChHadFrac);
		hMETDPhivsMass_etacut->Fill(fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]),fMT2tree->jet[i].lv.M());
		if(fMT2tree->misc.MT2>200.){
			hMT2ChvsNeu_etacut->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].NeuMult);
			hMT2ChvsAll_etacut->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].NConstituents);
			hMT2NeuvsAll_etacut->Fill(fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].NConstituents);
			hMT2ChMNeuvsAll_etacut->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].NConstituents);
			hMT2ChvsMass_etacut->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].lv.M());
			hMT2AllvsMass_etacut->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].lv.M());
			hMT2ChMNeuvsMass_etacut->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].lv.M());
			hMT2ChvsMETDPhi_etacut->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]));
			hMT2AllvsMETDPhi_etacut->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]));
			hMT2ChMNeuvsMETDPhi_etacut->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]));
			hMT2ChvsCHF_etacut->Fill(fMT2tree->jet[i].ChMult,fMT2tree->jet[i].ChHadFrac);
			hMT2NeuvsCHF_etacut->Fill(fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].ChHadFrac);
			hMT2AllvsCHF_etacut->Fill(fMT2tree->jet[i].NConstituents,fMT2tree->jet[i].ChHadFrac);
			hMT2ChMNeuvsCHF_etacut->Fill(fMT2tree->jet[i].ChMult-fMT2tree->jet[i].NeuMult,fMT2tree->jet[i].ChHadFrac);
			hMT2MassvsCHF_etacut->Fill(fMT2tree->jet[i].lv.M(),fMT2tree->jet[i].ChHadFrac);
			hMT2METDPhivsCHF_etacut->Fill(fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]),fMT2tree->jet[i].ChHadFrac);
			hMT2METDPhivsMass_etacut->Fill(fMT2tree->jet[i].lv.DeltaPhi(fMT2tree->pfmet[0]),fMT2tree->jet[i].lv.M());
		}
	}
}//while

gROOT->cd();
TFile *f = new TFile("TOBTECtesthistos.root","RECREATE");
f->cd();

hPFMEToverCaloMET->Write();
hPFMEToverCaloMETMuCorr->Write();
hPFMETvsCaloMET->Write();
hPFMETvsCaloMETMuCorr->Write();
hMT2PFMEToverCaloMET->Write();
hMT2PFMEToverCaloMETMuCorr->Write();
hMT2PFMETvsCaloMET->Write();
hMT2PFMETvsCaloMETMuCorr->Write();

hChvsNeu->Write();
hChvsAll->Write();
hNeuvsAll->Write();
hChMNeuvsAll->Write();
hChvsMass->Write();
hAllvsMass->Write();
hChMNeuvsMass->Write();
hChvsMETDPhi->Write();
hAllvsMETDPhi->Write();
hChMNeuvsMETDPhi->Write();
hChvsCHF->Write();
hNeuvsCHF->Write();
hAllvsCHF->Write();
hChMNeuvsCHF->Write();
hMassvsCHF->Write();
hMETDPhivsCHF->Write();
hMETDPhivsMass->Write();

hChvsEta->Write();
hNeuvsEta->Write();
hAllvsEta->Write();
hChMNeuvsEta->Write();
hMassvsEta->Write();
hMETDPhivsEta->Write();
hCHFvsEta->Write();

hMT2ChvsNeu->Write();
hMT2ChvsAll->Write();
hMT2NeuvsAll->Write();
hMT2ChMNeuvsAll->Write();
hMT2ChvsMass->Write();
hMT2AllvsMass->Write();
hMT2ChMNeuvsMass->Write();
hMT2ChvsMETDPhi->Write();
hMT2AllvsMETDPhi->Write();
hMT2ChMNeuvsMETDPhi->Write();
hMT2ChvsCHF->Write();
hMT2NeuvsCHF->Write();
hMT2AllvsCHF->Write();
hMT2ChMNeuvsCHF->Write();
hMT2MassvsCHF->Write();
hMT2METDPhivsCHF->Write();
hMT2METDPhivsMass->Write();

hMT2ChvsEta->Write();
hMT2NeuvsEta->Write();
hMT2AllvsEta->Write();
hMT2ChMNeuvsEta->Write();
hMT2MassvsEta->Write();
hMT2METDPhivsEta->Write();
hMT2CHFvsEta->Write();

hChvsNeu_etacut->Write();
hChvsAll_etacut->Write();
hNeuvsAll_etacut->Write();
hChMNeuvsAll_etacut->Write();
hChvsMass_etacut->Write();
hAllvsMass_etacut->Write();
hChMNeuvsMass_etacut->Write();
hChvsMETDPhi_etacut->Write();
hAllvsMETDPhi_etacut->Write();
hChMNeuvsMETDPhi_etacut->Write();
hChvsCHF_etacut->Write();
hNeuvsCHF_etacut->Write();
hAllvsCHF_etacut->Write();
hChMNeuvsCHF_etacut->Write();
hMassvsCHF_etacut->Write();
hMETDPhivsCHF_etacut->Write();
hMETDPhivsMass_etacut->Write();
hMT2ChvsNeu_etacut->Write();
hMT2ChvsAll_etacut->Write();
hMT2NeuvsAll_etacut->Write();
hMT2ChMNeuvsAll_etacut->Write();
hMT2ChvsMass_etacut->Write();
hMT2AllvsMass_etacut->Write();
hMT2ChMNeuvsMass_etacut->Write();
hMT2ChvsMETDPhi_etacut->Write();
hMT2AllvsMETDPhi_etacut->Write();
hMT2ChMNeuvsMETDPhi_etacut->Write();
hMT2ChvsCHF_etacut->Write();
hMT2NeuvsCHF_etacut->Write();
hMT2AllvsCHF_etacut->Write();
hMT2ChMNeuvsCHF_etacut->Write();
hMT2MassvsCHF_etacut->Write();
hMT2METDPhivsCHF_etacut->Write();
hMT2METDPhivsMass_etacut->Write();

   TCanvas *c1 = new TCanvas("c1", "c1",485,220,700,504);
//   c1->Range(82.71719,-0.4425771,532.9945,2.212885);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
  // c1->SetLeftMargin(0.1494253);
  // c1->SetRightMargin(0.07327586);
  // c1->SetTopMargin(0.08016878);
  // c1->SetBottomMargin(0.1666667);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);
gPad->SetLogy();

string base = "TOBTECPlots/";
string totname = base + hPFMEToverCaloMET->GetName() + ".eps"; hPFMEToverCaloMET->Draw(); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hPFMEToverCaloMETMuCorr->GetName() + ".eps"; hPFMEToverCaloMETMuCorr->Draw(); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2PFMEToverCaloMET->GetName() + ".eps"; hMT2PFMEToverCaloMET->Draw(); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2PFMEToverCaloMETMuCorr->GetName() + ".eps"; hMT2PFMEToverCaloMETMuCorr->Draw(); c1->SaveAs(totname.c_str()); c1->Clear();

c1->cd(); totname = base + hPFMEToverCaloMET->GetName() + ".pdf"; hPFMEToverCaloMET->Draw(); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hPFMEToverCaloMETMuCorr->GetName() + ".pdf"; hPFMEToverCaloMETMuCorr->Draw(); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2PFMEToverCaloMET->GetName() + ".pdf"; hMT2PFMEToverCaloMET->Draw(); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2PFMEToverCaloMETMuCorr->GetName() + ".pdf"; hMT2PFMEToverCaloMETMuCorr->Draw(); c1->SaveAs(totname.c_str()); c1->Clear();

gPad->SetLogy(0);
gPad->SetLogz();

c1->cd(); totname = base + hPFMETvsCaloMET->GetName() + ".eps"; hPFMETvsCaloMET->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hPFMETvsCaloMETMuCorr->GetName() + ".eps"; hPFMETvsCaloMETMuCorr->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2PFMETvsCaloMET->GetName() + ".eps"; hMT2PFMETvsCaloMET->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2PFMETvsCaloMETMuCorr->GetName() + ".eps"; hMT2PFMETvsCaloMETMuCorr->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();

c1->cd(); totname = base + hChvsNeu->GetName() + ".eps"; hChvsNeu->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsAll->GetName() + ".eps"; hChvsAll->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hNeuvsAll->GetName() + ".eps"; hNeuvsAll->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsAll->GetName() + ".eps"; hChMNeuvsAll->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsMass->GetName() + ".eps"; hChvsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsMass->GetName() + ".eps"; hAllvsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsMass->GetName() + ".eps"; hChMNeuvsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsMETDPhi->GetName() + ".eps"; hChvsMETDPhi->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsMETDPhi->GetName() + ".eps"; hAllvsMETDPhi->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsMETDPhi->GetName() + ".eps"; hChMNeuvsMETDPhi->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsCHF->GetName() + ".eps"; hChvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hNeuvsCHF->GetName() + ".eps"; hNeuvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsCHF->GetName() + ".eps"; hAllvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsCHF->GetName() + ".eps"; hChMNeuvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMassvsCHF->GetName() + ".eps"; hMassvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMETDPhivsCHF->GetName() + ".eps"; hMETDPhivsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMETDPhivsMass->GetName() + ".eps"; hMETDPhivsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsNeu->GetName() + ".eps"; hMT2ChvsNeu->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsAll->GetName() + ".eps"; hMT2ChvsAll->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2NeuvsAll->GetName() + ".eps"; hMT2NeuvsAll->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsAll->GetName() + ".eps"; hMT2ChMNeuvsAll->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsMass->GetName() + ".eps"; hMT2ChvsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsMass->GetName() + ".eps"; hMT2AllvsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsMass->GetName() + ".eps"; hMT2ChMNeuvsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsMETDPhi->GetName() + ".eps"; hMT2ChvsMETDPhi->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsMETDPhi->GetName() + ".eps"; hMT2AllvsMETDPhi->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsMETDPhi->GetName() + ".eps"; hMT2ChMNeuvsMETDPhi->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsCHF->GetName() + ".eps"; hMT2ChvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2NeuvsCHF->GetName() + ".eps"; hMT2NeuvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsCHF->GetName() + ".eps"; hMT2AllvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsCHF->GetName() + ".eps"; hMT2ChMNeuvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2MassvsCHF->GetName() + ".eps"; hMT2MassvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2METDPhivsCHF->GetName() + ".eps"; hMT2METDPhivsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2METDPhivsMass->GetName() + ".eps"; hMT2METDPhivsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();

c1->cd(); totname = base + hChvsNeu_etacut->GetName() + ".eps"; hChvsNeu_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsAll_etacut->GetName() + ".eps"; hChvsAll_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hNeuvsAll_etacut->GetName() + ".eps"; hNeuvsAll_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsAll_etacut->GetName() + ".eps"; hChMNeuvsAll_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsMass_etacut->GetName() + ".eps"; hChvsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsMass_etacut->GetName() + ".eps"; hAllvsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsMass_etacut->GetName() + ".eps"; hChMNeuvsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsMETDPhi_etacut->GetName() + ".eps"; hChvsMETDPhi_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsMETDPhi_etacut->GetName() + ".eps"; hAllvsMETDPhi_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsMETDPhi_etacut->GetName() + ".eps"; hChMNeuvsMETDPhi_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsCHF_etacut->GetName() + ".eps"; hChvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hNeuvsCHF_etacut->GetName() + ".eps"; hNeuvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsCHF_etacut->GetName() + ".eps"; hAllvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsCHF_etacut->GetName() + ".eps"; hChMNeuvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMassvsCHF_etacut->GetName() + ".eps"; hMassvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMETDPhivsCHF_etacut->GetName() + ".eps"; hMETDPhivsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMETDPhivsMass_etacut->GetName() + ".eps"; hMETDPhivsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsNeu_etacut->GetName() + ".eps"; hMT2ChvsNeu_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsAll_etacut->GetName() + ".eps"; hMT2ChvsAll_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2NeuvsAll_etacut->GetName() + ".eps"; hMT2NeuvsAll_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsAll_etacut->GetName() + ".eps"; hMT2ChMNeuvsAll_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsMass_etacut->GetName() + ".eps"; hMT2ChvsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsMass_etacut->GetName() + ".eps"; hMT2AllvsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsMass_etacut->GetName() + ".eps"; hMT2ChMNeuvsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsMETDPhi_etacut->GetName() + ".eps"; hMT2ChvsMETDPhi_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsMETDPhi_etacut->GetName() + ".eps"; hMT2AllvsMETDPhi_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsMETDPhi_etacut->GetName() + ".eps"; hMT2ChMNeuvsMETDPhi_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsCHF_etacut->GetName() + ".eps"; hMT2ChvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2NeuvsCHF_etacut->GetName() + ".eps"; hMT2NeuvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsCHF_etacut->GetName() + ".eps"; hMT2AllvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsCHF_etacut->GetName() + ".eps"; hMT2ChMNeuvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2MassvsCHF_etacut->GetName() + ".eps"; hMT2MassvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2METDPhivsCHF_etacut->GetName() + ".eps"; hMT2METDPhivsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2METDPhivsMass_etacut->GetName() + ".eps"; hMT2METDPhivsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();

c1->cd(); totname = base + hChvsEta->GetName() + ".eps"; hChvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hNeuvsEta->GetName() + ".eps"; hNeuvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsEta->GetName() + ".eps"; hAllvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsEta->GetName() + ".eps"; hChMNeuvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMassvsEta->GetName() + ".eps"; hMassvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMETDPhivsEta->GetName() + ".eps"; hMETDPhivsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hCHFvsEta->GetName() + ".eps"; hCHFvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsEta->GetName() + ".eps"; hMT2ChvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2NeuvsEta->GetName() + ".eps"; hMT2NeuvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsEta->GetName() + ".eps"; hMT2AllvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsEta->GetName() + ".eps"; hMT2ChMNeuvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2MassvsEta->GetName() + ".eps"; hMT2MassvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2METDPhivsEta->GetName() + ".eps"; hMT2METDPhivsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2CHFvsEta->GetName() + ".eps"; hMT2CHFvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();


c1->cd(); totname = base + hPFMETvsCaloMET->GetName() + ".pdf"; hPFMETvsCaloMET->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hPFMETvsCaloMETMuCorr->GetName() + ".pdf"; hPFMETvsCaloMETMuCorr->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2PFMETvsCaloMET->GetName() + ".pdf"; hMT2PFMETvsCaloMET->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2PFMETvsCaloMETMuCorr->GetName() + ".pdf"; hMT2PFMETvsCaloMETMuCorr->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsNeu->GetName() + ".pdf"; hChvsNeu->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsAll->GetName() + ".pdf"; hChvsAll->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hNeuvsAll->GetName() + ".pdf"; hNeuvsAll->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsAll->GetName() + ".pdf"; hChMNeuvsAll->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsMass->GetName() + ".pdf"; hChvsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsMass->GetName() + ".pdf"; hAllvsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsMass->GetName() + ".pdf"; hChMNeuvsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsMETDPhi->GetName() + ".pdf"; hChvsMETDPhi->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsMETDPhi->GetName() + ".pdf"; hAllvsMETDPhi->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsMETDPhi->GetName() + ".pdf"; hChMNeuvsMETDPhi->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsCHF->GetName() + ".pdf"; hChvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hNeuvsCHF->GetName() + ".pdf"; hNeuvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsCHF->GetName() + ".pdf"; hAllvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsCHF->GetName() + ".pdf"; hChMNeuvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMassvsCHF->GetName() + ".pdf"; hMassvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMETDPhivsCHF->GetName() + ".pdf"; hMETDPhivsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMETDPhivsMass->GetName() + ".pdf"; hMETDPhivsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsNeu->GetName() + ".pdf"; hMT2ChvsNeu->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsAll->GetName() + ".pdf"; hMT2ChvsAll->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2NeuvsAll->GetName() + ".pdf"; hMT2NeuvsAll->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsAll->GetName() + ".pdf"; hMT2ChMNeuvsAll->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsMass->GetName() + ".pdf"; hMT2ChvsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsMass->GetName() + ".pdf"; hMT2AllvsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsMass->GetName() + ".pdf"; hMT2ChMNeuvsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsMETDPhi->GetName() + ".pdf"; hMT2ChvsMETDPhi->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsMETDPhi->GetName() + ".pdf"; hMT2AllvsMETDPhi->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsMETDPhi->GetName() + ".pdf"; hMT2ChMNeuvsMETDPhi->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsCHF->GetName() + ".pdf"; hMT2ChvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2NeuvsCHF->GetName() + ".pdf"; hMT2NeuvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsCHF->GetName() + ".pdf"; hMT2AllvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsCHF->GetName() + ".pdf"; hMT2ChMNeuvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2MassvsCHF->GetName() + ".pdf"; hMT2MassvsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2METDPhivsCHF->GetName() + ".pdf"; hMT2METDPhivsCHF->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2METDPhivsMass->GetName() + ".pdf"; hMT2METDPhivsMass->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsNeu_etacut->GetName() + ".pdf"; hChvsNeu_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsAll_etacut->GetName() + ".pdf"; hChvsAll_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hNeuvsAll_etacut->GetName() + ".pdf"; hNeuvsAll_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsAll_etacut->GetName() + ".pdf"; hChMNeuvsAll_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsMass_etacut->GetName() + ".pdf"; hChvsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsMass_etacut->GetName() + ".pdf"; hAllvsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsMass_etacut->GetName() + ".pdf"; hChMNeuvsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsMETDPhi_etacut->GetName() + ".pdf"; hChvsMETDPhi_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsMETDPhi_etacut->GetName() + ".pdf"; hAllvsMETDPhi_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsMETDPhi_etacut->GetName() + ".pdf"; hChMNeuvsMETDPhi_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsCHF_etacut->GetName() + ".pdf"; hChvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hNeuvsCHF_etacut->GetName() + ".pdf"; hNeuvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsCHF_etacut->GetName() + ".pdf"; hAllvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsCHF_etacut->GetName() + ".pdf"; hChMNeuvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMassvsCHF_etacut->GetName() + ".pdf"; hMassvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMETDPhivsCHF_etacut->GetName() + ".pdf"; hMETDPhivsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMETDPhivsMass_etacut->GetName() + ".pdf"; hMETDPhivsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsNeu_etacut->GetName() + ".pdf"; hMT2ChvsNeu_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsAll_etacut->GetName() + ".pdf"; hMT2ChvsAll_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2NeuvsAll_etacut->GetName() + ".pdf"; hMT2NeuvsAll_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsAll_etacut->GetName() + ".pdf"; hMT2ChMNeuvsAll_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsMass_etacut->GetName() + ".pdf"; hMT2ChvsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsMass_etacut->GetName() + ".pdf"; hMT2AllvsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsMass_etacut->GetName() + ".pdf"; hMT2ChMNeuvsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsMETDPhi_etacut->GetName() + ".pdf"; hMT2ChvsMETDPhi_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsMETDPhi_etacut->GetName() + ".pdf"; hMT2AllvsMETDPhi_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsMETDPhi_etacut->GetName() + ".pdf"; hMT2ChMNeuvsMETDPhi_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsCHF_etacut->GetName() + ".pdf"; hMT2ChvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2NeuvsCHF_etacut->GetName() + ".pdf"; hMT2NeuvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsCHF_etacut->GetName() + ".pdf"; hMT2AllvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsCHF_etacut->GetName() + ".pdf"; hMT2ChMNeuvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2MassvsCHF_etacut->GetName() + ".pdf"; hMT2MassvsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2METDPhivsCHF_etacut->GetName() + ".pdf"; hMT2METDPhivsCHF_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2METDPhivsMass_etacut->GetName() + ".pdf"; hMT2METDPhivsMass_etacut->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChvsEta->GetName() + ".pdf"; hChvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hNeuvsEta->GetName() + ".pdf"; hNeuvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hAllvsEta->GetName() + ".pdf"; hAllvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hChMNeuvsEta->GetName() + ".pdf"; hChMNeuvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMassvsEta->GetName() + ".pdf"; hMassvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMETDPhivsEta->GetName() + ".pdf"; hMETDPhivsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hCHFvsEta->GetName() + ".pdf"; hCHFvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChvsEta->GetName() + ".pdf"; hMT2ChvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2NeuvsEta->GetName() + ".pdf"; hMT2NeuvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2AllvsEta->GetName() + ".pdf"; hMT2AllvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2ChMNeuvsEta->GetName() + ".pdf"; hMT2ChMNeuvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2MassvsEta->GetName() + ".pdf"; hMT2MassvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2METDPhivsEta->GetName() + ".pdf"; hMT2METDPhivsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();
c1->cd(); totname = base + hMT2CHFvsEta->GetName() + ".pdf"; hMT2CHFvsEta->Draw("colz"); c1->SaveAs(totname.c_str()); c1->Clear();

}