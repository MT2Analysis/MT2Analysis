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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path of MT2tree.hh

//use via root -l -b -q PassingEventStudy2.C++

using namespace std;

//plots TOBTEC noise filter efficiency (actually passing efficiency = 1-filter efficiency) via passing events / all events
//for several variables
//and optionally fit this efficiency to get a number
//also 'performs reweighting', i.e.:
//select one jet most likely cause of effect (see more details in AN2013/215)
//and saves passing efficiency w.r.t. this jet
//and thus creates map that will be used for the reweighting
void PassingEventStudy2(){

bool fit = false;//perform fit, default  = false
float scale = 1.;//jet pt scaled, default = 1, use for roughly seeing influence of JES

gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetPaintTextFormat(".2f");

TChain *c = new TChain("MassTree");
//datafiles
int cnum1 = c->Add("~mangano/filterTOBTEC/codeHJ/CMSSW_5_3_11/myTest/MT2Code/output/bloc1/*.root");
int cnum2 = c->Add("~mangano/filterTOBTEC/codeHJ/CMSSW_5_3_11/myTest/MT2Code/output/bloc2/*.root");
int cnum3 = c->Add("~mangano/filterTOBTEC/codeHJ/CMSSW_5_3_11/myTest/MT2Code/output/bloc3/*.root");
int cnum4 = c->Add("~mangano/filterTOBTEC/codeHJ/CMSSW_5_3_11/myTest/MT2Code/output/bloc4/*.root");
int cnum5 = c->Add("~/MT2Analysis/MT2trees/MT2_V02-03-02/20130914_8TeV_1g_removed/lowHT/SinglePhoton-Run2012C-PromptReco-v2-2_Filter.root");
cout << cnum3 << " " << cnum4 << " " << cnum5 << endl;
cout << "files in chain: " << cnum1 + cnum2 << endl;

//signal selection, including TOBTEC filter (all cuts, all cuts but MinDPhi(jets,MET), leptons, and/or photons)
TString cutallnomindphi = "((misc.MT2>200&&misc.HT>=450&&misc.HT<750&&misc.MET>200&&((trigger.HLT_PFMET150_v2==1||trigger.HLT_PFMET150_v3==1||trigger.HLT_PFMET150_v4==1||trigger.HLT_PFMET150_v5==1||trigger.HLT_PFMET150_v6==1||trigger.HLT_PFMET150_v7==1)||(trigger.HLT_PFHT350_PFMET100_v3==1||trigger.HLT_PFHT350_PFMET100_v4==1||trigger.HLT_PFHT350_PFMET100_v5==1||trigger.HLT_PFHT350_PFMET100_v6==1||trigger.HLT_PFHT350_PFMET100_v7==1||trigger.HLT_PFNoPUHT350_PFMET100_v1==1||trigger.HLT_PFNoPUHT350_PFMET100_v3==1||trigger.HLT_PFNoPUHT350_PFMET100_v4==1)))||(misc.MT2>100&misc.HT>=750&&misc.MET>30&&(trigger.HLT_PFHT650_v5==1||trigger.HLT_PFHT650_v6==1||trigger.HLT_PFHT650_v7==1||trigger.HLT_PFHT650_v8==1||trigger.HLT_PFHT650_v9==1||trigger.HLT_PFNoPUHT650_v1==1||trigger.HLT_PFNoPUHT650_v3==1||trigger.HLT_PFNoPUHT650_v4==1)))&&NMuons==0 && NEles==0 && NTausIDLoose3Hits==0&&misc.Jet0Pass ==1&&misc.Jet1Pass ==1&&misc.SecondJPt  > 100&&misc.PassJet40ID ==1&&misc.Vectorsumpt < 70&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&misc.MET/misc.CaloMETRaw<=2.&&ExtraBeamHaloFilter==0&&TOBTECTagger>=0";
TString cutallnomindphinolepveto = "((misc.MT2>200&&misc.HT>=450&&misc.HT<750&&misc.MET>200&&((trigger.HLT_PFMET150_v2==1||trigger.HLT_PFMET150_v3==1||trigger.HLT_PFMET150_v4==1||trigger.HLT_PFMET150_v5==1||trigger.HLT_PFMET150_v6==1||trigger.HLT_PFMET150_v7==1)||(trigger.HLT_PFHT350_PFMET100_v3==1||trigger.HLT_PFHT350_PFMET100_v4==1||trigger.HLT_PFHT350_PFMET100_v5==1||trigger.HLT_PFHT350_PFMET100_v6==1||trigger.HLT_PFHT350_PFMET100_v7==1||trigger.HLT_PFNoPUHT350_PFMET100_v1==1||trigger.HLT_PFNoPUHT350_PFMET100_v3==1||trigger.HLT_PFNoPUHT350_PFMET100_v4==1)))||(misc.MT2>100&misc.HT>=750&&misc.MET>30&&(trigger.HLT_PFHT650_v5==1||trigger.HLT_PFHT650_v6==1||trigger.HLT_PFHT650_v7==1||trigger.HLT_PFHT650_v8==1||trigger.HLT_PFHT650_v9==1||trigger.HLT_PFNoPUHT650_v1==1||trigger.HLT_PFNoPUHT650_v3==1||trigger.HLT_PFNoPUHT650_v4==1)))&&misc.Jet0Pass ==1&&misc.Jet1Pass ==1&&misc.SecondJPt  > 100&&misc.PassJet40ID ==1&&misc.Vectorsumpt < 70&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&misc.MET/misc.CaloMETRaw<=2.&&ExtraBeamHaloFilter==0&&TOBTECTagger>=0";
TString cutallnomindphinolepvetoandphotons = "((misc.MT2>200&&misc.HT>=450&&misc.HT<750&&misc.MET>200&&((trigger.HLT_PFMET150_v2==1||trigger.HLT_PFMET150_v3==1||trigger.HLT_PFMET150_v4==1||trigger.HLT_PFMET150_v5==1||trigger.HLT_PFMET150_v6==1||trigger.HLT_PFMET150_v7==1)||(trigger.HLT_PFHT350_PFMET100_v3==1||trigger.HLT_PFHT350_PFMET100_v4==1||trigger.HLT_PFHT350_PFMET100_v5==1||trigger.HLT_PFHT350_PFMET100_v6==1||trigger.HLT_PFHT350_PFMET100_v7==1||trigger.HLT_PFNoPUHT350_PFMET100_v1==1||trigger.HLT_PFNoPUHT350_PFMET100_v3==1||trigger.HLT_PFNoPUHT350_PFMET100_v4==1)))||(misc.MT2>100&misc.HT>=750&&misc.MET>30&&(trigger.HLT_PFHT650_v5==1||trigger.HLT_PFHT650_v6==1||trigger.HLT_PFHT650_v7==1||trigger.HLT_PFHT650_v8==1||trigger.HLT_PFHT650_v9==1||trigger.HLT_PFNoPUHT650_v1==1||trigger.HLT_PFNoPUHT650_v3==1||trigger.HLT_PFNoPUHT650_v4==1))||trigger.HLT_SinglePhotons)&&misc.Jet0Pass ==1&&misc.Jet1Pass ==1&&misc.SecondJPt  > 100&&misc.PassJet40ID ==1&&misc.Vectorsumpt < 70&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&type1pfmet.Pt()/misc.CaloMETRaw<=2.&&ExtraBeamHaloFilter==0&&TOBTECTagger>=0";

//signal selection, excluding TOBTEC filter (all cuts, all cuts but MinDPhi(jets,MET), leptons, and/or photons)
TString cutall = "((misc.MT2>200&&misc.HT>=450&&misc.HT<750&&misc.MET>200&&((trigger.HLT_PFMET150_v2==1||trigger.HLT_PFMET150_v3==1||trigger.HLT_PFMET150_v4==1||trigger.HLT_PFMET150_v5==1||trigger.HLT_PFMET150_v6==1||trigger.HLT_PFMET150_v7==1)||(trigger.HLT_PFHT350_PFMET100_v3==1||trigger.HLT_PFHT350_PFMET100_v4==1||trigger.HLT_PFHT350_PFMET100_v5==1||trigger.HLT_PFHT350_PFMET100_v6==1||trigger.HLT_PFHT350_PFMET100_v7==1||trigger.HLT_PFNoPUHT350_PFMET100_v1==1||trigger.HLT_PFNoPUHT350_PFMET100_v3==1||trigger.HLT_PFNoPUHT350_PFMET100_v4==1)))||(misc.MT2>100&misc.HT>=750&&misc.MET>30&&(trigger.HLT_PFHT650_v5==1||trigger.HLT_PFHT650_v6==1||trigger.HLT_PFHT650_v7==1||trigger.HLT_PFHT650_v8==1||trigger.HLT_PFHT650_v9==1||trigger.HLT_PFNoPUHT650_v1==1||trigger.HLT_PFNoPUHT650_v3==1||trigger.HLT_PFNoPUHT650_v4==1)))&&NMuons==0 && NEles==0 && NTausIDLoose3Hits==0&&misc.Jet0Pass ==1&&misc.Jet1Pass ==1&&misc.SecondJPt  > 100&&misc.PassJet40ID ==1&&misc.MinMetJetDPhi4Pt40 >0.3&&misc.Vectorsumpt < 70&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&misc.MET/misc.CaloMETRaw<=2.&&ExtraBeamHaloFilter==0&&TOBTECTagger>=0";
TString cutpass = "((misc.MT2>200&&misc.HT>=450&&misc.HT<750&&misc.MET>200&&((trigger.HLT_PFMET150_v2==1||trigger.HLT_PFMET150_v3==1||trigger.HLT_PFMET150_v4==1||trigger.HLT_PFMET150_v5==1||trigger.HLT_PFMET150_v6==1||trigger.HLT_PFMET150_v7==1)||(trigger.HLT_PFHT350_PFMET100_v3==1||trigger.HLT_PFHT350_PFMET100_v4==1||trigger.HLT_PFHT350_PFMET100_v5==1||trigger.HLT_PFHT350_PFMET100_v6==1||trigger.HLT_PFHT350_PFMET100_v7==1||trigger.HLT_PFNoPUHT350_PFMET100_v1==1||trigger.HLT_PFNoPUHT350_PFMET100_v3==1||trigger.HLT_PFNoPUHT350_PFMET100_v4==1)))||(misc.MT2>100&misc.HT>=750&&misc.MET>30&&(trigger.HLT_PFHT650_v5==1||trigger.HLT_PFHT650_v6==1||trigger.HLT_PFHT650_v7==1||trigger.HLT_PFHT650_v8==1||trigger.HLT_PFHT650_v9==1||trigger.HLT_PFNoPUHT650_v1==1||trigger.HLT_PFNoPUHT650_v3==1||trigger.HLT_PFNoPUHT650_v4==1)))&&NMuons==0 && NEles==0 && NTausIDLoose3Hits==0&&misc.Jet0Pass ==1&&misc.Jet1Pass ==1&&misc.SecondJPt  > 100&&misc.PassJet40ID ==1&&misc.MinMetJetDPhi4Pt40 >0.3&&misc.Vectorsumpt < 70&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&misc.MET/misc.CaloMETRaw<=2.&&ExtraBeamHaloFilter==0&&TOBTECTagger>=0&&TOBTECTagger<=8";
TString cutallHT = "((misc.MT2>100&&misc.HT>=750&&misc.MET>30&&(trigger.HLT_PFHT650_v5==1||trigger.HLT_PFHT650_v6==1||trigger.HLT_PFHT650_v7==1||trigger.HLT_PFHT650_v8==1||trigger.HLT_PFHT650_v9==1||trigger.HLT_PFNoPUHT650_v1==1||trigger.HLT_PFNoPUHT650_v3==1||trigger.HLT_PFNoPUHT650_v4==1)))&&NMuons==0 && NEles==0 && NTausIDLoose3Hits==0&&misc.Jet0Pass ==1&&misc.Jet1Pass ==1&&misc.SecondJPt  > 100&&misc.PassJet40ID ==1&&misc.MinMetJetDPhi4Pt40 >0.3&&misc.Vectorsumpt < 70&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&misc.MET/misc.CaloMETRaw<=2.&&ExtraBeamHaloFilter==0&&TOBTECTagger>=0";
TString cutpassHT = "((misc.MT2>100&&misc.HT>=750&&misc.MET>30&&(trigger.HLT_PFHT650_v5==1||trigger.HLT_PFHT650_v6==1||trigger.HLT_PFHT650_v7==1||trigger.HLT_PFHT650_v8==1||trigger.HLT_PFHT650_v9==1||trigger.HLT_PFNoPUHT650_v1==1||trigger.HLT_PFNoPUHT650_v3==1||trigger.HLT_PFNoPUHT650_v4==1)))&&NMuons==0 && NEles==0 && NTausIDLoose3Hits==0&&misc.Jet0Pass ==1&&misc.Jet1Pass ==1&&misc.SecondJPt  > 100&&misc.PassJet40ID ==1&&misc.MinMetJetDPhi4Pt40 >0.3&&misc.Vectorsumpt < 70&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&misc.MET/misc.CaloMETRaw<=2.&&ExtraBeamHaloFilter==0&&TOBTECTagger>=0&&TOBTECTagger<=8";
TString cutallMET = "((misc.MT2>200&&misc.HT>=450&&misc.HT<750&&misc.MET>200&&((trigger.HLT_PFMET150_v2==1||trigger.HLT_PFMET150_v3==1||trigger.HLT_PFMET150_v4==1||trigger.HLT_PFMET150_v5==1||trigger.HLT_PFMET150_v6==1||trigger.HLT_PFMET150_v7==1)||(trigger.HLT_PFHT350_PFMET100_v3==1||trigger.HLT_PFHT350_PFMET100_v4==1||trigger.HLT_PFHT350_PFMET100_v5==1||trigger.HLT_PFHT350_PFMET100_v6==1||trigger.HLT_PFHT350_PFMET100_v7==1||trigger.HLT_PFNoPUHT350_PFMET100_v1==1||trigger.HLT_PFNoPUHT350_PFMET100_v3==1||trigger.HLT_PFNoPUHT350_PFMET100_v4==1))))&&NMuons==0 && NEles==0 && NTausIDLoose3Hits==0&&misc.Jet0Pass ==1&&misc.Jet1Pass ==1&&misc.SecondJPt  > 100&&misc.PassJet40ID ==1&&misc.MinMetJetDPhi4Pt40 >0.3&&misc.Vectorsumpt < 70&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&misc.MET/misc.CaloMETRaw<=2.&&ExtraBeamHaloFilter==0&&TOBTECTagger>=0";
TString cutpassMET = "((misc.MT2>200&&misc.HT>=450&&misc.HT<750&&misc.MET>200&&((trigger.HLT_PFMET150_v2==1||trigger.HLT_PFMET150_v3==1||trigger.HLT_PFMET150_v4==1||trigger.HLT_PFMET150_v5==1||trigger.HLT_PFMET150_v6==1||trigger.HLT_PFMET150_v7==1)||(trigger.HLT_PFHT350_PFMET100_v3==1||trigger.HLT_PFHT350_PFMET100_v4==1||trigger.HLT_PFHT350_PFMET100_v5==1||trigger.HLT_PFHT350_PFMET100_v6==1||trigger.HLT_PFHT350_PFMET100_v7==1||trigger.HLT_PFNoPUHT350_PFMET100_v1==1||trigger.HLT_PFNoPUHT350_PFMET100_v3==1||trigger.HLT_PFNoPUHT350_PFMET100_v4==1))))&&NMuons==0 && NEles==0 && NTausIDLoose3Hits==0&&misc.Jet0Pass ==1&&misc.Jet1Pass ==1&&misc.SecondJPt  > 100&&misc.PassJet40ID ==1&&misc.MinMetJetDPhi4Pt40 >0.3&&misc.Vectorsumpt < 70&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&misc.MET/misc.CaloMETRaw<=2.&&ExtraBeamHaloFilter==0&&TOBTECTagger>=0&&TOBTECTagger<=8";

//all histogram
TH1D *MT2all = new TH1D("MT2all","",17,150,1000);
TH1D *MT2allHT = new TH1D("MT2allHT","",17,150,1000);
TH1D *MT2allMET = new TH1D("MT2allMET","",16,200,1000);
TH1D *MT2pass = new TH1D("MT2pass","",17,150,1000);
TH1D *MT2passHT = new TH1D("MT2passHT","",17,150,1000);
TH1D *MT2passMET = new TH1D("MT2passMET","",16,200,1000);
MT2all->GetXaxis()->SetTitle("M_{T2} [GeV]"); MT2all->GetYaxis()->SetTitle("all events / 50 GeV");
MT2allHT->GetXaxis()->SetTitle("M_{T2} [GeV]"); MT2allHT->GetYaxis()->SetTitle("all events / 50 GeV");
MT2allMET->GetXaxis()->SetTitle("M_{T2} [GeV]"); MT2allMET->GetYaxis()->SetTitle("all events / 50 GeV");
MT2pass->GetXaxis()->SetTitle("M_{T2} [GeV]"); MT2pass->GetYaxis()->SetTitle("passing events / 50 GeV");
MT2passHT->GetXaxis()->SetTitle("M_{T2} [GeV]"); MT2passHT->GetYaxis()->SetTitle("passing events / 50 GeV");
MT2passMET->GetXaxis()->SetTitle("M_{T2} [GeV]"); MT2passMET->GetYaxis()->SetTitle("passing events / 50 GeV");
TH1D *HTall = new TH1D("HTall","",16,450,2050);
TH1D *HTallHT = new TH1D("HTallHT","", 13,750,2050);
TH1D *HTallMET = new TH1D("HTallMET","",6,450,750);
TH1D *HTpass = new TH1D("HTpass","", 16,450,2050);
TH1D *HTpassHT = new TH1D("HTpassHT","", 13,750,2050);
TH1D *HTpassMET = new TH1D("HTpassMET","",6,450,750);
HTall->GetXaxis()->SetTitle("H_{T} [GeV]"); HTall->GetYaxis()->SetTitle("all events / 100 GeV");
HTallHT->GetXaxis()->SetTitle("H_{T} [GeV]"); HTallHT->GetYaxis()->SetTitle("all events / 100 GeV");
HTallMET->GetXaxis()->SetTitle("H_{T} [GeV]"); HTallMET->GetYaxis()->SetTitle("all events / 50 GeV");
HTpass->GetXaxis()->SetTitle("H_{T} [GeV]"); HTpass->GetYaxis()->SetTitle("passing events / 100 GeV");
HTpassHT->GetXaxis()->SetTitle("H_{T} [GeV]"); HTpassHT->GetYaxis()->SetTitle("passing events / 100 GeV");
HTpassMET->GetXaxis()->SetTitle("H_{T} [GeV]"); HTpassMET->GetYaxis()->SetTitle("passing events / 50 GeV");
TH1D *METall = new TH1D("METall","",17,150,1000);
TH1D *METallHT = new TH1D("METallHT","",17,150,1000);
TH1D *METallMET = new TH1D("METallMET","",16,200,1000);
TH1D *METpass = new TH1D("METpass","",17,150,1000);
TH1D *METpassHT = new TH1D("METpassHT","",17,150,1000);
TH1D *METpassMET = new TH1D("METpassMET","",16,200,1000);
METall->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]"); METall->GetYaxis()->SetTitle("all events / 50 GeV");
METallHT->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]"); METallHT->GetYaxis()->SetTitle("all events / 50 GeV");
METallMET->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]"); METallMET->GetYaxis()->SetTitle("all events / 50 GeV");
METpass->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]"); METpass->GetYaxis()->SetTitle("passing events / 50 GeV");
METpassHT->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]"); METpassHT->GetYaxis()->SetTitle("passing events / 50 GeV");
METpassMET->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]"); METpassMET->GetYaxis()->SetTitle("passing events / 50 GeV");
TH1D *NJetsall = new TH1D("NJetsall","",8,2,10);
TH1D *NJetsallHT = new TH1D("NJetsallHT","", 8,2,10);
TH1D *NJetsallMET = new TH1D("NJetsallMET","", 8,2,10);
TH1D *NJetspass = new TH1D("NJetspass","", 8,2,10);
TH1D *NJetspassHT = new TH1D("NJetspassHT","", 8,2,10);
TH1D *NJetspassMET = new TH1D("NJetspassMET","", 8,2,10);
NJetsall->GetXaxis()->SetTitle("NJets"); NJetsall->GetYaxis()->SetTitle("all events / 1");
NJetsallHT->GetXaxis()->SetTitle("NJets"); NJetsallHT->GetYaxis()->SetTitle("all events / 1");
NJetsallMET->GetXaxis()->SetTitle("NJets"); NJetsallMET->GetYaxis()->SetTitle("all events / 1");
NJetspass->GetXaxis()->SetTitle("NJets"); NJetspass->GetYaxis()->SetTitle("passing events / 1");
NJetspassHT->GetXaxis()->SetTitle("NJets"); NJetspassHT->GetYaxis()->SetTitle("passing events / 1");
NJetspassMET->GetXaxis()->SetTitle("NJets"); NJetspassMET->GetYaxis()->SetTitle("passing events / 1");
TH1D *NBJetsall = new TH1D("NBJetsall","",5,0,5);
TH1D *NBJetsallHT = new TH1D("NBJetsallHT","", 5,0,5);
TH1D *NBJetsallMET = new TH1D("NBJetsallMET","", 5,0,5);
TH1D *NBJetspass = new TH1D("NBJetspass","", 5,0,5);
TH1D *NBJetspassHT = new TH1D("NBJetspassHT","", 5,0,5);
TH1D *NBJetspassMET = new TH1D("NBJetspassMET","", 5,0,5);
NBJetsall->GetXaxis()->SetTitle("NBJets"); NBJetsall->GetYaxis()->SetTitle("all events / 1");
NBJetsallHT->GetXaxis()->SetTitle("NBJets"); NBJetsallHT->GetYaxis()->SetTitle("all events / 1");
NBJetsallMET->GetXaxis()->SetTitle("NBJets"); NBJetsallMET->GetYaxis()->SetTitle("all events / 1");
NBJetspass->GetXaxis()->SetTitle("NBJets"); NBJetspass->GetYaxis()->SetTitle("passing events / 1");
NBJetspassHT->GetXaxis()->SetTitle("NBJets"); NBJetspassHT->GetYaxis()->SetTitle("passing events / 1");
NBJetspassMET->GetXaxis()->SetTitle("NBJets"); NBJetspassMET->GetYaxis()->SetTitle("passing events / 1");
TH1D *Jet1Etaall = new TH1D("Jet1Etaall","",12,0.,2.4);
TH1D *Jet1EtaallHT = new TH1D("Jet1EtaallHT","", 12,0.,2.4);
TH1D *Jet1EtaallMET = new TH1D("Jet1EtaallMET","", 12,0.,2.4);
TH1D *Jet1Etapass = new TH1D("Jet1Etapass","", 12,0.,2.4);
TH1D *Jet1EtapassHT = new TH1D("Jet1EtapassHT","", 12,0.,2.4);
TH1D *Jet1EtapassMET = new TH1D("Jet1EtapassMET","", 12,0.,2.4);
Jet1Etaall->GetXaxis()->SetTitle("jet 1 #eta"); Jet1Etaall->GetYaxis()->SetTitle("all events / 0.2");
Jet1EtaallHT->GetXaxis()->SetTitle("jet 1 #eta"); Jet1EtaallHT->GetYaxis()->SetTitle("all events / 0.2");
Jet1EtaallMET->GetXaxis()->SetTitle("jet 1 #eta"); Jet1EtaallMET->GetYaxis()->SetTitle("all events / 0.2");
Jet1Etapass->GetXaxis()->SetTitle("jet 1 #eta"); Jet1Etapass->GetYaxis()->SetTitle("passing events / 0.2");
Jet1EtapassHT->GetXaxis()->SetTitle("jet 1 #eta"); Jet1EtapassHT->GetYaxis()->SetTitle("passing events / 0.2");
Jet1EtapassMET->GetXaxis()->SetTitle("jet 1 #eta"); Jet1EtapassMET->GetYaxis()->SetTitle("passing events / 0.2");
TH1D *Jet2Etaall = new TH1D("Jet2Etaall","",12,0.,2.4);
TH1D *Jet2EtaallHT = new TH1D("Jet2EtaallHT","", 12,0.,2.4);
TH1D *Jet2EtaallMET = new TH1D("Jet2EtaallMET","", 12,0.,2.4);
TH1D *Jet2Etapass = new TH1D("Jet2Etapass","", 12,0.,2.4);
TH1D *Jet2EtapassHT = new TH1D("Jet2EtapassHT","", 12,0.,2.4);
TH1D *Jet2EtapassMET = new TH1D("Jet2EtapassMET","", 12,0.,2.4);
Jet2Etaall->GetXaxis()->SetTitle("jet 2 #eta"); Jet2Etaall->GetYaxis()->SetTitle("all events / 0.2");
Jet2EtaallHT->GetXaxis()->SetTitle("jet 2 #eta"); Jet2EtaallHT->GetYaxis()->SetTitle("all events / 0.2");
Jet2EtaallMET->GetXaxis()->SetTitle("jet 2 #eta"); Jet2EtaallMET->GetYaxis()->SetTitle("all events / 0.2");
Jet2Etapass->GetXaxis()->SetTitle("jet 2 #eta"); Jet2Etapass->GetYaxis()->SetTitle("passing events / 0.2");
Jet2EtapassHT->GetXaxis()->SetTitle("jet 2 #eta"); Jet2EtapassHT->GetYaxis()->SetTitle("passing events / 0.2");
Jet2EtapassMET->GetXaxis()->SetTitle("jet 2 #eta"); Jet2EtapassMET->GetYaxis()->SetTitle("passing events / 0.2");
TH1D *Jet3Etaall = new TH1D("Jet3Etaall","",12,0.,2.4);
TH1D *Jet3EtaallHT = new TH1D("Jet3EtaallHT","", 12,0.,2.4);
TH1D *Jet3EtaallMET = new TH1D("Jet3EtaallMET","", 12,0.,2.4);
TH1D *Jet3Etapass = new TH1D("Jet3Etapass","", 12,0.,2.4);
TH1D *Jet3EtapassHT = new TH1D("Jet3EtapassHT","", 12,0.,2.4);
TH1D *Jet3EtapassMET = new TH1D("Jet3EtapassMET","", 12,0.,2.4);
Jet3Etaall->GetXaxis()->SetTitle("jet 3 #eta"); Jet3Etaall->GetYaxis()->SetTitle("all events / 0.2");
Jet3EtaallHT->GetXaxis()->SetTitle("jet 3 #eta"); Jet3EtaallHT->GetYaxis()->SetTitle("all events / 0.2");
Jet3EtaallMET->GetXaxis()->SetTitle("jet 3 #eta"); Jet3EtaallMET->GetYaxis()->SetTitle("all events / 0.2");
Jet3Etapass->GetXaxis()->SetTitle("jet 3 #eta"); Jet3Etapass->GetYaxis()->SetTitle("passing events / 0.2");
Jet3EtapassHT->GetXaxis()->SetTitle("jet 3 #eta"); Jet3EtapassHT->GetYaxis()->SetTitle("passing events / 0.2");
Jet3EtapassMET->GetXaxis()->SetTitle("jet 3 #eta"); Jet3EtapassMET->GetYaxis()->SetTitle("passing events / 0.2");
TH1D *Jet4Etaall = new TH1D("Jet4Etaall","",12,0.,2.4);
TH1D *Jet4EtaallHT = new TH1D("Jet4EtaallHT","", 12,0.,2.4);
TH1D *Jet4EtaallMET = new TH1D("Jet4EtaallMET","", 12,0.,2.4);
TH1D *Jet4Etapass = new TH1D("Jet4Etapass","", 12,0.,2.4);
TH1D *Jet4EtapassHT = new TH1D("Jet4EtapassHT","", 12,0.,2.4);
TH1D *Jet4EtapassMET = new TH1D("Jet4EtapassMET","", 12,0.,2.4);
Jet4Etaall->GetXaxis()->SetTitle("jet 4 #eta"); Jet4Etaall->GetYaxis()->SetTitle("all events / 0.2");
Jet4EtaallHT->GetXaxis()->SetTitle("jet 4 #eta"); Jet4EtaallHT->GetYaxis()->SetTitle("all events / 0.2");
Jet4EtaallMET->GetXaxis()->SetTitle("jet 4 #eta"); Jet4EtaallMET->GetYaxis()->SetTitle("all events / 0.2");
Jet4Etapass->GetXaxis()->SetTitle("jet 4 #eta"); Jet4Etapass->GetYaxis()->SetTitle("passing events / 0.2");
Jet4EtapassHT->GetXaxis()->SetTitle("jet 4 #eta"); Jet4EtapassHT->GetYaxis()->SetTitle("passing events / 0.2");
Jet4EtapassMET->GetXaxis()->SetTitle("jet 4 #eta"); Jet4EtapassMET->GetYaxis()->SetTitle("passing events / 0.2");
const int netabin = 12;
double etabin[netabin+1];
for(int n = 0; n<=netabin;++n) etabin[n] = (double)n*0.2;
const int nptbin = 14;
double ptbin[nptbin+1] = {20,30,40,50,75,100,125,150,175,200,250,300,400,500,1000};
TH2D *JetPtEtaall = new TH2D("JetPtEtaall","",nptbin, ptbin,netabin,etabin);
TH2D *JetPtEtaallHT = new TH2D("JetPtEtaallHT","",nptbin, ptbin,netabin,etabin);
TH2D *JetPtEtaallMET = new TH2D("JetPtEtaallMET","",nptbin, ptbin,netabin,etabin);
TH2D *JetPtEtapass = new TH2D("JetPtEtapass","",nptbin, ptbin,netabin,etabin);
TH2D *JetPtEtapassHT = new TH2D("JetPtEtapassHT","",nptbin, ptbin,netabin,etabin);
TH2D *JetPtEtapassMET = new TH2D("JetPtEtapassMET","",nptbin, ptbin,netabin,etabin);
JetPtEtaall->GetXaxis()->SetTitle("jet p_{T} [GeV]"); JetPtEtaall->GetYaxis()->SetTitle("jet #eta"); JetPtEtaall->GetZaxis()->SetTitle("all events");
JetPtEtaallHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); JetPtEtaallHT->GetYaxis()->SetTitle("jet #eta"); JetPtEtaallHT->GetZaxis()->SetTitle("all events");
JetPtEtaallMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); JetPtEtaallMET->GetYaxis()->SetTitle("jet #eta"); JetPtEtaallMET->GetZaxis()->SetTitle("all events");
JetPtEtapass->GetXaxis()->SetTitle("jet p_{T} [GeV]"); JetPtEtapass->GetYaxis()->SetTitle("jet #eta"); JetPtEtapass->GetZaxis()->SetTitle("passing events");
JetPtEtapassHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); JetPtEtapassHT->GetYaxis()->SetTitle("jet #eta"); JetPtEtapassHT->GetZaxis()->SetTitle("passing events");
JetPtEtapassMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); JetPtEtapassMET->GetYaxis()->SetTitle("jet #eta"); JetPtEtapassMET->GetZaxis()->SetTitle("passing events");
TH1D *JetEtaall = new TH1D("JetEtaall","",netabin,etabin);
TH1D *JetEtaallHT = new TH1D("JetEtaallHT","",netabin,etabin);
TH1D *JetEtaallMET = new TH1D("JetEtaallMET","",netabin,etabin);
TH1D *JetEtapass = new TH1D("JetEtapass","",netabin,etabin);
TH1D *JetEtapassHT = new TH1D("JetEtapassHT","",netabin,etabin);
TH1D *JetEtapassMET = new TH1D("JetEtapassMET","",netabin,etabin);
JetEtaall->GetXaxis()->SetTitle("jet #eta");    JetEtaall->GetYaxis()->SetTitle("all events");
JetEtaallHT->GetXaxis()->SetTitle("jet #eta");  JetEtaallHT->GetYaxis()->SetTitle("all events");
JetEtaallMET->GetXaxis()->SetTitle("jet #eta"); JetEtaallMET->GetYaxis()->SetTitle("all events");
JetEtapass->GetXaxis()->SetTitle("jet #eta");   JetEtapass->GetYaxis()->SetTitle("passing events");
JetEtapassHT->GetXaxis()->SetTitle("jet #eta"); JetEtapassHT->GetYaxis()->SetTitle("passing events");
JetEtapassMET->GetXaxis()->SetTitle("jet #eta");JetEtapassMET->GetYaxis()->SetTitle("passing events");
TH1D *JetPtall = new TH1D("JetPtall","",nptbin, ptbin);
TH1D *JetPtallHT = new TH1D("JetPtallHT","",nptbin, ptbin);
TH1D *JetPtallMET = new TH1D("JetPtallMET","",nptbin, ptbin);
TH1D *JetPtpass = new TH1D("JetPtpass","",nptbin, ptbin);
TH1D *JetPtpassHT = new TH1D("JetPtpassHT","",nptbin, ptbin);
TH1D *JetPtpassMET = new TH1D("JetPtpassMET","",nptbin, ptbin);
JetPtall->GetXaxis()->SetTitle("jet p_{T} [GeV]");     JetPtall->GetYaxis()->SetTitle("all events");
JetPtallHT->GetXaxis()->SetTitle("jet p_{T} [GeV]");   JetPtallHT->GetYaxis()->SetTitle("all events");
JetPtallMET->GetXaxis()->SetTitle("jet p_{T} [GeV]");  JetPtallMET->GetYaxis()->SetTitle("all events");
JetPtpass->GetXaxis()->SetTitle("jet p_{T} [GeV]");    JetPtpass->GetYaxis()->SetTitle("passing events");
JetPtpassHT->GetXaxis()->SetTitle("jet p_{T} [GeV]");  JetPtpassHT->GetYaxis()->SetTitle("passing events");
JetPtpassMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); JetPtpassMET->GetYaxis()->SetTitle("passing events");
TH2D *SJetPt40Etaall = new TH2D("SJetPt40Etaall","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt40EtaallHT = new TH2D("SJetPt40EtaallHT","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt40EtaallMET = new TH2D("SJetPt40EtaallMET","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt40Etapass = new TH2D("SJetPt40Etapass","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt40EtapassHT = new TH2D("SJetPt40EtapassHT","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt40EtapassMET = new TH2D("SJetPt40EtapassMET","",nptbin, ptbin,netabin,etabin);
SJetPt40Etaall->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt40Etaall->GetYaxis()->SetTitle("jet #eta"); SJetPt40Etaall->GetZaxis()->SetTitle("all events");
SJetPt40EtaallHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt40EtaallHT->GetYaxis()->SetTitle("jet #eta"); SJetPt40EtaallHT->GetZaxis()->SetTitle("all events");
SJetPt40EtaallMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt40EtaallMET->GetYaxis()->SetTitle("jet #eta"); SJetPt40EtaallMET->GetZaxis()->SetTitle("all events");
SJetPt40Etapass->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt40Etapass->GetYaxis()->SetTitle("jet #eta"); SJetPt40Etapass->GetZaxis()->SetTitle("passing events");
SJetPt40EtapassHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt40EtapassHT->GetYaxis()->SetTitle("jet #eta"); SJetPt40EtapassHT->GetZaxis()->SetTitle("passing events");
SJetPt40EtapassMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt40EtapassMET->GetYaxis()->SetTitle("jet #eta"); SJetPt40EtapassMET->GetZaxis()->SetTitle("passing events");
TH2D *SJetPt60Etaall = new TH2D("SJetPt60Etaall","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt60EtaallHT = new TH2D("SJetPt60EtaallHT","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt60EtaallMET = new TH2D("SJetPt60EtaallMET","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt60Etapass = new TH2D("SJetPt60Etapass","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt60EtapassHT = new TH2D("SJetPt60EtapassHT","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt60EtapassMET = new TH2D("SJetPt60EtapassMET","",nptbin, ptbin,netabin,etabin);
SJetPt60Etaall->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt60Etaall->GetYaxis()->SetTitle("jet #eta"); SJetPt60Etaall->GetZaxis()->SetTitle("all events");
SJetPt60EtaallHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt60EtaallHT->GetYaxis()->SetTitle("jet #eta"); SJetPt60EtaallHT->GetZaxis()->SetTitle("all events");
SJetPt60EtaallMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt60EtaallMET->GetYaxis()->SetTitle("jet #eta"); SJetPt60EtaallMET->GetZaxis()->SetTitle("all events");
SJetPt60Etapass->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt60Etapass->GetYaxis()->SetTitle("jet #eta"); SJetPt60Etapass->GetZaxis()->SetTitle("passing events");
SJetPt60EtapassHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt60EtapassHT->GetYaxis()->SetTitle("jet #eta"); SJetPt60EtapassHT->GetZaxis()->SetTitle("passing events");
SJetPt60EtapassMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt60EtapassMET->GetYaxis()->SetTitle("jet #eta"); SJetPt60EtapassMET->GetZaxis()->SetTitle("passing events");
TH2D *SJetPt80Etaall = new TH2D("SJetPt80Etaall","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt80EtaallHT = new TH2D("SJetPt80EtaallHT","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt80EtaallMET = new TH2D("SJetPt80EtaallMET","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt80Etapass = new TH2D("SJetPt80Etapass","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt80EtapassHT = new TH2D("SJetPt80EtapassHT","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt80EtapassMET = new TH2D("SJetPt80EtapassMET","",nptbin, ptbin,netabin,etabin);
SJetPt80Etaall->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt80Etaall->GetYaxis()->SetTitle("jet #eta"); SJetPt80Etaall->GetZaxis()->SetTitle("all events");
SJetPt80EtaallHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt80EtaallHT->GetYaxis()->SetTitle("jet #eta"); SJetPt80EtaallHT->GetZaxis()->SetTitle("all events");
SJetPt80EtaallMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt80EtaallMET->GetYaxis()->SetTitle("jet #eta"); SJetPt80EtaallMET->GetZaxis()->SetTitle("all events");
SJetPt80Etapass->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt80Etapass->GetYaxis()->SetTitle("jet #eta"); SJetPt80Etapass->GetZaxis()->SetTitle("passing events");
SJetPt80EtapassHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt80EtapassHT->GetYaxis()->SetTitle("jet #eta"); SJetPt80EtapassHT->GetZaxis()->SetTitle("passing events");
SJetPt80EtapassMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt80EtapassMET->GetYaxis()->SetTitle("jet #eta"); SJetPt80EtapassMET->GetZaxis()->SetTitle("passing events");
TH2D *SJetPt100Etaall = new TH2D("SJetPt100Etaall","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt100EtaallHT = new TH2D("SJetPt100EtaallHT","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt100EtaallMET = new TH2D("SJetPt100EtaallMET","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt100Etapass = new TH2D("SJetPt100Etapass","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt100EtapassHT = new TH2D("SJetPt100EtapassHT","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt100EtapassMET = new TH2D("SJetPt100EtapassMET","",nptbin, ptbin,netabin,etabin);
SJetPt100Etaall->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt100Etaall->GetYaxis()->SetTitle("jet #eta"); SJetPt100Etaall->GetZaxis()->SetTitle("all events");
SJetPt100EtaallHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt100EtaallHT->GetYaxis()->SetTitle("jet #eta"); SJetPt100EtaallHT->GetZaxis()->SetTitle("all events");
SJetPt100EtaallMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt100EtaallMET->GetYaxis()->SetTitle("jet #eta"); SJetPt100EtaallMET->GetZaxis()->SetTitle("all events");
SJetPt100Etapass->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt100Etapass->GetYaxis()->SetTitle("jet #eta"); SJetPt100Etapass->GetZaxis()->SetTitle("passing events");
SJetPt100EtapassHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt100EtapassHT->GetYaxis()->SetTitle("jet #eta"); SJetPt100EtapassHT->GetZaxis()->SetTitle("passing events");
SJetPt100EtapassMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt100EtapassMET->GetYaxis()->SetTitle("jet #eta"); SJetPt100EtapassMET->GetZaxis()->SetTitle("passing events");
TH2D *SJetPt150Etaall = new TH2D("SJetPt150Etaall","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt150EtaallHT = new TH2D("SJetPt150EtaallHT","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt150EtaallMET = new TH2D("SJetPt150EtaallMET","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt150Etapass = new TH2D("SJetPt150Etapass","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt150EtapassHT = new TH2D("SJetPt150EtapassHT","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt150EtapassMET = new TH2D("SJetPt150EtapassMET","",nptbin, ptbin,netabin,etabin);
SJetPt150Etaall->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt150Etaall->GetYaxis()->SetTitle("jet #eta"); SJetPt150Etaall->GetZaxis()->SetTitle("all events");
SJetPt150EtaallHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt150EtaallHT->GetYaxis()->SetTitle("jet #eta"); SJetPt150EtaallHT->GetZaxis()->SetTitle("all events");
SJetPt150EtaallMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt150EtaallMET->GetYaxis()->SetTitle("jet #eta"); SJetPt150EtaallMET->GetZaxis()->SetTitle("all events");
SJetPt150Etapass->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt150Etapass->GetYaxis()->SetTitle("jet #eta"); SJetPt150Etapass->GetZaxis()->SetTitle("passing events");
SJetPt150EtapassHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt150EtapassHT->GetYaxis()->SetTitle("jet #eta"); SJetPt150EtapassHT->GetZaxis()->SetTitle("passing events");
SJetPt150EtapassMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt150EtapassMET->GetYaxis()->SetTitle("jet #eta"); SJetPt150EtapassMET->GetZaxis()->SetTitle("passing events");
TH2D *SJetPt200Etaall = new TH2D("SJetPt200Etaall","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt200EtaallHT = new TH2D("SJetPt200EtaallHT","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt200EtaallMET = new TH2D("SJetPt200EtaallMET","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt200Etapass = new TH2D("SJetPt200Etapass","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt200EtapassHT = new TH2D("SJetPt200EtapassHT","",nptbin, ptbin,netabin,etabin);
TH2D *SJetPt200EtapassMET = new TH2D("SJetPt200EtapassMET","",nptbin, ptbin,netabin,etabin);
SJetPt200Etaall->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt200Etaall->GetYaxis()->SetTitle("jet #eta"); SJetPt200Etaall->GetZaxis()->SetTitle("all events");
SJetPt200EtaallHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt200EtaallHT->GetYaxis()->SetTitle("jet #eta"); SJetPt200EtaallHT->GetZaxis()->SetTitle("all events");
SJetPt200EtaallMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt200EtaallMET->GetYaxis()->SetTitle("jet #eta"); SJetPt200EtaallMET->GetZaxis()->SetTitle("all events");
SJetPt200Etapass->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt200Etapass->GetYaxis()->SetTitle("jet #eta"); SJetPt200Etapass->GetZaxis()->SetTitle("passing events");
SJetPt200EtapassHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt200EtapassHT->GetYaxis()->SetTitle("jet #eta"); SJetPt200EtapassHT->GetZaxis()->SetTitle("passing events");
SJetPt200EtapassMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt200EtapassMET->GetYaxis()->SetTitle("jet #eta"); SJetPt200EtapassMET->GetZaxis()->SetTitle("passing events");
TH1D *SJetPt40all = new TH1D("SJetPt40all","",nptbin, ptbin);
TH1D *SJetPt40allHT = new TH1D("SJetPt40allHT","",nptbin, ptbin);
TH1D *SJetPt40allMET = new TH1D("SJetPt40allMET","",nptbin, ptbin);
TH1D *SJetPt40pass = new TH1D("SJetPt40pass","",nptbin, ptbin);
TH1D *SJetPt40passHT = new TH1D("SJetPt40passHT","",nptbin, ptbin);
TH1D *SJetPt40passMET = new TH1D("SJetPt40passMET","",nptbin, ptbin);
SJetPt40all->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt40all->GetYaxis()->SetTitle("all events");
SJetPt40allHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt40allHT->GetYaxis()->SetTitle("all events");
SJetPt40allMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt40allMET->GetYaxis()->SetTitle("all events");
SJetPt40pass->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt40pass->GetYaxis()->SetTitle("passing events");
SJetPt40passHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt40passHT->GetYaxis()->SetTitle("passing events");
SJetPt40passMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt40passMET->GetYaxis()->SetTitle("passing events");
TH1D *SJetPt60all = new TH1D("SJetPt60all","",nptbin, ptbin);
TH1D *SJetPt60allHT = new TH1D("SJetPt60allHT","",nptbin, ptbin);
TH1D *SJetPt60allMET = new TH1D("SJetPt60allMET","",nptbin, ptbin);
TH1D *SJetPt60pass = new TH1D("SJetPt60pass","",nptbin, ptbin);
TH1D *SJetPt60passHT = new TH1D("SJetPt60passHT","",nptbin, ptbin);
TH1D *SJetPt60passMET = new TH1D("SJetPt60passMET","",nptbin, ptbin);
SJetPt60all->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt60all->GetYaxis()->SetTitle("all events");
SJetPt60allHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt60allHT->GetYaxis()->SetTitle("all events");
SJetPt60allMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt60allMET->GetYaxis()->SetTitle("all events");
SJetPt60pass->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt60pass->GetYaxis()->SetTitle("passing events");
SJetPt60passHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt60passHT->GetYaxis()->SetTitle("passing events");
SJetPt60passMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt60passMET->GetYaxis()->SetTitle("passing events");
TH1D *SJetPt80all = new TH1D("SJetPt80all","",nptbin, ptbin);
TH1D *SJetPt80allHT = new TH1D("SJetPt80allHT","",nptbin, ptbin);
TH1D *SJetPt80allMET = new TH1D("SJetPt80allMET","",nptbin, ptbin);
TH1D *SJetPt80pass = new TH1D("SJetPt80pass","",nptbin, ptbin);
TH1D *SJetPt80passHT = new TH1D("SJetPt80passHT","",nptbin, ptbin);
TH1D *SJetPt80passMET = new TH1D("SJetPt80passMET","",nptbin, ptbin);
SJetPt80all->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt80all->GetYaxis()->SetTitle("all events");
SJetPt80allHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt80allHT->GetYaxis()->SetTitle("all events");
SJetPt80allMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt80allMET->GetYaxis()->SetTitle("all events");
SJetPt80pass->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt80pass->GetYaxis()->SetTitle("passing events");
SJetPt80passHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt80passHT->GetYaxis()->SetTitle("passing events");
SJetPt80passMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt80passMET->GetYaxis()->SetTitle("passing events");
TH1D *SJetPt100all = new TH1D("SJetPt100all","",nptbin, ptbin);
TH1D *SJetPt100allHT = new TH1D("SJetPt100allHT","",nptbin, ptbin);
TH1D *SJetPt100allMET = new TH1D("SJetPt100allMET","",nptbin, ptbin);
TH1D *SJetPt100pass = new TH1D("SJetPt100pass","",nptbin, ptbin);
TH1D *SJetPt100passHT = new TH1D("SJetPt100passHT","",nptbin, ptbin);
TH1D *SJetPt100passMET = new TH1D("SJetPt100passMET","",nptbin, ptbin);
SJetPt100all->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt100all->GetYaxis()->SetTitle("all events");
SJetPt100allHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt100allHT->GetYaxis()->SetTitle("all events");
SJetPt100allMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt100allMET->GetYaxis()->SetTitle("all events");
SJetPt100pass->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt100pass->GetYaxis()->SetTitle("passing events");
SJetPt100passHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt100passHT->GetYaxis()->SetTitle("passing events");
SJetPt100passMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt100passMET->GetYaxis()->SetTitle("passing events");
TH1D *SJetPt150all = new TH1D("SJetPt150all","",nptbin, ptbin);
TH1D *SJetPt150allHT = new TH1D("SJetPt150allHT","",nptbin, ptbin);
TH1D *SJetPt150allMET = new TH1D("SJetPt150allMET","",nptbin, ptbin);
TH1D *SJetPt150pass = new TH1D("SJetPt150pass","",nptbin, ptbin);
TH1D *SJetPt150passHT = new TH1D("SJetPt150passHT","",nptbin, ptbin);
TH1D *SJetPt150passMET = new TH1D("SJetPt150passMET","",nptbin, ptbin);
SJetPt150all->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt150all->GetYaxis()->SetTitle("all events");
SJetPt150allHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt150allHT->GetYaxis()->SetTitle("all events");
SJetPt150allMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt150allMET->GetYaxis()->SetTitle("all events");
SJetPt150pass->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt150pass->GetYaxis()->SetTitle("passing events");
SJetPt150passHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt150passHT->GetYaxis()->SetTitle("passing events");
SJetPt150passMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt150passMET->GetYaxis()->SetTitle("passing events");
TH1D *SJetPt200all = new TH1D("SJetPt200all","",nptbin, ptbin);
TH1D *SJetPt200allHT = new TH1D("SJetPt200allHT","",nptbin, ptbin);
TH1D *SJetPt200allMET = new TH1D("SJetPt200allMET","",nptbin, ptbin);
TH1D *SJetPt200pass = new TH1D("SJetPt200pass","",nptbin, ptbin);
TH1D *SJetPt200passHT = new TH1D("SJetPt200passHT","",nptbin, ptbin);
TH1D *SJetPt200passMET = new TH1D("SJetPt200passMET","",nptbin, ptbin);
SJetPt200all->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt200all->GetYaxis()->SetTitle("all events");
SJetPt200allHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt200allHT->GetYaxis()->SetTitle("all events");
SJetPt200allMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt200allMET->GetYaxis()->SetTitle("all events");
SJetPt200pass->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt200pass->GetYaxis()->SetTitle("passing events");
SJetPt200passHT->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt200passHT->GetYaxis()->SetTitle("passing events");
SJetPt200passMET->GetXaxis()->SetTitle("jet p_{T} [GeV]"); SJetPt200passMET->GetYaxis()->SetTitle("passing events");
TH1D *SJet40Etaall = new TH1D("SJet40Etaall","",netabin, etabin);
TH1D *SJet40EtaallHT = new TH1D("SJet40EtaallHT","",netabin,etabin);
TH1D *SJet40EtaallMET = new TH1D("SJet40EtaallMET","",netabin,etabin);
TH1D *SJet40Etapass = new TH1D("SJet40Etapass","",netabin,etabin);
TH1D *SJet40EtapassHT = new TH1D("SJet40EtapassHT","",netabin,etabin);
TH1D *SJet40EtapassMET = new TH1D("SJet40EtapassMET","",netabin,etabin);
SJet40Etaall->GetXaxis()->SetTitle("jet #eta"); SJet40Etaall->GetYaxis()->SetTitle("all events");
SJet40EtaallHT->GetXaxis()->SetTitle("jet #eta"); SJet40EtaallHT->GetYaxis()->SetTitle("all events");
SJet40EtaallMET->GetXaxis()->SetTitle("jet #eta"); SJet40EtaallMET->GetYaxis()->SetTitle("all events");
SJet40Etapass->GetXaxis()->SetTitle("jet #eta"); SJet40Etapass->GetYaxis()->SetTitle("passing events");
SJet40EtapassHT->GetXaxis()->SetTitle("jet #eta"); SJet40EtapassHT->GetYaxis()->SetTitle("passing events");
SJet40EtapassMET->GetXaxis()->SetTitle("jet #eta"); SJet40EtapassMET->GetYaxis()->SetTitle("passing events");
TH1D *SJet60Etaall = new TH1D("SJet60Etaall","",netabin,etabin);
TH1D *SJet60EtaallHT = new TH1D("SJet60EtaallHT","",netabin,etabin);
TH1D *SJet60EtaallMET = new TH1D("SJet60EtaallMET","",netabin,etabin);
TH1D *SJet60Etapass = new TH1D("SJet60Etapass","",netabin,etabin);
TH1D *SJet60EtapassHT = new TH1D("SJet60EtapassHT","",netabin,etabin);
TH1D *SJet60EtapassMET = new TH1D("SJet60EtapassMET","",netabin,etabin);
SJet60Etaall->GetXaxis()->SetTitle("jet #eta"); SJet60Etaall->GetYaxis()->SetTitle("all events");
SJet60EtaallHT->GetXaxis()->SetTitle("jet #eta"); SJet60EtaallHT->GetYaxis()->SetTitle("all events");
SJet60EtaallMET->GetXaxis()->SetTitle("jet #eta"); SJet60EtaallMET->GetYaxis()->SetTitle("all events");
SJet60Etapass->GetXaxis()->SetTitle("jet #eta"); SJet60Etapass->GetYaxis()->SetTitle("passing events");
SJet60EtapassHT->GetXaxis()->SetTitle("jet #eta"); SJet60EtapassHT->GetYaxis()->SetTitle("passing events");
SJet60EtapassMET->GetXaxis()->SetTitle("jet #eta"); SJet60EtapassMET->GetYaxis()->SetTitle("passing events");
TH1D *SJet80Etaall = new TH1D("SJet80Etaall","",netabin,etabin);
TH1D *SJet80EtaallHT = new TH1D("SJet80EtaallHT","",netabin,etabin);
TH1D *SJet80EtaallMET = new TH1D("SJet80EtaallMET","",netabin,etabin);
TH1D *SJet80Etapass = new TH1D("SJet80Etapass","",netabin,etabin);
TH1D *SJet80EtapassHT = new TH1D("SJet80EtapassHT","",netabin,etabin);
TH1D *SJet80EtapassMET = new TH1D("SJet80EtapassMET","",netabin,etabin);
SJet80Etaall->GetXaxis()->SetTitle("jet #eta"); SJet80Etaall->GetYaxis()->SetTitle("all events");
SJet80EtaallHT->GetXaxis()->SetTitle("jet #eta"); SJet80EtaallHT->GetYaxis()->SetTitle("all events");
SJet80EtaallMET->GetXaxis()->SetTitle("jet #eta"); SJet80EtaallMET->GetYaxis()->SetTitle("all events");
SJet80Etapass->GetXaxis()->SetTitle("jet #eta"); SJet80Etapass->GetYaxis()->SetTitle("passing events");
SJet80EtapassHT->GetXaxis()->SetTitle("jet #eta"); SJet80EtapassHT->GetYaxis()->SetTitle("passing events");
SJet80EtapassMET->GetXaxis()->SetTitle("jet #eta"); SJet80EtapassMET->GetYaxis()->SetTitle("passing events");
TH1D *SJet100Etaall = new TH1D("SJet100Etaall","",netabin,etabin);
TH1D *SJet100EtaallHT = new TH1D("SJet100EtaallHT","",netabin,etabin);
TH1D *SJet100EtaallMET = new TH1D("SJet100EtaallMET","",netabin,etabin);
TH1D *SJet100Etapass = new TH1D("SJet100Etapass","",netabin,etabin);
TH1D *SJet100EtapassHT = new TH1D("SJet100EtapassHT","",netabin,etabin);
TH1D *SJet100EtapassMET = new TH1D("SJet100EtapassMET","",netabin,etabin);
SJet100Etaall->GetXaxis()->SetTitle("jet #eta"); SJet100Etaall->GetYaxis()->SetTitle("all events");
SJet100EtaallHT->GetXaxis()->SetTitle("jet #eta"); SJet100EtaallHT->GetYaxis()->SetTitle("all events");
SJet100EtaallMET->GetXaxis()->SetTitle("jet #eta"); SJet100EtaallMET->GetYaxis()->SetTitle("all events");
SJet100Etapass->GetXaxis()->SetTitle("jet #eta"); SJet100Etapass->GetYaxis()->SetTitle("passing events");
SJet100EtapassHT->GetXaxis()->SetTitle("jet #eta"); SJet100EtapassHT->GetYaxis()->SetTitle("passing events");
SJet100EtapassMET->GetXaxis()->SetTitle("jet #eta"); SJet100EtapassMET->GetYaxis()->SetTitle("passing events");
TH1D *SJet150Etaall = new TH1D("SJet150Etaall","",netabin,etabin);
TH1D *SJet150EtaallHT = new TH1D("SJet150EtaallHT","",netabin,etabin);
TH1D *SJet150EtaallMET = new TH1D("SJet150EtaallMET","",netabin,etabin);
TH1D *SJet150Etapass = new TH1D("SJet150Etapass","",netabin,etabin);
TH1D *SJet150EtapassHT = new TH1D("SJet150EtapassHT","",netabin,etabin);
TH1D *SJet150EtapassMET = new TH1D("SJet150EtapassMET","",netabin,etabin);
SJet150Etaall->GetXaxis()->SetTitle("jet #eta"); SJet150Etaall->GetYaxis()->SetTitle("all events");
SJet150EtaallHT->GetXaxis()->SetTitle("jet #eta"); SJet150EtaallHT->GetYaxis()->SetTitle("all events");
SJet150EtaallMET->GetXaxis()->SetTitle("jet #eta"); SJet150EtaallMET->GetYaxis()->SetTitle("all events");
SJet150Etapass->GetXaxis()->SetTitle("jet #eta"); SJet150Etapass->GetYaxis()->SetTitle("passing events");
SJet150EtapassHT->GetXaxis()->SetTitle("jet #eta"); SJet150EtapassHT->GetYaxis()->SetTitle("passing events");
SJet150EtapassMET->GetXaxis()->SetTitle("jet #eta"); SJet150EtapassMET->GetYaxis()->SetTitle("passing events");
TH1D *SJet200Etaall = new TH1D("SJet200Etaall","",netabin,etabin);
TH1D *SJet200EtaallHT = new TH1D("SJet200EtaallHT","",netabin,etabin);
TH1D *SJet200EtaallMET = new TH1D("SJet200EtaallMET","",netabin,etabin);
TH1D *SJet200Etapass = new TH1D("SJet200Etapass","",netabin,etabin);
TH1D *SJet200EtapassHT = new TH1D("SJet200EtapassHT","",netabin,etabin);
TH1D *SJet200EtapassMET = new TH1D("SJet200EtapassMET","",netabin,etabin);
SJet200Etaall->GetXaxis()->SetTitle("jet #eta"); SJet200Etaall->GetYaxis()->SetTitle("all events");
SJet200EtaallHT->GetXaxis()->SetTitle("jet #eta"); SJet200EtaallHT->GetYaxis()->SetTitle("all events");
SJet200EtaallMET->GetXaxis()->SetTitle("jet #eta"); SJet200EtaallMET->GetYaxis()->SetTitle("all events");
SJet200Etapass->GetXaxis()->SetTitle("jet #eta"); SJet200Etapass->GetYaxis()->SetTitle("passing events");
SJet200EtapassHT->GetXaxis()->SetTitle("jet #eta"); SJet200EtapassHT->GetYaxis()->SetTitle("passing events");
SJet200EtapassMET->GetXaxis()->SetTitle("jet #eta"); SJet200EtapassMET->GetYaxis()->SetTitle("passing events");


gROOT->cd();
MT2tree* fMT2tree = new MT2tree();
c->SetBranchAddress("MT2tree", &fMT2tree);
Float_t fTOBTECTagger;// = new Float_t;
c->SetBranchAddress("TOBTECTagger", &fTOBTECTagger);
Long64_t nentries =  c->GetEntries();
Long64_t nbytes = 0, nb = 0;
int nev =0;
//otional cuts
//cutallnomindphinolepvetoandphotons = "misc.HT<750&&" + cutallnomindphinolepvetoandphotons;
//cutallnomindphinolepvetoandphotons = "misc.HT>750&&" + cutallnomindphinolepvetoandphotons;
//cutallnomindphinolepvetoandphotons = "misc.HT>750&&misc.HT<1200&&" + cutallnomindphinolepvetoandphotons;
//cutallnomindphinolepvetoandphotons = "misc.HT>1200&&" + cutallnomindphinolepvetoandphotons;
c->Draw(">>selList", cutallnomindphinolepvetoandphotons);//cuts HERE
TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
c->SetEventList(myEvtList);
int counter=0;
cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
//run over loosest selection
while(myEvtList->GetEntry(counter++) !=-1){	
	int jentry = myEvtList->GetEntry(counter-1);
	nb =  c->GetEntry(jentry);   nbytes += nb;
	if ( counter % 500 == 0  )  cout << "+++ Proccessing event " << counter << " " << fTOBTECTagger << endl;

	bool HT = false;
	bool MET = false;
	bool tagged = false;
	if(fTOBTECTagger>8) tagged = true;
	if(((fMT2tree->misc.MT2>200&&fMT2tree->misc.HT>=450&&fMT2tree->misc.HT<750&&fMT2tree->misc.MET>200&&((fMT2tree->trigger.HLT_PFMET150_v2==1||fMT2tree->trigger.HLT_PFMET150_v3==1||fMT2tree->trigger.HLT_PFMET150_v4==1||fMT2tree->trigger.HLT_PFMET150_v5==1||fMT2tree->trigger.HLT_PFMET150_v6==1||fMT2tree->trigger.HLT_PFMET150_v7==1)||(fMT2tree->trigger.HLT_PFHT350_PFMET100_v3==1||fMT2tree->trigger.HLT_PFHT350_PFMET100_v4==1||fMT2tree->trigger.HLT_PFHT350_PFMET100_v5==1||fMT2tree->trigger.HLT_PFHT350_PFMET100_v6==1||fMT2tree->trigger.HLT_PFHT350_PFMET100_v7==1||fMT2tree->trigger.HLT_PFNoPUHT350_PFMET100_v1==1||fMT2tree->trigger.HLT_PFNoPUHT350_PFMET100_v3==1||fMT2tree->trigger.HLT_PFNoPUHT350_PFMET100_v4==1))))) MET = true;
	if(((fMT2tree->misc.MT2>100&&fMT2tree->misc.HT>=750&&fMT2tree->misc.MET>30&&(fMT2tree->trigger.HLT_PFHT650_v5==1||fMT2tree->trigger.HLT_PFHT650_v6==1||fMT2tree->trigger.HLT_PFHT650_v7==1||fMT2tree->trigger.HLT_PFHT650_v8==1||fMT2tree->trigger.HLT_PFHT650_v9==1||fMT2tree->trigger.HLT_PFNoPUHT650_v1==1||fMT2tree->trigger.HLT_PFNoPUHT650_v3==1||fMT2tree->trigger.HLT_PFNoPUHT650_v4==1)))) HT = true;

	//fill all histograms
	MT2all   ->Fill(fMT2tree->misc.MT2);
	HTall    ->Fill(fMT2tree->misc.HT);
	METall   ->Fill(fMT2tree->misc.MET);
	NJetsall ->Fill(fMT2tree->NJetsIDLoose40);
	NBJetsall->Fill(fMT2tree->NBJets40CSVM);
	if(fMT2tree->jet[0].isPFIDLoose) Jet1Etaall->Fill(fabs(fMT2tree->jet[0].lv.Eta()));
	if(fMT2tree->jet[1].isPFIDLoose) Jet2Etaall->Fill(fabs(fMT2tree->jet[1].lv.Eta()));
	if(fMT2tree->jet[2].isPFIDLoose) Jet3Etaall->Fill(fabs(fMT2tree->jet[2].lv.Eta()));
	if(fMT2tree->jet[3].isPFIDLoose) Jet4Etaall->Fill(fabs(fMT2tree->jet[3].lv.Eta()));
	if(HT){
		MT2allHT   ->Fill(fMT2tree->misc.MT2);
		HTallHT    ->Fill(fMT2tree->misc.HT);
		METallHT   ->Fill(fMT2tree->misc.MET);
		NJetsallHT ->Fill(fMT2tree->NJetsIDLoose40);
		NBJetsallHT->Fill(fMT2tree->NBJets40CSVM);
		if(fMT2tree->jet[0].isPFIDLoose) Jet1EtaallHT->Fill(fabs(fMT2tree->jet[0].lv.Eta()));
		if(fMT2tree->jet[1].isPFIDLoose) Jet2EtaallHT->Fill(fabs(fMT2tree->jet[1].lv.Eta()));
		if(fMT2tree->jet[2].isPFIDLoose) Jet3EtaallHT->Fill(fabs(fMT2tree->jet[2].lv.Eta()));
		if(fMT2tree->jet[3].isPFIDLoose) Jet4EtaallHT->Fill(fabs(fMT2tree->jet[3].lv.Eta()));	
	} if(MET){
		MT2allMET   ->Fill(fMT2tree->misc.MT2);
		HTallMET    ->Fill(fMT2tree->misc.HT);
		METallMET   ->Fill(fMT2tree->misc.MET);
		NJetsallMET ->Fill(fMT2tree->NJetsIDLoose40);
		NBJetsallMET->Fill(fMT2tree->NBJets40CSVM);
		if(fMT2tree->jet[0].isPFIDLoose) Jet1EtaallMET->Fill(fabs(fMT2tree->jet[0].lv.Eta()));
		if(fMT2tree->jet[1].isPFIDLoose) Jet2EtaallMET->Fill(fabs(fMT2tree->jet[1].lv.Eta()));
		if(fMT2tree->jet[2].isPFIDLoose) Jet3EtaallMET->Fill(fabs(fMT2tree->jet[2].lv.Eta()));
		if(fMT2tree->jet[3].isPFIDLoose) Jet4EtaallMET->Fill(fabs(fMT2tree->jet[3].lv.Eta()));
	}
	if(!tagged){
		MT2pass   ->Fill(fMT2tree->misc.MT2);
		HTpass    ->Fill(fMT2tree->misc.HT);
		METpass   ->Fill(fMT2tree->misc.MET);
		NJetspass ->Fill(fMT2tree->NJetsIDLoose40);
		NBJetspass->Fill(fMT2tree->NBJets40CSVM);
		if(fMT2tree->jet[0].isPFIDLoose) Jet1Etapass->Fill(fabs(fMT2tree->jet[0].lv.Eta()));
		if(fMT2tree->jet[1].isPFIDLoose) Jet2Etapass->Fill(fabs(fMT2tree->jet[1].lv.Eta()));
		if(fMT2tree->jet[2].isPFIDLoose) Jet3Etapass->Fill(fabs(fMT2tree->jet[2].lv.Eta()));
		if(fMT2tree->jet[3].isPFIDLoose) Jet4Etapass->Fill(fabs(fMT2tree->jet[3].lv.Eta()));
		if(HT){
			MT2passHT   ->Fill(fMT2tree->misc.MT2);
			HTpassHT    ->Fill(fMT2tree->misc.HT);
			METpassHT   ->Fill(fMT2tree->misc.MET);
			NJetspassHT ->Fill(fMT2tree->NJetsIDLoose40);
			NBJetspassHT->Fill(fMT2tree->NBJets40CSVM);
			if(fMT2tree->jet[0].isPFIDLoose) Jet1EtapassHT->Fill(fabs(fMT2tree->jet[0].lv.Eta()));
			if(fMT2tree->jet[1].isPFIDLoose) Jet2EtapassHT->Fill(fabs(fMT2tree->jet[1].lv.Eta()));
			if(fMT2tree->jet[2].isPFIDLoose) Jet3EtapassHT->Fill(fabs(fMT2tree->jet[2].lv.Eta()));
			if(fMT2tree->jet[3].isPFIDLoose) Jet4EtapassHT->Fill(fabs(fMT2tree->jet[3].lv.Eta()));	
		} if(MET){
			MT2passMET   ->Fill(fMT2tree->misc.MT2);
			HTpassMET    ->Fill(fMT2tree->misc.HT);
			METpassMET   ->Fill(fMT2tree->misc.MET);
			NJetspassMET ->Fill(fMT2tree->NJetsIDLoose40);
			NBJetspassMET->Fill(fMT2tree->NBJets40CSVM);
			if(fMT2tree->jet[0].isPFIDLoose) Jet1EtapassMET->Fill(fabs(fMT2tree->jet[0].lv.Eta()));
			if(fMT2tree->jet[1].isPFIDLoose) Jet2EtapassMET->Fill(fabs(fMT2tree->jet[1].lv.Eta()));
			if(fMT2tree->jet[2].isPFIDLoose) Jet3EtapassMET->Fill(fabs(fMT2tree->jet[2].lv.Eta()));
			if(fMT2tree->jet[3].isPFIDLoose) Jet4EtapassMET->Fill(fabs(fMT2tree->jet[3].lv.Eta()));
		}
	}
	//get jet for creating 'reweighting map' - done for several jet pt thresholds
	float pt200(-1), pt150(-1), pt100(-1), pt80(-1), pt60(-1), pt40(-1);
	float eta200(99), eta150(99), eta100(99), eta80(99), eta60(99), eta40(99);
	int n200(-1), n150(-1), n100(-1), n80(-1), n60(-1), n40(-1);
	for(int n = 0; n<fMT2tree->NJets; ++n){
		if(!(fMT2tree->jet[n].isPFIDLoose))     continue;
		if(fMT2tree->jet[n].lv.Pt()*scale<40)         continue;
		if(fabs(fMT2tree->jet[n].lv.Eta())>2.4) continue;
		JetPtEtaall->Fill(fMT2tree->jet[n].lv.Pt()*scale, fabs(fMT2tree->jet[n].lv.Eta()));
		JetEtaall->Fill(fabs(fMT2tree->jet[n].lv.Eta()));
		JetPtall->Fill(fMT2tree->jet[n].lv.Pt()*scale);
		if(HT){
			JetPtEtaallHT->Fill(fMT2tree->jet[n].lv.Pt()*scale, fabs(fMT2tree->jet[n].lv.Eta()));
			JetEtaallHT->Fill(fabs(fMT2tree->jet[n].lv.Eta()));
			JetPtallHT->Fill(fMT2tree->jet[n].lv.Pt()*scale);
		} if(MET){
			JetPtEtaallMET->Fill(fMT2tree->jet[n].lv.Pt()*scale, fabs(fMT2tree->jet[n].lv.Eta()));
			JetEtaallMET->Fill(fabs(fMT2tree->jet[n].lv.Eta()));
			JetPtallMET->Fill(fMT2tree->jet[n].lv.Pt()*scale);
		}
		if(!tagged){
			JetPtEtapass->Fill(fMT2tree->jet[n].lv.Pt()*scale, fabs(fMT2tree->jet[n].lv.Eta()));
			JetEtapass->Fill(fabs(fMT2tree->jet[n].lv.Eta()));
			JetPtpass->Fill(fMT2tree->jet[n].lv.Pt()*scale);
			if(HT){
				JetPtEtapassHT->Fill(fMT2tree->jet[n].lv.Pt()*scale, fabs(fMT2tree->jet[n].lv.Eta()));
				JetEtapassHT->Fill(fabs(fMT2tree->jet[n].lv.Eta()));
				JetPtpassHT->Fill(fMT2tree->jet[n].lv.Pt()*scale);
			} if(MET){
				JetPtEtapassMET->Fill(fMT2tree->jet[n].lv.Pt()*scale, fabs(fMT2tree->jet[n].lv.Eta()));
				JetEtapassMET->Fill(fabs(fMT2tree->jet[n].lv.Eta()));
				JetPtpassMET->Fill(fMT2tree->jet[n].lv.Pt()*scale);
			}
		}
		//create the map by first selecting correct jet
		float deltaeta40 = fabs(1.5-eta40);
		float deltaeta = fabs(1.5-fabs(fMT2tree->jet[n].lv.Eta()));
		float eta = fabs(fMT2tree->jet[n].lv.Eta());
		float pt = fMT2tree->jet[n].lv.Pt()*scale;
		//cout << "pt " << pt << " eta " << eta << " deltaeta " << deltaeta << " deltaeta40 " << deltaeta40 << " pt40 " << pt40 << endl;
		//40
		if(     deltaeta< 0.2&&deltaeta40>=0.2) {pt40 = pt; eta40 = eta; n40 = n;}
		else if(deltaeta< 0.2&&deltaeta40< 0.2) { if(pt>pt40) {pt40 = pt; eta40 = eta; n40 = n;} }
		else if(deltaeta>=0.2&&deltaeta40>=0.2){
		if(     deltaeta< 0.4&&deltaeta40>=0.4) {pt40 = pt; eta40 = eta; n40 = n;}
		else if(deltaeta< 0.4&&deltaeta40< 0.4) { if(pt>pt40) {pt40 = pt; eta40 = eta; n40 = n;} }  
		else if(deltaeta>=0.4&&deltaeta40>=0.4){
		if(     deltaeta< 0.6&&deltaeta40>=0.6) {pt40 = pt; eta40 = eta; n40 = n;}
		else if(deltaeta< 0.6&&deltaeta40< 0.6) { if(pt>pt40) {pt40 = pt; eta40 = eta; n40 = n;} }  
		else if(deltaeta>=0.6&&deltaeta40>=0.6){
		if(     deltaeta< 0.8&&deltaeta40>=0.8) {pt40 = pt; eta40 = eta; n40 = n;}
		else if(deltaeta< 0.8&&deltaeta40< 0.8) { if(pt>pt40) {pt40 = pt; eta40 = eta; n40 = n;} }  
		else if(deltaeta>=0.8&&deltaeta40>=0.8){
		if(     deltaeta< 1.0&&deltaeta40>=1.0) {pt40 = pt; eta40 = eta; n40 = n;}
		else if(deltaeta< 1.0&&deltaeta40< 1.0) { if(pt>pt40) {pt40 = pt; eta40 = eta; n40 = n;} }  
		else if(deltaeta>=1.0&&deltaeta40>=1.0){
		if(     deltaeta< 1.2&&deltaeta40>=1.2) {pt40 = pt; eta40 = eta; n40 = n;}
		else if(deltaeta< 1.2&&deltaeta40< 1.2) { if(pt>pt40) {pt40 = pt; eta40 = eta; n40 = n;} }  
		else if(deltaeta>=1.2&&deltaeta40>=1.2){
		if(     deltaeta< 1.4&&deltaeta40>=1.4) {pt40 = pt; eta40 = eta; n40 = n;}
		else if(deltaeta< 1.4&&deltaeta40< 1.4) { if(pt>pt40) {pt40 = pt; eta40 = eta; n40 = n;} }  
		else if(deltaeta>=1.4&&deltaeta40>=1.4){
		if(     deltaeta< 1.6&&deltaeta40>=1.6) {pt40 = pt; eta40 = eta; n40 = n;}
		else if(deltaeta< 1.6&&deltaeta40< 1.6) { if(pt>pt40) {pt40 = pt; eta40 = eta; n40 = n;} }  }}}}}}}
		if(pt<60) continue;
		float deltaeta60 = fabs(1.5-eta60);
		if(     deltaeta< 0.2&&deltaeta60>=0.2) {pt60 = pt; eta60 = eta; n60 = n;}
		else if(deltaeta< 0.2&&deltaeta60< 0.2) { if(pt>pt60) {pt60 = pt; eta60 = eta; n60 = n;} }
		else if(deltaeta>=0.2&&deltaeta60>=0.2){
		if(     deltaeta< 0.4&&deltaeta60>=0.4) {pt60 = pt; eta60 = eta; n60 = n;}
		else if(deltaeta< 0.4&&deltaeta60< 0.4) { if(pt>pt60) {pt60 = pt; eta60 = eta; n60 = n;} }  
		else if(deltaeta>=0.4&&deltaeta60>=0.4){
		if(     deltaeta< 0.6&&deltaeta60>=0.6) {pt60 = pt; eta60 = eta; n60 = n;}
		else if(deltaeta< 0.6&&deltaeta60< 0.6) { if(pt>pt60) {pt60 = pt; eta60 = eta; n60 = n;} }  
		else if(deltaeta>=0.6&&deltaeta60>=0.6){
		if(     deltaeta< 0.8&&deltaeta60>=0.8) {pt60 = pt; eta60 = eta; n60 = n;}
		else if(deltaeta< 0.8&&deltaeta60< 0.8) { if(pt>pt60) {pt60 = pt; eta60 = eta; n60 = n;} }  
		else if(deltaeta>=0.8&&deltaeta60>=0.8){
		if(     deltaeta< 1.0&&deltaeta60>=1.0) {pt60 = pt; eta60 = eta; n60 = n;}
		else if(deltaeta< 1.0&&deltaeta60< 1.0) { if(pt>pt60) {pt60 = pt; eta60 = eta; n60 = n;} }  
		else if(deltaeta>=1.0&&deltaeta60>=1.0){
		if(     deltaeta< 1.2&&deltaeta60>=1.2) {pt60 = pt; eta60 = eta; n60 = n;}
		else if(deltaeta< 1.2&&deltaeta60< 1.2) { if(pt>pt60) {pt60 = pt; eta60 = eta; n60 = n;} }  
		else if(deltaeta>=1.2&&deltaeta60>=1.2){
		if(     deltaeta< 1.4&&deltaeta60>=1.4) {pt60 = pt; eta60 = eta; n60 = n;}
		else if(deltaeta< 1.4&&deltaeta60< 1.4) { if(pt>pt60) {pt60 = pt; eta60 = eta; n60 = n;} }  
		else if(deltaeta>=1.4&&deltaeta60>=1.4){
		if(     deltaeta< 1.6&&deltaeta60>=1.6) {pt60 = pt; eta60 = eta; n60 = n;}
		else if(deltaeta< 1.6&&deltaeta60< 1.6) { if(pt>pt60) {pt60 = pt; eta60 = eta; n60 = n;} }  }}}}}}}
		if(pt<80) continue;
		float deltaeta80 = fabs(1.5-eta80);
		if(     deltaeta< 0.2&&deltaeta80>=0.2) {pt80 = pt; eta80 = eta; n80 = n;}
		else if(deltaeta< 0.2&&deltaeta80< 0.2) { if(pt>pt80) {pt80 = pt; eta80 = eta; n80 = n;} }
		else if(deltaeta>=0.2&&deltaeta80>=0.2){
		if(     deltaeta< 0.4&&deltaeta80>=0.4) {pt80 = pt; eta80 = eta; n80 = n;}
		else if(deltaeta< 0.4&&deltaeta80< 0.4) { if(pt>pt80) {pt80 = pt; eta80 = eta; n80 = n;} }  
		else if(deltaeta>=0.4&&deltaeta80>=0.4){
		if(     deltaeta< 0.6&&deltaeta80>=0.6) {pt80 = pt; eta80 = eta; n80 = n;}
		else if(deltaeta< 0.6&&deltaeta80< 0.6) { if(pt>pt80) {pt80 = pt; eta80 = eta; n80 = n;} }  
		else if(deltaeta>=0.6&&deltaeta80>=0.6){
		if(     deltaeta< 0.8&&deltaeta80>=0.8) {pt80 = pt; eta80 = eta; n80 = n;}
		else if(deltaeta< 0.8&&deltaeta80< 0.8) { if(pt>pt80) {pt80 = pt; eta80 = eta; n80 = n;} }  
		else if(deltaeta>=0.8&&deltaeta80>=0.8){
		if(     deltaeta< 1.0&&deltaeta80>=1.0) {pt80 = pt; eta80 = eta; n80 = n;}
		else if(deltaeta< 1.0&&deltaeta80< 1.0) { if(pt>pt80) {pt80 = pt; eta80 = eta; n80 = n;} }  
		else if(deltaeta>=1.0&&deltaeta80>=1.0){
		if(     deltaeta< 1.2&&deltaeta80>=1.2) {pt80 = pt; eta80 = eta; n80 = n;}
		else if(deltaeta< 1.2&&deltaeta80< 1.2) { if(pt>pt80) {pt80 = pt; eta80 = eta; n80 = n;} }  
		else if(deltaeta>=1.2&&deltaeta80>=1.2){
		if(     deltaeta< 1.4&&deltaeta80>=1.4) {pt80 = pt; eta80 = eta; n80 = n;}
		else if(deltaeta< 1.4&&deltaeta80< 1.4) { if(pt>pt80) {pt80 = pt; eta80 = eta; n80 = n;} }  
		else if(deltaeta>=1.4&&deltaeta80>=1.4){
		if(     deltaeta< 1.6&&deltaeta80>=1.6) {pt80 = pt; eta80 = eta; n80 = n;}
		else if(deltaeta< 1.6&&deltaeta80< 1.6) { if(pt>pt80) {pt80 = pt; eta80 = eta; n80 = n;} }  }}}}}}}
		if(pt<100) continue;
		float deltaeta100 = fabs(1.5-eta100);
		//cout << "pt " << pt << " eta " << eta << " deltaeta " << deltaeta << " deltaeta100 " << deltaeta100 << " pt100 " << pt100 << endl;
		if(     deltaeta< 0.2&&deltaeta100>=0.2) {pt100 = pt; eta100 = eta; n100 = n;}
		else if(deltaeta< 0.2&&deltaeta100< 0.2) { if(pt>pt100) {pt100 = pt; eta100 = eta; n100 = n;} }
		else if(deltaeta>=0.2&&deltaeta100>=0.2){
		if(     deltaeta< 0.4&&deltaeta100>=0.4) {pt100 = pt; eta100 = eta; n100 = n;}
		else if(deltaeta< 0.4&&deltaeta100< 0.4) { if(pt>pt100) {pt100 = pt; eta100 = eta; n100 = n;} }  
		else if(deltaeta>=0.4&&deltaeta100>=0.4){
		if(     deltaeta< 0.6&&deltaeta100>=0.6) {pt100 = pt; eta100 = eta; n100 = n;}
		else if(deltaeta< 0.6&&deltaeta100< 0.6) { if(pt>pt100) {pt100 = pt; eta100 = eta; n100 = n;} }  
		else if(deltaeta>=0.6&&deltaeta100>=0.6){
		if(     deltaeta< 0.8&&deltaeta100>=0.8) {pt100 = pt; eta100 = eta; n100 = n;}
		else if(deltaeta< 0.8&&deltaeta100< 0.8) { if(pt>pt100) {pt100 = pt; eta100 = eta; n100 = n;} }  
		else if(deltaeta>=0.8&&deltaeta100>=0.8){
		if(     deltaeta< 1.0&&deltaeta100>=1.0) {pt100 = pt; eta100 = eta; n100 = n;}
		else if(deltaeta< 1.0&&deltaeta100< 1.0) { if(pt>pt100) {pt100 = pt; eta100 = eta; n100 = n;} }  
		else if(deltaeta>=1.0&&deltaeta100>=1.0){
		if(     deltaeta< 1.2&&deltaeta100>=1.2) {pt100 = pt; eta100 = eta; n100 = n;}
		else if(deltaeta< 1.2&&deltaeta100< 1.2) { if(pt>pt100) {pt100 = pt; eta100 = eta; n100 = n;} }  
		else if(deltaeta>=1.2&&deltaeta100>=1.2){
		if(     deltaeta< 1.4&&deltaeta100>=1.4) {pt100 = pt; eta100 = eta; n100 = n;}
		else if(deltaeta< 1.4&&deltaeta100< 1.4) { if(pt>pt100) {pt100 = pt; eta100 = eta; n100 = n;} }  
		else if(deltaeta>=1.4&&deltaeta100>=1.4){
		if(     deltaeta< 1.6&&deltaeta100>=1.6) {pt100 = pt; eta100 = eta; n100 = n;}
		else if(deltaeta< 1.6&&deltaeta100< 1.6) { if(pt>pt100) {pt100 = pt; eta100 = eta; n100 = n;} }  }}}}}}}
		if(pt<150) continue;
		float deltaeta150 = fabs(1.5-eta150);
		if(     deltaeta< 0.2&&deltaeta150>=0.2) {pt150 = pt; eta150 = eta; n150 = n;}
		else if(deltaeta< 0.2&&deltaeta150< 0.2) { if(pt>pt150) {pt150 = pt; eta150 = eta; n150 = n;} }
		else if(deltaeta>=0.2&&deltaeta150>=0.2){
		if(     deltaeta< 0.4&&deltaeta150>=0.4) {pt150 = pt; eta150 = eta; n150 = n;}
		else if(deltaeta< 0.4&&deltaeta150< 0.4) { if(pt>pt150) {pt150 = pt; eta150 = eta; n150 = n;} }  
		else if(deltaeta>=0.4&&deltaeta150>=0.4){
		if(     deltaeta< 0.6&&deltaeta150>=0.6) {pt150 = pt; eta150 = eta; n150 = n;}
		else if(deltaeta< 0.6&&deltaeta150< 0.6) { if(pt>pt150) {pt150 = pt; eta150 = eta; n150 = n;} }  
		else if(deltaeta>=0.6&&deltaeta150>=0.6){
		if(     deltaeta< 0.8&&deltaeta150>=0.8) {pt150 = pt; eta150 = eta; n150 = n;}
		else if(deltaeta< 0.8&&deltaeta150< 0.8) { if(pt>pt150) {pt150 = pt; eta150 = eta; n150 = n;} }  
		else if(deltaeta>=0.8&&deltaeta150>=0.8){
		if(     deltaeta< 1.0&&deltaeta150>=1.0) {pt150 = pt; eta150 = eta; n150 = n;}
		else if(deltaeta< 1.0&&deltaeta150< 1.0) { if(pt>pt150) {pt150 = pt; eta150 = eta; n150 = n;} }  
		else if(deltaeta>=1.0&&deltaeta150>=1.0){
		if(     deltaeta< 1.2&&deltaeta150>=1.2) {pt150 = pt; eta150 = eta; n150 = n;}
		else if(deltaeta< 1.2&&deltaeta150< 1.2) { if(pt>pt150) {pt150 = pt; eta150 = eta; n150 = n;} }  
		else if(deltaeta>=1.2&&deltaeta150>=1.2){
		if(     deltaeta< 1.4&&deltaeta150>=1.4) {pt150 = pt; eta150 = eta; n150 = n;}
		else if(deltaeta< 1.4&&deltaeta150< 1.4) { if(pt>pt150) {pt150 = pt; eta150 = eta; n150 = n;} }  
		else if(deltaeta>=1.4&&deltaeta150>=1.4){
		if(     deltaeta< 1.6&&deltaeta150>=1.6) {pt150 = pt; eta150 = eta; n150 = n;}
		else if(deltaeta< 1.6&&deltaeta150< 1.6) { if(pt>pt150) {pt150 = pt; eta150 = eta; n150 = n;} }  }}}}}}}
		if(pt<200) continue;
		float deltaeta200 = fabs(1.5-eta200);
		if(     deltaeta< 0.2&&deltaeta200>=0.2) {pt200 = pt; eta200 = eta; n200 = n;}
		else if(deltaeta< 0.2&&deltaeta200< 0.2) { if(pt>pt200) {pt200 = pt; eta200 = eta; n200 = n;} }
		else if(deltaeta>=0.2&&deltaeta200>=0.2){
		if(     deltaeta< 0.4&&deltaeta200>=0.4) {pt200 = pt; eta200 = eta; n200 = n;}
		else if(deltaeta< 0.4&&deltaeta200< 0.4) { if(pt>pt200) {pt200 = pt; eta200 = eta; n200 = n;} }  
		else if(deltaeta>=0.4&&deltaeta200>=0.4){
		if(     deltaeta< 0.6&&deltaeta200>=0.6) {pt200 = pt; eta200 = eta; n200 = n;}
		else if(deltaeta< 0.6&&deltaeta200< 0.6) { if(pt>pt200) {pt200 = pt; eta200 = eta; n200 = n;} }  
		else if(deltaeta>=0.6&&deltaeta200>=0.6){
		if(     deltaeta< 0.8&&deltaeta200>=0.8) {pt200 = pt; eta200 = eta; n200 = n;}
		else if(deltaeta< 0.8&&deltaeta200< 0.8) { if(pt>pt200) {pt200 = pt; eta200 = eta; n200 = n;} }  
		else if(deltaeta>=0.8&&deltaeta200>=0.8){
		if(     deltaeta< 1.0&&deltaeta200>=1.0) {pt200 = pt; eta200 = eta; n200 = n;}
		else if(deltaeta< 1.0&&deltaeta200< 1.0) { if(pt>pt200) {pt200 = pt; eta200 = eta; n200 = n;} }  
		else if(deltaeta>=1.0&&deltaeta200>=1.0){
		if(     deltaeta< 1.2&&deltaeta200>=1.2) {pt200 = pt; eta200 = eta; n200 = n;}
		else if(deltaeta< 1.2&&deltaeta200< 1.2) { if(pt>pt200) {pt200 = pt; eta200 = eta; n200 = n;} }  
		else if(deltaeta>=1.2&&deltaeta200>=1.2){
		if(     deltaeta< 1.4&&deltaeta200>=1.4) {pt200 = pt; eta200 = eta; n200 = n;}
		else if(deltaeta< 1.4&&deltaeta200< 1.4) { if(pt>pt200) {pt200 = pt; eta200 = eta; n200 = n;} }  
		else if(deltaeta>=1.4&&deltaeta200>=1.4){
		if(     deltaeta< 1.6&&deltaeta200>=1.6) {pt200 = pt; eta200 = eta; n200 = n;}
		else if(deltaeta< 1.6&&deltaeta200< 1.6) { if(pt>pt200) {pt200 = pt; eta200 = eta; n200 = n;} }  }}}}}}}
	}
	if(n40<0) cout << "ERROR THIS SHOULD NOT HAPPEN " << __LINE__ << endl;
	if(n60<0) cout << "ERROR THIS SHOULD NOT HAPPEN " << __LINE__ << endl;
	if(n80<0) cout << "ERROR THIS SHOULD NOT HAPPEN " << __LINE__ << endl;
	if(n100<0)cout << "ERROR THIS SHOULD NOT HAPPEN " << __LINE__ << endl;
	//create the map by filling histograms pass and all
	if(n40 >=0) SJetPt40Etaall ->Fill(fMT2tree->jet[n40 ].lv.Pt()*scale, fabs(fMT2tree->jet[n40 ].lv.Eta()));
	if(n60 >=0) SJetPt60Etaall ->Fill(fMT2tree->jet[n60 ].lv.Pt()*scale, fabs(fMT2tree->jet[n60 ].lv.Eta()));
	if(n80 >=0) SJetPt80Etaall ->Fill(fMT2tree->jet[n80 ].lv.Pt()*scale, fabs(fMT2tree->jet[n80 ].lv.Eta()));
	if(n100>=0) SJetPt100Etaall->Fill(fMT2tree->jet[n100].lv.Pt()*scale, fabs(fMT2tree->jet[n100].lv.Eta()));
	if(n150>=0) SJetPt150Etaall->Fill(fMT2tree->jet[n150].lv.Pt()*scale, fabs(fMT2tree->jet[n150].lv.Eta()));
	else        SJetPt150Etaall->Fill(40.                        , 0.0                            );
	if(n200>=0) SJetPt200Etaall->Fill(fMT2tree->jet[n200].lv.Pt()*scale, fabs(fMT2tree->jet[n200].lv.Eta()));
	else        SJetPt200Etaall->Fill(40.                        , 0.0                            );
	if(HT){
		if(n40 >=0) SJetPt40EtaallHT ->Fill(fMT2tree->jet[n40 ].lv.Pt()*scale, fabs(fMT2tree->jet[n40 ].lv.Eta()));
		if(n60 >=0) SJetPt60EtaallHT ->Fill(fMT2tree->jet[n60 ].lv.Pt()*scale, fabs(fMT2tree->jet[n60 ].lv.Eta()));
		if(n80 >=0) SJetPt80EtaallHT ->Fill(fMT2tree->jet[n80 ].lv.Pt()*scale, fabs(fMT2tree->jet[n80 ].lv.Eta()));
		if(n100>=0) SJetPt100EtaallHT->Fill(fMT2tree->jet[n100].lv.Pt()*scale, fabs(fMT2tree->jet[n100].lv.Eta()));
		if(n150>=0) SJetPt150EtaallHT->Fill(fMT2tree->jet[n150].lv.Pt()*scale, fabs(fMT2tree->jet[n150].lv.Eta()));
		else        SJetPt150EtaallHT->Fill(40.                        , 0.0                            );
		if(n200>=0) SJetPt200EtaallHT->Fill(fMT2tree->jet[n200].lv.Pt()*scale, fabs(fMT2tree->jet[n200].lv.Eta()));
		else        SJetPt200EtaallHT->Fill(40.                        , 0.0                            );
	} if(MET){
		if(n40 >=0) SJetPt40EtaallMET ->Fill(fMT2tree->jet[n40 ].lv.Pt()*scale, fabs(fMT2tree->jet[n40 ].lv.Eta()));
		if(n60 >=0) SJetPt60EtaallMET ->Fill(fMT2tree->jet[n60 ].lv.Pt()*scale, fabs(fMT2tree->jet[n60 ].lv.Eta()));
		if(n80 >=0) SJetPt80EtaallMET ->Fill(fMT2tree->jet[n80 ].lv.Pt()*scale, fabs(fMT2tree->jet[n80 ].lv.Eta()));
		if(n100>=0) SJetPt100EtaallMET->Fill(fMT2tree->jet[n100].lv.Pt()*scale, fabs(fMT2tree->jet[n100].lv.Eta()));
		if(n150>=0) SJetPt150EtaallMET->Fill(fMT2tree->jet[n150].lv.Pt()*scale, fabs(fMT2tree->jet[n150].lv.Eta()));
		else        SJetPt150EtaallMET->Fill(40.                        , 0.0                            );
		if(n200>=0) SJetPt200EtaallMET->Fill(fMT2tree->jet[n200].lv.Pt()*scale, fabs(fMT2tree->jet[n200].lv.Eta()));
		else        SJetPt200EtaallMET->Fill(40.                        , 0.0                            );
	}
	if(!tagged){
		if(n40 >=0) SJetPt40Etapass ->Fill(fMT2tree->jet[n40 ].lv.Pt()*scale, fabs(fMT2tree->jet[n40 ].lv.Eta()));
		if(n60 >=0) SJetPt60Etapass ->Fill(fMT2tree->jet[n60 ].lv.Pt()*scale, fabs(fMT2tree->jet[n60 ].lv.Eta()));
		if(n80 >=0) SJetPt80Etapass ->Fill(fMT2tree->jet[n80 ].lv.Pt()*scale, fabs(fMT2tree->jet[n80 ].lv.Eta()));
		if(n100>=0) SJetPt100Etapass->Fill(fMT2tree->jet[n100].lv.Pt()*scale, fabs(fMT2tree->jet[n100].lv.Eta()));
		if(n150>=0) SJetPt150Etapass->Fill(fMT2tree->jet[n150].lv.Pt()*scale, fabs(fMT2tree->jet[n150].lv.Eta()));
		else        SJetPt150Etapass->Fill(40.                        , 0.0                            );
		if(n200>=0) SJetPt200Etapass->Fill(fMT2tree->jet[n200].lv.Pt()*scale, fabs(fMT2tree->jet[n200].lv.Eta()));
		else        SJetPt200Etapass->Fill(40.                        , 0.0                            );
		if(HT){
			if(n40 >=0) SJetPt40EtapassHT ->Fill(fMT2tree->jet[n40 ].lv.Pt()*scale, fabs(fMT2tree->jet[n40 ].lv.Eta()));
			if(n60 >=0) SJetPt60EtapassHT ->Fill(fMT2tree->jet[n60 ].lv.Pt()*scale, fabs(fMT2tree->jet[n60 ].lv.Eta()));
			if(n80 >=0) SJetPt80EtapassHT ->Fill(fMT2tree->jet[n80 ].lv.Pt()*scale, fabs(fMT2tree->jet[n80 ].lv.Eta()));
			if(n100>=0) SJetPt100EtapassHT->Fill(fMT2tree->jet[n100].lv.Pt()*scale, fabs(fMT2tree->jet[n100].lv.Eta()));
			if(n150>=0) SJetPt150EtapassHT->Fill(fMT2tree->jet[n150].lv.Pt()*scale, fabs(fMT2tree->jet[n150].lv.Eta()));
			else        SJetPt150EtapassHT->Fill(40.                        , 0.0                            );
			if(n200>=0) SJetPt200EtapassHT->Fill(fMT2tree->jet[n200].lv.Pt()*scale, fabs(fMT2tree->jet[n200].lv.Eta()));
			else        SJetPt200EtapassHT->Fill(40.                        , 0.0                            );
		} if(MET){
			if(n40 >=0) SJetPt40EtapassMET ->Fill(fMT2tree->jet[n40 ].lv.Pt()*scale, fabs(fMT2tree->jet[n40 ].lv.Eta()));
			if(n60 >=0) SJetPt60EtapassMET ->Fill(fMT2tree->jet[n60 ].lv.Pt()*scale, fabs(fMT2tree->jet[n60 ].lv.Eta()));
			if(n80 >=0) SJetPt80EtapassMET ->Fill(fMT2tree->jet[n80 ].lv.Pt()*scale, fabs(fMT2tree->jet[n80 ].lv.Eta()));
			if(n100>=0) SJetPt100EtapassMET->Fill(fMT2tree->jet[n100].lv.Pt()*scale, fabs(fMT2tree->jet[n100].lv.Eta()));
			if(n150>=0) SJetPt150EtapassMET->Fill(fMT2tree->jet[n150].lv.Pt()*scale, fabs(fMT2tree->jet[n150].lv.Eta()));
			else        SJetPt150EtapassMET->Fill(40.                        , 0.0                            );
			if(n200>=0) SJetPt200EtapassMET->Fill(fMT2tree->jet[n200].lv.Pt()*scale, fabs(fMT2tree->jet[n200].lv.Eta()));
			else        SJetPt200EtapassMET->Fill(40.                        , 0.0                            );
		}
	}
	if(n40 >=0) SJetPt40all ->Fill(fMT2tree->jet[n40 ].lv.Pt()*scale);
	if(n60 >=0) SJetPt60all ->Fill(fMT2tree->jet[n60 ].lv.Pt()*scale);
	if(n80 >=0) SJetPt80all ->Fill(fMT2tree->jet[n80 ].lv.Pt()*scale);
	if(n100>=0) SJetPt100all->Fill(fMT2tree->jet[n100].lv.Pt()*scale);
	if(n150>=0) SJetPt150all->Fill(fMT2tree->jet[n150].lv.Pt()*scale);
	else        SJetPt150all->Fill(40.                        );
	if(n200>=0) SJetPt200all->Fill(fMT2tree->jet[n200].lv.Pt()*scale);
	else        SJetPt200all->Fill(40.                        );
	if(HT){
		if(n40 >=0) SJetPt40allHT ->Fill(fMT2tree->jet[n40 ].lv.Pt()*scale);
		if(n60 >=0) SJetPt60allHT ->Fill(fMT2tree->jet[n60 ].lv.Pt()*scale);
		if(n80 >=0) SJetPt80allHT ->Fill(fMT2tree->jet[n80 ].lv.Pt()*scale);
		if(n100>=0) SJetPt100allHT->Fill(fMT2tree->jet[n100].lv.Pt()*scale);
		if(n150>=0) SJetPt150allHT->Fill(fMT2tree->jet[n150].lv.Pt()*scale);
		else        SJetPt150allHT->Fill(40.                        );
		if(n200>=0) SJetPt200allHT->Fill(fMT2tree->jet[n200].lv.Pt()*scale);
		else        SJetPt200allHT->Fill(40.                        );
	} if(MET){
		if(n40 >=0) SJetPt40allMET ->Fill(fMT2tree->jet[n40 ].lv.Pt()*scale);
		if(n60 >=0) SJetPt60allMET ->Fill(fMT2tree->jet[n60 ].lv.Pt()*scale);
		if(n80 >=0) SJetPt80allMET ->Fill(fMT2tree->jet[n80 ].lv.Pt()*scale);
		if(n100>=0) SJetPt100allMET->Fill(fMT2tree->jet[n100].lv.Pt()*scale);
		if(n150>=0) SJetPt150allMET->Fill(fMT2tree->jet[n150].lv.Pt()*scale);
		else        SJetPt150allMET->Fill(40.                        );
		if(n200>=0) SJetPt200allMET->Fill(fMT2tree->jet[n200].lv.Pt()*scale);
		else        SJetPt200allMET->Fill(40.                        );
	}
	if(!tagged){
		if(n40 >=0) SJetPt40pass ->Fill(fMT2tree->jet[n40 ].lv.Pt()*scale);
		if(n60 >=0) SJetPt60pass ->Fill(fMT2tree->jet[n60 ].lv.Pt()*scale);
		if(n80 >=0) SJetPt80pass ->Fill(fMT2tree->jet[n80 ].lv.Pt()*scale);
		if(n100>=0) SJetPt100pass->Fill(fMT2tree->jet[n100].lv.Pt()*scale);
		if(n150>=0) SJetPt150pass->Fill(fMT2tree->jet[n150].lv.Pt()*scale);
		else        SJetPt150pass->Fill(40.                        );
		if(n200>=0) SJetPt200pass->Fill(fMT2tree->jet[n200].lv.Pt()*scale);
		else        SJetPt200pass->Fill(40.                        );
		if(HT){
			if(n40 >=0) SJetPt40passHT ->Fill(fMT2tree->jet[n40 ].lv.Pt()*scale);
			if(n60 >=0) SJetPt60passHT ->Fill(fMT2tree->jet[n60 ].lv.Pt()*scale);
			if(n80 >=0) SJetPt80passHT ->Fill(fMT2tree->jet[n80 ].lv.Pt()*scale);
			if(n100>=0) SJetPt100passHT->Fill(fMT2tree->jet[n100].lv.Pt()*scale);
			if(n150>=0) SJetPt150passHT->Fill(fMT2tree->jet[n150].lv.Pt()*scale);
			else        SJetPt150passHT->Fill(40.                        );
			if(n200>=0) SJetPt200passHT->Fill(fMT2tree->jet[n200].lv.Pt()*scale);
			else        SJetPt200passHT->Fill(40.                        );
		} if(MET){
			if(n40 >=0) SJetPt40passMET ->Fill(fMT2tree->jet[n40 ].lv.Pt()*scale);
			if(n60 >=0) SJetPt60passMET ->Fill(fMT2tree->jet[n60 ].lv.Pt()*scale);
			if(n80 >=0) SJetPt80passMET ->Fill(fMT2tree->jet[n80 ].lv.Pt()*scale);
			if(n100>=0) SJetPt100passMET->Fill(fMT2tree->jet[n100].lv.Pt()*scale);
			if(n150>=0) SJetPt150passMET->Fill(fMT2tree->jet[n150].lv.Pt()*scale);
			else        SJetPt150passMET->Fill(40.                        );
			if(n200>=0) SJetPt200passMET->Fill(fMT2tree->jet[n200].lv.Pt()*scale);
			else        SJetPt200passMET->Fill(40.                        );
		}
	}
	if(n40 >=0) SJet40Etaall ->Fill(fabs(fMT2tree->jet[n40 ].lv.Eta()));
	if(n60 >=0) SJet60Etaall ->Fill(fabs(fMT2tree->jet[n60 ].lv.Eta()));
	if(n80 >=0) SJet80Etaall ->Fill(fabs(fMT2tree->jet[n80 ].lv.Eta()));
	if(n100>=0) SJet100Etaall->Fill(fabs(fMT2tree->jet[n100].lv.Eta()));
	if(n150>=0) SJet150Etaall->Fill(fabs(fMT2tree->jet[n150].lv.Eta()));
	else        SJet150Etaall->Fill(0.0                            );
	if(n200>=0) SJet200Etaall->Fill(fabs(fMT2tree->jet[n200].lv.Eta()));
	else        SJet200Etaall->Fill(0.0                            );
	if(HT){
		if(n40 >=0) SJet40EtaallHT ->Fill(fabs(fMT2tree->jet[n40 ].lv.Eta()));
		if(n60 >=0) SJet60EtaallHT ->Fill(fabs(fMT2tree->jet[n60 ].lv.Eta()));
		if(n80 >=0) SJet80EtaallHT ->Fill(fabs(fMT2tree->jet[n80 ].lv.Eta()));
		if(n100>=0) SJet100EtaallHT->Fill(fabs(fMT2tree->jet[n100].lv.Eta()));
		if(n150>=0) SJet150EtaallHT->Fill(fabs(fMT2tree->jet[n150].lv.Eta()));
		else        SJet150EtaallHT->Fill(0.0                            );
		if(n200>=0) SJet200EtaallHT->Fill(fabs(fMT2tree->jet[n200].lv.Eta()));
		else        SJet200EtaallHT->Fill(0.0                            );
	} if(MET){
		if(n40 >=0) SJet40EtaallMET ->Fill(fabs(fMT2tree->jet[n40 ].lv.Eta()));
		if(n60 >=0) SJet60EtaallMET ->Fill(fabs(fMT2tree->jet[n60 ].lv.Eta()));
		if(n80 >=0) SJet80EtaallMET ->Fill(fabs(fMT2tree->jet[n80 ].lv.Eta()));
		if(n100>=0) SJet100EtaallMET->Fill(fabs(fMT2tree->jet[n100].lv.Eta()));
		if(n150>=0) SJet150EtaallMET->Fill(fabs(fMT2tree->jet[n150].lv.Eta()));
		else        SJet150EtaallMET->Fill(0.0                            );
		if(n200>=0) SJet200EtaallMET->Fill(fabs(fMT2tree->jet[n200].lv.Eta()));
		else        SJet200EtaallMET->Fill(0.0                            );
	}
	if(!tagged){
		if(n40 >=0) SJet40Etapass ->Fill(fabs(fMT2tree->jet[n40 ].lv.Eta()));
		if(n60 >=0) SJet60Etapass ->Fill(fabs(fMT2tree->jet[n60 ].lv.Eta()));
		if(n80 >=0) SJet80Etapass ->Fill(fabs(fMT2tree->jet[n80 ].lv.Eta()));
		if(n100>=0) SJet100Etapass->Fill(fabs(fMT2tree->jet[n100].lv.Eta()));
		if(n150>=0) SJet150Etapass->Fill(fabs(fMT2tree->jet[n150].lv.Eta()));
		else        SJet150Etapass->Fill(0.0                            );
		if(n200>=0) SJet200Etapass->Fill(fabs(fMT2tree->jet[n200].lv.Eta()));
		else        SJet200Etapass->Fill(0.0                            );
		if(HT){
			if(n40 >=0) SJet40EtapassHT ->Fill(fabs(fMT2tree->jet[n40 ].lv.Eta()));
			if(n60 >=0) SJet60EtapassHT ->Fill(fabs(fMT2tree->jet[n60 ].lv.Eta()));
			if(n80 >=0) SJet80EtapassHT ->Fill(fabs(fMT2tree->jet[n80 ].lv.Eta()));
			if(n100>=0) SJet100EtapassHT->Fill(fabs(fMT2tree->jet[n100].lv.Eta()));
			if(n150>=0) SJet150EtapassHT->Fill(fabs(fMT2tree->jet[n150].lv.Eta()));
			else        SJet150EtapassHT->Fill(0.0                            );
			if(n200>=0) SJet200EtapassHT->Fill(fabs(fMT2tree->jet[n200].lv.Eta()));
			else        SJet200EtapassHT->Fill(0.0                            );
		} if(MET){
			if(n40 >=0) SJet40EtapassMET ->Fill(fabs(fMT2tree->jet[n40 ].lv.Eta()));
			if(n60 >=0) SJet60EtapassMET ->Fill(fabs(fMT2tree->jet[n60 ].lv.Eta()));
			if(n80 >=0) SJet80EtapassMET ->Fill(fabs(fMT2tree->jet[n80 ].lv.Eta()));
			if(n100>=0) SJet100EtapassMET->Fill(fabs(fMT2tree->jet[n100].lv.Eta()));
			if(n150>=0) SJet150EtapassMET->Fill(fabs(fMT2tree->jet[n150].lv.Eta()));
			else        SJet150EtapassMET->Fill(0.0                            );
			if(n200>=0) SJet200EtapassMET->Fill(fabs(fMT2tree->jet[n200].lv.Eta()));
			else        SJet200EtapassMET->Fill(0.0                            );
		}
	}

}
//get all efficiencies
TH1D *MT2alleff = (TH1D*)MT2all->Clone("MT2alleff");
TH1D *MT2passeff = (TH1D*)MT2pass->Clone("MT2passeff");
TH1D *MT2allHTeff = (TH1D*)MT2allHT->Clone("MT2allHTeff");
TH1D *MT2passHTeff = (TH1D*)MT2passHT->Clone("MT2passHTeff");
TH1D *MT2allMETeff = (TH1D*)MT2allMET->Clone("MT2allMETeff");
TH1D *MT2passMETeff = (TH1D*)MT2passMET->Clone("MT2passMETeff");
TEfficiency *MT2eff = new TEfficiency((*MT2passeff),(*MT2alleff));
MT2eff->SetName("MT2eff");
TEfficiency *MT2effHT = new TEfficiency((*MT2passHTeff),(*MT2allHTeff));
MT2effHT->SetName("MT2effHT");
TEfficiency *MT2effMET = new TEfficiency((*MT2passMETeff),(*MT2allMETeff));
MT2effMET->SetName("MT2effMET");
TF1 *MT2fit = new TF1("MT2fit","[0]",150,1000);
TF1 *MT2fitHT = new TF1("MT2fitHT","[0]",150,1000);
TF1 *MT2fitMET = new TF1("MT2fitMET","[0]",200,1000);
if(fit) MT2eff->Fit(MT2fit,"R");
if(fit) MT2effHT->Fit(MT2fitHT,"R");
if(fit) MT2effMET->Fit(MT2fitMET,"R");
TH1D *HTalleff = (TH1D*)HTall->Clone("HTalleff");
TH1D *HTpasseff = (TH1D*)HTpass->Clone("HTpasseff");
TH1D *HTallHTeff = (TH1D*)HTallHT->Clone("HTallHTeff");
TH1D *HTpassHTeff = (TH1D*)HTpassHT->Clone("HTpassHTeff");
TH1D *HTallMETeff = (TH1D*)HTallMET->Clone("HTallMETeff");
TH1D *HTpassMETeff = (TH1D*)HTpassMET->Clone("HTpassMETeff");
TEfficiency *HTeff = new TEfficiency((*HTpasseff),(*HTalleff));
HTeff->SetName("HTeff");
TEfficiency *HTeffHT = new TEfficiency((*HTpassHTeff),(*HTallHTeff));
HTeffHT->SetName("HTeffHT");
TEfficiency *HTeffMET = new TEfficiency((*HTpassMETeff),(*HTallMETeff));
HTeffMET->SetName("HTeffMET");
TF1 *HTfit1 = new TF1("HTfit1","[0]+x*[1]",450,750);
TF1 *HTfit2 = new TF1("HTfit2","[0]",750,2050);
TF1 *HTfitHT = new TF1("HTfitHT","[0]", 750,2050);
TF1 *HTfitMET = new TF1("HTfitMET","[0]+x*[1]", 450,750);
if(fit) HTeff->Fit(HTfit1,"R");
if(fit) HTeff->Fit(HTfit2,"R");
if(fit) HTeffHT->Fit(HTfitHT,"R");
if(fit) HTeffMET->Fit(HTfitMET,"R");
TH1D *METalleff = (TH1D*)METall->Clone("METalleff");
TH1D *METpasseff = (TH1D*)METpass->Clone("METpasseff");
TH1D *METallHTeff = (TH1D*)METallHT->Clone("METallHTeff");
TH1D *METpassHTeff = (TH1D*)METpassHT->Clone("METpassHTeff");
TH1D *METallMETeff = (TH1D*)METallMET->Clone("METallMETeff");
TH1D *METpassMETeff = (TH1D*)METpassMET->Clone("METpassMETeff");
TEfficiency *METeff = new TEfficiency((*METpasseff),(*METalleff));
METeff->SetName("METeff");
TEfficiency *METeffHT = new TEfficiency((*METpassHTeff),(*METallHTeff));
METeffHT->SetName("METeffHT");
TEfficiency *METeffMET = new TEfficiency((*METpassMETeff),(*METallMETeff));
METeffMET->SetName("METeffMET");
TF1 *METfit = new TF1("METfit","[0]",150,1000);
TF1 *METfitHT = new TF1("METfitHT","[0]",150,1000);
TF1 *METfitMET = new TF1("METfitMET","[0]",200,1000);
if(fit) METeff->Fit(METfit,"R");
if(fit) METeffHT->Fit(METfitHT,"R");
if(fit) METeffMET->Fit(METfitMET,"R");
TH1D *NJetsalleff = (TH1D*)NJetsall->Clone("NJetsalleff");
TH1D *NJetspasseff = (TH1D*)NJetspass->Clone("NJetspasseff");
TH1D *NJetsallHTeff = (TH1D*)NJetsallHT->Clone("NJetsallHTeff");
TH1D *NJetspassHTeff = (TH1D*)NJetspassHT->Clone("NJetspassHTeff");
TH1D *NJetsallMETeff = (TH1D*)NJetsallMET->Clone("NJetsallMETeff");
TH1D *NJetspassMETeff = (TH1D*)NJetspassMET->Clone("NJetspassMETeff");
TEfficiency *NJetseff = new TEfficiency((*NJetspasseff),(*NJetsalleff));
NJetseff->SetName("NJetseff");
TEfficiency *NJetseffHT = new TEfficiency((*NJetspassHTeff),(*NJetsallHTeff));
NJetseffHT->SetName("NJetseffHT");
TEfficiency *NJetseffMET = new TEfficiency((*NJetspassMETeff),(*NJetsallMETeff));
NJetseffMET->SetName("NJetseffMET");
TF1 *NJetsfit = new TF1("NJetsfit","[0]",2,10);
TF1 *NJetsfitHT = new TF1("NJetsfitHT","[0]",2,10);
TF1 *NJetsfitMET = new TF1("NJetsfitMET","[0]",2,10);
if(fit) NJetseff->Fit(NJetsfit,"R");
if(fit) NJetseffHT->Fit(NJetsfitHT,"R");
if(fit) NJetseffMET->Fit(NJetsfitMET,"R");
TH1D *NBJetsalleff = (TH1D*)NBJetsall->Clone("NBJetsalleff");
TH1D *NBJetspasseff = (TH1D*)NBJetspass->Clone("NBJetspasseff");
TH1D *NBJetsallHTeff = (TH1D*)NBJetsallHT->Clone("NBJetsallHTeff");
TH1D *NBJetspassHTeff = (TH1D*)NBJetspassHT->Clone("NBJetspassHTeff");
TH1D *NBJetsallMETeff = (TH1D*)NBJetsallMET->Clone("NBJetsallMETeff");
TH1D *NBJetspassMETeff = (TH1D*)NBJetspassMET->Clone("NBJetspassMETeff");
TEfficiency *NBJetseff = new TEfficiency((*NBJetspasseff),(*NBJetsalleff));
NBJetseff->SetName("NBJetseff");
TEfficiency *NBJetseffHT = new TEfficiency((*NBJetspassHTeff),(*NBJetsallHTeff));
NBJetseffHT->SetName("NBJetseffHT");
TEfficiency *NBJetseffMET = new TEfficiency((*NBJetspassMETeff),(*NBJetsallMETeff));
NBJetseffMET->SetName("NBJetseffMET");
TF1 *NBJetsfit = new TF1("NBJetsfit","[0]",0,5);
TF1 *NBJetsfitHT = new TF1("NBJetsfitHT","[0]",0,5);
TF1 *NBJetsfitMET = new TF1("NBJetsfitMET","[0]",0,5);
if(fit) NBJetseff->Fit(NBJetsfit,"R");
if(fit) NBJetseffHT->Fit(NBJetsfitHT,"R");
if(fit) NBJetseffMET->Fit(NBJetsfitMET,"R");
TH1D *Jet1Etaalleff = (TH1D*)Jet1Etaall->Clone("Jet1Etaalleff");
TH1D *Jet1Etapasseff = (TH1D*)Jet1Etapass->Clone("Jet1Etapasseff");
TH1D *Jet1EtaallHTeff = (TH1D*)Jet1EtaallHT->Clone("Jet1EtaallHTeff");
TH1D *Jet1EtapassHTeff = (TH1D*)Jet1EtapassHT->Clone("Jet1EtapassHTeff");
TH1D *Jet1EtaallMETeff = (TH1D*)Jet1EtaallMET->Clone("Jet1EtaallMETeff");
TH1D *Jet1EtapassMETeff = (TH1D*)Jet1EtapassMET->Clone("Jet1EtapassMETeff");
TEfficiency *Jet1Etaeff = new TEfficiency((*Jet1Etapasseff),(*Jet1Etaalleff));
Jet1Etaeff->SetName("Jet1Etaeff");
TEfficiency *Jet1EtaeffHT = new TEfficiency((*Jet1EtapassHTeff),(*Jet1EtaallHTeff));
Jet1EtaeffHT->SetName("Jet1EtaeffHT");
TEfficiency *Jet1EtaeffMET = new TEfficiency((*Jet1EtapassMETeff),(*Jet1EtaallMETeff));
Jet1EtaeffMET->SetName("Jet1EtaeffMET");
TH1D *Jet2Etaalleff = (TH1D*)Jet2Etaall->Clone("Jet2Etaalleff");
TH1D *Jet2Etapasseff = (TH1D*)Jet2Etapass->Clone("Jet2Etapasseff");
TH1D *Jet2EtaallHTeff = (TH1D*)Jet2EtaallHT->Clone("Jet2EtaallHTeff");
TH1D *Jet2EtapassHTeff = (TH1D*)Jet2EtapassHT->Clone("Jet2EtapassHTeff");
TH1D *Jet2EtaallMETeff = (TH1D*)Jet2EtaallMET->Clone("Jet2EtaallMETeff");
TH1D *Jet2EtapassMETeff = (TH1D*)Jet2EtapassMET->Clone("Jet2EtapassMETeff");
TEfficiency *Jet2Etaeff = new TEfficiency((*Jet2Etapasseff),(*Jet2Etaalleff));
Jet2Etaeff->SetName("Jet2Etaeff");
TEfficiency *Jet2EtaeffHT = new TEfficiency((*Jet2EtapassHTeff),(*Jet2EtaallHTeff));
Jet2EtaeffHT->SetName("Jet2EtaeffHT");
TEfficiency *Jet2EtaeffMET = new TEfficiency((*Jet2EtapassMETeff),(*Jet2EtaallMETeff));
Jet2EtaeffMET->SetName("Jet2EtaeffMET");
TH1D *Jet3Etaalleff = (TH1D*)Jet3Etaall->Clone("Jet3Etaalleff");
TH1D *Jet3Etapasseff = (TH1D*)Jet3Etapass->Clone("Jet3Etapasseff");
TH1D *Jet3EtaallHTeff = (TH1D*)Jet3EtaallHT->Clone("Jet3EtaallHTeff");
TH1D *Jet3EtapassHTeff = (TH1D*)Jet3EtapassHT->Clone("Jet3EtapassHTeff");
TH1D *Jet3EtaallMETeff = (TH1D*)Jet3EtaallMET->Clone("Jet3EtaallMETeff");
TH1D *Jet3EtapassMETeff = (TH1D*)Jet3EtapassMET->Clone("Jet3EtapassMETeff");
TEfficiency *Jet3Etaeff = new TEfficiency((*Jet3Etapasseff),(*Jet3Etaalleff));
Jet3Etaeff->SetName("Jet3Etaeff");
TEfficiency *Jet3EtaeffHT = new TEfficiency((*Jet3EtapassHTeff),(*Jet3EtaallHTeff));
Jet3EtaeffHT->SetName("Jet3EtaeffHT");
TEfficiency *Jet3EtaeffMET = new TEfficiency((*Jet3EtapassMETeff),(*Jet3EtaallMETeff));
Jet3EtaeffMET->SetName("Jet3EtaeffMET");
TH1D *Jet4Etaalleff = (TH1D*)Jet4Etaall->Clone("Jet4Etaalleff");
TH1D *Jet4Etapasseff = (TH1D*)Jet4Etapass->Clone("Jet4Etapasseff");
TH1D *Jet4EtaallHTeff = (TH1D*)Jet4EtaallHT->Clone("Jet4EtaallHTeff");
TH1D *Jet4EtapassHTeff = (TH1D*)Jet4EtapassHT->Clone("Jet4EtapassHTeff");
TH1D *Jet4EtaallMETeff = (TH1D*)Jet4EtaallMET->Clone("Jet4EtaallMETeff");
TH1D *Jet4EtapassMETeff = (TH1D*)Jet4EtapassMET->Clone("Jet4EtapassMETeff");
TEfficiency *Jet4Etaeff = new TEfficiency((*Jet4Etapasseff),(*Jet4Etaalleff));
Jet4Etaeff->SetName("Jet4Etaeff");
TEfficiency *Jet4EtaeffHT = new TEfficiency((*Jet4EtapassHTeff),(*Jet4EtaallHTeff));
Jet4EtaeffHT->SetName("Jet4EtaeffHT");
TEfficiency *Jet4EtaeffMET = new TEfficiency((*Jet4EtapassMETeff),(*Jet4EtaallMETeff));
Jet4EtaeffMET->SetName("Jet4EtaeffMET");
TH2D *JetPtEtaalleff = (TH2D*)JetPtEtaall->Clone("JetPtEtaalleff");
TH2D *JetPtEtapasseff = (TH2D*)JetPtEtapass->Clone("JetPtEtapasseff");
TH2D *JetPtEtaallHTeff = (TH2D*)JetPtEtaallHT->Clone("JetPtEtaallHTeff");
TH2D *JetPtEtapassHTeff = (TH2D*)JetPtEtapassHT->Clone("JetPtEtapassHTeff");
TH2D *JetPtEtaallMETeff = (TH2D*)JetPtEtaallMET->Clone("JetPtEtaallMETeff");
TH2D *JetPtEtapassMETeff = (TH2D*)JetPtEtapassMET->Clone("JetPtEtapassMETeff");
TEfficiency *JetPtEtaeff = new TEfficiency((*JetPtEtapasseff),(*JetPtEtaalleff));
JetPtEtaeff->SetName("JetPtEtaeff");
TEfficiency *JetPtEtaeffHT = new TEfficiency((*JetPtEtapassHTeff),(*JetPtEtaallHTeff));
JetPtEtaeffHT->SetName("JetPtEtaeffHT");
TEfficiency *JetPtEtaeffMET = new TEfficiency((*JetPtEtapassMETeff),(*JetPtEtaallMETeff));
JetPtEtaeffMET->SetName("JetPtEtaeffMET");
TH1D *JetEtaalleff = (TH1D*)JetEtaall->Clone("JetEtaalleff");
TH1D *JetEtapasseff = (TH1D*)JetEtapass->Clone("JetEtapasseff");
TH1D *JetEtaallHTeff = (TH1D*)JetEtaallHT->Clone("JetEtaallHTeff");
TH1D *JetEtapassHTeff = (TH1D*)JetEtapassHT->Clone("JetEtapassHTeff");
TH1D *JetEtaallMETeff = (TH1D*)JetEtaallMET->Clone("JetEtaallMETeff");
TH1D *JetEtapassMETeff = (TH1D*)JetEtapassMET->Clone("JetEtapassMETeff");
TEfficiency *JetEtaeff = new TEfficiency((*JetEtapasseff),(*JetEtaalleff));
JetEtaeff->SetName("JetEtaeff");
TEfficiency *JetEtaeffHT = new TEfficiency((*JetEtapassHTeff),(*JetEtaallHTeff));
JetEtaeffHT->SetName("JetEtaeffHT");
TEfficiency *JetEtaeffMET = new TEfficiency((*JetEtapassMETeff),(*JetEtaallMETeff));
JetEtaeffMET->SetName("JetEtaeffMET");
TH1D *JetPtalleff = (TH1D*)JetPtall->Clone("JetPtalleff");
TH1D *JetPtpasseff = (TH1D*)JetPtpass->Clone("JetPtpasseff");
TH1D *JetPtallHTeff = (TH1D*)JetPtallHT->Clone("JetPtallHTeff");
TH1D *JetPtpassHTeff = (TH1D*)JetPtpassHT->Clone("JetPtpassHTeff");
TH1D *JetPtallMETeff = (TH1D*)JetPtallMET->Clone("JetPtallMETeff");
TH1D *JetPtpassMETeff = (TH1D*)JetPtpassMET->Clone("JetPtpassMETeff");
TEfficiency *JetPteff = new TEfficiency((*JetPtpasseff),(*JetPtalleff));
JetPteff->SetName("JetPteff");
TEfficiency *JetPteffHT = new TEfficiency((*JetPtpassHTeff),(*JetPtallHTeff));
JetPteffHT->SetName("JetPteffHT");
TEfficiency *JetPteffMET = new TEfficiency((*JetPtpassMETeff),(*JetPtallMETeff));
JetPteffMET->SetName("JetPteffMET");
TH2D *SJetPt40Etaalleff = (TH2D*)SJetPt40Etaall->Clone("SJetPt40Etaalleff");
TH2D *SJetPt40Etapasseff = (TH2D*)SJetPt40Etapass->Clone("SJetPt40Etapasseff");
TH2D *SJetPt40EtaallHTeff = (TH2D*)SJetPt40EtaallHT->Clone("SJetPt40EtaallHTeff");
TH2D *SJetPt40EtapassHTeff = (TH2D*)SJetPt40EtapassHT->Clone("SJetPt40EtapassHTeff");
TH2D *SJetPt40EtaallMETeff = (TH2D*)SJetPt40EtaallMET->Clone("SJetPt40EtaallMETeff");
TH2D *SJetPt40EtapassMETeff = (TH2D*)SJetPt40EtapassMET->Clone("SJetPt40EtapassMETeff");
TEfficiency *SJetPt40Etaeff = new TEfficiency((*SJetPt40Etapasseff),(*SJetPt40Etaalleff));
SJetPt40Etaeff->SetName("SJetPt40Etaeff");
TEfficiency *SJetPt40EtaeffHT = new TEfficiency((*SJetPt40EtapassHTeff),(*SJetPt40EtaallHTeff));
SJetPt40EtaeffHT->SetName("SJetPt40EtaeffHT");
TEfficiency *SJetPt40EtaeffMET = new TEfficiency((*SJetPt40EtapassMETeff),(*SJetPt40EtaallMETeff));
SJetPt40EtaeffMET->SetName("SJetPt40EtaeffMET");
TH2D *SJetPt60Etaalleff = (TH2D*)SJetPt60Etaall->Clone("SJetPt60Etaalleff");
TH2D *SJetPt60Etapasseff = (TH2D*)SJetPt60Etapass->Clone("SJetPt60Etapasseff");
TH2D *SJetPt60EtaallHTeff = (TH2D*)SJetPt60EtaallHT->Clone("SJetPt60EtaallHTeff");
TH2D *SJetPt60EtapassHTeff = (TH2D*)SJetPt60EtapassHT->Clone("SJetPt60EtapassHTeff");
TH2D *SJetPt60EtaallMETeff = (TH2D*)SJetPt60EtaallMET->Clone("SJetPt60EtaallMETeff");
TH2D *SJetPt60EtapassMETeff = (TH2D*)SJetPt60EtapassMET->Clone("SJetPt60EtapassMETeff");
TEfficiency *SJetPt60Etaeff = new TEfficiency((*SJetPt60Etapasseff),(*SJetPt60Etaalleff));
SJetPt60Etaeff->SetName("SJetPt60Etaeff");
TEfficiency *SJetPt60EtaeffHT = new TEfficiency((*SJetPt60EtapassHTeff),(*SJetPt60EtaallHTeff));
SJetPt60EtaeffHT->SetName("SJetPt60EtaeffHT");
TEfficiency *SJetPt60EtaeffMET = new TEfficiency((*SJetPt60EtapassMETeff),(*SJetPt60EtaallMETeff));
SJetPt60EtaeffMET->SetName("SJetPt60EtaeffMET");
TH2D *SJetPt80Etaalleff = (TH2D*)SJetPt80Etaall->Clone("SJetPt80Etaalleff");
TH2D *SJetPt80Etapasseff = (TH2D*)SJetPt80Etapass->Clone("SJetPt80Etapasseff");
TH2D *SJetPt80EtaallHTeff = (TH2D*)SJetPt80EtaallHT->Clone("SJetPt80EtaallHTeff");
TH2D *SJetPt80EtapassHTeff = (TH2D*)SJetPt80EtapassHT->Clone("SJetPt80EtapassHTeff");
TH2D *SJetPt80EtaallMETeff = (TH2D*)SJetPt80EtaallMET->Clone("SJetPt80EtaallMETeff");
TH2D *SJetPt80EtapassMETeff = (TH2D*)SJetPt80EtapassMET->Clone("SJetPt80EtapassMETeff");
TEfficiency *SJetPt80Etaeff = new TEfficiency((*SJetPt80Etapasseff),(*SJetPt80Etaalleff));
SJetPt80Etaeff->SetName("SJetPt80Etaeff");
TEfficiency *SJetPt80EtaeffHT = new TEfficiency((*SJetPt80EtapassHTeff),(*SJetPt80EtaallHTeff));
SJetPt80EtaeffHT->SetName("SJetPt80EtaeffHT");
TEfficiency *SJetPt80EtaeffMET = new TEfficiency((*SJetPt80EtapassMETeff),(*SJetPt80EtaallMETeff));
SJetPt80EtaeffMET->SetName("SJetPt80EtaeffMET");
TH2D *SJetPt100Etaalleff = (TH2D*)SJetPt100Etaall->Clone("SJetPt100Etaalleff");
TH2D *SJetPt100Etapasseff = (TH2D*)SJetPt100Etapass->Clone("SJetPt100Etapasseff");
TH2D *SJetPt100EtaallHTeff = (TH2D*)SJetPt100EtaallHT->Clone("SJetPt100EtaallHTeff");
TH2D *SJetPt100EtapassHTeff = (TH2D*)SJetPt100EtapassHT->Clone("SJetPt100EtapassHTeff");
TH2D *SJetPt100EtaallMETeff = (TH2D*)SJetPt100EtaallMET->Clone("SJetPt100EtaallMETeff");
TH2D *SJetPt100EtapassMETeff = (TH2D*)SJetPt100EtapassMET->Clone("SJetPt100EtapassMETeff");
TEfficiency *SJetPt100Etaeff = new TEfficiency((*SJetPt100Etapasseff),(*SJetPt100Etaalleff));
SJetPt100Etaeff->SetName("SJetPt100Etaeff");
TEfficiency *SJetPt100EtaeffHT = new TEfficiency((*SJetPt100EtapassHTeff),(*SJetPt100EtaallHTeff));
SJetPt100EtaeffHT->SetName("SJetPt100EtaeffHT");
TEfficiency *SJetPt100EtaeffMET = new TEfficiency((*SJetPt100EtapassMETeff),(*SJetPt100EtaallMETeff));
SJetPt100EtaeffMET->SetName("SJetPt100EtaeffMET");
TH2D *SJetPt150Etaalleff = (TH2D*)SJetPt150Etaall->Clone("SJetPt150Etaalleff");
TH2D *SJetPt150Etapasseff = (TH2D*)SJetPt150Etapass->Clone("SJetPt150Etapasseff");
TH2D *SJetPt150EtaallHTeff = (TH2D*)SJetPt150EtaallHT->Clone("SJetPt150EtaallHTeff");
TH2D *SJetPt150EtapassHTeff = (TH2D*)SJetPt150EtapassHT->Clone("SJetPt150EtapassHTeff");
TH2D *SJetPt150EtaallMETeff = (TH2D*)SJetPt150EtaallMET->Clone("SJetPt150EtaallMETeff");
TH2D *SJetPt150EtapassMETeff = (TH2D*)SJetPt150EtapassMET->Clone("SJetPt150EtapassMETeff");
TEfficiency *SJetPt150Etaeff = new TEfficiency((*SJetPt150Etapasseff),(*SJetPt150Etaalleff));
SJetPt150Etaeff->SetName("SJetPt150Etaeff");
TEfficiency *SJetPt150EtaeffHT = new TEfficiency((*SJetPt150EtapassHTeff),(*SJetPt150EtaallHTeff));
SJetPt150EtaeffHT->SetName("SJetPt150EtaeffHT");
TEfficiency *SJetPt150EtaeffMET = new TEfficiency((*SJetPt150EtapassMETeff),(*SJetPt150EtaallMETeff));
SJetPt150EtaeffMET->SetName("SJetPt150EtaeffMET");
TH2D *SJetPt200Etaalleff = (TH2D*)SJetPt200Etaall->Clone("SJetPt200Etaalleff");
TH2D *SJetPt200Etapasseff = (TH2D*)SJetPt200Etapass->Clone("SJetPt200Etapasseff");
TH2D *SJetPt200EtaallHTeff = (TH2D*)SJetPt200EtaallHT->Clone("SJetPt200EtaallHTeff");
TH2D *SJetPt200EtapassHTeff = (TH2D*)SJetPt200EtapassHT->Clone("SJetPt200EtapassHTeff");
TH2D *SJetPt200EtaallMETeff = (TH2D*)SJetPt200EtaallMET->Clone("SJetPt200EtaallMETeff");
TH2D *SJetPt200EtapassMETeff = (TH2D*)SJetPt200EtapassMET->Clone("SJetPt200EtapassMETeff");
TEfficiency *SJetPt200Etaeff = new TEfficiency((*SJetPt200Etapasseff),(*SJetPt200Etaalleff));
SJetPt200Etaeff->SetName("SJetPt200Etaeff");
TEfficiency *SJetPt200EtaeffHT = new TEfficiency((*SJetPt200EtapassHTeff),(*SJetPt200EtaallHTeff));
SJetPt200EtaeffHT->SetName("SJetPt200EtaeffHT");
TEfficiency *SJetPt200EtaeffMET = new TEfficiency((*SJetPt200EtapassMETeff),(*SJetPt200EtaallMETeff));
SJetPt200EtaeffMET->SetName("SJetPt200EtaeffMET");
TH1D *SJet40Etaalleff = (TH1D*)SJet40Etaall->Clone("SJet40Etaalleff");
TH1D *SJet40Etapasseff = (TH1D*)SJet40Etapass->Clone("SJet40Etapasseff");
TH1D *SJet40EtaallHTeff = (TH1D*)SJet40EtaallHT->Clone("SJet40EtaallHTeff");
TH1D *SJet40EtapassHTeff = (TH1D*)SJet40EtapassHT->Clone("SJet40EtapassHTeff");
TH1D *SJet40EtaallMETeff = (TH1D*)SJet40EtaallMET->Clone("SJet40EtaallMETeff");
TH1D *SJet40EtapassMETeff = (TH1D*)SJet40EtapassMET->Clone("SJet40EtapassMETeff");
TEfficiency *SJet40Etaeff = new TEfficiency((*SJet40Etapasseff),(*SJet40Etaalleff));
SJet40Etaeff->SetName("SJet40Etaeff");
TEfficiency *SJet40EtaeffHT = new TEfficiency((*SJet40EtapassHTeff),(*SJet40EtaallHTeff));
SJet40EtaeffHT->SetName("SJet40EtaeffHT");
TEfficiency *SJet40EtaeffMET = new TEfficiency((*SJet40EtapassMETeff),(*SJet40EtaallMETeff));
SJet40EtaeffMET->SetName("SJet40EtaeffMET");
TH1D *SJet60Etaalleff = (TH1D*)SJet60Etaall->Clone("SJet60Etaalleff");
TH1D *SJet60Etapasseff = (TH1D*)SJet60Etapass->Clone("SJet60Etapasseff");
TH1D *SJet60EtaallHTeff = (TH1D*)SJet60EtaallHT->Clone("SJet60EtaallHTeff");
TH1D *SJet60EtapassHTeff = (TH1D*)SJet60EtapassHT->Clone("SJet60EtapassHTeff");
TH1D *SJet60EtaallMETeff = (TH1D*)SJet60EtaallMET->Clone("SJet60EtaallMETeff");
TH1D *SJet60EtapassMETeff = (TH1D*)SJet60EtapassMET->Clone("SJet60EtapassMETeff");
TEfficiency *SJet60Etaeff = new TEfficiency((*SJet60Etapasseff),(*SJet60Etaalleff));
SJet60Etaeff->SetName("SJet60Etaeff");
TEfficiency *SJet60EtaeffHT = new TEfficiency((*SJet60EtapassHTeff),(*SJet60EtaallHTeff));
SJet60EtaeffHT->SetName("SJet60EtaeffHT");
TEfficiency *SJet60EtaeffMET = new TEfficiency((*SJet60EtapassMETeff),(*SJet60EtaallMETeff));
SJet60EtaeffMET->SetName("SJet60EtaeffMET");
TH1D *SJet80Etaalleff = (TH1D*)SJet80Etaall->Clone("SJet80Etaalleff");
TH1D *SJet80Etapasseff = (TH1D*)SJet80Etapass->Clone("SJet80Etapasseff");
TH1D *SJet80EtaallHTeff = (TH1D*)SJet80EtaallHT->Clone("SJet80EtaallHTeff");
TH1D *SJet80EtapassHTeff = (TH1D*)SJet80EtapassHT->Clone("SJet80EtapassHTeff");
TH1D *SJet80EtaallMETeff = (TH1D*)SJet80EtaallMET->Clone("SJet80EtaallMETeff");
TH1D *SJet80EtapassMETeff = (TH1D*)SJet80EtapassMET->Clone("SJet80EtapassMETeff");
TEfficiency *SJet80Etaeff = new TEfficiency((*SJet80Etapasseff),(*SJet80Etaalleff));
SJet80Etaeff->SetName("SJet80Etaeff");
TEfficiency *SJet80EtaeffHT = new TEfficiency((*SJet80EtapassHTeff),(*SJet80EtaallHTeff));
SJet80EtaeffHT->SetName("SJet80EtaeffHT");
TEfficiency *SJet80EtaeffMET = new TEfficiency((*SJet80EtapassMETeff),(*SJet80EtaallMETeff));
SJet80EtaeffMET->SetName("SJet80EtaeffMET");
TH1D *SJet100Etaalleff = (TH1D*)SJet100Etaall->Clone("SJet100Etaalleff");
TH1D *SJet100Etapasseff = (TH1D*)SJet100Etapass->Clone("SJet100Etapasseff");
TH1D *SJet100EtaallHTeff = (TH1D*)SJet100EtaallHT->Clone("SJet100EtaallHTeff");
TH1D *SJet100EtapassHTeff = (TH1D*)SJet100EtapassHT->Clone("SJet100EtapassHTeff");
TH1D *SJet100EtaallMETeff = (TH1D*)SJet100EtaallMET->Clone("SJet100EtaallMETeff");
TH1D *SJet100EtapassMETeff = (TH1D*)SJet100EtapassMET->Clone("SJet100EtapassMETeff");
TEfficiency *SJet100Etaeff = new TEfficiency((*SJet100Etapasseff),(*SJet100Etaalleff));
SJet100Etaeff->SetName("SJet100Etaeff");
TEfficiency *SJet100EtaeffHT = new TEfficiency((*SJet100EtapassHTeff),(*SJet100EtaallHTeff));
SJet100EtaeffHT->SetName("SJet100EtaeffHT");
TEfficiency *SJet100EtaeffMET = new TEfficiency((*SJet100EtapassMETeff),(*SJet100EtaallMETeff));
SJet100EtaeffMET->SetName("SJet100EtaeffMET");
TH1D *SJet150Etaalleff = (TH1D*)SJet150Etaall->Clone("SJet150Etaalleff");
TH1D *SJet150Etapasseff = (TH1D*)SJet150Etapass->Clone("SJet150Etapasseff");
TH1D *SJet150EtaallHTeff = (TH1D*)SJet150EtaallHT->Clone("SJet150EtaallHTeff");
TH1D *SJet150EtapassHTeff = (TH1D*)SJet150EtapassHT->Clone("SJet150EtapassHTeff");
TH1D *SJet150EtaallMETeff = (TH1D*)SJet150EtaallMET->Clone("SJet150EtaallMETeff");
TH1D *SJet150EtapassMETeff = (TH1D*)SJet150EtapassMET->Clone("SJet150EtapassMETeff");
TEfficiency *SJet150Etaeff = new TEfficiency((*SJet150Etapasseff),(*SJet150Etaalleff));
SJet150Etaeff->SetName("SJet150Etaeff");
TEfficiency *SJet150EtaeffHT = new TEfficiency((*SJet150EtapassHTeff),(*SJet150EtaallHTeff));
SJet150EtaeffHT->SetName("SJet150EtaeffHT");
TEfficiency *SJet150EtaeffMET = new TEfficiency((*SJet150EtapassMETeff),(*SJet150EtaallMETeff));
SJet150EtaeffMET->SetName("SJet150EtaeffMET");
TH1D *SJet200Etaalleff = (TH1D*)SJet200Etaall->Clone("SJet200Etaalleff");
TH1D *SJet200Etapasseff = (TH1D*)SJet200Etapass->Clone("SJet200Etapasseff");
TH1D *SJet200EtaallHTeff = (TH1D*)SJet200EtaallHT->Clone("SJet200EtaallHTeff");
TH1D *SJet200EtapassHTeff = (TH1D*)SJet200EtapassHT->Clone("SJet200EtapassHTeff");
TH1D *SJet200EtaallMETeff = (TH1D*)SJet200EtaallMET->Clone("SJet200EtaallMETeff");
TH1D *SJet200EtapassMETeff = (TH1D*)SJet200EtapassMET->Clone("SJet200EtapassMETeff");
TEfficiency *SJet200Etaeff = new TEfficiency((*SJet200Etapasseff),(*SJet200Etaalleff));
SJet200Etaeff->SetName("SJet200Etaeff");
TEfficiency *SJet200EtaeffHT = new TEfficiency((*SJet200EtapassHTeff),(*SJet200EtaallHTeff));
SJet200EtaeffHT->SetName("SJet200EtaeffHT");
TEfficiency *SJet200EtaeffMET = new TEfficiency((*SJet200EtapassMETeff),(*SJet200EtaallMETeff));
SJet200EtaeffMET->SetName("SJet200EtaeffMET");
TH1D *SJetPt40alleff = (TH1D*)SJetPt40all->Clone("SJetPt40alleff");
TH1D *SJetPt40passeff = (TH1D*)SJetPt40pass->Clone("SJetPt40passeff");
TH1D *SJetPt40allHTeff = (TH1D*)SJetPt40allHT->Clone("SJetPt40allHTeff");
TH1D *SJetPt40passHTeff = (TH1D*)SJetPt40passHT->Clone("SJetPt40passHTeff");
TH1D *SJetPt40allMETeff = (TH1D*)SJetPt40allMET->Clone("SJetPt40allMETeff");
TH1D *SJetPt40passMETeff = (TH1D*)SJetPt40passMET->Clone("SJetPt40passMETeff");
TEfficiency *SJetPt40eff = new TEfficiency((*SJetPt40passeff),(*SJetPt40alleff));
SJetPt40eff->SetName("SJetPt40eff");
TEfficiency *SJetPt40effHT = new TEfficiency((*SJetPt40passHTeff),(*SJetPt40allHTeff));
SJetPt40effHT->SetName("SJetPt40effHT");
TEfficiency *SJetPt40effMET = new TEfficiency((*SJetPt40passMETeff),(*SJetPt40allMETeff));
SJetPt40effMET->SetName("SJetPt40effMET");
TH1D *SJetPt60alleff = (TH1D*)SJetPt60all->Clone("SJetPt60alleff");
TH1D *SJetPt60passeff = (TH1D*)SJetPt60pass->Clone("SJetPt60passeff");
TH1D *SJetPt60allHTeff = (TH1D*)SJetPt60allHT->Clone("SJetPt60allHTeff");
TH1D *SJetPt60passHTeff = (TH1D*)SJetPt60passHT->Clone("SJetPt60passHTeff");
TH1D *SJetPt60allMETeff = (TH1D*)SJetPt60allMET->Clone("SJetPt60allMETeff");
TH1D *SJetPt60passMETeff = (TH1D*)SJetPt60passMET->Clone("SJetPt60passMETeff");
TEfficiency *SJetPt60eff = new TEfficiency((*SJetPt60passeff),(*SJetPt60alleff));
SJetPt60eff->SetName("SJetPt60eff");
TEfficiency *SJetPt60effHT = new TEfficiency((*SJetPt60passHTeff),(*SJetPt60allHTeff));
SJetPt60effHT->SetName("SJetPt60effHT");
TEfficiency *SJetPt60effMET = new TEfficiency((*SJetPt60passMETeff),(*SJetPt60allMETeff));
SJetPt60effMET->SetName("SJetPt60effMET");
TH1D *SJetPt80alleff = (TH1D*)SJetPt80all->Clone("SJetPt80alleff");
TH1D *SJetPt80passeff = (TH1D*)SJetPt80pass->Clone("SJetPt80passeff");
TH1D *SJetPt80allHTeff = (TH1D*)SJetPt80allHT->Clone("SJetPt80allHTeff");
TH1D *SJetPt80passHTeff = (TH1D*)SJetPt80passHT->Clone("SJetPt80passHTeff");
TH1D *SJetPt80allMETeff = (TH1D*)SJetPt80allMET->Clone("SJetPt80allMETeff");
TH1D *SJetPt80passMETeff = (TH1D*)SJetPt80passMET->Clone("SJetPt80passMETeff");
TEfficiency *SJetPt80eff = new TEfficiency((*SJetPt80passeff),(*SJetPt80alleff));
SJetPt80eff->SetName("SJetPt80eff");
TEfficiency *SJetPt80effHT = new TEfficiency((*SJetPt80passHTeff),(*SJetPt80allHTeff));
SJetPt80effHT->SetName("SJetPt80effHT");
TEfficiency *SJetPt80effMET = new TEfficiency((*SJetPt80passMETeff),(*SJetPt80allMETeff));
SJetPt80effMET->SetName("SJetPt80effMET");
TH1D *SJetPt100alleff = (TH1D*)SJetPt100all->Clone("SJetPt100alleff");
TH1D *SJetPt100passeff = (TH1D*)SJetPt100pass->Clone("SJetPt100passeff");
TH1D *SJetPt100allHTeff = (TH1D*)SJetPt100allHT->Clone("SJetPt100allHTeff");
TH1D *SJetPt100passHTeff = (TH1D*)SJetPt100passHT->Clone("SJetPt100passHTeff");
TH1D *SJetPt100allMETeff = (TH1D*)SJetPt100allMET->Clone("SJetPt100allMETeff");
TH1D *SJetPt100passMETeff = (TH1D*)SJetPt100passMET->Clone("SJetPt100passMETeff");
TEfficiency *SJetPt100eff = new TEfficiency((*SJetPt100passeff),(*SJetPt100alleff));
SJetPt100eff->SetName("SJetPt100eff");
TEfficiency *SJetPt100effHT = new TEfficiency((*SJetPt100passHTeff),(*SJetPt100allHTeff));
SJetPt100effHT->SetName("SJetPt100effHT");
TEfficiency *SJetPt100effMET = new TEfficiency((*SJetPt100passMETeff),(*SJetPt100allMETeff));
SJetPt100effMET->SetName("SJetPt100effMET");
TH1D *SJetPt150alleff = (TH1D*)SJetPt150all->Clone("SJetPt150alleff");
TH1D *SJetPt150passeff = (TH1D*)SJetPt150pass->Clone("SJetPt150passeff");
TH1D *SJetPt150allHTeff = (TH1D*)SJetPt150allHT->Clone("SJetPt150allHTeff");
TH1D *SJetPt150passHTeff = (TH1D*)SJetPt150passHT->Clone("SJetPt150passHTeff");
TH1D *SJetPt150allMETeff = (TH1D*)SJetPt150allMET->Clone("SJetPt150allMETeff");
TH1D *SJetPt150passMETeff = (TH1D*)SJetPt150passMET->Clone("SJetPt150passMETeff");
TEfficiency *SJetPt150eff = new TEfficiency((*SJetPt150passeff),(*SJetPt150alleff));
SJetPt150eff->SetName("SJetPt150eff");
TEfficiency *SJetPt150effHT = new TEfficiency((*SJetPt150passHTeff),(*SJetPt150allHTeff));
SJetPt150effHT->SetName("SJetPt150effHT");
TEfficiency *SJetPt150effMET = new TEfficiency((*SJetPt150passMETeff),(*SJetPt150allMETeff));
SJetPt150effMET->SetName("SJetPt150effMET");
TH1D *SJetPt200alleff = (TH1D*)SJetPt200all->Clone("SJetPt200alleff");
TH1D *SJetPt200passeff = (TH1D*)SJetPt200pass->Clone("SJetPt200passeff");
TH1D *SJetPt200allHTeff = (TH1D*)SJetPt200allHT->Clone("SJetPt200allHTeff");
TH1D *SJetPt200passHTeff = (TH1D*)SJetPt200passHT->Clone("SJetPt200passHTeff");
TH1D *SJetPt200allMETeff = (TH1D*)SJetPt200allMET->Clone("SJetPt200allMETeff");
TH1D *SJetPt200passMETeff = (TH1D*)SJetPt200passMET->Clone("SJetPt200passMETeff");
TEfficiency *SJetPt200eff = new TEfficiency((*SJetPt200passeff),(*SJetPt200alleff));
SJetPt200eff->SetName("SJetPt200eff");
TEfficiency *SJetPt200effHT = new TEfficiency((*SJetPt200passHTeff),(*SJetPt200allHTeff));
SJetPt200effHT->SetName("SJetPt200effHT");
TEfficiency *SJetPt200effMET = new TEfficiency((*SJetPt200passMETeff),(*SJetPt200allMETeff));
SJetPt200effMET->SetName("SJetPt200effMET");


//save everything
TFile *newfile;
     if(scale>1.01)      newfile = new TFile("~/TOBTECfiles/PassingEventsStudy2_NoMinDPhiNoLeptonVetoIncludingPhotons_JetScaledUp.root","RECREATE");
else if(scale<0.99)      newfile = new TFile("~/TOBTECfiles/PassingEventsStudy2_NoMinDPhiNoLeptonVetoIncludingPhotons_JetScaledDown.root","RECREATE");
else                     newfile = new TFile("~/TOBTECfiles/PassingEventsStudy2_NoMinDPhiNoLeptonVetoIncludingPhotons.root","RECREATE");
newfile ->cd();
MT2all     ->Write();
MT2allHT   ->Write();
MT2allMET  ->Write();
MT2pass    ->Write();
MT2passHT  ->Write();
MT2passMET ->Write();
MT2eff     ->Write();
MT2effHT   ->Write();
MT2effMET  ->Write();
MT2fit     ->Write();
MT2fitHT   ->Write();
MT2fitMET  ->Write();
HTall     ->Write();
HTallHT   ->Write();
HTallMET  ->Write();
HTpass    ->Write();
HTpassHT  ->Write();
HTpassMET ->Write();
HTeff     ->Write();
HTeffHT   ->Write();
HTeffMET  ->Write();
HTfit1    ->Write();
HTfit2    ->Write();
HTfitHT   ->Write();
HTfitMET  ->Write();
METall     ->Write();
METallHT   ->Write();
METallMET  ->Write();
METpass    ->Write();
METpassHT  ->Write();
METpassMET ->Write();
METeff     ->Write();
METeffHT   ->Write();
METeffMET  ->Write();
METfit     ->Write();
METfitHT   ->Write();
METfitMET  ->Write();
NJetsall     ->Write();
NJetsallHT   ->Write();
NJetsallMET  ->Write();
NJetspass    ->Write();
NJetspassHT  ->Write();
NJetspassMET ->Write();
NJetseff     ->Write();
NJetseffHT   ->Write();
NJetseffMET  ->Write();
NJetsfit     ->Write();
NJetsfitHT   ->Write();
NJetsfitMET  ->Write();
NBJetsall     ->Write();
NBJetsallHT   ->Write();
NBJetsallMET  ->Write();
NBJetspass    ->Write();
NBJetspassHT  ->Write();
NBJetspassMET ->Write();
NBJetseff     ->Write();
NBJetseffHT   ->Write();
NBJetseffMET  ->Write();
NBJetsfit     ->Write();
NBJetsfitHT   ->Write();
NBJetsfitMET  ->Write();
Jet1Etaall     ->Write();
Jet1EtaallHT   ->Write();
Jet1EtaallMET  ->Write();
Jet1Etapass    ->Write();
Jet1EtapassHT  ->Write();
Jet1EtapassMET ->Write();
Jet1Etaeff     ->Write();
Jet1EtaeffHT   ->Write();
Jet1EtaeffMET  ->Write();
Jet2Etaall     ->Write();
Jet2EtaallHT   ->Write();
Jet2EtaallMET  ->Write();
Jet2Etapass    ->Write();
Jet2EtapassHT  ->Write();
Jet2EtapassMET ->Write();
Jet2Etaeff     ->Write();
Jet2EtaeffHT   ->Write();
Jet2EtaeffMET  ->Write();
Jet3Etaall     ->Write();
Jet3EtaallHT   ->Write();
Jet3EtaallMET  ->Write();
Jet3Etapass    ->Write();
Jet3EtapassHT  ->Write();
Jet3EtapassMET ->Write();
Jet3Etaeff     ->Write();
Jet3EtaeffHT   ->Write();
Jet3EtaeffMET  ->Write();
Jet4Etaall     ->Write();
Jet4EtaallHT   ->Write();
Jet4EtaallMET  ->Write();
Jet4Etapass    ->Write();
Jet4EtapassHT  ->Write();
Jet4EtapassMET ->Write();
Jet4Etaeff     ->Write();
Jet4EtaeffHT   ->Write();
Jet4EtaeffMET  ->Write();
JetPtEtaall     ->Write();
JetPtEtaallHT   ->Write();
JetPtEtaallMET  ->Write();
JetPtEtapass    ->Write();
JetPtEtapassHT  ->Write();
JetPtEtapassMET ->Write();
JetPtEtaeff     ->Write();
JetPtEtaeffHT   ->Write();
JetPtEtaeffMET  ->Write();
JetEtaall     ->Write();
JetEtaallHT   ->Write();
JetEtaallMET  ->Write();
JetEtapass    ->Write();
JetEtapassHT  ->Write();
JetEtapassMET ->Write();
JetEtaeff     ->Write();
JetEtaeffHT   ->Write();
JetEtaeffMET  ->Write();
JetPtall     ->Write();
JetPtallHT   ->Write();
JetPtallMET  ->Write();
JetPtpass    ->Write();
JetPtpassHT  ->Write();
JetPtpassMET ->Write();
JetPteff     ->Write();
JetPteffHT   ->Write();
JetPteffMET  ->Write();
SJetPt40all     ->Write();
SJetPt40allHT   ->Write();
SJetPt40allMET  ->Write();
SJetPt40pass    ->Write();
SJetPt40passHT  ->Write();
SJetPt40passMET ->Write();
SJetPt40eff     ->Write();
SJetPt40effHT   ->Write();
SJetPt40effMET  ->Write();
SJetPt60all     ->Write();
SJetPt60allHT   ->Write();
SJetPt60allMET  ->Write();
SJetPt60pass    ->Write();
SJetPt60passHT  ->Write();
SJetPt60passMET ->Write();
SJetPt60eff     ->Write();
SJetPt60effHT   ->Write();
SJetPt60effMET  ->Write();
SJetPt80all     ->Write();
SJetPt80allHT   ->Write();
SJetPt80allMET  ->Write();
SJetPt80pass    ->Write();
SJetPt80passHT  ->Write();
SJetPt80passMET ->Write();
SJetPt80eff     ->Write();
SJetPt80effHT   ->Write();
SJetPt80effMET  ->Write();
SJetPt100all    ->Write();
SJetPt100allHT  ->Write();
SJetPt100allMET ->Write();
SJetPt100pass   ->Write();
SJetPt100passHT ->Write();
SJetPt100passMET->Write();
SJetPt100eff    ->Write();
SJetPt100effHT  ->Write();
SJetPt100effMET ->Write();
SJetPt150all    ->Write();
SJetPt150allHT  ->Write();
SJetPt150allMET ->Write();
SJetPt150pass   ->Write();
SJetPt150passHT ->Write();
SJetPt150passMET->Write();
SJetPt150eff    ->Write();
SJetPt150effHT  ->Write();
SJetPt150effMET ->Write();
SJetPt200all    ->Write();
SJetPt200allHT  ->Write();
SJetPt200allMET ->Write();
SJetPt200pass   ->Write();
SJetPt200passHT ->Write();
SJetPt200passMET->Write();
SJetPt200eff    ->Write();
SJetPt200effHT  ->Write();
SJetPt200effMET ->Write();
SJet40Etaall     ->Write();
SJet40EtaallHT   ->Write();
SJet40EtaallMET  ->Write();
SJet40Etapass    ->Write();
SJet40EtapassHT  ->Write();
SJet40EtapassMET ->Write();
SJet40Etaeff     ->Write();
SJet40EtaeffHT   ->Write();
SJet40EtaeffMET  ->Write();
SJet60Etaall     ->Write();
SJet60EtaallHT   ->Write();
SJet60EtaallMET  ->Write();
SJet60Etapass    ->Write();
SJet60EtapassHT  ->Write();
SJet60EtapassMET ->Write();
SJet60Etaeff     ->Write();
SJet60EtaeffHT   ->Write();
SJet60EtaeffMET  ->Write();
SJet80Etaall     ->Write();
SJet80EtaallHT   ->Write();
SJet80EtaallMET  ->Write();
SJet80Etapass    ->Write();
SJet80EtapassHT  ->Write();
SJet80EtapassMET ->Write();
SJet80Etaeff     ->Write();
SJet80EtaeffHT   ->Write();
SJet80EtaeffMET  ->Write();
SJet100Etaall    ->Write();
SJet100EtaallHT  ->Write();
SJet100EtaallMET ->Write();
SJet100Etapass   ->Write();
SJet100EtapassHT ->Write();
SJet100EtapassMET->Write();
SJet100Etaeff    ->Write();
SJet100EtaeffHT  ->Write();
SJet100EtaeffMET ->Write();
SJet150Etaall    ->Write();
SJet150EtaallHT  ->Write();
SJet150EtaallMET ->Write();
SJet150Etapass   ->Write();
SJet150EtapassHT ->Write();
SJet150EtapassMET->Write();
SJet150Etaeff    ->Write();
SJet150EtaeffHT  ->Write();
SJet150EtaeffMET ->Write();
SJet200Etaall    ->Write();
SJet200EtaallHT  ->Write();
SJet200EtaallMET ->Write();
SJet200Etapass   ->Write();
SJet200EtapassHT ->Write();
SJet200EtapassMET->Write();
SJet200Etaeff    ->Write();
SJet200EtaeffHT  ->Write();
SJet200EtaeffMET ->Write();
SJetPt40Etaall     ->Write();
SJetPt40EtaallHT   ->Write();
SJetPt40EtaallMET  ->Write();
SJetPt40Etapass    ->Write();
SJetPt40EtapassHT  ->Write();
SJetPt40EtapassMET ->Write();
SJetPt40Etaeff     ->Write();
SJetPt40EtaeffHT   ->Write();
SJetPt40EtaeffMET  ->Write();
SJetPt60Etaall     ->Write();
SJetPt60EtaallHT   ->Write();
SJetPt60EtaallMET  ->Write();
SJetPt60Etapass    ->Write();
SJetPt60EtapassHT  ->Write();
SJetPt60EtapassMET ->Write();
SJetPt60Etaeff     ->Write();
SJetPt60EtaeffHT   ->Write();
SJetPt60EtaeffMET  ->Write();
SJetPt80Etaall     ->Write();
SJetPt80EtaallHT   ->Write();
SJetPt80EtaallMET  ->Write();
SJetPt80Etapass    ->Write();
SJetPt80EtapassHT  ->Write();
SJetPt80EtapassMET ->Write();
SJetPt80Etaeff     ->Write();
SJetPt80EtaeffHT   ->Write();
SJetPt80EtaeffMET  ->Write();
SJetPt100Etaall    ->Write();
SJetPt100EtaallHT  ->Write();
SJetPt100EtaallMET ->Write();
SJetPt100Etapass   ->Write();
SJetPt100EtapassHT ->Write();
SJetPt100EtapassMET->Write();
SJetPt100Etaeff    ->Write();
SJetPt100EtaeffHT  ->Write();
SJetPt100EtaeffMET ->Write();
SJetPt150Etaall    ->Write();
SJetPt150EtaallHT  ->Write();
SJetPt150EtaallMET ->Write();
SJetPt150Etapass   ->Write();
SJetPt150EtapassHT ->Write();
SJetPt150EtapassMET->Write();
SJetPt150Etaeff    ->Write();
SJetPt150EtaeffHT  ->Write();
SJetPt150EtaeffMET ->Write();
SJetPt200Etaall    ->Write();
SJetPt200EtaallHT  ->Write();
SJetPt200EtaallMET ->Write();
SJetPt200Etapass   ->Write();
SJetPt200EtapassHT ->Write();
SJetPt200EtapassMET->Write();
SJetPt200Etaeff    ->Write();
SJetPt200EtaeffHT  ->Write();
SJetPt200EtaeffMET ->Write();

newfile ->Close();

cout << "histograms stored in " << newfile->GetName() << endl;
//gStyle->SetPaintTextFormat("2.2f");

}