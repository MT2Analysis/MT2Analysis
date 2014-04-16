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
#include "TChain.h"
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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//user your own path of MT2tree.hh
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"

using namespace std;

//run via root -l -b -q MT2vsMET.C++
//makes 2d correlation plots of MT2 vs MET.
void MT2vsMET(){

//TString samples = "samples/samples_HTandMET_filter.dat";

gROOT->cd();
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

TChain *cLM6 = new TChain("MassTree");
TChain *cWJets = new TChain("MassTree");
TChain *cZJets = new TChain("MassTree");
TChain *cTTbar = new TChain("MassTree");
TChain *QCD300 = new TChain("MassTree");
TChain *QCD470 = new TChain("MassTree");
TChain *QCD600 = new TChain("MassTree");
TChain *QCD800 = new TChain("MassTree");
TChain *QCD1000 = new TChain("MassTree");
TChain *QCD1400 = new TChain("MassTree");
TChain *QCD1800 = new TChain("MassTree");


cLM6->Add("~mmasciov/MT2Analysis/MT2trees/MT2_V02-03-02/20130430_8TeV/highHT/SUSY-LM6-sftsht-8TeV-pythia6-Summer12-DR53X-PU-S10-START53-V7A-v1.root");
cLM6->Add("~mmasciov/MT2Analysis/MT2trees/MT2_V02-03-02/20130430_8TeV/lowHT/SUSY-LM6-sftsht-8TeV-pythia6-Summer12-DR53X-PU-S10-START53-V7A-v1.root");
cWJets->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/highHT/WJetsToLNu-HT-400ToInf*.root");
cWJets->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/WJetsToLNu-HT-400ToInf*.root");
cZJets->Add("~haweber/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/highHT/ZJetsToNuNu-400-HT-inf*.root");
cZJets->Add("~haweber/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/ZJetsToNuNu-400-HT-inf*.root");
cTTbar->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130503_8TeV/highHT/TTJets-SemiLeptMGDecays-8TeV-madgraph-tauola-Summer12-DR53X-PU-S10-START53-V7C-v1.root");
cTTbar->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130503_8TeV/lowHT/TTJets-SemiLeptMGDecays-8TeV-madgraph-tauola-Summer12-DR53X-PU-S10-START53-V7C-v1.root");
QCD300->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/QCD-Pt-300to470*.root");
QCD300->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/highHT/QCD-Pt-300to470*.root");
QCD470->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/QCD-Pt-470to600*.root");
QCD470->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/highHT/QCD-Pt-470to600*.root");
QCD600->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/QCD-Pt-600to800*.root");
QCD600->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/highHT/QCD-Pt-600to800*.root");
QCD800->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/QCD-Pt-800to1000*.root");
QCD800->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/highHT/QCD-Pt-800to1000*.root");
QCD1000->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/QCD-Pt-1000to1400*.root");
QCD1000->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/highHT/QCD-Pt-1000to1400*.root");
QCD1400->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/QCD-Pt-1400to1800*.root");
QCD1400->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/highHT/QCD-Pt-1400to1800*.root");
QCD1800->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/QCD-Pt-1800*.root");
QCD1800->Add("~casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/highHT/QCD-Pt-1800*.root");

TH2D *LM6      = new TH2D("LM6"     ,"", 60, 0, 1500, 60, 0, 1500); LM6     ->Sumw2();
TH2D *LM6lHT   = new TH2D("LM6lHT"  ,"", 60, 0, 1500, 60, 0, 1500); LM6lHT  ->Sumw2();
TH2D *LM6mHT   = new TH2D("LM6mHT"  ,"", 60, 0, 1500, 60, 0, 1500); LM6mHT  ->Sumw2();
TH2D *LM6hHT   = new TH2D("LM6hHT"  ,"", 60, 0, 1500, 60, 0, 1500); LM6hHT  ->Sumw2();

TH2D *QCD      = new TH2D("QCD"     ,"", 28, 30, 450, 30, 0, 450); QCD     ->Sumw2();
TH2D *QCDlHT   = new TH2D("QCDlHT"  ,"", 28, 30, 450, 30, 0, 450); QCDlHT  ->Sumw2();
TH2D *QCDmHT   = new TH2D("QCDmHT"  ,"", 28, 30, 450, 30, 0, 450); QCDmHT  ->Sumw2();
TH2D *QCDhHT   = new TH2D("QCDhHT"  ,"", 28, 30, 450, 30, 0, 450); QCDhHT  ->Sumw2();

TH2D *WJets      = new TH2D("WJets"     ,"", 39, 30, 1005, 40, 0, 1000); WJets       ->Sumw2();
TH2D *WJetslHT   = new TH2D("WJetslHT"  ,"", 39, 30, 1005, 40, 0, 1000); WJetslHT    ->Sumw2();
TH2D *WJetsmHT   = new TH2D("WJetsmHT"  ,"", 39, 30, 1005, 40, 0, 1000); WJetsmHT    ->Sumw2();
TH2D *WJetshHT   = new TH2D("WJetshHT"  ,"", 39, 30, 1005, 40, 0, 1000); WJetshHT    ->Sumw2();

TH2D *ZJets      = new TH2D("ZJets"     ,"", 59, 30, 1505, 60, 0, 1500); ZJets     ->Sumw2();
TH2D *ZJetslHT   = new TH2D("ZJetslHT"  ,"", 59, 30, 1505, 60, 0, 1500); ZJetslHT  ->Sumw2();
TH2D *ZJetsmHT   = new TH2D("ZJetsmHT"  ,"", 59, 30, 1505, 60, 0, 1500); ZJetsmHT  ->Sumw2();
TH2D *ZJetshHT   = new TH2D("ZJetshHT"  ,"", 59, 30, 1505, 60, 0, 1500); ZJetshHT  ->Sumw2();

TH2D *TTbar      = new TH2D("TTbar"     ,"", 39, 30, 1005, 40, 0, 1000); TTbar     ->Sumw2();
TH2D *TTbarlHT   = new TH2D("TTbarlHT"  ,"", 39, 30, 1005, 40, 0, 1000); TTbarlHT  ->Sumw2();
TH2D *TTbarmHT   = new TH2D("TTbarmHT"  ,"", 39, 30, 1005, 40, 0, 1000); TTbarmHT  ->Sumw2();
TH2D *TTbarhHT   = new TH2D("TTbarhHT"  ,"", 39, 30, 1005, 40, 0, 1000); TTbarhHT  ->Sumw2();

gROOT->cd();
cout << "LM6" << endl;
cLM6->Draw("misc.MT2:misc.MET>>LM6",      "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30","goff");
cout << "LM6lHT" << endl;
cLM6->Draw("misc.MT2:misc.MET>>LM6lHT",   "NJetsIDLoose40>=2&&misc.HT>450&&misc.HT<750&&misc.MET>30","goff");
cout << "LM6mHT" << endl;
cLM6->Draw("misc.MT2:misc.MET>>LM6mHT",   "NJetsIDLoose40>=2&&misc.HT>750&&misc.HT<1200&&misc.MET>30","goff");
cout << "LM6hHT" << endl;
cLM6->Draw("misc.MT2:misc.MET>>LM6hHT",   "NJetsIDLoose40>=2&&misc.HT>1200&&misc.MET>30","goff");

cout << "WJets" << endl;
cWJets->Draw("misc.MT2:misc.MET>>WJets",      "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30","goff");
cout << "WJetslHT" << endl;
cWJets->Draw("misc.MT2:misc.MET>>WJetslHT",   "NJetsIDLoose40>=2&&misc.HT>450&&misc.HT<750&&misc.MET>30","goff");
cout << "WJetsmHT" << endl;
cWJets->Draw("misc.MT2:misc.MET>>WJetsmHT",   "NJetsIDLoose40>=2&&misc.HT>750&&misc.HT<1200&&misc.MET>30","goff");
cout << "WJetshHT" << endl;
cWJets->Draw("misc.MT2:misc.MET>>WJetshHT",   "NJetsIDLoose40>=2&&misc.HT>1200&&misc.MET>30","goff");

cout << "ZJets" << endl;
cZJets->Draw("misc.MT2:misc.MET>>ZJets",      "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30","goff");
cout << "ZJetslHT" << endl;
cZJets->Draw("misc.MT2:misc.MET>>ZJetslHT",   "NJetsIDLoose40>=2&&misc.HT>450&&misc.HT<750&&misc.MET>30","goff");
cout << "ZJetsmHT" << endl;
cZJets->Draw("misc.MT2:misc.MET>>ZJetsmHT",   "NJetsIDLoose40>=2&&misc.HT>750&&misc.HT<1200&&misc.MET>30","goff");
cout << "ZJetshHT" << endl;
cZJets->Draw("misc.MT2:misc.MET>>ZJetshHT",   "NJetsIDLoose40>=2&&misc.HT>1200&&misc.MET>30","goff");

cout << "TTbar" << endl;
cTTbar->Draw("misc.MT2:misc.MET>>TTbar",      "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30","goff");
cout << "TTbarlHT" << endl;
cTTbar->Draw("misc.MT2:misc.MET>>TTbarlHT",   "NJetsIDLoose40>=2&&misc.HT>450&&misc.HT<750&&misc.MET>30","goff");
cout << "TTbarmHT" << endl;
cTTbar->Draw("misc.MT2:misc.MET>>TTbarmHT",   "NJetsIDLoose40>=2&&misc.HT>750&&misc.HT<1200&&misc.MET>30","goff");
cout << "TTbarhHT" << endl;
cTTbar->Draw("misc.MT2:misc.MET>>TTbarhHT",   "NJetsIDLoose40>=2&&misc.HT>1200&&misc.MET>30","goff");

cout << "QCD300" << endl;
QCD300->Draw("misc.MT2:misc.MET>>+QCD",      "(1.172) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30)","goff");
cout << "QCDlHT" << endl;
QCD300->Draw("misc.MT2:misc.MET>>+QCDlHT",   "(1.172) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.HT<750&&misc.MET>30)","goff");
cout << "QCDmHT" << endl;
QCD300->Draw("misc.MT2:misc.MET>>+QCDmHT",   "(1.172) * (NJetsIDLoose40>=2&&misc.HT>750&&misc.HT<1200&&misc.MET>30)","goff");
cout << "QCDhHT" << endl;
QCD300->Draw("misc.MT2:misc.MET>>+QCDhHT",   "(1.172) * (NJetsIDLoose40>=2&&misc.HT>1200&&misc.MET>30)","goff");
cout << "QCD470" << endl;
QCD470->Draw("misc.MT2:misc.MET>>+QCD",      "(0.5568) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30)","goff");
cout << "QCDlHT" << endl;
QCD470->Draw("misc.MT2:misc.MET>>+QCDlHT",   "(0.5568) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.HT<750&&misc.MET>30)","goff");
cout << "QCDmHT" << endl;
QCD470->Draw("misc.MT2:misc.MET>>+QCDmHT",   "(0.5568) * (NJetsIDLoose40>=2&&misc.HT>750&&misc.HT<1200&&misc.MET>30)","goff");
cout << "QCDhHT" << endl;
QCD470->Draw("misc.MT2:misc.MET>>+QCDhHT",   "(0.5568) * (NJetsIDLoose40>=2&&misc.HT>1200&&misc.MET>30)","goff");
cout << "QCD600" << endl;
QCD600->Draw("misc.MT2:misc.MET>>+QCD",      "(0.1319) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30)","goff");
cout << "QCDlHT" << endl;
QCD600->Draw("misc.MT2:misc.MET>>+QCDlHT",   "(0.1319) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.HT<750&&misc.MET>30)","goff");
cout << "QCDmHT" << endl;
QCD600->Draw("misc.MT2:misc.MET>>+QCDmHT",   "(0.1319) * (NJetsIDLoose40>=2&&misc.HT>750&&misc.HT<1200&&misc.MET>30)","goff");
cout << "QCDhHT" << endl;
QCD600->Draw("misc.MT2:misc.MET>>+QCDhHT",   "(0.1319) * (NJetsIDLoose40>=2&&misc.HT>1200&&misc.MET>30)","goff");
cout << "QCD800" << endl;
QCD800->Draw("misc.MT2:misc.MET>>+QCD",      "(0.01733) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30)","goff");
cout << "QCDlHT" << endl;
QCD800->Draw("misc.MT2:misc.MET>>+QCDlHT",   "(0.01733) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.HT<750&&misc.MET>30)","goff");
cout << "QCDmHT" << endl;
QCD800->Draw("misc.MT2:misc.MET>>+QCDmHT",   "(0.01733) * (NJetsIDLoose40>=2&&misc.HT>750&&misc.HT<1200&&misc.MET>30)","goff");
cout << "QCDhHT" << endl;
QCD800->Draw("misc.MT2:misc.MET>>+QCDhHT",   "(0.01733) * (NJetsIDLoose40>=2&&misc.HT>1200&&misc.MET>30)","goff");
cout << "QCD1000" << endl;
QCD1000->Draw("misc.MT2:misc.MET>>+QCD",      "(0.007326) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30)","goff");
cout << "QCDlHT" << endl;
QCD1000->Draw("misc.MT2:misc.MET>>+QCDlHT",   "(0.007326) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.HT<750&&misc.MET>30)","goff");
cout << "QCDmHT" << endl;
QCD1000->Draw("misc.MT2:misc.MET>>+QCDmHT",   "(0.007326) * (NJetsIDLoose40>=2&&misc.HT>750&&misc.HT<1200&&misc.MET>30)","goff");
cout << "QCDhHT" << endl;
QCD1000->Draw("misc.MT2:misc.MET>>+QCDhHT",   "(0.007326) * (NJetsIDLoose40>=2&&misc.HT>1200&&misc.MET>30)","goff");
cout << "QCD1400" << endl;
QCD1400->Draw("misc.MT2:misc.MET>>+QCD",      "(0.0003271) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30)","goff");
cout << "QCDlHT" << endl;
QCD1400->Draw("misc.MT2:misc.MET>>+QCDlHT",   "(0.0003271) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.HT<750&&misc.MET>30)","goff");
cout << "QCDmHT" << endl;
QCD1400->Draw("misc.MT2:misc.MET>>+QCDmHT",   "(0.0003271) * (NJetsIDLoose40>=2&&misc.HT>750&&misc.HT<1200&&misc.MET>30)","goff");
cout << "QCDhHT" << endl;
QCD1400->Draw("misc.MT2:misc.MET>>+QCDhHT",   "(0.0003271) * (NJetsIDLoose40>=2&&misc.HT>1200&&misc.MET>30)","goff");
cout << "QCD1800" << endl;
QCD1800->Draw("misc.MT2:misc.MET>>+QCD",      "(0.00003648) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30)","goff");
cout << "QCDlHT" << endl;
QCD1800->Draw("misc.MT2:misc.MET>>+QCDlHT",   "(0.00003648) * (NJetsIDLoose40>=2&&misc.HT>450&&misc.HT<750&&misc.MET>30)","goff");
cout << "QCDmHT" << endl;
QCD1800->Draw("misc.MT2:misc.MET>>+QCDmHT",   "(0.00003648) * (NJetsIDLoose40>=2&&misc.HT>750&&misc.HT<1200&&misc.MET>30)","goff");
cout << "QCDhHT" << endl;
QCD1800->Draw("misc.MT2:misc.MET>>+QCDhHT",   "(0.00003648) * (NJetsIDLoose40>=2&&misc.HT>1200&&misc.MET>30)","goff");


TFile *file = new TFile("MT2vsMETfile.root","RECREATE");
file->cd();
LM6   ->Write();
LM6lHT->Write();
LM6mHT->Write();
LM6hHT->Write();
QCD   ->Write();
QCDlHT->Write();
QCDmHT->Write();
QCDhHT->Write();
ZJets   ->Write();
ZJetslHT->Write();
ZJetsmHT->Write();
ZJetshHT->Write();
WJets   ->Write();
WJetslHT->Write();
WJetsmHT->Write();
WJetshHT->Write();
TTbar   ->Write();
TTbarlHT->Write();
TTbarmHT->Write();
TTbarhHT->Write();


file->Close();

cout << "File saved: " << file->GetName() << endl;

}