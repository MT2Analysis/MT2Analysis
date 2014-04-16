#include "TEfficiency.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TLegend.h"
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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path of MT2tree.hh
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"

//run via root -l -b -q UTMplots.C++

using namespace std;

//this makes shape plots of MT2 distributions for several VSPT bins for various SM types (W,Z,ttbar,susy,data)
void UTMplots(){

gROOT->cd();
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

TChain *LM6 = new TChain("MassTree");
TChain *WJets = new TChain("MassTree");
TChain *ZJets = new TChain("MassTree");
TChain *TTbar = new TChain("MassTree");
TChain *Data = new TChain("MassTree");

//read in trees
LM6->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/SUSY/MassTrees/MT2_V02-03-02/20130430_8TeV/SUSY-LM6-sftsht-8TeV-pythia6-Summer12-DR53X-PU-S10-START53-V7A-v1/*.root");
WJets->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/WJetsToLNu-HT-400ToInf-8TeV-madgraph-Summer12-DR53X-PU-S10-START53-V7A-v1/*.root");
WJets->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/WJetsToLNu-HT-400ToInf-8TeV-madgraph-v2-Summer12-DR53X-PU-S10-START53-V7A-v1/*.root");
ZJets->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/ZJetsToNuNu-400-HT-inf-TuneZ2Star-8TeV-madgraph-Summer12-DR53X-PU-S10-START53-V7A-v1/*.root");
ZJets->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/ZJetsToNuNu-400-HT-inf-TuneZ2Star-8TeV-madgraph-ext-Summer12-DR53X-PU-S10-START53-V7A-v1/*.root");
TTbar->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/SUSY/MassTrees/MT2_V02-03-02/20130503_8TeV/TTJets-SemiLeptMGDecays-8TeV-madgraph-tauola-Summer12-DR53X-PU-S10-START53-V7C-v1/*.root");
int nd = Data->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/HT-Run2012A-13Jul2012-v1-2/*.root");
nd += Data->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/HT-Run2012A-recover-06Aug2012-v1-2/*.root");
nd += Data->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/JetHT-Run2012B-13Jul2012-v1-2/*.root");
nd += Data->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/JetHT-Run2012C-24Aug2012-v2-2/*.root");
nd += Data->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/JetHT-Run2012C-EcalRecover-11Dec2012-v1-2/*.root");
nd += Data->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/JetHT-Run2012C-PromptReco-v2-3/*.root");
nd += Data->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/JetHT-Run2012D-PromptReco-v1-3/*.root");
cout << nd << endl;

//define histograms
TH1D *LM6VSPT0   = new TH1D("LM6VSPT0"  ,"", 75, 0, 750); LM6VSPT0  ->Sumw2(); LM6VSPT0  ->SetLineColor(kBlack);  LM6VSPT0  ->SetLineWidth(2);
TH1D *LM6VSPT20  = new TH1D("LM6VSPT20" ,"", 75, 0, 750); LM6VSPT20 ->Sumw2(); LM6VSPT20 ->SetLineColor(kRed);    LM6VSPT20 ->SetLineWidth(2);
TH1D *LM6VSPT30  = new TH1D("LM6VSPT30" ,"", 75, 0, 750); LM6VSPT30 ->Sumw2(); LM6VSPT30 ->SetLineColor(kGreen);  LM6VSPT30 ->SetLineWidth(2);
TH1D *LM6VSPT40  = new TH1D("LM6VSPT40" ,"", 75, 0, 750); LM6VSPT40 ->Sumw2(); LM6VSPT40 ->SetLineColor(kBlue);   LM6VSPT40 ->SetLineWidth(2);
TH1D *LM6VSPT50  = new TH1D("LM6VSPT50" ,"", 75, 0, 750); LM6VSPT50 ->Sumw2(); LM6VSPT50 ->SetLineColor(kYellow); LM6VSPT50 ->SetLineWidth(2);
TH1D *LM6VSPT70  = new TH1D("LM6VSPT70" ,"", 75, 0, 750); LM6VSPT70 ->Sumw2(); LM6VSPT70 ->SetLineColor(kMagenta);LM6VSPT70 ->SetLineWidth(2);
TH1D *LM6VSPT100 = new TH1D("LM6VSPT100","", 75, 0, 750); LM6VSPT100->Sumw2(); LM6VSPT100->SetLineColor(kCyan);   LM6VSPT100->SetLineWidth(2);
TH1D *WJetsVSPT0   = new TH1D("WJetsVSPT0"  ,"", 75, 0, 750); WJetsVSPT0  ->Sumw2(); WJetsVSPT0  ->SetLineColor(kBlack);   WJetsVSPT0  ->SetLineWidth(2);
TH1D *WJetsVSPT20  = new TH1D("WJetsVSPT20" ,"", 75, 0, 750); WJetsVSPT20 ->Sumw2(); WJetsVSPT20 ->SetLineColor(kRed);     WJetsVSPT20 ->SetLineWidth(2);
TH1D *WJetsVSPT30  = new TH1D("WJetsVSPT30" ,"", 75, 0, 750); WJetsVSPT30 ->Sumw2(); WJetsVSPT30 ->SetLineColor(kGreen);   WJetsVSPT30 ->SetLineWidth(2);
TH1D *WJetsVSPT40  = new TH1D("WJetsVSPT40" ,"", 75, 0, 750); WJetsVSPT40 ->Sumw2(); WJetsVSPT40 ->SetLineColor(kBlue);    WJetsVSPT40 ->SetLineWidth(2);
TH1D *WJetsVSPT50  = new TH1D("WJetsVSPT50" ,"", 75, 0, 750); WJetsVSPT50 ->Sumw2(); WJetsVSPT50 ->SetLineColor(kYellow);  WJetsVSPT50 ->SetLineWidth(2);
TH1D *WJetsVSPT70  = new TH1D("WJetsVSPT70" ,"", 75, 0, 750); WJetsVSPT70 ->Sumw2(); WJetsVSPT70 ->SetLineColor(kMagenta); WJetsVSPT70 ->SetLineWidth(2);
TH1D *WJetsVSPT100 = new TH1D("WJetsVSPT100","", 75, 0, 750); WJetsVSPT100->Sumw2(); WJetsVSPT100->SetLineColor(kCyan);    WJetsVSPT100->SetLineWidth(2);
TH1D *ZJetsVSPT0   = new TH1D("ZJetsVSPT0"  ,"", 75, 0, 750); ZJetsVSPT0  ->Sumw2(); ZJetsVSPT0  ->SetLineColor(kBlack);   ZJetsVSPT0  ->SetLineWidth(2);
TH1D *ZJetsVSPT20  = new TH1D("ZJetsVSPT20" ,"", 75, 0, 750); ZJetsVSPT20 ->Sumw2(); ZJetsVSPT20 ->SetLineColor(kRed);     ZJetsVSPT20 ->SetLineWidth(2);
TH1D *ZJetsVSPT30  = new TH1D("ZJetsVSPT30" ,"", 75, 0, 750); ZJetsVSPT30 ->Sumw2(); ZJetsVSPT30 ->SetLineColor(kGreen);   ZJetsVSPT30 ->SetLineWidth(2);
TH1D *ZJetsVSPT40  = new TH1D("ZJetsVSPT40" ,"", 75, 0, 750); ZJetsVSPT40 ->Sumw2(); ZJetsVSPT40 ->SetLineColor(kBlue);    ZJetsVSPT40 ->SetLineWidth(2);
TH1D *ZJetsVSPT50  = new TH1D("ZJetsVSPT50" ,"", 75, 0, 750); ZJetsVSPT50 ->Sumw2(); ZJetsVSPT50 ->SetLineColor(kYellow);  ZJetsVSPT50 ->SetLineWidth(2);
TH1D *ZJetsVSPT70  = new TH1D("ZJetsVSPT70" ,"", 75, 0, 750); ZJetsVSPT70 ->Sumw2(); ZJetsVSPT70 ->SetLineColor(kMagenta); ZJetsVSPT70 ->SetLineWidth(2);
TH1D *ZJetsVSPT100 = new TH1D("ZJetsVSPT100","", 75, 0, 750); ZJetsVSPT100->Sumw2(); ZJetsVSPT100->SetLineColor(kCyan);    ZJetsVSPT100->SetLineWidth(2);
TH1D *TTbarVSPT0   = new TH1D("TTbarVSPT0"  ,"", 75, 0, 750); TTbarVSPT0  ->Sumw2(); TTbarVSPT0  ->SetLineColor(kBlack);   TTbarVSPT0  ->SetLineWidth(2);
TH1D *TTbarVSPT20  = new TH1D("TTbarVSPT20" ,"", 75, 0, 750); TTbarVSPT20 ->Sumw2(); TTbarVSPT20 ->SetLineColor(kRed);     TTbarVSPT20 ->SetLineWidth(2);
TH1D *TTbarVSPT30  = new TH1D("TTbarVSPT30" ,"", 75, 0, 750); TTbarVSPT30 ->Sumw2(); TTbarVSPT30 ->SetLineColor(kGreen);   TTbarVSPT30 ->SetLineWidth(2);
TH1D *TTbarVSPT40  = new TH1D("TTbarVSPT40" ,"", 75, 0, 750); TTbarVSPT40 ->Sumw2(); TTbarVSPT40 ->SetLineColor(kBlue);    TTbarVSPT40 ->SetLineWidth(2);
TH1D *TTbarVSPT50  = new TH1D("TTbarVSPT50" ,"", 75, 0, 750); TTbarVSPT50 ->Sumw2(); TTbarVSPT50 ->SetLineColor(kYellow);  TTbarVSPT50 ->SetLineWidth(2);
TH1D *TTbarVSPT70  = new TH1D("TTbarVSPT70" ,"", 75, 0, 750); TTbarVSPT70 ->Sumw2(); TTbarVSPT70 ->SetLineColor(kMagenta); TTbarVSPT70 ->SetLineWidth(2);
TH1D *TTbarVSPT100 = new TH1D("TTbarVSPT100","", 75, 0, 750); TTbarVSPT100->Sumw2(); TTbarVSPT100->SetLineColor(kCyan);    TTbarVSPT100->SetLineWidth(2);
//small MinDPhi
TH1D *DataQCDVSPT0   = new TH1D("DataQCDVSPT0"  ,"", 75, 0, 750); DataQCDVSPT0  ->Sumw2(); DataQCDVSPT0  ->SetLineColor(kBlack);  DataQCDVSPT0  ->SetLineWidth(2);
TH1D *DataQCDVSPT20  = new TH1D("DataQCDVSPT20" ,"", 75, 0, 750); DataQCDVSPT20 ->Sumw2(); DataQCDVSPT20 ->SetLineColor(kRed);    DataQCDVSPT20 ->SetLineWidth(2);
TH1D *DataQCDVSPT30  = new TH1D("DataQCDVSPT30" ,"", 75, 0, 750); DataQCDVSPT30 ->Sumw2(); DataQCDVSPT30 ->SetLineColor(kGreen);  DataQCDVSPT30 ->SetLineWidth(2);
TH1D *DataQCDVSPT40  = new TH1D("DataQCDVSPT40" ,"", 75, 0, 750); DataQCDVSPT40 ->Sumw2(); DataQCDVSPT40 ->SetLineColor(kBlue);   DataQCDVSPT40 ->SetLineWidth(2);
TH1D *DataQCDVSPT50  = new TH1D("DataQCDVSPT50" ,"", 75, 0, 750); DataQCDVSPT50 ->Sumw2(); DataQCDVSPT50 ->SetLineColor(kYellow); DataQCDVSPT50 ->SetLineWidth(2);
TH1D *DataQCDVSPT70  = new TH1D("DataQCDVSPT70" ,"", 75, 0, 750); DataQCDVSPT70 ->Sumw2(); DataQCDVSPT70 ->SetLineColor(kMagenta);DataQCDVSPT70 ->SetLineWidth(2);
TH1D *DataQCDVSPT100 = new TH1D("DataQCDVSPT100","", 75, 0, 750); DataQCDVSPT100->Sumw2(); DataQCDVSPT100->SetLineColor(kCyan);   DataQCDVSPT100->SetLineWidth(2);
//1lep,0b
TH1D *DataWVSPT0   = new TH1D("DataWVSPT0"  ,"", 75, 0, 750); DataWVSPT0  ->Sumw2(); DataWVSPT0  ->SetLineColor(kBlack);  DataWVSPT0  ->SetLineWidth(2);
TH1D *DataWVSPT20  = new TH1D("DataWVSPT20" ,"", 75, 0, 750); DataWVSPT20 ->Sumw2(); DataWVSPT20 ->SetLineColor(kRed);    DataWVSPT20 ->SetLineWidth(2);
TH1D *DataWVSPT30  = new TH1D("DataWVSPT30" ,"", 75, 0, 750); DataWVSPT30 ->Sumw2(); DataWVSPT30 ->SetLineColor(kGreen);  DataWVSPT30 ->SetLineWidth(2);
TH1D *DataWVSPT40  = new TH1D("DataWVSPT40" ,"", 75, 0, 750); DataWVSPT40 ->Sumw2(); DataWVSPT40 ->SetLineColor(kBlue);   DataWVSPT40 ->SetLineWidth(2);
TH1D *DataWVSPT50  = new TH1D("DataWVSPT50" ,"", 75, 0, 750); DataWVSPT50 ->Sumw2(); DataWVSPT50 ->SetLineColor(kYellow); DataWVSPT50 ->SetLineWidth(2);
TH1D *DataWVSPT70  = new TH1D("DataWVSPT70" ,"", 75, 0, 750); DataWVSPT70 ->Sumw2(); DataWVSPT70 ->SetLineColor(kMagenta);DataWVSPT70 ->SetLineWidth(2);
TH1D *DataWVSPT100 = new TH1D("DataWVSPT100","", 75, 0, 750); DataWVSPT100->Sumw2(); DataWVSPT100->SetLineColor(kCyan);   DataWVSPT100->SetLineWidth(2);
//1lep,1b 
TH1D *DataTopVSPT0   = new TH1D("DataTopVSPT0"  ,"", 75, 0, 750); DataTopVSPT0  ->Sumw2(); DataTopVSPT0  ->SetLineColor(kBlack);  DataTopVSPT0  ->SetLineWidth(2);
TH1D *DataTopVSPT20  = new TH1D("DataTopVSPT20" ,"", 75, 0, 750); DataTopVSPT20 ->Sumw2(); DataTopVSPT20 ->SetLineColor(kRed);    DataTopVSPT20 ->SetLineWidth(2);
TH1D *DataTopVSPT30  = new TH1D("DataTopVSPT30" ,"", 75, 0, 750); DataTopVSPT30 ->Sumw2(); DataTopVSPT30 ->SetLineColor(kGreen);  DataTopVSPT30 ->SetLineWidth(2);
TH1D *DataTopVSPT40  = new TH1D("DataTopVSPT40" ,"", 75, 0, 750); DataTopVSPT40 ->Sumw2(); DataTopVSPT40 ->SetLineColor(kBlue);   DataTopVSPT40 ->SetLineWidth(2);
TH1D *DataTopVSPT50  = new TH1D("DataTopVSPT50" ,"", 75, 0, 750); DataTopVSPT50 ->Sumw2(); DataTopVSPT50 ->SetLineColor(kYellow); DataTopVSPT50 ->SetLineWidth(2);
TH1D *DataTopVSPT70  = new TH1D("DataTopVSPT70" ,"", 75, 0, 750); DataTopVSPT70 ->Sumw2(); DataTopVSPT70 ->SetLineColor(kMagenta);DataTopVSPT70 ->SetLineWidth(2);
TH1D *DataTopVSPT100 = new TH1D("DataTopVSPT100","", 75, 0, 750); DataTopVSPT100->Sumw2(); DataTopVSPT100->SetLineColor(kCyan);   DataTopVSPT100->SetLineWidth(2);
//highDeltaPhi,0b 
TH1D *DataZVSPT0   = new TH1D("DataZVSPT0"  ,"", 75, 0, 750); DataZVSPT0  ->Sumw2(); DataZVSPT0  ->SetLineColor(kBlack);  DataZVSPT0  ->SetLineWidth(2);
TH1D *DataZVSPT20  = new TH1D("DataZVSPT20" ,"", 75, 0, 750); DataZVSPT20 ->Sumw2(); DataZVSPT20 ->SetLineColor(kRed);    DataZVSPT20 ->SetLineWidth(2);
TH1D *DataZVSPT30  = new TH1D("DataZVSPT30" ,"", 75, 0, 750); DataZVSPT30 ->Sumw2(); DataZVSPT30 ->SetLineColor(kGreen);  DataZVSPT30 ->SetLineWidth(2);
TH1D *DataZVSPT40  = new TH1D("DataZVSPT40" ,"", 75, 0, 750); DataZVSPT40 ->Sumw2(); DataZVSPT40 ->SetLineColor(kBlue);   DataZVSPT40 ->SetLineWidth(2);
TH1D *DataZVSPT50  = new TH1D("DataZVSPT50" ,"", 75, 0, 750); DataZVSPT50 ->Sumw2(); DataZVSPT50 ->SetLineColor(kYellow); DataZVSPT50 ->SetLineWidth(2);
TH1D *DataZVSPT70  = new TH1D("DataZVSPT70" ,"", 75, 0, 750); DataZVSPT70 ->Sumw2(); DataZVSPT70 ->SetLineColor(kMagenta);DataZVSPT70 ->SetLineWidth(2);
TH1D *DataZVSPT100 = new TH1D("DataZVSPT100","", 75, 0, 750); DataZVSPT100->Sumw2(); DataZVSPT100->SetLineColor(kCyan);   DataZVSPT100->SetLineWidth(2);

//fill histograms
gROOT->cd();
cout << "LM6VSPT0" << endl;
LM6->Draw("misc.MT2>>LM6VSPT0",  "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=0&&misc.Vectorsumpt<20","goff");
cout << "LM6VSPT20" << endl;
LM6->Draw("misc.MT2>>LM6VSPT20", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=20&&misc.Vectorsumpt<30","goff");
cout << "LM6VSPT30" << endl;
LM6->Draw("misc.MT2>>LM6VSPT30", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=30&&misc.Vectorsumpt<40","goff");
cout << "LM6VSPT40" << endl;
LM6->Draw("misc.MT2>>LM6VSPT40", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=40&&misc.Vectorsumpt<50","goff");
cout << "LM6VSPT50" << endl;
LM6->Draw("misc.MT2>>LM6VSPT50", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=50&&misc.Vectorsumpt<70","goff");
cout << "LM6VSPT70" << endl;
LM6->Draw("misc.MT2>>LM6VSPT70", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=70&&misc.Vectorsumpt<100","goff");
cout << "LM6VSPT100" << endl;
LM6->Draw("misc.MT2>>LM6VSPT100","NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=100","goff");

cout << "WJetsVSPT0" << endl;
WJets->Draw("misc.MT2>>WJetsVSPT0",  "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=0&&misc.Vectorsumpt<20","goff");
cout << "WJetsVSPT20" << endl;
WJets->Draw("misc.MT2>>WJetsVSPT20", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=20&&misc.Vectorsumpt<30","goff");
cout << "WJetsVSPT30" << endl;
WJets->Draw("misc.MT2>>WJetsVSPT30", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=30&&misc.Vectorsumpt<40","goff");
cout << "WJetsVSPT40" << endl;
WJets->Draw("misc.MT2>>WJetsVSPT40", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=40&&misc.Vectorsumpt<50","goff");
cout << "WJetsVSPT50" << endl;
WJets->Draw("misc.MT2>>WJetsVSPT50", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=50&&misc.Vectorsumpt<70","goff");
cout << "WJetsVSPT70" << endl;
WJets->Draw("misc.MT2>>WJetsVSPT70", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=70&&misc.Vectorsumpt<100","goff");
cout << "WJetsVSPT100" << endl;
WJets->Draw("misc.MT2>>WJetsVSPT100","NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=100","goff");

cout << "ZJetsVSPT0" << endl;
ZJets->Draw("misc.MT2>>ZJetsVSPT0",  "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=0&&misc.Vectorsumpt<20","goff");
cout << "ZJetsVSPT20" << endl;
ZJets->Draw("misc.MT2>>ZJetsVSPT20", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=20&&misc.Vectorsumpt<30","goff");
cout << "ZJetsVSPT30" << endl;
ZJets->Draw("misc.MT2>>ZJetsVSPT30", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=30&&misc.Vectorsumpt<40","goff");
cout << "ZJetsVSPT40" << endl;
ZJets->Draw("misc.MT2>>ZJetsVSPT40", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=40&&misc.Vectorsumpt<50","goff");
cout << "ZJetsVSPT50" << endl;
ZJets->Draw("misc.MT2>>ZJetsVSPT50", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=50&&misc.Vectorsumpt<70","goff");
cout << "ZJetsVSPT70" << endl;
ZJets->Draw("misc.MT2>>ZJetsVSPT70", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=70&&misc.Vectorsumpt<100","goff");
cout << "ZJetsVSPT100" << endl;
ZJets->Draw("misc.MT2>>ZJetsVSPT100","NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=100","goff");

cout << "TTbarVSPT0" << endl;
TTbar->Draw("misc.MT2>>TTbarVSPT0",  "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=0&&misc.Vectorsumpt<20","goff");
cout << "TTbarVSPT20" << endl;
TTbar->Draw("misc.MT2>>TTbarVSPT20", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=20&&misc.Vectorsumpt<30","goff");
cout << "TTbarVSPT30" << endl;
TTbar->Draw("misc.MT2>>TTbarVSPT30", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=30&&misc.Vectorsumpt<40","goff");
cout << "TTbarVSPT40" << endl;
TTbar->Draw("misc.MT2>>TTbarVSPT40", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=40&&misc.Vectorsumpt<50","goff");
cout << "TTbarVSPT50" << endl;
TTbar->Draw("misc.MT2>>TTbarVSPT50", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=50&&misc.Vectorsumpt<70","goff");
cout << "TTbarVSPT70" << endl;
TTbar->Draw("misc.MT2>>TTbarVSPT70", "NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=70&&misc.Vectorsumpt<100","goff");
cout << "TTbarVSPT100" << endl;
TTbar->Draw("misc.MT2>>TTbarVSPT100","NJetsIDLoose40>=2&&misc.HT>450&&misc.MET>30&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.Vectorsumpt>=100","goff");

cout << "DataQCDVSPT0" << endl;
Data->Draw("misc.MT2>>DataQCDVSPT0",  "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40<0.2&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=0&&misc.Vectorsumpt<20","goff");
cout << "DataQCDVSPT20" << endl;
Data->Draw("misc.MT2>>DataQCDVSPT20", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40<0.2&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=20&&misc.Vectorsumpt<30","goff");
cout << "DataQCDVSPT30" << endl;
Data->Draw("misc.MT2>>DataQCDVSPT30", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40<0.2&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=30&&misc.Vectorsumpt<40","goff");
cout << "DataQCDVSPT40" << endl;
Data->Draw("misc.MT2>>DataQCDVSPT40", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40<0.2&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=40&&misc.Vectorsumpt<50","goff");
cout << "DataQCDVSPT50" << endl;
Data->Draw("misc.MT2>>DataQCDVSPT50", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40<0.2&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=50&&misc.Vectorsumpt<70","goff");
cout << "DataQCDVSPT70" << endl;
Data->Draw("misc.MT2>>DataQCDVSPT70", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40<0.2&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=70&&misc.Vectorsumpt<100","goff");
cout << "DataQCDVSPT100" << endl;
Data->Draw("misc.MT2>>DataQCDVSPT100","NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40<0.2&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=100","goff");

cout << "DataZVSPT0" << endl;
Data->Draw("misc.MT2>>DataZVSPT0",  "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=0&&misc.Vectorsumpt<20","goff");
cout << "DataZVSPT20" << endl;
Data->Draw("misc.MT2>>DataZVSPT20", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=20&&misc.Vectorsumpt<30","goff");
cout << "DataZVSPT30" << endl;
Data->Draw("misc.MT2>>DataZVSPT30", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=30&&misc.Vectorsumpt<40","goff");
cout << "DataZVSPT40" << endl;
Data->Draw("misc.MT2>>DataZVSPT40", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=40&&misc.Vectorsumpt<50","goff");
cout << "DataZVSPT50" << endl;
Data->Draw("misc.MT2>>DataZVSPT50", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=50&&misc.Vectorsumpt<70","goff");
cout << "DataZVSPT70" << endl;
Data->Draw("misc.MT2>>DataZVSPT70", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=70&&misc.Vectorsumpt<100","goff");
cout << "DataZVSPT100" << endl;
Data->Draw("misc.MT2>>DataZVSPT100","NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&NEles==0&&NMuons==0&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=100","goff");

cout << "DataWVSPT0" << endl;
Data->Draw("misc.MT2>>DataWVSPT0",  "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=0&&misc.Vectorsumpt<20","goff");
cout << "DataWVSPT20" << endl;
Data->Draw("misc.MT2>>DataWVSPT20", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=20&&misc.Vectorsumpt<30","goff");
cout << "DataWVSPT30" << endl;
Data->Draw("misc.MT2>>DataWVSPT30", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=30&&misc.Vectorsumpt<40","goff");
cout << "DataWVSPT40" << endl;
Data->Draw("misc.MT2>>DataWVSPT40", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=40&&misc.Vectorsumpt<50","goff");
cout << "DataWVSPT50" << endl;
Data->Draw("misc.MT2>>DataWVSPT50", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=50&&misc.Vectorsumpt<70","goff");
cout << "DataWVSPT70" << endl;
Data->Draw("misc.MT2>>DataWVSPT70", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=70&&misc.Vectorsumpt<100","goff");
cout << "DataWVSPT100" << endl;
Data->Draw("misc.MT2>>DataWVSPT100","NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM==0&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=100","goff");

cout << "DataWVSPT0" << endl;
Data->Draw("misc.MT2>>DataTopVSPT0",  "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM>=1&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=0&&misc.Vectorsumpt<20","goff");
cout << "DataTopVSPT20" << endl;
Data->Draw("misc.MT2>>DataTopVSPT20", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM>=1&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=20&&misc.Vectorsumpt<30","goff");
cout << "DataTopVSPT30" << endl;
Data->Draw("misc.MT2>>DataTopVSPT30", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM>=1&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=30&&misc.Vectorsumpt<40","goff");
cout << "DataTopVSPT40" << endl;
Data->Draw("misc.MT2>>DataTopVSPT40", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM>=1&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=40&&misc.Vectorsumpt<50","goff");
cout << "DataTopVSPT50" << endl;
Data->Draw("misc.MT2>>DataTopVSPT50", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM>=1&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=50&&misc.Vectorsumpt<70","goff");
cout << "DataTopVSPT70" << endl;
Data->Draw("misc.MT2>>DataTopVSPT70", "NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM>=1&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=70&&misc.Vectorsumpt<100","goff");
cout << "DataTopVSPT100" << endl;
Data->Draw("misc.MT2>>DataTopVSPT100","NJetsIDLoose40>=2&&misc.HT>750&&misc.MET>10&&(NEles+NMuons)==1&&NTausIDLoose3Hits==0&&misc.MinMetJetDPhi4Pt40>0.5&&NBJets40CSVM>=1&&misc.PassJet40ID ==1&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(misc.MET>30||misc.MET/misc.CaloMETRaw<=2.)&&misc.Vectorsumpt>=100","goff");

TLegend *leg = new TLegend(.6,.6,.90,.90);
leg->SetName("leg");
leg -> SetFillColor(0);
leg -> SetBorderSize(0);
leg->AddEntry(LM6VSPT0,"0 #leq VSPT < 20 GeV", "l");
leg->AddEntry(LM6VSPT20,"20 #leq VSPT < 30 GeV", "l");
leg->AddEntry(LM6VSPT30,"30 #leq VSPT < 40 GeV", "l");
leg->AddEntry(LM6VSPT40,"40 #leq VSPT < 50 GeV", "l");
leg->AddEntry(LM6VSPT50,"50 #leq VSPT < 70 GeV", "l");
leg->AddEntry(LM6VSPT70,"70 #leq VSPT < 100 GeV", "l");
leg->AddEntry(LM6VSPT100,"100 GeV #leq VSPT", "l");

//store the file
TFile *file = new TFile("UTMfile.root","RECREATE");
file->cd();
LM6VSPT0  ->Write();
LM6VSPT20 ->Write();
LM6VSPT30 ->Write();
LM6VSPT40 ->Write();
LM6VSPT50 ->Write();
LM6VSPT70 ->Write();
LM6VSPT100->Write();
WJetsVSPT0  ->Write();
WJetsVSPT20 ->Write();
WJetsVSPT30 ->Write();
WJetsVSPT40 ->Write();
WJetsVSPT50 ->Write();
WJetsVSPT70 ->Write();
WJetsVSPT100->Write();
ZJetsVSPT0  ->Write();
ZJetsVSPT20 ->Write();
ZJetsVSPT30 ->Write();
ZJetsVSPT40 ->Write();
ZJetsVSPT50 ->Write();
ZJetsVSPT70 ->Write();
ZJetsVSPT100->Write();
TTbarVSPT0  ->Write();
TTbarVSPT20 ->Write();
TTbarVSPT30 ->Write();
TTbarVSPT40 ->Write();
TTbarVSPT50 ->Write();
TTbarVSPT70 ->Write();
TTbarVSPT100->Write();
DataQCDVSPT0  ->Write();
DataQCDVSPT20 ->Write();
DataQCDVSPT30 ->Write();
DataQCDVSPT40 ->Write();
DataQCDVSPT50 ->Write();
DataQCDVSPT70 ->Write();
DataQCDVSPT100->Write();
DataWVSPT0  ->Write();
DataWVSPT20 ->Write();
DataWVSPT30 ->Write();
DataWVSPT40 ->Write();
DataWVSPT50 ->Write();
DataWVSPT70 ->Write();
DataWVSPT100->Write();
DataZVSPT0  ->Write();
DataZVSPT20 ->Write();
DataZVSPT30 ->Write();
DataZVSPT40 ->Write();
DataZVSPT50 ->Write();
DataZVSPT70 ->Write();
DataZVSPT100->Write();
DataTopVSPT0  ->Write();
DataTopVSPT20 ->Write();
DataTopVSPT30 ->Write();
DataTopVSPT40 ->Write();
DataTopVSPT50 ->Write();
DataTopVSPT70 ->Write();
DataTopVSPT100->Write();
leg->Write();
file->Close();

cout << "File saved: " << file->GetName() << endl;

}