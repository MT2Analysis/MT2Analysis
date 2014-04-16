#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

//run via root -l -b -q TauEffHistos.C++
//just dump the numbers of the TauPOG into histograms
void TauEffHistos(){

const int neffbins = 16;
double effbins[neffbins+1] = {  15,  18,  20,  22,  24,  26,  28,  30,  34,  38,  42,  46,  50,  56,  62,  70,  80};
double effvalues[neffbins] = {0.51,0.54,0.56,0.58,0.59,0.60,0.60,0.60,0.60,0.60,0.60,0.61,0.60,0.59,0.61,0.60};

const int nfakebins = 11;
double fakebins[nfakebins+1] = {   20,   25,   30,   35,   40,   45,   50,   60,   75,  100,  140,  200};
double fakevalues[nfakebins] = {0.013,0.027,0.036,0.039,0.037,0.033,0.027,0.020,0.013,0.008,0.005};

const int nelebins = 2;
double elebins[nelebins+1] = {   0.,   1.4,2.5};
//double elevalues[nelebins] = {0.146,0.318};//loose
//double elevalues[nelebins] = {0.019,0.063};//medium
double elevalues[nelebins] = {0.011,0.036};//tight
//double elevalues[nelebins] = {0.005,0.031};//mva

const int nmuobins = 3;
double muobins[nmuobins+1] = {     0.,    1.2,    1.7,  2.5};
//double muovalues[nmuobins] = {0.0058,0.0037,0.0011};//loose
//double muovalues[nmuobins] = {0.00052,0.00040,0.00041};//medium
double muovalues[nmuobins] = {0.00024,0.00018,0.00015};//tight

TH1F *taueff  = new TH1F("HadTauTauEfficiency_mc", "HadTauTauEfficiency_mc", neffbins, effbins);
TH1F *jetfake = new TH1F("JetTauEfficiency_mc", "JetTauEfficiency_mc", nfakebins, fakebins);
TH1F *elefake = new TH1F("ElectronTauEfficiency_mc", "ElectronTauEfficiency_mc", nelebins, elebins);
TH1F *muofake = new TH1F("MuonTauEfficiency_mc", "MuonTauEfficiency_mc", nmuobins, muobins);

for(int i = 1; i<= neffbins; ++i) taueff ->SetBinContent(i,  effvalues[i-1]-0.075-0.10);//0.075 due to muorej, 0.10 due to elerej
for(int i = 1; i<=nfakebins; ++i) jetfake->SetBinContent(i, fakevalues[i-1]);
for(int i = 1; i<= nelebins; ++i) elefake->SetBinContent(i,  elevalues[i-1]);
for(int i = 1; i<= nmuobins; ++i) muofake->SetBinContent(i,  muovalues[i-1]);

TFile *file = new TFile("/shome/haweber/MT2Analysis_8TeV/Code/Efficiencies/TauEfficiencies_LooseTauCommonObject_fromTauPOGplots.root", "RECREATE");
file->cd();
taueff ->Write();
jetfake->Write();
elefake->Write();
muofake->Write();

cout << "TauEfficiencies_Histos saved in " << file->GetName() << endl;
}

