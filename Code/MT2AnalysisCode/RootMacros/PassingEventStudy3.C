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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path for MT2tree.hh

//run via root -l -b -q PassingEventStudy3.C++

using namespace std;

//test TOBTEC reweighting
//and observe its influence
//plotting was done interactively
void PassingEventStudy3(){

gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

//get reweighting map
//TFile *effmapfile = TFile::Open("/shome/haweber/TOBTECfiles/PassingEventsStudy2.root");
TFile *effmapfile = TFile::Open("/shome/haweber/TOBTECfiles/PassingEventsStudy2_NoMinDPhiNoLeptonVetoIncludingPhotons.root");
effmapfile->cd();
TEfficiency *effmap = (TEfficiency*)effmapfile->Get("SJetPt100Etaeff");

//get a test sample
TChain *c = new TChain("MassTree");
int cnum1 = c->Add("/shome/haweber/MT2Analysis/MT2trees/MT2_V02-03-02/20130717_8TeV_sms/skimmed_T2tt_400_150/SMS-T2tt-mStop-375to475-mLSP-0to375-8TeV-Pythia6Z-Summer12-START52-V9-FSIM-v1-2.root");

//define 'reweighted' and notreweighted ('plain') histograms
TString cut = "NJetsIDLoose40>0";
TH1D *MT2plain    = new TH1D("MT2plain",   "",50,  0, 800); MT2plain   ->Sumw2();
TH1D *MT2reweight = new TH1D("MT2reweight","",50,  0, 800); MT2reweight->Sumw2();
TH1D *METplain    = new TH1D("METplain",   "",50,  0, 800); METplain   ->Sumw2();
TH1D *METreweight = new TH1D("METreweight","",50,  0, 800); METreweight->Sumw2();
TH1D *plain       = new TH1D("plain",      "", 1,  0,   1); plain      ->Sumw2();
TH1D *reweight    = new TH1D("reweight",   "", 1,  0,   1); reweight   ->Sumw2();
TH1D *HTplain     = new TH1D("HTplain",    "",21,450,1500); HTplain    ->Sumw2();
TH1D *HTreweight  = new TH1D("HTreweight", "",21,450,1500); HTreweight ->Sumw2();

const int netabin = 12;
double etabin[netabin+1];
for(int n = 0; n<=netabin;++n) etabin[n] = (double)n*0.2;
const int nptbin = 14;
double ptbin[nptbin+1] = {20,30,40,50,75,100,125,150,175,200,250,300,400,500,1000};

TH1D *SJetPtplain     = new TH1D("SJetPtplain",     "",nptbin,  ptbin ); SJetPtplain     ->Sumw2();
TH1D *SJetPtreweight  = new TH1D("SJetPtreweight",  "",nptbin,  ptbin ); SJetPtreweight  ->Sumw2();
TH1D *SJetEtaplain    = new TH1D("SJetEtaplain",    "",netabin, etabin); SJetEtaplain    ->Sumw2();
TH1D *SJetEtareweight = new TH1D("SJetEtareweight", "",netabin, etabin); SJetEtareweight ->Sumw2();
TH1D *Jet1Ptplain     = new TH1D("Jet1Ptplain",     "",nptbin,  ptbin ); Jet1Ptplain     ->Sumw2();
TH1D *Jet1Ptreweight  = new TH1D("Jet1Ptreweight",  "",nptbin,  ptbin ); Jet1Ptreweight  ->Sumw2();
TH1D *Jet1Etaplain    = new TH1D("Jet1Etaplain",    "",netabin, etabin); Jet1Etaplain    ->Sumw2();
TH1D *Jet1Etareweight = new TH1D("Jet1Etareweight", "",netabin, etabin); Jet1Etareweight ->Sumw2();
TH1D *Jet2Ptplain     = new TH1D("Jet2Ptplain",     "",nptbin,  ptbin ); Jet2Ptplain     ->Sumw2();
TH1D *Jet2Ptreweight  = new TH1D("Jet2Ptreweight",  "",nptbin,  ptbin ); Jet2Ptreweight  ->Sumw2();
TH1D *Jet2Etaplain    = new TH1D("Jet2Etaplain",    "",netabin, etabin); Jet2Etaplain    ->Sumw2();
TH1D *Jet2Etareweight = new TH1D("Jet2Etareweight", "",netabin, etabin); Jet2Etareweight ->Sumw2();

gROOT->cd();
MT2tree* fMT2tree = new MT2tree();
c->SetBranchAddress("MT2tree", &fMT2tree);
Float_t fTOBTECTagger;// = new Float_t;
c->SetBranchAddress("TOBTECTagger", &fTOBTECTagger);
Long64_t nentries =  c->GetEntries();
Long64_t nbytes = 0, nb = 0;
int nev =0;
c->Draw(">>selList", cut);//cuts HERE
TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
c->SetEventList(myEvtList);
int counter=0; int skipped=0;
cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
//run over all events
while(myEvtList->GetEntry(counter++) !=-1){	
	int jentry = myEvtList->GetEntry(counter-1);
	nb =  c->GetEntry(jentry);   nbytes += nb;
	if ( counter % 500 == 0  )  cout << "+++ Proccessing event " << counter << " " << fTOBTECTagger << endl;

	//get 'selected jet' on which we will reweight
	float ptmax(-1), etamin(99); int nind(-1);
	vector<int> ninds; ninds.clear();
	int njets = 0;
	for(int n = 0; n<fMT2tree->NJets; ++n){
		if(!(fMT2tree->jet[n].isPFIDLoose))     continue;
		if(fMT2tree->jet[n].lv.Pt()<100)         continue;
		if(fabs(fMT2tree->jet[n].lv.Eta())>2.4) continue; 
		float deltaetamin = fabs(1.5-etamin);
		float deltaeta = fabs(1.5-fabs(fMT2tree->jet[n].lv.Eta()));
		float eta = fabs(fMT2tree->jet[n].lv.Eta());
		float pt = fMT2tree->jet[n].lv.Pt();
		if(     deltaeta< 0.2&&deltaetamin>=0.2) {ptmax = pt; etamin = eta; nind = n; ninds.push_back(n); }
		else if(deltaeta< 0.2&&deltaetamin< 0.2) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} ninds.push_back(n);  }
		else if(deltaeta>=0.2&&deltaetamin>=0.2){
		if(     deltaeta< 0.4&&deltaetamin>=0.4) {ptmax = pt; etamin = eta; nind = n;}
		else if(deltaeta< 0.4&&deltaetamin< 0.4) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }  
		else if(deltaeta>=0.4&&deltaetamin>=0.4){
		if(     deltaeta< 0.6&&deltaetamin>=0.6) {ptmax = pt; etamin = eta; nind = n;}
		else if(deltaeta< 0.6&&deltaetamin< 0.6) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }  
		else if(deltaeta>=0.6&&deltaetamin>=0.6){
		if(     deltaeta< 0.8&&deltaetamin>=0.8) {ptmax = pt; etamin = eta; nind = n;}
		else if(deltaeta< 0.8&&deltaetamin< 0.8) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }  
		else if(deltaeta>=0.8&&deltaetamin>=0.8){
		if(     deltaeta< 1.0&&deltaetamin>=1.0) {ptmax = pt; etamin = eta; nind = n;}
		else if(deltaeta< 1.0&&deltaetamin< 1.0) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }  
		else if(deltaeta>=1.0&&deltaetamin>=1.0){
		if(     deltaeta< 1.2&&deltaetamin>=1.2) {ptmax = pt; etamin = eta; nind = n;}
		else if(deltaeta< 1.2&&deltaetamin< 1.2) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }  
		else if(deltaeta>=1.2&&deltaetamin>=1.2){
		if(     deltaeta< 1.4&&deltaetamin>=1.4) {ptmax = pt; etamin = eta; nind = n;}
		else if(deltaeta< 1.4&&deltaetamin< 1.4) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }  
		else if(deltaeta>=1.4&&deltaetamin>=1.4){
		if(     deltaeta< 1.6&&deltaetamin>=1.6) {ptmax = pt; etamin = eta; nind = n;}
		else if(deltaeta< 1.6&&deltaetamin< 1.6) { if(pt>ptmax) {ptmax = pt; etamin = eta; nind = n;} }  }}}}}}}
	}
	if(nind<0) {++skipped; continue;}
	//get TOBTEC 'reweighting' efficiency
	int bin = effmap->FindFixBin(ptmax,etamin);
	float eff = effmap->GetEfficiency(bin);
	float efflow = effmap->GetEfficiencyErrorLow(bin);
	float effup = effmap->GetEfficiencyErrorUp(bin);
	if( counter % 500 == 0  ) cout << "pt " << ptmax << " eta " << etamin << " --> binnumber " << bin << " ==> eff " << eff << endl;
	//fill 'reweighted' and 'plain' histogram
	MT2plain   ->Fill(fMT2tree->misc.MT2);
	MT2reweight->Fill(fMT2tree->misc.MT2, eff);
	METplain   ->Fill(fMT2tree->misc.MET);
	METreweight->Fill(fMT2tree->misc.MET, eff);
	HTplain    ->Fill(fMT2tree->misc.HT);
	HTreweight ->Fill(fMT2tree->misc.HT, eff);
	plain   ->Fill(0.5);
	reweight->Fill(0.5, eff);
	SJetPtplain     ->Fill(     fMT2tree->jet[nind].lv.Pt()   );
	SJetPtreweight  ->Fill(     fMT2tree->jet[nind].lv.Pt()  , eff);
	SJetEtaplain    ->Fill(fabs(fMT2tree->jet[nind].lv.Eta()) );
	SJetEtareweight ->Fill(fabs(fMT2tree->jet[nind].lv.Eta()), eff);
	Jet1Ptplain     ->Fill(     fMT2tree->jet[ 0  ].lv.Pt()   );
	Jet1Ptreweight  ->Fill(     fMT2tree->jet[ 0  ].lv.Pt()  , eff);
	Jet1Etaplain    ->Fill(fabs(fMT2tree->jet[ 0  ].lv.Eta()) );
	Jet1Etareweight ->Fill(fabs(fMT2tree->jet[ 0  ].lv.Eta()), eff);
	Jet2Ptplain     ->Fill(     fMT2tree->jet[  1 ].lv.Pt()   );
	Jet2Ptreweight  ->Fill(     fMT2tree->jet[  1 ].lv.Pt()  , eff);
	Jet2Etaplain    ->Fill(fabs(fMT2tree->jet[  1 ].lv.Eta()) );
	Jet2Etareweight ->Fill(fabs(fMT2tree->jet[  1 ].lv.Eta()), eff);
}

cout << "the FastSim(reweighted) / FastSim(original) value is " << reweight->GetBinContent(1)/plain->GetBinContent(1) << endl;
cout << "the FastSim(reweighted) / FastSim(original) value is " << reweight->Integral()/plain->Integral() << endl;
cout << "the FastSim(reweighted) / FastSim(original) value is " << MT2reweight->Integral()/MT2plain->Integral() << endl;
cout << "the FastSim(reweighted) / FastSim(original) value is " << METreweight->Integral()/METplain->Integral() << endl;

//save all histograms
TFile *newfile = new TFile("/shome/haweber/TOBTECfiles/PassingEventsStudy3_NoMinDPhiNoLeptonVetoIncludingPhotons.root","RECREATE");
//TFile *newfile = new TFile("/shome/haweber/TOBTECfiles/PassingEventsStudy3_NoMinDPhiNoLeptonVetoIncludingPhotons_eff2.root","RECREATE");
//TFile *newfile = new TFile("/shome/haweber/TOBTECfiles/PassingEventsStudy3_multiplied.root","RECREATE");//if using ninds
newfile    ->cd();
MT2plain   ->Write();
MT2reweight->Write();
METplain   ->Write();
METreweight->Write();
HTplain    ->Write();
HTreweight ->Write();
plain      ->Write();
reweight   ->Write();
SJetPtplain    ->Write();
SJetPtreweight ->Write();
SJetEtaplain   ->Write();
SJetEtareweight->Write();
Jet1Ptplain    ->Write();
Jet1Ptreweight ->Write();
Jet1Etaplain   ->Write();
Jet1Etareweight->Write();
Jet2Ptplain    ->Write();
Jet2Ptreweight ->Write();
Jet2Etaplain   ->Write();
Jet2Etareweight->Write();
newfile    ->Close();
cout << "saved histos in " << newfile->GetName() << endl;

}
