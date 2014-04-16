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

//run via root -l quickanalysis.C

using namespace std;

//this one takes events of the lowest MT2 signal bin in 3-5j,1b,high HT
//and plots several standard variables, just to check if there is some kind of noise or QCD there
//nothing conclusive was found.
void quickanalysis(){

TChain *c = new TChain("MassTree");
int cnum1 = c->Add("~mangano/filterTOBTEC/codeHJ/CMSSW_5_3_11/myTest/MT2Code/output/bloc1/*.root");
//int cnum2 = c->Add("~mangano/filterTOBTEC/codeHJ/CMSSW_5_3_11/myTest/MT2Code/output/bloc2/*.root");

TString cutpass = "((misc.MT2>200&&misc.HT>=450&&misc.HT<750&&misc.MET>200&&((trigger.HLT_PFMET150_v2==1||trigger.HLT_PFMET150_v3==1||trigger.HLT_PFMET150_v4==1||trigger.HLT_PFMET150_v5==1||trigger.HLT_PFMET150_v6==1||trigger.HLT_PFMET150_v7==1)||(trigger.HLT_PFHT350_PFMET100_v3==1||trigger.HLT_PFHT350_PFMET100_v4==1||trigger.HLT_PFHT350_PFMET100_v5==1||trigger.HLT_PFHT350_PFMET100_v6==1||trigger.HLT_PFHT350_PFMET100_v7==1||trigger.HLT_PFNoPUHT350_PFMET100_v1==1||trigger.HLT_PFNoPUHT350_PFMET100_v3==1||trigger.HLT_PFNoPUHT350_PFMET100_v4==1)))||(misc.MT2>100&misc.HT>=750&&misc.MET>30&&(trigger.HLT_PFHT650_v5==1||trigger.HLT_PFHT650_v6==1||trigger.HLT_PFHT650_v7==1||trigger.HLT_PFHT650_v8==1||trigger.HLT_PFHT650_v9==1||trigger.HLT_PFNoPUHT650_v1==1||trigger.HLT_PFNoPUHT650_v3==1||trigger.HLT_PFNoPUHT650_v4==1)))&&NMuons==0 && NEles==0 && NTausIDLoose3Hits==0&&misc.Jet0Pass ==1&&misc.Jet1Pass ==1&&misc.SecondJPt  > 100&&misc.PassJet40ID ==1&&misc.MinMetJetDPhi4Pt40 >0.3&&misc.Vectorsumpt < 70&&misc.HBHENoiseFlag == 0&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&misc.MET/misc.CaloMETRaw<=2.&&ExtraBeamHaloFilter==0&&TOBTECTagger>=0&&TOBTECTagger<=8";


TString cutbin = "NJetsIDLoose40>=3&&NJetsIDLoose40<=5&&NBJets40CSVM==1&&misc.HT>=1200&&misc.MT2>=150&&misc.MT2<=180&&"+cutpass;

TH1D *hDeltaPhiStar = new TH1D("hDeltaPhiStar","",32,0,3.2); hDeltaPhiStar->Sumw2();
TH1D *hDeltaPhiStarTest = new TH1D("hDeltaPhiStarTest","",32,0,3.2); hDeltaPhiStarTest->Sumw2();
TH1D *hMinDPhiMHT = new TH1D("hMinDPhiMHT","",32,0,3.2); hMinDPhiMHT->Sumw2();
TH1D *hMinDPhiMET = new TH1D("hMinDPhiMET","",32,0,3.2); hMinDPhiMET->Sumw2();

TH1D *halphaT = new TH1D("halphaT","", 25, 0, 1.); halphaT->Sumw2();

TH1D *hDeltaPhiStarAllJets20 = new TH1D("hDeltaPhiStarAllJets20","",32,0,3.2); hDeltaPhiStarAllJets20->Sumw2();
TH1D *hMinDPhiMHTAllJets20 = new TH1D("hMinDPhiMHTAllJets20","",32,0,3.2); hMinDPhiMHTAllJets20->Sumw2();
TH1D *hMinDPhiMETAllJets20 = new TH1D("hMinDPhiMETAllJets20","",32,0,3.2); hMinDPhiMETAllJets20->Sumw2();

TH1D *hMETPhi = new TH1D("hMETPhi","",96,-3.2,6.4); hMETPhi->Sumw2();
TH1D *hNVertices = new TH1D("hNVertices","",50,0,50); hNVertices->Sumw2();
TH1D *hHemiDPhi = new TH1D("hHemiDPhi","",32,0,3.2); hHemiDPhi->Sumw2();
TH1D *hNJets20 = new TH1D("hNJets20","",10,0,10); hNJets20->Sumw2();
TH1D *hNJets40 = new TH1D("hNJets40","",10,0,10); hNJets40->Sumw2();
TH1D *hNJets20alleta = new TH1D("hNJets20alleta","",10,0,10); hNJets20alleta->Sumw2();
TH1D *hNJets40alleta = new TH1D("hNJets40alleta","",10,0,10); hNJets40alleta->Sumw2();
TH1D *hPseudoJetMetDPhi = new TH1D("hPseudoJetMetDPhi","",32,0,3.2); hPseudoJetMetDPhi->Sumw2();


gROOT->cd();
MT2tree* fMT2tree = new MT2tree();
c->SetBranchAddress("MT2tree", &fMT2tree);
Float_t fTOBTECTagger;// = new Float_t;
c->SetBranchAddress("TOBTECTagger", &fTOBTECTagger);
Long64_t nentries =  c->GetEntries();
Long64_t nbytes = 0, nb = 0;
int nev =0;

c->Draw(">>selList", cutbin);//cuts HERE
TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
c->SetEventList(myEvtList);
int counter=0;
cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
while(myEvtList->GetEntry(counter++) !=-1){	
	int jentry = myEvtList->GetEntry(counter-1);
	nb =  c->GetEntry(jentry);   nbytes += nb;
	if ( counter % 5 == 0  )  cout << "+++ Proccessing event " << counter << " " << fTOBTECTagger << endl;
	float mindphi(99.);
	float mindphiMET(99.);
	float mindphiMHT(99.);
	TLorentzVector MHTminusJet(0.,0.,0.,0.);
	TLorentzVector Jetlv(0.,0.,0.,0.);
	int count(0);
	for(int i = 0; i<=fMT2tree->NJets; ++i){
		if(!(fMT2tree->jet[i].isPFIDLoose))     continue;
		if(fMT2tree->jet[i].lv.Pt()<40)         continue;
		if(fabs(fMT2tree->jet[i].lv.Eta())>2.4) continue;
		++count;
		MHTminusJet.SetPtEtaPhiE(0.,0.,0.,0.);
		Jetlv.SetPtEtaPhiE(0.,0.,0.,0.);
		MHTminusJet = fMT2tree->MHT[0] + fMT2tree->jet[i].lv;
		Jetlv = fMT2tree->jet[i].lv;
		float deltaphi = fabs(MHTminusJet.DeltaPhi(Jetlv));
		if(deltaphi<mindphi) mindphi = deltaphi;
		deltaphi = fabs(Jetlv.DeltaPhi(fMT2tree->MHT[0]));
		if(deltaphi<mindphiMHT) mindphiMHT = deltaphi;
		deltaphi = fabs(Jetlv.DeltaPhi(fMT2tree->pfmet[0]));
		if(deltaphi<mindphiMET) mindphiMET = deltaphi;
		//cout << deltaphi << " " << mindphi << endl;
	}
	if(count<3||count>5) cout << "ERROR " << count << endl;
	//cout << fMT2tree->misc.Run<<":"<<fMT2tree->misc.LumiSection<<":"<<fMT2tree->misc.Event<<"," << endl;
	hDeltaPhiStar->Fill(mindphi);
	hMinDPhiMHT->Fill(mindphiMHT);
	hMinDPhiMET->Fill(mindphiMET);

	mindphi= (99.);
	mindphiMET = (99.);
	mindphiMHT = (99.);
	for(int i = 0; i<=fMT2tree->NJets; ++i){
		if(!(fMT2tree->jet[i].isPFIDLoose))     continue;
		if(fMT2tree->jet[i].lv.Pt()<20)         continue;
		if(fabs(fMT2tree->jet[i].lv.Eta())>2.4) continue;
		++count;
		MHTminusJet.SetPtEtaPhiE(0.,0.,0.,0.);
		Jetlv.SetPtEtaPhiE(0.,0.,0.,0.);
		MHTminusJet = fMT2tree->MHT[0] + fMT2tree->jet[i].lv;
		Jetlv = fMT2tree->jet[i].lv;
		float deltaphi = fabs(MHTminusJet.DeltaPhi(Jetlv));
		if(deltaphi<mindphi) mindphi = deltaphi;
		deltaphi = fabs(Jetlv.DeltaPhi(fMT2tree->MHT[0]));
		if(deltaphi<mindphiMHT) mindphiMHT = deltaphi;
		deltaphi = fabs(Jetlv.DeltaPhi(fMT2tree->pfmet[0]));
		if(deltaphi<mindphiMET) mindphiMET = deltaphi;
		//cout << deltaphi << " " << mindphi << endl;
	}
	hDeltaPhiStarAllJets20->Fill(mindphi);
	hMinDPhiMHTAllJets20->Fill(mindphiMHT);
	hMinDPhiMETAllJets20->Fill(mindphiMET);

	double alphaT = fMT2tree->hemi[1].AlphaT;
	if(alphaT>2.5) cout << "max alphaT " << alphaT << endl;
	if(alphaT<=0) cout << "min alphaT " << alphaT << endl;
	halphaT->Fill(alphaT);
	hMETPhi->Fill(fMT2tree->misc.METPhi);
	hNVertices->Fill(fMT2tree->pileUp.NVertices);
	hHemiDPhi->Fill(fMT2tree->hemi[0].dPhi);
	hNJets40alleta->Fill(fMT2tree->GetNjets(40.,10.,1));
	hNJets20alleta->Fill(fMT2tree->GetNjets(20.,10.,1));
	hNJets40->Fill(fMT2tree->GetNjets(40.,2.4,1));
	hNJets20->Fill(fMT2tree->GetNjets(20.,2.4,1));
	hDeltaPhiStarTest->Fill(fMT2tree->BiasedDPhi(1,40.,2.4));
	hPseudoJetMetDPhi->Fill(fMT2tree->PseudoJetMetDPhi());
	cout << "R:LS:E " << fMT2tree->misc.Run << ":"<<fMT2tree->misc.LumiSection<<":"<<fMT2tree->misc.Event<<" MT2 " << fMT2tree->misc.MT2 << " MET " << fMT2tree->misc.MET << " HT " << fMT2tree->misc.HT << " NJets " << fMT2tree->NJetsIDLoose40 << " TOBTEC " << fTOBTECTagger << " HO " << fMT2tree->misc.MET/fMT2tree->misc.CaloMETRaw << endl;
}

TFile *newfile = new TFile("DeltaPhiStar.root","RECREATE");
newfile->cd();
hDeltaPhiStar->Write();
hDeltaPhiStarTest->Write();
hMinDPhiMHT->Write();
hMinDPhiMET->Write();
hDeltaPhiStarAllJets20->Write();
hMinDPhiMHTAllJets20->Write();
hMinDPhiMETAllJets20->Write();
halphaT->Write();
hMETPhi->Write();
hNVertices->Write();
hHemiDPhi->Write();
hNJets20->Write();
hNJets40->Write();
hNJets20alleta->Write();
hNJets40alleta->Write();
hPseudoJetMetDPhi->Write();
newfile->Close();
cout << "written " << newfile->GetName() << endl;
}


