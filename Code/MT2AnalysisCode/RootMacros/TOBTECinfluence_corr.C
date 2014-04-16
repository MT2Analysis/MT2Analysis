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
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TMap.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <cmath>
#include <limits.h>
#include <utility>
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path of MT2tree.hh

//run via root -l -b -q TOBTECinfluence_corr.C++

using namespace std;

//this is not the newest code - the newest is TOBTECinfluence.C
//but has some features the TOBTECinfluence.C has not
//basically this code looks at most important variables for original RECO and RERECO with fixed tracking
//compares the relative difference with respect to the TOBTEC tagger (this feature is not there anymore in newest version TOBTECinfluence.C)
//and also compares absolute distributions of RECO and RERECO version with the TOBTEC tagger value
void TOBTECinfluence_corr() {

   //these are three different versions of Boris TOBTEC tagger - I actually forgot the number, which of these is our default
   int version = 0;//0,1,2

   gROOT->SetStyle("Plain");
   gStyle->SetPalette(1);

   //histograms - reco version
   TH2D *hMT2 = new TH2D("hMT2", "", 50,0,1000, 10,0,20); hMT2->GetXaxis()->SetTitle("M_{T2} (reco) [GeV]");  hMT2->GetYaxis()->SetTitle("tagger"); hMT2->Sumw2();
   TH2D *hHT = new TH2D("hHT", "", 100,0,2000, 10,0,20); hHT->GetXaxis()->SetTitle("H_{T} (reco) [GeV]");  hHT->GetYaxis()->SetTitle("tagger"); hHT->Sumw2();
   TH2D *hMET = new TH2D("hMET", "", 50,0,1000, 10,0,20); hMET->GetXaxis()->SetTitle("E_{T}^{miss} (reco) [GeV]");  hMET->GetYaxis()->SetTitle("tagger"); hMET->Sumw2();
   TH2D *hVSPT = new TH2D("hVSPT", "", 30,0,150, 10,0,20); hVSPT->GetXaxis()->SetTitle("VSPT (reco) [GeV]");  hVSPT->GetYaxis()->SetTitle("tagger"); hVSPT->Sumw2();
   TH2D *hMinDPhi = new TH2D("hMinDPhi", "", 32,0,3.2, 10,0,20); hMinDPhi->GetXaxis()->SetTitle("min#Delta#phi (reco)");  hMinDPhi->GetYaxis()->SetTitle("tagger"); hMinDPhi->Sumw2();
   TH2D *hJ1Pt = new TH2D("hJ1Pt", "", 50,0,1000, 10,0,20); hJ1Pt->GetXaxis()->SetTitle("j_{1}-p_{T} (reco) [GeV]");  hJ1Pt->GetYaxis()->SetTitle("tagger"); hJ1Pt->Sumw2();
   TH2D *hJ2Pt = new TH2D("hJ2Pt", "", 50,0,1000, 10,0,20); hJ2Pt->GetXaxis()->SetTitle("j_{2}-p_{T} (reco) [GeV]");  hJ2Pt->GetYaxis()->SetTitle("tagger"); hJ2Pt->Sumw2();
   TH2D *hNJets = new TH2D("hNJets", "", 10,0,10, 10,0,20); hNJets->GetXaxis()->SetTitle("# jets (reco)");  hNJets->GetYaxis()->SetTitle("tagger"); hNJets->Sumw2();
   TH2D *hNBJets = new TH2D("hNBJets", "", 5,0,5, 10,0,20); hNBJets->GetXaxis()->SetTitle("# b-jets (reco)");  hNBJets->GetYaxis()->SetTitle("tagger"); hNBJets->Sumw2();
   TH2D *hNLeps = new TH2D("hNLeps", "", 5,0,5, 10,0,20); hNLeps->GetXaxis()->SetTitle("# leptons (reco)");  hNLeps->GetYaxis()->SetTitle("tagger"); hNLeps->Sumw2();
   //histogram-  rereco version
   TH2D *hMT2RE = new TH2D("hMT2RE", "", 50,0,1000, 10,0,20); hMT2RE->GetXaxis()->SetTitle("M_{T2} (rereco) [GeV]");  hMT2RE->GetYaxis()->SetTitle("tagger"); hMT2RE->Sumw2();
   TH2D *hHTRE = new TH2D("hHTRE", "", 100,0,2000, 10,0,20); hHTRE->GetXaxis()->SetTitle("H_{T} (rereco) [GeV]");  hHTRE->GetYaxis()->SetTitle("tagger"); hHTRE->Sumw2();
   TH2D *hMETRE = new TH2D("hMETRE", "", 50,0,1000, 10,0,20); hMETRE->GetXaxis()->SetTitle("E_{T}^{miss} (rereco) [GeV]");  hMETRE->GetYaxis()->SetTitle("tagger"); hMETRE->Sumw2();
   TH2D *hVSPTRE = new TH2D("hVSPTRE", "", 30,0,150, 10,0,20); hVSPTRE->GetXaxis()->SetTitle("VSPT (rereco) [GeV]");  hVSPTRE->GetYaxis()->SetTitle("tagger"); hVSPTRE->Sumw2();
   TH2D *hMinDPhiRE = new TH2D("hMinDPhiRE", "", 32,0,3.2, 10,0,20); hMinDPhiRE->GetXaxis()->SetTitle("min#Delta#phi (rereco)");  hMinDPhiRE->GetYaxis()->SetTitle("tagger"); hMinDPhiRE->Sumw2();
   TH2D *hJ1PtRE = new TH2D("hJ1PtRE", "", 50,0,1000, 10,0,20); hJ1PtRE->GetXaxis()->SetTitle("j_{1}-p_{T} (rereco) [GeV]");  hJ1PtRE->GetYaxis()->SetTitle("tagger"); hJ1PtRE->Sumw2();
   TH2D *hJ2PtRE = new TH2D("hJ2PtRE", "", 50,0,1000, 10,0,20); hJ2PtRE->GetXaxis()->SetTitle("j_{2}-p_{T} (rereco) [GeV]");  hJ2PtRE->GetYaxis()->SetTitle("tagger"); hJ2PtRE->Sumw2();
   TH2D *hNJetsRE = new TH2D("hNJetsRE", "", 10,0,10, 10,0,20); hNJetsRE->GetXaxis()->SetTitle("# jets (rereco)");  hNJetsRE->GetYaxis()->SetTitle("tagger"); hNJetsRE->Sumw2();
   TH2D *hNBJetsRE = new TH2D("hNBJetsRE", "", 5,0,5, 10,0,20); hNBJetsRE->GetXaxis()->SetTitle("# b-jets (rereco)");  hNBJetsRE->GetYaxis()->SetTitle("tagger"); hNBJetsRE->Sumw2();
   TH2D *hNLepsRE = new TH2D("hNLepsRE", "", 5,0,5, 10,0,20); hNLepsRE->GetXaxis()->SetTitle("# leptons (rereco)");  hNLepsRE->GetYaxis()->SetTitle("tagger"); hNLepsRE->Sumw2();
   //histogram - relative difference (reco-rereco)/reco
   TH2D *hMT2Tag = new TH2D("hMT2Tag", "", 60,-3,3, 10,0,20); hMT2Tag->GetXaxis()->SetTitle("M_{T2} rel.change");  hMT2Tag->GetYaxis()->SetTitle("tagger"); hMT2Tag->Sumw2();
   TH2D *hHTTag = new TH2D("hHTTag", "", 60,-3,3, 10,0,20); hHTTag->GetXaxis()->SetTitle("H_{T} rel.change");  hHTTag->GetYaxis()->SetTitle("tagger"); hHTTag->Sumw2();
   TH2D *hMETTag = new TH2D("hMETTag", "", 60,-3,3, 10,0,20); hMETTag->GetXaxis()->SetTitle("E_{T}^{miss} rel.change");  hMETTag->GetYaxis()->SetTitle("tagger"); hMETTag->Sumw2();
   TH2D *hVSPTTag = new TH2D("hVSPTTag", "", 60,-3,3, 10,0,20); hVSPTTag->GetXaxis()->SetTitle("VSPT rel.change");  hVSPTTag->GetYaxis()->SetTitle("tagger"); hVSPTTag->Sumw2();
   TH2D *hMinDPhiTag = new TH2D("hMinDPhiTag", "", 60,-3,3, 10,0,20); hMinDPhiTag->GetXaxis()->SetTitle("min#Delta#phi rel.change");  hMinDPhiTag->GetYaxis()->SetTitle("tagger"); hMinDPhiTag->Sumw2();
   TH2D *hJ1PtTag = new TH2D("hJ1PtTag", "", 60,-3,3, 10,0,20); hJ1PtTag->GetXaxis()->SetTitle("j_{1}-p_{T} rel.change");  hJ1PtTag->GetYaxis()->SetTitle("tagger"); hJ1PtTag->Sumw2();
   TH2D *hJ2PtTag = new TH2D("hJ2PtTag", "", 60,-3,3, 10,0,20); hJ2PtTag->GetXaxis()->SetTitle("j_{2}-p_{T} rel.change");  hJ2PtTag->GetYaxis()->SetTitle("tagger"); hJ2PtTag->Sumw2();
   TH2D *hNJetsTag = new TH2D("hNJetsTag", "", 60,-3,3, 10,0,20); hNJetsTag->GetXaxis()->SetTitle("# jets rel.change");  hNJetsTag->GetYaxis()->SetTitle("tagger"); hNJetsTag->Sumw2();
   TH2D *hNBJetsTag = new TH2D("hNBJetsTag", "", 60,-3,3, 10,0,20); hNBJetsTag->GetXaxis()->SetTitle("# b-jets rel.change");  hNBJetsTag->GetYaxis()->SetTitle("tagger"); hNBJetsTag->Sumw2();
   TH2D *hNLepsTag = new TH2D("hNLepsTag", "", 60,-3,3, 10,0,20); hNLepsTag->GetXaxis()->SetTitle("# leptons rel.change");  hNLepsTag->GetYaxis()->SetTitle("tagger"); hNLepsTag->Sumw2();

   //inputfiles
   TFile *recofile = TFile::Open("/shome/haweber/TOBTECfiles/MT2trees_TOBTEC_RECO_MarcsFiles.root");
   TFile *rerecofile = TFile::Open("/shome/haweber/TOBTECfiles/MT2trees_TOBTEC_RERECO_MarcsFiles.root");
   TFile *tobtecfilterfile;
   if(version==0) tobtecfilterfile = TFile::Open("/shome/haweber/TOBTECfiles/tobtec.root");
   if(version==1) tobtecfilterfile = TFile::Open("/shome/haweber/TOBTECfiles/tobtec_v2b.root");
   if(version==2) tobtecfilterfile = TFile::Open("/shome/haweber/TOBTECfiles/tobtec_v2c.root");

   TTree *recotree = (TTree*)recofile->Get("MassTree");
   TTree *rerecotree = (TTree*)rerecofile->Get("MassTree");
   TTree *tobtectree = (TTree*)tobtecfilterfile->Get("tobtec");

   MT2tree* fMT2treeRECO = new MT2tree();
   recotree->SetBranchAddress("MT2tree", &fMT2treeRECO);
   MT2tree* fMT2treeRERECO = new MT2tree();
   rerecotree->SetBranchAddress("MT2tree", &fMT2treeRERECO);
   UInt_t evt;
   Int_t LS, run;
   Float_t borisRECO, borisRERECO;
   Bool_t  RA2bRECO, RA2bRERECO;
   //taggers to test are the first ones, match with MT2tree using run:lumisection:event
   tobtectree->SetBranchAddress("BorisFilterRECO", &borisRECO);
   tobtectree->SetBranchAddress("BorisFilterRERECO", &borisRERECO);
   tobtectree->SetBranchAddress("RA2bFilterRECO", &RA2bRECO);
   tobtectree->SetBranchAddress("RA2bFilterRERECO", &RA2bRERECO);
   tobtectree->SetBranchAddress("Run", &run);
   tobtectree->SetBranchAddress("LumiSection", &LS);
   tobtectree->SetBranchAddress("Event", &evt);

   int counter =0, counter2;
   int cgf(0),cgp(0),cbf(0),cbp(0);
   Long64_t nentries =  recotree->GetEntries();
   Long64_t nentriesRR =  rerecotree->GetEntries();
   Long64_t nentriesTT =  tobtectree->GetEntries();
   Long64_t nbytesRECO = 0, nbRECO = 0;
   Long64_t nbytesRERECO = 0, nbRERECO = 0;
   Long64_t nbytesTT = 0, nbTT = 0;
   int nev =0;
   for (Long64_t jentryRECO=0; jentryRECO<nentries;jentryRECO++) {
	nbRECO =  recotree->GetEntry(jentryRECO);   nbytesRECO += nbRECO;
	recotree->SetBranchAddress("MT2tree", &fMT2treeRECO);
	double MT2RECO = fMT2treeRECO->misc.MT2;
	double HTRECO = fMT2treeRECO->misc.HT;
	double METRECO = fMT2treeRECO->misc.MET;
	double METPhiRECO = fMT2treeRECO->misc.METPhi;
	double MHTRECO = fMT2treeRECO->MHT[0].Pt();
	int NLepsRECO = fMT2treeRECO->NEles+fMT2treeRECO->NMuons+fMT2treeRECO->NTausIDLoose3Hits;
	int NJetsRECO = fMT2treeRECO->NJetsIDLoose40;
	int NJets20RECO = fMT2treeRECO->NJetsIDLoose;
	int NBJetsRECO = fMT2treeRECO->NBJets40CSVM;
	double VSPTRECO = fMT2treeRECO->GetMHTminusMET(0, 20, 2.4, true);
	double MinDPhiRECO = fMT2treeRECO->misc.MinMetJetDPhi4Pt40;
	double J1PtRECO = fMT2treeRECO->misc.LeadingJPt;
	double J2PtRECO = fMT2treeRECO->misc.SecondJPt;
	int eventRECO = fMT2treeRECO->misc.Event;
	int runRECO = fMT2treeRECO->misc.Run;
	int LSRECO = fMT2treeRECO->misc.LumiSection;
	for (Long64_t jentryRERECO=0; jentryRERECO<nentries;jentryRERECO++) {
		nbRERECO =  rerecotree->GetEntry(jentryRECO);   nbytesRERECO += nbRERECO;
		rerecotree->SetBranchAddress("MT2tree", &fMT2treeRERECO);
		int eventRERECO = fMT2treeRERECO->misc.Event;
		int runRERECO = fMT2treeRERECO->misc.Run;
		int LSRERECO = fMT2treeRERECO->misc.LumiSection;
		if(eventRECO!=eventRERECO || LSRECO!=LSRERECO || runRECO!=runRERECO){//match rereco tree to reco tree
			nbytesRERECO -= nbRERECO;
			nbRERECO =  rerecotree->GetEntry(jentryRERECO);   nbytesRERECO += nbRERECO;
			rerecotree->SetBranchAddress("MT2tree", &fMT2treeRERECO);
		}
		double MT2RERECO = fMT2treeRERECO->misc.MT2;
		double HTRERECO = fMT2treeRERECO->misc.HT;
		double METRERECO = fMT2treeRERECO->misc.MET;
		double METPhiRERECO = fMT2treeRERECO->misc.METPhi;
		double MHTRERECO = fMT2treeRERECO->MHT[0].Pt();
		int NLepsRERECO = fMT2treeRERECO->NEles+fMT2treeRERECO->NMuons+fMT2treeRERECO->NTausIDLoose3Hits;
		int NJetsRERECO = fMT2treeRERECO->NJetsIDLoose40;
		int NJets20RERECO = fMT2treeRERECO->NJetsIDLoose;
		int NBJetsRERECO = fMT2treeRERECO->NBJets40CSVM;
		double VSPTRERECO = fMT2treeRERECO->GetMHTminusMET(0, 20, 2.4, true);
		double MinDPhiRERECO = fMT2treeRERECO->misc.MinMetJetDPhi4Pt40;
		double J1PtRERECO = fMT2treeRERECO->misc.LeadingJPt;
		double J2PtRERECO = fMT2treeRERECO->misc.SecondJPt;
		eventRERECO = fMT2treeRERECO->misc.Event;
		runRERECO = fMT2treeRERECO->misc.Run;
		LSRERECO = fMT2treeRERECO->misc.LumiSection;
		if(eventRECO!=eventRERECO || LSRECO!=LSRERECO || runRECO!=runRERECO) continue;

		for (Long64_t jentryTT=0; jentryTT<nentries;jentryTT++) {
			nbTT =  tobtectree->GetEntry(jentryRECO);   nbytesTT += nbTT;
			tobtectree->SetBranchAddress("Run", &run);
			tobtectree->SetBranchAddress("LumiSection", &LS);
			tobtectree->SetBranchAddress("Event", &evt);
			if(evt!=eventRECO || LS!=LSRECO || run!=runRERECO){//match tobtectree to reco tree
				nbytesTT -= nbTT;
				nbTT =  tobtectree->GetEntry(jentryTT);   nbytesTT += nbTT;
				tobtectree->SetBranchAddress("Run", &run);
				tobtectree->SetBranchAddress("LumiSection", &LS);
				tobtectree->SetBranchAddress("Event", &evt);
			}
			if(evt!=eventRECO || LS!=LSRECO || run!=runRERECO) continue;
			tobtectree->SetBranchAddress("BorisFilterRECO", &borisRECO);
			tobtectree->SetBranchAddress("BorisFilterRERECO", &borisRERECO);
			tobtectree->SetBranchAddress("RA2bFilterRECO", &RA2bRECO);
			tobtectree->SetBranchAddress("RA2bFilterRERECO", &RA2bRERECO);
			++counter;
			//fill histograms
			hMT2->Fill(MT2RECO,borisRECO);
			hHT->Fill(HTRECO,borisRECO);
			hMET->Fill(METRECO,borisRECO);
			hVSPT->Fill(VSPTRECO,borisRECO);
			hMinDPhi->Fill(MinDPhiRECO,borisRECO);
			hJ1Pt->Fill(J1PtRECO,borisRECO);
			hJ2Pt->Fill(J2PtRECO,borisRECO);
			hNJets->Fill(NJetsRECO,borisRECO);
			hNBJets->Fill(NBJetsRECO,borisRECO);
			hNLeps->Fill(NLepsRECO,borisRECO);
			hMT2RE->Fill(MT2RERECO,borisRECO);
			hHTRE->Fill(HTRERECO,borisRECO);
			hMETRE->Fill(METRERECO,borisRECO);
			hVSPTRE->Fill(VSPTRERECO,borisRECO);
			hMinDPhiRE->Fill(MinDPhiRERECO,borisRECO);
			hJ1PtRE->Fill(J1PtRERECO,borisRECO);
			hJ2PtRE->Fill(J2PtRERECO,borisRECO);
			hNJetsRE->Fill(NJetsRERECO,borisRECO);
			hNBJetsRE->Fill(NBJetsRERECO,borisRECO);
			hNLepsRE->Fill(NLepsRERECO,borisRECO);
			if(MT2RECO!=0) hMT2Tag->Fill((MT2RECO-MT2RERECO)/MT2RECO,borisRECO);
			if(HTRECO!=0) hHTTag->Fill((HTRECO-HTRERECO)/HTRECO,borisRECO);
			if(METRECO!=0) hMETTag->Fill((METRECO-METRERECO)/METRECO,borisRECO);
			if(VSPTRECO!=0) hVSPTTag->Fill((VSPTRECO-VSPTRERECO)/VSPTRECO,borisRECO);
			if(MinDPhiRECO!=0) hMinDPhiTag->Fill((MinDPhiRECO-MinDPhiRERECO)/MinDPhiRECO,borisRECO);
			if(J1PtRECO!=0) hJ1PtTag->Fill((J1PtRECO-J1PtRERECO)/J1PtRECO,borisRECO);
			if(J2PtRECO!=0) hJ2PtTag->Fill((J2PtRECO-J2PtRERECO)/J2PtRECO,borisRECO);
			if(NJetsRECO!=0) hNJetsTag->Fill((NJetsRECO-NJetsRERECO)/NJetsRECO,borisRECO);
			if(NBJetsRECO!=0) hNBJetsTag->Fill((NBJetsRECO-NBJetsRERECO)/NBJetsRECO,borisRECO);
			if(NLepsRECO!=0) hNLepsTag->Fill((NLepsRECO-NLepsRERECO)/NLepsRECO,borisRECO);

			break;//loop finished
		}//TOBTEC
		break;//loop finished
	}//RERECO
   }//RECO
   //add over- and underflow!
	cout << "good events failing tagger " << cgf << ", good events passing tagger " << cgp << ", bad events failing tagger " << cbf << ", bad events passing tagger " << cbp << endl;

 
   //save the histograms
   TFile *newfile;
   if(version==0) newfile = new TFile("/shome/haweber/TOBTECfiles/CorrelationsHistograms.root","RECREATE");
   if(version==1) newfile = new TFile("/shome/haweber/TOBTECfiles/CorrelationsHistograms_v2b.root","RECREATE");
   if(version==2) newfile = new TFile("/shome/haweber/TOBTECfiles/CorrelationsHistograms_v2c.root","RECREATE");
   newfile->cd();

   hMT2        ->Write();
   hHT         ->Write();
   hMET        ->Write();
   hVSPT       ->Write();
   hMinDPhi    ->Write();
   hJ1Pt       ->Write();
   hJ2Pt       ->Write();
   hNJets      ->Write();
   hNBJets     ->Write();
   hNLeps      ->Write();
   hMT2RE        ->Write();
   hHTRE         ->Write();
   hMETRE        ->Write();
   hVSPTRE       ->Write();
   hMinDPhiRE    ->Write();
   hJ1PtRE       ->Write();
   hJ2PtRE       ->Write();
   hNJetsRE      ->Write();
   hNBJetsRE     ->Write();
   hNLepsRE      ->Write();
   hMT2Tag     ->Write();
   hHTTag      ->Write();
   hMETTag     ->Write();
   hVSPTTag    ->Write();
   hMinDPhiTag ->Write();
   hJ1PtTag    ->Write();
   hJ2PtTag    ->Write();
   hNJetsTag   ->Write();
   hNBJetsTag  ->Write();
   hNLepsTag   ->Write();

   newfile->Close();
   cout << "histograms save in " << newfile->GetName() << endl;

}//void