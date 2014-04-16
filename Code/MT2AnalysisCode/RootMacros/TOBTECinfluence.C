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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path for MT2tree.hh

//use via root -l -b -q TOBTECinfluence.C++

using namespace std;

//this code checks the influence of the TOBTEC tagger by comparing the tagger value with 
//the original reco main quantities and the same quantities after fixed tracking (rereco)
//checks done:
//- simply plot tagger distribution before and after rereco
//- reco quantity vs. rereco quantity (before and after a TOBTEC tagger cut)
//- of tagged events (=bad) and also all events: what variables change a lot from reco --> rereco (change set relative by lvl variable)
//- plot the curve of 'good events' tagged/all vs. 'bad events' tagged/all (for various tagger values) where
// 'good' means relative change is smaller than lvl, and 'bad' means relative larger than lvl for main variables (like MT2, HT, ...)
//  - from this curve get also overall 'good' and 'bad' efficiency
void TOBTECinfluence() {

   int version = 3;//0,1,2,3 // version of Boris' TOBTEC tagger implementation
   double lvl = 0.1;//10% difference allowed in fabs(reco-rereco) > lvl * reco

   gROOT->SetStyle("Plain");
   gStyle->SetPalette(1);

   //define all histograms
   TH2D *hHTtest = new TH2D("hHTtest", "", 100,0,2000, 100,0,2000); hHTtest->GetXaxis()->SetTitle("H_{T} (rereco) [GeV]");  hHTtest->GetYaxis()->SetTitle("H_{T} (reco) [GeV]"); hHTtest->Sumw2();
   TH2D *hNJets20 = new TH2D("hNJets20", "", 10,0,10, 10,0,10); hNJets20->GetXaxis()->SetTitle("# jets (rereco)");  hNJets20->GetYaxis()->SetTitle("# jets (reco)"); hNJets20->Sumw2();
   TH2D *hMHT = new TH2D("hMHT", "", 50,0,1000, 50,0,1000); hMHT->GetXaxis()->SetTitle("H_{T}^{miss} (rereco) [GeV]");  hMHT->GetYaxis()->SetTitle("H_{T}^{miss} (reco) [GeV]"); hMHT->Sumw2();

   TH2D *hMT2 = new TH2D("hMT2", "", 50,0,1000, 50,0,1000); hMT2->GetXaxis()->SetTitle("M_{T2} (rereco) [GeV]");  hMT2->GetYaxis()->SetTitle("M_{T2} (reco) [GeV]"); hMT2->Sumw2();
   TH2D *hHT = new TH2D("hHT", "", 100,0,2000, 100,0,2000); hHT->GetXaxis()->SetTitle("H_{T} (rereco) [GeV]");  hHT->GetYaxis()->SetTitle("H_{T} (reco) [GeV]"); hHT->Sumw2();
   TH2D *hMET = new TH2D("hMET", "", 50,0,1000, 50,0,1000); hMET->GetXaxis()->SetTitle("E_{T}^{miss} (rereco) [GeV]");  hMET->GetYaxis()->SetTitle("E_{T}^{miss} (reco) [GeV]"); hMET->Sumw2();
   TH2D *hVSPT = new TH2D("hVSPT", "", 30,0,150, 30,0,150); hVSPT->GetXaxis()->SetTitle("VSPT (rereco) [GeV]");  hVSPT->GetYaxis()->SetTitle("VSPT (reco) [GeV]"); hVSPT->Sumw2();
   TH2D *hMinDPhi = new TH2D("hMinDPhi", "", 32,0,3.2, 32,0,3.2); hMinDPhi->GetXaxis()->SetTitle("min#Delta#phi (rereco)");  hMinDPhi->GetYaxis()->SetTitle("min#Delta#phi (reco)"); hMinDPhi->Sumw2();
   TH2D *hJ1Pt = new TH2D("hJ1Pt", "", 50,0,1000, 50,0,1000); hJ1Pt->GetXaxis()->SetTitle("j_{1}-p_{T} (rereco) [GeV]");  hJ1Pt->GetYaxis()->SetTitle("j_{1}-p_{T} (reco) [GeV]"); hJ1Pt->Sumw2();
   TH2D *hJ2Pt = new TH2D("hJ2Pt", "", 50,0,1000, 50,0,1000); hJ2Pt->GetXaxis()->SetTitle("j_{2}-p_{T} (rereco) [GeV]");  hJ2Pt->GetYaxis()->SetTitle("j_{2}-p_{T} (reco) [GeV]"); hJ2Pt->Sumw2();
   TH2D *hNJets = new TH2D("hNJets", "", 10,0,10, 10,0,10); hNJets->GetXaxis()->SetTitle("# jets (rereco)");  hNJets->GetYaxis()->SetTitle("# jets (reco)"); hNJets->Sumw2();
   TH2D *hNBJets = new TH2D("hNBJets", "", 5,0,5, 5,0,5); hNBJets->GetXaxis()->SetTitle("# b-jets (rereco)");  hNBJets->GetYaxis()->SetTitle("# b-jets (reco)"); hNBJets->Sumw2();
   TH2D *hNLeps = new TH2D("hNLeps", "", 5,0,5, 5,0,5); hNLeps->GetXaxis()->SetTitle("# leptons (rereco)");  hNLeps->GetYaxis()->SetTitle("# leptons (reco)"); hNLeps->Sumw2();

   TH2D *hMT2Tag = new TH2D("hMT2Tag", "", 50,0,1000, 50,0,1000); hMT2Tag->GetXaxis()->SetTitle("M_{T2} (rereco) [GeV]");  hMT2Tag->GetYaxis()->SetTitle("M_{T2} (reco) [GeV]"); hMT2Tag->Sumw2();
   TH2D *hHTTag = new TH2D("hHTTag", "", 100,0,2000, 100,0,2000); hHTTag->GetXaxis()->SetTitle("H_{T} (rereco) [GeV]");  hHTTag->GetYaxis()->SetTitle("H_{T} (reco) [GeV]"); hHTTag->Sumw2();
   TH2D *hMETTag = new TH2D("hMETTag", "", 50,0,1000, 50,0,1000); hMETTag->GetXaxis()->SetTitle("E_{T}^{miss} (rereco) [GeV]");  hMETTag->GetYaxis()->SetTitle("E_{T}^{miss} (reco) [GeV]"); hMETTag->Sumw2();
   TH2D *hVSPTTag = new TH2D("hVSPTTag", "", 30,0,150, 30,0,150); hVSPTTag->GetXaxis()->SetTitle("VSPT (rereco) [GeV]");  hVSPTTag->GetYaxis()->SetTitle("VSPT (reco) [GeV]"); hVSPTTag->Sumw2();
   TH2D *hMinDPhiTag = new TH2D("hMinDPhiTag", "", 32,0,3.2, 32,0,3.2); hMinDPhiTag->GetXaxis()->SetTitle("min#Delta#phi (rereco)");  hMinDPhiTag->GetYaxis()->SetTitle("min#Delta#phi (reco)"); hMinDPhiTag->Sumw2();
   TH2D *hJ1PtTag = new TH2D("hJ1PtTag", "", 50,0,1000, 50,0,1000); hJ1PtTag->GetXaxis()->SetTitle("j_{1}-p_{T} (rereco) [GeV]");  hJ1PtTag->GetYaxis()->SetTitle("j_{1}-p_{T} (reco) [GeV]"); hJ1PtTag->Sumw2();
   TH2D *hJ2PtTag = new TH2D("hJ2PtTag", "", 50,0,1000, 50,0,1000); hJ2PtTag->GetXaxis()->SetTitle("j_{2}-p_{T} (rereco) [GeV]");  hJ2PtTag->GetYaxis()->SetTitle("j_{2}-p_{T} (reco) [GeV]"); hJ2PtTag->Sumw2();
   TH2D *hNJetsTag = new TH2D("hNJetsTag", "", 10,0,10, 10,0,10); hNJetsTag->GetXaxis()->SetTitle("# jets (rereco)");  hNJetsTag->GetYaxis()->SetTitle("# jets (reco)"); hNJetsTag->Sumw2();
   TH2D *hNBJetsTag = new TH2D("hNBJetsTag", "", 5,0,5, 5,0,5); hNBJetsTag->GetXaxis()->SetTitle("# b-jets (rereco)");  hNBJetsTag->GetYaxis()->SetTitle("# b-jets (reco)"); hNBJetsTag->Sumw2();
   TH2D *hNLepsTag = new TH2D("hNLepsTag", "", 5,0,5, 5,0,5); hNLepsTag->GetXaxis()->SetTitle("# leptons (rereco)");  hNLepsTag->GetYaxis()->SetTitle("# leptons (reco)"); hNLepsTag->Sumw2();

   TH1D *hFailAbsCutFailure = new TH1D("hFailAbsCutFailure","",11,0,11); hFailAbsCutFailure->GetXaxis()->SetTitle("failed cut"); hFailAbsCutFailure->GetYaxis()->SetTitle("Events"); hFailAbsCutFailure->Sumw2();
   hFailAbsCutFailure->GetXaxis()->SetBinLabel(1,"no cut"); hFailAbsCutFailure->GetXaxis()->SetBinLabel(2, "NLeps"); hFailAbsCutFailure->GetXaxis()->SetBinLabel(3,"NJets"); hFailAbsCutFailure->GetXaxis()->SetBinLabel(4,"HT"); hFailAbsCutFailure->GetXaxis()->SetBinLabel(5, "MET"); hFailAbsCutFailure->GetXaxis()->SetBinLabel(6, "j_{1} p_{T}"); hFailAbsCutFailure->GetXaxis()->SetBinLabel(7, "j_{2} p_{T}"); hFailAbsCutFailure->GetXaxis()->SetBinLabel(8, "min#Delta#phi"); hFailAbsCutFailure->GetXaxis()->SetBinLabel(9, "VSPT"); hFailAbsCutFailure->GetXaxis()->SetBinLabel(10, "M_{T2}"); hFailAbsCutFailure->GetXaxis()->SetBinLabel(11,"Pass JID");

   TH1D *hPassAbsCutFailure = new TH1D("hPassAbsCutFailure","",11,0,11); hPassAbsCutFailure->GetXaxis()->SetTitle("failed cut"); hPassAbsCutFailure->GetYaxis()->SetTitle("Events"); hPassAbsCutFailure->Sumw2();
   hPassAbsCutFailure->GetXaxis()->SetBinLabel(1,"no cut"); hPassAbsCutFailure->GetXaxis()->SetBinLabel(2, "NLeps"); hPassAbsCutFailure->GetXaxis()->SetBinLabel(3,"NJets"); hPassAbsCutFailure->GetXaxis()->SetBinLabel(4,"HT"); hPassAbsCutFailure->GetXaxis()->SetBinLabel(5, "MET"); hPassAbsCutFailure->GetXaxis()->SetBinLabel(6, "j_{1} p_{T}"); hPassAbsCutFailure->GetXaxis()->SetBinLabel(7, "j_{2} p_{T}"); hPassAbsCutFailure->GetXaxis()->SetBinLabel(8, "min#Delta#phi"); hPassAbsCutFailure->GetXaxis()->SetBinLabel(9, "VSPT"); hPassAbsCutFailure->GetXaxis()->SetBinLabel(10, "M_{T2}"); hPassAbsCutFailure->GetXaxis()->SetBinLabel(11,"Pass JID");

   TH1D *hFailCutTaggedGood = new TH1D("hFailCutTaggedGood","",6,0,6); hFailCutTaggedGood->GetXaxis()->SetTitle("failed cut"); hFailCutTaggedGood->GetYaxis()->SetTitle("Events"); hFailCutTaggedGood->Sumw2();
   hFailCutTaggedGood->GetXaxis()->SetBinLabel(1,"good"); hFailCutTaggedGood->GetXaxis()->SetBinLabel(2,"fail H_{T}"); hFailCutTaggedGood->GetXaxis()->SetBinLabel(3,"fail E_{T}^{miss}"); hFailCutTaggedGood->GetXaxis()->SetBinLabel(4,"fail E_{T}^{miss}-#varphi");  hFailCutTaggedGood->GetXaxis()->SetBinLabel(5,"fail j_{1}-p_{T}"); hFailCutTaggedGood->GetXaxis()->SetBinLabel(6,"fail j_{2}-p_{T}");
   TH1D *hFailCutTaggedBad = new TH1D("hFailCutTaggedBad","",6,0,6); hFailCutTaggedBad->GetXaxis()->SetTitle("failed cut"); hFailCutTaggedBad->GetYaxis()->SetTitle("Events"); hFailCutTaggedBad->Sumw2();
   hFailCutTaggedBad->GetXaxis()->SetBinLabel(1,"good"); hFailCutTaggedBad->GetXaxis()->SetBinLabel(2,"fail H_{T}"); hFailCutTaggedBad->GetXaxis()->SetBinLabel(3,"fail E_{T}^{miss}"); hFailCutTaggedBad->GetXaxis()->SetBinLabel(4,"fail E_{T}^{miss}-#varphi");  hFailCutTaggedBad->GetXaxis()->SetBinLabel(5,"fail j_{1}-p_{T}"); hFailCutTaggedBad->GetXaxis()->SetBinLabel(6,"fail j_{2}-p_{T}");

   TH1D *hFailCutTaggedGoodAll = new TH1D("hFailCutTaggedGoodAll","",6,0,6); hFailCutTaggedGoodAll->GetXaxis()->SetTitle("failed cut"); hFailCutTaggedGoodAll->GetYaxis()->SetTitle("Events"); hFailCutTaggedGoodAll->Sumw2();
   hFailCutTaggedGoodAll->GetXaxis()->SetBinLabel(1,"good"); hFailCutTaggedGoodAll->GetXaxis()->SetBinLabel(2,"fail H_{T}"); hFailCutTaggedGoodAll->GetXaxis()->SetBinLabel(3,"fail E_{T}^{miss}"); hFailCutTaggedGoodAll->GetXaxis()->SetBinLabel(4,"fail E_{T}^{miss}-#varphi");  hFailCutTaggedGoodAll->GetXaxis()->SetBinLabel(5,"fail j_{1}-p_{T}"); hFailCutTaggedGoodAll->GetXaxis()->SetBinLabel(6,"fail j_{2}-p_{T}");
   TH1D *hFailCutTaggedBadAll = new TH1D("hFailCutTaggedBadAll","",6,0,6); hFailCutTaggedBadAll->GetXaxis()->SetTitle("failed cut"); hFailCutTaggedBadAll->GetYaxis()->SetTitle("Events"); hFailCutTaggedBadAll->Sumw2();
   hFailCutTaggedBadAll->GetXaxis()->SetBinLabel(1,"good"); hFailCutTaggedBadAll->GetXaxis()->SetBinLabel(2,"fail H_{T}"); hFailCutTaggedBadAll->GetXaxis()->SetBinLabel(3,"fail E_{T}^{miss}"); hFailCutTaggedBadAll->GetXaxis()->SetBinLabel(4,"fail E_{T}^{miss}-#varphi");  hFailCutTaggedBadAll->GetXaxis()->SetBinLabel(5,"fail j_{1}-p_{T}"); hFailCutTaggedBadAll->GetXaxis()->SetBinLabel(6,"fail j_{2}-p_{T}");

   TGraphAsymmErrors  *epass_vs_efail = new TGraphAsymmErrors(); epass_vs_efail ->SetName("epass_vs_efail"); epass_vs_efail->SetMarkerStyle(20);
   TGraphAsymmErrors  *epass_vs_efailRA2b = new TGraphAsymmErrors(); epass_vs_efailRA2b ->SetName("epass_vs_efailRA2b");epass_vs_efailRA2b->SetLineColor(kGreen-2); epass_vs_efailRA2b->SetMarkerColor(kGreen-2); epass_vs_efailRA2b->SetMarkerStyle(20);
   TGraphAsymmErrors  *egood_vs_ebad = new TGraphAsymmErrors(); egood_vs_ebad ->SetName("egood_vs_ebad"); egood_vs_ebad->SetMarkerStyle(20);
   TGraphAsymmErrors *epass_vs_tagger = new TGraphAsymmErrors(); epass_vs_tagger->SetName("epass_vs_tagger");
   TH1D *tagger_goodevents = new TH1D("tagger_goodevents","",60,0,20); tagger_goodevents->Sumw2();
   TH1D *tagger_badevents = new TH1D("tagger_badevents" ,"",60,0,20); tagger_badevents->Sumw2();
   TH1D *tagger_allevents  = new TH1D("tagger_allevents" , "",60,0,20);tagger_allevents ->Sumw2();
   TH1D *dpass = new TH1D("dpass","",2,0,2);dpass->Sumw2();
   TH1D *dtot  = new TH1D("dtot", "",2,0,2);dtot ->Sumw2();

   TH1D *taggerRECO   = new TH1D("taggerRECO"  ,"",40,0,20); taggerRECO  ->Sumw2();
   TH1D *taggerRERECO = new TH1D("taggerRERECO","",40,0,20); taggerRERECO->Sumw2();

   TFile *recofile = TFile::Open("/shome/haweber/TOBTECfiles/MT2trees_TOBTEC_RECO_MarcsFiles.root");
   TFile *rerecofile = TFile::Open("/shome/haweber/TOBTECfiles/MT2trees_TOBTEC_RERECO_MarcsFiles_noHO_fixed.root");
   TFile *tobtecfilterfile;
   if(version==0) tobtecfilterfile = TFile::Open("/shome/haweber/TOBTECfiles/tobtecFixed.root");
   if(version==1) tobtecfilterfile = TFile::Open("/shome/haweber/TOBTECfiles/tobtec_v2b.root");
   if(version==2) tobtecfilterfile = TFile::Open("/shome/haweber/TOBTECfiles/tobtec_v2c.root");
   if(version==3) tobtecfilterfile = TFile::Open("/shome/haweber/TOBTECfiles/tobtec_v2d.root");

   //input files
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
	double CaloMETRECO = fMT2treeRECO->misc.CaloMETRaw;
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
		if(eventRECO!=eventRERECO || LSRECO!=LSRERECO || runRECO!=runRERECO){//rereco tree matched to reco tree by run:LS:event
			nbytesRERECO -= nbRERECO;
			nbRERECO =  rerecotree->GetEntry(jentryRERECO);   nbytesRERECO += nbRERECO;
			rerecotree->SetBranchAddress("MT2tree", &fMT2treeRERECO);
		}
		double MT2RERECO = fMT2treeRERECO->misc.MT2;
		double HTRERECO = fMT2treeRERECO->misc.HT;
		double METRERECO = fMT2treeRERECO->misc.MET;
		double CaloMETRERECO = fMT2treeRERECO->misc.CaloMETRaw;
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
			if(evt!=eventRECO || LS!=LSRECO || run!=runRECO){//tobtec tree matched to reco tree by run:LS:event
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

			//fill histogram - tagger distribution
			taggerRECO->Fill(borisRECO);
			taggerRERECO->Fill(borisRERECO);
			//fill histogram - 2d correlation reco vs. rereco (also after tagger applied)
			hMT2->Fill(MT2RERECO,MT2RECO);
			hHT->Fill(HTRERECO,HTRECO);
			hHTtest->Fill(HTRERECO,HTRECO);
			hMET->Fill(METRERECO,METRECO);
			hMHT->Fill(MHTRERECO,MHTRECO);
			hVSPT->Fill(VSPTRERECO,VSPTRECO);
			hMinDPhi->Fill(MinDPhiRERECO,MinDPhiRECO);
			hJ1Pt->Fill(J1PtRERECO,J1PtRECO);
			hJ2Pt->Fill(J2PtRERECO,J2PtRECO);
			hNJets->Fill(NJetsRERECO,NJetsRECO);
			hNJets20->Fill(NJets20RERECO,NJets20RECO);
			hNBJets->Fill(NBJetsRERECO,NBJetsRECO);
			hNLeps->Fill(NLepsRERECO,NLepsRECO);
			if(borisRECO>=8){
				hMT2Tag->Fill(MT2RERECO,MT2RECO);
				hHTTag->Fill(HTRERECO,HTRECO);
				hMETTag->Fill(METRERECO,METRECO);
				hVSPTTag->Fill(VSPTRERECO,VSPTRECO);
				hMinDPhiTag->Fill(MinDPhiRERECO,MinDPhiRECO);
				hJ1PtTag->Fill(J1PtRERECO,J1PtRECO);
				hJ2PtTag->Fill(J2PtRERECO,J2PtRECO);
				hNJetsTag->Fill(NJetsRERECO,NJetsRECO);
				hNBJetsTag->Fill(NBJetsRERECO,NBJetsRECO);
				hNLepsTag->Fill(NLepsRERECO,NLepsRECO);
			}
			if(borisRECO>=8){
				//if tagged, what what cut does it fail (relatively)
				bool pass = true;
				if(     fabs(HTRERECO  -HTRECO  )>lvl*HTRECO  ) { pass = false; hFailCutTaggedBad->Fill(1); }
				else if(fabs(METRERECO -METRECO )>lvl*METRECO ) { pass = false; hFailCutTaggedBad->Fill(2); }
				else if(fabs(METPhiRERECO-METPhiRECO)>lvl*3.  ) { pass = false; hFailCutTaggedBad->Fill(3); }
				else if(fabs(J1PtRERECO-J1PtRECO)>lvl*J1PtRECO) { pass = false; hFailCutTaggedBad->Fill(4); }
				else if(fabs(J2PtRERECO-J2PtRECO)>lvl*J2PtRECO) { pass = false; hFailCutTaggedBad->Fill(5); }
				else hFailCutTaggedBad->Fill(0);
				if(fabs(HTRERECO  -HTRECO  )>lvl*HTRECO  ) { hFailCutTaggedBadAll->Fill(1); }
				if(fabs(METRERECO -METRECO )>lvl*METRECO ) { hFailCutTaggedBadAll->Fill(2); }
				if(fabs(METPhiRERECO-METPhiRECO)>lvl*3.  ) { hFailCutTaggedBadAll->Fill(3); }
				if(fabs(J1PtRERECO-J1PtRECO)>lvl*J1PtRECO) { hFailCutTaggedBadAll->Fill(4); }
				if(fabs(J2PtRERECO-J2PtRECO)>lvl*J2PtRECO) { hFailCutTaggedBadAll->Fill(5); }
				if(pass) hFailCutTaggedBadAll->Fill(0);
				if(     NLepsRERECO  >0  ) { hFailAbsCutFailure->Fill(1); }
				else if(NJetsRERECO  <2  ) { hFailAbsCutFailure->Fill(2); }
				else if(HTRERECO     <450) { hFailAbsCutFailure->Fill(3); }
				else if(HTRERECO<750&&METRERECO<200) { hFailAbsCutFailure->Fill(4); }
				else if(METRERECO    <30 ) { hFailAbsCutFailure->Fill(4); }
				else if(J1PtRERECO   <100) { hFailAbsCutFailure->Fill(5); }
				else if(J2PtRERECO   <100) { hFailAbsCutFailure->Fill(6); }
				else if(MinDPhiRERECO<0.3) { hFailAbsCutFailure->Fill(7); }
				else if(VSPTRERECO   > 70) { hFailAbsCutFailure->Fill(8); }
				else if(HTRERECO<750&&METRERECO>200&&MT2RERECO<200) { hFailAbsCutFailure->Fill(9); }
				else if(MT2RERECO    <100) { hFailAbsCutFailure->Fill(9); }
				else if(fMT2treeRERECO->misc.PassJet40ID==false) { hFailAbsCutFailure->Fill(10); }
				else hFailAbsCutFailure->Fill(0);
				if(pass) ++cgf;
				else     ++cbf;
			} else {
				//if not tagged, what what cut does it fail (relatively)
				bool pass = true;
				if(     fabs(HTRERECO  -HTRECO  )>lvl*HTRECO  ) { pass = false; hFailCutTaggedGood->Fill(1); }
				else if(fabs(METRERECO -METRECO )>lvl*METRECO ) { pass = false; hFailCutTaggedGood->Fill(2); }
				else if(fabs(METPhiRERECO-METPhiRECO)>lvl*3.  ) { pass = false; hFailCutTaggedGood->Fill(3); }
				else if(fabs(J1PtRERECO-J1PtRECO)>lvl*J1PtRECO) { pass = false; hFailCutTaggedGood->Fill(4); }
				else if(fabs(J2PtRERECO-J2PtRECO)>lvl*J2PtRECO) { pass = false; hFailCutTaggedGood->Fill(5); }
				else hFailCutTaggedGood->Fill(0);
				if(fabs(HTRERECO  -HTRECO  )>lvl*HTRECO  ) { hFailCutTaggedGoodAll->Fill(1); }
				if(fabs(METRERECO -METRECO )>lvl*METRECO ) { hFailCutTaggedGoodAll->Fill(2); }
				if(fabs(METPhiRERECO-METPhiRECO)>lvl*3.  ) { hFailCutTaggedGoodAll->Fill(3); }
				if(fabs(J1PtRERECO-J1PtRECO)>lvl*J1PtRECO) { hFailCutTaggedGoodAll->Fill(4); }
				if(fabs(J2PtRERECO-J2PtRECO)>lvl*J2PtRECO) { hFailCutTaggedGoodAll->Fill(5); }
				if(pass) hFailCutTaggedGoodAll->Fill(0);
				if(     NLepsRERECO  >0  ) { hPassAbsCutFailure->Fill(1); }
				else if(NJetsRERECO  <2  ) { hPassAbsCutFailure->Fill(2); }
				else if(HTRERECO     <450) { hPassAbsCutFailure->Fill(3); }
				else if(HTRERECO<750&&METRERECO<200) { hPassAbsCutFailure->Fill(4); }
				else if(METRERECO    <30 ) { hPassAbsCutFailure->Fill(4); }
				else if(J1PtRERECO   <100) { hPassAbsCutFailure->Fill(5); }
				else if(J2PtRERECO   <100) { hPassAbsCutFailure->Fill(6); }
				else if(MinDPhiRERECO<0.3) { hPassAbsCutFailure->Fill(7); }
				else if(VSPTRERECO   > 70) { hPassAbsCutFailure->Fill(8); }
				else if(HTRERECO<750&&METRERECO>200&&MT2RERECO<200) { hPassAbsCutFailure->Fill(9); }
				else if(MT2RERECO    <100) { hPassAbsCutFailure->Fill(9); }
				else if(fMT2treeRERECO->misc.PassJet40ID==false) { hPassAbsCutFailure->Fill(10); }
				else  hPassAbsCutFailure->Fill(0); 
				if(pass) ++cgp;
				else     ++cbp;
			}
			bool pass = true;
			if(     fabs(HTRERECO  -HTRECO  )>lvl*HTRECO  ) { pass = false; }
			else if(fabs(METRERECO -METRECO )>lvl*METRECO ) { pass = false; }
			else if(fabs(METPhiRERECO-METPhiRECO)>lvl*3.  ) { pass = false; }
			else if(fabs(J1PtRERECO-J1PtRECO)>lvl*J1PtRECO) { pass = false; }
			else if(fabs(J2PtRERECO-J2PtRECO)>lvl*J2PtRECO) { pass = false; }

			//tagger distribution of good and bad events
			if(pass) tagger_goodevents->Fill(borisRECO);
			else     tagger_badevents->Fill(borisRECO);
			if(!pass) cout << run << ":" << LS << ":" << evt << endl;
			tagger_allevents->Fill(borisRECO);
			int RA2b = -1;
			if(RA2bRECO) RA2b=0; else RA2b = 1;
			if(pass) dpass->Fill(RA2b);
			dtot->Fill(RA2b);

			break;//loop finished
		}//TOBTEC
		break;//loop finished
	}//RERECO
   }//RECO

   //add over- and underflow!
   tagger_goodevents->SetBinContent(tagger_goodevents->GetNbinsX(), tagger_goodevents->GetBinContent(tagger_goodevents->GetNbinsX())+tagger_goodevents->GetBinContent(tagger_goodevents->GetNbinsX()+1));
   tagger_allevents->SetBinContent(tagger_allevents->GetNbinsX(), tagger_allevents->GetBinContent(tagger_allevents->GetNbinsX())+tagger_allevents->GetBinContent(tagger_allevents->GetNbinsX()+1));
   for(int x = 1; x<=100; ++x){
//   for(int y = 1; y<=100; ++y){
	if(hVSPT->GetBinContent(x,hVSPT->GetNbinsX()+1)>0) hVSPT->SetBinContent(x,hVSPT->GetNbinsX(), hVSPT->GetBinContent(x,hVSPT->GetNbinsX())+hVSPT->GetBinContent(x,hVSPT->GetNbinsX()+1));
	if(hVSPT->GetBinContent(hVSPT->GetNbinsX()+1,x)>0) hVSPT->SetBinContent(hVSPT->GetNbinsX(),x, hVSPT->GetBinContent(hVSPT->GetNbinsX(),x)+hVSPT->GetBinContent(hVSPT->GetNbinsX()+1,x));
	if(hVSPT->GetBinContent(x,0)>0) hVSPT->SetBinContent(x,1, hVSPT->GetBinContent(x,1)+hVSPT->GetBinContent(x,0));
	if(hVSPT->GetBinContent(0,x)>0) hVSPT->SetBinContent(1,x, hVSPT->GetBinContent(1,x)+hVSPT->GetBinContent(0,x));
	if(hVSPTTag->GetBinContent(x,hVSPT->GetNbinsX()+1)>0) hVSPTTag->SetBinContent(x,hVSPTTag->GetNbinsX(), hVSPTTag->GetBinContent(x,hVSPT->GetNbinsX())+hVSPTTag->GetBinContent(x,hVSPT->GetNbinsX()+1));
	if(hVSPTTag->GetBinContent(hVSPT->GetNbinsX()+1,x)>0) hVSPTTag->SetBinContent(hVSPTTag->GetNbinsX(),x, hVSPTTag->GetBinContent(hVSPT->GetNbinsX(),x)+hVSPTTag->GetBinContent(hVSPT->GetNbinsX()+1,x));
	if(hVSPTTag->GetBinContent(x,0)>0) hVSPTTag->SetBinContent(x,1, hVSPTTag->GetBinContent(x,1)+hVSPTTag->GetBinContent(x,0));
	if(hVSPTTag->GetBinContent(0,x)>0) hVSPTTag->SetBinContent(1,x, hVSPTTag->GetBinContent(1,x)+hVSPTTag->GetBinContent(0,x));
	if(x>90||x>90) continue;
	if(hHT->GetBinContent(x,hHT->GetNbinsX()+1)>0) hHT->SetBinContent(x,hHT->GetNbinsX(), hHT->GetBinContent(x,hHT->GetNbinsX())+hHT->GetBinContent(x,hHT->GetNbinsX()+1));
	if(hHT->GetBinContent(hHT->GetNbinsX()+1,x)>0) hHT->SetBinContent(hHT->GetNbinsX(),x, hHT->GetBinContent(hHT->GetNbinsX(),x)+hHT->GetBinContent(hHT->GetNbinsX()+1,x));
	if(hHT->GetBinContent(x,0)>0) hHT->SetBinContent(x,1, hHT->GetBinContent(x,1)+hHT->GetBinContent(x,0));
	if(hHT->GetBinContent(0,x)>0) hHT->SetBinContent(1,x, hHT->GetBinContent(1,x)+hHT->GetBinContent(0,x));
	if(hHTTag->GetBinContent(x,hHT->GetNbinsX()+1)>0) hHTTag->SetBinContent(x,hHTTag->GetNbinsX(), hHTTag->GetBinContent(x,hHT->GetNbinsX())+hHTTag->GetBinContent(x,hHT->GetNbinsX()+1));
	if(hHTTag->GetBinContent(hHT->GetNbinsX()+1,x)>0) hHTTag->SetBinContent(hHTTag->GetNbinsX(),x, hHTTag->GetBinContent(hHT->GetNbinsX(),x)+hHTTag->GetBinContent(hHT->GetNbinsX()+1,x));
	if(hHTTag->GetBinContent(x,0)>0) hHTTag->SetBinContent(x,1, hHTTag->GetBinContent(x,1)+hHTTag->GetBinContent(x,0));
	if(hHTTag->GetBinContent(0,x)>0) hHTTag->SetBinContent(1,x, hHTTag->GetBinContent(1,x)+hHTTag->GetBinContent(0,x));
	if(x>50||x>50) continue;
	if(hMT2->GetBinContent(x,hMT2->GetNbinsX()+1)>0) hMT2->SetBinContent(x,hMT2->GetNbinsX(), hMT2->GetBinContent(x,hMT2->GetNbinsX())+hMT2->GetBinContent(x,hMT2->GetNbinsX()+1));
	if(hMT2->GetBinContent(hMT2->GetNbinsX()+1,x)>0) hMT2->SetBinContent(hMT2->GetNbinsX(),x, hMT2->GetBinContent(hMT2->GetNbinsX(),x)+hMT2->GetBinContent(hMT2->GetNbinsX()+1,x));
	if(hMT2->GetBinContent(x,0)>0) hMT2->SetBinContent(x,1, hMT2->GetBinContent(x,1)+hMT2->GetBinContent(x,0));
	if(hMT2->GetBinContent(0,x)>0) hMT2->SetBinContent(1,x, hMT2->GetBinContent(1,x)+hMT2->GetBinContent(0,x));
	if(hMT2Tag->GetBinContent(x,hMT2Tag->GetNbinsX()+1)>0) hMT2Tag->SetBinContent(x,hMT2Tag->GetNbinsX(), hMT2Tag->GetBinContent(x,hMT2Tag->GetNbinsX())+hMT2Tag->GetBinContent(x,hMT2Tag->GetNbinsX()+1));
	if(hMT2Tag->GetBinContent(hMT2Tag->GetNbinsX()+1,x)>0) hMT2Tag->SetBinContent(hMT2Tag->GetNbinsX(),x, hMT2Tag->GetBinContent(hMT2Tag->GetNbinsX(),x)+hMT2Tag->GetBinContent(hMT2Tag->GetNbinsX()+1,x));
	if(hMT2Tag->GetBinContent(x,0)>0) hMT2Tag->SetBinContent(x,1, hMT2Tag->GetBinContent(x,1)+hMT2Tag->GetBinContent(x,0));
	if(hMT2Tag->GetBinContent(0,x)>0) hMT2Tag->SetBinContent(1,x, hMT2Tag->GetBinContent(1,x)+hMT2Tag->GetBinContent(0,x));
	if(hMET->GetBinContent(x,hMET->GetNbinsX()+1)>0) hMET->SetBinContent(x,hMET->GetNbinsX(), hMET->GetBinContent(x,hMET->GetNbinsX())+hMET->GetBinContent(x,hMET->GetNbinsX()+1));
	if(hMET->GetBinContent(hMET->GetNbinsX()+1,x)>0) hMET->SetBinContent(hMET->GetNbinsX(),x, hMET->GetBinContent(hMET->GetNbinsX(),x)+hMET->GetBinContent(hMET->GetNbinsX()+1,x));
	if(hMET->GetBinContent(x,0)>0) hMET->SetBinContent(x,1, hMET->GetBinContent(x,1)+hMET->GetBinContent(x,0));
	if(hMET->GetBinContent(0,x)>0) hMET->SetBinContent(1,x, hMET->GetBinContent(1,x)+hMET->GetBinContent(0,x));
	if(hMETTag->GetBinContent(x,hMET->GetNbinsX()+1)>0) hMETTag->SetBinContent(x,hMETTag->GetNbinsX(), hMETTag->GetBinContent(x,hMET->GetNbinsX())+hMETTag->GetBinContent(x,hMET->GetNbinsX()+1));
	if(hMETTag->GetBinContent(hMET->GetNbinsX()+1,x)>0) hMETTag->SetBinContent(hMETTag->GetNbinsX(),x, hMETTag->GetBinContent(hMET->GetNbinsX(),x)+hMETTag->GetBinContent(hMET->GetNbinsX()+1,x));
	if(hMETTag->GetBinContent(x,0)>0) hMETTag->SetBinContent(x,1, hMETTag->GetBinContent(x,1)+hMETTag->GetBinContent(x,0));
	if(hMETTag->GetBinContent(0,x)>0) hMETTag->SetBinContent(1,x, hMETTag->GetBinContent(1,x)+hMETTag->GetBinContent(0,x));
	if(hJ1Pt->GetBinContent(x,hJ1Pt->GetNbinsX()+1)>0) hJ1Pt->SetBinContent(x,hJ1Pt->GetNbinsX(), hJ1Pt->GetBinContent(x,hJ1Pt->GetNbinsX())+hJ1Pt->GetBinContent(x,hJ1Pt->GetNbinsX()+1));
	if(hJ1Pt->GetBinContent(hJ1Pt->GetNbinsX()+1,x)>0) hJ1Pt->SetBinContent(hJ1Pt->GetNbinsX(),x, hJ1Pt->GetBinContent(hJ1Pt->GetNbinsX(),x)+hJ1Pt->GetBinContent(hJ1Pt->GetNbinsX()+1,x));
	if(hJ1Pt->GetBinContent(x,0)>0) hJ1Pt->SetBinContent(x,1, hJ1Pt->GetBinContent(x,1)+hJ1Pt->GetBinContent(x,0));
	if(hJ1Pt->GetBinContent(0,x)>0) hJ1Pt->SetBinContent(1,x, hJ1Pt->GetBinContent(1,x)+hJ1Pt->GetBinContent(0,x));
	if(hJ1PtTag->GetBinContent(x,hJ1Pt->GetNbinsX()+1)>0) hJ1PtTag->SetBinContent(x,hJ1PtTag->GetNbinsX(), hJ1PtTag->GetBinContent(x,hJ1Pt->GetNbinsX())+hJ1PtTag->GetBinContent(x,hJ1Pt->GetNbinsX()+1));
	if(hJ1PtTag->GetBinContent(hJ1Pt->GetNbinsX()+1,x)>0) hJ1PtTag->SetBinContent(hJ1PtTag->GetNbinsX(),x, hJ1PtTag->GetBinContent(hJ1Pt->GetNbinsX(),x)+hJ1PtTag->GetBinContent(hJ1Pt->GetNbinsX()+1,x));
	if(hJ1PtTag->GetBinContent(x,0)>0) hJ1PtTag->SetBinContent(x,1, hJ1PtTag->GetBinContent(x,1)+hJ1PtTag->GetBinContent(x,0));
	if(hJ1PtTag->GetBinContent(0,x)>0) hJ1PtTag->SetBinContent(1,x, hJ1PtTag->GetBinContent(1,x)+hJ1PtTag->GetBinContent(0,x));
	if(hJ2Pt->GetBinContent(x,hJ2Pt->GetNbinsX()+1)>0) hJ2Pt->SetBinContent(x,hJ2Pt->GetNbinsX(), hJ2Pt->GetBinContent(x,hJ2Pt->GetNbinsX())+hJ2Pt->GetBinContent(x,hJ2Pt->GetNbinsX()+1));
	if(hJ2Pt->GetBinContent(hJ2Pt->GetNbinsX()+1,x)>0) hJ2Pt->SetBinContent(hJ2Pt->GetNbinsX(),x, hJ2Pt->GetBinContent(hJ2Pt->GetNbinsX(),x)+hJ2Pt->GetBinContent(hJ2Pt->GetNbinsX()+1,x));
	if(hJ2Pt->GetBinContent(x,0)>0) hJ2Pt->SetBinContent(x,1, hJ2Pt->GetBinContent(x,1)+hJ2Pt->GetBinContent(x,0));
	if(hJ2Pt->GetBinContent(0,x)>0) hJ2Pt->SetBinContent(1,x, hJ2Pt->GetBinContent(1,x)+hJ2Pt->GetBinContent(0,x));
	if(hJ2PtTag->GetBinContent(x,hJ2Pt->GetNbinsX()+1)>0) hJ2PtTag->SetBinContent(x,hJ2PtTag->GetNbinsX(), hJ2PtTag->GetBinContent(x,hJ2Pt->GetNbinsX())+hJ2PtTag->GetBinContent(x,hJ2Pt->GetNbinsX()+1));
	if(hJ2PtTag->GetBinContent(hJ2Pt->GetNbinsX()+1,x)>0) hJ2PtTag->SetBinContent(hJ2PtTag->GetNbinsX(),x, hJ2PtTag->GetBinContent(hJ2Pt->GetNbinsX(),x)+hJ2PtTag->GetBinContent(hJ2Pt->GetNbinsX()+1,x));
	if(hJ2PtTag->GetBinContent(x,0)>0) hJ2PtTag->SetBinContent(x,1, hJ2PtTag->GetBinContent(x,1)+hJ2PtTag->GetBinContent(x,0));
	if(hJ2PtTag->GetBinContent(0,x)>0) hJ2PtTag->SetBinContent(1,x, hJ2PtTag->GetBinContent(1,x)+hJ2PtTag->GetBinContent(0,x));
	if(x>32||x>32) continue;
	if(hMinDPhi->GetBinContent(x,hMinDPhi->GetNbinsX()+1)>0) hMinDPhi->SetBinContent(x,hMinDPhi->GetNbinsX(), hMinDPhi->GetBinContent(x,hMinDPhi->GetNbinsX())+hMinDPhi->GetBinContent(x,hMinDPhi->GetNbinsX()+1));
	if(hMinDPhi->GetBinContent(hMinDPhi->GetNbinsX()+1,x)>0) hMinDPhi->SetBinContent(hMinDPhi->GetNbinsX(),x, hMinDPhi->GetBinContent(hMinDPhi->GetNbinsX(),x)+hMinDPhi->GetBinContent(hMinDPhi->GetNbinsX()+1,x));
	if(hMinDPhi->GetBinContent(x,0)>0) hMinDPhi->SetBinContent(x,1, hMinDPhi->GetBinContent(x,1)+hMinDPhi->GetBinContent(x,0));
	if(hMinDPhi->GetBinContent(0,x)>0) hMinDPhi->SetBinContent(1,x, hMinDPhi->GetBinContent(1,x)+hMinDPhi->GetBinContent(0,x));
	if(hMinDPhiTag->GetBinContent(x,hMinDPhi->GetNbinsX()+1)>0) hMinDPhiTag->SetBinContent(x,hMinDPhiTag->GetNbinsX(), hMinDPhiTag->GetBinContent(x,hMinDPhi->GetNbinsX())+hMinDPhiTag->GetBinContent(x,hMinDPhi->GetNbinsX()+1));
	if(hMinDPhiTag->GetBinContent(hMinDPhi->GetNbinsX()+1,x)>0) hMinDPhiTag->SetBinContent(hMinDPhiTag->GetNbinsX(),x, hMinDPhiTag->GetBinContent(hMinDPhi->GetNbinsX(),x)+hMinDPhiTag->GetBinContent(hMinDPhi->GetNbinsX()+1,x));
	if(hMinDPhiTag->GetBinContent(x,0)>0) hMinDPhiTag->SetBinContent(x,1, hMinDPhiTag->GetBinContent(x,1)+hMinDPhiTag->GetBinContent(x,0));
	if(hMinDPhiTag->GetBinContent(0,x)>0) hMinDPhiTag->SetBinContent(1,x, hMinDPhiTag->GetBinContent(1,x)+hMinDPhiTag->GetBinContent(0,x));
	if(x>5||x>5) continue;
	if(hNJets->GetBinContent(x,hNJets->GetNbinsX()+1)>0) hNJets->SetBinContent(x,hNJets->GetNbinsX(), hNJets->GetBinContent(x,hNJets->GetNbinsX())+hNJets->GetBinContent(x,hNJets->GetNbinsX()+1));
	if(hNJets->GetBinContent(hNJets->GetNbinsX()+1,x)>0) hNJets->SetBinContent(hNJets->GetNbinsX(),x, hNJets->GetBinContent(hNJets->GetNbinsX(),x)+hNJets->GetBinContent(hNJets->GetNbinsX()+1,x));
	if(hNJets->GetBinContent(x,0)>0) hNJets->SetBinContent(x,1, hNJets->GetBinContent(x,1)+hNJets->GetBinContent(x,0));
	if(hNJets->GetBinContent(0,x)>0) hNJets->SetBinContent(1,x, hNJets->GetBinContent(1,x)+hNJets->GetBinContent(0,x));
	if(hNJetsTag->GetBinContent(x,hNJets->GetNbinsX()+1)>0) hNJetsTag->SetBinContent(x,hNJetsTag->GetNbinsX(), hNJetsTag->GetBinContent(x,hNJets->GetNbinsX())+hNJetsTag->GetBinContent(x,hNJets->GetNbinsX()+1));
	if(hNJetsTag->GetBinContent(hNJets->GetNbinsX()+1,x)>0) hNJetsTag->SetBinContent(hNJetsTag->GetNbinsX(),x, hNJetsTag->GetBinContent(hNJets->GetNbinsX(),x)+hNJetsTag->GetBinContent(hNJets->GetNbinsX()+1,x));
	if(hNJetsTag->GetBinContent(x,0)>0) hNJetsTag->SetBinContent(x,1, hNJetsTag->GetBinContent(x,1)+hNJetsTag->GetBinContent(x,0));
	if(hNJetsTag->GetBinContent(0,x)>0) hNJetsTag->SetBinContent(1,x, hNJetsTag->GetBinContent(1,x)+hNJetsTag->GetBinContent(0,x));
	if(+hNBJets->GetBinContent(x,hNBJets->GetNbinsX()+1)>0) hNBJets->SetBinContent(x,hNBJets->GetNbinsX(), hNBJets->GetBinContent(x,hNBJets->GetNbinsX())+hNBJets->GetBinContent(x,hNBJets->GetNbinsX()+1));
	if(hNBJets->GetBinContent(hNBJets->GetNbinsX()+1,x)>0) hNBJets->SetBinContent(hNBJets->GetNbinsX(),x, hNBJets->GetBinContent(hNBJets->GetNbinsX(),x)+hNBJets->GetBinContent(hNBJets->GetNbinsX()+1,x));
	if(hNBJets->GetBinContent(x,0)>0) hNBJets->SetBinContent(x,1, hNBJets->GetBinContent(x,1)+hNBJets->GetBinContent(x,0));
	if(hNBJets->GetBinContent(0,x)>0) hNBJets->SetBinContent(1,x, hNBJets->GetBinContent(1,x)+hNBJets->GetBinContent(0,x));
	if(hNBJetsTag->GetBinContent(x,hNBJets->GetNbinsX()+1)>0) hNBJetsTag->SetBinContent(x,hNBJetsTag->GetNbinsX(), hNBJetsTag->GetBinContent(x,hNBJets->GetNbinsX())+hNBJetsTag->GetBinContent(x,hNBJets->GetNbinsX()+1));
	if(hNBJetsTag->GetBinContent(hNBJets->GetNbinsX()+1,x)>0) hNBJetsTag->SetBinContent(hNBJetsTag->GetNbinsX(),x, hNBJetsTag->GetBinContent(hNBJets->GetNbinsX(),x)+hNBJetsTag->GetBinContent(hNBJets->GetNbinsX()+1,x));
	if(hNBJetsTag->GetBinContent(x,0)>0) hNBJetsTag->SetBinContent(x,1, hNBJetsTag->GetBinContent(x,1)+hNBJetsTag->GetBinContent(x,0));
	if(hNBJetsTag->GetBinContent(0,x)>0) hNBJetsTag->SetBinContent(1,x, hNBJetsTag->GetBinContent(1,x)+hNBJetsTag->GetBinContent(0,x));
	if(hNLeps->GetBinContent(x,hNLeps->GetNbinsX()+1)>0) hNLeps->SetBinContent(x,hNLeps->GetNbinsX(), hNLeps->GetBinContent(x,hNLeps->GetNbinsX())+hNLeps->GetBinContent(x,hNLeps->GetNbinsX()+1));
	if(hNLeps->GetBinContent(hNLeps->GetNbinsX()+1,x)>0) hNLeps->SetBinContent(hNLeps->GetNbinsX(),x, hNLeps->GetBinContent(hNLeps->GetNbinsX(),x)+hNLeps->GetBinContent(hNLeps->GetNbinsX()+1,x));
	if(hNLeps->GetBinContent(x,0)>0) hNLeps->SetBinContent(x,1, hNLeps->GetBinContent(x,1)+hNLeps->GetBinContent(x,0));
	if(hNLeps->GetBinContent(0,x)>0) hNLeps->SetBinContent(1,x, hNLeps->GetBinContent(1,x)+hNLeps->GetBinContent(0,x));
	if(hNLepsTag->GetBinContent(x,hNLeps->GetNbinsX()+1)>0) hNLepsTag->SetBinContent(x,hNLepsTag->GetNbinsX(), hNLepsTag->GetBinContent(x,hNLeps->GetNbinsX())+hNLepsTag->GetBinContent(x,hNLeps->GetNbinsX()+1));
	if(hNLepsTag->GetBinContent(hNLeps->GetNbinsX()+1,x)>0) hNLepsTag->SetBinContent(hNLepsTag->GetNbinsX(),x, hNLepsTag->GetBinContent(hNLeps->GetNbinsX(),x)+hNLepsTag->GetBinContent(hNLeps->GetNbinsX()+1,x));
	if(hNLepsTag->GetBinContent(x,0)>0) hNLepsTag->SetBinContent(x,1, hNLepsTag->GetBinContent(x,1)+hNLepsTag->GetBinContent(x,0));
	if(hNLepsTag->GetBinContent(0,x)>0) hNLepsTag->SetBinContent(1,x, hNLepsTag->GetBinContent(1,x)+hNLepsTag->GetBinContent(0,x));
   }//}
	cout << "good events failing tagger " << cgf << ", good events passing tagger " << cgp << ", bad events failing tagger " << cbf << ", bad events passing tagger " << cbp << endl;

   //get passing/failing efficiency vs. tagger value
   for(int i = 1; i<=tagger_goodevents->GetNbinsX(); ++i){
	TH1D *dp = new TH1D("dp","",1,0,1);
	TH1D *dt = new TH1D("dt","",1,0,1);
	dp->SetBinContent(1,tagger_goodevents->Integral(1,i));
	dt->SetBinContent(1,tagger_allevents ->Integral(1,i));
	TEfficiency *teff = new TEfficiency((*dp),(*dt));
	double eff = teff->GetEfficiency(1);//efficiency of tagger to tag good events
	double efferrl = teff->GetEfficiencyErrorLow(1);
	double efferru = teff->GetEfficiencyErrorUp(1);
	int n = epass_vs_tagger->GetN();
	epass_vs_tagger->SetPoint(n,tagger_goodevents->GetBinLowEdge(i),eff);
	epass_vs_tagger->SetPointEYhigh(n, efferru);
	epass_vs_tagger->SetPointEYlow(n, efferrl);

	TH1D *dpe = new TH1D("dpe","",1,0,1);
	TH1D *dte = new TH1D("dte","",1,0,1);
	TH1D *dpf = new TH1D("dpf","",1,0,1);
	TH1D *dtf = new TH1D("dtf","",1,0,1);
	dpe->SetBinContent(1,tagger_goodevents->Integral(1,i));
	dte->SetBinContent(1,tagger_goodevents->Integral()   );
	TEfficiency *teffe = new TEfficiency((*dpe),(*dte));
	double effe = teffe->GetEfficiency(1);//efficiency of selection events of passing the tagger
	double efferrle = teffe->GetEfficiencyErrorLow(1);
	double efferrue = teffe->GetEfficiencyErrorUp(1);
	dpf->SetBinContent(1,tagger_badevents->Integral(1,i));
	dtf->SetBinContent(1,tagger_badevents->Integral()   );
	TEfficiency *tefff = new TEfficiency((*dpf),(*dtf));//effieciency of selection failing events to pass the tagger
	double efff = tefff->GetEfficiency(1);
	double efferrlf = tefff->GetEfficiencyErrorLow(1);
	double efferruf = tefff->GetEfficiencyErrorUp(1);
	n = epass_vs_efail->GetN();
	epass_vs_efail->SetPoint(n, effe, efff);
	epass_vs_efail->SetPointError(n,efferrle,efferrue,efferrlf,efferruf);
	double eee = efferrle; if(efferrle<efferrue) eee = efferrue;
	double eef = efferrlf; if(efferrlf<efferruf) eef = efferruf;
	//cout << "point " << n << " with tagger cut " << tagger_goodevents->GetBinLowEdge(i)+tagger_goodevents->GetBinWidth(i) << " has tagger efficiency for passing/failing events " << effe << "+/-" << eee << " / " << efff << "+/-" << eef << endl;

	delete dpe;delete dte;delete dpf;delete dtf;delete dt;delete dp;
	delete teff;delete teffe;delete tefff;

        //get 'good' efficiency vs. 'bad' efficiency
	TH1D *dpg = new TH1D("dpg","",1,0,1);
	TH1D *dtg = new TH1D("dtg","",1,0,1);
	TH1D *dpb = new TH1D("dpb","",1,0,1);
	TH1D *dtb = new TH1D("dtb","",1,0,1);
	dpg->SetBinContent(1,tagger_goodevents->Integral(1,i));
	dtg->SetBinContent(1,tagger_allevents ->Integral(1,i));
	TEfficiency *teffg = new TEfficiency((*dpg),(*dtg));
	double effg = teffg->GetEfficiency(1);//efficiency of events passing the tagger to be selected
	double efferrlg = teffg->GetEfficiencyErrorLow(1);
	double efferrug = teffg->GetEfficiencyErrorUp(1);
	dpb->SetBinContent(1,tagger_goodevents->Integral(i+1,tagger_goodevents->GetNbinsX()));
	dtb->SetBinContent(1,tagger_allevents ->Integral(i+1,tagger_allevents ->GetNbinsX()));
	TEfficiency *teffb = new TEfficiency((*dpb),(*dtb));//effieciency of events failing the tagger to be selected
	double effb = teffb->GetEfficiency(1);
	double efferrlb = teffb->GetEfficiencyErrorLow(1);
	double efferrub = teffb->GetEfficiencyErrorUp(1);
	n = egood_vs_ebad->GetN();
	egood_vs_ebad->SetPoint(n, effg, effb);
	egood_vs_ebad->SetPointError(n,efferrlg,efferrug,efferrlb,efferrub);
	delete dpg;delete dtg;delete dpb;delete dtb;
	delete teffg;delete teffb;

   }
	//cross check of pass/fail efficiency for RA2b tagger
	TH1D *dpe = new TH1D("dpe","",1,0,1);
	TH1D *dte = new TH1D("dte","",1,0,1);
	TH1D *dpf = new TH1D("dpf","",1,0,1);
	TH1D *dtf = new TH1D("dtf","",1,0,1);
	dpe->SetBinContent(1,dpass->GetBinContent(2));
	dte->SetBinContent(1,dpass->Integral()   );
	TEfficiency *teffe = new TEfficiency((*dpe),(*dte));
	double effe = teffe->GetEfficiency(1);//efficiency of good events of passing the tagger
	double efferrle = teffe->GetEfficiencyErrorLow(1);
	double efferrue = teffe->GetEfficiencyErrorUp(1);
	dpf->SetBinContent(1,dtot->GetBinContent(2)-dpass->GetBinContent(2));
	dtf->SetBinContent(1,dtot->Integral()      -dpass->Integral()   );
	TEfficiency *tefff = new TEfficiency((*dpf),(*dtf));//effieciency of failing events to pass the tagger
	double efff = tefff->GetEfficiency(1);
	double efferrlf = tefff->GetEfficiencyErrorLow(1);
	double efferruf = tefff->GetEfficiencyErrorUp(1);
	int n = epass_vs_efailRA2b->GetN();
	epass_vs_efailRA2b->SetPoint(n, effe, efff);
	epass_vs_efailRA2b->SetPointError(n,efferrle,efferrue,efferrlf,efferruf);

	delete dpe;delete dte;delete dpf;delete dtf;
	delete teffe;delete tefff;

   //save all histograms, efficiencies, and curves...
   TFile *newfile;
   if(version==0) newfile = new TFile("/shome/haweber/TOBTECfiles/ChecksnoHOHistogramsFixed.root","RECREATE");
   if(version==1) newfile = new TFile("/shome/haweber/TOBTECfiles/ChecksnoHOHistograms_v2b.root","RECREATE");
   if(version==2) newfile = new TFile("/shome/haweber/TOBTECfiles/ChecksnoHOHistograms_v2c.root","RECREATE");
   if(version==3) newfile = new TFile("/shome/haweber/TOBTECfiles/ChecksnoHOHistograms_v2d.root","RECREATE");
   newfile->cd();

   taggerRECO  ->Write();
   taggerRERECO->Write();

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

   hFailAbsCutFailure   ->Write();
   hPassAbsCutFailure   ->Write();
   hFailCutTaggedGood   ->Write();
   hFailCutTaggedBad    ->Write();
   hFailCutTaggedGoodAll->Write();
   hFailCutTaggedBadAll ->Write();

   tagger_badevents ->Write();
   tagger_goodevents ->Write();
   tagger_allevents  ->Write();

   hHTtest     ->Write();
   hNJets20    ->Write();
   hMHT        ->Write();

   epass_vs_tagger   ->Write();
   epass_vs_efail    ->Write();
   epass_vs_efailRA2b->Write();
   egood_vs_ebad     ->Write();

   newfile->Close();
   cout << "histograms save in " << newfile->GetName() << endl;

}//void