#include "TEfficiency.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLine.h"
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
#include "THStack.h"
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

//run via root -l -b -q UTMcheck.C++

using namespace std;

void UTMcheck();
void printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, vector<TH1D*> h_sig, int nsig, TLegend* legend,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale);
void plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, TH1* h3, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale);
void load(const char* filename);

//struct that combines the MT2trees with important information like cross section
struct sample{
	TString name;
	TString sname;
	TString shapename;
	TString type;
	TFile *file;
	TTree *tree;
	float xsection;
	float nevents;
	float kfact;
	float PU_avg_weight;
	float lumi;
	int color;
};

vector < sample >  fSamples;
const int fVerbose = 3;
TString fPath;

//this function was used to do a check of the full UTM/VSPT distribution (incl. VSPT>70 GeV)
//in leptonic control region, e.g. for >=2b region
//needed as in original skims we did not have VSPT>70 GeV
//therefore this code should be obsolete
void UTMcheck(){

gROOT->cd();
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

//for legend
TH1D *QCDdummy = new TH1D("QCD","",1,0,1); QCDdummy->SetFillColor(401);QCDdummy->SetLineColor(401);QCDdummy->SetMarkerColor(401);
TH1D *Zdummy = new TH1D("Z","",1,0,1); Zdummy->SetFillColor(419);Zdummy->SetLineColor(419);Zdummy->SetMarkerColor(419);
TH1D *Otherdummy = new TH1D("Other","",1,0,1); Otherdummy->SetFillColor(603);Otherdummy->SetLineColor(603);Otherdummy->SetMarkerColor(603);
TH1D *Wincl = new TH1D("Wincl","", 15,0,150); Wincl->SetFillColor(417);Wincl->SetLineColor(417);Wincl->SetMarkerColor(417); Wincl->Sumw2();
TH1D *Wsignal = new TH1D("Wsignal","", 15,0,150); Wsignal->SetFillColor(417);Wsignal->SetLineColor(417);Wsignal->SetMarkerColor(417); Wsignal->Sumw2();
TH1D *TTincl = new TH1D("TTincl","", 15,0,150); TTincl->SetFillColor(855);TTincl->SetLineColor(855);TTincl->SetMarkerColor(855); TTincl->Sumw2();
TH1D *TTsignal = new TH1D("TTsignal","", 15,0,150); TTsignal->SetFillColor(855);TTsignal->SetLineColor(855);TTsignal->SetMarkerColor(855); TTsignal->Sumw2();
TH1D *Dataincl = new TH1D("Dataincl","",15,0,150); Dataincl->Sumw2(); Dataincl->SetMarkerColor(kBlack); Dataincl->SetLineColor(kBlack); Dataincl->SetMarkerStyle(20);
TH1D *Datasignal = new TH1D("Datasignal","",15,0,150); Datasignal->Sumw2(); Datasignal->SetMarkerColor(kBlack); Datasignal->SetLineColor(kBlack); Datasignal->SetMarkerStyle(20);

gROOT->cd();

TLegend* Legend1 = new TLegend(.71,.68,.91,.92);
Legend1->SetName("legend");
Legend1 -> SetFillColor(0);
Legend1 -> SetBorderSize(0);
Legend1 ->AddEntry(QCDdummy, "QCD"  , "f");
Legend1 ->AddEntry(Wincl, "W+jets"  , "f");
Legend1 ->AddEntry(Zdummy, "Z+jets"  , "f");
Legend1 ->AddEntry(TTincl, "Top"  , "f");
Legend1 ->AddEntry(Otherdummy, "Other"  , "f");
Legend1 ->AddEntry(Dataincl, "data"  , "p");

//TOBTECvariable
//needed to reload the tagger
   vector<pair<pair<int,pair<int,int> >,float> > rls; rls.clear();
   char buffer[200];
   ifstream filterdat("/shome/haweber/AODFiles_MT2SR/taggerlist/ControlRegionHTMHTMET.dat");
   while( filterdat.getline(buffer, 200, '\n') ){
	int rrun(-1), lls(-1), eevent(-1), d1(-1); float ttagger(-1);
	sscanf(buffer, "*\t%d\t*\t%d\t*\t%d\t*\t%d\t*\t%f", &d1, &eevent, &lls, &rrun, &ttagger);
	if(eevent<0||eevent>INT_MAX) cout << rrun<<":"<<lls<<":"<<eevent << " - be careful for this event number, tagger: " << ttagger << endl;
	pair<int,int> t1(lls,eevent);
	pair<int,pair<int,int> > t2(rrun,t1);
        pair<pair<int,pair<int,int> >, float> t3(t2,ttagger);
	rls.push_back(t3);
   }

   vector<pair<pair<int,int >,int> > rlsused; rlsused.clear();//to avoid doublecounting of data events

	// if you want to apply an additional event selection, define cuts here
	std::ostringstream cutStream;
	std::ostringstream cutStreamBase;
	cutStream << " " 
		<< "NJetsIDLoose40>=2"                   << "&&"
		<< "NBJets40CSVM>=2"                   << "&&"
		<< "(NTausIDLoose3Hits+NMuons+NEles)==1"                   << "&&"
		<< "((NTausIDLoose3Hits==1&&tau[0].MT<100)"                   << "||"
		<< "(NMuons==1&&muo[0].MT<100)"                              << "||"
		<< "(NEles==1&&ele[0].MT<100))"                               << "&&"
		<< "misc.Jet0Pass ==1"                      << "&&"
		<< "misc.Jet1Pass ==1"                      << "&&"
		<< "misc.MinMetJetDPhi4Pt40 >0.3";
		cutStream << "&&((misc.MET>200&&misc.HT<=750&&misc.HT>=450))";
//	         cutStream << "&&misc.MT2>=100";//lowest border in MT2
	
	cutStreamBase << " " 
      << "misc.PassJet40ID ==1"                      << "&&"
      << "misc.HBHENoiseFlag == 0"  << "&&" // for rare SM samples
      << "misc.CSCTightHaloIDFlag == 0"             << "&&"
      << "misc.trackingFailureFlag==0"              << "&&"
      << "misc.eeBadScFlag==0"                      << "&&"
      << "misc.EcalDeadCellTriggerPrimitiveFlag==0" << "&&"
      << "misc.TrackingManyStripClusFlag==0"             << "&&"
      << "misc.TrackingTooManyStripClusFlag==0"          << "&&"
      << "misc.TrackingLogErrorTooManyClustersFlag==0"   << "&&"
      << "misc.CrazyHCAL==0";
cutStreamBase << "&& misc.MET/misc.CaloMETRaw<=2.";

  std::ostringstream triggerStream;
  triggerStream << "( ( ( "
		<< "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
  		<< "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )"
		<< "||("
		<< "trigger.HLT_PFHT350_PFMET100_v3==1 || trigger.HLT_PFHT350_PFMET100_v4==1 || trigger.HLT_PFHT350_PFMET100_v5==1 || "
		<< "trigger.HLT_PFHT350_PFMET100_v6==1 || trigger.HLT_PFHT350_PFMET100_v7==1 || trigger.HLT_PFNoPUHT350_PFMET100_v1==1 || "
		<< "trigger.HLT_PFNoPUHT350_PFMET100_v3==1 || trigger.HLT_PFNoPUHT350_PFMET100_v4==1 ) ) )";

TString trigger = triggerStream.str().c_str();
  
  TString cuts = cutStream.str().c_str();
  TString basecuts = cutStreamBase.str().c_str();

TString samples = "samples/samples_dummy.dat";//this contained samples with VSPT<70 GeV as well as VSPT>=70 GeV
load(samples.Data());
  	for(size_t i = 0; i < fSamples.size(); ++i){
		//define sample
	    string sampletype = (string)fSamples[i].type;
	    if(sampletype==(string)"mc"){
		if(fSamples[i].sname=="QCD")         sampletype = (string)"QCD";
		else if(fSamples[i].sname=="Wtolnu") sampletype = (string)"WJets";
		else if(fSamples[i].sname=="DY")     sampletype = (string)"ZJets";
		else if(fSamples[i].name=="TTbar")   sampletype = (string)"TTbar";
		else if(fSamples[i].name=="TTbar_Madgraph0l")   sampletype = (string)"TTbar";
		else if(fSamples[i].name=="TTbar_Madgraph1l")   sampletype = (string)"TTbar";
		else if(fSamples[i].name=="TTbar_Madgraph2l")   sampletype = (string)"TTbar";
		else if(fSamples[i].sname=="Top")    sampletype = (string)"SingleTop";//no ttbar, includes TTZ, TTW
		else sampletype = (string)"Other";
	    }

	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "ZnunuNumbers: looping over " << fSamples[i].name << " added in " << sampletype << endl;
	    if(fVerbose>2) cout << "              sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts + "&&" + basecuts;
        
	    if( fSamples[i].type=="data") myCuts += " && " + trigger; //cuts to be aplied only on data
        
   	    cout << "Cuts for Flow: " << myCuts << endl;
   	    fSamples[i].tree->Draw(">>selList", myCuts);
        
	    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
	    fSamples[i].tree->SetEventList(myEvtList);
	    int counter=0;
	    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
	    if(myEvtList->GetSize()==0) continue;
	//run over all selected events
        while(myEvtList->GetEntry(counter++) !=-1){	
		int jentry = myEvtList->GetEntry(counter-1);
            
		nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
		fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
		
		if ( fVerbose>2 && counter % 100 == 0  )  cout << "+++ Proccessing event " << counter << endl;
		
		Double_t weight = sample_weight;
		if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 
		if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->SFWeight.BTagCSV40ge2; // pile-up reweighting for MC 
		bool beamhalo = false;//bad events = true;
		if(fMT2tree->jet[0].lv.DeltaPhi(fMT2tree->jet[1].lv)<0.2 && !(fMT2tree->jet[0].isPFIDMedium)) beamhalo = true;
		if(beamhalo==true) continue;
		if(sampletype=="WJets"){
		//separately fill WJets, TTbar, Data (rest is negligible for ==1 lepton, >=2b
		Wincl->Fill(fMT2tree->misc.Vectorsumpt, weight);
		if(fMT2tree->misc.MT2>200.) Wsignal->Fill(fMT2tree->misc.Vectorsumpt, weight);
		}	
		if(sampletype=="TTbar"){
		TTincl->Fill(fMT2tree->misc.Vectorsumpt, weight);
		if(fMT2tree->misc.MT2>200.) TTsignal->Fill(fMT2tree->misc.Vectorsumpt, weight);
		}
		if(fMT2tree->misc.isData ){
			//redo tobtec tagger
			bool keep = true;
			int p=-1;
			float TOBTECTagger = -1;
			for(unsigned int nn=0; nn<rls.size();++nn){
				if((rls[nn].first).first!=fMT2tree->misc.Run) continue; // --> run over matching 
				if(((rls[nn].first).second).first!=fMT2tree->misc.LumiSection) continue; // --> LS over matching 
				if(fMT2tree->misc.Event>0){
					if(((rls[nn].first).second).second==fMT2tree->misc.Event) { 
						TOBTECTagger = rls[nn].second;
						keep=false; p=nn; break;
					}
				}
				else{	//if misc.Event is negative due to int overflow
					int evtt = ((rls[nn].first).second).second + INT_MAX;
					evtt = evtt - INT_MAX;
					if(evtt==fMT2tree->misc.Event) {keep=false; p=nn; break; TOBTECTagger = rls[nn].second;}
				}
			}
			if(!keep){
				int tt = p;
				rls.erase(rls.begin()+p);
			}
			if(TOBTECTagger>8.) continue;
			pair<int,int> t1(fMT2tree->misc.Run,fMT2tree->misc.LumiSection);
			pair<pair<int,int>, int > t2(t1,fMT2tree->misc.Event);
			bool additional = true;
			for(unsigned int nn=0; nn<rlsused.size();++nn){
				if((rlsused[nn].first).first!=fMT2tree->misc.Run) continue; // --> run over matching 
				if( (rlsused[nn].first).second!=fMT2tree->misc.LumiSection) continue; // --> LS over matching 
				if(fMT2tree->misc.Event>0){
					if(rlsused[nn].second==fMT2tree->misc.Event) { 
						additional = false; break;
					}
				}
				else{
					int evtt = rlsused[nn].second + INT_MAX;
					evtt = evtt - INT_MAX;
					if(evtt==fMT2tree->misc.Event) {additional = false; break;}
				}
			}
			if(additional){
				rlsused.push_back(t2);
				Dataincl->Fill(fMT2tree->misc.Vectorsumpt, weight);
				if(fMT2tree->misc.MT2>200.) Datasignal->Fill(fMT2tree->misc.Vectorsumpt, weight);
			}
		}
	}//while
	delete fMT2tree;
	delete fSamples[i].tree;

	}//for samples


//add overflow bin
Dataincl->SetBinContent(15, Dataincl->GetBinContent(15)+Dataincl->GetBinContent(16));
Dataincl->SetBinError(15, sqrt(pow(Dataincl->GetBinError(15),2)+pow(Dataincl->GetBinError(16),2) ) );
Datasignal->SetBinContent(15, Datasignal->GetBinContent(15)+Datasignal->GetBinContent(16));
Datasignal->SetBinError(15, sqrt(pow(Datasignal->GetBinError(15),2)+pow(Datasignal->GetBinError(16),2) ) );
Wincl->SetBinContent(15, Wincl->GetBinContent(15)+Wincl->GetBinContent(16));
Wincl->SetBinError(15, sqrt(pow(Wincl->GetBinError(15),2)+pow(Wincl->GetBinError(16),2) ) );
Wsignal->SetBinContent(15, Wsignal->GetBinContent(15)+Wsignal->GetBinContent(16));
Wsignal->SetBinError(15, sqrt(pow(Wsignal->GetBinError(15),2)+pow(Wsignal->GetBinError(16),2) ) );
TTincl->SetBinContent(15, TTincl->GetBinContent(15)+TTincl->GetBinContent(16));
TTincl->SetBinError(15, sqrt(pow(TTincl->GetBinError(15),2)+pow(TTincl->GetBinError(16),2) ) );
TTsignal->SetBinContent(15, TTsignal->GetBinContent(15)+TTsignal->GetBinContent(16));
TTsignal->SetBinError(15, sqrt(pow(TTsignal->GetBinError(15),2)+pow(TTsignal->GetBinError(16),2) ) );

TH1D*    mcincl   = new TH1D("mcincl","", 15,0,150);mcincl-> SetFillStyle(3004);mcincl-> SetFillColor(kBlack);mcincl->Sumw2();
TH1D*    mcsignal = new TH1D("mcsignal","", 15,0,150);mcsignal-> SetFillStyle(3004);mcsignal-> SetFillColor(kBlack);mcsignal->Sumw2();

THStack* stackincl       = new THStack("stackincl", "");
THStack* stacksignal     = new THStack("stacksignal", "");
stacksignal->Add(Wsignal);
stacksignal->Add(TTsignal);
stackincl->Add(Wincl);
stackincl->Add(TTincl);
mcsignal->Add(Wsignal);
mcsignal->Add(TTsignal);
mcincl->Add(Wincl);
mcincl->Add(TTincl);

TH1D* dummySUSY = new TH1D("dummySUSY","", 15,0,150);
vector <TH1D*> hsig; hsig.push_back(dummySUSY);
TString ytitle = "Events / 10 GeV";
TString xtitle = "|MHT-MET| [GeV]";
TString outname = "VSPT_ge2jge2b_1lepCR_lowHT_noMT2cut";
printHisto(stackincl, Dataincl, mcincl, hsig, 1, Legend1 , outname, "hist", false, xtitle, ytitle, -2, -2, -10, 1.);
plotRatioStack(stackincl,  mcincl, Dataincl, mcincl,  false, false, outname, Legend1, xtitle, ytitle, -2, -2, -10, 1.);
outname  = "VSPT_ge2jge2b_1lepCR_lowHT_MT2ge200";
printHisto(stacksignal, Datasignal, mcsignal, hsig, 1, Legend1 , outname, "hist", false, xtitle, ytitle, -2, -2, -10, 1.);
plotRatioStack(stacksignal,  mcsignal, Datasignal, mcsignal,  false, false, outname, Legend1, xtitle, ytitle, -2, -2, -10, 1.);

//store all the histograms
TFile *file = new TFile("UTMcheck.root","RECREATE");
file->cd();
Dataincl->Write();
Datasignal->Write();
Wincl->Write();
Wsignal->Write();
TTincl->Write();
TTsignal->Write();
QCDdummy->Write();
Otherdummy->Write();
Zdummy->Write();
mcincl->Write();
mcsignal->Write();
stackincl->Write();
stacksignal->Write();
Legend1->Write();
file->Close();
cout << "File saved: " << file->GetName() << endl;

}

//this function is a copy from the MassPlotter.cc function - no detailed comments
void printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, vector<TH1D*> h_sig, int nsig, TLegend* legend,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 600, 600);
	col->SetRightMargin(0.08);
	col->SetLeftMargin(0.18);
	col->SetTopMargin(0.07);
	col->SetBottomMargin(0.17);

	TLegend *leg = (TLegend*) legend->Clone("leg");

	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h        -> SetMinimum(0.05);
		h_mc_sum -> SetMinimum(0.05);
		h_data   -> SetMinimum(0.05);
		for (int i=0; i<nsig; i++)  h_sig[i]-> SetMinimum(0.05);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h_mc_sum->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h_mc_sum->SetMaximum(max);
	h       ->SetMaximum(max);

	h->Draw(drawopt);
	if(h_data->Integral()>0) {
		h_data       ->Draw("sameE");
	}
	for (int i=0; i<nsig; i++) {
	  h_sig[i]->Scale(overlayScale ? overlayScale : h_data->Integral() / h_sig[i]->Integral());
	  h_sig[i]->SetLineStyle(2);
	  h_sig[i]->SetFillColor(0);
	  h_sig[i]->Draw("samehist");
	}
	if(leg != NULL ){
		leg -> SetY1NDC(0.68);
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.0275);
	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d jets",abs(njets)) : TString::Format("%d jets",abs(njets));
	text += nbjets==-10 ? "" : nbjets < 0 ? TString::Format(", #geq %d b-tag",abs(nbjets)) : TString::Format(", %d b-tag",abs(nbjets));
	text += nleps == 1 ? ", 1 lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());//standard
	TLatex LumiBox;
	LumiBox.SetNDC();
	LumiBox.SetTextSize(0.0275);
	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
	LumiBox.DrawLatex(0.43,0.943,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//standard

	h->GetXaxis()->SetTitle(xtitle);
	h->GetXaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetTitleOffset(1.1);
	
	stringstream yTitle;
		yTitle << ytitle.Data();
	h->GetYaxis()->SetTitle(yTitle.str().c_str());
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleOffset(1.47);

	gPad->RedrawAxis();
//this part has been commented for UTMcheck.C
//	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
//	if(fSave)Util::PrintEPS(col, canvname, fOutputDir);
// 	delete col;

}

//function is a copy from a MassPlotter.cc function - no detailed comments
void plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, TH1* h3, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale){

	// define canvas and pads 
	TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
	TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	h1->SetMarkerStyle(20);
	h2->SetMarkerStyle(20);

	TCanvas* c1 = new TCanvas(name+"c_ratio","",0,0,600,600);
	c1->SetFrameLineWidth(1);
	c1 -> cd();
	
	float border = 0.2;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad(name+"_plotpad",  "Pad containing the overlay plot", 0,0.211838,1,1);
	p_plot->SetLeftMargin(0.131579);
	p_plot->SetRightMargin(0.08);
	p_plot->SetTopMargin(0.06895515);
	p_plot->SetBottomMargin(0.07206074);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad(name+"_ratiopad", "Pad containing the ratio",   0,0.01863354,0.9967105,0.2189441);
	p_ratio->SetLeftMargin(0.1336634);	
	p_ratio->SetRightMargin(0.075);
	p_ratio->SetTopMargin(0.06976745);
	p_ratio->SetBottomMargin(0.2790698);

	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	if(logflag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2.5*max;
	else max = 1.5*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	hstack->SetMaximum(max);
	stringstream yTitle;
	yTitle << ytitle.Data();
	hstack->GetYaxis()->SetTitle(yTitle.str().c_str());
	hstack->GetYaxis()->SetLabelSize(0.05);
	hstack->GetYaxis()->SetTitleSize(0.05);
	hstack->GetYaxis()->SetTitleOffset(1.3);

	hstack->SetMinimum(0.02);
	hstack->Draw("hist");
	h2    ->Draw("sameE");
	h3->Scale(overlayScale ? overlayScale : h2->Integral() / h3->Integral());
	h3->SetFillColor(0);
	h3->SetLineStyle(kDotted);
	h3->SetLineWidth(4);
	h3->Draw("samehist");

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.03);

	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d jets",abs(njets)) : TString::Format("%d jets",abs(njets));
	text += nbjets==-10 ? "" : nbjets < 0 ? TString::Format(", #geq %d b-tag",abs(nbjets)) : TString::Format(", %d b-tag",abs(nbjets));
	text += nleps == 1 ? ", 1 lepton" : "";

	TitleBox.DrawLatex(0.13,0.943,text.Data());
	TLatex LumiBox;
	LumiBox.SetNDC();
	LumiBox.SetTextSize(0.0305);
	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
	LumiBox.DrawLatex(0.49,0.943,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//for CMS Preliminary

 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
	//gPad->SetLogy();
 
	TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");	
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);

 	h_ratio ->Divide(h2, h1);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	//MC with errors
	TH1D*h_ratio_mc = (TH1D*)h1_orig->Clone("h1_copy");
	h_ratio_mc->Divide(h1);
	h_ratio_mc->GetYaxis()->SetRangeUser(0,2);
	h_ratio_mc->GetXaxis()->SetLabelSize( 0.);
	h_ratio_mc->GetYaxis()->SetTitle("Data / MC");
	h_ratio_mc->GetXaxis()->SetTitle(xtitle);
	h_ratio_mc->GetXaxis()->SetTitleSize(0.2);
	h_ratio_mc->GetXaxis()->SetTitleOffset(0.5);
	h_ratio_mc->GetYaxis()->SetLabelSize(0.19);
	h_ratio_mc->GetXaxis()->SetTickLength(0.09);
	h_ratio_mc	->GetYaxis()->SetTitleSize(0.18);
	h_ratio_mc->GetYaxis()->SetTitleOffset(0.36);
	h_ratio_mc->GetYaxis()->SetNdivisions(509);

	
	h_ratio_mc->SetFillStyle(3001);
	h_ratio_mc->Draw("E2");
	h_ratio ->DrawCopy("Esame");//LEO MOD
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
//	if(fSave)Util::Print(c1, save, fOutputDir, fOutputFile);//commented for UTMcheck.C

}

//reads in the samples.dat
void load(const char* filename){
	
	fSamples.clear();
	char buffer[200];
	ifstream IN(filename);

	char ParName[100], StringValue[1000];
	float ParValue;

	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Sample File  " << filename << endl;
	int counter(0);
	
	while( IN.getline(buffer, 200, '\n') ){
		// ok = false;
		if (buffer[0] == '#') {
			continue; // Skip lines commented with '#'
		}
		if( !strcmp(buffer, "GENERAL") ){
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Path\t%s", StringValue);
			fPath = StringValue;	
			cout << fPath << endl;
			
			if(fVerbose >0){
				cout << " ----  " << endl;
				cout << "  Path " << fPath << endl;
			}

		}
		if( !strcmp(buffer, "SAMPLE")){

			sample s;
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Name\t%s", StringValue);
			s.name = TString(StringValue);

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "SName\t%s", StringValue);
			s.sname = TString(StringValue);

			IN.getline(buffer, 200, '\n');			//new
			sscanf(buffer, "ShapeName\t%s", StringValue);	//new
			s.shapename = TString(StringValue);		//new

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "File\t%s", StringValue);
			TString file =fPath+StringValue;
			TFile *f = TFile::Open(file);
			s.file = f;
			s.tree = (TTree*)f->Get("MassTree");
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Xsection\t%f", &ParValue);
			s.xsection = ParValue;
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Kfact\t%f", &ParValue);
			s.kfact = ParValue;
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Lumi\t%f", &ParValue);
			s.lumi = ParValue;

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Type\t%s", StringValue);
			s.type = StringValue;
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Color\t%f", &ParValue);
			s.color = ParValue;

			if(s.type!="data"){
			TH1F *h_PUWeights = (TH1F*) s.file->Get("h_PUWeights");
			TH1F *h_Events    = (TH1F*) s.file->Get("h_Events");
			if(h_PUWeights==0 || h_Events==0){
				cout << "ERROR: sample " << (s.file)->GetName() << " does not have PU and NEvents histos! " << endl;
				exit(1);
			}
			s.type!="data" ? s.PU_avg_weight = h_PUWeights->GetMean()    : s.PU_avg_weight =1;
			s.type!="data" ? s.nevents       = h_Events   ->GetEntries() : s.nevents       =1;
			delete h_PUWeights;
			delete h_Events;
			} else {
			s.PU_avg_weight =1;
			s.nevents       =1;
			}
			if(fVerbose > 0){
				cout << " ---- " << endl;
				cout << "  New sample added: " << s.name << endl;
				cout << "   Sample no.      " << counter << endl;
				cout << "   Short name:     " << s.sname << endl;
				cout << "   File:           " << (s.file)->GetName() << endl;
				cout << "   Events:         " << s.nevents  << endl;
				cout << "   Events in tree: " << s.tree->GetEntries() << endl; 
				cout << "   Xsection:       " << s.xsection << endl;
				cout << "   Lumi:           " << s.lumi << endl;
				cout << "   kfactor:        " << s.kfact << endl;
				cout << "   avg PU weight:  " << s.PU_avg_weight << endl;
				cout << "   type:           " << s.type << endl;
				cout << "   Color:          " << s.color << endl;
			}
			fSamples.push_back(s);
			counter++;
		}
	}
	if(fVerbose > 0) cout << "------------------------------------" << endl;
}
