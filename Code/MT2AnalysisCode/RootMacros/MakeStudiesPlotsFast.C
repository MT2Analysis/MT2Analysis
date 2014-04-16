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
#include "TLegend.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLine.h"
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

//run via root -l -b -q MakeStudiesPlotsFast.C++

using namespace std;

void MakeStudiesPlotsFast();
void MakePlots(map<string, TH1D*> histos, const int signalregionsize, string *signalregions, const int controlsize, string *control, const int mindphisize, string *minmetjetdphi, map<string, THStack*> stacks, TLegend *Legend1, vector<string> histonames, TString outputdirs);
void Make1DPlotsRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, TH1D* histsusy, bool logFlag, bool normalize, TString outname, TLegend *leg, TString xtitle, TString ytitle, float overlayScale);
void Make1DPlotsNoRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, TH1D* histsusy, bool logFlag, bool normalize, TString outname, TLegend *leg, TString xtitle, TString ytitle, float overlayScale);

	//defines directory where plots are stored in, and in which the root file with the histograms is
	TString outputdir = "MT2Studies/CSVM40/new/HT/";
//	TString outputdir = "MT2Studies/CSVM40/new/MET/";
	TString outputdirold = outputdir;

TString outputname                = "HistogramsPVMinDPhiNJets.root";//this is the file with the histograms from MakeStudies.C

Bool_t  logflag                   = true;//plot y-axis in logstyle, default  = true
Bool_t  fSave                     = true;//save plots, default = true
Bool_t  plotonlywithratio         = true;//do only the plot with the ratio, alternative there is also a plot without ratio, default = true

//takes the root file from MakeStudies.C and creates plots
//you can choose plots by (un)commenting the histogram definition and 'minmetjetdphi' definitions below
void MakeStudiesPlotsFast(){

gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

    map<string, TH1D*> histos;
    map<string, THStack*> stacks;
    TLegend *Legend1 = new TLegend(.71,.54,.91,.92);
    vector<string> histonames; histonames.clear();
    Legend1 -> SetFillColor(0);
    Legend1 -> SetBorderSize(0);

    TFile *oldfile = TFile::Open(outputdirold + outputname);


	//defines all regions and cutvariables you want to study
	const int sampletypesize = 9;
	string sample_type[sampletypesize] = {"QCD", "WJets", "ZJets", "TTbar", "SingleTop", "Other", "mc", "susy", "data"};//same as type in samples.dat

	const int signalregionsize = 9;
	string signalregions[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b"};
	int NMT2bins[signalregionsize]     = { 40,  30,  40,  30,  25,  30,  25,  25,  25};
	double MT2upedge[signalregionsize] = {800, 600, 800, 600, 500, 600, 500, 500, 500};

	const int controlsize = 2;
	string control[controlsize] = {"had", "1l"};

	const int mindphisize = 5;//4;//5;//3;
	//string minmetjetdphi[mindphisize] = {"0p3", "0p4", "0p5"};//MinPhi
	//string minmetjetdphi[mindphisize] = {"le50", "le70", "le90", "le120", "le150"};//VSPT
	//string minmetjetdphi[mindphisize] = {"NoPass", "Pass20", "Pass40", "Pass50"};//PassJetID
	string minmetjetdphi[mindphisize] = {"2", "3", "4", "5", "6"};//MinPhiNJets

	const int nverticessize = 4;
	string numvert[nverticessize] = {"0to12V","13to19V","20upV","all"};

	//new
	const int jetthresholdsize = 3;//4;
	string jetthreshold[jetthresholdsize] = {"j20"/*, "j30"*/, "j40", "j50"};
	string mapname;

	//uncomment the variable you want to study
	for(int js = 0; js<signalregionsize; ++js){
	for(int ks = 0; ks<controlsize;      ++ks){
	for(int ls = 0; ls<mindphisize;      ++ls){
	for(int ms = 0; ms<nverticessize;    ++ms){
	for(int ns = 0; ns<jetthresholdsize; ++ns){
		if(ns!=0) continue;//use only MT2j20
		if(ms!=3) continue;//use only all vertices
		if(ks!=0) continue;//use only hadronic sample
		string hs = "_" + minmetjetdphi[ls] + "_" + signalregions[js] + "_" + control[ks] + "_" + numvert[ms] + "_";// + "_" + sample_type[is];
	/*	//MinDPhi
		mapname = "MT2" + jetthreshold[ns] + "_NoMinDPhi"; if(av) vs.push_back(mapname+hs);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhi"; if(av) vs.push_back(mapname+hs);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhi4"; if(av) vs.push_back(mapname+hs);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhiPt40"; if(av) vs.push_back(mapname+hs);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhi4Pt40"; if(av) vs.push_back(mapname+hs);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhiPt50"; if(av) vs.push_back(mapname+hs);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhi4Pt50"; if(av) vs.push_back(mapname+hs);
	*/
	/*	//VSPT
		mapname = "MT2" + jetthreshold[ns] + "_noVSPT"; if(av) vs.push_back(mapname+hs);
		mapname = "MT2" + jetthreshold[ns] + "_VSPTj20"; if(av) vs.push_back(mapname+hs);
		mapname = "MT2" + jetthreshold[ns] + "_VSPTj40"; if(av) vs.push_back(mapname+hs);
		mapname = "MT2" + jetthreshold[ns] + "_VSPTj50"; if(av) vs.push_back(mapname+hs);
	*/
	/*	//PassJetID
		mapname = "MT2" + jetthreshold[ns] + "_JetIDPt"; if(av) vs.push_back(mapname+hs);
	*/
		//MinDPhiNJets
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhi20NJets"; histonames.push_back(mapname+hs);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhi40NJets"; histonames.push_back(mapname+hs);

	}}}}}

	//load histograms
	for(int is = 0; is<sampletypesize; ++is){
		string hs = sample_type[is];
		for(unsigned int in = 0; in<histonames.size(); ++in){
			string name = histonames[in] + hs;
			if(histos.count(name)==0) histos[name ] = (TH1D*)oldfile->Get((name).c_str() );
		}
	}
	//stack for plotting
	for(unsigned int n = 0; n<histonames.size(); ++n){
			string h = "data";
			if(stacks.count(histonames[n]+h)==0) stacks[(histonames[n])+h ] = new THStack((histonames[n]+h).c_str(),(histonames[n]+h).c_str());
			for(int is = 0; is<sampletypesize;   ++is){
				if(is==8 || is==7 || is==6) continue;
				stacks[histonames[n]+h]->Add(histos[histonames[n]+sample_type[is] ]);
			}
	}
	bool leggy = true;
	for(unsigned int n = 0; n<histonames.size(); ++n){
		if(leggy){
	   	for(int is = 0; is<sampletypesize; ++is){
			string hs = sample_type[is];
			if(sample_type[is]=="mc")             continue;
			if(histos.count(histonames[n]+hs)==0) continue;
			if(sample_type[is]!="data") Legend1->AddEntry(histos[histonames[n]+hs], (sample_type[is]).c_str(), "f");
			else Legend1->AddEntry(histos[histonames[n]+hs], (sample_type[is]).c_str(), "l");
			leggy = false;
		}
		}
	}

	cout << "plotting histograms ..." << endl;
	MakePlots(histos, signalregionsize, signalregions, controlsize, control, mindphisize, minmetjetdphi, stacks, Legend1, histonames, outputdir);

}

//function to invoke Make1DPlotsRatio(...)
void MakePlots(map<string, TH1D*> histos, const int signalregionsize, string *signalregions, const int controlsize, string *control, const int mindphisize, string *minmetjetdphi, map<string, THStack*> stacks, TLegend *Legend1, vector<string> histonames, TString outputdirs){


	for(unsigned int n = 0; n<histonames.size(); ++n){//plot 1d
			string name = histonames[n];
			string h   = string("data");
			string hs1 = string("mc");
			string hs2 = string("data");
			string hs3 = string("susy");
			if(histos.count(name+hs1)==0) continue;
			TString ytitle = "Events / 20 GeV";//histos[name+hs1]->GetYaxis()->GetTitle();
			TString xtitle = histos[name+hs1]->GetXaxis()->GetTitle();

			TString outname = name + hs3 + (logflag ? "_log" : "") + "_overlay";

		if(histos[name+hs2]->Integral()>0 || histos[name+hs1]->Integral()>0){
			if(!plotonlywithratio) Make1DPlotsNoRatio(stacks[name+h], histos[name+hs1], histos[name+hs2], histos[name+hs3], logflag, false, outname, Legend1, xtitle, ytitle, 1.);

			Make1DPlotsRatio(stacks[name+h], histos[name+hs1], histos[name+hs2], histos[name+hs3], logflag, false, outname, Legend1, xtitle, ytitle, 1.);

		}
	}

}

//a copy of a MassPlotter.cc function - no detailed comments
void Make1DPlotsRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, TH1D* histsusy, bool logFlag, bool normalize, TString outname, TLegend *leg, TString xtitle, TString ytitle, float overlayScale){
	// define canvas and pads 
	TH1D *h1 = (TH1D*)histmc->Clone("h1_copy");
	TH1D *h2 = (TH1D*)histdata->Clone("h2_copy");
	TH1D *h3 = (TH1D*)histsusy->Clone("h3_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	h1->SetMarkerStyle(20);
	h2->SetMarkerStyle(20);

	TCanvas* c1 = new TCanvas(outname+"_ratio","",0,0,600,600);
	c1->SetFrameLineWidth(1);
	c1 -> cd();

	float border = 0.2;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad(outname+"_plotpad",  "Pad containing the overlay plot", 0,0.211838,1,1);
	p_plot->SetLeftMargin(0.131579);
	p_plot->SetRightMargin(0.08);
	p_plot->SetTopMargin(0.06895515);
	p_plot->SetBottomMargin(0.07206074);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad(outname+"_ratiopad", "Pad containing the ratio",   0,0.01863354,0.9967105,0.2189441);
	p_ratio->SetLeftMargin(0.1336634);	
	p_ratio->SetRightMargin(0.075);
	p_ratio->SetTopMargin(0.06976745);
	p_ratio->SetBottomMargin(0.2790698);
	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	if(logFlag) gPad->SetLogy(1);
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
	hstack->SetTitle("");

	h1    ->SetTitle("");
	h2    ->SetTitle("");
	h3    ->SetTitle("");
	hstack->SetMinimum(0.02);
	hstack->Draw("hist");

	h2    ->Draw("sameE");
	h3->Scale(overlayScale ? overlayScale : h2->Integral() / h3->Integral());
	h3->SetFillColor(0);
	h3->SetLineStyle(kDotted);
	h3->Draw("samehist");
	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.03);
	TString text ="";
	text = outname;
	TitleBox.DrawLatex(0.13,0.943,text.Data());
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
 
	TH1D *h_ratio = (TH1D*)histdata->Clone("h2_copy_2");
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);
 	h_ratio ->Divide(h2, h1);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->SetTitle("");	
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	//MC with errors
	TH1D*h_ratio_mc = (TH1D*)histmc->Clone("h1_copy_2");
	h_ratio_mc->Divide(h1);
	h_ratio_mc->SetTitle("");	
	h_ratio_mc->GetYaxis()->SetRangeUser(0,2);
	h_ratio_mc->GetXaxis()->SetLabelSize( 0.);
	h_ratio_mc->GetYaxis()->SetTitle("Data / MC");
	h_ratio_mc->GetXaxis()->SetTitle(xtitle);
	h_ratio_mc->GetXaxis()->SetTitleSize(0.2);
	h_ratio_mc->GetXaxis()->SetTitleOffset(0.5);
	h_ratio_mc->GetYaxis()->SetLabelSize(0.19);
	h_ratio_mc->GetXaxis()->SetTickLength(0.09);
	h_ratio_mc->GetYaxis()->SetTitleSize(0.18);
	h_ratio_mc->GetYaxis()->SetTitleOffset(0.36);
	h_ratio_mc->GetYaxis()->SetNdivisions(509);//xxxnew
	
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

	TString save=outname+"_ratio";
	if(fSave)Util::Print(c1, save, outputdir);


}

//a copy from a MassPlotter.cc function - no detailed comments
void Make1DPlotsNoRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, TH1D* histsusy, bool logFlag, bool normalize, TString outname, TLegend *leg, TString xtitle, TString ytitle, float overlayScale){

	TCanvas *col = new TCanvas(outname, "", 0, 0, 600, 600);
	col->SetRightMargin(0.08);
	col->SetLeftMargin(0.18);
	col->SetTopMargin(0.07);
	col->SetBottomMargin(0.17);

	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logFlag) {
		gPad->SetLogy(1);
		hstack     -> SetMinimum(0.05);
		histmc     -> SetMinimum(0.05);
		histdata   -> SetMinimum(0.05);
		histsusy   -> SetMinimum(0.05);
	}else{
		hstack->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = histdata->GetMaximum();
	double max2 = histmc->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	histdata  ->SetMaximum(max);
	histmc    ->SetMaximum(max);
	hstack    ->SetMaximum(max);

	hstack->Draw("hist");
	if(histdata->Integral()>0) {
		histdata       ->Draw("sameE");
	}
	histsusy->Scale(overlayScale ? overlayScale : histdata->Integral() / histsusy->Integral());
	histsusy->SetLineStyle(kDotted);
	histsusy->SetFillColor(0);
	histsusy->Draw("samehist");
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
	TitleBox.SetTextSize(0.0305);
	TString text = outname;
	TitleBox.DrawLatex(0.18,0.943,text.Data());


	hstack->GetXaxis()->SetTitle(xtitle);
	hstack->GetXaxis()->SetLabelSize(0.05);
	hstack->GetXaxis()->SetTitleSize(0.05);
	hstack->GetXaxis()->SetTitleOffset(1.1);
	hstack->GetYaxis()->SetTitle(ytitle);
	hstack->GetYaxis()->SetLabelSize(0.05);
	hstack->GetYaxis()->SetTitleSize(0.05);
	hstack->GetYaxis()->SetTitleOffset(1.47);

	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, outname, outputdir);
	if(fSave)Util::PrintEPS(col, outname, outputdir);

}
