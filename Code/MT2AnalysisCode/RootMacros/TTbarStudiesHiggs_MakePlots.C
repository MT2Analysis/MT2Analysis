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
#include "TPad.h"
#include "TLine.h"
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
#include "TLatex.h"
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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/Utilities.hh"

//run via root -l -b -q TTbarStudiesHiggs_MakePlots.C++

//this function is basically a copy of TTbarStudies_MakePlots.C
//but shapes are plotted for Mbb instead of MT2, event selection is the MT2 Higgs selection instead of inclusive MT2 analysis
//therefore no additional comments are given


using namespace std;

void MakePlot(int numhistos, TH1D *h0_orig, TH1D *h1_orig, TH1D *h2_orig, TH1D *h3_orig, TLegend *leg, bool logflag, TString titlebox, TString name, TString outputdir, bool normalized);
void TTbarStudiesHiggs_MakePlots();

const int fVerbose = 3;

const int sampletypesize = 34;
string sample_type[sampletypesize] = {"Nominal", "MatchingUp", "MatchingDown", "ScaleUp", "ScaleDown", "MassUp", "MassDown", "JESUp", "JESDown", "METUp", "METDown", "BSFUp", "BSFDown", "ISRUp", "ISR", "ISRDown", "PUUp", "PUDown", "TopPtReweighted", "nominal0l", "nominal1l", "nominal2l", "nominal1Sample", "nominalPowheg", "WUp", "WDown", "TopUp", "TopDown", "SingleTopUp", "SingleTopDown", "TTVUp", "TTVDown", "MT2relaxedcut", "MT2fullcut"};//maybe not all filled  - ! include TTbar and WJets
const int samplekindsize = 6;
string sample_kind[samplekindsize] = {"allMC","allTop", "TTbar", "SingleTop", "TTV", "WJets"};
const int HTregionsize = 2;
string HT_region[HTregionsize] = {"lowHT", "highHT"};

bool fSave  = true;//save the plots

void TTbarStudiesHiggs_MakePlots(){

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	bool onlyttbar = false;//plot shape only for ttbar, default = false, if true: onlytop = false, onlyW = false
	bool onlytop   = false;//plot shape only for top(ttbar+single top), default = false, if true: onlyW = false
	bool onlyW     = false;//plot shape only for W, default = false

	bool normalized = true; //plot normalized distribution, default = true
	bool notnorm    = false;//plot not-normalized distribution, default = false
	bool logflag    = false;//make y-axis log style, default = false (as only signal region is considered)

	bool WnoScaleMatching = true;//see in TTbarStudiesHiggs.C
	bool ISRdefault       = true;//see in TTbarStudiesHiggs.C (fISRreweight)

	//number of systematics, and which systematics you want to plot
	const int numberdifferentplots = 13;
	string differentplots[numberdifferentplots] = {"Nominal", "differentMC", "JES", "MET", "Scale", "Matching", "BSF", "PU" , "ISR", "Wxs", "Topxs", "MT2cuts", "bla"};
	bool   whattoplot[numberdifferentplots]     = {false,     true,          true,  true,  true,    true,       true,  true,   true, true,  true,    true,          false};

	TString   inputdir = "Filtered/TTbarStudies/Higgs/";
	TString  inputname = "TTbarStudiesHistograms_all_ISR_limitMbb.root";
	if(WnoScaleMatching)  inputname = "TTbarStudiesHistograms_all_noWscaleupdown_ISR_limitMbb.root";
	if(!ISRdefault){
		inputname = "TTbarStudiesHistograms_all_limitMbb.root";
		if(WnoScaleMatching)  inputname = "TTbarStudiesHistograms_all_noWscaleupdown_limitMbb.root";
	}

	TString outputdir  = "Filtered/TTbarStudies/Higgs/Plots/all/ISR/Mbblimited/";
	if(WnoScaleMatching) outputdir  = "Filtered/TTbarStudies/Higgs/Plots/all_noWscaleupdown/ISR/Mbblimited/";
	if(!ISRdefault){
		outputdir  = "Filtered/TTbarStudies/Higgs/Plots/all/Mbblimited/";
		if(WnoScaleMatching) outputdir  = "Filtered/TTbarStudies/Higgs/Plots/all_noWscaleupdown/Mbblimited/";
	}
	if(onlyttbar)    outputdir = outputdir + "ttbar/";
	else if(onlytop) outputdir = outputdir + "top/";
	else if(onlyW)   outputdir = outputdir + "W/";
	if(logflag) outputdir = outputdir + "log/";
	else        outputdir = outputdir + "linear/";
    	Util::MakeOutputDir(outputdir);

	TFile *infile = TFile::Open(inputdir + inputname);

	map<string, TH1D*>    histos;
	map<string, TH1D*>    histosnormalized;

	string samplekind;
	if(onlyttbar)    samplekind = "TTbar";
	else if(onlytop) samplekind = "allTop";
	else if(onlyW)   samplekind = "WJets";
	else             samplekind = "allMC";

	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		string hs = string("_") + sample_type[i1] + string("_") + samplekind + string("_") + HT_region[i3];// + string("_") + signal_region[i2];
		string mapname = "MT2" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)infile->Get(mapname.c_str());
		string mapnamenorm = "Norm"+mapname;
		if(histosnormalized.count(mapnamenorm) == 0 ) histosnormalized[mapnamenorm] = (TH1D*)infile->Get(mapnamenorm.c_str());
	}}

	for(int i1 = 0; i1<numberdifferentplots; ++i1){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		if(whattoplot[i1] == false) continue;
		string hs0, hs1, hs2, hs3; int numhistos;
		string hs = string("_") + samplekind + string("_") + HT_region[i3];
		hs0 = "_Nominal"+hs;
		if(!WnoScaleMatching){
			if(differentplots[i1]=="Nominal"||differentplots[i1]=="JES"||differentplots[i1]=="MET"||differentplots[i1]=="BSF"||differentplots[i1]=="PU"||differentplots[i1]=="ISR"||differentplots[i1]=="Wxs"||differentplots[i1]=="Topxs"||differentplots[i1]=="SingleTopxs"||differentplots[i1]=="Restxs") continue;
		}
		if(differentplots[i1]=="Nominal")     { numhistos = 0; }
		if(differentplots[i1]=="differentMC") { numhistos = 2; hs1 = "_nominal1Sample"+hs; hs2 = "_nominalPowheg"+hs; }
		if(differentplots[i1]=="JES")         { numhistos = 2; hs1 = "_JESUp"+hs;          hs2 = "_JESDown"+hs;       }
		if(differentplots[i1]=="MET")         { numhistos = 2; hs1 = "_METUp"+hs;          hs2 = "_METDown"+hs;       }
		if(differentplots[i1]=="Scale")       { numhistos = 2; hs1 = "_ScaleUp"+hs;        hs2 = "_ScaleDown"+hs;     }
		if(differentplots[i1]=="Matching")    { numhistos = 2; hs1 = "_MatchingUp"+hs;     hs2 = "_MatchingDown"+hs;  }
		if(differentplots[i1]=="BSF")         { numhistos = 2; hs1 = "_BSFUp"+hs;          hs2 = "_BSFDown"+hs;       }
		if(differentplots[i1]=="PU")          { numhistos = 2; hs1 = "_PUUp"+hs;           hs2 = "_PUDown"+hs;        }
		if(differentplots[i1]=="ISR")         { numhistos = 3; hs1 = "_ISRUp"+hs;          hs2 = "_ISRDown"+hs;       hs3 = "_ISR"+hs; }
		if(differentplots[i1]=="Wxs")         { numhistos = 2; hs1 = "_WUp"+hs;            hs2 = "_WDown"+hs;         }
		if(differentplots[i1]=="Topxs")       { numhistos = 2; hs1 = "_TopUp"+hs;          hs2 = "_TopDown"+hs;       }
		if(differentplots[i1]=="MT2cuts")     { numhistos = 2; hs1 = "_MT2relaxedcut"+hs;    hs2 = "_MT2fullcut"+hs; }
		if(differentplots[i1]=="Restxs")      { numhistos = 2; hs1 = "_TTVUp"+hs;          hs2 = "_TTVDown"+hs;       }
		if(numhistos<3) hs3 = "_Nominal"+hs;//dummy
		if(numhistos<2) hs2 = "_Nominal"+hs;//dummy
		if(numhistos<1) hs1 = "_Nominal"+hs;//dummy

		//only thing to be set is the color
		histos["MT2"+hs0]->SetLineColor(kBlack); histos["MT2"+hs0]->SetMarkerColor(kBlack);
		if(numhistos>0){
		histos["MT2"+hs1]->SetLineColor(kRed);   histos["MT2"+hs1]->SetMarkerColor(kRed);
		histos["MT2"+hs2]->SetLineColor(kBlue);  histos["MT2"+hs2]->SetMarkerColor(kBlue);
		}
		if(numhistos>2) histos["MT2"+hs3]->SetLineColor(kGreen+2);   histos["MT2"+hs3]->SetMarkerColor(kGreen+2);
		histosnormalized["NormMT2"+hs0]->SetLineColor(kBlack); histosnormalized["NormMT2"+hs0]->SetMarkerColor(kBlack); histosnormalized["NormMT2"+hs0]->GetYaxis()->SetTitle("Events (norm.)");
		if(numhistos>0){
		histosnormalized["NormMT2"+hs1]->SetLineColor(kRed);   histosnormalized["NormMT2"+hs1]->SetMarkerColor(kRed); histosnormalized["NormMT2"+hs1]->GetYaxis()->SetTitle("Events (norm.)");
		histosnormalized["NormMT2"+hs2]->SetLineColor(kBlue);  histosnormalized["NormMT2"+hs2]->SetMarkerColor(kBlue); histosnormalized["NormMT2"+hs2]->GetYaxis()->SetTitle("Events (norm.)");
		}
		if(numhistos>2) histosnormalized["NormMT2"+hs3]->SetLineColor(kGreen+2);   histosnormalized["NormMT2"+hs3]->SetMarkerColor(kGreen+2); histosnormalized["NormMT2"+hs3]->GetYaxis()->SetTitle("Events (norm.)");

		TLegend* Legend1 = new TLegend(.71,.68,.91,.92);
		string legendname = differentplots[i1]+hs+"_legend";
		Legend1->SetName(legendname.c_str());
		                                                     Legend1->AddEntry(histos["MT2"+hs0], "Nominal",               "lp");
		if(differentplots[i1]!="MT2cuts"){
		if(numhistos>0 && differentplots[i1]!="differentMC") Legend1->AddEntry(histos["MT2"+hs1], "+ 1 #sigma",            "lp");
		else if(numhistos>0)                                 Legend1->AddEntry(histos["MT2"+hs1], "Top MC Madgraph Incl.", "lp");
		if(numhistos>1 && differentplots[i1]!="differentMC") Legend1->AddEntry(histos["MT2"+hs2], "- 1 #sigma",            "lp");
		else if(numhistos>1)                                 Legend1->AddEntry(histos["MT2"+hs2], "Top MC Powheg",         "lp");
		if(numhistos>2)                                      Legend1->AddEntry(histos["MT2"+hs3], "std reweighting",       "lp");
		} else {
		if(i3==1) Legend1->AddEntry(histos["MT2"+hs1], "M_{T2} > 50 GeV",            "lp");
		if(i3==1) Legend1->AddEntry(histos["MT2"+hs2], "M_{T2} > 125 GeV",            "lp");
		if(i3==0) Legend1->AddEntry(histos["MT2"+hs1], "M_{T2} > 75 GeV",            "lp");
		if(i3==0) Legend1->AddEntry(histos["MT2"+hs2], "M_{T2} > 200 GeV",            "lp");
		}

		TString name = differentplots[i1];
		name.ReplaceAll("xs", " cross section");
		if(differentplots[i1]!="Nominal") name = name + " systematics";
		if(samplekind!="allMC") name = name + ", " + samplekind;
		name = name + ", " + HT_region[i3];

		TString filename = differentplots[i1]+"_"+samplekind+"_"+ HT_region[i3];

		if(notnorm){
			TString thisoutputdir = outputdir + "notnorm/";
    			Util::MakeOutputDir(thisoutputdir);
			MakePlot(numhistos, histos["MT2"+hs0], histos["MT2"+hs1], histos["MT2"+hs2], histos["MT2"+hs3], Legend1, logflag, name, filename, thisoutputdir, false);
		}
		if(normalized){
			TString thisoutputdir = outputdir + "normalized/";
    			Util::MakeOutputDir(thisoutputdir);
			MakePlot(numhistos, histosnormalized["NormMT2"+hs0], histosnormalized["NormMT2"+hs1], histosnormalized["NormMT2"+hs2], histosnormalized["NormMT2"+hs3], Legend1, logflag, name, filename, thisoutputdir, true);
		}
	}}//}

}//TTbarStudiesPlots()

void MakePlot(int numhistos, TH1D *h0_orig, TH1D *h1_orig, TH1D *h2_orig, TH1D *h3_orig, TLegend *leg, bool logflag, TString titlebox, TString name, TString outputdir, bool normalized){

	TH1D *h0 = (TH1D*)h0_orig->Clone("h0");
	TH1D *h1 = (TH1D*)h1_orig->Clone("h1");
	TH1D *h2 = (TH1D*)h2_orig->Clone("h2");
	TH1D *h3 = (TH1D*)h3_orig->Clone("h3");

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
 
 	p_plot ->cd();

	if(logflag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);

	double max1 = h0->GetMaximum();
	double max2 = h1->GetMaximum();
	double max3 = h2->GetMaximum();
	double max4 = h3->GetMaximum();
	double maxA  = (max1>max2)?max1:max2;
	double maxB  = (max3>max4)?max3:max4;
	double max   = (maxA>maxB)?maxA:maxB;
	if(logflag) max = 2.5*max;
	else max = 1.5*max;

	h0->SetMaximum(max);
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	h3->SetMaximum(max);

	if(normalized&&logflag) h0->SetMinimum(0.005);
	else if(       logflag) h0->SetMinimum(0.05);
	else                    h0->SetMinimum(0.0);

	h0->GetXaxis()->SetTitle("");
	h0->GetXaxis()->SetLabelSize(0.05);
	h0->GetYaxis()->SetLabelSize(0.05);
	h0->GetYaxis()->SetTitleSize(0.05);
	h0->GetYaxis()->SetTitleOffset(1.3);
	h0    ->Draw("");//nominal shape without error bars
	if(numhistos>0) h1    ->Draw("sameE");//up
	if(numhistos>1) h2    ->Draw("sameE");//down
	if(numhistos>2) h3    ->Draw("same");//ISRnominal / or notISR

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.03);
	TitleBox.DrawLatex(0.13,0.943,titlebox.Data());

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

	//MC with errors
	TH1D*h_ratio_mc = (TH1D*)h0_orig->Clone("h0_copy");
	h_ratio_mc->Divide(h0);

	if(name.Contains("Scale")||name.Contains("Matching")) h_ratio_mc->GetYaxis()->SetRangeUser(0.0,2.0);
	else h_ratio_mc->GetYaxis()->SetRangeUser(0.5,1.5);

	h_ratio_mc->GetXaxis()->SetLabelSize( 0.);
	h_ratio_mc->GetYaxis()->SetTitle("ratios");
	h_ratio_mc->GetXaxis()->SetTitleSize(0.2);
	h_ratio_mc->GetXaxis()->SetTitleOffset(0.5);
	h_ratio_mc->GetYaxis()->SetLabelSize(0.19);
	h_ratio_mc->GetXaxis()->SetTickLength(0.09);
	h_ratio_mc->GetYaxis()->SetTitleSize(0.18);
	h_ratio_mc->GetYaxis()->SetTitleOffset(0.36);
	h_ratio_mc->GetYaxis()->SetNdivisions(509);	
	h_ratio_mc->SetFillStyle(3001);
	h_ratio_mc->Draw("E2");

	if(numhistos>0){
		TH1D *h_ratio = (TH1D*)h1_orig->Clone("h1_copy");	
		h_ratio ->SetStats(0);
		h_ratio ->SetMarkerStyle(20);
		h_ratio ->Divide(h1, h0);
		h_ratio ->SetMinimum(0.4);
		h_ratio ->SetMaximum(3.0);
		h_ratio ->GetYaxis()->SetTitleOffset(h0->GetYaxis()->GetTitleOffset());
		h_ratio ->DrawCopy("Esame");//LEO MOD
	}
	if(numhistos>1){
		TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");	
		h_ratio ->SetStats(0);
		h_ratio ->SetMarkerStyle(20);
		h_ratio ->Divide(h2, h0);
		h_ratio ->SetMinimum(0.4);
		h_ratio ->SetMaximum(3.0);
		h_ratio ->GetYaxis()->SetTitleOffset(h0->GetYaxis()->GetTitleOffset());
		h_ratio ->DrawCopy("Esame");//LEO MOD
	} 
	if(numhistos>2){
		TH1D *h_ratio = (TH1D*)h3_orig->Clone("h3_copy");	
		h_ratio ->SetStats(0);
		h_ratio ->SetMarkerStyle(20);
		h_ratio ->Divide(h3, h0);
		h_ratio ->SetMinimum(0.4);
		h_ratio ->SetMaximum(3.0);
		h_ratio ->GetYaxis()->SetTitleOffset(h0->GetYaxis()->GetTitleOffset());
		h_ratio ->DrawCopy("same");//LEO MOD
	} 

	TLine *l3 = new TLine(h0->GetXaxis()->GetXmin(), 1.00, h0->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name;
	if(fSave)Util::Print(c1, save, outputdir, 0);	
	if(fSave)Util::PrintPDF(c1, save, outputdir);

}
