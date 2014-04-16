#include "TEventList.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
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
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"//use your own path

//run via root -l -b -q fit_lowercut_uppercut_ratios.C++

void fit_lowercut_uppercut_ratios();

const int gNbins_lowHT_lowMT2                   = 13;
const double  gbins_lowHT_lowMT2[gNbins_lowHT_lowMT2+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 100, 125, 150, 200, 300, 500};

const int gNbins_highHT_lowMT2                   = 13;
const double  gbins_highHT_lowMT2[gNbins_highHT_lowMT2+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 100, 125, 150, 180, 260, 500};


const int gNbins_highMT2                   = 15;
const double  gbins_highMT2[gNbins_highMT2+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 100, 125, 150, 200, 275, 375, 500, 700};

//this code fits the ratio obtained from files of save_lowercut_uppercut_histos.C using TEfficiencies, that allow for correct uncertainty covering (as binominal uncertainties are 0, if e.g. two events out of two pass a selection cut)
//this code was needed as at this time this piece of code could not be used together with save_lowercut_uppercut_histos.C due to conflicting ROOT copies needed.
//flags common with save_lowercut_uppercut_histos.C are not commented
//later this macro was combined with save_lowercut_uppercut_histos.C ==> see run_draw_ratio.C
void fit_lowercut_uppercut_ratios(){
	gROOT->ProcessLine(".x SetStyle_PRD.C");
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);

	bool lowMT2  = false;
	bool highMT2 = true;
	
	bool lowHT   = false;
	bool highHT  = false;
	bool allHT   = true;

	bool ele     = true;
	bool muo     = false;

	bool fitting = true; //perform the fitting, if false - only a plot will be done, default = true

	TString outputdir = "LostLepton/fits";//newfits if lepton-corrected quantities, fits if normal quantities
	TString outputfile = "MinDPhiRatio_";
	outputfile = "MinDPhiRatio_";//new
	outputdir = "LostLepton/lepcortest";//new
	if(ele) outputfile = outputfile + (TString)"ele_";
	if(muo) outputfile = outputfile + (TString)"muo_";
	if(lowMT2)  outputfile = outputfile + (TString)"MT2b_";
	if(highMT2) outputfile = outputfile + (TString)"MT2_";
	if(lowHT)  outputfile = outputfile + (TString)"lowHT";
	if(highHT) outputfile = outputfile + (TString)"highHT";
	if(allHT)  outputfile = outputfile + (TString)"allHT";

	TString outputfile2 = outputfile + (TString)"_fitted";
	TString outputdir2 = outputdir +(TString)"/newfits";
  	Util::MakeOutputDir(outputdir2);


	double lowerfitborder = 125.;
	if(lowMT2)  lowerfitborder = 125.;
	if(highMT2) lowerfitborder = 150.;
	double upperfitborder = 1000.;
	if(lowMT2 && highHT) upperfitborder = gbins_highHT_lowMT2[gNbins_highHT_lowMT2];//high HT, MT2b
	if(lowMT2 && (lowHT||allHT)) upperfitborder = gbins_lowHT_lowMT2[gNbins_lowHT_lowMT2];//low HT, MT2b
	if(highMT2) upperfitborder = gbins_highMT2[gNbins_highMT2];//both HT, MT2

       TFile *oldfile = TFile::Open(outputdir + (TString)"/" + outputfile + (TString)".root");
	//load the file/histograms from save_lowercut_uppercut_histos.C
  	TH1D* h_lower_band_mc_old	= (TH1D*)oldfile->Get("h_lower_band_mc");
  	TH1D* h_upper_band_mc_old	= (TH1D*)oldfile->Get("h_upper_band_mc");
  	TH1D* h_ratio_mc_old		= (TH1D*)oldfile->Get("h_ratio_mc");
  	TH1D* h_lower_band_susy_old	= (TH1D*)oldfile->Get("h_lower_band_susy");
  	TH1D* h_upper_band_susy_old	= (TH1D*)oldfile->Get("h_upper_band_susy");
  	TH1D* h_ratio_susy_old		= (TH1D*)oldfile->Get("h_ratio_susy");
  	TH1D* h_lower_band_data_old	= (TH1D*)oldfile->Get("h_lower_band_data");
  	TH1D* h_upper_band_data_old	= (TH1D*)oldfile->Get("h_upper_band_data");
  	TH1D* h_ratio_data_old		= (TH1D*)oldfile->Get("h_ratio_data");
  	TH1D* h_axis_old		= (TH1D*)oldfile->Get("h_axis");

	TString title = "MinMetJetDPhi ratio";
	if(lowHT)  title = title + (TString)", low HT";
	if(highHT) title = title + (TString)", high HT";
	if(ele)    title = title + (TString)", 1 Electron";
	if(muo)    title = title + (TString)", 1 Muon";

	TString title2 = "MT2 with or w/o MinMetJetDPhi";
	if(lowHT)  title2 = title2 + (TString)", low HT";
	if(highHT) title2 = title2 + (TString)", high HT";
	if(ele)    title2 = title2 + (TString)", 1 Electron";
	if(muo)    title2 = title2 + (TString)", 1 Muon";

	h_axis_old->SetTitle(title);

	h_lower_band_mc_old->SetTitle(title);
	h_upper_band_mc_old->SetTitle(title);
	h_lower_band_susy_old->SetTitle(title);
	h_upper_band_susy_old->SetTitle(title);
	h_lower_band_data_old->SetTitle(title);
	h_upper_band_data_old->SetTitle(title);

	//make clones used for TEfficiencies
	TH1D *h_low_mc   = (TH1D*)h_lower_band_mc_old  ->Clone("h_low_mc"  );
	TH1D *h_up_mc    = (TH1D*)h_upper_band_mc_old  ->Clone("h_up_mc"   );
	TH1D *h_low_susy = (TH1D*)h_lower_band_susy_old->Clone("h_low_susy");
	TH1D *h_up_susy  = (TH1D*)h_upper_band_susy_old->Clone("h_up_susy" );
	TH1D *h_low_data = (TH1D*)h_lower_band_data_old->Clone("h_low_data");
	TH1D *h_up_data  = (TH1D*)h_upper_band_data_old->Clone("h_up_data" );

	//get the efficiency
  TEfficiency* e_ratio_mc   = new TEfficiency((*h_low_mc),   (*h_up_mc)  ); 
  TEfficiency* e_ratio_susy = new TEfficiency((*h_low_susy), (*h_up_susy)); 
  TEfficiency* e_ratio_data = new TEfficiency((*h_low_data), (*h_up_data)); 

  TH1D *h_ratio_mcnew;
  if(lowMT2 && !(highHT)) h_ratio_mcnew = new TH1D("h_ratio_mcnew",        outputfile,        gNbins_lowHT_lowMT2, gbins_lowHT_lowMT2);
  if(lowMT2 && highHT)    h_ratio_mcnew = new TH1D("h_ratio_mcnew",        outputfile,        gNbins_highHT_lowMT2, gbins_highHT_lowMT2);
  if(highMT2)             h_ratio_mcnew = new TH1D("h_ratio_mcnew",        outputfile,        gNbins_highMT2, gbins_highMT2);
  TGraphAsymmErrors *g_ratio = new TGraphAsymmErrors(h_ratio_mcnew);
  g_ratio->SetNameTitle("g_ratio", "g_ratio");

	//for plotting purposes (and at that time also fitting purposes) store TEfficiency results into a TGraphAsymmErrors
  for(int n =1; n<=h_ratio_mcnew->GetNbinsX(); ++n){
	double x,y;
	g_ratio->GetPoint(n-1, x, y);
	g_ratio->SetPoint(n-1, x, e_ratio_mc->GetEfficiency(n));
	g_ratio->SetPointEYhigh(n-1, e_ratio_mc->GetEfficiencyErrorUp(n));
	g_ratio->SetPointEYlow(n-1, e_ratio_mc->GetEfficiencyErrorLow(n));
  }

  e_ratio_mc  ->SetNameTitle("e_ratio_mc",   outputfile);
  e_ratio_data->SetNameTitle("e_ratio_data", outputfile);
  e_ratio_susy->SetNameTitle("e_ratio_susy", outputfile);
  e_ratio_mc              ->SetLineColor(1);
  e_ratio_susy            ->SetLineColor(4);
  e_ratio_data            ->SetLineColor(2);

  //try different fit functions
  TF1 *fitfunc = new TF1("fit", "pol0", lowerfitborder, upperfitborder);//defined 110
  TF1 *fitfunc1 = new TF1("fit1", "pol1", lowerfitborder, upperfitborder);//defined 110
  TF1 *fitfunc2 = new TF1("fit2", "[0]*(1 - exp(-(x*[1])))", lowerfitborder, upperfitborder);//defined 110
  TF1 *fitfunc3 = new TF1("fit3", "[0]*(1 - exp(-(x*[1])))", 0, upperfitborder);//defined 110
  TF1 *fitfunc2b = new TF1("fit2b", "[0]*(1 - [2]*exp(-(x*[1])))", lowerfitborder, upperfitborder);//defined 110
  TF1 *fitfunc3b = new TF1("fit3b", "[0]*(1 - [2]*exp(-(x*[1])))", 0, upperfitborder);//defined 110

  if(fitting){
	//perform the fitting
	double mean;
	int lbin = h_ratio_mc_old->FindBin(lowerfitborder+0.01);//+0.01 to exclude border effects
	int ubin = h_ratio_mc_old->FindBin(upperfitborder-0.01);
	mean = (h_ratio_mc_old->Integral(lbin, ubin))/(ubin-lbin+1);//mean value of ratio for MT2>100
	fitfunc->SetParameter(0, mean);
	fitfunc->SetParLimits(0, 0.,1.);
	fitfunc->SetLineColor(6);
	fitfunc->SetLineWidth(2);
	fitfunc1->SetParameter(0, mean);
	fitfunc2->SetParameter(0, mean);
	fitfunc3->SetParameter(0, mean);
	fitfunc2b->SetParameter(0, mean);
	fitfunc3b->SetParameter(0, mean);
	fitfunc1->SetParameter(1, 0);
	fitfunc2->SetParameter(1, 0);
	fitfunc3->SetParameter(1, 0);
	fitfunc2b->SetParameter(1, 0);
	fitfunc3b->SetParameter(1, 0);
	fitfunc2b->SetParameter(2, 1);
	fitfunc3b->SetParameter(2, 1);
	cout << endl << "fit horizontal line" << endl;
	g_ratio->Fit(fitfunc, "R");
	cout << endl << "fit slopped line" << endl;
	g_ratio->Fit(fitfunc1, "R0+");
	cout << endl << "fit exponential a*(1-exp(-b*MT2))" << endl;
	g_ratio->Fit(fitfunc2, "R0+");
	cout << endl << "fit full-range exponential line a*(1-exp(-b*MT2))" << endl;
	g_ratio->Fit(fitfunc3, "R0+");
	cout << endl << "fit exponential a*(1-c*exp(-b*MT2))" << endl;
	g_ratio->Fit(fitfunc2b, "R0+");
	cout << endl << "fit full-range exponential line a*(1-c*exp(-b*MT2))" << endl;
	g_ratio->Fit(fitfunc3b, "R0+");
	cout << "fitting done" << endl << endl;
	cout << "Fit " << outputfile << " in range MT2(" << lowerfitborder << "," << upperfitborder << ")" << " (for MC)" << endl;
	cout << "Fit is a constant: c" << endl;
	cout << "Fit constant = " << fitfunc->GetParameter(0) << " +/- " << fitfunc->GetParError(0) << endl;
	cout << "Goodness of fit: chi^2/NDF = " << fitfunc->GetChisquare() << "/" << fitfunc->GetNDF() << " = " << (fitfunc->GetChisquare())/(fitfunc->GetNDF()) << endl;
	cout << "Fit probability: " << fitfunc->GetProb() << endl;
	cout << "Next Fit: Fit is a tilted line: a + c*MT2" << endl;
	cout << "Fit constant = " << fitfunc1->GetParameter(0) << " +/- " << fitfunc1->GetParError(0) << endl;
	cout << "Fit slope    = " << fitfunc1->GetParameter(1) << " +/- " << fitfunc1->GetParError(1) << endl;
	cout << "Goodness of fit: chi^2/NDF = " << fitfunc1->GetChisquare() << "/" << fitfunc1->GetNDF() << " = " << (fitfunc1->GetChisquare())/(fitfunc1->GetNDF()) << endl;
	cout << "Fit probability: " << fitfunc1->GetProb() << endl;
	cout << "Next Fit: Fit is a exponential: a*(1-exp(-b*MT2))" << endl;
	cout << "Fit constant = " << fitfunc2->GetParameter(0) << " +/- " << fitfunc2->GetParError(0) << endl;
	cout << "Fit exp.fact = " << fitfunc2->GetParameter(1) << " +/- " << fitfunc2->GetParError(1) << endl;
	cout << "Goodness of fit: chi^2/NDF = " << fitfunc2->GetChisquare() << "/" << fitfunc2->GetNDF() << " = " << (fitfunc2->GetChisquare())/(fitfunc2->GetNDF()) << endl;
	cout << "Fit probability: " << fitfunc2->GetProb() << endl;
	cout << "Next Fit: Fit is a exponential in complete range: a*(1-exp(-b*MT2))" << endl;
	cout << "Fit constant = " << fitfunc3->GetParameter(0) << " +/- " << fitfunc3->GetParError(0) << endl;
	cout << "Fit exp.fact = " << fitfunc3->GetParameter(1) << " +/- " << fitfunc3->GetParError(1) << endl;
	cout << "Goodness of fit: chi^2/NDF = " << fitfunc3->GetChisquare() << "/" << fitfunc3->GetNDF() << " = " << (fitfunc3->GetChisquare())/(fitfunc3->GetNDF()) << endl;
	cout << "Fit probability: " << fitfunc3->GetProb() << endl;
	cout << "Next Fit: Fit is a exponential: a*(1-c*exp(-b*MT2))" << endl;
	cout << "Fit constant = " << fitfunc2b->GetParameter(0) << " +/- " << fitfunc2b->GetParError(0) << endl;
	cout << "Fit exp.fact = " << fitfunc2b->GetParameter(1) << " +/- " << fitfunc2b->GetParError(1) << endl;
	cout << "Fit expconst = " << fitfunc2b->GetParameter(2) << " +/- " << fitfunc2b->GetParError(2) << endl;
	cout << "Goodness of fit: chi^2/NDF = " << fitfunc2b->GetChisquare() << "/" << fitfunc2b->GetNDF() << " = " << (fitfunc2b->GetChisquare())/(fitfunc2b->GetNDF()) << endl;
	cout << "Fit probability: " << fitfunc2b->GetProb() << endl;
	cout << "Next Fit: Fit is a exponential in complete range: a*(1-c*exp(-b*MT2))" << endl;
	cout << "Fit constant = " << fitfunc3b->GetParameter(0) << " +/- " << fitfunc3b->GetParError(0) << endl;
	cout << "Fit exp.fact = " << fitfunc3b->GetParameter(1) << " +/- " << fitfunc3b->GetParError(1) << endl;
	cout << "Fit expconst = " << fitfunc3b->GetParameter(2) << " +/- " << fitfunc3b->GetParError(2) << endl;
	cout << "Goodness of fit: chi^2/NDF = " << fitfunc3b->GetChisquare() << "/" << fitfunc3b->GetNDF() << " = " << (fitfunc3b->GetChisquare())/(fitfunc3b->GetNDF()) << endl;
	cout << "Fit probability: " << fitfunc3b->GetProb() << endl;
  }
  //save everything
  TFile *ratio_file = new TFile(outputdir2 + (TString)"/" +  outputfile2 +TString(".root"), "RECREATE");
  ratio_file		->cd();
  h_lower_band_mc_old	->Write();
  h_upper_band_mc_old	->Write();
  h_ratio_mc_old	->Write();
  h_lower_band_susy_old	->Write();
  h_upper_band_susy_old	->Write();
  h_ratio_susy_old	->Write();
  h_lower_band_data_old	->Write();
  h_upper_band_data_old	->Write();
  h_ratio_data_old	->Write();
  e_ratio_data		->Write();
  e_ratio_susy		->Write();
  e_ratio_mc		->Write();
  fitfunc		->Write();
  fitfunc1		->Write();
  fitfunc2		->Write();
  fitfunc3		->Write();
  fitfunc2b		->Write();
  fitfunc3b		->Write();
  g_ratio		->Write();
  cout << endl << "saved in " << outputdir2 + (TString)"/" +  outputfile2 +TString(".root") << endl;


	//make the plots now - for now the saving option has been commented.
	double min = 0.;
	double max = 1.5;

	h_ratio_data_old  ->SetMaximum(max);
	h_ratio_mc_old->SetMaximum(max);
	h_ratio_data_old  ->SetMinimum(min);
	h_ratio_mc_old->SetMinimum(min);
	h_axis_old->SetMaximum(max);
	h_axis_old->SetMinimum(min);

	TLegend* Legend2 = new TLegend(.55,.6,.85,.85);
	Legend2     ->AddEntry(h_upper_band_mc_old, "MC relaxed cuts", "lp");
	Legend2     ->AddEntry(h_lower_band_mc_old, "MC full cuts",    "lp");
	Legend2     -> SetFillColor(0);
	Legend2     -> SetBorderSize(0);

  	TCanvas *col2 = new TCanvas("col2", outputfile+TString("_mc"), 0, 0, 900, 700);
	gPad->SetLogy(1);
  	h_upper_band_mc_old->Draw("E1");
  	h_lower_band_mc_old->Draw("E1same");
	Legend2->Draw();
	gPad->RedrawAxis();
	col2 ->Update();
//	Util::PrintNoEPS(col2, TString("Z_") + outputfile+TString("_mc"), outputdir2, false);
//	Util::PrintEPS(col2, TString("Z_") + outputfile+TString("_mc"), outputdir2);

	TLegend* Legend3 = new TLegend(.55,.6,.85,.85);
	Legend3     ->AddEntry(h_upper_band_data_old, "Data relaxed cuts", "lp");
	Legend3     ->AddEntry(h_lower_band_data_old, "Data full cuts",    "lp");
	Legend3     -> SetFillColor(0);
	Legend3     -> SetBorderSize(0);

  	TCanvas *col3 = new TCanvas("col3", outputfile+TString("_data"), 0, 0, 900, 700);
	gPad->SetLogy(1);
  	h_upper_band_data_old->Draw("E1");
  	h_lower_band_data_old->Draw("E1same");
	Legend3->Draw();
	gPad->RedrawAxis();
	col3 ->Update();
//	Util::PrintNoEPS(col3, TString("Z_") + outputfile+TString("_data"), outputdir2, false);
//	Util::PrintEPS(col3, TString("Z_") + outputfile+TString("_data"), outputdir2);

	TLegend* Legend1 = new TLegend(.6,.75,.85,0.85);
	Legend1     ->AddEntry(e_ratio_mc,   "ratio mc",   "lp");
	Legend1     ->AddEntry(e_ratio_data, "ratio data", "lp");
	Legend1     -> SetFillColor(0);
	Legend1     -> SetBorderSize(0);

  	TCanvas *col = new TCanvas("col", outputfile, 0, 0, 900, 700);
	gPad->SetLogy(0);
	h_axis_old->Draw("");
  	e_ratio_data->Draw("same");
  	e_ratio_mc->Draw("same");
	if(fitting) fitfunc->Draw("same");
	Legend1->Draw();
	gPad->RedrawAxis();
	col ->Update();
//	Util::PrintNoEPS(col, outputfile, outputdir2, false);
//	Util::PrintEPS(col, outputfile, outputdir2);

}