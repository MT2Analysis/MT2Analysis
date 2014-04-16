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
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"



using namespace std;


void load(const char* filename="samples/samples_20110606_n90.dat");
void run_draw_ratio();
void draw_ratio(TString lower_cut, TString upper_cut, TString var, TString variablelabel, TString outputdir, TString outputfile, TString basecut, TString trigger, TString samples, const int nbins, const double *bins, bool fitting, double lowbord, double upbord);

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
const int fVerbose = 3; // already defined below
TString fPath;

const int gNbins_lowHT_lowMT2                   = 13;
const double  gbins_lowHT_lowMT2[gNbins_lowHT_lowMT2+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 100, 125, 150, 200, 300, 500};

const int gNbins_highHT_lowMT2                   = 13;
const double  gbins_highHT_lowMT2[gNbins_highHT_lowMT2+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 100, 125, 150, 180, 260, 500};


const int gNbins_highMT2                   = 15;
const double  gbins_highMT2[gNbins_highMT2+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 100, 125, 150, 200, 275, 375, 500, 700};

//this function combines the functionality of save_lowercut_uppercut_histos.C and fit_lowercut_uppercut_ratios.C
//for comments please see this other two functions - this function is left without much comment
//note that this is a 7 TeV macro - which is not used for 8 TeV
void run_draw_ratio(){

	gROOT->ProcessLine(".x SetStyle_PRD.C");
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);

	bool lowMT2  = true;
	bool highMT2 = false;
	
	bool lowHT   = false;
	bool highHT  = false;
	bool allHT   = true;

	bool ele     = false;
	bool muo     = false;

	TString samples;

     if(highMT2 && highHT)
	samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4400_noMinDPhi_MT2_highHT_2.dat";
     if(highMT2 && (lowHT || allHT))
	samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4400_noMinDPhi_MT2_2.dat";
     if(lowMT2 && (lowHT || allHT))
	samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4400_noMinDPhi_MT2bSF_2.dat";
     if(lowMT2 && highHT)
	samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4400_noMinDPhi_MT2bSF_highHT_2.dat";

	std::ostringstream triggerStream;
	triggerStream << "( "
	<< "(trigger.HLT_HT440_v2 ==1 && misc.Run<161216)" << "||"
	<< "(trigger.HLT_HT450_v2 ==1 && (misc.Run>=161216 && misc.Run< 163269))" << "||"
	<< "(trigger.HLT_HT500_v3 ==1 && (misc.Run>=163269 && misc.Run<=163869))" << "||"
	<< "(trigger.HLT_HT500_v4 ==1 && (misc.Run>=165088 && misc.Run< 165970))" << "||"
	<< "(trigger.HLT_HT550_v5 ==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" << "||"
	<< "(trigger.HLT_HT550_v6 ==1 && (misc.Run==166346))" << "||"
	<< "(trigger.HLT_HT550_v7 ==1 && (misc.Run>=167078 && misc.Run< 170249))" << "||"
	<< "(trigger.HLT_HT550_v8 ==1 && (misc.Run>=170249 && misc.Run< 173236))" << "||"
	<< "(trigger.HLT_HT600_v1 ==1 && (misc.Run>=173236 && misc.Run< 178420))" << "||"
	<< "(trigger.HLT_HT650_v4 ==1 && (misc.Run>=178420 && misc.Run< 179959))" << " )";//lowHT
//	<< "(trigger.HLT_HT700_v2 ==1 && (misc.Run>=178420 && misc.Run< 179959))" << "||"//highHT
//	<< "(trigger.HLT_PFHT650_v1==1 && misc.Run>=179959)" << " )";//highHT
	TString trigger = triggerStream.str().c_str();

	std::ostringstream cutStream;
	cutStream << " " 
    << "misc.MET>=30"                                                  << "&&"
    << "misc.Jet0Pass==1"                                              << "&&"
    << "misc.Jet1Pass==1"                                              << "&&"
    << "misc.SecondJPt  >100"                                          << "&&"
    << "misc.PassJetID ==1"                                            << "&&"
    << "misc.Vectorsumpt < 70"                                         << "&&"
    << "(misc.ProcessID!=6 || (misc.Event!=814918 && misc.Event!=5500089 && misc.Event!=1934425))" << "&&"
    // Noise
    << "misc.HBHENoiseFlagIso==0"                                      << "&&"
    << "misc.CSCTightHaloID==0"                                        << "&&"
    << "misc.CrazyHCAL==0";
    if(ele){
        cutStream << "&&"
	<< "(NMuons==0 || muo[0].lv.Pt()<10)"                        << "&&"
	<< "NEles>0 && ele[0].lv.Pt()>10&&ele[1].lv.Pt()<10";/*	     << "&&"
	<< "ele[0].MT<100";*/
    }
    else if(muo){
        cutStream << "&&"
	<< "(NEles==0  || ele[0].lv.Pt()<10)"                        << "&&"
	<< "NMuons>0 && muo[0].lv.Pt()>10&&muo[1].lv.Pt()<10";/*	     << "&&"
	<< "muo[0].MT<100";*/
    }
    else{
        cutStream << "&&"
	<< "(NEles==0  || ele[0].lv.Pt()<10)"                        << "&&"
	<< "(NMuons==0 || muo[0].lv.Pt()<10)";
    }
    if(lowMT2){
        cutStream << "&&"
        << "misc.LeadingJPt  >150"                                     << "&&"
        << "NJetsIDLoose40 >=4"                                        << "&&"
        << "NBJets >0"                                                 << "&&"
        << "Sum$(jet.lv.Pt()>20&&jet.lv.Pt()<500&&jet.bTagProbSSVHP>2&&abs(jet.lv.Eta())<2.4)>=1";
    }
    if(highMT2){
        cutStream << "&&"
        << "NJetsIDLoose40 >=3";
    }
    if(highHT){
		cutStream << "&&" << "misc.HT>=950";
    }
    if(lowHT){
		cutStream << "&&" << "misc.HT<950&&misc.HT>750";
    }
    if(allHT){
		cutStream << "&&" << "misc.HT>750";
    }
  	TString cuts = cutStream.str().c_str();


  	TString upper_cut = "NJetsIDLoose40>=0";//dummy
	TString lower_cut;
	if(highMT2) lower_cut = "misc.MinMetJetDPhi >0.3";
	if(lowMT2)  lower_cut = "misc.MinMetJetDPhi4 >0.3";

	TString outputdir = "LostLepton/newfits";
	TString outputfile = "MinDPhiRatio_";
	outputdir  = "kk/dummy";
	outputfile = "MTcut_";
	if(ele) outputfile = outputfile + (TString)"ele_";
	if(muo) outputfile = outputfile + (TString)"muo_";
	if(lowMT2)  outputfile = outputfile + (TString)"MT2b_";
	if(highMT2) outputfile = outputfile + (TString)"MT2_";
	if(lowHT)  outputfile = outputfile + (TString)"lowHT";
	if(highHT) outputfile = outputfile + (TString)"highHT";
	if(allHT)  outputfile = outputfile + (TString)"allHT";

	TString variable = "misc.MT2";
	TString variablelabel = "MT2";

	double lowerfitborder = 125.;
	if(lowMT2)  lowerfitborder = 125.;
	if(highMT2) lowerfitborder = 150.;
	double upperfitborder = 1000.;
	if(lowMT2 && highHT) upperfitborder = gbins_highHT_lowMT2[gNbins_highHT_lowMT2];//high HT, MT2b
	if(lowMT2 && (lowHT||allHT)) upperfitborder = gbins_lowHT_lowMT2[gNbins_lowHT_lowMT2];//low HT, MT2b
	if(highMT2) upperfitborder = gbins_highMT2[gNbins_highMT2];//both HT, MT2

	if(highMT2) draw_ratio(lower_cut, upper_cut, variable, variablelabel, outputdir, outputfile, cuts, trigger, samples, gNbins_highMT2, gbins_highMT2, true, lowerfitborder, upperfitborder); //MT2
	if(lowMT2 && highHT) draw_ratio(lower_cut, upper_cut, variable, variablelabel, outputdir, outputfile, cuts, trigger, samples, gNbins_highHT_lowMT2, gbins_highHT_lowMT2, true, lowerfitborder, upperfitborder); //MT2b high HT
	if(lowMT2 && (lowHT||allHT)) draw_ratio(lower_cut, upper_cut, variable, variablelabel, outputdir, outputfile, cuts, trigger, samples, gNbins_lowHT_lowMT2, gbins_lowHT_lowMT2, true, lowerfitborder, upperfitborder); //MT2b low HT

}

void draw_ratio(TString lower_cut, TString upper_cut, TString var, TString variablelabel, TString outputdir, TString outputfile, TString basecut, TString trigger, TString samples, const int nbins, const double *bins, bool fitting, double lowbord, double upbord){


  load(samples.Data());

  //histos
  TH1D *h_axis        = new TH1D("h_axis",        outputfile,        nbins, bins);
  h_axis   ->SetXTitle(variablelabel); h_axis   ->SetYTitle("ratio");

  TH1D *h_lower_band_mc   = new TH1D("h_lower_band_mc",   outputfile,   nbins, bins);
  h_lower_band_mc         ->SetLineColor(1);               h_lower_band_mc     ->Sumw2();
  h_lower_band_mc 	  ->SetMarkerStyle(22);
  TH1D *h_upper_band_mc   = new TH1D("h_upper_band_mc",   outputfile,   nbins, bins);
  h_upper_band_mc         ->SetLineColor(2);               h_upper_band_mc     ->Sumw2();
  h_upper_band_mc 	  ->SetMarkerStyle(23);
  TH1D *h_ratio_mc        = new TH1D("h_ratio_mc",        outputfile,        nbins, bins);
  h_ratio_mc              ->SetLineColor(1);               h_ratio_mc          ->Sumw2();
  h_ratio_mc 	  	  ->SetMarkerStyle(23);
  TGraphAsymmErrors* g_ratio_mc = new TGraphAsymmErrors();
  g_ratio_mc              ->SetLineColor(1);
  g_ratio_mc 	  	  ->SetMarkerStyle(23);

  TH1D *h_lower_band_susy = new TH1D("h_lower_band_susy", outputfile, nbins, bins);
  h_lower_band_susy       ->SetLineColor(1);               h_lower_band_susy   ->Sumw2();
  TH1D *h_upper_band_susy = new TH1D("h_upper_band_susy", outputfile, nbins, bins);
  h_upper_band_susy       ->SetLineColor(2);               h_upper_band_susy   ->Sumw2();
  TH1D *h_ratio_susy      = new TH1D("h_ratio_susy",      outputfile,      nbins, bins);
  h_ratio_susy            ->SetLineColor(4);               h_ratio_susy        ->Sumw2();
  TGraphAsymmErrors* g_ratio_susy = new TGraphAsymmErrors();
  g_ratio_susy              ->SetLineColor(4);

  TH1D *h_lower_band_data = new TH1D("h_lower_band_data", outputfile, nbins, bins);
  h_lower_band_data       ->SetLineColor(1);               h_lower_band_data   ->Sumw2();
  h_lower_band_data   	  ->SetMarkerStyle(22);
  TH1D *h_upper_band_data = new TH1D("h_upper_band_data", outputfile, nbins, bins);
  h_upper_band_data       ->SetLineColor(2);               h_upper_band_data   ->Sumw2();
  h_upper_band_data   	  ->SetMarkerStyle(23);
  TH1D *h_ratio_data      = new TH1D("h_ratio_data",      outputfile,      nbins, bins);
  h_ratio_data            ->SetLineColor(2);               h_ratio_data        ->Sumw2();
  h_ratio_data 	  	  ->SetMarkerStyle(22);
  TGraphAsymmErrors* g_ratio_data = new TGraphAsymmErrors();
  g_ratio_data              ->SetLineColor(2);
  g_ratio_data 	  	  ->SetMarkerStyle(22);

  h_lower_band_mc   ->SetXTitle(variablelabel); h_lower_band_mc   ->SetYTitle("events");
  h_upper_band_mc   ->SetXTitle(variablelabel); h_upper_band_mc   ->SetYTitle("events");
  h_lower_band_susy ->SetXTitle(variablelabel); h_lower_band_susy ->SetYTitle("events");
  h_upper_band_susy ->SetXTitle(variablelabel); h_upper_band_susy ->SetYTitle("events");
  h_upper_band_data ->SetXTitle(variablelabel); h_upper_band_data ->SetYTitle("events");
  h_lower_band_data ->SetXTitle(variablelabel); h_lower_band_data ->SetYTitle("events");

  h_ratio_mc   ->SetXTitle(variablelabel); h_ratio_mc   ->SetYTitle("ratio");
  h_ratio_susy ->SetXTitle(variablelabel); h_ratio_susy ->SetYTitle("ratio");
  h_ratio_data ->SetXTitle(variablelabel); h_ratio_data ->SetYTitle("ratio");

  g_ratio_mc   ->SetNameTitle("g_ratio_mc", outputfile); g_ratio_mc   ->GetXaxis()->SetTitle(variablelabel); g_ratio_mc   ->GetYaxis()->SetTitle("ratio");
  g_ratio_susy ->SetNameTitle("g_ratio_susy", outputfile); g_ratio_susy ->GetXaxis()->SetTitle(variablelabel); g_ratio_susy ->GetYaxis()->SetTitle("ratio");
  g_ratio_data ->SetNameTitle("g_ratio_data", outputfile); g_ratio_data ->GetXaxis()->SetTitle(variablelabel); g_ratio_data ->GetYaxis()->SetTitle("ratio");

  for(size_t i = 0; i < fSamples.size(); ++i){
    
    if(fSamples[i].sname=="QCD") continue;//for debugging
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    
//    Double_t weight = fSamples[i].weight();//old version
    Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
    cout << "----------------------------------------" << endl;
    cout << "++ looping over " << fSamples[i].name << endl;
    cout << "   sample has weight "   << weight << " and " << fSamples[i].nevents << "(" << fSamples[i].tree->GetEntries() << ") entries" << endl; 

    TString cuts = basecut;
    if(fSamples[i].type=="data" && trigger!="") cuts += " && ("+trigger+")";

    TString weights   = (fSamples[i].type!="data" ?  TString::Format("(%.15f*pileUp.Weight)",weight) : TString::Format("(%.15f)",weight));
    TString selection = TString::Format("(%s) * (%s)"      ,weights.Data(),cuts.Data());
    TString sel_up    = TString::Format("(%s) * (%s && %s)",weights.Data(),cuts.Data(),upper_cut.Data());
    TString sel_lo    = TString::Format("(%s) * (%s && %s)",weights.Data(),cuts.Data(),lower_cut.Data());
      
    TString variable;

    if (fSamples[i].type == "mc"){
      variable  = TString::Format("%s>>+%s",var.Data(),h_lower_band_mc->GetName());
      cout << " drawing: " << variable << ", sel: "<< sel_lo << endl;
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      cout << "Events found: " << h_lower_band_mc->Integral() << endl;
      variable  = TString::Format("%s>>+%s",var.Data(),h_upper_band_mc->GetName());
      cout << " drawing: " << variable << ", sel: "<< sel_up << endl;
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
      cout << "Events found: " << h_upper_band_mc->Integral() << endl;
    }
    else if (fSamples[i].type == "susy"){
      variable  = TString::Format("%s>>+%s",var.Data(),h_lower_band_susy->GetName());
      cout << " drawing: " << variable << ", sel: "<< sel_lo << endl;
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      cout << "Events found: " << h_lower_band_data->Integral() << endl;
      variable  = TString::Format("%s>>+%s",var.Data(),h_upper_band_susy->GetName());
      cout << " drawing: " << variable << ", sel: "<< sel_up << endl;
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
      cout << "Events found: " << h_upper_band_susy->Integral() << endl;
    }
    else if (fSamples[i].type == "data"){
      variable  = TString::Format("%s>>+%s",var.Data(),h_lower_band_data->GetName());
      cout << " drawing: " << variable << ", sel: "<< sel_lo << endl;
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      cout << "Events found: " << h_lower_band_data->Integral() << endl;
      variable  = TString::Format("%s>>+%s",var.Data(),h_upper_band_data->GetName());
      cout << " drawing: " << variable << ", sel: "<< sel_up << endl;
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
      cout << "Events found: " << h_upper_band_data->Integral() << endl;
    }

  }//for size_t =i


  h_ratio_mc   ->Divide(h_lower_band_mc  , h_upper_band_mc,   1, 1, "B");
  h_ratio_data ->Divide(h_lower_band_data, h_upper_band_data, 1, 1, "B");
  h_ratio_susy ->Divide(h_lower_band_susy, h_upper_band_susy, 1, 1, "B");

  TH1D *h_low_mc;
  TH1D *h_up_mc;
  TH1D *h_low_susy;
  TH1D *h_up_susy;
  TH1D *h_low_data;
  TH1D *h_up_data;


  h_low_mc   = (TH1D*)h_lower_band_mc  ->Clone("h_low_mc"  );
  h_up_mc    = (TH1D*)h_upper_band_mc  ->Clone("h_up_mc"   );
  h_low_susy = (TH1D*)h_lower_band_susy->Clone("h_low_susy");
  h_up_susy  = (TH1D*)h_upper_band_susy->Clone("h_up_susy" );
  h_low_data = (TH1D*)h_lower_band_data->Clone("h_low_data");
  h_up_data  = (TH1D*)h_upper_band_data->Clone("h_up_data" );

  TEfficiency* e_ratio_mc   = new TEfficiency((*h_low_mc),   (*h_up_mc)  );
  TEfficiency* e_ratio_susy = new TEfficiency((*h_low_susy), (*h_up_susy));
  TEfficiency* e_ratio_data = new TEfficiency((*h_low_data), (*h_up_data));
  TEfficiency* e_fit        = new TEfficiency((*h_low_mc),   (*h_up_mc)  );

  e_ratio_mc  ->SetNameTitle("e_ratio_mc",   outputfile);
  e_ratio_data->SetNameTitle("e_ratio_data", outputfile);
  e_ratio_susy->SetNameTitle("e_ratio_susy", outputfile);
  e_ratio_mc              ->SetLineColor(1);
  e_ratio_mc 	  	  ->SetMarkerStyle(23);
  e_ratio_susy            ->SetLineColor(4);
  e_ratio_data            ->SetLineColor(2);
  e_ratio_data 	  	  ->SetMarkerStyle(22);

  g_ratio_mc   ->BayesDivide(h_lower_band_mc  , h_upper_band_mc  );
  g_ratio_data ->BayesDivide(h_lower_band_data, h_upper_band_data);
  g_ratio_susy ->BayesDivide(h_lower_band_susy, h_upper_band_susy);

  TF1 *fitfunc = new TF1("fit", "[0]", lowbord, upbord);//defined 110
  TF1 *fitfunc1 = new TF1("fit1", "[0] + x*[1]", lowbord, upbord);//defined 110
  TF1 *fitfunc2 = new TF1("fit2", "[0]*(1 - exp(-(x*[1])))", lowbord, upbord);//defined 110
  TF1 *fitfunc3 = new TF1("fit3", "[0]*(1 - exp(-(x*[1])))", 0, upbord);//defined 110

  if(fitting){
	double mean;
	int lbin = h_ratio_mc->FindBin(lowbord+0.01);//+0.01 to exclude border effects
	int ubin = h_ratio_mc->FindBin(upbord-0.01);
	mean = (h_ratio_mc->Integral(lbin, ubin))/(ubin-lbin+1);//mean value of ratio for MT2>100
	fitfunc->SetParameter(0, mean);
	fitfunc->SetLineColor(6);
	fitfunc->SetLineWidth(2);
	fitfunc1->SetParameter(0, mean);
	fitfunc2->SetParameter(0, mean);
	fitfunc3->SetParameter(0, mean);
	fitfunc1->SetParameter(1, 0);
	fitfunc2->SetParameter(1, 0);
	fitfunc3->SetParameter(1, 0);
	fitfunc1->SetLineColor(7);
	fitfunc1->SetLineWidth(1);
	fitfunc2->SetLineColor(8);
	fitfunc2->SetLineWidth(1);
	fitfunc3->SetLineColor(9);
	fitfunc3->SetLineWidth(1);
	cout << endl << "fit horizontal line" << endl;
	h_ratio_mc->Fit(fitfunc, "R");
	cout << endl << "fit slopped line" << endl;
	h_ratio_mc->Fit(fitfunc1, "R0+");
	cout << endl << "fit exponential" << endl;
	h_ratio_mc->Fit(fitfunc2, "R0+");
	cout << endl << "fit full-range exponential line" << endl;
	h_ratio_mc->Fit(fitfunc3, "R0+");
	cout << "fitting done" << endl;
	cout << "Fit " << outputfile << " in range MT2(" << lowbord << "," << upbord << ")" << " (for MC)" << endl;
	cout << "Fit is a straight line." << endl;
	cout << "Fit constant = " << fitfunc->GetParameter(0) << " +/- " << fitfunc->GetParError(0) << endl;
	cout << "Goodness of fit: chi^2/NDF = " << fitfunc->GetChisquare() << "/" << fitfunc->GetNDF() << " = " << (fitfunc->GetChisquare())/(fitfunc->GetNDF()) << endl;
	cout << "Fit is a tilted line." << endl;
	cout << "Fit constant = " << fitfunc1->GetParameter(0) << " +/- " << fitfunc1->GetParError(0) << endl;
	cout << "Fit slope    = " << fitfunc1->GetParameter(1) << " +/- " << fitfunc1->GetParError(1) << endl;
	cout << "Goodness of fit: chi^2/NDF = " << fitfunc1->GetChisquare() << "/" << fitfunc1->GetNDF() << " = " << (fitfunc1->GetChisquare())/(fitfunc1->GetNDF()) << endl;
	cout << "Fit is a exponential." << endl;
	cout << "Fit constant = " << fitfunc2->GetParameter(0) << " +/- " << fitfunc2->GetParError(0) << endl;
	cout << "Fit exp.fact = " << fitfunc2->GetParameter(1) << " +/- " << fitfunc2->GetParError(1) << endl;
	cout << "Goodness of fit: chi^2/NDF = " << fitfunc2->GetChisquare() << "/" << fitfunc2->GetNDF() << " = " << (fitfunc2->GetChisquare())/(fitfunc2->GetNDF()) << endl;
	cout << "Fit is a exponential in complete range." << endl;
	cout << "Fit constant = " << fitfunc3->GetParameter(0) << " +/- " << fitfunc3->GetParError(0) << endl;
	cout << "Fit exp.fact = " << fitfunc3->GetParameter(1) << " +/- " << fitfunc3->GetParError(1) << endl;
	cout << "Goodness of fit: chi^2/NDF = " << fitfunc3->GetChisquare() << "/" << fitfunc3->GetNDF() << " = " << (fitfunc3->GetChisquare())/(fitfunc3->GetNDF()) << endl;

  }

  TCanvas *can1 = new TCanvas("can1", "mc", 0, 0, 900, 700);
  can1->SetLogy(1);
  h_lower_band_mc->Draw();  
  h_upper_band_mc->Draw("same");

  TCanvas *can2 = new TCanvas("can2", "ratio mc", 0, 0, 900, 700);
  //can2->SetLogy(1);
  h_ratio_mc->Draw("E1");

  TCanvas *can3 = new TCanvas("can3", "susy", 0, 0, 900, 700);
  can3->SetLogy(1);
  h_lower_band_susy->Draw();  
  h_upper_band_susy->Draw("same");

  TCanvas *can4 = new TCanvas("can4", "ratio susy", 0, 0, 900, 700);
  //can4->SetLogy(1);
  h_ratio_susy->Draw("E1");

  TCanvas *can5 = new TCanvas("can5", "data", 0, 0, 900, 700);
  can5->SetLogy(1);
  h_lower_band_data->Draw();  
  h_upper_band_data->Draw("same");

  TCanvas *can6 = new TCanvas("can6", "ratio data", 0, 0, 900, 700);
  //can4->SetLogy(1);
  h_ratio_data->Draw("E1");

  Util::MakeOutputDir(outputdir);
  TFile *ratio_file = new TFile(outputdir + "/" +  outputfile+TString(".root"), "RECREATE");
  ratio_file		->cd();
  h_lower_band_mc	->Write();
  h_upper_band_mc	->Write();
  h_ratio_mc		->Write();
  h_lower_band_susy	->Write();
  h_upper_band_susy	->Write();
  h_ratio_susy		->Write();
  h_lower_band_data	->Write();
  h_upper_band_data	->Write();
  h_ratio_data		->Write();
  g_ratio_data		->Write();
  g_ratio_susy		->Write();
  g_ratio_mc		->Write();
  e_ratio_data		->Write();
  e_ratio_susy		->Write();
  e_ratio_mc		->Write();
  if(fitting) fitfunc	->Write();
  if(fitting) fitfunc1	->Write();
  if(fitting) fitfunc2	->Write();
  if(fitting) fitfunc3	->Write();
  cout << "saved in " << outputfile << endl;

	double min = 0;
	double max = 1.5;

	h_ratio_data  ->SetMaximum(max);
	h_ratio_mc->SetMaximum(max);
	h_ratio_data  ->SetMinimum(min);
	h_ratio_mc->SetMinimum(min);
	h_axis->SetMaximum(max);
	h_axis->SetMinimum(min);

	TLegend* Legend2 = new TLegend(.55,.6,.85,.85);
	Legend2     ->AddEntry(h_upper_band_mc, "MC relaxed cuts", "lp");
	Legend2     ->AddEntry(h_lower_band_mc, "MC full cuts",    "lp");
	Legend2     -> SetFillColor(0);
	Legend2     -> SetBorderSize(0);

  	TCanvas *col2 = new TCanvas("col2", outputfile+TString("_mc"), 0, 0, 900, 700);
	gPad->SetLogy(1);
  	h_upper_band_mc->Draw("E1");
  	h_lower_band_mc->Draw("E1same");
	Legend2->Draw();
	gPad->RedrawAxis();
	col2 ->Update();
	Util::PrintNoEPS(col2, outputfile+TString("_mc"), outputdir+TString("normalized/"), false);
	Util::PrintEPS(col2, outputfile+TString("_mc"), outputdir+TString("normalized/"));

	TLegend* Legend3 = new TLegend(.55,.6,.85,.85);
	Legend3     ->AddEntry(h_upper_band_data, "Data relaxed cuts", "lp");
	Legend3     ->AddEntry(h_lower_band_data, "Data full cuts",    "lp");
	Legend3     -> SetFillColor(0);
	Legend3     -> SetBorderSize(0);

  	TCanvas *col3 = new TCanvas("col3", outputfile+TString("_data"), 0, 0, 900, 700);
	gPad->SetLogy(1);
  	h_upper_band_data->Draw("E1");
  	h_lower_band_data->Draw("E1same");
	Legend3->Draw();
	gPad->RedrawAxis();
	col3 ->Update();
	Util::PrintNoEPS(col3, outputfile+TString("_data"), outputdir+TString("normalized/"), false);
	Util::PrintEPS(col3, outputfile+TString("_data"), outputdir+TString("normalized/"));

	TLegend* Legend1 = new TLegend(.6,.7,.88,.88);
	Legend1     ->AddEntry(h_ratio_mc,   "MT cut efficiency",   "lp");
//	Legend1     ->AddEntry(h_ratio_data, "ratio data", "lp");
	Legend1     -> SetFillColor(0);
	Legend1     -> SetBorderSize(0);

  	TCanvas *col = new TCanvas("col", outputfile, 0, 0, 900, 700);
	gPad->SetLogy(0);
	h_ratio_mc->Draw("AXIS");
  	h_ratio_mc->Draw("same");
  //	h_ratio_data->Draw("same");
	fitfunc->Draw("same");
	Legend1->Draw();
	gPad->RedrawAxis();
	col ->Update();
	Util::PrintNoEPS(col, outputfile + TString("new"), outputdir, false);
	Util::PrintEPS(col, outputfile + TString("new"), outputdir);

}

//reads in samples.dat
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