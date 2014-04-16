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
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path of MT2tree.hh
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"

//note that this is a 7 TeV code
//run via root -l -b -q save_lowercut_uppercut_histos.C++

using namespace std;


void load(const char* filename="samples/samples_20110606_n90.dat");
void save_lowercut_uppercut_histos();
void draw_ratio(TString lower_cut, TString upper_cut, TString title, TString var, TString variablelabel, TString outputdir, TString outputfile, TString basecut, TString trigger, TString samples, const int nbins, const double *bins, bool fitting, double lowbord, double upbord);

//samples combining MT2trees with necessary information like cross section
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


	bool lowMT2  = false;//MT2b (7 TeV)
	bool highMT2 = true; //MT2 (7 TeV)
	
	bool lowHT   = false;//low HT (7 TeV)
	bool highHT  = false;//high HT (7 TeV)
	bool allHT   = true; //low+high HT

	bool ele     = false;//1 electron selection, else hadronic or 1 muon selection
	bool muo     = true; //1 muon selection, else hadronic or 1 electron selection

//this cut saves ratios (for data, SM simulation, and SUSY signal) for a loose cut and a tight cut (can be overlapping cuts or not,
//e.g. tight >=1b, loose 0b, but also loose >=0b possible).
//note that this is a 7 TeV code
//it was used for the 1 fb-1 result (PAS SUS-11-005) for Multijet background prediction for the 'low MT2' Analysis
//later this macro was combined with fit_lowercut_uppercut_ratios.C ==> see run_draw_ratio.C
void save_lowercut_uppercut_histos(){

	gROOT->ProcessLine(".x SetStyle_PRD.C");
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);

	//get the samples
	TString samples;
     if(highMT2 && highHT)
	samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4400_noMinDPhi_MT2_highHT_3.dat";
     if(highMT2 && (lowHT || allHT))
	samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4400_noMinDPhi_MT2_3.dat";
     if(lowMT2 && (lowHT || allHT))
	samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4400_noMinDPhi_MT2bSF_2.dat";
     if(lowMT2 && highHT)
	samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4400_noMinDPhi_MT2bSF_highHT_2.dat";

	//get the (loose) event selection
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
	<< "NEles>0 && ele[0].lv.Pt()>10&&ele[1].lv.Pt()<10"	     << "&&"
	<< "ele[0].MT<100";
    }
    else if(muo){
        cutStream << "&&"
	<< "(NEles==0  || ele[0].lv.Pt()<10)"                        << "&&"
	<< "NMuons>0 && muo[0].lv.Pt()>10&&muo[1].lv.Pt()<10"	     << "&&"
	<< "muo[0].MT<100";
    }
    else{
        cutStream << "&&"
	<< "(NEles==0  || ele[0].lv.Pt()<10)"                        << "&&"
	<< "(NMuons==0 || muo[0].lv.Pt()<10)";
    }
    if(lowMT2){
        cutStream << "&&"
        << "misc.LeadingJPt  >150"                                     << "&&"
        << "( (NJetsIDLoose40 >=4) || (NJetsIDLoose40 >=3 && ele[0].lv.Pt()>40) || (NJetsIDLoose40 >=3 && muo[0].lv.Pt()>40) )" << "&&"
        << "NBJets >0";//                                                 << "&&"
    }
    if(highMT2){
        cutStream << "&&"
        << "( (NJetsIDLoose40 >=3) || (NJetsIDLoose40 >=2 && ele[0].lv.Pt()>40) || (NJetsIDLoose40 >=2 && muo[0].lv.Pt()>40) )";
    }
    if(highHT){
	if(ele)
		cutStream << "&&" << "(((misc.HT + ele[0].lv.Pt())>=950&&ele[0].lv.Pt()>50&&misc.HT>400)||(misc.HT>=950&&ele[0].lv.Pt()<=50))";
	else if(muo)
		cutStream << "&&" << "(((misc.HT + muo[0].lv.Pt())>=950&&muo[0].lv.Pt()>50&&misc.HT>400)||(misc.HT>=950&&muo[0].lv.Pt()<=50))";
	else
		cutStream << "&&" << "misc.HT>=950";
    }
    if(lowHT){
	if(ele)
		cutStream << "&&" << "(((misc.HT + ele[0].lv.Pt())>=750 && (misc.HT + ele[0].lv.Pt())<950&&ele[0].lv.Pt()>50&&misc.HT>400) || (misc.HT>=750&&misc.HT<950&&ele[0].lv.Pt()<=50))";
	else if(muo)
		cutStream << "&&" << "(((misc.HT + muo[0].lv.Pt())>=750 && (misc.HT + muo[0].lv.Pt())<950&&muo[0].lv.Pt()>50&&misc.HT>400) || (misc.HT>=750&&misc.HT<950&&muo[0].lv.Pt()<=50))";
	else
		cutStream << "&&" << "misc.HT<950&&misc.HT>750";
    }
    if(allHT){
	if(ele)
		cutStream << "&&" << "(((misc.HT + ele[0].lv.Pt())>=750&&ele[0].lv.Pt()>50&&misc.HT>400)||(misc.HT>=750&&ele[0].lv.Pt()<=50))";
	else if(muo)
		cutStream << "&&" << "(((misc.HT + muo[0].lv.Pt())>=750&&muo[0].lv.Pt()>50&&misc.HT>400)||(misc.HT>=750&&muo[0].lv.Pt()<=50))";
	else
		cutStream << "&&" << "misc.HT>750";
    }
  	TString cuts = cutStream.str().c_str();


	//get loose/tight selection
  	TString upper_cut = "NJetsIDLoose40>=0";//dummy
	TString lower_cut;//veto lowercut in cutstream
	if(highMT2) lower_cut = "misc.MinMetJetDPhi >0.3";
	if(lowMT2)  lower_cut = "misc.MinMetJetDPhi4 >0.3";

	//get titles and names
	TString outputdir = "LostLepton/fits";//newfits if lepton-corrected quantities, fits if normal quantities
	TString outputfile = "MinDPhiRatio_";
	if(ele) outputfile = outputfile + (TString)"ele_";
	if(muo) outputfile = outputfile + (TString)"muo_";
	if(lowMT2)  outputfile = outputfile + (TString)"MT2b_";
	if(highMT2) outputfile = outputfile + (TString)"MT2_";
	if(lowHT)  outputfile = outputfile + (TString)"lowHT";
	if(highHT) outputfile = outputfile + (TString)"highHT";
	if(allHT)  outputfile = outputfile + (TString)"allHT";

	TString title = "MinMetJetDPhi ratio";
	if(lowHT)  title = title + (TString)", low HT";
	if(highHT) title = title + (TString)", high HT";
	if(ele)    title = title + (TString)", 1 Electron";
	if(muo)    title = title + (TString)", 1 Muon";


	TString variable = "misc.MT2";
	TString variablelabel = "MT2";

	double lowerfitborder = 125.;
	if(lowMT2)  lowerfitborder = 125.;
	if(highMT2) lowerfitborder = 150.;
	double upperfitborder = 1000.;
	if(lowMT2 && highHT) upperfitborder = gbins_highHT_lowMT2[gNbins_highHT_lowMT2];//high HT, MT2b
	if(lowMT2 && (lowHT||allHT)) upperfitborder = gbins_lowHT_lowMT2[gNbins_lowHT_lowMT2];//low HT, MT2b
	if(highMT2) upperfitborder = gbins_highMT2[gNbins_highMT2];//both HT, MT2
	//draw ratios
	if(highMT2) draw_ratio(lower_cut, upper_cut, title, variable, variablelabel, outputdir, outputfile, cuts, trigger, samples, gNbins_highMT2, gbins_highMT2, true, lowerfitborder, upperfitborder); //MT2
	if(lowMT2 && highHT) draw_ratio(lower_cut, upper_cut, title, variable, variablelabel, outputdir, outputfile, cuts, trigger, samples, gNbins_highHT_lowMT2, gbins_highHT_lowMT2, true, lowerfitborder, upperfitborder); //MT2b high HT
	if(lowMT2 && (lowHT||allHT)) draw_ratio(lower_cut, upper_cut, title, variable, variablelabel, outputdir, outputfile, cuts, trigger, samples, gNbins_lowHT_lowMT2, gbins_lowHT_lowMT2, true, lowerfitborder, upperfitborder); //MT2b low HT

}

void draw_ratio(TString lower_cut, TString upper_cut, TString title, TString var, TString variablelabel, TString outputdir, TString outputfile, TString basecut, TString trigger, TString samples, const int nbins, const double *bins, bool fitting, double lowbord, double upbord){


  load(samples.Data());

  //histos
  TH1D *h_axis        = new TH1D("h_axis",        title,        nbins, bins);
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

  TH1D *h_lower_band_susy = new TH1D("h_lower_band_susy", outputfile, nbins, bins);
  h_lower_band_susy       ->SetLineColor(1);               h_lower_band_susy   ->Sumw2();
  TH1D *h_upper_band_susy = new TH1D("h_upper_band_susy", outputfile, nbins, bins);
  h_upper_band_susy       ->SetLineColor(2);               h_upper_band_susy   ->Sumw2();
  TH1D *h_ratio_susy      = new TH1D("h_ratio_susy",      outputfile,      nbins, bins);
  h_ratio_susy            ->SetLineColor(4);               h_ratio_susy        ->Sumw2();
  h_ratio_susy	  	  ->SetMarkerStyle(23);

  TH1D *h_lower_band_data = new TH1D("h_lower_band_data", outputfile, nbins, bins);
  h_lower_band_data       ->SetLineColor(1);               h_lower_band_data   ->Sumw2();
  h_lower_band_data   	  ->SetMarkerStyle(22);
  TH1D *h_upper_band_data = new TH1D("h_upper_band_data", outputfile, nbins, bins);
  h_upper_band_data       ->SetLineColor(2);               h_upper_band_data   ->Sumw2();
  h_upper_band_data   	  ->SetMarkerStyle(23);
  TH1D *h_ratio_data      = new TH1D("h_ratio_data",      outputfile,      nbins, bins);
  h_ratio_data            ->SetLineColor(2);               h_ratio_data        ->Sumw2();
  h_ratio_data 	  	  ->SetMarkerStyle(22);

  h_lower_band_mc   ->SetXTitle(variablelabel); h_lower_band_mc   ->SetYTitle("events");
  h_upper_band_mc   ->SetXTitle(variablelabel); h_upper_band_mc   ->SetYTitle("events");
  h_lower_band_susy ->SetXTitle(variablelabel); h_lower_band_susy ->SetYTitle("events");
  h_upper_band_susy ->SetXTitle(variablelabel); h_upper_band_susy ->SetYTitle("events");
  h_upper_band_data ->SetXTitle(variablelabel); h_upper_band_data ->SetYTitle("events");
  h_lower_band_data ->SetXTitle(variablelabel); h_lower_band_data ->SetYTitle("events");

  h_ratio_mc   ->SetXTitle(variablelabel); h_ratio_mc   ->SetYTitle("efficiency");
  h_ratio_susy ->SetXTitle(variablelabel); h_ratio_susy ->SetYTitle("efficiency");
  h_ratio_data ->SetXTitle(variablelabel); h_ratio_data ->SetYTitle("efficiency");

  for(size_t i = 0; i < fSamples.size(); ++i){
    
//    if(fSamples[i].sname=="QCD") continue;//for debugging
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    
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
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");//draw loose/lower selection
      cout << "Events found: " << h_lower_band_mc->Integral() << endl;
      variable  = TString::Format("%s>>+%s",var.Data(),h_upper_band_mc->GetName());
      cout << " drawing: " << variable << ", sel: "<< sel_up << endl;
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");//draw tight/upper selection
      cout << "Events found: " << h_upper_band_mc->Integral() << endl;
    }
    else if (fSamples[i].type == "susy"){
      variable  = TString::Format("%s>>+%s",var.Data(),h_lower_band_susy->GetName());
      cout << " drawing: " << variable << ", sel: "<< sel_lo << endl;
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");//draw loose/lower selection
      cout << "Events found: " << h_lower_band_data->Integral() << endl;
      variable  = TString::Format("%s>>+%s",var.Data(),h_upper_band_susy->GetName());
      cout << " drawing: " << variable << ", sel: "<< sel_up << endl;
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");//draw tight/upper selection
      cout << "Events found: " << h_upper_band_susy->Integral() << endl;
    }
    else if (fSamples[i].type == "data"){
      variable  = TString::Format("%s>>+%s",var.Data(),h_lower_band_data->GetName());
      cout << " drawing: " << variable << ", sel: "<< sel_lo << endl;
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");//draw loose/lower selection
      cout << "Events found: " << h_lower_band_data->Integral() << endl;
      variable  = TString::Format("%s>>+%s",var.Data(),h_upper_band_data->GetName());
      cout << " drawing: " << variable << ", sel: "<< sel_up << endl;
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");//draw tight/upper selection
      cout << "Events found: " << h_upper_band_data->Integral() << endl;
    }

  }//for size_t =i

  //make the ratios
  h_ratio_mc   ->Divide(h_lower_band_mc  , h_upper_band_mc,   1, 1, "B");//binomial uncertainties as it is fail/pass
  h_ratio_data ->Divide(h_lower_band_data, h_upper_band_data, 1, 1, "B");
  h_ratio_susy ->Divide(h_lower_band_susy, h_upper_band_susy, 1, 1, "");

  Util::MakeOutputDir(outputdir);//create directory, if not existent
  //save the files
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
  h_axis		->Write();
  cout << "saved in " << outputfile << endl;

}

//function to read in samples.dat
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