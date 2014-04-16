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
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path of MT2tree.hh
#include "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"

//run via root -l make2Dplot.C++

using namespace std;

void load(const char* filename = "datasamples/samples_2141_dataonly.dat");

//sample combines MT2trees with necessary information like cross section
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
Int_t   fVerbose                  = 3;
TString fPath;

Bool_t  runMC                     = true; //if false, false for all MC events, default = true (same for all MC samples)
Bool_t  runQCD                    = true; // set to false to save some time
Bool_t  runPhotons                = true;
Bool_t  runZll                    = true;
Bool_t  runZnunu                  = true;
Bool_t  runWJets                  = true;
Bool_t  runTTbar                  = true; // includes TTZ,TTW
Bool_t  runSingleTop              = true;
Bool_t  runData                   = true; //default = true
Bool_t  runSUSY                   = false;//default = false

Bool_t  do3D                      = false; //instead of 2d color plot, 3d plot, does not really work with overlaying to plots, default = false
Bool_t  compareMCDATA             = true; //in plotting compare MC to data, MC in color, data as number, default = true
Bool_t  compareMCSUSY             = false;//in plotting compare MC to SUSY, MC in color, SUSY as number, default = false
Bool_t  compareSUSYDATA           = false;//in plotting compare SUSY to data, SUSY in color, data as number, default = false
//a way to make only MC plot or SUSY/data is by setting flags runMC, runData, runSUSY accordingly
Bool_t  logflag			  = true; //set logflag, usually true is a good choice

Bool_t  fSave                     = true;//save plots

const int nxbins = 15; double xbins[nxbins+1] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15.};
const int nybins =  5; double ybins[nybins+1] = {0., 1., 2., 3., 4., 5.};
TString xtitle = "NJets40";//title of x axis in plot
TString ytitle = "NBJets"; //title of y axis in plot
TString xvar = "NJetsIDLoose40";//the variable name of x axis as in the MT2tree branch/leaf
TString yvar = "NBJets";        //the variable name of y axis as in the MT2tree branch/leaf


TString samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4400_noMinDPhi.dat";
TString fOutputDir = "../MassPlots/2Dplots";//output directory
TFile*  fOutputFile = Util::MakeOutputFile(fOutputDir + "MassPlots.root");

//this function makes a 2d plot for chosen variables
//this function comes from 7 TeV analysis!!!
void make2Dplot(){

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	//event selection
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
	cutStream  
		<< "misc.MT2>=125"					     << "&&"

		<< "misc.MET>=30"                                            << "&&"
//		<< "misc.HT<950&&misc.HT>750"                                << "&&"
//		<< "misc.HT>950 "                                            << "&&"
		<< "misc.HT>750 "                                            << "&&"
	//	<< "NJetsIDLoose40 >=3"                                      << "&&"
	//	<< "NJetsIDLoose40 >=4"                                      << "&&"
	//	<< "NBJets >0"                                               << "&&"
		<< "(NEles==0  || ele[0].lv.Pt()<10)"                        << "&&"
		<< "(NMuons==0 || muo[0].lv.Pt()<10)"                        << "&&"

		<< "misc.Jet0Pass==1"                                        << "&&"
		<< "misc.Jet1Pass==1"                                        << "&&"
		<< "misc.SecondJPt >100"                                     << "&&"
		<< "misc.LeadingJPt >150"                                    << "&&"
		<< "misc.PassJetID ==1"                                      << "&&"
		<< "misc.Vectorsumpt < 70"                                   << "&&"
		<< "(misc.MinMetJetDPhi >0.3||misc.MinMetJetDPhiIndex>3)"    << "&&"
	//	<< "misc.MinMetJetDPhi >0.3"                                 << "&&"
		<< "misc.HBHENoiseFlagIso==0"                                << "&&"
		<< "misc.CSCTightHaloID==0"                                  << "&&"
		<< "misc.CrazyHCAL==0";

  	TString cuts = cutStream.str().c_str();

	//***********************************************************************************************************************
	// The part below is basically a copy from MassPlotter.cc, just with the difference of 2 drawing variables instead of one
	// therefore no detailed comments are done here
	//***********************************************************************************************************************

	TString xvarname = Util::removeFunnyChar(xvar.Data());
	TString yvarname = Util::removeFunnyChar(yvar.Data());

  	TH2D*    h_data      = new TH2D   (xvarname+yvarname+"data"  , "", nxbins, xbins, nybins, ybins );
	TH2D*    h_mc_sum    = new TH2D   (xvarname+yvarname+"mc_sum", "", nxbins, xbins, nybins, ybins );
	TH2D*    h_susy      = new TH2D   (xvarname+yvarname+"susy"  , "", nxbins, xbins, nybins, ybins );	

	// h_data
	h_data -> Sumw2();
	h_data -> SetMarkerSize(1);
	h_data -> SetMarkerStyle(20);
	h_data -> SetMarkerColor(kBlack);
	h_data -> SetLineColor(kBlack);
	h_data -> SetStats(false);

	h_mc_sum -> SetFillStyle(3004);
	h_mc_sum -> SetFillColor(kBlack);
	h_mc_sum -> SetLineColor(0);
	h_mc_sum -> SetStats(0);
	h_mc_sum -> Sumw2();

	h_susy -> Sumw2();
	h_susy -> SetStats(0);

	load(samples.Data());

	vector<TH2D*> h_samples;

	for(size_t i = 0; i < fSamples.size(); ++i){

		h_samples.push_back(new TH2D(xvarname+yvarname+"_"+fSamples[i].name, "", nxbins, xbins, nybins, ybins ));
		h_samples[i] -> Sumw2();
		h_samples[i] -> SetFillColor(0);
		h_samples[i] -> SetLineColor(fSamples[i].color);
		h_samples[i] -> SetMarkerColor(fSamples[i].color);
		h_samples[i] -> SetStats(false);

		Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
		if(fVerbose>2) cout << "MakePlot: looping over " << fSamples[i].sname << "-----------------------------------" <<  endl;
		if(fVerbose>2) cout << "  +++++++ xsection:    "    << fSamples[i].xsection << " k-fact " << fSamples[i].kfact << endl;
		if(fVerbose>2) cout << "  +++++++ tot events:  "  << fSamples[i].nevents  << " avg pu weight " << fSamples[i].PU_avg_weight << endl;
		if(fVerbose>2) cout << "  +++++++ weight:      " << weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	
		TString variable  = TString::Format("%s:%s>>+%s",yvar.Data(), xvar.Data(),h_samples[i]->GetName());
		TString theCuts = cuts;
		if(fSamples[i].type=="data" && trigger!="") theCuts += " &&("+trigger+")"; // triggers for data

		TString selection;
		if(fSamples[i].type!="data") selection      = TString::Format("(%.15f*pileUp.Weight) * (%s)",weight,theCuts.Data());
		else                        selection      = TString::Format("(%.15f) * (%s)"              ,weight,theCuts.Data()); 
		  
		if(fVerbose>2) cout << "  +++++++ Drawing      " << variable  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;

		int nev = fSamples[i].tree->Draw(variable.Data(),selection.Data(),"goff");

		if(fVerbose>2) cout << "  +++++++ MC   events : "  <<  nev << endl;
	
	
		//over and underflow bins
		h_samples[i]->SetBinContent(1,1,
				h_samples[i]->GetBinContent(0,0) + h_samples[i]->GetBinContent(1,1) + h_samples[i]->GetBinContent(0,1) + h_samples[i]->GetBinContent(0,1));
		h_samples[i]->SetBinContent(h_samples[i]->GetNbinsX(),h_samples[i]->GetNbinsY(),
				h_samples[i]->GetBinContent(h_samples[i]->GetNbinsX(),h_samples[i]->GetNbinsY()) + h_samples[i]->GetBinContent(h_samples[i]->GetNbinsX()+1,h_samples[i]->GetNbinsY()+1) + 
				h_samples[i]->GetBinContent(h_samples[i]->GetNbinsX()+1,h_samples[i]->GetNbinsY()) + h_samples[i]->GetBinContent(h_samples[i]->GetNbinsX(),h_samples[i]->GetNbinsY()+1));
		h_samples[i]->SetBinError(1,1,
				sqrt(pow(h_samples[i]->GetBinError(0,0),2) + pow(h_samples[i]->GetBinError(1,1),2) + pow(h_samples[i]->GetBinError(0,1),2) + pow(h_samples[i]->GetBinError(0,1),2)) );
		h_samples[i]->SetBinError(h_samples[i]->GetNbinsX(),h_samples[i]->GetNbinsY(),
				sqrt(pow(h_samples[i]->GetBinError(h_samples[i]->GetNbinsX(),h_samples[i]->GetNbinsY()),2) + pow(h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()+1,h_samples[i]->GetNbinsY()+1),2) + 
					pow(h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()+1,h_samples[i]->GetNbinsY()),2) + pow(h_samples[i]->GetBinError(h_samples[i]->GetNbinsX(),h_samples[i]->GetNbinsY()+1),2)) );
		for(int ix = 2; ix< h_samples[i]->GetNbinsX(); ++ix){
			h_samples[i]->SetBinContent(ix,h_samples[i]->GetNbinsY(),
				h_samples[i]->GetBinContent(ix,h_samples[i]->GetNbinsY()) + h_samples[i]->GetBinContent(ix,h_samples[i]->GetNbinsY()+1));
			h_samples[i]->SetBinContent(ix,1,
				h_samples[i]->GetBinContent(ix,0) + h_samples[i]->GetBinContent(ix,1));
			h_samples[i]->SetBinError(ix,h_samples[i]->GetNbinsY(),
				sqrt(pow(h_samples[i]->GetBinError(ix,h_samples[i]->GetNbinsY()),2) + pow(h_samples[i]->GetBinError(ix,h_samples[i]->GetNbinsY()+1),2)) );
			h_samples[i]->SetBinError(ix,1,
				sqrt(pow(h_samples[i]->GetBinError(ix,0),2) + pow(h_samples[i]->GetBinError(ix,1),2)) );
		}
		for(int iy = 2; iy< h_samples[i]->GetNbinsY(); ++iy){
			h_samples[i]->SetBinContent(h_samples[i]->GetNbinsX(), iy, 
				h_samples[i]->GetBinContent(h_samples[i]->GetNbinsX(), iy) + h_samples[i]->GetBinContent(h_samples[i]->GetNbinsX()+1, iy));
			h_samples[i]->SetBinContent(1,iy,
				h_samples[i]->GetBinContent(0,iy) + h_samples[i]->GetBinContent(1,iy));
			h_samples[i]->SetBinError(h_samples[i]->GetNbinsX(),iy,
				sqrt(pow(h_samples[i]->GetBinError(h_samples[i]->GetNbinsX(),iy),2) + pow(h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()+1,iy),2)) );
			h_samples[i]->SetBinError(1,iy,
				sqrt(pow(h_samples[i]->GetBinError(0,iy),2) + pow(h_samples[i]->GetBinError(1,iy),2)) );
		}
		if(fSamples[i].type=="mc") h_mc_sum->Add(h_samples[i]);
		if(fSamples[i].type=="data") h_data->Add(h_samples[i]);
		if(fSamples[i].type=="susy") h_susy->Add(h_samples[i]);

	}

  	if(fVerbose > 2) {
	           cout << "------------------------------------"                << endl
			<< "TOTAL BG:          " << h_mc_sum->Integral()         << endl
		        << "SUSY:              " << h_susy->Integral()           << endl
	        << "Data:              " << h_data->Integral()           << endl;
	}

	TString title = "Events";
	TString oname = xvar+"_"+yvar+ "_CUT_";
	oname += cuts;

	oname.ReplaceAll("misc.BadEcalTP==0&&misc.BadEcalBE==0&&misc.TrackingFailurePVtx>=0.1&&misc.RecovRecHitFilterFlag==0&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloID==0&&misc.CrazyHCAL==0", "AllFilters");
	oname.ReplaceAll("misc.TrackingFailurePVtx>=0.1&&misc.RecovRecHitFilterFlag==0&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloID==0&&misc.CrazyHCAL==0", "AllButEBTPFilters");
	oname.ReplaceAll("misc.RecovRecHitFilterFlag==0&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloID==0&&misc.CrazyHCAL==0", "MT2FiltersInclRecovHits");
	oname.ReplaceAll("misc.HBHENoiseFlag==0&&misc.CSCTightHaloID==0&&misc.CrazyHCAL==0", "MT2Filters");
	oname.ReplaceAll("misc.BadEcalTP==0&&misc.BadEcalBE==0&&misc.TrackingFailurePVtx>=0.1&&misc.RecovRecHitFilterFlag==0&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloID==0&&misc.CrazyHCAL==0", "AllFilters");
	oname.ReplaceAll("misc.TrackingFailurePVtx>=0.1&&misc.RecovRecHitFilterFlag==0&&misc.HBHENoiseFlagIso==0&&misc.CSCTightHaloID==0&&misc.CrazyHCAL==0", "AllButEBTPFiltersIso");
	oname.ReplaceAll("misc.RecovRecHitFilterFlag==0&&misc.HBHENoiseFlagIso==0&&misc.CSCTightHaloID==0&&misc.CrazyHCAL==0", "MT2FiltersIsoInclRecovHits");
	oname.ReplaceAll("misc.HBHENoiseFlagIso==0&&misc.CSCTightHaloID==0&&misc.CrazyHCAL==0", "MT2FiltersIso");
	oname.ReplaceAll("misc.BadEcalTP==0&&misc.BadEcalBE==0", "ECALFlg");

	oname.ReplaceAll("((NEles==0||ele[0].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10))", "NoLeps");
	oname.ReplaceAll("(((NEles>0&&ele[0].lv.Pt()>=10&&ele[1].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10))||((NMuons>0&&muo[0].lv.Pt()>=10&&muo[1].lv.Pt()<10)&&(NEles==0||ele[0].lv.Pt()<10)))", "1Lep");
	oname.ReplaceAll("((NEles>0&&ele[0].lv.Pt()>=10&&ele[1].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10))", "1Ele");
	oname.ReplaceAll("((NMuons>0&&muo[0].lv.Pt()>=10&&muo[1].lv.Pt()<10)&&(NEles==0||ele[0].lv.Pt()<10))", "1Muo");	
	oname.ReplaceAll("misc.Jet0Pass==1&&misc.Jet1Pass==1", "J1J1Pass");
	oname.ReplaceAll("(misc.ProcessID!=6 || (misc.Event!=814918 && misc.Event!=5500089 && misc.Event!=1934425))", "ProcessIDcut");
	oname.ReplaceAll("(misc.MinMetJetDPhi >0.3||misc.MinMetJetDPhiIndex>3)", "MinDPhiCutMT2b");
	oname.ReplaceAll("(NEles==0  || ele[0].lv.Pt()<10)&&(NMuons==0 || muo[0].lv.Pt()<10)", "NoLeps");

	oname.ReplaceAll(">=" ,".ge");
	oname.ReplaceAll("<=" ,".le");
	oname.ReplaceAll(">" ,".gt");
	oname.ReplaceAll("<" ,".lt");
	oname.ReplaceAll("==",".eq");
	oname.ReplaceAll("!=",".ne");
	oname.ReplaceAll("&&","_");
	oname.ReplaceAll("||","_");
	oname.ReplaceAll("misc.","");
	oname.ReplaceAll("LeptConfig","LepCfg");
	oname.ReplaceAll("Vectorsumpt","VSPT");
	oname.ReplaceAll("EcalDeadCellBEFlag","BEFlg");
	oname.ReplaceAll("HBHENoiseFlag","HBHEFlg");
	oname.ReplaceAll("NJetsIDLoose","NJIDLoose");
	oname.ReplaceAll("isPFIDLoose","isJLoose");
	oname.ReplaceAll("IsGoodPFJet","IsGoodPFJ");
	oname.ReplaceAll("MinMetJetDPhi","MinDPhi");
	oname.ReplaceAll("Znunu.METplusLeptsPtReco","METLeptReco");
	oname.ReplaceAll("Znunu.MinMetplusLeptJetDPhiReco","MinMetLeptJetDPhi");
	oname.ReplaceAll("Znunu.caloHT50_matchedReco","caloHTmaReco");
	oname.ReplaceAll("Znunu.caloMHT30_matchedReco","caloMHTmReco");
	oname.ReplaceAll("Znunu.RecoOSmumu_mll","ROSmm_mll");
	oname.ReplaceAll("Znunu.RecoOSee_mll","ROSee_mll");
	oname.ReplaceAll("trigger.HLT_HT260_MHT60_v2","HLT_HT260_MHT60_v2");
	oname.ReplaceAll("trigger.HLT_HT250_MHT60_v3","HLT_HT250_MHT60_v3");
	oname.ReplaceAll("trigger.HLT_HT250_MHT60_v2","HLT_HT250_MHT60_v2");
	oname.ReplaceAll("trigger.HLT_HT440_v2","HLT_HT440_v2");
	oname.ReplaceAll("trigger.HLT_HT450_v2","HLT_HT450_v2");
	oname.ReplaceAll("trigger.HLT_HT500_v3","HLT_HT500_v3");
	oname.ReplaceAll(",","-");

	oname.ReplaceAll("trigger.", "");
	oname.ReplaceAll("LeadingJPt","J1Pt");
	oname.ReplaceAll("SecondJPt", "J2Pt");
	oname.ReplaceAll("Jet", "J");
	oname.ReplaceAll("GetNBtags", "NBtags");
	oname.ReplaceAll("bTagProb", "");
	oname.ReplaceAll("calo", "c");
	oname.ReplaceAll("CSCTightHaloID", "CSCHalo");
	oname.ReplaceAll("RecovRecHitFilterFlag", "RecovRecHit");
	oname.ReplaceAll("TrackingFailurePVtx", "TrkFailPVtx");
	oname.ReplaceAll("BadEcalBE", "BEFlg");
	oname.ReplaceAll("BadEcalTP", "TPFlg");

        TString outname = Util::removeFunnyChar(oname.Data());
	outname =  outname + (logflag ? "_log" : "");

	//do plotting (not TDR style)
	TCanvas *col = new TCanvas(outname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogz(1);
		h_mc_sum -> SetMinimum(0.05);
		h_data   -> SetMinimum(0.05);
		h_susy   -> SetMinimum(0.05);
	}else{
		h_mc_sum -> SetMinimum(0.0);
		h_data   -> SetMinimum(0.0);
		h_susy   -> SetMinimum(0.0);	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h_mc_sum->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h_mc_sum->SetMaximum(max);
	h_susy  ->SetMaximum(max);

	if(compareMCDATA){
		if(!(do3D)){
			h_mc_sum->Draw("COLZ");
			if(h_data->Integral()>0) {
				h_data       ->Draw("same TEXT");
			}
		} else{
			h_mc_sum->Draw("LEGO2Z");
			if(h_data->Integral()>0) {
				h_data       ->Draw("same LEGO");
			}
		}
	} else if(compareMCSUSY){
		if(!(do3D)){
			h_mc_sum->Draw("COLZ");
			if(h_susy->Integral()>0) {
				h_susy       ->Draw("same TEXT");
			}
		} else{
			h_mc_sum->Draw("LEGO2Z");
			if(h_susy->Integral()>0) {
				h_susy       ->Draw("same LEGO");
			}
		}
	} else if(compareSUSYDATA){
		if(!(do3D)){
			h_susy->Draw("COLZ");
			if(h_data->Integral()>0) {
				h_data       ->Draw("same TEXT");
			}
		} else{
			h_susy->Draw("LEGO2Z");
			if(h_data->Integral()>0) {
				h_data       ->Draw("same LEGO");
			}
		}
	}

	gPad->RedrawAxis();
	col ->Update();

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.05);
	TString text = title;
	TitleBox.DrawLatex(0.18,0.943,text.Data());

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.07);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.01, 0.9, ytitle);
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.1 : 0.06;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
	//Util::MakeOutputDir(fOutputDir)
	if(fSave)Util::PrintNoEPS(col, outname, fOutputDir, fOutputFile);
	if(fSave)Util::PrintEPS(col, outname, fOutputDir);

 	for(int i=0; i<h_samples.size(); ++i){
 		delete h_samples[i];
 	}
	h_samples.clear();

}

//read in samples.dat
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

			if(s.type=="mc" && runMC==false) continue;
			if(s.type=="data" && runData==false) continue;
			if(s.type=="susy" && runSUSY==false) continue;
			if(s.sname=="QCD" && runQCD==false) continue;
			if(s.sname=="Photons" && runPhotons==false) continue;
			if(s.sname=="Wtolnu" && runWJets==false) continue;
			if(s.shapename=="ZJetsToLL" && runZll==false) continue;
			if(s.shapename=="ZJetsToNuNu" && runZnunu==false) continue;
			if((s.name=="TTbar" || s.name=="TTZ" || s.name=="TTW") && runTTbar==false) continue;
			if((s.name=="Tbar_s" || s.name=="Tbar_t" || s.name=="Tbar_tW" || s.name=="T_s" || s.name=="T_t" || s.name=="T_tW") && runSingleTop==false) continue;

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