#include "TEventList.h"
#include "TLorentzVector.h"
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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MassPlotter.hh"//use printHisto macro

using namespace std;

//run via root -l -b -q TauStudies_Efficiency.C++

//subset of TauStudies.C
//study the tau efficiency using histograms from root file done with TauStudies.C
void TauStudies_Efficiency(){

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

  	gROOT->ProcessLine(".x SetStyle_PRD.C");
	
	//this flag are identical to the ones in TauStudies.C
	bool MET  = true;
	bool HT   = false;
	bool T2bb = true;
	bool T2tt = false;
	bool correcttau = true;
	
	//output directory
	TString outputdir = "dummy";
	if(HT){
		if(T2tt) outputdir = "data_11fb-1_HT/T2tt/TauStudies/";
		if(T2bb) outputdir = "data_11fb-1_HT/T2bb/TauStudies/";
	}
	if(MET){
		if(T2tt) outputdir = "data_11fb-1_MET/T2tt/TauStudies/";
		if(T2bb) outputdir = "data_11fb-1_MET/T2bb/TauStudies/";
	}
	if(correcttau) outputdir = outputdir + "CorrectTauPt/";

    	Util::MakeOutputDir(outputdir);


	map<string, TH1D*>    histos;


    TFile *fsavefile = TFile::Open(outputdir + "Histograms.root");//do not delete this one


    const int signalregionsize = 9;
    string signalregions[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b"};

	for(int js = 0; js<signalregionsize; ++js){
		string hs = string("_") + signalregions[js] + string("_mc");
		string mapname;
		mapname = "MT2_RecoTausge1" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)fsavefile->Get(mapname.c_str());
		mapname = "MT2_GenTausge1" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)fsavefile->Get(mapname.c_str());
		mapname = "MT2_GenHadTausge1" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)fsavefile->Get(mapname.c_str());
		mapname = "MT2_GenTausge1_acc" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)fsavefile->Get(mapname.c_str());
		mapname = "MT2_GenHadTausge1_acc" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)fsavefile->Get(mapname.c_str());
		mapname = "MT2_RecoTau_matchGenTau" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)fsavefile->Get(mapname.c_str());
		mapname = "MT2_RecoTau_matchGenHadTau" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)fsavefile->Get(mapname.c_str());
	}

	//do the efficiencies
	TH1D *passgentaus[signalregionsize];
	TH1D  *totgentaus[signalregionsize];
	TH1D *passgentaus_acc[signalregionsize];
	TH1D  *totgentaus_acc[signalregionsize];
	TH1D *passgenhadtaus[signalregionsize];
	TH1D  *totgenhadtaus[signalregionsize];
	TH1D *passgenhadtaus_acc[signalregionsize];
	TH1D  *totgenhadtaus_acc[signalregionsize];

	TH1D *matchedpassgentaus[signalregionsize];
	TH1D  *matchedtotgentaus[signalregionsize];
	TH1D *matchedpassgentaus_acc[signalregionsize];
	TH1D  *matchedtotgentaus_acc[signalregionsize];
	TH1D *matchedpassgenhadtaus[signalregionsize];
	TH1D  *matchedtotgenhadtaus[signalregionsize];
	TH1D *matchedpassgenhadtaus_acc[signalregionsize];
	TH1D  *matchedtotgenhadtaus_acc[signalregionsize];

	TEfficiency *taueff[signalregionsize];
	TEfficiency *hadtaueff[signalregionsize];
	TEfficiency *tauwithinacceff[signalregionsize];
	TEfficiency *hadtauwithinacceff[signalregionsize];
	TEfficiency *matchedtaueff[signalregionsize];
	TEfficiency *matchedhadtaueff[signalregionsize];
	TEfficiency *matchedtauwithinacceff[signalregionsize];
	TEfficiency *matchedhadtauwithinacceff[signalregionsize];

	for(int il = 0; il<signalregionsize; ++il){
		string pgt    = "MT2_RecoTausge1_"            + signalregions[il] + "_mc"; passgentaus[il] = (TH1D*)histos[pgt]->Clone((pgt+(string)"_copy").c_str());
		string tgt    = "MT2_GenTausge1_"             + signalregions[il] + "_mc"; totgentaus[il] = (TH1D*)histos[tgt]->Clone((tgt+(string)"_copy").c_str());
		string pgta   = "MT2_RecoTausge1_"            + signalregions[il] + "_mc"; passgentaus_acc[il] = (TH1D*)histos[pgta]->Clone((pgta+(string)"_copy").c_str());
		string tgta   = "MT2_GenTausge1_acc_"         + signalregions[il] + "_mc"; totgentaus_acc[il] = (TH1D*)histos[tgta]->Clone((tgta+(string)"_copy").c_str());
		string pght   = "MT2_RecoTausge1_"            + signalregions[il] + "_mc"; passgenhadtaus[il] = (TH1D*)histos[pght]->Clone((pght+(string)"_copy").c_str());
		string tght   = "MT2_GenHadTausge1_"          + signalregions[il] + "_mc"; totgenhadtaus[il] = (TH1D*)histos[tght]->Clone((tght+(string)"_copy").c_str());
		string pghta  = "MT2_RecoTausge1_"            + signalregions[il] + "_mc"; passgenhadtaus_acc[il] = (TH1D*)histos[pghta]->Clone((pghta+(string)"_copy").c_str());
		string tghta  = "MT2_GenHadTausge1_acc_"      + signalregions[il] + "_mc"; totgenhadtaus_acc[il] = (TH1D*)histos[tghta]->Clone((tghta+(string)"_copy").c_str());

		string mpgt   = "MT2_RecoTau_matchGenTau_"    + signalregions[il] + "_mc"; matchedpassgentaus[il] = (TH1D*)histos[mpgt]->Clone((mpgt+(string)"_copy").c_str());
		string mtgt   = "MT2_GenTausge1_"             + signalregions[il] + "_mc"; matchedtotgentaus[il] = (TH1D*)histos[mtgt]->Clone((mtgt+(string)"_copy").c_str());
		string mpgta  = "MT2_RecoTau_matchGenTau_"    + signalregions[il] + "_mc"; matchedpassgentaus_acc[il] = (TH1D*)histos[mpgta]->Clone((mpgta+(string)"_copy").c_str());
		string mtgta  = "MT2_GenTausge1_acc_"         + signalregions[il] + "_mc"; matchedtotgentaus_acc[il] = (TH1D*)histos[mtgta]->Clone((mtgta+(string)"_copy").c_str());
		string mpght  = "MT2_RecoTau_matchGenHadTau_" + signalregions[il] + "_mc"; matchedpassgenhadtaus[il] = (TH1D*)histos[mpght]->Clone((mpght+(string)"_copy").c_str());
		string mtght  = "MT2_GenHadTausge1_"          + signalregions[il] + "_mc"; matchedtotgenhadtaus[il] = (TH1D*)histos[mtght]->Clone((mtght+(string)"_copy").c_str());
		string mpghta = "MT2_RecoTau_matchGenHadTau_" + signalregions[il] + "_mc"; matchedpassgenhadtaus_acc[il] = (TH1D*)histos[mpghta]->Clone((mpghta+(string)"_copy").c_str());
		string mtghta = "MT2_GenHadTausge1_acc_"      + signalregions[il] + "_mc"; matchedtotgenhadtaus_acc[il] = (TH1D*)histos[mtghta]->Clone((mtghta+(string)"_copy").c_str());

		taueff[il] = new TEfficiency((*(passgentaus[il])), (*(totgentaus[il])));
		taueff[il]->SetNameTitle(("taueff"+signalregions[il]).c_str(), ("taueff"+signalregions[il]).c_str());
		tauwithinacceff[il] = new TEfficiency((*(passgentaus_acc[il])), (*(totgentaus_acc[il])));
		tauwithinacceff[il]->SetNameTitle(("tauwithinacceff"+signalregions[il]).c_str(), ("tauwithinacceff"+signalregions[il]).c_str());
		hadtaueff[il] = new TEfficiency((*(passgenhadtaus[il])), (*(totgenhadtaus[il])));
		hadtaueff[il]->SetNameTitle(("hadtaueff"+signalregions[il]).c_str(), ("hadtaueff"+signalregions[il]).c_str());
		hadtauwithinacceff[il] = new TEfficiency((*(passgenhadtaus_acc[il])), (*(totgenhadtaus_acc[il])));
		hadtauwithinacceff[il]->SetNameTitle(("hadtauwithinacceff"+signalregions[il]).c_str(), ("hadtauwithinacceff"+signalregions[il]).c_str());

		matchedtaueff[il] = new TEfficiency((*(matchedpassgentaus[il])), (*(matchedtotgentaus[il])));
		matchedtaueff[il]->SetNameTitle(("matchedtaueff"+signalregions[il]).c_str(), ("matchedtaueff"+signalregions[il]).c_str());
		matchedtauwithinacceff[il] = new TEfficiency((*(matchedpassgentaus_acc[il])), (*(matchedtotgentaus_acc[il])));
		matchedtauwithinacceff[il]->SetNameTitle(("matchedtauwithinacceff"+signalregions[il]).c_str(), ("matchedtauwithinacceff"+signalregions[il]).c_str());
		matchedhadtaueff[il] = new TEfficiency((*(matchedpassgenhadtaus[il])), (*(matchedtotgenhadtaus[il])));
		matchedhadtaueff[il]->SetNameTitle(("matchedhadtaueff"+signalregions[il]).c_str(), ("matchedhadtaueff"+signalregions[il]).c_str());
		matchedhadtauwithinacceff[il] = new TEfficiency((*(matchedpassgenhadtaus_acc[il])), (*(matchedtotgenhadtaus_acc[il])));
		matchedhadtauwithinacceff[il]->SetNameTitle(("matchedhadtauwithinacceff"+signalregions[il]).c_str(), ("matchedhadtauwithinacceff"+signalregions[il]).c_str());
	}
	//and plot the efficiencies
	TCanvas* c1 = new TCanvas("dummy","",0,0,600,600 /*37, 60,636,670*/);
	for(int il = 0; il<signalregionsize; ++il){
		TString canvname;
		c1->cd();
		taueff[il]->Draw();
		c1->Update();
		canvname = taueff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		hadtaueff[il]->Draw();
		c1->Update();
		canvname = hadtaueff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		tauwithinacceff[il]->Draw();
		c1->Update();
		canvname = tauwithinacceff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		hadtauwithinacceff[il]->Draw();
		c1->Update();
		canvname = hadtauwithinacceff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		matchedtaueff[il]->Draw();
		c1->Update();
		canvname = matchedtaueff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		matchedhadtaueff[il]->Draw();
		c1->Update();
		canvname = matchedhadtaueff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		matchedtauwithinacceff[il]->Draw();
		c1->Update();
		canvname = matchedtauwithinacceff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
		c1->cd();
		matchedhadtauwithinacceff[il]->Draw();
		c1->Update();
		canvname = matchedhadtauwithinacceff[il]->GetName();
		Util::PrintNoEPS(c1, canvname, outputdir);
		Util::PrintEPS(c1, canvname, outputdir);
		c1->Clear();
	}
}
