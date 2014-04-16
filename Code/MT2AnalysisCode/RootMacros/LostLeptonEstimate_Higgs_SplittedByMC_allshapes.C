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
#include "TLegend.h"
#include "TLegendEntry.h"
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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/Utilities.hh"

//use via root -l LostLeptonEstimate_Higgs_SplittedByMC_allshapes.C++

using namespace std;

void LostLeptonEstimate_Higgs_SplittedByMC_allshapes();



const int sampletypesize = 34;
string sample_type[sampletypesize] = {"Nominal", "MatchingUp", "MatchingDown", "ScaleUp", "ScaleDown", "MassUp", "MassDown", "JESUp", "JESDown", "METUp", "METDown", "BSFUp", "BSFDown", "ISRUp", "ISR", "ISRDown", "PUUp", "PUDown", "TopPtReweighted", "nominal0l", "nominal1l", "nominal2l", "nominal1Sample", "nominalPowheg", "WUp", "WDown", "TopUp", "TopDown", "SingleTopUp", "SingleTopDown", "TTVUp", "TTVDown", "MT2relaxedcut", "MT2fullcut"};//maybe not all filled  - ! include TTbar and WJets
const int samplekindsize = 6;
string sample_kind[samplekindsize] = {"allMC","allTop", "TTbar", "SingleTop", "TTV", "WJets"};
const int HTregionsize = 2;
string HT_region[HTregionsize] = {"lowHT", "highHT"};
//const int signalregionsize = 9;
//string signal_region[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b"};

bool fSave  = true;
std::ostringstream* fLogStream     = 0;

//this is an old copy (i.e. not cleaned up) of LostLeptonEstimate_Higgs_SplittedByMC.C
//this code does not reweight only the nominal shape of Mbb to the LostLeptonEstimation yield
//but also the modified shapes used for assessing the shape uncertainties of the LostLeptonEstimate.
void LostLeptonEstimate_Higgs_SplittedByMC_allshapes(){

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	bool fixemptybin = true; //if true use MC +/- 100% if data bin is empty.

	bool onlyttbar = false;//if true: onlytop = false, onlyW = false
	bool onlytop   = false;//if true: onlyW = false
	bool onlyW     = false;

	bool normalized = true;
	bool notnorm    = false;
	bool logflag    = false;

	fLogStream = new std::ostringstream();
	Bool_t  fWriteToFile               = false; // writes couts to file
	Bool_t  fAppend                    = true; // append at end of file (if existing), otherwise delete previous content


	bool WnoScaleMatching = true;//true if want W matchingup/down scaleup/down, Wincl sample in Nominal1Sample (default == false due to statistics of those W samples) // here only in plotting
	bool ISRusage = true;

	const int numberdifferentplots = 13;
	string differentplots[numberdifferentplots] = {"Nominal", "differentMC", "JES", "MET", "Scale", "Matching", "BSF", "PU" , "ISR", "Wxs", "Topxs", "MT2cuts", "bla"};
	bool   whattoplot[numberdifferentplots]     = {true,     false,          true,  true,  true,    true,       true,  true,   true, true,  true,    true,          false};

	TString   inputdir = "Filtered/TTbarStudies/Higgs/";
	TString  inputname = "TTbarStudiesHistograms_all_ISR_limitMbb.root";
	if(WnoScaleMatching)  inputname = "TTbarStudiesHistograms_all_noWscaleupdown_ISR_limitMbb.root";
	TString  outputname = "LostLeptonHiggs_ShapesScaledByEstimate_all_ISR_limitMbb.root";
	if(WnoScaleMatching)  outputname = "LostLeptonHiggs_ShapesScaledByEstimate_all_noWscaleupdown_ISR_limitMbb.root";

	TString outputdir  = "Filtered/TTbarStudies/Higgs/Plots/all/ISR/";
	if(WnoScaleMatching) outputdir  = "Filtered/TTbarStudies/Higgs/Plots/all_noWscaleupdown/ISR/";
	if(onlyttbar)    outputdir = outputdir + "ttbar/";
	else if(onlytop) outputdir = outputdir + "top/";
	else if(onlyW)   outputdir = outputdir + "W/";
	if(logflag) outputdir = outputdir + "log/";
	else        outputdir = outputdir + "linear/";
    	Util::MakeOutputDir(outputdir);

	TFile *infile = TFile::Open(inputdir + inputname);

	map<string, TH1D*>    histos;
	map<string, TH1D*>    histosnormalized;
	vector<string> histonamestosave; histonamestosave.clear();

	string samplekind;
	if(onlyttbar)    samplekind = "TTbar";
	else if(onlytop) samplekind = "allTop";
	else if(onlyW)   samplekind = "WJets";
	else             samplekind = "allMC";
	for(int i1 = 0; i1<sampletypesize;   ++i1){
//	for(int i2 = 0; i2<signalregionsize; ++i2){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		string hs = string("_") + sample_type[i1] + string("_") + samplekind + string("_") + HT_region[i3];// + string("_") + signal_region[i2];
		string mapname = "MT2" + hs;
		if(histos.count(mapname) == 0 ) histos[mapname] = (TH1D*)infile->Get(mapname.c_str());
		string mapnamenorm = "Norm"+mapname;
		if(histosnormalized.count(mapnamenorm) == 0 ) histosnormalized[mapnamenorm] = (TH1D*)infile->Get(mapnamenorm.c_str());
	}}//}
    TH1D *hLLest    = new TH1D("hLLest"   ,"",6,0,6); hLLest   ->SetMarkerStyle(20), hLLest   ->SetMarkerColor(kBlack); hLLest->SetLineWidth(3); hLLest->SetLineColor(kBlack);
    TH1D *hLLest_MC = new TH1D("hLLest_MC","",6,0,6); hLLest_MC->SetFillStyle(3001); hLLest_MC->SetFillColor(kBlue);

    TLegend *leg = new TLegend(0.6551724,0.7299578,0.8706897,0.8987342,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(62);
    leg->SetTextSize(0.04575163);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->AddEntry(hLLest_MC, "MC truth", "f");
    leg->AddEntry(hLLest, "data prediction", "lp");
    
    TH1D *haxis = new TH1D("haxis","",6,0,6);

    haxis->GetXaxis()->SetBinLabel(1, "e, low H_{T}");
    haxis->GetXaxis()->SetBinLabel(2, "#mu, low H_{T}");
    haxis->GetXaxis()->SetBinLabel(3, "#tau, low H_{T}");
    haxis->GetXaxis()->SetBinLabel(4, "e, high H_{T}");
    haxis->GetXaxis()->SetBinLabel(5, "#mu, high H_{T}");
    haxis->GetXaxis()->SetBinLabel(6, "#tau, high H_{T}");
    haxis->GetXaxis()->SetTitle("signal channel");
    haxis->GetXaxis()->SetLabelFont(42);    haxis->GetXaxis()->SetLabelSize(0.06);    haxis->GetXaxis()->SetTitleSize(0.05);    haxis->GetXaxis()->SetTitleFont(42);
    haxis->GetYaxis()->SetTitle("Lost Lepton yield");
    haxis->GetYaxis()->SetLabelFont(42);    haxis->GetYaxis()->SetLabelSize(0.05);    haxis->GetYaxis()->SetTitleSize(0.05);    haxis->GetYaxis()->SetTitleFont(42);
    haxis->GetZaxis()->SetLabelFont(42);    haxis->GetZaxis()->SetLabelSize(0.035);   haxis->GetZaxis()->SetTitleSize(0.035);   haxis->GetZaxis()->SetTitleFont(42);

    hLLest_MC->SetBinContent(1,14.30); hLLest_MC->SetBinError(1,1.65);// e,   lHT
    hLLest   ->SetBinContent(1,15.54); hLLest   ->SetBinError(1,sqrt(pow(4.33,2)+pow(2.43,2)+pow(0.25,2)+pow(0.,2)));
    hLLest_MC->SetBinContent(2,11.15); hLLest_MC->SetBinError(2,1.41);// mu,  lHT
    hLLest   ->SetBinContent(2,12.49); hLLest   ->SetBinError(2,sqrt(pow(3.05,2)+pow(2.19,2)+pow(0.12,2)+pow(0.,2)));
    hLLest_MC->SetBinContent(3,18.31); hLLest_MC->SetBinError(3,1.81);// tau, lHT
//    hLLest   ->SetBinContent(3,6.83); hLLest   ->SetBinError(3,sqrt(pow(6.30,2)+pow(1.46,2)+pow(0.72,2)+pow(0.,2)));
    hLLest   ->SetBinContent(3,9.05); hLLest   ->SetBinError(3,sqrt(pow(6.26,2)+pow(1.60,2)+pow(0.11,2)+pow(0.,2)));
    hLLest_MC->SetBinContent(4,18.98); hLLest_MC->SetBinError(4,1.64);// e,   hHT
    hLLest   ->SetBinContent(4,11.86); hLLest   ->SetBinError(4,sqrt(pow(3.79,2)+pow(1.67,2)+pow(0.12,2)+pow(0.,2)));
    hLLest_MC->SetBinContent(5,17.96); hLLest_MC->SetBinError(5,2.44);// mu,  hHT
    hLLest   ->SetBinContent(5,22.08); hLLest   ->SetBinError(5,sqrt(pow(4.04,2)+pow(3.40,2)+pow(0.11,2)+pow(0.,2)));
    hLLest_MC->SetBinContent(6,29.98); hLLest_MC->SetBinError(6,3.78);// tau, hHT
    hLLest   ->SetBinContent(6,30.87); hLLest   ->SetBinError(6,sqrt(pow(14.04,2)+pow(5.24,2)+pow(0.32,2)+pow(0.,2)));

	for(int i3 = 0; i3<HTregionsize;     ++i3){
	for(int i1 = 0; i1<numberdifferentplots; ++i1){
		string hs0, hs1, hs2, hs3; int numhistos;
		string hs = string("_") + samplekind + string("_") + HT_region[i3];// + string("_") + signal_region[i2];
		hs0 = "_Nominal"+hs;
		if(whattoplot[i1]==false) continue;
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
		if(differentplots[i1]=="Restxs")      { numhistos = 2; hs1 = "_TTVUp"+hs;          hs2 = "_TTVDown"+hs;      }
		if(numhistos<3) hs3 = "_Nominal"+hs;//dummy
		if(numhistos<2) hs2 = "_Nominal"+hs;//dummy
		if(numhistos<1) hs1 = "_Nominal"+hs;//dummy
		double scale(1.), scaleerr(0.);
		double scaleMC(1.), scaleerrMC(0.);
		if(HT_region[i3]=="lowHT"){
			double ele    = hLLest->GetBinContent(1);
			double eleerr = hLLest->GetBinError(  1);
			double muo    = hLLest->GetBinContent(2);
			double muoerr = hLLest->GetBinError(  2);
			double tau    = hLLest->GetBinContent(3);
			double tauerr = hLLest->GetBinError(  3);
			double MCele    = hLLest_MC->GetBinContent(1);
			double MCeleerr = hLLest_MC->GetBinError(  1);
			double MCmuo    = hLLest_MC->GetBinContent(2);
			double MCmuoerr = hLLest_MC->GetBinError(  2);
			double MCtau    = hLLest_MC->GetBinContent(3);
			double MCtauerr = hLLest_MC->GetBinError(  3);
			if(fixemptybin&&ele==0) { ele = MCele; eleerr = MCele;}
			if(fixemptybin&&muo==0) { muo = MCmuo; eleerr = MCmuo;}
			if(fixemptybin&&tau==0) { tau = MCtau; eleerr = MCtau;}
			scaleMC = MCele+MCmuo+MCtau;
			scaleerrMC = sqrt(MCeleerr*MCeleerr+MCmuoerr*MCmuoerr+MCtauerr*MCtauerr);
			scale = ele+muo+tau;
			scaleerr = sqrt(eleerr*eleerr+muoerr*muoerr+tauerr*tauerr);
		}
		if(HT_region[i3]=="highHT"){
			double ele    = hLLest->GetBinContent(4);
			double eleerr = hLLest->GetBinError(  4);
			double muo    = hLLest->GetBinContent(5);
			double muoerr = hLLest->GetBinError(  5);
			double tau    = hLLest->GetBinContent(6);
			double tauerr = hLLest->GetBinError(  6);
			double MCele    = hLLest_MC->GetBinContent(4);
			double MCeleerr = hLLest_MC->GetBinError(  4);
			double MCmuo    = hLLest_MC->GetBinContent(5);
			double MCmuoerr = hLLest_MC->GetBinError(  5);
			double MCtau    = hLLest_MC->GetBinContent(6);
			double MCtauerr = hLLest_MC->GetBinError(  6);
			if(fixemptybin&&ele==0) { ele = MCele; eleerr = MCele;}
			if(fixemptybin&&muo==0) { muo = MCmuo; eleerr = MCmuo;}
			if(fixemptybin&&tau==0) { tau = MCtau; eleerr = MCtau;}
			scaleMC = MCele+MCmuo+MCtau;
			scaleerrMC = sqrt(MCeleerr*MCeleerr+MCmuoerr*MCmuoerr+MCtauerr*MCtauerr);
			scale = ele+muo+tau;
			scaleerr = sqrt(eleerr*eleerr+muoerr*muoerr+tauerr*tauerr);
		}
		double squaredsum = 0;
		if(i1==0){
			for(int n = 1; n<= histosnormalized["NormMT2"+hs0]->GetNbinsX(); ++n){
				double content = histosnormalized["NormMT2"+hs0]->GetBinContent(n);
				histosnormalized["NormMT2"+hs0]->SetBinContent(n, content*scale);
			//	squaredsum += pow(content*scale,2);
				//correct error1:
				// Delta(a)/a = Delta(b)/b = Delta(c)/c = Delta(tot)/tot;
				histosnormalized["NormMT2"+hs0]->SetBinError(  n, content*scaleerr);
				bool existing = false;
				for(int m = 0; m<histonamestosave.size();++m){
			//		cout << "NormMT2"+hs0 << endl;
					if(histonamestosave[m]=="NormMT2"+hs0) existing = true;
			//		if(!existing) cout << "not yet there" << endl; else cout << "is there" << endl;
				}
				if(!existing) histonamestosave.push_back("NormMT2"+hs0);
			}
		}
		if(numhistos>=1){
			for(int n = 1; n<= histosnormalized["NormMT2"+hs1]->GetNbinsX(); ++n){
				double content = histosnormalized["NormMT2"+hs1]->GetBinContent(n);
				histosnormalized["NormMT2"+hs1]->SetBinContent(n, content*scale);
			//	squaredsum += pow(content*scale,2);
				histosnormalized["NormMT2"+hs1]->SetBinError(  n, content*scaleerr);
				bool existing = false;
				for(int m = 0; m<histonamestosave.size();++m){
					if(histonamestosave[m]=="NormMT2"+hs1) existing = true;
				}
				if(!existing) histonamestosave.push_back("NormMT2"+hs1);
			}
		}
		if(numhistos>=2){
			for(int n = 1; n<= histosnormalized["NormMT2"+hs2]->GetNbinsX(); ++n){
				double content = histosnormalized["NormMT2"+hs2]->GetBinContent(n);
				histosnormalized["NormMT2"+hs2]->SetBinContent(n, content*scale);
			//	squaredsum += pow(content*scale,2);
				histosnormalized["NormMT2"+hs2]->SetBinError(  n, content*scaleerr);
				bool existing = false;
				for(int m = 0; m<histonamestosave.size();++m){
					if(histonamestosave[m]=="NormMT2"+hs2) existing = true;
				}
				if(!existing) histonamestosave.push_back("NormMT2"+hs2);
			}
		}
		if(numhistos>=3){
			for(int n = 1; n<= histosnormalized["NormMT2"+hs3]->GetNbinsX(); ++n){
				double content = histosnormalized["NormMT2"+hs3]->GetBinContent(n);
				histosnormalized["NormMT2"+hs3]->SetBinContent(n, content*scale);
			//	squaredsum += pow(content*scale,2);
				histosnormalized["NormMT2"+hs3]->SetBinError(  n, content*scaleerr);
				bool existing = false;
				for(int m = 0; m<histonamestosave.size();++m){
					if(histonamestosave[m]=="NormMT2"+hs3) existing = true;
				}
				if(!existing) histonamestosave.push_back("NormMT2"+hs3);
			}
		}
	}}
    	TFile *fsavefile = new TFile(inputdir + outputname,"RECREATE");
	fsavefile->cd();
	for(int i2 = 0; i2<histonamestosave.size(); ++i2){
		histosnormalized[histonamestosave[i2] ]->Write();
	}
	hLLest->Write();
	hLLest_MC->Write();
	cout << "files saved in " << fsavefile->GetName() << endl;
}
