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
#include "TChain.h"
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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"// use your own path of MT2tree.hh
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"

using namespace std;

//call this macro via root -l -b -q GammaVsZnunuStudies.C++
void GammaVsZnunuStudies();
void load(const char* filename = "samples/datasamples/samples_2141_dataonly.dat");

//struct combines the MT2tree with necesarry information like sample type, xsection,etc.
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

//This code is a modification of GammaVsZllStudies - see that macro first!
//Here we get the histograms for Znunu/gamma MC ratios
void GammaVsZnunuStudies(){
gROOT->cd();
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

//for explanation of all these inputs see GammaVsZllStudies
double minphotonpt = 20;
double minzpt      = 20;
double maxphotonpt = -1;
double maxzpt      = -1;
double minHT       = 1200;
double maxHT       = -1;
double minMT2      = -1;
double maxMT2      = -1;
double METcutval   = -1;
int    njets       = -2;
int    nbjets      = 0;

TString plotvariable = "misc.MT2";//choose "VPt" to plot all VPts!
TString varname      = "MT2";//VPt // this is for saving the histograms!
TString varlabel     = "M_{T2} [GeV]";//V-p_{T} [GeV]
int    nbins       = 40;
double binlow      = 0;
double binup       = 1000;
//6 histograms: Z, Z_bg, G, G_bg, Z_data, G_data
//VPt, MT2: 30, 0, 750;   HT: 20, 450,1450;   NJets: 8, 2, 10;   NBJets: 5, 0, 5

int    HTselection = 1;//1:allHT - only for ll, 2: highPt, for ll and g, 3: highHT
bool   addoverflow = true;
bool   addunderflow= false;//keep that false
bool   removebkg   = false;//keep that false, can do it later

TString samples_2l; TString samples_1g;//files to load for photon or dilepton selection
if(HTselection==2)             samples_2l = "samples/samples_Znunu_HTMET.dat"; samples_1g = "samples/samples_1g_GEst_noHTdata.dat";
if(HTselection==3)             samples_2l = "samples/samples_Znunu_HTMET.dat"; samples_1g = "samples/samples_1g_GEst_HTdata.dat";
if(HTselection==1)             samples_2l = "samples/samples_Znunu_HTMET.dat"; samples_1g = "samples/samples_1g_gRemoved_HTMET.dat";//VPt>80 for g-data // VPt>125 for Z-MC
if(HTselection==2&&maxHT==750) samples_2l = "samples/samples_Znunu_HTMET.dat"; samples_1g = "samples/samples_1g_GEst_MET.dat";

TString outputdir  = "GammaVsZnunuStudies/";

TString preselectioncuts = "NJetsIDLoose40>=2&&misc.HT>450&&NTausIDLoose3Hits==0&&NEles==0&&NMuons==0&&misc.Vectorsumpt<70";
TString filters = "misc.PassJet40ID ==1&&(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&(type1pfmet[0].Pt()>30||type1pfmet[0].Pt()/misc.CaloMETRaw<=2.)";
TString METcut           = TString::Format("misc.MET>=%d", (int)METcutval);
TString MinDPhi_g        = "misc.MinMetJetDPhi4Pt40>0.3";
TString MinDPhi_Z        = "MinMetJetDPhi(1,40.,2.4,1,4)>0.3";
TString Znunusel         = TString::Format("GenZ[0].Pt()>0");
TString Znunuvar         = TString::Format("GenZ[0].Pt()>=%d", (int)minzpt);
TString upper_Znunuvar   = TString::Format("GenZ[0].Pt()<%d", (int)maxzpt);
TString Gsel             = TString::Format("GenPhoton[0].Pt()>0");
TString Gvar             = TString::Format("GenPhoton[0].Pt()>=%d", (int)minzpt);
TString upper_Gvar       = TString::Format("GenPhoton[0].Pt()<%d", (int)maxzpt);
TString HTlow            = TString::Format("misc.HT>=%d", (int)minHT);
TString HTup             = TString::Format("misc.HT<=%d", (int)maxHT);
TString MT2low           = TString::Format("misc.MT2>=%d", (int)minMT2);
TString MT2up            = TString::Format("misc.MT2<=%d", (int)maxMT2);
TString NJetsCut         = "NJetsIDLoose40>=2";
if (njets>=10) NJetsCut  = TString::Format("NJetsIDLoose40>=%d",njets/10)+"&&"+TString::Format("NJetsIDLoose40<=%d",njets%10);
else           NJetsCut  = TString::Format("NJetsIDLoose40") + (njets < 0 ? TString::Format(">=") : TString::Format("==")) + TString::Format("%d",abs(njets));
TString NBJetsCut        = "NBJets40CSVM";
NBJetsCut               += nbjets < 0 ? ">=" : "==";
NBJetsCut               += nbjets==-10 ? "0" : TString::Format("%d",abs(nbjets));

//change varname according to cuts
                        varname += "_SF_";//SF is a dummy
if(njets>=0)            varname += TString::Format("%dj"  ,   abs(njets));
else if(njets!=-10)     varname += TString::Format("ge%dj",   abs(njets));
if(njets!=-10)          varname += "_";
if(nbjets>=0)           varname += TString::Format("%db"  ,  abs(nbjets));
else if(nbjets!=-10)    varname += TString::Format("ge%db",  abs(nbjets));
if(nbjets!=-10)         varname += "_";
if(minzpt>0)            varname += TString::Format("VPtge%d",(int)minzpt);
if(minzpt>0&&maxzpt>0)  varname += TString::Format(   "le%d",(int)maxzpt);
else if(     maxzpt>0)  varname += TString::Format("VPtle%d",(int)maxzpt);
if(minzpt>0||maxzpt>0)  varname += "_";
if(minHT >0)            varname += TString::Format( "HTge%d",(int)minHT );
if(minHT >0&&maxHT >0)  varname += TString::Format(   "le%d",(int)maxHT );
else if(     maxHT >0)  varname += TString::Format( "HTle%d",(int)maxHT );
if(minHT >0||maxHT >0)  varname += "_";
if(minMT2>0)            varname += TString::Format("MT2ge%d",(int)minMT2);
if(minMT2>0&&maxMT2>0)  varname += TString::Format(   "le%d",(int)maxMT2);
else if(     maxMT2>0)  varname += TString::Format("MT2le%d",(int)maxMT2);
if(minMT2>0||maxMT2>0)  varname += "_";
if(HTselection==1)      varname += "all";
if(HTselection==2)      varname += "highVPt";
if(HTselection==3)      varname += "highHT";

//define the histograms
//make ratio data/MC (2), make ratio G/Z (also 2), rest in postprocessing of stupid code (e.g. 0b/1b etc.)
TString histoname = varname + "_Z";
TH1D *hZ            = new TH1D(histoname.Data(), "", nbins, binlow, binup); hZ              ->Sumw2();
TString ytitle = TString::Format("Events / %d GeV", (int)hZ->GetBinWidth(1));
  hZ           ->GetXaxis()->SetTitle(varlabel.Data()); hZ           ->GetYaxis()->SetTitle(ytitle.Data());
  hZ           ->SetLineColor(kOrange  ); hZ           ->SetMarkerColor(kOrange  );   hZ           ->SetMarkerStyle(22);
histoname = varname + "_G";
TH1D *hG            = new TH1D(histoname.Data(), "", nbins, binlow, binup); hG           ->Sumw2();
  hG           ->GetXaxis()->SetTitle(varlabel.Data()); hG           ->GetYaxis()->SetTitle(ytitle.Data());
  hG           ->SetLineColor(kViolet-3); hG           ->SetMarkerColor(kViolet-3);   hG           ->SetMarkerStyle(20);
histoname = varname + "_ZG_MCratio";
TH1D *hZG_MCratio   = new TH1D(histoname.Data(), "", nbins, binlow, binup); hZG_MCratio  ->Sumw2();
ytitle = TString::Format("#frac{Z(ll)}{#gamma} / %d GeV", (int)hZ->GetBinWidth(1));
  hZG_MCratio  ->GetXaxis()->SetTitle(varlabel.Data()); hZG_MCratio  ->GetYaxis()->SetTitle(ytitle.Data());
  hZG_MCratio  ->SetLineColor(kViolet-3); hZG_MCratio  ->SetMarkerColor(kViolet-3);   hZG_MCratio  ->SetMarkerStyle(20);

//make two forloops, for Znunu (called Zll here, but it is Znunu) and Gamma
load(samples_2l.Data());
for(size_t i = 0; i < fSamples.size(); ++i){ //ZLL
   string sampletype = "";

   if(fSamples[i].sname=="ZJetsToNuNu") sampletype = "Z_inv";
   else if(fSamples[i].shapename=="ZJetsToNuNu") sampletype = "Z";
   else { cout << "sample " << fSamples[i].name << " is not considered" << endl; continue; }

    MT2tree* fMT2tree = new MT2tree();
    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    int nev =0;
    TString myCuts = preselectioncuts + "&&" + filters + "&&" + NJetsCut + "&&" + NBJetsCut;

    if((int)minHT >0) myCuts = myCuts+ "&&" + HTlow;
    if((int)maxHT >0 && maxHT >minHT ) myCuts = myCuts+ "&&" + HTup;
    if(HTselection==3 && minHT<750.  ) myCuts = myCuts+ "&&" + "misc.HT>=750";
    if((int)minMT2>0) myCuts = myCuts+ "&&" + MT2low;
    if((int)maxMT2>0 && maxMT2>minMT2) myCuts = myCuts+ "&&" + MT2up;
       if(minzpt>0)                myCuts = myCuts+ "&&" + Znunuvar;
       else                        myCuts = myCuts+ "&&" + Znunusel;
       if(maxzpt>0&&maxzpt>minzpt) myCuts = myCuts+ "&&" + upper_Znunuvar;
 

    double weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
    if(fVerbose>2) cout << "===========================================================" <<  endl;
    if(fVerbose>2) cout << "MakePlot: looping over " << fSamples[i].name << "-----------------------------------" <<  endl;
    if(fVerbose>2) cout << "  +++++++ weight:      " << weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 

    TString btagweight = "1.00"; //stored btag weights up to >=3, default is weight==1 to avoid non-existing weights
    if(nbjets>=0 && nbjets<=3) btagweight = TString::Format("SFWeight.BTagCSV40eq%d",abs(nbjets));
    else if(nbjets>=-3)        btagweight = TString::Format("SFWeight.BTagCSV40ge%d",abs(nbjets));

    TString selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), myCuts.Data());
    TString variable = TString::Format("%s>>+%s",plotvariable.Data(),hZ     ->GetName());
    else continue;
    if(plotvariable=="VPt"){
       TString tempvar = "GenZ[0].Pt()";
       variable = TString::Format("%s>>+%s",tempvar.Data(),hZ     ->GetName());
       else continue;
    }
    gROOT->cd();
    if(fVerbose>2) cout << "  +++++++ Drawing      " << variable  << endl
                        << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
    nev = fSamples[i].tree->Draw(variable.Data(),selection.Data(),"goff");//this fills the histograms
    if(fVerbose>2) cout << "  +++++++ MC   events : "  <<  nev << endl;

}//for(size_t i = 0; i < fSamples.size(); ++i) //ZLL

load(samples_1g.Data());
gROOT->cd();
for(size_t i = 0; i < fSamples.size(); ++i){ //gamma
   string sampletype = "";
   if(fSamples[i].sname=="Photons")          sampletype = "G";
   else { cout << "sample " << fSamples[i].name << " is not considered" << endl; continue; }

    MT2tree* fMT2tree = new MT2tree();
    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    int nev =0;
    TString myCuts = preselectioncuts + "&&" + filters + "&&" + NJetsCut + "&&" + NBJetsCut;
    if((int)minHT >0) myCuts = myCuts+ "&&" + HTlow;
    if((int)maxHT >0 && maxHT >minHT ) myCuts = myCuts+ "&&" + HTup;
    if(HTselection==3 && minHT<750.  ) myCuts = myCuts+ "&&" + "misc.HT>=750";
    if((int)minMT2>0) myCuts = myCuts+ "&&" + MT2low;
    if((int)maxMT2>0 && maxMT2>minMT2) myCuts = myCuts+ "&&" + MT2up;
    if(minphotonpt>0) myCuts = myCuts+ "&&" + Gvar;
    else              myCuts = myCuts+ "&&" + Gsel;
    if(maxphotonpt>0&&maxphotonpt>minphotonpt) myCuts = myCuts+ "&&" + upper_Gvar;

    double weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
    if(fVerbose>2) cout << "===========================================================" <<  endl;
    if(fVerbose>2) cout << "MakePlot: looping over " << fSamples[i].name << "-----------------------------------" <<  endl;
    if(fVerbose>2) cout << "  +++++++ weight:      " << weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 

    TString btagweight = "1.00"; //stored btag weights up to >=3, default is weight==1 to avoid non-existing weights
    if(nbjets>=0 && nbjets<=3) btagweight = TString::Format("SFWeight.BTagCSV40eq%d",abs(nbjets));
    else if(nbjets>=-3)        btagweight = TString::Format("SFWeight.BTagCSV40ge%d",abs(nbjets));

    TString selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), myCuts.Data());

    TString variable = TString::Format("%s>>+%s",plotvariable.Data(),hG     ->GetName());
    else { cout << "WTF " << fSamples[i].name << endl; continue;}
    if(plotvariable=="VPt"){
       TString tempvar = "GenPhoton[0].Pt()";
       variable = TString::Format("%s>>+%s",tempvar.Data(),hG     ->GetName());
    }
    gROOT->cd();
    if(fVerbose>2) cout << "  +++++++ Drawing      " << variable  << endl
                        << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
    nev = fSamples[i].tree->Draw(variable.Data(),selection.Data(),"goff");//fills histograms
    if(fVerbose>2) cout << "  +++++++ MC   events : "  <<  nev << endl;

}//for(size_t i = 0; i < fSamples.size(); ++i) //gamma

//add over/underflow
if(addunderflow) {
   hZ->SetBinContent(1, hZ->GetBinContent(0) + hZ->GetBinContent(1));
   hZ->SetBinError(  1, sqrt(pow(hZ->GetBinError(0),2) +  pow(hZ->GetBinError(1),2)));
   hG->SetBinContent(1, hG->GetBinContent(0) + hG->GetBinContent(1));
   hG->SetBinError(  1, sqrt(pow(hG->GetBinError(0),2) +  pow(hG->GetBinError(1),2)));
 } if(addoverflow){
   int nbinsX = hZ->GetNbinsX();
   hZ->SetBinContent(nbinsX, hZ->GetBinContent(nbinsX  ) + hZ->GetBinContent(nbinsX+1) );
   hZ->SetBinError(nbinsX, sqrt(pow(hZ->GetBinError(nbinsX  ),2) + pow(hZ->GetBinError(nbinsX+1),2) ) );
   hG->SetBinContent(nbinsX, hG->GetBinContent(nbinsX  ) + hG->GetBinContent(nbinsX+1) );
   hG->SetBinError(nbinsX, sqrt(pow(hG->GetBinError(nbinsX  ),2) + pow(hG->GetBinError(nbinsX+1),2) ) );
}

//make ratio and save histograms
hZG_MCratio->Divide(hZ,hG);

TFile *fnew = new TFile(outputdir+varname+".root","RECREATE");
fnew->cd();
hZ            ->Write();
hG            ->Write();
hZG_MCratio   ->Write();
fnew->Close();
cout << "saved histos in " << fnew->GetName() << endl;

}

//standard load function to load the MT2trees from the samples.dat
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
