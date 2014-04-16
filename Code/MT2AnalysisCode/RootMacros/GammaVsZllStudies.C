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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path of MT2tree.hh
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"

using namespace std;

//call this macro via root -l -b -q GammaVsZllStudies.C++
void GammaVsZllStudies();
void load(const char* filename = "samples/datasamples/samples_2141_dataonly.dat");

//standard sample struct that combines the MT2trees with the nessecary information of the tree, like sample type, cross-section, etc.
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

//this code only stores histograms of
//dileptonic Z data and MC (and the MC background for the dileptonic selection)
//single photon data and MC (and the MC background for the photon selection)
//In the beginning, the important cut variables are defined
//the variable you want to plot (including number of bins, and histogram borders)
//the variables possible are boson-pT, MT2, HT, NJets, NBJets)
//Furthermore you need to set HT-selection (defines datasample used), and lepton selection
void GammaVsZllStudies(){
gROOT->cd();
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

//apply weights for high HT: G: 0.956509, G_bg: 1.407555 - use this if using samples_1g_GEst_HT.dat, not for samples_1g_GEst_HTdata.dat // see l.392/393

//cut selection
// -1 means no cut, except for
//njets,nbjets, there it is -10
// -X for njets,nbjets means njets>=X
//while X means njets==X
double minphotonpt = 20;
double minzpt      = 20;
double maxphotonpt = -1;
double maxzpt      = -1;
double minHT       = 450;
double maxHT       = -1;
double minMT2      = -1;
double maxMT2      = -1;
double METcutval   = -1;
int    njets       = -2;
int    nbjets      = -10;

//plotvariable is the name within the MT2trees, they are
//misc.MT2, misc.HT, NBJets40CSVM, NJetsIDLoose40, VPt
//varname are the names the histogram is stored with:
//MT2, HT, NBJets, NJets, VPt
//varlabel is the label of the x-axis
//M_{T2} [GeV], H_{T} [GeV], number of jets, number of b jets, boson-p_{T} [GeV]
TString plotvariable = "VPt";//choose "VPt" to plot all VPts!
TString varname      = "VPt";//VPt // this is for saving the histograms!
TString varlabel     = "boson-p_{T} [GeV]";//boson-p_{T} [GeV]
int    nbins       = 30;
double binlow      = 0;
double binup       = 750;
//6 histograms: Z, Z_bg, G, G_bg, Z_data, G_data
//values for nbins, binlow, binup
//VPt or MT2: 30, 0, 750;   HT: 20, 450,1450;   NJets: 8, 2, 10;   NBJets: 5, 0, 5

int    HTselection = 1;//1:allHT - only for dileptons possible, for photons use PhotonHad trigger, 2: highPt (VPt>180 GeV due to trigger), 3: highHT (HT>750 GeV due to trigger)
int    LLselection = 1;//loose selection on leptons mean negative number, otherwise positive (default, due to dileptonic trigger)
                       //1:  ll(ee+mumu), 2: ee, 3: mumu, 4: emu, 5: all(ee+mumu+emu)
bool   addoverflow = true;//default = true
bool   addunderflow= false;//keep that false, default = false
bool   removebkg   = false;//keep that false, can do it later, default = false

TString samples_2l; TString samples_1g;//files to load for photon or dilepton selection
if(HTselection==2)             samples_2l = "samples/samples_2l_noHTdata.dat"; samples_1g = "samples/samples_1g_GEst_noHTdata.dat";
if(HTselection==3)             samples_2l = "samples/samples_2l_highHT.dat";   samples_1g = "samples/samples_1g_GEst_HTdata.dat";
if(HTselection==1)             samples_2l = "samples/samples_2l_noHTdata.dat"; samples_1g = "samples/samples_1g_GEst_PhotonHad.dat";//VPt>80 for g-data // VPt>125 for Z-MC
if(HTselection==2&&maxHT==750) samples_2l = "samples/samples_2l_lowHT.dat";    samples_1g = "samples/samples_1g_GEst_MET.dat";

TString outputdir  = "GammaVsZllStudies/";//if this directory does not exist, you need to create it before.

//Cuts should be self-explanatory
TString preselectioncuts = "NJetsIDLoose40>=2&&misc.HT>450&&NTausIDLoose3Hits==0&&misc.Vectorsumpt<70";
TString filters = "misc.PassJet40ID ==1&&(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)&&misc.CSCTightHaloIDFlag == 0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0";//&&(type1pfmet[0].Pt()>30||type1pfmet[0].Pt()/misc.CaloMETRaw<=2.)";
TString METcut           = TString::Format("misc.MET>=%d", (int)METcutval);
TString MinDPhi_g        = "misc.MinMetJetDPhi4Pt40>0.3";
TString MinDPhi_Z        = "MinMetJetDPhi(1,40.,2.4,1,4)>0.3";
TString eecutsloose      = "NEles==2&&NMuons==0";
TString eecuts           = "NEles==2&&NMuons==0&&ele[0].lv.Pt()>20&&ele[1].lv.Pt()>20&&ele[0].IDLoose&&ele[1].IDLoose";
TString mumucutsloose    = "NEles==0&&NMuons==2";
TString mumucuts         = "NEles==0&&NMuons==2&&muo[0].lv.Pt()>20&&muo[1].lv.Pt()>20";
TString llcutsloose      = "(("+eecutsloose+")||("+mumucutsloose+"))";
TString llcuts           = "(("+eecuts     +")||("+mumucuts     +"))";
TString emucutsloose     = "NEles==1&&NMuons==1";
TString emucuts          = "NEles==1&&NMuons==1&&muo[0].lv.Pt()>20&&ele[0].lv.Pt()>20&&ele[0].IDLoose";
TString gselloose        = "NEles==0&&NMuons==0&&NPhotons==1&&photon[0].isLooseID";
TString gselvar          = TString::Format("NEles==0&&NMuons==0&&NPhotons==1&&photon[0].isLooseID&&photon[0].lv.Pt()>=%d",(int)minphotonpt);
TString upper_gcutvar    = TString::Format("photon[0].lv.Pt()<%d",(int)maxphotonpt);
TString Zsel             = "GetDiLeptonPt(0,1,0,10,76,106)>=0";
TString Zselvar          = TString::Format("GetDiLeptonPt(0,1,0,10,76,106)>=%d", (int)minzpt);
TString upper_Zcutvar    = TString::Format("GetDiLeptonPt(0,1,0,10,76,106)<%d", (int)maxzpt);
TString Zemusel          = "GetDiLeptonPt(0,0,0,10,76,106)>=0";
TString Zemuselvar       = TString::Format("GetDiLeptonPt(0,0,0,10,76,106)>=%d", (int)minzpt);
TString upper_Zemucutvar = TString::Format("GetDiLeptonPt(0,0,0,10,76,106)<%d", (int)maxzpt);
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
TString ZRecalculate     = "ZllRecalculate()";//this means that dileptons are removed from event and added to MET; MT2, etc. get recalculated

TString triggerHT        = "( (trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1 || trigger.HLT_PFNoPUHT650_v4 == 1))";
TString triggerEE        = "(trigger.HLT_DiElectrons==1)";
TString triggerMM        = "(trigger.HLT_DiMuons==1)";
TString triggerEM        = "(trigger.HLT_EMu==1)";
TString triggerG         = "(trigger.HLT_SinglePhotons==1)";
TString triggerGHad      = "(trigger.HLT_SinglePhoton70_HT400==1)";

//change varname according to cuts, this will be the name of the TFile
                        varname += "_";
if(LLselection<0)       varname += "loose";
if(abs(LLselection)==1) varname += "SF";
if(abs(LLselection)==2) varname += "EE";
if(abs(LLselection)==3) varname += "MM";
if(abs(LLselection)==4) varname += "OF";
if(abs(LLselection)==5) varname += "SaOF";//same&opp.flavour
                        varname += "_";
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
//make ratio data/MC, make ratio G/Z <-- actually not important, will be done later in GammaVsZllStudiesRatios.C
TString histoname = varname + "_Z";
TH1D *hZ            = new TH1D(histoname.Data(), "", nbins, binlow, binup); hZ              ->Sumw2();
TString ytitle = TString::Format("Events / %d GeV", (int)hZ->GetBinWidth(1));
  hZ           ->GetXaxis()->SetTitle(varlabel.Data()); hZ           ->GetYaxis()->SetTitle(ytitle.Data());
  hZ           ->SetLineColor(kOrange  ); hZ           ->SetMarkerColor(kOrange  );   hZ           ->SetMarkerStyle(22);
histoname = varname + "_Zbg";
TH1D *hZ_bg         = new TH1D(histoname.Data(), "", nbins, binlow, binup); hZ_bg        ->Sumw2();
  hZ_bg        ->GetXaxis()->SetTitle(varlabel.Data()); hZ_bg        ->GetYaxis()->SetTitle(ytitle.Data());
  hZ_bg        ->SetLineColor(kOrange-3); hZ_bg        ->SetMarkerColor(kOrange-3);   hZ_bg        ->SetMarkerStyle(22);
histoname = varname + "_G";
TH1D *hG            = new TH1D(histoname.Data(), "", nbins, binlow, binup); hG           ->Sumw2();
  hG           ->GetXaxis()->SetTitle(varlabel.Data()); hG           ->GetYaxis()->SetTitle(ytitle.Data());
  hG           ->SetLineColor(kViolet-3); hG           ->SetMarkerColor(kViolet-3);   hG           ->SetMarkerStyle(20);
histoname = varname + "_Gbg";
TH1D *hG_bg         = new TH1D(histoname.Data(), "", nbins, binlow, binup); hG_bg        ->Sumw2();
  hG_bg        ->GetXaxis()->SetTitle(varlabel.Data()); hG_bg        ->GetYaxis()->SetTitle(ytitle.Data());
  hG_bg        ->SetLineColor(kPink    ); hG_bg        ->SetMarkerColor(kPink    );   hG_bg        ->SetMarkerStyle(20);
histoname = varname + "_Zdata";
TH1D *hZ_data       = new TH1D(histoname.Data(), "", nbins, binlow, binup); hZ_data      ->Sumw2();
  hZ_data      ->GetXaxis()->SetTitle(varlabel.Data()); hZ_data      ->GetYaxis()->SetTitle(ytitle.Data());  hZ_data      ->SetMarkerStyle(26);
histoname = varname + "_Gdata";
TH1D *hG_data       = new TH1D(histoname.Data(), "", nbins, binlow, binup); hG_data      ->Sumw2();
  hG_data      ->GetXaxis()->SetTitle(varlabel.Data()); hG_data      ->GetYaxis()->SetTitle(ytitle.Data());  hG_data      ->SetMarkerStyle(24);
histoname = varname + "_ZG_MCratio";
TH1D *hZG_MCratio   = new TH1D(histoname.Data(), "", nbins, binlow, binup); hZG_MCratio  ->Sumw2();
ytitle = TString::Format("#frac{Z(ll)}{#gamma} / %d GeV", (int)hZ->GetBinWidth(1));
  hZG_MCratio  ->GetXaxis()->SetTitle(varlabel.Data()); hZG_MCratio  ->GetYaxis()->SetTitle(ytitle.Data());
  hZG_MCratio  ->SetLineColor(kViolet-3); hZG_MCratio  ->SetMarkerColor(kViolet-3);   hZG_MCratio  ->SetMarkerStyle(20);
histoname = varname + "_ZG_Dataratio";
TH1D *hZG_Dataratio = new TH1D(histoname.Data(), "", nbins, binlow, binup); hZG_Dataratio->Sumw2();
  hZG_Dataratio->GetXaxis()->SetTitle(varlabel.Data()); hZG_Dataratio->GetYaxis()->SetTitle(ytitle.Data());  hZG_Dataratio->SetMarkerStyle(26);
histoname = varname + "_Z_DataVsMC";
TH1D *hZ_DataVsMC   = new TH1D(histoname.Data(), "", nbins, binlow, binup); hZ_DataVsMC  ->Sumw2();
ytitle = TString::Format("#frac{data}{MC} / %d GeV", (int)hZ->GetBinWidth(1));
  hZ_DataVsMC  ->GetXaxis()->SetTitle(varlabel.Data()); hZ_DataVsMC  ->GetYaxis()->SetTitle(ytitle.Data());
  hZ_DataVsMC  ->SetLineColor(kOrange  ); hZ_DataVsMC  ->SetMarkerColor(kOrange  );   hZ_DataVsMC  ->SetMarkerStyle(22);
histoname = varname + "_G_DataVsMC";
TH1D *hG_DataVsMC   = new TH1D(histoname.Data(), "", nbins, binlow, binup); hG_DataVsMC  ->Sumw2();
 hG_DataVsMC   ->GetXaxis()->SetTitle(varlabel.Data()); hG_DataVsMC  ->GetYaxis()->SetTitle(ytitle.Data());
 hG_DataVsMC   ->SetLineColor(kViolet-3); hG_DataVsMC  ->SetMarkerColor(kViolet-3);   hG_DataVsMC  ->SetMarkerStyle(20);

//define cuts depending on sample /also Z vs ZhighPt
//make two for-loops, first for dileptons, then for photons
load(samples_2l.Data());
for(size_t i = 0; i < fSamples.size(); ++i){ //ZLL
   string sampletype = "";
   //use only samples according to HT, VPt selection
   if(fSamples[i].sname=="ZJetsToLL_PtZ100" && minzpt< 180) continue;
   if(fSamples[i].sname=="ZJetsToLL"        && minzpt>=180) continue;
   if(HTselection==3 && fSamples[i].sname=="EE-Data"   ) continue;
   if(HTselection==3 && fSamples[i].sname=="EMu-Data"  ) continue;
   if(HTselection==3 && fSamples[i].sname=="MuMu-Data" ) continue;
   if(HTselection!=3 && fSamples[i].sname=="HT-Data_2l") continue;
   if(abs(LLselection)==1 && fSamples[i].sname=="EMu-Data" ) continue;
   if(abs(LLselection)==2 && fSamples[i].sname=="EMu-Data" ) continue;
   if(abs(LLselection)==3 && fSamples[i].sname=="EMu-Data" ) continue;
   if(abs(LLselection)==2 && fSamples[i].sname=="MuMu-Data") continue;
   if(abs(LLselection)==4 && fSamples[i].sname=="MuMu-Data") continue;
   if(abs(LLselection)==3 && fSamples[i].sname=="EE-Data"  ) continue;
   if(abs(LLselection)==4 && fSamples[i].sname=="EE-Data"  ) continue;

   //assign samples
   if(fSamples[i].sname=="ZJetsToLL_PtZ100") sampletype = "Z";
   else if(fSamples[i].sname=="ZJetsToLL")   sampletype = "Z";
   else if(fSamples[i].sname=="Top")         sampletype = "Z_bg";
   else if(fSamples[i].sname=="VV")          sampletype = "Z_bg";
   else if(fSamples[i].sname=="Wtolnu")      sampletype = "Z_bg";
   else if(fSamples[i].sname=="EE-Data")     sampletype = "Z_data";
   else if(fSamples[i].sname=="MuMu-Data")   sampletype = "Z_data";
   else if(fSamples[i].sname=="EMu-Data")    sampletype = "Z_data";
   else if(fSamples[i].sname=="HT-Data")     sampletype = "Z_data";
   else if(fSamples[i].sname=="ZJetsToNuNu") sampletype = "Z_inv";
   else { cout << "sample " << fSamples[i].name << " is not considered" << endl; continue; }

   //load MT2trees
    MT2tree* fMT2tree = new MT2tree();
    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    int nev =0;

    //these lines define the cuts to select only the events within previously defined selection, a bit cumbersome
    TString myCuts = preselectioncuts + "&&" + filters + "&&" + NJetsCut + "&&" + NBJetsCut + "&&" + ZRecalculate;
    if((int)minHT >0) myCuts = myCuts+ "&&" + HTlow;
    if((int)maxHT >0 && maxHT >minHT ) myCuts = myCuts+ "&&" + HTup;
    if(HTselection==3 && minHT<750.  ) myCuts = myCuts+ "&&" + "misc.HT>=750";
    if((int)minMT2>0) myCuts = myCuts+ "&&" + MT2low;
    if((int)maxMT2>0 && maxMT2>minMT2) myCuts = myCuts+ "&&" + MT2up;
    if(abs(LLselection)==5){
       if(LLselection<0) myCuts = myCuts+ "&&("+llcutsloose + "||(" + emucutsloose + "))";
       else              myCuts = myCuts+ "&&("+llcuts      + "||(" + emucuts      + "))";
    } else if(abs(LLselection)==4){
       if(LLselection<0) myCuts = myCuts+ "&&" + emucutsloose;
       else              myCuts = myCuts+ "&&" + emucuts;
    } else if(abs(LLselection)==3){
       if(LLselection<0) myCuts = myCuts+ "&&" + mumucutsloose;
       else              myCuts = myCuts+ "&&" + mumucuts;
    } else if(abs(LLselection)==2){
       if(LLselection<0) myCuts = myCuts+ "&&" + eecutsloose;
       else              myCuts = myCuts+ "&&" + eecuts;
    } else if(abs(LLselection)==1){
       if(LLselection<0) myCuts = myCuts+ "&&" + llcutsloose;
       else              myCuts = myCuts+ "&&" + llcuts;
    }
    if(abs(LLselection)==5){
       if(minzpt>0)                myCuts = myCuts+ "&&(" + Zselvar + "||" + Zemuselvar + ")";
       else                        myCuts = myCuts+ "&&(" + Zsel    + "||" + Zemusel    + ")";
       if(maxzpt>0&&maxzpt>minzpt) myCuts = myCuts+ "&&(" + upper_Zcutvar + "||" + upper_Zemucutvar + ")";
    } else if(abs(LLselection)==4){
       if(minzpt>0)                myCuts = myCuts+ "&&" + Zemuselvar;
       else                        myCuts = myCuts+ "&&" + Zemusel;
       if(maxzpt>0&&maxzpt>minzpt) myCuts = myCuts+ "&&" + upper_Zemucutvar;
    } else {
       if(minzpt>0)                myCuts = myCuts+ "&&" + Zselvar;
       else                        myCuts = myCuts+ "&&" + Zsel;
       if(maxzpt>0&&maxzpt>minzpt) myCuts = myCuts+ "&&" + upper_Zcutvar;
    }

    if(fSamples[i].sname=="EE-Data")     myCuts = myCuts+ "&&" + eecuts   + "&&" + triggerEE;
    if(fSamples[i].sname=="MuMu-Data")   myCuts = myCuts+ "&&" + mumucuts + "&&" + triggerMM;
    if(fSamples[i].sname=="EMu-Data")    myCuts = myCuts+ "&&" + emucuts  + "&&" + triggerEM;
    if(fSamples[i].sname=="HT-Data_2l")  myCuts = myCuts+ "&&" + "misc.HT>=750" + "&&" + triggerHT;


   //global event weight
    double weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
    if(fVerbose>2) cout << "===========================================================" <<  endl;
    if(fVerbose>2) cout << "MakePlot: looping over " << fSamples[i].name << "-----------------------------------" <<  endl;
    if(fVerbose>2) cout << "  +++++++ weight:      " << weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 

    TString btagweight = "1.00"; //stored btag weights up to >=3, default is weight==1 to avoid non-existing weights
    //load correct BTV SF weight
    if(nbjets>=0 && nbjets<=3) btagweight = TString::Format("SFWeight.BTagCSV40eq%d",abs(nbjets));
    else if(nbjets>=-3)        btagweight = TString::Format("SFWeight.BTagCSV40ge%d",abs(nbjets));

    //Define selection and variable as it is used in the tree->Draw() function
    TString selection;
    if(    fSamples[i].type!="data") selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), myCuts.Data());
    else                             selection = TString::Format("(%.15f) * (%s)",                 weight,                    myCuts.Data()); 

    TString variable;
    if(     sampletype=="Z"     ) variable = TString::Format("%s>>+%s",plotvariable.Data(),hZ     ->GetName());
    else if(sampletype=="Z_bg"  ) variable = TString::Format("%s>>+%s",plotvariable.Data(),hZ_bg  ->GetName());
    else if(sampletype=="Z_data") variable = TString::Format("%s>>+%s",plotvariable.Data(),hZ_data->GetName());
    else { cout << "WTF " << fSamples[i].name << endl; continue;}
    if(plotvariable=="VPt"){
       TString tempvar = "";
       if(abs(LLselection)==4) tempvar = "GetDiLeptonPt(0,0,0,10,76,106)";
       else                    tempvar = "GetDiLeptonPt(0,1,0,10,76,106)";
       if(     sampletype=="Z"     ) variable = TString::Format("%s>>+%s",tempvar.Data(),hZ     ->GetName());
       else if(sampletype=="Z_bg"  ) variable = TString::Format("%s>>+%s",tempvar.Data(),hZ_bg  ->GetName());
       else if(sampletype=="Z_data") variable = TString::Format("%s>>+%s",tempvar.Data(),hZ_data->GetName());
    }
    gROOT->cd();//needed, before I had problems that the histograms were not filled.
    if(fVerbose>2) cout << "  +++++++ Drawing      " << variable  << endl
                        << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
    nev = fSamples[i].tree->Draw(variable.Data(),selection.Data(),"goff");//this fills the histograms
    if(fVerbose>2) cout << "  +++++++ MC   events : "  <<  nev << endl;

    if(plotvariable=="VPt" && abs(LLselection)==5){//did draw already ll, now emu
       TString tempvar = "GetDiLeptonPt(0,0,0,10,76,106)";
       if(     sampletype=="Z"     ) variable = TString::Format("%s>>+%s",tempvar.Data(),hZ     ->GetName());
       else if(sampletype=="Z_bg"  ) variable = TString::Format("%s>>+%s",tempvar.Data(),hZ_bg  ->GetName());
       else if(sampletype=="Z_data") variable = TString::Format("%s>>+%s",tempvar.Data(),hZ_data->GetName());
       if(fVerbose>2) cout << "  +++++++ Drawing      " << variable  << endl
                           << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
       nev = fSamples[i].tree->Draw(variable.Data(),selection.Data(),"goff");
       if(fVerbose>2) cout << "  +++++++ MC   events : "  <<  nev << endl;
    }

}//for(size_t i = 0; i < fSamples.size(); ++i) //ZLL

//this for loop is equivalent to the previous one but for photon selection
load(samples_1g.Data());
gROOT->cd();
for(size_t i = 0; i < fSamples.size(); ++i){ //gamma
   string sampletype = "";
   if(HTselection==3 && fSamples[i].sname=="G-Data"    ) continue;
   if(HTselection!=3 && fSamples[i].sname=="HT-Data"   ) continue;
   if(fSamples[i].sname=="Photons")          sampletype = "G";
   else if(fSamples[i].sname=="QCD")         sampletype = "G_bg";
   else if(fSamples[i].sname=="G-Data")      sampletype = "G_data";
   else if(fSamples[i].sname=="HT-Data")     sampletype = "G_data";
   else { cout << "sample " << fSamples[i].name << " is not considered" << endl; continue; }

    MT2tree* fMT2tree = new MT2tree();
    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    int nev =0;
    TString myCuts = preselectioncuts + "&&" + filters + "&&" + NJetsCut + "&&" + NBJetsCut;
//optional
//    myCuts = myCuts + "&&" + METcut + "&&" + MinDPhi_g;
    if((int)minHT >0) myCuts = myCuts+ "&&" + HTlow;
    if((int)maxHT >0 && maxHT >minHT ) myCuts = myCuts+ "&&" + HTup;
    if(HTselection==3 && minHT<750.  ) myCuts = myCuts+ "&&" + "misc.HT>=750";
    if((int)minMT2>0) myCuts = myCuts+ "&&" + MT2low;
    if((int)maxMT2>0 && maxMT2>minMT2) myCuts = myCuts+ "&&" + MT2up;
    if(minphotonpt>0) myCuts = myCuts+ "&&" + gselvar;
    else              myCuts = myCuts+ "&&" + gselloose;
    if(maxphotonpt>0&&maxphotonpt>minphotonpt) myCuts = myCuts+ "&&" + upper_gcutvar;

    if(fSamples[i].sname=="G-Data" && fSamples[i].name.Contains("PhotonHad"))    myCuts = myCuts+ "&&" + "photon[0].lv.Pt()>80"  + "&&" + triggerGHad;
    else if(fSamples[i].sname=="G-Data")                                         myCuts = myCuts+ "&&" + "photon[0].lv.Pt()>180" + "&&" + triggerG;
    if(fSamples[i].sname=="HT-Data_1g")                                          myCuts = myCuts+ "&&" + "misc.HT>=750" + "&&" + triggerHT;

    double weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
    if(fVerbose>2) cout << "===========================================================" <<  endl;
    if(fVerbose>2) cout << "MakePlot: looping over " << fSamples[i].name << "-----------------------------------" <<  endl;
    if(fVerbose>2) cout << "  +++++++ weight:      " << weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 

    TString btagweight = "1.00"; //stored btag weights up to >=3, default is weight==1 to avoid non-existing weights
    if(nbjets>=0 && nbjets<=3) btagweight = TString::Format("SFWeight.BTagCSV40eq%d",abs(nbjets));
    else if(nbjets>=-3)        btagweight = TString::Format("SFWeight.BTagCSV40ge%d",abs(nbjets));

    TString selection;
    if(    fSamples[i].type!="data") selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), myCuts.Data());
    else                             selection = TString::Format("(%.15f) * (%s)",                 weight,                    myCuts.Data()); 

    TString variable;
    if(     sampletype=="G"     ) variable = TString::Format("%s>>+%s",plotvariable.Data(),hG     ->GetName());
    else if(sampletype=="G_bg"  ) variable = TString::Format("%s>>+%s",plotvariable.Data(),hG_bg  ->GetName());
    else if(sampletype=="G_data") variable = TString::Format("%s>>+%s",plotvariable.Data(),hG_data->GetName());
    else { cout << "WTF " << fSamples[i].name << endl; continue;}
    if(plotvariable=="VPt"){
       TString tempvar = "photon[0].lv.Pt()";
       if(     sampletype=="G"     ) variable = TString::Format("%s>>+%s",tempvar.Data(),hG     ->GetName());
       else if(sampletype=="G_bg"  ) variable = TString::Format("%s>>+%s",tempvar.Data(),hG_bg  ->GetName());
       else if(sampletype=="G_data") variable = TString::Format("%s>>+%s",tempvar.Data(),hG_data->GetName());
    }
    gROOT->cd();
    if(fVerbose>2) cout << "  +++++++ Drawing      " << variable  << endl
                        << "  +++++++ with cuts:   " << setw(40)  << selection << endl;
    nev = fSamples[i].tree->Draw(variable.Data(),selection.Data(),"goff");
    if(fVerbose>2) cout << "  +++++++ MC   events : "  <<  nev << endl;

}//for(size_t i = 0; i < fSamples.size(); ++i) //gamma

//add over/underflow
if(addunderflow) {
   hZ->SetBinContent(1, hZ->GetBinContent(0) + hZ->GetBinContent(1));
   hZ->SetBinError(  1, sqrt(pow(hZ->GetBinError(0),2) +  pow(hZ->GetBinError(1),2)));
   hZ_bg->SetBinContent(1, hZ_bg->GetBinContent(0) + hZ_bg->GetBinContent(1));
   hZ_bg->SetBinError(  1, sqrt(pow(hZ_bg->GetBinError(0),2) +  pow(hZ_bg->GetBinError(1),2)));
   hZ_data->SetBinContent(1, hZ_data->GetBinContent(0) + hZ_data->GetBinContent(1));
   hZ_data->SetBinError(  1, sqrt(pow(hZ_data->GetBinError(0),2) +  pow(hZ_data->GetBinError(1),2)));

   hG->SetBinContent(1, hG->GetBinContent(0) + hG->GetBinContent(1));
   hG->SetBinError(  1, sqrt(pow(hG->GetBinError(0),2) +  pow(hG->GetBinError(1),2)));
   hG_bg->SetBinContent(1, hG_bg->GetBinContent(0) + hG_bg->GetBinContent(1));
   hG_bg->SetBinError(  1, sqrt(pow(hG_bg->GetBinError(0),2) +  pow(hG_bg->GetBinError(1),2)));
   hG_data->SetBinContent(1, hG_data->GetBinContent(0) + hG_data->GetBinContent(1));
   hG_data->SetBinError(  1, sqrt(pow(hG_data->GetBinError(0),2) +  pow(hG_data->GetBinError(1),2)));
} if(addoverflow){
   int nbinsX = hZ->GetNbinsX();
   hZ->SetBinContent(nbinsX, hZ->GetBinContent(nbinsX  ) + hZ->GetBinContent(nbinsX+1) );
   hZ->SetBinError(nbinsX, sqrt(pow(hZ->GetBinError(nbinsX  ),2) + pow(hZ->GetBinError(nbinsX+1),2) ) );
   hZ_bg->SetBinContent(nbinsX, hZ_bg->GetBinContent(nbinsX  ) + hZ_bg->GetBinContent(nbinsX+1) );
   hZ_bg->SetBinError(nbinsX, sqrt(pow(hZ_bg->GetBinError(nbinsX  ),2) + pow(hZ_bg->GetBinError(nbinsX+1),2) ) );
   hZ_data->SetBinContent(nbinsX, hZ_data->GetBinContent(nbinsX  ) + hZ_data->GetBinContent(nbinsX+1) );
   hZ_data->SetBinError(nbinsX, sqrt(pow(hZ_data->GetBinError(nbinsX  ),2) + pow(hZ_data->GetBinError(nbinsX+1),2) ) );
   hG->SetBinContent(nbinsX, hG->GetBinContent(nbinsX  ) + hG->GetBinContent(nbinsX+1) );
   hG->SetBinError(nbinsX, sqrt(pow(hG->GetBinError(nbinsX  ),2) + pow(hG->GetBinError(nbinsX+1),2) ) );
   hG_bg->SetBinContent(nbinsX, hG_bg->GetBinContent(nbinsX  ) + hG_bg->GetBinContent(nbinsX+1) );
   hG_bg->SetBinError(nbinsX, sqrt(pow(hG_bg->GetBinError(nbinsX  ),2) + pow(hG_bg->GetBinError(nbinsX+1),2) ) );
   hG_data->SetBinContent(nbinsX, hG_data->GetBinContent(nbinsX  ) + hG_data->GetBinContent(nbinsX+1) );
   hG_data->SetBinError(nbinsX, sqrt(pow(hG_data->GetBinError(nbinsX  ),2) + pow(hG_data->GetBinError(nbinsX+1),2) ) );
}

//make ratios - this is not needed, because the ratio analysis is done in a separate macro
//called GammaVsZllStudiesRatios.C
if(removebkg){
    hZG_MCratio->Divide(hZ,hG);
  //optional
   TH1D *tempZ = (TH1D*)hZ_data->Clone("cloneZ");
   tempZ->Add(hZ_bg,-1.);
   TH1D *tempG = (TH1D*)hG_data->Clone("cloneG");
   tempG->Add(hG_bg,-1.);
   hZG_Dataratio->Divide(tempZ,tempG);
} else {
   hZG_MCratio->Divide(hZ,hG);
   hZG_Dataratio->Divide(hZ_data,hG_data);
} if(removebkg){ 
   //optional
   TH1D *tempZ = (TH1D*)hZ_data->Clone("cloneZ");
   tempZ->Add(hZ_bg,-1.);
   hZ_DataVsMC->Divide(tempZ,hZ);
} else {
   hZ_DataVsMC->Divide(hZ_data,hZ);
} if(removebkg){
   //optional
   TH1D *tempG = (TH1D*)hG_data->Clone("cloneG");
   tempG->Add(hG_bg,-1.);
   hG_DataVsMC->Divide(tempG,hG);
} else {
   hG_DataVsMC->Divide(hG_data,hG);
}

//store the histograms in a default named TFile
TFile *fnew = new TFile(outputdir+varname+".root","RECREATE");
fnew->cd();
hZ            ->Write();
hZ_bg         ->Write();
hZ_data       ->Write();
hG            ->Write();
hG_bg         ->Write();
hG_data       ->Write();
hZG_MCratio   ->Write();
hZG_Dataratio ->Write();
hZ_DataVsMC   ->Write();
hG_DataVsMC   ->Write();
fnew->Close();
cout << "saved histos in " << fnew->GetName() << endl;


}

//standard function to load the MT2trees in the samples.dat
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
