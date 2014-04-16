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

//run via root -l -b -q makeCutFlow.C++
//note that this is a code from the 7 TeV analysis and has not been updated
//it is very clumsy and could be done better

using namespace std;

void load(const char* filename = "/shome/haweber/CMSSW_4_2_3/src/DiLeptonAnalysis/NTupleProducer/MT2Analysis/Code/MT2AnalysisCode/MT2Code/samples/datasamples/samples_2141_dataonly.dat");
void GetContent2DRange(TH2D* num, double &eff, double &err, double xlow, double xup, double ylow, double yup);
void GetContent2DRangeBin(TH2D* num, double &eff, double &err, int xlow, int xup, int ylow, int yup);
void PrintCutflowTableLine(double HTlow, double HTup, double MT2low, double MT2up, bool susy, map<string, TH2D*> hists );
void makeCutFlow();

//struct that combines MT2trees with necessary information like cross section
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


// USER INPUT -------------------------------------------
Bool_t  lowMT2                    = true; //MT2b analysis selection
Bool_t  highMT2                   = false;//MT2 analysis selection
Bool_t  lowHT                     = false;//low HT selection
Bool_t  highHT                    = true; //high HT selection
Bool_t  dofastestimate            = false;//do an MT2 cut from the beginning, otherwise first line in table is without MT2 cut
TString samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/_4400_noMinDPhi.dat";
TString outputdir= "CutFlows/";//default
TString outputname = "CutflowHistograms.root";//default
// --------------------------------------------------------

//produces a cutflow table (along MT2) for data and simulation, and possibly a SUSY signal
void makeCutFlow(){

    //safe histograms used for the cutflow table
    if(highMT2)               outputname = "CutflowHistograms_MT2.root";
    if(lowMT2)                outputname = "CutflowHistograms_MT2b.root";
    if(highMT2 && highHT)     outputname = "CutflowHistograms_highHT_MT2.root";
    else if(lowMT2 && highHT) outputname = "CutflowHistograms_highHT_MT2b.root";
    else if(highMT2 && lowHT) outputname = "CutflowHistograms_lowHT_MT2.root";
    else if(lowMT2 && lowHT)  outputname = "CutflowHistograms_lowHT_MT2b.root";

    //hard-coded samples.dat
    if(lowHT)  samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4400_noMinDPhi.dat";
    if(highHT) samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4600_noMinDPhi.dat";
    if(highMT2 && lowHT)  samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4400_noMinDPhi_MT2_LM.dat";
    if(highMT2 && highHT) samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4600_noMinDPhi_MT2_LM.dat";
    if(lowMT2 && lowHT)  samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4400_noMinDPhi_MT2b_LM.dat";
    if(lowMT2 && highHT) samples = "/shome/haweber/MT2Analysis/Code/MT2AnalysisCode/RootMacros/samples/samples_4600_noMinDPhi_MT2b_LM.dat";
    
    //load histograms for output
    map<string, TH2D*> histos;

    const int NHTbins  =   5;
    double  HTbins[NHTbins+1] = {400., 750., 950., 1500., 3000., 7000.};
    
    const int NMT2bins =  12;//must contain all bin borders
    double MT2bins[NMT2bins+1] = {0., 80., 125., 150., 180., 200., 260., 275., 300., 375., 500., 1000., 7000.};
    const int NMT2bins2 = 19;
    double MT2bins2[NMT2bins2+1] = {0.,10.,20.,30.,40.,50., 60., 70., 80.,100.,125.,150.,180.,200.,260.,275.,300.,375.,500.,1000.};
    const int NMT2bins3 = 12;
    double MT2bins3[NMT2bins3+1] = {0.,125.,150.,175.,200.,225.,250.,275.,300.,325.,350.,375.,500.};
    
    const int sampletypesize = 16;
    string sample_type[sampletypesize] = {"QCD", "WJets", "ZJets", "Top", "Other", "mc", "LM1", "LM2", "LM3", "LM4", "LM5", "LM6", "LM7", "LM8", "LM9", "data"};

	for(int is = 0; is<sampletypesize; ++is){
		string hs = string("_") + sample_type[is];
		string mapname;
        	mapname = "Events" + hs;
        	if(histos.count(mapname) == 0 ) histos[mapname] = new TH2D(mapname.c_str(), "eff", NMT2bins2, MT2bins2, NHTbins, HTbins);
	}
    
	for(map<string,TH2D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Sumw2();}
    
	cout << "Histograms loaded" << endl;

	//event selection
    std::ostringstream cutStream;
	cutStream << " " 
    << "misc.MET>=30"                                                  << "&&"
    << "misc.HT > 750 "                                                << "&&"//trigger cut
    << "NJetsIDLoose40 >=3"                                            << "&&"
    << "(NEles==0  || ele[0].lv.Pt()<10)"                              << "&&"
    << "(NMuons==0 || muo[0].lv.Pt()<10)"                              << "&&"
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
    if(dofastestimate){
        cutStream << "&&"
        << "misc.MT2>=125";//lowest MT2 for both MT2 and MT2b
    }
    if(lowMT2){
        cutStream << "&&"
        << "misc.LeadingJPt  >150"                                     << "&&"
        << "NJetsIDLoose40 >=4"                                        << "&&"
        << "NBJets >0"                                                 << "&&"
        << "(misc.MinMetJetDPhi >0.3||misc.MinMetJetDPhiIndex>3)";
	if(dofastestimate) cutStream << "&&" << "misc.MT2>=125";

    }
    if(highMT2){
        cutStream << "&&"
        << "NJetsIDLoose40 >=3"                                        << "&&"
        << "misc.MinMetJetDPhi >0.3";
	if(dofastestimate) cutStream << "&&" << "misc.MT2>=150";
    }
    if(highHT && dofastestimate){
        cutStream << "&&"
        << "misc.HT>=950";
    }
    if(lowHT && dofastestimate){
        cutStream << "&&"
        << "misc.HT<950&&misc.HT>=750";
    }

	TString cuts = cutStream.str().c_str();
   
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
	<< "(trigger.HLT_HT600_v1 ==1 && (misc.Run>=173236 && misc.Run< 178420))" << "||";
	if(lowHT || !(highHT)){// !(highHT) if both lowHT and highHT are false
	   triggerStream << "(trigger.HLT_HT650_v4 ==1 && (misc.Run>=178420 && misc.Run< 179959))" << " )";
        }
	if(highHT){
	   triggerStream << "(trigger.HLT_HT700_v2 ==1 && (misc.Run>=178420 && misc.Run< 179959))" << "||"
	                 << "(trigger.HLT_PFHT650_v1==1 && misc.Run>=179959)" << " )";
        }
	TString trigger = triggerStream.str().c_str();

	//load samples
	load(samples.Data());

  for(size_t i = 0; i < fSamples.size(); ++i){

     string sampletype = (string)fSamples[i].type;
     if(sampletype==(string)"mc"){
	if(fSamples[i].sname=="QCD") sampletype = (string)"QCD";
	else if(fSamples[i].sname=="Wtolnu") sampletype = (string)"WJets";
	else if(fSamples[i].sname=="DY") sampletype = (string)"ZJets";
	else if(fSamples[i].sname=="Top") sampletype = (string)"Top";//no ttbar, includes TTZ, TTW
	else sampletype = (string)"Other";
    }
    if(sampletype==(string)"susy") sampletype = (string)fSamples[i].sname;

    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents * fSamples[i].PU_avg_weight );
    if(fVerbose>2) cout << "PrintCutFlow: looping over " << fSamples[i].name << endl;
    if(fVerbose>2) cout << "              sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
    if(fVerbose>2) cout << "              Average PU weight: " << fSamples[i].PU_avg_weight << endl;
    if(fVerbose>2) cout << "              Original Entries: " << fSamples[i].nevents << endl;
    
    MT2tree *fMT2tree = new MT2tree();
    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    int nev =0;

    //leo tweak - filtering out the TTree
    TString myCuts = cuts;

    if( fSamples[i].type=="data") myCuts += " && " + trigger; //cuts to be aplied only on data
    cout << "Cuts for Flow: " << myCuts << endl;
    fSamples[i].tree->Draw(">>selList", myCuts);


    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
    fSamples[i].tree->SetEventList(myEvtList);
    int counter=0;
    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
    
    if(myEvtList->GetSize()==0) continue;
    //run over selected events - note that one could do this also using tree->Draw() function!
    while(myEvtList->GetEntry(counter++) !=-1){
      
      int jentry = myEvtList->GetEntry(counter-1);
      
      //for (Long64_t jentry=0; jentry<nentries;jentry++) {
      nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
      fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
      
      if ( fVerbose>2 && jentry % 50000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;
      
      Double_t weight = sample_weight;
      if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 

      histos["Events_" + sampletype]->Fill(fMT2tree->misc.MT2, fMT2tree->misc.HT, weight);//fill histogram
   }//while event
  }//for sample

  //sum up mc samples to mcsum histogram
  for(int ists = 0; ists<5; ++ists){//first five samples are mc
	if(sample_type[ists]==string("data") || sample_type[ists]==string("susy") || sample_type[ists]==string("mc") ) continue;//have all mc subsamples only
	string helptype = string("_") + sample_type[ists];
	histos[string("Events_mc")]->Add(histos[string("Events") + helptype], 1.);
  }

      TFile *fsavefile = new TFile(outputdir + outputname,"RECREATE");
	fsavefile->cd();
	for(map<string,TH2D*>::iterator h=histos.begin(); h!=histos.end();++h){
		if(h->second->GetEntries()>0) h->second->Write();
	}
	fsavefile->Close();
	cout << "Saved histograms in " << outputdir << outputname << endl;
    
	cout << endl << endl;


	//output numbers in table
	//first mc + data
	cout << "\%BEGINLATEX\%"             << endl;
	cout << "\\begin{table}"             << endl
	     << "\\begin{center}"            << endl
	//     << "\\small"                    << endl
             << "\\begin{tabular}{lccccccc}" << endl	     
	     << "\\hline\\hline"             << endl;
	cout << " Process                  " << " & " << "QCD" << " & " << "W+Jets" << " & " << "Z+Jets" << " & " << "Top" << " & " << "Other" << " & " << "Total Bgk." << " & " << "Data" << "\\\\" << endl;
	cout << "\\hline"                    << endl;
	if(highMT2 && !(lowHT || highHT)){//high MT2, both HT;
		PrintCutflowTableLine(750., 7000., 0.,   7000., false, histos);
		PrintCutflowTableLine(750., 7000., 150., 200.,  false, histos);
		PrintCutflowTableLine(750., 7000., 200., 275.,  false, histos);
		PrintCutflowTableLine(750., 7000., 275., 375.,  false, histos);
		PrintCutflowTableLine(750., 7000., 375., 500.,  false, histos);
		PrintCutflowTableLine(750., 7000., 500., 7000., false, histos);
	}
	if(highMT2 && lowHT){
		PrintCutflowTableLine(750., 950.,  0.,   7000., false, histos);
		PrintCutflowTableLine(750., 950.,  150., 200.,  false, histos);
		PrintCutflowTableLine(750., 950.,  200., 275.,  false, histos);
		PrintCutflowTableLine(750., 950.,  275., 375.,  false, histos);
		PrintCutflowTableLine(750., 950.,  375., 500.,  false, histos);
		PrintCutflowTableLine(750., 950.,  500., 7000., false, histos);
	}
	if(highMT2 && highHT){
		PrintCutflowTableLine(950., 7000., 0.,   7000., false, histos);
		PrintCutflowTableLine(950., 7000., 150., 200.,  false, histos);
		PrintCutflowTableLine(950., 7000., 200., 275.,  false, histos);
		PrintCutflowTableLine(950., 7000., 275., 375.,  false, histos);
		PrintCutflowTableLine(950., 7000., 375., 500.,  false, histos);
		PrintCutflowTableLine(950., 7000., 500., 7000., false, histos);
	}
	if(lowMT2 && !(lowHT || highHT)){//low MT2, both HT, use lowHT binning
		PrintCutflowTableLine(750., 7000., 0.,   7000., false, histos);
		PrintCutflowTableLine(750., 7000., 125., 150.,  false, histos);
		PrintCutflowTableLine(750., 7000., 150., 200.,  false, histos);
		PrintCutflowTableLine(750., 7000., 200., 300.,  false, histos);
		PrintCutflowTableLine(750., 7000., 300., 7000., false, histos);
	}
	if(lowMT2 && lowHT){
		PrintCutflowTableLine(750., 950.,  0.,   7000., false, histos);
		PrintCutflowTableLine(750., 950.,  125., 150.,  false, histos);
		PrintCutflowTableLine(750., 950.,  150., 200.,  false, histos);
		PrintCutflowTableLine(750., 950.,  200., 300.,  false, histos);
		PrintCutflowTableLine(750., 950.,  300., 7000., false, histos);
	}
	if(lowMT2 && highHT){
		PrintCutflowTableLine(950., 7000., 0.,   7000., false, histos);
		PrintCutflowTableLine(950., 7000., 125., 150.,  false, histos);
		PrintCutflowTableLine(950., 7000., 150., 180.,  false, histos);
		PrintCutflowTableLine(950., 7000., 180., 260.,  false, histos);
		PrintCutflowTableLine(950., 7000., 260., 7000., false, histos);
	}
	cout << "\\hline\\hline"             << endl
	     << "\\end{tabular}"             << endl
	     << "\\end{center}"              << endl
	     << "\\end{table}"               << endl
	     << "\%ENDLATEX\%"               << endl
	     << endl << endl;

	//then SUSY
	cout << "\%BEGINLATEX\%"               << endl;
	cout << "\\begin{table}"               << endl
	     << "\\begin{center}"              << endl
	//     << "\\small"                      << endl
             << "\\begin{tabular}{lccccccccc}" << endl	     
	     << "\\hline\\hline"               << endl;
	cout << " Process                  "   << " & " << "LM1" << " & " << "LM2" << " & " << "LM3" << " & " << "LM4" << " & " << "LM5" << " & " << "LM6" << " & " << "LM7" << " & " << "LM8" << " & " << "LM9" << "\\\\" << endl;
	cout << "\\hline"                      << endl;
	if(highMT2 && !(lowHT || highHT)){//high MT2, both HT;
		PrintCutflowTableLine(750., 7000., 0.,   7000., true, histos);
		PrintCutflowTableLine(750., 7000., 150., 200.,  true, histos);
		PrintCutflowTableLine(750., 7000., 200., 275.,  true, histos);
		PrintCutflowTableLine(750., 7000., 275., 375.,  true, histos);
		PrintCutflowTableLine(750., 7000., 375., 500.,  true, histos);
		PrintCutflowTableLine(750., 7000., 500., 7000., true, histos);
	}
	if(highMT2 && lowHT){
		PrintCutflowTableLine(750., 950.,  0.,   7000., true, histos);
		PrintCutflowTableLine(750., 950.,  150., 200.,  true, histos);
		PrintCutflowTableLine(750., 950.,  200., 275.,  true, histos);
		PrintCutflowTableLine(750., 950.,  275., 375.,  true, histos);
		PrintCutflowTableLine(750., 950.,  375., 500.,  true, histos);
		PrintCutflowTableLine(750., 950.,  500., 7000., true, histos);
	}
	if(highMT2 && highHT){
		PrintCutflowTableLine(950., 7000., 0.,   7000., true, histos);
		PrintCutflowTableLine(950., 7000., 150., 200.,  true, histos);
		PrintCutflowTableLine(950., 7000., 200., 275.,  true, histos);
		PrintCutflowTableLine(950., 7000., 275., 375.,  true, histos);
		PrintCutflowTableLine(950., 7000., 375., 500.,  true, histos);
		PrintCutflowTableLine(950., 7000., 500., 7000., true, histos);
	}
	if(lowMT2 && !(lowHT || highHT)){//low MT2, both HT, use lowHT binning
		PrintCutflowTableLine(750., 7000., 0.,   7000., true, histos);
		PrintCutflowTableLine(750., 7000., 125., 150.,  true, histos);
		PrintCutflowTableLine(750., 7000., 150., 200.,  true, histos);
		PrintCutflowTableLine(750., 7000., 200., 300.,  true, histos);
		PrintCutflowTableLine(750., 7000., 300., 7000., true, histos);
	}
	if(lowMT2 && lowHT){
		PrintCutflowTableLine(750., 950.,  0.,   7000., true, histos);
		PrintCutflowTableLine(750., 950.,  125., 150.,  true, histos);
		PrintCutflowTableLine(750., 950.,  150., 200.,  true, histos);
		PrintCutflowTableLine(750., 950.,  200., 300.,  true, histos);
		PrintCutflowTableLine(750., 950.,  300., 7000., true, histos);
	}
	if(lowMT2 && highHT){
		PrintCutflowTableLine(950., 7000., 0.,   7000., true, histos);
		PrintCutflowTableLine(950., 7000., 125., 150.,  true, histos);
		PrintCutflowTableLine(950., 7000., 150., 180.,  true, histos);
		PrintCutflowTableLine(950., 7000., 180., 260.,  true, histos);
		PrintCutflowTableLine(950., 7000., 260., 7000., true, histos);
	}
	cout << "\\hline\\hline"               << endl
	     << "\\end{tabular}"               << endl
	     << "\\end{center}"                << endl
	     << "\\end{table}"                 << endl
	     << "\%ENDLATEX\%"                 << endl
	     << endl << endl;

}//makeCutFlow

//this prints a single line of the cutflow table
void PrintCutflowTableLine(double HTlow, double HTup, double MT2low, double MT2up, bool susy, map<string, TH2D*> hists){

	double  qcd,  wjets,  zjets,  top,  other,  mc,  data,  lm1,  lm2,  lm3,  lm4,  lm5,  lm6,  lm7,  lm8,  lm9;
	double eqcd, ewjets, ezjets, etop, eother, emc, edata, elm1, elm2, elm3, elm4, elm5, elm6, elm7, elm8, elm9;//errors if one wants to use them
	if(susy==false){
		GetContent2DRange(hists[string("Events_") + string("QCD")],   qcd,   eqcd,   MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("WJets")], wjets, ewjets, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("ZJets")], zjets, ezjets, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("Top")],   top,   etop,   MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("Other")], other, eother, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("mc")],    mc,    emc,    MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("data")],  data,  edata,  MT2low, MT2up, HTlow, HTup);
	} else {
		GetContent2DRange(hists[string("Events_") + string("LM1")],     lm1,   elm1, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("LM2")],     lm2,   elm2, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("LM3")],     lm3,   elm3, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("LM4")],     lm4,   elm4, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("LM5")],     lm5,   elm5, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("LM6")],     lm6,   elm6, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("LM7")],     lm7,   elm7, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("LM8")],     lm8,   elm8, MT2low, MT2up, HTlow, HTup);
		GetContent2DRange(hists[string("Events_") + string("LM9")],     lm9,   elm9, MT2low, MT2up, HTlow, HTup);
	}

	if(MT2up==7000. && MT2low==0.) cout << " full selection &" ;
	else if(MT2up==7000.) cout  << " MT2 $\\ge$ " << int(MT2low) << " GeV & ";
	else            cout << " " << int(MT2low) << " GeV $\\le$ MT2 $<$ " << int(MT2up) << " GeV & ";
	if(susy==false){
		cout << fixed << setprecision(1)
		<< qcd << " & " << wjets << " & " << zjets << " & " << top << " & " << other << " & " << mc << " & "
		<< setprecision(0)
		<< data << "\\\\\n";
	}
	else {
		cout << fixed << setprecision(1)
		<< lm1 << " & " << lm2 << " & " << lm3 << " & " << lm4 << " & " << lm5 << " & " << lm6 << " & " << lm7 << " & " << lm8 << " & " << lm9
		<< "\\\\\n";
	}

}


//obtain the content sum of a range within a 2d histogram by defining the value range in x-y
void GetContent2DRange(TH2D* num, double &eff, double &err, double xlow, double xup, double ylow, double yup){
    int xbinlow, xbinup, ybinlow, ybinup;
    //change by 10e-6 to avoid border effects where one bin is included just as it is the lower border of a bin, so it in fact should not be included
    xbinlow = num->GetXaxis()->FindBin(xlow+10e-6);
    xbinup  = num->GetXaxis()->FindBin(xup -10e-6);
    ybinlow = num->GetYaxis()->FindBin(ylow+10e-6);
    ybinup  = num->GetYaxis()->FindBin(yup -10e-6);
    GetContent2DRangeBin(num, eff, err, xbinlow, xbinup, ybinlow, ybinup);
}
//obtain the content sum of a range within a 2d histogram by defining the bin range in x-y
void GetContent2DRangeBin(TH2D* num, double &eff, double &err, int xlow, int xup, int ylow, int yup){
    TH2D *num_c;
    num_c = (TH2D*) num->Clone();
    vector<double> xx, yy;
    for(int i = xlow; i<=xup; ++i){
        xx.push_back(num_c->GetXaxis()->GetBinLowEdge(i));}
    xx.push_back(num_c->GetXaxis()->GetBinLowEdge(xup)+num_c->GetXaxis()->GetBinWidth(xup));//upper border of last bin
    for(int i = ylow; i<=yup; ++i){
        yy.push_back(num_c->GetYaxis()->GetBinLowEdge(i));}
    yy.push_back(num_c->GetYaxis()->GetBinLowEdge(yup)+num_c->GetYaxis()->GetBinWidth(yup));//upper border of last bin
    int sizex = xx.size();
    int sizey = yy.size();
    double x[sizex], y[sizey];
    for(int i = 0; i<sizex; ++i){ x[i] = xx[i]; }
    for(int i = 0; i<sizey; ++i){ y[i] = yy[i]; }
    TH2D *num_sm = new TH2D("num_sm", "num_sm", sizex-1, x, sizey-1, y); num_sm->Sumw2();
    for(int nx=1; nx<=num_sm->GetNbinsX(); ++nx){
        for(int ny=1; ny<=num_sm->GetNbinsY(); ++ny){
            num_sm->SetBinContent(nx, ny, num_c->GetBinContent(num_c->FindBin(num_sm->GetXaxis()->GetBinCenter(nx), num_sm->GetYaxis()->GetBinCenter(ny) )));
            num_sm->SetBinError(nx, ny, num_c->GetBinError(num_c->FindBin(num_sm->GetXaxis()->GetBinCenter(nx), num_sm->GetYaxis()->GetBinCenter(ny) )));
        }
    }
    num_sm->RebinX( num_sm->GetNbinsX() );
    num_sm->RebinY( num_sm->GetNbinsY() );
    (eff) = (double)num_sm->GetBinContent(1,1);
    (err) = (double)num_sm->GetBinError(1,1);
    delete num_sm;
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

