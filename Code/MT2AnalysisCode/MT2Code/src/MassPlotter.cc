/*****************************************************************************
*   Small Class to make plots for MassAnalysis                               *
*****************************************************************************/

#include "MassPlotter.hh"

#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "helper/Utilities.hh"
#include "helper/Monitor.hh"

#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"
#include "TMath.h"

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TEventList.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TVirtualPad.h"
#include "TClonesArray.h"

#include "TGraph2D.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TVectorD.h"
#include "TPaveStats.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <time.h> // access to date/time


using namespace std;

//____________________________________________________________________________
MassPlotter::MassPlotter(){
	fSave=true;
	fisPhoton=false;
// Default constructor, no samples are set
}

//____________________________________________________________________________
MassPlotter::MassPlotter(TString outputdir){
// Explicit constructor with output directory
	setOutputDir(outputdir);
	fSave=true;
	fisPhoton=false;
}

//____________________________________________________________________________
MassPlotter::MassPlotter(TString outputdir, TString outputfile){
// Explicit constructor with output directory and output file
	setOutputDir(outputdir);
	setOutputFile(outputfile);
	fSave=true;
	fisPhoton=false;
}

//____________________________________________________________________________
MassPlotter::~MassPlotter(){
	fOutputFile->Close();
	delete fOutputFile;
}

//____________________________________________________________________________
void MassPlotter::init(TString filename){
	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Initializing MassPlotter ... " << endl;
	//Util::SetStyle();
	loadSamples(filename);

	// map for Z->nunu plots
	// to be used for Z->nunu           to be used for Z->ll
	RemoveLeptMap["misc.MET"]       = "Znunu.METplusLeptsPtReco";
}

//___________________________________________________________________________
void MassPlotter::makeSmallCopy(int nevents, int sample){
	if(fSamples.size()<sample) return;
	cout << "making small copy of tree " << fSamples[sample].name << endl;
	TFile *newfile = new TFile(fOutputDir+"/"+fSamples[sample].name+"_small.root","recreate");
	TTree *newtree = fSamples[sample].tree->CloneTree(0);
	Int_t nentries = (Int_t)fSamples[sample].tree->GetEntries();

	for (Int_t i=0;i<nentries; i++) {
		if(i > nevents) continue;
		fSamples[sample].tree->GetEntry(i);
		newtree->Fill();
	}


	newfile->Write();

	delete newfile;

}
//___________________________________________________________________________
void MassPlotter::makePlots(){
	double dPhisplit[4]={0., 1.0, 1.5, 3.142};
	MakeMT2PredictionAndPlots(false, dPhisplit, 2.5);  // not cleaned, fudgefactor 2
}


//__________________________________________________________________________
void MassPlotter::MakeMT2PredictionAndPlots(bool cleaned , double dPhisplit[], double fudgefactor){
	TString cleanflag;
	if(cleaned) cleanflag="cleaned";
	else        cleanflag="notcleaned";
	std::ostringstream o0, o1, o2, o3;
	o0 << dPhisplit[0];
	o1 << dPhisplit[1];
	o2 << dPhisplit[2];
	o3 << dPhisplit[3];
	TString dPhirange1 = (TString) o0.str() + "to" +(TString) o1.str();
	TString dPhirange2 = (TString) o2.str() + "to" +(TString) o3.str();
	TString dPhirange  = (TString) o0.str()+(TString) o1.str() + (TString) o2.str() +(TString) o3.str();
	
	std::vector<sample>  QCDSamples;
	std::vector<sample>  Samples_NOTData;
	std::vector<sample>  QCDandDataSamples;
	for(int i=0; i< fSamples.size(); ++i){
		if(fSamples[i].sname=="QCD"){
			QCDSamples.push_back(fSamples[i]);
			QCDandDataSamples.push_back(fSamples[i]);
		}
		if(fSamples[i].type!="data"){
			Samples_NOTData.push_back(fSamples[i]);
		}
		if(fSamples[i].type=="data"){
			QCDandDataSamples.push_back(fSamples[i]);
		}
	}

	double ZerotoPi[4]    ={0., 3.142, 3.142, 3.142};
	double controlsplit[4]={dPhisplit[2], dPhisplit[3], dPhisplit[3], dPhisplit[3]};

	std::ostringstream cutStream;
	cutStream  
		  << "misc.MET > 30"                                << "&&"
		  << "misc.HT > 300"                                << "&&"
		  << "misc.Jet0Pass == 1"                           << "&&"
		  << "misc.Jet1Pass == 1"                           << "&&"
		  << "misc.PassJetID == 1"                          << "&&"
		  << "misc.Vectorsumpt<70"                          << "&&"
		  << "misc.MinMetJetDPhi>0.3"                       << "&&" 
		  << "misc.EcalDeadCellBEFlag==1"                   << "&&"
		  << "misc.HBHENoiseFlag == 1"                   ;

	TString cuts = cutStream.str().c_str();
	
	//                 variable                            cuts njets  nlepts HLT title     bins               flip_order  log  composite    ratio  stacked overlay 
	MakePlot(fSamples,"misc.MT2" ,                         cuts, -3,   0,    "",  "MT2" , gNMT2bins, gMT2bins , false,  true ,  true,      true,  true,  false);
//	MakePlot(fSamples,"hemi[0].lv1.M()" ,                  cuts, -3,   0,    "",  "h1.M", 30,       0, 1000  , false,  true ,  true,      true,  true,  false);
//	MakePlot(fSamples,"hemi[0].lv2.M()" ,                  cuts, -3,   0,    "",  "h1.M", 30,       0, 1000  , false,  true ,  true,      true,  true,  false);

}

//________________________________________________________________________

void MassPlotter::makePlot(TString var, TString cuts, int njets, int nleps, TString HLT,  TString xtitle,
			   const int nbins, const double min, const double max,
			   bool flip_order, bool logflag, bool composited, bool ratio, bool stacked, 
			   bool overlaySUSY, float overlayScale, bool add_underflow){

  MakePlot(fSamples, var, cuts, njets, nleps, HLT, xtitle, nbins, min, max, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow);
  

}

//________________________________________________________________________

void MassPlotter::MakePlot(TString var, TString cuts, int njets, int nleps, TString HLT, TString xtitle, 
			   const int nbins,  const double *bins,
			   bool flip_order, bool logflag, bool composited, bool ratio, bool stacked, 
			   bool overlaySUSY, float overlayScale, bool add_underflow){

  MakePlot(fSamples, var, cuts, njets, nleps, HLT, xtitle, nbins, bins, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow);


}

// ________________________________________________________________________
void MassPlotter::PrintWEfficiency(int sample_index ,TString process,  std::string lept, Long64_t nevents, bool includeTaus){
	sample Sample = fSamples[sample_index];

      	std::cout << setfill('=') << std::setw(70) << "" << std::endl;
	cout << "printing W efficiency: \n"
	     << Sample.name << endl;	
      	std::cout << setfill('-') << std::setw(70) << "" << std::endl;

        enum counters_t { count_begin, all=count_begin, presel, NJetsIDLoose , PassJetID, MinMetJetDPhi, VectorSumPt, MT2,  count_end };
 	Monitor counters[count_end];
	TString lablx[count_end] = {"all events", "presel", "NJetsIDLoose" , "PassJetID", "MinMetJetDPhi", "VectorSumPt", "MT2"};

	string to_measure;
	if     (lept == "ele" && process=="W" )       {to_measure = "W->enu  recoed";} 
	else if(lept == "muo" && process=="W" )       {to_measure = "W->munu recoed";} 
	else if(lept == "ele" && process=="Top")      {to_measure = "Top W->enu  recoed";}
	else if(lept == "muo" && process=="Top")      {to_measure = "Top W->munu  recoed";}

    	fMT2tree = new MT2tree();
    	Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    	Long64_t nentries =  Sample.tree->GetEntries();
    	Long64_t nbytes = 0, nb = 0;
    	for (Long64_t jentry=0; jentry<min(nentries, nevents);jentry++) {
      		nb = Sample.tree->GetEntry(jentry);   nbytes += nb;
      		Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
      
      		if ( fVerbose>2 && jentry % 50000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;

		bool leptfound(false);
		bool eventgood(false);
		bool acceptance(false);
		if(process=="W"){
			if(lept=="ele" && fMT2tree->GenLeptFromW(11, 0 , 1000,includeTaus)==true)                                              eventgood =true;
			if(lept=="muo" && fMT2tree->GenLeptFromW(13, 0 , 1000,includeTaus)==true)                                              eventgood =true;
			if(lept=="ele" && fMT2tree->GenLeptFromW(11,10,  2.4 ,includeTaus)==true)                                              acceptance=true;
			if(lept=="muo" && fMT2tree->GenLeptFromW(13,10,  2.4 ,includeTaus)==true)                                              acceptance=true;
			if(lept=="ele" && fMT2tree->NEles==1 && fMT2tree->NMuons==0 && fMT2tree->GenLeptFromW(11, 0, 1000,includeTaus)==true && fMT2tree->ele[0].lv.Pt()>10)  leptfound =true;
			if(lept=="muo" && fMT2tree->NMuons==1&& fMT2tree->NEles== 0 && fMT2tree->GenLeptFromW(13, 0, 1000,includeTaus)==true && fMT2tree->muo[0].lv.Pt()>10)  leptfound =true;
		}else if(process=="Top"){ 
			if(lept=="ele" && fMT2tree->TopDecayModeResult(11)==true)                                                 eventgood =true;
			if(lept=="muo" && fMT2tree->TopDecayModeResult(13)==true)                                                 eventgood =true;
			if(lept=="ele" && fMT2tree->TopDecayModeResult(11)==true && fMT2tree->SLTopAccept(10, 2.4)==true )        acceptance =true;
			if(lept=="muo" && fMT2tree->TopDecayModeResult(13)==true && fMT2tree->SLTopAccept(10, 2.4)==true )        acceptance =true;
			if(lept=="ele" && fMT2tree->NEles==1 && fMT2tree->NMuons==0 && eventgood &&fMT2tree->ele[0].lv.Pt()>10)             leptfound =true;
			if(lept=="muo" && fMT2tree->NMuons==1&& fMT2tree->NEles== 0 && eventgood &&fMT2tree->muo[0].lv.Pt()>10)             leptfound =true;
		}

		Double_t weight = fMT2tree->pileUp.Weight;

		// all events
		if(eventgood)   counters[all].fill("all events", weight);
		if(acceptance)  counters[all].fill("acceptance", weight);
		if(leptfound)   counters[all].fill(to_measure, weight);
		
		// presel
		if(fMT2tree->misc.MET                     < 30    )  continue;
		if(fMT2tree->misc.HT                      < 300   )  continue;
		//if(fMT2tree->misc.caloHT50_ID             < 600   )  continue;
		if(fMT2tree->misc.Jet0Pass                ==0     )  continue;
		if(fMT2tree->misc.Jet1Pass                ==0     )  continue;
//		if(fMT2tree->misc.LeadingJPt              < 150   )  continue;
		if(fMT2tree->misc.SecondJPt               < 100   )  continue;
//		if(fMT2tree->NBJets                       < 1     )  continue;
		if(fMT2tree->misc.HBHENoiseFlag           ==0     )  continue;
		if(fMT2tree->misc.CrazyHCAL               ==1     )  continue;
		if(eventgood)   counters[presel].fill("presel", weight);
		if(acceptance)  counters[presel].fill("acceptance", weight);
		if(leptfound)   counters[presel].fill(to_measure, weight);
		
		
		// NJetsIDLoose
		if(fMT2tree->NJetsIDLoose                 <3     ) continue;  
		if(eventgood)   counters[NJetsIDLoose].fill("NJetsIDLoose", weight);
		if(acceptance)  counters[NJetsIDLoose].fill("acceptance", weight);
		if(leptfound)   counters[NJetsIDLoose].fill(to_measure, weight);

		// PassJetID
		if(fMT2tree->misc.PassJetID               ==0     ) continue;
		if(eventgood)   counters[PassJetID].fill("PassJetID", weight);
		if(acceptance)  counters[PassJetID].fill("acceptance", weight);
		if(leptfound)   counters[PassJetID].fill(to_measure, weight);
		
		// MinMetJetDPhi
		if(fMT2tree->misc.MinMetJetDPhi <0.3 )       continue;
		if(eventgood)   counters[MinMetJetDPhi].fill("MinMetJetDPhi", weight);
		if(acceptance)  counters[MinMetJetDPhi].fill("acceptance", weight);
		if(leptfound)   counters[MinMetJetDPhi].fill(to_measure, weight);
		
		// VectorSumPt
		if(fMT2tree->misc.Vectorsumpt             >70     )  continue;
		if(eventgood)   counters[VectorSumPt].fill("VectorSumPt", weight);
		if(acceptance)  counters[VectorSumPt].fill("acceptance", weight);
		if(leptfound)   counters[VectorSumPt].fill(to_measure, weight);
		
		// MT2
//		if(fMT2tree->misc.MT2             <400     )  continue;
//		if(fMT2tree->misc.MT2             <150     )  continue;
		if(fMT2tree->misc.MT2 <200 || fMT2tree->misc.MT2>400 )  continue;
//		if(fMT2tree->misc.MT2 <150 || fMT2tree->misc.MT2>300 )  continue;
//		if(fMT2tree->misc.MT2 <100 || fMT2tree->misc.MT2>150  )  continue;
		if(eventgood)   counters[MT2].fill("MT2", weight);
		if(acceptance)  counters[MT2].fill("acceptance", weight);
		if(leptfound)   counters[MT2].fill(to_measure, weight);
		
      	}

	// print stats     
      	std::cout << setfill('=') << std::setw(70) << "" << std::endl;
       	std::cout << "Statistics" << std::endl;
       	std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;
       	for ( counters_t iCount=count_begin; iCount<count_end; iCount = counters_t(iCount+1) ) {
         	counters[iCount].print();
         	std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;  
       	}
  	std::cout << setfill('=') << std::setw(70) << "" << std::endl;

	// fill histo
	TH1D* h_all = new TH1D("h_all", "", count_end, 0., (double) count_end );
	TH1D* h_acc = new TH1D("h_acc", "", count_end, 0., (double) count_end );
	TH1D* h_rec = new TH1D("h_rec", "", count_end, 0., (double) count_end );
	TH1D* h_1   = new TH1D("h_1"  , "", count_end, 0., (double) count_end );
	TH1D* h_2   = new TH1D("h_2"  , "", count_end, 0., (double) count_end );
	TH1D* h_3   = new TH1D("h_3"  , "", count_end, 0., (double) count_end );
	for(int i=0; i<count_end; ++i){
		h_1    ->GetXaxis()->SetBinLabel(i+1, lablx[i]);
		h_2    ->GetXaxis()->SetBinLabel(i+1, lablx[i]);
		h_3    ->GetXaxis()->SetBinLabel(i+1, lablx[i]);
		h_all->SetBinContent(i+1,counters[i].counts((string) lablx[i]));
		h_acc->SetBinContent(i+1,counters[i].counts("acceptance"));
		h_rec->SetBinContent(i+1,counters[i].counts(to_measure));
	}
	h_1    ->Sumw2();
	h_2    ->Sumw2();
	h_3    ->Sumw2();
	h_1    ->Divide(h_acc,h_all);	
	h_2    ->Divide(h_rec,h_acc);	
	h_3    ->Divide(h_rec,h_all);	
	TCanvas *col  = new TCanvas("Wevents", "", 0, 0, 900, 700);
	col -> cd();
	h_1    ->SetLineColor(kBlue);
	h_1    ->SetMarkerStyle(20);
	h_1    ->SetMinimum(0);
	h_1    ->SetMaximum(1);
	h_1    ->SetDrawOption("E");
	h_1    ->Draw();
	string name  ="Wevents_"+lept+"_acceptance_"+(string) Sample.name; 
	h_2    ->SetLineColor(kBlue);
	h_2    ->SetMinimum(0);
	h_2    ->SetMaximum(1);
	h_2    ->SetMarkerStyle(20);
	h_2    ->SetDrawOption("E");
	h_2    ->Draw();
	name  ="Wevents_"+lept+"_reco_"+(string) Sample.name; 
	h_3    ->SetLineColor(kBlue);
	h_3    ->SetMinimum(0);
	h_3    ->SetMaximum(1);
	h_3    ->SetMarkerStyle(20);
	h_3    ->SetDrawOption("E");
	h_3    ->Draw();
	name  ="Wevents_"+lept+"_prob_"+(string) Sample.name; 

	cout << "lept " << lept << " process " << process << endl;	
	//__________________________________
	if        (lept=="ele" && process=="W") {
		fWpred.Wenu_acc        =h_1->GetBinContent(count_end);
		fWpred.Wenu_acc_err    =h_1->GetBinError(count_end); 
		fWpred.Wenu_rec        =h_2->GetBinContent(count_end);
		fWpred.Wenu_rec_err    =h_2->GetBinError(count_end);  
		fWpred.Wenu_prob       =h_3->GetBinContent(count_end);
		fWpred.Wenu_prob_err   =h_3->GetBinError(count_end);  
	}else if  (lept=="ele" && process=="Top"){
		fWpred.TopWenu_acc     =h_1->GetBinContent(count_end);
		fWpred.TopWenu_acc_err =h_1->GetBinError(count_end); 
		fWpred.TopWenu_rec     =h_2->GetBinContent(count_end);
		fWpred.TopWenu_rec_err =h_2->GetBinError(count_end);  
		fWpred.TopWenu_prob    =h_3->GetBinContent(count_end);
		fWpred.TopWenu_prob_err=h_3->GetBinError(count_end);  
	}else if  (lept=="muo" && process=="W") {
		fWpred.Wmunu_acc       =h_1->GetBinContent(count_end);
		fWpred.Wmunu_acc_err   =h_1->GetBinError(count_end); 
		fWpred.Wmunu_rec       =h_2->GetBinContent(count_end);
		fWpred.Wmunu_rec_err   =h_2->GetBinError(count_end); 
		fWpred.Wmunu_prob      =h_3->GetBinContent(count_end);
		fWpred.Wmunu_prob_err  =h_3->GetBinError(count_end);  
	}else if  (lept=="muo" && process=="Top") {
		fWpred.TopWmunu_acc    =h_1->GetBinContent(count_end);
		fWpred.TopWmunu_acc_err=h_1->GetBinError(count_end); 
		fWpred.TopWmunu_rec    =h_2->GetBinContent(count_end);
		fWpred.TopWmunu_rec_err=h_2->GetBinError(count_end); 
		fWpred.TopWmunu_prob    =h_3->GetBinContent(count_end);
		fWpred.TopWmunu_prob_err=h_3->GetBinError(count_end);  
	}

	//_________________________________
	delete h_1;
	delete h_2;
	delete h_3;
	delete h_all;
	delete h_acc;
	delete h_rec;
	delete col;
}

// ________________________________________________________________________
void MassPlotter::PrintZllEfficiency(int sample_index , bool data, std::string lept, Long64_t nevents, double lower_mass, double upper_mass, bool pileUp_reweight){
//void MassPlotter::PrintZllEfficiency(int sample_index , bool data, std::string lept, Long64_t nevents, double lower_mass, double upper_mass){
	sample Sample = fSamples[sample_index];

      	std::cout << setfill('=') << std::setw(70) << "" << std::endl;
	cout << "printing kinetic & geom. acceptance for Z->ll for sample: \n"
	     << Sample.name << endl;	
      	std::cout << setfill('-') << std::setw(70) << "" << std::endl;

        enum counters_t { count_begin, all=count_begin, presel,  HCAL_ECAL_noise, VectorSumPt ,PassJetID, MinMetJetDPhi, MET, count_end };
 	Monitor counters[count_end];
	TString lablx[count_end] = {"all events", "presel", "Cal_Noise", "VectorSumPt", "PassJetID", "MinMetJetDPhi", "MET"};

	int pid=-1, flavour=-1; 
	string to_measure;
	if(lept == "ele" )            {pid = 11; flavour = 1; to_measure = "Zee recoed";} 
	else if(lept == "muo" )       {pid = 13; flavour = 2; to_measure = "Zmumu recoed";} 
	else if(lept == "neutrinos")  {pid = 12; flavour = 0; to_measure = "Znunu within acceptance";}

    	fMT2tree = new MT2tree();
    	Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    	Long64_t nentries =  Sample.tree->GetEntries();
    	Long64_t nbytes = 0, nb = 0;
    	for (Long64_t jentry=0; jentry<min(nentries, nevents);jentry++) {
      		nb = Sample.tree->GetEntry(jentry);   nbytes += nb;
      		Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
      
      		if ( fVerbose>2 && jentry % 50000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;

		// check if event has Z->e+e- within acceptance	
		
		bool Zll(false);         
		bool ZllSelect(false); 
		bool ZllGood(false); 
		if(  data) Zll              = (fMT2tree->GetDiLeptonInvMass(0,1,flavour,5,1) > 71 && fMT2tree->GetDiLeptonInvMass(0,1,flavour,5,1) < 111);
		if(! data){ 
			if     (pid ==11) {
				ZllSelect =   (fMT2tree->Znunu.GenZee_mll_acc   > lower_mass && fMT2tree->Znunu.GenZee_mll_acc   < upper_mass  );
				ZllGood     = (fMT2tree->Znunu.RecoOSee_mll     > lower_mass && fMT2tree->Znunu.RecoOSee_mll     < upper_mass  && ZllSelect );
			}else if(pid ==13) {
				ZllSelect =   (fMT2tree->Znunu.GenZmumu_mll_acc > lower_mass && fMT2tree->Znunu.GenZmumu_mll_acc < upper_mass  );
				ZllGood     = (fMT2tree->Znunu.RecoOSmumu_mll   > lower_mass && fMT2tree->Znunu.RecoOSmumu_mll   < upper_mass  && ZllSelect  );
			}else if(lept=="neutrinos"){
				ZllSelect   =  ((fMT2tree->Znunu.GenZnunu_e_mll        > lower_mass && fMT2tree->Znunu.GenZnunu_e_mll    < upper_mass)  || 
					 	(fMT2tree->Znunu.GenZnunu_mu_mll       > lower_mass && fMT2tree->Znunu.GenZnunu_mu_mll   < upper_mass)  ||  
					        (fMT2tree->Znunu.GenZnunu_tau_mll      > lower_mass && fMT2tree->Znunu.GenZnunu_tau_mll  < upper_mass));	
				ZllGood     =  ((fMT2tree->Znunu.GenZnunu_e_mll_acc    > lower_mass && fMT2tree->Znunu.GenZnunu_e_mll_acc    < upper_mass)  || 
					 	(fMT2tree->Znunu.GenZnunu_mu_mll_acc   > lower_mass && fMT2tree->Znunu.GenZnunu_mu_mll_acc   < upper_mass)  ||  
					        (fMT2tree->Znunu.GenZnunu_tau_mll_acc  > lower_mass && fMT2tree->Znunu.GenZnunu_tau_mll_acc  < upper_mass));	
			}else {cout << "choose appropriate lepton" << endl; return;}

		}

		Double_t weight =(pileUp_reweight)? fMT2tree->pileUp.Weight:1;
			
		// all events
		if(ZllSelect)   counters[all].fill("all events", weight);
		if(ZllGood    ) counters[all].fill(to_measure  , weight);
		
		// presel
		if(pid == 11 || pid == 13){
			if(fMT2tree->Znunu.HTmatched              <300    )  continue;
			if(fMT2tree->Znunu.caloHT50ID_matched     <600    )  continue;
			if(fMT2tree->Znunu.Jet0Pass_matched       ==0     )  continue;
			if(fMT2tree->Znunu.Jet1Pass_matched       ==0     )  continue;
			if(fMT2tree->Znunu.SecondJPt_matched      <100     )  continue;
			if(fMT2tree->Znunu.NJetsIDLoose_matched   <3      )  continue;
		} else if(lept == "neutrinos"){
			if(fMT2tree->misc.HT                      <300    )  continue;
			//if(fMT2tree->misc.caloHT50_ID             <600    )  continue;
			if(fMT2tree->misc.Jet0Pass                ==0     )  continue;
			if(fMT2tree->misc.Jet1Pass                ==0     )  continue;
			if(fMT2tree->misc.SecondJPt               <100     )  continue;
			if(fMT2tree->NJetsIDLoose                  <3     )  continue;
		}
		if(ZllSelect)   counters[presel].fill("presel"  , weight);
		if(ZllGood    ) counters[presel].fill(to_measure, weight);
		
		// ECAL HCAL Noise
		if(fMT2tree->misc.HBHENoiseFlag                   ==0     )  continue;
		if(fMT2tree->misc.CrazyHCAL                       ==1     )  continue;
		if(ZllSelect  ) counters[HCAL_ECAL_noise].fill("Cal_Noise",weight);
		if(ZllGood    ) counters[HCAL_ECAL_noise].fill(to_measure ,weight);
		
		// VectorSumPt
		if(pid == 11 || pid == 13){
			if(fMT2tree->Znunu.Vectorsumpt_matched    >70     )  continue;
		} else if(lept == "neutrinos"){
			if(fMT2tree->misc.Vectorsumpt             >70     )  continue;
		}
		if(ZllSelect)   counters[VectorSumPt].fill("VectorSumPt", weight);
		if(ZllGood    ) counters[VectorSumPt].fill(to_measure,    weight);
		
		// PassJetID
		if(pid ==11 || pid == 13){
			if(fMT2tree->Znunu.PassJetID_matched      ==0     ) continue;
		}else if (lept == "neutrinos"){
			if(fMT2tree->misc.PassJetID               ==0     ) continue;
		}
		if(ZllSelect)   counters[PassJetID].fill("PassJetID"   ,weight);
		if(ZllGood    ) counters[PassJetID].fill(to_measure    ,weight);
		
		// MinMetJetDPhi
		if(pid ==11 || pid == 13){
			if(fMT2tree->Znunu.MinMetplusLeptJetDPhi  <0.3    )  continue;
		}else if (lept == "neutrinos"){
			if(fMT2tree->misc.MinMetJetDPhi           <0.3    )  continue;
		}
		if(ZllSelect)   counters[MinMetJetDPhi].fill("MinMetJetDPhi",weight);
		if(ZllGood    ) counters[MinMetJetDPhi].fill(to_measure     ,weight);
		
		// MET
		if(pid ==11 || pid == 13){
			if(fMT2tree->Znunu.METplusLeptsPt         < 30    ) continue;
		}else if (lept == "neutrinos"){
			if(fMT2tree->misc.MET                     < 30    )  continue;
		}
		if(ZllSelect)   counters[MET].fill("MET"       ,weight);
		if(ZllGood    ) counters[MET].fill(to_measure  ,weight);
      	}

	// print stats     
      	std::cout << setfill('=') << std::setw(70) << "" << std::endl;
       	std::cout << "Statistics" << std::endl;
       	std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;
       	for ( counters_t iCount=count_begin; iCount<count_end; iCount = counters_t(iCount+1) ) {
         	counters[iCount].print();
         	std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;  
       	}
  	std::cout << setfill('=') << std::setw(70) << "" << std::endl;

	if(data) return;

	// fill histo
	TH1D* h_Z_num = new TH1D("h_Z_num", "", count_end, 0., (double) count_end );
	TH1D* h_Z_den = new TH1D("h_Z_den", "", count_end, 0., (double) count_end );
	TH1D* h_Z     = new TH1D("h_Z"    , "", count_end, 0., (double) count_end );
	for(int i=0; i<count_end; ++i){
		h_Z    ->GetXaxis()->SetBinLabel(i+1, lablx[i]);
		h_Z_num->SetBinContent(i+1,counters[i].counts(to_measure));
		h_Z_den->SetBinContent(i+1,counters[i].counts((string) lablx[i]));
	}
	h_Z    ->Sumw2();
	h_Z    ->Divide(h_Z_num,h_Z_den);	
	TCanvas *col  = new TCanvas("GenDYToLL", "", 0, 0, 900, 700);
	col -> cd();
	h_Z    ->SetLineColor(kBlue);
	h_Z    ->SetMarkerStyle(20);
	h_Z    ->SetDrawOption("E");
	h_Z    ->Draw();
	string name  ="GenDYToLL_"+lept+"_"+(string) Sample.name; 
	if(fSave) Util::PrintEPS(col, name, fOutputDir);
	

	//__________________________________
	if(lept=="ele") {
		fZpred.ele_reco     =h_Z  ->GetBinContent(count_end); 
		fZpred.ele_reco_err =h_Z  ->GetBinError  (count_end);
	}
	if(lept=="muo") {
		fZpred.muo_reco     =h_Z  ->GetBinContent(count_end); 
		fZpred.muo_reco_err =h_Z  ->GetBinError  (count_end); 
	}
	if(lept=="neutrinos") {
		fZpred.nu_acc       =h_Z  ->GetBinContent(count_end);
		fZpred.nu_acc_err   =h_Z  ->GetBinError  (count_end);
	}
	//_________________________________
	delete h_Z;
	delete h_Z_den;
	delete h_Z_num;
	delete col;
}

//________________________________________________________________________

void MassPlotter::PrintCutFlow(int njets, int nleps, TString trigger, TString cuts){
  
  Monitor counters[fSamples.size()];

  const int   nProc      = 17;  
  TString  cnames[nProc] = {"QCD", "W+jets", "Z+jets", "Top","Other", "Total Bkg.", "data", "LM1", "LM2", "LM3", "LM4", "LM5", "LM8", "LM9", "LM11", "LM12", "LM13"};
  Monitor  ccount[nProc], ccount_100[nProc];

  for(size_t i = 0; i < fSamples.size(); ++i){

    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents * fSamples[i].PU_avg_weight );
    if(fVerbose>2) cout << "PrintCutFlow: looping over " << fSamples[i].name << endl;
    if(fVerbose>2) cout << "              sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
    if(fVerbose>2) cout << "              Average PU weight: " << fSamples[i].PU_avg_weight << endl;
    if(fVerbose>2) cout << "              Original Entries: " << fSamples[i].nevents << endl;
    
    fMT2tree = new MT2tree();
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
    
    while(myEvtList->GetEntry(counter++) !=-1){
      
      int jentry = myEvtList->GetEntry(counter-1);
      
      //for (Long64_t jentry=0; jentry<nentries;jentry++) {
      nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
      fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
      
      if ( fVerbose>2 && jentry % 50000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;
      
      Double_t weight = sample_weight;
      if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 

       bool isMT2gt100 = fMT2tree->misc.MT2 > 100.;

      if( fMT2tree->NJetsIDLoose < 1 || !fMT2tree->misc.Jet0Pass                                      )  continue;
      if( fMT2tree->NJetsIDLoose < 2 || !fMT2tree->misc.Jet1Pass    || fMT2tree->misc.SecondJPt < 100 )  continue;
      if( !(!fMT2tree->misc.isData   || (fMT2tree->misc.Run <162803 || fMT2tree->misc.Run >162909  )) )  continue;
      if(trigger=="HT"){
	      if( fMT2tree->misc.isData ==1 
		  && fMT2tree->trigger.HLT_PFHT650_v5 ==0 
		  && fMT2tree->trigger.HLT_PFHT650_v6 ==0 
		  && fMT2tree->trigger.HLT_PFHT650_v7 ==0 ) continue;
              //if( fMT2tree->misc.caloHT50_ID< 550)        continue;
              if( fMT2tree->misc.MET     < 30)  continue;
      }else if(trigger=="MHT_HT"){
      	      if( fMT2tree->misc.isData ==1 
		  && fMT2tree->trigger.HLT_PFHT350_PFMET100_v3 ==0 
		  && fMT2tree->trigger.HLT_PFHT350_PFMET100_v4 ==0 
		  && fMT2tree->trigger.HLT_PFHT350_PFMET100_v5 ==0 ) continue;
              if( fMT2tree->misc.caloHT50  <= 320) continue;
              if( fMT2tree->misc.caloMHT30 <  130) continue;
              if( fMT2tree->misc.MET       <   30) continue;
	      if( fMT2tree->misc.HT        <  450) continue;
      }

      if     (njets>1) {
	if(fMT2tree->NJetsIDLoose != njets) continue;
	TString text = TString::Format("All events (jets == %d)",njets);
	counters[i].fill(text.Data(),weight);
	FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, text, weight);
	if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, text, weight);
      }
      else if(njets <-1) {
	if(fMT2tree->NJetsIDLoose < abs(njets)) continue;
	TString text = TString::Format("All events(jets >= %d)",abs(njets));
	counters[i].fill(text.Data(),weight);
	FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, text, weight);
	if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, text, weight);
      }

      if(fMT2tree->PassMinMetJetDPhi03() ==0)      continue;      
      counters[i].fill("Minimum DPhi(MET,jet) > 0.3",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "Minimum DPhi(MET,jet) > 0.3", weight);
      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "Minimum DPhi(MET,jet) > 0.3", weight);

      if( fMT2tree->misc.HBHENoiseFlag != 1 )  continue;
      if( fMT2tree->misc.CrazyHCAL     != 0 )  continue;
      counters[i].fill("HBHE noise veto",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "HBHE noise veto", weight);
      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "HBHE noise veto", weight);

//      if( fMT2tree->misc.EcalDeadCellBEFlag != 1 )  continue;
//      counters[i].fill("Boundary energy veto (dead ecal)",weight);
//      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "Boundary energy veto (dead ecal)", weight);
//      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "Boundary energy veto (dead ecal)", weight);

      if( fMT2tree->misc.Vectorsumpt > 70. )  continue;
      counters[i].fill("VectorSumPt < 70",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "VectorSumPt < 70", weight);
      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "VectorSumPt < 70", weight);

      if( !fMT2tree->misc.PassJetID )  continue;
      counters[i].fill("jets > 50GeV failing PFID event veto",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "jets > 50GeV failing PFID event veto", weight);
      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "jets > 50GeV failing PFID event veto", weight);

      // Only considering (so far): no lepton requirement (<0);  1 lepton (==1), lepton veto (otherwise)
      string nLeps = (std::string) TString::Format("%d",abs(nleps));
      if(nleps !=-10){ 
	      if(nleps>0){ // exactly nleps lepton
		string cut_name = "NLeps == "+nLeps;
		if( (fMT2tree->NEles + fMT2tree->NMuons) != abs(nleps))  continue;
		counters[i].fill(cut_name,weight);
		FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, cut_name, weight);
		if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, cut_name, weight);
	      }
	      else if( nleps==0 ) { // lepton veto
		string cut_name = "Lepton veto";
		if( (fMT2tree->NEles + fMT2tree->NMuons) != 0 )  continue;
		counters[i].fill(cut_name,weight);
		FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, cut_name, weight);
		if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, cut_name, weight);
	      }
	      else if(nleps<0){ // at least nleps lepton
		string cut_name = "NLeps >= "+nLeps;
		if( (fMT2tree->NEles + fMT2tree->NMuons) < abs(nleps))  continue;
		counters[i].fill(cut_name ,weight);
		FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, cut_name, weight);
		if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, cut_name, weight);
	      }
      }
      
//      if( !(fMT2tree->GetNBtags(2,1.74) >= 1) ) continue;
//      TString text = "b-tags >=1";
//      counters[i].fill(text.Data(),weight);
//      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, text, weight);
//      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, text, weight);
      
      if( fMT2tree->misc.MT2 < 80. )  continue;
      counters[i].fill("MT2 > 80 GeV" ,weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 80 GeV", weight);
      if( fMT2tree->misc.MT2 < 100. )  continue;
      counters[i].fill("MT2 > 100 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 100 GeV", weight);
      if( fMT2tree->misc.MT2 < 120. )  continue;
      counters[i].fill("MT2 > 120 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 120 GeV", weight);
      if( fMT2tree->misc.MT2 < 135. )  continue;
      counters[i].fill("MT2 > 135 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 135 GeV", weight);
      if( fMT2tree->misc.MT2 < 150. )  continue;
      counters[i].fill("MT2 > 150 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 150 GeV", weight);
      if( fMT2tree->misc.MT2 < 165. )  continue;
      counters[i].fill("MT2 > 165 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 165 GeV", weight);
      if( fMT2tree->misc.MT2 < 180. )  continue;
      counters[i].fill("MT2 > 180 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 180 GeV", weight);
      if( fMT2tree->misc.MT2 < 200. )  continue;
      counters[i].fill("MT2 > 200 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 200 GeV", weight);
      if( fMT2tree->misc.MT2 < 225. )  continue;
      counters[i].fill("MT2 > 225 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 225 GeV", weight);
      if( fMT2tree->misc.MT2 < 250. )  continue;
      counters[i].fill("MT2 > 250 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 250 GeV", weight);
      if( fMT2tree->misc.MT2 < 275. )  continue;
      counters[i].fill("MT2 > 275 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 275 GeV", weight);
      if( fMT2tree->misc.MT2 < 300. )  continue;
      counters[i].fill("MT2 > 300 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 300 GeV", weight);
      if( fMT2tree->misc.MT2 < 325. )  continue;
      counters[i].fill("MT2 > 325 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 325 GeV", weight);
      if( fMT2tree->misc.MT2 < 350. )  continue;
      counters[i].fill("MT2 > 350 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 350 GeV", weight);
      if( fMT2tree->misc.MT2 < 375. )  continue;
      counters[i].fill("MT2 > 375 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 375 GeV", weight);
      if( fMT2tree->misc.MT2 < 400. )  continue;
      counters[i].fill("MT2 > 400 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 400 GeV", weight);
      if( fMT2tree->misc.MT2 < 425. )  continue;
      counters[i].fill("MT2 > 425 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 425 GeV", weight);
      if( fMT2tree->misc.MT2 < 450. )  continue;
      counters[i].fill("MT2 > 450 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 450 GeV", weight);
      if( fMT2tree->misc.MT2 < 475. )  continue;
      counters[i].fill("MT2 > 475 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 475 GeV", weight);
      if( fMT2tree->misc.MT2 < 500. )  continue;
      counters[i].fill("MT2 > 500 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 500 GeV", weight);
      
    }
    delete fMT2tree;
  }

  std::cout << setfill('=') << std::setw(70) << "" << std::endl;
  std::cout << "Statistics - by sample" << std::endl;
  std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
  for ( size_t i = 0; i < fSamples.size(); ++i){
    std::cout << "++++  " << fSamples[i].name << std::endl;  
    std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
    counters[i].print();
    std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;    
  }  

  std::cout << setfill('=') << std::setw(70) << "" << std::endl;
  std::cout << "Statistics - by process" << std::endl;
  std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
  for ( size_t i = 0; i < nProc; ++i){
    std::cout << "++++  " << cnames[i] << std::endl;  
    std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
    ccount[i].print();
    std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;    
  }  

  std::cout << setfill('=') << std::setw(70) << "" << std::endl;
  std::cout << "Statistics - by process (MT2 > 100GeV)" << std::endl;
  std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
  for ( size_t i = 0; i < nProc; ++i){
    std::cout << "++++  " << cnames[i] << std::endl;  
    std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
    ccount_100[i].print();
    std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;    
  }  

}

// MT2 vs HT cutflot  --------------------------------------------------------------------
void MassPlotter::PrintCutFlowMT2vsHT(TString trigger, TString cuts){
  
  // define here the HT & MT2 sampling to be used
  int   numHTbins   =  15;
  float minHT       = 400;
  float HTgridsize  =  50;
  float minMT2cut   = 100;
  float maxMT2cut   = 700;
  float MT2gridsize =  25;
  
  const int   nProc      = 17;  
  TString  cnames[nProc] = {"QCD", "W+jets", "Z+jets", "Top", "Other", "Total Bkg.", "data", "LM1", "LM2", "LM3", "LM4", "LM5", "LM8", "LM9", "LM11", "LM12", "LM13"};
  Monitor Monitors[numHTbins][nProc];
  Monitor counters[fSamples.size()];

  // Loop over samples
  for(size_t i = 0; i < fSamples.size(); ++i){

    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents * fSamples[i].PU_avg_weight );
    if(fVerbose>2) cout << "PrintCutFlow: looping over " << fSamples[i].name << endl;
    if(fVerbose>2) cout << "              sample has weight " << sample_weight << " and input tree has " << fSamples[i].tree->GetEntries() << " entries" << endl; 
    if(fVerbose>2) cout << "              Average PU weight: " << fSamples[i].PU_avg_weight << endl;
    if(fVerbose>2) cout << "              Original Entries: " << fSamples[i].nevents << endl;

    
    fMT2tree = new MT2tree();
    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    int nev =0;

    // don't do anything for empty trees: (maybe tree already skimmed)
    if (fSamples[i].tree->GetEntries()==0) continue;

    //leo tweak - filtering out the TTree
    TString myCuts = cuts;

    if( fSamples[i].type=="data" && trigger!="") myCuts += " && " + trigger; //cuts to be aplied only on data
    cout << "              Cuts for Flow: " << myCuts << endl;
    fSamples[i].tree->Draw(">>selList", myCuts);


    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
    fSamples[i].tree->SetEventList(myEvtList);
    int counter=0;
    cout << "              Filtering done, size=" <<myEvtList->GetN()  << endl;
    
    if(myEvtList->GetSize()==0) continue;
    
    // loop over entried in MT2trees
    while(myEvtList->GetEntry(counter++) !=-1){
      
      int jentry = myEvtList->GetEntry(counter-1);
      
      nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
      fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
      
      if ( fVerbose>2 && jentry % 50000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;
      
      Double_t weight = sample_weight;
      if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 

	counters[i].fill("All events",weight); // fill all events that pass selection, without MT2 cut
	float MT2cut =minMT2cut; // initial cut on MT2
	while(MT2cut <=maxMT2cut ){
		if(fMT2tree->misc.MT2 >MT2cut){
			TString cutname = TString::Format("MT2 > %.0f GeV",MT2cut);
			counters[i].fill(cutname.Data(),weight);
			for(int HTbin=0; HTbin< numHTbins; ++HTbin){
				float HTcut = minHT + HTbin*HTgridsize;
				if( fMT2tree->misc.HT > HTcut)   FillMonitor(Monitors[HTbin], fSamples[i].sname, fSamples[i].type,cutname.Data() , weight);
			
			}
		}
		MT2cut+=MT2gridsize;
	}
    }
    delete fMT2tree;
  }

  std::ostringstream logStream;
  // by sample, no HT cuts
  logStream << setfill('=') << std::setw(70) << "" << std::endl;
  logStream << "Statistics - by sample, no HT cut" << std::endl;
  logStream << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
  for ( size_t i = 0; i < fSamples.size(); ++i){
    logStream << "++++  " << fSamples[i].name << std::endl;  
    logStream << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
    logStream << counters[i];
    logStream << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;    
  }  
  ofstream f_log ("Cutflow_bySample_defaultcuts.log", ios::app);
  f_log << logStream.str();
  logStream.str(""); // clear the ostringstream

  for ( int HTbin =0; HTbin < numHTbins; ++HTbin){
	  logStream << setfill('=') << std::setw(70) << "" << std::endl;
	  TString HTcut= TString::Format("%.0f", minHT+HTbin*HTgridsize);
	  logStream << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
	  logStream << "the following cuts were applied: " << std::endl;
	  logStream << cuts                                << std::endl;
 	  logStream << "with pfHT > " << HTcut	<< " GeV " << std::endl;  
	  logStream << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
	  logStream << "Statistics - by process (HT > " << HTcut << " GeV)" << std::endl;
	  logStream << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
	  for ( size_t i = 0; i < nProc; ++i){
	    logStream << "++++  " << cnames[i] << std::endl;  
	    logStream << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
	    logStream << Monitors[HTbin][i];
	    logStream << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;    
	  }  

	TString cutflowlog = "Cutflow_HT" +HTcut + "_byProcess.log";
	ofstream f_log (cutflowlog.Data(), ios::app);
	f_log << logStream.str();
	logStream.str(""); // clear the ostringstream
  }

}

//________________________________________________________________________

void MassPlotter::FillMonitor(Monitor *count, TString sname, TString type, TString cut, double weight){
  if     (sname == "QCD"    ) 	count[ 0].fill(cut.Data(),weight);
  else if(sname == "Wtolnu" ) 	count[ 1].fill(cut.Data(),weight);
  else if(sname == "DY"     ) 	count[ 2].fill(cut.Data(),weight);
  else if(sname == "Top"    )	count[ 3].fill(cut.Data(),weight);
  else if(sname == "VV" || sname == "VGamma" || sname == "Photons") 	
                                count[ 4].fill(cut.Data(),weight);
  if     (type  == "mc"     )   count[ 5].fill(cut.Data(),weight);
  else if(type  == "data"   )	count[ 6].fill(cut.Data(),weight);
  else if(sname == "LM1"    )	count[ 7].fill(cut.Data(),weight);
  else if(sname == "LM2"    )	count[ 8].fill(cut.Data(),weight);
  else if(sname == "LM3"    )	count[ 9].fill(cut.Data(),weight);
  else if(sname == "LM4"    )	count[10].fill(cut.Data(),weight);
  else if(sname == "LM5"    )	count[11].fill(cut.Data(),weight);
  else if(sname == "LM8"    )	count[12].fill(cut.Data(),weight);
  else if(sname == "LM9"    )	count[13].fill(cut.Data(),weight);
  else if(sname == "LM11"    )	count[14].fill(cut.Data(),weight);
  else if(sname == "LM12"    )	count[15].fill(cut.Data(),weight);
  else if(sname == "LM13"    )	count[16].fill(cut.Data(),weight);
}

//________________________________________________________________________

void MassPlotter::plotSig(TString var, TString cuts, TString xtitle, int nbins, double min, double max, bool flip_order, int type ){

        TString varname = Util::removeFunnyChar(var.Data());

	TH1D*    h_mc_sum    = new TH1D   (varname+"mc_sum", "", nbins, min, max );
	TH1D*    h_susy      = new TH1D   (varname+"susy"  , "", nbins, min, max );	
	h_mc_sum    ->Sumw2();
	h_susy      ->Sumw2();
	// vector of all histos
	vector<TH1D*> h_samples;

	for(size_t i = 0; i < fSamples.size(); ++i){
	        if(!(fSamples[i].type=="susy") && !(fSamples[i].type=="mc")) continue;

		h_samples.push_back(new TH1D(varname+"_"+fSamples[i].name, "", nbins, min, max));
		h_samples[i] -> Sumw2();

		Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);

		if(fVerbose>2) cout << "MakePlot: looping over " << fSamples[i].sname << endl;
		if(fVerbose>2) cout << "           sample has weight " << weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	
		TString variable  = TString::Format("%s>>%s",var.Data(),h_samples[i]->GetName());
		TString selection = TString::Format("(%f) * (%s)",weight,cuts.Data());
		  
		if(fVerbose>2) cout << "+++++ Drawing " << variable  << endl
				    << "\twith cuts: "  << selection << endl;

		int nev = fSamples[i].tree->Draw(variable.Data(),selection.Data(),"goff");
		
		if(fVerbose>2) cout << "\tevents found : "  <<  nev << endl
				    << "\t->Integral() : "  <<  h_samples[i]->Integral() << endl;
		
		if( fSamples[i].type=="mc"        )   h_mc_sum->Add(h_samples[i]);
		else if( fSamples[i].type=="susy" )   h_susy  ->Add(h_samples[i]);
	}
	Float_t  x[nbins], y[nbins];
	for (int i = 1; i <=nbins+1; i++){
	  x[i-1] = h_mc_sum->GetBinLowEdge(i);
	  float s = h_susy  ->Integral(i,nbins+1);
	  float b = h_mc_sum->Integral(i,nbins+1);
	  switch (type){
	  case 0:
	    y[i-1] = s/sqrt(b);
	    break;
	  case 1:
	    y[i-1] = s/sqrt(s+b);
	    break;
	  case 2:
	    y[i-1] = s/b;
	    break;
	  }
	}
	TGraph *sig = new TGraph(nbins+1,x,y);
	sig->SetTitle("");
	sig->GetXaxis()->SetTitle(xtitle);
	sig->SetMarkerStyle(20);
	switch (type){
	case 0:
	  sig->GetYaxis()->SetTitle("S/#sqrt{B}");
	  break;
	case 1:
	  sig->GetYaxis()->SetTitle("S/#sqrt{S+B}");
	  break;
	case 2:
	  sig->GetYaxis()->SetTitle("S/B");
	  break;
	}
	sig->Draw("ACP");
}
//__________________________________________________________________________
void MassPlotter::CompSamples(TString var, TString cuts, TString optcut, bool RemoveLepts, 
		              TString xtitle, const int nbins, const double *bins, bool add_underflow, bool logflag, double scale_factor, bool normalize){

  	CompSamples(fSamples, var, cuts, optcut, RemoveLepts, xtitle, nbins, bins, add_underflow, logflag, scale_factor, normalize);
}
//___________________________________________________________________________________________________________________
void MassPlotter::compSamples(TString var, TString cuts, TString optcut, bool RemoveLepts,
			   TString xtitle, const int nbins, const double min, const double max, bool add_underflow, bool logflag, double scale_factor, bool normalize){

  	CompSamples(fSamples, var, cuts, optcut, RemoveLepts, xtitle, nbins, min, max, add_underflow, logflag, scale_factor, normalize);

}
// _______________________________________________________________________
void MassPlotter::CompSamples(std::vector<sample> Samples, TString var, TString cuts, TString optcut, bool RemoveLepts,
			   TString xtitle, const int nbins, const double min, const double max, bool add_underflow, bool logflag, double scale_factor, bool normalize){

  	double bins[nbins];
  	bins[0] = min;
  	for(int i=1; i<=nbins; i++)
    	bins[i] = min+i*(max-min)/nbins;
  	CompSamples(Samples, var, cuts, optcut, RemoveLepts, xtitle, nbins, bins, add_underflow, logflag, scale_factor, normalize);

}
// __________________________________________________________________________
void MassPlotter::CompSamples(std::vector<sample> Samples, TString var, TString cuts, TString optcut, bool RemoveLepts,
			   TString xtitle, const int nbins, const double *bins, bool add_underflow, bool logflag, double scale_factor, bool normalize){

        TString varname = Util::removeFunnyChar(var.Data());
	TH1D*    histos[2];
	cout << Samples.size() << endl;
	histos[0] = new TH1D(varname+Samples[0].sname  , "", nbins, bins);
	histos[1] = new TH1D(varname+Samples[1].sname  , "", nbins, bins);

	histos[0] -> Sumw2(); 
	histos[1] -> Sumw2();
	histos[0] -> SetLineColor(kBlue); 
	histos[1] -> SetLineColor(kRed); 
	histos[0] -> SetLineWidth(4);
	histos[1] -> SetLineWidth(4);

	// legend
	TLegend* Legend1 = new TLegend(.71,.54,.91,.92);

	for(int i=0; i<2; ++i){
		Double_t weight = scale_factor * Samples[i].xsection * Samples[i].kfact * Samples[i].lumi / (Samples[i].nevents);
		if(fVerbose>2) cout << "GetHisto: looping over " << Samples[i].sname << endl;
		if(fVerbose>2) cout << "           sample has weight " << weight << " and " << Samples[i].tree->GetEntries() << " entries" << endl; 
		// exchange variable to plot with the corresponding one with removed leptons;	
		TString variable; 
		TString selection; 
		if(RemoveLepts && (Samples[i].sname == "DYToLL" || Samples[i].sname=="QCD_2")){
			MapType::iterator iter = RemoveLeptMap.begin();
			iter = RemoveLeptMap.find(var);
			if (iter != RemoveLeptMap.end() ){
				variable = iter -> second;
			} else {cout << "found RemoveLepts==true, but no corresponding map" << endl; variable = var;}
		}else {
			variable = var;
		}
		if(optcut!="_" && Samples[i].sname == "DYToLL"){
			TString newcuts = cuts + "&&" + optcut;  
			if(Samples[i].type!="data") selection      = TString::Format("(%.15f*pileUp.Weight) * (%s)",weight,newcuts.Data());
			else                        selection      = TString::Format("(%.15f) * (%s)"              ,weight,newcuts.Data()); 
		} else {
			if(Samples[i].type!="data") selection      = TString::Format("(%.15f*pileUp.Weight) * (%s)",weight,cuts.Data());
			else                        selection      = TString::Format("(%.15f) * (%s)"              ,weight,cuts.Data()); 
		}

		variable          = TString::Format("%s>>%s",variable.Data(),histos[i]->GetName());
		
			  
		if(fVerbose>2) cout << "+++++ Drawing " << variable  << endl
				    << "\twith cuts: "  << selection << endl;

		int nev = Samples[i].tree->Draw(variable.Data(),selection.Data(),"goff");

		// Add underflow & overflow bins
		// This failed for older ROOT version when the first(last) bin is empty
		// and there are underflow (overflow) events --- must check whether this 
		// is still the case
		if(add_underflow) {
			histos[i]->SetBinContent(1,
				    histos[i]->GetBinContent(0) + histos[i]->GetBinContent(1));
			histos[i]->SetBinError(1,
				  sqrt(histos[i]->GetBinError(0)*histos[i]->GetBinError(0)+
				       histos[i]->GetBinError(1)*histos[i]->GetBinError(1) ));
		}
		histos[i]->SetBinContent(histos[i]->GetNbinsX(),
					    histos[i]->GetBinContent(histos[i]->GetNbinsX()  )+ 
					    histos[i]->GetBinContent(histos[i]->GetNbinsX()+1) );
		histos[i]->SetBinError(histos[i]->GetNbinsX(),
					  sqrt(histos[i]->GetBinError(histos[i]->GetNbinsX()  )*
					       histos[i]->GetBinError(histos[i]->GetNbinsX()  )+
					       histos[i]->GetBinError(histos[i]->GetNbinsX()+1)*
					       histos[i]->GetBinError(histos[i]->GetNbinsX()+1)  ));
		
		if(fVerbose>2) cout << "\tevents found : "  <<  nev << endl
				    << "\t->Integral() : "  <<  histos[i]->Integral() << endl;
	     	Legend1 ->AddEntry(histos[i], Samples[i].sname, "l"); 
	}
	
	plotRatio(histos[0], histos[1], logflag, normalize, varname+"_"+Samples[0].sname+"_"+Samples[1].sname, Legend1 , xtitle, "");
	
	
}

//________________________________________________________________________

void MassPlotter::MakePlot(std::vector<sample> Samples, TString var, TString cuts, int njets, int nleps, TString HLT,
			   TString xtitle, const int nbins, const double min, const double max,
			   bool flip_order, bool logflag, bool composited, bool ratio, 
			   bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow){

  double bins[nbins];
  bins[0] = min;
  for(int i=1; i<=nbins; i++)
    bins[i] = min+i*(max-min)/nbins;
  MakePlot(Samples, var, cuts, njets, nleps, HLT, xtitle, nbins, bins, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow);

}

//________________________________________________________________________

void MassPlotter::MakePlot(std::vector<sample> Samples, TString var, TString cuts, int njets, int nleps, TString HLT,
			   TString xtitle, const int nbins, const double *bins, 
			   bool flip_order, bool logflag, bool composited, bool ratio, 
			   bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow){

        TString varname = Util::removeFunnyChar(var.Data());

	TString nJets = "NJetsIDLoose40";
	nJets += njets < 0 ? ">=" : "==";
	nJets += TString::Format("%d",abs(njets));

	TString  nLeps;
	if     (nleps < 0 )  nLeps = " && (NEles + NMuons) >=";
	else if(nleps >=0  ) nLeps = " && (NEles + NMuons) ==";
	nLeps += TString::Format("%d",abs(nleps));
	if     (nleps ==-10) nLeps = " "; 
	if     (nleps ==-11) nLeps = " && NEles ==1 && NMuons ==0"; 
	if     (nleps ==-13) nLeps = " && NEles ==0 && NMuons ==1"; 

	THStack* h_stack     = new THStack(varname, "");
  	TH1D*    h_data      = new TH1D   (varname+"data"  , "", nbins, bins );
	TH1D*    h_mc_sum    = new TH1D   (varname+"mc_sum", "", nbins, bins );
	TH1D*    h_susy      = new TH1D   (varname+"susy"  , "", nbins, bins );	

	// h_data
	h_data -> Sumw2();
	h_data -> SetMarkerStyle(20);
	h_data -> SetMarkerColor(kBlack);
	h_data -> SetLineColor(kBlack);
	h_data -> SetStats(false);

	h_mc_sum -> SetFillStyle(3004);
	h_mc_sum -> SetFillColor(kBlack);
	h_mc_sum -> SetStats(0);
	
	h_mc_sum    ->Sumw2();
	h_susy      ->Sumw2();


	// vector of all histos
	vector<TH1D*> h_samples;

	// arrays of composited stuff
	TH1D    *h_composited[7];
	//TString  cnames[5] = {"QCD", "W/Z/#gamma production", "Top production", "susy", "data"};
	//int      ccolor[5] = {401, 418, 602, 0, 632};
	TString  cnames[7] = {"QCD", "W+jets", "Z+jets", "Top", "Other", "susy", "data"};
	int      ccolor[7] = { 401,   417,       419,      600,    603,     0,     632};
	for (int i=0; i<7; i++){
	  h_composited[i] = new TH1D(varname+"_"+cnames[i], "", nbins, bins);
	  h_composited[i] -> Sumw2();
	  h_composited[i] -> SetFillColor  (stacked ? ccolor[i] : 0);
	  h_composited[i] -> SetLineColor  (ccolor[i]);
	  if (!stacked) {
	    h_composited[i] -> SetLineWidth(4);
	  }
	  h_composited[i] -> SetMarkerColor(ccolor[i]);
	  h_composited[i] -> SetStats(false);
	}
	if(fisPhoton){
		h_composited[4] -> SetFillColor(kViolet-3); 
		h_composited[4] -> SetLineColor(kViolet-3);
		cnames[4]="Photons";
	}

	// legend
	TLegend* Legend1 = new TLegend(.71,.68,.91,.92);
	//TLegend* Legend1 = new TLegend(.3,.5,.6,.88);

	for(size_t i = 0; i < Samples.size(); ++i){
		h_samples.push_back(new TH1D(varname+"_"+Samples[i].name, "", nbins, bins));
		h_samples[i] -> Sumw2();
		h_samples[i] -> SetFillColor(stacked ? Samples[i].color : 0);
		h_samples[i] -> SetLineColor(Samples[i].color);
		if (!stacked) {
		  h_samples[i] -> SetLineWidth(4);
		}
		h_samples[i] -> SetMarkerColor(Samples[i].color);
		h_samples[i] -> SetStats(false);
		if(Samples[i].type == "susy" ){
			h_samples[i] -> SetLineColor(kBlack);
			h_samples[i] -> SetLineStyle(kDotted);
			h_composited[5] -> SetLineColor(kBlack);
			h_composited[5] -> SetLineStyle(kDotted);
		}
	
		Double_t weight = Samples[i].xsection * Samples[i].kfact * Samples[i].lumi / (Samples[i].nevents*Samples[i].PU_avg_weight);
		if(fVerbose>2) cout << "MakePlot: looping over " << Samples[i].sname << "-----------------------------------" <<  endl;
		if(fVerbose>2) cout << "  +++++++ xsection:    "    << Samples[i].xsection << " k-fact " << Samples[i].kfact << endl;
		if(fVerbose>2) cout << "  +++++++ tot events:  "  << Samples[i].nevents  << " avg pu weight " << Samples[i].PU_avg_weight << endl;
		if(fVerbose>2) cout << "  +++++++ weight:      " << weight << " and " << Samples[i].tree->GetEntries() << " entries" << endl; 
	
		TString variable  = TString::Format("%s>>%s",var.Data(),h_samples[i]->GetName());
		TString theCuts = nJets + nLeps + "&&" + cuts;
		if(Samples[i].type=="data" && HLT!="") theCuts += " &&("+HLT+")"; // triggers for data

		TString selection;
		if(Samples[i].type!="data") selection      = TString::Format("(%.15f*pileUp.Weight) * (%s)",weight,theCuts.Data());
		else                        selection      = TString::Format("(%.15f) * (%s)"              ,weight,theCuts.Data()); 
		  
		if(fVerbose>2) cout << "  +++++++ Drawing      " << variable  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;

		int nev = Samples[i].tree->Draw(variable.Data(),selection.Data(),"goff");

		// Add underflow & overflow bins
		// This failed for older ROOT version when the first(last) bin is empty
		// and there are underflow (overflow) events --- must check whether this 
		// is still the case
		if(add_underflow) {
		  h_samples[i]->SetBinContent(1,
					      h_samples[i]->GetBinContent(0) + h_samples[i]->GetBinContent(1));
		  h_samples[i]->SetBinError(1,
					    sqrt(h_samples[i]->GetBinError(0)*h_samples[i]->GetBinError(0)+
						 h_samples[i]->GetBinError(1)*h_samples[i]->GetBinError(1) ));
		}
		h_samples[i]->SetBinContent(h_samples[i]->GetNbinsX(),
					    h_samples[i]->GetBinContent(h_samples[i]->GetNbinsX()  )+ 
					    h_samples[i]->GetBinContent(h_samples[i]->GetNbinsX()+1) );
		h_samples[i]->SetBinError(h_samples[i]->GetNbinsX(),
					  sqrt(h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()  )*
					       h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()  )+
					       h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()+1)*
					       h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()+1)  ));
		
		if(fVerbose>2) cout << "  +++++++ MC   events : "  <<  nev << endl;

		/// event count with errors
		TH1F * clone = (TH1F*)h_samples[i]->Clone();
		clone->Rebin(clone->GetNbinsX());
		if(fVerbose>2) cout << "  +++++++ Events:       " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;
		
		if (Samples[i].sname.Contains("QCD")) {
		  h_composited[0]->Add(h_samples[i]);
		}
		else if(Samples[i].sname.Contains("Top") ){
		  h_composited[3]->Add(h_samples[i]);
		}
		else if(Samples[i].sname=="Wtolnu") {
		  h_composited[1]->Add(h_samples[i]);
		}
		else if(Samples[i].sname=="DY") {
		  h_composited[2]->Add(h_samples[i]);
		}
		else if(Samples[i].sname=="VV" || Samples[i].sname=="VGamma" || Samples[i].sname=="Photons") {
		  h_composited[4]->Add(h_samples[i]);
		}
		else if(Samples[i].type.Contains("susy")){
		  h_composited[5]->Add(h_samples[i]);
		  cnames[5] = Samples[i].sname;
		}
		else if(Samples[i].type.Contains("data")){
		  h_composited[6]->Add(h_samples[i]);
		}

	}

	if (composited){
	  if(flip_order){
	    for (int i=5; i>=0; --i){
	      if (!stacked)  {
	        h_composited[i]->Scale(1./h_composited[i]->Integral());
	        h_composited[i] -> SetMinimum(0.00005);
	        //if (fVerbose>2) cout << " h_composited[" << i<< "]->Integral = " << h_composited[i]->Integral()  << endl;
	      }
	      if( i!=5 || !overlaySUSY) h_stack  -> Add(h_composited[i]);
	      if( i==5)                 h_susy   -> Add(h_composited[i]);
	      else        	        h_mc_sum -> Add(h_composited[i]);
	      if( i==5 && overlaySUSY)
	        Legend1 ->AddEntry(h_composited[i], cnames[i]  , "f");  
	      else
	        Legend1 ->AddEntry(h_composited[i], cnames[i], stacked ? "f" : "l");
	    }
	  }
	  else{
	    for (int i=0; i<6; ++i){
	      if (!stacked)  {
	        h_composited[i]->Scale(1./h_composited[i]->Integral());
	        h_composited[i] -> SetMinimum(0.00005);
	        //if (fVerbose>2) cout << " h_composited[" << i<< "]->Integral = " << h_composited[i]->Integral()  << endl;
	      }
	      if( i!=5 || !overlaySUSY) h_stack  -> Add(h_composited[i]);
	      if( i==5)                 h_susy   -> Add(h_composited[i]);
	      else        	      h_mc_sum -> Add(h_composited[i]);
	      if( i==5 && overlaySUSY){
	        if(h_composited[i]->GetEntries()>0) Legend1 ->AddEntry(h_composited[i], cnames[i] + (overlayScale ? overlayScale==1? "":
		  					      TString::Format(" x %.0f",overlayScale) :
		  					      " scaled to data")   , "f");  
	      }else{
	        if(h_composited[i]->GetEntries()>0) Legend1 ->AddEntry(h_composited[i], cnames[i], stacked ? "f" : "l");
	      }
	    }
	  }
	  h_data->Add(h_composited[6]);
	  if(h_data->Integral()>0 && stacked){
	    Legend1     ->AddEntry(h_data, "data", "l");
	  }
	  if(fVerbose > 2 && composited) {
		           cout << "------------------------------------"                << endl
	                        << "QCD Integral:      " << h_composited[0]->Integral()  << endl
			        << "W+jets Integral:   " << h_composited[1]->Integral()  << endl
			        << "Z+jets Integral:   " << h_composited[2]->Integral()  << endl
			        << "Top Integral:      " << h_composited[3]->Integral()  << endl
			        << "Other Integral:    " << h_composited[4]->Integral()  << endl
				<< "TOTAL BG:          " << h_composited[0]->Integral()+h_composited[1]->Integral()+h_composited[2]->Integral()+h_composited[3]->Integral()+h_composited[4]->Integral() <<endl
			        << "SUSY:              " << h_composited[5]->Integral()  << endl
			        << "Data:              " << h_composited[6]->Integral()  << endl;

			   //Same, but with errors
			   string OverallSamples[8] = {"QCD","W+jets","Z+Jets","Top","Other","SUSY","Data","TOTAL BG"};
			   for(int os=0; os<8; os++){
			     string tSample = OverallSamples[os];
			     if(tSample!="TOTAL BG"){
			       TH1F* clone = (TH1F*) h_composited[os]->Clone();
			       clone->Rebin(clone->GetNbinsX());
			       if(fVerbose>2) cout << tSample << ": " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;
			     }
			     else{
			       TH1F* clone = (TH1F*) h_composited[0]->Clone();
			       clone->Add(h_composited[1]);
			       clone->Add(h_composited[2]);
                               clone->Add(h_composited[3]);
			       clone->Add(h_composited[4]);
			       clone->Rebin(clone->GetNbinsX());
                               if(fVerbose>2) cout << tSample << ": " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;
			     }
			   }
			  // save nevents for W background prediction
			  if(nleps==-11){
			  	fWpred.QCD_bg_e   = h_composited[0]->Integral();	
				fWpred.W_bg_e     = h_composited[1]->Integral();
				fWpred.Z_bg_e     = h_composited[2]->Integral();
				fWpred.Top_bg_e   = h_composited[3]->Integral();
				fWpred.Other_bg_e = h_composited[4]->Integral();
			  }else if(nleps==-13){
			  	fWpred.QCD_bg_mu   = h_composited[0]->Integral();	
				fWpred.W_bg_mu     = h_composited[1]->Integral();
				fWpred.Z_bg_mu     = h_composited[2]->Integral();
				fWpred.Top_bg_mu   = h_composited[3]->Integral();
				fWpred.Other_bg_mu = h_composited[4]->Integral();
			  } 
			  if(var=="Znunu.RecoOSee_mll"){
			  	fZpred.QCD_bg_e   = h_composited[0]->Integral();
				fZpred.W_bg_e     = h_composited[1]->Integral();
				fZpred.Top_bg_e   = h_composited[3]->Integral();
				fZpred.Other_bg_e = h_composited[4]->Integral();
			  }else if(var=="Znunu.RecoOSmumu_mll"){
			  	fZpred.QCD_bg_mu  = h_composited[0]->Integral();
				fZpred.W_bg_mu    = h_composited[1]->Integral();
				fZpred.Top_bg_mu  = h_composited[3]->Integral();
				fZpred.Other_bg_mu= h_composited[4]->Integral();
			  }
	  }
	}
	else{
	  for(int i=0; i<Samples.size(); ++i){
	    if (!stacked) {
	      h_samples[i]->Scale(1./h_samples[i]->Integral());
	      h_samples[i] -> SetMinimum(0.0001);
	    }
	    if((Samples[i].sname=="QCD" )){
	      h_samples[i] -> SetMinimum(0.01);
	      h_stack->Add(h_samples[i]);
	      h_mc_sum->Add(h_samples[i]);
	      Legend1 ->AddEntry(h_samples[i], Samples[i].sname, "f");
	    }
	    if(Samples[i].sname!="QCD" && Samples[i].type!="data"){
	      if(Samples[i].type=="susy") h_susy->Add(h_samples[i]);
	      if(Samples[i].type!="susy" || !overlaySUSY) h_stack->Add(h_samples[i]);
	      if(Samples[i].type!="susy") h_mc_sum->Add(h_samples[i]);
	      if(Samples[i].type=="susy" && overlaySUSY){
	         if(h_samples[i]->Integral()>0) Legend1 ->AddEntry(h_samples[i], Samples[i].sname + (overlayScale ? 
								       TString::Format(" x %.0f",overlayScale) :
								       " scaled to data")   , "f") ;
	      }else{
		 if(h_samples[i]->Integral()>0) Legend1 ->AddEntry(h_samples[i], Samples[i].sname, stacked ? "f" : "l");
	      }	
	    }
	    if(Samples[i].type == "data"){
	      h_data -> Add(h_samples[i]);
	    }
	    
	  }
	  if(h_data->Integral()>0){
	    Legend1     ->AddEntry(h_data, "data", "l");
	  }
	}

	TString ytitle = "Events";

	TString oname = var+ "_CUT_";
	oname += njets < 0 ? TString::Format("ge%dJets",abs(njets)) : TString::Format("%dJets",abs(njets));
	oname += nleps == 0 ? "_0Lep_" : nleps == -10 ? "_anyLep_" : nleps == -11 ? "_1Ele_" : nleps == -13 ? "_1Muo_" :
			     nleps < 0 ? TString::Format("_ge%dLeps_",abs(nleps)) : TString::Format("_%dLeps_",abs(nleps));
	oname += cuts;
	oname.ReplaceAll(">=" ,".ge");
	oname.ReplaceAll("<=" ,".le");
	oname.ReplaceAll(">" ,".gt");
	oname.ReplaceAll("<" ,".lt");
	oname.ReplaceAll("==",".eq");
	oname.ReplaceAll("!=",".ne");
	oname.ReplaceAll("&&","_");
	oname.ReplaceAll("||","_");
	oname.ReplaceAll("misc.","");
	oname.ReplaceAll("LeadingJPt","J1Pt");
	oname.ReplaceAll("SecondJPt","J2Pt");
	oname.ReplaceAll("LeptConfig","LepCfg");
	oname.ReplaceAll("Vectorsumpt","VSPT");
	oname.ReplaceAll("EcalDeadCellBEFlag","BEFlg");
	oname.ReplaceAll("HBHENoiseFlag","HBHEFlg");
	oname.ReplaceAll("CSCTightHaloID","CSCFlg");
	oname.ReplaceAll("NJetsIDLoose40","NJIDLoose40");
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
        TString outname = Util::removeFunnyChar(oname.Data());
	outname.ReplaceAll("NMuons.eq0_muo0.lv.Pt.lt10_NEles.eq0_ele0.lv.Pt.lt10","noLepPt10");
	outname.ReplaceAll("NMuons.gt0_muo0.lv.Pt.gt10_NEles.gt0_ele0.lv.Pt.gt10","gt1lepPt10");

	if(!stacked) {
	  outname = outname + (flip_order ? "_flipped" : "") + (logflag ? "_log" : "") + "_shape";
	  printHisto(h_stack, h_data, h_mc_sum, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nleps, false);

	}
	else if (!overlaySUSY){
	  outname = outname + (flip_order ? "_flipped" : "") + (logflag ? "_log" : "") + (composited ? "_comp" : "");
	  printHisto(h_stack, h_data, h_mc_sum, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nleps);
	  if (ratio) 
	    plotRatioStack(h_stack,  h_mc_sum, h_data,  logflag, false, outname, Legend1, xtitle, ytitle, njets, nleps);
	}
	else {
	  outname =  outname + (flip_order ? "_flipped" : "") + (logflag ? "_log" : "") + (composited ? "_comp" : "") + "_overlay";	  
	  printHisto(h_stack, h_data, h_mc_sum, h_susy, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nleps, overlayScale);
	  if (ratio) 
	    plotRatioStack(h_stack,  h_mc_sum, h_data, h_susy,  logflag, false, outname, Legend1, xtitle, ytitle, njets, nleps, overlayScale);

	}
	
// 	for(int i=0; i<h_samples.size(); ++i){
// 		delete h_samples[i];
// 	}
// 	delete h_mc_sum;
// 	delete Legend1;
// 	delete h_data;
// 	delete h_susy;
	
}




//__________________________________________________________________________
void MassPlotter::abcd_MT2(TString var, TString basecut, TString upper_cut, TString lower_cut, const int nbins,const double min, const double max, double fit_min, double fit_max){

  double bins[nbins];
  bins[0] = min;
  for(int i=1; i<=nbins; i++)
    bins[i] = min+i*(max-min)/nbins;
  ABCD_MT2(var, basecut, upper_cut, lower_cut, nbins, bins, fit_min, fit_max);
  
}

//__________________________________________________________________________
void MassPlotter::ABCD_MT2(TString var, TString basecut, TString upper_cut, TString lower_cut, const int nbins, const double *bins, double fit_min, double fit_max){
  TH2D *h_ABCD_MT2_qcd           = new TH2D("ABCD_MT2__qcd"           , "",  100, 0, 600, 180, 0, TMath::Pi());  h_ABCD_MT2_qcd          ->Sumw2();
  TH2D *h_ABCD_MT2_susy          = new TH2D("ABCD_MT2__susy"          , "",  100, 0, 600, 180, 0, TMath::Pi());  h_ABCD_MT2_susy         ->Sumw2();
  TH1D *h_ABCD_lower_y_band_susy = new TH1D("ABCD_lower_y_band__susy" , "",  nbins, bins);		        h_ABCD_lower_y_band_susy->Sumw2();
  TH1D *h_ABCD_upper_y_band_susy = new TH1D("ABCD_upper_y_band__susy" , "",  nbins, bins);		        h_ABCD_upper_y_band_susy->Sumw2();
  TH1D *h_ABCD_lower_y_band_data = new TH1D("ABCD_lower_y_band__data" , "",  nbins, bins);		        h_ABCD_lower_y_band_data->Sumw2();
  TH1D *h_ABCD_upper_y_band_data = new TH1D("ABCD_upper_y_band__data" , "",  nbins, bins);		        h_ABCD_upper_y_band_data->Sumw2();
  TH1D *h_ABCD_lower_y_band_qcd  = new TH1D("ABCD_lower_y_band__qcd"  , "",  nbins, bins);		        h_ABCD_lower_y_band_qcd ->Sumw2();
  TH1D *h_ABCD_upper_y_band_qcd  = new TH1D("ABCD_upper_y_band__qcd"  , "",  nbins, bins);		        h_ABCD_upper_y_band_qcd ->Sumw2();
  TH1D *h_ABCD_lower_y_band_mc   = new TH1D("ABCD_lower_y_band__mc"   , "",  nbins, bins);		        h_ABCD_lower_y_band_mc  ->Sumw2();
  TH1D *h_ABCD_upper_y_band_mc   = new TH1D("ABCD_upper_y_band__mc"   , "",  nbins, bins);		        h_ABCD_upper_y_band_mc  ->Sumw2();
  TH1D *h_ABCD_lower_y_band_sub  = new TH1D("ABCD_lower_y_band__sub"  , "",  nbins, bins);		        h_ABCD_lower_y_band_sub ->Sumw2();
  TH1D *h_ABCD_upper_y_band_sub  = new TH1D("ABCD_upper_y_band__sub"  , "",  nbins, bins);		        h_ABCD_upper_y_band_sub ->Sumw2();
  TH1D *ratio_data               = new TH1D("ratio__data"             , "",  nbins, bins);		        ratio_data              ->Sumw2();
  TH1D *ratio_qcd                = new TH1D("ratio__qcd"              , "",  nbins, bins);		        ratio_qcd               ->Sumw2();
  TH1D *ratio_sub                = new TH1D("ratio__sub"              , "",  nbins, bins);		        ratio_sub               ->Sumw2();
  
  for(size_t i = 0; i < fSamples.size(); ++i){
    
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    
    Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
    if(fVerbose>2) cout << "ABCD: looping over " << fSamples[i].name << endl;
    if(fVerbose>2) cout << "      sample has weight " << weight << endl; 
    
//     TEntryList *elist = new TEntryList("elist","elist");
//     int nev = fSamples[i].tree->Draw(">>elist",basecut.Data(),"entrylist");
//     fSamples[i].tree->SetEntryList(elist);
    
    TString  weights   = (fSamples[i].type!="data" ?  TString::Format("(%.15f*pileUp.Weight)",weight) : TString::Format("(%.15f)",weight));
    TString selection = TString::Format("(%s) * (%s)"      ,weights.Data(),basecut.Data());
    TString sel_up    = TString::Format("(%s) * (%s && %s)",weights.Data(),basecut.Data(),upper_cut.Data());
    TString sel_lo    = TString::Format("(%s) * (%s && %s)",weights.Data(),basecut.Data(),lower_cut.Data());

    // Remove jets intentionaly!!!
    // RemoveAndRecalcMT2() has to be called before any selection to make the internal changes effective!
    TString removeJet = "1";
    //TString removeJet = "SmearAndRecalcMT2()>-5.";
    //TString removeJet = "RemoveAndRecalcMT2(0,0.000001)>-5.";
    //TString removeJet = "RemoveAndRecalcMT2(1,0.005)>-5.";
    TString selection_QCD = TString::Format("(%s) * (%s && %s)"      ,weights.Data(),removeJet.Data(),basecut.Data());
    TString sel_up_QCD    = TString::Format("(%s) * (%s && %s && %s)",weights.Data(),removeJet.Data(),basecut.Data(),upper_cut.Data());
    TString sel_lo_QCD    = TString::Format("(%s) * (%s && %s && %s)",weights.Data(),removeJet.Data(),basecut.Data(),lower_cut.Data());

    TString variable;
    TString var2 = "misc.MT2";
    if (fSamples[i].type == "data"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_data->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_data->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
      if(fVerbose>2) cout << "\tData lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_data->Integral() + h_ABCD_lower_y_band_data->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\tData upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_data->Integral() + h_ABCD_upper_y_band_data->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_120to170"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      variable  = TString::Format("%s:misc.MT2>>+%s",var.Data(),h_ABCD_MT2_qcd->GetName());
      TString sel_u    = /*sel_up_QCD; /*/TString::Format("(%s) * (%s && misc.MT2<80&& %s)",weights.Data(),basecut.Data(),upper_cut.Data());
      TString sel_l    = /*sel_lo_QCD; /*/TString::Format("(%s) * (%s && misc.MT2<80&& %s)",weights.Data(),basecut.Data(),lower_cut.Data());
      int nev2d= fSamples[i].tree->Draw(variable,selection_QCD,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_qcd->Integral() + h_ABCD_lower_y_band_qcd->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_qcd->Integral() + h_ABCD_upper_y_band_qcd->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_170to300"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      variable  = TString::Format("%s:misc.MT2>>+%s",var.Data(),h_ABCD_MT2_qcd->GetName());
      TString sel_u    = /*sel_up_QCD; /*/TString::Format("(%s) * (%s && misc.MT2<140&& %s)",weights.Data(),basecut.Data(),upper_cut.Data());
      TString sel_l    = /*sel_lo_QCD; /*/TString::Format("(%s) * (%s && misc.MT2<170&& %s)",weights.Data(),basecut.Data(),lower_cut.Data());
      int nev2d= fSamples[i].tree->Draw(variable,selection_QCD,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_qcd->Integral() + h_ABCD_lower_y_band_qcd->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_qcd->Integral() + h_ABCD_upper_y_band_qcd->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_300to470"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      variable  = TString::Format("%s:misc.MT2>>+%s",var.Data(),h_ABCD_MT2_qcd->GetName());
      TString sel_u    = /*sel_up_QCD; /*/TString::Format("(%s) * (%s && misc.MT2<200&& %s)",weights.Data(),basecut.Data(),upper_cut.Data());
      TString sel_l    = /*sel_lo_QCD; /*/TString::Format("(%s) * (%s && %s)",weights.Data(),basecut.Data(),lower_cut.Data());
      int nev2d= fSamples[i].tree->Draw(variable,selection_QCD,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_qcd->Integral() + h_ABCD_lower_y_band_qcd->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_qcd->Integral() + h_ABCD_upper_y_band_qcd->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "mc" && fSamples[i].sname == "QCD"){
      variable  = TString::Format("%s:misc.MT2>>+%s",var.Data(),h_ABCD_MT2_qcd->GetName());
      int nev2d= fSamples[i].tree->Draw(variable,selection_QCD,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo_QCD,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up_QCD,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_qcd->Integral() + h_ABCD_lower_y_band_qcd->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_qcd->Integral() + h_ABCD_upper_y_band_qcd->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "mc"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_mc->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_mc->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_mc->Integral() + h_ABCD_lower_y_band_mc->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_mc->Integral() + h_ABCD_upper_y_band_mc->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "susy"){
      variable  = TString::Format("%s:misc.MT2>>+%s",var.Data(),h_ABCD_MT2_susy->GetName());
      int nev2d= fSamples[i].tree->Draw(variable,selection,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_susy->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_susy->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_susy->Integral() + h_ABCD_lower_y_band_susy->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_susy->Integral() + h_ABCD_upper_y_band_susy->GetBinContent(nbins+1) << endl;
    }
  }
  
  
  h_ABCD_lower_y_band_sub->Add(h_ABCD_lower_y_band_data,h_ABCD_lower_y_band_mc,1,-1);
  h_ABCD_upper_y_band_sub->Add(h_ABCD_upper_y_band_data,h_ABCD_upper_y_band_mc,1,-1);

  for (int i= 0; i<=nbins; i++){
    if (h_ABCD_lower_y_band_sub->GetBinContent(i)<0)  h_ABCD_lower_y_band_sub->SetBinContent(i, 0.);
    if (h_ABCD_upper_y_band_sub->GetBinContent(i)<0)  h_ABCD_upper_y_band_sub->SetBinContent(i, 0.);
  }

  ratio_qcd ->Divide(h_ABCD_upper_y_band_qcd ,h_ABCD_lower_y_band_qcd );
  ratio_data->Divide(h_ABCD_upper_y_band_data,h_ABCD_lower_y_band_data);
  ratio_sub ->Divide(h_ABCD_upper_y_band_sub ,h_ABCD_lower_y_band_sub );
  
  h_ABCD_lower_y_band_data->SetLineColor(1);  h_ABCD_upper_y_band_data->SetLineColor(2);
  h_ABCD_lower_y_band_qcd ->SetLineColor(1);  h_ABCD_upper_y_band_qcd ->SetLineColor(2);
  h_ABCD_lower_y_band_sub ->SetLineColor(1);  h_ABCD_upper_y_band_sub ->SetLineColor(2);
  h_ABCD_lower_y_band_sub ->SetLineStyle(3);  h_ABCD_upper_y_band_sub ->SetLineStyle(3);   ratio_sub->SetLineStyle(3);
  ratio_qcd ->SetMarkerStyle(24);   ratio_qcd ->SetMarkerSize(0.8);
  ratio_data->SetMarkerStyle(24);   ratio_data->SetMarkerSize(0.8);
  ratio_sub ->SetMarkerStyle(24);   ratio_sub ->SetMarkerSize(0.8);
  ratio_qcd ->SetXTitle("M_{T2}");   ratio_qcd ->SetYTitle("ratio");
  ratio_sub ->SetXTitle("M_{T2}");   ratio_sub ->SetYTitle("ratio");

  //TF1 *f_qcd = new TF1("f_qcd","pol0(0)+expo(1)",bins[0], bins[nbins]);   f_qcd->SetLineColor(8);
  //TF1 *f_sub = new TF1("f_sub","pol0(0)+expo(1)",bins[0], bins[nbins]);   f_sub->SetLineColor(8);
  //TF1 *f_qcd = new TF1("f_qcd","exp([0]+300.*[1])+expo(0)",bins[0], bins[nbins]);   f_qcd->SetLineColor(8);
  //TF1 *f_sub = new TF1("f_sub","exp([0]+300.*[1])+expo(0)",bins[0], bins[nbins]);   f_sub->SetLineColor(8);
  TF1 *f_qcd1 = new TF1("f_qcd1","expo(0)",bins[0], bins[nbins]);   f_qcd1->SetLineColor(8);
  TF1 *f_sub1 = new TF1("f_sub1","expo(0)",bins[0], bins[nbins]);   f_sub1->SetLineColor(8);
  TF1 *f_qcd2 = new TF1("f_qcd2","exp([0]+200.*[1])+expo(0)",bins[0], bins[nbins]);   f_qcd2->SetLineColor(9);
  TF1 *f_sub2 = new TF1("f_sub2","exp([0]+200.*[1])+expo(0)",bins[0], bins[nbins]);   f_sub2->SetLineColor(8);
  TF1 *f_qcd3 = new TF1("f_qcd3","expo(0)+[2]"              ,bins[0], bins[nbins]);   f_qcd3->SetLineColor(50);

  gStyle->SetPalette(1);
  gStyle->SetOptFit (0);
  gStyle->SetOptStat(0);

  TCanvas *can = new TCanvas("can", "2d qcd", 0, 0, 900, 700);
  can->SetLogz(1);
  h_ABCD_MT2_qcd->SetMinimum(0.0001);
  h_ABCD_MT2_qcd->Draw("colz");
  TCanvas *ccan = new TCanvas("ccan", "2d susy", 0, 0, 900, 700);
  ccan->SetLogz(1);
  h_ABCD_MT2_susy->SetMinimum(0.0001);
  h_ABCD_MT2_susy->Draw("colz");
  TCanvas *can2 = new TCanvas("can2", "qcd", 0, 0, 900, 700);
  can2->SetLogy(1);
  h_ABCD_lower_y_band_qcd ->Draw();  
  h_ABCD_upper_y_band_qcd ->Draw("same");
  TCanvas *can3 = new TCanvas("can3", "data w/ and w/o substraction", 0, 0, 900, 700);
  can3->SetLogy(1);
  h_ABCD_lower_y_band_data->Draw();  
  h_ABCD_upper_y_band_data->Draw("same");
  h_ABCD_lower_y_band_sub ->Draw("same");  
  h_ABCD_upper_y_band_sub ->Draw("same");
  TCanvas *can4 = new TCanvas("can4", "ratio qcd", 0, 0, 900, 700);
  can4->SetLogy(1);
  ratio_qcd->Draw("E1");
  //f_qcd->FixParameter(0,0.006);
  ratio_qcd->SetTitle("QCD");
  ratio_qcd->Fit("f_qcd1","0","",fit_min,fit_max);
  f_qcd1->Draw("same");
  f_qcd2->SetParameter(0,f_qcd1->GetParameter(0));
  f_qcd2->SetParameter(1,f_qcd1->GetParameter(1));
  f_qcd2->Draw("same");
  f_qcd3->SetParameter(0,f_qcd1->GetParameter(0));
  f_qcd3->SetParameter(1,f_qcd1->GetParameter(1));
  f_qcd3->SetParameter(2,0.002);
  ratio_qcd->Fit("f_qcd3","0","",fit_min,bins[nbins]);
  f_qcd3->Draw("same");
  TLine *l_qcd_min = new TLine(fit_min,f_qcd1->GetYaxis()->GetXmin(),fit_min,f_qcd1->GetMaximum()); l_qcd_min->SetLineStyle(2); l_qcd_min->Draw("same");
  TLine *l_qcd_max = new TLine(fit_max,f_qcd1->GetYaxis()->GetXmin(),fit_max,f_qcd1->GetMaximum()); l_qcd_max->SetLineStyle(2); l_qcd_max->Draw("same");

  TCanvas *can5 = new TCanvas("can5", "ratio data", 0, 0, 900, 700);
  can5->SetLogy(1);
  ratio_sub ->Draw("PE1");
  ratio_data->Draw("sameE1");
  //float p0 = f_qcd->GetParameter(0);
  //f_sub->FixParameter(0,0.006);
  //ratio_data->Fit("f_sub1","0","",fit_min,fit_max);  // fit w/o subtracting background
  ratio_sub->Fit("f_sub1","0","",fit_min,fit_max);   // background subtracted fit
  f_sub1->Draw("same");
  f_sub2->SetParameter(0,f_sub1->GetParameter(0));
  f_sub2->SetParameter(1,f_sub1->GetParameter(1));
  f_sub2->Draw("same");
  TLine *l_sub_min = new TLine(fit_min,f_sub1->GetYaxis()->GetXmin(),fit_min,f_sub1->GetMaximum()); l_sub_min->SetLineStyle(2); l_sub_min->Draw("same");
  TLine *l_sub_max = new TLine(fit_max,f_sub1->GetYaxis()->GetXmin(),fit_max,f_sub1->GetMaximum()); l_sub_max->SetLineStyle(2); l_sub_max->Draw("same");
  TLegend *leg = new TLegend(.55,.5,.85,.75);
  leg->SetFillColor(0);
  leg->AddEntry(ratio_data,"data","le");
  leg->AddEntry(ratio_sub ,"data (EWK subtracted)","le");
  leg->Draw("same");

  PrintABCDPredictions(var, basecut, upper_cut, lower_cut, f_qcd1,f_sub1,f_qcd3);
  
}

//_________________________________________________________________________________
void MassPlotter::PrintABCDPredictions(TString var, TString basecut, TString upper_cut, TString lower_cut, TF1* func_qcd, TF1* func_sub, TF1* func_qcd_model){  
  int nbins=180;
  float min=100., max=1000.;
  TH1D *h_pred_lower_y_band_data = new TH1D("pred_lower_y_band__data" , "", nbins,min,max);		        h_pred_lower_y_band_data->Sumw2();
  TH1D *h_pred_upper_y_band_data = new TH1D("pred_upper_y_band__data" , "", nbins,min,max);		        h_pred_upper_y_band_data->Sumw2();
  TH1D *h_pred_lower_y_band_qcd  = new TH1D("pred_lower_y_band__qcd"  , "", nbins,min,max);		        h_pred_lower_y_band_qcd ->Sumw2();
  TH1D *h_pred_upper_y_band_qcd  = new TH1D("pred_upper_y_band__qcd"  , "", nbins,min,max);		        h_pred_upper_y_band_qcd ->Sumw2();
  TH1D *h_pred_lower_y_band_mc   = new TH1D("pred_lower_y_band__mc"   , "", nbins,min,max);		        h_pred_lower_y_band_mc  ->Sumw2();
  TH1D *h_pred_upper_y_band_mc   = new TH1D("pred_upper_y_band__mc"   , "", nbins,min,max);		        h_pred_upper_y_band_mc  ->Sumw2();
  TH1D *h_pred_lower_y_band_sub  = new TH1D("pred_lower_y_band__sub"  , "", nbins,min,max);		        h_pred_lower_y_band_sub ->Sumw2();
  TH1D *h_pred_upper_y_band_sub  = new TH1D("pred_upper_y_band__sub"  , "", nbins,min,max);		        h_pred_upper_y_band_sub ->Sumw2();
  TH1D *h_pred_lower_y_band_susy = new TH1D("pred_lower_y_band__susy" , "", nbins,min,max);		        h_pred_lower_y_band_susy->Sumw2();
  TH1D *h_pred_upper_y_band_susy = new TH1D("pred_upper_y_band__susy" , "", nbins,min,max);		        h_pred_upper_y_band_susy->Sumw2();
  //TF1 *f_qcd = new TF1("f_qcd","pol0(0)+expo(1)",100,1000);
  //TF1 *f_sub = new TF1("f_sub","pol0(0)+expo(1)",100,1000);
//   TF1 *f_qcd = new TF1("f_qcd","exp([0]+300.*[1])+expo(0)",100,1000);
//   TF1 *f_sub = new TF1("f_sub","exp([0]+300.*[1])+expo(0)",100,1000);
  TF1 *f_qcd = new TF1("f_qcd","expo(0)",min,max);
  TF1 *f_sub = new TF1("f_sub","expo(0)",min,max);
  f_qcd->SetParameters(func_qcd->GetParameter(0),func_qcd->GetParameter(1));//,func_qcd->GetParameter(2));
  f_sub->SetParameters(func_sub->GetParameter(0),func_sub->GetParameter(1));//,func_sub->GetParameter(2));
  TF1 *f_qcd2 = new TF1("f_qcd2","exp([0]+200.*[1])+expo(0)",min,max);
  TF1 *f_sub2 = new TF1("f_sub2","exp([0]+200.*[1])+expo(0)",min,max);
  f_qcd2->SetParameters(func_qcd->GetParameter(0),func_qcd->GetParameter(1));//,func_qcd->GetParameter(2));
  f_sub2->SetParameters(func_sub->GetParameter(0),func_sub->GetParameter(1));//,func_sub->GetParameter(2));
  TF1 *f_qcd_model = new TF1("f_qcd_model","expo(0)+[2]",min,max);
  f_qcd_model->SetParameters(func_qcd_model->GetParameter(0),func_qcd_model->GetParameter(1),func_qcd_model->GetParameter(2));

  for(size_t i = 0; i < fSamples.size(); ++i){
   
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    
    Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
    if(fVerbose>2) cout << "ABCD: looping over " << fSamples[i].name << endl;
    if(fVerbose>2) cout << "      sample has weight " << weight << endl; 

    TString weights   = fSamples[i].type!="data" ?  TString::Format("(%.15f*pileUp.Weight)",weight) : TString::Format("(%.15f)",weight);
    TString selection = TString::Format("(%s) * (%s)"      ,weights.Data(),basecut.Data());
    TString sel_up    = TString::Format("(%s) * (%s && %s)",weights.Data(),basecut.Data(),upper_cut.Data());
    TString sel_lo    = TString::Format("(%s) * (%s && %s)",weights.Data(),basecut.Data(),lower_cut.Data());
    if(fVerbose>2) cout << "      sel_lo = " << weights << endl; 
//     TString selection = TString::Format("(%f) * (%s)"      ,weight,basecut.Data());
//     TString sel_up    = TString::Format("(%f) * (%s && %s)",weight,basecut.Data(),upper_cut.Data());
//     TString sel_lo    = TString::Format("(%f) * (%s && %s)",weight,basecut.Data(),lower_cut.Data());

    TString variable;
    TString var2 = "misc.MT2";
    if (fSamples[i].type == "data"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_data->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_data->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_120to170"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      TString sel_u    = TString::Format("(%s) * (%s && misc.MT2<80&& %s)",weights.Data(),basecut.Data(),upper_cut.Data());
      TString sel_l    = TString::Format("(%s) * (%s && misc.MT2<80&& %s)",weights.Data(),basecut.Data(),lower_cut.Data());
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_170to300"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      TString sel_u    = TString::Format("(%s) * (%s && misc.MT2<140&& %s)",weights.Data(),basecut.Data(),upper_cut.Data());
      TString sel_l    = TString::Format("(%s) * (%s && misc.MT2<170&& %s)",weights.Data(),basecut.Data(),lower_cut.Data());
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_300to470"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      TString sel_u    = TString::Format("(%s) * (%s && misc.MT2<200&& %s)",weights.Data(),basecut.Data(),upper_cut.Data());
      TString sel_l    = TString::Format("(%s) * (%s && %s)",weights.Data(),basecut.Data(),lower_cut.Data());
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
    }
    else if (fSamples[i].type == "mc" && fSamples[i].sname == "QCD"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
    }
    else if (fSamples[i].type == "mc"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_mc->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_mc->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
    }
    else if (fSamples[i].type == "susy"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_susy->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_susy->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
    }
  }
  
  
  h_pred_lower_y_band_sub ->Add(h_pred_lower_y_band_data,h_pred_lower_y_band_mc,1,-1);
  h_pred_upper_y_band_sub ->Add(h_pred_upper_y_band_data,h_pred_upper_y_band_mc,1,-1);
  h_pred_lower_y_band_susy->Add(h_pred_lower_y_band_sub); // adding susy contamination to subtracted data

  for (int i= 0; i<=nbins+1; i++){
    if (h_pred_lower_y_band_sub->GetBinContent(i)<0)  h_pred_lower_y_band_sub->SetBinContent(i, 0.);
    if (h_pred_upper_y_band_sub->GetBinContent(i)<0)  h_pred_upper_y_band_sub->SetBinContent(i, 0.);
  }

  TH1D *h_pred_lower_y_band_data_2   = (TH1D*)h_pred_lower_y_band_data->Clone("h_pred_lower_y_band_data_2");         h_pred_lower_y_band_data_2   ->Sumw2();
  TH1D *h_pred_lower_y_band_qcd_2    = (TH1D*)h_pred_lower_y_band_qcd ->Clone("h_pred_lower_y_band_qcd_2" );         h_pred_lower_y_band_qcd_2    ->Sumw2();
  TH1D *h_pred_lower_y_band_sub_2    = (TH1D*)h_pred_lower_y_band_sub ->Clone("h_pred_lower_y_band_sub_2" );         h_pred_lower_y_band_sub_2    ->Sumw2();
  TH1D *h_pred_lower_y_band_data_c   = (TH1D*)h_pred_lower_y_band_data->Clone("h_pred_lower_y_band_data_c");         h_pred_lower_y_band_data_c   ->Sumw2();
  TH1D *h_pred_lower_y_band_qcd_c    = (TH1D*)h_pred_lower_y_band_qcd ->Clone("h_pred_lower_y_band_qcd_c" );         h_pred_lower_y_band_qcd_c    ->Sumw2();
  TH1D *h_pred_lower_y_band_sub_c    = (TH1D*)h_pred_lower_y_band_sub ->Clone("h_pred_lower_y_band_sub_c" );         h_pred_lower_y_band_sub_c    ->Sumw2();
  TH1D *h_pred_lower_y_band_susy_c   = (TH1D*)h_pred_lower_y_band_susy->Clone("h_pred_lower_y_band_susy_c" );        h_pred_lower_y_band_susy_c  ->Sumw2();
  TH1D *h_pred_lower_y_band_data_2_c = (TH1D*)h_pred_lower_y_band_data->Clone("h_pred_lower_y_band_data_2_c");       h_pred_lower_y_band_data_2_c ->Sumw2();
  TH1D *h_pred_lower_y_band_qcd_2_c  = (TH1D*)h_pred_lower_y_band_qcd ->Clone("h_pred_lower_y_band_qcd_2_c" );       h_pred_lower_y_band_qcd_2_c  ->Sumw2();
  TH1D *h_pred_lower_y_band_sub_2_c  = (TH1D*)h_pred_lower_y_band_sub ->Clone("h_pred_lower_y_band_sub_2_c" );       h_pred_lower_y_band_sub_2_c  ->Sumw2();
  TH1D *h_pred_lower_y_band_qcd_model= (TH1D*)h_pred_lower_y_band_qcd ->Clone("h_pred_lower_y_band_qcd_model");      h_pred_lower_y_band_qcd_model->Sumw2();

  // from qcd fit (exp+const) in whole range (50-500)
  h_pred_lower_y_band_qcd_model->Multiply(f_qcd_model);
  // from qcd fit

  h_pred_lower_y_band_data  ->Multiply(f_qcd);
  h_pred_lower_y_band_qcd   ->Multiply(f_qcd);
  h_pred_lower_y_band_sub   ->Multiply(f_qcd);
  // from data fit
  h_pred_lower_y_band_data_2->Multiply(f_sub);
  h_pred_lower_y_band_qcd_2 ->Multiply(f_sub);
  h_pred_lower_y_band_sub_2 ->Multiply(f_sub);
  // from qcd fit (exp+const)
  h_pred_lower_y_band_data_c  ->Multiply(f_qcd2);
  h_pred_lower_y_band_qcd_c   ->Multiply(f_qcd2);
  h_pred_lower_y_band_sub_c   ->Multiply(f_qcd2);
  // from data fit (exp+const)
  h_pred_lower_y_band_data_2_c->Multiply(f_sub2);
  h_pred_lower_y_band_qcd_2_c ->Multiply(f_sub2);
  h_pred_lower_y_band_sub_2_c ->Multiply(f_sub2);
  // with susy contamination from data fit (exp+const)
  h_pred_lower_y_band_susy  ->Multiply(f_sub);
  h_pred_lower_y_band_susy_c->Multiply(f_sub2);

  cout << "SUSY prediction:" << endl;
  printEstimation(h_pred_upper_y_band_susy,h_pred_upper_y_band_susy, nbins, min, max);
  cout << "QCD prediction:" << endl;
  printEstimation(h_pred_upper_y_band_qcd,h_pred_upper_y_band_qcd, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from QCD model (exp+const):" << endl;
  printEstimation(h_pred_lower_y_band_qcd_model,h_pred_lower_y_band_qcd_model, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from QCD (fit to QCD):" << endl;
  printEstimation(h_pred_lower_y_band_qcd,h_pred_lower_y_band_qcd_c, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from QCD (fit to data):" << endl;
  printEstimation(h_pred_lower_y_band_qcd_2,h_pred_lower_y_band_qcd_2_c, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from data (fit to QCD):" << endl;
  printEstimation(h_pred_lower_y_band_data,h_pred_lower_y_band_data_c, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from data (fit to data):" << endl;
  printEstimation(h_pred_lower_y_band_data_2,h_pred_lower_y_band_data_2_c, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from EWK subtracted data (fit to QCD):" << endl;
  printEstimation(h_pred_lower_y_band_sub,h_pred_lower_y_band_sub_c, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from EWK subtracted data (fit to data):" << endl;
  printEstimation(h_pred_lower_y_band_sub_2,h_pred_lower_y_band_sub_2_c, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from QCD with susy contamination (fit to data):" << endl;
  printEstimation(h_pred_lower_y_band_susy,h_pred_lower_y_band_susy_c, nbins, min, max);
  cout << "------------------------------------" << endl;

}

void MassPlotter::printEstimation(TH1D* h_pred, TH1D* h_pred_c, int nbins, float min, float max){
  float delta = (max-min)/((float)nbins);
  double yield[5][2], error[5][2];
  yield[0][0] = Util::IntegralAndError(h_pred,  (int)((100.-min)/delta)+1, nbins+1, error[0][0]);
  cout << "\tMT2>100 = " << yield[0][0]  << " +/- " << error[0][0];
  yield[0][1] = Util::IntegralAndError(h_pred_c,(int)((100.-min)/delta)+1, nbins+1, error[0][1]);
  cout << "\tMT2>100 = " << yield[0][1]  << " +/- " << error[0][1] << endl;
  yield[1][0] = Util::IntegralAndError(h_pred,  (int)((150.-min)/delta)+1, nbins+1, error[1][0]);
  cout << "\tMT2>150 = " << yield[1][0]  << " +/- " << error[1][0];
  yield[1][1] = Util::IntegralAndError(h_pred_c,(int)((150.-min)/delta)+1, nbins+1, error[1][1]);
  cout << "\tMT2>150 = " << yield[1][1]  << " +/- " << error[1][1] << endl;
  yield[2][0] = Util::IntegralAndError(h_pred,  (int)((200.-min)/delta)+1, nbins+1, error[2][0]);
  cout << "\tMT2>200 = " << yield[2][0]  << " +/- " << error[2][0];
  yield[2][1] = Util::IntegralAndError(h_pred_c,(int)((200.-min)/delta)+1, nbins+1, error[2][1]);
  cout << "\tMT2>200 = " << yield[2][1]  << " +/- " << error[2][1] << endl;
  yield[3][0] = Util::IntegralAndError(h_pred,  (int)((250.-min)/delta)+1, nbins+1, error[3][0]);
  cout << "\tMT2>250 = " << yield[3][0]  << " +/- " << error[3][0];
  yield[3][1] = Util::IntegralAndError(h_pred_c,(int)((250.-min)/delta)+1, nbins+1, error[3][1]);
  cout << "\tMT2>250 = " << yield[3][1]  << " +/- " << error[3][1] << endl;
  yield[4][0] = Util::IntegralAndError(h_pred,  (int)((300.-min)/delta)+1, nbins+1, error[4][0]);
  cout << "\tMT2>300 = " << yield[4][0]  << " +/- " << error[4][0];
  yield[4][1] = Util::IntegralAndError(h_pred_c,(int)((300.-min)/delta)+1, nbins+1, error[4][1]);
  cout << "\tMT2>300 = " << yield[4][1]  << " +/- " << error[4][1] << endl;

  for (int i=0; i<2; i++){
    cout << "yield_" << (i==0 ? "opt" : "pes") << " = {" << endl;
    for (int j=0; j<5; j++){
      cout << yield[j][i] << (j<4 ? ", " : "\n}\n");
    }
    cout << "error_" << (i==0 ? "opt" : "pes") << " = {" << endl;
    for (int j=0; j<5; j++){
      cout << error[j][i] << (j<4 ? ", " : "\n}\n");
    }
  }

}

//_________________________________________________________________________________
void MassPlotter::plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets, int nleps){
	// define canvas and pads 
	TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
	TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	
	TCanvas* c1 = new TCanvas(name,"", 20,100,1000,700);
	cout << name << endl;
	c1 -> cd();
	
	float border = 0.3;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad(name+"_plotpad",  "Pad containing the overlay plot", 0.00, border, 1.00, 1.00, 0, 0);
 	p_plot->SetBottomMargin(0);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad(name+"_ratiopad", "Pad containing the ratio",        0.00, 0.00, 1.00, border, 0, 0);
 	p_ratio->SetTopMargin(0);
 	p_ratio->SetBottomMargin(0.35);
 	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	if(logflag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	hstack->SetMaximum(max);
//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);

	hstack->SetMinimum(0.02);
	hstack ->Draw("hist");
	//h2     ->Draw("sameEX0");
        h2     ->Draw("sameE1"); //LEO MOD

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.1);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 0.9, ytitle);
		
	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.05);
	TString text = njets < 0 ? TString::Format("#geq %d Jets",abs(njets)) : TString::Format("%d Jets",abs(njets));
	text += nleps == 1 ? ", 1 Lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());

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
 
	TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");	

	h_ratio->GetXaxis()->SetTitleSize(scale * h1->GetXaxis()->GetTitleSize());
	h_ratio->GetYaxis()->SetTitleSize(scale * h1->GetYaxis()->GetTitleSize());
	h_ratio->GetXaxis()->SetLabelSize(scale * h1->GetXaxis()->GetLabelSize());
	h_ratio->GetYaxis()->SetLabelSize(scale * h1->GetYaxis()->GetLabelSize());
	h_ratio->GetXaxis()->SetTickLength(scale * h1->GetXaxis()->GetTickLength());
	h_ratio->GetYaxis()->SetTickLength(h1->GetYaxis()->GetTickLength());
	h_ratio->SetLineWidth(2);
	//h_ratio->SetFillColor(kBlue);//leo
	//h_ratio->SetLineColor(kBlue);
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);
	//h_ratio ->SetMarkerSize(0.1);

 	h_ratio ->Divide(h2, h1);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	h_ratio ->DrawCopy("E"/*2"*/);
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	// ratio title
	lat.SetTextAlign(33);
	lat.SetTextColor(1);
	lat.SetTextSize(0.15);
	lat.SetNDC(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 1.0, "data / MC");
	
	// x axis title
	lat.SetTextSize(0.2);
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.23 : 0.16;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
	if(fSave) Util::Print(c1, save, fOutputDir, fOutputFile);	

// 	delete h1;
// 	delete h2;
// 	delete h_ratio;
// 	delete p_plot;
// 	delete p_ratio;
// 	delete c1;

}
//_________________________________________________________________________________
void MassPlotter::plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, TH1* h3, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets, int nleps, float overlayScale){
  //LEO TRUE USE THIS

	// define canvas and pads 
	TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
	TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	h1->SetMarkerStyle(20);
	h2->SetMarkerStyle(20);

	
	//TCanvas* c1 = new TCanvas(name,"", 20,100,1000,700);
	TCanvas* c1 = new TCanvas("c_ratio","",0,0,600,600 /*37, 60,636,670*/);
	c1->SetFrameLineWidth(1);
	c1 -> cd();
	
	float border = 0.2;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad(name+"_plotpad",  "Pad containing the overlay plot", 0,0.211838,1,1 /*0.00, border, 1.00, 1.00, 0, 0*/);
 	//p_plot->SetBottomMargin(0.05);
	//p_plot->SetTopMargin(0.09);
	//p_plot->SetLeftMargin(0.1669107);
      	//p_plot->SetRightMargin(0.02);
	p_plot->SetLeftMargin(0.131579);
	p_plot->SetRightMargin(0.08);
	p_plot->SetTopMargin(0.06895515);
	p_plot->SetBottomMargin(0.07206074);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad(name+"_ratiopad", "Pad containing the ratio",   0,0.01863354,0.9967105,0.2189441/*     0.00, 0.05, 1.00, border, 0, 0*/);
 	//p_ratio->SetTopMargin(0.03);
 	//p_ratio->SetBottomMargin(0.05/*5*/);
	//p_ratio->SetRightMargin(0.02);
	p_ratio->SetLeftMargin(0.1336634);	
	p_ratio->SetRightMargin(0.075);
	p_ratio->SetTopMargin(0.06976745);
	p_ratio->SetBottomMargin(0.2790698);

	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	if(logflag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2.5*max;
	else max = 1.5*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	hstack->SetMaximum(max);
//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);
	stringstream yTitle;
	if(fEventsPerGeV){
		if(fabs(h1_orig->GetBinWidth(1) -h1_orig->GetBinWidth(h1_orig->GetNbinsX()-1))<0.01){
		double binwidth = h1_orig->GetBinWidth(1);
		yTitle.precision(3);
		yTitle << ytitle.Data();
		yTitle << " / ";
		yTitle << binwidth;
		yTitle << " GeV";
		} else{
		cout << h1_orig->GetBinWidth(1) << " " << h1_orig->GetBinWidth(h1_orig->GetNbinsX()-1) << endl;
		}
	}else{
		yTitle << ytitle.Data();
	}
	hstack->GetYaxis()->SetTitle(yTitle.str().c_str());
	hstack->GetYaxis()->SetLabelSize(0.05);
	hstack->GetYaxis()->SetTitleSize(0.05);
	hstack->GetYaxis()->SetTitleOffset(1.3);

	
	//MT2_bSel[0]->SetTitleSize(0.03);
	///MT2_bSel[0]->SetTitleOffset(1.);
	hstack->SetMinimum(0.02);
	hstack->Draw("hist");
	h2    ->Draw("sameE");
	h3->Scale(overlayScale ? overlayScale : h2->Integral() / h3->Integral());
	h3->SetFillColor(0);
	h3->SetLineStyle(kDotted);
	h3->Draw("samehist");

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.03);
	TString text ="";
	text = fMT2Analysis?  "M_{T2} Analysis                                          ":"";
	text +=fMT2bAnalysis? "M_{T2}b Analysis                                         ":"";
	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
	text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";
	TitleBox.DrawLatex(0.13,0.943,text.Data());

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
 
	TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");	
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);

 	h_ratio ->Divide(h2, h1);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	//MC with errors
	TH1D*h_ratio_mc = (TH1D*)h1_orig->Clone("h1_copy");
	h_ratio_mc->Divide(h1);
	h_ratio_mc->GetYaxis()->SetRangeUser(0,2);
	h_ratio_mc->GetXaxis()->SetLabelSize( 0.);
	h_ratio_mc->GetYaxis()->SetTitle("Data / MC");
	h_ratio_mc->GetXaxis()->SetTitle(xtitle);
	h_ratio_mc->GetXaxis()->SetTitleSize(0.2);
	h_ratio_mc->GetXaxis()->SetTitleOffset(0.5);
	h_ratio_mc->GetYaxis()->SetLabelSize(0.19);
	h_ratio_mc->GetXaxis()->SetTickLength(0.09);
	h_ratio_mc	->GetYaxis()->SetTitleSize(0.18);
	h_ratio_mc->GetYaxis()->SetTitleOffset(0.36);

	
	h_ratio_mc->SetFillStyle(3001);
	h_ratio_mc->Draw("E2");
	h_ratio ->DrawCopy("Esame");//LEO MOD
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
	if(fSave)Util::Print(c1, save, fOutputDir, fOutputFile);	


}
//_________________________________________________________________________________
void MassPlotter::plotRatio(TH1* h1_orig, TH1* h2_orig, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle){
	// define canvas and pads
	TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
	TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	
	TCanvas* c1 = new TCanvas("c1","", 20,100,1000,700);
	c1 -> cd();
	
	float border = 0.3;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad("plotpad",  "Pad containing the overlay plot", 0.00, border, 1.00, 1.00, 0, 0);
 	p_plot->SetBottomMargin(0);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio",        0.00, 0.00, 1.00, border, 0, 0);
 	p_ratio->SetTopMargin(0);
 	p_ratio->SetBottomMargin(0.35);
 	p_ratio->Draw();
 
	gPad->SetFillStyle(0);
	// draw overlay plot
 	p_plot ->cd();

	if(logflag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);
	h2     ->Draw();
	h1     ->Draw("same");
	
	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.07);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 0.9, ytitle);
		
 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
 
	TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");	

	h_ratio->GetXaxis()->SetTitleSize(scale * h1->GetXaxis()->GetTitleSize());
	h_ratio->GetYaxis()->SetTitleSize(scale * h1->GetYaxis()->GetTitleSize());
	h_ratio->GetXaxis()->SetLabelSize(scale * h1->GetXaxis()->GetLabelSize());
	h_ratio->GetYaxis()->SetLabelSize(scale * h1->GetYaxis()->GetLabelSize());
	h_ratio->GetXaxis()->SetTickLength(scale * h1->GetXaxis()->GetTickLength());
	h_ratio->GetYaxis()->SetTickLength(h1->GetYaxis()->GetTickLength());
	h_ratio->SetLineWidth(2);
	h_ratio->SetFillColor(kBlue);
	h_ratio->SetLineColor(kBlue);
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);
	h_ratio ->SetMarkerSize(0.1);

 	h_ratio ->Divide(h1, h2);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	gPad->SetFillStyle(0);
	h_ratio ->DrawCopy("E2");
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	// ratio title
	lat.SetTextAlign(33);
	lat.SetTextColor(1);
	lat.SetTextSize(0.07);
	lat.SetNDC(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 1.0, "ratio");
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.23 : 0.16;
	lat.DrawLatex(0.9, ypos, xtitle);
	//gPad->SetLogy(1);
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
	if(fSave)Util::PrintEPS(c1, save, fOutputDir);	

// 	delete h1;
//	delete h2;
//	delete h_ratio;
// 	delete p_plot;
// 	delete p_ratio;
// 	delete c1;

}
//___________________________________________________________________________
void MassPlotter::printHisto(TH1* h_orig, TString canvname,  Option_t *drawopt, bool logflag){
	TH1D* h = (TH1D*)h_orig->Clone("h1_copy");
	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	gStyle->SetPalette(1);
	if(logflag) gPad->SetLogy(1);
	h->SetMinimum(0.01);
	gPad->SetRightMargin(0.15);
	gPad->SetLogz();
	h->SetStats(false);
	h->DrawCopy(drawopt);
	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	h->Write();
	delete col;
	delete h;

}

//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TString canvname, Option_t *drawopt, bool logflag){
	
	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) gPad->SetLogy(1);
	h->Draw(drawopt);
	gPad->RedrawAxis();
	col ->Update();
	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	delete col;

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets, int nleps, bool stacked){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h        -> SetMinimum(stacked ? 0.05 : 0.0001);
		h_mc_sum -> SetMinimum(0.05);
		h_data   -> SetMinimum(0.05);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h_mc_sum->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	if(stacked){
	  h_data  ->SetMaximum(max);
	  h_mc_sum->SetMaximum(max);
	  h       ->SetMaximum(max);
	}

	h->Draw(stacked ? drawopt : "histnostack");
	//h_mc_sum -> Draw("same, E2");
	if(h_data->Integral()>0 && stacked) {
		h_data       ->Draw("same");
	}
	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.05);
	TString text = njets < 0 ? TString::Format("#geq %d Jets",abs(njets)) : TString::Format("%d Jets",abs(njets));
	text += nleps == 1 ? ", 1 Lepton" : "";
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
	lat.DrawLatex(0.03, 0.9, ytitle);
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.1 : 0.06;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	if(fSave)Util::PrintEPS(col, canvname, fOutputDir);
// 	delete col;

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, TH1* h_susy, TLegend* legend,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets, int nleps, float overlayScale){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 600, 600);
	col->SetRightMargin(0.08);
	col->SetLeftMargin(0.18);
	col->SetTopMargin(0.07);
	col->SetBottomMargin(0.17);

	TLegend *leg = (TLegend*) legend->Clone("leg");

	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h        -> SetMinimum(0.05);
		h_mc_sum -> SetMinimum(0.05);
		h_data   -> SetMinimum(0.05);
		h_susy   -> SetMinimum(0.05);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h_mc_sum->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h_mc_sum->SetMaximum(max);
	h       ->SetMaximum(max);

	h->Draw(drawopt);
	//h_mc_sum -> Draw("same, E2");
	if(h_data->Integral()>0) {
		h_data       ->Draw("sameE");
	}
	h_susy->Scale(overlayScale ? overlayScale : h_data->Integral() / h_susy->Integral());
	h_susy->SetLineStyle(kDotted);
	h_susy->SetFillColor(0);
	h_susy->Draw("samehist");
	if(leg != NULL ){
		leg -> SetY1NDC(0.68);
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.0305);
	TString text = fMT2Analysis? "M_{T2} Analysis        ":"";
	text +=fMT2bAnalysis? "M_{T2}b Analysis     ":"";
	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
	text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";
	TitleBox.DrawLatex(0.18,0.943,text.Data());


	h->GetXaxis()->SetTitle(xtitle);
	h->GetXaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetTitleOffset(1.1);
	
	stringstream yTitle;
	if(fEventsPerGeV){
		if(fabs(h_data->GetBinWidth(1) -h_data->GetBinWidth(h_data->GetNbinsX()-1))<0.01){
		double binwidth = h_data->GetBinWidth(1);
		yTitle.precision(3);
		yTitle << ytitle.Data();
		yTitle << " / ";
		yTitle << binwidth;
		yTitle << " GeV";
		} else{
		cout << h_data->GetBinWidth(1) << " " << h_data->GetBinWidth(h_data->GetNbinsX()-1) << endl;
		}
	}else{
		yTitle << ytitle.Data();
	}
	h->GetYaxis()->SetTitle(yTitle.str().c_str());
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleOffset(1.47);

	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	if(fSave)Util::PrintEPS(col, canvname, fOutputDir);
// 	delete col;

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag , TString xtitle, TString ytitle){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h->SetMinimum(0.01);
		h_data->SetMinimum(0.01);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h     ->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h       ->SetMaximum(max);

	h->Draw(drawopt);
	if(h_data->Integral()>0) {
		h_data  ->Draw("same, E1");
	}
	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.07);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.03, 0.9, ytitle);
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.1 : 0.06;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
	col ->Update();
	if(fSave) Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	delete col;

}
//____________________________________________________________________________
void MassPlotter::loadSamples(const char* filename){
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
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "ShapeName\t%s", StringValue);
			s.shapename = TString(StringValue);

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


//_________________________________________________________________________
void MassPlotter::TauContamination(int sample_index, Long64_t nevents, int flag){
  float sys = 0.05; //systematic uncertanty on the fake rate and efficiency
  float IsoCut = 0.5;//cut on the isolation 0.5:VLoose, 1.5:Loose...
  float maxDPhiTauMET = 11.5;//To cut on the DeltaPhi between the tau and MET

  setFlags(flag); 

  TH1::SetDefaultSumw2();
  TString  cnames[NumberOfSamples] = {"QCD", "Wjets", "Zjets", "Top", "MC", "data"};
  int      ccolor[NumberOfSamples] = {  401,     417,     419,   600,  500,  632};
  TString varname = "MT2";
  for (int i=0; i<NumberOfSamples; i++){
    MT2[i] = new TH1F(varname+"_"+cnames[i], "", 150, 0, 750);
    MT2[i] -> SetFillColor  (ccolor[i]);
    MT2[i] -> SetLineColor  (ccolor[i]);
    MT2[i] -> SetLineWidth  (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }

  MT2[5] -> SetMarkerStyle(20);
  MT2[5] -> SetMarkerColor(kBlack);
  MT2[5] -> SetLineColor(kBlack);
  
  MT2[4] -> SetFillStyle(3004);
  MT2[4] -> SetFillColor(kBlack);

  TH2F* hEtaPtJet           = new TH2F("hEtaPtJet",           "hEtaPtJet",          1,0,1000.0,1,0,2.5);//tauable jets in MC in QCD dominated region (MT2 < 70)
  TH2F* hEtaPtTauFake       = new TH2F("hEtaPtTauFake",       "hEtaPtTauFake",      1,0,1000.0,1,0,2.5);//taus in MC in QCD dominated region  (MT2 < 70)
  TH2F* hEtaPtJetData       = new TH2F("hEtaPtJetData",       "hEtaPtJetData",      1,0,1000.0,1,0,2.5);//tauable jets in data in QCD dominated region (MT2 < 70)
  TH2F* hEtaPtTauFakeData   = new TH2F("hEtaPtTauFakeData",   "hEtaPtTauFakeData",  1,0,1000.0,1,0,2.5);//taus in data in QCD dominated region  (MT2 < 70)

  TH1F* hMT2TauFound        = new TH1F("hMT2TauFound",   "hMT2TauFound",   NumberOfMT2Bins - 1, xMT2bin);//Rec/Id taus
  TH1F* hMT2JetSig          = new TH1F("hMT2JetSig",     "hMT2JetSig",     NumberOfMT2Bins - 1, xMT2bin);//tauable jets in signal region (MC)
  TH1F* hMT2TauSig          = new TH1F("hMT2TauSig",     "hMT2TauSig",     NumberOfMT2Bins - 1, xMT2bin);//taus in the signal region (MC)
  TH1F* hMT2JetSigData      = new TH1F("hMT2JetSigData", "hMT2JetSigData", NumberOfMT2Bins - 1, xMT2bin);//tauable jets in signal region (data)
  TH1F* hMT2TauSigData      = new TH1F("hMT2TauSigData", "hMT2TauSigData", NumberOfMT2Bins - 1, xMT2bin);//taus in the signal region (data)
  TH1F* hRatioTight         = new TH1F("hRatioTight",       "hRatioTight",       NumberOfMT2Bins - 1, xMT2bin);//Used for the final rescaling the estimation from relaxed to appplied cut
  TH1F* hRatioRelaxed       = new TH1F("hRatioRelaxed",       "hRatioRelaxed",       NumberOfMT2Bins - 1, xMT2bin);//Used for the final rescaling the estimation from relaxed to appplied cut

  TH1F* hNAllTaus           = new TH1F("hNAllTaus",      "hNAllTaus",      NumberOfMT2Bins - 1, xMT2bin);//All of the generated hadronic taus
  TH1F* hNAllTauJets        = new TH1F("hNAllTauJets",   "hNAllTauJets",   NumberOfMT2Bins - 1, xMT2bin);//All of the generated hadronic taus matched with a rec jet


  for(int i = 0; i < fSamples.size(); i++){

    if(sample_index != -1 && i != sample_index)
      continue;

    int data = 0;
    sample Sample = fSamples[i];
    
    if(Sample.type == "data")
      data = 1;

    if(Sample.type == "susy") 
      continue;
     
    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  Sample.tree->GetEntries();

    float Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents * Sample.PU_avg_weight );

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    for (Long64_t jentry=0; jentry<min(nentries, nevents);jentry++) {
      Sample.tree->GetEntry(jentry); 
      if ( fVerbose>2 && jentry % 100000 == 0 ){  cout << "+++ Proccessing event " << jentry << endl;}
 
      float weight = Weight;
      if(data == 0)
	weight *= fMT2tree->pileUp.Weight;
      else
	weight = 1.0;

      //For 2011 needs to be updated for 2012
// 	  if (data == 1){
// 	    if(!((fMT2tree->trigger.HLT_HT440_v2 ==1 && fMT2tree->misc.Run<161216)                                                               ||
// 		 (fMT2tree->trigger.HLT_HT450_v2 ==1 && fMT2tree->misc.Run>=161216 && fMT2tree->misc.Run< 163269)                                ||
// 		 (fMT2tree->trigger.HLT_HT500_v3 ==1 && fMT2tree->misc.Run>=163269 && fMT2tree->misc.Run<=163869)                                ||
// 		 (fMT2tree->trigger.HLT_HT500_v4 ==1 && fMT2tree->misc.Run>=165088 && fMT2tree->misc.Run< 165970)                                ||
// 		 (fMT2tree->trigger.HLT_HT550_v5 ==1 && fMT2tree->misc.Run>=165970 && fMT2tree->misc.Run<=167043 && fMT2tree->misc.Run!=166346)  ||
// 		 (fMT2tree->trigger.HLT_HT550_v6 ==1 && fMT2tree->misc.Run==166346)                                                              ||
// 		 (fMT2tree->trigger.HLT_HT550_v7 ==1 && fMT2tree->misc.Run>=167078 && fMT2tree->misc.Run< 170249)                                ||
// 		 (fMT2tree->trigger.HLT_HT550_v8 ==1 && fMT2tree->misc.Run>=170249 && fMT2tree->misc.Run< 173236)                                ||
// 		 (fMT2tree->trigger.HLT_HT600_v1 ==1 && fMT2tree->misc.Run>=173236 && fMT2tree->misc.Run< 178420)                                || 
// 		 (fMT2tree->trigger.HLT_HT650_v4 ==1 && fMT2tree->misc.Run>=178420 && fMT2tree->misc.Run< 179959)))
// 	      continue;
// 	  }


      if (!(fMT2tree->misc.MET             >  30                                   &&
	    fMT2tree->misc.HT              >  750                                  &&
	    (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<10)                  &&
	    (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<10)                  &&
	    fMT2tree->misc.Jet0Pass        == 1                                    &&
	    fMT2tree->misc.Jet1Pass        == 1                                    &&
	    fMT2tree->misc.SecondJPt       >  100                                  &&
	    fMT2tree->misc.PassJetID       == 1                                    &&
	    fMT2tree->misc.Vectorsumpt     <  70                                   &&
            fMT2tree->misc.HBHENoiseFlagIso == 0                                   &&
	    fMT2tree->misc.CSCTightHaloID   == 0                                   &&
            fMT2tree->misc.CrazyHCAL        == 0))
	  continue;

 	  if(High && fMT2tree->misc.HT             <  950)
 	    continue;

 	  if(Low && fMT2tree->misc.HT              >  950)
	    continue;

	  if(mT2){
	    if(!( fMT2tree->misc.LeadingJPt      >  100                                  &&
//  		  fMT2tree->misc.MinMetJetDPhi   >  0.3                                  &&
		  fMT2tree->NJetsIDLoose40       >= 3 
		  ))
	      continue;
	  }else{
	    if(!( fMT2tree->NJetsIDLoose40       >= 4                                    &&
// 		  fMT2tree->NBJets               >  0                                    &&
//  		  fMT2tree->NBJetsHE             >  0                                    &&
//  		  fMT2tree->misc.MinMetJetDPhi4  >  0.3                                  &&
		  fMT2tree->misc.LeadingJPt      >  150                                  
 		  ))
 	    continue;
	    
	    if(!((BTagRelaxed == 0 && fMT2tree->NBJets > 0) || (BTagRelaxed == 1 && fMT2tree->NBJetsHE > 0)))
	      continue;	
	  }

	
	   
	  float METPhi = fMT2tree->misc.METPhi; //To be used if a cone close to met is used to look for taus
	  
   for(int t = 0; t < fMT2tree->NTaus; t++){
 	if(fMT2tree->tau[t].MT > 100 || fMT2tree->tau[t].Isolation < IsoCut)
 	  continue;

	float minDeltaR = 100.0;

 	float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[t].lv.Phi(), fMT2tree->misc.METPhi);

	if(DeltaPhi > maxDPhiTauMET)
	  continue;

	  if(data == 1){
	    
	    MT2[5]->Fill(fMT2tree->misc.MT2, weight);//data
	    
	  }else{
	    MT2[4]->Fill(fMT2tree->misc.MT2, weight);
	    
	    if(Sample.sname == "Top")  
	      MT2[3]->Fill(fMT2tree->misc.MT2, weight);
	    else
	      if(Sample.sname == "DY")	
		MT2[2]->Fill(fMT2tree->misc.MT2, weight);
	      else
		if(Sample.name == "WJetsToLNu")
		  MT2[1]->Fill(fMT2tree->misc.MT2, weight);
		else
		  if(Sample.sname == "QCD")
		    MT2[0]->Fill(fMT2tree->misc.MT2, weight);
	  }}

   for(int j = 0; j < fMT2tree->NGenLepts; j++){
     if(abs(fMT2tree->genlept[j].ID) != 15)
       continue;
	
	    int HadTau = 1;
	    
	    for(int i = 0; i < fMT2tree->NGenLepts; i++){
	      if((abs(fMT2tree->genlept[i].ID) == 11 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID) || 
		 (abs(fMT2tree->genlept[i].ID) == 13 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID)){
		HadTau = 0;
	      }
	    }

	    if(HadTau == 0)
	      continue;
        
	if((mT2 && fMT2tree->misc.MinMetJetDPhi   >  0.3) || (mT2b && fMT2tree->NBJets >  0) || (mT2b && fMT2tree->misc.MinMetJetDPhi4   >  0.3 && fMT2tree->NBJets >  0))
	  hRatioTight->Fill(fMT2tree->misc.MT2, weight);
	else
	  hRatioRelaxed->Fill(fMT2tree->misc.MT2, weight);	  
	
	hNAllTaus->Fill(fMT2tree->misc.MT2, weight);
	

	//To find the jet matched to this gen had tau
	int jetIndex = -1;
	float minDR = 100.0;

	float genleptEta = fMT2tree->genlept[j].lv.Eta();

	float genleptPhi = fMT2tree->genlept[j].lv.Phi();

	for(int k = 0; k < fMT2tree->NJets; k++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->jet[k].lv.Eta(), genleptPhi, fMT2tree->jet[k].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    jetIndex = k;
	  }
	}	

	if(minDR < 0.1){
	  hNAllTauJets->Fill(fMT2tree->misc.MT2, weight);
	}

	jetIndex = -1;
	minDR = 100.0;

	//To find the rec/id tau matched to this gen had tau
	for(int k = 0; k < fMT2tree->NTaus; k++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->tau[k].lv.Eta(), genleptPhi, fMT2tree->tau[k].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    jetIndex = k;
	  }
	}

	float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[jetIndex].lv.Phi(), METPhi);

	if(DeltaPhi > maxDPhiTauMET)
	  continue;
	
	if(minDR < 0.1){
	  if(fMT2tree->tau[jetIndex].Isolation > IsoCut && fMT2tree->tau[jetIndex].MT < 100)
	    hMT2TauFound->Fill(fMT2tree->misc.MT2, weight);
	}
      }
      
      //QCD dominated region, use these events to find the fake rate
      if(fMT2tree->misc.MT2 < 70){
	if(data == 1){
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    
	    if(fMT2tree->jet[j].isTauMatch < 0 )// || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	
	    hEtaPtJetData->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);

	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation > IsoCut && fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT < 100){  
	      hEtaPtTauFakeData->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	    }
	  }
	}else{
	  for(int j = 0; j < fMT2tree->NJets; j++){

	    if(fMT2tree->jet[j].isTauMatch < 0 )// || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;	    

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;

	    hEtaPtJet->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	    
	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation > IsoCut && fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT < 100){ 
	      
	      hEtaPtTauFake->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	    }
	  }
	}
      }
      //Fill The number of Rec/Id taus and tauable jets in the signal region
      if(fMT2tree->misc.MT2 >= MT2Min){
	if(data == 1){
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    if(fMT2tree->jet[j].isTauMatch < 0)// || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	    hMT2JetSigData->Fill(fMT2tree->misc.MT2, weight);
	  }
	
	  for(int j = 0; j < fMT2tree->NTaus; j++){

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	    
	    if(fMT2tree->tau[j].MT < 100 && fMT2tree->tau[j].Isolation > IsoCut){

	      hMT2TauSigData->Fill(fMT2tree->misc.MT2, weight);
	    }
	  }
	}else{
	  
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    if(fMT2tree->jet[j].isTauMatch < 0)// || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100) 
	      continue;

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	
	    hMT2JetSig->Fill(fMT2tree->misc.MT2, weight);
	  }
	
	  for(int j = 0; j < fMT2tree->NTaus; j++){

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	
	    if(fMT2tree->tau[j].MT < 100 && fMT2tree->tau[j].Isolation > IsoCut){
	      hMT2TauSig->Fill(fMT2tree->misc.MT2, weight);

	      float minDR = 100.0;
	      
	      float tauEta = fMT2tree->tau[j].lv.Eta();

	      float tauPhi = fMT2tree->tau[j].lv.Phi();
	    }
	  }
	}
      }
    }//for (Long64_t jentry=0; jentry<min(nentries, nevents);jentry++) 
  }// for(int i = 0; i < fSamples.size(); i++){


 for(int j = 1; j < NumberOfMT2Bins; j++){
   float id = hMT2TauFound  ->GetBinContent(j);
   float jet= hNAllTauJets->GetBinContent(j);
   float gen= hNAllTaus   ->GetBinContent(j);
   float EffTotal = -1.0;
   if(gen != 0) EffTotal = id/gen;
   float Eff = -1.0;
   if(jet != 0) Eff = id/jet;
   float Acc = -1.0;
   if(gen != 0) Acc = jet/gen;

   cout<<"$"<<xMT2bin[j-1]<<"-"<<xMT2bin[j]<<" Eff = "<<EffTotal<<", Eff * Acc = "<<Eff<<" * "<<Acc<<endl;
 }

 int p = 0;
 hMT2TauFound->Divide(hNAllTaus);
 Cout(p++, hMT2TauFound);

 //To add the syatematic uncertainty
 for(int j = 1; j < NumberOfMT2Bins; j++){
   float val = hMT2TauFound->GetBinContent(j);
   float Err = hMT2TauFound->GetBinError(j);
   hMT2TauFound->SetBinError(j, sqrt(Err * Err + sys * sys * val * val));
 }

 Cout(p++, hEtaPtTauFake);
 Cout(p++, hEtaPtJet);
 float TauFake  = hEtaPtTauFake->Integral();
  float Jet      = hEtaPtJet    ->Integral();
  float Fake     = TauFake/Jet;
  cout<<" Fake "<<Fake<<endl;
  hEtaPtTauFake->Divide(hEtaPtJet);
  Cout(p++, hEtaPtTauFake);

  Cout(p++, hMT2JetSig);

  hMT2JetSig->Scale(Fake);

  Cout(p++, hMT2JetSig);
  Cout(p++, hMT2TauSig);
  hMT2JetSig->Add(hMT2TauSig, -1);
  hMT2JetSig->Scale(-1.0);
  Cout(p++, hMT2JetSig);

  hMT2JetSig->Divide(hMT2TauFound);
  Cout(p++, hMT2JetSig);

  float TauFakeData  = hEtaPtTauFakeData->Integral();
  float JetData      = hEtaPtJetData    ->Integral();
  float FakeData     = TauFakeData/JetData;
  cout<<" FakeData "<<FakeData<<endl;
  hEtaPtTauFakeData->Divide(hEtaPtJetData);
  
  Cout(p++, hMT2JetSigData);
  hMT2JetSigData->Scale(FakeData);
  
  Cout(p++, hMT2JetSigData);
  Cout(p++, hMT2TauSigData);
  hMT2JetSigData->Add(hMT2TauSigData, -1);
  hMT2JetSigData->Scale(-1.0);
  Cout(p++, hMT2JetSigData);
    

  hMT2JetSigData->Divide(hMT2TauFound);
  
  Cout(p++, hMT2JetSigData);
  
  int k = 1;

  int Chi2  = 0;

  float Sim = 0;

  float Pre = 0;

  float PreData    = 0;

  float SimErr     = 0;

  float PreErr     = 0;

  float PreErrData = 0;

  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);

  Cout(p++, hRatioRelaxed);
  Cout(p++, hRatioTight);

  hRatioRelaxed->Add(hRatioTight);
  
  hRatioTight->Divide(hRatioRelaxed);

  Cout(p++, hRatioTight);

  hMT2JetSig->Multiply(hRatioTight);

  hMT2JetSigData->Multiply(hRatioTight);

  hNAllTaus->Multiply(hRatioTight);
 
  for(int j = 1; j < NumberOfMT2Bins; j++){
    
    float TotalErrSim = hMT2JetSig    ->GetBinError(j);
    float StatErrSim  = (hMT2TauSig    ->GetBinError(j)/hMT2TauFound->GetBinContent(j))*hRatioTight->GetBinContent(j);
    float SysErrSim   = sqrt(fabs(TotalErrSim * TotalErrSim - StatErrSim * StatErrSim));

    float TotalErrData = hMT2JetSigData    ->GetBinError(j);
    float StatErrData  = (hMT2TauSigData    ->GetBinError(j)/hMT2TauFound->GetBinContent(j))*hRatioTight->GetBinContent(j);
    float SysErrData   = sqrt(fabs(TotalErrData * TotalErrData - StatErrData * StatErrData));
    
    cout<<"$"<<xMT2bin[j-1]<<"-"<<xMT2bin[j]<<"$: Sim  "<<hNAllTaus     ->GetBinContent(j)<<" +/- "<<hNAllTaus        ->GetBinError(j)*hRatioTight->GetBinContent(j)<<
                                                 " Pre  "<<hMT2JetSig    ->GetBinContent(j)<<" +/- "<<StatErrSim <<" +/- "<< SysErrSim <<                              
                                                 " Data "<<hMT2JetSigData->GetBinContent(j)<<" +/- "<<StatErrData<<" +/- "<< SysErrData<<endl;
  }

  for(int j = 0; j < NumberOfSamples; j++){
     AddOverAndUnderFlow(MT2[j]);
  }
  printYield();

THStack* h_stack     = new THStack(varname, "");
  for(int j = 0; j < NumberOfSamples; j++){
    MT2[j]->Rebin(3);
    TH1F* mt2 = (TH1F*)MT2[j]->Clone();
    mt2->SetName("mt2");
    if(j < (NumberOfSamples - 2))
      h_stack  -> Add(MT2[j]);
    delete mt2;
  }  

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W+jets", "f");
  Legend1->AddEntry(MT2[2], "Z+jets", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[5], "data", "l");
  
  //  TLegend *Legend1;
  printHisto(h_stack, MT2[5], MT2[4], Legend1 , "MTC", "hist", true, "MT2", "Events", -4, 0);

  plotRatioStack(h_stack,  MT2[4], MT2[5], true, false, "MT2_ratio", Legend1, "MT2", "Events", -4, 0);

  TCanvas *CanvasMT2 = new TCanvas("myMT2", "myMT2");
  hNAllTaus      ->SetLineColor(kBlue);
  hMT2JetSig     ->SetLineColor(kRed);
  hMT2JetSigData ->SetLineColor(kGreen);
  hNAllTaus      ->Draw();
  hMT2JetSig     ->Draw("sames");
  hMT2JetSigData ->Draw("sames");
}

void MassPlotter::setFlags(int flag){
    mT2  = 0;
    mT2b = 0;
    High = 0;
    Low  = 0;
    BTagRelaxed = 0;

    if(flag < 10)  mT2 = 1;
    if(flag > 10)  mT2b= 1;
    if(flag == 5 || flag == 15)  Low = 1;
    if(flag == 7 || flag == 17)  High = 1;

    if(flag == 19){
      High = 1;
      BTagRelaxed = 1;}


    if(mT2){
      MT2Min = 150;

      xMT2bin[1] = 200.0;
      xMT2bin[2] = 275.0;
      xMT2bin[3] = 375.0;
      xMT2bin[4] = 500.0;
      xMT2bin[5] = 2000.0;
   
      
      newxMT2bin[0] = 0.0;
      newxMT2bin[1] = 30.0;
      newxMT2bin[2] = 90.0;
      newxMT2bin[3] = 120.0;
      newxMT2bin[4] = 150.0;
      newxMT2bin[5] = 200.0;
      newxMT2bin[6] = 275.0;
      newxMT2bin[7] = 375.0;
      newxMT2bin[8] = 500.0;
      newxMT2bin[9] = 750.0;
    }


    if(mT2b){
      MT2Min = 125;

      xMT2bin[1] = 150.0;
      xMT2bin[2] = 200.0;
      xMT2bin[3] = 300.0;
      xMT2bin[4] = 2000.0;
      xMT2bin[5] = 3000.0;
  
      newxMT2bin[0] = 0.0;
      newxMT2bin[1] = 25.0;
      newxMT2bin[2] = 50.0;
      newxMT2bin[3] = 75.0;
      newxMT2bin[4] = 100.0;
      newxMT2bin[5] = 125.0;
      newxMT2bin[6] = 150.0;
      newxMT2bin[7] = 200.0;
      newxMT2bin[8] = 300.0;
      newxMT2bin[9] = 750.0;
    }

    if(mT2b && High){
      MT2Min = 125;

      xMT2bin[1] = 150.0;
      xMT2bin[2] = 180.0;
      xMT2bin[3] = 260.0;
      xMT2bin[4] = 2000.0;
      xMT2bin[5] = 3000.0;


      newxMT2bin[0] = 0.0;
      newxMT2bin[1] = 25.0;
      newxMT2bin[2] = 50.0;
      newxMT2bin[3] = 75.0;
      newxMT2bin[4] = 100.0;
      newxMT2bin[5] = 125.0;
      newxMT2bin[6] = 150.0;
      newxMT2bin[7] = 180.0;
      newxMT2bin[8] = 260.0;
      newxMT2bin[9] = 750.0;
 
    }
    xMT2bin[0] = MT2Min;
  }


void MassPlotter::AddOverAndUnderFlow(TH1 * Histo){
  // Add underflow & overflow bins
  // This failed for older ROOT version when the first(last) bin is empty
  // and there are underflow (overflow) events --- must check whether this 
  // is still the case
  Histo->SetBinContent(1,
		       Histo->GetBinContent(0) + Histo->GetBinContent(1));
  Histo->SetBinError(1,
		     sqrt(Histo->GetBinError(0)*Histo->GetBinError(0)+
			  Histo->GetBinError(1)*Histo->GetBinError(1) ));
  Histo->SetBinContent(Histo->GetNbinsX(),
		       Histo->GetBinContent(Histo->GetNbinsX()  )+ 
		       Histo->GetBinContent(Histo->GetNbinsX()+1) );
  Histo->SetBinError(Histo->GetNbinsX(),
		     sqrt(Histo->GetBinError(Histo->GetNbinsX()  )*
			  Histo->GetBinError(Histo->GetNbinsX()  )+
			  Histo->GetBinError(Histo->GetNbinsX()+1)*
			  Histo->GetBinError(Histo->GetNbinsX()+1)  ));
  Histo->SetBinContent(Histo->GetNbinsX() + 1, 0.0);
  Histo->SetBinContent(0, 0.0);
}

void MassPlotter::Cout(int k, TH1F * Histo){
  cout<<k<<" "<<Histo->GetName()<<endl;
  for(int j = 1; j <= Histo->GetNbinsX(); j++)
    cout<<"Bin :: "<<j<<" "<<Histo->GetBinContent(j)<<" +/- "<<Histo->GetBinError(j)<<endl;
}

void MassPlotter::Cout(int k, TH2F * Histo){
  cout<<k<<" "<<Histo->GetName()<<endl;
  int nYbins = Histo->GetNbinsY();
  int nXbins = Histo->GetNbinsX();

  for(int i = 1; i <= nYbins; i++){
    for(int j = 1; j <= nXbins; j++){
      cout<<Histo->GetBinContent((nXbins + 2) * i + j)<<" +/- "<<Histo->GetBinError((nXbins + 2) * i + j)<<endl;
    }
  }
}

void MassPlotter::printYield(){
  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);

    cout<<"           & "<<MT2[0]->GetName()     << " & "<<MT2[2]->GetName()      <<" & "<<MT2[1]->GetName()      <<" & "<<MT2[3]->GetName()     <<" & "<<MT2[4]->GetName()<<"& "<<MT2[5]->GetName()     <<"\\"<<endl;
  
    int tmpxMT2bin[10];
  
    for(int i = 0; i < NumberOfMT2Bins; i++)
      tmpxMT2bin[i] = xMT2bin[i]/5;

    double err;

    MT2[4]->IntegralAndError(0,tmpxMT2bin[NumberOfMT2Bins],err);
    cout<<" full selection & "<<MT2[0]->Integral()     << " & "<<MT2[2]->Integral()      <<" & "<<MT2[1]->Integral()      <<" & "<<MT2[3]->Integral()     <<" & "<<MT2[4]->Integral(0,tmpxMT2bin[NumberOfMT2Bins])<<" $pm$ "<<err<<" & "<<MT2[5]->Integral()     <<"\\"<<endl;


    for(int i = 0; i < NumberOfMT2Bins-1; i++){
      MT2[4]->IntegralAndError(tmpxMT2bin[i]+1,tmpxMT2bin[i+1],err);
      cout<<"$"<<xMT2bin[i]<<"-"<<xMT2bin[i+1]<<"$ &      "<<MT2[0]->Integral(tmpxMT2bin[i]+1,tmpxMT2bin[i+1])<< " & "<<MT2[2]->Integral(tmpxMT2bin[i]+1,tmpxMT2bin[i+1]) <<" & "<<MT2[1]->Integral(tmpxMT2bin[i]+1,tmpxMT2bin[i+1]) <<" & "<<MT2[3]->Integral(tmpxMT2bin[i]+1,tmpxMT2bin[i+1])<<" & "<<MT2[4]->Integral(tmpxMT2bin[i]+1,tmpxMT2bin[i+1])<<" $pm$ "<<err<<" & "<<MT2[5]->Integral(tmpxMT2bin[i]+1,tmpxMT2bin[i+1])<<"\\"<<endl;
    }
  }

