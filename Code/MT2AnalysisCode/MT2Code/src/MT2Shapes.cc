/*****************************************************************************
*   Small Class to make Shape-Templated for MT2Analysis                      *
*****************************************************************************/

#include "MT2Shapes.hh"

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
MT2Shapes::MT2Shapes(){
// Default constructor, no samples are set
}

//____________________________________________________________________________
MT2Shapes::MT2Shapes(TString outputdir){
// Explicit constructor with output directory
	setOutputDir(outputdir);
}

//____________________________________________________________________________
MT2Shapes::MT2Shapes(TString outputdir, TString outputfile){
// Explicit constructor with output directory and output file
	setOutputDir(outputdir);
	setOutputFile(outputfile);
}

//____________________________________________________________________________
MT2Shapes::~MT2Shapes(){
	fOutputFile->Close();
	delete fOutputFile;
}

//____________________________________________________________________________
void MT2Shapes::init(TString filename){
	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Initializing MT2Shapes ... " << endl;
	loadSamples(filename);
}

//____________________________________________________________________________
void MT2Shapes::GetShapes( TString var, TString cuts, int njets, int nleps, TString selection_name, TString HLT,
		          TString xtitle, const int nbins, const double min, const double max){
	
	double bins[nbins];
	bins[0] = min;
	for(int i=1; i<=nbins; i++) bins[i] = min+i*(max-min)/nbins;
	GetShapes(var, cuts, njets, nleps, selection_name, HLT, xtitle, nbins, bins);
}

//________________________________________________________________________

void MT2Shapes::GetShapes(TString var, TString cuts, int njets, int nleps, TString selection_name, TString HLT,
			  TString xtitle, const int nbins, const double *bins){
	if(fVerbose > 1) {
		cout << "....................................................................................................." << endl;
		cout << "........ GetShapes for " << selection_name << "..............................................................." << endl;
	}

	// parsiong cuts
	TString nJets = "NJetsIDLoose";
	nJets += njets < 0 ? ">=" : "==";
	nJets += TString::Format("%d",abs(njets));
	TString  nLeps;
	if     (nleps < 0 )  nLeps = " && (NEles + NMuons) >=";
	else if(nleps >=0  ) nLeps = " && (NEles + NMuons) ==";
	nLeps += TString::Format("%d",abs(nleps));
	if     (nleps ==-10) nLeps = " "; 
	if     (nleps ==-11) nLeps = " && NEles ==1 && NMuons ==0"; 
	if     (nleps ==-13) nLeps = " && NEles ==0 && NMuons ==1";

	TString basecuts = nJets + nLeps + "&&" + cuts;
	if(fVerbose > 1){
		cout << "\n---> applied cuts:     " << basecuts           << endl;
		cout << "\n---> trigger for data: " << HLT        << "\n" << endl;
	}


	// vector of all histos
	vector<TH1D*>   h_shapes;
	vector<TString> ShapeNames;
	int nShapes=0;

	for(size_t i = 0; i < fSamples.size(); ++i){
		// see if sample should be added to existing shape or not
		bool newshape(true);
		for (int shape=0; shape<ShapeNames.size(); ++shape){
			if(ShapeNames[shape]==fSamples[i].shapename) newshape=false;
		}
		if(newshape) {
			ShapeNames.push_back(fSamples[i].shapename);
			h_shapes.  push_back(new TH1D(fSamples[i].shapename+"_"+selection_name, fSamples[i].shapename+"_"+selection_name, nbins, bins));
			nShapes++;
			h_shapes[nShapes-1]->SetLineColor(fSamples[i].color);
			h_shapes[nShapes-1]-> Sumw2();
		}

		// curr histo
		TH1D* hcurr= new TH1D("hcurr", "", nbins, bins);
		hcurr      ->SetStats(0);
		hcurr      ->Sumw2();
		
		// calculate weight and print some info
		Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
		if(fVerbose>2) cout << "---> looping over " << fSamples[i].sname << "--------------------------------------------"<< endl;
		if(fVerbose>2) cout << "          sample belongs to shape " << fSamples[i].shapename << endl;
		if(fVerbose>2) cout << "          sample has weight " << weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	
		TString variable  = TString::Format("%s>>%s",var.Data(),hcurr->GetName());
		TString theCuts = basecuts;
		if(fSamples[i].type=="data" && HLT!="") theCuts += " &&("+HLT+")"; // triggers for data

		TString selection;
		if(fSamples[i].type!="data") selection      = TString::Format("(%.15f*pileUp.Weight) * (%s)",weight,theCuts.Data());
		else                         selection      = TString::Format("(%.15f) * (%s)"              ,weight,theCuts.Data()); 
		  
		int nev = fSamples[i].tree->Draw(variable.Data(),selection.Data(),"goff");

		// fix over and underflow bins
		FixOverAndUnderflowBins(hcurr);

		/// n MC events passing the cuts and weigthed event count with errors
		TH1F * clone = (TH1F*)hcurr->Clone();
		clone->Rebin(clone->GetNbinsX());
		if(fVerbose>2) cout << "          MC events found : "  <<  nev << endl;
		if(fVerbose>2) cout << "          Events: " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;
		delete clone;

		// draw clone of hist
		for (int iShape=0; iShape<ShapeNames.size(); ++iShape){
			if(ShapeNames[iShape]==fSamples[i].shapename) {
				h_shapes[iShape]->Add(hcurr);
				break;
			}
		}
		delete hcurr;
	}

	// plot clone of hostos and save histos in fOutputFile
	fOutputFile->cd();
	for (int iShape=0; iShape<h_shapes.size(); ++iShape){
		h_shapes[iShape]->SetXTitle(xtitle);
		if(ShapeNames[iShape]=="Data") {h_shapes[iShape]->SetMarkerStyle(20);h_shapes[iShape]->SetMarkerColor(kBlack);}
		DrawHisto(h_shapes[iShape], h_shapes[iShape]->GetName(),ShapeNames[iShape]=="Data"?"EX0":"hist");
		h_shapes[iShape]->Write();
	}

	// cleanup
	for(unsigned int iShape=0; iShape<h_shapes.size(); ++iShape){
		delete h_shapes[iShape];
	}

}

//____________________________________________________________________________
void MT2Shapes::loadSamples(const char* filename){
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
			sscanf(buffer, "Nevents\t%f", &ParValue);
			s.nevents = ParValue;
			
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

			if(fVerbose > 0){
				cout << " ---- " << endl;
				cout << "  New sample added: " << s.name << endl;
				cout << "   Sample no.      " << counter << endl;
				cout << "   Short name:     " << s.sname << endl;
				cout << "   ShapeName:      " << s.shapename << endl;
				cout << "   File:           " << (s.file)->GetName() << endl;
				cout << "   Events:         " << s.nevents  << endl;
				cout << "   Events in tree: " << s.tree->GetEntries() << endl; 
				cout << "   Xsection:       " << s.xsection << endl;
				cout << "   Lumi:           " << s.lumi << endl;
				cout << "   kfactor:        " << s.kfact << endl;
				cout << "   type:           " << s.type << endl;
				cout << "   Color:          " << s.color << endl;
			}
			fSamples.push_back(s);
			counter++;
		}
	}
	if(fVerbose > 0) cout << "------------------------------------" << endl;
}

//______________________________________________________________________________
void MT2Shapes::FixOverAndUnderflowBins(TH1D*h ){
		// Add underflow & overflow bins
		// This failed for older ROOT version when the first(last) bin is empty
		// and there are underflow (overflow) events --- must check whether this 
		// is still the case
		h->SetBinContent(1,                h->GetBinContent(0) + h->GetBinContent(1));
		h->SetBinError  (1,                sqrt(h->GetBinError(0)*h->GetBinError(0)+h->GetBinError(1)*h->GetBinError(1) ));
		h->SetBinContent(h->GetNbinsX(),h->GetBinContent(h->GetNbinsX()  )+ h->GetBinContent(h->GetNbinsX()+1) );
		h->SetBinError  (h->GetNbinsX(),sqrt(h->GetBinError(h->GetNbinsX()  )*h->GetBinError(h->GetNbinsX()  )+h->GetBinError(h->GetNbinsX()+1)*h->GetBinError(h->GetNbinsX()+1)  ));
}

//_______________________________________________________________________________
void MT2Shapes::DrawHisto(TH1* h_orig, TString name,  Option_t *drawopt){
	TH1D* h = (TH1D*)h_orig->Clone(name);
	TCanvas *col = new TCanvas(name, "", 0, 0, 600, 500);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	gStyle->SetPalette(1);
	gPad->SetRightMargin(0.15);
	h->DrawCopy(drawopt);
	gPad->RedrawAxis();
}
