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
#include "TRandom.h"
#include "TRandom3.h"
#include "TMap.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <cmath>
#include <limits.h>
#include <utility>
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use your own path of MT2tree.hh
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"

//use via root -l -b -q PrintOutEvents.C++

using namespace std;

void load(const char* filename = "samples/datasamples/samples_2141_dataonly.dat");
void PrintOutEvents();

//struct combining MT2trees with important information like cross section
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


//this code create the so-called 'PrintOuts' as done with the MT2tree::PrintOut() function
void PrintOutEvents(){

	bool savetofile = true;//save the printout to the file PrintOut.log, default = true
	bool onlyData = true;  //run only over data, default = true

//	TString  samples = "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/samples/samples_MET_filter.dat";//only dummy
//	TString  samples = "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/samples/samples_HT_filter.dat";//only dummy
	TString samples = "samples/samples_HTandMET_filter.dat";

	//event selection for which you want to have the printout - complicated stuff can also be done later
	std::ostringstream cutStream;
	std::ostringstream cutStreamBase;
	cutStream << " " 
	//	<< "NBJets40CSVM >= 3"		<< "&&"
		<< "NJetsIDLoose40 >= 2"		<< "&&"
		<< "NEles+NMuons+NTausIDLoose3Hits  ==0"<< "&&"
	//	<< "misc.MT2>=250"                      << "&&"
		<< "misc.Jet0Pass ==1"                      << "&&"
		<< "misc.Jet1Pass ==1"                      << "&&"
		<< "misc.Vectorsumpt < 70";
		cutStream  << "&& misc.MinMetJetDPhi4Pt40 >0.3";
//	cutStream << "&&misc.HT<750&&misc.HT>=450&&misc.MET>200";
//	cutStream << "&&misc.HT>=750&&misc.MET>30";
	cutStream << "&&((misc.HT<750&&misc.HT>=450&&misc.MET>200)||(misc.HT>=750&&misc.MET>30))";

	//Higgs
//	cutStream <<"&&NJetsIDLoose40 >= 4"<<"&&"
  //  		  <<"NBJetsCSVM>=2"<<"&&"
	//	  <<"GetSelBBMinv()>250";
	
	cutStreamBase << " " 
      << "misc.PassJet40ID ==1"                      << "&&"
      << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"  << "&&" // for rare SM samples
      << "misc.CSCTightHaloIDFlag == 0"             << "&&"
      << "misc.trackingFailureFlag==0"              << "&&"
      << "misc.eeBadScFlag==0"                      << "&&"
      << "misc.EcalDeadCellTriggerPrimitiveFlag==0" << "&&"
      << "misc.TrackingManyStripClusFlag==0"             << "&&"
      << "misc.TrackingTooManyStripClusFlag==0"          << "&&"
      << "misc.TrackingLogErrorTooManyClustersFlag==0"   << "&&"
      << "misc.CrazyHCAL==0";
	cutStreamBase << "&&misc.MET/misc.CaloMETRaw<=2.";
	cutStreamBase << "&&TOBTECTagger<=8&&ExtraBeamHaloFilter==0";// ONLY ON DATA
	TString cuts = cutStream.str().c_str();
	TString basecuts = cutStreamBase.str().c_str();

	load(samples.Data());

	//if have extremely many events can only select a random subset of those
	TRandom3 *random = new TRandom3(0);//somehow truly random with this seed
	vector<int> randvec; randvec.clear();
	for(int i =0; i<40;++i)
		randvec.push_back(random->Integer(7750));
	delete random;
	sort(randvec.begin(), randvec.end());
	//for(int i=0; i<randvec.size();++i) cout << "randvec["<<i<<"]="<<randvec[i] << endl;
	int counter1 = 0; int counter2 = 0;

	//don't run events twice - this vector makes that sure (in case event is in both MET and HT trigger stream)
	   vector<pair<int,pair<int,int> > > rls; rls.clear();


   	for(size_t i = 0; i < fSamples.size(); ++i){
        
   	    if(onlyData && fSamples[i].type!="data") continue;

	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts + "&&" + basecuts;
        
	    //if( fSamples[i].type=="data") myCuts += " && " + trigger; //cuts to be aplied only on data
        
	    cout << "sample " << fSamples[i].name << endl;
   	    cout << "Cuts for Flow: " << myCuts << endl;
   	    fSamples[i].tree->Draw(">>selList", myCuts);
        
	    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
	    fSamples[i].tree->SetEventList(myEvtList);
	    int counter=0;
	    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
	    if(myEvtList->GetSize()==0) continue;
		//run over selected events
        while(myEvtList->GetEntry(counter++) !=-1){	
		int jentry = myEvtList->GetEntry(counter-1);
            
		nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
		fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
		//implementation of running over random subset
//		if(counter1==randvec[counter2]){
//			++counter2;
//			fMT2tree->PrintOut(savetofile);
//		}

		pair <int,int> dummy1(fMT2tree->misc.LumiSection,fMT2tree->misc.Event);
		pair <int,pair<int,int> > dummy(fMT2tree->misc.Run, dummy1);
		bool newevt = true;
		for(int nn = 0; nn<rls.size(); ++nn){
			if(rls[nn].first==dummy.first && rls[nn].second.first==dummy.second.first && rls[nn].second.second==dummy.second.second) {newevt = false; break;}
		}
		if(newevt){
			rls.push_back(dummy);
			fMT2tree->PrintOut(savetofile);//safe the printout
			++counter1;
		}

	}
	delete fMT2tree;
	delete fSamples[i].tree;
	}
	cout << "Printed out " << counter1 << " events." << endl;
}

//read the samples.dat
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
			counter++;
			fSamples.push_back(s);
		}
	}
	if(fVerbose > 0) cout << "------------------------------------" << endl;
}
