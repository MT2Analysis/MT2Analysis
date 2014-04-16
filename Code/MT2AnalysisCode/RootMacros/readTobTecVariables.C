#include "TEventList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCut.h"
#include "TROOT.h"
#include "TSystem.h"
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

//run via root -l -b -q readTobTecVariables.C++

using namespace std;

//from the AOD files + TobTec Tagger variables that Boris produced
//a text file/table has been produced by using tree->Scan() function
//this short macro uses this table/text file, reads in all variables
//and stores them to a more usuable root file
void readTobTecVariables(){

TFile *newfile = new TFile("tobtec.root","RECREATE");
TTree *newtree = new TTree("tobtec","tobtec");

   Bool_t  RA2bFilterRECO;
   Bool_t  RA2bFilterRERECO;
   Float_t BorisFilterRECO;
   Float_t BorisFilterRERECO;
   UInt_t   Event;
   Int_t   LumiSection;
   Int_t   Run;

   newtree->Branch("BorisFilterRECO", &BorisFilterRECO, "BorisFilterRECO/F");
   newtree->Branch("BorisFilterRERECO", &BorisFilterRERECO, "BorisFilterRERECO/F");
   newtree->Branch("RA2bFilterRECO", &RA2bFilterRECO, "RA2bFilterRECO/O");
   newtree->Branch("RA2bFilterRERECO", &RA2bFilterRERECO, "RA2bFilterRERECO/O");
   newtree->Branch("Run", &Run, "Run/I");
   newtree->Branch("LumiSection", &LumiSection, "LumiSection/I");
   newtree->Branch("Event", &Event, "Event/i");

   ifstream filterdat("TobTecVariables.dat");
   char buffer[200];
   //read in file
   while( filterdat.getline(buffer, 200, '\n') ){
	int rrun(-1), lls(-1), RA2bold(-1),RA2bnew(-1);
	unsigned int eevent(-1);
	bool oldRA2b(false), newRA2b(false);
	float Borisnew(-1), Borisold(-1);
	int d1(-1);
	sscanf(buffer, "*\t%d\t*\t%u\t*\t%d\t*\t%d\t*\t%f\t*\t%f\t*\t%d\t*\t%d\t*", &d1, &eevent, &lls, &rrun, &Borisnew, &Borisold,&RA2bnew,&RA2bold);
	if(eevent<0||eevent>INT_MAX) cout << rrun<<":"<<lls<<":"<<eevent << " - be careful for this event number" << endl;
	oldRA2b = bool(RA2bold); newRA2b = bool(RA2bnew);
	RA2bFilterRECO = oldRA2b;
	RA2bFilterRERECO = newRA2b;
	BorisFilterRECO = Borisold;
	BorisFilterRERECO = Borisnew;
	Event = eevent;
	LumiSection = lls;
	Run = rrun;
	newtree->Fill();//save variables to the tree
   }
   newtree->Print();//save the tree
   newtree->AutoSave();

}