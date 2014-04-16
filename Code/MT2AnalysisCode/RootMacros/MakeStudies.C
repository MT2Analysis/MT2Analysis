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

//run via root -l -b -q MakeStudies.C++

using namespace std;

void MakeStudies();
void load(const char* filename = "/shome/haweber/CMSSW_4_2_3/src/DiLeptonAnalysis/NTupleProducer/MT2Analysis/Code/MT2AnalysisCode/MT2Code/samples/datasamples/samples_2141_dataonly.dat");
void Make1DPlotsRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, TH1D* histsusy, bool logFlag, bool normalize, TString outname, TString outputdir, TLegend *leg, TString xtitle, TString ytitle, float overlayScale);

//struct that combines MT2trees with necessary information
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
const int fVerbose = 3;
TString fPath;


Bool_t  met    = false;//use MET triggers (low HT region)
Bool_t  ht     = true; //use HT triggers (medium+high HT region), one of the two should be true, other false

Bool_t  fSave                     = true; //save the plots, needs doplotting = true
Bool_t  logflag                   = true; //make y-axis in plot logstyle
Bool_t  doplotting                = false;//make plots, default = false (because otherwise you create really tons of plots, use MakeStudiesPlotsFast.C)

Bool_t runData                    = true; //run over data - default = false, as these are object studies for sensitivity
Bool_t calcsusy                   = true; //run over susy sample - default = true
Bool_t runQCD                     = true; //run over QCD sample - default = true (unless study only high MT2)

//note that this code does not have any kind of reweighting - not needed as this is a sensitivity study for object/cut definitions
//this code makes plots/histograms for several configurations (i.e. definitions of certain high level objects)
//to find most sensitive one (be careful as sensitivity might be defined by a certain susy model point and not be valid for others)
//and stores them into a file (or plots them directly)
//uncomment below the variable you want to study (3x - in 'variables to study', in histogram definition, and in event loop)
void MakeStudies(){

	//samples to load
	TString samples;
	if(ht)       samples = "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/samples/samples_type1MET_CHS_53X_HT.dat";//HT
	else if(met) samples = "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/samples/samples_type1MET_CHS_53X_MET.dat";//MET
	TString outputdir;// directory of output
	if(ht)       outputdir = "MT2Studies/CSVM40/new/HT/";
	else if(met) outputdir = "MT2Studies/CSVM40/new/MET/";
	Util::MakeOutputDir(outputdir);

	TString outputname = "HistogramsPVMinDPhiNJets.root";//file where histograms are stored

	map<string, TH1D*> histograms;
        map<string, THStack*> stacks;
        TLegend *Legend1 = new TLegend(.71,.54,.91,.92);
	vector<string> histonames; histonames.clear();


	//define signal region
	const int sampletypesize = 9;
	string sample_type[sampletypesize] = {"QCD", "WJets", "ZJets", "TTbar", "SingleTop", "Other", "mc", "susy", "data"};//same as type in samples.dat

	const int signalregionsize = 9;
	string signalregions[signalregionsize] = {"2j0b", "2j1to2b", "3to5j0b", "3to5j1b", "3to5j2b", "6j0b", "6j1b", "6j2b", "3b"};
	int NMT2bins[signalregionsize]     = { 40,  30,  40,  30,  25,  30,  25,  25,  25};
	double MT2upedge[signalregionsize] = {800, 600, 800, 600, 500, 600, 500, 500, 500};

	const int controlsize = 2;
	string control[controlsize] = {"had", "1l"};

	//variables to study
	const int mindphisize = 5;//5;//3;
	//string minmetjetdphi[mindphisize] = {"0p3", "0p4", "0p5"};//MinPhi
	//string minmetjetdphi[mindphisize] = {"le50", "le70", "le90", "le120", "le150"};//VSPT
	//string minmetjetdphi[mindphisize] = {"NoPass", "Pass20", "Pass40", "Pass50"};//PassJetID
	string minmetjetdphi[mindphisize] = {"2", "3", "4", "5", "6"};//MinPhiNJets

	const int nverticessize = 4;
	string numvert[nverticessize] = {"0to12V","13to19V","20upV", "all"};

	//MT2 input jet threshold
	const int jetthresholdsize = 3;//2;//3;//4;
						// 1.      4.        2.     3.
	string jetthreshold[jetthresholdsize] = {"j20"/*, "j30"*/, "j40", "j50"};

	const int histokinds = 7;

	//here I do not use a map - to many name definitions
	//create histograms
	TH1D   *histos[sampletypesize][signalregionsize][controlsize][mindphisize][nverticessize][jetthresholdsize][histokinds];
	THStack *stack[signalregionsize][controlsize][mindphisize][nverticessize][jetthresholdsize][histokinds];
	cout << "histos["<<"sampletypesize"<<"]["<<"signalregionsize"<<"]["<<"controlsize"<<"]["<<"mindphisize"<<"]["<<"nverticessize"<<"]["<<"jetthresholdsize"<<"]["<<histokinds<<"];" << endl;
	int numhistos = 0;
	for(int is = 0; is<sampletypesize;   ++is){
	for(int js = 0; js<signalregionsize; ++js){
	for(int ks = 0; ks<controlsize;      ++ks){
	for(int ls = 0; ls<mindphisize;      ++ls){
	for(int ms = 0; ms<nverticessize;    ++ms){
	for(int ns = 0; ns<jetthresholdsize; ++ns){
		string hs = "_" + minmetjetdphi[ls] + "_" + signalregions[js] + "_" + control[ks] + "_" + numvert[ms] + "_" + sample_type[is];
		string mapname;
		int i = 0;
/*		//MinDPhi study
		mapname = "MT2" + jetthreshold[ns] + "_NoMinDPhi"; histonames.push_back(mapname); mapname = mapname + hs;//i==0
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhi"; histonames.push_back(mapname); mapname = mapname + hs;//i==1
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhi4"; histonames.push_back(mapname); mapname = mapname + hs;//i==2
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhiPt40"; histonames.push_back(mapname); mapname = mapname + hs;//i==3
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhi4Pt40"; histonames.push_back(mapname); mapname = mapname + hs;//i==4
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhiPt50"; histonames.push_back(mapname); mapname = mapname + hs;//i==5
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhi4Pt50"; histonames.push_back(mapname); mapname = mapname + hs;//i==6
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
*/
/*		//VSPT study
		mapname = "MT2" + jetthreshold[ns] + "_noVSPT";histonames.push_back(mapname); mapname = mapname + hs;//i==0
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
		mapname = "MT2" + jetthreshold[ns] + "_VSPTj20";histonames.push_back(mapname); mapname = mapname + hs;//i==1
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
		mapname = "MT2" + jetthreshold[ns] + "_VSPTj40";histonames.push_back(mapname); mapname = mapname + hs;//i==2
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
		mapname = "MT2" + jetthreshold[ns] + "_VSPTj50";histonames.push_back(mapname); mapname = mapname + hs;//i==3
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
*/
/*		//PassJetID study
		mapname = "MT2" + jetthreshold[ns] + "_JetIDPt";histonames.push_back(mapname); mapname = mapname + hs;//i==0
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
*/
		//MinDPhi40NJets study
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhi20NJets"; histonames.push_back(mapname); mapname = mapname + hs;//i==0
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
		mapname = "MT2" + jetthreshold[ns] + "_MinDPhi40NJets"; histonames.push_back(mapname); mapname = mapname + hs;//i==1
		if(is==7) stack[js][ks][ls][ms][ns][i] = new THStack( mapname.c_str(), mapname.c_str());
		histos[is][js][ks][ls][ms][ns][i++] = new TH1D( mapname.c_str(), mapname.c_str(), NMT2bins[js], 0., MT2upedge[js]);
	

		numhistos = i;
		for(int j = 0; j<numhistos; ++j) {
			histos[is][js][ks][ls][ms][ns][j]->Sumw2();
			histos[is][js][ks][ls][ms][ns][j]->SetStats(false);
			histos[is][js][ks][ls][ms][ns][j]->GetXaxis()->SetTitle("M_{T2} [GeV]");
			histos[is][js][ks][ls][ms][ns][j]->GetYaxis()->SetTitle("Events / 20 GeV");
		}
	}}}}}}
	
	std::ostringstream cutStream;
	cutStream       << " " 	  
	//define susy signal point you want to study
     	<< "(misc.ProcessID!=10||(Susy.MassChi==350&&Susy.MassLSP==50))" << "&&" // T2bb M_stop = 350, M_lsp = 50
//    	<< "(misc.ProcessID!=10||(Susy.MassChi==600&&Susy.MassLSP==50))" << "&&" // T2bb M_stop = 600, M_lsp = 50

	//define event selection
	  << "misc.MET>=30"                                      << "&&"
	  << "misc.Jet0Pass ==1"                                 << "&&"
	  << "misc.Jet1Pass ==1"                                 << "&&"
	  << "misc.PassJetID ==1"                                << "&&"
	  << "misc.Vectorsumpt < 70"                             << "&&"
//        JetSkim
	  << "NJetsIDLoose40 >=2"                                << "&&"
// Noise 
	  << "( misc.ProcessID==10 || misc.HBHENoiseFlag == 0)"  << "&&"
	  << "misc.CSCTightHaloIDFlag == 0"                      << "&&"
	  << "(misc.ProcessID==10||misc.hcalLaserEventFlag==0)"  << "&&"
	  << "misc.trackingFailureFlag == 0"                     << "&&"
	  << "misc.eeBadScFlag == 0"                             << "&&"
	  << "misc.EcalDeadCellTriggerPrimitiveFlag == 0";

	if(ht) cutStream << "&& misc.HT >=750";
	else if(met) cutStream << "&& misc.MET >=200";
	TString cuts = cutStream.str().c_str();

	std::ostringstream triggerStream;
	if(ht){
	triggerStream << "( "
			<< "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || "
			<< "trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
			<< "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1)";
	} else if(met){
	triggerStream << "( "
			<< "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
			<< "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )";
	}
	TString trigger = triggerStream.str().c_str();

	load(samples.Data());//load samples

	for(size_t i = 0; i < fSamples.size(); ++i){

   	    if(runData==false && fSamples[i].type=="data") continue;
	    if(calcsusy==false && fSamples[i].type=="susy") continue;
	    if(runQCD==false   && fSamples[i].sname=="QCD") continue;

		//define sampletype
	    string sampletype = (string)fSamples[i].type;
	    if(sampletype==(string)"mc"){
		if(fSamples[i].sname=="QCD")         sampletype = (string)"QCD";
		else if(fSamples[i].sname=="Wtolnu") sampletype = (string)"WJets";
		else if(fSamples[i].sname=="DY")     sampletype = (string)"ZJets";
		else if(fSamples[i].name=="TTbar")   sampletype = (string)"TTbar";
		else if(fSamples[i].sname=="TTbarV") sampletype = (string)"TTbar";
		else if(fSamples[i].sname=="Top")    sampletype = (string)"SingleTop";//no ttbar
		else sampletype = (string)"Other";
	    }

		//define global histogram bins within histogram-array
	    int ctrl;
	    int sr;
	    int st;

	    if(     sampletype=="QCD"      ) st = 0;
	    else if(sampletype=="WJets"    ) st = 1;
	    else if(sampletype=="ZJets"    ) st = 2;
	    else if(sampletype=="TTbar"    ) st = 3;
	    else if(sampletype=="SingleTop") st = 4;
	    else if(sampletype=="Other"    ) st = 5;
	    else if(sampletype=="susy"     ) st = 7;
	    else if(sampletype=="data"     ) st = 8;
	    else {st = 5; cout << "sample ? " << sampletype << endl;}

	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
	    if(fVerbose>2) cout << "MakeStudies: looping over " << fSamples[i].name << endl;
	    if(fVerbose>2) cout << "              sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;
	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;

	    TString myCuts = cuts;
	    if( fSamples[i].type=="data") { myCuts += " && " + trigger; }

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
            
            if ( fVerbose>2 && counter % 50000 == 0 )  cout << "+++ Proccessing event " << counter << endl;

            Double_t weight = sample_weight;
            if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 

		//define event-by-event histogram bins within histogram-array
	    int NumBJets = fMT2tree->NBJets40CSVM;
	    int NumJets  = fMT2tree->NJetsIDLoose40;
	    int NumTags = 0;
	    if(NumBJets == 0) NumTags = 0;
	    else if(NumBJets == 1 && fMT2tree->NJetsIDLoose40==2) NumTags = -1;
	    else if(NumBJets == 1) NumTags = 1;
	    else if(NumBJets == 2 && fMT2tree->NJetsIDLoose40==2) NumTags = -1;
	    else if(NumBJets == 2) NumTags = 2;
	    else if(NumBJets >= 3) NumTags = -3;
	    else cout << __LINE__ << " WTF NBJets " << NumBJets << "  NJets " << fMT2tree->NJetsIDLoose40 << endl;

	    if(     (fMT2tree->NEles + fMT2tree->NMuons)==0) ctrl = 0;//"had";	
	    else if((fMT2tree->NEles + fMT2tree->NMuons)==1) ctrl = 1;//"1l";
	    else continue;//skip all events with more than one lepton for now
	    if(     NumJets == 2 &&               NumBJets == 0) sr = 0;//"2j0b";
	    else if(NumJets == 2 &&               NumBJets >= 1) sr = 1;//"2j1to2b";
	    else if(NumJets >= 3 && NumJets<=5 && NumBJets == 0) sr = 2;//"3to5j0b";
	    else if(NumJets >= 3 && NumJets<=5 && NumBJets == 1) sr = 3;//"3to5j1b";
	    else if(NumJets >= 3 && NumJets<=5 && NumBJets == 2) sr = 4;//"3to5j2b";
	    else if(NumJets >= 6 &&               NumBJets == 0) sr = 5;//"6j0b";
	    else if(NumJets >= 6 &&               NumBJets == 1) sr = 6;//"6j1b";
	    else if(NumJets >= 6 &&               NumBJets == 2) sr = 7;//"6j2b";
	    else if(                              NumBJets >= 3) sr = 8;//"3b";
	    else cout << __LINE__ << " WTF NBJets " << NumBJets << "  NJets " << NumJets << endl;

	    int nv;
	    if(     fMT2tree->pileUp.NVertices<=12) nv = 0;
	    else if(fMT2tree->pileUp.NVertices<=20) nv = 1;
	    else                                    nv = 2;
	/*
	//MinDPhiStudy
	    histos[st][sr][ctrl][0][nv][0][0]->Fill(fMT2tree->misc.MT2     , weight);//_0p3 is dummy
	    histos[st][sr][ctrl][0][nv][1][0]->Fill(fMT2tree->misc.MT2jet40, weight);//_0p3 is dummy
	    histos[st][sr][ctrl][0][nv][2][0]->Fill(fMT2tree->misc.MT2jet50, weight);//_0p3 is dummy
	    if(fMT2tree->misc.MinMetJetDPhi>0.3){
		histos[st][sr][ctrl][0][nv][0][1]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][0][nv][1][1]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][0][nv][2][1]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhi>0.4){
		histos[st][sr][ctrl][1][nv][0][1]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][1][nv][1][1]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][1][nv][2][1]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhi>0.5){
		histos[st][sr][ctrl][2][nv][0][1]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][2][nv][1][1]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][2][nv][2][1]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhi4>0.3){
		histos[st][sr][ctrl][0][nv][0][2]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][0][nv][1][2]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][0][nv][2][2]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhi4>0.4){
		histos[st][sr][ctrl][1][nv][0][2]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][1][nv][1][2]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][1][nv][2][2]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhi4>0.5){
		histos[st][sr][ctrl][2][nv][0][2]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][2][nv][1][2]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][2][nv][2][2]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhiPt40>0.3){
		histos[st][sr][ctrl][0][nv][0][3]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][0][nv][1][3]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][0][nv][2][3]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhiPt40>0.4){
		histos[st][sr][ctrl][1][nv][0][3]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][1][nv][1][3]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][1][nv][2][3]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhiPt40>0.5){
		histos[st][sr][ctrl][2][nv][0][3]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][2][nv][1][3]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][2][nv][2][3]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhi4Pt40>0.3){
		histos[st][sr][ctrl][0][nv][0][4]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][0][nv][1][4]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][0][nv][2][4]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhi4Pt40>0.4){
		histos[st][sr][ctrl][1][nv][0][4]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][1][nv][1][4]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][1][nv][2][4]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhi4Pt40>0.5){
		histos[st][sr][ctrl][2][nv][0][4]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][2][nv][1][4]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][2][nv][2][4]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhiPt50>0.3){
		histos[st][sr][ctrl][0][nv][0][5]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][0][nv][1][5]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][0][nv][2][5]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhiPt50>0.4){
		histos[st][sr][ctrl][1][nv][0][5]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][1][nv][1][5]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][1][nv][2][5]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhiPt50>0.5){
		histos[st][sr][ctrl][2][nv][0][5]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][2][nv][1][5]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][2][nv][2][5]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhi4Pt50>0.3){
		histos[st][sr][ctrl][0][nv][0][6]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][0][nv][1][6]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][0][nv][2][6]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhi4Pt50>0.4){
		histos[st][sr][ctrl][1][nv][0][6]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][1][nv][1][6]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][1][nv][2][6]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } if(fMT2tree->misc.MinMetJetDPhi4Pt50>0.5){
		histos[st][sr][ctrl][2][nv][0][6]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][2][nv][1][6]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][2][nv][2][6]->Fill(fMT2tree->misc.MT2jet50, weight);
	    } 
	*/
/*	//VSPT study
	    double VSPTj20 = fMT2tree->misc.Vectorsumpt;
	    double VSPTj40 = fMT2tree->GetMHTminusMET(1, 40, 2.4, true);
	    double VSPTj50 = fMT2tree->GetMHTminusMET(1, 50, 2.4, true);
	    histos[st][sr][ctrl][0][nv][0][0]->Fill(fMT2tree->misc.MT2     , weight);//_le50 is dummy
	    histos[st][sr][ctrl][0][nv][1][0]->Fill(fMT2tree->misc.MT2jet40, weight);//_le50 is dummy
	      if(VSPTj20<50){
	    histos[st][sr][ctrl][0][nv][0][1]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][0][nv][1][1]->Fill(fMT2tree->misc.MT2jet40, weight);
	    } if(VSPTj20<70){
	    histos[st][sr][ctrl][1][nv][0][1]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][1][nv][1][1]->Fill(fMT2tree->misc.MT2jet40, weight);
            } if(VSPTj20<90){
	    histos[st][sr][ctrl][2][nv][0][1]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][2][nv][1][1]->Fill(fMT2tree->misc.MT2jet40, weight);
            } if(VSPTj20<120){
	    histos[st][sr][ctrl][3][nv][0][1]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][3][nv][1][1]->Fill(fMT2tree->misc.MT2jet40, weight);
            } if(VSPTj20<150){
	    histos[st][sr][ctrl][4][nv][0][1]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][4][nv][1][1]->Fill(fMT2tree->misc.MT2jet40, weight);
            } if(VSPTj40<50){
	    histos[st][sr][ctrl][0][nv][0][2]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][0][nv][1][2]->Fill(fMT2tree->misc.MT2jet40, weight);
	    } if(VSPTj40<70){
	    histos[st][sr][ctrl][1][nv][0][2]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][1][nv][1][2]->Fill(fMT2tree->misc.MT2jet40, weight);
            } if(VSPTj40<90){
	    histos[st][sr][ctrl][2][nv][0][2]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][2][nv][1][2]->Fill(fMT2tree->misc.MT2jet40, weight);
            } if(VSPTj40<120){
	    histos[st][sr][ctrl][3][nv][0][2]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][3][nv][1][2]->Fill(fMT2tree->misc.MT2jet40, weight);
            } if(VSPTj40<150){
	    histos[st][sr][ctrl][4][nv][0][2]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][4][nv][1][2]->Fill(fMT2tree->misc.MT2jet40, weight);
            } if(VSPTj50<50){
	    histos[st][sr][ctrl][0][nv][0][3]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][0][nv][1][3]->Fill(fMT2tree->misc.MT2jet40, weight);
	    } if(VSPTj50<70){
	    histos[st][sr][ctrl][1][nv][0][3]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][1][nv][1][3]->Fill(fMT2tree->misc.MT2jet40, weight);
            } if(VSPTj50<90){
	    histos[st][sr][ctrl][2][nv][0][3]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][2][nv][1][3]->Fill(fMT2tree->misc.MT2jet40, weight);
            } if(VSPTj50<120){
	    histos[st][sr][ctrl][3][nv][0][3]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][3][nv][1][3]->Fill(fMT2tree->misc.MT2jet40, weight);
            } if(VSPTj50<150){
	    histos[st][sr][ctrl][4][nv][0][3]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][4][nv][1][3]->Fill(fMT2tree->misc.MT2jet40, weight);
            }
	*/
	/*	//PassJetId
	    histos[st][sr][ctrl][0][nv][0][0]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][0][nv][1][0]->Fill(fMT2tree->misc.MT2jet40, weight);
	    if(fMT2tree->PassJetID(20,2.4,1)){
	    histos[st][sr][ctrl][1][nv][0][0]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][1][nv][1][0]->Fill(fMT2tree->misc.MT2jet40, weight);
	    } if(fMT2tree->misc.PassJet40ID){
	    histos[st][sr][ctrl][2][nv][0][0]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][2][nv][1][0]->Fill(fMT2tree->misc.MT2jet40, weight);
	    } if(fMT2tree->misc.PassJetID){
	    histos[st][sr][ctrl][3][nv][0][0]->Fill(fMT2tree->misc.MT2     , weight);
	    histos[st][sr][ctrl][3][nv][1][0]->Fill(fMT2tree->misc.MT2jet40, weight);
	    }
	*/
		//MinDPhiNJetsstudy
		double mindphi20no4 = fMT2tree->misc.MinMetJetDPhi4;
		double mindphi40no4 = fMT2tree->misc.MinMetJetDPhi4Pt40;
		double mindphi20no2 = fMT2tree->MinMetJetDPhi(0,20,5.0,1,2);
		double mindphi20no3 = fMT2tree->MinMetJetDPhi(0,20,5.0,1,3);
		double mindphi20no5 = fMT2tree->MinMetJetDPhi(0,20,5.0,1,5);
		double mindphi20no6 = fMT2tree->MinMetJetDPhi(0,20,5.0,1,6);
		double mindphi40no2 = fMT2tree->MinMetJetDPhi(0,40,5.0,1,2);
		double mindphi40no3 = fMT2tree->MinMetJetDPhi(0,40,5.0,1,3);
		double mindphi40no5 = fMT2tree->MinMetJetDPhi(0,40,5.0,1,5);
		double mindphi40no6 = fMT2tree->MinMetJetDPhi(0,40,5.0,1,6);
	     if(mindphi20no2>0.3){
		histos[st][sr][ctrl][0][nv][0][0]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][0][nv][1][0]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][0][nv][2][0]->Fill(fMT2tree->misc.MT2jet50, weight);
	   } if(mindphi20no3>0.3){
		histos[st][sr][ctrl][1][nv][0][0]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][1][nv][1][0]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][1][nv][2][0]->Fill(fMT2tree->misc.MT2jet50, weight);
	   } if(mindphi20no4>0.3){
		histos[st][sr][ctrl][2][nv][0][0]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][2][nv][1][0]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][2][nv][2][0]->Fill(fMT2tree->misc.MT2jet50, weight);
	   } if(mindphi20no5>0.3){
		histos[st][sr][ctrl][3][nv][0][0]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][3][nv][1][0]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][3][nv][2][0]->Fill(fMT2tree->misc.MT2jet50, weight);
	   } if(mindphi20no6>0.3){
		histos[st][sr][ctrl][4][nv][0][0]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][4][nv][1][0]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][4][nv][2][0]->Fill(fMT2tree->misc.MT2jet50, weight);
	   } if(mindphi40no2>0.3){
		histos[st][sr][ctrl][0][nv][0][1]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][0][nv][1][1]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][0][nv][2][1]->Fill(fMT2tree->misc.MT2jet50, weight);
	   } if(mindphi40no3>0.3){
		histos[st][sr][ctrl][1][nv][0][1]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][1][nv][1][1]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][1][nv][2][1]->Fill(fMT2tree->misc.MT2jet50, weight);
	   } if(mindphi40no4>0.3){
		histos[st][sr][ctrl][2][nv][0][1]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][2][nv][1][1]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][2][nv][2][1]->Fill(fMT2tree->misc.MT2jet50, weight);
	   } if(mindphi40no5>0.3){
		histos[st][sr][ctrl][3][nv][0][1]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][3][nv][1][1]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][3][nv][2][1]->Fill(fMT2tree->misc.MT2jet50, weight);
	   } if(mindphi40no6>0.3){
		histos[st][sr][ctrl][4][nv][0][1]->Fill(fMT2tree->misc.MT2     , weight);
		histos[st][sr][ctrl][4][nv][1][1]->Fill(fMT2tree->misc.MT2jet40, weight);
		histos[st][sr][ctrl][4][nv][2][1]->Fill(fMT2tree->misc.MT2jet50, weight);
	   }

	}//while

	delete fMT2tree;

	}//for sample

	cout << "adding overflow/underflow to 1st, last bin" << endl;
	for(int is = 0; is<sampletypesize;   ++is){
	for(int js = 0; js<signalregionsize; ++js){
	for(int ks = 0; ks<controlsize;      ++ks){
	for(int ls = 0; ls<mindphisize;      ++ls){
	for(int ms = 0; ms<nverticessize;    ++ms){
	for(int ns = 0; ns<jetthresholdsize; ++ns){
	for(int zs = 0; zs<numhistos;        ++zs){
		histos[is][js][ks][ls][ms][ns][zs]->SetBinContent(1,
				    histos[is][js][ks][ls][ms][ns][zs]->GetBinContent(0) + histos[is][js][ks][ls][ms][ns][zs]->GetBinContent(1));
		histos[is][js][ks][ls][ms][ns][zs]->SetBinError(1,
				  sqrt(histos[is][js][ks][ls][ms][ns][zs]->GetBinError(0)*histos[is][js][ks][ls][ms][ns][zs]->GetBinError(0)+
				       histos[is][js][ks][ls][ms][ns][zs]->GetBinError(1)*histos[is][js][ks][ls][ms][ns][zs]->GetBinError(1) ));
		histos[is][js][ks][ls][ms][ns][zs]->SetBinContent(histos[is][js][ks][ls][ms][ns][zs]->GetNbinsX(),
					    histos[is][js][ks][ls][ms][ns][zs]->GetBinContent(histos[is][js][ks][ls][ms][ns][zs]->GetNbinsX()  )+ 
					    histos[is][js][ks][ls][ms][ns][zs]->GetBinContent(histos[is][js][ks][ls][ms][ns][zs]->GetNbinsX()+1) );
		histos[is][js][ks][ls][ms][ns][zs]->SetBinError(histos[is][js][ks][ls][ms][ns][zs]->GetNbinsX(),
					  sqrt(histos[is][js][ks][ls][ms][ns][zs]->GetBinError(histos[is][js][ks][ls][ms][ns][zs]->GetNbinsX()  )*
					       histos[is][js][ks][ls][ms][ns][zs]->GetBinError(histos[is][js][ks][ls][ms][ns][zs]->GetNbinsX()  )+
					       histos[is][js][ks][ls][ms][ns][zs]->GetBinError(histos[is][js][ks][ls][ms][ns][zs]->GetNbinsX()+1)*
					       histos[is][js][ks][ls][ms][ns][zs]->GetBinError(histos[is][js][ks][ls][ms][ns][zs]->GetNbinsX()+1)  ));
	}}}}}}}


	cout << "Adding all mc samples to mc histo" << endl;
	for(int is = 0; is<sampletypesize;   ++is){
	for(int js = 0; js<signalregionsize; ++js){
	for(int ks = 0; ks<controlsize;      ++ks){
	for(int ls = 0; ls<mindphisize;      ++ls){
	for(int ms = 0; ms<nverticessize;    ++ms){
	for(int ns = 0; ns<jetthresholdsize; ++ns){
	for(int zs = 0; zs<numhistos;        ++zs){
		if(is==8 || is==7 || is==6) continue;
		histos[6][js][ks][ls][ms][ns][zs]->Add(histos[is][js][ks][ls][ms][ns][zs],1);
	}}}}}}}

	cout << "Adding pv samples " << endl;
	for(int is = 0; is<sampletypesize;   ++is){
	for(int js = 0; js<signalregionsize; ++js){
	for(int ks = 0; ks<controlsize;      ++ks){
	for(int ls = 0; ls<mindphisize;      ++ls){
	for(int ms = 0; ms<nverticessize;    ++ms){
	for(int ns = 0; ns<jetthresholdsize; ++ns){
	for(int zs = 0; zs<numhistos;        ++zs){
		if(ms==3) continue;
		histos[is][js][ks][ls][3][ns][zs]->Add(histos[is][js][ks][ls][ms][ns][zs],1);
	}}}}}}}

	cout << "setting stack color" << endl;
	Legend1 -> SetFillColor(0);
   	Legend1 -> SetBorderSize(0);
	for(int is = 0; is<sampletypesize;   ++is){
	for(int js = 0; js<signalregionsize; ++js){
	for(int ks = 0; ks<controlsize;      ++ks){
	for(int ls = 0; ls<mindphisize;      ++ls){
	for(int ms = 0; ms<nverticessize;    ++ms){
	for(int ns = 0; ns<jetthresholdsize; ++ns){
	for(int zs = 0; zs<numhistos;        ++zs){
		Color_t colour = 603;//set default others if we have some failure
		     if(is==0) colour = 401;
		else if(is==1) colour = 417;
		else if(is==2) colour = 419;
		else if(is==3) colour = 855;
		else if(is==4) colour = 600;
		else if(is==5) colour = 603;
		else if(is==6) colour = 603;
		if(is<=6){
			histos[is][js][ks][ls][ms][ns][zs]->SetFillColor(colour);
			histos[is][js][ks][ls][ms][ns][zs]->SetLineColor(colour);
		}
		else if(is==7){
			histos[is][js][ks][ls][ms][ns][zs]->SetLineColor(kBlack);
			histos[is][js][ks][ls][ms][ns][zs]->SetLineStyle(kDotted);
			histos[is][js][ks][ls][ms][ns][zs]->SetLineWidth(4);
		}
		else if(is==8){
			histos[is][js][ks][ls][ms][ns][zs]->SetLineColor(kBlack);
			histos[is][js][ks][ls][ms][ns][zs]->SetMarkerStyle(20);
			histos[is][js][ks][ls][ms][ns][zs]->SetLineWidth(1);
		}
	}}}}}}}

	cout << "setting stacks and legend" << endl;
	bool leggy = true;
	for(int is = 0; is<sampletypesize;   ++is){
	for(int js = 0; js<signalregionsize; ++js){
	for(int ks = 0; ks<controlsize;      ++ks){
	for(int ls = 0; ls<mindphisize;      ++ls){
	for(int ms = 0; ms<nverticessize;    ++ms){
	for(int ns = 0; ns<jetthresholdsize; ++ns){
	for(int zs = 0; zs<numhistos;        ++zs){
		if(is==8 || is==7 || is==6) continue;
		stack[js][ks][ls][ms][ns][zs]->Add(histos[is][js][ks][ls][ms][ns][zs]);
	}}}}}}}
	for(int js = 0; js<signalregionsize; ++js){
	for(int ks = 0; ks<controlsize;      ++ks){
	for(int ls = 0; ls<mindphisize;      ++ls){
	for(int ms = 0; ms<nverticessize;    ++ms){
	for(int ns = 0; ns<jetthresholdsize; ++ns){
	for(int zs = 0; zs<numhistos;        ++zs){
	if(leggy){
	for(int is = 0; is<sampletypesize;   ++is){
		if(is==6) continue;
		if(is!=6) Legend1->AddEntry(histos[is][js][ks][ls][ms][ns][zs], sample_type[is].c_str(), "f");
		else Legend1->AddEntry(histos[is][js][ks][ls][ms][ns][zs], (sample_type[is]).c_str(), "l");
		leggy = false;
	}}}}}}}}

	cout << "saving " << endl;
    TFile *fsavefile = new TFile(outputdir + outputname,"RECREATE");//do not delete this one
	fsavefile->cd();
	for(int is = 0; is<sampletypesize;   ++is){
	for(int js = 0; js<signalregionsize; ++js){
	for(int ks = 0; ks<controlsize;      ++ks){
	for(int ls = 0; ls<mindphisize;      ++ls){
	for(int ms = 0; ms<nverticessize;    ++ms){
	for(int ns = 0; ns<jetthresholdsize; ++ns){
	for(int zs = 0; zs<numhistos;        ++zs){
		histos[is][js][ks][ls][ms][ns][zs]->Write();
		if(is==0) stack[js][ks][ls][ms][ns][zs]->Write();
	}}}}}}}

	fsavefile->Close();
	cout << "Saved histograms in " << outputdir << outputname << endl;


	if(doplotting){
		cout << "plotting histograms ..." << endl;

		for(int js = 0; js<signalregionsize; ++js){
		for(int ks = 0; ks<controlsize;      ++ks){
		for(int ls = 0; ls<mindphisize;      ++ls){
		for(int ms = 0; ms<nverticessize;    ++ms){
		for(int ns = 0; ns<jetthresholdsize; ++ns){
		for(int zs = 0; zs<numhistos;        ++zs){
			string name = histos[8][js][ks][ls][ms][ns][zs]->GetName();
			TString ytitle = "Events / 20 GeV";//histos[name+hs1]->GetYaxis()->GetTitle();
			TString xtitle = "M_{T2} [GeV]";
			TString outname = name + (logflag ? "_log" : "") + "_overlay";
			TH1D *hmc   = (TH1D*)histos[6][js][ks][ls][ms][ns][zs]->Clone("hmc");
			TH1D *hsusy = (TH1D*)histos[7][js][ks][ls][ms][ns][zs]->Clone("hsusy");
			TH1D *hdata = (TH1D*)histos[8][js][ks][ls][ms][ns][zs]->Clone("hdata");
			THStack *st = (THStack*)stack[js][ks][ls][ms][ns][zs]->Clone("stack");
			if(histos[6][js][ks][ls][ms][ns][zs]->Integral()>0 || histos[8][js][ks][ls][ms][ns][zs]->Integral()>0){
				Make1DPlotsRatio(st, hmc, hdata, hsusy, logflag, false, outname, outputdir, Legend1, xtitle, ytitle, 1.);
			}
	
		}}}}}}

	}

	delete Legend1;
	delete fsavefile;
}//function

//copy from a MassPlotter.cc function - no detail comments
void Make1DPlotsRatio(THStack *hstack, TH1D* histmc, TH1D *histdata, TH1D* histsusy, bool logFlag, bool normalize, TString outname, TString outputdir, TLegend *leg, TString xtitle, TString ytitle, float overlayScale){


	// define canvas and pads 
	TH1D *h1 = (TH1D*)histmc->Clone("h1_copy");
	TH1D *h2 = (TH1D*)histdata->Clone("h2_copy");
	TH1D *h3 = (TH1D*)histsusy->Clone("h3_copy");

	h1->SetTitle("");
	h2->SetTitle("");
	h3->SetTitle("");
	hstack->SetTitle("");

	h1->SetStats(0);	
	h2->SetStats(0);	
	h1->SetMarkerStyle(20);
	h2->SetMarkerStyle(20);

	TCanvas* c1 = new TCanvas(outname+"_ratio","",0,0,600,600);
	c1->SetFrameLineWidth(1);
	c1 -> cd();

 	TPad *p_plot  = new TPad(outname+"_plotpad",  "Pad containing the overlay plot", 0,0.211838,1,1);
	p_plot->SetLeftMargin(0.131579);
	p_plot->SetRightMargin(0.08);
	p_plot->SetTopMargin(0.06895515);
	p_plot->SetBottomMargin(0.07206074);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad(outname+"_ratiopad", "Pad containing the ratio",   0,0.01863354,0.9967105,0.2189441);
	p_ratio->SetLeftMargin(0.1336634);	
	p_ratio->SetRightMargin(0.075);
	p_ratio->SetTopMargin(0.06976745);
	p_ratio->SetBottomMargin(0.2790698);

	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	if(logFlag) gPad->SetLogy(1);
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
	hstack->SetMinimum(0.02);
	hstack->Draw("hist");
	hstack->GetYaxis()->SetTitle(ytitle);
	hstack->GetYaxis()->SetLabelSize(0.05);
	hstack->GetYaxis()->SetTitleSize(0.05);
	hstack->GetYaxis()->SetTitleOffset(1.3);

	h2    ->Draw("sameE");
	h3->Scale(overlayScale ? overlayScale : h2->Integral() / h3->Integral());
	h3->SetFillColor(0);
	h3->SetLineStyle(kDotted);
	h3->Draw("samehist");

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.03);
	TString text ="";
	text = outname;
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
 
	TH1D *h_ratio = (TH1D*)histdata->Clone("h2_copy_2");	
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);

 	h_ratio ->Divide(h2, h1);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	//MC with errors
	TH1D*h_ratio_mc = (TH1D*)histmc->Clone("h1_copy_2");
	h_ratio_mc->Divide(h1);
	h_ratio_mc->GetYaxis()->SetRangeUser(0,2);
	h_ratio_mc->GetXaxis()->SetLabelSize( 0.);
	h_ratio_mc->GetYaxis()->SetTitle("Data / MC");
	h_ratio_mc->GetXaxis()->SetTitle(xtitle);
	h_ratio_mc->GetXaxis()->SetTitleSize(0.2);
	h_ratio_mc->GetXaxis()->SetTitleOffset(0.5);
	h_ratio_mc->GetYaxis()->SetLabelSize(0.19);
	h_ratio_mc->GetXaxis()->SetTickLength(0.09);
	h_ratio_mc->GetYaxis()->SetTitleSize(0.18);
	h_ratio_mc->GetYaxis()->SetTitleOffset(0.36);
	h_ratio_mc->GetYaxis()->SetNdivisions(509);//xxxnew

	h_ratio   ->SetTitle("");
	h_ratio_mc->SetTitle("");
	
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

	TString save=outname+"_ratio";
	if(fSave)Util::Print(c1, save, outputdir);
}

//read samples.dat
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
