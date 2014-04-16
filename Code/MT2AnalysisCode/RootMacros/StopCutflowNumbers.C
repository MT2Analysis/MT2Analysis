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

using namespace std;

//run via root -l -b -q StopCutflowNumbers.C++

void StopCutflowNumbers();
void GetContent3DRange(TH3D* num, double &eff, double &err, double xlow, double xup, double ylow, double yup, double zlow, double zup);
void GetContent3DRangeBin(TH3D* num, double &eff, double &err, int xlow, int xup, int ylow, int yup, int zlow, int zup);
void GetContent2DRange(TH2D* num, double &eff, double &err, double xlow, double xup, double ylow, double yup);
void GetContent2DRangeBin(TH2D* num, double &eff, double &err, int xlow, int xup, int ylow, int yup);
void GetContent1DRange(TH1D* num, double &eff, double &err, double xlow, double xup);
void GetContent1DRangeBin(TH1D* num, double &eff, double &err, int xlow, int xup);
void GetCutflowFootHead(bool head, bool twoD, bool mc, bool ttbar, bool stop, bool data);
void GetCutflowLine(string histoname, string title, double xlow, double xup, double ylow, double yup, bool twoD, bool witherror, bool mc, bool ttbar, bool stop, bool data, map<string, TH1D*> histos, map<string, TH2D*> histos2);
void GetCutflowLine3D(string histoname, string title, double xlow, double xup, double ylow, double yup, double zlow, double zup, bool witherror, bool mc, bool ttbar, bool stop, bool data, map<string, TH3D*> histos);
void GetCutflowFootHead3D(bool head, bool mc, bool ttbar, bool stop, bool data);

const int fVerbose = 3; // already defined below
TString outputdir                 = "Stops/test53X_corr_Stop400bChi300_SlepSnu200/bjet40_jpt40";
TString outputname                = "DiLeptonicStops.root";//default
Bool_t MC                         = true;
Bool_t TTbar                      = false;
Bool_t Stop                       = true;
Bool_t Data                       = true;
Bool_t WithErr                    = true;
Bool_t printonlyLL                = true;
Bool_t printSF                    = false;//only works if printonlyLL==false
Bool_t printEMu                   = false;//only works if printonlyLL==false and printSF==false

//cut used for dileptonic stop feasibility study.
//makes cut flow tables for different cuts ont MT2(l), MT2(b), MT2(lb)
//too get most sensitive cut (for chosen test point)
void StopCutflowNumbers(){

    map<string, TH1D*> histos;
    map<string, TH2D*> histos2;
    const unsigned int sampletypesize = 10;//11;
    string sample_type[sampletypesize] = {/*"QCD", */"WJets", "ZJets", "TTbar", "SingleTop", "TTbarV", "VV/VVV", "Other", "mc", "Stop", "data"};

    const unsigned int leptontypesize = 5;//4;
    string lepton_type[leptontypesize] = {"MuMu", "EMu", "EE", "SF", "LL"};

    vector<string> histonames; histonames.clear();
    vector<string> vs1; vs1.clear();
    vector<string> vs2; vs2.clear();
    string mapname;

		//load histograms from StopFeasibility.C
		mapname = "Mll_0b"; histonames.push_back(mapname);//0-700
		mapname = "MT2ll_eq1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_eq1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_eq1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_withMlbcut_eq1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_withMlbcut_eq1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2ll_1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_withMlbcut_1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_withMlbcut_1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2ll_2b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_2b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_2b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_withMlbcut_2b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_withMlbcut_2b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_withMlbcut_afterMT2llge85_1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_afterMT2llge85_1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85_1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_withMlbcut_afterMT2llge85_eq1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_afterMT2llge85_eq1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85_eq1b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_withMlbcut_afterMT2llge85_2b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_afterMT2llge85_2b"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85_2b"; histonames.push_back(mapname);//0-700
		mapname = "D2"; histonames.push_back(mapname);//0-300
		mapname = "MT2ll"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_withMlbcut"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_withMlbcut"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_withMlbcut_afterMT2llge85"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_withMlbcut_afterMT2llge85"; histonames.push_back(mapname);//0-700
		mapname = "MT2lb_massless_afterMT2llge85"; histonames.push_back(mapname);//0-700
		vs1 = vs;

		mapname = "MT2ll_vs_D2"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2ll_vs_MT2bb"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2ll_vs_MT2bb_lInMET"; vs2.push_back(mapname); histonames.push_back(mapname);
		mapname = "MT2ll_vs_MT2bb_lInMET_testmass80"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_massless_vs_D2"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_massless_vs_MET"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_massless_vs_MT2bb"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_massless_vs_MT2bb_lInMET"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_massless_vs_MT2bb_lInMET_testmass80"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_massless_vs_MT2ll"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_vs_D2"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_vs_MET"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_vs_MT2bb"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_vs_MT2bb_lInMET"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_vs_MT2bb_lInMET_testmass80"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_vs_MT2ll"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_withMlbcut_vs_MT2ll"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_withMlbcut_vs_MT2bb"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_withMlbcut_vs_MT2bb_lInMET"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_withMlbcut_vs_MT2bb_lInMET_mW"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_massless_withMlb_vs_MT2ll"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_massless_vs_MT2ll_1b"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_massless_withMlb_vs_MT2ll_1b"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_withMlb_vs_MT2ll_1b"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_massless_vs_MT2ll_2b"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_massless_withMlb_vs_MT2ll_2b"; vs2.push_back(mapname); histonames.push_back(mapname); 
		mapname = "MT2lb_withMlb_vs_MT2ll_2b"; vs2.push_back(mapname); histonames.push_back(mapname); 

    TFile *oldfile = TFile::Open(outputdir + outputname);//comment above

	for(unsigned int is = 0; is<sampletypesize; ++is){
           for(unsigned int il = 0; il<leptontypesize; ++il){
		string hs = string("_") + lepton_type[il] + string("_") + sample_type[is];
		for(unsigned int in = 0; in<vs1.size(); ++in){
			string name = vs1[in] + hs;
			if(histos.count(name)==0) {histos[name ] = (TH1D*)oldfile->Get((name).c_str() ); }//1D histograms
		}
		for(unsigned int in = 0; in<vs2.size(); ++in){
			string name = vs2[in] + hs;
			if(histos2.count(name)==0) {histos2[name ] = (TH2D*)oldfile->Get((name).c_str() ); }//2D histograms
		}
	   }
	}

	cout << endl << endl << outputdir << outputname << endl << endl;

	//make the cutflow tables
	cout << endl << ">=2b >=2b >=2b" << endl << endl;
	GetCutflowFootHead(true,false,MC,TTbar,Stop,Data);
	GetCutflowLine("MT2lb_massless_2b",                     "massless $M_{T2}^{min}(lb)$",                                          0.,700.,-999.,-999.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2lb_massless_2b",                     "massless $M_{T2}^{min}(lb)$",                                        200.,700.,-999.,-999.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2lb_withMlbcut_2b",                     "$M_{T2}^{min}(lb)$ with $M_{lb}$ cut",                                          0.,700.,-999.,-999.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2lb_withMlbcut_2b",                     "$M_{T2}^{min}(lb)$ with $M_{lb}$ cut",                                        230.,700.,-999.,-999.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2lb_withMlbcut_2b",                     "$M_{T2}^{min}(lb)$ with $M_{lb}$ cut",                                        240.,700.,-999.,-999.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2lb_withMlbcut_afterMT2llge85_2b",       "$M_{T2}^{min}(lb)$ with $M_{lb}$ cut histonames. $M_{T2}(l)$",           0.,700.,85.,700.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2lb_withMlbcut_afterMT2llge85_2b",       "$M_{T2}^{min}(lb)$ with $M_{lb}$ cut histonames. $M_{T2}(l)$",         220.,700.,85.,700.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2ll_2b",       "$M_{T2}(l)$",           0.,700.,85.,700.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2ll_2b",       "$M_{T2}(l)$",         130.,700.,85.,700.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowFootHead(false,false,MC,TTbar,Stop,Data);

	cout << endl << "=1b =1b =1b" << endl << endl;
	GetCutflowFootHead(true,false,MC,TTbar,Stop,Data);
	GetCutflowLine("MT2lb_massless_eq1b",                     "massless $M_{T2}^{min}(lb)$",                                          0.,700.,-999.,-999.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2lb_massless_eq1b",                     "massless $M_{T2}^{min}(lb)$",                                        200.,700.,-999.,-999.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2lb_withMlbcut_eq1b",                     "$M_{T2}^{min}(lb)$ with $M_{lb}$ cut",                                          0.,700.,-999.,-999.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2lb_withMlbcut_eq1b",                     "$M_{T2}^{min}(lb)$ with $M_{lb}$ cut",                                        230.,700.,-999.,-999.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2lb_withMlbcut_eq1b",                     "$M_{T2}^{min}(lb)$ with $M_{lb}$ cut",                                        240.,700.,-999.,-999.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2lb_withMlbcut_afterMT2llge85_eq1b",       "$M_{T2}^{min}(lb)$ with $M_{lb}$ cut histonames. $M_{T2}(l)$",           0.,700.,85.,700.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2lb_withMlbcut_afterMT2llge85_eq1b",       "$M_{T2}^{min}(lb)$ with $M_{lb}$ cut histonames. $M_{T2}(l)$",         220.,700.,85.,700.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2ll_eq1b",       "$M_{T2}(l)$",           0.,700.,85.,700.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowLine("MT2ll_eq1b",       "$M_{T2}(l)$",         130.,700.,85.,700.,false,WithErr,MC,TTbar,Stop,Data,histos,histos2);
	GetCutflowFootHead(false,false,MC,TTbar,Stop,Data);

}//void StopCutflowNumbers

//get a single line of the cut flow table
void GetCutflowLine(string histoname, string title, double xlow, double xup, double ylow, double yup, bool twoD, bool witherror, bool mc, bool ttbar, bool stop, bool data, map<string, TH1D*> histos, map<string, TH2D*> histos2){

	string lept = "_LL";//default
	if(printonlyLL)    lept = "_LL";
	else if (printSF ) lept = "_SF";
	else if (printEMu) lept = "_EMu";
	string hmc  = lept+"_mc"; string httbar = lept+"_TTbar"; string hstop = lept+"_Stop"; string hdata = lept+"_data";
	//string hmc  = lept+"_mc"; string httbar = lept+"_TTbar"; string hstop = lept+"_ZJets"; string hdata = lept+"_data";


	double mcsum, mcerr, ttsum, tterr, stopsum, stoperr, datasum, dataerr;
	if(twoD){
		GetContent2DRange(histos2[histoname+hmc],    mcsum,  mcerr,  xlow,xup,ylow,yup);
		GetContent2DRange(histos2[histoname+httbar], ttsum,  tterr,  xlow,xup,ylow,yup);
		GetContent2DRange(histos2[histoname+hstop],  stopsum,stoperr,xlow,xup,ylow,yup);
		GetContent2DRange(histos2[histoname+hdata],  datasum,dataerr,xlow,xup,ylow,yup);
	} else{
		GetContent1DRange(histos[histoname+hmc],     mcsum,  mcerr,  xlow,xup);
		GetContent1DRange(histos[histoname+httbar],  ttsum,  tterr,  xlow,xup);
		GetContent1DRange(histos[histoname+hstop],   stopsum,stoperr,xlow,xup);
		GetContent1DRange(histos[histoname+hdata],   datasum,dataerr,xlow,xup);
	}
	cout << " " << title << "    & " << fixed << setprecision(1) << (xlow) << " - " << (xup);
	if(twoD)     cout << fixed << setprecision(2) << "  & " << int(ylow) << " - " << int(yup);
	if(stop)     cout << fixed << setprecision(2) << "  & " << stopsum;
	if(stop && witherror) cout << fixed << setprecision(2) << " $\\pm$ " << stoperr;
	if(ttbar)    cout << fixed << setprecision(2) << "  & " << ttsum ;//<< "$\\pm$"<<tterr;
	if(ttbar && witherror) cout << fixed << setprecision(2) << " $\\pm$ " << tterr;
	if(mc)       cout << fixed << setprecision(2) << "  & " << mcsum;
	if(mc && witherror) cout << fixed << setprecision(2) << " $\\pm$ " << mcerr;
	if(data)     cout << fixed << setprecision(2) << "  & " << int(datasum);
	if(stop&&mc) { if(mcsum>0) cout << fixed << setprecision(3) << "  & " << stopsum/mcsum;
                       else        cout << fixed << setprecision(2) << "  & " << "-"; }
	if(stop&&mc) { if(mcsum>0) cout << fixed << setprecision(3) << "  & " << stopsum/sqrt(mcsum);
                       else        cout << fixed << setprecision(2) << "  & " << "-"; }
	cout << " \\\\" << endl;

}

//get top code or bottom code of a latex table
void GetCutflowFootHead(bool head, bool twoD, bool mc, bool ttbar, bool stop, bool data){
	if(head){
		int ncolumns = 1;//first column subtracted, its a l column;
		if(twoD) ++ncolumns;
		if(mc) ++ncolumns;
		if(ttbar) ++ncolumns;
		if(stop) ++ncolumns;
		if(data) ++ncolumns;
		cout << "\%BEGINLATEX\%"             << endl;
		cout << "\\begin{table}"             << endl
		<< "\\begin{center}"            << endl;
		cout << "\\tiny"                    << endl;
		cout << "\\begin{tabular}{l";
		for(int n = 0; n< ncolumns; ++n) cout << "c";
		if(mc&&stop) cout << "cc";
		             cout << "}"  << endl;
		             cout << "\\hline\\hline"             << endl;
		             cout << "  cut         & x-range (GeV) ";
		if(twoD)     cout << "  & y-range (GeV) ";
		if(stop)     cout << "  & $\\tilde{t}\\bar{\\tilde{t}}$  ";
		if(ttbar)    cout << " & $t\\bar{t}$ ";
		if(mc)       cout << " & Total MC ";
		if(data)     cout << " & data  ";
		if(stop&&mc) cout << " & S/B  ";
		if(stop&&mc) cout << " & $S/\\sqrt{B}$  ";
		             cout << "\\\\" << endl;
		             cout << "\\hline\\hline"             << endl;
	}
	else{
		cout << "\\hline\\hline"                                                                                                << endl
		<< "\\end{tabular}"                                                                                                << endl
		<< "\\end{center}"                                                                                                 << endl
		<< "\\end{table}"                                                                                                  << endl
		<< "\%ENDLATEX\%"                                                                                                  << endl
		<< endl;
	}
}

//get a line of a cutflow table from a 3D histogram
void GetCutflowLine3D(string histoname, string title, double xlow, double xup, double ylow, double yup, double zlow, double zup, bool witherror, bool mc, bool ttbar, bool stop, bool data, map<string, TH3D*> histos){

	string lept = "_LL";
	string hmc  = lept+"_mc"; string httbar = lept+"_TTbar"; string hstop = lept+"_Stop"; string hdata = lept+"_data";

	double mcsum, mcerr, ttsum, tterr, stopsum, stoperr, datasum, dataerr;
	GetContent3DRange(histos[histoname+hmc],    mcsum,  mcerr,  xlow,xup,ylow,yup,zlow,zup);
	GetContent3DRange(histos[histoname+httbar], ttsum,  tterr,  xlow,xup,ylow,yup,zlow,zup);
	GetContent3DRange(histos[histoname+hstop],  stopsum,stoperr,xlow,xup,ylow,yup,zlow,zup);
	GetContent3DRange(histos[histoname+hdata],  datasum,dataerr,xlow,xup,ylow,yup,zlow,zup);
	cout << " " << title << "    & " << fixed << setprecision(2) << int(xlow) << "-" << int(xup);
	cout << fixed << setprecision(2) << "  & " << int(ylow) << "-" << int(yup) << "  & " << int(zlow) << "-" << int(zup);
	if(stop)     cout << fixed << setprecision(2) << "  & " << stopsum;
	if(ttbar)    cout << fixed << setprecision(2) << "  & " << ttsum;
	if(mc)       cout << fixed << setprecision(2) << "  & " << mcsum;
	if(data)     cout << fixed << setprecision(2) << "  & " << int(datasum);
	if(stop&&mc) { if(mcsum>0) cout << fixed << setprecision(3) << "  & " << stopsum/mcsum;
                       else        cout << fixed << setprecision(2) << "  & " << "-"; }
	cout << " \\\\" << endl;

}

//get top code or bottom code of a latex table for 3D histogram cuts
void GetCutflowFootHead3D(bool head, bool mc, bool ttbar, bool stop, bool data){
	if(head){
		int ncolumns = 3;//first column subtracted, its a l column;
		if(mc) ++ncolumns;
		if(ttbar) ++ncolumns;
		if(stop) ++ncolumns;
		if(data) ++ncolumns;
		cout << "\%BEGINLATEX\%"             << endl;
		cout << "\\begin{table}"             << endl
		<< "\\begin{center}"            << endl;
		//<< "\\small"                    << endl;
		cout << "\\begin{tabular}{l";
		for(int n = 0; n< ncolumns; ++n) cout << "c";
		if(mc&&stop) cout << "c";
		             cout << "}"  << endl;
		             cout << "\\hline\\hline"             << endl;
		             cout << "  cut         & x-range (GeV) ";
		cout << "  & y-range (GeV)   & z-range (GeV)";
		if(stop)     cout << "  & $\\tilde{t}\\bar{\\tilde{t}}$  ";
		if(ttbar)    cout << " & $t\\bar{t}$ ";
		if(mc)       cout << " & Total MC ";
		if(data)     cout << " & data  ";
		if(stop&&mc) cout << " & S/B  ";
		             cout << "\\\\" << endl;
		             cout << "\\hline\\hline"             << endl;
	}
	else{
		cout << "\\hline\\hline"                                                                                                << endl
		<< "\\end{tabular}"                                                                                                << endl
		<< "\\end{center}"                                                                                                 << endl
		<< "\\end{table}"                                                                                                  << endl
		<< "\%ENDLATEX\%"                                                                                                  << endl
		<< endl;
	}
}

//get chosen bins of a 3D histogram by axis borders
void GetContent3DRange(TH3D* num, double &eff, double &err, double xlow, double xup, double ylow, double yup, double zlow, double zup){
    int xbinlow, xbinup, ybinlow, ybinup, zbinlow, zbinup;
    //change by 10e-6 to avoid border effects where one bin is included just as it is the lower border of a bin, so it in fact should not be included
    xbinlow = num->GetXaxis()->FindBin(xlow+10e-6);
    xbinup  = num->GetXaxis()->FindBin(xup -10e-6);
    ybinlow = num->GetYaxis()->FindBin(ylow+10e-6);
    ybinup  = num->GetYaxis()->FindBin(yup -10e-6);
    zbinlow = num->GetZaxis()->FindBin(zlow+10e-6);
    zbinup  = num->GetZaxis()->FindBin(zup -10e-6);
    GetContent3DRangeBin(num, eff, err, xbinlow, xbinup, ybinlow, ybinup, zbinlow, zbinup);
}

//get sum of the contents (and quadratic error) of chosen bins of a 3D histogram
void GetContent3DRangeBin(TH3D* num, double &eff, double &err, int xlow, int xup, int ylow, int yup, int zlow, int zup){
    TH3D *num_c;
    num_c = (TH3D*) num->Clone();
    vector<double> xx, yy, zz;
    for(int i = xlow; i<=xup; ++i){
        xx.push_back(num_c->GetXaxis()->GetBinLowEdge(i));}
    xx.push_back(num_c->GetXaxis()->GetBinLowEdge(xup)+num_c->GetXaxis()->GetBinWidth(xup));//upper border of last bin
    for(int i = ylow; i<=yup; ++i){
        yy.push_back(num_c->GetYaxis()->GetBinLowEdge(i));}
    yy.push_back(num_c->GetYaxis()->GetBinLowEdge(yup)+num_c->GetYaxis()->GetBinWidth(yup));//upper border of last bin
    for(int i = zlow; i<=zup; ++i){
        zz.push_back(num_c->GetZaxis()->GetBinLowEdge(i));}
    zz.push_back(num_c->GetZaxis()->GetBinLowEdge(zup)+num_c->GetZaxis()->GetBinWidth(zup));//upper border of last bin
    int sizex = xx.size();
    int sizey = yy.size();
    int sizez = zz.size();
    double x[sizex], y[sizey], z[sizez];
    for(int i = 0; i<sizex; ++i){ x[i] = xx[i]; }
    for(int i = 0; i<sizey; ++i){ y[i] = yy[i]; }
    for(int i = 0; i<sizez; ++i){ z[i] = zz[i]; }
    TH3D *num_sm = new TH3D("num_sm", "num_sm", sizex-1, x, sizey-1, y, sizez-1, z); num_sm->Sumw2();
    for(int nx=1; nx<=num_sm->GetNbinsX(); ++nx){
        for(int ny=1; ny<=num_sm->GetNbinsY(); ++ny){
            for(int nz=1; nz<=num_sm->GetNbinsZ(); ++nz){
                num_sm->SetBinContent(nx, ny, nz, num_c->GetBinContent(num_c->FindBin(num_sm->GetXaxis()->GetBinCenter(nx), num_sm->GetYaxis()->GetBinCenter(ny), num_sm->GetZaxis()->GetBinCenter(nz) )));
                num_sm->SetBinError(nx, ny, nz, num_c->GetBinError(num_c->FindBin(num_sm->GetXaxis()->GetBinCenter(nx), num_sm->GetYaxis()->GetBinCenter(ny), num_sm->GetZaxis()->GetBinCenter(nz) )));
            }
        }
    }
    //NOTE: why is this crashing
    //num_sm->Rebin3D(num_sm->GetNbinsX(),num_sm->GetNbinsY(),num_sm->GetNbinsZ() );//only in ROOT for >= 5.30
    (eff) = (double)num_sm->GetBinContent(1,1,1);
    (err) = (double)num_sm->GetBinError(1,1,1);
    //do the dirty way --> errors might be wrong due to weights - only temporary
    eff=0.;err=0.;
    for(int nx=1; nx<=num_sm->GetNbinsX(); ++nx){
        for(int ny=1; ny<=num_sm->GetNbinsY(); ++ny){
            for(int nz=1; nz<=num_sm->GetNbinsZ(); ++nz){
                eff += num_sm->GetBinContent(nx, ny, nz);
                err += pow(num_sm->GetBinError(nx, ny, nz),2);
            }
        }
    }
    err = sqrt(err);
    delete num_sm;
}

//get chosen bins of a 2D histogram by axis borders
void GetContent2DRange(TH2D* num, double &eff, double &err, double xlow, double xup, double ylow, double yup){
    int xbinlow, xbinup, ybinlow, ybinup;
    //change by 10e-6 to avoid border effects where one bin is included just as it is the lower border of a bin, so it in fact should not be included
    xbinlow = num->GetXaxis()->FindBin(xlow+10e-6);
    xbinup  = num->GetXaxis()->FindBin(xup -10e-6);
    ybinlow = num->GetYaxis()->FindBin(ylow+10e-6);
    ybinup  = num->GetYaxis()->FindBin(yup -10e-6);
    GetContent2DRangeBin(num, eff, err, xbinlow, xbinup, ybinlow, ybinup);
}
//get sum of the contents (and quadratic error) of chosen bins of a 2D histogram
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

//get chosen bins of a 1D histogram by axis borders
void GetContent1DRange(TH1D* num, double &eff, double &err, double xlow, double xup){
    int xbinlow, xbinup;
    //change by 10e-6 to avoid border effects where one bin is included just as it is the lower border of a bin, so it in fact should not be included
    xbinlow = num->GetXaxis()->FindBin(xlow+10e-6);
    xbinup  = num->GetXaxis()->FindBin(xup -10e-6);
    GetContent1DRangeBin(num, eff, err, xbinlow, xbinup);
}
//get sum of the contents (and quadratic error) of chosen bins of a 1D histogram
void GetContent1DRangeBin(TH1D* num, double &eff, double &err, int xlow, int xup){
    TH1D *num_c;
    num_c = (TH1D*) num->Clone();
    vector<double> xx;
    for(int i = xlow; i<=xup; ++i){
        xx.push_back(num_c->GetXaxis()->GetBinLowEdge(i));}
    xx.push_back(num_c->GetXaxis()->GetBinLowEdge(xup)+num_c->GetXaxis()->GetBinWidth(xup));//upper border of last bin
    int sizex = xx.size();
    double x[sizex];
    for(int i = 0; i<sizex; ++i){ x[i] = xx[i]; }
    TH1D *num_sm = new TH1D("num_sm", "num_sm", sizex-1, x); num_sm->Sumw2();
    for(int nx=1; nx<=num_sm->GetNbinsX(); ++nx){
            num_sm->SetBinContent(nx, num_c->GetBinContent(num_c->FindBin(num_sm->GetXaxis()->GetBinCenter(nx) )));
            num_sm->SetBinError(nx, num_c->GetBinError(num_c->FindBin(num_sm->GetXaxis()->GetBinCenter(nx) )));
    }
    num_sm->Rebin( num_sm->GetNbinsX() );
    (eff) = (double)num_sm->GetBinContent(1);
    (err) = (double)num_sm->GetBinError(1);
    delete num_sm;
}