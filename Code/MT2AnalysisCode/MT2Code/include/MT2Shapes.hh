/*****************************************************************************
*   Small Class to make Shape-Templated for MT2Analysis                      *
*****************************************************************************/

#ifndef MT2Shapes_HH
#define MT2Shapes_HH

#include "MT2tree.hh"
#include "helper/Utilities.hh"
#include "helper/Monitor.hh"
#include "THStack.h"
#include "TTree.h"
#include <map>

//________________________________________________________________________________
class MT2Shapes  {

public:
	MT2Shapes();
	MT2Shapes(TString);
	MT2Shapes(TString, TString);
	virtual ~MT2Shapes();


	void init(TString filename = "samples.dat");
	void loadSamples(const char* filename = "samples.dat");

	struct sample{
		TString name;
		TString sname;
		TString type;
		TString shapename;
		TFile *file;
		TTree *tree;
		float xsection;
		float nevents;
		float kfact;
		float lumi;
		int color;
	};
	std::vector<sample>  fSamples;
	

	void setVerbose(int v){ fVerbose = v;};
	void setOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
	void setOutputFile(TString filename){ fOutputFile = Util::MakeOutputFile(fOutputDir + filename); };
	void GetShapes(TString var, TString cuts, int njets, int nleps, TString selection_name, TString HLT,
			         TString xtitle, const int nbins, const double *bins);
	void GetShapes(TString var, TString cuts, int njets, int nleps, TString selection_name, TString HLT,
			         TString xtitle, const int nbins, const double min, const double max);

private:

	TString fOutputDir;
	TFile *fOutputFile;
	int fVerbose;
	TString fPath;

	MT2tree* fMT2tree;
	TTree*   fTree;

	void DrawHisto(TH1* h_orig, TString canvname,  Option_t *drawopt);
	void FixOverAndUnderflowBins(TH1D* h);

};

#endif
