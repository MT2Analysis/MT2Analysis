#ifndef MT2Analyzer_hh
#define MT2Analyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "MT2Analysis.hh"


class MT2Analyzer : public TreeAnalyzerBase {
public:
	MT2Analyzer(TTree *tree = 0);
	virtual ~MT2Analyzer();
	void BeginJob(TString filename="MassTree.root" , TString setofcuts="default",
	              bool isData=false, string data_PileUp="", string mc_PileUp="", string JEC="");
	void EndJob();
	void Loop();
	void SetMaxEvents(int a){fMaxEvents=a;};
	void SetProcessID(int ID){fID=ID;};
  	bool isS3;
  	bool noPU;  
	bool removePhoton;
        bool doPDF;
  bool isScan;
private:
	MT2Analysis             *fMT2Analysis;
  	int fMaxEvents;   
	int fID;
};
#endif
