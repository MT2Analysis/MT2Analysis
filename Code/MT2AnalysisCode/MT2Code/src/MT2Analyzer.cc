#include "MT2Analyzer.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "MT2Analysis.hh"

using namespace std;

MT2Analyzer::MT2Analyzer(TTree *tree) : TreeAnalyzerBase(tree) {
	fMT2Analysis             = new MT2Analysis(fTR);
	Util::SetStyle();
}

MT2Analyzer::~MT2Analyzer(){
	delete fMT2Analysis;
	if(!fTR->fChain) cout << "MT2Analyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void MT2Analyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	
	if(fMaxEvents==-1)nentries=fTR->GetEntries();
	if(fMaxEvents > 0){
		nentries=min((Long64_t)fMaxEvents, fTR->GetEntries());
		cout << " only running on first " << nentries << " events" << endl;
	}
	// loop over all ntuple entries
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		fTR->GetEntry(jentry);
                if ( fCurRun != fTR->Run ) {
        		fCurRun = fTR->Run;
			fMT2Analysis        ->BeginRun(fCurRun);
                  	skipRun = false; // re-initialize
                  	if ( !CheckRun() ) skipRun = true;
                }
                // Check if new lumi is in JSON file
                if ( !skipRun && fCurLumi != fTR->LumiSection ) {
                  	fCurLumi = fTR->LumiSection;
                  	skipLumi = false; // Re-initialise
                  	if ( !CheckRunLumi() ) skipLumi = true;
                }
		if ( !(skipRun || skipLumi) ) {
			fMT2Analysis        ->Analyze();
		}
	}
}

// Method called before starting the event loop
void MT2Analyzer::BeginJob(TString filename, TString setofcuts, bool isData, string data_PileUp, string mc_PileUp){
	fMT2Analysis             ->ReadCuts(setofcuts);
	fMT2Analysis             ->SetType(isData);
	fMT2Analysis             ->SetPileUpSrc(data_PileUp, mc_PileUp);
	fMT2Analysis             ->SetOutputDir(fOutputDir);
	fMT2Analysis             ->fVerbose        = fVerbose;
	fMT2Analysis             ->Begin(filename);

        fMT2Analysis             ->isS3        = isS3;
	fMT2Analysis             ->noPU        = noPU;
}

// Method called after finishing the event loop
void MT2Analyzer::EndJob(){
  fMT2Analysis         ->End();
  cout << " MT2Analyzer::End()                                             " << endl;
  
}
