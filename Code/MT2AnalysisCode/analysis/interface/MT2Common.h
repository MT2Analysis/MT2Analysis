#ifndef MT2Common_h
#define MT2Common_h

#include <string>

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"


struct MT2Sample {

  std::string name;
  std::string sname;
  std::string shapename;
  std::string type;
  TFile *file;
  TTree *tree;
  float xsection;
  float nevents;
  float kfact;
  float PU_avg_weight;
  float lumi;
  int color;

};





class MT2Common {


 public:

  MT2Common();
  ~MT2Common();

  static void getBins( float ht_min, int njet_min, int njet_max, int nbjet_min, int nbjet_max, int &nBins, Double_t* bins);
  static void getBins( const std::string& signal_region, int &nBins, Double_t* bins);

  static std::string getSignalRegion( float ht_min, int njet_min, int njet_max, int nbjet_min, int nbjet_max);
  static std::string getSingleSignalRegionString( int n_min , int n_max, const std::string& suffix );

  static std::vector<MT2Sample> loadSamples(const std::string& filename);
  static std::string getSampleType( MT2Sample sample );

 private:



};

#endif
