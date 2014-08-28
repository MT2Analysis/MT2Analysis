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
  std::string file;
  float xsection;
  float nevents;
  float kfact;
  float PU_avg_weight;
  float lumi;
  int color;

};


struct MT2SampleBaby {

  std::string name;
  std::string sname;
  std::string dir;
  float nevents;
  int   id;
  float xsection;
  float PU_avg_weight;

};


struct MT2SampleBaby_basic {

  std::string name;
  std::string sname;
  std::string dir;
  std::string file;
  float xsection;
  float nevents;
  float kfact;
  float PU_avg_weight;
  float lumi;
  //int color;

};


struct MT2SamplePostBaby {

  std::string name;
  std::string sname;
  std::string dir;
  std::string file;
  int id;
  float nevents;
  float xsection;
  float lumi;
  float kfact;
  float filter;
  float scale1fb;
  float PU_avg_weight;

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

  static std::vector<MT2SampleBaby> loadSamplesBaby(const std::string& dir, const std::string& outputdir);
  static float getXSection( int sampleId );

  static std::vector<MT2SampleBaby_basic> loadSamplesBaby_basic(const std::string& filename);
  
  static std::vector<MT2SamplePostBaby> loadSamplesPostBaby(const std::string& filename);

 private:


};

#endif
