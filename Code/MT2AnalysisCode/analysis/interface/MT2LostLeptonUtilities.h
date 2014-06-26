#ifndef MT2LostLeptonUtilities_h
#define MT2LostLeptonUtilities_h


#include "MT2Region.h"
#include "TEfficiency.h"
#include "TH1D.h"
#include <map>



class MT2SingleLLEstimate {

 public:

  MT2SingleLLEstimate( const MT2SingleLLEstimate& rhs );
  MT2SingleLLEstimate( const std::string& aname, const MT2Region& aregion );
  ~MT2SingleLLEstimate() {};


  // this is just a name to differentiate different
  // instances of the same class
  std::string name;

  void getBins( int& nBins, double* bins ) const {
    return region->getBins(nBins, bins);
  }

  MT2HTRegion* htRegion() const {
    return region->htRegion();
  }
  MT2SignalRegion* sigRegion() const {
    return region->sigRegion();
  }

  // this is a univocal identifier of the region
  // regions with the same definition (jet numbers and HT cuts)
  // have the same name
  std::string regionName() const {
    return region->getName();
  }

  MT2Region* region;

  TH1D* yield;

  TH1D* effLept_pass;
  TH1D* effLept_tot;
  TH1D* effMT_pass;
  TH1D* effMT_tot;


  MT2SingleLLEstimate operator+( const MT2SingleLLEstimate& rhs ) const;

  TEfficiency* effLept() {
    return new TEfficiency(*effLept_pass, *effLept_tot);
  }

  TEfficiency* effMT() {
    return new TEfficiency(*effMT_pass, *effMT_tot);
  }

 private:

};



class MT2LeptonTypeLLEstimate {

 public:

  MT2LeptonTypeLLEstimate( const std::string& aname, const std::string& aSName) {
    name = aname;
    SName = aSName;
  }
  MT2LeptonTypeLLEstimate( const std::string& aname, const std::string& aSName, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions );
  ~MT2LeptonTypeLLEstimate() {};


  MT2SingleLLEstimate* getRecoRegion( const std::string& regionName ) const;
  MT2SingleLLEstimate* getGenRegion( const std::string& regionName ) const;

  MT2LeptonTypeLLEstimate operator+( const MT2LeptonTypeLLEstimate& rhs ) const;

  std::string name;
  std::string SName;

  std::vector< MT2SingleLLEstimate* > reco;
  std::vector< MT2SingleLLEstimate* > gen;

 private:

};



class MT2LostLeptonEstimate {

 public:

  MT2LostLeptonEstimate( const std::string& aname, const std::string& aSName="") {
    name = aname;
    SName = aSName;
  }
  ~MT2LostLeptonEstimate() {};

  void add(const MT2LostLeptonEstimate& rhs);

  MT2LostLeptonEstimate operator+(const MT2LostLeptonEstimate& rhs) const;
  const MT2LostLeptonEstimate& operator=(const MT2LostLeptonEstimate& rhs);

  std::string name;
  std::string SName;
  std::map< std::string, MT2LeptonTypeLLEstimate*> l; // "Ele" and "Muo"

 private:

};




#endif
