#ifndef MT2LostLeptonUtilities_h
#define MT2LostLeptonUtilities_h


#include "MT2Region.h"
#include "TEfficiency.h"
#include "TH1D.h"
#include <map>






// this class for simtruth estimates (has only yield)

class MT2SingleLLEstimateBase {

 public:

  MT2SingleLLEstimateBase( const MT2SingleLLEstimateBase& rhs );  
  MT2SingleLLEstimateBase( const std::string& aname, const MT2Region& aregion );
  ~MT2SingleLLEstimateBase();
 
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

  MT2SingleLLEstimateBase operator+( const MT2SingleLLEstimateBase& rhs ) const;
  void addOverflow();

  void writeTable( const std::string& fileName ) const;

 private:

};




// this class for reco estimates (has also efficiencies)

class MT2SingleLLEstimate : public MT2SingleLLEstimateBase {

 public:

  MT2SingleLLEstimate( const MT2SingleLLEstimate& rhs );
  MT2SingleLLEstimate( const std::string& aname, const MT2Region& aregion );
  ~MT2SingleLLEstimate();


  TH1D* yield_btagUp;
  TH1D* yield_btagDown;

  TH1D* effLept_pass;
  TH1D* effLept_tot;
  TH1D* effMT_pass;
  TH1D* effMT_tot;


  MT2SingleLLEstimate operator+( const MT2SingleLLEstimate& rhs ) const;

  void addOverflow();

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
  ~MT2LeptonTypeLLEstimate();


  MT2SingleLLEstimate* getRegion( const MT2Region& region ) const {
    return this->getRegion( region.getName() );
  }
  MT2SingleLLEstimateBase* getRegionGen( const MT2Region& region ) const {
    return this->getRegionGen( region.getName() );
  }
  MT2SingleLLEstimate* getRegion( const std::string& regionName ) const;
  MT2SingleLLEstimateBase* getRegionGen( const std::string& regionName ) const;

  MT2LeptonTypeLLEstimate operator+( const MT2LeptonTypeLLEstimate& rhs ) const;

  void addOverflow();

  std::string name;
  std::string SName;

  std::vector< MT2SingleLLEstimate* > pred;
  std::vector< MT2SingleLLEstimateBase* > simtruth;

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
