#ifndef MT2RegionAnalysis_h
#define MT2RegionAnalysis_h



#include "MT2Region.h"
#include "TEfficiency.h"
#include "TH1D.h"
#include <map>





// this class for simtruth estimates (has only yield)

class MT2RegionAnalysisBase {

 public:

  MT2RegionAnalysisBase( const MT2RegionAnalysisBase& rhs );  
  MT2RegionAnalysisBase( const std::string& aname, const MT2Region& aregion );
  ~MT2RegionAnalysisBase();
 
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

  TH1F* yield;

  MT2RegionAnalysisBase operator+( const MT2RegionAnalysisBase& rhs ) const;
  void addOverflow();

  void writeTable( const std::string& fileName ) const;

 private:

};




// this class for b-tag up/down yields

class MT2RegionAnalysis : public MT2RegionAnalysisBase {

 public:

  MT2RegionAnalysis( const MT2RegionAnalysis& rhs );
  MT2RegionAnalysis( const std::string& aname, const MT2Region& aregion );
  ~MT2RegionAnalysis();


  TH1D* yield_btagUp;
  TH1D* yield_btagDown;


  MT2RegionAnalysis operator+( const MT2RegionAnalysis& rhs ) const;

  void addOverflow();
  
 private:

};



class MT2YieldPreAnalysis {

 public:

  MT2YieldPreAnalysis( const std::string& aSName ) {
    SName = aSName;
  }
  MT2YieldPreAnalysis( const std::string& aSName, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions );
  ~MT2YieldPreAnalysis();
  
  
  MT2RegionAnalysis* getRegion( const MT2Region& region ) const {
    return this->getRegion( region.getName() );
  }
  MT2RegionAnalysisBase* getRegionGen( const MT2Region& region ) const {
    return this->getRegionGen( region.getName() );
  }
  MT2RegionAnalysis* getRegion( const std::string& regionName ) const;
  MT2RegionAnalysisBase* getRegionGen( const std::string& regionName ) const;



  MT2YieldPreAnalysis operator+( const MT2YieldPreAnalysis& rhs ) const;
  
  void addOverflow();

  
  std::string SName;

  std::vector< MT2RegionAnalysis* > EventYield;
  std::vector< MT2RegionAnalysisBase* > simtruthEventYield;

 private:

};



class MT2YieldAnalysis {

 public:

  MT2YieldAnalysis( const std::string& aSName="") {
    SName = aSName;
  }
  ~MT2YieldAnalysis() {};

  void add(const MT2YieldAnalysis& rhs);

  MT2YieldAnalysis operator+(const MT2YieldAnalysis& rhs) const;
  const MT2YieldAnalysis& operator=(const MT2YieldAnalysis& rhs);


  std::string SName;
  std::map< std::string, MT2YieldPreAnalysis*> f; // "fake" map


 private:

};


#endif
