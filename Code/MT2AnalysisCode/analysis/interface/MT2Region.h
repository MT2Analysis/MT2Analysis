#ifndef MT2Region_h
#define MT2Region_h

#include <string>





class MT2HTRegion {

 public:

  MT2HTRegion( const MT2HTRegion& rhs );
  MT2HTRegion( const std::string& aname, float ahtMin, float ahtMax, float ametMin, const std::string& aHLT_selection="" );
  ~MT2HTRegion() {};


  // univocal identifier:
  std::string getName() const;


  std::string name;
  float htMin;
  float htMax;
  float metMin;

  std::string HLT_selection;

 private:

};





class MT2SignalRegion {

 public:

  MT2SignalRegion( int njmin, int njmax, int nbmin, int nbmax );
  MT2SignalRegion( const MT2SignalRegion& rhs );

  ~MT2SignalRegion() {};

  // univocal identifier:
  std::string getName() const;

  std::string getSingleSignalRegionString( const std::string& suffix, int n_min , int n_max=-1 ) const;
  
  int nJetsMin; 
  int nJetsMax;
  int nBJetsMin;
  int nBJetsMax;


 private:

};





class MT2Region {

 public:

  MT2Region( const MT2Region& region ) {
    htRegion_ = new MT2HTRegion(*(region.htRegion()));
    sigRegion_ = new MT2SignalRegion(*(region.sigRegion()));
  }

  MT2Region( MT2HTRegion* htRegion, MT2SignalRegion* sigRegion ) {
    htRegion_ = new MT2HTRegion(*htRegion);
    sigRegion_ = new MT2SignalRegion(*sigRegion);
  }
  ~MT2Region() {};


  // this is a univocal identifier (same jet/HT thresholds always correspond to same name)
  std::string getName() const {
    return htRegion_->getName() + "_" + sigRegion_->getName();
  }

  void getBins( int& nBins, double*& bins ) const;

  MT2HTRegion* htRegion() const {
    return htRegion_;
  }

  MT2SignalRegion* sigRegion() const {
    return sigRegion_;
  }

  float htMin() { return htRegion_->htMin; };
  float htMax() { return htRegion_->htMax; };
  int nJetsMin() { return sigRegion_->nJetsMin; };
  int nJetsMax() { return sigRegion_->nJetsMax; };
  int nBJetsMin() { return sigRegion_->nBJetsMin; };
  int nBJetsMax() { return sigRegion_->nBJetsMax; };



 private:
  
  MT2HTRegion* htRegion_;
  MT2SignalRegion* sigRegion_;

};





#endif
