#include "../interface/MT2RegionAnalysisUtilities.h"

#include <iostream>
#include <cstdlib>
#include <cmath>





// *******************************
//                               | 
//                               | 
//     MT2RegionAnalysisBase     |
//                               |
//                               |
// *******************************



MT2RegionAnalysisBase::MT2RegionAnalysisBase( const MT2RegionAnalysisBase& rhs ) {

  MT2RegionAnalysisBase( rhs.name, *(rhs.region) );

}



MT2RegionAnalysisBase::MT2RegionAnalysisBase( const std::string& aname, const MT2Region& aregion ) {

  region = new MT2Region(aregion);
  name = aname;

  int nBins;
  double* bins;
  region->getBins(nBins, bins);

  yield = new TH1D(Form("yield_%s", name.c_str()), "", nBins, bins);
  yield->Sumw2();


}



MT2RegionAnalysisBase::~MT2RegionAnalysisBase() {

  delete region;
  delete yield;

}



void MT2RegionAnalysisBase::addOverflow() {

  yield->SetBinContent(yield->GetNbinsX(),
      yield->GetBinContent(yield->GetNbinsX()  )+
      yield->GetBinContent(yield->GetNbinsX()+1)  );
  yield->SetBinError(  yield->GetNbinsX(),
      sqrt(yield->GetBinError(yield->GetNbinsX()  )*
           yield->GetBinError(yield->GetNbinsX()  )+
           yield->GetBinError(yield->GetNbinsX()+1)*
           yield->GetBinError(yield->GetNbinsX()+1)  ));

}




MT2RegionAnalysisBase MT2RegionAnalysisBase::operator+( const MT2RegionAnalysisBase& rhs ) const {


  if( this->region->getName() != rhs.region->getName() ) {
    std::cout << "[MT2RegionAnalysisBase::operator+] ERROR! Can't add MT2RegionAnalysis with different MT2Regions!" << std::endl;
    exit(113);
  }

  std::string newname = this->name + "_" + rhs.name;

  MT2RegionAnalysisBase result(newname, *(this->region) );

  result.yield->Add(this->yield);
  result.yield->Add(rhs.yield);

  return result;

}



// ***************************
//                           | 
//                           | 
//     MT2RegionAnalysis     |
//                           |
//                           |
// ***************************




MT2RegionAnalysis::MT2RegionAnalysis( const MT2RegionAnalysis& rhs ) : MT2RegionAnalysisBase( rhs.name, *(rhs.region) ) {

  MT2RegionAnalysis( rhs.name, *(rhs.region) );

}



MT2RegionAnalysis::MT2RegionAnalysis( const std::string& aname, const MT2Region& aregion ) : MT2RegionAnalysisBase( aname, aregion ) {

  int nBins;
  double* bins;
  region->getBins(nBins, bins);


  yield_btagUp = new TH1D(Form("yield_btagUp_%s", name.c_str()), "", nBins, bins);
  yield_btagUp->Sumw2();

  yield_btagDown = new TH1D(Form("yield_btagDown_%s", name.c_str()), "", nBins, bins);
  yield_btagDown->Sumw2();


}




MT2RegionAnalysis::~MT2RegionAnalysis() {

  delete yield_btagUp;
  delete yield_btagDown;

}




MT2RegionAnalysis MT2RegionAnalysis::operator+( const MT2RegionAnalysis& rhs ) const {


  if( this->region->getName() != rhs.region->getName() ) {
    std::cout << "[MT2RegionAnalysis::operator+] ERROR! Can't add MT2RegionAnalysis with different MT2Regions!" << std::endl;
    exit(113);
  }

  std::string newname = this->name + "_" + rhs.name;

  MT2RegionAnalysis result(newname, *(this->region) );

  result.yield->Add(this->yield);
  result.yield->Add(rhs.yield);

  result.yield_btagUp->Add(this->yield_btagUp);
  result.yield_btagUp->Add(rhs.yield_btagUp);

  result.yield_btagDown->Add(this->yield_btagDown);
  result.yield_btagDown->Add(rhs.yield_btagDown);

  return result;

}




void MT2RegionAnalysis::addOverflow() {

  yield->SetBinContent(yield->GetNbinsX(),
      yield->GetBinContent(yield->GetNbinsX()  )+
      yield->GetBinContent(yield->GetNbinsX()+1)  );
  yield->SetBinError(  yield->GetNbinsX(),
      sqrt(yield->GetBinError(yield->GetNbinsX()  )*
           yield->GetBinError(yield->GetNbinsX()  )+
           yield->GetBinError(yield->GetNbinsX()+1)*
           yield->GetBinError(yield->GetNbinsX()+1)  ));

  yield_btagUp->SetBinContent(yield_btagUp->GetNbinsX(),
      yield_btagUp->GetBinContent(yield_btagUp->GetNbinsX()  )+
      yield_btagUp->GetBinContent(yield_btagUp->GetNbinsX()+1)  );
  yield_btagUp->SetBinError(  yield_btagUp->GetNbinsX(),
      sqrt(yield_btagUp->GetBinError(yield_btagUp->GetNbinsX()  )*
           yield_btagUp->GetBinError(yield_btagUp->GetNbinsX()  )+
           yield_btagUp->GetBinError(yield_btagUp->GetNbinsX()+1)*
           yield_btagUp->GetBinError(yield_btagUp->GetNbinsX()+1)  ));

  yield_btagDown->SetBinContent(yield_btagDown->GetNbinsX(),
      yield_btagDown->GetBinContent(yield_btagDown->GetNbinsX()  )+
      yield_btagDown->GetBinContent(yield_btagDown->GetNbinsX()+1)  );
  yield_btagDown->SetBinError(  yield_btagDown->GetNbinsX(),
      sqrt(yield_btagDown->GetBinError(yield_btagDown->GetNbinsX()  )*
           yield_btagDown->GetBinError(yield_btagDown->GetNbinsX()  )+
           yield_btagDown->GetBinError(yield_btagDown->GetNbinsX()+1)*
           yield_btagDown->GetBinError(yield_btagDown->GetNbinsX()+1)  ));

}






// *****************************
//                             | 
//                             | 
//       MT2YieldPreAnalysis         |
//                             |
//                             |
// *****************************



MT2YieldPreAnalysis::MT2YieldPreAnalysis( const std::string& aSName, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions ) {

  SName = aSName;

  for( unsigned iHR=0; iHR<HTRegions.size(); ++iHR ) {

    for( unsigned iSR=0; iSR<signalRegions.size(); ++iSR ) {

      std::string suffix = SName + "_" + HTRegions[iHR].name + "_" + signalRegions[iSR].getName();

      MT2Region thisRegion(&(HTRegions[iHR]), &(signalRegions[iSR]));

      MT2RegionAnalysis* newEventYield = new MT2RegionAnalysis( "EventYield_" + suffix, thisRegion );
      EventYield.push_back(newEventYield);

      MT2RegionAnalysisBase* newsimtruthEventYield = new MT2RegionAnalysisBase( "simtruthEventYield_" + suffix, thisRegion );
      simtruthEventYield.push_back(newsimtruthEventYield);

    } // for SR

  } // for HR

}



MT2YieldPreAnalysis::~MT2YieldPreAnalysis() {

  for( unsigned i=0; i<EventYield.size(); ++i ) {
    delete EventYield[i];
    EventYield[i] = 0;
  }

  for( unsigned i=0; i<simtruthEventYield.size(); ++i ) {
    delete simtruthEventYield[i];
    simtruthEventYield[i] = 0;
  }

}


MT2RegionAnalysis* MT2YieldPreAnalysis::getRegion( const std::string& regionName ) const {


  MT2RegionAnalysis* theRegion = 0;

  for( unsigned i=0; i<EventYield.size(); ++i ) {

    if( EventYield[i]->regionName() == regionName ) {

      theRegion = EventYield[i];
      break;

    }

  }


  return theRegion;

}



MT2RegionAnalysisBase* MT2YieldPreAnalysis::getRegionGen( const std::string& regionName ) const {


  MT2RegionAnalysisBase* theRegion = 0;

  for( unsigned i=0; i<simtruthEventYield.size(); ++i ) {

    if( simtruthEventYield[i]->regionName() == regionName ) {

      theRegion = simtruthEventYield[i];
      break;

    }

  }


  return theRegion;

}



MT2YieldPreAnalysis MT2YieldPreAnalysis::operator+( const MT2YieldPreAnalysis& rhs ) const {

  
  std::string newSName = this->SName;
  if( rhs.SName != this->SName ) newSName += "_" + rhs.SName;

  MT2YieldPreAnalysis result( newSName );

  int nRegions = this->EventYield.size();

  for( unsigned i=0; i<nRegions; ++i ) {

    if( (this->EventYield[i]->regionName() != rhs.EventYield[i]->regionName()) ) {
      std::cout << "[MT2YieldPreAnalysis::operator+] ERROR! Can't add yields with different regions. Exiting." << std::endl;
      exit(313);
    }

    MT2RegionAnalysis* reco_sum = new MT2RegionAnalysis( *(this->EventYield[i]) + *(rhs.EventYield[i]) );
    result.EventYield.push_back( reco_sum );

    MT2RegionAnalysisBase* gen_sum = new MT2RegionAnalysisBase( *(this->simtruthEventYield[i]) + *(rhs.simtruthEventYield[i]) );
    result.simtruthEventYield.push_back( gen_sum );

  }


  return result;

}



void MT2YieldPreAnalysis::addOverflow() {

  for( unsigned i=0; i<EventYield.size(); ++i )
    EventYield[i]->addOverflow();
 
  for( unsigned i=0; i<simtruthEventYield.size(); ++i )
    simtruthEventYield[i]->addOverflow();
 
}




// ****************************
//                            |
//                            |
//      MT2YieldAnalysis      |
//                            |
//                            |
// ****************************                                                                                                                                                                                                               


const MT2YieldAnalysis& MT2YieldAnalysis::operator=( const MT2YieldAnalysis& rhs ) {

  SName = rhs.SName;

  for(std::map< std::string, MT2YieldPreAnalysis*>::const_iterator j=rhs.f.begin(); j!=rhs.f.end(); ++j)
    f[j->first] = new MT2YieldPreAnalysis(*(j->second));
  
  return *this;
  
}


void MT2YieldAnalysis::add(const MT2YieldAnalysis& rhs) {

  (*this) = (*this) + rhs;

}


MT2YieldAnalysis MT2YieldAnalysis::operator+( const MT2YieldAnalysis& rhs ) const {


  MT2YieldAnalysis result("");

  if( f.size()==0 ) {

    result = rhs;

  } else if( rhs.f.size()==0 ) {

    result = *this;

  } else {

    std::string newSName = this->SName;
    if( rhs.SName != this->SName ) {
      if( newSName!="" && rhs.SName!="" ) newSName += "_" + rhs.SName;
      else if( rhs.SName!="" ) newSName = rhs.SName;
    }
    result.SName = newSName;

    for(std::map< std::string, MT2YieldPreAnalysis*>::const_iterator i=f.begin(); i!=f.end(); ++i) {
      for(std::map< std::string, MT2YieldPreAnalysis*>::const_iterator j=rhs.f.begin(); j!=rhs.f.end(); ++j) {
        if( i->first != j->first ) continue;
        result.f[i->first] = new MT2YieldPreAnalysis(*(i->second) + *(j->second));
      }
    }

  } // else                                                                                                                                                                                                                                   



  return result;

}
