#include "../interface/MT2LostLeptonUtilities.h"

#include <iostream>
#include <cstdlib>





// ***************************
//                           | 
//                           | 
//   MT2SingleLLEstimate     |
//                           |
//                           |
// ***************************


MT2SingleLLEstimate::MT2SingleLLEstimate( const MT2SingleLLEstimate& rhs ) {

  MT2SingleLLEstimate( rhs.name, *(rhs.region) );

}



MT2SingleLLEstimate::MT2SingleLLEstimate( const std::string& aname, const MT2Region& aregion ) {

  region = new MT2Region(aregion);
  name = aname;

  int nBins;
  double* bins;
  region->getBins(nBins, bins);

  yield = new TH1D(Form("yield_%s", name.c_str()), "", nBins, bins);
  yield->Sumw2();

  effLept_pass = new TH1D(Form("effLept_%s_pass", name.c_str()), "", 1, 0., 1000000.);
  effLept_pass->Sumw2();
  effLept_tot  = new TH1D(Form("effLept_%s_tot", name.c_str()), "", 1, 0., 1000000.);
  effLept_tot->Sumw2();

  effMT_pass = new TH1D(Form("effMT_%s_pass", name.c_str()), "", 1, 0., 1000000.);
  effMT_pass->Sumw2();
  effMT_tot  = new TH1D(Form("effMT_%s_tot", name.c_str()), "", 1, 0., 1000000.);
  effMT_tot->Sumw2();

}



MT2SingleLLEstimate MT2SingleLLEstimate::operator+( const MT2SingleLLEstimate& rhs ) const {


  if( this->region->getName() != rhs.region->getName() ) {
    std::cout << "[MT2SingleLLEstimate::operator+] ERROR! Can't add MT2SingleLLEstimate with different MT2Regions!" << std::endl;
    exit(113);
  }

  std::string newname = this->name + "_" + rhs.name;

  MT2SingleLLEstimate result(newname, *(this->region) );

  result.yield->Add(this->yield);
  result.yield->Add(rhs.yield);

  result.effLept_pass->Add(this->effLept_pass);
  result.effLept_pass->Add(rhs.effLept_pass);

  result.effLept_tot->Add(this->effLept_tot);
  result.effLept_tot->Add(rhs.effLept_tot);

  result.effMT_pass->Add(this->effMT_pass);
  result.effMT_pass->Add(rhs.effMT_pass);

  result.effMT_tot->Add(this->effMT_tot);
  result.effMT_tot->Add(rhs.effMT_tot);


  return result;

}







// *****************************
//                             | 
//                             | 
//   MT2LeptonTypeLLEstimate   |
//                             |
//                             |
// *****************************



MT2LeptonTypeLLEstimate::MT2LeptonTypeLLEstimate( const std::string& aname, const std::string& aSName, std::vector<MT2HTRegion> HTRegions, std::vector<MT2SignalRegion> signalRegions ) {

  name = aname;
  SName = aSName;

  for( unsigned iHR=0; iHR<HTRegions.size(); ++iHR ) {

    for( unsigned iSR=0; iSR<signalRegions.size(); ++iSR ) {

      std::string suffix = name + "_" + SName + "_" + HTRegions[iHR].name + "_" + signalRegions[iSR].getName();

      MT2Region thisRegion(&(HTRegions[iHR]), &(signalRegions[iSR]));

      MT2SingleLLEstimate* recoEst = new MT2SingleLLEstimate( "reco_" + suffix, thisRegion );
      MT2SingleLLEstimate* genEst  = new MT2SingleLLEstimate( "gen_"  + suffix, thisRegion );

      reco.push_back(recoEst);
       gen.push_back(genEst);

    } // for SR

  } // for HR

}




MT2SingleLLEstimate* MT2LeptonTypeLLEstimate::getRecoRegion( const std::string& regionName ) const {


  MT2SingleLLEstimate* theRegion = 0;

  for( unsigned i=0; i<reco.size(); ++i ) {

    if( reco[i]->regionName() == regionName ) {

      theRegion = reco[i];
      break;

    }

  }


  return theRegion;

}





MT2SingleLLEstimate* MT2LeptonTypeLLEstimate::getGenRegion( const std::string& regionName ) const {


  MT2SingleLLEstimate* theRegion;

  for( unsigned i=0; i<gen.size(); ++i ) {

    if( gen[i]->regionName() == regionName ) {

      theRegion = gen[i];
      break;

    }

  }


  return theRegion;

}




MT2LeptonTypeLLEstimate MT2LeptonTypeLLEstimate::operator+( const MT2LeptonTypeLLEstimate& rhs ) const {

  std::string newname = this->name + "_" + rhs.name;

  std::string newSName = this->SName;
  if( rhs.SName != this->SName ) newSName += "_" + rhs.SName;

  MT2LeptonTypeLLEstimate result( newname, newSName );

  int nRegions = this->reco.size();

  for( unsigned i=0; i<nRegions; ++i ) {

    if( (this->reco[i]->regionName() != rhs.reco[i]->regionName()) || (this->gen[i]->regionName() != rhs.gen[i]->regionName())) {
      std::cout << "[MT2LeptonTypeLLEstimate::operator+] ERROR! Can't add estimates with different regions. Exiting." << std::endl;
      exit(313);
    }

    MT2SingleLLEstimate* reco_sum = new MT2SingleLLEstimate( *(this->reco[i]) + *(rhs.reco[i]) );
    MT2SingleLLEstimate* gen_sum  = new MT2SingleLLEstimate( *(this->gen[i])  + *(rhs.gen[i]) );

    result.reco.push_back( reco_sum );
    result.gen.push_back( gen_sum );

  }


  return result;

}







// ****************************
//                            | 
//                            | 
//   MT2LostLeptonEstimate    |
//                            |
//                            |
// ****************************




const MT2LostLeptonEstimate& MT2LostLeptonEstimate::operator=( const MT2LostLeptonEstimate& rhs ) {

  name = rhs.name;
  SName = rhs.SName;

  for(std::map< std::string, MT2LeptonTypeLLEstimate*>::const_iterator j=rhs.l.begin(); j!=rhs.l.end(); ++j) 
    l[j->first] = new MT2LeptonTypeLLEstimate(*(j->second));
 

  return *this;

}


void MT2LostLeptonEstimate::add(const MT2LostLeptonEstimate& rhs) {

  (*this) = (*this) + rhs;

}


MT2LostLeptonEstimate MT2LostLeptonEstimate::operator+( const MT2LostLeptonEstimate& rhs ) const {


  MT2LostLeptonEstimate result("tmp");

  if( l.size()==0 ) {

    result = rhs;

  } else if( rhs.l.size()==0 ) {

    result = *this;

  } else {

    std::string newName = this->name + "_" + rhs.name;
    result.name = newName;

    std::string newSName = this->SName;
    if( rhs.SName != this->SName ) {
      if( newSName!="" && rhs.SName!="" ) newSName += "_" + rhs.SName;
      else if( rhs.SName!="" ) newSName = rhs.SName;
    }
    result.SName = newSName;

    for(std::map< std::string, MT2LeptonTypeLLEstimate*>::const_iterator i=l.begin(); i!=l.end(); ++i) {
      for(std::map< std::string, MT2LeptonTypeLLEstimate*>::const_iterator j=rhs.l.begin(); j!=rhs.l.end(); ++j) {
        if( i->first != j->first ) continue;
        result.l[i->first] = new MT2LeptonTypeLLEstimate(*(i->second) + *(j->second));
      }
    }

  } // else



  return result;

}

