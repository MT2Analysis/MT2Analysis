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



MT2SingleLLEstimate MT2SingleLLEstimate::operator+( const MT2SingleLLEstimate& rhs ) {


  if( this->region->getName() != rhs.region->getName() ) {
    std::cout << "[MT2SingleLLEstimate::operator+] ERROR! Can't add MT2SingleLLEstimate with different MT2Regions!" << std::endl;
    exit(113);
  }

  std::string newname = this->name;
  if( rhs.name != this->name ) newname += "_" + rhs.name;


  MT2SingleLLEstimate result(name, *(this->region) );

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




MT2LeptonTypeLLEstimate MT2LeptonTypeLLEstimate::operator+( const MT2LeptonTypeLLEstimate& rhs ) {

  std::string newname = this->name;
  if( rhs.name != this->name ) newname += "_" + rhs.name;

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
    MT2SingleLLEstimate* gen_sum = new MT2SingleLLEstimate( *(this->gen[i]) + *(rhs.gen[i]) );

    reco.push_back( reco_sum );
    gen.push_back( gen_sum );

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




MT2LostLeptonEstimate MT2LostLeptonEstimate::operator+( const MT2LostLeptonEstimate& rhs ) {

  std::string newName = this->SName;
  if( rhs.SName != this->SName ) newName += "_" + rhs.SName;

  MT2LostLeptonEstimate result(newName);

  //for(std::map< std::string, MT2LeptonTypeLLEstimate*>::iterator i=l.begin(); i!=l.end(); ++i) 
  //  *(result.l[i->first]) = *(i->second) + *(rhs.l[i->first]);

  return result;

}

