#include "../interface/MT2LostLeptonUtilities.h"

#include <iostream>
#include <cstdlib>
#include <cmath>





// *******************************
//                               | 
//                               | 
//   MT2SingleLLEstimateBase     |
//                               |
//                               |
// *******************************



MT2SingleLLEstimateBase::MT2SingleLLEstimateBase( const MT2SingleLLEstimateBase& rhs ) {

  MT2SingleLLEstimateBase( rhs.name, *(rhs.region) );

}



MT2SingleLLEstimateBase::MT2SingleLLEstimateBase( const std::string& aname, const MT2Region& aregion ) {

  region = new MT2Region(aregion);
  name = aname;

  int nBins;
  double* bins;
  region->getBins(nBins, bins);

  yield = new TH1D(Form("yield_%s", name.c_str()), "", nBins, bins);
  yield->Sumw2();


}



MT2SingleLLEstimateBase::~MT2SingleLLEstimateBase() {

  delete region;
  delete yield;

}



void MT2SingleLLEstimateBase::addOverflow() {

  yield->SetBinContent(yield->GetNbinsX(),
      yield->GetBinContent(yield->GetNbinsX()  )+
      yield->GetBinContent(yield->GetNbinsX()+1)  );
  yield->SetBinError(  yield->GetNbinsX(),
      sqrt(yield->GetBinError(yield->GetNbinsX()  )*
           yield->GetBinError(yield->GetNbinsX()  )+
           yield->GetBinError(yield->GetNbinsX()+1)*
           yield->GetBinError(yield->GetNbinsX()+1)  ));

}




MT2SingleLLEstimateBase MT2SingleLLEstimateBase::operator+( const MT2SingleLLEstimateBase& rhs ) const {


  if( this->region->getName() != rhs.region->getName() ) {
    std::cout << "[MT2SingleLLEstimateBase::operator+] ERROR! Can't add MT2SingleLLEstimate with different MT2Regions!" << std::endl;
    exit(113);
  }

  std::string newname = this->name + "_" + rhs.name;

  MT2SingleLLEstimateBase result(newname, *(this->region) );

  result.yield->Add(this->yield);
  result.yield->Add(rhs.yield);

  return result;

}



// ***************************
//                           | 
//                           | 
//   MT2SingleLLEstimate     |
//                           |
//                           |
// ***************************




MT2SingleLLEstimate::MT2SingleLLEstimate( const MT2SingleLLEstimate& rhs ) : MT2SingleLLEstimateBase( rhs.name, *(rhs.region) ) {

  MT2SingleLLEstimate( rhs.name, *(rhs.region) );

}



MT2SingleLLEstimate::MT2SingleLLEstimate( const std::string& aname, const MT2Region& aregion ) : MT2SingleLLEstimateBase( aname, aregion ) {

  int nBins;
  double* bins;
  region->getBins(nBins, bins);


  yield_btagUp = new TH1D(Form("yield_btagUp_%s", name.c_str()), "", nBins, bins);
  yield_btagUp->Sumw2();

  yield_btagDown = new TH1D(Form("yield_btagDown_%s", name.c_str()), "", nBins, bins);
  yield_btagDown->Sumw2();


  float xMin = bins[0];
  float xMax = bins[nBins];


  // efficiencies are inclusive (one big bin):
  effLept_pass = new TH1D(Form("effLept_%s_pass", name.c_str()), "", 1, xMin, xMax );
  effLept_tot  = new TH1D(Form("effLept_%s_tot", name.c_str()), "", 1, xMin, xMax );

  effMT_pass = new TH1D(Form("effMT_%s_pass", name.c_str()), "", 1, xMin, xMax );
  effMT_tot  = new TH1D(Form("effMT_%s_tot", name.c_str()), "", 1, xMin, xMax );


}




MT2SingleLLEstimate::~MT2SingleLLEstimate() {

  delete yield_btagUp;
  delete yield_btagDown;
  delete effLept_pass;
  delete effLept_tot;
  delete effMT_pass;
  delete effMT_tot;

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

  result.yield_btagUp->Add(this->yield_btagUp);
  result.yield_btagUp->Add(rhs.yield_btagUp);

  result.yield_btagDown->Add(this->yield_btagDown);
  result.yield_btagDown->Add(rhs.yield_btagDown);

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




void MT2SingleLLEstimate::addOverflow() {

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

      MT2SingleLLEstimate* newPred = new MT2SingleLLEstimate( "pred_" + suffix, thisRegion );
      pred.push_back(newPred);

      MT2SingleLLEstimateBase* newsimtruth = new MT2SingleLLEstimateBase( "simtruth_" + suffix, thisRegion );
      simtruth.push_back(newsimtruth);

    } // for SR

  } // for HR

}




MT2SingleLLEstimate* MT2LeptonTypeLLEstimate::getRegion( const std::string& regionName ) const {


  MT2SingleLLEstimate* theRegion = 0;

  for( unsigned i=0; i<pred.size(); ++i ) {

    if( pred[i]->regionName() == regionName ) {

      theRegion = pred[i];
      break;

    }

  }


  return theRegion;

}




MT2SingleLLEstimateBase* MT2LeptonTypeLLEstimate::getRegionGen( const std::string& regionName ) const {


  MT2SingleLLEstimateBase* theRegion = 0;

  for( unsigned i=0; i<simtruth.size(); ++i ) {

    if( simtruth[i]->regionName() == regionName ) {

      theRegion = simtruth[i];
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

  int nRegions = this->pred.size();

  for( unsigned i=0; i<nRegions; ++i ) {

    if( (this->pred[i]->regionName() != rhs.pred[i]->regionName()) ) {
      std::cout << "[MT2LeptonTypeLLEstimate::operator+] ERROR! Can't add estimates with different regions. Exiting." << std::endl;
      exit(313);
    }

    MT2SingleLLEstimate* reco_sum = new MT2SingleLLEstimate( *(this->pred[i]) + *(rhs.pred[i]) );
    result.pred.push_back( reco_sum );

    MT2SingleLLEstimateBase* gen_sum = new MT2SingleLLEstimateBase( *(this->simtruth[i]) + *(rhs.simtruth[i]) );
    result.simtruth.push_back( gen_sum );

  }


  return result;

}




void MT2LeptonTypeLLEstimate::addOverflow() {

  for( unsigned i=0; i<pred.size(); ++i )
    pred[i]->addOverflow();
 
  for( unsigned i=0; i<simtruth.size(); ++i )
    simtruth[i]->addOverflow();
 
}



// ****************************
//                            | 
//                            | 
//   MT2LostLeptonEstimate    |
//                            |
//                            |
// ****************************


MT2LeptonTypeLLEstimate::~MT2LeptonTypeLLEstimate() {

  for( unsigned i=0; i<pred.size(); ++i ) {
    delete pred[i];
    pred[i] = 0;
  }

  for( unsigned i=0; i<simtruth.size(); ++i ) {
    delete simtruth[i];
    simtruth[i] = 0;
  }

}



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


  MT2LostLeptonEstimate result("");

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



/*
void MT2LostLeptonEstimate::writeTable( const std::string& fileName ) const {


  ofstream ofs(fileName.c_str());

  ofs << "\\begin{table}"              << endl
       << "\\begin{center}"             << endl
       << "\\small"                     << endl
       << "\\begin{tabular}{lccccccc}"  << endl
       << "\\hline\\hline"              << endl;
  
  else        ofs << "signal region       ";
  ofs  << " & $" << "N^{QCD}" << "$ & $" << "N^{Z}" << "$  & $" << "N^{W}" << "$  & $" << "N^{Top}" << "$  & $" << "N^{Other}" << "$ &           $"  << "N^{MC}" << "$      & $" << "N^{data}" << "$ ";
  ofs << "\\\\" << endl << "\\hline\\hline"             << endl;
  string oldlep = "dummy";
  string oldsr = "dummy";
  string oldHT = "dummy";

  for(unsigned int n = 0; n<sr.size(); ++n){

    string lep;
    if(lepr[n]==0) lep="Muo";
    else if(lepr[n]==1) lep="Ele";
    else continue;
    string sigreg;
    if(sr[n]==0) sigreg="2j, 0b";
    if(sr[n]==1) sigreg="2j, $\\ge 1$b";
    if(sr[n]==2) sigreg="$3-5$j, 0b";
    if(sr[n]==3) sigreg="$3-5$j, 1b";
    if(sr[n]==4) sigreg="$3-5$j, 2b";
    if(sr[n]==5) sigreg="$\\ge 6$j, 0b";
    if(sr[n]==6) sigreg="$\\ge 6$j, 1b";
    if(sr[n]==7) sigreg="$\\ge 6$j, 2b";
    if(sr[n]==8) sigreg="$\\ge 3$j, $\\ge 3$b";
    string htreg;
    if(htr[n]==0 && fMET) htreg = "450 GeV $\\leq H_{T} < 750$ GeV";
    if(htr[n]==1 && fHT) htreg = "750 GeV $\\leq H_{T} < 1200$ GeV";
    if(htr[n]==2 && fHT) htreg = "$H_{T}\\ge 1200$ GeV";
    if(htreg!=oldHT){
      ofs << " \\hline " << endl << "\\multicolumn{7}{l}{" <<  htreg << "} \\\\ " << endl << "\\hline" << endl;
      oldHT = htreg;
    }
    if(lep!=oldlep){
      ofs << " \\hline " << endl << lep << "\\\\ " << endl << "\\hline" << endl;
      oldlep = lep;
    }
    if(!fRebin){
    if(sigreg!=oldsr){
      ofs << " \\hline  " << endl << sigreg << "\\\\" << endl;
      oldsr = sigreg;
    }
    if(MT2up[n]==10000.) ofs << "$" << int(MT2low[n]) << "-" << "\\infty$" << " " << setw(4) << " & ";
    else            ofs << "$" << int(MT2low[n]) << "-" << int(MT2up[n]) << "$" << " " << setw(7) << " & ";
    }
    else ofs << " " << setw(18) << sigreg << " & ";
    ofs << fixed << setprecision(2)
    << " " << setw(7) << numQCD[n] << " & " << " " << setw(7) << numZ[n] << " & " << " " << setw(7) << numW[n] << " &  " << " " << setw(8) << numT[n] << " & " << " " << setw(10) << numOther[n] << " & " << " " << setw(10) 
    << numMC[n] << "$\\pm" << " " << setw(7) << BGerr[n] << "$ & ";
    ofs << " " << setw(10) << int(numData[n]) << " \\\\" << endl;
  }
  ofs << "\\hline\\hline"                                                                                            << endl
      << "\\end{tabular}"                                                                                                << endl
      << "\\end{center}"                                                                                                 << endl
      << "\\end{table}"                                                                                                  << endl;
  ofs << endl << endl;

  ofs.close();

}

*/
