#include "../interface/MT2Region.h"

#include <iostream>
#include <stdio.h>
#include <cstdlib>




MT2HTRegion::MT2HTRegion( const MT2HTRegion& rhs ) {

  name = rhs.name;
  htMin = rhs.htMin;
  htMax = rhs.htMax;
  metMin = rhs.metMin;
  HLT_selection = rhs.HLT_selection;
  
}



MT2HTRegion::MT2HTRegion( const std::string& aname, float ahtMin, float ahtMax, float ametMin, const std::string& aHLT_selection ) {

  name = aname;
  htMin = ahtMin;
  htMax = ahtMax;
  metMin = ametMin;
  HLT_selection = aHLT_selection;

}


std::string MT2HTRegion::getName() const {

  char n[512];
  sprintf( n, "HTge%.0f", htMin );
  std::string n_str(n);

  return n_str;

}




MT2SignalRegion::MT2SignalRegion(int njmin, int njmax, int nbmin, int nbmax ) {

  nJetsMin = njmin;
  nJetsMax = njmax;
  nBJetsMin = nbmin;
  nBJetsMax = nbmax;

}




MT2SignalRegion::MT2SignalRegion( const MT2SignalRegion& rhs ) {

  nJetsMin = rhs.nJetsMin;
  nJetsMax = rhs.nJetsMax;
  nBJetsMin = rhs.nBJetsMin;
  nBJetsMax = rhs.nBJetsMax;

}



std::string MT2SignalRegion::getName() const {
 
  std::string jString = getSingleSignalRegionString( "j", nJetsMin, nJetsMax   );
  std::string bString = getSingleSignalRegionString( "b", nBJetsMin, nBJetsMax );

  std::string signal_region = jString + bString;

  return signal_region;

}



std::string MT2SignalRegion::getSingleSignalRegionString( const std::string& suffix, int n_min , int n_max ) const {

  if( n_min==-1 && n_max==-1 ) return std::string("");
  if( n_max==-1 ) n_max=n_min;

  char signal_region_ch[64];
  if( n_max!=n_min )
    sprintf( signal_region_ch, "%dto%d%s", n_min, n_max, suffix.c_str() );
  else 
    sprintf( signal_region_ch, "%d%s", n_min, suffix.c_str() );
  
  std::string signal_region(signal_region_ch);

  return signal_region;

}

 


void MT2Region::getBins( int &nBins, double*& bins) const {

  std::string regionName = this->getName();

  if( regionName == "HTge1200_2j0b" ) {
    const int nBins_tmp                      = 6;
    bins = new double[nBins_tmp+1]{120, 150, 200, 260, 350, 550, 900};
    nBins = nBins_tmp;
  } else if( regionName == "HTge1200_2j1to2b" ) {
    const int nBins_tmp                      = 2;
    bins = new double[nBins_tmp+1]{100, 180, 350};
    nBins = nBins_tmp;
  } else if( regionName == "HTge1200_3to5j0b" ) {
    const int nBins_tmp                      = 7;
    bins = new double[nBins_tmp+1]{160, 185, 220, 270, 350, 450, 650, 1000};
    nBins = nBins_tmp;
  } else if( regionName == "HTge1200_3to5j1b" ) {
    const int nBins_tmp                      = 4;
    bins = new double[nBins_tmp+1]{150, 180, 230, 350, 550};
    nBins = nBins_tmp;
  } else if( regionName == "HTge1200_3to5j2b" ) {
    const int nBins_tmp                      = 2;
    bins = new double[nBins_tmp+1]{130, 200, 350};
    nBins = nBins_tmp;
  } else if( regionName == "HTge1200_6j0b" ) {
    const int nBins_tmp                      = 3;
    bins = new double[nBins_tmp+1]{160, 200, 300, 500};
    nBins = nBins_tmp;
  } else if( regionName == "HTge1200_6j1b" ) {
    const int nBins_tmp                      = 3;
    bins = new double[nBins_tmp+1]{150, 200, 300, 500};
    nBins = nBins_tmp;
  } else if( regionName == "HTge1200_6j2b" ) {
    const int nBins_tmp                      = 2;
    bins = new double[nBins_tmp+1]{130, 200, 350};
    nBins = nBins_tmp;
  } else if( regionName == "HTge1200_3b" ) {
    const int nBins_tmp                      = 1;
    bins = new double[nBins_tmp+1]{125, 300};
    nBins = nBins_tmp;

  // HT > 750 && HT < 1200
  } else if( regionName == "HTge750_2j0b" ) {
    const int nBins_tmp                      = 9;
    bins = new double[nBins_tmp+1]{125, 150, 180, 220, 270, 325, 425, 580, 780, 1000};
    nBins = nBins_tmp;
  } else if( regionName == "HTge750_2j1to2b" ) {
    const int nBins_tmp                      = 5;
    bins = new double[nBins_tmp+1]{100, 135, 170, 260, 450, 700};
    nBins = nBins_tmp;
  } else if( regionName == "HTge750_3to5j0b" ) {
    const int nBins_tmp                      = 9;
    bins = new double[nBins_tmp+1]{160, 185, 215, 250, 300, 370, 480, 640, 800, 1000};
    nBins = nBins_tmp;
  } else if( regionName == "HTge750_3to5j1b" ) {
    const int nBins_tmp                      = 6;
    bins = new double[nBins_tmp+1]{150, 175, 210, 270, 380, 600, 900};
    nBins = nBins_tmp;
  } else if( regionName == "HTge750_3to5j2b" ) {
    const int nBins_tmp                      = 5;
    bins = new double[nBins_tmp+1]{130, 160, 200, 270, 370, 500};
    nBins = nBins_tmp;
  } else if( regionName == "HTge750_6j0b" ) {
    const int nBins_tmp                      = 5;
    bins = new double[nBins_tmp+1]{160, 200, 250, 325, 425, 600};
    nBins = nBins_tmp;
  } else if( regionName == "HTge750_6j1b" ) {
    const int nBins_tmp                      = 4;
    bins = new double[nBins_tmp+1]{150, 190, 250, 350, 500};
    nBins = nBins_tmp;
  } else if( regionName == "HTge750_6j2b" ) {
    const int nBins_tmp                      = 4;
    bins = new double[nBins_tmp+1]{130, 170, 220, 300, 450};
    nBins = nBins_tmp;
  } else if( regionName == "HTge750_3b" ) {
    const int nBins_tmp                        = 3;
    bins = new double[nBins_tmp+1]{125, 175, 275, 450};
    nBins = nBins_tmp;

  // HT > 450 && HT < 750
  } else if( regionName == "HTge450_2j0b" ) {
    const int nBins_tmp                      = 8;
    bins = new double[nBins_tmp+1]{200, 240, 290, 350, 420, 490, 570, 650, 750};
    nBins = nBins_tmp;
  } else if( regionName == "HTge450_2j1to2b" ) {
    const int nBins_tmp                      = 6;
    bins = new double[nBins_tmp+1]{200, 250, 310, 380, 450, 550, 700};
    nBins = nBins_tmp;
  } else if( regionName == "HTge450_3to5j0b" ) {
  const int nBins_tmp                      = 8;
  bins = new double[nBins_tmp+1]{200, 240, 290, 350, 420, 490, 570, 650, 750};
    nBins = nBins_tmp;
  } else if( regionName == "HTge450_3to5j1b" ) {
    const int nBins_tmp                      = 6;
    bins = new double[nBins_tmp+1]{200, 250, 310, 380, 460, 550, 700};
    nBins = nBins_tmp;
  } else if( regionName == "HTge450_3to5j2b" ) {
    const int nBins_tmp                      = 4;
    bins = new double[nBins_tmp+1]{200, 250, 325, 425, 550};
    nBins = nBins_tmp;
  } else if( regionName == "HTge450_6j0b" ) {
    const int nBins_tmp                      = 3;
    bins = new double[nBins_tmp+1]{200, 280, 380, 520};
    nBins = nBins_tmp;
  } else if( regionName == "HTge450_6j1b" ) {
    const int nBins_tmp                      = 3;
    bins = new double[nBins_tmp+1]{200, 250, 325, 450};
    nBins = nBins_tmp;
  } else if( regionName == "HTge450_6j2b" ) {
    const int nBins_tmp                      = 3;
    bins = new double[nBins_tmp+1]{200, 250, 300, 400};
    nBins = nBins_tmp;
  } else if( regionName == "HTge450_3b" ) {
    const int nBins_tmp                        = 2;
    bins = new double[nBins_tmp+1]{200, 280, 400};
    nBins = nBins_tmp;


  } else {

    std::cout << "[MT2Region::getBins] WARNING! Unknown regionName: " << regionName << std::endl;
    std::cout << "-> Aborting." << std::endl;
    exit(981);

  }


}




