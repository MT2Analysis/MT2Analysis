#include "../interface/MT2Common.h"

#include <iostream>
#include <cstdlib>
#include <fstream>



MT2Common::MT2Common() {};



MT2Common::~MT2Common() {};


void MT2Common::getBins( float ht_min, int njet_min, int njet_max, int nbjet_min, int nbjet_max, int &nBins, Double_t* bins) {

  std::string signal_region = MT2Common::getSignalRegion( ht_min, njet_min, njet_max, nbjet_min, nbjet_max );

  return MT2Common::getBins( signal_region, nBins, bins );

}




void MT2Common::getBins( const std::string& signal_region, int &nBins, Double_t* bins) {

  if( signal_region == "HTge1200_2j0b" ) {
    const int nBins_tmp                      = 6;
    bins = new Double_t[nBins_tmp+1]{120, 150, 200, 260, 350, 550, 900};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge1200_2j1b" ) {
    const int nBins_tmp                      = 2;
    bins = new Double_t[nBins_tmp+1]{100, 180, 350};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge1200_3j0b" ) {
    const int nBins_tmp                      = 7;
    bins = new Double_t[nBins_tmp+1]{160, 185, 220, 270, 350, 450, 650, 1000};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge1200_3j1b" ) {
    const int nBins_tmp                      = 4;
    bins = new Double_t[nBins_tmp+1]{150, 180, 230, 350, 550};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge1200_3j2b" ) {
    const int nBins_tmp                      = 2;
    bins = new Double_t[nBins_tmp+1]{130, 200, 350};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge1200_6j0b" ) {
    const int nBins_tmp                      = 3;
    bins = new Double_t[nBins_tmp+1]{160, 200, 300, 500};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge1200_6j1b" ) {
    const int nBins_tmp                      = 3;
    bins = new Double_t[nBins_tmp+1]{150, 200, 300, 500};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge1200_6j2b" ) {
    const int nBins_tmp                      = 2;
    bins = new Double_t[nBins_tmp+1]{130, 200, 350};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge1200_3b" ) {
    const int nBins_tmp                      = 1;
    bins = new Double_t[nBins_tmp+1]{125, 300};
    nBins = nBins_tmp;

  // HT > 750 && HT < 1200
  } else if( signal_region == "HTge750_2j0b" ) {
    const int nBins_tmp                      = 9;
    bins = new Double_t[nBins_tmp+1]{125, 150, 180, 220, 270, 325, 425, 580, 780, 1000};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge750_2j1b" ) {
    const int nBins_tmp                      = 5;
    bins = new Double_t[nBins_tmp+1]{100, 135, 170, 260, 450, 700};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge750_3j0b" ) {
    const int nBins_tmp                      = 9;
    bins = new Double_t[nBins_tmp+1]{160, 185, 215, 250, 300, 370, 480, 640, 800, 1000};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge750_3j1b" ) {
    const int nBins_tmp                      = 6;
    bins = new Double_t[nBins_tmp+1]{150, 175, 210, 270, 380, 600, 900};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge750_3j2b" ) {
    const int nBins_tmp                      = 5;
    bins = new Double_t[nBins_tmp+1]{130, 160, 200, 270, 370, 500};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge750_6j0b" ) {
    const int nBins_tmp                      = 5;
    bins = new Double_t[nBins_tmp+1]{160, 200, 250, 325, 425, 600};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge750_6j1b" ) {
    const int nBins_tmp                      = 4;
    bins = new Double_t[nBins_tmp+1]{150, 190, 250, 350, 500};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge750_6j2b" ) {
    const int nBins_tmp                      = 4;
    bins = new Double_t[nBins_tmp+1]{130, 170, 220, 300, 450};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge750_3b" ) {
    const int nBins_tmp                        = 3;
    bins = new Double_t[nBins_tmp+1]{125, 175, 275, 450};
    nBins = nBins_tmp;

  // HT > 450 && HT < 750
  } else if( signal_region == "HTge450_2j0b" ) {
    const int nBins_tmp                      = 8;
    bins = new Double_t[nBins_tmp+1]{200, 240, 290, 350, 420, 490, 570, 650, 750};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge450_2j1b" ) {
    const int nBins_tmp                      = 6;
    bins = new Double_t[nBins_tmp+1]{200, 250, 310, 380, 450, 550, 700};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge450_3j0b" ) {
    const int nBins_tmp                      = 8;
    bins = new Double_t[nBins_tmp+1]{200, 240, 290, 350, 420, 490, 570, 650, 750};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge450_3j1b" ) {
    const int nBins_tmp                      = 6;
    bins = new Double_t[nBins_tmp+1]{200, 250, 310, 380, 460, 550, 700};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge450_3j2b" ) {
    const int nBins_tmp                      = 4;
    bins = new Double_t[nBins_tmp+1]{200, 250, 325, 425, 550};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge450_6j0b" ) {
    const int nBins_tmp                      = 3;
    bins = new Double_t[nBins_tmp+1]{200, 280, 380, 520};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge450_6j1b" ) {
    const int nBins_tmp                      = 3;
    bins = new Double_t[nBins_tmp+1]{200, 250, 325, 450};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge450_6j2b" ) {
    const int nBins_tmp                      = 3;
    bins = new Double_t[nBins_tmp+1]{200, 250, 300, 400};
    nBins = nBins_tmp;
  } else if( signal_region == "HTge450_3b" ) {
    const int nBins_tmp                        = 2;
    bins = new Double_t[nBins_tmp+1]{200, 280, 400};
    nBins = nBins_tmp;


  } else {

    std::cout << "WARNING! Unknown signal_region: " << signal_region << std::endl;
    std::cout << "-> Aborting." << std::endl;
    exit(981);

  }

}



std::string MT2Common::getSignalRegion( float ht_min, int njet_min, int njet_max, int nbjet_min, int nbjet_max) {

  char htString_ch[64];
  sprintf( htString_ch, "HTge%.0f", ht_min );
  std::string htString(htString_ch);

  std::string jString = MT2Common::getSingleSignalRegionString( njet_min , njet_max  , "j");
  std::string bString = MT2Common::getSingleSignalRegionString( nbjet_min, nbjet_max , "b");

  std::string signal_region = htString + "_" + jString + bString;

  return signal_region;

}


std::string MT2Common::getSingleSignalRegionString( int n_min , int n_max, const std::string& suffix ) {

  if( n_max==999 ) n_max=n_min;

  char signal_region_ch[64];
  if( n_max!=n_min )
    sprintf( signal_region_ch, "%dto%d%s", n_min, n_max, suffix.c_str() );
  else 
    sprintf( signal_region_ch, "%d%s", n_min, suffix.c_str() );
  
  std::string signal_region(signal_region_ch);

  return signal_region;

}




std::vector<MT2Sample> MT2Common::loadSamples(const std::string& filename) {

  std::vector<MT2Sample> fSamples;

  char buffer[200];
  ifstream IN(filename.c_str());
  
  //char ParName[100];
  char StringValue[1000];
  float ParValue;
  
  std::cout << "------------------------------------" << std::endl;
  std::cout << "Sample File  " << filename << std::endl;

  int counter(0);
  TString fPath; 

  
  while( IN.getline(buffer, 200, '\n') ) {
    
    // ok = false;
    if (buffer[0] == '#') {
      continue; // Skip lines commented with '#'
    }

    if( !strcmp(buffer, "GENERAL") ) {

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Path\t%s", StringValue);
      fPath = StringValue;
      
      std::cout << " ----  " << std::endl;
      std::cout << "  Path " << fPath << std::endl;
    
    }

    if( !strcmp(buffer, "SAMPLE")){
    
      MT2Sample s;
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Name\t%s", StringValue);
      s.name = TString(StringValue);
      
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "SName\t%s", StringValue);
      s.sname = TString(StringValue);
      
      IN.getline(buffer, 200, '\n');			//new
      sscanf(buffer, "ShapeName\t%s", StringValue);	//new
      s.shapename = TString(StringValue);		//new
      
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "File\t%s", StringValue);
      TString fileName =fPath+StringValue;
      //TFile *f = TFile::Open(file);
      s.file = fileName;
      //s.tree = (TTree*)f->Get("MassTree");
      
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Xsection\t%f", &ParValue);
      s.xsection = ParValue;
      
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Kfact\t%f", &ParValue);
      s.kfact = ParValue;
      
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Lumi\t%f", &ParValue);
      s.lumi = ParValue;
      
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Type\t%s", StringValue);
      s.type = StringValue;
      
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Color\t%f", &ParValue);
      s.color = ParValue;
      

      TFile* file = TFile::Open(s.file.c_str());
      TTree* tree = (TTree*)file->Get("MassTree");

      if(s.type!="data"){

        TH1F *h_PUWeights = (TH1F*) file->Get("h_PUWeights");
        TH1F *h_Events    = (TH1F*) file->Get("h_Events");
        if(h_PUWeights==0 || h_Events==0){
          std::cout << "ERROR: sample " << (file)->GetName() << " does not have PU and NEvents histos! " << std::endl;
          exit(1);
        }
        s.PU_avg_weight = h_PUWeights->GetMean();
        s.nevents       = h_Events   ->GetEntries();
        delete h_PUWeights;
        delete h_Events;

      } else {

        s.PU_avg_weight =1;
        s.nevents       =1;

      }

      std::cout << " ---- " << std::endl;
      std::cout << "  New sample added: " << s.name << std::endl;
      std::cout << "   Sample no.      " << counter << std::endl;
      std::cout << "   Short name:     " << s.sname << std::endl;
      std::cout << "   File:           " << s.file << std::endl;
      std::cout << "   Events:         " << s.nevents  << std::endl;
      std::cout << "   Events in tree: " << tree->GetEntries() << std::endl; 
      std::cout << "   Xsection:       " << s.xsection << std::endl;
      std::cout << "   Lumi:           " << s.lumi << std::endl;
      std::cout << "   kfactor:        " << s.kfact << std::endl;
      std::cout << "   avg PU weight:  " << s.PU_avg_weight << std::endl;
      std::cout << "   type:           " << s.type << std::endl;
      std::cout << "   Color:          " << s.color << std::endl;
      fSamples.push_back(s);
      file->Close();
      tree = 0;
      counter++;
    }

  }

  std::cout << "------------------------------------" << std::endl;

  return fSamples;

}



std::string MT2Common::getSampleType( MT2Sample sample ) {

  std::string sampletype = sample.type;
  if(sampletype=="mc"){
    if(sample.sname=="QCD")         sampletype = "QCD";
    else if(sample.sname=="Wtolnu") sampletype = "WJets";
    else if(sample.sname=="DY")     sampletype = "ZJets";
    else if(sample.name=="TTbar")   sampletype = "TTbar";
    else if(sample.name=="TTbar_Madgraph0l")   sampletype = "TTbar";
    else if(sample.name=="TTbar_Madgraph1l")   sampletype = "TTbar";
    else if(sample.name=="TTbar_Madgraph2l")   sampletype = "TTbar";
    else if(sample.sname=="Top")    sampletype = "SingleTop";//no ttbar, includes TTZ, TTW
    else sampletype = "Other";
  }

  return sampletype;

}
