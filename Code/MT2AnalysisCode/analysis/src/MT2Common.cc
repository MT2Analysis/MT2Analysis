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
      s.file = fileName;
      
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


std::vector<MT2SampleBaby> MT2Common::loadSamplesBaby(const std::string& dir, const std::string& outputdir) {

  std::vector<MT2SampleBaby> fSamples;

  char buffer[200];
  ifstream IN(dir.c_str());

  //char ParName[100];

  char StringValue[1000];
  float ParValue;

  std::cout << "------------------------------------" << std::endl;
  std::cout << "Sample Directory  " << dir << std::endl;
  
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
      std::cout << " Path  " << fPath << std::endl;

    }

    if( !strcmp(buffer, "SAMPLE")){

      MT2SampleBaby s;

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Name\t%s", StringValue);
      s.name = TString(StringValue);

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "SName\t%s", StringValue);
      s.sname = TString(StringValue);

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Dir\t%s", StringValue);
      TString dirName =fPath+StringValue;
      s.dir = dirName;

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "SampleID\t%f", &ParValue);
      s.id = ParValue;

      s.xsection = MT2Common::getXSection( s.id );
    

      ////// hadd of all files in dir -> Result in 'outputdir'
      system(Form("mkdir -p %s", outputdir.c_str()));

      system(Form("ls %s/*.root > listOfFiles_%s", dir.c_str(), s.name.c_str()));
      system(Form("hadd -f %s/%s.root @listOfFiles_%s", outputdir.c_str(), s.name.c_str(), s.name.c_str()));
      
      TString fileName = outputdir+"/"+s.name+".root";
      
      std::cout << "File " << fileName << " will be processed..." << std::endl;

      TFile* file = TFile::Open(fileName);
      TTree* tree = (TTree*)file->Get("MassTree");

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

      std::cout << "   ---- " << std::endl;
      std::cout << "   New sample added: " << s.name << std::endl;
      std::cout << "   Sample no.        " << counter << std::endl;
      std::cout << "   Short name:       " << s.sname << std::endl;
      std::cout << "   Directory:        " << s.dir << std::endl;
      std::cout << "   Events:           " << s.nevents  << std::endl;
      std::cout << "   Events in tree:   " << tree->GetEntries() << std::endl;
      std::cout << "   ID:               " << s.id << std::endl;
      std::cout << "   Xsection:         " << s.xsection << std::endl;
      std::cout << "   avg PU weight:    " << s.PU_avg_weight << std::endl;

      fSamples.push_back(s);
      file->Close();
      tree = 0;
      counter++;

    }

  }

  std::cout << "------------------------------------" << std::endl;

  return fSamples;

}


//////
std::vector<MT2SamplePostBaby> MT2Common::loadSamplesPostBaby(const std::string& filename) {

  std::vector<MT2SamplePostBaby> fSamples;

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

      MT2SamplePostBaby s;
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Name\t%s", StringValue);
      s.name = TString(StringValue);

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "SName\t%s", StringValue);
      s.sname = TString(StringValue);

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "File\t%s", StringValue);
      TString fileName =fPath+StringValue;
      s.file = fileName;



      TFile* file = TFile::Open(s.file.c_str());
      TTree* tree = (TTree*)file->Get("treeProducerSusyFullHad");

      int evt_id;
      float evt_nEvts, evt_filter, evt_kfactor, evt_xsec, evt_scale1fb;

      tree->SetBranchAddress("evt_id",&evt_id);
      tree->SetBranchAddress("evt_nEvts",&evt_nEvts);
      tree->SetBranchAddress("evt_filter",&evt_filter);
      tree->SetBranchAddress("evt_kfactor",&evt_kfactor);
      tree->SetBranchAddress("evt_xsec",&evt_xsec);
      tree->SetBranchAddress("evt_scale1fb",&evt_scale1fb);
      
      tree->GetEntry(0);

      s.id       = evt_id;
      s.nevents  = evt_nEvts;
      s.xsection = evt_xsec;
      s.filter   = evt_filter;
      s.kfact    = evt_kfactor;
      s.scale1fb = evt_scale1fb;

      // Manually set variables
      s.PU_avg_weight = 1;
      s.lumi = 1000; // 1/fb
      
      std::cout << " ---- " << std::endl;
      std::cout << "  New sample added: " << s.name << std::endl;
      std::cout << "   Sample no.      " << counter << std::endl;
      std::cout << "   Short name:     " << s.sname << std::endl;
      std::cout << "   File:           " << s.file << std::endl;
      std::cout << "   ID:             " << s.id << std::endl;
      std::cout << "   Events:         " << s.nevents  << std::endl;
      std::cout << "   Xsection:       " << s.xsection << std::endl;
      std::cout << "   Lumi:           " << s.lumi << std::endl;
      std::cout << "   Kfactor:        " << s.kfact << std::endl;
      std::cout << "   Filter:         " << s.filter << std::endl;
      std::cout << "   Scale to 1/fb:  " << s.scale1fb << std::endl;
      std::cout << "   Avg PU weight:  " << s.PU_avg_weight << std::endl;
      //std::cout << "   type:           " << s.type << std::endl;
      //std::cout << "   Color:          " << s.color << std::endl;
      fSamples.push_back(s);
      file->Close();
      tree = 0;
      counter++;
    
    }

  }

  std::cout << "------------------------------------" << std::endl;

  return fSamples;

}
//////
//////
std::vector<MT2SampleBaby_basic> MT2Common::loadSamplesBaby_basic(const std::string& filename) {

  std::vector<MT2SampleBaby_basic> fSamples;

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

      MT2SampleBaby_basic s;
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Name\t%s", StringValue);
      s.name = TString(StringValue);

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "SName\t%s", StringValue);
      s.sname = TString(StringValue);

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "File\t%s", StringValue);
      TString fileName =fPath+StringValue;
      s.file = fileName;

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Xsection\t%f", &ParValue);
      s.xsection = ParValue;

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Kfact\t%f", &ParValue);
      s.kfact = ParValue;

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Lumi\t%f", &ParValue);
      s.lumi = ParValue;


      TFile* file = TFile::Open(s.file.c_str());
      TTree* tree = (TTree*)file->Get("treeProducerSusyFullHad");

      s.PU_avg_weight = 1;
      s.nevents       = tree->GetEntries();


      std::cout << " ---- " << std::endl;
      std::cout << "  New sample added: " << s.name << std::endl;
      std::cout << "   Sample no.      " << counter << std::endl;
      std::cout << "   Short name:     " << s.sname << std::endl;
      std::cout << "   File:           " << s.file << std::endl;
      std::cout << "   Events:         " << s.nevents  << std::endl;
      //std::cout << "   Events in tree: " << tree->GetEntries() << std::endl;
      std::cout << "   Xsection:       " << s.xsection << std::endl;
      std::cout << "   Lumi:           " << s.lumi << std::endl;
      std::cout << "   kfactor:        " << s.kfact << std::endl;
      std::cout << "   avg PU weight:  " << s.PU_avg_weight << std::endl;
      //std::cout << "   type:           " << s.type << std::endl;
      //std::cout << "   Color:          " << s.color << std::endl;
      fSamples.push_back(s);
      file->Close();
      tree = 0;
      counter++;
    }

  }

  std::cout << "------------------------------------" << std::endl;

  return fSamples;

}     
//////

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
    else if(sample.sname=="Top")    sampletype = "SingleTop"; // no ttbar, includes TTZ, TTW
    else sampletype = "Other";
  }

  return sampletype;

}



float MT2Common::getXSection( int sampleId ) {

  float xsection=0.; //[pb]


  ////if      ( sampleId >= 0   && sampleId < 100 )   xsection = 1.0; // data
  
  ////else if ( sampleId >= 100 && sampleId < 200 )   xsecti0n = 1.1; // QCD
  //else if ( sampleId == 100 )   xsection =           1.0000; // QCD_Pt-15to8000 Flat
  //else if ( sampleId == 101 )   xsection =           1.0000; // QCD_Pt-0to5
  /*else*/ if ( sampleId == 102 )   xsection = 80710000000.0000; // QCD_Pt-5to10_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1[v2]
  else if ( sampleId == 103 )   xsection =  7528000000.0000; // QCD_Pt-10to15_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 104 )   xsection =  2237000000.0000; // QCD_Pt-15to30_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 105 )   xsection =   161500000.0000; // QCD_Pt-30to50_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v2 
  else if ( sampleId == 106 )   xsection =    22110000.0000; // QCD_Pt-50to80_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 107 )   xsection =     3000114.3000; // QCD_Pt-80to120_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 108 )   xsection =      493200.0000; // QCD_Pt-120to170_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 109 )   xsection =      120300.0000; // QCD_Pt-170to300_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 110 )   xsection =        7475.0000; // QCD_Pt-300to470_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 111 )   xsection =         587.1000; // QCD_Pt-470to600_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 112 )   xsection =         167.0000; // QCD_Pt-600to800_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 113 )   xsection =          28.2500; // QCD_Pt-800to1000_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 114 )   xsection =           8.1950; // QCD_Pt-1000to1400_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 115 )   xsection =           0.7346; // QCD_Pt-1400to1800_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 116 )   xsection =           0.1091; // QCD_Pt-1800_Tune4C_13TeV_pythia8/Spring14dr-castor_PU_S14_POSTLS170_V6-v1

  ////else if ( sampleId >= 200 && sampleId < 300 )   xsection = 1.2;  // GJets
  //else if ( sampleId == 201 )   xsection =           1.2000; // GJets_Pt15to8000 Flat
  else if ( sampleId == 201 )   xsection =        1534.0000; // GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 202 )   xsection =         489.9000; // GJets_HT-200to400_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 203 )   xsection =          62.0500; // GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU_S14_POSTLS170_V6-v1
  //else if ( sampleId == 203 )   xsection =           1.2000; // GJets_HT-600toInf

  ////else if ( sampleId >= 300 && sampleId < 400 )   xsection = 1.3; // TTbar
  else if ( sampleId == 300 )   xsection =         424.5000; // TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU_S14_POSTLS170_V6-v1 [inclusive]
  //else if ( sampleId == 301 )   xsection =           1.3000; // TTJets-HadronicMGDecays
  //else if ( sampleId == 302 )   xsection =           1.3000; // TTJets-SemiLeptMGDecays
  //else if ( sampleId == 303 )   xsection =           1.3000; // TTJets-FullLeptMGDecays

  ////else if ( sampleId >= 400 && sampleId < 500 )   xsection = 1.4; // Top
  //else if ( sampleId == 400 )   xsection =           1.4000; // T-s-channel
  //else if ( sampleId == 401 )   xsection =           1.4000; // T-t-channel
  //else if ( sampleId == 402 )   xsection =           1.4000; // T-tW-channel
  //else if ( sampleId == 403 )   xsection =           1.4000; // Tbar-s-channel
  //else if ( sampleId == 404 )   xsection =           1.4000; // Tbar-t-channel
  //else if ( sampleId == 405 )   xsection =           1.4000; // Tbar-tW-channel

  ////else if ( sampleId >= 500 && sampleId < 600 )   xsection = 1.5; // WJets
  //else if ( sampleId == 500 )   xsection =           1.5000; // WJetsToLNu inclusive
  else if ( sampleId == 501 )   xsection =        1817.0000; // WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 502 )   xsection =         471.6000; // WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 503 )   xsection =          55.6100; // WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 504 )   xsection =          18.8100; // WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU_S14_POSTLS170_V6-v1

  ////else if ( sampleId >= 600 && sampleId < 700 )   xsection = 1.6; // ZJets & DY
  else if ( sampleId == 600 )   xsection =          52.2400; // DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 601 )   xsection =           6.5460; // DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU_S14_POSTLS170_V6-v1
  else if ( sampleId == 602 )   xsection =           2.1790; // DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU_S14_POSTLS170_V6-v1
  //else if ( sampleId == 603 )   xsection =           1.6000; // ZJetsToNuNu-HT-50to100 
  //else if ( sampleId == 604 )   xsection =           1.6000; // ZJetsToNuNu-HT-100to200
  else if ( sampleId == 605 )   xsection =         100.8000; // ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU_S14_POSTLS170_V6-v1
  //else if ( sampleId == 606 )   xsection =           1.6000; // ZJetsToNuNu-HT-400toInf
  

  ////else if ( sampleId >= 700 && sampleId < 800 )   xsection = 1.7; // VV 
  
  ////else if ( sampleId >= 1000                  )   xsection = 1.8; // SUSY

  else {
    
    std::cout << "ERROR: Id " << sampleId << " does not correspond to any sample! " << std::endl;
    exit(1);
  
  }

  return xsection;

}
