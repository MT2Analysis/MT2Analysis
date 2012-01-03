// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>

#include "MT2Analyzer.hh"

using namespace std;

//________________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	cout << "Usage: RunMT2Analyzer [-d dir] [-o filename] [-v verbose] [-j json]                  " << endl;
	cout << "                      [-m set_of_cuts] [-n maxEvents] [-t type]                      " << endl;
	cout << "                      [-p data_PileUp] [-P mc_PileUP]                                " << endl; 
        cout << "                      [-s S3,noPU,3D] [-C JEC]                                          " << endl;
	cout << "                      [-r photon ] [-i ID ]                                          " << endl;
	cout << "                      [-l] file1 [... filen]"                                          << endl;
	cout << "  where:"                                                                              << endl;
	cout << "     dir           is the output directory                                           " << endl;
	cout << "                   default is TempOutput/                                            " << endl;
	cout << "     filename      is the output filename for the MassAnalysis                       " << endl;
	cout << "     verbose       sets the verbose level                                            " << endl;
	cout << "                   default is 0 (quiet mode)                                         " << endl;
	cout << "     json          json file to be read                                              " << endl;
	cout << "     set_of_cuts   optional cuts for MT2Analysis                                     " << endl;
	cout << "     data_PileUp   root file from which the expected # pile-up                       " << endl;
	cout << "                   interactions is read                                              " << endl;
	cout << "     mc_PileUP     root file from which the generated # pile up                      " << endl;
	cout << "                   interactions is read                                              " << endl;
	cout << "     type          data or mc=default                                                " << endl;
	cout << "     photon        add Photon to MET and remove jets/ele matched to photon           " << endl;
	cout << "     ID            ProcessID: 0=data, 1=Znunu, 2=Zll, 3=WJets, 4=Top,                " << endl;
	cout << "                              5=Gamma+jets, 6=QCD ,[empty], 9=Other, 10=Signal       " << endl;
	cout << "     JEC           redo JEC: dir of JEC to be used                                   " << endl;
	cout << "                   /shome/pnef/MT2Analysis/Code/JetEnergyCorrection/[JEC]            " << endl;
	cout << "                   ak5 pf-jets will be corrected with L1FastL2L3 (+RES if type=data) " << endl;
	cout << "     filen         are the input files (by default: ROOT files)                      " << endl;
	cout << "                   with option -l, these are read as text files                      " << endl;
	cout << "                   with one ROOT file name per line                                  " << endl;
	cout << endl;
	exit(status);
}

//________________________________________________________________________________________
int main(int argc, char* argv[]) {
// Default options
	bool isList = false;
	TString outputdir = "TempOutput/";
	TString filename  = "MassTree.root";
	TString setofcuts = "default";
	string  puScenario = "";
  	string  jsonFileName = "";
	string  data_PileUp = "";
	string  mc_PileUp = "";
	string  type = "mc";
	string  JEC  ="";
	bool isData  = false;
	int verbose  = 0;
	int maxEvents=-1;
	int ID       =-1;
	bool isS3    = false;
	bool noPU    = false;
    	bool removePhoton = false;
	string photon = "";

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "s:d:o:v:j:m:n:p:P:t:r:i:C:lh?")) != -1 ) {
	  switch (ch) {
	  case 'd': outputdir       = TString(optarg);break;
	  case 'o': filename        = TString(optarg);break;
	  case 'v': verbose         = atoi(optarg);   break;
	  case 'j': jsonFileName    = string(optarg); break;
	  case 'm': setofcuts       = TString(optarg);break;
	  case 'n': maxEvents       = atoi(optarg);   break;
	  case 'p': data_PileUp     = string(optarg); break;
	  case 'P': mc_PileUp       = string(optarg); break;
	  case 't': type            = string(optarg); break;
	  case 'r': photon          = string(optarg); break;
	  case 'i': ID              = atoi(optarg);   break;
	  case 's': puScenario      = string(optarg); break;
          case 'C': JEC             = "/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/"+string(optarg)+"/"; break;
	  case 'l': isList          = true; break;
	    //case 'noPU': noPU = true; break;  
	  case '?':
	  case 'h': usage(0); break;
	  default:
	    cerr << "*** Error: unknown option " << optarg << std::endl;
	    usage(-1);
	  }
	}

	argc -= optind;
	argv += optind;

// Check arguments
	if( argc<1 ) {
		usage(-1);
	}
	if      (type=="data") isData =true;
	else if (type=="mc"  ) isData =false;
	else    usage(-1);
	if      (photon == "photon") removePhoton=true;
	if      (!isData && data_PileUp.length()==0  ) {
		cout << "                        WARNING: need data_PileUp to run on MC " << endl;
	}
	if      ( isData && (data_PileUp.length() >0 || mc_PileUp.length() >0)  ) {
		cout << "ERROR: you are running on data, no reweighting needed... " << endl; exit(-1);
	}

	setofcuts   ="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/MT2_cuts/"+setofcuts+".dat";
	if(data_PileUp.length()!=0){data_PileUp ="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/Certification/pileUp_data/"+data_PileUp;}
	if(mc_PileUp.length()  !=0){mc_PileUp   ="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/Certification/pileUp_mc/"  + mc_PileUp;}
	//setofcuts   ="/shome/leo/Analysis/MT2_cuts/"+setofcuts+".dat";
	//if(data_PileUp.length()!=0){data_PileUp ="/shome/leo/Analysis/Certification/pileUp_data/"+data_PileUp;}
        //if(mc_PileUp.length()  !=0){mc_PileUp   ="/shome/leo/Analysis/Certification/pileUp_mc/"  + mc_PileUp;}

	if(jsonFileName.length() !=0){jsonFileName="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/Certification/"           +jsonFileName;}

	if(puScenario=="3D"){ isS3=true; noPU=true; } // THIS IS A DIRTY TRICK TO TEST 3D REWEIGHT WITHOUT ADD A NEW VAR
	if(puScenario=="S3") isS3=true;
	else if(puScenario=="noPU") noPU=true;

	TChain *theChain = new TChain("analyze/Analysis");
	for(int i = 0; i < argc; i++){
		if( !isList ){
			theChain->Add(argv[i]);
			printf(" Adding file: %s\n",argv[i]);
		} else {
			TString rootFile;
			ifstream is(argv[i]);
			while(rootFile.ReadLine(is) && (!rootFile.IsNull())){
				if(rootFile[0] == '#') continue;
				theChain->Add(rootFile);
				printf(" Adding file: %s\n", rootFile.Data());
			}
		}
	}


	cout << "--------------" << endl;
	cout << "OutputDir is:                   " << outputdir << endl;
	cout << "Type is:                        " << type << endl;
	cout << "Verbose level is:               " << verbose << endl;
  	cout << "JSON file is:                   " << (jsonFileName.length()>0?jsonFileName:"empty") << endl;
  	cout << "MC_PileUp file:                 " << (mc_PileUp.length()>0?mc_PileUp:"empty") << endl;
  	cout << "Data_PileUp file:               " << (data_PileUp.length()>0?data_PileUp:"empty") << endl;
	if(noPU && !isS3) cout << "WARNING: NoPU option set, all the PU weights will be set to 1" << endl;
	cout << "Set of Cuts is:                 " << setofcuts << endl;
	cout << "Number of events:               " << theChain->GetEntries() << endl;
	if(JEC.length()!=0){
	cout << "Redo JEC with set               " << JEC << endl;
	}
	if(removePhoton){
	cout << "WARNING: Photon is added to MET and jet/ele match to photon is removed!!" << endl;
	}
	cout << "--------------" << endl;

	MT2Analyzer *tA = new MT2Analyzer(theChain);
	tA->SetOutputDir(outputdir);
	tA->SetVerbose(verbose);
	tA->SetMaxEvents(maxEvents);
	tA->SetProcessID(ID);
	tA->isS3 = isS3;
	tA->noPU = noPU;
	tA->removePhoton = removePhoton;
  	if (jsonFileName!="") tA->ReadJSON(jsonFileName.c_str());
	tA->BeginJob(filename, setofcuts, isData, data_PileUp, mc_PileUp, JEC);
	tA->Loop();
	tA->EndJob();

	delete tA;
	return 0;
}

