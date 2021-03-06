#include "helper/Utilities.hh"
#include "MT2Analysis.hh"
#include "TLorentzVector.h"
#include <sstream>


using namespace std;

MT2Analysis::MT2Analysis(TreeReader *tr) : UserAnalysisBase(tr){
	Util::SetStyle();	
	fCut_PFMET_min                      = 0;
	fCut_HT_min                         = 0;
	fCut_caloHT50_min                   = 0;
	fCut_caloHT50ID_min                 = 0;
	fCut_caloMHT30_min                  = 0;
	fCut_caloMHT30ID_min                = 0;
	fCut_JPt_hardest_min                = 0;
	fCut_JPt_second_min                 = 0;
        fCut_PtHat_max                      = 999999.;
	fCut_Run_min                        = 0;
	fCut_Run_max                        = 9999999;
	fDoJESUncertainty                   = false;
	fJESUpDown                          = 0;

	fRemovePhoton                       = 0;
	fID                                 = -1;
	fbtagFileName                       = "";

	fRequiredHLT.clear();
	fVetoedHLT.clear();
	
	fBEfiles.clear();
	fTPfiles.clear();
}

MT2Analysis::~MT2Analysis(){
	delete fJecUncPF;    
	delete fJecUncCalo;  
	
	delete fJecCaloL2;
	delete fJecCaloL3; 
	delete fJecPFL1;   
	delete fJecPFL2;   
	delete fJecPFL3;   
	delete fJecPFRES;  

	delete fJetCorrectorPF;
	delete fJetCorrectorCalo;

}

void MT2Analysis::End(){
	cout << " *************************************************************** " << endl;
	cout << " MT2Analysis::End()                                             " << endl;
	cout << " *************************************************************** " << endl;
	
	fHistFile->cd();	

	// write tree
	fH_PUWeights->Write();
	fH_Events->Write();

	if(isScan){
	  fH2_mSugraEvents->Write();
	  fH2_SMSEvents->Write();
	  for(int s=1;s<11; s++){
	    fH_mSugraSubProcEvents[s]->Write();
	  }
	}
	fATree->Write();
	fHistFile                ->Close();


	cout << " MT2Analysis::RealEnd()                                             " << endl;
	cout << " *************************************************************** " << endl;
}

void MT2Analysis::Begin(const char* filename){
	// Define the output file of histograms
	fHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");
	TDirectory *dir = gDirectory;

	//define btagging files
	bool existing=true;
	std::ifstream ifile(fbtagFileName.c_str() );
	if(!(ifile)) existing = false;
	if(existing){
		btagfile = TFile::Open(fbtagFileName.c_str() );
		hbeff = (TH1D*)btagfile->Get("h_beff"); hbeff->SetDirectory(dir);
		hceff = (TH1D*)btagfile->Get("h_ceff"); hceff->SetDirectory(dir);
		hleff = (TH1D*)btagfile->Get("h_leff"); hleff->SetDirectory(dir);
		btagfile->Close();
		dir->cd();
	}else{
		cout << "No btagfile existing: use b-eff = 0.5, c-eff = 0.08, l-eff = 0.001" << endl;
	}


	// book tree
	fMT2tree = new MT2tree();
	BookTree();
	
	// define which triggers to fill
	if(fisData){
		// HT with dPhi
		fTriggerMap["HLT_HT500_JetPt60_DPhi2p94_v1"] = &fMT2tree->trigger.HLT_HT500_JetPt60_DPhi2p94_v1;
		fTriggerMap["HLT_HT550_JetPt60_DPhi2p94_v1"] = &fMT2tree->trigger.HLT_HT550_JetPt60_DPhi2p94_v1;

		// HT
		fTriggerMap["HLT_HT150_v2"]                 = &fMT2tree->trigger.HLT_HT150_v2;
		fTriggerMap["HLT_HT150_v3"]                 = &fMT2tree->trigger.HLT_HT150_v3;
		fTriggerMap["HLT_HT160_v2"]                 = &fMT2tree->trigger.HLT_HT160_v2;
		fTriggerMap["HLT_HT200_v2"]                 = &fMT2tree->trigger.HLT_HT200_v2;
		fTriggerMap["HLT_HT200_v3"]                 = &fMT2tree->trigger.HLT_HT200_v3;
		fTriggerMap["HLT_HT240_v2"]                 = &fMT2tree->trigger.HLT_HT240_v2;
		fTriggerMap["HLT_HT250_v2"]                 = &fMT2tree->trigger.HLT_HT250_v2;
		fTriggerMap["HLT_HT250_v3"]                 = &fMT2tree->trigger.HLT_HT250_v3;
		fTriggerMap["HLT_HT260_v2"]                 = &fMT2tree->trigger.HLT_HT260_v2;
		fTriggerMap["HLT_HT300_v2"]                 = &fMT2tree->trigger.HLT_HT300_v2;
		fTriggerMap["HLT_HT300_v3"]                 = &fMT2tree->trigger.HLT_HT300_v3;
		fTriggerMap["HLT_HT300_v4"]                 = &fMT2tree->trigger.HLT_HT300_v4;
		fTriggerMap["HLT_HT300_v5"]                 = &fMT2tree->trigger.HLT_HT300_v5;
		fTriggerMap["HLT_HT350_v2"]                 = &fMT2tree->trigger.HLT_HT350_v2;
		fTriggerMap["HLT_HT350_v3"]                 = &fMT2tree->trigger.HLT_HT350_v3;
		fTriggerMap["HLT_HT350_v4"]                 = &fMT2tree->trigger.HLT_HT350_v4;
		fTriggerMap["HLT_HT360_v2"]                 = &fMT2tree->trigger.HLT_HT360_v2;
		fTriggerMap["HLT_HT400_v10"]                = &fMT2tree->trigger.HLT_HT400_v10;
		fTriggerMap["HLT_HT400_v2"]                 = &fMT2tree->trigger.HLT_HT400_v2;
		fTriggerMap["HLT_HT400_v3"]                 = &fMT2tree->trigger.HLT_HT400_v3;
		fTriggerMap["HLT_HT400_v4"]                 = &fMT2tree->trigger.HLT_HT400_v4;
		fTriggerMap["HLT_HT400_v5"]                 = &fMT2tree->trigger.HLT_HT400_v5;
		fTriggerMap["HLT_HT400_v6"]                 = &fMT2tree->trigger.HLT_HT400_v6;
		fTriggerMap["HLT_HT400_v7"]                 = &fMT2tree->trigger.HLT_HT400_v7;
		fTriggerMap["HLT_HT400_v8"]                 = &fMT2tree->trigger.HLT_HT400_v8;
		fTriggerMap["HLT_HT400_v9"]                 = &fMT2tree->trigger.HLT_HT400_v9;
		fTriggerMap["HLT_HT440_v2"]                 = &fMT2tree->trigger.HLT_HT440_v2;
		fTriggerMap["HLT_HT450_v10"]                = &fMT2tree->trigger.HLT_HT450_v10;
		fTriggerMap["HLT_HT450_v2"]                 = &fMT2tree->trigger.HLT_HT450_v2;
		fTriggerMap["HLT_HT450_v3"]                 = &fMT2tree->trigger.HLT_HT450_v3;
		fTriggerMap["HLT_HT450_v4"]                 = &fMT2tree->trigger.HLT_HT450_v4;
		fTriggerMap["HLT_HT450_v5"]                 = &fMT2tree->trigger.HLT_HT450_v5;
		fTriggerMap["HLT_HT450_v6"]                 = &fMT2tree->trigger.HLT_HT450_v6;
		fTriggerMap["HLT_HT450_v7"]                 = &fMT2tree->trigger.HLT_HT450_v7;
		fTriggerMap["HLT_HT450_v8"]                 = &fMT2tree->trigger.HLT_HT450_v8;
		fTriggerMap["HLT_HT450_v9"]                 = &fMT2tree->trigger.HLT_HT450_v9;
		fTriggerMap["HLT_HT500_v10"]                = &fMT2tree->trigger.HLT_HT500_v10;
		fTriggerMap["HLT_HT500_v11"]                = &fMT2tree->trigger.HLT_HT500_v11;
		fTriggerMap["HLT_HT500_v2"]                 = &fMT2tree->trigger.HLT_HT500_v2;
		fTriggerMap["HLT_HT500_v3"]                 = &fMT2tree->trigger.HLT_HT500_v3;
		fTriggerMap["HLT_HT500_v4"]                 = &fMT2tree->trigger.HLT_HT500_v4;
		fTriggerMap["HLT_HT500_v5"]                 = &fMT2tree->trigger.HLT_HT500_v5;
		fTriggerMap["HLT_HT500_v6"]                 = &fMT2tree->trigger.HLT_HT500_v6;
		fTriggerMap["HLT_HT500_v7"]                 = &fMT2tree->trigger.HLT_HT500_v7;
		fTriggerMap["HLT_HT500_v8"]                 = &fMT2tree->trigger.HLT_HT500_v8;
		fTriggerMap["HLT_HT500_v9"]                 = &fMT2tree->trigger.HLT_HT500_v9;
		fTriggerMap["HLT_HT550_v10"]                = &fMT2tree->trigger.HLT_HT550_v10;
		fTriggerMap["HLT_HT550_v11"]                = &fMT2tree->trigger.HLT_HT550_v11;
		fTriggerMap["HLT_HT550_v2"]                 = &fMT2tree->trigger.HLT_HT550_v2;
		fTriggerMap["HLT_HT550_v3"]                 = &fMT2tree->trigger.HLT_HT550_v3;
		fTriggerMap["HLT_HT550_v4"]                 = &fMT2tree->trigger.HLT_HT550_v4;
		fTriggerMap["HLT_HT550_v5"]                 = &fMT2tree->trigger.HLT_HT550_v5;
		fTriggerMap["HLT_HT550_v6"]                 = &fMT2tree->trigger.HLT_HT550_v6;
		fTriggerMap["HLT_HT550_v7"]                 = &fMT2tree->trigger.HLT_HT550_v7;
		fTriggerMap["HLT_HT550_v8"]                 = &fMT2tree->trigger.HLT_HT550_v8;
		fTriggerMap["HLT_HT550_v9"]                 = &fMT2tree->trigger.HLT_HT550_v9;
		fTriggerMap["HLT_HT600_v1"]                 = &fMT2tree->trigger.HLT_HT600_v1;
		fTriggerMap["HLT_HT600_v2"]                 = &fMT2tree->trigger.HLT_HT600_v2;
		fTriggerMap["HLT_HT600_v3"]                 = &fMT2tree->trigger.HLT_HT600_v3;
		fTriggerMap["HLT_HT600_v4"]                 = &fMT2tree->trigger.HLT_HT600_v4;
		fTriggerMap["HLT_HT650_v1"]                 = &fMT2tree->trigger.HLT_HT650_v1;
		fTriggerMap["HLT_HT650_v2"]                 = &fMT2tree->trigger.HLT_HT650_v2;
		fTriggerMap["HLT_HT650_v3"]                 = &fMT2tree->trigger.HLT_HT650_v3;
		fTriggerMap["HLT_HT650_v4"]                 = &fMT2tree->trigger.HLT_HT650_v4;
		fTriggerMap["HLT_HT700_v2"]                 = &fMT2tree->trigger.HLT_HT700_v2;
		fTriggerMap["HLT_HT750_v3"]                 = &fMT2tree->trigger.HLT_HT750_v3;
		fTriggerMap["HLT_HT750_L1FastJet_v3"]       = &fMT2tree->trigger.HLT_HT750_L1FastJet_v3;
		fTriggerMap["HLT_PFHT650_v1"]               = &fMT2tree->trigger.HLT_PFHT650_v1;
						            
		// MHT_HT			            
		fTriggerMap["HLT_HT250_MHT60_v2"]           = &fMT2tree->trigger.HLT_HT250_MHT60_v2;
		fTriggerMap["HLT_HT250_MHT60_v3"]           = &fMT2tree->trigger.HLT_HT250_MHT60_v3;
		fTriggerMap["HLT_HT250_MHT60_v4"]           = &fMT2tree->trigger.HLT_HT250_MHT60_v4;
		fTriggerMap["HLT_HT250_MHT60_v5"]           = &fMT2tree->trigger.HLT_HT250_MHT60_v5;
		fTriggerMap["HLT_HT250_MHT60_v6"]           = &fMT2tree->trigger.HLT_HT250_MHT60_v6;
		fTriggerMap["HLT_HT250_MHT70_v1"]           = &fMT2tree->trigger.HLT_HT250_MHT70_v1;
		fTriggerMap["HLT_HT250_MHT70_v2"]           = &fMT2tree->trigger.HLT_HT250_MHT70_v2;
		fTriggerMap["HLT_HT250_MHT70_v3"]           = &fMT2tree->trigger.HLT_HT250_MHT70_v3;
		fTriggerMap["HLT_HT250_MHT70_v4"]           = &fMT2tree->trigger.HLT_HT250_MHT70_v4;
		fTriggerMap["HLT_HT250_MHT90_v1"]           = &fMT2tree->trigger.HLT_HT250_MHT90_v1;
		fTriggerMap["HLT_HT250_MHT90_v2"]           = &fMT2tree->trigger.HLT_HT250_MHT90_v2;
		fTriggerMap["HLT_HT250_MHT100_v2"]          = &fMT2tree->trigger.HLT_HT250_MHT100_v2;
		fTriggerMap["HLT_HT260_MHT60_v2"]           = &fMT2tree->trigger.HLT_HT260_MHT60_v2;
		fTriggerMap["HLT_HT300_MHT75_v4"]           = &fMT2tree->trigger.HLT_HT300_MHT75_v4;
		fTriggerMap["HLT_HT300_MHT75_v5"]           = &fMT2tree->trigger.HLT_HT300_MHT75_v5;
		fTriggerMap["HLT_HT300_MHT75_v7"]           = &fMT2tree->trigger.HLT_HT300_MHT75_v7;
		fTriggerMap["HLT_HT300_MHT75_v8"]           = &fMT2tree->trigger.HLT_HT300_MHT75_v8;
		fTriggerMap["HLT_HT300_MHT80_v1"]           = &fMT2tree->trigger.HLT_HT300_MHT80_v1;
		fTriggerMap["HLT_HT300_MHT80_v2"]           = &fMT2tree->trigger.HLT_HT300_MHT80_v2;
		fTriggerMap["HLT_HT300_MHT90_v1"]           = &fMT2tree->trigger.HLT_HT300_MHT90_v1;
		fTriggerMap["HLT_HT300_MHT90_v2"]           = &fMT2tree->trigger.HLT_HT300_MHT90_v2;
		fTriggerMap["HLT_HT350_MHT70_v1"]           = &fMT2tree->trigger.HLT_HT350_MHT70_v1;
		fTriggerMap["HLT_HT350_MHT70_v2"]           = &fMT2tree->trigger.HLT_HT350_MHT70_v2;
		fTriggerMap["HLT_HT350_MHT80_v1"]           = &fMT2tree->trigger.HLT_HT350_MHT80_v1;
		fTriggerMap["HLT_HT350_MHT80_v2"]           = &fMT2tree->trigger.HLT_HT350_MHT80_v2;
		fTriggerMap["HLT_HT350_MHT90_v1"]           = &fMT2tree->trigger.HLT_HT350_MHT90_v1;
		fTriggerMap["HLT_HT350_MHT100_v3"]          = &fMT2tree->trigger.HLT_HT350_MHT100_v3;
		fTriggerMap["HLT_HT350_L1FastJet_MHT100_v1"]= &fMT2tree->trigger.HLT_HT350_L1FastJet_MHT100_v1;
		fTriggerMap["HLT_HT350_MHT110_v3"]          = &fMT2tree->trigger.HLT_HT350_MHT110_v3;
		fTriggerMap["HLT_HT400_MHT80_v1"]           = &fMT2tree->trigger.HLT_HT400_MHT80_v1;
		fTriggerMap["HLT_HT400_MHT90_v3"]           = &fMT2tree->trigger.HLT_HT400_MHT90_v3;
		fTriggerMap["HLT_PFHT350_PFMHT90_v1"]       = &fMT2tree->trigger.HLT_PFHT350_PFMHT90_v1;
		fTriggerMap["HLT_PFHT350_PFMHT100_v1"]      = &fMT2tree->trigger.HLT_PFHT350_PFMHT100_v1;
		fTriggerMap["HLT_PFHT400_PFMHT80_v1"]       = &fMT2tree->trigger.HLT_PFHT400_PFMHT80_v1;
		fTriggerMap["HLT_PFHT400_PFMHT90_v1"]       = &fMT2tree->trigger.HLT_PFHT400_PFMHT90_v1;
		// Muons			            
		fTriggerMap["HLT_DoubleMu3_HT160_v2"]       = &fMT2tree->trigger.HLT_DoubleMu3_HT160_v2;
		fTriggerMap["HLT_DoubleMu3_v3"]             = &fMT2tree->trigger.HLT_DoubleMu3_v3;
		fTriggerMap["HLT_Mu8_Jet40_v2"]             = &fMT2tree->trigger.HLT_Mu8_Jet40_v2;

                // **** MET Dataset ****
                // CentralJet+MET
                fTriggerMap["HLT_PFMHT150_v1" ]   = &fMT2tree->trigger.HLT_PFMHT150_v1 ;
                fTriggerMap["HLT_PFMHT150_v2" ]   = &fMT2tree->trigger.HLT_PFMHT150_v2 ;
                fTriggerMap["HLT_PFMHT150_v4" ]   = &fMT2tree->trigger.HLT_PFMHT150_v4 ;
                fTriggerMap["HLT_PFMHT150_v6" ]   = &fMT2tree->trigger.HLT_PFMHT150_v6 ;
                fTriggerMap["HLT_PFMHT150_v7" ]   = &fMT2tree->trigger.HLT_PFMHT150_v7 ;
                fTriggerMap["HLT_PFMHT150_v8" ]   = &fMT2tree->trigger.HLT_PFMHT150_v8 ;
                fTriggerMap["HLT_PFMHT150_v9" ]   = &fMT2tree->trigger.HLT_PFMHT150_v9 ;
                fTriggerMap["HLT_PFMHT150_v11"]   = &fMT2tree->trigger.HLT_PFMHT150_v11;
                fTriggerMap["HLT_PFMHT150_v12"]   = &fMT2tree->trigger.HLT_PFMHT150_v12;
                fTriggerMap["HLT_PFMHT150_v16"]   = &fMT2tree->trigger.HLT_PFMHT150_v16;
                fTriggerMap["HLT_CentralJet80_MET65_v1"]   = &fMT2tree->trigger.HLT_CentralJet80_MET65_v1;
                fTriggerMap["HLT_CentralJet80_MET65_v2"]   = &fMT2tree->trigger.HLT_CentralJet80_MET65_v2;
                fTriggerMap["HLT_CentralJet80_MET65_v3"]   = &fMT2tree->trigger.HLT_CentralJet80_MET65_v3;
                fTriggerMap["HLT_CentralJet80_MET65_v4"]   = &fMT2tree->trigger.HLT_CentralJet80_MET65_v4;
                fTriggerMap["HLT_CentralJet80_MET65_v5"]   = &fMT2tree->trigger.HLT_CentralJet80_MET65_v5;
                fTriggerMap["HLT_CentralJet80_MET65_v6"]   = &fMT2tree->trigger.HLT_CentralJet80_MET65_v6;
                fTriggerMap["HLT_CentralJet80_MET65_v7"]   = &fMT2tree->trigger.HLT_CentralJet80_MET65_v7;
                fTriggerMap["HLT_CentralJet80_MET65_v10"]  = &fMT2tree->trigger.HLT_CentralJet80_MET65_v10;
                fTriggerMap["HLT_CentralJet80_MET80_v1"]   = &fMT2tree->trigger.HLT_CentralJet80_MET80_v1;
                fTriggerMap["HLT_CentralJet80_MET80_v2"]   = &fMT2tree->trigger.HLT_CentralJet80_MET80_v2;
                fTriggerMap["HLT_CentralJet80_MET80HF_v2"] = &fMT2tree->trigger.HLT_CentralJet80_MET80HF_v2;
                fTriggerMap["HLT_CentralJet80_MET80HF_v3"] = &fMT2tree->trigger.HLT_CentralJet80_MET80HF_v3;
                fTriggerMap["HLT_CentralJet80_MET80HF_v4"] = &fMT2tree->trigger.HLT_CentralJet80_MET80HF_v4;
                fTriggerMap["HLT_CentralJet80_MET80HF_v5"] = &fMT2tree->trigger.HLT_CentralJet80_MET80HF_v5;
                fTriggerMap["HLT_CentralJet80_MET80_v6"]   = &fMT2tree->trigger.HLT_CentralJet80_MET80_v6;
                fTriggerMap["HLT_CentralJet80_MET80_v9"]   = &fMT2tree->trigger.HLT_CentralJet80_MET80_v9;
                fTriggerMap["HLT_CentralJet80_MET95_v3"]   = &fMT2tree->trigger.HLT_CentralJet80_MET95_v3;
                fTriggerMap["HLT_CentralJet80_MET110_v3"]  = &fMT2tree->trigger.HLT_CentralJet80_MET110_v3;
                fTriggerMap["HLT_CentralJet80_MET100_v1"]  = &fMT2tree->trigger.HLT_CentralJet80_MET100_v1;
                fTriggerMap["HLT_CentralJet80_MET100_v2"]  = &fMT2tree->trigger.HLT_CentralJet80_MET100_v2;
                fTriggerMap["HLT_CentralJet80_MET100_v3"]  = &fMT2tree->trigger.HLT_CentralJet80_MET100_v3;
                fTriggerMap["HLT_CentralJet80_MET100_v4"]  = &fMT2tree->trigger.HLT_CentralJet80_MET100_v4;
                fTriggerMap["HLT_CentralJet80_MET100_v5"]  = &fMT2tree->trigger.HLT_CentralJet80_MET100_v5;
                fTriggerMap["HLT_CentralJet80_MET100_v6"]  = &fMT2tree->trigger.HLT_CentralJet80_MET100_v6;
                fTriggerMap["HLT_CentralJet80_MET100_v7"]  = &fMT2tree->trigger.HLT_CentralJet80_MET100_v7;
                // DiCentralJet+MET
                fTriggerMap["HLT_DiCentralJet20_MET80_v1"] = &fMT2tree->trigger.HLT_DiCentralJet20_MET80_v1;
                fTriggerMap["HLT_DiCentralJet20_MET80_v2"] = &fMT2tree->trigger.HLT_DiCentralJet20_MET80_v2;
                fTriggerMap["HLT_DiCentralJet20_MET80_v3"] = &fMT2tree->trigger.HLT_DiCentralJet20_MET80_v3;
                fTriggerMap["HLT_DiCentralJet20_MET80_v4"] = &fMT2tree->trigger.HLT_DiCentralJet20_MET80_v4;
                fTriggerMap["HLT_DiCentralJet20_MET80_v5"] = &fMT2tree->trigger.HLT_DiCentralJet20_MET80_v5;
                fTriggerMap["HLT_DiCentralJet20_MET80_v8"] = &fMT2tree->trigger.HLT_DiCentralJet20_MET80_v8;
                fTriggerMap["HLT_DiCentralJet20_BTagIP_MET65_v2"]  = &fMT2tree->trigger.HLT_DiCentralJet20_BTagIP_MET65_v2;
                fTriggerMap["HLT_DiCentralJet20_BTagIP_MET65_v3"]  = &fMT2tree->trigger.HLT_DiCentralJet20_BTagIP_MET65_v3;
                fTriggerMap["HLT_DiCentralJet20_BTagIP_MET65_v4"]  = &fMT2tree->trigger.HLT_DiCentralJet20_BTagIP_MET65_v4;
                fTriggerMap["HLT_DiCentralJet20_BTagIP_MET65_v5"]  = &fMT2tree->trigger.HLT_DiCentralJet20_BTagIP_MET65_v5;
                fTriggerMap["HLT_DiCentralJet20_BTagIP_MET65_v6"]  = &fMT2tree->trigger.HLT_DiCentralJet20_BTagIP_MET65_v6;
                fTriggerMap["HLT_DiCentralJet20_BTagIP_MET65_v7"]  = &fMT2tree->trigger.HLT_DiCentralJet20_BTagIP_MET65_v7;
                fTriggerMap["HLT_DiCentralJet20_BTagIP_MET65_v10"] = &fMT2tree->trigger.HLT_DiCentralJet20_BTagIP_MET65_v10;
                fTriggerMap["HLT_DiCentralJet20_BTagIP_MET65_v11"] = &fMT2tree->trigger.HLT_DiCentralJet20_BTagIP_MET65_v11;
                fTriggerMap["HLT_DiCentralJet20_MET100_HBHENoiseFiltered_v1"] = &fMT2tree->trigger.HLT_DiCentralJet20_MET100_HBHENoiseFiltered_v1;
                fTriggerMap["HLT_DiCentralJet20_MET100_HBHENoiseFiltered_v4"] = &fMT2tree->trigger.HLT_DiCentralJet20_MET100_HBHENoiseFiltered_v4;
                fTriggerMap["HLT_DiCentralPFJet30_PFMHT80_v1"] = &fMT2tree->trigger.HLT_DiCentralPFJet30_PFMHT80_v1;
                fTriggerMap["HLT_DiCentralPFJet50_PFMHT80_v1"] = &fMT2tree->trigger.HLT_DiCentralPFJet50_PFMHT80_v1;
                fTriggerMap["HLT_DiJet60_MET45_v1" ]  = &fMT2tree->trigger.HLT_DiJet60_MET45_v1 ;
                fTriggerMap["HLT_DiJet60_MET45_v2" ]  = &fMT2tree->trigger.HLT_DiJet60_MET45_v2 ;
                fTriggerMap["HLT_DiJet60_MET45_v3" ]  = &fMT2tree->trigger.HLT_DiJet60_MET45_v3 ;
                fTriggerMap["HLT_DiJet60_MET45_v4" ]  = &fMT2tree->trigger.HLT_DiJet60_MET45_v4 ;
                fTriggerMap["HLT_DiJet60_MET45_v5" ]  = &fMT2tree->trigger.HLT_DiJet60_MET45_v5 ;
                fTriggerMap["HLT_DiJet60_MET45_v6" ]  = &fMT2tree->trigger.HLT_DiJet60_MET45_v6 ;
                fTriggerMap["HLT_DiJet60_MET45_v7" ]  = &fMT2tree->trigger.HLT_DiJet60_MET45_v7 ;
                fTriggerMap["HLT_DiJet60_MET45_v10"]  = &fMT2tree->trigger.HLT_DiJet60_MET45_v10;
	}


	// initialize fDeadCellFilterBE and fDeadCellFilterTP
	fDeadCellFilterBE.Reset();
	fDeadCellFilterTP.Reset();
	for(int i=0; i<fTPfiles.size(); ++i){
		MT2Analysis::DeadCellParser(fDeadCellFilterTP, fTPfiles[i]);
	}
	for(int i=0; i<fBEfiles.size(); ++i){
		MT2Analysis::DeadCellParser(fDeadCellFilterBE, fBEfiles[i]);
	}

	if(fJEC.length()!=0){
	// initialize JEC and JESuncertainty
	cout << "--------------------- " << endl;
	cout << " -> initialize JetEnergyCorrection & JetCorrectionUncertainty " << endl;
	Initialize_JetEnergyCorrection();
	Initialize_JetCorrectionUncertainty();
	cout << "--------------------- " << endl;
	}else{
		fJecUncPF         =NULL;
		fJecUncCalo       =NULL;
	        fJecCaloL2        =NULL;
	        fJecCaloL3        =NULL;
	        fJecPFL1          =NULL;
	        fJecPFL2          =NULL; 
	        fJecPFL3          =NULL;
	        fJecPFRES         =NULL;
	        fJetCorrectorPF   =NULL;
	        fJetCorrectorCalo =NULL;
	}
	if(fDoJESUncertainty && fJEC.length()==0){
		cout << "ERROR: need to know JEC set to upscale/downscale the Jet Energy" << endl;
		exit(-1);
	}

	//LHAPDF init
        string PDF_SET="cteq66";
        string PDF_PATH = "/shome/leo/Installations/LHAPDF/lhapdf-5.8.4/";
	if(doPDF){
	  LHAPDF::initPDFSet(PDF_PATH+"/share/lhapdf/"+PDF_SET, LHAPDF::LHGRID);
	  nPDFs = LHAPDF::numberPDF();
	  cout << "nPDF: " << nPDFs << endl;
	}




}

// ***********************************************************************
// Analyze called for each event
void MT2Analysis::Analyze(){	

	// ---------------------------------------------------
	// Initialize fElecs, fJetsLoose, fBJets, fMuons, fLeptConfig 
	InitializeEvent();

	// ------------------------------------------------------------------
	// check if event passes selection cuts and trigger requirements
	if(! IsSelectedEvent()){return;}

	// -------------------------------------------------------------------
	// calculate MT2tree variables and store it in output TTree
	
	// reset tree variables
	ResetTree();

	// Fill Tree !! has to be called at the very end !!
	bool basicfilled = FillMT2TreeBasics();

	// Fill complicated MT2tree variables;
	if(basicfilled) FillMT2treeCalculations();

	// FillTree
	FillTree();
}

// ***********************************************************************************************
// fill MT2 tree
void MT2Analysis::BookTree(){
	fATree = new TTree("MassTree", "MassTree");
	fATree->Branch("MT2tree" , "MT2tree" , &fMT2tree);
}

void MT2Analysis::ResetTree(){
        fMT2tree->Reset();
}

void MT2Analysis::FillTree(){
	// fill output tree
	fATree            ->Fill();
}

bool MT2Analysis::FillMT2TreeBasics(){
	// check size of jets electrons muons and genleptons
	if(fJets.size()     > 25) {cout << "ERROR: fJets.size()   > 25: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return false;}
	if(fElecs.size()    > 5 ) {cout << "ERROR: fElecs.size()  >  5: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return false;}
	if(fMuons.size()    > 5 ) {cout << "ERROR: fMuons.size()  >  5: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return false;}
	if(fPhotons.size()  > 5 ) {cout << "ERROR: fPhotons.size()>  5: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return false;}

	//pdf weights
	if(doPDF){
          fMT2tree->NPdfs = nPDFs;
          fMT2tree->pdfW[0]=1;
	  LHAPDF::initPDF(0);

          float pdf01 = LHAPDF::xfx(fTR->PDFx1, fTR->PDFScalePDF, fTR->PDFID1)/fTR->PDFx1 ;
          float pdf02 = LHAPDF::xfx(fTR->PDFx2, fTR->PDFScalePDF, fTR->PDFID2)/fTR->PDFx2 ;

          for(int pdf=1; pdf<= nPDFs; pdf++){
	    LHAPDF::initPDF(pdf);
            float pdf1 = LHAPDF::xfx(fTR->PDFx1, fTR->PDFScalePDF, fTR->PDFID1)/fTR->PDFx1 ;
            float pdf2 = LHAPDF::xfx(fTR->PDFx2, fTR->PDFScalePDF, fTR->PDFID2)/fTR->PDFx2 ;
            fMT2tree->pdfW[pdf] = pdf1/pdf01*pdf2/pdf02;
          }
        }

	//MC info
	if(!fisData){
	  fMT2tree->GenProcessID = fTR->process;
	  fMT2tree->GenWeight = fTR->GenWeight;
	}
	if(isScan){
	  fMT2tree->Susy.MassGlu = fTR->MassGlu;
	  fMT2tree->Susy.MassChi= fTR->MassChi;
	  fMT2tree->Susy.MassLSP= fTR->MassLSP;
	  fMT2tree->Susy.M0= fTR->M0;
	  fMT2tree->Susy.M12= fTR->M12;
	  fMT2tree->Susy.A0= fTR->A0;
	  fMT2tree->Susy.Mu= fTR->signMu;
	  fMT2tree->Susy.XSec = fTR->IntXSec;
	}

	// ---------------------------------------------------------------
	// Fill jets 4-momenta & ID's 
	for(int i=0; i<fJets.size(); ++i) {
		fMT2tree->jet[i].lv.SetPtEtaPhiE( Jet(fJets[i]).Pt(),  Jet(fJets[i]).Eta(), 
						  Jet(fJets[i]).Phi(), Jet(fJets[i]).E()  ); //SetLV(GetJet4Momenta(fJets[i]));
		// b-tag info now should be available
		fMT2tree->jet[i].bTagProbTCHE  =  fTR->PF2PAT3JbTagProbTkCntHighEff [fJets[i]];
		fMT2tree->jet[i].bTagProbTCHP  =  fTR->PF2PAT3JbTagProbTkCntHighPur [fJets[i]];
		fMT2tree->jet[i].bTagProbSSVHE =  fTR->PF2PAT3JbTagProbSimpSVHighEff[fJets[i]];
		fMT2tree->jet[i].bTagProbSSVHP =  fTR->PF2PAT3JbTagProbSimpSVHighPur[fJets[i]];
	  
		// Jet id variables
		if(IsGoodBasicPFJetPAT3 (fJets[i], 20., 2.4))  fMT2tree->jet[i].isPFIDLoose   =true;
		if(IsGoodPFJetMediumPAT3(fJets[i], 20., 2.4))  fMT2tree->jet[i].isPFIDMedium  =true;
		if(IsGoodPFJetTightPAT3 (fJets[i], 20., 2.4))  fMT2tree->jet[i].isPFIDTight   =true;
		if(fTR->PF2PAT3JIDLoose [fJets[i]]==1       )  fMT2tree->jet[i].isPATPFIDLoose=true;;  // PF jetID regardless of eta and pt
		fMT2tree->jet[i].ChHadFrac      = fTR->PF2PAT3JChHadfrac     [fJets[i]];	
		fMT2tree->jet[i].NeuHadFrac     = fTR->PF2PAT3JNeuHadfrac    [fJets[i]];
		fMT2tree->jet[i].ChEmFrac       = fTR->PF2PAT3JChEmfrac      [fJets[i]];
		fMT2tree->jet[i].NeuEmFrac      = fTR->PF2PAT3JNeuEmfrac     [fJets[i]];
		fMT2tree->jet[i].ChMult         = fTR->PF2PAT3JChMult        [fJets[i]];
		fMT2tree->jet[i].NeuMult        = fTR->PF2PAT3JNeuMult       [fJets[i]];
		fMT2tree->jet[i].NConstituents  = fTR->PF2PAT3JNConstituents [fJets[i]];
		fMT2tree->jet[i].Area           = fTR->PF2PAT3JArea          [fJets[i]];
		fMT2tree->jet[i].Flavour        = (fTR->fChain->GetBranch("PF2PAT3JFlavour")==NULL || fisData)? -77777.77: fTR->PF2PAT3JFlavour[fJets[i]];        // Branch added in  V02-04-01
		// JET energy correction factors
		// full correction factor
		if      (fJEC.length()!=0) fMT2tree->jet[i].Scale = GetPFJEC(fTR->PF2PAT3JPt   [fJets[i]], fTR->PF2PAT3JScale[fJets[i]], fTR->PF2PAT3JEta[fJets[i]], fTR->PF2PAT3JArea[fJets[i]], fTR->Rho, fisData?1230:123);
		else                       fMT2tree->jet[i].Scale = fTR->PF2PAT3JScale[fJets[i]];
		// L1FastScale
		if      (fTR->fChain->GetBranch("PF2PAT3JL1FastJetScale")==NULL) fMT2tree->jet[i].L1FastJetScale=-77777.77;  // Branch added in  V02-04-01
		else if (fJEC.length()!=0) fMT2tree->jet[i].L1FastJetScale=GetPFJEC(fTR->PF2PAT3JPt   [fJets[i]], fTR->PF2PAT3JScale[fJets[i]], fTR->PF2PAT3JEta[fJets[i]], fTR->PF2PAT3JArea[fJets[i]], fTR->Rho, 1);
		else                       fMT2tree->jet[i].L1FastJetScale=fTR->PF2PAT3JL1FastJetScale[fJets[i]];
	}

	// --------------------------------------------------------------
	// match taus to jets
	for(int i=0; i<fTaus.size(); ++i){
		TLorentzVector tau;
		tau.SetPxPyPzE(fTR->PfTau3Px[fTaus[i]],fTR->PfTau3Py[fTaus[i]],fTR->PfTau3Pz[fTaus[i]],fTR->PfTau3E[fTaus[i]]);
		float mindR =10000; int jindex =-1;
		for(int j=0; j<fJets.size(); ++j){
			float dR = tau.DeltaR(fMT2tree->jet[j].lv);
			if(dR < mindR && dR < 0.5) {mindR=dR; jindex = j;} 	
		}
		if(jindex ==-1) continue;
		fMT2tree->jet[jindex].NTauMatch++;
		if(fMT2tree->jet[jindex].isTauMatch && fMT2tree->jet[jindex].TauDR < mindR) continue;
		fMT2tree->jet[jindex].isTauMatch = 1;
		fMT2tree->jet[jindex].TauDR      = mindR;
		fMT2tree->jet[jindex].TauDPt     = fMT2tree->jet[jindex].lv.Pt()-tau.Pt();
	}

	// -----------------------------------------------------------------
	// Photons
	for(int i=0; i<fPhotons.size(); ++i){
		fMT2tree->photon[i].lv.SetPtEtaPhiM(fTR->PhoPt[fPhotons[i]], fTR->PhoEta[fPhotons[i]], fTR->PhoPhi[fPhotons[i]], 0.);
		fMT2tree->photon[i].TrkIso              =fTR->PhoIso04TrkHollow[fPhotons[i]]; 
		fMT2tree->photon[i].EcalIso             =fTR->PhoIso04Ecal[fPhotons[i]];
		fMT2tree->photon[i].HcalIso             =fTR->PhoIso04Hcal[fPhotons[i]];
		fMT2tree->photon[i].SigmaIEtaIEta       =fTR->PhoSigmaIetaIeta[fPhotons[i]];
		fMT2tree->photon[i].HoverE              =fTR->PhoHoverE[fPhotons[i]];
		fMT2tree->photon[i].isEGMlooseID        =IsGoodPhotonEGMLooseID(fPhotons[i]);
		fMT2tree->photon[i].isEGMlooseIso       =IsGoodPhotonEGMLooseISO(fPhotons[i]);
		fMT2tree->photon[i].isEGMlooseRelIso    =IsGoodPhotonEGMLooseRelISO(fPhotons[i]);
		fMT2tree->photon[i].isEGMtightID        =IsGoodPhotonEGMTightID(fPhotons[i]);
		fMT2tree->photon[i].isEGMtightIso       =IsGoodPhotonEGMTightISO(fPhotons[i]);
		fMT2tree->photon[i].JetRemoved          =fPhotonJetOverlapRemoved[i];
		fMT2tree->photon[i].MCmatchexitcode     =fTR->PhoMCmatchexitcode[fPhotons[i]]; 
		if(fTR->PhoMCmatchindex[fPhotons[i]]>=0){
		int index=fTR->PhoMCmatchindex[fPhotons[i]];
		fTR->fChain->GetBranch("GenPhotonPartonMindR")==NULL ? fMT2tree->photon[i].GenJetMinDR= -777.77 :
			                                               fMT2tree->photon[i].GenJetMinDR= fTR->GenPhotonPartonMindR[index];
		}
		// exit code meaning:  
		//                     0 = matched to particle that is not a photon -> fake
		//                     1 = photon  with status 3 quark or gluon as mother (hard scatter) (gluon happens, though gamma does not have color..)
		//                     2 = PromptPhoton with status 3 photon as mother,  
		//                     3 = other particle, i.e. showering.. 
		//                     negative values: no MC match. 
	}
	
	// ---------------------------------------------------------------
	// Set NJets, NElecs, NMuons
	fMT2tree->SetNJets         (fJets.size());
	fMT2tree->SetNGenJets      (fTR->NGenJets > gNGenJets ? gNGenJets: fTR->NGenJets);
	fMT2tree->SetNJetsIDLoose  (fMT2tree->GetNjets(20, 2.4, 1));
	fMT2tree->SetNBJets        (fMT2tree->GetNBtags(3,2.0,20,2.4,1));
	fMT2tree->SetNEles         ((Int_t)fElecs.size());
	fMT2tree->SetNMuons        ((Int_t)fMuons.size());
	fMT2tree->SetNTaus         ((Int_t)fTaus.size());
	fMT2tree->SetNPhotons      ((Int_t)fPhotons.size());
	fMT2tree->NJetsIDLoose40 = fMT2tree->GetNjets(40, 2.4, 1);
	fMT2tree->NJetsIDLoose50 = fMT2tree->GetNjets(50, 2.4, 1);
	
	// --------------------------------------------------------------
	// Fill GenJets
	for(int i=0; i<fTR->NGenJets; ++i){
		if(i >= gNGenJets) {
			cout << "WARNING: Event " << fTR->Event << " Run " << fTR->Run << " has more than " << gNGenJets << " GenJets " << endl;
			continue;
	       	}
		fMT2tree->genjet[i].lv.SetPtEtaPhiE(fTR->GenJetPt[i], fTR->GenJetEta[i], fTR->GenJetPhi[i], fTR->GenJetE[i]);
		float mindR=999.99;
		int    index=-1;
		for(int j=0; j<fMT2tree->NJets; ++j){
			float dR=fMT2tree->jet[j].lv.DeltaR(fMT2tree->genjet[i].lv);
			if(dR < mindR) {
				mindR = dR;
				index = j;
			}
		}
		fMT2tree->genjet[i].DeltaR        = mindR;
		fMT2tree->genjet[i].JetMatchIndex = index;
	}
	
	
	// -----------------------------------------------------------------
	// Fill leptons 4-momenta
	TLorentzVector METlv;
	METlv.SetPtEtaPhiM(MET().Pt(), 0., MET().Phi(), 0.);

	for(int i=0; i<fElecs.size(); ++i) {
	  	fMT2tree->ele[i].lv.SetPtEtaPhiE(fTR->PfEl3Pt [fElecs[i]], fTR->PfEl3Eta[fElecs[i]], 
				          fTR->PfEl3Phi[fElecs[i]], fTR->PfEl3E  [fElecs[i]]); // = GetEle4Momenta(fElecs[i]);
		fMT2tree->ele[i].MT     = fMT2tree->GetMT(fMT2tree->ele[i].lv, 0., METlv, 0.); 
		fMT2tree->ele[i].Charge = fTR->PfEl3Charge[fElecs[i]];
		fMT2tree->ele[i].Iso    = (fTR->PfEl3ChargedHadronIso[fElecs[i]] + fTR->PfEl3NeutralHadronIso[fElecs[i]] + fTR->PfEl3PhotonIso[fElecs[i]])/fTR->PfEl3Pt[fElecs[i]];
		fMT2tree->ele[i].ID90   = fTR->PfEl3ID90[fElecs[i]];
		fMT2tree->ele[i].ID95   = fTR->PfEl3ID95[fElecs[i]];
		
		// match to caloJets
		for(int j=0; j<fTR->CANJets; ++j){
			if(fMT2tree->ele[i].lv.DeltaR(CAJet(j)) <0.4) {
				fMT2tree->ele[i].CAJ_n90     = fTR->CAJn90[j];
				fMT2tree->ele[i].CAJ_n90Hits = fTR->CAJID_n90Hits[j];
			}
		}
	}
	for(int i=0; i<fMuons.size(); ++i) {
	  	fMT2tree->muo[i].lv.SetPtEtaPhiM(fTR->PfMu3Pt [fMuons[i]], fTR->PfMu3Eta[fMuons[i]], 
				          fTR->PfMu3Phi[fMuons[i]], 0.106);                     // = GetMuo4Momenta(fMuons[i]);
		fMT2tree->muo[i].MT       = fMT2tree->GetMT(fMT2tree->muo[i].lv, fMT2tree->muo[i].lv.M(), METlv, 0.); 
		fMT2tree->muo[i].Charge   = fTR->PfMu3Charge[fMuons[i]];	
		fMT2tree->muo[i].Iso      = (fTR->PfMu3ChargedHadronIso[fMuons[i]] + fTR->PfMu3NeutralHadronIso[fMuons[i]] + fTR->PfMu3PhotonIso[fMuons[i]])/fTR->PfMu3Pt[fMuons[i]];
		fMT2tree->muo[i].NMatches = fTR->PfMu3NMatches[fMuons[i]];
		fMT2tree->muo[i].PtErr    = fTR->PfMu3PtErr[fMuons[i]];
	}
	
	// ---------------------------------------------------------------
	// GenMET	
	fMT2tree->genmet[0].SetPtEtaPhiM(fTR->GenMET, 0., fTR->GenMETphi, 0.);
	// -------------------------------------------------------------------
	// Genleptons
	int NGenLepts=0;
	for(int i=0; i<fTR->NGenLeptons; ++i){
		int ID=abs(fTR->GenLeptonID[i]);
		int MID=abs(fTR->GenLeptonMID[i]);
		if( (ID == 11 || ID == 12 || ID == 13 || ID == 14 || ID == 15 || ID == 16) && (MID > 24) ) continue;
		float mass=0;
		if     (abs(fTR->GenLeptonID[i]) == 15)   mass=1.776; // tau
		else if(abs(fTR->GenLeptonID[i]) == 13)   mass=0.106; // mu
		else if(abs(fTR->GenLeptonID[i]) == 11)   mass=0.;    // el 
		else if(abs(fTR->GenLeptonID[i]) == 12 || 
			abs(fTR->GenLeptonID[i]) == 14 || 
			abs(fTR->GenLeptonID[i]) == 16)   mass=0.;    // nu 
		else if(abs(fTR->GenLeptonID[i]) == 5 )   mass=4.2;   // bottom-quark
		else if(abs(fTR->GenLeptonID[i]) == 22 )  mass=0;     // photon
		else if(abs(fTR->GenLeptonID[i]) == 23 )  mass=91.2;   //  Z
		else if(abs(fTR->GenLeptonID[i]) == 24 )  mass=80.4;   //  W
		else   continue;
		NGenLepts++;
		if(NGenLepts >= 20 ) {cout << "ERROR: NGenLepts >=20: skipping remaining genlepts for event " << fTR->Event << endl; continue;}
		fMT2tree->genlept[NGenLepts-1].lv.SetPtEtaPhiM(fTR->GenLeptonPt[i], fTR->GenLeptonEta[i], fTR->GenLeptonPhi[i], mass);
		fMT2tree->genlept[NGenLepts-1].ID       = fTR->GenLeptonID[i];
		fMT2tree->genlept[NGenLepts-1].MID      = fTR->GenLeptonMID[i];
		fMT2tree->genlept[NGenLepts-1].MStatus  = fTR->GenLeptonMStatus[i];
		fMT2tree->genlept[NGenLepts-1].GMID     = fTR->GenLeptonGMID[i];
		fMT2tree->genlept[NGenLepts-1].GMStatus = fTR->GenLeptonGMStatus[i];
		if(abs(fMT2tree->genlept[NGenLepts-1].ID) == 11 || abs(fMT2tree->genlept[NGenLepts-1].ID) == 13 || abs(fMT2tree->genlept[NGenLepts-1].ID) == 15   ){
			fMT2tree->genlept[NGenLepts-1].MT = fMT2tree->GetMT(fMT2tree->genlept[NGenLepts-1].lv, fMT2tree->genlept[NGenLepts-1].lv.M(), fMT2tree->genmet[0], 0.);
		}
		if(abs(fTR->GenLeptonID[NGenLepts-1]) == 11){// match to caloJets
			for(int j=0; j<fTR->CANJets; ++j){
				if(fMT2tree->genlept[NGenLepts-1].lv.DeltaR(CAJet(j)) <0.4) {
					fMT2tree->genlept[NGenLepts-1].CAJ_n90     = fTR->CAJn90[j];
					fMT2tree->genlept[NGenLepts-1].CAJ_n90Hits = fTR->CAJID_n90Hits[j];
				}
			}
		}
	}
	// add GenPhotons
	fMT2tree->NGenLepts = NGenLepts;


	// --------------------------------------------------------------------
	// MET 
	fMT2tree->pfmet[0]=MET();
	// raw, uncorrected and not modified pfmet
	fMT2tree->rawpfmet[0].SetPtEtaPhiM(fTR->PFMETPAT, 0, fTR->PFMETPATphi, 0);
	// raw calomet
	fMT2tree->misc.CaloMETRaw = fTR->RawMET;
	fMT2tree->misc.CaloMETMuJesCorr = fTR->MuJESCorrMET;

	// Pile UP info and reco vertices
	if(!fisData){
		fMT2tree->pileUp.PUnumInt          = fTR->PUnumInteractions;        
		fMT2tree->pileUp.PUnumIntLate      = (fTR->fChain->GetBranch("PUOOTnumInteractionsLate") ==NULL) ? -77777.77: fTR->PUOOTnumInteractionsLate;   // branch added in V02-03-01 
		fMT2tree->pileUp.PUnumIntEarly     = (fTR->fChain->GetBranch("PUOOTnumInteractionsEarly")==NULL) ? -77777.77: fTR->PUOOTnumInteractionsEarly;  // branch added in V02-03-01 
		fMT2tree->pileUp.PtHat             = fTR->PtHat;
		fMT2tree->pileUp.isS3              = (int) isS3;
		//////////// S3 vs S4
		if(noPU && isS3){
		  fMT2tree->pileUp.Weight=GetPUWeight3D(fTR->PUOOTnumInteractionsEarly ,fTR->PUnumInteractions , fTR->PUOOTnumInteractionsLate);
		}
		else if(noPU)
		  fMT2tree->pileUp.Weight            = 1;
	        else if(isS3)
		  fMT2tree->pileUp.Weight            = (fTR->fChain->GetBranch("PUOOTnumInteractionsLate")==NULL) ? 1: GetPUWeight(fTR->PUnumInteractions, fTR->PUOOTnumInteractionsLate); // branch added in V02-03-01 
		else
		  fMT2tree->pileUp.Weight            = GetPUWeight(fTR->PUnumInteractions);
		
		if(fVerbose > 3) {
			cout << "fTR->PUnumInteractions " <<  fTR->PUnumInteractions << " weight "  
		     	     << " GetPUWeight() " << GetPUWeight(fTR->PUnumInteractions) << endl; 
		}
	}
	fMT2tree->pileUp.Rho               = fTR->Rho; // ATTENTION: this rho is from KT6 PF jets without pf-CHS
	//cout << fTR->Rho << " " <<  fMT2tree->pileUp.Rho << endl;

	int nvertex=0;
	for(int i=0; i<fTR->NVrtx; ++i){
		if(fabs(fTR->VrtxZ[i]) > 24) continue;
		if(sqrt( (fTR->VrtxX[i])*(fTR->VrtxX[i]) + (fTR->VrtxY[i])*(fTR->VrtxY[i])) > 2) continue;
		if(fTR->VrtxNdof[i]<=4) continue;
		nvertex++;
	}

	fMT2tree->pileUp.NVertices=nvertex;
	// _________

	// _________
	// HLT triggers --------------------------------------------------------------------------------
	for (StringBoolMap::iterator iter = fTriggerMap.begin(); iter != fTriggerMap.end(); ++iter){
		if(GetHLTResult(iter->first)) *iter->second =1;
	}
	if(fisData){
		// single photon triggers ---------------------------------------------------------------------
		string singlePhotonTriggers[100];
		int photontriggernumber=0;
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_IsoL_v1";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_IsoL_v2";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_IsoL_v3";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_IsoL_v4";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_IsoL_v5";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_IsoL_v6";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_IsoL_v7";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_IsoL_v8";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_IsoL_v9";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_IsoL_v10";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_v1";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_v2";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_v3";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_v4";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_v5";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_v6";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon75_CaloIdVL_v7";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon90_CaloIdVL_IsoL_v1";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon90_CaloIdVL_IsoL_v2";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon90_CaloIdVL_IsoL_v3";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon90_CaloIdVL_IsoL_v4";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon90_CaloIdVL_IsoL_v5";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon90_CaloIdVL_IsoL_v6";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon90_CaloIdVL_IsoL_v7";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon90_CaloIdVL_v1";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon90_CaloIdVL_v2";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon90_CaloIdVL_v3";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon90_CaloIdVL_v4";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon125_v1";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon125_v1";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon135_v1";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon135_v2";

		bool SiglePhotFired(false);
		for(int i=0; i<photontriggernumber; ++i){
			if( GetHLTResult(singlePhotonTriggers[i])) SiglePhotFired=true;
		}
		if(SiglePhotFired) fMT2tree->trigger.HLT_SinglePhotons = true;
		//NOTE NEW START
		string MuHadTriggers[100];
		int muhadtriggernumber=0;
		MuHadTriggers[muhadtriggernumber++] = "HLT_Mu8_HT200_v2";
		MuHadTriggers[muhadtriggernumber++] = "HLT_Mu8_HT200_v3";
		MuHadTriggers[muhadtriggernumber++] = "HLT_Mu8_HT200_v4";
		MuHadTriggers[muhadtriggernumber++] = "HLT_Mu15_HT200_v2";
		MuHadTriggers[muhadtriggernumber++] = "HLT_Mu15_HT200_v3";
		MuHadTriggers[muhadtriggernumber++] = "HLT_Mu15_HT200_v4";
		MuHadTriggers[muhadtriggernumber++] = "HLT_Mu30_HT200_v1";
		MuHadTriggers[muhadtriggernumber++] = "HLT_Mu30_HT200_v3";
		MuHadTriggers[muhadtriggernumber++] = "HLT_Mu40_HT200_v4";
		MuHadTriggers[muhadtriggernumber++] = "HLT_Mu40_HT300_v4";
		MuHadTriggers[muhadtriggernumber++] = "HLT_Mu40_HT300_v5";

		bool MuHadFired(false);
		for(int i=0; i<muhadtriggernumber; ++i){
			if( GetHLTResult(MuHadTriggers[i])) MuHadFired=true;
		}
		if(MuHadFired) fMT2tree->trigger.HLT_MuHad = true;
		//NOTE NEW END

		// di-electron triggers ------------------------------------------------------------------------------------------------------
		string diElectronTriggers[100];
		int diElectronTriggernumber=0;
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v7";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v8";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v1";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v5";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10";
		bool DiElectronFired(false);
		for(int i=0; i<diElectronTriggernumber; ++i){
			if(GetHLTResult(diElectronTriggers[i])) DiElectronFired=true;
		}
		if(DiElectronFired) fMT2tree->trigger.HLT_DiElectrons =true;

		// DiMuon triggers
		string diMuonTriggers[100];
		int diMuonTriggernumber=0;
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu6_v1";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu6_v2";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu6_v3";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu6_v4";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu6_v5";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu6_v6";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu6_v7";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu6_v8";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu7_v1";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu7_v2";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu7_v3";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu7_v4";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu7_v5";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu7_v6";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu7_v7";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu7_v8";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu7_v9";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu7_v10";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu7_v11";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_DoubleMu7_v12";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v1";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v2";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v3";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v4";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v5";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v6";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v7";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v8";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v9";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v10";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v11";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v1";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v2";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v3";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v4";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v5";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v6";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v7";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v8";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v9";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v10";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v11";
		bool DiMuonFired(false);
		for(int i=0; i<diMuonTriggernumber; ++i){
			if(GetHLTResult(diMuonTriggers[i])) DiMuonFired=true;
		}
		if(DiMuonFired) fMT2tree->trigger.HLT_DiMuons =true;

		//EMu
		string EMuTriggers[100];
		int EMuTriggernumber=0;
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdL_v1";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdL_v2";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdL_v3";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdL_v4";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdL_v5";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdL_v6";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdL_v7";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdL_v8";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdL_v9";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdL_v1";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdL_v2";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdL_v3";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdL_v4";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdL_v5";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdL_v6";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdL_v7";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdL_v8";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdL_v9";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v1";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v2";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v5";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v6";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v1";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v2";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v3";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v5";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v6";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8";
		bool EMuFired(false);
		for(int i=0; i<EMuTriggernumber; ++i){
			if(GetHLTResult(EMuTriggers[i])) EMuFired=true;
		}
		if(EMuFired) fMT2tree->trigger.HLT_EMu =true;
	}



	// ___________________________________________________________________________
	
	//__________
	// ECAL Dead Cell
	for(int i=0;i<fDeadCellFilterBE.event.size(); ++i){
		if(fTR->Run         !=fDeadCellFilterBE.run[i] )  continue;
		if(fTR->LumiSection !=fDeadCellFilterBE.lumi[i])  continue;
		if(fTR->Event       !=fDeadCellFilterBE.event[i]) continue;
		fMT2tree->misc.BadEcalBE = 1;
	}
	for(int i=0;i<fDeadCellFilterTP.event.size(); ++i){
		if(fTR->Run         !=fDeadCellFilterTP.run[i] )  continue;
		if(fTR->LumiSection !=fDeadCellFilterTP.lumi[i])  continue;
		if(fTR->Event       !=fDeadCellFilterTP.event[i]) continue;
		fMT2tree->misc.BadEcalTP = 1;
	}
	if (fTR->fChain->GetBranch("EcalDeadTPFilterFlag")!=NULL){            // Branch added in  V02-04-00
	       	fMT2tree->misc.BadEcalTP = (fTR->EcalDeadTPFilterFlag==1) ? 0:1; // store bad events as "true" 
	}
	
	// ------------------------------------------------------------------
	// fill misc 
	fMT2tree->misc.isData              = fisData;
	fMT2tree->misc.ProcessID           = fID;
	fMT2tree->misc.Run                 = fTR->Run;
	fMT2tree->misc.Event		   = fTR->Event;
	fMT2tree->misc.LumiSection	   = fTR->LumiSection;
	fMT2tree->misc.LeptConfig          = (int) fLeptConfig;
	fMT2tree->misc.HBHENoiseFlag	   = (fTR->HBHENoiseFlag==1) ? 0:1;  // store bad events as "true"
	fMT2tree->misc.CrazyHCAL           = fCrazyHCAL;                     // store bad events as "true"
	fMT2tree->misc.NegativeJEC         = fNegativeJEC;
	fMT2tree->misc.HT                  = fHT;
	
	fMT2tree->misc.MET                 = MET().Pt();
	fMT2tree->misc.METPhi              = MET().Phi();

	fMT2tree->misc.LeadingJPt          = (fMT2tree->NJets > 0) ? fMT2tree->jet[0].lv.Pt() : 0;
	fMT2tree->misc.SecondJPt           = (fMT2tree->NJets > 1) ? fMT2tree->jet[1].lv.Pt() : 0;
	
	// RA2 tracking failure
	fMT2tree->misc.TrackingFailure     = fTR->TrkPtSum/fHT;
	fMT2tree->misc.TrackingFailurePVtx = fTR->PrimVtxPtSum/fHT;

	if(fTR->fChain->GetBranch("QCDPartonicHT")        !=NULL) fMT2tree->misc.QCDPartonicHT         = fTR->QCDPartonicHT;                   // Branch added in ntuple V02-04-07
	if(fTR->fChain->GetBranch("CSCTightHaloID")       !=NULL) fMT2tree->misc.CSCTightHaloID        = fTR->CSCTightHaloID;                  // Branch added in ntuple V02-04-00
	if(fTR->fChain->GetBranch("RecovRecHitFilterFlag")!=NULL) fMT2tree->misc.RecovRecHitFilterFlag = (fTR->RecovRecHitFilterFlag==1)? 0:1; // Branch added in ntuple V02-04-06	
	if(fTR->fChain->GetBranch("HBHENoiseFlagIso")     !=NULL) fMT2tree->misc.HBHENoiseFlagIso      = (fTR->HBHENoiseFlagIso==1)     ? 0:1; // Branch added in ntuple V02-04-06 
	
	// add gen photon
	if(fTR->fChain->GetBranch("NGenPhotons") !=NULL && !fisData){
		TLorentzVector photon(0,0,0,0);			
		for(int i=0; i<fTR->NGenPhotons; ++i){
			// mother status ==3 means mother takes part in hard statter: 
			// this means the mother can be photon itself or e.g. quark
			if(fTR->GenPhotonMotherStatus[i]==3) { 
				photon.SetPtEtaPhiM(fTR->GenPhotonPt[i], fTR->GenPhotonEta[i],fTR->GenPhotonPhi[i],0);
				break;
			} 
		}
		fMT2tree->GenPhoton[0]=photon;
	}
	// add gen Z
	vector<int> indices;
	for(int i=0; i<fTR->NGenLeptons; ++i){
		int ID  =abs(fTR->GenLeptonID[i]);
		int MID =fTR->GenLeptonMID[i];
		if((ID == 11 || ID == 12 || ID == 13 || ID == 14 || ID == 15 || ID == 16) && MID==23){
			indices.push_back(i);
		} 
	}
	if(indices.size()==2){
		TLorentzVector l1(0,0,0,0), l2(0,0,0,0);
		l1.SetPtEtaPhiM(fTR->GenLeptonPt[indices[0]], fTR->GenLeptonEta[indices[0]],fTR->GenLeptonPhi[indices[0]],0);
		l2.SetPtEtaPhiM(fTR->GenLeptonPt[indices[1]], fTR->GenLeptonEta[indices[1]],fTR->GenLeptonPhi[indices[1]],0);
		fMT2tree->GenZ[0]=l1+l2;
	}

	// ----------------------------------------------------------------
	// that was it
	return true;
}

void MT2Analysis::FillMT2treeCalculations(){
	// -----------------------------------------------------------
	// fill MT2hemi 
	
	// hemi 0
	// testmass 0, massless pseudojets, PF-JID, JPt > 20, |jet-eta|<2.4, hemi seed =2 (max inv mass), hemi assoc =3, pf-met, hemi-index 1
	fMT2tree->FillMT2Hemi(0,0,1,20,2.4,2,3,1,0);  
	
	// testmass 0, massless pseudojets, PF-JID, JPt > 20, |jet-eta|<2.4, minimizing Delta_HT, pf-met, hemi-index 1  -> AlphaT version
	// fMT2tree->FillMT2HemiMinDHT(0,0,1,20,2.4,1,2);  
	
	// store MT2 misc variables
	fMT2tree->misc.MT2                 = fMT2tree->hemi[0].MT2;    // note: this is a bit dangerous, 
	fMT2tree->misc.MCT                 = fMT2tree->hemi[0].MCT;

	// other variables to be computed based on MT2tree
	// ----------------------------------------------------------------------------------
	fMT2tree->misc.Vectorsumpt	   = fMT2tree->GetMHTminusMET(1, 20, 2.4, true); // including leptons, ID jets only
	fMT2tree->misc.MinMetJetDPhi       = fMT2tree->MinMetJetDPhi(0,20,5.0,1);
	fMT2tree->misc.MinMetJetDPhi4      = fMT2tree->MinMetJetDPhi(0,20,5.0,1,4);  // use first 4 jets
	fMT2tree->misc.MinMetJetDPhiIndex  = fMT2tree->MinMetJetDPhiIndex(0,20,5.0,1);
	fMT2tree->misc.MinMetBJetDPhi      = fMT2tree->BJetMETdPhi(2,1.74,20,5,1);  // minmetjet dPhi w.r.t SSVHEM bjets with pt > 20 and no eta restriction. 
	fMT2tree->misc.PassJetID           = fMT2tree->PassJetID(50,2.4,1);
	if(fMT2tree->NJets > 0) {
		fMT2tree->misc.Jet0Pass      = (Int_t) fMT2tree->jet[0].IsGoodPFJet(100,2.4,1);
	} else  fMT2tree->misc.Jet0Pass      = 0; 
	if(fMT2tree->NJets > 1) {
		fMT2tree->misc.Jet1Pass      = (Int_t) fMT2tree->jet[1].IsGoodPFJet( 60,2.4,1);
	} else  fMT2tree->misc.Jet1Pass      = 0;
	
	// MHT from jets and taus
	fMT2tree->MHT[0]=fMT2tree->GetMHTlv(1, 20, 2.4, true); // only jets satisfying the loose PF-ID and leptons
	

	// ---------------------------	
	// calo HT and MHT
	float caHT40=0;
	TLorentzVector mht40(0,0,0,0);
	for(int j=0; j<fTR->CANJets; ++j){
	  	if( CAJet(j).Pt()<40 || fabs(CAJet(j).Eta())>3.0 ) continue;
	    	caHT40   += CAJet(j).Pt();
	    	mht40    -= CAJet(j);
	}
	fMT2tree->misc.caloHT40     = caHT40;
	fMT2tree->misc.caloMHT40    = mht40.Pt();
	fMT2tree->misc.caloHT50     = fCaloHT50;
	fMT2tree->misc.caloHT50_ID  = fCaloHT50_ID;
	fMT2tree->misc.caloMHT30    = fCaloMHT30;

	// W and Top decay modes
	fMT2tree->misc.WDecayMode   = fMT2tree->WDecayMode  ();
	fMT2tree->misc.TopDecayMode = fMT2tree->TopDecayMode();

	// _________
	// stuff for Z->nunu (close your eyes...)	
	vector<int>    jindi;
	vector<float>  jpt;
	bool PassJetID_matched(true);
	float HTmatched=0, vectorsumpt_matched_px=0, vectorsumpt_matched_py=0;
	int    NJetsIDLoose_matched=0;
	for(int i=0; i<fTR->PF2PAT3NJets; ++i){
		if(Jet(i).Pt() < 20) continue;  
		bool jet(true);
		for(int gen=0; gen<fMT2tree->NGenLepts; ++gen){
			if( ((abs(fMT2tree->genlept[gen].ID) == 11 || abs(fMT2tree->genlept[gen].ID)==13) && fMT2tree->genlept[gen].MID ==23) ) {
				float deltaR = Util::GetDeltaR(Jet(i).Eta(), fMT2tree->genlept[gen].lv.Eta(), Jet(i).Phi(), fMT2tree->genlept[gen].lv.Phi());
				if(deltaR < 0.4) jet=false;
			}
		}
		if(  jet == false) continue;	
		if(Jet(i).Pt() > 50 && fabs(Jet(i).Eta())<2.4 && IsGoodBasicPFJetPAT3(i,  50., 2.4)==false ){
			PassJetID_matched  = false;
		}
		jindi.push_back(i); jpt.push_back(Jet(i).Pt());
		if(Jet(i).Pt() >50 && fabs(Jet(i).Eta())<3 ){
			HTmatched += Jet(i).Pt();
		}
		if(! IsGoodBasicPFJetPAT3(i,  20., 2.4) ) continue;
		vectorsumpt_matched_px+=Jet(i).Px();
		vectorsumpt_matched_py+=Jet(i).Py();
		NJetsIDLoose_matched++;
	}
	TLorentzVector met     = fMT2tree->pfmet[0];
	for(int gen=0; gen<fMT2tree->NGenLepts; ++gen){
		if( ((abs(fMT2tree->genlept[gen].ID) == 11 || abs(fMT2tree->genlept[gen].ID)==13) && fMT2tree->genlept[gen].MID ==23) ) {
			vectorsumpt_matched_px+=fMT2tree->genlept[gen].lv.Px();
			vectorsumpt_matched_py+=fMT2tree->genlept[gen].lv.Py();
			met +=fMT2tree->genlept[gen].lv;
		}
	}
	
	float           caHT50_matched      =0.0; 
	float           caHT50ID_matched    =0.0; 
	float           caHT50_matchedReco  =0.0;
	float           caHT50ID_matchedReco=0.0;
	TLorentzVector  mht30_matched       (0,0,0,0); 
	TLorentzVector  mht30ID_matched     (0,0,0,0); 
	TLorentzVector  mht30_matchedReco   (0,0,0,0);
	TLorentzVector  mht30ID_matchedReco (0,0,0,0);
	for(int j=0; j<fTR->CANJets; ++j){
	  	if( CAJet(j).Pt()<30 || fabs(CAJet(j).Eta())>3.0 ) continue;
		bool jet(true);
		for(int gen=0; gen<fMT2tree->NGenLepts; ++gen){
			if( ((abs(fMT2tree->genlept[gen].ID) == 11 || abs(fMT2tree->genlept[gen].ID)==13) && fMT2tree->genlept[gen].MID ==23) ) {
				float deltaR = Util::GetDeltaR(CAJet(j).Eta(), fMT2tree->genlept[gen].lv.Eta(), CAJet(j).Phi(), fMT2tree->genlept[gen].lv.Phi());
				if(deltaR < 0.4) jet=false;
			}
		}
		if(  jet == false) continue;	
	  	mht30_matched    -= CAJet(j);  // MHT30
		if(fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001){
		mht30ID_matched  -= CAJet(j); // MHT30_ID
		} 
	  	if(CAJet(j).Pt()>50) {
			caHT50_matched    += CAJet(j).Pt(); //HT50
			if(fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001) {
			caHT50ID_matched  += CAJet(j).Pt(); //HT50
			}
		}
	}
	for(int j=0; j<fTR->CANJets; ++j){
	  	if( CAJet(j).Pt()<30 || fabs(CAJet(j).Eta())>3.0 ) continue;
		bool jet(true);
		// remove overlap from reco electrons and muons
		for(int e=0; e<fMT2tree->NEles; ++e){
			float dR=Util::GetDeltaR(CAJet(j).Eta(), fMT2tree->ele[e].lv.Eta(), CAJet(j).Phi(), fMT2tree->ele[e].lv.Phi());
			if(dR < 0.4) jet=false;
		}
		for(int m=0; m<fMT2tree->NMuons; ++m){
			float dR=Util::GetDeltaR(CAJet(j).Eta(), fMT2tree->muo[m].lv.Eta(), CAJet(j).Phi(), fMT2tree->muo[m].lv.Phi());
			if(dR < 0.4) jet=false;
		}
		if(jet==false) continue;
	  	mht30_matchedReco   -= CAJet(j);  // MHT30
		if(fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001){
		mht30ID_matchedReco -= CAJet(j);  // MHT30_ID
		}
	  	if(CAJet(j).Pt()>50) {
			caHT50_matchedReco    += CAJet(j).Pt(); //HT50
			if(fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001){
			caHT50ID_matchedReco  += CAJet(j).Pt(); //HT50ID
			}
		}
	}


	fMT2tree->Znunu.caloMHT30_matched      =mht30_matched.Pt();
	fMT2tree->Znunu.caloMHT30ID_matched    =mht30ID_matched.Pt();
	fMT2tree->Znunu.caloMHT30_matchedReco  =mht30_matchedReco.Pt();
	fMT2tree->Znunu.caloMHT30ID_matchedReco=mht30ID_matchedReco.Pt();
	fMT2tree->Znunu.caloHT50_matched       =caHT50_matched;
	fMT2tree->Znunu.caloHT50ID_matched     =caHT50ID_matched;
	fMT2tree->Znunu.caloHT50_matchedReco   =caHT50_matchedReco;
	fMT2tree->Znunu.caloHT50ID_matchedReco =caHT50ID_matchedReco;

	float mindPhi=10;
	if(jindi.size()<1){mindPhi = -999.99;}
	else{
		for(int i=0; i<jindi.size(); ++i){
			if(Jet(jindi[i]).Pt()       < 20 ) continue;
			if(fabs(Jet(jindi[i]).Eta())> 5.0) continue;
			float dphi = TMath::Abs(Jet(jindi[i]).DeltaPhi(met));
			if(dphi < mindPhi){ 
				mindPhi = dphi;
			}
		}
	}
	if(mindPhi==10.){
		fMT2tree->Znunu.MinMetplusLeptJetDPhi      = -999.99;
	} else  fMT2tree->Znunu.MinMetplusLeptJetDPhi      =  mindPhi;
	fMT2tree->Znunu.MinMetplusLeptJetDPhiReco          =  fMT2tree->MinMetJetDPhi(0, 20, 5.0 ,3);

	jindi   = Util::VSort(jindi, jpt);
	if(jindi.size() >0){
		fMT2tree->Znunu.Jet0Pass_matched   = (Int_t) IsGoodBasicPFJetPAT3(jindi[0],  100., 2.4);
		fMT2tree->Znunu.LeadingJPt_matched = jpt[0];
	} else  fMT2tree->Znunu.Jet0Pass_matched   =0; 
	if(jindi.size() >1){
		fMT2tree->Znunu.Jet1Pass_matched   = (Int_t) IsGoodBasicPFJetPAT3(jindi[1],   60., 2.4);
		fMT2tree->Znunu.SecondJPt_matched  = jpt[1];
	} else  fMT2tree->Znunu.Jet1Pass_matched   =0;
	fMT2tree->Znunu.PassJetID_matched          = (Int_t) PassJetID_matched;
	fMT2tree->Znunu.Vectorsumpt_matched        = sqrt( pow(vectorsumpt_matched_px+MET().Px(),2) + pow(vectorsumpt_matched_py+MET().Py(),2));
	
	fMT2tree->Znunu.HTmatched                  = HTmatched;
	fMT2tree->Znunu.NJetsIDLoose_matched       = NJetsIDLoose_matched;

	fMT2tree->Znunu.RecoOSee_mll               = fMT2tree->GetDiLeptonInvMass(0, 1, 1, 20, 1); 
	fMT2tree->Znunu.RecoOSmumu_mll             = fMT2tree->GetDiLeptonInvMass(0, 1, 2, 20, 1); 

	fMT2tree->Znunu.GenZee_mll                 = fMT2tree->GenOSDiLeptonInvMass(11,23,0,100);
	fMT2tree->Znunu.GenZee_mll_acc             = fMT2tree->GenOSDiLeptonInvMass(11,23,20,2.4);
	fMT2tree->Znunu.GenZmumu_mll               = fMT2tree->GenOSDiLeptonInvMass(13,23,0,100);
	fMT2tree->Znunu.GenZmumu_mll_acc           = fMT2tree->GenOSDiLeptonInvMass(13,23,20,2.4);

	fMT2tree->Znunu.GenZnunu_e_mll             = fMT2tree->GenOSDiLeptonInvMass(12,23,0,100);
	fMT2tree->Znunu.GenZnunu_e_mll_acc         = fMT2tree->GenOSDiLeptonInvMass(12,23,20,2.4);
	fMT2tree->Znunu.GenZnunu_mu_mll            = fMT2tree->GenOSDiLeptonInvMass(14,23,0,100);
	fMT2tree->Znunu.GenZnunu_mu_mll_acc        = fMT2tree->GenOSDiLeptonInvMass(14,23,20,2.4);
	fMT2tree->Znunu.GenZnunu_tau_mll           = fMT2tree->GenOSDiLeptonInvMass(16,23,0,100);
	fMT2tree->Znunu.GenZnunu_tau_mll_acc       = fMT2tree->GenOSDiLeptonInvMass(16,23,20,2.4);

	fMT2tree->Znunu.METplusLeptsPt             = fMT2tree->GetMETPlusGenLepts(0, 1, 1,  1113, 23, 0, 100, 0, 10000);
	fMT2tree->Znunu.METplusLeptsPtReco         = fMT2tree->GetMETPlusLepts(1);

	//btag SF
	if(!fMT2tree->misc.isData && fbtagFileName.length() !=0){
		bool existing = true;
		if(hceff==0 || hbeff==0 || hleff ==0) existing=false;
		//implementation for >=1 btag only, SSVHPT
		float SFweightErr = 0;
		float SFweight = 1;//outside, since need it there later
		string tagger = "SSVHPT";
		vector<float> jetEff;
		vector<float> jetEffErr;
		vector<float> jetSF;
		vector<float> jetSFErr;
		int njetsusuable = 0;
		for(int n = 0; n<fMT2tree->NJets; ++n){
			if(fMT2tree->jet[n].isPFIDLoose==false) continue;
			if(fMT2tree->jet[n].Flavour<=-7777) continue;//if samples has no flavour information like qcd
			float effPt  = fMT2tree->jet[n].lv.Pt();
			float effEta = fMT2tree->jet[n].lv.Eta();
			++njetsusuable;
			if(abs(fMT2tree->jet[n].Flavour)==5){
				if(existing){
					jetEff.push_back( float(hbeff->GetBinContent(hbeff->FindBin(effPt))) );
					jetEffErr.push_back( float(hbeff->GetBinError(hbeff->FindBin(effPt))) );//here only dummy
				}
				else{
					jetEff.push_back( 0.5 );
					jetEffErr.push_back( 0.05 );//here only dummy
				}
				float SFErr;
				float SF = getBTagSF(SFErr, tagger, effPt, effEta);
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr );//here only dummy
			}
			else if(abs(fMT2tree->jet[n].Flavour)==4){
				if(existing){
					jetEff.push_back( float(hceff->GetBinContent(hceff->FindBin(effPt))) );
					jetEffErr.push_back( float(hceff->GetBinError(hceff->FindBin(effPt))) );//here only dummy
				}
				else{
					jetEff.push_back( 0.08 );
					jetEffErr.push_back( 0.01 );//here only dummy
				}
				float SFErr;
				float SF = getBTagSF(SFErr, tagger, effPt, effEta );
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr*2.);//here only dummy
			}
			else {
				if(existing){
					jetEff.push_back( float(hleff->GetBinContent(hleff->FindBin(effPt))) );
					jetEffErr.push_back( float(hleff->GetBinError(hleff->FindBin(effPt))) );//here only dummy
				}
				else{
					jetEff.push_back( 0.001 );
					jetEffErr.push_back( 0.0003 );//here only dummy
				}
				float SFErr;
				float SF = getMistagSF(SFErr, tagger, effPt, effEta);
				jetSF.push_back(   SF    );
				jetSFErr.push_back(SFErr );//here only dummy
			}
		}
		SFweight = getBTagEventWeightError(SFweightErr, jetEff, jetEffErr, jetSF, jetSFErr, -1);
		if(SFweight==0){
			if(njetsusuable!=0){
				std::cout << "Event has zero weight, set BTagWeight to 0" << std::endl;
				fMT2tree->misc.BTagWeight = SFweight;
			}
			else { //event has no flavour information, use an average event weight
				SFweight = 0.945572;
				SFweightErr = sqrt(0.0257166*0.0257166 + 0.0370919+0.0370919);//here only dummy
			}
		}
		fHistFile->cd();
		fMT2tree->misc.BTagWeight = SFweight;
	} else fMT2tree->misc.BTagWeight = -999.99;
}

// *****************************************************************************
// initialize event 
void MT2Analysis::InitializeEvent(){
	GetLeptonJetIndices();
	FindLeptonConfig();
}

void MT2Analysis::GetLeptonJetIndices(){
	fIsNANObj = false; 

	fElecs.clear();
	fMuons.clear();
	fTaus.clear();
	fJets.clear();
	fPhotons.clear();
	fPhotonJetOverlapRemoved.clear();

	vector<float> mutight;
	for(int i=0; i< fTR->PfMu3NObjs; ++i){
		if(std::isnan(fTR->PfMu3Pt[i]))  {fIsNANObj = true; continue;} //protection against objects with NAN-Pt
		fMuons.push_back(i);
		mutight.push_back(fTR->PfMu3Pt[i]);
	}
	fMuons      = Util::VSort(fMuons     , mutight);
	
	vector<float> eltight;
	for(int i=0; i< fTR->PfEl3NObjs; ++i){
		if(std::isnan(fTR->PfEl3Pt[i]))  {fIsNANObj = true; continue;} //protection against objects with NAN-Pt
		fElecs.push_back(i);
		eltight.push_back(fTR->PfEl3Pt[i]);
	}
	fElecs      = Util::VSort(fElecs     , eltight);
	
	vector<float> pt1; 
	for(int ij=0; ij < fTR->PF2PAT3NJets; ++ij){
		if(std::isnan(fTR->PF2PAT3JPt[ij])){ fIsNANObj = true; continue; } //protection against objects with NAN-Pt
		if(fTR->PF2PAT3JScale[ij]<0) continue;  // ignore jets with negative JEC
		if(Jet(ij).Pt() < 20)        continue;  // note: ETH ntuple only stores PF2PAT3Jets > 15 GeV (defualt config)
		fJets.push_back(ij);                    // fJets has all jets except for duplicates with selected leptons
		pt1.push_back(Jet(ij).Pt());
	}
	fJets        = Util::VSort(fJets,       pt1);

	// Warning: taus are also contained in the jet-collection.
	vector<float> taus;
	for(int i=0; i< fTR->PfTau3NObjs; ++i){
		if(std::isnan(fTR->PfTau3Pt[i])) { fIsNANObj = true; continue;} //protection against objects with NAN-Pt
		if(fTR->PfTau3Pt[i]   < 20    ) continue; // note: taus go up to 2.5 in Eta
		fTaus.push_back(i);
		taus.push_back(fTR->PfTau3Pt[i]);
	}
	fTaus          = Util::VSort(fTaus     , taus);

	// Photons ---------------------------------------------------------------------------------
	if(fTR->fChain->GetBranch("NPhotons")     ==NULL){return;} 	     // protection for TESCO ntuples without photon info
	vector<float> photon_pts;
	for(int i=0; i<fTR->NPhotons; ++i){
		if(std::isnan(fTR->PhoPt[i]))  {fIsNANObj = true; continue;} // protection against objects with NAN-Pt	
		if(! IsGoodPhoton(i))                             continue; 
		if(! IsGoodPhotonEGMLooseRelISO(i))               continue;   // preselection: use only photons passing the loose RelIso
		fPhotons.push_back(i);
		photon_pts.push_back(fTR->PhoPt[i]);
		fPhotonJetOverlapRemoved.push_back(false);
	}
	fPhotons     = Util::VSort(fPhotons, photon_pts);

	// remove jet or ele matched to photon ----------------------------------------------------
	if(! fRemovePhoton)   return;

	vector<int> removeEleIndices;
	vector<int> removeJetIndices;
	
	if(fRemovePhoton){
		// find jet or ele index to be removed
		for(int iphoton=0; iphoton<fPhotons.size(); ++iphoton){
			float minDR=10; int index=-1; bool jetmatch(false);
			TLorentzVector phot(0,0,0,0);
			phot.SetPtEtaPhiM(fTR->PhoPt[fPhotons[iphoton]], fTR->PhoEta[fPhotons[iphoton]], fTR->PhoPhi[fPhotons[iphoton]], 0.);
			for(int el=0; el<fElecs.size(); ++el){
				TLorentzVector ele(0,0,0,0);
				ele.SetPtEtaPhiM(fTR->PfEl3Pt[fElecs[el]], fTR->PfEl3Eta[fElecs[el]], fTR->PfEl3Phi[fElecs[el]], 0.);
				float dR = phot.DeltaR(ele);
				if (dR< minDR) {minDR= dR; index =fElecs[el];} 
			}
			for(int ijet=0; ijet<fJets.size(); ++ijet){
				TLorentzVector jet(0,0,0,0);
				jet.SetPtEtaPhiE(fTR->PF2PAT3JPt[fJets[ijet]], fTR->PF2PAT3JEta[fJets[ijet]], fTR->PF2PAT3JPhi[fJets[ijet]], fTR->PF2PAT3JE[fJets[ijet]]);
				float dR = phot.DeltaR(jet);
				if (dR< minDR) {minDR= dR; index =fJets[ijet]; jetmatch=true;} 
			}
			if(minDR<0.2){ // require match within dR=0.2
				if(jetmatch){removeJetIndices.push_back(index);}
				else        {removeEleIndices.push_back(index);}
				fPhotonJetOverlapRemoved[iphoton]=true; // set mark that for photon "iphoton", the corresponding jet was found (will be) removed 
			}	
		}
	}
	if(fRemovePhoton && fPhotons.size()>0){	
		for(int i=0; i<removeEleIndices.size(); ++i){
			for(int j=0; j<fElecs.size(); ++j){
				if(fElecs[j]==removeEleIndices[i]) {
					fElecs    .erase(fElecs.begin()   +j);
					eltight   .erase(eltight.begin()  +j);
				}
			}
		}
		for(int i=0; i<removeJetIndices.size(); ++i){
			for(int j=0; j<fJets.size(); ++j){
				if(fJets[j]==removeJetIndices[i]) {
					fJets .erase(fJets.begin()   +j);
					pt1   .erase(pt1.begin()     +j);
				}
			}
		}
	}
	// -----------
}

void MT2Analysis::FindLeptonConfig(){
	fLeptConfig = null;
	if(fElecs.size()==1 && fMuons.size()==0){
		fLeptConfig=e;
	}else if(fElecs.size()==0 && fMuons.size()==1) {
		fLeptConfig=mu;
	}else if(fElecs.size() + fMuons.size() ==2){
		int charge1, charge2;
		if(fElecs.size()==2){
			charge1=fTR->PfEl3Charge[fElecs[0]];
			charge2=fTR->PfEl3Charge[fElecs[1]];
			if(charge1*charge2==1){fLeptConfig=SS_ee;}
			else{fLeptConfig=OS_ee;}
		} else if(fMuons.size()==2){
			charge1=fTR->PfMu3Charge[fMuons[0]];
			charge2=fTR->PfMu3Charge[fMuons[1]];
			if(charge1*charge2==1){fLeptConfig=SS_mumu;}
			else{fLeptConfig=OS_mumu;}
		} else{
			charge1=fTR->PfEl3Charge[fElecs[0]];
			charge2=fTR->PfMu3Charge[fMuons[0]];			
			if(charge1*charge2==1){fLeptConfig=SS_emu;}
			else{fLeptConfig=OS_emu;}
		}
	}else if(fElecs.size() + fMuons.size() >2){
		fLeptConfig=multilept;
	}

}

// *****************************************************************************
// event selector
bool MT2Analysis::IsSelectedEvent(){
	// goodevt from UserAnalysisBase
	if(!IsGoodEvent()) {return false;}


	// Run
	if(fTR->Run < fCut_Run_min ) {return false;}
	if(fTR->Run > fCut_Run_max ) {return false;}

	// Protection against events with NAN-Pt objects
	if(fIsNANObj) {
		cout << "WARNING: Event " << fTR->Event << " Lumi " << fTR->LumiSection << " Run " << fTR->Run << " has NAN-Pt Objects!! Skipping Event" << endl;
		return false;
	}
	
	//PtHat
	if(fTR->PtHat > fCut_PtHat_max ){return false;}
	
	// MET
	if(MET().Pt() < fCut_PFMET_min){return false;}
	
	// HLT triggers
	if(fRequiredHLT.size() !=0 ){
		bool HLT_fire(false);
		for(int i=0; i<fRequiredHLT.size(); ++i){
			if( GetHLTResult(fRequiredHLT[i]) ){  // require HLT bit
				HLT_fire = true;
			} 
		}
		if(! HLT_fire) return false;
	}
	if(fVetoedHLT.size() !=0){
		for(int i=0; i<fVetoedHLT.size(); ++i){
			if( GetHLTResult(fVetoedHLT[i]) ){   // veto HLT bit 
				return false;
			} 
		}
	}

	// HT from jets + taus
	float HT=0;
	for(int j=0; j<fJets.size(); ++j){
	  if(Jet(fJets[j]).Pt() > 50 && fabs(Jet(fJets[j]).Eta())<3.0){
	    HT += Jet(fJets[j]).Pt();
	  }
	}
	fHT = HT;
	if(HT<fCut_HT_min){return false;}
	
	//caloHT and calo MHT
	TLorentzVector mht30(0,0,0,0);
	fCaloHT50 =0.0, fCaloHT50_ID =0.0;
	for(int j=0; j<fTR->CANJets; ++j){
		if( CAJet(j).Pt()<30 || fabs(CAJet(j).Eta())>3.0 ) continue;
	  	// MHT
		mht30 -= CAJet(j);
		// HT
		if( CAJet(j).Pt()<50 ) continue;
		fCaloHT50  += CAJet(j).Pt();
		// HT_ID
		if(fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001) fCaloHT50_ID += CAJet(j).Pt();
	}
	fCaloMHT30   =mht30.Pt();
	if(fCaloHT50      < fCut_caloHT50_min   ) return false;
	if(fCaloHT50_ID   < fCut_caloHT50ID_min ) return false;
	if(fCaloMHT30     < fCut_caloMHT30_min  ) return false;

	// leading jets including JID for jets
	bool leadingjets(true);
	if(fCut_JPt_hardest_min > 0){
		if(fJets.size() <1) leadingjets=false;
		if(! IsGoodBasicPFJetPAT3(fJets[0], fCut_JPt_hardest_min, 2.4)) {leadingjets=false;}
	}
	if(fCut_JPt_second_min > 0){
		if(fJets.size() <2) leadingjets=false;
		if(! IsGoodBasicPFJetPAT3(fJets[1], fCut_JPt_second_min, 2.4)) {leadingjets=false;}
	}
	if(leadingjets == false) return false;
	
	// flag crazy events (dedicated to filter HCAL "laser" events)
	if(fTR->HCALSumEt > 10000 && fTR->PF2PAT3NJets > 25 ){
	  fCrazyHCAL =1;
		cout << "WARNING: crazy HCAL event: Run " << fTR->Run << " LS " << fTR->LumiSection  << " Event " << fTR->Event 
		     << " has HCALSumET " << fTR->HCALSumEt << " njets " << fTR->PF2PAT3NJets << endl;
	}
	else    fCrazyHCAL =0;	

	// flag events with Jets with correction factor from JE correction 
	fNegativeJEC =0; 
	for( int i=0; i<fTR->PF2PAT3NJets; ++i){
		if(fTR->PF2PAT3JScale[i]>=0) continue;
		fNegativeJEC =1;
	}
	if(fNegativeJEC && fVerbose >1) {
		cout << "WARNING: Event with Jets with negative JEC: Run " << fTR->Run << " LS " << fTR->LumiSection << " Event " << fTR->Event
		     << " N PFJets " << fTR->PF2PAT3NJets << " NCaloJets " << fTR->CANJets << endl; 
	}

	// ------------------------------------------------------------------------------------------	
	return true;	
}


// ----------------------------------------------------------------------------
// Parsing and reading of cuts
void MT2Analysis::ReadCuts(const char* SetofCuts="MT2_cuts/default.dat"){
	
	ifstream IN(SetofCuts);
	char buffer[200];
	char ParName[100];
	char StringValue[100];
	float ParValue;
	int   FlagValue;
	int   IntValue;

	bool verbose(true);
	bool ok(false);
	
	
	while( IN.getline(buffer, 200, '\n') ){
		ok = false;
		if (buffer[0] == '#') {continue;} // Skip lines commented with '#'

		// strings
		sscanf(buffer, "%s %s", ParName, StringValue);
		if( !strcmp(ParName, "SetName") ){
			fSetName = TString(StringValue); ok = true;
			if(verbose){cout << "Reading cut parameters for set: " << fSetName << endl; }
		} else if( !strcmp(ParName, "HLT_required") ){
			fRequiredHLT.push_back(StringValue); ok = true;
		} else if( !strcmp(ParName, "HLT_vetoed") ){
			fVetoedHLT.push_back(StringValue); ok = true;
		} else if( !strcmp(ParName, "ECALBEfile") ){
			fBEfiles.push_back(StringValue); ok = true;
		} else if( !strcmp(ParName, "ECALTPfile") ){
			fTPfiles.push_back(StringValue); ok = true;
		}	

		// ints
		sscanf(buffer, "%s %i", ParName, &IntValue);
		if( !strcmp(ParName, "Run_min") ){
			fCut_Run_min = int(IntValue); ok = true;
		} else if( !strcmp(ParName, "Run_max") ){
			fCut_Run_max = int(IntValue); ok = true;
		} else if( !strcmp(ParName, "JESUpDown")){
			fJESUpDown   = int(IntValue); ok = true;
		}

		// floats 
		sscanf(buffer, "%s %f", ParName, &ParValue);
		if( !strcmp(ParName, "PFMET_min") ){
			fCut_PFMET_min            = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "HT_min") ){
			fCut_HT_min               = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "caloHT50_min") ){
			fCut_caloHT50_min         = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "caloHT50ID_min") ){
			fCut_caloHT50ID_min       = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "caloMHT30_min") ){
			fCut_caloMHT30_min        = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "caloMHT30ID_min") ){
			fCut_caloMHT30ID_min      = float(ParValue); ok = true;
		} else if( !strcmp(ParName, "JPt_hardest_min") ){
			fCut_JPt_hardest_min      = float(ParValue); ok = true;			
		} else if( !strcmp(ParName, "JPt_second_min") ){
			fCut_JPt_second_min      = float(ParValue); ok = true;			
		} else if( !strcmp(ParName, "PtHat_max")){
			fCut_PtHat_max            = float(ParValue); ok = true;
		} 

		if(!ok) cout << "%% MT2Analysis::ReadCuts ==> ERROR: Unknown variable " << ParName << endl;
	}	
	if(verbose){
		cout << "setting cuts to: " << endl;
		cout << "  PFMET_min                   " << fCut_PFMET_min                  <<endl;
		cout << "  HT_min                      " << fCut_HT_min                     <<endl;
		cout << "  caloHT50_min                " << fCut_caloHT50_min               <<endl;
		cout << "  caloHT50ID_min              " << fCut_caloHT50ID_min             <<endl;
		cout << "  caloMHT30_min               " << fCut_caloMHT30_min              <<endl;
		cout << "  caloMHT30ID_min             " << fCut_caloMHT30ID_min            <<endl;
		cout << "  JPt_hardest_min             " << fCut_JPt_hardest_min            <<endl;
		cout << "  JPt_second_min              " << fCut_JPt_second_min             <<endl;
		cout << "  PtHat_max                   " << fCut_PtHat_max                  <<endl;
		cout << "  Run_min                     " << fCut_Run_min                    <<endl;
		cout << "  Run_max                     " << fCut_Run_max                    <<endl;
		cout << "  JESUpDown                   " << fJESUpDown                      <<endl;
		if(fJESUpDown ==-1 || fJESUpDown==1){
		if(fJESUpDown ==1) cout << "  Do up  -scaling of Jets and MET "             <<endl;
		else               cout << "  Do down-scaling of Jets and MET "             <<endl;
		fDoJESUncertainty = true;
		}

		for(int i=0; i<fRequiredHLT.size(); ++i){
			cout << "  HLTRequired (logic OR)      " << fRequiredHLT[i]                  <<endl;
		}
		for(int i=0; i<fVetoedHLT.size(); ++i){
			cout << "  HLTVetoed                   " << fVetoedHLT[i]                    <<endl;
		}
		for(int i=0; i<fBEfiles.size(); ++i){
			cout << "  ECALBEfile                  " << fBEfiles[i]                      << endl;
		}
		for(int i=0; i<fTPfiles.size(); ++i){
			cout << "  ECALTPfile                  " << fTPfiles[i]                      << endl;
		}
		cout << "--------------"    << endl;	
	}			
}

// BE file parser
void MT2Analysis::DeadCellParser(DeadCellFilter &DeadCellFilter_, string file_){
	string line;	
	string path="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/ECALDeadCell/";
	string file=path+file_;
	ifstream IN(file.c_str());
	if (!IN.is_open()) {cout << "ERROR: cannot open dead cell file " << file << endl; exit(1);}
	else{
		if(fVerbose>0) cout << "--------------------------"          << endl;
		if(fVerbose>0) cout << "DeadCellParser: read file " << file  << endl;
		while ( ! IN.eof() ){
			getline (IN, line);
			TString Line=  line.c_str();
			if(Line.EndsWith("\",") && Line.BeginsWith("\"")){
				Line.ReplaceAll("\"", "");
				Line.ReplaceAll(",", "");
				TObjArray *p= (TObjArray*) Line.Tokenize(":");
				TString run  =((TObjString*)p->At(0))->GetString();
				TString lumi =((TObjString*)p->At(1))->GetString();
				TString event=((TObjString*)p->At(2))->GetString();
				DeadCellFilter_.run.  push_back(run.Atoi());
				DeadCellFilter_.lumi. push_back(lumi.Atoi());
				DeadCellFilter_.event.push_back(event.Atoi());
			}
		}
		if(fVerbose >0) cout << "DeadCellParser: read " <<  DeadCellFilter_.run.size() << " events to be filtered. " << endl; 
		if(fVerbose >0) cout << "--------------------------"          << endl;
	}
}

//****************************************************************************************************
// Photon Selectors

bool MT2Analysis::IsGoodPhotonEGMLooseID(int i){
	// EGM-10-006 Loose Photon selection
	bool isGood(true);
	if(!IsGoodPhoton(i)                                                    ) isGood=false;
	if( fabs(fTR->PhoEta[i]) < 1.4442 && fTR->PhoSigmaIetaIeta[i] > 0.01   ) isGood=false;
	if( fabs(fTR->PhoEta[i]) > 1.566  && fTR->PhoSigmaIetaIeta[i] > 0.03   ) isGood=false;
	if( fTR->PhoHoverE[i]   > 0.05                                         ) isGood=false; // H/E cut
	return isGood;
}
bool MT2Analysis::IsGoodPhotonEGMLooseISO(int i){
	// EGM-10-006 Loose Photon selection
	bool isGood(true);
	if(!IsGoodPhoton(i)                                                    ) isGood=false;
	if( fTR->PhoIso04TrkHollow[i] > 2.0                                    ) isGood=false; // trk Iso
	if( fTR->PhoIso04Ecal[i] > 4.2                                         ) isGood=false; // Ecal ISo
	if( fTR->PhoIso04Hcal[i] > 2.2                                         ) isGood=false; // Hcal Iso
	return isGood;
}
bool MT2Analysis::IsGoodPhotonEGMLooseRelISO(int i){
	// e.g. CMS AN AN-2011/033
	bool isGood(true);
	float pt = fTR->PhoPt[i];
	if(!IsGoodPhoton(i)                                                     ) isGood=false;
	if( fTR->PhoIso04TrkHollow[i] > (2.0 + 0.001*pt)                        ) isGood=false; // trk Iso
	if( fTR->PhoIso04Ecal[i]      > (4.2 + 0.003*pt)                        ) isGood=false; // Ecal ISo
	if( fTR->PhoIso04Hcal[i]      > (2.2 + 0.001*pt)                        ) isGood=false; // Hcal Iso
	return isGood;
}
bool MT2Analysis::IsGoodPhotonEGMTightID(int i){
	// EGM-10-006 Loose Photon selection
	bool isGood(true);
	if(!IsGoodPhoton(i)                                                    ) isGood=false;
	if( fabs(fTR->PhoEta[i]) < 1.4442 && fTR->PhoSigmaIetaIeta[i] > 0.01   ) isGood=false;
	if( fabs(fTR->PhoEta[i]) > 1.566  && fTR->PhoSigmaIetaIeta[i] > 0.028  ) isGood=false;
	if( fTR->PhoHoverE[i]   > 0.03                                         ) isGood=false; // H/E cut
	return isGood;
}
bool MT2Analysis::IsGoodPhotonEGMTightISO(int i){
	// EGM-10-006 Loose Photon selection
	bool isGood(true);
	if(!IsGoodPhoton(i)                                                    ) isGood=false;
	if( fTR->PhoIso04TrkHollow[i] > 0.9                                    ) isGood=false; // trk Iso
	if( fTR->PhoIso04Ecal[i] > 2.4                                         ) isGood=false; // Ecal ISo
	if( fTR->PhoIso04Hcal[i] > 1.0                                         ) isGood=false; // Hcal Iso
	return isGood;
}

bool MT2Analysis::IsGoodPhoton(int i){
	if( fTR->fChain->GetBranch("PhoPt")             ==NULL   ) {cout << "WARNING: Photon variable not found" << endl; return false; }
	if( fTR->fChain->GetBranch("PhoEta")            ==NULL   ) {cout << "WARNING: Photon variable not found" << endl; return false; } 
	if( fTR->fChain->GetBranch("PhoEta")            ==NULL   ) {cout << "WARNING: Photon variable not found" << endl; return false; }  
	if( fTR->fChain->GetBranch("PhoHoverE")         ==NULL   ) {cout << "WARNING: Photon variable not found" << endl; return false; }  
	if( fTR->fChain->GetBranch("PhoSigmaIetaIeta")  ==NULL   ) {cout << "WARNING: Photon variable not found" << endl; return false; }  
	if( fTR->fChain->GetBranch("PhoIso04TrkHollow") ==NULL   ) {cout << "WARNING: Photon variable not found" << endl; return false; }  
	if( fTR->fChain->GetBranch("PhoIso04Ecal")      ==NULL   ) {cout << "WARNING: Photon variable not found" << endl; return false; }  
	if( fTR->fChain->GetBranch("PhoIso04Hcal")      ==NULL   ) {cout << "WARNING: Photon variable not found" << endl; return false; }  
	if( fTR->fChain->GetBranch("PhoHasPixSeed")     ==NULL   ) {cout << "WARNING: Photon variable not found" << endl; return false; }  

	if( fTR->PhoPt[i] < 20                                                 ) return false; // pt cut
	if( fabs(fTR->PhoEta[i])> 2.4                                          ) return false;
	if( fabs(fTR->PhoEta[i])> 1.4442 && fabs(fTR->PhoEta[i])<1.566         ) return false; // veto EB-EE gap
	if( fTR->PhoHoverE[i]   > 0.05                                         ) return false; // H/E cut
	if( fTR->PhoHasPixSeed[i]==1                                           ) return false; // veto pixel seed for electron rejection 
	return true;
}

// ****************************************************************************************************
// Jet Selectors

bool MT2Analysis::IsGoodBasicPFJetPAT3(int index, float ptcut, float absetacut){
	// Basic PF jet cleaning and ID cuts
	// cut at pt of ptcut (default = 30 GeV)
	// cut at abs(eta) of absetacut (default = 2.5)
	if(Jet(index).Pt() < ptcut                  ) return false;
	if(fabs(fTR->PF2PAT3JEta[index]) > absetacut) return false;
	if(fTR->PF2PAT3JIDLoose[index]    ==0       ) return false;
	if(fTR->PF2PAT3JScale[index]     < 0        ) return false;
	return true;
}

bool MT2Analysis::IsGoodPFJetMediumPAT3(int index, float ptcut, float absetacut) {
	// Medium PF JID
	if ( ! IsGoodBasicPFJetPAT3(index, ptcut, absetacut)  ) return false;
	if ( !(fTR->PF2PAT3JNeuHadfrac[index] < 0.95)         ) return false;
	if ( !(fTR->PF2PAT3JNeuEmfrac[index]  < 0.95)         ) return false;
	return true;
}

bool MT2Analysis::IsGoodPFJetTightPAT3(int index, float ptcut, float absetacut) {
	// Tight PF JID
	if ( ! IsGoodBasicPFJetPAT3(index, ptcut, absetacut)  ) return false;
	if ( !(fTR->PF2PAT3JNeuHadfrac[index] < 0.90)         ) return false;
	if ( !(fTR->PF2PAT3JNeuEmfrac[index]  < 0.90)         ) return false;
	return true;
}



// Jets and JES uncertainty
void MT2Analysis::Initialize_JetCorrectionUncertainty(){
	string Calo=fJEC+"/AK5PF_Uncertainty.txt";
	string PF  =fJEC+"/AK5Calo_Uncertainty.txt";

	ifstream fileCalo(Calo.c_str());
	ifstream filePF  (PF.c_str());
	if(! filePF   ) {cout << "ERROR: cannot open file " << PF     << endl; exit(1); }
	else            {cout << "  using file "            << PF     << endl;          }
	if(! fileCalo ) {cout << "ERROR: cannot open file " << Calo   << endl; exit(1); }
	else            {cout << "  using file "            << Calo   << endl;          } 

	fJecUncCalo = new JetCorrectionUncertainty(Calo);
	fJecUncPF   = new JetCorrectionUncertainty(PF); 
}

void MT2Analysis::Initialize_JetEnergyCorrection(){
	string Calo_L2  =fJEC+"/AK5Calo_L2Relative.txt";
	string Calo_L3  =fJEC+"/AK5Calo_L3Absolute.txt";
	string Calo_RES =fJEC+"/AK5Calo_L2L3Residual.txt";
	string PF_L1    =fJEC+"/AK5PF_L1FastJet.txt";
	string PF_L2    =fJEC+"/AK5PF_L2Relative.txt";
	string PF_L3    =fJEC+"/AK5PF_L3Absolute.txt";
	string PF_RES   =fJEC+"/AK5PF_L2L3Residual.txt";

	ifstream fileCaloL2  (Calo_L2.c_str());
	ifstream fileCaloL3  (Calo_L3.c_str());
	ifstream fileCaloRES (Calo_RES.c_str());
	ifstream filePFL1    (PF_L1.c_str());
	ifstream filePFL2    (PF_L2.c_str());
	ifstream filePFL3    (PF_L3.c_str());
	ifstream filePFRES   (PF_RES.c_str());
	if(! filePFL1   ) {cout << "ERROR: cannot open file " << PF_L1     << endl; exit(1); }
	else              {cout << "  using file "            << PF_L1     << endl;          }
	if(! filePFL2   ) {cout << "ERROR: cannot open file " << PF_L2     << endl; exit(1); }
	else              {cout << "  using file "            << PF_L2     << endl;          }
	if(! filePFL3   ) {cout << "ERROR: cannot open file " << PF_L3     << endl; exit(1); }
	else              {cout << "  using file "            << PF_L3     << endl;          }
	if(! filePFRES  ) {cout << "ERROR: cannot open file " << PF_RES    << endl; exit(1); }
	else              {cout << "  using file "            << PF_RES    << endl;          }
	if(! fileCaloL2 ) {cout << "ERROR: cannot open file " << Calo_L2   << endl; exit(1); }
	else              {cout << "  using file "            << Calo_L2   << endl;          } 
	if(! fileCaloL3 ) {cout << "ERROR: cannot open file " << Calo_L3   << endl; exit(1); }
	else              {cout << "  using file "            << Calo_L3   << endl;          } 
	if(! fileCaloRES) {cout << "ERROR: cannot open file " << Calo_RES   << endl; exit(1); }
	else              {cout << "  using file "            << Calo_RES   << endl;          } 

	fJecCaloL2  = new JetCorrectorParameters(Calo_L2);
	fJecCaloL3  = new JetCorrectorParameters(Calo_L3);
	fJecCaloRES = new JetCorrectorParameters(Calo_RES);
	fJecPFL1    = new JetCorrectorParameters(PF_L1); 
	fJecPFL2    = new JetCorrectorParameters(PF_L2); 
	fJecPFL3    = new JetCorrectorParameters(PF_L3); 
	fJecPFRES   = new JetCorrectorParameters(PF_RES); 

	vector<JetCorrectorParameters> JecPF;
	vector<JetCorrectorParameters> JecCalo;
	// ATTENTION: order matters here!
	JecPF.push_back(*fJecPFL1);     JecPF.push_back(*fJecPFL2);     JecPF.push_back(*fJecPFL3);     JecPF.push_back(*fJecPFRES);
	JecCalo.push_back(*fJecCaloL2); JecCalo.push_back(*fJecCaloL3); JecCalo.push_back(*fJecCaloRES); 
	
	fJetCorrectorPF   = new FactorizedJetCorrector(JecPF);
	fJetCorrectorCalo = new FactorizedJetCorrector(JecCalo);
}

// pf
TLorentzVector MT2Analysis::Jet(int index){
	TLorentzVector j(0,0,0,0);
	j.SetPtEtaPhiE(fTR->PF2PAT3JPt[index], fTR->PF2PAT3JEta[index], fTR->PF2PAT3JPhi[index], fTR->PF2PAT3JE[index]);
	if      (fJEC.length()==0) return j;
	else                       return MT2Analysis::PFJetScaled(j, fTR->PF2PAT3JScale[index], fTR->PF2PAT3JArea[index], fTR->Rho);
}

TLorentzVector MT2Analysis::PFJetScaled(TLorentzVector j, float old_scale, float area, float rho){
	// get new correction factor from DB dumpled txt files
	float new_scale = GetPFJEC(j.Pt(), old_scale, j.Eta(), area, rho, (fisData)?1230:123); // get new scale factor: L1FastL2L3 for MC, L1fastL2L3+RES for data

	TLorentzVector j_scaled(0.,0.,0.,0);
	j_scaled.SetPtEtaPhiM(j.Perp()*new_scale/old_scale, j.Eta(), j.Phi(), j.M());

	// optionally up or downscale jet: JES uncertainty is a function of corrected Pt
	if(fJESUpDown== 1) j_scaled.SetPtEtaPhiM(j_scaled.Perp()*(1+GetJECUncertPF(j_scaled.Pt(),j_scaled.Eta())), j_scaled.Eta(), j_scaled.Phi(), j_scaled.M());
	if(fJESUpDown==-1) j_scaled.SetPtEtaPhiM(j_scaled.Perp()*(1-GetJECUncertPF(j_scaled.Pt(),j_scaled.Eta())), j_scaled.Eta(), j_scaled.Phi(), j_scaled.M());
	return j_scaled;
}

float MT2Analysis::GetJECUncertPF(float pt, float eta){
	// pt must be corrected! not raw  	

	// eta cannot be greater than 5! 
	if      (eta> 5.0) eta = 5.0;
	else if (eta<-5.0) eta =-5.0;

	fJecUncPF->setJetPt(pt);   
	fJecUncPF->setJetEta(eta); 
	float uncert= fJecUncPF->getUncertainty(true);
	return uncert; 
}

float MT2Analysis::GetPFJEC(float corrpt, float scale, float eta, float area, float rho, int level){

	fJetCorrectorPF->setJetPt(corrpt/scale); // WARNING: JEC is a function of RAW pt
	fJetCorrectorPF->setJetEta(eta);
	fJetCorrectorPF->setJetA(area);
	fJetCorrectorPF->setRho(rho);
	
	vector<float> factors = fJetCorrectorPF->getSubCorrections();
	// convention for correction factors as follows: 0->L1, 1->L1L2, 2->L1L2L3, 3->L1L2L3Res
	// see MT2Analysis::Initialize_JetEnergyCorrection()
	
	if(!fisData && level==1230){
		cout << "ERROR: residual corrections are only applied on data, this is MC!"<<endl;
		exit(-1);
	}

	float scalefactor=-1;
	if     (level==1)    scalefactor=factors[0];
	else if(level==12)   scalefactor=factors[1];
	else if(level==123)  scalefactor=factors[2];
	else if(level==1230) scalefactor=factors[3];
	else {
		cout << "ERROR: convention for correction factors as follows: 0->L1, 1->L1L2, 2->L1L2L3, 3->L1L2L3Res" << endl;
		exit(-1);	
	}
	return scalefactor;
}

// Calo
TLorentzVector MT2Analysis::CAJet(int index){
	TLorentzVector j(0,0,0,0);
	j.SetPtEtaPhiE(fTR->CAJPt[index], fTR->CAJEta[index], fTR->CAJPhi[index], fTR->CAJE[index]);
	if  (fJEC.length()==0)    return j;
	else                      return MT2Analysis::CAJetScaled(j, fTR->CAJScale[index]);
}


TLorentzVector MT2Analysis::CAJetScaled(TLorentzVector j, float old_scale){
	// get new correction factor from DB dumpled txt files
	float new_scale = GetCaloJEC(j.Pt(), old_scale, j.Eta(), (fisData)?230:23); // get new scale factor: L2L3 for MC, L2L3+RES for data

	TLorentzVector j_scaled(0.,0.,0.,0);
	j_scaled.SetPtEtaPhiM(j.Perp()*new_scale/old_scale, j.Eta(), j.Phi(), j.M());

	// optionally up or downscale jet: JES uncertainty is a function of corrected Pt
	if(fJESUpDown== 1) j_scaled.SetPtEtaPhiM(j_scaled.Perp()*(1+GetJECUncertCalo(j_scaled.Pt(),j_scaled.Eta())), j_scaled.Eta(), j_scaled.Phi(), j_scaled.M());
	if(fJESUpDown==-1) j_scaled.SetPtEtaPhiM(j_scaled.Perp()*(1-GetJECUncertCalo(j_scaled.Pt(),j_scaled.Eta())), j_scaled.Eta(), j_scaled.Phi(), j_scaled.M());
	return j_scaled;
}

float MT2Analysis::GetJECUncertCalo(float pt, float eta){
	// pt must be corrected! not raw  	
	
	// eta cannot be greater than 5! 
	if      (eta> 5.0) eta = 5.0;
	else if (eta<-5.0) eta =-5.0;
	fJecUncCalo->setJetPt(pt);
	fJecUncCalo->setJetEta(eta);
	float uncert= fJecUncCalo->getUncertainty(true);
	return uncert; 
}

float MT2Analysis::GetCaloJEC(float corrpt, float scale, float eta, int level){
	fJetCorrectorCalo->setJetEta(eta);
	fJetCorrectorCalo->setJetPt(corrpt/scale); // WARNING: JEC is a function of RAW pt
	
	vector<float> factors = fJetCorrectorCalo->getSubCorrections();
	// convention for correction factors as follows: 0->L2, 1->L2L3, 2->L2L3+RES
	// see MT2Analysis::Initialize_JetEnergyCorrection()
	
	if     (level==2)     return factors[0];
	else if(level==23)    return factors[1];
	else if(level==230)   return factors[2];
	else {
		cout << "ERROR: convention for correction factors as follows: 0->L2, 1->L2L3, 2->L2L3+RES" << endl;
		exit(-1);	
	}
}

// pfMET 
TLorentzVector MT2Analysis::MET(){
	TLorentzVector MET(0,0,0,0);
	if(!fDoJESUncertainty){
		if(fRemovePhoton && fPhotons.size()==1){
			TLorentzVector trueMET(0,0,0,0), photon(0,0,0,0);
			trueMET.SetPtEtaPhiM(fTR->PFMETPAT, 0., fTR->PFMETPATphi, 0);
			photon .SetPtEtaPhiM(fTR->PhoPt[fPhotons[0]], 0., fTR->PhoPhi[fPhotons[0]],   0);
			MET    = photon + trueMET; 
		}else{
			MET.SetPtEtaPhiM(fTR->PFMETPAT, 0., fTR->PFMETPATphi, 0);
		}
		return MET;
	} else{
		if      (fJESUpDown==1){
    			MET.SetPtEtaPhiM(fTR->PFMETPAT*1.05, 0., fTR->PFMETPATphi, 0);
		}else if(fJESUpDown==-1){
			MET.SetPtEtaPhiM(fTR->PFMETPAT*0.95, 0., fTR->PFMETPATphi, 0);
		}else{
			cout << " something wrong in met scaling" << endl; exit(1);
		}		
      	return MET;
	}
}

