#include "helper/Utilities.hh"
#include "MT2Analysis.hh"
#include "TLorentzVector.h"
#include "MT2tree.hh"
#include <sstream>


using namespace std;

MT2Analysis::MT2Analysis(TreeReader *tr) : UserAnalysisBase(tr){
	Util::SetStyle();	
	fCut_PFMET_min                      = 0;
	fCut_HT_min                         = 0;
	fCut_JPt_hardest_min                = 0;
	fCut_JPt_second_min                 = 0;
        fCut_PtHat_max                      = 999999.;
	fCut_Run_min                        = 0;
	fCut_Run_max                        = 9999999;
	fDoJESUncertainty                   = false;
	fJESUpDown                          = 0;
	fCut_NJets40_min                    = 0;
	fCut_NLeptons_min                   = 0;

	fRemovePhoton                       = 0;
	fID                                 = -1;
	fbtagFileName                       = "";

	fisType1MET                         = false;
	fisCHSJets                          = false;
	fisFastSim                          = false;

	fRequiredHLT.clear();
	fVetoedHLT.clear();
	
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

	  // HT/jetHT dataset
	  fTriggerMap["HLT_HT500_v1"]       = &fMT2tree->trigger.HLT_HT500_v1;
	  fTriggerMap["HLT_HT500_v2"]       = &fMT2tree->trigger.HLT_HT500_v2;
	  fTriggerMap["HLT_HT500_v3"]       = &fMT2tree->trigger.HLT_HT500_v3;
	  fTriggerMap["HLT_HT500_v4"]       = &fMT2tree->trigger.HLT_HT500_v4;
	  fTriggerMap["HLT_HT550_v1"]       = &fMT2tree->trigger.HLT_HT550_v1;
	  fTriggerMap["HLT_HT550_v2"]       = &fMT2tree->trigger.HLT_HT550_v2;
	  fTriggerMap["HLT_HT550_v3"]       = &fMT2tree->trigger.HLT_HT550_v3;
	  fTriggerMap["HLT_HT550_v4"]       = &fMT2tree->trigger.HLT_HT550_v4;
	  fTriggerMap["HLT_HT600_v1"]       = &fMT2tree->trigger.HLT_HT600_v1;
	  fTriggerMap["HLT_HT600_v2"]       = &fMT2tree->trigger.HLT_HT600_v2;
	  fTriggerMap["HLT_HT600_v3"]       = &fMT2tree->trigger.HLT_HT600_v3;
	  fTriggerMap["HLT_HT600_v4"]       = &fMT2tree->trigger.HLT_HT600_v4;
	  fTriggerMap["HLT_HT650_v1"]       = &fMT2tree->trigger.HLT_HT650_v1;
	  fTriggerMap["HLT_HT650_v2"]       = &fMT2tree->trigger.HLT_HT650_v2;
	  fTriggerMap["HLT_HT650_v3"]       = &fMT2tree->trigger.HLT_HT650_v3;
	  fTriggerMap["HLT_HT650_v4"]       = &fMT2tree->trigger.HLT_HT650_v4;
	  fTriggerMap["HLT_PFHT350_v3"]     = &fMT2tree->trigger.HLT_PFHT350_v3;
	  fTriggerMap["HLT_PFHT350_v4"]     = &fMT2tree->trigger.HLT_PFHT350_v4;
	  fTriggerMap["HLT_PFHT350_v5"]     = &fMT2tree->trigger.HLT_PFHT350_v5;
	  fTriggerMap["HLT_PFHT350_v6"]     = &fMT2tree->trigger.HLT_PFHT350_v6;
	  fTriggerMap["HLT_PFHT350_v7"]     = &fMT2tree->trigger.HLT_PFHT350_v7;
	  fTriggerMap["HLT_PFHT650_v5"]     = &fMT2tree->trigger.HLT_PFHT650_v5;
	  fTriggerMap["HLT_PFHT650_v6"]     = &fMT2tree->trigger.HLT_PFHT650_v6;
	  fTriggerMap["HLT_PFHT650_v7"]     = &fMT2tree->trigger.HLT_PFHT650_v7;
	  fTriggerMap["HLT_PFHT650_v8"]     = &fMT2tree->trigger.HLT_PFHT650_v8;
	  fTriggerMap["HLT_PFHT650_v9"]     = &fMT2tree->trigger.HLT_PFHT650_v9;
	  fTriggerMap["HLT_PFHT700_v3"]     = &fMT2tree->trigger.HLT_PFHT700_v3;
	  fTriggerMap["HLT_PFHT700_v4"]     = &fMT2tree->trigger.HLT_PFHT700_v4;
	  fTriggerMap["HLT_PFHT700_v5"]     = &fMT2tree->trigger.HLT_PFHT700_v5;
	  fTriggerMap["HLT_PFHT700_v6"]     = &fMT2tree->trigger.HLT_PFHT700_v6;
	  fTriggerMap["HLT_PFHT700_v7"]     = &fMT2tree->trigger.HLT_PFHT700_v7;
	  fTriggerMap["HLT_PFHT750_v3"]     = &fMT2tree->trigger.HLT_PFHT750_v3;
	  fTriggerMap["HLT_PFHT750_v4"]     = &fMT2tree->trigger.HLT_PFHT750_v4;
	  fTriggerMap["HLT_PFHT750_v5"]     = &fMT2tree->trigger.HLT_PFHT750_v5;
	  fTriggerMap["HLT_PFHT750_v6"]     = &fMT2tree->trigger.HLT_PFHT750_v6;
	  fTriggerMap["HLT_PFHT750_v7"]     = &fMT2tree->trigger.HLT_PFHT750_v7;
	  
	  // HT/HTMHT dataset
	  fTriggerMap["HLT_PFHT350_PFMET100_v3"]     = &fMT2tree->trigger.HLT_PFHT350_PFMET100_v3;
	  fTriggerMap["HLT_PFHT350_PFMET100_v4"]     = &fMT2tree->trigger.HLT_PFHT350_PFMET100_v4;
	  fTriggerMap["HLT_PFHT350_PFMET100_v5"]     = &fMT2tree->trigger.HLT_PFHT350_PFMET100_v5;
	  fTriggerMap["HLT_PFHT350_PFMET100_v6"]     = &fMT2tree->trigger.HLT_PFHT350_PFMET100_v6;
	  fTriggerMap["HLT_PFHT350_PFMET100_v7"]     = &fMT2tree->trigger.HLT_PFHT350_PFMET100_v7;
	  fTriggerMap["HLT_PFHT400_PFMET100_v3"]     = &fMT2tree->trigger.HLT_PFHT400_PFMET100_v3;
	  fTriggerMap["HLT_PFHT400_PFMET100_v4"]     = &fMT2tree->trigger.HLT_PFHT400_PFMET100_v4;
	  fTriggerMap["HLT_PFHT400_PFMET100_v5"]     = &fMT2tree->trigger.HLT_PFHT400_PFMET100_v5;
	  fTriggerMap["HLT_PFHT400_PFMET100_v6"]     = &fMT2tree->trigger.HLT_PFHT400_PFMET100_v6;
	  fTriggerMap["HLT_PFHT400_PFMET100_v7"]     = &fMT2tree->trigger.HLT_PFHT400_PFMET100_v7;
	  
	  // MET Dataset
	  fTriggerMap["HLT_MET120_v9"]                       = &fMT2tree->trigger.HLT_MET120_v9;
	  fTriggerMap["HLT_MET120_v10"]                      = &fMT2tree->trigger.HLT_MET120_v10;
	  fTriggerMap["HLT_MET200_v9"]                       = &fMT2tree->trigger.HLT_MET200_v9;
	  fTriggerMap["HLT_MET200_v10"]                      = &fMT2tree->trigger.HLT_MET200_v10;
	  fTriggerMap["HLT_MET120_HBHENoiseCleaned_v2"]      = &fMT2tree->trigger.HLT_MET120_HBHENoiseCleaned_v2;
	  fTriggerMap["HLT_MET120_HBHENoiseCleaned_v3"]      = &fMT2tree->trigger.HLT_MET120_HBHENoiseCleaned_v3;
	  fTriggerMap["HLT_MET120_HBHENoiseCleaned_v4"]      = &fMT2tree->trigger.HLT_MET120_HBHENoiseCleaned_v4;
	  fTriggerMap["HLT_MET200_HBHENoiseCleaned_v2"]      = &fMT2tree->trigger.HLT_MET200_HBHENoiseCleaned_v2;
	  fTriggerMap["HLT_MET200_HBHENoiseCleaned_v3"]      = &fMT2tree->trigger.HLT_MET200_HBHENoiseCleaned_v3;
	  fTriggerMap["HLT_MET200_HBHENoiseCleaned_v4"]      = &fMT2tree->trigger.HLT_MET200_HBHENoiseCleaned_v4;
	  fTriggerMap["HLT_PFMET150_v2"]                     = &fMT2tree->trigger.HLT_PFMET150_v2;
	  fTriggerMap["HLT_PFMET150_v3"]                     = &fMT2tree->trigger.HLT_PFMET150_v3;
	  fTriggerMap["HLT_PFMET150_v4"]                     = &fMT2tree->trigger.HLT_PFMET150_v4;
	  fTriggerMap["HLT_PFMET180_v2"]                     = &fMT2tree->trigger.HLT_PFMET180_v2;
	  fTriggerMap["HLT_PFMET180_v3"]                     = &fMT2tree->trigger.HLT_PFMET180_v3;
	  fTriggerMap["HLT_PFMET180_v4"]                     = &fMT2tree->trigger.HLT_PFMET180_v4;
	  fTriggerMap["HLT_DiCentralPFJet30_PFMHT80_v5"]     = &fMT2tree->trigger.HLT_DiCentralPFJet30_PFMHT80_v5;
	  fTriggerMap["HLT_DiCentralPFJet30_PFMHT80_v6"]     = &fMT2tree->trigger.HLT_DiCentralPFJet30_PFMHT80_v6;
	  fTriggerMap["HLT_DiCentralPFJet30_PFMHT80_v7"]     = &fMT2tree->trigger.HLT_DiCentralPFJet30_PFMHT80_v7;
	  fTriggerMap["HLT_DiCentralPFJet50_PFMET80_v3"]     = &fMT2tree->trigger.HLT_DiCentralPFJet50_PFMET80_v3;
	  fTriggerMap["HLT_DiCentralPFJet50_PFMET80_v4"]     = &fMT2tree->trigger.HLT_DiCentralPFJet50_PFMET80_v4;
	  fTriggerMap["HLT_DiCentralPFJet50_PFMET80_v5"]     = &fMT2tree->trigger.HLT_DiCentralPFJet50_PFMET80_v5;
	  fTriggerMap["HLT_DiCentralPFJet50_PFMET80_v6"]     = &fMT2tree->trigger.HLT_DiCentralPFJet50_PFMET80_v6;
	  
	  // Multijet dataset
	  fTriggerMap["HLT_DiJet80_DiJet60_DiJet20_v1"]     = &fMT2tree->trigger.HLT_DiJet80_DiJet60_DiJet20_v1;
	  fTriggerMap["HLT_DiJet80_DiJet60_DiJet20_v2"]     = &fMT2tree->trigger.HLT_DiJet80_DiJet60_DiJet20_v2;
	  fTriggerMap["HLT_QuadJet60_DiJet20_v1"]           = &fMT2tree->trigger.HLT_QuadJet60_DiJet20_v1;
	  fTriggerMap["HLT_QuadJet60_DiJet20_v2"]           = &fMT2tree->trigger.HLT_QuadJet60_DiJet20_v2;
	  fTriggerMap["HLT_QuadJet50_v1"]                   = &fMT2tree->trigger.HLT_QuadJet50_v1;
	  fTriggerMap["HLT_QuadJet50_v2"]                   = &fMT2tree->trigger.HLT_QuadJet50_v2;
	  fTriggerMap["HLT_QuadJet70_v1"]                   = &fMT2tree->trigger.HLT_QuadJet70_v1;
	  fTriggerMap["HLT_QuadJet70_v2"]                   = &fMT2tree->trigger.HLT_QuadJet70_v2;
	  fTriggerMap["HLT_QuadJet70_v3"]                   = &fMT2tree->trigger.HLT_QuadJet70_v3;
	  fTriggerMap["HLT_QuadJet80_v1"]                   = &fMT2tree->trigger.HLT_QuadJet80_v1;
	  fTriggerMap["HLT_QuadJet80_v2"]                   = &fMT2tree->trigger.HLT_QuadJet80_v2;
	  fTriggerMap["HLT_QuadJet80_v3"]                   = &fMT2tree->trigger.HLT_QuadJet80_v3;
	  fTriggerMap["HLT_SixJet35_v1"]                    = &fMT2tree->trigger.HLT_SixJet35_v1;
	  fTriggerMap["HLT_SixJet35_v2"]                    = &fMT2tree->trigger.HLT_SixJet35_v2;
	  fTriggerMap["HLT_SixJet35_v3"]                    = &fMT2tree->trigger.HLT_SixJet35_v3;
	  fTriggerMap["HLT_SixJet45_v1"]                    = &fMT2tree->trigger.HLT_SixJet45_v1;
	  fTriggerMap["HLT_SixJet45_v2"]                    = &fMT2tree->trigger.HLT_SixJet45_v2;
	  fTriggerMap["HLT_SixJet45_v3"]                    = &fMT2tree->trigger.HLT_SixJet45_v3;

          // Jet/JetMon dataset
          fTriggerMap["HLT_PFJet320_v3"]                     = &fMT2tree->trigger.HLT_PFJet320_v3;
          fTriggerMap["HLT_PFJet320_v4"]                     = &fMT2tree->trigger.HLT_PFJet320_v4;
          fTriggerMap["HLT_PFJet320_v5"]                     = &fMT2tree->trigger.HLT_PFJet320_v5;
          fTriggerMap["HLT_PFJet260_v3"]                     = &fMT2tree->trigger.HLT_PFJet260_v3;
          fTriggerMap["HLT_PFJet260_v4"]                     = &fMT2tree->trigger.HLT_PFJet260_v4;
          fTriggerMap["HLT_PFJet260_v5"]                     = &fMT2tree->trigger.HLT_PFJet260_v5;
          fTriggerMap["HLT_PFJet200_v3"]                     = &fMT2tree->trigger.HLT_PFJet200_v3;
          fTriggerMap["HLT_PFJet200_v4"]                     = &fMT2tree->trigger.HLT_PFJet200_v4;
          fTriggerMap["HLT_PFJet200_v5"]                     = &fMT2tree->trigger.HLT_PFJet200_v5;
          fTriggerMap["HLT_PFJet140_v3"]                     = &fMT2tree->trigger.HLT_PFJet140_v3;
          fTriggerMap["HLT_PFJet140_v4"]                     = &fMT2tree->trigger.HLT_PFJet140_v4;
          fTriggerMap["HLT_PFJet140_v5"]                     = &fMT2tree->trigger.HLT_PFJet140_v5;
          fTriggerMap["HLT_DiPFJetAve200_v3"]                = &fMT2tree->trigger.HLT_DiPFJetAve200_v3;
          fTriggerMap["HLT_DiPFJetAve200_v4"]                = &fMT2tree->trigger.HLT_DiPFJetAve200_v4;
          fTriggerMap["HLT_DiPFJetAve200_v5"]                = &fMT2tree->trigger.HLT_DiPFJetAve200_v5;
          fTriggerMap["HLT_DiPFJetAve200_v6"]                = &fMT2tree->trigger.HLT_DiPFJetAve200_v6;
          fTriggerMap["HLT_DiPFJetAve140_v3"]                = &fMT2tree->trigger.HLT_DiPFJetAve140_v3;
          fTriggerMap["HLT_DiPFJetAve140_v4"]                = &fMT2tree->trigger.HLT_DiPFJetAve140_v4;
          fTriggerMap["HLT_DiPFJetAve140_v5"]                = &fMT2tree->trigger.HLT_DiPFJetAve140_v5;
          fTriggerMap["HLT_DiPFJetAve140_v6"]                = &fMT2tree->trigger.HLT_DiPFJetAve140_v6;
          fTriggerMap["HLT_DiPFJetAve80_v3"]                 = &fMT2tree->trigger.HLT_DiPFJetAve80_v3;
          fTriggerMap["HLT_DiPFJetAve80_v4"]                 = &fMT2tree->trigger.HLT_DiPFJetAve80_v4;
          fTriggerMap["HLT_DiPFJetAve80_v5"]                 = &fMT2tree->trigger.HLT_DiPFJetAve80_v5;
          fTriggerMap["HLT_DiPFJetAve80_v6"]                 = &fMT2tree->trigger.HLT_DiPFJetAve80_v6;
          fTriggerMap["HLT_DiPFJetAve40_v3"]                 = &fMT2tree->trigger.HLT_DiPFJetAve40_v3;
          fTriggerMap["HLT_DiPFJetAve40_v4"]                 = &fMT2tree->trigger.HLT_DiPFJetAve40_v4;
          fTriggerMap["HLT_DiPFJetAve40_v5"]                 = &fMT2tree->trigger.HLT_DiPFJetAve40_v5;
          fTriggerMap["HLT_DiPFJetAve40_v6"]                 = &fMT2tree->trigger.HLT_DiPFJetAve40_v6;
	  
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

//	//LHAPDF init   // FIXME
//        string PDF_SET="cteq66";
//        string PDF_PATH = "/shome/leo/Installations/LHAPDF/lhapdf-5.8.4/";
//	if(doPDF){
//	  LHAPDF::initPDFSet(PDF_PATH+"/share/lhapdf/"+PDF_SET, LHAPDF::LHGRID);
//	  nPDFs = LHAPDF::numberPDF();
//	  cout << "nPDF: " << nPDFs << endl;
//	}




}

// ***********************************************************************
// Analyze called for each event
void MT2Analysis::Analyze(){	

	// ---------------------------------------------------
	// Initialize fElecs, fBJets, fMuons, 
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
	if(basicfilled) FillTree();
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
	if(Jets.size()      > 25) {cout << "ERROR: Jets.size()    > 25: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return false;}
	if(fTaus.size()     > 8 ) {cout << "ERROR: fTaus.size()   >  8: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return false;}
	if(fElecs.size()    > 8 ) {cout << "ERROR: fElecs.size()  >  8: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return false;}
	if(fMuons.size()    > 8 ) {cout << "ERROR: fMuons.size()  >  8: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return false;}
	if(fPhotons.size()  > 8 ) {cout << "ERROR: fPhotons.size()>  8: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return false;}


	// ---------------------------------------------------------------
	// Fill jets 4-momenta & ID's 
	for(int i=0; i<Jets.size(); ++i) {
		fMT2tree->jet[i] = (MT2Jet) Jets[i];
	}
	// ---------------------------------------------------------------
	// Set NJets, NElecs, NMuons
	fMT2tree->SetNJets         (Jets.size());
	fMT2tree->NJetsIDLoose40 = fMT2tree->GetNjets(40, 2.4, 1);
	fMT2tree->NJetsIDLoose50 = fMT2tree->GetNjets(50, 2.4, 1);
	fMT2tree->SetNGenJets      (fTR->NGenJets > gNGenJets ? gNGenJets: fTR->NGenJets);
	fMT2tree->SetNJetsIDLoose  (fMT2tree->GetNjets(20, 2.4, 1));
	fMT2tree->SetNBJets        (fMT2tree->GetNBtags(3,2.0  ,20,2.4,1));
	fMT2tree->SetNBJetsHE      (fMT2tree->GetNBtags(2,1.74 ,20,2.4,1));
	fMT2tree->SetNBJetsCSVM    (fMT2tree->GetNBtags(4,0.679,20,2.4,1));
	fMT2tree->SetNBJetsCSVT    (fMT2tree->GetNBtags(4,0.898,20,2.4,1));
	fMT2tree->SetNBJets40CSVM  (fMT2tree->GetNBtags(4,0.679,40,2.4,1));
	fMT2tree->SetNBJets40CSVT  (fMT2tree->GetNBtags(4,0.898,40,2.4,1));
	fMT2tree->SetNEles         ((Int_t)fElecs.size());
	fMT2tree->SetNMuons        ((Int_t)fMuons.size());
	fMT2tree->SetNTaus         ((Int_t)fTaus.size());
	fMT2tree->SetNPhotons      ((Int_t)fPhotons.size());

	// --------------------------------------------------------------
	// match taus to jets
	for(int i=0; i<fTaus.size(); ++i){
		TLorentzVector tau;
		tau.SetPxPyPzE(fTR->TauPx[fTaus[i]],fTR->TauPy[fTaus[i]],fTR->TauPz[fTaus[i]],fTR->TauE[fTaus[i]]);
		float mindR =10000; int jindex =-1;
		for(int j=0; j<Jets.size(); ++j){
			float dR = tau.DeltaR(fMT2tree->jet[j].lv);
			if(dR < mindR && dR < 0.5) {mindR=dR; jindex = j;} 	
		}
		if(jindex ==-1) continue;
		fMT2tree->jet[jindex].NTauMatch++;
		if(fMT2tree->jet[jindex].isTauMatch > 0 && fMT2tree->jet[jindex].TauDR < mindR) continue;
		fMT2tree->jet[jindex].isTauMatch = i; // tells you if pf-jet is matched to a tau.  -1 not matched, >= 0 gives the tau index in the tau collection
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
		if(!fisData){
			fMT2tree->photon[i].MCmatchexitcode     =fTR->PhoMCmatchexitcode[fPhotons[i]]; 
			if(fTR->PhoMCmatchindex[fPhotons[i]]>=0){
				int index=fTR->PhoMCmatchindex[fPhotons[i]];
				fMT2tree->photon[i].GenJetMinDR= fTR->GenPhotonPartonMindR[index];
			}
		}
		// exit code meaning:  
		//                     0 = matched to particle that is not a photon -> fake
		//                     1 = photon  with status 3 quark or gluon as mother (hard scatter) (gluon happens, though gamma does not have color..)
		//                     2 = PromptPhoton with status 3 photon as mother,  
		//                     3 = other particle, i.e. showering.. 
		//                     negative values: no MC match. 
	}
	
	
	
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
	  	fMT2tree->ele[i].lv.SetPtEtaPhiE(fTR->ElPt [fElecs[i]], fTR->ElEta[fElecs[i]], fTR->ElPhi[fElecs[i]], fTR->ElE  [fElecs[i]]); 
		fMT2tree->ele[i].MT       = fMT2tree->GetMT(fMT2tree->ele[i].lv, 0., METlv, 0.); 
		fMT2tree->ele[i].Charge   = fTR->ElCharge[fElecs[i]];
		fMT2tree->ele[i].Iso      = ElePFIso(fElecs[i]);
		fMT2tree->ele[i].IDMedium = IsGoodMT2ElectronMediumID(fElecs[i]);
		fMT2tree->ele[i].IDLoose  = IsGoodMT2ElectronLooseID(fElecs[i]);
		fMT2tree->ele[i].IDVeto   = IsGoodMT2ElectronVetoID(fElecs[i]);
	}
	for(int i=0; i<fMuons.size(); ++i) {
	  	fMT2tree->muo[i].lv.SetPtEtaPhiM(fTR->MuPt [fMuons[i]], fTR->MuEta[fMuons[i]], fTR->MuPhi[fMuons[i]], 0.106); 
		fMT2tree->muo[i].MT       = fMT2tree->GetMT(fMT2tree->muo[i].lv, fMT2tree->muo[i].lv.M(), METlv, 0.); 
		fMT2tree->muo[i].Charge   = fTR->MuCharge[fMuons[i]];	
		fMT2tree->muo[i].Iso      = MuPFIso(fMuons[i]);
	}
	for(int i=0; i<fTaus.size(); ++i) {
	  	fMT2tree->tau[i].lv.SetPtEtaPhiE(fTR->TauPt [fTaus[i]], fTR->TauEta[fTaus[i]], 
				          fTR->TauPhi[fTaus[i]], fTR->TauE  [fTaus[i]]); 
		fMT2tree->tau[i].MT       = fMT2tree->GetMT(fMT2tree->tau[i].lv, 0., METlv, 0.); 
		fMT2tree->tau[i].Charge   = fTR->TauCharge[fTaus[i]];
		fMT2tree->tau[i].JetPt    = fTR->TauJetPt[fTaus[i]];
		fMT2tree->tau[i].JetEta   = fTR->TauJetEta[fTaus[i]];
		fMT2tree->tau[i].JetPhi   = fTR->TauJetPhi[fTaus[i]];
		fMT2tree->tau[i].JetMass  = fTR->TauJetMass[fTaus[i]];

		fMT2tree->tau[i].Isolation= 0;
		if(fTR->TauVLooseCombinedIsoDBSumPtCorr[fTaus[i]] > 0.5)
		  fMT2tree->tau[i].Isolation = 1;
		if(fTR->TauLooseCombinedIsoDBSumPtCorr[fTaus[i]] > 0.5)
		  fMT2tree->tau[i].Isolation = 2;
		if(fTR->TauMediumCombinedIsoDBSumPtCorr[fTaus[i]] > 0.5)
		  fMT2tree->tau[i].Isolation = 3;
		if(fTR->TauTightCombinedIsoDBSumPtCorr[fTaus[i]] > 0.5)
		  fMT2tree->tau[i].Isolation = 4;
		
		fMT2tree->tau[i].ElectronRej= 0;
		if(fTR->TauLooseElectronRejection[fTaus[i]]  > 0.5)
		fMT2tree->tau[i].ElectronRej= 1;
		if(fTR->TauMediumElectronRejection[fTaus[i]] > 0.5)
		fMT2tree->tau[i].ElectronRej= 2;
		if(fTR->TauTightElectronRejection[fTaus[i]]  > 0.5)
		fMT2tree->tau[i].ElectronRej= 3;
	
		fMT2tree->tau[i].MuonRej= 0;
		if(fTR->TauLooseMuonRejection[fTaus[i]]  > 0.5)
		fMT2tree->tau[i].MuonRej= 1;
		if(fTR->TauTightMuonRejection[fTaus[i]]  > 0.5)
		fMT2tree->tau[i].MuonRej= 3;
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
		//allow for genleptons from chargino or slepton decays
		//allow for genleptons from Higgs decays --> changed (MID>24) to the thing below
		if( (ID == 11 || ID == 12 || ID == 13 || ID == 14 || ID == 15 || ID == 16) && 
		    (MID > 37 && ( MID<1000001||MID>1000039) && (MID<2000001||MID>2000015) ) ) continue;
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
		fMT2tree->genlept[NGenLepts-1].Mlv.SetPtEtaPhiM(fTR->GenLeptonMPt[i], fTR->GenLeptonMEta[i], fTR->GenLeptonMPhi[i], 0);//how to get the mass
		fMT2tree->genlept[NGenLepts-1].GMlv.SetPtEtaPhiM(fTR->GenLeptonGMPt[i], fTR->GenLeptonGMEta[i], fTR->GenLeptonGMPhi[i], 0);//how to get mass
		fMT2tree->genlept[NGenLepts-1].ID       = fTR->GenLeptonID[i];
		fMT2tree->genlept[NGenLepts-1].MID      = fTR->GenLeptonMID[i];
		fMT2tree->genlept[NGenLepts-1].MStatus  = fTR->GenLeptonMStatus[i];
		fMT2tree->genlept[NGenLepts-1].GMID     = fTR->GenLeptonGMID[i];
		fMT2tree->genlept[NGenLepts-1].GMStatus = fTR->GenLeptonGMStatus[i];
		if(abs(fMT2tree->genlept[NGenLepts-1].ID) == 11 || abs(fMT2tree->genlept[NGenLepts-1].ID) == 13 || abs(fMT2tree->genlept[NGenLepts-1].ID) == 15   ){
			fMT2tree->genlept[NGenLepts-1].MT = fMT2tree->GetMT(fMT2tree->genlept[NGenLepts-1].lv, fMT2tree->genlept[NGenLepts-1].lv.M(), fMT2tree->genmet[0], 0.);
		}
	}
	fMT2tree->NGenLepts      = NGenLepts;


	// --------------------------------------------------------------------
	// MET 
	fMT2tree->pfmet[0]=MET();
	// raw, uncorrected and not modified pfmet
	fMT2tree->rawpfmet[0].SetPtEtaPhiM(fTR->PFMET, 0, fTR->PFMETphi, 0);
	// type1corrected pfMET
	fMT2tree->type1pfmet[0].SetPtEtaPhiM(fTR->PFType1MET, 0, fTR->PFType1METphi, 0);
	// raw calomet
	fMT2tree->misc.CaloMETRaw = fTR->RawMET;
	fMT2tree->misc.CaloMETMuJesCorr = fTR->MuJESCorrMET;

//	//pdf weights // FIXME
//	if(doPDF){
//          fMT2tree->NPdfs = nPDFs;
//          fMT2tree->pdfW[0]=1;
//	  LHAPDF::initPDF(0);
//
//          float pdf01 = LHAPDF::xfx(fTR->PDFx1, fTR->PDFScalePDF, fTR->PDFID1)/fTR->PDFx1 ;
//          float pdf02 = LHAPDF::xfx(fTR->PDFx2, fTR->PDFScalePDF, fTR->PDFID2)/fTR->PDFx2 ;
//
//          for(int pdf=1; pdf<= nPDFs; pdf++){
//	    LHAPDF::initPDF(pdf);
//            float pdf1 = LHAPDF::xfx(fTR->PDFx1, fTR->PDFScalePDF, fTR->PDFID1)/fTR->PDFx1 ;
//            float pdf2 = LHAPDF::xfx(fTR->PDFx2, fTR->PDFScalePDF, fTR->PDFID2)/fTR->PDFx2 ;
//            fMT2tree->pdfW[pdf] = pdf1/pdf01*pdf2/pdf02;
//          }
//        }

	//MC info ---------------------------------------------------------------
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

	// Pile UP info and reco vertices -------------------------------------------------------
	if(!fisData){
		fMT2tree->pileUp.PUnumInt          = fTR->PUnumInteractions;        
		fMT2tree->pileUp.PUtrueNumInt      = fTR->PUnumTrueInteractions;
		fMT2tree->pileUp.PUnumIntLate      = fTR->PUOOTnumInteractionsLate;   // branch added in V02-03-01 
		fMT2tree->pileUp.PUnumIntEarly     = fTR->PUOOTnumInteractionsEarly;  // branch added in V02-03-01 
		fMT2tree->pileUp.PtHat             = fTR->PtHat;
		fMT2tree->pileUp.PUScenario        = (int) fPUScenario;

		if       (fPUScenario==noPU  )  {fMT2tree->pileUp.Weight            = 1;}
		else if  (fPUScenario==MC2012)  {fMT2tree->pileUp.Weight            = GetPUWeight(fTR->PUnumTrueInteractions);}
		if(fVerbose > 3) {
			cout << "fPUScenario " << fPUScenario <<  " fTR->PUnumInteractions " <<  fTR->PUnumInteractions << " weight "  
		     	     << " fMT2tree->pileUp.Weight "         << fMT2tree->pileUp.Weight << endl; 
		}
	}
	fMT2tree->pileUp.Rho               = fTR->Rho; // ATTENTION: this rho is from KT6 PF jets without pf-CHS
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
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon150_v1";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon150_v2";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon150_v3";
		singlePhotonTriggers[photontriggernumber++] = "HLT_Photon150_v4";

		bool SinglePhotFired(false);
		for(int i=0; i<photontriggernumber; ++i){
			if( GetHLTResult(singlePhotonTriggers[i])) SinglePhotFired=true;
		}
		if(SinglePhotFired) fMT2tree->trigger.HLT_SinglePhotons = true;

		// di-electron triggers ------------------------------------------------------------------------------------------------------
		string diElectronTriggers[100];
		int diElectronTriggernumber=0;
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17";
		diElectronTriggers[diElectronTriggernumber++] = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18";
		bool DiElectronFired(false);
		for(int i=0; i<diElectronTriggernumber; ++i){
			if(GetHLTResult(diElectronTriggers[i])) DiElectronFired=true;
		}
		if(DiElectronFired) fMT2tree->trigger.HLT_DiElectrons =true;

		// DiMuon triggers
		string diMuonTriggers[100]; 
		int diMuonTriggernumber=0;
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v16";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v17";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v18";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v19";
	//	diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v20";//is not there?
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_Mu8_v21";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_TkMu8_v9";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_TkMu8_v10";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_TkMu8_v11";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_TkMu8_v12";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu17_TkMu8_v13";
/*
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v16";//prescaled
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v17";//prescaled
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v18";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v19";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v20";//is not there?
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu13_Mu8_v21";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu22_TkMu8_v4";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu22_TkMu8_v5";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu22_TkMu8_v6";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu22_TkMu8_v7";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu22_TkMu8_v8";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu22_TkMu22_v4";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu22_TkMu22_v5";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu22_TkMu22_v6";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu22_TkMu22_v7";
		diMuonTriggers[diMuonTriggernumber++] = "HLT_Mu22_TkMu22_v8";
*/
		bool DiMuonFired(false);
		for(int i=0; i<diMuonTriggernumber; ++i){
			if(GetHLTResult(diMuonTriggers[i])) DiMuonFired=true;
		}
		if(DiMuonFired) fMT2tree->trigger.HLT_DiMuons =true;

		//EMu
		string EMuTriggers[100];
		int EMuTriggernumber=0;
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8";
		EMuTriggers[EMuTriggernumber++] = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9";
		bool EMuFired(false);
		for(int i=0; i<EMuTriggernumber; ++i){
			if(GetHLTResult(EMuTriggers[i])) EMuFired=true;
		}
		if(EMuFired) fMT2tree->trigger.HLT_EMu =true;

		//SingleMu
		string SingleMuTriggers[100];
		int SingleMuTriggerNumber=0;
		SingleMuTriggers[SingleMuTriggerNumber++] = "HLT_IsoMu24_eta2p1_v11";
		SingleMuTriggers[SingleMuTriggerNumber++] = "HLT_IsoMu24_eta2p1_v12";
		SingleMuTriggers[SingleMuTriggerNumber++] = "HLT_IsoMu24_eta2p1_v13";
		SingleMuTriggers[SingleMuTriggerNumber++] = "HLT_IsoMu24_eta2p1_v14";
		SingleMuTriggers[SingleMuTriggerNumber++] = "HLT_IsoMu24_eta2p1_v15";
		/*
		SingleMuTriggers[SingleMuTriggerNumber++] = "HLT_IsoMu24_v15"; // only from run 193834 on
		SingleMuTriggers[SingleMuTriggerNumber++] = "HLT_IsoMu24_v16";
		SingleMuTriggers[SingleMuTriggerNumber++] = "HLT_IsoMu24_v17";
		*/
		bool SingleMuFired(false);
		for(int i=0; i<SingleMuTriggerNumber; ++i){
			if(GetHLTResult(SingleMuTriggers[i])) SingleMuFired=true;
		}
		if(SingleMuFired) fMT2tree->trigger.HLT_SingleMu =true;

		//SingleMuJet
		string SingleMuJetTriggers[100];
		int SingleMuJetTriggerNumber=0;
		SingleMuJetTriggers[SingleMuJetTriggerNumber++] = "HLT_IsoMu20_eta2p1_CentralPFJet80_v3";
		SingleMuJetTriggers[SingleMuJetTriggerNumber++] = "HLT_IsoMu20_eta2p1_CentralPFJet80_v4";
		SingleMuJetTriggers[SingleMuJetTriggerNumber++] = "HLT_IsoMu20_eta2p1_CentralPFJet80_v5";
		SingleMuJetTriggers[SingleMuJetTriggerNumber++] = "HLT_IsoMu20_eta2p1_CentralPFJet80_v6";
		SingleMuJetTriggers[SingleMuJetTriggerNumber++] = "HLT_IsoMu20_eta2p1_CentralPFJet80_v7";
		SingleMuJetTriggers[SingleMuJetTriggerNumber++] = "HLT_IsoMu20_eta2p1_CentralPFJet80_v8";
		bool SingleMuJetFired(false);
		for(int i=0; i<SingleMuJetTriggerNumber; ++i){
			if(GetHLTResult(SingleMuJetTriggers[i])) SingleMuJetFired=true;
		}
		if(SingleMuJetFired) fMT2tree->trigger.HLT_SingleMu_Jet =true;

		//SingleMuDiJet // only from run 193834 on -->there is a similar pass including MET
		string SingleMuDiJetTriggers[100];
		int SingleMuDiJetTriggerNumber=0;
		SingleMuDiJetTriggers[SingleMuDiJetTriggerNumber++] = "HLT_IsoMu24_CentralPFJet30_CentralPFJet25_v1";
		SingleMuDiJetTriggers[SingleMuDiJetTriggerNumber++] = "HLT_IsoMu24_CentralPFJet30_CentralPFJet25_v2";
		SingleMuDiJetTriggers[SingleMuDiJetTriggerNumber++] = "HLT_IsoMu24_CentralPFJet30_CentralPFJet25_v3";
		SingleMuDiJetTriggers[SingleMuDiJetTriggerNumber++] = "HLT_IsoMu24_CentralPFJet30_CentralPFJet25_v4";
		bool SingleMuDiJetFired(false);
		for(int i=0; i<SingleMuDiJetTriggerNumber; ++i){
			if(GetHLTResult(SingleMuDiJetTriggers[i])) SingleMuDiJetFired=true;
		}
		if(SingleMuDiJetFired) fMT2tree->trigger.HLT_SingleMu_DiJet =true;

		//SingleMuTriJet
		string SingleMuTriJetTriggers[100];
		int SingleMuTriJetTriggerNumber=0;
		SingleMuTriJetTriggers[SingleMuTriJetTriggerNumber++] = "HLT_IsoMu17_eta2p1_TriCentralPFJet30_v2";
		SingleMuTriJetTriggers[SingleMuTriJetTriggerNumber++] = "HLT_IsoMu17_eta2p1_TriCentralPFJet30_v3";
		SingleMuTriJetTriggers[SingleMuTriJetTriggerNumber++] = "HLT_IsoMu17_eta2p1_TriCentralPFJet30_v4";
		SingleMuTriJetTriggers[SingleMuTriJetTriggerNumber++] = "HLT_IsoMu17_eta2p1_TriCentralPFJet30_v5";
		SingleMuTriJetTriggers[SingleMuTriJetTriggerNumber++] = "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_v1";
		SingleMuTriJetTriggers[SingleMuTriJetTriggerNumber++] = "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1";
		SingleMuTriJetTriggers[SingleMuTriJetTriggerNumber++] = "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v2";
		SingleMuTriJetTriggers[SingleMuTriJetTriggerNumber++] = "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_v3";
		SingleMuTriJetTriggers[SingleMuTriJetTriggerNumber++] = "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_v4";
		bool SingleMuTriJetFired(false);
		for(int i=0; i<SingleMuTriJetTriggerNumber; ++i){
			if(GetHLTResult(SingleMuTriJetTriggers[i])) SingleMuTriJetFired=true;
		}
		if(SingleMuTriJetFired) fMT2tree->trigger.HLT_SingleMu_TriJet =true;

		//SingleEleDiJetMET
		string SingleEleDiJetMETTriggers[100];
		int SingleEleDiJetMETTriggerNumber=0;
		SingleEleDiJetMETTriggers[SingleEleDiJetMETTriggerNumber++] = "HLT_Ele27_WP80_CentralPFJet30_CentralPFJet25_PFMET20_v2";
		SingleEleDiJetMETTriggers[SingleEleDiJetMETTriggerNumber++] = "HLT_Ele27_WP80_CentralPFJet30_CentralPFJet25_PFMET20_v3";
		SingleEleDiJetMETTriggers[SingleEleDiJetMETTriggerNumber++] = "HLT_Ele27_WP80_CentralPFJet30_CentralPFJet25_PFMET20_v4";
		SingleEleDiJetMETTriggers[SingleEleDiJetMETTriggerNumber++] = "HLT_Ele32_WP80_CentralPFJet35_CentralPFJet25_PFMET20_v1";
		SingleEleDiJetMETTriggers[SingleEleDiJetMETTriggerNumber++] = "HLT_Ele32_WP80_CentralPFJet35_CentralPFJet25_PFMET20_v2";
		SingleEleDiJetMETTriggers[SingleEleDiJetMETTriggerNumber++] = "HLT_Ele32_WP80_CentralPFJet35_CentralPFJet25_PFMET20_v3";
		SingleEleDiJetMETTriggers[SingleEleDiJetMETTriggerNumber++] = "HLT_Ele32_WP80_CentralPFJet35_CentralPFJet25_PFMET20_v4";
		bool SingleEleDiJetMETFired(false);
		for(int i=0; i<SingleEleDiJetMETTriggerNumber; ++i){
			if(GetHLTResult(SingleEleDiJetMETTriggers[i])) SingleEleDiJetMETFired=true;
		}
		if(SingleEleDiJetMETFired) fMT2tree->trigger.HLT_SingleEle_DiJet_MET =true;

		//SingleEleMETMT
		string SingleEleMETMTTriggers[100];
		int SingleEleMETMTTriggerNumber=0;
		SingleEleMETMTTriggers[SingleEleMETMTTriggerNumber++] = "HLT_Ele27_WP80_PFMET_MT50_v2";
		SingleEleMETMTTriggers[SingleEleMETMTTriggerNumber++] = "HLT_Ele27_WP80_PFMET_MT50_v3";
		SingleEleMETMTTriggers[SingleEleMETMTTriggerNumber++] = "HLT_Ele27_WP80_PFMET_MT50_v4";
		SingleEleMETMTTriggers[SingleEleMETMTTriggerNumber++] = "HLT_Ele27_WP80_PFMET_MT50_v5";
		SingleEleMETMTTriggers[SingleEleMETMTTriggerNumber++] = "HLT_Ele27_WP80_PFMET_MT50_v6";
		SingleEleMETMTTriggers[SingleEleMETMTTriggerNumber++] = "HLT_Ele27_WP80_PFMET_MT50_v7";
		bool SingleEleMETMTFired(false);
		for(int i=0; i<SingleEleMETMTTriggerNumber; ++i){
			if(GetHLTResult(SingleEleMETMTTriggers[i])) SingleEleMETMTFired=true;
		}
		if(SingleEleDiJetMETFired) fMT2tree->trigger.HLT_SingleEle_MET_MT =true;
	}
	// ___________________________________________________________________________
	
	// ------------------------------------------------------------------
	// fill misc 
	fMT2tree->misc.isData              = fisData;
	fMT2tree->misc.isType1MET          = fisType1MET;
	fMT2tree->misc.isCHSJets           = fisCHSJets;
	fMT2tree->misc.isFastSim           = fisFastSim;
	fMT2tree->misc.QCDPartonicHT       = fTR->QCDPartonicHT;     
	fMT2tree->misc.ProcessID           = fID;
	fMT2tree->misc.Run                 = fTR->Run;
	fMT2tree->misc.Event		   = fTR->Event;
	fMT2tree->misc.LumiSection	   = fTR->LumiSection;
	fMT2tree->misc.HT                  = fHT;
	
	fMT2tree->misc.MET                 = MET().Pt();
	fMT2tree->misc.METPhi              = MET().Phi();

	fMT2tree->misc.LeadingJPt          = (fMT2tree->NJets > 0) ? fMT2tree->jet[0].lv.Pt() : 0;
	fMT2tree->misc.SecondJPt           = (fMT2tree->NJets > 1) ? fMT2tree->jet[1].lv.Pt() : 0;
	
	// noise filters    // store bad events as "true"
	
	fMT2tree->misc.CSCTightHaloIDFlag                  = fTR->CSCTightHaloID                    ? 0:1;                  
	fMT2tree->misc.HBHENoiseFlag                       = fTR->HBHENoiseFilterResult             ? 0:1;
	fMT2tree->misc.hcalLaserEventFlag                  = fTR->hcalLaserEventFilter              ? 0:1;
	fMT2tree->misc.trackingFailureFlag                 = fTR->trackingFailureFilter             ? 0:1;
	fMT2tree->misc.eeBadScFlag                         = fTR->eeBadScFilter                     ? 0:1;
	fMT2tree->misc.EcalDeadCellTriggerPrimitiveFlag    = fTR->EcalDeadCellTriggerPrimitiveFilter? 0:1;
	fMT2tree->misc.CrazyHCAL                           = fCrazyHCAL;                 
	fMT2tree->misc.NegativeJEC                         = fNegativeJEC;
	
	// add gen photon
	if(!fisData){
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
	fMT2tree->FillMT2HemiMinDHT(0,0,1,20,2.4,1,1);  
	
	// store MT2 misc variables
	fMT2tree->misc.MT2                 = fMT2tree->hemi[0].MT2;    // note: this is a bit dangerous, 
	fMT2tree->misc.MCT                 = fMT2tree->hemi[0].MCT;
	fMT2tree->misc.MT2jet40            = fMT2tree->GetMT2(0,false,1,40,2.4,2,3,1);
	fMT2tree->misc.MT2jet50            = fMT2tree->GetMT2(0,false,1,50,2.4,2,3,1);


	// other variables to be computed based on MT2tree
	// ----------------------------------------------------------------------------------
	fMT2tree->misc.Vectorsumpt	   = fMT2tree->GetMHTminusMET(1, 20, 2.4, true); // including leptons, ID jets only
	fMT2tree->misc.MinMetJetDPhi       = fMT2tree->MinMetJetDPhi(0,20,5.0,1);
	fMT2tree->misc.MinMetJetDPhi4      = fMT2tree->MinMetJetDPhi(0,20,5.0,1,4);  // use first 4 jets
	fMT2tree->misc.MinMetJetDPhiPt40   = fMT2tree->MinMetJetDPhi(0,40,5.0,1);    // use jets having pT>40 GeV
	fMT2tree->misc.MinMetJetDPhiPt50   = fMT2tree->MinMetJetDPhi(0,50,5.0,1);    // use jets having pT>50 GeV
	fMT2tree->misc.MinMetJetDPhi4Pt40  = fMT2tree->MinMetJetDPhi(0,40,5.0,1,4);  // use first 4 jets having pT>40 GeV
	fMT2tree->misc.MinMetJetDPhi4Pt50  = fMT2tree->MinMetJetDPhi(0,50,5.0,1,4);  // use first 4 jets having pT>50 GeV
	fMT2tree->misc.MinMetJetDPhiInCrack= fMT2tree->MinMetJetDPhi(0,20,5.0,0,99,true);
	fMT2tree->misc.MinMetJetDPhiIndex  = fMT2tree->MinMetJetDPhiIndex(0,20,5.0,1);
	fMT2tree->misc.MinMetBJetDPhi      = fMT2tree->BJetMETdPhi(2,1.74,20,5,1);  // minmetjet dPhi w.r.t SSVHEM bjets with pt > 20 and no eta restriction. 
	fMT2tree->misc.PassJetID           = fMT2tree->PassJetID(50,2.4,1);
	fMT2tree->misc.PassJet40ID         = fMT2tree->PassJetID(40,2.4,1);
	if(fMT2tree->NJets > 0) {
		fMT2tree->misc.Jet0Pass      = (Int_t) fMT2tree->jet[0].IsGoodPFJet(100,2.4,1);
	} else  fMT2tree->misc.Jet0Pass      = 0; 
	if(fMT2tree->NJets > 1) {
		fMT2tree->misc.Jet1Pass      = (Int_t) fMT2tree->jet[1].IsGoodPFJet(100,2.4,1);
	} else  fMT2tree->misc.Jet1Pass      = 0;
	
	// MHT from jets and taus
	fMT2tree->MHT[0]=fMT2tree->GetMHTlv(1, 20, 2.4, true); // only jets satisfying the loose PF-ID and leptons
	

	// pf HT and MHT
	float pfHT30=0,pfHT35=0,pfHT40=0,pfHT45=0,pfHT50=0;
	for(int j=0; j<fTR->NJets; ++j){ // pfHT from Jets without CHS 
	  if(fTR->JEcorr[j] <0) continue;
	  if(fTR->JPt[j] > 30 && fabs(fTR->JEta[j])<3.0)   pfHT30 += fTR->JPt[j];	  
	  if(fTR->JPt[j] > 35 && fabs(fTR->JEta[j])<3.0)   pfHT35 += fTR->JPt[j];	  
	  if(fTR->JPt[j] > 40 && fabs(fTR->JEta[j])<3.0)   pfHT40 += fTR->JPt[j];	  
	  if(fTR->JPt[j] > 45 && fabs(fTR->JEta[j])<3.0)   pfHT45 += fTR->JPt[j];	  
	  if(fTR->JPt[j] > 50 && fabs(fTR->JEta[j])<3.0)   pfHT50 += fTR->JPt[j];	  
	}
	fMT2tree->misc.pfHT30       = pfHT30;
	fMT2tree->misc.pfHT35       = pfHT35;
	fMT2tree->misc.pfHT40       = pfHT40;
	fMT2tree->misc.pfHT45       = pfHT45;
	fMT2tree->misc.pfHT50       = pfHT50;

	// W and Top decay modes
	fMT2tree->misc.WDecayMode   = fMT2tree->WDecayMode  ();
	fMT2tree->misc.TopDecayMode = fMT2tree->TopDecayMode();

	//btag SF --------------------------------------------------------------------------------------------------------
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
	// -------------------------------------------------------------------------------------------------
}

// *****************************************************************************
// initialize event 
void MT2Analysis::InitializeEvent(){
	GetLeptonJetIndices();
}

void MT2Analysis::GetLeptonJetIndices(){
	fIsNANObj = false; 

	fElecs.clear();
	fMuons.clear();
	fTaus.clear();
	fPhotons.clear();
	fPhotonJetOverlapRemoved.clear();
	Jets.clear();
	
	// Photons -----------------
	vector<float> photon_pts;
	for(int i=0; i<fTR->NPhotons; ++i){
		if(! IsGoodPhoton(i))                             continue; 
		if(! IsGoodPhotonEGMLooseRelISO(i))               continue;   // preselection: use only photons passing the loose RelIso
		fPhotons.push_back(i);
		photon_pts.push_back(fTR->PhoPt[i]);
		fPhotonJetOverlapRemoved.push_back(false);
	}
	fPhotons     = Util::VSort(fPhotons, photon_pts);

  	// #--- muon loop
	vector<float> muloose;
  	for(int muIndex=0;muIndex<fTR->NMus;muIndex++){
		if(! IsGoodMT2Muon(muIndex)) continue;
		fMuons.push_back(muIndex);
		muloose.push_back(fTR->MuPt[muIndex]);
	}
	fMuons      = Util::VSort(fMuons     , muloose);
	
  	// #--- electron loop
	vector<float> eltight;
  	for(int elIndex=0;elIndex<fTR->NEles;elIndex++){
		if(! IsGoodMT2ElectronVetoID(elIndex)) continue;
		fElecs.push_back(elIndex);
		eltight.push_back(fTR->ElPt[elIndex]);
	}
	fElecs      = Util::VSort(fElecs     , eltight);
	
	// Warning: taus are also contained in the jet-collection.: overlap is not removed! 
	vector<float> taus;
	for(int i=0; i< fTR->TauNObjs; ++i){
		if(std::isnan(fTR->TauPt[i])) { fIsNANObj = true; continue;} //protection against objects with NAN-Pt
		if(!(IsGoodTau(i))          ) continue;
		fTaus.push_back(i);
		taus.push_back(fTR->TauPt[i]);
	}
	fTaus          = Util::VSort(fTaus     , taus);
	

	// Jets loop ------------------------------------------------------------------
	for(int ij=0; ij < (fisCHSJets?fTR->PFCHSNJets:fTR->NJets); ++ij){
		MT2AnalysisJet* jet = new MT2AnalysisJet(ij, "PF", this);
		if( jet->Pt()    < 20)  continue; 
		if( jet->Scale   < 0 )  continue; 

		// Delta R (jet-lepton cleaning)
		Bool_t jGood(true);
		for (int elIndex=0; elIndex<fElecs.size(); ++elIndex){
			TLorentzVector el;
			el.SetPtEtaPhiE(fTR->ElPt[fElecs[elIndex]], fTR->ElEta[fElecs[elIndex]], fTR->ElPhi[fElecs[elIndex]], fTR->ElE[fElecs[elIndex]]);
			if(jet->lv.DeltaR(el)<0.4) {jGood=false;}
		}
		for (int muIndex=0; muIndex<fMuons.size(); ++muIndex){
			TLorentzVector mu;
			mu.SetPtEtaPhiE(fTR->MuPt[fMuons[muIndex]], fTR->MuEta[fMuons[muIndex]], fTR->MuPhi[fMuons[muIndex]], fTR->MuE[fMuons[muIndex]]);
			if(jet->lv.DeltaR(mu)<0.4) {jGood=false;}
		}
		for (int iphoton=0; iphoton<fPhotons.size(); ++iphoton){
			TLorentzVector phot(0,0,0,0);
			phot.SetPtEtaPhiM(fTR->PhoPt[fPhotons[iphoton]], fTR->PhoEta[fPhotons[iphoton]], fTR->PhoPhi[fPhotons[iphoton]], 0.);
			if(jet->lv.DeltaR(phot)<0.2 && fRemovePhoton ) {jGood=false; fPhotonJetOverlapRemoved[iphoton]=true; }
		}
		if (! jGood) continue; 

		Jets.push_back(*jet); delete jet;
	}
	// -----------
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
	for(int j=0; j<Jets.size(); ++j){
	  if(Jets[j].Pt() > 50 && fabs(Jets[j].Eta())<3.0){
	    HT += Jets[j].Pt();
	  }
	}
	fHT = HT;
	if(HT<fCut_HT_min){return false;}
	
	// leading jets including JID for jets
	bool leadingjets(true);
	if(fCut_JPt_hardest_min > 0){
		if(Jets.size() <1) leadingjets=false;
		if(! Jets[0].IsGoodPFJet(fCut_JPt_hardest_min, 2.4,1)       ) {leadingjets=false;} 
	}
	if(fCut_JPt_second_min > 0){
		if(Jets.size() <2) leadingjets=false;
		if(! Jets[1].IsGoodPFJet(fCut_JPt_second_min, 2.4,1)       ) {leadingjets=false;} 
	}
	if(leadingjets == false) return false;
	
	// flag crazy events (dedicated to filter HCAL "laser" events)
	if(fTR->HCALSumEt > 10000 && fTR->NJets > 25 ){
	  fCrazyHCAL =1;
		cout << "WARNING: crazy HCAL event: Run " << fTR->Run << " LS " << fTR->LumiSection  << " Event " << fTR->Event 
		     << " has HCALSumET " << fTR->HCALSumEt << " njets " << fTR->NJets << endl;
	}
	else    fCrazyHCAL =0;	

	// flag events with Jets with correction factor from JE correction 
	fNegativeJEC =0; 
	for( int i=0; i<fTR->NJets; ++i){
		if(fTR->JEcorr[i]>=0) continue;
		fNegativeJEC =1;
	}
	if(fNegativeJEC && fVerbose >1) {
		cout << "WARNING: Event with Jets with negative JEC: Run " << fTR->Run << " LS " << fTR->LumiSection << " Event " << fTR->Event
		     << " N PFJets " << fTR->NJets << " NCaloJets " << fTR->CANJets << endl; 
	}
	//only electrons and muons
	if(fCut_NLeptons_min > 0){
		if(fElecs.size() + fMuons.size()<fCut_NLeptons_min) {return false;}
	}
	if(fCut_NJets40_min > 0){
	int NJets40_min=0;
	for(int j=0; j<Jets.size(); ++j){
	  if(Jets[j].IsGoodPFJet(40, 2.4,1) ) {++NJets40_min;}
	}
		if(NJets40_min<fCut_NJets40_min) {return false;}
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
		} else if( !strcmp(ParName, "JPt_hardest_min") ){
			fCut_JPt_hardest_min      = float(ParValue); ok = true;			
		} else if( !strcmp(ParName, "JPt_second_min") ){
			fCut_JPt_second_min      = float(ParValue); ok = true;			
		} else if( !strcmp(ParName, "PtHat_max")){
			fCut_PtHat_max            = float(ParValue); ok = true;
		} 
		else if( !strcmp(ParName, "NJets40_min")){
			fCut_NJets40_min         = int(ParValue); ok = true;
		} 
		//only electrons and muons
		else if( !strcmp(ParName, "NLeptons_min")){
			fCut_NLeptons_min         = int(ParValue); ok = true;
		} 

		if(!ok) cout << "%% MT2Analysis::ReadCuts ==> ERROR: Unknown variable " << ParName << endl;
	}	
	if(verbose){
		cout << "setting cuts to: " << endl;
		cout << "  PFMET_min                   " << fCut_PFMET_min                  <<endl;
		cout << "  HT_min                      " << fCut_HT_min                     <<endl;
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
		cout << "--------------"    << endl;	
	}			
}

//***************************************************************************************************
// Muon Selector
bool MT2Analysis::IsGoodMT2Muon(const int index){
	// Acceptance
	if (!(fTR->MuPt[index] > 10) )             return false;
	if (!(fabs(fTR->MuEta[index])<2.4) )       return false;
	// Quality cuts
	if ( !fTR->MuIsGlobalMuon[index] )         return false;
//	if ( !fTR->MuIsPFMuon[index] )             return false;  // optional cut 
	// Hits
	if ( !(fTR->MuNChi2[index] < 10) )         return false;
	if ( !(fTR->MuNMuHits[index] > 0) )        return false;   //  muon.outerTrack()->hitPattern().numberOfValidHits()
	if ( !(fTR->MuNPxHits[index] > 0) )        return false;   //  muon.innerTrack()->hitPattern().numberOfValidPixelHits()
	if ( !(fTR->MuNMatches[index] > 1) )       return false;   //  numberOfMatches()
	if ( !(fTR->MuNSiLayers[index] > 5) )      return false;   //  track()->hitPattern().trackerLayersWithMeasurement() 
	// Vertex compatibility
	if ( !(fabs(fTR->MuD0PV[index]) < 0.2 ) )  return false;
	if ( !(fabs(fTR->MuDzPV[index]) < 0.5 ) )  return false;
	// Isolation
	if ( !(MuPFIso(index) < 0.2))              return false;

	return true;
}
float MT2Analysis::MuPFIso(const int index){
	  double neutral = (fTR->MuPfIsoR03NeHad[index] + fTR->MuPfIsoR03Photon[index] - 0.5*fTR->MuPfIsoR03SumPUPt[index] );
	  float iso      = (fTR->MuPfIsoR03ChHad[index] + TMath::Max(0.0, neutral) ) / fTR->MuPt[index];
	  return iso;
}

//***************************************************************************************************
// Electron Selector
bool MT2Analysis::IsGoodMT2ElectronVetoID(const int index){
  	if(!(fabs(fTR->ElEta[index]) < 2.4) ) return false;
  	if(!(fabs(fTR->ElPt[index]) > 10.0 ) ) return false;

  	// ECAL gap veto
  	if ( fabs(fTR->ElSCEta[index]) > 1.4442 && fabs(fTR->ElSCEta[index]) < 1.566 )  return false;  
	
	// Medium Working Point
	if ( fabs(fTR->ElEta[index]) < 1.479 ) { // Barrel
		if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.007)) return false;
		if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.8))   return false;
		if(!(fTR->ElSigmaIetaIeta[index]<0.01))                    return false;
		if(!(fTR->ElHcalOverEcal[index]<0.15))                     return false;
	} else { // Endcap
		if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.01)) return false;
		if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.7))  return false;
		if(!(fTR->ElSigmaIetaIeta[index]<0.03))                   return false;
	}
  
	// Vertex
	if(!(abs(fTR->ElD0PV[index])<0.04))                               return false;
	if(!(abs(fTR->ElDzPV[index])<0.2))                                return false;


	// Iso
	float pfIso = ElePFIso(index);
	if ( !(pfIso  < 0.15 ) )                                          return false;

	return true;
}

bool MT2Analysis::IsGoodMT2ElectronLooseID(const int index){
  	if(!(fabs(fTR->ElEta[index]) < 2.4) ) return false;
  	if(!(fabs(fTR->ElPt[index]) > 10.0 ) ) return false;

  	// ECAL gap veto
  	if ( fabs(fTR->ElSCEta[index]) > 1.4442 && fabs(fTR->ElSCEta[index]) < 1.566 )  return false;  
	
	// Medium Working Point
	if ( fabs(fTR->ElEta[index]) < 1.479 ) { // Barrel
		if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.007)) return false;
		if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.15))  return false;
		if(!(fTR->ElSigmaIetaIeta[index]<0.01))                    return false;
		if(!(fTR->ElHcalOverEcal[index]<0.12)) return false;
	} else { // Endcap
		if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.009)) return false;
		if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.10))  return false;
		if(!(fTR->ElSigmaIetaIeta[index]<0.03)) return false;
		if(!(fTR->ElHcalOverEcal[index]<0.10)) return false;
	}
  
	// Vertex
	if(!(abs(fTR->ElD0PV[index])<0.02)) return false;
	if(!(abs(fTR->ElDzPV[index])<0.2)) return false;

	// Conversion rejection
  	if(!(fTR->ElPassConversionVeto[index])) return false;
 	if(!(fTR->ElNumberOfMissingInnerHits[index]<=1)) return false;

	// |1/e-1/p|<0.05
	float e=fTR->ElCaloEnergy[index];
	float p=fTR->ElCaloEnergy[index]/fTR->ElESuperClusterOverP[index];
	if(!(fabs(1./e-1./p)<0.05)) return false;
  

  	float pfIso = ElePFIso(index);
  	if ( fabs(fTR->ElEta[index]) < 1.479 || fTR->ElPt[index]>20.0) { // Barrel
		if ( !((pfIso  < 0.15) ) ) return false;
	} else {
    	//Endcap with pt<20
    		if ( !((pfIso  < 0.10) ) ) return false;
  	}
	return true;
}

bool MT2Analysis::IsGoodMT2ElectronMediumID(const int index){
  	if(!(fabs(fTR->ElEta[index]) < 2.4) ) return false;
  	if(!(fabs(fTR->ElPt[index]) > 10.0 ) ) return false;

  	// ECAL gap veto
  	if ( fabs(fTR->ElSCEta[index]) > 1.4442 && fabs(fTR->ElSCEta[index]) < 1.566 )  return false;  
	
	// Medium Working Point
	if ( fabs(fTR->ElEta[index]) < 1.479 ) { // Barrel
		if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.004)) return false;
		if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.06))  return false;
		if(!(fTR->ElSigmaIetaIeta[index]<0.01))                    return false;
		if(!(fTR->ElHcalOverEcal[index]<0.12)) return false;
	} else { // Endcap
		if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.007)) return false;
		if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.03))  return false;
		if(!(fTR->ElSigmaIetaIeta[index]<0.03)) return false;
		if(!(fTR->ElHcalOverEcal[index]<0.10)) return false;
	}
  
	// Vertex
	if(!(abs(fTR->ElD0PV[index])<0.02)) return false;
	if(!(abs(fTR->ElDzPV[index])<0.1)) return false;

	// Conversion rejection
  	if(!(fTR->ElPassConversionVeto[index])) return false;
 	if(!(fTR->ElNumberOfMissingInnerHits[index]<=1)) return false;

	// |1/e-1/p|<0.05
	float e=fTR->ElCaloEnergy[index];
	float p=fTR->ElCaloEnergy[index]/fTR->ElESuperClusterOverP[index];
	if(!(fabs(1./e-1./p)<0.05)) return false;
  

  	float pfIso = ElePFIso(index);
  	if ( fabs(fTR->ElEta[index]) < 1.479 || fTR->ElPt[index]>20.0) { // Barrel
		if ( !((pfIso  < 0.15) ) ) return false;
	} else {
    	//Endcap with pt<20
    		if ( !((pfIso  < 0.10) ) ) return false;
  	}
	return true;
}

float MT2Analysis::ElePFIso(const int index){
	// Isolation: EffArea corrected, with cone size 03
	// see here: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolation
	double  neutral = fTR->ElEventelPFIsoValueNeutral03PFIdStandard[index] + fTR->ElEventelPFIsoValueGamma03PFIdStandard[index];
	double  rhocorr = fTR->RhoForIso * EffArea(fTR->ElEta[index]);
	float   iso = ( fTR->ElEventelPFIsoValueCharged03PFIdStandard[index] + TMath::Max(0., neutral - rhocorr) )/ fTR->ElPt[index];
	return iso;
}

const float MT2Analysis::EffArea(float abseta) {
	// see here:  https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection 
	abseta=fabs(abseta); // making sure we're looking at |eta|
	if(abseta<1.0)   return 0.10;
	if(abseta<1.479) return 0.12;
	if(abseta<2.0)   return 0.085;
	if(abseta<2.2)   return 0.11;
	if(abseta<2.3)   return 0.12;
	if(abseta<2.4)   return 0.12;
	return 0.13;
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
	if( fTR->PhoPt[i] < 20                                                 ) return false; // pt cut
	if( fabs(fTR->PhoEta[i])> 2.4                                          ) return false;
	if( fabs(fTR->PhoEta[i])> 1.4442 && fabs(fTR->PhoEta[i])<1.566         ) return false; // veto EB-EE gap
	if( fTR->PhoHoverE[i]   > 0.05                                         ) return false; // H/E cut
	if( fTR->PhoHasPixSeed[i]==1                                           ) return false; // veto pixel seed for electron rejection 
	return true;
}

bool MT2Analysis::IsGoodTau(int i){
       if(fTR->TauPt[i]                          < 15.0         ) return false;    
       if(fabs(fTR->TauEta[i])                   > 2.3          ) return false;
       if(fTR->TauDecayModeFinding[i]            < 0.5          ) return false;  
       if(fTR->TauLooseElectronRejection[i]      < 0.5          ) return false;  
       if(fTR->TauLooseMuonRejection[i]          < 0.5          ) return false;  
       if(fTR->TauVLooseCombinedIsoDBSumPtCorr[i]< 0.5          ) return false;
       return true;
}


// ****************************************************************************************************
// Jet Helper Class

// Jets and JES uncertainty
void MT2Analysis::Initialize_JetCorrectionUncertainty(){
	string PF  ="/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/"+fJEC+"/"+fJEC+ (fisCHSJets?"_Uncertainty_AK5PFchs.txt":"_Uncertainty_AK5PF.txt");
	string Calo="/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/"+fJEC+"/"+fJEC+"_Uncertainty_AK5Calo.txt";

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
	string Calo_L2  ="/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/"+fJEC+"/"+fJEC+"_L2Relative_AK5Calo.txt";
	string Calo_L3  ="/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/"+fJEC+"/"+fJEC+"_L3Absolute_AK5Calo.txt";
	string Calo_RES ="/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/"+fJEC+"/"+fJEC+"_L2L3Residual_AK5Calo.txt";
	string PF_L1    ="/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/"+fJEC+"/"+fJEC+(fisCHSJets?"_L1FastJet_AK5PFchs.txt"    :"_L1FastJet_AK5PF.txt");
	string PF_L2    ="/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/"+fJEC+"/"+fJEC+(fisCHSJets?"_L2Relative_AK5PFchs.txt"   :"_L2Relative_AK5PF.txt");
	string PF_L3    ="/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/"+fJEC+"/"+fJEC+(fisCHSJets?"_L3Absolute_AK5PFchs.txt"   :"_L3Absolute_AK5PF.txt");
	string PF_RES   ="/shome/pnef/MT2Analysis/Code/JetEnergyCorrection/"+fJEC+"/"+fJEC+(fisCHSJets?"_L2L3Residual_AK5PFchs.txt" :"_L2L3Residual_AK5PF.txt");

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


// pfMET 
TLorentzVector MT2Analysis::MET(){
	TLorentzVector MET(0,0,0,0);
	if(!fDoJESUncertainty){
		if(fRemovePhoton && fPhotons.size()==1){
			TLorentzVector trueMET(0,0,0,0), photon(0,0,0,0);
			trueMET.SetPtEtaPhiM((fisType1MET?fTR->PFType1MET:fTR->PFMET), 0., (fisType1MET?fTR->PFType1METphi:fTR->PFMETphi), 0);
			photon .SetPtEtaPhiM(fTR->PhoPt[fPhotons[0]], 0., fTR->PhoPhi[fPhotons[0]],   0);
			MET    = photon + trueMET; 
		}else{
			MET.SetPtEtaPhiM((fisType1MET?fTR->PFType1MET:fTR->PFMET), 0., (fisType1MET?fTR->PFType1METphi:fTR->PFMETphi), 0);
		}
		return MET;
	} else{
		if      (fJESUpDown==1){
    			MET.SetPtEtaPhiM((fisType1MET?fTR->PFType1MET:fTR->PFMET)*1.05, 0., (fisType1MET?fTR->PFType1METphi:fTR->PFMETphi), 0);
		}else if(fJESUpDown==-1){
			MET.SetPtEtaPhiM((fisType1MET?fTR->PFType1MET:fTR->PFMET)*0.95, 0., (fisType1MET?fTR->PFType1METphi:fTR->PFMETphi), 0);
		}else{
			cout << " something wrong in met scaling" << endl; exit(1);
		}		
      	return MET;
	}
}


// MT2AnalysisJet Helper Class ----------------------------------------------------------
MT2AnalysisJet::MT2AnalysisJet(int index, TString type, MT2Analysis *ana): MT2Jet(){
	fAna  =ana;
	fType =type; 
	TLorentzVector jraw(0,0,0,0);
	TLorentzVector jorig(0,0,0,0);
	// PFJets, NO-CHS -------------------------------------------
	if(fType=="PF" && !fAna->fisCHSJets){ 
		jorig.SetPtEtaPhiE(fAna->fTR->JPt[index],               fAna->fTR->JEta[index], fAna->fTR->JPhi[index], fAna->fTR->JE[index]);
		jraw .SetPtEtaPhiM(jorig.Pt()/fAna->fTR->JEcorr[index], jorig.Eta(),            jorig.Phi(),            jorig.M());
		if(fAna->fJEC.length()==0){
			lv    = jorig;
			Scale = fAna->fTR->JEcorr[index];
		}else{                       // REDO JEC
			Scale=        GetPFJEC   (fAna->fTR->JPt[index]/fAna->fTR->JEcorr[index], 
					          fAna->fTR->JEta[index], 
						  fAna->fTR->JArea[index], 
						  fAna->fTR->Rho, 
						  fAna->fisData?1230:123);// L1FastL2L3 + Res(data)
			L1FastJetScale= GetPFJEC (fAna->fTR->JPt[index]/fAna->fTR->JEcorr[index], 
					          fAna->fTR->JEta[index], 
						  fAna->fTR->JArea[index], 
						  fAna->fTR->Rho, 1);               // L1Fast
			lv          = PFJetScaled(jraw, fAna->fTR->JArea[index], fAna->fTR->Rho, fAna->fisData?1230:123);
		}
		// JArea
		Area          =  fAna->fTR->JArea[index];
		// btags
		bTagProbTCHE  =  fAna->fTR->JnewPFTrackCountingHighEffBJetTags[index];
		bTagProbTCHP  =  fAna->fTR->JnewPFTrackCountingHighPurBJetTags[index];
		bTagProbSSVHE =  fAna->fTR->JnewPFSimpleSecondaryVertexHighEffBJetTags[index];
		bTagProbSSVHP =  fAna->fTR->JnewPFSimpleSecondaryVertexHighPurBJetTags[index];
		bTagProbJProb =  fAna->fTR->JnewPFJetProbabilityBPFJetTags[index];
		bTagProbCSV   =  fAna->fTR->JnewPFCombinedSecondaryVertexBPFJetTags[index];
		// PartonFlavour
		if(!fAna->fisData) Flavour = fAna->fTR->JPartonFlavour[index]; // only from V03-08-00 onwards
		// Jet energy fractions
		ChHadFrac     =  fAna->fTR->JChargedHadFrac       [index];	
		NeuHadFrac    =  fAna->fTR->JNeutralHadFrac       [index];
		ChEmFrac      =  fAna->fTR->JChargedEmFrac        [index];
		NeuEmFrac     =  fAna->fTR->JNeutralEmFrac        [index];
		ChMuFrac      =  fAna->fTR->JChargedMuEnergyFrac  [index];
		ChMult        =  fAna->fTR->JNAssoTracks          [index];
		NeuMult       =  fAna->fTR->JNNeutrals            [index];
		NConstituents =  fAna->fTR->JNConstituents        [index];
		isPFIDLoose   =  IsGoodPFJet(10, 5, 1);
		isPFIDMedium  =  IsGoodPFJet(10, 5, 2);
		isPFIDTight   =  IsGoodPFJet(10, 5, 3);

		if(fAna->fVerbose>4) cout << fAna->fTR->Event << " lv " << lv.Pt() << " fTR->JPt[index] " <<  fAna->fTR->JPt[index] << endl;
	} else if (fType == "PF" && fAna->fisCHSJets){
		jorig.SetPtEtaPhiE(fAna->fTR->PFCHSJPt[index],               fAna->fTR->PFCHSJEta[index], fAna->fTR->PFCHSJPhi[index], fAna->fTR->PFCHSJE[index]);
		jraw .SetPtEtaPhiM(jorig.Pt()/fAna->fTR->PFCHSJScale[index], jorig.Eta(),                 jorig.Phi(),                 jorig.M());
		if(fAna->fJEC.length()==0){
			lv             = jorig;
			Scale          = fAna->fTR->PFCHSJScale[index];
			L1FastJetScale = fAna->fTR->PFCHSJL1FastJetScale[index];
		}else{                       // REDO JEC
			Scale=        GetPFJEC   (fAna->fTR->PFCHSJPt[index]/fAna->fTR->PFCHSJScale[index], 
					          fAna->fTR->PFCHSJEta[index], 
						  fAna->fTR->PFCHSJArea[index], 
						  fAna->fTR->Rho, 
						  fAna->fisData?1230:123);// L1FastL2L3 + Res(data)
			L1FastJetScale= GetPFJEC (fAna->fTR->PFCHSJPt[index]/fAna->fTR->PFCHSJScale[index], 
					          fAna->fTR->PFCHSJEta[index], 
						  fAna->fTR->PFCHSJArea[index], 
						  fAna->fTR->Rho, 1);               // L1Fast
			lv          = PFJetScaled(jraw, fAna->fTR->PFCHSJArea[index], fAna->fTR->Rho, fAna->fisData?1230:123);
		}
		// JArea
		Area          =  fAna->fTR->PFCHSJArea[index];
		// btags
		bTagProbTCHE  =  fAna->fTR->PFCHSJtrackCountingHighEffBJetTags[index];
		bTagProbTCHP  =  fAna->fTR->PFCHSJtrackCountingHighPurBJetTags[index];
		bTagProbSSVHE =  fAna->fTR->PFCHSJsimpleSecondaryVertexHighEffBJetTags[index];
		bTagProbSSVHP =  fAna->fTR->PFCHSJsimpleSecondaryVertexHighPurBJetTags[index];
		bTagProbJProb =  fAna->fTR->PFCHSJjetProbabilityBJetTags[index];
		bTagProbCSV   =  fAna->fTR->PFCHSJcombinedSecondaryVertexBJetTags[index];
		// PartonFlavour
		if(!fAna->fisData) Flavour = fAna->fTR->PFCHSJFlavour[index];
		// Jet energy fractions
		ChHadFrac     =  fAna->fTR->PFCHSJChHadfrac       [index];	
		NeuHadFrac    =  fAna->fTR->PFCHSJNeuHadfrac      [index];
		ChEmFrac      =  fAna->fTR->PFCHSJChEmfrac        [index];
		NeuEmFrac     =  fAna->fTR->PFCHSJNeuEmfrac       [index];
		ChMuFrac      =  fAna->fTR->PFCHSJChMufrac        [index];
		ChMult        =  fAna->fTR->PFCHSJChMult          [index];
		NeuMult       =  fAna->fTR->PFCHSJNeuMult         [index];
		NConstituents =  fAna->fTR->PFCHSJNConstituents   [index];
		isPFIDLoose   =  IsGoodPFJet(10, 5, 1);
		isPFIDMedium  =  IsGoodPFJet(10, 5, 2);
		isPFIDTight   =  IsGoodPFJet(10, 5, 3);
		if(fAna->fVerbose>4) cout << fAna->fTR->Event << " lv " << lv.Pt() << " fTR->PFCHSJPt[index] " <<  fAna->fTR->PFCHSJPt[index] << endl;
	}else if (fType == "Calo"){
		jorig.SetPtEtaPhiE(fAna->fTR->CAJPt[index],               fAna->fTR->CAJEta[index],   fAna->fTR->CAJPhi[index],     fAna->fTR->CAJE[index]);
		jraw .SetPtEtaPhiM(jorig.Pt()/fAna->fTR->CAJScale[index], jorig.Eta(),                 jorig.Phi(),                 jorig.M());
		if(fAna->fJEC.length()==0){
			lv    = jorig;
			Scale = fAna->fTR->CAJScale[index];
		}else{                       // REDO JEC
			Scale = GetCaloJEC (fAna->fTR->CAJPt[index]/fAna->fTR->PFCHSJScale[index], 
				            fAna->fTR->CAJEta[index], 
					    fAna->fisData?230:23);// L1FastL2L3 + Res(data)
			lv    = CAJetScaled(jraw, fAna->fisData?230:23);
		}
	}else{
		cout << "ERROR: MT2AnalysisJet::MT2AnalysisJet. exiting.." <<endl;
		exit(-1);
	}
} 
float MT2AnalysisJet::Pt() {return lv.Pt();  }
float MT2AnalysisJet::Eta(){return lv.Eta(); }
float MT2AnalysisJet::Phi(){return lv.Phi(); }
float MT2AnalysisJet::E()  {return lv.E();   }
float MT2AnalysisJet::M()  {return lv.M();   }

TLorentzVector MT2AnalysisJet::PFJetScaled(TLorentzVector jraw, float area, float rho, int level){
	// get new correction factor from DB dumpled txt files
	float scale = GetPFJEC(jraw.Pt(), jraw.Eta(), area, rho, level); 

	TLorentzVector j_scaled(0.,0.,0.,0);
	j_scaled.SetPtEtaPhiM(jraw.Pt()*scale, jraw.Eta(), jraw.Phi(), jraw.M());

	// optionally up or downscale jet: JES uncertainty is a function of corrected Pt
	if(fAna->fJESUpDown== 1) j_scaled.SetPtEtaPhiM(j_scaled.Perp()*(1+GetJECUncertPF(j_scaled.Pt(),j_scaled.Eta())), j_scaled.Eta(), j_scaled.Phi(), j_scaled.M());
	if(fAna->fJESUpDown==-1) j_scaled.SetPtEtaPhiM(j_scaled.Perp()*(1-GetJECUncertPF(j_scaled.Pt(),j_scaled.Eta())), j_scaled.Eta(), j_scaled.Phi(), j_scaled.M());
	return j_scaled;
}

TLorentzVector MT2AnalysisJet::CAJetScaled(TLorentzVector jraw, int level){
	// get new correction factor from DB dumpled txt files
	float scale = GetCaloJEC(jraw.Pt(), jraw.Eta(), level);

	TLorentzVector j_scaled(0.,0.,0.,0);
	j_scaled.SetPtEtaPhiM(jraw.Pt()*scale, jraw.Eta(), jraw.Phi(), jraw.M());

	// optionally up or downscale jet: JES uncertainty is a function of corrected Pt
	if(fAna->fJESUpDown== 1) j_scaled.SetPtEtaPhiM(j_scaled.Perp()*(1+GetJECUncertCalo(j_scaled.Pt(),j_scaled.Eta())), j_scaled.Eta(), j_scaled.Phi(), j_scaled.M());
	if(fAna->fJESUpDown==-1) j_scaled.SetPtEtaPhiM(j_scaled.Perp()*(1-GetJECUncertCalo(j_scaled.Pt(),j_scaled.Eta())), j_scaled.Eta(), j_scaled.Phi(), j_scaled.M());
	return j_scaled;
}

float MT2AnalysisJet::GetJECUncertPF(float pt, float eta){
	// pt must be corrected! not raw  	

	// eta cannot be greater than 5! 
	if      (eta> 5.0) eta = 5.0;
	else if (eta<-5.0) eta =-5.0;

	fAna->fJecUncPF->setJetPt(pt);   
	fAna->fJecUncPF->setJetEta(eta); 
	float uncert= fAna->fJecUncPF->getUncertainty(true);
	return uncert; 
}

float MT2AnalysisJet::GetJECUncertCalo(float pt, float eta){
	// pt must be corrected! not raw  	
	
	// eta cannot be greater than 5! 
	if      (eta> 5.0) eta = 5.0;
	else if (eta<-5.0) eta =-5.0;
	fAna->fJecUncCalo->setJetPt(pt);
	fAna->fJecUncCalo->setJetEta(eta);
	float uncert= fAna->fJecUncCalo->getUncertainty(true);
	return uncert; 
}

float MT2AnalysisJet::GetPFJEC(float rawpt, float eta, float area, float rho, int level){

	fAna->fJetCorrectorPF->setJetPt(rawpt); // WARNING: JEC is a function of RAW pt
	fAna->fJetCorrectorPF->setJetEta(eta);
	fAna->fJetCorrectorPF->setJetA(area);
	fAna->fJetCorrectorPF->setRho(rho);
	
	vector<float> factors = fAna->fJetCorrectorPF->getSubCorrections();
	// convention for correction factors as follows: 0->L1, 1->L1L2, 2->L1L2L3, 3->L1L2L3Res
	// see MT2Analysis::Initialize_JetEnergyCorrection()
	
	if( !(fAna->fisData) && level==1230){
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

float MT2AnalysisJet::GetCaloJEC(float rawpt, float eta, int level){
	fAna->fJetCorrectorCalo->setJetEta(eta);
	fAna->fJetCorrectorCalo->setJetPt(rawpt); // WARNING: JEC is a function of RAW pt
	
	vector<float> factors = fAna->fJetCorrectorCalo->getSubCorrections();
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




