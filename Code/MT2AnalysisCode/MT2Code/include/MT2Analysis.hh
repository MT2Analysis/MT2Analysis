#ifndef MT2Analysis_hh
#define MT2Analysis_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <numeric>

#include "base/UserAnalysisBase.hh"
#include "base/TreeReader.hh"
#include "helper/Davismt2.h"
#include "helper/TMctLib.h"
#include "helper/Hemisphere.hh"
#include "BTagWeight.hh"
#include "MT2tree.hh"
#include <TLorentzVector.h>
#include "JetCorrectionUncertainty.h" 
#include "JetCorrectorParameters.h" 
#include "FactorizedJetCorrector.h" 

#include "TH2F.h"
#include "TTree.h"

// #include </shome/leo/Installations/LHAPDF/lhapdf-5.8.4/include/LHAPDF/LHAPDF.h>


static const int gNHemispheres = 4;
static const int gNGenJets     = 20;

class MT2Analysis: public UserAnalysisBase{
public:
	MT2Analysis(TreeReader *tr = NULL);
	virtual ~MT2Analysis();

	void Begin(const char* filename = "Mass_histos.root");
	void Analyze();
	void End();
	void SetType(bool isData=false){fisData=isData;};
	void SetProcessID(int ID){fID=ID;};
	void SetBTagEfficiency(string btagFileName){ fbtagFileName = btagFileName;};


	// redo JEC
	void SetJEC(string JEC){fJEC=JEC;};
	// parser
	void ReadCuts(const char* SetofCuts);

	// PU
  	enum PUScenario {noPU, MC2012};
	PUScenario fPUScenario;
	void SetPUReweighting(string PU, string data_PileUp, string mc_PileUp){
		if      (PU =="MC2012") {fPUScenario=MC2012;   SetPileUpSrc  (data_PileUp, mc_PileUp);}
		else                     fPUScenario=noPU;
	};

 	 //create pdf weights
        bool doPDF;

  	//is a susy scan?
  	bool isScan;

  	// remove Photon
  	bool fRemovePhoton;
  
  	//Control histos
  	TH1F *fH_PUWeights, *fH_Events ;
  	TH2F *fH2_SMSEvents, *fH2_mSugraEvents;
  	//SUSY subprocess histos
  	TH2F *fH_mSugraSubProcEvents[11];


private:
	// files and trees ----------------------------------------------------------------
	// file for histograms:
	TFile* fHistFile;

        // MT2 mini-tree
        MT2tree* fMT2tree;
        TTree* fATree;

        //Control Histograms
  //  TH1F *fH_PUWeights ;

	//btagging histograms and files
	TFile *btagfile;
	TH1D *hbeff;
	TH1D *hceff;
	TH1D *hleff;

	// data members-----------------------------------------
	bool fisData;
	int  fID;
	int  fCounter;
	bool fBasicMT2treeFilled;
	string fJEC;
	string fbtagFileName;


	// vectors for jets and lepton indices
	vector<int> fElecs;
	vector<int> fPhotons;
	vector<int> fTaus;
	vector<int> fMuons;
	vector<int> fJets;
	vector<bool>fPhotonJetOverlapRemoved;
	
	// MT2 and hemisphere
	Davismt2 *fMT2;
	TMctLib  *fMCT;
	Hemisphere *fHemisphere;

  //PDFs
  int nPDFs;

	// cut variables
	float fHT;
	float fCaloHT50;
	float fCaloHT50_ID;
	float fCaloMHT30;
	float fCaloMHT30_ID;
	bool  fCrazyHCAL;
	bool  fIsNANObj;
	bool  fNegativeJEC;

	
	//  ---- set of cuts ---
	TString fSetName;
	float fCut_PFMET_min;
	float fCut_HT_min;
	float fCut_caloHT50_min;
	float fCut_caloHT50ID_min;
	float fCut_caloMHT30_min;
	float fCut_caloMHT30ID_min;
	float fCut_JPt_hardest_min;
	float fCut_JPt_second_min;
	float fCut_PtHat_max;
	int   fCut_Run_min;
	int   fCut_Run_max;
	bool  fDoJESUncertainty;
	int   fJESUpDown;
  int   fCut_NJets40_min;


	// ---- required and vetoed triggers ----
	std::vector<std::string> fRequiredHLT; 
	std::vector<std::string> fVetoedHLT;
	
	// ---- file for EcalDeadCell veto
	std::vector<std::string> fTPfiles;
	std::vector<std::string> fBEfiles;
	
	typedef std::map <string, bool*> StringBoolMap;
	StringBoolMap fTriggerMap;
	
	struct DeadCellFilter{
		vector<int> run;
		vector<int> lumi;
		vector<int> event;
		void Reset(){
			run.clear();
			lumi.clear();
			event.clear();
		}
	} fDeadCellFilterBE, fDeadCellFilterTP;



	// member functions -------------------------------------------------------------
	void BookTree();
	void FillTree();
	void ResetTree(); 
	bool FillMT2TreeBasics();
	void FillMT2treeCalculations();
	void GetLeptonJetIndices();
	bool IsSelectedEvent();
	void InitializeEvent();
	void DeadCellParser(DeadCellFilter &DeadCellFilter_, string file_);
	// photons
	bool IsGoodPhotonEGMLooseISO(int index);
	bool IsGoodPhotonEGMLooseRelISO(int index);
	bool IsGoodPhotonEGMLooseID(int index);
	bool IsGoodPhotonEGMTightISO(int index);
	bool IsGoodPhotonEGMTightID(int index);
	bool IsGoodPhoton(int index);

        //Taus
        bool IsGoodTau(int index);

	//Muons 
	bool IsGoodMT2Muon(const int index);
	float MuPFIso(const int index);
	// Electrons
	bool IsGoodMT2ElectronVetoID(const int index);
	bool IsGoodMT2ElectronMediumID(const int index);
	const float EffArea(float abseta);
	float ElePFIso(const int index);

	//pfJetID
	bool IsGoodMT2PFJetIDLoose(int index, float ptcut, float absetacut);
	bool IsGoodMT2PFJetIDMedium(int index, float ptcut, float absetacut);
	bool IsGoodMT2PFJetIDTight(int index, float ptcut, float absetacut);
	// jets
	TLorentzVector Jet(int index);
	TLorentzVector CAJet(int index);
	TLorentzVector PFJetScaled(TLorentzVector j, float old_scale, float area, float rho);
	TLorentzVector CAJetScaled(TLorentzVector j, float old_scale);
	TLorentzVector MET();
	float GetJECUncertPF(float pt, float eta);
	float GetJECUncertCalo(float pt, float eta);
	float GetPFJEC  (float corrpt, float scale, float eta, float area, float rho, int level);
	float GetCaloJEC(float corrpt, float scale, float eta, int level);
	void Initialize_JetCorrectionUncertainty();
	void Initialize_JetEnergyCorrection();
	JetCorrectionUncertainty *fJecUncPF;   
	JetCorrectionUncertainty *fJecUncCalo; 

	JetCorrectorParameters* fJecCaloL2;
	JetCorrectorParameters* fJecCaloL3; 
	JetCorrectorParameters* fJecCaloRES; 
	JetCorrectorParameters* fJecPFL1;   
	JetCorrectorParameters* fJecPFL2;   
	JetCorrectorParameters* fJecPFL3;   
	JetCorrectorParameters* fJecPFRES;  

	FactorizedJetCorrector* fJetCorrectorPF;
	FactorizedJetCorrector* fJetCorrectorCalo;

};
#endif
