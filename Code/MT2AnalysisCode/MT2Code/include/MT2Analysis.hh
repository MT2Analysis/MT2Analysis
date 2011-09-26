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
#include "MT2tree.hh"
#include <TLorentzVector.h>
#include "JetCorrectionUncertainty.h" 

static const int gNHemispheres = 4;
static const int gNGenJets     = 20;

class MT2Analysis: public UserAnalysisBase{
public:
	MT2Analysis(TreeReader *tr = NULL);
	virtual ~MT2Analysis();

	void Begin(const char* filename = "Mass_histos.root");
	void Analyze();
	void End();
	void SetType(bool isData=false){
		fisData=isData;
	};
	// parser
	void ReadCuts(const char* SetofCuts);

  	//  reweight::LumiReWeighting *LumiReW;
  	bool isS3;
  	bool noPU;
private:
	// files and trees ----------------------------------------------------------------
	// file for histograms:
	TFile* fHistFile;

        // MT2 mini-tree
        MT2tree* fMT2tree;
        TTree* fATree;

	// data members-----------------------------------------
	bool fisData;
	int   fCounter;
	bool fBasicMT2treeFilled;

	// vectors for jets and lepton indices
	vector<int> fElecs;
	vector<int> fTaus;
	vector<int> fMuons;
	vector<int> fJets;
	
	// MT2 and hemisphere
	Davismt2 *fMT2;
	TMctLib  *fMCT;
	Hemisphere *fHemisphere;


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


	// ---- required and vetoed triggers ----
	std::vector<std::string> fRequiredHLT; 
	std::vector<std::string> fVetoedHLT;
	
	// ---- file for EcalDeadCell veto
	std::vector<std::string> fTPfiles;
	std::vector<std::string> fBEfiles;
	
	typedef std::map <string, bool*> StringBoolMap;
	StringBoolMap fTriggerMap;
	
	// structs
	enum LeptConfig {
	 	e, mu, OS_emu, OS_ee, OS_mumu, SS_emu, SS_ee, SS_mumu, multilept, null
  	};
	LeptConfig fLeptConfig;
	
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
	void FindLeptonConfig();
	void GetLeptonJetIndices();
	bool IsSelectedEvent();
	void InitializeEvent();
	void DeadCellParser(DeadCellFilter &DeadCellFilter_, string file_);
	//pfJetID
	bool IsGoodBasicPFJetPAT3(int index, float ptcut, float absetacut);
	bool IsGoodPFJetMediumPAT3(int index, float ptcut, float absetacut);
	bool IsGoodPFJetTightPAT3(int index, float ptcut, float absetacut);
	// jets
	TLorentzVector Jet(int index);
	TLorentzVector CAJet(int index);
	TLorentzVector PFJetJESScaled(TLorentzVector j);
	TLorentzVector CAJetJESScaled(TLorentzVector j);
	TLorentzVector MET();
	float GetJECUncertPF(float pt, float eta);
	float GetJECUncertCalo(float pt, float eta);
	void Initialize_JetCorrectionUncertainty();
	JetCorrectionUncertainty *fJecUncPF;   
	JetCorrectionUncertainty *fJecUncCalo; 

};
#endif
