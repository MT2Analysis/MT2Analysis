#ifndef MT2tree_hh
#define MT2tree_hh

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"

enum {m_jetSize = 25, m_genjetSize = 20,  m_tauSize = 8, m_eleSize = 8, m_muoSize = 8, m_phoSize = 8, m_genleptSize=20, m_hemiSize=2};



// MT2Misc ----------------------------------
class MT2Misc : public TObject {

public:
  MT2Misc();
  virtual ~MT2Misc();

  void Reset();
  
  Bool_t   isData;
  Bool_t   isCHSJets;
  Bool_t   isType1MET;

  Bool_t   isFastSim;
  Int_t    Run;
  Int_t    Event;
  Int_t    LumiSection;
  Int_t    ProcessID;             // 0=data, 1=Znunu, 2=Zll, 3=WJets, 4=ttbar, 5=Gamma+jets, 6=QCD ,[empty], 9=Other, 10=Signal
  Int_t    Jet0Pass;
  Int_t    Jet1Pass;
  Int_t    PassJetID;
  Int_t    PassJet40ID;
  Int_t    PassJet30ID;//can probably be deleted soon, again? (common object definition hadronic analyses)
  Float_t  MT2;
  Float_t  MT2jet40;
  Float_t  MCT;
  Float_t  MET;
  Float_t  METPhi;
  Float_t  CaloMETRaw;
  Float_t  CaloMETMuJesCorr;

  Float_t  LeadingJPt;
  Float_t  SecondJPt;
  Float_t  J3Pt;
  Float_t  J4Pt;
  Float_t  Vectorsumpt;
  Float_t  MinMetJetDPhi;       // use all jets
  Float_t  MinMetJetDPhi4;      // use first 4 jets
  Float_t  MinMetJetDPhiPt40;   // use jets that have pT>40GeV
  Float_t  MinMetJetDPhi4Pt40;  // use first 4 jets that have pT>40GeV
  Int_t    MinMetJetDPhiIndex;
  Float_t  MinMetBJetDPhi;
  Float_t  QCDPartonicHT;
  Float_t  HT;
  Float_t  pfHT30;//do we need all pfHTs
  Float_t  pfHT35;
  Float_t  pfHT40;
  Float_t  pfHT45;
  Float_t  pfHT50;
  Int_t    WDecayMode;
  Int_t    TopDecayMode;

  //Noise Filters
  Bool_t   CrazyHCAL;             // store bad events as "true"
  Bool_t   NegativeJEC;           
  Bool_t   CSCTightHaloIDFlag;    
  Bool_t   HBHENoiseFlag;
  Bool_t   hcalLaserEventFlag;
  Bool_t   trackingFailureFlag;
  Bool_t   eeBadScFlag;
  Bool_t   EcalDeadCellTriggerPrimitiveFlag;
  Bool_t   EcalLaserCorrFlag;
  Bool_t   TrackingManyStripClusFlag;
  Bool_t   TrackingTooManyStripClusFlag;
  Bool_t   TrackingLogErrorTooManyClustersFlag;
  
  ClassDef(MT2Misc, 37)
};



// ----------------------------------------
class MT2Susy : public TObject {
  
public:
  MT2Susy();
  virtual ~MT2Susy();
  void Reset();
  
  float MassGlu;
  float MassChi;
  float MassLSP;
  float M0;
  float M12;
  float A0;
  float Mu;
  float XSec;
  float TanBeta;

  ClassDef(MT2Susy, 1);
};


// ----------------------------------------
class MT2PileUp : public TObject {

public:
	MT2PileUp();
	virtual ~MT2PileUp();
	void Reset();

	Int_t    PUnumInt;
	Int_t    PUnumIntEarly;
	Int_t    PUnumIntLate;
	Int_t	 PUtrueNumInt;
	Int_t    PUScenario;
	Float_t  PtHat;
	Float_t  Weight;
  	Int_t    NVertices;  // good reco vertices
	Float_t  Rho;


	ClassDef(MT2PileUp, 6);
};

// --------------------------------
class MT2Trigger : public TObject {

public:
	MT2Trigger();
	virtual ~MT2Trigger();
	void Reset();

        // 2012 datasets

	// HT/jetHT dataset
	Bool_t HLT_HT500_v1;
	Bool_t HLT_HT500_v2;
	Bool_t HLT_HT500_v3;
	Bool_t HLT_HT500_v4;
	Bool_t HLT_HT550_v1;
	Bool_t HLT_HT550_v2;
	Bool_t HLT_HT550_v3;
	Bool_t HLT_HT550_v4;
	Bool_t HLT_HT600_v1;
	Bool_t HLT_HT600_v2;
	Bool_t HLT_HT600_v3;
	Bool_t HLT_HT600_v4;
	Bool_t HLT_HT650_v1;
	Bool_t HLT_HT650_v2;
	Bool_t HLT_HT650_v3;
	Bool_t HLT_HT650_v4;
	Bool_t HLT_PFHT350_v3;
	Bool_t HLT_PFHT350_v4;
	Bool_t HLT_PFHT350_v5;
	Bool_t HLT_PFHT350_v6;
	Bool_t HLT_PFHT350_v7;
	Bool_t HLT_PFHT650_v5;
	Bool_t HLT_PFHT650_v6;
	Bool_t HLT_PFHT650_v7;
	Bool_t HLT_PFHT650_v8;
	Bool_t HLT_PFHT650_v9;
	Bool_t HLT_PFHT700_v3;
	Bool_t HLT_PFHT700_v4;
	Bool_t HLT_PFHT700_v5;
	Bool_t HLT_PFHT700_v6;
	Bool_t HLT_PFHT700_v7;
	Bool_t HLT_PFHT750_v3;
	Bool_t HLT_PFHT750_v4;
	Bool_t HLT_PFHT750_v5;
	Bool_t HLT_PFHT750_v6;
	Bool_t HLT_PFHT750_v7;
	Bool_t HLT_PFNoPUHT350_v1;
	Bool_t HLT_PFNoPUHT350_v2;
	Bool_t HLT_PFNoPUHT350_v3;
	Bool_t HLT_PFNoPUHT350_v4;
	Bool_t HLT_PFNoPUHT650_v1;
	Bool_t HLT_PFNoPUHT650_v2;
	Bool_t HLT_PFNoPUHT650_v3;
	Bool_t HLT_PFNoPUHT650_v4;
	Bool_t HLT_PFNoPUHT700_v1;
	Bool_t HLT_PFNoPUHT700_v2;
	Bool_t HLT_PFNoPUHT700_v3;
	Bool_t HLT_PFNoPUHT700_v4;
	Bool_t HLT_PFNoPUHT750_v1;
	Bool_t HLT_PFNoPUHT750_v2;
	Bool_t HLT_PFNoPUHT750_v3;
	Bool_t HLT_PFNoPUHT750_v4;

	// HT/HTMHT dataset
	Bool_t HLT_PFHT350_PFMET100_v3;
	Bool_t HLT_PFHT350_PFMET100_v4;
	Bool_t HLT_PFHT350_PFMET100_v5;
	Bool_t HLT_PFHT350_PFMET100_v6;
	Bool_t HLT_PFHT350_PFMET100_v7;
	Bool_t HLT_PFHT400_PFMET100_v3;
	Bool_t HLT_PFHT400_PFMET100_v4;
	Bool_t HLT_PFHT400_PFMET100_v5;
	Bool_t HLT_PFHT400_PFMET100_v6;
	Bool_t HLT_PFHT400_PFMET100_v7;
	Bool_t HLT_PFNoPUHT350_PFMET100_v1;
	Bool_t HLT_PFNoPUHT350_PFMET100_v3;
	Bool_t HLT_PFNoPUHT350_PFMET100_v4;
	Bool_t HLT_PFNoPUHT400_PFMET100_v1;
	Bool_t HLT_PFNoPUHT400_PFMET100_v3;
	Bool_t HLT_PFNoPUHT400_PFMET100_v4;

        // MET Dataset
        Bool_t HLT_MET120_v9;
        Bool_t HLT_MET120_v10;
        Bool_t HLT_MET200_v9;
        Bool_t HLT_MET200_v10;
        Bool_t HLT_MET120_HBHENoiseCleaned_v2;
        Bool_t HLT_MET120_HBHENoiseCleaned_v3;
        Bool_t HLT_MET120_HBHENoiseCleaned_v4;
        Bool_t HLT_MET200_HBHENoiseCleaned_v2;
        Bool_t HLT_MET200_HBHENoiseCleaned_v3;
        Bool_t HLT_MET200_HBHENoiseCleaned_v4;
        Bool_t HLT_PFMET150_v2;
        Bool_t HLT_PFMET150_v3;
        Bool_t HLT_PFMET150_v4;
        Bool_t HLT_PFMET150_v5;
        Bool_t HLT_PFMET150_v6;
        Bool_t HLT_PFMET150_v7;
        Bool_t HLT_PFMET180_v2;
        Bool_t HLT_PFMET180_v3;
        Bool_t HLT_PFMET180_v4;
        Bool_t HLT_PFMET180_v5;
        Bool_t HLT_PFMET180_v6;
        Bool_t HLT_PFMET180_v7;
        Bool_t HLT_DiCentralPFJet30_PFMHT80_v5;
        Bool_t HLT_DiCentralPFJet30_PFMHT80_v6;
        Bool_t HLT_DiCentralPFJet30_PFMHT80_v7;
        Bool_t HLT_DiCentralPFJet50_PFMET80_v3;
        Bool_t HLT_DiCentralPFJet50_PFMET80_v4;
        Bool_t HLT_DiCentralPFJet50_PFMET80_v5;
        Bool_t HLT_DiCentralPFJet50_PFMET80_v6;
        Bool_t HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v1;
        Bool_t HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v3;
        Bool_t HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4;

        // Multijet dataset
        Bool_t HLT_DiJet80_DiJet60_DiJet20_v1;
        Bool_t HLT_DiJet80_DiJet60_DiJet20_v2;
        Bool_t HLT_DiJet80_DiJet60_DiJet20_v3;
        Bool_t HLT_DiJet80_DiJet60_DiJet20_v4;
        Bool_t HLT_DiJet80_DiJet60_DiJet20_v5;
        Bool_t HLT_DiJet80_DiJet60_DiJet20_v6;
        Bool_t HLT_QuadJet60_DiJet20_v1;
        Bool_t HLT_QuadJet60_DiJet20_v2;
        Bool_t HLT_QuadJet60_DiJet20_v3;
        Bool_t HLT_QuadJet60_DiJet20_v5;
        Bool_t HLT_QuadJet60_DiJet20_v6;
        Bool_t HLT_QuadJet50_v1;
        Bool_t HLT_QuadJet50_v2;
        Bool_t HLT_QuadJet50_v3;
        Bool_t HLT_QuadJet50_v5;
        Bool_t HLT_QuadJet70_v1;
        Bool_t HLT_QuadJet70_v2;
        Bool_t HLT_QuadJet70_v3;
        Bool_t HLT_QuadJet70_v4;
        Bool_t HLT_QuadJet70_v6;
        Bool_t HLT_QuadJet80_v1;
        Bool_t HLT_QuadJet80_v2;
        Bool_t HLT_QuadJet80_v3;
        Bool_t HLT_QuadJet80_v4;
        Bool_t HLT_QuadJet80_v6;
        Bool_t HLT_SixJet35_v1;
        Bool_t HLT_SixJet35_v2;
        Bool_t HLT_SixJet35_v3;
        Bool_t HLT_SixJet35_v4;
        Bool_t HLT_SixJet35_v6;
        Bool_t HLT_SixJet45_v1;
        Bool_t HLT_SixJet45_v2;
        Bool_t HLT_SixJet45_v3;
        Bool_t HLT_SixJet45_v4;
        Bool_t HLT_SixJet45_v6;

        // Jet/JetMon dataset
        Bool_t HLT_PFJet320_v3; //JetHT
        Bool_t HLT_PFJet320_v4; //JetHT
        Bool_t HLT_PFJet320_v5; //JetHT
        Bool_t HLT_PFJet260_v3;
        Bool_t HLT_PFJet260_v4;
        Bool_t HLT_PFJet260_v5;
        Bool_t HLT_PFJet200_v3;
        Bool_t HLT_PFJet200_v4;
        Bool_t HLT_PFJet200_v5;
        Bool_t HLT_PFJet140_v3;
        Bool_t HLT_PFJet140_v4;
        Bool_t HLT_PFJet140_v5;
        Bool_t HLT_DiPFJetAve200_v3;
        Bool_t HLT_DiPFJetAve200_v4;
        Bool_t HLT_DiPFJetAve200_v5;
        Bool_t HLT_DiPFJetAve200_v6;
        Bool_t HLT_DiPFJetAve140_v3;
        Bool_t HLT_DiPFJetAve140_v4;
        Bool_t HLT_DiPFJetAve140_v5;
        Bool_t HLT_DiPFJetAve140_v6;
        Bool_t HLT_DiPFJetAve80_v3;
        Bool_t HLT_DiPFJetAve80_v4;
        Bool_t HLT_DiPFJetAve80_v5;
        Bool_t HLT_DiPFJetAve80_v6;
        Bool_t HLT_DiPFJetAve40_v3;
        Bool_t HLT_DiPFJetAve40_v4;
        Bool_t HLT_DiPFJetAve40_v5;
        Bool_t HLT_DiPFJetAve40_v6;

	// Photon
	Bool_t HLT_SinglePhotons;
	Bool_t HLT_SinglePhoton70_HT400;
	Bool_t HLT_SinglePhoton70_MET100;
	// Dileptons
	Bool_t HLT_DiElectrons;
	Bool_t HLT_DiMuons;
	Bool_t HLT_EMu;
	// Single Leptons
	// please look at exact triggers (some are pfnopu, others not)
	Bool_t HLT_SingleMu;
	Bool_t HLT_SingleMu_Jet;
	Bool_t HLT_SingleMu_DiJet; // only from run 193834 on
	Bool_t HLT_SingleMu_TriJet;
	Bool_t HLT_SingleEle_DiJet_MET;
	Bool_t HLT_SingleEle_MET_MT;

	ClassDef(MT2Trigger, 22);
};

// MT2Jet ----------------------------------
class MT2Jet : public TObject {

public:
  MT2Jet();
  virtual ~MT2Jet();

  void Reset();
  void SetLV(const TLorentzVector v);
  Bool_t IsGoodPFJet(float minJPt=20., float maxJEta=2.4, int PFJID=1); // PFJID: 1 - loose, 2 - medium, 3 - tight
  Bool_t IsBJet(Int_t algo=4); // algo 4 = CSV, algo 3 = SSVHP, algo 2 = SSVHE
  Bool_t IsBJet(Int_t algo, Int_t WP);
  TLorentzVector lv;

  Float_t bTagProbTCHE;
  Float_t bTagProbTCHP;
  Float_t bTagProbSSVHE;
  Float_t bTagProbSSVHP;
  Float_t bTagProbJProb;
  Float_t bTagProbCSV;

  Bool_t isPFIDLoose;
  Bool_t isPFIDMedium;
  Bool_t isPFIDTight;

  Float_t ChHadFrac;
  Float_t NeuHadFrac;
  Float_t ChEmFrac;
  Float_t NeuEmFrac;
  Float_t ChMuFrac;
  Int_t   ChMult;
  Int_t   NeuMult;
  Int_t   NConstituents;

  Float_t Scale;          // scale factor from JE correction
  Float_t L1FastJetScale; // correction factor from raw to L1FastJetcorrected
  Float_t Area;

  Int_t   Flavour;   // JetFlavour for MC
  
  //  Bool_t  isTauMatch; // tells you if pf-jet is matched to a tau
  Int_t  isTauMatch; // tells you if pf-jet is matched to a tau, >= 0 gives the tau index in the tau collection, -1 not matched
  Float_t TauDR;
  Float_t TauDPt;
  Int_t   NTauMatch;

  ClassDef(MT2Jet, 15)
};

// MT2GenJet -------------------------
class MT2GenJet : public TObject {

public:
  MT2GenJet();
  virtual ~MT2GenJet();

  void Reset();
  TLorentzVector lv;
  
  Int_t   JetMatchIndex;
  Float_t DeltaR;

  ClassDef(MT2GenJet, 2)
};


// MT2Hemi ---------------------------
class MT2Hemi : public TObject {

public:
  MT2Hemi();
  virtual ~MT2Hemi();

  void Reset();
  Int_t         seed_method;
  Int_t         assoc_method;
  Int_t         NHemiIter;
  Float_t       MT2;
  Float_t       MCT;
  Float_t       AlphaT;
  Float_t       minDHT;
  Float_t       maxDR;
  Float_t       dPhi;

  Int_t          jindices1  [m_jetSize];
  Int_t          jindices2  [m_jetSize];
  Int_t          eleindices1[m_eleSize];
  Int_t          eleindices2[m_eleSize];
  Int_t          muoindices1[m_muoSize];
  Int_t          muoindices2[m_muoSize];
  TLorentzVector lv1;
  TLorentzVector lv2;
  TLorentzVector UTM;

  ClassDef(MT2Hemi, 6)
};

// MT2Tau  ----------------------------------
class MT2Tau : public TObject {

public:
  MT2Tau();
  virtual ~MT2Tau();

  void Reset();
  void SetLV(const TLorentzVector v);
  Bool_t IsGoodTau(float minJPt, float maxJEta, int iso, int elerej, int muorej);

  TLorentzVector lv;

  Float_t  MT;
  Int_t    Charge;
  Float_t  JetPt;
  Float_t  JetEta;
  Float_t  JetPhi;
  Float_t  JetMass;
  Int_t    Isolation;
  Int_t    ElectronRej;
  Int_t    MuonRej;
  Int_t    Isolation3Hits;
  Int_t    IsolationMVA2;
  Int_t    ElectronRejMVA3;
  Int_t    MuonRej2;
  Bool_t   isLooseID;
  Bool_t   isLooseID3Hits;
  Bool_t   isLooseIDMVA;
  ClassDef(MT2Tau, 5)
};

// MT2Elec ----------------------------------
class MT2Elec : public TObject {

public:
  MT2Elec();
  virtual ~MT2Elec();

  void Reset();
  void SetLV(const TLorentzVector v);

  TLorentzVector lv;

  Float_t  MT;
  Float_t  Iso;
  Float_t  Iso04;
  Int_t    Charge;
  Int_t    IDMedium;
  Int_t    IDLoose;
  Int_t    IDVeto;

  ClassDef(MT2Elec, 12)
};

// MT2Muon ----------------------------------
class MT2Muon : public TObject {

public:
  MT2Muon();
  virtual ~MT2Muon();

  void Reset();
  void SetLV(const TLorentzVector v);

  TLorentzVector lv;

  Float_t  MT;
  Float_t  Iso;
  Float_t  Iso04;
  Int_t    Charge;

  ClassDef(MT2Muon, 9)
};

// MT2Photon  ----------------------------------
class MT2Photon : public TObject {

public:
  MT2Photon();
  virtual ~MT2Photon();

  void Reset();
  void SetLV(const TLorentzVector v);
  Bool_t IsIsolated(int i);
  Bool_t IsID(int i, float minpt, float maxeta);

  TLorentzVector lv;
  Float_t TrkIso;
  Float_t EcalIso; 
  Float_t HcalIso;
  Float_t SigmaIEtaIEta;
  Float_t HoverE;
  Float_t HoverE2012;
  Float_t ChargedHadIso;
  Float_t NeutralHadIso;
  Float_t PhotonIso;
  Int_t   MCmatchexitcode;
  Float_t GenJetMinDR;
  Bool_t  JetRemoved;
  Bool_t  isLooseID;//contains Iso
  Bool_t  isMediumID;//contains Iso
  Bool_t  isTightID;//contains Iso
  Bool_t  isLooseIso;
  Bool_t  isMediumIso;
  Bool_t  isTightIso;

  ClassDef(MT2Photon, 8)
};


// MT2GenLept ----------------------------------
class MT2GenLept : public TObject {

public:
  MT2GenLept();
  virtual ~MT2GenLept();

  void Reset();

  TLorentzVector lv;
  TLorentzVector Mlv;
  TLorentzVector GMlv;
  Int_t          ID;
  Int_t          MID;
  Int_t          MStatus;
  Int_t          GMID;
  Int_t          GMStatus;
  Float_t        MT;


  ClassDef(MT2GenLept, 6)
};

// MT2SFWeight -----------------------------
class MT2SFWeight : public TObject {

public:
  MT2SFWeight();
  virtual ~MT2SFWeight();

  void Reset();

  Float_t BTagCSV40eq0;
  Float_t BTagCSV40eq1;
  Float_t BTagCSV40eq2;
  Float_t BTagCSV40eq3;
  Float_t BTagCSV40ge1;
  Float_t BTagCSV40ge2;
  Float_t BTagCSV40ge3;
  Float_t BTagCSV40eq0Error;
  Float_t BTagCSV40eq1Error;
  Float_t BTagCSV40eq2Error;
  Float_t BTagCSV40eq3Error;
  Float_t BTagCSV40ge1Error;
  Float_t BTagCSV40ge2Error;
  Float_t BTagCSV40ge3Error;
  Float_t TauTageq0;
  Float_t TauTageq1;
  Float_t TauTagge1;
  Float_t TauTageq0Error;
  Float_t TauTageq1Error;
  Float_t TauTagge1Error;


  ClassDef(MT2SFWeight, 1);
};

// MT2tree ----------------------------------
class MT2tree : public TObject {

public:
  MT2tree();
  virtual ~MT2tree();

  void Reset();

  void SetNJets            (int n);
  void SetNGenJets         (int n);
  void SetNJetsIDLoose     (int n);
  void SetNBJetsCSVM       (int n);
  void SetNBJetsCSVT       (int n);
  void SetNBJets40CSVM     (int n);
  void SetNBJets40CSVT     (int n);
  void SetNEles            (int n);
  void SetNMuons           (int n);
  void SetNMuonsCommonIso  (int n);
  void SetNPhotons         (int n);
  void SetNPhotonsIDLoose  (int n);
  void SetNPhotonsIDLoose25(int n);
  void SetNTaus            (int n);
  void SetNTausIDLoose     (int n);
  void SetNTausIDLoose3Hits(int n);
  void SetNTausIDLooseMVA  (int n);
  void SetNTausIDLoose2    (int n);
  
  // My functions here
  // NJets
  Int_t    GetNjets   (float minJPt=20, float maxJEta=5., int PFJID=0);  // PFJETID not depends on pt and eta
  Int_t    GetJetIndex(int ijet=0, int PFJID=0, float minJPt=20, float maxJEta=2.4);
  Int_t    GetJetIndexByEta(int ijet=0, int PFJID=0, float minJPt=20, float maxJEta=2.4);
  Int_t    GetNBtags  (int algo=3, float value=2., float minJPt=20, float maxJEta=2.4, int PFJID=0);  // algo - 0:TCHE, 1:TCHP, 2:SSVHE, 3:SSVHP
  Float_t  JetPt      (int ijet=0, int PFJID=0, float minJPt=20, float maxJEta=2.4);
  Int_t    GetNTaus   (float minJPt=20, float maxJEta=2.4, int iso=2, int elerej=1, int muorej=3);

  // HT, MHT, ...
  Float_t GetHT         (int PFJID=0, float minJPt=50, float maxJEta=2.4);
  TLorentzVector GetMHTlv(int PFJID=0, float minJPt=20, float maxJEta=2.4, bool inclLepts=false);
  Float_t GetMHT        (int PFJID=0, float minJPt=20, float maxJEta=2.4, bool inclLepts=false);
  Float_t GetMHTPhi     (int PFJID=0, float minJPt=20, float maxJEta=2.4, bool inclLepts=false);
  Float_t GetMHTminusMET(int PFJID=0, float minJPt=20, float maxJEta=2.4, bool inclLepts=false);
  // dPhi and friends
  Bool_t  PassJetID(float minJPt=50, float maxJEta=5.0, int PFJID=1);
  Float_t JetsDPhi(int j1=1, int j2=0, int PFJID=0);
  Float_t JetsInvMass(int j1=0, int j2=1);
  Float_t MetJetDPhi(int ijet = 0, int PFJID=0, int met=1);
  Bool_t  PassMinMetJetDPhi03();
  Float_t GetMinR12R21      (int PFJID=0, float minJPt=20, float maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Float_t MinMetJetDPhi     (int PFJID=0, float minJPt=20, float maxJEta=6., int met=1, int njets=0, bool onlyCrack=false); // electrons and muons not considered for minDPhi, njets=0->all jets
  Int_t   MinMetJetDPhiIndex(int PFJID=0, float minJPt=20, float maxJEta=6., int met=1, int njets=0, bool onlyCrack=false); // electrons and muons not considered for minDPhi, njets=0->all jets
  Int_t   BiasedDPhiIndex   (int PFJID, float minJPt, float maxJEta);
  Float_t BiasedDPhi        (int PFJID, float minJPt, float maxJEta);
  Float_t MaxMetJetDPhi     (int PFJID=0, float minJPt=20, float maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Int_t   MaxMetJetDPhiIndex(int PFJID=0, float minJPt=20, float maxJEta=6., int met=1); // electrons and muons not considered for minDPhi
  Float_t MinMetJetDPhiL2L3 ();
  Float_t PseudoJetMetDPhi();
  Float_t GetPseudoJetMetDPhi(int hemi_index=1, int pj=1, int whichmet=1, float met=30);
  Float_t GetPseudoJetsdPhi(int hemi_seed=2, int hemi_association=3, int PFJID=0, double minJPt=20, double maxJEta=2.4);
  Float_t GetPseudoJetsdPhiOnline(int PFJID=0, double minJPt=20, double maxJEta=2.4);
  Float_t PseudoJetPtRatio(Bool_t inclMET, Bool_t vsHT);
  Float_t GetBJetDR(int algo, float value, float minJPt, float maxJEta, int PFJID);
  Float_t BJetMETdPhi(int algo, float value, float minJPt, float maxJEta, int PFJID);
  Float_t MHTMETDPhi();

  // MT2 & friends
  Float_t GetMT2            (Float_t testmass=0, bool massive=false,       Int_t PFJID=1, Float_t minJPt=20, Float_t maxJEta=2.4, 
		              Int_t hemi_seed=2,   Int_t hemi_association=3, Int_t met=1 );
  Float_t GetMT2MinDHT      (Float_t testmass=0, bool massive=false,       Int_t PFJID=1, Float_t minJPt=20, Float_t maxJEta=2.4, Int_t met=1);
  Bool_t  FillMT2Hemi       (Float_t testmass=0, bool massive=false,       Int_t PFJID=1, Float_t minJPt=20, Float_t maxJEta=2.4, 
		              Int_t hemi_seed=2,   Int_t hemi_association=3, Int_t met=1,   Int_t hemi_nr=-1);
  Bool_t  FillMT2HemiMinDHT (Float_t testmass=0, bool massive=false,       Int_t PFJID=1, Float_t minJPt=20, Float_t maxJEta=2.4, Int_t met=1, Int_t hemi_nr=-1);
  Float_t GetMCTcorr(TLorentzVector p1, TLorentzVector p2, TLorentzVector DTM, TVector2 pmiss);
  Float_t CalcMT2(float testmass, bool massive, 
		  TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET );
  Float_t GetMT(TLorentzVector lv1, float m1, TLorentzVector lv2, float m2);
  Float_t GetMT(TLorentzVector lv1, TLorentzVector lv2);
  Float_t SimpleMT2(bool pseudo=true, int heminr=1);
  Float_t GetSqrtS(float testmass=0, bool massive=true,int PFJID=0, float minJPt=20, float maxJEta=2.4,int met=1);
  Float_t GetMaxHemiMass(int hemi_index=0);

  // Leptons
  Float_t GenOSDiLeptonInvMass(unsigned int pid=11, unsigned int mother=23, float pt=10, float eta=2.4);
  Float_t GetDiLeptonInvMass(int same_sign=0, int same_flavour=1, int flavour=0, float pt=10, bool exclDiLept=false);
  Bool_t  IsDiLeptonMll(int same_sign=0, int same_flavour=1, int flavour=0, float pt=10, bool exclDiLept=false, float lower_mass=71, float upper_mass=111);
  Bool_t  IsGenOSDiLepton(unsigned int pid=11, unsigned int mother=23, float pt=10, float eta=2.4, float lower_mass=71, float upper_mass=111);
  TLorentzVector GetMETPlusLeptsLV(int OSDiLeptFromZ =1);
  Float_t GetMETPlusLepts(int OSDiLeptFromZ =1);
  Float_t GetMETPlusGenLepts(int met, int RemoveOSSFDiLepts=0, int require_cuts=1, unsigned int pid=11, 
		              unsigned int mother=23, float pt=10, float eta=2.4, float lower_mass=71, float upper_mass=111);
  Float_t GetDiLeptonPt(int same_sign=0, int same_flavour=1, int flavour=0, float pt=10, float lower_mass=71, float upper_mass=111);
  Float_t GetGenLeptPt(int which, int pid, int mother, float pt, float eta);
  Float_t GetGenLeptEta(int which, int pid, int mother, float pt, float eta);
  Int_t   GetGenLeptIndex(int which, int pid, int mother, float pt, float eta);
  Float_t GetGenLeptPt2(int which, int pid, int mother, int gmother, float pt, float eta);
  Float_t GetGenLeptEta2(int which, int pid, int mother, int gmother, float pt, float eta);
  Int_t   GetGenLeptIndex2(int which, int pid, int mother, int gmother, float pt, float eta);
  Bool_t  GenLeptFromW(int pid, float pt, float eta, bool includeTaus);
  Int_t   GenNumLeptFromW(int pid, float pt, float eta, bool includeTaus);
  Int_t   GenNumLeptFromMID(int mid, int pid, float pt, float eta, bool includeTaus);
  Float_t GetLeptPt(int index);
  Float_t ElClosestJet();
  Int_t   WDecayMode();
  Int_t   TopDecayMode();
  Bool_t  TopDecayModeResult(Int_t nlepts);
  Bool_t  SLTopAccept(float pt, float eta);
  Float_t SLTopEta(float pt);
  Float_t LeptJetDR(int pid, int index, bool bjet, int ID);
  Float_t GenDiLeptPt(Float_t l_ptmin, Float_t l_etamax, Float_t mll_min, Float_t mll_max, Bool_t charged);
  Float_t GenDiLeptRapidity(Float_t l_ptmin, Float_t l_etamax, Float_t mll_min, Float_t mll_max, Bool_t charged);
  TLorentzVector GenDiLeptLv(Float_t l_ptmin, Float_t l_etamax, Float_t mll_min, Float_t mll_max, Bool_t charged);
  Float_t GenZPt();
  TLorentzVector RecoOSDiLeptLv(Float_t l_ptmin, Float_t l_etamax, Float_t mll_min, Float_t mll_max);
  Float_t RecoOSDiLeptPt(Float_t l_ptmin, Float_t l_etamax, Float_t mll_min, Float_t mll_max);
  Float_t RecoOSDiLeptM(Float_t l_ptmin, Float_t l_etamax, Float_t mll_min, Float_t mll_max);
  Float_t RecoOSDiLeptRapidity(Float_t l_ptmin, Float_t l_etamax, Float_t mll_min, Float_t mll_max);
  Bool_t  ZllRecalculate();
  Bool_t  ZAcceptance(float minleptpt, float maxleptpt, float minlepteta, float maxlepteta);

  // Photons
  Int_t   GenPhotonGenJetDRJIndex(float minJPt, float maxJEta, int PFJID );
  Float_t GenPhotonGenJetDR(float minJPt, float maxJEta, int PFJID );
  Int_t   PhotonJetDRJIndex(int ph_index, float minJPt, float maxJEta, int PFJID );
  Float_t PhotonJetDR(int ph_index, float minJPt, float maxJEta, int PFJID );
  Int_t   PhotonEleDREIndex(int ph_index, float minEPt, float maxEEta);
  Float_t PhotonEleDR(int ph_index, float minEPt, float maxEEta);
  Float_t GenPhotonAndLeadingJetPt();
  Float_t PhotonToJetPtRatio();
  Float_t MinGenBosonJetsDR();

  // PrintOut 
  Bool_t   PrintOut(Bool_t logfile);

  //Bosons
  Float_t GetGenVPt(int pid);

  Int_t   NJets;
  Int_t   NGenJets;
  Int_t   NJetsIDLoose;
  Int_t   NJetsIDLoose40;
  Int_t   NJetsIDLoose50;
  Int_t   NBJetsCSVM;
  Int_t   NBJetsCSVT;
  Int_t   NBJets40CSVM;
  Int_t   NBJets40CSVT;
  Int_t   NEles;
  Int_t   NMuons;
  Int_t   NMuonsCommonIso;//can be probably be deleted soon, again? (common object definition hadronic analyses)
  Int_t   NPhotons;
  Int_t   NPhotonsIDLoose;
  Int_t   NPhotonsIDLoose25;//can be probably be deleted soon, again? (common object definition hadronic analyses)
  Int_t   NTaus;
  Int_t   NTausIDLoose;
  Int_t   NTausIDLoose3Hits;
  Int_t   NTausIDLooseMVA;
  Int_t   NTausIDLoose2;//can be probably be deleted soon, again? (common object definition hadronic analyses)
  Int_t   NGenLepts;
  Int_t   NPdfs;
  Int_t GenProcessID;
  Double_t GenWeight;

  MT2Susy        Susy;
  MT2Misc        misc;
  MT2PileUp      pileUp;
  MT2Trigger     trigger;
  MT2SFWeight    SFWeight;
  MT2Jet         jet[m_jetSize];
  MT2GenJet      genjet[m_genjetSize];
  MT2Hemi        hemi[m_hemiSize];
  MT2Tau         tau[m_tauSize];
  MT2Elec        ele[m_eleSize];
  MT2Muon        muo[m_muoSize];
  MT2Photon      photon[m_phoSize];
  MT2GenLept     genlept[m_genleptSize];
  TLorentzVector pfmet[2];
  TLorentzVector genmet[2];
  TLorentzVector MHT[2];
  TLorentzVector GenPhoton[2];
  TLorentzVector GenZ[2];
  TLorentzVector type1pfmet[2]; // type1 corrected PFMET
  TLorentzVector rawpfmet[2]; // this is raw, uncorrected and not scaled pf-met. 
                              // identical to pfmet unless met is JEScales or modified in data-driven estimates
  double pdfW[100];
  
  ClassDef(MT2tree, 29)
};

#endif
