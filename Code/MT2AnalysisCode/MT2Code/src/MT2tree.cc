#include "MT2tree.hh"

#include "helper/Davismt2.h"
#include "helper/TMctLib.h"

#include <vector>
#include "helper/Hemisphere.hh"
#include "helper/Utilities.hh"
#include "TLorentzVector.h"
#include "TMath.h"

#include <cmath>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

// MT2Misc -----------------------------------
MT2Misc::MT2Misc(){
  Reset();
}

MT2Misc::~MT2Misc(){
}

void MT2Misc::Reset() {
  HBHENoiseFlag           =  0;
  HBHENoiseFlagIso        =  0;
  RecovRecHitFilterFlag   =  0;
  BadEcalTP               =  0;
  BadEcalBE               =  0;
  CSCTightHaloID          =  0;
  CrazyHCAL               =  0;
  NegativeJEC             =  0;
  isData                  =  0;
  Run                     = -1;	  
  Event		  	  = -1;	  
  LumiSection		  = -1;	  
  LeptConfig		  = -1;	  
  PassJetID               = -1;
  Jet0Pass                = -1;
  Jet1Pass                = -1;
  MT2                     = -99999.99;
  MCT                     = -99999.99;
  MET                     = -99999.99;
  METPhi                  = -99999.99;
  LeadingJPt              = -99999.99;
  SecondJPt               = -99999.99;
  Vectorsumpt		  = -99999.99;
  HT			  = -99999.99;
  QCDPartonicHT		  = -99999.99;
  caloHT40       	  = -99999.99;
  caloHT50       	  = -99999.99;
  caloHT50_ID             = -99999.99;
  caloMHT30       	  = -99999.99;
  caloMHT40       	  = -99999.99;
  MinMetJetDPhi           = -99999.99;
  MinMetJetDPhiIndex      = -1;
  MinMetBJetDPhi          = -99999.99;
  TrackingFailure         = -99999.99;
  TrackingFailurePVtx     = -99999.99;

}

// ------------------------------------------------------------
// MT2PileUp
MT2PileUp::MT2PileUp(){
	Reset();
}

MT2PileUp::~MT2PileUp(){
}

void MT2PileUp::Reset(){
	PUnumInt       = -999;
	PUnumIntEarly  = -999;
	PUnumIntLate   = -999;
	isS3           = -1;
	PtHat          = -999.99;
	Weight         = -999.99;
	Rho            = -999.99;
  	NVertices      = -1;
}

// ------------------------------------------------------
// MT2trigger
MT2Trigger::MT2Trigger(){
	Reset();
}
MT2Trigger::~MT2Trigger(){
}

void MT2Trigger::Reset(){
	
	// HT with DPhi
	HLT_HT500_JetPt60_DPhi2p94_v1 = false;
	HLT_HT550_JetPt60_DPhi2p94_v1 = false;
	
	// HT
	HLT_HT150_v2            = false;
	HLT_HT150_v3            = false;
	HLT_HT160_v2            = false;
	HLT_HT200_v2            = false;
	HLT_HT200_v3            = false;
	HLT_HT240_v2            = false;
	HLT_HT250_v2            = false;
	HLT_HT250_v3            = false;
	HLT_HT260_v2            = false;
	HLT_HT300_v2            = false;
	HLT_HT300_v3            = false;
	HLT_HT300_v4            = false;
	HLT_HT300_v5            = false;
	HLT_HT350_v2            = false;
	HLT_HT350_v3            = false;
	HLT_HT350_v4            = false;
	HLT_HT360_v2            = false;
	HLT_HT400_v10           = false;
	HLT_HT400_v2            = false;
	HLT_HT400_v3            = false;
	HLT_HT400_v4            = false;
	HLT_HT400_v5            = false;
	HLT_HT400_v6            = false;
	HLT_HT400_v7            = false;
	HLT_HT400_v8            = false;
	HLT_HT400_v9            = false;
	HLT_HT440_v2            = false;
	HLT_HT450_v10           = false;
	HLT_HT450_v2            = false;
	HLT_HT450_v3            = false;
	HLT_HT450_v4            = false;
	HLT_HT450_v5            = false;
	HLT_HT450_v6            = false;
	HLT_HT450_v7            = false;
	HLT_HT450_v8            = false;
	HLT_HT450_v9            = false;
	HLT_HT500_v10           = false;
	HLT_HT500_v2            = false;
	HLT_HT500_v3            = false;
	HLT_HT500_v4            = false;
	HLT_HT500_v5            = false;
	HLT_HT500_v6            = false;
	HLT_HT500_v7            = false;
	HLT_HT500_v8            = false;
	HLT_HT500_v9            = false;
	HLT_HT550_v10           = false;
	HLT_HT550_v2            = false;
	HLT_HT550_v3            = false;
	HLT_HT550_v4            = false;
	HLT_HT550_v5            = false;
	HLT_HT550_v6            = false;
	HLT_HT550_v7            = false;
	HLT_HT550_v8            = false;
	HLT_HT550_v9            = false;
	HLT_HT600_v1            = false;
	HLT_HT600_v2            = false;
	HLT_HT600_v3            = false;
	HLT_HT650_v1            = false;
	HLT_HT650_v2            = false;
	HLT_HT650_v3            = false;
	// HT_MHT
	HLT_HT250_MHT60_v2      = false;
	HLT_HT250_MHT60_v3      = false;
	HLT_HT250_MHT60_v4      = false;
	HLT_HT250_MHT60_v5      = false;
	HLT_HT250_MHT60_v6      = false;
	HLT_HT250_MHT70_v1      = false;
	HLT_HT250_MHT70_v2      = false;
	HLT_HT250_MHT70_v3      = false;
	HLT_HT250_MHT70_v4      = false;
	HLT_HT250_MHT90_v1      = false;
	HLT_HT250_MHT90_v2      = false;
	HLT_HT260_MHT60_v2      = false;
	HLT_HT300_MHT75_v4      = false;
	HLT_HT300_MHT75_v5      = false;
	HLT_HT300_MHT75_v7      = false;
	HLT_HT300_MHT75_v8      = false;
	HLT_HT300_MHT80_v1      = false;
	HLT_HT300_MHT80_v2      = false;
	HLT_HT300_MHT90_v1      = false;
	HLT_HT300_MHT90_v2      = false;
	HLT_HT350_MHT70_v1      = false;
	HLT_HT350_MHT70_v2      = false;
	HLT_HT350_MHT80_v1      = false;
	HLT_HT350_MHT80_v2      = false;
	// Muons
	HLT_DoubleMu3_HT160_v2  = false;
	HLT_Mu8_Jet40_v2        = false;
	HLT_DoubleMu3_v3        = false;
}


// MT2Znunu ------------------------------------
MT2Znunu::MT2Znunu(){
	Reset();
}

MT2Znunu::~MT2Znunu(){
}

void MT2Znunu::Reset(){
	NJetsIDLoose_matched      = -999;
	PassJetID_matched         = -999;
	Jet1Pass_matched          = -999;
	Jet0Pass_matched          = -999;
	LeadingJPt_matched        = -99999.99;
	SecondJPt_matched         = -99999.99;
	HTmatched                 = -99999.99;
	caloMHT30_matched         = -99999.99;
	caloMHT30ID_matched       = -99999.99;
	caloMHT30_matchedReco     = -99999.99;
	caloMHT30ID_matchedReco   = -99999.99;
	caloHT50_matched          = -99999.99;
	caloHT50ID_matched        = -99999.99;
	caloHT50_matchedReco      = -99999.99;
	caloHT50ID_matchedReco    = -99999.99;
	GenZmumu_mll              = -99999.99;
	GenZmumu_mll_acc          = -99999.99;
	GenZee_mll                = -99999.99;
	GenZee_mll_acc            = -99999.99;
	GenZnunu_e_mll            = -99999.99;
	GenZnunu_e_mll_acc        = -99999.99;
	GenZnunu_mu_mll           = -99999.99;
	GenZnunu_mu_mll_acc       = -99999.99;
	GenZnunu_tau_mll          = -99999.99;
	GenZnunu_tau_mll_acc      = -99999.99;
	RecoOSee_mll              = -99999.99;
	RecoOSmumu_mll            = -99999.99;
	METplusLeptsPt            = -99999.99;
	METplusLeptsPtReco        = -99999.99;
	MinMetplusLeptJetDPhi     = -99999.99;
	MinMetplusLeptJetDPhiReco = -99999.99;
	Vectorsumpt_matched       = -99999.99;
}


// MT2Jet -----------------------------------
MT2Jet::MT2Jet(){
  Reset();
}

MT2Jet::~MT2Jet(){
}

void MT2Jet::Reset() {
  lv.SetPxPyPzE(0, 0, 0, 0);
  bTagProbTCHE  = -99999.99;
  bTagProbTCHP  = -99999.99;
  bTagProbSSVHE = -99999.99;
  bTagProbSSVHP = -99999.99;

  isPATPFIDLoose= 0; // PAT PFJetIDLoose regardless of eta and pt
  isPFIDLoose   = 0;
  isPFIDMedium  = 0;
  isPFIDTight   = 0;
  ChHadFrac     = -99999.99; 
  NeuHadFrac    = -99999.99; 
  ChEmFrac      = -99999.99;
  NeuEmFrac     = -99999.99; 
  ChMult        = -1; 
  NeuMult       = -1; 
  NConstituents = -1;

  Scale         = -99999.99; // correction factor
  L1FastJetScale= -99999.99; // correction factor from raw to L1FastJetcorrected
  Area          = -99999.99;

  Flavour       = -9999;
  
  isTauMatch    = 0;  // tell you if the jet is matched to a tau.
  TauDR         = -99999.99;
  TauDPt        = -99999.99;
  NTauMatch     = 0;
}

void MT2Jet::SetLV(const TLorentzVector v) {
  lv = v;
}

Bool_t MT2Jet::IsGoodPFJet(float minJPt, float maxJEta, int PFJID) {

  float pt = lv.Pt();
  float eta = lv.Eta();
  if ( pt < minJPt || fabs(eta) > maxJEta )     return false;
  if ( Scale <0)                                return false; // jet has negative Scale from JE correction
  
  switch (PFJID) {
  case 3:               // TIGHT
    if ( !(NeuEmFrac     < 0.90) ) return false;
    if ( !(NeuHadFrac    < 0.90) ) return false;
    // break;   // No break: medium contains tight -> check common cuts
  case 2:               // MEDIUM
    if ( !(NeuEmFrac     < 0.95) ) return false;
    if ( !(NeuHadFrac    < 0.95) ) return false;
    // break;   // No break: loose contains medium -> check common cuts
  case 1:               // LOOSE
    if ( ! isPATPFIDLoose     )    return false; // loose PF-jet ID from PAT
    break;
  default:
    // None of the above. Do we want any default cut?
    break;
  }

  return true;
}

Bool_t MT2Jet::IsBJet(Int_t algo) {
  if     (algo==3 && bTagProbSSVHP >2.0 ) return true;
  else if(algo==2 && bTagProbSSVHE >1.74) return true;
  else                                    return false;
}

// MT2GenJet --------------------------------------
MT2GenJet::MT2GenJet(){
  Reset();
}

MT2GenJet::~MT2GenJet(){
}

void MT2GenJet::Reset(){
  DeltaR             = -9999.99;
  JetMatchIndex      = -1;
  
  lv.SetPxPyPzE(0, 0, 0, 0);
}

// MT2Hemisphere -----------------------------------
MT2Hemi::MT2Hemi(){
  Reset();
}

MT2Hemi::~MT2Hemi(){
}

void MT2Hemi::Reset(){
  seed_method    = -1;
  assoc_method   = -1;
  NHemiIter      = -1;
  MT2            = -99999.99;
  MCT            = -99999.99;
  AlphaT         = -99999.99;
  minDHT         = -99999.99;
  maxDR          = -99999.99;
  dPhi           = -99999.99;
  for(int i=0; i<m_jetSize; ++i){
  	jindices1[i]  =-1;
  	jindices2[i]  =-1;
  }
  for(int i=0; i<m_eleSize; ++i){
  	eleindices1[i]=-1;
  	eleindices2[i]=-1;
  }
  for(int i=0; i<m_muoSize; ++i){
  	muoindices1[i]=-1;
  	muoindices2[i]=-1;
  }
  
  lv1. SetPxPyPzE(0, 0, 0, 0);
  lv2. SetPxPyPzE(0, 0, 0, 0);
  UTM. SetPxPyPzE(0, 0, 0, 0); 
}


// MT2GenLept -----------------------------------
MT2GenLept::MT2GenLept(){
  Reset();
}

MT2GenLept::~MT2GenLept(){
}

void MT2GenLept::Reset(){
  lv.  SetPxPyPzE(0, 0, 0, 0);

  ID          = -999;
  MID         = -999;
  MStatus     = -999;
  GMID        = -999;
  GMStatus    = -999;
  MT          = -9999.99;
  CAJ_n90     = -9999.99;
  CAJ_n90Hits = -9999.99;
}

// MT2Muon -----------------------------------
MT2Muon::MT2Muon(){
  Reset();
}

MT2Muon::~MT2Muon(){
}

void MT2Muon::Reset() {
  lv.SetPxPyPzE(0, 0, 0, 0);
  MT            = -9999.99;
  Iso           = -9999.99;
  Charge        = -999;
  NMatches      = -999;
  PtErr         = -999.99;
}

void MT2Muon::SetLV(const TLorentzVector v) {
  lv = v;
}

// MT2Elec -----------------------------------
MT2Elec::MT2Elec(){
  Reset();
}

MT2Elec::~MT2Elec(){
}

void MT2Elec::Reset() {
  lv.SetPxPyPzE(0, 0, 0, 0);
  MT            = -9999.99;
  Iso           = -9999.99;
  Charge        = -999;
  ID95          = -999;
  ID90          = -999;
  CAJ_n90       = -999;
  CAJ_n90Hits   = -999;
}

void MT2Elec::SetLV(const TLorentzVector v) {
  lv = v;
}

// MT2tree ----------------------------------
MT2tree::MT2tree(){
  Reset();
}

MT2tree::~MT2tree(){
}

void MT2tree::Reset() {
  NJets            = 0;
  NTaus            = 0;
  NJetsIDLoose     = 0;
  NJetsIDLoose40   = 0;
  NJetsIDLoose50   = 0;
  NEles            = 0;
  NMuons           = 0;
  NGenLepts        = 0;
  NGenJets         = 0;

  misc.Reset();
  Znunu.Reset();
  pileUp.Reset();
  trigger.Reset();

  for (int i = 0; i < m_jetSize; ++i) {
    jet[i].Reset();
  }
  for (int i = 0; i< m_genjetSize; ++i){
    genjet[i].Reset();
  }
  for (int i = 0; i < m_eleSize; ++i) {
    ele[i].Reset();
  }
  for (int i = 0; i < m_muoSize; ++i) {
    muo[i].Reset();
  }
  for (int i = 0; i < m_genleptSize; ++i) {
    genlept[i].Reset();
  }
  for (int i = 0; i < m_hemiSize; ++i) {
    hemi[i].Reset();
  }
  pfmet     [0].SetPxPyPzE(0., 0., 0., 0.);
  genmet    [0].SetPxPyPzE(0., 0., 0., 0.);
  MHT       [0].SetPxPyPzE(0., 0., 0., 0.);
}

void MT2tree::SetNJets(int n) {
  NJets = n;
}

void MT2tree::SetNGenJets(int n) {
  NGenJets = n;
}

void MT2tree::SetNJetsIDLoose(int n) {
  NJetsIDLoose = n;
}

void MT2tree::SetNBJets(int n) {
  NBJets = n;
}

void MT2tree::SetNEles(int n) {
  NEles = n;
}

void MT2tree::SetNMuons(int n) {
  NMuons = n;
}

void MT2tree::SetNTaus(int n) {
  NTaus = n;
}

// --------------------------------------------------------
// NJets and friends
Int_t MT2tree::GetNjets(float minJPt, float maxJEta, int PFJID){
  int njets=0;
  for(int i=0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
	njets++;
  }
  return njets;
}

Int_t MT2tree::GetNBtags (int algo, float value, float minJPt, float maxJEta, int PFJID){  // algo - 0:TCHE, 1:TCHP, 2:SSVHE, 3:SSVHP
  int nbjets=0;
  for(int i=0; i<NJets; ++i){
    if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false ) continue; 
    switch(algo){
    case 0: 
      if( jet[i].bTagProbTCHE < value ) continue;
      break;
    case 1: 
      if( jet[i].bTagProbTCHP < value ) continue;
      break;
    case 2: 
      if( jet[i].bTagProbSSVHE < value ) continue;
      break;
    case 3: 
      if( jet[i].bTagProbSSVHP < value ) continue;
      break;
    default:
      continue;
      break;
    }
    nbjets++; 
  }
  return nbjets;
}

Float_t MT2tree::JetPt(int ijet, int PFJID, float minJPt, float maxJEta) {
  int index = GetJetIndex(ijet, PFJID,minJPt,maxJEta);
  if ( index < 0 )       return -999.;
  return jet[index].lv.Pt();
}

Int_t MT2tree::GetJetIndex(int ijet, int PFJID, float minJPt, float maxJEta) {
  std::vector<int> indices;
  for(int i = 0; i<NJets; ++i){
    if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
    indices.push_back(i);
  }
  if (ijet>=indices.size())  return -9;
  return indices[ijet];
}

Int_t MT2tree::GetJetIndexByEta(int ijet, int PFJID, float minJPt, float maxJEta) {
  std::vector<int> indices, indx;
  std::vector<float> etas;
  for(int i = 0; i<NJets; ++i){
    if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
    indices.push_back(i);
    etas.push_back(fabs(jet[i].lv.Eta()));
  }
  if (ijet>=indices.size())  return -9;
  indx = Util::VSort(indices,etas,true);
  return indx[ijet];
}

// ----------------------------------------------------------------------------
// dPhi and friends
Float_t MT2tree::GetMinR12R21(int PFJID, float minJPt, float maxJEta, int met){
	TLorentzVector MET(0., 0., 0., 0.);
	if(met==1)      MET=pfmet[0];
        else if(met==2) MET=MHT[0];
	else            return -900;

	vector<int> indices;
	for(int i=0; i<NJets; ++i){
		if( jet[i].IsGoodPFJet(minJPt, maxJEta, PFJID)==false) continue;
		indices.push_back(i);
	}
	if(indices.size()<2) return -999;
	float d_phi1 = fabs(Util::DeltaPhi(jet[indices[0]].lv.Phi(), MET.Phi()));	
	float d_phi2 = fabs(Util::DeltaPhi(jet[indices[1]].lv.Phi(), MET.Phi()));	
	float R12    = sqrt( d_phi1*d_phi1 + (TMath::Pi()-d_phi2)*(TMath::Pi()-d_phi2) );
	float R21    = sqrt( d_phi2*d_phi2 + (TMath::Pi()-d_phi1)*(TMath::Pi()-d_phi1) );

	if(R12<R21) return R12;
	else        return R21;
}

Bool_t MT2tree::PassJetID(float minJPt, float maxJEta, int PFJID) {
	int njets=0;
	for(int i=0; i<NJets; ++i){
		if(jet[i].lv.Pt() >= minJPt && fabs(jet[i].lv.Eta()) <= maxJEta && 
		   jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false)   return false;
	}
	return true;
}

Float_t MT2tree::JetsInvMass(int j1, int j2){
	if(NJets < 2) return -999.99;
	TLorentzVector sum = jet[j1].lv + jet[j2].lv;
	return sum.M();
}

Float_t MT2tree::JetsDPhi(int j1, int j2, int PFJID){
  std::vector<int> indices;
  for(int i = 0; i<NJets; ++i){
  	if(jet[i].isPFIDLoose ==false && PFJID==1  )continue;
  	if(jet[i].isPFIDMedium==false && PFJID==2  )continue;
  	if(jet[i].isPFIDTight ==false && PFJID==3  )continue;
	indices.push_back(i);
  }
  if (j1>=indices.size() || j2>=indices.size())  return -999;
  return TMath::Abs(jet[indices[j1]].lv.DeltaPhi(jet[indices[j2]].lv));
}

Float_t MT2tree::MetJetDPhi(int ijet, int PFJID, int met) {
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHT[0];
  else            return -999;

  std::vector<int> indices;
  for(int i = 0; i<NJets; ++i){
  	if(jet[i].isPFIDLoose ==false && PFJID==1  )continue;
  	if(jet[i].isPFIDMedium==false && PFJID==2  )continue;
  	if(jet[i].isPFIDTight ==false && PFJID==3  )continue;
	indices.push_back(i);
  }
  if (ijet>=indices.size())  return -999;
  return TMath::Abs(jet[indices[ijet]].lv.DeltaPhi(MET));
}

Float_t MT2tree::MinMetJetDPhi(int PFJID, float minJPt, float maxJEta, int met) {
// Attention: electrons and muons are not considered for minDPhi
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = GetMHTlv(PFJID, minJPt, maxJEta);
  else if(met==3) {
	  float mass = GetDiLeptonInvMass(0,1,0,5.0,1);
	  if(mass > 71 && mass <111) MET = GetMETPlusLeptsLV(1) + pfmet[0] ;
	  else return -888.;
  }
  else         return -999;


  int index = MinMetJetDPhiIndex(PFJID, minJPt, maxJEta, met);
  if(index >=0) {
	  return TMath::Abs(jet[index].lv.DeltaPhi(MET));
  } else  return -999.99;
}

Int_t MT2tree::MinMetJetDPhiIndex(int PFJID, float minJPt, float maxJEta, int met) {
// Attention: electrons and muons are not considered for minDPhi
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = GetMHTlv(PFJID, minJPt, maxJEta);
  else if(met==3) {
	  float mass = GetDiLeptonInvMass(0,1,0,5.0,1);
	  if(mass > 71 && mass <111) MET = GetMETPlusLeptsLV(1) + pfmet[0] ;
	  else return -888.;
  }
  else            return -999;

  std::vector<int> indices;
  for(int i = 0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
	indices.push_back(i);
  }
  if(indices.size()<1)  return -999;
  Float_t minDPhi=10;
  Int_t    imin;
  for (int i=0; i<indices.size(); i++){
    Float_t dphi = TMath::Abs(jet[indices[i]].lv.DeltaPhi(MET));
    if(dphi<minDPhi)  {
      minDPhi=dphi;
      imin = indices[i];
    }
  }
  return imin;
}

Float_t MT2tree::MinMetJetDPhiL2L3(){
	vector<TLorentzVector> jets;
	for (int i=0; i<NJets; ++i){
		TLorentzVector j(0,0,0,0);
		if(jet[i].lv.Pt()/jet[i].L1FastJetScale < 20) continue;
		if(fabs(jet[i].lv.Eta())                > 5 ) continue;
		j.SetPtEtaPhiM(jet[i].lv.Pt()/jet[i].L1FastJetScale, jet[i].lv.Eta(), jet[i].lv.Phi(), jet[i].lv.M());
		jets.push_back(j);
	}
	Float_t minDPhi=10;
	for(int i=0; i<jets.size(); ++i){
		Float_t dphi = TMath::Abs(jets[i].DeltaPhi(pfmet[0]));
		if(dphi<minDPhi) minDPhi=dphi;
	}
	return minDPhi;
}

Float_t MT2tree::MaxMetJetDPhi(int PFJID, float minJPt, float maxJEta, int met) {
// Attention: electrons and muons are not considered for minDPhi
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHT[0];
  else            return -999;

  int index = MaxMetJetDPhiIndex(PFJID, minJPt, maxJEta, met);
  return TMath::Abs(jet[index].lv.DeltaPhi(MET));
}

Int_t MT2tree::MaxMetJetDPhiIndex(int PFJID, float minJPt, float maxJEta, int met) {
// Attention: electrons and muons are not considered for minDPhi
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHT[0];
  else            return -999;

  std::vector<int> indices;
  for(int i = 0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
	indices.push_back(i);
  }
  if(indices.size()<1)  return -999;
  Float_t maxDPhi=0.;
  Int_t    imax;
  for (int i=0; i<indices.size(); i++){
    Float_t dphi = TMath::Abs(jet[indices[i]].lv.DeltaPhi(MET));
    if(dphi>maxDPhi)  {
      maxDPhi=dphi;
      imax = indices[i];
    }
  }
  return imax;
}

Int_t MT2tree::BiasedDPhiIndex(int PFJID, float minJPt, float maxJEta){
  std::vector<int> indices;
  for(int i = 0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
	indices.push_back(i);
  }
  if(indices.size()<2)  return -999;
  Float_t minDPhi=10.;
  Int_t    imin;
  for (int i=0; i<indices.size(); i++){
    TLorentzVector recoil(0,0,0,0);
    for(int j=0; j<indices.size(); j++){
    	if(j==i) continue;
	recoil -= jet[j].lv;
    }
    Float_t dphi = TMath::Abs(jet[i].lv.DeltaPhi(recoil));
    if(dphi<minDPhi)  {
      minDPhi=dphi;
      imin = indices[i];
    }
  }
  return imin;
}

Float_t MT2tree::BiasedDPhi(int PFJID, float minJPt, float maxJEta) {
  int index = BiasedDPhiIndex(PFJID, minJPt, maxJEta);
  TLorentzVector recoil(0,0,0,0);
  if(index <0) return -999.99;
  for(int i = 0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
	if(i== index) continue;
	recoil -= jet[i].lv;
  }
  return TMath::Abs(jet[index].lv.DeltaPhi(recoil));
}


Bool_t MT2tree::PassMinMetJetDPhi03(){
	if( NJetsIDLoose < 3)                              return true;
	if( NJetsIDLoose >=3  && misc.MinMetJetDPhi > 0.3) return true;
	return false;
}

// Maximum hemi mass ---------------------------------------------------------------
Float_t MT2tree::GetMaxHemiMass(int hemi_index){
	return TMath::Max(hemi[hemi_index].lv1.M(),hemi[hemi_index].lv2.M());
}

// Pseudojet-met-DPhi ------------------------------------------------------------
Float_t MT2tree::GetPseudoJetMetDPhi(int hemi_index, int pj, int whichmet, float met){
	if(whichmet==1 && pfmet[0].Pt() < met)                                    return -777;
	if(whichmet==2 && (hemi[hemi_index].lv1+hemi[hemi_index].lv2).Pt() < met) return -777;
	if(pj ==1){
		if(whichmet==1)        return TMath::Abs(hemi[hemi_index].lv1.DeltaPhi(pfmet[0]));
		else if (whichmet ==2) return TMath::Abs(hemi[hemi_index].lv1.DeltaPhi(-hemi[hemi_index].lv1 -hemi[hemi_index].lv2));
		else                   return -111;
	}else if (pj ==2){
		if(whichmet==1)        return TMath::Abs(hemi[hemi_index].lv2.DeltaPhi(pfmet[0]));
		else if (whichmet ==2) return TMath::Abs(hemi[hemi_index].lv2.DeltaPhi(-hemi[hemi_index].lv1 -hemi[hemi_index].lv2));
		else                   return -111;       
	}else {return -999.99;}
}

// Pseudojet-met-DPhi ------------------------------------------------------------
Float_t MT2tree::PseudoJetMetDPhi(){
	if(hemi[0].MT2 <0) return -999;
	float dPhi1 = hemi[0].lv1.DeltaPhi(pfmet[0]);
	float dPhi2 = hemi[0].lv2.DeltaPhi(pfmet[0]);
	return (fabs(dPhi1) < fabs(dPhi2)) ? fabs(dPhi1) : fabs(dPhi2);
}

// PseudoJetPtRatio -----------------------------------------------
Float_t MT2tree::PseudoJetPtRatio(Bool_t inclMET, Bool_t vsHT){
	if(! inclMET && !vsHT){
		if(hemi[0].lv1.Pt() > hemi[0].lv2.Pt()) return hemi[0].lv2.Pt()/hemi[0].lv1.Pt();
		else return hemi[0].lv1.Pt()/hemi[0].lv2.Pt();
	}else if(inclMET && !vsHT){
		if(hemi[0].lv1.Pt() > hemi[0].lv2.Pt()) return (hemi[0].lv2.Pt()+misc.MET)/hemi[0].lv1.Pt();
		else return (hemi[0].lv1.Pt()+misc.MET)/hemi[0].lv2.Pt();
	}else if(inclMET && vsHT){
		if(hemi[0].lv1.Pt() > hemi[0].lv2.Pt()) return (hemi[0].lv2.Pt()+misc.MET)/misc.HT;
		else return (hemi[0].lv1.Pt()+misc.MET)/misc.HT;
	}else if(!inclMET && vsHT){
		if(hemi[0].lv1.Pt() > hemi[0].lv2.Pt()) return (hemi[0].lv2.Pt())/misc.HT;
		else return (hemi[0].lv1.Pt())/misc.HT;
	}
}


// BJetDR --------------------------------------------------------
Float_t MT2tree::GetBJetDR(int algo, float value, float minJPt, float maxJEta, int PFJID){

	// get highest pt b-jet index
	float jpt=0; int jindex=-1;
	for(int i=0; i<NJets; ++i){
		if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue; 
		switch(algo){
			case 0: 
			if( jet[i].bTagProbTCHE < value ) continue;
			break;
			case 1: 
			if( jet[i].bTagProbTCHP < value ) continue;
			break;
			case 2: 
			if( jet[i].bTagProbSSVHE < value ) continue;
			break;
			case 3: 
			if( jet[i].bTagProbSSVHP < value ) continue;
			break;
			default:
			continue;
			break;
		}
		if(jet[i].lv.Pt()>jpt) jindex=i;
	}
	if(jindex==-1) return -999.99;

	// get dR to closest Jet
	float minDR=100;
	for(int i=0; i<NJets; ++i){
		if(jet[i].IsGoodPFJet(20,5,0)==false) continue;
		if(i==jindex)                         continue;
		float dR=TMath::Abs(jet[jindex].lv.DeltaR(jet[i].lv));
		if(dR < minDR) minDR=dR;
	}
	if(minDR==100) return -999.99;
	else           return minDR;
}

Float_t MT2tree::BJetMETdPhi(int algo, float value, float minJPt, float maxJEta, int PFJID){
	Float_t minDPhi =10;
	for(int i=0; i<NJets; ++i){
		if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue; 
		switch(algo){
			case 0: 
			if( jet[i].bTagProbTCHE < value ) continue;
			break;
			case 1: 
			if( jet[i].bTagProbTCHP < value ) continue;
			break;
			case 2: 
			if( jet[i].bTagProbSSVHE < value ) continue;
			break;
			case 3: 
			if( jet[i].bTagProbSSVHP < value ) continue;
			break;
			default:
			continue;
			break;
		}
		Float_t dPhi =TMath::Abs(jet[i].lv.DeltaPhi(pfmet[0]));
		if(dPhi < minDPhi) minDPhi = dPhi;
	}
	if(minDPhi==10) return -999.99;
	else            return minDPhi;


}
// ----------------------------------------------------------------
// HT and friends
Float_t MT2tree::GetHT(int PFJID, float minJPt, float maxJEta){
  Float_t ht=0;
  for(int i = 0; i<NJets; ++i){
    if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
    ht+=jet[i].lv.Pt();
  }
  return ht;
}

TLorentzVector MT2tree::GetMHTlv(int PFJID, float minJPt, float maxJEta, bool inclLepts){
  TLorentzVector mht(0,0,0,0);
  TLorentzVector j(0,0,0,0);
  for(int i = 0; i<NJets; ++i){
    if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
    j.SetPtEtaPhiM(jet[i].lv.Pt(),0,jet[i].lv.Phi(),0);
    mht-=j;
  }
  if(!inclLepts) return mht;
  // add leptons
  for(int i=0; i<NEles; ++i){
    j.SetPtEtaPhiM(ele[i].lv.Pt(),0, ele[i].lv.Phi(), 0);
    mht-=j;
  }
  for(int i=0; i<NMuons; ++i){
    j.SetPtEtaPhiM(muo[i].lv.Pt(),0, muo[i].lv.Phi(), 0);
    mht-=j;
  }
  return mht;
}

Float_t MT2tree::GetMHT(int PFJID, float minJPt, float maxJEta, bool inclLepts){
  TLorentzVector mht = GetMHTlv(PFJID,minJPt,maxJEta, inclLepts);
  return mht.Pt();
}

Float_t MT2tree::GetMHTPhi(int PFJID, float minJPt, float maxJEta, bool inclLepts){
  TLorentzVector mht = GetMHTlv(PFJID,minJPt,maxJEta, inclLepts);
  return mht.Phi();
}

Float_t MT2tree::GetMHTminusMET(int PFJID, float minJPt, float maxJEta, bool inclLepts){
  TLorentzVector mht = GetMHTlv(PFJID,minJPt,maxJEta, inclLepts);
  return (mht-pfmet[0]).Pt();
}

// -------------------------------------------------------------------------------------------
// MT2 and friends

// calculate MT2 with hemispheres  ---------------------------------------------------------------------------------------------------------
Float_t MT2tree::GetMT2(Float_t testmass, bool massive, Int_t PFJID, Float_t minJPt, Float_t maxJEta, Int_t hemi_seed, Int_t hemi_association, Int_t met) {

  if(misc.CrazyHCAL ) return false; // protection against crazy HCAL events
  
  // fill Pseudojets with selected objects
  vector<float> px, py, pz, E;
  for(int i=0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt, maxJEta, PFJID) ==false) continue;
  	px.push_back(jet[i].lv.Px());
	py.push_back(jet[i].lv.Py());
	pz.push_back(jet[i].lv.Pz());
	E .push_back(jet[i].lv.E ());
  }
  for(int i=0; i<NEles; ++i){
	if(met==1113) continue; 
  	px.push_back(ele[i].lv.Px());
	py.push_back(ele[i].lv.Py());
	pz.push_back(ele[i].lv.Pz());
	E .push_back(ele[i].lv.E ());
  }
  for(int i=0; i<NMuons; ++i){
	if(met==1113) continue; 
  	px.push_back(muo[i].lv.Px());
	py.push_back(muo[i].lv.Py());
	pz.push_back(muo[i].lv.Pz());
	E .push_back(muo[i].lv.E ());
  }
  if (px.size()<2) return -999.99;
  

  // get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)
  Hemisphere* hemisp = new Hemisphere(px, py, pz, E, hemi_seed, hemi_association);
  vector<int> grouping = hemisp->getGrouping();

  TLorentzVector pseudojet1(0.,0.,0.,0.);
  TLorentzVector pseudojet2(0.,0.,0.,0.);
	
  for(int i=0; i<px.size(); ++i){
	if(grouping[i]==1){
		pseudojet1.SetPx(pseudojet1.Px() + px[i]);
		pseudojet1.SetPy(pseudojet1.Py() + py[i]);
		pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
		pseudojet1.SetE( pseudojet1.E()  + E[i]);	
	}else if(grouping[i] == 2){
		pseudojet2.SetPx(pseudojet2.Px() + px[i]);
		pseudojet2.SetPy(pseudojet2.Py() + py[i]);
		pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
		pseudojet2.SetE( pseudojet2.E()  + E[i]);
	}
  }
  delete hemisp;

  // define MET
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1){
	MET = pfmet[0];
  } else if(met==1113) { // add leptons to MET
	MET = pfmet[0]; 
	for(int i=0; i<NEles; ++i){
		MET = MET + ele[i].lv;
	}
	for(int i=0; i<NMuons; ++i){
		MET = MET + muo[i].lv;
	}
  } else if(met==4) {MET = -pseudojet1 - pseudojet2;}
  else return -999.99;

  Float_t MT2 = CalcMT2(testmass, massive, pseudojet1, pseudojet2, MET);

  // return MT2
  return MT2;
}

// Calculate MT2 with minimized DHT  ----------------------------------------------------------------------------------------------------
Float_t MT2tree::GetMT2MinDHT(Float_t testmass, bool massive, Int_t PFJID, Float_t minJPt, Float_t maxJEta, Int_t met) {

  if(misc.CrazyHCAL ) return false; // protection against crazy HCAL events
  
  // fill Pseudojets with selected objects
  vector<TLorentzVector> p4s;
  for(int i=0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt, maxJEta, PFJID) ==false) continue;
	p4s.push_back(jet[i].lv);
  }
  for(int i=0; i<NEles; ++i){
	if(met==1113) continue; 
	p4s.push_back(ele[i].lv);
  }
  for(int i=0; i<NMuons; ++i){
	if(met==1113) continue; 
	p4s.push_back(muo[i].lv);
  }
  if (p4s.size()<2) return -999.99;
  

  // grouping according to minimal dHT
  std::vector<std::vector<float> > ht( 1<<(p4s.size()-1) , std::vector<float>( 2, 0.) ); 
      // initializiert einen vector ht der size 1<<(p4s.size()-1), wobei jeder eintrag ein vector (0,0) ist

  std::vector<std::vector<float> > px( 1<<(p4s.size()-1) , std::vector<float>( 2, 0.) ); 
  std::vector<std::vector<float> > py( 1<<(p4s.size()-1) , std::vector<float>( 2, 0.) ); 
  std::vector<std::vector<float> > pz( 1<<(p4s.size()-1) , std::vector<float>( 2, 0.) ); 
  std::vector<std::vector<float> > E ( 1<<(p4s.size()-1) , std::vector<float>( 2, 0.) ); 

  for(unsigned i=0; i < ht.size(); i++) {                                             
  for(unsigned j=0; j < p4s.size(); j++) {
		ht [i] [(i/(1<<j))%2] += p4s[j].Pt();
		px [i] [(i/(1<<j))%2] += p4s[j].Px();
		py [i] [(i/(1<<j))%2] += p4s[j].Py();
		pz [i] [(i/(1<<j))%2] += p4s[j].Pz();
		E  [i] [(i/(1<<j))%2] += p4s[j].E();
  }  
  }
  std::vector<float> deltaHT; for(unsigned i=0; i<ht.size(); i++) deltaHT.push_back(fabs(ht[i][0]-ht[i][1]));
  const float mDHT = *(std::min_element( deltaHT.begin(), deltaHT.end() ));
  int pos=distance(deltaHT.begin(), min_element(deltaHT.begin(), deltaHT.end()));
  TLorentzVector pseudojet1, pseudojet2;
  pseudojet1.SetPxPyPzE(px[pos][0], py[pos][0], pz[pos][0], E[pos][0]);
  pseudojet2.SetPxPyPzE(px[pos][1], py[pos][1], pz[pos][1], E[pos][1]);

  // define MET
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1){
	MET = pfmet[0];
  } else if(met==1113) { // add leptons to MET
	MET = pfmet[0]; 
	for(int i=0; i<NEles; ++i){
		MET = MET + ele[i].lv;
	}
	for(int i=0; i<NMuons; ++i){
		MET = MET + muo[i].lv;
	}
  } else if(met==4) {MET = -pseudojet1 - pseudojet2;}
  else return -999.99;

  Float_t MT2 = CalcMT2(testmass, massive, pseudojet1, pseudojet2, MET);

  // return MT2
  return MT2;
}

// Fill MT2 with hemispheres  ---------------------------------------------------------------------------------------------------------
Bool_t MT2tree::FillMT2Hemi(Float_t testmass, bool massive, Int_t PFJID, Float_t minJPt, Float_t maxJEta, Int_t hemi_seed, Int_t hemi_association, Int_t met, Int_t hemi_nr) {
  
  if(misc.CrazyHCAL ) return false; // protection against crazy HCAL events
  
  // struct to keep track of indices in hemi association
  struct HemiObj {vector<TString> type; vector<int> index; vector<int> hemisphere;
  } hemiobjs;

  // fill Pseudojets with selected objects
  vector<float> px, py, pz, E;
  for(int i=0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt, maxJEta, PFJID) ==false) continue;
  	px.push_back(jet[i].lv.Px());
	py.push_back(jet[i].lv.Py());
	pz.push_back(jet[i].lv.Pz());
	E .push_back(jet[i].lv.E ());
	hemiobjs.index.push_back(i); hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
  }
  for(int i=0; i<NEles; ++i){
	if(met==1113) continue; 
  	px.push_back(ele[i].lv.Px());
	py.push_back(ele[i].lv.Py());
	pz.push_back(ele[i].lv.Pz());
	E .push_back(ele[i].lv.E ());
	hemiobjs.index.push_back(i); hemiobjs.type.push_back("ele"), hemiobjs.hemisphere.push_back(0);
  }
  for(int i=0; i<NMuons; ++i){
	if(met==1113) continue; 
  	px.push_back(muo[i].lv.Px());
	py.push_back(muo[i].lv.Py());
	pz.push_back(muo[i].lv.Pz());
	E .push_back(muo[i].lv.E ());
	hemiobjs.index.push_back(i); hemiobjs.type.push_back("muo"), hemiobjs.hemisphere.push_back(0);
  }
  if (px.size()<2) return false;
  

  // get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)
  Hemisphere* hemisp = new Hemisphere(px, py, pz, E, hemi_seed, hemi_association);
  vector<int> grouping = hemisp->getGrouping();
  int NHemiIterations  = hemisp->GetNumLoop();

  TLorentzVector pseudojet1(0.,0.,0.,0.);
  TLorentzVector pseudojet2(0.,0.,0.,0.);
	
  float dHT=0;
  for(int i=0; i<px.size(); ++i){
	if(grouping[i]==1){
		pseudojet1.SetPx(pseudojet1.Px() + px[i]);
		pseudojet1.SetPy(pseudojet1.Py() + py[i]);
		pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
		pseudojet1.SetE( pseudojet1.E()  + E[i]);	
		dHT += sqrt(px[i]*px[i] + py[i]*py[i]);
		hemiobjs.hemisphere[i]=1;
	}else if(grouping[i] == 2){
		pseudojet2.SetPx(pseudojet2.Px() + px[i]);
		pseudojet2.SetPy(pseudojet2.Py() + py[i]);
		pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
		pseudojet2.SetE( pseudojet2.E()  + E[i]);
		dHT -= sqrt(px[i]*px[i] + py[i]*py[i]);
		hemiobjs.hemisphere[i]=2;
	}
  }
  delete hemisp;

  // define MET
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1){
	MET = pfmet[0];
  } else if(met==1113) { // add leptons to MET
	MET = pfmet[0]; 
	for(int i=0; i<NEles; ++i){
		MET = MET + ele[i].lv;
	}
	for(int i=0; i<NMuons; ++i){
		MET = MET + muo[i].lv;
	}
  } else if(met==4) {MET = -pseudojet1 - pseudojet2;}
  else return false;

  Float_t MT2 = CalcMT2(testmass, massive, pseudojet1, pseudojet2, MET);

  // fill MT2Hemi
  if(hemi_nr!=-1 && hemi_nr < m_hemiSize){
	  hemi[hemi_nr].seed_method  = hemi_seed;
	  hemi[hemi_nr].assoc_method = hemi_association;
	  hemi[hemi_nr].NHemiIter    = NHemiIterations;
	  hemi[hemi_nr].MT2          = MT2;
	  hemi[hemi_nr].dPhi         = Util::DeltaPhi(pseudojet1.Phi(), pseudojet2.Phi());
	  hemi[hemi_nr].minDHT       = dHT;
	  hemi[hemi_nr].maxDR        = -999.99; // option not implemented
	  hemi[hemi_nr].AlphaT       = -999.99; // only implemented for minimized DHT
	  hemi[hemi_nr].UTM          = - MET - pseudojet1 - pseudojet2;
	  hemi[hemi_nr].lv1          = pseudojet1;
	  hemi[hemi_nr].lv2          = pseudojet2;
	  // fill MCT
	  TVector2 pmiss_vector2;
	  pmiss_vector2.Set(MET.Px(), MET.Py());
	  TLorentzVector downstream(0.,0.,0.,0.); // no particles are downstream, i.e. not selected jets are upstream. 
	  hemi[hemi_nr].MCT          = GetMCTcorr(pseudojet1, pseudojet2, downstream, pmiss_vector2);
	  // fill indices for hemisphere association
	  int jcount1=0, jcount2=0, elecount1=0, elecount2=0, muocount1=0, muocount2=0; 
	  for(int i=0; i<hemiobjs.index.size();++i){
	  	if(hemiobjs.type[i]=="jet" && hemiobjs.hemisphere[i]==1){
			hemi[hemi_nr].jindices1[jcount1]=hemiobjs.index[i];
			jcount1++;
		}else if(hemiobjs.type[i]=="jet" && hemiobjs.hemisphere[i]==2){
			hemi[hemi_nr].jindices2[jcount2]=hemiobjs.index[i];
			jcount2++;
		}else if(hemiobjs.type[i]=="ele" && hemiobjs.hemisphere[i]==1){
			hemi[hemi_nr].eleindices1[elecount1]=hemiobjs.index[i];
			elecount1++;
		}else if(hemiobjs.type[i]=="ele" && hemiobjs.hemisphere[i]==2){
			hemi[hemi_nr].eleindices2[elecount2]=hemiobjs.index[i];
			elecount2++;
		}else if(hemiobjs.type[i]=="muo" && hemiobjs.hemisphere[i]==1){
			hemi[hemi_nr].muoindices1[muocount1]=hemiobjs.index[i];
			muocount1++;
		}else if(hemiobjs.type[i]=="muo" && hemiobjs.hemisphere[i]==2){
			hemi[hemi_nr].muoindices2[muocount2]=hemiobjs.index[i];
			muocount2++;
		}
	  }
	  return true;
  }else return false;
}

// fill MT2 for minimized Delta_HT (RA1 way of doing things)
Bool_t MT2tree::FillMT2HemiMinDHT(Float_t testmass, bool massive, Int_t PFJID, Float_t minJPt, Float_t maxJEta, Int_t met, Int_t hemi_nr) {

  if(misc.CrazyHCAL ) return false; // protection against crazy HCAL events
  
  // fill Pseudojets with selected objects
  vector<TLorentzVector> p4s;
  for(int i=0; i<NJets; ++i){
	if(jet[i].IsGoodPFJet(minJPt, maxJEta, PFJID) ==false) continue;
	p4s.push_back(jet[i].lv);
  }
  for(int i=0; i<NEles; ++i){
	if(met==1113) continue; 
	p4s.push_back(ele[i].lv);
  }
  for(int i=0; i<NMuons; ++i){
	if(met==1113) continue; 
	p4s.push_back(muo[i].lv);
  }
  if (p4s.size()<2) return false;
  
  // grouping according to minimal dHT
  std::vector<std::vector<float> > ht( 1<<(p4s.size()-1) , std::vector<float>( 2, 0.) ); 
      // initializiert einen vector ht der size 1<<(p4s.size()-1), wobei jeder eintrag ein vector (0,0) ist

  std::vector<std::vector<float> > px( 1<<(p4s.size()-1) , std::vector<float>( 2, 0.) ); 
  std::vector<std::vector<float> > py( 1<<(p4s.size()-1) , std::vector<float>( 2, 0.) ); 
  std::vector<std::vector<float> > pz( 1<<(p4s.size()-1) , std::vector<float>( 2, 0.) ); 
  std::vector<std::vector<float> > E ( 1<<(p4s.size()-1) , std::vector<float>( 2, 0.) ); 

  for(unsigned i=0; i < ht.size(); i++) {                                             
  for(unsigned j=0; j < p4s.size(); j++) {
		ht [i] [(i/(1<<j))%2] += p4s[j].Pt();
		px [i] [(i/(1<<j))%2] += p4s[j].Px();
		py [i] [(i/(1<<j))%2] += p4s[j].Py();
		pz [i] [(i/(1<<j))%2] += p4s[j].Pz();
		E  [i] [(i/(1<<j))%2] += p4s[j].E();
  }  
  }
  std::vector<float> deltaHT; for(unsigned i=0; i<ht.size(); i++) deltaHT.push_back(fabs(ht[i][0]-ht[i][1]));
  const float mDHT = *(std::min_element( deltaHT.begin(), deltaHT.end() )); // get minimal configuration
  int pos=distance(deltaHT.begin(), min_element(deltaHT.begin(), deltaHT.end()));
  TLorentzVector pseudojet1, pseudojet2;
  pseudojet1.SetPxPyPzE(px[pos][0], py[pos][0], pz[pos][0], E[pos][0]);
  pseudojet2.SetPxPyPzE(px[pos][1], py[pos][1], pz[pos][1], E[pos][1]);

  // calculate alpha_T  
  std::vector<float> pTs; for(unsigned i=0; i<p4s.size(); i++) pTs.push_back(p4s[i].Pt());
  const float sumPT = accumulate( pTs.begin(), pTs.end(), float(0) );
  const TLorentzVector sumP4 = accumulate( p4s.begin(), p4s.end(), TLorentzVector() );
  float alphaT  = 0.5 * ( sumPT - mDHT ) / sqrt( sumPT*sumPT - sumP4.Perp2() );
  
  // define MET
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1){
	MET = pfmet[0];
  } else if(met==1113) { // add leptons to MET
	MET = pfmet[0]; 
	for(int i=0; i<NEles; ++i){
		MET = MET + ele[i].lv;
	}
	for(int i=0; i<NMuons; ++i){
		MET = MET + muo[i].lv;
	}
  } else if(met==4) {MET = -pseudojet1 - pseudojet2;}
  else return false;

  // fill MT2Hemi
  if(hemi_nr!=-1 && hemi_nr < m_hemiSize){
	  hemi[hemi_nr].seed_method  = -1;
	  hemi[hemi_nr].assoc_method = -1;
	  hemi[hemi_nr].NHemiIter    = -1; 
	  hemi[hemi_nr].MT2          = CalcMT2(testmass, massive, pseudojet1, pseudojet2, MET);
	  hemi[hemi_nr].dPhi         = Util::DeltaPhi(pseudojet1.Phi(), pseudojet2.Phi());
	  hemi[hemi_nr].minDHT       = mDHT;
	  hemi[hemi_nr].maxDR        = -999.99; // option not implemented
	  hemi[hemi_nr].AlphaT       = alphaT;  // only implemented for minimized DHT
	  hemi[hemi_nr].UTM          = - MET - pseudojet1 - pseudojet2;
	  hemi[hemi_nr].lv1          = pseudojet1;
	  hemi[hemi_nr].lv2          = pseudojet2;
	  // fill MCT
	  TVector2 pmiss_vector2;
	  pmiss_vector2.Set(MET.Px(), MET.Py());
	  TLorentzVector downstream(0.,0.,0.,0.); // no particles are downstream, i.e. not selected jets are upstream. 
	  hemi[hemi_nr].MCT          = GetMCTcorr(pseudojet1, pseudojet2, downstream, pmiss_vector2);
	  return true;
  }else return false;
}

// calculate MT2 ---------------------------------------------------------------------------------------------
Float_t MT2tree::CalcMT2(float testmass, bool massive, TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET ){
  
  double pa[3];
  double pb[3];
  double pmiss[3];
  
  pmiss[0] = 0;
  pmiss[1] = static_cast<double> (MET.Px());
  pmiss[2] = static_cast<double> (MET.Py());
  
  pa[0] = static_cast<double> (massive ? visible1.M() : 0);
  pa[1] = static_cast<double> (visible1.Px());
  pa[2] = static_cast<double> (visible1.Py());
  
  pb[0] = static_cast<double> (massive ? visible2.M() : 0);
  pb[1] = static_cast<double> (visible2.Px());
  pb[2] = static_cast<double> (visible2.Py());
  
  Davismt2 *mt2 = new Davismt2();
  mt2->set_momenta(pa, pb, pmiss);
  mt2->set_mn(testmass);
  Float_t MT2=mt2->get_mt2();
  delete mt2;
  return MT2;

}

Float_t MT2tree::SimpleMT2(bool pseudo, int heminr){
  if (!pseudo){
    if(NJets < 2) return -999.99;
    return sqrt(2*jet[0].lv.Pt()*jet[1].lv.Pt()*(1+TMath::Cos(jet[0].lv.DeltaPhi(jet[1].lv))));
  }
  else{
    return sqrt(2*hemi[heminr].lv1.Pt()*hemi[heminr].lv2.Pt()*(1+TMath::Cos(hemi[heminr].lv1.DeltaPhi(hemi[heminr].lv2))));
  }
}


// transverse mass MT ----------------------------------------------------------------------------------
Float_t MT2tree::GetMT(TLorentzVector lv1, float m1, TLorentzVector lv2, float m2){
	// returns Mt. Not not rely on Pz measurement -> also suitable if one lv == MET  
	float ET_1 = sqrt(m1*m1 + lv1.Perp2());
	float ET_2 = sqrt(m2*m2 + lv2.Perp2());
	float MTsquared = m1*m1 + m2*m2 + 2*ET_1*ET_2 - 2*lv1.Px()*lv2.Px() - 2*lv1.Py()*lv2.Py();
	if(MTsquared < 0) return -sqrt(MTsquared);
	else              return  sqrt(MTsquared);
}

Float_t MT2tree::GetMT(TLorentzVector lv1, TLorentzVector lv2){
	return GetMT(lv1, 0., lv2, 0.);
}

// cotransverse mass MCT  -----------------------------------------------------------------------------------
Float_t MT2tree::GetMCTcorr(TLorentzVector p1, TLorentzVector p2, TLorentzVector DTM, TVector2 pmiss){
	// Tovey MCT corrected
	TMctLib* MCT = new TMctLib;
	float MCTcorr = MCT -> mctcorr(p1, p2, DTM, pmiss, 7000., 0.);
	delete MCT;
	return MCTcorr;
}

// sqrt s-hat min ---------------------------------------------------------------------------
Float_t MT2tree::GetSqrtS(float testmass, bool massive, int PFJID, float minJPt, float maxJEta,int met){
  TLorentzVector MET(0., 0., 0., 0.);
  if(met==1)      MET = pfmet[0];
  else if(met==2) MET = MHT[0];
  else            return -999;

  TLorentzVector obj(0,0,0,0);
  TLorentzVector vis(0,0,0,0);
  for(int i = 0; i<NJets; ++i){
    if(jet[i].IsGoodPFJet(minJPt,maxJEta,PFJID)==false) continue;
    obj.SetPtEtaPhiM(jet[i].lv.Pt(),jet[i].lv.Eta(),jet[i].lv.Phi(),massive ? jet[i].lv.M() : 0);
    vis+=obj;
  }
  for(int i = 0; i<NEles; ++i){
    obj.SetPtEtaPhiM(ele[i].lv.Pt(),ele[i].lv.Eta(),ele[i].lv.Phi(),ele[i].lv.M());
    vis+=obj;
  }
  for(int i = 0; i<NMuons; ++i){
    obj.SetPtEtaPhiM(muo[i].lv.Pt(),muo[i].lv.Eta(),muo[i].lv.Phi(),muo[i].lv.M());
    vis+=obj;
  }
  return sqrt(TMath::Power(vis.M(),2)+TMath::Power(vis.Pt(),2)) + sqrt(TMath::Power(2*testmass,2)+TMath::Power(MET.Pt(),2));
}



// -----------------------------------------------------------------------------------------------------------
// Leptons and friends

Float_t MT2tree::GenOSDiLeptonInvMass(unsigned int pid, unsigned int mother, float pt, float eta){
	 // returns true if there are 
	 //  - two particles with abs(ID)=pid
	 //  - they have opposite charge
	 //  - they come from "mother"
	 //  - they have Pt() > pt and |Eta()| < eta
	 //  - their inv mass is in (lower_Minv, upper_Minv)
	vector<int> indices;	
	for(int i=0; i<NGenLepts; ++i){
		if(abs(genlept[i].ID)       !=pid   ) continue;
		if(abs(genlept[i].MID)      !=mother) continue;
		if(    genlept[i].lv.Pt()   < pt    ) continue;
		if(fabs(genlept[i].lv.Eta())> eta   ) continue;
		indices.push_back(i);
	}
	if(indices.size()!=2)                                        return -100;
	if( (genlept[indices[0]].ID) * (genlept[indices[1]].ID) >=0 && pid !=12 && pid != 14 && pid !=16) return -150;
	
	float Mass= (genlept[indices[0]].lv + genlept[indices[1]].lv).M();
	if(Mass<0) Mass = -Mass;
	return Mass;

}


Bool_t   MT2tree::IsGenOSDiLepton(unsigned int pid, unsigned int mother, float pt, float eta, float lower_mass, float upper_mass){
	 float Mass = GenOSDiLeptonInvMass(pid, mother, pt, eta);
	 if(Mass < 0 ) return false;
	 if(Mass < lower_mass || Mass > upper_mass) return false;
	 return true;
}


Float_t MT2tree::GetDiLeptonInvMass(int same_sign, int same_flavour, int flavour, float pt, bool exclDiLept){
	// flavour == 0 : don't care if el or mu
	// flavour == 1 : only electrons
	// flavour == 2 : only muons
	
	struct lepton {TLorentzVector lv; int charge; string flavour;} lepts[NEles + NMuons];	
	int counter = 0;
	for(int i=0; i<NEles; ++i){
		lepts[counter].lv=ele[i].lv; lepts[counter].charge =ele[i].Charge; lepts[counter].flavour ="ele";
		counter++;
	}
	for(int i=0; i<NMuons; ++i){
		lepts[counter].lv=muo[i].lv; lepts[counter].charge =muo[i].Charge; lepts[counter].flavour ="muo"; 
		counter++;
	}
	
	int index1=-1, index2=-1;
	if( (NEles + NMuons) <= 1)                                               return -10;
	if( (NEles + NMuons) >  2 && exclDiLept)                                 return -100;
	else if( (NEles + NMuons) >  2 && !exclDiLept){
		float pt1=0, pt2=0;
		for(int i=0; i<(NEles + NMuons); ++i){
			if(lepts[i].lv.Pt()>pt1){
				pt2 =pt1;
				index2=index1;
				pt1 =lepts[i].lv.Pt();
				index1=i;
			} else if(lepts[i].lv.Pt()>pt2){
				index2=i;
				pt2   =lepts[i].lv.Pt();
			}
		}	
	
	}else if ((NEles + NMuons)==2) {
		index1=0; index2=1; 
	}

	if(  lepts[index1].charge * lepts[index2].charge ==  1 && same_sign ==0)             return -200;
	if(  lepts[index1].charge * lepts[index2].charge == -1 && same_sign ==1)             return -300;
	if(  lepts[index1].lv.Pt() < pt || lepts[index2].lv.Pt() < pt )                      return -400;
	if(  lepts[index1].flavour != lepts[index2].flavour && same_flavour == 1 )           return -500;
	if(  lepts[index1].flavour == lepts[index2].flavour && same_flavour == 0 )           return -600;
	if( (lepts[index1].flavour =="muo" || lepts[index2].flavour=="muo") && flavour == 1) return -910;
	if( (lepts[index1].flavour =="ele" || lepts[index2].flavour=="ele") && flavour == 2) return -920;

	TLorentzVector sum = lepts[index1].lv + lepts[index2].lv;
	float mass = sum.M();
	if(mass < 0) return -mass;
	else return mass;
}

Bool_t MT2tree::IsDiLeptonMll(int same_sign, int same_flavour, int flavour, float pt, bool exclDiLept, float lower_mass, float upper_mass){
	float mass = GetDiLeptonInvMass(same_sign, same_flavour, flavour, pt, exclDiLept);
	if(mass < 0) return false;
	if(mass < lower_mass ||  mass > upper_mass) return false;
	return true;
}


TLorentzVector MT2tree::GetMETPlusLeptsLV(int OSDiLeptFromZ){
	TLorentzVector lv = pfmet[0];
	if(OSDiLeptFromZ ==1) {
		float mass = GetDiLeptonInvMass(0,1,0,5.0,1);
		if(mass < 0 )               return lv;
		if(mass < 71 || mass > 111) return lv; 
	}
	for(int i=0; i<NEles; ++i){
		lv+=ele[i].lv;
	}
	for(int i=0; i<NMuons; ++i){
		lv+=muo[i].lv;
	}
	return lv;
}

Float_t MT2tree::GetMETPlusLepts(int OSDiLeptFromZ){
	if(OSDiLeptFromZ ==1) {
		float mass = GetDiLeptonInvMass(0,1,0,5.0,1);
		if(mass < 0 )               return -2000;
		if(mass < 71 || mass > 111) return -1000; 
	}
	TLorentzVector lv = pfmet[0];
	for(int i=0; i<NEles; ++i){
		lv+=ele[i].lv;
	}
	for(int i=0; i<NMuons; ++i){
		lv+=muo[i].lv;
	}
	return lv.Pt();
}

Float_t MT2tree::GetDiLeptonPt(int same_sign, int same_flavour, int flavour, float pt, float lower_mass, float upper_mass){
	float mass = GetDiLeptonInvMass(same_sign, same_flavour, flavour, pt, 1);
	if(mass < 0 )                              return -2000;
	if(mass < lower_mass || mass > upper_mass) return -1000; 
	TLorentzVector lv(0., 0., 0., 0.);
	for(int i=0; i<NEles; ++i){
		lv+=ele[i].lv;
	}
	for(int i=0; i<NMuons; ++i){
		lv+=muo[i].lv;
	}
	return lv.Pt();
}

Float_t MT2tree::GetMETPlusGenLepts(int met, int RemoveOSSFDiLepts, int require_cuts,  unsigned int pid, unsigned int mother, float pt, float eta, float lower_mass, float upper_mass ){

	TLorentzVector lv;
        if(met==1) lv= genmet[0];
        if(met==0) lv= pfmet[0];


	vector<int> indices;	
	if(RemoveOSSFDiLepts ==1 || require_cuts ==1 ) {
		for(int i=0; i<NGenLepts; ++i){
			if(pid==1113) {if(abs(genlept[i].ID) !=11 && abs(genlept[i].ID)!=13  ) continue;}
			if(pid!=1113) {if(abs(genlept[i].ID) !=pid   ) continue;}
			if(abs(genlept[i].MID)               !=mother) continue;
			if(    genlept[i].lv.Pt()            < pt    ) continue;
			if(fabs(genlept[i].lv.Eta())         > eta   ) continue;
			indices.push_back(i);
		}
		if(indices.size()!=2)                                         return -100;
		if( (genlept[indices[0]].ID) * (genlept[indices[1]].ID) >=0 ) return -150;
		
		float Mass= (genlept[indices[0]].lv + genlept[indices[1]].lv).M();
		if(Mass<0) Mass = -Mass;
		if(Mass > upper_mass || Mass < lower_mass)                    return -200;
		if(RemoveOSSFDiLepts ==1) {
			lv+=genlept[indices[0]].lv;
			lv+=genlept[indices[1]].lv;
		}
	}
	
	return lv.Pt();
}

Int_t   MT2tree::GetGenLeptIndex(int which, int pid, int mother, float pt, float eta){
	vector<int>    indices;
	vector<float> pts, etas;
	for(int i=0; i<NGenLepts; ++i){
		if(pid==1113)        {if(abs(genlept[i].ID) !=11 && abs(genlept[i].ID)!=13  ) continue;}
		else if(pid==121416) {if(abs(genlept[i].ID) !=12 && abs(genlept[i].ID)!=14 && abs(genlept[i].ID)!=16) continue;}
		else if(abs(genlept[i].ID) !=pid   ) continue;
		if(abs(genlept[i].MID)               !=mother) continue;
		if(    genlept[i].lv.Pt()            < pt    ) continue;
		if(fabs(genlept[i].lv.Eta())         > eta   ) continue;
		indices.push_back(i);
		pts.push_back(genlept[i].lv.Pt());
	}
	if(indices.size() < 1 || which >= indices.size()) return -1;
	else indices = Util::VSort(indices, pts);

	return (Int_t) indices[which]; 	
}

Float_t MT2tree::GetGenLeptEta(int which, int pid, int mother, float pt, float eta){
	Int_t index = GetGenLeptIndex(which, pid, mother, pt, eta);	
	if(index ==-1) return -999.99;
	else           return genlept[index].lv.Eta();	
}

Float_t MT2tree::GetGenLeptPt(int which, int pid, int mother, float pt, float eta){
	Int_t index = GetGenLeptIndex(which, pid, mother, pt, eta);	
	if(index ==-1) return -999.99;
	else           return genlept[index].lv.Pt();	
}

Bool_t MT2tree::GenLeptFromW(int pid, float pt, float eta, bool includeTau){
	bool good(false);
	for(int i=0; i<NGenLepts; ++i){
		if(abs(genlept[i].ID) !=pid                         ) continue;
		if( (!includeTau) && abs(genlept[i].MID) !=24       ) continue;
		if(   includeTau  && !((abs(genlept[i].MID)==15 && abs(genlept[i].GMID)==24 ) || abs(genlept[i].MID)==24)) continue;
		if(genlept[i].lv.Pt()       < pt                    ) continue;
		if(fabs(genlept[i].lv.Eta())>eta                    ) continue;
		good=true;
	}
	return good;	
}

Float_t MT2tree::GetLeptPt(int index){
	if(NEles + NMuons ==0) return -999;
	float pt_0 =0;
	float pt_1 =0;
	for(int i=0; i<NEles; ++i){
		if(ele[i].lv.Pt() > pt_0)                          {pt_1 = pt_0; pt_0 = ele[i].lv.Pt();}
		if(ele[i].lv.Pt() < pt_0 && ele[i].lv.Pt() > pt_1) {pt_1 = ele[i].lv.Pt();}
	}
	for(int i=0; i<NMuons; ++i){
		if(muo[i].lv.Pt() > pt_0)                          {pt_1 = pt_0; pt_0 = muo[i].lv.Pt();}
		if(muo[i].lv.Pt() < pt_0 && muo[i].lv.Pt() > pt_1) {pt_1 = muo[i].lv.Pt();}
	}

	if     (index==0) return pt_0;
	else if(index==1) return pt_1;
	else return -1;
}

Float_t MT2tree::ElClosestJet(){
	float dR=1000;
	for(int i=0; i<NEles; ++i){
		for(int j=0; j<NJets; ++j){
			if(ele[i].lv.DeltaR(jet[j].lv) < dR) {dR=ele[i].lv.DeltaR(jet[j].lv);}
		}	
	}
	return dR;
}

Int_t MT2tree::TopDecayMode(){
	// bit map: 
	// 1 = electron 1
	// 2 = electron 2
	// 4 = muon 1
	// 8 = muon 2
	// 16 = tau 1
	// 32 = tau 2
	// 64 = leptonic tau1
	// 128= leptonic tau2
	Int_t bit=0;
	Bool_t acceptance(true);
	for(int i=0; i<NGenLepts; ++i){
		if( abs(genlept[i].ID)==11 && abs(genlept[i].MID)==24 && abs(genlept[i].GMID)==6 ) {
			if( (bit & 1)==0) bit = bit | 1;
			else              bit = bit | 2;
		} 
		if( abs(genlept[i].ID)==13 && abs(genlept[i].MID)==24 && abs(genlept[i].GMID)==6 ) {
			if( (bit & 4)==0) bit = bit | 4;
			else              bit = bit | 8;
		} 
		if( abs(genlept[i].ID)==11 && abs(genlept[i].MID)==15 && abs(genlept[i].GMID)==24 ) {
			if     ( (bit & 16)==0) bit = bit | 16;
			else                    bit = bit | 32;
			if     ( (bit & 1 )==0) bit = bit |  1;
			else                    bit = bit |  2;
			if     ( (bit & 64)==0) bit = bit | 64;
			else                    bit = bit |128;
		} 
		if( abs(genlept[i].ID)==13 && abs(genlept[i].MID)==15 && abs(genlept[i].GMID)==24 ) {
			if     ( (bit & 16)==0) bit = bit | 16;
			else                    bit = bit | 32;
			if     ( (bit & 4 )==0) bit = bit |  4;
			else                    bit = bit |  8;
			if     ( (bit & 64)==0) bit = bit | 64;
			else                    bit = bit |128;
		}
		if( abs(genlept[i].ID)==16 && abs(genlept[i].MID)==24 && abs(genlept[i].GMID)==6 ){
			if     ( (bit & 16)==0) bit = bit | 16;
			else                    bit = bit | 32;
		}
	}
	return bit;
}


Bool_t MT2tree::TopDecayModeResult(Int_t nlepts){
	Int_t bit =TopDecayMode();
	if(nlepts == 1){ // semileptonic without leptonic tau
		if     ( (bit & 2 )==2 || (bit & 8)==8) return false; // more than one e/mu
		if     ( (bit & 64)==64)                return false; // at least one leptonic tau
		if     ( (bit & 1 )==1 && (bit & 4)==0) return true;
		else if( (bit & 1 )==0 && (bit & 4)==4) return true;
		else                                    return false;
	}else if(nlepts == 115){ // semileptonic with leptonic tau
		if     ( (bit & 2 )==2 || (bit & 8)==8) return false; // more than one e/mu
		if     ( (bit & 1 )==1 && (bit & 4)==0) return true;
		else if( (bit & 1 )==0 && (bit & 4)==4) return true;
		else                                    return false;
	}else if(nlepts == 215){ // fully leptonic with leptonic tau
		if     ( (bit & 2 )==2 || (bit & 8)==8) return true; // two eles or two muons
		if     ( (bit & 1 )==1 && (bit & 4)==4) return true; // one ele and one muo
		else                                    return false;
	}else if(nlepts == 2){ // fully leptonic without leptonic tau
		if     ( (bit & 64)==64               ) return false; // leptonic tau
		if     ( (bit & 2 )==2 || (bit & 8)==8) return true; // two eles or two muons
		if     ( (bit & 1 )==1 && (bit & 4)==4) return true; // one ele and one muo
		else                                    return false;
	}else if(nlepts == 0){ // fully hadronic without hadronic tau
		if     ( (bit & 1 )==1 || (bit & 4)==4) return false; // ele or muo
		if     ( (bit & 16)==16               ) return false; // tau
		else                                    return true;
	}else if(nlepts == 15){ // fully hadronic with hadronic tau
		if     ( (bit & 1 )==1 || (bit & 4)==4) return false; // ele or muo
		else                                    return true;
	}else if(nlepts ==11){
		if     ( (bit & 4 )==4 ) return false; // muon
		if     ( (bit & 2 )==2 ) return false; // two electron
		if     ( (bit & 1 )==1 ) return true;  // electron
		else                     return false;
	}else if(nlepts ==13){
		if     ( (bit & 1 )==1 ) return false; // ele
		if     ( (bit & 8 )==8 ) return false; // two muons
		if     ( (bit & 4 )==4 ) return true;  // muon
		else                     return false;
	}
	else                                           return false;
}

Bool_t MT2tree::SLTopAccept(float pt, float eta){
	for(int i=0; i<NGenLepts; ++i){
		if(    abs(genlept[i].ID)  !=11 && abs(genlept[i].ID) !=13                               ) continue;
		if( ! (abs(genlept[i].MID) ==15 && abs(genlept[i].GMID)==24 || abs(genlept[i].MID) ==24 )) continue;
		if(genlept[i].lv.Pt()>pt && fabs(genlept[i].lv.Eta()) < eta)                              return true;	
	}
	return false;
}

Float_t MT2tree::SLTopEta(float pt){
	for(int i=0; i<NGenLepts; ++i){
		if(    abs(genlept[i].ID)  !=11 && abs(genlept[i].ID) !=13                               ) continue;
		if( ! (abs(genlept[i].MID) ==15 && abs(genlept[i].GMID)==24 || abs(genlept[i].MID) ==24 )) continue;
		if(genlept[i].lv.Pt()>pt)      return genlept[i].lv.Eta();	
	}
	return -999.99;
}


Int_t MT2tree::WDecayMode(){
	// bit map:
	// 0 not recognized
	// 1= ele
	// 2= muo
	// 4= tau
	// 8= tau stable (i.e. problem in MC sample)
	Int_t result =0;
	for(int i=0; i<NGenLepts; ++i){
		if( abs(genlept[i].ID)==11 && abs(genlept[i].MID)==24 )                               {result = result | 1;} 
		if( abs(genlept[i].ID)==13 && abs(genlept[i].MID)==24 )                               {result = result | 2;} 
		if( abs(genlept[i].ID)==16 && abs(genlept[i].MID)==24 )                               {result = result | 4;}  // tau neutrino
		if( abs(genlept[i].ID)==11 && abs(genlept[i].MID)==15 && abs(genlept[i].GMID)==24 )   {result = result | 5;} 
		if( abs(genlept[i].ID)==13 && abs(genlept[i].MID)==15 && abs(genlept[i].GMID)==24 )   {result = result | 6;} 
		if( abs(genlept[i].ID)==15 && abs(genlept[i].MID)==24 )                               {result = result | 8;} // stable tau
	}
	return result;
}

Float_t MT2tree::LeptJetDR( int pid, int index, bool bjet, int ID){
        if(pid==11 && NEles  < index+1) return -1;
	if(pid==13 && NMuons < index+1) return -1;
	Float_t dRmin=999.99;
	for(int i=0; i<NJets; ++i){
		if(! jet[i].IsGoodPFJet(20,2.8,ID))         continue;
		if(bjet && jet[i].bTagProbSSVHP < 2.0)      continue;
		Float_t dR=999.99;
		if(pid==13) dR = jet[i].lv.DeltaR(muo[index].lv);
		if(pid==11) dR = jet[i].lv.DeltaR(ele[index].lv);
		if(dR < dRmin) dRmin=dR;
	}
	return dRmin;
}

Float_t MT2tree::GetGenVPt(int pid){
  TLorentzVector V_p;
  int countNLepts=0;
  for(int i=0; i<NGenLepts; ++i){
    if( abs(genlept[i].MID) != pid       ) continue;
    V_p += genlept[i].lv;
    countNLepts++;
  }
  //if(countNLepts!=2)cout << "[WARNING]: " << countNLepts << " leptons for this boson" << endl;
  return V_p.Pt();
}



// Print-Outs ---------------------------------------------------------------------------------------------------------
Bool_t MT2tree::PrintOut(Bool_t logfile){
	
	std::ostringstream logStream;
	// detailed Event Printout
	logStream << "********************************************************************************************* " << endl;
	logStream << "Event " << misc.Event   << " Lumi " << misc.LumiSection << " run " <<  misc.Run                 << endl;
	logStream << "  NJetsIDLoose (pT > 20, |eta| < 2.4, Pf-JID) " << NJetsIDLoose                                 << endl;
	logStream << "  NEles " << NEles  << ", NMuons "<< NMuons                                                     << endl;
	logStream << "  NVertices " << pileUp.NVertices                                                               << endl;
        logStream << "  pf-HT " << misc.HT << ", caloHT50_ID " << misc.caloHT50_ID                                    << endl;
	logStream << "  pf-MET Pt:" << misc.MET << " Phi " << pfmet[0].Phi()                                          << endl;	
	if(!misc.isData){
	logStream << "  gen-met : " << genmet[0].Pt() << " phi " << genmet[0].Phi()                                   << endl;
	}
	logStream << " Data quality ------------------------------------------------------------------------------- "  << endl;
	logStream << "  HBHENoiseFlag " << misc.HBHENoiseFlag << " (1=good),  CrazyHCAL " << misc.CrazyHCAL 
	     << ", BadEcalTP " << misc.BadEcalTP         << ", BadEcalBE " << misc.BadEcalBE 
	     << ", CSCTightHaloID " << misc.CSCTightHaloID                                                        << endl;
        logStream << "  Jet0Pass " << misc.Jet0Pass <<" Jet1Pass " << misc.Jet1Pass << " PassJetID " << misc.PassJetID   << endl;	
	logStream << "  MinMetJetDPhi " << misc.MinMetJetDPhi                                                          << endl;
	logStream << " Jet-Info -------------------------------------------------------------------------------------" << endl;
	for(int i=0; i<NJets; ++i){
	logStream << "  jet " << i << ":\n"
	     << "   pt " << jet[i].lv.Pt() << ", eta " << jet[i].lv.Eta() << ", phi " <<   jet[i].lv.Phi() << ", E " << jet[i].lv.E() << ", Mass " << jet[i].lv.M() << "\n"	
	     << "   px " << jet[i].lv.Px() << ", py "  << jet[i].lv.Py()  << ", pz "  <<   jet[i].lv.Pz() << "\n"
	     << "   JID " << jet[i].isPATPFIDLoose << ", isTau " << jet[i].isTauMatch    << " Flavour " << jet[i].Flavour <<"\n"
	     << "   CHF " << jet[i].ChHadFrac << ", NHF " << jet[i].NeuHadFrac << ", CEF " << jet[i].ChEmFrac 
	                  << ", NEF " << jet[i].NeuEmFrac << ", NConstituents " << jet[i].NConstituents  
	                  << ", Ch Mult " << jet[i].ChMult << ", Neu Mul " << jet[i].NeuMult << "\n"
	     << "   isBtag SSVHP: " << (jet[i].IsBJet(3)? "true (":"false (") << jet[i].bTagProbSSVHP 
	     << "), isBtag SSVHE: " << (jet[i].IsBJet(2)? "true (":"false (") << jet[i].bTagProbSSVHE << ")\n"
	     << "   isBtag TCHEM: " << (jet[i].IsBJet(0)? "true (":"false (") << jet[i].bTagProbTCHE 
	     << "), isBtag TCHPM: " << (jet[i].IsBJet(1)? "true (":"false (") << jet[i].bTagProbTCHP << ")\n"
	     << "   L1FastL2L3 JE corr factor " << jet[i].Scale << ", L1Fast factor " << jet[i].L1FastJetScale << ", jet Area " <<  jet[i].Area << "\n"
	     << "   jet-MET-dPhi " << MetJetDPhi(i, 0, 1)  
	     << endl; 
	}	
	if(NEles >0){
	logStream << " Eles Info ------------------------------------------------------------------------------------" << endl;
	for (int i=0; i<NEles; ++i){
	logStream << "   Ele " << i << ":\n";
	logStream << "    Pt " << ele[i].lv.Pt() << " Eta " << ele[i].lv.Eta() << " Phi " << ele[i].lv.Phi() << " E " << ele[i].lv.E() << endl;
	logStream << "    Isolation " << ele[i].Iso                                                                     << endl;
	logStream << "    Charge    " << ele[i].Charge                                                                  << endl;
	logStream << "    VBTF ID 95 " << ele[i].ID95  << ", ID90 " << ele[i].ID90                                      << endl;
	logStream << "    transverse Mass with MET " << ele[i].MT                                                       << endl; 
	}
	logStream << "   Ele-Jet combinations -------------------------------------------------------------------------"   << endl;
	for (int i=0; i<NEles; ++i){
	for (int j=0; j<NJets;  ++j){
	TLorentzVector pj=jet[j].lv + ele[i].lv;
	logStream << "     Ele " << i << " and jet " << j << ": Pt " << pj.Pt() << " Px " << pj.Px() << " Py " << pj.Py() << " Pz " << pj.Pz()
		  << " Eta " << pj.Eta() << " Phi " << pj.Phi() << " E " << pj.E() << " M " << pj.M() << endl;	
	}
	}
	if(NMuons>1){
	logStream << "   Ele-Ele combinations --------------------------------------------------------------------------"   << endl;
	for(int i=0; i<NEles; ++i){
	for(int j=i+1; j<NEles; ++j){
	TLorentzVector pj=ele[i].lv +ele[j].lv;
	logStream << "     Ele " << i << " and ele " << j << ": Pt " << pj.Pt() << " Px " << pj.Px() << " Py " << pj.Py() << " Pz " << pj.Pz()
		  << " Eta " << pj.Eta() << " Phi " << pj.Phi() << " E " << pj.E() << " M " << pj.M() << endl;	
	}
	}
	}
	}
	if(NMuons >0){
	logStream << " Muons Info -----------------------------------------------------------------------------------"   << endl;
	for (int i=0; i<NMuons; ++i){
	logStream << "   Muon " << i << ":\n";
	logStream << "    Pt " << muo[i].lv.Pt() << " Eta " << muo[i].lv.Eta() << " Phi " << muo[i].lv.Phi() << " E " << muo[i].lv.E() << endl;
	logStream << "    Isolation " << muo[i].Iso                                                                      << endl;
	logStream << "    Charge    " << muo[i].Charge                                                                   << endl;
	logStream << "    transverse Mass with MET " << muo[i].MT                                                        << endl; 
	}
	logStream << "   Muon-Jet combinations ------------------------------------------------------------------------"   << endl;
	for (int i=0; i<NMuons; ++i){
	for (int j=0; j<NJets;  ++j){
	TLorentzVector pj=jet[j].lv + muo[i].lv;
	logStream << "     Muon " << i << " and jet " << j << ": Pt " << pj.Pt() << " Px " << pj.Px() << " Py " << pj.Py() << " Pz " << pj.Pz()
		  << " Eta " << pj.Eta() << " Phi " << pj.Phi() << " E " << pj.E() << " M " << pj.M() << endl;	
	}
	}
	if(NMuons>1){
	logStream << "   Muon-Muon combinations -----------------------------------------------------------------------"   << endl;
	for(int i=0; i<NMuons; ++i){
	for(int j=i+1; j<NMuons; ++j){
	TLorentzVector pj=muo[i].lv +muo[j].lv;
	logStream << "     Muon " << i << " and muon " << j << ": Pt " << pj.Pt() << " Px " << pj.Px() << " Py " << pj.Py() << " Pz " << pj.Pz()
		  << " Eta " << pj.Eta() << " Phi " << pj.Phi() << " E " << pj.E() << " M " << pj.M() << endl;	
	}
	}
	}
	}
	if(NMuons>0 && NEles>0){
	logStream << " Muon-Electron combinations --------------------------------------------------------------------"   << endl;
	TLorentzVector pj=ele[0].lv+muo[0].lv;
	logStream << "   Muon " << 0 << " and ele " << 0 << ": Pt " << pj.Pt() << " Px " << pj.Px() << " Py " << pj.Py() << " Pz " << pj.Pz()
		  << " Eta " << pj.Eta() << " Phi " << pj.Phi() << " E " << pj.E() << " M " << pj.M() << endl;	
	}
	logStream << " 2-jet combinations ----------------------------------------------------------------------------" << endl;
	for (int i=0;   i<NJets; ++i){
	for (int j=i+1; j<NJets; ++j){
	TLorentzVector pj = jet[i].lv +jet[j].lv;
	logStream << "  Jet " << i << " and jet " << j << ": Pt " << pj.Pt() << " Px " << pj.Px() << " Py " << pj.Py() << " Pz " << pj.Pz()
		  << " Eta " << pj.Eta() << " Phi " << pj.Phi() << " E " << pj.E() << " M " << pj.M() << endl;	
	}
	}
	logStream << " 3-jet combinations ----------------------------------------------------------------------------" << endl;
	for (int i=0;   i<NJets; ++i){
	for (int j=i+1; j<NJets; ++j){
	for (int k=j+1; k<NJets; ++k){
	TLorentzVector pj = jet[i].lv +jet[j].lv + jet[k].lv;
	logStream << "  Jet " << i << ", jet " << j << " and jet " << k << ": Pt " << pj.Pt() << " Px " << pj.Px() << " Py " << pj.Py() << " Pz " << pj.Pz()
		  << " Eta " << pj.Eta() << " Phi " << pj.Phi() << " E " << pj.E() << " M " << pj.M() << endl;	
	}
	}
	}

	logStream << " MT2 Info ---------------------------------------------------------------------------------------" << endl;
	logStream << "  MT2 " << misc.MT2  << " simple MT2 (sqrt{pt1*pt2*(1+cos phi)}) " << SimpleMT2(true, 1) << " UTM " << hemi[0].UTM.Pt()         << endl;
	logStream << "  MET " << misc.MET  << " HT " << misc.HT << " SqrtSmin " << GetSqrtS(0,true,1,20,2.4,1)          << endl; 
	logStream << "  Hemi-DPhi " << hemi[0].dPhi                                                                     << endl;
	logStream << "  PseudoJet1 "                                                                                    << endl;
	logStream << "   Pt " << hemi[0].lv1.Pt() << " Eta " << hemi[0].lv1.Eta() << " Phi " << hemi[0].lv1.Phi() << " E " << hemi[0].lv1.E() << " M " << hemi[0].lv1.M() << endl; 
	logStream << "   Contributing Objects "                                                                             << endl;
	for (int i=0; i<NJets; ++i){
		if(hemi[0].jindices1[i]>=0) logStream << "    jet " << hemi[0].jindices1[i] << " with Pt " << jet[hemi[0].jindices1[i]].lv.Pt() << endl; 
	}
	for (int i=0; i<NEles; ++i){
		if(hemi[0].eleindices1[i]>=0) logStream << "    ele " << hemi[0].eleindices1[i] << " with Pt " << ele[hemi[0].eleindices1[i]].lv.Pt() << endl; 
	}
	for (int i=0; i<NMuons; ++i){
		if(hemi[0].muoindices1[i]>=0) logStream << "    muo " << hemi[0].muoindices1[i] << " with Pt " << muo[hemi[0].muoindices1[i]].lv.Pt() << endl; 
	}
	logStream << "  PseudoJet2 "                                                                                    << endl;
	logStream << "   Pt " << hemi[0].lv2.Pt() << " Eta " << hemi[0].lv2.Eta() << " Phi " << hemi[0].lv2.Phi() << " E " << hemi[0].lv2.E() << " M " << hemi[0].lv2.M() << endl; 
	logStream << "   Contributing Objects "                                                                             << endl;
	for (int i=0; i<NJets; ++i){
		if(hemi[0].jindices2[i]>=0) logStream << "    jet " << hemi[0].jindices2[i] << " with Pt " << jet[hemi[0].jindices2[i]].lv.Pt() << endl; 
	}
	for (int i=0; i<NEles; ++i){
		if(hemi[0].eleindices2[i]>=0) logStream << "    ele " << hemi[0].eleindices2[i] << " with Pt " << ele[hemi[0].eleindices2[i]].lv.Pt() << endl; 
	}
	for (int i=0; i<NMuons; ++i){
		if(hemi[0].muoindices2[i]>=0) logStream << "    muo " << hemi[0].muoindices2[i] << " with Pt " << muo[hemi[0].muoindices2[i]].lv.Pt() << endl; 
	}
	if(!misc.isData){
	logStream << " GenLevel Info --------------------------------------------------------------------------------"<< endl;
	for(int i=0; i<NGenLepts; ++i){
	logStream << " genParticle ID " << genlept[i].ID  << ", with Mother " << genlept[i].MID << ", and GM " << genlept[i].GMID  
	     << ",  Pt " << genlept[i].lv.Pt() << " Eta " << genlept[i].lv.Eta() << " Phi " << genlept[i].lv.Phi() << " Mass " << genlept[i].lv.M() << endl;	
	}
	logStream << " GenJets       ---------------------------------------------------------------------------------"<< endl;
	for(int i=0; i<NGenJets; ++i){
	logStream << " GenJet Pt " << genjet[i].lv.Pt() << " Eta " << genjet[i].lv.Eta() << " Phi " << genjet[i].lv.Phi() <<  " E " << genjet[i].lv.E() <<" M " << genjet[i].lv.E() << endl;
	}
	}

	if(logfile){
		ofstream f_log ("PrintOut.log", ios::app);
		f_log << logStream.str();
	} else{
		cout << logStream.str();
	}

	return true;
}

// ----------------------------------------------------------------------------------------------------------
ClassImp(MT2Misc)
ClassImp(MT2Znunu)
ClassImp(MT2PileUp)
ClassImp(MT2Trigger)
ClassImp(MT2Jet)
ClassImp(MT2Elec)
ClassImp(MT2Muon)
ClassImp(MT2Hemi)
ClassImp(MT2GenLept)
ClassImp(MT2tree)
