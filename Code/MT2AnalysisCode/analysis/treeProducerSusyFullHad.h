//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 19 14:54:54 2014 by ROOT version 5.32/00
// from TTree treeProducerSusyFullHad/treeProducerSusyFullHad
// found on file: /scratch/mmasciov/treeProducerSusyFullHad_tree_small.root
//////////////////////////////////////////////////////////

#ifndef treeProducerSusyFullHad_h
#define treeProducerSusyFullHad_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class treeProducerSusyFullHad {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           evt;
   Int_t           isData;
   Int_t           HLT_ht350met100;
   Int_t           HLT_HT650;
   Int_t           HLT_MET150;
   Float_t         puWeight;
   Int_t           nTrueInt;
   Float_t         rho;
   Int_t           nVert;
   Int_t           nJet25;
   Int_t           nBJetLoose25;
   Int_t           nBJetMedium25;
   Int_t           nBJetTight25;
   Int_t           nJet40;
   Int_t           nJet40a;
   Int_t           nBJetLoose40;
   Int_t           nBJetMedium40;
   Int_t           nBJetTight40;
   Int_t           nLepGood20;
   Int_t           nLepGood15;
   Int_t           nLepGood10;
   Int_t           GenHeaviestQCDFlavour;
   Float_t         LepEff_1lep;
   Float_t         LepEff_2lep;
   Int_t           GenSusyMScan1;
   Int_t           GenSusyMScan2;
   Int_t           GenSusyMScan3;
   Int_t           GenSusyMScan4;
   Int_t           GenSusyMGluino;
   Int_t           GenSusyMGravitino;
   Int_t           GenSusyMStop;
   Int_t           GenSusyMSbottom;
   Int_t           GenSusyMStop2;
   Int_t           GenSusyMSbottom2;
   Int_t           GenSusyMSquark;
   Int_t           GenSusyMNeutralino;
   Int_t           GenSusyMNeutralino2;
   Int_t           GenSusyMNeutralino3;
   Int_t           GenSusyMNeutralino4;
   Int_t           GenSusyMChargino;
   Int_t           GenSusyMChargino2;
   Int_t           nBJet40;
   Float_t         ht;
   Float_t         deltaPhiMin;
   Float_t         diffMetMht;
   Float_t         mht_pt;
   Float_t         mht_phi;
   Int_t           nMuons10;
   Int_t           nElectrons10;
   Int_t           nTaus20;
   Int_t           nGammas20;
   Float_t         mt2;
   Float_t         mt2_gen;
   Float_t         met_pt;
   Float_t         met_eta;
   Float_t         met_phi;
   Float_t         met_mass;
   Float_t         met_sumEt;
   Float_t         met_genPt;
   Float_t         met_genPhi;
   Float_t         met_genEta;
   Float_t         metNoPU_pt;
   Float_t         metNoPU_eta;
   Float_t         metNoPU_phi;
   Float_t         metNoPU_mass;
   Float_t         pseudoJet1_pt;
   Float_t         pseudoJet1_eta;
   Float_t         pseudoJet1_phi;
   Float_t         pseudoJet1_mass;
   Float_t         pseudoJet2_pt;
   Float_t         pseudoJet2_eta;
   Float_t         pseudoJet2_phi;
   Float_t         pseudoJet2_mass;
   Int_t           nGenQuark;
   Float_t         GenQuark_pt[1];   //[nGenQuark]
   Float_t         GenQuark_eta[1];   //[nGenQuark]
   Float_t         GenQuark_phi[1];   //[nGenQuark]
   Float_t         GenQuark_mass[1];   //[nGenQuark]
   Int_t           GenQuark_pdgId[1];   //[nGenQuark]
   Float_t         GenQuark_charge[1];   //[nGenQuark]
   Int_t           GenQuark_sourceId[1];   //[nGenQuark]
   Int_t           ngenPart;
   Float_t         genPart_pt[127];   //[ngenPart]
   Float_t         genPart_eta[127];   //[ngenPart]
   Float_t         genPart_phi[127];   //[ngenPart]
   Float_t         genPart_mass[127];   //[ngenPart]
   Int_t           genPart_pdgId[127];   //[ngenPart]
   Float_t         genPart_charge[127];   //[ngenPart]
   Int_t           genPart_motherId[127];   //[ngenPart]
   Int_t           genPart_grandmaId[127];   //[ngenPart]
   Int_t           nGenTop;
   Float_t         GenTop_pt[1];   //[nGenTop]
   Float_t         GenTop_eta[1];   //[nGenTop]
   Float_t         GenTop_phi[1];   //[nGenTop]
   Float_t         GenTop_mass[1];   //[nGenTop]
   Int_t           GenTop_pdgId[1];   //[nGenTop]
   Float_t         GenTop_charge[1];   //[nGenTop]
   Int_t           nisoTrack;
   Float_t         isoTrack_pt[5];   //[nisoTrack]
   Float_t         isoTrack_eta[5];   //[nisoTrack]
   Float_t         isoTrack_phi[5];   //[nisoTrack]
   Float_t         isoTrack_mass[5];   //[nisoTrack]
   Int_t           isoTrack_pdgId[5];   //[nisoTrack]
   Int_t           isoTrack_charge[5];   //[nisoTrack]
   Float_t         isoTrack_dz[5];   //[nisoTrack]
   Float_t         isoTrack_absIso[5];   //[nisoTrack]
   Int_t           isoTrack_mcMatchId[5];   //[nisoTrack]
   Int_t           njet;
   Float_t         jet_pt[100];   //[njet]
   Float_t         jet_eta[100];   //[njet]
   Float_t         jet_phi[100];   //[njet]
   Float_t         jet_mass[100];   //[njet]
   Float_t         jet_btagCSV[100];   //[njet]
   Float_t         jet_rawPt[100];   //[njet]
   Float_t         jet_mcPt[100];   //[njet]
   Int_t           jet_mcFlavour[100];   //[njet]
   Int_t           jet_mcMatchId[100];   //[njet]
   Int_t           jet_mcMatchFlav[100];   //[njet]
   Float_t         jet_area[100];   //[njet]
   Int_t           jet_puId[100];   //[njet]
   Int_t           jet_id[100];   //[njet]
   Int_t           nlep;
   Float_t         lep_pt[2];   //[nlep]
   Float_t         lep_eta[2];   //[nlep]
   Float_t         lep_phi[2];   //[nlep]
   Float_t         lep_mass[2];   //[nlep]
   Int_t           lep_pdgId[2];   //[nlep]
   Int_t           lep_charge[2];   //[nlep]
   Float_t         lep_dxy[2];   //[nlep]
   Float_t         lep_dz[2];   //[nlep]
   Float_t         lep_edxy[2];   //[nlep]
   Float_t         lep_edz[2];   //[nlep]
   Float_t         lep_ip3d[2];   //[nlep]
   Float_t         lep_sip3d[2];   //[nlep]
   Int_t           lep_tightId[2];   //[nlep]
   Int_t           lep_convVeto[2];   //[nlep]
   Int_t           lep_lostHits[2];   //[nlep]
   Int_t           lep_looseIdSusy[2];   //[nlep]
   Float_t         lep_relIso03[2];   //[nlep]
   Float_t         lep_relIso04[2];   //[nlep]
   Float_t         lep_chargedHadRelIso03[2];   //[nlep]
   Float_t         lep_chargedHadRelIso04[2];   //[nlep]
   Int_t           lep_convVetoFull[2];   //[nlep]
   Int_t           lep_eleCutIdCSA14_25ns_v1[2];   //[nlep]
   Int_t           lep_eleCutIdCSA14_50ns_v1[2];   //[nlep]
   Int_t           lep_eleMVAId[2];   //[nlep]
   Int_t           lep_tightCharge[2];   //[nlep]
   Float_t         lep_mvaId[2];   //[nlep]
   Float_t         lep_mvaIdTrig[2];   //[nlep]
   Float_t         lep_nStations[2];   //[nlep]
   Float_t         lep_trkKink[2];   //[nlep]
   Int_t           lep_trackerLayers[2];   //[nlep]
   Int_t           lep_pixelLayers[2];   //[nlep]
   Float_t         lep_mvaTTH[2];   //[nlep]
   Float_t         lep_jetPtRatio[2];   //[nlep]
   Float_t         lep_jetBTagCSV[2];   //[nlep]
   Float_t         lep_jetDR[2];   //[nlep]
   Int_t           lep_mcMatchId[2];   //[nlep]
   Int_t           lep_mcMatchAny[2];   //[nlep]
   Int_t           lep_mcMatchTau[2];   //[nlep]
   Int_t           nGenLepFromTau;
   Float_t         GenLepFromTau_pt[3];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_eta[3];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_phi[3];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_mass[3];   //[nGenLepFromTau]
   Int_t           GenLepFromTau_pdgId[3];   //[nGenLepFromTau]
   Float_t         GenLepFromTau_charge[3];   //[nGenLepFromTau]
   Int_t           GenLepFromTau_sourceId[3];   //[nGenLepFromTau]
   Int_t           ngamma;
   Float_t         gamma_pt[1];   //[ngamma]
   Float_t         gamma_eta[1];   //[ngamma]
   Float_t         gamma_phi[1];   //[ngamma]
   Float_t         gamma_mass[1];   //[ngamma]
   Int_t           gamma_pdgId[1];   //[ngamma]
   Int_t           gamma_idCutBased[1];   //[ngamma]
   Float_t         gamma_hOverE[1];   //[ngamma]
   Float_t         gamma_r9[1];   //[ngamma]
   Float_t         gamma_sigmaIetaIeta[1];   //[ngamma]
   Float_t         gamma_chHadIso[1];   //[ngamma]
   Float_t         gamma_neuHadIso[1];   //[ngamma]
   Float_t         gamma_phIso[1];   //[ngamma]
   Int_t           gamma_mcMatchId[1];   //[ngamma]
   Int_t           nGenP6StatusThree;
   Float_t         GenP6StatusThree_pt[1];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_eta[1];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_phi[1];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_mass[1];   //[nGenP6StatusThree]
   Int_t           GenP6StatusThree_pdgId[1];   //[nGenP6StatusThree]
   Float_t         GenP6StatusThree_charge[1];   //[nGenP6StatusThree]
   Int_t           GenP6StatusThree_motherId[1];   //[nGenP6StatusThree]
   Int_t           GenP6StatusThree_grandmaId[1];   //[nGenP6StatusThree]
   Int_t           nGenLep;
   Float_t         GenLep_pt[1];   //[nGenLep]
   Float_t         GenLep_eta[1];   //[nGenLep]
   Float_t         GenLep_phi[1];   //[nGenLep]
   Float_t         GenLep_mass[1];   //[nGenLep]
   Int_t           GenLep_pdgId[1];   //[nGenLep]
   Float_t         GenLep_charge[1];   //[nGenLep]
   Int_t           GenLep_sourceId[1];   //[nGenLep]
   Int_t           ntau;
   Float_t         tau_pt[2];   //[ntau]
   Float_t         tau_eta[2];   //[ntau]
   Float_t         tau_phi[2];   //[ntau]
   Float_t         tau_mass[2];   //[ntau]
   Int_t           tau_pdgId[2];   //[ntau]
   Int_t           tau_charge[2];   //[ntau]
   Float_t         tau_dxy[2];   //[ntau]
   Float_t         tau_dz[2];   //[ntau]
   Int_t           tau_idCI3hit[2];   //[ntau]
   Float_t         tau_isoCI3hit[2];   //[ntau]
   Int_t           tau_mcMatchId[2];   //[ntau]
   Int_t           nGenBQuark;
   Float_t         GenBQuark_pt[1];   //[nGenBQuark]
   Float_t         GenBQuark_eta[1];   //[nGenBQuark]
   Float_t         GenBQuark_phi[1];   //[nGenBQuark]
   Float_t         GenBQuark_mass[1];   //[nGenBQuark]
   Int_t           GenBQuark_pdgId[1];   //[nGenBQuark]
   Float_t         GenBQuark_charge[1];   //[nGenBQuark]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_HLT_ht350met100;   //!
   TBranch        *b_HLT_HT650;   //!
   TBranch        *b_HLT_MET150;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_nTrueInt;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_nVert;   //!
   TBranch        *b_nJet25;   //!
   TBranch        *b_nBJetLoose25;   //!
   TBranch        *b_nBJetMedium25;   //!
   TBranch        *b_nBJetTight25;   //!
   TBranch        *b_nJet40;   //!
   TBranch        *b_nJet40a;   //!
   TBranch        *b_nBJetLoose40;   //!
   TBranch        *b_nBJetMedium40;   //!
   TBranch        *b_nBJetTight40;   //!
   TBranch        *b_nLepGood20;   //!
   TBranch        *b_nLepGood15;   //!
   TBranch        *b_nLepGood10;   //!
   TBranch        *b_GenHeaviestQCDFlavour;   //!
   TBranch        *b_LepEff_1lep;   //!
   TBranch        *b_LepEff_2lep;   //!
   TBranch        *b_GenSusyMScan1;   //!
   TBranch        *b_GenSusyMScan2;   //!
   TBranch        *b_GenSusyMScan3;   //!
   TBranch        *b_GenSusyMScan4;   //!
   TBranch        *b_GenSusyMGluino;   //!
   TBranch        *b_GenSusyMGravitino;   //!
   TBranch        *b_GenSusyMStop;   //!
   TBranch        *b_GenSusyMSbottom;   //!
   TBranch        *b_GenSusyMStop2;   //!
   TBranch        *b_GenSusyMSbottom2;   //!
   TBranch        *b_GenSusyMSquark;   //!
   TBranch        *b_GenSusyMNeutralino;   //!
   TBranch        *b_GenSusyMNeutralino2;   //!
   TBranch        *b_GenSusyMNeutralino3;   //!
   TBranch        *b_GenSusyMNeutralino4;   //!
   TBranch        *b_GenSusyMChargino;   //!
   TBranch        *b_GenSusyMChargino2;   //!
   TBranch        *b_nBJet40;   //!
   TBranch        *b_ht;   //!
   TBranch        *b_deltaPhiMin;   //!
   TBranch        *b_diffMetMht;   //!
   TBranch        *b_mht_pt;   //!
   TBranch        *b_mht_phi;   //!
   TBranch        *b_nMuons10;   //!
   TBranch        *b_nElectrons10;   //!
   TBranch        *b_nTaus20;   //!
   TBranch        *b_nGammas20;   //!
   TBranch        *b_mt2;   //!
   TBranch        *b_mt2_gen;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_eta;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_mass;   //!
   TBranch        *b_met_sumEt;   //!
   TBranch        *b_met_genPt;   //!
   TBranch        *b_met_genPhi;   //!
   TBranch        *b_met_genEta;   //!
   TBranch        *b_metNoPU_pt;   //!
   TBranch        *b_metNoPU_eta;   //!
   TBranch        *b_metNoPU_phi;   //!
   TBranch        *b_metNoPU_mass;   //!
   TBranch        *b_pseudoJet1_pt;   //!
   TBranch        *b_pseudoJet1_eta;   //!
   TBranch        *b_pseudoJet1_phi;   //!
   TBranch        *b_pseudoJet1_mass;   //!
   TBranch        *b_pseudoJet2_pt;   //!
   TBranch        *b_pseudoJet2_eta;   //!
   TBranch        *b_pseudoJet2_phi;   //!
   TBranch        *b_pseudoJet2_mass;   //!
   TBranch        *b_nGenQuark;   //!
   TBranch        *b_GenQuark_pt;   //!
   TBranch        *b_GenQuark_eta;   //!
   TBranch        *b_GenQuark_phi;   //!
   TBranch        *b_GenQuark_mass;   //!
   TBranch        *b_GenQuark_pdgId;   //!
   TBranch        *b_GenQuark_charge;   //!
   TBranch        *b_GenQuark_sourceId;   //!
   TBranch        *b_ngenPart;   //!
   TBranch        *b_genPart_pt;   //!
   TBranch        *b_genPart_eta;   //!
   TBranch        *b_genPart_phi;   //!
   TBranch        *b_genPart_mass;   //!
   TBranch        *b_genPart_pdgId;   //!
   TBranch        *b_genPart_charge;   //!
   TBranch        *b_genPart_motherId;   //!
   TBranch        *b_genPart_grandmaId;   //!
   TBranch        *b_nGenTop;   //!
   TBranch        *b_GenTop_pt;   //!
   TBranch        *b_GenTop_eta;   //!
   TBranch        *b_GenTop_phi;   //!
   TBranch        *b_GenTop_mass;   //!
   TBranch        *b_GenTop_pdgId;   //!
   TBranch        *b_GenTop_charge;   //!
   TBranch        *b_nisoTrack;   //!
   TBranch        *b_isoTrack_pt;   //!
   TBranch        *b_isoTrack_eta;   //!
   TBranch        *b_isoTrack_phi;   //!
   TBranch        *b_isoTrack_mass;   //!
   TBranch        *b_isoTrack_pdgId;   //!
   TBranch        *b_isoTrack_charge;   //!
   TBranch        *b_isoTrack_dz;   //!
   TBranch        *b_isoTrack_absIso;   //!
   TBranch        *b_isoTrack_mcMatchId;   //!
   TBranch        *b_njet;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_mass;   //!
   TBranch        *b_jet_btagCSV;   //!
   TBranch        *b_jet_rawPt;   //!
   TBranch        *b_jet_mcPt;   //!
   TBranch        *b_jet_mcFlavour;   //!
   TBranch        *b_jet_mcMatchId;   //!
   TBranch        *b_jet_mcMatchFlav;   //!
   TBranch        *b_jet_area;   //!
   TBranch        *b_jet_puId;   //!
   TBranch        *b_jet_id;   //!
   TBranch        *b_nlep;   //!
   TBranch        *b_lep_pt;   //!
   TBranch        *b_lep_eta;   //!
   TBranch        *b_lep_phi;   //!
   TBranch        *b_lep_mass;   //!
   TBranch        *b_lep_pdgId;   //!
   TBranch        *b_lep_charge;   //!
   TBranch        *b_lep_dxy;   //!
   TBranch        *b_lep_dz;   //!
   TBranch        *b_lep_edxy;   //!
   TBranch        *b_lep_edz;   //!
   TBranch        *b_lep_ip3d;   //!
   TBranch        *b_lep_sip3d;   //!
   TBranch        *b_lep_tightId;   //!
   TBranch        *b_lep_convVeto;   //!
   TBranch        *b_lep_lostHits;   //!
   TBranch        *b_lep_looseIdSusy;   //!
   TBranch        *b_lep_relIso03;   //!
   TBranch        *b_lep_relIso04;   //!
   TBranch        *b_lep_chargedHadRelIso03;   //!
   TBranch        *b_lep_chargedHadRelIso04;   //!
   TBranch        *b_lep_convVetoFull;   //!
   TBranch        *b_lep_eleCutIdCSA14_25ns_v1;   //!
   TBranch        *b_lep_eleCutIdCSA14_50ns_v1;   //!
   TBranch        *b_lep_eleMVAId;   //!
   TBranch        *b_lep_tightCharge;   //!
   TBranch        *b_lep_mvaId;   //!
   TBranch        *b_lep_mvaIdTrig;   //!
   TBranch        *b_lep_nStations;   //!
   TBranch        *b_lep_trkKink;   //!
   TBranch        *b_lep_trackerLayers;   //!
   TBranch        *b_lep_pixelLayers;   //!
   TBranch        *b_lep_mvaTTH;   //!
   TBranch        *b_lep_jetPtRatio;   //!
   TBranch        *b_lep_jetBTagCSV;   //!
   TBranch        *b_lep_jetDR;   //!
   TBranch        *b_lep_mcMatchId;   //!
   TBranch        *b_lep_mcMatchAny;   //!
   TBranch        *b_lep_mcMatchTau;   //!
   TBranch        *b_nGenLepFromTau;   //!
   TBranch        *b_GenLepFromTau_pt;   //!
   TBranch        *b_GenLepFromTau_eta;   //!
   TBranch        *b_GenLepFromTau_phi;   //!
   TBranch        *b_GenLepFromTau_mass;   //!
   TBranch        *b_GenLepFromTau_pdgId;   //!
   TBranch        *b_GenLepFromTau_charge;   //!
   TBranch        *b_GenLepFromTau_sourceId;   //!
   TBranch        *b_ngamma;   //!
   TBranch        *b_gamma_pt;   //!
   TBranch        *b_gamma_eta;   //!
   TBranch        *b_gamma_phi;   //!
   TBranch        *b_gamma_mass;   //!
   TBranch        *b_gamma_pdgId;   //!
   TBranch        *b_gamma_idCutBased;   //!
   TBranch        *b_gamma_hOverE;   //!
   TBranch        *b_gamma_r9;   //!
   TBranch        *b_gamma_sigmaIetaIeta;   //!
   TBranch        *b_gamma_chHadIso;   //!
   TBranch        *b_gamma_neuHadIso;   //!
   TBranch        *b_gamma_phIso;   //!
   TBranch        *b_gamma_mcMatchId;   //!
   TBranch        *b_nGenP6StatusThree;   //!
   TBranch        *b_GenP6StatusThree_pt;   //!
   TBranch        *b_GenP6StatusThree_eta;   //!
   TBranch        *b_GenP6StatusThree_phi;   //!
   TBranch        *b_GenP6StatusThree_mass;   //!
   TBranch        *b_GenP6StatusThree_pdgId;   //!
   TBranch        *b_GenP6StatusThree_charge;   //!
   TBranch        *b_GenP6StatusThree_motherId;   //!
   TBranch        *b_GenP6StatusThree_grandmaId;   //!
   TBranch        *b_nGenLep;   //!
   TBranch        *b_GenLep_pt;   //!
   TBranch        *b_GenLep_eta;   //!
   TBranch        *b_GenLep_phi;   //!
   TBranch        *b_GenLep_mass;   //!
   TBranch        *b_GenLep_pdgId;   //!
   TBranch        *b_GenLep_charge;   //!
   TBranch        *b_GenLep_sourceId;   //!
   TBranch        *b_ntau;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_eta;   //!
   TBranch        *b_tau_phi;   //!
   TBranch        *b_tau_mass;   //!
   TBranch        *b_tau_pdgId;   //!
   TBranch        *b_tau_charge;   //!
   TBranch        *b_tau_dxy;   //!
   TBranch        *b_tau_dz;   //!
   TBranch        *b_tau_idCI3hit;   //!
   TBranch        *b_tau_isoCI3hit;   //!
   TBranch        *b_tau_mcMatchId;   //!
   TBranch        *b_nGenBQuark;   //!
   TBranch        *b_GenBQuark_pt;   //!
   TBranch        *b_GenBQuark_eta;   //!
   TBranch        *b_GenBQuark_phi;   //!
   TBranch        *b_GenBQuark_mass;   //!
   TBranch        *b_GenBQuark_pdgId;   //!
   TBranch        *b_GenBQuark_charge;   //!

   treeProducerSusyFullHad(TTree *tree=0);
   virtual ~treeProducerSusyFullHad();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef treeProducerSusyFullHad_cxx
treeProducerSusyFullHad::treeProducerSusyFullHad(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/scratch/mmasciov/treeProducerSusyFullHad_tree_small.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("/scratch/mmasciov/treeProducerSusyFullHad_tree_small.root");
//      }
//      f->GetObject("treeProducerSusyFullHad",tree);
//
//   }
//   Init(tree);
}

treeProducerSusyFullHad::~treeProducerSusyFullHad()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t treeProducerSusyFullHad::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t treeProducerSusyFullHad::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void treeProducerSusyFullHad::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("HLT_ht350met100", &HLT_ht350met100, &b_HLT_ht350met100);
   fChain->SetBranchAddress("HLT_HT650", &HLT_HT650, &b_HLT_HT650);
   fChain->SetBranchAddress("HLT_MET150", &HLT_MET150, &b_HLT_MET150);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("nTrueInt", &nTrueInt, &b_nTrueInt);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nVert", &nVert, &b_nVert);
   fChain->SetBranchAddress("nJet25", &nJet25, &b_nJet25);
   fChain->SetBranchAddress("nBJetLoose25", &nBJetLoose25, &b_nBJetLoose25);
   fChain->SetBranchAddress("nBJetMedium25", &nBJetMedium25, &b_nBJetMedium25);
   fChain->SetBranchAddress("nBJetTight25", &nBJetTight25, &b_nBJetTight25);
   fChain->SetBranchAddress("nJet40", &nJet40, &b_nJet40);
   fChain->SetBranchAddress("nJet40a", &nJet40a, &b_nJet40a);
   fChain->SetBranchAddress("nBJetLoose40", &nBJetLoose40, &b_nBJetLoose40);
   fChain->SetBranchAddress("nBJetMedium40", &nBJetMedium40, &b_nBJetMedium40);
   fChain->SetBranchAddress("nBJetTight40", &nBJetTight40, &b_nBJetTight40);
   fChain->SetBranchAddress("nLepGood20", &nLepGood20, &b_nLepGood20);
   fChain->SetBranchAddress("nLepGood15", &nLepGood15, &b_nLepGood15);
   fChain->SetBranchAddress("nLepGood10", &nLepGood10, &b_nLepGood10);
   fChain->SetBranchAddress("GenHeaviestQCDFlavour", &GenHeaviestQCDFlavour, &b_GenHeaviestQCDFlavour);
   fChain->SetBranchAddress("LepEff_1lep", &LepEff_1lep, &b_LepEff_1lep);
   fChain->SetBranchAddress("LepEff_2lep", &LepEff_2lep, &b_LepEff_2lep);
   fChain->SetBranchAddress("GenSusyMScan1", &GenSusyMScan1, &b_GenSusyMScan1);
   fChain->SetBranchAddress("GenSusyMScan2", &GenSusyMScan2, &b_GenSusyMScan2);
   fChain->SetBranchAddress("GenSusyMScan3", &GenSusyMScan3, &b_GenSusyMScan3);
   fChain->SetBranchAddress("GenSusyMScan4", &GenSusyMScan4, &b_GenSusyMScan4);
   fChain->SetBranchAddress("GenSusyMGluino", &GenSusyMGluino, &b_GenSusyMGluino);
   fChain->SetBranchAddress("GenSusyMGravitino", &GenSusyMGravitino, &b_GenSusyMGravitino);
   fChain->SetBranchAddress("GenSusyMStop", &GenSusyMStop, &b_GenSusyMStop);
   fChain->SetBranchAddress("GenSusyMSbottom", &GenSusyMSbottom, &b_GenSusyMSbottom);
   fChain->SetBranchAddress("GenSusyMStop2", &GenSusyMStop2, &b_GenSusyMStop2);
   fChain->SetBranchAddress("GenSusyMSbottom2", &GenSusyMSbottom2, &b_GenSusyMSbottom2);
   fChain->SetBranchAddress("GenSusyMSquark", &GenSusyMSquark, &b_GenSusyMSquark);
   fChain->SetBranchAddress("GenSusyMNeutralino", &GenSusyMNeutralino, &b_GenSusyMNeutralino);
   fChain->SetBranchAddress("GenSusyMNeutralino2", &GenSusyMNeutralino2, &b_GenSusyMNeutralino2);
   fChain->SetBranchAddress("GenSusyMNeutralino3", &GenSusyMNeutralino3, &b_GenSusyMNeutralino3);
   fChain->SetBranchAddress("GenSusyMNeutralino4", &GenSusyMNeutralino4, &b_GenSusyMNeutralino4);
   fChain->SetBranchAddress("GenSusyMChargino", &GenSusyMChargino, &b_GenSusyMChargino);
   fChain->SetBranchAddress("GenSusyMChargino2", &GenSusyMChargino2, &b_GenSusyMChargino2);
   fChain->SetBranchAddress("nBJet40", &nBJet40, &b_nBJet40);
   fChain->SetBranchAddress("ht", &ht, &b_ht);
   fChain->SetBranchAddress("deltaPhiMin", &deltaPhiMin, &b_deltaPhiMin);
   fChain->SetBranchAddress("diffMetMht", &diffMetMht, &b_diffMetMht);
   fChain->SetBranchAddress("mht_pt", &mht_pt, &b_mht_pt);
   fChain->SetBranchAddress("mht_phi", &mht_phi, &b_mht_phi);
   fChain->SetBranchAddress("nMuons10", &nMuons10, &b_nMuons10);
   fChain->SetBranchAddress("nElectrons10", &nElectrons10, &b_nElectrons10);
   fChain->SetBranchAddress("nTaus20", &nTaus20, &b_nTaus20);
   fChain->SetBranchAddress("nGammas20", &nGammas20, &b_nGammas20);
   fChain->SetBranchAddress("mt2", &mt2, &b_mt2);
   fChain->SetBranchAddress("mt2_gen", &mt2_gen, &b_mt2_gen);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_mass", &met_mass, &b_met_mass);
   fChain->SetBranchAddress("met_sumEt", &met_sumEt, &b_met_sumEt);
   fChain->SetBranchAddress("met_genPt", &met_genPt, &b_met_genPt);
   fChain->SetBranchAddress("met_genPhi", &met_genPhi, &b_met_genPhi);
   fChain->SetBranchAddress("met_genEta", &met_genEta, &b_met_genEta);
   fChain->SetBranchAddress("metNoPU_pt", &metNoPU_pt, &b_metNoPU_pt);
   fChain->SetBranchAddress("metNoPU_eta", &metNoPU_eta, &b_metNoPU_eta);
   fChain->SetBranchAddress("metNoPU_phi", &metNoPU_phi, &b_metNoPU_phi);
   fChain->SetBranchAddress("metNoPU_mass", &metNoPU_mass, &b_metNoPU_mass);
   fChain->SetBranchAddress("pseudoJet1_pt", &pseudoJet1_pt, &b_pseudoJet1_pt);
   fChain->SetBranchAddress("pseudoJet1_eta", &pseudoJet1_eta, &b_pseudoJet1_eta);
   fChain->SetBranchAddress("pseudoJet1_phi", &pseudoJet1_phi, &b_pseudoJet1_phi);
   fChain->SetBranchAddress("pseudoJet1_mass", &pseudoJet1_mass, &b_pseudoJet1_mass);
   fChain->SetBranchAddress("pseudoJet2_pt", &pseudoJet2_pt, &b_pseudoJet2_pt);
   fChain->SetBranchAddress("pseudoJet2_eta", &pseudoJet2_eta, &b_pseudoJet2_eta);
   fChain->SetBranchAddress("pseudoJet2_phi", &pseudoJet2_phi, &b_pseudoJet2_phi);
   fChain->SetBranchAddress("pseudoJet2_mass", &pseudoJet2_mass, &b_pseudoJet2_mass);
   fChain->SetBranchAddress("nGenQuark", &nGenQuark, &b_nGenQuark);
   fChain->SetBranchAddress("GenQuark_pt", &GenQuark_pt, &b_GenQuark_pt);
   fChain->SetBranchAddress("GenQuark_eta", &GenQuark_eta, &b_GenQuark_eta);
   fChain->SetBranchAddress("GenQuark_phi", &GenQuark_phi, &b_GenQuark_phi);
   fChain->SetBranchAddress("GenQuark_mass", &GenQuark_mass, &b_GenQuark_mass);
   fChain->SetBranchAddress("GenQuark_pdgId", &GenQuark_pdgId, &b_GenQuark_pdgId);
   fChain->SetBranchAddress("GenQuark_charge", &GenQuark_charge, &b_GenQuark_charge);
   fChain->SetBranchAddress("GenQuark_sourceId", &GenQuark_sourceId, &b_GenQuark_sourceId);
   fChain->SetBranchAddress("ngenPart", &ngenPart, &b_ngenPart);
   fChain->SetBranchAddress("genPart_pt", genPart_pt, &b_genPart_pt);
   fChain->SetBranchAddress("genPart_eta", genPart_eta, &b_genPart_eta);
   fChain->SetBranchAddress("genPart_phi", genPart_phi, &b_genPart_phi);
   fChain->SetBranchAddress("genPart_mass", genPart_mass, &b_genPart_mass);
   fChain->SetBranchAddress("genPart_pdgId", genPart_pdgId, &b_genPart_pdgId);
   fChain->SetBranchAddress("genPart_charge", genPart_charge, &b_genPart_charge);
   fChain->SetBranchAddress("genPart_motherId", genPart_motherId, &b_genPart_motherId);
   fChain->SetBranchAddress("genPart_grandmaId", genPart_grandmaId, &b_genPart_grandmaId);
   fChain->SetBranchAddress("nGenTop", &nGenTop, &b_nGenTop);
   fChain->SetBranchAddress("GenTop_pt", &GenTop_pt, &b_GenTop_pt);
   fChain->SetBranchAddress("GenTop_eta", &GenTop_eta, &b_GenTop_eta);
   fChain->SetBranchAddress("GenTop_phi", &GenTop_phi, &b_GenTop_phi);
   fChain->SetBranchAddress("GenTop_mass", &GenTop_mass, &b_GenTop_mass);
   fChain->SetBranchAddress("GenTop_pdgId", &GenTop_pdgId, &b_GenTop_pdgId);
   fChain->SetBranchAddress("GenTop_charge", &GenTop_charge, &b_GenTop_charge);
   fChain->SetBranchAddress("nisoTrack", &nisoTrack, &b_nisoTrack);
   fChain->SetBranchAddress("isoTrack_pt", isoTrack_pt, &b_isoTrack_pt);
   fChain->SetBranchAddress("isoTrack_eta", isoTrack_eta, &b_isoTrack_eta);
   fChain->SetBranchAddress("isoTrack_phi", isoTrack_phi, &b_isoTrack_phi);
   fChain->SetBranchAddress("isoTrack_mass", isoTrack_mass, &b_isoTrack_mass);
   fChain->SetBranchAddress("isoTrack_pdgId", isoTrack_pdgId, &b_isoTrack_pdgId);
   fChain->SetBranchAddress("isoTrack_charge", isoTrack_charge, &b_isoTrack_charge);
   fChain->SetBranchAddress("isoTrack_dz", isoTrack_dz, &b_isoTrack_dz);
   fChain->SetBranchAddress("isoTrack_absIso", isoTrack_absIso, &b_isoTrack_absIso);
   fChain->SetBranchAddress("isoTrack_mcMatchId", isoTrack_mcMatchId, &b_isoTrack_mcMatchId);
   fChain->SetBranchAddress("njet", &njet, &b_njet);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_mass", jet_mass, &b_jet_mass);
   fChain->SetBranchAddress("jet_btagCSV", jet_btagCSV, &b_jet_btagCSV);
   fChain->SetBranchAddress("jet_rawPt", jet_rawPt, &b_jet_rawPt);
   fChain->SetBranchAddress("jet_mcPt", jet_mcPt, &b_jet_mcPt);
   fChain->SetBranchAddress("jet_mcFlavour", jet_mcFlavour, &b_jet_mcFlavour);
   fChain->SetBranchAddress("jet_mcMatchId", jet_mcMatchId, &b_jet_mcMatchId);
   fChain->SetBranchAddress("jet_mcMatchFlav", jet_mcMatchFlav, &b_jet_mcMatchFlav);
   fChain->SetBranchAddress("jet_area", jet_area, &b_jet_area);
   fChain->SetBranchAddress("jet_puId", jet_puId, &b_jet_puId);
   fChain->SetBranchAddress("jet_id", jet_id, &b_jet_id);
   fChain->SetBranchAddress("nlep", &nlep, &b_nlep);
   fChain->SetBranchAddress("lep_pt", lep_pt, &b_lep_pt);
   fChain->SetBranchAddress("lep_eta", lep_eta, &b_lep_eta);
   fChain->SetBranchAddress("lep_phi", lep_phi, &b_lep_phi);
   fChain->SetBranchAddress("lep_mass", lep_mass, &b_lep_mass);
   fChain->SetBranchAddress("lep_pdgId", lep_pdgId, &b_lep_pdgId);
   fChain->SetBranchAddress("lep_charge", lep_charge, &b_lep_charge);
   fChain->SetBranchAddress("lep_dxy", lep_dxy, &b_lep_dxy);
   fChain->SetBranchAddress("lep_dz", lep_dz, &b_lep_dz);
   fChain->SetBranchAddress("lep_edxy", lep_edxy, &b_lep_edxy);
   fChain->SetBranchAddress("lep_edz", lep_edz, &b_lep_edz);
   fChain->SetBranchAddress("lep_ip3d", lep_ip3d, &b_lep_ip3d);
   fChain->SetBranchAddress("lep_sip3d", lep_sip3d, &b_lep_sip3d);
   fChain->SetBranchAddress("lep_tightId", lep_tightId, &b_lep_tightId);
   fChain->SetBranchAddress("lep_convVeto", lep_convVeto, &b_lep_convVeto);
   fChain->SetBranchAddress("lep_lostHits", lep_lostHits, &b_lep_lostHits);
   fChain->SetBranchAddress("lep_looseIdSusy", lep_looseIdSusy, &b_lep_looseIdSusy);
   fChain->SetBranchAddress("lep_relIso03", lep_relIso03, &b_lep_relIso03);
   fChain->SetBranchAddress("lep_relIso04", lep_relIso04, &b_lep_relIso04);
   fChain->SetBranchAddress("lep_chargedHadRelIso03", lep_chargedHadRelIso03, &b_lep_chargedHadRelIso03);
   fChain->SetBranchAddress("lep_chargedHadRelIso04", lep_chargedHadRelIso04, &b_lep_chargedHadRelIso04);
   fChain->SetBranchAddress("lep_convVetoFull", lep_convVetoFull, &b_lep_convVetoFull);
   fChain->SetBranchAddress("lep_eleCutIdCSA14_25ns_v1", lep_eleCutIdCSA14_25ns_v1, &b_lep_eleCutIdCSA14_25ns_v1);
   fChain->SetBranchAddress("lep_eleCutIdCSA14_50ns_v1", lep_eleCutIdCSA14_50ns_v1, &b_lep_eleCutIdCSA14_50ns_v1);
   fChain->SetBranchAddress("lep_eleMVAId", lep_eleMVAId, &b_lep_eleMVAId);
   fChain->SetBranchAddress("lep_tightCharge", lep_tightCharge, &b_lep_tightCharge);
   fChain->SetBranchAddress("lep_mvaId", lep_mvaId, &b_lep_mvaId);
   fChain->SetBranchAddress("lep_mvaIdTrig", lep_mvaIdTrig, &b_lep_mvaIdTrig);
   fChain->SetBranchAddress("lep_nStations", lep_nStations, &b_lep_nStations);
   fChain->SetBranchAddress("lep_trkKink", lep_trkKink, &b_lep_trkKink);
   fChain->SetBranchAddress("lep_trackerLayers", lep_trackerLayers, &b_lep_trackerLayers);
   fChain->SetBranchAddress("lep_pixelLayers", lep_pixelLayers, &b_lep_pixelLayers);
   fChain->SetBranchAddress("lep_mvaTTH", lep_mvaTTH, &b_lep_mvaTTH);
   fChain->SetBranchAddress("lep_jetPtRatio", lep_jetPtRatio, &b_lep_jetPtRatio);
   fChain->SetBranchAddress("lep_jetBTagCSV", lep_jetBTagCSV, &b_lep_jetBTagCSV);
   fChain->SetBranchAddress("lep_jetDR", lep_jetDR, &b_lep_jetDR);
   fChain->SetBranchAddress("lep_mcMatchId", lep_mcMatchId, &b_lep_mcMatchId);
   fChain->SetBranchAddress("lep_mcMatchAny", lep_mcMatchAny, &b_lep_mcMatchAny);
   fChain->SetBranchAddress("lep_mcMatchTau", lep_mcMatchTau, &b_lep_mcMatchTau);
   fChain->SetBranchAddress("nGenLepFromTau", &nGenLepFromTau, &b_nGenLepFromTau);
   fChain->SetBranchAddress("GenLepFromTau_pt", GenLepFromTau_pt, &b_GenLepFromTau_pt);
   fChain->SetBranchAddress("GenLepFromTau_eta", GenLepFromTau_eta, &b_GenLepFromTau_eta);
   fChain->SetBranchAddress("GenLepFromTau_phi", GenLepFromTau_phi, &b_GenLepFromTau_phi);
   fChain->SetBranchAddress("GenLepFromTau_mass", GenLepFromTau_mass, &b_GenLepFromTau_mass);
   fChain->SetBranchAddress("GenLepFromTau_pdgId", GenLepFromTau_pdgId, &b_GenLepFromTau_pdgId);
   fChain->SetBranchAddress("GenLepFromTau_charge", GenLepFromTau_charge, &b_GenLepFromTau_charge);
   fChain->SetBranchAddress("GenLepFromTau_sourceId", GenLepFromTau_sourceId, &b_GenLepFromTau_sourceId);
   fChain->SetBranchAddress("ngamma", &ngamma, &b_ngamma);
   fChain->SetBranchAddress("gamma_pt", gamma_pt, &b_gamma_pt);
   fChain->SetBranchAddress("gamma_eta", gamma_eta, &b_gamma_eta);
   fChain->SetBranchAddress("gamma_phi", gamma_phi, &b_gamma_phi);
   fChain->SetBranchAddress("gamma_mass", gamma_mass, &b_gamma_mass);
   fChain->SetBranchAddress("gamma_pdgId", gamma_pdgId, &b_gamma_pdgId);
   fChain->SetBranchAddress("gamma_idCutBased", gamma_idCutBased, &b_gamma_idCutBased);
   fChain->SetBranchAddress("gamma_hOverE", gamma_hOverE, &b_gamma_hOverE);
   fChain->SetBranchAddress("gamma_r9", gamma_r9, &b_gamma_r9);
   fChain->SetBranchAddress("gamma_sigmaIetaIeta", gamma_sigmaIetaIeta, &b_gamma_sigmaIetaIeta);
   fChain->SetBranchAddress("gamma_chHadIso", gamma_chHadIso, &b_gamma_chHadIso);
   fChain->SetBranchAddress("gamma_neuHadIso", gamma_neuHadIso, &b_gamma_neuHadIso);
   fChain->SetBranchAddress("gamma_phIso", gamma_phIso, &b_gamma_phIso);
   fChain->SetBranchAddress("gamma_mcMatchId", gamma_mcMatchId, &b_gamma_mcMatchId);
   fChain->SetBranchAddress("nGenP6StatusThree", &nGenP6StatusThree, &b_nGenP6StatusThree);
   fChain->SetBranchAddress("GenP6StatusThree_pt", &GenP6StatusThree_pt, &b_GenP6StatusThree_pt);
   fChain->SetBranchAddress("GenP6StatusThree_eta", &GenP6StatusThree_eta, &b_GenP6StatusThree_eta);
   fChain->SetBranchAddress("GenP6StatusThree_phi", &GenP6StatusThree_phi, &b_GenP6StatusThree_phi);
   fChain->SetBranchAddress("GenP6StatusThree_mass", &GenP6StatusThree_mass, &b_GenP6StatusThree_mass);
   fChain->SetBranchAddress("GenP6StatusThree_pdgId", &GenP6StatusThree_pdgId, &b_GenP6StatusThree_pdgId);
   fChain->SetBranchAddress("GenP6StatusThree_charge", &GenP6StatusThree_charge, &b_GenP6StatusThree_charge);
   fChain->SetBranchAddress("GenP6StatusThree_motherId", &GenP6StatusThree_motherId, &b_GenP6StatusThree_motherId);
   fChain->SetBranchAddress("GenP6StatusThree_grandmaId", &GenP6StatusThree_grandmaId, &b_GenP6StatusThree_grandmaId);
   fChain->SetBranchAddress("nGenLep", &nGenLep, &b_nGenLep);
   fChain->SetBranchAddress("GenLep_pt", GenLep_pt, &b_GenLep_pt);
   fChain->SetBranchAddress("GenLep_eta", GenLep_eta, &b_GenLep_eta);
   fChain->SetBranchAddress("GenLep_phi", GenLep_phi, &b_GenLep_phi);
   fChain->SetBranchAddress("GenLep_mass", GenLep_mass, &b_GenLep_mass);
   fChain->SetBranchAddress("GenLep_pdgId", GenLep_pdgId, &b_GenLep_pdgId);
   fChain->SetBranchAddress("GenLep_charge", GenLep_charge, &b_GenLep_charge);
   fChain->SetBranchAddress("GenLep_sourceId", GenLep_sourceId, &b_GenLep_sourceId);
   fChain->SetBranchAddress("ntau", &ntau, &b_ntau);
   fChain->SetBranchAddress("tau_pt", tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("tau_eta", tau_eta, &b_tau_eta);
   fChain->SetBranchAddress("tau_phi", tau_phi, &b_tau_phi);
   fChain->SetBranchAddress("tau_mass", tau_mass, &b_tau_mass);
   fChain->SetBranchAddress("tau_pdgId", tau_pdgId, &b_tau_pdgId);
   fChain->SetBranchAddress("tau_charge", tau_charge, &b_tau_charge);
   fChain->SetBranchAddress("tau_dxy", tau_dxy, &b_tau_dxy);
   fChain->SetBranchAddress("tau_dz", tau_dz, &b_tau_dz);
   fChain->SetBranchAddress("tau_idCI3hit", tau_idCI3hit, &b_tau_idCI3hit);
   fChain->SetBranchAddress("tau_isoCI3hit", tau_isoCI3hit, &b_tau_isoCI3hit);
   fChain->SetBranchAddress("tau_mcMatchId", tau_mcMatchId, &b_tau_mcMatchId);
   fChain->SetBranchAddress("nGenBQuark", &nGenBQuark, &b_nGenBQuark);
   fChain->SetBranchAddress("GenBQuark_pt", &GenBQuark_pt, &b_GenBQuark_pt);
   fChain->SetBranchAddress("GenBQuark_eta", &GenBQuark_eta, &b_GenBQuark_eta);
   fChain->SetBranchAddress("GenBQuark_phi", &GenBQuark_phi, &b_GenBQuark_phi);
   fChain->SetBranchAddress("GenBQuark_mass", &GenBQuark_mass, &b_GenBQuark_mass);
   fChain->SetBranchAddress("GenBQuark_pdgId", &GenBQuark_pdgId, &b_GenBQuark_pdgId);
   fChain->SetBranchAddress("GenBQuark_charge", &GenBQuark_charge, &b_GenBQuark_charge);
   Notify();
}

Bool_t treeProducerSusyFullHad::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void treeProducerSusyFullHad::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t treeProducerSusyFullHad::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef treeProducerSusyFullHad_cxx
