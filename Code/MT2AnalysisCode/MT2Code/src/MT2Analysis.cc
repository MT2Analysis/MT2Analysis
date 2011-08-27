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

	fRequiredHLT.clear();
	fVetoedHLT.clear();
	
	fBEfiles.clear();
	fTPfiles.clear();
}

MT2Analysis::~MT2Analysis(){
	delete fJecUncPF;    
	delete fJecUncCalo;  
}

void MT2Analysis::End(){
	cout << " *************************************************************** " << endl;
	cout << " MT2Analysis::End()                                             " << endl;
	cout << " *************************************************************** " << endl;
	
	fHistFile->cd();	

	// write tree
	fATree->Write();
	fHistFile                ->Close();

	cout << " MT2Analysis::RealEnd()                                             " << endl;
	cout << " *************************************************************** " << endl;
}

void MT2Analysis::Begin(const char* filename){
	// Define the output file of histograms
	fHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");

	// book tree
	fMT2tree = new MT2tree();
	BookTree();

	// define which triggers to fill
	if(fisData){
		// HT with dPhi
		fTriggerMap["HLT_HT500_JetPt60_DPhi2p94_v1"] = &fMT2tree->trigger.HLT_HT500_JetPt60_DPhi2p94_v1;
		fTriggerMap["HLT_HT550_JetPt60_DPhi2p94_v1"] = &fMT2tree->trigger.HLT_HT550_JetPt60_DPhi2p94_v1;

		// HT
		fTriggerMap["HLT_HT150_v2"]            = &fMT2tree->trigger.HLT_HT150_v2;
		fTriggerMap["HLT_HT150_v3"]            = &fMT2tree->trigger.HLT_HT150_v3;
		fTriggerMap["HLT_HT160_v2"]            = &fMT2tree->trigger.HLT_HT160_v2;
		fTriggerMap["HLT_HT200_v2"]            = &fMT2tree->trigger.HLT_HT200_v2;
		fTriggerMap["HLT_HT200_v3"]            = &fMT2tree->trigger.HLT_HT200_v3;
		fTriggerMap["HLT_HT240_v2"]            = &fMT2tree->trigger.HLT_HT240_v2;
		fTriggerMap["HLT_HT250_v2"]            = &fMT2tree->trigger.HLT_HT250_v2;
		fTriggerMap["HLT_HT250_v3"]            = &fMT2tree->trigger.HLT_HT250_v3;
		fTriggerMap["HLT_HT260_v2"]            = &fMT2tree->trigger.HLT_HT260_v2;
		fTriggerMap["HLT_HT300_v2"]            = &fMT2tree->trigger.HLT_HT300_v2;
		fTriggerMap["HLT_HT300_v3"]            = &fMT2tree->trigger.HLT_HT300_v3;
		fTriggerMap["HLT_HT300_v4"]            = &fMT2tree->trigger.HLT_HT300_v4;
		fTriggerMap["HLT_HT300_v5"]            = &fMT2tree->trigger.HLT_HT300_v5;
		fTriggerMap["HLT_HT350_v2"]            = &fMT2tree->trigger.HLT_HT350_v2;
		fTriggerMap["HLT_HT350_v3"]            = &fMT2tree->trigger.HLT_HT350_v3;
		fTriggerMap["HLT_HT350_v4"]            = &fMT2tree->trigger.HLT_HT350_v4;
		fTriggerMap["HLT_HT360_v2"]            = &fMT2tree->trigger.HLT_HT360_v2;
		fTriggerMap["HLT_HT400_v2"]            = &fMT2tree->trigger.HLT_HT400_v2;
		fTriggerMap["HLT_HT400_v3"]            = &fMT2tree->trigger.HLT_HT400_v3;
		fTriggerMap["HLT_HT400_v4"]            = &fMT2tree->trigger.HLT_HT400_v4;
		fTriggerMap["HLT_HT400_v5"]            = &fMT2tree->trigger.HLT_HT400_v5;
		fTriggerMap["HLT_HT400_v6"]            = &fMT2tree->trigger.HLT_HT400_v6;
		fTriggerMap["HLT_HT400_v7"]            = &fMT2tree->trigger.HLT_HT400_v7;
		fTriggerMap["HLT_HT400_v8"]            = &fMT2tree->trigger.HLT_HT400_v8;
		fTriggerMap["HLT_HT440_v2"]            = &fMT2tree->trigger.HLT_HT440_v2;
		fTriggerMap["HLT_HT450_v2"]            = &fMT2tree->trigger.HLT_HT450_v2;
		fTriggerMap["HLT_HT450_v3"]            = &fMT2tree->trigger.HLT_HT450_v3;
		fTriggerMap["HLT_HT450_v4"]            = &fMT2tree->trigger.HLT_HT450_v4;
		fTriggerMap["HLT_HT450_v5"]            = &fMT2tree->trigger.HLT_HT450_v5;
		fTriggerMap["HLT_HT450_v6"]            = &fMT2tree->trigger.HLT_HT450_v6;
		fTriggerMap["HLT_HT450_v7"]            = &fMT2tree->trigger.HLT_HT450_v7;
		fTriggerMap["HLT_HT450_v8"]            = &fMT2tree->trigger.HLT_HT450_v8;
		fTriggerMap["HLT_HT500_v2"]            = &fMT2tree->trigger.HLT_HT500_v2;
		fTriggerMap["HLT_HT500_v3"]            = &fMT2tree->trigger.HLT_HT500_v3;
		fTriggerMap["HLT_HT500_v4"]            = &fMT2tree->trigger.HLT_HT500_v4;
		fTriggerMap["HLT_HT500_v5"]            = &fMT2tree->trigger.HLT_HT500_v5;
		fTriggerMap["HLT_HT500_v6"]            = &fMT2tree->trigger.HLT_HT500_v6;
		fTriggerMap["HLT_HT500_v7"]            = &fMT2tree->trigger.HLT_HT500_v7;
		fTriggerMap["HLT_HT500_v8"]            = &fMT2tree->trigger.HLT_HT500_v8;
		fTriggerMap["HLT_HT550_v2"]            = &fMT2tree->trigger.HLT_HT550_v2;
		fTriggerMap["HLT_HT550_v3"]            = &fMT2tree->trigger.HLT_HT550_v3;
		fTriggerMap["HLT_HT550_v4"]            = &fMT2tree->trigger.HLT_HT550_v4;
		fTriggerMap["HLT_HT550_v5"]            = &fMT2tree->trigger.HLT_HT550_v5;
		fTriggerMap["HLT_HT550_v6"]            = &fMT2tree->trigger.HLT_HT550_v6;
		fTriggerMap["HLT_HT550_v7"]            = &fMT2tree->trigger.HLT_HT550_v7;
		fTriggerMap["HLT_HT550_v8"]            = &fMT2tree->trigger.HLT_HT550_v8;
		fTriggerMap["HLT_HT600_v1"]            = &fMT2tree->trigger.HLT_HT600_v1;
		// MHT_HT
		fTriggerMap["HLT_HT250_MHT60_v2"]      = &fMT2tree->trigger.HLT_HT250_MHT60_v2;
		fTriggerMap["HLT_HT250_MHT60_v3"]      = &fMT2tree->trigger.HLT_HT250_MHT60_v3;
		fTriggerMap["HLT_HT250_MHT60_v4"]      = &fMT2tree->trigger.HLT_HT250_MHT60_v4;
		fTriggerMap["HLT_HT250_MHT60_v5"]      = &fMT2tree->trigger.HLT_HT250_MHT60_v5;
		fTriggerMap["HLT_HT250_MHT60_v6"]      = &fMT2tree->trigger.HLT_HT250_MHT60_v6;
		fTriggerMap["HLT_HT250_MHT70_v1"]      = &fMT2tree->trigger.HLT_HT250_MHT70_v1;
		fTriggerMap["HLT_HT250_MHT70_v2"]      = &fMT2tree->trigger.HLT_HT250_MHT70_v2;
		fTriggerMap["HLT_HT250_MHT70_v3"]      = &fMT2tree->trigger.HLT_HT250_MHT70_v3;
		fTriggerMap["HLT_HT250_MHT70_v4"]      = &fMT2tree->trigger.HLT_HT250_MHT70_v4;
		fTriggerMap["HLT_HT250_MHT90_v1"]      = &fMT2tree->trigger.HLT_HT250_MHT90_v1;
		fTriggerMap["HLT_HT250_MHT90_v2"]      = &fMT2tree->trigger.HLT_HT250_MHT90_v2;
		fTriggerMap["HLT_HT260_MHT60_v2"]      = &fMT2tree->trigger.HLT_HT260_MHT60_v2;
		fTriggerMap["HLT_HT300_MHT75_v4"]      = &fMT2tree->trigger.HLT_HT300_MHT75_v4;
		fTriggerMap["HLT_HT300_MHT75_v5"]      = &fMT2tree->trigger.HLT_HT300_MHT75_v5;
		fTriggerMap["HLT_HT300_MHT75_v7"]      = &fMT2tree->trigger.HLT_HT300_MHT75_v7;
		fTriggerMap["HLT_HT300_MHT75_v8"]      = &fMT2tree->trigger.HLT_HT300_MHT75_v8;
		fTriggerMap["HLT_HT300_MHT80_v1"]      = &fMT2tree->trigger.HLT_HT300_MHT80_v1;
		fTriggerMap["HLT_HT300_MHT80_v2"]      = &fMT2tree->trigger.HLT_HT300_MHT80_v2;
		fTriggerMap["HLT_HT300_MHT90_v1"]      = &fMT2tree->trigger.HLT_HT300_MHT90_v1;
		fTriggerMap["HLT_HT300_MHT90_v2"]      = &fMT2tree->trigger.HLT_HT300_MHT90_v2;
		fTriggerMap["HLT_HT350_MHT70_v1"]      = &fMT2tree->trigger.HLT_HT350_MHT70_v1;
		fTriggerMap["HLT_HT350_MHT70_v2"]      = &fMT2tree->trigger.HLT_HT350_MHT70_v2;
		fTriggerMap["HLT_HT350_MHT80_v1"]      = &fMT2tree->trigger.HLT_HT350_MHT80_v1;
		fTriggerMap["HLT_HT350_MHT80_v2"]      = &fMT2tree->trigger.HLT_HT350_MHT80_v2;
		// Muons
		fTriggerMap["HLT_DoubleMu3_HT160_v2"]  = &fMT2tree->trigger.HLT_DoubleMu3_HT160_v2;
		fTriggerMap["HLT_DoubleMu3_v3"]        = &fMT2tree->trigger.HLT_DoubleMu3_v3;
		fTriggerMap["HLT_Mu8_Jet40_v2"]        = &fMT2tree->trigger.HLT_Mu8_Jet40_v2;
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
	if(fJets.size()   > 40) {cout << "ERROR: fJets.size()   > 40: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return false;}
	if(fElecs.size()  > 5 ) {cout << "ERROR: fElecs.size()  >  5: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return false;}
	if(fMuons.size()  > 5 ) {cout << "ERROR: fMuons.size()  >  5: " << "run " << fTR->Run << " Event " << fTR->Event << " skip event" << endl; return false;}

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
		fMT2tree->jet[i].Scale          = fTR->PF2PAT3JScale         [fJets[i]];
		fMT2tree->jet[i].L1FastJetScale = fTR->PF2PAT3JL1FastJetScale[fJets[i]];
		fMT2tree->jet[i].Area           = fTR->PF2PAT3JArea          [fJets[i]];
		fMT2tree->jet[i].isTau          = false;
		if(!fisData){
		fMT2tree->jet[i].Flavour        = fTR->PF2PAT3JFlavour       [fJets[i]];
		}
	}

	// --------------------------------------------------------------
	// match taus to jets
	for(int i=0; i<fTaus.size(); ++i){
		TLorentzVector tau;
		tau.SetPxPyPzE(fTR->PfTau3Px[fTaus[i]],fTR->PfTau3Py[fTaus[i]],fTR->PfTau3Pz[fTaus[i]],fTR->PfTau3E[fTaus[i]]);
		double mindR =10000; int jindex =-1;
		for(int j=0; j<fJets.size(); ++j){
			double dR = tau.DeltaR(fMT2tree->jet[j].lv);
			if(dR < mindR && dR < 0.5) {mindR=dR; jindex = j;} 	
		}
		if(jindex ==-1) continue;
		fMT2tree->jet[jindex].NTauMatch++;
		if(fMT2tree->jet[jindex].isTauMatch && fMT2tree->jet[jindex].TauDR < mindR) continue;
		fMT2tree->jet[jindex].isTauMatch = 1;
		fMT2tree->jet[jindex].TauDR      = mindR;
		fMT2tree->jet[jindex].TauDPt     = fMT2tree->jet[jindex].lv.Pt()-tau.Pt();
	}

	
	// ---------------------------------------------------------------
	// Set NJets, NElecs, NMuons
	fMT2tree->SetNJets         (fJets.size());
	fMT2tree->SetNGenJets      (fTR->NGenJets > gNGenJets ? gNGenJets: fTR->NGenJets);
	fMT2tree->SetNJetsIDLoose  (fMT2tree->GetNjets(20, 2.4, 1));
	fMT2tree->SetNJetsAcc      (fMT2tree->GetNjets(20, 2.4, 0));
	fMT2tree->SetNBJets        (fMT2tree->GetNBtags(3,2.0,20,2.4,1));
	fMT2tree->SetNEles         ((Int_t)fElecs.size());
	fMT2tree->SetNMuons        ((Int_t)fMuons.size());
	fMT2tree->SetNTaus         ((Int_t)fTaus.size());
	
	// --------------------------------------------------------------
	// Fill GenJets
	for(int i=0; i<fTR->NGenJets; ++i){
		if(i >= gNGenJets) {
			cout << "WARNING: Event " << fTR->Event << " Run " << fTR->Run << " has more than " << gNGenJets << " GenJets " << endl;
			continue;
	       	}
		fMT2tree->genjet[i].lv.SetPtEtaPhiE(fTR->GenJetPt[i], fTR->GenJetEta[i], fTR->GenJetPhi[i], fTR->GenJetE[i]);
		double mindR=999.99;
		int    index=-1;
		for(int j=0; j<fMT2tree->NJets; ++j){
			double dR=fMT2tree->jet[j].lv.DeltaR(fMT2tree->genjet[i].lv);
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
		double mass=0;
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
		if(NGenLepts >= 30 ) {cout << "ERROR: NGenLepts >=30: skipping remaining genlepts for event " << fTR->Event << endl; continue;}
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
	fMT2tree->NGenLepts = NGenLepts;

	// --------------------------------------------------------------------
	// MET and MPT
	fMT2tree->pfmet[0]=MET();

	// Fill vector sum of tracks
	TVector3 tracks(0.,0.,0.);
	for(int i=0; i< fTR->NTracks; ++i){
		TVector3 track;
		track.SetPtEtaPhi(fabs(fTR->TrkPt[i]), fTR->TrkEta[i], fTR->TrkPhi[i]);
		tracks += track;
	}	
	fMT2tree->MPT[0].SetXYZM(-tracks.Px(), -tracks.Py(), 0, 0);

	// Pile UP info and reco vertices
	if(!fisData){
		fMT2tree->pileUp.PUnumInt          = fTR->PUnumInteractions;
		fMT2tree->pileUp.PUnumIntLate      = fTR->PUOOTnumInteractionsLate;
		fMT2tree->pileUp.PUnumIntEarly     = fTR->PUOOTnumInteractionsEarly;
		fMT2tree->pileUp.PtHat             = fTR->PtHat;
		fMT2tree->pileUp.isS3              = (int) isS3;
		//////////// S3 vs S4
		if(noPU)
		  fMT2tree->pileUp.Weight            = 1;
	        else if(isS3)
		  fMT2tree->pileUp.Weight            = GetPUWeight(fTR->PUnumInteractions, fTR->PUOOTnumInteractionsLate);
		else
		  fMT2tree->pileUp.Weight            = GetPUWeight(fTR->PUnumInteractions);
		
		fMT2tree->pileUp.Rho               = fTR->Rho;
		
		if(fVerbose > 3) {
			cout << "fTR->PUnumInteractions " <<  fTR->PUnumInteractions << " weight "  
		     	     << " GetPUWeight() " << GetPUWeight(fTR->PUnumInteractions) << endl; 
		}
	}
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
	// HLT triggers
	for (StringBoolMap::iterator iter = fTriggerMap.begin(); iter != fTriggerMap.end(); ++iter){
		if(GetHLTResult(iter->first)) *iter->second =1;
	}
	// _________
	
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
	if(fTR->EcalDeadTPFilterFlag==0) fMT2tree->misc.BadEcalTP=1;
	
	// ------------------------------------------------------------------
	// fill misc 
	fMT2tree->misc.isData              = fisData;
	fMT2tree->misc.Run                 = fTR->Run;
	fMT2tree->misc.Event		   = fTR->Event;
	fMT2tree->misc.LumiSection	   = fTR->LumiSection;
	fMT2tree->misc.LeptConfig          = (int) fLeptConfig;
	fMT2tree->misc.HBHENoiseFlag	   = fTR->HBHENoiseFlag;
	fMT2tree->misc.CrazyHCAL           = fCrazyHCAL;
	fMT2tree->misc.NegativeJEC         = fNegativeJEC;
	fMT2tree->misc.CSCTightHaloID      = fTR->CSCTightHaloID;
	fMT2tree->misc.HT                  = fHT;
	fMT2tree->misc.PFMETsign	   = (MET().Pt())/sqrt(fTR->SumEt);
	
	fMT2tree->misc.MET                 = MET().Pt();
	fMT2tree->misc.METPhi              = MET().Phi();

	fMT2tree->misc.LeadingJPt          = (fMT2tree->NJets > 0) ? fMT2tree->jet[0].lv.Pt() : 0;
	fMT2tree->misc.SecondJPt           = (fMT2tree->NJets > 1) ? fMT2tree->jet[1].lv.Pt() : 0;
	
	// RA2 tracking failure
	fMT2tree->misc.TrackingFailure     = fTR->TrkPtSum/fHT;
	fMT2tree->misc.TrackingFailurePVtx = fTR->PrimVtxPtSum/fHT;

	// ----------------------------------------------------------------
	// that was it
	return true;
}

void MT2Analysis::FillMT2treeCalculations(){
	// -----------------------------------------------------------
	// fill MT2hemi 
	// hemi 0
	// testmass 0, massless pseudojets, no PF-JID, JPt > 20, |jet-eta|<2.4, hemi seed =2 (max inv mass), hemi assoc =3, pf-met, hemi-index 0
	fMT2tree->FillMT2Hemi(0,0,0,20,2.4,2,3,1,0);  
	
	// hemi 1
	// testmass 0, massless pseudojets, PF-JID, JPt > 20, |jet-eta|<2.4, hemi seed =2 (max inv mass), hemi assoc =3, pf-met, hemi-index 1
	fMT2tree->FillMT2Hemi(0,0,1,20,2.4,2,3,1,1);  
	
	// hemi 2
	// testmass 0, massless pseudojets, PF-JID, JPt > 20, |jet-eta|<2.4, minimizing Delta_HT, pf-met, hemi-index 1
	fMT2tree->FillMT2HemiMinDHT(0,0,1,20,2.4,1,2);  
	
	// store MT2 misc variables
	fMT2tree->misc.MT2                 = fMT2tree->hemi[1].MT2;    // note: this is a bit dangerous, 
	fMT2tree->misc.MT2all              = fMT2tree->hemi[0].MT2;    
	fMT2tree->misc.MCT                 = fMT2tree->hemi[1].MCT;
	fMT2tree->misc.AlphaT              = fMT2tree->hemi[2].AlphaT;

	// other variables to be computed based on MT2tree
	// ----------------------------------------------------------------------------------
	fMT2tree->misc.Vectorsumpt	   = fMT2tree->GetMHTminusMET(1, 20, 2.4, true); // including leptons, ID jets only
	fMT2tree->misc.VectorsumptAll	   = fMT2tree->GetMHTminusMET(0, 20, 2.4, true); // including leptons
	fMT2tree->misc.MinMetJetDPhi       = fMT2tree->MinMetJetDPhi(0,20,5.0,1);
	fMT2tree->misc.PassJetID           = fMT2tree->PassJetID(50,2.4,1);
	fMT2tree->misc.PassJetID20         = fMT2tree->PassJetID(20,2.4,1);
	if(fMT2tree->NJets > 0) {
		fMT2tree->misc.Jet0Pass      = (Int_t) fMT2tree->jet[0].IsGoodPFJet(100,2.4,1);
	} else  fMT2tree->misc.Jet0Pass      = 0; 
	if(fMT2tree->NJets > 1) {
		fMT2tree->misc.Jet1Pass      = (Int_t) fMT2tree->jet[1].IsGoodPFJet( 60,2.4,1);
	} else  fMT2tree->misc.Jet1Pass      = 0;
	
	// MHT from jets and taus
	fMT2tree->MHTloose[0]=fMT2tree->GetMHTlv(1, 20, 2.4, true); // only jets satisfying the loose PF-ID and leptons
	fMT2tree->MHT     [0]=fMT2tree->GetMHTlv(0, 20, 2.4, true); // jets and leptons
	
	// DPhiMhtMpt: should be replaced by method to calculate it on the fly. 
	fMT2tree->misc.DPhiMhtMpt=Util::DeltaPhi(fMT2tree->MHT[0].Phi(), fMT2tree->MPT[0].Phi());	
	

	// ---------------------------	
	// calo HT and MHT
	float caHT40=0, caHT40_ID=0;
	TLorentzVector mht40(0,0,0,0), mht40_ID(0,0,0,0);
	for(int j=0; j<fTR->CANJets; ++j){
	  	if( CAJet(j).Pt()<40 || fabs(CAJet(j).Eta())>3.0 ) continue;

	    	caHT40   += CAJet(j).Pt();
	    	mht40    -= CAJet(j);
	  	
	  	if(! (fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001)) continue;
		caHT40_ID += CAJet(j).Pt();
		mht40_ID  -= CAJet(j);
	}
	fMT2tree->misc.caloHT40     = caHT40;
	fMT2tree->misc.caloHT40_ID  = caHT40_ID;
	fMT2tree->misc.caloMHT40    = mht40.Pt();
	fMT2tree->misc.caloMHT40_ID = mht40_ID.Pt();
	fMT2tree->misc.caloHT50     = fCaloHT50;
	fMT2tree->misc.caloHT50_ID  = fCaloHT50_ID;
	fMT2tree->misc.caloMHT30    = fCaloMHT30;
	fMT2tree->misc.caloMHT30_ID = fCaloMHT30_ID;

	// _________
	// stuff for Z->nunu (close your eyes...)	
	vector<int>    jindi;
	vector<double> jpt;
	bool PassJetID_matched(true);
	double HTmatched=0, vectorsumpt_matched_px=0, vectorsumpt_matched_py=0;
	int    NJetsIDLoose_matched=0;
	for(int i=0; i<fTR->PF2PAT3NJets; ++i){
		if(Jet(i).Pt() < 20) continue;  
		bool jet(true);
		for(int gen=0; gen<fMT2tree->NGenLepts; ++gen){
			if( ((abs(fMT2tree->genlept[gen].ID) == 11 || abs(fMT2tree->genlept[gen].ID)==13) && fMT2tree->genlept[gen].MID ==23) ) {
				double deltaR = Util::GetDeltaR(Jet(i).Eta(), fMT2tree->genlept[gen].lv.Eta(), Jet(i).Phi(), fMT2tree->genlept[gen].lv.Phi());
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
				double deltaR = Util::GetDeltaR(CAJet(j).Eta(), fMT2tree->genlept[gen].lv.Eta(), CAJet(j).Phi(), fMT2tree->genlept[gen].lv.Phi());
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
			double dR=Util::GetDeltaR(CAJet(j).Eta(), fMT2tree->ele[e].lv.Eta(), CAJet(j).Phi(), fMT2tree->ele[e].lv.Phi());
			if(dR < 0.4) jet=false;
		}
		for(int m=0; m<fMT2tree->NMuons; ++m){
			double dR=Util::GetDeltaR(CAJet(j).Eta(), fMT2tree->muo[m].lv.Eta(), CAJet(j).Phi(), fMT2tree->muo[m].lv.Phi());
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

	double mindPhi=10;
	if(jindi.size()<1){mindPhi = -999.99;}
	else{
		for(int i=0; i<jindi.size(); ++i){
			if(Jet(jindi[i]).Pt()       < 20 ) continue;
			if(fabs(Jet(jindi[i]).Eta())> 5.0) continue;
			double dphi = TMath::Abs(Jet(jindi[i]).DeltaPhi(met));
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

	fMT2tree->Znunu.NJetsToRemoveEle           = fNJets_toremove_ele;
	fMT2tree->Znunu.NJetsToRemoveMuo           = fNJets_toremove_muo;

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

	vector<double> mutight;
	for(int i=0; i< fTR->PfMu3NObjs; ++i){
		if(std::isnan(fTR->PfMu3Pt[i]))  {fIsNANObj = true; continue;} //protection against objects with NAN-Pt
		fMuons.push_back(i);
		mutight.push_back(fTR->PfMu3Pt[i]);
	}
	fMuons      = Util::VSort(fMuons     , mutight);
	
	vector<double> eltight;
	for(int i=0; i< fTR->PfEl3NObjs; ++i){
		if(std::isnan(fTR->PfEl3Pt[i]))  {fIsNANObj = true; continue;} //protection against objects with NAN-Pt
		fElecs.push_back(i);
		eltight.push_back(fTR->PfEl3Pt[i]);
	}
	fElecs      = Util::VSort(fElecs     , eltight);
	
	vector<double> pt1; 
	for(int ij=0; ij < fTR->PF2PAT3NJets; ++ij){
		if(std::isnan(fTR->PF2PAT3JPt[ij])){ fIsNANObj = true; continue; } //protection against objects with NAN-Pt
		if(Jet(ij).Pt() < 20)       continue;  // note: ETH ntuple only stores PF2PAT3Jets > 15 GeV (defualt config)
		fJets.push_back(ij);                   // fJets has all jets except for duplicates with selected leptons
		pt1.push_back(Jet(ij).Pt());
	}
	fJets        = Util::VSort(fJets,       pt1);
	pt1.clear();

	vector<double> taus;
	for(int i=0; i< fTR->PfTau3NObjs; ++i){
		if(std::isnan(fTR->PfTau3Pt[i])) { fIsNANObj = true; continue;} //protection against objects with NAN-Pt
		if(fTR->PfTau3Pt[i]   < 20    ) continue; // note: taus go up to 2.5 in Eta
		fTaus.push_back(i);
		taus.push_back(fTR->PfTau3Pt[i]);
	}
	fTaus          = Util::VSort(fTaus     , taus);
	
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
	double HT=0;
	for(int j=0; j<fJets.size(); ++j){
		if(Jet(fJets[j]).Pt() > 50 && fabs(Jet(fJets[j]).Eta())<3.0){
			HT += Jet(fJets[j]).Pt();
		}
	}
	fHT = HT;
	if(HT<fCut_HT_min){return false;}
	
	//caloHT and calo MHT
	TLorentzVector mht30(0,0,0,0), mht30_ID(0,0,0,0);
	fCaloHT50 =0.0, fCaloHT50_ID =0.0;
	for(int j=0; j<fTR->CANJets; ++j){
		if( CAJet(j).Pt()<30 || fabs(CAJet(j).Eta())>3.0 ) continue;
	  	// MHT
		mht30 -= CAJet(j);
		// MHT_ID
		if(fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001) mht30_ID -= CAJet(j);	
		// HT
		if( CAJet(j).Pt()<50 ) continue;
		fCaloHT50  += CAJet(j).Pt();
		// HT_ID
		if(fTR->CAJn90[j]>=2 && fTR->CAJEMfrac[j]>=0.000001) fCaloHT50_ID += CAJet(j).Pt();
	}
	fCaloMHT30   =mht30.Pt();
	fCaloMHT30_ID=mht30_ID.Pt();
	if(fCaloHT50      < fCut_caloHT50_min   ) return false;
	if(fCaloHT50_ID   < fCut_caloHT50ID_min ) return false;
	if(fCaloMHT30     < fCut_caloMHT30_min  ) return false;
	if(fCaloMHT30_ID  < fCut_caloMHT30ID_min) return false;

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
	if(fNegativeJEC) {
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
	if(fDoJESUncertainty){
		cout << "--------------------- " << endl;
		cout << " -> initialize JetCorrectionUncertainty " << endl;
		Initialize_JetCorrectionUncertainty();
		cout << "--------------------- " << endl;
	}else{
		fJecUncPF   =NULL; 
 		fJecUncCalo =NULL; 
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

// ****************************************************************************************************
// Jet Selectors

bool MT2Analysis::IsGoodBasicPFJetPAT3(int index, double ptcut, double absetacut){
	// Basic PF jet cleaning and ID cuts
	// cut at pt of ptcut (default = 30 GeV)
	// cut at abs(eta) of absetacut (default = 2.5)
	if(Jet(index).Pt() < ptcut                  ) return false;
	if(fabs(fTR->PF2PAT3JEta[index]) > absetacut) return false;
	if(fTR->PF2PAT3JIDLoose[index]    ==0       ) return false;
	if(fTR->PF2PAT3JScale[index]     < 0        ) return false;
	return true;
}

bool MT2Analysis::IsGoodPFJetMediumPAT3(int index, double ptcut, double absetacut) {
	// Medium PF JID
	if ( ! IsGoodBasicPFJetPAT3(index, ptcut, absetacut)  ) return false;
	if ( !(fTR->PF2PAT3JNeuHadfrac[index] < 0.95)         ) return false;
	if ( !(fTR->PF2PAT3JNeuEmfrac[index]  < 0.95)         ) return false;
	return true;
}

bool MT2Analysis::IsGoodPFJetTightPAT3(int index, double ptcut, double absetacut) {
	// Tight PF JID
	if ( ! IsGoodBasicPFJetPAT3(index, ptcut, absetacut)  ) return false;
	if ( !(fTR->PF2PAT3JNeuHadfrac[index] < 0.90)         ) return false;
	if ( !(fTR->PF2PAT3JNeuEmfrac[index]  < 0.90)         ) return false;
	return true;
}

// Jets and JES uncertainty
void MT2Analysis::Initialize_JetCorrectionUncertainty(){
	string Calo="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/JESuncertainty/Fall10_AK5Calo_Uncertainty.txt";
	string PF  ="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/JESuncertainty/Fall10_AK5PF_Uncertainty.txt";

	ifstream fileCalo(Calo.c_str());
	ifstream filePF  (PF.c_str());
	if(! filePF   ) {cout << "ERROR: cannot open file " << PF     << endl; exit(1); }
	else            {cout << "  using file "            << PF     << endl;          }
	if(! fileCalo ) {cout << "ERROR: cannot open file " << Calo   << endl; exit(1); }
	else            {cout << "  using file "            << Calo   << endl;          } 

	fJecUncCalo = new JetCorrectionUncertainty(Calo);
	fJecUncPF   = new JetCorrectionUncertainty(PF); 
}
// pf
TLorentzVector MT2Analysis::Jet(int index){
	TLorentzVector j(0,0,0,0);
	j.SetPtEtaPhiE(fTR->PF2PAT3JPt[index], fTR->PF2PAT3JEta[index], fTR->PF2PAT3JPhi[index], fTR->PF2PAT3JE[index]);
	if  (! fDoJESUncertainty) return j;
	else                      return MT2Analysis::PFJetJESScaled(j);
}

TLorentzVector MT2Analysis::PFJetJESScaled(TLorentzVector j){
	TLorentzVector j_scaled(0.,0.,0.,0);
	if(fJESUpDown== 1) j_scaled.SetPtEtaPhiM(j.Perp()*(1+GetJECUncertPF(j.Pt(),j.Eta())), j.Eta(), j.Phi(), j.M());
	if(fJESUpDown==-1) j_scaled.SetPtEtaPhiM(j.Perp()*(1-GetJECUncertPF(j.Pt(),j.Eta())), j.Eta(), j.Phi(), j.M());
	return j_scaled;
}

double MT2Analysis::GetJECUncertPF(double pt, double eta){
	// pt must be corrected! not raw  	

	// eta cannot be greater than 5! 
	if      (eta> 5.0) eta = 5.0;
	else if (eta<-5.0) eta =-5.0;

	fJecUncPF->setJetPt(pt);   
	fJecUncPF->setJetEta(eta); 
	double uncert= fJecUncPF->getUncertainty(true);
	if     (fabs(eta) <=1.5)  uncert = sqrt(uncert*uncert + 0.02*0.02);
	else if(fabs(eta) <=3.0)  uncert = sqrt(uncert*uncert + 0.06*0.06);
	else if(fabs(eta) <=5.0)  uncert = sqrt(uncert*uncert + 0.20*0.20);
	return uncert; 
}

// Calo
TLorentzVector MT2Analysis::CAJet(int index){
	TLorentzVector j(0,0,0,0);
	j.SetPtEtaPhiE(fTR->CAJPt[index], fTR->CAJEta[index], fTR->CAJPhi[index], fTR->CAJE[index]);
	if  (! fDoJESUncertainty) return j;
	else                      return MT2Analysis::CAJetJESScaled(j);
}


TLorentzVector MT2Analysis::CAJetJESScaled(TLorentzVector j){
	TLorentzVector j_scaled(0.,0.,0.,0);
	if(fJESUpDown== 1) j_scaled.SetPtEtaPhiM(j.Perp()*(1+GetJECUncertCalo(j.Pt(),j.Eta())), j.Eta(), j.Phi(), j.M());
	if(fJESUpDown==-1) j_scaled.SetPtEtaPhiM(j.Perp()*(1-GetJECUncertCalo(j.Pt(),j.Eta())), j.Eta(), j.Phi(), j.M());
	return j_scaled;
}

double MT2Analysis::GetJECUncertCalo(double pt, double eta){
	// pt must be corrected! not raw  	
	
	// eta cannot be greater than 5! 
	if      (eta> 5.0) eta = 5.0;
	else if (eta<-5.0) eta =-5.0;
	fJecUncCalo->setJetPt(pt);
	fJecUncCalo->setJetEta(eta);
	double uncert= fJecUncCalo->getUncertainty(true);
	if     (fabs(eta) <=1.5)  uncert = sqrt(uncert*uncert + 0.02*0.02);
	else if(fabs(eta) <=3.0)  uncert = sqrt(uncert*uncert + 0.06*0.06);
	else if(fabs(eta) <=5.0)  uncert = sqrt(uncert*uncert + 0.20*0.20);
	return uncert; 
}

// pfMET 
TLorentzVector MT2Analysis::MET(){
	TLorentzVector MET(0,0,0,0);
	if(!fDoJESUncertainty){
		MET.SetPtEtaPhiM(fTR->PFMETPAT, 0., fTR->PFMETPATphi, 0);
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

