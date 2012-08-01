/*********************************************************************************
*  macros to make plots and get reco & acceptance efficiencies                   *
*  for Z->nunu backgorund estimate with Z->ll events                             *
*                                                                                *
*  Pascal, MT2Analysis, Feb. 8th 2011                                            *
*                                                                                *
*  updated: 20110227                                                             *
*********************************************************************************/
#include <iomanip>
#include "../include/helper/Utilities.hh"
#include "../include/MT2tree.hh"

using namespace std;

// USER INPUT *******************************
// D6T or Z2 for Z->ll
int     fTune      = 1;                    // 1=Z2, 2=D6
bool    fSummer11  = true;                 // set to true to use Summer11 for DYtoLL 
bool    fFastSim   = true;                 // set to true to use Summer11 fastsim for Zll
TString fHLTflag   = "HT";                 // HT or MHT_HT
TString fSelection = "whatever";           // only used if fMakeEfficienciesFast==true; 
int     fEffStats  = 10000000;
bool    fMCrescale(false);                  // set to true if using MHT_HT path and need to rescale MC for electrons
bool    fMakeData(true);                    // set to false for MC closure only
bool    fMakePrediction(true);           
bool    fMakeEfficiencies(true);            // this must be set to true 
bool    fMakeEfficienciesPablo(true);      // to use Pablo efficiencies: not fully implemented yet
bool    fUsePabloEffMC(true);              // set to true to use pablo's MC T&P efficiencies
int     fTypeOfEffErr =0;                   // Pablo's errors: not fully implemented yet
bool    fMakeEfficienciesPresel(true);      // to calculate efficiencies based on samples with preselection: always use this!!
bool    fMakeEfficienciesFast(false);       // set to true to take hard coded efficiencies
bool    fMakeTable(true);
bool    fMCStatError(false);
bool    fPileUp                     = false;  // put to true to take into account pileup weights
TString outputdir = "MassPlots/Zinvisible/HT/20110725/highMT2/"; // this is where the plots will go
TString fPablofile_ele_data = "/shome/pnef/SUSY/SUSY_macros/LeptJetMult/efficiencies/effElectron.root"; // pablo's files, currently not used
TString fPablofile_muo_data = "/shome/pnef/SUSY/SUSY_macros/LeptJetMult/efficiencies/effMuon.root";     // pablo's files, currently not used
TString fPablofile_ele_mc   = "/shome/pnef/SUSY/SUSY_macros/LeptJetMult/efficiencies/effElectronMC.root"; // pablo's files, currently not used
TString fPablofile_muo_mc   = "/shome/pnef/SUSY/SUSY_macros/LeptJetMult/efficiencies/effMuonMC.root";     // pablo's files, currently not used
int fVerbose =2;


// ******************************************

// --------------------------------------------------------------------------------
struct leptonicZ{
	double nZ_data;
	double nZ;
	double nZ_bg;
	double nZ_leptveto;
	double R;
	double R_MC;
	double R_err;
	double R_MC_err;
	double rec;
	double rec_MC;
	double rec_err;
	double rec_MC_err;
	double Top_bg;
	double W_bg;
	double QCD_bg;
	double Other_bg;
	double reweight;
	vector<double> pt1;
	vector<double> pt2;
	vector<double> eta1;
	vector<double> eta2;
	vector<double> dR1;
	vector<double> dR2;
	vector<double> eff1;
	vector<double> eff2;
	vector<double> err1;
	vector<double> err2;
	vector<double> weight;
	TH2F* EffHistdR;
	TH2F* EffHistEta;
	TH2F* EffHistdRMC;
	TH2F* EffHistEtaMC;
	static const double eff_err =0.1; // systematic uncertainty on ele_R and muo_R  // fixme
	double P_mll;
	double P_mll_err(){
		return P_mll*0.03;
	}
	double R_err_tot(){ 
		// sys : corresponds to a "eff_err" uncertainty on both reco and acceptance efficiency
		// stat: statistical error from MC statistics
		// tot: squared sum
		return sqrt( pow(eff_err*sqrt(2.)*R,2) +  pow(R_err,2));	
	}
	double R_err_tot_MC(){ 
		// sys : corresponds to a "eff_err" uncertainty on both reco and acceptance efficiency
		// stat: statistical error from MC statistics
		// tot: squared sum
		return sqrt( pow(eff_err*sqrt(2.)*R_MC,2) +  pow(R_MC_err,2));	
	}
	double nZ_scaled(){
		return nZ*global.sigma_Zll_m50 *global.lumi/global.nZllEvents();	
	}
	double nZ_bg_scaled(){
		return nZ_bg*global.sigma_Zll_m50 *global.lumi/global.nZllEvents();
	}
	double nZinv_pred_mll_acc(){
		return ((nZ_scaled()-nZ_bg_scaled())/rec_MC)*global.sigma_Znunu_mll/(global.sigma_Zll_m50*(1./3.)*P_mll)/reweight;
	}
	double nBG_scaled(){
		return nZ_bg_scaled()+Top_bg+W_bg+QCD_bg+Other_bg;
	}
	double pred(bool data){
		if(! data) return (nZ_scaled()-nZ_bg_scaled())*R_MC   /P_mll*global.sigma_Znunu/(global.sigma_Zll_m50*1./3.)/reweight;
		else       return (nZ_data    -nBG_scaled()  )*R      /P_mll*global.sigma_Znunu/(global.sigma_Zll_m50*1./3.)/reweight;
	}
	double pred_err_stat(bool data){
		if(!data) {
			if   (!fMCStatError) return 0.0;
			else                 return 1./P_mll*global.sigma_Znunu/(global.sigma_Zll_m50*1./3.)*sqrt( pow(R_MC*sqrt(nZ-nZ_bg),2) )/reweight;
		}
		else      return 1./P_mll*global.sigma_Znunu/(global.sigma_Zll_m50*1./3.)* 
				 sqrt( pow(R*sqrt(nZ_data),2) )/reweight;
	}
	double pred_err_sys(bool data){
		if(! data) return (nZ_scaled()-nZ_bg_scaled())/(global.sigma_Zll_m50*1./3.)* 
			          sqrt( pow(R_err_tot_MC()/P_mll*global.sigma_Znunu,2) + pow(1/(P_mll*P_mll)*P_mll_err()*R_MC*global.sigma_Znunu,2)+ pow(R_MC/P_mll*global.sigma_Znunu_err,2))/reweight;
		else       return (nZ_data    -nBG_scaled()  )/(global.sigma_Zll_m50*1./3.)*
			          sqrt( pow(R_err_tot()/P_mll*global.sigma_Znunu,2) + pow(1/(P_mll*P_mll)*P_mll_err()*R*global.sigma_Znunu,2) + pow(R/P_mll*global.sigma_Znunu_err,2))/reweight;
	}
} Zee, Zmm;

struct invisibleZ{
	double acc;
	double acc_err;
	double nZ;
	double nZ_scaled(){
		if(fFastSim) return nZ       *33.8*1.06          *global.lumi/global.nZnunuEvents();
		else         return nZ       *global.sigma_Znunu  *global.lumi/global.nZnunuEvents();
	}
	double nZacc_mll;
	double nZacc_mll_scaled(){
		if(fFastSim) return nZacc_mll*33.8*1.06          *global.lumi/global.nZnunuEvents();
		else         return nZacc_mll*global.sigma_Znunu  *global.lumi/global.nZnunuEvents();
	}
} Zinv;

struct global{
	static const double ele_mc_scale    = 0.43; //make sure this is synchronized with the samples.dat
	static const double lower_mass      = 71.;
	static const double upper_mass      = 111.;
	static const double lumi            = 1080.0;
	static const double sigma_Zll_m50   = 3048.;
	static const double sigma_Znunu     = 5760.;
	static const double sigma_Znunu_mll = 5760.*0.9670;
	static const double sigma_Znunu_err = 5760.*0.0707;  // fixme
	double nZnunuEvents(){
		if(fFastSim) return  408999;
		else         return 2165002;
	}
	double nZllEvents(){
		if     (fTune ==1 &&! fSummer11)   return  2329439; // Z2
		else if(fTune ==2 &&! fSummer11)   return  2518269; // D6T
		else if(fTune ==1 &&  fSummer11)   return 34940743; // Summer11 Z2
		else                               return 0;
	}
	double pred_combined(bool data){
		return  (Zee.pred(data)+Zmm.pred(data))/2.;
	}
	double pred_combined_stat_err(bool data){
		double nZee_err=0; double nZmm_err=0; 
		if(data) {nZee_err=sqrt(Zee.nZ_data); nZmm_err=sqrt(Zmm.nZ_data);}
		if(! data && fMCStatError) {nZee_err=sqrt(Zee.nZ-Zee.nZ_bg); nZmm_err=sqrt(Zmm.nZ-Zmm.nZ_bg)}
		double nZee=0;     double nZmm=0;
		if(data) {nZee=Zee.nZ_data/Zee.reweight;     nZmm=Zmm.nZ_data;    }
		else     {nZee=Zee.nZ_scaled()/Zee.reweight; nZmm=Zmm.nZ_scaled();}
		double Rele=(data)?Zee.R:Zee.R_MC;
		double Rmuo=(data)?Zmm.R:Zmm.R_MC;
		return sigma_Znunu/(2./3.*sigma_Zll_m50)*
		       sqrt(pow(nZee_err  *Rele/Zee.P_mll ,2) +
		            pow(nZmm_err  *Rmuo/Zmm.P_mll ,2));
	}
	double pred_combined_sys_err(bool data){
		double nZee_bg=0; double nZmm_bg=0;
		if(data){ nZee_bg = Zee.nBG_scaled()/Zee.reweight  ; nZmm_bg = Zmm.nBG_scaled();   }
		else    { nZee_bg = Zee.nZ_bg_scaled()/Zee.reweight; nZmm_bg = Zmm.nZ_bg_scaled(); }
		double nZee=0; double nZmm=0;
		if(data){ nZee    = Zee.nZ_data/Zee.reweight;        nZmm    = Zmm.nZ_data;        }	
		else    { nZee    = Zee.nZ_scaled()/Zee.reweight;    nZmm    = Zmm.nZ_scaled();    }	
	
		double Rele_err_tot = (data)?Zee.R_err_tot():Zee.R_err_tot_MC();
		double Rmuo_err_tot = (data)?Zee.R_err_tot():Zee.R_err_tot_MC();
		double Rele         = (data)?Zee.R:Zee.R_MC;
		double Rmuo         = (data)?Zmm.R:Zmm.R_MC;
		return 	1/(2./3*sigma_Zll_m50)*
		       	sqrt( pow((nZee-nZee_bg)    *sigma_Znunu*Rele_err_tot/Zee.P_mll ,2) + 
			      pow((nZmm-nZmm_bg)    *sigma_Znunu*Rmuo_err_tot/Zmm.P_mll ,2) +
			      pow(1/(Zee.P_mll*Zee.P_mll)  *Zee.P_mll_err() *(nZee-nZee_bg)    *sigma_Znunu*Rele,2) +
			      pow(1/(Zmm.P_mll*Zmm.P_mll)  *Zmm.P_mll_err() *(nZmm-nZmm_bg)    *sigma_Znunu*Rmuo,2) +
                              pow((nZee-nZee_bg)    *sigma_Znunu_err*Rele/Zee.P_mll ,2) +
                              pow((nZmm-nZmm_bg)    *sigma_Znunu_err*Rmuo/Zmm.P_mll ,2)
			      );
	}
} global;

//__________________________________
const int gNMT2bins                   = 19;
double    gMT2bins[gNMT2bins+1]   = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 135, 180, 260, 360, 500};

const int gNMT2bins_l                   = 14;
double    gMT2bins_l[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500};
//__________________________________



void makeZtoInv_from_Ztoll(){
	gSystem->Load("libPhysics");
	gSystem->CompileMacro("src/MassPlotter.cc", "k");


	if     (fTune ==1) {
		Zmm.P_mll        = 0.9099;  // z2 measured in sample (takes into account bremsstrahlung)
		Zee.P_mll        = 0.8819;  // z2 measured in sample (takes into account bremsstrahlung)
	}
	else if(fTune ==2) {
		Zmm.P_mll        = 0.90684;  // D6T measured in sample (takes into account bremsstrahlung)
		Zee.P_mll        = 0.87737;  // D6T measured in sample (takes into account bremsstrahlung)
	}
	if(fMCrescale){ // to rescale the electron MC yield to account for trigger inefficiency for MHT_HHT path
		if(fHLTflag!="MHT_HT") {cout << "ERROR: you've chosen to reweight the MC yield for Z->ee without HT_HLT path" << endl; exit(1);}
		Zee.reweight=global.ele_mc_scale;
		Zmm.reweight=1.0;
	}else {
		Zee.reweight=1.0;
		Zmm.reweight=1.0;
	}
	if(fTune !=1 && fSummer11) {cout << "ERROR: Summer11 DY is Z2 tune" << endl; exit(1);}

        // ------------------------
	if(fMakeEfficiencies )                                        getEfficiencies(fEffStats);
	if(fMakeEfficiencies &&   fMakeEfficienciesPablo)             readEfficienciesPablo();
	if(fMakeData)                                                 makePlots();
	if(fMakePrediction && fMakeEfficiencies)                      makePrediction(fMakeData);
	if(fMakeTable && fMakePrediction && fMakeEfficiencies)        PrintTable();
        // ------------------------
}

void getEfficiencies(Long64_t nevents){

	// assign stored efficiencies (hard coded)	
	if(fMakeEfficienciesFast){
		AssignEfficiency();	
		return;
	}


	int verbose = 3;

	MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
	tA->setVerbose(verbose);
	if(fMakeEfficienciesPresel){
		if     (fTune ==1 && !fSummer11) tA->init("samples_Zinv/samples_DYJetsToLL_presel.dat");	
		else if(fTune ==1 &&  fSummer11) tA->init("samples_Zinv/samples_DYJetsToLL_presel_Summer11.dat");	
		else if(fTune ==2 && !fSummer11) tA->init("samples_Zinv/samples_DYJetsToLL_D6T_presel.dat");
	}
	else{
		if     (fTune ==1 && !fSummer11) tA->init("samples_Zinv/samples_DYJetsToLL.dat");	
		else if(fTune ==1 &&  fSummer11) tA->init("samples_Zinv/samples_DYJetsToLL_Summer11.dat");	
		else if(fTune ==2 && !fSummer11) tA->init("samples_Zinv/samples_DYJetsToLL_D6T.dat");
	}
		
	tA->PrintZllEfficiency(0, false, "ele", nevents, global.lower_mass, global.upper_mass, fPileUp);	
	tA->PrintZllEfficiency(0, false, "muo", nevents, global.lower_mass, global.upper_mass, fPileUp);	

	if(fMakeEfficienciesPresel){
		if(fFastSim) tA->init("samples_Zinv/samples_ZinvisibleJets_presel_FastSim.dat");	
		else         tA->init("samples_Zinv/samples_ZinvisibleJets_presel.dat");	
	}else{
		if(fFastSim) tA->init("samples_Zinv/samples_ZinvisibleJets.dat");	
		else         tA->init("samples_Zinv/samples_ZinvisibleJets_FastSim.dat");	
	}
	tA->PrintZllEfficiency(0, false, "neutrinos", nevents, global.lower_mass, global.upper_mass, fPileUp);	


	// efficiencies and errors
	Zinv.acc       = tA->fZpred.nu_acc;
	Zinv.acc_err   = tA->fZpred.nu_acc_err;
	
	Zee.rec        = tA->fZpred.ele_reco;     Zee.rec_MC     =Zee.rec;
	Zee.rec_err    = tA->fZpred.ele_reco_err; Zee.rec_MC_err =Zee.rec_err;
	Zee.R          = tA->fZpred.R("ele");     Zee.R_MC       =Zee.R;
	Zee.R_err      = tA->fZpred.R_err("ele"); Zee.R_MC_err   =Zee.R_err;
	Zmm.rec        = tA->fZpred.muo_reco;     Zmm.rec_MC     =Zmm.rec;
	Zmm.rec_err    = tA->fZpred.muo_reco_err; Zmm.rec_MC_err =Zmm.rec_err;
	Zmm.R          = tA->fZpred.R("muo");     Zmm.R_MC       =Zmm.R;
	Zmm.R_err      = tA->fZpred.R_err("muo"); Zmm.R_MC_err   =Zmm.R_err;
	
	cout << "-----------------------------------" << endl;
	cout << "acc efficiency for Z->nunu         " << endl;
	cout << "  acc  efficiency = " << Zinv.acc    << " pm " << Zinv.acc_err  << endl;
	cout << "reco efficiency for Z->ee          " << endl;
	cout << "  reco efficiency = " << Zee.rec     << " pm " << Zee.rec_err << endl;
	cout << "  R=1/(reco*acc)  = " << Zee.R       << " pm " << Zee.R_err    << endl;
	cout << "reco efficiency for Z->mumu        " << endl;
	cout << "  reco efficiency = " << Zmm.rec     << " pm " << Zmm.rec_err << endl;
	cout << "  R=1/(reco*acc)  = " << Zmm.R       << " pm " << Zmm.R_err    << endl;
	cout << "-----------------------------------" << endl;

	delete tA;
}

void makePlots(){
	int verbose = 3;

	MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
	tA->setVerbose(verbose);

	std::ostringstream cutStream;
	cutStream << " " 
	  << "misc.Jet0Pass ==1"                     << "&&"
	  << "misc.Jet1Pass ==1"                     << "&&"
	  << "misc.SecondJPt >100"                   << "&&"
	  << "misc.Vectorsumpt < 70"                 << "&&"
	  << "misc.PassJetID ==1"                    << "&&"
	  << "Znunu.MinMetplusLeptJetDPhiReco > 0.3" << "&&"
	  << "Znunu.METplusLeptsPtReco >30"          << "&&"
	  << "misc.HT > 300"                         << "&&"
	  << "Znunu.caloHT50ID_matchedReco > 600 "   << "&&"
	  << "misc.HBHENoiseFlag == 1"               << "&&"
	  << "misc.CrazyHCAL==0" ;
		

	TString cuts = cutStream.str().c_str();
	TString theCuts;
  	
	std::ostringstream triggerStream;
	triggerStream << "( " << "(trigger.HLT_HT440_v2==1 && misc.Run<=161176)" << "||"
	      << "(trigger.HLT_HT450_v2==1 && (misc.Run>=161217 && misc.Run<=163261))" << "||"
	      << "(trigger.HLT_HT500_v3==1 && (misc.Run>=163270 && misc.Run<=163869))" << "||"
	      << "(trigger.HLT_HT550_v4==1 && (misc.Run>=165088 && misc.Run<=165633))" << "||"
	      << "(trigger.HLT_HT550_v5==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" 
	      << " || (trigger.HLT_HT550_v6==1 && (misc.Run==166346))"                     << "||"
	      << "(trigger.HLT_HT550_v7==1 && (misc.Run>=167078 && misc.Run<=167913))" 
	      << " )";

	TString HLT=""; // NOTE: not actually used as data events are calculated from get_n_events
	if(fHLTflag=="HT") {HLT=triggerStream.str().c_str();}
	//__________________________________________________________________________________________
	// ATTENTION! makePlot adds the underflow and the overflow to the fist and last bin
	//            this is probably not what you want. 
	//            comment out appropriate lines in MassPlotter::makePlot()
	//__________________________________________________________________________________________

	// Z->ee ----------------------------------------------------------------------------------
	TString samples;  // for electrons we need a MC scale factor to account for the trigger inefficiency when using caloMHT. 
	if     (fTune ==1 ) {
		if     (  fSummer11 && ! fFastSim)  samples = "samples_Zinv/samples_MG_V02-01-02_Z2_SUMMER11.dat";
		else if(! fSummer11 &&   fFastSim)  samples = "samples_Zinv/samples_MG_V02-01-02_Z2_FastSim.dat";
		else if(  fSummer11 &&   fFastSim)  samples = "samples_Zinv/samples_MG_V02-01-02_Z2_SUMMER11_FastSim.dat";
		else{
			if(fMCrescale) samples = "samples_Zinv/samples_MG_V02-01-02_Z2_rescaled.dat"; 
			else           samples = "samples_Zinv/samples_MG_V02-01-02_Z2.dat";
		}
	}
	else if(fTune ==2 ) {
		if(fMCrescale) samples = "samples_Zinv/samples_MG_V02-01-02_D6T_rescaled.dat"; 
		else           samples = "samples_Zinv/samples_MG_V02-01-02_D6T.dat";
	}
	tA->init(samples);
	theCuts = cuts+ "&& Znunu.RecoOSee_mll >71 && Znunu.RecoOSee_mll <111";
	//            samples , variable,            cuts,   njet, nlep, HLT,  xtitle       , nbins     , bins     , cleaned, log  , comp ,  ratio, stack, overlay
	tA->makePlot("Znunu.RecoOSee_mll",           theCuts,  -3,  2  , HLT , "Mee"        , 10,  71, 111,          false,   false , true,   true,  true,  false , 1);
	
	// Z->mumu ------------------------------------------------------------------------------------
	TString samples;
	if     (fTune ==1 ) {
		if     (  fSummer11 && ! fFastSim)  samples = "samples_Zinv/samples_MG_V02-01-02_Z2_SUMMER11.dat";
		else if(! fSummer11 &&   fFastSim)  samples = "samples_Zinv/samples_MG_V02-01-02_Z2_FastSim.dat";
		else if(  fSummer11 &&   fFastSim)  samples = "samples_Zinv/samples_MG_V02-01-02_Z2_SUMMER11_FastSim.dat";
		else{
			if(fMCrescale) samples = "samples_Zinv/samples_MG_V02-01-02_Z2_rescaled.dat"; 
			else           samples = "samples_Zinv/samples_MG_V02-01-02_Z2.dat";
		}
	}
	else if(fTune ==2 ) {samples = "samples_Zinv/samples_MG_V02-01-02_D6T.dat";}
	
	tA->init(samples);
	theCuts = cuts+ "&& Znunu.RecoOSmumu_mll >71 && Znunu.RecoOSmumu_mll <111";
	//            samples , variable,            cuts,   njet, nlep, HLT,  xtitle       , nbins     , bins     , cleaned, log  , comp ,  ratio, stack, overlay
	tA->makePlot("Znunu.RecoOSmumu_mll",         theCuts,  -3,  2   ,HLT, "Mmumu"       , 10,  71, 111,          false,   false , true,   true,  true,  false , 1);


	// store results
	Zee.Top_bg   = tA->fZpred.Top_bg_e;			
	Zee.W_bg     = tA->fZpred.W_bg_e;			
	Zee.QCD_bg   = tA->fZpred.QCD_bg_e;			
	Zee.Other_bg = tA->fZpred.Other_bg_e;			
	
	Zmm.Top_bg   = tA->fZpred.Top_bg_mu;			
	Zmm.W_bg     = tA->fZpred.W_bg_mu;			
	Zmm.QCD_bg   = tA->fZpred.QCD_bg_mu;			
	Zmm.Other_bg = tA->fZpred.Other_bg_mu;			

	delete tA;
}

void makePrediction(bool data){
	int counter=0;
	Zinv.nZacc_mll             = get_n_events("Znunu_within_Acc_mll" , false);
	Zinv.nZ                    = get_n_events("Znunu"                , false);
	Zee.nZ                     = get_n_events("Zee"                  , false);
	Zee.nZ_bg                  = get_n_events("Zee"                  , true);
	Zmm.nZ                     = get_n_events("Zmumu"                , false);
	Zmm.nZ_bg                  = get_n_events("Zmumu"                , true);
	if(data){
		double Zee.nZ_data = get_n_events("data_ee"              , false);
		double Zmm.nZ_data = get_n_events("data_mm"              , false);
	}else {
		double Zee.nZ_data =0; 
		double Zmm.nZ_data =0; 
	}
	makePrediction();
}

void makePrediction(){
	

	cout << "_____________________________________________________"              << endl;
	cout << "not rescaled number of MC events                     "              << endl;
	cout << "  the bg corrected number of events divided by the   "              << endl;
	cout << "  reco_efficiency is equal to the number of events   "              << endl;
	cout << "  within acceptance and Mll window.                  "              << endl;
	cout << "  this should be similar for Z->ee and Z->mumu       "              << endl;
	cout << "  (up to differences because of the mll efficiency)  "              << endl;
	cout << "-----------------------------------------------------"              << endl;
	cout << "(nZee-nZee_bg)    /ele_reco = " << (Zee.nZ-Zee.nZ_bg  )  /Zee.rec_MC   << endl;
	cout << "nZee              /ele_reco = " << Zee.nZ /Zee.rec_MC               << endl;
	cout << "nZee                        = " << Zee.nZ                           << endl;
	cout << "nZee_bg                     = " << Zee.nZ_bg                        << endl;
	cout << "(nZmumu-nZmumu_bg)/muo_reco = " << (Zmm.nZ-Zmm.nZ_bg)/Zmm.rec_MC    << endl;
	cout << "nZmumu            /muo_reco = " << Zmm.nZ            /Zmm.rec_MC    << endl;
	cout << "nZmumu                      = " << Zmm.nZ                           << endl;
	cout << "nZmumu_bg                   = " << Zmm.nZ_bg                        << endl;
	cout << "-----------------------------------------------------"              << endl;
	cout << "nZnunu                      = " << Zinv.nZ                          << endl;
	cout << "nZnunu_acc_mll              = " << Zinv.nZacc_mll                   << endl;
	cout << "----------------------------------------------------"               << endl;
	cout << "same numbers as above but rescaled to "<<global.lumi<<"pb-1"        << endl;
	cout << " (# events within mll and accept, corr for reco eff)"               << endl;
	cout << "Zee  : " << (Zee.nZ_scaled()-Zee.nZ_bg_scaled())/Zee.rec_MC         << endl;
	cout << "Zmm  : " << (Zmm.nZ_scaled()-Zmm.nZ_bg_scaled())/Zmm.rec_MC         << endl;
	cout << "____________________________________________________"               << endl;
	cout << "number of Z->nunu events within acceptance and Mll  "               << endl;
	cout << " for " << global.lumi  << " pb -1                          "        << endl;
	cout << "predicted from ee:   " << Zee.nZinv_pred_mll_acc()                  << endl;
	cout << "predicted from mumu: " << Zmm.nZinv_pred_mll_acc()                  << endl;
	cout << "truth nunu:          " << Zinv.nZacc_mll_scaled()                   << endl;
	cout << "____________________________________________________"               << endl;
	cout << "MC N Z->ll events for "<<global.lumi<< " pb-1              "        << endl;
	cout << " this should be compared with the data event yield"                 << endl;
	cout << "nZee     " << Zee.nZ_scaled()     << " nZee_bg   " << Zee.nZ_bg_scaled()    << endl;
	cout << "nZmumu   " << Zmm.nZ_scaled()     << " nZmumu_bg " << Zmm.nZ_bg_scaled()  << endl;
	cout << "nZnunu   " << Zinv.nZ_scaled()                                      << endl;
	cout << "____________________________________________________"               << endl;
	if(fMakeData){
	cout << "Data N Z->ee events in " << global.lumi << " pb-1         "         << endl;
	cout << "nZee_data   "<<Zee.nZ_data  << " bg " << Zee.nBG_scaled()           << endl;
	cout << " bg Z(!ee)  "<<Zee.nZ_bg_scaled()                                   << endl;
	cout << " bg Top     "<<Zee.Top_bg                                           << endl;
	cout << " bg W       "<<Zee.W_bg                                             << endl;
	cout << " bg QCD     "<<Zee.QCD_bg                                           << endl;
	cout << " bg Other   "<<Zee.Other_bg                                         << endl;
	cout << "Data N Z->mm events in " << global.lumi << " pb-1          "        << endl;
	cout << "nZmumu_data "<<Zmm.nZ_data  << " bg " << Zmm.nBG_scaled()           << endl;
	cout << " bg Z(!mm)  "<<Zmm.nZ_bg_scaled()                                   << endl;
	cout << " bg Top     "<<Zmm.Top_bg                                           << endl;
	cout << " bg W       "<<Zmm.W_bg                                             << endl;
	cout << " bg QCD     "<<Zmm.QCD_bg                                           << endl;
	cout << " bg Other   "<<Zmm.Other_bg                                         << endl;
	cout << "____________________________________________________"               << endl;
	}
	cout << "data reco and acceptance efficiencies                    "               << endl;
	cout << "ele_reco " << Zee.rec << " err " << Zee.rec_err                     << endl;
	cout << "muo_reco " << Zmm.rec << " err " << Zmm.rec_err                     << endl;
	cout << "nu_acc   " << Zinv.acc   << " err " << Zinv.acc_err                 << endl;
	cout << "ele_R    " << Zee.R << " ele_R_err " << Zee.R_err                   << endl;
	cout << "muo_R    " << Zmm.R << " muo_R_err " << Zmm.R_err                   << endl;
	cout << "MC reco and acceptance efficiencies                    "            << endl;
	cout << "ele_reco " << Zee.rec_MC << " err " << Zee.rec_MC_err               << endl;
	cout << "muo_reco " << Zmm.rec_MC << " err " << Zmm.rec_MC_err               << endl;
	cout << "nu_acc   " << Zinv.acc   << " err " << Zinv.acc_err                 << endl;
	cout << "ele_R    " << Zee.R_MC << " ele_R_err " << Zee.R_MC_err            << endl;
	cout << "muo_R    " << Zmm.R_MC << " muo_R_err " << Zmm.R_MC_err             << endl;
	cout << "____________________________________________________"               << endl;
	cout << "MC predicted N Z->nunu without acc and mll cuts     "               << endl;
	cout << "  N = " << global.pred_combined(false)            << " pm ";
	cout <<             global.pred_combined_sys_err(false)    << " (sys)  "     << endl;
	cout << " ee channel                                         "               << endl;
	cout << "  N = " << Zee.pred(false)            << " pm ";
	cout <<             Zee.pred_err_sys(false)     << " (sys)  "                << endl;
	cout << " mumu channel                                       "               << endl;
	cout << "  N = " << Zmm.pred(false)           << " pm ";
	cout <<             Zmm.pred_err_sys(false)   << " (sys)  "                  << endl;
	cout << " -------------------------------------------------  "               << endl;
	if(fMakeData){
	cout << "DATA predicted N Z->nunu without acc and mll cuts   "               << endl;
	cout << "  N = " << global.pred_combined(true)            << " pm ";
	cout <<             global.pred_combined_stat_err(true)   << " (stat) pm ";
	cout <<             global.pred_combined_sys_err(true)    << " (sys)  "      << endl;
	cout << " ee channel                                         "               << endl;
	cout << "  N = " << Zee.pred(true)            << " pm ";
	cout <<             Zee.pred_err_stat(true)    << " (stat) pm ";
	cout <<             Zee.pred_err_sys(true)     << " (sys)  "                 << endl;
	cout << " mumu channel                                       "               << endl;
	cout << "  N = " << Zmm.pred(true)           << " pm ";
	cout <<             Zmm.pred_err_stat(true)  << " (stat) pm ";
	cout <<             Zmm.pred_err_sys(true)   << " (sys)  "                   << endl;
	cout << " -------------------------------------------------  "               << endl;
	}
	cout << "expected N Z->nunu from MC                          " << endl;
	cout << "  N = " << Zinv.nZ_scaled()                           << endl;
	cout << "____________________________________________________" << endl;

}


double get_n_events(string process, bool background){
	setStyle();
	TString      basecut   =  " ";
	TChain       *chain    = new TChain("MassTree");
	double       weight    = 1.;
	TString lower_mass_str = TString::Format("%f", global.lower_mass);
	TString upper_mass_str = TString::Format("%f", global.upper_mass);
	TString path_fastsim="/shome/casal/pascual/MT2_V00-06-10/20110730_MC_HT300_data_nocuts/";
	TString path        ="/shome/casal/pascual/MT2_V00-06-10/20110730_MC_HT300_data_nocuts/";
	TString path_data   ="/shome/casal/pascual/MT2_V00-06-10/20110730_MC_HT300_data_nocuts/";

	if(process == "Zee"){
		basecut +="Znunu.RecoOSee_mll > "+lower_mass_str+" && Znunu.RecoOSee_mll < "+upper_mass_str+" &&";
		if(background){
			basecut +="! (Znunu.GenZee_mll_acc > "  +lower_mass_str+" && Znunu.GenZee_mll_acc < "+upper_mass_str+") &&"; 	
		}
		if     (fTune ==1 &&! fSummer11) chain   ->Add(path+"DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root");
		else if(fTune ==1 &&  fSummer11) chain   ->Add(path+"DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11.root");
		else if(fTune ==2)               chain   ->Add(path+"DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola.root");
	} else if ( process == "Zmumu"){
		basecut +="Znunu.RecoOSmumu_mll > "+lower_mass_str+" && Znunu.RecoOSmumu_mll < "+upper_mass_str+" &&";
		if(background){
			basecut +="! (Znunu.GenZmumu_mll_acc > "+lower_mass_str+" && Znunu.GenZmumu_mll_acc < "+upper_mass_str+") &&"; 	
		}
		if     (fTune ==1 &&! fSummer11) chain   ->Add(path+"DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root");
		else if(fTune ==1 &&  fSummer11) chain   ->Add(path+"DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11.root");
		else if(fTune ==2)               chain   ->Add(path+"DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola.root");
	} else if ( process == "Znunu"){
		if(fFastSim) chain   ->Add(path_fastsim+"ZInvisible-HT200toinf-Madgraph-Z2-Fastsim-4_2_7_p1-v4_FlatDist10_2011EarlyData.root");
		else         chain   ->Add(path        +"ZinvisibleJets_7TeV-madgraph.root");
	} else if ( process == "Znunu_within_Acc_mll"){
		if(fFastSim) chain   ->Add(path_fastsim+"ZInvisible-HT200toinf-Madgraph-Z2-Fastsim-4_2_7_p1-v4_FlatDist10_2011EarlyData.root");
		else         chain   ->Add(path        +"ZinvisibleJets_7TeV-madgraph.root");
		basecut +="((Znunu.GenZnunu_e_mll_acc   > "+lower_mass_str+" && Znunu.GenZnunu_e_mll_acc   < "+upper_mass_str+") ||"; 	
		basecut +="( Znunu.GenZnunu_mu_mll_acc  > "+lower_mass_str+" && Znunu.GenZnunu_mu_mll_acc  < "+upper_mass_str+") ||"; 	
		basecut +="( Znunu.GenZnunu_tau_mll_acc > "+lower_mass_str+" && Znunu.GenZnunu_tau_mll_acc < "+upper_mass_str+")) &&"; 	
	} else if ( process == "data_ee" || process == "data_mm"){
		chain   ->Add(path_data+"HT-Run2011A-May10ReReco-v1-AOD.root");
		chain   ->Add(path_data+"HT-Run2011A-PromptReco-v4-AOD_4.root");
	} else {return -1;}


	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// define cuts here
	if(process == "Zee" || process == "Zmumu"){
		basecut +="Znunu.METplusLeptsPt >30 && ";
		basecut +="Znunu.MinMetplusLeptJetDPhi >0.3 && ";
		basecut +="Znunu.caloHT50ID_matched > 600 && ";
	}else if (process == "Znunu" || process == "Znunu_within_Acc_mll"){
		basecut +="misc.MET > 30 &&";
		basecut +="misc.MinMetJetDPhi > 0.3 &&";
		basecut +="misc.caloHT50_ID > 600 &&";
	}else if(process == "data_ee" ){
		basecut +="Znunu.RecoOSee_mll > "+lower_mass_str+" && Znunu.RecoOSee_mll < "+upper_mass_str+" &&";
		basecut +="Znunu.METplusLeptsPtReco > 30 &&";
		basecut +="Znunu.MinMetplusLeptJetDPhiReco >0.3 &&";
		basecut +="Znunu.caloHT50ID_matchedReco > 600 && ";
	}else if(process == "data_mm" ){
		basecut +="Znunu.RecoOSmumu_mll > "+lower_mass_str+" && Znunu.RecoOSmumu_mll < "+upper_mass_str+" &&";
		basecut +="Znunu.METplusLeptsPtReco > 30 &&";
		basecut +="Znunu.MinMetplusLeptJetDPhiReco >0.3 &&";
		basecut +="Znunu.caloHT50ID_matchedReco > 600 && ";
	}
	if(process =="data_mm" || process == "data_ee"){
		std::ostringstream triggerStream;
		triggerStream << "( " << "(trigger.HLT_HT440_v2==1 && misc.Run<=161176)" << "||"
		      << "(trigger.HLT_HT450_v2==1 && (misc.Run>=161217 && misc.Run<=163261))" << "||"
		      << "(trigger.HLT_HT500_v3==1 && (misc.Run>=163270 && misc.Run<=163869))" << "||"
		      << "(trigger.HLT_HT550_v4==1 && (misc.Run>=165088 && misc.Run<=165633))" << "||"
		      << "(trigger.HLT_HT550_v5==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" 
		      << " || (trigger.HLT_HT550_v6==1 && (misc.Run==166346))"                     << "||"
		      << "(trigger.HLT_HT550_v7==1 && (misc.Run>=167078 && misc.Run<=167913))" 
		      << " )";

		TString HLT=""; // NOTE: not actually used as data events are calculated from get_n_events
		if(fHLTflag=="HT") {HLT=triggerStream.str().c_str();}
	}

	basecut  += "misc.Jet0Pass==1 && ";
	basecut  += "misc.Jet1Pass==1 && ";
	basecut  += "misc.SecondJPt >100 &&";
	basecut  += "misc.Vectorsumpt <70 && ";
	basecut  += "misc.HT >300 && ";
	basecut  += "NJetsIDLoose >=3 &&";
	basecut  += "misc.HBHENoiseFlag == 1 &&";
	basecut  += "misc.CrazyHCAL ==0";

	// rescale MC for ee channel
	if(fMCrescale && process == "Zee"){
		weight=weight*global.ele_mc_scale;
	}

	TString selection;
	if(process == "data_ee" || process == "data_mm" || fPileUp==false){
		selection = TString::Format("(%f) * (%s)",weight, basecut.Data());
	}else{
		selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",weight, basecut.Data());
	}

 	double nevents=0;
	TH1D* h;
	if(process=="data_ee" || process=="data_mm" || (process=="Zee"&&!background && fUsePabloEffMC) || (process=="Zmumu" && !background && fUsePabloEffMC)){
		nevents = GetHistoAndEff(chain, selection, process); 
	} else {
		h       = GetHisto(chain, "misc.MET           >>", selection, "h",       0, 1000  , process) ;
	    	cout << "+++ get_n_events: " << process << " background=" << background << " integral " << h->Integral() << endl;	
		nevents = h->Integral();
	}
	
	delete chain;
	delete h;

	return nevents;
}

double   GetHistoAndEff(TChain* chainorig, TString basecut, TString process){
	if(fVerbose >0){
	cout << " +++ Getting n-events and efficiencies for process " << process << endl;
	cout << " +++  applying cuts: " << endl;
	cout << basecut << endl;
	}
	TChain *chain=chainorig->Clone();
	MT2tree *fMT2tree = new MT2tree();
	chain->SetBranchAddress("MT2tree", &fMT2tree);
	Long64_t nentries =  chain->GetEntries();
	Long64_t nbytes = 0, nb = 0;
	int nev =0;
	chain->Draw(">>selList", basecut);

	TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
	chain->SetEventList(myEvtList);
	int counter=0;
	cout << "+++ Filtering done, size=" <<myEvtList->GetN()  << endl;

	if(myEvtList->GetSize()==0) continue;
	cout << "+++ Passing events for " << process << endl;	    
	Zmm.pt1.clear(); Zmm.pt2.clear(); Zmm.eta1.clear(); Zmm.eta2.clear(); Zmm.dR1.clear(); Zmm.dR2.clear(); Zmm.eff1.clear(); Zmm.eff2.clear(); Zmm.err1.clear(); Zmm.err2.clear(); Zmm.weight.clear();
	Zee.pt1.clear(); Zee.pt2.clear(); Zee.eta1.clear(); Zee.eta2.clear(); Zee.dR1.clear(); Zee.dR2.clear(); Zee.eff1.clear(); Zee.eff2.clear(); Zee.err1.clear(); Zee.err2.clear(); Zee.weight.clear();
	while(myEvtList->GetEntry(counter++) !=-1){
		int jentry = myEvtList->GetEntry(counter-1);
		nb =  chain->GetEntry(jentry);   nbytes += nb;
		chain->SetBranchAddress("MT2tree", &fMT2tree);
                  
		if(process=="data_ee" || (process=="Zee" )){
			float eff1=0; float eff2=0;
			float err1=0; float err2=0;
			bool data= (process=="Zee")? false:true;
			if(fMakeEfficienciesPablo) PabloEffEle(fMT2tree->ele[0].lv.Pt(), fMT2tree->ele[0].lv.Eta(), fMT2tree->LeptJetDR(11,0,0,1), fTypeOfEffErr, eff1, err1, data);
			if(fMakeEfficienciesPablo) PabloEffEle(fMT2tree->ele[1].lv.Pt(), fMT2tree->ele[1].lv.Eta(), fMT2tree->LeptJetDR(11,1,0,1), fTypeOfEffErr, eff2, err2, data);
			Zee.pt1 .push_back(fMT2tree->ele[0].lv.Pt());
			Zee.pt2 .push_back(fMT2tree->ele[1].lv.Pt());
			Zee.eta1.push_back(fMT2tree->ele[0].lv.Eta());
			Zee.eta2.push_back(fMT2tree->ele[1].lv.Eta());
			Zee.dR1 .push_back(fMT2tree->LeptJetDR(11,0,0,1));
			Zee.dR2 .push_back(fMT2tree->LeptJetDR(11,1,0,1));
			Zee.eff1.push_back((double) eff1);
			Zee.eff2.push_back((double) eff2);
			Zee.err1.push_back((double) err1);
			Zee.err2.push_back((double) err2);
			Zee.weight.push_back((fPileUp)? fMT2tree->pileUp.Weight:1);
			if(fVerbose >1){
			cout << "++event    " << fMT2tree->misc.Event                                                      << endl;
			cout << "  ele pt1  " << fMT2tree->ele[0].lv.Pt()      << " pt2  " << fMT2tree->ele[1].lv.Pt()      << endl;
			cout << "  ele eta1 " << fMT2tree->ele[0].lv.Eta()     << " eta2 " << fMT2tree->ele[1].lv.Eta()     << endl;
			cout << "  ele dR1  " << fMT2tree->LeptJetDR(11,0,0,1) << " dR2  " << fMT2tree->LeptJetDR(11,1,0,1) << endl;
			cout << "  ele eff1 " << eff1 << " pm " << err1        << " eff2 " << eff2   << " pm " << err2      << endl;
			}
		}else if(process=="data_mm" || (process=="Zmumu")){
			float eff1=0; float eff2=0;
			float err1=0; float err2=0;
			bool data= (process=="Zmumu")? false:true;
			if(fMakeEfficienciesPablo) PabloEffMuo(fMT2tree->muo[0].lv.Pt(), fMT2tree->muo[0].lv.Eta(), fMT2tree->LeptJetDR(13,0,0,1), fTypeOfEffErr, eff1, err1, data);
			if(fMakeEfficienciesPablo) PabloEffMuo(fMT2tree->muo[1].lv.Pt(), fMT2tree->muo[1].lv.Eta(), fMT2tree->LeptJetDR(13,1,0,1), fTypeOfEffErr, eff2, err2, data);
			Zmm.pt1 .push_back(fMT2tree->muo[0].lv.Pt());
			Zmm.pt2 .push_back(fMT2tree->muo[1].lv.Pt());
			Zmm.eta1.push_back(fMT2tree->muo[0].lv.Eta());
			Zmm.eta2.push_back(fMT2tree->muo[1].lv.Eta());
			Zmm.dR1 .push_back(fMT2tree->LeptJetDR(13,0,0,1));
			Zmm.dR2 .push_back(fMT2tree->LeptJetDR(13,1,0,1));
			Zmm.eff1.push_back((double) eff1);
			Zmm.eff2.push_back((double) eff2);
			Zmm.err1.push_back((double) err1);
			Zmm.err2.push_back((double) err2);
			Zmm.weight.push_back((fPileUp)? fMT2tree->pileUp.Weight:1);
			if(fVerbose >1){
			cout << "++event    " << fMT2tree->misc.Event                                                      << endl;
			cout << "  muo pt1  " << fMT2tree->muo[0].lv.Pt()      << " pt2  " << fMT2tree->muo[1].lv.Pt()      << endl;
			cout << "  muo eta1 " << fMT2tree->muo[0].lv.Eta()     << " eta2 " << fMT2tree->muo[1].lv.Eta()     << endl;
			cout << "  muo dR1  " << fMT2tree->LeptJetDR(13,0,0,1) << " dR2  " << fMT2tree->LeptJetDR(13,1,0,1) << endl;
			cout << "  muo eff1 " << eff1 << " pm " << err1        << " muo2 " << eff2   << " pm " << err2      << endl;
			}
		}
	}

	// make new averaged efficienies
	if(fMakeEfficienciesPablo){
		double reco_eff=0;
		double reco_eff1=0;
		double reco_eff2=0;
		double reco_eff_err_squared=0;
		double integral_weight=0;
		for(int i=0; i<myEvtList->GetN(); ++i){
			if(process=="data_ee"){
				reco_eff            +=Zee.eff1[i]*Zee.eff2[i];
				reco_eff_err_squared+= pow(Zee.err1[i]*Zee.eff2[i],2)+pow(Zee.eff1[i]*Zee.err2[i],2);
			}else if(process=="data_mm"){
				reco_eff            +=Zmm.eff1[i]*Zmm.eff2[i];
				reco_eff_err_squared+= pow(Zmm.err1[i]*Zmm.eff2[i],2)+pow(Zmm.eff1[i]*Zmm.err2[i],2);
			}else if(process=="Zee" ){
				reco_eff            +=Zee.weight[i]*Zee.eff1[i]*Zee.eff2[i];
				reco_eff_err_squared+= pow(Zee.weight[i]*Zee.err1[i]*Zee.eff2[i],2)+pow(Zee.weight[i]*Zee.eff1[i]*Zee.err2[i],2);
				integral_weight     +=Zee.weight[i];
				reco_eff1           +=Zee.eff1[i];
				reco_eff2           +=Zee.eff2[i];
			}else if(process=="Zmumu" ){
				reco_eff            +=Zmm.weight[i]*Zmm.eff1[i]*Zmm.eff2[i];
				reco_eff_err_squared+= pow(Zmm.weight[i]*Zmm.err1[i]*Zmm.eff2[i],2)+pow(Zmm.weight[i]*Zmm.eff1[i]*Zmm.err2[i],2);
				integral_weight     +=Zmm.weight[i];
				reco_eff1           +=Zmm.eff1[i];
				reco_eff2           +=Zmm.eff2[i];
			}
		}
		if(process=="data_ee" ){
			cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
			Zee.rec     = reco_eff/myEvtList->GetN();		
			Zee.rec_err = sqrt(reco_eff_err_squared)/myEvtList->GetN();		
			Zee.R       = 1/(Zee.rec*Zinv.acc);
			Zee.R_err   =Zee.R*sqrt(pow(Zee.rec_err/Zee.rec,2)+pow(Zinv.acc_err/Zinv.acc,2));
			cout << " + tag&probe data efficiencies               " << endl;
			cout << " Zee.rec " << Zee.rec  << " pm " << Zee.rec_err << endl;
			cout << " Zee.R "   << Zee.R    << " pm " << Zee.R_err   << endl;
			cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
		}else if(process=="data_mm"){
			cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
			Zmm.rec     = reco_eff/myEvtList->GetN();
			Zmm.rec_err = sqrt(reco_eff_err_squared)/myEvtList->GetN();		
			Zmm.R       = 1/(Zmm.rec*Zinv.acc);
			Zmm.R_err   =Zmm.R*sqrt(pow(Zmm.rec_err/Zmm.rec,2)+pow(Zinv.acc_err/Zinv.acc,2));
			cout << " + tag&probe data efficiencies               " << endl;
			cout << " Zmm.rec " << Zmm.rec  << " pm " << Zmm.rec_err << endl;
			cout << " Zmm.R "   << Zmm.R    << " pm " << Zmm.R_err   << endl;
			cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
		}else if(process=="Zee" ){
			cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
			reco_eff1      = reco_eff1/myEvtList->GetN();
			reco_eff2      = reco_eff2/myEvtList->GetN();
			Zee.rec_MC     = reco_eff1*reco_eff2;
//			Zee.rec_MC     = reco_eff/myEvtList->GetN();
			Zee.rec_MC_err = sqrt(reco_eff_err_squared)/myEvtList->GetN();		
			Zee.R_MC       = 1/(Zee.rec_MC*Zinv.acc);
			Zee.R_MC_err   = Zee.R_MC*sqrt(pow(Zee.rec_MC_err/Zee.rec_MC,2)+pow(Zinv.acc_err/Zinv.acc,2));
			cout << " + tag&probe MC efficiencies               " << endl;
			cout << " Zee.rec_MC " << Zee.rec_MC  << " pm " << Zee.rec_MC_err << endl;
			cout << " Zee.R_MC "   << Zee.R_MC    << " pm " << Zee.R_MC_err   << endl;
			cout << " average pileup weight " << integral_weight/myEvtList->GetN();
			cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
		}else if(process=="Zmumu"){
			cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
			reco_eff1      = reco_eff1/myEvtList->GetN();
			reco_eff2      = reco_eff2/myEvtList->GetN();
			Zmm.rec_MC     = reco_eff1*reco_eff2;
//			Zmm.rec_MC     = reco_eff/myEvtList->GetN();
			Zmm.rec_MC_err = sqrt(reco_eff_err_squared)/myEvtList->GetN();		
			Zmm.R_MC       = 1/(Zmm.rec_MC*Zinv.acc);
			Zmm.R_MC_err   = Zmm.R_MC*sqrt(pow(Zmm.rec_MC_err/Zmm.rec_MC,2)+pow(Zinv.acc_err/Zinv.acc,2));
			cout << " + tag&probe MC efficiencies               " << endl;
			cout << " Zmm.rec_MC " << Zmm.rec_MC  << " pm " << Zmm.rec_MC_err << endl;
			cout << " Zmm.R_MC "   << Zmm.R_MC    << " pm " << Zmm.R_MC_err   << endl;
			cout << " average pileup weight " << integral_weight/myEvtList->GetN();
			cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
		}
	}

	delete fMT2tree;
	delete chain;
	cout << "done " << endl;
	return myEvtList->GetN();
}

TH1D* GetHisto(TChain* chain, TString var, TString basecut, TString name, double binmin, double binmax, TString process){
	TH1D* h1 = new TH1D(name,"", 100, binmin, binmax);
	h1->SetXTitle(var.Data());
	h1->SetMarkerSize(4);
	h1->SetMarkerColor(kBlue);
	
	TString varname = var + h1->GetName();
	cout << "+++ Getting N-Events for process:  " << process << "\n"
	     << "+++  applying cuts" << "\n"
	     << basecut << endl;
	int n = chain->Draw(varname,basecut,"");
	
	return h1;
}

void setStyle(){
//	gROOT->ProcessLine(".x ~casal/SetStyle_PRD.C");
	gStyle->SetLegendBorderSize(0);
	gStyle->SetFillStyle(0);
	gStyle->SetTextFont(62);
	gStyle->SetTextSize(0.045);
	gStyle->SetPalette(1);
}

void AssignEfficiency(){
	if(fTune==2 && fSelection == "nocuts"){
		Zee.rec   =0.656516,   Zee.rec_err  =0.00165482,  Zee.R=2.94105,    Zee.R_err=0.00818553; // no cuts D6, pile up reweighted
		Zmm.rec   =0.783602,   Zmm.rec_err  =0.00184442,  Zmm.R=2.46407,    Zmm.R_err=0.00648804;  // no cuts D6
		Zinv.acc  =0.517907,   Zinv.acc_err =0.00061121;                                         // no cuts D6
	}
	else if(fTune==2 && fSelection == "HT300"){
		Zee.rec   =0.586667,   Zee.rec_err  =0.0454812,  Zee.R=2.07783,    Zee.R_err=0.177532; // no cuts D6
		Zmm.rec   =0.660173,   Zmm.rec_err  =0.0487062,  Zmm.R=1.84648,    Zmm.R_err=0.151516;  // no cuts D6
		Zinv.acc  =0.820349,   Zinv.acc_err =0.0294654;                                         // no cuts D6
	}
	else if(fTune==2 && fSelection == "HT300_allcuts"){
		Zee.rec   =0.637387,   Zee.rec_err  =0.0592846,  Zee.R=1.96201,    Zee.R_err=0.203794; // no cuts D6
		Zmm.rec   =0.702365,   Zmm.rec_err  =0.0653043,  Zmm.R=1.78049,    Zmm.R_err=0.184886;  // no cuts D6
		Zinv.acc  =0.799644,   Zinv.acc_err =0.036973;                                         // no cuts D6
	}
	else if(fTune==2 && fSelection == "caloMHT130_caloHT320_HT300_allcuts"){
		Zee.rec   =0.762341,   Zee.rec_err  =0.107947,  Zee.R=1.49725,    Zee.R_err=0.236916; // no cuts D6
		Zmm.rec   =0.698309,   Zmm.rec_err  =0.0949222, Zmm.R=1.63454,    Zmm.R_err=0.250384;  // no cuts D6
		Zinv.acc  =0.876107,   Zinv.acc_err =0.061874;                                         // no cuts D6
	}
	else{ cout << "ERROR: AssignEfficiency" << endl; exit(1);}
}

void PrintTable(){
	int precision=4;
	if(fSelection=="nocuts") precision=6;
	cout << "*********************************************************************" << endl;
	cout << "\%BEGINLATEX\%"                                                                                                                      << endl;
	cout << "\\begin{table}"                                                                                                                      << endl
	     << "\\begin{center}"                                                                                                                     << endl
             << "\\begin{tabular}{lcc}"                                                                                                               << endl	     
	     << "\\hline\\hline"                                                                                                                      << endl
	     << "                               &  $"     << "Z\\to e \ e"          << "$ & $"          << " Z\\to \\mu \\mu$          \\\\"          << endl
	     << "\\hline"                                                                                                                             << endl;
	if(fMakeData){
	cout << " $N_{ll}^{\\mathrm{data}}$      &  $"    << Zee.nZ_data            << "$ & $"          << Zmm.nZ_data            << "$ \\\\"         << endl
	     << "\\hline"                                                                                                                             << endl;
	}
	cout << setprecision(precision) 
	     << " $N_{ll}^{\\mathrm{MC}}(Z)$     &  $"    << Zee.nZ_scaled()-Zee.nZ_bg_scaled()
     	                                                                            << "$ & $" << Zmm.nZ_scaled()-Zmm.nZ_bg_scaled()         
									                                                          << "$ \\\\"         << endl
	     << " $N_{ll}^{\\mathrm{MC}}(BG)$    &  $"    << Zee.nBG_scaled()       << "$ & $" << Zmm.nBG_scaled()                << "$ \\\\"         << endl
	     << "\\hline"                                                                                                                             << endl
	     << setprecision(4) 
	     << " $\\varepsilon_{acc}$           &        \\multicolumn{2}{c}{"     << Zinv.acc                                   <<  "}   \\\\ "     << endl
	     << " $\\varepsilon_{reco}$          &  $"    << Zee.rec                << "$ & $" << Zmm.rec                         << "$ \\\\"         << endl
	     << "\\hline"                                                                                                                             << endl
	     << setprecision(precision) 
	     << " $Z_{\\nu \ \\nu}^{\\mathrm{pred}} $ MC    &  $" << Zee.pred(false)        
	                                                  << " \\pm "               << Zee.pred_err_stat(false)   
						          << " \\ (stat) \\ \\pm "  << Zee.pred_err_sys(false)                     << " \\ (sys) \\"                 
	                                                  << " $ & $"               << Zmm.pred(false)        
	                                                  << " \\ \\pm "            << Zmm.pred_err_stat(false)   
						          << " \\ (stat) \\ \\pm "  << Zmm.pred_err_sys(false)                     << "\\ (sys) $ \\\\" << endl;
     	if(fMakeData){	
	cout << " $Z_{\\nu \ \\nu}^{\\mathrm{pred}} $ data  &  $" << Zee.pred(true)        
	                                                  << " \\pm "               << Zee.pred_err_stat(true)   
						          << " \\ (stat) \\ \\pm "  << Zee.pred_err_sys(true)                     << " \\ (sys) \\"                 
	                                                  << " $ & $"               << Zmm.pred(true)        
	                                                  << " \\ \\pm "            << Zmm.pred_err_stat(true)   
						          << " \\ (stat) \\ \\pm "  << Zmm.pred_err_sys(true)                     << "\\ (sys) $ \\\\" << endl   
	     << "\\hline"                                                                                                                              << endl
	     << " $Z_{\\nu \ \\nu}^{\\mathrm{pred}} $ data,\\ combined  & \\multicolumn{2}{c}{$" << global.pred_combined(true) 
	                                                  << " \\pm "               << global.pred_combined_stat_err(true)
	                                                  << " \\ (stat) \\ \\pm "  << global.pred_combined_sys_err(true)         << "\\ (sys) $} \\\\"<< endl;
	}
	cout << "\\hline"                                                                                                                              << endl;
	if(fFastSim){
	cout << setprecision(precision) 
	     << " $Z_{\\nu \\nu} $ MC (FastSim)         &        \\multicolumn{2}{c}{"    << Zinv.nZ_scaled()                           <<  "}   \\\\ "      << endl;
	}else {
	cout << setprecision(precision) 
	     << " $Z_{\\nu \\nu} $ MC (D6T)         &        \\multicolumn{2}{c}{"    << Zinv.nZ_scaled()                           <<  "}   \\\\ "      << endl;
	}
	cout << "\\hline\\hline"                                                                                                                       << endl
	     << "\\end{tabular}"                                                                                                                       << endl
	     << "\\end{center}"                                                                                                                        << endl
	     << "\\end{table}"                                                                                                                         << endl
	     << "\%ENDLATEX\%"                                                                                                                         << endl
	     << endl;

	cout << "*********************************************************************" << endl;


}

void readEfficienciesPablo(){	
	// electros

	cout << "--------------------------------------" << endl;
	cout << "Read Pablo EffTrees " << endl;	
  	TFile* fele   = TFile::Open(fPablofile_ele_data);
  	TFile* feleMC = TFile::Open(fPablofile_ele_mc);
  	if (!fele || !feleMC) {
		cout << "ERROR : Could not open file " << fPablofile_ele_data << " or " fPablofile_ele_mc<< "!"  << endl;
	    	exit(1);
  	}
  
  	fele->GetObject("histo_electron_var1pt_var2dr_data",Zee.EffHistdR) ;
  	fele->GetObject("histo_electron_var1pt_var2eta_data",Zee.EffHistEta);
  	feleMC->GetObject("histo_electron_var1pt_var2dr_mc",Zee.EffHistdRMC) ;
  	feleMC->GetObject("histo_electron_var1pt_var2eta_mc",Zee.EffHistEtaMC);
	if(Zee.EffHistdR==NULL || Zee.EffHistEta==NULL || Zee.EffHistdRMC==NULL || Zee.EffHistEtaMC==NULL) {cout << " cannot get hist" << endl; exit(-1);}

	// muons
	TFile* fmuo   = TFile::Open(fPablofile_muo_data);
	TFile* fmuoMC = TFile::Open(fPablofile_muo_mc);
  	if (!fmuo || !fmuoMC) {
		cout << "ERROR : Could not open file " << fPablofile_muo_data << " or " << fPablofile_muo_mc << "!"  << endl;
	    	exit(1);
  	}
  
  	fmuo  ->GetObject("histo_muon_var1pt_var2dr_data" ,Zmm.EffHistdR) ;
  	fmuo  ->GetObject("histo_muon_var1pt_var2eta_data",Zmm.EffHistEta);
  	fmuoMC->GetObject("histo_muon_var1pt_var2dr_mc" ,Zmm.EffHistdRMC) ;
  	fmuoMC->GetObject("histo_muon_var1pt_var2eta_mc",Zmm.EffHistEtaMC);
	if(Zmm.EffHistdR==NULL || Zmm.EffHistEta==NULL || Zmm.EffHistdRMC==NULL || Zmm.EffHistdRMC==NULL) {cout << " cannot get hist" << endl; exit(-1);}

	cout << " done"                                  << endl;
	cout << "--------------------------------------" << endl;
}

void PabloEffEle(float pt_, float eta_, float dr_, int typeOfError, float &eff, float &err, bool data) {
	if(! fMakeEfficienciesPablo) return;

	float efficiency,  e_efficiency1, e_efficiency2;
	float pt, eta, dr, njets;
	float ptwidth, etawidth, drwidth, njetswidth;
	TH2F* H_eta=(data)? Zee.EffHistEta->Clone(): Zee.EffHistEtaMC->Clone();
	TH2F* H_dR =(data)? Zee.EffHistdR ->Clone(): Zee.EffHistdRMC ->Clone();
	if(pt_ < 20){
		int bin = 	H_eta->FindBin(pt_, eta_);
		if(bin!=1){ 
			eff = H_eta->GetBinContent( bin  );
			err = H_eta->GetBinError( bin  );
		}
		else exit(-1);
	}
	else{
		int bin = -1;
		float tmpDr=dr_, tmpPt=pt_;
		// for leptons with dR<0.5, I'm taking dR=0.5. Same for pt
		if(dr_ <0.5) tmpDr = 0.5;
		if(pt_ >200) tmpPt = 199;
		bin = H_dR->FindBin(tmpPt, tmpDr );
		if(bin!=-1){ 
			eff = H_dR->GetBinContent( bin  );
			err = H_dR->GetBinError( bin  );
		}
		else exit(-1);
	}

	if(eff==0) {cout << "eff not found " << endl; exit(1);}

}

void PabloEffMuo(float pt_, float eta_, float dr_, int typeOfError, float &eff, float &err, bool data) {
	if(! fMakeEfficienciesPablo) return;

	float efficiency,  e_efficiency1, e_efficiency2;
	float pt, eta, dr, njets;
	float ptwidth, etawidth, drwidth, njetswidth;
	TH2F* H_eta=(data)? Zmm.EffHistEta->Clone(): Zmm.EffHistEtaMC->Clone();
	TH2F* H_dR =(data)? Zmm.EffHistdR ->Clone(): Zmm.EffHistdRMC ->Clone();


	if(pt_ < 20){
		int bin = 	H_eta->FindBin(pt_, eta_);
		if(bin!=1){ 
			eff = H_eta->GetBinContent( bin  );
			err = H_eta->GetBinError( bin  );
		}
		else exit(-1);
	}
	else{
		int bin = -1;
		float tmpDr=dr_, tmpPt=pt_;
		// for leptons with dR<0.5, I'm taking dR=0.5. Same for pt
		if(dr_ <0.5) tmpDr = 0.5;
		if(pt_ >200) tmpPt = 199;
		bin = H_dR->FindBin(tmpPt, tmpDr );
		if(bin!=-1){ 
			eff = H_dR->GetBinContent( bin  );
			err = H_dR->GetBinError( bin  );
		}
		else exit(-1);
	}

	if(eff==0) {cout << "eff not found " << endl; exit(1);}

}
