#include "TEfficiency.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TTree.h"
#include "TCut.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TMap.h"
#include "TLegend.h"
#include "TLatex.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string>
#include <cmath>
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"//use own path of MT2tree.hh
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/Utilities.hh"

//run via root -l -b -q TTbarStudiesHiggs.C++

//this function is basically a copy of TTbarStudies.C
//but shapes are done for Mbb instead of MT2, event selection is the MT2 Higgs selection instead of inclusive MT2 analysis
//therefore no additional comments are given


using namespace std;

void load(const char* filename = "samples/datasamples/samples_2141_dataonly.dat");
void TTbarStudiesHiggs();

struct sample{
	TString name;
	TString sname;
	TString shapename;
	TString type;
	TFile *file;
	TTree *tree;
	float xsection;
	float nevents;
	float kfact;
	float PU_avg_weight;
	float lumi;
	int color;
};



vector < sample >  fSamples;
const int fVerbose = 3; // already defined below
TString fPath;

Bool_t fbTagReweight    = true;//for all shapes (except BTV SF systematics) apply BTV SF weights
Bool_t fISRreweight     = true;//for nominal shape (and successive other shapes except ISR systematics) apply the 'ISR reweighting', default = true, but consider also false
Bool_t WnoScaleMatching = true;//do not compute scale/matching up/down (default == true due to statistics of those W samples)

double ISRweights[4] = {1.00,0.95,0.90,0.80};
double ISupperpts[4] = {120.,150.,250.,99999.};

const int sampletypesize = 34;
string sample_type[sampletypesize] = {"Nominal", "MatchingUp", "MatchingDown", "ScaleUp", "ScaleDown", "MassUp", "MassDown", "JESUp", "JESDown", "METUp", "METDown", "BSFUp", "BSFDown", "ISRUp", "ISR", "ISRDown", "PUUp", "PUDown", "TopPtReweighted", "nominal0l", "nominal1l", "nominal2l", "nominal1Sample", "nominalPowheg", "WUp", "WDown", "TopUp", "TopDown", "SingleTopUp", "SingleTopDown", "TTVUp", "TTVDown", "MT2relaxedcut", "MT2fullcut"};//maybe not all filled  - ! include TTbar and WJets
const int samplekindsize = 6;
string sample_kind[samplekindsize] = {"allMC","allTop", "TTbar", "SingleTop", "TTV", "WJets"};
const int HTregionsize = 2;
string HT_region[HTregionsize] = {"lowHT", "highHT"};

void TTbarStudiesHiggs(){


	TFile *pu = TFile::Open("/shome/haweber/MT2Analysis_8TeV/Code/Certification/pileup_weights.root");
	TH1D *pu_up = (TH1D*)pu->Get("pileup_weight_xsecPlusFivePercent");
	TH1D *pu_down = (TH1D*)pu->Get("pileup_weight_xsecMinusFivePercent");

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	TString  outputdir = "Filtered/TTbarStudies/Higgs/";
    	Util::MakeOutputDir(outputdir);

	TString outputname = "TTbarStudiesHistograms_all_noWscaleupdown_limitMbb.root";
	if(!WnoScaleMatching) outputname = "TTbarStudiesHistograms_all_limitMbb.root";
	if(fISRreweight){
		outputname = "TTbarStudiesHistograms_all_noWscaleupdown_ISR_limitMbb.root";
		if(!WnoScaleMatching) outputname = "TTbarStudiesHistograms_all_ISR_limitMbb.root";
	}

	TString  samples = "samples/samples_TTbar_filter.dat";
	map<string, TH1D*>    histos;
	map<string, TH1D*>    histosnormalized;
	for(int i0 = 0; i0<samplekindsize;   ++i0){
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i3 = 0; i3<HTregionsize;     ++i3){

		string hs = string("_") + sample_type[i1] + string("_") + sample_kind[i0] + string("_") + HT_region[i3];
		string mapname = "MT2" + hs;//called MT2, but is Mbb
		if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", 33,5,500);
		Color_t colour = 1;
		histos[mapname]->SetLineColor(colour); histos[mapname]->SetMarkerColor(colour); histos[mapname]->SetLineWidth(2);
		histos[mapname]->GetXaxis()->SetTitle("M_{bb} [GeV]"); histos[mapname]->GetYaxis()->SetTitle("Events"); 

	}}}//}

	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Sumw2();}

	std::ostringstream cutStream;
	std::ostringstream cutStreamBase;
	cutStream << " " 
		<< "NJetsIDLoose>=2"                        << "&&"//preselection only
		<< "NTausIDLoose3Hits==0"                   << "&&"
		<< "NEles==0"                               << "&&"
		<< "NMuons==0"                              << "&&"
		<< "misc.HT>400" << "&&"//;//preselection only
		<< "misc.Jet0Pass ==1"                      << "&&"
		<< "misc.Jet1Pass ==1"                      << "&&"
		<< "misc.Vectorsumpt < 95";//precut
	cutStream  << "&& misc.MinMetJetDPhi4Pt40 >0.3";//have nowhere angular uncertainties
	cutStream << "&&((misc.HT>450&&misc.HT<750&&misc.MET>10)||(misc.HT>750&&misc.MET>10))";
	cutStream <<"&&NJetsIDLoose40 >= 4"<<"&&"
		  <<"NBJetsCSVM>=2";
	cutStream << "&&GetSelBBMinv()>=20&&GetSelBBMinv()<=200";

	cutStreamBase << " " 
      << "misc.PassJet40ID ==1"                      << "&&"
      << "(misc.HBHENoiseFlag == 0 || misc.ProcessID==10)"  << "&&" // for rare SM samples
      << "misc.CSCTightHaloIDFlag == 0"             << "&&"
      << "misc.trackingFailureFlag==0"              << "&&"
      << "misc.eeBadScFlag==0"                      << "&&"
      << "misc.EcalDeadCellTriggerPrimitiveFlag==0" << "&&"
      << "misc.TrackingManyStripClusFlag==0"             << "&&"
      << "misc.TrackingTooManyStripClusFlag==0"          << "&&"
      << "misc.TrackingLogErrorTooManyClustersFlag==0"   << "&&"
      << "misc.CrazyHCAL==0";
      cutStreamBase << "&&misc.MET/misc.CaloMETRaw<=2.";

	TString cuts = cutStream.str().c_str();
	TString basecuts = cutStreamBase.str().c_str();

	load(samples.Data());

   	for(size_t i = 0; i < fSamples.size(); ++i){
        
		TString kindsample = fSamples[i].sname;
		kindsample.ReplaceAll("_lowHT","");
		kindsample.ReplaceAll("_highHT","");

		TString namesample = fSamples[i].name;
		namesample.ReplaceAll("_lowHT","");
		namesample.ReplaceAll("_highHT","");
		namesample.ReplaceAll("v1","");
		namesample.ReplaceAll("v2","");

	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
	    if(fVerbose>2) cout << "SampleName: looping over " << fSamples[i].name << " (" << namesample << ") added in " << kindsample << endl;
	    if(fVerbose>2) cout << "            sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;

		string nominal = " ";
		if(namesample=="TTbar_Madgraph0l" || namesample=="TTbar_Madgraph1l" || namesample=="TTbar_Madgraph2l")                          nominal = "TTbar";
		if(namesample=="WJetsHT200"||namesample=="WJetsHT250"||namesample=="WJetsHT300"||namesample=="WJetsHT400")                      nominal = "WJets";
		if(namesample=="T_s"||namesample=="Tbar_s"||namesample=="T_t"||namesample=="Tbar_t"||namesample=="T_tW"||namesample=="Tbar_tW") nominal = "SingleTop";
		if(namesample=="WW"||namesample=="WZ"||namesample=="ZZ"||namesample=="TTW"||namesample=="TTZ"||namesample=="TTG")               nominal = "TTV";
		cout << "this sample is " << nominal << " of kind " << kindsample << " with name " << namesample << endl;

	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;
        
	    TString myCuts = cuts + "&&" + basecuts;

   	    cout << "Cuts for Flow: " << myCuts << endl;
   	    fSamples[i].tree->Draw(">>selList", myCuts);
        
	    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
	    fSamples[i].tree->SetEventList(myEvtList);
	    int counter=0;
	    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
	    if(myEvtList->GetSize()==0) continue;
        while(myEvtList->GetEntry(counter++) !=-1){	
		int jentry = myEvtList->GetEntry(counter-1);
            
		nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
		fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
		
		if ( fVerbose>2 && counter % 50000 == 0 )  cout << "+++ Proccessing event " << counter << endl;
		
		Double_t weight = sample_weight;
		if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 

		if(fMT2tree->misc.HT>=750.&&fMT2tree->misc.HT<1200.&&fMT2tree->NJetsIDLoose40==2&&fMT2tree->NBJets40CSVM>=1&&fMT2tree->misc.MT2>=135.&&fMT2tree->misc.MT2<=170.&&nominal=="SingleTop"&&fMT2tree->misc.LumiSection==149&&fMT2tree->misc.Event==44699) continue;//sample T_t_highHT

		//these are the values that can change
		int NumBJets  = fMT2tree->NBJets40CSVM;
		int NumJets   = fMT2tree->NJetsIDLoose40;
		float HTval   = fMT2tree->misc.HT;
		float METval  = fMT2tree->misc.MET;
		float MinDPhi = fMT2tree->misc.MinMetJetDPhi4Pt40;
		float VSPT    = fMT2tree->misc.Vectorsumpt;
		float J1Pt    = fMT2tree->misc.LeadingJPt;
		float J2Pt    = fMT2tree->misc.SecondJPt;
		float MT2val  = fMT2tree->misc.MT2;
		float Mbbval  = fMT2tree->GetSelBBMinv();
		if(Mbbval<0) continue;

            Double_t ISRweight(1.); Double_t weightnoISR = weight;
            if(!fMT2tree->misc.isData){
                TLorentzVector hardgenlv; hardgenlv.SetPtEtaPhiM(0,0,0,0);
                if(fSamples[i].sname=="Wtolnu"||namesample.Contains("WJets")||namesample.Contains("W-")||nominal=="WJets"){
                    bool foundW(false);
                    for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
                        int ID  =abs(fMT2tree->genlept[ngl].ID);
                        int MID =abs(fMT2tree->genlept[ngl].MID);
                        if((ID == 11 || ID == 12 || ID == 13 || ID == 14 || ID == 15 || ID == 16) && MID==24){
                            hardgenlv = fMT2tree->genlept[ngl].Mlv; foundW = true;
                        }
                        if(foundW) break;
                    }
                    if(!foundW){
                        for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
                            int ID  =abs(fMT2tree->genlept[ngl].ID);
                            int MID =abs(fMT2tree->genlept[ngl].MID);
                            int GMID=abs(fMT2tree->genlept[ngl].GMID);
                            if((ID == 11 || ID == 12 || ID == 13 || ID == 14 || ID == 16) && MID==15 && GMID==24){
                                hardgenlv = fMT2tree->genlept[ngl].Mlv; foundW = true;
                            }
                            if(foundW) break;
                        }
                    }
                } if(namesample.Contains("ZJets")){
                    hardgenlv = fMT2tree->GenZ[0];
                } if(fSamples[i].sname=="Top"||namesample.Contains("TTbar")||namesample.Contains("TT-")||nominal=="TTbar"){
                    TLorentzVector top1(0.,0.,0.,0.), top2(0.,0.,0.,0.);
                    bool top1f(false), top2f(false);
                    for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
                        int id   = abs(fMT2tree->genlept[ngl].ID);
                        if(id!=5) continue;
                        int mid  = fMT2tree->genlept[ngl].MID;//from b
                        if(mid==6&&top1f) continue;
                        else if(mid==6) { top1 = fMT2tree->genlept[ngl].Mlv; top1f = true; }
                        if(mid==-6&&top2f) continue;
                        else if(mid==-6){ top2 = fMT2tree->genlept[ngl].Mlv; top2f = true; }
                        if(top1f&&top2f) {
                            hardgenlv = top1+top2;
                            break;
                        }
                    }
                } if(namesample.Contains("T_t")||namesample.Contains("T_s")||namesample.Contains("Tbar_t")||namesample.Contains("Tbar_s")||namesample.Contains("T-t")||namesample.Contains("T-s")||namesample.Contains("Tbar-t")||namesample.Contains("Tbar-s")||nominal=="SingleTop"){
                    if((fSamples[i].name).Contains("tW")){//t + W
                        TLorentzVector top(0.,0.,0.,0.), W(0.,0.,0.,0.);
                        bool topf(false), Wf(false);
                        for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
                            int id    = abs(fMT2tree->genlept[ngl].ID);
                            int mid   = abs(fMT2tree->genlept[ngl].MID);//from b
                            int gmid  = abs(fMT2tree->genlept[ngl].GMID);//from b
                            if(mid==6&&topf) continue;
                            else if(mid==6) { top = fMT2tree->genlept[ngl].Mlv; topf = true; }
                            if(mid==24&&gmid!=6&&Wf) continue;
                            if((id == 11 || id == 12 || id == 13 || id == 14 || id == 15 || id == 16) && mid==24 && gmid!=6){
                                W = fMT2tree->genlept[ngl].Mlv; Wf = true;
                            }
                            if(topf&&Wf){
                                hardgenlv = top+W;
                                break;
                            }
                        }
                        if(!Wf){//this might be wrong
                            for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
                                int id    = abs(fMT2tree->genlept[ngl].ID);
                                int mid   = abs(fMT2tree->genlept[ngl].MID);//from b
                                int gmid  = abs(fMT2tree->genlept[ngl].GMID);//from b
                                if(mid==6||gmid==6) continue;
                                if(gmid==24&&gmid==15&&Wf) continue;
                                if((id == 11 || id == 12 || id == 13 || id == 14 || id == 15 || id == 16) && mid==15 && gmid==24){
                                    W = fMT2tree->genlept[ngl].Mlv; Wf = true;
                                }
                                if(topf&&Wf){
                                    hardgenlv = top+W;
                                    break;
                                }
                            }
                        }
                        
                    } else {
                        TLorentzVector top(0.,0.,0.,0.);
                        bool topf(false), Wf(false);
                        for(int ngl = 0; ngl<fMT2tree->NGenLepts; ++ngl){
                            int id    = abs(fMT2tree->genlept[ngl].ID);
                            int mid   = abs(fMT2tree->genlept[ngl].MID);//from b
                            int gmid  = abs(fMT2tree->genlept[ngl].GMID);//from b
                            if(mid==6&&topf) continue;
                            else if(mid==6) { top = fMT2tree->genlept[ngl].Mlv; topf = true; }
                            if(topf){
                                hardgenlv = top;
                                break;
                            }
                        }
                    }
                }
                if(hardgenlv.Pt()>250.) ISRweight = 0.8;
                else if(hardgenlv.Pt()>150.) ISRweight = 0.9;
                else if(hardgenlv.Pt()>120.) ISRweight = 0.95;
                else                         ISRweight = 1.;
                if(fISRreweight_ weight = weight * ISRweight;
            }

		double btagSF(1.), btagSFerr(0.);
		if(!fMT2tree->misc.isData && fbTagReweight){
			btagSF = fMT2tree->SFWeight.BTagCSV40ge2; btagSFerr = fMT2tree->SFWeight.BTagCSV40ge2Error;
			weight = weight * btagSF;
			//if(!fbTagError) btagSFerr = 0;//don't store btag error is equivalent to no error at all
		}
		int HTbin(-1); string sHT = "";
		if(HTval>=450.  && HTval<750.  /*&& METval>200.*/) { HTbin = 0; sHT = "_lowHT";    }
		if(HTval>=750.                 && METval>30. ) { HTbin = 1; sHT = "_highHT";   }
		string ssignal = ""; int signalbin = 1;

		string hh = string("_") + namesample.Data() + string("_") + kindsample.Data() + sHT;
		string helptype = " ";

		//nominal TTbar
		if(signalbin>=0 && HTbin>=0 && VSPT<70. && J1Pt>100. && J2Pt>100. && MinDPhi>0.3/* && MT2val>100. && ((HTbin==0&&MT2val>200.)||(HTbin==1&&MT2val>125.))*/ ){
		//can put lower cut on MT2 here, as no event modification is made here; only weights change
		//nominal shape, xsec up/down for each process, ISR, TopMass, Matching, Scale
		if(nominal=="TTbar"){
			if(namesample=="TTbar_Madgraph0l") helptype = "_nominal0l";
			if(namesample=="TTbar_Madgraph1l") helptype = "_nominal1l";
			if(namesample=="TTbar_Madgraph2l") helptype = "_nominal2l";
			hh = "_TTbar"+sHT+ssignal;
			if (histos.count("MT2"+helptype     +hh)==0) cout << "MT2"+helptype     +hh << endl;

			histos["MT2"+helptype     +hh]->Fill(Mbbval, weight);
			histos["MT2_Nominal"      +hh]->Fill(Mbbval, weight);
			if((HTbin==0&&MT2val>75. )||(HTbin==1&&MT2val>50. )) histos["MT2_MT2relaxedcut" +hh]->Fill(Mbbval, weight);
			if((HTbin==0&&MT2val>200.)||(HTbin==1&&MT2val>125.)) histos["MT2_MT2fullcut"    +hh]->Fill(Mbbval, weight);
			histos["MT2_TopUp"        +hh]->Fill(Mbbval, weight*1.07);
			histos["MT2_TopDown"      +hh]->Fill(Mbbval, weight*0.93);
			histos["MT2_SingleTopUp"  +hh]->Fill(Mbbval, weight);
			histos["MT2_SingleTopDown"+hh]->Fill(Mbbval, weight);
			histos["MT2_WUp"          +hh]->Fill(Mbbval, weight);
			histos["MT2_WDown"        +hh]->Fill(Mbbval, weight);
			histos["MT2_TTVUp"        +hh]->Fill(Mbbval, weight);
			histos["MT2_TTVDown"      +hh]->Fill(Mbbval, weight);

			int t1weight(1), t2weight(1);
			bool t1(false), t2(false);
			for(int n = 0; n<20; ++n){
				int mid  = fMT2tree->genlept[n].MID;
				if(mid==6&&t1) continue;
				else if(mid==6){
					double t1pt = fMT2tree->genlept[n].Mlv.Pt();
					if(t1pt>200.&&t1pt<600.) t1weight = 1.-(t1pt-200.)/800.;
					if(t1pt>600.)            t1weight = 0.25;
					t1 = true;
				}
				if(mid==-6&&t2) continue;
				else if(mid==-6){
					double t2pt = fMT2tree->genlept[n].Mlv.Pt();
					if(t2pt>200.&&t2pt<600.) t2weight = 1.-(t2pt-200.)/800.;
					if(t2pt>600.)            t2weight = 0.25;
					t2 = true;
				}
				if(t1&&t2) break;
			}
			histos["MT2_TopPtReweighted"+hh]->Fill(Mbbval, weight*t1weight*t2weight);
			histos["MT2_ISR"+hh]->Fill(Mbbval, weightnoISR*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(Mbbval, weightnoISR);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(Mbbval, weightnoISR*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
		}
		//nominal W
		if(nominal=="WJets"){
			hh = "_WJets"+sHT+ssignal;
			histos["MT2_Nominal"      +hh]->Fill(Mbbval, weight);
			if((HTbin==0&&MT2val>75. )||(HTbin==1&&MT2val>50. )) histos["MT2_MT2relaxedcut" +hh]->Fill(Mbbval, weight);
			if((HTbin==0&&MT2val>200.)||(HTbin==1&&MT2val>125.)) histos["MT2_MT2fullcut"    +hh]->Fill(Mbbval, weight);
			histos["MT2_WUp"          +hh]->Fill(Mbbval, weight*1.5);
			histos["MT2_WDown"        +hh]->Fill(Mbbval, weight*0.5);
			histos["MT2_TopUp"        +hh]->Fill(Mbbval, weight);
			histos["MT2_TopDown"      +hh]->Fill(Mbbval, weight);
			histos["MT2_SingleTopUp"  +hh]->Fill(Mbbval, weight);
			histos["MT2_SingleTopDown"+hh]->Fill(Mbbval, weight);
			histos["MT2_TTVUp"        +hh]->Fill(Mbbval, weight);
			histos["MT2_TTVDown"      +hh]->Fill(Mbbval, weight);
			histos["MT2_MassUp"       +hh]->Fill(Mbbval, weight);
			histos["MT2_MassDown"     +hh]->Fill(Mbbval, weight);
			histos["MT2_ISR"+hh]->Fill(Mbbval, weightnoISR*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(Mbbval, weightnoISR);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(Mbbval, weightnoISR*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
		}
		//nominal singleTop
		if(nominal=="SingleTop"){
			hh = "_SingleTop"+sHT+ssignal;
			histos["MT2_Nominal"      +hh]->Fill(Mbbval, weight);
			if((HTbin==0&&MT2val>75. &&METval>200.)||(HTbin==1&&MT2val>50. )) histos["MT2_MT2relaxedcut" +hh]->Fill(Mbbval, weight);
			if((HTbin==0&&MT2val>200.&&METval>200.)||(HTbin==1&&MT2val>125.)) histos["MT2_MT2fullcut"    +hh]->Fill(Mbbval, weight);
			histos["MT2_WUp"          +hh]->Fill(Mbbval, weight);
			histos["MT2_WDown"        +hh]->Fill(Mbbval, weight);
			histos["MT2_TopUp"        +hh]->Fill(Mbbval, weight);
			histos["MT2_TopDown"      +hh]->Fill(Mbbval, weight);
			histos["MT2_SingleTopUp"  +hh]->Fill(Mbbval, weight*1.1);
			histos["MT2_SingleTopDown"+hh]->Fill(Mbbval, weight*0.9);
			histos["MT2_TTVUp"        +hh]->Fill(Mbbval, weight);
			histos["MT2_TTVDown"      +hh]->Fill(Mbbval, weight);
			histos["MT2_MassUp"       +hh]->Fill(Mbbval, weight);
			histos["MT2_MassDown"     +hh]->Fill(Mbbval, weight);
			histos["MT2_MatchingUp"   +hh]->Fill(Mbbval, weight);
			histos["MT2_MatchingDown" +hh]->Fill(Mbbval, weight);
			histos["MT2_ScaleUp"      +hh]->Fill(Mbbval, weight);
			histos["MT2_ScaleDown"    +hh]->Fill(Mbbval, weight);
			histos["MT2_ISR"+hh]->Fill(Mbbval, weightnoISR*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(Mbbval, weightnoISR);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(Mbbval, weightnoISR*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)

		}
		//nominal rest
		if(nominal=="TTV"){
			hh = "_TTV"+sHT+ssignal;
			histos["MT2_Nominal"      +hh]->Fill(Mbbval, weight);
			if((HTbin==0&&MT2val>75. )||(HTbin==1&&MT2val>50. )) histos["MT2_MT2relaxedcut" +hh]->Fill(Mbbval, weight);
			if((HTbin==0&&MT2val>200.)||(HTbin==1&&MT2val>125.)) histos["MT2_MT2fullcut"    +hh]->Fill(Mbbval, weight);
			histos["MT2_WUp"          +hh]->Fill(Mbbval, weight);
			histos["MT2_WDown"        +hh]->Fill(Mbbval, weight);
			histos["MT2_TopUp"        +hh]->Fill(Mbbval, weight);
			histos["MT2_TopDown"      +hh]->Fill(Mbbval, weight);
			histos["MT2_SingleTopUp"  +hh]->Fill(Mbbval, weight);
			histos["MT2_SingleTopDown"+hh]->Fill(Mbbval, weight);
			histos["MT2_TTVUp"        +hh]->Fill(Mbbval, weight*1.3);
			histos["MT2_TTVDown"      +hh]->Fill(Mbbval, weight*0.7);
			histos["MT2_MassUp"       +hh]->Fill(Mbbval, weight);
			histos["MT2_MassDown"     +hh]->Fill(Mbbval, weight);
			histos["MT2_MatchingUp"   +hh]->Fill(Mbbval, weight);
			histos["MT2_MatchingDown" +hh]->Fill(Mbbval, weight);
			histos["MT2_ScaleUp"      +hh]->Fill(Mbbval, weight);
			histos["MT2_ScaleDown"    +hh]->Fill(Mbbval, weight);
			histos["MT2_ISR"+hh]->Fill(Mbbval, weightnoISR*ISRweight);//normalization done by Nominal
			histos["MT2_ISRUp"+hh]->Fill(Mbbval, weightnoISR);//ISRUp is Nominal
			histos["MT2_ISRDown"+hh]->Fill(Mbbval, weightnoISR*(2.*ISRweight - 1.));//1.-2.*(1-ISRweight)
		}
		//scale up/down
		if(namesample=="TT-scaleup"||namesample=="TT-scaledown"||namesample=="W-scaleup"||namesample=="W-scaledown"){
			hh = sHT+ssignal;
			if(namesample=="TT-scaleup"  ) histos["MT2_ScaleUp_TTbar"  +hh]->Fill(Mbbval,weight);
			if(namesample=="TT-scaledown") histos["MT2_ScaleDown_TTbar"+hh]->Fill(Mbbval,weight);
			if(!WnoScaleMatching && namesample== "W-scaleup"  ) histos["MT2_ScaleUp_WJets"  +hh]->Fill(Mbbval,weight);
			if(!WnoScaleMatching && namesample== "W-scaledown") histos["MT2_ScaleDown_WJets"+hh]->Fill(Mbbval,weight);
		}
		//mass up/down
		if(namesample=="TT-massup"||namesample=="TT-massdown"){
			hh = "_TTbar"+sHT+ssignal;
			if(namesample=="TT-massup"  ) histos["MT2_MassUp"  +hh]->Fill(Mbbval,weight);
			if(namesample=="TT-massdown") histos["MT2_MassDown"+hh]->Fill(Mbbval,weight);
		}
		//matching up/down
		if(namesample=="TT-matchingup"||namesample=="TT-matchingdown"||namesample=="W-matchingup"||namesample=="W-matchingdown"){
			hh = sHT+ssignal;
			if(namesample=="TT-matchingup"  ) histos["MT2_MatchingUp_TTbar"  +hh]->Fill(Mbbval,weight);
			if(namesample=="TT-matchingdown") histos["MT2_MatchingDown_TTbar"+hh]->Fill(Mbbval,weight);
			if(!WnoScaleMatching && namesample== "W-matchingup"  ) histos["MT2_MatchingUp_WJets"  +hh]->Fill(Mbbval,weight);
			if(!WnoScaleMatching && namesample== "W-matchingdown") histos["MT2_MatchingDown_WJets"+hh]->Fill(Mbbval,weight);
		}
		if(namesample=="WJets_Incl"){//moved up
			hh = "_WJets"+sHT+ssignal;
			if(!WnoScaleMatching) histos["MT2_nominal1Sample"+hh]->Fill(Mbbval, weight);
		}
		if(nominal=="WJets" && WnoScaleMatching){//only if not using scale up/down for WJets above due to too small statistics
			hh = sHT+ssignal;
			histos["MT2_ScaleUp_WJets"  +hh]->Fill(Mbbval,weight);
			histos["MT2_ScaleDown_WJets"+hh]->Fill(Mbbval,weight);
			histos["MT2_MatchingUp_WJets"  +hh]->Fill(Mbbval,weight);
			histos["MT2_MatchingDown_WJets"+hh]->Fill(Mbbval,weight);
			histos["MT2_nominal1Sample_WJets"+hh]->Fill(Mbbval, weight);
		}
		//JES up/down
		if(namesample=="TTbar_Madgraph0l_JESdown" || namesample=="TTbar_Madgraph1l_JESdown" || namesample=="TTbar_Madgraph2l_JESdown") helptype = "TTbar_JESdown";
		if(namesample=="TTbar_Madgraph0l_JESup" || namesample=="TTbar_Madgraph1l_JESup" || namesample=="TTbar_Madgraph2l_JESup") helptype = "TTbar_JESup";
		if(namesample=="WJetsHT200_JESdown"||namesample=="WJetsHT250_JESdown"||namesample=="WJetsHT300_JESdown"||namesample=="WJetsHT400_JESdown") helptype = "WJets_JESdown";
		if(namesample=="WJetsHT200_JESup"||namesample=="WJetsHT250_JESup"||namesample=="WJetsHT300_JESup"||namesample=="WJetsHT400_JESup") helptype = "WJets_JESup";
		if(namesample=="T-s_JESdown"||namesample=="Tbar-s_JESdown"||namesample=="T-t_JESdown"||namesample=="Tbar-t_JESdown"||namesample=="T-tW_JESdown"||namesample=="Tbar-tW_JESdown") helptype = "SingleTop_JESdown";
		if(namesample=="T-s_JESup"||namesample=="Tbar-s_JESup"||namesample=="T-t_JESup"||namesample=="Tbar-t_JESup"||namesample=="T-tW_JESup"||namesample=="Tbar-tW_JESup") helptype = "SingleTop_JESup";
		if(namesample=="WW_JESdown"||namesample=="WZ_JESdown"||namesample=="ZZ_JESdown"||namesample=="TTW_JESdown"||namesample=="TTZ_JESdown"||namesample=="TTG_JESdown") helptype = "TTV_JESdown";
		if(namesample=="WW_JESup"||namesample=="WZ_JESup"||namesample=="ZZ_JESup"||namesample=="TTW_JESup"||namesample=="TTZ_JESup"||namesample=="TTG_JESup") helptype = "TTV_JESup";

		if(helptype=="TTbar_JESdown"||helptype=="TTbar_JESup"||helptype=="WJets_JESdown"||helptype=="WJets_JESup"||helptype=="SingleTop_JESdown"||helptype=="SingleTop_JESup"||helptype=="TTV_JESdown"||helptype=="TTV_JESup"){
			hh = sHT+ssignal;
			if(helptype=="TTbar_JESdown")     histos["MT2_JESDown_TTbar"    +hh]->Fill(Mbbval,weight);
			if(helptype=="TTbar_JESup"  )     histos["MT2_JESUp_TTbar"      +hh]->Fill(Mbbval,weight);
			if(helptype=="WJets_JESdown")     histos["MT2_JESDown_WJets"    +hh]->Fill(Mbbval,weight);
			if(helptype=="WJets_JESup"  )     histos["MT2_JESUp_WJets"      +hh]->Fill(Mbbval,weight);
			if(helptype=="SingleTop_JESdown") histos["MT2_JESDown_SingleTop"+hh]->Fill(Mbbval,weight);
			if(helptype=="SingleTop_JESup"  ) histos["MT2_JESUp_SingleTop"  +hh]->Fill(Mbbval,weight);
			if(helptype=="TTV_JESdown")       histos["MT2_JESDown_TTV"      +hh]->Fill(Mbbval,weight);
			if(helptype=="TTV_JESup"  )       histos["MT2_JESUp_TTV"        +hh]->Fill(Mbbval,weight);
		}
		//BTV, PU, different TTbar MC
		if(nominal=="TTbar"||nominal=="WJets"||nominal=="SingleTop"||nominal=="TTV"){
			hh = "_"+nominal+sHT+ssignal;
			//BTV
			if(btagSFerr>0 && btagSF!=0){
			histos["MT2_BSFUp"         +hh]->Fill(Mbbval, weight*(btagSF+btagSFerr)/btagSF);
			histos["MT2_BSFDown"       +hh]->Fill(Mbbval, weight*(btagSF-btagSFerr)/btagSF);
			} else {
			histos["MT2_BSFUp"         +hh]->Fill(Mbbval, weight);
			histos["MT2_BSFDown"       +hh]->Fill(Mbbval, weight);
			}
			//PUWeight
			int pubin = pu_up->FindBin(fMT2tree->pileUp.PUtrueNumInt);
			double weightup = sample_weight * pu_up->GetBinContent(pubin);
			double weightdown = sample_weight * pu_down->GetBinContent(pubin);
			histos["MT2_PUUp"          +hh]->Fill(Mbbval, weightup);
			histos["MT2_PUDown"        +hh]->Fill(Mbbval, weightdown);
			if(nominal!="TTbar"){
			if(nominal!="WJets") histos["MT2_nominal1Sample"+hh]->Fill(Mbbval, weight);
			histos["MT2_nominalPowheg" +hh]->Fill(Mbbval, weight);
			}
		}
		if(namesample=="TTbar_MadgraphIncl"){
			hh = "_TTbar"+sHT+ssignal;
			histos["MT2_nominal1Sample"+hh]->Fill(Mbbval, weight);
		}
		if(namesample=="TTbar_Powheg"){
			hh = "_TTbar"+sHT+ssignal;
			histos["MT2_nominalPowheg" +hh]->Fill(Mbbval, weight);
		}
		}//if(signalbin>=0 && HTbin>=0 && VSPT<70. && J1Pt>100. && J2Pt>100. && MinDPhi<0.3)
		//unclustered MET
		if(nominal=="TTbar"||nominal=="WJets"||nominal=="SingleTop"||nominal=="TTV"){
			hh = "_"+nominal+sHT+ssignal;
			//unclustered MET
			if(signalbin>=0 && J1Pt>100. && J2Pt>100.&& MinDPhi>0.3){
				//no angular uncertainty
				TLorentzVector METunclusteredup = fMT2tree->pfmet[0] - fMT2tree->MHT[0];
				METunclusteredup.SetPtEtaPhiM(1.1*METunclusteredup.Pt(), METunclusteredup.Eta(), METunclusteredup.Phi(), METunclusteredup.M());
				if(METunclusteredup.Pt()<70){//VSPT cut
					METunclusteredup += fMT2tree->MHT[0];
					int uHTbin(-1); string usHT;
					if(HTval>=450.  && HTval<750.  /*&& METunclusteredup.Pt()>200.*/) { uHTbin = 0; usHT = "_lowHT";    }
					if(HTval>=750.                 && METunclusteredup.Pt()>30. ) { uHTbin = 1; usHT = "_highHT";   }
					if(uHTbin>=0){
						hh = "_"+nominal+usHT+ssignal;
						histos["MT2_METUp"+hh]->Fill(Mbbval, weight);
					}
				}
				TLorentzVector METunclustereddown = fMT2tree->pfmet[0] - fMT2tree->MHT[0];
				METunclustereddown.SetPtEtaPhiM(0.9*METunclustereddown.Pt(), METunclustereddown.Eta(), METunclustereddown.Phi(), METunclustereddown.M());
				if(METunclustereddown.Pt()<70){//VSPT cut
					METunclustereddown += fMT2tree->MHT[0];
					int dHTbin(-1); string dsHT;
					if(HTval>=450.  && HTval<750.  /*&& METunclustereddown.Pt()>200.*/) { dHTbin = 0; dsHT = "_lowHT";    }
					if(HTval>=750.                 && METunclustereddown.Pt()>30. ) { dHTbin = 1; dsHT = "_highHT";   }
					if(dHTbin>=0){
						hh = "_"+nominal+dsHT+ssignal;
						histos["MT2_METDown"+hh]->Fill(Mbbval, weight);
					}
				}
			}
		}
		//all direct changes

 	}//while(myEvtList->GetEntry(counter++) !=-1)
	delete fMT2tree;
	delete fSamples[i].tree;
	}//for(size_t i = 0; i < fSamples.size(); ++i)
	cout << "add overflow to last bin" << endl;
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->SetBinContent(h->second->GetNbinsX(),
					 h->second->GetBinContent(h->second->GetNbinsX()  )+ 
					 h->second->GetBinContent(h->second->GetNbinsX()+1)  );
		h->second->SetBinError(  h->second->GetNbinsX(),
					 sqrt(h->second->GetBinError(h->second->GetNbinsX()  )*
					      h->second->GetBinError(h->second->GetNbinsX()  )+
					      h->second->GetBinError(h->second->GetNbinsX()+1)*
					      h->second->GetBinError(h->second->GetNbinsX()+1)  ));
		h->second->SetBinContent(0, h->second->GetBinContent(0 )+  h->second->GetBinContent(1)  );
		h->second->SetBinError(  0, sqrt(pow(h->second->GetBinError(0  ),2)+pow(h->second->GetBinError(1),2)));
	}
	cout << "add samples to allMC" << endl;
	for(int i1 = 0; i1<sampletypesize;   ++i1){
	for(int i3 = 0; i3<HTregionsize;     ++i3){
		string ht = string("MT2_") + sample_type[i1];
		string hs = string("_") + HT_region[i3];// + string("_") + signal_region[i2];
		string mapnameallMC   = ht + string("_allMC")     + hs;
		string mapnameTop     = ht + string("_allTop")    + hs;
		string mapnameTTbar   = ht + string("_TTbar")     + hs;
		string mapnameSingleT = ht + string("_SingleTop") + hs;
		string mapnameTTV     = ht + string("_TTV")       + hs;
		string mapnameWJets   = ht + string("_WJets")     + hs;
		histos[mapnameTop]  ->Add(histos[mapnameTTbar]);
		histos[mapnameTop]  ->Add(histos[mapnameSingleT]);
		histos[mapnameallMC]->Add(histos[mapnameTTbar]);
		histos[mapnameallMC]->Add(histos[mapnameSingleT]);
		histos[mapnameallMC]->Add(histos[mapnameTTV]);
		histos[mapnameallMC]->Add(histos[mapnameWJets]);
	}}
	cout << "create normalized distribution" << endl;
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		string mapname    = h->first;
		string mapnamenorm = "Norm"+mapname;
		TH1D *hcopy = (TH1D*)h->second->Clone(mapnamenorm.c_str());
		if(hcopy->Integral()>0) hcopy->Scale(1./hcopy->Integral());
		if(histosnormalized.count(mapnamenorm) == 0 ) histosnormalized[mapnamenorm] = hcopy;
	}
	cout << "Saving." << endl;
    	TFile *fsavefile = new TFile(outputdir + outputname,"RECREATE");
	fsavefile->cd();
	for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Write();
	}
	for(map<string,TH1D*>::iterator h=histosnormalized.begin(); h!=histosnormalized.end();++h){
		h->second->Write();
	}
	fsavefile->Close();
	cout << "Saved histograms in " << outputdir << outputname << endl;

	//plotting is done in TTbarStudiesHiggs_MakePlots.C

}//void TTbarStudiesHiggs



void load(const char* filename){
	
	fSamples.clear();
	char buffer[200];
	ifstream IN(filename);

	char ParName[100], StringValue[1000];
	float ParValue;

	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Sample File  " << filename << endl;
	int counter(0);
	
	while( IN.getline(buffer, 200, '\n') ){
		// ok = false;
		if (buffer[0] == '#') {
			continue; // Skip lines commented with '#'
		}
		if( !strcmp(buffer, "GENERAL") ){
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Path\t%s", StringValue);
			fPath = StringValue;	
			cout << fPath << endl;
			
			if(fVerbose >0){
				cout << " ----  " << endl;
				cout << "  Path " << fPath << endl;
			}

		}
		if( !strcmp(buffer, "SAMPLE")){

			sample s;
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Name\t%s", StringValue);
			s.name = TString(StringValue);

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "SName\t%s", StringValue);
			s.sname = TString(StringValue);

			IN.getline(buffer, 200, '\n');			//new
			sscanf(buffer, "ShapeName\t%s", StringValue);	//new
			s.shapename = TString(StringValue);		//new

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "File\t%s", StringValue);
			TString file =fPath+StringValue;
			TFile *f = TFile::Open(file);
			s.file = f;
			s.tree = (TTree*)f->Get("MassTree");
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Xsection\t%f", &ParValue);
			s.xsection = ParValue;
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Kfact\t%f", &ParValue);
			s.kfact = ParValue;
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Lumi\t%f", &ParValue);
			s.lumi = ParValue;

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Type\t%s", StringValue);
			s.type = StringValue;
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Color\t%f", &ParValue);
			s.color = ParValue;

			if(s.type!="data"){
			TH1F *h_PUWeights = (TH1F*) s.file->Get("h_PUWeights");
			TH1F *h_Events    = (TH1F*) s.file->Get("h_Events");
			if(h_PUWeights==0 || h_Events==0){
				cout << "ERROR: sample " << (s.file)->GetName() << " does not have PU and NEvents histos! " << endl;
				exit(1);
			}
			s.type!="data" ? s.PU_avg_weight = h_PUWeights->GetMean()    : s.PU_avg_weight =1;
			s.type!="data" ? s.nevents       = h_Events   ->GetEntries() : s.nevents       =1;
			delete h_PUWeights;
			delete h_Events;
			} else {
			s.PU_avg_weight =1;
			s.nevents       =1;
			}
			if(fVerbose > 0){
				cout << " ---- " << endl;
				cout << "  New sample added: " << s.name << endl;
				cout << "   Sample no.      " << counter << endl;
				cout << "   Short name:     " << s.sname << endl;
				cout << "   File:           " << (s.file)->GetName() << endl;
				cout << "   Events:         " << s.nevents  << endl;
				cout << "   Events in tree: " << s.tree->GetEntries() << endl; 
				cout << "   Xsection:       " << s.xsection << endl;
				cout << "   Lumi:           " << s.lumi << endl;
				cout << "   kfactor:        " << s.kfact << endl;
				cout << "   avg PU weight:  " << s.PU_avg_weight << endl;
				cout << "   type:           " << s.type << endl;
				cout << "   Color:          " << s.color << endl;
			}
			fSamples.push_back(s);
			counter++;
		}
	}
	if(fVerbose > 0) cout << "------------------------------------" << endl;
}
