{//run via root -l testStitching.C

//checks - using partonicHT - if stitching works (basically if partonic HT distribution is smooth of if some cross section is wrong)
//for WJets,Znunu,QCD,GJets Madgraph - uncomment accordingly
//cross sections are hardcoded


gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
/*
//WJets
TFile *_file0 = TFile::Open("/shome/casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/WJetsToLNu-HT-200To250-8TeV-madgraph-Summer12-DR53X-PU-S10-START53-V7C-v1-2.root");
_file0->cd();
TTree *tree0 = (TTree*)_file0->Get("MassTree");
TFile *_file1 = TFile::Open("/shome/casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/WJetsToLNu-HT-250To300-8TeV-madgraph-v2-Summer12-DR53X-PU-S10-START53-V7A-v1.root");
_file1->cd();
TTree *tree1 = (TTree*)_file1->Get("MassTree");
TFile *_file2 = TFile::Open("/shome/casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/WJetsToLNu-HT-300To400-8TeV-madgraph-v2-Summer12-DR53X-PU-S10-START53-V7A-v1.root");
_file2->cd();
TTree *tree2 = (TTree*)_file2->Get("MassTree");
TFile *_file3 = TFile::Open("/shome/casal/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/WJetsToLNu-HT-400ToInf-8TeV-madgraph-v2-Summer12-DR53X-PU-S10-START53-V7A-v1.root");
_file3->cd();
TTree *tree3 = (TTree*)_file3->Get("MassTree");
*/
/*
//WJets
TChain *tree0 = new TChain("MassTree");
tree0->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/WJetsToLNu-HT-200To250-8TeV-madgraph-Summer12-DR53X-PU-S10-START53-V7C-v1-2/output_*.root");
TChain *tree1 = new TChain("MassTree");
tree1->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/WJetsToLNu-HT-250To300-8TeV-madgraph-v2-Summer12-DR53X-PU-S10-START53-V7A-v1/output_*.root");
TChain *tree2 = new TChain("MassTree");
tree2->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/WJetsToLNu-HT-300To400-8TeV-madgraph-v2-Summer12-DR53X-PU-S10-START53-V7A-v1/output_*.root");
TChain *tree3 = new TChain("MassTree");
tree3->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/mmasciov/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/WJetsToLNu-HT-400ToInf-8TeV-madgraph-v2-Summer12-DR53X-PU-S10-START53-V7A-v1/output_*.root");
*/
/*
//Znunu
TFile *_file0 = TFile::Open("/shome/haweber/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/ZJetsToNuNu-50-HT-100-TuneZ2Star-8TeV-madgraph-ext-Summer12-DR53X-PU-S10-START53-V7A-v1.root");
_file0->cd();
TTree *tree0 = (TTree*)_file0->Get("MassTree");
TFile *_file1 = TFile::Open("/shome/haweber/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/ZJetsToNuNu-100-HT-200-TuneZ2Star-8TeV-madgraph-ext-Summer12-DR53X-PU-S10-START53-V7C-v1.root");
_file1->cd();
TTree *tree1 = (TTree*)_file1->Get("MassTree");
TFile *_file2 = TFile::Open("/shome/haweber/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/ZJetsToNuNu-200-HT-400-TuneZ2Star-8TeV-madgraph-ext-Summer12-DR53X-PU-S10-START53-V7A-v1.root");
_file2->cd();
TTree *tree2 = (TTree*)_file2->Get("MassTree");
TFile *_file3 = TFile::Open("/shome/haweber/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/ZJetsToNuNu-400-HT-inf-TuneZ2Star-8TeV-madgraph-ext-Summer12-DR53X-PU-S10-START53-V7A-v1.root");
_file3->cd();
TTree *tree3 = (TTree*)_file3->Get("MassTree");
*/
/*
TChain *tree0 = new TChain("MassTree");
tree0->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/ZJetsToNuNu-50-HT-100-TuneZ2Star-8TeV-madgraph-ext-Summer12-DR53X-PU-S10-START53-V7A-v1/output_*.root");
TChain *tree1 = new TChain("MassTree");
tree1->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/ZJetsToNuNu-100-HT-200-TuneZ2Star-8TeV-madgraph-ext-Summer12-DR53X-PU-S10-START53-V7C-v1/output_*.root");
TChain *tree2 = new TChain("MassTree");
tree2->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/ZJetsToNuNu-200-HT-400-TuneZ2Star-8TeV-madgraph-ext-Summer12-DR53X-PU-S10-START53-V7A-v1/output_*.root");
TChain *tree3 = new TChain("MassTree");
tree3->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/ZJetsToNuNu-400-HT-inf-TuneZ2Star-8TeV-madgraph-ext-Summer12-DR53X-PU-S10-START53-V7A-v1/output_*.root");
*/

//QCD
/*
TFile *_file0 = TFile::Open("/shome/haweber/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/QCD-HT-100To250-TuneZ2star-8TeV-madgraph-pythia-Summer12-DR53X-PU-S10-START53-V7A-v1.root");
_file0->cd();
TTree *tree0 = (TTree*)_file0->Get("MassTree");
TFile *_file1 = TFile::Open("/shome/haweber/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/QCD-HT-250To500-TuneZ2star-8TeV-madgraph-pythia6-Summer12-DR53X-PU-S10-START53-V7A-v1-2.root");
_file1->cd();
TTree *tree1 = (TTree*)_file1->Get("MassTree");
TFile *_file2 = TFile::Open("/shome/haweber/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/QCD-HT-500To1000-TuneZ2star-8TeV-madgraph-pythia6-Summer12-DR53X-PU-S10-START53-V7A-v1-2.root");
_file2->cd();
TTree *tree2 = (TTree*)_file2->Get("MassTree");
TFile *_file3 = TFile::Open("/shome/haweber/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/QCD-HT-1000ToInf-TuneZ2star-8TeV-madgraph-pythia6-Summer12-DR53X-PU-S10-START53-V7A-v1-2.root");
_file3->cd();
TTree *tree3 = (TTree*)_file3->Get("MassTree");
*/
/*
TChain *tree0 = new TChain("MassTree");
tree0->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/QCD-HT-100To250-TuneZ2star-8TeV-madgraph-pythia-Summer12-DR53X-PU-S10-START53-V7A-v1/output_*.root");
TChain *tree1 = new TChain("MassTree");
tree1->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/QCD-HT-250To500-TuneZ2star-8TeV-madgraph-pythia6-Summer12-DR53X-PU-S10-START53-V7A-v1-2/output_*.root");
TChain *tree2 = new TChain("MassTree");
tree2->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/QCD-HT-500To1000-TuneZ2star-8TeV-madgraph-pythia6-Summer12-DR53X-PU-S10-START53-V7A-v1-2/output_*.root");
TChain *tree3 = new TChain("MassTree");
tree3->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/casal/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV/QCD-HT-1000ToInf-TuneZ2star-8TeV-madgraph-pythia6-Summer12-DR53X-PU-S10-START53-V7A-v1-2/output_*.root");
*/

//GJets
TChain *tree0 = new TChain("MassTree");
tree0->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV_1g_Gremoved/GJets-HT-200To400-8TeV-madgraph-v2-Summer12-DR53X-PU-S10-START53-V7A-v1-2ndShot/output_*.root");
TChain *tree1 = new TChain("MassTree");
tree1->Add("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/SUSY/MassTrees/MT2_V02-03-02/20130318_8TeV_1g_Gremoved/GJets-HT-400ToInf-8TeV-madgraph-v3-Summer12-DR53X-PU-S10-START53-V7C-v1-2ndShot/output_*.root");

TH1D *h0_HT = new TH1D("h0_HT", "", 100, 0., 1000.); h0_HT->Sumw2(); h0_HT->SetLineColor(kBlue); h0_HT->SetMinimum(0.002); 
TH1D *h1_HT = new TH1D("h1_HT", "", 100, 0., 1000.); h1_HT->Sumw2(); h1_HT->SetLineColor(kRed); h1_HT->SetMinimum(0.002);
TH1D *h2_HT = new TH1D("h2_HT", "", 100, 0., 1000.); h2_HT->Sumw2(); h2_HT->SetLineColor(kCyan); h2_HT->SetMinimum(0.002);
TH1D *h3_HT = new TH1D("h3_HT", "", 100, 0., 1000.); h3_HT->Sumw2(); h3_HT->SetLineColor(kBlack); h3_HT->SetMinimum(0.002);

//WJets
/*
tree0->Draw("misc.QCDPartonicHT>>h0_HT","(0.1800) * (NJets>0)","goff");
tree1->Draw("misc.QCDPartonicHT>>h1_HT","(0.1947) * (NJets>0)","goff");
tree2->Draw("misc.QCDPartonicHT>>h2_HT","(0.1490) * (NJets>0)","goff");
tree3->Draw("misc.QCDPartonicHT>>h3_HT","(0.1018) * (NJets>0)","goff");
*/
/*
//Znunu
tree0->Draw("misc.QCDPartonicHT>>h0_HT","(0.3808) * (NJets>0)","goff");
tree1->Draw("misc.QCDPartonicHT>>h1_HT","(0.5754) * (NJets>0)","goff");
tree2->Draw("misc.QCDPartonicHT>>h2_HT","(0.1769) * (NJets>0)","goff");
tree3->Draw("misc.QCDPartonicHT>>h3_HT","(0.0258) * (NJets>0)","goff");
*/
/*
//QCD
tree0->Draw("misc.QCDPartonicHT>>h0_HT","(4133.29) * (NJets>0)","goff");
tree1->Draw("misc.QCDPartonicHT>>h1_HT","(204.17) * (NJets>0)","goff");
tree2->Draw("misc.QCDPartonicHT>>h2_HT","(5.512) * (NJets>0)","goff");
tree3->Draw("misc.QCDPartonicHT>>h3_HT","(0.2951) * (NJets>0)","goff");
*/
//GJets
tree0->Draw("misc.QCDPartonicHT>>h0_HT","(0.34561) * (NJets>0)","goff");
tree1->Draw("misc.QCDPartonicHT>>h1_HT","(0.0530207) * (NJets>0)","goff");


double max(0.), max0(0.), max1(0.), max2(0.), max3(0.);
max0=h0_HT->GetMaximum();
max1=h1_HT->GetMaximum();
max2=h2_HT->GetMaximum();
max3=h3_HT->GetMaximum();
max  = (max0>max1)?max0:max1;
max  = (max2>max)?max2:max;
max  = (max3>max)?max3:max;
max = 3.*max;
h0_HT->SetMaximum(max);
h1_HT->SetMaximum(max);
h2_HT->SetMaximum(max);
h3_HT->SetMaximum(max);

TCanvas* c_HT = new TCanvas("c_HT","",0,0,1000,1000);
c_HT->SetFrameLineWidth(1);
c_HT -> cd();
gPad->SetLogy(1);
gPad->SetFillStyle(0);
// h0_HT->DrawNormalized("");
// h1_HT->DrawNormalized("same");
h0_HT->Draw("");
h1_HT->Draw("same");
h2_HT->Draw("same");
h3_HT->Draw("same");

}