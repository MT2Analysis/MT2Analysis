#include <iostream>
#include <sstream>

#include <TH2.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "helper/Utilities.hh"
#include "Utilities.hh"

#include <TTreeFormula.h>
#include <TTreePerfStats.h>

#include "interface/MT2Common.h"
#include "interface/MT2Region.h"
#include "interface/MT2RegionAnalysisUtilities.h"

#define treeProducerSusyFullHad_cxx
#include "treeProducerSusyFullHad.h"

int main( int argc, char* argv[] ) {


  if( argc!=2 ) {
    std::cout << "USAGE: ./inclusivePlots_BabyTrees [samplesFileName]" << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  
  bool isReduceTree = false;

  std::string sampleName(argv[1]);
  
  std::string samplesFileName = "samples/samples_" + sampleName + ".dat";
  std::cout << std::endl << std::endl;
  std::cout << "-> Loading samples from file: " << samplesFileName << std::endl;
  
  std::vector<MT2SampleBaby_basic> fSamples = MT2Common::loadSamplesBaby_basic(samplesFileName);

  
  std::string outputdir = "InclusivePlots/T1bbbb_1500-100";
  system(Form("mkdir -p %s", outputdir.c_str()));


  std::cout << "-> Starting declaration of stack histograms..." << std::endl;

  // Declaring stack histogram - one per variable
  THStack* h_njets_stack = new THStack("h_njets_stack", "");

  THStack* h_nbjets_stack = new THStack("h_nbjets_stack", "");

  THStack* h_nleptons_stack = new THStack("h_nleptons_stack", "");

  THStack* h_ht_stack = new THStack("h_ht_stack", "");

  THStack* h_met_stack = new THStack("h_met_stack", "");

  THStack* h_mt2_stack = new THStack("h_mt2_stack", "");

  THStack* h_jetpt_stack = new THStack("h_jetpt_stack", "");


  // Declaring histograms - one per sample
  int nSamples = fSamples.size(); 
  
  cout<< "nSamples " << nSamples << endl;

  TH1F* h_njets[nSamples+1];
  
  TH1F* h_nbjets[nSamples+1];
  
  TH1F* h_nleptons[nSamples+1];

  TH1F* h_ht[nSamples+1];

  TH1F* h_met[nSamples+1];

  TH1F* h_mt2[nSamples+1];

  TH1F* h_jetpt[nSamples+1];

  
  for( unsigned i=0; i<fSamples.size(); ++i ){
    
    //Initializing histograms:
    
    string name_njets="njets"+fSamples[i].name;
    h_njets[i] = new TH1F(name_njets.c_str(), "number of jets", 12, 0, 12 );
    h_njets[i]->GetXaxis()->SetTitle("N(jets)");
    h_njets[i]->GetYaxis()->SetTitle("Events");
    if(fSamples[i].sname == "WJets"){
      h_njets[i]->SetLineColor(3);
      h_njets[i]->SetFillColor(3);
    }
    else if(fSamples[i].sname == "Top"){
      h_njets[i]->SetLineColor(4);
      h_njets[i]->SetFillColor(4);
    }
    else if(fSamples[i].sname == "QCD"){
      h_njets[i]->SetLineColor(5);
      h_njets[i]->SetFillColor(5);
    }
      
    string name_nbjets="nbjets"+fSamples[i].name;
    h_nbjets[i] = new TH1F(name_nbjets.c_str(), "number of b-jets", 6, 0, 6 );
    h_nbjets[i]->GetXaxis()->SetTitle("N(b-jets)");
    h_nbjets[i]->GetYaxis()->SetTitle("Events");
    if(fSamples[i].sname == "WJets"){
      h_nbjets[i]->SetLineColor(3);
      h_nbjets[i]->SetFillColor(3);
    }
    else if(fSamples[i].sname == "Top"){
      h_nbjets[i]->SetLineColor(4);
      h_nbjets[i]->SetFillColor(4);
    }
    else if(fSamples[i].sname == "QCD"){
      h_nbjets[i]->SetLineColor(5);
      h_nbjets[i]->SetFillColor(5);
    }

    string name_nleptons="nleptons"+fSamples[i].name;
    h_nleptons[i] = new TH1F(name_nleptons.c_str(), "number of leptons", 10, 0, 10 );
    h_nleptons[i]->GetXaxis()->SetTitle("N(leptons)");
    h_nleptons[i]->GetYaxis()->SetTitle("Events");
    if(fSamples[i].sname == "WJets"){
      h_nleptons[i]->SetLineColor(3);
      h_nleptons[i]->SetFillColor(3);
    }
    else if(fSamples[i].sname == "Top"){
      h_nleptons[i]->SetLineColor(4);
      h_nleptons[i]->SetFillColor(4);
    }
    else if(fSamples[i].sname == "QCD"){
      h_nleptons[i]->SetLineColor(5);
      h_nleptons[i]->SetFillColor(5);
    }
    
    string name_ht="ht"+fSamples[i].name;
    h_ht[i] = new TH1F(name_ht.c_str(), "H_{T}", 300, 0, 3000 );
    h_ht[i]->GetXaxis()->SetTitle("H_{T} [GeV]");
    h_ht[i]->GetYaxis()->SetTitle("Events/10 GeV");
    if(fSamples[i].sname == "WJets"){
      h_ht[i]->SetLineColor(3);
      h_ht[i]->SetFillColor(3);
    }
    else if(fSamples[i].sname == "Top"){
      h_ht[i]->SetLineColor(4);
      h_ht[i]->SetFillColor(4);
    }
    else if(fSamples[i].sname == "QCD"){
      h_ht[i]->SetLineColor(5);
      h_ht[i]->SetFillColor(5);
    }

    string name_met="met"+fSamples[i].name;
    h_met[i] = new TH1F(name_met.c_str(), "E_{T}^{miss}", 100, 0, 1000 );
    h_met[i]->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
    h_met[i]->GetYaxis()->SetTitle("Events/10 GeV");
    if(fSamples[i].sname == "WJets"){
      h_met[i]->SetLineColor(3);
      h_met[i]->SetFillColor(3);
    }
    else if(fSamples[i].sname == "Top"){
      h_met[i]->SetLineColor(4);
      h_met[i]->SetFillColor(4);
    }
    else if(fSamples[i].sname == "QCD"){
      h_met[i]->SetLineColor(5);
      h_met[i]->SetFillColor(5);
    }

    string name_mt2="mt2"+fSamples[i].name;;
    h_mt2[i] = new TH1F(name_mt2.c_str(), "M_{T2}", 100, 0, 1000 );
    h_mt2[i]->GetXaxis()->SetTitle("M_{T2} [GeV]");
    h_mt2[i]->GetYaxis()->SetTitle("Events/10 GeV");
    if(fSamples[i].sname == "WJets"){
      h_mt2[i]->SetLineColor(3);
      h_mt2[i]->SetFillColor(3);
    }
    else if(fSamples[i].sname == "Top"){
      h_mt2[i]->SetLineColor(4);
      h_mt2[i]->SetFillColor(4);
    }
    else if(fSamples[i].sname == "QCD"){
      h_mt2[i]->SetLineColor(5);
      h_mt2[i]->SetFillColor(5);
    }

    string name_jetpt="jetpt"+fSamples[i].name;;
    h_jetpt[i] = new TH1F(name_jetpt.c_str(), "p_{T} (jet)", 300, 0, 3000 );
    h_jetpt[i]->GetXaxis()->SetTitle("p_{T}(jet) [GeV]");
    h_jetpt[i]->GetYaxis()->SetTitle("Events/10 GeV");
    if(fSamples[i].sname == "WJets"){
      h_jetpt[i]->SetLineColor(3);
      h_jetpt[i]->SetFillColor(3);
    }
    else if(fSamples[i].sname == "Top"){
      h_jetpt[i]->SetLineColor(4);
      h_jetpt[i]->SetFillColor(4);
    }
    else if(fSamples[i].sname == "QCD"){
      h_jetpt[i]->SetLineColor(5);
      h_jetpt[i]->SetFillColor(5);
    }


    std::cout << std::endl << std::endl;
    std::cout << "-> Starting filling histograms for sample: " << fSamples[i].name << std::endl;

    TFile* file = TFile::Open(fSamples[i].file.c_str());
    TTree* tree = (TTree*)file->Get("treeProducerSusyFullHad");

    TFile* tmpFile = TFile::Open("tmp.root", "recreate");
    TTree* tree_reduced;

    if(isReduceTree){
    
      // Define selection if reducing the tree before looping over entries:

      std::ostringstream preselectionStream;
      preselectionStream << " "
	//<< "(nTaus20==0 && nMuons10==0 && nElectrons10==0)"                   << " && "
			 << "(nVert > 0)";                      
        //<< " && "
        //		       << "(nJet40 > 1)"                     << " && "
        //		       << "(jet_pt[1] > 100)"                << " && "
        //		       << "(deltaPhiMin > 0.3)"              << " && "
        //		       << "(diffMetMht < 70)";
      
      TString preselection = preselectionStream.str().c_str();
      TString cuts = preselection;
      
      /*
	TFile* tmpFile = TFile::Open("tmp.root", "recreate");
	tmpFile->cd();
	TTree* tree_reduced = tree->CopyTree(cuts);
      */

      tmpFile->cd();
      tree_reduced = tree->CopyTree(cuts);
      
    }// reduceTree

    treeProducerSusyFullHad myTree;
    if(isReduceTree)
      myTree.Init(tree_reduced);
    else
      myTree.Init(tree);

    // global sample weight:
    //Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents*fSamples[i].PU_avg_weight);
    Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
    
    //std::cout << "Weight = " << weight << std::endl;
    //cout<<fSamples[i].xsection<<endl;
    //cout<<fSamples[i].kfact<<endl;
    //cout<<fSamples[i].lumi<<endl;
    //cout<<fSamples[i].nevents<<endl;
    //cout<<fSamples[i].PU_avg_weight<<endl;

    int nentries;

    // Performance statistics
    TTreePerfStats *ps;
    
    if(isReduceTree){      

      nentries = tree_reduced->GetEntries();
      ps = new TTreePerfStats("ioperf", tree_reduced);
      
    }
    else {

      nentries = tree->GetEntries(); 
      ps = new TTreePerfStats("ioperf", tree);
    
    }
    
    std::cout << std::endl << std::endl;
    std::cout << "-> Starting loop over " << nentries << " entries..." << std::endl;

    for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

      if( iEntry % 50000 == 0 ) std::cout << "    Entry: " << iEntry << " / " << nentries << std::endl;

      myTree.GetEntry(iEntry);
      
      if(!isReduceTree){
	
	//Define event selection if NOT reducing the input tree:
	if ( myTree.nVert < 1 ) continue;
	
      }
      
      float ht       = myTree.ht;
      float met      = myTree.met_pt;
      float mt2      = myTree.mt2;
      int njets      = myTree.nJet40;
      int nbjets     = myTree.nBJet40;
      int nmuons     = myTree.nMuons10;
      int nelectrons = myTree.nElectrons10;
      int ntaus      = myTree.nTaus20;
      int nleptons   = nmuons+nelectrons+ntaus;
      int njet_all   = myTree.njet;
      
      h_njets[i]    ->Fill(njets, weight);
      h_nbjets[i]   ->Fill(nbjets, weight);
      h_nleptons[i] ->Fill(nleptons, weight);
      h_ht[i]       ->Fill(ht, weight);
      h_met[i]      ->Fill(met, weight);
      h_mt2[i]      ->Fill(mt2, weight);

      for (int jeti=0; jeti<njet_all; ++jeti)
	h_jetpt[i]  ->Fill(myTree.jet_pt[jeti]);
      
      
    }// entries
    
    std::cout << "Done with sample " << fSamples[i].name << std::endl << std::endl; 
 
    h_njets_stack    ->Add(h_njets[i]);
    h_nbjets_stack   ->Add(h_nbjets[i]);
    h_nleptons_stack ->Add(h_nleptons[i]);
    h_ht_stack       ->Add(h_ht[i]);
    h_met_stack      ->Add(h_met[i]);
    h_mt2_stack      ->Add(h_mt2[i]);
    h_jetpt_stack    ->Add(h_jetpt[i]);

    // Performance statistics
    std::string ps_name = outputdir+"/perfstat"+fSamples[i].name+".root";
    ps->SaveAs(ps_name.c_str());
    
    delete ps;

    if(isReduceTree)
      delete tree_reduced;
    delete tree;
    
    tmpFile->Close();
    delete tmpFile;
       
    file->Close();
    delete file;

  }// samples

  system( "rm tmp.root" );

  // Declaring empty histos for legend definition
  TH1F* hWJets = new TH1F("WJets", "WJets", 1, 0, 1);
  hWJets->SetFillColor(3);

  TH1F* hTop = new TH1F("Top", "Top", 1, 0, 1);
  hTop->SetFillColor(4);

  TH1F* hQCD = new TH1F("QCD", "QCD", 1, 0, 1);
  hQCD->SetFillColor(5);

  TH1F* hSignal = new TH1F("Signal", "Signal", 1, 0, 1);

  TLegend* Leg = new TLegend(.65,.65,.85,.85, "");
  Leg->SetFillColor(0);
  //Leg->AddEntry(hWJets, "W + jets", "F");
  //Leg->AddEntry(hTop, "Top", "F");
  //Leg->AddEntry(hQCD, "QCD", "F");
  Leg->AddEntry(hSignal, "T1bbbb - 1500, 100", "F");

  TCanvas* c1=new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  h_njets_stack ->Draw("hist");
  h_njets_stack ->GetXaxis()->SetTitle("N(jets)");
  h_njets_stack ->GetYaxis()->SetTitle("Events");
  h_njets_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_njets_stack ->Draw("hist");

  Leg->Draw("same");
  
  c1->Update();

  c1->SaveAs(Form("%s/inclusiveNJets.eps", outputdir.c_str()));
  
  TCanvas* c2=new TCanvas("c2", "c2", 600, 600);
  c2->cd();
  h_nbjets_stack ->Draw("hist");
  h_nbjets_stack ->GetXaxis()->SetTitle("N(b-jets)");
  h_nbjets_stack ->GetYaxis()->SetTitle("Events");
  h_nbjets_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_nbjets_stack ->Draw("hist");

  Leg->Draw("same");
  
  c2->Update();
  c2->SaveAs(Form("%s/inclusiveNbjets.eps", outputdir.c_str()));

  TCanvas* c3=new TCanvas("c3", "c3", 600, 600);
  c3->cd();
  h_nleptons_stack ->Draw("hist");
  h_nleptons_stack ->GetXaxis()->SetTitle("N(leptons)");
  h_nleptons_stack ->GetYaxis()->SetTitle("Events");
  h_nleptons_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_nleptons_stack ->Draw("hist");
  
  Leg->Draw("same");

  c3->Update();
  c3->SaveAs(Form("%s/inclusiveNleptons.eps", outputdir.c_str()));

  TCanvas* c4=new TCanvas("c4", "c4", 600, 600);
  c4->cd();
  gPad->SetLogy();

  h_ht_stack ->Draw("hist");
  h_ht_stack ->GetXaxis()->SetTitle("H_{T} [GeV]");
  h_ht_stack ->GetYaxis()->SetTitle("Events/10 GeV");
  h_ht_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_ht_stack ->SetMinimum(0.1);
  h_ht_stack ->Draw("hist");

  Leg->Draw("same");

  c4->Update();
  c4->SaveAs(Form("%s/HTinclusive.eps", outputdir.c_str()));

  TCanvas* c5=new TCanvas("c5", "c5", 600, 600);
  c5->cd();
  gPad->SetLogy();

  h_met_stack ->Draw("hist");
  h_met_stack ->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  h_met_stack ->GetYaxis()->SetTitle("Events/10 GeV");
  h_met_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_met_stack ->SetMinimum(0.1);
  h_met_stack ->Draw("hist");
  
  Leg->Draw("same");
  
  c5->Update();
  c5->SaveAs(Form("%s/METinclusive.eps", outputdir.c_str()));

  TCanvas* c6=new TCanvas("c6", "c6", 600, 600);
  c6->cd();
  gPad->SetLogy();

  h_mt2_stack ->Draw("hist");
  h_mt2_stack ->GetXaxis()->SetTitle("M_{T2} [GeV]");
  h_mt2_stack ->GetYaxis()->SetTitle("Events/10 GeV");
  h_mt2_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_mt2_stack ->SetMinimum(0.1);
  h_mt2_stack ->Draw("hist");

  Leg->Draw("same");

  c6->Update();
  c6->SaveAs(Form("%s/MT2inclusive.eps", outputdir.c_str()));

  TCanvas* c7=new TCanvas("c7", "c7", 600, 600);
  c7->cd();
  gPad->SetLogy();

  h_jetpt_stack ->Draw("hist");
  h_jetpt_stack ->GetXaxis()->SetTitle("p_{T}(jet) [GeV]");
  h_jetpt_stack ->GetYaxis()->SetTitle("Events/10 GeV");
  h_jetpt_stack ->GetYaxis()->SetTitleOffset(1.5);
  h_jetpt_stack ->SetMinimum(0.1);
  h_jetpt_stack ->Draw("hist");

  Leg->Draw("same");

  c7->Update();
  c7->SaveAs(Form("%s/jetpt_inclusive.eps", outputdir.c_str()));


  return 0;

}// main
