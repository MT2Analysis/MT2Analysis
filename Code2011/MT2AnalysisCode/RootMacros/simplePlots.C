{

  TString q = "Sum$(jet[].lv.Pt()>20)"; TString bin="(20,0,20)";
  TString label ="NJetsPt20";

  //TString q = "misc.MT2"; TString bin="(20,0,100)";
  //TString label ="MT2";

  //TString F = "/shome/pnef/MT2Analysis/MT2trees/MT2_V01-01-01/20111004_MC_HT400_data_nocuts/skimmed_v1/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11_2.root"; TString FLab="TTBar_LowMT2";
  //TString F = "/shome/pnef/MT2Analysis/MT2trees/MT2_V01-01-01/20111004_MC_HT400_data_nocuts/skimmed_v1/ZJetsToNuNu_200_HT_inf_7TeV-madgraph_Summer11.root";

  //TString F = "/shome/leo/MT2Analysis/MT2trees/MT2_V01-01-02/20111017_Pileupstudies_2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11_2.root"; TString FLab="TTBar";  
  TString F = "/shome/leo/MT2Analysis/MT2trees/MT2_V01-01-02/20111017_Pileupstudies_PFnoPU/TTJets_TuneZ2_7TeV-madgraph-tauola-Summer11.root"; TString FLab="TTBar_PFnoPU";

  //TString F = "/scratch/leo/DYToEE_M-20_TuneZ2_7TeV-pythia6_Fall11-Peak32PU_PFnoPU.root"; TString FLab="DY_PFnoPU";


  int REBIN=1;
  bool ratio = false;

  TChain *t = new TChain("MassTree");
  t->Add(F);
  
  
  std::ostringstream cutStream;
  cutStream       << " " 
    //	  << "NBJets >0"                               << "&&"
    //<< "Sum$(ele[].lv.Pt()>10) == 0" << " && "
    //	  << "Sum$(muo[].lv.Pt()>10) == 0" << " && "
    << "misc.caloHT50_ID > 600"             << "&&"
    //	  << "misc.HT  > 650"                     << "&&"
    << "misc.Jet0Pass ==1"                  << "&&"
		  << "misc.Jet1Pass ==1"                  << "&&"
    //<< "NJetsIDLoose >=3"                   << "&&"
    	  << "misc.SecondJPt > 100"               << "&&"
    	  << "misc.MinMetJetDPhi >0.3"            << "&&"
		  << "misc.PassJetID ==1"                 << "&&"
		  << "misc.Vectorsumpt < 70"              << "&&"
    		  << "misc.MET>=30"                       << "&&"
    //	<< "misc.MT2 > 300"                     << "&&"
		  << "misc.HBHENoiseFlag == 0"            << "&&"
		  << "misc.CrazyHCAL==0";
  
  
  TString cuts    = cutStream.str().c_str();
  cuts="1==1 ";
  cout << cuts << endl;
  /*
    t   ->Draw("misc.MT2>>h1(70, 0, 700)", cuts+" && pileUp.NVertices==1");
    h1->SetName("h1");
    t   ->Draw("misc.MT2>>h2(70, 0, 700)", cuts+" && pileUp.NVertices>7");
    h2->SetName("h2");
  */
  
  t   ->Draw(q+">>h1"+bin, cuts+" && pileUp.NVertices<=10");
  h1->SetName("h1");
  t   ->Draw(q+">>h2"+bin, cuts+" && pileUp.NVertices>12 && pileUp.NVertices<18");
  h2->SetName("h2");
  t   ->Draw(q+">>h3"+bin, cuts+" && pileUp.NVertices>20");
  h3->SetName("h3");


  q.ReplaceAll("[",1,"",1);
  q.ReplaceAll("]",1,"",1);
  q.ReplaceAll("(",1,"",1);
  q.ReplaceAll(")",1,"",1);
  q.ReplaceAll(">",1,"",1);

  TCanvas *c1 = new TCanvas("histo_"+FLab+"_"+label,"histo_"+FLab+"_"+label);
  if(ratio){
    c1->Divide(2);
    c1->cd(1);
  }
  h2 ->SetLineColor(kBlue);
  //  h2 ->SetMarkerStyle(20);
  h2 ->SetMarkerColor(kBlue);
  h3 ->SetLineColor(kRed);
  h3 ->SetMarkerColor(kRed);
  h1 ->SetLineColor(kBlack);
  

  cout << " h1 " << h1->Integral()<< endl;
  cout << " h2 " << h2->Integral()<< endl;
  cout << " h3 " << h3->Integral()<< endl;

  float i1 = h1->Integral();
  float i2 = h2->Integral();
  float i3 = h3->Integral();

  h1->Scale(1./ i1 );
  h2->Scale(1./ i2 );
  h3->Scale(1./ i3 );

  h1->GetXaxis()->SetTitle(label);
  h1->Draw("");
  h2->Draw("hsames");
  h3->Draw("hsames");
  
  
  TLegend* Legend1 = new TLegend(.71,.60,.91,.92);
  Legend1 ->AddEntry(h1, "nV<=10"  , "l");
  Legend1 ->AddEntry(h2, "nV>12 && nV<18"   , "l");
  Legend1 ->AddEntry(h3, "nV>20"   , "l");
  Legend1 ->SetBorderSize(0.00001);
  Legend1 ->SetFillColor(kWhite);
  Legend1 ->Draw();
  
  
  if(ratio){
    c1->cd(2);
    
    TH1F *hc = (TH1F*) h2->Clone();
    hc->GetYaxis()->SetTitle("ratio");
    hc->Divide(h1);
    hc->Draw();
    TH1F *hc2 = (TH1F*) h3->Clone();
    hc2->Divide(h1);
    hc2->Draw("sames");

    
  }
  
}
