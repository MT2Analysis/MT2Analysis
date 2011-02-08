
void comp_dist_from_single_tree(){
	setStyle();
	TString      basecut =  " misc.LeptConfig==9 &&";
	  //         basecut += " misc.Jet0Pass==1 && misc.Jet1Pass==1 &&";
	             basecut += " GetNjets(20,2.0,0)>=2 && ";
	             basecut += " misc.MET >= 0 ";
	  //         basecut += " misc.PassJetID==1";

	TChain *chain = new TChain("MassTree");
//	chain ->Add("testQCDverysmall/QCD_Pt_170to300_TuneZ2_7TeV_pythia6_V01-11-02_small.root");
//	chain ->Add("testDataverysmall/MultiJet-Run2010B-Nov4ReReco_v1_RECO_V01-11-02_small.root");
	chain ->Add("testLM1verysmall/SUSY_LM1_small.root");
//	chain ->Add("testQCD/QCD_Pt_170to300_TuneZ2_7TeV_pythia6_V01-11-02_small.root");
  	
	TH1D* h1 = GetHisto(chain, "GetMT2Hemi(0,false,0,20,2.4,3,1)     >>", basecut, "MT2_0_20_24_3_MET");	
	TH1D* h2 = GetHisto(chain, "GetMT2Hemi(0,false,0,20,3.0,3,1)     >>", basecut, "MT2_0_20_30_3_MET");	
	TH1D* h3 = GetHisto(chain, "GetMT2HemiMinDHT(0,false,0,20,2.4,1) >>", basecut, "MT2_0_20_24_minDHT_MET");	
	TH1D* h4 = GetHisto(chain, "GetMT2HemiMinDHT(0,false,0,20,3.0,1) >>", basecut, "MT2_0_20_30_minDHT_MET");	
  
	h1->SetLineColor(kGreen);
	h2->SetLineColor(kGreen+2);
	h3->SetLineColor(kBlue);
	h4->SetLineColor(kBlue+2);

	h2->SetLineStyle(kDashed);
	h4->SetLineStyle(kDashed);

	TLegend* leg = new TLegend(0.2, 0.6, 0.5, 0.9);	
	leg->AddEntry(h1,h1->GetName(),"l" );
	leg->AddEntry(h2,h2->GetName(),"l" );
	leg->AddEntry(h3,h3->GetName(),"l" );
	leg->AddEntry(h4,h4->GetName(),"l" );

	TCanvas *c1 = new TCanvas();
  	gPad->SetLogy(1);
	h1->SetXTitle("MT2");
  	h1->Draw("");
  	h2->Draw("same");
	h3->Draw("same");
	h4->Draw("same");
	leg->Draw();
}

TH1D* GetHisto(TChain* chain, TString var, TString basecut, TString name){
	TH1D* h = new TH1D(name,"", 100, 0, 800);
	TString varname = var + h->GetName();
	int n = chain->Draw(varname,basecut,"goff");
	cout << name << " " << h->Integral() << endl;
	h->Scale(1./h->Integral());
	return h;
}

void setStyle(){
	gROOT->ProcessLine(".x ~casal/SetStyle_PRD.C");
	gStyle->SetLegendBorderSize(0);
	gStyle->SetFillStyle(0);
	gStyle->SetTextFont(62);
	gStyle->SetTextSize(0.045);
	gStyle->SetPalette(1);
}
