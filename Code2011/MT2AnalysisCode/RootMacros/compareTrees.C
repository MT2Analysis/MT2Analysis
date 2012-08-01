{
	TChain *std = new TChain("MassTree");
	TChain *pfnopu = new TChain("MassTree");
	
	pfnopu->Add("~/MT2Analysis/MT2trees/MT2_V01-01-02/20111014_Pileupstudies/DYToEE_M-20_TuneZ2_7TeV-pythia6_Fall11-Peak32PU_PFnoPU.root");
	std->Add("~/MT2Analysis/MT2trees/MT2_V01-01-02/20111014_Pileupstudies/DYToEE_M-20_TuneZ2_7TeV-pythia6_Fall11-Peak32PU.root");

	TH1D* hstd  = new TH1D("hstd" , "", 30, 0, 100);
	TH1D* hnopu = new TH1D("hnopu", "", 30, 0, 100);

	hstd ->Sumw2();
	hnopu->Sumw2();

	pfnopu->Draw("misc.MT2>>hnopu", "(NJetsIDLoose+NMuons+NEles)>=2");
	std   ->Draw("misc.MT2>>hstd",  "(NJetsIDLoose+NMuons+NEles)>=2");

	hstd->GetXaxis()->SetTitle("misc.MT2");

	hstd ->SetLineColor(kBlue);
	hnopu->SetMarkerStyle(20);

	gPad->SetLogy(1);
	hstd ->DrawNormalized("hist");
	hnopu->DrawNormalized("sameEX0");



	TLegend *leg = new TLegend(.55,.5,.85,.75);
	leg->SetFillColor(0);
	leg->AddEntry(hstd ,"DYToEE Peak32 - L1FastJet","le");
	leg->AddEntry(hnopu,"DYToEE Peak32 - hybrid (chs + L1FastJet)","le");
	leg->Draw("same");


}
