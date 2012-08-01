{
	TChain *std = new TChain("MassTree");
	TChain *pfnopu = new TChain("MassTree");
	
	pfnopu->Add("~/MT2Analysis/MT2trees/MT2_V01-01-02/20111014_Pileupstudies/DYToEE_M-20_TuneZ2_7TeV-pythia6_Fall11-Peak32PU_PFnoPU.root");
	std->Add("~/MT2Analysis/MT2trees/MT2_V01-01-02/20111014_Pileupstudies/DYToEE_M-20_TuneZ2_7TeV-pythia6_Fall11-Peak32PU.root");

	TH2D* hstd  = new TH2D("hstd" , "", 60, 0, 60, 20, 0, 100);
	TH2D* hnopu = new TH2D("hnopu", "", 60, 0, 60, 20, 0, 100);

	hstd ->Sumw2();
	hnopu->Sumw2();

	pfnopu->Draw("misc.HT:pileUp.PUnumInt>>hnopu", "(NJetsIDLoose+NMuons+NEles)>=2");
	std   ->Draw("misc.HT:pileUp.PUnumInt>>hstd",  "(NJetsIDLoose+NMuons+NEles)>=2");

	hstd->GetXaxis()->SetTitle("Num InTime PU interactions");
	hstd->GetYaxis()->SetTitle("M_{T2}");

	hstd ->SetLineColor(kBlue);
	hnopu->SetMarkerStyle(20);

	gPad->SetLogy(1);
	hstd->ProfileX()->Draw();
	hnopu->ProfileX()->Draw("sameEX0");



	TLegend *leg = new TLegend(.55,.5,.85,.75);
	leg->SetFillColor(0);
	leg->AddEntry(hstd ,"DYToEE Peak32 - L1FastJet","le");
	leg->AddEntry(hnopu,"DYToEE Peak32 - hybrid (chs + L1FastJet)","le");
	leg->Draw("same");


}
