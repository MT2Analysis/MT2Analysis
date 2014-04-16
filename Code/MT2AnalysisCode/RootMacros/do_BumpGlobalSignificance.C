{
//run via root -l do_bumpGlobalSignificance.C
//this code checks the global significance of a local bump/dip by running toy experiments

	int nTrails=250000;//number of toy experiments
	int counter=0;
	//choose your bump
	//model is either prediction or simulation
	//the bump/dip is between dump_low - dump_high
/*
	double bump_low=38.99;
	double bump_high=41.01;
	TFile * file = TFile::Open("/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/Results/Filtered/DataPredOneHistogram.root");
	TH1D *model = (TH1D*)file->Get("pred");
	TH1D *data = (TH1D*)file->Get("data");
*//*
	double bump_low=6.99;
	double bump_high=7.01;
	TFile * file = TFile::Open("/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/Results/Filtered/PredictionFile_fixemptybin_ISR_TOBTECreweight.root");
	TH1D *model = (TH1D*)file->Get("predSum");
	TH1D *data = (TH1D*)file->Get("predData");
*/
	double bump_low=16.99;
	double bump_high=17.01;
	TFile *_file0 = TFile::Open("/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/LostLepton/tryout/LLonehistogram.root");
	TH1D* model   = (TH1D*) _file0->Get("LLMC");
	TH1D* data    = (TH1D*) _file0->Get("LLdata");

	//get data and model count within bump window and their uncertainties
	double data_count=0;
	double mc_count  =0;
	double mc_err2   =0;
	int    binnr     = -1;
	for(int i=1;i<data->GetNbinsX(); ++i){
		if(data->GetBinLowEdge(i)<bump_low || data->GetBinLowEdge(i)>=bump_high) continue;
		data_count += data->GetBinContent(i);
		mc_count   += model->GetBinContent(i);
		mc_err2    += model->GetBinError(i)*model->GetBinError(i);
		binnr       = i;
		cout << data->GetXaxis()->GetBinLabel(i) << endl;
	}
	mc_err2 += pow(0.1*mc_count,2);//for mc add 10% uncertainty due to luminosity and cross section uncertainties
	cout << "data in bin " << binnr << "  " << data_count << endl;
	cout << "mc   in bin " << binnr << "  " << mc_count   << endl;
	cout << "err  in bin " << binnr << "  " << sqrt(mc_err2)    << endl;

	//double local_sig_bump = (data_count-mc_count)/sqrt(mc_count); //bump - simple
	//double local_sig_bump = (data_count-mc_count)/sqrt(mc_err2); //bump - simple
	//double local_sig_bump = (mc_count-data_count)/sqrt(mc_count);  //dip - simple
	//get local significance by using Bob Cousin's et al. Z_bi (arXiv:physics/0702156)
	double tau = mc_count/mc_err2;
	double n_off = tau*mc_count;
	double P_Bi = TMath::BetaIncomplete(1./(1.+tau),data_count,n_off+1);
	double Z_Bi = sqrt(2.)*TMath::ErfInverse(1.-2.*P_Bi);
	double Z_Bi = sqrt(2.)*TMath::ErfInverse(1.-2.*P_Bi);
	double local_sig_bump = Z_Bi;
       	cout << "bump local significance " << local_sig_bump  << endl;	
	
	TCanvas * c1 = new TCanvas("c1", "", 0, 0, 500, 500);
	model->Draw("hist");

	TH1D* smeared1 = (TH1D*) model->Clone("smeared1");//just plot the smeared distributions for all trials
	TH1D* sigma    = new TH1D("sigma", "", 50, -5, 5); //saves local signifcance for toys for each bin
	TH1D* sigmamax    = new TH1D("sigmamax", "", 80, -8, 8);//saves local signifcance for toys for most extreme bin


	for(int trial=0; trial<nTrails; ++trial){
	double localmax(99.);
	double bg(0), bgerr(0), pseudodata(0),bini(0); string binn;
	//double localmax(-99.);
		for (int i=1; i<=smeared1->GetNbinsX(); ++i){
			double orig = model->GetBinContent(i);
			double origerr = model->GetBinError(i);
			origerr = sqrt(origerr*origerr+pow(0.1*orig,2));
			smeared1->SetBinContent(i, gRandom->PoissonD(orig+gRandom->Gaus(0,origerr)));
			double tn_on = gRandom->PoissonD(orig+gRandom->Gaus(0,origerr));//pseudodata from smeared model + poissonian statistical smearing
			double ttau = orig/origerr;
			double tn_off = ttau*orig;
			double tP_Bi = TMath::BetaIncomplete(1./(1.+ttau),tn_on,tn_off+1);
			double tZ_Bi = sqrt(2.)*TMath::ErfInverse(1.-2.*tP_Bi);
			double local_significance = tZ_Bi;
			sigma->Fill(local_significance);
			//if(local_significance>local_sig_bump) {counter++; break;}//do not do this in loop
			//if(local_significance<local_sig_bump) {counter++; break;} //do not do this in loop
			//for dip (negative significance) do
			if(local_significance<localmax) {localmax = local_significance; bini=i;binn=smeared1->GetXaxis()->GetBinLabel(i); bg = orig; bgerr = origerr; pseudodata=tn_on;}
			//for bump (positive significance) do
			//if(local_significance>localmax) {localmax = local_significance; bini=i;binn=smeared1->GetXaxis()->GetBinLabel(i); bg = orig; bgerr = origerr; pseudodata=tn_on;}
		}
		sigmamax->Fill(localmax);
		if(localmax<local_sig_bump) ++counter;//this saves the probability of such a bump/dip
	}

	smeared1->SetMarkerStyle(4);
	smeared1->SetLineColor(1);
	smeared1->SetMarkerColor(1);
	smeared1->Draw("same");

	TCanvas * c2 = new TCanvas("c2", "", 0, 0, 500, 500);
	sigma->SetMarkerStyle(4);
	sigma->Draw();

	TCanvas * c3 = new TCanvas("c3", "", 0, 0, 500, 500);
	sigmamax->SetMarkerStyle(4);
	sigmamax->Draw();

	//result
	cout << nTrails << " Trials gave " << counter  << " times a local significance > " << local_sig_bump << " i.e. " << double(counter)/double(nTrails)*100. << " percent for " << smeared1->GetNbinsX() << " regions" << endl;

}