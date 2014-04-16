{
	//run via root -l do_BumpGlobalSignificanceComplex.C

	//this does not directly do a 'global significance' calculation,
	//but rather checks how likely a 'global distribution' of all signal bins in our analysis is
	//using pull mean
	//as this is a copy of do_BumpGlobalSignificance.C less comments are made here

	int nTrials=25000;
	int counter=0;

	//hPullData has mean 0.26 +/- 0.07 and rms 0.81 +/- 0.05
	double bump_low=38.99;
	double bump_high=41.01;
//	TFile * file = TFile::Open("/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/Results/Filtered/DataPredOneHistogram.root");
//	TFile * file = TFile::Open("/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/Results/Filtered/PredictionFile_fixemptybin_ISR_TOBTECreweight_PoissonErrors_AllInOne.root");
	TFile * file = TFile::Open("/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/Results/Filtered/PredictionFile_fixemptybin_ISR_AllInOne.root");
	TH1D *model = (TH1D*)file->Get("hPredAll");
	TH1D *data = (TH1D*)file->Get("hDataAll");
	TH1D *zinv = (TH1D*)file->Get("hZinvAll");
	TH1D *qcd = (TH1D*)file->Get("hQCDAll");
	TH1D *LL = (TH1D*)file->Get("hLLAll");
	//need all the histograms to cover correlations

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
	mc_err2 += pow(0.1*mc_count,2);
	cout << "data in bin " << binnr << "  " << data_count << endl;
	cout << "mc   in bin " << binnr << "  " << mc_count   << endl;
	cout << "err  in bin " << binnr << "  " << sqrt(mc_err2)    << endl;

	//actually this is not needed here!!!
	//double local_sig_bump = (data_count-mc_count)/sqrt(mc_count); //bump
	//double local_sig_bump = (data_count-mc_count)/sqrt(mc_err2); //bump
	//double local_sig_bump = (mc_count-data_count)/sqrt(mc_count);  //dip
	double tau = mc_count/mc_err2;
	double n_off = tau*mc_count;
	double P_Bi = TMath::BetaIncomplete(1./(1.+tau),data_count,n_off+1);
	double Z_Bi = sqrt(2.)*TMath::ErfInverse(1.-2.*P_Bi);
	double Z_Bi = sqrt(2.)*TMath::ErfInverse(1.-2.*P_Bi);
	double local_sig_bump = Z_Bi;
       	cout << "bump local significance " << local_sig_bump  << endl;	
	
	TCanvas * c1 = new TCanvas("c1", "", 0, 0, 500, 500);
	model->Draw("hist");

	TH1D* smeared1 = (TH1D*) model->Clone("smeared1");
	TH1D* smeared2 = (TH1D*) model->Clone("smeared2");
	TH1D* sigma    = new TH1D("sigma", "", 50, -5, 5);
	TH1D* sigmamax    = new TH1D("sigmamax", "", 80, -8, 8);

	TH1D* pullmean    = new TH1D("pullmean", "", 400, -1., 1.); pullmean->Sumw2(); pullmean->GetXaxis()->SetTitle("pull mean"); pullmean->GetYaxis()->SetTitle("trials");
	TH1D* pullrms     = new TH1D("pullrms" , "", 400,  0., 2.); pullrms ->Sumw2(); pullrms ->GetXaxis()->SetTitle("pull rms" ); pullrms ->GetYaxis()->SetTitle("trials");

	for(int trial=0; trial<nTrials; ++trial){
	if(trial%(nTrials/25)==0) cout << "trial " << trial << "/" << nTrials << endl;
	TH1D *pulltemp = new TH1D("pulltemp","",50,-5.,5.); pulltemp->Sumw2();
	double localmax(99.);
	double bg(0), bgerr(0), pseudodata(0),bini(0); string binn;
	double LLerr(0), Zinv20err(0), Zinv1b0berr(0);
		for (int i=1; i<=smeared1->GetNbinsX(); ++i){
			//get correlations
			double zinv20err = 0.2;
			if((i>=29&&i<33)||(i>=39&&i<44)||(i>=73&&i<78)||(i>=87&&i<94)||(i>=113&&i<115)||(i>=121&&i<=123)) zinv20err = 0;
			double LLrelerr  = LL->GetBinError(i)/LL->GetBinContent(i);
			double zinv1b0berr=0;
			if(i>=  9&&i< 15) zinv1b0berr = 0.022/0.093;
			if(i>= 23&&i< 29) zinv1b0berr = 0.024/0.162;
			if(i>= 36&&i< 39) zinv1b0berr = 0.168/0.269;
			if(i>= 53&&i< 58) zinv1b0berr = 0.050/0.093;
			if(i>= 67&&i< 73) zinv1b0berr = 0.018/0.162;
			if(i>= 83&&i< 87) zinv1b0berr = 0.143/0.269;
			if(i>=100&&i<102) zinv1b0berr = 0.060/0.093;
			if(i>=109&&i<113) zinv1b0berr = 0.030/0.162;
			if(i>=118&&i<121) zinv1b0berr = 0.201/0.269;
			double orig = model->GetBinContent(i);
			double origerr = model->GetBinError(i);
			double origuncorrerr2 = origerr*origerr-LL->GetBinError(i)*LL->GetBinError(i)-zinv20err*zinv20err*zinv->GetBinContent(i)*zinv->GetBinContent(i)-zinv1b0berr*zinv1b0berr*zinv->GetBinContent(i)*zinv->GetBinContent(i);
			if(origuncorrerr2<0) cout << "ERROR!!! " << origuncorrerr2 << endl;
			double origuncorrerr = sqrt(origuncorrerr2);
			if(i==1) Zinv20err = gRandom->Gaus(0,zinv20err);
			double Zinv20errUsed = Zinv20err;
			if((i>=29&&i<33)||(i>=39&&i<44)||(i>=73&&i<78)||(i>=87&&i<94)||(i>=113&&i<115)||(i>=121&&i<=123)) Zinv20errUsed = 0;
			if(i==  1) Zinv1b0berr = 0; if(i==  9) Zinv1b0berr = gRandom->Gaus(0,zinv1b0berr);
			if(i== 15) Zinv1b0berr = 0; if(i== 23) Zinv1b0berr = gRandom->Gaus(0,zinv1b0berr);
			if(i== 29) Zinv1b0berr = 0; if(i== 36) Zinv1b0berr = gRandom->Gaus(0,zinv1b0berr);
			if(i== 39) Zinv1b0berr = 0; if(i== 53) Zinv1b0berr = gRandom->Gaus(0,zinv1b0berr);
			if(i== 58) Zinv1b0berr = 0; if(i== 67) Zinv1b0berr = gRandom->Gaus(0,zinv1b0berr);
			if(i== 73) Zinv1b0berr = 0; if(i== 83) Zinv1b0berr = gRandom->Gaus(0,zinv1b0berr);
			if(i== 87) Zinv1b0berr = 0; if(i==100) Zinv1b0berr = gRandom->Gaus(0,zinv1b0berr);
			if(i==102) Zinv1b0berr = 0; if(i==109) Zinv1b0berr = gRandom->Gaus(0,zinv1b0berr);
			if(i==113) Zinv1b0berr = 0; if(i==118) Zinv1b0berr = gRandom->Gaus(0,zinv1b0berr); if(i==121) Zinv1b0berr = 0;
			if(i==  1) LLerr = gRandom->Gaus(0,LLrelerr); if(i==  9) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 15) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i== 23) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 29) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 33) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i== 36) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 39) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 42) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i== 44) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 53) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 58) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i== 67) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 73) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 78) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i== 83) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 87) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 91) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i== 94) LLerr = gRandom->Gaus(0,LLrelerr); if(i==100) LLerr = gRandom->Gaus(0,LLrelerr); if(i==102) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i==109) LLerr = gRandom->Gaus(0,LLrelerr); if(i==113) LLerr = gRandom->Gaus(0,LLrelerr); if(i==115) LLerr = gRandom->Gaus(0,LLrelerr); 

			double origuncorrsmearederr = gRandom->Gaus(0,origuncorrerr);
			//tn_on considers most important correlations
			double tn_on = gRandom->PoissonD(orig+origuncorrsmearederr+LLerr*LL->GetBinContent(i)+Zinv1b0berr*zinv->GetBinContent(i)+Zinv20errUsed*zinv->GetBinContent(i));
			double pullpseudodata = (tn_on - orig) / sqrt(pow(origerr,2) + tn_on);//pull with pseudodata
			pulltemp->Fill(pullpseudodata);

		}
		//FILL HERE PULL MEAN HISTO, PULL RMS HISTO
		double mean = pulltemp->GetMean(); double meanerr = pulltemp->GetMeanError();
		double rms = pulltemp->GetRMS(); double rmserr = pulltemp->GetRMSError();
		if(trial%(nTrials/25)==0) cout << "mean " << mean << "+/-" << meanerr << ", rms " << rms << "+/-" << rmserr << endl;
		delete pulltemp;
		pullmean->Fill(mean);
		pullrms ->Fill(rms);
		if(0.26<mean) ++counter;//mean of our analysis is 0.26
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

	TCanvas * c4 = new TCanvas("c4", "", 0, 0, 500, 500);
	pullmean->SetMarkerStyle(4);
	pullmean->Draw();

	TCanvas * c5 = new TCanvas("c5", "", 0, 0, 500, 500);
	pullrms->SetMarkerStyle(4);
	pullrms->Draw();


	cout << nTrials << " Trials gave " << counter  << " times a local significance > " << local_sig_bump << " i.e. " << double(counter)/double(nTrials)*100. << " percent for " << smeared1->GetNbinsX() << " regions" << endl;


}