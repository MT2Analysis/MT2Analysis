{
	//run via root -l do_BumpGlobalSignificanceReallyComplex.C
	//this code is the same asdo_BumpGlobalSignificanceComplex.C
	//but takes into account all correlations that can be implemented by hand.
	//no additional comments are made here

	int nTrials=25000;
	int counter=0;

	double bump_low=38.99;
	double bump_high=41.01;
//	TFile * file = TFile::Open("/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/Results/Filtered/DataPredOneHistogram.root");
	//TFile * file = TFile::Open("/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/Results/Filtered/PredictionFile_fixemptybin_ISR_TOBTECreweight_PoissonErrors_AllInOne.root");
	TFile * file = TFile::Open("/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/Results/Filtered/PredictionFile_fixemptybin_ISR_AllInOne.root");
	TH1D *model = (TH1D*)file->Get("hPredAll");
	TH1D *data = (TH1D*)file->Get("hDataAll");
	TH1D *zinv = (TH1D*)file->Get("hZinvAll");
	TH1D *qcd = (TH1D*)file->Get("hQCDAll");
	TH1D *LL = (TH1D*)file->Get("hLLAll");

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
	cout << "data in bin " << binnr << "  " << data_count << endl;
	cout << "mc   in bin " << binnr << "  " << mc_count   << endl;
	cout << "err  in bin " << binnr << "  " << sqrt(mc_err2)    << endl;
	//double local_sig_bump = (data_count-mc_count)/sqrt(mc_err2); //bump
	double tau = mc_count/mc_err2;
	double n_off = tau*mc_count;
	double P_Bi = TMath::BetaIncomplete(1./(1.+tau),data_count,n_off+1);
	double Z_Bi = sqrt(2.)*TMath::ErfInverse(1.-2.*P_Bi);
	double Z_Bi = sqrt(2.)*TMath::ErfInverse(1.-2.*P_Bi);
	double local_sig_bump = Z_Bi;
       	cout << "bump local significance " << local_sig_bump  << endl;	
	
	TCanvas * c1 = new TCanvas("c1", "", 0, 0, 500, 500);
	model->Draw("hist");

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
	double Zinv1b0berr2(0.), Zinv1b0berr3(0.), Zinv1b0berr6(0.);
		for (int i=1; i<=model->GetNbinsX(); ++i){
			bool thirty = false; bool zero = false;//thirty: 30% instead of 20% for R/Z gamma ratio, same for zero
			if(i>=4&&i<9) thirty = true; else if(i>=18&&i<23) thirty = true; else if(i==35||i==38||i==86||i==112||i==120) thirty = true; else if(i>=49&&i<53) thirty = true; else if(i>=63&&i<67) thirty = true; else if(i>=81&&i<83) thirty = true; else if(i>=98&&i<100) thirty = true; else if(i>=106&&i<109) thirty = true; else if(i>=116&&i<118) thirty = true; else if(i>=11&&i<15) thirty = true; else if(i>=25&&i<=29) thirty = true; else if(i>=56&&i<58) thirty = true; else if(i>=71&&i<73) thirty = true;
			else if(i>=29&&i<33) zero = true; else if(i>=39&&i<44) zero = true; else if(i>=73&&i<78) zero = true; else if(i>=87&&i<94) zero = true; else if(i>=113&&i<115) zero = true; else if(i>=121&&i<=123) zero = true;
			double zinv20err = 0.2;
			if(thirty) zinv20err = 0.3;
			else if(zero) zinv20err = 0;
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
			if(origuncorrerr2<0) cout << "ERROR!!! " << i << " " <<  origuncorrerr2 << " " << zinv20err*zinv->GetBinContent(i) << endl;
			double origuncorrerr = sqrt(origuncorrerr2);
			if(i==1) Zinv20err = gRandom->Gaus(0,zinv20err);
			double Zinv20errUsed = Zinv20err;
			if(thirty) Zinv20errUsed = 1.5 * Zinv20err;
			if(zero) Zinv20errUsed = 0;
			if(i==  1) Zinv1b0berr = 0; if(i==  9) { Zinv1b0berr2 = gRandom->Gaus(0,zinv1b0berr); Zinv1b0berr = Zinv1b0berr2;}
			if(i== 15) Zinv1b0berr = 0; if(i== 23) { Zinv1b0berr3 = gRandom->Gaus(0,zinv1b0berr); Zinv1b0berr = Zinv1b0berr3;}
			if(i== 29) Zinv1b0berr = 0; if(i== 36) { Zinv1b0berr6 = gRandom->Gaus(0,zinv1b0berr); Zinv1b0berr = Zinv1b0berr6;}
			if(i== 39) Zinv1b0berr = 0; if(i== 53) { Zinv1b0berr2 = Zinv1b0berr2*0.050/0.022; Zinv1b0berr = Zinv1b0berr2;}
			if(i== 58) Zinv1b0berr = 0; if(i== 67) { Zinv1b0berr3 = Zinv1b0berr3*0.018/0.024; Zinv1b0berr = Zinv1b0berr3;}
			if(i== 73) Zinv1b0berr = 0; if(i== 83) { Zinv1b0berr6 = Zinv1b0berr6*0.143/0.168; Zinv1b0berr = Zinv1b0berr6;}
			if(i== 87) Zinv1b0berr = 0; if(i==100) { Zinv1b0berr2 = Zinv1b0berr2*0.060/0.050; Zinv1b0berr = Zinv1b0berr2;}
			if(i==102) Zinv1b0berr = 0; if(i==109) { Zinv1b0berr3 = Zinv1b0berr3*0.030/0.018; Zinv1b0berr = Zinv1b0berr3;}
			if(i==113) Zinv1b0berr = 0; if(i==118) { Zinv1b0berr6 = Zinv1b0berr6*0.201/0.143; Zinv1b0berr = Zinv1b0berr6;} if(i==121) Zinv1b0berr = 0;
			if(i==  1) LLerr = gRandom->Gaus(0,LLrelerr); if(i==  9) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 15) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i== 23) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 29) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 33) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i== 36) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 39) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 42) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i== 44) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 53) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 58) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i== 67) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 73) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 78) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i== 83) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 87) LLerr = gRandom->Gaus(0,LLrelerr); if(i== 91) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i== 94) LLerr = gRandom->Gaus(0,LLrelerr); if(i==100) LLerr = gRandom->Gaus(0,LLrelerr); if(i==102) LLerr = gRandom->Gaus(0,LLrelerr); 
			if(i==109) LLerr = gRandom->Gaus(0,LLrelerr); if(i==113) LLerr = gRandom->Gaus(0,LLrelerr); if(i==115) LLerr = gRandom->Gaus(0,LLrelerr); 

			double origuncorrsmearederr = gRandom->Gaus(0,origuncorrerr);
			double tn_on = gRandom->PoissonD(orig+origuncorrsmearederr+LLerr*LL->GetBinContent(i)+Zinv1b0berr*zinv->GetBinContent(i)+Zinv20errUsed*zinv->GetBinContent(i));
			double pullpseudodata = (tn_on - orig) / sqrt(pow(origerr,2) + tn_on);
			pulltemp->Fill(pullpseudodata);

		}
		//FILL HERE PULL MEAN HISTO, PULL RMS HISTO
		double mean = pulltemp->GetMean(); double meanerr = pulltemp->GetMeanError();
		double rms = pulltemp->GetRMS(); double rmserr = pulltemp->GetRMSError();
		if(trial%(nTrials/25)==0) cout << "mean " << mean << "+/-" << meanerr << ", rms " << rms << "+/-" << rmserr << endl;
		delete pulltemp;
		pullmean->Fill(mean);
		pullrms ->Fill(rms);
		if(0.26<mean) ++counter;//pull mean of observed experiment is 0.26
	}

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