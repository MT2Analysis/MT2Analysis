{
  TString SAMPLE="histoScan_T1tttt_v2_Moriond2_JES42V23_NLL_ISR";
  //  TString SYST = "jes";

  TString SYST2[3] = {"jes", "btag", "th"};

  for(int j=0; j<3; j++){
    TString SYST = SYST2[j];

  cout << "Sample: " << SAMPLE << endl;
  cout << "Syst: " << SYST << endl;
  cout << "Outfile: " << SAMPLE+"-"+SYST+"Unc.root" <<  endl;


  TFile *fNom;
  TFile *fUp;
  TFile *fDown;
  
  if(SYST=="btag"){
    fNom =  TFile::Open(SAMPLE+"_wBWeight.root");
    fUp = TFile::Open(SAMPLE+"_b_up_wBWeight.root");
    fDown = TFile::Open(SAMPLE+"_b_down_wBWeight.root");  
  }
  else if(SYST=="th"){
    fNom =  TFile::Open(SAMPLE+"_wBWeight.root");
    fUp = TFile::Open(SAMPLE+"_th_up_wBWeight.root");
    fDown = TFile::Open(SAMPLE+"_th_down_wBWeight.root");  
  }
  else if(SYST=="jes"){
    fNom =  TFile::Open(SAMPLE+"_wBWeight.root");
    fUp = TFile::Open(SAMPLE+"_jes_up_wBWeight.root");
    fDown = TFile::Open(SAMPLE+"_jes_down_wBWeight.root");  
  }



  TH2F * hNom, *hSUp, *hSDown;

  TH2F *hRes[200];
  TH2F *hUnc[200], *hUncAvg[200];
  TH2F *hUncDiff[200];

  TIter next(fNom->GetListOfKeys());
  TKey *key;
  int  count2D=0;
  while ((key=(TKey*)next())) {
    TString className =  key->GetClassName();
    if(className=="TH2F"){

      TString hname = key->GetName();
      if(!hname.Contains("NLO")) continue;
      //cout << "--------------------------" << key->GetName() << endl;


      hNom   = (TH2F*) fNom->Get( hname );
      hSUp   = (TH2F*) fUp->Get( hname );
      hSDown = (TH2F*) fDown->Get( hname );
      hUnc[count2D]   = (TH2F*) hNom->Clone(hname+"_Unc");
      hUncDiff[count2D]   = (TH2F*) hNom->Clone(hname+"_UncDiff");
      hUncAvg[count2D]   = (TH2F*) hNom->Clone(hname+"_UncAvg");


      for(int x=1; x<hNom->GetNbinsX(); x++){
	for(int y=1; y<hNom->GetNbinsY(); y++){
	  float nom = hNom->GetBinContent(x,y);
	  hUnc[count2D]->SetBinContent(x,y,0);
	  if(nom==0) continue;

	  float sup = hSUp->GetBinContent(x,y);
	  float sdo = hSDown->GetBinContent(x,y);
	  
	  //Theo perc error
	  float sup_p = fabs(sup - nom)/nom;
	  float sdo_p = fabs(sdo - nom)/nom;

	  if(hname=="MT2-HT_950to9999-MT2_150to200_NLO"){
	    /*	    cout << hNom->GetXaxis()->GetBinLowEdge(x) << " " <<  hNom->GetYaxis()->GetBinLowEdge(y)<< " " << nom << " " 
		 << " " << sup << " " << sdo
                 << " +-% " << hSUp->GetBinError(x,y)/sup  << " " << hSDown->GetBinError(x,y)/sdo
		 << " % " << sup_p << " " << sdo_p<< endl;*/
	  }
	
	  hUncDiff[count2D]->SetBinContent(x,y, sup_p - sdo_p);
	  if(sdo_p>sup_p) hUnc[count2D]->SetBinContent(x,y,sdo_p);
	  else  hUnc[count2D]->SetBinContent(x,y,sup_p);
	  hUncAvg[count2D]->SetBinContent(x,y,fabs(sup_p+sdo_p)/2.);
	}
      }

      count2D++;
    }
    else{
      cout <<"[WARNING] What is this? Skipping " << key->GetName()<<endl;
    }
  }

  
  TFile *fUnc = TFile::Open(SAMPLE+"-"+SYST+"Unc.root","RECREATE");
  for(int i=0;i<count2D; i++){
    hUnc[i]->Write();
    hUncAvg[i]->Write();
    hUncDiff[i]->Write();
  }
  fUnc->Write();
  fUnc->Close();
  }    

}
