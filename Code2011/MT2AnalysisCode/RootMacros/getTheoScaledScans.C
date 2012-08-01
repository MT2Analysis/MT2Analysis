{

  TFile *fNom =  TFile::Open("histoScan_mSugra_tanb10_v3.root");
  TFile *fPDF = TFile::Open("histoScan_mSugra_tanb10_v4_NLOPdfUp.root");
  TFile *fScaleUp = TFile::Open("histoScan_mSugra_tanb10_v1_NLO2.root");
  TFile *fScaleDown = TFile::Open("histoScan_mSugra_tanb10_v1_NLO05.root");
  
  TH2F * hNom, *hPdf, *hSUp, *hSDown;

  TH2F *hRes[200];
  TH2F *hPdfUnc[200], *hScaleUnc[200];

  TH2F* hEvUp[200], *hEvDown[200];

  TIter next(fNom->GetListOfKeys());
  TKey *key;
  int  count2D=0;
  while ((key=(TKey*)next())) {
    TString className =  key->GetClassName();
    if(className=="TH2F"){

      TString hname = key->GetName();
      if(!hname.Contains("NLO")) continue;
      cout << "--------------------------" << key->GetName() << endl;


      hNom   = (TH2F*) fNom->Get( hname );
      hPdf   = (TH2F*) fPDF->Get( hname );
      hSUp   = (TH2F*) fScaleUp->Get( hname );
      hSDown = (TH2F*) fScaleDown->Get( hname );
      hPdfUnc[count2D]   = (TH2F*) hNom->Clone(hname+"_pdfUnc");
      hScaleUnc[count2D] = (TH2F*) hNom->Clone(hname+"_scaleUnc");

      hEvUp[count2D]   = (TH2F*) hNom->Clone(hname+"_EvUp");
      hEvDown[count2D]   = (TH2F*) hNom->Clone(hname+"_EvDown");

      for(int x=1; x<hNom->GetNbinsX(); x++){
	for(int y=1; y<hNom->GetNbinsY(); y++){
	  float nom = hNom->GetBinContent(x,y);
	  hScaleUnc[count2D]->SetBinContent(x,y,0);
	  if(nom==0) continue;

	  float pdf = hPdf->GetBinContent(x,y);
	  float sup = hSUp->GetBinContent(x,y);
	  float sdo = hSDown->GetBinContent(x,y);
	  
	  //Theo perc error
	  float pdf_p = fabs( pdf - nom)/nom;
	  float sup_p = fabs(sup - nom)/nom;
	  float sdo_p = fabs(sdo - nom)/nom;

	  float tot_plus = sqrt( sup_p*sup_p + pdf_p*pdf_p  );
	  float tot_minus = sqrt( sdo_p*sdo_p + pdf_p*pdf_p  );

	  if(hname=="MT2-HT_950to9999-MT2_500to9999_NLO"){
	    //cout << hNom->GetXaxis()->GetBinLowEdge(x) << " " <<  hNom->GetYaxis()->GetBinLowEdge(y)<< " " << nom << " " << pdf << " " << sup << " " << sdo << endl;
	    cout << hNom->GetXaxis()->GetBinLowEdge(x) << " " <<  hNom->GetYaxis()->GetBinLowEdge(y)<< " " << nom << " " << 
	      100.*pdf_p << " " << 100.*sup_p << " " << 100.*sdo_p << endl;

	  }
	  
	  if(sdo_p>sup_p) hScaleUnc[count2D]->SetBinContent(x,y,sdo_p);
	  else  hScaleUnc[count2D]->SetBinContent(x,y,sup_p);
	  hPdfUnc[count2D]->SetBinContent(x,y,pdf_p);  

	  hEvUp[count2D]->SetBinContent(x,y, nom*(1.+tot_plus));  
	  hEvDown[count2D]->SetBinContent(x,y, nom*(1.-tot_minus));  

	}
      }

      count2D++;
    }
    else{
      cout <<"[WARNING] What is this? Skipping " << key->GetName()<<endl;
    }
  }

  
  TFile *fUnc = TFile::Open("mSugra_uncertainties_v4.root","RECREATE");
  for(int i=0;i<count2D; i++){
    hScaleUnc[i]->Write();
    hPdfUnc[i]->Write();
  }
  fUnc->Write();
  fUnc->Close();
  
  
  TFile *fUp = TFile::Open("histoScan_mSugra_tanb10_v4_Plus1Sigma.root","RECREATE");
  for(int i=0;i<count2D; i++){
    TString hname = hEvUp[i]->GetName(   );
    hEvUp[i]->SetName( hname.ReplaceAll("_EvUp","")  );
    cout << hname.ReplaceAll("_EvUp","")  << endl;
    hEvUp[i]->Write();
  }
  fUp->Write();
  fUp->Close();


  TFile *fDown = TFile::Open("histoScan_mSugra_tanb10_v4_Minus1Sigma.root","RECREATE");
  for(int i=0;i<count2D; i++){
    TString hname = hEvDown[i]->GetName(   );
    hEvDown[i]->SetName( hname.ReplaceAll("_EvDown","")  );
    cout << hname.ReplaceAll("_EvDown","")  << endl;
    hEvDown[i]->Write();
  }
  fDown->Write();
  fDown->Close(); 
  

}
