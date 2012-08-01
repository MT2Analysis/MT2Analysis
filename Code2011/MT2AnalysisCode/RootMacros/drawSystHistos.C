{
  gROOT->ProcessLine(".x SetStyle_PRD.C");


  TString SYST="ISRUnc";
  TString SAMPLE="histoScan_T1_v1_Moriond2_JES42V23-"+SYST;
  TString MT2="MT2-";
  TString HT="HT_950";

    cout << "Sample: " << SAMPLE << endl;

  TFile * fNom =  TFile::Open(SAMPLE+".root");


  gStyle->SetPalette(1);
  TH2F * hNom[200];

  TIter next(fNom->GetListOfKeys());
  TKey *key;
  int  count2D=0;
  while ((key=(TKey*)next())) {
    TString className =  key->GetClassName();
    if(className=="TH2F"){

      TString hname = key->GetName();
      if(!hname.Contains(MT2)) continue;
      if(!hname.Contains(HT)) continue;
      if(hname.Contains("Avg")) continue;
      if(hname.Contains("Diff")) continue;


      cout << hname << endl;

      hNom[count2D]   = (TH2F*) fNom->Get( hname );
      count2D++;
    }
    else{
      cout <<"[WARNING] What is this? Skipping " << key->GetName()<<endl;
    }
  }

  cout << count2D << endl;
  
  TString cName = "histosUnc_"+SAMPLE+"-"+MT2+"-"+HT;;
  TCanvas *c = new TCanvas(cName, cName);

  c->Divide(3,2);
  for(int i=0; i< count2D; i++){
    c->cd(i+1);
    hNom[i]->Rebin2D(5,5);
    hNom[i]->GetXaxis()->SetRangeUser(0,1000);
    hNom[i]->GetYaxis()->SetRangeUser(0,1000);
    hNom[i]->GetZaxis()->SetRangeUser(0,0.5);
    cout << hNom[i]->GetTitle() << endl;
    hNom[i]->Draw("col");
  }

  //  fNom->Close();
  

}
