{
TChain *t50  = new TChain("MassTree");
TChain *t100 = new TChain("MassTree");
TChain *t200 = new TChain("MassTree");
TChain *tinc = new TChain("MassTree");

t50->Add("~/MT2Analysis/MT2trees/MT2_V01-00-00/20110919_MC_nocuts_2/ZJetsToNuNu_50_HT_100_Summer11.root");
t100->Add("~/MT2Analysis/MT2trees/MT2_V01-00-00/20110919_MC_nocuts_2/ZJetsToNuNu_100_HT_200_Summer11.root");
t200->Add("~/MT2Analysis/MT2trees/MT2_V01-00-00/20110919_MC_nocuts_2/ZJetsToNuNu_200_HT_inf_Summer11.root");
tinc->Add("~/MT2Analysis/MT2trees/MT2_V01-00-00/20110919_MC_nocuts_2/ZinvisibleJets_7TeV-madgraph.root");

cout << t50->GetEntries() << " " << t100->GetEntries() << " " << t200->GetEntries() << endl;

//t50 ->Draw("GenJetHT()>>h50(200, 0, 2000)");
//cout << " h50 done" << endl;
//t100->Draw("GenJetHT()>>h100(200, 0, 2000)");
//cout << " h100 done" << endl;
//t200->Draw("GenJetHT()>>h200(200, 0, 2000)");
//cout << " h200 done" << endl;

t50 ->Draw("misc.HT>>h50(200, 0, 2000)");
cout << " h50 done" << endl;
t100->Draw("misc.HT>>h100(200, 0, 2000)");
cout << " h100 done" << endl;
t200->Draw("misc.HT>>h200(200, 0, 2000)");
cout << " h200 done" << endl;
tinc->Draw("misc.HT>>hinc(200, 0, 2000)");
cout << " hinc done" << endl;


Int_t nice_blue1  = TColor::GetColor("#0489B1");
Int_t nice_blue2  = TColor::GetColor("#0174DF");
Int_t nice_blue3  = TColor::GetColor("#013ADF");

h50 ->SetLineColor(nice_blue1);
h100->SetLineColor(nice_blue2);
h200->SetLineColor(nice_blue3);
hinc->SetLineColor(kMagenta);

h50 ->SetLineWidth(2);
h100->SetLineWidth(2);
h200->SetLineWidth(2);
hinc->SetLineWidth(2);

// ---------------
h50 ->Scale(309.5/8163840.);
h100->Scale(125.2/3837779.);
h200->Scale(28.64/3067017.);
hinc->Scale(4500./2106977);
//hinc->Scale(5760./2106977);

TH1F* hcombined = new TH1F("hcombined", "",200, 0, 2000);
hcombined->Add(h50);
hcombined->Add(h100);
hcombined->Add(h200);

h50 ->Draw();
h100->Draw("same");
h200->Draw("same");
hinc->Draw("same");
hcombined->Draw("same");

}
