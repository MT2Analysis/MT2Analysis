//called via root -l -b -q HO.C
//this is a very stupid and extrmely slow macro studying several pfmet and calomet inputs tested for RA2b's HO tagger.
//I hope that this code will be obsolete for 13 TeV.
//Macro should be easy to understand.
{

gROOT->cd();
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

TChain *c = new TChain("MassTree");
c->Add("~mmasciov/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/highHT/HT*.root");
c->Add("~mmasciov/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/highHT/JetHT*.root");
c->Add("~mmasciov/MT2Analysis/MT2trees/MT2_V02-03-02/20130318_8TeV/lowHT/MET*.root");

TH2D *calovspf = new TH2D("calovspf", "", 75, 0, 1500, 75, 0 ,1500);
TH2D *calomuvspf = new TH2D("calomuvspf", "", 75, 0, 1500, 75, 0 ,1500);
TH1D *pfdivcalo = new TH1D("pfdivcalo", "", 75,0,15);
TH1D *pfdivcalomu = new TH1D("pfdivcalomu", "", 75,0,15);

TH2D *j0calovspf = new TH2D("j0calovspf", "", 75, 0, 1500, 75, 0 ,1500);
TH2D *j0calomuvspf = new TH2D("j0calomuvspf", "", 75, 0, 1500, 75, 0 ,1500);
TH1D *j0pfdivcalo = new TH1D("j0pfdivcalo", "", 75,0,15);
TH1D *j0pfdivcalomu = new TH1D("j0pfdivcalomu", "", 75,0,15);

TH2D *MT2calovspf = new TH2D("MT2calovspf", "", 75, 0, 1500, 75, 0 ,1500);
TH2D *MT2calomuvspf = new TH2D("MT2calomuvspf", "", 75, 0, 1500, 75, 0 ,1500);
TH1D *MT2pfdivcalo = new TH1D("MT2pfdivcalo", "", 75,0,15);
TH1D *MT2pfdivcalomu = new TH1D("MT2pfdivcalomu", "", 75,0,15);

TH2D *MT2j0calovspf = new TH2D("MT2j0calovspf", "", 75, 0, 1500, 75, 0 ,1500);
TH2D *MT2j0calomuvspf = new TH2D("MT2j0calomuvspf", "", 75, 0, 1500, 75, 0 ,1500);
TH1D *MT2j0pfdivcalo = new TH1D("MT2j0pfdivcalo", "", 75,0,15);
TH1D *MT2j0pfdivcalomu = new TH1D("MT2j0pfdivcalomu", "", 75,0,15);

calovspf->GetXaxis()->SetTitle("PFMET [GeV]"); calovspf->GetYaxis()->SetTitle("CaloMET [GeV]");
calomuvspf->GetXaxis()->SetTitle("PFMET [GeV]"); calomuvspf->GetYaxis()->SetTitle("CaloMET (#mu-corr) [GeV]");
j0calovspf->GetXaxis()->SetTitle("PFMET [GeV]"); j0calovspf->GetYaxis()->SetTitle("CaloMET [GeV]");
j0calomuvspf->GetXaxis()->SetTitle("PFMET [GeV]"); j0calomuvspf->GetYaxis()->SetTitle("CaloMET (#mu-corr) [GeV]");
pfdivcalo->GetXaxis()->SetTitle("PFMET/CaloMET"); pfdivcalo->GetYaxis()->SetTitle("Events / 0.2");
pfdivcalomu->GetXaxis()->SetTitle("PFMET/CaloMET"); pfdivcalomu->GetYaxis()->SetTitle("Events / 0.2");
j0pfdivcalo->GetXaxis()->SetTitle("PFMET/CaloMET (#mu-corr)"); j0pfdivcalo->GetYaxis()->SetTitle("Events / 0.2");
j0pfdivcalomu->GetXaxis()->SetTitle("PFMET/CaloMET (#mu-corr)"); j0pfdivcalomu->GetYaxis()->SetTitle("Events / 0.2");

MT2calovspf->GetXaxis()->SetTitle("PFMET [GeV]"); MT2calovspf->GetYaxis()->SetTitle("CaloMET [GeV]");
MT2calomuvspf->GetXaxis()->SetTitle("PFMET [GeV]"); MT2calomuvspf->GetYaxis()->SetTitle("CaloMET (#mu-corr) [GeV]");
MT2j0calovspf->GetXaxis()->SetTitle("PFMET [GeV]"); MT2j0calovspf->GetYaxis()->SetTitle("CaloMET [GeV]");
MT2j0calomuvspf->GetXaxis()->SetTitle("PFMET [GeV]"); MT2j0calomuvspf->GetYaxis()->SetTitle("CaloMET (#mu-corr) [GeV]");
MT2pfdivcalo->GetXaxis()->SetTitle("PFMET/CaloMET"); MT2pfdivcalo->GetYaxis()->SetTitle("Events / 0.2");
MT2pfdivcalomu->GetXaxis()->SetTitle("PFMET/CaloMET"); MT2pfdivcalomu->GetYaxis()->SetTitle("Events / 0.2");
MT2j0pfdivcalo->GetXaxis()->SetTitle("PFMET/CaloMET (#mu-corr)"); MT2j0pfdivcalo->GetYaxis()->SetTitle("Events / 0.2");
MT2j0pfdivcalomu->GetXaxis()->SetTitle("PFMET/CaloMET (#mu-corr)"); MT2j0pfdivcalomu->GetYaxis()->SetTitle("Events / 0.2");

TCanvas *c1 = new TCanvas();

int nev = 0;

gROOT->cd();
nev = c->Draw("misc.CaloMETRaw:misc.MET>>+j0calovspf","((misc.HT>750&&misc.MET>30)||(misc.HT>450&&misc.HT<750&&misc.MET>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&Sum$(jet.isPFIDLoose&&jet.lv.Pt()>40&&abs(jet.lv.Eta())<0.5)>0","goff");
std::cout << "misc.CaloMETRaw:misc.MET>>+j0calovspf" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.CaloMETMuJesCorr:misc.MET>>+j0calomuvspf","((misc.HT>750&&misc.MET>30)||(misc.HT>450&&misc.HT<750&&misc.MET>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&Sum$(jet.isPFIDLoose&&jet.lv.Pt()>40&&abs(jet.lv.Eta())<0.5)>0","goff");
std::cout << "misc.CaloMETMuJesCorr:misc.MET>>+j0calomuvspf" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.CaloMETRaw:misc.MET>>+calovspf","((misc.HT>750&&misc.MET>30)||(misc.HT>450&&misc.HT<750&&misc.MET>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0","goff");
std::cout << "misc.CaloMETRaw:misc.MET>>+calovspf" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.CaloMETMuJesCorr:misc.MET>>+calomuvspf","((misc.HT>750&&misc.MET>30)||(misc.HT>450&&misc.HT<750&&misc.MET>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0","goff");
std::cout << "misc.CaloMETMuJesCorr:misc.MET>>+calomuvspf" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.CaloMETRaw:misc.MET>>+MT2j0calovspf","((misc.HT>750&&misc.MET>30&&misc.MT2>125)||(misc.HT>450&&misc.HT<750&&misc.MET>200&&misc.MT2>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&Sum$(jet.isPFIDLoose&&jet.lv.Pt()>40&&abs(jet.lv.Eta())<0.5)>0","goff");
std::cout << "misc.CaloMETRaw:misc.MET>>+MT2j0calovspf" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.CaloMETMuJesCorr:misc.MET>>+MT2j0calomuvspf","((misc.HT>750&&misc.MET>30&&misc.MT2>125)||(misc.HT>450&&misc.HT<750&&misc.MET>200&&misc.MT2>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&Sum$(jet.isPFIDLoose&&jet.lv.Pt()>40&&abs(jet.lv.Eta())<0.5)>0","goff");
std::cout << "misc.CaloMETMuJesCorr:misc.MET>>+MT2j0calomuvspf" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.CaloMETRaw:misc.MET>>+MT2calovspf","((misc.HT>750&&misc.MET>30&&misc.MT2>125)||(misc.HT>450&&misc.HT<750&&misc.MET>200&&misc.MT2>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0","goff");
std::cout << "misc.CaloMETRaw:misc.MET>>+MT2calovspf" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.CaloMETMuJesCorr:misc.MET>>+MT2calomuvspf","((misc.HT>750&&misc.MET>30&&misc.MT2>125)||(misc.HT>450&&misc.HT<750&&misc.MET>200&&misc.MT2>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0","goff");
std::cout << "misc.CaloMETMuJesCorr:misc.MET>>+MT2calomuvspf" << " --> events: " << nev << std::endl;

gROOT->cd();
nev = c->Draw("misc.MET/misc.CaloMETRaw>>+j0pfdivcalo","((misc.HT>750&&misc.MET>30)||(misc.HT>450&&misc.HT<750&&misc.MET>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&Sum$(jet.isPFIDLoose&&jet.lv.Pt()>40&&abs(jet.lv.Eta())<0.5)>0","goff");
std::cout << "misc.MET/misc.CaloMETRaw>>+j0pfdivcalo" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.MET/misc.CaloMETMuJesCorr>>+j0pfdivcalomu","((misc.HT>750&&misc.MET>30)||(misc.HT>450&&misc.HT<750&&misc.MET>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&Sum$(jet.isPFIDLoose&&jet.lv.Pt()>40&&abs(jet.lv.Eta())<0.5)>0","goff");
std::cout << "misc.MET/misc.CaloMETMuJesCorr>>+j0pfdivcalomu" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.MET/misc.CaloMETRaw>>+pfdivcalo","((misc.HT>750&&misc.MET>30)||(misc.HT>450&&misc.HT<750&&misc.MET>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0","goff");
std::cout << "misc.MET/misc.CaloMETRaw>>+pfdivcalo" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.MET/misc.CaloMETMuJesCorr>>+pfdivcalomu","((misc.HT>750&&misc.MET>30)||(misc.HT>450&&misc.HT<750&&misc.MET>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0","goff");
std::cout << "misc.MET/misc.CaloMETMuJesCorr>>+pfdivcalomu" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.MET/misc.CaloMETRaw>>+MT2j0pfdivcalo","((misc.HT>750&&misc.MET>30&&misc.MT2>125)||(misc.HT>450&&misc.HT<750&&misc.MET>200&&misc.MT2>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&Sum$(jet.isPFIDLoose&&jet.lv.Pt()>40&&abs(jet.lv.Eta())<0.5)>0","goff");
std::cout << "misc.MET/misc.CaloMETRaw>>+MT2j0pfdivcalo" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.MET/misc.CaloMETMuJesCorr>>+MT2j0pfdivcalomu","((misc.HT>750&&misc.MET>30&&misc.MT2>125)||(misc.HT>450&&misc.HT<750&&misc.MET>200&&misc.MT2>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0&&Sum$(jet.isPFIDLoose&&jet.lv.Pt()>40&&abs(jet.lv.Eta())<0.5)>0","goff");
std::cout << "misc.MET/misc.CaloMETMuJesCorr>>+MT2j0pfdivcalomu" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.MET/misc.CaloMETRaw>>+MT2pfdivcalo","((misc.HT>750&&misc.MET>30&&misc.MT2>125)||(misc.HT>450&&misc.HT<750&&misc.MET>200&&misc.MT2>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0","goff");
std::cout << "misc.MET/misc.CaloMETRaw>>+MT2pfdivcalo" << " --> events: " << nev << std::endl;
gROOT->cd();
nev = c->Draw("misc.MET/misc.CaloMETMuJesCorr>>+MT2pfdivcalomu","((misc.HT>750&&misc.MET>30&&misc.MT2>125)||(misc.HT>450&&misc.HT<750&&misc.MET>200&&misc.MT2>200))&&misc.Jet0Pass==1&&misc.Jet1Pass==1&&misc.Vectorsumpt<70&&misc.PassJet40ID&&(NEles+NMuons+NTausIDLoose3Hits)==0&&misc.MinMetJetDPhi4Pt40>0.3&&misc.HBHENoiseFlag==0&&misc.CSCTightHaloIDFlag==0&&misc.trackingFailureFlag==0&&misc.eeBadScFlag==0&&misc.EcalDeadCellTriggerPrimitiveFlag==0&&misc.TrackingManyStripClusFlag==0&&misc.TrackingTooManyStripClusFlag==0&&misc.TrackingLogErrorTooManyClustersFlag==0&&misc.CrazyHCAL==0","goff");
std::cout << "misc.MET/misc.CaloMETMuJesCorr>>+MT2pfdivcalomu" << " --> events: " << nev << std::endl;


int nbinsx = calovspf->GetNbinsX();
pfdivcalo->SetBinContent(nbinsx, pfdivcalo->GetBinContent(nbinsx)+pfdivcalo->GetBinContent(nbinsx+1));
pfdivcalomu->SetBinContent(nbinsx, pfdivcalomu->GetBinContent(nbinsx)+pfdivcalomu->GetBinContent(nbinsx+1));
j0pfdivcalo->SetBinContent(nbinsx, j0pfdivcalo->GetBinContent(nbinsx)+j0pfdivcalo->GetBinContent(nbinsx+1));
j0pfdivcalomu->SetBinContent(nbinsx, j0pfdivcalomu->GetBinContent(nbinsx)+j0pfdivcalomu->GetBinContent(nbinsx+1));
MT2pfdivcalo->SetBinContent(nbinsx, MT2pfdivcalo->GetBinContent(nbinsx)+MT2pfdivcalo->GetBinContent(nbinsx+1));
MT2pfdivcalomu->SetBinContent(nbinsx, MT2pfdivcalomu->GetBinContent(nbinsx)+MT2pfdivcalomu->GetBinContent(nbinsx+1));
MT2j0pfdivcalo->SetBinContent(nbinsx, MT2j0pfdivcalo->GetBinContent(nbinsx)+MT2j0pfdivcalo->GetBinContent(nbinsx+1));
MT2j0pfdivcalomu->SetBinContent(nbinsx, MT2j0pfdivcalomu->GetBinContent(nbinsx)+MT2j0pfdivcalomu->GetBinContent(nbinsx+1));
for(int i=1; i<=nbinsx;++i){
calovspf->SetBinContent(i, nbinsx, calovspf->GetBinContent(i, nbinsx)+calovspf->GetBinContent(i, nbinsx+1));
calovspf->SetBinContent(nbinsx, i, calovspf->GetBinContent(nbinsx, i)+calovspf->GetBinContent(nbinsx+1, i));
calomuvspf->SetBinContent(i, nbinsx, calomuvspf->GetBinContent(i, nbinsx)+calomuvspf->GetBinContent(i, nbinsx+1));
calomuvspf->SetBinContent(nbinsx, i, calomuvspf->GetBinContent(nbinsx, i)+calomuvspf->GetBinContent(nbinsx+1, i));
j0calovspf->SetBinContent(i, nbinsx, j0calovspf->GetBinContent(i, nbinsx)+j0calovspf->GetBinContent(i, nbinsx+1));
j0calovspf->SetBinContent(nbinsx, i, j0calovspf->GetBinContent(nbinsx, i)+j0calovspf->GetBinContent(nbinsx+1, i));
j0calomuvspf->SetBinContent(i, nbinsx, j0calomuvspf->GetBinContent(i, nbinsx)+j0calomuvspf->GetBinContent(i, nbinsx+1));
j0calomuvspf->SetBinContent(nbinsx, i, j0calomuvspf->GetBinContent(nbinsx, i)+j0calomuvspf->GetBinContent(nbinsx+1, i));
MT2calovspf->SetBinContent(i, nbinsx, MT2calovspf->GetBinContent(i, nbinsx)+MT2calovspf->GetBinContent(i, nbinsx+1));
MT2calovspf->SetBinContent(nbinsx, i, MT2calovspf->GetBinContent(nbinsx, i)+MT2calovspf->GetBinContent(nbinsx+1, i));
MT2calomuvspf->SetBinContent(i, nbinsx, MT2calomuvspf->GetBinContent(i, nbinsx)+MT2calomuvspf->GetBinContent(i, nbinsx+1));
MT2calomuvspf->SetBinContent(nbinsx, i, MT2calomuvspf->GetBinContent(nbinsx, i)+MT2calomuvspf->GetBinContent(nbinsx+1, i));
MT2j0calovspf->SetBinContent(i, nbinsx, MT2j0calovspf->GetBinContent(i, nbinsx)+MT2j0calovspf->GetBinContent(i, nbinsx+1));
MT2j0calovspf->SetBinContent(nbinsx, i, MT2j0calovspf->GetBinContent(nbinsx, i)+MT2j0calovspf->GetBinContent(nbinsx+1, i));
MT2j0calomuvspf->SetBinContent(i, nbinsx, MT2j0calomuvspf->GetBinContent(i, nbinsx)+MT2j0calomuvspf->GetBinContent(i, nbinsx+1));
MT2j0calomuvspf->SetBinContent(nbinsx, i, MT2j0calomuvspf->GetBinContent(nbinsx, i)+MT2j0calomuvspf->GetBinContent(nbinsx+1, i));
}

TFile *newfile = new TFile("CaloMET/HO.root","RECREATE");
newfile->cd();
calovspf->Write();
calomuvspf->Write();
pfdivcalo->Write();
pfdivcalomu->Write();
j0calovspf->Write();
j0calomuvspf->Write();
j0pfdivcalo->Write();
j0pfdivcalomu->Write();
MT2calovspf->Write();
MT2calomuvspf->Write();
MT2pfdivcalo->Write();
MT2pfdivcalomu->Write();
MT2j0calovspf->Write();
MT2j0calomuvspf->Write();
MT2j0pfdivcalo->Write();
MT2j0pfdivcalomu->Write();
newfile->Close();

string savename;
c1->Clear();
c1->cd();
calovspf->Draw("COLZ");
savename = "CaloMET/calovspf.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
calomuvspf->Draw("COLZ");
savename = "CaloMET/calomuvspf.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
j0calovspf->Draw("COLZ");
savename = "CaloMET/j0calovspf.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
j0calomuvspf->Draw("COLZ");
savename = "CaloMET/j0calomuvspf.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
MT2calovspf->Draw("COLZ");
savename = "CaloMET/MT2calovspf.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
MT2calomuvspf->Draw("COLZ");
savename = "CaloMET/MT2calomuvspf.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
MT2j0calovspf->Draw("COLZ");
savename = "CaloMET/MT2j0calovspf.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
MT2j0calomuvspf->Draw("COLZ");
savename = "CaloMET/MT2j0calomuvspf.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
pfdivcalo->Draw("");
savename = "CaloMET/pfdivcalo.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
pfdivcalomu->Draw("");
savename = "CaloMET/pfdivcalomu.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
j0pfdivcalo->Draw("");
savename = "CaloMET/j0pfdivcalo.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
j0pfdivcalomu->Draw("");
savename = "CaloMET/j0pfdivcalomu.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
MT2pfdivcalo->Draw("");
savename = "CaloMET/MT2pfdivcalo.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
MT2pfdivcalomu->Draw("");
savename = "CaloMET/MT2pfdivcalomu.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
MT2j0pfdivcalo->Draw("");
savename = "CaloMET/MT2j0pfdivcalo.eps";
c1->SaveAs(savename.c_str());

c1->Clear();
c1->cd();
MT2j0pfdivcalomu->Draw("");
savename = "CaloMET/MT2j0pfdivcalomu.eps";
c1->SaveAs(savename.c_str());

}