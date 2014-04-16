//call macro via root -l GammaVsZnunuStudiesPlots.C
//this macro plots the gamma/Z(nunu) ratios (and distributions)
//this code is a modified clone of GammaVsZllStudiesRatios.C - see that macro first!
void GammaVsZnunuStudiesPlots(){

gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

int minvpt      = 20;
int minHT       = 450;
int maxHT       = -1;
int njets       = -2;
int nbjets       = 0;

double binlow      = -1;//define x-axis range (lower border), default = -1 (use border of histogram


TString outputdir  = "GammaVsZnunuStudies/Plots/";
TString inputdir   = "GammaVsZnunuStudies/";

TString varname = "MT2";//define variable name here
                        varname += "_SF_";
if(njets>=0)            varname += TString::Format("%dj"  ,   abs(njets));
else if(njets!=-10)     varname += TString::Format("ge%dj",   abs(njets));
if(njets!=-10)          varname += "_";
if(nbjets>=0)           varname += TString::Format("%db"  ,  abs(nbjets));
else if(nbjets!=-10)    varname += TString::Format("ge%db",  abs(nbjets));
if(nbjets!=-10)         varname += "_";
if(minvpt>0)            varname += TString::Format("VPtge%d",(int)minvpt);
if(minvpt>0)            varname += "_";
if(minHT >0)            varname += TString::Format( "HTge%d",(int)minHT );
if(minHT >0&&maxHT >0)  varname += TString::Format(   "le%d",(int)maxHT );
else if(     maxHT >0)  varname += TString::Format( "HTle%d",(int)maxHT );
if(minHT >0||maxHT >0)  varname += "_";
                        varname += "all";

TString inputfile = inputdir + varname + ".root";
TString hZname    = varname + "_Z";
TString hGname    = varname + "_G";
TString ratname   = varname + "_ZG_MCratio";

TFile *f = TFile::Open(inputfile.Data());
f->cd();
TH1D *Z  = (TH1D*)f->Get( hZname.Data());
TH1D *G  = (TH1D*)f->Get( hGname.Data());
TH1D *ZG = (TH1D*)f->Get(ratname.Data());

TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",334,93,600,600);
gStyle->SetOptFit(1);           gStyle->SetOptStat(0);           gStyle->SetOptTitle(0);
Canvas_1->Range(-82.79222,-0.9532603,793.8312,3.394941);         Canvas_1->SetTickx(1);          Canvas_1->SetTicky(1);
Canvas_1->SetFillColor(0);      Canvas_1->SetBorderMode(0);      Canvas_1->SetBorderSize(2);
Canvas_1->SetRightMargin(0.05); Canvas_1->SetTopMargin(0.07);    Canvas_1->SetBottomMargin(0.15);Canvas_1->SetLeftMargin(0.18);
Canvas_1->SetFrameFillStyle(0); Canvas_1->SetFrameBorderMode(0); Canvas_1->SetFrameBorderMode(0);

   double lowbin = (binlow>0)?binlow:Z->GetBinLowEdge(1);

   TString xtitle = Z->GetXaxis()->GetTitle();
   TString ytitle = Z->GetYaxis()->GetTitle();
   int binwidth = (int)Z->GetBinWidth(Z->GetNbinsX() );
   stringstream yTitle;
		yTitle << "Events / ";
		yTitle << binwidth;
		yTitle << " GeV";
   stringstream yTitle2;
		yTitle2 << "#frac{Z(#nu#bar{#nu})}{#gamma} / ";
		yTitle2 << binwidth;
		yTitle2 << " GeV";
   TH1D *haxis = new TH1D("haxis", "", Z->GetNbinsX(), lowbin, Z->GetBinLowEdge(Z->GetNbinsX() )+Z->GetBinWidth(Z->GetNbinsX() ) );
   haxis->GetXaxis()->SetTitle(xtitle.Data());
   haxis->GetYaxis()->SetTitle(yTitle.str().c_str());
haxis->GetXaxis()->SetLabelFont(42);
haxis->GetXaxis()->SetLabelSize(0.05);
haxis->GetXaxis()->SetLabelOffset(0.007);
haxis->GetXaxis()->SetTitleFont(42);
haxis->GetXaxis()->SetTitleSize(0.06);
haxis->GetXaxis()->SetTitleOffset(0.9);
haxis->GetYaxis()->SetLabelFont(42);
haxis->GetYaxis()->SetLabelSize(0.05);
haxis->GetYaxis()->SetLabelOffset(0.007);
haxis->GetYaxis()->SetTitleFont(42);
haxis->GetYaxis()->SetTitleSize(0.06);
haxis->GetYaxis()->SetTitleOffset(1.25);

   double max1 = G->GetMaximum();
   double max2 = Z->GetMaximum();
   double max = (max1>max2)?max1:max2; max = 2.5*max;
   double min1 = G->GetMinimum();
   double min2 = Z->GetMinimum();
   double min = (min1<min2)?min1:min2; min = min/1.5;
   if(min<=0.) min = 0.02;
   haxis->SetMaximum(max);
   haxis->SetMinimum(min);
   Canvas_1->cd();
   haxis->Draw();
   Z->Draw("same");
   G->Draw("same");

TLegend *leg = new TLegend(0.6,0.68,0.8,0.88,NULL,"brNDC");
leg->SetBorderSize(0); leg->SetTextFont(42); leg->SetTextSize(0.04181185); leg->SetLineColor(1);
leg->SetLineStyle(1);  leg->SetLineWidth(2); leg->SetFillColor(0);         leg->SetFillStyle(1001);
   leg->AddEntry(G,"#gamma + jets","lp");
   leg->AddEntry(Z,"Z(#nu#bar{#nu}) + jets","lp");
   TLatex *   tex = new TLatex(0.328859,0.95,"CMS Simulation, #sqrt{s} = 8 TeV");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04181185);
   tex->SetLineWidth(2);
   tex->Draw();
    TString text = "";
    if(njets== 2) text += "2 jets";
    if(njets==-2) text += "#geq 2 jets";
    if(njets==35) text += "3-5 jets";
    if(njets==-6) text += "#geq 6 jets";
    if(nbjets==0) text += ", 0 b-jets";
   TLatex *   tex2 = new TLatex(0.9,0.803987,text.Data());
   tex2->SetNDC();
   tex2->SetTextAlign(31);
   tex2->SetTextFont(42);
   tex2->SetTextSize(0.04181185);
   tex2->SetLineWidth(2);
   tex2->Draw();
    if(minHT==1200)             text = "H_{T} #geq 1200 GeV";
    if(minHT==750&&maxHT==1200) text = "750 GeV #leq H_{T} < 1200 GeV";
    else if(minHT==750)         text = "H_{T} #geq 750 GeV";
    if(minHT==450&&maxHT==750)  text = "450 GeV #leq H_{T} < 750 GeV";
    else if(minHT==450)         text = "H_{T} #geq 450 GeV";
   TLatex *   tex3 = new TLatex(0.9,0.86711,text.Data());
   tex3->SetNDC();
   tex3->SetTextAlign(31);
   tex3->SetTextFont(42);
   tex3->SetTextSize(0.04181185);
   tex3->SetLineWidth(2);
   tex3->Draw();
    leg->Draw();
   Canvas_1->Modified();
   Canvas_1->cd();
   TString outnamepdf = outputdir + varname + ".pdf";
   TString outnameeps = outputdir + varname + ".eps";
   Canvas_1->SaveAs(outnameeps.Data());
   Canvas_1->SaveAs(outnamepdf.Data());

TCanvas *Canvas_2 = new TCanvas("Canvas_2", "Canvas_2",334,93,600,600);
gStyle->SetOptFit(1);           gStyle->SetOptStat(0);           gStyle->SetOptTitle(0);
Canvas_2->Range(-165.9091,-0.1923077,1061.364,1.089744);         Canvas_2->SetTickx(1);          Canvas_2->SetTicky(1);
Canvas_2->SetFillColor(0);      Canvas_2->SetBorderMode(0);      Canvas_2->SetBorderSize(2);
Canvas_2->SetRightMargin(0.05); Canvas_2->SetTopMargin(0.07);    Canvas_2->SetBottomMargin(0.15);Canvas_2->SetLeftMargin(0.18);
Canvas_2->SetFrameFillStyle(0); Canvas_2->SetFrameBorderMode(0); Canvas_2->SetFrameBorderMode(0);

   lowbin = (binlow>0)?binlow:ZG->GetBinLowEdge(1);
   TH1D *haxisZG = new TH1D("haxisZG", "", ZG->GetNbinsX(), lowbin, ZG->GetBinLowEdge(ZG->GetNbinsX() )+ZG->GetBinWidth(ZG->GetNbinsX() ) );
   haxisZG->GetXaxis()->SetTitle(xtitle);
   haxisZG->GetYaxis()->SetTitle(yTitle2.str().c_str());
haxis->GetXaxis()->SetLabelFont(42);
haxis->GetXaxis()->SetLabelSize(0.05);
haxis->GetXaxis()->SetLabelOffset(0.007);
haxis->GetXaxis()->SetTitleFont(42);
haxis->GetXaxis()->SetTitleSize(0.06);
haxis->GetXaxis()->SetTitleOffset(0.9);
haxis->GetYaxis()->SetLabelFont(42);
haxis->GetYaxis()->SetLabelSize(0.05);
haxis->GetYaxis()->SetLabelOffset(0.007);
haxis->GetYaxis()->SetTitleFont(42);
haxis->GetYaxis()->SetTitleSize(0.06);
haxis->GetYaxis()->SetTitleOffset(1.25);
   if(njets==-6) haxisZG->SetMaximum(2.);
   else          haxisZG->SetMaximum(1.);
   haxisZG->SetMinimum(0.);
   Canvas_2->cd();
   haxisZG->Draw("");
   ZG->Draw("same");
   tex->Draw();
   tex2->Draw();
   tex3->Draw();
   Canvas_2->Modified();
   Canvas_2->cd();
   outnamepdf = outputdir + varname + "_ratio.pdf";
   outnameeps = outputdir + varname + "_ratio.eps";
   Canvas_2->SaveAs(outnameeps.Data());
   Canvas_2->SaveAs(outnamepdf.Data());
}