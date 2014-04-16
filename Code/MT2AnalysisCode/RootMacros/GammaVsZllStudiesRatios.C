{//call this macro via root -l GammaVsZllStudiesRatios.C

//cuts and names equivalent to GammaVsZllStudies.C
//as you might also build ratio between Z events (like Z(1b)/Z(0b) ) there are two selections here
double minVPt1 = 20;
double minVPt2 = 20;
double minHT1  = 450;
double minHT2  = 450;
double maxHT1  = -1;
double maxHT2  = -1;
int    njets1  = -2;
int    njets2  = -2;
int    nbjets1 = 1;
int    nbjets2 = 0;  
int    HTselection = -1;//automatically set correctly depending on cuts
TString varname      = "MT2";//VPt //MT2 //HT //NJets //NBJets
    

bool   removebkg   = true;//remove background from selection using simulation (default = true)
double  rescaleZ   = -1.;//negative: no rescaling, otherwise rescale Z(ll) signal by this factor
double  rescaleG   = -1.;//negative: no rescaling, otherwise rescale 1 photon signal by this factor
double  rescaleZbg = -1.;//negative: no rescaling, otherwise rescale Z(ll) background by this factor
double  rescaleGbg = -1.;//negative: no rescaling, otherwise rescale 1 photon background by this factor
double  cutoff     = 0.01;//default 0.01 - kills bins that have negative entries due to background subtraction (expected background > observed data)

bool printDistZvsZorGvsG    = true; //default = true,  prints the distributions of Z(1st selection) vs. Z(2nd selection) (and same for photons)
bool printDistZvsG          = false;//default = true,  prints the distributions of Z vs. photons (for both selections)
bool printRatiosDataVsMC    = false;//default = false, prints both data/MC ratios (of the input distributions, i.e. this does not mean data/MC ratio of ratios like Z/G)
bool printRatiosZVsG        = false;//default = true,  prints ratio of Z/G for both selections (both in data and MC)
bool printRatiosZvsZdivGvsG = false;//default = false, prints double ratio [Z(1st sel.)/Z(2nd sel.)] / [G(1st sel.)/G(2nd sel.)] both in data and MC
bool printRatiosZvsZorGvsG  = true; //default = true,  prints ratio Z(1st sel.)/Z(2nd sel.) and ratio G(1st sel.)/G(2nd sel.) both in data and MC
bool dataonly               = false;//in distributions do not overlay data with MC (for printDistZvsZorGvsG, printDistZvsG)
bool dataandMC              = true;//in distributions overlay data with MC (for printDistZvsZorGvsG, printDistZvsG)
bool saveC      = false;//save canvases as root .C macro
bool saveeps    = false;//saves canvases as eps files
bool savepng    = false;//saves canvases as png files
bool makeclones = true;//create clones of canvases --> allows to modify each canvas afterwards

//I used these fits to get the Z(1b)/Z(0b) scaling factors and uncertainties on them
bool dofit      = false;//default = false, fit the several ratios with a constant or 2 parameters to check stability of the ratios
bool fit2par    = true;//default = false: false means fit constant, true, means fit a linear function (needed for NJets fit only)
int  fit300     = 0;//for MT2 and VPt ratio of Z/G we state that the ratio is constant for MT2/VPt > 300 GeV. This int defines the lower border of the fit range, defaults = 0, or 300
bool fitMC      = false;//fit the MC ratio on top of the data ratio

//x-axis range for plotting, if set to -1 use the borders of the histogram
double lowbin   = -1.;
double upbin    = -1.;
//rebin the histogram (needed if there are too low statistics, default = -1 (means no rebinning)
int    rebin    = -1.;
    
if(     (int)minHT1 >=750.) HTselection = 3;
else if((int)minVPt1>=180.) HTselection = 2;
else                        HTselection = 1;
    
//load the file names
TString fname1 = varname + "_SF_";
if(njets1>=0)            fname1 += TString::Format("%dj"  ,   abs(njets1));
else if(njets1!=-10)     fname1 += TString::Format("ge%dj",   abs(njets1));
if(njets1!=-10)          fname1 += "_";
if(nbjets1>=0)           fname1 += TString::Format("%db"  ,  abs(nbjets1));
else if(nbjets1!=-10)    fname1 += TString::Format("ge%db",  abs(nbjets1));
if(nbjets1!=-10)         fname1 += "_";
if(minVPt1>0)            fname1 += TString::Format("VPtge%d",(int)minVPt1);
if(minVPt1>0)            fname1 += "_";
if(minHT1 >0)            fname1 += TString::Format( "HTge%d",(int)minHT1);
if(minHT1 >0&&maxHT1 >0) fname1 += TString::Format(   "le%d",(int)maxHT1);
else if(      maxHT1 >0) fname1 += TString::Format( "HTle%d",(int)maxHT1);
if(minHT1 >0||maxHT1 >0) fname1 += "_";
if(HTselection==1)       fname1 += "all";
if(HTselection==2)       fname1 += "highVPt";
if(HTselection==3)       fname1 += "highHT";

TString fname2 = varname + "_SF_";
if(njets2>=0)            fname2 += TString::Format("%dj"  ,   abs(njets2));
else if(njets2!=-10)     fname2 += TString::Format("ge%dj",   abs(njets2));
if(njets2!=-10)          fname2 += "_";
if(nbjets2>=0)           fname2 += TString::Format("%db"  ,  abs(nbjets2));
else if(nbjets2!=-10)    fname2 += TString::Format("ge%db",  abs(nbjets2));
if(nbjets2!=-10)         fname2 += "_";
if(minVPt2>0)            fname2 += TString::Format("VPtge%d",(int)minVPt2);
if(minVPt2>0)            fname2 += "_";
if(minHT2 >0)            fname2 += TString::Format( "HTge%d",(int)minHT2);
if(minHT2 >0&&maxHT2 >0) fname2 += TString::Format(   "le%d",(int)maxHT2);
else if(      maxHT2 >0) fname2 += TString::Format( "HTle%d",(int)maxHT2);
if(minHT2 >0||maxHT2 >0) fname2 += "_";
if(HTselection==1)       fname2 += "all";
if(HTselection==2)       fname2 += "highVPt";
if(HTselection==3)       fname2 += "highHT";
 
//defines the names of the output file
TString  tmp_oname = varname + "_SF_";
TString diff_oname = "Ratio_";
if(njets1==njets2){
if(njets1>=0)            tmp_oname += TString::Format("%dj"  ,   abs(njets1));
else if(njets1!=-10)     tmp_oname += TString::Format("ge%dj",   abs(njets1));
if(njets1!=-10)          tmp_oname += "_";
} else{
TString j1="";
if(njets1>=0)            j1 = TString::Format("%dj"  ,   abs(njets1));
else if(njets1!=-10)     j1 = TString::Format("ge%dj",   abs(njets1));
TString j2="";
if(njets2>=0)            j2 = TString::Format("%dj"  ,   abs(njets2));
else if(njets2!=-10)     j2 = TString::Format("ge%dj",   abs(njets2));
diff_oname += j1 + "vs" + j2 + "_";
}
if(nbjets1==nbjets2){
if(nbjets1>=0)           tmp_oname += TString::Format("%db"  ,  abs(nbjets1));
else if(nbjets1!=-10)    tmp_oname += TString::Format("ge%db",  abs(nbjets1));
if(nbjets1!=-10)         tmp_oname += "_";
} else{
TString b1 = "";
if(nbjets1>=0)           b1 = TString::Format("%db"  ,  abs(nbjets1));
else if(nbjets1!=-10)    b1 = TString::Format("ge%db",  abs(nbjets1));
TString b2 = "";
if(nbjets2>=0)           b2 = TString::Format("%db"  ,  abs(nbjets2));
else if(nbjets2!=-10)    b2 = TString::Format("ge%db",  abs(nbjets2));
diff_oname += b1 + "vs" + b2 + "_";
}
if(minVPt1==minVPt2){
if(minVPt1>0)            tmp_oname += TString::Format("VPtge%d",(int)minVPt1);
if(minVPt1>0)            tmp_oname += "_";
} else{
TString v1 = "";
if(minVPt1>0)            v1 = TString::Format("VPtge%d",(int)minVPt1);
TString v2 = "";
if(minVPt2>0)            v2 = TString::Format("VPtge%d",(int)minVPt2);
diff_oname += v1 + "vs" + v2 + "_";
}
if(maxHT1==maxHT2&&minHT1==minHT2){
if(minHT1 >0)            tmp_oname += TString::Format( "HTge%d",(int)minHT1);
if(minHT1 >0&&maxHT1 >0) tmp_oname += TString::Format(   "le%d",(int)maxHT1);
else if(      maxHT1 >0) tmp_oname += TString::Format( "HTle%d",(int)maxHT1);
if(minHT1 >0||maxHT1 >0) tmp_oname += "_";
} else {
TString h1 = "";
if(minHT1 >0)            h1 = TString::Format( "HTge%d",(int)minHT1);
if(minHT1 >0&&maxHT1 >0) h1 = TString::Format(   "le%d",(int)maxHT1);
else if(      maxHT1 >0) h1 = TString::Format( "HTle%d",(int)maxHT1);
TString h2 = "";
if(minHT2 >0)            h2 = TString::Format( "HTge%d",(int)minHT2);
if(minHT2 >0&&maxHT2 >0) h2 = TString::Format(   "le%d",(int)maxHT2);
else if(      maxHT2 >0) h2 = TString::Format( "HTle%d",(int)maxHT2);
diff_oname += h1 + "vs" + h2 + "_";
}
if(HTselection==1)       tmp_oname += "all";
if(HTselection==2)       tmp_oname += "highVPt";
if(HTselection==3)       tmp_oname += "highHT";

diff_oname += tmp_oname;
    

string filename1 = fname1.Data();
string filename2 = fname2.Data();


//load the files and histograms
string outname   = diff_oname.Data();
string fn1 = "../RootMacros/GammaVsZllStudies/"+filename1 + ".root";
string fn2 = "../RootMacros/GammaVsZllStudies/"+filename2 + ".root";
TFile *f1 = TFile::Open(fn1.c_str()); 
TFile *f2 = TFile::Open(fn2.c_str()); 
cout << "loaded files " << fn1 << " and " << fn2 << endl;

string hname;
f1->cd();
hname = filename1+"_Z";     TH1D *hZ_1     = (TH1D*)f1->Get(hname.c_str());
hname = filename1+"_Zbg";   TH1D *hZbg_1   = (TH1D*)f1->Get(hname.c_str());
hname = filename1+"_Zdata"; TH1D *hZdata_1 = (TH1D*)f1->Get(hname.c_str());
hname = filename1+"_G";     TH1D *hG_1     = (TH1D*)f1->Get(hname.c_str());
hname = filename1+"_Gbg";   TH1D *hGbg_1   = (TH1D*)f1->Get(hname.c_str());
hname = filename1+"_Gdata"; TH1D *hGdata_1 = (TH1D*)f1->Get(hname.c_str());
f2->cd();
hname = filename2+"_Z";     TH1D *hZ_2     = (TH1D*)f2->Get(hname.c_str());
hname = filename2+"_Zbg";   TH1D *hZbg_2   = (TH1D*)f2->Get(hname.c_str());
hname = filename2+"_Zdata"; TH1D *hZdata_2 = (TH1D*)f2->Get(hname.c_str());
hname = filename2+"_G";     TH1D *hG_2     = (TH1D*)f2->Get(hname.c_str());
hname = filename2+"_Gbg";   TH1D *hGbg_2   = (TH1D*)f2->Get(hname.c_str());
hname = filename2+"_Gdata"; TH1D *hGdata_2 = (TH1D*)f2->Get(hname.c_str());

if(rescaleZ>0){
cout << "rescale Z MC by " << rescaleZ << endl;
hZ_1->Scale(rescaleZ);
hZ_2->Scale(rescaleZ);
}
if(rescaleZbg>0){
cout << "rescale Zbg MC by " << rescaleZbg << endl;
hZbg_1->Scale(rescaleZbg);
hZbg_2->Scale(rescaleZbg);
}
if(rescaleG>0){
cout << "rescale Gamma MC by " << rescaleG << endl;
hG_1->Scale(rescaleG);
hG_2->Scale(rescaleG);
}
if(rescaleGbg>0){
cout << "rescale GammaBg MC by " << rescaleGbg << endl;
hGbg_1->Scale(rescaleGbg);
hGbg_2->Scale(rescaleGbg);
}
if(removebkg){
hZdata_1->Add(hZbg_1,-1);
hGdata_1->Add(hGbg_1,-1);
hZdata_2->Add(hZbg_2,-1);
hGdata_2->Add(hGbg_2,-1);
}

if(rebin>1){
hZ_1->Rebin(rebin);
hZ_2->Rebin(rebin);
hZbg_1->Rebin(rebin);
hZbg_2->Rebin(rebin);
hG_1->Rebin(rebin);
hG_2->Rebin(rebin);
hGbg_1->Rebin(rebin);
hGbg_2->Rebin(rebin);
hZdata_1->Rebin(rebin);
hGdata_1->Rebin(rebin);
hZdata_2->Rebin(rebin);
hGdata_2->Rebin(rebin);
}

if(removebkg) outname   += "__bgRemoved";
if(removebkg) filename1 += "__bgRemoved";
if(removebkg) filename2 += "__bgRemoved";

//define possible y-axis label
TString yTitle;
int binwidth = hZdata_1->GetBinWidth(1);
if(varname=="VPt"||varname=="MT2"||varname=="HT") yTitle = TString::Format(" / %d GeV", binwidth);
else                                              yTitle = "";//NJets, NBJets;


TLatex LumiBox;
LumiBox.SetNDC();
LumiBox.SetTextSize(0.04181185);
TString lumi = TString::Format("%1.2f",19.5);//define luminosity

//build all the ratios
//ratio Z_1/G_1
TH1D *MCratioZG_1 = (TH1D*)hZ_1->Clone("MCratioZG_1");
MCratioZG_1->Divide(hG_1);
//ratio Zdata_1/Gdata_1
TH1D *DataratioZG_1 = (TH1D*)hZdata_1->Clone("DataratioZG_1");
DataratioZG_1->Divide(hGdata_1);
//ratio Zdata_1/Z_1, Gdata_1/G_1
TH1D *DataMCratioZ_1 = (TH1D*)hZdata_1->Clone("DataMCratioZ_1");
DataMCratioZ_1->Divide(hZ_1);
TH1D *DataMCratioG_1 = (TH1D*)hGdata_1->Clone("DataMCratioG_1");
DataMCratioG_1->Divide(hG_1);
//double ratio Zdata_1/Z_1 / Gdata_1/G_1
TH1D *DoubleRatioDataMCZG_1 = (TH1D*)DataMCratioZ_1->Clone("DoubleRatioDataMCZG_1");
DoubleRatioDataMCZG_1->Divide(DataMCratioG_1);
//ratio Z_2/G_2, Zdata_2/Gdata_2; Zdata_2/Z_2, Gdata_2/G_2
TH1D *MCratioZG_2 = (TH1D*)hZ_2->Clone("MCratioZG_2");
MCratioZG_2->Divide(hG_2);
TH1D *DataratioZG_2 = (TH1D*)hZdata_2->Clone("DataratioZG_2");
DataratioZG_2->Divide(hGdata_2);
TH1D *DataMCratioZ_2 = (TH1D*)hZdata_2->Clone("DataMCratioZ_2");
DataMCratioZ_2->Divide(hZ_2);
TH1D *DataMCratioG_2 = (TH1D*)hGdata_2->Clone("DataMCratioG_2");
DataMCratioG_2->Divide(hG_2);
//double ratio Zdata_2/Z_2 / Gdata_2/G_2
TH1D *DoubleRatioDataMCZG_2 = (TH1D*)DataMCratioZ_2->Clone("DoubleRatioDataMCZG_2");
DoubleRatioDataMCZG_2->Divide(DataMCratioG_2);
//ratio Z_1/Z_2, Zdata_1/Zdata_2
TH1D *MCratioZ1Z2 = (TH1D*)hZ_1->Clone("MCratioZ1Z2");
MCratioZ1Z2->Divide(hZ_2);
TH1D *DataratioZ1Z2 = (TH1D*)hZdata_1->Clone("DataratioZ1Z2");
DataratioZ1Z2->Divide(hZdata_2);
//ratio G_1/G_2, Gdata_1/Gdata_2
TH1D *MCratioG1G2 = (TH1D*)hG_1->Clone("MCratioG1G2");
MCratioG1G2->Divide(hG_2);
TH1D *DataratioG1G2 = (TH1D*)hGdata_1->Clone("DataratioG1G2");
DataratioG1G2->Divide(hGdata_2);
//ratio Z_1/Z_2 / G_1/G_2; Zdata_1/Zdata_2 / Gdata_1/Gdata_2
TH1D *MCratioZZGG = (TH1D*)MCratioZ1Z2->Clone("MCratioZZGG");
MCratioZZGG->Divide(MCratioG1G2);
TH1D *DataratioZZGG = (TH1D*)DataratioZ1Z2->Clone("DataratioZZGG");
DataratioZZGG->Divide(DataratioG1G2);

TF1 *fitZG1, *fitZG2, *fitZZ, *fitGG, *fitZZGG;
TF1 *MCfitZG1, *MCfitZG2, *MCfitZZ, *MCfitGG, *MCfitZZGG;
if(dofit){
if(fit300<0) fitZG1  = new TF1("fitZG1","[0]", DataratioZG_1->GetBinLowEdge(1), DataratioZG_1->GetBinLowEdge(DataratioZG_1->GetNbinsX()) + DataratioZG_1->GetBinWidth(DataratioZG_1->GetNbinsX()) );
else         fitZG1  = new TF1("fitZG1","[0]", fit300, DataratioZG_1->GetBinLowEdge(DataratioZG_1->GetNbinsX()) + DataratioZG_1->GetBinWidth(DataratioZG_1->GetNbinsX()) );
if(fit300<0) fitZG2  = new TF1("fitZG2","[0]", DataratioZG_2->GetBinLowEdge(1), DataratioZG_2->GetBinLowEdge(DataratioZG_2->GetNbinsX()) + DataratioZG_2->GetBinWidth(DataratioZG_2->GetNbinsX()) );
else         fitZG2  = new TF1("fitZG2","[0]", fit300, DataratioZG_2->GetBinLowEdge(DataratioZG_2->GetNbinsX()) + DataratioZG_2->GetBinWidth(DataratioZG_2->GetNbinsX()) );
if(!fit2par) fitZZ   = new TF1("fitZZ","[0]", DataratioZ1Z2->GetBinLowEdge(1), DataratioZ1Z2->GetBinLowEdge(DataratioZ1Z2->GetNbinsX()) + DataratioZ1Z2->GetBinWidth(DataratioZ1Z2->GetNbinsX()) );
else         fitZZ   = new TF1("fitZZ","[0]+[1]*x", DataratioZ1Z2->GetBinLowEdge(1), DataratioZ1Z2->GetBinLowEdge(DataratioZ1Z2->GetNbinsX()) + DataratioZ1Z2->GetBinWidth(DataratioZ1Z2->GetNbinsX()) );
if(!fit2par) fitGG   = new TF1("fitGG","[0]", DataratioG1G2->GetBinLowEdge(1), DataratioG1G2->GetBinLowEdge(DataratioG1G2->GetNbinsX()) + DataratioG1G2->GetBinWidth(DataratioG1G2->GetNbinsX()) );
else         fitGG   = new TF1("fitGG","[0]+[1]*x", DataratioG1G2->GetBinLowEdge(1), DataratioG1G2->GetBinLowEdge(DataratioG1G2->GetNbinsX()) + DataratioG1G2->GetBinWidth(DataratioG1G2->GetNbinsX()) );
             fitZZGG = new TF1("fitZZGG","[0]", DataratioZZGG->GetBinLowEdge(1), DataratioZZGG->GetBinLowEdge(DataratioZZGG->GetNbinsX()) + DataratioZZGG->GetBinWidth(DataratioZZGG->GetNbinsX()) );
//copy histograms - remove bins with <=0 content
TH1D *DataratioZG_1_cp = new TH1D("DataratioZG_1_cp", "", DataratioZG_1->GetNbinsX(), DataratioZG_1->GetBinLowEdge(1), DataratioZG_1->GetBinLowEdge(DataratioZG_1->GetNbinsX()) + DataratioZG_1->GetBinWidth(DataratioZG_1->GetNbinsX()));
TH1D *DataratioZG_2_cp = new TH1D("DataratioZG_2_cp", "", DataratioZG_2->GetNbinsX(), DataratioZG_2->GetBinLowEdge(1), DataratioZG_2->GetBinLowEdge(DataratioZG_2->GetNbinsX()) + DataratioZG_2->GetBinWidth(DataratioZG_2->GetNbinsX()));
TH1D *DataratioZ1Z2_cp = new TH1D("DataratioZ1Z2_cp", "", DataratioZ1Z2->GetNbinsX(), DataratioZ1Z2->GetBinLowEdge(1), DataratioZ1Z2->GetBinLowEdge(DataratioZ1Z2->GetNbinsX()) + DataratioZ1Z2->GetBinWidth(DataratioZ1Z2->GetNbinsX()));
TH1D *DataratioG1G2_cp = new TH1D("DataratioG1G2_cp", "", DataratioG1G2->GetNbinsX(), DataratioG1G2->GetBinLowEdge(1), DataratioG1G2->GetBinLowEdge(DataratioG1G2->GetNbinsX()) + DataratioG1G2->GetBinWidth(DataratioG1G2->GetNbinsX()));
TH1D *DataratioZZGG_cp = new TH1D("DataratioZZGG_cp", "", DataratioZZGG->GetNbinsX(), DataratioZZGG->GetBinLowEdge(1), DataratioZZGG->GetBinLowEdge(DataratioZZGG->GetNbinsX()) + DataratioZZGG->GetBinWidth(DataratioZZGG->GetNbinsX()));
for(int ty = 1; ty<=DataratioZZGG_cp->GetNbinsX(); ++ty){
if(DataratioZG_1->GetBinContent(ty)>cutoff) { DataratioZG_1_cp->SetBinContent(ty, DataratioZG_1->GetBinContent(ty)); DataratioZG_1_cp->SetBinError(ty, DataratioZG_1->GetBinError(ty)); }
if(DataratioZG_2->GetBinContent(ty)>cutoff) { DataratioZG_2_cp->SetBinContent(ty, DataratioZG_2->GetBinContent(ty)); DataratioZG_2_cp->SetBinError(ty, DataratioZG_2->GetBinError(ty)); }
if(DataratioZ1Z2->GetBinContent(ty)>cutoff) { DataratioZ1Z2_cp->SetBinContent(ty, DataratioZ1Z2->GetBinContent(ty)); DataratioZ1Z2_cp->SetBinError(ty, DataratioZ1Z2->GetBinError(ty)); }
if(DataratioG1G2->GetBinContent(ty)>cutoff) { DataratioG1G2_cp->SetBinContent(ty, DataratioG1G2->GetBinContent(ty)); DataratioG1G2_cp->SetBinError(ty, DataratioG1G2->GetBinError(ty)); }
if(DataratioZZGG->GetBinContent(ty)>cutoff) { DataratioZZGG_cp->SetBinContent(ty, DataratioZZGG->GetBinContent(ty)); DataratioZZGG_cp->SetBinError(ty, DataratioZZGG->GetBinError(ty)); }
}
if(fitMC){
if(fit300<0) MCfitZG1  = new TF1("MCfitZG1","[0]", MCratioZG_1->GetBinLowEdge(1), MCratioZG_1->GetBinLowEdge(MCratioZG_1->GetNbinsX()) + MCratioZG_1->GetBinWidth(MCratioZG_1->GetNbinsX()) );
else         MCfitZG1  = new TF1("MCfitZG1","[0]", fit300, MCratioZG_1->GetBinLowEdge(MCratioZG_1->GetNbinsX()) + MCratioZG_1->GetBinWidth(MCratioZG_1->GetNbinsX()) );
if(fit300<0) MCfitZG2  = new TF1("MCfitZG2","[0]", MCratioZG_2->GetBinLowEdge(1), MCratioZG_2->GetBinLowEdge(MCratioZG_2->GetNbinsX()) + MCratioZG_2->GetBinWidth(MCratioZG_2->GetNbinsX()) );
else         MCfitZG2  = new TF1("MCfitZG2","[0]", fit300, MCratioZG_2->GetBinLowEdge(MCratioZG_2->GetNbinsX()) + MCratioZG_2->GetBinWidth(MCratioZG_2->GetNbinsX()) );
if(!fit2par) MCfitZZ   = new TF1("MCfitZZ","[0]", MCratioZ1Z2->GetBinLowEdge(1), MCratioZ1Z2->GetBinLowEdge(MCratioZ1Z2->GetNbinsX()) + MCratioZ1Z2->GetBinWidth(MCratioZ1Z2->GetNbinsX()) );
else         MCfitZZ   = new TF1("MCfitZZ","[0]+[1]*x", MCratioZ1Z2->GetBinLowEdge(1), MCratioZ1Z2->GetBinLowEdge(MCratioZ1Z2->GetNbinsX()) + MCratioZ1Z2->GetBinWidth(MCratioZ1Z2->GetNbinsX()) );
if(!fit2par) MCfitGG   = new TF1("MCfitGG","[0]", MCratioG1G2->GetBinLowEdge(1), MCratioG1G2->GetBinLowEdge(MCratioG1G2->GetNbinsX()) + MCratioG1G2->GetBinWidth(MCratioG1G2->GetNbinsX()) );
else         MCfitGG   = new TF1("MCfitGG","[0]+[1]*x", MCratioG1G2->GetBinLowEdge(1), MCratioG1G2->GetBinLowEdge(MCratioG1G2->GetNbinsX()) + MCratioG1G2->GetBinWidth(MCratioG1G2->GetNbinsX()) );
             MCfitZZGG = new TF1("MCfitZZGG","[0]", MCratioZZGG->GetBinLowEdge(1), MCratioZZGG->GetBinLowEdge(MCratioZZGG->GetNbinsX()) + MCratioZZGG->GetBinWidth(MCratioZZGG->GetNbinsX()) );
//copy histograms - remove bins with <=0 content
TH1D *MCratioZG_1_cp = new TH1D("MCratioZG_1_cp", "", MCratioZG_1->GetNbinsX(), MCratioZG_1->GetBinLowEdge(1), MCratioZG_1->GetBinLowEdge(MCratioZG_1->GetNbinsX()) + MCratioZG_1->GetBinWidth(MCratioZG_1->GetNbinsX()));
TH1D *MCratioZG_2_cp = new TH1D("MCratioZG_2_cp", "", MCratioZG_2->GetNbinsX(), MCratioZG_2->GetBinLowEdge(1), MCratioZG_2->GetBinLowEdge(MCratioZG_2->GetNbinsX()) + MCratioZG_2->GetBinWidth(MCratioZG_2->GetNbinsX()));
TH1D *MCratioZ1Z2_cp = new TH1D("MCratioZ1Z2_cp", "", MCratioZ1Z2->GetNbinsX(), MCratioZ1Z2->GetBinLowEdge(1), MCratioZ1Z2->GetBinLowEdge(MCratioZ1Z2->GetNbinsX()) + MCratioZ1Z2->GetBinWidth(MCratioZ1Z2->GetNbinsX()));
TH1D *MCratioG1G2_cp = new TH1D("MCratioG1G2_cp", "", MCratioG1G2->GetNbinsX(), MCratioG1G2->GetBinLowEdge(1), MCratioG1G2->GetBinLowEdge(MCratioG1G2->GetNbinsX()) + MCratioG1G2->GetBinWidth(MCratioG1G2->GetNbinsX()));
TH1D *MCratioZZGG_cp = new TH1D("MCratioZZGG_cp", "", MCratioZZGG->GetNbinsX(), MCratioZZGG->GetBinLowEdge(1), MCratioZZGG->GetBinLowEdge(MCratioZZGG->GetNbinsX()) + MCratioZZGG->GetBinWidth(MCratioZZGG->GetNbinsX()));
for(int ty = 1; ty<=MCratioZZGG_cp->GetNbinsX(); ++ty){
if(MCratioZG_1->GetBinContent(ty)>cutoff) { MCratioZG_1_cp->SetBinContent(ty, MCratioZG_1->GetBinContent(ty)); MCratioZG_1_cp->SetBinError(ty, MCratioZG_1->GetBinError(ty)); }
if(MCratioZG_2->GetBinContent(ty)>cutoff) { MCratioZG_2_cp->SetBinContent(ty, MCratioZG_2->GetBinContent(ty)); MCratioZG_2_cp->SetBinError(ty, MCratioZG_2->GetBinError(ty)); }
if(MCratioZ1Z2->GetBinContent(ty)>cutoff) { MCratioZ1Z2_cp->SetBinContent(ty, MCratioZ1Z2->GetBinContent(ty)); MCratioZ1Z2_cp->SetBinError(ty, MCratioZ1Z2->GetBinError(ty)); }
if(MCratioG1G2->GetBinContent(ty)>cutoff) { MCratioG1G2_cp->SetBinContent(ty, MCratioG1G2->GetBinContent(ty)); MCratioG1G2_cp->SetBinError(ty, MCratioG1G2->GetBinError(ty)); }
if(MCratioZZGG->GetBinContent(ty)>cutoff) { MCratioZZGG_cp->SetBinContent(ty, MCratioZZGG->GetBinContent(ty)); MCratioZZGG_cp->SetBinError(ty, MCratioZZGG->GetBinError(ty)); }
}
}
//plot out all the fit falues
DataratioZG_1_cp->Fit(fitZG1 ,"R");
cout << "fitZG1:  " << "constant = " << fitZG1->GetParameter(0) << " +/- " << fitZG1->GetParError(0);
if(fit2par) cout << ", slope = " << fitZG1->GetParameter(1) << " +/- " << fitZG1->GetParError(1);
cout << ", chi2/N = " << fitZG1->GetChisquare() << "/" << fitZG1->GetNDF() << " = " << fitZG1->GetChisquare()/fitZG1->GetNDF() << endl << endl;
DataratioZG_2_cp->Fit(fitZG2 ,"R");
cout << "fitZG2:  " << "constant = " << fitZG2->GetParameter(0) << " +/- " << fitZG2->GetParError(0);
if(fit2par) cout << ", slope = " << fitZG2->GetParameter(1) << " +/- " << fitZG2->GetParError(1);
cout << ", chi2/N = " << fitZG2->GetChisquare() << "/" << fitZG2->GetNDF() << " = " << fitZG2->GetChisquare()/fitZG2->GetNDF() << endl << endl;
DataratioZ1Z2_cp->Fit(fitZZ  ,"R");
cout << "fitZZ:   " << "constant = " << fitZZ->GetParameter(0) << " +/- " << fitZZ->GetParError(0);
if(fit2par) cout << ", slope = " << fitZZ->GetParameter(1) << " +/- " << fitZZ->GetParError(1);
cout << ", chi2/N = " << fitZZ->GetChisquare() << "/" << fitZZ->GetNDF() << " = " << fitZZ->GetChisquare()/fitZZ->GetNDF() << endl << endl;
DataratioG1G2_cp->Fit(fitGG  ,"R");
cout << "fitGG:   " << "constant = " << fitGG->GetParameter(0) << " +/- " << fitGG->GetParError(0);
if(fit2par) cout << ", slope = " << fitGG->GetParameter(1) << " +/- " << fitGG->GetParError(1);
cout << ", chi2/N = " << fitGG->GetChisquare() << "/" << fitGG->GetNDF() << " = " << fitGG->GetChisquare()/fitGG->GetNDF() << endl << endl;
DataratioZZGG_cp->Fit(fitZZGG,"R");
cout << "fitZZGG: " << "constant = " << fitZZGG->GetParameter(0) << " +/- " << fitZZGG->GetParError(0);
cout << ", chi2/N = " << fitZZGG->GetChisquare() << "/" << fitZZGG->GetNDF() << " = " << fitZZGG->GetChisquare()/fitZZGG->GetNDF() << endl << endl;
if(fitMC){
MCratioZG_1_cp->Fit(MCfitZG1 ,"RQ0");
cout << "MCfitZG1:  " << "constant = " << MCfitZG1->GetParameter(0) << " +/- " << MCfitZG1->GetParError(0);
if(fit2par) cout << ", slope = " << MCfitZG1->GetParameter(1) << " +/- " << MCfitZG1->GetParError(1);
cout << ", chi2/N = " << MCfitZG1->GetChisquare() << "/" << MCfitZG1->GetNDF() << " = " << MCfitZG1->GetChisquare()/MCfitZG1->GetNDF() << endl << endl;
MCratioZG_2_cp->Fit(MCfitZG2 ,"RQ0");
cout << "MCfitZG2:  " << "constant = " << MCfitZG2->GetParameter(0) << " +/- " << MCfitZG2->GetParError(0);
if(fit2par) cout << ", slope = " << MCfitZG2->GetParameter(1) << " +/- " << MCfitZG2->GetParError(1);
cout << ", chi2/N = " << MCfitZG2->GetChisquare() << "/" << MCfitZG2->GetNDF() << " = " << MCfitZG2->GetChisquare()/MCfitZG2->GetNDF() << endl << endl;
MCratioZ1Z2_cp->Fit(MCfitZZ  ,"RQ0");
cout << "MCfitZZ:   " << "constant = " << MCfitZZ->GetParameter(0) << " +/- " << MCfitZZ->GetParError(0);
if(fit2par) cout << ", slope = " << MCfitZZ->GetParameter(1) << " +/- " << MCfitZZ->GetParError(1);
cout << ", chi2/N = " << MCfitZZ->GetChisquare() << "/" << MCfitZZ->GetNDF() << " = " << MCfitZZ->GetChisquare()/MCfitZZ->GetNDF() << endl << endl;
MCratioG1G2_cp->Fit(MCfitGG  ,"RQ0");
cout << "MCfitGG:   " << "constant = " << MCfitGG->GetParameter(0) << " +/- " << MCfitGG->GetParError(0);
if(fit2par) cout << ", slope = " << MCfitGG->GetParameter(1) << " +/- " << MCfitGG->GetParError(1);
cout << ", chi2/N = " << MCfitGG->GetChisquare() << "/" << MCfitGG->GetNDF() << " = " << MCfitGG->GetChisquare()/MCfitGG->GetNDF() << endl << endl;
MCratioZZGG_cp->Fit(MCfitZZGG,"RQ0");
cout << "MCfitZZGG: " << "constant = " << MCfitZZGG->GetParameter(0) << " +/- " << MCfitZZGG->GetParError(0);
cout << ", chi2/N = " << MCfitZZGG->GetChisquare() << "/" << MCfitZZGG->GetNDF() << " = " << MCfitZZGG->GetChisquare()/MCfitZZGG->GetNDF() << endl << endl;
}
}

if(dofit) outname   = "Fit_"+outname;
if(dofit) filename1 = "Fit_"+filename1;
if(dofit) filename2 = "Fit_"+filename2;

//save all histograms, ratios, and fits
string on = "GammaVsZllStudies/RatioHistos/"+outname + ".root";
TFile *nf = new TFile(on.c_str(), "RECREATE"); 
nf->cd();
//save original histos
hZ_1                 ->Write();
hZbg_1               ->Write();
hZdata_1             ->Write();
hG_1                 ->Write();
hGbg_1               ->Write();
hGdata_1             ->Write();
hZ_2                 ->Write();
hZbg_2               ->Write();
hZdata_2             ->Write();
hG_2                 ->Write();
hGbg_2               ->Write();
hGdata_2             ->Write();
//save ratios
MCratioZG_1          ->Write();
DataratioZG_1        ->Write();
DataMCratioZ_1       ->Write();
DataMCratioG_1       ->Write();
MCratioZG_2          ->Write();
DataratioZG_2        ->Write();
DataMCratioZ_2       ->Write();
DataMCratioG_2       ->Write();
MCratioZ1Z2          ->Write();
DataratioZ1Z2        ->Write();
MCratioG1G2          ->Write();
DataratioG1G2        ->Write();
MCratioZZGG          ->Write();
DataratioZZGG        ->Write();
DoubleRatioDataMCZG_1->Write();
DoubleRatioDataMCZG_2->Write();
if(dofit){
fitZG1               ->Write();
fitZG2               ->Write();
fitZZ                ->Write();
fitGG                ->Write();
fitZZGG              ->Write();
}
cout << "saved histos in " << nf->GetName() << endl;

if(rebin==2){
outname   += "_doublebinsize";
filename1 += "_doublebinsize";
filename2 += "_doublebinsize";
} else if(rebin>2){
outname   += TString::Format("_%dtimesbinsize"  ,  rebin);
filename1 += TString::Format("_%dtimesbinsize"  ,  rebin);
filename2 += TString::Format("_%dtimesbinsize"  ,  rebin);
}

//now all the print out everything you want
if(printDistZvsZorGvsG || printDistZvsG || printRatiosDataVsMC || printRatiosZVsG || printRatiosZvsZdivGvsG || printRatiosZvsZorGvsG){

TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",334,93,600,600);
gStyle->SetOptFit(1);           gStyle->SetOptStat(0);           gStyle->SetOptTitle(0);
Canvas_1->Range(-82.79222,-0.9532603,793.8312,3.394941);         Canvas_1->SetTickx(1);          Canvas_1->SetTicky(1);
Canvas_1->SetFillColor(0);      Canvas_1->SetBorderMode(0);      Canvas_1->SetBorderSize(2);
Canvas_1->SetRightMargin(0.05); Canvas_1->SetTopMargin(0.07);    Canvas_1->SetBottomMargin(0.15);Canvas_1->SetLeftMargin(0.18);
Canvas_1->SetFrameFillStyle(0); Canvas_1->SetFrameBorderMode(0); Canvas_1->SetFrameBorderMode(0);

TLegend *leg4_dist = new TLegend(0.6,0.68,0.8,0.88,NULL,"brNDC");
leg4_dist->SetBorderSize(0); leg4_dist->SetTextFont(42); leg4_dist->SetTextSize(0.04181185); leg4_dist->SetLineColor(1);
leg4_dist->SetLineStyle(1);  leg4_dist->SetLineWidth(2); leg4_dist->SetFillColor(0);         leg4_dist->SetFillStyle(1001);

TLegend *leg2_dist = new TLegend(0.6,0.68,0.8,0.88,NULL,"brNDC");
leg2_dist->SetBorderSize(0); leg2_dist->SetTextFont(42); leg2_dist->SetTextSize(0.04181185); leg2_dist->SetLineColor(1);
leg2_dist->SetLineStyle(1);  leg2_dist->SetLineWidth(2); leg2_dist->SetFillColor(0);         leg2_dist->SetFillStyle(1001);

TLegend *leg_rat   = new TLegend(0.25,0.78,0.45,0.88,NULL,"brNDC");
leg_rat  ->SetBorderSize(0); leg_rat  ->SetTextFont(42); leg_rat  ->SetTextSize(0.04181185); leg_rat  ->SetLineColor(1);
leg_rat  ->SetLineStyle(1);  leg_rat  ->SetLineWidth(2); leg_rat  ->SetFillColor(0);         leg_rat  ->SetFillStyle(1001);

//copy
TString haxistitle = "Events" + yTitle;
TH1D *haxis = (TH1D*)hGdata_2->Clone("haxis");
haxis->GetYaxis()->SetTitle(haxistitle.Data());
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
double lowborder = haxis->GetBinLowEdge(1);
double upborder  = haxis->GetBinLowEdge(haxis->GetNbinsX()) + haxis->GetBinWidth(haxis->GetNbinsX());
if(lowbin>0 && lowbin>lowborder) lowborder = lowbin;
if(upbin >0 && upbin <upborder ) upborder  = upbin;
haxis->GetXaxis()->SetRangeUser(lowborder, upborder);
//haxis->Draw("axis");

string savename = "GammaVsZllStudies/Plots/";


if(printDistZvsZorGvsG){
Canvas_1->Clear();
TH1D *Gd2 = (TH1D*)hGdata_2->Clone("Gd2"); Gd2->SetMarkerStyle(26);
TH1D *Gd1 = (TH1D*)hGdata_1->Clone("Gd1"); Gd1->SetMarkerStyle(32);
TH1D *Gm2 = (TH1D*)hG_2    ->Clone("Gm2"); Gm2->SetMarkerStyle(22);
TH1D *Gm1 = (TH1D*)hG_1    ->Clone("Gm1"); Gm1->SetMarkerStyle(23);
leg4_dist->AddEntry(Gd2, "#gamma(0b) data",       "lp");
leg4_dist->AddEntry(Gd1, "#gamma(1b) data",       "lp");
leg4_dist->AddEntry(Gm2, "#gamma(0b) simulation", "lp");
leg4_dist->AddEntry(Gm1, "#gamma(1b) simulation", "lp");
leg2_dist->AddEntry(Gd2, "#gamma(0b) data",       "lp");
leg2_dist->AddEntry(Gd1, "#gamma(1b) data",       "lp");
double maxd2 = Gd2->GetMaximum(); double maxd1 = Gd1->GetMaximum(); double maxm2 = Gm2->GetMaximum(); double maxm1 = Gm1->GetMaximum();
double max = (maxd1>maxd2)?maxd1:maxd2; max = (max>maxm2)?max:maxm2; max = (max>maxm1)?max:maxm1;
max = 3.*max;

double mind2 = Gd2->GetMinimum(1E-8); double mind1 = Gd1->GetMinimum(1E-8); double minm2 = Gm2->GetMinimum(1E-8); double minm1 = Gm1->GetMinimum(1E-8);
double min = (mind1<mind2)?mind1:mind2; min = (min<minm2)?min:minm2; min = (min<minm1)?min:minm1;
min = min/2.5;

Canvas_1->cd();
haxis->GetYaxis()->SetRangeUser(min,max);
gPad->SetLogy(1);
if(dataonly)(
haxis->Draw("axis");
Gd1->Draw("sameE");
Gd2->Draw("sameE");
leg2_dist->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/DistZvsZorGvsG/GvsG/Data/"+outname+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/DistZvsZorGvsG/GvsG/Data/"+outname+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/DistZvsZorGvsG/GvsG/Data/"+outname+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c1 = (TCanvas*)Canvas_1->DrawClone();
Canvas_1->Clear();
)
if(dataandMC){
Canvas_1->cd();
haxis->Draw("axis");
Gm1->Draw("sameE");
Gm2->Draw("sameE");
Gd1->Draw("sameE");
Gd2->Draw("sameE");
leg4_dist->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/DistZvsZorGvsG/GvsG/DataMC/"+outname+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/DistZvsZorGvsG/GvsG/DataMC/"+outname+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/DistZvsZorGvsG/GvsG/DataMC/"+outname+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c2 = (TCanvas*)Canvas_1->DrawClone();
Canvas_1->Clear();
}

leg4_dist->Clear();
leg2_dist->Clear();

Canvas_1->cd();
TH1D *Zd2 = (TH1D*)hZdata_2->Clone("Zd2"); Zd2->SetMarkerStyle(26);
TH1D *Zd1 = (TH1D*)hZdata_1->Clone("Zd1"); Zd1->SetMarkerStyle(32); 
TH1D *Zm2 = (TH1D*)hZ_2    ->Clone("Zm2"); Zm2->SetMarkerStyle(22);
TH1D *Zm1 = (TH1D*)hZ_1    ->Clone("Zm1"); Zm1->SetMarkerStyle(23); Zm1->SetMarkerColor(kOrange+7), Zm1->SetLineColor(kOrange+7);
leg4_dist->AddEntry(Zd2, "Z_{ll}(0b) data",       "lp");
leg4_dist->AddEntry(Zd1, "Z_{ll}(1b) data",       "lp");
leg4_dist->AddEntry(Zm2, "Z_{ll}(0b) simulation", "lp");
leg4_dist->AddEntry(Zm1, "Z_{ll}(1b) simulation", "lp");
leg2_dist->AddEntry(Zd2, "Z_{ll}(0b) data",       "lp");
leg2_dist->AddEntry(Zd1, "Z_{ll}(1b) data",       "lp");
double maxd2 = Zd2->GetMaximum(); double maxd1 = Zd1->GetMaximum(); double maxm2 = Zm2->GetMaximum(); double maxm1 = Zm1->GetMaximum();
double max = (maxd1>maxd2)?maxd1:maxd2; max = (max>maxm2)?max:maxm2; max = (max>maxm1)?max:maxm1;
max = 3.*max;

double mind2 = Zd2->GetMinimum(1E-8); double mind1 = Zd1->GetMinimum(1E-8); double minm2 = Zm2->GetMinimum(1E-8); double minm1 = Zm1->GetMinimum(1E-8);
double min = (mind1<mind2)?mind1:mind2; min = (min<minm2)?min:minm2; min = (min<minm1)?min:minm1;
min = min/2.5;

haxis->GetYaxis()->SetRangeUser(min,max);
gPad->SetLogy(1);
if(dataonly){
haxis->Draw("axis");
Zd1->Draw("sameE");
Zd2->Draw("sameE");
leg2_dist->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/DistZvsZorGvsG/ZvsZ/Data/"+outname+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/DistZvsZorGvsG/ZvsZ/Data/"+outname+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/DistZvsZorGvsG/ZvsZ/Data/"+outname+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c3 = (TCanvas*)Canvas_1->DrawClone();
Canvas_1->Clear();
}
if(dataandMC){
Canvas_1->cd();
haxis->Draw("axis");
Zm1->Draw("sameE");
Zm2->Draw("sameE");
Zd1->Draw("sameE");
Zd2->Draw("sameE");
leg4_dist->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/DistZvsZorGvsG/ZvsZ/DataMC/"+outname+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/DistZvsZorGvsG/ZvsZ/DataMC/"+outname+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/DistZvsZorGvsG/ZvsZ/DataMC/"+outname+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c4 = (TCanvas*)Canvas_1->DrawClone();
}
leg4_dist->Clear();
leg2_dist->Clear();

}


if(printDistZvsG){
Canvas_1->Clear();
TH1D *Gd2 = (TH1D*)hGdata_2->Clone("Gd2"); Gd2->SetMarkerStyle(24);
TH1D *Gd1 = (TH1D*)hGdata_1->Clone("Gd1"); Gd1->SetMarkerStyle(24);
TH1D *Gm2 = (TH1D*)hG_2    ->Clone("Gm2"); Gm2->SetMarkerStyle(20);
TH1D *Gm1 = (TH1D*)hG_1    ->Clone("Gm1"); Gm1->SetMarkerStyle(20);
TH1D *Zd2 = (TH1D*)hZdata_2->Clone("Zd2"); Zd2->SetMarkerStyle(32);
TH1D *Zd1 = (TH1D*)hZdata_1->Clone("Zd1"); Zd1->SetMarkerStyle(32);
TH1D *Zm2 = (TH1D*)hZ_2    ->Clone("Zm2"); Zm2->SetMarkerStyle(23);
TH1D *Zm1 = (TH1D*)hZ_1    ->Clone("Zm1"); Zm1->SetMarkerStyle(23);
leg2_dist->AddEntry(Gd2, "#gamma data",       "lp");
leg2_dist->AddEntry(Zd2, "Z_{ll} data",       "lp");
leg4_dist->AddEntry(Gd2, "#gamma data",       "lp");
leg4_dist->AddEntry(Zd2, "Z_{ll} data",       "lp");
leg4_dist->AddEntry(Gm2, "#gamma simulation", "lp");
leg4_dist->AddEntry(Zm2, "Z_{ll} simulation", "lp");
double maxd2 = Gd2->GetMaximum(); double maxd1 = Zd2->GetMaximum(); double maxm2 = Gm2->GetMaximum(); double maxm1 = Zm2->GetMaximum();
double max = (maxd1>maxd2)?maxd1:maxd2; max = (max>maxm2)?max:maxm2; max = (max>maxm1)?max:maxm1;
max = 3.*max;
double mind2 = Gd2->GetMinimum(1E-8); double mind1 = Zd2->GetMinimum(1E-8); double minm2 = Gm2->GetMinimum(1E-8); double minm1 = Zm2->GetMinimum(1E-8);
double min = (mind1<mind2)?mind1:mind2; min = (min<minm2)?min:minm2; min = (min<minm1)?min:minm1;
min = min/2.5;

Canvas_1->cd();
haxis->GetYaxis()->SetRangeUser(min,max);
gPad->SetLogy(1);
if(dataonly){
haxis->Draw("axis");
Gd2->Draw("sameE");
Zd2->Draw("sameE");
leg2_dist->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/DistZvsG/Data/"+filename2+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/DistZvsG/Data/"+filename2+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/DistZvsG/Data/"+filename2+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c7 = (TCanvas*)Canvas_1->DrawClone();
Canvas_1->Clear();
}
if(dataandMC){
Canvas_1->cd();
haxis->Draw("axis");
Gm2->Draw("sameE");
Zm2->Draw("sameE");
Gd2->Draw("sameE");
Zd2->Draw("sameE");
leg4_dist->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/DistZvsG/DataMC/"+filename2+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/DistZvsG/DataMC/"+filename2+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/DistZvsG/DataMC/"+filename2+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c8 = (TCanvas*)Canvas_1->DrawClone();
Canvas_1->Clear();
}

leg2_dist->Clear();
leg4_dist->Clear();

leg2_dist->AddEntry(Gd1, "#gamma data",       "lp");
leg2_dist->AddEntry(Zd1, "Z_{ll} data",       "lp");
leg4_dist->AddEntry(Gd1, "#gamma data",       "lp");
leg4_dist->AddEntry(Zd1, "Z_{ll} data",       "lp");
leg4_dist->AddEntry(Gm1, "#gamma simulation", "lp");
leg4_dist->AddEntry(Zm1, "Z_{ll} simulation", "lp");
double maxd2 = Gd1->GetMaximum(); double maxd1 = Zd1->GetMaximum(); double maxm2 = Gm1->GetMaximum(); double maxm1 = Zm1->GetMaximum();
double max = (maxd1>maxd2)?maxd1:maxd2; max = (max>maxm2)?max:maxm2; max = (max>maxm1)?max:maxm1;
max = 3.*max;
double mind2 = Gd1->GetMinimum(1E-8); double mind1 = Zd1->GetMinimum(1E-8); double minm2 = Gm1->GetMinimum(1E-8); double minm1 = Zm1->GetMinimum(1E-8);
double min = (mind1<mind2)?mind1:mind2; min = (min<minm2)?min:minm2; min = (min<minm1)?min:minm1;
min = min/2.5;

Canvas_1->cd();
haxis->GetYaxis()->SetRangeUser(min,max);
gPad->SetLogy(1);
if(dataonly){
haxis->Draw("axis");
Gd1->Draw("sameE");
Zd1->Draw("sameE");
leg2_dist->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/DistZvsG/Data/"+filename1+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/DistZvsG/Data/"+filename1+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/DistZvsG/Data/"+filename1+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c5 = (TCanvas*)Canvas_1->DrawClone();
Canvas_1->Clear();
}
if(dataandMC){
Canvas_1->cd();
haxis->Draw("axis");
Gm1->Draw("sameE");
Zm1->Draw("sameE");
Gd1->Draw("sameE");
Zd1->Draw("sameE");
leg4_dist->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/DistZvsG/DataMC/"+filename1+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/DistZvsG/DataMC/"+filename1+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/DistZvsG/DataMC/"+filename1+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c6 = (TCanvas*)Canvas_1->DrawClone();
}

leg4_dist->Clear();
leg2_dist->Clear();

}
if(printRatiosDataVsMC){
Canvas_1->Clear();

TH1D *RZ1 = (TH1D*)DataMCratioZ_1->Clone("RZ1"); RZ1->SetMarkerStyle(22); RZ1->SetMarkerColor(kOrange  ); RZ1->SetLineColor(kOrange  );
TH1D *RZ2 = (TH1D*)DataMCratioZ_2->Clone("RZ2"); RZ2->SetMarkerStyle(22); RZ2->SetMarkerColor(kOrange  ); RZ2->SetLineColor(kOrange  );
TH1D *RG1 = (TH1D*)DataMCratioG_1->Clone("RG1"); RG1->SetMarkerStyle(20); RG1->SetMarkerColor(kViolet-3), RG1->SetLineColor(kViolet-3);
TH1D *RG2 = (TH1D*)DataMCratioG_2->Clone("RG2"); RG2->SetMarkerStyle(20); RG2->SetMarkerColor(kViolet-3), RG2->SetLineColor(kViolet-3);
leg_rat->AddEntry(RG1, "#gamma data/simulation", "lp");
leg_rat->AddEntry(RZ2, "Z_{ll} data/simulation", "lp");
double maxd2 = RZ1->GetMaximum(); double maxd1 = RG1->GetMaximum(); 
double max = (maxd1>maxd2)?maxd1:maxd2; 
max = 1.5*max;

Canvas_1->cd();
TString ytitle = "ratio" + yTitle;
haxis->GetYaxis()->SetTitle(ytitle.Data());
haxis->GetYaxis()->SetRangeUser(0.,max);
gPad->SetLogy(0);
haxis->Draw("axis");
RZ1->Draw("sameE");
RG1->Draw("sameE");
leg_rat->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/RatiosDataVsMC/"+filename1+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/RatiosDataVsMC/"+filename1+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/RatiosDataVsMC/"+filename1+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c9 = (TCanvas*)Canvas_1->DrawClone();

Canvas_1->Clear();

double maxd2 = RZ2->GetMaximum(); double maxd1 = RG2->GetMaximum(); 
double max = (maxd1>maxd2)?maxd1:maxd2; 
max = 1.5*max;

Canvas_1->cd();
haxis->GetYaxis()->SetRangeUser(0.,max);
haxis->Draw("axis");
RZ2->Draw("sameE");
RG2->Draw("sameE");
leg_rat->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/RatiosDataVsMC/"+filename2+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/RatiosDataVsMC/"+filename2+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/RatiosDataVsMC/"+filename2+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c10 = (TCanvas*)Canvas_1->DrawClone();

leg_rat->Clear();
}
if(printRatiosZVsG){
Canvas_1->Clear();
TH1D *ZGm1 = (TH1D*)MCratioZG_1  ->Clone("ZGm1"); ZGm1->SetMarkerStyle(20); ZGm1->SetMarkerColor(kViolet-3); ZGm1->SetLineColor(kViolet-3);
TH1D *ZGm2 = (TH1D*)MCratioZG_2  ->Clone("ZGm2"); ZGm2->SetMarkerStyle(20); ZGm2->SetMarkerColor(kViolet-3); ZGm2->SetLineColor(kViolet-3);
TH1D *ZGd1 = (TH1D*)DataratioZG_1->Clone("ZGd1"); ZGd1->SetMarkerStyle(26); ZGd1->SetMarkerColor(kBlack   ), ZGd1->SetLineColor(kBlack   );
TH1D *ZGd2 = (TH1D*)DataratioZG_2->Clone("ZGd2"); ZGd2->SetMarkerStyle(26); ZGd2->SetMarkerColor(kBlack   ), ZGd2->SetLineColor(kBlack   );
leg_rat->AddEntry(ZGd1, "Z_{ll}/#gamma data",       "lp");
leg_rat->AddEntry(ZGm1, "Z_{ll}/#gamma simulation", "lp");
double maxd2 = ZGm1->GetMaximum(); double maxd1 = ZGd1->GetMaximum(); 
double max = (maxd1>maxd2)?maxd1:maxd2; 
max = 1.5*max;

Canvas_1->cd();
TString ytitle = "ratio" + yTitle;
haxis->GetYaxis()->SetTitle(ytitle.Data());
//haxis->GetYaxis()->SetTitle("Z_{ll} / #gamma");
haxis->GetYaxis()->SetRangeUser(0.,max);
gPad->SetLogy(0);
haxis->Draw("axis");
ZGm1->Draw("sameE");
ZGd1->Draw("sameE");
if(dofit) fitZG1->Draw("same");
leg_rat->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/RatiosZVsG/"+filename1+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/RatiosZVsG/"+filename1+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/RatiosZVsG/"+filename1+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c11 = (TCanvas*)Canvas_1->DrawClone();
Canvas_1->Clear();


double maxd2 = ZGm2->GetMaximum(); double maxd1 = ZGd2->GetMaximum(); 
double max = (maxd1>maxd2)?maxd1:maxd2; 
max = 1.5*max;

Canvas_1->cd();
haxis->GetYaxis()->SetRangeUser(0.,max);
haxis->Draw("axis");
ZGm2->Draw("sameE");
ZGd2->Draw("sameE");
if(dofit) fitZG2->Draw("same");
leg_rat->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/RatiosZVsG/"+filename2+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/RatiosZVsG/"+filename2+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/RatiosZVsG/"+filename2+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c12 = (TCanvas*)Canvas_1->DrawClone();

leg_rat->Clear();
}

if(printRatiosZvsZdivGvsG){
Canvas_1->Clear();

TH1D *ZGm = (TH1D*)MCratioZZGG  ->Clone("ZGm"); ZGm->SetMarkerStyle(20); ZGm->SetMarkerColor(kViolet-3); ZGm->SetLineColor(kViolet-3);
TH1D *ZGd = (TH1D*)DataratioZZGG->Clone("ZGd"); ZGd->SetMarkerStyle(26); ZGd->SetMarkerColor(kBlack   ); ZGd->SetLineColor(kBlack   );
leg_rat->AddEntry(ZGd, "R_{Z}/R_{#gamma} (1b/0b) data",       "lp");
leg_rat->AddEntry(ZGm, "R_{Z}/R_{#gamma} (1b/0b) simulation", "lp");
double maxd2 = ZGm->GetMaximum(); double maxd1 = ZGd->GetMaximum(); 
double max = (maxd1>maxd2)?maxd1:maxd2; 
max = 1.5*max;

Canvas_1->cd();
TString ytitle = "double ratio" + yTitle;
haxis->GetYaxis()->SetTitle(ytitle.Data());
haxis->GetYaxis()->SetRangeUser(0.,max);
gPad->SetLogy(0);haxis->Draw("axis");
ZGm->Draw("sameE");
ZGd->Draw("sameE");
if(dofit) fitZZGG->Draw("same");
leg_rat->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/RatiosZvsZdivGvsG/"+outname+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/RatiosZvsZdivGvsG/"+outname+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/RatiosZvsZdivGvsG/"+outname+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c13 = (TCanvas*)Canvas_1->DrawClone();
leg_rat->Clear();

}
if(printRatiosZvsZorGvsG){
Canvas_1->Clear();

TH1D *ZZm = (TH1D*)MCratioZ1Z2  ->Clone("ZZm"); ZZm->SetMarkerStyle(20); ZZm->SetMarkerColor(kOrange  ); ZZm->SetLineColor(kOrange  );
TH1D *ZZd = (TH1D*)DataratioZ1Z2->Clone("ZZd"); ZZd->SetMarkerStyle(26); ZZd->SetMarkerColor(kBlack   ); ZZd->SetLineColor(kBlack   );
TH1D *GGm = (TH1D*)MCratioG1G2  ->Clone("GGm"); GGm->SetMarkerStyle(20); GGm->SetMarkerColor(kViolet-3), GGm->SetLineColor(kViolet-3);
TH1D *GGd = (TH1D*)DataratioG1G2->Clone("GGd"); GGd->SetMarkerStyle(26); GGd->SetMarkerColor(kBlack   ), GGd->SetLineColor(kBlack   );
leg_rat->AddEntry(GGd, "#gamma(1b)/#gamma(0b) data",       "lp");
leg_rat->AddEntry(GGm, "#gamma(1b)/#gamma(0b) simulation", "lp");
double maxd2 = GGd->GetMaximum(); double maxd1 = GGm->GetMaximum(); 
double max = (maxd1>maxd2)?maxd1:maxd2; 
max = 1.5*max;

Canvas_1->cd();
TString ytitle = "ratio" + yTitle;
haxis->GetYaxis()->SetTitle(ytitle.Data());
haxis->GetYaxis()->SetRangeUser(0.,max);
gPad->SetLogy(0);
haxis->Draw("axis");
GGm->Draw("sameE");
GGd->Draw("sameE");
if(dofit) fitGG->Draw("same");
leg_rat->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/RatiosZvsZorGvsG/GvsG/"+outname+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/RatiosZvsZorGvsG/GvsG/"+outname+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/RatiosZvsZorGvsG/GvsG/"+outname+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c14 = (TCanvas*)Canvas_1->DrawClone();
Canvas_1->Clear();

leg_rat->Clear();

leg_rat->AddEntry(ZZd, "Z_{ll}(1b)/Z_{ll}(0b) data",       "lp");
leg_rat->AddEntry(ZZm, "Z_{ll}(1b)/Z_{ll}(0b) simulation", "lp");
double maxd2 = ZZd->GetMaximum(); double maxd1 = ZZm->GetMaximum(); 
double max = (maxd1>maxd2)?maxd1:maxd2; 
max = 1.5*max;
    

Canvas_1->cd();
TString ytitle = "ratio" + yTitle;
haxis->GetYaxis()->SetTitle(ytitle.Data());
haxis->GetYaxis()->SetRangeUser(0.,max);
gPad->SetLogy(0);
haxis->Draw("axis");
ZZm->Draw("sameE");
ZZd->Draw("sameE");
if(dofit) fitZZ->Draw("same");
leg_rat->Draw();
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/RatiosZvsZorGvsG/ZvsZ/"+outname+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/RatiosZvsZorGvsG/ZvsZ/"+outname+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/RatiosZvsZorGvsG/ZvsZ/"+outname+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c15 = (TCanvas*)Canvas_1->DrawClone();


}

Canvas_1->Clear();

TH1D *DR2 = (TH1D*)DoubleRatioDataMCZG_2  ->Clone("DR2"); DR2->SetMarkerStyle(26); DR2->SetMarkerColor(kBlack   ); DR2->SetLineColor(kBlack  );
TH1D *DR1 = (TH1D*)DoubleRatioDataMCZG_1  ->Clone("DR1"); DR1->SetMarkerStyle(26); DR1->SetMarkerColor(kBlack   ); DR1->SetLineColor(kBlack  );
double maxd2 = DR2->GetMaximum(); double maxd1 = DR1->GetMaximum(); 

Canvas_1->cd();
TString ytitle = "#frac{Z_{ll} data/simulation}{#gamma data/simulation} " + yTitle;
haxis->GetYaxis()->SetTitle(ytitle.Data());
haxis->GetYaxis()->SetRangeUser(0.,maxd2);
gPad->SetLogy(0);
haxis->Draw("axis");
DR2->Draw("sameE");
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/DoubleRatioDataVsMCZvsG/"+filename2+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/DoubleRatioDataVsMCZvsG/"+filename2+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/DoubleRatioDataVsMCZvsG/"+filename2+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c16 = (TCanvas*)Canvas_1->DrawClone();
Canvas_1->Clear();

Canvas_1->cd();
haxis->GetYaxis()->SetRangeUser(0.,maxd1);
gPad->SetLogy(0);
haxis->Draw("axis");
DR1->Draw("sameE");
LumiBox.DrawLatex(0.216443,0.95,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");
if(saveC  ){ savename = "GammaVsZllStudies/Plots/DoubleRatioDataVsMCZvsG/"+filename1+".C";   Canvas_1->SaveAs(savename.c_str()); }
if(saveeps){ savename = "GammaVsZllStudies/Plots/DoubleRatioDataVsMCZvsG/"+filename1+".eps"; Canvas_1->SaveAs(savename.c_str()); }
if(savepng){ savename = "GammaVsZllStudies/Plots/DoubleRatioDataVsMCZvsG/"+filename1+".png"; Canvas_1->SaveAs(savename.c_str()); }
if(makeclones) TCanvas *c17 = (TCanvas*)Canvas_1->DrawClone();



}//if(printDistGvsG || printDistZvsG || printRatiosDataVsMC || printRatiosZVsG || printRatiosZvsZdivGvsG || printRatiosZvsZorGvsG)

}
