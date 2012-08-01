void rootlogon() {
  gSystem->SetIncludePath(" -I./include/ -I../TESCO/include/");
  gSystem->Load("libPhysics");
  gSystem->Load("shlib/libDiLeptonAnalysis.so");
  gROOT->ProcessLine(".x SetStyle_PRD.C");
}
