void rootlogon() {
  gSystem->SetIncludePath(" -Iinclude/");
  gSystem->SetIncludePath(" -I../TESCO/include/");
  gSystem->Load("libPhysics");
  gSystem->Load("shlib/libDiLeptonAnalysis.so");
  gROOT->ProcessLine(".x SetStyle_PRD.C");
}
