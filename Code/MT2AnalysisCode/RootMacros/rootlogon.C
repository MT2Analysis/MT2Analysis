void rootlogon() {
  //I think this does not need any explanation
  //maybe you need to adjust the directory of MT2Code and TESCO
  gSystem->SetIncludePath(" -I../MT2Code/include/ -I../ASAnalysis/include/");
  gSystem->Load("libPhysics");
  gSystem->Load("libFWCoreFWLite.so");
  gSystem->Load("libRooFit") ;
  gSystem->Load("libRooFitCore") ;
  AutoLibraryLoader::enable();
  gSystem->Load("libDataFormatsFWLite.so");
  gSystem->Load("../MT2Code/shlib/libDiLeptonAnalysis.so");
  gROOT->ProcessLine(".x SetStyle_PRD.C");
}
