{
  TString outputdir = "../MassPlots/";
  TString samples   = "./samples/samples.dat";
  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "k");
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  int gNMT2bins                   = 17;
  double  gMT2bins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 400, 550, 800}; 	
  
  int gNMT2Bbins                   = 15;
  double  gMT2Bbins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500}; 	

  int gNMT2bins_l                   = 14;
  double  gMT2bins_l[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500}; 	

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
  tA->SetSave(false);
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
  
	std::ostringstream cutStream;
	cutStream << " " 
          << "misc.MET>=30"                                                << "&&"
//	  << "misc.HT<950&&misc.HT>750  "                                  << "&&"
//	  << "misc.HT>950"                                                 << "&&"
	  << "misc.HT>750"                                                 << "&&"
//	  << "misc.MT2<200&&misc.MT2>150  "                                << "&&"
//	  << "misc.MT2>200  "                                              << "&&"
	  << "NJetsIDLoose40  >=3"                                         << "&&"
//	  << "NBJets ==0"                                                  << "&&"
	  << "(NEles==0  || ele[0].lv.Pt()<10)"                            << "&&"
	  << "(NMuons==0 || muo[0].lv.Pt()<10)"                            << "&&"
//	  << "(muo[0].lv.Pt()>10 || ele[0].lv.Pt()>10)"                    << "&&"
	  << "misc.Jet0Pass ==1"                                           << "&&"
	  << "misc.Jet1Pass ==1"                                           << "&&"
	  << "misc.SecondJPt >100"                                         << "&&"
//	  << "misc.LeadingJPt >150"                                        << "&&"
	  << "misc.PassJetID ==1"                                          << "&&"
	  << "misc.Vectorsumpt < 70"                                       << "&&"
//	  << "(misc.MinMetJetDPhi >0.3||misc.MinMetJetDPhiIndex>3)"        << "&&"
	  << "(misc.MinMetJetDPhi >0.3)"                                   << "&&"
          << "misc.HBHENoiseFlagIso ==0"                                   << "&&"
          << "misc.CSCTightHaloID==0"                                      << "&&"
          << "misc.CrazyHCAL==0";

	std::ostringstream triggerStream;
	triggerStream << "( "
		<< "(trigger.HLT_HT440_v2 ==1 && misc.Run<161216)" << "||"
		<< "(trigger.HLT_HT450_v2 ==1 && (misc.Run>=161216 && misc.Run< 163269))" << "||"
		<< "(trigger.HLT_HT500_v3 ==1 && (misc.Run>=163269 && misc.Run<=163869))" << "||"
		<< "(trigger.HLT_HT500_v4 ==1 && (misc.Run>=165088 && misc.Run< 165970))" << "||"
		<< "(trigger.HLT_HT550_v5 ==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" << "||"
		<< "(trigger.HLT_HT550_v6 ==1 && (misc.Run==166346))" << "||"
		<< "(trigger.HLT_HT550_v7 ==1 && (misc.Run>=167078 && misc.Run< 170249))" << "||"
		<< "(trigger.HLT_HT550_v8 ==1 && (misc.Run>=170249 && misc.Run< 173236))" << "||"
		<< "(trigger.HLT_HT600_v1 ==1 && (misc.Run>=173236 && misc.Run< 178420))" << "||"
		<< "(trigger.HLT_HT650_v4 ==1 && (misc.Run>=178420 && misc.Run< 179959))" << ")";
	TString trigger = triggerStream.str().c_str();
	  TString cuts = cutStream.str().c_str();
    
//                     variable,              cuts,  njet, nlep,     HLT,     xtitle             nbins  bins      flip_order,    log  , comp ,  ratio, stack, overlay
//       tA->makePlot("misc.HT",              cuts,    -1,  -10 ,    trigger , "HT"      ,40,   700, 2000,        false,         true ,  true,   true,  true,  true, 1);
         tA->makePlot("misc.MT2",             cuts,    -1,  -10 ,    trigger , "MT2"      ,70,   0,  700,         false,         true ,  true,   true,  true,  true, 1);
}
