{
  //run via root -l -b -q run_CutFlow.C

  //this function calls the MassPlotter function PrintCutFlow(...)
  //this is from the 7 TeV analysis

  TString outputdir = "Temp/";
 TString samples   = "samples/samples_lowMT2_4400_1.dat";
  int verbose = 3;

//  gSystem->Load("libPhysics");
  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "k");//load MassPlotter.cc
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");//define MassPlotter
  tA->setVerbose(verbose);
  tA->init(samples);


//LowMT2
  //the event selection
  std::ostringstream cutStream;
  cutStream	//uncommented one due to cutflow

		<< "misc.MET>=30"                                            << "&&"
		<< "misc.HT<950&&misc.HT>750"                                << "&&"
	//	<< "misc.HT>950 "                                            << "&&"
		<< "NJetsIDLoose40 >=4"                                      << "&&"
		<< "NBJets >0"                                               << "&&"
		<< "(NEles==0  || ele[0].lv.Pt()<10)"                        << "&&"
		<< "(NMuons==0 || muo[0].lv.Pt()<10)"                        << "&&"
		<< "misc.Jet0Pass==1"                                        << "&&"
		<< "misc.Jet1Pass==1"                                        << "&&"
		<< "misc.SecondJPt >100"                                     << "&&"
		<< "misc.LeadingJPt >150"                                    << "&&"
		<< "misc.PassJetID ==1"                                      << "&&"
		<< "misc.Vectorsumpt < 70"                                   << "&&"
		<< "(misc.MinMetJetDPhi >0.3||misc.MinMetJetDPhiIndex>3)"    << "&&"
	//	<< "misc.MinMetJetDPhi >0.3"                                 << "&&"
		<< "misc.HBHENoiseFlagIso==0"                                << "&&"
		<< "misc.CSCTightHaloID==0"                                  << "&&"
		<< "misc.CrazyHCAL==0";

  TString cuts = cutStream.str().c_str();

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
	<< "(trigger.HLT_HT650_v4 ==1 && (misc.Run>=178420 && misc.Run< 179959))" << " )";
	TString trigger = triggerStream.str().c_str();

  //tA->PrintCutFlow(njets,nleps,trigger,cuts);
  //tA->PrintCutFlow(njets,nleps,nbtag,b-tag algo,disc. value,trigger,cuts);

  tA->PrintCutFlow(-4,0,trigger, cuts);   // no btag


}