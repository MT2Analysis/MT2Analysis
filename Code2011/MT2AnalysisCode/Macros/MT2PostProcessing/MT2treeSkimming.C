/*****************************************************************
* root macro to skim MT2trees                                    *
* -> sample and path to shlib given as arguments                 *
*                                                                *
* Pascal Nef                             September 18th, 2011    *  
*****************************************************************/
void MT2treeSkimming(string sample, string shlib, string prefix) {
  	gSystem->Load("libPhysics");
  	gSystem->Load(shlib.c_str());

	string LABEL  = "";
	string file   = sample;
	string outfile = prefix+"/"+sample;

	// log file 
	TString log=sample+".skim.log";
	ofstream f_log;
	f_log.open(log.Data());
	f_log << "Skimming file: " << sample << " with cuts: " << endl;


	// cuts --------------------------------------------
	  //GENERIC
  	std::ostringstream cutStream;
	cutStream       << " " 	  
	  //	  << "misc.MT2 >=50"                                     << "&&" 
			<< "misc.MET>=30"                                      << "&&"
			<< "misc.HT > 750 "                                    << "&&"
			<< "misc.Jet0Pass ==1"                                 << "&&"
			<< "misc.Jet1Pass ==1"                                 << "&&"
			<< "misc.SecondJPt  >100"                              << "&&"
			<< "misc.PassJetID ==1"                                << "&&"
			<< "misc.Vectorsumpt < 70"                             << "&&"
	  //<< "((misc.MinMetJetDPhi >0.3&&NBJets==0)||NBJets>=1)" << "&&"
	  << "misc.MinMetJetDPhi4 >0.3 " << "&&"
	  // Lepton Veto
//	  << "(NEles==0 || ele[0].lv.Pt()<10)"                   << "&&"
//	  << "(NMuons==0 || muo[0].lv.Pt()<10)"                  << "&&"
// Lepton Skim
//	  << "(ele[0].lv.Pt()>10 || muo[0].lv.Pt()>10)"          << "&&"
// LowMT2 ----------------------------
//	  << "misc.LeadingJPt >150"                              << "&&"
//	  << "NBJets >0"                                         << "&&"
//	  << "NJetsIDLoose >=4"                                  << "&&"
// -----------------------------------
	<< "NJetsIDLoose40 >=2"                              ;//  << "&&"
// Photons
//          << "(GenZ[0].Pt()>100)"                                << "&&"
//          << "(GenDiLeptPt(0,10,0,1000,true)>=100||GenPhoton[0].Pt()>=100)"  << "&&"
//          << "NPhotons >0"                                      << "&&"
// Noise
//	  << "misc.HBHENoiseFlag == 0"                           << "&&"
//	  << "misc.CSCTightHaloID==0"                            << "&&"
//	  << "misc.CrazyHCAL==0";



	////////////Hadronic
	/*
        std::ostringstream cutStream;
        cutStream       << " "
          //        << "misc.MT2 >=50"                                     << "&&"
	  //		<< "misc.MET>=30"                                      << "&&"
			<< "NJetsIDLoose40 >=3"                                  << "&&"
			<< "misc.HT > 750 "                                    << "&&"
                        << "misc.Jet0Pass ==1"                                 << "&&"
                        << "misc.Jet1Pass ==1"                                 << "&&"
                        << "misc.SecondJPt  >100"                              << "&&"
                        << "misc.PassJetID ==1";//                                << "&&"
	//<< "misc.Vectorsumpt < 70"                             << "&&"
	//<< "((misc.MinMetJetDPhi >0.3&&NBJets==0)||NBJets>=1)" << "&&"
          // Lepton Veto
	  //<< "(NEles==0 || ele[0].lv.Pt()<10)"                   << "&&"
	  //	<< "(NMuons==0 || muo[0].lv.Pt()<10)"                  << "&&"
          // Lepton Skim
          //        << "(ele[0].lv.Pt()>10 || muo[0].lv.Pt()>10)"          << "&&"
          // LowMT2 ----------------------------
          //        << "misc.LeadingJPt >150"                              << "&&"
          //        << "NBJets >0"                                         << "&&"
          //        << "NJetsIDLoose >=4"                                  << "&&"
          // -----------------------------------
	//                  << "NJetsIDLoose40 >=3"                                  << "&&"
          // Photons
          //        << "GenPhoton[0].Pt()>200"                             << "&&"
          //        << "NPhotons >0"                                       << "&&"
          // Noise
	  //		<< "misc.HBHENoiseFlagIso == 0"                           << "&&"
	  //			<< "misc.CSCTightHaloID==0"                            << "&&"
	  //		<< "misc.CrazyHCAL==0";
	  */

        ///////////// LEPTONIC, for ZInv
	/*
	std::ostringstream cutStream;
        cutStream       << " "
          //<< "misc.MT2 >= 150"                         << "&&"
          //<< "misc.MET>=30"                            << "&&"
                        << "misc.HT > 600 "                          << "&&"
          //      << "misc.caloHT50_ID >600 "                  << "&&"
          //<< "misc.Jet0Pass ==1"                       << "&&"
          //            << "misc.Jet1Pass ==1"                       << "&&"
          //    << "misc.SecondJPt  >100"                    << "&&"
          //<< "misc.PassJetID ==1"                      << "&&"
          //<< "misc.Vectorsumpt < 70"                   << "&&"
          //    << "misc.MinMetJetDPhi >0.3"                 << "&&"
                        << "misc.HBHENoiseFlag == 0"                 << "&&"
                        << "misc.CrazyHCAL==0"                       << "&&"
                        << "misc.CSCTightHaloID==0"                            << "&&"
                        << "( GetGenVPt(23)>1 || WDecayMode()!=0 || (NEles !=0 || NMuons!=0) )"                << "&&"
          // LowMT2 ----------------------------
          //        << "misc.LeadingJPt >150"                    << "&&"
          //        << "NBJets >0"                               << "&&"
          //        << "NJetsIDLoose >=4"                        << "&&"
          // -----------------------------------
                        << "NJetsIDLoose >=3";
	*/
	
	TString basecut = cutStream.str();
	string  SEL= "("+basecut+")";
	cout << "Skimming with sel: " << SEL << endl;
	TString cuts_log = basecut.ReplaceAll("&&", "\n");
	f_log << cuts_log << endl;
	//  --------------------------------------------


	// files ---------------------------------------
	TFile *_file0 = TFile::Open( (file).c_str()); 
	TTree * t = (TTree*) _file0->Get("MassTree");

	TH1F* hists[100];
	TH2F* hists2D[100]; 
	TIter next(_file0->GetListOfKeys());
	TKey *key;
	int count1D=0, count2D=0;
	while ((key=(TKey*)next())) { 
	  TString className =  key->GetClassName();
	  cout << className << "-" << endl;
	  if(className=="TH1F"){
	    cout << key->GetName() << endl;
	    hists[count1D++] = (TH1F*) _file0->Get( key->GetName() ); 
	  }
	  else if(className=="TH2F"){
	    cout << key->GetName() << endl;
	    hists2D[count2D++] = (TH2F*) _file0->Get( key->GetName() );
	  }
	}

	t->SetMaxTreeSize(19000000000);
	TFile*out = TFile::Open( (outfile).c_str(),"RECREATE");
	TTree *tc = t->CopyTree(SEL.c_str());   
	int nentries = tc->GetEntries();
	f_log << "skimmed tree has " << nentries << " entries." <<endl;
	f_log.close();
	out->Write();
	
	for(int h=0; h<count1D; h++){
	  cout << "Writing " << hists[h]->GetName() << endl;
	  hists[h]->Write();
	}
	for(int h=0; h<count2D; h++){
	  cout << "Writing " << hists2D[h]->GetName() << endl;
	  hists2D[h]->Write();
	}

	out->Close();
	_file0->Close();
	cout << "Result file: " << outfile << endl;
	// -------------------------------------------------

}
