#ifndef Utilities_HH
#define Utilities_HH

#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>

#include "TString.h"

using namespace std;

//general purpose utilities

//creates histos for counting events
// the map is: m[cut][sample]
inline void createCounterHistos(map<TString, map<TString, TH1F*> > *hm, map<TString, TString> cutLabels, TString sName){
  bool isThere=false;

  for(map<TString,TString>::iterator l=cutLabels.begin(); l!=cutLabels.end();l++){
    TString hName = "counter_"+l->first+"_"+sName;

    for(map<TString, TH1F*>::iterator h=(*hm)[ l->first].begin(); h!=(*hm)[ l->first].end(); h++){
      if(h->first == sName) isThere=true;
    }
    if(!isThere){
      (*hm)[l->first][sName] = new TH1F(hName,"",500,0,5000);
      (*hm)[l->first][sName]->Sumw2();
    }
    //cout << sName << " " << l->first << endl;
  }
  
}


//
inline void getEfficiency(TH1F* num, TH1F*den, float *eff, float *err){
  TH1F *num_c, *den_c;
  num_c = (TH1F*) num->Clone();
  TString hname = num->GetName();
  num_c->SetName( hname+"_clone");
  hname = den->GetName();
  den_c = (TH1F*) den->Clone();
  den_c->SetName( hname+"_clone");


  num_c->Rebin( num_c->GetNbinsX() );
  den_c->Rebin( den_c->GetNbinsX() );
  num_c->Divide(den_c);
  (*eff) = (float)num_c->GetBinContent(1);
  (*err) = (float)num_c->GetBinError(1);

  delete num_c;
  delete den_c;

}


inline void printLatexTable(  vector<TString> samplesName, map<TString, TString> cutLabels, vector<TString> orderedCutLabels, map<TString, map<TString, TH1F*> > counterH ){

  vector<TString> samplesList;
  cout <<"\\begin{table} \n\\begin{center} \n\\begin{tabular}{lcccccccccccccccc} \n\\hline\\hline" << endl;

  string header = "";
  for( vector<TString>::const_iterator c_i = orderedCutLabels.begin(); c_i != orderedCutLabels.end(); c_i++ ){
    TString myCut = *c_i;
    map< TString, TH1F*> myH = counterH[myCut];
    if(header==""){
      for( vector<TString>::const_iterator s_i=samplesName.begin(); s_i != samplesName.end(); s_i++)
	header += "& " + *s_i ;
      
      header += " \\\\ \\hline \\hline";
      cout << header << endl;
    }   
    string line = "";
    line += " "+ cutLabels[myCut];
    //for( map< TString, TH1F* >::const_iterator H_i = myH.begin(); H_i != myH.end(); H_i++){
    for( vector<TString>::const_iterator s_i=samplesName.begin(); s_i != samplesName.end(); s_i++){
      char tmp[50];
      TString mySample = *s_i;
      //cout << *s_i << endl;
  
      if(mySample=="MC" || mySample=="HT-Data"){
	double err;
	double integral = counterH[myCut][ mySample]->IntegralAndError(1 ,counterH[myCut][ mySample]->GetNbinsX() , err);
	sprintf( tmp, " & %.2f +- %.2f", integral, err );
      }
      else
	sprintf( tmp, " & %.2f ", counterH[myCut][ mySample]->Integral() );
      line += tmp;
    }
    cout << line <<  " \\\\ " << endl;
    
  }
  string footer = "\\hline\\hline \n\\end{tabular} \n\\end{center} \n\\end{table}";
  cout << footer << endl;

}


inline void printTwikiTable(  vector<TString> samplesName, map<TString, TString> cutLabels, vector<TString> orderedCutLabels, map<TString, map<TString, TH1F*> > counterH ){

  vector<TString> samplesList;
  //cout <<"\\begin{table} \n\\begin{center} \n\\begin{tabular}{lcccccccccccccccc} \n\\hline\\hline" << endl;

  string header = "|  ";
  for( vector<TString>::const_iterator c_i = orderedCutLabels.begin(); c_i != orderedCutLabels.end(); c_i++ ){
    TString myCut = *c_i;
    map< TString, TH1F*> myH = counterH[myCut];
    if(header=="|  "){
      for( vector<TString>::const_iterator s_i=samplesName.begin(); s_i != samplesName.end(); s_i++)
	header += " | *" + *s_i + "* ";
      
      header += " | ";
      cout << header << endl;
    }   
    string line = "";
    line += "| *"+ cutLabels[myCut]+"* ";
    //for( map< TString, TH1F* >::const_iterator H_i = myH.begin(); H_i != myH.end(); H_i++){
    for( vector<TString>::const_iterator s_i=samplesName.begin(); s_i != samplesName.end(); s_i++){
      char tmp[50];
      TString mySample = *s_i;
      //cout << *s_i << endl;
  
      if(mySample=="MC" || mySample=="HT-Data"){
	double err;
	double integral = counterH[myCut][ mySample]->IntegralAndError(1 ,counterH[myCut][ mySample]->GetNbinsX() , err);
	sprintf( tmp, " | %.2f +- %.2f", integral, err );
      }
      else
	sprintf( tmp, " | %.2f ", counterH[myCut][ mySample]->Integral() );
      line += tmp;
    }
    cout << line <<  " | " << endl;
    
  }
  //  string footer = "\\hline\\hline \n\\end{tabular} \n\\end{center} \n\\end{table}";
  //cout << footer << endl;

}


#endif
