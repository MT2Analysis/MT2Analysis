#include "helper/Hemisphere.hh" 
#include "TLorentzVector.h"

using namespace std;
void HemisphereDBG(){

	struct HemiObj {vector<TString> type; vector<int> index; vector<int> hemisphere;
	} hemiobjs;
	
	TLorentzVector p;
	std::vector<TLorentzVector> jets;
	
	TLorentzVector j(0,0,0,0);
// 	Event 1
//	j.SetPxPyPzE(  -105.87, 37.6926, -172.621, 205.979); jets.push_back(j);
//	j.SetPxPyPzE(  54.7914, -73.5522, -364.353, 375.72); jets.push_back(j);
//	j.SetPxPyPzE(  41.7905, 47.931, -33.6835, 71.9611); jets.push_back(j);
//	j.SetPxPyPzE(  -34.1219, -49.3399, -94.0637, 111.565); jets.push_back(j);

// 	Event 2
//	j.SetPxPyPzE( -68.6346, 115.285, -76.6308 ,154.511); jets.push_back(j);
//	j.SetPxPyPzE(  24.165, -104.687, -46.0788, 116.904); jets.push_back(j);
//	j.SetPxPyPzE(   68.8099, 35.8709, -14.5291, 78.947); jets.push_back(j);
//	j.SetPxPyPzE(  -35.7681, -20.7164, -23.391, 47.4939); jets.push_back(j);

// 	Event 3
//	j.SetPxPyPzE(  -0.876251, -69.699, -185.086, 197.776); jets.push_back(j);
//	j.SetPxPyPzE(  66.0559 ,  -9.59629, -102.383, 122.22); jets.push_back(j);
//	j.SetPxPyPzE(  -55.2165, 13.3919, -48.9457, 74.9926); jets.push_back(j);
//	j.SetPxPyPzE(    6.6639, 52.4509, -36.2669, 64.1154); jets.push_back(j);

// 	Event 4
///	j.SetPxPyPzE(    -35.415, -29.6071, 28.2148, 54.1006); jets.push_back(j);
///	j.SetPxPyPzE(    41.9645, -9.37466, -18.6909, 46.8856); jets.push_back(j);
///	j.SetPxPyPzE(    -4.79925, -42.5528, -36.6005, 56.3327); jets.push_back(j);
///	j.SetPxPyPzE(    36.3802, 21.9876, -46.5538, 63.0415); jets.push_back(j);
///	j.SetPxPyPzE(    -4.73678, 40.9546, 19.1024, 45.4381); jets.push_back(j);

// 	Event 5
	j.SetPxPyPzE(  97.1063, -5.30115, -74.5191, 122.519); jets.push_back(j);
	j.SetPxPyPzE(  -37.3242, -41.7401, -71.6686, 90.949); jets.push_back(j);
	j.SetPxPyPzE(  -22.4109, 49.7455, 8.75688, 55.2589); jets.push_back(j);
	j.SetPxPyPzE(  13.8803, -43.1, 53.5636, 70.1379); jets.push_back(j);


	cout << "hemisphere association for jets with:" << endl;
	for(int i=0; i<jets.size(); ++i){
		cout <<  " E "<< jets[i].E() << " Px " << jets[i].Px() << " Py " << jets[i].Py() << " Pz " << jets[i].Pz() << endl;
	}

	vector<float> px, py, pz, E;
	for(int i=0; i<jets.size(); ++i){
		px.push_back(jets[i].Px());
		py.push_back(jets[i].Py());
		pz.push_back(jets[i].Pz());
		E. push_back(jets[i].E());
		hemiobjs.index.push_back(i); hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
	}

 	Hemisphere* hemisp = new Hemisphere(px, py, pz, E, 2, 3);
	hemisp->SetDebug(1);
  	vector<int> grouping = hemisp->getGrouping();
	int nloop  = hemisp->GetNumLoop();

	
  	TLorentzVector pseudojet1(0.,0.,0.,0.);
  	TLorentzVector pseudojet2(0.,0.,0.,0.);
	
	for(int i=0; i<px.size(); ++i){
	if(grouping[i]==1){
		pseudojet1.SetPx(pseudojet1.Px() + px[i]);
		pseudojet1.SetPy(pseudojet1.Py() + py[i]);
		pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
		pseudojet1.SetE( pseudojet1.E()  + E[i]);	
		hemiobjs.hemisphere[i]=1;
	}else if(grouping[i] == 2){
		pseudojet2.SetPx(pseudojet2.Px() + px[i]);
		pseudojet2.SetPy(pseudojet2.Py() + py[i]);
		pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
		pseudojet2.SetE( pseudojet2.E()  + E[i]);
		hemiobjs.hemisphere[i]=2;
	}
	}
	delete hemisp;


	cout << "-------------------------------" << endl;
	cout << "nloops:      " << nloop << endl;
	cout << "association: (";
       	for(int i=0; i<hemiobjs.hemisphere.size(); ++i){
		cout << hemiobjs.hemisphere[i] << ", ";
	}
	cout << ")" << endl;	
	cout << "-------------------------------" << endl;

}
