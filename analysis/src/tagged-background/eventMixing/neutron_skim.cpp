#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unordered_map>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TVectorT.h"
#include "TFitResultPtr.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "../../deuteron_dis/include/constants.h"

using namespace std;

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
void LoadGlobalShift();
void LoadRunByRunShift();
void LoadGlobalShift2ndIter();
double FADC_GLOBSHIFT[600] = {0.};
double FADC_RUNBYRUNSHIFT[10000] = {0.};
double FADC_GLOBSHIFT2NDITER[600] = {0.};
char* getRunNumber( char* parse );
double getBeamEnergy( int runNum );


int main(int argc, char ** argv){
	if (argc < 2){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./w_plot [outputRootfile] [inputDatafiles]\n";
		return -1;
	}

	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("neutrons","skimmed neutron events");
	int	output_thisRun			 = 0;
	int 	output_nHits     		 = 0;
	int 	output_sector    		 = 0;
	int 	output_layer     		 = 0;
	int 	output_component 		 = 0;
	double 	output_adcLcorr  		 = 0;
	double 	output_adcRcorr  		 = 0;
	double 	output_meantimeFadc		 = 0;
	double 	output_meantimeTdc		 = 0;
	double 	output_difftimeFadc		 = 0;
	double 	output_difftimeTdc		 = 0;
	double 	output_x        		 = 0; 
	double 	output_y        		 = 0; 
	double 	output_z        		 = 0; 
	double 	output_dL       		 = 0; 
	double 	output_theta_n  		 = 0; 
	double 	output_phi_n    		 = 0; 
	double 	output_ToF      		 = 0; 
	double 	output_beta     		 = 0; 
	double 	output_p_n      		 = 0; 
	double 	output_E_n      		 = 0; 
	outTree->Branch("thisRun",	&output_thisRun			,"thisRun/I");
	outTree->Branch("nHits",	&output_nHits     		,"nHits/I");
	outTree->Branch("sector",	&output_sector    		,"sector/I");
	outTree->Branch("layer",	&output_layer     		,"layer/I");
	outTree->Branch("component",	&output_component 		,"component/I");
	outTree->Branch("adcLcorr",	&output_adcLcorr  		,"adcLcorr/D");
	outTree->Branch("adcRcorr",	&output_adcRcorr  		,"adcRcorr/D");
	outTree->Branch("meantimeFadc",	&output_meantimeFadc		,"meantimeFadc/D");
	outTree->Branch("meantimeTdc",	&output_meantimeTdc		,"meantimeTdc/D");
	outTree->Branch("difftimeFadc",	&output_difftimeFadc		,"difftimeFadc/D");
	outTree->Branch("difftimeTdc",	&output_difftimeTdc		,"difftimeTdc/D");
	outTree->Branch("x",		&output_x        		,"x/D");
	outTree->Branch("y",		&output_y        		,"y/D");
	outTree->Branch("z",		&output_z        		,"z/D");
	outTree->Branch("dL",		&output_dL       		,"dL/D");
	outTree->Branch("theta_n",	&output_theta_n  		,"theta_n/D");
	outTree->Branch("phi_n",	&output_phi_n    		,"phi_n/D");
	outTree->Branch("ToF",		&output_ToF      		,"ToF/D");
	outTree->Branch("beta",		&output_beta     		,"beta/D");	     
	outTree->Branch("p_n",		&output_p_n      		,"p_n/D");
	outTree->Branch("E_n",		&output_E_n      		,"E_n/D");

	// Define 4D binned ToF plots and bin edges and kin cuts
	double const cut_lo_MeVee	= 5;
	// Define background edges
	const double bkgrd_min = -50;
	const double bkgrd_max = 0;
	const double signal_min = 13;
	const double signal_max = 63;


	LoadGlobalShift();
	LoadRunByRunShift();
	LoadGlobalShift2ndIter();

	// Loop over all the files that are given to me
	for( int i = 2 ; i < argc ; i++ ){
		TFile * inFile = new TFile(argv[i]);
		if (!(inFile->GetListOfKeys()->Contains("skim"))){
			cerr << "File has no entries\n";
			return -2;
		}
		TTree * inTree = (TTree*)inFile->Get("skim");

		int thisRun = atoi(getRunNumber(argv[i]));
		// Get beam energy from this run number:
		double fixed_Ebeam = getBeamEnergy(thisRun);
		double STTime			 = 0;
		int nHits     			 = 0;
		int sector    			 = 0;
		int layer     			 = 0;
		int component 			 = 0;
		double adcLcorr  		 = 0;
		double adcRcorr  		 = 0;
		double meantimeFadc		 = 0;
		double meantimeTdc		 = 0;
		double difftimeFadc		 = 0;
		double difftimeTdc		 = 0;
		double x        		 = 0; 
		double y        		 = 0; 
		double z        		 = 0; 
		double dL       		 = 0; 
		double theta_n  		 = 0; 
		double phi_n    		 = 0; 
		double ToF      		 = 0; 
		double beta     		 = 0; 
		double p_n      		 = 0; 
		double E_n      		 = 0; 
		int nADC     			 = 0; 
		int nTDC     			 = 0; 
		double phaseCorr		 = 0; 
		double adcLraw  		 = 0; 
		double adcRraw  		 = 0; 
		double tTdcLraw 		 = 0; 
		double tTdcRraw 		 = 0; 
		double tFadcLraw		 = 0; 
		double tFadcRraw		 = 0; 
		double byHand_adcL		 = 0; 
		double byHand_adcR		 = 0; 
		double byHand_meantimeFadc 	 = 0;
		double byHand_meantimeTdc  	 = 0;
		double byHand_difftimeFadc 	 = 0;
		double byHand_difftimeTdc 	 = 0;
		double byHand_x 		 = 0; 
		double byHand_y 		 = 0; 
		double byHand_z 		 = 0; 
		double byHand_dL		 = 0;
		inTree->SetBranchAddress("STTime",	&STTime			);
		inTree->SetBranchAddress("nHits",	&nHits     		 );
		inTree->SetBranchAddress("sector",	&sector    		 );
		inTree->SetBranchAddress("layer",	&layer     		 );
		inTree->SetBranchAddress("component",	&component 		 );
		inTree->SetBranchAddress("adcLcorr",	&adcLcorr  		 );
		inTree->SetBranchAddress("adcRcorr",	&adcRcorr  		 );
		inTree->SetBranchAddress("meantimeFadc",&meantimeFadc		 );
		inTree->SetBranchAddress("meantimeTdc",	&meantimeTdc		 );
		inTree->SetBranchAddress("difftimeFadc",&difftimeFadc		 );
		inTree->SetBranchAddress("difftimeTdc",	&difftimeTdc		 );
		inTree->SetBranchAddress("x",		&x        		 ); 
		inTree->SetBranchAddress("y",		&y        		 ); 
		inTree->SetBranchAddress("z",		&z        		 ); 
		inTree->SetBranchAddress("dL",		&dL       		 ); 
		inTree->SetBranchAddress("theta_n",	&theta_n  		 ); 
		inTree->SetBranchAddress("phi_n",	&phi_n    		 ); 
		inTree->SetBranchAddress("ToF",		&ToF      		 ); 
		inTree->SetBranchAddress("beta",	&beta     		 ); 
		inTree->SetBranchAddress("p_n",		&p_n      		 ); 
		inTree->SetBranchAddress("E_n",		&E_n      		 ); 
		inTree->SetBranchAddress("nADC",	&nADC     		 ); 
		inTree->SetBranchAddress("nTDC",	&nTDC     		 ); 
		inTree->SetBranchAddress("phaseCorr",	&phaseCorr		 ); 
		inTree->SetBranchAddress("adcLraw",	&adcLraw  		 ); 
		inTree->SetBranchAddress("adcRraw",	&adcRraw  		 ); 
		inTree->SetBranchAddress("tTdcLraw",	&tTdcLraw 		 ); 
		inTree->SetBranchAddress("tTdcRraw",	&tTdcRraw 		 ); 
		inTree->SetBranchAddress("tFadcLraw",	&tFadcLraw		 ); 
		inTree->SetBranchAddress("tFadcRraw",	&tFadcRraw		 ); 
		inTree->SetBranchAddress("byHand_adcL",	&byHand_adcL		 ); 
		inTree->SetBranchAddress("byHand_adcR",	&byHand_adcR		 ); 
		inTree->SetBranchAddress("byHand_meantimeFadc",	&byHand_meantimeFadc 	 );
		inTree->SetBranchAddress("byHand_meantimeTdc",	&byHand_meantimeTdc  	 );
		inTree->SetBranchAddress("byHand_difftimeFadc",	&byHand_difftimeFadc 	 );
		inTree->SetBranchAddress("byHand_difftimeTdc",	&byHand_difftimeTdc 	 );
		inTree->SetBranchAddress("byHand_x",	&byHand_x 		 ); 
		inTree->SetBranchAddress("byHand_y",	&byHand_y 		 ); 
		inTree->SetBranchAddress("byHand_z",	&byHand_z 		 ); 
		inTree->SetBranchAddress("byHand_dL",	&byHand_dL		 ); 

		cout << "Working on file: " << argv[i] << "\n";

		for( int ev = 0 ; ev < inTree->GetEntries() ; ev++ ){
			if( ev % 1000000 == 0 ) cout << "\ton event " << ev << "\n";
			STTime			 = 0;
			nHits     		 = 0;
			sector    		 = 0;
			layer     		 = 0;
			component 		 = 0;
			adcLcorr  		 = 0;
			adcRcorr  		 = 0;
			meantimeFadc		 = 0;
			meantimeTdc		 = 0;
			difftimeFadc		 = 0;
			difftimeTdc		 = 0;
			x        		 = 0; 
			y        		 = 0; 
			z        		 = 0; 
			dL       		 = 0; 
			theta_n  		 = 0; 
			phi_n    		 = 0; 
			ToF      		 = 0; 
			beta     		 = 0; 
			p_n      		 = 0; 
			E_n      		 = 0; 
			nADC     		 = 0; 
			nTDC     		 = 0; 
			phaseCorr		 = 0; 
			adcLraw  		 = 0; 
			adcRraw  		 = 0; 
			tTdcLraw 		 = 0; 
			tTdcRraw 		 = 0; 
			tFadcLraw		 = 0; 
			tFadcRraw		 = 0; 
			byHand_adcL		 = 0; 
			byHand_adcR		 = 0; 
			byHand_meantimeFadc 	 = 0;
			byHand_meantimeTdc  	 = 0;
			byHand_difftimeFadc 	 = 0;
			byHand_difftimeTdc 	 = 0;
			byHand_x 		 = 0; 
			byHand_y 		 = 0; 
			byHand_z 		 = 0; 
			byHand_dL		 = 0; 

			inTree->GetEntry(ev);

				
			if( nHits == 1 && meantimeFadc!=0 && sqrt(adcLcorr*adcRcorr) > cut_lo_MeVee*3000){
				int barID = sector*100 + layer*10 + component;
				ToF = (meantimeFadc - STTime) - FADC_GLOBSHIFT[barID] - FADC_RUNBYRUNSHIFT[thisRun] - FADC_GLOBSHIFT2NDITER[barID];
				if( ToF > bkgrd_max || ToF < bkgrd_min ) continue;

				// use this neutron!!
				output_thisRun			 = thisRun	 ;
				output_nHits     		 = nHits     	 ;
				output_sector    		 = sector    	 ;
				output_layer     		 = layer     	 ;
				output_component 		 = component 	 ;
				output_adcLcorr  		 = adcLcorr  	 ;
				output_adcRcorr  		 = adcRcorr  	 ;
				output_meantimeFadc		 = meantimeFadc	 ;
				output_meantimeTdc		 = meantimeTdc	 ;
				output_difftimeFadc		 = difftimeFadc	 ;
				output_difftimeTdc		 = difftimeTdc	 ;
				output_x        		 = x        	 ; 
				output_y        		 = y        	 ; 
				output_z        		 = z        	 ; 
				output_dL       		 = dL       	 ; 
				output_theta_n  		 = theta_n  	 ; 
				output_phi_n    		 = phi_n    	 ; 
				output_ToF      		 = ToF      	 ; 
				output_beta     		 = beta     	 ; 
				output_p_n      		 = p_n      	 ; 
				output_E_n      		 = E_n      	 ; 
				outTree->Fill();
				
			}
			else{ continue; }
		} // end loop over events
		inFile->Close();
	}// end loop over files

	outFile->cd();
	outTree->Write();
	outFile->Close();



	return 0;
}

void LoadGlobalShift(){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open("global_offset_fadc-10082019.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		FADC_GLOBSHIFT[barId] = mean;
		f >> sig;
		f >> temp;
		f >> temp;
	}
	f.close();
}

char* getRunNumber( char* parse ){
	char * parse_copy = (char*) malloc( strlen(parse)+1 );
	char * parsed;

	strcpy( parse_copy, parse );
	char * loop = strtok( parse_copy, ".");
	while(loop){
		char* equals_sign = strchr(loop, '_');
		if (equals_sign){
			*equals_sign = 0;
			equals_sign++;
			parsed = (char*)malloc(strlen(equals_sign) + 1);
			strcpy(parsed, equals_sign);
		}
		loop = strtok(NULL, ".");
	}
	free(parse_copy);
	return parsed;
}

void LoadRunByRunShift(){
	ifstream f;
	int runnum;
	double pol0, height, mean, sig, temp;

	f.open("runByrun_offset_fadc-10082019.txt");
	while(!f.eof()){
		f >> runnum;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		FADC_RUNBYRUNSHIFT[runnum] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}

void LoadGlobalShift2ndIter(){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open("global_offset_fadc_2ndIter-10082019.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		FADC_GLOBSHIFT2NDITER[barId] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}

double getBeamEnergy( int runNum ){
        double thisEn = 0;

        if( runNum <= 6399 ) thisEn = 10.6;
        else{ thisEn = 10.2; }
        if( runNum == 6523 || runNum == 6524 || runNum == 6525 ) thisEn = 10.;

        return thisEn;
}
