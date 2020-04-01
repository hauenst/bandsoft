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
	TTree * outTree = new TTree("electrons","calibrated ToF events");
	double 	output_p_e       		 = 0;
	double 	output_theta_e   		 = 0;
	double 	output_phi_e     		 = 0;
	double 	output_q         		 = 0;
	double 	output_theta_q   		 = 0;
	double 	output_phi_q     		 = 0;
	double 	output_Q2        		 = 0;
	double 	output_nu        		 = 0;
	double 	output_xB        		 = 0;
	double 	output_W2        		 = 0;
	double 	output_EoP       		 = 0;
	double 	output_STTime    		 = 0;
	int 	output_sector_e	  		 = 0;
	double 	output_vrt_x_e   		 = 0;
	double 	output_vrt_y_e   		 = 0;
	double 	output_vrt_z_e   		 = 0;
	double 	output_t_e       		 = 0;
	double 	output_beta_e    		 = 0;
	double 	output_chi2pid_e 		 = 0;
	double	output_Ebeam			 = 0;
	outTree->Branch("p_e",		&output_p_e       		,"p_e/D");
	outTree->Branch("theta_e",	&output_theta_e   		,"theta_e/D");
	outTree->Branch("phi_e",	&output_phi_e     		,"phi_e/D");
	outTree->Branch("q",		&output_q         		,"q/D");
	outTree->Branch("theta_q",	&output_theta_q   		,"theta_q/D");
	outTree->Branch("phi_q",	&output_phi_q     		,"phi_q/D");
	outTree->Branch("Q2",		&output_Q2        		,"Q2/D");
	outTree->Branch("nu",		&output_nu        		,"nu/D");
	outTree->Branch("xB",		&output_xB        		,"xB/D");
	outTree->Branch("W2",		&output_W2        		,"W2/D");
	outTree->Branch("EoP",		&output_EoP       		,"EoP/D");
	outTree->Branch("STTime",	&output_STTime    		,"STTime/D");
	outTree->Branch("sector_e",	&output_sector_e	  	,"sector_e/I");
	outTree->Branch("vrt_x_e",	&output_vrt_x_e   		,"vrt_x_e/D");
	outTree->Branch("vrt_y_e",	&output_vrt_y_e   		,"vrt_y_e/D");
	outTree->Branch("vrt_z_e",	&output_vrt_z_e   		,"vrt_z_e/D");
	outTree->Branch("t_e",		&output_t_e       		,"t_e/D");
	outTree->Branch("beta_e",	&output_beta_e    		,"beta_e/D");
	outTree->Branch("chi2pid_e",	&output_chi2pid_e 		,"chi2pid_e/D");
	outTree->Branch("Ebeam",	&output_Ebeam			,"Ebeam/D");

	// Define 4D binned ToF plots and bin edges and kin cuts
	double const min_Q2 		= 2;
	double const max_Q2 		= 10;
 

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
		output_Ebeam = fixed_Ebeam;

		double p_e       		 = 0;
		double theta_e   		 = 0;
		double phi_e     		 = 0;
		double q         		 = 0;
		double theta_q   		 = 0;
		double phi_q     		 = 0;
		double Q2        		 = 0;
		double nu        		 = 0;
		double xB        		 = 0;
		double W2        		 = 0;
		double EoP       		 = 0;
		double STTime    		 = 0;
		int sector_e	  		 = 0;
		double vrt_x_e   		 = 0;
		double vrt_y_e   		 = 0;
		double vrt_z_e   		 = 0;
		double t_e       		 = 0;
		double beta_e    		 = 0;
		double chi2pid_e 		 = 0;
		inTree->SetBranchAddress("p_e",		&p_e       		 );
		inTree->SetBranchAddress("theta_e",	&theta_e   		 );
		inTree->SetBranchAddress("phi_e",	&phi_e     		 );
		inTree->SetBranchAddress("q",		&q         		 );
		inTree->SetBranchAddress("theta_q",	&theta_q   		 );
		inTree->SetBranchAddress("phi_q",	&phi_q     		 );
		inTree->SetBranchAddress("Q2",		&Q2        		 );
		inTree->SetBranchAddress("nu",		&nu        		 );
		inTree->SetBranchAddress("xB",		&xB        		 );
		inTree->SetBranchAddress("W2",		&W2        		 );
		inTree->SetBranchAddress("EoP",		&EoP       		 );
		inTree->SetBranchAddress("STTime",	&STTime    		 );
		inTree->SetBranchAddress("sector_e",	&sector_e  		 );
		inTree->SetBranchAddress("vrt_x_e",	&vrt_x_e   		 );
		inTree->SetBranchAddress("vrt_y_e",	&vrt_y_e   		 );
		inTree->SetBranchAddress("vrt_z_e",	&vrt_z_e   		 );
		inTree->SetBranchAddress("t_e",		&t_e       		 );
		inTree->SetBranchAddress("beta_e",	&beta_e    		 );
		inTree->SetBranchAddress("chi2pid_e",	&chi2pid_e 		 );

		cout << "Working on file: " << argv[i] << "\n";

		for( int ev = 0 ; ev < inTree->GetEntries() ; ev++ ){
			if( ev % 1000000 == 0 ) cout << "\ton event " << ev << "\n";

			p_e       		 = 0;
			theta_e   		 = 0;
			phi_e     		 = 0;
			q         		 = 0;
			theta_q   		 = 0;
			phi_q     		 = 0;
			Q2        		 = 0;
			nu        		 = 0;
			xB        		 = 0;
			W2        		 = 0;
			EoP       		 = 0;
			STTime    		 = 0;
			sector_e  		 = 0;
			vrt_x_e   		 = 0;
			vrt_y_e   		 = 0;
			vrt_z_e   		 = 0;
			t_e       		 = 0;
			beta_e    		 = 0;
			chi2pid_e 		 = 0;
			inTree->GetEntry(ev);

			// Recalculate quantities based on fixed beam energy
			TVector3 beamVec(0,0,fixed_Ebeam);
			TVector3 eVec;	eVec.SetMagThetaPhi(p_e,theta_e,phi_e);
			TVector3 qVec;	qVec = beamVec - eVec;
			q 	= qVec.Mag();
			theta_q = qVec.Theta();
			phi_q 	= qVec.Phi();
			nu 	= fixed_Ebeam - sqrt( p_e*p_e + mE*mE );
			Q2 	= q*q - nu*nu;
			xB	= Q2 / (2.*mP*nu);
			W2	= mP*mP - Q2 + 2*nu*mP;

			// Save and fill variables
			output_p_e       		 = p_e       	 ;
			output_theta_e   		 = theta_e   	 ;
			output_phi_e     		 = phi_e     	 ;
			output_q         		 = q         	 ;
			output_theta_q   		 = theta_q   	 ;
			output_phi_q     		 = phi_q     	 ;
			output_Q2        		 = Q2        	 ;
			output_nu        		 = nu        	 ;
			output_xB        		 = xB        	 ;
			output_W2        		 = W2        	 ;
			output_EoP       		 = EoP       	 ;
			output_STTime    		 = STTime    	 ;
			output_sector_e	  		 = sector_e	 ;
			output_vrt_x_e   		 = vrt_x_e   	 ;
			output_vrt_y_e   		 = vrt_y_e   	 ;
			output_vrt_z_e   		 = vrt_z_e   	 ;
			output_t_e       		 = t_e       	 ;
			output_beta_e    		 = beta_e    	 ;
			output_chi2pid_e 		 = chi2pid_e 	 ;
			outTree->Fill();
			
		} // end loop over inclusive e-
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
