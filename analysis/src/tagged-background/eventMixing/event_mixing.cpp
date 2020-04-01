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
	if (argc != 4){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./event_mixing [outputRootfile] [inputElectronSkim] [inputNeutronSkim]\n";
		return -1;
	}

	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("mixed","event-mixed background");
	double  output_fixed_Ebeam		 = 0;
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
	double 	output_dL       		 = 0; 
	double 	output_theta_n  		 = 0; 
	double 	output_phi_n    		 = 0; 
	double 	output_ToF      		 = 0; 
	double 	output_beta     		 = 0; 
	double 	output_p_n      		 = 0; 
	double 	output_phi_nq   		 = 0; 
	double 	output_theta_nq 		 = 0; 
	double 	output_E_n      		 = 0; 
	double 	output_CosTheta_nq 		 = 0;
	double 	output_Xp       		 = 0; 
	double 	output_Wp       		 = 0; 
	double 	output_As       		 = 0; 
	outTree->Branch("Ebeam",	&output_fixed_Ebeam		,"Ebeam/D");
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
	outTree->Branch("dL",		&output_dL       		,"dL/D");
	outTree->Branch("theta_n",	&output_theta_n  		,"theta_n/D");
	outTree->Branch("phi_n",	&output_phi_n    		,"phi_n/D");
	outTree->Branch("ToF",		&output_ToF      		,"ToF/D");
	outTree->Branch("beta",		&output_beta     		,"beta/D");	     
	outTree->Branch("p_n",		&output_p_n      		,"p_n/D");
	outTree->Branch("phi_nq",	&output_phi_nq   		,"phi_nq/D");
	outTree->Branch("theta_nq",	&output_theta_nq 		,"theta_nq/D");
	outTree->Branch("E_n",		&output_E_n      		,"E_n/D");
	outTree->Branch("CosTheta_nq",	&output_CosTheta_nq 		,"CosTheta_nq/D");
	outTree->Branch("Xp",		&output_Xp       		,"Xp/D");
	outTree->Branch("Wp",		&output_Wp       		,"Wp/D");
	outTree->Branch("As",		&output_As       		,"As/D");

	LoadGlobalShift();
	LoadRunByRunShift();
	LoadGlobalShift2ndIter();
	TRandom3 * myRand = new TRandom3(0);
		
	// Define 4D binned ToF plots and bin edges and kin cuts
	double const cut_lo_MeVee	= 5;
	double const min_Q2 		= 2;
	double const max_Q2 		= 10;
	double const min_CosTheta_nq 	= -1;
	double const max_CosTheta_nq 	= -0.8;
	int nBins_Q2		= 1;
	int nBins_CosTheta_nq 	= 1;
	TH1D *** hToF_4D	= new TH1D**[nBins_Q2];
	for( int bin_Q2 = 0; bin_Q2 < nBins_Q2 ; bin_Q2++){
		hToF_4D[bin_Q2]	= new TH1D*[nBins_CosTheta_nq];
		for( int bin_CosTheta_nq = 0; bin_CosTheta_nq < nBins_CosTheta_nq ; bin_CosTheta_nq++){
			hToF_4D[bin_Q2][bin_CosTheta_nq] = new TH1D(Form("hToF_%i_%i",bin_Q2,bin_CosTheta_nq),"",1600,-100,300);
		}
	}
	// Define background and signal edges
	const double bkgrd_min = -50;
	const double bkgrd_max = 0;
	const double signal_min = 13;
	const double signal_max = 63;
	
	// Loop over the neutron and electron skim files given
	TFile * inFile_e = new TFile(argv[2]);
	TFile * inFile_n = new TFile(argv[3]);
	TTree * inTree_e = (TTree*)inFile_e->Get("electrons");
	TTree * inTree_n = (TTree*)inFile_n->Get("neutrons");
	double input_p_e       		 = 0;
	double input_theta_e   		 = 0;
	double input_phi_e     		 = 0;
	double input_fixed_Ebeam		 = 0;
	double input_dL       		 = 0; 
	double input_theta_n  		 = 0; 
	double input_phi_n    		 = 0; 
	inTree_e->SetBranchAddress("p_e",	&input_p_e       		 );
	inTree_e->SetBranchAddress("theta_e",	&input_theta_e   		 );
	inTree_e->SetBranchAddress("phi_e",	&input_phi_e     		 );
	inTree_e->SetBranchAddress("Ebeam",	&input_fixed_Ebeam		 );
	inTree_n->SetBranchAddress("dL",	&input_dL       		 ); 
	inTree_n->SetBranchAddress("theta_n",	&input_theta_n  		 ); 
	inTree_n->SetBranchAddress("phi_n",	&input_phi_n    		 ); 

	double n_theta	[inTree_n->GetEntries()];
	double n_phi	[inTree_n->GetEntries()];
	double n_dL	[inTree_n->GetEntries()];
	double e_theta	[inTree_e->GetEntries()/100];
	double e_phi	[inTree_e->GetEntries()/100];
	double e_mom	[inTree_e->GetEntries()/100];
	double e_beam	[inTree_e->GetEntries()/100];
	// Read entire tree into memory due to issue with jumping around in tree as we read it..
	cout << "Saving neutrons to mem...\n";
	for( int neutron = 0 ; neutron < inTree_n->GetEntries() ; neutron++ ){
		input_theta_n  		 = 0; 
		input_phi_n    		 = 0; 
		input_dL		 = 0;
		
		inTree_n->GetEntry(neutron);
		n_theta		[neutron]	= input_theta_n		;
		n_phi		[neutron]	= input_phi_n		;
		n_dL		[neutron]	= input_dL		;
	}
	cout << "Saving electrons to mem...\n";
	for( int electron = 0 ; electron < inTree_e->GetEntries() ; electron++ ){
		input_p_e       		 = 0;
		input_theta_e   		 = 0;
		input_phi_e     		 = 0;
		input_fixed_Ebeam                = 0;
		
		inTree_e->GetEntry(electron);
		if( electron > inTree_e->GetEntries()/100 ) break;
		e_theta		[electron]	= input_theta_e		;
		e_phi		[electron]	= input_phi_e		;
		e_mom		[electron]	= input_p_e		;
		e_beam		[electron]	= input_fixed_Ebeam	;
	}
	
	// Loop over all neutron events
	for( int neutron = 0 ; neutron < inTree_n->GetEntries() ; neutron++ ){
		if( neutron % 1000 == 0 ) cout << "working on neutron " << neutron << "\n";

		int nPairs = 0; // counts how many pairs were made for this neutron
		int nFails = 0; // counts how many times a pair was tried for
		while( nPairs < 100 && nFails < 100 ){
			int electron = myRand->Rndm() * inTree_e->GetEntries()/100;
			
			// Read in the stored variables
			double fixed_Ebeam 	= e_beam	[electron];
			double p_e		= e_mom		[electron];
			double theta_e		= e_theta	[electron];
			double phi_e		= e_phi		[electron];
			double theta_n		= n_theta	[neutron];
			double phi_n		= n_phi		[neutron];
			double dL		= n_dL		[neutron];
		
			// Redraw my ToF and calculate all quantities to save
			double ToF = myRand->Rndm()*( signal_max - signal_min ) + signal_min ;
			double beta = dL / (ToF*cAir);
			double p_n = mN / sqrt( 1./pow(beta,2) - 1. );
			
			// Create vectors for calculating angles
			TVector3 beamVec(0,0,fixed_Ebeam);
			TVector3 eVec;	eVec.SetMagThetaPhi(p_e,theta_e,phi_e);
			TVector3 qVec;	qVec = beamVec - eVec;
			TVector3 nVec;	nVec.SetMagThetaPhi(p_n,theta_n,phi_n);

			double q 	= qVec.Mag();
			double theta_q  = qVec.Theta();
			double phi_q 	= qVec.Phi();
			double nu 	= fixed_Ebeam - sqrt( p_e*p_e + mE*mE );
			double Q2 	= q*q - nu*nu;
			double xB	= Q2 / (2.*mP*nu);
			double W2	= mP*mP - Q2 + 2*nu*mP;
			double E_n 	= sqrt( mN*mN + p_n*p_n );

			// Calculate the nq angles
			TVector3 norm_scatter = qVec.Cross( beamVec );
			norm_scatter 	= norm_scatter.Unit();
			TVector3 norm_reaction = qVec.Cross( nVec );
			norm_reaction 	= norm_reaction.Unit();
			double phi_nq 	= norm_scatter.Angle( norm_reaction );
			double theta_nq = nVec.Angle( qVec );
			double CosTheta_nq = cos(theta_nq);
			TVector3 direction = norm_scatter.Cross(norm_reaction);
			if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
			}
			else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
				phi_nq *= (-1);
			}
			double W_primeSq = mD*mD - Q2 + mN*mN + 2.*mD*(nu-E_n) - 2.*nu*E_n + 2.*q*p_n*cos(theta_nq);
			double Wp = sqrt(W_primeSq);
			double Xp = Q2/(2.*( nu*(mD-E_n) + p_n*q*CosTheta_nq));
			double As = (E_n - p_n*CosTheta_nq)/mN;

			// Ask if passes cuts, and if so, count valid pair and fill histogram/tree
			if( Q2 > min_Q2 && Q2 < max_Q2 && CosTheta_nq > min_CosTheta_nq && CosTheta_nq < max_CosTheta_nq ){
				// Valid pair
				nFails = 0;
				nPairs++;
				output_fixed_Ebeam 	= fixed_Ebeam 	;
				output_p_e		= p_e		;
				output_theta_e		= theta_e	;	
				output_phi_e		= phi_e		;
				output_theta_n		= theta_n	;	
				output_phi_n		= phi_n		;
				output_dL		= dL		;
				output_ToF 		= ToF 		;
				output_beta 		= beta 		;
				output_p_n 		= p_n 		;
				output_q 		= q 		;
				output_theta_q  	= theta_q  	;
				output_phi_q 		= phi_q 	;	
				output_nu 		= nu 		;
				output_Q2 		= Q2 		;
				output_xB		= xB		;
				output_W2		= W2		;
				output_E_n 		= E_n 		;
				output_phi_nq 		= phi_nq 	;	
				output_theta_nq 	= theta_nq 	;
				output_CosTheta_nq 	= CosTheta_nq 	;
				output_Wp 		= Wp 		;
				output_Xp 		= Xp 		;
				output_As 		= As 		;

				outTree->Fill();
				

			}
			else{ nFails++; }
			
			
		} // end loop over finding electrons

	}

	outFile->cd();
	for( int bin_Q2 = 0; bin_Q2 < nBins_Q2 ; bin_Q2++){
		for( int bin_CosTheta_nq = 0; bin_CosTheta_nq < nBins_CosTheta_nq ; bin_CosTheta_nq++){
			hToF_4D[bin_Q2][bin_CosTheta_nq] -> Write();
		}
	}
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
