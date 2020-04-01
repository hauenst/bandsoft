#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

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
void PrettyHist( TH1D * hist , bool sig , double weight , bool xp, int thisSlice );


int main(int argc, char ** argv){
	if (argc < 2){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./w_plot [outputRootfile] [inputDatafiles]\n";
		return -1;
	}

	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("calib","calibrated ToF events");
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
	double 	output_phi_nq   		 = 0; 
	double 	output_theta_nq 		 = 0; 
	double 	output_E_n      		 = 0; 
	double 	output_phi_en   		 = 0; 
	double 	output_CosTheta_nq 		 = 0;
	double 	output_Xp       		 = 0; 
	double 	output_Wp       		 = 0; 
	double 	output_As       		 = 0; 
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
	outTree->Branch("phi_nq",	&output_phi_nq   		,"phi_nq/D");
	outTree->Branch("theta_nq",	&output_theta_nq 		,"theta_nq/D");
	outTree->Branch("E_n",		&output_E_n      		,"E_n/D");
	outTree->Branch("phi_en",	&output_phi_en   		,"phi_en/D");
	outTree->Branch("CosTheta_nq",	&output_CosTheta_nq 		,"CosTheta_nq/D");
	outTree->Branch("Xp",		&output_Xp       		,"Xp/D");
	outTree->Branch("Wp",		&output_Wp       		,"Wp/D");
	outTree->Branch("As",		&output_As       		,"As/D");

	LoadGlobalShift();
	LoadRunByRunShift();
	LoadGlobalShift2ndIter();
	TRandom3 * myRand = new TRandom3(0);
		
	// Define histograms to be used:
	TH1D * hToF_spb = new TH1D("hToF_spb","hToF_spb",1600,-100,300);
	TH1D * hToF_b = new TH1D("hToF_b","hToF_b",1600,-100,300);
	TH1D * hToF_orig_spb = new TH1D("hToF_orig_spb","hToF_orig_spb",1600,-100,300);
	TH1D * hToF_orig_b = new TH1D("hToF_orig_b","hToF_orig_b",1600,-100,300);
	TH1D * hW_spb = new TH1D("hW_spb","hW_spb",500,0,5.0);
	TH1D * hWp_spb = new TH1D("hWp_spb","hXp_spb",50,0,1);
	TH1D * hW_b = new TH1D("hW_b","hW_b",500,0,5.0);
	TH1D * hWp_b = new TH1D("hWp_b","hXp_b",50,0,1);
	TH1D * p_n_loX = new TH1D("p_n_loX","p_n_loX",5,0.1,0.6);
	TH1D * hW_full = new TH1D("hW_full","hW_full",500,0,5.0);
	TH1D * hXp_basic = new TH1D("hXp_basic","hXp_basic",50,0,1);
	int x_slices = 10;
	TH1D ** ToF_spb_xslices 	= new TH1D*[x_slices];
	for( int thisSlice = 0 ; thisSlice < x_slices ; thisSlice++ )
		ToF_spb_xslices[thisSlice]	= new TH1D(Form("ToF_spb_xslices_%i",thisSlice)	,Form("ToF_spb_xslices_%i",thisSlice),		300,-75,225);
	

	// Define 4D binned ToF plots and bin edges and kin cuts
	double const cut_lo_Q2		= 2;
	double const cut_lo_MeVee	= 5;
	double const min_Q2 		= 2;
	double const max_Q2 		= 10;
	double const min_CosTheta_nq 	= -1;
	double const max_CosTheta_nq 	= -0.8;
	int nBins_Q2		= 1;
	int nBins_CosTheta_nq 	= 1;
	double avg_xp_cnts[200] = {0.};
	double avg_cnts[200] = {0.};
	TH1D *** hToF_4D	= new TH1D**[nBins_Q2];
	for( int bin_Q2 = 0; bin_Q2 < nBins_Q2 ; bin_Q2++){
		hToF_4D[bin_Q2]	= new TH1D*[nBins_CosTheta_nq];
		for( int bin_CosTheta_nq = 0; bin_CosTheta_nq < nBins_CosTheta_nq ; bin_CosTheta_nq++){
			hToF_4D[bin_Q2][bin_CosTheta_nq] = new TH1D(Form("hToF_%i_%i",bin_Q2,bin_CosTheta_nq),"",1600,-100,300);
		}
	}
	// Define background edges
	const double bkgrd_min = -50;
	const double bkgrd_max = 0;
	const double signal_min = 13;
	const double signal_max = 63;
	
 

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
		double phi_nq   		 = 0; 
		double theta_nq 		 = 0; 
		double E_n      		 = 0; 
		double phi_en   		 = 0; 
		double CosTheta_nq 		 = 0;
		double Xp       		 = 0; 
		double Wp       		 = 0; 
		double As       		 = 0; 
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
		inTree->SetBranchAddress("phi_nq",	&phi_nq   		 ); 
		inTree->SetBranchAddress("theta_nq",	&theta_nq 		 ); 
		inTree->SetBranchAddress("E_n",		&E_n      		 ); 
		inTree->SetBranchAddress("phi_en",	&phi_en   		 ); 
		inTree->SetBranchAddress("CosTheta_nq",	&CosTheta_nq 		 );
		inTree->SetBranchAddress("Xp",		&Xp       		 ); 
		inTree->SetBranchAddress("Wp",		&Wp       		 ); 
		inTree->SetBranchAddress("As",		&As       		 ); 
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
			phi_nq   		 = 0; 
			theta_nq 		 = 0; 
			E_n      		 = 0; 
			phi_en   		 = 0; 
			CosTheta_nq 		 = 0;
			Xp       		 = 0; 
			Wp       		 = 0; 
			As       		 = 0; 
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
			// Tagged skim -- look for 1 hit in BAND with hi energy dep
			if( nHits == 1 && meantimeFadc!=0 && sqrt(adcLcorr*adcRcorr) > cut_lo_MeVee*3000){
	
				// Recalculate electron quantities based on new beam energy
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

				// Recalculate momentum based on shifts
				int barID = sector*100 + layer*10 + component;
				ToF = (meantimeFadc - STTime) - FADC_GLOBSHIFT[barID] - FADC_RUNBYRUNSHIFT[thisRun] - FADC_GLOBSHIFT2NDITER[barID];
	
				// ONLY LOOK AT NEUTRON DIRECTION AS SOME ToF ARE UNPHYSICAL AND WILL GIVE NAN VALUES
				TVector3 nVec;	nVec.SetMagThetaPhi(1.,theta_n,phi_n);

				// Calculate the nq angles
				TVector3 norm_scatter = qVec.Cross( beamVec );
				norm_scatter 	= norm_scatter.Unit();

				TVector3 norm_reaction = qVec.Cross( nVec );
				norm_reaction 	= norm_reaction.Unit();

				phi_nq 	= norm_scatter.Angle( norm_reaction );
				theta_nq = nVec.Angle( qVec );
				CosTheta_nq = cos(theta_nq);

				TVector3 direction = norm_scatter.Cross(norm_reaction);
				if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
				}
				else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
					phi_nq *= (-1);
				}

				// Now fill our histogram to do background-subtraction
				if( Q2 > min_Q2					 &&
				    Q2 < max_Q2					 &&
				    CosTheta_nq > min_CosTheta_nq		 &&
				    CosTheta_nq	< max_CosTheta_nq		 ){
	
					hToF_4D[0][0] -> Fill( ToF );
					//int slice = (xB - 0.05)/0.1;
					//if( slice > -1 && slice < 10 ) ToF_spb_xslices[slice]->Fill( ToF );

					if( ToF > signal_min && ToF < signal_max ){
						
						beta = dL / (ToF*cAir);
						p_n = mN / sqrt( 1./pow(beta,2) - 1. );
						nVec.SetMagThetaPhi(p_n,theta_n,phi_n);
						E_n 	= sqrt( mN*mN + p_n*p_n );
						double W_primeSq = mD*mD - Q2 + mN*mN + 2.*mD*(nu-E_n) - 2.*nu*E_n + 2.*q*p_n*cos(theta_nq);
						Wp = sqrt(W_primeSq);
						Xp = Q2/(2.*( nu*(mD-E_n) + p_n*q*CosTheta_nq));
						As = (E_n - p_n*CosTheta_nq)/mN;

						hToF_spb->Fill( ToF );
						hWp_spb->Fill( Xp );

						int tof_bin = (ToF - 13)/0.25;
						avg_xp_cnts[tof_bin] += Xp;
						avg_cnts[tof_bin]++;
						//cout << ToF << " " << tof_bin << " " << avg_xp_cnts[tof_bin] << " " <<  avg_cnts[tof_bin] << "\n";

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
						output_dL       		 = dL       	 ; 
						output_theta_n  		 = theta_n  	 ; 
						output_phi_n    		 = phi_n    	 ; 
						output_ToF      		 = ToF      	 ; 
						output_beta     		 = beta     	 ; 
						output_p_n      		 = p_n      	 ; 
						output_phi_nq   		 = phi_nq   	 ; 
						output_theta_nq 		 = theta_nq 	 ; 
						output_E_n      		 = E_n      	 ; 
						output_CosTheta_nq 		 = CosTheta_nq 	 ;
						output_Xp       		 = Xp       	 ; 
						output_Wp       		 = Wp       	 ; 
						output_As       		 = As       	 ; 


						outTree->Fill();
					}
					else if( ToF > bkgrd_min && ToF < bkgrd_max ){
						hToF_orig_b->Fill( ToF );
						double fake_ToF = signal_min + myRand->Rndm()*(signal_max - signal_min);
						hToF_b->Fill( fake_ToF );

						// Calculate W and Wp based on this:
						beta = dL / (fake_ToF*cAir);
						p_n = mN / sqrt( 1./pow(beta,2) - 1. );
						TVector3 nVec;	nVec.SetMagThetaPhi(p_n,theta_n,phi_n);
						TVector3 norm_scatter = qVec.Cross( beamVec );
						norm_scatter 	= norm_scatter.Unit();

						TVector3 norm_reaction = qVec.Cross( nVec );
						norm_reaction 	= norm_reaction.Unit();

						phi_nq 	= norm_scatter.Angle( norm_reaction );
						theta_nq = nVec.Angle( qVec );
						TVector3 direction = norm_scatter.Cross(norm_reaction);
						if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
						}
						else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
							phi_nq *= (-1);
						}
						E_n = sqrt( mN*mN + p_n*p_n );
						double W_primeSq = mD*mD - Q2 + mN*mN + 2.*mD*(nu-E_n) - 2.*nu*E_n + 2.*q*p_n*cos(theta_nq);
						Wp = sqrt(W_primeSq);

						hWp_b->Fill( Xp );
						if( fabs(Xp-0.3)<0.04 ) p_n_loX->Fill(p_n);
			
					}
				}
				
				/*
				// Recalculate neutron quantities based on fixed ToF
				beta = dL / (ToF*cAir);
				p_n = mN / sqrt( 1./pow(beta,2) - 1. );
				TVector3 nVec;	nVec.SetMagThetaPhi(p_n,theta_n,phi_n);

				// Calculate the nq angles
				TVector3 norm_scatter = qVec.Cross( beamVec );
				norm_scatter 	= norm_scatter.Unit();

				TVector3 norm_reaction = qVec.Cross( nVec );
				norm_reaction 	= norm_reaction.Unit();

				phi_nq 	= norm_scatter.Angle( norm_reaction );
				theta_nq = nVec.Angle( qVec );

				TVector3 direction = norm_scatter.Cross(norm_reaction);
				if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
				}
				else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
					phi_nq *= (-1);
				}
				E_n = sqrt( mN*mN + p_n*p_n );
				double W_primeSq = mD*mD - Q2 + mN*mN + 2.*mD*(nu-E_n) - 2.*nu*E_n + 2.*q*p_n*cos(theta_nq);
				CosTheta_nq = cos(theta_nq);
				Wp = sqrt(W_primeSq);
				Xp = Q2/(2.*( nu*(mD-E_n) + p_n*q*CosTheta_nq));
				// Fill Signal+Background distributions
				if( ToF > sig_lower && ToF < sig_upper ){
					hToF_spb->Fill( ToF );
					hWp_spb->Fill( Wp );
					hW_spb->Fill( sqrt(W2) );
				}
				// Fill Backgorund distributions
				else if( ToF < bkg_upper && ToF > bkg_lower){
					hToF_orig_b->Fill( ToF );
					double fake_ToF = sig_lower + myRand->Rndm()*(sig_upper-sig_lower);
					hToF_b->Fill( fake_ToF );

					// Calculate W and Wp based on this:
					beta = dL / (fake_ToF*cAir);
					p_n = mN / sqrt( 1./pow(beta,2) - 1. );
					TVector3 nVec;	nVec.SetMagThetaPhi(p_n,theta_n,phi_n);
					TVector3 norm_scatter = qVec.Cross( beamVec );
					norm_scatter 	= norm_scatter.Unit();

					TVector3 norm_reaction = qVec.Cross( nVec );
					norm_reaction 	= norm_reaction.Unit();

					phi_nq 	= norm_scatter.Angle( norm_reaction );
					theta_nq = nVec.Angle( qVec );
					TVector3 direction = norm_scatter.Cross(norm_reaction);
					if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
					}
					else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
						phi_nq *= (-1);
					}
					E_n = sqrt( mN*mN + p_n*p_n );
					double W_primeSq = mD*mD - Q2 + mN*mN + 2.*mD*(nu-E_n) - 2.*nu*E_n + 2.*q*p_n*cos(theta_nq);
					Wp = sqrt(W_primeSq);

					hWp_b->Fill( Wp );
					hW_b->Fill( sqrt(W2) );
		
				}
				*/
				// Only fill if we have a neutron and electron
			} // end if on neutron info fill
		} // end loop over events
		inFile->Close();
	}// end loop over files

	// For each bin in Q2, Theta_nq of ToF spectrum, 
	// fit the S+B to get the normalization of total
	// background events under the S+B region
	outFile->cd();
	TF1 *** ToF_fits = new TF1**[nBins_Q2];
	for( int bin_Q2 = 0; bin_Q2 < nBins_Q2 ; bin_Q2++){
		ToF_fits[bin_Q2] = new TF1*[nBins_CosTheta_nq];
		for( int bin_CosTheta_nq = 0; bin_CosTheta_nq < nBins_CosTheta_nq ; bin_CosTheta_nq++){
			TVectorT<double> norm(2);
			hToF_4D[bin_Q2][bin_CosTheta_nq] -> Write();
			
			ToF_fits[bin_Q2][bin_CosTheta_nq] = new TF1(Form("ToF_fits_%i_%i",bin_Q2,bin_CosTheta_nq),"pol0",bkgrd_min,bkgrd_max);
			double start = hToF_4D[bin_Q2][bin_CosTheta_nq]->FindBin(bkgrd_min);
			double end = hToF_4D[bin_Q2][bin_CosTheta_nq]->FindBin(bkgrd_max);
			double nbins = end - start;
			double bkgrd_guess = hToF_4D[bin_Q2][bin_CosTheta_nq]->Integral(start,end) / nbins;
			ToF_fits[bin_Q2][bin_CosTheta_nq]->SetParameter(0,bkgrd_guess);
			hToF_4D[bin_Q2][bin_CosTheta_nq] -> Fit( ToF_fits[bin_Q2][bin_CosTheta_nq] , "QESR" );
			cout << "parameter: " << ToF_fits[bin_Q2][bin_CosTheta_nq]->GetParameter(0) << "\n";
			double signal_start = hToF_4D[bin_Q2][bin_CosTheta_nq]->FindBin(signal_min);
			double signal_end = hToF_4D[bin_Q2][bin_CosTheta_nq]->FindBin(signal_max);
			double signal_nbins = signal_end - signal_start;
			cout << "normalization: " << ToF_fits[bin_Q2][bin_CosTheta_nq]->GetParameter(0) * (signal_nbins) << "\n";
			norm[0] = ToF_fits[bin_Q2][bin_CosTheta_nq]->GetParameter(0);
			norm[1] = ToF_fits[bin_Q2][bin_CosTheta_nq]->GetParameter(0) * (signal_nbins);
			norm.Write(Form("norm_%i_%i",bin_Q2,bin_CosTheta_nq));
		}
	}
	for( int tof_bin = 0; tof_bin < 200; tof_bin++){
		avg_xp_cnts[tof_bin] /= avg_cnts[tof_bin];
		avg_cnts[tof_bin] -= 2.69112;
		hXp_basic->Fill( avg_xp_cnts[tof_bin] , avg_cnts[tof_bin] );
	}

	/*
	// Weight our background distributions correctly
	double binwidth_b = (300.-(-100.))/hToF_orig_b->GetNbinsX();
	double binwidth_spb = (300.-(-100.))/hToF_orig_spb->GetNbinsX();
	double avg_height = hToF_orig_b->Integral() / ((bkg_upper-bkg_lower)/binwidth_b);
	double weight = hToF_b->Integral() / ((sig_upper-sig_lower)/binwidth_spb) / avg_height;


	hToF_b->Scale(1./weight);
	hWp_b->Scale(1./weight);
	hW_b->Scale(1./weight);
	hToF_b->SetLineColor(2);
	hWp_b->SetLineColor(2);
	hW_b->SetLineColor(2);
	
	*/
	outFile->cd();
	for( int thisSlice = 0 ; thisSlice < x_slices ; thisSlice++ ){
		PrettyHist( ToF_spb_xslices[thisSlice], true, 0  ,false,thisSlice);
		ToF_spb_xslices[thisSlice] -> Write();
	}
	hToF_spb->Write();
	hToF_b->Write();
	hToF_orig_spb->Write();
	hToF_orig_b->Write();
	hWp_b->Write();
	hW_b->Write();
	hW_spb->Write();
	hWp_spb->Write();
	hW_full->Write();
	p_n_loX->Write();
	hXp_basic->Write();
	

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

void PrettyHist( TH1D * hist , bool sig , double weight , bool xp, int thisSlice ){
	hist->SetLineWidth(3);
	hist->SetMinimum(0);
	if( xp ) hist->SetTitle(Form("ToF for |x'-%.2f|<0.05, W' > 2",0.05+thisSlice*0.1 ));
	else{ hist->SetTitle(Form("ToF for |x-%.2f|<0.05, W > 2",0.05+thisSlice*0.1 )); }
	if( !sig){ hist->SetLineColor(2); hist->Scale( weight ); }
	return;
}
