#include <iostream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "TFile.h"
#include "TTree.h"
#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TRandom3.h"
#include "TVectorT.h"
#include "TF1.h"

#include "constants.h"
#include "spWf.h"
#include "sigmaep_bound.h"

using namespace std;

// Functions to be used
TVector3 makeVector( double mag, double theta, double phi );
double getAngleBetween( double theta1, double phi1, double theta2, double phi2);
double protonMag( double omega, double q2, double cosTheta_pq, double choice);

TF1 * momWeight = new TF1("momWeight", "(x<=0.02)*1 + (x>0.02 && x<=0.6)*exp( -11*(x-0.1) ) + (x>0.6)*1");
TF1 * thetaWeight = new TF1("thetaWeight","exp( -(x-1.4)*(x-1.4)/(2*0.5*0.5) )");

// Foam Integrand class for d(e,e'p) QE cross section
class deepCS : public TFoamIntegrand{
	public:
		double Density(int nDim, double * args);
};

// Ranges for generator
const double Ebeam 		= 2.2;
const double min_theta_e 	= 9.5	*M_PI/180.;	// 5 deg is out-bending min, 9.5 deg is in-bending min
const double max_theta_e 	= 40.	*M_PI/180.;
const double min_phi_e 		= 0.	*M_PI/180.;
const double max_phi_e 		= 360.	*M_PI/180.;
const double min_p_e 		= 1.		  ;
const double max_p_e 		= Ebeam		  ;
const double min_theta_p 	= 5	*M_PI/180.;
const double max_theta_p 	= 120.	*M_PI/180.;
const double min_phi_p		= 0.	*M_PI/180.;
const double max_phi_p 		= 360.	*M_PI/180.;

// Initialize wavefunction class
spWf * momDist = new spWf();
// Initialize CC1 class
sigmaep_bound * epCC1 = new sigmaep_bound();

// Main
int main(int argc, char ** argv){

	if (argc<3){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./gen_eep [nEvents] [outputFile]\n";
		return -1;
	}

	// Set up files
	TFile * outFile = new TFile(argv[2],"RECREATE");
	TTree * outTree = new TTree("MCout","Generator Output");
	
	// Set output branches
	double mom_pE,mom_pP,mom_pN,Ebeam_out;
	double theta_pE,theta_pP,theta_pN;
	double phi_pE,phi_pP,phi_pN;
	double xB, weight;
	outTree->Branch("mom_pE", &mom_pE);
	outTree->Branch("mom_pP", &mom_pP);
	outTree->Branch("mom_pN", &mom_pN);
	outTree->Branch("Ebeam",  &Ebeam_out);
	outTree->Branch("theta_pE", &theta_pE);
	outTree->Branch("theta_pP", &theta_pP);
	outTree->Branch("theta_pN", &theta_pN);
	outTree->Branch("phi_pE", &phi_pE);
	outTree->Branch("phi_pP", &phi_pP);
	outTree->Branch("phi_pN", &phi_pN);
	outTree->Branch("xB",&xB);
	outTree->Branch("weight",&weight);

	// Initialize random number seed
	TRandom3 * rand = new TRandom3(0);

	// Create TFoam 
	TFoam * csFoam = new TFoam("csFoam");

	// Initialize our cross section for d(e,e'p)n QE
	deepCS * csTotal = new deepCS();

	// Set up Foam
	csFoam->SetkDim(6);
	csFoam->SetRho(csTotal);
	csFoam->SetPseRan(rand);
	// optional
	csFoam->SetnCells(10000);
	csFoam->SetnSampl(1000);
	// initialize -- this will run the cells specified above, see DENSITY function later in this code
	csFoam->Initialize();
	
	// Create memory for each event
	double * eventData = new double[6];
	const int nEvents=atoi(argv[1]);
	for (int i=0 ; i<nEvents ; i++){
		if ( (i%10000) ==1) cerr << "Working on event " << i << "\n";

		// Grab an event from our CS foam
		csFoam->MakeEvent();
		csFoam->GetMCvect(eventData);

		double theta_e 	= min_theta_e + eventData[0]*(max_theta_e - min_theta_e);
		double phi_e	= min_phi_e + eventData[1]*(max_phi_e - min_phi_e);
		double p_e 	= min_p_e + eventData[2]*(max_p_e - min_p_e);
		
		double theta_p 	= min_theta_p + eventData[3]*(max_theta_p - min_theta_p);
		double phi_p 	= min_phi_p + eventData[4]*(max_phi_p - min_phi_p);

		double choice = eventData[5];

		// Create e' and e vectors
		TVector3 ePrime = makeVector( p_e, theta_e, phi_e );
		TVector3 e0 = makeVector( Ebeam, 0 , 0);
		// Create Momentum transfer
		TVector3 q = e0 - ePrime;
		double omega = Ebeam - p_e;
		double phi_q = phi_e - 180. * M_PI/180.;
		double theta_q = acos( q.Z() / q.Mag() );

		// Get angle between proton and q
		double cosTheta_pq = getAngleBetween( theta_q, phi_q, theta_p, phi_p);
			
		// Now solve for proton momentum given kinematics
		double p_p = protonMag( omega, q.Mag2() , cosTheta_pq, choice);

		// Now define pPrime vector and corresponding neutron vector
		TVector3 pPrime = makeVector( p_p, theta_p, phi_p );
		TVector3 pN_init = q - pPrime;
		TVector3 pP_init = pPrime - q;

		// Values to store in the generator tree
		mom_pE = p_e;
		mom_pP = p_p;
		mom_pN = pN_init.Mag();
		Ebeam_out = Ebeam;
		theta_pE = theta_e;
		theta_pP = theta_p;
		theta_pN = pN_init.Theta();
		phi_pE = phi_e;
		phi_pP = phi_p;
		phi_pN = pN_init.Phi();

		double Q2 = -( pow(omega,2) - q.Mag2() ); 	// Q2
		xB = Q2 / (2. * omega * mP);		// Bjorken x

		weight = (1./momWeight->Eval(mom_pN))*(1./thetaWeight->Eval(theta_pN));
		//weight = 1;
		outTree->Fill();
	}

	outFile->cd();
	outTree->Write();

	// Store CS from this file
	TVectorT<double> totalCS(2);
	csFoam->GetIntegMC(totalCS[0],totalCS[1]);
	totalCS.Write("totalCS");

	// Clean up
	outFile->Close();
	delete csFoam;
	delete csTotal;
	delete rand;

	return 0;
}

TVector3 makeVector( double mag, double theta, double phi ){
	TVector3 vec ( mag*sin(theta)*cos(phi) , mag*sin(theta)*sin(phi) , mag*cos(theta) );
	return vec;
}
double getAngleBetween( double theta1, double phi1, double theta2, double phi2){
	// returns cos(theta12) assuming that the norm of both vectors is 1
	return sin(theta1)*sin(theta2) * ( cos(phi1)*cos(phi2) + sin(phi1)*sin(phi2) ) + cos(theta1)*cos(theta2) ;
}

double protonMag( double omega, double q2, double cosTheta_pq, double choice){
	// Solve for proton momentum given kinematics
	
	double x = mD + omega;
	double y2 = x*x + mP*mP - mN*mN - q2;
	
	double A = 4. * q2 * pow(cosTheta_pq,2)  -  4.*x*x ;
	double B = 4. * y2 * sqrt(q2) * cosTheta_pq;
	double C = y2*y2 - 4. * x*x * mP*mP;

	// Make sure determinant is positive
	double p_p = 0.;
	if( (B*B - 4.*A*C) < 0 ){ return 0.; }
	else{
		double p_p1 = ( -B + sqrt(B*B-4.*A*C) ) / (2.*A);
		double p_p2 = ( -B - sqrt(B*B-4.*A*C) ) / (2.*A);

		// Now check if they are valid:
		bool validP1 = (p_p1 >= 0.) && ( (x - sqrt(p_p1*p_p1 + mP*mP) ) >= 0.) && ( (y2 + 2*p_p1*sqrt(q2) * cosTheta_pq ) >= 0.) ;
		bool validP2 = (p_p2 >= 0.) && ( (x - sqrt(p_p2*p_p2 + mP*mP) ) >= 0.) && ( (y2 + 2*p_p2*sqrt(q2) * cosTheta_pq ) >= 0.) ;

		if( (!validP1) && (!validP2) ) return 0.;
		if( (!validP1) ) p_p = p_p2;
		else if( (!validP2) ) p_p = p_p1;
		else{ // both valid, so choose one
			if( choice > 0.5 ) p_p = p_p1;
			else if( choice <= 0.5) p_p = p_p2;
		}
 
	}
	return p_p;
}


double deepCS::Density(int nDim, double *args){
	// Returns CS sample -- used to generate our foam surface
	
	// Variables (there should be 6 of them)
	double theta_e 	= min_theta_e + args[0]*(max_theta_e - min_theta_e);
	double phi_e	= min_phi_e + args[1]*(max_phi_e - min_phi_e);
	double p_e 	= min_p_e + args[2]*(max_p_e - min_p_e);
	
	double theta_p 	= min_theta_p + args[3]*(max_theta_p - min_theta_p);
	double phi_p 	= min_phi_p + args[4]*(max_phi_p - min_phi_p);

	double choice = args[5];

	// Create e' and e vectors
	TVector3 ePrime = makeVector( p_e, theta_e, phi_e );
	TVector3 e0 = makeVector( Ebeam, 0 , 0);
	// Create Momentum transfer
	TVector3 q = e0 - ePrime;
	double omega = Ebeam - p_e;
	double phi_q = phi_e - 180. * M_PI/180.;
	double theta_q = acos( q.Z() / q.Mag() );
	// Get angle between proton and q
	double cosTheta_pq = getAngleBetween( theta_q, phi_q, theta_p, phi_p);
	
	// Now solve for proton momentum given kinematics
	double p_p = protonMag( omega, q.Mag2() , cosTheta_pq, choice);

	// Now define pPrime vector and corresponding neutron vector
	TVector3 pPrime = makeVector( p_p, theta_p, phi_p );
	TVector3 pN_init = q - pPrime;
	TVector3 pP_init = pPrime - q;

	double En = sqrt( pN_init.Mag2() + mN*mN );
	double Ep = sqrt( pPrime.Mag2()  + mP*mP );

	double jacobian = En*Ep / fabs( En*p_p + Ep*(p_p - q.Mag()*cosTheta_pq) );
	double epCS = epCC1->sigmaCC1(Ebeam, ePrime, pPrime);

	double differential_e = (max_theta_e - min_theta_e)*(2.*M_PI)*sin(theta_e)*(max_p_e - min_p_e);
	double differential_p = (max_theta_p - min_theta_p)*(2.*M_PI)*sin(theta_p);

	double result = -1.;
	result = (epCS) * (p_p*p_p) * (jacobian) * (momDist->getDensity(pP_init.Mag())) * (differential_e) * (differential_p); 
	result *= 1./(momWeight->Eval(pP_init.Mag()));
	result *= 1./(thetaWeight->Eval(pN_init.Theta()) );
		// in [cm2] * [GeV2] * [1] * [1/GeV3] * [GeV] * [1]
	
	if (result < 0.) return 0.;
	return result;
}
