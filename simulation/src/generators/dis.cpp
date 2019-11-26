#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TRandom3.h"
#include "TVectorT.h"

#include "constants.h"
#include "deuteronwf.h"
#include "crossdis.h"

#include "gen_tree.h"

using namespace std;

// Settings for cross section: 
double sigmainput=40.;          //sigma parameter in rescattering amplitude [mb], default value
double betainput=8.;            //beta parameter in rescattering amplitude [GeV^-2], default value
double epsinput=-0.5;           //epsilon parameter in rescattering amplitude [], default value
double lambdainput=1.2;
double betaoffinput=8.;
int offshellset=3;              // 0=mass diff suppression, 1=dipole suppression, 2=beta parameterization, 3=no offshell, 4=full off-shell
int symm=0;                     //symmetric FSI in inclusive reaction
int phiavg=0;                   //average cross section over phi
int F_param=0;                  //structure functions parametrization 0=SLAC, 1=Christy&Bosted, 2=Alekhin et al leading twist

extern "C"{
extern struct{
	char homedir[128];
	char a09file1[128];
	char a09file2[128];
} dir_;
}//shared fortran variables

// Foam ranges
const double min_theta_e = 8.*M_PI/180.;
const double max_theta_e = 40.*M_PI/180.;
const double min_p_e = 2.;
const double max_p_e = Ebeam-1;

const double min_theta_r =155.*M_PI/180.;
const double max_theta_r =180.*M_PI/180.;
const double min_p_r =0.200;
const double max_p_r =0.650;

const double min_phi_r=0.;
const double max_phi_r=2.*M_PI;

// Foam Integrand
class disCS : public TFoamIntegrand
{
	public:
		double Density(int nDim, double * args);
};


int main(int argc, char *argv[])
{
	// Read in arguments
	if (argc != 3)
	{
		cerr << "Wrong number of arguments. Instead use: [nEvents] [/path/to/output/file]\n";
		exit(-1);
	}

	const int nEvents=atoi(argv[1]);
	TFile * outputFile = new TFile(argv[2],"RECREATE");

	// Set external variables
	if(offshellset==0) epsinput=-0.5;
	if(offshellset==1) lambdainput=0;
	if(offshellset==2) betaoffinput=16.8;
	if(offshellset==3) epsinput=-0.5;
	if(offshellset==4) epsinput=-0.5;

	strcpy(dir_.homedir,getenv("HOME"));
	strcat(dir_.a09file1,"./grids/a09.sfs_lNNC");
	strcat(dir_.a09file2,"./grids/a09.dsfs_lNNC");

	// Create a tree
	TTree * outputTree = new TTree("MCout","Generator Output");

	// Initialize the branches
	Gen_Event * thisEvent = new Gen_Event;
	outputTree->Branch("event",&thisEvent);

	// Random number generator
	TRandom3 * rand = new TRandom3(0);

	// Initialize the foam
	TFoam * csFoam = new TFoam("csFoam");
	disCS * csTotal = new disCS();
	csFoam->SetkDim(5);
	csFoam->SetRho(csTotal);
	csFoam->SetPseRan(rand);
	// optional
	csFoam->SetnCells(10000);
	csFoam->SetnSampl(10000);
	// initialize
	csFoam->Initialize();

	// Create memory for each event
	double * eventData = new double[5];

	for (int i=0 ; i<nEvents ; i++)
	{
		if (nEvents >= 10000)
			if (i%10000==0)
				cerr << "Working on event " << i << "\n";

		thisEvent->particles.clear();
		csFoam->MakeEvent();
		csFoam->GetMCvect(eventData);

		double theta_e 	= min_theta_e 	+ eventData[0]*(max_theta_e - min_theta_e);
		double p_e	= min_p_e 	+ eventData[1]*(max_p_e - min_p_e);

		double theta_r 	= min_theta_r 	+ eventData[2]*(max_theta_r - min_theta_r);
		double p_r 	= min_p_r 	+ eventData[3]*(max_p_r - min_p_r);
		double phi_r 	= min_phi_r 	+ eventData[4]*(max_phi_r - min_phi_r);

		// Pull a random phi_e and rotate all vectors by that phi_e to do the global phi rotation that
		// we didn't waste our foam dimension on
		double phi_e = 2.*M_PI * rand->Rndm();
		phi_r += phi_e;
		if( phi_r > 2.*M_PI ) phi_r -= 2.*M_PI;

		// Store electron
		Gen_Particle electron;
		electron.type="e-";
		electron.momentum.SetMagThetaPhi(p_e,theta_e,phi_e);
		electron.t0=0.;
		thisEvent->particles.push_back(electron);
		
		// Store neutron
		Gen_Particle neutron;
		neutron.type="neutron";
		neutron.momentum.SetMagThetaPhi(p_r,theta_r,phi_r);
		neutron.t0=0.;
		thisEvent->particles.push_back(neutron);

		// Fill tree
		outputTree->Fill();
	}

	// Write tree
	outputTree->Write();

	// Store the total cross section (and error) in the output file in a TVectorT
	TVectorT<double> totalCS(2);
	csFoam->GetIntegMC(totalCS[0],totalCS[1]);
	totalCS.Write("totalCS");

	outputFile->Close();

	// Clean up
	delete csFoam;
	delete csTotal;
	delete rand;

	return 0;
}

double disCS::Density(int nDim, double *args){
	// Parameters
	int proton = 1; // (0 for DIS on neutron with spectator proton, 1 for DIS on proton with spectator neutron).
	int which_wave =1;     // Non relativistic deuteron wf: 0 - Paris wf, 1 - AV18, 2 - CD Bonn, 3 - AV18*
	int decay=0;
	int num_res=1;
	
	///////////////////////////////////////////////////////////////////////
	// Electron Drawn Variables (there should be 2 of them)
	double theta_e 	= min_theta_e 	+ args[0]*(max_theta_e - min_theta_e);
	double p_e	= min_p_e 	+ args[1]*(max_p_e - min_p_e);

	// Define beam and scattered electron vector
	TVector3 beamVec(0,0,Ebeam);
	TVector3 eVec;	eVec.SetMagThetaPhi(p_e,theta_e,0.);

	// Define q vector
	TVector3 qVec;	qVec = beamVec - eVec;
	double q 	= qVec.Mag();
	double theta_q 	= qVec.Theta();
	double phi_q 	= qVec.Phi();
	if( phi_q < 0 ) phi_q += 2.*M_PI;
	
	// Define other kinematic variables
	double nu 	= Ebeam - sqrt( p_e*p_e + mE*mE );
	double Q2 	= q*q - nu*nu;
	double xB	= Q2 / (2.*mP*nu);

	
	///////////////////////////////////////////////////////////////////////
	// Neutron Drawn Variables (there should be 3 of them)
	double theta_r 	= min_theta_r 	+ args[2]*(max_theta_r - min_theta_r);
	double p_r 	= min_p_r 	+ args[3]*(max_p_r - min_p_r);
	double phi_r 	= min_phi_r 	+ args[4]*(max_phi_r - min_phi_r); // this is the lab frame phi of the neutron

	// Define neutron vector
	TVector3 nVec;	nVec.SetMagThetaPhi(p_r,theta_r,phi_r);
	
	// Calculate phi_nq, the angle between the scattering plane and reaction plane 
	// 	(or rotation angle of neutron around q vector)
	TVector3 norm_scatter = qVec.Cross( beamVec );
	norm_scatter 	= norm_scatter.Unit();

	TVector3 norm_reaction = qVec.Cross( nVec );
	norm_reaction 	= norm_reaction.Unit();

	double phi_rq 	= norm_scatter.Angle( norm_reaction );
	double theta_rq = nVec.Angle( qVec );
	TVector3 direction = norm_scatter.Cross(norm_reaction);
	if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
	}
	else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
		phi_rq *= (-1);
	}
	

	///////////////////////////////////////////////////////////////////////
	// Create quantities for cross section input and jacobian
	double E_r = sqrt( mN*mN + p_r*p_r );
	double W_primeSq = mD*mD - Q2 + mN*mN + 2.*mD*(nu-E_r) - 2.*nu*E_r + 2.*q*p_r*cos(theta_rq);
	if (W_primeSq < 0.) return 0.;
	double W_prime = sqrt(W_primeSq);
	sigmainput = (25.3+53*(W_prime-mN))/(Q2);
		// Wim's cross section:
	phi_rq = 0.; 	// let's turn off phi_rq dependence of CS, and always say it's 0, that way the foam has no
			// sensitivity to phi_rq
	double crosstotal1 = calc_cross(Ebeam, Q2, xB, p_r, theta_rq, phi_rq, proton, which_wave, decay, num_res, 0);	// [nb/GeV^4]
		// Jacobian for recoil vector:
	double jacobian_r = p_r*p_r*sin(theta_r);									// [GeV^2]
		// Jacobian for scattered electron:
	double jacobian_e = 2.*xB*Ebeam*p_e/nu * sin(theta_e);								// [GeV]
		// Additional factors
	double add_factors = (1./E_r) * (1./(2*M_PI) );									// [1/GeV]
		// Put them all together:
	double result = crosstotal1 * jacobian_r * jacobian_e * add_factors;	
		// and multiply by our integrating range
	double differential_e = (max_theta_e - min_theta_e) * (2.*M_PI)		      * (max_p_e - min_p_e);		// [GeV]
	double differential_p = (max_theta_r - min_theta_r) * (max_phi_r - min_phi_r) * (max_p_r - min_p_r);		// [GeV]
	result = result * differential_e * differential_p;

	if ( result < 0. || result!=result) return 0.;
	return result; 													// [nb]
}

