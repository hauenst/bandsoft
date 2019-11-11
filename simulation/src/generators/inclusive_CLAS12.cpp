#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "TFile.h"
#include "TTree.h"
#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TRandom3.h"
#include "TVectorT.h"

using namespace std;

#include "constants.h"
#include "deuteronwf.h"
#include "crossdis.h"
#include "crossincl.h"
#include "gen_tree.h"

// Settings for cross section: 
double sigmainput=40.;          //sigma parameter in rescattering amplitude [mb], default value
double betainput=8.;            //beta parameter in rescattering amplitude [GeV^-2], default value
double epsinput=-0.5;           //epsilon parameter in rescattering amplitude [], default value
double lambdainput=1.2;
double betaoffinput=8.;
int offshellset=3;              // 0=mass diff suppression, 1=dipole suppression, 2=beta parameterization, 3=no offshell, 4=full off-shell
int symm=0;                     //symmetric FSI in inclusive reaction
int phiavg=1;                   //average cross section over phi
int F_param=0;                  //structure functions parametrization 0=SLAC, 1=Christy&Bosted, 2=Alekhin et al leading twist

// Constant settings
const int calc=0;
const int which_wave=1;
const int this_decay=0;
const int num_res=1;

// Memory
double *c, *d, *m;
int c_length=0;

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

// Foam Integrand
class InclusiveCS : public TFoamIntegrand
{
	public:
		double Density(int nDim, double * args);
};

int main(int argc, char *argv[]){
	// Read in arguments
	if (argc != 3){
		cerr << "Wrong number of arguments. Instead use: [nEvents] [/path/to/output/file]\n";
		return -1;
	}
	const int nEvents=atoi(argv[1]);
	TFile * outputFile = new TFile(argv[2],"RECREATE");

	// Run sanity check on parameters
	sanitycheck(calc,0,F_param,which_wave,offshellset);

	// Get some parameters
	get_wf_param(&c, &d, &m, c_length, which_wave);

	// Read in data files
	strcpy(dir_.homedir,getenv("HOME"));
	strcpy(dir_.a09file1,getenv("HOME"));
	strcpy(dir_.a09file2,getenv("HOME"));
	strcat(dir_.a09file1,"/software/band/deuteron_dis/wim_lib/grids/a09.sfs_lNNC");
	strcat(dir_.a09file2,"/software/band/deuteron_dis/wim_lib/grids/a09.dsfs_lNNC");

	// Create a tree
	TTree * outputTree = new TTree("MCout","Generator Output");

	// Initialize the branches
	Gen_Event * thisEvent = new Gen_Event;
	outputTree->Branch("event",&thisEvent);

	// Random number generator
	TRandom3 * rand = new TRandom3(0);

	// Initialize the foam
	TFoam * csFoam = new TFoam("csFoam");
	InclusiveCS * csTotal = new InclusiveCS();
	csFoam->SetkDim(2);
	csFoam->SetRho(csTotal);
	csFoam->SetPseRan(rand);
	// optional
	csFoam->SetnCells(10000);
	csFoam->SetnSampl(10000);
	// initialize
	csFoam->Initialize();

	// Create memory for each event
	double * eventData = new double[2];

	for (int i=0 ; i<nEvents ; i++){
		if (nEvents >= 10000)
			if (i%10000==0)
				cerr << "Working on event " << i << "\n";

		csFoam->MakeEvent();
		csFoam->GetMCvect(eventData);

		// Extract useful quantities
		double theta_e = min_theta_e + eventData[0]*(max_theta_e - min_theta_e);
		double p_e = min_p_e + eventData[1]*(max_p_e - min_p_e);

		// Handle phi
		double phi_e = 2.*M_PI * rand->Rndm();

		// Write to tree
		thisEvent->particles.clear();
		Gen_Particle electron;
		electron.type="e-";
		electron.momentum.SetMagThetaPhi(p_e,theta_e,phi_e);
		electron.t0=0.;
		thisEvent->particles.push_back(electron);
		outputTree->Fill();
	}

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

inline double sq(double x) {return x*x;};

double InclusiveCS::Density(int nDim, double *args)
{
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

	// Some memory for the cross section calculation results
	double QEpw, DISpw, DISfsi, DISfsi2;

	// Calculate the cross section
	calc_inclusive2(QEpw, DISpw, DISfsi, DISfsi2, Ebeam, Q2,xB, c, d, m, c_length, which_wave, offshellset, this_decay, num_res, calc);
	double dsigma=DISpw+QEpw-DISfsi;
	double jacobian = xB*Ebeam*p_e/(M_PI*nu);
	double differential = (max_theta_e - min_theta_e)*(2.*M_PI)*(max_p_e - min_p_e)*sin(theta_e);

	return dsigma*jacobian*differential;
}
