#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"

#include "constants.h"
#include "gen_tree.h"

using namespace std;

const double minTrueKE = 0;
const double minMomR = 0.2;
const double maxMomR = 0.6;
const double maxThetaR = 178.*M_PI/180.;
const double minThetaR = 150.*M_PI/180.;
const double bandZ = -267; // cm

inline double sq(double x){ return x*x; };
inline double beta(double p){ return 1./sqrt(1. + sq(mN/p)); };
double getNeutronCS(double threshold); // in nb/sr
double getNeutronKE(double r);

int main(int argc, char ** argv)
{
	if (argc != 3)
	{
		cerr << "Wrong number of arguments. Instead use:\n\t/path/to/inclusive/file /path/to/outputfile\n";
		exit(-1);
	}

	// Prep
	TRandom3 * myRand = new TRandom3(0);
	const double maxBetaR = beta(maxMomR);
	const double minBetaR = beta(minMomR);
	const double minCosThetaR = cos(maxThetaR);
	const double maxCosThetaR = cos(minThetaR);
	const double minTime = -100;//-bandZ/(maxBetaR*cAir) - 5.;  		// Put in 5 ns for cushion
	const double maxTime = 100;//(7.2*5-bandZ)/(minBetaR*cAir) + 5.; 	// Put in another 5 ns for more cushion
	const double timeWindow = maxTime - minTime;
	const double BAND_center =  bandZ - (7.2*5)/2.;

	// Open the files
	TFile *infile = new TFile(argv[1]);
	TFile *outfile = new TFile(argv[2],"RECREATE");

	// Get the cross section information
	TVectorT<double> * csVec = (TVectorT<double>*)infile->Get("totalCS");

	// Get the tree and assign branches
	TTree * inTree = (TTree*) infile->Get("MCout");
	Gen_Event * inEvent = NULL;
	inTree->SetBranchAddress("event",&inEvent);

	// Produce an output tree
	TTree * outTree = new TTree("MCout","Random Coincidence Output");
	Gen_Event * outEvent = new Gen_Event;
	outTree->Branch("event",&outEvent);
	// Loop over the events
	const int nEvents = inTree->GetEntries();
	for (int i=0 ; i<nEvents ; i++)
	{
		if (i % 10000 == 0) cout << "Event " << i << endl;
		inTree->GetEntry(i);

		// Require that there be one particle in the event
		if (inEvent->particles.size() != 1){
			cerr << "Event " << i << " is not an inclusive event! It has " << inEvent->particles.size() << "particles!!!\n";
			continue;
		}
		
		// Generate random neutron
		double cosThetaR = minCosThetaR + myRand->Rndm()*(maxCosThetaR-minCosThetaR);
		double phiR = 2.*M_PI * myRand->Rndm();
		double keR = getNeutronKE(myRand->Rndm());

		// Calculate the derived neutron info
		double thetaR = acos(cosThetaR);
		double momR = sqrt(sq(keR + mN) - mN*mN);
		double betaR = momR / (keR + mN);

		// Calculate a random time for the hit to occur
		double hitTime = minTime + timeWindow*myRand->Rndm();

		// Work back the time it would take to reach the center of BAND, store to tree
		double t0 = hitTime - ( -BAND_center )/(betaR*cAir);

		// Write tree
		outEvent->particles.clear();
		outEvent->particles.push_back(inEvent->particles[0]);
		Gen_Particle neutron;
		neutron.type="neutron";
		neutron.momentum.SetMagThetaPhi(momR,thetaR,phiR);
		neutron.t0=t0;
		outEvent->particles.push_back(neutron);
		outTree->Fill();
	}
	// Update cross section info
	double neutronCS = getNeutronCS(minTrueKE) * (2.*M_PI)*(maxCosThetaR - minCosThetaR);
	TVectorT<double> csSqVec(3);
	csSqVec[0]=neutronCS*timeWindow*(*csVec)[0]; // Units of nb^2 * ns
	csSqVec[1]=neutronCS*timeWindow*(*csVec)[1];
	csSqVec[2]=timeWindow;
	csSqVec.Write("totalCSSq");

	// Clean-up
	//gFile = outfile;
	outTree->Write();
	outfile->Close();
}

double getNeutronCS(double threshold)
{
	double pavel_rate = 967451. * exp(-29.2539 * threshold) + 5.43632e+06 * exp(-548.598 * threshold);

	return pavel_rate / (3.E36) * (1.E33) / (0.1); // divide by Pavel's lumi, convert to nb, divide by 0.1 sr
}

double CDF(double Tr)
{
	return 1. - getNeutronCS(Tr)/getNeutronCS(minTrueKE);
}

double getNeutronKE(double r)
{
	double minKE = minTrueKE;
	double maxKE = E1;
	double testKE = 0.5*(minKE + maxKE);

	while (fabs(maxKE - minKE) > 1.E-6)
	{
		double testCDF = CDF(testKE);

		if (r < testCDF)
			maxKE = testKE;
		else
			minKE = testKE;

		testKE = 0.5*(minKE + maxKE);
	}

	return testKE;
}
