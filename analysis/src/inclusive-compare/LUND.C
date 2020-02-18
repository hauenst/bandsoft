#include "gen_tree.cpp"

void LUND(TString filename) {

	TFile* file = new TFile(Form("%s.root", filename.Data()));
	TTree* MCout=(TTree*)file->Get("MCout");

	double targetL = 5.;	// cm
	double rasterX = 0.04;	// cm
	double rasterY = 0.04; 	// cm

	const int nParticles = 1;
	TString particleType[nParticles] = {"e-"};
        int particleID[nParticles] = {11};
        double particleMass[nParticles] = {0.511e-3};   // GeV

	int targA = 2;
	int targZ = 1;
	double targP = 0.;
	double beamP = 0.;
	int interactN = 1;
	int beamType = 11;
	double beamE = 10.6;	// GeV
	double weight = 1.0;

	unsigned int seed = time(0) + (int) getpid();
	TRandom3 *r = new TRandom3(seed);

	Gen_Event* event = new Gen_Event();
	Gen_Particle* particle = new Gen_Particle(); 

	MCout->SetBranchAddress("event", &event);
	
	int nEvents = MCout->GetEntries();

	ofstream outfile;
//	TString outfilename = Form("lund_%s.dat", filename.Data());
	TString outfilename = "lund_wim_inclusive.dat";
	outfile.open(outfilename); 

	TString formatstring, outstring;

	for (int i = 0; i<nEvents; i++) {
		if( i == 10000 ) break;

		MCout->GetEntry(i);

		// LUND header for the event:
		formatstring = "%i \t %i \t %i \t %.3f \t %.3f \t %i \t %.1f \t %i \t %i \t %.3f \n";
		outstring = Form(formatstring, nParticles, targA, targZ, targP, beamP, beamType, beamE, interactN, i, weight);

		outfile << outstring; 

		double vx = r->Uniform(-rasterX/2., rasterX/2.);	
		double vy = r->Uniform(-rasterY/2., rasterY/2.);
		double vz = r->Uniform(-targetL/2., targetL/2.);	

		// LUND info for each particle in the event
		for (int j = 0; j<nParticles; j++) {

			particle = &event->particles[j];
			
			double px = particle->momentum.x();
			double py = particle->momentum.y();
			double pz = particle->momentum.z();

			double E = sqrt(px*px + py*py + pz*pz + particleMass[j]*particleMass[j] );			

			formatstring = "%i \t %.3f \t %i \t %i \t %i \t %i \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \n";
			outstring = Form(formatstring, j+1, 0.0, 1, particleID[0], 0, 0, px, py, pz, E, particleMass[0], vx, vy, vz);  
		
			outfile << outstring;
		
		}	


	}

	outfile.close();

}
