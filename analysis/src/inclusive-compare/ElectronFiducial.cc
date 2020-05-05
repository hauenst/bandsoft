#include "ElectronFiducial.hh"

#include <fstream>

ElectronFiducial::ElectronFiducial(TString upperParFile, TString lowerParFile) {

	ReadMomentumParameters(upperParFile, lowerParFile);
	
	InitMomentumFunctions();
	InitFiducialFunctions();

}

ElectronFiducial::~ElectronFiducial() {
}

void ElectronFiducial::ReadMomentumParameters(TString upperParFile, TString lowerParFile) {

	ifstream upperFile(upperParFile);
	ifstream lowerFile(lowerParFile);
	
	int sector, fidPar;
	double momPar0, momPar1, momPar2;

	while(upperFile >> sector >> fidPar >> momPar0 >> momPar1 >> momPar2) {

		momentumPar[0][sector-1][fidPar][0] = momPar0;
		momentumPar[0][sector-1][fidPar][1] = momPar1;
		momentumPar[0][sector-1][fidPar][2] = momPar2;

	}

	while(lowerFile >> sector >> fidPar >> momPar0 >> momPar1 >> momPar2) {

		momentumPar[1][sector-1][fidPar][0] = momPar0;
		momentumPar[1][sector-1][fidPar][1] = momPar1;
		momentumPar[1][sector-1][fidPar][2] = momPar2;

	}


}

void ElectronFiducial::InitMomentumFunctions() {

	TString momentum_functions[3] = {"[0] + exp([1]*x+[2])", "[0] - exp([1]*x+[2])", "[0] + exp([1]*x+[2])"};
	TString boundary[2] = {"upper","lower"};

	for(int sector = 0; sector < 6; sector++) {
		for(int fidPar = 0; fidPar < 3; fidPar++) {
			for(int i = 0; i < 2; i++) {
				momentumFunction[i][sector][fidPar] = new TF1(Form("%sMomFunc_%i_%i", boundary[i].Data(), sector, fidPar), momentum_functions[fidPar], 0., 12.);
				for(int momPar = 0; momPar < 3; momPar++) {
					momentumFunction[i][sector][fidPar]->SetParameter(momPar, momentumPar[i][sector][fidPar][momPar]);
				}	
			}
		}
	}
}


void ElectronFiducial::InitFiducialFunctions() {

	TString fiducial_functions[2] = {"[0] - exp([1]*x+[2])", "[0] + exp([1]*x+[2])"};
	TString boundary[2] = {"upper","lower"};

	for(int sector = 0; sector < 6; sector++) {
		for(int i = 0; i < 2; i++) {
			fiducialFunction[i][sector] = new TF1(Form("%sFidFunc_%i", boundary[i].Data(), sector), fiducial_functions[i], 0., 60.);
		}
	}
}

void ElectronFiducial::SetFiducialParameters(double pe) {

	for(int sector = 0; sector < 6; sector++) {
		for(int i = 0; i < 2; i++) {
			for(int fidPar = 0; fidPar < 3; fidPar++) {
				double par = GetFiducialParameter(i, sector, fidPar, pe);
				fiducialFunction[i][sector]->SetParameter(fidPar, par);
			}
		}
	}
}


double ElectronFiducial::GetFiducialParameter(int boundary, int sector, int par, double pe) {

	return momentumFunction[boundary][sector][par]->Eval(pe);

}


double ElectronFiducial::GetElectronAcceptance(double theta, double phi, double pe) {

	SetFiducialParameters(pe);

	for(int sector = 0; sector < 6; sector++) {

		double maxPhi = fiducialFunction[0][sector]->Eval(theta);
		double minPhi = fiducialFunction[1][sector]->Eval(theta);

		if (phi > minPhi && phi < maxPhi) {
			return 1.;
		}
	}

	return 0.;
}



