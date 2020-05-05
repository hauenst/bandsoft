#ifndef EFIDUCIAL_HH
#define EFIDUCIAL_HH

class TF1;

class ElectronFiducial{

public:
	ElectronFiducial(TString, TString);
	~ElectronFiducial();

	double GetElectronAcceptance(double, double, double);

private:

	void ReadMomentumParameters(TString, TString);

	void InitMomentumFunctions();
	void InitFiducialFunctions();
	void SetFiducialParameters(double);

	double GetFiducialParameter(int, int, int, double);
	
	double momentumPar[2][6][3][3];	

	TF1* fiducialFunction[2][6];
	TF1* momentumFunction[2][6][3];
};

#endif
