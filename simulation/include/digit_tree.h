#ifndef __DIGIT_TREE_H__
#define __DIGIT_TREE_H__

#include <string>
#include <vector>
#include "TObject.h"
#include "TVector3.h"

class Digit_Particle : public TObject
{
	public:
		Digit_Particle();
		~Digit_Particle();

		std::string type;
		TVector3 momTrue;
		TVector3 momRecon;

		TVector3 pos; // in mm
		double time; // in ns
		double E_dep; // in MeV
		int barNo; // indicating what bar

		ClassDef(Digit_Particle,1);
};

class BAND_Hit : public TObject
{
	public:
		BAND_Hit();
		~BAND_Hit();

		double timeL;
		double timeR;
		double adcL;
		double adcR;

		TVector3 pos;
};

class CLAS_Hit : public TObject
{
	public:
		CLAS_Hit();
		~CLAS_Hit();
	
		double time;
		double adc;
		TVector3 pos;
		TVector3 mom;
};

class Digit_Event : public TObject
{
	public:
		Digit_Event();
		~Digit_Event();
		std::vector<Digit_Particle> particles;

		std::vector<BAND_Hit> band_hits;
		std::vector<CLAS_Hit> clas_hits;

		ClassDef(Digit_Event,1);
};

#endif
