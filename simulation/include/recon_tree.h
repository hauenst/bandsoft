#ifndef __RECON_TREE_H__
#define __RECON_TREE_H__

#include <string>
#include <vector>
#include "TObject.h"
#include "TVector3.h"

class Recon_Particle : public TObject
{
	public:
		Recon_Particle();
		~Recon_Particle();

		std::string type;
		double tof; // in ns
		TVector3 pos; // in mm
		double E_dep; // in MeV
		TVector3 mom;

		// old so i don't break my code:
		TVector3 momTrue;
		TVector3 momRecon;
		double time; // in ns
		int barNo; // indicating what bar

		ClassDef(Recon_Particle,3);
};

class Recon_Event : public TObject
{
	public:
		Recon_Event();
		~Recon_Event();
		std::vector<Recon_Particle> particles;

		ClassDef(Recon_Event,3);
};

#endif
