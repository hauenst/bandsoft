#ifndef __PROP_TREE_H__
#define __PROP_TREE_H__

#include <vector>
#include <string>
#include "TObject.h"
#include "TVector3.h"

class Prop_Particle : public TObject
{
	public:
		Prop_Particle();
		~Prop_Particle();
		TVector3 pos; // in mm
		TVector3 mom; // in MeV/c -- use initial momentum or reconstructed momentum
		double time; // in ns
		double E_dep; // in MeV
		int barNo; // indicating what bar
		std::string type;
		

		ClassDef(Prop_Particle,1);
};

class Prop_Event : public TObject
{
	public:
		Prop_Event();
		~Prop_Event();
		std::vector<Prop_Particle> particles;

		ClassDef(Prop_Event,1);
};

#endif
