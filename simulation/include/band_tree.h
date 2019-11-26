#ifndef __BAND_TREE_H__
#define __BAND_TREE_H__

#include <vector>
#include <string>
#include "TObject.h"
#include "TVector3.h"

class BAND_Hits : public TObject
{
	public:
		BAND_Hits();
		~BAND_Hits();
		TVector3 pos; // in mm
		TVector3 mom; // in MeV/c -- use initial momentum or reconstructed momentum
		double time; // in ns
		double E_dep; // in MeV
		int barNo; // indicating what bar
		std::string type;
		

		//ClassDef(BAND_Hits,11);
};

class BAND_Event : public TObject
{
	public:
		BAND_Event();
		~BAND_Event();
		std::vector<BAND_Hits> hits;

		//ClassDef(BAND_Event,6);
};

#endif
