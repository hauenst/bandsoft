#include "digit_tree.h"

// For digit particle:
Digit_Particle::Digit_Particle()
{}

Digit_Particle::~Digit_Particle()
{}

ClassImp(Digit_Particle);

// For BAND hit:
BAND_Hit::BAND_Hit()
{}
BAND_Hit::~BAND_Hit()
{}

// For CLAS hit:
CLAS_Hit::CLAS_Hit()
{}
CLAS_Hit::~CLAS_Hit()
{}

// For digit event:
Digit_Event::Digit_Event()
{
	particles.clear();
	band_hits.clear();
	clas_hits.clear();
}

Digit_Event::~Digit_Event()
{}

ClassImp(Digit_Event);
