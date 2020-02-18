#include "gen_tree.h"

Gen_Particle::Gen_Particle()
{}

Gen_Particle::~Gen_Particle()
{}

ClassImp(Gen_Particle);

Gen_Event::Gen_Event()
{
  particles.clear();
}

Gen_Event::~Gen_Event()
{}

ClassImp(Gen_Event);
