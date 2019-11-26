#include "prop_tree.h"

Prop_Particle::Prop_Particle()
{}

Prop_Particle::~Prop_Particle()
{}

ClassImp(Prop_Particle);

Prop_Event::Prop_Event()
{
  particles.clear();
}

Prop_Event::~Prop_Event()
{}

ClassImp(Prop_Event);
