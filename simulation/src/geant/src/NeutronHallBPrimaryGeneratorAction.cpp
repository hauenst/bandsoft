// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBPrimaryGeneratorAction.cc
/// \brief Implementation of the NeutronHallBPrimaryGeneratorAction class

#include "NeutronHallBPrimaryGeneratorAction.h"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

#include "TFile.h"
#include "TTree.h"

#include "gen_tree.h"


	NeutronHallBPrimaryGeneratorAction::NeutronHallBPrimaryGeneratorAction(void * tp)
: G4VUserPrimaryGeneratorAction()
{
	// Initialize
	eventNum = 0;
	thisEvent = NULL;

	// Set up the tree and event branch
	genTree = (TTree*)tp;
	genTree->SetBranchAddress("event",&thisEvent);
}

NeutronHallBPrimaryGeneratorAction::~NeutronHallBPrimaryGeneratorAction()
{
}

void NeutronHallBPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	// This function is called at the begining of event

	// Load the event from the tree
	genTree->GetEvent(eventNum);

	// Loop over the particles
	for (unsigned int p=0 ; p<thisEvent->particles.size() ; p++)
	{
		// Only shoot neutrons in Geant
		G4String particleName = thisEvent->particles[p].type;
		G4String wantedParticle = "neutron";
		if (particleName != wantedParticle) continue;

		// Set vertex position
		G4PrimaryVertex* thisVertex=new G4PrimaryVertex(G4ThreeVector(0.,0.,0.),thisEvent->particles[p].t0*ns);
		//G4PrimaryVertex* thisVertex=new G4PrimaryVertex(G4ThreeVector(0.,0.,0.),0.*ns);
		// Create a new primary for the event
		//std::cout << "T0 for event: " << thisEvent->particles[p].t0 << "\n";
		G4PrimaryParticle *thisParticle=new G4PrimaryParticle(G4ParticleTable::GetParticleTable()->FindParticle(G4String(thisEvent->particles[p].type)));
		thisParticle->SetMomentum(thisEvent->particles[p].momentum.X() * GeV,
				thisEvent->particles[p].momentum.Y() * GeV,
				thisEvent->particles[p].momentum.Z() * GeV);
		thisVertex->SetPrimary(thisParticle);

		anEvent->AddPrimaryVertex(thisVertex);
	}

	// Update the event number
	eventNum += 1;
	if (eventNum >= genTree->GetEntries() ){
		eventNum = 0; // If we run off the end of the tree, reset to the beginning.
	}
}

G4int NeutronHallBPrimaryGeneratorAction::getCurrentEventNum()
{
	return eventNum;
}
