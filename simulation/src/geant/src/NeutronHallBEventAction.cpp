// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBEventAction.cc
/// \brief Implementation of the NeutronHallBEventAction class

#include "NeutronHallBEventAction.h"
#include "NeutronHallBRunAction.h"
#include "NeutronHallBAnalysis.h"
#include "NeutronHallBTrackerHit.h"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"


#include <G4VHit.hh>
#include <G4VHitsCollection.hh>
#include "Randomize.hh"

#include "G4UnitsTable.hh"
#include <iomanip>
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"

#include "TTree.h"

	NeutronHallBEventAction::NeutronHallBEventAction(void * treePtr)
: G4UserEventAction()
	//fEdep(0.)
{
	rand2  = new TRandom2();
	outTree = (TTree*) treePtr;
}

NeutronHallBEventAction::~NeutronHallBEventAction()
{}

void NeutronHallBEventAction::BeginOfEventAction(const G4Event*){
	//fEdep = 0.;
	//    G4cout << "starting event...." << G4endl;
}

void NeutronHallBEventAction::EndOfEventAction(const G4Event* event)
{
	//fRunAction->AddEdep(fEdep);
	G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
	G4int n_trajectories = 0;
	if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
	G4int eventID = event->GetEventID();
	//G4cout << ">>> Event: " << eventID  << G4endl;
	if ( eventID % 10000  == 0) {
		G4cout << ">>> Event: " << eventID  << G4endl;
		if ( trajectoryContainer ) {
			G4cout << "    " << n_trajectories
				<< " trajectories stored in this event." << G4endl;
		}
	}

	outTree->Fill();

}

