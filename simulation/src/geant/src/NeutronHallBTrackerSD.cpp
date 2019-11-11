// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBTrackerSD.cc
/// \brief Implementation of the NeutronHallBTrackerSD class

#include "NeutronHallBTrackerSD.h"
#include "NeutronHallBAnalysis.h"
#include "NeutronHallBPrimaryGeneratorAction.h"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

#include "TTree.h"


using namespace CLHEP;

	NeutronHallBTrackerSD::NeutronHallBTrackerSD(void * treePtr, const G4String& name)
: G4VSensitiveDetector(name)
	//fEventAction(eventAction)
{
	//fEventAction(eventAction)
	hitList = new BAND_Event;
	outTree = (TTree*) treePtr;
	outTree->Branch(name.data(),&hitList);
}

NeutronHallBTrackerSD::~NeutronHallBTrackerSD() 
{}

void NeutronHallBTrackerSD::Initialize(G4HCofThisEvent* hce)
{
	// Clear the last event
	hitList->hits.clear();
	barHits.clear();

}

G4bool NeutronHallBTrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{ 
	// Get the track:
	G4Track * thisTrack = aStep->GetTrack();
	hitTime = thisTrack->GetGlobalTime()/ns;
	
	G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable()); 
	hitBarNo = touchable->GetReplicaNumber(0);

	hitE = aStep->GetTotalEnergyDeposit()/MeV;
	//if (hitE<= 0.) return true;

	/*
	   G4String hitParticleName = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
	   G4String wantedParticle = "neutron";
	   if (hitParticleName != wantedParticle) return true;
	   if( std::abs(aStep->GetPreStepPoint()->GetPosition().z()/mm +2620.) > 0.001) return true;
	   if (hitE != 0.) return true;
	   */

	//hitTime = aStep->GetPreStepPoint()->GetGlobalTime()/ns;
	//std::cout << "\tGlobal step time " << hitTime << "\n";
	hitPos = aStep->GetPreStepPoint()->GetPosition();
	G4ThreeVector hitMom = aStep->GetPreStepPoint()->GetMomentum();

	// For each event with multiple hits, if there are multiple hits
	// IN THE SAME BAR, take only the earliest time / position in the bar.
	// But, for each event, with multiple bars fired, save all bars fired.
	if(barHits.count(hitBarNo) == 0){
		BAND_Hits newHit;

		newHit.E_dep = hitE;
		newHit.time = hitTime; //* hitE;
		newHit.pos = TVector3(hitPos.x() / mm,hitPos.y() / mm,hitPos.z() / mm); //* hitE;
		newHit.barNo = hitBarNo;
		newHit.mom = TVector3( hitMom.x() / GeV, hitMom.y() / GeV, hitMom.z() / GeV );
		//std::cout << "\t\tMomentum: " << newHit.mom.Mag() << "\n";
		//newHit.mom = aStep->GetPreStepPoint()->GetMomentum();

		barHits[hitBarNo] = newHit;
	}
	else{
		//return true;
		barHits[hitBarNo].E_dep += hitE;
		if (hitTime < barHits[hitBarNo].time){
			barHits[hitBarNo].time = hitTime;
			barHits[hitBarNo].pos = TVector3(hitPos.x() / mm,hitPos.y() / mm,hitPos.z() / mm); //* hitE;
			barHits[hitBarNo].mom = TVector3( hitMom.x() / GeV, hitMom.y() / GeV, hitMom.z() / GeV );
			//std::cout << "\t\tMomentum: " << barHits[hitBarNo].mom.Mag() << "\n";
		}
	}

	return true;
}

void NeutronHallBTrackerSD::EndOfEvent(G4HCofThisEvent*)
{

	for(std::map<int,BAND_Hits>::iterator it=barHits.begin() ; it != barHits.end() ; it++){
		hitList->hits.push_back(it->second);

	}


}
