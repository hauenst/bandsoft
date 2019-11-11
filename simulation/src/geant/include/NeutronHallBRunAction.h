// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBRunAction.hh
/// \brief Definition of the NeutronHallBRunAction class

#ifndef NeutronHallBRunAction_h
#define NeutronHallBRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Event.hh"
#include "globals.hh"
//#include <TFile.h>
//#include <TTree.h>
//#include <TLorentzVector.h>
//#include <TVector3.h>
//#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;

/// Run action class

class NeutronHallBRunAction : public G4UserRunAction
{
private:

public:
    
    NeutronHallBRunAction();
    virtual ~NeutronHallBRunAction();

    virtual void BeginOfRunAction(const G4Run* run);
    virtual void EndOfRunAction(const G4Run* run);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
