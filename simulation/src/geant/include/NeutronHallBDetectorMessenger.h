// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBDetectorMessenger.hh
/// \brief Definition of the NeutronHallBDetectorMessenger class

#ifndef NeutronHallBDetectorMessenger_h
#define NeutronHallBDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class NeutronHallBDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Messenger class that defines commands for NeutronHallBDetectorConstruction.
///
/// It implements commands:
/// - /NeutronHallB/det/setTargetMaterial name
/// - /NeutronHallB/det/setChamberMaterial name
/// - /NeutronHallB/det/stepMax value unit

class NeutronHallBDetectorMessenger: public G4UImessenger
{
  public:
    NeutronHallBDetectorMessenger(NeutronHallBDetectorConstruction* );
    virtual ~NeutronHallBDetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    NeutronHallBDetectorConstruction*  fDetectorConstruction;

    G4UIdirectory*           fNeutronHallBDirectory;
    G4UIdirectory*           fDetDirectory;

    G4UIcmdWithAString*      fTargMatCmd;
    G4UIcmdWithAString*      fChamMatCmd;

    G4UIcmdWithADoubleAndUnit* fStepMaxCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
