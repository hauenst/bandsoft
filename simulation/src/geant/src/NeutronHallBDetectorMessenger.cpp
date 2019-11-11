// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
// $Id: NeutronHallBDetectorMessenger.cc 69706 2013-05-13 09:12:40Z gcosmo $
// 
/// \file NeutronHallBDetectorMessenger.cc
/// \brief Implementation of the NeutronHallBDetectorMessenger class

#include "NeutronHallBDetectorMessenger.h"
#include "NeutronHallBDetectorConstruction.h"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHallBDetectorMessenger::NeutronHallBDetectorMessenger(NeutronHallBDetectorConstruction* Det)
 : G4UImessenger(),
   fDetectorConstruction(Det)
{
  fNeutronHallBDirectory = new G4UIdirectory("/NeutronHallB/");
  fNeutronHallBDirectory->SetGuidance("UI commands specific to this example.");

  fDetDirectory = new G4UIdirectory("/NeutronHallB/det/");
  fDetDirectory->SetGuidance("Detector construction control");

  fTargMatCmd = new G4UIcmdWithAString("/NeutronHallB/det/setTargetMaterial",this);
  fTargMatCmd->SetGuidance("Select Material of the Target.");
  fTargMatCmd->SetParameterName("choice",false);
  fTargMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fChamMatCmd = new G4UIcmdWithAString("/NeutronHallB/det/setChamberMaterial",this);
  fChamMatCmd->SetGuidance("Select Material of the Chamber.");
  fChamMatCmd->SetParameterName("choice",false);
  fChamMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fStepMaxCmd = new G4UIcmdWithADoubleAndUnit("/NeutronHallB/det/stepMax",this);
  fStepMaxCmd->SetGuidance("Define a step max");
  fStepMaxCmd->SetParameterName("stepMax",false);
  fStepMaxCmd->SetUnitCategory("Length");
  fStepMaxCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHallBDetectorMessenger::~NeutronHallBDetectorMessenger()
{
  delete fTargMatCmd;
  delete fChamMatCmd;
  delete fStepMaxCmd;
  delete fNeutronHallBDirectory;
  delete fDetDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHallBDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == fTargMatCmd )
   { fDetectorConstruction->SetTargetMaterial(newValue);}

  if( command == fChamMatCmd )
   { fDetectorConstruction->SetChamberMaterial(newValue);}

  if( command == fStepMaxCmd ) {
    fDetectorConstruction
      ->SetMaxStep(fStepMaxCmd->GetNewDoubleValue(newValue));
  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
