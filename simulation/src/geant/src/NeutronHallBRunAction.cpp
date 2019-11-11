// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
// $Id: NeutronHallBRunAction.cc 75214 2013-10-29 16:04:42Z gcosmo $
//
/// \file NeutronHallBRunAction.cc
/// \brief Implementation of the NeutronHallBRunAction class

#include "NeutronHallBRunAction.h"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "NeutronHallBAnalysis.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHallBRunAction::NeutronHallBRunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(1000);
    
    // Create Ntuple
    G4AnalysisManager* man = G4AnalysisManager::Instance(); // Open an output file
    man->CreateNtuple("Ntuple", "generated data");
    man->CreateNtupleIColumn("Particle_ID_g");        //branch 0
    man->CreateNtupleDColumn("KineticEnergy_g");      //branch 1
    man->CreateNtupleDColumn("Momentum_g");           //branch 2
    man->CreateNtupleDColumn("x_g");                  //branch 3
    man->CreateNtupleDColumn("y_g");                  //branch 4
    man->CreateNtupleDColumn("theta_g");              //branch 5
    man->CreateNtupleIColumn("Sector_g");             //branch 6
    
    man->CreateNtupleIColumn("NeutronReachedDet");    //branch 7
    man->CreateNtupleIColumn("Particle_ID");          //branch 8
    man->CreateNtupleDColumn("KineticEnergy");        //branch 9
    man->CreateNtupleDColumn("Momentum");             //branch 10
    man->CreateNtupleDColumn("Theta");                //branch 11
    man->CreateNtupleDColumn("x");                    //branch 12
    man->CreateNtupleDColumn("y");                    //branch 13
    man->CreateNtupleDColumn("z");                    //branch 14
    man->CreateNtupleIColumn("Sector");               //branch 15
    man->CreateNtupleDColumn("MomentumLoss");         //branch 16
    man->CreateNtupleDColumn("EnergyDeposition");     //branch 17
    man->CreateNtupleDColumn("TOF");                  //branch 18
    man->CreateNtupleDColumn("p_reconstructed");      //branch 19 - momentum determined by TOF
    man->CreateNtupleDColumn("delta_p");              //branch 20 - (p(rec)-p(gen) / p(gen)) [%]
    man->FinishNtuple();
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHallBRunAction::~NeutronHallBRunAction()
{
//    G4cout << "ending NeutronHallBRunAction" << G4endl;
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHallBRunAction::BeginOfRunAction(const G4Run*)
{ 
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

    
    // Get analysis manager
    G4AnalysisManager* man = G4AnalysisManager::Instance(); // Open an output file
    man -> OpenFile("NeutronHallB");
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void NeutronHallBRunAction::EndOfRunAction(const G4Run* ){
    
//    G4cout << "NeutronHallBRunAction: end of run" << G4endl;
    G4AnalysisManager * man = G4AnalysisManager::Instance();
    man -> Write();
    man -> CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
