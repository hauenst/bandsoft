// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBPrimaryGeneratorAction.hh
/// \brief Definition of the NeutronHallBPrimaryGeneratorAction class

#ifndef NeutronHallBPrimaryGeneratorAction_h
#define NeutronHallBPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4Event;
class Gen_Event;
class TTree;

class NeutronHallBPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
  NeutronHallBPrimaryGeneratorAction(void * tp);    
  virtual ~NeutronHallBPrimaryGeneratorAction();
  
  virtual void GeneratePrimaries(G4Event* );
  
  // Set methods
  void SetRandomFlag(G4bool );
  
  G4double E_kin_g;
  G4int particle_ID_g;
  G4double p_g;
  G4int getCurrentEventNum();
  
 private:
  TTree * genTree;
  G4int eventNum;
  Gen_Event *thisEvent;
};

#endif
