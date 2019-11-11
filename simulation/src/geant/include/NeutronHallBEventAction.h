// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBEventAction.hh
/// \brief Definition of the NeutronHallBEventAction class

#ifndef NeutronHallBEventAction_h
#define NeutronHallBEventAction_h 1

#include "G4UserEventAction.hh"
#include "TRandom2.h"

#include "globals.hh"

class TTree;

/// Event action class
class NeutronHallBEventAction : public G4UserEventAction
{
 public:
  TRandom2 * rand2;
  NeutronHallBEventAction(void * treePtr);
  virtual ~NeutronHallBEventAction();
  
  virtual void  BeginOfEventAction(const G4Event* );
  virtual void    EndOfEventAction(const G4Event* );

  //void AddEdep(G4double edep) { fEdep += edep; }

 private:
  TTree * outTree;
  //G4double fEdep;
};


#endif
