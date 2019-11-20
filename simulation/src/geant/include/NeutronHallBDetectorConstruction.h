// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBDetectorConstruction.hh
/// \brief Definition of the NeutronHallBDetectorConstruction class

#ifndef NeutronHallBDetectorConstruction_h
#define NeutronHallBDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "tls.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;
class G4VisAttributes;

class NeutronHallBDetectorMessenger;

/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.

class NeutronHallBDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  NeutronHallBDetectorConstruction(void * t);
  virtual ~NeutronHallBDetectorConstruction();
  
 public:
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  
  // Set methods
  void SetTargetMaterial (G4String );
  void SetChamberMaterial(G4String );
  void SetMaxStep (G4double );
  void SetCheckOverlaps(G4bool );
  
 private:
  // methods
  void DefineMaterials();
  G4VPhysicalVolume* DefineVolumes();
  
  // data members
  void * treePtr;
  G4int fNbOfChambers;
  
  // pointers to logical volumes
  G4LogicalVolume*   fLogicTarget;     
  G4LogicalVolume**  fLogicChamber;    
  G4LogicalVolume*   logicDetectorBarsE[5];
  G4LogicalVolume*   logicDetectorBarsD[5];
  G4LogicalVolume*   logicDetectorBarsC[5];
  G4LogicalVolume*   logicDetectorBarsB[5];
  G4LogicalVolume*   logicDetectorBarsA[5];
  //G4LogicalVolume*   wedge_logical[4];

  G4LogicalVolume*   logicLeadWall;
  G4LogicalVolume*   logic_solenoid;
  
  G4LogicalVolume*   LogicTOF[4];
  G4LogicalVolume*   Logic_slantedTOF_front[4];
  G4LogicalVolume*   logic_SSTube_P1;
  G4LogicalVolume*   logic_SSTube_P1_coating;
  G4LogicalVolume*   logic_SSTube_Connector1;
  //G4LogicalVolume*   logic_SSTube_Connector1_coating;
  G4LogicalVolume*   logic_SSTube_Connector2;
  G4LogicalVolume*   logic_SSTube_Connector2_coating;
  G4LogicalVolume*   logic_SSTube_P2;
  G4LogicalVolume*   logic_SSTube_P2_coating;
  G4LogicalVolume*   LogicElBox[6];

  
  G4Material*        fTargetMaterial;  // pointer to the target  material
  G4Material*        fChamberMaterial; // pointer to the chamber material


  
  G4UserLimits* fStepLimit;            // pointer to user step limits
  
  NeutronHallBDetectorMessenger*  fMessenger;   // messenger
  std::vector<G4VisAttributes*> fVisAttributes;

  static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
  // magnetic field messenger
  
  G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
