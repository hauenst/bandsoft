// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBDetectorConstruction.cc
/// \brief Implementation of the NeutronHallBDetectorConstruction class

#include "NeutronHallBDetectorConstruction.h"
#include "NeutronHallBDetectorMessenger.h"
#include "NeutronHallBTrackerSD.h"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4GenericMessenger.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"

#include "G4Sphere.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Trd.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"

#include "G4SystemOfUnits.hh"



G4ThreadLocal G4GlobalMagFieldMessenger* NeutronHallBDetectorConstruction::fMagFieldMessenger = 0;

NeutronHallBDetectorConstruction::NeutronHallBDetectorConstruction(void * t)
	:G4VUserDetectorConstruction(),
	fNbOfChambers(0),
	fLogicTarget(NULL), fLogicChamber(NULL),
	fTargetMaterial(NULL), fChamberMaterial(NULL),
	fStepLimit(NULL),
	fVisAttributes(),
	fCheckOverlaps(true){
		fMessenger = new NeutronHallBDetectorMessenger(this);

		fNbOfChambers = 5;
		fLogicChamber = new G4LogicalVolume*[fNbOfChambers];

		treePtr = t;
	}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NeutronHallBDetectorConstruction::~NeutronHallBDetectorConstruction(){
	delete [] fLogicChamber;
	delete fStepLimit;
	delete fMessenger;
	for (auto visAttributes: fVisAttributes) {
		delete visAttributes;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* NeutronHallBDetectorConstruction::Construct(){
	// Define materials
	DefineMaterials();

	// Define volumes
	return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHallBDetectorConstruction::DefineMaterials(){
	// Material definition

	G4NistManager* man = G4NistManager::Instance();

	// Define StainlessSteel not in NIST, following http://hypernews.slac.stanford.edu/HyperNews/geant4/get/geometry/915.html?inline=1
	G4Element* C  = man->FindOrBuildElement("C");
	G4Element* Si = man->FindOrBuildElement("Si");
	G4Element* Cr = man->FindOrBuildElement("Cr");
	G4Element* Mn = man->FindOrBuildElement("Mn");
	G4Element* Fe = man->FindOrBuildElement("Fe");
	G4Element* Ni = man->FindOrBuildElement("Ni");

	G4Material* StainlessSteel = new G4Material( "StainlessSteel", 8.06*g/cm3, 6 );
	StainlessSteel->AddElement(C, 0.001);
	StainlessSteel->AddElement(Si, 0.007);
	StainlessSteel->AddElement(Cr, 0.18);
	StainlessSteel->AddElement(Mn, 0.01);
	StainlessSteel->AddElement(Fe, 0.712);
	StainlessSteel->AddElement(Ni, 0.09);

	G4Material* Air = man->FindOrBuildMaterial("G4_AIR");
	G4Material* vacuum = new G4Material("Vacuum", 1.e-25*g/cm3, 1, kStateGas, 2.e-2*bar, CLHEP::STP_Temperature);
	vacuum -> AddMaterial(Air, 1.);



	// Print materials
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* NeutronHallBDetectorConstruction::DefineVolumes(){
	G4NistManager* nist = G4NistManager::Instance();
	G4Material* vacuum          = nist->FindOrBuildMaterial("Vacuum");
	G4Material* default_mat     = nist->FindOrBuildMaterial("G4_AIR");
	G4Material* SS_mat          = nist->FindOrBuildMaterial("StainlessSteel");
	G4Material* LeadWall_mat    = nist->FindOrBuildMaterial("G4_Pb");
	G4Material* BAND_mat        = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

	G4double length;
	G4double shift;

	// ------------------ Defining World ---------------- //
	G4double world_sizeX = 5000.*cm;
	G4double world_sizeY = 5000.*cm;
	G4double world_sizeZ = 5000.*cm;


	G4Box* solidWorld =    
		new G4Box("World",                                     // its name
				0.5*world_sizeX, 
				0.5*world_sizeY, 
				0.5*world_sizeZ);                            // its size

	G4LogicalVolume* logicWorld =                         
		new G4LogicalVolume(solidWorld,                        // its solid
				vacuum,                            // its world material
				"World");                          // its name

	G4VPhysicalVolume* physWorld = 
		new G4PVPlacement(0,                                   // no rotation
				G4ThreeVector(),                     // at (0,0,0)
				logicWorld,                          // its logical volume
				"World",                             // its name
				0,                                   // its mother  volume
				false,                               // no boolean operation
				0,                                   // copy number
				fCheckOverlaps);                     // overlaps checking
	// ------------------ End Defining World ---------------- //

	// ------------------ Defining Target - in World ---------------- //
	/*
	   G4double targetRadius = 1.*mm;

	   G4Sphere* solid_target = 
	   new G4Sphere("target",                                 // name
	   0.,                                      // inner radius 
	   targetRadius,                            // outer radius
	   0.,                                      // starting phi angle
	   CLHEP::twopi,                            // ending phi angle
	   0.,                                      // starting theta angle
	   CLHEP::twopi/2.);                        // ending theta angle

	   G4LogicalVolume* logic_target = 
	   new G4LogicalVolume(solid_target,                      // its solid
	   default_mat,                       // its material
	   "target");                         // its name

	   new G4PVPlacement(0,                                   // no rotation
	   G4ThreeVector(),                     // at 0,0,0
	   logic_target,                        // its logical volume
	   "target",                            // its name
	   logicWorld,                          // its mother volume
	   false,                               // no boolean operation
	   0,                                   // copy number
	   fCheckOverlaps);                     // overlap checkign

	   logic_target -> SetVisAttributes (G4Colour::G4Colour( 0.65 , 0.82 , 0.77 ));
	   */
	// ------------------ End Defining Target ---------------- //

	// ------------------ Defining BAND Detector - in World ---------------- //

	const int BAND_wallReplicas = 5;
	const int BAND_wallGroups = 5;
	//const int BAND_numRows = 18;
	//
	//double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};
	//

	G4double BAND_barCrossSection = 7.2*cm;
	G4double BAND_groupA_length = 1637.*mm; // temp
	G4double BAND_groupB_length = 2019.*mm; // temp
	G4double BAND_groupC_length = 2019.*mm; // temp
	G4double BAND_groupD_length = 512.*mm;
	G4double BAND_groupE_length = 2019.*mm; // temp
	G4double BAND_groupD_offset = (491.35 + BAND_groupD_length/2.)*mm; // temp  

	G4double BAND_Z_start = -2670.*mm;
	G4double BAND_Z_thickness = BAND_barCrossSection * BAND_wallReplicas;
	G4double BAND_Z_coord_center = BAND_Z_start - (BAND_Z_thickness/2.)*mm;
	/* old below - above is updated 4/25/2019
	   G4double BAND_barCrossSection = 7.4*cm;
	   G4double BAND_groupA_length = 1634.58*mm; // temp
	   G4double BAND_groupB_length = 1946.53*mm; // temp
	   G4double BAND_groupC_length = 2018.35*mm; // temp
	   G4double BAND_groupD_length = 500.*mm;
	   G4double BAND_groupE_length = 2018.35*mm; // temp
	   G4double BAND_groupD_offset = (509.17 + BAND_groupD_length/2.)*mm; // temp  

	   G4double BAND_Z_coord_center = -2805.*mm;
	   G4double BAND_Z_thickness = BAND_barCrossSection * BAND_wallReplicas;
	   G4double BAND_Z_start = BAND_Z_coord_center + (BAND_Z_thickness/2.)*mm;
	   */
	G4double BAND_Y_offset =  -(BAND_barCrossSection * 5)*mm;

	G4int numBars;
	G4double wall_Z_offset;
	// create 5 box walls that will hold the BAND array
	for (int wall = 0; wall < BAND_wallReplicas; wall++){
		wall_Z_offset = (BAND_Z_start - BAND_barCrossSection*(1./2+wall)*mm );

		G4Box* solidDetectorBarsE =
			new G4Box("DetectorBars_E",
					0.5*BAND_groupE_length,
					0.5*BAND_barCrossSection,
					0.5*BAND_barCrossSection);
		logicDetectorBarsE[wall] = 
			new G4LogicalVolume(solidDetectorBarsE,
					BAND_mat,
					"DetectorBarsLV_E"); 

		G4Box* solidDetectorBarsD =
			new G4Box("DetectorBars_D",
					0.5*BAND_groupD_length,
					0.5*BAND_barCrossSection,
					0.5*BAND_barCrossSection);
		logicDetectorBarsD[wall] = 
			new G4LogicalVolume(solidDetectorBarsD,
					BAND_mat,
					"DetectorBarsLV_D");   

		G4Box* solidDetectorBarsC =
			new G4Box("DetectorBars_C",
					0.5*BAND_groupC_length,
					0.5*BAND_barCrossSection,
					0.5*BAND_barCrossSection);
		logicDetectorBarsC[wall] = 
			new G4LogicalVolume(solidDetectorBarsC,
					BAND_mat,
					"DetectorBarsLV_C");

		G4Box* solidDetectorBarsB =
			new G4Box("DetectorBars_B",
					0.5*BAND_groupB_length,
					0.5*BAND_barCrossSection,
					0.5*BAND_barCrossSection);
		logicDetectorBarsB[wall] = 
			new G4LogicalVolume(solidDetectorBarsB,
					BAND_mat,
					"DetectorBarsLV_B");    

		G4Box* solidDetectorBarsA =
			new G4Box("DetectorBars_A",
					0.5*BAND_groupA_length,
					0.5*BAND_barCrossSection,
					0.5*BAND_barCrossSection);
		logicDetectorBarsA[wall] = 
			new G4LogicalVolume(solidDetectorBarsA,
					BAND_mat,
					"DetectorBarsLV_A");        

		// now for each wall, need to parameterize it into the 5 groups.
		G4double offset = 0.;
		G4double prevNumBars = 0.;
		for (int group = 0; group < BAND_wallGroups; group++){

			if    (group == 0){
				numBars = 2;
				for(int barNo = 0; barNo < numBars; barNo++){
					new G4PVPlacement(0, G4ThreeVector(0,BAND_Y_offset + (1/2. + barNo)*BAND_barCrossSection, wall_Z_offset), logicDetectorBarsE[wall],                 
							"BarsE", logicWorld, false, barNo+24*wall, fCheckOverlaps);
				}
				offset +=  numBars *BAND_barCrossSection*mm;
			}

			else if(group == 1){
				prevNumBars += numBars;
				numBars = 6;
				for(int barNo = 0; barNo < numBars; barNo++){
					new G4PVPlacement(0, G4ThreeVector(-BAND_groupD_offset,BAND_Y_offset + offset + (1/2. + barNo)*BAND_barCrossSection,wall_Z_offset ), logicDetectorBarsD[wall],                 
							"BarsD1", logicWorld, false, (prevNumBars + barNo)+24*wall , fCheckOverlaps);
					new G4PVPlacement(0, G4ThreeVector(+BAND_groupD_offset,BAND_Y_offset + offset + (1/2. + barNo)*BAND_barCrossSection,wall_Z_offset ), logicDetectorBarsD[wall],                 
							"BarsD2", logicWorld, false, (prevNumBars+(barNo+numBars))+24*wall, fCheckOverlaps);
				}
				offset +=  numBars *BAND_barCrossSection*mm;
			}

			else if(group == 2){
				prevNumBars += numBars*2 ;
				numBars = 4;
				for(int barNo = 0; barNo < numBars; barNo++){
					new G4PVPlacement(0, G4ThreeVector(0,BAND_Y_offset + offset+ (1/2. + barNo)*BAND_barCrossSection,wall_Z_offset ), logicDetectorBarsC[wall],                 
							"BarsC", logicWorld, false, (prevNumBars+barNo)+24*wall, fCheckOverlaps);
				}
				offset +=  numBars *BAND_barCrossSection*mm;
			}

			else if(group == 3){
				prevNumBars += numBars ;
				numBars = 3;
				for(int barNo = 0; barNo < numBars; barNo++){
					new G4PVPlacement(0, G4ThreeVector(0,BAND_Y_offset + offset+ (1/2. + barNo)*BAND_barCrossSection,wall_Z_offset ), logicDetectorBarsB[wall],                 
							"BarsB", logicWorld, false, (prevNumBars+barNo)+24*wall, fCheckOverlaps);
				}
				offset +=  numBars *BAND_barCrossSection*mm;
			}
			else if(group == 4){
				prevNumBars += numBars ;
				numBars = 3;
				for(int barNo = 0; barNo < numBars; barNo++){
					new G4PVPlacement(0, G4ThreeVector(0,BAND_Y_offset + offset+ (1/2. + barNo)*BAND_barCrossSection,wall_Z_offset ), logicDetectorBarsA[wall],                 
							"BarsA", logicWorld, false, (prevNumBars+barNo)+24*wall, fCheckOverlaps);
				}
			}
		}
	}
	// ------------------ END Defining BAND Detector  ---------------- //
	
	// ------------------ Defining Lead Wall  - in World ---------------- //
	/* OLD Pb WALL GEOMETRY -- UPDATED 4/29/2019
	G4double lead_thickness = 25.*mm/2;
	G4double lead_Z_coord = (BAND_Z_coord_center + BAND_Z_thickness/2.)*mm + 15*mm + lead_thickness;

	G4double RinDetector = 225.*mm;
	G4double RoutDetector = 1300.*mm;
	G4Tubs* solidLeadWall = 
		new G4Tubs("LeadWall",                                   // name
				RinDetector,                                 // inner radius
				RoutDetector,                                // outer raidus
				lead_thickness,                                    // z extent
				0.,                                          // starting phi angle
				CLHEP::twopi);                            // ending phi angle

	logicLeadWall = 
		new G4LogicalVolume(solidLeadWall,                       // its solid
				LeadWall_mat,                        // material
				"LeadWall");                         // its name

	new G4PVPlacement(0,                                     // no rotation
			G4ThreeVector(0, 0,lead_Z_coord),      // its (0,0,z)
			logicLeadWall,                         // logical volume
			"LeadWall",                            // name
			logicWorld,                            // mother volume
			false,                                 // no boolean operation
			0,                                     // copy number
			fCheckOverlaps);                       // overlap checking

	G4Colour LeadWall_Color( 0.25 , 0.38 , 0.25 );
	G4VisAttributes* LeadWallVisAttributes = new G4VisAttributes( true , LeadWall_Color );
	fVisAttributes.push_back(LeadWallVisAttributes);
	logicLeadWall -> SetVisAttributes(LeadWallVisAttributes);
	*/
	/*
	double BAND_pbCrossSection = 20.*mm;
	double Pbwall_Z_offset = (BAND_Z_start - BAND_barCrossSection*0.5*mm) + 174.63*mm + 73.0*mm; 
			// BAND_Z_start is face of layer 5, closest to target. Then subtract off
			// 1/2 layer thickness to get to the "PMT" on layer 5, and then we have
			// measurements that take us to the center of the Pb wall by Rey

	G4Box* PbBarsE =
		new G4Box("PbBars_E",
				0.5*BAND_groupE_length,
				0.5*BAND_barCrossSection,
				0.5*BAND_pbCrossSection);
	G4Box* PbBarsD =
		new G4Box("PbBars_D",
				0.5*BAND_groupD_length,
				0.5*BAND_barCrossSection,
				0.5*BAND_pbCrossSection);
	G4Box* PbBarsC =
		new G4Box("PbBars_C",
				0.5*BAND_groupC_length,
				0.5*BAND_barCrossSection,
				0.5*BAND_pbCrossSection);
	G4Box* PbBarsB =
		new G4Box("PbBars_B",
				0.5*BAND_groupB_length,
				0.5*BAND_barCrossSection,
				0.5*BAND_pbCrossSection);
	G4Box* PbBarsA =
		new G4Box("PbBars_A",
				0.5*BAND_groupA_length,
				0.5*BAND_barCrossSection,
				0.5*BAND_pbCrossSection);
	G4LogicalVolume * logicalPbBarsE = 
		new G4LogicalVolume(PbBarsE,
				LeadWall_mat,
				"PbBarsLV_E"); 
	G4LogicalVolume * logicalPbBarsD = 
		new G4LogicalVolume(PbBarsD,
				LeadWall_mat,
				"PbBarsLV_D"); 
	G4LogicalVolume * logicalPbBarsC = 
		new G4LogicalVolume(PbBarsC,
				LeadWall_mat,
				"PbBarsLV_C"); 
	G4LogicalVolume * logicalPbBarsB = 
		new G4LogicalVolume(PbBarsB,
				LeadWall_mat,
				"PbBarsLV_B"); 
	G4LogicalVolume * logicalPbBarsA = 
		new G4LogicalVolume(PbBarsA,
				LeadWall_mat,
				"PbBarsLV_A"); 

	G4double offset = 0;
	int prevNumBars = 0;
	for (int group = 0; group < BAND_wallGroups; group++){

		if    (group == 0){
			numBars = 2;
			for(int barNo = 0; barNo < numBars; barNo++){
				new G4PVPlacement(	0, 
							G4ThreeVector(0, BAND_Y_offset + (1/2. + barNo)*BAND_barCrossSection, Pbwall_Z_offset), 
							logicalPbBarsE,                 
							Form("PbBarsE_%i",barNo), 
							logicWorld, 
							false, 
							barNo, 
							fCheckOverlaps);
			}
			offset +=  numBars *BAND_barCrossSection*mm;
		}

		else if(group == 1){
			prevNumBars += numBars;
			numBars = 6;
			for(int barNo = 0; barNo < numBars; barNo++){
				new G4PVPlacement(0, 
						G4ThreeVector(-BAND_groupD_offset,BAND_Y_offset + offset + (1/2. + barNo)*BAND_barCrossSection,Pbwall_Z_offset ), 
						logicalPbBarsD,                 
						Form("PbBarsD1_%i",barNo), 
						logicWorld, 
						false, 
						(prevNumBars + barNo), 
						fCheckOverlaps);
				new G4PVPlacement(0, 
						G4ThreeVector(+BAND_groupD_offset,BAND_Y_offset + offset + (1/2. + barNo)*BAND_barCrossSection,Pbwall_Z_offset ), 
						logicalPbBarsD,                 
						Form("PbBarsD2_%i",barNo), 
						logicWorld, 
						false, 
						(prevNumBars+(barNo+numBars)), 
						fCheckOverlaps);
			}
			offset +=  numBars *BAND_barCrossSection*mm;
		}

		else if(group == 2){
			prevNumBars += numBars*2 ;
			numBars = 4;
			for(int barNo = 0; barNo < numBars; barNo++){
				new G4PVPlacement(	0, 
							G4ThreeVector(0,BAND_Y_offset + offset+ (1/2. + barNo)*BAND_barCrossSection,Pbwall_Z_offset ), 
							logicalPbBarsC,                 
							Form("PbBarsC_%i",barNo), 
							logicWorld, 
							false, 
							(prevNumBars+barNo), 
							fCheckOverlaps);
			}
			offset +=  numBars *BAND_barCrossSection*mm;
		}

		else if(group == 3){
			prevNumBars += numBars ;
			numBars = 3;
			for(int barNo = 0; barNo < numBars; barNo++){
				new G4PVPlacement(	0, 
							G4ThreeVector(0,BAND_Y_offset + offset+ (1/2. + barNo)*BAND_barCrossSection,Pbwall_Z_offset ), 
							logicalPbBarsB,                 
							Form("PbBarsB_%i",barNo), 
							logicWorld, 
							false, 
							(prevNumBars+barNo), 
							fCheckOverlaps);
			}
			offset +=  numBars *BAND_barCrossSection*mm;
		}
		else if(group == 4){
			prevNumBars += numBars ;
			numBars = 3;
			for(int barNo = 0; barNo < numBars; barNo++){
				new G4PVPlacement(	0, 
							G4ThreeVector(0,BAND_Y_offset + offset+ (1/2. + barNo)*BAND_barCrossSection,Pbwall_Z_offset ), 
							logicalPbBarsA,                 
							Form("PbBarsA_%i",barNo), 
							logicWorld, 
							false, (prevNumBars+barNo), 
							fCheckOverlaps);
			}
		}
	}

	*/
	// ------------------ END Defining Lead Wall  - in World ---------------- //

	/*
	// ------------------ Defining Solenoid ---------------- //

	length = 480.*mm;                                          // length of the first tube pipe

	G4Tubs* solid_solenoid =     
		new G4Tubs("Solenoid",                                 // name
				476.*mm ,                                  // inner radius
				1200.*mm,                                  // outer radius
				length ,                                   // z-extent
				0.,                                        // starting phi
				CLHEP::twopi);                             // ending phi

	logic_solenoid =      
		new G4LogicalVolume(solid_solenoid,                    // solid logic volume
				SS_mat,                          // material
				"Solenoid");                       // name

	new G4PVPlacement(0,                                   // rotation
			G4ThreeVector(),                     // at (0,0,0)
			logic_solenoid,                      // logical volume
			"Solenoid",                          // name
			logicWorld,                          // mother volume
			false,                               // no boolean operation
			0,                                   // copy number
			fCheckOverlaps);                     // checking overlap

	logic_solenoid -> SetVisAttributes (G4Colour::G4Colour( 0.99 , 0.88 , 0.66 ));


	const int NTOFs = 4; //CAUTION -->                         // if you edit this must edit .h file  
	G4Tubs                *SolidTOF[NTOFs];                    // solid 
	G4VPhysicalVolume     *PhysTOF[NTOFs];                     // the physical volume
	G4ThreeVector         posTOF[NTOFs];                       // the vector for the position of the boxes

	for (int i_tof = 0; i_tof < NTOFs; i_tof++) {              // for each TOF we want to add
		posTOF[i_tof]     = G4ThreeVector();                   // defining position vector for box

		SolidTOF[i_tof] =                                      // solid volume for each box
			new G4Tubs("TOF",                                  // name
					(252.+50.8*i_tof)*mm,                          // inner radius
					(252.+50.8*(i_tof+1))*mm,                      // outer radius
					length ,                                       // z-extent
					0.,                                            // starting phi
					CLHEP::twopi);                                 // ending phi 

		LogicTOF[i_tof] =                                      // logical volume for each
			new G4LogicalVolume(SolidTOF[i_tof],               // TOF is solid
					BAND_mat,                      // material for each TOF
					"TOF" );                       // name for each logical volume TOF

		LogicTOF[i_tof]->SetVisAttributes (G4Colour::G4Colour( 0.75+(i_tof/10.) , 0.6-(i_tof/10.) , 0.75 ));

		PhysTOF[i_tof]=                                        // physical volume for TOF
			new G4PVPlacement(0,                               // rotation
					G4ThreeVector(),                 // at 0,0,0
					LogicTOF[i_tof],                 // logical volume of each TOF
					"TOF",                           // name of each TOF
					logicWorld,                      // mother volume for TOF
					false,                           // no boolean operation
					0,                               // copy number
					fCheckOverlaps);                 // checking overlap
	}

	length = 1200/2.*mm;
	shift = 480*mm;

	// the 3 CND lightguides are 144.97 degrees (35.03 degrees)
	// the 1 CND TOF is 158.19 degrees (21.81 degrees)
	const int slanted_TOFs = 4; // CAUTION -->                    // if you edit this, must edit .h file

	G4Cons                *Solid_slantedTOF_front[slanted_TOFs];  // solid TOF
	G4VPhysicalVolume     *Phys_slantedTOF_front[slanted_TOFs];   // the physical volume
	G4ThreeVector         pos_slantedTOF_front[slanted_TOFs];     // the vector for the position of the TOF

	for (int i_tof = 0; i_tof < slanted_TOFs; i_tof++) {
		pos_slantedTOF_front[i_tof]     = G4ThreeVector();        // defining position vector for TOF

		if(i_tof==0){  // This is the CND TOF
			Solid_slantedTOF_front[i_tof] =                       // solid volume for CND
				new G4Cons("slanted_TOF_front",                   // name
						(445.84+252.)*mm,                         // inner radius 1
						(445.84+252.+50.8)*mm,                    // outer raidus 1
						(252.)*mm,                             // inner radius 2
						(252.+50.8)*mm,                        // outer radius 2
						length ,                               // z-extent
						0.,                                    // starting phi
						CLHEP::twopi);                         // ending phi             
		}
		else{           // These are the CND lightguides
			Solid_slantedTOF_front[i_tof] =                       // solid volume for each TOF
				new G4Cons("slanted_TOF_front",                   // name
						(688.8+252.+50.8*(i_tof))*mm,            // inner radius 1
						(688.8+252.+50.8*(i_tof+1))*mm,          // outer radius 1
						(252.+50.8*(i_tof))*mm,                // inner radius 2
						(252.+50.8*(i_tof+1))*mm,              // outer radius 2
						length ,                               // z-extent
						0.,                                    // starting phi
						CLHEP::twopi);                         // ending phi 
		}

		Logic_slantedTOF_front[i_tof] =                           // logical volume for each TOF
			new G4LogicalVolume(Solid_slantedTOF_front[i_tof],    // TOF is solid
					BAND_mat,                         // material for each TOF
					"slanted_TOF_front" );            // name for each logical volume TOF

		Logic_slantedTOF_front[i_tof]->SetVisAttributes(G4Colour::G4Colour( 0.75 , 0.6+(i_tof/10.) , 0.75-(i_tof/10.) ));   

		Phys_slantedTOF_front[i_tof]=                              // physical volume for TOF
			new G4PVPlacement(0,                                   // rotation
					G4ThreeVector(0,0,-1.*(shift+length)),// position for each TOF
					Logic_slantedTOF_front[i_tof],       // logical volume of each TOF
					"slanted_TOF_front",                 // name of each TOF
					logicWorld,                          // mother volume for TOF
					false,                               // no boolean operation
					0,                                   // copy number
					fCheckOverlaps);                     // checking overlap
	}

	// ------------------ END Defining Solenoid  - in World ---------------- //
	*/
	// ------------------ Defining SS Tubing around beam pipe ---------------- //
	/*
	length = 658.*mm; // length of the first tube pipe
	G4Tubs* solid_SSTube_P1 =     
		new G4Tubs("SSTube_P1",                                // name
				191*mm,                                  // inner radius
				194*mm,                                   // outer radius
				length ,                                   // z-extent
				0.,                                        // starting phi
				CLHEP::twopi);                             // ending phi

	logic_SSTube_P1 =      
		new G4LogicalVolume(solid_SSTube_P1,                   // solid logic volume
				SS_mat,                            // material
				"SSTube_P1");                      // name

	new G4PVPlacement(0,                                   // rotation
			G4ThreeVector(),                     // at (0,0,0)
			logic_SSTube_P1,                     // logical volume
			"SSTube_P1",                         // name
			logicWorld,                          // mother volume
			false,                               // no boolean operation
			0,                                   // copy number
			fCheckOverlaps);                     // checking overlap

	//logic_SSTube_P1 -> SetVisAttributes (G4Colour::G4Colour( 0.99 , 0.88 , 0.66 ));

	//  ----------1 cm plastic coating around tub to represent cables ---------- //
	G4Tubs* solid_SSTube_P1_coating = 
		new G4Tubs("SSTube_P1_coating",                        // name
				195*mm,                                  // inner radius -- starting at outer radius of 1st pipe
				205*mm,                                   // outer radius
				length ,                                   // z extent
				0.,                                        // starting phi
				CLHEP::twopi);                             // ending phi

	logic_SSTube_P1_coating =      
		new G4LogicalVolume(solid_SSTube_P1_coating,           // solid logic volume
				BAND_mat ,                         // material
				"SSTube_P1_coating");              // name

	new G4PVPlacement(0,                                   // no rotation
			G4ThreeVector(),                     // at (0,0,0)
			logic_SSTube_P1_coating,             // logical volume
			"SSTube_P1_coating",                 // name
			logicWorld,                          // mother volume
			false,                               // no boolean operation
			0,                                   // copy number
			fCheckOverlaps);                     // checking overlaps

	//logic_SSTube_P1_coating -> SetVisAttributes (G4Colour::G4Colour( 0.75 , 0.6 , 0.75 ));

	//  ---------- SS connector from target tube to beam tube ---------- //
	length = 24.*mm;
	shift = 658.*mm;

	G4Tubs* solid_SSTube_Connector1 =         
		new G4Tubs("SSTube_Connector1",                        // name
				194*mm,                                  // inner radius
				221*mm,                                   // outer radius
				length,                                   // z extent
				0.,                                        // starting phi
				CLHEP::twopi);                             // ending phi

	logic_SSTube_Connector1 = 
		new G4LogicalVolume(solid_SSTube_Connector1,           // solid object
				SS_mat ,                           // material
				"SSTube_Connector1");              // name

	new G4PVPlacement(0,                                   // no rotation
			G4ThreeVector(0,0,-(length+shift)),  // at (x,y,z)
			logic_SSTube_Connector1,             // logical volume
			"SSTube_Connector1",                 // name
			logicWorld,                          // mother volume
			false,                               // boolean operation  
			0,                                   // copy number
			fCheckOverlaps);                     // checking overlaps

	//logic_SSTube_Connector1 -> SetVisAttributes (G4Colour::G4Colour( 0.99 , 0.88 , 0.66 ));

	//  ---------- SS connector #2 from target tube to beam tube ---------- //
	length = 34.*mm;
	shift = 658.*mm;

	G4Tubs* solid_SSTube_Connector2 = 
		new G4Tubs("SSTube_Connector2",
				221*mm,
				239*mm,
				length,
				0.,
				CLHEP::twopi);

	logic_SSTube_Connector2 = 
		new G4LogicalVolume(solid_SSTube_Connector2, 
				SS_mat, 
				"SSTube_Connector2");

	new G4PVPlacement(0, 
			G4ThreeVector(0,0,-(length+shift)),
			logic_SSTube_Connector2, 
			"SSTube_Connector2",
			logicWorld, 
			false,  
			0,  
			fCheckOverlaps);

	//logic_SSTube_Connector2 -> SetVisAttributes (G4Colour::G4Colour( 0.99 , 0.38 , 0.2 ));

	//  ---------- SS connector #2 coating to represent cables ---------- //
	G4Tubs* solid_SSTube_Connector2_coating = 
		new G4Tubs("SSTube_Connector2_coating", 
				240*mm, 
				250*mm, 
				length, 
				0., 
				CLHEP::twopi);

	logic_SSTube_Connector2_coating = 
		new G4LogicalVolume(solid_SSTube_Connector2_coating, 
				BAND_mat, 
				"SSTube_Connector2_coating");

	new G4PVPlacement(0, 
			G4ThreeVector(0,0,-(length+shift)),
			logic_SSTube_Connector2_coating, 
			"SSTube_Connector2_coating", 
			logicWorld, 
			false,  
			0,  
			fCheckOverlaps);

	//logic_SSTube_Connector2_coating -> SetVisAttributes (G4Colour::G4Colour( 0.75 , 0.6 , 0.75 ));

	//  ---------- SS tube from connector #2 to BAND detector ---------- //
	length = 3304*mm+BAND_Z_thickness/2.*mm-1750*mm;
	shift = 706.*mm;
	G4Tubs* solid_SSTube_P2 = 
		new G4Tubs("SSTube_P2", 
				208*mm, 
				211*mm, 
				length, 
				0., 
				CLHEP::twopi);

	logic_SSTube_P2 = 
		new G4LogicalVolume(solid_SSTube_P2, 
				SS_mat, 
				"SSTube_P2");

	new G4PVPlacement(0, 
			G4ThreeVector(0,0,-(length+shift)),
			logic_SSTube_P2, 
			"SSTube_P2", 
			logicWorld, 
			false,  
			0,  
			fCheckOverlaps);

	//logic_SSTube_P2 -> SetVisAttributes (G4Colour::G4Colour( 0.99 , 0.88 , 0.66 ));

	//  ---------- SS tube coatings from connector #2 to BAND detector ---------- //
	G4Tubs* solid_SSTube_P2_coating = 
		new G4Tubs("SSTube_P2_coating", 
				211.*mm, 
				216*mm, 
				length, 
				0., 
				CLHEP::twopi);

	logic_SSTube_P2_coating = 
		new G4LogicalVolume(solid_SSTube_P2_coating, 
				BAND_mat, 
				"SSTube_P2_coating");

	new G4PVPlacement(0, 
			G4ThreeVector(0,0,-(length+shift)),  
			logic_SSTube_P2_coating, 
			"SSTube_P2_coating", 
			logicWorld, 
			false,  
			0,  
			fCheckOverlaps);

	//logic_SSTube_P2_coating -> SetVisAttributes (G4Colour::G4Colour( 0.75 , 0.6 , 0.75 ));
	// ------------------ END SS Tubing around beam pipe ---------------- //

	// ------------------ Defining Electronic Boxes  - in World ---------------- //

	G4double BoxLength = 20*cm;
	G4double BoxWidth = 25*cm;
	G4double BoxThickness = 10*cm;
	G4double BoxZ = -150*cm;
	G4double BoxDisFromAxis = 240.*mm + BoxWidth/2. + 5*mm;

	G4double BoxX;
	G4double BoxY;
	G4double BoxAngle;

	const int NBoxes = 6;      //if edit this, must edit .h file  // Creating 6 boxes that will go around the beam line
	G4Box                 *SolidElBox[NBoxes];                 // solid box
	//G4LogicalVolume       *LogicElBox[NBoxes];               // the logical volumes
	G4VPhysicalVolume     *PhysElBox[NBoxes];                  // the physical volume
	G4ThreeVector         posElBox[NBoxes];                    // the vector for the position of the boxes
	//G4VisAttributes       *ElBoxVisAtt = new G4VisAttributes( true , G4Colour::G4Colour( 0.75 , 0.6 , 0.75 ));   // setting color
	//fVisAttributes.push_back(ElBoxVisAtt);
	G4RotationMatrix      rot_mat;                             // rotation matrix to rotate each box
	G4Transform3D         transform;                           // the rotation that will be applied to the placement

	for (int i_box = 0; i_box < NBoxes; i_box++) {
	BoxAngle            = i_box * 60 * deg;                               // this box's angle around beam line
	BoxX                = BoxDisFromAxis * std::cos(BoxAngle);            // the x-placement based on angle
	BoxY                = BoxDisFromAxis * std::sin(BoxAngle);            // the y-placement based on angle
	rot_mat             = G4RotationMatrix();                             // initializing rotation matrix
	rot_mat.rotateZ(BoxAngle);                                            // defining rotation matrix around z-beamLine for box
	posElBox[i_box]     = G4ThreeVector( BoxX , BoxY , BoxZ);             // defining position vector for box
	transform           = G4Transform3D( rot_mat , posElBox[i_box] );     // creating rotation of each box 

	SolidElBox[i_box] =                                    // solid volume for each box
	new G4Box("ElBox",                                 // name
	BoxLength/2.,                             // length of box
	BoxWidth/2.,                              // width of box
	BoxThickness/2.);                         // thickness of box

	LogicElBox[i_box] =                                    // logical volume for each box
	new G4LogicalVolume(SolidElBox[i_box],             // box is solid
	BAND_mat,                        // material for each box
	"ElBoxLV" );                     // name for each logical volume box

	//LogicElBox[i_box]->SetVisAttributes (ElBoxVisAtt);     // setting color

	PhysElBox[i_box]=                                      // physical volume for box
	new G4PVPlacement(transform,                       // rotation for each box
	LogicElBox[i_box],                 // logical volume of each box
	"ElBox",                           // name of each box
	logicWorld,                        // mother volume for box
	false,                             // no boolean operation
	0,                                 // copy number
	fCheckOverlaps);                   // checking overlap
	}
	// ------------------ END Defining Electronic Boxes  ---------------- //
	*/

	// Always return the physical world
	return physWorld;
}

void NeutronHallBDetectorConstruction::ConstructSDandField()
{
	// Sensitive detectors
	auto sdManager = G4SDManager::GetSDMpointer();
	G4String SDname;

	NeutronHallBTrackerSD* bandSD = new NeutronHallBTrackerSD(treePtr,SDname="band");
	sdManager->AddNewDetector(bandSD);

	for (int wall = 0; wall < 5; wall++) {
		logicDetectorBarsE[wall]->SetSensitiveDetector(bandSD);
		logicDetectorBarsD[wall]->SetSensitiveDetector(bandSD);
		logicDetectorBarsC[wall]->SetSensitiveDetector(bandSD); 
		logicDetectorBarsB[wall]->SetSensitiveDetector(bandSD); 
		logicDetectorBarsA[wall]->SetSensitiveDetector(bandSD);            
	}

	// Electro-magnetic fields (currently zero)  
	G4ThreeVector fieldValue = G4ThreeVector();
	G4AutoDelete::Register(fMagFieldMessenger);
}

void NeutronHallBDetectorConstruction::SetTargetMaterial(G4String materialName)
{
	G4NistManager* nistManager = G4NistManager::Instance();

	G4Material* pttoMaterial =
		nistManager->FindOrBuildMaterial(materialName);

	if (fTargetMaterial != pttoMaterial) {
		if ( pttoMaterial ) {
			fTargetMaterial = pttoMaterial;
			if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
			//        G4cout << "\n----> The target is made of " << materialName << G4endl;
		} else {
			//        G4cout << "\n-->  WARNING from SetTargetMaterial : "
			//               << materialName << " not found" << G4endl;
		}
	}
}

void NeutronHallBDetectorConstruction::SetChamberMaterial(G4String materialName)
{
	G4NistManager* nistManager = G4NistManager::Instance();

	G4Material* pttoMaterial =
		nistManager->FindOrBuildMaterial(materialName);

	if (fChamberMaterial != pttoMaterial) {
		if ( pttoMaterial ) {
			fChamberMaterial = pttoMaterial;
			for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {
				if (fLogicChamber[copyNo]) fLogicChamber[copyNo]->
					SetMaterial(fChamberMaterial);
			}
			//        G4cout << "\n----> The chambers are made of " << materialName << G4endl;
		} else {
			//        G4cout << "\n-->  WARNING from SetChamberMaterial : "
			//               << materialName << " not found" << G4endl;
		}
	}
}


void NeutronHallBDetectorConstruction::SetMaxStep(G4double maxStep)
{
	if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

void NeutronHallBDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
	fCheckOverlaps = checkOverlaps;
}  
