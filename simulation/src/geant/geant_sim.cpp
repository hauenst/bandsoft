#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "Randomize.hh"
#include "G4RunManager.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"

#include "NeutronHallBDetectorConstruction.h"
#include "NeutronHallBActionInitialization.h"

#include "QGSP_BIC.hh"
#include "QGSP_BERT.hh"
#include "FTFP_BERT.hh"


using namespace std;

int main(int argc, char ** argv)
{
	if (argc != 5)
	{
		cerr << "Wrong number of arguments! Instead use:\n"
			<< "\tgeant_sim /path/to/generated/root/file /path/to/output/file /path/to/macro [1 for vis / 0 for batch]\n\n";
		exit(-1);
	}

	char * inFileName=argv[1];
	char * outFileName=argv[2];
	G4String macroName=argv[3];
	bool useVis = atoi(argv[4]);

	// Input tree
	TFile * inFile = new TFile(inFileName);
	TTree * inTree = (TTree*) inFile->Get("MCout");
	//const int nEvents = inTree->GetEntries();
	const int nEvents = 10000;
	// Output file + tree
	TFile * outFile = new TFile(outFileName,"RECREATE");
	TTree * outTree = new TTree("PropTree","Propagation Tree");

	cout << "Beginning geant_sim initialization...\n"
		<< "\tGenerator file: " << inFileName << "\n"
		<< "\tOutput file:    " << outFileName << "\n"
		<< "\tMacro:          " << macroName << "\n";
	if (useVis)
		cout << "\tVisualisation:  ON\n";
	else
		cout << "\tVisualization:  OFF\n";

	// Initialization
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4RunManager * runManager = new G4RunManager;

	// Detector Geometry
	runManager->SetUserInitialization(new NeutronHallBDetectorConstruction((void*)outTree));

	// Physics List
	G4VModularPhysicsList* physicsList = new QGSP_BERT;//FTFP_BERT_HP;//FTFP_BERT;
	physicsList->RegisterPhysics(new G4StepLimiterPhysics());
	runManager->SetUserInitialization(physicsList);

	// User Action
	runManager->SetUserInitialization(new NeutronHallBActionInitialization(inTree,outTree));

	// Initialize G4 kernel
	runManager->Initialize();
	cout << "\n\nGeant4 is initialized!\n\n";

	// Visualization and User-Interface
	G4UIExecutive* ui = NULL;
	if (useVis)
		ui = new G4UIExecutive(1,argv);
	G4VisManager *visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();
	G4UImanager* UImanager = G4UImanager::GetUIpointer();

	// Execute macro
	G4String command = "/control/execute ";
	UImanager->ApplyCommand(command+macroName);

	// Either start visualization or do the run.
	if (useVis)
		ui->SessionStart();
	else
	{
		// We're in batch mode so let's start the run
		char runCmd[100];
		sprintf(runCmd,"/run/beamOn %d",nEvents);
		UImanager->ApplyCommand(G4String(runCmd));	
	}

	// Clean Up
	outFile->Write();
	inFile->Close();
	outFile->Close();
	delete runManager;
	if (visManager) delete visManager;
	if (ui) delete ui;

	return 0;
}
