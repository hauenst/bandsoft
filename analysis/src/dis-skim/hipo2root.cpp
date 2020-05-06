#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH2.h"
#include "TH1.h"

#include "reader.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

#include "RCDB/Connection.h"

#include "constants.h"

using namespace std;

const int maxProtons	= 100;
const int maxNeutrons	= 1;

int getRunNumber( string filename );
void getEventInfo( BEvent eventInfo, double &integrated_charge, double &livetime, double &starttime );
void getElectronInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex, 
			double& time, int& charge, double& beta, double& chi2pid, int& status );
bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status,
			double lV, double lW , double E_tot);
void getProtonInfo( BParticle particles, double pid[maxProtons], TVector3 momentum[maxProtons], TVector3 vertex[maxProtons],
			double time[maxProtons], double charge[maxProtons], double beta[maxProtons], double chi2pid[maxProtons], double status[maxProtons] , int& multiplicity );
bool checkProton( int pid, TVector3 momentum, TVector3 del_vertex, double time, int charge, double beta, double chi2pid, int status, int mult );
void getNeutronInfo( BBand band_hits, int& mult, int id[maxNeutrons], double edep[maxNeutrons],
			double time[maxNeutrons], TVector3 path[maxNeutrons] , double starttime , int thisRun);
bool pointsToBand(double theta,double phi,double z_m);

void LoadGlobalShift();
void LoadRunByRunShift();

double FADC_GLOBSHIFT[600] = {0.};
double FADC_RUNBYRUNSHIFT[10000] = {0.};

double TDC_GLOBSHIFT[600] = {0.};
double TDC_RUNBYRUNSHIFT[10000] = {0.};

int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 3 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [inputFile]\n\n";
		return -1;
	}

	// Create output tree
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("skim","CLAS and BAND Physics");
	//	Event info:
	double Ebeam		= 0;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	//	Electron info:
	int ePid		= 0;
	int eCharge		= 0;
	int eStatus		= 0;
	double eTime		= 0;
	double eBeta 		= 0;
	double eChi2pid		= 0;
	double E_tot		= 0;
	double E_pcal		= 0;
	double t_e		= 0;
	double dL_e		= 0;
	double lU		= 0;
	double lV		= 0;
	double lW		= 0;
	double e_vtx		= 0;
	double e_vty		= 0;
	double e_vtz		= 0;
	double p_e		= 0;
	double theta_e		= 0;
	double phi_e		= 0;
	double q		= 0;
	double theta_q		= 0;
	double phi_q		= 0;
	double nu		= 0;
	double Q2		= 0;
	double xB		= 0;
	double W2		= 0;
	// 	Neutron info:
	int nMult		= 0;
	int barID		[maxNeutrons]= {0};
	double dL_n		[maxNeutrons]= {0.};
	double theta_n		[maxNeutrons]= {0.};
	double phi_n		[maxNeutrons]= {0.};
	double p_n		[maxNeutrons]= {0.};
	double nTime		[maxNeutrons]= {0.};
	double nEdep		[maxNeutrons]= {0.};
	outTree->Branch("Ebeam"		,&Ebeam			);
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("livetime"	,&livetime		);
	outTree->Branch("starttime"	,&starttime		);
	outTree->Branch("ePid"		,&ePid			);
	outTree->Branch("eCharge"	,&eCharge		);
	outTree->Branch("eStatus"	,&eStatus		);
	outTree->Branch("eTime"		,&eTime			);
	outTree->Branch("eBeta"		,&eBeta 		);
	outTree->Branch("eChi2pid"	,&eChi2pid		);
	outTree->Branch("E_tot"		,&E_tot			);
	outTree->Branch("E_pcal"	,&E_pcal		);
	outTree->Branch("t_e"		,&t_e			);
	outTree->Branch("dL_e"		,&dL_e			);
	outTree->Branch("lU"		,&lU			);
	outTree->Branch("lV"		,&lV			);
	outTree->Branch("lW"		,&lW			);
	outTree->Branch("e_vtx"		,&e_vtx			);
	outTree->Branch("e_vty"		,&e_vty			);
	outTree->Branch("e_vtz"		,&e_vtz			);
	outTree->Branch("p_e"		,&p_e			);
	outTree->Branch("theta_e"	,&theta_e		);
	outTree->Branch("phi_e"		,&phi_e			);
	outTree->Branch("q"		,&q			);
	outTree->Branch("theta_q"	,&theta_q		);
	outTree->Branch("phi_q"		,&phi_q			);
	outTree->Branch("nu"		,&nu			);
	outTree->Branch("Q2"		,&Q2			);
	outTree->Branch("xB"		,&xB			);
	outTree->Branch("W2"		,&W2			);
	outTree->Branch("nMult"		,&nMult			);
	outTree->Branch("barID"		,&barID			,"barID[nMult]/D"	);
	outTree->Branch("dL_n"		,&dL_n			,"dL_n[nMult]/D"	);
	outTree->Branch("theta_n"	,&theta_n		,"theta_n[nMult]/D"	);
	outTree->Branch("phi_n"		,&phi_n			,"phi_n[nMult]/D"	);
	outTree->Branch("p_n"		,&p_n			,"p_n[nMult]/D"		);
	outTree->Branch("nTime"		,&nTime			,"nTime[nMult]/D"	);
	outTree->Branch("nEdep"		,&nEdep			,"nEdep[nMult]/D"	);
	

	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

	// Load FADC offsets
	LoadGlobalShift();
	LoadRunByRunShift();

	// Load input file
	for( int i = 2 ; i < argc ; i++ ){
		// Using run number of current file, grab the beam energy from RCDB
		int runNum = getRunNumber(argv[i]);
		auto cnd = connection.GetCondition(runNum, "beam_energy");
		Ebeam = cnd->ToDouble() / 1000. * 1.018; // [GeV] -- conversion factor due to miscalibration in RCDB

		// Setup hipo reading for this file
		TString inputFile = argv[i];
		hipo::reader reader;
		reader.open(inputFile);
		hipo::dictionary  factory;      
		hipo::schema	  schema;
		reader.readDictionary(factory); 
		BEvent		event_info		(factory.getSchema("REC::Event"		));
		BParticle	particles		(factory.getSchema("REC::Particle"	));
		BCalorimeter	calorimeter		(factory.getSchema("REC::Calorimeter"	));
		BScintillator	scintillator		(factory.getSchema("REC::Scintillator"	));
		BBand		band_hits		(factory.getSchema("BAND::hits"		));
		hipo::bank	band_adcs		(factory.getSchema("BAND::adc"		));
		hipo::bank	band_tdcs		(factory.getSchema("BAND::tdc"		));
		hipo::bank	scaler			(factory.getSchema("RUN::scaler"	));
		hipo::event 	readevent;
		
		// Loop over all events in file
		int event_counter = 0;
		gated_charge = 0;
		livetime	= 0;
		while(reader.next()==true){
			// Clear all branches
			starttime 	= 0;
			ePid		= 0;
			eCharge		= 0;
			eStatus		= 0;
			eTime		= 0;
			eBeta 		= 0;
			eChi2pid	= 0;
			E_tot		= 0;
			E_pcal		= 0;
			t_e		= 0;
			dL_e		= 0;
			lU		= 0;
			lV		= 0;
			lW		= 0;
			e_vtx		= 0;
			e_vty		= 0;
			e_vtz		= 0;
			p_e		= 0;
			theta_e		= 0;
			phi_e		= 0;
			q		= 0;
			theta_q		= 0;
			phi_q		= 0;
			nu		= 0;
			Q2		= 0;
			xB		= 0;
			W2		= 0;
			nMult		= 0;
			memset( barID		,0	,sizeof(barID		)	);
			memset( dL_n		,0	,sizeof(dL_n		)	);
			memset( theta_n		,0	,sizeof(theta_n		)	);
			memset( phi_n		,0	,sizeof(phi_n		)	);
			memset( p_n		,0	,sizeof(p_n		)	);
			memset( nTime		,0	,sizeof(nTime		)	);
			memset( nEdep		,0	,sizeof(nEdep		)	);

			// Count events
			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			//if( event_counter > 100000 ) break;
			event_counter++;

			// Load data structure for this event:
			reader.read(readevent);
			readevent.getStructure(event_info);
			readevent.getStructure(particles);
			readevent.getStructure(calorimeter);
			readevent.getStructure(scintillator);
			readevent.getStructure(scaler);
			readevent.getStructure(band_hits);
			readevent.getStructure(band_adcs);
			readevent.getStructure(band_tdcs);
	
			// Currently, REC::Event has uncalibrated livetime / charge, so these will have to work
			livetime 		= 	scaler.getFloat(2,0);
			gated_charge 		= 	scaler.getFloat(0,0) * 0.001; // [microC] -- this seems to be ~10-20% accurate

			// Get integrated charge, livetime and start-time from REC::Event
			if( event_info.getRows() == 0 ) continue;
			getEventInfo( event_info, gated_charge, livetime, starttime );

			// Get electron from particle bank REC::Particle
			TVector3 eVertex, eMomentum;
			getElectronInfo( particles, ePid, eMomentum, eVertex, eTime ,eCharge, eBeta, eChi2pid, eStatus );
			e_vtx = eVertex.X(); e_vty = eVertex.Y(); e_vtz = eVertex.Z(); 
			//	get electron information from scint and calo banks:
			t_e 	= scintillator.getTime(0) - starttime;
			dL_e	= scintillator.getPath(0);
			E_tot 	= calorimeter.getTotE(0);
			E_pcal 	= calorimeter.getPcalE(0);
			lU	= calorimeter.getLU(0);
			lV	= calorimeter.getLV(0);
			lW	= calorimeter.getLW(0);
			//	Do electron PID cuts
			//		none implemented for the moment
			bool ePass = checkElectron( ePid, eMomentum, eVertex, eTime ,eCharge, eBeta, eChi2pid, eStatus , lV , lW , E_tot );
			if( !ePass ) continue;

			// From electron information and beam information, create kinematic variables
			TVector3 beamMomentum(0,0,Ebeam);
			TVector3 qMomentum; qMomentum = beamMomentum - eMomentum;
			p_e		= eMomentum.Mag();
			theta_e		= eMomentum.Theta();
			phi_e		= eMomentum.Phi();
			q		= qMomentum.Mag();
			theta_q		= qMomentum.Theta();
			phi_q		= qMomentum.Phi();
			nu		= Ebeam - p_e;
			Q2		= q*q - nu*nu;
			xB		= Q2 / (2.*mP*nu);
			W2		= mP*mP - Q2 + 2.*nu*mP;
			
			
			// Grab the neutron information:
			TVector3 nMomentum[maxNeutrons], nPath[maxNeutrons];
			getNeutronInfo( band_hits, nMult, barID, nEdep, nTime, nPath , starttime , runNum);
			for( int n = 0 ; n < nMult ; n++ ){
				dL_n[n]		= nPath[n].Mag();
				theta_n[n]	= nPath[n].Theta();
				phi_n[n]	= nPath[n].Phi();
				p_n[n]		= 0.;	// no conversion yet for ToF due to missing calibrations	
			}
			
			// Fill tree based on d(e,e'n)X
			if( nMult == 1 ) outTree->Fill();

		} // end loop over events
	}// end loop over files
	
	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}



void getNeutronInfo( BBand band_hits, int& mult, int id[maxNeutrons], double edep[maxNeutrons],
			double time[maxNeutrons], TVector3 path[maxNeutrons] , double starttime , int thisRun){
	
	if( band_hits.getRows() > maxNeutrons ) return; // not interested in events with more than 1 BAND hit for now
	for( int hit = 0 ; hit < band_hits.getRows() ; hit++ ){
		if( band_hits.getLayer(hit) == 6 ) continue; // not interested in a veto hit

		id[hit]		= band_hits.getBarKey(hit);
		edep[hit]	= sqrt( band_hits.getAdcLcorr(hit) * band_hits.getAdcRcorr(hit) );
		double tof_fix = (band_hits.getMeantimeTdc(hit) - starttime ) - TDC_GLOBSHIFT[band_hits.getBarKey(hit)] - TDC_RUNBYRUNSHIFT[thisRun];
		time[hit]	= tof_fix;
		path[hit].SetXYZ(	band_hits.getX(hit), band_hits.getY(hit), band_hits.getZ(hit) 	);

		mult++;
	}
	
}

int getRunNumber( string filename ){
	string parsed = filename.substr( filename.find("inc") );
	//string parsed = filename.substr( filename.find("_clas") );
	//string parsed = filename.substr( filename.find("_band") );
	string moreparse = parsed.substr(6,8);
	cout << filename << " " << parsed << "\n";
	cout << moreparse << " " << stoi(moreparse) << "\n\n";
        return stoi(moreparse);
}

void getEventInfo( BEvent eventInfo, double &integrated_charge, double &livetime, double &starttime ){
	if( eventInfo.getRows() != 1 ){ 
		cerr << "getEventInfo::NotImplementedFunction\n"; 
		exit(-1); 
	}
	//integrated_charge       += (double)eventInfo.getBCG(0); 	// not calibrated currently
	//livetime 		= (double)eventInfo.getLT(0);		// not calibrated currently
	starttime		= (double)eventInfo.getSTTime(0);
	return;
}
void getElectronInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex,
			double& time, int& charge, double& beta, double& chi2pid, int& status ){
	pid 		= particles.getPid(0);
	momentum 	= particles.getV3P(0);
	vertex		= particles.getV3v(0);
	time		= particles.getVt(0);
	charge		= particles.getCharge(0);
	beta		= particles.getBeta(0);
	chi2pid		= particles.getChi2pid(0);
	status		= particles.getStatus(0);
	return;
}
void getProtonInfo( BParticle particles, double pid[maxProtons], TVector3 momentum[maxProtons], TVector3 vertex[maxProtons],
			double time[maxProtons], double charge[maxProtons], double beta[maxProtons], double chi2pid[maxProtons], double status[maxProtons] , int& multiplicity ){
	// Takes the first proton in the bank
	multiplicity = 0;
	for( int row = 1 ; row < particles.getRows() ; row++ ){ // start after electron information
		pid[multiplicity] 		= particles.getPid(row);
		charge[multiplicity]		= particles.getCharge(row);
		if( charge[multiplicity] == 1 ){
			momentum 	[multiplicity]	= particles.getV3P(row);
			vertex		[multiplicity]	= particles.getV3v(row);
			time		[multiplicity]	= particles.getVt(row);
			beta		[multiplicity]	= particles.getBeta(row);
			chi2pid		[multiplicity]	= particles.getChi2pid(row);
			status		[multiplicity]	= particles.getStatus(row);
			multiplicity ++;
		}
	}
}
bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status,
			double lV, double lW , double E_tot){
			// Anymore fiducial cuts that are needed for the electron:
			// 	E/p cut
			// 	vertex cut
			// 	minimum momentum  cut / maximum  momentum  cut
			// 	electron ToF cut
			// 	minimum W cut
			//	chi2 cut
			//	pcal energy cut?
			//

	if( pid != 11 || charge != -1 ) return false;
	//if( lV < 2 || lW < 2 ) return false;
	//if( momentum.Mag() < 1 || momentum.Mag() > 4.2 ) return false;
	//if( chi2pid == 0 || chi2pid > 1 ) return false;
	//if( E_tot/momentum.Mag() > 0.4 || E_tot/momentum.Mag() < 0.1 ) return false;
	//if( vertex.X() < -2 || vertex.X() > 2) return false;
	//if( vertex.Y() < -2 || vertex.Y() > 2) return false;
	//if( vertex.Z() < -7 || vertex.Z() > 2) return false;
	//if( time < 15 ) return false;
	//if( momentum.Mag() < 2 || momentum.Mag() > 10.6 ) return false;
	//if( E_tot / momentum.Mag() < 0.15 || E_tot / momentum.Mag() > 0.3 ) return false;

	return true;
}
bool checkProton( int pid, TVector3 momentum, TVector3 del_vertex, double time, int charge, double beta, double chi2pid, int status, int mult ){
	//if( momentum.Mag() == 0 || momentum.Mag() > 4.2 ) return false;
	//if( beta <= 0 || beta > 1 ) return false;
	//if( mult != 1 ) return false;
	//if( momentum.Theta() == 0 ) return false;
	//if( beta <= 0 ) return false;
	//if( del_vertex.Z() < 0 || del_vertex.Z() > 5 ) return false;
	//if( chi2pid < -3 || chi2pid > 3 ) return false;
	
	
	return true;
}

bool pointsToBand(double theta,double phi,double z_m){
	//double z = z_m*100; // from m to cm
	double z = z_m;

	// Numbers taken from band/src/main/java/org/jlab/rec/band/constants/Parameters.java
	double thickness  = 7.2;                                // thickness of each bar (cm)
	double layerGap[] = {7.94, 7.62, 7.94, 7.62, 7.3};      // gap between center of neighbouring layers (cm), 1-2, 2-3, 3-4, 4-5, 5-6

	// Numbers taken from clas-band-calib/bin/src/org/clas/fcmon/band/BANDConstants.java
	double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};

	// Distance from ideal target to upstream end of BAND
	// (from BAND survey report, 02/18/2019)
	double zUpst = (-302.69-302.69-302.57-302.64)/4.; // [cm]

	// Distance from ideal target to downstream end of layer 5
	double zDown = (zUpst + 5*thickness) - z_m;

	double rho   = zDown/cos(theta);
	double xDown = rho*sin(theta)*cos(phi);
	double yDown = rho*sin(theta)*sin(phi);

	double globalX = (-240.5-240.5+241.0+243.7)/4.; // [cm] --> Not using this yet (need to make sure we have the right coordinate system)
	double globalY = (-211.0+228.1-210.6+228.1)/4.; // [cm]

	// Sector boundaries
	double topSec1  = globalY + 13*thickness;
	double topSec2  = globalY + 10*thickness;
	double topSec34 = globalY +  3*thickness;
	double topSec5  = globalY -  3*thickness;
	double downSec5 = globalY -  5*thickness;

	if( yDown >= topSec1 || yDown <= downSec5 ) return 0;

	if(             (yDown < topSec1  && yDown >= topSec2  && fabs(xDown) < bandlen[0]/2. )||
			(yDown < topSec2  && yDown >= topSec34 && fabs(xDown) < bandlen[1]/2. )||
			(yDown < topSec34 && yDown >= topSec5  && fabs(xDown) < bandlen[1]/2. && fabs(xDown) > bandlen[1]/2.-bandlen[2])||
			(yDown < topSec5  && yDown >= downSec5 && fabs(xDown) < bandlen[4]/2. )
	  )
		return 1;

	return 0;
}

void LoadGlobalShift(){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open("../include/global_offset_fadc-10082019.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		FADC_GLOBSHIFT[barId] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
	f.open("../include/global_offset_tdc_1stIter-04132020.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		TDC_GLOBSHIFT[barId] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}
void LoadRunByRunShift(){
	ifstream f;
	int runnum;
	double pol0, height, mean, sig, temp;

	f.open("../include/runByrun_offset_fadc-10082019.txt");
	while(!f.eof()){
		f >> runnum;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		FADC_RUNBYRUNSHIFT[runnum] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
	f.open("../include/runByrun_offset_tdc_firstIter-04132020.txt");
	while(!f.eof()){
		f >> runnum;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		TDC_RUNBYRUNSHIFT[runnum] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}
