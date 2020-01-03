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

int getRunNumber( string filename );
double getBeamEnergy( int runNum );
void getEventInfo( BEvent eventInfo, double &integrated_charge, double &livetime, double &starttime );
void getElectronInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex, 
			double& time, int& charge, double& beta, double& chi2pid, int& status );
bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status,
			double lV, double lW , double E_tot);
void getProtonInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex,
			double& time, int& charge, double& beta, double& chi2pid, int& status , int& multiplicity );
bool checkProton( int pid, TVector3 momentum, TVector3 del_vertex, double time, int charge, double beta, double chi2pid, int status, int mult );
void getNeutronInfo(BBand band_hits, int& nHits, vector<int>& barKey, vector<int>& layer, vector<double>& meanADC, 
			vector<double>& meanTimeFADC, vector<TVector3>& hitPos);
bool checkNeutron(int nHits, vector<int> barKey, vector<int> layer, vector<double> meanADC, vector<double> meanTimeFADC, 
			vector<TVector3> hitPos, double starttime, double* FADC_GLOB_SHIFT, double &ToF , int & nMult , 
			double &Edep, int &thisBar, double &dL, double &theta, double &phi ) ;
bool pointsToBand(double theta,double phi,double z_m);
void LoadGlobalShift(double* FADC_GLOB_SHIFT);

int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 3 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [inputFile]\n\n";
		return -1;
	}

	// Create output tree
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("skim","CLAS and BAND Physics");
	double gated_charge	= 0;
	int ePid		= 0;
	int eCharge		= 0;
	int eStatus		= 0;
	double eTime		= 0;
	double eBeta 		= 0;
	double eChi2pid		= 0;
	double E_tot		= 0;
	double E_pcal		= 0;
	double t_e		= 0;
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
	double pPid		= 0;
	double pCharge		= 0;
	double pStatus		= 0;
	double pTime		= 0;
	double pBeta		= 0;
	double pChi2pid		= 0;
	double p_vtx		= 0;
	double p_vty		= 0;
	double p_vtz		= 0;
	double p_p		= 0;
	double theta_p		= 0;
	double phi_p		= 0;
	double protons		= 0;
	double p_miss		= 0;
	double m_miss		= 0;
	double theta_miss	= 0;
        double phi_miss		= 0;
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("ePid"		,&ePid			);
	outTree->Branch("eCharge"	,&eCharge		);
	outTree->Branch("eStatus"	,&eStatus		);
	outTree->Branch("eTime"	,&eTime				);
	outTree->Branch("eBeta"	,&eBeta 			);
	outTree->Branch("eChi2pid"	,&eChi2pid			);
	outTree->Branch("E_tot"	,&E_tot				);
	outTree->Branch("E_pcal"	,&E_pcal			);
	outTree->Branch("t_e"	,&t_e				);
	outTree->Branch("lU"	,&lU				);
	outTree->Branch("lV"	,&lV				);
	outTree->Branch("lW"	,&lW				);
	outTree->Branch("e_vtx"	,&e_vtx				);
	outTree->Branch("e_vty"	,&e_vty				);
	outTree->Branch("e_vtz"	,&e_vtz				);
	outTree->Branch("p_e"	,&p_e				);
	outTree->Branch("theta_e"	,&theta_e			);
	outTree->Branch("phi_e"	,&phi_e				);
	outTree->Branch("q"	,&q				);
	outTree->Branch("theta_q"	,&theta_q			);
	outTree->Branch("phi_q"	,&phi_q				);
	outTree->Branch("nu"	,&nu				);
	outTree->Branch("Q2"	,&Q2				);
	outTree->Branch("xB"	,&xB				);
	outTree->Branch("W2"	,&W2				);
	outTree->Branch("pPid"	,&pPid				);
	outTree->Branch("pCharge"	,&pCharge			);
	outTree->Branch("pStatus"	,&pStatus			);
	outTree->Branch("pTime"	,&pTime				);
	outTree->Branch("pBeta"	,&pBeta				);
	outTree->Branch("pChi2pid"	,&pChi2pid			);
	outTree->Branch("p_vtx"	,&p_vtx				);
	outTree->Branch("p_vty"	,&p_vty				);
	outTree->Branch("p_vtz"	,&p_vtz				);
	outTree->Branch("p_p"	,&p_p				);
	outTree->Branch("theta_p"	,&theta_p			);
	outTree->Branch("phi_p"	,&phi_p				);
	outTree->Branch("protons"	,&protons			);
	outTree->Branch("p_miss"	,&p_miss			);
	outTree->Branch("m_miss"	,&m_miss			);
	outTree->Branch("theta_miss"	,&theta_miss			);
	outTree->Branch("phi_miss"	,&phi_miss			);
	double ToF		= 0;
	double Edep		= 0;
	int nMult		= 0;
	int thisBar		= 0;
	int pointing		= 0;
	double dL		= 0;
	double theta_n		= 0;
	double phi_n		= 0;
	outTree->Branch("ToF"		,&ToF			);
	outTree->Branch("Edep"		,&Edep			);
	outTree->Branch("nMult"		,&nMult			);
	outTree->Branch("thisBar"	,&thisBar		);
	outTree->Branch("pointing"	,&pointing		);
	outTree->Branch("dL"		,&dL			);
	outTree->Branch("theta_n"	,&theta_n		);
	outTree->Branch("phi_n"		,&phi_n			);

	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

	// Load FADC offsets
	double FADC_GLOB_SHIFT[600];
	LoadGlobalShift(FADC_GLOB_SHIFT);

	// Load input file
	for( int i = 2 ; i < argc ; i++ ){
		// Using run number of current file, grab the beam energy from RCDB
		int runNum = getRunNumber(argv[i]);
		auto cnd = connection.GetCondition(runNum, "beam_energy");
		double Ebeam = cnd->ToDouble() / 1000.; // [GeV]

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
		while(reader.next()==true){
			// electron clear
			ePid		= 0;
			eCharge		= 0;
			eStatus		= 0;
			eTime		= 0;
			eBeta 		= 0;
			eChi2pid	= 0;
			E_tot		= 0;
			E_pcal		= 0;
			t_e		= 0;
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
			pPid		= 0;
			pCharge		= 0;
			pStatus		= 0;
			pTime		= 0;
			pBeta		= 0;
			pChi2pid	= 0;
			p_vtx		= 0;
			p_vty		= 0;
			p_vtz		= 0;
			p_p		= 0;
			theta_p		= 0;
			phi_p		= 0;
			protons		= 0;
			p_miss		= 0;
			m_miss		= 0;
			theta_miss	= 0;
			phi_miss	= 0;
			ToF		= 0;
			Edep		= 0;
			nMult		= 0;
			thisBar		= 0;
			pointing	= 0;
			dL		= 0;
			theta_n		= 0;
			phi_n		= 0;


			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			//if(event_counter >= 1000000) break;
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
			double livetime 	= 	scaler.getFloat(2,0);
			gated_charge 		= 	scaler.getFloat(0,0) * 0.001; // [microC] -- this seems to be ~10-20% accurate

			// Get integrated charge, livetime and start-time from REC::Event
			double starttime = 0;
			if( event_info.getRows() == 0 ) continue;
			getEventInfo( event_info, gated_charge, livetime, starttime );

			// Get electron from particle bank REC::Particle
			TVector3 eVertex, eMomentum;
			getElectronInfo( particles, ePid, eMomentum, eVertex, eTime ,eCharge, eBeta, eChi2pid, eStatus );
			//	get electron information from scint and calo banks:
			t_e 	= scintillator.getTime(0) - starttime;
			E_tot 	= calorimeter.getTotE(0);
			E_pcal 	= calorimeter.getPcalE(0);
			lU	= calorimeter.getLU(0);
			lV	= calorimeter.getLV(0);
			lW	= calorimeter.getLW(0);
			//	Do electron PID cuts
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
			
			// Grab the proton information:
			TVector3 pVertex, pMomentum;
			int pPid = 0, pCharge = 0, pStatus = 0, pMult = 0;
			double pTime = 0, pBeta = 0, pChi2pid = 0;
			getProtonInfo( particles, pPid, pMomentum, pVertex, pTime ,pCharge, pBeta, pChi2pid, pStatus, pMult );
			//	do proton PID cuts		
			bool pPass = checkProton( pPid, pMomentum, pVertex-eVertex, pTime ,pCharge, pBeta, pChi2pid, pStatus , pMult );
			if( !pPass ) continue;
			// Create more kinematic variables
			p_p		= pMomentum.Mag();
			theta_p		= pMomentum.Theta();
			phi_p		= pMomentum.Phi();
			protons		= pMult;
			double E_p	= sqrt(p_p*p_p + mP*mP);
			// 	assuming d e -> e' p 
			TVector3 missMomentum; missMomentum = pMomentum - qMomentum;
			p_miss		= missMomentum.Mag();
			m_miss		= mP - mD + sqrt( pow(nu+mD-E_p,2) + pow(p_miss,2) );
			theta_miss	= missMomentum.Theta();
			phi_miss	= missMomentum.Phi();

			
			// Get the neutron information
			int nHits;
			vector<int> barKey;
			vector<int> layer;
			vector<double> meanADC;
			vector<double> meanTimeFADC;
			vector<TVector3> hitPos;
				
			getNeutronInfo(band_hits, nHits, barKey, layer, meanADC, meanTimeFADC, hitPos);
			
			bool nPass = checkNeutron(nHits, barKey, layer, meanADC, meanTimeFADC, hitPos, 
					starttime, FADC_GLOB_SHIFT, ToF, nMult, Edep, thisBar , dL, theta_n, phi_n );
			//cout << ToF << " " << Edep << " " << dL << " " << theta_n << " " << phi_n << "\n\n";
			
			TVector3 nMomentum = -missMomentum;
			pointing = pointsToBand(nMomentum.Theta() , nMomentum.Phi() , eVertex.Z() );
			//if (!nPass) continue;
		
			//band_adcs.show();
			//band_tdcs.show();


			// put more neutron stuff here 


			// Fill tree to do any more plots on the fly
			outTree->Fill();

		} // end loop over events
		cout << "Total charge collected in file: " << gated_charge << " [microC]\n";
		cout << "Total number of events (e,e'pn): " << outTree->GetEntries() << "\n";
	}// end loop over files
	

	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}



void getNeutronInfo(BBand band_hits, int& nHits, vector<int>& barKey, vector<int>& layer, vector<double>& meanADC, vector<double>& meanTimeFADC, vector<TVector3>& hitPos) {

	nHits = band_hits.getRows();

	for (int hit = 0; hit < nHits; hit++) {
	
		barKey.push_back(band_hits.getBarKey(hit));
		layer.push_back(band_hits.getLayer(hit));
		double ADCL = band_hits.getAdcLcorr(hit);
		double ADCR = band_hits.getAdcRcorr(hit);
		double ADCLR = sqrt(ADCL*ADCR);
		meanADC.push_back(ADCLR);
		meanTimeFADC.push_back(band_hits.getMeantimeFadc(hit));
		double hitx = band_hits.getX(hit);
		double hity = band_hits.getY(hit);
		double hitz = band_hits.getZ(hit);
		hitPos.push_back(TVector3(hitx, hity, hitz));	
		
	}
	
}

int getRunNumber( string filename ){
	//string parsed = filename.substr( filename.find("inc") );
	//string parsed = filename.substr( filename.find("_clas") );
	string parsed = filename.substr( filename.find("_band") );
	string moreparse = parsed.substr(6,8);
	cout << filename << " " << parsed << "\n";
	cout << moreparse << " " << stoi(moreparse) << "\n\n";
        return stoi(moreparse);
}

double getBeamEnergy( int runNum ){
        double thisEn = 0;

        if( runNum <= 6399 ) thisEn = 10.6;
        else{ thisEn = 10.2; }
        if( runNum == 6523 || runNum == 6524 || runNum == 6525 ) thisEn = 10.;
	

        return thisEn;
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
void getProtonInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex,
			double& time, int& charge, double& beta, double& chi2pid, int& status , int& multiplicity ){
	// Takes the first proton in the bank
	multiplicity = 0;
	for( int row = 0 ; row < particles.getRows() ; row++ ){
		pid 		= particles.getPid(row);
		charge		= particles.getCharge(row);
		if( pid == 2212 && charge == 1 ){
			if( multiplicity == 0 ){ // if this is the first proton, save the information
				momentum 	= particles.getV3P(row);
				vertex		= particles.getV3v(row);
				time		= particles.getVt(row);
				beta		= particles.getBeta(row);
				chi2pid		= particles.getChi2pid(row);
				status		= particles.getStatus(row);
			}
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

bool checkNeutron(int nHits, vector<int> barKey, vector<int> layer, vector<double> meanADC, vector<double> meanTimeFADC, 
			vector<TVector3> hitPos, double starttime, double* FADC_GLOB_SHIFT, double &ToF , int & nMult , 
			double &Edep, int &thisBar, double &dL, double &theta, double &phi ) {

	for (int hit = 0; hit < nHits; hit++) {
		if( nMult == 0 ){
			ToF = meanTimeFADC[hit] - starttime - FADC_GLOB_SHIFT[barKey[hit]];
			Edep = meanADC[hit];
			thisBar = barKey[hit];
			dL = hitPos[hit].Mag();
			theta = hitPos[hit].Theta();
			phi = hitPos[hit].Phi();
			
			//cout << ToF << " " << Edep << " " << dL << " " << theta << " " << phi << "\n";
		}
		nMult++;

	}

	return true;
}


void LoadGlobalShift(double* FADC_GLOB_SHIFT){

	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open("../include/global_offset_fadc.txt");

	while(!f.eof()){

		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		FADC_GLOB_SHIFT[barId] = mean;
		f >> sig;
		f >> temp;
		f >> temp;

	}

	f.close();

}



bool pointsToBand(double theta,double phi,double z_m){
	double z = z_m*100; // from m to cm

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
