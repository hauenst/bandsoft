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

const int maxPions	= 100;

int getRunNumber( string filename );
void getEventInfo( BEvent eventInfo, double &integrated_charge, double &livetime, double &starttime );
void getElectronInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex, 
			double& time, int& charge, double& beta, double& chi2pid, int& status );
bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status,
			double lV, double lW , double E_tot);
void getProtonInfo( BParticle particles, double pid[maxPions], TVector3 momentum[maxPions], TVector3 vertex[maxPions],
			double time[maxPions], double charge[maxPions], double beta[maxPions], double chi2pid[maxPions], double status[maxPions] , int& multiplicity );
bool checkProton( int pid, TVector3 momentum, TVector3 del_vertex, double time, int charge, double beta, double chi2pid, int status, int mult );

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
	// 	Proton info:
	int pMult		= 0;
	double pPid		[maxPions]= {0.};
	double pCharge		[maxPions]= {0.};
	double pStatus		[maxPions]= {0.};
	double pTime		[maxPions]= {0.};
	double pBeta		[maxPions]= {0.};
	double pChi2pid		[maxPions]= {0.};
	double p_vtx		[maxPions]= {0.};
	double p_vty		[maxPions]= {0.};
	double p_vtz		[maxPions]= {0.};
	double p_p		[maxPions]= {0.};
	double theta_p		[maxPions]= {0.};
	double phi_p		[maxPions]= {0.};
	double p_miss		[maxPions]= {0.};
	double m_miss		[maxPions]= {0.};
	double theta_miss	[maxPions]= {0.};
        double phi_miss		[maxPions]= {0.};
	double theta_pq		[maxPions]= {0.};
	double Wp		[maxPions]= {0.};
	double z		[maxPions]= {0.};
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
	outTree->Branch("pMult"		,&pMult			);
	outTree->Branch("pPid"		,pPid			,"pPid[pMult]/D"	);
	outTree->Branch("pCharge"	,pCharge		,"pCharge[pMult]/D"	);
	outTree->Branch("pStatus"	,pStatus		,"pStatus[pMult]/D"	);
	outTree->Branch("pTime"		,pTime			,"pTime[pMult]/D"	);
	outTree->Branch("pBeta"		,pBeta			,"pBeta[pMult]/D"	);
	outTree->Branch("pChi2pid"	,pChi2pid		,"pChi2pid[pMult]/D"	);
	outTree->Branch("p_vtx"		,p_vtx			,"p_vtx[pMult]/D"	);
	outTree->Branch("p_vty"		,p_vty			,"p_vty[pMult]/D"	);
	outTree->Branch("p_vtz"		,p_vtz			,"p_vtz[pMult]/D"	);
	outTree->Branch("p_p"		,p_p			,"p_p[pMult]/D"		);
	outTree->Branch("theta_p"	,theta_p		,"theta_p[pMult]/D"	);
	outTree->Branch("phi_p"		,phi_p			,"phi_p[pMult]/D"	);
	outTree->Branch("p_miss"	,p_miss			,"p_miss[pMult]/D"	);
	outTree->Branch("m_miss"	,m_miss			,"m_miss[pMult]/D"	);
	outTree->Branch("theta_miss"	,theta_miss		,"theta_miss[pMult]/D"	);
	outTree->Branch("phi_miss"	,phi_miss		,"phi_miss[pMult]/D"	);
	outTree->Branch("theta_pq"	,theta_pq		,"theta_pq[pMult]/D"	);
	outTree->Branch("Wp"		,Wp			,"Wp[pMult]/D"		);
	outTree->Branch("z"		,z			,"z[pMult]/D"		);
	

	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

	// Load input file
	for( int i = 2 ; i < argc ; i++ ){
		// Using run number of current file, grab the beam energy from RCDB
		int runNum = getRunNumber(argv[i]);
		auto cnd = connection.GetCondition(runNum, "beam_energy");
		Ebeam = cnd->ToDouble() / 1000.; // [GeV] -- conversion factor due to miscalibration in RCDB

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
		hipo::bank	scaler			(factory.getSchema("RUN::scaler"	));
		hipo::event 	readevent;
		
		// Loop over all events in file
		int event_counter = 0;
		gated_charge = 0;
		livetime	= 0;
		while(reader.next()==true){
			if( event_counter > 1000000 ) break;
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
			pMult		= 0;
			memset(	pPid		,0	,sizeof(pPid		)	);
			memset(	pCharge		,0	,sizeof(pCharge		)	);
			memset(	pStatus		,0	,sizeof(pStatus		)	);
			memset(	pTime		,0	,sizeof(pTime		)	);
			memset(	pBeta		,0	,sizeof(pBeta		)	);
			memset(	pChi2pid	,0	,sizeof(pChi2pid	)	);
			memset(	p_vtx		,0	,sizeof(p_vtx		)	);
			memset(	p_vty		,0	,sizeof(p_vty		)	);
			memset(	p_vtz		,0	,sizeof(p_vtz		)	);
			memset(	p_p		,0	,sizeof(p_p		)	);
			memset(	theta_p		,0	,sizeof(theta_p		)	);
			memset(	phi_p		,0	,sizeof(phi_p		)	);
			memset(	p_miss		,0	,sizeof(p_miss		)	);
			memset(	m_miss		,0	,sizeof(m_miss		)	);
			memset(	theta_miss	,0	,sizeof(theta_miss	)	);
			memset(	phi_miss	,0	,sizeof(phi_miss	)	);
			memset( theta_pq	,0	,sizeof(theta_pq	)	);
			memset( Wp		,0	,sizeof(Wp		)	);
			memset( z		,0	,sizeof(z		)	);

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
			
			
			// Grab the pion information:
			TVector3 pVertex[maxPions], pMomentum[maxPions];
			getProtonInfo( particles, pPid, pMomentum, pVertex, pTime ,pCharge, pBeta, pChi2pid, pStatus, pMult );
			for( int p = 0 ; p < pMult ; p++ ){
				p_vtx[p]	= pVertex[p].X();
				p_vty[p]	= pVertex[p].Y();
				p_vtz[p]	= pVertex[p].Z();
				p_p[p]		= pMomentum[p].Mag();
				theta_p[p]	= pMomentum[p].Theta();
				phi_p[p]	= pMomentum[p].Phi();

				TVector3 missMomentum; missMomentum = -(pMomentum[p] - qMomentum);
				p_miss[p]	= missMomentum.Mag();
				theta_miss[p]	= missMomentum.Theta();
				phi_miss[p]	= missMomentum.Phi();
				theta_pq[p]	= qMomentum.Angle(pMomentum[p]);

				double E_p = sqrt( p_p[p]*p_p[p] + mPi*mPi );
				m_miss[p] 	= sqrt( pow( nu + mD - E_p , 2 ) - ( q*q + p_p[p]*p_p[p] - 2*q*p_p[p]*cos(theta_pq[p]) ) );
				if( m_miss[p] != m_miss[p] ) m_miss[p] = 0.;

				double W_primeSq = mD*mD - Q2 + mPi*mPi + 2.*mD*(nu-E_p) - 2.*nu*E_p + 2.*q*p_p[p]*cos(theta_pq[p]);
				Wp[p] = sqrt(W_primeSq);
				if( Wp[p] != Wp[p] ) Wp[p] = 0.;

				z[p] = E_p / nu;

			}

			// Fill tree to do any more plots on the fly
			outTree->Fill();

		} // end loop over events
	}// end loop over files
	

	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
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
void getProtonInfo( BParticle particles, double pid[maxPions], TVector3 momentum[maxPions], TVector3 vertex[maxPions],
			double time[maxPions], double charge[maxPions], double beta[maxPions], double chi2pid[maxPions], double status[maxPions] , int& multiplicity ){
	// Takes the first pion in the bank
	multiplicity = 0;
	for( int row = 1 ; row < particles.getRows() ; row++ ){ // start after electron information
		pid[multiplicity] 		= particles.getPid(row);
		charge[multiplicity]		= particles.getCharge(row);
		if( charge[multiplicity] == 1 && pid[multiplicity] == 211 ){
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


