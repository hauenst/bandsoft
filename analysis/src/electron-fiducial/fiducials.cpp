#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TGraph.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TCutG.h"

using namespace std;

int main( int argc, char** argv){
	if( argc != 3 ){
		cerr << "Incorrect number of arguments. Please instead use:\n"
			<< "\t./code [InputInclusiveFile] [OutputName - used for pdf and root file]\n";
		return -1;
	}

	/////////////////////////////////////////////////////////////////////
	// Setup the input root file
	TFile * inFile = new TFile(argv[1]);
	TTree * inTree = (TTree*)inFile->Get("skim");
	double Ebeam		= 0;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	int    ePid		= 0;
	int    eCharge		= 0;
	int    eStatus		= 0;
	double eTime		= 0;
	double eBeta 		= 0;
	double eChi2pid		= 0;
	double E_tot		= 0;
	double E_pcal		= 0;
	int    eSector		= 0;
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
	double weight		= 0;
	inTree->SetBranchAddress("Ebeam"		,&Ebeam			);
	inTree->SetBranchAddress("gated_charge"		,&gated_charge		);
	inTree->SetBranchAddress("livetime"		,&livetime		);
	inTree->SetBranchAddress("starttime"		,&starttime		);
	inTree->SetBranchAddress("ePid"			,&ePid			);
	inTree->SetBranchAddress("eCharge"		,&eCharge		);
	inTree->SetBranchAddress("eStatus"		,&eStatus		);
	inTree->SetBranchAddress("eTime"		,&eTime			);
	inTree->SetBranchAddress("eBeta"		,&eBeta 		);
	inTree->SetBranchAddress("eChi2pid"		,&eChi2pid		);
	inTree->SetBranchAddress("E_tot"		,&E_tot			);
	inTree->SetBranchAddress("E_pcal"		,&E_pcal		);
	inTree->SetBranchAddress("eSector"		,&eSector		);
	inTree->SetBranchAddress("t_e"			,&t_e			);
	inTree->SetBranchAddress("dL_e"			,&dL_e			);
	inTree->SetBranchAddress("lU"			,&lU			);
	inTree->SetBranchAddress("lV"			,&lV			);
	inTree->SetBranchAddress("lW"			,&lW			);
	inTree->SetBranchAddress("e_vtx"		,&e_vtx			);
	inTree->SetBranchAddress("e_vty"		,&e_vty			);
	inTree->SetBranchAddress("e_vtz"		,&e_vtz			);
	inTree->SetBranchAddress("p_e"			,&p_e			);
	inTree->SetBranchAddress("theta_e"		,&theta_e		);
	inTree->SetBranchAddress("phi_e"		,&phi_e			);
	inTree->SetBranchAddress("q"			,&q			);
	inTree->SetBranchAddress("theta_q"		,&theta_q		);
	inTree->SetBranchAddress("phi_q"		,&phi_q			);
	inTree->SetBranchAddress("nu"			,&nu			);
	inTree->SetBranchAddress("Q2"			,&Q2			);
	inTree->SetBranchAddress("xB"			,&xB			);
	inTree->SetBranchAddress("W2"			,&W2			);
	inTree->SetBranchAddress("weight"		,&weight		);

	/////////////////////////////////////////////////////////////////////
	// Setup the output root file and histograms
	TString outname = Form("%s.root",argv[2]);
	TFile * outFile = new TFile(outname,"RECREATE");
	TH2D ** h2_theta_phi = new TH2D*[6];
	for( int sector = 0 ; sector < 6 ; sector++ )
		h2_theta_phi[sector] = new TH2D(Form("h2_theta_phi_%i",sector),"",420,-220,200,160,0,40);

	
	/////////////////////////////////////////////////////////////////////////////
	// Let's take the (theta,phi) coverage with some basic electron PID cuts
	// and then fit for each sector. Once we do that, take the (p,mom) curve
	// and verify that it follows some Q2 dependence minimum based on CLAS accetpance
	for( int event = 0 ; event < inTree->GetEntries() ; event++ ){
		if( event%10000==0) cout << "\tevent: " << event << "\n";
		//if( event > 10000000 ) break;
		// Clear out all the branch variables
		Ebeam		= 0;
		gated_charge	= 0;
		livetime	= 0;
		starttime	= 0;
		ePid		= 0;
		eCharge		= 0;
		eStatus		= 0;
		eTime		= 0;
		eBeta 		= 0;
		eChi2pid	= 0;
		E_tot		= 0;
		E_pcal		= 0;
		eSector		= 0;
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
		weight		= 0;
		// Grab the event from the tree
		inTree->GetEvent(event);
		
		// Ask if we pass some basic PID cuts:
		if(	E_tot/p_e 	< 0.17	||
			E_tot/p_e 	> 0.27	||
			E_tot		< 0.25  ||
			e_vtz 		< -6	||
			e_vtz		> 0	||
			p_e 		> Ebeam 	) continue;
		// Fill the theta-phi coverage for each sector:
		if( eSector == 4 ){
			if( phi_e > 0 ) phi_e -= 2*M_PI;
		}
		h2_theta_phi[eSector-1]	-> Fill( (phi_e)*180./M_PI, theta_e*180./M_PI	);
		
	}
	
	/////////////////////////////////////////////////////////////////////////////
	// For each of these histograms, take slices in theta (y-axis) and try to find
	// the edges in phi based on First/LastBinAbove:
	double 	sliceWidth 	= 0.5; // degrees
	int 	nTotBins 	= h2_theta_phi[0]->GetNbinsY();
	double  histXmin	= h2_theta_phi[0]->GetYaxis()->GetXmin();
	double 	histRange 	= h2_theta_phi[0]->GetYaxis()->GetXmax() - h2_theta_phi[0]->GetYaxis()->GetXmin();
	int 	numSlices 	= histRange/sliceWidth;
	
	outFile->cd();
	TGraph ** fiducialsThetaSlice = new TGraph*[6];
	for( int sector = 0 ; sector < 6 ; sector++ ){
	
		std::vector<double> left_thetas;
		std::vector<double> left_phis;
		std::vector<double> right_thetas;
		std::vector<double> right_phis;
		for( int slice = 0 ; slice < numSlices ; slice++ ){
		
			int firstBin 	= h2_theta_phi[sector]->GetYaxis()->FindBin(histXmin + slice*sliceWidth);
			int lastBin 	= h2_theta_phi[sector]->GetYaxis()->FindBin(histXmin + (slice+1)*sliceWidth);
			TH1D* proj_x 	= (TH1D*)h2_theta_phi[sector]->
						ProjectionX(Form("h2_theta_phi_%i_%i",sector,slice),firstBin,lastBin);

			if( proj_x->GetEntries() < 1000. ) continue;
			
			// Just try to do a simple maximum*0.1 
			double threshold = proj_x->GetMaximum()	* 0.1;
			double low_phi = proj_x->GetBinCenter(proj_x->FindFirstBinAbove(threshold));
			double high_phi = proj_x->GetBinCenter(proj_x->FindLastBinAbove(threshold));
			double mean_theta = histXmin + (slice+0.5)*sliceWidth;
	
			// Store to vectors
			left_thetas.push_back(mean_theta);
			left_phis.push_back(low_phi);
			right_thetas.push_back(mean_theta);
			right_phis.push_back(high_phi);
	
			// Save the projection
			//proj_x->Write();

		}
		//h2_theta_phi[sector]->Write();
		
		// Combine the left and right side of the fiducial into one 
		// vector, but we need to flip the order of the left side so that it will
		// be ordered in increasing phi
		std::vector<double> full_thetas;
		std::vector<double> full_phis;
		//	flip the left side:
		std::reverse(	left_phis.begin(), 	left_phis.end());
		std::reverse(	left_thetas.begin(), 	left_thetas.end());
		// 	Add the left thetas to the full thetas
		full_thetas.insert( full_thetas.end() , left_thetas.begin(),      left_thetas.end() );
		full_thetas.insert( full_thetas.end() , right_thetas.begin(),     right_thetas.end() );
		full_phis.insert( full_phis.end() , left_phis.begin(),		left_phis.end()	);
		full_phis.insert( full_phis.end() , right_phis.begin(), 	right_phis.end() );

		fiducialsThetaSlice[sector] = new TGraph( full_thetas.size(), &full_phis[0], &full_thetas[0]);
		fiducialsThetaSlice[sector]->SetName(Form("fiducialsThetaSlice_%i",sector));
		//fiducialsThetaSlice[sector]->Write();


		// Save them as a TCutG:
		TCutG *thisCut = new TCutG(Form("eCut_%i",sector),full_thetas.size()+1);
		thisCut->SetVarX("phi_e");
		thisCut->SetVarY("theta_e");
		for( int i = 0 ; i < full_thetas.size() ; i++ ){
			thisCut->SetPoint(i,full_phis[i]*M_PI/180.,full_thetas[i]*M_PI/180.);
		}
		thisCut->SetPoint(full_thetas.size(),full_phis[0]*M_PI/180.,full_thetas[0]*M_PI/180.);
		thisCut->Write();
	}
	
	/////////////////////////////////////////////////////////////////////////////
	// Now do the same but take slices in phi (x-axis)
	sliceWidth 	= 1; // degreesa
	nTotBins 	= h2_theta_phi[0]->GetNbinsX();
	histXmin	= h2_theta_phi[0]->GetXaxis()->GetXmin();
	histRange 	= h2_theta_phi[0]->GetXaxis()->GetXmax() - h2_theta_phi[0]->GetXaxis()->GetXmin();
	numSlices 	= histRange/sliceWidth;
	TGraph ** fiducialsPhiSlice = new TGraph*[6];
	for( int sector = 0 ; sector < 6 ; sector++ ){

		std::vector<double> thetas;
		std::vector<double> phis;
		for( int slice = 0 ; slice < numSlices ; slice++ ){
		
			int firstBin 	= h2_theta_phi[sector]->GetXaxis()->FindBin(histXmin + slice*sliceWidth);
			int lastBin 	= h2_theta_phi[sector]->GetXaxis()->FindBin(histXmin + (slice+1)*sliceWidth);
			TH1D* proj_y 	= (TH1D*)h2_theta_phi[sector]->
						ProjectionY(Form("h2_theta_phi_Y_%i_%i",sector,slice),firstBin,lastBin);

			if( proj_y->GetEntries() < 10000. ) continue;
			
			// For this, let's fit the histogram to linear and extrapolate to y=0:
			double xmin = proj_y->GetBinCenter(proj_y->FindFirstBinAbove( proj_y->GetMaximum()*0.1 ));
			double xmax = proj_y->GetBinCenter(proj_y->FindFirstBinAbove( proj_y->GetMaximum()*0.9 ));
			TFitResultPtr fit = (TFitResultPtr) proj_y->Fit("pol1","QESR","",xmin,xmax);
			double theta = abs(fit->Parameter(0) / fit->Parameter(1));
			double phi = (histXmin + (slice+0.5)*sliceWidth);	
			
			thetas.push_back(theta);
			phis.push_back(phi);

			// Save the projection
			//proj_y->Write();
		}
		fiducialsPhiSlice[sector] = new TGraph( thetas.size() , &phis[0], &thetas[0] );
		fiducialsPhiSlice[sector]->SetName(Form("fiducialsPhiSlice_%i",sector));
		//fiducialsPhiSlice[sector]->Write();
	}


	// Write the histograms
	outFile->Close();



	return 1;
}

	
