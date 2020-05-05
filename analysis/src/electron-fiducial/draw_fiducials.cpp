#include <iostream>
#include <cmath>
#include <vector>

#include "TFile.h"
#include "TCut.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TCutG.h"
#include "TChain.h"
#include "TEventList.h"
#include "TCanvas.h"

using namespace std;

int main( int argc, char** argv){
	if( argc != 2){
		cerr << "Incorrect number of arguments. Please instead use:\n"
			<< "\t./code [InclusiveDataDirectory] \n";
		return -1;
	}

	/////////////////////////////////////////////////////////////////////
	// Setup the input root file
//	TFile * inFile = new TFile(argv[1]);
//	TTree * inTree = (TTree*)inFile->Get("skim");
	
	// Use TChain instead
	TChain* inTree = new TChain("skim");
	inTree->Add(Form("%s/*.root",argv[1]));


	/////////////////////////////////////////////////////////////////////
	// Setup the output root file and histograms
//	TH1F* pHist = new TH1F("pHist","",100,pMin-0.5,pMax+0.5);
	
	TCut electron_cut = "E_tot/p_e > 0.17 && E_tot/p_e < 0.27 && E_tot > 0.25 && e_vtz > -6.0 && e_vtz < 0 && p_e < Ebeam";
	TCut sector_cut;
	TCut momentum_cut;
	
	TEventList* goodElectrons = new TEventList("goodElectrons");
	inTree->Draw(">>goodElectrons",electron_cut);
	inTree->SetEventList(goodElectrons);

	const int nSector = 6;
	const int nMom = 12;
	int pStart = 1250;	// MeV
	double pStep0 = 250;	// MeV
	double pScale = 1.1;	// rate of bin increase with momentum

	TH2F* acc_hist[nSector][nMom];
	TH1F* p_hist[nSector][nMom];

	TCanvas* c1 = new TCanvas("c1","",2400,1200);
	c1->Divide(3,2);

	int pMinNext = pStart;

	for(int p = 0; p < nMom; p++) {

		int pMin = pMinNext;
		int pMax = pMin + pStep0*pow(pScale, p);	

		cout << "Looking at electrons between " << pMin << " and " << pMax << " MeV" << endl;

		momentum_cut = Form("p_e > %f && p_e < %f", pMin/1000., pMax/1000.);

		for(int sector = 0; sector < nSector; sector++) {

			sector_cut = Form("eSector == %i", sector+1);

			TString acc_draw = Form("theta_e*TMath::RadToDeg():phi_e*TMath::RadToDeg()>>acc_sect_%i_p%i", sector,p);

			if(sector == 3) {
				acc_draw = Form("theta_e*TMath::RadToDeg():(phi_e*TMath::RadToDeg()-(phi_e > 0)*360.)>>acc_sect_%i_p%i", sector, p);
			}
			TString p_draw = Form("p_e>>p_sect_%i_p%i",sector,p);

			acc_hist[sector][p] = new TH2F(Form("acc_sect_%i_p%i",sector,p), "", 420, -220, 200, 160, 0, 40);
			p_hist[sector][p] = new TH1F(Form("p_sect_%i_p%i",sector,p) , "", 100, pMin/1000., pMax/1000.);

			inTree->Draw(p_draw, momentum_cut && sector_cut, "GOFF");
			p_hist[sector][p]->GetXaxis()->SetTitle("p_{e} (GeV)");
			p_hist[sector][p]->SetTitle(Form("Sector %i", sector+1));
			double pAvg = p_hist[sector][p]->GetMean();
			c1->cd(sector+1);
			inTree->Draw(acc_draw, momentum_cut && sector_cut, "COLZ");
			acc_hist[sector][p]->GetXaxis()->SetTitle("#phi (deg)");
			acc_hist[sector][p]->GetYaxis()->SetTitle("#theta (deg)");
			acc_hist[sector][p]->SetTitle(Form("Sector %i, p = %.2f GeV", sector+1, pAvg));
		}

		c1->SaveAs(Form("electron_acc_p%i.pdf", p));
		
		pMinNext = pMax;

	} // momentum for loop


	TString outname = "/Users/tylerkutz/research/BAND/fiducial/hist_out.root";
	TFile * outFile = new TFile(outname,"RECREATE");

	for(int sector = 0; sector < nSector; sector++) {
		for(int p = 0; p < nMom; p++) {
			acc_hist[sector][p]->Write();
			p_hist[sector][p]->Write();
		}
	}


	outFile->Close();

	return 1;
}

	
