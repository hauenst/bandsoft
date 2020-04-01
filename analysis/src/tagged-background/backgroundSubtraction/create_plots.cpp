#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TCut.h"


using namespace std;
void PrettyHist( TH1D * hist , bool sig , double weight , bool xp, int thisSlice );

int main( int argc, char** argv){
	if( argc != 4 ){
		cerr << "Wrong number of arguments. Instead use:\n\t./create_plots [outputRootFile] [SpB RootFile] [B RootFile]\n";
		return -1;
	}

	// Grab input and output files
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TFile * inFile_spb	= new TFile(argv[2]);
	TFile * inFile_b	= new TFile(argv[3]);
	
	// Load trees
	TTree * inTree_spb	= (TTree*)inFile_spb->Get("calib");
	TTree * inTree_b	= (TTree*)inFile_b->Get("mixed");
	TH1D  * hWp_offtime	= (TH1D*)inFile_spb->Get("hWp_b");
	TH1D  * hWp_spb		= (TH1D*)inFile_spb->Get("hWp_spb");
	TH1D  * hWp_mixed	= new TH1D("hWp_mixed","hWp_mixed",500,0,5);
	
	// Define background weights:
	double const background_weight = (16645.940) / inTree_b->GetEntries();
	
	PrettyHist( hWp_offtime, false, (16645.940)/hWp_offtime->Integral()  ,false, 0);
	hWp_offtime->SetTitle("");
	inTree_b->Draw("Wp >> hWp_mixed");
	PrettyHist( hWp_mixed, false, background_weight  ,false, 0);
	hWp_mixed->SetTitle("");
	outFile->cd();
	hWp_mixed->Write();
	hWp_offtime->Write();
	
	// Define histograms to use:
	int x_slices = 10;
	TH1D ** ToF_spb_xslices 	= new TH1D*[x_slices];
	TH1D ** ToF_spb_xpslices	= new TH1D*[x_slices];
	TH1D ** ToF_b_xslices 		= new TH1D*[x_slices];
	TH1D ** ToF_b_xpslices		= new TH1D*[x_slices];
	TCut * cuts_xslices	= new TCut[x_slices];
	TCut * cuts_xpslices	= new TCut[x_slices];
	for( int thisSlice = 0 ; thisSlice < x_slices ; thisSlice++ ){
		ToF_spb_xslices[thisSlice]	= new TH1D(Form("ToF_spb_xslices_%i",thisSlice)	,Form("ToF_spb_xslices_%i",thisSlice),		60,13,63);
		ToF_b_xslices[thisSlice]	= new TH1D(Form("ToF_b_xslices_%i",thisSlice)	,Form("ToF_b_xslices_%i",thisSlice),		60,13,63);
		ToF_spb_xpslices[thisSlice]	= new TH1D(Form("ToF_spb_xpslices_%i",thisSlice)	,Form("ToF_spb_xpslices_%i",thisSlice),	60,13,63);
		ToF_b_xpslices[thisSlice]	= new TH1D(Form("ToF_b_xpslices_%i",thisSlice)		,Form("ToF_b_xpslices_%i",thisSlice),	60,13,63);

		cuts_xslices[thisSlice] = Form("sqrt(W2)>2 && fabs(xB - %f)<0.05",0.05+thisSlice*0.1);
		cuts_xpslices[thisSlice] = Form("Wp>2 && fabs(Xp - %f)<0.05",0.05+thisSlice*0.1);

		// Draw signal + background x slice
		inTree_spb	->Draw(Form("ToF >> ToF_spb_xslices_%i",thisSlice) 	, cuts_xslices[thisSlice] );
		inTree_b	->Draw(Form("ToF >> ToF_b_xslices_%i",thisSlice) 	, cuts_xslices[thisSlice] );
		// Draw signal + background xp slice
		inTree_spb	->Draw(Form("ToF >> ToF_spb_xpslices_%i",thisSlice) 	, cuts_xpslices[thisSlice] );
		inTree_b	->Draw(Form("ToF >> ToF_b_xpslices_%i",thisSlice) 	, cuts_xpslices[thisSlice] );

		// Do weighting and coloring
		PrettyHist( ToF_spb_xslices[thisSlice], true, background_weight  ,false,thisSlice);
		PrettyHist( ToF_spb_xpslices[thisSlice], true, background_weight ,true,thisSlice);
		PrettyHist( ToF_b_xslices[thisSlice], false, background_weight   ,false, thisSlice);
		PrettyHist( ToF_b_xpslices[thisSlice], false, background_weight  ,true, thisSlice);

		

		outFile->cd();
		ToF_spb_xslices[thisSlice]->Write();
		ToF_b_xslices[thisSlice]->Write();
		ToF_spb_xpslices[thisSlice]->Write();
		ToF_b_xpslices[thisSlice]->Write();
	}




	return 0;
}

void PrettyHist( TH1D * hist , bool sig , double weight , bool xp, int thisSlice ){
	hist->SetLineWidth(3);
	hist->SetMinimum(0);
	if( xp ) hist->SetTitle(Form("ToF for |x'-%.2f|<0.05, W' > 2",0.05+thisSlice*0.1 ));
	else{ hist->SetTitle(Form("ToF for |x-%.2f|<0.05, W > 2",0.05+thisSlice*0.1 )); }
	if( !sig){ hist->SetLineColor(2); hist->Scale( weight ); }
	return;
}
