#include "helpers.h"

// Pretty histograming ----------------------------------------------
void zoomTH2F(TH2F * h2, double zoomAround, double zoomLower, double zoomUpper ) {
	h2 -> GetXaxis() -> CenterTitle();
	h2 -> GetYaxis() -> CenterTitle();

	h2 -> SetTitleSize(0.5);

	h2 -> GetXaxis() -> SetLabelSize(0.09);
	h2 -> GetYaxis() -> SetLabelSize(0.09);
	h2 -> GetXaxis() -> SetTitleSize(0.09);
	h2 -> GetYaxis() -> SetTitleSize(0.09);

	h2 -> GetYaxis() -> SetTitleOffset(0.80);

	h2 -> GetYaxis() -> SetNdivisions(109);
	h2 -> GetXaxis() -> SetNdivisions(107);

	int midbin = h2->GetXaxis()->FindBin(zoomAround);
	double val = h2->ProjectionY("trash",midbin,midbin+1)->GetMean();
	double min_val = val - zoomLower;
	double max_val = val + zoomUpper;
	h2->GetYaxis() -> SetRangeUser(min_val,max_val);
}

// Trigger phase function ----------------------------------------------
double getTriggerPhase( long timeStamp ) {
	double tPh = 0.;

	long Period = 4.0;
	long Cycles = 6.0;
	long Phase  = 3.0;

	if( timeStamp != -1 ) 
		tPh = (double)(Period *( (timeStamp + Phase) % Cycles ));

	return tPh;
}

void drawTH2F( TH2F * hist , TCanvas * c , int cd ){
	c -> cd(cd);
	gPad -> SetLeftMargin(0.16);
	gPad -> SetBottomMargin(0.20);
	gPad -> SetTopMargin(0);
	hist->Draw("COLZ");
	return;
}

// Making projection of a TH2F by doing slices in x and fitting y to a gauss ----------------------------------------------
int doProj( TH2F * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE, int thres, int lastBin ){
	int cnt = 0;
	int step = 0;
	TCanvas * trash = new TCanvas("trash");
	while ( cnt < thres ){
		char temp[100];
		sprintf(temp,"slice_%d_%d",flag,bin);
		TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
		cnt = (int) pj->Integral();
		delete pj;
		if( (bin+step) >= lastBin) break;
		step+=1;
	}
	
	if( cnt >= thres ){
		char temp[100];
		sprintf(temp,"slice_%d_%d",flag,bin);
		TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
		
		// Getting the mean x value in this range:
		hist->GetXaxis()->SetRange(bin,bin+step);
		x = hist->GetMean(1);
		hist->GetXaxis()->SetRange();
		TFitResultPtr f = pj->Fit("gaus","QESR","",-300,300);
		y = f->Parameter(1);
		yE = f->ParError(1);
		sig = f->Parameter(2);
		sigE = f->ParError(2);
		
		if( write ){
			pj->Write();
		}
		delete pj;
	}

	delete trash;

	return step;
}
