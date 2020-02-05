#include <cstdlib>
#include <iostream>

#include "TRint.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "reader.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

using namespace std;

double startADC = 2000;
double endADC = 40000;
double decayRate = 10000;
int cnt = 0;

// Forward-declaring functions
void PrettyTH1F(TH1F * h1);
TFitResultPtr changeFit(TH1F * h1, double land_amp, double land_mean, double land_sigma, double exp_amp, double exp_decay);
TFitResultPtr doFit(TH1F * h1);
int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
// ========================================================================================================================================
int main(int argc, char** argv) {
	// Declaring histograms
	const int nArray = 600;

	/*
#ifdef WITHRINT
TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
TApplication *myapp = new TApplication("myapp",0,0);
#endif
*/
	gStyle->SetTitleSize(0.3,"t");
	gStyle->SetOptStat(0);

	TString inputFile;
	TString outPDFFile;
	TString outParamFile;

	if(argc==4) {
		inputFile = argv[1];
		outPDFFile = argv[2];
		outParamFile = argv[3];
	}
	else {
		cout << "=========================\nRun this code as:\n./specFit path/to/input/Hippo /path/to/output/pdf /path/to/output/Param/text\n=========================" << endl;
		exit(0);
	}

	// ----------------------------------------------------------------------------------
	double parL[nArray][5] = {{0}};
	double parR[nArray][5] = {{0}};

	TH1F ** h1_adc_spec_L = new TH1F * [nArray];
	TH1F ** h1_adc_spec_R = new TH1F * [nArray];
	for(int i = 0 ; i < nArray ; i++){
		h1_adc_spec_L[i] = new TH1F(Form("h1_adc_spec_L_%i",i),";ADC",800,startADC,endADC);
		h1_adc_spec_R[i] = new TH1F(Form("h1_adc_spec_R_%i",i),";ADC",800,startADC,endADC);
		PrettyTH1F(h1_adc_spec_L[i]);
		PrettyTH1F(h1_adc_spec_R[i]);
	}

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	//Read Dictionary of Hipo File  // new hipo4
	hipo::dictionary  factory;      // new hipo4
	reader.readDictionary(factory); // new hipo4
	//factory.show();               // new hipo4

	hipo::bank BAND_ADC  (factory.getSchema("BAND::adc"  ));	// new hipo4
	//hipo::bank BAND_TDC  (factory.getSchema("BAND::tdc"  ));	// new hipo4
	//hipo::bank RUN_config(factory.getSchema("RUN::config"));	// new hipo4

	//One also needs a hipo::event object which is called from the reader for each event to get
	//the information for each bank
	hipo::event readevent;  // new hipo4

	int event_counter = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while((reader.next()==true)){

		//Reader has to load information about event in hipo::event class
		reader.read(readevent); // new hipo4

		//Load explicitly all information for each bank for the event
		readevent.getStructure(BAND_ADC   );   // new hipo4
		//readevent.getStructure(BAND_TDC   );   // new hipo4

		//Now everything is loaded and can be used as before with HIPO3 files. There is only one difference:
		//The number of hits in each bank is determined by the function "getRows()" and not by "getSize" as before.

		if(event_counter%1000000==0) cout << "event: " << event_counter << endl;

		event_counter++;	
		int nADC = BAND_ADC.getRows();
		//int nTDC = BAND_TDC.getRows();

		// Skip events with no entries
		if(nADC==0) continue;

		for(int aIdx = 0 ; aIdx < nADC ; aIdx++){

			int   ADC_sector    = BAND_ADC.getInt  (0,aIdx);
			int   ADC_layer     = BAND_ADC.getInt  (1,aIdx);
			int   ADC_component = BAND_ADC.getInt  (2,aIdx);
			int   ADC_order     = BAND_ADC.getInt  (3,aIdx);
			float ADC_adc       = (float)(BAND_ADC.getInt(4,aIdx));
			float ADC_amp       = (float)(BAND_ADC.getInt(5,aIdx));
			float ADC_time      = BAND_ADC.getFloat(6,aIdx);
			if( ADC_amp >= 4095) continue; //overflow event

			int barKey = 100*ADC_sector + 10*ADC_layer + ADC_component;

			// Left PMTs
			if(ADC_order==0){
				h1_adc_spec_L[barKey] -> SetTitle(Form("Left PMTs, Sector:%i, Layer:%i, Component:%i",ADC_sector,ADC_layer,ADC_component));
				h1_adc_spec_L[barKey] -> Fill(ADC_adc);
			}

			// Right PMTs
			if(ADC_order==1){
				h1_adc_spec_R[barKey] -> SetTitle(Form("Right PMTs, Sector:%i, Layer:%i, Component:%i",ADC_sector,ADC_layer,ADC_component));
				h1_adc_spec_R[barKey] -> Fill(ADC_adc);
			}			
		}
	}// end file
	// -------------------------------------------------------------------------------------------------
	// Printing results onto canvases
	TCanvas * c0 = new TCanvas("c0");
	c0 -> Divide(2,1);
	c0 -> Modified();
	c0 -> Update();

	TCanvas *** cSLC = new TCanvas**[5];

	for(int is = 0 ; is < 5 ; is++){
		cSLC[is] = new TCanvas*[6];
		for(int il = 0 ; il < 6 ; il++){
			cSLC[is][il] = new TCanvas(Form("S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);
			cSLC[is][il] -> Divide(2,7);

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);
				int notEmpty = (h1_adc_spec_L[identifier]->Integral()+h1_adc_spec_R[identifier]->Integral());
				if(notEmpty){

					int min, max;

					cSLC[is][il] -> cd(2*cIdx+1);
					gPad -> SetBottomMargin(0.26);
					h1_adc_spec_L[identifier] -> Draw("COLZ");

					cout << "Fitting slco: " << (is+1) << " " << (il+1) << " " << (cIdx+1) << "\n";
					if( h1_adc_spec_L[identifier]->Integral() ){
						TFitResultPtr Fit_L = doFit(h1_adc_spec_L[identifier]);
						for(int j = 0; j < 5; j++){
							if(Fit_L == 0){  parL[identifier][j] = Fit_L->Parameter(j); }
						}
					}

					cSLC[is][il] -> cd(2*cIdx+2);
					gPad -> SetBottomMargin(0.26);
					h1_adc_spec_R[identifier] -> Draw("COLZ");
					if( h1_adc_spec_R[identifier]->Integral() ){
						TFitResultPtr Fit_R = doFit(h1_adc_spec_R[identifier]);
						for(int j = 0; j < 5; j++){
							if(Fit_R == 0){  parR[identifier][j] = Fit_R->Parameter(j); }
						}
					}

				}
				cSLC[is][il] -> Modified();	cSLC[is][il] -> Update();
			}
		}
	}

	// -------------------------------------------------------------------------------------------------
	// Deleting empty histograms
	for(int i = 0 ; i < nArray ; i++){
		int notEmpty = (h1_adc_spec_L[i]->Integral()+h1_adc_spec_R[i]->Integral());
		if(!notEmpty){
			delete h1_adc_spec_L[i];
			delete h1_adc_spec_R[i];
		}
	}

	// -------------------------------------------------------------------------------------------------
	// Saving plots to a pdf file
	c0 -> Print(outPDFFile + "(");
	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 6 ; il++){
			cSLC[is][il] -> Print(outPDFFile);
		}
	}
	c0 -> Print(outPDFFile + ")");


	// -------------------------------------------------------------------------------------------------
	// Saving fit values to ccdb tables
	ofstream outText;
	outText.open(outParamFile);

	outText << "# Sector\t Layer\t Component\t Order\t Land_Amp\t Land_Mean\t Land_Sigma\t Exp_Amp\t Exp_Decay\t \n";

	for(int is = 1 ; is <= 5 ; is++){	
		for(int il = 1 ; il <= 6 ; il++){
			for(int ic = 1 ; ic <= slc[il-1][is-1] ; ic++){
				int idx = 100*is + 10*il + ic;

				outText << is << " " << il << " " << ic << " " << 0 << " ";
				outText << parL[idx][0] << " " << parL[idx][1] << " " << parL[idx][2] << " " << parL[idx][3] << " " << parL[idx][4] << " "  << endl;
				//outText << parL[idx][1] << " " << parL[idx][2] << endl;
				// ---
				outText << is << " " << il << " " << ic << " " << 1 << " ";
				outText << parR[idx][0] << " " << parR[idx][1] << " " << parR[idx][2] << " " << parR[idx][3] << " " << parR[idx][4] << " "  << endl;				
				//outText << parR[idx][1] << " " << parR[idx][2] << endl;
			}
		}
	}
	outText.close();

	return 0;
}
// ========================================================================================================================================
TFitResultPtr doFit(TH1F * h1) {

	double land_amp = (h1->GetMaximum());
	double land_mean = 10000;//h1->GetBinCenter(h1->GetMaximumBin());
	//You can use the mean or the mode for the landau guess
	//double land_mean = h1_adc_spec_L[i]->GetMean();
	double land_sigma = h1->GetStdDev();
	double exp_amp = (h1->GetBinContent(1)) * TMath::Exp(startADC / decayRate);
	double exp_decay = decayRate;

	return changeFit(h1, land_amp, land_mean, land_sigma, exp_amp, exp_decay);

}
// ========================================================================================================================================
TFitResultPtr changeFit(TH1F * h1, double land_amp, double land_mean, double land_sigma, double exp_amp, double exp_decay) {
	cnt++;
	TF1 * f_Spec = new TF1(Form("f_cosmic_spec_%i",cnt),"[&](double *x, double *p) { return p[0]*TMath::Landau(x[0],p[1],p[2]) + p[3]*TMath::Exp(-x[0]/p[4]) ; }",startADC, endADC,5);
	f_Spec->SetParameter(0,land_amp);
	f_Spec->SetParLimits(0,0,10*(h1->GetMaximum()));

	f_Spec->SetParameter(1,land_mean);
	f_Spec->SetParLimits(1,startADC,endADC);

	f_Spec->SetParameter(2,land_sigma);			
	f_Spec->SetParLimits(2,0,5000);

	f_Spec->SetParameter(3,exp_amp);			
	f_Spec->SetParameter(4,exp_decay);			

	TFitResultPtr fitStatus = h1 -> Fit(Form("f_cosmic_spec_%i",cnt),"Qesr");

	if(fitStatus == 0){
		f_Spec->Draw("same");
	}

	return fitStatus;

}
// ========================================================================================================================================
void PrettyTH1F(TH1F * h1) {
	h1 -> GetXaxis() -> CenterTitle();

	h1 -> SetTitleSize(0.5);

	h1 -> GetXaxis() -> SetLabelSize(0.13);
	h1 -> GetXaxis() -> SetTitleSize(0.13);

	h1 -> GetXaxis() -> SetNdivisions(509);
}

/*
   TF1 * f_timewalk_L = new TF1("f_cosmic_spec_L","[land_amp]*TMath::Landau(x,[land_mean],[land_sigma]) +[exp_amp]*exp(-x/[exp_decay])");
   f_timewalk_L->SetParameter(0,h1_adc_spec_L[i]->GetMaximum());
//You can use the mean or the mode for the landau guess
//f_timewalk_L->SetParameter(1,h1_adc_spec_L[i]->GetMean());
f_timewalk_L->SetParameter(1,h1_adc_spec_L[i]->GetBinCenter(h1_adc_spec_L[i]->GetMaximumBin()));
f_timewalk_L->SetParameter(2,h1_adc_spec_L[i]->GetStdDev());			
f_timewalk_L->SetParameter(3,0.5 * (h1_adc_spec_L[i]->GetMaximum()));			
f_timewalk_L->SetParameter(4,500);			
TFitResultPtr fitStatus_L = h1_adc_spec_L[i] -> Fit("f_cosmic_spec_L","Qes");

// ========================================================================================================================================
TFitResultPtr changeFitTest(TH1F * h1, double land_amp, double land_mean, double land_sigma, double exp_amp, double exp_decay) {

//TF1 * f_Spec = new TF1("f_cosmic_spec","(700*TMath::Landau(x,15000,1000)) + ([0]*exp(-x/[1]))",startADC,30000);
TF1 * f_Spec = new TF1("f_cosmic_spec_test",[&](double *x, double *par) { return 700*TMath::Landau(x[0],15000,5000) + par[0]*TMath::Exp(-x[0]/par[1]) ; },startADC, 30000,2);
//  TF1 * f_Spec = new TF1("f_cosmic_spec","[land_amp]*TMath::Landau(x,[land_mean],[land_sigma])",startADC,30000);
cout <<"Starting params:\n";
cout << land_amp << " " << land_mean << " " << land_sigma << " " << exp_amp << " " << exp_decay << " \n";
f_Spec->SetParameter(0,exp_amp);			
f_Spec->SetParameter(1,exp_decay);			

f_Spec->Draw("same");
//TFitResultPtr fitStatus = h1 -> Fit("f_cosmic_spec","Qesr");


for(int j = 0; j < 5; j++){
//if(fitStatus == 0){  cout<<fitStatus->Parameter(j)<<" "; }                                                                                                                  
} 
cout<<"\n";
return 1; //fitStatus;

}


//TF1 * f_Spec = new TF1("f_cosmic_spec","([0]*TMath::Landau(x,[1],[2])) + ([3]*exp(-x/[4]))",startADC,30000);
//TF1 * f_Spec = new TF1(Form("f_cosmic_spec_%i",cnt),"([0]*TMath::Landau(x,[1],[2])) + ([3]*exp(-x/[4]))",startADC,30000,5);
//TF1 * f_Spec = new TF1(Form("f_cosmic_spec_%i",cnt),"([0]*TMath::Landau(x,[1],[2])) + ([3]*exp(-x/[4]))",startADC,30000);

//  TF1 * f_Spec = new TF1("f_cosmic_spec","[land_amp]*TMath::Landau(x,[land_mean],[land_sigma])",startADC,30000);


*/
