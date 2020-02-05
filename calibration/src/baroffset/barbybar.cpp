#include <cstdlib>
#include <iostream>
// Root includes:
#include "TFile.h"
#include "TStyle.h"
#include "TGraphErrors.h"
// Dependencies includes:
#include "reader.h"
#include "bank.h"
// Calibration includes:
#include "calibclass.h"
#include "helpers.h"

using namespace std;

// Number of bars per sector-layer
int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
void offsetCorr(std::vector<double> *adcs		,
		std::vector<double> *adcsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		double widthCut				,
		TH2F * hist				);
void fitOffset(TH2F * hist , TCanvas * c, int cd, double &par1, double &par2 , double &par1_err, double &par2_err );

// main:
int main(int argc, char** argv) {
	// Check number of arguments provided:
	if( argc!= 3 ){
		cerr << "Wrong number of arguments. Please use instead:\n\t./code [InputFile] [OutputPDFName]\n";
		return -1;
	}
	gStyle->SetTitleSize(0.2,"t");
	gStyle->SetOptStat(0);

	// Read in input file
	TString inputFile = argv[1];
	hipo::reader reader;
	reader.open(inputFile);
        hipo::dictionary  factory;      
        reader.readDictionary(factory); 
	hipo::bank BAND_ADC  (factory.getSchema("BAND::adc"  ));
	hipo::bank BAND_TDC  (factory.getSchema("BAND::tdc"  ));
	hipo::bank RUN_config  (factory.getSchema("RUN::config"  ));
        hipo::event readevent;
	cout << "\n\tLoaded input file\n";

	// Load calibration class and parameters
	calibclass * BAND = new calibclass();
	BAND->LoadTimeWalk();
	BAND->LoadLROffsets();
	BAND->LoadVelocityMap();
	cout << "\n\tLoaded BAND calibrations\n\n";

	// Define histograms to be filled:
	const int nHistos = 600;
	TH2F ** h2_tdc_tmean_gmadc = new TH2F * [nHistos];
	TH2F ** h2_ftdc_tmean_gmadc = new TH2F * [nHistos];
	TH2F ** h2_tdc_tmean_gmamp = new TH2F * [nHistos];
	TH2F ** h2_ftdc_tmean_gmamp = new TH2F * [nHistos];
	for(int i = 0 ; i < nHistos ; i++){
		h2_tdc_tmean_gmadc[i] = new TH2F(Form("h2_tdc_tmean_gmadc_%i",i),"	;sqrt(ADC_{L}*ADC_{R});(t_{TDC,L}+t_{TDC,R})/2. - ref  [ns]	",		500,0,20000,800,-20,20);
		h2_tdc_tmean_gmamp[i] = new TH2F(Form("h2_tdc_tmean_gmamp_%i",i),"	;sqrt(AMP_{L}*AMP_{R});(t_{TDC,L}+t_{TDC,R})/2. - ref  [ns]	",		500,0,5000, 800,-20,20);
                                                                                                              
		h2_ftdc_tmean_gmadc[i] = new TH2F(Form("h2_ftdc_tmean_gmadc_%i",i),"	;sqrt(ADC_{L}*ADC_{R});(t_{FADC,L}+t_{FADC,R})/2. - ref  [ns]	",		500,0,20000,800,-20,20);
		h2_ftdc_tmean_gmamp[i] = new TH2F(Form("h2_ftdc_tmean_gmamp_%i",i),"	;sqrt(AMP_{L}*AMP_{R});(t_{FADC,L}+t_{FADC,R})/2. - ref  [ns]	",		500,0,5000, 800,-20,20);
	}

	// Loop over all events in file
	int event_counter = 0;
	while(reader.next()==true){
		if(event_counter%10000==0) cout << "event: " << event_counter << endl;
		if( event_counter > 150000 ) break;
		event_counter++;

		// Load data structure for this event
                reader.read(readevent); 
                readevent.getStructure(BAND_ADC );   
		readevent.getStructure(BAND_TDC );   
		readevent.getStructure(RUN_config );

		// Get phase correction for TDC-FADC difference
		long timestamp = RUN_config.getLong(4,0);
		double phaseCorr = getTriggerPhase(timestamp);

		// Get number of ADCs and TDCs to do event skimming:
		int nADC = BAND_ADC.getRows();
		int nTDC = BAND_TDC.getRows();
			// Skip events with no entries
		if(nADC==0||nTDC==0) continue;
			// Skip non-laser events
		if(nADC<90||nTDC<90) continue;

		// Create PMTs
		BAND->CreatePMTs( BAND_ADC, BAND_TDC, phaseCorr );
		// Create Bars from PMTs
		BAND->CreateBars();
		// Prepare maps to access information from bars
		BAND->InitGettersBars();
		
		// Access information we need from bars
		std::map<int,double> tdc_tmean = BAND->GetBars_TDC_tMean();
		std::map<int,double> ftdc_tmean = BAND->GetBars_FADC_tMean();
		std::map<int,double> adc = BAND->GetBars_FADC_gmAdc();
		std::map<int,double> amp = BAND->GetBars_FADC_gmAmp();

		// Loop through maps and fill our histograms:
		std::map<int,double>::iterator iter_bar;
		for( iter_bar = tdc_tmean.begin() ; iter_bar != tdc_tmean.end() ; iter_bar++ ){
			int barKey = iter_bar->first;
			int component = (barKey % 10);
			int layer = (barKey - component)/10 % 10;
			int sector = (barKey - component - layer*10)/100 % 10;

			int ref_bar = 200 + layer*10 + 1;
			if( 	tdc_tmean.count(ref_bar) == 0 || ftdc_tmean.count(ref_bar) == 0 ||
				adc.count(ref_bar) == 0 || amp.count(ref_bar) == 0 ) continue;

			h2_tdc_tmean_gmamp[barKey] -> Fill (  amp[barKey] , tdc_tmean[barKey]  - tdc_tmean[ref_bar]  );
			h2_tdc_tmean_gmadc[barKey] -> Fill (  adc[barKey] , tdc_tmean[barKey]  - tdc_tmean[ref_bar]  );
			h2_ftdc_tmean_gmamp[barKey] -> Fill ( amp[barKey] , ftdc_tmean[barKey] - ftdc_tmean[ref_bar] );
			h2_ftdc_tmean_gmadc[barKey] -> Fill ( adc[barKey] , ftdc_tmean[barKey] - ftdc_tmean[ref_bar] );
		}

	} // end loop over events

	// Create canvases to output
	TString outputName = argv[2];
	TCanvas * c0 = new TCanvas("c0","c0",700,900);
	TCanvas * c1 = new TCanvas("c1","c0",700,900);
	c0 -> Modified();
	c0 -> Update(); 
	c1 -> Modified();
	c1 -> Update(); 
	TCanvas *** cSLC_tdc = new TCanvas**[5];
	TCanvas *** cSLC_ftdc = new TCanvas**[5];
	c0 -> Print(Form("tdc_%s(" ,outputName.Data()));
	c1 -> Print(Form("ftdc_%s(",outputName.Data()));
	// 	Loop through all bars and draw
	for(int is = 0 ; is < 5 ; is++){
		cSLC_tdc[is] = new TCanvas*[6];
		cSLC_ftdc[is] = new TCanvas*[6];
		for(int il = 0 ; il < 6 ; il++){
			cSLC_tdc[is][il] = new TCanvas(Form("TDC_S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),700,900);
			cSLC_ftdc[is][il] = new TCanvas(Form("FADC_S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),700,900);
			cSLC_tdc[is][il] -> Divide(3,7);
			cSLC_ftdc[is][il] -> Divide(3,7);

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);
				if( h2_tdc_tmean_gmamp[identifier]->Integral()  ){
					zoomTH2F( h2_tdc_tmean_gmadc[identifier], 10000, 1.5, 1.5 );	drawTH2F( h2_tdc_tmean_gmadc[identifier], cSLC_tdc[is][il]   , 3*cIdx+1 );
					zoomTH2F( h2_tdc_tmean_gmamp[identifier], 2000, 1.5, 1.5 ); 	drawTH2F( h2_tdc_tmean_gmamp[identifier], cSLC_tdc[is][il]   , 3*cIdx+2 );
					double temp;
					fitOffset( h2_tdc_tmean_gmadc[identifier], cSLC_tdc[is][il], 3*cIdx+3, temp, temp, temp, temp );
				}
				if( h2_ftdc_tmean_gmamp[identifier]->Integral()  ){
					zoomTH2F( h2_ftdc_tmean_gmadc[identifier], 10000, 1.5, 1.5 );	drawTH2F( h2_ftdc_tmean_gmadc[identifier], cSLC_ftdc[is][il]   , 3*cIdx+1 );
					zoomTH2F( h2_ftdc_tmean_gmamp[identifier], 2000, 1.5, 1.5 ); 	drawTH2F( h2_ftdc_tmean_gmamp[identifier], cSLC_ftdc[is][il]   , 3*cIdx+2 );
					double temp;
					fitOffset( h2_ftdc_tmean_gmadc[identifier], cSLC_ftdc[is][il], 3*cIdx+3, temp, temp, temp, temp );
				}
			}
			cSLC_tdc[is][il] -> Modified();     
			cSLC_tdc[is][il] -> Update();
			cSLC_tdc[is][il] -> Print(Form("tdc_%s",outputName.Data()));
			cSLC_ftdc[is][il] -> Modified();     
			cSLC_ftdc[is][il] -> Update();
			cSLC_ftdc[is][il] -> Print(Form("ftdc_%s",outputName.Data()));
		}
	}
	c0 -> Print(Form("tdc_%s)",outputName.Data()));
	c1 -> Print(Form("ftdc_%s)",outputName.Data()));

	return 0;
}

void fitOffset(TH2F * hist , TCanvas * c, int cd, double &par1, double &par2 , double &par1_err, double &par2_err ){
	std::vector<double> xs;
	std::vector<double> ys;
	std::vector<double> xErrs;
	std::vector<double> yErrs;
	std::vector<double> res;
	std::vector<double> resErrs;
	offsetCorr( &xs, &xErrs, &ys, &yErrs, &res, &resErrs, 4000, hist );

	int dim = xs.size();
	TGraphErrors *g = new TGraphErrors(dim, &xs[0], &ys[0], &xErrs[0], &res[0]);
	//TF1 * model = new TF1("timeWalk",wlk,0,20000,2);
	//TFitResultPtr ptr = g->Fit(model,"QES");
	//par1 = ptr->Parameter(0);
	//par1_err = ptr->ParError(0);
	//par2 = ptr->Parameter(1);
	//par2_err = ptr->ParError(1);
		
	c->cd(cd);
	g->SetMarkerStyle(20);
	g->Draw("AP");
	
	return;
}

void offsetCorr(std::vector<double> *adcs		,
		std::vector<double> *adcsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		double widthCut				,
		TH2F * hist				){
	int currBin = 0;
	while( currBin < 500 ){
		double xPt, yPt, yEr, ySig, ySigEr;
		int step = doProj( hist , currBin , false , 0, xPt, yPt, yEr, ySig, ySigEr , 500, 500 );
		currBin += step ;
		if( xPt < widthCut) continue;
		adcs		->push_back(xPt);
		adcsErr		->push_back(0);
		times		->push_back(yPt);
		timesErr	->push_back(yEr);
		res		->push_back(ySig);
		resErr		->push_back(ySigEr);
	}
}
