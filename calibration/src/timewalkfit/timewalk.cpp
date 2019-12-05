#include <cstdlib>
#include <iostream>
// Root includes:
#include "TFile.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TF1.h"
// Dependencies includes:
#include "reader.h"
#include "bank.h"
// Calibration includes:
#include "calibclass.h"
#include "helpers.h"

using namespace std;

// Number of bars per sector-layer
int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
// Bar lengths in each sector
double BARLENGTHS[]  = {163.7,201.9,51.2,51.2,201.9};

void walkCorr(	std::vector<double> *adcs		,
		std::vector<double> *adcsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		double widthCut				,
		TH2F * hist				);
void fitTW(TH2F * hist , double &par1, double &par2 , double &par1_err, double &par2_err );
double wlk( double *x , double *p);

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
	//BAND->LoadTimeWalk();
	//BAND->LoadLROffsets();
	//BAND->LoadVelocityMap();
	//cout << "\n\tLoaded BAND calibrations\n\n";

	// Define histograms to be filled:
	const int nHistos = 600;
	double parA_L[nHistos][2] = {{0}};
	double parB_L[nHistos][2] = {{0}};
	double parA_R[nHistos][2] = {{0}};
	double parB_R[nHistos][2] = {{0}};
	TH2F ** h2_tdc_adc_L = new TH2F * [nHistos];
	TH2F ** h2_tdc_adc_R = new TH2F * [nHistos];
	TH2F ** h2_ftdc_adc_L = new TH2F * [nHistos];
	TH2F ** h2_ftdc_adc_R = new TH2F * [nHistos];
	for(int i = 0 ; i < nHistos ; i++){
		h2_tdc_adc_L[i] = new TH2F(Form("h2_tdc_adc_L_%i",i),";ADC_{L};t_{TDC,L} - ref [ns]",800,0,20000,3000,0,300);
		h2_tdc_adc_R[i] = new TH2F(Form("h2_tdc_adc_R_%i",i),";ADC_{R};t_{TDC,R} - ref [ns]",800,0,20000,3000,0,300);
		h2_ftdc_adc_L[i] = new TH2F(Form("h2_ftdc_adc_L_%i",i),";ADC_{L};t_{FADC,L} - ref [ns]",800,0,20000,3000,0,300);
		h2_ftdc_adc_R[i] = new TH2F(Form("h2_ftdc_adc_R_%i",i),";ADC_{R};t_{FADC,R} - ref [ns]",800,0,20000,3000,0,300);
	}

	// Loop over all events in file
	int event_counter = 0;
	while(reader.next()==true){
		if(event_counter%10000==0) cout << "event: " << event_counter << endl;
		if( event_counter > 100000 ) break;
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

		// Create PMTs
		BAND->CreatePMTs( BAND_ADC, BAND_TDC, phaseCorr );
		// Prepare maps to access information from pmts
		BAND->InitGettersPMTs();
		
		// Access information we need from pmts
		std::map<int,double> adc = BAND->GetPMTs_FADC_adc();
		std::map<int,double> amp = BAND->GetPMTs_FADC_amp();
		std::map<int,double> ftdc = BAND->GetPMTs_FADC_time();
		std::map<int,double> tdc = BAND->GetPMTs_TDC_time();

		// Loop through maps and fill our histograms:
		std::map<int,double>::iterator iter_pmt;
		for( iter_pmt = adc.begin() ; iter_pmt != adc.end() ; iter_pmt++ ){
			int pmtKey = iter_pmt->first;
			int order = (pmtKey % 10);
			int component = (pmtKey - order)/10 % 10;
			int layer = (pmtKey - order - component*10)/100 % 10;
			int sector  = (pmtKey - order - component*10 - layer*100)/1000 % 10;
			int barKey = sector * 100 + layer * 10 + component;
			// 	Reference is PMT in bar 116B_R (layer = 1, sector = 4, comp = 6, order = 3)
			//	( pmtKey == 4161 )
			int ref_photodiode = 4161;

			if( 	adc[pmtKey] == 0 		|| adc[ref_photodiode] == 0 	||
				amp[pmtKey] == 0 		|| amp[ref_photodiode] == 0 	||
				amp[pmtKey] >= 4095 		|| amp[ref_photodiode] >= 4095 	|| 
				tdc[pmtKey] == 0 		|| tdc[ref_photodiode] == 0	||
				ftdc[pmtKey] == 0		|| ftdc[ref_photodiode] == 0 	||
				tdc[pmtKey] == 1e10		|| tdc[ref_photodiode] == 1e10	||
				ftdc[pmtKey] == 1e10		|| ftdc[ref_photodiode] == 1e10  ) continue;
			
			if( order == 0 ){
				h2_tdc_adc_L[barKey] -> Fill(	adc[pmtKey]	,	tdc[pmtKey] - tdc[ref_photodiode] 	);
				h2_ftdc_adc_L[barKey] -> Fill(	adc[pmtKey]	,	ftdc[pmtKey] - ftdc[ref_photodiode] 	);
			}
			else if( order == 1){
				h2_tdc_adc_R[barKey] -> Fill(	adc[pmtKey]	,	tdc[pmtKey] - tdc[ref_photodiode] 	);
				h2_ftdc_adc_R[barKey] -> Fill(	adc[pmtKey]	,	ftdc[pmtKey] - ftdc[ref_photodiode] 	);
			}
		}

	} // end loop over events

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
	double par1_L[600] = {0};
	double par2_L[600] = {0};
	double par1_R[600] = {0};
	double par2_R[600] = {0};
	double par1e_L[600] = {0};
	double par2e_L[600] = {0};
	double par1e_R[600] = {0};
	double par2e_R[600] = {0};
	
	for(int is = 0 ; is < 5 ; is++){
		cSLC_tdc[is] = new TCanvas*[6];
		cSLC_ftdc[is] = new TCanvas*[6];
		for(int il = 0 ; il < 6 ; il++){
			cSLC_tdc[is][il] = new TCanvas(Form("TDC_S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),700,900);
			cSLC_ftdc[is][il] = new TCanvas(Form("FADC_S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),700,900);
			cSLC_tdc[is][il] -> Divide(2,7);
			cSLC_ftdc[is][il] -> Divide(2,7);

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);

				if( h2_tdc_adc_L[identifier]->Integral()  ){
					double par1, par2, par1e, par2e = 0;
					zoomTH2F( h2_tdc_adc_L[identifier] , 10000 , 2, 8 );
					drawTH2F( h2_tdc_adc_L[identifier] , cSLC_tdc[is][il] ,  2*cIdx+1 );
					fitTW( 	  h2_tdc_adc_L[identifier] , par1 , par2 , par1e, par2e );
					// Fit results:
					par1_L[identifier] = par1;
					par2_L[identifier] = par2;
					par1e_L[identifier] = par1e;
					par2e_L[identifier] = par2e;
				}
				if( h2_tdc_adc_R[identifier]->Integral()  ){
					double par1, par2, par1e, par2e = 0;
					zoomTH2F( h2_tdc_adc_R[identifier] , 10000 , 2, 8 );
					drawTH2F( h2_tdc_adc_R[identifier] , cSLC_tdc[is][il] ,  2*cIdx+2 );
					if( identifier!= 416 ) fitTW( 	  h2_tdc_adc_R[identifier] , par1 , par2 , par1e, par2e );
					// Fit results:
					par1_R[identifier] = par1;
					par2_R[identifier] = par2;
					par1e_R[identifier] = par1e;
					par2e_R[identifier] = par2e;
				}
				if( h2_ftdc_adc_L[identifier]->Integral()  ){
					zoomTH2F( h2_ftdc_adc_L[identifier] , 10000 , 4, 4 );
					drawTH2F( h2_ftdc_adc_L[identifier] , cSLC_ftdc[is][il] ,  2*cIdx+1 );
				}
				if( h2_ftdc_adc_R[identifier]->Integral()  ){
					zoomTH2F( h2_ftdc_adc_R[identifier] , 10000 , 4, 4 );
					drawTH2F( h2_ftdc_adc_R[identifier] , cSLC_ftdc[is][il] ,  2*cIdx+2 );
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
	// Print out the results
	ofstream tabL, tabR;
	tabL.open("test_timeWalkPar_L.txt");
	tabR.open("test_timeWalkPar_R.txt");
	for(int is = 1 ; is <= 5 ; is++){
		for(int il = 1 ; il <= 6 ; il++){
			for(int ic = 1 ; ic <= slc[il-1][is-1] ; ic++){
				int idx = 100*is + 10*il + ic;
				tabL << is << "\t" << il << "\t" << ic << "\t" 	<< par2_L[idx] << "\t" << par1_L[idx] << "\t" 
										<< par2e_L[idx]  << "\t" << par1e_L[idx] << "\n";
				tabR << is << "\t" << il << "\t" << ic << "\t" 	<< par2_R[idx] << "\t" << par1_R[idx] << "\t" 
										<< par2e_R[idx]  << "\t" << par1e_R[idx] << "\n";
			}
		}
	}
	tabL.close();
	tabR.close();

	return 0;
}


void fitTW(TH2F * hist , double &par1, double &par2 , double &par1_err, double &par2_err ){
	std::vector<double> xs;
	std::vector<double> ys;
	std::vector<double> xErrs;
	std::vector<double> yErrs;
	std::vector<double> res;
	std::vector<double> resErrs;
	walkCorr( &xs, &xErrs, &ys, &yErrs, &res, &resErrs, 4000, hist );

	int dim = xs.size();
	TGraphErrors *g = new TGraphErrors(dim, &xs[0], &ys[0], &xErrs[0], &yErrs[0]);
	TF1 * model = new TF1("timeWalk",wlk,0,20000,2);
	TFitResultPtr ptr = g->Fit(model,"QES");
	par1 = ptr->Parameter(0);
	par1_err = ptr->ParError(0);
	par2 = ptr->Parameter(1);
	par2_err = ptr->ParError(1);
	return;
}

void walkCorr(	std::vector<double> *adcs		,
		std::vector<double> *adcsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		double widthCut				,
		TH2F * hist				){
	int currBin = 0;
	while( currBin < 800 ){
		double xPt, yPt, yEr, ySig, ySigEr;
		int step = doProj( hist , currBin , false , 0, xPt, yPt, yEr, ySig, ySigEr , 800, 500 );
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

double wlk( double *x , double *p){
	double var = *x;
	return p[0] + p[1] / sqrt(var);
}




/*
#include <cstdlib>
#include <iostream>

#include "TLine.h"
#include "TProfile.h"
#include "TRint.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "reader.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

using namespace std;

// Forward-declaring functions
void PrettyTH1F(TH1F * h1, int color);
void PrettyTH2F(TH2F * h2);
void PrettyTGraphErrors(TGraphErrors * gP, int color);
double getTriggerPhase( long timeStamp );
void walkCorr(	std::vector<double> *adcs		,
		std::vector<double> *adcsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		double widthCut				,
		TH2F * hist				);
int doProj( TH2F * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE);
void fitTW(TH2F * hist , double &par1, double &par2 , double &par1_err, double &par2_err );
double wlk( double *x , double *p);

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

int main(int argc, char** argv) {
	if( argc!= 2 ){
		cerr << "Wrong number of arguments. Please use instead:\n\t./code [InputFile]\n";
		return -1;
	}
	gStyle->SetTitleSize(0.2,"t");
	gStyle->SetOptStat(0);

	TString inputFile = argv[1];


	const int nHistos = 600;
	double parA_L[nHistos][2] = {{0}};
	double parB_L[nHistos][2] = {{0}};
	double parA_R[nHistos][2] = {{0}};
	double parB_R[nHistos][2] = {{0}};
	TH2F ** h2_tdc_adc_L = new TH2F * [nHistos];
	TH2F ** h2_tdc_adc_R = new TH2F * [nHistos];
	TH2F ** h2_ftdc_adc_L = new TH2F * [nHistos];
	TH2F ** h2_ftdc_adc_R = new TH2F * [nHistos];
	for(int i = 0 ; i < nHistos ; i++){
		h2_tdc_adc_L[i] = new TH2F(Form("h2_tdc_adc_L_%i",i),";ADC_{L};t_{TDC,L} - ref [ns]",800,0,20000,3000,0,300);
		h2_tdc_adc_R[i] = new TH2F(Form("h2_tdc_adc_R_%i",i),";ADC_{R};t_{TDC,R} - ref [ns]",800,0,20000,3000,0,300);
		h2_ftdc_adc_L[i] = new TH2F(Form("h2_ftdc_adc_L_%i",i),";ADC_{L};t_{FADC,L} - ref [ns]",800,0,20000,3000,-1200,-900);
		h2_ftdc_adc_R[i] = new TH2F(Form("h2_ftdc_adc_R_%i",i),";ADC_{R};t_{FADC,R} - ref [ns]",800,0,20000,3000,-1200,-900);
		PrettyTH2F(h2_tdc_adc_L[i]);
		PrettyTH2F(h2_tdc_adc_R[i]);
		PrettyTH2F(h2_ftdc_adc_L[i]);
		PrettyTH2F(h2_ftdc_adc_R[i]);
	}

	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);
        hipo::dictionary  factory;      // new hipo4
        reader.readDictionary(factory); // new hipo4
	hipo::bank BAND_ADC  (factory.getSchema("BAND::adc"  ));
	hipo::bank BAND_TDC  (factory.getSchema("BAND::tdc"  ));
	hipo::bank RUN_config  (factory.getSchema("RUN::config"  ));
        hipo::event readevent;  // new hipo4

	int event_counter = 0;
	while(reader.next()==true){
		if(event_counter%10000==0) cout << "event: " << event_counter << endl;
		if( event_counter > 100000 ) break;
		event_counter++;

                reader.read(readevent); // new hipo4
                readevent.getStructure(BAND_ADC );   // new hipo4
		readevent.getStructure(BAND_TDC );   // new hipo4
		readevent.getStructure(RUN_config );   // new hipo4

		long timestamp = RUN_config.getLong(4,0);
		double phaseCorr = getTriggerPhase(timestamp);


		int nADC = BAND_ADC.getRows();
		int nTDC = BAND_TDC.getRows();
		double ref_time = -100000;


		// Skip events with no entries
		if(nADC==0||nTDC==0) continue;
		// Skip events with less than 50 entries (the laser lights all bars at the same time)
		if(nADC<50||nTDC<50) continue;

		std::map<int,std::vector<int>> ADCInd;
		std::map<int,std::vector<int>> TDCInd;

		for(int aIdx = 0 ; aIdx < nADC ; aIdx++){
			int   ADC_sector    = BAND_ADC.getInt  (0,aIdx);
                        int   ADC_layer     = BAND_ADC.getInt  (1,aIdx);
                        int   ADC_component = BAND_ADC.getInt  (2,aIdx);
                        int   ADC_order     = BAND_ADC.getInt  (3,aIdx);
			int pmtKey = 1000*ADC_sector + 100*ADC_layer + 10*ADC_component + ADC_order;
			
			ADCInd[pmtKey].push_back( aIdx );
		}
		for(int tIdx = 0 ; tIdx < nTDC ; tIdx++){
			int   TDC_sector    = BAND_TDC.getInt  (0,tIdx);
                        int   TDC_layer     = BAND_TDC.getInt  (1,tIdx);
                        int   TDC_component = BAND_TDC.getInt  (2,tIdx);
                        int   TDC_order     = BAND_TDC.getInt  (3,tIdx);
			TDC_order -= 2;
			int pmtKey = 1000*TDC_sector + 100*TDC_layer + 10*TDC_component + TDC_order;
				
			TDCInd[pmtKey].push_back( tIdx );

			// Get reference time
			// 	Reference is PMT in bar 116B_R (layer = 1, sector = 4, comp = 6, order = 3)
			if( pmtKey == 4161 ){
				float TDC_tdc       = (float)(BAND_TDC.getInt(4,tIdx));
				float TDC_time      = TDC_tdc*0.02345 - phaseCorr;
				if( ref_time > 0 && ref_time > TDC_time) ref_time = TDC_time;
				else if( ref_time < 0 ) ref_time = TDC_time;
			}
		}
		if(ref_time == -100000) continue; // i.e. no reference fired

		// ADC is always larger than TDC map:
		std::map<int,std::vector<int>>::iterator iter;
		for( iter = ADCInd.begin() ; iter != ADCInd.end() ; ++iter ){
			// Grab which pmt we are on:
			int pmtKey = iter->first;
			int order = (pmtKey % 10);
			int component = (pmtKey - order)/10 % 10;
			int layer = (pmtKey - order - component*10)/100 % 10;
			int sector  = (pmtKey - order - component*10 - layer*100)/1000 % 10;
			int barKey = sector * 100 + layer * 10 + component;
			
			// Get the ADC inds and TDC inds
			std::vector<int> adcs = iter->second;
			std::vector<int> tdcs = TDCInd[pmtKey];
			double tdiff,adc,fadctdiff;
			
			// More than 1 hit in BOTH ADC and TDC is highly unlikely (by looking  at data) so just throw those events
			// away, and then matching is easy...
			if( tdcs.size() > 1 && adcs.size() > 1 ) continue; // throw away events where I have more than 1 hit in TDC AND ADC
			else if( tdcs.size() == 1 && adcs.size() == 1 ){
				int aIdx = adcs.at(0);
				adc = (double)(BAND_ADC.getInt(4,aIdx));
				fadctdiff = (double)BAND_ADC.getFloat(6,aIdx) - ref_time;

				int tIdx = tdcs.at(0);
				tdiff = ((BAND_TDC.getInt(4,tIdx)) * 0.02345 - phaseCorr) - ref_time;
			}
			else if( tdcs.size() > 1 ){ // need to do time-matching:
				int aIdx = adcs.at(0);
				adc = (double)(BAND_ADC.getInt(4,aIdx));
				fadctdiff = (double)BAND_ADC.getFloat(6,aIdx) - ref_time;

				double adctime = (double)BAND_ADC.getFloat(6,aIdx);
				double prevdiff = 1e10;
				for( int i = 0; i < tdcs.size() ; i++){
					int tIdx = tdcs.at(i);
					double tdctime = (BAND_TDC.getInt(4,tIdx)) * 0.02345 - phaseCorr;
					double thisdiff = fabs(adctime-tdctime);
					if( thisdiff < prevdiff){
						prevdiff = thisdiff;
						tdiff = tdctime - ref_time;
					}
				}
			}
			else{ continue; }

			if( order == 0 ){
				h2_tdc_adc_L[barKey] -> SetTitle(Form("Left PMTs, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
				h2_tdc_adc_L[barKey] -> Fill(adc,tdiff);
				h2_ftdc_adc_L[barKey] -> SetTitle(Form("Left PMTs, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
				h2_ftdc_adc_L[barKey] -> Fill(adc,fadctdiff);
			}
			if( order ==1 ){
				h2_tdc_adc_R[barKey] -> SetTitle(Form("Right PMTs, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
				h2_tdc_adc_R[barKey] -> Fill(adc,tdiff);
				h2_ftdc_adc_R[barKey] -> SetTitle(Form("Right PMTs, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
				h2_ftdc_adc_R[barKey] -> Fill(adc,fadctdiff);
			}

		}

	} // end loop over events

	// Printing results onto canvases
	TCanvas * c0 = new TCanvas("c0","c0",700,900);
	TCanvas * c1 = new TCanvas("c1","c0",700,900);
	c0 -> Modified();
	c0 -> Update(); 
	c1 -> Modified();
	c1 -> Update(); 
	TCanvas *** cSLC_tdc = new TCanvas**[5];
	TCanvas *** cSLC_ftdc = new TCanvas**[5];
	c0 -> Print("tdc_results_timewalk_corr_perPMT.pdf(");
	c0 -> Print("ftdc_results_timewalk_corr_perPMT.pdf(");
	double par1_L[600] = {0};
	double par2_L[600] = {0};
	double par1_R[600] = {0};
	double par2_R[600] = {0};
	double par1e_L[600] = {0};
	double par2e_L[600] = {0};
	double par1e_R[600] = {0};
	double par2e_R[600] = {0};
	for(int is = 0 ; is < 5 ; is++){
		cSLC_tdc[is] = new TCanvas*[6];
		cSLC_ftdc[is] = new TCanvas*[6];
		for(int il = 0 ; il < 6 ; il++){
			cSLC_tdc[is][il] = new TCanvas(Form("TDC_S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),700,900);
			cSLC_ftdc[is][il] = new TCanvas(Form("FADC_S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),700,900);
			cSLC_tdc[is][il] -> Divide(2,7);
			cSLC_ftdc[is][il] -> Divide(2,7);

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);
				bool notEmpty = ( h2_ftdc_adc_L[identifier]->Integral() || h2_tdc_adc_L[identifier]->Integral() );
				if(notEmpty){
					cout << "\tDoing fit for L " << identifier << "\n";
					int midbin;
					double val, min_val, max_val;

					cSLC_tdc[is][il] -> cd(2*cIdx+1);
					gPad -> SetLeftMargin(0.16);
					gPad -> SetBottomMargin(0.20);
					gPad -> SetTopMargin(0);
					midbin = h2_tdc_adc_L[identifier]->GetXaxis()->FindBin(10000);
					val = h2_tdc_adc_L[identifier]->ProjectionY("trash",midbin,midbin+1)->GetMean();
					min_val = val - 2;
					max_val = val + 8;
					h2_tdc_adc_L[identifier] -> GetYaxis() -> SetRangeUser(min_val,max_val);
					h2_tdc_adc_L[identifier] -> Draw("COLZ");
					// Now do the fit:
					double par1 = 0, par2 = 0;
					double par1e = 0, par2e = 0;
					fitTW( h2_tdc_adc_L[identifier] , par1 , par2 , par1e, par2e );
					par1_L[identifier] = par1;
					par2_L[identifier] = par2;
					par1e_L[identifier] = par1e;
					par2e_L[identifier] = par2e;
					//tabL << (is+1) << "\t" << (il+1) << "\t" << (cIdx+1) << "\t" << par1 << "\t" <<par2 << "\t" << par1e  << "\t" << par2e << "\n";


					cSLC_ftdc[is][il] -> cd(2*cIdx+1);
					gPad -> SetLeftMargin(0.16);
					gPad -> SetBottomMargin(0.20);
					gPad -> SetTopMargin(0);
					midbin = h2_ftdc_adc_L[identifier]->GetXaxis()->FindBin(10000);
					val = h2_ftdc_adc_L[identifier]->ProjectionY("trash",midbin,midbin+1)->GetMean();
					min_val = val - 2;
					max_val = val + 8;
					h2_ftdc_adc_L[identifier] -> GetYaxis() -> SetRangeUser(min_val,max_val);
					h2_ftdc_adc_L[identifier] -> Draw("COLZ");
					// No fit for FADC 
				}
				notEmpty = (h2_ftdc_adc_R[identifier]->Integral() || h2_tdc_adc_R[identifier]->Integral() );
				if(notEmpty){
					cout << "\tDoing fit for R " << identifier << "\n";
					int midbin;
					double val, min_val, max_val;

					cSLC_tdc[is][il] -> cd(2*cIdx+2);
					gPad -> SetLeftMargin(0.16);
					gPad -> SetBottomMargin(0.20);
					gPad -> SetTopMargin(0);
					midbin = h2_tdc_adc_R[identifier]->GetXaxis()->FindBin(10000);
					val = h2_tdc_adc_R[identifier]->ProjectionY("trash",midbin,midbin+1)->GetMean();
					min_val = val - 2;
					max_val = val + 8;
					h2_tdc_adc_R[identifier] -> GetYaxis() -> SetRangeUser(min_val,max_val);
					h2_tdc_adc_R[identifier] -> Draw("COLZ");
					// Now do the fit:
					double par1 = 0, par2 = 0;
					double par1e = 0, par2e = 0;
					if( identifier != 416) fitTW( h2_tdc_adc_R[identifier] , par1 , par2 , par1e, par2e );
					par1_R[identifier] = par1;
					par2_R[identifier] = par2;
					par1e_R[identifier] = par1e;
					par2e_R[identifier] = par2e;
					//tabR << (is+1) << "\t" << (il+1) << "\t" << (cIdx+1) << "\t" << par1 << "\t" <<par2 << "\t" << par1e  << "\t" << par2e << "\n";
					cSLC_ftdc[is][il] -> cd(2*cIdx+2);
					gPad -> SetLeftMargin(0.16);
					gPad -> SetBottomMargin(0.20);
					gPad -> SetTopMargin(0);
					midbin = h2_ftdc_adc_R[identifier]->GetXaxis()->FindBin(10000);
					val = h2_ftdc_adc_R[identifier]->ProjectionY("trash",midbin,midbin+1)->GetMean();
					min_val = val - 2;
					max_val = val + 8;
					h2_ftdc_adc_R[identifier] -> GetYaxis() -> SetRangeUser(min_val,max_val);
					h2_ftdc_adc_R[identifier] -> Draw("COLZ");
					// No for for FADC
				}
			}
			cSLC_tdc[is][il] -> Modified();     
			cSLC_tdc[is][il] -> Update();
			cSLC_tdc[is][il] -> Print("tdc_results_timewalk_corr_perPMT.pdf");
			cSLC_ftdc[is][il] -> Modified();     
			cSLC_ftdc[is][il] -> Update();
			cSLC_ftdc[is][il] -> Print("ftdc_results_timewalk_corr_perPMT.pdf");
		}
	}
	c0 -> Print("tdc_results_timewalk_corr_perPMT.pdf)");
	c1 -> Print("ftdc_results_timewalk_corr_perPMT.pdf)");

	// -------------------------------------------------------------------------------------------------
	// Saving fit values to ccdb tables
	ofstream tabL, tabR;
	tabL.open("timeWalkPar_L.txt");
	tabR.open("timeWalkPar_R.txt");
	for(int is = 1 ; is <= 5 ; is++){
		for(int il = 1 ; il <= 6 ; il++){
			for(int ic = 1 ; ic <= slc[il-1][is-1] ; ic++){
				int idx = 100*is + 10*il + ic;
				tabL << is << "\t" << il << "\t" << ic << "\t" 	<< par2_L[idx] << "\t" << par1_L[idx] << "\t" 
										<< par2e_L[idx]  << "\t" << par1e_L[idx] << "\n";
				tabR << is << "\t" << il << "\t" << ic << "\t" 	<< par2_R[idx] << "\t" << par1_R[idx] << "\t" 
										<< par2e_R[idx]  << "\t" << par1e_R[idx] << "\n";
			}
		}
	}
	tabL.close();
	tabR.close();

	return 0;
}
// ========================================================================================================================================
void PrettyTH2F(TH2F * h2) {
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
}
// ========================================================================================================================================
void PrettyTH1F(TH1F * h1, int color) {
	h1 -> GetXaxis() -> CenterTitle();
	h1 -> GetYaxis() -> CenterTitle();

	h1 -> SetLineColor(color);
	h1 -> SetLineWidth(2);

	h1 -> SetTitleSize(0.5);

	h1 -> GetXaxis() -> SetLabelSize(0.09);
	h1 -> GetYaxis() -> SetLabelSize(0.09);
	h1 -> GetXaxis() -> SetTitleSize(0.09);
	h1 -> GetYaxis() -> SetTitleSize(0.09);

	h1 -> GetYaxis() -> SetTitleOffset(0.60);

	h1 -> GetYaxis() -> SetNdivisions(109);
	h1 -> GetXaxis() -> SetNdivisions(107);
}
double getTriggerPhase( long timeStamp ) {
	double tPh = 0.;

	long Period = 4.0;
	long Cycles = 6.0;
	long Phase  = 3.0;

	if( timeStamp != -1 ) 
		tPh = (double)(Period *( (timeStamp + Phase) % Cycles ));

	return tPh;
}

void fitTW(TH2F * hist , double &par1, double &par2 , double &par1_err, double &par2_err ){
	std::vector<double> xs;
	std::vector<double> ys;
	std::vector<double> xErrs;
	std::vector<double> yErrs;
	std::vector<double> res;
	std::vector<double> resErrs;
	walkCorr( &xs, &xErrs, &ys, &yErrs, &res, &resErrs, 4000, hist );

	int dim = xs.size();
	TGraphErrors *g = new TGraphErrors(dim, &xs[0], &ys[0], &xErrs[0], &yErrs[0]);
	TF1 * model = new TF1("timeWalk",wlk,0,20000,2);
	TFitResultPtr ptr = g->Fit(model,"QES");
	par1 = ptr->Parameter(0);
	par1_err = ptr->ParError(0);
	par2 = ptr->Parameter(1);
	par2_err = ptr->ParError(1);
	return;
}

void walkCorr(	std::vector<double> *adcs		,
		std::vector<double> *adcsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		double widthCut				,
		TH2F * hist				){
	int currBin = 0;
	while( currBin < 800 ){
		double xPt, yPt, yEr, ySig, ySigEr;
		int step = doProj( hist , currBin , false , 0, xPt, yPt, yEr, ySig, ySigEr );
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
int doProj( TH2F * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE){
	int cnt = 0;
	int step = 0;
	TCanvas * trash = new TCanvas("trash");
	//int thres = hist->GetEntries() / 50;
	int thres = 2000;
	while ( cnt < thres ){
		char temp[100];
		sprintf(temp,"slice_%d_%d",flag,bin);
		TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
		cnt = (int) pj->Integral();
		delete pj;
		if( (bin+step) >= 800) break;
		step+=1;
	}
	char temp[100];
	sprintf(temp,"slice_%d_%d",flag,bin);
	TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
	
	// Getting the mean x value in this range:
	hist->GetXaxis()->SetRange(bin,bin+step);
	x = hist->GetMean(1);
	hist->GetXaxis()->SetRange();
	TFitResultPtr f = pj->Fit("gaus","QES","",-1200,300);
	y = f->Parameter(1);
	yE = f->ParError(1);
	sig = f->Parameter(2);
	sigE = f->ParError(2);
	
	if( write ){
		pj->Write();
	}

	delete trash;
	delete pj;

	return step;
}


double wlk( double *x , double *p){
	double var = *x;
	return p[0] + p[1] / sqrt(var);
}

*/
