
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
void fitTW(TH2F * hist , TCanvas * c, int cd, int s, int l, int co, int o, double cut, double &par1, double &par2 , double &par1_err, double &par2_err );
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
	BAND->LoadTimeWalk();
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
	TH2F ** h2_tdc_amp_L = new TH2F * [nHistos];
	TH2F ** h2_tdc_amp_R = new TH2F * [nHistos];
	//TH2F ** h2_ftdc_adc_L = new TH2F * [nHistos];
	//TH2F ** h2_ftdc_adc_R = new TH2F * [nHistos];
	for(int i = 0 ; i < nHistos ; i++){
		h2_tdc_adc_L[i] = new TH2F(Form("h2_tdc_adc_L_%i",i),";ADC_{L};t_{TDC,L} - ref [ns]",800,0,20000,3000,-150,150);
		h2_tdc_adc_R[i] = new TH2F(Form("h2_tdc_adc_R_%i",i),";ADC_{R};t_{TDC,R} - ref [ns]",800,0,20000,3000,-150,150);
		h2_tdc_amp_L[i] = new TH2F(Form("h2_tdc_amp_L_%i",i),";AMP_{L};t_{TDC,L} - ref [ns]",600,0,4200,3000,-150,150);
		h2_tdc_amp_R[i] = new TH2F(Form("h2_tdc_amp_R_%i",i),";AMP_{R};t_{TDC,R} - ref [ns]",600,0,4200,3000,-150,150);
		//h2_ftdc_adc_L[i] = new TH2F(Form("h2_ftdc_adc_L_%i",i),";ADC_{L};t_{FADC,L} - ref [ns]",800,0,20000,3000,0,300);
		//h2_ftdc_adc_R[i] = new TH2F(Form("h2_ftdc_adc_R_%i",i),";ADC_{R};t_{FADC,R} - ref [ns]",800,0,20000,3000,0,300);
	}

	// Loop over all events in file
	int event_counter = 0;
	while(reader.next()==true){
		if(event_counter%10000==0) cout << "event: " << event_counter << endl;
		//if( event_counter > 100000 ) break;
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
			// Skip non-laser entries:
		if(nADC<90||nTDC<90) continue;

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
				h2_tdc_amp_L[barKey] -> Fill(	amp[pmtKey]	,	tdc[pmtKey] - tdc[ref_photodiode] 	);
				//h2_ftdc_adc_L[barKey] -> Fill(	adc[pmtKey]	,	ftdc[pmtKey] - ftdc[ref_photodiode] 	);
			}
			else if( order == 1){
				h2_tdc_adc_R[barKey] -> Fill(	adc[pmtKey]	,	tdc[pmtKey] - tdc[ref_photodiode] 	);
				h2_tdc_amp_R[barKey] -> Fill(	amp[pmtKey]	,	tdc[pmtKey] - tdc[ref_photodiode] 	);
				//h2_ftdc_adc_R[barKey] -> Fill(	adc[pmtKey]	,	ftdc[pmtKey] - ftdc[ref_photodiode] 	);
			}
		}

	} // end loop over events

	TString outputName = argv[2];
	TCanvas * c0 = new TCanvas("c0","c0",700,900);
	//TCanvas * c1 = new TCanvas("c1","c0",700,900);
	TCanvas * c2 = new TCanvas("c2","c0",700,900);
	c0 -> Modified();
	c0 -> Update(); 
	//c1 -> Modified();
	//c1 -> Update(); 
	c2 -> Modified();
	c2 -> Update(); 
	TCanvas *** cSLC_adc = new TCanvas**[5];
	//TCanvas *** cSLC_ftdc = new TCanvas**[5];
	TCanvas *** cSLC_amp = new TCanvas**[5];
	c0 -> Print(Form("adc_%s(" ,outputName.Data()));
	//c1 -> Print(Form("ftdc_%s(",outputName.Data()));
	c2 -> Print(Form("amp_%s(",outputName.Data()));
	double par1_L[600][2] 	= {0};
	double par2_L[600][2] 	= {0};
	double par1_R[600][2] 	= {0};
	double par2_R[600][2] 	= {0};
	double par1e_L[600][2] 	= {0};
	double par2e_L[600][2] 	= {0};
	double par1e_R[600][2] 	= {0};
	double par2e_R[600][2] 	= {0};
	
	for(int is = 0 ; is < 5 ; is++){
		cSLC_adc[is] = new TCanvas*[6];
		cSLC_amp[is] = new TCanvas*[6];
		//cSLC_ftdc[is] = new TCanvas*[6];
		for(int il = 0 ; il < 6 ; il++){
			cSLC_adc[is][il] = new TCanvas(Form("TDC_S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),700,900);
			cSLC_amp[is][il] = new TCanvas(Form("AMP_S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),700,900);
			//cSLC_ftdc[is][il] = new TCanvas(Form("FADC_S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),700,900);
			cSLC_adc[is][il] -> Divide(4,7);
			cSLC_amp[is][il] -> Divide(4,7);
			//cSLC_ftdc[is][il] -> Divide(2,7);

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);

				if( h2_tdc_adc_L[identifier]->Integral()  ){
					double par1 = 0, par2 = 0, par1e = 0, par2e = 0;
					//zoomTH2F( h2_tdc_adc_L[identifier] , 10000 , 4, 8 );
					zoomTH2F( h2_tdc_amp_L[identifier] , 1000 , 4, 8 );
					//drawTH2F( h2_tdc_adc_L[identifier] , cSLC_adc[is][il] ,  4*cIdx+1 );
					drawTH2F( h2_tdc_amp_L[identifier] , cSLC_amp[is][il] ,  4*cIdx+1 );
					cout << "\tfitting SLCO: " << (identifier*10 + 0) << "\n";
					//fitTW( 	  h2_tdc_adc_L[identifier] , cSLC_adc[is][il] ,  4*cIdx+2 , (is+1), (il+1), (cIdx+1), 0, 4000,par1 , par2 , par1e, par2e );
					// Fit results:
					par1_L[identifier][0] = par1;
					par2_L[identifier][0] = par2;
					par1e_L[identifier][0] = par1e;
					par2e_L[identifier][0] = par2e;
					par1 = 0, par2 = 0, par1e = 0, par2e = 0;
					fitTW( 	  h2_tdc_amp_L[identifier] , cSLC_amp[is][il] ,  4*cIdx+2 , (is+1), (il+1), (cIdx+1), 0, 500,par1 , par2 , par1e, par2e );
					par1_L[identifier][1] = par1;
					par2_L[identifier][1] = par2;
					par1e_L[identifier][1] = par1e;
					par2e_L[identifier][1] = par2e;
				}
				if( h2_tdc_adc_R[identifier]->Integral()  ){
					double par1 = 0, par2 = 0, par1e = 0, par2e = 0;
					//zoomTH2F( h2_tdc_adc_R[identifier] , 10000 , 4, 8 );
					zoomTH2F( h2_tdc_amp_R[identifier] , 1000 , 4, 8 );
					//drawTH2F( h2_tdc_adc_R[identifier] , cSLC_adc[is][il] ,  4*cIdx+3 );
					drawTH2F( h2_tdc_amp_R[identifier] , cSLC_amp[is][il] ,  4*cIdx+3 );
					cout << "\tfitting SLCO: " << (identifier*10 + 1) << "\n";
					//if( identifier!= 416 )
					//	fitTW( 	  h2_tdc_adc_R[identifier] , cSLC_adc[is][il] ,  4*cIdx+4 , (is+1), (il+1), (cIdx+1), 1, 4000,par1 , par2 , par1e, par2e );
					par1 = 0, par2 = 0, par1e = 0, par2e = 0;
					// Fit results:
					par1_R[identifier][1] = par1;
					par2_R[identifier][1] = par2;
					par1e_R[identifier][1] = par1e;
					par2e_R[identifier][1] = par2e;
					if( identifier!= 416 )
						fitTW( 	  h2_tdc_amp_R[identifier] , cSLC_amp[is][il] ,  4*cIdx+4 , (is+1), (il+1), (cIdx+1), 1, 500,par1 , par2 , par1e, par2e );
				}
				/*
				if( h2_ftdc_adc_L[identifier]->Integral()  ){
					zoomTH2F( h2_ftdc_adc_L[identifier] , 10000 , 4, 4 );
					drawTH2F( h2_ftdc_adc_L[identifier] , cSLC_ftdc[is][il] ,  2*cIdx+1 );
				}
				if( h2_ftdc_adc_R[identifier]->Integral()  ){
					zoomTH2F( h2_ftdc_adc_R[identifier] , 10000 , 4, 4 );
					drawTH2F( h2_ftdc_adc_R[identifier] , cSLC_ftdc[is][il] ,  2*cIdx+2 );
				}
				*/
			}
			cSLC_adc[is][il] -> Modified();     
			cSLC_adc[is][il] -> Update();
			cSLC_adc[is][il] -> Print(Form("adc_%s",outputName.Data()));
			//cSLC_ftdc[is][il] -> Modified();     
			//cSLC_ftdc[is][il] -> Update();
			//cSLC_ftdc[is][il] -> Print(Form("ftdc_%s",outputName.Data()));
			cSLC_amp[is][il] -> Modified();     
			cSLC_amp[is][il] -> Update();
			cSLC_amp[is][il] -> Print(Form("amp_%s",outputName.Data()));
		}
	}
	c0 -> Print(Form("adc_%s)",outputName.Data()));
	//c1 -> Print(Form("ftdc_%s)",outputName.Data()));
	c2 -> Print(Form("amp_%s)",outputName.Data()));
	
	/*
	// Print out the results
	ofstream tabL, tabR;
	tabL.open("secondIter_timeWalkPar_L.txt");
	tabR.open("secondIter_timeWalkPar_R.txt");
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
	*/

	return 0;
}


void fitTW(TH2F * hist , TCanvas * c, int cd, int s, int l, int co, int o, double cut, double &par1, double &par2 , double &par1_err, double &par2_err ){
	std::vector<double> xs;
	std::vector<double> ys;
	std::vector<double> xErrs;
	std::vector<double> yErrs;
	std::vector<double> res;
	std::vector<double> resErrs;
	walkCorr( &xs, &xErrs, &ys, &yErrs, &res, &resErrs, cut, hist );

	int dim = xs.size();
	TGraphErrors *g = new TGraphErrors(dim, &xs[0], &ys[0], &xErrs[0], &res[0]);
	TF1 * model = new TF1("timeWalk",wlk,0,20000,3);
	model->SetParameter(0,-15);
	model->SetParameter(1,80);
	model->SetParameter(2,0.5);
	TFitResultPtr ptr = g->Fit(model,"QES");
	//par1 = ptr->Parameter(0);
	//par1_err = ptr->ParError(0);
	//par2 = ptr->Parameter(1);
	//par2_err = ptr->ParError(1);

	c->cd(cd);
	gStyle->SetTitleW(0.6);
	g->SetTitle(Form("SLCO: %i %i %i %i",s,l,co,o));
	g->SetMarkerStyle(20);
	g->Draw("AP");

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
		int step = doProj( hist , currBin , false , 0, xPt, yPt, yEr, ySig, ySigEr , 500, 800 );
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
	//return p[0] + p[1] / sqrt(var);
	return p[0] + p[1] / pow(var,p[2]);
}
