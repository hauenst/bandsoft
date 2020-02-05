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

void attenFit( TH2F * hist , TCanvas * c , int cd , int is , int il, int ic, double &mu , double &mu_e);
void attenCorr(	std::vector<double> *adcs		,
		std::vector<double> *adcsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		TH2F * hist				,
		double xposcut				);
void reflectFit( TH2F * hist , TCanvas * c , int cd , int is , int il, int ic, double mu_guess , double &mu , double &mu_e);
double lnR_reflect_short( double *x , double *p);
double lnR_reflect_med( double *x , double *p);
double lnR_reflect_long( double *x , double *p);

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
	TH2F ** h2_tdc_tdiff_adc = new TH2F * [nHistos];
	TH2F ** h2_ftdc_tdiff_adc = new TH2F * [nHistos];
	TH2F ** h2_tdc_tdiff_amp = new TH2F * [nHistos];
	TH2F ** h2_ftdc_tdiff_amp = new TH2F * [nHistos];
	for(int i = 0 ; i < nHistos ; i++){
		h2_tdc_tdiff_adc[i] = new TH2F(Form("h2_tdc_tdiff_adc_%i",i),"	;x_{TDC}  [cm];ln(ADC_{L}/ADC_{R})",		500,-150,150,400,-10,10);
		h2_tdc_tdiff_amp[i] = new TH2F(Form("h2_tdc_tdiff_amp_%i",i),"	;x_{TDC}  [cm];ln(AMP_{L}/AMP_{R})",		500,-150,150,400,-10,10);

		h2_ftdc_tdiff_adc[i] = new TH2F(Form("h2_ftdc_tdiff_adc_%i",i),";x_{FTDC} [cm];ln(ADC_{L}/ADC_{R});",		500,-150,150,400,-10,10);
		h2_ftdc_tdiff_amp[i] = new TH2F(Form("h2_ftdc_tdiff_amp_%i",i),";x_{FTDC} [cm];ln(AMP_{L}/AMP_{R});",		500,-150,150,400,-10,10);
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
		// Create Bars from PMTs
		BAND->CreateBars();
		// Prepare maps to access information from bars
		BAND->InitGettersBars();
		
		// Access information we need from bars
		std::map<int,double> lnR_adc = BAND->GetBars_FADC_lnAdc();
		std::map<int,double> lnR_amp = BAND->GetBars_FADC_lnAmp();
		std::map<int,double> tdc_xpos = BAND->GetBars_TDC_xPos();
		std::map<int,double> ftdc_xpos = BAND->GetBars_FADC_xPos();

		// Loop through maps and fill our histograms:
		std::map<int,double>::iterator iter_bar;
		for( iter_bar = tdc_xpos.begin() ; iter_bar != tdc_xpos.end() ; iter_bar++ ){
			int barKey = iter_bar->first;
			int component = (barKey % 10);
			int layer = (barKey - component)/10 % 10;
			int sector = (barKey - component - layer*10)/100 % 10;

			h2_tdc_tdiff_adc[barKey] -> Fill( 	tdc_xpos[barKey],	lnR_adc[barKey]	);
			h2_ftdc_tdiff_adc[barKey] -> Fill( 	ftdc_xpos[barKey],	lnR_adc[barKey] );
			h2_tdc_tdiff_amp[barKey] -> Fill( 	tdc_xpos[barKey],	lnR_amp[barKey]	);
			h2_ftdc_tdiff_amp[barKey] -> Fill( 	ftdc_xpos[barKey],	lnR_amp[barKey] );
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
	
	for(int is = 0 ; is < 5 ; is++){
		cSLC_tdc[is] = new TCanvas*[6];
		cSLC_ftdc[is] = new TCanvas*[6];
		for(int il = 0 ; il < 6 ; il++){
			cSLC_tdc[is][il] = new TCanvas(Form("TDC_S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),700,900);
			cSLC_ftdc[is][il] = new TCanvas(Form("FADC_S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),700,900);
			cSLC_tdc[is][il] -> Divide(4,7);
			cSLC_ftdc[is][il] -> Divide(4,7);

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);

				if( h2_tdc_tdiff_amp[identifier]->Integral()  ){
					double mu = 0, mu_e = 0;
					zoomTH2F( h2_tdc_tdiff_amp[identifier], 0, 1, 1 ); drawTH2F( h2_tdc_tdiff_amp[identifier] , cSLC_tdc[is][il] , 4*cIdx+1 );
					attenFit( h2_tdc_tdiff_amp[identifier] , cSLC_tdc[is][il] , 4*cIdx+2 , is, il, cIdx , mu , mu_e );

					// Now we could look at the adc plot and try to extract the reflection coefficient
					zoomTH2F( h2_tdc_tdiff_adc[identifier], 0, 1, 1 ); drawTH2F( h2_tdc_tdiff_adc[identifier] , cSLC_tdc[is][il] , 4*cIdx+3 );
					double n = 0, n2 = 0;
					reflectFit( h2_tdc_tdiff_adc[identifier] , cSLC_tdc[is][il] , 4*cIdx+4 , is, il, cIdx , mu, n , n2 );
				}
				if( h2_ftdc_tdiff_amp[identifier]->Integral()  ){
					double mu = 0, mu_e = 0;
					zoomTH2F( h2_ftdc_tdiff_amp[identifier], 0, 1, 1 ); drawTH2F( h2_ftdc_tdiff_amp[identifier] , cSLC_ftdc[is][il] , 4*cIdx+1 );
					attenFit( h2_ftdc_tdiff_amp[identifier] , cSLC_ftdc[is][il] , 4*cIdx+2 , is, il, cIdx , mu , mu_e );

					// Now we could look at the adc plot and try to extract the reflection coefficient
					zoomTH2F( h2_ftdc_tdiff_adc[identifier], 0, 1, 1 ); drawTH2F( h2_ftdc_tdiff_adc[identifier] , cSLC_ftdc[is][il] , 4*cIdx+3 );
					double n = 0, n2 = 0;
					reflectFit( h2_ftdc_tdiff_adc[identifier] , cSLC_ftdc[is][il] , 4*cIdx+4 , is, il, cIdx , mu, n , n2 );
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

void attenFit( TH2F * hist , TCanvas * c , int cd , int is , int il, int ic, double &mu , double &mu_e){
	std::vector<double> xs;
	std::vector<double> ys;
	std::vector<double> xErrs;
	std::vector<double> yErrs;
	std::vector<double> res;
	std::vector<double> resErrs;

	double zoom = 25;
	if( (is+1)==3 || (is+1)==4 ) zoom = 12;
	attenCorr( &xs, &xErrs, &ys, &yErrs, &res, &resErrs, hist , zoom );

	int dim = xs.size();
	TGraphErrors *g = new TGraphErrors(dim, &xs[0], &ys[0], &xErrs[0], &yErrs[0]);

	TFitResultPtr ptr = g->Fit("pol1","QESR","",-zoom,zoom);
	mu = 2./ptr->Parameter(1);
	mu_e = 0.;
	char string[100];
	sprintf(string, "Mu: %.2g +/- %.2g", mu,mu_e);

	g->SetTitle(string);
	c->cd(cd);
	g->SetMarkerStyle(20);
	g->Draw("AP");	
	return;
}

void attenCorr(	std::vector<double> *adcs		,
		std::vector<double> *adcsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		TH2F * hist				,
		double xposcut				){
	int currBin = 0;
	while( currBin < 500 ){
		double xPt, yPt, yEr, ySig, ySigEr;
		int step = doProj( hist , currBin , false , 0, xPt, yPt, yEr, ySig, ySigEr , 500, 500 );
		currBin += step ;
		if( fabs(xPt) > xposcut ) continue;
		adcs		->push_back(xPt);
		adcsErr		->push_back(0);
		times		->push_back(yPt);
		timesErr	->push_back(yEr);
		res		->push_back(ySig);
		resErr		->push_back(ySigEr);
	}
}

void reflectFit( TH2F * hist , TCanvas * c , int cd , int is , int il, int ic, double mu_guess , double &mu , double &mu_e){
	std::vector<double> xs;
	std::vector<double> ys;
	std::vector<double> xErrs;
	std::vector<double> yErrs;
	std::vector<double> res;
	std::vector<double> resErrs;
	
	double zoom = 120;
	if( (is+1)==3 || (is+1)==4 ) zoom = 40;
	attenCorr( &xs, &xErrs, &ys, &yErrs, &res, &resErrs, hist , zoom );

	int dim = xs.size();
	TGraphErrors *g = new TGraphErrors(dim, &xs[0], &ys[0], &xErrs[0], &yErrs[0]);

	TF1 * model = NULL;
	if( (is+1)==3 || (is+1)==4 ) 		model = new TF1("atten_reflect",lnR_reflect_short,-150,150,3);
	else if( (is+1)==1 )	 		model = new TF1("atten_reflect",lnR_reflect_med,-150,150,3);
	else if( (is+1)==2 || (is+1)==5) 	model = new TF1("atten_reflect",lnR_reflect_long,-150,150,3);
	if( (is+1)==3 || (is+1)==4 ) mu_guess *= (-1);
	model->SetParameter(0,mu_guess);
	//model->SetParLimits(0,mu_guess*0.8,mu_guess*1.2);
	model->SetParameter(1,0.2);
	model->SetParLimits(1,0,1);
	model->SetParameter(2,0);
	//model->SetParLimits(2,-0.2,0.2);
	//model->SetParameter(2,0);
	cout << (is+1) << " " << (il+1) << " " << (ic+1) << "\n";
	TFitResultPtr ptr = g->Fit(model,"QESR","",-zoom/2.,zoom/2.);
	TMatrixDSym cov = ptr->GetCovarianceMatrix();
	//cout << ptr->Parameter(0) << " " << ptr->Parameter(1) << " " << ptr->Parameter(2) << "\n";
	cout << "Mu: " 		<< ptr->Parameter(0) << " +/- " << ptr->ParError(0) << "\n";
	cout << "Alpha: " 	<< ptr->Parameter(1) << " +/- " << ptr->ParError(1) << "\n";
	cout << "Offset: " 	<< ptr->Parameter(2) << " +/- " << ptr->ParError(2) << "\n\n";
		
	c->cd(cd);
	g->SetMarkerStyle(20);
	g->Draw("AP");	
	return;
}
double lnR_reflect_med( double *x , double *p){
	double var = (*x);
	double mu = p[0];
	double a = p[1];
	// p[0] = mu
	// p[1] = alpha
	double L = BARLENGTHS[0];
	
	return p[2] + 2.*var/mu + log(1. + a*exp( -(2*var+L)/mu )) - log(1. + a*exp( -(-2*var+L)/mu ));
	//return p[2] + 2.*var/mu + log(1. + a*exp( -(2*var+L)/mu ) + a*a*exp(-2*L/mu) + a*a*a*exp( -(2*var+3*L)/mu ) ) - log(1. + a*exp( -(-2*var+L)/mu ) + a*a*exp(-2.*L/mu) + a*a*a*exp( -(-2.*var+3*L)/mu  ) );
	//double A = exp(L/mu)*(a-1) - exp(2.*var/mu)*a;
	//double B = exp( (L+2*var)/mu )*(a-1) - a;
	//return log( A / B );
}
double lnR_reflect_long( double *x , double *p){
	double var = -(*x);
	double mu = p[0];
	double a = p[1];
	// p[0] = mu
	// p[1] = alpha
	double L = BARLENGTHS[1];
	return p[2] + 2.*var/mu + log(1. + a*exp( -(2*var+L)/mu )) - log(1. + a*exp( -(-2*var+L)/mu ));
	//return p[2] + 2.*var/mu + log(1. + a*exp( -(2*var+L)/mu ) + a*a*exp(-2*L/mu) + a*a*a*exp( -(2*var+3*L)/mu ) ) - log(1. + a*exp( -(-2*var+L)/mu ) + a*a*exp(-2.*L/mu) + a*a*a*exp( -(-2.*var+3*L)/mu  ) );
	//double A = exp(L/mu)*(a-1) - exp(2.*var/mu)*a;
	//double B = exp( (L+2*var)/mu )*(a-1) - a;
	//return log( A / B );
}
double lnR_reflect_short( double *x , double *p){
	double var = -(*x);
	double mu = p[0];
	double a = p[1];
	// p[0] = mu
	// p[1] = alpha
	double L = BARLENGTHS[2];
	return p[2] + 2.*var/mu + log(1. + a*exp( -(2*var+L)/mu )) - log(1. + a*exp( -(-2*var+L)/mu ));
	//return p[2] + 2.*var/mu + log(1. + a*exp( -(2*var+L)/mu ) + a*a*exp(-2*L/mu) + a*a*a*exp( -(2*var+3*L)/mu ) ) - log(1. + a*exp( -(-2*var+L)/mu ) + a*a*exp(-2.*L/mu) + a*a*a*exp( -(-2.*var+3*L)/mu  ) );
	//double A = exp(L/mu)*(a-1) - exp(2.*var/mu)*a;
	//double B = exp( (L+2*var)/mu )*(a-1) - a;
	//return log( A / B );
}
