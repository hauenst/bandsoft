
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

#include "calibclass.h"


using namespace std;

// Forward-declaring functions
void PrettyTH1F(TH1F * h1, int color);
void PrettyTH2F(TH2F * h2);
void PrettyTGraphErrors(TGraphErrors * gP, int color);
double getTriggerPhase( long timeStamp );
void LoadTimeWalk();
void LoadLROffsets();
void LoadVelocityMap();
double parA_L[600] = {0.};		// loaded
double parB_L[600] = {0.};		// loaded
double parA_R[600] = {0.};		// loaded
double parB_R[600] = {0.};		// loaded
double TDC_TDIFF[600] = {0.};		// loaded
double FADC_TDIFF[600] = {0.};		// loaded
double TDC_VEFF[600] = {0.};		// loaded
double FADC_VEFF[600] = {0.};		// loaded
double BARLENGTHS[]  = {163.7,201.9,51.2,51.2,201.9};
struct PMTHit{
	int ID;
	double adc, amp;
	double tdc=1e10, ftdc=1e10;
	bool processed = false;
};
struct BarHit{
	int ID;
	double adc, amp;
	double tdc_tdiff, ftdc_tdiff;
	double tdc_tmean, ftdc_tmean;
};
void drawHist( TH2F * hist , TCanvas * c , int cd );
void attenFit( TH2F * hist , TCanvas * c , int cd , int is , int il, int ic, double &mu , double &mu_e);

int doProj( TH2F * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE);
void attenCorr(	std::vector<double> *adcs		,
		std::vector<double> *adcsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		TH2F * hist				,
		double xposcut				);
double lnR_reflect_short( double *x , double *p);
double lnR_reflect_med( double *x , double *p);
double lnR_reflect_long( double *x , double *p);
void reflectFit( TH2F * hist , TCanvas * c , int cd , int is , int il, int ic, double mu_guess , double &mu , double &mu_e);

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

int main(int argc, char** argv) {
	if( argc!= 2 ){
		cerr << "Wrong number of arguments. Please use instead:\n\t./code [InputFile]\n";
		return -1;
	}
	gStyle->SetTitleSize(0.2,"t");
	gStyle->SetOptStat(0);

	TString inputFile = argv[1];

	// Load calibration class and parameters
	calibclass * BAND = new calibclass();
	BAND->LoadTimeWalk();
	BAND->LoadLROffsets();
	BAND->LoadVelocityMap();

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

		PrettyTH2F(h2_tdc_tmean_gmadc[i]);
		PrettyTH2F(h2_ftdc_tmean_gmadc[i]);
		PrettyTH2F(h2_tdc_tmean_gmamp[i]);
		PrettyTH2F(h2_ftdc_tmean_gmamp[i]);
	}

	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);
        hipo::dictionary  factory;      
        reader.readDictionary(factory); 
	hipo::bank BAND_ADC  (factory.getSchema("BAND::adc"  ));
	hipo::bank BAND_TDC  (factory.getSchema("BAND::tdc"  ));
	hipo::bank RUN_config  (factory.getSchema("RUN::config"  ));
        hipo::event readevent;

	int event_counter = 0;
	while(reader.next()==true){
		if(event_counter%10000==0) cout << "event: " << event_counter << endl;
		if( event_counter > 200000 ) break;
		event_counter++;

                reader.read(readevent); 
                readevent.getStructure(BAND_ADC );   
		readevent.getStructure(BAND_TDC );   
		readevent.getStructure(RUN_config );

		long timestamp = RUN_config.getLong(4,0);
		double phaseCorr = getTriggerPhase(timestamp);

		int nADC = BAND_ADC.getRows();
		int nTDC = BAND_TDC.getRows();
		// Skip events with no entries
		if(nADC==0||nTDC==0) continue;
		// Skip non-laser events
		if(nADC<90||nTDC<90) continue;

		BAND->CreatePMTs( BAND_ADC, BAND_TDC, phaseCorr );
		BAND->CreateBars();
		BAND->InitGetters();
		
		// Now we can access these bars however we need
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

	
	TCanvas * c0 = new TCanvas("c0","c0",700,900);
	TCanvas * c1 = new TCanvas("c1","c0",700,900);
	c0 -> Modified();
	c0 -> Update(); 
	c1 -> Modified();
	c1 -> Update(); 
	TCanvas *** cSLC_tdc = new TCanvas**[5];
	TCanvas *** cSLC_ftdc = new TCanvas**[5];
	c0 -> Print("tdc_results_barbybar.pdf(");
	c1 -> Print("ftdc_results_barbybar.pdf(");
	
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

				if( h2_tdc_tmean_gmamp[identifier]->Integral()  ){
					drawHist( h2_tdc_tmean_gmadc[identifier] , cSLC_tdc[is][il] , 2*cIdx+1 );
					drawHist( h2_tdc_tmean_gmamp[identifier] , cSLC_tdc[is][il] , 2*cIdx+2 );
				}
				if( h2_ftdc_tmean_gmamp[identifier]->Integral()  ){
					drawHist( h2_ftdc_tmean_gmadc[identifier] , cSLC_ftdc[is][il] , 2*cIdx+1 );
					drawHist( h2_ftdc_tmean_gmamp[identifier] , cSLC_ftdc[is][il] , 2*cIdx+2 );
				}
			}
			cSLC_tdc[is][il] -> Modified();     
			cSLC_tdc[is][il] -> Update();
			cSLC_tdc[is][il] -> Print("tdc_results_barbybar.pdf");
			cSLC_ftdc[is][il] -> Modified();     
			cSLC_ftdc[is][il] -> Update();
			cSLC_ftdc[is][il] -> Print("ftdc_results_barbybar.pdf");
		}
	}
	c0 -> Print("tdc_results_barbybar.pdf)");
	c1 -> Print("ftdc_results_barbybar.pdf)");

	return 0;
}

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


void LoadTimeWalk(){
	ifstream f;
	int sector, layer, component, barId;
	double parA, parB, temp;

	f.open("../../include/time_walk_corr_left.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> parA;
		f >> parB;
		f >> temp;
		f >> temp;
		parA_L[barId] = parA;
		parB_L[barId] = parB;
	}
	f.close();

	f.open("../../include/time_walk_corr_right.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> parA;
		f >> parB;
		f >> temp;
		f >> temp;
		parA_R[barId] = parA;
		parB_R[barId] = parB;
	}
	f.close();
	return;
}
void LoadLROffsets(){
	ifstream f;
	int sector, layer, component, barId;
	double tdc_off, fadc_off, temp;

	f.open("../../include/lr_offsets.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> tdc_off;
		f >> fadc_off;
		f >> temp;
		f >> temp;
		TDC_TDIFF[barId] = tdc_off;
		FADC_TDIFF[barId] = fadc_off;
	}
	f.close();
	return;
}
void LoadVelocityMap(){
	ifstream f;
	int sector, layer, component, barId;
	double veff_tdc, veff_fadc, temp;

	f.open("../../include/effective_velocity.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> veff_tdc;
		f >> veff_fadc;
		f >> temp;
		f >> temp;
		TDC_VEFF[barId] = veff_tdc;
		FADC_VEFF[barId] = veff_fadc;
	}
	f.close();

	return;
}

void drawHist( TH2F * hist , TCanvas * c , int cd ){
	c -> cd(cd);
	gPad -> SetLeftMargin(0.16);
	gPad -> SetBottomMargin(0.20);
	gPad -> SetTopMargin(0);
	int midbin = hist->GetXaxis()->FindBin(10000.);
	double val = hist->ProjectionY("trash",midbin,midbin+1)->GetMean();
	double min_val = val - 1.5;
	double max_val = val + 1.5;
	hist->GetYaxis() -> SetRangeUser(min_val,max_val);
	hist->Draw("COLZ");
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
		int step = doProj( hist , currBin , false , 0, xPt, yPt, yEr, ySig, ySigEr );
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
int doProj( TH2F * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE){
	int cnt = 0;
	int step = 0;
	TCanvas * trash = new TCanvas("trash");
	//int thres = hist->GetEntries() / 50;
	int thres = 500;
	while ( cnt < thres ){
		char temp[100];
		sprintf(temp,"slice_%d_%d",flag,bin);
		TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
		cnt = (int) pj->Integral();
		delete pj;
		if( (bin+step) >= 500) break;
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
		TFitResultPtr f = pj->Fit("gaus","QESR","",-10,10);
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
	model->FixParameter(0,mu_guess);
	//model->SetParLimits(0,mu_guess*0.8,mu_guess*1.2);
	model->SetParameter(1,0.2);
	model->SetParLimits(1,0,1);
	model->SetParameter(2,0);
	model->SetParLimits(2,-0.2,0.2);
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
	
	return p[2] + 2.*var/mu + log(1. + a*exp( -(2*var+L)/mu ) + a*a*exp(-2*L/mu) + a*a*a*exp( -(2*var+3*L)/mu ) ) - log(1. + a*exp( -(-2*var+L)/mu ) + a*a*exp(-2.*L/mu) + a*a*a*exp( -(-2.*var+3*L)/mu  ) );
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
	return p[2] + 2.*var/mu + log(1. + a*exp( -(2*var+L)/mu ) + a*a*exp(-2*L/mu) + a*a*a*exp( -(2*var+3*L)/mu ) ) - log(1. + a*exp( -(-2*var+L)/mu ) + a*a*exp(-2.*L/mu) + a*a*a*exp( -(-2.*var+3*L)/mu  ) );
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
	return p[2] + 2.*var/mu + log(1. + a*exp( -(2*var+L)/mu ) + a*a*exp(-2*L/mu) + a*a*a*exp( -(2*var+3*L)/mu ) ) - log(1. + a*exp( -(-2*var+L)/mu ) + a*a*exp(-2.*L/mu) + a*a*a*exp( -(-2.*var+3*L)/mu  ) );
	//double A = exp(L/mu)*(a-1) - exp(2.*var/mu)*a;
	//double B = exp( (L+2*var)/mu )*(a-1) - a;
	//return log( A / B );
}
