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
void drawHist( TH2F * hist , TCanvas * c , int cd );
void offsetFit( TH2F * hist , TCanvas * c , int cd , int is , int il, int ic, int flag ,double &offset, double&veff , double &width );
void fitAtten(TH2F * hist , TCanvas* c, int cd , int is, double mid, double &par1, double &par2 , double &par1_err, double &par2_err );
int doProj( TH2F * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE);
void attenCorr(	std::vector<double> *adcs		,
		std::vector<double> *adcsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		TH2F * hist				);
double lnR_reflect( double *x , double *p);

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

int main(int argc, char** argv) {
	if( argc!= 2 ){
		cerr << "Wrong number of arguments. Please use instead:\n\t./code [InputFile]\n";
		return -1;
	}
	gStyle->SetTitleSize(0.2,"t");
	gStyle->SetOptStat(0);

	TString inputFile = argv[1];

	// Load calibration parameters
	LoadTimeWalk();
	cout << "loaded 1\n";
	LoadLROffsets();
	cout << "loaded 1\n";
	LoadVelocityMap();
	cout << "loaded 1\n";

	const int nHistos = 600;
	TH2F ** h2_tdc_tdiff_adc = new TH2F * [nHistos];
	TH2F ** h2_ftdc_tdiff_adc = new TH2F * [nHistos];
	TH2F ** h2_tdc_tdiff_amp = new TH2F * [nHistos];
	TH2F ** h2_ftdc_tdiff_amp = new TH2F * [nHistos];
	for(int i = 0 ; i < nHistos ; i++){
		h2_tdc_tdiff_adc[i] = new TH2F(Form("h2_tdc_tdiff_adc_%i",i),";t_{TDC,L} - t_{TDC,R} [ns];ln(ADC_{L}/ADC_{R})",500,-25,25,400,-10,10);
		h2_ftdc_tdiff_adc[i] = new TH2F(Form("h2_ftdc_tdiff_adc_%i",i),";t_{FTDC,L} - t_{FTDC,R} [ns];ln(ADC_{L}/ADC_{R});",500,-25,25,400,-10,10);
		h2_tdc_tdiff_amp[i] = new TH2F(Form("h2_tdc_tdiff_amp_%i",i),";t_{TDC,L} - t_{TDC,R} [ns];ln(AMP_{L}/AMP_{R})",500,-25,25,400,-10,10);
		h2_ftdc_tdiff_amp[i] = new TH2F(Form("h2_ftdc_tdiff_amp_%i",i),";t_{FTDC,L} - t_{FTDC,R} [ns];ln(AMP_{L}/AMP_{R});",500,-25,25,400,-10,10);
		PrettyTH2F(h2_tdc_tdiff_adc[i]);
		PrettyTH2F(h2_ftdc_tdiff_adc[i]);
		PrettyTH2F(h2_tdc_tdiff_amp[i]);
		PrettyTH2F(h2_ftdc_tdiff_amp[i]);
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
		event_counter++;

                reader.read(readevent); // new hipo4
                readevent.getStructure(BAND_ADC );   // new hipo4
		readevent.getStructure(BAND_TDC );   // new hipo4
		readevent.getStructure(RUN_config );   // new hipo4

		long timestamp = RUN_config.getLong(4,0);
		double phaseCorr = getTriggerPhase(timestamp);

		int nADC = BAND_ADC.getRows();
		int nTDC = BAND_TDC.getRows();

		// Skip events with no entries
		if(nADC==0||nTDC==0) continue;
		std::map<int,PMTHit> PMTs;

		// First we need to get an ADC and a TDC for all PMTs:
		for(int aIdx = 0 ; aIdx < nADC ; aIdx++){
			int   ADC_sector    	= BAND_ADC.getInt  (0,aIdx);
                        int   ADC_layer     	= BAND_ADC.getInt  (1,aIdx);
                        int   ADC_component 	= BAND_ADC.getInt  (2,aIdx);
                        int   ADC_order     	= BAND_ADC.getInt  (3,aIdx);
			double ADC_adc = (double)(BAND_ADC.getInt  (4,aIdx));
			double ADC_amp = (double)(BAND_ADC.getInt  (5,aIdx));
			double ADC_tdc = (double)(BAND_ADC.getFloat  (6,aIdx));
			
			int pmtKey = 1000*ADC_sector + 100*ADC_layer + 10*ADC_component + ADC_order;
			if( PMTs.count(pmtKey) ){
				if( PMTs[pmtKey].ftdc > ADC_tdc ){
					// replace this PMT
					PMTs[pmtKey].adc = ADC_adc;
					PMTs[pmtKey].ftdc = ADC_tdc;
					PMTs[pmtKey].amp = ADC_amp;
				}
			}
			else{
				PMTHit newHit;
				newHit.adc = ADC_adc;
				newHit.ftdc = ADC_tdc;
				newHit.amp = ADC_amp;
				PMTs[pmtKey] = newHit;
			}
		}
		for(int tIdx = 0 ; tIdx < nTDC ; tIdx++){
			int   TDC_sector    = BAND_TDC.getInt  (0,tIdx);
                        int   TDC_layer     = BAND_TDC.getInt  (1,tIdx);
                        int   TDC_component = BAND_TDC.getInt  (2,tIdx);
                        int   TDC_order     = BAND_TDC.getInt  (3,tIdx);
			TDC_order -= 2;
			double TDC_tdc	    = (BAND_TDC.getInt(4,tIdx)) * 0.02345 - phaseCorr;

			int pmtKey = 1000*TDC_sector + 100*TDC_layer + 10*TDC_component + TDC_order;
			
			if( PMTs.count(pmtKey) == 0 ) continue; // if I don't have an ADC for this PMT already, skip it.
		
			// Do timewalk correction:
			int barKey = 100*TDC_sector + 10*TDC_layer + TDC_component;
			TDC_tdc = (TDC_order == 0) ? 
					TDC_tdc -  ( parB_L[barKey] + parA_L[barKey]/sqrt(PMTs[pmtKey].adc) )
				:	TDC_tdc -  ( parB_R[barKey] + parA_R[barKey]/sqrt(PMTs[pmtKey].adc) ) ;

				

			double prevdiff = fabs( PMTs[pmtKey].ftdc - PMTs[pmtKey].tdc );
			double currdiff = fabs( PMTs[pmtKey].ftdc - TDC_tdc );
			if( currdiff < prevdiff ) PMTs[pmtKey].tdc = TDC_tdc; // if we have a more similar time, save that TDC instead
		}
	
		// Now that we have single PMT hits, we just to to construct our bars:
		std::map<int,PMTHit>::iterator iter;
		for( iter = PMTs.begin() ; iter != PMTs.end() ; iter++ ){
			if( iter->second.processed == true ) continue;
			iter->second.processed = true;

			// Parse key into ID:
			int pmtKey = iter->first;
			int order = (pmtKey % 10);
			int component = (pmtKey - order)/10 % 10;
			int layer = (pmtKey - order - component*10)/100 % 10;
			int sector  = (pmtKey - order - component*10 - layer*100)/1000 % 10;
			int barKey = sector * 100 + layer * 10 + component;

			// Get this PMT and the matching PMT:
			PMTHit thisPMT = iter->second;
			int pmtMate = (order == 0) ? (pmtKey+1) : (pmtKey-1);
			if( PMTs.count(pmtMate) == 0 ) continue; // cannot match this bar as it doesn't have partner

			PMTHit matePMT = PMTs[pmtMate];
			matePMT.processed = true; // also flag the matching PMT so we don't double-process

			// Plot bar ln(adcL/adcR) : tL - tR:
				// if order == 0, then thisPMT is L and if order == 1 then thisPMT is R
			int sign = (order == 0) ? 1 : -1;
			if( thisPMT.tdc == 0 || matePMT.tdc == 0 || thisPMT.ftdc == 0 || matePMT.ftdc == 0 || thisPMT.adc == 0 || matePMT.adc == 0 ||
				thisPMT.amp == 0 || matePMT.amp == 0 || thisPMT.amp >= 4095 || matePMT.amp >= 4095 ) continue;

			// Apply correction to tdiff:
			thisPMT.tdc = thisPMT.tdc - sign*TDC_TDIFF[barKey];
			thisPMT.ftdc = thisPMT.ftdc - sign*FADC_TDIFF[barKey];

			double tdcdiff 	= (sign)*(	thisPMT.tdc - matePMT.tdc	);
			double ftdcdiff = (sign)*(	thisPMT.ftdc - matePMT.ftdc	);
			double lnR_adc 	= (sign)*(	log(thisPMT.adc / matePMT.adc )	);
			double lnR_amp 	= (sign)*(	log(thisPMT.amp / matePMT.amp )	);

			h2_tdc_tdiff_adc[barKey] -> Fill( tdcdiff , lnR_adc );
			h2_ftdc_tdiff_adc[barKey] -> Fill( ftdcdiff , lnR_adc );
			h2_tdc_tdiff_amp[barKey] -> Fill( tdcdiff , lnR_amp );
			h2_ftdc_tdiff_amp[barKey] -> Fill( ftdcdiff , lnR_amp );

		}
		
	
	} // end loop over events

	
	TCanvas * c0 = new TCanvas("c0","c0",700,900);
	TCanvas * c1 = new TCanvas("c1","c0",700,900);
	double tdc_off		[600] = {0};
	double ftdc_off		[600] = {0};
	double tdc_veff		[600] = {0};
	double ftdc_veff	[600] = {0};
	c0 -> Modified();
	c0 -> Update(); 
	c1 -> Modified();
	c1 -> Update(); 
	TCanvas *** cSLC_tdc = new TCanvas**[5];
	TCanvas *** cSLC_ftdc = new TCanvas**[5];
	c0 -> Print("corr2nd_tdc_results_offsets.pdf(");
	c1 -> Print("corr2nd_ftdc_results_offsets.pdf(");
	
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

				if( h2_tdc_tdiff_amp[identifier]->Integral()  ) drawHist( h2_tdc_tdiff_amp[identifier] , cSLC_tdc[is][il] ,  3*cIdx+1 );
				if( h2_tdc_tdiff_adc[identifier]->Integral()  ){
					drawHist( h2_tdc_tdiff_adc[identifier] , cSLC_tdc[is][il] ,  3*cIdx+2 );
					// do projection for l-r and veff:
					double offset, veff, width;
					offsetFit( h2_tdc_tdiff_adc[identifier] , cSLC_tdc[is][il] , 3*cIdx+3 , is, il, cIdx , 0 ,offset, veff, width );
					tdc_off[identifier] = offset;
					tdc_veff[identifier] = veff;
				}
				if( h2_ftdc_tdiff_amp[identifier]->Integral()  ) drawHist( h2_ftdc_tdiff_amp[identifier] , cSLC_ftdc[is][il] ,  3*cIdx+1 );
				if( h2_ftdc_tdiff_adc[identifier]->Integral()  ){
					drawHist( h2_ftdc_tdiff_adc[identifier] , cSLC_ftdc[is][il] ,  3*cIdx+2 );
					// do projection for l-r and veff:
					double offset, veff, width;
					offsetFit( h2_ftdc_tdiff_adc[identifier] , cSLC_ftdc[is][il] , 3*cIdx+3 , is, il, cIdx , 1 ,offset, veff, width );
					ftdc_off[identifier] = offset;
					ftdc_veff[identifier] = veff;
				}
			}
			cSLC_tdc[is][il] -> Modified();     
			cSLC_tdc[is][il] -> Update();
			cSLC_tdc[is][il] -> Print("corr2nd_tdc_results_offsets.pdf");
			cSLC_ftdc[is][il] -> Modified();     
			cSLC_ftdc[is][il] -> Update();
			cSLC_ftdc[is][il] -> Print("corr2nd_ftdc_results_offsets.pdf");
		}
	}
	c0 -> Print("corr2nd_tdc_results_offsets.pdf)");
	c1 -> Print("corr2nd_ftdc_results_offsets.pdf)");

	ofstream tabL, tabR;
	tabL.open("corr2nd_effective_velocity.txt");
	tabR.open("corr2nd_lr_offsets.txt");
	for(int is = 1 ; is <= 5 ; is++){
		for(int il = 1 ; il <= 6 ; il++){
			for(int ic = 1 ; ic <= slc[il-1][is-1] ; ic++){
				int idx = 100*is + 10*il + ic;
				tabL << is << "\t" << il << "\t" << ic << "\t"  << tdc_veff[idx] << "\t" << ftdc_veff[idx] << "\n";
				tabR << is << "\t" << il << "\t" << ic << "\t"  << tdc_off[idx] << "\t" << ftdc_off[idx] << "\n";
				
			}
		}
	}
	tabL.close();
	tabR.close();
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

	f.open("../../include/improved_lr_offsets.txt");
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
	int midbin = hist->GetXaxis()->FindBin(0.);
	double val = hist->ProjectionY("trash",midbin,midbin+1)->GetMean();
	double min_val = val - 1;
	double max_val = val + 1;
	hist->GetYaxis() -> SetRangeUser(min_val,max_val);
	hist->Draw("COLZ");
	return;
}

void offsetFit( TH2F * hist , TCanvas * c , int cd , int is , int il, int ic, int flag ,double &offset, double&veff , double &width ){
	// Do projection on x-axis for L-R offset fit:
	TH1D* projX = hist->ProjectionX(Form("%i_tdiff_%i_%i_%i",flag,(is+1),(il+1),(ic+1)));
	gPad -> SetLeftMargin(0.16);
	gPad -> SetBottomMargin(0.20);
	gPad -> SetTopMargin(0);
	c->cd(cd);
	projX->SetTitle(Form("Sector %i, Layer %i, Comp %i",(is+1),(il+1),(ic+1)));
	projX->Draw();
	double percent = 0.5;
	double len = BARLENGTHS[is];
	double thres = projX->Integral()*0.1 / (len/15.*2.) * percent;
	double xlower = projX->GetXaxis()->GetBinCenter( projX->FindFirstBinAbove( thres ) );
	double xupper = projX->GetXaxis()->GetBinCenter( projX->FindLastBinAbove( thres ) );
	TLine *lineL = new TLine(xlower,0,xlower,thres/percent);
	TLine *lineU = new TLine(xupper,0,xupper,thres/percent);
	TLine *lineM = new TLine((xlower+xupper)/2.,0,(xlower+xupper)/2.,thres/percent);
	lineL->SetLineColor(2); lineU->SetLineColor(2); lineM->SetLineColor(2);
	lineL->SetLineWidth(2); lineU->SetLineWidth(2); lineM->SetLineWidth(2);
	lineL->Draw("same"); lineU->Draw("same"); lineM->Draw("same");
	c->Update(); c->Modified();

	offset = (xlower+xupper)/2.;
	width = xupper-xlower;
	veff = len/(width/2.);

	return;
}

void attenCorr(	std::vector<double> *adcs		,
		std::vector<double> *adcsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		TH2F * hist				){
	int currBin = 0;
	while( currBin < 500 ){
		double xPt, yPt, yEr, ySig, ySigEr;
		int step = doProj( hist , currBin , false , 0, xPt, yPt, yEr, ySig, ySigEr );
		currBin += step ;
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
	char temp[100];
	sprintf(temp,"slice_%d_%d",flag,bin);
	TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
	
	// Getting the mean x value in this range:
	hist->GetXaxis()->SetRange(bin,bin+step);
	x = hist->GetMean(1);
	hist->GetXaxis()->SetRange();
	TFitResultPtr f = pj->Fit("gaus","QES","",-10,10);
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

void fitAtten(TH2F * hist , TCanvas* c, int cd , int is, double mid, double &par1, double &par2 , double &par1_err, double &par2_err ){
	std::vector<double> xs;
	std::vector<double> ys;
	std::vector<double> xErrs;
	std::vector<double> yErrs;
	std::vector<double> res;
	std::vector<double> resErrs;
	attenCorr( &xs, &xErrs, &ys, &yErrs, &res, &resErrs, hist );

	int dim = xs.size();
	TGraphErrors *g = new TGraphErrors(dim, &xs[0], &ys[0], &xErrs[0], &yErrs[0]);
	double len = BARLENGTHS[is];
	//TString formula = Form("log( (exp(%f/[0])*([1]-1.) - [1]*exp(2.*x/[0]) )  / ( exp((%f+2.*x)/[0])*([1]-1.) -[1] ) )",len,len);
	//TF1 * model = new TF1("atten_reflect",formula.Data(),-25,25,2);
	TF1 * model = new TF1("atten_reflect",lnR_reflect,-25,25,3);
	//TF1 * model = new TF1("timeWalk",wlk,0,20000,2);
	double zoom = 6;
	if( (is+1)==3 || (is+1)==4 ) zoom = 1.5;
	//TFitResultPtr ptr = g->Fit("pol1","QESR","",mid-zoom,mid+zoom);
	TFitResultPtr ptr = g->Fit(model,"QESR","",mid-zoom,mid+zoom);
	
	//par1 = ptr->Parameter(0);
	//par1_err = ptr->ParError(0);
	//par2 = ptr->Parameter(1);
	//par2_err = ptr->ParError(1);
	c->cd(cd);
	g->SetMarkerStyle(20);
	g->Draw("AP");	
	return;
}


double lnR_reflect( double *x , double *p){
	double var = *x;
	double mu = p[0];
	double a = p[1];
	// p[0] = mu
	// p[1] = alpha
	// p[2] = L
	double L = BARLENGTHS[1];
	double A = exp(L/mu)*(a-1) - exp(2.*var/mu)*a;
	double B = exp( (L+2*var)/mu )*(a-1) - a;
	return p[2] + log( A / B );

}

