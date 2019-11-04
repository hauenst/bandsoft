#include <cmath>
#include <iostream>
#include <fstream>

#include "reader.h"
#include "bank.h"

#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace std;

int getKey( hipo::bank bandbank, int idx );
void readFile( TString inputFile, TH1D * hist );
double single_edge( double *x, double *p);
double gaus_edge( double *x, double *p);
TGraph * histDerivative( TH1D * hist , TH1D * back );
TCanvas * cobaltEdge( TH1D * cobalt);
TCanvas * sodiumEdge( TH1D * sodium);
void backgroundSubtract( TH1D * sig, TH1D * back);
void LoadLROffsets();
double TDC_TDIFF[600] = {0};
double FADC_TDIFF[600] = {0};

int main( int argc, char** argv){

	gStyle->SetOptFit(1);

	if( argc != 6 ){
		cerr << "Incorrect number of arguments. Please use instead:\n\t"
			<< "./singlebar [Co60-hipoFile] [Cs137-hipoFile] "
			<< "[Na22-hipoFile] [Background-hipoFile] [output-rootFile]\n";
		return -1;
	}
	LoadLROffsets();

	// Create strings for source files for this particular bar
	TString inputCo = argv[1];
	TString inputCs = argv[2];
	TString inputNa = argv[3];
	TString inputBa = argv[4];	

	// Create output file
	TFile * outFile = new TFile(argv[5],"RECREATE");

	// Prep histograms for the background spectrum, and source spectra
	TH1D * baHist = new TH1D("baHist","baHist",750,0,30000);
	//TH1D * coHist = new TH1D("coHist","coHist",750,0,30000);
	TH1D * coHist = new TH1D("coHist","coHist",5000,-5,5);
	TH1D * csHist = new TH1D("csHist","csHist",750,0,30000);
	TH1D * naHist = new TH1D("naHist","naHist",750,0,30000);
	

	// Read in files and store relevant histograms
	readFile( inputBa, baHist );
	readFile( inputCo, coHist );
	//readFile( inputCs, csHist );
	//readFile( inputNa, naHist );

	// Subtract background from all source files:
	//backgroundSubtract( coHist, baHist );
	//backgroundSubtract( naHist, baHist );
	//backgroundSubtract( csHist, baHist );

	//// Try both methods for cobalt edge:
	//TCanvas * coFitResult = (TCanvas *) cobaltEdge( coHist );
	//TCanvas * naFitResult = (TCanvas *) sodiumEdge( naHist );
	//TGraph * coDeri = (TGraph*) histDerivative( coHist , baHist );

	outFile->cd();
	//coDeri->Write();
	//coFitResult->Write();
	//naFitResult->Write();
	baHist->Write();
	coHist->Write();
	csHist->Write();
	naHist->Write();
	outFile->Close();

	return 0;
}

int getKey( hipo::bank bandbank, int idx ){
	int sector    = bandbank.getInt(0,idx);
	int layer     = bandbank.getInt(1,idx);
	int component = bandbank.getInt(2,idx);
	return sector*100 + layer*10 + component;
}

void readFile( TString inputFile, TH1D * hist ){

	// Initialize reader
	hipo::reader reader;
	hipo::dictionary factory;
	hipo::event readevent; 

	reader.open(inputFile);
	reader.readDictionary(factory); 
	hipo::bank	BAND_ADC	(factory.getSchema("BAND::adc"));
	hipo::bank	BAND_TDC	(factory.getSchema("BAND::tdc"));
	while(reader.next()==true){
		reader.read(readevent);

		//Load explicitly all information for each bank for the event
		readevent.getStructure(BAND_ADC);
		readevent.getStructure(BAND_TDC);

		int nADC = BAND_ADC.getRows();
		int nTDC = BAND_TDC.getRows();
		// We should expect 2 ADC and 2 TDC for a good source hit:
		if( !(nADC == 2 && nTDC == 2) ) continue;
		
		// Grab ADC info
		int adc_barKey = 0;
		double adcL, tFadcL, adcR, tFadcR = 0;
		// Read ADC information
		for(int aIdx = 0 ; aIdx < nADC ; aIdx++){
			adc_barKey = getKey( BAND_ADC, aIdx );
			int   ADC_order     = BAND_ADC.getInt  (3,aIdx);
			if( ADC_order == 0 ){
				adcL 	= 	BAND_ADC.getInt(4,aIdx);
				tFadcL	= 	BAND_ADC.getFloat(5,aIdx);
			}
			else if( ADC_order == 1 ){
				adcR	= 	BAND_ADC.getInt(4,aIdx);
				tFadcR 	= BAND_ADC.getFloat(5,aIdx);
			}
		}
	
		// Grab TDC info
		int tdc_barKey = 0;
		double tTdcL, tTdcR = 0;
		for(int tIdx = 0 ; tIdx < nTDC ; tIdx++){
			int   TDC_order     = BAND_TDC.getInt  (3,tIdx);
			TDC_order -= 2;
			tdc_barKey = getKey( BAND_TDC, tIdx);
			if( TDC_order == 0 ) tTdcL = BAND_TDC.getInt(4,tIdx)*0.02345;
			else if( TDC_order == 1 ) tTdcR = BAND_TDC.getInt(4,tIdx)*0.02345;
		}
		
		// Check that everything is non-zero for a good event:
		if( adcL == 0 || adcR == 0 || tFadcL == 0 || tFadcR == 0 || tTdcL == 0 || tTdcR == 0 ) continue;
		
		//if( fabs(tFadcL - tFadcR - FADC_TDIFF[adc_barKey] ) > 4 ) continue;
		//if( fabs(tTdcL - tTdcR - TDC_TDIFF[tdc_barKey] ) > 4 ) continue;
		if( sqrt(adcL*adcR) < 1500 ) continue;
		hist->Fill( log(adcL/adcR) );
		//hist->Fill( sqrt(adcL*adcR) );

	}

	return;
}


double single_edge( double *x, double *p){
	double var = *x;
	double alpha = 0.5*( p[0]*var + p[1] );
	double beta = p[0];
	//return p[0]*exp( -pow(var-p[1],2)/(2*p[2]*p[2]) ) + 
	//		alpha * erfc( (var - p[5]) / (sqrt(2)*p[6] ) ) + beta * exp( -pow(var-p[5],2)/(2*p[6]*p[6]));
	return 	alpha * erfc( (var - p[2]) / (sqrt(2)*p[3] ) ) + beta * exp( -pow(var-p[2],2)/(2*p[3]*p[3]));
	
}

double gaus_edge( double *x, double *p){
	double var = *x;
	return 	p[0] * exp( -pow(var-p[1],2)/(2*p[2]*p[2]));
	
}

TGraph * histDerivative( TH1D * hist , TH1D * back ){
	// Do background subtraction, normalizing the cosmic peaks
	double background = back->Integral( back->FindBin(13000),back->FindBin(22000) );
	double signal = hist->Integral( hist->FindBin(13000),hist->FindBin(22000) );
	back->Scale( signal / background );
	hist->Add( back, -1 );

	std::vector<double> xs;
	std::vector<double> derivs;
	int end = hist->FindLastBinAbove(10);
	double h = hist->GetXaxis()->GetBinWidth(1);
	for( int bin = 1 ; bin < end ; bin ++ ){
		// Get values of bins 2 away:
		double f2hm = hist->GetBinContent( bin - 2 );
		double f1hm = hist->GetBinContent( bin - 1 );
		double f1hp = hist->GetBinContent( bin + 1 );
		double f2hp = hist->GetBinContent( bin + 2 );
		if( f2hm == 0 || f1hm == 0 || f1hp == 0 || f2hp == 0 ) continue;
		
		double est = ( f2hm - 8.*f1hm + 8.*f1hp - f2hp ) / (12.*h );
		derivs.push_back( est );
		xs.push_back( hist->GetXaxis()->GetBinCenter(bin) );
	}
	int dim = xs.size();
	// Find local min at the end of the vector:
	double edgeAdc = 0;
	double errEst = 0;
	for( int i = dim-1; i > 0 ; i--){
		double curr = derivs.at(i);
		double next = derivs.at(i-1);
		if( next > (curr+0.75) && edgeAdc == 0.){
			edgeAdc = xs.at(i);
		}
		else if( next > (curr+2) ){
			errEst = fabs(edgeAdc - xs.at(i-1));
			break;
		}
	}
	cout << edgeAdc << " " << errEst << "\n";
	
	return new TGraph( dim, &xs[0], &derivs[0] );
        
}


TCanvas * cobaltEdge( TH1D * cobalt){
	// Try to determine where to start trying to fit for the edge by going away from peak
	double peakMean = cobalt->GetXaxis()->GetBinCenter(cobalt->GetMaximumBin());
	TF1 * gausFit = new TF1("gausFit","gaus",0,5000);
	gausFit->SetParameter(0,cobalt->GetMaximum()) ; gausFit->SetParameter(1,peakMean);
	TCanvas * junk = new TCanvas("junk");
	cobalt->Fit(gausFit,"QESR");

	// Fit the edge
	double startEdge = gausFit->GetParameter(1) + gausFit->GetParameter(2)*3.5;
	TF1 * trialFit = new TF1("trialFit",single_edge,startEdge,10000,4);
	trialFit->SetParameter(0,-1);
	trialFit->SetParameter(1,5000);
	trialFit->SetParameter(2,2300);
	trialFit->SetParameter(3,200);

	TCanvas * c1_co = new TCanvas("c1_co");
	cobalt->Draw("h");
	cobalt->Fit(trialFit,"QESR");
	trialFit->Draw("same");
	c1_co->Update();
	
	cout << trialFit->GetParameter(2) << " " << trialFit->GetParameter(3) << "\n";
	delete junk;
	return c1_co;
}

TCanvas * sodiumEdge( TH1D * sodium){
	// Try to determine where to start trying to fit for the edge by going away from peak
	double peakMean = sodium->GetXaxis()->GetBinCenter(sodium->GetMaximumBin());
	TF1 * initGaus = new TF1("initGaus","gaus",0,5000);
	initGaus->SetParameter(0,sodium->GetMaximum()) ; initGaus->SetParameter(1,peakMean);
	TCanvas * junk = new TCanvas("junk");
	sodium->Fit(initGaus,"QESR");

	// Fit the left edge
	double startEdge = initGaus->GetParameter(1) - initGaus->GetParameter(2)*5;
	double endEdge = initGaus->GetParameter(1) + initGaus->GetParameter(2)*0.5;
	TF1 * newFit = new TF1("newFit",gaus_edge,startEdge,endEdge,3);
	newFit->SetParameter(0,initGaus->GetParameter(0));
	newFit->SetParameter(1,initGaus->GetParameter(1));
	newFit->SetParameter(2,initGaus->GetParameter(2)*0.7);

	// Now create a new histogram and subtract off gaussian peak to see what is left over:
	TString tit = Form("%s_subtracted",sodium->GetTitle());
	int bins = sodium->GetXaxis()->GetNbins();
	TH1D * newHist = new TH1D(tit,tit,750,0,30000);
	for( int i = 1; i <= bins ; i++ ){
		double x = newHist->GetXaxis()->GetBinCenter(i);
		newHist->SetBinContent(i, sodium->GetBinContent(i) - newFit->Eval(x) );
		cout << x << " " << sodium->GetBinContent(i) << " " << sodium->GetBinContent(i) - newFit->Eval(x) << "\n";
	}
	newHist->SetLineColor(2);
	newHist->SetLineWidth(3);

	TCanvas * c1_na = new TCanvas("c1_na");
	sodium->Draw("h");
	sodium->Fit(newFit,"QESR");
	newFit->Draw("same");
	newHist->Draw("same");
	c1_na->Update(); c1_na->Modified();
	
	delete junk;
	//delete newHist;
	return c1_na;
}

void backgroundSubtract( TH1D * sig, TH1D * back){
	double background = back->Integral( back->FindBin(13000),back->FindBin(22000) );
	double signal = sig->Integral( sig->FindBin(13000),sig->FindBin(22000) );
	back->Scale( signal / background );
	sig->Add( back, -1 );
	back->Scale( background / signal );
	return;
}


void LoadLROffsets(){
	ifstream f;
	int sector, layer, component, barId;
	double tdc_off, fadc_off, temp;

	f.open("../include/lr_offsets.txt");
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
