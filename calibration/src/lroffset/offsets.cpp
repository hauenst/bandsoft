#include <cstdlib>
#include <iostream>
// Root includes:
#include "TFile.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLine.h"
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

void offsetFit( TH2F * hist , TCanvas * c , int cd , int is , int il, int ic, int flag ,double &offset, double&veff , double &width );

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
	cout << "\n\tLoaded BAND calibrations\n\n";

	// Define histograms to be filled:
	const int nHistos = 600;
	TH2F ** h2_tdc_tdiff_adc = new TH2F * [nHistos];
	TH2F ** h2_ftdc_tdiff_adc = new TH2F * [nHistos];
	TH2F ** h2_tdc_tdiff_amp = new TH2F * [nHistos];
	TH2F ** h2_ftdc_tdiff_amp = new TH2F * [nHistos];
	for(int i = 0 ; i < nHistos ; i++){
		h2_tdc_tdiff_adc[i] = new TH2F(Form("h2_tdc_tdiff_adc_%i",i),"	;t_{TDC,L} - t_{TDC,R} [ns];ln(ADC_{L}/ADC_{R})",500,-25,25,400,-10,10);
		h2_tdc_tdiff_amp[i] = new TH2F(Form("h2_tdc_tdiff_amp_%i",i),";t_{TDC,L} - t_{TDC,R} [ns];ln(AMP_{L}/AMP_{R})",500,-25,25,400,-10,10);

		h2_ftdc_tdiff_adc[i] = new TH2F(Form("h2_ftdc_tdiff_adc_%i",i),";t_{FTDC,L} - t_{FTDC,R} [ns];ln(ADC_{L}/ADC_{R});",500,-25,25,400,-10,10);
		h2_ftdc_tdiff_amp[i] = new TH2F(Form("h2_ftdc_tdiff_amp_%i",i),";t_{FTDC,L} - t_{FTDC,R} [ns];ln(AMP_{L}/AMP_{R});",500,-25,25,400,-10,10);
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

		// Create PMTs
		BAND->CreatePMTs( BAND_ADC, BAND_TDC, phaseCorr );
		// Create Bars from PMTs
		BAND->CreateBars();
		// Prepare maps to access information from bars
		BAND->InitGettersBars();
		
		// Access information we need from bars
		std::map<int,double> lnR_adc = BAND->GetBars_FADC_lnAdc();
		std::map<int,double> lnR_amp = BAND->GetBars_FADC_lnAmp();
		std::map<int,double> tdc_tdiff = BAND->GetBars_TDC_tDiff();
		std::map<int,double> ftdc_tdiff = BAND->GetBars_FADC_tDiff();

		// Loop through maps and fill our histograms:
		std::map<int,double>::iterator iter_bar;
		for( iter_bar = tdc_tdiff.begin() ; iter_bar != tdc_tdiff.end() ; iter_bar++ ){
			int barKey = iter_bar->first;
			int component = (barKey % 10);
			int layer = (barKey - component)/10 % 10;
			int sector = (barKey - component - layer*10)/100 % 10;

			h2_tdc_tdiff_adc[barKey] -> Fill( 	tdc_tdiff[barKey],	lnR_adc[barKey]	);
			h2_ftdc_tdiff_adc[barKey] -> Fill( 	ftdc_tdiff[barKey],	lnR_adc[barKey] );
			h2_tdc_tdiff_amp[barKey] -> Fill( 	tdc_tdiff[barKey],	lnR_amp[barKey]	);
			h2_ftdc_tdiff_amp[barKey] -> Fill( 	ftdc_tdiff[barKey],	lnR_amp[barKey] );
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
	double tdc_off		[600] = {0};
	double ftdc_off		[600] = {0};
	double tdc_veff		[600] = {0};
	double ftdc_veff	[600] = {0};
	
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

				if( h2_tdc_tdiff_amp[identifier]->Integral()  ){
					zoomTH2F( h2_tdc_tdiff_amp[identifier] , 0 , 1, 1 );
					drawTH2F( h2_tdc_tdiff_amp[identifier] , cSLC_tdc[is][il] ,  3*cIdx+1 );
				}
				if( h2_tdc_tdiff_adc[identifier]->Integral()  ){
					zoomTH2F( h2_tdc_tdiff_adc[identifier] , 0 , 1, 1 );
					drawTH2F( h2_tdc_tdiff_adc[identifier] , cSLC_tdc[is][il] ,  3*cIdx+2 );
					// do projection for l-r and veff:
					double offset = 0, veff = 0, width = 0;
					offsetFit( h2_tdc_tdiff_adc[identifier] , cSLC_tdc[is][il] , 3*cIdx+3 , is, il, cIdx , 0 ,offset, veff, width );
					tdc_off[identifier] = offset;
					tdc_veff[identifier] = veff;
				}
				if( h2_ftdc_tdiff_amp[identifier]->Integral()  ){
					zoomTH2F( h2_ftdc_tdiff_amp[identifier] , 0 , 1 , 1 );
					drawTH2F( h2_ftdc_tdiff_amp[identifier] , cSLC_ftdc[is][il] ,  3*cIdx+1 );
				}
				if( h2_ftdc_tdiff_adc[identifier]->Integral()  ){
					zoomTH2F( h2_ftdc_tdiff_adc[identifier] , 0 , 1 , 1 );
					drawTH2F( h2_ftdc_tdiff_adc[identifier] , cSLC_ftdc[is][il] ,  3*cIdx+2 );
					// do projection for l-r and veff:
					double offset = 0, veff = 0, width = 0;
					offsetFit( h2_ftdc_tdiff_adc[identifier] , cSLC_ftdc[is][il] , 3*cIdx+3 , is, il, cIdx , 1 ,offset, veff, width );
					ftdc_off[identifier] = offset;
					ftdc_veff[identifier] = veff;
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
	tabL.open("effective_velocity.txt");
	tabR.open("lr_offsets.txt");
	for(int is = 1 ; is <= 5 ; is++){
		for(int il = 1 ; il <= 6 ; il++){
			for(int ic = 1 ; ic <= slc[il-1][is-1] ; ic++){
				int idx = 100*is + 10*il + ic;
				tabL << is << "\t" << il << "\t" << ic << "\t"  << tdc_veff[idx] << "\t" << ftdc_veff[idx] << "\t" << 0 << "\t" << 0 << "\n";
				tabR << is << "\t" << il << "\t" << ic << "\t"  << tdc_off[idx] << "\t" << ftdc_off[idx] << "\t" << 0 << "\t" << 0 << "\n";
				
			}
		}
	}
	tabL.close();
	tabR.close();

	ofstream outfile;
	outfile.open("th2d_atten_lroffset.txt");
	for( int binx = 1; binx < h2_tdc_tdiff_adc[254]->GetXaxis()->GetNbins() ; binx++){
		for( int biny = 1; biny < h2_tdc_tdiff_adc[254]->GetYaxis()->GetNbins() ; biny++){
			outfile << h2_tdc_tdiff_adc[254]->GetXaxis()->GetBinCenter(binx) << " " 
				<< h2_tdc_tdiff_adc[254]->GetYaxis()->GetBinCenter(biny) << " "
				<< h2_tdc_tdiff_adc[254]->GetBinContent(binx,biny) 	 << " " 
				<< h2_tdc_tdiff_amp[254]->GetBinContent(binx,biny)	 << "\n";
		}
	}
	outfile.close();

	outfile.open("th1d_lroffset.txt");
	TH1D * proj = h2_tdc_tdiff_adc[254]->ProjectionX("proj");	
	for( int binx = 1; binx < proj->GetXaxis()->GetNbins() ; binx++){
			outfile << proj->GetXaxis()->GetBinCenter(binx) << " " 
				<< proj->GetBinContent(binx)	 << "\n";
	}
	outfile.close();
	
	

	
	return 0;
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
