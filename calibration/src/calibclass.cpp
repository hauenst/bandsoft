#include "calibclass.h"


// Constructor -----------------------------------------------
calibclass::calibclass(){
}
// Denstructor -----------------------------------------------
calibclass::~calibclass(){}

// Getting functions -----------------------------------------------
// 	Bar TDC:
std::map<int,double> calibclass::GetBars_TDC_tDiff	(){return Bar_TDC_tdiff;}
std::map<int,double> calibclass::GetBars_TDC_tMean	(){return Bar_TDC_tmean;}
std::map<int,double> calibclass::GetBars_TDC_xPos	(){return Bar_TDC_xpos;}
// 	Bar FADC:
std::map<int,double> calibclass::GetBars_FADC_tDiff	(){return Bar_FADC_tdiff;}
std::map<int,double> calibclass::GetBars_FADC_tMean	(){return Bar_FADC_tmean;}
std::map<int,double> calibclass::GetBars_FADC_xPos	(){return Bar_FADC_xpos;}
std::map<int,double> calibclass::GetBars_FADC_gmAdc	(){return Bar_FADC_gmadc;}
std::map<int,double> calibclass::GetBars_FADC_gmAmp	(){return Bar_FADC_gmamp;}
std::map<int,double> calibclass::GetBars_FADC_lnAdc	(){return Bar_FADC_lnadc;}
std::map<int,double> calibclass::GetBars_FADC_lnAmp	(){return Bar_FADC_lnamp;}
// 	PMT FADC:
std::map<int,double> calibclass::GetPMTs_FADC_adc	(){return PMT_FADC_adc;}
std::map<int,double> calibclass::GetPMTs_FADC_amp	(){return PMT_FADC_amp;}
std::map<int,double> calibclass::GetPMTs_FADC_time	(){return PMT_FADC_time;}
// 	PMT TDC:
std::map<int,double> calibclass::GetPMTs_TDC_time	(){return PMT_TDC_time;}

// Filling for getting functions PMTs -----------------------------------------------
void calibclass::InitGettersPMTs(){
	PMT_FADC_adc.clear();
	PMT_FADC_amp.clear();
	PMT_FADC_time.clear();
	PMT_TDC_time.clear();

	// Do PMT filling:
	for( pmt_iter = PMTs.begin() ; pmt_iter != PMTs.end() ; pmt_iter++ ){
		int pmtKey = pmt_iter->first;
		PMTHit thisPMT = pmt_iter->second;

		PMT_FADC_adc[pmtKey] 	= thisPMT.adc;
		PMT_FADC_amp[pmtKey] 	= thisPMT.amp;
		PMT_FADC_time[pmtKey]	= thisPMT.ftdc;
		PMT_TDC_time[pmtKey]	= thisPMT.tdc;
	}
	// After filling the memory, let's clear up the PMT struct:
	PMTs.clear();
}

// Filling for getting functions Bars -----------------------------------------------
void calibclass::InitGettersBars(){
	Bar_TDC_tdiff.clear();
	Bar_TDC_tmean.clear();
	Bar_TDC_xpos.clear();
	Bar_FADC_tdiff.clear();
	Bar_FADC_tmean.clear();
	Bar_FADC_xpos.clear();
	Bar_FADC_gmadc.clear();
	Bar_FADC_gmamp.clear();
	
	// Do Bar filling:
	for( bar_iter = Bars.begin() ; bar_iter != Bars.end() ; bar_iter++ ){
		int barKey = bar_iter->first;
		BarHit thisBar = bar_iter->second;
	
		Bar_TDC_tdiff[barKey] = thisBar.tdc_tdiff;
		Bar_TDC_tmean[barKey] = thisBar.tdc_tmean;
		Bar_TDC_xpos[barKey]  = thisBar.tdc_xpos;
		Bar_FADC_tdiff[barKey] = thisBar.ftdc_tdiff;
		Bar_FADC_tmean[barKey] = thisBar.ftdc_tmean;
		Bar_FADC_xpos[barKey]  = thisBar.ftdc_xpos;
		Bar_FADC_gmadc[barKey] = thisBar.adc;
		Bar_FADC_gmamp[barKey] = thisBar.amp;
		Bar_FADC_lnadc[barKey] = thisBar.lnR_adc;
		Bar_FADC_lnamp[barKey] = thisBar.lnR_amp;
	}
	// After filling the memory, let's clear up the Bar struct:
	Bars.clear();
}

// Parse hipo banks for PMT structure -----------------------------------------------
void calibclass::CreatePMTs(hipo::bank BAND_ADC, hipo::bank BAND_TDC , double phaseCorr ){
	PMTs.clear();

	int nADC = BAND_ADC.getRows();
	int nTDC = BAND_TDC.getRows();

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
				TDC_tdc -  ( parampB_L[barKey] + parampA_L[barKey]/pow(PMTs[pmtKey].amp,parampC_L[barKey]) )
			:	TDC_tdc -  ( parampB_R[barKey] + parampA_R[barKey]/pow(PMTs[pmtKey].amp,parampC_R[barKey]) ) ;
		//		TDC_tdc -  ( paradcB_L[barKey] + paradcA_L[barKey]/sqrt(PMTs[pmtKey].adc) )
		//	:	TDC_tdc -  ( paradcB_R[barKey] + paradcA_R[barKey]/sqrt(PMTs[pmtKey].adc) ) ;


		double prevdiff = fabs( PMTs[pmtKey].ftdc - PMTs[pmtKey].tdc );
		double currdiff = fabs( PMTs[pmtKey].ftdc - TDC_tdc );
		if( currdiff < prevdiff ) PMTs[pmtKey].tdc = TDC_tdc; // if we have a more similar time, save that TDC instead
	}
}

// Parse fired PMTs for valid Bars  -----------------------------------------------
void calibclass::CreateBars(){
	Bars.clear();

	// Loop through all pmts and do matching
	for( pmt_iter = PMTs.begin() ; pmt_iter != PMTs.end() ; pmt_iter++ ){
		if( pmt_iter->second.processed == true ) continue;
		pmt_iter->second.processed = true;

		// Parse key into ID:
		int pmtKey = pmt_iter->first;
		int order = (pmtKey % 10);
		int component = (pmtKey - order)/10 % 10;
		int layer = (pmtKey - order - component*10)/100 % 10;
		int sector  = (pmtKey - order - component*10 - layer*100)/1000 % 10;
		int barKey = sector * 100 + layer * 10 + component;

		// Get this PMT and the matching PMT:
		PMTHit thisPMT = pmt_iter->second;
		int pmtMate = (order == 0) ? (pmtKey+1) : (pmtKey-1);
		if( PMTs.count(pmtMate) == 0 ) continue; // cannot match this bar as it doesn't have partner

		PMTHit matePMT = PMTs[pmtMate];
		matePMT.processed = true; // also flag the matching PMT so we don't double-process

		// if order == 0, then thisPMT is L and if order == 1 then thisPMT is R
		int sign = (order == 0) ? 1 : -1;
		// Skip this match if the TDC time = 0 or ADC time = 0 
		// 		   or the ADC = 0 or AMP = 0	
		// 		   or the AMP is overflow
		// 		   or the TDC time was never set or the FADC time was never set
		// 		   or if the ADC is too low
		if( 	thisPMT.tdc == 0 || matePMT.tdc == 0 || thisPMT.ftdc == 0 || matePMT.ftdc == 0 || 
			thisPMT.adc == 0 || matePMT.adc == 0 ||thisPMT.amp == 0 || matePMT.amp == 0 || 
			thisPMT.amp >= 4095 || matePMT.amp >= 4095 || 
			thisPMT.tdc == 1e10 || matePMT.tdc == 1e10 || thisPMT.ftdc == 1e10 || matePMT.ftdc == 1e10 || 
			thisPMT.adc < 4000 || matePMT.adc < 4000 ) continue;

		// Apply correction to tdiff:
		thisPMT.tdc = thisPMT.tdc 	- (sign)*TDC_TDIFF[barKey];
		thisPMT.ftdc = thisPMT.ftdc 	- (sign)*FADC_TDIFF[barKey];

		// Create the variables we are interested in saving for a bar
		double tdcdiff 	= (sign)*(	thisPMT.tdc - matePMT.tdc	);
		double ftdcdiff = (sign)*(	thisPMT.ftdc - matePMT.ftdc	);
		double tdcxpos 	= (sign)*(	thisPMT.tdc - matePMT.tdc	)*(-TDC_VEFF[barKey]/2.);
		double ftdcxpos = (sign)*(	thisPMT.ftdc - matePMT.ftdc	)*(-FADC_VEFF[barKey]/2.);
		double lnR_adc 	= (sign)*(	log(thisPMT.adc / matePMT.adc )	);
		double lnR_amp 	= (sign)*(	log(thisPMT.amp / matePMT.amp )	);
		double tdcmean = (thisPMT.tdc + matePMT.tdc)/2.;
		double ftdcmean = (thisPMT.ftdc + matePMT.ftdc)/2.;
		double gm_adc = sqrt( thisPMT.adc * matePMT.adc );
		double gm_amp = sqrt( thisPMT.amp * matePMT.amp );

		// Store the bar
		BarHit newBar;
		newBar.tdc_tdiff = tdcdiff;
		newBar.ftdc_tdiff = ftdcdiff;
		newBar.tdc_xpos = tdcxpos;
		newBar.ftdc_xpos = ftdcxpos;
		newBar.tdc_tmean = tdcmean;
		newBar.ftdc_tmean = ftdcmean;
		newBar.adc = gm_adc;
		newBar.amp = gm_amp;
		newBar.lnR_adc = lnR_adc;
		newBar.lnR_amp = lnR_amp;
		Bars[barKey] = newBar;
	} // end loop over pmts
}

// Load timewalk calibration -----------------------------------------------
void calibclass::LoadTimeWalk(){
	ifstream f;
	int sector, layer, component, barId;
	double parA, parB, parC, temp;
	f.open("../../include/time_walk_corr_adc_left.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> parA;
		f >> parB;
		f >> parC;
		f >> temp;
		f >> temp;
		f >> temp;
		paradcA_L[barId] = parA;
		paradcB_L[barId] = parB;
	}
	f.close();

	f.open("../../include/time_walk_corr_adc_right.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> parA;
		f >> parB;
		f >> parC;
		f >> temp;
		f >> temp;
		f >> temp;
		paradcA_R[barId] = parA;
		paradcB_R[barId] = parB;
	}
	f.close();

	f.open("../../include/time_walk_corr_amp_left.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> parA;
		f >> parB;
		f >> parC;
		f >> temp;
		f >> temp;
		f >> temp;
		parampA_L[barId] = parA;
		parampB_L[barId] = parB;
		parampC_L[barId] = parC;
	}
	f.close();

	f.open("../../include/time_walk_corr_amp_right.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> parA;
		f >> parB;
		f >> parC;
		f >> temp;
		f >> temp;
		f >> temp;
		parampA_R[barId] = parA;
		parampB_R[barId] = parB;
		parampC_R[barId] = parC;
	}
	f.close();
	return;
}

// Load lr offset calibration -----------------------------------------------
void calibclass::LoadLROffsets(){
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

// Load effective velocity calibration -----------------------------------------------
void calibclass::LoadVelocityMap(){
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
