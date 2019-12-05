#ifndef __CALIBCLASS_H__
#define __CALIBCLASS_H__

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "reader.h"
#include "bank.h"

using namespace std;

class calibclass
{
	// Structure for PMT hits
	struct PMTHit{
		int ID;
		double adc, amp;
		double tdc=1e10, ftdc=1e10;
		bool processed = false;
	};
	// Structure for Bar hits
	struct BarHit{
		int ID;
		double adc, amp;
		double tdc_tdiff, ftdc_tdiff;
		double tdc_tmean, ftdc_tmean;
		double tdc_xpos, ftdc_xpos;
		double lnR_adc, lnR_amp;
	};


	public:
		calibclass();
		~calibclass();
		
		void LoadTimeWalk();
		void LoadLROffsets();
		void LoadVelocityMap();

		void CreatePMTs(hipo::bank BAND_ADC, hipo::bank BAND_TDC, double phaseCorr );
		void CreateBars();
		void InitGettersPMTs();
		void InitGettersBars();

		// Getter functions:
		// 	Bar TDC:
		std::map<int,double> GetBars_TDC_tDiff	();
		std::map<int,double> GetBars_TDC_tMean	();
		std::map<int,double> GetBars_TDC_xPos	();
		// 	Bar FADC:
		std::map<int,double> GetBars_FADC_tDiff	();
		std::map<int,double> GetBars_FADC_tMean	();
		std::map<int,double> GetBars_FADC_xPos	();
		std::map<int,double> GetBars_FADC_gmAdc	();
		std::map<int,double> GetBars_FADC_gmAmp	();
		std::map<int,double> GetBars_FADC_lnAdc	();
		std::map<int,double> GetBars_FADC_lnAmp	();
		// 	PMT TDC:
		std::map<int,double> GetPMTs_TDC_time	();
		// 	PMT FADC:
		std::map<int,double> GetPMTs_FADC_adc	();
		std::map<int,double> GetPMTs_FADC_amp	();
		std::map<int,double> GetPMTs_FADC_time	();

	private:

		// Holders for timewalk:
		double parA_L[600] = {0.};		
		double parB_L[600] = {0.};		
		double parA_R[600] = {0.};		
		double parB_R[600] = {0.};		
		// Holders for LR offsets:
		double TDC_TDIFF[600] = {0.};	
		double FADC_TDIFF[600] = {0.};	
		// Holders for effective vel:
		double TDC_VEFF[600] = {0.};	
		double FADC_VEFF[600] = {0.};	

		// Storage container for PMT hits and iterator:
		std::map<int,PMTHit> PMTs;
		std::map<int,PMTHit>::iterator pmt_iter;
		
		// Storage container for Bar hits and iterator:
		std::map<int,BarHit> Bars;
		std::map<int,BarHit>::iterator bar_iter;

		// Storage for getters:
		std::map<int,double> Bar_TDC_tdiff;
		std::map<int,double> Bar_TDC_tmean;
		std::map<int,double> Bar_TDC_xpos;
		std::map<int,double> Bar_FADC_tdiff;
		std::map<int,double> Bar_FADC_tmean;
		std::map<int,double> Bar_FADC_xpos;
		std::map<int,double> Bar_FADC_gmadc;
		std::map<int,double> Bar_FADC_gmamp;
		std::map<int,double> Bar_FADC_lnadc;
		std::map<int,double> Bar_FADC_lnamp;
		std::map<int,double> PMT_FADC_adc;
		std::map<int,double> PMT_FADC_amp;
		std::map<int,double> PMT_FADC_time;
		std::map<int,double> PMT_TDC_time;

};

#endif
