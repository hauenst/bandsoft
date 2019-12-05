#ifndef __HELPERS_H__
#define __HELPERS_H__

#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
using namespace std;

void zoomTH2F(TH2F * h2, double zoomAround, double zoomLower, double zoomUpper ) ;
double getTriggerPhase( long timeStamp );
int doProj( TH2F * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE, int thres, int lastBin );
void drawTH2F( TH2F * hist , TCanvas * c , int cd );

#endif
