#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <cmath>

const double alpha=0.00729735253;

// Units of GeV
const double mP = 0.93827231;
const double mN = 0.93956536;
const double mD = 1.8756;
const double Ebeam = 10.3;
const double mE = 0.0005109989;

// Units of cm / ns
const double cAir = 29.9792458;

// LAD constants
const double lad_max_phi=60.*M_PI/180.;
const double lad_min_theta_deg=90.;
const double lad_max_theta_deg=165.;
const double lad_target_z = 25.;
const double lad_angles[3]={102.*M_PI/180.,127.*M_PI/180.,152.*M_PI/180.};
const double lad_radii[3]={564.,472.,564.}; // cm
const double lad_gem_radii[2]={100.,150.};
const double hallc_phi_range=24.*M_PI/180.;
const double hallc_min_theta=8.*M_PI/180.;
const double hallc_max_theta=24.*M_PI/180.;
const double hms_acc_theta=1.6*M_PI/180.;
const double shms_acc_theta=1.4*M_PI/180.;
const double hms_acc_mom=1.1;
const double shms_acc_mom_hi=1.22;
const double shms_acc_mom_lo=1.1;
const double hms_acc=0.006; // 6msr
const double shms_acc=0.005; // 5msr


#endif
