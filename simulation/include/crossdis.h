#ifndef CROSS_H
#define CROSS_H

void neutronstructurefit(double Ein, double Q2, double Wprime, double pr,
		  int proton, int which_wave, int decay, int num_res, int which_calc, int offshell, int varvalue, double xref);
double calc_deeps(double Ein, double Q2, double Wprime, double pr, double thetar, 
		  int proton, int which_wave, int decay, int num_res, int FSI);
double calc_cross(double Ein, double Q2, double x, double pr, double thetar, double phir, 
		  int proton, int which_wave, int decay, int num_res, int FSI);
double structfunct(double x, double Q2,double nu,double qvec,double massi,double tan,
		   double pr,double thetar,double phir,double Einoff, double Er, int proton);
double F2Dincl(double x, double Q2,double nu,double qvec,double massi,double tan2,
		   double pr,double thetar,double phir,double Einoff, double Er, int proton);
void calc_F2D(double x, double Q2,double nu,double qvec,double massi,double tan2,
		   double pr,double thetar,double Einoff, double Er, int proton, int offshell, double &ff, double &F2, double &F2D);
		   double f2p_b(double massi, double xp, double q2, double wstar2);
double f2n_b(double massi, double xp, double q2, double wstar2);
double bodek(double wm, double qsq);
void init_param(double *scattfactor, double *przprime, double prz,  
		double **scatterparam, double *Gamma, int num_res,
		double nu, double qvec, double Q2, double Wx2, double Wxprime2, double Wxprime2_0, double massi, double Er, double s);
double intcheck(double p, va_list ap);
double phiint(double phir, va_list ap);
void readin_deeps(double ******pdeepsarray, std::string dir);
void maint_deepsarray(double *****deepsarray);
void getF_deeps(double *****deepsarray,int proton, double Ein, double residu, int num_res, int which_wave, int offshell, int decay, int varvalue);
//void neutronstructurefit_old(double Ein, double Q2, double x, 
//		  int proton, int which_wave, int decay, int num_res, int offshell, int fitvar,int varvalue);
#endif
