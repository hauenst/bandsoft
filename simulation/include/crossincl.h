#ifndef CROSSINCL_H
#define CROSSINCL_H

void calc_inclusive2(double &QEpw, double &DISpw, double &fsi1, double &fsi2, double Ein, double Q2,double x,
			double *c, double *d, double *m, int c_length, int which_wave,  int offshellset, int decay, int num_res, int calc);
void calc_inclusive_DIS(double &pw, double &fsi, double &fsi2, double Q2,double nu,double qvec, double tan2, 
			double *c, double *d, double *m, int c_length, int which_wave, int offshellset, int decay, int num_res);
void int_p_incl(double prnorm, double* results, va_list ap);
void int_costheta_incl(double costhetar, double* results, va_list ap);
void int_qt(double qt, double* qresults, va_list ap);
void int_qphi(double qphi, double* qresults, va_list ap);
double get_przprimebis(double pt, double massr, double nu, double qvec, double Er, double W1sq, double prz, double& W2sq);
void init_scattparam(double *scattparam, double W1sq, double Q2);
void getfitscattparam(double *scattparam,double W2sq, double Q2);

double calc_inclusive_QE_pw(double Q2,double nu,double qvec, double tan2, double *c, double *d, double *m, int c_length, int which_wave);
double int_QE_pw(double costh1, va_list ap);
void calc_currentD(double &RL, double &RT, double pr, double costhr, double pf, double costhf, 
		   double massi, double massr, complex<double>*** J, int proton,
		   double *c, double *d, double *m, int c_length, int which_wave);

void Initialise_Current_Operator(double* q, double Q2,complex<double>**** pJ,float luGE,float luGM,float mup, int proton, double tau);
double Get_GpE_bba(double QmuQmu,float luGE,float GE_pnull);
double Get_GpM_bba(double QmuQmu,float luGM,float GM_pnull);
double Get_GnE_bba(double QmuQmu,float luGE,float GE_nnull, double tau);
double Get_GnM_bba(double QmuQmu,float luGM,float GM_nnull);
double Get_F1_bba(double QmuQmu,float luGE,float luGM,float GE_pnull,float GM_pnull, double tau);
double Get_F2_bba(double QmuQmu,float luGE,float luGM,float GE_pnull,float GM_pnull, double tau);
double Get_F1_bba_n(double QmuQmu,float luGE,float luGM,float GE_nnull,float GM_nnull, double tau);
double Get_F2_bba_n(double QmuQmu,float luGE,float luGM,float GE_nnull,float GM_nnull, double tau);

#endif 