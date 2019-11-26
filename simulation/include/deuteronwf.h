#ifndef DEUT_H
#define DEUT_H

#include <complex>

std::complex<double> deuteronwf_on(int dspin, int proton, int spinp, int spinn, double p, double theta, double phi,
			   int decay, double Gamma, int which_wave);
std::complex<double> deuteronwf_on(int dspin, int proton, int spinp, int spinn, double p, double theta, double phi,
			    int decay, double Gamma, int which_wave, double *c, double *d, double *m, int c_length);
std::complex<double> deuteronwf_off(int dspin, int proton, int spinp, int spinn, double p, double theta, double phi,
			       int decay, double Gamma, int which_wave);
std::complex<double> deuteronwf_off(int dspin, int proton, int spinp, int spinn, double p, double theta, double phi,
			    int decay, double Gamma, int which_wave, double *c, double *d, double *m, int c_length);			       
double Ufront(int M, int spinp, int spinn);
std::complex<double> get_stensor(int dspin, int spinp, int spinn, double theta, double phi);
void get_wf_param(double **c, double **d, double **bm, int &c_length, int which_wave);
double uu(double x, double xt, double xz, int decay, double Gamma,
	  double *c, double *m, int c_length);
double wd(double x, double xt, double xz, int decay, double Gamma,
	  double *d, double *m, int c_length);
double uu1(double x, double xt, int decay, double Gamma,
	  double *c, double *m, int c_length);
double wd1(double x, double xt, int decay, double Gamma,
	  double *d, double *m, int c_length);
double wd2(double xt,
	  double *d, double *m, int c_length);
void error(const char * msg);
double total_dens(double p, double theta, double phi, double* scattfactor, double s, double nu, double *przprime, double qvec,
		  double massr, double massi, double Er, double Wx2, double *c, double *d, double *m, int c_length,
		  int decay, double *Gamma, int which_wave, int num_res, double **scattparam, int FSI);
std::complex<double> totdens_pt(double phi, va_list ap);
std::complex<double> totdens_phi(double phi, va_list ap);
std::complex<double> totdens_qt(double phi, va_list ap);
std::complex<double> totdens_qphi(double phi, va_list ap);
double get_przprime(double pt, double massr, double nu, double qvec, double Er, double Wx2, double prz, double &Wxprime2);
double born_dens(double p, double Er, int which_wave, double *c, double *d, double *m, int c_length);
void cross_fsi_dens(double p, double theta, double* scattfactor, double s, double nu, double *przprime, double qvec,
		  double massr, double massi, double Er,  double *c, double *d, double *m, int c_length,
		  int decay, double *Gamma, int which_wave, int num_res, double **scattparam, double &crossdens, double &fsidens);
std::complex<double> densint1(double pt, va_list ap);
std::complex<double> densint2(double pt, va_list ap);
std::complex<double> densint3(double pt, va_list ap);
std::complex<double> densint4(double pt, va_list ap);
std::complex<double> scatter(double t, double *scattparam);
double offshell(double B, double offshellm);
double offshell2(double lambda, double massX2, double offshellm2);
double offshell3(double Bdiff, double t);
double normfactor(int which_wave);
std::complex<double> complexromberg(std::complex<double> (*function)(double, va_list), double a, double b, double acc, int min, int max, ...);
double romberger(double (*function)(double, va_list), double a, double b, double acc, int min, int max, ...);
double power(double x, int y);
double BESSI(int N, double X);
double BESSI0(double X);
double BESSI1(double X);
void sanitycheck(int calc, int proton, int F_param, int which_wave,int offshellset);

template <class T> void rombergerN(void (*function)(double, T*,  va_list), double a, double b, 
		       int N, T* results, double acc, int min, int max, double *estimate, ...)
{

  if(abs(a-b)<1e-06){for(int i=0;i<N;i++) results[i]=0.;return;}
  T **DN1= new T*[N];
  T **DN2= new T*[N];
  double x,h=b-a;
  T sum[N], value[N];
  double dev[N];
  
  va_list ap;
  va_start(ap,estimate);
  
  for(int i=0;i<N;i++) {
    DN1[i] = new T[1];
  }

  function(a,value,ap);
  va_end(ap);
  for(int i=0;i<N;i++) DN1[i][0] = value[i]; 
  va_start(ap,estimate);
  function(b,value,ap);
  va_end(ap);
  for(int i=0;i<N;i++){
    DN1[i][0] += value[i]; 
  DN1[i][0] *= 0.5*(b-a);
  }
  

  
  for (int n = 1; n <= max ; n++) {
    for(int i=0;i<N;i++) {
      sum[i] = 0.;
      //store new results
      DN2[i] = new T[n+1];
    }
    int ceiling = power(2,n);
    // Trapezium rule recursive rule
    for (int k = 1; k < ceiling; k += 2) {
      x = a + power(0.5,n)*k*h;
      va_start(ap,estimate);
      function(x,value,ap);
      va_end(ap);
      for(int i=0;i<N;i++) sum[i] += value[i];
    }
    for(int i=0;i<N;i++) DN2[i][0] = 0.5*DN1[i][0] + power(0.5,n)*h*sum[i];
    int p=4;
    for (int m = 1; m <= n; m++) {
      for(int i=0;i<N;i++) DN2[i][m] = (double(p)*DN2[i][m-1] - DN1[i][m-1])/(p - 1.);
      p*=4;
    }
    if (n >= min) {
      double deviation=0.;
      for(int i=0;i<N;i++) {
	dev[i] = (DN2[i][n] == std::complex<double>(0.,0.)) ? 0. : abs((DN2[i][n]-DN1[i][n-1]))/abs(DN2[i][n]);
	if(dev[i]>deviation) deviation=dev[i];
      }
     if (((deviation < acc ) ) || ((abs(DN2[0][n]-DN1[0][n-1]) <acc*1e-07 ))) {
       for(int i=0;i<N;i++) {
	 results[i] = DN2[i][n];
	 delete [] DN2[i];
	 delete [] DN1[i];
	 if(abs(results[i])>(*estimate)) (*estimate)=abs(results[i]);
       }
       delete [] DN1; delete [] DN2;
       return;
     }
     double dummy = 0.;
     for(int i=0;i<N;i++) dummy +=abs(DN2[i][n]);
     if(dummy <(*estimate)*1e-05){
       for(int i=0;i<N;i++){
	 results[i] = DN2[i][n];
	 delete [] DN2[i];
	 delete [] DN1[i];	
       }
       delete [] DN1; delete [] DN2;
       return;
     }
    }
    for(int i=0;i<N;i++){
      delete [] DN1[i];
      DN1[i] = DN2[i];
    }
  }
    
  for(int i=0;i<N;i++){
    results[i] = DN2[i][max];
    delete [] DN2[i];
  }
  delete [] DN1; delete [] DN2;
  va_end(ap);
  return;  
  
}


// double cross_dens(double p, double theta, double *scattfactor, double *przprime, double qvec,
// 		  double massr, double massi, double Er,  double *c, double *d, double *m, int c_length,
// 		  int decay, double *Gamma, int which_wave, int num_res, double **scattparam);
// std::complex<double> crossterm1(double pt, va_list ap);
// std::complex<double> crossterm2(double pt, va_list ap);
// std::complex<double> crossterm3(double pt, va_list ap);
// std::complex<double> crossterm4(double pt, va_list ap);
// double FSI_dens(double p, double theta, double* scattfactor, double *przprime, double qvec,
// 		double massr, double massi, double Er,  double *c, double *d, double *m, int c_length,
// 		int decay, double *Gamma, int which_wave, int num_res, double **scattparam);
// std::complex<double> intU(double pt, va_list ap);
// std::complex<double> intW1(double pt, va_list ap);
// std::complex<double> intW2(double pt, va_list ap);
// std::complex<double> intW3(double pt, va_list ap);
#endif
