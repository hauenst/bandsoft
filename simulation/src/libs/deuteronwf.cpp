#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <complex>
#include <fstream>
#include <cstdarg>
#include <time.h>
#include <signal.h>
#include <iomanip>

using namespace std;

#include "constants.h"
#include "deuteronwf.h"
#include "crossincl.h"

#define NORMFACTOR 11.4081948 //1/HBARC^{3/2}
#define UFACTOR 1./sqrt(4.*PI)  
#define PREC 1.E-03
#define INTLIMIT 10

extern double lambdainput;
extern double betaoffinput;
extern int offshellset;
extern double W20;
extern int symm;

/*********************************************************************
*       Subroutine calculates the deuteron wave function by spin and isospin 
*       components for given deuteron spin projection
*       dspin  = 1,0,-1 - deuteron spin projection
*       proton  = 1,0 - struck out nucleon isospin 1-proton 0 neutron
*       spinp = 1,0 - spin projection of proton 1 up, 0 down
*       spinn = 1,0 - spin projection of neutron 1 up, 0 down
*       p, theta, phi - momentum(GeV/c), polar (rad) and azimuthal(rad) angles of 
*                       relative p-n momentum
*
*	which_wave: which deuteron wave function used: 
*
*       Based on M. Sargsian's fortran code, translated to c++
*	
*       W. Cosyn, Miami 15 march '10
**********************************************************************/

//general deuteron wf [GeV^3/2]
complex<double> deuteronwf_on(int dspin, int proton, int spinp, int spinn, double p, double theta, double phi,
			   int decay, double Gamma, int which_wave){
  
  complex<double> stensor=get_stensor(dspin,spinp,spinn,theta,phi);
  double invsqrt8 = 1./sqrt(8.);
  double x = p*INVHBARC;
  double xt = x*sin(theta);
  double xz = x*cos(theta);
  Gamma*=INVHBARC;
  double *c, *d, *m;
  int c_length=0;  //parametriz array dim
  get_wf_param(&c, &d, &m,c_length, which_wave);
  complex<double> result = Ufront(dspin,spinp,spinn)*uu(x,xt,xz,decay,Gamma,c,m,c_length) + invsqrt8*wd(x,xt,xz,decay,Gamma,d,m,c_length)*stensor;
      
  delete [] c; delete [] d; delete [] m;
  result*=normfactor(which_wave);
  return (proton? result:-result);
}

//if arrays are already present [GeV^3/2]
complex<double> deuteronwf_on(int dspin, int proton, int spinp, int spinn, double p, double theta, double phi,
			    int decay, double Gamma, int which_wave, double *c, double *d, double *m, int c_length){
  complex<double> stensor=get_stensor(dspin,spinp,spinn,theta,phi);
  double invsqrt8 = 1./sqrt(8.);
  double x = p*INVHBARC;
  double xt = x*sin(theta);
  double xz = x*cos(theta);
  Gamma*=INVHBARC;
  complex<double> result = Ufront(dspin,spinp,spinn)*uu(x,xt,xz,decay,Gamma,c,m,c_length) + invsqrt8*wd(x,xt,xz,decay,Gamma,d,m,c_length)*stensor;
  result*=normfactor(which_wave);
  return (proton? result:-result);

}

//off-shell part of the wave function used in FSI amplitude [GeV^3/2]
complex<double> deuteronwf_off(int dspin, int proton, int spinp, int spinn, double p, double theta, double phi,
			       int decay, double Gamma, int which_wave){

  complex<double> stensor=get_stensor(dspin,spinp,spinn,theta,phi);
  double invsqrt8 = 1./sqrt(8.);
  double x = p*INVHBARC;
  double xz = x*cos(theta);
  double xt = x*sin(theta);
  Gamma*=INVHBARC;
  double *c, *d, *m;
  int c_length=0;  //parametriz array dim
  get_wf_param(&c, &d, &m,c_length, which_wave);
  
  double tan2theta=pow(tan(theta),2.);
  complex<double> result= -I*xz*(Ufront(dspin,spinp,spinn)*uu1(x,xt,decay,Gamma,c,m,c_length)+invsqrt8*(wd1(x,xt,decay,Gamma,d,m,c_length)*stensor
					+tan2theta*wd2(xt,d,m,c_length)*(stensor-get_stensor(dspin,spinp,spinn,PI/2.,phi))));
					
  delete [] c; delete [] d; delete [] m;
  result*=normfactor(which_wave);
  return (proton? result:-result);
}

//off-shell part of the wave function used in FSI amplitude, arrays already read in [GeV^3/2]
complex<double> deuteronwf_off(int dspin, int proton, int spinp, int spinn, double p, double theta, double phi,
			    int decay, double Gamma, int which_wave, double *c, double *d, double *m, int c_length){
  
  complex<double> stensor=get_stensor(dspin,spinp,spinn,theta,phi);
  double invsqrt8 = 1./sqrt(8.);
  double x = p*INVHBARC;
  double xz = x*cos(theta);
  if(abs(xz)<1.e-10) return 0.;  
  double xt = x*sin(theta);
  Gamma*=INVHBARC;
  double tan2theta=pow(tan(theta),2.);
  complex<double> result= -I*xz*(Ufront(dspin,spinp,spinn)*uu1(x,xt,decay,Gamma,c,m,c_length)+invsqrt8*(wd1(x,xt,decay,Gamma,d,m,c_length)*stensor
					+tan2theta*wd2(xt,d,m,c_length)*(stensor-get_stensor(dspin,spinp,spinn,PI/2.,phi))));
  result*=normfactor(which_wave);
  return (proton? result:-result);			  
}

//front factor S-wave
double Ufront(int M, int spinp, int spinn){
  double front=0.;
  if(M==1&&spinp==1&&spinn==1) front=1.;
  if(M==-1&&spinp==0&&spinn==0) front=1.;
  if(M==0&&spinp==1&&spinn==0) front=1./sqrt(2.);
  if(M==0&&spinp==0&&spinn==1) front=1./sqrt(2.);
  
  return front;
}


//s-tensor part that goes witht the W part of the wf
complex<double> get_stensor(int dspin, int spinp, int spinn, double theta, double phi){
  
   complex<double> stensor;
   double costheta = cos(theta);
   double sintheta = sin(theta);
   double cosphi = cos(phi);
   double sinphi = sin(phi);
   double cos2phi = cos(2.*phi);
   double sin2phi = sin(2.*phi);
   
   switch(dspin){
    case -1:
      switch(spinp){
	case 0:
	  switch(spinn){
	    case 0:
	      stensor=complex<double>(3.0*costheta*costheta-1.0,0.);
	      break;
	    case 1:
	      stensor=complex<double>(-3.0*costheta*sintheta*cosphi,3.0*costheta*sintheta*sinphi);
	      break;
	    default:
	      cout << "invalid neutron spin" << endl;
	      exit(1);
	  }
	  break;
	case 1:
	  switch(spinn){
	    case 0:
	      stensor=complex<double>(-3.0*costheta*sintheta*cosphi,3.0*costheta*sintheta*sinphi);
	      break;
	    case 1:
	      stensor=complex<double>(3.0*sintheta*sintheta*cos2phi,-3.0*sintheta*sintheta*sin2phi);
	      break;
	    default:
	      cout << "invalid neutron spin" << endl;
	      exit(1);
	  }
	  break;
	default:
	  cout << "invalid proton spin" << endl;
	  exit(1);
      }
      break;
    case 0:
      switch(spinp){
	case 0:
	  switch(spinn){
	    case 0:
	      stensor=complex<double>(-3.0*costheta*sintheta*cosphi,-3.0*costheta*sintheta*sinphi);
	      break;
	    case 1:
	      stensor=complex<double>(-(3.0*costheta*costheta-1.0),0.);
	      break;
	    default:
	      cout << "invalid neutron spin" << endl;
	      exit(1);
	  }
	  break;
	case 1:
	  switch(spinn){
	    case 0:
	      stensor=complex<double>(-(3.0*costheta*costheta-1.0),0.);
	      break;
	    case 1:
	      stensor=complex<double>(3.0*costheta*sintheta*cosphi,3.0*costheta*sintheta*sinphi);
	      break;
	    default:
	      cout << "invalid neutron spin" << endl;
	      exit(1);
	  }
	  break;
	default:
	  cout << "invalid proton spin" << endl;
	  exit(1);
      }
      stensor*=sqrt(2.);
      break;
    case 1:
      switch(spinp){
	case 0:
	  switch(spinn){
	    case 0:
	      stensor=complex<double>(3.0*sintheta*sintheta*cos2phi,3.0*sintheta*sintheta*sin2phi);
	      break;
	    case 1:
	      stensor=complex<double>(3.0*costheta*sintheta*cosphi,3.0*costheta*sintheta*sinphi);
	      break;
	    default:
	      cout << "invalid neutron spin" << endl;
	      exit(1);
	  }
	  break;
	case 1:
	  switch(spinn){
	    case 0:
	      stensor=complex<double>(3.0*costheta*sintheta*cosphi,3.0*costheta*sintheta*sinphi);
	      break;
	    case 1:
	      stensor=complex<double>(3.0*costheta*costheta-1.0,0.);
	      break;
	    default:
	      cout << "invalid neutron spin" << endl;
	      exit(1);
	  }
	  break;
	default:
	  cout << "invalid proton spin" << endl;
	  exit(1);
      }
      break;
    default:
      cout << "invalid deuteron spin" << endl;
      exit(1);
  } 
  return stensor;
  
}

//get parameters for the radial part
//which_wave: 4 diff parametrisations
//0 - Paris wf
//1 - AV18
//2 - CD Bonn
//3 - AV18*

void get_wf_param(double **c, double **d, double **bm, int &c_length, int which_wave){
  
  switch(which_wave){
    case 0:  //Paris wf
    {
      c_length=13;
      (*c) = new double[c_length];
      (*d) = new double[c_length];
      (*bm) = new double[c_length];
      (*c)[0]=0.88688076;                                                   
      (*c)[1]=-0.34717093;                                                  
      (*c)[2]=-3.050238;                                                   
      (*c)[3]=56.207766;                                                    
      (*c)[4]=-749.57334;                                                   
      (*c)[5]=5336.5279;                                                    
      (*c)[6]=-22706.863;                                                   
      (*c)[7]=60434.4690;                                                  
      (*c)[8]=-102920.58 ; 
      (*c)[9]=112233.57;                                                
      (*c)[10]=-75925.226;                                       
      (*c)[11]=29059.715; 
      (*c)[12]=0.;
      for(int i=0;i<c_length-1;i++) (*c)[12]-=(*c)[i];
      
      (*d)[0]=0.023135193;
      (*d)[1]=-0.85604572;
      (*d)[2]=5.6068193  ;
      (*d)[3]=-69.462922 ;
      (*d)[4]=416.31118  ;
      (*d)[5]=-1254.6621 ;
      (*d)[6]=1238.783   ;
      (*d)[7]=3373.9172  ;
      (*d)[8]=-13041.151 ;
      (*d)[9]=19512.524 ;
      
      for(int i=0;i<c_length;i++) (*bm)[i]=0.23162461+i;
      
      double a=0.,b=0.,cc=0.;
      for(int i=0;i<10;i++){
	a+=(*d)[i]/((*bm)[i]*(*bm)[i]);
	b+=(*d)[i];
	cc+=(*d)[i]*((*bm)[i]*(*bm)[i]);
      }
      (*d)[10]=(*bm)[10]*(*bm)[10]/((*bm)[12]*(*bm)[12]-(*bm)[10]*(*bm)[10])/((*bm)[11]*(*bm)[11]-(*bm)[10]*(*bm)[10])
	      *(-(*bm)[11]*(*bm)[11]*(*bm)[12]*(*bm)[12]*a+((*bm)[11]*(*bm)[11]+(*bm)[12]*(*bm)[12])*b-cc);  
      (*d)[11]=(*bm)[11]*(*bm)[11]/((*bm)[10]*(*bm)[10]-(*bm)[11]*(*bm)[11])/((*bm)[12]*(*bm)[12]-(*bm)[11]*(*bm)[11])
	      *(-(*bm)[12]*(*bm)[12]*(*bm)[10]*(*bm)[10]*a+((*bm)[12]*(*bm)[12]+(*bm)[10]*(*bm)[10])*b-cc);  
      (*d)[12]=(*bm)[12]*(*bm)[12]/((*bm)[11]*(*bm)[11]-(*bm)[12]*(*bm)[12])/((*bm)[10]*(*bm)[10]-(*bm)[12]*(*bm)[12])
	      *(-(*bm)[10]*(*bm)[10]*(*bm)[11]*(*bm)[11]*a+((*bm)[10]*(*bm)[10]+(*bm)[11]*(*bm)[11])*b-cc);  
    }
    //for(int i=0;i<13;i++) cout << setprecision(20) << (*c)[i] << " " << (*d)[i] << " " << (*bm)[i] << endl;
      break;
      
    case 1:  //AV18      
    {
      c_length=12;
      (*c) = new double[c_length];
      (*d) = new double[c_length];
      (*bm) = new double[c_length];
      (*c)[0]  =  0.706699E+00;                                                
      (*c)[1]  = -0.169743E+00;                                             
      (*c)[2]  =  0.112368E+01;                                             
      (*c)[3]  = -0.852995E+01;                                             
      (*c)[4]  =  0.195033E+02;                                             
      (*c)[5]  = -0.757831E+02;                                             
      (*c)[6]  =  0.283739E+03;                                             
      (*c)[7]  = -0.694734E+03;                                             
      (*c)[8]  =  0.885257E+03;                                             
      (*c)[9] = -0.720739E+03 ;                                            
      (*c)[10] =  0.412969E+03;
      (*c)[11] = -0.103336E+03;                                             

      (*d)[0]  =  0.176655E-01;                                               
      (*d)[1]  = -0.124551E+00;                                            
      (*d)[2]  = -0.108815E+01;                                            
      (*d)[3]  =  0.384848E+01;                                            
      (*d)[4]  = -0.852442E+01;                                            
      (*d)[5]  =  0.209435E+02;                                            
      (*d)[6]  = -0.490728E+02;                                            
      (*d)[7]  =  0.577382E+02;                                            
      (*d)[8]  = -0.127114E+01;                                            
      (*d)[9] = -0.628361E+02;
      (*d)[10] =  0.581016E+02;
      (*d)[11] = -0.177062E+02;                                            

      (*bm)[0]  = 0.2316; 
      (*bm)[1]  = 1.0;
      (*bm)[2]  = 1.5;
      (*bm)[3]  = 2.0;
      (*bm)[4]  = 2.5;
      (*bm)[5]  = 3.5;
      (*bm)[6]  = 4.5;
      (*bm)[7]  = 5.5;
      (*bm)[8]  = 6.5;
      (*bm)[9] = 8.0;
      (*bm)[10] = 9.5;
      (*bm)[11] = 11.0;
    }
      break;
      
    case 2:  //CD Bonn
    {
      c_length=11;
      (*c) = new double[c_length];
      (*d) = new double[c_length];
      (*bm) = new double[c_length];
      (*c)[0] =   0.88472985;
      (*c)[1] = - 0.26408759;
      (*c)[2] = - 0.44114404e-01;
      (*c)[3] = - 0.14397512e+02;
      (*c)[4] =   0.85591256e+02;
      (*c)[5] = - 0.31876761e+03;
      (*c)[6] =   0.70336701e+03;
      (*c)[7] = - 0.90049586e+03;
      (*c)[8] =   0.66145441e+03;
      (*c)[9]= - 0.25958894e+03;
      (*c)[10]=0.;
      for(int i=0;i<10;i++) (*c)[10]-=(*c)[i];
      
      (*d)[0] =   0.22623762e-01;
      (*d)[1] = - 0.50471056e+00;
      (*d)[2] =   0.56278897e+00;
      (*d)[3] = - 0.16079764e+02;
      (*d)[4] =   0.11126803e+03;
      (*d)[5] = - 0.44667490e+03;
      (*d)[6] =   0.10985907e+04;
      (*d)[7] = - 0.16114995e+04;

      for(int i=0;i<c_length;i++) (*bm)[i]=0.2315380 + 0.9*i;
      
      double tm0=0.,tm1=0.,tm2=0.;
      for(int i=0;i<8;i++){
	tm0+=(*d)[i];
	tm1+=(*d)[i]/((*bm)[i]*(*bm)[i]);
	tm2+=(*d)[i]*((*bm)[i]*(*bm)[i]);
      }
      (*d)[8] = (*bm)[8]*(*bm)[8]/((*bm)[10]*(*bm)[10]-(*bm)[8]*(*bm)[8])/((*bm)[9]*(*bm)[9]-(*bm)[8]*(*bm)[8]) *
              ( -(*bm)[9]*(*bm)[9]*(*bm)[10]*(*bm)[10]*tm1 + ((*bm)[9]*(*bm)[9]+(*bm)[10]*(*bm)[10])*tm0-tm2);
      (*d)[9] = (*bm)[9]*(*bm)[9]/((*bm)[10]*(*bm)[10]-(*bm)[9]*(*bm)[9])/((*bm)[8]*(*bm)[8]-(*bm)[9]*(*bm)[9]) *
              ( -(*bm)[10]*(*bm)[10]*(*bm)[8]*(*bm)[8]*tm1 + ((*bm)[10]*(*bm)[10]+(*bm)[8]*(*bm)[8])*tm0-tm2);
      (*d)[10] = (*bm)[10]*(*bm)[10]/((*bm)[9]*(*bm)[9]-(*bm)[10]*(*bm)[10])/((*bm)[8]*(*bm)[8]-(*bm)[10]*(*bm)[10]) *
              ( -(*bm)[8]*(*bm)[8]*(*bm)[9]*(*bm)[9]*tm1 + ((*bm)[8]*(*bm)[8]+(*bm)[9]*(*bm)[9])*tm0-tm2);
    }
      break;
     
    case 3: // AV18b
    {
      c_length=12;
      (*c) = new double[c_length];
      (*d) = new double[c_length];
      (*bm) = new double[c_length];
      (*bm)[0] = 0.232500e+00;
      (*bm)[1] = 0.500000e+00;
      (*bm)[2] = 0.800000e+00;
      (*bm)[3] = 0.120000e+01;
      (*bm)[4] = 0.160000e+01;
      (*bm)[5] = 0.200000e+01;
      (*bm)[6] = 0.400000e+01;
      (*bm)[7] = 0.600000e+01;
      (*bm)[8] = 0.100000e+02;
      (*bm)[9] = 0.140000e+02;
      (*bm)[10] = 0.180000e+02;
      (*bm)[11] = 0.220000e+02;

      (*c)[0] =  0.105252223e+02;
      (*c)[1] =  0.124352529e+02;
      (*c)[2] = -0.687541641e+02;
      (*c)[3] =  0.239111042e+03;
      (*c)[4] = -0.441014422e+03;
      (*c)[5] =  0.300140328e+03;
      (*c)[6] = -0.230639939e+03;
      (*c)[7] =  0.409671540e+03;
      (*c)[8] = -0.733453611e+03;
      (*c)[9]=  0.123506081e+04;
      (*c)[10]= -0.120520606e+04;
      (*c)[11]=0.;
      for(int i=0;i<11;i++) (*c)[11]-=(*c)[i];

      (*d)[0] =  0.280995496e+00;
      (*d)[1] =  0.334117629e-01;
      (*d)[2] = -0.727192237e+00;
      (*d)[3] = -0.302809607e+01;
      (*d)[4] = -0.903824982e+01;
      (*d)[5] =  0.496045967e+01;
      (*d)[6] = -0.271985613e+02;
      (*d)[7] =  0.125334598e+03;
      (*d)[8] = -0.346742235e+03;
      
      double sp2=0.,sm=0.,sm2=0.;
      for(int i=0;i<9;i++){
	sp2+=(*d)[i]/((*bm)[i]*(*bm)[i]);
	sm+=(*d)[i];
	sm2+=(*d)[i]*((*bm)[i]*(*bm)[i]);
      }
      double a,b,cc;
      a = (*bm)[11]*(*bm)[11];
      b = (*bm)[10]*(*bm)[10];
      cc = (*bm)[9]*(*bm)[9];
      (*d)[9] = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2);
 
      a = (*bm)[9]*(*bm)[9];
      b = (*bm)[11]*(*bm)[11];
      cc = (*bm)[10]*(*bm)[10];
      (*d)[10] = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2);
 
      a = (*bm)[10]*(*bm)[10];
      b = (*bm)[9]*(*bm)[9];
      cc = (*bm)[11]*(*bm)[11];
      (*d)[11] = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2);
      
      double fact=4.*PI*sqrt(PI/2.);
      for(int i=0;i<12;i++){
	(*c)[i]/=fact;
	(*d)[i]/=fact;
      }
    } 
      break;
    default:
      cout << "invalid wf set" << endl;
      exit(1);
  }
  
  
}

//all wf arguments are in fm^-1 (x,xt,xz,Gamma)
//wf is in fm^(3/2)
//Decay part unused for now

//regular U part
double uu(double x, double xt, double xz, int decay, double Gamma,
	  double *c, double *m, int c_length){
 
  double result=0.;
  for(int i=0;i<c_length;i++) result+=c[i]/(x*x+m[i]*m[i]);
  
  if(decay){ 
    double decaypart=0.;
    for(int i=0;i<c_length;i++) decaypart+=c[i]/pow(x*x+m[i]*m[i],2.)*(m[i]*m[i]+xt*xt-xz*xz)/sqrt(m[i]*m[i]+xt*xt);
    result-=Gamma*decaypart;
  }  
  return result;
}

//regular W part
double wd(double x, double xt, double xz, int decay, double Gamma,
	  double *d, double *m, int c_length){
 
  double result=0.;
  for(int i=0;i<c_length;i++) result+=d[i]/(x*x+m[i]*m[i]);
  
  if(decay){ 
    double decaypart=0.;
    for(int i=0;i<c_length;i++) decaypart+=d[i]/pow(x*x+m[i]*m[i],2.)*(m[i]*m[i]+xt*xt-xz*xz)/sqrt(m[i]*m[i]+xt*xt);
    result-=Gamma*decaypart;
  }  
  return result;
}

//off shell U part in FSI
double uu1(double x, double xt, int decay, double Gamma,
	  double *c, double *m, int c_length){
 
  double result=0.;
  for(int i=0;i<c_length;i++) result+=c[i]/(x*x+m[i]*m[i])/sqrt(xt*xt+m[i]*m[i]);
  
  if(decay){ 
    double decaypart=0.;
    for(int i=0;i<c_length;i++) decaypart+=c[i]/pow(x*x+m[i]*m[i],2.);
    result-=2.*Gamma*decaypart;
  }
  return result;
}

//off shell W part in FSI that goes with S(p)
double wd1(double x, double xt, int decay, double Gamma,
	  double *d, double *m, int c_length){
 
  double result=0.;
  for(int i=0;i<c_length;i++) result+=d[i]/(x*x+m[i]*m[i])/sqrt(xt*xt+m[i]*m[i]);
  
  if(decay){ 
    double decaypart=0.;
    for(int i=0;i<c_length;i++) decaypart+=d[i]/pow(x*x+m[i]*m[i],2.);
    result-=2.*Gamma*decaypart;
  }  
  return result;
}

//off-shell W part in FSI that goes wisth S(p)-S(p_perp)
double wd2(double xt, double *d, double *m, int c_length){
 
  double result=0.;
  for(int i=0;i<c_length;i++) result+=d[i]/(m[i]*m[i])/sqrt(xt*xt+m[i]*m[i]);
  return result;
}
  
void error(const char * msg){
  cout << "Error!: " << msg << endl;
  exit(1);
}

//calculates the distorted spectral function
double total_dens(double p, double theta, double phi, double* scattfactor, double s, double nu, double *przprime, double qvec,
		  double massr, double massi, double Er,  double Wx2, double *c, double *d, double *m, int c_length,
		  int decay, double* Gamma, int which_wave, int num_res, double **scattparam, int FSI){
  
  double resultpw=0.;  //plane wave contribution
  double resultcrossq=0.; //cross term contribution
  double resultfsiq=0.; //FSI contribution
  double resulttotalq=0.; //total contribution
  for(int M=-1;M<=1;M++){ //deuteron spin
    for(int spinr=0;spinr<=1;spinr++){ //spectator spin
      complex<double> wave=deuteronwf_on(M, 1, 0, spinr, p, theta, phi, 0, 0., which_wave,c,d,m,c_length); //plane-wave contribution
      complex<double> wavetotalq=wave;
      complex<double> wave2q;
      if(FSI){
	for(int i=0;i<num_res;i++){
	  //fsi contribution, integrated over vec{q_t} = vec{ps_t}-vec{ps'_t}
	  wave2q =I/(32.*PI*PI*qvec*sqrt(Er))*
		  complexromberg(totdens_qt,0.,1.,PREC,3,INTLIMIT,massr,massi,przprime[i],decay,
				  Gamma[i],c,d,m,c_length,scattparam[i], s, nu, 
				  M, spinr, which_wave,Er,p*cos(theta),p*sin(theta),phi,qvec, Wx2)/sqrt(MASSD/(2.*(MASSD-Er)));
	  wavetotalq+=wave2q; //total wave function
	  
	}
      }
      //add squares of wf to spectral function
      resultpw+=real(wave*conj(wave));
      resulttotalq+=real(wavetotalq*conj(wavetotalq));
      resultcrossq+=real(wave*conj(wave2q)+wave2q*conj(wave));
      resultfsiq+=real(wave2q*conj(wave2q));
    }
  }
  //spin factor
  resultpw*=2./3.; 
  resulttotalq*=2./3.;
  resultcrossq*=2./3.;
  resultfsiq*=2./3.;
  return resulttotalq*MASSD/(2.*(MASSD-Er)); //add factor for relativistic normalization
}

//radial integration of FSI amplitude
complex<double> totdens_qt(double qt, va_list ap){
  
  double massr = va_arg(ap,double);
  double massi = va_arg(ap,double);
  double przprime = va_arg(ap,double);
  int decay = va_arg(ap,int);
  double Gamma = va_arg(ap,double);
  double* c = va_arg(ap,double*);
  double* d = va_arg(ap,double*);
  double* m = va_arg(ap,double*);
  int c_length = va_arg(ap,int);
  double *scattparam = va_arg(ap,double*);
  double s = va_arg(ap,double);
  double nu = va_arg(ap,double);
  int M = va_arg(ap,int);
  int spinr = va_arg(ap,int);
  int which_wave = va_arg(ap,int);
  double Er = va_arg(ap,double);
  double prz = va_arg(ap,double);
  double prt = va_arg(ap,double);
  double phi = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double Wx2 = va_arg(ap,double);
 

  return qt*complexromberg(totdens_qphi,0.,2*PI,PREC,3,INTLIMIT,qt,massr,massi,przprime,decay,
			   Gamma,c,d,m,c_length,M,spinr,which_wave,Er,prz,prt,phi,scattparam, s,nu,qvec, Wx2);
  
}

complex<double> totdens_qphi(double qphi, va_list ap){
  
  double qt = va_arg(ap,double);
  double massr = va_arg(ap,double);
  double massi = va_arg(ap,double);
  double przprime = va_arg(ap,double);
  int decay = va_arg(ap,int);
  double Gamma = va_arg(ap,double);
  double* c = va_arg(ap,double*);
  double* d = va_arg(ap,double*);
  double* m = va_arg(ap,double*);
  int c_length = va_arg(ap,int);
  int M = va_arg(ap,int);
  int spinr = va_arg(ap,int);
  int which_wave = va_arg(ap,int);
  double Er = va_arg(ap,double);
  double prz = va_arg(ap,double);
  double prt = va_arg(ap,double);
  double phi = va_arg(ap,double);
  double *scattparam = va_arg(ap,double*);
  double s = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double Wx2 = va_arg(ap,double);
 
  double pt=0.;
  pt=sqrt(prt*prt+qt*qt-2.*prt*qt*cos(qphi-phi));  
  
  //find true pole of the z-integration
  double Wxprime2=-1.;
  W20=nu*nu-qvec*qvec+pow(MASSD-massi,2.)+2*nu*(MASSD-massi); 
  przprime=get_przprime(pt,massr,nu,qvec,Er,Wx2,prz,Wxprime2);
  if(Wxprime2<0.) return 0.;
  
  //construct kinematical variables
  double pprime = sqrt(pt*pt+przprime*przprime);
  double thetaprime=acos(przprime/pprime);
  double phiprime=atan2(prt*sin(phi)-qt*sin(qphi),prt*cos(phi)-qt*cos(qphi));
  double Erprime=sqrt(massr*massr+pprime*pprime);
  double t=(Er-Erprime)*(Er-Erprime)-(prz-przprime)*(prz-przprime)-qt*qt;
   
  //kin factor that goes with rescattering amplitude
  double chi=sqrt(s*s-2.*s*(Wxprime2+massr*massr)+pow(massr*massr-Wxprime2,2.));
  //different off-shell description of rescattering
  double offshellness=0.;
  //mass suppression factor
  if(offshellset==0){
    double massdiff=(MASSD-Erprime)*(MASSD-Erprime)-pprime*pprime-massi*massi+2.*nu*(MASSD-Erprime-massi)+2.*qvec*przprime;
    offshellness=offshell(scattparam[2], massdiff);
  }
  //dipole suppression 
  if(offshellset==1) {
    double onshellm=nu*nu-qvec*qvec+massi*massi+2.*nu*massi;
    offshellness=offshell2(lambdainput, Wxprime2, onshellm);
  }
  //beta_offshell suppression
  if(offshellset==2) offshellness=offshell3(betaoffinput-scattparam[2],t);
  //no off-shell
  if(offshellset==3) offshellness=0.;
  //full off-shell, equal to on-shell amplitude
  if(offshellset==4) offshellness=1.;
  return chi*scatter(t,scattparam)*(deuteronwf_on(M, 1, 0, spinr, pprime, thetaprime, phiprime, decay, Gamma, which_wave,c,d,m,c_length)
	 +(offshellness*deuteronwf_off(M, 1, 0, spinr, pprime, thetaprime, phiprime, decay, Gamma, which_wave,c,d,m,c_length)))
	*sqrt(MASSD/(2.*(MASSD-Erprime)*Erprime));

  
}

//find the true pole, through recursion
double get_przprime(double pt, double massr, double nu, double qvec, double Er, double Wx2, double prz, double& Wxprime2){
  double przprime=0.;
  for(int i=0;i<50;i++){
    double f_Erprime=sqrt(massr*massr+pt*pt+przprime*przprime);
    double f_Wxprime2=nu*nu-qvec*qvec+pow(MASSD-f_Erprime,2.)-pt*pt-przprime*przprime+2.*nu*(MASSD-f_Erprime)+2.*qvec*przprime;
    double f_przprime=prz-(nu+MASSD)/qvec*(Er-f_Erprime);
    if(!symm) if(Wx2>f_Wxprime2) f_przprime-=(Wx2-f_Wxprime2)/(2.*qvec);  //inclusive FSI case
    if((abs((przprime-f_przprime)/f_przprime)<1e-05)) {Wxprime2= f_Wxprime2; return f_przprime;}
    przprime=f_przprime;
    Wxprime2=f_Wxprime2;  
  }
  return przprime;
}

//Born Term density
double born_dens(double p, double Er, int which_wave, double *c, double *d, double *m, int c_length){
  
  double x=p*INVHBARC;
  double result = pow(uu(x,0.,0.,0,0.,c,m,c_length),2.)+pow(wd(x,0.,0.,0,0.,d,m,c_length),2.);
  return result*pow(normfactor(which_wave),2.)*(MASSD/(2.*(MASSD-Er)));
}

//return rescattering amplitude sigma*(I+epsilon)*e^{-beta*t/2}
complex<double> scatter(double t, double *scattparam){

   return scattparam[0]*(I+scattparam[1])*exp(scattparam[2]*t/2.);

}

//different off-shell suppression factors
double offshell(double B, double offshellm){
  
  return exp(B*offshellm);
}
double offshell2(double lambda, double offshellm2, double massX2){
  
  return pow(lambda*lambda-massX2,2.)/(pow(offshellm2-massX2,2.)+pow(lambda*lambda-massX2,2.));
}

double offshell3(double Bdiff, double t){

  return exp(Bdiff*t/2.);

}


//normalization of the different wave function parametrizations
double normfactor(int which_wave){
  double result=1.;
  switch(which_wave){
    case 0:
      result*=0.79788456*UFACTOR;
      break;
    case 1:
      result *=UFACTOR;
      break;
    case 2:
      result /=sqrt(2.)*PI;
      break;
    case 3:
      result *=UFACTOR;
      break;
    default:
      error("invalid wf set");
  }
  result*=NORMFACTOR;

  return result;
}

//integration function, based on Romberg algorithm
double romberger(double (*function)(double, va_list), double a, double b, double acc, int min, int max, ...)
{
  if (abs(a-b)<1e-05) return 0.;
  double *D1=NULL,*D2=NULL;
  double x,h;
  double sum, deviation, solution;
  
  va_list ap;

  D1 = new double[1];
  va_start(ap,max);
  double fa=function(a,ap);
  va_end(ap);
  va_start(ap,max);
  double fb=function(b,ap);
  va_end(ap);

  D1[0] = 0.5*(b-a)*(fa + fb);
  h = b - a;

  for (int n = 1; n <= max ; n++) {
    sum = 0.;
    //store new results
    D2 = new double[n+1];
    int ceiling = power(2,n);
    // Trapezium rule recursive rule
    for (int k = 1; k < ceiling; k += 2) {
      x = a + power(0.5,n)*k*h;
      va_start(ap,max);
      sum += function(x,ap);
      va_end(ap);
    }
    D2[0] = 0.5*D1[0] + power(0.5,n)*h*sum;
    
    // Richarson Interpolation
    int p=4;
    for (int m = 1; m <= n; m++) {
      D2[m] = (p*D2[m-1] - D1[m-1])/(p - 1.);
      p*=4;
    }

    if (n >= min) {
      deviation = (D2[n] == 0.) ? 0. : abs((D2[n]-D1[n-1])/D2[n]);
      if ((deviation < acc )||( abs(D2[n]-D1[n-1]) < acc*1e-07)) {
	solution = D2[n];
	// Free memory
	delete [] D2;
	delete [] D1;
	//cout << "steps " << n << endl;
	//if(abs(solution)>(*estimate)) (*estimate)=abs(solution);
	return solution;
      }
//       if(abs(D2[n])<(*estimate)*1e-05){
// 	solution = D2[n];
// 	delete [] D2;
// 	delete [] D1;
// 	return solution;
//       }
    }
    //D1 will now contain intermediate results
    delete [] D1;
    D1 = D2;
  }
  
  solution = D2[max];
  // Free memory
  delete [] D2;

  //cout << "Too much integrations needed" << endl;
  //cout << deviation << " " << D2[max] << endl << endl;
  va_end(ap);
  return solution;
}


//another integration function, for complex functions
complex<double> complexromberg(complex<double> (*function)(double, va_list), double a, double b, double acc, int min, int max, ...)
{
  if(abs(a-b)<1e-06) return complex<double>(0.,0.);
  complex<double> *D1=NULL,*D2=NULL;
  double x,h;
  complex<double> sum, solution;
  double deviation1;
  
  va_list ap;

  D1 = new complex<double>[1];
  va_start(ap,max);
  complex<double> fa=function(a,ap);
  va_end(ap);
  va_start(ap,max);
  complex<double> fb=function(b,ap);
  va_end(ap);
            
  D1[0] = 0.5*(b-a)*(fa + fb);
  h = b - a;

  for (int n = 1; n <= max ; n++) {
    sum = complex<double>(0.,0.);
    //store new results
    D2 = new complex<double>[n+1];
    int ceiling = power(2,n);
    // Trapezium rule recursive rule
    for (int k = 1; k < ceiling; k += 2) {
      x = a + power(0.5,n)*k*h;
      va_start(ap,max);
      sum += function(x,ap);
      va_end(ap);
    }
    D2[0] = 0.5*D1[0] + power(0.5,n)*h*sum;
    
    // Richarson Interpolation
    int p=4;
    for (int m = 1; m <= n; m++) {
      D2[m] = (double(p)*D2[m-1] - D1[m-1])/(p - 1.);
      p*=4;
    }

    if (n >= min) {
      deviation1 = (D2[n] == complex<double>(0.,0.)) ? 0. : abs((D2[n])-(D1[n-1]))/abs((D2[n]));
      //deviation2 = abs(imag(D2[n])-imag(D1[n-1])/imag(D2[n]));
      if (deviation1 < acc|| abs(D2[n]-D1[n-1]) < acc*1e-05) {
	solution = D2[n];
	// Free memory
	delete [] D2;
	delete [] D1;
	return solution;
      }
    }
    //D1 will now contain intermediate results
    delete [] D1;
    D1 = D2;
  }

  
  solution = D2[max];
  // Free memory
  delete [] D2;

  //cout << "Too much integrations needed" << endl;
  //cout << deviation << " " << D2[max] << endl << endl;
  va_end(ap);
  return solution;
}

//some faster functions than the regular C++ ones
double power(double x, int y){
  switch(y){
  case 0:
    return 1.;
    break;
  case 1:
    return x;
    break;
  case 2:
    return x*x;
    break;
  case 3:
    return x*x*x;
    break;
  default:
    if(y%2==0) {
      double temp = power(x, y/2);
      return temp*temp;
    }
    else{
      double temp = power(x, (y-1)/2);
      return x*temp*temp;
    }
  }
}

/* --------------------------------------------------------------------
! Program to calculate the first kind modified Bessel function
! of integer order N, for any REAL X, using the function BESSI(N,X).
! ---------------------------------------------------------------------
! SAMPLE RUN:
!
! (Calculate Bessel function for N=2, X=0.75).
! 
! Bessel function of order 2 for X =  0.7500:
!
!      Y = 0.073667
!
! ---------------------------------------------------------------------
! Reference: From Numath Library By Tuan Dang Trong in Fortran 77
!            [BIBLI 18].
!
!                         Visual C++ Release 1.0 By J-P Moreau, Paris
! --------------------------------------------------------------------*/


 // ---------------------------------------------------------------------
  double BESSI(int N, double X) {
/*
!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows. 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
*/
      int IACC = 40; 
	  double BIGNO = 1e10, BIGNI = 1e-10;
      double TOX, BIM, BI, BIP, BSI;
      int J, M;

      if (N==0)  return (BESSI0(X));
      if (N==1)  return (BESSI1(X));
      if (X==0.0) return 0.0;

      TOX = 2.0/X;
      BIP = 0.0;
      BI  = 1.0;
      BSI = 0.0;
      M = (int) (2*((N+floor(sqrt(IACC*N)))));
      for (J = M; J>0; J--) {
        BIM = BIP+J*TOX*BI;
        BIP = BI;
        BI  = BIM;
        if (fabs(BI) > BIGNO) {
          BI  = BI*BIGNI;
          BIP = BIP*BIGNI;
          BSI = BSI*BIGNI;
        }
        if (J==N)  BSI = BIP;
      }
      return (BSI*BESSI0(X)/BI);
  }

// ----------------------------------------------------------------------
//  Auxiliary Bessel functions for N=0, N=1
  double BESSI0(double X) {
      double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
      P1=1.0; P2=3.5156229; P3=3.0899424; P4=1.2067429;
      P5=0.2659732; P6=0.360768e-1; P7=0.45813e-2;
      Q1=0.39894228; Q2=0.1328592e-1; Q3=0.225319e-2;
      Q4=-0.157565e-2; Q5=0.916281e-2; Q6=-0.2057706e-1;
      Q7=0.2635537e-1; Q8=-0.1647633e-1; Q9=0.392377e-2;
      if (fabs(X) < 3.75) {
        Y=(X/3.75)*(X/3.75);
        return (P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
      }
      else {
        AX=fabs(X);
        Y=3.75/AX;
        BX=exp(AX)/sqrt(AX);
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        return (AX*BX);
      }
  }

// ---------------------------------------------------------------------
  double BESSI1(double X) {
      double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
      P1=0.5; P2=0.87890594; P3=0.51498869; P4=0.15084934;
      P5=0.2658733e-1; P6=0.301532e-2; P7=0.32411e-3;
      Q1=0.39894228; Q2=-0.3988024e-1; Q3=-0.362018e-2;
      Q4=0.163801e-2; Q5=-0.1031555e-1; Q6=0.2282967e-1;
      Q7=-0.2895312e-1; Q8=0.1787654e-1; Q9=-0.420059e-2;
      if (fabs(X) < 3.75) {
        Y=(X/3.75)*(X/3.75);
        return(X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))));
      }
      else {
        AX=fabs(X);
        Y=3.75/AX;
        BX=exp(AX)/sqrt(AX);
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        return (AX*BX);
      }
  }

// ---------------------------------------------------------------------



void sanitycheck(int calc, int proton, int F_param, int which_wave,int offshellset){
  if((calc<0)||(calc>1)){cout << "calc variable invalid value, should be [0,1]" << endl; exit(1);} 
  if((proton<0)||(proton>1)){cout << "proton variable invalid value, should be [0,1]" << endl; exit(1);} 
  if((F_param<0)||(F_param>2)){cout << "structure function variable invalid value, should be [0,2]" << endl; exit(1);} 
  if((which_wave<0)||(which_wave>3)){cout << "wave function variable invalid value, should be [0,3]" << endl; exit(1);} 
  if((offshellset<0)||(offshellset>4)){cout << "offshellset variable invalid value, should be [0,1]" << endl; exit(1);} 
  return;
}
  










