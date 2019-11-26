#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <complex>
#include <fstream>
#include <cstdarg>
#include <time.h>

using namespace std;

#include "constants.h"
#include "deuteronwf.h"
#include "crossdis.h"
#include "crossincl.h"

//global variables
extern double sigmainput;
extern double betainput;
extern double epsinput;
extern double lambdainput;
extern double betaoffinput;
extern int symm;
extern int F_param;


//fortran prototypes for structure functions
extern"C"{
  void f1f2in09_(double *Z, double *A, double *QSQ, double *Wsq, double *F1, double *F2, double *rc);
}
extern"C"{
  void alekhin_(double *xb,double *q2,double PDFS[22],int *IPAR,int *ICOL);
}

double PREC=1.e-03;
double W20=0.;

//calculates the ratio of F_2D/(F_2P+F_2N) splits contributions in quasi-elastic plane-wave, DIS pw, DIS FSI (disabled for now)
void calc_inclusive2(double &QEpw, double &DISpw, double &fsi1, double &fsi2, double Ein, double Q2,double x,
			double *c, double *d, double *m, int c_length, int which_wave,  int offshellset, int decay, int num_res, int calc){
  
  QEpw=DISpw=fsi1=fsi2=0.;
  double *results=new double[3];
    
  double massi=MASSn;  //average nucleon mass
  double massr=MASSn;
  //electron kinematics
  double nu=Q2/(2.*massi*x);
  double qvec=sqrt(Q2+nu*nu);
  double Eout=Ein-nu;
  double thetae=asin(sqrt(Q2/(4.*Ein*Eout)))*2.;
  double tan2=pow(tan(thetae/2.),2.);
      
    
  double s=MASSD*MASSD-Q2+2.*nu*MASSD; //mandelstam variable
  W20= nu*nu-qvec*qvec+pow(MASSD-massr,2.)+2.*nu*(MASSD-massr); //mass of X (part of it)
  //double x=Q2/(2.*massi*nu);
  double pestimate=0.,thetaestimate=0.; //integration control variables
  double F2N=0.;
  //neutron and proton structure functions
  if(calc!=0){
    double fm2=-Q2+massi*massi+2.*massi*nu;
    if(F_param==0){
      F2N=f2p_b(massi, x,Q2, sqrt(fm2))+f2n_b(massi, x,Q2, sqrt(fm2));
    }
    else{
      if(F_param==1){
	if(fm2>25.) F2N=f2p_b(massi, x,Q2, sqrt(fm2))+f2n_b(massi, x,Q2, sqrt(fm2));
	else{
	  double ONE=1.;
	  double ZERO=0.;
	  double F1=0.,F2=0.,rc=0.;	
	  f1f2in09_(&ONE, &ZERO, &Q2, &fm2,&F1,&F2,&rc);
	  F2N+=F2;
	  f1f2in09_(&ZERO, &ONE, &Q2, &fm2,&F1,&F2,&rc);
	  F2N+=F2;	
	}	
      }
      else{
	int ONE=1;
	int ZERO=0;
	double SF[22];	
	alekhin_(&x,&Q2,SF,&ZERO,&ONE);
	F2N+=SF[1]+SF[4];
      }
    }
  }
  //quasi-elastic plane wave
  QEpw=calc_inclusive_QE_pw(Q2,nu,qvec,tan2, c, d, m, c_length, which_wave)*nu*(calc==0? (1./(2.*MASSD*nu)
	      *(2.*PI*Eout*ALPHA*ALPHA*4.*Eout*Ein*pow(cos(thetae/2.),2.))/(Ein*Q2*Q2)
	      *HBARC*HBARC*1.e07) :2.*PI/F2N);
  if(std::isnan(QEpw)) QEpw=0.;	      
  //loop over proton and neutron
  for(int proton=0;proton<2;proton++){
    //DIS plane-wave and FSI part, integration over spectator momentum
    rombergerN(int_p_incl,1.e-03,1.,3,results,PREC,3,8,&pestimate,Q2,nu,qvec, tan2, x, massi, massr, s, proton,
		      c, d, m, c_length, which_wave, offshellset, decay, num_res, calc, &thetaestimate);
    //frontfactors
    DISpw+=results[0]*2.*massi/MASSD*(calc==0? (1./nu*(2.*PI*Eout*ALPHA*ALPHA*4.*Eout*Ein*pow(cos(thetae/2.),2.))/(Ein*Q2*Q2)
		  *HBARC*HBARC*1.e07):2.*PI/F2N);
    fsi1+=results[1]*2.*massi/MASSD*(calc==0? (1./nu*(2.*PI*Eout*ALPHA*ALPHA*4.*Eout*Ein*pow(cos(thetae/2.),2.))/(Ein*Q2*Q2)
		  *HBARC*HBARC*1.e07):2.*PI/F2N);
    fsi2+=results[2]*2.*massi/MASSD*(calc==0? (1./nu*(2.*PI*Eout*ALPHA*ALPHA*4.*Eout*Ein*pow(cos(thetae/2.),2.))/(Ein*Q2*Q2)
		  *HBARC*HBARC*1.e07):2.*PI/F2N);
    
  }
  //memory management
  delete [] results;
  return;
  
}




//integration over spectator momentum norm
void int_p_incl(double prnorm, double* results, va_list ap){
  
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double tan2 = va_arg(ap,double);
  double x = va_arg(ap,double);
  double massi = va_arg(ap,double);
  double massr = va_arg(ap,double);
  double s = va_arg(ap,double);
  int proton = va_arg(ap,int);
  double* c = va_arg(ap,double*);
  double* d = va_arg(ap,double*);
  double* m = va_arg(ap,double*);
  int c_length = va_arg(ap,int);
  int which_wave = va_arg(ap,int);
  int offshellset = va_arg(ap,int);
  int decay = va_arg(ap,int);
  int num_res = va_arg(ap,int);
  int calc = va_arg(ap,int);
  double* pthetaestimate = va_arg(ap,double*);
  
  //on-shell spectator s
  double Er=sqrt(massr*massr+prnorm*prnorm);
  
  //integration limits for cos(theta_r)
  double lo=-1;
  double hi=1.;
  double lowerlimit= -(-Q2+pow(MASSD-Er,2.)-prnorm*prnorm+2.*(MASSD-Er)*nu)/(2.*qvec*prnorm);
  if(lowerlimit>lo) {lo=lowerlimit;}
  if(hi<lo){ results[0]=results[1]=results[2]=0.; return;}
  //integration over cos(theta_r)
  rombergerN(int_costheta_incl,lo+1e-04,hi-1.e-04,3,results,PREC,3,7,pthetaestimate,prnorm, Er, Q2,nu,qvec, tan2, x, 
				      massi, massr, s, proton, c, d, m, c_length, which_wave,  offshellset, decay, num_res, calc);
  //phase space					
  results[0]*= prnorm*prnorm;
  results[1]*= prnorm*prnorm;
  results[2]*= prnorm*prnorm;
  return;
}

void int_costheta_incl(double costhetar, double* results, va_list ap){
  
  double prnorm = va_arg(ap,double);
  double Er = va_arg(ap,double);
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double tan2 = va_arg(ap,double);  
  double x = va_arg(ap,double);
  double massi = va_arg(ap,double);
  double massr = va_arg(ap,double);
  double s = va_arg(ap,double);
  int proton = va_arg(ap,int);
  double* c = va_arg(ap,double*);
  double* d = va_arg(ap,double*);
  double* m = va_arg(ap,double*);
  int c_length = va_arg(ap,int);
  int which_wave = va_arg(ap,int);
  int offshellset = va_arg(ap,int);
  int decay = va_arg(ap,int);
  int num_res = va_arg(ap,int);
  int calc = va_arg(ap,int);
  
  double thetar=acos(costhetar);
  double structfactor=0.;
  //get structure function factor
  if(calc==0) structfactor = structfunct(x,Q2,nu,qvec,massi,tan2,prnorm,thetar,0., MASSD-Er, Er, proton);
  else structfactor = F2Dincl(x,Q2,nu,qvec,massi,tan2,prnorm,thetar,0., MASSD-Er, Er, proton);
  if(abs(structfactor<1E-09)||std::isnan(structfactor)){results[0]=results[1]=results[2]=0.; return;}  //sanity check

  //plane-wave result
  results[0]= structfactor*born_dens(prnorm,Er,which_wave,c,d,m,c_length);
  //FSI results, disabled for now
  results[1]= 0.;
  results[2]= results[1];
  return;
  
}


//quasi-elastic plane-wave contribution
double calc_inclusive_QE_pw(double Q2,double nu,double qvec, double tan2, double *c, double *d, double *m, int c_length, int which_wave){

  double result=0.;
  double q[4] = {nu,0.,0.,qvec};
    
  
  
  //loop over proton and neutron
  for(int proton=0;proton<2;proton++){
   
    double massi=(proton?MASSP:MASSN);
    double massr=(proton?MASSN:MASSP);
    double  tau = Q2/(4*massi*massi);
    complex<double> ***J;  //Dirac matrix, contains the current operator
    Initialise_Current_Operator(q, Q2, &J, 1.,1.,MUP, proton,tau);
    
    //determine integration limits
    double A=nu+MASSD;
    double B=massi*massi-massr*massr+qvec*qvec-A*A;
    double C=2.*qvec;
    double limit= (4.*A*A*massr*massr-B*B)/(C*C*massr*massr);
    if(limit<0.)  result+=romberger(int_QE_pw,-1.,1.,1.e-2,3,10,Q2,nu,qvec, tan2, massi, massr, proton,
		      c, d, m, c_length, which_wave,J);
    else{ 
      if(limit>1.) result+=0.;
      else{
	double bound=sqrt(limit);
	result+=romberger(int_QE_pw,-1.,-bound-1.e-04,PREC,3,10,Q2,nu,qvec, tan2, massi, massr, proton,
		      c, d, m, c_length, which_wave,J);
	result+=romberger(int_QE_pw,bound+1e-04,1.,PREC,3,10,Q2,nu,qvec, tan2, massi, massr, proton,
		      c, d, m, c_length, which_wave,J);
      }
    }
    //memory management
    for(int i=0;i<4;i++){
      for(int j=0;j<4;j++) delete [] J[i][j];
      delete [] J[i];
    }
    delete [] J;	  
    
  }
 return result; 
}


double int_QE_pw(double costhr, va_list ap){
  
  double Q2 = va_arg(ap,double);
  double nu = va_arg(ap,double);
  double qvec = va_arg(ap,double);
  double tan2 = va_arg(ap,double);
  double massi = va_arg(ap,double);
  double massr = va_arg(ap,double);
  int proton = va_arg(ap,int);
  double* c = va_arg(ap,double*);
  double* d = va_arg(ap,double*);
  double* m = va_arg(ap,double*);
  int c_length = va_arg(ap,int);
  int which_wave = va_arg(ap,int);
  complex<double>*** J = va_arg(ap,complex<double>***);

  //determine norm of pr through quasi-elastic condition
  double A=nu+MASSD;
  double B=massi*massi-massr*massr+qvec*qvec-A*A;
  double C=2.*qvec*costhr;
  double aa=C*C-4.*A*A;
  double bb=2.*B*C;
  double cc=B*B-4.*A*A*massr*massr;
  double discr=4.*abs(A)*sqrt(C*C*massr*massr+cc);
  double pr1=(bb+discr)/(2.*aa);
  double pr2=(bb-discr)/(2.*aa);
  double result=0.,RL=0.,RT=0.;  
  if((pr1>0.)&&(pr1<1.)){ 
    double pf1=sqrt(pr1*pr1+qvec*qvec-C*pr1);
    double costhf1=(qvec*qvec+pf1*pf1-pr1*pr1)/(2.*pf1*qvec);
    double Er1=sqrt(massr*massr+pr1*pr1);
    calc_currentD(RL,RT,pr1,costhr,pf1,costhf1,massi,massr, J, proton, c,d,m,c_length,which_wave);
    result+=pr1*pr1/(Er1*abs(pr1/Er1+(pr1-qvec*costhr)/sqrt(massi*massi+pf1*pf1)))
	*(pow(Q2/(qvec*qvec),2.)*RL+(tan2+Q2/(2.*qvec*qvec))*RT)*(MASSD/(2.*(MASSD-Er1)));
  }
    
  if((pr2>0.)&&(pr2<1.)){ 
    double pf2=sqrt(pr2*pr2+qvec*qvec-C*pr2);
    double Er2=sqrt(massr*massr+pr2*pr2);
    double costhf2=(qvec*qvec+pf2*pf2-pr2*pr2)/(2.*pf2*qvec);
    calc_currentD(RL,RT,pr2,costhr,pf2,costhf2,massi,massr, J,proton, c,d,m,c_length,which_wave);
    result+=pr2*pr2/(Er2*abs(pr2/Er2+(pr2-qvec*costhr)/sqrt(massi*massi+pf2*pf2)))
	*(pow(Q2/(qvec*qvec),2.)*RL+(tan2+Q2/(2.*qvec*qvec))*RT)*(MASSD/(2.*(MASSD-Er2)));
  }
  return result;
  
}

//calculate the quasi-elastic current matrix element
void calc_currentD(double &RL, double &RT, double pr, double costhr, double pf, double costhf, 
		   double massi, double massr, complex<double>*** J, int proton,
		   double *c, double *d, double *m, int c_length, int which_wave){
  
  double sinthr=sin(acos(costhr));
  double thr=acos(costhr);
  double sinthf=-sin(acos(costhr));
  double Er=sqrt(massr*massr+pr*pr);
  double Ei=sqrt(massi*massi+pr*pr);
  double Ef=sqrt(massi*massi+pf*pf);
  double Eioff=MASSD-Er;
  double offshellfactor=(Eioff-Ei)/(2.*massi);
  
  complex<double> uf_up[4] = {1.,0.,-pf*costhf/(Ef+massi),-pf*sinthf/(Ef+massi)};  // \bar{uf}!! spninors without normalization factor sqrt(E+m) (multiplied in the end)
  complex<double> uf_down[4] = {0.,1.,-pf*sinthf/(Ef+massi),pf*costhf/(Ef+massi)};
  
  complex<double> ui_up[4] = {1.,0.,-pr*costhr/(Ei+massi),-pr*sinthr/(Ei+massi)};  // u_i spninors without normalization factor sqrt(E+m) (multiplied in the end)
  complex<double> ui_down[4] = {0.,1.,-pr*sinthr/(Ei+massi),pr*costhr/(Ei+massi)};  // \vec{p}_i = -\vec{p}_r !!!
  
  RL=RT=0.;
  for(int spinf=0; spinf<2; spinf++){
    complex<double> temp_mu0[4] = {0.,0.,0.,0.};
    complex<double> temp_mu1[4] = {0.,0.,0.,0.};
    complex<double> temp_mu2[4] = {0.,0.,0.,0.};
    if(spinf==0){
      for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
	  temp_mu0[i]+=uf_down[j]*J[0][j][i];
	  temp_mu1[i]+=uf_down[j]*J[1][j][i];
	  temp_mu2[i]+=uf_down[j]*J[2][j][i];
	}
      }
    }
    else{
      for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
	  temp_mu0[i]+=uf_up[j]*J[0][j][i];
	  temp_mu1[i]+=uf_up[j]*J[1][j][i];
	  temp_mu2[i]+=uf_up[j]*J[2][j][i];
	}
      }
    }    
    for(int spinr=0; spinr<2; spinr++){
      for(int spinD=-1; spinD<2; spinD++){
	  complex<double> temptensor_mu0=0.,temptensor_mu1=0.,temptensor_mu2=0.;
	  for(int spini=0; spini<2; spini++){
	    complex<double> tempcurrent_mu0=0.,tempcurrent_mu1=0.,tempcurrent_mu2=0.;
	    if(spini==0){
	      for(int i=0;i<4;i++){
		tempcurrent_mu0+=(i<2?-1.:1.)*temp_mu0[i]*ui_down[i];
		tempcurrent_mu1+=(i<2?-1.:1.)*temp_mu1[i]*ui_down[i];
		tempcurrent_mu2+=(i<2?-1.:1.)*temp_mu2[i]*ui_down[i];
	      }
	    }
	    else{
	      for(int i=0;i<4;i++){
		tempcurrent_mu0+=(i<2?-1.:1.)*temp_mu0[i]*ui_up[i];
		tempcurrent_mu1+=(i<2?-1.:1.)*temp_mu1[i]*ui_up[i];
		tempcurrent_mu2+=(i<2?-1.:1.)*temp_mu2[i]*ui_up[i];
	      }
	    }
	    tempcurrent_mu0*=offshellfactor;
	    tempcurrent_mu0*=offshellfactor;
	    tempcurrent_mu0*=offshellfactor;
	    if(spini==0){
	      for(int i=0;i<4;i++){
		tempcurrent_mu0+=temp_mu0[i]*ui_down[i];
		tempcurrent_mu1+=temp_mu1[i]*ui_down[i];
		tempcurrent_mu2+=temp_mu2[i]*ui_down[i];
	      }
	    }
	    else{
	      for(int i=0;i<4;i++){
		tempcurrent_mu0+=temp_mu0[i]*ui_up[i];
		tempcurrent_mu1+=temp_mu1[i]*ui_up[i];
		tempcurrent_mu2+=temp_mu2[i]*ui_up[i];
	      }
	    }
	    complex<double> wf=deuteronwf_on(spinD, proton, proton?spini:spinr, proton?spinr:spini, 
							  pr, thr, 0.,
							  0, 0., which_wave, c, d, m, c_length);
	    temptensor_mu0+=tempcurrent_mu0*wf;
	    temptensor_mu1+=tempcurrent_mu1*wf;
	    temptensor_mu2+=tempcurrent_mu2*wf;
	  }
	  RL+=pow(abs(temptensor_mu0),2.);
	  RT+=pow(abs(temptensor_mu1),2.)+pow(abs(temptensor_mu2),2.);    
      }
    }
  }
  double normfactor=(Ef+massi)*(Ei+massi)/3.;
  RT*=normfactor;
  RL*=normfactor;  
}





//construct the current operator
void Initialise_Current_Operator(double* q, double Q2, 				 
				 complex<double>**** pJ,
				 float luGE,float luGM,float mup, int proton, double tau){
			 //int relat, int isospin,double rs2,double mu_s)

  (*pJ) = new complex<double>**[4];
  for(int i=0;i<4;i++){
    (*pJ)[i] = new complex<double>*[4];
    for(int j=0;j<4;j++) (*pJ)[i][j] = new complex<double>[4];
  }
      
  

  //Q2 in GeV^2
  //q in MeV
  
  //double sin2Weinberg = 0.2224; // sin squared of the Weinberg angle


  double GE_null, GM_null;
  if(proton==1){
    GE_null = 1; 
    GM_null = mup;
  }
  else{
    GE_null = 1.91304; 
    GM_null = -1.91304;
  }
  double  F1, F2, F1acc, F2acc;
  
  if(proton){
    F1 = Get_F1_bba(Q2,luGE,luGM,GE_null,GM_null, tau);
    F2 = Get_F2_bba(Q2,luGE,luGM,GE_null,GM_null, tau);
    F1acc = F1/(2*MASSP);
    F2acc = F2/(2*MASSP);
  }

  else{
    F1 = Get_F1_bba_n(Q2,luGE,luGM,GE_null,GM_null, tau);
    F2 = Get_F2_bba_n(Q2,luGE,luGM,GE_null,GM_null, tau);      
    F1acc = F1/(2*MASSN);
    F2acc = F2/(2*MASSN);
  }
  
  (*pJ)[0][0][0] = F1;
  (*pJ)[0][0][1] = 0;
  (*pJ)[0][0][2] = F2acc*q[3];
  (*pJ)[0][0][3] = F2acc*(q[1]-I*q[2]);
  (*pJ)[0][1][0] = 0;
  (*pJ)[0][1][1] = F1;
  (*pJ)[0][1][2] = F2acc*(q[1]+I*q[2]);
  (*pJ)[0][1][3] = -F2acc*q[3];
  (*pJ)[0][2][0] = (*pJ)[0][0][2];
  (*pJ)[0][2][1] = (*pJ)[0][0][3];
  (*pJ)[0][2][2] = -F1;
  (*pJ)[0][2][3] = 0;
  (*pJ)[0][3][0] = (*pJ)[0][1][2];
  (*pJ)[0][3][1] = (*pJ)[0][1][3];
  (*pJ)[0][3][2] = 0;
  (*pJ)[0][3][3] = -F1;
  
  (*pJ)[1][0][0] = -I*F2acc*q[2];
  (*pJ)[1][0][1] = F2acc*q[3];
  (*pJ)[1][0][2] = 0;
  (*pJ)[1][0][3] = F1+F2acc*q[0];
  (*pJ)[1][1][0] = -F2acc*q[3];
  (*pJ)[1][1][1] = I*F2acc*q[2];
  (*pJ)[1][1][2] = F1+F2acc*q[0];
  (*pJ)[1][1][3] = 0;
  (*pJ)[1][2][0] = 0;
  (*pJ)[1][2][1] = -F1+F2acc*q[0];
  (*pJ)[1][2][2] = (*pJ)[1][0][0];
  (*pJ)[1][2][3] = (*pJ)[1][0][1];
  (*pJ)[1][3][0] = -F1+F2acc*q[0];
  (*pJ)[1][3][1] = 0;
  (*pJ)[1][3][2] = (*pJ)[1][1][0];
  (*pJ)[1][3][3] = (*pJ)[1][1][1];
  
  (*pJ)[2][0][0] = I*F2acc*q[1];
  (*pJ)[2][0][1] = -I*F2acc*q[3];
  (*pJ)[2][0][2] = 0;
  (*pJ)[2][0][3] = -I*(F1+F2acc*q[0]);
  (*pJ)[2][1][0] = -I*F2acc*q[3];
  (*pJ)[2][1][1] = -I*F2acc*q[1];
  (*pJ)[2][1][2] = I*(F1+F2acc*q[0]);
  (*pJ)[2][1][3] = 0;
  (*pJ)[2][2][0] = 0;
  (*pJ)[2][2][1] = I*(F1-F2acc*q[0]);
  (*pJ)[2][2][2] = (*pJ)[2][0][0];
  (*pJ)[2][2][3] = (*pJ)[2][0][1];
  (*pJ)[2][3][0] = -I*(F1-F2acc*q[0]);
  (*pJ)[2][3][1] = 0;
  (*pJ)[2][3][2] = (*pJ)[2][1][0];
  (*pJ)[2][3][3] = (*pJ)[2][1][1];
  
  (*pJ)[3][0][0] = 0;
  (*pJ)[3][0][1] = -F2acc*(q[1]-I*q[2]);
  (*pJ)[3][0][2] = F1+F2acc*q[0];
  (*pJ)[3][0][3] = 0;
  (*pJ)[3][1][0] = F2acc*(q[1]+I*q[2]);
  (*pJ)[3][1][1] = 0;
  (*pJ)[3][1][2] = 0;
  (*pJ)[3][1][3] = -F1-F2acc*q[0];
  (*pJ)[3][2][0] = -F1+F2acc*q[0];
  (*pJ)[3][2][1] = 0;
  (*pJ)[3][2][2] = (*pJ)[3][0][0];
  (*pJ)[3][2][3] = (*pJ)[3][0][1];
  (*pJ)[3][3][0] = 0;
  (*pJ)[3][3][1] = F1-F2acc*q[0];
  (*pJ)[3][3][2] = (*pJ)[3][1][0];
  (*pJ)[3][3][3] = (*pJ)[3][1][1];
  
  
}

//////////////////////////////////////////////////////////////////////////////
/// proton electric form factor of BBA-2003
//////////////////////////////////////////////////////////////////////////////


double Get_GpE_bba(double QmuQmu,float luGE,float GE_pnull)
{
  return GE_pnull*luGE/(1 + 3.253 * QmuQmu + 1.422 * pow(QmuQmu,2) + 0.08582 * pow(QmuQmu,3) 
	    + 0.3318 * pow(QmuQmu,4) +  (-0.09371)*pow(QmuQmu,5) + 
	    0.01076*pow(QmuQmu,6));
}

//////////////////////////////////////////////////////////////////////////////
/// proton magnetic form factor of BBA-2003
//////////////////////////////////////////////////////////////////////////////

double Get_GpM_bba(double QmuQmu,float luGM,float GM_pnull)
{
  return GM_pnull*luGM/(1 + 3.104 *QmuQmu + 1.428 * pow(QmuQmu,2) + 0.1112 * pow(QmuQmu,3)
	       +(-0.006981)*pow(QmuQmu,4) + 0.0003705*pow(QmuQmu,5) + 
	       (-0.000007063)*pow(QmuQmu,6));
}

//////////////////////////////////////////////////////////////////////////////
/// neutron electric form factor of BBA-2003
//////////////////////////////////////////////////////////////////////////////


double Get_GnE_bba(double QmuQmu,float luGE,float GE_nnull, double tau)
{
  return luGE*GE_nnull*pow(1+QmuQmu/.71,-2)*tau*0.942/(1+4.61*tau);
  
}

//////////////////////////////////////////////////////////////////////////////
/// neutron magnetic form factor of BBA-2003
//////////////////////////////////////////////////////////////////////////////

double Get_GnM_bba(double QmuQmu,float luGM,float GM_nnull)
{
  return GM_nnull*luGM/(1 + 3.043 *QmuQmu + 0.8548 * pow(QmuQmu,2) +
		0.6806 * pow(QmuQmu,3)
		+(-0.1287)*pow(QmuQmu,4) + 0.008912*pow(QmuQmu,5));
}


///////////////////////////////////////////////////////////////////////////////
/// Dirac form factor for the proton for BBA-2003
///////////////////////////////////////////////////////////////////////////////

double Get_F1_bba(double QmuQmu,float luGE,float luGM,float GE_pnull,float GM_pnull, double tau)

  
{
  return (tau*Get_GpM_bba(QmuQmu,luGM,GM_pnull)+Get_GpE_bba(QmuQmu,luGE,GE_pnull))/(tau+1);
}

///////////////////////////////////////////////////////////////////////////////
/// Pauli form factor for the proton for BBA-2003
///////////////////////////////////////////////////////////////////////////////

double Get_F2_bba(double QmuQmu,float luGE,float luGM,float GE_pnull,float GM_pnull, double tau)

  
{
  return (Get_GpM_bba(QmuQmu,luGM,GM_pnull)-Get_GpE_bba(QmuQmu,luGE,GE_pnull))/(tau+1);
}

///////////////////////////////////////////////////////////////////////////////
/// Dirac form factor for the neuton for BBA-2003
///////////////////////////////////////////////////////////////////////////////

double Get_F1_bba_n(double QmuQmu,float luGE,float luGM,float GE_nnull,float GM_nnull, double tau)

  
{
  return (tau*Get_GnM_bba(QmuQmu,luGM,GM_nnull)+Get_GnE_bba(QmuQmu,luGE,GE_nnull, tau))/(tau+1);
}

///////////////////////////////////////////////////////////////////////////////
/// Pauli form factor for the neutron for BBA-2003
///////////////////////////////////////////////////////////////////////////////

double Get_F2_bba_n(double QmuQmu,float luGE,float luGM,float GE_nnull,float GM_nnull, double tau)

  
{
  return (Get_GnM_bba(QmuQmu,luGM,GM_nnull)-Get_GnE_bba(QmuQmu,luGE,GE_nnull, tau))/(tau+1);
}

