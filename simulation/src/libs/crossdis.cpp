#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <complex>
#include <fstream>
#include <cstdarg>
#include <time.h>
#include <signal.h>
#include <string>

using namespace std;

#include "constants.h"
#include "deuteronwf.h"
#include "crossdis.h"

//global variables
extern double sigmainput;
extern double betainput;
extern double epsinput;

extern int phiavg;
extern int F_param;

//fortran linked functions for structure functions
extern"C"{
    void f1f2in09_(double *Z, double *A, double *QSQ, double *Wsq, double *F1, double *F2, double *rc);
}

extern"C"{
  void alekhin_(double *xb,double *q2,double PDFS[22],int *IPAR,int *ICOL);
}


//Deeps data calculations, binned in Q^2, pr and Wprime
//output is F_2N*P along the lines of the Deeps article.  We calc cross section averaged over phi and divide by the Deeps prefactor.
double calc_deeps(double Ein, double Q2, double Wprime, double pr, double costhetar, 
		  int proton, int which_wave, int decay, int num_res, int FSI){
  
  //mass of nucleon interacting with photon
  double massi=proton? MASSP:MASSN;
  //mass of spectator nucleon
  double massr=proton? MASSN:MASSP;
  
  double Er=sqrt(massr*massr+pr*pr); //spectator on-shell energy
  double Einoff=MASSD-Er; //DIS nucleon off-shell energy
  double massoff=sqrt(Einoff*Einoff-pr*pr); //off-shell mass DIS nucleon
  
  double xprime=Q2/(Wprime*Wprime-massoff*massoff+Q2); //x'=Q^2/(2p_i*q)
  
  //calc photon energy
  double prz=pr*costhetar;
  double aaa=Einoff*Einoff-prz*prz;
  double bbb=-Einoff*Q2/xprime;
  double ccc=Q2*Q2/(4.*xprime*xprime)-Q2*prz*prz;
  
  double discr=sqrt(bbb*bbb-4.*aaa*ccc);
  double nu1=(-bbb+discr)/(2.*aaa);
  double nu2=(-bbb-discr)/(2.*aaa);
  double nu=nu2;
  if(costhetar<0.) nu=nu1;
  double qvec=sqrt(Q2+nu*nu);
  double x=Q2/(2.*massi*nu); //Bjorken x
  double Eout=Ein-nu;
  double thetae=asin(sqrt(Q2/(4.*Ein*Eout)))*2.; //electron angle
  double y=nu/Ein;
  
  
  //variables needed for deeps prefactor
  double thetain=acos((Ein*Ein+qvec*qvec-Eout*Eout)/(2.*Ein*qvec));
  double yprime=(Einoff*nu+prz*qvec)/(Ein*Einoff+Ein*cos(thetain)*pr*costhetar);
   
  double tan2=pow(tan(thetae/2.),2.); //tan^2(theta_e/2)
  
  double front=(4.*PI*ALPHA*ALPHA)/(x*Q2*Q2)*(1-y-(x*x*y*y*massi*massi)/Q2); //kinematical front factor
  double R=0.18; 
  double frontdeeps=(4.*PI*ALPHA*ALPHA)/(xprime*Q2*Q2)*
      (yprime*yprime/(2.*(1.+R))+(1.-yprime)+(pow(massoff*xprime*yprime,2.)*(1.-R))/(Q2*(1+R))); //deeps kin front factor
  double dxprimedx=-2.*xprime*xprime*nu/(x*qvec)*((Einoff+prz)/(nu-qvec)+1./(2.*xprime)); //(dx'/dx) jacobian
  double Wxprime2=massi*massi+2.*nu*massi-Q2; //(m_i+q)^2 invariant mass of X' created on stationary nucleon
  double Wxprime2_0=massr*massr+MASSD*MASSD-Q2+2.*nu*MASSD; //part of real invariant mass 
  double Wx2=Wprime*Wprime;
  double s=MASSD*MASSD-Q2+2.*nu*MASSD; //mandelstam variable
  
  
  double thetar=acos(costhetar);
  //get structure function part
  double structfactor=structfunct(x,Q2,nu,qvec,massi,tan2,pr,thetar,0., Einoff, Er, proton);
  if(abs(structfactor<1E-09)||std::isnan(structfactor)) return structfactor;  //sanity check
  
  //read in arrays with wave-functions
  double *c, *d, *m;
  int c_length=0;  //parametriz array dim
  get_wf_param(&c, &d, &m,c_length, which_wave);

  //construct some parameters for the FSI part
  double scattfactor[num_res];
  double przprime[num_res];
  double Gamma[num_res];
  double **scatterparam=new double*[num_res];
  for(int i=0;i<num_res;i++){
    scatterparam[i]=new double[3];
  }  
  init_param(scattfactor,przprime,pr*cos(thetar),scatterparam,Gamma,num_res,
	     nu, qvec, Q2, Wx2, Wxprime2, Wxprime2_0, massr, Er, s);
    
  //calc spectral function	      
  double dens=total_dens(pr,thetar,0.,scattfactor,s,nu,przprime,qvec,massr,massi,Er, Wx2,
			  c,d,m,c_length,decay,Gamma,which_wave,num_res,scatterparam, FSI);

  //memory cleanup			  
  delete [] c; delete [] d; delete [] m;
  for(int i=0;i<num_res;i++){
   delete [] scatterparam[i];
  }
  delete [] scatterparam;
  return front*structfactor*dens/frontdeeps/dxprimedx;  //F_2D*P [GeV-3]
}


double calc_cross(double Ein, double Q2, double x, double pr, double thetar, double phir, 
		  int proton, int which_wave, int decay, int num_res, int FSI){
  
  
  //mass of nucleon interacting with photon
  double massi=proton? MASSP:MASSN;
  //mass of spectator nucleon
  double massr=proton? MASSN:MASSP;
  
  
  //electron kinematics
  double nu=Q2/(2.*massi*x);
  double Eout=Ein-nu;
  double thetae=asin(sqrt(Q2/(4.*Ein*Eout)))*2.;
  double qvec=sqrt(Q2+nu*nu);
  //double x=Q2/(2.*massi*nu);
  double y=nu/Ein;
  
   
  double tan2=pow(tan(thetae/2.),2.);
  double Er=sqrt(massr*massr+pr*pr); //spectator on-shell energy
  double Einoff=MASSD-Er; //DIS nucleon off-shell energy
  
  double front=(4.*PI*ALPHA*ALPHA)/(x*Q2*Q2)*(1-y-(x*x*y*y*massi*massi)/Q2); //kinematical front factor
  
  
   double Wxprime2=massi*massi+2.*nu*massi-Q2;//initial X mass (quasi-free kinematics)
  double Wxprime2_0=massr*massr+MASSD*MASSD-Q2+2.*nu*MASSD; //initial X mass (real kinematics), part of it
  double Wx2=MASSD*MASSD-Q2+massr*massr+2.*MASSD*(nu-Er)-2.*nu*Er+2.*qvec*pr*cos(thetar); //final X mass
  //cout << Wxprime2 << " " << Wx2 << endl;
  double s=MASSD*MASSD-Q2+2.*nu*MASSD; //mandelstam variable
  
  
  
  double structfactor=structfunct(x,Q2,nu,qvec,massi,tan2,pr,thetar,phir, Einoff, Er, proton);  //calc structure functions
  
  //wave function arrays read-in
  double *c, *d, *m;
  int c_length=0;  //parametriz array dim
  get_wf_param(&c, &d, &m,c_length, which_wave);

  //fsi parameters initialization
  double scattfactor[num_res];
  double przprime[num_res];
  double Gamma[num_res];
  //double **offshellparam=new double*[num_res];
  double **scatterparam=new double*[num_res];
  for(int i=0;i<num_res;i++){
    //offshellparam[i]=new double[2];
    scatterparam[i]=new double[3];
  }
  
  init_param(scattfactor,przprime,pr*cos(thetar),scatterparam,Gamma,num_res,
	     nu, qvec, Q2, Wx2, Wxprime2, Wxprime2_0, massr, Er, s);
  
  //calc spectral function	     
  double dens=total_dens(pr,thetar,phir,scattfactor,s,nu,przprime,qvec,massr,massi,Er,Wx2,
			  c,d,m,c_length,decay,Gamma,which_wave,num_res,scatterparam, FSI);

  //memory management			  
  delete [] c; delete [] d; delete [] m;
  for(int i=0;i<num_res;i++){
   delete [] scatterparam[i];
  }
  delete [] scatterparam;
//   cout << thetar*RADTODEGR << " " << dens*Er << endl;
  return front*structfactor*dens*Er*HBARC*HBARC*1.e07;  //nanobarn/GeV^4
}


//calculates the structure function part of the diff cross section in the semi-inclusive case
double structfunct(double x, double Q2,double nu,double qvec,double massi,double tan2,
		   double pr,double thetar,double phir,double Einoff, double Er, int proton){
  
  double alphaq=(nu-qvec)/(MASSD/2.);
  double alphai=(Einoff+pr*cos(thetar))/(MASSD/2.);  //vec{pi}=-vec{pr} in PW!
  double pt=pr*sin(thetar);
  double cosdelta=nu/qvec;
  double sindelta2=Q2/(qvec*qvec);
  double piq=(Einoff*nu+qvec*pr*cos(thetar));  //vec{pi}=-vec{pr} in PW!
  double nutilde=piq/massi;
  double wstar2 = Einoff*Einoff-pr*pr; //effective mass off-shell nucleon
  double xtilde=Q2/(2*piq);
  double nuoffshell=(wstar2-massi*massi+2.*piq)/(2.*massi); //(m_i+qoffshell)^2=(p_i+q)^2
  double xoffshell=Q2/(2.*massi*nuoffshell);
  double fm=sqrt(wstar2+2.*piq-Q2);//sqrt((p_i+q)^2)
  //if(std::isnan(fm)) cout << wstar2 << " " << wstar2+2.*piq-Q2 << endl;
  double F2=proton? f2p_b(massi, xoffshell,Q2, fm):f2n_b(massi, xoffshell,Q2, fm);
  double R=0.18;
  double F1=F2*2.*xtilde/(1+R)*(pow(alphai/alphaq+1/(2.*xtilde),2.)-pt*pt/(2.*Q2)*R);  //for moving nucleon
      
  //Christy&Bosted F2,F1
  if(F_param==1&&fm<5.){
    double F1_CB,F2_CB, R_CB;
    double ONE=1.;
    double ZERO=0.;
    double Wsq=fm*fm;
    if(proton) f1f2in09_(&ONE, &ZERO, &Q2, &Wsq, &F1_CB, &F2_CB, &R_CB);
    else f1f2in09_(&ZERO, &ONE, &Q2, &Wsq, &F1_CB, &F2_CB, &R_CB);
    F1=F1_CB;
    F2=F2_CB;
  }
  //leading twist F2, F1
  if(F_param==2){
    int ONE=1;
    int ZERO=0;
    double SF[22];	
    alekhin_(&xtilde,&Q2,SF,&ZERO,&ONE);
    if(proton){
      F1=SF[7]/(2.*xtilde);
      F2=SF[1];
    }
    else{
      F1=SF[10]/(2.*xtilde);
      F2=SF[4];
    }
  }    
  //cout << F_param << " " << F1 << " " << F2 << endl;
  
  double FL=pow((alphai+alphaq*piq/Q2)*(1+cosdelta),2.)*nu/nutilde*F2-nu/massi*sindelta2*F1;
  double FT=2.*F1+pt*pt/(massi*nutilde)*F2;
  double FTT=nu*pt*pt*sindelta2/(nutilde*massi*massi*2.)*F2;
  double FTL=2.*(1.+cosdelta)*pt*nu/(massi*nutilde)*(alphai+alphaq*piq/Q2)*F2;
  //not averaged over phi, all structure functions survive
  if(phiavg==0) return (FL+(Q2/(2.*qvec*qvec)+tan2)*nu/massi*FT)+sqrt(sindelta2+tan2)*cos(phir)*FTL+cos(2.*phir)*FTT;
  //averaged over phi
  else{
    if(!F_param) return F2*(pow(alphai/alphaq+1/(2.*xtilde),2.)+pt*pt/(2.*Q2))*2.*xtilde*nu/massi*(sindelta2+2.*tan2/(1+R));
    else return FL+(Q2/(2.*qvec*qvec)+tan2)*nu/massi*FT;
  }
}

//calculates the structure function part of the diff cross section in the inclusive case
double F2Dincl(double x, double Q2,double nu,double qvec,double massi,double tan2,
		   double pr,double thetar,double phir,double Einoff, double Er, int proton){
  
  double alphaq=(nu-qvec)/(MASSD/2.);
  double alphai=(Einoff+pr*cos(thetar))/(MASSD/2.);  //vec{pi}=-vec{pr} in PW!
  double pt=pr*sin(thetar);
  double cosdelta=nu/qvec;
  double sindelta2=Q2/(qvec*qvec);
  double piq=(Einoff*nu+qvec*pr*cos(thetar));  //vec{pi}=-vec{pr} in PW!
  double nutilde=piq/massi;
  double wstar2 = Einoff*Einoff-pr*pr; //effective mass off-shell nucleon
  double xtilde=Q2/(2*piq);
  double nuoffshell=(wstar2-massi*massi+2.*piq)/(2.*massi); //(m_i+qoffshell)^2=(p_i+q)^2
  double xoffshell=Q2/(2.*massi*nuoffshell);
  double fm=sqrt(wstar2+2.*piq-Q2);//sqrt((p_i+q)^2)
  double F2=proton? f2p_b(massi, xoffshell,Q2, fm):f2n_b(massi, xoffshell,Q2, fm);
  double R=0.18;
  double F1=F2*2.*xtilde/(1+R)*(pow(alphai/alphaq+1/(2.*xtilde),2.)-pt*pt/(2.*Q2)*R);
  
  //Christy&Bosted F2,F1
  if(F_param&&fm<5.){
    double F1_CB,F2_CB, R_CB;
    double ONE=1.;
    double ZERO=0.;
    double Wsq=fm*fm;
    if(proton) f1f2in09_(&ONE, &ZERO, &Q2, &Wsq, &F1_CB, &F2_CB, &R_CB);
    else f1f2in09_(&ZERO, &ONE, &Q2, &Wsq, &F1_CB, &F2_CB, &R_CB);
    F1=F1_CB;
    F2=F2_CB;
  }
  //leading twist F2, F1
  if(F_param==2){
    int ONE=1;
    int ZERO=0;
    double SF[22];	
    alekhin_(&xtilde,&Q2,SF,&ZERO,&ONE);
    if(proton){
      F1=SF[7]/(2.*xtilde);
      F2=SF[1];
    }
    else{
      F1=SF[10]/(2.*xtilde);
      F2=SF[4];
    }    
  }    

  
  double FL=pow((alphai+alphaq*piq/Q2)*(1+cosdelta),2.)*nu/nutilde*F2-nu/massi*sindelta2*F1;
  double FT=2.*F1+pt*pt/(massi*nutilde)*F2;
  
   return FL+Q2/(2.*qvec*qvec)*nu/massi*FT;
  
}



//	Bodek's Parameterization of Proton's F2
double f2p_b(double massi, double xp, double q2, double fm){

  double pmo=massi;
  double f2p,yn,wwi,t,gw;
  //f2p =  0.0;
  yn  =  q2/(2*massi*xp);//nuoffshell
  //fm2 = -q2 +2.0*massi*yn + wstar2;
  if ((xp > 1.0) | (fm < massi)){return 0;}
  //	fm=pow(fm2,.5);
  //fm  = sqrt(fm2);
  //cout << fm;
  wwi = (2.*yn*pmo+1.642)/(q2 +0.376);
  t = 1.0 - 1.0/wwi;
  gw=0.256*pow(t,3)+2.178*pow(t,4)+
        0.898*pow(t,5)-6.716*pow(t,6)+3.756*pow(t,7);
  f2p = bodek(fm,q2)*gw*wwi*xp;
  return f2p;
}

// Bodek's parameterization of neutron's F2 
double f2n_b(double massi, double xp, double q2, double fm){
  
  
  double pmo=massi;
  double f2n,yn,wwi,t,gw;
  //f2n =  0.0;
  yn  =  q2/(2*massi*xp);
  //fm2 = -q2 + 2.0*massi*yn + wstar2;
  if ((xp > 1.0) | (fm < massi)){return 0;}
  //  fm=pow(fm2,2);
  //fm = sqrt(fm2);
  //cout << fm;
  wwi = (2*yn*pmo + 1.642)/(q2 + 0.376);
  t   = 1.0-1.0/wwi;
  gw  = 0.064*pow(t,3) + 0.225*pow(t,4)+
  4.106*pow(t,5) - 7.079*pow(t,6)+3.055*pow(t,7);
  f2n = bodek(fm,q2)*gw*wwi*xp;
  return f2n;
}


/************************************************/
/*		BODEK PARAMETRIZATION		*/
/************************************************/

double bodek(double wm, double qsq)
{
 const double c[25]={0., 1.0741163,  0.75531124, 3.3506491   , 1.7447015  ,
                    3.5102405,  1.040004  , 1.2299128   , 0.10625394 ,
                    0.48132786, 1.5101467 , 0.081661975 , 0.65587179 ,
                    1.7176216 , 0.12551987, 0.7473379   , 1.953819   ,
                    0.19891522,-0.17498537, 0.0096701919,-0.035256748, 
		    3.5185207 ,-0.59993696, 4.7615828   , 0.41167589}; 
        const int lspin[5]={0,1,2,3,2};
	      int nres=4,nbkg=5,index,i;
	
	double pmsq=.880324,pm2=1.876512,pm=0.938256;
	double b,b1,b2,wsq,omega,x,xpx,piemsq;
	double eb1,eb2,bbkg,bres,ressum;
	double ram,rma,rwd,qstarn,qstaro,j,k;
	double term,termo,gamres,brwig,res;
	b=0;

	//  cout <<  "c2 " << c[2]  <<  endl;

	if (wm < .94){ return b;}
	wsq=wm*wm;
	omega=1+(wsq-pmsq)/qsq;
	x=1/omega;
	xpx=c[22]+c[23]*(x-c[24])*(x-c[24]);
	//  cout <<  "xpx " << xpx  <<  endl;

	piemsq=(c[1]-pm)*(c[1]-pm);
	b1=0;
	if (wm == c[1]){goto label11;}
	b1=max(0. , (wm - c[1]))/(wm-c[1])*c[2];
	//  cout <<  "b1 " << b1  <<  endl;
	
label11:
	eb1=c[3]*(wm-c[1]);
	if (eb1 > 25.){goto label1;}
	b1=b1*(1-exp(-eb1));
	b2=0;
	if (wm == c[4]){goto label12;}
label1:
	b2=max(0., (wm - c[4]))/(wm-c[4])*(1.-c[2]);
	//  cout <<  "b2 " << b2  <<  endl;

label12:
	eb2=c[5]*(wsq-c[4]*c[4]);
	if (eb2 > 25.){goto label2;}
	b2=b2*(1.-exp(-eb2));
	//  cout <<  "b2 " << b2  <<  endl;

label2:
	bbkg=b1+b2;
	bres=c[2]+b2;
	ressum=0.;
	for(i=1;i<=nres;i++)
	{
		index=(i-1)*3+1+nbkg;
		ram=c[index];
		if (i==1) {ram=c[index]+c[18]*qsq+c[19]*qsq*qsq;}
		rma=c[index+1];
		if (i==3) {rma=rma*(1+c[20]/(1+c[21]*qsq));}
		rwd=c[index+2];
		qstarn=sqrt(max(0,pow((wsq+pmsq-piemsq)/(2*wm),2)-pmsq));
		qstaro=sqrt(max(0,pow((rma*rma-pmsq+piemsq)/(2*rma),2)-piemsq));
		if (qstaro <= 1E-10){goto label40;}
		term=6.08974*qstarn;
		termo=6.08974*qstaro;
		j=2*lspin[i];
		k=j+1;
		gamres=rwd*pow((term/termo),k)*(1+pow(termo,j))/(1+pow(term,j));
		gamres=gamres/2;
		brwig=gamres/((pow(wm-rma,2)+pow(gamres,2))*3.1415926);
		res=ram*brwig/pm2;
		goto label30;
	label40:
		res=0;
	label30:
		ressum=ressum+res;
	}
	b=bbkg*(1.+(1.-bbkg)*xpx)+ressum*(1.-bres);
	return b;
}

//construct some parameters used in the FSI part of the code
void init_param(double *scattfactor, double *przprime, double prz,  
		double **scatterparam, double *Gamma, int num_res, 
		double nu, double qvec, double Q2, double Wx2, double Wxprime2, double Wxprime2_0, double massr, double Er, double s){
  
  for(int i=0;i<num_res;i++){
//    scattfactor[i]=sqrt((s-pow(massr-sqrt(Wxprime2),2.))*(s-pow(massr+sqrt(Wxprime2),2.)));
   przprime[i]=prz-(Er-massr)*(MASSD+nu)/qvec;
   if(Wx2>Wxprime2){ przprime[i]-=(Wx2-Wxprime2)/(2.*qvec);}
   Gamma[i]=0.;
   //for(int j=0;j<2;j++) offshellparam[i][j]=1.;
   scatterparam[i][0]=sigmainput/10.*INVHBARC*INVHBARC; //sigma_tot in GeV-2
   scatterparam[i][1]=epsinput; //alpha no units
   scatterparam[i][2]=betainput; //beta in GeV-2
   scattfactor[i]=Wxprime2_0/*+2.*qvec*przprime[i]*/;
  }
  
}


