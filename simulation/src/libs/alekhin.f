	SUBROUTINE alekhin (xb,q2,PDFS,IPAR,ICOL)
c      subroutine SFinclDIS_abkm09(xb,q2,PDFS,DPDFS,IPAR,ICOL)
c
c     This is a code accessing the grid for DIS inclusive structure 
c     functions (SFs). At the momentum transfer squared Q2>1 GeV^2 
c     the SFs are calculated in pQCD using the NNLO PDFs taken in the 
c     3-flavour scheme with account of the QCD corrections to the 
c     massless Wilson coefficients up to the NNLO. The c- and b-quark 
c     contributions are calculated in the 3-flavour scheme with account of 
c     the full NLO corrections and partial NNLO corrections due to the 
c     threshold resummation (cf. Ref.[PLB 672, 166 (2009)]).
c     The target-mass corrections by Georgi-Politzer and the 
c     dynamical high-twist (HT) contribution are taken into account. 
c     The PDFs and the HT terms for the electromagnetic SFs 
c     were fitted to the inclusive DIS non-resonant data down to 
c     W=1.8 GeV and Q2=2.5 GeV^2 (cf. Ref.[PRD 81, 014032 (2010)]). 
c     The SFs are given in the kinematics, which corresponds to the 
c     c.m.s. energy less than 2e5 GeV and momentum transferred squared 
c     bigger than 0.6 GeV^2, outside of this region they are set to 0. 
c     The full electromagnetic SFs, which include the leading-twist 
c     terms corrected for the target-mass effects and the HT terms 
c     come from a smooth interpolation between the pQCD calculations
c     at Q2>1 GeV^2 and the constraints defined by the current 
c     conservation at Q2=0 (cf. Ref.[arXiv:0710.0124]). The leading twist 
c     SFs, with and without account of the target-mass corrections, 
c     calculated in the pQCD through the whole grid kinematics are 
c     also provided, for completeness.  

c     The corrections on nuclear effects in deuterium used in the 
c     analysis of Ref.[PRD 81, 014032 (2010)] are also provided, for 
c     completeness. They are calculated with the model of 
c     Ref.[NPA 765,126 (2006)] using the PDFs stored in this grid.   

c  Input parameters:
c    XB is the parton momentum fraction 
c    Q2 is the momentum transferred. The renormalization/factorization scale
c       is set to Q2 for the massless case and sqrt(Q2+4m_h^2) for the 
c       heavy-quark contribution, where m_h is the heavy-quark mass.
c    IPAR controls output of the PDFs and \alpha_s uncertainties (see 
c         description of the output parameters below)
c    ICOL selects the collision type (the beam, final state, and current)
c         ICOL=1  --  the charged-lepton neutral current (photon exchange only)
c         ICOL=2  --  the neutrino charged current 
c         ICOL=3  --  the anti-neutrino charged current 
c         ICOL=4  --  the neutrino/anti-neutrino neutral current 

c  Output parameters:
c     The array PDFS contains fitted value of the strong coupling constant 
c     and the SFs at given x and Q:

c         PDFS(0) -- \alpha_s

c         PDFS(1) -- F_2 (p), leading twist  
c         PDFS(2) -- F_2 (p), leading twist + TMC 
c         PDFS(3) -- F_2 (p), leading twist + TMC + twist-4

c         PDFS(4) -- F_2 (n), leading twist
c         PDFS(5) -- F_2 (n), leading twist + TMC
c         PDFS(6) -- F_2 (n), leading twist + TMC + twist-4

c         PDFS(7) -- F_T (p), leading twist
c         PDFS(8) -- F_T (p), leading twist + TMC
c         PDFS(9) -- F_T (p), leading twist + TMC + twist-4

c         PDFS(10) -- F_T (n), leading twist
c         PDFS(11) -- F_T (n), leading twist + TMC
c         PDFS(12) -- F_T (n), leading twist + TMC + twist-4

c         PDFS(13) -- xF_3 (p), leading twist
c         PDFS(14) -- xF_3 (p), leading twist + TMC
c         PDFS(15) -- xF_3 (p), leading twist + TMC + twist-4

c         PDFS(16) -- xF_3 (n), leading twist
c         PDFS(17) -- xF_3 (n), leading twist + TMC
c         PDFS(18) -- xF_3 (n), leading twist + TMC + twist-4

c         PDFS(19) -- the correction on nuclear effects in deuterium for F_2
c         PDFS(20) -- the correction on nuclear effects in deuterium for F_T
c         PDFS(21) -- the correction on nuclear effects in deuterium for xF_3

c     where (p) means proton, (n) -- neutron, and F_T=2xF_1. 
c     At the moment the structure functions xF_3 are set to 0.  

c     Output array DPDFS(0:NP,NVAR) contains derivatives of \alpha_s and
c     the SFs on the parameters corresponding to the independent 
c     sources of the uncertainties, NVAR is the number of these sources, 
c     equal to 46 in the current version. The input parameter IPAR is used to 
c     optimize performance of the code. If IPAR=0, no 
c     uncertainties are returned in DPDFS; if 0<IPAR<=NVAR, only 
c     the uncertainty due to the IPAR-th source is returned; if IPAR<0
c     all uncertainties for the sources from 1 to NVAR are returned. 
c     Using derivatives stored in DPDFS one can take into account the 
c     correlations between different SFs and between SFs and \alpha_s. 
c     All derivatives are transformed to the orthonormal basis of eigenvectors 
c     of the parameters error matrix therefore variation of the SFs by 
c     the values of DPDFS is performed independently. For example, 
c     after the call with IPAR=-1 the dispersion of the i-th SF can 
c     be stored in DELPDF using the code 
c
c-----------------
c          DELPDF=0.
c          do k=1,nvar
c            DELPDF=DELPDF+dpdfs(i,k)**2
c          end do
c-----------------
c     and its random value can be stored in RPDF using the code 
c-----------------
c          RPDF=pdfs(i)          
c          do k=1,nvar
c            s=0.
c            do l=1,96
c              s=s+(2*rndm(xxx)-1)/sqrt(32.)
c            end do
c            RPDF=RPDF+s*dpdfs(i,k)
c          end do
c-----------------
c         Comments: Sergey.Alekhin@ihep.ru                      
c                                                               
c     Initial version:                                              Feb 2010    
c     - The grid was optimized to provide better interpolation 
c     accuracy at large x                                           Mar 2010  
c     - The deuteron corrections added to the grid                  May 2010

      implicit none

      character locdir*128,inputfile*128,devfile*128
      common/dir/ locdir,inputfile,devfile
      integer nxb,nq,np,nvar
      integer k,i,n,m,kx,nxbb,nplus,nminus
      parameter(nxb=99,nq=60,np=21,nvar=46)

      real*8 f(nxb,nq+1,0:np),xx(nxb)
      real*8 fsp(nxb),bs(nxb),cs(nxb),ds(nxb)
      real*8 bsp(nxb,nq+1,0:np),csp(nxb,nq+1,0:np),dsp(nxb,nq+1,0:np)
      real*8 bspd(nvar,nxb,nq+1,0:np),cspd(nvar,nxb,nq+1,0:np)
     ,      ,dspd(nvar,nxb,nq+1,0:np)
      real*8 pdfs(0:np),dpdfs(0:np,nvar)
      real*8 df(nvar,0:np,nxb,nq+1)
      real*8 x,qsq,dels,delx,x1,delx1,xlog1,xd,b,aa,ss,f0,fp,fm
      real*8 xb,q2,df0,dfp,dfm

      character suffix*5
      dimension suffix(4)
      integer icol,icols

      real*8 xmin,xmax,qsqmin,qsqmax
      integer npdf,npar1,npar2,ipar
      integer lnblnk

c I/O channel to read the data
      integer nport
      data nport/1/

      data npdf /21/

      data suffix/'lNNC ','nNNC ','nbNNC','nNNC' /

      data xmin,xmax,qsqmin,qsqmax/1d-7,1d0,0.6d0,2d5/

      save f,df,dels,delx,x1,delx1,xlog1,nxbb,xx,icols

c put in your local address of the PDF files in LOCDIR 
!       data locdir /'/home/wim/DeuteronDIS/source'/
c or in the system variable GRIDS
c      CALL GETENV( 'GRIDS', locdir ) 
c      locdir=locdir(:LNBLNK(locdir))//'pdfs/a09/'

      if (ipar.gt.nvar) print *,'Wrong call of the uncertainties' 

      if (ipar.eq.0) then 
        npar1=0
        npar2=0
      end if
      if (ipar.ge.0) then 
        npar1=ipar
        npar2=ipar
      end if
      if (ipar.lt.0) then 
        npar1=1
        npar2=nvar
      end if

      if (icols.eq.icol) goto 10

      icols=icol

      dels=(dlog(dlog(qsqmax/0.04d0))-
     +      dlog(dlog(qsqmin/0.04d0)))/dble(nq-1)

*...Read input tables
!       print *,'***** Reading SFs from tables *****'

      open(unit=nport,status='old',err=199
     , ,file=inputfile)
c     ,    ,file=locdir(:LNBLNK(locdir))
c     /     //'/a09.sfs_'//suffix(icol))

*...X GRID
      read(nport,*) (xx(i),i=1,nxb)

      do n=1,nxb
        do m=1,nq
          read(nport,*) (f(n,m,i),i=0,npdf)
        end do
      end do
      do i=0,npdf
        do m=1,nq
          do n=1,nxb
            fsp(n)=f(n,m,i)
          end do
          call spline (nxb,xx,fsp,bs,cs,ds)
          do n=1,nxb
            bsp(n,m,i)=bs(n)
            csp(n,m,i)=cs(n)
            dsp(n,m,i)=ds(n)
          end do
        end do
      end do
      close(unit=nport)

      open(unit=nport,status='old'
     , ,file=devfile)
c     ,    ,file=locdir(:LNBLNK(locdir))
c     /     //'/a09.dsfs_'//suffix(icol))
      do n=1,nxb
        do m=1,nq
          do i=0,npdf 
            read (nport,*) (df(k,i,n,m),k=1,nvar)
          end do
        end do
      end do
      close(unit=nport)

      do k=1,nvar
        do i=0,npdf
          do m=1,nq
            do n=1,nxb
              fsp(n)=df(k,i,n,m)
            end do
            call spline (nxb,xx,fsp,bs,cs,ds)
            do n=1,nxb
              bspd(k,n,m,i)=bs(n)
              cspd(k,n,m,i)=cs(n)
              dspd(k,n,m,i)=ds(n)
            end do
          end do
        end do
      end do
  10  continue

      if((q2.lt.qsqmin).or.(q2.gt.qsqmax)) then
         print 99,q2,qsqmin,qsqmax
         return
      end if
      if((xb.lt.xmin).or.(xb.gt.xmax)) then
!          print 98,xb,xmin,xmax
         return
      end if
  99  format('  AGRIDS WARNING:  Q^2 VALUE IS OUT OF RANGE   ',3g12.3)
  98  format('  AGRIDS WARNING:   X  VALUE IS OUT OF RANGE   ',3g12.3)

      x=max(xb,xmin)
      x=min(xb,xmax)
      qsq=max(q2,qsqmin)
      qsq=min(q2,qsqmax)

      do n=1,nxb
        if (x.lt.xx(n+1)) goto 300
      end do
 300  aa=x-xx(n)

      ss=dlog(dlog(qsq/0.04d0))-dlog(dlog(qsqmin/0.04d0))
      m=int(ss/dels)+1
      b=ss/dels-dble(m)+1.d0

      do i=0,npdf
        f0=f(n,m,i) + aa*bsp(n,m,i) + aa**2*csp(n,m,i) 
     +              + aa**3*dsp(n,m,i)
        fp=f(n,m+1,i) + aa*bsp(n,m+1,i) + aa**2*csp(n,m+1,i)
     +                + aa**3*dsp(n,m+1,i)
        
        if (m.ge.2) then 
          fm=f(n,m-1,i) + aa*bsp(n,m-1,i) + aa**2*csp(n,m-1,i)
     +                   +aa**3*dsp(n,m-1,i)
          pdfs(i)=(fm*b*(b-1d0)/2d0 + f0*(1d0-b**2) + fp*b*(b+1d0)/2d0)
        else 
          pdfs(i)= (f0*(1d0-b) + fp*b)              
        end if
        if (npar1.gt.0) then 
          do k=npar1,npar2
            df0=df(k,i,n,m) + aa*bspd(k,n,m,i) + aa**2*cspd(k,n,m,i) 
     +                      + aa**3*dspd(k,n,m,i)
            dfp=df(k,i,n,m+1)+aa*bspd(k,n,m+1,i)+aa**2*cspd(k,n,m+1,i)
     +                        + aa**3*dspd(k,n,m+1,i)
            if (m.ge.2) then 
              dfm=df(k,i,n,m-1)+aa*bspd(k,n,m-1,i)+aa**2*cspd(k,n,m-1,i)
     +                          + aa**3*dspd(k,n,m-1,i)
              dpdfs(i,k)=dfm*b*(b-1d0)/2d0 
     +              + df0*(1d0-b**2) +dfp*b*(b+1d0)/2d0    
            else 
              dpdfs(i,k) = df0*(1d0-b) + dfp*b            
            end if
          end do
        end if
        pdfs(i)=max(pdfs(i),0d0)
      end do

      return

 199  print *,'The grid is inavailable'

      return
      end
* ---------------------------------------------------------------------
      SUBROUTINE SPLINE(N,X,Y,B,C,D)
* ---------------------------------------------------------------------
* CALCULATE THE COEFFICIENTS B,C,D IN A CUBIC SPLINE INTERPOLATION.
* INTERPOLATION SUBROUTINES ARE TAKEN FROM
* G.E. FORSYTHE, M.A. MALCOLM AND C.B. MOLER,
* COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS (PRENTICE-HALL, 1977).
*
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION X(N), Y(N), B(N), C(N), D(N)
*
      NM1=N-1
      IF(N.LT.2) RETURN
      IF(N.LT.3) GO TO 250
      D(1)=X(2)-X(1)
      C(2)=(Y(2)-Y(1))/D(1)
      DO 210 K=2,NM1
         D(K)=X(K+1)-X(K)
         B(K)=2.0D0*(D(K-1)+D(K))
         C(K+1)=(Y(K+1)-Y(K))/D(K)
         C(K)=C(K+1)-C(K)
  210 CONTINUE
      B(1)=-D(1)
      B(N)=-D(N-1)
      C(1)=0.0D0
      C(N)=0.0D0
      IF(N.EQ.3) GO TO 215
      C(1)=C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
      C(N)=C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
      C(1)=C(1)*D(1)**2.0D0/(X(4)-X(1))
      C(N)=-C(N)*D(N-1)**2.0D0/(X(N)-X(N-3))
 215  CONTINUE
      DO 220 K=2,N
         T=D(K-1)/B(K-1)
         B(K)=B(K)-T*D(K-1)
         C(K)=C(K)-T*C(K-1)
 220  CONTINUE
      C(N)=C(N)/B(N)
      DO 230 IB=1,NM1
         K=N-IB
         C(K)=(C(K)-D(K)*C(K+1))/B(K)
 230  CONTINUE
      B(N)=(Y(N)-Y(NM1))/D(NM1)
     1     +D(NM1)*(C(NM1)+2.0D0*C(N))
      DO 240 K=1,NM1
         B(K)=(Y(K+1)-Y(K))/D(K)
     1        -D(K)*(C(K+1)+2.0D0*C(K))
         D(K)=(C(K+1)-C(K))/D(K)
         C(K)=3.0D0*C(K)
 240  CONTINUE
      C(N)=3.0D0*C(N)
      D(N)=D(N-1)
      RETURN
 250  CONTINUE
      B(1)=(Y(2)-Y(1))/(X(2)-X(1))
      C(1)=0.0D0
      D(1)=0.0D0
      B(2)=B(1)
      C(2)=0.0D0
      D(2)=0.0D0
      RETURN
      END
