C======================================================

      REAL FUNCTION DVMPX(del2,xb,Q2, Phi_g,E,heli,MESONMASS)
C
C     dvmpx(-0.493,0.187,1.61,3.92699,5.75,0.,0.135)      
C
c  dsigma/dQ2 dX dt dphi for ep-->ep pi0
C
C exc pi0 x-section
c
cinput:
cdel2=t (negative GeV^2)           NEGATIVE !!!!
cxb,Q2 x and Q^2 (GeV^2)
c Phi_g angle in the photon frame (radians)
c E energy of the electron in GeV
c heli electron helicity -1.0 or +1.0
c MESONMASS is the mass of th epi0 or eta

      IMPLICIT NONE
      REAL del2,xb,Q2, Phi_g,E,heli,MESONMASS
      real Mp, mele, pi
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)

      REAL     S_T, S_L, S_LT, S_TT,S_LTP
      REAL     XSIGMA_T, XSIGMA_L, XSIGMA_LT, XSIGMA_TT,XSIGMA_LTP
      REAL     EPS, EPSILON, FLUXW, SIGMA_TOT
      EXTERNAL EPSILON,XSIGMA_T,XSIGMA_L, XSIGMA_LT, XSIGMA_TT
      EXTERNAL XSIGMA_LTP
      REAL     FL
      LOGICAL  XCHECK_KINE
      
      INCLUDE 'pi0eta.par'

      CALL XSINIT(MESONMASS)
      
      DVMPX = 0.0
      
      IF(.NOT.XCHECK_KINE(del2,xb,Q2,E)) RETURN
       EPS  =  EPSILON(XB,q2,e)
       S_T  =  XSIGMA_T   (del2,xb,Q2,E)
       S_L  =  XSIGMA_L   (del2,xb,Q2,E)
       S_LT =  XSIGMA_LT  (del2,xb,Q2,E)
       S_TT =  XSIGMA_TT  (del2,xb,Q2,E)
       S_LTP = XSIGMA_LTP (del2,xb,Q2,E)

       DVMPX =FLUXW(xb,Q2,E)/(2.*PI)*(
     *                       S_T   + EPS*S_L          + 
     * EPS                  *S_TT  * COS(2*PHI_G)       + 
     * SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G)     +
     * HELI*SQRT(2.*EPS*(1.-EPS))*S_LTP * SIN(2*PHI_G)     
     * ) 
      if(DVMPX.lt.0.) DVMPX=0.

      RETURN
      END
c-------------------------------------------------------------------------

      SUBROUTINE XSINIT(MESONMASS)
C
      IMPLICIT NONE
      REAL MESONMASS
c      INTEGER K
c      COMMON/PIETA/K

      INCLUDE 'pi0eta.par'

      IF(MESONMASS.GT.0.140) THEN
         K=2
      ELSE
         K=1
      ENDIF

      RETURN
      END


C=============================================================================
      REAL FUNCTION XRXSB(del2,x,Q2, Phi_g,E,heli,MESONMASS)
      IMPLICIT NONE
      REAL       del2,x,Q2,Phi_g,E,heli,MESONMASS
      REAL       DVMPX,fluxw, pi
      PARAMETER (PI=3.1415926536)

      CALL XSINIT(MESONMASS)

      XRXSB=DVMPX(del2,x,Q2,Phi_g,E,heli,MESONMASS)/FLUXW(x,Q2,E)
C     *      *(2.*PI)

      RETURN
      END

C====================================================================================c
      REAL FUNCTION HT(del2,X,Q2)
      IMPLICIT NONE
      REAL          del2,x,Q2,E,KSI
      EXTERNAL      tminq
      REAL          tminq,T0,SLOPE,T,PI,PHASE
      REAL P1,P2,P3,P4,P5
      REAL MP,ALPHA
      DATA ALPHA/0.00729927/
      DATA PI/3.14159265/
      DATA MP/0.938272/

C      COMMON /PAWPAR/ P(10)
C      REAL P

      REAL W2,XNORM

      INCLUDE 'pi0eta.par'

      P1=ABS(P(1,K))
      P2=ABS(P(2,K))
      P3=ABS(P(3,K))
      P4=ABS(P(4,K))
      P5=ABS(P(5,K))

      T0=tminq(Q2,X)
      T=-DEL2

      XNORM=0.15
      SLOPE=(P2/2.)*ALOG(X)/ALOG(XNORM)
      HT=P1*EXP(-SLOPE*T)
     *     /Q2**(P3/4.)
     *     * (X**P4*(1-X)**P5)

C     PRINT *,'HT= ',HT,'  SLOPE=',slope
c      print *,'P1 ',p1,p2,p3,p4,p5
      RETURN
      END
C====================================================================================c
      REAL FUNCTION ET(del2,X,Q2)
      IMPLICIT NONE
      REAL          del2,x,Q2,E,KSI
      EXTERNAL      tminq,PHASE
      REAL          tminq,T0,SLOPE,T,PI,PHASE
      REAL P1,P2,P3,P4,P5,P6,P7,P8,P9,P10
      LOGICAL       XCHECK_KINE
      REAL MP,ALPHA
      DATA ALPHA/0.00729927/
      DATA PI/3.14159265/
      DATA MP/0.938272/
C      COMMON /PAWPAR/ P(10)
C      REAL P

     INCLUDE 'pi0eta.par'
       
      REAL W2,xnorm
      P1=ABS(P(6,K))
      P2=ABS(P(7,K))
      P3=ABS(P(8,K))
      P4=ABS(P(9,K))
      P5=ABS(P(10,K))

      T0=tminq(Q2,X)
      T=-DEL2
      XNORM=0.15
      SLOPE=(P2/2.)*ALOG(X)/ALOG(XNORM)
      ET=P1*EXP(-SLOPE*T)
     *     /Q2**(P3/4.)
     *     *(X**P4*(1-X)**P5)

c     PRINT *,'ET= ',ET,' SLOPE=',slope
c      print *,'P6 ',p1,p2,p3,p4,p5

      RETURN
      END

C====================================================================================c

      REAL FUNCTION XSIGMA_U(del2,x,Q2,E)
      IMPLICIT NONE
      REAL          del2,x,Q2,E
      EXTERNAL      XSIGMA_T,XSIGMA_L,EPSILON
      REAL          tminq,T0,SLOPE,T,XSIGMA_T,
     *               XSIGMA_L,EPSILON
      LOGICAL       XCHECK_KINE

       REAL AM(2)
       DATA AM/0.134976,0.547300/
       INTEGER K
       COMMON/PIETA/K
c      INCLUDE 'pi0eta.par'
      
      IF(XCHECK_KINE(del2,x,Q2,E)) THEN
        XSIGMA_U=XSIGMA_T(DEL2,X,Q2,E)+
     *  EPSILON(x,Q2,E)*XSIGMA_L(DEL2,X,Q2,E)
      ELSE
         XSIGMA_U=0.0
      ENDIF
      RETURN
      END

C====================================================================================c

      REAL FUNCTION XSIGMA_T(del2,x,Q2,E)
      IMPLICIT NONE
      REAL          del2,x,Q2,E,ALPHA,KSI
      EXTERNAL      tminq,PHASE,HT,ET
      REAL          tminq,T0,SLOPE,T,PI,PHASE
      REAL          HT,ET,MP,HC2
      LOGICAL       XCHECK_KINE
      DATA ALPHA/0.00729927/
      DATA PI/3.14159265/
      DATA MP/0.938272/
      DATA HC2/389379.36/
      
      IF(XCHECK_KINE(del2,x,Q2,E)) THEN
        T0=tminq(Q2,X)
        T=-DEL2
        KSI=X/(2-X)*(1+MP**2/Q2)
        XSIGMA_T=4.*PI*ALPHA/2./PHASE(Q2,X)*(
     *          (1-KSI**2)*    HT(del2,X,Q2)**2
     *       +  (T-T0)/8/MP**2*ET(del2,X,Q2)**2
     *   )
      ELSE
         XSIGMA_T=0.0
      ENDIF
      XSIGMA_T=XSIGMA_T*HC2

c     print *,'XSIGMA_T=', XSIGMA_T

      RETURN
      END
      

C====================================================================================c

       REAL FUNCTION XSIGMA_TT(del2,x,Q2,E)
       IMPLICIT NONE
       REAL          del2,x,Q2,E
       EXTERNAL      tminq
       REAL          tminq,T0,SLOPE,T,HC2
       REAL ALPHA,PI,MP,PHASE,ET
       LOGICAL XCHECK_KINE
       DATA ALPHA/0.00729927/
       DATA PI/3.14159265/
       DATA MP/0.938272/
       DATA HC2/389379.36/

       IF(XCHECK_KINE(del2,x,Q2,E)) THEN
         T0=tminq(Q2,X)
         T=-DEL2
        XSIGMA_TT=-4.*PI*ALPHA/2./PHASE(Q2,X)
     *         *  (T-T0)/8/MP**2*ET(del2,X,Q2)**2
      ELSE
         XSIGMA_TT=0.0
      ENDIF

      XSIGMA_TT=XSIGMA_TT*HC2

c     print *,' t=',-del2,' XSIGMA_TT=', xsigma_tt
      
       RETURN
       END

C====================================================================================c

      REAL FUNCTION XSIGMA_L(del2,x,Q2,E)
      IMPLICIT NONE
      REAL          del2,x,Q2,E
      EXTERNAL      tminq
      REAL          tminq,T0,SLOPE,T
      LOGICAL XCHECK_KINE

c      INCLUDE 'pi0eta.par'

      XSIGMA_L=0.0

       RETURN
       END
C====================================================================================c

       REAL FUNCTION XSIGMA_LT(del2,x,Q2,E)
       IMPLICIT NONE
       REAL          del2,x,Q2,E
       EXTERNAL      tminq
       REAL          tminq,T0,SLOPE,T
       LOGICAL XCHECK_KINE

       XSIGMA_LT=0.0

       RETURN
       END


C====================================================================================c

       REAL FUNCTION XSIGMA_LTP(del2,x,Q2,E)
       IMPLICIT NONE
       REAL del2,x,Q2,E
       REAL EPS, EPSILON,A2_PHI,FLUXW,SIGMA_TOT
       EXTERNAL EPSILON,FLUXW


       XSIGMA_LTP=0.

       RETURN
       END

C====================================================================================c

      LOGICAL FUNCTION XCHECK_KINE(del2,xb,Q2,E)
      IMPLICIT NONE
      REAL del2,xb,Q2,E
      real Mp, mele, pi,mpi0
      real fluxw, rxs
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)
           

      real nu,W2,W,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
      real  xmin1,xmax1
      real y, e1, epsilon,tminq

       REAL AM(2)
       DATA AM/0.134976,0.547300/
       INTEGER K
       COMMON/PIETA/K

      MPI0=AM(K)
      XCHECK_KINE=.FALSE.
      xmin1 = Q2/(2.0*Mp*E)
      xmax1 = 1.0
      nu  = Q2/(2D0*Mp*xb)
      W2  = Mp**2 + 2.0*Mp*nu - Q2
C      print *,'W2=',W2
      IF(W2.LT.(Mp+Mpi0)**2)                  RETURN
      W   = sqrt(W2)
      qmod = sqrt(nu**2 + Q2)

      E1cm = Mp*(Mp + nu)/W
      P1cm = Mp*qmod/W
      E2cm = (W2 + Mp**2-Mpi0**2)/(2.*W)
C      print *,'E2cm=',E2cm,Mp
      IF(E2cm.LE.Mp)                             RETURN
      P2cm = SQRT(E2CM**2 - Mp**2)
      del2max = 2.0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
      del2min = 2.0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)
C      print *,'xB=',xb,xmin1,xmax1
      IF( xb.le.xmin1 .or. xb.gt.xmax1 )         return !    x  out of range
C      print *,'t=',del2,del2min,del2max,tminq(Q2,xb)
      IF( del2.ge.del2min .or. del2.le.del2max ) return   ! delta out of range

      y=Q2/(2*Mp*xb*E)
      e1=(y*xb*Mp)**2/Q2
      EPSILON=(1.0-y-e1)/(1-y+y**2/2+e1)
C      print *,'Epsilon',epsilon
      IF(EPSILON.LT.0. .OR.EPSILON .GT.1.)         RETURN

      XCHECK_KINE=.TRUE.

      RETURN
      END

C====================================================================================c

      real function fluxw(x,Q2,E)
      implicit none
      real x,Q2,E,y,eps,e1
      real alpha,Mp,PI
      parameter (alpha=1.0/137.036,Mp=0.93827231,PI=3.14151926)

c
      y=Q2/(2*Mp*x*E)
      e1=(y*x*Mp)**2/Q2
      eps=(1.0-y-e1)/(1-y+y**2/2+e1)
      fluxw=alpha/(2*PI)*y*y/(1-eps)*(1-x)/x/Q2
      return
      end

C====================================================================================c

      real function EPSILON(x,Q2,E)
      implicit none
      real x,Q2,E,y,eps,e1
      real alpha,Mp,PI
      parameter (alpha=1.0/137.036,Mp=0.93827231,PI=3.14151926)
c
      y=Q2/(2*Mp*x*E)
      e1=(y*x*Mp)**2/Q2
      EPSILON=(1.0-y-e1)/(1-y+y**2/2+e1)
      return
      end


C======================================================================================c

      REAL FUNCTION tminq(Q2,X)
      IMPLICIT NONE
      REAL Q2,X
      REAL W2,W,S,E1CM,E3CM,P1CM,P3CM,TMAX
      INTEGER I
      real alpha,Mp,PI,Me
      REAL MPI0

      INCLUDE 'pi0eta.par'

      MPI0=AM(K)
      tminq=0.

      Mp=0.9382723
      Me=0.0051
      PI=3.14151926

      IF(X.LE.0. .OR. X.GE.1.)        RETURN
      W2 = Q2*(1./X-1.)+Mp**2
      W=SQRT(W2)
      IF(W.LT.Mp+Mpi0)                RETURN

      E1CM=(W2+Q2+Mp**2)/(2*W)
      P1CM=SQRT(E1CM**2-MP**2)

      E3CM=(W2-MPI0**2+Mp**2)/(2*W)
      P3CM=SQRT(E3CM**2-MP**2)
      
      TMINQ=-((Q2+MPI0**2)**2/4./W2-(P1CM-P3CM)**2)

      RETURN
      END

C====================================================================================c
      REAL FUNCTION PHASE(Q2,X)
      IMPLICIT NONE
      REAL Q2,X
      REAL W,W2,W4,Q4,MP,MP2,MP4,PI,LAMBDA
      DATA MP/0.938272/
      DATA PI/3.14159265/
      PHASE=0.0
      
      IF(X.GE.1.) RETURN
      W2=Q2*(1/X-1.)+MP**2
      W4=W2*W2
      MP2=MP*MP
      MP4=MP2*MP2

      LAMBDA=W4+Q4+MP4+2*W2*Q2-2*W2*MP2+2*Q2*MP2
      IF(LAMBDA.LE.0) RETURN

      PHASE=16.*PI*(W2-MP2)*SQRT(LAMBDA)
      RETURN
      END
      
      
