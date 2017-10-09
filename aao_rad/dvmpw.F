c*******************Valery Model****************************************
      subroutine DVMPW(COST,W,Q2, Phi_g,E,heli,MESONMASS,
     1                 sigma0,S_T,S_L,S_TT,S_LT,S_LTP)
C
c  dsigma/dQ2 dphi dCosTheta* dW for ep-->ep pi0
C
C  exc pi0 x-section
c
cinput:
c COST = Cos(Theta*), Theta* is the angle (pi0-gamma* or  p-p')in the CM system
c W is W
c xb,Q2 x and Q^2 (GeV^2)
c Phi_g angle in the photon frame (radians)
c E energy of the electron in GeV
c heli electron helicity -1.0 or +1.0
c
cdel2=t (negative GeV^2)           NEGATIVE !!!!
C
C MESONMASS is the mass of the pi0 or eta.
C The actual masses that will be used for the calculations are in pi0eta.par file


      IMPLICIT NONE
      REAL    COST,W,del2,xb,Q2, Phi_g,E,heli, mesonmass
      REAL    Mp, mele, pi
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)

      REAL     sigma0,S_T, S_L, S_LT, S_TT,S_LTP
      REAL     SIGMA_T, SIGMA_L, SIGMA_LT, SIGMA_TT,SIGMA_LTP
      REAL     EPS, EPSILON, FLUXW, SIGMA_TOT
      EXTERNAL EPSILON,SIGMA_T,SIGMA_L, SIGMA_LT, SIGMA_TT
      EXTERNAL SIGMA_LTP
      REAL     JACG,JACR
      LOGICAL  CHECK_KINE
c      REAL DVMPX
      INCLUDE 'pi0eta.par'
c
      
      CALL XSINIT(MESONMASS)

c     
      S_T = 0.0
      S_L = 0.0
      S_LT = 0.0
      S_TT = 0.0
      S_LTP = 0.0
      IF(.NOT.CHECK_KINE(COST,W,Q2,E,del2,xb,JACG,JACR)) RETURN 
       EPS=EPSILON(XB,q2,e)
       S_T  =  0.001/(2.*PI)*SIGMA_T   (COST,W,Q2,E)
       S_L  =  0.001/(2.*PI)*SIGMA_L   (COST,W,Q2,E)
       S_LT =  0.001/(2.*PI)*2.0*SIGMA_LT  (COST,W,Q2,E)
       S_TT =  0.001/(2.*PI)*SIGMA_TT  (COST,W,Q2,E)
       S_LTP = 0.001/(2.*PI)*2.0*SIGMA_LTP (COST,W,Q2,E)

c
c       sigma0 =0.000001*FLUXW(xb,Q2,E)/(2.*PI)*JACG*(
       sigma0 =1/(2.*PI)*(                                                                                                                             

     *                       S_T   + EPS*S_L          + 
     * EPS                  *S_TT  * COS(2*PHI_G)       + 
     * SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G)     +
     * HELI*SQRT(2.*EPS*(1.-EPS))*S_LTP * SIN(2*PHI_G)     
     * )

       IF (sigma0.LT.0.0) then
         sigma0=0.0
         S_T = 0.0 
         S_L = 0.0 
         S_LT = 0.0 
         S_TT = 0.0 
         S_LTP = 0.0 
       ENDIF
c       PRINT *,DVMPW,COST,W,Q2,del2,xb,phi_g,JACG,JACR
c      PRINT *,DVMPX(del2,xb,Q2, Phi_g,E,heli,MESONMASS)*JACG*JACR
c      PRINT *,DVMPX(del2,xb,Q2, Phi_g,E,heli,MESONMASS)
c      PRINT *,del2,xb,Q2, Phi_g,E,heli,MESONMASS

c       PRINT *,DVMPW,COST,W,Q2,phi_g,eps
c       print *,'DVMP',sigma0, S_T,S_L,S_TT,S_LT
c       print *,'S_TT  ',  S_TT, EPS     *S_TT  * COS(2*PHI_G) 
c       print *, 'S_LT ',S_LT,SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G) 

      RETURN
      END

C====================================================================================c

      REAL FUNCTION SIGMA_T(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      EXTERNAL      tminq
      REAL          tminq,T0,SLOPE,T,JACG,JACR
      LOGICAL       CHECK_KINE
      EXTERNAL      XSIGMA_T
      REAL          XSIGMA_T
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_T=JACR*XSIGMA_T(del2,x,Q2,E)
      ELSE
        SIGMA_T=0.0
      ENDIF
      RETURN
      END
C=============================================================      
      REAL FUNCTION SIGMA_L(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      EXTERNAL      tminq
      REAL          tminq,T0,SLOPE,T,JACG,JACR
      LOGICAL CHECK_KINE
      EXTERNAL      XSIGMA_L
      REAL          XSIGMA_L
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_L=JACR*XSIGMA_L(del2,x,Q2,E)
      ELSE
        SIGMA_L=0.0
      ENDIF
      RETURN
      END
      
C=============================================================      
      REAL FUNCTION SIGMA_TT(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      EXTERNAL      tminq
      REAL          tminq,T0,SLOPE,T,JACG,JACR
      LOGICAL CHECK_KINE
      EXTERNAL      XSIGMA_TT
      REAL          XSIGMA_TT
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_TT=JACR*XSIGMA_TT(del2,x,Q2,E)
      ELSE
        SIGMA_TT=0.0
      ENDIF
      RETURN
      END
      
C=============================================================      
      REAL FUNCTION SIGMA_LT(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      EXTERNAL      tminq
      REAL          tminq,T0,SLOPE,T,JACG,JACR
      LOGICAL CHECK_KINE
      EXTERNAL      XSIGMA_LT
      REAL          XSIGMA_LT
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_LT=JACR*XSIGMA_LT(del2,x,Q2,E)
      ELSE
        SIGMA_LT=0.0
      ENDIF
      RETURN
      END
      
C====================================================================================c

       REAL FUNCTION SIGMA_LTP(COST,W,Q2,E)
       IMPLICIT NONE
       REAL COST,W,del2,x,Q2,E
       REAL EPS, EPSILON,A2_PHI,FLUXW,SIGMA_TOT
       EXTERNAL EPSILON,FLUXW

       SIGMA_LTP=0.

       RETURN
       END

C====================================================================================c

      LOGICAL FUNCTION CHECK_KINE(COST,W,Q2,E,del2,xb,JACG,JACR)
      IMPLICIT NONE
      REAL COST,W,del2,xb,Q2,E,JACG,JACR
      real Mp, mele, pi,mpi0
      real fluxw, rxs
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)
           

      real nu,W2,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
      real  xmin1,xmax1
      real y, e1, epsilon
      INCLUDE 'pi0eta.par'

      MPI0=AM(K)

      CHECK_KINE=.FALSE.
      W2=W*W
      XB=Q2/(W2+Q2-Mp*Mp)
      nu  = (W2+Q2-Mp*Mp)/(2*Mp)
      qmod = sqrt(nu**2 + Q2)
      IF(W         .LT. Mp+Mpi0)              RETURN
      IF(W         .GT. -Q2+2*Mp*E+Mp*Mp)     RETURN
      IF(Q2/(4*E*(E-NU)) .GE. 1.0)            RETURN
      IF(XB        .LT. Q2/(2*Mp*E))          RETURN
      IF(ABS(COST) .GT. 1.0)                  RETURN
      xmin1 = Q2/(2.0*Mp*E)
      xmax1 = 1.0

      E1cm = Mp*(Mp + nu)/W
      E2cm = (W2 + Mp**2-Mpi0**2)/(2D0*W)
      IF(E1cm.LE.Mp .OR. E2cm.LE. Mp) RETURN
      P1cm = Mp*qmod/W
      P2cm = SQRT(E2CM**2 - Mp**2)
      del2max = 2.0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
      del2min = 2.0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)
      del2=del2min-2.*P1cm*P2cm*(1-COST)
      IF( xb.le.xmin1 .or. xb.gt.xmax1 )         return   !    x  out of range
      IF( del2.ge.del2min .or. del2.le.del2max ) return   ! delta out of range

      y=Q2/(2*Mp*xb*E)
      e1=(y*xb*Mp)**2/Q2
      EPSILON=(1.0-y-e1)/(1-y+y**2/2+e1)

      IF(EPSILON.LT.0. .OR.EPSILON .GT.1.)         RETURN

C      jacobian
C      JAC=4*P1cm*P2CM*W*XB/(W2+Q2-Mp*Mp)
      JACG=2*W*XB/(W2+Q2-Mp*Mp)
      JACR=2*P1cm*P2CM

      CHECK_KINE=.TRUE.

      RETURN
      END

C=======================================================================================C

