c 	Multipole and Helicity Amplitudes

      integer wave_L		! max pi-N ang. mom. in table
      integer mwave_L		! max pi-N ang. mom. calculated
      integer max_L		! max hel wave summed
      integer method_helicity 	! =1, for Tiator's
	                      	! =2, for Burkert's

      parameter(max_L = 5)
      parameter(wave_L = 5)
	  
      complex sp(0:max_L+1), sm(0:max_L+1)
      complex ep(0:max_L+1), em(0:max_L+1)
      complex mp(0:max_L+1), mm(0:max_L+1)
      complex ap(0:max_L+1), am(0:max_L+1)
      complex bp(0:max_L+1), bm(0:max_L+1)
      complex cp(0:max_L+1), cm(0:max_L+1)

      common /mp_amplitudes/sp, sm, ep, em, mp, mm, 
     1                      ap, am, bp, bm, cp, cm
      common /method/method_helicity,mwave_L

      complex hh1,hh2,hh3,hh4,hh5,hh6 ! helicity amplitudes
      complex ff1,ff2,ff3,ff4,ff5,ff6 ! CGLN amplitudes

      common /hel_amplitudes/hh1,hh2,hh3,hh4,hh5,hh6
      common /cgln_amplitudes/ff1,ff2,ff3,ff4,ff5,ff6
      
      real pol(0:max_L+2,2)

      common /legend/ pol

C Parameters and options

	  character*132 data_file

	  common /input_file/data_file      

C*******************************************************************************
CEV The following is for the response functions.
C*******************************************************************************
	  
	  real sigma_0, sigma_l, sigma_t, sigma_tt, sigma_lt, sigma_ltp, asym_ltp
	  real sigma_0_5fold
	  
	  common /reponse_fcns/asym_ltp, sigma_0, sigma_l, sigma_t, 
     1                     sigma_tt, sigma_lt, sigma_ltp, sigma_0_5fold

C*******************************************************************************
C The following is for the interpolation routines.
C*******************************************************************************



C Parameters

      INTEGER MAX_MP
      INTEGER NVAR1, NVAR2 

      PARAMETER( MAX_MP = 62 )
C for OOPS:   NVAR1 = 16, NVAR2 = 26
C for SatoLee:   NVAR1 = 77, NVAR2 = 26
C for NSTAR:  NVAR1 = 77, NVAR2 = 61
c for NSTAR 2000 format:  NVAR1 = 81, NVAR2 = 91

C 	  PARAMETER( NVAR1  = 77,   NVAR2 = 26 )	! Q2,W 
c      PARAMETER( NVAR1  = 77 , NVAR2 = 61 )	! Q2,W 
C      PARAMETER( NVAR1  = 16 , NVAR2 = 26 )	! Q2,W 
c      PARAMETER( NVAR1  = 81 , NVAR2 = 91 )	! Q2,W 
      PARAMETER( NVAR1  = 101 , NVAR2 = 93 )	! Q2,W maid 07

C Vectors of Kinematic Variables
C
      REAL VAR1(NVAR1),VAR2(NVAR2)
      REAL VAR(NVAR1+NVAR2)
      integer NA(2)
      REAL var1_min, var1_max
      REAL var2_min, var2_max
C
	  integer method_spline
C
      data NA /nvar1, nvar2/
C
      COMMON /var_maxmin/var1_min,var1_max,var2_min,var2_max,
     1                   method_spline

C Labeling Scheme for multipoles:
C
C                     SF1,2   = Real, Imag S0+
C                     SF3,4   = Real, Imag S1+
C                     SF5,6   = Real, Imag S2+
C                     SF7,8   = Real, Imag S3+
C                     SF9,10  = Real, Imag S4+
C                     SF11,12 = Real, Imag S5+
C
C                     SF13,14 = Real, Imag S1-
C                     SF15,16 = Real, Imag S2-
C                     SF17,18 = Real, Imag S3-
C                     SF19,20 = Real, Imag S4-
C                     SF21,22 = Real, Imag S5-
C
C                     SF23,24 = Real, Imag E0+
C                     SF25,26 = Real, Imag E1+
C                     SF27,28 = Real, Imag E2+
C                     SF29,30 = Real, Imag E3+
C                     SF31,32 = Real, Imag E4+
C                     SF33,34 = Real, Imag E5+
C
C                     SF35,36 = Real, Imag E2-
C                     SF37,38 = Real, Imag E3-
C                     SF39,40 = Real, Imag E4-
C                     SF41,42 = Real, Imag E5-
C
C                     SF43,44 = Real, Imag M1+
C                     SF45,46 = Real, Imag M2+
C                     SF47,48 = Real, Imag M3+
C                     SF49,50 = Real, Imag M4+
C                     SF51,52 = Real, Imag M5+
C
C                     SF53,54 = Real, Imag M1-
C                     SF55,56 = Real, Imag M2-
C                     SF57,58 = Real, Imag M3-
C                     SF59,60 = Real, Imag M4-
C                     SF61,62 = Real, Imag M5-

CEV Multipoles:

      REAL    SF1(NVAR1,NVAR2),  SF2(NVAR1,NVAR2),
     &        SF3(NVAR1,NVAR2),  SF4(NVAR1,NVAR2),
     &        SF5(NVAR1,NVAR2),  SF6(NVAR1,NVAR2),
     &        SF7(NVAR1,NVAR2),  SF8(NVAR1,NVAR2),
     &        SF9(NVAR1,NVAR2), SF10(NVAR1,NVAR2),
     &       SF11(NVAR1,NVAR2), SF12(NVAR1,NVAR2),
     &       SF13(NVAR1,NVAR2), SF14(NVAR1,NVAR2),
     &       SF15(NVAR1,NVAR2), SF16(NVAR1,NVAR2),
     &       SF17(NVAR1,NVAR2), SF18(NVAR1,NVAR2),
     &       SF19(NVAR1,NVAR2), SF20(NVAR1,NVAR2),
     &       SF21(NVAR1,NVAR2), SF22(NVAR1,NVAR2),
     &       SF23(NVAR1,NVAR2), SF24(NVAR1,NVAR2),
     &       SF25(NVAR1,NVAR2), SF26(NVAR1,NVAR2),
     &       SF27(NVAR1,NVAR2), SF28(NVAR1,NVAR2),
     &       SF29(NVAR1,NVAR2), SF30(NVAR1,NVAR2),
     &       SF31(NVAR1,NVAR2), SF32(NVAR1,NVAR2),
     &       SF33(NVAR1,NVAR2), SF34(NVAR1,NVAR2),
     &       SF35(NVAR1,NVAR2), SF36(NVAR1,NVAR2),
     &       SF37(NVAR1,NVAR2), SF38(NVAR1,NVAR2)
      real   SF39(NVAR1,NVAR2), SF40(NVAR1,NVAR2),
     &       SF41(NVAR1,NVAR2), SF42(NVAR1,NVAR2),
     &       SF43(NVAR1,NVAR2), SF44(NVAR1,NVAR2),
     &       SF45(NVAR1,NVAR2), SF46(NVAR1,NVAR2),
     &       SF47(NVAR1,NVAR2), SF48(NVAR1,NVAR2),
     &       SF49(NVAR1,NVAR2), SF50(NVAR1,NVAR2),
     &       SF51(NVAR1,NVAR2), SF52(NVAR1,NVAR2),
     &       SF53(NVAR1,NVAR2), SF54(NVAR1,NVAR2),
     &       SF55(NVAR1,NVAR2), SF56(NVAR1,NVAR2),
     &       SF57(NVAR1,NVAR2), SF58(NVAR1,NVAR2),
     &       SF59(NVAR1,NVAR2), SF60(NVAR1,NVAR2),
     &       SF61(NVAR1,NVAR2), SF62(NVAR1,NVAR2)

CEV Second derivatives for cubic spline:

      REAL    D2SF1(NVAR1,NVAR2),  D2SF2(NVAR1,NVAR2),
     &        D2SF3(NVAR1,NVAR2),  D2SF4(NVAR1,NVAR2),
     &        D2SF5(NVAR1,NVAR2),  D2SF6(NVAR1,NVAR2),
     &        D2SF7(NVAR1,NVAR2),  D2SF8(NVAR1,NVAR2),
     &        D2SF9(NVAR1,NVAR2), D2SF10(NVAR1,NVAR2),
     &       D2SF11(NVAR1,NVAR2), D2SF12(NVAR1,NVAR2),
     &       D2SF13(NVAR1,NVAR2), D2SF14(NVAR1,NVAR2),
     &       D2SF15(NVAR1,NVAR2), D2SF16(NVAR1,NVAR2),
     &       D2SF17(NVAR1,NVAR2), D2SF18(NVAR1,NVAR2),
     &       D2SF19(NVAR1,NVAR2), D2SF20(NVAR1,NVAR2),
     &       D2SF21(NVAR1,NVAR2), D2SF22(NVAR1,NVAR2),
     &       D2SF23(NVAR1,NVAR2), D2SF24(NVAR1,NVAR2),
     &       D2SF25(NVAR1,NVAR2), D2SF26(NVAR1,NVAR2),
     &       D2SF27(NVAR1,NVAR2), D2SF28(NVAR1,NVAR2),
     &       D2SF29(NVAR1,NVAR2), D2SF30(NVAR1,NVAR2),
     &       D2SF31(NVAR1,NVAR2), D2SF32(NVAR1,NVAR2),
     &       D2SF33(NVAR1,NVAR2), D2SF34(NVAR1,NVAR2),
     &       D2SF35(NVAR1,NVAR2), D2SF36(NVAR1,NVAR2),
     &       D2SF37(NVAR1,NVAR2), D2SF38(NVAR1,NVAR2)
      REAL   D2SF39(NVAR1,NVAR2), D2SF40(NVAR1,NVAR2),
     &       D2SF41(NVAR1,NVAR2), D2SF42(NVAR1,NVAR2),
     &       D2SF43(NVAR1,NVAR2), D2SF44(NVAR1,NVAR2),
     &       D2SF45(NVAR1,NVAR2), D2SF46(NVAR1,NVAR2),
     &       D2SF47(NVAR1,NVAR2), D2SF48(NVAR1,NVAR2),
     &       D2SF49(NVAR1,NVAR2), D2SF50(NVAR1,NVAR2),
     &       D2SF51(NVAR1,NVAR2), D2SF52(NVAR1,NVAR2),
     &       D2SF53(NVAR1,NVAR2), D2SF54(NVAR1,NVAR2),
     &       D2SF55(NVAR1,NVAR2), D2SF56(NVAR1,NVAR2),
     &       D2SF57(NVAR1,NVAR2), D2SF58(NVAR1,NVAR2),
     &       D2SF59(NVAR1,NVAR2), D2SF60(NVAR1,NVAR2),
     &       D2SF61(NVAR1,NVAR2), D2SF62(NVAR1,NVAR2)

      COMMON/THEORY_MESH/VAR1,VAR2, VAR,
     &                    SF1, SF2, SF3, SF4, SF5, SF6, SF7, SF8,
     &                    SF9,SF10,SF11,SF12,SF13,SF14,SF15,SF16,
     &                   SF17,SF18,SF19,SF20,SF21,SF22,SF23,SF24,
     &                   SF25,SF26,SF27,SF28,SF29,SF30,SF31,SF32,
     &                   SF33,SF34,SF35,SF36,SF37,SF38,SF39,SF40,
     &                   SF41,SF42,SF43,SF44,SF45,SF46,SF47,SF48,
     &                   SF49,SF50,SF51,SF52,SF53,SF54,SF55,SF56,
     &                   SF57,SF58,SF59,SF60,SF61,SF62,
     &                    D2SF1, D2SF2, D2SF3, D2SF4, D2SF5, D2SF6,
     &                    D2SF7, D2SF8, D2SF9,D2SF10,D2SF11,D2SF12,
     &                   D2SF13,D2SF14,D2SF15,D2SF16,D2SF17,D2SF18,
     &                   D2SF19,D2SF20,D2SF21,D2SF22,D2SF23,D2SF24,
     &                   D2SF25,D2SF26,D2SF27,D2SF28,D2SF29,D2SF30,
     &                   D2SF31,D2SF32,D2SF33,D2SF34,D2SF35,D2SF36,
     &                   D2SF37,D2SF38,D2SF39,D2SF40,D2SF41,D2SF42,
     &                   D2SF43,D2SF44,D2SF45,D2SF46,D2SF47,D2SF48,
     &                   D2SF49,D2SF50,D2SF51,D2SF52,D2SF53,D2SF54,
     &                   D2SF55,D2SF56,D2SF57,D2SF58,D2SF59,D2SF60,
     &                   D2SF61,D2SF62

















