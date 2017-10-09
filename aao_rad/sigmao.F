*************************************************************************                                                            
         REAL FUNCTION SIGMAO(EPW,EPQ2,EPEPS,EPCOS,EPPHI,EPIREA,
     *   POL_ELEC,POL_TARG,POL_TARG_THETA,POL_TARG_PHI,sig0
     *   ,sigu,sigt,sigl,sigi,fkt)

         implicit none
         real epw,epq2,epeps,epcos,EPPHI,pol_elec,pol_targ
         real pol_targ_theta,pol_targ_phi,sig0,sigu,sigt,sigl,sigi
         real fkt
         real sindd,cosdd,pol_x,pol_y,pol_z,sigpt,siget
         real xila,sige,theta_gam,enue,theta_elec,e_elec_scatt
         real e_elec_prim
         real c2,c3
         real q2,sinx,cosx,sinx2,cosx2
         real eph1b,eph2b,eph3b,eph4b,eph5b,eph6b
         real xmpot
         real sq2,pi,xmprot,xmpip,xmpi0,xmpi,dum
         real epsin,epsin2,epcos2,pener,pimp,gener,gimp,gamk
         real p1s,p2s,p3s,p4s,p5s,p6s,p7s,p2s2,p3s2,p4s2,p5s2
         real p6s2,p7s2
         real x,x0,f
         real bpipy0,bpipy1,bpipy2,bpipy3
         integer irea
c	The folowing were not declared in the original AO (RCM)
c	I think they are supposed to be complex
         complex xa4m,xb4m,xc4m,xa4p,xb4p,xc4p,xa5m,xb5m,xc5m,xa5p
         complex xb5p,xc5p,xa6m,xb6m,xc6m

c      CHARACTER*3 GON
c      LOGICAL TEST
      COMMON /IP/IT,IB,NF,IBORN,CUT,IP11
c      COMMON /GOCOM/TEST,GON
      real cut
      integer it,ib,nf,iborn,ip11

         COMMON /BACKG/      
     * BGA0P0,BGA1P0,BGA2P0,BGA3P0,  
     * BGA1M0,BGA2M0,BGA3M0,
     * BGB1P0,BGB2P0,BGB3P0, 
     * BGB2M0,BGB3M0,  
     * BGC0P0,BGC1P0,BGC2P0,BGC3P0,
     * BGC1M0,BGC2M0,BGC3M0,
     * BGA0PP,BGA1PP,BGA2PP,BGA3PP,
     * BGA1MP,BGA2MP,BGA3MP,
     * BGB1PP,BGB2PP,BGB3PP,
     * BGB2MP,BGB3MP,
     * BGC0PP,BGC1PP,BGC2PP,BGC3PP,
     * BGC1MP,BGC2MP,BGC3MP
       COMMON/RCONST/      
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710, 
     * WP13_1720, WD13_1520, WD13_1700, 
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600, 
     * WS31_1620, WS31_1900, WP31_1910,   
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710, 
     * LP13_1720, LD13_1520, LD13_1700, 
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600, 
     * LS31_1620, LS31_1900, LP31_1910,  
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710, 
     * PIBP13_1720, PIBD13_1520, PIBD13_1700, 
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600, 
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,  
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975
         COMMON/RAMP/
     * RAS11_1535,             RCS11_1535, 
     * RAS11_1650,             RCS11_1650, 
     * RAP11_1440,             RCP11_1440, 
     * RAP11_1710,             RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,         
     * RAP13_1910, RBP13_1910, RCP13_1910,         
     * RAD13_1520, RBD13_1520, RCD13_1520, 
     * RAD13_1700, RBD13_1700, RCD13_1700, 
     * RAD15_1675, RBD15_1675, RCD15_1675, 
     * RAF15_1680, RBF15_1680, RCF15_1680, 
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,             RCS31_1620, 
     * RAS31_1900,             RCS31_1900, 
     * RAP31_1910,             RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232, 
     * RAP33_1920, RBP33_1920, RCP33_1920, 
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420,
     * RAP33_1600, RBP33_1600, RCP33_1600,
     * RAF17_1990, RBF17_1990, RCF17_1990,
     * RAF15_2000, RBF15_2000, RCF15_2000,
     * RAP11_2100, RBP11_2100, RCP11_2100,
     * RAF35_2000, RBF35_2000, RCF35_2000,
     * RAP13_1870, RBP13_1870, RCP13_1870,
     * RAP31_1925, RBP31_1925, RCP31_1925,
     * RAP13_1980, RBP13_1980, RCP13_1980,
     * RAF15_1955, RBF15_1955, RCF15_1955,
     * RAP13_1955, RBP13_1955, RCP13_1955,
     * RAP33_1975, RBP33_1975, RCP33_1975

         real
     * BGA0P0,BGA1P0,BGA2P0,BGA3P0,  
     * BGA1M0,BGA2M0,BGA3M0,
     * BGB1P0,BGB2P0,BGB3P0, 
     * BGB2M0,BGB3M0,  
     * BGC0P0,BGC1P0,BGC2P0,BGC3P0,
     * BGC1M0,BGC2M0,BGC3M0,
     * BGA0PP,BGA1PP,BGA2PP,BGA3PP,
     * BGA1MP,BGA2MP,BGA3MP,
     * BGB1PP,BGB2PP,BGB3PP,
     * BGB2MP,BGB3MP,
     * BGC0PP,BGC1PP,BGC2PP,BGC3PP,
     * BGC1MP,BGC2MP,BGC3MP

       real
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710, 
     * WP13_1720, WD13_1520, WD13_1700, 
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600, 
     * WS31_1620, WS31_1900, WP31_1910,   
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710, 
     * LP13_1720, LD13_1520, LD13_1700, 
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600, 
     * LS31_1620, LS31_1900, LP31_1910,  
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710, 
     * PIBP13_1720, PIBD13_1520, PIBD13_1700, 
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600, 
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,  
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975

       real
     * RAS11_1535,             RCS11_1535, 
     * RAS11_1650,             RCS11_1650, 
     * RAP11_1440,             RCP11_1440, 
     * RAP11_1710,             RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,         
     * RAP13_1910, RBP13_1910, RCP13_1910,         
     * RAD13_1520, RBD13_1520, RCD13_1520, 
     * RAD13_1700, RBD13_1700, RCD13_1700, 
     * RAD15_1675, RBD15_1675, RCD15_1675, 
     * RAF15_1680, RBF15_1680, RCF15_1680, 
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,             RCS31_1620, 
     * RAS31_1900,             RCS31_1900, 
     * RAP31_1910,             RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232, 
     * RAP33_1920, RBP33_1920, RCP33_1920, 
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420,
     * RAP33_1600, RBP33_1600, RCP33_1600,
     * RAF17_1990, RBF17_1990, RCF17_1990,
     * RAF15_2000, RBF15_2000, RCF15_2000,
     * RAP11_2100, RBP11_2100, RCP11_2100,
     * RAF35_2000, RBF35_2000, RCF35_2000,
     * RAP13_1870, RBP13_1870, RCP13_1870,
     * RAP31_1925, RBP31_1925, RCP31_1925,
     * RAP13_1980, RBP13_1980, RCP13_1980,
     * RAF15_1955, RBF15_1955, RCF15_1955,
     * RAP13_1955, RBP13_1955, RCP13_1955,
     * RAP33_1975, RBP33_1975, RCP33_1975

    
      COMPLEX                                                           
     * QAS11_1535, QCS11_1535, 
     * QAS11_1650, QCS11_1650, 
     * QAP11_1440, QCP11_1440,
     * QAP11_1710, QCP11_1710,
     * QAP13_1720, qbp13_1720, QCP13_1720,
     * QAD13_1520, QBD13_1520, QCD13_1520, 
     * QAD13_1700, QBD13_1700, QCD13_1700, 
     * QAD15_1675, QBD15_1675, QCD15_1675, 
     * QAF15_1680, QBF15_1680, QCF15_1680,
     * QAG17_2190, QBG17_2190, QCG17_2190, 
     * QAG19_2250, QBG19_2250, QCG19_2250, 
     * QAH19_2220, QBH19_2220, QCH19_2220,
     * QAI111_2600,QBI111_2600,QCI111_2600, 
     * QAS31_1620, QCS31_1620,
     * QAP31_1910, QCP31_1910,
     * QAP33_1232, QBP33_1232, QCP33_1232,
     * QAP33_1920, QBP33_1920, QCP33_1920,
     * QAD33_1700, QBD33_1700, QCD33_1700,
     * QAD35_1930, QBD35_1930, QCD35_1930,
     * QAF35_1905, QBF35_1905, QCF35_1905,
     * QAF37_1950, QBF37_1950, QCF37_1950,
     * QAH311_2420,QBH311_2420,QCH311_2420,
     * QAP33_1600, QBP33_1600, QCP33_1600,
     * QAF17_1990, QBF17_1990, QCF17_1990,
     * QAF15_2000, QBF15_2000, QCF15_2000,
     * QAP11_2100, QBP11_2100, QCP11_2100,
     * QAF35_2000, QBF35_2000, QCF35_2000,
     * QAP13_1870, QBP13_1870, QCP13_1870,
     * QAP31_1925, QBP31_1925, QCP31_1925,
     * QAP13_1980, QBP13_1980, QCP13_1980,
     * QAF15_1955, QBF15_1955, QCF15_1955,
     * QAP13_1955, QBP13_1955, QCP13_1955,
     * QAP33_1975, QBP33_1975, QCP33_1975,
     * qas31_1900, qcs31_1900,
     * XA0P  ,XA1P  ,XA2P  ,XA3P  ,                                     
     * XA1M  ,XA2M  ,XA3M  ,                                            
     * XB1P  ,XB2P  ,XB3P  ,                                            
     * XB2M  ,XB3M  ,                                                   
     * XC0P  ,XC1P  ,XC2P  ,XC3P  ,                                     
     * XC1M  ,XC2M  ,XC3M                                               
      COMPLEX                                                           
     * RS11_1535, RS11_1650, RP11_1440, RP11_1710, 
     * RP13_1720, RD13_1520, RD13_1700, 
     * RD15_1675, RF15_1680,
     * RG17_2190,
     * RG19_2250, RH19_2220,
     * RI111_2600, 
     * RS31_1620, RS31_1900, RP31_1910,  
     * RP33_1232, RP33_1920, RD33_1700,
     * RD35_1930, RF35_1905,
     * RF37_1950,
     * RH311_2420,
     * RP33_1600, RF17_1990, RF15_2000, RP11_2100, RF35_2000,
     * RP13_1870,RP31_1925,RP13_1980,RF15_1955,RP13_1955,RP33_1975
     
      COMPLEX EPH1,EPH2,EPH3,EPH4,EPH5,EPH6                             
      COMPLEX ADDA0,ADDA1,ADDA2,ADDA3,ADDA4,ADDA5
      COMPLEX ADDB1,ADDB2,ADDB3,ADDB4,ADDB5
      COMPLEX ADDC0,ADDC1,ADDC2,ADDC3,ADDC4,ADDC5
      COMPLEX SUBA0,SUBA1,SUBA2,SUBA3,SUBA4,SUBA5
      COMPLEX SUBB1,SUBB2,SUBB3,SUBB4,SUBB5
      COMPLEX SUBC0,SUBC1,SUBC2,SUBC3,SUBC4,SUBC5
      COMPLEX HNP,HNM,HFP,HFM,HN0,HF0                                   
      COMPLEX X1,X2,Y1,Y2,Y3,Y4,Z1,Z2                                   
      COMPLEX R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11
      REAL HB(6)
                                      
         COMPLEX ZERO                 
         INTEGER EPIREA               
         INTEGER IR                   
         COMPLEX EPRES                
         REAL YMPI,W,PIM,GK           
                                      
         DATA ZERO/(0.,0.)/           
      DATA SQ2        /1.41421  /     
      DATA PI         /3.14159  /     
      DATA XMPROT     /.938     /     
      DATA XMPIP      /.1395    /     
      DATA XMPI0      /.1349    /     
c      DATA XMETA      /.5488    /     
        

        sig0=0.

         DUM=ACOS(EPCOS)                                                
         EPSIN=SIN(DUM)                                                 
         EPSIN2=SIN(DUM/2.)                                             
         EPCOS2=COS(DUM/2.)                                             
         IF(EPIREA.EQ.1.OR.EPIREA.EQ.2)XMPI=XMPI0                                      
         IF(EPIREA.EQ.3.OR.EPIREA.EQ.4)XMPI=XMPIP                                      
                                                                                 
           PENER  =(EPW**2-XMPROT**2+XMPI   **2)/(2.*EPW   )            
           PIMP   =SQRT(PENER   **2 - XMPI   **2)                       
           GENER  =(EPW   **2-XMPROT**2-EPQ2   )/(2.*EPW   )            
           GIMP   =SQRT(GENER**2 + EPQ2   )                             
           GAMK   =(EPW**2-XMPROT**2)/(2.*EPW   )                       
           X=EPCOS                                                      
           P1S  =1.   
           P2S  =3.*X 
           P3S  =(15.*X*X-3.)/2.
           P4S  =(35.*X*X*X-15.*X)/2.
           P5S  =(315.*X*X*X*X-210.*X*X+15.)/8. 
           P6S  =(693.*X*X*X*X*X-630.*X*X*X)/8.
           P7S  =(3003.*X*X*X*X*X*X-3465.*X*X*X*X)/16.
           P2S2 =3.
           P3S2 =15.*X 
           P4S2 =(105.*X*X-15.)/2.
           P5S2 =(315.*X*X*X-105.*X)/2.
           P6S2 =(3465.*X*X*X*X-1890.*X*X)/8.                             
           P7S2 =(18018.*X*X*X*X*X-13860.*X*X*X)/16.
           X0=1./SQRT((PIMP   **2)+(.35**2))                            
           F=(1.+(PIMP   **2)/.71)                                      
           BPIPY0   =1./F                                               
           BPIPY1   =X0*PIMP/F                                       
           BPIPY2   =X0*X0*PIMP*PIMP/F                            
           BPIPY3   =X0*X0*X0*PIMP*PIMP*PIMP/F  

           if (gamk .gt. 0.)then
              FKT=PIMP/GAMK
           else
              FKT=0.
           endif


           Q2=EPQ2                                                      
           W =EPW                                                       
           SINX =EPSIN                                                  
           COSX =EPCOS                                                  
           SINX2=EPSIN2                                                 
           COSX2=EPCOS2                                                 
           IREA=EPIREA                                                  
           IF(IBORN.EQ.1) THEN
            CALL BORNT(Q2,W,SINX,COSX,SINX2,COSX2,HB,IREA,CUT)               
            EPH1B   =HB(1)
            EPH2B   =HB(2)
            EPH3B   =HB(3)
            EPH4B   =HB(4)
            EPH5B   =HB(5)
            EPH6B   =HB(6)                                           
           ELSE
            EPH1B=0.
            EPH2B=0.
            EPH3B=0.
            EPH4B=0.
            EPH5B=0.
            EPH6B=0.
           END IF
                                                                        
           IR=EPIREA         
           YMPI  =  XMPI                                                
           W     =  EPW                                                 
           PIM   =  PIMP                                                
           GK    =  GAMK                                                
      
      RS11_1535 = EPRES(IR,XMPROT,YMPI,W,LS11_1535,WS11_1535,
     * .35,0,1,1,PIM,GK,.TRUE.)
       RS11_1650  = EPRES(IR,XMPROT,YMPI,W,LS11_1650,WS11_1650,
     * .35,0,1,1,PIM,GK,.FALSE.)
       RP11_1440 = EPRES(IR,XMPROT,YMPI,W,LP11_1440,WP11_1440,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP11_1710 = EPRES(IR,XMPROT,YMPI,W,LP11_1710,WP11_1710,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP13_1720 = EPRES(IR,XMPROT,YMPI,W,LP13_1720,WP13_1720,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RD13_1520 = EPRES(IR,XMPROT,YMPI,W,LD13_1520,WD13_1520,
     * .35,2,1,1,PIM,GK,.FALSE.)
       RD13_1700 = EPRES(IR,XMPROT,YMPI,W,LD13_1700,WD13_1700,
     * .35,2,1,1,PIM,GK,.FALSE.)
       RD15_1675 = EPRES(IR,XMPROT,YMPI,W,LD15_1675,WD15_1675,
     * .35,2,2,1,PIM,GK,.FALSE.)
       RF15_1680 = EPRES(IR,XMPROT,YMPI,W,LF15_1680,WF15_1680,
     * .35,3,2,1,PIM,GK,.FALSE.)
       RG17_2190 = EPRES(IR,XMPROT,YMPI,W,LG17_2190,WG17_2190,
     * .35,4,3,1,PIM,GK,.FALSE.)
       RG19_2250 = EPRES(IR,XMPROT,YMPI,W,LG19_2250,WG19_2250,
     * .35,4,4,1,PIM,GK,.FALSE.)
       RH19_2220 = EPRES(IR,XMPROT,YMPI,W,LH19_2220,WH19_2220,
     * .35,5,4,1,PIM,GK,.FALSE.)
       RI111_2600 = EPRES(IR,XMPROT,YMPI,W,LI111_2600,WI111_2600,
     * .35,6,5,1,PIM,GK,.FALSE.)


       RS31_1620 = EPRES(IR,XMPROT,YMPI,W,LS31_1620,WS31_1620,
     * .35,0,1,3,PIM,GK,.FALSE.)
       RS31_1900 = EPRES(IR,XMPROT,YMPI,W,LS31_1900,WS31_1900,
     * .35,0,1,3,PIM,GK,.FALSE.)
       RP31_1910 = EPRES(IR,XMPROT,YMPI,W,LP31_1910,WP31_1910,
     * .35,1,1,3,PIM,GK,.FALSE.)
       RP33_1232 = EPRES(IR,XMPROT,YMPI,W,LP33_1232,WP33_1232,
     * .185,1,1,3,PIM,GK,.FALSE.)
       RP33_1920 = EPRES(IR,XMPROT,YMPI,W,LP33_1920,WP33_1920,
     * .35,1,1,3,PIM,GK,.FALSE.)
       RD33_1700 = EPRES(IR,XMPROT,YMPI,W,LD33_1700,WD33_1700,
     * .35,2,1,3,PIM,GK,.FALSE.)
       RD35_1930 = EPRES(IR,XMPROT,YMPI,W,LD35_1930,WD35_1930,
     * .35,2,2,3,PIM,GK,.FALSE.)
       RF35_1905 = EPRES(IR,XMPROT,YMPI,W,LF35_1905,WF35_1905,
     * .35,3,2,3,PIM,GK,.FALSE.)
       RF37_1950 = EPRES(IR,XMPROT,YMPI,W,LF37_1950,WF37_1950,
     * .35,3,3,3,PIM,GK,.FALSE.)
       RH311_2420 = EPRES(IR,XMPROT,YMPI,W,LH311_2420,WH311_2420,
     * .35,5,5,3,PIM,GK,.FALSE.)
                                                                         
C* THE FOLOWING ARE ONE OR TWO STAR RESONANCES NOT INCLUDED IN THE ABOVE
       RP33_1600 = EPRES(IR,XMPROT,YMPI,W,LP33_1600,WP33_1600,
     * .35,1,1,3,PIM,GK,.FALSE.)
       RF17_1990 = EPRES(IR,XMPOT,YMPI,W,LF17_1990,WF17_1990,
     * .35,3,3,1,PIM,GK,.FALSE.)
       RF15_2000 = EPRES(IR,XMPOT,YMPI,W,LF15_2000,WF15_2000,
     * .35,3,2,1,PIM,GK,.FALSE.)
       RP11_2100 = EPRES(IR,XMPROT,YMPI,W,LP11_2100,WP11_2100,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RF35_2000 = EPRES(IR,XMPROT,YMPI,W,LF35_2000,WF35_2000,
     * .35,3,2,3,PIM,GK,.FALSE.)

C*  THE FOLOWING ARE RESONANCES PREDICTED BY THE QUARK MODEL ONLY
       RP13_1870 = EPRES(IR,XMPROT,YMPI,W,LP13_1870,WP13_1870,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP31_1925 = EPRES(IR,XMPROT,YMPI,W,LP31_1925,WP31_1925,
     * .35,1,1,3,PIM,GK,.FALSE.)
       RP13_1980 = EPRES(IR,XMPROT,YMPI,W,LP13_1980,WP13_1980,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RF15_1955 = EPRES(IR,XMPROT,YMPI,W,LF15_1955,WF15_1955,
     * .35,3,2,1,PIM,GK,.FALSE.)
       RP13_1955 = EPRES(IR,XMPROT,YMPI,W,LP13_1955,WP13_1955,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP33_1975 = EPRES(IR,XMPROT,YMPI,W,LP33_1975,WP33_1975,
     * .35,1,1,3,PIM,GK,.FALSE.)
                                                     
c           CALL BACK(IB)
                                                                        
           IF(EPIREA   .EQ.1)THEN                                       
             XA0P    = BGA0P0*BPIPY0   +ZERO                            
             XA1P    = BGA1P0*BPIPY1   +ZERO                            
             XA2P    = BGA2P0*BPIPY2   +ZERO                            
	     XA3P    = BGA3P0*BPIPY3   +ZERO
	     xa4p    = zero
	     xa5p    = zero
	     
             XB1P    = BGB1P0*BPIPY1   +ZERO 
             XB2P    = BGB2P0*BPIPY2   +ZERO 
             XB3P    = BGB3P0*BPIPY3   +ZERO
	     xb4p    = zero
	     xb5p    = zero
	     
             XC0P    = BGC0P0*BPIPY0   +ZERO                            
             XC1P    = BGC1P0*BPIPY1   +ZERO                            
             XC2P    = BGC2P0*BPIPY2   +ZERO
             XC3P    = BGC3P0*BPIPY3   +ZERO
             xc4p     =  zero
             xc5p     = zero
                                            
             XA1M    = BGA1M0*BPIPY1   +ZERO
             XA2M    = BGA2M0*BPIPY2   +ZERO
             XA3M    = BGA3M0*BPIPY3   +ZERO
	     xa4m    = zero
	     xa5m     = zero
	     xa6m     = zero
	     	
             XB2M    = BGB2M0*BPIPY2   +ZERO
             XB3M    = BGB3M0*BPIPY3   +ZERO                            
	     xb4m     = zero
	     xb5m    = zero
	     xb6m     = zero
	     
             XC1M    = BGC1M0*BPIPY1   +ZERO                            
             XC2M    = BGC2M0*BPIPY2   +ZERO
             XC3M    = BGC3M0*BPIPY3   +ZERO
             xc4m    =  zero
             xc5m     = zero
             xc6m    = zero
                            
           ELSE 

             XA0P    = BGA0PP*BPIPY0   +ZERO 
             XA1P    = BGA1PP*BPIPY1   +ZERO 
             XA2P    = BGA2PP*BPIPY2   +ZERO 
             XA3P    = BGA3PP*BPIPY3   +ZERO 

             XB1P    = BGB1PP*BPIPY1   +ZERO 
             XB2P    = BGB2PP*BPIPY2   +ZERO 
             XB3P    = BGB3PP*BPIPY3   +ZERO 

             XC0P    = BGC0PP*BPIPY0   +ZERO 
             XC1P    = BGC1PP*BPIPY1   +ZERO 
             XC2P    = BGC2PP*BPIPY2   +ZERO 
             XC3P    = BGC3PP*BPIPY3   +ZERO 
                                    
             XA1M    = BGA1MP*BPIPY1   +ZERO 
             XA2M    = BGA2MP*BPIPY2   +ZERO 
             XA3M    = BGA3MP*BPIPY3   +ZERO 

             XB2M    = BGB2MP*BPIPY2   +ZERO 
             XB3M    = BGB3MP*BPIPY3   +ZERO 

             XC1M    = BGC1MP*BPIPY1   +ZERO 
             XC2M    = BGC2MP*BPIPY2   +ZERO 
             XC3M    = BGC3MP*BPIPY3   +ZERO 

           END IF
 
C* IT=1 FOR EXPT ONLY, IT=2 FOR EXPT+ONE OR TWO STAR FROM QKM
C* IT=2 FOR EXPT+ALL THE OTHER FROM QKM
C* IT=4 FOR ALL FROM QKM
C* C2,C3 ARE USED TO TURN THE QKM PREDICTION OFF AND ON

         IF(IT.EQ.1) THEN
           CALL EXPA(EPIREA,EPQ2)
         IF(IP11.EQ.2) THEN
           RAP11_1440=0.0
         ENDIF
         IF(IP11.EQ.3.OR.IP11.EQ.4) THEN
c           CALL QKMA(EPIREA,EPQ2)
         ENDIF
         C2=0.
         C3=0.
         GO TO 35
         ENDIF

         IF(IT.EQ.2) THEN
         C2=1.
         C3=0.
c         CALL QKMA(EPIREA,EPQ2)
         CALL EXPA(EPIREA,EPQ2)
         GO TO 25
         ENDIF

         IF (IT.EQ.3) THEN
c         CALL QKMA(EPIREA,EPQ2)
         CALL EXPA(EPIREA,EPQ2)
         C2=1.
         C3=1.
         ENDIF

         IF(IT.EQ.4) THEN
         C2=0.
         C3=0.
C* here I get ride of the terms which is not used in
c* the EXPT
c         CALL QKMA(EPIREA,EPQ2)
C*         GO TO 35
         ENDIF

C* THIS FOLLOWING RESONANCES ARE PREDICTED BY QK MODEL ONLY
15       QAP13_1870   =RP13_1870   * RAP13_1870
         QBP13_1870   =RP13_1870   * RBP13_1870
         QCP13_1870   =RP13_1870   * 0.0

         QAP31_1925   =RP31_1925   * RAP31_1925
         QBP31_1925   =RP31_1925   * RBP31_1925
         QCP31_1925   =RP31_1925   * 0.0
 
         QAP13_1980   =RP13_1980   * RAP13_1980
         QBP13_1980   =RP13_1980   * RBP13_1980
         QCP13_1980   =RP13_1980   * 0.0
        
         QAF15_1955   =RF15_1955   * RAF15_1955
         QBF15_1955   =RF15_1955   * RBF15_1955
         QCF15_1955   =RF15_1955   * 0.0
         
         QAP13_1955   =RP13_1955   * RAP13_1955
         QBP13_1955   =RP13_1955   * RBP13_1955
         QCP13_1955   =RP13_1955   * 0.0
     
         QAP33_1975   =RP33_1975   * RAP33_1975
         QBP33_1975   =RP33_1975   * RBP33_1975
         QCP33_1975   =RP33_1975   * 0.0
         

C*  THE FOLLOWING RESONANCES ARE SEEN AS ONE OR TWO STAR
25       QAP33_1600   =RP33_1600   * RAP33_1600
         QBP33_1600   =RP33_1600   * RBP33_1600
         QCP33_1600   =RP33_1600   * 0.0

         QAF17_1990   =RF17_1990   * RAF17_1990
         QBF17_1990   =RF17_1990   * RBF17_1990
         QCF17_1990   =RF17_1990   * 0.0
        
         QAF15_2000   =RF15_2000   * RAF15_2000
         QBF15_2000   =RF15_2000   * RBF15_2000
         QCF15_2000   =RF15_2000   * 0.0
    
         QAP11_2100   =RP11_2100   * RAP11_2100
         QBP11_2100   =RP11_2100   * RBP11_2100
         QCP11_2100   =RP11_2100   * 0.0

         QAF35_2000   =RF35_2000   * RAF35_2000
         QBF35_2000   =RF35_2000   * RBF35_2000
         QCF35_2000   =RF35_2000   * RCF35_2000

         
C* THE FOLLWING RESONANCES ARE SEEN BY EXPERIMENT                                                                        
35       QAS11_1535   =RS11_1535   * RAS11_1535 
         QCS11_1535   =RS11_1535   * RCS11_1535

         QAS11_1650   =RS11_1650   * RAS11_1650
         QCS11_1650   =RS11_1650   * RCS11_1650

         QAP11_1440   =RP11_1440   * RAP11_1440  
         QCP11_1440   =RP11_1440   * RCP11_1440  

         QAP13_1720   =RP13_1720   * RAP13_1720
         QBP13_1720   =RP13_1720   * RBP13_1720
         QCP13_1720   =RP13_1720   * RCP13_1720

         QAD13_1520   =RD13_1520   * RAD13_1520
         QBD13_1520   =RD13_1520   * RBD13_1520
         QCD13_1520   =RD13_1520   * RCD13_1520

         QAD13_1700   =RD13_1700   * RAD13_1700
         QBD13_1700   =RD13_1700   * RBD13_1700
         QCD13_1700   =RD13_1700   * RCD13_1700

         QAD15_1675   =RD15_1675   * RAD15_1675
         QBD15_1675   =RD15_1675   * RBD15_1675
         QCD15_1675   =RD15_1675   * RCD15_1675

         QAF15_1680   =RF15_1680   * RAF15_1680  
         QBF15_1680   =RF15_1680   * RBF15_1680
         QCF15_1680   =RF15_1680   * RCF15_1680

         QAG17_2190   =RG17_2190   * RAG17_2190
         QBG17_2190   =RG17_2190   * RBG17_2190
         QCG17_2190   =RG17_2190   * RCG17_2190

         QAG19_2250   =RG19_2250   * RAG19_2250
         QBG19_2250   =RG19_2250   * RBG19_2250
         QCG19_2250   =RG19_2250   * RCG19_2250

         QAH19_2220   =RH19_2220   * RAH19_2220
         QBH19_2220   =RH19_2220   * RBH19_2220
         QCH19_2220   =RH19_2220   * RCH19_2220

         QAI111_2600  =RI111_2600  * RAI111_2600
         QBI111_2600  =RI111_2600  * RBI111_2600
         QCI111_2600  =RI111_2600  * RCI111_2600

         QAS31_1620   =RS31_1620   * RAS31_1620  
         QCS31_1620   =RS31_1620   * RCS31_1620  

         QAS31_1900   =RS31_1900   * RAS31_1900  
         QCS31_1900   =RS31_1900   * RCS31_1900  

         QAP31_1910   =RP31_1910   * RAP31_1910  
         QCP31_1910   =RP31_1910   * RCP31_1910  

         QAP33_1232   =RP33_1232   * RAP33_1232  
         QBP33_1232   =RP33_1232   * RBP33_1232  
         QCP33_1232   =RP33_1232   * RCP33_1232

         QAP33_1920   =RP33_1920   * RAP33_1920  
         QBP33_1920   =RP33_1920   * RBP33_1920  
         QCP33_1920   =RP33_1920   * RCP33_1920  


         QAD33_1700   =RD33_1700   * RAD33_1700  
         QBD33_1700   =RD33_1700   * RBD33_1700  
         QCD33_1700   =RD33_1700   * RCD33_1700  

         QAD35_1930   =RD35_1930   * RAD35_1930
         QBD35_1930   =RD35_1930   * RBD35_1930
         QCD35_1930   =RD35_1930   * RCD35_1930

         QAF35_1905   =RF35_1905   * RAF35_1905
         QBF35_1905   =RF35_1905   * RBF35_1905
         QCF35_1905   =RF35_1905   * RCF35_1905

         QAF37_1950   =RF37_1950   * RAF37_1950
         QBF37_1950   =RF37_1950   * RBF37_1950
         QCF37_1950   =RF37_1950   * RCF37_1950

         QAH311_2420  =RH311_2420  * RAH311_2420
         QBH311_2420  =RH311_2420  * RBH311_2420
         QCH311_2420  =RH311_2420  * RCH311_2420

c	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c	The following were not initialized in my copy of AO (RCM)
	qap11_1710=zero
    	qcp11_1710=zero
c	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
    				
         XA0P=XA0P + QAS11_1535 + QAS11_1650 + QAS31_1620 + QAS31_1900
         XC0P=XC0P + QCS11_1535 + QCS11_1650 + QCS31_1620 + QCS31_1900
                                                       
         XA1M=XA1M + QAP11_1440 + QAP11_1710 + QAP31_1910
     *       +C2*(QAP11_2100)+C3*(QAP31_1925)
         XC1M=XC1M + QCP11_1440 + QCP11_1710 + QCP31_1910
     *       +C2*(QCP11_2100)+C3*(QCP31_1925)
                                                                   
         XA1P=XA1P + QAP13_1720 + QAP33_1232 + QAP33_1920
     *+C2*(QAP33_1600)+C3*(QAP13_1870+QAP13_1980+QAP13_1955+QAP33_1975)
         XB1P=XB1P + QBP13_1720 + QBP33_1232 + QBP33_1920
     *+C2*(QBP33_1600)+C3*(QBP13_1870+QBP13_1980+QBP13_1955+QBP33_1975)
         XC1P=XC1P + QCP13_1720 + QCP33_1232 + QCP33_1920
     *+C2*(QCP33_1600)+C3*(QCP13_1870+QCP13_1980+QCP13_1955+QCP33_1975)
                                                                        
         XA2M = XA2M + QAD13_1520 + QAD13_1700 + QAD33_1700
         XB2M = XB2M + QBD13_1520 + QBD13_1700 + QBD33_1700
         XC2M = XC2M + QCD13_1520 + QCD13_1700 + QCD33_1700
         
         XA2P = XA2P + QAD15_1675 + QAD35_1930
         XB2P = XB2P + QBD15_1675 + QBD35_1930
         XC2P = XC2P + QCD15_1675 + QCD35_1930

         XA3M = XA3M + QAF15_1680 + QAF35_1905
     *        +C2*(QAF15_2000+QAF35_2000)+C3*(QAF15_1955)    
         XB3M = XB3M + QBF15_1680 + QBF35_1905
     *        +C2*(QBF15_2000+QBF35_2000)+C3*(QBF15_1955)    
         XC3M = XC3M + QCF15_1680 + QCF35_1905
     *        +C2*(QCF15_2000+QCF35_2000)+C3*(QCF15_1955)    
                                      
         XA3P = XA3P + QAF37_1950
     *        +C2*(QAF17_1990)
         XB3P = XB3P + QBF37_1950
     *        +C2*(QBF17_1990)
         XC3P = XC3P + QCF37_1950
     *        +C2*(QCF17_1990)

         XA4M = XA4M + QAG17_2190
         XB4M = XB4M + QBG17_2190
         XC4M = XC4M + QCG17_2190
                                  
         XA4P = XA4P + QAG19_2250
         XB4P = XB4P + QBG19_2250
         XC4P = XC4P + QCG19_2250

         XA5M = XA5M + QAH19_2220
         XB5M = XB5M + QBH19_2220
         XC5M = XC5M + QCH19_2220

         XA5P = XA5P + QAH311_2420
         XB5P = XB5P + QBH311_2420
         XC5P = XC5P + QCH311_2420

         XA6M = XA6M + QAI111_2600
         XB6M = XB6M + QBI111_2600
         XC6M = XC6M + QCI111_2600

         SUBA0=XA0P-XA1M 
         SUBA1=XA1P-XA2M 
         SUBA2=XA2P-XA3M 
         SUBA3=XA3P-XA4M 
         SUBA4=XA4P-XA5M 
         SUBA5=XA5P-XA6M 
                         
         SUBC0=XC0P-XC1M 
         SUBC1=XC1P-XC2M 
         SUBC2=XC2P-XC3M 
         SUBC3=XC3P-XC4M 
         SUBC4=XC4P-XC5M 
         SUBC5=XC5P-XC6M 
                         
         SUBB1=XB1P-XB2M 
         SUBB2=XB2P-XB3M 
         SUBB3=XB3P-XB4M 
         SUBB4=XB4P-XB5M 
         SUBB5=XB5P-XB6M 
                         
         ADDA0=XA0P+XA1M 
         ADDA1=XA1P+XA2M 
         ADDA2=XA2P+XA3M 
         ADDA3=XA3P+XA4M 
         ADDA4=XA4P+XA5M 
         ADDA5=XA5P+XA6M 
                         
         ADDC0=XC0P+XC1M 
         ADDC1=XC1P+XC2M 
         ADDC2=XC2P+XC3M 
         ADDC3=XC3P+XC4M 
         ADDC4=XC4P+XC5M 
         ADDC5=XC5P+XC6M 
                         
         ADDB1=XB1P+XB2M 
         ADDB2=XB2P+XB3M 
         ADDB3=XB3P+XB4M 
         ADDB4=XB4P+XB5M 
         ADDB5=XB5P+XB6M 
c
C   Now do the partial wave expansion for the Walker amplitudes
C
       EPH1   =SUBB1*(       -P2S2   )                                  
     *       + SUBB2*(P2S2   -P3S2   )                                  
     *       + SUBB3*(P3S2   -P4S2   )                                  
     *       + SUBB4*(P4S2   -P5S2   )                                  
     *       + SUBB5*(P5S2   -P6S2   )
       EPH1   = EPH1   * EPSIN   * EPCOS2   /SQ2                           
                                                                        
                                                                        
       EPH2  = SUBA0   *(       -P1S   )                                
     *       + SUBA1   *(P1S    -P2S   )                                
     *       + SUBA2   *(P2S    -P3S   )                                
     *       + SUBA3   *(P3S    -P4S   )
     *       + SUBA4   *(P4S    -P5S   )
     *       + SUBA5   *(P5S    -P6S   )
       EPH2  = EPH2    * EPCOS2   * SQ2

                                                                        
       EPH3   = ADDB1    *(       +P2S2   )                             
     *       +  ADDB2    *(P2S2   +P3S2   )                             
     *       +  ADDB3    *(P3S2   +P4S2   )
     *       +  ADDB4    *(P4S2   +P5S2   )
     *       +  ADDB5    *(P5S2   +P6S2   )
       EPH3   = EPH3   *  EPSIN   * EPSIN2   /SQ2     
                                                                        
                                                                        
       EPH4  = ADDA0     * (       +P1S   )                              
     *       + ADDA1     * (P1S    +P2S   )                             
     *       + ADDA2     * (P2S    +P3S   )                             
     *       + ADDA3     * (P3S    +P4S   )
     *       + ADDA4     * (P4S    +P5S   )
     *       + ADDA5     * (P5S    +P6S   )
       EPH4   = EPH4   * EPSIN2   * SQ2              
                                                                        
                                                                        
       EPH5   =SUBC0     *(       -P1S   )                              
     *       + SUBC1     *(P1S    -P2S   )                              
     *       + SUBC2     *(P2S    -P3S   )                              
     *       + SUBC3     *(P3S    -P4S   )
     *       + SUBC4     *(P4S    -P5S   )
     *       + SUBC5     *(P5S    -P6S   )
       EPH5   = EPH5   * EPCOS2   * SQ2                
                                                                        
                                                                        
       EPH6  = ADDC0     * (       +P1S   )                              
     *       + ADDC1     * (P1S    +P2S   )                             
     *       + ADDC2     * (P2S    +P3S   )                             
     *       + ADDC3     * (P3S    +P4S   )
     *       + ADDC4     * (P4S    +P5S   )
     *       + ADDC5     * (P5S    +P6S   )
        EPH6   =  EPH6   * EPSIN2   * SQ2              
                                                                        
C
C      Now add the Born term contributions. Remember, these are real functions!
C                                                                        
         EPH1   =EPH1   +CMPLX(EPH1B   ,0.)                             
         EPH2   =EPH2   +CMPLX(EPH2B   ,0.)                             
         EPH3   =EPH3   +CMPLX(EPH3B   ,0.)                             
         EPH4   =EPH4   +CMPLX(EPH4B   ,0.)                             
         EPH5   =EPH5   +CMPLX(EPH5B   ,0.)                             
         EPH6   =EPH6   +CMPLX(EPH6B   ,0.)               
   
         
C
C       Now we are ready to calculate all the observables
C
C          First the unpolarized transverse term
C              
                                                                        
           SIGU   =(CABS(EPH1   )**2                                    
     *            +CABS(EPH2   )**2                                     
     *            +CABS(EPH3   )**2                                     
     *            +CABS(EPH4   )**2)/2                                  
C
C            Next the transverse polarized term
C                                                                       
           SIGT   =REAL(EPH2   *CONJG(EPH3   )                          
     *            -     EPH1   *CONJG(EPH4   ))                         
C
C             And the longitudinal term
C                                                                        
           SIGL   =CABS(EPH5   )**2                                     
     *            +CABS(EPH6   )**2                                     
C
C             Finally, the longitudinal-transverse interference term
C
           SIGI   =SQ2*REAL(EPH5   *CONJG(EPH4   -EPH1   )              
     *            +         EPH6   *CONJG(EPH3   +EPH2   ))             
C
C              Sum up all the unpolarized cross section terms
C
                                                                        
             sig0  = SIGU       + EPEPS*SIGL
           SIGMAO   = sig0
     *             +EPEPS       *SIGT       *COS(EPPHI   *PI/90.)       
     *             +SQRT(EPEPS       *(1+EPEPS       )/2.)*SIGI         
     *             *COS(EPPHI   *PI/180.)                               



C
         if(POL_ELEC .eq. 0. .and. POL_TARG .eq. 0.)then
            sigmao=sigmao*FKT

            return
         endif



C
C              Convert Walker amplitudes to Bartl & Majerotto spin-flip and non spin-flip amplitudes
C
           HNP=(EPH4+EPH1)/SQ2                                          
           HNM=(EPH4-EPH1)/SQ2                                          
           HFP=(EPH3-EPH2)/SQ2                                          
           HFM=(EPH3+EPH2)/SQ2                                          
           HN0=EPH5                                                     
           HF0=EPH6                                                     
C
C               Define the appreviations for the various polarized cross section terms
C               Polarized beam, polarized target, polarized beam - polarized target
C
           X1=HF0*CONJG(HNP) + HN0*CONJG(HFP)                           
           X2=HFM*CONJG(HNP) + HNM*CONJG(HFP)                           
           Y1=HNP*CONJG(HFP) + HNM*CONJG(HFM)                           
           Y2=HNM*CONJG(HFM) - HNP*CONJG(HFP)                           
           Y3=HN0*CONJG(HF0)                                            
           Y4=HN0*CONJG(HFM) - HF0*CONJG(HNM)                           
	   Z1=HN0*CONJG(HNP) - HF0*CONJG(HFP)
           Z2=HNM*CONJG(HNP) - HFM*CONJG(HFP)                           
C
C               These are the terms for the recoil polarizations
C
	   R1=CONJG(HN0)*HF0
           R2=CONJG(HNP)*HFP - CONJG(HNM)*HFM
           R3=CONJG(HN0)*HFM + CONJG(HF0)*HNP
           R4=CONJG(HF0)*HNM - CONJG(HN0)*HFM
           R5=CONJG(HN0)*HNP + CONJG(HF0)*HFP
           R6=CONJG(HNP)*HFP + CONJG(HNM)*HFM
           R7=CONJG(HFM)*HNP + CONJG(HFP)*HNM
           R8=CONJG(HFP)*HFM + CONJG(HNP)*HNM
           R9=CONJG(HNP)*HFM - CONJG(HNM)*HFP
           R10=CONJG(HNP)*HNM + CONJG(HFP)*HFM
           R11=CONJG(HN0)*HNM + CONJG(HF0)*HFM
           

           XILA=SQRT(2.*EPEPS*(EPEPS+1))                                
           SIGE=-POL_ELEC*SQRT(2*EPEPS*(1-EPEPS))*SIN(EPPHI*PI/180.)*         
     *     AIMAG(HN0*CONJG(HNM)+HF0*CONJG(HFM))

                         
C           SIGPT=PX*(-XILA*SIN(EPPHI*PI/180)*AIMAG(X1)-EPEPS*            
C     *         SIN(EPPHI*PI/90.)*AIMAG(X2))                             
C     *       -PY*(AIMAG(Y1)+EPEPS*COS(EPPHI*PI/90.)*AIMAG(Y2)           
C     *          +2.*EPEPS*AIMAG(Y3)+XILA*COS(EPPHI*PI/180)*AIMAG(Y4))   
C     *         +PZ*(EPEPS*SIN(EPPHI*PI/90.)*AIMAG(Z2)+                  
C     *              XILA*SIN(EPPHI*PI/180.)*AIMAG(Z1))                  
C       SIGET=-POL_ELEC*(-PX*(SQRT(2.*EPEPS*(1-EPEPS))*COS(EPPHI*PI/180)  
C     *               *REAL(X1)+SQRT(1-EPEPS**2)*REAL(X2))               
C     *       +PY*SQRT(2.*EPEPS*(1-EPEPS))*SIN(EPPHI*PI/180.)*REAL(Y4)   
C     *       +PZ*(SQRT(1-EPEPS**2)*REAL(Z2)+SQRT(2*EPEPS*(1-EPEPS))*    
C     *        COS(EPPHI*PI/180.)*REAL(Z1)))                             
C           SIGMA = SIGMA + SIGE + SIGPT + SIGET
C           SIGMA=SIGMA*FKT                                              
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C         The target is polarized along the beam axis.
c                                       PS = -1 anti-parallel to beam
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
           IF(EPQ2.LT.1.E-06.OR.EPEPS.LT.1.E-06)THEN
                THETA_GAM = 0.
           ELSE
           ENUE = (EPW**2 - XMPROT**2 + EPQ2)/2./XMPROT
           THETA_ELEC = 2.*ATAN(SQRT(EPQ2/2./(EPQ2+ENUE**2)*
     *                  (1.-EPEPS)/EPEPS))
           E_ELEC_SCATT = -ENUE/2 + SQRT((ENUE/2.)**2 + 
     *                   EPQ2/4./SIN(THETA_ELEC/2.)**2)
           E_ELEC_PRIM = E_ELEC_SCATT + ENUE
           THETA_GAM = ATAN(SIN(THETA_ELEC)/
     *                 (E_ELEC_PRIM/E_ELEC_SCATT - COS(THETA_ELEC)))
           THETA_GAM = THETA_GAM*180./PI
           ENDIF

           POL_X = POL_TARG*SINDD(THETA_GAM)*COSDD(EPPHI)
           POL_Y = POL_TARG*SINDD(THETA_GAM)*SINDD(EPPHI)
           POL_Z = POL_TARG*COSDD(THETA_GAM)

           SIGPT = - POL_X*(SINDD(EPPHI)*
     1  SQRT(2.*EPEPS*(1. + EPEPS))*AIMAG(X1) + 
     2  EPEPS*SINDD(2.*EPPHI)*AIMAG(X2)) - 
     3  POL_Y*(AIMAG(Y1) + 
     4  EPEPS*COSDD(2.*EPPHI)*AIMAG(Y2) + 
     5  2.*EPEPS*AIMAG(Y3) + 
     6  SQRT(2.*EPEPS*(1.+EPEPS))*COSDD(EPPHI)*AIMAG(Y4)) + 
     7  POL_Z*(EPEPS*SINDD(2.*EPPHI)*AIMAG(Z2) + 
     8  SQRT(2.*EPEPS*(1. + EPEPS))*SINDD(EPPHI)*AIMAG(Z1))

        SIGET = POL_ELEC*POL_TARG*(SINDD(THETA_GAM)*(COSDD(EPPHI))**2*
     1  (REAL(Y4) + SQRT(2.*EPEPS*(1.-EPEPS))*REAL(X1))+
     2  COSDD(EPPHI)*(SINDD(THETA_GAM)*SQRT(1.-EPEPS**2)*REAL(X2) - 
     3  COSDD(THETA_GAM)*SQRT(2.*EPEPS*(1.-EPEPS))*REAL(Z1)) -
     4  SINDD(THETA_GAM)*REAL(Y4) - COSDD(THETA_GAM)*
     5  SQRT(1.-EPEPS**2)*REAL(Z2))
   
           SIGMAO = SIGMAO + SIGE + SIGPT + SIGET 

 99        SIGMAO = SIGMAO*FKT
           
           RETURN                                                       
           END                                                          

*****************************************************************************

      real function sindd(x)
      implicit none

      real pi,x
      pi=4.*atan(1.)
      sindd=sin(pi*x/180.)
      return
      end
*****************************************************************************
*****************************************************************************


*****************************************************************************
*****************************************************************************

      real function cosdd(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      cosdd=cos(pi*x/180.)
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************

      real function tandd(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      tandd=tan(pi*x/180.)
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************
      
      real function asind(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      asind=180.*asin(x)/pi
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************
      
      real function acosd(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      acosd=180.*acos(x)/pi
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************
      
      real function atand(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      atand=180.*atan(x)/pi
      return
      end
*****************************************************************************
