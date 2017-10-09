
      program aaorad_gen

c     This program makes an n-tuple that can be used with Paw to
c     make distributions of energies, angles, resonance
c     mass resulting from internal bremmstrahlung associated with pion
c     production on a proton. The exact integration formula of Mo and Tsai
c     is used.
c     The n-tuple contains the photon energy(EG), the true hadronic invariant
c     mass (W), the components of the proton momentum (PPX, PPY, PPZ),
c     the proton energy (EP), the pion momentum (PPIX, PPIY, PPIZ) and
c     pion energy (EPI), the  angles for the hadronic
c     decay in the hadronic frame (CSTCM, PHICM), the missing mass (MM),
c     and the photon angles relative to the q vector, (CSTHK, PHIK).
c
c     This program forces the monte carlo to concentrate on the regions
c     of photon emission along the directions of the incident and
c     scattered electrons.
c
c     The electrons are radiated as they pass through the target. Resolution
c     of detectors is not folded into the results.  If this is desired it
c     should be done with a second program that can operate on the n-tuple
c     and make a new version.

      implicit none
      
#include "bcs.inc"
#include "mc.inc"
#include "names.inc"

      COMMON/ALPHA/ ALPHA,PI,MP,MPI,MEL,WG,EPIREA,TH_OPT,RES_OPT
      common /radcal/T0,es,ep,ps,pp,rs,rp,u0,pu,uu,cst0,snt0,csths,csthp
     * ,snths,snthp,pdotk,sdotk
      common /random/idum

      real lundmass(4)
      integer timevals(8)
      integer sseed(64)
      integer ssize
      real*8 ek,Tk,delta
      real*8 alpha,pi,mp,mpi,mel,wg,T0
      real*8 es,ep,ps,pp,rs,rp,u0,pu,uu,pdotk,sdotk
      real*8 cst0,snt0,csths,csthp,snths,snthp

      real csran,csrng,csrnge,csrngb,delphi
      real csthcm
      real csthcm_max
      real cstk
      real cstk1,cstk2
      real cstmp
      real csdotk,cpdotk,cqdotk
      real delinf
      real deltar
      real ek_max
      real ekmax
      real ekx,eky,ekz
      real epeps
      real epi
      real epmax
      real eprng
      real eprot
      real epw
      real ep_min,ep_max,ep_test
      real ep_sav
      real events
      real f
      real fkt
      real fmcall
      real g
      real jacob
      real kexp
      real kfac
      real mcfac
      real mm2
      real mm_exp,mm_cut
      real mpfac
      real mpi0
      real mpip
      real meta
      real mpi_s
      real myran
      real nu
      real phik
      real ppx,ppy,ppz
      real ppix,ppiy,ppiz
      real phicm
      real phicm_max
      real phir
      real px,py
      real q0
      real q2
      real q2_min,q2_max,q2max
      real qsq
      real qvecx
      real qvecz
      real ran
      real reg1,reg2,reg3,reg4
      real rn1,rn2
      real rotc,rots
      real rtest
      real s
      real s1
      real s2
      real sigi
      real sigl
      real sigma
      real sigma0
      real signr
      real sigr,sig_ratio
      real sigr_max
      real sigr1
      real*8 sig_tot,sig_sum
      real sigt
      real sigu
      real sigip, asym_p
      real sntk
      real sp
      real spence
      real stest
      real t_elapse
      real th0
      real theta
      real tk_max
      real tp
      real tries
      real itime1, itime2
      real ts
      real uek
      real uq2, uq2_min,uq2_max,uq2rng
      real w2
      real wreal, w_min_input
      real x1
      real x2
      real targs,targp,xs,eloss,gxs,xtest,ebeam,t_targ,bfac,r_targ,temp
      real sig_int, hydrogen_rad, vertex_x, vertex_y, vertex_z, vz

      integer dismc(6,100)
      integer intreg

      integer epirea, th_opt, res_opt
      integer i
      integer ir1
c      integer iext
      integer j
      integer jj
      integer mcall
      integer mcall_max
      integer nprint,ntell,ntold
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer get_spin
      integer flag_ehel
      integer ehel
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer*4 idum
      integer*4 ntries
      
      real cfac,asig

c     Parameters for the n-tuple, which is named func1 and contains
c     15 elements per event.

      common /pawc/h(5000000)
      integer h,n,nevent,nmax,lrecl,istat,icycle
      parameter (n=32)
      real*4 ntp(n)
      
c     tag is the an array of names for the variables in the n-tuple.

c     character*1 ich1
      character*3 month
      character*2 day
      character*2 year
      character*5 tag(n)
      character*13 filerz
      character*13 file_out
      character*13 lund_out
      character*13 file_sum
      character*13 file_bos
      character*8 recname
      character*28  ctime

      data tag /' ES  ',' EP  ','THETE','  W  ','WREAL',
     * ' PPX ',' PPY ',' PPZ ','EPROT','PPIX ',
     *  'PPIY ','PPIZ ',' EPI ','CSTCM','PHICM',' MM  ','  EG ',
     * 'CSTHK','PHIK ', ' QX  ',' QZ  ',' Q0  ','CSTHE',' EGX ',
     * ' EGY ',' EGZ ', ' VX  ', ' VY  ', ' VZ  ', ' Q2  ',
     * ' HEL ','ASYM '/


      DATA PI   /3.1415926/
      DATA MPIP /.1395/
      DATA MPI0 /.1349/
      DATA META /.547853/
      DATA MP   /.938/
      DATA MEL  /.511E-3/

      data filerz   /'aao_rad.hbook'/
      data file_out /'aao_rad.out'/
      data lund_out /'aao_rad.lund'/
      data file_sum /'aao_rad.sum'/
      data file_bos /'aao_rad.bos'/
      data ctime    /'                            '/

      do j=1,6
        do i=1,100
         dismc(j,i)=0
        enddo
      enddo

c     set up parameters for breaking the monte-carlo integration region
c     over csthk into 5 parts:

      csrng=.04

c      Region sizes suggested for 4 GeV: reg1=.23, reg2=.14, reg3=.11, reg4=.10

 14   print *, 'Enter theory_opt,resonance_opt: (1=AO,4=MAID98,5=MAID2000)'
      read(5,*) th_opt
c
      print *, th_opt

      print *, 'Enter 1 for polarized electron, 0 for unpolarized electron'
      read(5,*) flag_ehel
c      
      print *, flag_ehel
      res_opt = 0
      
      print *, 'Input the sizes of the integration regions'
      read(5,*)reg1,reg2,reg3,reg4
      
      reg2 = reg1 + reg2
      reg3 = reg2 + reg3
      reg4 = reg3 + reg4
      
      if (reg4 .gt. .95) then
        write(6,*)' The sum of the region sizes must be less than .95'
        go to 14
      endif
      
      alpha = 1/137.

c     set up parameters for bos bank input to GSIM

      print *, ' Input 2 for two charged particles in the bos bank'
      print *, ' Input 4 to get the neutral hadron and photon as well'              
      read(5,*)npart
      print *, ' Input=',npart
      
      if (npart .ne. 4) npart=2

      q(1) 	= -1
      id(1)	= 3	!Geant ID, e-
      pdgid(1)  = 11	!PDG ID, e-
      
      if (npart .eq. 4)then
         q(4)	   = 0
         id(4)     = 1	!Geant ID, photon
         pdgid(4)  = 22	!PDG ID, photon
      endif
c   
c     Choose whether a neutral or charged pion is made in the reaction

 5    print *, ' Input epirea (1 for pi0, 3 for pi+, 5 for eta)'
      read(5,*)epirea
      print *, 'epiarea=', epirea
      
      print *, ' Input a limit on the error in (mm)**2'
      read(5,*)mm_cut
      print *, 'mm_cut=', mm_cut
      
      IF(epirea.eq.1)then
         MPI	= MPI0
         mm_exp = MPI0**2
         id(2)	= 14		!Geant ID, proton
         q(2) 	= 1
         pdgid(2)  = 2212		!PDG ID, proton
         if (npart .eq. 4)then
            id(3)     =  7              !Geant ID, pi-zero
            pdgid(3)  = 111		!PDG ID, pi-zero
            q(3)      = 0
         endif
      elseif(epirea.eq.3)then
         MPI	= MPIP
         mm_exp	= mp**2
         id(2)     = 8                !Geant ID, pi-plus
         pdgid(2)  = 211		!PDG ID, pi-plus
         q(2) 	   = 1
         if (npart .eq. 4)then
            id(3)	= 13		!Geant ID, neutron
            pdgid(3)    = 2112		!PDG ID, neutron
            q(3) 	= 0
         endif
      elseif(epirea.eq.5)then
         MPI	= META
         mm_exp = META**2
         id(2)	= 14		!Geant ID, proton
         q(2) 	= 1
         pdgid(2)  = 2212		!PDG ID, proton
         if (npart .eq. 4)then
            id(3)     = 17              !Geant ID, eta-zero
            pdgid(3)  = 221		!PDG ID, pi-zero
            q(3)      = 0
         endif
      else
         go to 5
      endif
      
      if (epirea .eq. 3)filerz(13:13)='p'

c     Set single precision version of pion mass
      mpi_s=mpi

c     Calculate the minimum hadronic mass for pion production:
      wg	= mp + mpi + .0005

      write(6,*)' Input the target thickness (cm)'
      read(5,*)t_targ            
      print *, 't_target=', t_targ

      write(6,*)' Input the target radius, (cm)'
      read(5,*)r_targ
      print *, 'r_target=', r_targ

      write(6,*)' Input the x-coordinate of the beam position, (cm)'
      read(5,*)vertex_x
      write(6,*)' Input the y-coordinate of the beam position, (cm)'
      read(5,*)vertex_y
      write(6,*)' Input the z-coordinate of the beam position, (cm)'
      read(5,*)vz
      print *, 'VX,VY,VZ=', vertex_x,vertex_y,vertex_z

      bfac		= 4./3.
      hydrogen_rad 	= 865 		! hydrogen radiation length (cm)
      t_targ 		= bfac * t_targ / hydrogen_rad 

      write(6,*)' Input the incident electron energy (GeV)'
      read(5,*)ebeam
      print *, 'E beam=', ebeam
      
c     calculate the incident momentum

      es	= ebeam
      ps	= sqrt(es**2-mel**2)
      rs	= ps/es
      
c     cut off q2 at the value for 90 degree elastic scattering

      s		= 0.5
      q2max	= 4.*ebeam**2*s/(1.+2.*ebeam*s/mp)
      
c     Choose two limits for Q**2

      write(6,*)'Input lower and upper limit for Q**2'
      read(5,*)q2_min,q2_max
      print *, ' Q2 min max=', q2_min,q2_max
      
      if (q2_max .gt. q2max) q2_max = q2max
      uq2_min	= 1/q2_max
      uq2_max	= 1/q2_min
      uq2rng	= uq2_max-uq2_min

c     Set the limits on the range of scattered electron energies

      write(6,*)'Input lower and upper limit for scat e- energy(GeV).'
      read(5,*)ep_min,ep_max
      print *, 'Eelectron min max', ep_min,ep_max
      
      epmax	= es - (wg**2 + q2_min - mp**2)/2./mp
      
      if (ep_max .lt. epmax)epmax=ep_max
      eprng	= epmax - ep_min

c     Choose a maximum value for the range of photon energies

      write(6,*)' Input minimum photon energy for integration'
      read(5,*)delta
      print *, ' Min photon energy=', delta

c     Select the number of events desired in the rz file

      write(6,*)' Input the desired number of events in the ntuple'
      read(5,*)nmax
      print *, 'Number of events in the ntuple=', nmax
      
      nprint	= nmax/25

      write(6,*)' Input a multiplication factor for sigr_max'
      read(5,*)fmcall
      print *, ' Factor=', fmcall
      
      if (fmcall .eq. 0.)then
         write(6,*)' Input sigr_max'
         read(5,*)sigr_max
         print *, 'sigr_max=', sigr_max
      endif
c
      if (th_opt .gt. 10) then
         write(6,*)'W mininum input'
         read(5,*) w_min_input
         print *, 'W min=', w_min_input
      endif
c
 1    mcall_max	= 0
      ntold	= 0
      events	= 0

c     Use the internal clock to initialize the random number generator

      call timex(itime1)
      call getunixtime(idum)
      call getasciitime(idum,ctime)

      idum	= -idum
      month 	= ctime(5:7)
      day   	= ctime(9:10)
      year  	= ctime(23:24)
      
      if (day(1:1).eq. ' ')then
        ir1=48
        day(1:1)=char(ir1)
      endif

      call date_and_time(values=timevals)
      call random_seed(size=ssize)
      write(6,*)'seed:',ssize
      sseed(:) = timevals(8)
      call random_seed(put=sseed)

c      call random_number(rnd)
c      write(6,*)'seed:',rnd
c      call random_number(rnd)
c      write(6,*)'seed:',rnd
c      call random_number(rnd)
c      write(6,*)'seed:',rnd
c      write(6,*)'seed:',idum,' from start time ',ctime


      nevent	= 0
      t_elapse	= 0.
      itime2	= itime1
      ntries	= 0
      sig_int	= 0.
      sig_tot	= 0.

c     Name the output rz file according to meson type and beam energy.
c     filerz=aaoradgen-pi0-1.6-0811.rz.0, for example.

      open(unit=12,file=file_out)
      open(unit=13,file=lund_out)

c     Initialize BOS
     
      bosout 	= file_bos
      recname 	= 'MCEVENT'
     
Cvpk      call bos_ini(recname)
  
c     set up the ntuple file
      print *,'     set up the ntuple file'

      lrecl	= 1024
      
      call hlimit(5000000)
      call hropen(1,'aaoradgen',filerz,'n',lrecl,istat)
      call hbookn(10,'func1',n,'aaoradgen',1000,tag)

      open(unit=12,file=file_out)

      write(12,*)' AO Calculation of Single Pion Production'
      write(12,*)' Starting time:', ctime
      write(12,*)' Epirea (1 for pi0, 3 for pi+, 5 for eta) =',epirea
      write(12,*)' Target thickness =',t_targ*3./4.,' (r.l.)'
      write(12,*)' Incident electron energy =',ebeam,' GeV'

      write(12,*)' Electron Q**2 limits:',q2_min,q2_max
      write(12,*)' Lower and upper limit for scattered electron',
     * ' energy(GeV):',ep_min,epmax
      write(12,*)' Minimum photon energy for integration (delta):',delta

c     Use a new variable in place of ek. Let uek=exp(-kek*ek)
c     ek=-(1/kexp)alog(uek).  The factor of 5. was chosen empirically for kexp
c     by looking at the ek spectrum for E0=1.6 GeV.
c     Let uek range from 0 to 1. Then ek will range from 0 and infinity.
c     This requires a jacobian. Jacobian=1./(kexp*uek)=(1./kexp)exp(kexp*ek)

      kexp	= 5.

      if (fmcall .eq. 0.)then
         write(6,*)' sigr_max from input data =',sigr_max
         write(12,*)' sigr_max from input data =',sigr_max
         go to 20
      endif

c     Do a preliminary calculation to estimate the maximum value
c     of the integrand

c     calculate the energy and momentum of the scattered electron,
c     and calculate Q**2 at the delta mass, 1.232 GeV

 10   q2	= q2_min
      q0	= (1.232**2 - mp**2 + q2)/2./mp
      if (th_opt.gt.10) then
        q0	= (w_min_input**2 - mp**2 + q2)/2./mp   
      endif
      ep	= es-q0
      pp	= sqrt(ep**2 - mel**2)
      rp	= pp/ep
      s		= q2/4/es/ep
      th0	= 2.*asin(sqrt(s))
      theta	= th0*180./pi
      T0	= th0
      snt0	= sin(th0)
      cst0	= cos(th0)

c     calculate kinematic quantities needed for the Mo and Tsai calculation

      u0	= es - ep + mp
      pu	= sqrt(ps**2 + pp**2 - 2*ps*pp*cst0)
      uu	= u0**2 - pu**2
      csths	= (ps - pp*cst0)/pu
      csthp	= (ps*cst0 - pp)/pu
      snths	= sqrt(1. - csths**2)
      snthp	= sqrt(1. - csthp**2)
      ts	= acos(csths)
      tp	= acos(csthp)
      qsq	= q2
      sp	= es*ep - ps*pp*cst0

      sigr_max	= 0.
      cstk1	= (es - ep*cst0)/(sqrt(es**2 + ep**2 - 2*es*ep*cst0))
      cstk2	= (es*cst0 - ep)/(sqrt(es**2 + ep**2 - 2*es*ep*cst0))
      csrnge	= csrng
      
      if ((1.-cstk1) .lt. csrnge)       csrnge = 1.-cstk1
      if ((cstk1-cstk2) .lt. 2.*csrnge) csrnge = 0.5*(cstk1-cstk2)
      
      csrngb	= csrng/40.
      
      if (csrngb .gt. csrnge/5.) csrngb = csrnge/5.
      
      delphi	= pi/9.
      phik	= (myran()-0.5)*delphi
      mpfac	= delphi/2./pi

      ek	= delta
      
      jacob	= exp(kexp*ek)/kexp*(q2**2/(2*es*ep))
c
      ehel = 0   
      if(flag_ehel.eq.1) ehel = 1    
c    
      do i=1,10000
        csthcm	= 2.*(myran()-0.5)
        phicm	= 360.*myran()

        csran	= myran()
        
        if (csran .gt. .3) then
            cstk = 2.*csrngb*(myran()-0.5) + cstk1
        else
            cstk = 2.*csrngb*(myran()-0.5) + cstk2
        endif
        
        mcfac	= csrngb /reg1
        phik	= (myran()-0.5)*delphi
        mpfac	= delphi/2./pi

        Tk	= acos(cstk)
        sntk	= sin(Tk)
        sdotk	= es*ek - ps*ek*cstk*csths - ps*ek*sntk*snths*cos(phik)
        pdotk	= ep*ek - pp*ek*cstk*csthp - pp*ek*sntk*snthp*cos(phik)

        sigr	= sigma(ek,Tk,csthcm,phicm,ehel)
        sigr	= sigr*mcfac*mpfac
        sigr	= sigr*jacob      
        if (sigr .gt. sigr_max) then
           sigr_max	= sigr
           ek_max	= ek
           tk_max	= Tk
           csthcm_max	= csthcm
           phicm_max	= phicm
        endif
        
      enddo

      write(6,*)'sigr_max,ek_max,tk_max,csthcm_max,phicm_max',
     * sigr_max,ek_max,tk_max,csthcm_max,phicm_max
      write(12,*)'sigr_max,ek_max,tk_max,csthcm_max,phicm_max',
     * sigr_max,ek_max,tk_max,csthcm_max,phicm_max
     
      sigr_max	= sigr_max*fmcall
      
 25   write(6,*) ' sigr_max changed to ',sigr_max
      write(12,*)' sigr_max changed to ',sigr_max

c   %%%%%%%%%%%%%%%%%%% Main Calculation  %%%%%%%%%%%%%%%%%%%%%%%
c     Use a Monte-Carlo to calculate a distribution of nmax events
c     distributed according to the Mo-Tsai integrand.
      ehel = 0      
 20   continue
      ntries	= ntries+1
      if(flag_ehel.eq.1) ehel = get_spin(idum)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      hel       = ehel
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c     calculate the energy of the electron at the scattering point
c     after making its way through the target.  First, randomly
c     choose the interaction point.
      
c     Change from units of r.l. to cm
      t_targ 	=  hydrogen_rad * t_targ / bfac 
      targs	= t_targ*myran()

c     Change into proper coordinate system
      vertex_z	= vz + targs - t_targ / 2.0
 
c     calculates the distance that the electron stays in the target
      targp	= r_targ/sin(T0)
      temp	= (t_targ - targs)/cos(T0)
      
      if (temp .lt. targp) targp = temp
      
c     Back to units of r.l.      
      t_targ 	=  bfac * t_targ / hydrogen_rad 
      targs  	=  bfac * targs / hydrogen_rad
      targp  	=  bfac * targp / hydrogen_rad      

c     Now calculate the radiation loss

 22   xs	= myran()
      eloss	= xs**(1./targs)
      gxs	= 1.-eloss
      xtest	= myran()
      
      if (xtest.gt.gxs) go to 22
      
      es	= ebeam*(1.-eloss)
      
c     Cut off the incident energy at e_s=ebeam/4.

      if (es .lt. ebeam/4.) go to 20
      
      ps	= sqrt(es**2-mel**2)
      rs	= ps/es

      uq2	= uq2_min+uq2rng*myran()
      q2	= 1./uq2

c     calculate the energy and momentum of the scattered electron,
c     and calculate Q**2

      ep	= epmax-eprng*myran()
      
c     check to see if the scattered electron energy is below the
c     detector threshold

      if (ep .lt. ep_min)go to 20

      q0	= es-ep
      s		= q2/4/es/ep
      
c     cut off scattering at 90 degree

      if (s .gt. .5)go to 20
      
      th0	= 2.*asin(sqrt(s))
      theta	= th0*180./pi
      T0	= th0
      snt0	= sin(th0)
      cst0	= cos(th0)
      
c     check to see if the scattered electron energy is above
c     the pion threshold for this angle.

      ep_test	= (mp**2 + 2*mp*es - wg**2)/2./(mp + 2.*es*s)
      
      if (ep .gt. ep_test)go to 20
      
      pp	= sqrt(ep**2-mel**2)
      rp	= pp/ep
      qsq	= q2
      
      if (qsq .le. 0.)then
         write(6,*)' Main-1:, qsq =',qsq
         go to 20
      endif
      
      qvecx	=     -pp*sin(th0)
      qvecz	= ps - pp*cos(th0)

      w2	= mp**2 + 2*mp*q0 - q2
      
      if (w2 .lt. mp**2)go to 20
      
      epw	= sqrt(w2)
      
      if (th_opt.gt.10) then
        if (epw .lt. W_min_input) go to 20
      else
        if (epw .lt. wg+0.002)go to 20
      endif
      
c     calculate kinematic quantities needed for the Mo and Tsai calculation

      u0	= es - ep + mp
      pu	= ps**2 + pp**2 - 2*ps*pp*cst0
      
      if (pu .le. 0.)then
         write(6,*)' Main-2, pu**2 =',pu
         go to 20
      endif
      
      pu	= sqrt(pu)
      uu	= u0**2 - pu**2
      csths	= (ps - pp*cst0)/pu
      csthp	= (ps*cst0 - pp)/pu
      snths	= 1.-csths**2
      
      if (snths**2 .le. 0.)then
         write(6,*)' Main-3: snths =',snths
         go to 20
      endif
      
      snths	= sqrt(snths)
      snthp	= 1.-csthp**2
      
      if (snthp .le. 0.)then
         write(6,*)' Main-4: snthp**2 =',snthp
         go to 20
      endif
      
      snthp	= sqrt(snthp)
      ts	= acos(csths)
      tp	= acos(csthp)
      sp	= es*ep - ps*pp*cst0

      cstk1	= csths
      cstk2	= csthp
      csrnge 	= csrng
      
      if (cstk1 .lt. cstk2)then
         cstmp=cstk1
         cstk1=cstk2
         cstk2=cstmp
      endif
      
      if ((1.-cstk1) .lt. csrnge)       csrnge=1.-cstk1
      if ((cstk1-cstk2) .lt. 2.*csrnge) csrnge=0.5*(cstk1-cstk2)
      
      csrngb	= csrng/40.
      
      if (csrngb .gt. csrnge/5.) csrngb=csrnge/5.
      
      csran	= myran()
      rn1	= myran()
      
      if (rn1 .gt. 0.5)then
        rn1=1.
      else
        rn1=-1.
      endif
      
      rn2	= myran()
      delphi	= ps*cstk1*csrngb/pp/sqrt(1.-cst0**2)/sqrt(1.-cstk1**2)
      
      if (delphi .lt. pi/9.) delphi=pi/9.
      if (delphi .gt. 2.*pi) delphi=2.*pi
      if (cstk1 .gt. .995)   delphi=2.*pi
      
      delphi	= pi/9.
       
      if (csran .lt. reg1)then
         intreg	= 1
         cstk	= cstk1+(2.*rn2-1.)*csrngb
         mcfac	= csrngb /reg1
         phik	= (myran()-0.5)*delphi
         mpfac	= delphi/2./pi
      elseif(csran .lt. reg2)then
         intreg	= 2
         cstk	= cstk2+(2.*rn2-1.)*csrngb
         mcfac	= csrngb /(reg2-reg1)
         phik	= (myran()-0.5)*delphi
         mpfac	= delphi/2./pi
      elseif(csran .lt. reg3)then
         intreg	= 3
         cstk	= cstk1+rn1*(csrngb+rn2*(csrnge-csrngb))
         mcfac	= (csrnge-csrngb) /(reg3-reg2)
         phik	= (myran()-0.5)*delphi
         mpfac	= delphi/2./pi
      elseif(csran .lt. reg4)then
         intreg	= 4
         cstk	= cstk2+rn1*(csrngb+rn2*(csrnge-csrngb))
         mcfac	= (csrnge-csrngb) /(reg4-reg3)
         phik	= (myran()-0.5)*delphi
         mpfac	= delphi/2./pi
      else
         intreg	= 5
 45      cstk	= 2.*rn2-1.
         phik	= 2.*pi*(myran()-0.5)
         if (abs(cstk-cstk1) .lt. csrnge .or.
     &       abs(cstk-cstk2) .lt. csrnge) then
         if(abs(phik).lt. delphi/2.) go to 45
        endif
           
c     combine mcfac and mpfac into one factor and set mpfac=1

        mcfac=(1.-csrnge*delphi/pi)/(1.-reg4)
        mpfac=1.
      endif

      Tk	= acos(cstk)
      sntk	= sin(Tk)
      
c     change the following on Jan. 19, 1999
c      phran=myran()
c      if (phran .lt. .2)then
c      else
c 48      phik=2*pi*(myran()-0.5)
c         if(abs(phik).lt. pi/180.)go to 48
c         mpfac=1.25
c      endif
c     end of jan 19, 1999 correction

      ekmax	= 0.5*(uu - wg**2)/(u0 - pu*cstk)
      
      if (ekmax .gt. ebeam)then
        write(6,*)' Main-5: ekmax =',ekmax
        ekmax=ebeam
      endif
      
c     choose ek by making a change of variables

 78   uek	= myran()
 
      if (uek .lt. 0.1E-20)then
        write(6,*)' Main-6: uek =',uek
        go to 20
      endif
      
      ek	= -alog(uek)/kexp
      
      if (ek .gt. ekmax)go to 20

      csthcm	= -1.+2.*myran()
      phicm	= 360.*myran()
      
c********************************** Ek < Delta ****************************************      

      if (ek .lt. delta)then
        intreg=6
        
c      print *, th0,qsq,epw,csthcm,phicm
c     calculate the non-radiative cross section

        call dsigma(th0,qsq,epw,csthcm,phicm,th_opt,epirea,res_opt,sigma0
     * ,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)
     
        if (sigma0.le.0.) go to 20
c       PRINT *,'DSIGMA=',
c     *  th0,qsq,epw,csthcm,phicm,th_opt,epirea,res_opt,sigma0
c     * ,sigu,sigt,sigl,sigi,sigip,asym_p,ehel
 
c       print *, 'AAO_RAD s0,u,t,l,i:',sigma0,sigu,sigt,sigl,sigi
     
c      if (abs(epw-1.232).lt.0.010.and.abs(qsq-0.35).lt.0.1) then 
c      print *, qsq,epw,csthcm,phicm
c      print *, 'AO:',sigma0,sigu,sigl,sigt,sigi
c      call dsigma(th0,qsq,epw,csthcm,phicm,4,epirea,res_opt,sigma0
c     * ,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)
c      print *, 'MAID:',sigma0,sigu,sigl,sigt,sigi
c      print *, ''
c     endif

c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     The AO program calculates sigu,sigl,sigt and sigi with
c     the kinematic factors included. To convert them to the
c     form needed for the above problem we must use

        epeps	= 1. + 2.*(1 + q0**2/qsq)*s/(1-s)
      
        if (epeps .le. 1.)then
           write(6,*)' Main-6: epeps-inv =',epeps
           go to 20
        endif
        
        epeps	= 1./epeps
        kfac	= (w2 - mp**2)/2./mp
        
c        cfac=1/Gamma_T
c        cfac=2*pi**2*qsq*(es/ep)*(1-epeps)/kfac/alpha
c        sigu=sigu*cfac
c        sigl=sigl*cfac
c        sigt=sigt*cfac
c        sigi=sigi*cfac

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Mo and Tsai give a formula for the cross section in terms
c     of amplitudes, f and g.  Here we try to extract f and g
c     from the AO cross section.  I may not have this part right.

c     sig0=4.*pi**2*alpha/mp**2

c     The following would work if there were no polarization terms

c     f=2.*kfac/(mp**3)*qsq/(qsq+q0**2)*(sigu+sigl)/sig0
c     g=(2.*kfac/mp)*(sigu)/sig0

c     I don\'t know what to do with the polarization terms: these give
c     the dependence on the center of mass phi decay angle.  I try
c     the following:

        nu	= (w2 + qsq - mp**2)/2/mp
        
        f = (1./(2*pi**2*alpha*mp))*(kfac/(1.+nu**2/qsq))*(sigu+sigl)
        g = mp/(2*pi**2*alpha) * kfac * sigu
        
        fkt = (w2-mp**2+mpi**2)/2./epw
        fkt = sqrt(fkt**2-mpi**2)*2.*epw/(w2-mp**2)
        f   = f*fkt
        g   = g*fkt

c     Calculate the non-radiative cross section

        signr = 2*(alpha*ep/qsq)**2
     &   * (f*mp*cos(th0/2)**2+2*g*sin(th0/2)**2/mp)
     
c        asig = fkt*(sigu+epeps*sigl)/cfac
c        print *, 'Compare x-sections',signr,asig 
     
c     Modulate the cross section with the polarization and interference terms

        signr = signr * (1.+(epeps*sigt*cos(phicm*pi/90.)
     &  + sqrt(epeps*(1.+epeps)/2)*sigi*cos(phicm*pi/180.)
     &  + ehel*sqrt(epeps*(1.-epeps)/2)*sigip*sin(phicm*pi/180.))
     &  /(sigu+epeps*sigl))

        if (signr .le. 0.)then
           write(6,*)'signr,sigu,sigl,sigt,sigi,epeps',
     &     signr,sigu,sigl,sigt,sigi,epeps
           write(6,*)'qsq,epw,th0,csthcm,phicm,epirea',
     &      qsq,epw,th0,csthcm,phicm,epirea
        endif

c     Calculate the radiative correction factor deltar for the cross
c      section.  This includes vertex corrections and the integration
c     up to photon energies of delta.

        x1	= (ep-es)/ep
        x2	= (es-ep)/es
        s1	= spence(x1)
        s2	= spence(x2)
        
        deltar	= -(alpha/pi)*(28./9. -(13./6.)*dlog(2.*sp/mel**2)
     &   -s1-s2)
        delinf	= -(alpha/pi)*dlog(es*ep/delta**2)*(dlog(2.*sp/mel**2)-1.)
      
        sigr1	= signr*(1.+deltar)*exp(delinf)
        
c	calculate average differential cross section in the region
c	from ek=0 to delta, and from cos(thetak)=-1 to 1.

	sigr	= sigr1/delta/4./pi
	
        if (sigr .gt. 0.) then
          go to 28
        else
          go to 20
        endif
        
c       end of section for calculation with ek < delta.
c       Normally, go to statement 28

      endif

      sdotk = es*ek - ps*ek*cstk*csths - ps*ek*sntk*snths*cos(phik)
      pdotk = ep*ek - pp*ek*cstk*csthp - pp*ek*sntk*snthp*cos(phik)

      sigr	= sigma(ek,Tk,csthcm,phicm,ehel)
      
      if (sigr .le. 0.) go to 20
      
28    jacob	= exp(kexp*ek)/kexp/(2*es*ep)*q2**2
      sigr	= sigr*jacob
      
      call missm(epirea,ebeam,es,ep,th0,ek,cstk,phik,mpi_s,
     *  ppx,ppy,ppz,eprot,ppix,ppiy,ppiz,epi,ekx,eky,ekz,
     *  csthcm,phicm,wreal,mm2)
      
      if (abs(mm2-mm_exp) .gt. mm_cut)go to 20

c     Compare sigr to the sigr_max to determine whether to generate
c     an event.

      sigr	= mcfac*mpfac*sigr
      sig_ratio	= sigr/sigr_max
      sig_tot	= sig_tot + sigr
      
c     Choose the number of times, mcall, to call the routine used
c     to calculate kinematic quantities for the n-tuple.

      rtest	= myran()
      mcall	= sig_ratio
      stest	= sig_ratio - mcall
      
      if (stest .gt. rtest)mcall = mcall + 1
       
      if (mcall .gt. mcall_max) mcall_max=mcall
      if (mcall .gt. 10)then
c        write(6,*)' mcall =',mcall,' intreg=',intreg
c        write(6,*)es,ep,th0,ek
c        write(6,*)'cstk1,cstk2,cstk,phik',cstk1,cstk2,cstk,phik
c        write(6,*)' csrnge=',csrnge
c        write(6,*) ppx,ppy,ppz,eprot,ekx,eky,ekz
c        write(6,*)csthcm,phicm,wreal,mm2

         csdotk	= cstk*csths+sntk*snths*cos(phik)
         cpdotk	= cstk*csthp+sntk*snthp*cos(phik)
         cqdotk	= (ps*csdotk-pp*cpdotk)/sqrt(ps**2+pp**2-2*ps*pp*cos(th0))
         
c         write(6,*)'es,ep,ek',es,ep,ek
c         write(6,*)'cstk,phik,csdotk,cpdotk'
c     +   ,cstk,phik,csdotk,cpdotk
c         write(6,*)'cqdotk',cqdotk

         write(12,*)' mcall =',mcall,' intreg=',intreg
         write(12,*)'es,ep,ek',es,ep,ek
         write(12,*)'cstk,phik,csdotk,cpdotk'
     +   ,cstk,phik,csdotk,cpdotk
         write(12,*)'cqdotk',cqdotk
      endif
      
c     If mcall .gt. 0 generate mcall n-tuple events

      if (mcall .eq. 0)go to 30

      ep_sav	= ep

      if (mcall .lt. 100)then
         dismc(intreg,mcall)=dismc(intreg,mcall)+1
      else
         dismc(intreg,100)=dismc(intreg,100)+1
      endif

      do j=1,mcall

c     Calculate the radiation loss for the electron leaving the target

 222     xs	= myran()
         eloss	= xs**(1./targp)
         gxs	= 1.-eloss
         xtest	= myran()
         
         if (xtest.gt.gxs)go to 222
         
         ep	= ep_sav*(1.-eloss)

c     correct the following section on Jan. 23, 1999

         if (ep .lt. ep_min)then
            sig_tot = sig_tot-sigr_max
            go to 24
         endif
         
c     end of correction

         call missm(epirea,ebeam,es,ep,th0,ek,cstk,phik,mpi_s,
     *  ppx,ppy,ppz,eprot,ppix,ppiy,ppiz,epi,ekx,eky,ekz,
     *  csthcm,phicm,wreal,mm2)

         if(mm2 .eq. 0. .or. ep .lt. ep_min)go to 24

         w2	= mp**2+2*mp*(ebeam-ep)-2*es*ep*(1-cos(th0))+2.*mel**2
         epw	= sqrt(w2)

c     Calculate the members of the n-tuple and ouput it to the rz file.

         ntp(1)	= es
         ntp(2)	= ep
         ntp(3) = theta
         ntp(4) = epw
         ntp(5) = wreal
         ntp(6) = ppx
         ntp(7) = ppy
         ntp(8) = ppz
         ntp(9) = eprot
         ntp(10)= ppix
         ntp(11)= ppiy
         ntp(12)= ppiz
         ntp(13)= epi
         ntp(14)= csthcm
         ntp(15)= phicm
         ntp(16)= mm2
         ntp(17)= ek
         ntp(18)= cstk
         ntp(19)= phik*180./pi
         ntp(20)= qvecx
         ntp(21)= qvecz
         ntp(22)= q0
         ntp(23)= cst0
         ntp(24)= ekx
         ntp(25)= eky
         ntp(26)= ekz
         ntp(27)= vertex_x
         ntp(28)= vertex_y
         ntp(29)= vertex_z
         ntp(30)= qsq                  
         ntp(31)= ehel
         ntp(32)= asym_p

         call hfn(10,ntp)
c         nevent	= nevent+1
        
         do jj = 1,npart
           v(jj,1) = vertex_x
           v(jj,2) = vertex_y
           v(jj,3) = vertex_z
         enddo
         
c     rotate all the momenentum by a random angle around the beam line

         phir	= 2.*pi*myran()
         rotc	= cos(phir)
         rots	= sin(phir)
          
c     momentum of scattered electron:

         px 	=  ep*sin(theta*pi/180.)
         py 	=  0.
         p(1,1) = px*rotc+py*rots
         p(1,2) = py*rotc-px*rots
         p(1,3) = ep*cos(theta*pi/180.)
         p(1,4) = ep
         
c     momentum of charged hadron

         if ((epirea .eq. 1).or.(epirea.eq.5)) then
             p(2,1) = ppx*rotc+ppy*rots
             p(2,2) = ppy*rotc-ppx*rots
             p(2,3) = ppz
             p(2,4) = eprot
             
c     momentum of scattered pion

             if (npart .eq. 4)then
                p(3,1) = ppix*rotc+ppiy*rots
                p(3,2) = ppiy*rotc-ppix*rots
                p(3,3) = ppiz
                p(3,4) = epi
c                print *,sqrt(p(3,4)**2-p(3,3)**2-p(3,2)**2-p(3,1)**2)
             endif
             
         else
         
c     momentum of scattered pion

             p(2,1) = ppix*rotc+ppiy*rots
             p(2,2) = ppiy*rotc-ppix*rots
             p(2,3) = ppiz
             p(2,4) = epi
             
             if (npart .eq. 4)then
                p(3,1) = ppx*rotc+ppy*rots
                p(3,2) = ppy*rotc-ppx*rots
                p(3,3) = ppz
                p(3,4) = eprot
             endif
             
         endif

c     momentum of radiated photon

         if (npart .eq. 4)then
         
c     suppress radiated photons when E-gamma is less than delta

           if (ek .le. delta)then
             p(4,1)=0.
             p(4,2)=0.
             p(4,3)=1.e-5
             p(4,4)=1.e-5
           else
             p(4,1) = ekx*rotc+eky*rots
             p(4,2) = eky*rotc-ekx*rots
             p(4,3) = ekz
             p(4,4) = ek
             endif
         endif
*
*   Cut of theta>16 deg applied. For e1-dvcs eta meson production.
*   12/12/2011, Ivan Bedlinskiy
*         
c         if (theta.lt.16) goto 24

         lundmass(1) = MEL
         lundmass(4) = 0
         if ((epirea .eq. 1).or.(epirea.eq.5)) then
            lundmass(2) = MP
            lundmass(3) = MPI
         else if (epirea.eq.3) then
            lundmass(2) = MPIP
* neutron mass 940 MeV
            lundmass(3) = 0.940 
         endif

         write(13,*) "4 1 1 0 0 ",q2/(2.0*MP*(ebeam-p(1,4))), 0, epw, 
     *   q2, (ebeam-p(1,4))
         write(13,*) 1, q(1), 1, pdgid(1), 0, 0,
     *               p(1,1), p(1,2), p(1,3), p(1,4),
     *               lundmass(1), vertex_x, vertex_y, vertex_z
         write(13,*) 2, q(2), 1, pdgid(2), 0, 0,
     *               p(2,1), p(2,2), p(2,3), p(2,4),
     *               lundmass(2), vertex_x, vertex_y, vertex_z
         write(13,*) 3, q(3), 1, pdgid(3), 0, 0,
     *               p(3,1), p(3,2), p(3,3), p(3,4),
     *               lundmass(3), vertex_x, vertex_y, vertex_z
         write(13,*) 4, q(4), 1, pdgid(4), 0, 0,
     *               p(4,1), p(4,2), p(4,3), p(4,4),
     *               lundmass(4), vertex_x, vertex_y, vertex_z

C         print *,'Event to be written'
Cvpk         call bos_out	! Pack the BOS banks and write out to file   
         nevent	= nevent+1
C         print *,'Events=',nevent,nprint
         goto 124
 24      continue
       enddo

c     Talk to the user every now and then
 124    continue
       ntell	= nevent/nprint-ntold
       
       if (ntell .gt. 0)then
          write(6,*)' ntries, nevent, mcall_max: '
     *    ,ntries,nevent,mcall_max
          write(12,*)' ntries, nevent, mcall_max: '
     *    ,ntries,nevent,mcall_max
     
         events	= nevent
         tries	= ntries
         sig_int= events/tries
         
c     photon phase space = 4*pi*delta-omega=4*pi (after change of variables)
c     hadron phase space = 4*pi
c     electron phase space =2*pi*ucrng*eprng

         sig_int = sig_int*sigr_max*(4.*pi)**2*(2.*pi*uq2rng*eprng)
         sig_sum = sig_tot*(4.*pi)**2*(2.*pi*uq2rng*eprng)/tries

         write(6,*)' Integrated cross section (MC, numerical) ='
     +   ,sig_int,sig_sum, ' mu-barns'
         write(6,*)' Beam time at Lum=1.0E34 =',events/sig_sum*1.E-4
     *    ,' seconds'
         write(12,*)' Integrated cross section =',sig_int,sig_sum
     *    ,' micro-barns'
         write(12,*)' Beam time at Lum=1.0E34 =',events/sig_sum*1.E-4
     *    ,' seconds'
         write(12,*)' Luminosity =',events/sig_sum
     *    ,' inverse micro-barns'

          call timex(itime2)
          
          WRITE(12,*)' t_elapse,itime1,iteme2=',
     *                t_elapse,itime1,itime2          
          WRITE(6,*)' t_elapse,itime1,iteme2=',
     *                t_elapse,itime1,itime2          
C          t_elapse	= t_elapse+INT(itime2-itime1)
          t_elapse	= t_elapse+itime2-itime1
          itime1	= itime2
          write(6,*)' Elapsed CPU time = ',t_elapse/60,' minutes'
          write(12,*)' Elapsed CPU time = ',t_elapse/60,' minutes'
          
          ntold	= ntold+1
          
       endif

c     Do we have enough events in the n-tuple?

       if (nevent .gt. nmax)go to 50
       
 30    go to 20

c     Close out the n-tuple file

 50   call hrout(0,icycle,' ')
      call hrend('aaoradgen')


      close(12)
      close(13)

      open(unit=14,file=file_sum)
      
      write(14,*)' AO Calculation of Single Pion Production'
      write(14,*)' AO Calculation of Single Pion Production'

      write(14,*)' Starting time: ',ctime
      write(14,*)' Epirea (1 for pi0, 3 for pi+) =',epirea
      write(14,*)'Target thickness =',t_targ*3./4.,' (r.l.)'
      write(14,*)' Incident electron energy =',ebeam,' GeV'

      write(14,*)'Electron Q**2 limits:',q2_min,q2_max
      write(14,*)'lower and upper limit for scattered electron',
     * ' energy(GeV):',ep_min,epmax
      write(14,*)' Minimum photon energy for integration (delta):',delta
      write(14,*)'sigr_max,ek_max,tk_max,csthcm_max,phicm_max',
     * sigr_max,ek_max,tk_max,csthcm_max,phicm_max
      write(14,*)' ntries, nevent, mcall_max: '
     *    ,ntries,nevent,mcall_max
      write(14,*)' Missing-mass squared cut at:',mm_cut
      write(14,*)' Integrated cross section =',sig_int,sig_sum
     *    ,' micro-barns'
      write(14,*)' Beam time at Lum=1.0E34 =',events/sig_sum*1.E-4
     *    ,' seconds'
       write(14,*)' Luminosity =', nevent/sig_sum
     *    ,' inverse micro-barns'
      write(14,*)' Elapsed CPU time = ',t_elapse/60,' minutes'
      write(14,*)' CPU time/event = ', t_elapse/nevent,' sec'

      write(14,1425)reg1,reg2-reg1,reg3-reg2, reg4-reg3,1.-reg4
      write(14,*)'csrng =',csrng
 1425 format(' size of cosine regions:',5f5.2)
      do intreg=1,6
        write(14,1410)intreg,(dismc(intreg,j),j=1,100)
      enddo

 1410 format(' Distribution of mcall values, region', i2 /10(1x,10i7/))

      close(14)
Cvpk      CALL bos_end(recname)

 99   continue

       write(6,*)' END OF LOOP'
      stop
      end

      real function sigma(ek,Tk,epcos,epphi,ehel)

c     Calculate the single pion electroproduction cross section with
c      radiative tail, according to the prescription in Mo and Tsai.
c     Mo and Tsai  calculate dsigma/d_omega dp for (omega>delta)
c     omega is the energy of the radiated photon.

c     The Mo and Tsai cross section for the 3-3 resonance is replaced with
c     the AO cross section for single pion production from the proton.


      implicit none
      COMMON/ALPHA/ ALPHA,PI,MP,MPI,MEL,WG,EPIREA,TH_OPT,RES_OPT
      common /radcal/T0,es,ep,ps,pp,rs,rp,u0,pu,uu,cst0,snt0,csths,csthp
     * ,snths,snthp,pdotk,sdotk

      real*8 ek,Tk
      real*8 alpha,pi,mp,mpi,mel,wg
      real*8 T0,es,ep,ps,pp,rs,rp,u0,pu,uu,pdotk,sdotk
      real*8 cst0,snt0,csths,csthp,snths,snthp
      real*8 csthk,snthk,qq,mf2
      real*8 sp
      real*4 ffac1,ffac2,ffac3,ffac4,ffac5,ffac6,ffac
      real*4 gfac1,gfac2,gfac3,gfac4,gfac

      real*4 f,g,fkt
      real*4 sig_r,sigf
      real nu
c     arguments for aaosub_1:
      real*4 qsq,epw,th0,sigma0,sigu,sigt,sigl,sigi
      real*4 sigip, asym_p
      real epcos,epphi
      real*4 q0,kfac
      real s2,epeps
      integer epirea, th_opt, res_opt
      integer ehel
c
      th0=T0
      csthk=cos(Tk)
      snthk=sin(Tk)
      qq=2*mel**2-2*es*ep+2*ps*pp*cst0-2*ek*(es-ep)+2*ek*pu*csthk
      mf2=uu-2*ek*(u0-pu*csthk)
      if (mf2 .lt. wg**2 .or. qq .ge. 0.)then
         sigma=0.
         return
      endif

      epw=sqrt(mf2)

      sp=es*ep-ps*pp*cst0

      ffac1=-(mel/pdotk)**2*(2.*es*(ep+ek)+qq/2)
      ffac2=-(mel/sdotk)**2*(2.*ep*(es-ek)+qq/2)
      ffac3=-2.
      ffac4=2/sdotk/pdotk*(mel**2*(sp-ek**2)+sp*(2*es*ep-sp+ek*(es-ep)))       
      ffac5=(2*(es*ep+es*ek+ep*ep)+qq/2-sp-mel**2)/pdotk
      ffac6=-(2*(es*ep-ep*ek+es*es)+qq/2-sp-mel**2)/sdotk
      ffac=ffac1+ffac2+ffac3+ffac4+ffac5+ffac6
      if (ffac .le. 0.)then
         sigma=0.1e-30
         return
      endif

      gfac1=mel**2*(2*mel**2+qq)* (1./(pdotk**2)+1./(sdotk**2) )
      gfac2=4.
      gfac3=4.*sp*(sp-2*mel**2)/pdotk/sdotk
      gfac4=(2*sp+2*mel**2-qq)*(1./pdotk-1./sdotk)

      gfac=gfac1+gfac2+gfac3+gfac4
      if (gfac .le. 0.)then
         sigma=0.1e-30
         return
      endif
      qsq=-qq
      q0=es-ep

c     Use the AO cross section to fake up the values of f and g
c     needed for the Mo and Tsai integration.
c     This should be compared to the formulas given in Mo and Tsai
c     for the 3-3 resonance.  We ought to agree pretty well in the
c     3-3 resonance region.

     
      call dsigma(th0,qsq,epw,epcos,epphi,th_opt,epirea,res_opt,sigma0
     * ,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)
c      print *, 'AO:',sigma0,sigu,sigt,sigl,sigi
c      call dsigma(th0,qsq,epw,epcos,epphi,4,epirea,res_opt,sigma0
c     * ,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)
c      print *, 'MAID:',sigma0,sigu,sigt,sigl,sigi
c      print *, ''
     
      if (sigma0 .le. 0.)then
         sigma = 0.0
         return
      endif

c     According to my calculations
c     F=(2K/Mp**3)[q2/(q2+nu**2)][(sigt+sigl)/sig0]
c     G=(2K/Mp)(sigt/sig0)
c     where K=(W**2-Mp**2)/(2Mp)
c     sig0=4*pi**2*alpha/mp**2=127 microbarns
C     In terms of the sigt and sigl, the cross section,
c     differential in omega_e and E_e, is
c     sigma=[alpha/(4 pi**2][2 * K * E_s]/[q2 * E_s * (1-epsilon)]
c     *  (sigt+epsilon*sigl)
c     where
c     1/epsilon = 1+2tan(theta_e/2)**2(1+nu**2/q2)

c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     The AO program calculates sigu,sigl,sigt and sigi with
c     the kinematic factors included. To convert them to the
c     form needed for the above problem we must use
      s2=(1.-cst0)/2.
      epeps=1+2*(1+q0**2/qsq)*s2/(1-s2)

      if (epeps .gt. 1.)go to 120
        write(6,*)' sigma-1: epeps-inv =',epeps
        sigma=0.1e-30
        return
 120  epeps=1./epeps

      Kfac=(mf2-mp**2)/2./mp
c      cfac=2*pi**2*qsq*es*(1-epeps)/kfac/ep/alpha
c      sigu=sigu*cfac
c      sigl=sigl*cfac
c      sigt=sigt*cfac
c      sigi=sigi*cfac
c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c      sig0=4.*pi**2*alpha/mp**2
c     The following would work if there were no polarization terms
c      f=2.*kfac/(mp**3)*qsq/(qsq+q0**2)*(sigu+sigl)/sig0
c      g=(2.*kfac/mp)*(sigu)/sig0
c     I don't know what to do with the polarization terms: these give
c     the dependence on the center of mass phi decay angle.  I try
c     the following and modulate the cross section later:
      nu=(mf2-mp**2+qsq)/2./mp
      if (nu .gt. 0.)go to 121
       sigma=0.1e-30
       write(6,*)' sigma-2: nu = ',nu
       return

 121  f=kfac/(2.*pi**2*alpha*mp)/(1.+nu**2/qsq)*(sigu+sigl)
      g=mp/(2.*pi**2*alpha)*kfac*sigu
        fkt=(epw**2-mp**2+mpi**2)/2./epw
        fkt=sqrt(fkt**2-mpi**2)*2.*epw/(epw**2-mp**2)
        f=f*fkt
        g=g*fkt

c      f_old=2.*kfac/(mp**3)*qsq/(qsq+q0**2)*(sigu+sigl)/sig0
c      g_old=(2.*kfac/mp)*(sigu)/sig0

c     ??????????????????????????
c     The following formula is the same as B.8 in Mo and Tsai,
c     except I have divided by 2pi.  It seems to me that the
c     result doesn't make sense otherwise.
c     ?????????????????????????


      sig_r=((alpha**3/(2*pi*qq)**2)/mp)*(ep/es)*ek
      sigf=mp**2*f*ffac+g*gfac
      if (sigf .gt. 0.)go to 122
c        print *, 'SIGMA:',th0,qsq,epw,epcos,epphi,theory_opt
c        write(6,*)' sigma-3: ffac, gfac, f, g,sigu,sigl, sigf ='
c     *  ,ffac,gfac,f,g,sigu,sigl,sigf       
        sigma=0.1e-30
        return

c     Modulate the cross section with the polarization and interference
c     terms
 122   if (sigu+epeps*sigl .gt. 0.)go to 123
         write(6,*)' sigma-4:, sigu+epeps*sigl =',sigu+epeps*sigl
         sigma=0.1e-30
         return

 123   sigf=sigf * (1.+(epeps*sigt*cos(epphi*pi/90.)
     &  +sqrt(epeps*(1.+epeps)/2)*sigi*cos(epphi*pi/180.)
     &  +ehel*sqrt(epeps*(1.-epeps)/2)*sigip*sin(epphi*pi/180.))
c     &  /(sigu+epeps*sigl))
     &  /sigf)

       if (sigf .gt. 0.)go to 124
         write(6,*)' sigma-5: sigf =',sigf
         
         sigma=0.1e-30
         return

 124   sigma=sig_r*sigf
       
      return
      end


      function myran()
c     Random number generator used because I can't find one in the
c     library.

      implicit none
      real myran

      call random_number(myran)

      return
      end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine missm(type,ebeam,es,ep,th0,ephot,cstk,phik,mpi,ppx,
     * ppy,ppz,eprot,ppix,ppiy,ppiz,epi,ephotx,ephoty,ephotz,
     * csthcm,phicm,wreal,mm2)

c
c     Input:
c         type= 1- pi0 , 3 - p+, 5 -eta
c         ebeam = incident electron beam energy
c         es  = incident electron energy at interaction point
c         ep = scattered electron energy
c         th0 = electron scattering angle
c         ephot = energy of radiated photon
c         cstk = cosine of the photon angle (relative to the q vector)
c         phik = azimuthal angle of photon
c         csthcm = cosine of hadronic decay angle in the hadronic frame
c         phicm = phi angle in the hadronic frame
c     Output:
c         ppx, ppy, ppz = proton momentum components
c         eprot = proton energy
c         ppix,ppiy,ppiz = pion momentum components
c         epi = pion energy
c         wreal = true hadronic invariant mass
c         mm2 = experimental missing mass

         
c     Choose the hadronic decay angles randomly and calculate the
c     missing mass, the proton momenta and pion momenta.

      implicit none
      common /random/idum

      real*8 es,ep,ps,pp,ephot

      real beta
      real cstk
      real csthe
      real cstq
      real csthcm
      real cxx,cxy,cxz,cyx,cyy,cyz,czx,czy,czz
      real csphk
      real csphi
      real ebeam
      real ephotx,ephoty,ephotz,eph_dot_q
      real epcm
      real epi
      real epicm
      real eprot
      real ewreal
      real gamma
      real mel
      real mp
      real mpi
      real mm2
      real nu
      real pfac
      real phicm
      real phik
      real pi
      real ppix,ppiy,ppiz
      real ppx,ppy,ppz
      real ppiwx,ppiwy,ppiwz
      real ppwx,ppwy,ppwz
      real pstar
      real pwrx,pwry,pwrz,pwr
      real q2,qvec
      real qx,qz
      real q_dot_pp
      real snphk
      real snphi
      real snthcm
      real snthe
      real sntk
      real sntq
      real th0
      real w2
      real wmin
      real wreal
      integer type
      integer*4 idum
c

      mp=.938
      mel=0.511E-3
      pi=3.14159
      wmin=mp+mpi
      csthe=cos(th0)
      snthe=sin(th0)
      nu=es-ep
      ps=abs(es**2-mel**2)
      pp=abs(ep**2-mel**2)
      ps=sqrt(ps)
      pp=sqrt(pp)
      q2=2.*es*ep-2.*ps*pp*csthe-2.*mel**2
      w2=mp**2-q2+2.*mp*nu
c     get components of the q vector
      qx=-pp*snthe
      qz=ps-pp*csthe
      qvec=sqrt(qx**2+qz**2)
      q2=2.*es*ep-2.*ps*pp*csthe-2.*mel**2
      w2=mp**2-q2+2.*mp*nu
c     get components of the q vector
      qx=-pp*snthe
      qz=ps-pp*csthe
      qvec=sqrt(qx**2+qz**2)
c     get components of the photon vector
      if (abs(cstk) .gt. 1.)then
         write(6,*)' missm-1: cstk =',cstk
         cstk=cstk/abs(cstk)
      endif
      sntk=sqrt(1.-cstk**2)
      csphk=cos(phik)
      snphk=sin(phik)

      cstq=qz/qvec
      sntq=sqrt(1.-cstq**2)
      ephotx=ephot*(sntk*csphk*cstq-cstk*sntq)
      ephoty=ephot*sntk*snphk
      ephotz=ephot*(cstk*cstq+sntk*csphk*sntq)

c     calculate the dot product of the photon vector and the q-vector
      eph_dot_q=ephotx*qx+ephotz*qz
c     calculate the mass of the actual hadronic system for the two
c     photon directions.

      wreal=(w2-2.*ephot*(nu+mp)+2.*eph_dot_q)
      if (wreal .le. wmin**2)go to 30
      wreal=sqrt(wreal)
c     go to the end of the loop if the hadronic mass is below the pion
c     threshold.

c     calculate the energy of the actual hadronic system, with radiation
      ewreal=nu+mp-ephot


c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Calculate the momentum components of the nucleon and pion in the
c     Lab frame

      snthcm=(1.-csthcm**2)
      if (snthcm .le. 0.)then
         write(6,*)' missm-2: snthcm =',snthcm
         snthcm=0.0000001
      endif
      snthcm=sqrt(snthcm)
      

c     Calculate cos and sin of phicm
      csphi=cos(phicm*pi/180.)
      snphi=sin(phicm*pi/180.)

c     calculate laboratory components of the resonance momentum vector
         pwrz=qz-ephotz
         pwrx=qx-ephotx
         pwry=-ephoty
         pwr=sqrt(pwrz**2+pwrx**2+pwry**2)
c     calculate the relativistic factors, gamma and beta, for the resonance
         beta=pwr/ewreal
         gamma=ewreal/wreal

c     define angle cosines for the laboratory system with the resonance
c     as the z axis.  Choose the y axis of the resonance frame perpendicular
c     to the laboratory x axis.
         pfac=sqrt(pwry**2+pwrz**2)
         cxx=pfac/pwr
         cxy=-pwrx*pwry/pfac/pwr
         cxz=-pwrx*pwrz/pfac/pwr
         cyx=0
         cyy=pwrz/pfac
         cyz=-pwry/pfac
         czx=pwrx/pwr
         czy=pwry/pwr
         czz=pwrz/pwr

c     calculate the momentum of the pion and proton in the resonance frame
         pstar=( (wreal**2-mp**2-mpi**2)**2/4.-(mp*mpi)**2 )/wreal**2
         if (pstar .le. 0.)then
            write(6,*)'missm-3: pstar=',pstar
            pstar=0.000001
         endif
         pstar=sqrt(pstar)
c     Calculate the energy of the proton and pion in the resonance center
c     of mass frame.
         epcm  = sqrt(pstar**2+mp**2)
         epicm = sqrt(pstar**2+mpi**2)

c     Calculate the pion momentum components and energy in the lab frame
c     where the z axis is the direction of the momentum of the resonance.
c     The center of mass angles are pion angles.
         ppiwx=pstar*snthcm*csphi
         ppiwy=pstar*snthcm*snphi
         ppiwz=gamma*(pstar*csthcm+beta*epicm)
         epi=gamma*(epicm+beta*pstar*csthcm)

c         print *,'M(eta)=',sqrt(epi**2-ppiwz**2-ppiwy**2-ppiwx**2)
c     Calculate momentum components of the pion in
c     the lab frame where the z axis is along the incident beam
         ppix=ppiwx*cxx+ppiwy*cyx+ppiwz*czx
         ppiy=ppiwx*cxy+ppiwy*cyy+ppiwz*czy
         ppiz=ppiwx*cxz+ppiwy*cyz+ppiwz*czz

c     Calculate the proton energy and momentum components
c     in the lab resonance frame.
c     The center of mass angles are pion angles.
         ppwx=-ppiwx
         ppwy=-ppiwy

c         ppwz=-ppiwz+gamma*beta*(epicm+epcm)
c     epicm+epcm=wreal
         ppwz=gamma*beta*wreal-ppiwz

         eprot=gamma*wreal-epi

c     Rotate the lab momentum components of the proton into the frame
c     where the z axis is along the incident electron beam.
         ppx=ppwx*cxx+ppwy*cyx+ppwz*czx
         ppy=ppwx*cxy+ppwy*cyy+ppwz*czy
         ppz=ppwx*cxz+ppwy*cyz+ppwz*czz


c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c     Calculate the square of the missing mass, associated with the
c     momentum components of the charged hadron:
      q2=2.*ebeam*ep*(1.-csthe)
      qz=ebeam-pp*csthe
      qvec=sqrt(qx**2+qz**2)
      nu=ebeam-ep
         if ((type .eq. 1).or.(type.eq.5)) then
         q_dot_pp=qx*ppx+qz*ppz
         mm2=-q2+2*mp**2+2*mp*(nu-eprot)-2*nu*eprot+2*q_dot_pp
         else       
         q_dot_pp=qx*ppix+qz*ppiz
         mm2=-q2+mp**2+mpi**2+2*mp*(nu-epi)-2*nu*epi+2*q_dot_pp
         endif

         return
 30      mm2=0.
         return
         end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      real function spence(x)
c     Calculate the Spence function needed for the Mo and Tsai formula.

      implicit none
      real x
      real pi
      real sintp,sintn
      pi=3.14159

      if (abs(x) .lt. 0.1)then
         spence=x+x**2/4.
c         write(6,*)' spence: abs(x) .lt. 0.1'
         return
      endif

      if (x .gt. 0.99 .and. x .lt. 1.01)then
         spence=pi**2/6.
c         write(6,*)' spence: x=1.'
         return
      endif

      if (x .gt. -1.01 .and. x .lt. -0.99)then
         spence=-pi**2/12.
c         write(6,*)' spence: x= -1.'
         return
      endif

      if (x .gt. 0.)then
         spence=.1025+sintp(x)
c         write(6,*)' x .gt. 0.'
         return
      endif
      spence=-0.0975+sintn(x)
c      write(6,*)' spence: x .lt. 0.'
      return
      end

c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      real function sintp(x)
      implicit none
      real x
      real xstep,sum,y,arg
      integer i

      xstep=(x-.1)/100.
      sum=0.
      y=.1-xstep/2.
      do i=1,100
        y=y+xstep
        arg=abs(1.-y)
        sum=sum-alog(arg)/y
      enddo
      sintp=sum*xstep
      return
      end

c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      real function sintn(x)
      implicit none
      real x,xa,ystep,y,sum
      integer i

      xa=abs(x)
      ystep=(xa-0.1)/100.
      sum=0.
      y=.1-ystep/2.
      do i=1,100
        y=y+ystep
        sum=sum-alog(1.+y)/y
      enddo
      sintn=sum*ystep
      return
      end

c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine gauss(x,y,sigma_x,sigma_y)

c     calculate two random numbers, x, y,  for gaussian distributions
c     with s.d. of sigma_x and sigma_y.

      implicit none

      common /random/idum

      real x,y,sigma_x,sigma_y
      real r1,r2,pi
      real myran

      integer*4 idum

         pi=3.14159
         r1=myran()
         r2=myran()
         r1=sqrt(-2.*alog(r1))
         r2=2.*pi*r2
         x=sigma_x*r1*cos(r2)
         y=sigma_y*r1*sin(r2)
      return
      end
C======================================================================
      function get_spin(iseed)
C----------------------------------------------------------------------
C-
C-   Purpose : Get spin (1 or -1)
C-
C-   Inputs  : random seed
C-
C-   Outputs : get_spin
C----------------------------------------------------------------------
      implicit none
	  integer get_spin
	  integer iseed
      real    random, myran 

	  random = myran()
      random = 0.5 - myran()
      get_spin = 1
      if(random.lt.0) get_spin = -1

	  return
	  end
C======================================================================
