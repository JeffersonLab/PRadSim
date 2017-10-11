! Subroutines copied from ELRADGEN 1.0

!===============================================================================
      subroutine elrad_init(elab, vmin)
     &bind(C, name = "elrad_init")
      use, intrinsic :: ISO_C_BINDING
!-------------------------------------------------------------------------------
      implicit none
      real(C_DOUBLE), intent(IN), VALUE :: elab, vmin
      include 'elrad_const.inc'
      include 'elrad_par.inc'

      aml = 0.5109989461d-3
      aml2 = aml * aml
      amp = 938.272046d-3
      amp2 = amp * amp
      alpha = 0.72973525664d-2
      pi = 3.14159265358979323846d0

      barn = 1d0
      
      ss = 2d0 * elab * amp
      vminn = vmin

      call grid_init

      end

!===============================================================================
      real(C_DOUBLE) function elrad_sigfs(tmin, tmax, q2)
     &bind(C, name = "elrad_sigfs")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      real(C_DOUBLE), intent(IN), VALUE :: tmin, tmax, q2
      real*8 fsirs, xs, ttmin, ttmax
      include 'elrad_par.inc'

      q22 = q2

      xs = fsirs(tmin + 1d-8)
      call simpsx(tmin + 1d-8, tmax - 1d-8, 10000, 1d-4, fsirs, xs)
      elrad_sigfs = xs

      end

!===============================================================================
      real*8 function fsirs(t)
c finite part of radiative cross section integrated over 
c v and phik ds/dQ2/dphi/dt
      implicit none
      real*8 t, fsigt
      include 'elrad_par.inc'

      fsirs = fsigt(ss, q22, 0d0, 0d0, vminn, t, 2)
      end

!===============================================================================
      real(C_DOUBLE) function elrad_sigfh(tmin, tmax, q2, vcut)
     &bind(C, name = "elrad_sigfh")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      real(C_DOUBLE), intent(IN), VALUE :: tmin, tmax, q2, vcut
      real*8 fsirh, xs
      include 'elrad_par.inc'

      q22 = q2
      vmaxx = vcut
      
      xs = fsirh(tmin)
      call simpsx(tmin, tmax, 10000, 1d-4, fsirh, xs)
      elrad_sigfh = xs

      end
   
!===============================================================================
      real*8 function fsirh(t)
c radiative cross section integrated over phik ds/dQ2/dphi/dt/dv
c a helper function to use the simpson integration subroutine
c kinematics t, v are shared through a common block
      implicit none
      real*8 t, fsigt
      include 'elrad_par.inc'

      fsirh = fsigt(ss, q22, 0d0, vminn, vmaxx, t, 0)

      end

!===============================================================================
      subroutine elrad_gentvp(s, q2, vmin, vmax, tgen, vgen, pgen)
     &bind(C, name = "elrad_gentvp")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      real(C_DOUBLE), intent(IN), VALUE :: s, q2, vmin, vmax
      real(C_DOUBLE), intent(OUT) :: tgen, vgen, pgen
      real*8 fsigt, fsigtv, urand
      real*8 sig, sigrand
      real*8 t1, t2, ttn, tto, sin, sio, epsq2
      real*8 v1, vp, vvn, vvo
      real*8 phn, pho
      integer it1, it2, itt, ivv, iph
      include 'elrad_const.inc'
      include 'elrad_grid.inc'
      include 'elrad_par.inc'

      epsq2 = 1d-13

      t1 = (2d0 * amp2 * q2 + vmax * (q2 + vmax -
     .    dsqrt((q2 + vmax)**2 + 4d0 * amp2 * q2))) /
     .    2d0 / (amp2 + vmax)
      t2 = (2d0 * amp2 * q2 + vmax * (q2 + vmax +
     .    dsqrt((q2 + vmax)**2 + 4d0 * amp2 * q2))) /
     .    2d0 / (amp2 + vmax)

      ttn = t1
      sin = 0d0
      distsit(0) = 0d0
      distart(0) = ttn
      do it1 = 1, nt1
        tto = ttn
        sio = sin
        ttn = t1 + (q2 - epsq2 - t1) * grt1(it1)
        sin = fsigt(s, q2, 0d0, vmin, vmax, ttn, 0)
        if (it1 .eq. 1) then
          distsit(it1) = sin * (ttn - tto)
        else
          distsit(it1) = distsit(it1 - 1) +
     .        (sin + sio) * (ttn - tto) / 2d0
        endif
        distart(it1) = ttn
      enddo

      ttn = q2 + epsq2
c      distart(nt1) = q2
      do it2 = 1, nt2
        tto = ttn
        sio = sin
        ttn = q2 + epsq2 + (t2 - q2 - epsq2) * grt2(it2)
        if (it2 .ne. nt2) then
          sin = fsigt(s, q2, 0d0, vmin, vmax, ttn, 0)
          distsit(nt1 + it2) = distsit(nt1 + it2 - 1) +
     .        (sin + sio) * (ttn - tto) / 2d0
        else
          distsit(nt1 + it2) = distsit(nt1 + it2 - 1) +
     .        sio * (ttn - tto)
        endif 
        distart(nt1 + it2) = ttn
      enddo

      sig = distsit(nt1 + nt2)
      sigrand = urand() * sig
      do itt = 1, nt1 + nt2
        if (distsit(itt) .gt. sigrand) then
          tgen = distart(itt - 1) + (distart(itt) -
     .        distart(itt - 1)) * (sigrand - distsit(itt - 1)) /
     .        (distsit(itt) - distsit(itt - 1))
          exit
        endif
      enddo

      v1 = max(vmin,
     .    (tgen - q2) * (dsqrt(tgen) + dsqrt(4d0 * amp2 + tgen)) /
     .    2d0 / dsqrt(tgen),
     .    (tgen - q2) * (dsqrt(tgen) - dsqrt(4d0 * amp2 + tgen)) /
     .    2d0 / dsqrt(tgen))
      vp = max(s - q2 * s / tgen, s + tgen - q2 - s * tgen / q2)

      vvn = v1
      sin = 0d0
      distsiv(0) = 0d0
      distarv(0) = vvn
      if (vp .lt. v1 .or. vp .gt. vmax) then
        do ivv = 1, nvv
          vvo = vvn
          sio = sin
          vvn = v1 + (vmax - v1) * grv(ivv, 0)
          sin = fsigtv(s, q2, tgen, vvn, 0d0, 0d0, 0)
          distsiv(ivv) = distsiv(ivv - 1) + (sin + sio) *
     .        (vvn - vvo) / 2d0
          distarv(ivv) = vvn
        enddo
      else
        do ivv = 1, nvv / 2
          vvo = vvn
          sio = sin
          vvn = v1 + (vp - v1) * grv0(ivv)
          sin = fsigtv(s, q2, tgen, vvn, 0d0, 0d0, 0)
          distsiv(ivv) = distsiv(ivv - 1) + (sin + sio) *
     .        (vvn - vvo) / 2d0
          distarv(ivv) = vvn
        enddo
        do ivv = nvv / 2 + 1, nvv
          vvo = vvn
          sio = sin
          vvn = vp + (vmax - vp) * grv0(ivv)
          sin = fsigtv(s, q2, tgen, vvn, 0d0, 0d0, 0)
          distsiv(ivv) = distsiv(ivv - 1) + (sin + sio) *
     .        (vvn - vvo) / 2d0
          distarv(ivv) = vvn
        enddo
      endif
     
      sig = distsiv(nvv)
      sigrand = urand() * sig
      do ivv = 1, nvv
        if (distsiv(ivv) .gt. sigrand) then
          vgen = distarv(ivv - 1) +
     .        (distarv(ivv) - distarv(ivv - 1)) *
     .        (sigrand - distsiv(ivv - 1)) /
     .        (distsiv(ivv) - distsiv(ivv - 1))
          exit
        endif
      enddo

      phn = 0d0
      sin = fsigtv(s, q2, tgen, vgen, 0d0, phn, 1)
      distsip(0) = 0d0
      distarp(0) = phn
      do iph = 1, nph
        pho = phn
        sio = sin
        phn = pi * grphi(iph)
        sin = fsigtv(s, q2, tgen, vgen, 0d0, phn, 1)
        distsip(iph) = distsip(iph - 1) + (sin + sio) *
     .      (phn - pho) / 2d0
        distarp(iph) = phn
      enddo

      sig = distsip(nph)
      sigrand = urand() * sig
      do iph = 1, nph
        if (distsip(iph) .gt. sigrand) then
          pgen = distarp(iph - 1) + (distarp(iph) -
     .        distarp(iph - 1)) * (sigrand - distsip(iph - 1)) /
     .        (distsip(iph) - distsip(iph - 1))
          return
        endif
      enddo

      end

!===============================================================================
      real*8 function fsigt(s,q2,phi,vmin,vmax,t,nad)
c the cross section of hard real photon emission ds/dQ2/dphi/dt
c numerical integration over v and analytical integration over phik
      implicit real*8(a-h,o-z)
      include 'elrad_grid.inc'
      include 'elrad_const.inc'

      v1=max(vmin,
     .    (t-q2)*(dsqrt(t)+dsqrt(4d0*amp2+t))/2d0/dsqrt(t),
     .    (t-q2)*(dsqrt(t)-dsqrt(4d0*amp2+t))/2d0/dsqrt(t))
      vp=max(s-q2*s/t,s+t-q2-s*t/q2)
      vvn=v1
      sin=0d0
      sum=0d0
      if(vp.lt.v1.or.vp.gt.vmax)then
        do ivv=1,nvv
          vvo=vvn
          sio=sin
          vvn=v1+(vmax-v1)*grv(ivv,0)
          sin=fsigtv(s,q2,t,vvn,phi,0d0,nad)
          sum=sum+(sin+sio)*(vvn-vvo)/2d0
        enddo
      else
        do ivv=1,nvv/2
          vvo=vvn
          sio=sin
          vvn=v1+(vp-v1)*grv0(ivv)
          sin=fsigtv(s,q2,t,vvn,phi,0d0,nad)
          sum=sum+(sin+sio)*(vvn-vvo)/2d0
        enddo
        do ivv=nvv/2+1,nvv
          vvo=vvn
          sio=sin
          vvn=vp+(vmax-vp)*grv0(ivv)
          sin=fsigtv(s,q2,t,vvn,phi,0d0,nad)
          sum=sum+(sin+sio)*(vvn-vvo)/2d0
        enddo
      endif
      fsigt=sum
      return

      end

!===============================================================================
      real*8 function fsigtv(s,q2,t,v,phi,phik,iphik)
c The cross section of hard real photon emission
c iphik=0 radiative cross section integrated over phik ds/dQ2/dphi/dt/dv
c iphik=1 radiative cross section ds/dQ2/dphi/dt/dv/dphik
c iphik=2 finite part of radiative cross section integrated over
c phik ds/dQ2/dphi/dt/dv
      implicit real*8(a-h,o-z)
      include 'elrad_const.inc'

      r=q2+v-t
      x=s-r-t
      sx=s-x
      sxp=s+x
      als=s**2-4d0*aml2*amp2
      ta=(t-q2)/r
      alq=sx**2+4d0*amp2*q2
      sqlq=dsqrt(alq)

      if(iphik.eq.0.or.iphik.eq.2)then
        b2=(-alq*ta+sxp*sx*ta+2.*sxp*q2)/2.
        b1=(-alq*ta-sxp*sx*ta-2.*sxp*q2)/2.
        c1=-(4.*(amp2*ta**2-sx*ta-q2)*aml2-(s*ta+q2)**2)
        c2=-(4.*(amp2*ta**2-sx*ta-q2)*aml2-(ta*x-q2)**2)
        bb=1./sqlq
        sc1=dsqrt(c1)
        sc2=dsqrt(c2)
        bis=(-b1/sc1/c1+b2/sc2/c2)
        bir=b2/sc2/c2+b1/sc1/c1
        bip=1d0/sc1+1d0/sc2
        bi12=(sxp*(sx*ta+2.*q2))/(sc1*sc2*(sc1+sc2))
      else
        tamin=(sx-sqlq)/2d0/amp2
        tamax=(sx+sqlq)/2d0/amp2
        sqrtmb=sqrt((ta-tamin)*(tamax-ta)*(s*x*q2-q2**2*amp2-aml2*alq))
        z1=(q2*sxp+ta*(s*sx+2d0*amp2*q2)-2d0*amp*cos(phik)*sqrtmb)/alq
        z2=(q2*sxp+ta*(s*sx-2d0*amp2*q2)-2d0*amp*cos(phik)*sqrtmb)/alq
        bb=1./sqlq/2.0/pi
        bi12=bb/(z1*z2)
        bip=bb*(1d0/z1+1d0/z2)
        bis=bb/z2**2+bb/z1**2
        bir=bb/z2**2-bb/z1**2
      endif

      bfir=aml2*bis-(q2+2.0d0*aml2)*bi12

c***** Beyond URA **************************************************************
      t11=4d0*(q2-2*aml2)*bfir
      t12=4d0*(ta*bfir)
      t13=-4d0*bb-2d0*ta**2*bi12
      th1n=t11/r**2+t12/r+t13
      t21=2d0*(s*x-amp2*q2)*bfir/amp2
      t22=(2d0*aml2*sxp*bir+sxp*sx*bip+2d0*(sx-2d0*amp2*ta)*bfir
     .    -ta*sxp**2*bi12)/2d0/amp2
      t23=(4d0*amp2*bb+(2d0*amp2*ta**2-sx*ta+4d0*aml2)*bi12-sxp*bip)
     .    /2d0/amp2
      th2n=t21/r**2+t22/r+t23
c*******************************************************************************
      call sffun(t,f1,f2)
      fsigtv=-alpha**3*barn/4d0/pi/als*(th1n*f1+th2n*f2)/t**2

      if(iphik.eq.2)then
        call sffun(q2,f10,f20)
c***** Beyond URA **************************************************************
        tm1=q2-2*aml2
c*******************************************************************************
        tm2=(s*(s-q2)-amp2*q2)/2d0/amp2
        fsigtv=fsigtv+alpha**3*barn/pi/als*bfir*(tm1*f10+tm2*f20)
     .      /(r*q2)**2
      endif

      end

!===============================================================================
      subroutine sffun(t,f1,f2)
c unpolarized (f1, f2) hadronic structure functions
      implicit real*8(a-h,l,m,o-z)
      include 'elrad_const.inc'

      call ffpro(t,gep,gmp)
      tau=t/4d0/amp2
      f1=4d0*tau*amp2*gmp**2
      f2=4d0*amp2*(gep**2+tau*gmp**2)/(1d0+tau)

      end

!===============================================================================
      subroutine ffpro(t, gep, gmp)
c electric and magnetic form factors of the proton
      implicit real*8(a-h, l, m, o-z)
      include 'elrad_const.inc'

      tau = t / 4d0 / amp2
      x1 = tau
      x2 = x1 * tau
      x3 = x2 * tau
      x4 = x3 * tau
      x5 = x4 * tau
      gep = (1d0 + 2.90966d0 * x1 - 1.11542229d0 * x2
     .    + 3.866171d-2 * x3) / (1.0 + 14.5187212d0 * x1
     .    + 40.88333d0 * x2 + 99.999998d0 * x3 + 4.579d-5 * x4
     .    + 10.3580447d0 * x5)
      gmp = 2.792782d0 * (1d0 - 1.43573d0 * x1 + 1.19052066d0 * x2
     .    + 2.5455841d-1 * x3) / (1d0 + 9.70703681d0 * x1
     .    + 3.7357d-4 * x2 + 6.0d-8 * x3 + 9.9527277d0 * x4
     .    + 12.7977739d0 * x5)

      end
