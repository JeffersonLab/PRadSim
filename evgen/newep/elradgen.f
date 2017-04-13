c      program main
      subroutine mainf(ranin,ranout,pebeam,ppq2,pvmin,pvcut,pphi,pppx,pppy
     .,pppz,ppp0,eppx,eppy,eppz,epp0,pphx,pphy,pphz,pph0,ppp,epp,ppph
     .,pq2,phige,sigmaborn,sigmatot,sigmarad,sigmabsv,pweight,ppt,ppv)
      implicit real*8(a-h,o-z)
      include 'const.inc'
      include 'output.inc'
      include 'pol.inc'
      include 'test.inc'
      real*4 ebeam,vmin,k1,k2
      dimension vprad(4),phrad(4),eprad(4)


cccccccccc random seed in input and outpout of the program to feed in the C++ code cccccccccccccc
      integer ranout,ranin
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
       pppx=0.0
       pppy=0.0
       pppz=0.0
       ppp0=0.0
       eppx=0.0
       eppy=0.0
       eppz=0.0
       epp0=0.0
       pphx=0.0
       pphy=0.0
       pphz=0.0
       pph0=0.0
       ppp=0.0
       epp=0.0
       ppph=0.0
       pq2=0.0
       phige=0.0
       sigmaborn=0.0
       ppt=0.0
       ppv=0.0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call elrad_init
c      print *,amp2
        
        ebeam=pebeam
	q2=ppq2
        vmin=pvmin
        phi=pphi
c        print *,ebeam,q2,vmin,phi
        s=2d0*ebeam*amp
c        vmaxkin=s-q2-(0.938272d0)**2*q2/s

c*************Beyond URA*************************************************
        vmaxkin=(2d0*q2*(s**2-4*aml2*amp2-q2*(s+aml2+amp2)))/
     .  (q2*(s+2*aml2)+sqrt(q2*(s**2-4d0*amp2*aml2)*(q2+3*aml2)))
c************************************************************************

c vcut -- cut on inelasticity. If vcut=0d0 no cut!!!
	vcut=pvcut
 	plrun=0d0
 	pnrun=0d0
	thetapn=pi/2d0!0d0!48d0*pi/180d0
        phipn=0d0!pi/4d0
c        phipn=pi/2d0

cccccc iy is initialize in gene.cppcccccccccc
c	 open(12,file='rnd.dat')
c	 read(12,*)iy 	
c	 close(12)
c         iy=urand(iy)	

         iy=ranin
c	write(*,'(4g12.4)')vpgen

c      call elradgen(ebeam,q2,phi,vmin,vcut)

c      do i=1,10

cccccccccccccccccccc        CALL of ELRADGEN         cccccccccccccccccccc
c      print *,2,vmin,vcut
      call elradgen(ebeam,q2,phi,vmin,vcut,psigmaborn,psigmatot,
     .psigmarad,psigmabsv,pt,pv)
c     scattered proton 4-momentum
       pppx=vprad(1)
       pppy=vprad(2)
       pppz=vprad(3)
       ppp0=vprad(4)
c     scattered electron 4-momentum
       eppx=eprad(1)
       eppy=eprad(2)
       eppz=eprad(3)
       epp0=eprad(4)
c     radiated photon 4-momentum
       pphx=phrad(1)
       pphy=phrad(2)
       pphz=phrad(3)
       pph0=phrad(4)

c      magnitude of the 4-momenta
       ppp=dsqrt(vprad(1)**2+vprad(2)**2+vprad(3)**2)
       epp=dsqrt(eprad(1)**2+eprad(2)**2+eprad(3)**2)
       ppph=dsqrt(phrad(1)**2+phrad(2)**2+phrad(3)**2)

c      4-momentum transfer
       pq2=-((phrad(4)+vprad(4))**2-(phrad(1)+vprad(1))**2
     .-(phrad(2)+vprad(2))**2-(phrad(3)+vprad(3))**2)
       phige=phigen*57.296

c      cross sections
       sigmaborn=psigmaborn
       sigmarad=psigmarad
       sigmabsv=psigmabsv
c      total cross section = rad + bsv
       sigmatot=psigmatot
       ranout=iy
c      weight = sigmatot/sigma/born
       pweight=weight
c      t and v generated variables
       ppt=pt
       ppv=pv
c       print *,ranout
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
c      print *,tgen,vgen,phigen
c      print *,vprad
c      print *,phrad,weight
c      print *,
c     .phrad(4)**2-phrad(1)**2-phrad(2)**2-phrad(3)**2
c      print *,
c     .(vprad(4)+amp)**2
c     .-vprad(1)**2-vprad(2)**2-vprad(3)**2,amp**2
c      print *,
c     .-((phrad(4)+vprad(4))**2
c     .-(phrad(1)+vprad(1))**2-(phrad(2)+vprad(2))**2
c     .-(phrad(3)+vprad(3))**2),q2


c      enddo
c      pause
      return
      end

      subroutine elrad_init
c subroutine with definition all necessary constant for generation
      implicit real*8(a-h,o-z)
      include 'const.inc'
      include 'pol.inc'
	data itest/0/
	aml2=.261112d-6
      aml=dsqrt(aml2)
	amp=.938272d0
	amp2=amp*amp
	alpha=.729735d-2
	pi=atan(1d0)*4d0
	barn=.389379d6
	call grid_init
      end



      subroutine elradgen(ebeam,q2,phi,vvmin,vcut,psigmaborn,psigmatot,
     .psigmarad,psigmabsv,pt,pv)
c      subroutine elradgen(ebeam,q2,phi,vvmin,vcut)
c basic subroutine for generation of one radiative event
c and reconstruction of four momenta of final particles
      implicit real*8(a-h,o-z)
      real*4 ebeam,vvmin,urand

      real*8 lambs,lambm,lamb0x,l0x,lm,ls,qm,aa,bb,qmbis,zxc
      real*8 sqlm,sqls,sql0x,logm,logsq2,jy0,s,q2
      
      include 'const.inc'
      include 'output.inc'
      include 'test.inc'
      include 'par.inc'
      external fsadd,fs
      include 'grid.inc'
      dimension vprad(4),phrad(4),eprad(4)
	dimension distsit(0:nt1+nt2),distarg(0:nt1+nt2)
      dimension distsiv(0:nvv),distarv(0:nvv)
      dimension distsip(0:nph),distarp(0:nph)
	data ikey/0/
	phib = phi

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ppx=0.0
       ppy=0.0
       ppz=0.0
       pp0=0.0
       epx=0.0
       epy=0.0
       epz=0.0
       ep0=0.0
       pphigen=0.0
       ranout=0
       pt=0.0
       pv=0.0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      s=2d0*ebeam*amp
      vmin=vvmin
      ss=s
      q22=q2
c      vmaxkin=s-q2-(0.938272d0)**2*q2/s
      vmaxkin=2*q2*(s*s-4*amp2*aml2-q2*(s+aml2+amp2))/
     .(q2*(s+2*aml2)+dsqrt(q2*((s*s-4*amp2*aml2)*(q2+4*aml2))))
c      print *,q2,s,vmaxkin,vmaxkin2
	if(vcut.lt.1d-9)then
        vmax=vmaxkin-1d-4
	else
	vmax=vcut
	endif
c       print *,1,vmin,vmax
      vminn=vmin
      vmaxx=vmax
      t2=(2d0*amp2*q2+vmax*(q2+vmax+
     .dsqrt((q2+vmax)**2+4d0*amp2*q2)))/2d0/(amp2+vmax)
c      t1=amp2*q2**2/(amp2+vmax)/t2
      t1=(2d0*amp2*q2+vmax*(q2+vmax-
     .dsqrt((q2+vmax)**2+4d0*amp2*q2)))/2d0/(amp2+vmax)
        if(itest.eq.2)goto 881
        if(itest.eq.3)goto 882

      t22=(2d0*amp2*q2+vmin*(q2+vmin+
     .dsqrt((q2+vmin)**2+4d0*amp2*q2)))/2d0/(amp2+vmin)
c      t11=amp2*q2**2/(amp2+vmin)/t22
      t11=(2d0*amp2*q2+vmin*(q2+vmin-
     .dsqrt((q2+vmin)**2+4d0*amp2*q2)))/2d0/(amp2+vmin)
        if(ikey.eq.0)then
        if(itest.eq.1)then
	  ikey=1

          step=(t2-t1)/dble(nbin)

          do i=0,nbin
          bin(i)=t1+dble(i)*step
          enddo
          endif

	 epsq2=1d-13
	 ttn=t1
	 sin=0d0
	 distsit(0)=0d0
	 distarg(0)=ttn

	 do it1=1,nt1
	   tto=ttn
	   sio=sin
	   ttn=t1+(q2-t1-epsq2)*grt1(it1)
	 sin=fsigt(s,q2,phi,vmin,vmax,ttn,0)
	   distsit(it1)=distsit(it1-1)+(sin+sio)*(ttn-tto)/2d0
	   distarg(it1)=ttn
	 enddo

	 ttn=q2+epsq2
	 	
	 distarg(nt1)=q2
	 do it2=1,nt2
	   tto=ttn
	   sio=sin
	   ttn=q2+epsq2+(t2-q2-epsq2)*grt2(it2)
	   sin=0d0
	  if(it2.ne.nt2)sin=fsigt(s,q2,phi,vmin,vmax,ttn,0)
	   distsit(nt1+it2)=distsit(nt1+it2-1)+(sin+sio)*(ttn-tto)/2d0
	   distarg(nt1+it2)=ttn
	 enddo
	 siborn=fsib(s,q2,phi)
         siborn2=fsib2(s,q2,phi)
c         print *,siborn,siborn2
c****************Beyond URA********************************************
c       q2=0.0002
c       s=2.064198838
       qm=q2+2d0*aml2
       lambm=q2**2+4d0*aml2*q2
       vmaxkin=(2d0*q2*(s**2-4*aml2*amp2-q2*(s+aml2+amp2)))/
     .  (q2*(s+2*aml2)+sqrt(q2*(s**2-4d0*amp2*aml2)*(q2+3*aml2)))
       

       lambs=s**2-4d0*aml2*amp2
       lamb0x=(s-q2)**2-4d0*aml2*amp2
       sqlm=dsqrt(lambm)
       sqls=dsqrt(lambs)
       sql0x=dsqrt(lamb0x)

       lm=(dlog((sqlm+q2)/(sqlm-q2)))/sqlm
       ls=(dlog((s+sqls)/(s-sqls)))/sqls
       l0x=(dlog((s-q2+sql0x)/(s-q2-sql0x)))/sql0x
       jy0=2.0*(qm*lm-1)

       aa=(s*(s-q2)-2d0*amp2*q2-4d0*amp2*aml2)/(2d0*amp2)
       bb=(q2*(s*(s-q2)-amp2*q2)-aml2*(q2*q2+4d0*amp2*q2))/amp2
       xpp=s-q2-2d0*amp2

       call sphi(qm,lambm,aa,bb,sphi1)
c       call sphibis(qm,lambm,s,s-q2-2*amp2,lamb0x,sphi3)
c       print *,sphi1

c********* S_phi in URA**************************
       logm=dlog(q2/aml2)
       logsq2=dlog(s/(s-q2))

       sphi2=0.5d0*logm*logm-logm*dlog(s*(s-q2)/(aml2*amp2))
     . - 0.5d0*logsq2*logsq2
     . +fspen((s*(s-q2)-q2*amp2)/(s*(s-q2)))-pi*pi/3d0
c*************************************************

       delinfx=(dlog(q2/aml2)-1d0)*dlog(vmax**2/s/(s-q2))
       deltaS=jy0*log(2d0)+.5d0*(s*ls+(s-q2)*l0x)
c       deltaS=2d0*jy0*dlog(vmax/aml/amp)+.5d0*(s*ls+(s-q2)*l0x)
       delinf=jy0*log(vmax/2/aml/amp)
c       delinf=jy0*dlog(vmax**2/s/(s-q2))
       deltavert=-2d0+(1.5d0*q2+4*aml2)*lm-(qm/sqlm)*(.5d0*lambm*lm**2+
     . 2d0*fspen((2d0*sqlm)/(q2+sqlm))-pi**2/4d0)
c       deltavert=-2d0+(1.5d0*q2+4*aml2)*lm-(qm/sqlm)*(.5d0*lambm*lm**2+
c     . 2d0*fspen((2d0*sqlm)/(q2+sqlm))-pi**2/2d0)
       delta=alpha/pi*(deltaS+sphi1+deltavert+delinf-delinfx+vacpol(q2))
c       delta=alpha/pi*(deltaS+sphi1+deltavert-delinf+vacpol(q2))
c       print *,delinf,delinfx
c       print *,delta,deltaS+sphi1+deltavert,vacpol(q2)
       
c*********MASCARAD delta**************************************************
c       deltax=alpha/pi*(1.5d0*dlog(q2/aml2)-2d0-0.5d0*dlog(s/(s-q2))**2
c     . +fspen(1d0-amp2*q2/s/(s-q2))-pi**2/6d0+vacpol(q2))

      deladd=2d0*alpha/pi*dlog(vmax/vmin)*(1d0-dlog(q2/aml2))
c       print *,deladd
       call simpxx(t11+1d-8,t22-1d-8,5000,1d-2,fsadd,sigadd)
c       print *,sigadd
       extai1=exp(alpha/pi*delinfx)
c       extai1=exp(alpha/pi*delinf)
c       print *,extai1
c BSV-part of the total cross section	
	 sibsv=siborn*extai1*(1d0+delta)+siborn*deladd+sigadd
c part of the total cross section with
c hard photon emission	
	 sirad=distsit(nt1+nt2)
c         print *,sirad
	 sitot=sibsv+sirad
c       print *,siborn,sitot,sirad,sibsv


        endif
       weight=sngl(sitot/siborn)

	if(itest.eq.1)weight=sngl(sirad)
      sirand=urand(iy)*sitot

cccccccccccccccccccccccccccccccccccccccccccccccccccc
      ranout=iy
c      print *,ranout
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c choosing "success" events	
      if(sirand.le.sibsv)then
	ich=0
	tgen=0d0
	vgen=0d0
	phigen=0d0	
      vprad(1)=sngl(-dsqrt(q2*vmaxkin/s)*dcos(phib))
      vprad(2)=sngl(-dsqrt(q2*vmaxkin/s)*dsin(phib))
      vprad(3)=sngl(q2*(2.*amp2+s)/2./amp/s)
      vprad(4)=sngl(q2/2./amp)
	do i=1,4
	phrad(i)=0.
	enddo
      psigmaborn=siborn
      psigmatot=sitot
      psigmarad=sirad
      psigmabsv=sibsv
	return
      endif
	ich=1
c generation of t-variable
	do itt=1,nt1+nt2
	  if(distsit(itt).gt.sirand-sibsv)then
	    tgen=distarg(itt-1)+(distarg(itt)-distarg(itt-1))*
     .	   (sirand-sibsv
     .	     -distsit(itt-1))/(distsit(itt)-distsit(itt-1))	
      if(itest.eq.1)return
	  goto 881
	  endif
      enddo
881   continue
        if(itest.eq.2)ich=1

	v1=max(vmin,
     .	  (tgen-q2)*(sqrt(tgen)+sqrt(4d0*amp2+tgen))/2d0/sqrt(tgen),
     .	  (tgen-q2)*(sqrt(tgen)-sqrt(4d0*amp2+tgen))/2d0/sqrt(tgen))
	vp=max(s-q2*s/tgen,s+tgen-q2-s*tgen/q2)
        if(ikey.eq.0)then
        if(itest.eq.2)then
	  ikey=1
        step=(vmax-v1)/dble(nbin)
        do i=0,nbin
        bin(i)=v1+dble(i)*step
        enddo
        endif

	 vvn=v1
	 sin=0d0
	 distsiv(0)=0d0
	 distarv(0)=vvn
	if(vp.lt.v1.or.vp.gt.vmax)then
	 do ivv=1,nvv
	   vvo=vvn
	   sio=sin
	   vvn=v1+(vmax-v1)*grv(ivv,0)
	   sin=fsigtv(s,q2,tgen,vvn,phi,0d0,0)
	   distsiv(ivv)=distsiv(ivv-1)+(sin+sio)*(vvn-vvo)/2d0
	   distarv(ivv)=vvn
	 enddo
	else
	 do ivv=1,nvv/2
	   vvo=vvn
	   sio=sin
	   vvn=v1+(vp-v1)*grv0(ivv)
	   sin=fsigtv(s,q2,tgen,vvn,phi,0d0,0)
	   distsiv(ivv)=distsiv(ivv-1)+(sin+sio)*(vvn-vvo)/2d0
	   distarv(ivv)=vvn
	 enddo
	 do ivv=nvv/2+1,nvv
	   vvo=vvn
	   sio=sin
	   vvn=vp+(vmax-vp)*grv0(ivv)
	   sin=fsigtv(s,q2,tgen,vvn,phi,0d0,0)
	   distsiv(ivv)=distsiv(ivv-1)+(sin+sio)*(vvn-vvo)/2d0
	   distarv(ivv)=vvn
	 enddo
	endif
	goto 44
	if(vp.lt.v1.or.vp.gt.vmax)iivn=0
	xivp=(vp-v1)/(vmax-v1)*40d0
	if(xivp.lt.7.5d0)then
	  iivn=1
	elseif(xivp.lt.12.5d0)then
	  iivn=2
	elseif(xivp.lt.17.5d0)then
	  iivn=3
	elseif(xivp.lt.22.5d0)then
	  iivn=4
	elseif(xivp.lt.27.5d0)then
	  iivn=5
	elseif(xivp.lt.32.5d0)then
	  iivn=6
	else
	  iivn=7
	endif
	 vvn=v1
	 sin=0d0
	 distsiv(0)=0d0
	 distarv(0)=vvn
	 do ivv=1,nvv
	   vvo=vvn
	   sio=sin
	   vvn=v1+(vmax-v1)*grv(ivv,iivn)
	   sin=fsigtv(s,q2,tgen,vvn,phi,0d0,0)
	   distsiv(ivv)=distsiv(ivv-1)+(sin+sio)*(vvn-vvo)/2d0
	   distarv(ivv)=vvn
	 enddo

44        endif
	if(itest.eq.2)weight=sngl(distsiv(nvv))
	 sirrvv=urand(iy)*distsiv(nvv)
c generation of v-variable
	 do ivv=1,nvv
	   if(distsiv(ivv).gt.sirrvv)then
	     vgen=distarv(ivv-1)+(distarv(ivv)-distarv(ivv-1))*
     .	   (sirrvv-distsiv(ivv-1))/(distsiv(ivv)-distsiv(ivv-1))
	if(itest.eq.2)return

	   goto 882
	   endif
	 enddo
882	 continue
        if(itest.eq.3)ich=1
        if(ikey.eq.0)then
        if(itest.eq.3)then
	  ikey=1
        step=2d0*pi/dble(nbin)
        do i=0,nbin
        bin(i)=-pi+dble(i)*step
c       print *,i,argbin(i),vmin,vmax
        enddo
c       stop
        endif

	 phn=0d0
	 sin=0d0
	 distsip(0)=0d0
	 distarp(0)=phn
	 do iph=1,nph
	   pho=phn
	   sio=sin
	   phn=pi*grphi(iph)
	   sin=fsigtv(s,q2,tgen,vgen,phi,phn,1)
	   distsip(iph)=distsip(iph-1)+(sin+sio)*(phn-pho)/2d0
	   distarp(iph)=phn
	 enddo
        endif
        if(itest.eq.3)weight=sngl(distsip(nph))
	 sirrph=urand(iy)*distsip(nph)
c generation of phigen-variable
	 do iph=1,nph
	   if(distsip(iph).gt.sirrph)then
	     phigen=distarp(iph-1)+(distarp(iph)-distarp(iph-1))*
     .	   (sirrph-distsip(iph-1))/(distsip(iph)-distsip(iph-1))
	   if(urand(iy).gt.0.5)phigen=-phigen
	if(itest.eq.3)return
	   goto 883
	   endif
	 enddo
883	 continue
	alq=(q2+vgen)**2+4.*amp2*q2
	slq=sqrt(alq)
	al1=s*(q2+vgen)+2.*amp2*q2
	al2=tgen*(q2+vgen)+2.*amp2*(q2+tgen)
	sl3=sqrt(tgen*vgen*(q2-tgen+vgen)-amp2*(q2-tgen)**2)
	sl4=sqrt(q2*s*(vmaxkin-vgen))
	cph=dcos(phi)
	sph=dsin(phi)
	cpg=dcos(phigen)
	spg=dsin(phigen)
      if(urand(iy).gt.0.5)sl3=-sl3
c four momenta reconstruction	
      vprad(1)=
     .-sngl((sl3*(al1*cph*cpg-slq*s*sph*spg)+al2*sl4*cph)/alq/s)
      vprad(2)=
     .-sngl((sl3*(al1*sph*cpg+slq*s*cph*spg)+al2*sl4*sph)/alq/s)
      vprad(3)=sngl((al1*al2-4*amp2*sl3*sl4*cpg)/2/alq/s/amp)
      vprad(4)=sngl((tgen+2*amp*amp)/2./amp)
      phrad(1)=-sngl((-sl3*(al1*cph*cpg-slq*s*sph*spg)
     . +(alq-al2)*sl4*cph)/alq/s)
      phrad(2)=-sngl((-sl3*(al1*sph*cpg+slq*s*cph*spg)
     . +(alq-al2)*sl4*sph)/alq/s)
      phrad(3)=sngl((al1*(alq-al2)+4*amp2*sl3*sl4*cpg)/2/alq/s/amp)
      phrad(4)=sngl((q2+vgen-tgen)/2./amp)
      
      eprad(1)=sngl(sl4*cpg/s)
      eprad(2)=sngl(sl4*spg/s)
      eprad(3)=sngl((s*s-al1)/(2*amp*s))
      eprad(4)=sngl((s-q2-vgen)/(2*amp))
c      print *,eprad(1),eprad(2),eprad(3),eprad(4),vprad(1),vprad(2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      psigmaborn=siborn
      psigmatot=sitot
      psigmarad=sirad
      psigmabsv=sibsv
      pt=tgen
      pv=vgen
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      end



      double precision function fsadd(t)
c finite part of radiative cross section integrated over 
c v and phik ds/dQ2/dphi/dt
      implicit real*8(a-h,o-z)
      include 'par.inc'
      fsadd=fsigt(ss,q22,phib,0d0,vminn,t,2)
      end

      double precision function fsigt(s,q2,phi,vmin,vmax,t,nad)
c the cross section of hard real photon emission ds/dQ2/dphi/dt
      implicit real*8(a-h,o-z)
      include 'pol.inc'
c	if(abs(plrun*pnrun).lt.1d-12)then
c	fsigt=fsigtan(s,q2,vmin,vmax,t,nad)
c	else	
	fsigt=fsigtnu(s,q2,phi,vmin,vmax,t,nad)	
c	endif
      end

      double precision function fsigtnu(s,q2,phi,vmin,vmax,t,nad)
c the cross section of hard real photon emission ds/dQ2/dphi/dt
c numerical integration over v and analytical integration over phik
      implicit real*8(a-h,o-z)
        include 'grid.inc'
      include 'const.inc'
	v1=max(vmin,
     .	  (t-q2)*(sqrt(t)+sqrt(4d0*amp2+t))/2d0/sqrt(t),
     .	  (t-q2)*(sqrt(t)-sqrt(4d0*amp2+t))/2d0/sqrt(t))
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
       fsigtnu=sum
	 return
      end


      double precision function fsigtan(s,q2,vmin,vmax,t,nad)
c the cross section of hard real photon emission ds/dQ2/dphi/dt
c analytical integration over v and phik
      implicit real*8(a-h,o-z)
      include 'const.inc'
      v1=max(vmin,
     .(t-q2)*(dsqrt(t)+dsqrt(4d0*amp2+t))/2d0/dsqrt(t),
c     .2d0*amp2*(q2-t)/(dsqrt(t)+dsqrt(4d0*amp2+t))/dsqrt(t))!,
     .(t-q2)*(dsqrt(t)-dsqrt(4d0*amp2+t))/2d0/dsqrt(t))
      v2=vmax
      c11=(t-q2)*(s-q2)+q2*v1
      c21=S*(q2-t)+t*v1
      c31=v1*t*(v1-t+q2)-amp2*(q2-t)**2
      c12=(t-q2)*(s-q2)+q2*v2
      c22=S*(q2-t)+t*v2
      c32=v2*t*(v2-t+q2)-amp2*(q2-t)**2
      c41=dsqrt(c11**2+4d0*aml2*c31)
      c51=dsqrt(c21**2+4d0*aml2*c31)
      c42=dsqrt(c12**2+4d0*aml2*c32)
      c52=dsqrt(c22**2+4d0*aml2*c32)
	al11=2d0*t*(q2-t+2d0*v1)*aml2+c11*q2+c41*sqrt(4d0*t*aml2+q2**2)
	al12=2d0*t*(q2-t+2d0*v2)*aml2+c12*q2+c42*sqrt(4d0*t*aml2+q2**2)
	al21=2d0*dsqrt(t)*(q2-t+2d0*v1)*aml2+c21*dsqrt(t)+c51
     .*sqrt(4d0*aml2+t)
	al22=2d0*dsqrt(t)*(q2-t+2d0*v2)*aml2+c22*dsqrt(t)+c52
     .*sqrt(4d0*aml2+t)
      dl1=log(al12)-log(al11)
      dl2=log(al22)-log(al21)
      d12=c42*c12-c41*c11
      d11=c42-c41
      d12=abs(c12)*c12-abs(c11)*c11
      d11=abs(c12)-abs(c11)
      d10=c12/c42-c11/c41
      d22=abs(c22)*c22-abs(c21)*c21
      d21=abs(c22)-abs(c21)
      d20=c22/c52-c21/c51
      bif=dlog(
     .(v2+q2+dsqrt((v2+q2)**2+4d0*amp2*q2))/
     .(v1+q2+dsqrt((v1+q2)**2+4d0*amp2*q2)))
      bip0=(q2-t)*(s/q2**2*dl1-(s-t)/t**2*dl2)+d11/q2**2+d21/t**2
      bip1=dl1/q2+dl2/t
      bi210=(q2-t)/2d0*s**2/q2**2*d10
      bi211=(s/q2*d10)/2d0
      bi2p2=(d10-q2/t*d20)/2d0/(q2-t)
      bid0=(q2-t)*(s**2/q2**3*dl1-(s-t)**2/t**3*dl2)
     .+2d0*(s*d11/q2**3+(s-t)*d21/t**3)
     .+(d12/q2**3-d22/t**3)/2d0/(q2-t)
      bid1=s/q2**2*dl1+(s-t)/t**2*dl2+(d11/q2**2-d21/t**2)/(q2-t)
      bid2=(dl1/q2-dl2/t)/(q2-t)

      th1a=-4d0*bif+4d0*t*bi2p2-2d0*(Q2**2+t**2)*bid2
      th2a=(4d0*amp2*bif-t*(bip0+2d0*bid0)
     .+t*(2d0*s-t)*bip1+4d0*(bi210)
     .-4d0*(2d0*s-t)*bi211+4d0*(s**2-(amp2+s)*t)
     .*bi2p2+t*(q2+4d0*s-3d0*t)*bid1+(t*(q2*t-(2d0*s-t)**2)
     .+2d0*amp2*(t**2+q2**2))*bid2)/2d0/amp2
      call sffun(t,f1,f2,f3,f4)
      fsigtan=-alpha**3*barn/4d0/pi/s**2*(th1a*f1+th2a*f2)/t**2
      if(nad.eq.2)then
      tm1=q2
      tm2=(s*(s-q2)-amp2*q2)/(2.d0*amp2)
      call sffun(q2,f10,f20,f30,f40)
      fir=bi2p2-q2*bid2
      fsigtan=fsigtan+fir*alpha**3*barn/pi/s**2*(tm1*f10+tm2*f20)/q2**2
      endif

      end



      double precision function fsigtv(s,q2,t,v,phi,phik,iphik)
c The cross section of hard real photon emission 
c
c iphik=0 radiative cross section integrated over phik ds/dQ2/dphi/dt/dv
c
c iphik=1 radiative cross section ds/dQ2/dphi/dt/dv/dphik
c
c iphik=2 finite part of radiative cross section integrated over
c phik ds/dQ2/dphi/dt/dv

      implicit real*8(a-h,o-z)
      include 'const.inc'
      include 'pol.inc'

      r=q2+v-t

      x=s-r-t
      sx=s-x
      sxp=s+x
      ta=(t-q2)/r
      alq=sx**2+4d0*amp2*q2
      sqlq=dsqrt(alq)
      sqls=sqrt(s**2-4d0*amp2*aml2)

      ym=q2+2d0*aml**2


       call polpar(s,q2,x,phi,ta,sps,spe,ccps,ccpe,de,apq,apn,dk2ks,
     .dapks,dksp1)
      if(iphik.eq.0.or.iphik.eq.2)then

       b2=(-alq*ta+sxp*sx*ta+2.*sxp*q2)/2.
       b1=(-alq*ta-sxp*sx*ta-2.*sxp*q2)/2.
       c1=-(4.*(amp2*ta**2-sx*ta-q2)*aml2-(s*ta+q2)**2)
       c2=-(4.*(amp2*ta**2-sx*ta-q2)*aml2-(ta*x-q2)**2)
       bb=1/sqlq
       sc1=dsqrt(c1)
       sc2=dsqrt(c2)
       bis=(-b1/sc1/c1+b2/sc2/c2)
       bir=b2/sc2/c2+b1/sc1/c1
       bip=1d0/sc1+1d0/sc2
       bi12=(sxp*(sx*ta+2.*q2))/(sc1*sc2*(sc1+sc2))
	sl3=0d0
      else
	sl3=sqrt(t*v*(q2-t+v)-amp2*(q2-t)**2)
       tamin=(sx-sqlq)/2d0/amp2
       tamax=(sx+sqlq)/2d0/amp2
       sqrtmb=sqrt((ta-tamin)*(tamax-ta)*(s*x*q2-q2**2*amp2-aml2*alq))
       z1=(q2*sxp+ta*(s*sx+2d0*amp2*q2)-2d0*amp*cos(phik)*sqrtmb)/alq
       z2=(q2*sxp+ta*(x*sx-2d0*amp2*q2)-2d0*amp*cos(phik)*sqrtmb)/alq
       bb=1./sqlq/pi
       bi12=bb/(z1*z2)
       bip=bb*(1d0/z1+1d0/z2)
       bis=bb/z2**2+bb/z1**2
       bir=bb/z2**2-bb/z1**2

c************New bim (correspond to a new F_{1-}=F*(1/z1 -1/z2) in MASCARD paper)*********
       bim=bb*(1d0/z1-1d0/z2)
c*****************************************************************************************
      endif

      bfir=aml2*bis-q2*bi12
      sis=((2.*bip+bir*ta)*sps+bis*ccps)/2.
      eis=((2.*bip+bir*ta)*spe+bis*ccpe+4.*sin(phik)*de*sl3/r
     ./sqlq*bis)/2.
      sir=((2.*bi12+bis)*sps*ta+bir*ccps)/2.
      si12=(bi12*ccps+bip*sps)/2.
      ei12=(bi12*ccpe+bip*spe+4.*sin(phi_k)*de*sl3/r/sqlq*bi12)/2.
      ois=((2.*bip+bir*ta)*(ccpe*sps+ccps*spe)+(ccpe*ccps+
     . spe*sps*ta**2)*bis+4.*(2.*bb+bi12*ta**2)*spe*sps+8.*sin(phi_k)
     .*de*sl3/r/sqlq*sis)/4.
      oir=( ((2.*bi12+bis)*(ccpe*sps+ccps*spe)*ta+(ccpe*ccps+
     . spe*sps*ta**2)*bir+4.*bip*spe*sps*ta)+8.*sin(phi_k)*de
     .*sl3/r/sqlq*sir)/4.
      oi12=((ccpe*sps+ccps*spe)*bip+(ccpe*ccps+spe*sps*ta**2)
     .*bi12+4.*bb*spe*sps+8.*sin(phi_k)*de*sl3/r/sqlq*si12)/4.
      sfir=aml2*sis-q2*si12
      efir=aml2*eis-q2*ei12
      ofir=aml2*ois-q2*oi12

c***********Modification of the thetas beyon Ultra-relativistic limit*********
      t11x=4d0*q2*bfir
      t12x=4d0*ta*bfir
      t13x=-4d0*bb-2d0*ta**2*bi12
      th1nx=t11x/r**2+t12x/r+t13x
      t21x=2d0*(s*x-amp2*q2)*bfir/amp2
      t22x=(2d0*aml2*sxp*bir+sxp*sx*bip+2d0*(sx-2d0*amp2*ta)*bfir
     .-ta*sxp**2*bi12)/2d0/amp2
      t23x=(4d0*amp2*bb+(2d0*amp2*ta**2-sx*ta)*bi12-sxp*bip)/2d0/amp2
      th2nx=t21x/r**2+t22x/r+t23x

      t11=4d0*(q2-2*aml2)*bfir
c-2d0*bil2*aml2**2
c      print *,t11k,t11l,t11
      t12=4d0*(ta*bfir)
c+aml2*bim)
      t13=-4d0*bb-2d0*ta**2*bi12
      th1n=t11/r**2+t12/r+t13
      t21=2d0*(s*x-amp2*q2)*bfir/amp2
      t22=(2d0*aml2*sxp*bir+sxp*sx*bip+2d0*(sx-2d0*amp2*ta)*bfir
     .-ta*sxp**2*bi12)/2d0/amp2
      t23=(4d0*amp2*bb+(2d0*amp2*ta**2-sx*ta+4d0*aml2)*bi12-sxp*bip)
     ./2d0/amp2
      th2n=t21/r**2+t22/r+t23
c****************************************************************************

      t31=8.*(apq*dk2ks-dapks*q2)*aml*bfir/amp
      t32=-2.*(2.*(bi12*dk2ks*ta-2.*sfir)*apq+(sis-sir-2d0*si12)*q2*apn
     .+4.*dapks*bfir*ta-4.*((2.*ei12-eis)*dk2ks)*aml2)*aml/amp
      t33=2.*(((2.*si12+sir-sis)*apn*ta-2.*dk2ks*ei12*ta
     . -6.*ofir-oir*q2+ois*ym)-4.*aml2*oi12)*aml/amp
      t34=2.*(2.*oi12-oir+ois)*aml*ta/amp
      th3n=t31/r**2+t32/r+t33+t34*r
      t4t1=4.*(2.*dksp1*q2-dk2ks*sx)*aml*bfir/amp2
      t4e1=4.*(2.*dksp1*q2-dk2ks*sx)*aml*efir/amp2
      t4t2=(2.*(sxp-2.*sx)*sfir+2.*bi12*dk2ks*sx*ta+8.*
     . dksp1*bfir*ta+sxp*(sis*ym-sir*q2)-4.*((2.*bi12-bis)*dk2ks+
     .(sis-si12)*sxp)*aml2)*aml/amp2
      t4e2=(2.*(sxp-2.*sx)*ofir+2.*ei12*dk2ks*sx*ta+8.*
     . dksp1*efir*ta+sxp*(ois*ym-oir*q2)-4.*((2.*ei12-eis)*dk2ks+
     .(ois-oi12)*sxp)*aml2)*aml/amp2
      t4t3=(((sxp*ta-ym)*sis-(sxp*ta-q2)*sir+2.*bi12*dk2ks
     . *ta+6.*sfir-2.*si12*sxp*ta)+4.*aml2*si12)*aml/amp2	
      t4e3=(((sxp*ta-ym)*ois-(sxp*ta-q2)*oir+2.*ei12*dk2ks
     . *ta+6.*ofir-2.*oi12*sxp*ta)+4.*aml2*oi12)*aml/amp2	
      t4t4=-(2.*si12-sir+sis)*aml*ta/amp2
      t4e4=-(2.*oi12-oir+ois)*aml*ta/amp2
      t41=apq*t4t1/amp	 
      t42=(apq*t4t2-t4e1)/amp	 
      t43=(apq*t4t3-t4e2)/amp	 
      t44=(apq*t4t4-t4e3)/amp	 
      t45=-t4e4/amp	 
      th4n=t41/r**2+t42/r+t43+t44*r+t45*r**2
      call sffun(t,f1,f2,f3,f4)
c      plrun=0.0
      fsigtv=-alpha**3*barn/4d0/pi/s**2*(th1n*f1+th2n*f2)/t**2
c     .+plrun*pnrun*(th3n*f3+th4n*f4))/t**2

      fsigtvx=-alpha**3*barn/4d0/pi/s**2*(th1nx*f1+th2nx*f2)/t**2
c     .+plrun*pnrun*(th3n*f3+th4n*f4))/t**2
c      print *,fsigtv,fsigtvx
      if(iphik.eq.2)then
      call sffun(q2,f10,f20,f30,f40)
       call polpar(s,q2,s-q2,phi,0d0,sps,spe,ccps,ccpe,de,apq0,apn0,
     .dk2ks0,dapks0,dksp10) 
c************* Beyond URA *****************************
c        tm1=q2
        tm1=q2-2*aml2
c******************************************************
	tm2=(s*(s-q2)-amp2*q2)/2d0/amp2
        tm3=2.*(apq0*dk2ks0-dapks0*q2)*aml/amp
        tm4=-apq0*(dk2ks0-2.*dksp10)*q2*aml/amp**3
	fsigtv=fsigtv+alpha**3*barn*bfir*(tm1*f10+tm2*f20
     .+plrun*pnrun*(tm3*f30+tm4*f40))/pi/(s*q2*r)**2
        fsigtvx=fsigtvx+alpha**3*barn*bfir*(tm1*f10+tm2*f20
     .+plrun*pnrun*(tm3*f30+tm4*f40))/pi/(s*q2*r)**2
c       print *,fsigtv,fsigtvx
      endif
      end



      subroutine  polpar(s,q2,x,phi,ta,sps,spe,ccps,ccpe,de,apq,apn,
     .dk2ks,dapks,dksp1)
c some parametrization for polarization vectors
      implicit real*8(a-h,l,m,o-z)
      include 'const.inc'
      include 'pol.inc'

	sx=s-x
      alq=sx**2+4d0*amp2*q2
      sqls=sqrt(s**2-4d0*amp2*aml2)
      as=s/2./aml/sqls
      bs=0.
      cs=-aml/sqls
      ym=q2+2d0*aml**2

	 ael=amp/sqls
	 bel=0.
	 cel=-s/2d0/amp/sqls

	 sqn=sqrt(s*x*q2-alq*aml2-amp2*q2**2)
	 aet=(-s*x+2d0*amp2*ym)/sqls/sqn/2.
	 bet=sqls/sqn/2.
	 cet=-(s*q2+2d0*aml2*sx)/sqls/sqn/2.
      ae=cos(thetapn)*ael+sin(thetapn)*cos(phi-phipn)*aet
      be=cos(thetapn)*bel+sin(thetapn)*cos(phi-phipn)*bet
      ce=cos(thetapn)*cel+sin(thetapn)*cos(phi-phipn)*cet
	de=sin(thetapn)*sin(phi-phipn)
      apq=-q2*(ae-be)+ce*sx
      apn=(q2+4.*aml2)*(ae+be)+ce*(s+x)
      dk2ks=as*ym+2d0*aml2*bs+cs*x
      dapks=2.*(2d0*aml2*(as*ae+bs*be)+2d0*amp2*cs*ce+ym*(as*be+bs*ae)
     .+s*(as*ce+cs*ae)+x*(bs*ce+cs*be))
      dksp1=as*s+bs*x+cs*2d0*amp2
	 sps=as+bs
	 spe=ae+be
	 ccpe=(ae-be)*ta+2.*ce
	 ccps=(as-bs)*ta+2.*cs
      end

      subroutine sffun(t,f1,f2,f3,f4)
c unpolarized (f1, f2)  and polarized (f3, f4) hadronic structure functions
      implicit real*8(a-h,l,m,o-z)
      include 'const.inc'
       call ffpro(t,gep,gmp)
       tau=t/4d0/amp2
       f1=4d0*tau*amp2*gmp**2
       f2=4d0*amp2*(gep**2+tau*gmp**2)/(1d0+tau)
       f3=-2d0*amp2*gep*gmp
       f4=-amp2*gmp*(Gep-gmp)/(1d0+tau)
      end

      real*8 function fsib(s,q2,phi)
c Born cross section ds/dQ2/dphi
      implicit real*8(a-h,o-z)
      include 'const.inc'
      include 'pol.inc'
       call polpar(s,q2,s-q2,phi,0d0,sps,spe,ccps,ccpe,de,apq0,apn0,
     .dk2ks0,dapks0,dksp10)
C************ Beyond URA*************************************** 
c	tm1=q2
        tm1=q2-2*aml2
c**************************************************************
	tm2=(s*(s-q2)-amp2*q2)/2d0/amp2
        tm3=2.*(apq0*dk2ks0-dapks0*q2)*aml/amp
        tm4=-apq0*(dk2ks0-2.*dksp10)*q2*aml/amp**3
      call sffun(q2,f1,f2,f3,f4)
      fsib=alpha**2*barn/s**2/q2**2*(f1*tm1+f2*tm2
     .+plrun*pnrun*(tm3*f3+tm4*f4))
      end

      real*8 function fsib2(s,q2,phi)
c Born cross section ds/dQ2/dphi
      implicit real*8(a-h,o-z)
      include 'const.inc'
      include 'pol.inc'
       call polpar(s,q2,s-q2,phi,0d0,sps,spe,ccps,ccpe,de,apq0,apn0,
     .dk2ks0,dapks0,dksp10)
C************ Beyond URA*************************************** 
c	tm1=q2
        tm1=q2-2*aml2
c**************************************************************
	tm2=(s*(s-q2)-amp2*q2)/2d0/amp2
        tm3=2.*(apq0*dk2ks0-dapks0*q2)*aml/amp
        tm4=-apq0*(dk2ks0-2.*dksp10)*q2*aml/amp**3
      call sffun(q2,f1,f2,f3,f4)
      fsib2=alpha**2*barn/dsqrt(s*s-2*s*(aml2+amp2)+(aml2-amp2)*
     .(aml2-amp2))/q2**2*(f1*tm1+f2*tm2
     .+plrun*pnrun*(tm3*f3+tm4*f4))
      end

      subroutine grid_init
c grids for generation of photonic variable
      implicit real*8(a-h,o-z)
        include 'grid.inc'
      do it1=1,nt1/2
	 grt1(it1)=0.5*(2d0*dble(it1)/dble(nt1))**6
	 grt1(nt1/2+it1)=0.5*(2d0-(2d0*dble(nt1/2-it1)/dble(nt1))**5)
      enddo
      do it2=1,nt2
	 grt2(it2)=(dble(it2)/dble(nt2))**7
      enddo
      do ivv=1,30
	 grv0(ivv)=(1d0-(dble(30-ivv)/29d0)**5)
	 grv0(30+ivv)=(dble(ivv)/30d0)**9
      enddo
      do ivv=1,nvv
	 grv(ivv,0)=dble(ivv)/dble(nvv)
      enddo

      do ivv=1,nvv/2
	 grv0(ivv)=(1d0-(dble(nvv/2-ivv)/dble(nvv/2-1))**5)
	 grv0(nvv/2+ivv)=(2d0*dble(ivv)/dble(nvv))**9
      enddo
      do ivv=1,nvv
	 grv(ivv,0)=dble(ivv)/dble(nvv)
      enddo
      do ii=1,7
	 ii15=(ii-1)*5
	 do ivv=1,ii15
	   grv(ivv,ii)=dble(ivv)/40d0
	 enddo
	 if(ii.eq.1)then
	   ggrr=0d0
	 else
	   ggrr=grv(ii15,ii)
	 endif
	 do ivv=1,30
	   grv(ii15+ivv,ii)=ggrr+dble(ivv)/120d0
	 enddo
	 do ivv=ii15+11,40
	   grv(20+ivv,ii)=dble(ivv)/40d0
	 enddo
      enddo
      do iph=1,15
	 grphi(iph)=0.01d0*dble(iph)/15d0
      enddo
      do iph=16,nph
	 grphi(iph)=0.01d0+0.99d0*(dble(iph)-15d0)/(dble(nph)-15d0)
      enddo
      do iph=1,nph
	 grphi(iph)=(dble(iph)/dble(nph))**3
c	print*,iph,grphi(iph)
      enddo
c      write(22,'(7f9.3)')((grv(iiv,ii),ii=1,7),iiv=1,60)
      end


       subroutine ffpro(t,gep,gmp)
c electric and magnetic form factors of the proton

      implicit real*8(a-h,l,m,o-z)
      include 'const.inc'
      tau=t/4d0/amp2
      x1=tau
      x2=x1*tau
      x3=x2*tau
      x4=x3*tau
      x5=x4*tau
      gep=(1.0+2.90966*x1-1.11542229*x2+3.866171d-2*x3)/(1+14.5187212*
     .x1+40.88333*x2+99.999998*x3+4.579d-5*x4+10.3580447*x5)
      gmp=2.792782*(1.0-1.43573*x1+1.19052066*x2+2.5455841d-1*x3)/(1+
     .9.70703681*x1+3.7357d-4*x2+6.0d-8*x3+9.9527277*x4+12.7977739*x5)
c     gep=1.2742/(1.+t/0.6394**2)-.2742/(1.+t/1.582**2)
c     gmp=(1.3262/(1.+t/0.6397**2)-.3262/(1.+t/1.3137**2))*2.7921
c     gep=1./((1.+.61*t)*(1.+2.31*t)*(1.+.04*t))
c     gmp=amm*gep
      end





c$nodebug
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C *						      *
      FUNCTION URAND(IY)
C *						      *
C *   This is a standard pseudo-random generator      *
C *   that work on IBM-370 and IBM-PC. We don't       *
C *   know does it work on SUN? 		      *
C *						      *
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
      REAL*4  URAND,S
      INTEGER*4 IY
c      INTEGER*4 A,C,MASK
c      PARAMETER (A  = 843314861)
c      PARAMETER (C  = 453816693)
c      PARAMETER (S  = 4.6566128E-10)

c      IY=IAND(A*IY+C,Z'7FFFFFFF')
c      URAND=FLOAT(IY)*S
      call RANDOM_NUMBER(S)
      URAND=S
      END



****************** vacpol *************************************

      double precision function vacpol(t)
c contribution from vacuum polarization by leptons (suml) and hadrons (sumh)
      implicit real*8(a-h,l,m,o-z)
c      common/cmp/pi,alpha,amp,amp2,aml,aml2,barn
      include 'const.inc'
      dimension am2(3)
c
c    am2 : squared masses of charge leptons
c
      data am2/.26110d-6,.111637d-1,3.18301d0/

      suml=0.
      do 10 i=1,3
	 a2=2.*am2(i)
	 sqlmi=dsqrt(t*t+2.*a2*t)
	 allmi=dlog((sqlmi+t)/(sqlmi-t))/sqlmi
  10  suml=suml+2.*(t+a2)*allmi/3.-10./9.+4.*a2*(1.-a2*allmi)/3./t
      if(t.lt.1.d0)then
	aaa = -1.345d-9
	bbb = -2.302d-3
	ccc = 4.091
      elseif(t.lt.64d0)then
	aaa = -1.512d-3
	bbb =  -2.822d-3
	ccc = 1.218
      else
	aaa = -1.1344d-3
	bbb = -3.0680d-3
	ccc = 9.9992d-1
      endif
      sumh = -(aaa+bbb*log(1.+ccc*t)) *2*pi/alpha

      vacpol=suml+sumh

      end
****************** fspens *************************************

      double precision function fspens(x)
c
c    spence function
c
      implicit real*8(a-h,o-z)
      f=0.d0
      a=1.d0
      an=0.d0
      tch=1.d-16
  1   an=an+1.d0
      a=a*x
      b=a/an**2
      f=f+b
      if(b-tch)2,2,1
  2   fspens=f
      return
      end

      double precision function fspen(x)
      implicit real*8(a-h,o-z)
      data f1/1.644934d0/
      if(x)8,1,1
  1   if(x-.5d0)2,2,3
    2 fspen=fspens(x)
      return
    3 if(x-1d0)4,4,5
    4 fspen=f1-dlog(x)*dlog(1d0-x+1d-10)-fspens(1d0-x)
      return
    5 if(x-2d0)6,6,7
    6 fspen=f1-.5*dlog(x)*dlog((x-1d0)**2/x)+fspens(1d0-1d0/x)
      return
    7 fspen=2d0*f1-.5d0*dlog(x)**2-fspens(1d0/x)
      return
    8 if(x+1d0)10,9,9
   9  fspen=-.5d0*dlog(1d0-x)**2-fspens(x/(x-1d0))
      return
  10  fspen=-.5*dlog(1.-x)*dlog(x**2/(1d0-x))-f1+fspens(1d0/(1d0-x))
      return
      end

      subroutine simps(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
c simps
c a1,b1 -the limits of integration
c h1 -an initial step of integration
c reps1,aeps1 - relative and absolute precision of integration
c funct -a name of function subprogram for calculation of integrand +
c x - an argument of the integrand
c ai - the value of integral
c aih- the value of integral with the step of integration
c aiabs- the value of integral for module of the integrand
c this subrogram calculates the definite integral with the relative or
c absolute precision by simpson+s method with the automatical choice
c of the step of integration
c if aeps1    is very small(like 1.e-17),then calculation of integral
c with reps1,and if reps1 is very small (like 1.e-10),then calculation
c of integral with aeps1
c when aeps1=reps1=0. then calculation with the constant step h1
c
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1)
      s=dsign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1)
      aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

      subroutine simpxx(a,b,np,ep,func,res)
      implicit real*8 (a-h,o-z)
      external func
      step=(b-a)/np
      call simps(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end

c****************** S_phi FUNCTION (no approximation) *************************************
      subroutine sphibis(s1,l1,s2,s3,l3,ssphibis)
      implicit real*8(a-h,l,m,o-z)
      include 'const.inc'
      real*8 tt(4),delij(4,4),sj(4),logtitju(4,4),logtitjl(4,4)
      real*8 logtitu(4,4),logtitl(4,4)
      real*8 tto,tti,sln1,sq1,snl2,snl3

      do i=1,4
        sj(i)=0d0
        do j=1,4
        delij(i,j)=0d0
        logtitju(i,j)=0d0
        logtitjl(i,j)=0d0
        logtitu(i,j)=0d0
        logtitl(i,j)=0d0
       enddo
      enddo

      l2=s2**2-4*aml2*amp2
      sql1=dsqrt(l1)
      sql2=dsqrt(l2)
      sql3=dsqrt(l3)
      
c      print *,s1,l1
      tt(1)=2d0*amp2*(s1+sql1)/(s2-sql2)-s3-2d0*amp2-sql2
      tt(2)=2d0*amp2*(s1-sql1)/(s2-sql2)-s3-2d0*amp2-sql2
      tt(3)=2d0*amp2*(-s1+sql1)/(s2+sql2)+s3+2d0*amp2-sql2
      tt(4)=2d0*amp2*(-s1-sql1)/(s2+sql2)+s3+2d0*amp2-sql2
      tl=-sql2+(s2*s3+2d0*amp2*(s2-s1))/sql2
      tu=sql3-sql2

      tto=(tu-tt(1))*(tu-tt(3))/(tu-tt(2))/(tu-tt(4))
      tti=(s2-sql2)/(s2+sql2)
      sln1=dlog(tti)*dlog(tto)
      sq1=s1/(2d0*sql1)

c      print *,sln1,sq1
      snl2=0.0d0
      snl3=0.0d0
 
      do i=1,4
        do j=1,4
        if(j.le.2.0) then
        sj(j)=1d0
        endif
        if(j.gt.2.0) then
        sj(j)=-1d0
        endif
        enddo
      enddo

      do i=1,4
        do j=1,4
        if(i.EQ.j)then
        logtitju(i,j)=0d0
        logtitjl(i,j)=0d0
        logtitu(i,j)=.5d0*dlog(abs(tu-tt(i)))*dlog(abs(tu-tt(i)))
        logtitl(i,j)=.5d0*dlog(abs(tl-tt(i)))*dlog(abs(tl-tt(i)))
        else
        logtitju(i,j)=dlog(abs(tt(i)-tt(j)))*dlog(abs(tu-tt(j)))
     .  -fspen((tu-tt(i))/(tt(j)-tt(i)))
        logtitjl(i,j)=dlog(abs(tt(i)-tt(j)))*dlog(abs(tl-tt(j)))
     .  -fspen((tl-tt(i))/(tt(j)-tt(i)))
        logtitu(i,j)=0d0
        logtitl(i,j)=0d0
     
        endif
c        print *,logtitju(i,j),logtitu(i,j),logtitjl(i,j),logtitl(i,j)
        enddo
      enddo
     

      do i=1,4
        do j=1,4
        snl2=snl2+(-1d0)**(i+1)*sj(j)*(logtitu(i,j)+logtitju(i,j))
        snl3=snl3+(-1d0)**(i+1)*sj(j)*(logtitl(i,j)+logtitjl(i,j))
        enddo
c       print *,snl2,snl3,snl2-snl3
      enddo
c      print *,snl2,snl3

      ssphibis=(snl1+snl2-snl3)*sq1

c      print *,ssphibis
      end


c****************** S_phi FUNCTION (no approximation) *************************************
      subroutine sphi(sss,lamb,a,b,ssphi)
      implicit real*8(a-h,l,m,o-z)
      include 'const.inc'
      real*8 sj(4),aj(4),tau(4),gam12(2),gam(2,4)
      real*8 slam,sqlam,d,sqd,sqb,gamu,q22,logms,logsq2s
      real*8 f1u,f2u,f11,f21
      integer*4 i,j,k
c      print *,sss,lamb,a,b
      do i=1,4
        sj(i)=0d0
        aj(i)=0d0
        tau(i)=0d0
      enddo

      do i=1,2
        gam12(i)=0d0
      enddo

      do i=1,2
       do j=1,4
       gam(i,j)=0d0
       enddo
      enddo

      sj(1)=1d0
      sj(2)=1d0
      sj(3)=-1d0
      sj(4)=-1d0
      sqlam=dsqrt(lamb)
      d=((lamb+b)*(lamb+b)/4d0)+(sss+a)*(lamb*a-sss*b)      
      sqd=dsqrt(d)
      sqb=dsqrt(b)
      gamu=(dsqrt(lamb+b)-sqb)/sqlam
c      print *,gamu

      do j=1,4
      aj(j)=sss-sj(j)*sqlam
      tau(j)=-a*sqlam+0.5*sj(j)*(b-lamb)+((-1)**(j))*sqd
c      print *,aj(j),tau(j)
      enddo

      gam12(1)=-((sqb-sqlam)**2)/(b-lamb)
      gam12(2)=((sqb+sqlam)**2)/(b-lamb)
 
      do i=1,2
        do j=1,4
        gam(i,j)=-(aj(j)*sqb+((-1)**(i+1))*dsqrt(b*aj(j)*aj(j)
     . +tau(j)*tau(j)))/tau(j)
        enddo
        
      enddo
      
      slam=sss/(2d0*sqlam)
      f1u=0.0d0
      f2u=0.0d0
      f11=0.0d0
      f21=0.0d0

      do i=1,2
       do k=1,2
        do j=1,4
       f2u=f2u+((-1)**i)*sj(j)*
     . fspen((gamu+(-1)**i)/(gam(k,j)+(-1)**i))
       f11=f11+((-1)**i)*sj(j)*
     . fspen((gam12(i)-gam12(1))/(gam12(i)-gam(k,j)))
       f21=f21+((-1)**i)*sj(j)*
     . fspen((gam12(1)+(-1)**i)/(gam(k,j)+(-1)**i))
       f1u=f1u+((-1)**i)*sj(j)*
     . fspen((gam12(i)-gamu)/(gam12(i)-gam(k,j)))
         enddo
       enddo
      enddo

      ssphi=slam*(f1u+f2u-f11-f21)
 
      end

