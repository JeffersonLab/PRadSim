      subroutine dsigma(the,q2,w,cscm,phicm,opt1,opt2,opt3
     *,sig0,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)

      implicit none

      real the,ki_mag,q2,w,cscm,phicm,kf_mag,s2
      real fkt
      real sig0,sigu,sigt,sigl,sigi,sigip,asym_p
      real nu,eps,eps1
      integer ehel
      real cthe     
      integer opt1,opt2,opt3
c variables for dvmp
      real q2_dummy, W_dummy, hcs_dummy, W2_dummy
      real Mp, mpi0, amp2, amu2, mesonmass, meta, mpip
      real xb_dummy, E_pi_cm, ppi_mag_cm, qv_mag_cm
      real nu_cm, t_dummy
      parameter (Mp=0.93827)
      PARAMETER (MPI0=0.134976,MPIP=0.13957018,META=0.547853)
      
      logical test1,test2,test3
c
      if (opt2.eq.1) then
        mesonmass = mpi0
      elseif (opt2.eq.3) then
        mesonmass = mpip
      elseif (opt2.eq.5) then
        mesonmass = meta
      endif

      nu 	= 0.5*(w**2 + q2 - Mp**2)/Mp
      s2	= sin(0.5*the)**2
      ki_mag	= (nu+sqrt(q2/s2+nu**2))*0.5

      kf_mag 	= ki_mag-nu
      eps 	= 1. /(1+2.0*(1+nu*nu/q2)*tan(0.5*the)**2)

c      print *, W_dummy, xb_dummy, hcs_dummy, t_dummy, amp2, amu2
c      test1	= ki_mag.lt.0.1.or.kf_mag.lt.0.01
c      test3	= opt1.ge.4.and.(w.lt.1.1.or.w.gt.1.7)
c      
      if (kf_mag.lt.0.1) then
        print *, 'low electron momentum',kf_mag,q2,w        
        sig0 = 0.
        sigu = 0.
        sigt = 0.
        sigl = 0.
        sigi = 0.
        sigip = 0.
        asym_p = 0.
c        return
      endif
c
c      if (opt1.gt.10.and.w.lt.1.7) then
c        print *, 'DVMP Low W',w,q2        
c        sig0 = 0.
c        sigu = 0.
c        sigt = 0.
c        sigl = 0.
c        sigi = 0.
c        sigip = 0.
c        asym_p = 0.
c        return
c      endif
c
      if(opt1.eq.1) call aao(q2,w,eps,cscm,phicm,1,sig0,
     1             sigu,sigt,sigl,sigi)
      if(opt1.eq.2) 
     1   call daresbury(q2,w,eps,cscm,phicm,1,sig0)
      if(opt1.eq.4.or.opt1.eq.5) then
c
        if(w.gt.1.82 .and. opt2.eq.1) then

c
        amp2 = 0.93827**2
        amu2 = mesonmass**2
        q2_dummy = q2
        W_dummy  = w
        hcs_dummy = cscm
        W2_dummy = W_dummy**2
        xb_dummy = q2_dummy/(W2_dummy - amp2 + q2_dummy)
        E_pi_cm    = 0.5*(w2_dummy+amu2-amp2)/W_dummy
        ppi_mag_cm = E_pi_cm**2 - amu2
        ppi_mag_cm = sqrt(ppi_mag_cm)
        qv_mag_cm  = ((W2_dummy+q2_dummy+amp2)/2.0/W_dummy)**2
     1                   -amp2
        qv_mag_cm  = sqrt(qv_mag_cm)
        nu_cm  = (W2_dummy-amp2-q2_dummy)/(2*W_dummy)
        t_dummy =  -(-q2_dummy) - amu2 + 2*nu_cm*E_pi_cm
     1             - hcs_dummy * (2*qv_mag_cm*ppi_mag_cm)
        call dvmpw(hcs_dummy,w_dummy,q2_dummy,phicm,ki_mag,ehel,mesonmass,
     1           sig0,sigu,sigl,sigt,sigi,sigip) 
c	    call dvmp(-t_dummy,xb_dummy,q2_dummy,phicm,ki_mag,0,
c     *   sig0,sigu,sigl,sigt,sigi,sigip,asym_p,ehel)
c         print *, 'DSIGMA: ',sig0,sigu,sigl,sigt,sigi,sigip,asym_p,ehel 
c        print *, '**************************************************'
        fkt  = 2*W_dummy*ppi_mag_cm/(W2_dummy - amp2)
        sig0 = fkt*sig0
      else
         call maid_lee(q2,w,eps,cscm,phicm,opt1,opt2,opt3,
     1        sig0,sigu,sigt,sigl,sigi,sigip,asym_p,ehel)
      endif


c     print *, 'dsigma.F:',sig0, fkt, opt1
      endif
c     
      end
