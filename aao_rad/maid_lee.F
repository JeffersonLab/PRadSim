      subroutine maid_lee(a1,a2,a3,a4,a5,theory_opt,channel_opt,
     1           resonance_opt,sigma0,sigu,sigt,sigl,sigi,sigip,
     1           asym_p,ehel)
c
      implicit none
c      
      real a1,a2,a3,a4,a5
c
      include 'mpintp.inc'
      include 'spp.inc'
c
      integer theory_opt, channel_opt, resonance_opt, n_call
      real E_pi_cm
      real m_pi, alpi
      real sigma0,sigu,sigt,sigl,sigi,sigip,asym_p
      integer ehel
c
      data n_call /0/
c
      m_pi  = m_pip
      if (channel_opt.eq.1) m_pi  = m_pi0
c
      q2          = a1
      w           = a2
      epsilon     = a3
      csthcm      = a4
      phicm       = a5
c
      if (n_call.ne.0) goto 100
c
      method_helicity = 1
c
      if (theory_opt.eq.3) then
        call revinm('CLAS_PARMS',
     1  'spp_tbl/hlee_pi0.tbl',data_file)
      endif      

      if (theory_opt.eq.4) then
        if (channel_opt.eq.1) then
          if (resonance_opt.eq.0) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.1) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_p33_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.2) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_p11_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.3) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_d13_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.4) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_s11_1535_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.5) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_s11_1635_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.6) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_f15_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.7) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_d33_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.8) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_born_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.9) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_rho_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.10) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_omega_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.11) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_p11_off-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.12) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_d13_off-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.13) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_s11_1535_off-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.14) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_all_res_off-PPpi.tbl',data_file)
          endif
        endif
        if (channel_opt.eq.2) then
          if (resonance_opt.eq.0) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.1) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_p33_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.2) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_p11_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.3) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_d13_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.4) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_s11_1535_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.5) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_s11_1635_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.6) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_f15_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.7) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_d33_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.8) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_born_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.9) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_rho_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.10) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_omega_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.11) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_p11_off-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.12) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_d13_off-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.13) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_s11_1535_off-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.14) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_all_res_off-NPpi.tbl',data_file)
          endif
        endif
        if (channel_opt.eq.3) then
          if (resonance_opt.eq.0) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.1) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_p33_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.2) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_p11_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.3) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_d13_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.4) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_s11_1535_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.5) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_s11_1635_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.6) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_f15_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.7) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_d33_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.8) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_born_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.9) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_rho_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.10) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_omega_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.11) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_p11_off-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.12) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_d13_off-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.13) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_s11_1535_off-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.14) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid98_all_res_off-PNpi.tbl',data_file)
          endif
        endif
      endif
c
      if (theory_opt.eq.5) then
        if (channel_opt.eq.1) then
          if (resonance_opt.eq.0) then
            data_file = '/group/clas/parms/spp_tbl/maid07-PPpi.tbl';
c            call revinm('CLAS_PARMS',
c     1      'spp_tbl/maid2000-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.1) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_p33_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.2) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_p11_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.3) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_d13_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.4) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_s11_1535_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.5) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_s11_1635_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.6) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_f15_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.7) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_d33_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.8) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_born_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.9) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_rho_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.10) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_omega_on-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.11) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_p11_off-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.12) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_d13_off-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.13) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_s11_1535_off-PPpi.tbl',data_file)
          elseif (resonance_opt.eq.14) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_all_res_off-PPpi.tbl',data_file)
          endif
        endif
        if (channel_opt.eq.2) then
          if (resonance_opt.eq.0) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.1) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_p33_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.2) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_p11_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.3) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_d13_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.4) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_s11_1535_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.5) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_s11_1635_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.6) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_f15_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.7) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_d33_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.8) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_born_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.9) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_rho_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.10) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_omega_on-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.11) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_p11_off-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.12) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_d13_off-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.13) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_s11_1535_off-NPpi.tbl',data_file)
          elseif (resonance_opt.eq.14) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_all_res_off-NPpi.tbl',data_file)
          endif
        endif
        if (channel_opt.eq.3) then
          if (resonance_opt.eq.0) then
            data_file = '/group/clas/parms/spp_tbl/maid07-PNpi.tbl';
c            call revinm('CLAS_PARMS',
c     1      'spp_tbl/maid2000-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.1) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_p33_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.2) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_p11_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.3) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_d13_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.4) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_s11_1535_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.5) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_s11_1635_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.6) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_f15_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.7) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_d33_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.8) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_born_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.9) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_rho_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.10) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_omega_on-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.11) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_p11_off-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.12) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_d13_off-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.13) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_s11_1535_off-PNpi.tbl',data_file)
          elseif (resonance_opt.eq.14) then
            call revinm('CLAS_PARMS',
     1      'spp_tbl/maid2000_all_res_off-PNpi.tbl',data_file)
          endif
        endif
      endif
c
c      write(6,*) 'Enter max pi-N angular momentum (0-5): '
c
c      read(5,*) mwave_L
      mwave_L = 5   
c  
      write(6,*) 'Reading multipoles from ',data_file
c
      call read_sf_file(data_file,2)
c
      n_call = n_call + 1
 
 100  E_pi_cm    = 0.5*(W*W+m_pi**2-m_p**2)/W
      ppi_mag_cm = E_pi_cm**2 - m_pi**2
      ppi_mag_cm = sqrt(ppi_mag_cm)
      qv_mag_cm  = ((W*W+Q2+m_p**2)/2.0/W)**2-m_p**2
      qv_mag_cm  = sqrt(qv_mag_cm)
c      nu_cm      = sqrt(qv_mag_cm**2-Q2)
c     this is the  right calculation for nu_cm
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      nu_cm  = (W*W-m_p**2-Q2)/(2*W)    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      fkt = 2.0*W*ppi_mag_cm/(W**2 - m_p**2)   
      e_hel = ehel
c   
      call xsection
c
      sigma0 = sigma_0
      sigu   = sigma_t
      sigl   = sigma_l
      sigt   = sigma_tt
      sigi   = sigma_lt
      sigip  = sigma_ltp
      if (e_hel.gt.-0.5.and.e_hel.lt.0.5) then
        asym_p = 0.0
      else
        asym_p = asym_ltp
      endif
      return
      end
      subroutine revinm(A,B,C)
      Character A*(*),B*(*),C*(*)
      return
      end
