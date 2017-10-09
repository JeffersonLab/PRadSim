      subroutine helicity_amps

      implicit none

      include 'mpintp.inc'
      include 'spp.inc'

      real   root2,s,s2,c,c2
      real   theta_cm
      integer l
      
      data root2/1.41421/

      call cgln_amps

      theta_cm = acos(csthcm)
      
      s   = sin(theta_cm)
      c   = cos(theta_cm)
      s2  = sin(theta_cm/2.0)
      c2  = cos(theta_cm/2.0)
      
      hh1 = 0.0
      hh2 = 0.0
      hh3 = 0.0
      hh4 = 0.0
      hh5 = 0.0
      hh6 = 0.0
      
      if(method_helicity.eq.1) then
c        hh1 = -s*(ff3 + ff4*c)/root2
c        hh2 = -(2*ff1 - 2*ff2*c + ff4*s2)/root2
c        hh3 = -ff4*s2/root2
c        hh4 = s*(2*ff2 + ff3 + ff4*c)/root2
c        hh5 = ff5 + ff6*c
c        hh6 = ff6*s
        hh1 = -s*c2*(ff3 + ff4)/root2
        hh2 =  c2*((ff2 - ff1) + 0.5*(1 - c)*(ff3 - ff4))*root2
        hh3 =  s*s2*(ff3 - ff4)/root2
        hh4 =  s2*((ff2 + ff1) + 0.5*(1 + c)*(ff3 + ff4))*root2
        hh5 =  c2*(ff5 + ff6)
        hh6 =  s2*(ff6 - ff5)
        
      else

        do l = 0, wave_L
          hh1 = hh1 + s*c2*(bp(l)-bm(l+1))*(pol(l,2)-pol(l+1,2))/root2
          hh2 = hh2 + c2*(ap(l)-am(l+1))*(pol(l,1)-pol(l+1,1))*root2
          hh3 = hh3 + s*s2*(bp(l)+bm(l+1))*(pol(l,2)+pol(l+1,2))/root2
          hh4 = hh4 + s2*(ap(l)+am(l+1))*(pol(l,1)+pol(l+1,1))*root2
          hh5 = hh5 - c2*(cp(l)-cm(l+1))*(pol(l,1)-pol(l+1,1))
          hh6 = hh6 - s2*(cp(l)+cm(l+1))*(pol(l,1)+pol(l+1,1))
        enddo
        
      endif
      
      end



