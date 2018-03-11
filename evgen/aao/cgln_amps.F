        subroutine cgln_amps
      
        implicit none

        include 'mpintp.inc'
        include 'spp.inc'

        integer l

        call multipole_amps

        call legendre
        
        if (method_helicity.eq.1) then

        ff1 = 0.0
        ff2 = 0.0
        ff3 = 0.0
        ff4 = 0.0
        ff5 = 0.0
        ff6 = 0.0
        
        do l = 0, wave_L

          if(l.lt.2) then
             ff1 = ff1 + (    l*mp(l) +  ep(l))*pol(l+1,1)
          else
             ff1 = ff1 + (    l*mp(l) +  ep(l))*pol(l+1,1)
     1                 + ((l+1)*mm(l) +  em(l))*pol(l-1,1)
          endif
        enddo

        do l = 1, wave_L
          ff2 = ff2 + ((l+1)*mp(l) + l*mm(l))*pol(l,1)
          if(l.lt.2) then
             ff3 = ff3 + (      ep(l) -   mp(l))*pol(l+1,2)
          else
             ff3 = ff3 + (      ep(l) -   mp(l))*pol(l+1,2) +
     1                   (      em(l) +   mm(l))*pol(l-1,2)
          endif
        enddo

        do l = 2, wave_L
          ff4 = ff4 + (mp(l) - ep(l) -mm(l) - em(l))* pol(l,2)
        enddo

        do l = 0, wave_L
          if(l.lt.2) then
             ff5 = ff5 + (l+1)*sp(l)*pol(l+1,1)
          else
             ff5 = ff5 + ((l+1)*sp(l)*pol(l+1,1) - l*sm(l)*pol(l-1,1))
          endif
        enddo

        do l = 1, wave_L
          ff6 = ff6 + (l*sm(l) - (l+1)*sp(l))*pol(l,1)
        enddo
        
        else
        do l = 0, max_L+1
          ap(l)  = 0.0
          am(l)  = 0.0
          bp(l)  = 0.0
          bm(l)  = 0.0
          cp(l)  = 0.0
          cm(l)  = 0.0
        enddo

        do l = 0, max_L
          ap(l)  = 0.5*(l*mp(l) + (l+2)*ep(l))
          am(l+1)= 0.5*((l+2)*mm(l+1) - l*em(l+1))
          bp(l)  = ep(l)   - mp(l)
          bm(l+1)= em(l+1) + mm(l+1)
          cp(l)  =  (l+1)*sp(l)
          cm(l+1)= -(l+1)*sm(l+1)
        enddo
        
        endif

      return
      end
