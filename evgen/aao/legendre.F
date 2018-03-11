c---------------------------------------------------------------------------
c       Calculatess the first and second derivatives of 
c       Legendre Polynomials.
c---------------------------------------------------------------------------
      subroutine legendre

      implicit none

      include 'mpintp.inc'
      include 'spp.inc'

      real X
      integer l
c
      X = csthcm
c
      pol(0,1) = 0.
      pol(1,1) = 1.
      pol(2,1) = 3.*X 
      pol(3,1) = (15.*X*X-3.)/2.
      pol(4,1) = (35.*X*X*X-15.*X)/2.
      pol(5,1) = (315.*X*X*X*X-210.*X*X+15.)/8.
      pol(6,1) = (693.*X*X*X*X*X-630.*X*X*X +105.*X)/8.
      pol(7,1) = (3003.*X*X*X*X*X*X-3465.*X*X*X*X+945*X*X-35.)/16.
      pol(0,2) = 0.
      pol(1,2) = 0.
      pol(2,2) = 3.
      pol(3,2) = 15.*X 
      pol(4,2) = (105.*X*X-15.)/2.
      pol(5,2) = (315.*X*X*X-105.*X)/2.
      pol(6,2) = (3465.*X*X*X*X-1890.*X*X+105.)/8.  
      pol(7,2) = (9009.*X*X*X*X*X-6930.*X*X*X +945.*X)/8.
c
c      do l = 0, 6
c          print *, 'legedre', 'l=',l
c          print *, pol(l,1), pol(l,2) 
c      enddo
      return
      end

