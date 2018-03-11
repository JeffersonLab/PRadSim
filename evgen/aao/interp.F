C------------------------------------------------------------------------------
C INTERP
C
C Does interpolartion of database given two variables V1 and V2

      subroutine interp(v1,v2,vec_f)
C
      implicit none

      include 'mpintp.inc'

C Passed variables:
               
      real  vec_f(max_mp)
      real  v1,v2
               
C Local variables:
      real X(2)  
      real fint
      integer mpp
C
      if (v1.gt.var1_max) then
c	write(6,'(a,2(1x,e15.8))')' FATAL ERROR: var1 .gt. max - ',v1,var1_max
c	STOP

      endif

      if (v1.lt.var1_min) then
c	write(6,'(a,2(1x,e15.8))')' FATAL ERROR: var1 .lt. min - ',v1,var1_min
c	STOP
	     goto 1
      endif

      if (v2.gt.var2_max) then
c	write(6,'(a,2(1x,e15.8))')' FATAL ERROR: var2 .gt. max - ',v2,var2_max
c	STOP
         goto 1
      endif

      if (v2.lt.var2_min) then
c	write(6,'(a,2(1x,e15.8))')' FATAL ERROR: var2 .lt. min - ',v2,var2_min
c	STOP
	     goto 1
      endif

      if (method_spline.eq.1) then
C Get the 62 multipoles (F's) at the grid point (V1,V2)
C by 2-dimensional natural cubic spline interpolation.

        call splin2(var1,var2,sf1 ,d2sf1 ,nvar1,nvar2,v1,v2,vec_f(1) )
        call splin2(var1,var2,sf2 ,d2sf2 ,nvar1,nvar2,v1,v2,vec_f(2) )
        call splin2(var1,var2,sf3 ,d2sf3 ,nvar1,nvar2,v1,v2,vec_f(3) )
        call splin2(var1,var2,sf4 ,d2sf4 ,nvar1,nvar2,v1,v2,vec_f(4) )
        call splin2(var1,var2,sf5 ,d2sf5 ,nvar1,nvar2,v1,v2,vec_f(5) )
        call splin2(var1,var2,sf6 ,d2sf6 ,nvar1,nvar2,v1,v2,vec_f(6) )
        call splin2(var1,var2,sf7 ,d2sf7 ,nvar1,nvar2,v1,v2,vec_f(7) )
        call splin2(var1,var2,sf8 ,d2sf8 ,nvar1,nvar2,v1,v2,vec_f(8) )
        call splin2(var1,var2,sf9 ,d2sf9 ,nvar1,nvar2,v1,v2,vec_f(9) )
        call splin2(var1,var2,sf10,d2sf10,nvar1,nvar2,v1,v2,vec_f(10))
        call splin2(var1,var2,sf11,d2sf11,nvar1,nvar2,v1,v2,vec_f(11))
        call splin2(var1,var2,sf12,d2sf12,nvar1,nvar2,v1,v2,vec_f(12))
        call splin2(var1,var2,sf13,d2sf13,nvar1,nvar2,v1,v2,vec_f(13))
        call splin2(var1,var2,sf14,d2sf14,nvar1,nvar2,v1,v2,vec_f(14))
        call splin2(var1,var2,sf15,d2sf15,nvar1,nvar2,v1,v2,vec_f(15))
        call splin2(var1,var2,sf16,d2sf16,nvar1,nvar2,v1,v2,vec_f(16))
        call splin2(var1,var2,sf17,d2sf17,nvar1,nvar2,v1,v2,vec_f(17))
        call splin2(var1,var2,sf18,d2sf18,nvar1,nvar2,v1,v2,vec_f(18))
        call splin2(var1,var2,sf19,d2sf19,nvar1,nvar2,v1,v2,vec_f(19))
        call splin2(var1,var2,sf20,d2sf20,nvar1,nvar2,v1,v2,vec_f(20))
        call splin2(var1,var2,sf21,d2sf21,nvar1,nvar2,v1,v2,vec_f(21))
        call splin2(var1,var2,sf22,d2sf22,nvar1,nvar2,v1,v2,vec_f(22))
        call splin2(var1,var2,sf23,d2sf23,nvar1,nvar2,v1,v2,vec_f(23))
        call splin2(var1,var2,sf24,d2sf24,nvar1,nvar2,v1,v2,vec_f(24))
        call splin2(var1,var2,sf25,d2sf25,nvar1,nvar2,v1,v2,vec_f(25))
        call splin2(var1,var2,sf26,d2sf26,nvar1,nvar2,v1,v2,vec_f(26))
        call splin2(var1,var2,sf27,d2sf27,nvar1,nvar2,v1,v2,vec_f(27))
        call splin2(var1,var2,sf28,d2sf28,nvar1,nvar2,v1,v2,vec_f(28))
        call splin2(var1,var2,sf29,d2sf29,nvar1,nvar2,v1,v2,vec_f(29))
        call splin2(var1,var2,sf30,d2sf30,nvar1,nvar2,v1,v2,vec_f(30))
        call splin2(var1,var2,sf31,d2sf31,nvar1,nvar2,v1,v2,vec_f(31))
        call splin2(var1,var2,sf32,d2sf32,nvar1,nvar2,v1,v2,vec_f(32))
        call splin2(var1,var2,sf33,d2sf33,nvar1,nvar2,v1,v2,vec_f(33))
        call splin2(var1,var2,sf34,d2sf34,nvar1,nvar2,v1,v2,vec_f(34))
        call splin2(var1,var2,sf35,d2sf35,nvar1,nvar2,v1,v2,vec_f(35))
        call splin2(var1,var2,sf36,d2sf36,nvar1,nvar2,v1,v2,vec_f(36))
        call splin2(var1,var2,sf37,d2sf37,nvar1,nvar2,v1,v2,vec_f(37))
        call splin2(var1,var2,sf38,d2sf38,nvar1,nvar2,v1,v2,vec_f(38))
        call splin2(var1,var2,sf39,d2sf39,nvar1,nvar2,v1,v2,vec_f(39))
        call splin2(var1,var2,sf40,d2sf40,nvar1,nvar2,v1,v2,vec_f(40))
        call splin2(var1,var2,sf41,d2sf41,nvar1,nvar2,v1,v2,vec_f(41))
        call splin2(var1,var2,sf42,d2sf42,nvar1,nvar2,v1,v2,vec_f(42))
        call splin2(var1,var2,sf43,d2sf43,nvar1,nvar2,v1,v2,vec_f(43))
        call splin2(var1,var2,sf44,d2sf44,nvar1,nvar2,v1,v2,vec_f(44))
        call splin2(var1,var2,sf45,d2sf45,nvar1,nvar2,v1,v2,vec_f(45))
        call splin2(var1,var2,sf46,d2sf46,nvar1,nvar2,v1,v2,vec_f(46))
        call splin2(var1,var2,sf47,d2sf47,nvar1,nvar2,v1,v2,vec_f(47))
        call splin2(var1,var2,sf48,d2sf48,nvar1,nvar2,v1,v2,vec_f(48))
        call splin2(var1,var2,sf49,d2sf49,nvar1,nvar2,v1,v2,vec_f(49))
        call splin2(var1,var2,sf50,d2sf50,nvar1,nvar2,v1,v2,vec_f(50))
        call splin2(var1,var2,sf51,d2sf51,nvar1,nvar2,v1,v2,vec_f(51))
        call splin2(var1,var2,sf52,d2sf52,nvar1,nvar2,v1,v2,vec_f(52))
        call splin2(var1,var2,sf53,d2sf53,nvar1,nvar2,v1,v2,vec_f(53))
        call splin2(var1,var2,sf54,d2sf54,nvar1,nvar2,v1,v2,vec_f(54))
        call splin2(var1,var2,sf55,d2sf55,nvar1,nvar2,v1,v2,vec_f(55))
        call splin2(var1,var2,sf56,d2sf56,nvar1,nvar2,v1,v2,vec_f(56))
        call splin2(var1,var2,sf57,d2sf57,nvar1,nvar2,v1,v2,vec_f(57))
        call splin2(var1,var2,sf58,d2sf58,nvar1,nvar2,v1,v2,vec_f(58))
        call splin2(var1,var2,sf59,d2sf59,nvar1,nvar2,v1,v2,vec_f(59))
        call splin2(var1,var2,sf60,d2sf60,nvar1,nvar2,v1,v2,vec_f(60))
        call splin2(var1,var2,sf61,d2sf61,nvar1,nvar2,v1,v2,vec_f(61))
        call splin2(var1,var2,sf62,d2sf62,nvar1,nvar2,v1,v2,vec_f(62))
C
      elseif (method_spline.eq.2) then
C Get the 62 multipoles (F's) at the grid point (V1,V2)
C by 2-dimensional natural cubic spline interpolation.
C
        X(1) = v1
        X(2) = v2
C
        vec_f(1 ) = fint(2,X,NA,VAR,sf1 )
        vec_f(2 ) = fint(2,X,NA,VAR,sf2 )
        vec_f(3 ) = fint(2,X,NA,VAR,sf3 )
        vec_f(4 ) = fint(2,X,NA,VAR,sf4 )
        vec_f(5 ) = fint(2,X,NA,VAR,sf5 )
        vec_f(6 ) = fint(2,X,NA,VAR,sf6 )
        vec_f(7 ) = fint(2,X,NA,VAR,sf7 )
        vec_f(8 ) = fint(2,X,NA,VAR,sf8 )
        vec_f(9 ) = fint(2,X,NA,VAR,sf9 )
        vec_f(10) = fint(2,X,NA,VAR,sf10)
        vec_f(11) = fint(2,X,NA,VAR,sf11)
        vec_f(12) = fint(2,X,NA,VAR,sf12)
        vec_f(13) = fint(2,X,NA,VAR,sf13)
        vec_f(14) = fint(2,X,NA,VAR,sf14)
        vec_f(15) = fint(2,X,NA,VAR,sf15)
        vec_f(16) = fint(2,X,NA,VAR,sf16)
        vec_f(17) = fint(2,X,NA,VAR,sf17)
        vec_f(18) = fint(2,X,NA,VAR,sf18)
        vec_f(19) = fint(2,X,NA,VAR,sf19)
        vec_f(20) = fint(2,X,NA,VAR,sf20)
        vec_f(21) = fint(2,X,NA,VAR,sf21)
        vec_f(22) = fint(2,X,NA,VAR,sf22)
        vec_f(23) = fint(2,X,NA,VAR,sf23)
        vec_f(24) = fint(2,X,NA,VAR,sf24)
        vec_f(25) = fint(2,X,NA,VAR,sf25)
        vec_f(26) = fint(2,X,NA,VAR,sf26)
        vec_f(27) = fint(2,X,NA,VAR,sf27)
        vec_f(28) = fint(2,X,NA,VAR,sf28)
        vec_f(29) = fint(2,X,NA,VAR,sf29)
        vec_f(30) = fint(2,X,NA,VAR,sf30)
        vec_f(31) = fint(2,X,NA,VAR,sf31)
        vec_f(32) = fint(2,X,NA,VAR,sf32)
        vec_f(33) = fint(2,X,NA,VAR,sf33)
        vec_f(34) = fint(2,X,NA,VAR,sf34)
        vec_f(35) = fint(2,X,NA,VAR,sf35)
        vec_f(36) = fint(2,X,NA,VAR,sf36)
        vec_f(37) = fint(2,X,NA,VAR,sf37)
        vec_f(38) = fint(2,X,NA,VAR,sf38)
        vec_f(39) = fint(2,X,NA,VAR,sf39)
        vec_f(40) = fint(2,X,NA,VAR,sf40)
        vec_f(41) = fint(2,X,NA,VAR,sf41)
        vec_f(42) = fint(2,X,NA,VAR,sf42)
        vec_f(43) = fint(2,X,NA,VAR,sf43)
        vec_f(44) = fint(2,X,NA,VAR,sf44)
        vec_f(45) = fint(2,X,NA,VAR,sf45)
        vec_f(46) = fint(2,X,NA,VAR,sf46)
        vec_f(47) = fint(2,X,NA,VAR,sf47)
        vec_f(48) = fint(2,X,NA,VAR,sf48)
        vec_f(49) = fint(2,X,NA,VAR,sf49)
        vec_f(50) = fint(2,X,NA,VAR,sf50)
        vec_f(51) = fint(2,X,NA,VAR,sf51)
        vec_f(52) = fint(2,X,NA,VAR,sf52)
        vec_f(53) = fint(2,X,NA,VAR,sf53)
        vec_f(54) = fint(2,X,NA,VAR,sf54)
        vec_f(55) = fint(2,X,NA,VAR,sf55)
        vec_f(56) = fint(2,X,NA,VAR,sf56)
        vec_f(57) = fint(2,X,NA,VAR,sf57)
        vec_f(58) = fint(2,X,NA,VAR,sf58)
        vec_f(59) = fint(2,X,NA,VAR,sf59)
        vec_f(60) = fint(2,X,NA,VAR,sf60)
        vec_f(61) = fint(2,X,NA,VAR,sf61)
        vec_f(62) = fint(2,X,NA,VAR,sf62)
      endif
C
      return
C
 1    do mpp = 1,max_mp
         vec_f(mpp) = 0.0
      enddo
      print *, 'v1= ',v1,' Q2 min= ',var1_min,' Q2 max= ',var1_max
      print *, 'v2= ',v2,' W  min= ',var2_min,' W  max= ',var2_max
      stop
C
      return
      end
