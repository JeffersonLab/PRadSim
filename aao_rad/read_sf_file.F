C-------------------------------------------------------------------------------
C READ_DATA_FILE
C
C Reads data file for use in interpolation
C
CEV Both these routines operate in double precision arithmetic
C   to suppress round off errors in the interpolation procedure.
C   But INTERP returns single precision results.

      subroutine read_sf_file(rfile,iunit)

      implicit none

      include 'mpintp.inc'
      include 'spp.inc'
C Passed variables:

      character rfile*(*)
      integer iunit

C Local variables:

      integer jvar1,jvar2
      real dumvar1,dumvar2,dumvar3,dumvar4
      character*80 dummy
      character*2 dummyX
      character*5 dummyY
      real var1tmp1,var1tmp2
      real var2tmp1,var2tmp2

c      write(6,*) 'Input method of spline fitting (1=cubic, 2=linear)'
c      read(5,*) method_spline 
      method_spline = 2

C First open the file.
      open(unit=iunit,file=rfile,form='formatted',status='old')

C Begin reading the file.

      do jvar1=1,nvar1         ! Q2
         do jvar2=1,nvar2      ! W
           read(iunit,FMT=15,err=1000) dummyX,  var2(jvar2), 
     1                                 dummyY,  var1(jvar1)   ! w,q2
 15        format(A8,f4.2,A7,f7.5)
	   if (jvar1.eq.1)     var1tmp1=var1(jvar1)
	   if (jvar1.eq.nvar1) var1tmp2=var1(jvar1)
	   if (jvar2.eq.1)     var2tmp1=var2(jvar2)
	   if (jvar2.eq.nvar2) var2tmp2=var2(jvar2)
C SL+
           read(iunit,'(a)',err=1000) dummy
           read(iunit,*,err=1000)  sf1(jvar1,jvar2), sf2(jvar1,jvar2),  
     &                             sf3(jvar1,jvar2), sf4(jvar1,jvar2),
     &                             sf5(jvar1,jvar2), sf6(jvar1,jvar2) 
           read(iunit,*,err=1000)  sf7(jvar1,jvar2), sf8(jvar1,jvar2),
     &                             sf9(jvar1,jvar2),sf10(jvar1,jvar2),
     &                            sf11(jvar1,jvar2),sf12(jvar1,jvar2)
C SL-
           read(iunit,'(a)',err=1000) dummy
           read(iunit,*,err=1000) dumvar1          ,dumvar2,
     &                            sf13(jvar1,jvar2),sf14(jvar1,jvar2),
     &                            sf15(jvar1,jvar2),sf16(jvar1,jvar2) 
           read(iunit,*,err=1000) sf17(jvar1,jvar2),sf18(jvar1,jvar2),
     &                            sf19(jvar1,jvar2),sf20(jvar1,jvar2),
     &                            sf21(jvar1,jvar2),sf22(jvar1,jvar2)
C EL+
           read(iunit,'(a)',err=1000) dummy
           read(iunit,*,err=1000) sf23(jvar1,jvar2),sf24(jvar1,jvar2),
     &                            sf25(jvar1,jvar2),sf26(jvar1,jvar2),
     &                            sf27(jvar1,jvar2),sf28(jvar1,jvar2) 
           read(iunit,*,err=1000) sf29(jvar1,jvar2),sf30(jvar1,jvar2),
     &                            sf31(jvar1,jvar2),sf32(jvar1,jvar2),
     &                            sf33(jvar1,jvar2),sf34(jvar1,jvar2)
C EL-
           read(iunit,'(a)',err=1000) dummy
           read(iunit,*,err=1000) dumvar1         ,dumvar2,
     &                            dumvar3         ,dumvar4,
     &                            sf35(jvar1,jvar2),sf36(jvar1,jvar2) 
           read(iunit,*,err=1000) sf37(jvar1,jvar2),sf38(jvar1,jvar2),
     &                            sf39(jvar1,jvar2),sf40(jvar1,jvar2),
     &                            sf41(jvar1,jvar2),sf42(jvar1,jvar2)
C ML+
           read(iunit,'(a)',err=1000) dummy
           read(iunit,*,err=1000) dumvar1          ,dumvar2,
     &                            sf43(jvar1,jvar2),sf44(jvar1,jvar2),
     &                            sf45(jvar1,jvar2),sf46(jvar1,jvar2) 
           read(iunit,*,err=1000) sf47(jvar1,jvar2),sf48(jvar1,jvar2),
     &                            sf49(jvar1,jvar2),sf50(jvar1,jvar2),
     &                            sf51(jvar1,jvar2),sf52(jvar1,jvar2)
C ML-
           read(iunit,'(a)',err=1000) dummy
           read(iunit,*,err=1000) dumvar1          ,dumvar2,
     &                            sf53(jvar1,jvar2),sf54(jvar1,jvar2),
     &                            sf55(jvar1,jvar2),sf56(jvar1,jvar2)
           read(iunit,*,err=1000) sf57(jvar1,jvar2),sf58(jvar1,jvar2),
     &                            sf59(jvar1,jvar2),sf60(jvar1,jvar2),
     &                            sf61(jvar1,jvar2),sf62(jvar1,jvar2)
         enddo
      enddo

C Close the input file.

      close(iunit)

C Determine minimum/maximum of variables

      if (var1tmp1.gt.var1tmp2) then
	     var1_max = var1tmp1
	     var1_min = var1tmp2
      else
	     var1_max = var1tmp2
	     var1_min = var1tmp1
      endif

      if (var2tmp1.gt.var2tmp2) then
         var2_max = var2tmp1
         var2_min = var2tmp2
      else
         var2_max = var2tmp2
         var2_min = var2tmp1
      endif

C Now construct natural cubic splines in the 2nd dimension (VAR2) and
C calculate the 2nd derivatives wrt. VAR2 (2nd entry of each D2SF) of
C each structure function (SF). This is done only once, so do it here.
      if(method_spline.eq.1) then
      call splie2(var2,sf1, nvar1,nvar2,d2sf1 )
      call splie2(var2,sf2, nvar1,nvar2,d2sf2 )
      call splie2(var2,sf3, nvar1,nvar2,d2sf3 )
      call splie2(var2,sf4, nvar1,nvar2,d2sf4 )
      call splie2(var2,sf5, nvar1,nvar2,d2sf5 )
      call splie2(var2,sf6, nvar1,nvar2,d2sf6 )
      call splie2(var2,sf7, nvar1,nvar2,d2sf7 )
      call splie2(var2,sf8, nvar1,nvar2,d2sf8 )
      call splie2(var2,sf9, nvar1,nvar2,d2sf9 )
      call splie2(var2,sf10,nvar1,nvar2,d2sf10)
      call splie2(var2,sf11,nvar1,nvar2,d2sf11)
      call splie2(var2,sf12,nvar1,nvar2,d2sf12)
      call splie2(var2,sf13,nvar1,nvar2,d2sf13)
      call splie2(var2,sf14,nvar1,nvar2,d2sf14)
      call splie2(var2,sf15,nvar1,nvar2,d2sf15)
      call splie2(var2,sf16,nvar1,nvar2,d2sf16)
      call splie2(var2,sf17,nvar1,nvar2,d2sf17)
      call splie2(var2,sf18,nvar1,nvar2,d2sf18)
      call splie2(var2,sf19,nvar1,nvar2,d2sf19)
      call splie2(var2,sf20,nvar1,nvar2,d2sf20)
      call splie2(var2,sf21,nvar1,nvar2,d2sf21)
      call splie2(var2,sf22,nvar1,nvar2,d2sf22)
      call splie2(var2,sf23,nvar1,nvar2,d2sf23)
      call splie2(var2,sf24,nvar1,nvar2,d2sf24)
      call splie2(var2,sf25,nvar1,nvar2,d2sf25)
      call splie2(var2,sf26,nvar1,nvar2,d2sf26)
      call splie2(var2,sf27,nvar1,nvar2,d2sf27)
      call splie2(var2,sf28,nvar1,nvar2,d2sf28)
      call splie2(var2,sf29,nvar1,nvar2,d2sf29)
      call splie2(var2,sf30,nvar1,nvar2,d2sf30)
      call splie2(var2,sf31,nvar1,nvar2,d2sf31)
      call splie2(var2,sf32,nvar1,nvar2,d2sf32)
      call splie2(var2,sf33,nvar1,nvar2,d2sf33)
      call splie2(var2,sf34,nvar1,nvar2,d2sf34)
      call splie2(var2,sf35,nvar1,nvar2,d2sf35)
      call splie2(var2,sf36,nvar1,nvar2,d2sf36)
      call splie2(var2,sf37,nvar1,nvar2,d2sf37)
      call splie2(var2,sf38,nvar1,nvar2,d2sf38)
      call splie2(var2,sf39,nvar1,nvar2,d2sf39)
      call splie2(var2,sf40,nvar1,nvar2,d2sf40)
      call splie2(var2,sf41,nvar1,nvar2,d2sf41)
      call splie2(var2,sf42,nvar1,nvar2,d2sf42)
      call splie2(var2,sf43,nvar1,nvar2,d2sf43)
      call splie2(var2,sf44,nvar1,nvar2,d2sf44)
      call splie2(var2,sf45,nvar1,nvar2,d2sf45)
      call splie2(var2,sf46,nvar1,nvar2,d2sf46)
      call splie2(var2,sf47,nvar1,nvar2,d2sf47)
      call splie2(var2,sf48,nvar1,nvar2,d2sf48)
      call splie2(var2,sf49,nvar1,nvar2,d2sf49)
      call splie2(var2,sf50,nvar1,nvar2,d2sf50)
      call splie2(var2,sf51,nvar1,nvar2,d2sf51)
      call splie2(var2,sf52,nvar1,nvar2,d2sf52)
      call splie2(var2,sf53,nvar1,nvar2,d2sf53)
      call splie2(var2,sf54,nvar1,nvar2,d2sf54)
      call splie2(var2,sf55,nvar1,nvar2,d2sf55)
      call splie2(var2,sf56,nvar1,nvar2,d2sf56)
      call splie2(var2,sf57,nvar1,nvar2,d2sf57)
      call splie2(var2,sf58,nvar1,nvar2,d2sf58)
      call splie2(var2,sf59,nvar1,nvar2,d2sf59)
      call splie2(var2,sf60,nvar1,nvar2,d2sf60)
      call splie2(var2,sf61,nvar1,nvar2,d2sf61)
      call splie2(var2,sf62,nvar1,nvar2,d2sf62)
      elseif(method_spline.eq.2) then
C****************************************************************
C the preparation for linear spline fitting
C
      do jvar1=1,nvar1         ! Q2
        var(jvar1) = var1(jvar1)
      enddo
      do jvar2= nvar1+1,nvar1+nvar2         ! W
        var(jvar2) = var2(jvar2-nvar1)
      enddo
C****************************************************************
      endif
      return
     
c 1000 STOP ' Error reading the structure function input file.'
 1000 print *, 'Error reading ', 'Q2= ',(0.2+0.05*(jvar1-1)), 
     1         'W= ', (1.1+0.01*(jvar2-1))
      stop
      end
