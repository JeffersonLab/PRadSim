C------------------------------------------------------------------------------
C multipole_amps
C
C Fills multipole amplitudes using results from interpolation
C------------------------------------------------------------------------------
        subroutine multipole_amps

        implicit none

        include 'mpintp.inc'
        include 'spp.inc'

        real vec_f(max_mp)
        real factor
        integer i
        if (w.le.1.1) w = 1.1 
        if (q2.le.5.and.w.le.2) then
        call interp(Q2,W,vec_f)
        endif
        if (q2.gt.5.and.w.le.2) then
        call interp(5.0,W,vec_f)
        endif
        if (q2.le.5.and.w.gt.2) then
        call interp(Q2,2.0,vec_f)
        endif
        if (q2.gt.5.and.w.gt.2) then
        call interp(5.0,2.0,vec_f)
        endif
c        call interp(Q2,W,vec_f)
        factor = 0.141383e0     ! Convert from 10**-3/m_pi+ to sqrt(microbarns)
c
c	nu_cm * S = qv_mag_cm * L
c
c We have to multiply nu_cm / qv_mag_cm because we are following Tiator's notation 
c which uses L-multipole rather than S-multipole.
c Note that Sato-Lee uses S-multipole.
c
	    factor = 0.141383d0 * nu_cm / qv_mag_cm
c        factor = nu_cm / qv_mag_cm
c        factor = 0.141383d0 
        sp(0) = cmplx(vec_F(1), vec_F(2) )*factor
        sp(1) = cmplx(vec_F(3), vec_F(4) )*factor
        sp(2) = cmplx(vec_F(5), vec_F(6) )*factor
        sp(3) = cmplx(vec_F(7), vec_F(8) )*factor
        sp(4) = cmplx(vec_F(9), vec_F(10))*factor
        sp(5) = cmplx(vec_F(11),vec_F(12))*factor

	sm(0) = 0.0
        sm(1) = cmplx(vec_F(13),vec_F(14))*factor
        sm(2) = cmplx(vec_F(15),vec_F(16))*factor
        sm(3) = cmplx(vec_F(17),vec_F(18))*factor
        sm(4) = cmplx(vec_F(19),vec_F(20))*factor
        sm(5) = cmplx(vec_F(21),vec_F(22))*factor

	factor = 0.141383d0
c
c        factor = 1.0
        ep(0) = cmplx(vec_F(23),vec_F(24))*factor
        ep(1) = cmplx(vec_F(25),vec_F(26))*factor
        ep(2) = cmplx(vec_F(27),vec_F(28))*factor
        ep(3) = cmplx(vec_F(29),vec_F(30))*factor
        ep(4) = cmplx(vec_F(31),vec_F(32))*factor
        ep(5) = cmplx(vec_F(33),vec_F(34))*factor

	    em(0) = 0.0
		em(1) = 0.0
        em(2) = cmplx(vec_F(35),vec_F(36))*factor
        em(3) = cmplx(vec_F(37),vec_F(38))*factor
        em(4) = cmplx(vec_F(39),vec_F(40))*factor
        em(5) = cmplx(vec_F(41),vec_F(42))*factor

		mp(0) = 0.0
        mp(1) = cmplx(vec_F(43),vec_F(44))*factor
        mp(2) = cmplx(vec_F(45),vec_F(46))*factor
        mp(3) = cmplx(vec_F(47),vec_F(48))*factor
        mp(4) = cmplx(vec_F(49),vec_F(50))*factor
        mp(5) = cmplx(vec_F(51),vec_F(52))*factor

		mm(0) = 0.0
        mm(1) = cmplx(vec_F(53),vec_F(54))*factor
        mm(2) = cmplx(vec_F(55),vec_F(56))*factor
        mm(3) = cmplx(vec_F(57),vec_F(58))*factor
        mm(4) = cmplx(vec_F(59),vec_F(60))*factor
        mm(5) = cmplx(vec_F(61),vec_F(62))*factor
c
c        print *, 'multipoles ------'
c        print *, (sp(i),sm(i),ep(i),em(i),mp(i),mm(i),i=0,5)
        return
        end
