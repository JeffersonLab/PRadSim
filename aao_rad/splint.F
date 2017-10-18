C------------------------------------------------------------------------

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
C
C Given the arrays ya(1:n) and xa(1:n) of length n, which tabulate a function
C and its independent variable (with the xa(i) in order, as in SPLINE),
C respectively, and given the array y2a(1:n), which is the output from SPLINE,
C and given a value of x, this routine returns a cubic-spline interpolated
C value y.  From "Numerical Recipes", Ch. 3.3, p. 110.
      IMPLICIT NONE
C Passed variables:
      INTEGER n
      REAL xa(n),ya(n),y2a(n),x,y
C Local variables:
      INTEGER k,klo,khi
      REAL a,b,h
C
      klo=1
      khi=n
C
C Locate the interpolation point by means of bisection.
C This is optimal if sequential calls to this routine are at random values
C of x. If sequential calls are in order, and closely spaced, it would be
C better to store previous values of klo and khi and test if they remain
C appropriate on the next call.
 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if (xa(k).gt.x) then
            khi=k
         else
            klo=k
         endif
         go to 1
      endif
C klo,khi now bracket x.
      h=xa(khi)-xa(klo)
C The xa's must be distinct:
      if (h.eq.0.e0) stop ' Bad xa input in splint'
C Evaluate now the cubic spline polynomial:
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     &  ((a*a*a-a)*y2a(klo)+(b*b*b-b)*y2a(khi))*(h*h)/6.e0
C
      return
      END
