      subroutine daresbury(q2,win,ep,cscm,phcm,opt,sigma0)
      
c     cscm=cos(theta_pi_cm)
c     phcm=phi_pi_cm
c     win=invariant hadronic mass, .ge. 1.205
c     q2=Q squared
c     ep=polarization, epsilon
c     opt=1, calculate cross section as function of angles

      real a(4,11),c(2,11),d(3,11),wdat(11),pi,sigma0
      integer opt
      
      data a/16.17,-0.57,-7.50,1.87,19.3,1.68,-9.75,0.62,
     1        15.3,2.72,-9.38,-1.64,9.65,0.35,-7.24,0.19,
     1        6.65,0.93,-4.71,-0.57,5.34,1.16,-4.14,-0.98,
     1        3.90,1.16,-2.77,-0.95,
     1        3.83,0.51,-3.29,-0.29,3.31,0.21,-3.16,-0.34,
     1        3.83,0.52,-3.28,-0.68,4.30,0.83,-5.24,-2.73/
     
      data c/-8.00,-5.08,-10.37,0.16,-6.37,0.74,-5.32,-0.58,
     1       -2.39,-0.23,-1.44,0.34,-1.18,0.6,-1.03,-0.22,
     1       -1.23,-0.57,-1.28,-1.88,-2.44,-2.23/
     
      data d/-1.14,-2.72,2.87,-0.36,-4.43,-0.48,-0.03,-2.09,0.26,
     1        0.87,-1.06,-1.21,0.02,-0.37,-0.21,-0.18,-0.2,-0.7,
     1        0.01,-0.51,-1.18,-0.43,-0.17,0.07,-0.25,0.18,-0.26,
     1       -0.57,0.2,-0.37,-0.38,-0.3,-1.22/
     
      data wdat/1.205,1.235,1.265,1.295,1.325,1.355,1.385,1.415,
     11.445,1.475,1.505/
     
      data pi/3.14159/,d2r/0.0174532/
      
      w = max(win,1.205)
     
      i = max(1,min(10,int((w-1.205)/0.03)+1))
      
      aa1 = (w-wdat(i))*(a(1,i+1)-a(1,i))/0.035 + a(1,i)
      aa2 = (w-wdat(i))*(a(2,i+1)-a(2,i))/0.035 + a(2,i)
      aa3 = (w-wdat(i))*(a(3,i+1)-a(3,i))/0.035 + a(3,i)
      aa4 = (w-wdat(i))*(a(4,i+1)-a(4,i))/0.035 + a(4,i)
      cc1 = (w-wdat(i))*(c(1,i+1)-c(1,i))/0.035 + c(1,i)
      cc2 = (w-wdat(i))*(c(2,i+1)-c(2,i))/0.035 + c(2,i)
      dd1 = (w-wdat(i))*(d(1,i+1)-d(1,i))/0.035 + d(1,i)
      dd2 = (w-wdat(i))*(d(2,i+1)-d(2,i))/0.035 + d(2,i)
      dd3 = (w-wdat(i))*(d(3,i+1)-d(3,i))/0.035 + d(3,i)
      
      if (opt.eq.1) then      
        snsm  = sqrt(1-cscm**2)
        c2phi = cos(2*phcm*d2r)
        cphi  = cos(phcm*d2r)
        aa    = aa1+aa2*cscm+aa3*cscm*cscm+aa4*cscm*cscm*cscm
        cc    = cc1+cc2*cscm
        dd    = dd1+dd2*cscm+dd3*cscm*cscm
        dsdo  = aa+ep*cc*snsm**2*c2phi+sqrt(2*ep*(ep+1))*dd*snsm*cphi
      elseif (opt.eq.2) then
        dsdo = 4*pi*(aa1+0.333*aa3)
      endif
      
      sigma0 = dsdo*(q2/0.57)*((1+(0.57/0.71))/(1+(q2/0.71)))**4
     
        
      end
      














