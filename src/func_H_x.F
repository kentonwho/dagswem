      double precision function func_H_x(x,y,t)
c
      implicit none
      real*8 x,y,t,x1,x2,omega,xi0,tau,H0,y1,y2
      omega=.0001405
      tau=3456.
      x1=40000.
      x2=151800.4503056
      y1=10000.
      y2=54720.1801222
      xi0=.25
      H0=2.
c
      func_H_x=-2.*omega*xi0*dsin(omega*(x-x1))/dcos(omega*(x2-x1))
     $     *dcos(omega*(y-y1))/dcos(omega*(y2-y1))
     $     *dcos(omega*(t+tau))
      return
      end

