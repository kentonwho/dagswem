      double precision function func_sy(x,y,t)
c
      implicit none
      real*8 x,y,t,x1,x2,omega,xi0,tau,H0,func_qxdH,func_H,func_qx,
     $     func_qy,func_qy_t,
     $     func_qx_t,func_qx_x,func_H_x,g,taubf,func_qydH,func_qx_y,
     $     func_qy_x,func_qy_y,func_H_y,y1,y2
      omega=.0001405
      tau=3456
      x1=40000
      x2=151800.4503056866
      y1=10000
      y2=54720.18012227464
      xi0=.25
      H0=2
      g=9.81
c
      func_qxdH=func_qx(x,y,t)/func_H(x,y,t)
      func_qydH=func_qy(x,y,t)/func_H(x,y,t)
      func_sy=func_qy_t(x,y,t)+2*func_qydH*func_qy_y(x,y,t)
c     $     +taubf*func_qy(x,y,t)
     $     -func_qydH**2*func_H_y(x,y,t)+g*func_H(x,y,t)*func_H_y(x,y,t)
     $     +func_qy_x(x,y,t)*func_qxdH+func_qx_x(x,y,t)*func_qydH
     $     -func_qydH*func_qxdH*func_H_x(x,y,t)
c      func_sy =0.
      return
      end
