C***********************************************************************
C
C     SUBROUTINE LAND_EDGE_HYDRO( )
C
C     This subroutine does the following:
C
C       1.  Calculates the values of the necessary variables at the edge
C           gauss points for NO-NORMAL FLOW edges
C       2.  Calls the appropriate subroutine to compute the flux at
C           these points.
C       3.  Calls the appropriate subroutine to compute the boundary
C           integrals.
C
C     Written by Ethan Kubatko (06-11-2004)
C
C-----------------------------------------------------------------------
C
C     01-02-2007, sb, Modified for LDG
C     08-xx-2005, sb, Modifications for wetting/drying
C
C***********************************************************************

      SUBROUTINE LAND_EDGE_HYDRO(lvl)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED
      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED, lvl
      REAL(SZ) TX, TY, W_IN, DEN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN
!CND
      real(sz) yedg1,xedg1,func_qx,func_qy,func_H

      DO 1000 L = 1,NLEDS
      
C.....Retrieve the global and local edge number

        GED = NLEDN(L)
        LED = NEDSD(1,GED)
!CND
        n1=nedno(1,ged)
        n2=nedno(2,ged)
        x1=x(n1)
        y1=y(n1)
        x2=x(n2)
        y2=y(n2)
!

C.....Retrieve the elements which share the edge

        EL_IN = NEDEL(1,GED)

C.....If edge element is not assign this ts level, skip CJT
        IF (TsGroup(EL_IN) /= lvl) GOTO 1000
        
C.....If the element is dry then skip the edge calculation

        IF (WDFLG(EL_IN).EQ.0) GOTO 1000

C.....Retrieve the components of the normal vector to the edge
        NX = COSNX(GED)
        NY = SINNX(GED)

           
      
C.....Set the components for the tangential vector to the edge
      
        TX = -NY
        TY =  NX

C.....Compute ZE, QX, QY, and HB at each Gauss point

        DO I = 1,NEGP

          ZE_IN = ZE(1,EL_IN,IRK)
          QX_IN = QX(1,EL_IN,IRK)
          QY_IN = QY(1,EL_IN,IRK)
          HB_IN = BATHED(I,LED,EL_IN)
CND
          SFAC_IN = SFACED(I,LED,EL_IN)
C
          LZ_XX_IN = LZ(1,1,1,EL_IN)
          LZ_XY_IN = LZ(1,1,2,EL_IN)
          LZ_YX_IN = LZ(1,2,1,EL_IN)
          LZ_YY_IN = LZ(1,2,2,EL_IN)

!CND analytic
c        if (nx==1.0.or.nx==-1.0) then
!          yedg1=xegp(i)*(y2-y1)/2.0d0 + (y1+y2)/2.0d0
c           write(*,*) x2,yedg1
c           write(*,*) func_qx(x2,yedg1,timedg)
c        endif
c        if (ny==1.0.or.ny==-1.0) then
!          xedg1=xegp(i)*(x2-x1)/2.0d0 + (x1+x2)/2.0d0
c           write(*,*) xedg1,y1
c           write(*,*) func_qy(xedg1,y1,timedg)
c        endif
!
C.....Compute the solution at the interior state
C.....(modified for wetting and drying)
c           write(*,*) xedg1,yedg1,nx,ny

          DO K = 2,DOF
          
            ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)
            QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)
            QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)

            ! LDG terms
            LZ_XX_IN = LZ_XX_IN + LZ(K,1,1,EL_IN)*PHI_EDGE(K,I,LED)
            LZ_XY_IN = LZ_XY_IN + LZ(K,1,2,EL_IN)*PHI_EDGE(K,I,LED)
            LZ_YX_IN = LZ_YX_IN + LZ(K,2,1,EL_IN)*PHI_EDGE(K,I,LED)
            LZ_YY_IN = LZ_YY_IN + LZ(K,2,2,EL_IN)*PHI_EDGE(K,I,LED)

          ENDDO
! CND analytic
c          ze_in=func_H(xedg1,yedg1,timedg)-2.0
c          qx_in=func_qx(xedg1,yedg1,timedg)
c          qy_in=func_qy(xedg1,yedg1,timedg)
c          if (nx==1.0.or.nx==-1.0) then
c             q_n_ext=func_qx(x1,yedg1,timedg)
c             q_n_int=q_n_ext
c          endif
c          if (ny==1.0.or.ny==-1.0) then
c             q_n_ext=func_qy(xedg1,y1,timedg)
c             q_n_int=q_n_ext
c             ze_in=func_H(xedg1,y1,timedg)-2.0
c          endif

        
C.....Compute the velocity in the normal and tangental direction
        
          Q_N_INT = QX_IN*NX + QY_IN*NY
          Q_T_INT = QX_IN*TX + QY_IN*TY

C.....Reflect the velocity in the normal direction

          Q_N_EXT = -Q_N_INT
          Q_T_EXT =  Q_T_INT

        
C.....Compute the x and y components of the external state flow

          DEN = 1.D0/(NX*TY - NY*TX)
          QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)*DEN
          QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)*DEN

          ZE_EX = ZE_IN
          HB_EX = HB_IN
          SFAC_EX = SFAC_IN
C.....Compute the flux
c          write(*,*) 'qx_in, qx_ex, qy_in, qy_ex ',
c     $         qx_in,qx_ex,qy_in,qy_ex
          CALL NUMERICAL_FLUX()
          
C          IF (EL_IN.EQ.242) THEN
C              PRINT*,F_HAT,G_HAT,H_HAT
C              STOP
C            ENDIF

C.....Add LDG terms
          
          G_HAT = G_HAT + LZ_XX_IN*NX*SFAC_IN + LZ_XY_IN*NY
          H_HAT = H_HAT + LZ_YX_IN*NX*SFAC_IN + LZ_YY_IN*NY

C.....Compute the edge integral
C.....(modified for wetting and drying)

C.....Compute the edge integral
          DO K = 1,DOF

            W_IN = 2.0*M_INV(K)/AREAS(EL_IN)*XLEN(GED)*
     &             PHI_EDGE(K,I,LED)*WEGP(I)

            RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT
            RHS_QX(K,EL_IN,IRK) = RHS_QX(K,EL_IN,IRK) - W_IN*G_HAT
            RHS_QY(K,EL_IN,IRK) = RHS_QY(K,EL_IN,IRK) - W_IN*H_HAT

          ENDDO

        ENDDO

1000  CONTINUE
      
      RETURN
      END SUBROUTINE


C***********************************************************************
C
C     SUBROUTINE LAND_EDGE_HYDRO_POST( )
C
C     This subroutine does the following:
C
C       1.  Direct velocity at each node on land edges toward
C           the tangential direction of the edge.
C           (This procedure needs to be called after the wetting and 
C            drying post-process.)
C
C     Written by Shintaro Bunya (02-13-2006)
C
C***********************************************************************

      SUBROUTINE LAND_EDGE_HYDRO_POST()

C.....Use appropriate modules

      USE SIZES,ONLY : SZ,MYPROC
      USE GLOBAL,ONLY : X,Y
      USE DG,ONLY : NLEDS,NLEDN,NEDSD,NEDEL,COSNX,SINNX,DOF,QX,QY,
     &              PHI_CORNER, IRK, EL_UPDATED, NEDNO, XEGP, YEGP

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED, EL, K, KK, NOD1, NOD2, NOD3, NEDGES
      REAL(SZ) NX, NY
      REAL(SZ) TX, TY
      REAL(SZ) DIR
      REAL(SZ) QXP(3),QXPN(3)
      REAL(SZ) QYP(3),QYPN(3)
      REAL(SZ) QP(3)
!CND
      integer n1,n2
      real(sz) yed1,yed2,x1,x2,y1,y2
!
      write(*,*) 'in land edge hydro'
      
      DO 1000 L=1, NLEDS

C.....Retrieve the global and local edge number
      GED = NLEDN(L)
      LED = NEDSD(1,GED)
! CND
      n1=nedno(1,ged)
      n2=nedno(2,ged)
      x1=x(n1)
      x2=x(n2)
      y1=y(n1)
      y2=y(n2)
!

C.....Retrieve the elements which share the edge

      EL = NEDEL(1,GED)
C.....Apply the following sequences only on dry, drying, or wetting elements
      IF(EL_UPDATED(EL).EQ.0) CYCLE

C.....Retrieve the components of the normal vector to the edge
        
      NX = COSNX(GED)
      NY = SINNX(GED)
! CND
      write(*,*) 'nx ,ny ',nx,ny
      if (nx==1.0d0) then
         yed1=xegp(1)/2.*(y2-y1)+.5d0*(y1+y2)
         yed2=xegp(2)/2.*(y2-y1)+.5d0*(y1+y2)
         write(*,*) x1,x2,yed1,yed2
      else
         yed1=0.0
         yed2=0.0
      endif
      
C.....Set the components for the tangential vector to the edge
      
      TX = -NY


C.....Compute nodal QX and QY
      QXP(:) = 0.D0
      QYP(:) = 0.D0

      DO K = 1,3
        DO KK = 1,DOF
          QXP(K) = QXP(K) + PHI_CORNER(KK,K)*QX(KK,EL,IRK+1)
          QYP(K) = QYP(K) + PHI_CORNER(KK,K)*QY(KK,EL,IRK+1)
        ENDDO
        QP(K) = SQRT(QXP(K)*QXP(K)+QYP(K)*QYP(K))
        
      ENDDO

C.....Retrieve the end nodes of the edge and the other node
      NOD1 = MOD(LED+0,3)+1
      NOD2 = MOD(LED+1,3)+1
      NOD3 = MOD(LED+2,3)+1
      
C.....Modify QX and QY at the end nodes

C.......Node 1
      IF(ABS( QXP(NOD1)*TX + QYP(NOD1)*TY ).GT.(10.D0**(-20.D0))) THEN
        DIR = ( QXP(NOD1)*TX + QYP(NOD1)*TY )
     &       / ABS( QXP(NOD1)*TX + QYP(NOD1)*TY )

#if 1
        QXPN(NOD1) = DIR*QP(NOD1)*TX
        QYPN(NOD1) = DIR*QP(NOD1)*TY
#else
        QXPN(NOD1) = (QXP(NOD1)*TX + QYP(NOD1)*TY)*TX
        QYPN(NOD1) = (QXP(NOD1)*TX + QYP(NOD1)*TY)*TY
#endif
      ELSE
        QXPN(NOD1) = 0.D0
        QYPN(NOD1) = 0.D0
      ENDIF

C.......Node 2
      IF(ABS( QXP(NOD2)*TX + QYP(NOD2)*TY ).GT.(10.D0**(-20.D0))) THEN
        DIR = ( QXP(NOD2)*TX + QYP(NOD2)*TY )
     &       / ABS( QXP(NOD2)*TX + QYP(NOD2)*TY )

#if 1
        QXPN(NOD2) = DIR*QP(NOD2)*TX
        QYPN(NOD2) = DIR*QP(NOD2)*TY
#else
        QXPN(NOD2) = (QXP(NOD2)*TX + QYP(NOD2)*TY)*TX
        QYPN(NOD2) = (QXP(NOD2)*TX + QYP(NOD2)*TY)*TY
#endif
      ELSE
        QXPN(NOD2) = 0.D0
        QYPN(NOD2) = 0.D0
      ENDIF

C.......Node 3
      QXPN(NOD3) = QXP(NOD3)
      QYPN(NOD3) = QYP(NOD3)

C.....Inverse the modified QX and QY
      QX(1,EL,IRK+1) = 1.D0/3.D0*(QXPN(1)+QXPN(2)+QXPN(3))
      QX(2,EL,IRK+1) =-1.D0/6.D0*(QXPN(1)+QXPN(2))+1.D0/3.D0*QXPN(3)
      QX(3,EL,IRK+1) =-0.5D0*QXPN(1)+0.5D0*QXPN(2)

      QY(1,EL,IRK+1) = 1.D0/3.D0*(QYPN(1)+QYPN(2)+QYPN(3))
      QY(2,EL,IRK+1) =-1.D0/6.D0*(QYPN(1)+QYPN(2))+1.D0/3.D0*QYPN(3)
      QY(3,EL,IRK+1) =-0.5D0*QYPN(1)+0.5D0*QYPN(2)

C      print *,'T=',TX,TY
C      print *,'QXP =',QXP (1),QXP (2),QXP (3)
C      print *,'QYP =',QYP (1),QYP (2),QYP (3)
C      print *,'QXPN=',QXPN(1),QXPN(2),QXPN(3)
C      print *,'QYPN=',QYPN(1),QYPN(2),QYPN(3)
C      print *,'EL=', EL,',N1=',NOD1,',N2=',NOD2,',N3=',NOD3

1000  CONTINUE
      RETURN
      END SUBROUTINE
