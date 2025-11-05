C***********************************************************************
C
C     SUBROUTINE RHS_DG_HYDRO()
C
C     This subroutine computes the area integrals for the DG hydro and
C     adds them into the RHS.
C
C     Written by Ethan Kubatko (06-11-2004)
C
C-----------------------------------------------------------------------
C
C     Feb 23, 2007, sb, Modified for better performance
C     Jan 02, 2007, sb, Modified for LDG
C     Aug xx, 2005, sb, Modified for wetting/drying
C
C***********************************************************************
      SUBROUTINE RHS_DG_HYDRO(LVL)
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED
      USE NodalAttributes, ONLY : TAU, IFLINBF, IFHYBF, HBREAK, FTHETA,
     &                            FGAMMA,LoadManningsN,ManningsN,CF

      USE sizes, ONLY: myproc



      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER L, LVL, elem  ! CJT
      REAL(SZ) DPSIDX(3), DPSIDY(3)
      REAL(SZ) RainInt, XR_L, XR_R, YR_L, YR_R
      REAL(SZ) AREA, IMASS, TKX, TKY
      REAL(SZ) PHI_AREA_KI,MN_IN
      REAL(SZ) LZ_XX, LZ_XY, LZ_YX, LZ_YY
      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL, FW_NL
      REAL(SZ) HUU, HVV, HUV, GH2, DPHIDX, DPHIDY
      REAL(SZ) DEPTH_C, FH_NL_C, UX_C, UY_C, UMAG_C, DTDPH
c
      real(sz) SFACQUAD,rkeep,rcheck1
c
! CND  analytic test case
      real(sz)  func_H,func_qx,func_qy,func_qx_x,func_qx_y,func_qy_x,
     $     func_qy_y,func_qx_t,func_qy_t,func_H_x,func_H_y,
     $     func_sx,func_sy,xintp,yintp


      DTDPH = 1.D0/DTDom(LVL)
      DO 1000 L = 1,Nelem(lvl)

         elem = ElemNum(lvl,L)
c
cnd
         advectqx(elem)=0.0
         advectqx(elem)=0.0
         sourceqx(elem)=0.0
         sourceqy(elem)=0.0
cnd
c         if (dabs(rhs_qx(2,l,irk)+rtemp2(l)).ge.1.d-6) then
c            write(779,*) l,irk,rhs_qx(2,l,irk),rtemp2(l),
c     $           ze(1,l,irk),ze(2,l,irk),ze(3,l,irk),
c     $           qx(1,l,irk),qx(2,l,irk),qx(3,l,irk)
c         endif
cnd
      
C.......If element is dry then skip calculations
      
        IF (WDFLG(elem).EQ.0) then
           GOTO 1000
        endif

C.......Retrieve the global node numbers for the element

        N1 = NM(elem,1)
        N2 = NM(elem,2)
        N3 = NM(elem,3)

C.......Compute avaraged values
C.......These will be used later when bottom friction is computed

        DEPTH_C = HB(1,elem,1) + ZE(1,elem,IRK)
        FH_NL_C = 1.D0/(NLEQ*DEPTH_C + LEQ)
        UX_C = QX(1,elem,IRK)*FH_NL_C
        UY_C = QY(1,elem,IRK)*FH_NL_C
        UMAG_C = SQRT(UX_C*UX_C + UY_C*UY_C)
        
C.......Compute derivatives of Lagrange basis functions at nodes
        
        IF ((NWS.NE.0).OR.(NTIP.NE.0)) THEN
          DPSIDX(1) = DRPSI(1)*DRDX(elem) + DSPSI(1)*DSDX(elem)
          DPSIDX(2) = DRPSI(2)*DRDX(elem) + DSPSI(2)*DSDX(elem)
          DPSIDX(3) = DRPSI(3)*DRDX(elem) + DSPSI(3)*DSDX(elem)
          DPSIDY(1) = DRPSI(1)*DRDY(elem) + DSPSI(1)*DSDY(elem)
          DPSIDY(2) = DRPSI(2)*DRDY(elem) + DSPSI(2)*DSDY(elem)
          DPSIDY(3) = DRPSI(3)*DRDY(elem) + DSPSI(3)*DSDY(elem)
        ENDIF

C.......Compute ZE, QX, QY, and HB at each area Gauss quadrature point

        DO I = 1,NAGP
      
          ZE_IN = ZE(1,elem,IRK)
          QX_IN = QX(1,elem,IRK)
          QY_IN = QY(1,elem,IRK)
        
          HB_IN = BATH(I,elem)
          DHB_X = DBATHDX(I,elem)
          DHB_Y = DBATHDY(I,elem)

          LZ_XX = LZ(1,1,1,elem)
          LZ_XY = LZ(1,1,2,elem)
          LZ_YX = LZ(1,2,1,elem)
          LZ_YY = LZ(1,2,2,elem)
CND
          SFACQUAD = SFAC_ELEM(I,elem)

          DO K = 2,DOF

            ZE_IN = ZE_IN + ZE(K,elem,IRK)*PHI_AREA(K,I)
            QX_IN = QX_IN + QX(K,elem,IRK)*PHI_AREA(K,I)
            QY_IN = QY_IN + QY(K,elem,IRK)*PHI_AREA(K,I)
c
            LZ_XX = LZ_XX + LZ(K,1,1,elem)*PHI_AREA(K,I)
            LZ_XY = LZ_XY + LZ(K,1,2,elem)*PHI_AREA(K,I)
            LZ_YX = LZ_YX + LZ(K,2,1,elem)*PHI_AREA(K,I)
            LZ_YY = LZ_YY + LZ(K,2,2,elem)*PHI_AREA(K,I)

          ENDDO

          DEPTH = ZE_IN + HB_IN

C.........Compute continuity fluxes

          F1_NL = NLEQ + LEQ*HB_IN
          FX_IN = QX_IN*F1_NL*SFACQUAD
c          FX_IN = QX_IN*F1_NL
          FY_IN = QY_IN*F1_NL

C.........Compute momentum flux terms

          FU_NL = NLEQ*QX_IN
          FV_NL = NLEQ*QY_IN
          FG_NL = NLEQG*ZE_IN*WDFLG(elem)
          FH_NL = 1.D0/(NLEQ*DEPTH + LEQ)
          U_IN  = QX_IN*FH_NL
          V_IN  = QY_IN*FH_NL

          HUU = FU_NL*U_IN
          HVV = FV_NL*V_IN
          HUV = FU_NL*V_IN
          GH2 = FG_NL*(0.5D0*ZE_IN + HB_IN) + FG_L*ZE_IN

C.........Compute x momentum fluxes

          GX_IN = (HUU + GH2 + LZ_XX)*SFACQUAD
          GY_IN = HUV + LZ_XY
cnd
          advectqx(elem)=advectqx(elem)+gx_in+gy_in
cnd

C.........Compute y momentum fluxes

          HX_IN = (HUV + LZ_YX)*SFACQUAD
          HY_IN = HVV + GH2 + LZ_YY

cnd
          advectqy(elem)=advectqy(elem)+hx_in+hy_in
cnd

C.........Compute the friction factor
          if (LoadManningsN) then
c             MN_IN=MANN(1,L)
c             do k=2,dof
c                MN_IN=MN_IN + MANN(K,L)*PHI_AREA(K,I)
c             enddo
             fric_el(elem)=G*
     $            ((ManningsN(n1)+ManningsN(n2)+ManningsN(n3))/3.)**2
c     $            MN_IN**2
     $            /(DEPTH**(1.d0/3.d0))
             if (fric_el(elem).lt.CF) fric_el(elem)=CF
          endif
          TAU = FRIC_EL(elem)
C          IF (IFLINBF.EQ.0) THEN
C            UMAG = SQRT( U_IN*U_IN + V_IN*V_IN )
C            TAU  = TAU*UMAG*FH_NL
C            IF (IFHYBF.EQ.1) TAU = TAU*
C     &             (1.D0  + (HBREAK*FH_NL)**FTHETA)**(FGAMMA/FTHETA)
C          ENDIF
C Modified to compute TAU using elemental averages.
C This seems necessary to avoid exessive bottom friction
C at wetting-drying fronts where the total column height is very
C small. S.B. 9-Feb-2008
          IF (IFLINBF.EQ.0) THEN
             UMAG = SQRT( U_IN*U_IN + V_IN*V_IN )
c cnd modified 4/23/10 to test friction 
c             TAU  = TAU*UMAG_C*FH_NL_C
             TAU  = TAU*UMAG*FH_NL
             IF (IFHYBF.EQ.1) TAU = TAU*
     &            (1.D0  + (HBREAK*FH_NL_C)**FTHETA)**(FGAMMA/FTHETA)
c     &            (1.D0  + (HBREAK*FH_NL)**FTHETA)**(FGAMMA/FTHETA)
C It is numerically probable that the bottom friction is large enoght
C to reverse the direction of currents backward due to a too small column
C height even though it does not happen in reality. To avoid this, the MIN
C function bellow is added. It is expected that this MIN function upper-limits
C TAU so the bottom friction force does not reverse the currents within
C half a time step.  S.B. 9-Feb-2008
c             IF(TAU.GT.2.D0*DTDPH) PRINT *, "TAU = ", TAU, 2.D0*DTDPH,
c     $           u_in,v_in,myproc,L
c             TAU = MIN(TAU, 2.D0*DTDPH)
             TAU = MIN(TAU, .9D0*DTDPH)
          ENDIF
c          IF (RAMPDG.LT.1.D0) TAU = MAX(TAU,0.001)

C.........Computer sources and sinks for mass conservation
          SOURCE_Z = 0.D0 ! cjt

C.........Compute the x momentum source/sink terms

          SOURCE_X =

C.........1.) Friction term

     &    - TAU*QX_IN

C.........2.) Bathymetric slope term

     &    + FG_NL*DHB_X*SFACQUAD

C.........3.) Coriolis force

     &    + CORI_EL(elem)*QY_IN

C.....Compute the y momentum source/sink terms

          SOURCE_Y =

C.........1.) Friction term

     &    - TAU*QY_IN

C.........2) Bathymetric slope term

     &    + FG_NL*DHB_Y

C.........3.) Coriolis force

     &    - CORI_EL(elem)*QX_IN

     
C.........4.) Wind and pressure forcing (in x and y)

C!!!!!!!!!!!!!!!! HARWIRE RAIN FOR NOW (CJT)
            ! Hardwire rainfall CJT
             RainInt = 0.d0
!             RainInt = 0.001d0
!             RI2(N1) = RainInt; RI2(N2) = RainInt; RI2(N3) = RainInt;
!
!!             XR_L = 400.d0; YR_L = 400.d0;
!!             XR_R = 600.d0; YR_R = 600.d0;
!!             RI2(N1) = 0.d0; RI2(N2) = 0.d0; RI2(N3) = 0.d0;
!!             IF(X(N1)>=XR_L.and.X(N1)<=XR_R) THEN
!!               IF (Y(N1)>=YR_L.and.Y(N1)<=YR_R) THEN
!!                 RI2(N1) = RainInt; 
!!               ENDIF
!!             ENDIF
!!             IF(X(N2)>=XR_L.and.X(N2)<=XR_R) THEN
!!               IF (Y(N2)>=YR_L.and.Y(N2)<=YR_R) THEN
!!                 RI2(N2) = RainInt; 
!!               ENDIF
!!             ENDIF
!!             IF(X(N3)>=XR_L.and.X(N3)<=XR_R) THEN
!!               IF (Y(N3)>=YR_L.and.Y(N3)<=YR_R) THEN
!!                 RI2(N3) = RainInt; 
!!               ENDIF
!!             ENDIF
!
!            ! Water source from rainfall
!             SOURCE_Z = SOURCE_Z + RI2(N1)*PSI1(I) + RI2(N2)*PSI2(I)
!     &                          + RI2(N3)*PSI3(I)
!
!c            ! Rain forcing (cjt)  [Rain drop shearing of water surface]
!c            SOURCE_X = SOURCE_X - (RI2(N1)*PSI1(I) + RI2(N2)*PSI2(I)
!c     &                          +  RI2(N3)*PSI3(I)) * QX_IN/DEPTH  ! cjt
!c            SOURCE_Y = SOURCE_Y - (RI2(N1)*PSI1(I) + RI2(N2)*PSI2(I)
!c     &                          +  RI2(N3)*PSI3(I)) * QY_IN/DEPTH  ! cjt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  CND analytic test case
!             xintp=x(n1)*psi1(i)+x(n2)*psi2(i)+x(n3)*psi3(i)
!             yintp=y(n1)*psi1(i)+y(n2)*psi2(i)+y(n3)*psi3(i)
!             source_x=source_x+func_sx(xintp,yintp,timedg)
c             source_x=source_x+func_sx(x(n1),y(n1),timedg)*psi1(i)
c     $            +func_sx(x(n2),y(n2),timedg)*psi2(i)
c     $            +func_sx(x(n3),y(n3),timedg)*psi3(i)
!             source_y=source_y+func_sy(xintp,yintp,timedg)
c             source_y=source_y+func_sy(x(n1),y(n1),timedg)*psi1(i)
c     $            +func_sy(x(n2),y(n2),timedg)*psi2(i)
c     $            +func_sy(x(n3),y(n3),timedg)*psi3(i)
c             write(79,*) x(n1),y(n1),func_sx(x(n1),y(n1),timedg),
c     $            func_sy(x(n1),y(n1),timedg),psi1(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          IF (NWS.NE.0) THEN
            FW_NL = 1.D0/F1_NL

            ! Water source from rainfall
            SOURCE_Z = SOURCE_Z + RI2(N1)*PSI1(I) + RI2(N2)*PSI2(I)
     &                          + RI2(N3)*PSI3(I)

            SOURCE_X = SOURCE_X + FW_NL*( WSX2(N1)*PSI1(I)
     &           + WSX2(N2)*PSI2(I)  + WSX2(N3)*PSI3(I) )
c     &                          - G*DEPTH*( PR2(N1)*DPSIDX(1)
cnd
     &           - G*SFACQUAD*DEPTH
     $           *( PR2(N1)*DPSIDX(1)
     &           + PR2(N2)*DPSIDX(2) + PR2(N3)*DPSIDX(3))
            SOURCE_Y = SOURCE_Y + FW_NL*( WSY2(N1)*PSI1(I)
     &           + WSY2(N2)*PSI2(I)  + WSY2(N3)*PSI3(I) )
     &           - G*DEPTH*( PR2(N1)*DPSIDY(1)
     &           + PR2(N2)*DPSIDY(2) + PR2(N3)*DPSIDY(3))

c            ! Rain forcing (cjt)  [Rain drop shearing of water surface]
c            SOURCE_X = SOURCE_X - (RI2(N1)*PSI1(I) + RI2(N2)*PSI2(I)
c     &                          +  RI2(N3)*PSI3(I)) * QX_IN/DEPTH  ! cjt
c            SOURCE_Y = SOURCE_Y - (RI2(N1)*PSI1(I) + RI2(N2)*PSI2(I)
c     &                          +  RI2(N3)*PSI3(I)) * QY_IN/DEPTH  ! cjt



          ENDIF


C.........5) Tidal potential forcing (in x and y)

          IF (NTIP.NE.0) THEN
            SOURCE_X = SOURCE_X + 
c     $            RAMPDG*G*DEPTH*( DPSIDX(1)*TIP2(N1)
cnd
     $            RAMPDG*G*DEPTH*SFACQUAD*( DPSIDX(1)*TIP2(N1)
     &            + DPSIDX(2)*TIP2(N2) + DPSIDX(3)*TIP2(N3) )
            SOURCE_Y = SOURCE_Y + RAMPDG*G*DEPTH*( DPSIDY(1)*TIP2(N1)
     &                      + DPSIDY(2)*TIP2(N2) + DPSIDY(3)*TIP2(N3) )
          ENDIF
cnd
          sourceqx(elem)=sourceqx(elem)+source_x
          sourceqy(elem)=sourceqy(elem)+source_y
C.........Build the rhs

          RHS_ZE(1,elem,IRK) = RHS_ZE(1,elem,IRK) + SRFAC(1,I,elem)*SOURCE_Z ! cjt
          RHS_QX(1,elem,IRK) = RHS_QX(1,elem,IRK) + SRFAC(1,I,elem)*SOURCE_X
          RHS_QY(1,elem,IRK) = RHS_QY(1,elem,IRK) + SRFAC(1,I,elem)*SOURCE_Y

          DO K = 2,DOF
            RHS_ZE(K,elem,IRK) = RHS_ZE(K,elem,IRK) + XFAC(K,I,elem)*FX_IN
     &                    + YFAC(K,I,elem)*FY_IN + SRFAC(K,I,elem)*SOURCE_Z
            RHS_QX(K,elem,IRK) = RHS_QX(K,elem,IRK) + XFAC(K,I,elem)*GX_IN
     &                    + YFAC(K,I,elem)*GY_IN + SRFAC(K,I,elem)*SOURCE_X
            RHS_QY(K,elem,IRK) = RHS_QY(K,elem,IRK) + XFAC(K,I,elem)*HX_IN
     &                    + YFAC(K,I,elem)*HY_IN + SRFAC(K,I,elem)*SOURCE_Y
          ENDDO


      ENDDO     !  I  loop

1000  CONTINUE      
      RETURN
      END SUBROUTINE
