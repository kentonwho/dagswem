C***********************************************************************
C
C     SUBROUTINE NUMERICAL_FLUX
C
C     This subroutine calls an appropriate numrical flux subroutine
C     based on FLUXTYPE variable.
C
C     Written by Shintaro Bunya (11-08-2006)
C
C***********************************************************************

      SUBROUTINE NUMERICAL_FLUX()
      
      USE DG,ONLY : FLUXTYPE

      IF(FLUXTYPE.EQ.1) THEN
        CALL ROE_FLUX()
      ELSEIF(FLUXTYPE.EQ.2) THEN
        CALL LLF_FLUX()
      ELSEIF (FLUXTYPE.EQ.3) THEN
        CALL HLLC_FLUX()
      ELSE
        STOP 'INVALID FLUXTYPE'
      ENDIF

      RETURN
      END SUBROUTINE


***********************************************************************
C
C     SUBROUTINE ROE_FLUX
C
C     This subroutine computes the Roe flux for the shallow water
C     equations.
C
C     Written by Ethan Kubatko (06-11-2004)
C
C***********************************************************************
      SUBROUTINE ROE_FLUX()
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE
      
C.....Declare local variables
  
      INTEGER II, GED
      REAL(SZ) DEN
      REAL T1, T2, T3, T4
      
C.....Compute the total height of the water column

      HT_IN = ZE_IN*IFNLFA + HB_IN
      HT_EX = ZE_EX*IFNLFA + HB_EX

C.....Compute the velocities in the x and y directions

      U_IN = QX_IN/HT_IN
      U_EX = QX_EX/HT_EX

      V_IN = QY_IN/HT_IN
      V_EX = QY_EX/HT_EX

C.....Compute the Roe averaged variables

      T1    = SQRT(HT_IN)
      T2    = SQRT(HT_EX)
      T3    = SQRT(HT_IN + HT_EX)
      T4    = 1.0D0/(T1 + T2)
      C_ROE = G2ROOT*T3
      U_ROE = T4*(U_IN*T1 + U_EX*T2)
      V_ROE = T4*(V_IN*T1 + V_EX*T2)

C.....Evaluate the eigenvalues at the Roe averaged variables

      EIGVAL(2) = (U_ROE*NX + V_ROE*NY)*IFNLCT
      EIGVAL(1) = EIGVAL(2) + C_ROE
      EIGVAL(3) = EIGVAL(2) - C_ROE
      
C.....Evaluate right eigenvectors at Roe averaged variables

      RI(1,1) = 1.D0
      RI(2,1) = U_ROE*IFNLCT + (IFNLCT*C_ROE + (1-IFNLCT)*G/C_ROE)*NX
      RI(3,1) = V_ROE*IFNLCT + (IFNLCT*C_ROE + (1-IFNLCT)*G/C_ROE)*NY

      RI(1,2) = 0.D0
      RI(2,2) = -NY
      RI(3,2) =  NX

      RI(1,3) = 1.D0
      RI(2,3) = U_ROE*IFNLCT - (IFNLCT*C_ROE + (1-IFNLCT)*G/C_ROE)*NX
      RI(3,3) = V_ROE*IFNLCT - (IFNLCT*C_ROE + (1-IFNLCT)*G/C_ROE)*NY

C.....Evaluate left eigenvectors at Roe averaged variables

      DEN = 1.0d0/(RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2)+RI(2,1)*RI(3,2)
     &                                                 -RI(3,1)*RI(2,2))

      LE(1,1) =  DEN*( RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2) )
      LE(1,2) =  DEN*RI(3,2)
      LE(1,3) = -DEN*RI(2,2)

      LE(2,1) = -DEN*( RI(2,1)*RI(3,3)-RI(2,3)*RI(3,1) )
      LE(2,2) =  DEN*( RI(3,3)-RI(3,1) )
      LE(2,3) = -DEN*( RI(2,3)-RI(2,1) )

      LE(3,1) =  DEN*( RI(2,1)*RI(3,2) - RI(2,2)*RI(3,1) )
      LE(3,2) = -DEN*RI(3,2)
      LE(3,3) =  DEN*RI(2,2)

C.....Compute the Roe matrix

      DO II = 1,3

        A_ROE(II,1) = RI(II,1)*ABS(EIGVAL(1))*LE(1,1) +
     &                RI(II,2)*ABS(EIGVAL(2))*LE(2,1) +
     &                RI(II,3)*ABS(EIGVAL(3))*LE(3,1)

        A_ROE(II,2) = RI(II,1)*ABS(EIGVAL(1))*LE(1,2) +
     &                RI(II,2)*ABS(EIGVAL(2))*LE(2,2) +
     &                RI(II,3)*ABS(EIGVAL(3))*LE(3,2)

        A_ROE(II,3) = RI(II,1)*ABS(EIGVAL(1))*LE(1,3) +
     &                RI(II,2)*ABS(EIGVAL(2))*LE(2,3) +
     &                RI(II,3)*ABS(EIGVAL(3))*LE(3,3)

      ENDDO

C.....Compute continuity fluxes at exterior state

      FX_EX = (IFNLCT*QX_EX + (1-IFNLCT)*QX_EX*HT_EX)*SFAC_EX
      FY_EX = IFNLCT*QY_EX + (1-IFNLCT)*QY_EX*HT_EX

C.....Compute x momentum fluxes at exterior state

      GX_EX = IFNLCT*(QX_EX*QX_EX/(HB_EX + ZE_EX) +
     &    0.5D0*G*(ZE_EX*ZE_EX + 2.D0*ZE_EX*HB_EX)) + (1-IFNLCT)*G*ZE_EX
      GX_EX = GX_EX*SFAC_EX
      GY_EX = IFNLCT*QX_EX*QY_EX/(HB_EX + ZE_EX)

C.....Compute y momentum fluxes at exterior state

      HX_EX = (IFNLCT*QX_EX*QY_EX/(HB_EX + ZE_EX))*SFAC_EX
      HY_EX = IFNLCT*(QY_EX*QY_EX/(HB_EX + ZE_EX) +
     &    0.5D0*G*(ZE_EX*ZE_EX + 2.D0*ZE_EX*HB_EX)) + (1-IFNLCT)*G*ZE_EX

C.....Compute continuity fluxes at interior state

      FX_IN = (IFNLCT*QX_IN + (1-IFNLCT)*QX_IN*HT_IN)*SFAC_IN
      FY_IN = IFNLCT*QY_IN + (1-IFNLCT)*QY_IN*HT_IN

C.....Compute x momentum fluxes at interior state

      GX_IN = IFNLCT*(QX_IN*QX_IN/(HB_IN + ZE_IN) +
     &    0.5D0*G*(ZE_IN*ZE_IN + 2.D0*ZE_IN*HB_IN)) + (1-IFNLCT)*G*ZE_IN
      GX_IN = GX_IN*SFAC_IN
      GY_IN = IFNLCT*QX_IN*QY_IN/(HB_IN + ZE_IN)

C.....Compute y momentum fluxes at interior state

      HX_IN = (IFNLCT*QX_IN*QY_IN/(HB_IN + ZE_IN))*SFAC_IN
      HY_IN = IFNLCT*(QY_IN*QY_IN/(HB_IN + ZE_IN) +
     &    0.5D0*G*(ZE_IN*ZE_IN + 2.D0*ZE_IN*HB_IN)) + (1-IFNLCT)*G*ZE_IN

C.....Compute the average flux in the normal direction

      F_AVG = 0.5D0*((FX_IN + FX_EX)*NX + (FY_IN + FY_EX)*NY)
      G_AVG = 0.5D0*((GX_IN + GX_EX)*NX + (GY_IN + GY_EX)*NY)
      H_AVG = 0.5D0*((HX_IN + HX_EX)*NX + (HY_IN + HY_EX)*NY)

C.....Compute the jump in the variables

      JUMP(1) = ZE_EX - ZE_IN
      JUMP(2) = QX_EX - QX_IN
      JUMP(3) = QY_EX - QY_IN

C.....Compute the Roe flux

      F_HAT = F_AVG - 0.5D0*(A_ROE(1,1)*JUMP(1) + A_ROE(1,2)*JUMP(2) +
     &                                            A_ROE(1,3)*JUMP(3) )
      G_HAT = G_AVG - 0.5D0*(A_ROE(2,1)*JUMP(1) + A_ROE(2,2)*JUMP(2) +
     &                                            A_ROE(2,3)*JUMP(3) )
      H_HAT = H_AVG - 0.5D0*(A_ROE(3,1)*JUMP(1) + A_ROE(3,2)*JUMP(2) +
     &                                            A_ROE(3,3)*JUMP(3) )
C      PRINT*,F_AVG,G_AVG,H_AVG
C      PRINT*,GX_IN, GX_EX
C      PRINT*,GY_IN, GY_EX
C      PRINT*,HX_IN, HX_EX
C      PRINT*,HY_IN, HY_EX
C      PRINT*,'---------------------------------------'

      RETURN
      END SUBROUTINE


C***********************************************************************
C
C     SUBROUTINE LLF_FLUX
C
C     This subroutine computes the Local Lax Friedrichs flux for the
C     shallow water equations.
C
C     Written by Ethan Kubatko (07-15-2005)
C
C***********************************************************************

      SUBROUTINE LLF_FLUX()

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE


C.....Declare local variables

      INTEGER II
      REAL(SZ) EIGVALS(6)
      REAL(SZ) EIGMAX, IDEPTH, HUU, HUV, HVV, GH2
      REAL(SZ) C_EX, C_IN, UN_EX, UN_IN
      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL
      REAL(SZ) FX1_IN, FX2_IN, FX3_IN, FY1_IN, FY2_IN, FY3_IN
      REAL(SZ) FX1_EX, FX2_EX, FX3_EX, FY1_EX, FY2_EX, FY3_EX
      REAL(SZ) F1_AVG, F2_AVG, F3_AVG, NX2, NY2, NXY
      REAL(SZ) TOL
      REAL(SZ) DEN

C.....Compute the jump in the variables

      JUMP(1) = ZE_EX - ZE_IN
      JUMP(2) = QX_EX - QX_IN
      JUMP(3) = QY_EX - QY_IN

C.....Compute the total height of the water column

      HT_IN = ZE_IN*NLEQ + HB_IN

C.....Compute continuity fluxes at interior state

      F1_NL  = NLEQ + LEQ*HB_IN
      FX1_IN = QX_IN*F1_NL*SFAC_IN
      FY1_IN = QY_IN*F1_NL

C.....Compute momentum flux terms at interior state

      FU_NL = NLEQ*QX_IN
      FV_NL = NLEQ*QY_IN
      FG_NL = NLEQG*ZE_IN
      FH_NL = 1.D0/(NLEQ*HT_IN + LEQ)
      U_IN  = QX_IN*FH_NL
      V_IN  = QY_IN*FH_NL

      HUU = FU_NL*U_IN
      HVV = FV_NL*V_IN
      HUV = FU_NL*V_IN
      GH2 = FG_NL*(0.5D0*ZE_IN + HB_IN) + FG_L*ZE_IN

C.....Compute x momentum fluxes at interior state

      FX2_IN = (HUU + GH2)*SFAC_IN
      FY2_IN = HUV

C.....Compute y momentum fluxes at interior state

      FX3_IN = HUV*SFAC_IN
      FY3_IN = HVV + GH2

C.....Compute the total height of the water column

      HT_EX = ZE_EX*NLEQ + HB_EX

C.....Compute continuity fluxes at exterior state

      F1_NL  = NLEQ + LEQ*HB_EX
      FX1_EX = QX_EX*F1_NL*SFAC_EX
      FY1_EX = QY_EX*F1_NL

C.....Compute momentum flux terms at interior state

      FU_NL = NLEQ*QX_EX
      FV_NL = NLEQ*QY_EX
      FG_NL = NLEQG*ZE_EX
      FH_NL = 1.D0/(NLEQ*HT_EX + LEQ)
      U_EX  = QX_EX*FH_NL
      V_EX  = QY_EX*FH_NL

      HUU = FU_NL*U_EX
      HVV = FV_NL*V_EX
      HUV = FU_NL*V_EX
      GH2 = FG_NL*(0.5D0*ZE_EX + HB_EX) + FG_L*ZE_EX

C.....Compute x momentum fluxes at exterior state

      FX2_EX = (HUU + GH2)*SFAC_EX
      FY2_EX = HUV

C.....Compute y momentum fluxes at exterior state

      FX3_EX = HUV*SFAC_EX
      FY3_EX = HVV + GH2

C.....Compute the average flux function

      F1_AVG = 0.5D0*((FX1_IN + FX1_EX)*NX + (FY1_IN + FY1_EX)*NY)
      F2_AVG = 0.5D0*((FX2_IN + FX2_EX)*NX + (FY2_IN + FY2_EX)*NY)
      F3_AVG = 0.5D0*((FX3_IN + FX3_EX)*NX + (FY3_IN + FY3_EX)*NY)

C.....Evaluate the eigenvalues at the interior and exterior states

      UN_IN = (U_IN*NX + V_IN*NY)*NLEQ
      C_IN = SQRT(G*HT_IN)
      EIGVALS(1) = ABS(UN_IN + C_IN)
      EIGVALS(2) = ABS(UN_IN)
      EIGVALS(3) = ABS(UN_IN - C_IN)

      UN_EX = (U_EX*NX + V_EX*NY)*NLEQ
      C_EX = SQRT(G*HT_EX)
      EIGVALS(4) = ABS(UN_EX + C_EX)
      EIGVALS(5) = ABS(UN_EX)
      EIGVALS(6) = ABS(UN_EX - C_EX)

C.....Find the maximum eigenvalue (in absolute value)

      EIGMAX = MAX( EIGVALS(1), EIGVALS(2), EIGVALS(3),
     &              EIGVALS(4), EIGVALS(5), EIGVALS(6) )

C.....Compute the Local Lax Friedrichs Fluxes

      NX2 = NLEQ+LEQ*NX*NX
      NY2 = NLEQ+LEQ*NY*NY
      NXY = NX*NY


      F_HAT = F1_AVG - 0.5D0*EIGMAX*(JUMP(1))
      G_HAT = F2_AVG - 0.5D0*EIGMAX*(JUMP(2)*NX2 + LEQ*JUMP(3)*NXY)
      H_HAT = F3_AVG - 0.5D0*EIGMAX*(JUMP(3)*NY2 + LEQ*JUMP(2)*NXY)

      RETURN
      END SUBROUTINE
      
C***********************************************************************
C
C     SUBROUTINE HLLC_FLUX()
C
C     This subroutine computes the Harten-Lax-van Leer Contact (HLLC)
C     flux for the shallow water equations.
C
C     Written by Ethan Kubatko
C
C***********************************************************************

      SUBROUTINE HLLC_FLUX()

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

C.....Declare local variables

      INTEGER II
      REAL(SZ) Q_STAR(3), Q(3)
      REAL(SZ) HUU, HUV, HVV, GH2
      REAL(SZ) C_EX, C_IN, UN_EX, UN_IN, UT_EX, UT_IN
      REAL(SZ) LS_EX, LS_IN, S_EX, S_IN, S_STAR
      REAL(SZ) H_STAR
      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL
      REAL(SZ) FX1_IN, FX2_IN, FX3_IN, FY1_IN, FY2_IN, FY3_IN
      REAL(SZ) FX1_EX, FX2_EX, FX3_EX, FY1_EX, FY2_EX, FY3_EX

C.....Compute the interior and exterior water column heights

      HT_IN = ZE_IN*NLEQ + HB_IN
      HT_EX = ZE_EX*NLEQ + HB_EX

C.....Compute the interior and exterior normal velocities

      FH_NL = 1.D0/(NLEQ*HT_IN + LEQ)
      U_IN  = QX_IN*FH_NL
      V_IN  = QY_IN*FH_NL
      UN_IN =  U_IN*NX + V_IN*NY
      UT_IN = -U_IN*NY + V_IN*NX

      FH_NL = 1.D0/(NLEQ*HT_EX + LEQ)
      U_EX  = QX_EX*FH_NL
      V_EX  = QY_EX*FH_NL
      UN_EX =  U_EX*NX + V_EX*NY
      UT_EX = -U_EX*NY + V_EX*NX

C.....Compute H*

      C_EX = SQRT(G*HT_EX)
      C_IN = SQRT(G*HT_IN)
      H_STAR = 0.5D0*(C_IN + C_EX) + 0.25D0*(UN_IN - UN_EX)
      H_STAR = 1.D0/G*H_STAR*H_STAR

C.....Compute the wave speeds

      LS_EX = 1.D0
      LS_IN = 1.D0
      IF (H_STAR.GT.HT_EX) LS_EX = SQRT( (H_STAR*H_STAR + H_STAR*HT_EX)
     &                                             /(2.D0*HT_EX*HT_EX) )
      IF (H_STAR.GT.HT_IN) LS_IN = SQRT( (H_STAR*H_STAR + H_STAR*HT_IN)
     &                                             /(2.D0*HT_IN*HT_IN) )

      S_IN   = UN_IN - C_IN*LS_IN
      S_EX   = UN_EX + C_EX*LS_EX
      S_STAR = ( S_IN*HT_EX*(UN_EX - S_EX) - S_EX*HT_IN*(UN_IN - S_IN) )
     &                  /( HT_EX*(UN_EX - S_EX) - HT_IN*(UN_IN - S_IN) )

C.....Compute the numerical flux based on the wave speeds

      IF (S_IN.GE.0) THEN

C.......Compute continuity fluxes at exterior state

        F1_NL  = NLEQ + LEQ*HB_IN
        FX1_IN = QX_IN*F1_NL
        FY1_IN = QY_IN*F1_NL

C.......Compute momentum flux terms at exterior state

        FU_NL = NLEQ*QX_IN
        FV_NL = NLEQ*QY_IN
        FG_NL = NLEQG*ZE_IN

        HUU = FU_NL*U_IN
        HVV = FV_NL*V_IN
        HUV = FU_NL*V_IN
        GH2 = FG_NL*(0.5D0*ZE_IN + HB_IN) + FG_L*ZE_IN

C.......Compute x momentum fluxes at exterior state

        FX2_IN = HUU + GH2
        FY2_IN = HUV

C.......Compute y momentum fluxes at exterior state

        FX3_IN = HUV
        FY3_IN = HVV + GH2

C.......Set the numerical fluxes equal to the exterior fluxes

        F_HAT = FX1_IN*NX + FY1_IN*NY
        G_HAT = FX2_IN*NX + FY2_IN*NY
        H_HAT = FX3_IN*NX + FY3_IN*NY

      ELSEIF ((S_IN.LE.0).AND.(S_STAR.GE.0)) THEN

C.......Compute continuity fluxes at exterior state

        F1_NL  = NLEQ + LEQ*HB_IN
        FX1_IN = QX_IN*F1_NL
        FY1_IN = QY_IN*F1_NL

C.......Compute momentum flux terms at exterior state

        FU_NL = NLEQ*QX_IN
        FV_NL = NLEQ*QY_IN
        FG_NL = NLEQG*ZE_IN

        HUU = FU_NL*U_IN
        HVV = FV_NL*V_IN
        HUV = FU_NL*V_IN
        GH2 = FG_NL*(0.5D0*ZE_IN + HB_IN) + FG_L*ZE_IN

C.......Compute x momentum fluxes at exterior state

        FX2_IN = HUU + GH2
        FY2_IN = HUV

C.......Compute y momentum fluxes at exterior state

        FX3_IN = HUV
        FY3_IN = HVV + GH2

C.......Compute Q*

        Q_STAR(1) = HT_IN*((S_IN - UN_IN)/(S_IN - S_STAR))
        Q_STAR(2) = Q_STAR(1)*S_STAR
        Q_STAR(3) = Q_STAR(1)*UT_IN
        Q(1) =  HT_IN
        Q(2) =  QX_IN*NX + QY_IN*NY
        Q(3) = -QX_IN*NY + QY_IN*NX

C.......Compute the numerical flux

        F_HAT = FX1_IN*NX + FY1_IN*NY + S_IN*(  Q_STAR(1) - Q(1))
        G_HAT = FX2_IN*NX + FY2_IN*NY + S_IN*( (Q_STAR(2) - Q(2))*NX -
     &                                         (Q_STAR(3) - Q(3))*NY )
        H_HAT = FX3_IN*NX + FY3_IN*NY + S_IN*( (Q_STAR(2) - Q(2))*NY +
     &                                         (Q_STAR(3) - Q(3))*NX )

      ELSEIF ((S_STAR.LE.0).AND.(S_EX.GE.0)) THEN

C.......Compute continuity fluxes at exterior state

        F1_NL  = NLEQ + LEQ*HB_EX
        FX1_EX = QX_EX*F1_NL
        FY1_EX = QY_EX*F1_NL

C.......Compute momentum flux terms at exterior state

        FU_NL = NLEQ*QX_EX
        FV_NL = NLEQ*QY_EX
        FG_NL = NLEQG*ZE_EX

        HUU = FU_NL*U_EX
        HVV = FV_NL*V_EX
        HUV = FU_NL*V_EX
        GH2 = FG_NL*(0.5D0*ZE_EX + HB_EX) + FG_L*ZE_EX

C.......Compute x momentum fluxes at exterior state

        FX2_EX = HUU + GH2
        FY2_EX = HUV

C.......Compute y momentum fluxes at exterior state

        FX3_EX = HUV
        FY3_EX = HVV + GH2

C.......Compute Q*

        Q_STAR(1) = HT_EX*((S_EX - UN_EX)/(S_EX - S_STAR))
        Q_STAR(2) = Q_STAR(1)*S_STAR
        Q_STAR(3) = Q_STAR(1)*UT_EX
        Q(1) =  HT_EX
        Q(2) =  QX_EX*NX + QY_EX*NY
        Q(3) = -QX_EX*NY + QY_EX*NX

C.......Compute the numerical flux

        F_HAT = FX1_EX*NX + FY1_EX*NY + S_EX*(  Q_STAR(1) - Q(1))
        G_HAT = FX2_EX*NX + FY2_EX*NY + S_EX*( (Q_STAR(2) - Q(2))*NX -
     &                                         (Q_STAR(3) - Q(3))*NY )
        H_HAT = FX3_EX*NX + FY3_EX*NY + S_EX*( (Q_STAR(2) - Q(2))*NY +
     &                                         (Q_STAR(3) - Q(3))*NX )

      ELSEIF (S_EX.LE.0) THEN

C.......Compute continuity fluxes at exterior state

        F1_NL  = NLEQ + LEQ*HB_EX
        FX1_EX = QX_EX*F1_NL
        FY1_EX = QY_EX*F1_NL

C.......Compute momentum flux terms at exterior state

        FU_NL = NLEQ*QX_EX
        FV_NL = NLEQ*QY_EX
        FG_NL = NLEQG*ZE_EX

        HUU = FU_NL*U_EX
        HVV = FV_NL*V_EX
        HUV = FU_NL*V_EX
        GH2 = FG_NL*(0.5D0*ZE_EX + HB_EX) + FG_L*ZE_EX

C.......Compute x momentum fluxes at exterior state

        FX2_EX = HUU + GH2
        FY2_EX = HUV

C.......Compute y momentum fluxes at exterior state

        FX3_EX = HUV
        FY3_EX = HVV + GH2

C.......Set the numerical fluxes equal to the exterior fluxes

        F_HAT = FX1_EX*NX + FY1_EX*NY
        G_HAT = FX2_EX*NX + FY2_EX*NY
        H_HAT = FX3_EX*NX + FY3_EX*NY

      ENDIF

      RETURN
      END SUBROUTINE HLLC_FLUX
      
      
