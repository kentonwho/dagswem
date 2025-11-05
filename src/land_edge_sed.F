C***********************************************************************
C
C     SUBROUTINE LAND_EDGE_SED( )
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
C     Written by Ethan Kubatko (03-05-2005)
C
C***********************************************************************


      SUBROUTINE LAND_EDGE_SED(L)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED
      REAL(SZ) TX,TY
      
C.....Retrieve the global and local edge number

      GED = NLEDN(L)
      LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge

      EL_IN = NEDEL(1,GED)

C.....Retrieve the components of the normal vector to the edge
        
      NX = COSNX(GED)
      NY = SINNX(GED)
      
C.....Set the components for the tangential vector to the edge
      
      TX = -NY
      TY =  NX

C.....Compute ZE, QX, QY, and HB at each Gauss point

      DO 111 I = 1,NEGP

        ZE_IN = 0.D0
        QX_IN = 0.D0
        QY_IN = 0.D0
        
        HB_IN = 0.D0

C.....Compute the solution at the interior state

        DO K = 1,DOF
          
          ZE_IN = ZE_IN + ZE(K,EL_IN,1  )*PHI_EDGE(K,I,LED)
          QX_IN = QX_IN + QX(K,EL_IN,1  )*PHI_EDGE(K,I,LED)
          QY_IN = QY_IN + QY(K,EL_IN,1  )*PHI_EDGE(K,I,LED)
          HB_IN = HB_IN + HB(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)

        ENDDO
        
C.....If sediment tranpsort due to waves is on compute the wave period,
C.....height, and angle at the given interior gauss point

       IF (SEDFLAG.EQ.2) THEN

          WAVE_RATIO = (XEGP(I) - (-1.D0))/2.D0

          T_WAVE = WAVE_T1(N1) + WAVE_RATIO*(WAVE_T2(N2) - WAVE_T1(N1))
          A_WAVE = WAVE_A1(N1) + WAVE_RATIO*(WAVE_A2(N2) - WAVE_A1(N1))
          H_WAVE = RAMPDG*(WAVE_H1(N1)
     &                         + WAVE_RATIO*(WAVE_H2(N2) - WAVE_H1(N1)))

        ENDIF
        
C.....Compute the velocity in the normal and tangental direction
        
C        Q_N_INT = QX_IN*NX + QY_IN*NY
C        Q_T_INT = QX_IN*TX + QY_IN*TY

C.....Reflect the velocity in the normal direction

C        Q_N_EXT = -Q_N_INT
C        Q_T_EXT =  Q_T_INT
        
C.....Compute the x and y components of the external state flow

C        QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)/(NX*TY - NY*TX)
C        QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)/(NX*TY - NY*TX)

        ZE_EX = ZE_IN
        HB_EX = HB_IN
        QX_EX = 0.D0
        QY_EX = 0.D0

C.....Compute the total height of the water column

        HT_IN = ZE_IN*IFNLFA + HB_IN
        HT_EX = ZE_EX*IFNLFA + HB_EX

C.....Compute the velocities in the x and y directions

        U_IN = QX_IN/HT_IN
        U_EX = QX_EX/HT_EX

        V_IN = QY_IN/HT_IN
        V_EX = QY_EX/HT_EX

C.....Compute the magnitudes of the velocities

        UMAG_IN = SQRT(U_IN**(2.D0) + V_IN**(2.D0))
        UMAG_EX = SQRT(U_EX**(2.D0) + V_EX**(2.D0))

C.....Compute the Roe averaged velocities

        U_ROE = (U_IN*SQRT(HT_IN) + U_EX*SQRT(HT_EX))/
     &                                     (SQRT(HT_IN) + SQRT(HT_EX))
        V_ROE = (V_IN*SQRT(HT_IN) + V_EX*SQRT(HT_EX))/
     &                                     (SQRT(HT_IN) + SQRT(HT_EX))

C.....If the roe velocity vector and normal > 0 use the interior element
C.....values to compute the sediment transport

        IF (((U_ROE*NX + V_ROE*NY).GT.0).AND.(UMAG_IN.NE.0)) THEN

          CALL LUND_FORMULA(HT_IN, U_IN, V_IN, HB_IN, QSX, QSY)

C.....If the roe velocity vector and normal < 0 use the interior element
C.....values to compute the sediment transport

        ELSEIF (((U_ROE*NX + V_ROE*NY).LT.0).AND.(UMAG_EX.NE.0)) THEN

          CALL LUND_FORMULA(HT_EX, U_EX, V_EX, HB_EX, QSX, QSY)

C.....Else skip the sediment transport and edge integral calculations

        ELSE
          GOTO 111
        ENDIF

C.....Compute the sediment flux normal to the edge

        QS_HAT = QSX*NX + QSY*NY

C.....Compute the edge integral

        DO K = 1,DOF
          CALL EDGE_INT_SED(EL_IN, LED, GED, I, QS_HAT)
        ENDDO

111   CONTINUE

      RETURN
      END SUBROUTINE


