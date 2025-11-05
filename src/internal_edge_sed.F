C***********************************************************************
C
C     SUBROUTINE INTERNAL_EDGE_SED( )
C
C     This subroutine does the following:
C
C       1.  Calculates the values of the necessary variables at the edge
C           gauss points for INTERNAL edges
C       2.  Calls the appropriate subroutine to compute the flux at
C           these points.
C       3.  Calls the appropriate subroutine to compute the boundary
C           integrals.
C
C     Written by Ethan Kubatko (03-05-2005)
C
C-----------------------------------------------------------------------
C
C     Modification history on hp_DG_ADCIRC after v9_sb2.2.1
C     
C     v9_sb2.2.1 - 08/16/05 - sb - Hard bottom implementation
C
C***********************************************************************

      SUBROUTINE INTERNAL_EDGE_SED(L)
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED
      
      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER L, LED_IN, LED_EX, GED, GP_IN, GP_EX
Csb-HB..added for hard bottom
      REAL(SZ) LDP_HB(2), GDP_HB, W1, W2

C.....Retrieve the global and local edge number

      GED = NIEDN(L)
      LED_IN = NEDSD(1,GED)
      LED_EX = NEDSD(2,GED)

C.....Retrieve the elements which share the edge

      EL_IN = NEDEL(1,GED)
      EL_EX = NEDEL(2,GED)
      
C.....Retrieve the nodes of the given element

      N1 = NEDNO(1,GED)
      N2 = NEDNO(2,GED)
      
C.....Retrieve the components of the normal vector to the edge
        
      NX = COSNX(GED)
      NY = SINNX(GED)

Csb-HB Retrieve the hard bottom depth at the end nodes of the edge
      
      IF(HBFLAG.GT.0) THEN
        LDP_HB(1) = DP_HB(N1)
        LDP_HB(2) = DP_HB(N2)
      ENDIF
      
C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

      DO 111 I = 1,NEGP

        GP_IN = I
        GP_EX = NEGP - I + 1

        ZE_IN = 0.D0
        QX_IN = 0.D0
        QY_IN = 0.D0
        HB_IN = 0.D0

        ZE_EX = 0.D0
        QX_EX = 0.D0
        QY_EX = 0.D0
        HB_EX = 0.D0

        DO K = 1,DOF

          ZE_IN = ZE_IN + ZE(K,EL_IN,1  )*PHI_EDGE(K,GP_IN,LED_IN)
          QX_IN = QX_IN + QX(K,EL_IN,1  )*PHI_EDGE(K,GP_IN,LED_IN)
          QY_IN = QY_IN + QY(K,EL_IN,1  )*PHI_EDGE(K,GP_IN,LED_IN)
          HB_IN = HB_IN + HB(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN)

          ZE_EX = ZE_EX + ZE(K,EL_EX,1  )*PHI_EDGE(K,GP_EX,LED_EX)
          QX_EX = QX_EX + QX(K,EL_EX,1  )*PHI_EDGE(K,GP_EX,LED_EX)
          QY_EX = QY_EX + QY(K,EL_EX,1  )*PHI_EDGE(K,GP_EX,LED_EX)
          HB_EX = HB_EX + HB(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX)

        ENDDO
        
C.....If sediment tranpsort due to waves is on compute the wave period,
C.....height, and angle at the given interior gauss point

       IF (SEDFLAG.EQ.2) THEN
          WAVE_RATIO = (XEGP(I) - (-1.D0))/2.D0
          
          T_WAVE = WAVE_T(N1) + WAVE_RATIO*(WAVE_T(N2) - WAVE_T(N1))
          A_WAVE = WAVE_A(N1) + WAVE_RATIO*(WAVE_A(N2) - WAVE_A(N1))
          H_WAVE = RAMPDG*(WAVE_H(N1)
     &                           + WAVE_RATIO*(WAVE_H(N2) - WAVE_H(N1)))
        ENDIF

C.....Compute the total height of the water column

        HT_IN = ZE_IN*IFNLFA + HB_IN
        HT_EX = ZE_EX*IFNLFA + HB_EX

C.....Compute the velocities in the x and y directions

        U_IN = QX_IN/HT_IN
        U_EX = QX_EX/HT_EX

        V_IN = QY_IN/HT_IN
        V_EX = QY_EX/HT_EX

C.....Compute the magnitudes of the velocities

        UMAG_IN = SQRT(U_IN**2.D0 + V_IN**2.D0)
        UMAG_EX = SQRT(U_EX**2.D0 + V_EX**2.D0)

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

C.....Limit the outgoing flux for the hard bottom

        IF(HBFLAG.GT.0) THEN
          W2 = (XEGP(I) + 1.D0)/2.D0
          W1 = 1.D0 - W2
          GDP_HB = W1*LDP_HB(1) + W2*LDP_HB(2)

          IF(QS_HAT.GE.0.D0) THEN
            IF(HB_IN.GE.GDP_HB) QS_HAT = 0.D0
          ELSE 
            IF(HB_EX.GE.GDP_HB) QS_HAT = 0.D0
          ENDIF

          IF(WDFLG(EL_IN).EQ.0.OR.WDFLG(EL_EX).EQ.0) THEN
            QS_HAT = 0.D0
          ENDIF
        ENDIF

C.....Compute the edge integral

        DO K = 1,DOF
          CALL EDGE_INT_SED(EL_IN, LED_IN, GED, GP_IN,  QS_HAT)
          CALL EDGE_INT_SED(EL_EX, LED_EX, GED, GP_EX, -QS_HAT)
        ENDDO

111   CONTINUE
       
      RETURN
      END SUBROUTINE
