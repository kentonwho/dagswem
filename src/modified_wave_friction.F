C***********************************************************************
C
C     SUBROUTINE MODIFIED_WAVE_FRICTION()
C
C     UNDER CONSTRUCTION!!!!
C
C     Written by Ethan Kubatko (06-11-2004)
C
C***********************************************************************

      SUBROUTINE MODIFIED_WAVE_FRICTION()
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

C.....Compute the wave length using a Pade approximation (this part is
C.....from Magnus' code)

      Y_WAVE = (2.D0*PI/T_WAVE)**2.D0*(HB_IN + ZE_IN)/G
      
      F_WAVE = Y_WAVE + 1.D0/(1.D0 + Y_WAVE*(0.66667D0 +
     &         Y_WAVE*(0.3555D0    + Y_WAVE*(0.16084D0 +
     &         Y_WAVE*(0.06320D0   + Y_WAVE*(0.02174D0 +
     &         Y_WAVE*(0.00654D0   + Y_WAVE*(0.00171D0 +
     &         Y_WAVE*(0.00039D0   + Y_WAVE*0.000111D0)))))))))

      L_WAVE = 2.D0*PI/SQRT(Y_WAVE*F_WAVE/(DP(I) + ETA1(I))**2.D0)

C.....Compute the orbital wave velocity

      WB = 2.D0*PI/WAVE_T*WAVE_H/(PI*SINH(2.D0*PI/
     &                                          L_WAVE*(HB_IN + ZE_IN)))

C.....Compute the modified friction factors due to waves

      COSA1 = COS(WAVE_A(I))
      SINA1 = SIN(WAVE_A(I))
      UVW1  = UU1(I)*UU1(I) + VV1(I)*VV1(I) + WB(I)*WB(I)
      UVW2  = 2.D0*(UU1(I)*COSA1 + VV1(I)*SINA1)*WB(I)
      WD1 = 0.5D0*((UVW1 + UVW2)**0.5D0 + (UVW1 - UVW2)**0.5D0)

      WD = 0
      IF (ABS(WD1).GE.(0.001)) THEN
        WD = WB(I)*WB(I)/WD1
      ENDIF

C              WDXX = WD*COSA1*COSA1
C              WDYY = WD*SINA1*SINA1
C              WDXY = WD*SINA1*COSA1
C              CHYBR = IFNLBF + IFHYBF*
C     &                         (1+(HBREAK/HH1)**FTHETA)**(FGAMMA/FTHETA)
C
C              TKXX(I) = FRIC(I)*(IFLINBF + ((WD1 + WDXX)/HH1)*CHYBR)
C              TKYY(I) = FRIC(I)*(IFLINBF + ((WD1 + WDYY)/HH1)*CHYBR)
C              TKXY(I) = FRIC(I)*(IFLINBF + WDXY/HH1*CHYBR)
C
C            ENDIF
C
C          endif

      RETURN
      END SUBROUTINE
