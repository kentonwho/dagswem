C***********************************************************************
C
C     SUBROUTINE LUND_FORMULA( )
C
C     This subroutine calculates the sediment flux based on the Lund
C     formula.
C
C     Written by Ethan Kubatko (03-05-2005)
C
C***********************************************************************

      SUBROUTINE LUND_FORMULA(HT_GP, U_GP, V_GP, HB_GP, QSTOT_X,QSTOT_Y)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

C.....Declare local variables

      INTEGER II
      REAL(SZ) HT_GP, U_GP, V_GP, HB_GP, QSTOT_X, QSTOT_Y

C.....Compute friction factor for current

      FC = 2.D0*(KAPPA/(1.D0 + LOG(Z0/HT_GP)))**2.D0
      
C.....Compute the magnitude of the current velocity

      UC = SQRT(U_GP**2.D0 + V_GP**2.D0)
      
C.....Compute shear stress from current

      TAUC = 0.5D0*RHOWAT0*FC*UC**2.D0

C.....If waves are included calculate sediment tranpsort due to waves

      IF ((SEDFLAG.EQ.2).AND.(H_WAVE.GT.1.E-10)) THEN

C.....Compute the wave length using a Pade approximation

        Y_WAVE = (2.D0*PI/T_WAVE)**2.D0*HT_GP/G

        F_WAVE = Y_WAVE + 1.D0/(1.D0 + Y_WAVE*(0.66667D0 +
     &           Y_WAVE*(0.3555D0    + Y_WAVE*(0.16084D0 +
     &           Y_WAVE*(0.06320D0   + Y_WAVE*(0.02174D0 +
     &           Y_WAVE*(0.00654D0   + Y_WAVE*(0.00171D0 +
     &           Y_WAVE*(0.00039D0   + Y_WAVE*0.000111D0)))))))))

        L_WAVE = 2*PI/SQRT(Y_WAVE*F_WAVE/HT_GP**2.D0)

C.....Compute the near-bed wave semi-excursion distance

        AW = 0.5D0*H_WAVE/SINH(2.D0*PI*HT_GP/L_WAVE)

C.....Compute the magnitude of the orbital wave speed

        UW = 0.5D0*H_WAVE*G*T_WAVE/(L_WAVE*COSH(2.D0*PI*HT_GP/L_WAVE))

C.....Compute the bottom friction due to waves

        IF (AW.EQ.0.D0) THEN
          FW = 0.D0
        ELSEIF ((AW.NE.0).AND.(AW.LE.AW_LIM)) THEN
          FW = 0.30
        ELSE
          FW = EXP(5.2D0*(KS/AW)**0.194D0 - 5.98D0)
        ENDIF

C.....Compute the bottom shear stress from waves (mean value)

        TAUW = 0.25D0*RHOWAT0*FW*UW**2.D0
        
C.....Compute the maximum wave shear stress value during a cycle

        TAUWM = 0.5D0*RHOWAT0*FW*UW**2.D0

C.....Compute the angle (in radians) between the wave and current

        WAVE_X = COS(A_WAVE*DEG2RAD)
        WAVE_Y = SIN(A_WAVE*DEG2RAD)

        DOT = WAVE_X*U_GP + WAVE_Y*V_GP

        DOT_U_MAG = DOT/UC

        IF (DOT_U_MAG.GT.+1) DOT_U_MAG = +1.D0
        IF (DOT_U_MAG.LT.-1) DOT_U_MAG = -1.D0

        PHI0 = ACOS(DOT_U_MAG)

C.....Compute combined shear stress from waves and currents (mean value)

        TAUCW = SQRT((TAUC + TAUW*COS(PHI0))**2.D0
     &                    + (TAUW*SIN(PHI0))**2.D0)
     
C.....Compute the maximum combined shear stress during a cycle

        TAUCWM = SQRT((TAUC + TAUWM*COS(PHI0))**2.D0
     &                     + (TAUWM*SIN(PHI0))**2.D0)
     
C.....Compute the energy dissipation due to the waves

        DW = 2.D0/3.D0/PI*RHOWAT0*FW*UW**3.D0

C.....If waves are not included set parameters accordingly
        
      ELSE
      
        TAUCW   = TAUC
        TAUCWM  = TAUC
        DW = 0.D0
        DB = 0.D0
        
      ENDIF
     
C.....Compute the energy dissipation due to the current

      DC = 0.5D0*RHOWAT0*FC*UC**3.D0

C.....Compute the total energy dissipation (including wave breaking
C.....energy dissipation, DB from ST-WAVE)
        
      DTOT = DC + DW + DB

C.....Compute the sediment decay parameter (non-dimensional)

      IF (DTOT.LE.DTOT_MIN) THEN
        ADIFF = 20.D0
      ELSE
        ADIFF = WFV/CD/(DTOT/RHOWAT0)**(1.D0/3.D0)
      ENDIF


C.....Compute reduction coefficients for suspended load

      CEXP = 1.D0 - EXP(-ADIFF)

      IF (ADIFF.LT.1.D0) THEN
        CVEL = 1.D0 - 0.0343D0*ADIFF
      ELSE
        FIG = ADIFF**CM*EXP(-ADIFF)*
     &                    (1.D0 + CM/ADIFF + CM*(CM - 1.D0)/ADIFF**2.D0)
        CVEL = (CM + 1.D0)*(FG - FIG)/(1.D0 - EXP(-ADIFF))/ADIFF**CM
      ENDIF
        
C.....Compute the reference concentration

      IF (TAUCW.GT.TAULIM) THEN
        CRM = CSL*TAULIM/TAUCR*EXP(-B_EMP*TAUCR/TAULIM)
        CR  = CRM*EXP(-ADEC*(TAUCW - TAULIM))
      ELSE
        CR = CSL*TAUCW/TAUCR*EXP(-B_EMP*TAUCR/TAUCWM)
      ENDIF
      
C.....Compute the magnitude of the bed load transport

      QSBED = CBL*SQRT(TAUC/RHOWAT0)*TAUCW*EXP(-B_EMP*TAUCR/TAUCWM)

C.....Compute the magnitude of the suspended load transport

      QSSUS = CR*HT_GP*UC/ADIFF*CEXP*CVEL

C.....Compute the magnitude of the total sediment transport

      QSTOT = MAG_SED*(QSBED + QSSUS)
      
C.....Compute the X and Y components of the sediment flux
      
      QSTOT_X = QSTOT*U_GP/UC
      QSTOT_Y = QSTOT*V_GP/UC
      
      RETURN
      END SUBROUTINE
      
      
      
      
        

      
