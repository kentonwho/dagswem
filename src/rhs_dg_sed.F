C***********************************************************************
C
C     SUBROUTINE RHS_DG_SED()
C
C     This subroutine computes the area integrals for the DG sed and
C     adds them into the RHS.
C
C     Written by Ethan Kubatko (03-05-2005)
C
C***********************************************************************

      SUBROUTINE RHS_DG_SED(L)
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER L
      REAL(SZ) AREA, IMASS
      
C.....Compute the area of the element
      
      AREA = 0.5D0*AREAS(L)
      
C.....Retrieve the global node numbers for the element

      N1 = NM(L,1)
      N2 = NM(L,2)
      N3 = NM(L,3)
      
C.....Compute ZE, QX, QY, and HB at each area Gauss quadrature point

      DO 111 I = 1,NAGP
      
        ZE_IN = 0.D0
        QX_IN = 0.D0
        QY_IN = 0.D0
        HB_IN = 0.D0

        DO K = 1,DOF                 

          ZE_IN = ZE_IN + ZE(K,L,1  )*PHI_AREA(K,I)
          QX_IN = QX_IN + QX(K,L,1  )*PHI_AREA(K,I)
          QY_IN = QY_IN + QY(K,L,1  )*PHI_AREA(K,I)
          HB_IN = HB_IN + HB(K,L,IRK)*PHI_AREA(K,I)

        ENDDO
        
C.....Compute the total height of the water column

        HT_IN = ZE_IN*IFNLFA + HB_IN
        
C.....Compute the velocities

        U_IN = QX_IN/HT_IN
        V_IN = QY_IN/HT_IN
        
C.....Compute the magnitude of the velocity

        UMAG_IN = SQRT(U_IN**2.D0 + V_IN**2.D0)
        
C.....Compute the sediment fluxes

        IF (UMAG_IN.NE.0) THEN
          CALL LUND_FORMULA(HT_IN, U_IN, V_IN, HB_IN, QSX_IN, QSY_IN)
        ELSE
          GOTO 111
        ENDIF

C.....Build the rhs
        
        DO K = 1,DOF
        
          IMASS = M_INV(K)/(0.5D0*AREA)

          RHS_HB(K,L,IRK) = RHS_HB(K,L,IRK) + IMASS*0.5D0*AREA
     &       *(QSX_IN*(DRPHI(K,I)*DRDX(L) + DSPHI(K,I)*DSDX(L))
     &       + QSY_IN*(DRPHI(K,I)*DRDY(L) + DSPHI(K,I)*DSDY(L)))*WAGP(I)

        ENDDO
        
111   CONTINUE
      
      RETURN
      END SUBROUTINE
