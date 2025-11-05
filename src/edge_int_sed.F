C***********************************************************************
C
C     SUBROUTINE EDGE_INT_SED()
C
C     This subroutine computes the edge integrals for the DG sed using
C     Gauss quadrature and adds them to the RHS.
C
C     Written by Ethan Kubatko (03-05-2005)
C
C***********************************************************************

      SUBROUTINE EDGE_INT_SED(EL_I, LED, GED, GP, FLUX1)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER EL_I, LED, GED, GP
      REAL(SZ) AREA, IMASS, FLUX1
      
C.....Retrieve the element area

      AREA = 0.5D0*AREAS(EL_I)

C.....Comput the edge integral

      IMASS = M_INV(K)/(0.5D0*AREA)

      RHS_HB(K,EL_I,IRK) = RHS_HB(K,EL_I,IRK)
     &          - IMASS*XLEN(GED)/2.D0*FLUX1*PHI_EDGE(K,GP,LED)*WEGP(GP)
      
      RETURN
      END SUBROUTINE
