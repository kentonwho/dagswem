C***********************************************************************
C
C     SUBROUTINE EDGE_INT_HYDRO()
C
C     This subroutine computes the edge integrals for the DG hydro using
C     Gauss quadrature and adds them to the RHS.
C
C     Written by Ethan Kubatko (03-05-2005)
C
C     02-23-2007 - sb - rewritten for better performance
C
C***********************************************************************

#if 0
      SUBROUTINE EDGE_INT_HYDRO(EL_I,LED,GED,GP,FLUX1,FLUX2,FLUX3)
      ! <ezpp-noinst>
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER EL_I, LED, GED, GP
      REAL(SZ) FLUX1, FLUX2, FLUX3, W
      
C.....Comput the edge integral

      W = 2.0*M_INV(K)/AREAS(EL_I)*XLEN(GED)*PHI_EDGE(K,GP,LED)
     &     *WEGP(GP)

      RHS_ZE(K,EL_I,IRK) = RHS_ZE(K,EL_I,IRK) - W*FLUX1
      RHS_QX(K,EL_I,IRK) = RHS_QX(K,EL_I,IRK) - W*FLUX2
      RHS_QY(K,EL_I,IRK) = RHS_QY(K,EL_I,IRK) - W*FLUX3
     
      RETURN
      END SUBROUTINE
#else
      SUBROUTINE EDGE_INT_HYDRO(EL_I,LED,GED,GP,FLUX1,FLUX2,FLUX3)
      ! <ezpp-noinst>
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER EL_I, LED, GED, GP
      REAL(SZ) AREA, IMASS, FLUX1, FLUX2, FLUX3, W
      
C.....Retrieve the element area
      
      AREA = 0.5D0*AREAS(EL_I)

C.....Comput the edge integral

      IMASS = M_INV(K)/(0.5D0*AREA)

      RHS_ZE(K,EL_I,IRK) = RHS_ZE(K,EL_I,IRK)
     &          - IMASS*XLEN(GED)/2.D0*FLUX1*PHI_EDGE(K,GP,LED)*WEGP(GP)
      RHS_QX(K,EL_I,IRK) = RHS_QX(K,EL_I,IRK)
     &          - IMASS*XLEN(GED)/2.D0*FLUX2*PHI_EDGE(K,GP,LED)*WEGP(GP)
      RHS_QY(K,EL_I,IRK) = RHS_QY(K,EL_I,IRK)
     &          - IMASS*XLEN(GED)/2.D0*FLUX3*PHI_EDGE(K,GP,LED)*WEGP(GP)
     
      RETURN
      END SUBROUTINE
#endif
