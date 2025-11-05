C***********************************************************************
C
C     SUBROUTINE COMPUTE_LDG_Z(IT)
C
C     Compute variable Z to be used in the LDG terms.
C
C     Written by Shintaro Bunya (01-01-2007)
C
C***********************************************************************

      SUBROUTINE COMPUTE_LDG_Z(IT)

C.....Use appropriate modules
      
      USE GLOBAL
      USE DG
      USE SED

c#ifdef CMPI
c      USE MESSENGER_ELEM
c#endif

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER IT,L,GED,NBOREL,NNBORS,NDRYNBORS
      INTEGER Detected
      REAL(SZ) QBCT1,QBCT2,QBC1,QBC2,ZP(3)
      REAL(SZ) DPAVG

C.....Initialize

      DO L=1,NE
         DO K=1,DOF
            LZ(K,1,1,L) = 0.d0
            LZ(K,1,2,L) = 0.d0
            LZ(K,2,1,L) = 0.d0
            LZ(K,2,2,L) = 0.d0
         ENDDO
      ENDDO

C.....Loop over edges

      

C.....Loop over elements

      

C.....Solve for Z      

      RETURN
      END SUBROUTINE


C***********************************************************************
C
C     SUBROUTINE INTERNAL_EDGE_LDG( )
C
C     This subroutine does the following:
C
C       1.  Calculates the values of the necessary variables at the edge
C           gauss points for INTERNAL edges
C       2.  Take the average of these values
C       3.  Perform boundary integration
C
C     Written by Shintaro Bunya (01-02-2007)
C
C***********************************************************************

      SUBROUTINE INTERNAL_EDGE_LDG()

C.....Use appropriate modules
      
      USE SIZE, ONLY : SZ
      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED_IN, LED_EX, GED, GP_IN, GP_EX
      REAL(SZ) QX_AVG, QY_AVG

      DO 1000 L = 1,NIEDS

C.....Retrieve the global and local edge number

        GED = NIEDN(L)
        LED_IN = NEDSD(1,GED)
        LED_EX = NEDSD(2,GED)

C.....Retrieve the elements which share the edge

        EL_IN = NEDEL(1,GED)
        EL_EX = NEDEL(2,GED)

C.....Retrieve the components of the normal vector to the edge

        NX = COSNX(GED)
        NY = SINNX(GED)

C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

        DO I = 1,NEGP

          GP_IN = I
          GP_EX = NEGP - I + 1

          QX_IN = 0.D0
          QY_IN = 0.D0

          QX_EX = 0.D0
          QY_EX = 0.D0

          DO K = 1,DOF

            QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN)
            QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN)

            QX_EX = QX_EX + QX(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX)
            QY_EX = QY_EX + QY(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX)

          ENDDO

C.....Take the average

          QX_AVG = 0.5*(QX_IN + QX_EX)
          QY_AVG = 0.5*(QY_IN + QY_EX)

C.....Compute the edge integral
          DO K = 1,DOF
            CALL EDGE_INT_HYDRO(EL_IN, LED_IN, GED, GP_IN,
     &                           F_HAT,  G_HAT,  H_HAT)
            CALL EDGE_INT_HYDRO(EL_EX, LED_EX, GED, GP_EX,
     &                          -F_HAT, -G_HAT, -H_HAT)
          ENDDO

        ENDDO

1000  CONTINUE
      RETURN
      END SUBROUTINE


C***********************************************************************
C
C     SUBROUTINE EDGE_INT_LDG()
C
C     This subroutine computes the edge integrals for the LDG boundary
C     terms using Gauss quadrature and adds them to the LZ.
C
C     Written by Shintaro Bunya (01-02-2007)
C
C***********************************************************************

      SUBROUTINE EDGE_INT_LDG(K,EL_I,LED,GED,GP,QX,QY,NX,NY)
      ! <ezpp-noinst>
      
C.....Use appropriate modules

      USE GLOBAL,ONLY : AREAS
      USE DG,ONLY : M_INV,XLEN,PHI_EDGE,WEGP,LZ

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER K, EL_I, LED, GED, GP
      REAL(SZ) AREA, IMASS, QX, QY
      
C.....Retrieve the element area
      
      AREA = 0.5D0*AREAS(EL_I)

C.....Comput the edge integral

      IMASS = M_INV(K)/AREA

      LZ(K,1,1,EL_I) = LZ(K,1,1,EL_I)
     &          - IMASS*XLEN(GED)*QX*PHI_EDGE(K,GP,LED)*NX*WEGP(GP)
      LZ(K,1,2,EL_I) = LZ(K,1,2,EL_I)
     &          - IMASS*XLEN(GED)*QX*PHI_EDGE(K,GP,LED)*NY*WEGP(GP)
      LZ(K,2,1,EL_I) = LZ(K,2,1,EL_I)
     &          - IMASS*XLEN(GED)*QY*PHI_EDGE(K,GP,LED)*NX*WEGP(GP)
      LZ(K,2,2,EL_I) = LZ(K,2,2,EL_I)
     &          - IMASS*XLEN(GED)*QY*PHI_EDGE(K,GP,LED)*NY*WEGP(GP)
     
      RETURN
      END SUBROUTINE
