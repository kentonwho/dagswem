C***********************************************************************
C
C     SUBROUTINE LDG_HYDRO(IT)
C
C     Compute variable Z to be used in the LDG terms.
C
C     Written by Shintaro Bunya (01-01-2007)
C
C***********************************************************************

      SUBROUTINE LDG_HYDRO(IT)

C.....Use appropriate modules
      
      USE SIZES,ONLY : SZ, myproc
      USE GLOBAL
      USE DG
      USE NodalAttributes, ONLY : EVM

c#ifdef CMPI
c      USE MESSENGER_ELEM,ONLY : UPDATELZ_ELEM
c#endif

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER IT,L,GED,NBOREL,NNBORS,NDRYNBORS
      INTEGER Detected
      REAL(SZ) QBCT1,QBCT2,QBC1,QBC2,ZP(3)
      REAL(SZ) DPAVG,EVMAvg

C.....Initialize

      DO L=1,NE
        DO K=1,DOF
          LZ(K,1,1,L) = 0.d0
          LZ(K,1,2,L) = 0.d0
          LZ(K,2,1,L) = 0.d0
          LZ(K,2,2,L) = 0.d0
        ENDDO
      ENDDO

C.....Compute elevation specified edges

      IF (NEEDS.GT.0)  CALL OCEAN_EDGE_LDG_HYDRO()

C.....Compute no-normal flow edges

      IF (NLEDS.GT.0)  CALL LAND_EDGE_LDG_HYDRO()

C.....Compute non-zero flow edges

      IF (NFEDS.GT.0)  CALL FLOW_EDGE_LDG_HYDRO()
        
C.....Compute radiation edges

      IF (NREDS.GT.0)  CALL RADIATION_EDGE_LDG_HYDRO()

C.....Compute internal edges
      CALL INTERNAL_EDGE_LDG_HYDRO()

C.....Loop over elements

      CALL RHS_LDG_HYDRO()

C.....Solve for Z
      DO L=1,NE
        N1 = NM(L,1)
        N2 = NM(L,2)
        N3 = NM(L,3)
        EVMAvg = (EVM(N1)+EVM(N2)+EVM(N3))/3.d0

        DO K=1,DOF
          LZ(K,1,1,L) = LZ(K,1,1,L)*EVMAvg
          LZ(K,1,2,L) = LZ(K,1,2,L)*EVMAvg
          LZ(K,2,1,L) = LZ(K,2,1,L)*EVMAvg
          LZ(K,2,2,L) = LZ(K,2,2,L)*EVMAvg
        ENDDO

      ENDDO
c
c#ifdef CMPI
c      CALL UPDATELZ_ELEM(LZ)
c#endif

      RETURN
      END SUBROUTINE


C***********************************************************************
C
C     SUBROUTINE OCEAN_EDGE_LDG_HYDRO( )
C
C     This subroutine does the following:
C
C       1.  Calculates the values of the necessary variables at the edge
C           gauss points for NO-NORMAL FLOW edges
C       2.  Compute the LDG flux at these points.
C       3.  Compute the boundary integrals.
C
C     Written by Shintaro Bunya (01-04-2007)
C
C***********************************************************************
      SUBROUTINE OCEAN_EDGE_LDG_HYDRO()

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED
      REAL(SZ) QX_AVG, QY_AVG

      DO 1000 L = 1, needs
      
C.....Retrieve the global and local edge number

        GED = NEEDN(L)
        LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge

        EL_IN = NEDEL(1,GED)
        
C.....Retrieve the components of the normal vector to the edge
        NX = COSNX(GED)
        NY = SINNX(GED)
        
C.....Retrieve the nodes of the edge
      
        N1 = NEDNO(1,GED)
        N2 = NEDNO(2,GED)
      
C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

        DO I = 1,NEGP

          QX_IN = 0.D0
          QY_IN = 0.D0

          QX_EX = 0.D0
          QY_EX = 0.D0
c
        
C.....Compute the solution at the interior state

          DO K = 1,DOF

            QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)
            QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)

          ENDDO
        
C.....Set the exterior state flows equal to the interior state flows

          QX_EX = QX_IN
          QY_EX = QY_IN

C.....Take the average

          QX_AVG = 0.5*(QX_IN + QX_EX)*SFACED(I,LED,EL_IN)
          QY_AVG = 0.5*(QY_IN + QY_EX)

C.....Compute the edge integral

          DO K = 1,DOF
            CALL EDGE_INT_LDG_HYDRO
     &           (K,EL_IN,LED,GED,I,QX_AVG,QY_AVG,NX,NY)
          ENDDO
        
        ENDDO

1000  CONTINUE
      
      RETURN
      END SUBROUTINE


C***********************************************************************
C
C     SUBROUTINE LAND_EDGE_HYDRO_LDG( )
C
C     This subroutine does the following:
C
C       1.  Calculates the values of the necessary variables at the edge
C           gauss points for NO-NORMAL FLOW edges
C       2.  Compute the LDG flux at these points.
C       3.  Compute the boundary integrals.
C
C     Written by Shintaro Bunya (01-02-2007)
C
C***********************************************************************

      SUBROUTINE LAND_EDGE_LDG_HYDRO()

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED
      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED
      REAL(SZ) TX, TY, QX_AVG, QY_AVG

      DO 1000 L = 1,NLEDS
      
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

        DO I = 1,NEGP

          QX_IN = 0.D0
          QY_IN = 0.D0
c

C.....Compute the solution at the interior state

          DO K = 1,DOF
          
            QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)
            QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)

          ENDDO
        
C.....Compute the velocity in the normal and tangental direction
        
          Q_N_INT = QX_IN*NX + QY_IN*NY
          Q_T_INT = QX_IN*TX + QY_IN*TY

C.....Reflect the velocity in the normal direction

          Q_N_EXT = -Q_N_INT
          Q_T_EXT =  Q_T_INT
        
C.....Compute the x and y components of the external state flow

          QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)/(NX*TY - NY*TX)
          QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)/(NX*TY - NY*TX)

C.....Take the average

          QX_AVG = 0.5*(QX_IN + QX_EX)*SFACED(I,LED,EL_IN)
          QY_AVG = 0.5*(QY_IN + QY_EX)

C.....Compute the edge integral
          DO K = 1,DOF
            CALL EDGE_INT_LDG_HYDRO
     &           (K,EL_IN,LED,GED,I,QX_AVG,QY_AVG,NX,NY)
         ENDDO

        ENDDO

1000  CONTINUE
      
      RETURN
      END SUBROUTINE


C***********************************************************************
C
C     SUBROUTINE FLOW_EDGE_LDG_HYDRO( )
C
C     This subroutine does the following:
C
C       1.  Calculates the values of the necessary variables at the edge
C           gauss points for INTERNAL edges
C       2.  Take the average of these values
C       3.  Perform boundary integration
C
C     Written by Shintaro Bunya (01-05-2007)
C
C***********************************************************************

      SUBROUTINE FLOW_EDGE_LDG_HYDRO()

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED
      REAL(SZ) TX, TY, QX_AVG, QY_AVG

      DO 1000 L = 1,NFEDS
      
C.....Retrieve the global and local edge number

        GED = NFEDN(L)
        LED = NEDSD(1,GED)

C.....Retrieve the element to which the edge belongs

        EL_IN = NEDEL(1,GED)
        
C.....Retrieve the components of the normal vector to the edge
        NX = COSNX(GED)
        NY = SINNX(GED)
      
C.....Set the components for the tangential vector to the edge

        TX = -NY
        TY =  NX
      
C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

        DO I = 1,NEGP

          QX_IN = 0.D0
          QY_IN = 0.D0

          QX_EX = 0.D0
          QY_EX = 0.D0
        
C.....Compute the specified flow boundaries for the exterior state

          Q_N_EXT = 0.D0
          DO JJ = 1,NFFR
        
            IF(FPER(JJ).EQ.0.D0) THEN
              NCYC = 0.D0
            ELSE
              NCYC = INT(TIMEDG/FPER(JJ))
            ENDIF

            ARGJ = FAMIG(JJ)*(TIMEDG - NCYC*FPER(JJ)) + FFACE(JJ)
            RFF  = FFF(JJ)*RAMPDG
          
            QNAM_GP = 0.5D0*(QNAM_DG(JJ,L,1) + QNAM_DG(JJ,L,2))
     &              + 0.5D0*(QNAM_DG(JJ,L,2) - QNAM_DG(JJ,L,1))*XEGP(I)
            QNPH_GP = 0.5D0*(QNPH_DG(JJ,L,1) + QNPH_DG(JJ,L,2))
     &              + 0.5D0*(QNPH_DG(JJ,L,2) - QNPH_DG(JJ,L,1))*XEGP(I)
     
            ARG = ARGJ - QNPH_GP
          
            Q_N_EXT = Q_N_EXT + QNAM_GP*RFF*COS(ARG)
            Q_T_EXT =  0.D0

            QX_EX = -( TY*Q_N_EXT - NY*Q_T_EXT)/(NX*TY - NY*TX)
            QY_EX = -(-TX*Q_N_EXT + NX*Q_T_EXT)/(NX*TY - NY*TX)

          ENDDO
        
C.....Compute the solution at the interior state

          DO K = 1,DOF

            QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)
            QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)

          ENDDO
        
C.....Take the average

          QX_AVG = 0.5*(QX_IN + QX_EX)*SFACED(I,LED,EL_IN)
          QY_AVG = 0.5*(QY_IN + QY_EX)

C.....Compute the edge integral
          DO K = 1,DOF
            CALL EDGE_INT_LDG_HYDRO
     &           (K,EL_IN,LED,GED,I,QX_AVG,QY_AVG,NX,NY)
         ENDDO

        ENDDO

1000  CONTINUE
      
      RETURN
      END SUBROUTINE


C***********************************************************************
C
C     SUBROUTINE LAND_EDGE_LDG_HYDRO( )
C
C     This subroutine does the following:
C
C       1.  Calculates the values of the necessary variables at the edge
C           gauss points for INTERNAL edges
C       2.  Take the average of these values
C       3.  Perform boundary integration
C
C     Written by Shintaro Bunya (01-10-2007)
C
C***********************************************************************

      SUBROUTINE RADIATION_EDGE_LDG_HYDRO()

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED
      REAL(SZ) TX, TY

      DO 1000 L = 1,NREDS
      
C.....Retrieve the global and local edge number

        GED = NREDN(L)
        LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge

        EL_IN = NEDEL(1,GED)

C.....Retrieve the components of the normal vector to the edge
        NX = COSNX(GED)
        NY = SINNX(GED)

C.....Compute ZE, QX, QY, and HB at each Gauss point

        DO I = 1,NEGP

          QX_IN = 0.D0
          QY_IN = 0.D0

C.....Compute the solution at the interior state

          DO K = 1,DOF
          
            QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)
            QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)

          ENDDO
CND
          QX_IN = QX_IN*SFACED(I,LED,EL_IN)
C.....Compute the edge integral
          DO K = 1,DOF
            CALL EDGE_INT_LDG_HYDRO
     &           (K,EL_IN,LED,GED,I,QX_IN,QY_IN,NX,NY)
         ENDDO

        ENDDO

1000  CONTINUE      
      RETURN
      END SUBROUTINE
C***********************************************************************
C
C     SUBROUTINE INTERNAL_EDGE_LDG_HYDRO( )
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

      SUBROUTINE INTERNAL_EDGE_LDG_HYDRO()

C.....Use appropriate modules
      
      USE SIZES,ONLY : SZ
      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED_IN, LED_EX, GED, GP_IN, GP_EX
      REAL(SZ) QX_AVG, QY_AVG, W_IN, W_EX

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

          QX_AVG = 0.5*(QX_IN*SFACED(I,LED_IN,EL_IN)
     $         + QX_EX*SFACED(I,LED_EX,EL_EX))
          QY_AVG = 0.5*(QY_IN + QY_EX)

C.....Compute the edge integral
          DO K = 1,DOF
#if 0
            CALL EDGE_INT_LDG_HYDRO
     &           (K,EL_IN,LED_IN,GED,GP_IN,QX_AVG,QY_AVG, NX, NY)
            CALL EDGE_INT_LDG_HYDRO
     &           (K,EL_EX,LED_EX,GED,GP_EX,QX_AVG,QY_AVG,-NX,-NY)
#else
            W_IN = 2.0*M_INV(K)/AREAS(EL_IN)*XLEN(GED)*
     &             PHI_EDGE(K,GP_IN,LED_IN)*WEGP(GP_IN)
            W_EX = 2.0*M_INV(K)/AREAS(EL_EX)*XLEN(GED)*
     &             PHI_EDGE(K,GP_EX,LED_EX)*WEGP(GP_EX)

            LZ(K,1,1,EL_IN) = LZ(K,1,1,EL_IN) - QX_AVG*NX*W_IN
            LZ(K,1,2,EL_IN) = LZ(K,1,2,EL_IN) - QX_AVG*NY*W_IN
            LZ(K,2,1,EL_IN) = LZ(K,2,1,EL_IN) - QY_AVG*NX*W_IN
            LZ(K,2,2,EL_IN) = LZ(K,2,2,EL_IN) - QY_AVG*NY*W_IN

            LZ(K,1,1,EL_EX) = LZ(K,1,1,EL_EX) + QX_AVG*NX*W_EX
            LZ(K,1,2,EL_EX) = LZ(K,1,2,EL_EX) + QX_AVG*NY*W_EX
            LZ(K,2,1,EL_EX) = LZ(K,2,1,EL_EX) + QY_AVG*NX*W_EX
            LZ(K,2,2,EL_EX) = LZ(K,2,2,EL_EX) + QY_AVG*NY*W_EX
#endif
          ENDDO

        ENDDO

1000  CONTINUE
      RETURN
      END SUBROUTINE


C***********************************************************************
C
C     SUBROUTINE EDGE_INT_LDG_HYDRO()
C
C     This subroutine computes the edge integrals for the LDG boundary
C     terms using Gauss quadrature and adds them to the LZ.
C
C     Written by Shintaro Bunya (01-02-2007)
C
C***********************************************************************

      SUBROUTINE EDGE_INT_LDG_HYDRO(K,EL,LED,GED,GP,QX,QY,NX,NY)
      ! <ezpp-noinst>
      
C.....Use appropriate modules

      USE SIZES,ONLY : SZ
      USE GLOBAL,ONLY : AREAS
      USE DG,ONLY : M_INV,XLEN,PHI_EDGE,WEGP,LZ

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER K, EL, LED, GED, GP
      REAL(SZ) AREA, IMASS, QX, QY, NX, NY
      
C.....Retrieve the element area
      
      AREA = 0.5D0*AREAS(EL)

C.....Comput the edge integral

      IMASS = M_INV(K)/(0.5D0*AREA)

      LZ(K,1,1,EL) = LZ(K,1,1,EL)
     &          - IMASS*XLEN(GED)*0.5*QX*PHI_EDGE(K,GP,LED)*NX*WEGP(GP)
      LZ(K,1,2,EL) = LZ(K,1,2,EL)
     &          - IMASS*XLEN(GED)*0.5*QX*PHI_EDGE(K,GP,LED)*NY*WEGP(GP)
      LZ(K,2,1,EL) = LZ(K,2,1,EL)
     &          - IMASS*XLEN(GED)*0.5*QY*PHI_EDGE(K,GP,LED)*NX*WEGP(GP)
      LZ(K,2,2,EL) = LZ(K,2,2,EL)
     &          - IMASS*XLEN(GED)*0.5*QY*PHI_EDGE(K,GP,LED)*NY*WEGP(GP)
     
      RETURN
      END SUBROUTINE


C***********************************************************************
C
C     SUBROUTINE RHS_LDG_HYDRO()
C
C     This subroutine computes the area integrals for the LDG hydro and
C     adds them into the LZ.
C
C     Written by Shintaro Bunya (01-02-2007)
C
C***********************************************************************

      SUBROUTINE RHS_LDG_HYDRO()
      
C.....Use appropriate modules
      USE SIZES,ONLY : SZ, myproc
      USE GLOBAL,ONLY : NE,NM,N1,N2,N3,sfacdub
      USE DG,ONLY : QX_IN,QY_IN,LZ,M_INV,
     &     DRPHI,DSPHI,DRDX,DRDY,DSDX,DSDY,
     &     WAGP,NAGP,DOF,QX,QY,PHI_AREA,IRK,SFAC_ELEM

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER L,K,I
      real(sz) SFACQUAD

      DO L=1,NE
      
C.....Retrieve the global node numbers for the element

      N1 = NM(L,1)
      N2 = NM(L,2)
      N3 = NM(L,3)
      
C.....Compute ZE, QX, QY, and HB at each area Gauss quadrature point

      DO I = 1,NAGP
      
        QX_IN = 0.D0
        QY_IN = 0.D0
        
        DO K = 1,DOF

          QX_IN = QX_IN + QX(K,L,IRK)*PHI_AREA(K,I)
          QY_IN = QY_IN + QY(K,L,IRK)*PHI_AREA(K,I)

        ENDDO
CND
        SFACQUAD = SFAC_ELEM(I,L)

C.....Build the rhs
        DO K = 1,DOF
          LZ(K,1,1,L) = LZ(K,1,1,L)
     &       + M_INV(K)*QX_IN*SFACQUAD*                  ! <--- QX/Mk
     &         ( DRPHI(K,I)*DRDX(L) + DSPHI(K,I)*DSDX(L) )*  ! <--- dphi/dx
     &         WAGP(I)                                       ! <--- weight

          LZ(K,1,2,L) = LZ(K,1,2,L)
     &       + M_INV(K)*QX_IN*                               ! <--- QX/Mk
     &         ( DRPHI(K,I)*DRDY(L) + DSPHI(K,I)*DSDY(L) )*  ! <--- dphi/dy
     &         WAGP(I)                                       ! <--- weight

          LZ(K,2,1,L) = LZ(K,2,1,L)
     &       + M_INV(K)*QY_IN*SFACQUAD*                  ! <--- QY/Mk
     &         ( DRPHI(K,I)*DRDX(L) + DSPHI(K,I)*DSDX(L) )*  ! <--- dphi/dx
     &         WAGP(I)                                       ! <--- weight

          LZ(K,2,2,L) = LZ(K,2,2,L)
     &       + M_INV(K)*QY_IN*                               ! <--- QY/Mk
     &         ( DRPHI(K,I)*DRDY(L) + DSPHI(K,I)*DSDY(L) )*  ! <--- dphi/dy
     &         WAGP(I)                                       ! <--- weight
        ENDDO

      ENDDO

      ENDDO
      RETURN
      END SUBROUTINE


