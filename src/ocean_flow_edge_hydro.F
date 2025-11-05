C***********************************************************************
C
C     SUBROUTINE FLOW_EDGE_HYDRO( )
C
C     This subroutine does the following:
C
C       1.  Calculates the values of the necessary variables at the edge
C           gauss points for NON-ZERO FLUX edges
C       2.  Calls the appropriate subroutine to compute the flux at
C           these points.
C       3.  Calls the appropriate subroutine to compute the boundary
C           integrals.
C
C     Written by Ethan Kubatko (06-11-2004)
C
C***********************************************************************

      SUBROUTINE OCEAN_FLOW_EDGE_HYDRO(L,QBC)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED
      REAL(SZ) TX, TY, QBC
      
C.....Retrieve the global and local edge number

      GED = NEEDN(L)
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

        ZE_IN = 0.D0
        QX_IN = 0.D0
        QY_IN = 0.D0

        ZE_EX = 0.D0
        QX_EX = 0.D0
        QY_EX = 0.D0
        
        HB_IN = 0.D0
        
C.....Compute the specified flow boundaries for the exterior state

        Q_N_EXT = QBC
        Q_T_EXT = 0.D0

        QX_EX = -( TY*Q_N_EXT - NY*Q_T_EXT)/(NX*TY - NY*TX)
        QY_EX = -(-TX*Q_N_EXT + NX*Q_T_EXT)/(NX*TY - NY*TX)

C.....Compute the solution at the interior state

        DO K = 1,DOF

          ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)
          QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)
          QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)
          
          HB_IN = HB_IN + HB(K,EL_IN,1  )*PHI_EDGE(K,I,LED)

        ENDDO
        
        IF (DOF.EQ.1) THEN
          DO K = 2,3
            HB_IN = HB_IN + HB(K,EL_IN,1  )*PHI_EDGE(K,I,LED)
          ENDDO
        ENDIF
        
C.....Set the exterior bed and surface elevation equal to the interior

        ZE_EX = ZE_IN
        HB_EX = HB_IN
        
C.....Compute the Roe flux

        CALL NUMERICAL_FLUX()
        
C.....Compute the edge integral

        DO K = 1,DOF
          CALL EDGE_INT_HYDRO(EL_IN, LED, GED, I, F_HAT, G_HAT, H_HAT)
        ENDDO
        
      ENDDO
      
      RETURN
      END SUBROUTINE
