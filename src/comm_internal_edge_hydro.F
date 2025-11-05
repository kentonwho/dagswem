C***********************************************************************
C
C     SUBROUTINE INTERNAL_EDGE_HYDRO( )
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
C     Written by Ethan Kubatko (06-11-2004)
C
C***********************************************************************

      SUBROUTINE INTERNAL_EDGE_HYDRO(L)
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED
      
      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER L, LED_IN, LED_EX, GED, GP_IN, GP_EX
       
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
      
        ZE_IN = 0.D0
        QX_IN = 0.D0
        QY_IN = 0.D0
        HB_IN = 0.D0
        
        ZE_EX = 0.D0
        QX_EX = 0.D0
        QY_EX = 0.D0
        HB_EX = 0.D0

        DO K = 1,DOF

          ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN)
          QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN)
          QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN)
          HB_IN = HB_IN + HB(K,EL_IN,1  )*PHI_EDGE(K,GP_IN,LED_IN)
          
          ZE_EX = ZE_EX + ZE(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX)
          QX_EX = QX_EX + QX(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX)
          QY_EX = QY_EX + QY(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX)
          HB_EX = HB_EX + HB(K,EL_EX,1  )*PHI_EDGE(K,GP_EX,LED_EX)

        ENDDO
        
        IF (DOF.EQ.1) THEN
          DO K = 2,3
            HB_IN = HB_IN + HB(K,EL_IN,1  )*PHI_EDGE(K,GP_IN,LED_IN)
            HB_EX = HB_EX + HB(K,EL_EX,1  )*PHI_EDGE(K,GP_EX,LED_EX)
          ENDDO
        ENDIF
        
C.....Compute the Roe flux

        CALL NUMERICAL_FLUX()

C.....Compute the edge integral
       
        DO K = 1,DOF
          CALL EDGE_INT_HYDRO(EL_IN, LED_IN, GED, GP_IN,
     &                                            F_HAT,  G_HAT,  H_HAT)
          CALL EDGE_INT_HYDRO(EL_EX, LED_EX, GED, GP_EX,
     &                                           -F_HAT, -G_HAT, -H_HAT)
        ENDDO

      ENDDO
       
      RETURN
      END SUBROUTINE
      
      

