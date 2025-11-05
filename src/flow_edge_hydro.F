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

      SUBROUTINE FLOW_EDGE_HYDRO(lvl)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED, lvl
      REAL(SZ) TX, TY

      DO 1000 L = 1,NFEDS
      
C.....Retrieve the global and local edge number

        GED = NFEDN(L)
        LED = NEDSD(1,GED)

C.....Retrieve the element to which the edge belongs

        EL_IN = NEDEL(1,GED)
 
C.....If edge element is not assign this ts level, skip CJT
        IF (TsGroup(EL_IN) /= lvl) GOTO 1000
       
C.....If the element is dry then skip the edge calculation

        IF (WDFLG(EL_IN).EQ.0) GOTO 1000
        
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
CND
          SFAC_IN = SFACED(I,LED,EL_IN)
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
C
          SFAC_EX = SFAC_IN
        
C.....Compute the Roe flux

          CALL NUMERICAL_FLUX()


        
C.....Compute the edge integral

          DO K = 1,DOF
            CALL EDGE_INT_HYDRO(EL_IN, LED, GED, I, F_HAT, G_HAT, H_HAT)
          ENDDO
        
        ENDDO

1000  CONTINUE
      
      RETURN
      END SUBROUTINE
