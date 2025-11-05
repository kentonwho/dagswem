      SUBROUTINE EBARRIER_EDGE_HYDRO()
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

C.....Declare local variables

      INTEGER GED, LED, L
      REAL(SZ) ABOVE, SUBSUP, W_IN, TX, TY

C.....Loop over the external barrier segments

      DO 1000 L = 1,NEBSEG
      
C.......Obtain the global and local edge numbers

        GED = NEBSEGN(L)
        LED = NEDSD(1,GED)

C.......Obtain the element

        EL = NEDEL(1,GED)
        
        IF (WDFLG(EL).EQ.0) GOTO 1000

C.....Retrieve the components of the normal vector to the edge
        NX = COSNX(GED)
        NY = SINNX(GED)
        
C.....Set the components for the tangential vector to the edge

        TX = -NY
        TY =  NX
         
C.......Compute the variables at the quadrature points

        DO I = 1,NEGP
        
C.........Obtain the height of the barrier at the quadrature point

          ZE_IN = ZE(1,EL,IRK)
          QX_IN = QX(1,EL,IRK)
          QY_IN = QY(1,EL,IRK)
          HB_IN = HB(1,EL,1  )

          DO K = 2,DOF
            ZE_IN = ZE_IN + ZE(K,EL,IRK)*PHI_EDGE(K,I,LED)
            QX_IN = QX_IN + QX(K,EL,IRK)*PHI_EDGE(K,I,LED)
            QY_IN = QY_IN + QY(K,EL,IRK)*PHI_EDGE(K,I,LED)
            HB_IN = HB_IN + HB(K,EL,1  )*PHI_EDGE(K,I,LED)
          ENDDO
C
CND
          SFAC_IN = SFACED(I,LED,EL)

          ABOVE  = ZE_IN - EBHT(L)
          SUBSUP = 2.D0*ABOVE/3.D0
          
C.........Case 1:  Water is below barrier
C         ---------------------------------------------
          
          IF (ABOVE.LT.BARMIN) THEN
            Q_N_EXT = -(QX_IN*NX + QY_IN*NY)
            Q_T_EXT =   QX_IN*TX + QY_IN*TY
            
C.........Case 2:  Water is above barrier
C         ------------------------------------------------------------
            
          ELSE
            Q_N_EXT = RAMPDG*EBCFSP(L)*SUBSUP*SQRT(SUBSUP*G)
            Q_T_EXT = 0.D0
          ENDIF
          
C.....Set the exterior bed and surface elevation equal to the interior

          ZE_EX = ZE_IN
          HB_EX = HB_IN
          SFAC_EX = SFAC_IN
C.........Compute the x and y components of the external state flow

          QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)/(NX*TY - NY*TX)
          QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)/(NX*TY - NY*TX)
          
C.........Compute numerical flux
          
          CALL NUMERICAL_FLUX()

          DO K = 1,DOF

            W_IN = 2.0*M_INV(K)/AREAS(EL)*XLEN(GED)
     &                                   *PHI_EDGE(K,I,LED)*WEGP(I)

            RHS_ZE(K,EL,IRK) = RHS_ZE(K,EL,IRK) - W_IN*F_HAT
            RHS_QX(K,EL,IRK) = RHS_QX(K,EL,IRK) - W_IN*G_HAT
            RHS_QY(K,EL,IRK) = RHS_QY(K,EL,IRK) - W_IN*H_HAT

          ENDDO
        ENDDO
1000  CONTINUE

      RETURN
      END SUBROUTINE EBARRIER_EDGE_HYDRO
