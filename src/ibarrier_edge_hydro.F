C***********************************************************************
C
C     SUBROUTINE IBARRIER_EDGE_HYDRO()
C
C     This subroutine does the following:
C
C       1.  Calculates the values of the necessary variables at the edge
C           gauss points for INTERNAL BARRIER edges
C       2.  Calls the appropriate subroutine to compute the flux at
C           these points.
C       3.  Computes the boundary integrals.
C
C     Written by Ethan Kubatko (08-11-2008)
C
C***********************************************************************

      SUBROUTINE IBARRIER_EDGE_HYDRO(IT,LVL)
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED
      USE SIZES, ONLY:  MYPROC

      IMPLICIT NONE

C.....Declare local variables

      INTEGER GEDB, GEDF, LEDB, LEDF, ELB, ELF, GPB, GPF, L, WEIR_FLOW
      INTEGER IT, NB1, NB2, NF1, NF2, LVL
      REAL(SZ) ZEB, QXB, QYB, HBB, ZEF, QXF, QYF, HBF
      REAL(SZ) SFACB,SFACF
      REAL(SZ) FB_HAT, GB_HAT, HB_HAT, FF_HAT, GF_HAT, HF_HAT
      REAL(SZ) ABOVEB, ABOVEF
      REAL(SZ) QNIB(2)
      REAL(SZ) NLEQG_TMP, G_TMP
      REAL(SZ) SUBSUPB, SUBSUPF, WEGPB, WEGPF
      REAL(SZ) NXB, NYB, NXF, NYF, TXB, TYB, TXF, TYF
      REAL(SZ) QB_N_INT, QF_N_INT, QB_N_EXT, QF_N_EXT
      REAL(SZ) QB_T_INT, QF_T_INT, QB_T_EXT, QF_T_EXT

C.....Loop over the internal barrier segments (note: an internal barrier
C.....segment consists of two internal barrier edges -- a "front" edge &
C.....a "back" edge.
      DO 1000 L = 1,NIBSEG
      
C.......Obtain the global and local edges of the back and front sides

        GEDB = NIBSEGN(1,L)
        GEDF = NIBSEGN(2,L)
        if (gedf.eq.0) go to 1000
        
        LEDB = NEDSD(1,GEDB)
        LEDF = NEDSD(1,GEDF)
         
C.......Obtain the elements of the back and front sides

        ELB = NEDEL(1,GEDB)
        ELF = NEDEL(1,GEDF)



                
        NB1 = NEDNO(1,GEDB)
        NB2 = NEDNO(2,GEDB)
        
        NF1 = NEDNO(1,GEDF)
        NF2 = NEDNO(2,GEDF)
        if (tsgroup(elb).ne.tsgroup(elf) ) then
           write(*,*) 'error in ibarrier ',myproc,elf,elb,nb1,nb2,nf1,nf2
           stop
        endif
        if (tsgroup(elb).ne.lvl) go to 1000
        IF((WDFLG(ELB).EQ.0).AND.(WDFLG(ELF).EQ.0)) GOTO 1000

C.......Retrieve the components of the normal vector to the edge

        NXB = COSNX(GEDB)
        NYB = SINNX(GEDB)
        
        NXF = COSNX(GEDF)
        NYF = SINNX(GEDF)
        
C.....Set the components for the tangential vector to the edge

        TXB = -NYB
        TYB =  NXB

        TXF = -NYF
        TYF =  NXF
         
C.......Compute the variables at the quadrature points

        DO I = 1,NEGP
        
          GPB = I
          GPF = NEGP - I + 1
        
C.........Obtain the height of the barrier at the quadrature point

          ZEB = ZE(1,ELB,IRK)
          QXB = QX(1,ELB,IRK)
          QYB = QY(1,ELB,IRK)
          HBB = HB(1,ELB,1  )
          
          ZEF = ZE(1,ELF,IRK)
          QXF = QX(1,ELF,IRK)
          QYF = QY(1,ELF,IRK)
          HBF = HB(1,ELF,1  )

          DO K = 2,DOF
            ZEB = ZEB + ZE(K,ELB,IRK)*PHI_EDGE(K,GPB,LEDB)
            QXB = QXB + QX(K,ELB,IRK)*PHI_EDGE(K,GPB,LEDB)
            QYB = QYB + QY(K,ELB,IRK)*PHI_EDGE(K,GPB,LEDB)
            HBB = HBB + HB(K,ELB,1  )*PHI_EDGE(K,GPB,LEDB)

            ZEF = ZEF + ZE(K,ELF,IRK)*PHI_EDGE(K,GPF,LEDF)
            QXF = QXF + QX(K,ELF,IRK)*PHI_EDGE(K,GPF,LEDF)
            QYF = QYF + QY(K,ELF,IRK)*PHI_EDGE(K,GPF,LEDF)
            HBF = HBF + HB(K,ELF,1  )*PHI_EDGE(K,GPF,LEDF)
          ENDDO
CND
          SFACB = SFACED(I,LEDB,ELB)
          SFACF = SFACED(I,LEDF,ELF)
C
          ABOVEB = 0
          ABOVEF = 0 
          IF (WDFLG(ELB).EQ.1) ABOVEB = ZEB - IBHT(L)
          IF (WDFLG(ELF).EQ.1) ABOVEF = ZEF - IBHT(L)

C.........Case 1:  Water is below barrier on both sides
C         ---------------------------------------------
          
          IF ((ABOVEF.LE.BARMIN).AND.(ABOVEB.LE.BARMIN)) THEN
            QB_N_INT = -(QXB*NXB + QYB*NYB)
            QF_N_INT = -(QXF*NXF + QYF*NYF)
            QB_T_INT = QXB*TXB + QYB*TYB
            QF_T_INT = QXF*TXF + QYF*TYF
            WEIR_FLOW = 0
            GOTO 100
            
C.........Case 2:  Water is above on both sides and equal (within tol)
C         ------------------------------------------------------------
            
          ELSEIF (ABS(ABOVEF-ABOVEB).LT.BARMIN) THEN
            QB_N_INT = -(QXB*NXB + QYB*NYB)
            QF_N_INT = -(QXF*NXF + QYF*NYF)
            QB_T_INT = QXB*TXB + QYB*TYB
            QF_T_INT = QXF*TXF + QYF*TYF
            WEIR_FLOW = 0
            GOTO 100
          ENDIF

          SUBSUPF = 2.D0*ABOVEF/3.D0
          SUBSUPB = 2.D0*ABOVEB/3.D0

C.........Case 3: Overtopping of barrier with water higher on front side
C                 Flow from front to back side
C         -------------------------------------------------------------

          IF ((ABOVEF.GT.ABOVEB).AND.(ABOVEF.GT.BARMIN)) THEN
          
            WEIR_FLOW = 1
          
C...........Case 3a) Subcritical flow
C           -------------------------

            IF (ABOVEB.GT.SUBSUPF) THEN
              QF_N_INT = RAMPDG*IBCFSB(L)*ABOVEB
     &                                   *SQRT((2.D0*G*(ABOVEF-ABOVEB)))
              QF_T_INT = 0.D0
              
C...........Case 3b) Supercritical flow
C           ---------------------------

            ELSE
              QF_N_INT = RAMPDG*IBCFSP(L)*SUBSUPF*SQRT(SUBSUPF*G)
              QF_T_INT = 0.D0
            ENDIF
            GOTO 100

          ENDIF
          
C.........Case 4: Overtopping of barrier with water higher on back side
C                 Flow from back to front side
C         --------------------------------------------------------------

          IF ((ABOVEB.GT.ABOVEF).AND.(ABOVEB.GT.BARMIN)) THEN
          
            WEIR_FLOW = -1
            
C...........Case 4a) Subcritical flow
C           -------------------------

            IF (ABOVEF.GT.SUBSUPB) THEN
              QB_N_INT = RAMPDG*IBCFSB(L)*ABOVEF
     &                                   *SQRT((2.D0*G*(ABOVEB-ABOVEF)))
              QB_T_INT = 0.D0
              
C...........Case 4b) Supercritical flow
C           ---------------------------
            
            ELSE
              QB_N_INT = RAMPDG*IBCFSP(L)*SUBSUPB*SQRT(SUBSUPB*G)
              QB_T_INT = 0.D0
            ENDIF
            GOTO 100
          ENDIF
          
100       CONTINUE

          IF (WEIR_FLOW.LE.0) THEN

C...........Compute the numerical flux for the back side edge

            ZE_IN = ZEB
            QX_IN = QXB
            QY_IN = QYB
            HB_IN = HBB
          
            ZE_EX = ZEB
            HB_EX = HBB
CND
            SFAC_IN = SFACB
            SFAC_EX = SFACB
C          
            NX = NXB
            NY = NYB

C...........Reflect the velocity in the normal direction

            Q_N_EXT = QB_N_INT
            Q_T_EXT = QB_T_INT

C...........Compute the x and y components of the external state flow

            QX_EX = ( TYB*Q_N_EXT - NYB*Q_T_EXT)/(NXB*TYB - NYB*TXB)
            QY_EX = (-TXB*Q_N_EXT + NXB*Q_T_EXT)/(NXB*TYB - NYB*TXB)

C...........Compute the numerical flux
          
            CALL NUMERICAL_FLUX()
            FB_HAT = F_HAT
            GB_HAT = G_HAT
            HB_HAT = H_HAT
            IF (WEIR_FLOW.LT.0) THEN
              ZE_IN = ZEF
              QX_IN = QXF
              QY_IN = QYF
              HB_IN = HBF
              ZE_EX = ZEF
              HB_EX = HBF
CND
              SFAC_IN = SFACF
              SFAC_EX = SFACF
C
              NX = NXF
              NY = NYF
              IF (WDFLG(ELF).EQ.0) THEN
                NLEQG_TMP = NLEQG
                NLEQG = 0.D0
                G_TMP = G
                G = 0.D0
              ENDIF
              CALL NUMERICAL_FLUX()
              FF_HAT = F_HAT
              GF_HAT = G_HAT
              HF_HAT = H_HAT
              IF (WDFLG(ELF).EQ.0) THEN
                NLEQG = NLEQG_TMP
                G = G_TMP
              ENDIF
              GOTO 200
            ENDIF
         ENDIF

         IF (WEIR_FLOW.GE.0) THEN
          
C...........Compute the numerical flux for the front side edge

            ZE_IN = ZEF
            QX_IN = QXF
            QY_IN = QYF
            HB_IN = HBF

            ZE_EX = ZEF
            HB_EX = HBF
CND
            SFAC_IN = SFACF
            SFAC_EX = SFACF
C          
            NX = NXF
            NY = NYF
          
C...........Reflect the velocity in the normal direction

            Q_N_EXT = QF_N_INT
            Q_T_EXT = QF_T_INT

C...........Compute the x and y components of the external state flow

            QX_EX = ( TYF*Q_N_EXT - NYF*Q_T_EXT)/(NXF*TYF - NYF*TXF)
            QY_EX = (-TXF*Q_N_EXT + NXF*Q_T_EXT)/(NXF*TYF - NYF*TXF)

            CALL NUMERICAL_FLUX()
            FF_HAT = F_HAT
            GF_HAT = G_HAT
            HF_HAT = H_HAT
            
            IF (WEIR_FLOW.GT.0) THEN
              ZE_IN = ZEB
              QX_IN = QXB
              QY_IN = QYB
              HB_IN = HBB
              ZE_EX = ZEB
              HB_EX = HBB
CND
              SFAC_IN = SFACB
              SFAC_EX = SFACB
C
              NX = NXB
              NY = NYB
              IF (WDFLG(ELB).EQ.0) THEN
                NLEQG_TMP = NLEQG
                NLEQG = 0.D0
                G_TMP = G
                G = 0.D0
              ENDIF
              CALL NUMERICAL_FLUX()
              FB_HAT = F_HAT
              GB_HAT = G_HAT
              IF (WDFLG(ELB).EQ.0) THEN
                NLEQG = NLEQG_TMP
                G = G_TMP
              ENDIF
              HB_HAT = H_HAT
              GOTO 200
            ENDIF
          ENDIF
          
200       CONTINUE
c
          DO K = 1,DOF

            WEGPB = 2.0*M_INV(K)/AREAS(ELB)*XLEN(GEDB)
     &                                   *PHI_EDGE(K,GPB,LEDB)*WEGP(GPB)
            WEGPF = 2.0*M_INV(K)/AREAS(ELF)*XLEN(GEDF)
     &                                   *PHI_EDGE(K,GPF,LEDF)*WEGP(GPF)

            RHS_ZE(K,ELB,IRK) = RHS_ZE(K,ELB,IRK) - WEGPB*FB_HAT
            RHS_QX(K,ELB,IRK) = RHS_QX(K,ELB,IRK) - WEGPB*GB_HAT
            RHS_QY(K,ELB,IRK) = RHS_QY(K,ELB,IRK) - WEGPB*HB_HAT

            RHS_ZE(K,ELF,IRK) = RHS_ZE(K,ELF,IRK) - WEGPF*FF_HAT
            RHS_QX(K,ELF,IRK) = RHS_QX(K,ELF,IRK) - WEGPF*GF_HAT
            RHS_QY(K,ELF,IRK) = RHS_QY(K,ELF,IRK) - WEGPF*HF_HAT

          ENDDO
        ENDDO
1000  CONTINUE
      RETURN
      END SUBROUTINE IBARRIER_EDGE_HYDRO
