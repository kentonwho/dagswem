C***********************************************************************
C
C     SUBROUTINE WETDRY()
C
C     This subroutine performs wetting and drying -- for details see:
C
C     "A Wetting and Drying Treatment for the Runge-Kutta Discontinuous
C     Galerkin Solution to the Shallow Water Equations", Shintaro Bunya
C     et al, Computer Methods in Applied Mechanics and Eng.  
C
C     Written by Shintaro Bunya and Ethan Kubatko
C
C***********************************************************************

      SUBROUTINE WETDRY(lvl)
      
C.....Use appropriate modules
      
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : NODECODE, X, Y, NM, DP, NE, NP, H0, AREAS
      USE DG
      
      IMPLICIT NONE

C.....Declare local variables

      REAL(SZ), PARAMETER :: ZERO = 1.D-6
      REAL(SZ), PARAMETER :: HUGE = 1.D+6

      INTEGER I, J, K1, K2, K, NM1, NM2, NM3, ISTOP,lvl,ii
      INTEGER EL1, EL2, EL3, DRYALL
      INTEGER ELEMENT_CHECK, HT_MAX_I, HT_MIN_I, NDRYNODE, TMP
      INTEGER L2H(3)
      
      REAL(SZ) AREA
      REAL(SZ) DEN, DP_AVG, DP_MIN
      REAL(SZ) HT_AVG, HT_MAX, HT_MIN, HT_VOL, H0C
      REAL(SZ) NUM1, NUM2
      REAL(SZ) QX_SUM, QY_SUM
      REAL(SZ) ZE_VOL

      REAL(SZ) DH_HAT(3)
      REAL(SZ) HT_HAT(NCHECK), HT_MODE(3)
      REAL(SZ) HT_NODE(NCHECK)
      REAL(SZ) QX_HAT(3), QX_NODE(3), QY_HAT(3), QY_NODE(3)
      REAL(SZ) ZE_HAT(NCHECK), ZE_NODE(NCHECK)
      
C.....Initialize node codes to dry

c      DO I = 1,NP
c        NODECODE(I) = 0
c      ENDDO
      
C.....Loop over the elements

      DO 100 ii = 1,nelem(lvl)
         j=elemnum(lvl,ii)
         wdflg_tmp(j)=wdflg(j)
      
C.......Retrieve node numbers for element

        NM1 = NM(J,1)
        NM2 = NM(J,2)
        NM3 = NM(J,3)
        nodecode(nm1)=0
        nodecode(nm2)=0
        nodecode(nm3)=0
        
C.......Compute value of ZE and HT at each node, find HT max and min,
C.......and count number of dry nodes

        NDRYNODE = 0
        HT_MIN = HUGE
        HT_MAX = ZERO
        DO I = 1,NCHECK
          ZE_NODE(I) = ZE(1,J,IRK+1)
          DO K = 2,DOF
            ZE_NODE(I) = ZE_NODE(I) + PHI_CHECK(K,I)*ZE(K,J,IRK+1)
          ENDDO
          HT_NODE(I) = ZE_NODE(I) + DP_NODE(I,J)
          IF (HT_NODE(I).LT.HT_MIN) THEN
            HT_MIN_I = I
            HT_MIN = HT_NODE(I)
          ENDIF
          IF (HT_NODE(I).GT.HT_MAX) THEN
            HT_MAX_I = I
            HT_MAX = HT_NODE(I)
          ENDIF
          IF (HT_NODE(I).LE.(H0+ZERO)) NDRYNODE = NDRYNODE + 1
        ENDDO
        
C.......Set "average" bathymetry based on p

        IF (P.EQ.0) THEN
          DP_AVG = DPE_MIN(J)
        ELSE
          DP_AVG = C13*(DP_NODE(1,J) + DP_NODE(2,J) + DP_NODE(3,J))
        ENDIF

C-----------------------------------------------------------------------
C.......Case 1: All nodes are wet
C-----------------------------------------------------------------------

        IF (NDRYNODE.EQ.0) THEN

C.........Case 1A:  Element was previously wet
C----------------------------------------------

          IF (WDFLG(J).EQ.1) THEN
        
C...........Make no adjustments and set element node codes to wet

            ELEMENT_CHECK = 0
            NODECODE(NM1) = 1
            NODECODE(NM2) = 1
            NODECODE(NM3) = 1
            
C.........Case 1B:  Element was previously dry
C----------------------------------------------
            
          ELSE
          
C...........Make no adjustments now

            ZE_HAT = ZE_NODE

C...........But set flag to do element check

            ELEMENT_CHECK = 1
          ENDIF
          
C-----------------------------------------------------------------------
C.......Case 2: All nodes are dry
C-----------------------------------------------------------------------

        ELSEIF (NDRYNODE.EQ.NCHECK) THEN
        
C.........Set surface elevation parallel to bathymetry such that total
C..........water depth at all points = average water depth

          ZE(:,J,IRK+1) = 0.D0
          ZE(1,J,IRK+1) = H0 - DP_AVG
          ZE(2,J,IRK+1) = -HB(2,J,1)
          ZE(3,J,IRK+1) = -HB(3,J,1)

C.........Zero out the fluxes

          QX(:,J,IRK+1) = 0.D0
          QY(:,J,IRK+1) = 0.D0
            
C.........Set flags to dry and skip element check

          WDFLG(J) = 0
          ELEMENT_CHECK = 0
          
C-----------------------------------------------------------------------
C.......Case 4: At least one node is dry and at least one node is wet
C-----------------------------------------------------------------------

        ELSE

C.........Set averages based on p

          IF (P.EQ.0) THEN
            HT_AVG = H0
          ELSEIF (P.EQ.1) THEN
            HT_AVG = ZE(1,J,IRK+1) + DP_AVG
          ELSE
            ZE_VOL = 0.D0
            DO K = 1,DOF
              ZE_VOL = ZE_VOL + ZE(K,J,IRK+1)*PHI_INTEGRATED(K)
            ENDDO
            HT_VOL = DP_VOL(J) + 0.25D0*AREAS(J)*ZE_VOL
            HT_AVG = HT_VOL/(0.5D0*AREAS(J))
          ENDIF
          
C         --------------------------------------------------------------
C.........Case 3A: Average water depth < H0
C         --------------------------------------------------------------

          IF (HT_AVG.LE.(H0+ZERO)) THEN

C...........Set surface elevation parallel to bathymetry such that total
C...........water depth at all points = average water depth

            ZE(:,J,IRK+1) = 0.D0
            ZE(1,J,IRK+1) = HT_AVG - DP_AVG
            ZE(2,J,IRK+1) = -HB(2,J,1)
            ZE(3,J,IRK+1) = -HB(3,J,1)
            ZE_HAT = ZE(1,J,IRK+1)

C...........Zero out the fluxes

            QX(:,J,IRK+1) = 0.D0
            QY(:,J,IRK+1) = 0.D0

C         --------------------------------------------------------------
C.........Case 3B: Aaverage water depth > H0
C         --------------------------------------------------------------
        
          ELSE

C...........If p>1 then find a linear function that has the same average
C...........water depth and the same linear components

            IF (P.GT.1) THEN
              HT_MODE(1) = HT_AVG
              HT_MODE(2) = -(ZE(2,J,1) + HB(2,J,1))
              HT_MODE(3) = -(ZE(3,J,1) + HB(3,J,1))
              HT_NODE = 0.D0
              DO K = 1,3
                HT_NODE(1) = HT_NODE(1) + PHI_CORNER(K,1)*HT_MODE(K)
                HT_NODE(2) = HT_NODE(2) + PHI_CORNER(K,2)*HT_MODE(K)
                HT_NODE(3) = HT_NODE(3) + PHI_CORNER(K,3)*HT_MODE(K)
              ENDDO
            ENDIF

C...........Reorder the nodes from shallowest to deepest

            L2H(1) = 1
            L2H(2) = 2
            L2H(3) = 3
            DO K1 = 1,3
              DO K2 = K1,3
                IF ( HT_NODE(L2H(K1)).GT.HT_NODE(L2H(K2)) ) THEN
                  TMP = L2H(K1)
                  L2H(K1) = L2H(K2)
                  L2H(K2) = TMP
                ENDIF
              ENDDO
            ENDDO

C...........Compute the factors needed for nodal ajustment

            DH_HAT(L2H(1)) = MAX(H0, HT_NODE(L2H(1))) - HT_NODE(L2H(1))
            DH_HAT(L2H(2)) = MAX(H0, HT_NODE(L2H(2))-0.5*DH_HAT(L2H(1)))
     &                     - HT_NODE(L2H(2))
            DH_HAT(L2H(3)) = - DH_HAT(L2H(1)) - DH_HAT(L2H(2))

C...........Adjust nodal heights, find location of maximum, sum nodal
C...........fluxes on dry nodes, and count new # of dry nodes

            NDRYNODE = 0
            QX_SUM = 0
            QY_SUM = 0
            HT_MAX = 0
            DO I = 1,3
              HT_HAT(I)  = HT_NODE(I) + DH_HAT(I)
              QX_NODE(I) = QX(1,J,IRK+1)
              QY_NODE(I) = QY(1,J,IRK+1)
              IF (HT_HAT(I).GT.HT_MAX) THEN
                HT_MAX_I = I
                HT_MAX = HT_HAT(I)
              ENDIF
              DO K = 2,DOF
                QX_NODE(I) = QX_NODE(I) + PHI_CORNER(K,I)*QX(K,J,IRK+1)
                QY_NODE(I) = QY_NODE(I) + PHI_CORNER(K,I)*QY(K,J,IRK+1)
              ENDDO
              IF (HT_HAT(I).LE.(H0+1.D-10)) THEN
                NDRYNODE = NDRYNODE + 1
                QX_SUM = QX_SUM + QX_NODE(I)
                QY_SUM = QY_SUM + QY_NODE(I)
              ENDIF
            ENDDO

C...........If all adjusted nodes are dry then zero out fluxes

            IF (NDRYNODE.EQ.3) THEN
              QX(:,J,IRK+1) = 0.D0
              QY(:,J,IRK+1) = 0.D0

C...........Else re-distribute nodal fluxes

            ELSE
              QX_SUM = QX_SUM/(3.D0-REAL(NDRYNODE))
              QY_SUM = QY_SUM/(3.D0-REAL(NDRYNODE))
              DO I = 1,3
                IF(HT_HAT(I).LE.(H0+1.D-10)) THEN
                  QX_HAT(I) = 0.D0
                  QY_HAT(I) = 0.D0
                ELSE
                  IF(QX_SUM*QX_NODE(I).LT.0.D0) THEN
                    QX_HAT(I) = QX_NODE(I) + QX_SUM
                  ELSE
                    QX_HAT(I) = QX_NODE(I)
                  ENDIF
                  IF(QY_SUM*QY_NODE(I).LT.0.D0) THEN
                    QY_HAT(I) = QY_NODE(I) + QY_SUM
                  ELSE
                    QY_HAT(I) = QY_NODE(I)
                  ENDIF
               ENDIF
              ENDDO

C.............Compute new linear modal dofs for x-direction flux
c              write(878,*) 'before ',j,qx(1,j,irk+1),qx(2,j,irk+1),
c     $             qx(3,j,irk+1)
              QX(1,J,IRK+1) =  C13*(QX_HAT(1) + QX_HAT(2) + QX_HAT(3))
              QX(2,J,IRK+1) = -C16*(QX_HAT(1) + QX_HAT(2))+C13*QX_HAT(3)
              QX(3,J,IRK+1) = -0.5D0*QX_HAT(1) + 0.5D0*QX_HAT(2)
c              write(878,*) 'after ',j,qx(1,j,irk+1),qx(2,j,irk+1),
c     $             qx(3,j,irk+1)

C.............Compute new linear modal dofs for y-direction flux

              QY(1,J,IRK+1) =  C13*(QY_HAT(1) + QY_HAT(2) + QY_HAT(3))
              QY(2,J,IRK+1) = -C16*(QY_HAT(1) + QY_HAT(2))+C13*QY_HAT(3)
              QY(3,J,IRK+1) = -0.5D0*QY_HAT(1) + 0.5D0*QY_HAT(2)
            ENDIF

C...........Compute new nodal surface elevations

            ZE_HAT(1:3) = HT_HAT(1:3) - DP_NODE(1:3,J)

C...........Compute the new linear modal dofs for the surface elevation

            ZE(2,J,IRK+1) = -C16*(ZE_HAT(1) + ZE_HAT(2)) + C13*ZE_HAT(3)
            ZE(3,J,IRK+1) = -0.5D0*ZE_HAT(1) + 0.5D0*ZE_HAT(2)
          ENDIF
          
C.........If applicable zero out higher dofs

          DO K = 4,DOF
            ZE(K,J,IRK+1) = 0.D0
            QX(K,J,IRK+1) = 0.D0
            QY(K,J,IRK+1) = 0.D0
          ENDDO
          
C.........Set flag to do element check

          ELEMENT_CHECK = 1
          
        ENDIF

C-----------------------------------------------------------------------
C     End of cases
C-----------------------------------------------------------------------

C.......Element check:
C.......Not clear to me what Shintaro is checking here (EJK)

        IF (ELEMENT_CHECK.EQ.1) THEN
          IF (ZE_HAT(HT_MAX_I).GT.(-DPE_MIN(J)+H0)) THEN

C.........Make no adjustments and set element and node flags to wet

            WDFLG(J) = 1
cnd  local time-stepping check
c  if element was previously dry and is now wet, set its lts to smallest dt
c
            if (LTSModflag.eq.2) then
               if (wdflg_tmp(j).eq.0.and.tsgroup(j).ne.1) then
                  tsgroupflag(j)=.true.
                  callresetflag2=1
               endif
            endif
cnd
            NODECODE(NM1) = 1
            NODECODE(NM2) = 1
            NODECODE(NM3) = 1
          ELSE

C...........Make adjustments and set element node codes to dry

            QX(:,J,IRK+1) = 0.D0
            QY(:,J,IRK+1) = 0.D0
            WDFLG(J) = 0
          ENDIF
        ENDIF
100   CONTINUE

      RETURN
      END SUBROUTINE WETDRY
