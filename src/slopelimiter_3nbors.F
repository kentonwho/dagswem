      SUBROUTINE SLOPELIMITER3_3NBORS(L)

C.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG


      IMPLICIT NONE

C.....Declare local variables

      INTEGER,INTENT(IN) ::  L
      INTEGER GED, II
      INTEGER EL2, EL3, EL_N, EL_N1, EL_N2
      
      REAL(SZ) DTM, DZEDR, DZEDS, DQXDR, DQXDS, DQYDR, DQYDS
      REAL(SZ) DZEDX(4), DZEDY(4)
      REAL(SZ) DQXDX(4), DQXDY(4)
      REAL(SZ) DQYDX(4), DQYDY(4)
      REAL(SZ) GRADZE(4), GRADQX(4), GRADQY(4)
      REAL(SZ) SL1, SL2
      REAL(SZ) ZELIM, ZELIMX, ZELIMY
      REAL(SZ) QXLIM, QXLIMX, QXLIMY
      REAL(SZ) QYLIM, QYLIMX, QYLIMY
      REAL(SZ) XB(0:3), YB(0:3)
      REAL(SZ) XM(1:3), YM(1:3)
      REAL(SZ) ZB(0:3), DB(0:3), QXB(0:3), QYB(0:3)    ! Values at the barycenters
      REAL(SZ) ZC(1:3), DC(1:3), QXC(1:3), QYC(1:3)  ! Values at the corner nodes
      REAL(SZ) ZM(1:3), DM(1:3), QXM(1:3), QYM(1:3)    ! Values at the edge midpoints
      REAL(SZ) HTB ! Total height at the barycenter
      REAL(SZ) NXX(3), NYY(3)
      REAL(SZ) ALPHA1(3), ALPHA2(3)
      REAL(SZ) A(2,2), B(2)
      REAL(SZ) TE(3,3)
      REAL(SZ) WB_O(3,0:2), WM_O(3)     ! Variable vectors in original space
      REAL(SZ) WB_C(3,0:2), WM_C(3) ! Variable vectors in characteristic space
      REAL(SZ) WC_O(3,3), WC_C(3,3) ! Variable vectors in original space at the corner nodes
      REAL(SZ) W_TILDA(3), DELTA_W(3)
      REAL(SZ) R_SQ ! Squared r
      REAL(SZ) HMIN, HMAX, HAVG, HLIMIT
      REAL(SZ) S
      REAL(SZ) DELTA_O(3,3), DELTA_C(3,3), DELTA_HAT(3,3)
      INTEGER N(3), NOT_LIMITED(3)
      REAL(SZ) POS, NEG, THETA_P, THETA_N
      LOGICAL LIMITIT(3)
      REAL(SZ), PARAMETER :: ZERO = 1.D-15

C.....Do nothing if this or a neighboring element is dry.
C      IF(WDFLG(L).EQ.0) RETURN
C      DO I = 1,3
C        EL_N = EL_NBORS(I,L)
C        IF(WDFLG(EL_N).EQ.0) RETURN
C      ENDDO
      
C.....Retrieve the barycenter info of this element
      XB(0) = XBC(L)
      YB(0) = YBC(L)
      ZB(0) = ZE(1,L,IRK+1)
      DB(0) = HB(1,L,1)
      QXB(0) = QX(1,L,IRK+1)
      QYB(0) = QY(1,L,IRK+1)
      
C.....Retrieve the nodal info

      N(1) = NM(L,1)
      N(2) = NM(L,2)
      N(3) = NM(L,3)

      DO K = 1,3
        ZC(K) = ZE(1,L,IRK+1)
        DC(K) = HB(1,L,1)
        QXC(K) = QX(1,L,IRK+1)
        QYC(K) = QY(1,L,IRK+1)
        DO KK = 2,3  ! This slope limiter assumes the functions to limit as a linear function.
          ZC(K) = ZC(K) + PHI_CORNER(KK,K)*ZE(KK,L,IRK+1)
          DC(K) = DC(K) + PHI_CORNER(KK,K)*HB(KK,L,1)
          QXC(K) = QXC(K) + PHI_CORNER(KK,K)*QX(KK,L,IRK+1)
          QYC(K) = QYC(K) + PHI_CORNER(KK,K)*QY(KK,L,IRK+1)
        ENDDO
      ENDDO

      HMIN = MAX(0.D0,MIN(ZC(1)+DC(1),ZC(2)+DC(2),ZC(3)+DC(3)))
      HMAX = MAX(ZC(1)+DC(1),ZC(2)+DC(2),ZC(3)+DC(3))
      HAVG = ZE(1,L,IRK+1) + HB(1,L,1)
      HLIMIT = MIN(1.D0,(HMIN*HMIN)/(SL3_MD*SL3_MD))
      HLIMIT = 1.D0
      IF(nolifa.eq.2.and.HMIN.LE.H0) RETURN
c      IF(HMIN.LE.0.01D0) RETURN

C.....Retrieve info at barycenters and edge midpoints
      DO I = 1,3
        EL_N = EL_NBORS(I,L)

C.......Edge midpoints
        N1 = N(MOD(I+0,3)+1)
        N2 = N(MOD(I+1,3)+1)
        XM(I) = 0.5D0*(X(N1) + X(N2))
        YM(I) = 0.5D0*(Y(N1) + Y(N2))

        N1 = MOD(I+0,3)+1
        N2 = MOD(I+1,3)+1
        ZM(I) = 0.5D0*(ZC(N1) + ZC(N2))
        DM(I) = 0.5D0*(DC(N1) + DC(N2))
        QXM(I) = 0.5D0*(QXC(N1) + QXC(N2))
        QYM(I) = 0.5D0*(QYC(N1) + QYC(N2))
        
C.......Barycenters
C        IF(EL_N == 0) THEN
C If there is no neighboring element on this side (e.g. domain boundary),
C locate a barycenter at the mirroring point and use extrapolated values.
C but all elements in this subroutine will have 3 neighbors
C          XB(I) = 2.D0*XM(I) - XB(0)
C          YB(I) = 2.D0*YM(I) - YB(0)
C          ZB(I) = 2.D0*ZM(I) - ZB(0)
C          DB(I) = 2.D0*DM(I) - DB(0)
C          QXB(I) = 2.D0*QXM(I) - QXB(0)
C          QYB(I) = 2.D0*QYM(I) - QYB(0)
C        ELSE
          XB(I) = XBC(EL_N)
          YB(I) = YBC(EL_N)
          ZB(I) = ZE(1,EL_N,IRK+1)
          DB(I) = HB(1,EL_N,1)
          QXB(I) = QX(1,EL_N,IRK+1)
          QYB(I) = QY(1,EL_N,IRK+1)
C        ENDIF

C.......Edge normal vector
        GED = NELED(I,L)
C.......Shouldn't be normal vector (EJK) see Cockburn and Shu reference
C        NXX(I) = COSNX(GED)
C        NYY(I) = SINNX(GED)
        NXX(I) = XM(I)-XB(I)
        IF (NXX(I).NE.0) NXX(I) = NXX(I)/ABS(NXX(I))
        NYY(I) = YM(I)-YB(I)
        IF (NYY(I).NE.0) NYY(I) = NYY(I)/ABS(NYY(I))
      ENDDO

C.....Compute alpha1 and alpha2      

      DO I = 1,3
        EL_N1 = I
        EL_N2 = MOD(I,3)+1
        A(1,1) = XB(EL_N1) - XB(0)
        A(1,2) = XB(EL_N2) - XB(0)
        A(2,1) = YB(EL_N1) - YB(0)
        A(2,2) = YB(EL_N2) - YB(0)
        B(1) = XM(I) - XB(0)
        B(2) = YM(I) - YB(0)
        DTM = A(1,1)*A(2,2) - A(1,2)*A(2,1)  ! Compute the determinant
        
C.......Compute the geometric factors alpha1 and alpha2. If either one
C.......is negative that element must be set to p = 0.

        ALPHA1(I) = ( A(2,2)*B(1) - A(1,2)*B(2))/DTM
        IF (ALPHA1(I).LT.0) THEN
          ZE(2:DOF,L,IRK+1) = 0.D0
          QX(2:DOF,L,IRK+1) = 0.D0
          QY(2:DOF,L,IRK+1) = 0.D0
          RETURN
        ENDIF
        ALPHA2(I) = (-A(2,1)*B(1) + A(1,1)*B(2))/DTM
        IF (ALPHA2(I).LT.0) THEN
          ZE(2:DOF,L,IRK+1) = 0.D0
          QX(2:DOF,L,IRK+1) = 0.D0
          QY(2:DOF,L,IRK+1) = 0.D0
          RETURN
        ENDIF
      ENDDO
      
C.....Set the variable vectors at the baricenter of element 0

      WB_O(1,0) = ZB(0)
      WB_O(2,0) = QXB(0)
      WB_O(3,0) = QYB(0)

C.....Start computing deltas for each edge

      DO I = 1,3

C.......Compute eigen values and vectors at the midpoint using the variables at the barycenter

        HTB = ZB(0)+DB(0)
        CALL HYDRO_EIGEN_VALUES(HTB, QXB(0)/HTB, QYB(0)/HTB,
     &       NXX(I), NYY(I), EIGVAL, RI, LE)

C.......Set variable vectors

        WM_O(1) = ZM(I)
        WM_O(2) = QXM(I)
        WM_O(3) = QYM(I)
        EL_N1 = I
        EL_N2 = MOD(I,3)+1
        WB_O(1,1) = ZB(EL_N1)
        WB_O(2,1) = QXB(EL_N1)
        WB_O(3,1) = QYB(EL_N1)
        
        WB_O(1,2) = ZB(EL_N2)
        WB_O(2,2) = QXB(EL_N2)
        WB_O(3,2) = QYB(EL_N2)

C.......Transform original variables into the characteristic space

        DO J = 0,2
          WB_C(1,J)
     &       = LE(1,1)*WB_O(1,J) + LE(1,2)*WB_O(2,J) + LE(1,3)*WB_O(3,J)
          WB_C(2,J)
     &       = LE(2,1)*WB_O(1,J) + LE(2,2)*WB_O(2,J) + LE(2,3)*WB_O(3,J)
          WB_C(3,J)
     &       = LE(3,1)*WB_O(1,J) + LE(3,2)*WB_O(2,J) + LE(3,3)*WB_O(3,J)
        ENDDO
        WM_C(1)= LE(1,1)*WM_O(1) + LE(1,2)*WM_O(2) + LE(1,3)*WM_O(3)
        WM_C(2)= LE(2,1)*WM_O(1) + LE(2,2)*WM_O(2) + LE(2,3)*WM_O(3)
        WM_C(3)= LE(3,1)*WM_O(1) + LE(3,2)*WM_O(2) + LE(3,3)*WM_O(3)

C.......Compute W_TILDA
        W_TILDA(1) = WM_C(1) - WB_C(1,0)
        W_TILDA(2) = WM_C(2) - WB_C(2,0)
        W_TILDA(3) = WM_C(3) - WB_C(3,0)
        
C.......Compute DELTA_W
        DELTA_W(1)
     &       = ALPHA1(I)*(WB_C(1,1) - WB_C(1,0))
     &       + ALPHA2(I)*(WB_C(1,2) - WB_C(1,0))
        DELTA_W(2)
     &       = ALPHA1(I)*(WB_C(2,1) - WB_C(2,0))
     &       + ALPHA2(I)*(WB_C(2,2) - WB_C(2,0))
        DELTA_W(3)
     &       = ALPHA1(I)*(WB_C(3,1) - WB_C(3,0))
     &       + ALPHA2(I)*(WB_C(3,2) - WB_C(3,0))

C.......Apply the TVB modified minmod function
C.......Switched to minmod function (EJK)

C        R_SQ = (XM(I) - XB(0))*(XM(I) - XB(0))
C     &       + (YM(I) - YB(0))*(YM(I) - YB(0))

        NOT_LIMITED = 0
        DO K = 1,3 ! Loop for components of variable vectors
C          IF(ABS(W_TILDA(K)) < SL2_M*HLIMIT*R_SQ) THEN
C            DELTA_C(K,I) = W_TILDA(K)
C          ELSE
            IF((W_TILDA(K)*DELTA_W(K) >= 0.D0)
     &                               .AND.(ABS(W_TILDA(K)) > ZERO)) THEN
              S = W_TILDA(K)/ABS(W_TILDA(K))
C.......Hard-wired SL2_NYU = 1.5 as taken in Cockburn and Shu reference
              DELTA_C(K,I) = S*MIN(
     &                          ABS(W_TILDA(K)),ABS(1.5*DELTA_W(K)))
              IF (DELTA_C(K,I).EQ.(S*ABS(W_TILDA(K))))
     &          NOT_LIMITED(K) = NOT_LIMITED(K) + 1
            ELSE
              DELTA_C(K,I) = 0.D0
            ENDIF
C          ENDIF
        ENDDO
      ENDDO ! Loop for edges ends here
      
C.....If none of the variables were limited then end
      
      IF ((NOT_LIMITED(1).EQ.3).AND.(NOT_LIMITED(2).EQ.3)
     &                         .AND.(NOT_LIMITED(3).EQ.3)) THEN
        RETURN
        
C.....Else perform limiting on characteristic variables

      ELSE
        DO I = 1,3 ! Loop for variable components
          POS = 0.D0
          NEG = 0.D0
          DO J = 1,3 ! Loop for edges
            POS = POS + MAX(0.D0, DELTA_C(I,J))
            NEG = NEG + MAX(0.D0,-DELTA_C(I,J))
          ENDDO
          IF((POS < ZERO).OR.(NEG < ZERO)) THEN
            DELTA_HAT(I,1) = 0.D0
            DELTA_HAT(I,2) = 0.D0
            DELTA_HAT(I,3) = 0.D0
          ELSE
            THETA_P = MIN(1.D0,NEG/POS)
            THETA_N = MIN(1.D0,POS/NEG)
            DELTA_HAT(I,1) = THETA_P*MAX(0.D0, DELTA_C(I,1))
     &                     - THETA_N*MAX(0.D0,-DELTA_C(I,1))
            DELTA_HAT(I,2) = THETA_P*MAX(0.D0, DELTA_C(I,2))
     &                     - THETA_N*MAX(0.D0,-DELTA_C(I,2))
            DELTA_HAT(I,3) = THETA_P*MAX(0.D0, DELTA_C(I,3))
     &                     - THETA_N*MAX(0.D0,-DELTA_C(I,3))
          ENDIF

C.........Computed limited characteristic variables at the corner nodes

          WC_C(I,1) = WB_C(I,0) - DELTA_HAT(I,1) + DELTA_HAT(I,2)
     &                                           + DELTA_HAT(I,3)
          WC_C(I,2) = WB_C(I,0) + DELTA_HAT(I,1) - DELTA_HAT(I,2)
     &                                           + DELTA_HAT(I,3)
          WC_C(I,3) = WB_C(I,0) + DELTA_HAT(I,1) + DELTA_HAT(I,2)
     &                                           - DELTA_HAT(I,3)

C.........Transform back to the original space

          WC_O(1,I) = RI(1,1)*WC_C(1,I) + RI(1,2)*WC_C(2,I)
     &                                  + RI(1,3)*WC_C(3,I)
          WC_O(2,I) = RI(2,1)*WC_C(1,I) + RI(2,2)*WC_C(2,I)
     &                                  + RI(2,3)*WC_C(3,I)
          WC_O(3,I) = RI(3,1)*WC_C(1,I) + RI(3,2)*WC_C(2,I)
     &                                  + RI(3,3)*WC_C(3,I)
        ENDDO ! End variable component loop
        
C.......If the limited water depth < 0, don't apply this slope limiting

        IF((WC_O(1,1)+DC(1)).LE.0.D0.OR.(WC_O(1,2)+DC(2)).LE.0.D0.OR.
     &     (WC_O(1,2)+DC(2)).LE.0.D0) RETURN
     
C.......Else compute new modal dofs

        ZE(2,L,IRK+1) = -1.D0/6.D0*(WC_O(1,1) + WC_O(1,2))
     &                +  1.D0/3.D0*WC_O(1,3)
        ZE(3,L,IRK+1) = -0.5D0*WC_O(1,1) + 0.5D0*WC_O(1,2)
        QX(2,L,IRK+1) = -1.D0/6.D0*(WC_O(2,1) + WC_O(2,2))
     &                +  1.D0/3.D0*WC_O(2,3)
        QX(3,L,IRK+1) = -0.5D0*WC_O(2,1) + 0.5D0*WC_O(2,2)
        QY(2,L,IRK+1) = -1.D0/6.D0*(WC_O(3,1) + WC_O(3,2))
     &                +  1.D0/3.D0*WC_O(3,3)
        QY(3,L,IRK+1) = -0.5D0*WC_O(3,1) + 0.5D0*WC_O(3,2)

      ENDIF

      RETURN
      END SUBROUTINE
