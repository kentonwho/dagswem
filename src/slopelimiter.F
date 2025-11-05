C***********************************************************************
C
C     SUBROUTINE SLOPELIMITER()
C
C     This subroutine selects either SLOPELIMITER1() or SLOPELIMITER2()
C     according to SLOPEFLAG. SLOPEFLAG is specified in fort.dg.
C
C***********************************************************************

      SUBROUTINE SLOPELIMITER(lvl)

      USE DG, ONLY : SLOPEFLAG
      
      IMPLICIT NONE

! CND
      integer lvl

      IF (SLOPEFLAG.EQ.1) THEN
        CALL SLOPELIMITER1()
      ELSE IF (SLOPEFLAG.EQ.2) THEN
        CALL SLOPELIMITER2()
      ELSE IF (SLOPEFLAG.EQ.3) THEN
        CALL SLOPELIMITER3()
      ELSE IF (SLOPEFLAG.EQ.5) THEN
        CALL SLOPELIMITER5(lvl)
      ENDIF

      RETURN
      END SUBROUTINE

C***********************************************************************
C
C     SUBROUTINE SLOPELIMITER1()
C
C     Written by Ethan Kubatko
C
C     08-Feb-2008
C       - Modified to apply this slopelimiter on domain 
C         boundaries. S.B.
C
C***********************************************************************

      SUBROUTINE SLOPELIMITER1()

C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER EL2, EL3, L
      
      REAL(SZ) DEN, DZEDR, DZEDS, DQXDR, DQXDS, DQYDR, DQYDS
      REAL(SZ) DZEDX(4), DZEDY(4)
      REAL(SZ) DQXDX(4), DQXDY(4)
      REAL(SZ) DQYDX(4), DQYDY(4)
      REAL(SZ) GRADZE(4), GRADQX(4), GRADQY(4)
      REAL(SZ) SL1, SL2
      REAL(SZ) ZELIM, ZELIMX, ZELIMY
      REAL(SZ) QXLIM, QXLIMX, QXLIMY
      REAL(SZ) QYLIM, QYLIMX, QYLIMY

C.....Save the original values  02/28/2007 sb
      DO L=1,NE
        DO K = 1,DOF
          ZE(K,L,NRK+2) = ZE(K,L,IRK+1)
          QX(K,L,NRK+2) = QX(K,L,IRK+1)
          QY(K,L,NRK+2) = QY(K,L,IRK+1)
        ENDDO
      ENDDO      

      DO 1000 L=1, ne

C.....Retrieve the barycenter coordinates of the element
      
      X1 = XBC(L)
      Y1 = YBC(L)

C.....Compute the x and y derivatives of the P1 part of ZE

      DZEDR = ZE(3,L,NRK+2)
      DZEDS = 3.D0/2.D0*ZE(2,L,NRK+2) + 1.D0/2.D0*ZE(3,L,NRK+2)

      DZEDX(4) = DZEDR*DRDX(L) + DZEDS*DSDX(L)
      DZEDY(4) = DZEDR*DRDY(L) + DZEDS*DSDY(L)
      
      GRADZE(4) = SQRT(DZEDX(4)**2.D0 + DZEDY(4)**2.D0)
      
C.....Compute the x and y derivatives of the P1 part of QX
      
      DQXDR = QX(3,L,NRK+2)
      DQXDS = 3.D0/2.D0*QX(2,L,NRK+2) + 1.D0/2.D0*QX(3,L,NRK+2)

      DQXDX(4) = DQXDR*DRDX(L) + DQXDS*DSDX(L)
      DQXDY(4) = DQXDR*DRDY(L) + DQXDS*DSDY(L)
      
      GRADQX(4) = SQRT(DQXDX(4)**2.D0 + DQXDY(4)**2.D0)
      
C.....Compute the x and y derivatives of the P1 part of QY
      
      DQYDR = QY(3,L,NRK+2)
      DQYDS = 3.D0/2.D0*QY(2,L,NRK+2) + 1.D0/2.D0*QY(3,L,NRK+2)

      DQYDX(4) = DQYDR*DRDX(L) + DQYDS*DSDX(L)
      DQYDY(4) = DQYDR*DRDY(L) + DQYDS*DSDY(L)
      
      GRADQY(4) = SQRT(DQYDX(4)**2.D0 + DQYDY(4)**2.D0)
      
      DO I = 1,3
      
        EL2 = EL_NBORS(I,L)
        EL3 = EL_NBORS(I+1,L)

C.......Special treatment for elements that do not have three neighboring elements
        IF((EL2.EQ.0).OR.(EL3.EQ.0).OR.(SL3(I,L).LE.0.D0)) THEN
          GRADZE(I) = 1.D10
          GRADQX(I) = 1.D10
          GRADQY(I) = 1.D10
          CYCLE
        ENDIF
      
        X2 = XBC(EL2)
        Y2 = YBC(EL2)
        
        X3 = XBC(EL3)
        Y3 = YBC(EL3)
        
C.....Compute limiting slopes for the surface elevation

        SL1 = ZE(1,L,NRK+2)*(Y3 - Y2) + ZE(1,EL2,NRK+2)*(Y1 - Y3)
     &                                + ZE(1,EL3,NRK+2)*(Y2 - Y1)
     
        SL2 = ZE(1,L,NRK+2)*(X2 - X3) + ZE(1,EL2,NRK+2)*(X3 - X1)
     &                                + ZE(1,EL3,NRK+2)*(X1 - X2)
     
        DZEDX(I) = -SL1/SL3(I,L)
        
        DZEDY(I) = -SL2/SL3(I,L)
        
        GRADZE(I) = SQRT(DZEDX(I)**2.D0 + DZEDY(I)**2.D0)
        
C.....Compute limiting slopes for the X flow

        SL1 = QX(1,L,NRK+2)*(Y3 - Y2) + QX(1,EL2,NRK+2)*(Y1 - Y3)
     &                                + QX(1,EL3,NRK+2)*(Y2 - Y1)

        SL2 = QX(1,L,NRK+2)*(X2 - X3) + QX(1,EL2,NRK+2)*(X3 - X1)
     &                                + QX(1,EL3,NRK+2)*(X1 - X2)

        DQXDX(I) = -SL1/SL3(I,L)

        DQXDY(I) = -SL2/SL3(I,L)
        
        GRADQX(I) = SQRT(DQXDX(I)**2.D0 + DQXDY(I)**2.D0)
        
C.....Compute limiting slopes for the Y flow

        SL1 = QY(1,L,NRK+2)*(Y3 - Y2) + QY(1,EL2,NRK+2)*(Y1 - Y3)
     &                                + QY(1,EL3,NRK+2)*(Y2 - Y1)

        SL2 = QY(1,L,NRK+2)*(X2 - X3) + QY(1,EL2,NRK+2)*(X3 - X1)
     &                                + QY(1,EL3,NRK+2)*(X1 - X2)

        DQYDX(I) = -SL1/SL3(I,L)

        DQYDY(I) = -SL2/SL3(I,L)
        
        GRADQY(I) = SQRT(DQYDX(I)**2.D0 + DQYDY(I)**2.D0)
        
      ENDDO
      
C.....Compute the (possibly) limited gradient

      ZELIM = MIN( GRADZE(1), GRADZE(2), GRADZE(3), GRADZE(4) )
      QXLIM = MIN( GRADQX(1), GRADQX(2), GRADQX(3), GRADQX(4) )
      QYLIM = MIN( GRADQY(1), GRADQY(2), GRADQY(3), GRADQY(4) )
      
      DO I = 1,4
      
        IF (ZELIM.EQ.GRADZE(I)) THEN
          ZELIMX = DZEDX(I)
          ZELIMY = DZEDY(I)
        ENDIF
        
        IF (QXLIM.EQ.GRADQX(I)) THEN
          QXLIMX = DQXDX(I)
          QXLIMY = DQXDY(I)
        ENDIF
        
        IF (QYLIM.EQ.GRADQY(I)) THEN
          QYLIMX = DQYDX(I)
          QYLIMY = DQYDY(I)
        ENDIF
        
      ENDDO

      DEN = DRDX(L)*DSDY(L) - DRDY(L)*DSDX(L)

      ZE(2,L,IRK+1) = 1.D0/3.D0*(2.D0*ZELIMY*DRDX(L)
     &      - 2.D0*DRDY(L)*ZELIMX - ZELIMX*DSDY(L) + DSDX(L)*ZELIMY)/DEN
      QX(2,L,IRK+1) = 1.D0/3.D0*(2.D0*QXLIMY*DRDX(L)
     &      - 2.D0*DRDY(L)*QXLIMX - QXLIMX*DSDY(L) + DSDX(L)*QXLIMY)/DEN
      QY(2,L,IRK+1) = 1.D0/3.D0*(2.D0*QYLIMY*DRDX(L)
     &      - 2.D0*DRDY(L)*QYLIMX - QYLIMX*DSDY(L) + DSDX(L)*QYLIMY)/DEN

      ZE(3,L,IRK+1) = (ZELIMX*DSDY(L) - DSDX(L)*ZELIMY)/DEN
      QX(3,L,IRK+1) = (QXLIMX*DSDY(L) - DSDX(L)*QXLIMY)/DEN
      QY(3,L,IRK+1) = (QYLIMX*DSDY(L) - DSDX(L)*QYLIMY)/DEN
      
C      ZE(2:DOF,L,IRK+1) = 0.D0
C      QX(2:DOF,L,IRK+1) = 0.D0
C      QY(2:DOF,L,IRK+1) = 0.D0
     
C      IF (ZELIM.EQ.GRADZE(4)) THEN
C        WRITE(85,*) 'DID NOT LIMIT ZE ON ELEMENT #',L
C      ENDIF
      
C      IF (QXLIM.EQ.GRADQX(4)) THEN
C        WRITE(85,*) 'DO NOT LIMIT QX ON ELEMENT #',L
C      ENDIF
      
C      IF (QYLIM.EQ.GRADQY(4)) THEN
C        WRITE(85,*) 'DID NOT LIMIT QY ON ELEMENT #',L
C      ENDIF

111   continue

1000  CONTINUE
      RETURN
      END SUBROUTINE

      
C***********************************************************************
C
C     SUBROUTINE SLOPELIMITER2()
C
C     Written by Shintaro Bunya - 08 Feb 2008
C
C     This subroutine is an implementation of the slope limiter used in
C     B. Cockburn and C-W. Shu, "The Runge-Kutta Discontinuous Galerkin
C     Method for Conservation Laws V, Multidimensional Systems,"
C     Journal of Computational Physics 141, 199-224 (1998).
C
C***********************************************************************

      SUBROUTINE SLOPELIMITER2()

C.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, GED, II
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
      REAL(SZ) WC_O(3,3) ! Variable vectors in original space at the corner nodes
      REAL(SZ) W_TILDA(3), DELTA_W(3)
      REAL(SZ) R_SQ ! Squared r
      REAL(SZ) S
      REAL(SZ) DELTA_O(3,3), DELTA_C(3,3), DELTA_HAT(3,3)
      INTEGER N(3)
      REAL(SZ) POS, NEG, THETA_P, THETA_N
      LOGICAL LIMITIT(3)
      REAL(SZ), PARAMETER :: ZERO = 1.D-15

C.....Save the original values  02/28/2007 sb
      DO L=1,NE
        DO K = 1,DOF
          ZE(K,L,NRK+2) = ZE(K,L,IRK+1)
          QX(K,L,NRK+2) = QX(K,L,IRK+1)
          QY(K,L,NRK+2) = QY(K,L,IRK+1)
        ENDDO
      ENDDO      

      DO 1000 L=1, NE

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
        IF(EL_N == 0) THEN
C If there is no neighboring element on this side (e.g. domain boundary),
C locate a barycenter at the mirroring point and use the values at the 
C baricenter of this element.
          XB(I) = 2.D0*XM(I) - XB(0)
          YB(I) = 2.D0*YM(I) - YB(0)
          ZB(I) = ZB(0)
          DB(I) = DB(0)
          QXB(I) = QXB(0)
          QYB(I) = QYB(0)
        ELSE
          XB(I) = XBC(EL_N)
          YB(I) = YBC(EL_N)
          ZB(I) = ZE(1,EL_N,IRK+1)
          DB(I) = HB(1,EL_N,1)
          QXB(I) = QX(1,EL_N,IRK+1)
          QYB(I) = QY(1,EL_N,IRK+1)
        ENDIF

C.......Edge normal vector
        GED = NELED(I,L)
        NXX(I) = COSNX(GED)
        NYY(I) = SINNX(GED)
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
        ALPHA1(I) = ( A(2,2)*B(1) - A(1,2)*B(2))/DTM
        ALPHA2(I) = (-A(2,1)*B(1) + A(1,1)*B(2))/DTM
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
        R_SQ = (XM(I) - XB(0))*(XM(I) - XB(0))
     &       + (YM(I) - YB(0))*(YM(I) - YB(0))

        DO K = 1,3 ! Loop for components of variable vectors
          IF(ABS(W_TILDA(K)) < SL2_M*R_SQ) THEN
            DELTA_C(K,I) = W_TILDA(K)
          ELSE
            IF((W_TILDA(K)*DELTA_W(K) >= 0.D0)
     &                               .AND.(ABS(W_TILDA(K)) > ZERO)) THEN
              S = W_TILDA(K)/ABS(W_TILDA(K))
              DELTA_C(K,I) = S*MIN(
     &                          ABS(W_TILDA(K)),ABS(SL2_NYU*DELTA_W(K)))
            ELSE
              DELTA_C(K,I) = 0.D0
            ENDIF
          ENDIF
        ENDDO

C.......Transform back to the original space
        DELTA_O(1,I) = RI(1,1)*DELTA_C(1,I) + RI(1,2)*DELTA_C(2,I)
     &                                      + RI(1,3)*DELTA_C(3,I)
        DELTA_O(2,I) = RI(2,1)*DELTA_C(1,I) + RI(2,2)*DELTA_C(2,I)
     &                                      + RI(2,3)*DELTA_C(3,I)
        DELTA_O(3,I) = RI(3,1)*DELTA_C(1,I) + RI(3,2)*DELTA_C(2,I)
     &                                      + RI(3,3)*DELTA_C(3,I)

      ENDDO ! Loop for edges ends here

C Set the variable vector at the midpoint 1. This is used later to judge
C whether the variables should be limited.
      WM_O(1) = ZM(1)
      WM_O(2) = QXM(1)
      WM_O(3) = QYM(1)

C Find the possibly limited DELTA and variable vectors at the corner nodes
      DO I = 1,3 ! Loop for variable components
C Find the possibly limited DELTA
        POS = 0.D0
        NEG = 0.D0
        DO J = 1,3 ! Loop for edges
          POS = POS + MAX(0.D0, DELTA_O(I,J))
          NEG = NEG + MAX(0.D0,-DELTA_O(I,J))
        ENDDO
        IF((POS < ZERO).OR.(NEG < ZERO)) THEN 
          DELTA_HAT(I,1) = 0.D0
          DELTA_HAT(I,2) = 0.D0
          DELTA_HAT(I,3) = 0.D0
        ELSE
          THETA_P = MIN(1.D0,NEG/POS)
          THETA_N = MIN(1.D0,POS/NEG)
          
          DELTA_HAT(I,1)
     &         = THETA_P*MAX(0.D0, DELTA_O(I,1))
     &         - THETA_N*MAX(0.D0,-DELTA_O(I,1))
          DELTA_HAT(I,2)
     &         = THETA_P*MAX(0.D0, DELTA_O(I,2))
     &         - THETA_N*MAX(0.D0,-DELTA_O(I,2))
          DELTA_HAT(I,3)
     &         = THETA_P*MAX(0.D0, DELTA_O(I,3))
     &         - THETA_N*MAX(0.D0,-DELTA_O(I,3))
        ENDIF

C Compute the possibly limited variable vectors at the corner nodes
        IF(( DELTA_HAT(I,1) < (WM_O(I) - WB_O(I,0) - ZERO) ).OR.
     &       ( DELTA_HAT(I,1) > (WM_O(I) - WB_O(I,0) + ZERO) )) THEN
C Update the variables only if they need to be limited.
C If the condition in the if statement is false, the following 
C isn't applied and therefore the quradratic and higher order 
C coefficients are preserved.

C Computing limited variables at the corner nodes
          WC_O(I,1) = WB_O(I,0)
     &         - DELTA_HAT(I,1)
     &         + DELTA_HAT(I,2)
     &         + DELTA_HAT(I,3)
          WC_O(I,2) = WB_O(I,0)
     &         + DELTA_HAT(I,1)
     &         - DELTA_HAT(I,2)
     &         + DELTA_HAT(I,3)
          WC_O(I,3) = WB_O(I,0)
     &         + DELTA_HAT(I,1)
     &         + DELTA_HAT(I,2)
     &         - DELTA_HAT(I,3)
          LIMITIT(I) = .TRUE.
        ELSE
          LIMITIT(I) = .FALSE.
        ENDIF
      ENDDO

      IF(LIMITIT(1)) THEN
        ZE(2,L,IRK+1) = -1.D0/6.D0*(WC_O(1,1) + WC_O(1,2))
     &                + 1.D0/3.D0*WC_O(1,3)
        ZE(3,L,IRK+1) = -0.5D0*WC_O(1,1) + 0.5D0*WC_O(1,2)
      ENDIF
      IF(LIMITIT(2)) THEN
        QX(2,L,IRK+1) = -1.D0/6.D0*(WC_O(2,1) + WC_O(2,2))
     &                + 1.D0/3.D0*WC_O(2,3)
        QX(3,L,IRK+1) = -0.5D0*WC_O(2,1) + 0.5D0*WC_O(2,2)
      ENDIF
      IF(LIMITIT(3)) THEN
        QY(2,L,IRK+1) = -1.D0/6.D0*(WC_O(3,1) + WC_O(3,2))
     &                + 1.D0/3.D0*WC_O(3,3)
        QY(3,L,IRK+1) = -0.5D0*WC_O(3,1) + 0.5D0*WC_O(3,2)
      ENDIF

1000  CONTINUE

      RETURN
      END SUBROUTINE

C***********************************************************************
C
C     SUBROUTINE SLOPELIMITER3()
C
C     Written by Shintaro Bunya - 27 Feb 2008
C
C     This subroutine is an extentioin of the slope limiter used in
C     B. Cockburn and C-W. Shu, "The Runge-Kutta Discontinuous Galerkin
C     Method for Conservation Laws V, Multidimensional Systems,"
C     Journal of Computational Physics 141, 199-224 (1998).
C
C     An additional parameter SL3_MD is used in this slope limiter to
C     make sure that slopes are limited in elements whose water depth
C     is small. (A special treatment for wetting and drying.
C
C***********************************************************************

      SUBROUTINE SLOPELIMITER3()

C.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, NNBORS

      DO L= 1,NE
c         IF (WDFLG(L).NE.0) THEN
c
            NNBORS = 0
            DO I=1,3
               IF(EL_NBORS(I,L).NE.0) THEN
                  NNBORS = NNBORS + 1
               ENDIF
            ENDDO
            IF(NNBORS.EQ.3) THEN
               CALL SLOPELIMITER3_3NBORS(L)
            ELSE IF(NNBORS.EQ.2) THEN
               CALL SLOPELIMITER3_2NBORS(L)
            ELSE IF(NNBORS.EQ.1) THEN
               CALL SLOPELIMITER3_1NBORS(L)
#ifdef CMPI
C     CYCLE
#else
C     STOP 'STILL UNDER CONSTRUCTION!!!'
#endif
            ENDIF
c         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE

      SUBROUTINE SLOPELIMITER3_3NBORS(L)

C.....Use appropriate modules

      USE SIZES, ONLY : SZ,myproc
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
c      if (myproc.eq.42.and.l.eq.4229) then
c         write(42,*) 'slope limiter ',i,el_n,qxm(i),qxb(i)
c      endif
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
c        if (myproc.eq.42.and.l.eq.4229) then
c           write(42,*) 'slope limiter ',delta_w(1),delta_w(2),
c     $          delta_w(3),w_tilda(1),w_tilda(2),w_tilda(3)
c        endif
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
       enddo
C.........Transform back to the original space
       do i=1,3   !loop over corners
          WC_O(1,I) = RI(1,1)*WC_C(1,I) + RI(1,2)*WC_C(2,I)
     &                                  + RI(1,3)*WC_C(3,I)
          WC_O(2,I) = RI(2,1)*WC_C(1,I) + RI(2,2)*WC_C(2,I)
     &                                  + RI(2,3)*WC_C(3,I)
          WC_O(3,I) = RI(3,1)*WC_C(1,I) + RI(3,2)*WC_C(2,I)
     &                                  + RI(3,3)*WC_C(3,I)

c        if (myproc.eq.42.and.l.eq.4229) then
c           write(42,*) 'slope limiter ',wc_c(i,1),wc_c(i,2),wc_c(i,3)
c        endif
       ENDDO                 

        
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

      SUBROUTINE SLOPELIMITER3_2NBORS(L)

C.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER,INTENT(IN) ::  L
      INTEGER GED, II
      INTEGER EL2, EL3, EL_N, EL_N1, EL_N2, EL_NN(2), NNBORS
      
      REAL(SZ) DTM, DZEDR, DZEDS, DQXDR, DQXDS, DQYDR, DQYDS
      REAL(SZ) DZEDX(4), DZEDY(4)
      REAL(SZ) DQXDX(4), DQXDY(4)
      REAL(SZ) DQYDX(4), DQYDY(4)
      REAL(SZ) GRADZE(4), GRADQX(4), GRADQY(4)
      REAL(SZ) SL1, SL2
      REAL(SZ) ZELIM, ZELIMX, ZELIMY
      REAL(SZ) QXLIM, QXLIMX, QXLIMY
      REAL(SZ) QYLIM, QYLIMX, QYLIMY
      REAL(SZ) XB(0:2), YB(0:2)
      REAL(SZ) XM(1:3), YM(1:3)
      REAL(SZ) ZB(0:2), DB(0:2), QXB(0:2), QYB(0:2)    ! Values at the barycenters
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


C.....Find two neighboring elements
      NNBORS = 0
      DO I = 1,3
        EL_N = EL_NBORS(I,L)
        IF(EL_N.NE.0) THEN
          NNBORS = NNBORS + 1
          EL_NN(NNBORS) = EL_N
        ENDIF
      ENDDO

C.....Do nothing if a neighboring element is dry.
C      IF(WDFLG(L).EQ.0) RETURN
C      DO I = 1,NNBORS
C        EL_N = EL_NN(I)
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

C.....Retrieve info at edge midpoints
      DO I = 1,3
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
        
C.......Edge normal vector
        GED = NELED(I,L)
c        NXX(I) = -COSNX(GED)
c        NYY(I) = -SINNX(GED)
        NXX(I) = XM(I)-XB(I)
        IF (NXX(I).NE.0) NXX(I) = NXX(I)/ABS(NXX(I))
        NYY(I) = YM(I)-YB(I)
        IF (NYY(I).NE.0) NYY(I) = NYY(I)/ABS(NYY(I))
      ENDDO

C.....Retreive barycenter info
      DO I = 1,NNBORS
        EL_N = EL_NN(I)
        XB(I) = XBC(EL_N)
        YB(I) = YBC(EL_N)
        ZB(I) = ZE(1,EL_N,IRK+1)
        DB(I) = HB(1,EL_N,1)
        QXB(I) = QX(1,EL_N,IRK+1)
        QYB(I) = QY(1,EL_N,IRK+1)
      ENDDO

C.....Compute alpha1 and alpha2      
      DO I = 1,3
        EL_N1 = 1
        EL_N2 = 2
        A(1,1) = XB(EL_N1) - XB(0)
        A(1,2) = XB(EL_N2) - XB(0)
        A(2,1) = YB(EL_N1) - YB(0)
        A(2,2) = YB(EL_N2) - YB(0)
        B(1) = XM(I) - XB(0)
        B(2) = YM(I) - YB(0)
        DTM = A(1,1)*A(2,2) - A(1,2)*A(2,1)  ! Compute the determinant
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

        EL_N1 = 1
        EL_N2 = 2

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
c        R_SQ = (XM(I) - XB(0))*(XM(I) - XB(0))
c     &       + (YM(I) - YB(0))*(YM(I) - YB(0))
        NOT_LIMITED=0

        DO K = 1,3 ! Loop for components of variable vectors
c          IF(ABS(W_TILDA(K)) < SL2_M*HLIMIT*R_SQ) THEN
c            DELTA_C(K,I) = W_TILDA(K)
c          ELSE
            IF((W_TILDA(K)*DELTA_W(K) >= 0.D0)
     &                               .AND.(ABS(W_TILDA(K)) > ZERO)) THEN
              S = W_TILDA(K)/ABS(W_TILDA(K))
              DELTA_C(K,I) = S*MIN(
     &                          ABS(W_TILDA(K)),ABS(SL2_NYU*DELTA_W(K)))
              IF (DELTA_C(K,I).EQ.(S*ABS(W_TILDA(K))))
     &          NOT_LIMITED(K) = NOT_LIMITED(K) + 1
            ELSE
              DELTA_C(K,I) = 0.D0
            ENDIF
c          ENDIF
        ENDDO
      ENDDO
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
       enddo! End variable component loop
C.........Transform back to the original space
       do i=1,3

          WC_O(1,I) = RI(1,1)*WC_C(1,I) + RI(1,2)*WC_C(2,I)
     &                                  + RI(1,3)*WC_C(3,I)
          WC_O(2,I) = RI(2,1)*WC_C(1,I) + RI(2,2)*WC_C(2,I)
     &                                  + RI(2,3)*WC_C(3,I)
          WC_O(3,I) = RI(3,1)*WC_C(1,I) + RI(3,2)*WC_C(2,I)
     &                                  + RI(3,3)*WC_C(3,I)

       ENDDO 
        
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

      SUBROUTINE SLOPELIMITER3_1NBORS(L)

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
      REAL(SZ) WC_O(3,3) ! Variable vectors in original space at the corner nodes
      REAL(SZ) W_TILDA(3), DELTA_W(3)
      REAL(SZ) R_SQ ! Squared r
      REAL(SZ) HMIN, HMAX, HAVG, HLIMIT
      REAL(SZ) S
      REAL(SZ) DELTA_O(3,3), DELTA_C(3,3), DELTA_HAT(3,3)
      INTEGER N(3)
      REAL(SZ) POS, NEG, THETA_P, THETA_N
      LOGICAL LIMITIT(3)
      REAL(SZ), PARAMETER :: ZERO = 1.D-15

C.....Retrieve the barycenter info of this element
      XB(0) = XBC(L)
      YB(0) = YBC(L)
      ZB(0) = ZE(1,L,IRK+1)
      DB(0) = HB(1,L,1)
      QXB(0) = QX(1,L,IRK+1)
      QYB(0) = QY(1,L,IRK+1)
      
C.....Retrieve the nodal info
      limitit(1)=.true.
      limitit(2)=.true.
      limitit(3)=.true.
      IF(LIMITIT(1)) THEN
        ZE(2,L,IRK+1) = 0.0
        ZE(3,L,IRK+1) = 0.0
      ENDIF
      IF(LIMITIT(2)) THEN
        QX(2,L,IRK+1) = 0.0
        QX(3,L,IRK+1) = 0.0
      ENDIF
      IF(LIMITIT(3)) THEN
        QY(2,L,IRK+1) = 0.0
        QY(3,L,IRK+1) = 0.0
        return
      ENDIF

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
      IF(NOLIFA.eq.2.and.HMIN.LE.H0) RETURN

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
        IF(EL_N == 0) THEN
C If there is no neighboring element on this side (e.g. domain boundary),
C locate a barycenter at the mirroring point and use extrapolated values.
          XB(I) = 2.D0*XM(I) - XB(0)
          YB(I) = 2.D0*YM(I) - YB(0)
          ZB(I) = 2.D0*ZM(I) - ZB(0)
          DB(I) = 2.D0*DM(I) - DB(0)
          QXB(I) = 2.D0*QXM(I) - QXB(0)
          QYB(I) = 2.D0*QYM(I) - QYB(0)
        ELSE
          XB(I) = XBC(EL_N)
          YB(I) = YBC(EL_N)
          ZB(I) = ZE(1,EL_N,IRK+1)
          DB(I) = HB(1,EL_N,1)
          QXB(I) = QX(1,EL_N,IRK+1)
          QYB(I) = QY(1,EL_N,IRK+1)
        ENDIF

C.......Edge normal vector
        GED = NELED(I,L)
        NXX(I) = COSNX(GED)
        NYY(I) = SINNX(GED)
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
        ALPHA1(I) = ( A(2,2)*B(1) - A(1,2)*B(2))/DTM
        ALPHA2(I) = (-A(2,1)*B(1) + A(1,1)*B(2))/DTM
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
        R_SQ = (XM(I) - XB(0))*(XM(I) - XB(0))
     &       + (YM(I) - YB(0))*(YM(I) - YB(0))

        DO K = 1,3 ! Loop for components of variable vectors
          IF(ABS(W_TILDA(K)) < SL2_M*HLIMIT*R_SQ) THEN
            DELTA_C(K,I) = W_TILDA(K)
          ELSE
            IF((W_TILDA(K)*DELTA_W(K) >= 0.D0)
     &                               .AND.(ABS(W_TILDA(K)) > ZERO)) THEN
              S = W_TILDA(K)/ABS(W_TILDA(K))
              DELTA_C(K,I) = S*MIN(
     &                          ABS(W_TILDA(K)),ABS(SL2_NYU*DELTA_W(K)))
            ELSE
              DELTA_C(K,I) = 0.D0
            ENDIF
          ENDIF
        ENDDO

C.......Transform back to the original space
        DELTA_O(1,I)
     &       = RI(1,1)*DELTA_C(1,I)
     &       + RI(1,2)*DELTA_C(2,I)
     &       + RI(1,3)*DELTA_C(3,I)
        DELTA_O(2,I)
     &       = RI(2,1)*DELTA_C(1,I)
     &       + RI(2,2)*DELTA_C(2,I)
     &       + RI(2,3)*DELTA_C(3,I)
        DELTA_O(3,I)
     &       = RI(3,1)*DELTA_C(1,I)
     &       + RI(3,2)*DELTA_C(2,I)
     &       + RI(3,3)*DELTA_C(3,I)
        

      ENDDO ! Loop for edges ends here

C Set the variable vector at the midpoint 1. This is used later to judge
C whether the variables should be limited.
      WM_O(1) = ZM(1)
      WM_O(2) = QXM(1)
      WM_O(3) = QYM(1)

C Find the possibly limited DELTA and variable vectors at the corner nodes
      DO I = 1,3 ! Loop for variable components
C Find the possibly limited DELTA
        POS = 0.D0
        NEG = 0.D0
        DO J = 1,3 ! Loop for edges
          POS = POS + MAX(0.D0, DELTA_O(I,J))
          NEG = NEG + MAX(0.D0,-DELTA_O(I,J))
        ENDDO
        IF((POS < ZERO).OR.(NEG < ZERO)) THEN 
          DELTA_HAT(I,1) = 0.D0
          DELTA_HAT(I,2) = 0.D0
          DELTA_HAT(I,3) = 0.D0
        ELSE
          THETA_P = MIN(1.D0,NEG/POS)
          THETA_N = MIN(1.D0,POS/NEG)
          
          DELTA_HAT(I,1)
     &         = THETA_P*MAX(0.D0, DELTA_O(I,1))
     &         - THETA_N*MAX(0.D0,-DELTA_O(I,1))
          DELTA_HAT(I,2)
     &         = THETA_P*MAX(0.D0, DELTA_O(I,2))
     &         - THETA_N*MAX(0.D0,-DELTA_O(I,2))
          DELTA_HAT(I,3)
     &         = THETA_P*MAX(0.D0, DELTA_O(I,3))
     &         - THETA_N*MAX(0.D0,-DELTA_O(I,3))
        ENDIF

C Compute the possibly limited variable vectors at the corner nodes
        IF(( DELTA_HAT(I,1) < (WM_O(I) - WB_O(I,0) - ZERO) ).OR.
     &       ( DELTA_HAT(I,1) > (WM_O(I) - WB_O(I,0) + ZERO) )) THEN
C Update the variables only if they need to be limited.
C If the condition in the if statement is false, the following 
C isn't applied and therefore the quradratic and higher order 
C coefficients are preserved.

C Computing limited variables at the corner nodes
          WC_O(I,1) = WB_O(I,0)
     &         - DELTA_HAT(I,1)
     &         + DELTA_HAT(I,2)
     &         + DELTA_HAT(I,3)
          WC_O(I,2) = WB_O(I,0)
     &         + DELTA_HAT(I,1)
     &         - DELTA_HAT(I,2)
     &         + DELTA_HAT(I,3)
          WC_O(I,3) = WB_O(I,0)
     &         + DELTA_HAT(I,1)
     &         + DELTA_HAT(I,2)
     &         - DELTA_HAT(I,3)
          LIMITIT(I) = .TRUE.
        ELSE
          LIMITIT(I) = .FALSE.
        ENDIF
      ENDDO

C If the limited water depth is less 0, don't apply this slope limiting.
      IF((WC_O(1,1)+DC(1)).LE.0.D0.OR.
     &   (WC_O(1,2)+DC(2)).LE.0.D0.OR.
     &   (WC_O(1,2)+DC(2)).LE.0.D0) THEN
        GOTO 10000
      ENDIF

      IF(LIMITIT(1)) THEN
        ZE(2,L,IRK+1) = -1.D0/6.D0*(WC_O(1,1) + WC_O(1,2))
     &                + 1.D0/3.D0*WC_O(1,3)
        ZE(3,L,IRK+1) = -0.5D0*WC_O(1,1) + 0.5D0*WC_O(1,2)
      ENDIF
      IF(LIMITIT(2)) THEN
        QX(2,L,IRK+1) = -1.D0/6.D0*(WC_O(2,1) + WC_O(2,2))
     &                + 1.D0/3.D0*WC_O(2,3)
        QX(3,L,IRK+1) = -0.5D0*WC_O(2,1) + 0.5D0*WC_O(2,2)
      ENDIF
      IF(LIMITIT(3)) THEN
        QY(2,L,IRK+1) = -1.D0/6.D0*(WC_O(3,1) + WC_O(3,2))
     &                + 1.D0/3.D0*WC_O(3,3)
        QY(3,L,IRK+1) = -0.5D0*WC_O(3,1) + 0.5D0*WC_O(3,2)
      ENDIF

10000 CONTINUE

      RETURN
      END SUBROUTINE

C***********************************************************************
C
C     SUBROUTINE SLOPELIMITER3_OLD()
C
C     Written by Shintaro Bunya - 27 Feb 2008
C
C     This subroutine is an extentioin of the slope limiter used in
C     B. Cockburn and C-W. Shu, "The Runge-Kutta Discontinuous Galerkin
C     Method for Conservation Laws V, Multidimensional Systems,"
C     Journal of Computational Physics 141, 199-224 (1998).
C
C     An additional parameter SL3_MD is used in this slope limiter to
C     make sure that slopes are limited in elements whose water depth
C     is small. (A special treatment for wetting and drying.
C
C***********************************************************************

      SUBROUTINE SLOPELIMITER3_OLD()

C.....Use appropriate modules

      USE SIZES, ONLY : SZ
      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, GED, II
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
      REAL(SZ) WC_O(3,3) ! Variable vectors in original space at the corner nodes
      REAL(SZ) W_TILDA(3), DELTA_W(3)
      REAL(SZ) R_SQ ! Squared r
      REAL(SZ) HMIN, HMAX, HAVG, HLIMIT
      REAL(SZ) S
      REAL(SZ) DELTA_O(3,3), DELTA_C(3,3), DELTA_HAT(3,3)
      INTEGER N(3)
      REAL(SZ) POS, NEG, THETA_P, THETA_N
      LOGICAL LIMITIT(3)
      REAL(SZ), PARAMETER :: ZERO = 1.D-15

C.....Save the original values  02/28/2007 sb
      DO L=1,NE
        DO K = 1,DOF
          ZE(K,L,NRK+2) = ZE(K,L,IRK+1)
          QX(K,L,NRK+2) = QX(K,L,IRK+1)
          QY(K,L,NRK+2) = QY(K,L,IRK+1)
        ENDDO
      ENDDO      

      DO 1000 L=1, NE

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
      IF(HMIN.LE.H0) RETURN

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
        IF(EL_N == 0) THEN
C If there is no neighboring element on this side (e.g. domain boundary),
C locate a barycenter at the mirroring point and use extrapolated values.
          XB(I) = 2.D0*XM(I) - XB(0)
          YB(I) = 2.D0*YM(I) - YB(0)
          ZB(I) = 2.D0*ZM(I) - ZB(0)
          DB(I) = 2.D0*DM(I) - DB(0)
          QXB(I) = 2.D0*QXM(I) - QXB(0)
          QYB(I) = 2.D0*QYM(I) - QYB(0)
        ELSE
          XB(I) = XBC(EL_N)
          YB(I) = YBC(EL_N)
          ZB(I) = ZE(1,EL_N,IRK+1)
          DB(I) = HB(1,EL_N,1)
          QXB(I) = QX(1,EL_N,IRK+1)
          QYB(I) = QY(1,EL_N,IRK+1)
        ENDIF

C.......Edge normal vector
        GED = NELED(I,L)
        NXX(I) = COSNX(GED)
        NYY(I) = SINNX(GED)
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
        ALPHA1(I) = ( A(2,2)*B(1) - A(1,2)*B(2))/DTM
        ALPHA2(I) = (-A(2,1)*B(1) + A(1,1)*B(2))/DTM
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
        R_SQ = (XM(I) - XB(0))*(XM(I) - XB(0))
     &       + (YM(I) - YB(0))*(YM(I) - YB(0))

        DO K = 1,3 ! Loop for components of variable vectors
          IF(ABS(W_TILDA(K)) < SL2_M*HLIMIT*R_SQ) THEN
            DELTA_C(K,I) = W_TILDA(K)
          ELSE
            IF((W_TILDA(K)*DELTA_W(K) >= 0.D0)
     &                               .AND.(ABS(W_TILDA(K)) > ZERO)) THEN
              S = W_TILDA(K)/ABS(W_TILDA(K))
              DELTA_C(K,I) = S*MIN(
     &                          ABS(W_TILDA(K)),ABS(SL2_NYU*DELTA_W(K)))
            ELSE
              DELTA_C(K,I) = 0.D0
            ENDIF
          ENDIF
        ENDDO

C.......Transform back to the original space
        DELTA_O(1,I)
     &       = RI(1,1)*DELTA_C(1,I)
     &       + RI(1,2)*DELTA_C(2,I)
     &       + RI(1,3)*DELTA_C(3,I)
        DELTA_O(2,I)
     &       = RI(2,1)*DELTA_C(1,I)
     &       + RI(2,2)*DELTA_C(2,I)
     &       + RI(2,3)*DELTA_C(3,I)
        DELTA_O(3,I)
     &       = RI(3,1)*DELTA_C(1,I)
     &       + RI(3,2)*DELTA_C(2,I)
     &       + RI(3,3)*DELTA_C(3,I)
        

      ENDDO ! Loop for edges ends here

C Set the variable vector at the midpoint 1. This is used later to judge
C whether the variables should be limited.
      WM_O(1) = ZM(1)
      WM_O(2) = QXM(1)
      WM_O(3) = QYM(1)

C Find the possibly limited DELTA and variable vectors at the corner nodes
      DO I = 1,3 ! Loop for variable components
C Find the possibly limited DELTA
        POS = 0.D0
        NEG = 0.D0
        DO J = 1,3 ! Loop for edges
          POS = POS + MAX(0.D0, DELTA_O(I,J))
          NEG = NEG + MAX(0.D0,-DELTA_O(I,J))
        ENDDO
        IF((POS < ZERO).OR.(NEG < ZERO)) THEN 
          DELTA_HAT(I,1) = 0.D0
          DELTA_HAT(I,2) = 0.D0
          DELTA_HAT(I,3) = 0.D0
        ELSE
          THETA_P = MIN(1.D0,NEG/POS)
          THETA_N = MIN(1.D0,POS/NEG)
          
          DELTA_HAT(I,1)
     &         = THETA_P*MAX(0.D0, DELTA_O(I,1))
     &         - THETA_N*MAX(0.D0,-DELTA_O(I,1))
          DELTA_HAT(I,2)
     &         = THETA_P*MAX(0.D0, DELTA_O(I,2))
     &         - THETA_N*MAX(0.D0,-DELTA_O(I,2))
          DELTA_HAT(I,3)
     &         = THETA_P*MAX(0.D0, DELTA_O(I,3))
     &         - THETA_N*MAX(0.D0,-DELTA_O(I,3))
        ENDIF

C Compute the possibly limited variable vectors at the corner nodes
        IF(( DELTA_HAT(I,1) < (WM_O(I) - WB_O(I,0) - ZERO) ).OR.
     &       ( DELTA_HAT(I,1) > (WM_O(I) - WB_O(I,0) + ZERO) )) THEN
C Update the variables only if they need to be limited.
C If the condition in the if statement is false, the following 
C isn't applied and therefore the quradratic and higher order 
C coefficients are preserved.

C Computing limited variables at the corner nodes
          WC_O(I,1) = WB_O(I,0)
     &         - DELTA_HAT(I,1)
     &         + DELTA_HAT(I,2)
     &         + DELTA_HAT(I,3)
          WC_O(I,2) = WB_O(I,0)
     &         + DELTA_HAT(I,1)
     &         - DELTA_HAT(I,2)
     &         + DELTA_HAT(I,3)
          WC_O(I,3) = WB_O(I,0)
     &         + DELTA_HAT(I,1)
     &         + DELTA_HAT(I,2)
     &         - DELTA_HAT(I,3)
          LIMITIT(I) = .TRUE.
        ELSE
          LIMITIT(I) = .FALSE.
        ENDIF
      ENDDO

      IF(LIMITIT(1)) THEN
        ZE(2,L,IRK+1) = -1.D0/6.D0*(WC_O(1,1) + WC_O(1,2))
     &                + 1.D0/3.D0*WC_O(1,3)
        ZE(3,L,IRK+1) = -0.5D0*WC_O(1,1) + 0.5D0*WC_O(1,2)
      ENDIF
      IF(LIMITIT(2)) THEN
        QX(2,L,IRK+1) = -1.D0/6.D0*(WC_O(2,1) + WC_O(2,2))
     &                + 1.D0/3.D0*WC_O(2,3)
        QX(3,L,IRK+1) = -0.5D0*WC_O(2,1) + 0.5D0*WC_O(2,2)
      ENDIF
      IF(LIMITIT(3)) THEN
        QY(2,L,IRK+1) = -1.D0/6.D0*(WC_O(3,1) + WC_O(3,2))
     &                + 1.D0/3.D0*WC_O(3,3)
        QY(3,L,IRK+1) = -0.5D0*WC_O(3,1) + 0.5D0*WC_O(3,2)
      ENDIF

1000  CONTINUE

      RETURN
      END SUBROUTINE

***********************************************************************
C
C     SUBROUTINE HYDRO_EIGEN_VALUES
C
C     This subroutine computes the eigen values and vectors of hydro
C     dynamics equations.
C
C     Borrowed from ROE_FLUX() - 06 Feb 2008, S.B.
C
C***********************************************************************
      SUBROUTINE HYDRO_EIGEN_VALUES(H, U, V, NX, NY, EIGVAL, RI, LE)
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : IFNLCT, G

      IMPLICIT NONE

      REAL(SZ) H,U,V,NX,NY,C,DTM
      REAL(SZ) EIGVAL(3) ! eigen values
      REAL(SZ) RI(3,3)   ! columns are the right eigen vectors
      REAL(SZ) LE(3,3)   ! columns are the left eigen vectors

      C = SQRT(G*H)

C.....Evaluate the eigenvalues at the Roe averaged variables

      EIGVAL(2) = (U*NX + V*NY)*IFNLCT
      EIGVAL(1) = EIGVAL(2) + C
      EIGVAL(3) = EIGVAL(2) - C
      
C.....Evaluate right eigenvectors at Roe averaged variables

      RI(1,1) = 1.D0
      RI(2,1) = U*IFNLCT + (IFNLCT*C + (1-IFNLCT)*G/C)*NX
      RI(3,1) = V*IFNLCT + (IFNLCT*C + (1-IFNLCT)*G/C)*NY

      RI(1,2) = 0.D0
      RI(2,2) = -NY
      RI(3,2) =  NX

      RI(1,3) = 1.D0
      RI(2,3) = U*IFNLCT - (IFNLCT*C + (1-IFNLCT)*G/C)*NX
      RI(3,3) = V*IFNLCT - (IFNLCT*C + (1-IFNLCT)*G/C)*NY

C.....Evaluate left eigenvectors at Roe averaged variables

C.....Determinant
      DTM = 1.0d0/(RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2)+RI(2,1)*RI(3,2)
     &                                                 -RI(3,1)*RI(2,2))

      LE(1,1) =  DTM*( RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2) )
      LE(1,2) =  DTM*RI(3,2)
      LE(1,3) = -DTM*RI(2,2)

      LE(2,1) = -DTM*( RI(2,1)*RI(3,3)-RI(2,3)*RI(3,1) )
      LE(2,2) =  DTM*( RI(3,3)-RI(3,1) )
      LE(2,3) = -DTM*( RI(2,3)-RI(2,1) )

      LE(3,1) =  DTM*( RI(2,1)*RI(3,2) - RI(2,2)*RI(3,1) )
      LE(3,2) = -DTM*RI(3,2)
      LE(3,3) =  DTM*RI(2,2)

      RETURN
      END SUBROUTINE

c$$$C invert a 3x3 matrix
c$$$      subroutine invert(a,ainv,succeed)
c$$$      sizes,only : sz
c$$$      implicit none
c$$$      real(sz),intent(in) ::  a(3,3)
c$$$      real(sz),intent(out) ::  ainv(3,3)
c$$$      logical,intent(out) :: succeed
c$$$      integer :: i, j, k
c$$$      real(sz) :: det
c$$$      real(sz) :: e(3,3)
c$$$
c$$$      det = 
c$$$     &     a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) +
c$$$     &     a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) +
c$$$     &     a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
c$$$
c$$$      if(det.eq.0.D0) then
c$$$        succeed = .false.
c$$$        return
c$$$      endif
c$$$
c$$$      Ainv(1,1) = (A(1,1)*A(2,2) - A(1,2)*A(2,1))/det
c$$$      Ainv(1,2) = (A(0,2)*A(2,1) - A(0,1)*A(2,2))/det
c$$$      Ainv(1,3) = (A(0,1)*A(1,2) - A(0,2)*A(1,1))/det
c$$$
c$$$      Ainv(1,0) = (A(1,2)*A(2,0) - A(1,0)*A(2,2))/det
c$$$      Ainv(1,1) = (A(0,0)*A(2,2) - A(0,2)*A(2,0))/det
c$$$      Ainv(1,2) = (A(0,2)*A(1,0) - A(0,0)*A(1,2))/det
c$$$  
c$$$      Ainv(2,0) = (A(1,0)*A(2,1) - A(1,1)*A(2,0))/det
c$$$      Ainv(2,1) = (A(0,1)*A(2,0) - A(0,0)*A(2,1))/det
c$$$      Ainv(2,2) = (A(0,0)*A(1,1) - A(0,1)*A(1,0))/det
c$$$
c$$$      do i = 1,3
c$$$        do j = 1,3
c$$$          E(i,j) = 0.0
c$$$          do k = 1,3
c$$$            E(i,j) += A(i,k)*Ainv(k,j)
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$      end subroutine
c$$$
c$$$/* compute projected velocity u */
c$$$void compute_projection_u
c$$$(double *H, double *q, 
c$$$ int num_pts, double *area_coords[3], double *weights,
c$$$ double *u)
c$$${
c$$$  double A[3][3], Ainv[3][3];
c$$$  double b[3];
c$$$
c$$$  int i, j, k, p;
c$$$
c$$$  /* initialize */
c$$$  for(j = 0; j < 3; j++) {
c$$$    for(i = 0; i < 3; i++) {
c$$$      A[j][i] = 0.0;
c$$$    }
c$$$    b[j] = 0.0;
c$$$  }
c$$$
c$$$  /* determine matrix A */
c$$$  for(j = 0; j < 3; j++) {
c$$$    for(i = 0; i < 3; i++) {
c$$$      for(k = 0; k < 3; k++) {
c$$$	for(p = 0; p < num_pts; p++) {
c$$$	  A[j][i] += 
c$$$	    H[k]*
c$$$	    area_coords[i][p]*
c$$$	    area_coords[j][p]*
c$$$	    area_coords[k][p]*
c$$$	    weights[p];
c$$$	}
c$$$      }
c$$$    }
c$$$  }
c$$$
c$$$  /* determine vector b */
c$$$  for(j = 0; j < 3; j++) {
c$$$    for(i = 0; i < 3; i++) {
c$$$      for(p = 0; p < num_pts; p++) {
c$$$	b[j] += 
c$$$	  q[i]*
c$$$	  area_coords[i][p]*
c$$$	  area_coords[j][p]*
c$$$	  weights[p];
c$$$      }
c$$$    }
c$$$  }
c$$$
c$$$  /* obtain the inverse of A */
c$$$  invert(A, Ainv);
c$$$
c$$$  /* obtain projection u */
c$$$  for(j = 0; j < 3; j++) {
c$$$    u[j] = 0;
c$$$    for(i = 0; i < 3; i++) {
c$$$      u[j] += Ainv[j][i]*b[i];
c$$$    }
c$$$  }
c$$$}

C***********************************************************************
C
C     SUBROUTINE SLOPELIMITER5()
C
C     Written by Clint Dawson - 30 June 2010
C
C
C***********************************************************************

      SUBROUTINE SLOPELIMITER5(lvl)

C.....Use appropriate modules

      USE SIZES, ONLY : SZ,myproc
      USE GLOBAL
      USE DG

#ifdef CMPI
      USE MESSENGER_node_new
#endif

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LL, INC1,INC2,INC3,KDP,NN,IVAR,lvl
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF,
     $     DIV,REDFAC,REDMAX,ZEMAX(3),ZEMIN(3),HMIN
      real*4 tarray1(2),tarray2(2)

C  FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS 
C  SHARING A NODE
      DO N1=1,NP
         VISIT_FLAG(N1)=.FALSE.
      ENDDO
      DO J=1,Nelem(lvl)
! CND
         I=ElemNum(lvl,J)
c         IF(WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
         DO KK=1,3
            N1=NM(I,KK)
C   IF THIS NODE HAS ALREADY BEEN DONE SKIP IT
            IF (VISIT_FLAG(N1).EQ..TRUE.) GO TO 1590
            VISIT_FLAG(N1)=.TRUE.

c     
            ZE_MIN(N1)=99999.
            ZE_MAX(N1)=-99999.
            QX_MIN(N1)=99999.
            QX_MAX(N1)=-99999.
            QY_MIN(N1)=99999.
            QY_MAX(N1)=-99999.
            NO_NBORS = EL_COUNT(N1)

            DO JJ = 1,NO_NBORS
               NBOR_EL = ELETAB(N1,1+JJ)

c     IF(WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
               if (tsgroup(nbor_el).eq.lvl) then
                  ZE_DG(JJ) = ZE(1,NBOR_EL,IRK+1)
                  QX_DG(JJ) = QX(1,NBOR_EL,IRK+1)
                  QY_DG(JJ) = QY(1,NBOR_EL,IRK+1)
               else
                  ZE_DG(JJ) = ZE(1,NBOR_EL,1)
                  QX_DG(JJ) = QX(1,NBOR_EL,1)
                  QY_DG(JJ) = QY(1,NBOR_EL,1)
               endif
C     
               IF (ZE_DG(JJ).LT.ZE_MIN(N1))THEN
                  ZE_MIN(N1)=ZE_DG(JJ)
               ENDIF
               IF (ZE_DG(JJ).GT.ZE_MAX(N1)) THEN
                  ZE_MAX(N1)=ZE_DG(JJ)
               ENDIF
               IF (QX_DG(JJ).LT.QX_MIN(N1))THEN
                  QX_MIN(N1)=QX_DG(JJ)
               ENDIF
               IF (QX_DG(JJ).GT.QX_MAX(N1)) THEN
                  QX_MAX(N1)=QX_DG(JJ)
               ENDIF
               IF (QY_DG(JJ).LT.QY_MIN(N1))THEN
                  QY_MIN(N1)=QY_DG(JJ)
               ENDIF
               IF (QY_DG(JJ).GT.QY_MAX(N1)) THEN
                  QY_MAX(N1)=QY_DG(JJ)
               ENDIF
C  END JJ LOOP
            ENDDO
 1590       CONTINUE
C  END KK LOOP            
         ENDDO
C  END J LOOP
      ENDDO
c      if (myproc.eq.0) write(*,*) 'in slopelimiter',lvl
#ifdef CMPI
c      CALL UPDATER(ZE_MIN,ZE_MAX,QX_MIN,3)
c      CALL UPDATER(QX_MAX,QY_MIN,QY_MAX,3)
      call etime(tarray1)
      call updater_node_mod(ze_min,lvl)
c      write(*,*) myproc,'in slopelimiter after updater 1'
      call updater_node_mod(ze_max,lvl)
c      write(*,*) myproc,'in slopelimiter after updater 2'
      call updater_node_mod(qx_min,lvl)
c      write(*,*) myproc,'in slopelimiter after updater 3'
      call updater_node_mod(qx_max,lvl)
c      write(*,*) myproc,'in slopelimiter after updater 4'
      call updater_node_mod(qy_min,lvl)
c      write(*,*) myproc,'in slopelimiter after updater 5'
      call updater_node_mod(qy_max,lvl)
      call etime(tarray2)
      TIMINGS(4)=TIMINGS(4)
     $     +tarray2(1)+tarray2(2)-(tarray1(1)+tarray1(2))

#endif

C
C  LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
C
      DO J=1,Nelem(lvl)
! CND
         I=ElemNum(lvl,J)
c         IF(WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
         N1=NM(I,1)
         N2=NM(I,2)
         N3=NM(I,3)
         DO IVAR=1,3
            IF (IVAR.EQ.1) THEN
               ZEC(1)=ZE(1,I,IRK+1)
               ZEC(2)=ZE(2,I,IRK+1)
               ZEC(3)=ZE(3,I,IRK+1)
               ZEMAX(1)=ZE_MAX(N1)
               ZEMIN(1)=ZE_MIN(N1)
               ZEMAX(2)=ZE_MAX(N2)
               ZEMIN(2)=ZE_MIN(N2)
               ZEMAX(3)=ZE_MAX(N3)
               ZEMIN(3)=ZE_MIN(N3)
            ENDIF
            IF (IVAR.EQ.2) THEN
               ZEC(1)=QX(1,I,IRK+1)
               ZEC(2)=QX(2,I,IRK+1)
               ZEC(3)=QX(3,I,IRK+1)
               ZEMAX(1)=QX_MAX(N1)
               ZEMIN(1)=QX_MIN(N1)
               ZEMAX(2)=QX_MAX(N2)
               ZEMIN(2)=QX_MIN(N2)
               ZEMAX(3)=QX_MAX(N3)
               ZEMIN(3)=QX_MIN(N3)
            ENDIF

            IF (IVAR.EQ.3) THEN
               ZEC(1)=QY(1,I,IRK+1)
               ZEC(2)=QY(2,I,IRK+1)
               ZEC(3)=QY(3,I,IRK+1)
               ZEMAX(1)=QY_MAX(N1)
               ZEMIN(1)=QY_MIN(N1)
               ZEMAX(2)=QY_MAX(N2)
               ZEMIN(2)=QY_MIN(N2)
               ZEMAX(3)=QY_MAX(N3)
               ZEMIN(3)=QY_MIN(N3)
            ENDIF
C  COMPUTE THE VERTEX VALUES
            ZEVERTEX(1)=ZEC(1)
            ZEVERTEX(2)=ZEC(1)
            ZEVERTEX(3)=ZEC(1)
            DO KK=2,DOF
               ZEVERTEX(1)=ZEVERTEX(1)+ PHI_CORNER(KK,1)*ZEC(KK)
               ZEVERTEX(2)=ZEVERTEX(2)+ PHI_CORNER(KK,2)*ZEC(KK)
               ZEVERTEX(3)=ZEVERTEX(3)+ PHI_CORNER(KK,3)*ZEC(KK)
            ENDDO
c            IF (IVAR.EQ.1) THEN
c               HMIN = MAX(0.D0,
c     $              MIN(ZEVERTEX(1),ZEVERTEX(2),ZEVERTEX(3)))
c               IF (NOLIFA.EQ.2.AND.HMIN.LT.H0) RETURN
c            ENDIF
C
C  RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
C  GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
C
            ZEVERTEX(1)=DMAX1(DMIN1(ZEVERTEX(1),ZEMAX(1)),ZEMIN(1))
            ZEVERTEX(2)=DMAX1(DMIN1(ZEVERTEX(2),ZEMAX(2)),ZEMIN(2))
            ZEVERTEX(3)=DMAX1(DMIN1(ZEVERTEX(3),ZEMAX(3)),ZEMIN(3))
C  LOOP OVER THE VERTICES 3 TIMES
C  IF THE VALUE AT THE VERTEX IS ABOVE (BELOW) THE MAX (MIN) AT THAT
C  VERTEX THEN SUBTRACT OFF THE DIFFERENCE AND ADD IT TO THE OTHER
C  VERTICES
C  
            DO LL=1,3
               SUMLOC=(ZEVERTEX(1)+ZEVERTEX(2)+ZEVERTEX(3))/3.0D0
               SUMDIF=(SUMLOC-ZEC(1))*3.0D0
               SIGNDIF=DSIGN(1.D0,SUMDIF)
               DIF(1)=(ZEVERTEX(1)-ZEC(1))*SIGNDIF
               DIF(2)=(ZEVERTEX(2)-ZEC(1))*SIGNDIF
               DIF(3)=(ZEVERTEX(3)-ZEC(1))*SIGNDIF
               INC1=0
               IF (DIF(1).GT.0) INC1=1
               INC2=0
               IF (DIF(2).GT.0) INC2=1
               INC3=0
               IF (DIF(3).GT.0) INC3=1
               KDP=INC1+INC2+INC3
C     
               DO K=1,3
                  DIV=DMAX1(1.D0,DFLOAT(KDP))
                  IF (DIF(K).GT.0) THEN
                     REDFAC=SUMDIF*SIGNDIF/DIV
                     KDP=KDP-1
                  ELSE
                     REDFAC=0
                  ENDIF
                  IF (SIGNDIF.GT.0) THEN
                     REDMAX=ZEVERTEX(K)-ZEMIN(K)
                  ELSE
                     REDMAX=ZEMAX(K)-ZEVERTEX(K)
                  ENDIF
                  REDFAC=DMIN1(REDFAC,REDMAX)
                  SUMDIF=SUMDIF-REDFAC*SIGNDIF
                  ZEVERTEX(K)=ZEVERTEX(K)-REDFAC*SIGNDIF
               ENDDO
            ENDDO
            IF (IVAR.EQ.1) THEN
               ZE(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               ZE(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.2) THEN
               QX(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               QX(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.3) THEN
               QY(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               QY(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE 
