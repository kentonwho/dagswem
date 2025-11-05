C***********************************************************************
C
C     SUBROUTINE OCEAN_EDGE_HYDRO( )
C
C     This subroutine does the following:
C
C       1.  Calculates the values of the necessary variables at the edge
C           gauss points for ELEVATION SPECIFIED edges
C       2.  Calls the appropriate subroutine to compute the flux at
C           these points.
C       3.  Calls the appropriate subroutine to compute the boundary
C           integrals.
C
C     Written by Ethan Kubatko (06-11-2004)
C
C***********************************************************************
      SUBROUTINE OCEAN_EDGE_HYDRO(lvl)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED
c
      USE NodalAttributes, ONLY: GeoidOffset, LoadGeoidOffset

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED, lvl
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN
      REAL(SZ) UFA_GP, UMO_GP, VFA_GP, VMO_GP, W_IN

      DO 1000 L = 1, needs
      
C.....Retrieve the global and local edge number

        GED = NEEDN(L)
        LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge

        EL_IN = NEDEL(1,GED)
!CND
        if (tsgroup(el_in).ne.lvl) go to 1000
        
C.....If the element is dry then skip the edge calculation

        IF (WDFLG(EL_IN).EQ.0) GOTO 1000
        
C.....Retrieve the components of the normal vector to the edge
        NX = COSNX(GED)
        NY = SINNX(GED)
        
C.....Retrieve the nodes of the edge
      
        N1 = NEDNO(1,GED)
        N2 = NEDNO(2,GED)
      
C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

        DO I = 1,NEGP

          ZE_IN = ZE(1,EL_IN,IRK)
          QX_IN = QX(1,EL_IN,IRK)
          QY_IN = QY(1,EL_IN,IRK)

          ZE_EX = 0.D0
          QX_EX = 0.D0
          QY_EX = 0.D0
          U_EX  = 0.D0
          V_EX  = 0.D0
        

          HB_IN = BATHED(I,LED,EL_IN)
cnd
          SFAC_IN = SFACED(I,LED,EL_IN)
c
        
          LZ_XX_IN = LZ(1,1,1,EL_IN)
          LZ_XY_IN = LZ(1,1,2,EL_IN)
          LZ_YX_IN = LZ(1,2,1,EL_IN)
          LZ_YY_IN = LZ(1,2,2,EL_IN)

C.....Compute the specified open ocean elevation
        
          DO JJ=1,NBFR
        
            IF (PER(JJ).EQ.0.D0) THEN
              NCYC = 0.D0
            ELSE
              NCYC = INT(TIMEDG/PER(JJ))
            ENDIF
            
C...........Surface Elevation

            ARGJ = AMIG(JJ)*(TIMEDG - NCYC*PER(JJ)) + FACE(JJ)
            RFF = FF(JJ)*RAMPDG
          
            EFA_GP = 0.5D0*(EFA_DG(JJ,L,1) + EFA_DG(JJ,L,2))
     &             + 0.5D0*(EFA_DG(JJ,L,2) - EFA_DG(JJ,L,1))*XEGP(I)
            EMO_GP = 0.5D0*(EMO_DG(JJ,L,1) + EMO_DG(JJ,L,2))
     &             + 0.5D0*(EMO_DG(JJ,L,2) - EMO_DG(JJ,L,1))*XEGP(I)

            ARG = ARGJ - EFA_GP
          
            ZE_EX = ZE_EX + EMO_GP*RFF*COS(ARG)

C...........X Velocity
            
C            UFA_GP = 0.5D0*(UFA_DG(JJ,L,1) + UFA_DG(JJ,L,2))
C     &             + 0.5D0*(UFA_DG(JJ,L,2) - UFA_DG(JJ,L,1))*XEGP(I)
C            UMO_GP = 0.5D0*(UMO_DG(JJ,L,1) + UMO_DG(JJ,L,2))
C     &             + 0.5D0*(UMO_DG(JJ,L,2) - UMO_DG(JJ,L,1))*XEGP(I)
C
C            ARG = ARGJ - UFA_GP
C
C            U_EX = U_EX + UMO_GP*RFF*COS(ARG)
            
C...........Y Velocity

C            VFA_GP = 0.5D0*(VFA_DG(JJ,L,1) + VFA_DG(JJ,L,2))
C     &             + 0.5D0*(VFA_DG(JJ,L,2) - VFA_DG(JJ,L,1))*XEGP(I)
C            VMO_GP = 0.5D0*(VMO_DG(JJ,L,1) + VMO_DG(JJ,L,2))
C     &             + 0.5D0*(VMO_DG(JJ,L,2) - VMO_DG(JJ,L,1))*XEGP(I)
C
C            ARG = ARGJ - VFA_GP
C
C            V_EX = V_EX + VMO_GP*RFF*COS(ARG)

          ENDDO


C.....Compute the solution at the interior state

          DO K = 2,DOF

            ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)
            QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)
            QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)

            ! LDG terms
            LZ_XX_IN = LZ_XX_IN + LZ(K,1,1,EL_IN)*PHI_EDGE(K,I,LED)
            LZ_XY_IN = LZ_XY_IN + LZ(K,1,2,EL_IN)*PHI_EDGE(K,I,LED)
            LZ_YX_IN = LZ_YX_IN + LZ(K,2,1,EL_IN)*PHI_EDGE(K,I,LED)
            LZ_YY_IN = LZ_YY_IN + LZ(K,2,2,EL_IN)*PHI_EDGE(K,I,LED)

          ENDDO

C.....Set the exterior value of the bathymetry equal to the interior

          HB_EX = HB_IN
          SFAC_EX = SFAC_IN
          IF (LoadGeoidOffset) then
             ZE_EX = ZE_EX + .5*(GeoidOffset(N1)+GeoidOffset(N2))
          endif


C.....Set the exterior state flows equal to the interior state flows

           QX_EX = QX_IN
           QY_EX = QY_IN
C          QX_EX = U_EX*(ZE_EX + HB_EX)
C          QY_EX = V_EX*(ZE_EX + HB_EX)
c           write(*,*) 'ze_in,ze_ex,qx_in,qx_ex',ze_in,ze_ex,qx_in,qx_ex
C.....Compute the Roe flux

          CALL NUMERICAL_FLUX()
        
C.....Add LDG terms
          
          G_HAT = G_HAT + LZ_XX_IN*NX*SFAC_IN + LZ_XY_IN*NY
          H_HAT = H_HAT + LZ_YX_IN*NX*SFAC_IN + LZ_YY_IN*NY

C.....Compute the edge integral
          DO K = 1,DOF

            W_IN = 2.0*M_INV(K)/AREAS(EL_IN)*XLEN(GED)*
     &             PHI_EDGE(K,I,LED)*WEGP(I)

            RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT
            RHS_QX(K,EL_IN,IRK) = RHS_QX(K,EL_IN,IRK) - W_IN*G_HAT
            RHS_QY(K,EL_IN,IRK) = RHS_QY(K,EL_IN,IRK) - W_IN*H_HAT

          ENDDO
        ENDDO

1000  CONTINUE
      
      RETURN
      END SUBROUTINE
