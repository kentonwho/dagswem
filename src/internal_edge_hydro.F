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
C-----------------------------------------------------------------------
C
C     01-02-2007, sb, Modified for LDG
C
C***********************************************************************

      SUBROUTINE INTERNAL_EDGE_HYDRO(lvl)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED
      USE NodalAttributes, ONLY :  ESLM
      use sizes, only:  sz, myproc

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED_IN, LED_EX, GED, GP_IN, GP_EX
      INTEGER IRK_IN, IRK_EX, lvl, EL_SMALL, EL_BIG
      REAL(SZ), PARAMETER :: ZERO = 1.D-12
      REAL(SZ) TX, TY, W_IN, W_EX
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN
      REAL(SZ) LZ_XX_EX, LZ_XY_EX, LZ_YX_EX, LZ_YY_EX
      REAL(SZ) EDFAC_IN, EDFAC_EX, DEN
      REAL(SZ) XLEN_EL_IN, XLEN_EL_EX
      REAL(SZ) MASS_EL_IN, MASS_EL_EX
      REAL(SZ), SAVE, ALLOCATABLE :: 
     &     RHS_ZE_IN(:), RHS_QX_IN(:), RHS_QY_IN(:),
     &     RHS_ZE_EX(:), RHS_QX_EX(:), RHS_QY_EX(:)
      REAL(SZ) ARK, BRK
      REAL(SZ) MAX_BOA ! Maximum of beta_il/alpha_il for all l
      REAL(SZ) NLEQG_TMP, G_TMP
      REAL(SZ) F_HAT_O, G_HAT_O, H_HAT_O
      REAL(SZ) G_HAT_IN, H_HAT_IN
      REAL(SZ) G_HAT_EX, H_HAT_EX
      REAL(SZ) DT_Ratio_EX, DT_Ratio_IN ! CJT

      IF(.NOT.ALLOCATED(RHS_ZE_IN)) THEN
        ALLOCATE ( RHS_ZE_IN(DOF),RHS_QX_IN(DOF),RHS_QY_IN(DOF) )
        ALLOCATE ( RHS_ZE_EX(DOF),RHS_QX_EX(DOF),RHS_QY_EX(DOF) )
      ENDIF

      DO 1000 L = 1,NIEDS

C.......Retrieve the global and local edge number

        GED = NIEDN(L)
        LED_IN = NEDSD(1,GED)
        LED_EX = NEDSD(2,GED)

C.......Retrieve the elements which share the edge

!        ! Force EL_IN to be the smaller timestep
!        IF (TsGroup(NEDEL(1,GED)) < TsGroup(NEDEL(2,GED))) THEN
!          EL_IN = NEDEL(1,GED)
!          EL_EX = NEDEL(2,GED)
!        ENDIF

        EL_IN = NEDEL(1,GED)
        EL_EX = NEDEL(2,GED)

	!WRITE(669,*) L,EL_IN,TsGroup(EL_IN),EL_EX,TsGroup(EL_EX)
        !IF (TsGroup(EL_EX)<TsGroup(EL_IN)) WRITE(*,*) L,EL_IN,TsGroup(EL_IN),EL_EX,TsGroup(EL_EX) 

!K: Handle the interfaces between different heats
C.......Find which element is smaller, exit if smaller level is not = lvl
c.......Only want to sweep across this edge once ..
        IF (TsGroup(EL_IN) /= lvl .and. TsGroup(EL_EX) /= lvl) THEN
          GOTO 1000
        ELSE
          IF (TsGroup(EL_IN) == lvl) THEN 
            IF (TsGroup(EL_EX) < lvl) GOTO 1000  ! flux already computed and added to rhs
            IRK_IN = IRK
            IF (TsGroup(EL_EX) == TsGroup(EL_IN)) THEN
              IRK_EX = IRK
            ELSE
              IRK_EX = 1
            ENDIF
          ELSE
            IF (TsGroup(EL_IN) < lvl) GOTO 1000  ! flux already computed and added to rhs
            IRK_EX = IRK
            IF (TsGroup(EL_IN) == TsGroup(EL_EX)) THEN
              IRK_IN = IRK
            ELSE
              IRK_IN = 1
            ENDIF
          ENDIF
        ENDIF
        IF (TSGroup(EL_IN)==lvl) DT_Ratio_IN = 1.d0;
        IF (TsGroup(EL_EX)==lvl) DT_Ratio_EX = 1.d0;
        IF (TsGroup(EL_EX)>lvl)  DT_Ratio_EX=DtDom(TsGroup(EL_IN))/DtDom(TsGroup(EL_EX));
        IF (TsGroup(EL_IN)>lvl)  DT_Ratio_IN=DtDom(TsGroup(EL_EX))/DtDom(TsGroup(EL_IN));


C.......If both elements on either side of edge are dry then skip
c        if (myproc.eq.22.and.el_in.eq.8801) then
c           write(5555,*) 'in internal edge in ',el_ex,wdflg(el_in),
c     $          wdflg(el_ex),tsgroup(el_ex)
c        endif
c        if (myproc.eq.22.and.el_ex.eq.8801) then
c           write(5555,*) 'in internal edge ex ',el_in,wdflg(el_ex),
c     $          wdflg(el_in),tsgroup(el_in)
c        endif
        IF((WDFLG(EL_IN).EQ.0).AND.(WDFLG(EL_EX).EQ.0)) GOTO 1000

C.....Save the current RHS

        DO K = 1,DOF
          RHS_ZE_IN(K) = RHS_ZE(K,EL_IN,IRK_IN);
          RHS_QX_IN(K) = RHS_QX(K,EL_IN,IRK_IN);
          RHS_QY_IN(K) = RHS_QY(K,EL_IN,IRK_IN);
          RHS_ZE_EX(K) = RHS_ZE(K,EL_EX,IRK_EX);
          RHS_QX_EX(K) = RHS_QX(K,EL_EX,IRK_EX);
          RHS_QY_EX(K) = RHS_QY(K,EL_EX,IRK_EX);
        ENDDO

C.....Compute the sum of the lengths of three edges

        XLEN_EL_IN = XLEN(NELED(1,EL_IN))
        XLEN_EL_IN = XLEN_EL_IN + XLEN(NELED(2,EL_IN))
        XLEN_EL_IN = XLEN_EL_IN + XLEN(NELED(3,EL_IN))

        XLEN_EL_EX = XLEN(NELED(1,EL_EX))
        XLEN_EL_EX = XLEN_EL_EX + XLEN(NELED(2,EL_EX))
        XLEN_EL_EX = XLEN_EL_EX + XLEN(NELED(3,EL_EX))

C.....Compute the total mass in the elements

        MASS_EL_IN = (ZE(1,EL_IN,IRK_IN)+HB(1,EL_IN,1))*AREAS(EL_IN)*0.5D0
        MASS_EL_EX = (ZE(1,EL_EX,IRK_EX)+HB(1,EL_EX,1))*AREAS(EL_EX)*0.5D0

C.....Retrieve the components of the normal vector to the edge
c
        NX = COSNX(GED)
        NY = SINNX(GED)
        
C.....Set the components for the tangential vector to the edge

        TX = -NY
        TY =  NX

        
        EDFAC_IN = XLEN(GED)/AREAS(EL_IN)
        EDFAC_EX = XLEN(GED)/AREAS(EL_EX)

C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

        DO I = 1,NEGP
          GP_IN = I
          GP_EX = NEGP - I + 1

          ZE_IN = ZE(1,EL_IN,IRK_IN)
          QX_IN = QX(1,EL_IN,IRK_IN)
          QY_IN = QY(1,EL_IN,IRK_IN)
          HB_IN = BATHED(GP_IN,LED_IN,EL_IN)
cnd
          SFAC_IN = SFACED(GP_IN,LED_IN,EL_IN)

          ZE_EX = ZE(1,EL_EX,IRK_EX)
          QX_EX = QX(1,EL_EX,IRK_EX)
          QY_EX = QY(1,EL_EX,IRK_EX)
          HB_EX = BATHED(GP_EX,LED_EX,EL_EX)

cnd
          SFAC_EX = SFACED(GP_EX,LED_EX,EL_EX)

          LZ_XX_IN = LZ(1,1,1,EL_IN)
          LZ_XY_IN = LZ(1,1,2,EL_IN)
          LZ_YX_IN = LZ(1,2,1,EL_IN)
          LZ_YY_IN = LZ(1,2,2,EL_IN)

          LZ_XX_EX = LZ(1,1,1,EL_EX)
          LZ_XY_EX = LZ(1,1,2,EL_EX)
          LZ_YX_EX = LZ(1,2,1,EL_EX)
          LZ_YY_EX = LZ(1,2,2,EL_EX)

          DO K = 2,DOF

            ZE_IN = ZE_IN + ZE(K,EL_IN,IRK_IN)*PHI_EDGE(K,GP_IN,LED_IN)
            QX_IN = QX_IN + QX(K,EL_IN,IRK_IN)*PHI_EDGE(K,GP_IN,LED_IN)
            QY_IN = QY_IN + QY(K,EL_IN,IRK_IN)*PHI_EDGE(K,GP_IN,LED_IN)

            ZE_EX = ZE_EX + ZE(K,EL_EX,IRK_EX)*PHI_EDGE(K,GP_EX,LED_EX)
            QX_EX = QX_EX + QX(K,EL_EX,IRK_EX)*PHI_EDGE(K,GP_EX,LED_EX)
            QY_EX = QY_EX + QY(K,EL_EX,IRK_EX)*PHI_EDGE(K,GP_EX,LED_EX)

            ! LDG terms
            LZ_XX_IN = LZ_XX_IN 
     &               + LZ(K,1,1,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN)
            LZ_XY_IN = LZ_XY_IN 
     &               + LZ(K,1,2,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN)
            LZ_YX_IN = LZ_YX_IN 
     &               + LZ(K,2,1,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN)
            LZ_YY_IN = LZ_YY_IN 
     &               + LZ(K,2,2,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN)

            LZ_XX_EX = LZ_XX_EX 
     &               + LZ(K,1,1,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX)
            LZ_XY_EX = LZ_XY_EX 
     &               + LZ(K,1,2,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX)
            LZ_YX_EX = LZ_YX_EX 
     &               + LZ(K,2,1,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX)
            LZ_YY_EX = LZ_YY_EX 
     &               + LZ(K,2,2,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX)

          ENDDO

C.....Compute the numerical flux

c          iwrite=0
c          if (myproc.eq.41.and.el_ex.eq.13987) iwrite=1
c          if (myproc.eq.41.and.el_in.eq.13987) iwrite=1
c          if (myproc.eq.41.and.el_ex.eq.12362) iwrite=1
c          if (myproc.eq.41.and.el_in.eq.12362) iwrite=1

          CALL NUMERICAL_FLUX()
c          iwrite=0
          F_HAT_O  = F_HAT
          G_HAT_IN = G_HAT
          H_HAT_IN = H_HAT
          G_HAT_EX = G_HAT
          H_HAT_EX = H_HAT
c
c          if (myproc.eq.22.and.el_in.eq.8801) then
c             write(5555,*) 'el_in ',time_a,f_hat,ze_in,ze_ex,qx_in,qx_ex,
c     $            qy_in,qy_ex
c          endif
c          if (myproc.eq.22.and.el_ex.eq.8801) then
c             write(5555,*) 'el_ex ',time_a,f_hat,ze_in,ze_ex,qx_in,qx_ex,
c     $            qy_in,qy_ex
c          endif
C.....Check if the flux is large enough to dry up the elements
C.....1.01D0 is a safty factor.

          IF ( (1.01*F_HAT*XLEN_EL_IN*MAX_BOA_DT(IRK_IN).GE.MASS_EL_IN).OR.
     &         (1.01*F_HAT*XLEN_EL_EX*MAX_BOA_DT(IRK_EX)*(-1.D0).GE.
     &                                                 MASS_EL_EX)) THEN


C...........Put back the saved RHS

            DO K = 1,DOF
              RHS_ZE(K,EL_IN,IRK_IN) = RHS_ZE_IN(K)
              RHS_QX(K,EL_IN,IRK_IN) = RHS_QX_IN(K)
              RHS_QY(K,EL_IN,IRK_IN) = RHS_QY_IN(K)

              RHS_ZE(K,EL_EX,IRK_EX) = RHS_ZE_EX(K)
              RHS_QX(K,EL_EX,IRK_EX) = RHS_QX_EX(K)
              RHS_QY(K,EL_EX,IRK_EX) = RHS_QY_EX(K)
            ENDDO
cnd try something
            if (tsgroup(el_in).ne.1) then
               TSGROUPFLAG(EL_IN)=.true.
               CALLRESETFLAG1=1
c               write(*,*) myproc,el_in
            endif
            if (tsgroup(el_ex).ne.1) then
               TSGROUPFLAG(EL_EX)=.true.
               CALLRESETFLAG1=1
c               write(*,*) myproc,el_ex
            endif

cnd 
c          if (myproc.eq.22) then
c            write(5556,*) 'putting back rhs',myproc,el_in,el_ex
c          endif
            GOTO 100
          ENDIF
          
C........Check to make sure mass flux is not coming from a dry element
          IF (abs(F_HAT).gt.1.d-12) THEN
          IF (WDFLG(EL_IN).EQ.0) THEN
            IF (F_HAT.GT.0) THEN
              Q_N_EXT = QX_EX*NX + QY_EX*NY
              Q_T_EXT = QX_EX*TX + QY_EX*TY
              Q_N_INT = -Q_N_EXT
              Q_T_INT =  Q_T_EXT
              DEN = 1.D0/(NX*TY - NY*TX)
              QX_IN = ( TY*Q_N_INT - NY*Q_T_INT)*DEN
              QY_IN = (-TX*Q_N_INT + NX*Q_T_INT)*DEN
              ZE_IN = ZE_EX
              HB_IN = HB_EX
              SFAC_IN = SFAC_EX
              CALL NUMERICAL_FLUX()
              F_HAT_O  = F_HAT
              G_HAT_IN = G_HAT
              H_HAT_IN = H_HAT
              G_HAT_EX = G_HAT
              H_HAT_EX = H_HAT
            ELSE
              NLEQG_TMP = NLEQG
              NLEQG = 0.D0
              G_TMP = G
              G = 0.D0
              CALL NUMERICAL_FLUX()
              NLEQG = NLEQG_TMP
              G = G_TMP
              G_HAT_IN = G_HAT
              H_HAT_IN = H_HAT
            ENDIF
          ELSEIF (WDFLG(EL_EX).EQ.0) THEN
            IF (F_HAT.LT.0) THEN
              Q_N_INT = QX_IN*NX + QY_IN*NY
              Q_T_INT = QX_IN*TX + QY_IN*TY
              Q_N_EXT = -Q_N_INT
              Q_T_EXT =  Q_T_INT
              DEN = 1.D0/(NX*TY - NY*TX)
              QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)*DEN
              QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)*DEN
              ZE_EX = ZE_IN
              HB_EX = HB_IN
              SFAC_EX = SFAC_IN
              CALL NUMERICAL_FLUX()
              F_HAT_O  = F_HAT
              G_HAT_IN = G_HAT
              H_HAT_IN = H_HAT
              G_HAT_EX = G_HAT
              H_HAT_EX = H_HAT
            ELSE
              NLEQG_TMP = NLEQG
              NLEQG = 0.D0
              G_TMP = G
              G = 0.D0
              CALL NUMERICAL_FLUX()
              NLEQG = NLEQG_TMP
              G = G_TMP
              G_HAT_EX = G_HAT
              H_HAT_EX = H_HAT
            ENDIF
          ENDIF
          ENDIF

          G_HAT_IN = G_HAT_IN + 0.5*(LZ_XX_IN*SFAC_IN
     $         + LZ_XX_EX*SFAC_EX)*NX
     &         + 0.5*(LZ_XY_IN + LZ_XY_EX)*NY
cnd try adding the penalty term in the LDG method
     $         +ESLM*.01*(QX_IN-QX_EX)
          H_HAT_IN = H_HAT_IN + 0.5*(LZ_YX_IN*SFAC_IN
     $         + LZ_YX_EX*SFAC_EX)*NX
     &         + 0.5*(LZ_YY_IN + LZ_YY_EX)*NY
cnd try adding the penalty term in the LDG method
     $         +ESLM*.01*(QY_IN-QY_EX)
          G_HAT_EX = G_HAT_EX + 0.5*(LZ_XX_IN*SFAC_IN
     $         + LZ_XX_EX*SFAC_EX)*NX
     &         + 0.5*(LZ_XY_IN + LZ_XY_EX)*NY
cnd try adding the penalty term in the LDG method
     $         +ESLM*.01*(QX_IN-QX_EX)
          H_HAT_EX = H_HAT_EX + 0.5*(LZ_YX_IN*SFAC_IN
     $         + LZ_YX_EX*SFAC_EX)*NX
     &         + 0.5*(LZ_YY_IN + LZ_YY_EX)*NY
cnd try adding the penalty term in the LDG method
     $         +ESLM*.01*(QY_IN-QY_EX)
c



C.....Compute the edge integral
c          IF (TSGroup(EL_IN)==lvl) DT_Ratio_IN = 1.d0;
c          IF (TsGroup(EL_EX)==lvl) DT_Ratio_EX = 1.d0;
c          IF (TsGroup(EL_EX)>lvl)  DT_Ratio_EX=DtDom(TsGroup(EL_IN))/DtDom(TsGroup(EL_EX));
c          IF (TsGroup(EL_IN)>lvl)  DT_Ratio_IN=DtDom(TsGroup(EL_EX))/DtDom(TsGroup(EL_IN));

          DO K = 1,DOF
          
            W_IN = EDFAC_IN*EDGEQ(K,GP_IN,LED_IN)
            W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX)

            RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT_O  * DT_Ratio_IN
            RHS_QX(K,EL_IN,IRK) = RHS_QX(K,EL_IN,IRK) - W_IN*G_HAT_IN * DT_Ratio_IN
            RHS_QY(K,EL_IN,IRK) = RHS_QY(K,EL_IN,IRK) - W_IN*H_HAT_IN * DT_Ratio_IN
            RHS_ZE(K,EL_EX,IRK) = RHS_ZE(K,EL_EX,IRK) + W_EX*F_HAT_O  * DT_Ratio_EX
            RHS_QX(K,EL_EX,IRK) = RHS_QX(K,EL_EX,IRK) + W_EX*G_HAT_EX * DT_Ratio_EX
            RHS_QY(K,EL_EX,IRK) = RHS_QY(K,EL_EX,IRK) + W_EX*H_HAT_EX * DT_Ratio_Ex
        
          ENDDO

        ENDDO
        
        GOTO 1000

 100    CONTINUE

C--------------------------------------------------------------------
C Compute reflection flux
C--------------------------------------------------------------------

C.....Set the components for the tangential vector to the edge
c        write(*,*) 'computing reflection flux ',myproc,el_in,el_ex,
c     $       tsgroup(el_in),tsgroup(el_ex)
        TX = -NY
        TY =  NX

C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

        DO I = 1,NEGP

          GP_IN = I
          GP_EX = NEGP - I + 1

C ----------------------------- EL_IN ------------------------------

          ZE_IN = ZE(1,EL_IN,IRK_IN)
          QX_IN = QX(1,EL_IN,IRK_IN)
          QY_IN = QY(1,EL_IN,IRK_IN)
          HB_IN = BATHED(GP_IN,LED_IN,EL_IN)
CND
          SFAC_IN = SFACED(GP_IN,LED_IN,EL_IN)

          LZ_XX_IN = LZ(1,1,1,EL_IN)
          LZ_XY_IN = LZ(1,1,2,EL_IN)
          LZ_YX_IN = LZ(1,2,1,EL_IN)
          LZ_YY_IN = LZ(1,2,2,EL_IN)

          DO K = 2,DOF

            ZE_IN = ZE_IN + ZE(K,EL_IN,IRK_IN)*PHI_EDGE(K,GP_IN,LED_IN)
            QX_IN = QX_IN + QX(K,EL_IN,IRK_IN)*PHI_EDGE(K,GP_IN,LED_IN)
            QY_IN = QY_IN + QY(K,EL_IN,IRK_IN)*PHI_EDGE(K,GP_IN,LED_IN)

            ! LDG terms
            LZ_XX_IN = LZ_XX_IN
     &               + LZ(K,1,1,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN)
            LZ_XY_IN = LZ_XY_IN
     &               + LZ(K,1,2,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN)
            LZ_YX_IN = LZ_YX_IN
     &               + LZ(K,2,1,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN)
            LZ_YY_IN = LZ_YY_IN
     &               + LZ(K,2,2,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN)

          ENDDO

C.....Compute the velocity in the normal and tangental direction

          Q_N_INT = QX_IN*NX + QY_IN*NY
          Q_T_INT = QX_IN*TX + QY_IN*TY

C.....Reflect the velocity in the normal direction

          Q_N_EXT = -Q_N_INT
          Q_T_EXT =  Q_T_INT

C.....Compute the x and y components of the external state flow

          DEN = 1.D0/(NX*TY - NY*TX)
          QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)*DEN
          QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)*DEN

          ZE_EX = ZE_IN
          HB_EX = HB_IN
          SFAC_EX = SFAC_IN

C.....Compute the numerical flux

C          CALL NUMERICAL_FLUX()
          IF(WDFLG(EL_IN).EQ.0) THEN
            NLEQG_TMP = NLEQG
            NLEQG = 0.D0
            G_TMP = G
            G = 0.D0
            CALL NUMERICAL_FLUX()
            NLEQG = NLEQG_TMP
            G = G_TMP
          ELSE
            CALL NUMERICAL_FLUX()
          ENDIF
          F_HAT = 0.D0          ! Ensure zero mass flux

C.....Add LDG terms

          G_HAT = G_HAT + LZ_XX_IN*NX*SFAC_IN
     &                  + LZ_XY_IN*NY
          H_HAT = H_HAT + LZ_YX_IN*NX*SFAC_IN
     &                  + LZ_YY_IN*NY

C.....Compute the edge integral
          DO K = 1,DOF

            W_IN = EDFAC_IN*EDGEQ(K,GP_IN,LED_IN)

            RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT*DT_Ratio_IN
            RHS_QX(K,EL_IN,IRK) = RHS_QX(K,EL_IN,IRK) - W_IN*G_HAT*DT_Ratio_IN
            RHS_QY(K,EL_IN,IRK) = RHS_QY(K,EL_IN,IRK) - W_IN*H_HAT*DT_Ratio_IN

          ENDDO

C ----------------------------- EL_EX ------------------------------

          ZE_EX = ZE(1,EL_EX,IRK_EX)
          QX_EX = QX(1,EL_EX,IRK_EX)
          QY_EX = QY(1,EL_EX,IRK_EX)
          HB_EX = BATHED(GP_EX,LED_EX,EL_EX)
          SFAC_EX = SFACED(GP_EX,LED_EX,EL_EX)

          LZ_XX_EX = LZ(1,1,1,EL_EX)
          LZ_XY_EX = LZ(1,1,2,EL_EX)
          LZ_YX_EX = LZ(1,2,1,EL_EX)
          LZ_YY_EX = LZ(1,2,2,EL_EX)

          DO K = 2,DOF

            ZE_EX = ZE_EX + ZE(K,EL_EX,IRK_EX)*PHI_EDGE(K,GP_EX,LED_EX)
            QX_EX = QX_EX + QX(K,EL_EX,IRK_EX)*PHI_EDGE(K,GP_EX,LED_EX)
            QY_EX = QY_EX + QY(K,EL_EX,IRK_EX)*PHI_EDGE(K,GP_EX,LED_EX)

            ! LDG terms
            LZ_XX_EX = LZ_XX_EX
     &               + LZ(K,1,1,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX)
            LZ_XY_EX = LZ_XY_EX
     &               + LZ(K,1,2,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX)
            LZ_YX_EX = LZ_YX_EX
     &               + LZ(K,2,1,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX)
            LZ_YY_EX = LZ_YY_EX
     &               + LZ(K,2,2,EL_EX)*PHI_EDGE(K,GP_EX,LED_EX)

          ENDDO

C.....Compute the velocity in the normal and tangental direction

          Q_N_EXT = QX_EX*NX + QY_EX*NY
          Q_T_EXT = QX_EX*TX + QY_EX*TY

C.....Reflect the velocity in the normal direction

          Q_N_INT = -Q_N_EXT
          Q_T_INT =  Q_T_EXT

C.....Compute the x and y components of the external state flow

          DEN = 1.D0/(NX*TY - NY*TX)
          QX_IN = ( TY*Q_N_INT - NY*Q_T_INT)*DEN
          QY_IN = (-TX*Q_N_INT + NX*Q_T_INT)*DEN

          ZE_IN = ZE_EX
          HB_IN = HB_EX
          SFAC_IN = SFAC_EX

C.....Compute the numerical flux

C          CALL NUMERICAL_FLUX()
          IF(WDFLG(EL_EX).EQ.0) THEN
            NLEQG_TMP = NLEQG
            NLEQG = 0.D0
            G_TMP = G
            G = 0.D0
            CALL NUMERICAL_FLUX()
            NLEQG = NLEQG_TMP
            G = G_TMP
          ELSE
            CALL NUMERICAL_FLUX()
          ENDIF
          F_HAT = 0.D0 ! Ensure zero mass flux

C.....Add LDG terms

          G_HAT = G_HAT + LZ_XX_EX*NX*SFAC_EX
     &                  + LZ_XY_EX*NY
          H_HAT = H_HAT + LZ_YX_EX*NX*SFAC_EX
     &                  + LZ_YY_EX*NY

C.....Compute the edge integral
          DO K = 1,DOF

            W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX)

            RHS_ZE(K,EL_EX,IRK) = RHS_ZE(K,EL_EX,IRK) + W_EX*F_HAT*DT_Ratio_EX
            RHS_QX(K,EL_EX,IRK) = RHS_QX(K,EL_EX,IRK) + W_EX*G_HAT*DT_Ratio_EX
            RHS_QY(K,EL_EX,IRK) = RHS_QY(K,EL_EX,IRK) + W_EX*H_HAT*DT_Ratio_EX

          ENDDO
        ENDDO
        
1000  CONTINUE

      RETURN
      END SUBROUTINE
