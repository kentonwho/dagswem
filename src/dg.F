C***********************************************************************
C
C     MODULE DG
C
C     This module declares all the necessary variables for DG
C
C     Written by Ethan Kubatko (3-10-05)
C
C-----------------------------------------------------------------------
C
C  Modification history on hp_DG_ADCIRC after v9_sb1
C     
C  v9_sb1     - 08/05    - sb - Parallel runs
C  v9_sb2     - 08/05    - sb - Wetting and drying
C  v9_sb2.2.1 - 08/16/05 - sb - Hard bottom implementation
C  v10_sb5      - Oct 03 - sb - Consolidate Ethan's slope limiter
C
C***********************************************************************

      MODULE DG
      
C.....Use appropriate modules

      USE SIZES

C.....Declare integer variables

      INTEGER DGFLAG, DGHOT, DGHOTSPOOL, DOF
      INTEGER EL
      INTEGER IRK
      INTEGER J1, J2, J3
      INTEGER NAGP, NCHECK, NEGP, NEDGES, NRK
      INTEGER NIEDS, NLEDS, NEEDS, NFEDS, NREDS, NEBEDS, NIBEDS
      INTEGER NIBSEG, NEBSEG
      INTEGER MNED, MNLED, MNSED, MNRAED, MNRIED, MODAL_IC
      INTEGER P, P_READ
      INTEGER SLOPEFLAG
      INTEGER FLUXTYPE
      INTEGER RK_STAGE, RK_ORDER
c
      integer iwrite

C.....Declare real variables

      REAL(SZ) C13, C16
      REAL(SZ) DOT, DHB_X, DHB_Y
      REAL(SZ) EFA_GP, EMO_GP
      REAL(SZ) EL_ANG
      REAL(SZ) FG_L
      REAL(SZ) HB_IN, HB_EX, H_TRI
      REAL(SZ) MAG1, MAG2
      REAL(SZ) NX, NY
cnd
      REAL(SZ) SFAC_IN,SFAC_EX
      REAL(SZ) RAMPDG
      REAL(SZ) S1, S2, SAV, SOURCE_X, SOURCE_Y
      REAL(SZ) TIMEDG, TIMEH_DG, TK
      REAL(SZ) QX_EX, QX_IN, QY_EX, QY_IN
      REAL(SZ) QNAM_GP, QNPH_GP
      REAL(SZ) SL2_M, SL2_NYU
      REAL(SZ) SL3_MD
      REAL(SZ) UMAG
      REAL(SZ) WSX_GP, WSY_GP
      REAL(SZ) ZE_EX, ZE_IN
Csb...Wetting and drying
      INTEGER, ALLOCATABLE :: WDFLG(:)  ! = 1 if wet, =0 if dry
      INTEGER, ALLOCATABLE :: WDFLG_TMP(:)
      INTEGER, ALLOCATABLE :: DOFS(:)
      INTEGER, ALLOCATABLE :: EL_UPDATED(:)
      REAL(SZ), ALLOCATABLE :: LEDGE_NVEC(:,:,:)
                                ! LEDGE_NVEC(1,1:3,1:NELEM) = whether a node is on the land boundary
                                ! LEDGE_NVEC(2,1:3,1:NELEM) = x component of the normal vector
                                ! LEDGE_NVEC(3,1:3,1:NELEM) = y component of the normal vector

C.....Declare real variable arrays

      REAL(SZ) DRPSI(3), DSPSI(3)
      REAL(SZ) VEC1(2), VEC2(2)
cnd
      REAL(SZ) TIMINGS(20)
      
C.....Declare allocatable integer arrays

      INTEGER, ALLOCATABLE :: NCOUNT(:)
      INTEGER, ALLOCATABLE :: NEDEL(:,:), NEDNO(:,:), NEDSD(:,:)
      INTEGER, ALLOCATABLE :: NIEDN(:), NLEDN(:), NEEDN(:)
      INTEGER, ALLOCATABLE :: NFEDN(:), NREDN(:), NEBEDN(:), NIBEDN(:)
      INTEGER, ALLOCATABLE :: NIBSEGN(:,:)
      INTEGER, ALLOCATABLE :: NEBSEGN(:)
      INTEGER, ALLOCATABLE :: EL_NBORS(:,:)
c      INTEGER, ALLOCATABLE :: BACKNODES(:,:)
      INTEGER, ALLOCATABLE :: MARK(:)

C.....Declare allocatable real arrays
      
      REAL(SZ), ALLOCATABLE :: ATVD(:,:), BTVD(:,:), CTVD(:,:)
      REAL(SZ), ALLOCATABLE :: DTVD(:), MAX_BOA_DT(:)
      REAL(SZ), ALLOCATABLE :: BATH(:,:), DBATHDX(:,:), DBATHDY(:,:)
      REAL(SZ), ALLOCATABLE :: SFAC_ELEM(:,:)
      REAL(SZ), ALLOCATABLE :: BATHED(:,:,:),SFACED(:,:,:)
      REAL(SZ), ALLOCATABLE :: COSNX(:), SINNX(:)
      REAL(SZ), ALLOCATABLE :: DP_NODE(:,:)
      REAL(SZ), ALLOCATABLE :: DP_VOL(:)
      REAL(SZ), ALLOCATABLE :: DRPHI(:,:), DSPHI(:,:)
      REAL(SZ), ALLOCATABLE :: DRDX(:), DSDX(:), DRDY(:), DSDY(:)
      REAL(SZ), ALLOCATABLE :: DXPHI2(:,:,:), DYPHI2(:,:,:), PHI2(:,:,:)
      REAL(SZ), ALLOCATABLE :: EFA_DG(:,:,:), EMO_DG(:,:,:)
      REAL(SZ), ALLOCATABLE :: UFA_DG(:,:,:), UMO_DG(:,:,:)
      REAL(SZ), ALLOCATABLE :: VFA_DG(:,:,:), VMO_DG(:,:,:)
      REAL(SZ), ALLOCATABLE :: XLEN(:)
      REAL(SZ), ALLOCATABLE :: HB(:,:,:)
      REAL(SZ), ALLOCATABLE :: MANN(:,:)
      REAL(SZ), ALLOCATABLE :: IBHT(:), EBHT(:)
      REAL(SZ), ALLOCATABLE :: EBCFSP(:), IBCFSP(:), IBCFSB(:)
      REAL(SZ), ALLOCATABLE :: JACOBI(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: M_INV(:)
      REAL(SZ), ALLOCATABLE :: PHI_AREA(:,:), PHI_EDGE(:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_CENTER(:), PHI_CORNER(:,:)
      REAL(SZ), ALLOCATABLE :: PHI_CHECK(:,:)
      REAL(SZ), ALLOCATABLE :: PHI_CORNER1(:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_MID(:,:)
      REAL(SZ), ALLOCATABLE :: PHI_INTEGRATED(:)
      REAL(SZ), ALLOCATABLE :: PSI_CHECK(:,:)
      REAL(SZ), ALLOCATABLE :: PSI1(:), PSI2(:), PSI3(:)
      REAL(SZ), ALLOCATABLE :: Q_HAT(:)
      REAL(SZ), ALLOCATABLE :: QIB(:)
      REAL(SZ), ALLOCATABLE :: QX(:,:,:), QY(:,:,:), ZE(:,:,:)
      REAL(SZ), ALLOCATABLE :: CORI_EL(:), FRIC_EL(:)
c      REAL(SZ), ALLOCATABLE :: ZE_MAX(:),ZE_MIN(:),DPE_MIN(:)
      REAL(SZ), ALLOCATABLE :: DPE_MIN(:)
      REAL(SZ), ALLOCATABLE :: WATER_DEPTH_OLD(:,:),WATER_DEPTH(:,:)
      REAL(SZ), ALLOCATABLE :: ADVECTQX(:),ADVECTQY(:)
      REAL(SZ), ALLOCATABLE :: SOURCEQX(:),SOURCEQY(:)
      REAL(SZ), ALLOCATABLE :: LZ(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: QNAM_DG(:,:,:), QNPH_DG(:,:,:)
      REAL(SZ), ALLOCATABLE :: RHS_ZE(:,:,:)
      REAL(SZ), ALLOCATABLE :: RHS_QX(:,:,:),RHS_QY(:,:,:)
      REAL(SZ), ALLOCATABLE :: XAGP(:), YAGP(:), WAGP(:)
      REAL(SZ), ALLOCATABLE :: XEGP(:), YEGP(:), WEGP(:)
      REAL(SZ), ALLOCATABLE :: SL3(:,:)
      REAL(SZ), ALLOCATABLE :: XBC(:), YBC(:)
      REAL(SZ), ALLOCATABLE :: XFAC(:,:,:), YFAC(:,:,:), SRFAC(:,:,:)
      REAL(SZ), ALLOCATABLE :: EDGEQ(:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI(:), DPHIDZ1(:), DPHIDZ2(:)
      REAL(SZ), ALLOCATABLE :: PHI_STAE(:,:), PHI_STAV(:,:)

! CJT - Rainfall
      REAL(SZ) SOURCE_Z
      REAL(SZ), ALLOCATABLE :: RI1(:),RI2(:),RIN1(:),RIN2(:),Rain(:)  

! CJT - Multi-rate time-stepping
      INTEGER(4) Mint, Ndom, Lts_Flag, Mbar, MaxNdom, dtcheck,LTSmodflag,
     $     nmodtimes,TsGroupmod(20)
      REAL(SZ) tmod1(20),tmod2(20)
      REAL(SZ) dt_max, dt_min,dt_min0,cflfac,slammin,slammax,sfeamin,sfeamax
      INTEGER(4), ALLOCATABLE :: Msteps(:), Nelem(:), TsGroup(:), ElemNum(:,:)
      INTEGER(4), ALLOCATABLE :: TsGroupNode(:,:), TsGroupSave(:)
      REAL(SZ), ALLOCATABLE :: DtLoc(:), DtDom(:), tnkp1_a(:), tnkp1_h(:)
      REAL(SZ), ALLOCATABLE :: tn_a(:), tn_h(:)
      REAL(SZ), ALLOCATABLE :: CENTROID_X(:), CENTROID_Y(:)
      REAL(SZ), ALLOCATABLE :: UFLUX(:), UFLUXOLD(:)
      INTEGER(4),  ALLOCATABLE :: TsGroupNodetemp(:)
      LOGICAL, ALLOCATABLE:: TSGROUPFLAG(:)
      INTEGER  CALLRESETFLAG1,CALLRESETFLAG2,CALLRESETFLAG3 
c
C********************** END OF DATA DECLARATIONS ***********************

      CONTAINS

C...,.Rainfall (CJT)
      SUBROUTINE ALLOC_WATER_SOURCE()
      ALLOCATE ( RI1(MNP),RI2(MNP),RIN1(MNP),RIN2(MNP),Rain(MNP) ) 
      RETURN
      END SUBROUTINE

C.....Multi-rate (CJT)
      SUBROUTINE ALLOC_MULTI_RATE()
      ALLOCATE ( Msteps(MaxNdom), DtLoc(MNE), DtDom(MaxNdom) )
      ALLOCATE ( TsGroup(MNE), ElemNum(MaxNdom,MNE), Nelem(MaxNdom) )
      ALLOCATE (TsGroupNode(MNP,MaxNdom), TsGroupSave(MNE))
      ALLOCATE ( tnkp1_a(MaxNdom), tnkp1_h(MaxNdom), tn_a(MaxNdom), 
     $     tn_h(MaxNdom) )
      ALLOCATE ( CENTROID_X(MNE), CENTROID_Y(MNE) ) 
      ALLOCATE ( UFLUX(MNE), TsGroupNodetemp(MNP), UFLUXOLD(MNE) ) 
      ALLOCATE (TSGROUPFLAG(MNE))
      uflux(:)=0
      RETURN
      END SUBROUTINE
      
C.....Set edge array sizes

      SUBROUTINE ALLOC_EDGES0()
      ALLOCATE ( IBHT(3*MNE), EBHT(3*MNE) )
      ALLOCATE ( EBCFSP(3*MNE), IBCFSP(3*MNE), IBCFSB(3*MNE) )
c      ALLOCATE ( BACKNODES(2,3*MNE) )
c      ALLOCATE ( BACKNODES(mnbou,mnvel) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_EDGES1()
      ALLOCATE ( NEDNO(2,MNED), NEDEL(2,MNED), NEDSD(2,MNED) )
      ALLOCATE ( NIBSEGN(2,MNED) )
      ALLOCATE ( NEBSEGN(MNED) )
      ALLOCATE ( NIEDN(MNED), NLEDN(MNED), NEEDN(MNED) )
      ALLOCATE ( NFEDN(MNED),  NREDN(MNED), NIBEDN(MNED), NEBEDN(MNED) )
      ALLOCATE ( NCOUNT(MNED) )
      ALLOCATE ( COSNX(MNED), SINNX(MNED), XLEN(MNED) )
      ALLOCATE ( Q_HAT(MNED) )
      RETURN
      END SUBROUTINE
      
C.....Set DG SWE array sizes

      SUBROUTINE ALLOC_DG1(MNBFR)
      ALLOCATE ( EFA_DG(MNBFR,NEEDS,2), EMO_DG(MNBFR,NEEDS,2) )
      ALLOCATE ( UFA_DG(MNBFR,NEEDS,2), UMO_DG(MNBFR,NEEDS,2) )
      ALLOCATE ( VFA_DG(MNBFR,NEEDS,2), VMO_DG(MNBFR,NEEDS,2) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_DG2(MNFFR)
      ALLOCATE ( QNAM_DG(MNFFR,NFEDS,2), QNPH_DG(MNFFR,NFEDS,2) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_DG3(MNP)
      ALLOCATE ( QIB(MNP) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_DG4()
Csb-20070228 NRK+1-->NRK+2 --- XX(:,:,NRK+2) will be used by slope limiter
      ALLOCATE ( HB(DOF,MNE,NRK+2) )
      ALLOCATE ( MANN(DOF,MNE) )
      ALLOCATE ( QY(DOF,MNE,NRK+2), QX(DOF,MNE,NRK+2) )
      ALLOCATE ( ZE(DOF,MNE,NRK+2) )
Csb-20060711 For wet/dry
c      ALLOCATE ( ZE_MAX(MNE),ZE_MIN(MNE),DPE_MIN(MNE) )
      ALLOCATE (DPE_MIN(MNE) )
      ALLOCATE ( WATER_DEPTH(MNE,3), WATER_DEPTH_OLD(MNE,3))
      ALLOCATE ( ADVECTQX(MNE), ADVECTQY(MNE),
     $     SOURCEQX(MNE),SOURCEQY(MNE))
      ALLOCATE ( MARK(MNE))
C--
Csb-20070101
      ALLOCATE ( LZ(DOF,2,2,MNE) )
C--
      ALLOCATE ( RHS_QX(DOF,MNE,NRK), RHS_QY(DOF,MNE,NRK) )
      ALLOCATE ( RHS_ZE(DOF,MNE,NRK) )
      ALLOCATE ( DRDX(MNE), DSDX(MNE), DRDY(MNE), DSDY(MNE) )
      ALLOCATE ( CORI_EL(MNE), FRIC_EL(MNE) )
      ALLOCATE ( PHI(DOF), DPHIDZ1(DOF), DPHIDZ2(DOF) )
      RETURN
      END SUBROUTINE

C.....Set RK time scheme parameters array sizes

      SUBROUTINE ALLOC_RK()
      ALLOCATE( ATVD(NRK,NRK), BTVD(NRK,NRK), CTVD(NRK,NRK)  )
      ALLOCATE( DTVD(NRK), MAX_BOA_DT(NRK) )
      RETURN
      END SUBROUTINE ALLOC_RK

C.....Set sizes for arrays used in orthobasis

      SUBROUTINE ALLOC_JACOBI()
      ALLOCATE ( JACOBI(P+1,2*P+3,2,NAGP+1) )
      ALLOCATE ( DXPHI2(P+1,P+1,NAGP+1),DYPHI2(P+1,P+1,NAGP+1) )
      ALLOCATE ( PHI2(P+1,P+1,NAGP+1) )
      ALLOCATE ( PHI_CORNER1(P+1,P+1,3) )
      RETURN
      END SUBROUTINE
      
C.....Set sizes for arrays for area integrals
      
      SUBROUTINE ALLOC_AREA_GAUSS()
      ALLOCATE ( XAGP(NAGP), YAGP(NAGP), WAGP(NAGP) )
      ALLOCATE ( PHI_AREA(DOF,NAGP) )
      ALLOCATE ( DSPHI(DOF,NAGP),DRPHI(DOF,NAGP) )
      ALLOCATE ( PHI_CORNER(DOF,3),PHI_MID(DOF,3) )
      ALLOCATE ( PHI_CENTER(DOF) )
      ALLOCATE ( PSI1(NAGP),PSI2(NAGP),PSI3(NAGP) )
      ALLOCATE ( BATH(NAGP,MNE), DBATHDX(NAGP,MNE), DBATHDY(NAGP,MNE) )
      ALLOCATE ( SFAC_ELEM(NAGP,MNE) )
      ALLOCATE ( XFAC(DOF,NAGP,MNE), YFAC(DOF,NAGP,MNE) )
      ALLOCATE ( SRFAC(DOF,NAGP,MNE) )
      RETURN
      END SUBROUTINE
      
C.....Set sizes for arrays for edge integrals

      SUBROUTINE ALLOC_EDGE_GAUSS()
      ALLOCATE ( XEGP(NEGP), WEGP(NEGP) )
      ALLOCATE ( PHI_EDGE(DOF,NEGP,3) )
      ALLOCATE ( M_INV(DOF) )
      ALLOCATE ( BATHED(NEGP,3,MNE),SFACED(NEGP,3,MNE) )
      ALLOCATE ( EDGEQ(DOF,NEGP,3) )
      RETURN
      END SUBROUTINE

C.....Set sizes for the arrays for the slope limiter

      SUBROUTINE ALLOC_SLOPELIM()
      ALLOCATE ( XBC(MNE), YBC(MNE) )
      ALLOCATE ( EL_NBORS(4,MNE) )
      ALLOCATE ( SL3(3,MNE) )
      RETURN
      END SUBROUTINE

Csb...Set sizes for arrays for wetting and drying
      SUBROUTINE ALLOC_DG_WETDRY()
      ALLOCATE ( WDFLG(MNE), DOFS(MNE) )
      ALLOCATE ( EL_UPDATED(MNE) )
      ALLOCATE ( WDFLG_TMP(MNE) )
      ALLOCATE ( LEDGE_NVEC(3,3,MNE) )
      ALLOCATE ( DP_VOL(MNE) )
      ALLOCATE ( PHI_INTEGRATED(DOF) )
      ALLOCATE ( PHI_CHECK(DOF,NCHECK) )
      ALLOCATE ( DP_NODE(NCHECK,MNE) )
      ALLOCATE ( PSI_CHECK(3,NCHECK) )
      RETURN
      END SUBROUTINE
      

      SUBROUTINE ALLOC_STAE(L)
      ALLOCATE ( PHI_STAE(DOF,L) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_STAV(L)
      ALLOCATE ( PHI_STAV(DOF,L) )
      RETURN
      END SUBROUTINE
      

      END MODULE DG


