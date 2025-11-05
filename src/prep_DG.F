C***********************************************************************
C
C     SUBROUTINE DG_PREP()
C
C     This subroutine does preparatory stuff for DG
C
C     Written by Ethan Kubatko (03-07-2005)
C
C-----------------------------------------------------------------------
C
C  mod history of hp_ADCIRC since v9
C
C  v9_sb1       - Aug 05 - sb - parallelized
C  v9_sb2       - Aug 05 - sb - wetting/drying
C  v10_sb1      - Aug 05 - sb - reflect fort.12as an initial surface elevation
C  v10_sb5      - Oct 03 - sb - consolidate Ethan's slope limiter
C  v10_sb5      - Oct 24 - sb - DG.65 output is added
C               - Aug 29, 2007 - Fixed bug for startdry = 1
C
C***********************************************************************

      SUBROUTINE PREP_DG()

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE WIND
      USE NodalAttributes, ONLY : STARTDRY, FRIC, GeoidOffset,
     $     LoadGeoidOffset,LoadManningsN,ManningsN,
     &     LoadDirEffRLen,
     &     ApplyDirectionalWindReduction
      USE OWIWIND,ONLY : NWS12INIT,NWS12GET

#ifdef CMPI
      USE MESSENGER_NEW
      USE MESSENGER_node_new
#endif

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER II, L, P_0, DOF_0,lvl
      REAL(SZ) AREA, ANGLE_SUM, HBB(3), CASUM, DPHIDX, DPHIDY, DP_MIN
      REAL(SZ) XI, YI, ZE1, ZE2, ZE3
      REAL(SZ) ZP(3) 
      INTEGER IT
      REAL(SZ) WindDragLimit    !jgf46.01 Cap on max wind drag coefficient
      REAL(SZ) RampExtFlux1   ! Ramp for external flux b.c.s @ITHS-1
      REAL(SZ) RampExtFlux2   ! Ramp for external flux b.c.s @ITHS
      REAL(SZ) RampIntFlux1   ! Ramp for internal flux b.c.s @ITHS-1
      REAL(SZ) RampIntFlux2   ! Ramp for internal flux b.c.s @ITHS
      REAL(SZ) RampTip2    ! Ramp for tidal potential @ITHS
      REAL(SZ) RampMete2   ! Ramp for wind and atmospheric pressure @ITHS
      REAL(SZ) RampWRad2   ! Ramp for wave radiation stress @ITHS


cnd
      real(sz) checkarea,arint(2,2),rhsint(2),edgeint,dpsdx,psimid,
     $     determ,sfacdub2max,sfacdub3max,R,YELEM,YED
      integer i1,i2,sfac_flag,led,ged
      integer ifac2max,ifac3max
!  CND
      real(sz) func_H,func_qx,func_qy, qxedge,qyedge,qnormal,zedge,hedge
!!!!!!!!!!!
      
      C13 = 1.D0/3.D0
      C16 = 1.D0/6.D0
      R=6378206.4d0
Csb-PDG1 moved from other places
C.....Compute the degrees of freedom per element

      DOF = (P+1)*(P+2)/2
      P_0 = P
      DOF_0 = DOF

C.....Compute the number of gauss points needed for the edge integrals
      
      NEGP = P + 1
      IF (P.EQ.0) THEN
        P     = 1
        NEGP  = 2
        DOF   = 3
        P_0   = 0
        DOF_0 = 1
      ENDIF
cnd
      iwrite=0
      
C.....Initilization for parallel DG run

#ifdef CMPI
      CALL MSG_TYPES_ELEM()   ! Determine Word Sizes for Message-Passing
      CALL MSG_TABLE_ELEM()   ! Read Message-Passing Tables
      CALL MSG_TYPES()
CND
      IF (SLOPEFLAG.EQ.5) THEN
         CALL MSG_TABLE()
      ENDIF
CND
#endif

C.....Create the edge based data

      IF(MYPROC.EQ.0) THEN
         PRINT*, 'CREATING EDGE DATA...'
         PRINT*, ''
      ENDIF
      CALL CREATE_EDGE_DATA()
      IF(MYPROC.EQ.0) THEN
         print *, 'CREATING EDGE DATA DONE'
         print *, ''
      ENDIF
#ifdef CMPI
      CALL MESSAGE_START_ELEM()  ! Startup persistent message passing
      if (slopeflag.eq.5) CALL MESSAGE_START()
#endif

C.....Re-arrange elevation specified boundary segment data for DG

      IF (NEEDS.GT.0) THEN
        CALL ALLOC_DG1(MNBFR)
        II = 1
        JJ = 1
        DO I = 1,NBFR
          DO J = 1,NOPE
            DO K = 1,NVDLL(J)-1
              EMO_DG(I,II,1) = EMO(I,JJ)
              EMO_DG(I,II,2) = EMO(I,JJ+1)
              EFA_DG(I,II,1) = EFA(I,JJ)
              EFA_DG(I,II,2) = EFA(I,JJ+1)
              UMO_DG(I,II,1) = UMO(I,JJ)
              UMO_DG(I,II,2) = UMO(I,JJ+1)
              UFA_DG(I,II,1) = UFA(I,JJ)
              UFA_DG(I,II,2) = UFA(I,JJ+1)
              VMO_DG(I,II,1) = VMO(I,JJ)
              VMO_DG(I,II,2) = VMO(I,JJ+1)
              VFA_DG(I,II,1) = VFA(I,JJ)
              VFA_DG(I,II,2) = VFA(I,JJ+1)
              II = II + 1
              JJ = JJ + 1
            ENDDO
            JJ = JJ + 1
          ENDDO
          II = 1
          JJ = 1
        ENDDO
      ENDIF

C.....Re-arrange non-zero flow specified boundary segment data for DG

      IF (NFEDS.GT.0) THEN
        CALL ALLOC_DG2(MNFFR)
        II = 1
        JJ = 1
        DO I = 1,NFFR
          DO J = 1,NBOU
            IF ( (SEGTYPE(J).EQ.2 ).OR.(SEGTYPE(J).EQ.12)
     &                             .OR.(SEGTYPE(J).EQ.22) ) THEN
              DO K = 1,NVELL(J)-1
                QNAM_DG(I,II,1) = QNAM(I,JJ)
                QNAM_DG(I,II,2) = QNAM(I,JJ+1)
                QNPH_DG(I,II,1) = QNPH(I,JJ)
                QNPH_DG(I,II,2) = QNPH(I,JJ+1)
                II = II + 1
                JJ = JJ + 1
              ENDDO
              JJ = JJ + 1
            ENDIF
          ENDDO
          II = 1
          JJ = 1
        ENDDO
      ENDIF
      
C.....If there are internal barriers allocate some stuff

      IF (NIBEDS.NE.0) CALL ALLOC_DG3(MNP)

C.....Allocate the array for node to element table

      CALL ALLOC_NNOEL1(MNP)

C.....Determine the number of elements connected at each node

      EL_COUNT = 0
      MAXEL = 1
      DO K = 1,3
        DO J = 1,MNE
            N1 = NM(J,K)
            EL_COUNT(N1) = EL_COUNT(N1) + 1
        ENDDO
      ENDDO
      MAXEL = MAXVAL(EL_COUNT)

C.....Allocate the array for the node to element table

      CALL ALLOC_NNOEL2(MNP,MAXEL)
      
C.....Construct node to element table

      EL_COUNT = 0
      DO K = 1,3
        DO J = 1,MNE
          N1 = NM(J,K)
          NNOEL(N1,1+EL_COUNT(N1)) = J
          EL_COUNT(N1) = EL_COUNT(N1) + 1
        ENDDO
      ENDDO

C.....Construct node to element angle table

      DO I = 1,MNP
        ETAMAX(I) = -99999
        KK = 1
        ELETAB(I,1) = I
        S1 = SFAC(I)
        J1 = NEITAB(I,1)
        DO 111 K = 1,NNEIGH(I)-1
          J2 = NEITAB(I,1+K)
          IF (K.LT.(NNEIGH(I)-1)) THEN
            J3 = NEITAB(I,2+K)
          ELSE
            J3 = NEITAB(I,2)
          ENDIF
          DO J = 1,EL_COUNT(I)
            EL = NNOEL(I,J)
            N1 = NM(EL,1)
            N2 = NM(EL,2)
            N3 = NM(EL,3)
            IF ((J1.EQ.N1).OR.(J1.EQ.N2).OR.(J1.EQ.N3)) THEN
              IF ((J2.EQ.N1).OR.(J2.EQ.N2).OR.(J2.EQ.N3)) THEN
                IF ((J3.EQ.N1).OR.(J3.EQ.N2).OR.(J3.EQ.N3)) THEN
                  ELETAB(I,1+KK) = EL
                  S2  = SFAC(J2)
                  SAV = (S1 + S2)/2.D0
                  VEC1(1) =      X(J1) - X(J2)
                  VEC1(2) = SAV*(Y(J1) - Y(J2))
                  VEC2(1) =      X(J1) - X(J3)
                  VEC2(2) = SAV*(Y(J1) - Y(J3))
                  MAG1 = SQRT(VEC1(1)**2 + VEC1(2)**2)
                  MAG2 = SQRT(VEC2(1)**2 + VEC2(2)**2)
                  DOT = DOT_PRODUCT(VEC1,VEC2)
                  EL_ANG  = ACOS(DOT/(MAG1*MAG2))
                  ANGTAB(I,KK+1) = RAD2DEG*EL_ANG
                  KK = KK + 1
                  GOTO 111
                ENDIF
              ENDIF
            ENDIF
          ENDDO
111     CONTINUE
      ENDDO

C.....Obtain RK time scheme parameters

      CALL RK_TIME()
      
C.....Allocate some DG stuff

      CALL ALLOC_DG4()
      
C.....Initialize the DG arrays

      ZE = 0.D0
      QX = 0.D0
      QY = 0.D0
      HB = 0.D0
      LZ = 0.d0
      MARK = 0
      
      RHS_ZE(:,:,:) = 0.D0
      RHS_QX(:,:,:) = 0.D0
      RHS_QY(:,:,:) = 0.D0
      
      WSX2(:) = 0
      WSY2(:) = 0

C.....If using modal initial conditions transform the bathymetry from
C.....nodal coordinates to modal dof
      timedg=0.0

       DO J = 1,NE
         N1 = NM(J,1)
         N2 = NM(J,2)
         N3 = NM(J,3)
         x1=x(n1)
         y1=y(n1)
         x2=x(n2)
         y2=y(n2)
         x3=x(n3)
         y3=y(n3)
         HB(1,J,1) =  1.D0/3.D0*(DP(N1) + DP(N2) + DP(N3))
         HB(2,J,1) = -1.D0/6.D0*(DP(N1) + DP(N2)) + 1.D0/3.D0*DP(N3)
         HB(3,J,1) = -0.5D0*DP(N1) + 0.5D0*DP(N2)
!  CND  for analytic test case only
!         ze(1,J,1)=1.D0/3.D0*(func_H(x1,y1,timedg) 
!     $        + func_H(x2,y2,timedg) + func_H(x3,y3,timedg))-2.0
!         ze(1,J,1)=func_H(centroid_x(j),centroid_y(j),timedg)-2.0
!         ze(2,J,1) = -1.D0/6.D0*(func_H(x1,y1,timedg) 
!     $        + func_H(x2,y2,timedg)) + 1.D0/3.D0*func_H(x3,y3,timedg)
!         ze(3,J,1) = -0.5D0*func_H(x1,y1,timedg) + 0.5D0*func_H(x2,y2,timedg)
!
!         qx(1,J,1)=1.D0/3.D0*(func_qx(x1,y1,timedg) 
!     $        + func_qx(x2,y2,timedg) + func_qx(x3,y3,timedg))
!         qx(1,J,1)=func_qx(centroid_x(j),centroid_y(j),timedg) 
!         qx(2,J,1) = -1.D0/6.D0*(func_qx(x1,y1,timedg) 
!     $        + func_qx(x2,y2,timedg)) + 1.D0/3.D0*func_qx(x3,y3,timedg)
!         qx(3,J,1) = -0.5D0*func_qx(x1,y1,timedg) + 0.5D0*func_qx(x2,y2,timedg)
! 
!         qy(1,j,1)=1.D0/3.D0*(func_qy(x1,y1,timedg) 
!     $        + func_qy(x2,y2,timedg) + func_qy(x3,y3,timedg))
!         qy(1,j,1)=func_qy(centroid_x(j),centroid_y(j),timedg) 
!         qy(2,J,1) = -1.D0/6.D0*(func_qy(x1,y1,timedg) 
!     $        + func_qy(x2,y2,timedg)) + 1.D0/3.D0*func_qy(x3,y3,timedg)
!         qy(3,J,1) = -0.5D0*func_qy(x1,y1,timedg) + 0.5D0*func_qy(x2,y2,timedg)
cnd
c         SFACDUB(1,J)= 1.D0/3.D0*(SFAC(N1) + SFAC(N2) + SFAC(N3))
c         SFACDUB(2,J) = -1.D0/6.D0*(SFAC(N1) + SFAC(N2))
c     $        + 1.D0/3.D0*SFAC(N3)
c         SFACDUB(3,J) = -0.5D0*SFAC(N1) + 0.5D0*SFAC(N2)
cnd compute the function Y in the Dubiner basis so that we can
c   easily compute the edge and element gauss points
c
         YDUB(1,J)= 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))
         YDUB(2,J) = -1.D0/6.D0*(Y(N1) + Y(N2))
     $        + 1.D0/3.D0*Y(N3)
         YDUB(3,J) = -0.5D0*Y(N1) + 0.5D0*Y(N2)
       ENDDO
       if (LoadManningsN) then
          DO J = 1,NE
             N1 = NM(J,1)
             N2 = NM(J,2)
             N3 = NM(J,3)
             MANN(1,J) =  1.D0/3.D0*(ManningsN(N1)
     $            + ManningsN(N2) + ManningsN(N3))
             MANN(2,J) = -1.D0/6.D0*(ManningsN(N1) 
     $            + ManningsN(N2)) + 1.D0/3.D0*ManningsN(N3)
             MANN(3,J) = -0.5D0*ManningsN(N1) + 0.5D0*ManningsN(N2)
          ENDDO
       endif

      IF (MODAL_IC.EQ.0) THEN
         time_a=0.0
C  this assumes a cold start
        if (LoadGeoidOffset) then
           DO J = 1,NE
              N1 = NM(J,1)
              N2 = NM(J,2)
              N3 = NM(J,3)
              ZE(1,J,1)=1.d0/3.d0*(GeoidOffset(N1)+GeoidOffset(N2)+
     $             GeoidOffset(N3))
              IF (DOF_0.NE.1) THEN
                 ZE(2,J,1)=-1.d0/6.d0*(GeoidOffset(N1)+GeoidOffset(N2))
     $                +1.d0/3.d0*GeoidOffset(N3)
                 ZE(3,J,1)=-.5d0*GeoidOffset(N1)+.5d0*GeoidOffset(N2)
              ENDIF
           ENDDO
c  need to call wetdry in case geoidoffset gives negative water depth
        endif
      ENDIF

C.....Prep the slopelmiter

      IF (SLOPEFLAG.NE.0) THEN
        CALL ALLOC_SLOPELIM()
        CALL PREP_SLOPELIM()
      ENDIF
C--

      IF(MYPROC.EQ.0) THEN
         print *, 'PREP FOR WET/DRY BEGINS...'
      ENDIF

C.....1. Set initial surface elevation above the bed elevation
C.....   if wetting-and-drying is enabled and the initial water depth is
C.....   not specified by fort.12
C.....2. Set wet-and-dry elemental flags
C.....3. Set the DOF at dry elements = 1

      NCHECK = 3
      IF (P.GT.1) NCHECK = NCHECK + 3*NEGP
      CALL ALLOC_DG_WETDRY()
      PHI_CHECK = 0.D0
      PSI_CHECK = 0.D0
      H0L = H0
      H0H = H0 * 1.0
      HABSMIN=H0*1.0
      
C.....Wetting and drying is not turned on

      IF(NOLIFA.EQ.0.OR.NOLIFA.EQ.1) THEN
        DO J = 1,NE
          WDFLG(J) = 1
          DOFS(J) = DOF
        ENDDO
        
C.....Wetting and drying is turned on but there are no dry nodes below
C.....geoid

      ELSEIF (NOLIFA.EQ.2.AND.NSTARTDRY.EQ.0) THEN
        DO J = 1,NE
        
C.........Check to see if there are initially any dry nodes

c          ZE1 = 0.D0
c          ZE2 = 0.D0
c          ZE3 = 0.D0
          ZE1 = ZE(1,J,1)
          ZE2 = ZE(1,J,1)
          ZE3 = ZE(1,J,1)
cnd
          IF (DP(NM(J,1)).LT.H0.or.ze1<0) ZE1 = max(ze1,H0 - DP(NM(J,1)))
          IF (DP(NM(J,2)).LT.H0.or.ze2<0) ZE2 = max(ze2,H0 - DP(NM(J,2)))
          IF (DP(NM(J,3)).LT.H0.or.ze3<0) ZE3 = max(ze3,H0 - DP(NM(J,3)))
          
C.........If so set initial surface elevation values
          
          IF ((ZE1 + ZE2 + ZE3)/3.D0.NE.ZE(1,J,1)) THEN
            IF (P.EQ.0) THEN
              DP_MIN = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)))
              ZE(1,J,1) = max(ze(1,j,1),H0 - DP_MIN)
            ELSE
              IF (ZE(1,J,1).GT.(ZE1+ZE2+ZE3)/3.d0) THEN
                 IF (DOF_0.NE.1) THEN
                    ZE(2,J,1)=0.d0
                    ZE(3,J,1)=0.d0
                 ENDIF
              ELSE
                 ZE(1,J,1)=(ZE1+ZE2+ZE3)/3.D0
                 IF (DOF_0.NE.1) THEN
                    ZE(2,J,1) = -1.D0/6.D0*(ZE1 + ZE2) + 1.D0/3.D0*ZE3
                    ZE(3,J,1) = -0.5D0*ZE1 + 0.5D0*ZE2
                 ENDIF
              ENDIF
            ENDIF

C............and set wet/dry flag (0 = dry, 1 = wet)

             WDFLG(J) = 0
             DOFS(J) = DOF
          ELSE
            WDFLG(J) = 1
            DOFS(J) = DOF
          ENDIF

        ENDDO

C.....If there are dry nodes below geoid

      ELSEIF (NOLIFA.EQ.2.AND.NSTARTDRY.EQ.1) THEN
      
C.......Loop over elements
      
        DO J = 1,NE

C........Retrieve global node numbers for element

          N1 = NM(J,1)
          N2 = NM(J,2)
          N3 = NM(J,3)
          
C........Check to see if nodes are initially dry

          ZE1 = 0
          ZE2 = 0
          ZE3 = 0

C........Check to see if nodes are initially dry

          ZE1 = 0
          ZE2 = 0
          ZE3 = 0
          IF (STARTDRY(N1).EQ.1) ZE1 = H0 - DP(N1)
          IF (DP(N1).LT.H0) ZE1 = H0 - DP(N1)
          IF (STARTDRY(N2).EQ.1) ZE2 = H0 - DP(N2)
          IF (DP(N2).LT.H0) ZE2 = H0 - DP(N2)
          IF (STARTDRY(N3).EQ.1) ZE3 = H0 - DP(N3)
          IF (DP(N3).LT.H0) ZE3 = H0 - DP(N3)

          IF (MODAL_IC.EQ.3) THEN
            IF (STARTDRY(N1).EQ.-88888) then
               ZE1 = H0 - DP(N1)
            else
               ZE1 = STARTDRY(N1)
            endif
            IF (STARTDRY(N2).EQ.-88888) then
              ZE2 = H0 - DP(N2)
            else
              ZE2 = STARTDRY(N2)
            endif
            IF (STARTDRY(N3).EQ.-88888) then
              ZE3 = H0 - DP(N3)
            else
              ZE3 = STARTDRY(N3)
            endif
          ENDIF

C.........If so set initial surface elevation values

          IF ((ZE1 + ZE2 + ZE3).NE.0) THEN
            IF (P.EQ.0) THEN
              DP_MIN = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)))
              ZE(1,J,1) = max(ze(1,j,1),H0 - DP_MIN)
            ELSE
              IF (ZE(1,J,1).GT.(ZE1+ZE2+ZE3)/3.d0) THEN
                 IF (DOF_0.NE.1) THEN
                    ZE(2,J,1)=0.d0
                    ZE(3,J,1)=0.d0
                 ENDIF
              ELSE
                 ZE(1,J,1)=(ZE1+ZE2+ZE3)/3.D0
                 IF (DOF_0.NE.1) THEN
                    ZE(2,J,1) = -1.D0/6.D0*(ZE1 + ZE2) + 1.D0/3.D0*ZE3
                    ZE(3,J,1) = -0.5D0*ZE1 + 0.5D0*ZE2
                 ENDIF
              ENDIF
            ENDIF
            
C...........Set wet/dry flag (0 = dry, 1 = wet)

             WDFLG(J) = 0
             DOFS(J) = DOF
          ELSE
            WDFLG(J) = 1
            DOFS(J) = DOF
          ENDIF
        ENDDO
      ENDIF

      IF(MYPROC.EQ.0) THEN
         print *, 'DONE'
         print *, ''
      ENDIF
      
C.....Read in modal dof for initial conditions

      IF (MODAL_IC.EQ.1) THEN
        OPEN(163,FILE=DIRNAME//'/'//'Initial_Conditions.163')
        OPEN(164,FILE=DIRNAME//'/'//'Initial_Conditions.164')
        OPEN(114,FILE=DIRNAME//'/'//'Initial_Bathymetry.114')
        READ(163,*) P_READ
        READ(164,*) P_READ,P_READ
        READ(114,*) P_READ
        IF (P_READ.NE.P) THEN
          PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
          STOP
        ENDIF
        DO J = 1,MNE
          DO K = 1,DOF_0
            READ(163,*) ZE(K,J,1)
            READ(164,*) QX(K,J,1), QY(K,J,1)
            READ(114,*) HB(K,J,1)
          ENDDO
        ENDDO
        H_TRI = SQRT((X(1)-X(2))**2 + (Y(1)-Y(2))**2)
        CLOSE(163)
        CLOSE(164)
        CLOSE(114)
      ENDIF

C.....Read in modal dof for hot start conditions

      IF (MODAL_IC.EQ.2) THEN
        OPEN(263,FILE=DIRNAME//'/'//'Hot_start.263')
        OPEN(264,FILE=DIRNAME//'/'//'Hot_start.264')
        OPEN(214,FILE=DIRNAME//'/'//'Hot_start.214')
        READ(263,*) P_READ
        READ(264,*) P_READ,P_READ
        READ(214,*) ITHS, DTOLD
        time_a=iths*dtold
C  ADJUST ITHS TO ACCOUNT FOR POSSIBLE CHANGE IN DT
        ITHS=ITHS*(DTOLD/DTDP)
        if (myproc.eq.0) then
           write(*,*) 'iths ',iths
        endif
        IF (P_READ.NE.P_0) THEN
          PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
          STOP
        ENDIF
        DO J = 1,MNE
          DO K = 1,DOF_0
            READ(263,*) ZE(K,J,1)
            READ(264,*) QX(K,J,1), QY(K,J,1)
            READ(214,*) HB(K,J,1), WDFLG(J)
          ENDDO
        ENDDO
        CLOSE(263)
        CLOSE(264)
        CLOSE(214)
cnd
c$$$        IF (lts_flag.ne.0) then
c$$$           do j=1,ne
c$$$              uflux(j)=0
c$$$              do kk=1,3
c$$$                 ged=neled(kk,j)
c$$$                 led=nedsd(1,ged)
c$$$                 qnormal=0
c$$$                 zedge=0
c$$$                 do i=1,negp
c$$$                    qxedge=0
c$$$                    qyedge=0
c$$$                    do k=1,dof
c$$$                       qxedge=qxedge+qx(k,j,1)*phi_edge(k,i,led)
c$$$                       qyedge=qyedge+qy(k,j,1)*phi_edge(k,i,led)
c$$$                       zedge=zedge+ze(k,j,1)*phi_edge(k,i,led)
c$$$                    enddo
c$$$                    qnormal=qnormal+qxedge*cosnx(ged)+qyedge*sinnx(ged)
c$$$                 enddo
c$$$                 n1=nedno(1,ged)
c$$$                 n2=nedno(2,ged)
c$$$                 hedge=.5*(dp(n1)+dp(n2))+zedge/negp
c$$$                 hedge=max(hedge,h0)
c$$$                 uflux(j)=uflux(j)+(qnormal/negp)/hedge
c$$$              enddo
c$$$              uflux(j)=uflux(j)/3
c$$$           enddo
c$$$        endif
cnd
        IF (NWS.EQ.12) THEN
           WTIME1 = STATIM*86400.D0
           WTIME2 = WTIME1 + WTIMINC
           CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
           CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
           DO IT = 1,ITHS
              TIMEIT = it*dtold + STATIM*86400.D0
              IF (TIMEIT.GT.WTIME2) THEN
                 WTIME1 = WTIME2
                 WTIME2 = WTIME2+WTIMINC
                 DO I = 1,NP
                    WVNX1(I) = WVNX2(I)
                    WVNY1(I) = WVNY2(I)
                    PRN1(I)  = PRN2(I)
                 ENDDO
                 CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
              ENDIF
           ENDDO
           WTRATIO = (TIME_A - WTIME1)/WTIMINC
           DO I = 1,NP
              WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
              WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
              WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
              WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
              IF (LoadDirEffRLen) THEN
                 CALL ApplyDirectionalWindReduction(I, WDragCo,
     &                WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                 WindMag = SQRT(WindX*WindX+WindY*WindY)
                 WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
              ENDIF
              WSX2(I)    = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
              WSY2(I)    = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
              PR2(I)     = RampMete2*PRN1(I)
              WVNXOUT(I) = RampMete2*WINDX
              WVNYOUT(I) = RampMete2*WINDY
           ENDDO
        ENDIF
      ENDIF
      
C.....Retrieve the normals to the edges

      CALL CALC_NORMAL()
         

C.....Retrieve the area integral gauss quadrature points

      CALL QUAD_PTS_AREA()

C.....Retrieve the edge integral gauss quadrature points

      CALL QUAD_PTS_EDGE()

C.....Evaluate the orthogonal basis and its derivatives at the area
C.....gauss quadrature points

      CALL ORTHOBASIS_AREA()

C.....Evaluate the orthogonal basis at the edge gauss quadrature points

      CALL ORTHOBASIS_EDGE()

Csb-20060711 Added to update ZE_MAX,ZE_MIN,DPE_MIN (do we need this?)
      DO J = 1,NE
c        ZP(:) = 0.D0
c        DO K = 1,3
c          DO KK = 1,DOF
c            ZP(K)  = ZP(K)  + PHI_CORNER(KK,K)*ZE(KK,J,1)
c          ENDDO
c        ENDDO
c
c        ZE_MAX(J) = MAX(ZP(1),ZP(2),ZP(3))
c        ZE_MIN(J) = MIN(ZP(1),ZP(2),ZP(3))
        DPE_MIN(J) = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)))
      ENDDO

C.....Compute the values of the (linear) nodal basis functions at the
C.....area gauss quadrature points

      DO I = 1,NAGP
        PSI1(I) = -1.D0/2.D0*(XAGP(I) + YAGP(I))
        PSI2(I) =  1.D0/2.D0*(XAGP(I) + 1.D0)
        PSI3(I) =  1.D0/2.D0*(YAGP(I) + 1.D0)
      ENDDO
      
C.....Store the derivatives of the (linear) nodal basis functions

      DRPSI(1) = -1.D0/2.D0
      DRPSI(2) =  1.D0/2.D0
      DRPSI(3) =  0.D0
      DSPSI(1) = -1.D0/2.D0
      DSPSI(2) =  0.D0
      DSPSI(3) =  1.D0/2.D0
      
C.....Pre-compute the derivatives of the coordinate transformation for
C.....each element

      DO J = 1,MNE
      
C.....Compute the area of the element

c        AREA = 0.5D0*AREAS(J)

C.....Retrieve the global node numbers for the element

         N1 = NM(J,1)
         N2 = NM(J,2)
         N3 = NM(J,3)
         x1=x(n1)
         y1=y(n1)
         x2=x(n2)
         y2=y(n2)
         x3=x(n3)
         y3=y(n3)
         AREA = (X1 - X3)*(Y2 - Y3) + (X3 - X2)*(Y1 - Y3)
         area=area*.5d0

C.....Compute the derivatives of the coordinate transformation

         DRDX(J) = 1.D0/AREA*(Y(N3) - Y(N1))
         DSDX(J) = 1.D0/AREA*(Y(N1) - Y(N2))

         DRDY(J) = 1.D0/AREA*(X(N1) - X(N3))
         DSDY(J) = 1.D0/AREA*(X(N2) - X(N1))
        
C.......Compute elemental Coriolis and friction terms

         CORI_EL(J) = (CORIF(N1) + CORIF(N2) + CORIF(N3))/3.D0
         FRIC_EL(J) = (FRIC(N1) + FRIC(N2) + FRIC(N3))/3.D0

C.......Pre-compute the bathymetry and the gradient of the bathymetry at
C.......the quadrature points and compute volume of water

         DP_VOL(J) = 0.D0
         DO I = 1,NAGP          ! Area quadrature points
            BATH(I,J) = 0.D0
            DBATHDX(I,J) = 0.D0
            DBATHDY(I,J) = 0.D0
            YELEM = 0.d0
            DO K = 1,DOF
               DPHIDX = DRPHI(K,I)*DRDX(J) + DSPHI(K,I)*DSDX(J)
               DPHIDY = DRPHI(K,I)*DRDY(J) + DSPHI(K,I)*DSDY(J)
               XFAC(K,I,J)  = M_INV(K)*WAGP(I)*DPHIDX
               YFAC(K,I,J)  = M_INV(K)*WAGP(I)*DPHIDY
               SRFAC(K,I,J) = M_INV(K)*WAGP(I)*PHI_AREA(K,I)
               BATH(I,J) = BATH(I,J) + HB(K,J,1)*PHI_AREA(K,I)
               YELEM = YELEM + YDUB(K,J)*PHI_AREA(K,I)
cnd compute sfac at the element gauss point and store it
               IF (ICS.EQ.1) THEN
                  SFAC_ELEM(I,J)=1.0D0
               ELSE
                  SFAC_ELEM(I,J)=COS(SFEA0)/COS(YELEM/R)
               ENDIF
               DBATHDX(I,J) = DBATHDX(I,J) + HB(K,J,1)*DPHIDX
               DBATHDY(I,J) = DBATHDY(I,J) + HB(K,J,1)*DPHIDY
               DP_VOL(J) = DP_VOL(J) + WAGP(I)*HB(K,J,1)*PHI_AREA(K,I)
            ENDDO
         ENDDO
         DP_VOL(J) = 0.25D0*AREAS(J)*DP_VOL(J)
         DO L = 1,3
            DO I = 1,NEGP       ! Edge quadrature points
               BATHED(I,L,J) = 0.D0
               YED = 0.D0
               DO K = 1,DOF
                  BATHED(I,L,J)=BATHED(I,L,J)
     $                 + HB(K,J,1)*PHI_EDGE(K,I,L)
                  YED=YED + YDUB(K,J)*PHI_EDGE(K,I,L)
c                  SFACED(I,L,J)=SFACED(I,L,J)
c     $                 + SFACDUB(K,J)*PHI_EDGE(K,I,L)
cnd compute sfac at the edge gauss point and store it
                  IF (ICS.EQ.1) THEN
                     SFACED(I,L,J)=1.0d0
                  ELSE
                     SFACED(I,L,J)=COS(SFEA0)/COS(YED/R)
                  ENDIF
                  EDGEQ(K,I,L) = 2.0*M_INV(K)*PHI_EDGE(K,I,L)*WEGP(I)
               ENDDO
            ENDDO
         ENDDO
        
C........Store bathymetry at vertices and edge gauss points for wet-dry

         DO I = 1,3
            DP_NODE(I,J) = DP(NM(J,I))
         ENDDO
         IF (NCHECK.GT.3) THEN
            II = 4
            DO L = 1,3
               DO I = 1,NEGP
                  DP_NODE(II,J) = BATHED(I,L,J)
                  II = II + 1
               ENDDO
            ENDDO
         ENDIF
      ENDDO
c$$$cnd compute sfacdub
c$$$c
c$$$      if (ics.eq.1) then
c$$$         do i=1,ne
c$$$            sfacdub(1,i)=1.0d0
c$$$            if (dof.eq.3) then
c$$$               sfacdub(2,i)=0
c$$$               sfacdub(3,i)=0
c$$$            endif
c$$$         enddo
c$$$      endif
c$$$      if (ics.eq.2) then
c$$$         do i=1,ne
c$$$            checkarea=0.0d0
c$$$            edgeint=0.0
c$$$            if (dabs(dbathdx(1,i)).gt.1.d-10) then
c$$$               sfac_flag=1
c$$$               dpsdx=dbathdx(1,i)
c$$$            else
c$$$               sfac_flag=0
c$$$               dpsdx=1.0d0
c$$$            endif
c$$$c$$$            do j=1,3
c$$$c$$$               i1=nedno(1,neled(j,i))
c$$$c$$$               i2=nedno(2,neled(j,i))
c$$$c$$$               if (sfac_flag.eq.1) then
c$$$c$$$                  psimid=.5*(dp(i1)+dp(i2))
c$$$c$$$               else
c$$$c$$$                  psimid=.5*(x(i1)+x(i2))
c$$$c$$$               endif
c$$$c$$$               if (i.eq.nedel(1,neled(j,i))) then
c$$$c$$$                 edgeint=edgeint+xlen(neled(j,i))
c$$$c$$$     $                 *psimid*cosnx(neled(j,i))
c$$$c$$$               else
c$$$c$$$                 edgeint=edgeint-xlen(neled(j,i))
c$$$c$$$     $                 *psimid*cosnx(neled(j,i))
c$$$c$$$               endif
c$$$c$$$            enddo
c$$$            arint(1,1)=0.0d0
c$$$            rhsint(1)=0.0d0
c$$$            do j=1,3
c$$$               led=neled(j,i)
c$$$               do l=1,negp
c$$$                  hb_in=bathed(l,j,i)
c$$$                  if (i.eq.nedel(1,neled(j,i))) then
c$$$                     rhsint(1)=rhsint(1)+xlen(led)
c$$$     $                    *phi_edge(1,l,j)*wegp(l)*(.5+hb_in)*cosnx(led)
c$$$c                     if (i.eq.3232) then
c$$$c                        checkarea=checkarea+xlen(led)*phi_edge(1,l,j)
c$$$c     $                       *wegp(l)*(.5)*cosnx(led)
c$$$c                     endif
c$$$                  else
c$$$                     rhsint(1)=rhsint(1)-xlen(led)
c$$$     $                    *phi_edge(1,l,j)*wegp(l)*(.5+hb_in)*cosnx(led)
c$$$c                     if (i.eq.3232) then
c$$$c                        checkarea=checkarea-xlen(led)*phi_edge(1,l,j)
c$$$c     $                       *wegp(l)*(.5)*cosnx(led)
c$$$c                     endif
c$$$                  endif
c$$$               enddo
c$$$            enddo
c$$$
c$$$            do l=1,nagp
c$$$               arint(1,1)=arint(1,1)+wagp(l)*dbathdx(l,i)
c$$$            enddo
c$$$c            sfacdub(1,i)=(2.d0*edgeint)/(areas(i)*dpsdx)
c$$$            rhsint(1)=rhsint(1)*2*m_inv(1)/areas(i)
c$$$            arint(1,1)=arint(1,1)*m_inv(1)
c$$$            sfacdub(1,i)=rhsint(1)/arint(1,1)
c$$$c            if (i.eq.3232) then
c$$$c               write(1907,*) 2.*m_inv(1)*rhsint(1)/areas(i)
c$$$c               write(3232,*) rhsint(1), arint(1,1)
c$$$c            endif
c$$$c
c$$$c            sfac_flag=0
c$$$            if (dof.eq.3.and.sfac_flag.eq.0) then
c$$$               sfacdub(2,i)=0
c$$$               sfacdub(3,i)=0
c$$$            endif
c$$$c
c$$$c            if (i.eq.3232) then 
c$$$c               write(3232,*) 'checkarea ',checkarea,sfacdub(1,i)
c$$$c            endif
c$$$c
c$$$            if (dof.eq.3.and.sfac_flag.eq.1) then
c$$$               rhsint(1)=0.0
c$$$               rhsint(2)=0.0
c$$$               rtemp2(i)=0.0
c$$$               do j=1,3
c$$$                  led=neled(j,i)
c$$$                  do l=1,negp
c$$$                     hb_in=bathed(l,j,i)
c$$$                     do k=2,dof
c$$$                        if (i.eq.nedel(1,neled(j,i))) then
c$$$                           rhsint(k-1)=rhsint(k-1)+xlen(led)*G
c$$$     $                   *phi_edge(k,l,j)*wegp(l)*(.5+hb_in)*cosnx(led)
c$$$c                           rhsint(k-1)=rhsint(k-1)+xlen(led)/areas(i)
c$$$c     $                       *edgeq(k,l,j)*hb_in*cosnx(led)
c$$$c                           if (k.eq.2) then
c$$$c                              rtemp2(i)=rtemp2(i)+xlen(led)/areas(i)
c$$$c     $                             *edgeq(k,l,j)*G*(.5+hb_in)
c$$$c     $                             *cosnx(led)
c$$$c                           endif
c$$$                        else
c$$$                           rhsint(k-1)=rhsint(k-1)-xlen(led)*G
c$$$     $                   *phi_edge(k,l,j)*wegp(l)*(.5+hb_in)*cosnx(led)
c$$$c                           rhsint(k-1)=rhsint(k-1)-xlen(led)/areas(i)
c$$$c     $                       *edgeq(k,l,j)*hb_in*cosnx(led)
c$$$c                           if (k.eq.2) then
c$$$c                              rtemp2(i)=rtemp2(i)-xlen(led)/areas(i)
c$$$c     $                             *edgeq(k,l,j)*G*(.5+hb_in)
c$$$c     $                             *cosnx(led)
c$$$c                           endif
c$$$                        endif
c$$$c                        if (i.eq.173.and.k.eq.2) then
c$$$c                           write(783,*) rtemp2(i),
c$$$c     $                          j,xlen(led)/areas(i),cosnx(led),
c$$$c     $                          edgeq(2,l,j),G*(.5+hb_in)*cosnx(led)
c$$$c                        endif
c$$$                     enddo
c$$$c                     if (i.eq.1) then
c$$$c                        write(780,*) nedel(1,neled(j,i)),
c$$$c     $                   rhsint(1),xlen(led),phi_edge(k,l,j),wegp(l),
c$$$c     $                   hb_in,cosnx(led)
c$$$c                     endif
c$$$                  enddo
c$$$               enddo
c$$$               if (i.eq.3232) then
c$$$                  write(3232,*) 'in prep ',
c$$$     $                 rhsint(1)*2*m_inv(2)/areas(3232),
c$$$     $                 rhsint(2)*2*m_inv(3)/areas(3232)
c$$$               endif
c$$$               arint(1,1)=0.d0
c$$$               arint(1,2)=0.d0
c$$$               arint(2,1)=0.d0
c$$$               arint(2,2)=0.d0
c$$$               do l=1,nagp
c$$$                  do k=2,dof
c$$$                     DPHIDX = DRPHI(K,l)*DRDX(i) + DSPHI(K,l)*DSDX(i)
c$$$                     rhsint(k-1)=rhsint(k-1)
c$$$     $                    -sfacdub(1,i)*areas(i)/2.*wagp(l)*G*
c$$$     $                ((.5+bath(l,i))*dphidx+dbathdx(l,i)*phi_area(k,l))
c$$$                     do kk=2,dof
c$$$                        arint(k-1,kk-1)=arint(k-1,kk-1)+
c$$$     $                       phi_area(kk,l)*areas(i)/2.
c$$$     $                       *wagp(l)*G*
c$$$     $                ((.5+bath(l,i))*dphidx+dbathdx(l,i)*phi_area(k,l))
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$	       do k=1,2
c$$$	          rhsint(k)=2.*m_inv(k+1)*rhsint(k)/areas(i)
c$$$		  arint(k,1)=2.*m_inv(k+1)*arint(k,1)/areas(i)
c$$$		  arint(k,2)=2.*m_inv(k+1)*arint(k,2)/areas(i)
c$$$               enddo	
c$$$               determ=arint(1,1)*arint(2,2)-arint(1,2)*arint(2,1)
c$$$               if (i.eq.3232) then
c$$$                  write(3232,*) rhsint(1)
c$$$                  write(3232,*) rhsint(2)
c$$$                  write(3232,*) areas(3232),dbathdx(1,i),bath(1,i)
c$$$                  write(3232,*) arint(1,1),arint(1,2)
c$$$                  write(3232,*) arint(2,1),arint(2,2)
c$$$                  write(3232,*) determ,sfacdub(1,i)
c$$$               endif
c$$$
c$$$               if (determ.eq.0) then
c$$$                  write(*,*) 'problem in sfac ',determ,i,
c$$$     $                 arint(1,1),arint(1,2),arint(2,1),arint(2,2)
c$$$                  stop
c$$$               endif
c$$$
c$$$c
c$$$c project into Dubiner basis
c$$$               sfacdub(2,i)=(arint(2,2)*rhsint(1)-arint(1,2)*rhsint(2))
c$$$     $              /determ
c$$$               sfacdub(3,i)=(-arint(2,1)*rhsint(1)+arint(1,1)*rhsint(2))
c$$$     $              /determ
c$$$
c$$$            endif
c$$$         enddo
c$$$      endif
c$$$        write(781,*)
c$$$        sfacdub2max=0.0
c$$$        sfacdub3max=0.0
c$$$        do j=1,ne
c$$$           write(781,*) j,sfacdub(1,j),sfacdub(2,j),sfacdub(3,j)
c$$$           if (dabs(sfacdub(2,j)).gt.sfacdub2max) then
c$$$              sfacdub2max=dabs(sfacdub(2,j))
c$$$              ifac2max=j
c$$$           endif
c$$$           if (dabs(sfacdub(3,j)).gt.sfacdub3max) then
c$$$              sfacdub3max=dabs(sfacdub(3,j))
c$$$              ifac3max=j
c$$$           endif
c$$$        enddo
c$$$        write(*,*) 'sfacdub ',sfacdub(2,ifac2max),sfacdub(3,ifac3max),
c$$$     $       ifac2max,ifac3max,dbathdx(1,ifac2max)
c$$$c
      
      DO I = 1,3
        IF (I.EQ.1) THEN
          XI = -1.D0
          YI = -1.D0
        ELSEIF (I.EQ.2) THEN
          XI =  1.D0
          YI = -1.D0
        ELSE
          XI = -1.D0
          YI =  1.D0
        ENDIF
        PSI_CHECK(1,I) = -1.D0/2.D0*(XI + YI)
        PSI_CHECK(2,I) =  1.D0/2.D0*(XI + 1.D0)
        PSI_CHECK(3,I) =  1.D0/2.D0*(YI + 1.D0)
        DO K = 1,DOF
          PHI_CHECK(K,I) = PHI_CORNER(K,I)
        ENDDO
      ENDDO
      IF (NCHECK.GT.3) THEN
        II = 4
        DO L = 1,3
          DO I = 1,NEGP
            DO K = 1,DOF
              PHI_CHECK(K,II) = PHI_EDGE(K,I,L)
            ENDDO
            II = II + 1
          ENDDO
        ENDDO
      ENDIF
      
C.....Integrate the basis functions

      PHI_INTEGRATED = 0.D0
      DO I = 1,NAGP
        DO K = 1,DOF
          PHI_INTEGRATED(K) = PHI_INTEGRATED(K) + WAGP(I)*PHI_AREA(K,I)
        ENDDO
      ENDDO

C.....If computing sediment transport call the subroutine prep_sed

      IF (SEDFLAG.GE.1) CALL PREP_SED()
        
C.....Initialize the DG.63 output file

      IF (ABS(NOUTGE).EQ.1) THEN
        OPEN(631,FILE=DIRNAME//'/'//'DG.63')
        WRITE(631,3220) RUNDES, RUNID, AGRID
        WRITE(631,3645) NDSETSE, DOF, DTDP*NSPOOLGE, NSPOOLGE, 1
c        WRITE(631,2120) TIME_A, IT
c2120          FORMAT(2X,E20.10,5X,I10)
        DO J = 1,NE
           DO K = 1,DOF
c                  if (k.eq.1) then
c                     WRITE(631,*) J, ZE(K,J,1),
c     $                  func_H(centroid_x(j),centroid_y(j),timedg)-2.0;
c                  else
                     WRITE(631,*) J, ZE(K,J,1)
c                  endif
           ENDDO
        ENDDO
      ENDIF

C.....Initialize the DG.64 output file

      IF (ABS(NOUTGV).EQ.1) THEN
        OPEN(641,FILE=DIRNAME//'/'//'DG.64')
        WRITE(641,3220) RUNDES, RUNID, AGRID
        WRITE(641,3645) NDSETSV, DOF, DTDP*NSPOOLGV, NSPOOLGV, 2
        DO J = 1,NE
           DO K = 1,DOF
              WRITE(641,*) J, QX(K,J,1),QY(J,K,1)
           ENDDO
        ENDDO
      ENDIF

C.....Initialize the DG.65 output file (contains elemental statuses such
C.....as the wet/dry status.

      IF ((ABS(NOUTGE).EQ.1).AND.(NOLIFA.GE.2)) THEN
        OPEN(651,FILE=DIRNAME//'/'//'DG.65')
        WRITE(651,3220) RUNDES, RUNID, AGRID
        WRITE(651,3645) NDSETSE, DOF, DTDP*NSPOOLGE, NSPOOLGE, 1
      ENDIF
3220  FORMAT(1X,A32,2X,A24,2X,A24)
3645  FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)

C.....Set p back to original value if p = 0

      IF (P_0.NE.P) THEN
        DOF = 1
        P = 0
      ENDIF
      
C.....Compute basis functions at stations

      IF (NSTAE.GT.0) THEN   ! Elevation stations
        CALL ALLOC_STAE( NSTAE )
        DO I = 1,NSTAE
          CALL STA_BASIS( XEL(I), YEL(I),  NNE(I), PHI_STAE(:,I) )
        ENDDO
      ENDIF
      
      IF (NSTAV.GT.0) THEN   ! Velocity Stations
        CALL ALLOC_STAV( NSTAV )
        DO I = 1,NSTAV
          CALL STA_BASIS( XEV(I), YEV(I),  NNV(I), PHI_STAV(:,I) )
        ENDDO
      ENDIF

C.....Close files

      CLOSE(115)
      CLOSE(25)

      RETURN
      END SUBROUTINE PREP_DG

C***********************************************************************
C
C     SUBROUTINE RK_TIME()
C
C     This subroutine does preparatory stuff for DG
C
C     Written by Ethan Kubatko (03-07-2005)
C
C***********************************************************************

      SUBROUTINE RK_TIME()

      USE GLOBAL
      USE DG

      IMPLICIT NONE

      INTEGER L
      REAL(SZ) ARK, BRK, CASUM, MAX_BOA

C.....Allocate the time stepping arrays

      NRK = RK_STAGE
      CALL ALLOC_RK()


C.....The forward Euler method

      IF ((RK_STAGE.EQ.1).AND.(RK_ORDER.EQ.1)) THEN

        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0

        ATVD(1,1) = 1.D0
        BTVD(1,1) = 1.D0

C.....SSP(s,2) schemes

      ELSEIF (RK_ORDER.EQ.2) THEN

        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0

        DO I = 1,NRK
          DO J = 0,NRK-1

            IF ((J.EQ.(I-1)).AND.(I.LT.NRK)) THEN
              ATVD(I,J+1) = 1.D0
              BTVD(I,J+1) = 1.D0/(NRK-1)
            ELSEIF ((J.EQ.0).AND.(I.EQ.NRK)) THEN
              ATVD(I,J+1) = 1.D0/NRK
            ELSEIF ((J.EQ.(NRK-1)).AND.(I.EQ.NRK)) THEN
              ATVD(I,J+1) = (NRK-1.D0)/NRK
              BTVD(I,J+1) = 1.D0/NRK
            ENDIF

          ENDDO
        ENDDO

C.....SSP(3,3) scheme

      ELSEIF ((RK_STAGE.EQ.3).AND.(RK_ORDER.EQ.3)) THEN

        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0

        ATVD(1,1) = 1.D0
        ATVD(2,1) = 3.D0/4.D0
        ATVD(2,2) = 1.D0/4.D0
        ATVD(3,1) = 1.D0/3.D0
        ATVD(3,3) = 2.D0/3.D0

        BTVD(1,1) = 1.D0
        BTVD(2,2) = 1.D0/4.D0
        BTVD(3,3) = 2.D0/3.D0

C.....SSP(4,3) scheme

      ELSEIF ((RK_STAGE.EQ.4).AND.(RK_ORDER.EQ.3)) THEN

        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0

        ATVD(1,1) = 1.D0
        ATVD(2,2) = 1.D0
        ATVD(3,1) = 2.D0/3.D0
        ATVD(3,3) = 1.D0/3.D0
        ATVD(4,4) = 1.D0

        BTVD(1,1) = 1.D0/2.D0
        BTVD(2,2) = 1.D0/2.D0
        BTVD(3,3) = 1.D0/6.D0
        BTVD(4,4) = 1.D0/2.D0

C.....SSP(5,3) scheme

      ELSEIF ((RK_STAGE.EQ.5).AND.(RK_ORDER.EQ.3)) THEN

        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0

        ATVD(1,1) = 1.D0
        ATVD(2,2) = 1.D0
        ATVD(3,1) = 0.355909775063327D0
        ATVD(3,3) = 0.644090224936674D0
        ATVD(4,1) = 0.367933791638137D0
        ATVD(4,4) = 0.632066208361863D0
        ATVD(5,3) = 0.237593836598569D0
        ATVD(5,5) = 0.762406163401431D0

        BTVD(1,1) = 0.377268915331368D0
        BTVD(2,2) = 0.377268915331368D0
        BTVD(3,3) = 0.242995220537396D0
        BTVD(4,4) = 0.238458932846290D0
        BTVD(5,5) = 0.287632146308408D0

C.....SSP(6,3) scheme

      ELSEIF ((RK_STAGE.EQ.6).AND.(RK_ORDER.EQ.3)) THEN

        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0

        ATVD(1,1) = 1.D0
        ATVD(2,2) = 1.D0
        ATVD(3,3) = 1.D0
        ATVD(4,1) = 0.476769811285196D0
        ATVD(4,2) = 0.098511733286064D0
        ATVD(4,4) = 0.424718455428740D0
        ATVD(5,5) = 1.D0
        ATVD(6,3) = 0.155221702560091D0
        ATVD(6,6) = 0.844778297439909D0

        BTVD(1,1) = 0.284220721334261D0
        BTVD(2,2) = 0.284220721334261D0
        BTVD(3,3) = 0.284220721334261D0
        BTVD(4,4) = 0.120713785765930D0
        BTVD(5,5) = 0.284220721334261D0
        BTVD(6,6) = 0.240103497065900D0

C.....SSP(7,3) scheme

      ELSEIF ((RK_STAGE.EQ.7).AND.(RK_ORDER.EQ.3)) THEN

        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0

        ATVD(1,1) = 1.D0
        ATVD(2,2) = 1.D0
        ATVD(3,3) = 1.D0
        ATVD(4,1) = 0.184962588071072D0
        ATVD(4,4) = 0.815037411928928D0
        ATVD(5,1) = 0.180718656570380D0
        ATVD(5,2) = 0.314831034403793D0
        ATVD(5,5) = 0.504450309025826D0
        ATVD(6,6) = 1.D0
        ATVD(7,4) = 0.120199000000000D0
        ATVD(7,7) = 0.879801000000000D0

        BTVD(1,1) = 0.233213863663009D0
        BTVD(2,2) = 0.233213863663009D0
        BTVD(3,3) = 0.233213863663009D0
        BTVD(4,4) = 0.190078023865845D0
        BTVD(5,5) = 0.117644805593912D0
        BTVD(6,6) = 0.233213863663009D0
        BTVD(7,7) = 0.205181790464579D0

C.....SSP(8,3) scheme

      ELSEIF ((RK_STAGE.EQ.8).AND.(RK_ORDER.EQ.3)) THEN

        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0

        ATVD(1,1) = 1.D0
        ATVD(2,2) = 1.D0
        ATVD(3,3) = 1.D0
        ATVD(4,4) = 1.D0
        ATVD(5,1) = 0.421366967085359D0
        ATVD(5,2) = 0.005949401107575D0
        ATVD(5,5) = 0.572683631807067D0
        ATVD(6,2) = 0.004254010666365D0
        ATVD(6,6) = 0.995745989333635D0
        ATVD(7,3) = 0.104380143093325D0
        ATVD(7,4) = 0.243265240906726D0
        ATVD(7,7) = 0.652354615999950D0
        ATVD(8,8) = 1.D0

        BTVD(1,1) = 0.195804015330143D0
        BTVD(2,2) = 0.195804015330143D0
        BTVD(3,3) = 0.195804015330143D0
        BTVD(4,4) = 0.195804015330143D0
        BTVD(5,5) = 0.112133754621673D0
        BTVD(6,6) = 0.194971062960412D0
        BTVD(7,7) = 0.127733653231944D0
        BTVD(8,8) = 0.195804015330143D0

C.....SSP(5,4) scheme

      ELSEIF ((RK_STAGE.EQ.5).AND.(RK_ORDER.EQ.4)) THEN

        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0

        ATVD(1,1) = 1.D0
        ATVD(2,1) = 0.44437049406734D0
        ATVD(2,2) = 0.55562950593266D0
        ATVD(3,1) = 0.62010185138540D0
        ATVD(3,3) = 0.37989814861460D0
        ATVD(4,1) = 0.17807995410773D0
        ATVD(4,4) = 0.82192004589227D0
        ATVD(5,1) = 0.00683325884039D0
        ATVD(5,3) = 0.51723167208978D0
        ATVD(5,4) = 0.12759831133288D0
        ATVD(5,5) = 0.34833675773694D0

        BTVD(1,1) = 0.39175222700392D0
        BTVD(2,2) = 0.36841059262959D0
        BTVD(3,3) = 0.25189177424738D0
        BTVD(4,4) = 0.54497475021237D0
        BTVD(5,4) = 0.08460416338212D0
        BTVD(5,5) = 0.22600748319395D0

C.....SSP(6,4) scheme

      ELSEIF ((RK_STAGE.EQ.6).AND.(RK_ORDER.EQ.4)) THEN

        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0

        ATVD(1,1) = 1.00000000000000D0
        ATVD(2,1) = 0.30948026455053D0
        ATVD(2,2) = 0.69051973544947D0
        ATVD(3,1) = 0.54205244285557D0
        ATVD(3,3) = 0.45794755714443D0
        ATVD(4,1) = 0.35984960863377D0
        ATVD(4,4) = 0.64015039136623D0
        ATVD(5,5) = 1.00000000000000D0
        ATVD(6,1) = 0.05776282890116D0
        ATVD(6,3) = 0.44216432622405D0
        ATVD(6,5) = 0.10115567086469D0
        ATVD(6,6) = 0.39891717401009D0

        BTVD(1,1) = 0.39270746575722D0
        BTVD(2,2) = 0.30154043149172D0
        BTVD(3,3) = 0.19997937335132D0
        BTVD(4,4) = 0.27954483459696D0
        BTVD(5,5) = 0.43668618869443D0
        BTVD(6,3) = 0.09150931531680D0
        BTVD(6,5) = 0.04417328437472D0
        BTVD(6,6) = 0.14911300530736D0

C.....SSP(7,4) scheme

      ELSEIF ((RK_STAGE.EQ.7).AND.(RK_ORDER.EQ.4)) THEN

        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0

        ATVD(1,1) = 1.D0
        ATVD(2,1) = 0.20161507213829D0
        ATVD(2,2) = 0.79838492786171D0
        ATVD(3,1) = 0.19469598207921D0
        ATVD(3,3) = 0.80530401792079D0
        ATVD(4,1) = 0.58143386885601D0
        ATVD(4,4) = 0.41856613114399D0
        ATVD(5,1) = 0.01934367892154D0
        ATVD(5,5) = 0.98065632107846D0
        ATVD(6,6) = 1.D0
        ATVD(7,1) = 0.06006304558847D0
        ATVD(7,3) = 0.30152730794242D0
        ATVD(7,4) = 0.10518998496676D0
        ATVD(7,5) = 0.01483791154585D0
        ATVD(7,7) = 0.51838174995650D0

        BTVD(1,1) = 0.30111872706068D0
        BTVD(2,2) = 0.24040865318216D0
        BTVD(3,3) = 0.24249212077315D0
        BTVD(4,4) = 0.12603810060080D0
        BTVD(5,5) = 0.29529398308716D0
        BTVD(6,6) = 0.30111872706068D0
        BTVD(7,3) = 0.09079551914158D0
        BTVD(7,4) = 0.02888359354880D0
        BTVD(7,7) = 0.15609445267839D0

C.....SSP(8,4) scheme

      ELSEIF ((RK_STAGE.EQ.8).AND.(RK_ORDER.EQ.4)) THEN

        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0

        ATVD(1,1) = 1.D0
        ATVD(2,1) = 0.10645325745007D0
        ATVD(2,2) = 0.89354674254993D0
        ATVD(3,3) = 1.D0
        ATVD(4,1) = 0.57175518477257D0
        ATVD(4,4) = 0.42824481522743D0
        ATVD(5,1) = 0.19161667219044D0
        ATVD(5,5) = 0.80838332780956D0
        ATVD(6,6) = 1.D0
        ATVD(7,7) = 1.D0
        ATVD(8,1) = 0.02580435327923D0
        ATVD(8,3) = 0.03629901341774D0
        ATVD(8,4) = 0.31859181340256D0
        ATVD(8,5) = 0.05186768980103D0
        ATVD(8,6) = 0.03944076217320D0
        ATVD(8,7) = 0.00511633747411D0
        ATVD(8,8) = 0.52288003045213D0

        BTVD(1,1) = 0.24120020561311D0
        BTVD(2,2) = 0.21552365802797D0
        BTVD(3,3) = 0.24120020561311D0
        BTVD(4,4) = 0.10329273748560D0
        BTVD(5,5) = 0.19498222488188D0
        BTVD(6,6) = 0.24120020561311D0
        BTVD(7,7) = 0.24120020561311D0
        BTVD(8,3) = 0.00875532949991D0
        BTVD(8,4) = 0.06195575835101D0
        BTVD(8,6) = 0.00951311994571D0
        BTVD(8,8) = 0.12611877085604D0

      ENDIF

C.....Compute the time dependent parameters

      DO K = 0,RK_STAGE-1
        DO I = 1,RK_STAGE
          CASUM = 0.D0
          DO L = K+1,I-1
            CASUM = CASUM + CTVD(L,K+1)*ATVD(I,L+1)
          ENDDO
          CTVD(I,K+1) = BTVD(I,K+1) + CASUM
        ENDDO
      ENDDO

      DO K = 1,RK_STAGE-1
        DTVD(K+1) = 0.D0
        DO L = 0,K-1
          DTVD(K+1) = DTVD(K+1) + CTVD(K,L+1)
        ENDDO
      ENDDO

C.....Compute the maximum beta over alpha ratio at each stage

      DO IRK = 1,RK_STAGE
        MAX_BOA = 0.D0
        DO I = 1,IRK
          ARK = ATVD(IRK,I)
          BRK = BTVD(IRK,I)
          IF (ARK.NE.0.D0) THEN
            IF (MAX_BOA.LT.BRK/ARK) MAX_BOA = BRK/ARK
          ENDIF
        ENDDO
        MAX_BOA_DT(IRK) = MAX_BOA*DT
      ENDDO

      DEALLOCATE(CTVD)
      RETURN
      END SUBROUTINE RK_TIME
