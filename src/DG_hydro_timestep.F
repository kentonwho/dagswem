C***********************************************************************
C
C     SUBROUTINE DG_HYDRO_TIMESTEP(IT)
C
C     This is the main subroutine for the DG hydro
C
C     Written by Ethan Kubatko (06-11-2004)
C
C-----------------------------------------------------------------------
C
C  Modification history on hp_DG_ADCIRC since v9_sb1
C
C     v??     - 08/05    - sb - parallelized
C     v??     - 08/05    - sb - wetting/drying
C     v10_sb5 - 10/05    - sb - wetting/drying algorithm is improved
C     v10_sb5 - 10/05    - sb - Ethan's slope limiter is consolidated
C
C***********************************************************************

      SUBROUTINE DG_HYDRO_TIMESTEP(IT,lvl)
      USE GLOBAL
      USE DG
      USE SED

#ifdef CMPI
      USE MESSENGER_NEW
#endif

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER IT,L,GED,NBOREL,NNBORS,NDRYNBORS,Istop,lvl,elem ! cjt
      INTEGER Detected
      REAL(SZ) QBCT1,QBCT2,QBC1,QBC2,ZP(3),QXP(3),QXP_OLD(3),
     $     QYP(3),QYP_OLD(3),ZP_OLD(3)
      REAL(SZ) DPAVG, ITDT, ARK, BRK, DTl
      real*4 tarray1(2),tarray2(2)

      ! CJT
c      WRITE(99,*) "Taking one of ",Msteps(lvl)," steps with all 
c     & elements of level: ",lvl," time: ",tnkp1_a(lvl)

C.....Assign ts to local value CJT
      DTl = DtDom(lvl)

C.....Update ETA1 (shintaro: do we need this?)

      DO I = 1,MNP
        ETA1(I) = ETA2(I)
      ENDDO

C.....Begin Runge-Kutta time stepper

      DO 100 IRK = 1,NRK
      
C.......Compute the DG time and DG ramp

!        TIMEDG = tn_a(lvl) - DTl + DTVD(IRK)*DTl
        TIMEDG = tn_a(lvl) + DTVD(IRK)*DTl
        RAMPDG = 1.D0
!        IF ((NRAMP.EQ.1).AND.(IT.LE.IRAMPING)) THEN
!          RAMPDG = TANH((2.D0*(IT + DTVD(IRK))*DTl/86400.D0)/DRAMP)
        IF ((NRAMP.EQ.1).AND.(TIMEDG.LE.DRAMP*86400.d0)) THEN
          RAMPDG = TANH((2.D0*(TIMEDG + DTVD(IRK)*DTl)/86400.D0)/DRAMP)
        ENDIF
        RAMP = RAMPDG
        
C.......Obtain the meteorological forcing
c        IF (NWS.NE.0) CALL MET_FORCING(IT)
        
C.......Compute tidal potential terms
c        IF (NTIP.NE.0) CALL TIDAL_POTENTIAL()

C.......Compute LDG auxiliary equations
        IF (EVMSUM.NE.0.D0) CALL LDG_HYDRO(IT)
     
C.......Compute elevation specified edges
c        write(*,*) 'call ocean edge hydro',lvl,irk
        IF (NEEDS.GT.0)  CALL OCEAN_EDGE_HYDRO(lvl)

C.......Compute no-normal flow edges
c        write(*,*) 'nleds', nleds
        IF (NLEDS.GT.0)  CALL LAND_EDGE_HYDRO(lvl)

C.......Compute non-zero flow edges
        IF (NFEDS.GT.0)  CALL FLOW_EDGE_HYDRO(lvl)
        
C.......Compute radiation edges
        IF (NREDS.GT.0)  CALL RADIATION_EDGE_HYDRO(lvl)

C.......Compute internal barrier edges
c        if (lts_flag.eq.1) then
c        if (myproc.eq.37) then
c           write(*,*) 'nibeds ',nibeds
c        endif
c      endif
        IF (NIBEDS.GT.0) CALL IBARRIER_EDGE_HYDRO(IT,lvl)
c        if (myproc.eq.37) write(*,*) 'back from ibarrier '
c
C.......Compute external barrier edges
c        IF (NIBEDS.GT.0) CALL EBARRIER_EDGE_HYDRO()

C.......Compute internal edges
c        if (lts_flag.eq.1) then
c        if (myproc.eq.37.and.lvl.eq.1) then
c           write(37,*) 'before ie',timeDG,rhs_qx(1,20355,irk),
c     $          rhs_qx(2,20355,irk),rhs_qx(3,20355,irk)
c        endif
c        else
c           if (myproc.eq.133) then
c              write(133,*) 'before ie',timeDG,rhs_qx(1,17511,irk),
c     $          rhs_qx(2,17511,irk),rhs_qx(3,17511,irk)
c           endif
c        endif
        call etime(tarray1)
        CALL INTERNAL_EDGE_HYDRO(lvl)
        call etime(tarray2)
        TIMINGS(1)=TIMINGS(1)+tarray2(1)+tarray2(2)-(tarray1(1)+tarray1(2))

c        if (lts_flag.eq.1) then
c           if (myproc.eq.37.and.lvl.eq.1) then
c              write(37,*) 'after ie',timeDG,rhs_qx(1,20355,irk),
c     $             rhs_qx(2,20355,irk)
c           endif
c        else
c           if (myproc.eq.133) then
c              write(133,*) 'after  ie',timeDG,rhs_qx(1,17511,irk),
c     $          rhs_qx(2,17511,irk),rhs_qx(3,17511,irk)
c           endif
c        endif

C.......Compute elements to finish building the rhs
c        write(*,*) 'irk, timedge ',irk,timedg
        call etime(tarray1)
        CALL RHS_DG_HYDRO(lvl)
        call etime(tarray2)
        TIMINGS(2)=TIMINGS(2)+tarray2(1)+tarray2(2)-(tarray1(1)+tarray1(2))
c        if (lts_flag.eq.1) then
c        if (myproc.eq.37.and.lvl.eq.1) then
c           write(37,*) 'after rhs ',timeDG,rhs_qx(1,20355,irk),
c     $          rhs_qx(2,20355,irk)
c        endif
c        else
c        if (myproc.eq.133) then
c           write(133,*) 'after rhs ',timeDG,rhs_qx(1,17511,irk),
c     $          rhs_qx(2,17511,irk)
c        endif
c        endif
           

C.......SSP Runge-Kutta Time Scheme

        DO I = 1,IRK
          ARK = ATVD(IRK,I)
          BRK = BTVD(IRK,I)*DTl
          DO J = 1,Nelem(lvl)
            elem = ElemNum(lvl,J)
            DO K = 1,DOF
              ZE(K,elem,IRK+1) = ZE(K,elem,IRK+1) + ARK*ZE(K,elem,I)
     &                                      + BRK*RHS_ZE(K,elem,I)
              QX(K,elem,IRK+1) = QX(K,elem,IRK+1) + ARK*QX(K,elem,I)
     &                                      + BRK*RHS_QX(K,elem,I)
              QY(K,elem,IRK+1) = QY(K,elem,IRK+1) + ARK*QY(K,elem,I)
     &                                      + BRK*RHS_QY(K,elem,I)
            ENDDO
          ENDDO
        ENDDO

C.......Corey -  Debug CJT
        ! Check to make sure RHS one element up is being added to from internal edge flux
c        IF (lvl<Ndom) THEN
c          DO J = 1,Nelem(lvl+1)
c		elem = ElemNum(lvl+1,J)
c		DO K = 1,DOF
c			WRITE(101,'(2i10,3f20.10)') elem,TsGroup(elem),RHS_ZE(1,elem,1),RHS_ZE(2,elem,1),RHS_ZE(3,elem,1)
c		ENDDO
c          ENDDO
c        ENDIF
c        WRITE(101,*) "-------------------------------------------------"


#ifdef CMPI
        CALL UPDATER_ELEM_MOD(ZE,IRK+1,lvl)
        CALL UPDATER_ELEM_MOD(QX,IRK+1,lvl)
        CALL UPDATER_ELEM_MOD(QY,IRK+1,lvl)
#endif


C.......Apply the slopelimiter if appropriate

        IF (SLOPEFLAG.NE.0) then
           call etime(tarray1)
           CALL SLOPELIMITER(lvl)
           call etime(tarray2)
           TIMINGS(3)=TIMINGS(3)
     $          +tarray2(1)+tarray2(2)-(tarray1(1)+tarray1(2))
        endif

C.......Apply the wet-dry algorithm if appropriate

        IF (NOLIFA.GE.2) CALL WETDRY(lvl)


C.......For parallel run update dofs if slope limiter and/or wetting andvdrying is being used
#ifdef CMPI
        IF (SLOPEFLAG.NE.0.OR.NOLIFA.GE.2) THEN
          CALL UPDATER_ELEM_MOD(ZE,IRK+1,lvl)
          CALL UPDATER_ELEM_MOD(QX,IRK+1,lvl)
          CALL UPDATER_ELEM_MOD(QY,IRK+1,lvl)
        ENDIF
#endif


100   CONTINUE

C.....Prepare for next ts
      DO J = 1,Nelem(lvl)
        elem = ElemNum(lvl,J)
        DO K = 1,DOF
          ZE(K,elem,1) = ZE(K,elem,NRK+1)
          QX(K,elem,1) = QX(K,elem,NRK+1)
          QY(K,elem,1) = QY(K,elem,NRK+1)
        ENDDO
      ENDDO

      DO KK = 2,NRK+1
        DO J = 1,Nelem(lvl)
          elem = ElemNum(lvl,J)
          DO K = 1,DOF
            ZE(K,elem,KK) = 0.0d0
            QX(K,elem,KK) = 0.0d0
            QY(K,elem,KK) = 0.0d0
            RHS_ZE(K,elem,KK-1) = 0.D0
            RHS_QX(K,elem,KK-1) = 0.D0
            RHS_QY(K,elem,KK-1) = 0.D0
          ENDDO
        ENDDO
      ENDDO



      RETURN
      END SUBROUTINE
