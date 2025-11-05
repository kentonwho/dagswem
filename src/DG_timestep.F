C***********************************************************************
C
C     SUBROUTINE DG_TIMESTEP(IT)
C
C     This is the main subroutine for the DG hydro
C
C     Written by Shintaro Bunya (01-01-2007)
C
C***********************************************************************

      SUBROUTINE DG_TIMESTEP(IT)

C.....Use appropriate modules
      
      USE SIZES, ONLY : SZ
      USE GLOBAL, ONLY : SEDFLAG, TIME_A, TIMEH,NE,NWS,NTIP
      USE SED, ONLY : MIT
      USE DG !, ONLY : NRK,DOF,QX,QY,ZE,RHS_QX,RHS_QY,RHS_ZE,tn_h,tn_a,tnkp1_a,tnkp1_h,Ndom,DtDom
#ifdef CMPI
      USE MESSENGER_NEW
#endif

      IMPLICIT NONE

C.....Declare local variables

      INTEGER i,j,k,kk,IT,lvl ! cjt
      
C.....Hydrodynamics CJT ******************************************************
      tn_a(:) = time_a
      tn_h(:) = timeh

      ! Set the next time level on the domain with the largest dt
      tnkp1_a(Ndom) = tn_a(Ndom) + DtDom(Ndom)
      tnkp1_h(Ndom) = tn_h(Ndom) + DtDom(Ndom)      

C.......Obtain the meteorological forcing
      IF (NWS.NE.0) CALL MET_FORCING(IT)
        
C.......Compute tidal potential terms
      IF (NTIP.NE.0) CALL TIDAL_POTENTIAL()

      IF (Ndom>1) CALL FIRST_LOOP(IT,Ndom-1)

      CALL DG_HYDRO_TIMESTEP(IT,Ndom)

C.....Sediment transport *****************************************************

c      IF ((SEDFLAG.GE.1).AND.(MIT.EQ.IT)) THEN
c        CALL DG_SED_TIMESTEP(IT)
c      ENDIF

C.....Get Ready for next ts  CJT *********************************************

      time_a = time_a + DtDom(Ndom)
      
!      DO J = 1,MNE
!        DO K = 1,DOF
!          ZE(K,J,1) = ZE(K,J,NRK+1)
!          QX(K,J,1) = QX(K,J,NRK+1)
!          QY(K,J,1) = QY(K,J,NRK+1)
!        ENDDO
!      ENDDO
!
!      DO KK = 2,NRK+1
!        DO J = 1,MNE
!          DO K = 1,DOF
!            ZE(K,J,KK) = 0.D0
!            QX(K,J,KK) = 0.D0
!            QY(K,J,KK) = 0.D0
!            RHS_ZE(K,J,KK-1) = 0.D0
!            RHS_QX(K,J,KK-1) = 0.D0
!            RHS_QY(K,J,KK-1) = 0.D0
!          ENDDO
!        ENDDO
!      ENDDO



C.....Write out results

#ifdef CMPI
#ifdef BLKOUT
      CALL MSG_BLOCKSYNC_START()
#endif
#endif

      CALL WRITE_RESULTS(IT,.FALSE.)
#ifdef CMPI
#ifdef BLKOUT
      CALL MSG_BLOCKSYNC_FINISH()
#endif
#endif

      CALL SCRUTINIZE_SOLUTION(IT)

      RETURN
      END SUBROUTINE


C***********************************************************************
C
C     SUBROUTINE SCRUTINIZE_SOLUTION(IT)
C
C     Detect a negative water depth and a NaN in solution arrays and 
C     terminate if either of them is found.
C
C     Written by Shintaro Bunya (01-01-2007)
C
C***********************************************************************

      SUBROUTINE SCRUTINIZE_SOLUTION(IT)

C.....Use appropriate modules
      
#ifdef CMPI
      USE SIZES, ONLY : SZ,MYPROC
      USE GLOBAL, ONLY : TIME_A
      USE MESSENGER_NEW, ONLY:  MESSAGE_FINI
C      USE MESSENGER_ELEM
      USE MESSENGER_node_new
#else
      USE SIZES, ONLY : SZ
#endif
      USE GLOBAL, ONLY : NE,NM,DP,X,Y,SLAM,SFEA,DEG2RAD
      USE DG, ONLY : DOF,ZE,QX,QY,TsGroup,WDFLG

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER IT,J,K
      LOGICAL Detected
      REAL(SZ) DPAVG
C
      INTEGER ErrorElevExceeded

C.....Detect negative mass
      ErrorElevExceeded=0
      DO J = 1,NE
         DPAVG = (DP(NM(J,1))+DP(NM(J,2))+DP(NM(J,3)))/3.D0
         Detected = .FALSE.
         IF((ZE(1,J,1)+DPAVG).LE.0.D0) THEN
            PRINT *, ''
#ifdef CMPI
            PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
            PRINT *, '  ZE(',1,',',J,',1) = ',ZE(1,J,1)
            PRINT *, '  HE(',1,',',J,',1) = ',ZE(1,J,1)+DPAVG
            PRINT *, '  QX(',1,',',J,',1) = ',QX(1,J,1)
            PRINT *, '  QY(',1,',',J,',1) = ',QY(1,J,1)
            PRINT *, '  DP = ',DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3))
            PRINT *, '  DPAVG = ',DPAVG
            PRINT *, ' x,  y ',SLAM(NM(J,1))/DEG2RAD,
     $           SFEA(NM(J,1))/DEG2RAD
            Detected = .TRUE.
            ErrorElevExceeded=1            
         ENDIF
      ENDDO
#ifdef CMPI
      CALL ErrorElevSum(ErrorElevExceeded)
      IF(ErrorElevExceeded.NE.0) THEN
         if (myproc.eq.0) then
            PRINT *, ''
            PRINT *,'  PROGRAM WILL BE TERMINATED in DG_timestep'
            PRINT *, ''
         endif

         call message_fini ()
         STOP
      endif
#else
c$$$      IF(ErrorElevExceeded.NE.0) THEN
c$$$         PRINT *, ''
c$$$         PRINT *,'  PROGRAM WILL BE TERMINATED in DG_timestep'
c$$$         PRINT *, ''
c$$$         stop
c$$$      ENDIF
#endif


C.....Detect NaN
      ErrorElevExceeded=0
C      
      DO J = 1,NE
        DO K = 1,DOF
          Detected = .FALSE.
          IF(ZE(K,J,1).NE.ZE(K,J,1)) THEN
            PRINT *, ''
#ifdef CMPI
            PRINT *, 'IN SUBDOMAIN ',MYPROC,time_a,it
#endif
            PRINT *, '  ZE(',K,',',J,',1) IS ',ZE(K,J,1)
            PRINT *, ' x,  y ',SLAM(NM(J,1))/DEG2RAD,
     $           SFEA(NM(J,1))/DEG2RAD,dp(NM(j,1)),TsGroup(J),WDFLG(J)
            Detected = .TRUE.
            ErrorElevExceeded = 1
          ENDIF

          IF(QX(K,J,1).NE.QX(K,J,1)) THEN
#ifdef CMPI
            PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
            PRINT *, ''
            PRINT *, '  QX(',K,',',J,',1) IS ',QX(K,J,1)
            PRINT *, ' x,  y ',SLAM(NM(J,1))/DEG2RAD,
     $           SFEA(NM(J,1))/DEG2RAD,TsGroup(J),WDFLG(J)
            Detected = .TRUE.
            ErrorElevExceeded = 1
          ENDIF
          IF(QY(K,J,1).NE.QY(K,J,1)) THEN
#ifdef CMPI
            PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
            PRINT *, ''
            PRINT *, '  QY(',K,',',J,',1) IS ',QY(K,J,1)
            PRINT *, ' x,  y ',SLAM(NM(J,1))/DEG2RAD,
     $           SFEA(NM(J,1))/DEG2RAD,TsGroup(J),WDFLG(J)
            Detected = .TRUE.
            ErrorElevExceeded = 1
          ENDIF

!          IF(Detected) THEN
!            PRINT *, ''
!            PRINT *,'  WRITING OUT THE LAST-MOMENT RESULTS...'
!            PRINT *, ''
!            CALL WRITE_RESULTS(IT,.TRUE.)
!            PRINT *, ''
!            PRINT *,'  PROGRAM WILL BE TERMINATED.'
!            PRINT *, ''
!#ifdef CMPI
!            call message_fini()
!            STOP
!#else
!            stop
!#endif
!          ENDIF
        ENDDO
      ENDDO
#ifdef CMPI
      CALL ErrorElevSum(ErrorElevExceeded)
c      write(*,*) myproc,ErrorElevExceeded
      IF (ErrorElevExceeded.NE.0) THEN
        CALL MESSAGE_FINI()
         stop
      ENDIF
#else
!      write(*,*)ErrorElevExceeded
c$$$      IF (ErrorElevExceeded.NE.0) THEN
c$$$         STOP
c$$$      ENDIF
#endif

      RETURN
      END SUBROUTINE


