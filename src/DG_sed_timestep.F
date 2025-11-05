C***********************************************************************
C
C     SUBROUTINE DG_SED_TIMESTEP()
C
C     This is the main subroutine for the DG sed
C
C     Written by Ethan Kubatko (03-05-2005)
C
C     Modifications for parallel runs, Shintaro Bunya, Aug 2005 
C
C***********************************************************************

      SUBROUTINE DG_SED_TIMESTEP(IT)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED
Csb-PDG1
#ifdef CMPI
      USE MESSENGER_NEW
#endif
C--

      IMPLICIT NONE

C.....Declare local variables

      INTEGER IT, L, GED
      
C.....Write out header

#ifdef CMPI
      IF(MYPROC.EQ.0) THEN
        WRITE(6,*) ''
        WRITE(6,*) '*** UPDATING THE BATHYMETRY ***'
        WRITE(6,*) ''
      ENDIF
#else
      WRITE(6,*) ''
      WRITE(6,*) '*** UPDATING THE BATHYMETRY ***'
      WRITE(6,*) ''
#endif

      DO 100 IRK = 1,NRK

C.....Loop over internal edges

        DO L = 1,NIEDS
          CALL INTERNAL_EDGE_SED(L)
        ENDDO

C.....Loop over elevation edges

        DO L = 1,NEEDS
          CALL OCEAN_EDGE_SED(L)
        ENDDO
        
C.....Loop over no normal flow edges

C        DO L = 1,NLEDS
C          CALL LAND_EDGE_SED(L)
C        ENDDO

C.....Loop over non-zero flow edges

        DO L = 1,NFEDS
          CALL FLOW_EDGE_SED(L)
        ENDDO

C.....Loop over elements to finish building the rhs

        DO L = 1,NE
          CALL RHS_DG_SED(L)
        ENDDO

C.....SSP Runge-Kutta Time Scheme

        DO I = 1,IRK
          DO J = 1,NE
            DO K = 1,DOF
              HB(K,J,IRK+1) = HB(K,J,IRK+1) + ATVD(IRK,I)*HB(K,J,I)
     &                  + BTVD(IRK,I)*(MORPH*DT)*(-PORFAC*RHS_HB(K,J,I))
     
            ENDDO
          ENDDO
        ENDDO

Csb-PDG2
c#ifdef CMPI
C.....Update values on ghost edges
c        CALL UPDATER_ELEM(HB,QX,QY,IRK+1,1)
c#endif
C--        

C.....Apply the slopelimiter

        IF (SLOPEFLAG.GT.0) THEN
C          DO L = 1,NLEDS
C            GED = NLEDN(L)
C            EL_IN = NEDEL(1,GED)
C            HB(2:DOF,EL_IN,IRK+1) = 0.D0
C          ENDDO
          HB(2:DOF,:,IRK+1) = 0.D0
        ENDIF

100   CONTINUE

C.....Update bathymetry for next time step

      DO J = 1,MNE
        DO K = 1,DOF

          HB(K,J,1) = HB(K,J,NRK+1)

          DO KK = 2,NRK+1

            HB(K,J,KK) = 0.D0

            RHS_HB(K,J,KK-1) = 0.D0

          ENDDO
        ENDDO
      ENDDO
      
C.....Increment morphological timestep counter

      MIT = MIT + MORPH

      RETURN
      END SUBROUTINE
