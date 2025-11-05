C***********************************************************************
C
C     SUBROUTINE SED_PREP()
C
C     This subroutine reads in, sets, or computes the various time
C     independent sediment transport parameters for the Lund formula.
C
C     An effort has been made to follow the notation and "flow" (no pun
C     intended) of Larson's clsed_new.f and sedsim.f codes.
C
C     Written by Ethan Kubatko (03-05-2005)
C
C-----------------------------------------------------------------------
C
C     Modification history on hp_DG_ADCIRC after v9_sb2.2.1
C     
C     v9_sb2.2.1 - 08/16/05 - sb - Hard bottom implementation
C
C***********************************************************************

      SUBROUTINE PREP_SED()

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE
      
C.....Read in the morpholigical time step

      READ(25,*) MORPH
      
C.....Read in the sediment transport magnification factor

      READ(25,*) MAG_SED
      
C.....Read in the median grain size (in mm)

      READ(25,*) D50

C.....Convert median grain size from mm to m

      D50 = D50/1000.D0
      
C.....Read in the porosity
      
      READ(25,*) POROSITY

C.....Read in sediment and water density

      READ(25,*) RHOSED0
      READ(25,*) RHOWAT0
      
C.....Compute relative density

      RHORELO = RHOSED0/RHOWAT0
      
C.....Read in the water temperature

      READ(25,*) TEMPC

Csb-HB..Read in the hard bottom flag

      READ(25,*) HBFLAG

C.....Compute the temperature dependent kinematic viscosity (taken from
C.....Magnus' viscosity.f program with added error message)

      IF ((TEMPC.LT.0).OR.(TEMPC.GT.40)) THEN
        WRITE(16,*) 'Your specification of TEMPC (a UNIT dg input ',
     &              'parameter) is out of range (0-40 C)'
        STOP
      ENDIF
      
      DO I = 1,8
         IF ((TEMPC.GE.5.0*REAL(I-1)).AND.(TEMPC.LE.5.0*REAL(I))) THEN
            TEMP1 = 5.D0*REAL(I-1)
            TEMP2 = 5.D0*REAL(I)
            VNY = NU(I) + (NU(I+1) - NU(I))*(TEMPC-TEMP1)/(TEMP2-TEMP1)
            VNY = VNY*0.000001D0
            GOTO 111
         ENDIF
      ENDDO
      
C.....Compute conversion factor for non-dimension to dimensional

111   SNORM = G*D50*(RHOSED0 - RHOWAT0)

C.....Set the empirical bed load constants

      A_EMP = 12.D0
      B_EMP = 4.5D0
      
C.....Set the various suspended load constants

      GAMMA0  = 0.0005D0
      CB      = 0.65D0
      CD      = 0.02D0
      CM      = 1.D0/7.D0
      ADE     = 0.07D0
      TETALIM = 1.D0
      
C.....Convert to appropriate dimensions

      TAULIM = TETALIM*SNORM
      ADEC   = ADE/SNORM
      
C.....Compute grain size parameter

      DS = D50*(G*(RHORELO - 1.D0)/VNY**2.D0)**(1.D0/3.D0)
      
C.....Compute the sediment fall speed parameter

      BPAR = DS**3

C.....Comput the fall velocity (Hallermeier's method)

      IF (BPAR.LT.39) THEN
        WFV = (RHORELO - 1.D0)*G*D50**2.D0/18.D0/VNY
      ELSEIF ((BPAR.GT.39).AND.(BPAR.LT.10000)) THEN
        WFV = ((RHORELO - 1.D0)*G)**0.7D0*D50**1.1D0/6.D0/VNY**0.4D0
      ELSE
        WFV = ((RHORELO - 1.D0)*G*D50/0.91D0)**0.5D0
      ENDIF
      
C.....Compute the non-dimensional critical shear stress
C.....(Soulsby and Whitehouse)

      TETACR = 0.3D0/(1.D0 + 1.2D0*DS) + 0.055D0*(1.D0 -EXP(-0.02D0*DS))

C.....Convert to dimensional shear stress

      TAUCR = SNORM*TETACR

C.....Compute the time independent bed load parameter

      CBL = A_EMP/(G*(RHOSED0 - RHOWAT0))
      
C.....Compute the time independent suspended load parameter

      CSL = 2.D0/PI*GAMMA0*CB
      
C.....Compute complete gamma function (asymptotic solution)

      FGL = (CM + 0.5D0)*LOG(CM + 1.D0) - (CM + 1.D0)
     &    + 0.5D0*LOG(2.D0*PI) + 1.D0/12.D0/(CM + 1.D0)
     &    - 1.D0/360.D0/(CM + 1.D0)**2.D0

      FG = EXP(FGL)
      
C.....Set friction factor paramter

      KAPPA = 0.41D0
      
C.....Compute the skin friction

      KS = 2.D0*D50

C.....Compute the roughness length

      Z0 = KS/30.D0

C.....Compute the porosity factor

      PORFAC = 1.D0/(1.D0 - POROSITY)
      
C.....Compute lower limit of the near-bed wave semi-excursion distance
C.....(based on an upper limit friction wave factor of 0.30 - this
C.....prevents blow-up of the friction wave factor for very small values
C.....of the near-bed wave semi-excursion distance)
      
      AW_LIM = KS*((LOG(0.30D0) + 5.98D0)/5.2D0)**(-1.D0/0.194D0)
      
C.....Likewise compute lower limit of the total energy dissipation,
C.....DTOT_MIN so that the sediment decay parameter, ADIFF does not
C.....exceed 20, i.e. if DTOT < DTOT_MIN then ADIFF = 20

      DTOT_MIN = RHOWAT0*(WFV/CD/20.D0)**3.D0
      
C.....Create the allocatable sediment arrays

      CALL ALLOC_SED1()

C.....Intialize variables

      RHS_HB(:,:,:) = 0.D0
      HB0(:,:,:) = HB(:,:,:)

Csb-HB Read in DP_HB (Nodal bathymetric depth of the hard bottom)
      IF(HBFLAG.GT.0) CALL READ_HB()


C.....Set the morphology timestep counter

      MIT  = ITHS + 1
      
C.....Check to make sure that the specification of NWS and SEDFLAG are
C.....consistent

        IF ((SEDFLAG.EQ.2).AND.(ABS(NWS).LT.100)) THEN
          WRITE(6,*) 'IF WAVE EFFECTS ARE INCLUDED IN THE SEDIMENT',
     &               'TRANSPORT CALCULATIONS (SEDFLAG = 2) THEN',
     &               '|NWS|>=100 - CHECK YOUR INPUT PARAMETERS!!!'
          STOP
        ENDIF

      RETURN
      END SUBROUTINE

      
Csb-HB Read in DP_HB (Nodal bathymetric depth of the hard bottom)
C.....Read input file 'fort.hb'
C.....The format of 'fort.hb' is the same as that of 'fort.14'
C.....but the part after the node info is ignored

      SUBROUTINE READ_HB()

      USE SIZES,ONLY : SZ
      USE GLOBAL,ONLY : NP,NSCREEN,MYPROC,DIRNAME
      USE SED,ONLY : DP_HB

      IMPLICIT NONE

      INTEGER DUM1,DUM2,I,JKI
      REAL(SZ) RDUM1,RDUM2
      

      OPEN(26,FILE=DIRNAME//'/'//'fort.hb')

      READ(26,*)                ! SKIP AGRID
      READ(26,*) DUM1,DUM2      ! NE,NP

      IF(DUM2.NE.NP) THEN
        IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99900)
        WRITE(16,99900)
99900 FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &       //,1X,'NP IN FORT.HB IS DIFFERENT FROM THAT IN FORT.14',
     &       //,1X,'!!!!!! EXECUTION WILL BE TERMINATED !!!!!!',//)
        STOP
      ENDIF

      DO I=1,NP
        READ(26,*) JKI,RDUM1,RDUM2,DP_HB(I)
        IF(JKI.NE.I) THEN
          IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99901)
          WRITE(16,99901)
99901 FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &       'INPUT ERROR  !!!!!!!!!',
     &       //,1X,'YOUR NODE NUMBERING IS NOT SEQUENTIAL ',
     &       /,1X,'CHECK YOUR FORT.HB FILE CAREFULLY',//)
        ENDIF
      ENDDO

      CLOSE(26)

      END SUBROUTINE
