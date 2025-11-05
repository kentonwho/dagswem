C***********************************************************************
C
C     MODULE SED
C
C     This module contains various data declarations for the sediment
C     transport model.
C
C     Written by Ethan Kubatko (03-07-2005)
C
C***********************************************************************

      MODULE SED
      
C.....Use appropriate modules

      USE SIZES
      USE DG
      
C.....Declare integer variables

      INTEGER MIT

C.....Declare real variables

      REAL(SZ) A_EMP, ADE, ADEC, ADIFF, AW, AW_LIM, A_WAVE
      REAL(SZ) B_EMP, BPAR
      REAL(SZ) CB, CBL, CD, CEXP, CM, CR, CRM, CSL, CVEL
      REAL(SZ) D50, DB, DC, DOT_U_MAG, DS, DTOT, DTOT_MIN, DW
      REAL(SZ) FC, FG, FGL, FIG, F_WAVE, FW
      REAL(SZ) GAMMA0
      REAL(SZ) H_WAVE
      REAL(SZ) KAPPA, KS
      REAL(SZ) L_WAVE
      REAL(SZ) MAG_SED, MORPH
      REAL(SZ) PHI0, PORFAC, POROSITY
      REAL(SZ) QS_HAT, QSBED, QSSUS, QSTOT, QSX, QSX_IN, QSY, QSY_IN
      REAL(SZ) RHOSED0, RHORELO
      REAL(SZ) SNORM
      REAL(SZ) TAUC, TAUCR, TAUCW, TAUCWM, TAULIM, TAUW, TAUWM
      REAL(SZ) TEMP1, TEMP2, TEMPC, TETACR, TETALIM, T_WAVE
      REAL(SZ) UC, UMAG_EX, UMAG_IN, UW
      REAL(SZ) VNY
      REAL(SZ) WAVE_RATIO, WAVE_X, WAVE_Y, WFV
      REAL(SZ) Y_WAVE
      REAL(SZ) Z0
Csb-HB
      INTEGER HBFLAG                   ! =0...w/o hard bottom, =1...w/ hard bottom
      REAL(SZ),ALLOCATABLE :: DP_HB(:) ! Nodal bathymetric depth of the hard bottom
C--
      
C.....Declare real arrays

      REAL(SZ), DIMENSION(9) :: NU = (/ 1.790D0, 1.510D0, 1.310D0,
     &                                  1.140D0, 1.000D0, 0.894D0,
     &                                  0.799D0, 0.724D0, 0.658D0  /)
      
C.....Declare real allocatable arrays

      REAL(SZ), ALLOCATABLE :: RHS_HB(:,:,:), HB0(:,:,:)
      REAL(SZ), ALLOCATABLE :: Q_ROE(:,:)

C********************** END OF DATA DECLARATIONS ***********************

      CONTAINS

      SUBROUTINE ALLOC_SED1()
      ALLOCATE ( RHS_HB(DOF,MNE,NRK), HB0(DOF,MNE,1) )
Csb-HB
      IF(HBFLAG.GT.0) ALLOCATE ( DP_HB(MNP) )
C--
      RETURN
      END SUBROUTINE

      END MODULE


