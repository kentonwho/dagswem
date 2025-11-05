C***********************************************************************
C
C     SUBROUTINE PREP_SLOPELIM()
C
C     This subroutine does preparatory stuff for the slope limiter
C
C     Written by Ethan Kubatko (07-13-2005)
C
C***********************************************************************

      SUBROUTINE PREP_SLOPELIM()
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER GED
      
C.....Initialize SL3
      
      SL3 = 0

C.....Loop over the elements

      DO 111 J = 1,MNE
        
C.....Retrieve the nodal coordinates of the given element

        N1 = NM(J,1)
        N2 = NM(J,2)
        N3 = NM(J,3)

C.....Compute the barycenter coordinates of the element and store

        XBC(J) = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
        YBC(J) = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))

        X1 = XBC(J)
        Y1 = YBC(J)
          
C.....Loop over the edges to fnd the neighboring elements

        DO I = 1,3
          
C.....Retrieve the global edge number of the element for the given edge

          GED = NELED(I,J)

C.....Retrieve the neighboring element number and store into an array

          EL_NBORS(I,J) = NEDEL(1,GED)

          IF (EL_NBORS(I,J).EQ.J) EL_NBORS(I,J) = NEDEL(2,GED)

C.....If the element has an edge that is on boundary go to next element
Csb-2007/07/27 commented out

C          IF (EL_NBORS(I,J).EQ.0) GOTO 111
            
        ENDDO
        
C.....Set the 4th neighbroing element equal to the first

        EL_NBORS(4,J) = EL_NBORS(1,J)

C.....Now loop over the three edges of the element again

        DO I = 1,3

          IF(EL_NBORS(I,J).EQ.0.OR.EL_NBORS(I+1,J).EQ.0) GOTO 110
          
C.....Compute the barycenter coordinates of two neighboring elements

          N1 = NM(EL_NBORS(I,J),1)
          N2 = NM(EL_NBORS(I,J),2)
          N3 = NM(EL_NBORS(I,J),3)
            
          X2 = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
          Y2 = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))
            
          N1 = NM(EL_NBORS(I+1,J),1)
          N2 = NM(EL_NBORS(I+1,J),2)
          N3 = NM(EL_NBORS(I+1,J),3)

          X3 = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
          Y3 = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))
            
C.....Compute the time independent planar constant

          SL3(I,J) = X1*(Y2 - Y3) + X2*(Y3 - Y1) + X3*(Y1 - Y2)

          IF (SL3(I,J).LE.0.AND.SLOPEFLAG.NE.0) then
            WRITE(16,*) 'WARNING. SL3(',I,',',J,') =',SL3(I,J),' <= 0.',
     &           '    ELEMENT ',J,
     &           ' WILL NOT BE CONSIDERED IN SLOPE LIMITING.'
            WRITE(16,*)
            WRITE(*,*) 'WARNING. SL3(',I,',',J,') =',SL3(I,J),' <= 0.',
     &           '    ELEMENT ',J,
     &           ' WILL NOT BE CONSIDERED IN SLOPE LIMITING.'
            WRITE(*,*)
            SL3(I,J) = 0.D0
          ENDIF
            
 110      CONTINUE

        ENDDO
        
111   CONTINUE
            
      RETURN
      END SUBROUTINE PREP_SLOPELIM
