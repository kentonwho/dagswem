C----------------------------------------------------------------------------
C
C                           MODULE DECOMP
C
C----------------------------------------------------------------------------
C
C                  For use with ADCPREP Version 1.6 (  5/21/03 )
C
C                     current for ADCIRC v43.03   5/20/2003
C----------------------------------------------------------------------------
C
      SUBROUTINE DECOMP()
      USE pre_global 
      implicit none
C
C---------------------------------------------------------------------------C
C                     (  Serial Version 1.1  5/04/99  )                     C 
C                                                                           C
C  Decomposes the ADCIRC grid into NPROC subdomains.                        C
C  The Decomposition Variables are defined in the include file adcprep.inc  C 
C  This version is compatible with ADCIRC version 34.03                     C
C                                                                           C
C  12/14/98 vjp  Added interface to METIS 4.0                               C
C   3/10/99 vjp  Rewritten to allow Weir-node pairs to both be ghost nodes  C
C   4/05/99 vjp  Fixed bugs in metis interface routine                      C
C                                                                           C
C---------------------------------------------------------------------------C
C
      INTEGER N1, N2, N3, KMIN, VTMAX
      INTEGER I,J,JD,JG,JP,K,L,ITEMP,ITEMP2,IPR,IPR1,ICNT
      INTEGER ITOT,IEL,IELG,ILNODE,ILNODE2,IPROC,IPROC2
      INTEGER IG1,IG2,IG3,IL1,IL2,IL3,PE1,PE2,PE3
      INTEGER I1,DISC,BBN,IBP,NCOUNT
      INTEGER INDEX,INDEX2, KK
      INTEGER,ALLOCATABLE :: ITVECT(:)
      CHARACTER PE*6,GLOBAL*6

C     -----------------------------------
C     TACC Additions (k.schulz /9/13/08)
C     -----------------------------------

      integer elem_local
      integer ifound,ifound2,marked_node,node
      integer ierr,itemp_size
      integer, allocatable :: temp_loc_to_global(:,:)
      integer :: start, end, rate, max
      integer my_index,my_index2
      real(8) :: mytime
      real, dimension(0:2) :: position, velocity
      logical found
      integer, parameter :: MAX_MATCHES=5
      integer, parameter :: MAX_SUPPORTED_PROCS=2147483600
      real, parameter    :: LNODE_MULT=2.5
      type temp_node_info

         ! this data structure is simply used for
         ! searching during the creation of necessary
         ! data structures. Currently a 4-byte int is
         ! used which imposes a limit of MAX_SUPPORTED_PROCS
         ! above for the number of nweir and nprocs.

         integer (4) :: num_matches
         integer (4) :: matches(MAX_MATCHES)
      end type temp_node_info
      Type (temp_node_info), allocatable :: temp_nnodes(:)

      if(nproc .gt. MAX_SUPPORTED_PROCS)then
         print*,'** Error (decomp): NPROCS > MAX_SUPPORTED_PROCS'
         print*,'** Please use a larger integer data-type and recompile'
         stop
      endif

C     End TACC Additions -----------------------

C
      VTMAX = 24*MNP
      ALLOCATE ( ITVECT(VTMAX) )

C     TACC Allocation - Set aside some temp space 
C     that is Order(# of local nodes,NPROCS)

      itemp_size = LNODE_MULT*(MNP/NPROC)
      allocate ( temp_loc_to_global(itemp_size,NPROC),stat=ierr )
      if(ierr .ne. 0) then
         print*,'** Error (decomp): unable to allocate memory for temp_loc_to_global'
         stop
      endif

C   STEP 1:
C-- Use Partition of nodes to compute the number of Resident Nodes
C   to be assigned to each processor.
C-- Then construct Local-to-Global and Global-to-Local Node maps
C   for resident nodes: IMAP_NOD_LG(I,PE),  IMAP_NOD_GL(1:2,I) 

      DO I=1, NPROC             ! Use METIS 4.0  Partition
         NOD_RES_TOT(I) = 0
      ENDDO
      DO J=1, NNODG
         NCOUNT = NOD_RES_TOT(PROC(J))+1
         IMAP_NOD_GL(1,J) = PROC(J)
         IMAP_NOD_GL(2,J) = NCOUNT
         IMAP_NOD_LG(NCOUNT,PROC(J)) = J
         NOD_RES_TOT(PROC(J)) = NCOUNT
      ENDDO
C     DO I = 1, NNODG
C        print *, I, IMAP_NOD_GL(1,I)
C     ENDDO

C STEP 2: (Original Code)
C  Construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
C  Add an element to the map if it has an resident node
 
C       DO I = 1,NPROC
C          NELP(I) = 0
C          DO K = 1,NELG
C             N1 = NNEG(1,K)
C             N2 = NNEG(2,K)
C             N3 = NNEG(3,K)
C             PE1 = IMAP_NOD_GL(1,N1) ! Is any vertex a resident node?
C             PE2 = IMAP_NOD_GL(1,N2)
C             PE3 = IMAP_NOD_GL(1,N3)
C             IF ((PE1.EQ.I).OR.(PE2.EQ.I).OR.(PE3.EQ.I)) THEN
C               NELP(I) = NELP(I) + 1
C               IMAP_EL_LG(NELP(I),I) = K
C             ENDIF
C          ENDDO
C          IF (NELP(I).GT.MNEP) STOP 'NELP(I) > MNEP'
C       ENDDO


C STEP 2: Koomie Version
C
C  Construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
C  Add an element to the map if it has an resident node
C
C  k.schulz (9/12/08) - do it more efficiently using only a 
C                       single loop over the elements.

      print*,''
      print*,'Koomie (decomp): new local_to_global element map creation'
      call system_clock(count=start, count_rate=rate, count_max=max)

      do i = 1,nproc
         nelp(i) = 0
      enddo
 
      do k = 1,nelg

         n1 = nneg(1,k)
         n2 = nneg(2,k)
         n3 = nneg(3,k)
         pe1 = imap_nod_gl(1,n1) ! is any vertex a resident node?
         pe2 = imap_nod_gl(1,n2)
         pe3 = imap_nod_gl(1,n3)

         do i = 1,nproc
            if ((pe1.eq.i).or.(pe2.eq.i).or.(pe3.eq.i)) then
               nelp(i) = nelp(i) + 1
               imap_el_lg(nelp(i),i) = k
            endif
            if (nelp(i).gt.mnep) stop 'nelp(i) > mnep'
         enddo
      enddo

!      do i = 1,nproc
!         write(*,'1x,a,i5,a,i10') '  --> Proc: ',i,
!     $        ' # of elements = ',nelp(i)
!      enddo

      call system_clock(count=end, count_rate=rate, count_max=max)
      mytime=REAL(end - start)/REAL(rate)

      print*,'--> Initial loc_to_global element map complete'
      write(*,'(a,f6.2,a)') ' --> Time to complete = ',mytime,' (secs)'


C STEP 3:
C--Using Local-to-Global Element map
C  Construct Local-to-Global Node map:  IMAP_NOD_LG(I,PE)
C  and reconstruct Global-to-Local map for resident nodes
C

      print*,' '
      print*,'Orig (decomp): Original Step 3 code...'

      call system_clock(count=start, count_rate=rate, count_max=max)

      DO I = 1,NPROC
         ITOT = 0
         DO J = 1,NELP(I)
            IEL = IMAP_EL_LG(J,I)
            DO L=1, 3
               ITOT = ITOT + 1
               ITVECT(ITOT) = NNEG(L,IEL)
            ENDDO
         ENDDO
         ITEMP = ITOT     
         IF (ITOT.GT.VTMAX) stop 'step3 decomp'
         CALL SORT(ITEMP,ITVECT) ! Sort and remove multiple occurrences
         ITOT = 1
         IMAP_NOD_LG(1,I) = ITVECT(1)
         IF (IMAP_NOD_GL(1,ITVECT(1)).EQ.I) THEN
           IMAP_NOD_GL(2,ITVECT(1))=1
         ENDIF
         DO J = 2, ITEMP
           IF (ITVECT(J).NE.IMAP_NOD_LG(ITOT,I)) THEN
             ITOT = ITOT + 1
             IMAP_NOD_LG(ITOT,I) = ITVECT(J)
             IF (IMAP_NOD_GL(1,ITVECT(J)).EQ.I) THEN
               IMAP_NOD_GL(2,ITVECT(J))=ITOT
             ENDIF
           ENDIF
         ENDDO
         NNODP(I) = ITOT
         IF (NNODP(I).GT.MNPP) STOP 'NNODP > MNPP'
      ENDDO

      call system_clock(count=end, count_rate=rate, count_max=max)
      mytime=REAL(end - start)/REAL(rate)
      write(*,'(a,f6.2,a)') ' --> Time to complete = ',mytime,' (secs)'

c     print *, "Number of Nodes Assigned to PEs"
c     DO I=1, NPROC
c        print *, I-1, NNODP(I)
c        DO J=1, NNODP(I)
c           print *, J,IMAP_NOD_LG(J,I)
c        ENDDO
c     ENDDO

C STEP 4 (Original Code):
C--If there are any global Weir-node pairs, construct
C  Local-to-Global Weir Node maps: WEIRP_LG(:,PE), WEIRDP_LG(:,PE)
C  Rule: if a global Weir node is assigned ( either as a resident or ghost node )
C        then make it and its dual a local Weir-node pair

C       IF (NWEIR.GT.0) THEN

C         DO I=1, NPROC
C            ITOT = 0
C            DO J = 1,NNODP(I)
C               CALL SEARCH(WEIR,NWEIR,IMAP_NOD_LG(J,I),INDEX)
C               IF (INDEX.NE.0) THEN
C                 ITOT = ITOT+1
C                 ITVECT(ITOT) = INDEX
C               ENDIF
C               CALL SEARCH(WEIRD,NWEIR,IMAP_NOD_LG(J,I),INDEX2)
C               IF (INDEX2.NE.0) THEN
C                 ITOT = ITOT+1
C                 ITVECT(ITOT) = INDEX2
C               ENDIF
C            ENDDO
C            NWEIRP(I) = 0
C            ITEMP = ITOT
C            IF (ITOT.GT.VTMAX) stop 'step4 decomp'
C            IF (ITEMP.GT.1) THEN
C              CALL SORT(ITEMP,ITVECT)
C              ITOT=1
C              INDEX = ITVECT(1)
C              WEIRP_LG(ITOT,I)  = WEIR(INDEX)       
C              WEIRDP_LG(ITOT,I) = WEIRD(INDEX)
C              DO J = 2,ITEMP
C                 IF (ITVECT(J).NE.INDEX) THEN
C                   INDEX = ITVECT(J)
C                   ITOT = ITOT+1
C                   WEIRP_LG(ITOT,I)  = WEIR(INDEX)
C                   WEIRDP_LG(ITOT,I) = WEIRD(INDEX)
C                 ENDIF
C              ENDDO
C              NWEIRP(I) = ITOT
C            ENDIF
C c       DO J = 1, NWEIRP(I)
C c          print *, J, WEIRP_LG(J,I),WEIRDP_LG(J,I)
C c       ENDDO
C c       print *, "decomp: Number of WEIR node-pairs on PE",I-1,
C c    &           " = ",NWEIRP(I)
C         ENDDO
C       ELSE
C         DO I=1, NPROC
C            NWEIRP(I) = 0
C c          print *, "decomp: Number of WEIR node-pairs on PE",I-1,
C c    &              " = ",NWEIRP(I)
C         ENDDO

C       ENDIF 


C STEP 4 (Koomie Version):
C--If there are any global Weir-node pairs, construct
C  Local-to-Global Weir Node maps: WEIRP_LG(:,PE), WEIRDP_LG(:,PE)
C  Rule: if a global Weir node is assigned ( either as a resident or ghost node )
C        then make it and its dual a local Weir-node pair

C     ------------------------------------------------------------
C     Step4: TACC Rewrite; let's use a temporary array to flag
C     all Weir nodes to avoid oodles of calls to search()
C     
C     k.schulz 9/14/08
C     ------------------------------------------------------------

      if (nweir.gt.0) then

         print*,' '
         print*,'Koomie (decomp): Step 4 -> local-to-global maps with Weir nodes'
         call system_clock(count=start, count_rate=rate, count_max=max)

         if(nweir .gt. MAX_SUPPORTED_PROCS)then
            print*,'** Error (decomp): nweir > MAX_SUPPORTED_PROCS'
            print*,'** Please use a larger integer data-type and recompile'
            stop
         endif

         print*,'--> Initializing temporary search arrays'

         allocate(temp_nnodes(nnodg),stat=ierr )
         if(ierr .ne. 0) then
            print*,'** Error (decomp): unable to allocate memory for temp_nnodes'
            stop
         endif

C        Pre-processing, mark the indices of all Weir nodes to save
C        time on searches (note that the list of Weir nodes can
C        apparently have duplicates).  Only save the first index
C        for a given node (this matches the logic of the previous
C        Step 4 search code).

         print*,'--> Max global nodes     = ',nnodg
         print*,'--> Number of Weir nodes = ',nweir

         do i = 1,nnodg
            temp_nnodes(i)%matches(1) = 0
            temp_nnodes(i)%matches(2) = 0
         enddo
         
         do j = 1,nweir
            index  = weir(j)
            index2 = weird(j)
            if(temp_nnodes(index)%matches(1) .eq. 0)then
               temp_nnodes(index)%matches(1)  = j
            endif
            if(temp_nnodes(index2)%matches(2) .eq. 0)then
               temp_nnodes(index2)%matches(2) = j
            endif
         enddo

        DO I=1, NPROC
           ITOT = 0
           DO J = 1,NNODP(I)

              index  = temp_nnodes(imap_nod_lg(j,i))%matches(1)
              index2 = temp_nnodes(imap_nod_lg(j,i))%matches(2)

              if (index.ne.0) then
                 itot = itot+1
                 itvect(itot) = index
              endif

              if (index2.ne.0) then
                itot = itot+1
                itvect(itot) = index2
             endif
           enddo
           nweirp(i) = 0
           itemp = itot
           if (itot.gt.vtmax) stop 'step4 decomp'
           if (itemp.gt.1) then
             call sort(itemp,itvect)
             itot=1
             index = itvect(1)
             weirp_lg(itot,i)  = weir(index)       
             weirdp_lg(itot,i) = weird(index)
             do j = 2,itemp
                if (itvect(j).ne.index) then
                  index = itvect(j)
                  itot = itot+1
                  weirp_lg(itot,i)  = weir(index)
                  weirdp_lg(itot,i) = weird(index)
                endif
             enddo
             nweirp(i) = itot
           endif
        enddo

        call system_clock(count=end, count_rate=rate, count_max=max)
        mytime=REAL(end - start)/REAL(rate)
        write(*,'(a,f6.2,a)') ' --> Time to complete = ',mytime,' (secs)'

      else
        do i=1, nproc
           nweirp(i) = 0
        enddo

      endif 


C STEP 5 (Original Code):
C--If there are any global Weir-node pairs,
C  Re-construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
C  Rule:  Add an element if it has an resident node or 
C         has the dual Weir node of a resident or ghost node
 
C       IF (NWEIR.GT.0) THEN
C         DO I = 1,NPROC
C            NELP(I) = 0
C c          print *, "PE = ",I-1
C            DO K = 1,NELG
C               N1 = NNEG(1,K)
C               N2 = NNEG(2,K)
C               N3 = NNEG(3,K)
C               PE1 = IMAP_NOD_GL(1,N1)   ! Is any vertex a resident node?
C               PE2 = IMAP_NOD_GL(1,N2)
C               PE3 = IMAP_NOD_GL(1,N3)   ! belong to a Weir-node pair ?
C               CALL SEARCH3(WEIRP_LG(1,I),NWEIRP(I),N1,N2,N3,INDEX)
C               CALL SEARCH3(WEIRDP_LG(1,I),NWEIRP(I),N1,N2,N3,INDEX2)
C               IF ((PE1.EQ.I).OR.(PE2.EQ.I).OR.(PE3.EQ.I)
C      &             .OR.(INDEX.NE.0).OR.(INDEX2.NE.0)) THEN
C c                print *, K, PE1,PE2,PE3,INDEX,INDEX2
C                 NELP(I) = NELP(I) + 1
C                 KK = K
C                 if (EL_SHARE(NELP(I)) > -1) KK = -K
C                 IMAP_EL_LG(NELP(I),I) = KK
C                 EL_SHARE(NELP(I))     = I
C               ENDIF
C            ENDDO
C          IF (NELP(I).GT.MNEP) STOP 'NELP(I) > MNEP'
C c       print *, "Number of elements on PE",I-1," = ",NELP(I)
C c       DO J = 1, NELP(I)
C c          print *, J, IMAP_EL_LG(J,I)
C c       ENDDO
C         ENDDO
C       ENDIF

C STEP 5 (Koomie Version):
C--If there are any global Weir-node pairs,
C  Re-construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
C  Rule:  Add an element if it has an resident node or 
C         has the dual Weir node of a resident or ghost node

C     ------------------------------------------------------------
C     STEP5: TACC Rewrite; since NELG should be >> NPROC, 
C     let's do the looping a little more efficiently (ks 9/10/08).
C     ------------------------------------------------------------

        IF (NWEIR.GT.0) THEN

           print*,' '
           print*,'Koomie (decomp): new local_to_global element map',
     $          ' creation with Weir pairs'

           call system_clock(count=start, count_rate=rate, count_max=max)

           print*,'--> Initializing temporary search arrays'

           do i = 1,nproc
              nelp(i) = 0
           enddo

C          Pre-processing, mark all the Weir pairs to save searches (note that
C          temp_nnodes was allocated previously for Weir node search).

           do i = 1,nnodg
              temp_nnodes(i)%num_matches = 0
           enddo

           print*,'--> Max global nodes     = ',nnodg

           do i = 1,nproc
              
              do j = 1,nweirp(i)

                 index  = weirp_lg(j,i)
                 index2 = weirdp_lg(j,i)

                 ifound  = temp_nnodes(index)%num_matches  + 1
                 ifound2 = temp_nnodes(index2)%num_matches + 1

                 if( (ifound .gt. MAX_MATCHES).or. (ifound2 .gt. MAX_MATCHES) )then
                    print*,'** Error: Insufficient memory allocated in decomp()'
                    print*,'          --> Increase MAX_MATCHES; recompile, and rerun.'
                    stop
                 endif

                 temp_nnodes(index)%num_matches = ifound
                 temp_nnodes(index)%matches(ifound) = i

                 temp_nnodes(index2)%num_matches = ifound2
                 temp_nnodes(index2)%matches(ifound2) = i
              enddo
           enddo

           index = 0

           do i = 1,nnodg
              if (temp_nnodes(i)%num_matches  .gt. index)then
                 index = temp_nnodes(i)%num_matches
              endif
          enddo

          print*,'--> Max matches per node = ',index

C          To avoid multiple loops through the all elements, let's do this in two
C          steps.  Step 1 is to flag all the desired nodes and save the global 
C          numbering temporarily.

           print*,'--> Begin element search'
         
           DO K = 1,NELG
              N1 = NNEG(1,K)
              N2 = NNEG(2,K)
              N3 = NNEG(3,K)
              PE1 = IMAP_NOD_GL(1,N1) ! Is any vertex a resident node?
              PE2 = IMAP_NOD_GL(1,N2)
              PE3 = IMAP_NOD_GL(1,N3) ! belong to a Weir-node pair ?o

              marked_node = temp_nnodes(n1)%num_matches
     $                     +temp_nnodes(n2)%num_matches
     $                     +temp_nnodes(n3)%num_matches
            
              DO I = 1,NPROC

                 ! Cheap check #1

                 if ( (pe1.eq.i).or.(pe2.eq.i).or.(pe3.eq.i) ) then
                    nelp(i) = nelp(i) + 1
                    if(nelp(i) .gt. itemp_size ) then
                       print*,' '
                       print*,'** Error (decomp): itemp_size is too small, please increase.',itemp_size
                       print*,' '
                       stop
                    endif
                    
                    temp_loc_to_global(nelp(i),i) = k
                 elseif(marked_node .ne. 0)then

                    ! More expensive check #2

                    found = .false.
                    
                    do node = 1,temp_nnodes(n1)%num_matches
                       if(temp_nnodes(n1)%matches(node) == i) found=.true.
                    enddo

                    if( .not. found)then
                       do node = 1,temp_nnodes(n2)%num_matches
                          if(temp_nnodes(n2)%matches(node) == i) found=.true.
                       enddo
                    endif

                    if( .not. found)then
                       do node = 1,temp_nnodes(n3)%num_matches
                          if(temp_nnodes(n3)%matches(node) == i) found=.true.
                       enddo
                    endif

                    if(found)then
                       nelp(i) = nelp(i) + 1
                       temp_loc_to_global(nelp(i),i) = k
                    endif

                 endif

                 if (nelp(i).gt.mnep) stop 'nelp(i) > mnep'
              enddo
           enddo

!           do i = 1,nproc
!              write(*,'1x,a,i5,a,i10') 'Proc: ',i
!     $             ,' --> # of elements = ',nelp(i)
!           enddo

          ! Step 2: We've found all the desired elements, let's build the desired map.

           do i = 1,nproc
              do elem_local = 1,nelp(i)
                 kk = temp_loc_to_global(elem_local,i)
                 if (el_share(elem_local) > -1) then
                    kk = -kk
                 endif
               
                 imap_el_lg(elem_local,i) = kk
                 el_share(elem_local)     = i
              enddo
           enddo

           call system_clock(count=end, count_rate=rate, count_max=max)
           mytime=REAL(end - start)/REAL(rate)
           write(*,'(a,f6.2,a)') ' --> Time to complete = ',mytime,' (secs)'

           deallocate(temp_loc_to_global)
           deallocate(temp_nnodes)

        endif

C STEP 6:
C--Using Local-to-Global Element map
C  Construct Local-to-Global Node map:  IMAP_NOD_LG(I,PE)
C  and reconstruct Global-to-Local map for resident nodes
C
      DO I = 1,NPROC
         ITOT = 0
         DO J = 1,NELP(I)
            IEL = abs(IMAP_EL_LG(J,I))
            DO L=1, 3
               ITOT = ITOT + 1
               ITVECT(ITOT) = NNEG(L,IEL)
            ENDDO
         ENDDO
         ITEMP = ITOT     
         IF (ITOT.GT.VTMAX) stop 'step6 decomp'
         CALL SORT(ITEMP,ITVECT) ! Sort and remove multiple occurrences
         ITOT = 1
         IMAP_NOD_LG(1,I) = ITVECT(1)
         IF (IMAP_NOD_GL(1,ITVECT(1)).EQ.I) THEN
           IMAP_NOD_GL(2,ITVECT(1))=1
         ENDIF
         DO J = 2, ITEMP
           IF (ITVECT(J).NE.IMAP_NOD_LG(ITOT,I)) THEN
             ITOT = ITOT + 1
             IMAP_NOD_LG(ITOT,I) = ITVECT(J)
             IF (IMAP_NOD_GL(1,ITVECT(J)).EQ.I) THEN
               IMAP_NOD_GL(2,ITVECT(J))=ITOT
             ENDIF
           ENDIF
         ENDDO
         NNODP(I) = ITOT
         IF (NNODP(I).GT.MNPP) STOP 'NNODP > MNPP'
      ENDDO
c     print *, "Number of Nodes Assigned to PEs"
c     DO I=1, NPROC
c        print *, I-1, NNODP(I)
c        DO J=1, NNODP(I)
c           print *, J,IMAP_NOD_LG(J,I)
c        ENDDO
c     ENDDO

C STEP 7:
C--Construct Local Element Connectivity Table for each PE: NNEP(3,I,PE)
C
      DO I = 1,NPROC
         ITEMP = NNODP(I)
         DO J = 1,NNODP(I)
            ITVECT(J) = IMAP_NOD_LG(J,I)
         ENDDO
         DO J = 1,NELP(I)
            IELG = abs(IMAP_EL_LG(J,I))
            DO L = 1,3
               IG1 = NNEG(L,IELG)
               CALL LOCATE(ITVECT,ITEMP,IG1,K)
               IF (K.LE.0) THEN
                 IF (IMAP_NOD_LG(K+1,I).EQ.IG1) THEN
                   NNEP(L,J,I) = K+1
                 ELSE
                   STOP 'ERROR IN IMAP_NOD_LG'
                 ENDIF
               ELSEIF (K.GE.NNODP(I))THEN
                 IF (IMAP_NOD_LG(K,I).EQ.IG1) THEN
                   NNEP(L,J,I) = K
                 ELSE
                   STOP 'ERROR IN IMAP_NOD_LG'
                 ENDIF
               ELSE
                 IF (IMAP_NOD_LG(K,I).EQ.IG1) THEN
                   NNEP(L,J,I) = K
                 ELSEIF (IMAP_NOD_LG(K+1,I).EQ.IG1) THEN
                   NNEP(L,J,I) = K+1
                 ELSE 
                   STOP 'ERROR IN IMAP_NOD_LG'
                 ENDIF
               ENDIF
           ENDDO
         ENDDO
      ENDDO

C STEP 8:
C--Compute the number of communicating PEs and their 
C  list for each PE:  NUMM_COMM_PE(PE) and COMM_PE_NUM(IPE,PE)
C
      DO I = 1,NPROC
         NUM_COMM_PE(I) = 0
         ITEMP = 0
         DO J = 1,NNODP(I)
            INDEX = IMAP_NOD_LG(J,I)
            IPR = IMAP_NOD_GL(1,INDEX)
            IF (IPR.NE.I)THEN
              ITEMP = ITEMP + 1
              ITVECT(ITEMP) = IPR
            ENDIF
         ENDDO
         IF (ITEMP.EQ.0) THEN
           NUM_COMM_PE(I) = 0
         ELSE
           IF (ITEMP.GT.VTMAX) stop 'step8 decomp'
           CALL SORT(ITEMP,ITVECT)
           COMM_PE_NUM(1,I) = ITVECT(1)
           ITOT = 1
           DO J = 1,ITEMP
             IF (ITVECT(J).NE.COMM_PE_NUM(ITOT,I)) THEN
               ITOT = ITOT + 1
               COMM_PE_NUM(ITOT,I) = ITVECT(J)
             ENDIF
           ENDDO
           NUM_COMM_PE(I) = ITOT
           IF (NUM_COMM_PE(I).GT.MNPROC) STOP 'NUM_COMM_PE>MNPROC'
         ENDIF
      ENDDO

C STEP 9:
C--Construct a Global-to-Local node mapping: IMAP_NOD_GL2(*,J)
C  This is not a function, but is rather a relation
C  It works for both resident and ghost nodes

      DO I = 1,NNODG
         ITOTPROC(I) = 0
      ENDDO
      DO I = 1,NPROC
         DO J = 1,NNODP(I)
            INDEX = IMAP_NOD_LG(J,I)
            ITOTPROC(INDEX) = ITOTPROC(INDEX) + 1
            IF (ITOTPROC(INDEX).GT.MNPROC)THEN
              WRITE(6,*)'Some nodes belong to more processors',
     $                  ' than MNPROC'
              STOP 
            ENDIF
            ITEMP = (ITOTPROC(INDEX)-1)*2 + 1
            IMAP_NOD_GL2(ITEMP,INDEX) = I
            IMAP_NOD_GL2(ITEMP+1,INDEX) = J
         ENDDO
      ENDDO
c     print *, "Global Nodes assigned to more than one PE"
c     do J=1, NNODG
c        if (ITOTPROC(J).GT.1) print *, J, ITOTPROC(J)
c     enddo

C     -------------------------------------------------------------
C     Step 9a: Koomie: Set values for MNEP/MNEP
C     now that we have calulated final exact values.  This reduces
C     the over-estimate that was provided in domsize() which
C     were use for memory allocation
C     -------------------------------------------------------------

      print*,' '
      print*,'Koomie (decomp): Determine final sizing info'

      mnep=0
      mnpp=0

      do i = 1,nproc
         if ( nnodp(i) .gt. mnpp) then
            mnpp = nnodp(i)
         endif

         if ( nelp(i) .gt. mnep) then
            mnep = nelp(i)
         endif
      enddo

      write(*,'(a,i10)') '  --> MNEP = ',mnep
      write(*,'(a,i10)') '  --> MNPP = ',mnpp
      print*,' '
 
C STEP 10:
C--Print Summary of Decomposition
C
      print *, "Decomposition Data"
      print *, "DOMAIN  RES_NODES  GHOST_NODES  TOT_NODES  ELEMENTS"
      print *, "------  ---------  -----------  ---------  --------"
      GLOBAL = "GLOBAL"
      WRITE(*,90) GLOBAL,NNODG,NELG      
      DO I=1, NPROC
         PE(1:6) = 'PE0000'
         CALL IWRITE(PE,3,6,I-1)
         WRITE(6,92) PE, NOD_RES_TOT(I),NNODP(I)-NOD_RES_TOT(I),
     &                   NNODP(I),NELP(I)
      ENDDO
 90   FORMAT(1X,A6,25X,I9,2X,I9)
 92   FORMAT(1X,A6,1X,I9,2X,I9,4X,I9,2X,I9)
C
      RETURN
      END

      SUBROUTINE DOMSIZE()
      USE pre_global 
      implicit none

C     ---------------------------------------------------------------
C     The main purpose of this function appears to be to estimate 
C     the max number of nodes and elements assigned to any subdomain
C     to determine MNPP and MNEP.
C     
C     MNPP and MNEP will then be used to allocate appropriately sized 
C     arrays for the remaining data-structures to be created.
C     Consequently, we do not need exact numbers here and there is 
C     no point wasting lots of cycles to get an upper bounds.
C     This is an attempt at a more stripped-down version which 
C     computes the number of elements per processor and then
C     estimates the rest based on simple geometry arguments.
C     
C     k.schulz: (9/12/08)
C     ---------------------------------------------------------------

      integer i,j,k,ncount
      integer n1,n2,n3
      integer pe1,pe2,pe3
      integer mnpp_temp,mnep_temp
      integer ierr

      integer :: start, end, rate, max
      real(8) :: mytime

      integer,allocatable :: node_gl1(:)
      integer,allocatable :: node_lg(:)

      real :: PAD_FACTOR  = 0.05    ! Pad estimates by 5%
      real :: NODE_FACTOR = 1.4     ! Conversion factor from elems to nodes

      print*,''
      print*,'Koomie (domsize): estimating max number of elems/nodes per processor'

      call system_clock(count=start, count_rate=rate, count_max=max)

      allocate ( node_gl1(mnp),stat=ierr )
      if(ierr .ne. 0) then
         print*,'** Error (decomp): unable to allocate memory for node_gl1'
         stop
      endif 

C     Step 1: Compute # of elements per processor

      do i=1, nproc            
         nod_res_tot(i) = 0
      enddo

      do j=1, nnodg
         ncount = nod_res_tot(proc(j))+1
         node_gl1(j) = proc(j)
         nod_res_tot(proc(j)) = ncount
      enddo

      do i = 1,nproc
         nelp(i) = 0
      enddo
 
      do k = 1,nelg

         n1 = nneg(1,k)
         n2 = nneg(2,k)
         n3 = nneg(3,k)

         pe1 = node_gl1(n1) ! is any vertex a resident node?
         pe2 = node_gl1(n2)
         pe3 = node_gl1(n3)

         do i = 1,nproc
            if ((pe1.eq.i).or.(pe2.eq.i).or.(pe3.eq.i)) then
               nelp(i) = nelp(i) + 1
            endif
         enddo
      enddo

C     Find the proc with the largest number of local elements

      mnep_temp = 0

      do i = 1,nproc
         if(nelp(i) .gt. mnep_temp)then
            mnep_temp = nelp(i)
         endif
      enddo

C     Step 2: Estimate the rest; first, pad the number of 
C     elements (if there are Weir pairs).

      if(nweir .gt. 0)then
         mnep_temp = mnep_temp + (nweir/nproc)
      endif

C     Estimate # nodes for a 2d unstructured mesh

      mnpp_temp = mnep_temp/node_factor

C     Add a little additional padding for good measure.

      mnep = mnep_temp + PAD_FACTOR*mnep_temp
      mnpp = mnpp_temp + PAD_FACTOR*mnpp_temp

      call system_clock(count=end, count_rate=rate, count_max=max)
      mytime=REAL(end - start)/REAL(rate)

      print*,'--> Setting MNPP = ',MNPP
      print*,'--> Setting MNEP = ',MNEP
      write(*,'(a,f6.2,a)') ' --> Time to complete = ',mytime,' (secs)'
      print*,' '

      deallocate ( node_gl1 )

      return
      end

      SUBROUTINE SORT(N,RA)
      IMPLICIT NONE
      INTEGER N, L, IR, RRA, I, J
      INTEGER RA(N)
C 
C---------------------------------------------------------------------------        
C  Sorts array RA of length N into ascending order using Heapsort algorithm.
C  N is input; RA is replaced on its output by its sorted rearrangement.
C  Ref: Numerical Recipes
C--------------------------------------------------------------------------- 
C
      L = N/2 + 1
      IR = N
10    CONTINUE
      IF (L.GT.1)THEN
        L=L-1
        RRA = RA(L)
      ELSE
        RRA=RA(IR)
        RA(IR)=RA(1)
        IR=IR-1
        IF (IR.EQ.1) THEN
          RA(1)=RRA
          RETURN
        ENDIF
      ENDIF
      I=L
      J=L+L
20    IF (J.LE.IR) THEN
        IF (J.LT.IR) THEN
          IF(RA(J).LT.RA(J+1)) J=J+1
        ENDIF
        IF (RRA.LT.RA(J)) THEN
          RA(I)=RA(J)
          I=J
          J=J+J
        ELSE
          J=IR+1
        ENDIF
        GO TO 20
      ENDIF
      RA(I)=RRA
      GO TO 10
      END

      SUBROUTINE LOCATE(XX,N,X,J)
      IMPLICIT NONE
      INTEGER JM,JL,JU,J,N,X,XX(N)
C
C--Given an array XX of length N, and given a value X, returns a value J
C--such that X is between XX(J) and XX(J+1). XX must be monotonic, either
C--increasing or decreasing. J=0 or J=N is returned to indicate that X is out 
C--of range.
C--
C--NUMERICAL RECIPES - The Art of Scientific Computing [FORTRAN Version]

C--Initialize lower and upper limits
         JL = 0
         JU = N+1
C
C--If we are not done yet, compute a mid-point, and replace either the lower 
C--limit or the upper limit, as appropriate.
C
10       IF(JU-JL.GT.1)THEN
            JM = (JU+JL)/2
            IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
              JL = JM
            ELSE
              JU = JM
            ENDIF
C--Repeat until the test condition 10 is satisfied.
         GO TO 10
         ENDIF
C--Then set the output and return.
         J = JL
         RETURN
         END

      SUBROUTINE SEARCH3(MAP,LEN,N1,N2,N3,INDEX)
      USE pre_global
      INTEGER MAP(*),LEN,N1,N2,N3,INDEX,IP
cvjp  rewritten 5/3/99
      INDEX = 0
      DO I=1,LEN
         IP = MAP(I)
         IF (IP.EQ.N1.OR.IP.EQ.N2.OR.IP.EQ.N3) THEN
           INDEX = I
           GOTO 99
         ENDIF
      ENDDO
 99   RETURN
      END

      SUBROUTINE SEARCH(MAP,N,TARGET,INDEX)
      INTEGER MAP(*),N,TARGET,INDEX,I
C
      INDEX = 0
      IF (N.EQ.0) GOTO 99
      DO I=1,N
         IF (MAP(I).EQ.TARGET) THEN
            INDEX = I
            GOTO 99
         ENDIF
      ENDDO
99    RETURN
      END
