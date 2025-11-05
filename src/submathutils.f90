MODULE DEF_BITREE_NODE
  IMPLICIT NONE

  TYPE, public:: BiTreeNode
     INTEGER:: Val
     TYPE (BiTreeNode), pointer:: right
     TYPE (BiTreeNode), pointer:: left
  END TYPE BiTreeNode   

 CONTAINS

  FUNCTION CreateBiTreeNode( Number ) RESULT(Tail)
    IMPLICIT NONE
    
    INTEGER:: Number
    TYPE (BiTreeNode), pointer:: Tail

    ALLOCATE(Tail) ; 
    Tail%Val = Number ; 
    
    nullify(Tail%right) ;
    nullify(Tail%left)  ;
    
    RETURN ;
  END FUNCTION CreateBiTreeNode    
  
  SUBROUTINE AddNodeBiTree( Dummy_Root, number )
    IMPLICIT NONE
    
    TYPE (BiTreeNode), pointer:: Dummy_Root
    INTEGER:: number
    
    TYPE (BiTreeNode), pointer:: Current, node
    
    Current => Dummy_Root ;
    DO
       IF ( number > Current%val ) THEN
          IF ( .NOT. associated(current%right) ) THEN
             Current%right => CreateBiTreeNode( number ) ;
             
             nullify(node) ;
             return ;   
          END IF
          Current => Current%right ;
       ELSE IF ( number < Current%val ) THEN
          IF ( .NOT. associated(current%left) ) then
             Current%left => CreateBiTreeNode( number ) ;
             
             nullify(node) ;
             return ;
          end if
          Current => Current%left ; 
       ELSE
          return ;
       END IF
    END DO
    
    RETURN ;
  END SUBROUTINE AddNodeBiTree
  
  RECURSIVE SUBROUTINE rFetchBiTreeSortedNodes( Current, n, ivec )
    IMPLICIT NONE
    
    !c Dummy variables c!
    INTEGER:: n
    INTEGER, pointer:: ivec(:)
    TYPE (BiTreeNode), pointer:: Current
    
    !c Local variables c!
    INTEGER:: l, m, npad
    INTEGER, pointer:: ivecdummy(:)
    
    IF ( associated(Current%Left) ) THEN
       CALL rFetchBiTreeSortedNodes(Current%Left, n, ivec ) ;
    END IF
    n = n + 1 ;
    
    m = ubound(ivec,1) ;
    npad =  m/2 ; 
    IF ( n > m ) THEN
       allocate(ivecdummy(m + npad)) ;
       
       ivecdummy(1:m) = ivec(1:m) ;
       ivecdummy(m+1:) = 0 ;
       
       deallocate(ivec) ;
       nullify(ivec) ;
       
       ivec => ivecdummy ; 
       nullify(ivecdummy) ;
    END IF
    ivec(n) = Current%val ;
    
    IF ( associated(Current%Right) ) THEN
       CALL rFetchBiTreeSortedNodes(Current%Right, n, ivec ) ;
    END IF
    
    RETURN ;
  END SUBROUTINE rFetchBiTreeSortedNodes

  RECURSIVE SUBROUTINE RemoveBiTreeLeaves( Current )
    IMPLICIT NONE
    
    !c Dummy variables c!
    TYPE (BiTreeNode), pointer:: Current
    
    !c Local variables c!
    TYPE (BiTreeNode), pointer:: DummyPT 
    
    IF ( .NOT. associated(Current) ) THEN
       return ;
    END IF
    
    !c Remove all leaves c!
    IF ( IsBiTreeLeaf(Current) ) THEN
       deallocate(Current) ;
       nullify(Current) ;
       return ;
    ELSE
       IF ( associated(Current%Left) ) THEN
          DummyPT => Current%Left ;
          IF ( IsBiTreeLeaf(DummyPT) ) THEN
             deallocate(Current%Left) ;
             nullify(Current%Left)
          ELSE    
             CALL RemoveBiTreeLeaves( DummyPT ) ; 
          END IF
       END IF
       
       IF ( associated(Current%Right) ) THEN
          DummyPT => Current%Right ;
          IF ( IsBiTreeLeaf(DummyPT) ) THEN
             deallocate(Current%Right) ;
             nullify(Current%Right) ;  
          ELSE    
             CALL RemoveBiTreeLeaves( DummyPT ) ; 
          END IF
       END IF
       
       deallocate(Current) ;
       nullify(Current) ;
    END IF
    
    RETURN ;
  END SUBROUTINE RemoveBiTreeLeaves
  
  SUBROUTINE DeleteBiTree0( Root )
    IMPLICIT NONE
    
    !c Dummy variables c!
    TYPE (BiTreeNode), pointer:: Root
    
    !c Local variables c!
    TYPE (BiTreeNode), pointer:: Head
    INTEGER:: i 
    
    i = 0 ; 
    DO
       i = i + 1 ;
       IF ( .NOT. associated(Root) ) THEN
          EXIT ;
       ELSE
          Head => Root ;
          CALL RemoveBiTreeLeaves( Head ) ;
       END IF
       root => head ;
       
    END DO
    NULLIFY( Root ) ;
    
    RETURN ;
  END SUBROUTINE DeleteBiTree0
  
  FUNCTION IsBiTreeLeaf( Node )
    IMPLICIT NONE
    
    TYPE (BiTreeNode), pointer:: Node
    LOGICAL:: IsBiTreeLeaf
    
    IsBiTreeLeaf = .false. ;
    IF ( .not. associated(Node%Left) .AND. &
         & .not. associated(Node%Right) ) THEN
       IsBiTreeLeaf = .true. 
    END IF
    
    RETURN ;
  END FUNCTION IsBiTreeLeaf
  
END MODULE DEF_BITREE_NODE

MODULE ModBiTreeSort
  USE DEF_BITREE_NODE
  
  IMPLICIT NONE
  
  !c For a binary-tree search c!
  TYPE (BiTreeNode), private, pointer:: BiTreeRoot
  
CONTAINS  
  
  SUBROUTINE BuildBiTree( ivec )
    IMPLICIT NONE
    
    !c Dummy variables c!
    INTEGER:: ivec(:)
    
    INTEGER:: i, m
    LOGICAL, save:: first = .true. ;
    
    IF ( first ) THEN
       BiTreeRoot => CreateBiTreeNode( ivec(1) ) ;
       first = .false. ;
    ELSE
       IF ( .not. associated(BiTreeRoot) ) THEN 
          BiTreeRoot => CreateBiTreeNode( ivec(1) ) ;
       END IF
    END IF
    
    m = ubound(ivec, 1) ;
    DO i = 1, m
       CALL AddNodeBiTree( BiTreeRoot, ivec(i) ) ;
    END DO
    
    RETURN ;
  END SUBROUTINE BuildBiTree
  
  FUNCTION FetchBiTreeNodes( nnode ) RESULT( ivec )
    IMPLICIT NONE
    
    INTEGER:: nnode
    INTEGER, pointer:: ivec(:)
    
    INTEGER:: n, ip
    
    n = 40 ; 
    allocate( ivec(n) ) ;
    
    nnode = 0 ;
    CALL rFetchBiTreeSortedNodes( BiTreeRoot, nnode, ivec ) ;
    
    RETURN ;
  END FUNCTION FetchBiTreeNodes
  
  SUBROUTINE DestroyBiTree( )
    IMPLICIT NONE
    
    CALL DeleteBiTree0( BiTreeRoot ) ;
    
    RETURN ;
  END SUBROUTINE DestroyBiTree
  
!c========================= Binary tree sort =====================c!
  !c Work for a vector with equal elements
  !c The efficiency depends strong on the randomness of the data
  !c  - best case log_{2}(ubound(ivec,1))
  !c  - worst case n
  FUNCTION BiTreeSortIntVec( nval, ivec ) RESULT( sivec )
    IMPLICIT NONE
    
    INTEGER:: nval, ivec(:)
    INTEGER, pointer:: sivec(:) 
    
    CALL BuildBiTree( ivec ) ;

    sivec => FetchBiTreeNodes( nval ) ;

    CALL DestroyBiTree() ;
    
    RETURN ;
  END FUNCTION BiTreeSortIntVec
  
END MODULE ModBiTreeSort


MODULE DEF_BITREE_NODE1
  IMPLICIT NONE

  TYPE, public:: BiTreeNode1
     INTEGER:: idx
     INTEGER:: Val(2)
     TYPE (BiTreeNode1), pointer:: right
     TYPE (BiTreeNode1), pointer:: left
  END TYPE BiTreeNode1
  
 CONTAINS
  
  FUNCTION CreateBiTreeNode1( Number, idx ) RESULT(Tail)
    IMPLICIT NONE
    
    INTEGER:: Number(2), idx
    TYPE (BiTreeNode1), pointer:: Tail
    
    ALLOCATE(Tail) ; 
    Tail%Val = Number ; 
    Tail%idx = idx ;
    
    nullify(Tail%right) ;
    nullify(Tail%left)  ;
    
    RETURN ;
  END FUNCTION CreateBiTreeNode1
  
  SUBROUTINE FindNodeBiTree1( Dummy_Root, key, Found, idx )
    IMPLICIT NONE

    !c Dummy c!
    INTEGER:: idx
    LOGICAL:: Found
    
    INTEGER:: key(2)
    TYPE (BiTreeNode1), pointer:: Dummy_Root
    
    !c Local c!
    TYPE (BiTreeNode1), pointer:: Current, node
    
    idx = 0 ; 
    Found = .FALSE. ;
    Current => Dummy_Root ;
    DO
       IF ( key(1) > Current%val(1) ) THEN
          IF ( associated(current%right) ) THEN
             Current => Current%right ;
          ELSE
             Found = .FALSE. ;
             return ;
          END IF
       ELSE IF ( key(1) < Current%val(1) ) THEN
          IF ( associated(current%left) ) then
             Current => Current%left ; 
          ELSE
             Found = .FALSE. ;
             return ;
          END IF   
       ELSE
          IF ( key(2) > Current%val(2) ) THEN
             IF ( associated(current%right) ) THEN
                Current => Current%right ;
             ELSE
                Found = .FALSE. ;
                return ;
             END IF
          ELSE IF ( key(2) < Current%val(2) ) THEN
             IF ( associated(current%left) ) then
                Current => Current%left ; 
             ELSE
                Found = .FALSE. ;
                return ;
             END IF   
          ELSE
             idx = Current%idx ;
             Found = .TRUE. ;
             return ;
          END IF
       END IF
    END DO
 
  END SUBROUTINE FindNodeBiTree1  

  SUBROUTINE AddNodeBiTree1( Dummy_Root, number, idx )
    IMPLICIT NONE

    TYPE (BiTreeNode1), pointer:: Dummy_Root
    INTEGER:: number(2), idx

    TYPE (BiTreeNode1), pointer:: Current, node

    Current => Dummy_Root ;
    DO
       IF ( number(1) > Current%val(1) ) THEN
          IF ( .NOT. associated(current%right) ) THEN
             Current%right => CreateBiTreeNode1( number, idx ) ;
             
             nullify(node) ;
             return ;   
          END IF
          Current => Current%right ;
       ELSE IF ( number(1) < Current%val(1) ) THEN
          IF ( .NOT. associated(current%left) ) then
             Current%left => CreateBiTreeNode1( number, idx ) ;

             nullify(node) ;
             return ;
          END IF
          Current => Current%left ; 
       ELSE
          IF ( number(2) > Current%val(2) ) THEN
             IF ( .NOT. associated(current%right) ) THEN
                Current%right => CreateBiTreeNode1( number, idx ) ;
                
                nullify(node) ;
                return ;   
             END IF
             Current => Current%right ;
          ELSE IF ( number(2) < Current%val(2) ) THEN
             IF ( .NOT. associated(current%left) ) then
                Current%left => CreateBiTreeNode1( number, idx ) ;
                
                nullify(node) ;
                return ;
             END IF
             Current => Current%left ; 
          ELSE
             return ;
          END IF
       END IF
    END DO
    
    RETURN ;
  END SUBROUTINE AddNodeBiTree1
  
  RECURSIVE SUBROUTINE rFetchBiTreeSortedNodes1( Current, n, ivec )
    IMPLICIT NONE

    !c Dummy variables c!
    INTEGER:: n
    INTEGER, pointer:: ivec(:,:)
    TYPE (BiTreeNode1), pointer:: Current

    !c Local variables c!
    INTEGER:: l, m, npad
    INTEGER, pointer:: ivecdummy(:,:)
    
    IF ( associated(Current%Left) ) THEN
       CALL rFetchBiTreeSortedNodes1(Current%Left, n, ivec ) ;
    END IF
    n = n + 1 ;

    m = ubound(ivec,1) ;
    npad =  m/2 ; 
    IF ( n > m ) THEN
!!$       allocate(ivecdummy(m + npad,2)) ;
       allocate(ivecdummy(m + npad, 3)) ;
       
       ivecdummy(1:m,1:3) = ivec(1:m,1:3) ;
       ivecdummy(m+1:,:) = 0 ;
       
       deallocate(ivec) ;
       nullify(ivec) ;

       ivec => ivecdummy ; 
       nullify(ivecdummy) ;
    END IF
    ivec(n,1:2) = Current%val ;
    ivec(n,3) = Current%idx ;

    IF ( associated(Current%Right) ) THEN
       CALL rFetchBiTreeSortedNodes1(Current%Right, n, ivec ) ;
    END IF
    
    RETURN ;
  END SUBROUTINE rFetchBiTreeSortedNodes1

  RECURSIVE SUBROUTINE RemoveBiTreeLeaves1( Current )
    IMPLICIT NONE
    
    !c Dummy variables c!
    TYPE (BiTreeNode1), pointer:: Current
    
    !c Local variables c!
    TYPE (BiTreeNode1), pointer:: DummyPT 
    
    IF ( .NOT. associated(Current) ) THEN
       return ;
    END IF
    
    !c Remove all leaves c!
    IF ( IsBiTreeLeaf1(Current) ) THEN
       deallocate(Current) ;
       nullify(Current) ;
       return ;
    ELSE
       IF ( associated(Current%Left) ) THEN
          DummyPT => Current%Left ;
          IF ( IsBiTreeLeaf1(DummyPT) ) THEN
             deallocate(Current%Left) ;
             nullify(Current%Left)
          ELSE    
             CALL RemoveBiTreeLeaves1( DummyPT ) ; 
          END IF
       END IF
       
       IF ( associated(Current%Right) ) THEN
          DummyPT => Current%Right ;
          IF ( IsBiTreeLeaf1(DummyPT) ) THEN
             deallocate(Current%Right) ;
             nullify(Current%Right) ;  
          ELSE    
             CALL RemoveBiTreeLeaves1( DummyPT ) ; 
          END IF
       END IF
       
       deallocate(Current) ;
       nullify(Current) ;
    END IF
    
    RETURN ;
  END SUBROUTINE RemoveBiTreeLeaves1  

  SUBROUTINE DeleteBiTree0_1( Root )
    IMPLICIT NONE
    
    !c Dummy variables c!
    TYPE (BiTreeNode1), pointer:: Root
    
    !c Local variables c!
    TYPE (BiTreeNode1), pointer:: Head
    INTEGER:: i 
    
    i = 0 ; 
    DO
       i = i + 1 ;
       IF ( .NOT. associated(Root) ) THEN
          EXIT ;
       ELSE
          Head => Root ;
          CALL RemoveBiTreeLeaves1( Head ) ;
       END IF
       root => head ;
    END DO
    NULLIFY( Root ) ;
    
    RETURN ;
  END SUBROUTINE DeleteBiTree0_1
  
  FUNCTION IsBiTreeLeaf1( Node )
    IMPLICIT NONE
    
    TYPE (BiTreeNode1), pointer:: Node
    LOGICAL:: IsBiTreeLeaf1

    IsBiTreeLeaf1 = .false. ;
    IF ( .not. associated(Node%Left) .AND. &
         & .not. associated(Node%Right) ) THEN
       IsBiTreeLeaf1 = .true. 
    END IF
    
    RETURN ;
  END FUNCTION IsBiTreeLeaf1
  
END MODULE DEF_BITREE_NODE1

MODULE ModBiTreeSort1
  USE DEF_BITREE_NODE1
  
  IMPLICIT NONE
  
  TYPE (BiTreeNode1), private, pointer:: BiTreeRoot2N

 CONTAINS  

  FUNCTION SetPtrToBiTreeRoot2N( ) RESULT( Root )
    IMPLICIT NONE

    TYPE (BiTreeNode1), pointer:: Root

    Root => BiTreeRoot2N ;

    RETURN ;
  END FUNCTION SetPtrToBiTreeRoot2N

  SUBROUTINE SetBiTreeRoot2NToPtr( Root )
    IMPLICIT NONE
    
    TYPE (BiTreeNode1), pointer:: Root
    
    BiTreeRoot2N => Root ;
    
    RETURN ;
  END SUBROUTINE SetBiTreeRoot2NToPtr

  !c Alway build a tree c!
  SUBROUTINE BuildBiTree1_V1( ivec, idx )
    IMPLICIT NONE
    
    !c Dummy variables c!
    INTEGER:: ivec(:,:), idx(:)
    
    INTEGER:: i, m
    LOGICAL, save:: first = .true. ;
    
    IF ( first ) THEN
       BiTreeRoot2N => CreateBiTreeNode1( ivec(1,:), idx(1) ) ;
       first = .false. ;
    ELSE
       IF ( .not. associated(BiTreeRoot2N) ) THEN 
          BiTreeRoot2N => CreateBiTreeNode1( ivec(1,:), idx(1) ) ;
       ELSE
          nullify(BiTreeRoot2N) ;
          BiTreeRoot2N => CreateBiTreeNode1( ivec(1,:), idx(1) ) ;
       END IF   
    END IF
    
    m = ubound(ivec, 1) ;
    DO i = 1, m
       CALL AddNodeBiTree1( BiTreeRoot2N, ivec(i,:), idx(i) ) ;
    END DO
    
    RETURN ;
  END SUBROUTINE BuildBiTree1_V1  

  !c Add new data to the existing data c!
  SUBROUTINE BuildBiTree1_V2( ivec, idx )
    IMPLICIT NONE
    
    !c Dummy variables c!
    INTEGER:: ivec(:,:), idx(:)
    
    INTEGER:: i, m
    LOGICAL, save:: first = .true. ;
    
    IF ( first ) THEN
       BiTreeRoot2N => CreateBiTreeNode1( ivec(1,:), idx(1) ) ;
       first = .false. ;
    ELSE
       IF ( .not. associated(BiTreeRoot2N) ) THEN 
          BiTreeRoot2N => CreateBiTreeNode1( ivec(1,:), idx(1) ) ;
       END IF
    END IF
    
    m = ubound(ivec, 1) ;
    DO i = 1, m
       CALL AddNodeBiTree1( BiTreeRoot2N, ivec(i,:), idx(i) ) ;
    END DO
    
    RETURN ;
  END SUBROUTINE BuildBiTree1_V2 

  FUNCTION FetchBiTreeNodes1( nnode ) RESULT( ivec )
    IMPLICIT NONE
    
    INTEGER:: nnode
    INTEGER, pointer:: ivec(:,:)

    INTEGER:: n, ip
    
    n = 20 ; 
!!$    allocate( ivec(n,2) ) ;
    allocate( ivec(n,3) ) ;

    nnode = 0 ;
    CALL rFetchBiTreeSortedNodes1( BiTreeRoot2N, nnode, ivec ) ;

    RETURN ;
  END FUNCTION FetchBiTreeNodes1

  SUBROUTINE SearchBiTree1( key, Found, idx )
    IMPLICIT NONE

    LOGICAL:: Found
    INTEGER:: idx, Key(2)

    idx = 0 ;
    Found = .FALSE. ;
    CALL FindNodeBiTree1( BiTreeRoot2N, key, Found, idx ) ;

    RETURN ;
  END SUBROUTINE SearchBiTree1    

  SUBROUTINE DestroyBiTree1( )
    IMPLICIT NONE
    
    CALL DeleteBiTree0_1( BiTreeRoot2N ) ;

    RETURN ;
  END SUBROUTINE DestroyBiTree1
  
  FUNCTION BiTreeSortIntVec1( nval, ivec, idx ) RESULT( sivec )
    IMPLICIT NONE
    
    INTEGER:: nval, ivec(:,:), idx(:)
    INTEGER, pointer:: sivec(:,:) 

    CALL BuildBiTree1_V1( ivec, idx ) 
    sivec => FetchBiTreeNodes1( nval ) ;
    CALL DestroyBiTree1() ;

    RETURN ;
  END FUNCTION BiTreeSortIntVec1
  
END MODULE ModBiTreeSort1

MODULE AltMathUtils
  USE ModBiTreeSort
  USE ModBiTreeSort1

  IMPLICIT NONE
  
 CONTAINS

!c===================== Binary Search ======================c!
  !c Perform a binary search for a 'Key' in array A(1:N) c! 
  SUBROUTINE IBinarySearch( A, N, Key, Position, Found )
    IMPLICIT NONE
    
    INTEGER:: A(:)
    INTEGER:: Key, N, Position
    LOGICAL:: Found
    
    INTEGER:: L, R, Mid
    
    Found = .false. ;
    L = lbound(A, 1) ; 
    IF ( N < L ) THEN 
       return ;
    END IF
    
    R = N ; 
    DO
       Mid = (L + R)/2 ; 
       
       IF ( Key < A(Mid) ) THEN
          IF ( L >= R ) THEN
             Position = L ;
             EXIT ;
          END IF
          R = Mid ; 
       ELSE IF ( Key > A(Mid) ) THEN
          IF ( L >= R ) THEN
             Position = L ; 
             EXIT ;
          END IF
          L = Mid + 1 ; 
       ELSE
          Position = Mid ;
          Found = .true. ;
          EXIT ;
       END IF
    END DO
    
    RETURN ;
  END SUBROUTINE IBinarySearch
    
!c========================= Hash sort =====================c!
  !c NOTE:
  !c   - May not work for an array having identical entries
  !c Source: R.A. Vowels, Algoithms and Data structures in F and Fortran
  SUBROUTINE IHashSort( A, N )
    IMPLICIT NONE
    
    !c Dummy variables c!
    INTEGER:: N
    INTEGER, dimension(-N:,1:):: A
    
    !c Local variables c!
    INTEGER:: Temp(2)
    INTEGER:: H, J, K
    REAL(8), parameter:: Factor = 2.D0
    
    REAL(8):: Divisor
    INTEGER:: Largest, Smallest, idxsm, idxval

    CALL iminval1( Smallest, idxsm, A(1:N,1), N ) ; 
    idxval = A(idxsm,2) ;
    
    Largest = maxval(A(1:N,1)) ; 
    
    IF ( Largest == smallest ) THEN
       return ;
    END IF
    Divisor = dble(Largest) - dble(smallest) ; 

    A(-N:0,1) = smallest ;
    DO J = 1, N
       IF ( A(J,1) /= Smallest ) THEN
          H = Factor*N*(dble(A(J,1)) - Smallest)/Divisor ;
          H = H - N ; 
          
          IF ( H < 0 ) THEN
             DO 
                IF ( A(H,1) == Smallest ) THEN
                   EXIT ; 
                END IF
                H = H + 1 ; 
             END DO
             A(H,1:2) = A(J,1:2) ;
             A(J,1:2) = (/ Smallest, idxval /) ; 
          END IF
       END IF
    END DO

    K = -N ;
    DO J = -N, N
       IF ( A(J,1) /= Smallest ) THEN
          IF ( J /= K ) THEN
             A(K,1:2) = A(J,1:2) ;
             A(J,1:2) = (/ Smallest, idxval /) ; 
          END IF
          K = K + 1 ; 
       END IF
    END DO

    DO J = -N, -1
       IF ( A(J,1) /= Smallest ) THEN 
          H = Factor*N*(real(A(J,1),8)-Smallest)/Divisor ;
          H = H - N ;
          IF ( H >= 0 ) THEN
             DO 
                IF ( A(H,1) == Smallest ) THEN
                   EXIT ; 
                END IF
                H = H - 1 ;
             END DO
             A(H,1:2) = A(J,1:2) ;
             A(J,1:2) = (/ Smallest, idxval /) ; 
          END IF
       END IF
    END DO

    K = N ; 
    DO J = N, -N, -1
       IF ( A(J,1) /= Smallest ) THEN
          IF ( J /= K ) THEN
             A(K,1:2) = A(J,1:2) ; 
             A(J,1:2) = (/ Smallest, idxval /) ; 
          END IF
          K = K - 1 ; 
       END IF
    END DO

    DO J = 1, N - 1
       IF ( A(J,1) > A(J+1,1) ) THEN
          Temp = A(J+1,1:2) ;
          DO K = J, 1, -1
             IF ( A(K,1) < Temp(1) ) THEN
                EXIT ; 
             END IF
             A(K+1,1:2) = A(K,1:2) ; 
          END DO
          A(K+1,1:2) = Temp ; 
       END IF
    END DO

    RETURN ;
  END SUBROUTINE IHashSort


  !c Find the minimum value in A(1:N)
  SUBROUTINE IMinVal1( val, idx, A, N )
    IMPLICIT NONE
    
    INTEGER:: N
    INTEGER:: val, idx, A(:)
    
    INTEGER:: L, R, Mid, mval(2)
    
    IF ( N == 1 ) THEN
       val = A(1) ;
       return ; 
    END IF

    L = 1 ; 
    R = N ;
    DO
       Mid = (L + R)/2 ;
       mval(1) = minval(A(L:Mid)) ;

       IF ( Mid < R ) THEN
          mval(2) = minval(A(Mid+1:R)) ;
       ELSE
          mval(2) = minval(A(Mid:R)) ;
       END IF

       IF ( mval(1) < mval(2) ) THEN
          R = Mid ;        
          idx = L ;
       ELSE IF ( mval(2) < mval(1) ) THEN
          L = Mid + 1 ;        
          idx = L ;
       ELSE
          IF ( R == L ) THEN
             idx = L ; 
             val = mval(1)
             EXIT ; 
          ELSE
             print*, "Warning: minval1(): "
             print*, "  multiple entries with identical minimum value" ; 
             
             R = Mid ;
             idx = L ; 
             CONTINUE ;
          END IF
       END IF
    END DO

    RETURN ;
  END SUBROUTINE IMinVal1
  
  !c Find the maximum value in A(1:N)
  SUBROUTINE IMaxVal1( val, idx, A, N )
    IMPLICIT NONE

    INTEGER:: N
    INTEGER:: val, idx, A(N)

    INTEGER:: L, R, Mid, mval(2)

    IF ( N == 1 ) THEN
       val = A(1) ;
       return ; 
    END IF
    
    L = 1 ; 
    R = N ;
    DO
       Mid = (L + R)/2 ;
       mval(1) = minval(A(L:Mid)) ;
       
       IF ( Mid < R ) THEN
          mval(2) = minval(A(Mid+1:R)) ;
       ELSE
          mval(2) = minval(A(Mid:R)) ;
       END IF
       
       IF ( mval(1) > mval(2) ) THEN
          R = Mid ;        
          idx = L ;
       ELSE IF ( mval(2) > mval(1) ) THEN
          L = Mid + 1 ;        
          idx = L ;
       ELSE
          IF ( R == L ) THEN
             idx = L ; 
             val = mval(1)
             EXIT ; 
          ELSE
             print*, "Warning: minval1(): "
             print*, "  multiple entries with identical maximum value" ; 
             
             R = Mid ;
             idx = L ; 
             CONTINUE ;
          END IF
       END IF
    END DO
    
    RETURN ;
  END SUBROUTINE IMaxVal1
!c========================= Hash sort =====================c!
  
END MODULE AltMathUtils
