MODULE UTILS

  USE TypeDefs

  IMPLICIT NONE

  INTERFACE PointDataInit
     MODULE PROCEDURE PointDataInitNRXF, PointDataInitUT, PointDataInitUTNRXF
  END INTERFACE PointDataInit

  INTERFACE ResizePointData
     MODULE PROCEDURE ResizePointDataNRXF, ResizePointDataUT
  END INTERFACE ResizePointData

  INTERFACE ExpandIntarray
     MODULE PROCEDURE Expand1IntArray, Expand2IntArray
  END INTERFACE ExpandIntarray

  CONTAINS

    SUBROUTINE Expand1IntArray(intarr,newsize_in,fill_in)

      INTEGER, ALLOCATABLE :: intarr(:), workarr(:)
      INTEGER, OPTIONAL :: newsize_in,fill_in
      INTEGER :: newsize, oldsize, fill

      oldsize = SIZE(intarr)
      IF(PRESENT(newsize_in)) THEN
        newsize = newsize_in
      ELSE
        newsize =  oldsize * 2
      END IF

      fill = 0
      IF(PRESENT(fill_in)) fill = fill_in

      ! ALLOCATE(workarr(oldsize))
      ! workarr(1:oldsize) = intarr(1:oldsize)
      ! DEALLOCATE(intarr)
      ALLOCATE(workarr(newsize))
      workarr = fill
      workarr(1:oldsize) = intarr(1:oldsize)

      !FORTRAN 2003 feature...
      CALL MOVE_ALLOC(workarr, intarr)

      IF(DebugMode) PRINT *,'DEBUG: Done move alloc'

    END SUBROUTINE Expand1IntArray

    SUBROUTINE Expand2IntArray(intarr,newsize_in,fill_in)

      INTEGER, ALLOCATABLE :: intarr(:,:), workarr(:,:)
      INTEGER, OPTIONAL :: newsize_in,fill_in
      INTEGER :: newsize, oldsize, dim1size, fill

      oldsize = SIZE(intarr,2)

      IF(PRESENT(newsize_in)) THEN
        newsize = newsize_in
      ELSE
        newsize =  oldsize * 2
      END IF

      fill = 0
      IF(PRESENT(fill_in)) fill = fill_in

      dim1size = SIZE(intarr,1)

      ! ALLOCATE(workarr(oldsize))
      ! workarr(1:oldsize) = intarr(1:oldsize)
      ! DEALLOCATE(intarr)
      ALLOCATE(workarr(dim1size,newsize))
      workarr = fill

      workarr(:,1:oldsize) = intarr(:,1:oldsize)

      !FORTRAN 2003 feature...
      CALL MOVE_ALLOC(workarr, intarr)

      IF(DebugMode) PRINT *,'DEBUG: Done move alloc'

    END SUBROUTINE Expand2IntArray

    SUBROUTINE PointDataInitNRXF(NRXF, n, partexpand, arrsize)
      TYPE(NRXF_T), TARGET :: NRXF
      INTEGER :: n,n_tot
      REAL*8,OPTIONAL :: partexpand
      INTEGER,OPTIONAL :: arrsize

      IF(PRESENT(partexpand) .EQV. PRESENT(arrsize)) THEN
        PRINT *,myid,' PointDataInitNRXF: Programming error, &
             &should provide one of partexpand and arrsize!'
        STOP
      END IF

      IF(PRESENT(partexpand)) THEN
        n_tot = n + CEILING(n * partexpand)
      ELSE
        n_tot = arrsize
      END IF

      IF(ALLOCATED(NRXF%A)) THEN
        PRINT *, "Programming error: abuse of PointDataInitNRXF"
        STOP
      END IF

      NRXF%mstrt = 1
      NRXF%cstrt = 1 + n
      NRXF%pstrt = 1 + n !no connected particles, initially

      NRXF%NN = n
      NRXF%NC = 0
      NRXF%NP = 0

      ALLOCATE(NRXF%A(3,n_tot),NRXF%PartInfo(2,n_tot))

      NRXF%PartInfo(:,:) = -1 !-1 = no point

      NRXF%A = 0.0
      NRXF%M => NRXF%A(:,1:n)
      ! NRXF%P => NRXF%A(:,n+1:UBOUND(NRXF%A,2))

      IF(DebugMode) PRINT *,'Debug,nrxf init n, ntot: ',n,n_tot, SIZE(NRXF%A), NRXF%NN, UBOUND(NRXF%A,2)

    END SUBROUTINE PointDataInitNRXF

    SUBROUTINE PointDataInitUT(UT, n, partexpand, arrsize)
      TYPE(UT_T), TARGET :: UT
      INTEGER :: n,n_tot
      REAL*8, OPTIONAL :: partexpand
      INTEGER,OPTIONAL :: arrsize

      IF(PRESENT(partexpand) .EQV. PRESENT(arrsize)) THEN
        PRINT *,myid,' PointDataInitNRXF: Programming error, &
             &should provide one of partexpand and arrsize!'
        STOP
      END IF

      IF(PRESENT(partexpand)) THEN
        n_tot = n + CEILING(n * partexpand)
      ELSE
        n_tot = arrsize
      END IF

      IF(ALLOCATED(UT%A)) THEN
        PRINT *, "Programming error: abuse of PointDataInitUT"
        STOP
      END IF

      ALLOCATE(UT%A(6*n_tot))
      UT%A = 0.0
      UT%M => UT%A(1:6*n)
      ! UT%P => UT%A((6*n+1):UBOUND(UT%A,1)) <- not used

      IF(DebugMode) PRINT *,'Debug,ut init n, ntot: ',n,n_tot, SIZE(UT%A), SIZE(UT%M), SIZE(UT%P)

    END SUBROUTINE PointDataInitUT

    SUBROUTINE PointDataInitUTNRXF(UT, NRXF)
      TYPE(UT_T), TARGET :: UT
      TYPE(NRXF_T) :: NRXF
      INTEGER :: n_tot, NN

      n_tot = SIZE(NRXF%A,2)

      IF(ALLOCATED(UT%A)) THEN
        PRINT *, "Programming error: abuse of PointDataInitUT"
        STOP
      END IF

      ALLOCATE(UT%A(6*n_tot))
      UT%A = 0.0
      UT%M => UT%A(1:6*NRXF%NN)
      ! UT%P => UT%A((6*n+1):UBOUND(UT%A,1)) <- not used

      IF(DebugMode) PRINT *,myid,' Debug, UT init NN, ntot: ',NRXF%NN,n_tot, SIZE(UT%A), &
           SIZE(UT%M), SIZE(UT%P)

    END SUBROUTINE PointDataInitUTNRXF

    SUBROUTINE InvPartInfoInit(InvPartInfo,neighparts,initsize_in)
      TYPE(InvPartInfo_t),ALLOCATABLE :: InvPartInfo(:)
      LOGICAL :: neighparts(0:)
      INTEGER, OPTIONAL :: initsize_in
      !------------------------------
      INTEGER :: initsize,i,j,neighcount

      neighcount = COUNT(neighparts)

      initsize = 100
      IF(PRESENT(initsize_in)) initsize = initsize_in

      ALLOCATE(InvPartInfo(0:ntasks-1))
      DO i=0,ntasks-1
        InvPartInfo(i) % NID = i
        InvPartInfo(i) % ccount = 0
        InvPartInfo(i) % pcount = 0
        InvPartInfo(i) % sccount = 0
        InvPartInfo(i) % spcount = 0
      END DO

      DO j=0,ntasks-1
        IF(.NOT. neighparts(j)) CYCLE
        ALLOCATE(InvPartInfo(j) % ConnIDs(initsize),&
             InvPartInfo(j) % ConnLocs(initsize),&
             InvPartInfo(j) % ProxIDs(initsize),&
             InvPartInfo(j) % ProxLocs(initsize),&
             InvPartInfo(j) % SConnIDs(initsize),&
             InvPartInfo(j) % SProxIDs(initsize))
        InvPartInfo(j) % ConnIDs = -1
        InvPartInfo(j) % ConnLocs = -1
        InvPartInfo(j) % ProxIDs = -1
        InvPartInfo(j) % ProxLocs = -1
        InvPartInfo(j) % SProxIDs = -1
        InvPartInfo(j) % SConnIDs = -1
      END DO

    END SUBROUTINE InvPartInfoInit

    SUBROUTINE NewInvPartInfo(InvPartInfo, nid, initsize_in)
      TYPE(InvPartInfo_t) :: InvPartInfo(0:)
      INTEGER :: nid
      INTEGER, OPTIONAL :: initsize_in
      !--------------------
      INTEGER :: initsize

      initsize = 100
      IF(PRESENT(initsize_in)) initsize = initsize_in

      ALLOCATE(InvPartInfo(nid) % ConnIDs(initsize),&
           InvPartInfo(nid) % ConnLocs(initsize),&
           InvPartInfo(nid) % ProxIDs(initsize),&
           InvPartInfo(nid) % ProxLocs(initsize),&
           InvPartInfo(nid) % SProxIDs(initsize),&
           InvPartInfo(nid) % SConnIDs(initsize))

      InvPartInfo(nid) % ConnIDs = -1
      InvPartInfo(nid) % ConnLocs = -1
      InvPartInfo(nid) % ProxIDs = -1
      InvPartInfo(nid) % ProxLocs = -1
      InvPartInfo(nid) % SConnIDs = -1
      InvPartInfo(nid) % SProxIDs = -1

    END SUBROUTINE NewInvPartInfo

    SUBROUTINE ResizePointDataNRXF(NRXF,scale,UT,UTM,do_M,do_C,do_P)

      TYPE(NRXF_T), TARGET :: NRXF
      TYPE(UT_T), OPTIONAL, TARGET :: UT, UTM
      REAL*8 :: scale
      LOGICAL, OPTIONAL :: do_M, do_C, do_P
      !---------------------
      REAL*8, ALLOCATABLE :: work_arr(:,:),work_arr1(:)
      INTEGER, ALLOCATABLE :: work_int(:,:)
      INTEGER :: a_oldsize,m_oldsize,c_oldsize,p_oldsize
      INTEGER :: a_newsize,m_newsize,c_newsize,p_newsize
      INTEGER :: cstrt_new, pstrt_new,cstrt_old, pstrt_old
      LOGICAL :: doM,doC,doP
      
      IF(PRESENT(do_M)) THEN
        doM = do_M
      ELSE
        doM = .FALSE.
      END IF

      IF(PRESENT(do_C)) THEN
        doC = do_C
      ELSE
        doC = .FALSE.
      END IF

      IF(PRESENT(do_P)) THEN
        doP = do_P
      ELSE
        doP = .TRUE.
      END IF

      IF(scale < 1.0) THEN
        PRINT *,'ERROR: NRXF size reduction not yet implemented, sorry!'
        STOP
      END IF

      !ISSUE in here - repeated resizing produces wacky results - see LOG_5492906.sdb in test_metis
      a_oldsize = SIZE(NRXF%A,2)
      m_oldsize = NRXF%cstrt - NRXF%mstrt
      IF(NRXF%NN > m_oldsize) THEN
        IF(DebugMode) PRINT *,myid,' WARNING: ResizePointDataNRXF: NRXF%NN overlaps NC'
        m_oldsize = NRXF%NN
        NRXF%cstrt = NRXF%NN+1
      END IF

      c_oldsize = NRXF%pstrt - NRXF%cstrt
      IF(NRXF%NC > c_oldsize) THEN
        IF(DebugMode) PRINT *,myid,' WARNING: ResizePointDataNRXF: NRXF%NC overlaps NP'
        c_oldsize = NRXF%NC
        NRXF%pstrt = NRXF%cstrt + NRXF%NC
      END IF

      p_oldsize = a_oldsize - NRXF%pstrt + 1
      IF(NRXF%NP > p_oldsize) THEN
        IF(DebugMode) PRINT *,myid,' WARNING: ResizePointDataNRXF: NRXF%NP overlaps array end'
        p_oldsize = NRXF%NP
      END IF

      IF(SIZE(NRXF%M,2) /= m_oldsize) THEN
        PRINT *, "ERROR: NRXF%M wrong size in ResizePointData"
        STOP
      END IF

      IF(doM) THEN
        m_newsize = MAX(CEILING(m_oldsize*scale),100)
      ELSE
        m_newsize = m_oldsize
      END IF

      IF(doC) THEN
        c_newsize = MAX(CEILING(c_oldsize*scale),100)
      ELSE
        c_newsize = c_oldsize
      END IF

      IF(doP) THEN
        p_newsize = MAX(CEILING(p_oldsize*scale),100)
      ELSE
        p_newsize = p_oldsize
      END IF

      a_newsize = m_newsize + c_newsize + p_newsize
      
      cstrt_new = m_newsize + 1
      pstrt_new = m_newsize + c_newsize + 1

      IF(DebugMode) PRINT *,myid,' debug resize nrxf: ',a_oldsize, m_oldsize, c_oldsize, p_oldsize,&
            'strts: ', NRXF%cstrt, NRXF%pstrt, 'new: ',a_newsize, m_newsize, c_newsize, p_newsize,&
            'strts: ', cstrt_new, pstrt_new

      IF(a_oldsize < 0 .OR. m_oldsize < 0 .OR. c_oldsize < 0 .OR. p_oldsize < 0 .OR. &
           a_newsize < 0 .OR. m_newsize < 0 .OR. c_newsize < 0 .OR. p_newsize < 0) THEN
        PRINT *,'ResizePointDataNRXF: Programming error, negative sizes!'
      END IF

      ALLOCATE(work_arr(3,a_newsize),work_int(2,a_newsize))
      work_arr = 0.0
      work_int = -1

      work_arr(:,1:m_oldsize) = NRXF%A(:,1:m_oldsize)
      work_arr(:,cstrt_new : cstrt_new+c_oldsize-1) = NRXF % A(:,NRXF%cstrt:NRXF%cstrt+c_oldsize-1)
      work_arr(:,pstrt_new : pstrt_new+p_oldsize-1) = NRXF % A(:,NRXF%pstrt:NRXF%pstrt+p_oldsize-1)

      work_int(:,1:m_oldsize) = NRXF%PartInfo(:,1:m_oldsize)
      work_int(:,cstrt_new : cstrt_new+c_oldsize-1) = NRXF%PartInfo(:,NRXF%cstrt:NRXF%cstrt+c_oldsize-1)
      work_int(:,pstrt_new : pstrt_new+p_oldsize-1) = NRXF%PartInfo(:,NRXF%pstrt:NRXF%pstrt+p_oldsize-1)

      CALL MOVE_ALLOC(work_arr, NRXF%A)
      CALL MOVE_ALLOC(work_int, NRXF%PartInfo)

      cstrt_old = NRXF%cstrt
      pstrt_old = NRXF%pstrt

      NRXF%cstrt = cstrt_new
      NRXF%pstrt = pstrt_new
      
      NRXF%M => NRXF%A(:,1:m_newsize)
      ! NRXF%P => NRXF%A(:,m_newsize+1 : a_newsize)

      IF(PRESENT(UT)) THEN
        ALLOCATE(work_arr1(a_newsize*6))
        work_arr1 = 0.0

        work_arr1(1:m_oldsize*6) = UT%A(1:m_oldsize*6)
        work_arr1((cstrt_new-1)*6 - 1 : (cstrt_new-1)*6 + c_oldsize*6) = &
             UT%A((cstrt_old-1)*6 - 1 : (cstrt_old-1)*6 + c_oldsize*6)
        work_arr1((pstrt_new-1)*6 - 1 : (pstrt_new-1)*6 + p_oldsize*6) = &
             UT%A((pstrt_old-1)*6 - 1 : (pstrt_old-1)*6 + p_oldsize*6)

        CALL MOVE_ALLOC(work_arr1, UT%A)
        UT%M => UT%A(1:m_newsize*6)
      END IF

      IF(PRESENT(UTM)) THEN
        ALLOCATE(work_arr1(a_newsize*6))
        work_arr1 = 0.0

        work_arr1(1:m_oldsize*6) = UTM%A(1:m_oldsize*6)
        work_arr1((cstrt_new-1)*6 - 1 : (cstrt_new-1)*6 + c_oldsize*6) = &
             UTM%A((cstrt_old-1)*6 - 1 : (cstrt_old-1)*6 + c_oldsize*6)
        work_arr1((pstrt_new-1)*6 - 1 : (pstrt_new-1)*6 + p_oldsize*6) = &
             UTM%A((pstrt_old-1)*6 - 1 : (pstrt_old-1)*6 + p_oldsize*6)

        CALL MOVE_ALLOC(work_arr1, UTM%A)
        UTM%M => UTM%A(1:m_newsize*6)
      END IF

      IF(DebugMode) PRINT *,myid,' debug2 resize nrxf: ',SIZE(NRXF%A,2),NRXF%cstrt, NRXF%pstrt

    END SUBROUTINE ResizePointDataNRXF

    SUBROUTINE ResizePointDataUT(UT,scale,do_M,do_C,do_P)

      TYPE(UT_T), TARGET :: UT
      REAL*8 :: scale
      LOGICAL, OPTIONAL :: do_M, do_C, do_P
      !---------------------
      REAL*8, ALLOCATABLE :: work_arr(:)
      INTEGER :: a_oldsize,m_oldsize,p_oldsize
      INTEGER :: a_newsize,m_newsize,p_newsize
      LOGICAL :: doM,doC,doP
      
      IF(PRESENT(do_M)) THEN
        doM = do_M
      ELSE
        doM = .FALSE.
      END IF

      IF(PRESENT(do_C)) THEN
        doC = do_C
      ELSE
        doC = .FALSE.
      END IF

      IF(PRESENT(do_P)) THEN
        doP = do_P
      ELSE
        doP = .TRUE.
      END IF

      IF(scale < 1.0) THEN
        PRINT *,'UT size reduction not yet implemented, sorry!'
        STOP
      END IF

      a_oldsize = SIZE(UT%A)/6
      m_oldsize = SIZE(UT%M)/6
      p_oldsize = SIZE(UT%P)/6

      IF(doP) THEN
        p_newsize = MAX(CEILING(p_oldsize*scale),100)
      ELSE
        p_newsize = p_oldsize
      END IF

      IF(doM) THEN
        m_newsize = MAX(CEILING(m_oldsize*scale),100)
      ELSE
        m_newsize = m_oldsize
      END IF

      a_newsize = m_newsize + p_newsize

      IF(DebugMode) PRINT *,myid,' debug resize ut: ',a_oldsize, m_oldsize, p_oldsize,&
           'new: ',a_newsize, m_newsize, p_newsize

      a_newsize = a_newsize * 6
      m_newsize = m_newsize * 6
      p_newsize = p_newsize * 6

      ALLOCATE(work_arr(a_newsize))

      work_arr = 0.0
      work_arr(1:m_oldsize) = UT % M(1:m_oldsize)
      work_arr(m_newsize+1 : m_newsize+p_oldsize) = UT % P(1:p_oldsize)

      CALL MOVE_ALLOC(work_arr, UT%A)

      UT%M => UT%A(1:m_newsize)
      UT%P => UT%A(m_newsize+1 : a_newsize)

    END SUBROUTINE ResizePointDataUT

    !Subroutine to intermittently clear out old (presumably prox) points
    !which are no longer near our partition, therefore no longer sent
    !These are marked with NRXF%PartInfo(:,i) == -1
    SUBROUTINE ClearOldParticles(NRXF,UT,UTM,InvPartInfo)
      TYPE(NRXF_t) :: NRXF
      TYPE(UT_t) :: UT,UTM
      TYPE(InvPartInfo_t) :: InvPartInfo(0:)
      !-------------------------
      INTEGER :: i,j,pstrt, pend, NP,upper, newNP
      LOGICAL, ALLOCATABLE :: UTValid(:)

      NP = NRXF%NP
      pstrt = NRXF%pstrt
      pend = pstrt+NP-1

      upper = UBOUND(NRXF%A,2)
      newNP = COUNT(NRXF%PartInfo(1,pstrt:pend) /= -1)

      ALLOCATE(UTValid(SIZE(UT%A)))
      UTValid = .TRUE.
      DO i=pstrt, upper
        IF(NRXF%PartInfo(1,i) == -1) UTValid(i*6-5:i*6) = .FALSE.
      END DO

      IF(DebugMode) PRINT *,myid,' debug clear: pstrt,pend,np,upper,newNP :',pstrt,pend,np,upper,newNP

      UT%A(pstrt*6 - 5 : (pstrt + newNP - 1)*6) = PACK(UT%A(pstrt*6-5 : upper*6),&
           UTValid(pstrt*6-5:upper*6))

      UTM%A(pstrt*6 - 5 : (pstrt + newNP - 1)*6) = PACK(UTM%A(pstrt*6-5 : upper*6),&
           UTValid(pstrt*6-5:upper*6))

      DO i=1,3
        NRXF%A(i,pstrt : pstrt + newNP - 1) = PACK(NRXF%A(i,pstrt : upper),&
             NRXF%PartInfo(1,pstrt : upper) /= -1)
      END DO
      DO i=1,2
        NRXF%PartInfo(i,pstrt : pstrt + newNP - 1) = PACK(NRXF%PartInfo(i,pstrt : upper),&
             NRXF%PartInfo(1,pstrt : upper) /= -1)
      END DO

      NRXF%NP = newNP

    END SUBROUTINE ClearOldParticles

    !A subroutine to quicksort an int array (arr1)
    ! while also sorting arr2 by arr1
    SUBROUTINE sort_int2(arr1,arr2,n)
      IMPLICIT NONE
      INTEGER :: arr1(:), arr2(:),n
      IF(n <= 1) RETURN
      CALL sort_int2_r(arr1,arr2,1,n)
    END SUBROUTINE sort_int2

    !Does the actual sorting - needs a wrapper to avoid
    !having to pass the extents first
    RECURSIVE SUBROUTINE sort_int2_r(arr1,arr2,start,fin)
      IMPLICIT NONE
      INTEGER :: arr1(:), arr2(:)
      INTEGER :: start,fin,pivot,pivot_val,i,j,hold
      !-----------------------------------

      i=start
      j=fin
      pivot = (fin + start) / 2
      pivot_val = arr1(pivot)

      DO WHILE(.TRUE.)
        DO WHILE(arr1(i) < pivot_val)
          i=i+1
        END DO
        DO WHILE(arr1(j) > pivot_val)
          j=j-1
        END DO
        !    PRINT *,'debug ij: ',i,j,pivot
        IF(i>=j) EXIT

        hold = arr1(i)
        arr1(i) = arr1(j)
        arr1(j) = hold

        hold = arr2(i)
        arr2(i) = arr2(j)
        arr2(j) = hold

        i=i+1
        j=j-1
      END DO

      IF(j+1 < fin) THEN
        CALL sort_int2_r(arr1,arr2,j+1,fin)
      END IF
      IF(start < i-1) THEN
        CALL sort_int2_r(arr1,arr2,start, i-1)
      END IF

      RETURN
    END SUBROUTINE sort_int2_r

END MODULE UTILS
