MODULE UTILS

  USE TypeDefs
  USE INOUT

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

  INTERFACE ExpandRealarray
     MODULE PROCEDURE Expand1RealArray, Expand2RealArray
  END INTERFACE ExpandRealarray

  CONTAINS

    SUBROUTINE Expand1IntArray(intarr,newsize_in,fill_in)

      INTEGER, ALLOCATABLE :: intarr(:), workarr(:)
      INTEGER, OPTIONAL :: newsize_in,fill_in
      INTEGER :: newsize, oldsize, fill

      IF(.NOT. ALLOCATED(intarr)) THEN
        PRINT *,myid,' Error: Expand1IntArray received an unallocated array!'
        STOP
      END IF

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
      DEALLOCATE(intarr) !Cray compiler bug?
      CALL MOVE_ALLOC(workarr, intarr)

      IF(DebugMode) PRINT *,'DEBUG: Done move alloc'

    END SUBROUTINE Expand1IntArray

    SUBROUTINE Expand2IntArray(intarr,newsize_in,fill_in)

      INTEGER, ALLOCATABLE :: intarr(:,:), workarr(:,:)
      INTEGER, OPTIONAL :: newsize_in,fill_in
      INTEGER :: newsize, oldsize, dim1size, fill

      IF(.NOT. ALLOCATED(intarr)) THEN
        PRINT *,myid,' Error: Expand2IntArray received an unallocated array!'
        STOP
      END IF

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
      DEALLOCATE(intarr) !Cray compiler bug?
      CALL MOVE_ALLOC(workarr, intarr)

      IF(DebugMode) PRINT *,'DEBUG: Done move alloc'

    END SUBROUTINE Expand2IntArray

    SUBROUTINE Expand1RealArray(realarr,newsize_in,fill_in)

      REAL(KIND=dp), ALLOCATABLE :: realarr(:), workarr(:)
      INTEGER, OPTIONAL :: newsize_in,fill_in
      INTEGER :: newsize, oldsize, fill

      IF(.NOT. ALLOCATED(realarr)) THEN
        PRINT *,myid,' Error: Expand1RealArray received an unallocated array!'
        STOP
      END IF

      oldsize = SIZE(realarr)
      IF(PRESENT(newsize_in)) THEN
        newsize = newsize_in
      ELSE
        newsize =  oldsize * 2
      END IF

      fill = 0
      IF(PRESENT(fill_in)) fill = fill_in

      ! ALLOCATE(workarr(oldsize))
      ! workarr(1:oldsize) = realarr(1:oldsize)
      ! DEALLOCATE(realarr)
      ALLOCATE(workarr(newsize))
      workarr = fill
      workarr(1:oldsize) = realarr(1:oldsize)

      !FORTRAN 2003 feature...
      DEALLOCATE(realarr) !Cray compiler bug?
      CALL MOVE_ALLOC(workarr, realarr)

      IF(DebugMode) PRINT *,'DEBUG: Done move alloc'

    END SUBROUTINE Expand1RealArray


    SUBROUTINE Expand2RealArray(realarr,newsize_in,fill_in)

      REAL(KIND=dp), ALLOCATABLE :: realarr(:,:), workarr(:,:)
      INTEGER, OPTIONAL :: newsize_in,fill_in
      INTEGER :: newsize, oldsize, dim1size, fill

      IF(.NOT. ALLOCATED(realarr)) THEN
        PRINT *,myid,' Error: Expand2RealArray received an unallocated array!'
        STOP
      END IF

      oldsize = SIZE(realarr,2)

      IF(PRESENT(newsize_in)) THEN
        newsize = newsize_in
      ELSE
        newsize =  oldsize * 2
      END IF

      fill = 0
      IF(PRESENT(fill_in)) fill = fill_in

      dim1size = SIZE(realarr,1)

      ! ALLOCATE(workarr(oldsize))
      ! workarr(1:oldsize) = intarr(1:oldsize)
      ! DEALLOCATE(intarr)
      ALLOCATE(workarr(dim1size,newsize))
      workarr = fill

      workarr(:,1:oldsize) = realarr(:,1:oldsize)

      !FORTRAN 2003 feature...
      DEALLOCATE(realarr) !Cray compiler bug?
      CALL MOVE_ALLOC(workarr, realarr)

      IF(DebugMode) PRINT *,'DEBUG: Done move alloc'

    END SUBROUTINE Expand2RealArray

    SUBROUTINE PointDataInitNRXF(NRXF, n, partexpand, arrsize)
      TYPE(NRXF_T), TARGET :: NRXF
      INTEGER :: n,n_tot
      REAL(KIND=dp),OPTIONAL :: partexpand
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

      ALLOCATE(NRXF%A(3,n_tot),NRXF%PartInfo(2,n_tot),NRXF%GID(n_tot))

      NRXF%PartInfo(:,:) = -1 !-1 = no point

      NRXF%A = 0.0
      NRXF%M => NRXF%A(:,1:n)
      ! NRXF%P => NRXF%A(:,n+1:UBOUND(NRXF%A,2))

      IF(DebugMode) PRINT *,'Debug,nrxf init n, ntot: ',n,n_tot, SIZE(NRXF%A), NRXF%NN, UBOUND(NRXF%A,2)

    END SUBROUTINE PointDataInitNRXF

    SUBROUTINE PointDataInitUT(UT, n, partexpand, arrsize)
      TYPE(UT_T), TARGET :: UT
      INTEGER :: n,n_tot
      REAL(KIND=dp), OPTIONAL :: partexpand
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
      INTEGER :: initsize,i,j,neighcount,nopartitions

      neighcount = COUNT(neighparts)

      nopartitions = SIZE(neighparts)

      initsize = 100
      IF(PRESENT(initsize_in)) initsize = initsize_in

      ALLOCATE(InvPartInfo(0:nopartitions-1))
      DO i=0,nopartitions-1
        InvPartInfo(i) % NID = i
        InvPartInfo(i) % ccount = 0
        InvPartInfo(i) % pcount = 0
        InvPartInfo(i) % sccount = 0
        InvPartInfo(i) % spcount = 0
      END DO

      DO j=0,nopartitions-1
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
      REAL(KIND=dp) :: scale
      LOGICAL, OPTIONAL :: do_M, do_C, do_P
      !---------------------
      REAL(KIND=dp), ALLOCATABLE :: work_arr(:,:),work_arr1(:)
      INTEGER, ALLOCATABLE :: work_int(:,:),work_gid(:)
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
        !TODO - this would lead to invalid array reads, I think
        PRINT *,myid,' WARNING: ResizePointDataNRXF: NRXF%NP overlaps array end'
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

      ALLOCATE(work_arr(3,a_newsize),work_int(2,a_newsize), work_gid(a_newsize))
      work_arr = 0.0
      work_int = -1
      work_gid = 0

      work_arr(:,1:m_oldsize) = NRXF%A(:,1:m_oldsize)
      work_arr(:,cstrt_new : cstrt_new+c_oldsize-1) = NRXF % A(:,NRXF%cstrt:NRXF%cstrt+c_oldsize-1)
      work_arr(:,pstrt_new : pstrt_new+p_oldsize-1) = NRXF % A(:,NRXF%pstrt:NRXF%pstrt+p_oldsize-1)

      work_int(:,1:m_oldsize) = NRXF%PartInfo(:,1:m_oldsize)
      work_int(:,cstrt_new : cstrt_new+c_oldsize-1) = NRXF%PartInfo(:,NRXF%cstrt:NRXF%cstrt+c_oldsize-1)
      work_int(:,pstrt_new : pstrt_new+p_oldsize-1) = NRXF%PartInfo(:,NRXF%pstrt:NRXF%pstrt+p_oldsize-1)

      work_gid(1:m_oldsize) = NRXF%GID(1:m_oldsize)
      work_gid(cstrt_new : cstrt_new+c_oldsize-1) = NRXF%GID(NRXF%cstrt:NRXF%cstrt+c_oldsize-1)
      work_gid(pstrt_new : pstrt_new+p_oldsize-1) = NRXF%GID(NRXF%pstrt:NRXF%pstrt+p_oldsize-1)

      DEALLOCATE(NRXF%A) !Cray compiler bug?
      CALL MOVE_ALLOC(work_arr, NRXF%A)
      DEALLOCATE(NRXF%PartInfo) !Cray compiler bug?
      CALL MOVE_ALLOC(work_int, NRXF%PartInfo)
      DEALLOCATE(NRXF%GID)
      CALL MOVE_ALLOC(work_gid, NRXF%GID)

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

        DEALLOCATE(UT%A) !Cray compiler bug?
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

        DEALLOCATE(UTM%A) !Cray compiler bug?
        CALL MOVE_ALLOC(work_arr1, UTM%A)
        UTM%M => UTM%A(1:m_newsize*6)
      END IF

      IF(DebugMode) PRINT *,myid,' debug2 resize nrxf: ',SIZE(NRXF%A,2),NRXF%cstrt, NRXF%pstrt

    END SUBROUTINE ResizePointDataNRXF

    SUBROUTINE ResizePointDataUT(UT,scale,do_M,do_C,do_P)

      TYPE(UT_T), TARGET :: UT
      REAL(KIND=dp) :: scale
      LOGICAL, OPTIONAL :: do_M, do_C, do_P
      !---------------------
      REAL(KIND=dp), ALLOCATABLE :: work_arr(:)
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

      DEALLOCATE(UT%A) !Cray compiler bug?
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
      TYPE(InvPartInfo_t), TARGET :: InvPartInfo(0:)
      !-------------------------
      INTEGER :: i,j,pstrt, pend, newPend, NP, newNP, upper, rmcount
      INTEGER, ALLOCATABLE :: locLUT(:)
      LOGICAL, ALLOCATABLE :: UTValid(:)
      TYPE(InvPartInfo_t), POINTER :: IPI=>NULL()

      NP = NRXF%NP
      pstrt = NRXF%pstrt
      pend = pstrt+NP-1

      upper = UBOUND(NRXF%A,2)
      newNP = COUNT(NRXF%PartInfo(1,pstrt:pend) /= -1)
      newPend = pstrt + newNP - 1

      IF(DebugMode) PRINT *,myid,' debug clear: pstrt,pend,np,upper,newNP :',pstrt,pend,np,upper,newNP

      IF(newNP == NP) THEN
        IF(DebugMode) PRINT *,myid,' debug, didnt need to clear any removed points.'
        RETURN
      END IF

      !Create a logical array to enable PACKing of UT
      ALLOCATE(UTValid(SIZE(UT%A)))
      UTValid = .TRUE.
      DO i=pstrt, upper
        IF(NRXF%PartInfo(1,i) == -1) UTValid(i*6-5:i*6) = .FALSE.
      END DO

      !This lookup table gives the updated index (in NRXF%A,PartInfo)
      !once the invalid values have been removed
      ALLOCATE(locLUT(pstrt:pend))
      locLUT = 0
      rmcount = 0
      DO i=pstrt,pend
        IF(NRXF%PartInfo(1,i) == -1) THEN
          rmcount = rmcount + 1
          locLUT(i) = -1
        ELSE
          locLUT(i) = i-rmcount
        END IF
      END DO

      !TODO - most of these 'upper's should be newPend - but i'm scared it might 
      !break something...
      !shift, and delete (= 0.0) beyond newPend
      UT%A(pstrt*6 - 5 : newPend*6) = PACK(UT%A(pstrt*6-5 : upper*6),&
           UTValid(pstrt*6-5:upper*6))
      UT%A(newPend*6 + 1) = 0.0

      UTM%A(pstrt*6 - 5 : newPend*6) = PACK(UTM%A(pstrt*6-5 : upper*6),&
           UTValid(pstrt*6-5:upper*6))
      UTM%A(newPend*6 + 1) = 0.0

      DO i=1,3
        NRXF%A(i,pstrt : newPend) = PACK(NRXF%A(i,pstrt : upper),&
             NRXF%PartInfo(1,pstrt : upper) /= -1)
      END DO
      NRXF%A(:,newPend+1:upper) = 0.0

      !shift, and delete (= -1) beyond newPend
      NRXF%PartInfo(2,pstrt : newPend) = PACK(NRXF%PartInfo(2,pstrt : upper),&
           NRXF%PartInfo(1,pstrt : upper) /= -1)
      NRXF%PartInfo(1,pstrt : newPend) = PACK(NRXF%PartInfo(1,pstrt : upper),&
           NRXF%PartInfo(1,pstrt : upper) /= -1)

      NRXF%PartInfo(:,newPend+1:upper) = -1

      NRXF%NP = newNP


      !Correct locations in InvPartInfo
      DO i=0,ntasks-1
        IPI => InvPartInfo(i)
        DO j=1,IPI%Pcount
          IF(locLUT(IPI % ProxLocs(j)) <= 0) THEN
            PRINT *,myid," Programming Error in locLUT"
            STOP
          END IF
          IPI % ProxLocs(j) = locLUT(IPI % ProxLocs(j))
        END DO
      END DO
    END SUBROUTINE ClearOldParticles

    !Checks particle position and velocity for anything suspicious
    !Freezes particles which are outside the domain, determines outliers 
    !for the purpose of BBox calculation
    SUBROUTINE CheckSolution(NRXF,UT,UTP,NN,NTOT,NANS,EFS,grid_bbox,DT,MAXUT,Lost,Outlier)

      INTEGER ::  NN, NTOT
      INTEGER :: NANS(:,:)
      REAL(KIND=dp), ALLOCATABLE :: EFS(:),UTP(:)
      REAL(KIND=dp) :: DT,grid_bbox(4),MAXUT
      TYPE(UT_t) :: UT
      TYPE(NRXF_t) :: NRXF
      LOGICAL, ALLOCATABLE :: Lost(:), Outlier(:)
      !-----------------------------------
      INTEGER :: i,j,ierr
      REAL(KIND=dp) :: BBox(6),bbox_vol,bbox_vols(ntasks),med_vol,max_gap(3),gap(3)
      REAL(KIND=dp) :: outlier_gap_prop, outlier_arr_prop, freezer(3)
      REAL(KIND=dp) :: X,Y,Z,dx,dy,dz,disp, disp_limit
      REAL(KIND=dp) :: minx,maxx,miny,maxy,minz,maxz,midx,midy,midz,speed_limit, speed
      REAL(KIND=dp) :: dist(3,NN), pos(3,NN), sigma_dist(3,NN),mean_dist(3),std_dev(3)
      INTEGER :: dist_sort(3,NN), pos_sort(3,NN)
      LOGICAL :: d_outlier(3,NN), poutlier(NN),check_outliers,too_fast(NN)


      IF(.NOT. (ALLOCATED(Lost) .AND. ALLOCATED(Outlier))) THEN
        PRINT *,myid," IsLost or IsOutlier not allocated!"
        STOP
      END IF

      freezer(:) = -500.0
      speed_limit = 50.0 !m/s - TODO unharcode
      disp_limit = speed_limit * DT !save some calculations
 
      check_outliers = .FALSE.
      too_fast = .FALSE.

      !First check for particles which have left the domain or are travelling too fast <- TODO 
      !TODO - hook this up to the computations in effload etc
      DO i=1,NN
        IF(Lost(i)) CYCLE

        !Check speed (not necessarily lost)
        dx = UTP(6*I-5) - UT%M(6*I-5)
        dy = UTP(6*I-4) - UT%M(6*I-4)
        dz = UTP(6*I-3) - UT%M(6*I-3)
        IF(ABS(dx) > disp_limit/2.0 .OR. ABS(dy) > disp_limit/2.0 &
             .OR. ABS(dz) > disp_limit/2.0) THEN

          disp = SQRT(dx**2.0 + dy**2.0 + dz**2.0)
          IF(disp > disp_limit) THEN
            IF(DebugMode) PRINT *,myid,' debug, particle breaking the speed limit!'
            too_fast(i) = .TRUE.
          END IF
        END IF

        !Check location
        !TODO - check vertical coordinate too
        x = NRXF%M(1,I) + UTP(6*I-5)
        y = NRXF%M(2,I) + UTP(6*I-4)
        IF(x < grid_bbox(1) .OR. x > grid_bbox(2) .OR. y < grid_bbox(3) .OR. y > grid_bbox(4)) THEN
          Lost(i) = .TRUE.
        END IF

        !Check total displacement
        IF(ABS(UTP(6*I-5)) > MAXUT .OR. ABS(UTP(6*I-4)) > MAXUT .OR. &
             ABS(UTP(6*I-3)) > MAXUT) Lost(i) = .TRUE.

        IF(Lost(i)) THEN
          PRINT *, myid, " Lost a particle : ",I

          !Set as outlier
          Outlier(i) = .TRUE.

          !Put it in the freezer
          UTP(6*I-5) = freezer(1) - NRXF%M(1,I)
          UTP(6*I-4) = freezer(2) - NRXF%M(2,I)
          UTP(6*I-3) = freezer(3) - NRXF%M(3,I)
          UTP(6*I-2) = UT%M(6*I-2)
          UTP(6*I-1) = UT%M(6*I-1)
          UTP(6*I-0) = UT%M(6*I-0)
        END IF
      END DO


      minx = HUGE(minx)
      miny = HUGE(miny)
      minz = HUGE(minz)
      maxx = -HUGE(maxx)
      maxy = -HUGE(maxy)
      maxz = -HUGE(maxz)
      midx = 0.0
      midy = 0.0
      midz = 0.0

      !Find min, max & mean coords
      DO i=1,NN
        IF(Outlier(i)) CYCLE !either because Lost above, or previous iteration outlier
        X=NRXF%M(1,i)+UT%M(6*i-5)
        Y=NRXF%M(2,i)+UT%M(6*i-4)
        Z=NRXF%M(3,i)+UT%M(6*i-3)

        pos(1,i) = X
        pos(2,i) = Y
        pos(3,i) = Z
        pos_sort(:,i) = i

        minx = MIN(minx, x)
        miny = MIN(miny, y)
        minz = MIN(minz, z)
        maxx = MAX(maxx, x)
        maxy = MAX(maxy, y)
        maxz = MAX(maxz, z)
        midx = midx + X
        midy = midy + Y
        midz = midz + Z
      END DO

      midx = midx / NN
      midy = midy / NN
      midz = midz / NN

      BBox(1) = minx
      BBox(2) = maxx
      BBox(3) = miny
      BBox(4) = maxy
      BBox(5) = minz
      BBox(6) = maxz

      !Communicate initial BBox volume - any partition with a much larger
      !than median bbox volume will check for outlier particles
      bbox_vol = (maxx - minx) * (maxy - miny) * (maxz - minz)
      CALL MPI_AllGather(bbox_vol,1,MPI_DOUBLE_PRECISION,bbox_vols,1,&
           MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)

      CALL sort_real(bbox_vols,ntasks)
      med_vol = bbox_vols(ntasks/2)

      IF(DebugMode) THEN
        PRINT *,myid,' debug bbox vol: ',bbox_vol
        IF(myid==0) PRINT *,' median bbox vol: ',med_vol
      END IF

      IF(bbox_vol / med_vol > 2.0) THEN
        IF(DebugMode) PRINT *,myid,' might have outliers!'
        check_outliers = .TRUE.
      END IF

      IF(check_outliers) THEN

        !Compute the distance from mean in every dimension
        dist = 0.0
        dist(1,:) = pos(1,:) - midx
        dist(2,:) = pos(2,:) - midy
        dist(3,:) = pos(3,:) - midz

        mean_dist(1) = SUM(ABS(dist(1,:)))/NN
        mean_dist(2) = SUM(ABS(dist(2,:)))/NN
        mean_dist(3) = SUM(ABS(dist(3,:)))/NN

        !Compute std devs from mean position in each dimension
        std_dev(1) = SQRT(SUM(dist(1,:)**2.0)/NN)
        std_dev(2) = SQRT(SUM(dist(2,:)**2.0)/NN)
        std_dev(3) = SQRT(SUM(dist(3,:)**2.0)/NN)

        sigma_dist(1,:) = dist(1,:) / std_dev(1)
        sigma_dist(2,:) = dist(2,:) / std_dev(2)
        sigma_dist(3,:) = dist(3,:) / std_dev(3)

        poutlier = .FALSE.

        DO i=1,NN
          poutlier(i) = ANY(ABS(sigma_dist(:,i)) > 4.0)
        END DO

        ! !Compute the distance from mean in every dimension
        ! DO i=1,NN
        !   dist(1,i) = NRXF%M(1,i)+UT%M(6*i-5) - midx
        !   dist(2,i) = NRXF%M(2,i)+UT%M(6*i-4) - midy
        !   dist(3,i) = NRXF%M(3,i)+UT%M(6*i-3) - midz
        !   dist_sort(:,i) = i
        ! END DO

        ! !Sort the above
        ! CALL sort_real2(dist(1,:),dist_sort(1,:),NN)
        ! CALL sort_real2(dist(2,:),dist_sort(2,:),NN)
        ! CALL sort_real2(dist(3,:),dist_sort(3,:),NN)



        ! !Look for outliers
        ! max_gap(1) = (maxx - minx) * outlier_gap_prop
        ! max_gap(2) = (maxy - miny) * outlier_gap_prop
        ! max_gap(3) = (maxz - minz) * outlier_gap_prop

        ! idx_ulimit = NINT(NN * (1-outlier_arr_prop))
        ! idx_llimit = NINT(NN * outlier_arr_prop)

        ! outlier = .FALSE.
        ! DO i=1,NN-1
        !   gap(1) = dist(1,i+1) - dist(1,i)
        !   gap(2) = dist(2,i+1) - dist(2,i)
        !   gap(3) = dist(3,i+1) - dist(3,i)

        !   DO j=1,3
        !     IF(gap(j) > max_gap(j)) THEN
        !       IF(i < idx_llimit) THEN
        !         outlier(j,1:i) = .TRUE.
        !       ELSE IF(i > idx_ulimit) THEN
        !         outlier(j,i+1:NN) = .TRUE.
        !       END IF
        !     END IF
        !   END DO
        ! END DO

        ! DO i=1,NN
        !   DO j=1,3
        !     IF(outlier(j,i)) poutlier(dist_sort(j,i)) = .TRUE.
        !   END DO
        ! END DO



        IF(ANY(poutlier)) THEN
          IF(DebugMode) PRINT *,myid, ' has ',COUNT(poutlier),' STDEV outliers: '
          DO i=1,NN
            IF(poutlier(i)) THEN
              IF(DebugMode) PRINT *,myid, ' outlier: ',i,' stdev: ',sigma_dist(:,i),&
                   ' coords: ',NRXF%M(1,i)+UT%M(6*i-5),&
                   NRXF%M(2,i)+UT%M(6*i-4),NRXF%M(3,i)+UT%M(6*i-3)

              DO j=1,NTOT
                IF(ANY(NANS(:,j) == i)) THEN
                  IF(DebugMode) PRINT *,myid,' outlier ',i,j,' efs: ',EFS(j)
                  IF(EFS(j) > 0.0) poutlier(i) = .FALSE.
                END IF
              END DO
            END IF
          END DO
        ENDIF

        IF(DebugMode) THEN
          IF(ANY(poutlier)) THEN
            PRINT *,myid, ' has ',COUNT(poutlier),' actual outliers: '
            DO i=1,NN
              IF(poutlier(i)) THEN
                PRINT *,myid, ' actual outlier: ',i,' stdev: ',sigma_dist(:,i),&
                     ' coords: ',NRXF%M(1,i)+UT%M(6*i-5),&
                     NRXF%M(2,i)+UT%M(6*i-4),NRXF%M(3,i)+UT%M(6*i-3)
              END IF
            END DO
            PRINT *,myid,' old bbox: ',minx,maxx,miny,maxy,minz,maxz
          ENDIF
        END IF

        !Mark final outliers
        DO i=1,NN
          IF(poutlier(i)) Outlier(i) = .TRUE.
        END DO

      END IF !check_outliers

      !halt any particles breaking the speed limit
      DO i=1,NN
        IF(too_fast(i)) THEN
          UTP(6*I-5) = UT%M(6*I-5)
          UTP(6*I-4) = UT%M(6*I-4)
          UTP(6*I-3) = UT%M(6*I-3)
        END IF
      END DO

    END SUBROUTINE CheckSolution

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

    !A subroutine to quicksort a real array (arr1)
    SUBROUTINE sort_real(arr1,n)
      IMPLICIT NONE
      INTEGER :: n
      REAL(KIND=dp) :: arr1(:)
      IF(n <= 1) RETURN
      CALL sort_real_r(arr1,1,n)
    END SUBROUTINE sort_real

    !Does the actual sorting - needs a wrapper to avoid
    !having to pass the extents first
    RECURSIVE SUBROUTINE sort_real_r(arr1,start,fin)
      IMPLICIT NONE
      REAL(KIND=dp) :: arr1(:),hold,pivot_val
      INTEGER :: start,fin,pivot,i,j,hold_int
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

        i=i+1
        j=j-1
      END DO

      IF(j+1 < fin) THEN
        CALL sort_real_r(arr1,j+1,fin)
      END IF
      IF(start < i-1) THEN
        CALL sort_real_r(arr1,start, i-1)
      END IF

      RETURN
    END SUBROUTINE sort_real_r

    !A subroutine to quicksort a real array (arr1)
    ! while also sorting an index arrayarr2 by arr1
    SUBROUTINE sort_real2(arr1,arr2,n)
      IMPLICIT NONE
      INTEGER :: arr2(:),n
      REAL(KIND=dp) :: arr1(:)
      IF(n <= 1) RETURN
      CALL sort_real2_r(arr1,arr2,1,n)
    END SUBROUTINE sort_real2

    !Does the actual sorting - needs a wrapper to avoid
    !having to pass the extents first
    RECURSIVE SUBROUTINE sort_real2_r(arr1,arr2,start,fin)
      IMPLICIT NONE
      INTEGER :: arr2(:)
      REAL(KIND=dp) :: arr1(:),hold,pivot_val
      INTEGER :: start,fin,pivot,i,j,hold_int
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

        hold_int = arr2(i)
        arr2(i) = arr2(j)
        arr2(j) = hold_int

        i=i+1
        j=j-1
      END DO

      IF(j+1 < fin) THEN
        CALL sort_real2_r(arr1,arr2,j+1,fin)
      END IF
      IF(start < i-1) THEN
        CALL sort_real2_r(arr1,arr2,start, i-1)
      END IF

      RETURN
    END SUBROUTINE sort_real2_r

    FUNCTION InterpRast(x,y,RAST,grid,origin,miss_strategy,fill_val) RESULT(zval)
      REAL(KIND=dp) :: x,y,grid,origin(2)
      REAL(KIND=dp) :: RAST(0:,0:)
      INTEGER :: miss_strategy
      REAL(KIND=dp), OPTIONAL :: fill_val
      !--------------------
      REAL(KIND=dp) :: I1,I2,zval
      INTEGER :: XK,YK

      XK = FLOOR((x - origin(1))/grid)
      YK = FLOOR((y - origin(2))/grid)
      I1=(x-origin(1))/grid - XK
      I2=(y-origin(2))/grid - YK

      IF(ValidRasterIndex(xk,yk,RAST)) THEN
        CALL BIPINT(I1,I2,RAST(XK,YK),RAST(XK,YK+1),RAST(XK+1,YK),RAST(XK+1,YK+1),zval)
      ELSE
        SELECT CASE(miss_strategy)
        CASE(INTERP_MISS_FILL)

          IF(.NOT. PRESENT(fill_val)) CALL FatalError("Programming error: &
               &requested fill missing interp but no fill value specified!")
          zval = fill_val
          RETURN

        CASE(INTERP_MISS_NEAREST)

          IF(XK<0) THEN
            XK=0
          ELSE IF(XK>=UBOUND(RAST,1)) THEN
            XK = UBOUND(RAST,1)
          ELSE
            CALL FatalError("Programming Error: Missing interp doesn't make sense")
          END IF
          IF(YK<0) THEN
            YK=0
          ELSE IF(YK>=UBOUND(RAST,2)) THEN
            YK = UBOUND(RAST,2)
          ELSE
            CALL FatalError("Programming Error: Missing interp doesn't make sense")
          END IF

          zval = RAST(XK,YK)
          RETURN

        CASE (INTERP_MISS_ERR)
          PRINT *,myid,' Debug, missing value xy: ',x,y
          CALL FatalError("Unexpected missing value in raster interpolation.")
        CASE DEFAULT
          CALL FatalError("Programming Error: Missing Interp Strategy not understood.")
        END SELECT
      END IF
    END FUNCTION InterpRast

    !-----------------------------------------------------
    Subroutine BIPINT(x,y,f11,f12,f21,f22,fint)
      Implicit none
      Real(KIND=dp) :: x,y,f11,f12,f21,f22,fint
      fint=f11*(1.0-x)*(1.0-y)+f21*x*(1.0-y)+f12*(1.0-x)*y+f22*x*y
    End Subroutine BIPINT

    !-----------------------------------------------------
    Subroutine BIPINTN(x,y,f11,f12,f21,f22,dix,diy,diz,grid)
      Implicit none
      Real(KIND=dp) :: x,y,f11,f12,f21,f22,dix,diy,diz,grid
      REAL(KIND=dp) norm
      dix=(-f11*(1.0-y)+f21*(1.0-y)-f12*y+f22*y)/grid
      diy=(-f11*(1.0-x)-f21*x+f12*(1.0-x)+f22*x)/grid
      diz=1.0
      norm=SQRT(dix**2.0+diy**2.0+1.0)
      dix=-dix/norm
      diy=-diy/norm
      diz=diz/norm
      !if (dix.gt.0.2) dix=0.2
      !if (diy.gt.0.2) diy=0.2
      !if (dix.lt.-0.2) dix=-0.2
      !if (diy.lt.-0.2) diy=-0.2
    End Subroutine BIPINTN


    FUNCTION ValidRasterIndex(i,j,rast) RESULT(valid)
      INTEGER :: i,j
      REAL(KIND=dp) :: rast(0:,0:)
      LOGICAL :: valid

      !'<' rather than '<=' because we also require that i+1 and j+1 are valid
      valid = (i >= 0) .AND. (i < UBOUND(rast,1)) .AND. (j >= 0) .AND. (j < UBOUND(rast,2))
    END FUNCTION ValidRasterIndex
END MODULE UTILS
