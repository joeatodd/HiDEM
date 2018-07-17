MODULE UTILS

  USE TypeDefs

  IMPLICIT NONE

  INTERFACE PointDataInit
     MODULE PROCEDURE PointDataInitNRXF, PointDataInitUT
  END INTERFACE PointDataInit

  INTERFACE ResizePointData
     MODULE PROCEDURE ResizePointDataNRXF, ResizePointDataUT
  END INTERFACE ResizePointData

  INTERFACE ExpandIntarray
     MODULE PROCEDURE Expand1IntArray, Expand2IntArray
  END INTERFACE ExpandIntarray

  CONTAINS

    SUBROUTINE Expand1IntArray(intarr)

      INTEGER, ALLOCATABLE :: intarr(:), workarr(:)
      INTEGER :: newsize, oldsize

      oldsize = SIZE(intarr)
      newsize =  oldsize * 2

      ! ALLOCATE(workarr(oldsize))
      ! workarr(1:oldsize) = intarr(1:oldsize)
      ! DEALLOCATE(intarr)
      ALLOCATE(workarr(newsize))
      workarr(1:oldsize) = intarr(1:oldsize)

      !FORTRAN 2003 feature...
      CALL MOVE_ALLOC(workarr, intarr)

      IF(DebugMode) PRINT *,'DEBUG: Done move alloc'

    END SUBROUTINE Expand1IntArray

    SUBROUTINE Expand2IntArray(intarr)

      INTEGER, ALLOCATABLE :: intarr(:,:), workarr(:,:)
      INTEGER :: newsize, oldsize, dim1size

      oldsize = SIZE(intarr,2)
      newsize =  oldsize * 2
      dim1size = SIZE(intarr,1)

      ! ALLOCATE(workarr(oldsize))
      ! workarr(1:oldsize) = intarr(1:oldsize)
      ! DEALLOCATE(intarr)
      ALLOCATE(workarr(dim1size,newsize))
      workarr(:,1:oldsize) = intarr(:,1:oldsize)

      !FORTRAN 2003 feature...
      CALL MOVE_ALLOC(workarr, intarr)

      IF(DebugMode) PRINT *,'DEBUG: Done move alloc'

    END SUBROUTINE Expand2IntArray

    SUBROUTINE PointDataInitNRXF(NRXF, n, partexpand)
      TYPE(NRXF2_T), TARGET :: NRXF
      INTEGER :: n,n_tot
      REAL*8 :: partexpand

      n_tot = n + CEILING(n * partexpand)

      IF(ALLOCATED(NRXF%A)) THEN
        PRINT *, "Programming error: abuse of PointDataInitNRXF"
        STOP
      END IF

      NRXF%mstrt = 1
      NRXF%cstrt = 1 + n
      NRXF%mstrt = 1 + n !no connected particles, initially

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

    SUBROUTINE PointDataInitUT(UT, n, partexpand)
      TYPE(UT2_T), TARGET :: UT
      INTEGER :: n,n_tot
      REAL*8 :: partexpand

      n_tot = n + CEILING(n * partexpand)

      IF(ALLOCATED(UT%A)) THEN
        PRINT *, "Programming error: abuse of PointDataInitUT"
        STOP
      END IF

      ALLOCATE(UT%A(6*n_tot))
      UT%A = 0.0
      UT%M => UT%A(1:6*n)
      UT%P => UT%A((6*n+1):UBOUND(UT%A,1))

      IF(DebugMode) PRINT *,'Debug,ut init n, ntot: ',n,n_tot, SIZE(UT%A), SIZE(UT%M), SIZE(UT%P)

    END SUBROUTINE PointDataInitUT

    SUBROUTINE InvPartInfoInit(InvPartInfo,neighparts,initsize_in)
      TYPE(InvPartInfo_t),ALLOCATABLE :: InvPartInfo(:)
      INTEGER :: neighparts(:)
      INTEGER, OPTIONAL :: initsize_in
      !------------------------------
      INTEGER :: initsize,i,j,neighcount

      neighcount = COUNT(neighparts /= -1)

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

      DO i=1,neighcount
        j = neighparts(i)
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

    SUBROUTINE ResizePointDataNRXF(NRXF,scale,do_M,do_C,do_P)

      TYPE(NRXF2_T), TARGET :: NRXF
      REAL*8 :: scale
      LOGICAL, OPTIONAL :: do_M, do_C, do_P
      !---------------------
      REAL*8, ALLOCATABLE :: work_arr(:,:)
      INTEGER, ALLOCATABLE :: work_int(:,:)
      INTEGER :: a_oldsize,m_oldsize,c_oldsize,p_oldsize
      INTEGER :: a_newsize,m_newsize,c_newsize,p_newsize
      INTEGER :: cstrt_new, pstrt_new
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

      a_oldsize = SIZE(NRXF%A,2)
      m_oldsize = NRXF%cstrt - NRXF%mstrt

      IF(SIZE(NRXF%M,2) /= m_oldsize) THEN
        PRINT *, "ERROR: NRXF%M wrong size in ResizePointData"
        STOP
      END IF

      c_oldsize = NRXF%pstrt - NRXF%cstrt
      p_oldsize = a_oldsize - NRXF%cstrt + 1

      IF(doM) THEN
        m_newsize = CEILING(m_oldsize*scale)
      ELSE
        m_newsize = m_oldsize
      END IF

      IF(doC) THEN
        c_newsize = CEILING(c_oldsize*scale)
      ELSE
        c_newsize = c_oldsize
      END IF

      IF(doP) THEN
        p_newsize = CEILING(p_oldsize*scale)
      ELSE
        p_newsize = p_oldsize
      END IF

      a_newsize = m_newsize + c_newsize + p_newsize
      
      cstrt_new = m_newsize + 1
      pstrt_new = m_newsize + c_newsize + 1

      IF(DebugMode) PRINT *,myid,' debug resize nrxf: ',a_oldsize, m_oldsize, p_oldsize,&
           'new: ',a_newsize, m_newsize, p_newsize, 'strts: ', cstrt_new, pstrt_new

      ALLOCATE(work_arr(3,a_newsize),work_int(2,a_newsize))
      work_arr = 0.0
      work_int = 0

      work_arr(:,1:m_oldsize) = NRXF%A(:,1:m_oldsize)
      work_arr(:,cstrt_new : cstrt_new+c_oldsize) = NRXF % A(:,NRXF%cstrt:NRXF%cstrt+c_oldsize)
      work_arr(:,pstrt_new : pstrt_new+p_oldsize) = NRXF % A(:,NRXF%pstrt:NRXF%pstrt+p_oldsize)

      work_int(:,1:m_oldsize) = NRXF%PartInfo(:,1:m_oldsize)
      work_int(:,cstrt_new : cstrt_new+c_oldsize) = NRXF%PartInfo(:,NRXF%cstrt:NRXF%cstrt+c_oldsize)
      work_int(:,pstrt_new : pstrt_new+p_oldsize) = NRXF%PartInfo(:,NRXF%pstrt:NRXF%pstrt+p_oldsize)

      CALL MOVE_ALLOC(work_arr, NRXF%A)
      CALL MOVE_ALLOC(work_int, NRXF%PartInfo)

      NRXF%cstrt = cstrt_new
      NRXF%pstrt = pstrt_new
      
      NRXF%M => NRXF%A(:,1:m_newsize)
      ! NRXF%P => NRXF%A(:,m_newsize+1 : a_newsize)


      PRINT *,myid,' debug2 resize nrxf: ',SIZE(NRXF%A,2),NRXF%cstrt, NRXF%pstrt

    END SUBROUTINE ResizePointDataNRXF

    SUBROUTINE ResizePointDataUT(UT,scale,do_M,do_C,do_P)

      TYPE(UT2_T), TARGET :: UT
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
        p_newsize = CEILING(p_oldsize*scale)
      ELSE
        p_newsize = p_oldsize
      END IF

      IF(doM) THEN
        m_newsize = CEILING(m_oldsize*scale)
      ELSE
        m_newsize = m_oldsize
      END IF

      a_newsize = m_newsize + p_newsize

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

END MODULE UTILS
