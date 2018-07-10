MODULE UTILS

  USE TypeDefs

  IMPLICIT NONE

  INTERFACE PointDataInit
     MODULE PROCEDURE PointDataInitNRXF, PointDataInitUT
  END INTERFACE PointDataInit

  INTERFACE ResizePointData
     MODULE PROCEDURE ResizePointDataNRXF, ResizePointDataUT
  END INTERFACE ResizePointData

  CONTAINS

    SUBROUTINE ExpandIntArray(intarr)

      INTEGER, ALLOCATABLE :: intarr(:), workarr(:)
      INTEGER :: newsize, oldsize

      oldsize = SIZE(intarr)
      newsize =  oldsize * 2

      ! ALLOCATE(workarr(oldsize))
      ! workarr(1:oldsize) = intarr(1:oldsize)
      ! DEALLOCATE(intarr)
      ALLOCATE(workarr(newsize))
      workarr(1:oldsize) = intarr(1:oldsize)

      PRINT *,'DEBUG: Doing move alloc'

      !FORTRAN 2003 feature...
      CALL MOVE_ALLOC(workarr, intarr)

      PRINT *,'DEBUG: Done move alloc'

    END SUBROUTINE ExpandIntArray

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
      
      ALLOCATE(NRXF%A(3,n_tot),NRXF%PartInfo(2,n_tot))

      NRXF%PartInfo(:,:) = -1 !-1 = no point

      NRXF%A = 0.0
      NRXF%M => NRXF%A(:,1:n)
      NRXF%P => NRXF%A(:,n+1:UBOUND(NRXF%A,2))

      PRINT *,'Debug,nrxf n, ntot: ',n,n_tot, SIZE(NRXF%A), SIZE(NRXF%M), SIZE(NRXF%P), UBOUND(NRXF%A,2)

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
    END SUBROUTINE PointDataInitUT

    SUBROUTINE ResizePointDataNRXF(NRXF,scale,do_M,do_P)

      TYPE(NRXF2_T), TARGET :: NRXF
      REAL*8 :: scale
      LOGICAL, OPTIONAL :: do_M, do_P
      !---------------------
      REAL*8, ALLOCATABLE :: work_arr(:,:)
      INTEGER :: a_oldsize,m_oldsize,p_oldsize
      INTEGER :: a_newsize,m_newsize,p_newsize
      LOGICAL :: doM,doP
      
      IF(PRESENT(do_M)) THEN
        doM = do_M
      ELSE
        doM = .FALSE.
      END IF

      IF(PRESENT(do_P)) THEN
        doP = do_P
      ELSE
        doP = .TRUE.
      END IF

      IF(scale < 1.0) THEN
        PRINT *,'NRXF size reduction not yet implemented, sorry!'
        STOP
      END IF

      a_oldsize = SIZE(NRXF%A,2)
      m_oldsize = SIZE(NRXF%M,2)
      p_oldsize = SIZE(NRXF%P,2)

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

      PRINT *,myid,' debug resize nrxf: ',a_oldsize, m_oldsize, p_oldsize, a_newsize, m_newsize, p_newsize

      ALLOCATE(work_arr(3,a_newsize))
      work_arr = 0.0
      work_arr(:,1:m_oldsize) = NRXF % M(:,1:m_oldsize)
      work_arr(:,m_newsize+1 : m_newsize+p_oldsize) = NRXF % P(:,1:p_oldsize)

      CALL MOVE_ALLOC(work_arr, NRXF%A)

      NRXF%M => NRXF%A(:,1:m_newsize)
      NRXF%P => NRXF%A(:,m_newsize+1 : a_newsize)

      PRINT *,myid,' debug2 resize nrxf: ',SIZE(NRXF%A,2),SIZE(NRXF%M,2),SIZE(NRXF%P,2)

    END SUBROUTINE ResizePointDataNRXF

    SUBROUTINE ResizePointDataUT(UT,scale,do_M,do_P)

      TYPE(UT2_T), TARGET :: UT
      REAL*8 :: scale
      LOGICAL, OPTIONAL :: do_M, do_P
      !---------------------
      REAL*8, ALLOCATABLE :: work_arr(:)
      INTEGER :: a_oldsize,m_oldsize,p_oldsize
      INTEGER :: a_newsize,m_newsize,p_newsize
      LOGICAL :: doM,doP
      
      IF(PRESENT(do_M)) THEN
        doM = do_M
      ELSE
        doM = .FALSE.
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
