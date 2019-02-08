MODULE TypeDefs

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER :: INTERP_MISS_FILL=0,INTERP_MISS_NEAREST=1,INTERP_MISS_ERR=2

  REAL(KIND=dp) :: part_expand=0.5
  INTEGER :: myid, ntasks
  LOGICAL :: DebugMode,PrintTimes

  INCLUDE 'param.dat'
  INCLUDE 'mpif.h'
  INCLUDE 'na90.dat'

  !Structure for holding initial position of this partition's points, as well
  !as relevant neighbouring [C]onnected and [P]roximal points.
  !Note that %C & %P are not currently in use.
  !NB: The %NP member, in particular, represents the *size* 
  ! of the array (from first to last active member)
  ! which may contain proximal points, but NOT the total number of proximal 
  ! points at any given time (because some may be invalid) - check for PartInfo == -1
  TYPE NRXF_t
     REAL(KIND=dp), ALLOCATABLE :: A(:,:)
     REAL(KIND=dp), POINTER :: M(:,:)=>NULL(), C(:,:)=>NULL(), P(:,:)=>NULL()
     INTEGER :: mstrt, cstrt, pstrt,NN,NC,NP
     INTEGER, ALLOCATABLE :: PartInfo(:,:), GID(:)
  END TYPE NRXF_t

  !Type to hold information on particles shared between partitions
  !Our partition receives 'CCount' connected and 'PCount' proximal particles
  !from partition 'NID'. The other partition's particle numbers for these are
  !stored in ConnIDs and ProxIDs respectively, and the corresponding array values
  !in ConnLocs and ProxLocs give the array addresses of NRXF (and UT, UTM) to which
  !these shared particles correspond. NID is in fact redundant because HiDEM allocates
  !an array InvPartInfo(0:ntasks-1), where InvPartInfo(i) % NID = i.
  !Thus, if InvPartInfo(3) % ConnIDs(1) = 13, and % ConnLocs(1) = 2459, this means:
  !Partition 3's 13th particle is shared with this partition, and its initial position
  !and displacement is stored in NRXF(:,2459) and UT(2459) respectively.
  !
  !Inverse information (connected and proximal particles *sent* to partition i are
  !stored in SConnIDs(:), SProxIDs(:), SCCount, SPCount
  TYPE InvPartInfo_t
     INTEGER, ALLOCATABLE :: ConnIDs(:), ConnLocs(:), ProxIDs(:), ProxLocs(:)
     INTEGER, ALLOCATABLE :: SConnIDs(:), SProxIDs(:)
     INTEGER :: CCount, PCount, SCCount, SPCount, NID
  END TYPE InvPartInfo_t

  TYPE UT_t
     REAL(KIND=dp), ALLOCATABLE :: A(:)
     REAL(KIND=dp), POINTER :: M(:)=>NULL(), C(:)=>NULL(), P(:)=>NULL()
  END TYPE UT_t

  !Loading melange from a previous sim requires one partition reading *all* partitions info
  ! for flexiblity. Therefore need an allocatable datastructure to allow holding all that info
  TYPE MelangeDataHolder_t
     REAL(KIND=dp), ALLOCATABLE :: EFS(:)
     REAL(KIND=dp) :: BBox(6)
     INTEGER, ALLOCATABLE :: NANS(:,:),NANPart(:),NDL(:,:)
     TYPE(InvPartInfo_t), ALLOCATABLE :: InvPartInfo(:)
     TYPE(NRXF_t) :: NRXF
     TYPE(UT_t) ::  UT, UTM
     INTEGER :: NTOT,NN,ND
     LOGICAL, ALLOCATABLE :: IsOutlier(:), IsLost(:),PartIsNeighbour(:)
     LOGICAL :: Active=.FALSE.
  END TYPE MelangeDataHolder_t


  TYPE PointEx_t
     INTEGER :: partid=-1,scount=0,rcount=0
     INTEGER, ALLOCATABLE :: SendIDs(:), RecvIDs(:),SendGIDs(:), RecvGIDs(:)
     REAL(KIND=dp), ALLOCATABLE :: S(:),R(:)
  END TYPE PointEx_t

  !Memory-flexible alternative for CN/NCN
  TYPE Conn_t
     INTEGER, ALLOCATABLE :: Conn(:), Part(:)
     INTEGER :: NCN, ID
     REAL(KIND=dp), ALLOCATABLE :: Dists(:)
  END TYPE Conn_t

END MODULE TypeDefs
