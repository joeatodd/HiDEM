MODULE Melange

USE INOUT
USE UTILS
USE TypeDefs

IMPLICIT NONE
!INCLUDE 'na90.dat'

CONTAINS

SUBROUTINE LoadMelange(SI, melange_data)
  INTEGER :: i,j,k,mel_parts, mel_nn,ierr,counter,N1,N2,P1
  INTEGER :: cstrt,NC,pstrt,NP,NSIZE,NTOT,BCC,IX,Part,ID
  INTEGER, ALLOCATABLE :: pnn(:),mel_pstrt(:)
  INTEGER, POINTER :: countptr
  REAL(KIND=dp) :: x,y,z,m,X1,E
  CHARACTER(LEN=256) :: MelangeRunName,wrkdir
  LOGICAL :: FileExists
  TYPE(NRXF_t) :: melange_NRXF
  TYPE(UT_t) :: melange_UT, melange_UTM
  TYPE(MelangeDataHolder_t), ALLOCATABLE, TARGET :: mel_data(:)
  TYPE(MelangeDataHolder_t) :: melange_data
  TYPE(SimInfo_t) :: SI

  !Note: Need to load all the files - not just our partitions.

  wrkdir = SI%wrkdir
  MelangeRunName = SI%MelRunName

  !Count the REST0 files
  i = 0
  DO WHILE(.TRUE.)
    INQUIRE( FILE=TRIM(wrkdir)//'/'//TRIM(MelangeRunName)//'_REST0'//na(i), EXIST=FileExists )
    IF(.NOT. FileExists) THEN
      IF(i==0) CALL FatalError("Didn't find partition 0 file for melange restart")
      EXIT
    END IF
    i = i + 1
  END DO
  mel_parts = i-1+1
  PRINT *,'Debug, number of partitions in melange sim: ',mel_parts

  ALLOCATE(pnn(0:mel_parts-1), &
       mel_pstrt(0:mel_parts-1), &
       mel_data(0:mel_parts-1))
  pnn = 0

  DO i=0,mel_parts-1
    ALLOCATE(mel_data(i) % PartIsNeighbour(0:mel_parts))
    mel_data(i) % PartIsNeighbour = .FALSE.
  END DO

  !Read number of particles & initialize data structures (for each partition)
  mel_nn = 0
  DO i=0,mel_parts-1

    OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(MelangeRunName)//'_REST0'//na(i),&
         STATUS='OLD',ACTION='read')
    READ(117+myid,*) pnn(i),cstrt,NC,pstrt,NP,NSIZE,mel_data(i) % NTOT,BCC
    !    READ(117+myid,*) MAXX,MAXY,MAXZ,DMPEN,ENM0
    !    READ(117+myid,*) DPE,BCE,MGH0,GSUM0,PSUM,T,RY0
    CLOSE(117+myid)

    mel_nn = mel_nn + pnn(i)
    mel_data(i) % NN = pnn(i)

    CALL PointDataInit(mel_data(i)%NRXF, pnn(i), arrsize=NSIZE)
    mel_data(i)%NRXF % cstrt = cstrt
    mel_data(i)%NRXF % NC = NC
    mel_data(i)%NRXF % pstrt = pstrt
    CALL PointDataInit(mel_data(i)%UT,mel_data(i)%NRXF)
    CALL PointDataInit(mel_data(i)%UTM,mel_data(i)%NRXF)

  END DO
  PRINT *,'Debug, total number of particles in melange: ',mel_nn

  mel_pstrt(0) = 1
  DO i=1,mel_parts-1
    mel_pstrt(i) = mel_pstrt(i-1) + pnn(i-1)
  END DO

  !Read particle init pos (for each partition)
  DO i=0,mel_parts-1

    OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(MelangeRunName)//'_NODFIL2'//na(i),&
         STATUS='UNKNOWN',ACTION='read',IOSTAT=ierr)
    IF(ierr /= 0) CALL FatalError("IO error reading melange _NODFIL2")

    DO j=1,pnn(i)
      READ(117+myid,*) IX,X,Y,Z,M
      mel_data(i)%NRXF%M(1,j)=X
      mel_data(i)%NRXF%M(2,j)=Y
      mel_data(i)%NRXF%M(3,j)=Z
      mel_data(i)%NRXF%PartInfo(1,IX) = i
      mel_data(i)%NRXF%PartInfo(2,IX) = IX
    END DO
    CLOSE (117+myid)
  END DO

  !Read in the partition & id of particles we share from other partitions (for each partition)
  DO i=0,mel_parts-1

    OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(MelangeRunName)//'_ONODFIL2'//na(i),&
         STATUS='UNKNOWN',ACTION='read')
    DO j=1,mel_data(i)%NRXF % NC!+NP
      READ(117+myid,*) IX,Part,ID
      IF(IX >= mel_data(i)%NRXF%pstrt .OR. IX < mel_data(i)%NRXF%cstrt) THEN
        PRINT *,myid,' debug IX, pstrt, cstrt: ',IX,mel_data(i)%NRXF%pstrt, mel_data(i)%NRXF%cstrt
        CALL FatalError("Restart programming error - NRXF bounds")
      END IF
      mel_data(i)%NRXF%PartInfo(1,IX) = Part
      mel_data(i)%NRXF%PartInfo(2,IX) = ID
      mel_data(i) % PartIsNeighbour(Part) = .TRUE.
    END DO
    CLOSE (117+myid)

  END DO

  !Construct InvPartInfo (for each partition)
  DO i=0,mel_parts-1
    CALL InvPartInfoInit(mel_data(i) % InvPartInfo, mel_data(i) % PartIsNeighbour)

    DO j=mel_data(i)%NRXF%cstrt, mel_data(i)%NRXF%cstrt + mel_data(i)%NRXF%NC - 1
      Part = mel_data(i)%NRXF%PartInfo(1,j)
      ID = mel_data(i)%NRXF%PartInfo(2,j)
      countptr => mel_data(i) % InvPartInfo(Part) % ccount
      countptr = countptr + 1

      IF(countptr > SIZE(mel_data(i) % InvPartInfo(Part) % ConnIDs)) THEN
        CALL ExpandIntArray(mel_data(i) % InvPartInfo(Part) % ConnIDs)
        CALL ExpandIntArray(mel_data(i) % InvPartInfo(Part) % ConnLocs)
      END IF

      mel_data(i) % InvPartInfo(Part) % ConnIDs(countptr) = ID
      mel_data(i) % InvPartInfo(Part) % ConnLocs(countptr) = I
    END DO
  END DO

  !Read in bond info (for each partition)
  DO i=0,mel_parts-1
    ALLOCATE(mel_data(i) % EFS(mel_data(i) % NTOT), &
         mel_data(i) % NANS(2,mel_data(i) % NTOT), &
         mel_data(i) % NANPart(mel_data(i) % NTOT))
    OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(MelangeRunName)//'_FS'//na(i),&
         STATUS='UNKNOWN',ACTION='read')
    DO j=1,mel_data(i) % NTOT
      READ(117+myid,*) N1,N2,P1,X1,X1,X1,X1,X1,X1,E
      mel_data(i) % NANS(1,j)=N1
      mel_data(i) % NANS(2,j)=N2
      mel_data(i) % NANPart(j)=P1
      mel_data(i) % EFS(j)=E
    END DO
    CLOSE (117+myid)
  END DO

  DO i=0,mel_parts-1
    ALLOCATE(mel_data(i) % IsOutlier(mel_data(i) % NN), &
         mel_data(i) % IsLost(mel_data(i) % NN))

    OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(MelangeRunName)//'_REST2'//na(i),&
         STATUS='OLD',ACTION='read')
    DO j=1,pnn(i)
      READ(117+myid,*) mel_data(i)%UT%M(6*J-5),mel_data(i)%UT%M(6*J-4),mel_data(i)%UT%M(6*J-3),&
           mel_data(i)%UT%M(6*J-2),mel_data(i)%UT%M(6*J-1),mel_data(i)%UT%M(6*J-0)
      READ(117+myid,*) mel_data(i)%UTM%M(6*J-5),mel_data(i)%UTM%M(6*J-4),mel_data(i)%UTM%M(6*J-3),&
           mel_data(i)%UTM%M(6*J-2),mel_data(i)%UTM%M(6*J-1),mel_data(i)%UTM%M(6*J-0)
      READ(117+myid,*) mel_data(i) % IsOutlier(j), mel_data(i) % IsLost(j)
    END DO
    CLOSE (117+myid)

  END DO

  !Construct a list of all melange particle positions and connections
  CALL PointDataInit(melange_data%NRXF, mel_nn, arrsize=mel_nn)
  CALL PointDataInit(melange_data%UT,melange_data%NRXF)
  CALL PointDataInit(melange_data%UTM,melange_data%NRXF)
  ALLOCATE(melange_data % IsLost(mel_nn),&
       melange_data % IsOutlier(mel_nn))

  melange_data % NN = mel_nn

  counter = 0
  DO i=0,mel_parts-1
    DO j=1,mel_data(i)%NN
      counter = counter + 1
      melange_data%NRXF%M(:,counter) = mel_data(i)%nrxf % M(:,j)
      melange_data%UT%M(6*counter-5) = mel_data(i)%UT%M(6*j-5)
      melange_data%UT%M(6*counter-4) = mel_data(i)%UT%M(6*j-4)
      melange_data%UT%M(6*counter-3) = mel_data(i)%UT%M(6*j-3)
      melange_data%UT%M(6*counter-2) = mel_data(i)%UT%M(6*j-2)
      melange_data%UT%M(6*counter-1) = mel_data(i)%UT%M(6*j-1)
      melange_data%UT%M(6*counter-0) = mel_data(i)%UT%M(6*j-0)
      melange_data%UTM%M(6*counter-5) = mel_data(i)%UTM%M(6*j-5)
      melange_data%UTM%M(6*counter-4) = mel_data(i)%UTM%M(6*j-4)
      melange_data%UTM%M(6*counter-3) = mel_data(i)%UTM%M(6*j-3)
      melange_data%UTM%M(6*counter-2) = mel_data(i)%UTM%M(6*j-2)
      melange_data%UTM%M(6*counter-1) = mel_data(i)%UTM%M(6*j-1)
      melange_data%UTM%M(6*counter-0) = mel_data(i)%UTM%M(6*j-0)

      melange_data % IsLost(counter) = mel_data(i) % IsLost(j)
      melange_data % IsOutlier(counter) = mel_data(i) % IsOutlier(j)
    END DO
  END DO

  !Resolve global particle connections...
  !Count them:
  NTOT = 0
  DO i=0,mel_parts-1
    DO j=1,mel_data(i) % NTOT
      IF(mel_data(i) % NANPart(j) > i) CYCLE
      NTOT = NTOT + 1
    END DO
  END DO

  melange_data % NTOT = NTOT
  ALLOCATE(melange_data % NANS(2,NTOT),&
       melange_data % EFS(NTOT))

  !Resolve them:
  counter = 0
  DO i=0,mel_parts-1
    DO j=1,mel_data(i) % NTOT

      P1 = mel_data(i) % NANPart(j)
      IF(P1 > i) CYCLE

      counter = counter + 1

      N1 = mel_data(i) % NANS(1,j)
      N2 = mel_data(i) % NANS(2,j)
      melange_data % NANS(1,counter) = mel_pstrt(P1) + mel_data(i)%NRXF%PartInfo(2,N1) - 1 
      melange_data % NANS(2,counter) = mel_pstrt(i) + mel_data(i)%NRXF%PartInfo(2,N2) - 1
      melange_data % EFS(counter) = mel_data(i) % EFS(j)
    END DO
  END DO
  IF(counter /= NTOT) CALL FatalError("LoadMelange: Programming Error - NTOT =/= counter")

END SUBROUTINE LoadMelange

!Get rid of the melange particles which overlap w/ ice particles, or which are lost
SUBROUTINE Prune_Melange(melange_data,xo,ip,SCL)

  USE Octree

  TYPE(MelangeDataHolder_t), TARGET :: melange_data
  INTEGER :: i,ip
  REAL(KIND=dp), ALLOCATABLE :: xo(:,:)
  REAL(KIND=dp) :: SCL
  !-----------------------
  REAL(KIND=dp) :: oct_BBox(2,3),eps,x(3),xk,yk,zk,search_dist
  REAL(KIND=dp), ALLOCATABLE :: work_nrxf(:,:),work_ut(:),work_utm(:),work_efs(:)
  INTEGER, ALLOCATABLE :: ngb_ids(:),node_loc(:),work_nans(:,:)
  INTEGER :: counter,num_ngb,N1,N2,pruned_NN, pruned_NTOT
  LOGICAL, ALLOCATABLE :: nodeRemove(:),NANRemove(:)

  type(point_type), allocatable :: points(:)

  IF(SIZE(xo,2) < ip) CALL FatalError("Prune_Melange: xo size smaller than ip!")

  eps = 1.0
  search_dist = SCL*1.0

  PRINT *,myid,' debug pruning melange - particles: ',melange_data%NN, ' conns: ',&
       melange_data%NTOT,' ip: ',ip

  ALLOCATE(NodeRemove(melange_data%NN),&
       NANRemove(melange_data%NTOT),&
       Points(ip),&
       ngb_ids(100))

  PRINT *,myid,' survived allocation'


  nodeRemove = .FALSE.
  NANRemove = .FALSE.

  DO i=1,ip
    Points(i) % id = i
    Points(i) % x(1) = xo(1,i)
    Points(i) % x(2) = xo(2,i)
    Points(i) % x(3) = xo(3,i)
  END DO

  oct_bbox(1,1) = MINVAL(xo(1,1:ip)) - eps - search_dist
  oct_bbox(1,2) = MINVAL(xo(2,1:ip)) - eps - search_dist
  oct_bbox(1,3) = MINVAL(xo(3,1:ip)) - eps - search_dist
  oct_bbox(2,1) = MAXVAL(xo(1,1:ip)) + eps + search_dist
  oct_bbox(2,2) = MAXVAL(xo(2,1:ip)) + eps + search_dist
  oct_bbox(2,3) = MAXVAL(xo(3,1:ip)) + eps + search_dist

  PRINT *,myid,' initializing octree'
  CALL Octree_init(max_depth=10, max_num_point=6,bbox=oct_bbox)

  PRINT *,myid,' building octree'
  CALL Octree_build(Points)

  PRINT *,myid,' searching octree'

  DO i=1,melange_data%NN
    num_ngb = 0
    ngb_ids = 0

    x(1) = melange_data % NRXF%M(1,i) + melange_data % UT%M(6*i-5)
    x(2) = melange_data % NRXF%M(2,i) + melange_data % UT%M(6*i-4)
    x(3) = melange_data % NRXF%M(3,i) + melange_data % UT%M(6*i-3)

    CALL Octree_search(x, search_dist, num_ngb, ngb_ids)

    IF(num_ngb > 0) nodeRemove(i) = .TRUE.
  END DO

  !Remove the lost particles (IsLost)
  DO i=1,melange_data%NN
    IF(melange_data % IsLost(i)) nodeRemove(i) = .TRUE.
  END DO

  PRINT *,myid,' destroying octree'

  CALL Octree_final()

  PRINT *,'Debug, removing ',COUNT(nodeRemove),' of ',SIZE(nodeRemove),&
       ' melange particles which overlap, leaving: ',SIZE(nodeRemove) - COUNT(nodeRemove)

  DO i=1,melange_data%NTOT
    N1 = melange_data % NANS(1,i)
    N2 = melange_data % NANS(2,i)

    IF(nodeRemove(N1) .OR. nodeRemove(N2)) THEN
      NANRemove(i) = .TRUE.
    END IF
  END DO

  PRINT *,'Debug, removing ',COUNT(nanremove),' of ',SIZE(nanremove),&
       ' melange connections, leaving: ',SIZE(nanremove) - COUNT(nanremove)

  pruned_NN = COUNT(.not. nodeRemove)
  pruned_NTOT = COUNT(.not. NANRemove)

  ALLOCATE(node_loc(melange_data%NN),&
       work_nrxf(3,pruned_NN),&
       work_ut(6*pruned_NN),&
       work_utm(6*pruned_NN),&
       work_nans(2,pruned_NTOT),&
       work_efs(pruned_NTOT))

  work_nrxf = 0.0
  work_ut = 0.0
  work_utm = 0.0
  work_nans = 0
  work_efs = 0.0

  !fill node_loc - the updated array locations of all particles
  !and work_nrxf, ut, utm
  counter = 0
  DO i=1,melange_data%NN
    IF(nodeRemove(i)) THEN
      node_loc(i) = -1
    ELSE
      counter = counter + 1
      node_loc(i) = counter
      work_nrxf(:,counter) = melange_data%NRXF%M(:,i)
      work_ut(counter*6-5:counter*6) = melange_data%UT%M(i*6-5:i*6)
      work_utm(counter*6-5:counter*6) = melange_data%UTM%M(i*6-5:i*6)
    END IF
  END DO
  IF(counter /= pruned_NN) CALL FatalError("Programming Error: counter /= pruned_NN")

  !Create work_nans
  counter = 0
  DO i=1,melange_data%NTOT
    IF(NANRemove(i)) THEN
      CYCLE
    ELSE
      counter = counter + 1
      work_nans(1,counter) = node_loc(melange_data % NANS(1,i))
      work_nans(2,counter) = node_loc(melange_data % NANS(2,i))
      work_efs(counter) = melange_data % EFS(i)
    END IF
  END DO
  IF(counter /= pruned_NTOT) CALL FatalError("Programming Error: counter /= pruned_NTOT")
  IF(ANY(work_nans < 0) .OR. ANY(work_nans > pruned_NN)) &
       CALL FatalError("Programming Error: work_nans filled incorrectly")

  melange_data % NRXF % NN = pruned_NN
  melange_data % NN = pruned_NN
  melange_data % NTOT = pruned_NTOT

  !TODO - what else to update?

  DEALLOCATE(melange_data%NRXF%A)
  CALL MOVE_ALLOC(work_nrxf,melange_data%NRXF%A)
  melange_data%NRXF%M => melange_data%NRXF%A

  DEALLOCATE(melange_data%UT%A)
  CALL MOVE_ALLOC(work_UT,melange_data%UT%A)
  melange_data%UT%M => melange_data%UT%A

  DEALLOCATE(melange_data%UTM%A)
  CALL MOVE_ALLOC(work_UTM,melange_data%UTM%A)
  melange_data%UTM%M => melange_data%UTM%A

  DEALLOCATE(melange_data%nans)
  CALL MOVE_ALLOC(work_nans, melange_data%nans)

  DEALLOCATE(melange_data%EFS)
  CALL MOVE_ALLOC(work_efs, melange_data%EFS)

  !Set partinfo so that FindNearbyParticles works
  IF(ALLOCATED(melange_data%NRXF%partinfo)) DEALLOCATE(melange_data%NRXF%partinfo)
  ALLOCATE(melange_data%NRXF%partinfo(2,pruned_NN))
  melange_data%NRXF%partinfo(1,:) = 0
  DO i=1,pruned_NN
    melange_data%NRXF%partinfo(2,:) = i
  END DO

  !Compute the bounding box
  melange_data%BBox(1:5:2) = HUGE(melange_data%BBox(1))
  melange_data%BBox(2:6:2) = -HUGE(melange_data%BBox(1))
  DO i=1,pruned_NN
    xk = melange_data%NRXF%M(1,i) + melange_data%UT%M(i*6 - 5)
    yk = melange_data%NRXF%M(2,i) + melange_data%UT%M(i*6 - 4)
    zk = melange_data%NRXF%M(3,i) + melange_data%UT%M(i*6 - 3)

    melange_data%BBox(1) = MIN(xk,melange_data%BBox(1))
    melange_data%BBox(2) = MAX(xk,melange_data%BBox(2))
    melange_data%BBox(3) = MIN(yk,melange_data%BBox(3))
    melange_data%BBox(4) = MAX(yk,melange_data%BBox(4))
    melange_data%BBox(5) = MIN(zk,melange_data%BBox(5))
    melange_data%BBox(6) = MAX(zk,melange_data%BBox(6))
  END DO
  
END SUBROUTINE Prune_Melange

!Convert the 'NAN' info from previous simulation into 'NCN' and 'CN'
SUBROUTINE MelangeBonds(melange_data, glac_ip, CN, nbeams_mel)
  TYPE(Conn_t), ALLOCATABLE :: CN(:)
  TYPE(MelangeDataHolder_t) :: melange_data
  INTEGER :: nbeams_mel, glac_ip
  !--------------------------
  INTEGER :: i,NN,N1,N2

  NN = melange_data % NN
  nbeams_mel = 0

  ALLOCATE(CN(NN))
  DO i=1,NN
    ALLOCATE(CN(i) % Conn(12))
    CN(i) % NCN = 0
    CN(i) % Conn = 0
    CN(i) % ID = i
  END DO

  !Now generate NCN, CN
  DO i=1,melange_data % NTOT
    N1 = melange_data % NANS(1,i)
    N2 = melange_data % NANS(2,i)

    CN(N1) % NCN = CN(N1) % NCN + 1
    CN(N2) % NCN = CN(N2) % NCN + 1

    CN(N1) % Conn(CN(N1) % NCN) = N2 + glac_ip
    CN(N2) % Conn(CN(N2) % NCN) = N1 + glac_ip
  END DO

  nbeams_mel = melange_data % NTOT
END SUBROUTINE MelangeBonds

!Passes UT,UTM & EFS for melange particles from a previous sim.
!When melange is loaded, NRXF + UT (1:3) => NRXF
!We can then revert this back to NRXF + UT (and UTM if desired)
!NOTE: Cray ftn compiler may not like MPI self-send (0 to 0) in this subroutine
SUBROUTINE PassMelangeData(SI,melange_data,NRXF,UT,UTM,NANS,EFS,InvPartInfo)

  TYPE(MelangeDataHolder_t) :: melange_data
  TYPE(SimInfo_t) :: SI
  INTEGER :: NANS(:,:)
  REAL(KIND=dp), ALLOCATABLE :: EFS(:)
  TYPE(UT_t) :: UT,UTM
  TYPE(NRXF_t) :: NRXF
  TYPE(InvPartInfo_t) :: InvPartInfo(0:)
  !----------------------------
  REAL(KIND=dp), ALLOCATABLE :: recv_stream(:)
  INTEGER :: i,j,part,ppcnt,my_ppcnt, my_beamcnt, ierr,counter,seek,N1,N2
  INTEGER :: NTOT, stat(MPI_STATUS_SIZE), statcnt
  INTEGER, ALLOCATABLE :: stats(:),recv_GIDs(:)
  TYPE(PointEx_t), ALLOCATABLE :: PointEx(:)
  LOGICAL :: send, matched

  IF(.NOT. melange_data % active) RETURN

  NTOT = SIZE(NANS, 2)

  IF(myid==0) ALLOCATE(PointEx(0:ntasks-1))
  ALLOCATE(stats(ntasks*3))
  stats = MPI_REQUEST_NULL
  statcnt = 0

  !melange_data % owner(:) contains the partition ID of each melange particle
  !melange_data % GID(:) contains the GID of each melange particle

  statcnt = statcnt + 1
  CALL MPI_IRecv(my_ppcnt, 1, MPI_INTEGER, 0, 300, MPI_COMM_WORLD,stats(statcnt),ierr)

  !Determine & communicate particle counts to partitions
  IF(myid == 0) THEN
    DO part=0,ntasks-1
      ppcnt = COUNT(melange_data % owner == part)
      CALL MPI_Send(ppcnt, 1, MPI_INTEGER, part, 300, MPI_COMM_WORLD, ierr)
    END DO
  END IF

  !Wait for the previous non-blocking sends, then reset stats
  CALL MPI_Waitall(statcnt, stats(1:statcnt), MPI_STATUSES_IGNORE, ierr)
  stats = MPI_REQUEST_NULL
  statcnt = 0

  !Allocate space and set up non-blocking receive
  !Communicate: GID,UT*6,UTM*6
  IF(my_ppcnt > 0) THEN
    ALLOCATE(recv_stream(my_ppcnt * 12),recv_GIDs(my_ppcnt))
    statcnt = statcnt + 1
    CALL MPI_IRecv(recv_GIDs, my_ppcnt, MPI_INTEGER, 0, 301, MPI_COMM_WORLD,stats(statcnt),ierr)
    statcnt = statcnt + 1
    CALL MPI_IRecv(recv_stream, my_ppcnt*12, MPI_DOUBLE_PRECISION, 0, 302, MPI_COMM_WORLD,stats(statcnt),ierr)
  END IF

  !Root constructs & sends arrays
  IF(myid == 0) THEN
    DO part=0,ntasks-1

      ppcnt = COUNT(melange_data % owner == part)
      IF(ppcnt == 0) CYCLE
      PointEx(part) % scount = ppcnt
      ALLOCATE(PointEx(part) % SendGIDs(ppcnt),PointEx(part) % S(ppcnt*12))

      !Cycle melange_data, constructing streams
      counter = 0
      DO i = 1,melange_data % NN
        IF(melange_data % owner(i) /= part) CYCLE
        counter = counter + 1
        PointEx(part) % SendGIDs(counter) = melange_data % GID(i)
        PointEx(part) % S(counter*12-11 : counter*12-6) = melange_data % UT%M(i*6-5:i*6)
        PointEx(part) % S(counter*12-5 : counter*12) = melange_data % UTM%M(i*6-5:i*6)
      END DO

      statcnt = statcnt + 1
      CALL MPI_ISend(PointEx(part) % SendGIDs, ppcnt, MPI_INTEGER, part, 301, &
           MPI_COMM_WORLD, stats(statcnt),ierr)
      statcnt = statcnt + 1
      CALL MPI_ISend(PointEx(part) % S, ppcnt*12, MPI_DOUBLE_PRECISION, part, 302, &
           MPI_COMM_WORLD, stats(statcnt),ierr)
    END DO
  END IF

  !Wait for the previous non-blocking sends, then reset stats
  CALL MPI_Waitall(statcnt, stats(1:statcnt), MPI_STATUSES_IGNORE, ierr)
  stats = MPI_REQUEST_NULL
  statcnt = 0

  !Put received values into UT, UTM
  !and adjust NRXF accordingly (NOTE: O(N^2))
  IF(my_ppcnt > 0) THEN
    DO i=1,my_ppcnt
      seek = recv_GIDs(i)
      matched = .FALSE.
      DO j=1,NRXF%NN
        IF(NRXF%GID(j) == seek) THEN
          UT%M(6*j-5:6*j) = recv_stream(i*12-11 : i*12-6)
          UTM%M(6*j-5:6*j) = recv_stream(i*12-5 : i*12)
          NRXF%M(1:3,j) = NRXF%M(1:3,j) - UT%M(6*j-5:6*j-3)
          matched = .TRUE.
          EXIT
        END IF
      END DO
      IF(.NOT. matched) THEN
        PRINT *,myid,' match particle failed, GID: ',seek
        CALL FatalError("Failed to match GID for particle")
      END IF
    END DO
  END IF

  ! Exchange bond info:
  ! melange_data % NANs is globally resolved (i.e. not partition specific)
  ! particle numbering for each connection. We also have per particle  %owner
  ! thus we can construct a list of bonds to send to each partition

  !Root cycles partitions, constructing list of beams
  IF(myid==0) THEN
    DO part=0,ntasks-1

      ppcnt = COUNT(melange_data % owner == part)
      IF(ppcnt == 0) CYCLE

      PointEx(part) % scount = 0
      IF(ALLOCATED(PointEx(part) % SendGIDs)) DEALLOCATE(PointEx(part) % SendGIDs)
      IF(ALLOCATED(PointEx(part) % S)) DEALLOCATE(PointEx(part) % S)
      
      ALLOCATE(PointEx(part) % SendGIDs(ppcnt * 6 * 2),& !send both particle GIDs
           PointEx(part) % S(ppcnt * 6)) !EFS

      !Cycle bonds, determining which ones to send & storing ID & EFS
      DO i=1,melange_data % NTOT
        N1 = melange_data % NANS(1,i)
        N2 = melange_data % NANS(2,i)
        send = (melange_data % owner(N1) == part) .OR. (melange_data % owner(N2) == part)
        IF(.NOT. send) CYCLE

        PointEx(part) % scount = PointEx(part) % scount + 1

        !Grow arrays if necessary
        IF(PointEx(part) % scount > SIZE(PointEx(part) % S)) THEN
          CALL ExpandIntArray(PointEx(part) % SendGIDs)
          CALL ExpandRealArray(PointEx(part) % S)
        END IF

        PointEx(part) % SendGIDs(PointEx(part) % scount * 2 - 1) = melange_data % GID(N1)
        PointEx(part) % SendGIDs(PointEx(part) % scount * 2 - 0) = melange_data % GID(N2)
        PointEx(part) % S(PointEx(part) % scount) = melange_data % EFS(i)
        
      END DO

      !Send beam count
      statcnt = statcnt + 1
      CALL MPI_ISend(PointEx(part) % scount, 1, MPI_INTEGER, part, 303,&
           MPI_COMM_WORLD, stats(statcnt),ierr)

      !Send beam particle GIDs
      statcnt = statcnt + 1
      CALL MPI_ISend(PointEx(part) % SendGIDs, PointEx(part) % scount*2, MPI_INTEGER, part, 304, &
           MPI_COMM_WORLD, stats(statcnt),ierr)

      statcnt = statcnt + 1
      CALL MPI_ISend(PointEx(part) % S, PointEx(part) % scount, MPI_DOUBLE_PRECISION, &
           part, 305, MPI_COMM_WORLD, stats(statcnt),ierr)

    END DO
  END IF

  !Recvs
  IF(my_ppcnt > 0) THEN

    DEALLOCATE(recv_stream,recv_GIDs)

    CALL MPI_Recv(my_beamcnt, 1, MPI_INTEGER, 0, 303, MPI_COMM_WORLD, stat, ierr)
    ALLOCATE(recv_GIDs(my_beamcnt * 2), recv_stream(my_beamcnt))

    CALL MPI_Recv(recv_GIDs, my_beamcnt*2, MPI_INTEGER, 0, 304, &
         MPI_COMM_WORLD, stat, ierr)
    CALL MPI_Recv(recv_stream, my_beamcnt, MPI_DOUBLE_PRECISION, 0, 305, &
         MPI_COMM_WORLD, stat, ierr)

    DO i=1,my_beamcnt
      matched = .FALSE.

      DO j=1,NTOT

        IF (.NOT. ANY(NRXF%GID(NANS(:,j))  == recv_GIDs(i*2-1))) CYCLE
        IF (.NOT. ANY(NRXF%GID(NANS(:,j))  == recv_GIDs(i*2))) CYCLE

        EFS(j) = recv_stream(i)
        matched = .TRUE.
        EXIT
      END DO

      IF(.NOT. matched) THEN
        PRINT *,myid,' failed: ',i,recv_GIDs(i*2-1:i*2), recv_stream(i)
        CALL FatalError("Programming Error: failed to match melange EFS value.")
      END IF
    END DO
  END IF

  CALL MPI_Waitall(statcnt, stats(1:statcnt), MPI_STATUSES_IGNORE, ierr)
  stats = MPI_REQUEST_NULL
  statcnt = 0

  
END SUBROUTINE PassMelangeData

END MODULE Melange
