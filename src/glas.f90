! *************************************************************************
! * The subroutine Initializefcc is based on a legacy code that used to be 
! * freely distributed on the internet. It has been modified to fit HiDEM. 
! * We are not aware who created the original version.
! *************************************************************************

MODULE Lattice

USE INOUT
USE Utils
USE TypeDefs

IMPLICIT NONE

CONTAINS

!defines and makes the FCC lattice - dense packing

!ip - returns number of nodes (in this partition)
!ntasks - how many cores
!myid - this partition id
SUBROUTINE FIBG3(SI,base,surf,origin,NN,NTOT,NANS,NRXF,NANPart,InvPartInfo,&
     neighcount,melange_data)

  IMPLICIT NONE

REAL(KIND=dp) :: surf(0:,0:),base(0:,0:)
REAL(KIND=dp),ALLOCATABLE :: melt(:,:)
Real(KIND=dp) :: x,y
REAL(KIND=dp) :: box,b,origin(2)
INTEGER :: NN,i,j,mask
INTEGER :: N1,N2,xk,yk,neighcount,NTOT
INTEGER, ALLOCATABLE :: NANS(:,:),NANPart(:)
TYPE(NRXF_t) :: NRXF
TYPE(InvPartInfo_t), ALLOCATABLE :: InvPartInfo(:)
TYPE(MelangeDataHolder_t) :: melange_data
TYPE(SimInfo_t) :: SI

!Open(300,file='mass.dat',STATUS='OLD')

!Not really used
ALLOCATE(melt(0:UBOUND(surf,1),0:UBOUND(surf,2)))
melt = SI%melt*0.0_dp

!box is never actualy used...
!b is used, which is box/l, so L never actually enters into this, so it's only
!used for the number of vertical layers
box=2.0d0**(2.0d0/3.0d0)*REAL(SI%LS,8) ! box size equal to fcc ground state
CALL Initializefcc(SI,NN,NTOT,NANS,NRXF,NANPart,InvPartInfo,&
     neighcount,box,surf,base,melt,origin,melange_data)

CLOSE(400)

END SUBROUTINE FIBG3

!---------------------------------------------------------------!


SUBROUTINE Initializefcc(SI,NN,NTOT,NANS,NRXF,NANPart, &
     InvPartInfo,neighcount,box,surf,base,melt,origin,melange_data)

  USE Melange

  IMPLICIT NONE

REAL(KIND=dp) :: surf(0:,0:),base(0:,0:),melt(0:,0:)
REAL(KIND=dp) :: b,x0(3,4),box,SCL,origin(2)
REAL(KIND=dp) :: gridminx, gridmaxx, gridminy, gridmaxy,T1,T2
REAL(KIND=dp) :: minx, miny
REAL(KIND=dp) :: z,x,y,sint,bint,mint,grid,lc,efficiency
REAL(KIND=dp) :: X1,X2,Y1,Y2,Z1,Z2,RC,rc_min,rc_max
REAL(KIND=dp), ALLOCATABLE :: xo(:,:),work_arr(:,:)

INTEGER i,j,k,n,ix,K1,k2,ip,glac_ip,NN,nb,xk,yk,nx,ny,nbeams,nbeams_mel,nprox_metis,ierr,pown
INTEGER NTOT,stats(4), local, part, neighcount
INTEGER, ALLOCATABLE ::  particles_L(:),NANS(:,:),NANPart(:),PartNN(:),&
     SendGIDs(:),ConnStream(:),RConnStream(:)
LOGICAL, ALLOCATABLE :: neighparts(:)
TYPE(Conn_t), ALLOCATABLE :: CN(:),CN_Glac(:),CN_metis(:),CN_Melange(:)

!metis stuff
INTEGER :: objval,counter
INTEGER, ALLOCATABLE :: metis_options(:), particlepart(:),counters(:)
INTEGER, POINTER :: vwgt=>NULL(),vsize=>NULL(),adjwgt=>NULL(),xadj(:),adjncy(:),countptr
REAL(KIND=dp), POINTER :: tpwgts=>NULL(),ubvec=>NULL()

!TYPE(NTOT_t) :: NTOT
TYPE(NRXF_t), TARGET :: NRXF
TYPE(NRXF_t) :: NRXFold
TYPE(InvPartInfo_t), ALLOCATABLE, TARGET :: InvPartInfo(:)
TYPE(MelangeDataHolder_t) :: melange_data
TYPE(SimInfo_t) :: SI

SCL = SI%SCL
grid = SI%grid

b=SCL*box/REAL(SI%LS,8)  ! the size of the unit cell is the box length divided by l
x0(:,:)=b/2.0d0; x0(:,1)=0.0d0; x0(3,2)=0.0d0; x0(2,3)=0.0d0; x0(1,4)=0.0d0 

IF(myid==0) THEN
  !Use bed and surf to determine the extent of the domain
  !This is not perfect, just checks that there is SCL dist between surf and bed+melt
  ALLOCATE(xo(3,NOMA))

  !Find the ice extent in raster coordinates
  gridminx = HUGE(gridminx); gridminy = HUGE(gridminy)
  gridmaxx = -HUGE(gridmaxx); gridmaxy = -HUGE(gridmaxy)
  DO i=LBOUND(surf,1), UBOUND(surf,1)

    DO j=LBOUND(surf,2), UBOUND(surf,2)
      IF(.NOT. surf(i,j) > (base(i,j) + melt(i,j))) CYCLE
      !more CYCLE?
      if (i < gridminx) gridminx = i
      if (i > gridmaxx) gridmaxx = i
      if (j < gridminy) gridminy = j
      if (j > gridmaxy) gridmaxy = j
    END DO

  END DO
  gridminx = gridminx - 1; gridminy = gridminy - 1
  gridmaxx = gridmaxx + 1; gridmaxy = gridmaxy + 1

  !How many boxes (~1.5 SCL) in the x (nx) and y (ny) direction
  nx = CEILING(((gridmaxx - gridminx) * grid) / b)
  ny = CEILING(((gridmaxy - gridminy) * grid) / b)

  !Minimum real coordinate with ice (i.e. location of a corner particle)
  minx = (gridminx * grid) + origin(1)
  miny = (gridminy * grid) + origin(2)

  !nx (ny) is the number of boxes in the x (y) direction
  ip=0
  DO i=1,nx !x step
    DO j=1,ny !y step
      DO k=-2,SI%LS !vertical layer
        DO k1=1,4

          x=(x0(1,k1) + minx + REAL(i-1)*b)
          y=(x0(2,k1) + miny + REAL(j-1)*b)
          z=x0(3,k1) + REAL(k-1)*b

          xk = FLOOR((x-origin(1))/grid)
          yk = FLOOR((y-origin(2))/grid)

          !just beyond edge of geom def
          IF(xk < 0 .OR. yk < 0) CYCLE
          IF(ANY(surf(xk:xk+1,yk:yk+1) == base(xk:xk+1,yk:yk+1)) .AND. SI%StrictDomain) CYCLE

          bint = InterpRast(X,Y,base,grid,origin,INTERP_MISS_FILL,1.0_dp) ! 1 > 0
          sint = InterpRast(X,Y,surf,grid,origin,INTERP_MISS_FILL,0.0_dp) ! i.e. no particle here
          mint = InterpRast(X,Y,melt,grid,origin,INTERP_MISS_FILL,0.0_dp)

          !             write(1510+myid,13) 40.0*x,40.0*y,bint,sint
          !             If (base(xk,yk).ne.0.0.and.base(xk,yk+1).ne.0.0.and.base(xk+1,yk).ne.0.0&
          ! .AND.base(xk+1,yk+1).NE.0.0) THEN

          !TODO - unhardcode this
          ! wl = SI%WL
          ! If (z.ge.bint+mint.and.z.le.sint.and.((sint-bint).gt.4.0*SCL.or.&
          ! (ABS(z-wl).LT.4.0*SCL.AND.bint.LT.wl))) THEN
          IF (z.GE.bint+mint .AND. z.LT.sint .AND. (sint-(bint+mint)).GT.SCL) THEN

            ! undercut shape functions
            ! lc=4420.0+1.5e-04*(x-3300.0)**2+0.42*exp((x-3700.0)/2.0e+02)
            ! UCV=lc-1500.0*exp(-(x-3500.0)**2/50000.0)
            ! If (y.lt.lc-SI%UC.or.y.gt.lc.or.z.gt.bint+3.0*SQRT(y-(lc-SI%UC)).or.z.ge.WL-20.0) then
            ! If (y.lt.lc-SI%UC.or.z.gt.bint+3.0*SQRT(y-(lc-SI%UC)).or.z.ge.WL-40.0) then
            ! If (y.lt.UCV.or.(z.gt.bint+3.0*sqrt(y-UCV)).or.z.ge.WL-40.0) then

            ip=ip+1
            !More space if req'd - TODO, pass this size back (in place of NOMA)
            IF(ip > SIZE(xo,2)) THEN
              ALLOCATE(work_arr(3,ip-1))
              work_arr = xo
              DEALLOCATE(xo)
              ALLOCATE(xo(3,ip*2))
              xo(:,1:ip-1) = work_arr(:,1:ip-1)
              DEALLOCATE(work_arr)
              !CALL FatalError("NOMA not large enough (too many points)")
              IF(DebugMode) CALL Warn("Doubling size of xo")
            END IF

            xo(1,ip) = x
            xo(2,ip) = y
            xo(3,ip) = z
            ! EndIF
            ! EndIF
          END IF
          !EndIF
        END DO
      END DO
    END DO
  END DO

  IF(DebugMode) PRINT *,myid,' Done generating particles: ',ip
  IF(DebugMode) PRINT *,myid,' Finding connections...'

  IF(PrintTimes) CALL CPU_TIME(T1)
  CALL FindBeams(xo, ip, SCL, CN_Glac, nbeams)
  IF(PrintTimes) CALL CPU_TIME(T2)

  IF(PrintTimes) PRINT *,myid,' Done finding connections: ',T2-T1,' secs'

  
  !If we are adding melange from a previous simulation, the equivalence between
  !proximity and 'bonded particles' breaks down, so we need to separately determine
  !which particles are proximal (for metis partitioning) and which are actually bonded
  !(for the model)
  IF(melange_data % active)  THEN

    !Get rid of any melange_data particles which overlap w/ the glacier geometry (xo(1:ip))
    CALL Prune_Melange(melange_data, xo, ip,SCL)

    PRINT *,myid,' debug about to exchange NN and expand'

    !Expand xo to make room for melange particles
    IF(ip+melange_data%NN > SIZE(xo,2)) CALL ExpandRealArray(xo,ip+melange_data%NN)

    PRINT *,myid,' debug about to fill xo'

    !Append melange particles to particle list
    DO i=1,melange_data%NN
      xo(1,ip+i) = melange_data%NRXF%M(1,i) + melange_data%UT%M(6*i-5)
      xo(2,ip+i) = melange_data%NRXF%M(2,i) + melange_data%UT%M(6*i-4)
      xo(3,ip+i) = melange_data%NRXF%M(3,i) + melange_data%UT%M(6*i-3)
    END DO

    glac_ip = ip
    ip = ip + melange_data%NN

    !Find nearby particles - this is important for partitioning
    !Because melange (unlike glacier) begins largely broken, connection info for
    !metis needs *proximity* rather than beams
    CALL FindNNearest(xo,ip,12,SCL, CN_metis)
    nprox_metis = 12 * ip

    !Ensure every bond goes both ways
    CALL EnsureDualGraph(CN_Metis)

    !Detect disconnected groups of particles & connect them
    CALL ConnectRegions(xo, CN_metis, SCL)

    !Generate *actual* bond info for melange particles
    CALL MelangeBonds(melange_data, glac_ip, CN_melange, nbeams_mel)

  ELSE

    ALLOCATE(CN_Metis(ip))
    DO i=1,ip
      ALLOCATE(CN_metis(i) % Conn(CN_Glac(i) % NCN))
      CN_metis(i) % NCN = CN_Glac(i) % NCN
      CN_metis(i) % Conn = CN_Glac(i) % Conn
    END DO
    nprox_metis = nbeams

    glac_ip = ip

  END IF


  ALLOCATE(ParticlePart(ip),PartNN(ntasks))

  PRINT *,'About to metis'
  ALLOCATE(metis_options(40),xadj(ip+1),adjncy(nprox_metis*2))
  
  !Put particle connections into CRS format
  xadj = 0
  adjncy = 0
  counter = 0
  xadj(1) = 1
  DO i=1,ip
    DO j=1,CN_metis(i) % NCN
      counter = counter + 1
      adjncy(counter) = CN_metis(i) % Conn(j)
    END DO
    xadj(i+1) = counter+1
  END DO

  !Manually deallocate these here to reduce RAM high water mark
  DO i=1,ip
    DEALLOCATE(CN_metis(i) % Conn)
  END DO
  DEALLOCATE(CN_metis)


  !TODO - ensure compatibility (REAL and INT size)
  CALL METIS_SetDefaultOptions(metis_options)
  metis_options(18) = 1 !fortran array numbering

  !Use METIS to partition the particles
  CALL METIS_PartGraphKway(ip,1,xadj,adjncy,vwgt,vsize,adjwgt,ntasks,&
       tpwgts, ubvec, metis_options,objval,particlepart)

  particlepart = particlepart - 1 !mpi partitions are zero indexed

  DO i=1,ntasks
    PartNN(i) = COUNT(particlepart == i-1)
  END DO

  !Construct all local nodenums:
  ALLOCATE(counters(ntasks),&
       particles_L(ip))
  counters = 0
  particles_L = 0

  DO i=1,ip
    counters(particlepart(i)+1) = counters(particlepart(i)+1) + 1
    particles_L(i) = counters(particlepart(i)+1)
  END DO

  IF(DebugMode) THEN
    PRINT *,'obvjal: ',objval
    PRINT *,'particlepart: ',MINVAL(particlepart), MAXVAL(particlepart), &
         COUNT(particlepart==0),COUNT(particlepart==ntasks-1)
    PRINT *,'--- METIS SUCCESS ---'
  END IF
  DEALLOCATE(xadj,adjncy)

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!Root sends particle & connection info to other processes:

!Send total particle count
CALL MPI_BCast(ip,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!Send partition particle count
CALL MPI_Scatter(PartNN,1,MPI_INTEGER,NN,1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!Allocate the structure holding the point data
CALL PointDataInit(NRXF, NN, part_expand)
ALLOCATE(work_arr(NN,3),&
     neighparts(0:ntasks-1))

neighparts = .FALSE.

!Send/recv partition particle coords (xo) and GIDs
stats = MPI_REQUEST_NULL
CALL MPI_IRECV(NRXF%GID(1:NN),NN,MPI_INTEGER,0,121,MPI_COMM_WORLD,stats(1),ierr)
CALL MPI_IRECV(work_arr(1:NN,1),NN,MPI_DOUBLE_PRECISION,0,122,MPI_COMM_WORLD,stats(2),ierr)
CALL MPI_IRECV(work_arr(1:NN,2),NN,MPI_DOUBLE_PRECISION,0,123,MPI_COMM_WORLD,stats(3),ierr)
CALL MPI_IRECV(work_arr(1:NN,3),NN,MPI_DOUBLE_PRECISION,0,124,MPI_COMM_WORLD,stats(4),ierr)

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!Cycle partitions sending particles
IF(myid==0) THEN
  
  ALLOCATE(SendGIDs(MAXVAL(PartNN)))

  DO i=1,ntasks

    counter = 0
    SendGIDs = 0
    DO j=1,ip
      IF(ParticlePart(j) == i-1) THEN
        counter = counter + 1
        IF(counter > PartNN(i)) CALL FatalError("Programming error, sending too many points")
        SendGIDs(counter) = j
      END IF
    END DO

    CALL MPI_SEND(SendGIDs(1:PartNN(i)),PartNN(i),MPI_INTEGER, i-1,121,MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(xo(1,SendGIDs(1:PartNN(i))),PartNN(i),MPI_DOUBLE_PRECISION, i-1,122,&
         MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(xo(2,SendGIDs(1:PartNN(i))),PartNN(i),MPI_DOUBLE_PRECISION, i-1,123,&
         MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(xo(3,SendGIDs(1:PartNN(i))),PartNN(i),MPI_DOUBLE_PRECISION, i-1,124,&
         MPI_COMM_WORLD, ierr)
  END DO
END IF

CALL MPI_Waitall(4, stats, MPI_STATUSES_IGNORE, ierr)

NRXF%A(1,1:NN) = work_arr(1:NN,1)
NRXF%A(2,1:NN) = work_arr(1:NN,2)
NRXF%A(3,1:NN) = work_arr(1:NN,3)

NRXF%PartInfo(1,1:NN) = myid
DO i=1,NN
  NRXF%PartInfo(2,i) = i
END DO

DEALLOCATE(work_arr)

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


!Cycle partitions sending connection info - send local id and partition for other particle.

!Set up the non-blocking recv
ALLOCATE(RConnStream(NN*12*2))
stats = MPI_REQUEST_NULL
CALL MPI_IRECV(RConnStream, NN*12*2, MPI_INTEGER, 0, 125, MPI_COMM_WORLD, stats(1), ierr)

IF(myid==0) THEN

  !Need at least to send a list of local ID & partition (and probably GID too, why not?)
  !Fill ConnStream with neighbouring particle localID and Partition
  DO i=1,ntasks

    ALLOCATE(ConnStream(PartNN(i)*12*2))
    ConnStream = 0
    counter = 0

    DO j=1,ip

      IF(ParticlePart(j) /= i-1) CYCLE
      counter = counter + 1

      IF(j <= glac_ip) THEN
        DO k=1,CN_Glac(j) % NCN
          ConnStream((counter-1)*12*2 + ((k-1) * 2) + 1) = &
               particles_L(CN_Glac(j) % Conn(k))  !local
          ConnStream((counter-1)*12*2 + ((k-1) * 2) + 2) = &
               ParticlePart(CN_Glac(j) % Conn(k)) !part
        END DO
      ELSE
        ix = j - glac_ip
        DO k=1,CN_Melange(ix) % NCN
          ConnStream((counter-1)*12*2 + ((k-1) * 2) + 1) = &
               particles_L(CN_Melange(ix) % Conn(k))  !local
          ConnStream((counter-1)*12*2 + ((k-1) * 2) + 2) = &
               ParticlePart(CN_Melange(ix) % Conn(k)) !part
        END DO
      END IF

    END DO

    CALL MPI_SEND(ConnStream, PartNN(i)*2*12, MPI_INTEGER, i-1, 125, MPI_COMM_WORLD, ierr)
    DEALLOCATE(ConnStream)

  END DO

  DEALLOCATE(CN_Glac,ParticlePart)
END IF

CALL MPI_Waitall(1, stats, MPI_STATUSES_IGNORE, ierr)

IF(SIZE(RConnStream) /= 12*2*NN) CALL FatalError("Programming Error: RConnStream wrong size!")

!Construct an array of all connections:
!(this differs from CN, which is the per-node connection info
!Note - loop twice, first count, then allocate and fill
!NTOT = total beams this partition owns + shares
counter = 0
DO k=1,2
  IF(k==2) THEN
    ALLOCATE(NANS(2,counter),NANPart(counter))
    NTOT = counter
    counter = 0
    NANPart = -1
  END IF

  DO i=1,NN
    DO j=1,12

      !Array location of this neighbour
      n = ((i-1) * 12 * 2) + ((j-1) * 2) + 1
      local = RConnStream(n)
      part = RConnStream(n+1)

      IF(local == 0) EXIT
      IF(part /= myid) neighparts(part) = .TRUE.

      !Count each beam only once - except across boundaries, both need to count
      IF(local > i .OR. part /= myid) THEN
        counter = counter + 1
        IF(k==2) THEN
          NANS(1,counter) = local !Note - folows convention N1 = other part
          NANS(2,counter) = i  !NANS = their/ourNN, ourNN
          NANpart(counter) = part !NANPart is the otherPart
        END IF
      END IF
    END DO
  END DO
END DO

DEALLOCATE(RConnStream)

neighcount = COUNT(neighparts)

!Construct lookup arrays in NRXF for across partition beams
!Note this fills the %PartInfo but doesn't actually fill the points
!This is done by ExchangeConnPoints
CALL InvPartInfoInit(InvPartInfo, neighparts)

counter = NRXF%cstrt - 1
DO n=0,ntasks-1
  IF(.NOT. neighparts(n)) CYCLE
  IF(DebugMode) PRINT *,myid, 'neighbour is ',InvPartInfo(n) % NID
  DO i=1,NTOT
    IF (NANPart(i) /= InvPartInfo(n) % NID) CYCLE
    IF (ANY(InvPartInfo(n) % ConnIDs == NANS(1,i))) CYCLE

    countptr => InvPartInfo(n) % CCount
    countptr = countptr + 1
    counter = counter + 1

    IF(countptr > SIZE(InvPartInfo(n) % ConnIDs)) THEN
      CALL ExpandIntArray(InvPartInfo(n) % ConnIDs)
      CALL ExpandIntArray(InvPartInfo(n) % ConnLocs)
    END IF

    InvPartInfo(n) % ConnIDs(countptr) = NANS(1,i)
    InvPartInfo(n) % ConnLocs(countptr) = counter
    
    IF(counter > SIZE(NRXF%PartInfo,2)) CALL ResizePointData(NRXF,1.5_8, do_C=.TRUE.,do_P=.FALSE.)
    NRXF%PartInfo(1,counter) = InvPartInfo(n) % NID
    NRXF%PartInfo(2,counter) = NANS(1,i)
    NRXF%NC = NRXF%NC + 1
  END DO
END DO
NRXF%pstrt = NRXF%cstrt + NRXF%NC

!Change NANS to array loc:
DO i=1,NTOT
  IF(NANpart(i) == myid) CYCLE
  n = NANpart(i)

  IF(InvPartInfo(n)%CCount == 0) THEN
    PRINT *,'programming error'
    stop
  END IF

  DO k=1,InvPartInfo(n)%CCount
    IF(InvPartInfo(n) % ConnIDs(k) == NANS(1,i)) THEN
      ! IF(DebugMode) PRINT *,myid,' setings nans ',i,' to ',&
      !      InvPartInfo(n) % ConnLocs(k), NANPart(i), NANS(1,i)

      NANS(1,i) = InvPartInfo(n) % ConnLocs(k)
      EXIT
    END IF
  END DO
END DO

! !TODO - need comms here
! OPEN(UNIT=117+myid,FILE=TRIM(SI%wrkdir)//'/FS'//na(myid),STATUS='UNKNOWN')
! DO I=1,NTOT
!   WRITE(117+myid,*) NANS(1,I),NANS(2,I),NANS(3,I),&
!        NRXF%M(1,NANS(1,I)),NRXF%M(2,NANS(1,I)),&
!        NRXF%M(3,NANS(1,I)),NRXF%M(1,NANS(2,I)),&
!        NRXF%M(2,NANS(2,I)),NRXF%M(3,NANS(2,I)),EFS % M(I)
! END DO
! CLOSE (117+myid)



CALL ExchangeConnPoints(NANS, NRXF, InvPartInfo)

!Write out my particle to nodfil
OPEN(510+myid,file=TRIM(SI%wrkdir)//'/'//TRIM(SI%runname)//'_NODFIL2'//na(myid))
DO i=1,NN
  WRITE(510+myid,'(I8,4F16.8)') i,NRXF%A(:,i),1.0
END DO
! DO i=NRXF%cstrt, NRXF%cstrt + NRXF%NC - 1
!   WRITE(510+myid,12) i,NRXF%A(:,i),1.0,NRXF%PartInfo(:,i)
! END DO
CLOSE(510+myid)

IF(DebugMode .AND. .FALSE.) THEN
  rc_min = HUGE(rc_min)
  rc_max = -HUGE(rc_max)

  DO i=1,NTOT
    RC = ((NRXF%A(1,NANS(1,i)) - NRXF%A(1,NANS(2,i)))**2.0 + &
          (NRXF%A(2,NANS(1,i)) - NRXF%A(2,NANS(2,i)))**2.0 + &
          (NRXF%A(3,NANS(1,i)) - NRXF%A(3,NANS(2,i)))**2.0) ** 0.5
    rc_min = MIN(rc_min, RC)
    rc_max = MAX(rc_max, RC)

    IF(NRXF%A(1,NANS(1,i)) == 0.0 .AND. &
         NRXF%A(2,NANS(1,i)) == 0.0 .AND. &
         NRXF%A(3,NANS(1,i)) == 0.0) &
         PRINT *,myid,' uninit point 1: ',i,NANS(1,i),NRXF%PartInfo(1,NANS(1,i))

    IF(NRXF%A(1,NANS(2,i)) == 0.0 .AND. &
         NRXF%A(2,NANS(2,i)) == 0.0 .AND. &
         NRXF%A(3,NANS(2,i)) == 0.0) PRINT *,myid,' uninit point 2: ',i,NANS(2,i)

    IF(rc < 50.0) PRINT *,myid,' short conn: ',rc,i,NANPart(i),NANS(1,i),&
         NANS(2,i),NRXF%A(1,NANS(1,i)),NRXF%A(1,NANS(2,i))
    IF(rc > 150.0) PRINT *,myid,' long conn: ',rc,i,NANPart(i),NANS(1,i),&
         NANS(2,i),NRXF%A(1,NANS(1,i)),NRXF%A(1,NANS(2,i))
  END DO

  PRINT *,myid,' max NRXF rc: ',rc_max, rc_min
END IF

!We return:
! Our nodes initial locs  - NRXF
! Our node global numbers - particles_G
! Our connections - NCN, CN
! Our neighparts - neighparts

END SUBROUTINE Initializefcc

SUBROUTINE GetBBoxes(NRXF, UT, NN, IsOutlier, BBox, PBBox)

   IMPLICIT NONE

   INTEGER ::  NN
   TYPE(UT_t) :: UT
   TYPE(NRXF_t) :: NRXF
   REAL(KIND=dp) :: BBox(6),PBBox(6,0:ntasks)
   LOGICAL, ALLOCATABLE :: IsOutlier(:)
   !-----------------
   INTEGER :: i,ierr
   REAL(KIND=dp) :: workarr(6*ntasks),X,Y,Z
   REAL(KIND=dp) :: minx,maxx,miny,maxy,minz,maxz

   minx = HUGE(minx)
   miny = HUGE(miny)
   minz = HUGE(minz)
   maxx = -HUGE(maxx)
   maxy = -HUGE(maxy)
   maxz = -HUGE(maxz)
   DO i=1,NN
     IF(IsOutlier(i)) CYCLE
     X=NRXF%M(1,i)+UT%M(6*i-5)
     Y=NRXF%M(2,i)+UT%M(6*i-4)
     Z=NRXF%M(3,i)+UT%M(6*i-3)
     minx = MIN(minx, x)
     miny = MIN(miny, y)
     minz = MIN(minz, z)
     maxx = MAX(maxx, x)
     maxy = MAX(maxy, y)
     maxz = MAX(maxz, z)
   END DO

   BBox(1) = minx
   BBox(2) = maxx
   BBox(3) = miny
   BBox(4) = maxy
   BBox(5) = minz
   BBox(6) = maxz

   CALL MPI_AllGather(BBox,6,MPI_DOUBLE_PRECISION,workarr,6,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)

   DO i=0,ntasks-1
     PBBox(:,i) = workarr((i*6)+1:(i+1)*6)
   END DO

 END SUBROUTINE GetBBoxes

 SUBROUTINE FindNeighbours(PBBox,PartIsNeighbour,SCL)
   REAL(KIND=dp) :: PBBox(:,0:),SCL
   LOGICAL :: PartIsNeighbour(0:ntasks-1)
   !-------------------------------
   INTEGER :: i
   REAL(KIND=dp) :: BBox(6),Buffer

   BBox = PBBox(:,myid)
   Buffer = 4.0*SCL

   PartIsNeighbour = .FALSE.


   !Cycle each partition
   DO i=0,ntasks-1
     IF(i==myid) CYCLE
     IF(BBox(1) - Buffer > PBBox(2,i)) CYCLE
     IF(BBox(3) - Buffer > PBBox(4,i)) CYCLE
     IF(BBox(5) - Buffer > PBBox(6,i)) CYCLE
     IF(BBox(2) + Buffer < PBBox(1,i)) CYCLE
     IF(BBox(4) + Buffer < PBBox(3,i)) CYCLE
     IF(BBox(6) + Buffer < PBBox(5,i)) CYCLE

     PartIsNeighbour(i) = .TRUE.
     IF(DebugMode) PRINT *,myid,' potential neighbour: ',i
   END DO

 END SUBROUTINE FindNeighbours

!Determine and pass point information between partitions (based on beams/connections)
SUBROUTINE ExchangeConnPoints(NANS, NRXF, InvPartInfo, UT, UTM, passNRXF)

  
  INTEGER, ALLOCATABLE :: NANS(:,:)
  TYPE(NRXF_t) :: NRXF
  TYPE(InvPartInfo_t) :: InvPartInfo(0:)
  TYPE(UT_t), OPTIONAL :: UT, UTM
  LOGICAL, OPTIONAL :: passNRXF
  !--------------------
  REAL(KIND=dp) :: T1,T2
  INTEGER :: i,j,id,counter,loc,neigh,getcount,sendcount,ierr
  INTEGER, ALLOCATABLE :: stats(:)
  TYPE(PointEx_t), ALLOCATABLE :: PointEx(:)
  LOGICAL :: doNRXF=.TRUE.

  IF(PrintTimes) CALL CPU_TIME(T1)

  IF(PRESENT(passNRXF)) doNRXF = passNRXF

  ALLOCATE(PointEx(0:ntasks-1),stats(0:ntasks*2-1))
  stats = MPI_REQUEST_NULL

  IF(DebugMode) PRINT *,myid, ' ExchangeConnPoints Checkpoint:  0'

  !First exchange the number of particles we need from other partitions
  DO i=0,ntasks-1
    IF(InvPartInfo(i) % ccount == 0) CYCLE
    neigh = InvPartInfo(i) % NID
    PointEx(i) % rcount = InvPartInfo(i) % CCount

    CALL MPI_ISend(PointEx(i) % rcount,1,MPI_INTEGER,neigh,&
         198,MPI_COMM_WORLD,stats(i*2),ierr)

    CALL MPI_IRecv(PointEx(i) % scount,1,MPI_INTEGER,neigh,&
         198,MPI_COMM_WORLD,stats(i*2+1),ierr)

  END DO

  IF(DebugMode) PRINT *,myid, ' ExchangeConnPoints Checkpoint:  1'

  !Wait for the previous non-blocking sends, then reset stats
  CALL MPI_Waitall(ntasks*2, stats, MPI_STATUSES_IGNORE, ierr)
  stats = MPI_REQUEST_NULL

  !Now send the particle IDs we need to receive from other partitions (RecvIDs)
  DO i=0,ntasks-1
    IF(InvPartInfo(i) % ccount == 0) CYCLE
    neigh = InvPartInfo(i) % NID

    sendcount = PointEx(i) % scount
    getcount = PointEx(i) % rcount

    PointEx(i) % partid = neigh
    ALLOCATE(PointEx(i) % SendIDs(sendcount),&
         PointEx(i) % RecvIDs(getcount),&
         PointEx(i) % SendGIDs(sendcount),&
         PointEx(i) % RecvGIDs(getcount),&
         PointEx(i) % S(3*sendcount),&
         PointEx(i) % R(3*getcount))

    PointEx(i) % RecvIDs = InvPartInfo(i) % ConnIDs(1:getcount)

    CALL MPI_ISend(PointEx(i) % RecvIDs,getcount,MPI_INTEGER,neigh,&
         199,MPI_COMM_WORLD,stats(i*2),ierr)

    CALL MPI_IRecv(PointEx(i) % SendIDs,sendcount,MPI_INTEGER,neigh,&
         199,MPI_COMM_WORLD,stats(i*2+1),ierr)

  END DO

  IF(DebugMode)  PRINT *,myid, ' ExchangeConnPoints Checkpoint:  2'

  !Wait for the previous non-blocking sends, then reset stats
  CALL MPI_Waitall(ntasks*2, stats, MPI_STATUSES_IGNORE, ierr)
  stats = MPI_REQUEST_NULL

  !Store info about which points we will send to each partition
  DO i=0,ntasks-1
    IF(InvPartInfo(i) % ccount == 0) CYCLE 
    sendcount = PointEx(i) % scount

    IF(sendcount == 0) CALL FatalError('Programming Error: scount == 0')

    InvPartInfo(i) % sccount = sendcount
    IF(SIZE(InvPartInfo(i) % SConnIDs) < sendcount) &
         CALL ExpandIntArray(InvPartInfo(i) % SConnIDs,sendcount)
    InvPartInfo(i) % SConnIDs(1:sendcount) = PointEx(i) % SendIDs(1:sendcount)
  END DO

  IF(doNRXF) THEN
    !Now send the actual point locations
    DO i=0,ntasks-1
      IF(InvPartInfo(i) % ccount == 0) CYCLE
      neigh = PointEx(i) % partid
      sendcount = PointEx(i) % scount

      DO j=1,sendcount
        id = PointEx(i) % SendIDs(j)
        PointEx(i) % S(j*3 - 2) = NRXF%M(1,id)
        PointEx(i) % S(j*3 - 1) = NRXF%M(2,id)
        PointEx(i) % S(j*3 - 0) = NRXF%M(3,id)
      END DO

      CALL MPI_ISend(PointEx(i) % S, sendcount*3, MPI_DOUBLE_PRECISION,neigh,&
           200,MPI_COMM_WORLD,stats(i*2), ierr)

      getcount = PointEx(i) % rcount
      CALL MPI_IRecv(PointEx(i) % R, getcount*3, MPI_DOUBLE_PRECISION,neigh,&
           200, MPI_COMM_WORLD,stats(i*2+1), ierr)

    END DO

    IF(DebugMode) PRINT *,myid, ' ExchangeConnPoints Checkpoint:  3'

    !Wait for the previous non-blocking sends, then reset stats
    CALL MPI_Waitall(ntasks*2, stats, MPI_STATUSES_IGNORE, ierr)
    stats = MPI_REQUEST_NULL

    !Put the points in NRXF
    DO i=0,ntasks-1
      IF(InvPartInfo(i) % ccount == 0) CYCLE

      neigh = PointEx(i) % partid
      getcount = PointEx(i) % rcount

      DO j=1,getcount
        loc = InvPartInfo(i) % ConnLocs(j)
        NRXF%A(1,loc) = PointEx(i) % R(j*3 - 2)
        NRXF%A(2,loc) = PointEx(i) % R(j*3 - 1)
        NRXF%A(3,loc) = PointEx(i) % R(j*3 - 0)
        NRXF%GID(loc) = PointEx(i) % RecvGIDs(j) !TODO!!!!
        !IF(DebugMode) PRINT *,myid,' neigh: ',neigh,' loc: ',loc
      END DO
    END DO
  END IF

  IF(DebugMode) PRINT *,myid, ' ExchangeConnPoints Checkpoint:  4'

  IF(PRESENT(UT)) THEN
    
    !Send and receive UT
    DO i=0,ntasks-1
      IF(InvPartInfo(i) % ccount == 0) CYCLE
      neigh = PointEx(i) % partid
      sendcount = PointEx(i) % scount
      getcount = PointEx(i) % rcount


      DEALLOCATE(PointEx(i) % S, PointEx(i) % R)
      ALLOCATE(PointEx(i) % S(12*sendcount),&
           PointEx(i) % R(12*getcount))


      PointEx(i) % S = 0.0
      PointEx(i) % R = 0.0

      DO j=1,sendcount
        id = PointEx(i) % SendIDs(j)
        PointEx(i) % S(j*12-11 : j*12-6) = UT%M(6*id-5 : 6*id)
        PointEx(i) % S(j*12-5 : j*12) = UTM%M(6*id-5 : 6*id)
      END DO

      CALL MPI_ISend(PointEx(i) % S, sendcount*12, MPI_DOUBLE_PRECISION,neigh,&
           201,MPI_COMM_WORLD,stats(i*2), ierr)

      CALL MPI_IRecv(PointEx(i) % R, getcount*12, MPI_DOUBLE_PRECISION,neigh,&
           201,MPI_COMM_WORLD,stats(i*2+1), ierr)
    END DO

    !Wait for the previous non-blocking sends, then reset stats
    CALL MPI_Waitall(ntasks*2, stats, MPI_STATUSES_IGNORE, ierr)
    stats = MPI_REQUEST_NULL

    IF(DebugMode) PRINT *,myid, ' ExchangeConnPoints Checkpoint:  5'

    !Put the points in UT
    DO i=0,ntasks-1
      IF(InvPartInfo(i) % ccount == 0) CYCLE

      neigh = PointEx(i) % partid
      getcount = PointEx(i) % rcount
      
      DO j=1,getcount
        loc = InvPartInfo(i) % ConnLocs(j)
        IF(6*loc > UBOUND(UT%A,1)) CALL FatalError(&
             "ExchangeConnPoints: Insufficient space in UT array")

        UT%A(6*loc-5 : 6*loc) = PointEx(i) % R(j*12-11 : j*12-6)
        UTM%A(6*loc-5 : 6*loc) = PointEx(i) % R(j*12-5 : j*12)
      END DO
    END DO
  END IF

  IF(PrintTimes) THEN
    CALL CPU_TIME(T2)
    PRINT *,myid,'Exchange Conn Points took: ',T2-T1,' secs'
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


END SUBROUTINE ExchangeConnPoints

!Determine and pass point information between partitions (based on proximity)
!Note this differs from ExchangeConnPoints, which 1) knows a priori which particle info 
!to pass and 2) passes the same particles each time. In contrast, we must here
!1) determine which particles are req'd, and also 2) when new particles are passed
!we need to also pass NRXF (only once!)
SUBROUTINE ExchangeProxPoints(NRXF, UT, UTM, NN, SCL, PBBox, InvPartInfo, PartIsNeighbour)

  USE INOUT
  USE UTILS


  TYPE(NRXF_t) :: NRXF
  TYPE(UT_t) :: UT, UTM
  INTEGER :: NN
  REAL(KIND=dp) :: SCL, PBBox(6,0:ntasks-1)
  LOGICAL :: PartIsNeighbour(0:ntasks-1)
  TYPE(InvPartInfo_t) :: InvPartInfo(0:)
  !--------------------
  REAL(KIND=dp) :: T1, T2, tstrt,tend
  INTEGER :: i,j,k,id,new_id,put_loc,put_loc_init,cnt,cnt2,rmcnt,neigh,&
       getcount,rmcount,loc,stat(MPI_STATUS_SIZE),ierr
  INTEGER, ALLOCATABLE :: WorkInt(:),stats(:)
  TYPE(PointEx_t), ALLOCATABLE :: PointEx(:),NRXFPointEx(:)
  LOGICAL, ALLOCATABLE :: IsConnected(:),PrevProx(:)
  LOGICAL :: loss

  if(PrintTimes) CALL CPU_TIME(tstrt)

  ALLOCATE(PointEx(0:ntasks-1),NRXFPointEx(0:ntasks-1), &
       stats(0:ntasks*4-1),IsConnected(NN),PrevProx(NN))
  stats = MPI_REQUEST_NULL

  DO i=0,ntasks-1
    IF(.NOT. PartIsNeighbour(i)) CYCLE
    ALLOCATE(PointEx(i) % SendIDs(NN/2)) !TODO - better guess for this
    PointEx(i) % partid = i
    PointEx(i) % scount = 0
    PointEx(i) % rcount = 0

    ALLOCATE(NRXFPointEx(i) % SendIDs(NN/2)) !and this
    NRXFPointEx(i) % partid = i
    NRXFPointEx(i) % scount = 0
    NRXFPointEx(i) % rcount = 0
  END DO

  IF(PrintTimes) CALL CPU_TIME(T1)

  !Count points in each BB
  !Strategy here:
  !If point is connected - don't send
  !If point wasn't already sent - send NRXF
  !InvPartInfo % ProxIDs, ProxLocs should contain only *unconnected* points
  !   i.e. they're not in ConnLoc, ConnIDs
  DO i=0,ntasks-1
    IF(.NOT. PartIsNeighbour(i)) CYCLE
    IF(.NOT. ALLOCATED(InvPartInfo(i) % SProxIDs)) CALL NewInvPartInfo(InvPartInfo, i)

    IF(DebugMode) PRINT *,myid,' sends to ',i,'count ',InvPartInfo(i) % sccount,' sum: ',&
         SUM(InvPartInfo(i) % SConnIDs(1:InvPartInfo(i) % sccount))

    !Construct these temporary arrays to avoid repeated lookups
    IsConnected = .FALSE.
    DO j=1,InvPartInfo(i) % sccount
      IsConnected(InvPartInfo(i) % SConnIDs(j)) = .TRUE.
    END DO

    PrevProx = .FALSE.
    DO j=1,InvPartInfo(i) % spcount
      PrevProx(InvPartInfo(i) % SProxIDs(j)) = .TRUE.
    END DO

    cnt = 0
    cnt2 = 0
    rmcount = 0
    DO j=1,NN

      !If this particle is part of a beam shared w/ this partition, it's already been sent
      IF(IsConnected(j)) CYCLE

      !If the point is near the other partition's bbox, add to list
      IF(PInBBox(j,NRXF,UT,PBBox(:,i),SCL*2.0)) THEN
        PointEx(i) % scount = PointEx(i) % scount + 1
        cnt = PointEx(i) % scount

        IF(cnt > SIZE(PointEx(i) % SendIDs)) CALL ExpandIntArray(&
             PointEx(i) % SendIDs)

        PointEx(i) % SendIDs(cnt) = j

        !If point hasn't previously been sent, add to list
        IF(.NOT. PrevProx(j)) THEN
          NRXFPointEx(i) % scount = NRXFPointEx(i) % scount + 1
          cnt2 = NRXFPointEx(i) % scount
          IF(cnt2 > SIZE(NRXFPointEx(i) % SendIDs)) CALL ExpandIntArray(&
               NRXFPointEx(i) % SendIDs)
          NRXFPointEx(i) % SendIDs(cnt2) = j

          InvPartInfo(i) % spcount = InvPartInfo(i) % spcount + 1
          IF(InvPartInfo(i) % spcount > SIZE(InvPartInfo(i) % SProxIDs)) &
               CALL ExpandIntArray(InvPartInfo(i) % SProxIDs,fill_in=-1)
          InvPartInfo(i) % SProxIDs(InvPartInfo(i) % spcount) = j
          PrevProx(j) = .TRUE.
        END IF
      ELSE IF(PrevProx(j)) THEN
        IF(DebugMode) PRINT *,myid,' debug, point ',j,' no longer near partition ',i
        PrevProx(j) = .FALSE.
        rmcount = rmcount + 1
      END IF
    END DO

    IF(rmcount > 0 .AND. DebugMode) PRINT *,myid,' rmcount: ',rmcount,i

    !Mark any nodes previously sent but no longer needed
    rmcount = 0
    DO j=1,InvPartInfo(i) % spcount
      IF(.NOT. PrevProx(InvPartInfo(i) % SProxIDs(j))) THEN
        IF(DebugMode) PRINT *,myid,' debug particle ',InvPartInfo(i) % SProxIDs(j),&
             ' no longer neighbour of ',i,NRXF%A(:,InvPartInfo(i) % SProxIDs(j))
        InvPartInfo(i) % SProxIDs(j) = -1
        rmcount = rmcount + 1
      END IF
    END DO

    !Remove them from the array
    InvPartInfo(i) % spcount = InvPartInfo(i) % spcount - rmcount
    InvPartInfo(i) % SProxIDs(1:InvPartInfo(i) % spcount) = &
         PACK(InvPartInfo(i) % SProxIDs,InvPartInfo(i) % SProxIDs /= -1)
    InvPartInfo(i) % SProxIDs(InvPartInfo(i) % spcount+1:UBOUND(InvPartInfo(i) % SProxIDs,1)) = -1

    IF(DebugMode) PRINT *,myid,' sending to ',i,' prox: ', &
         cnt, ' already sent conn: ',COUNT(IsConnected), ' prev sent nrxf: ',NRXFPointEx(i) % scount

    ALLOCATE(PointEx(i) % S(12 * PointEx(i) % scount))
    CALL MPI_ISend(PointEx(i) % scount,1,MPI_INTEGER,i,&
         201,MPI_COMM_WORLD,stats(i*4),ierr)

    CALL MPI_ISend(PointEx(i) % SendIDs(1:cnt),cnt,MPI_INTEGER,i,&
         202,MPI_COMM_WORLD,stats(i*4+1),ierr)

    ALLOCATE(NRXFPointEx(i) % S(3 * NRXFPointEx(i) % scount))
    CALL MPI_ISend(NRXFPointEx(i) % scount,1,MPI_INTEGER,i,&
         203,MPI_COMM_WORLD,stats(i*4+2),ierr)
    CALL MPI_ISend(NRXFPointEx(i) % SendIDs(1:cnt2),cnt2,MPI_INTEGER,i,&
         204,MPI_COMM_WORLD,stats(i*4+3),ierr)
  END DO


  !Receive count and IDs sent above
  !Note: this code can throw up array bounds checking errors when receiving zero points, but
  ! its not an issue.
  DO i=0,ntasks-1
    IF(.NOT. PartIsNeighbour(i)) CYCLE

    CALL MPI_Recv(PointEx(i) % rcount, 1, MPI_INTEGER, i,201,MPI_COMM_WORLD, stat, ierr)
    cnt = PointEx(i) % rcount

    ALLOCATE(PointEx(i) % R(12 * cnt),&
         PointEx(i) % RecvIDs(cnt))

    CALL MPI_Recv(PointEx(i) % RecvIDs(1:cnt),cnt, MPI_INTEGER, i,202,MPI_COMM_WORLD, stat, ierr)

    CALL MPI_Recv(NRXFPointEx(i) % rcount, 1, MPI_INTEGER, i,203,MPI_COMM_WORLD, stat, ierr)
    cnt = NRXFPointEx(i) % rcount

    ALLOCATE(NRXFPointEx(i) % R(3 * cnt),&
         NRXFPointEx(i) % RecvIDs(cnt))

    CALL MPI_Recv(NRXFPointEx(i) % RecvIDs(1:cnt),cnt, MPI_INTEGER, i,204,MPI_COMM_WORLD, stat, ierr)
  END DO

  !Wait for the previous non-blocking sends, then reset stats
  CALL MPI_Waitall(ntasks*4, stats, MPI_STATUSES_IGNORE, ierr)
  stats = MPI_REQUEST_NULL

  IF(DebugMode) PRINT *,myid,' NRXF pstart, np init: ',NRXF%pstrt, NRXF%NP

  !For new prox points received, put them into InvPartInfo and NRXF
  put_loc = NRXF%pstrt + NRXF%NP
  put_loc_init = put_loc

  DO i=0,ntasks-1
    IF(.NOT. PartIsNeighbour(i)) CYCLE

    DO j=1,NRXFPointEx(i) % rcount
      new_id = NRXFPointEx(i) % RecvIDs(j)

      IF(put_loc > SIZE(NRXF%A,2)) THEN
        IF(DebugMode) PRINT *,myid,' ExchangeProxPoints: resizing NRXF & UT: ',&
             SIZE(NRXF%A,2),SIZE(UT%A)/6.0, SIZE(UTM%A)/6.0
        CALL ResizePointData(NRXF,2.0_8,UT,UTM,.FALSE.,.FALSE.,.TRUE.)
      END IF

      NRXF%PartInfo(1,put_loc) = i
      NRXF%PartInfo(2,put_loc) = new_id

      InvPartInfo(i) % PCount = InvPartInfo(i) % PCount + 1
      cnt = InvPartInfo(i) % PCount

      IF(cnt > SIZE(InvPartInfo(i) % ProxIDs)) THEN
        CALL ExpandIntArray(InvPartInfo(i) % ProxIDs, fill_in=-1)
        CALL ExpandIntArray(InvPartInfo(i) % ProxLocs, fill_in=-1)
      END IF
      InvPartInfo(i) % ProxIDs(cnt) = new_id
      InvPartInfo(i) % ProxLocs(cnt) = put_loc
      NRXF%NP = NRXF%NP + 1
      put_loc = put_loc + 1
    END DO
  END DO

  IF(DebugMode) PRINT *,myid,' NRXF pstart, np putloc post: ',NRXF%pstrt, NRXF%NP,put_loc

  !Pass NRXF
  DO i=0,ntasks-1
    IF(.NOT. PartIsNeighbour(i)) CYCLE

    DO j=1,NRXFPointEx(i) % scount
      id = NRXFPointEx(i) % SendIDs(j)
      NRXFPointEx(i) % S(j*3 - 2) = NRXF%M(1,id)
      NRXFPointEx(i) % S(j*3 - 1) = NRXF%M(2,id)
      NRXFPointEx(i) % S(j*3 - 0) = NRXF%M(3,id)
    END DO

    CALL MPI_ISend(NRXFPointEx(i) % S, NRXFPointEx(i) % scount*3, MPI_DOUBLE_PRECISION, &
         i, 203, MPI_COMM_WORLD, stats(i*2), ierr)
    CALL MPI_IRecv(NRXFPointEx(i) % R, NRXFPointEx(i) % rcount*3, MPI_DOUBLE_PRECISION, &
         i, 203, MPI_COMM_WORLD, stats(i*2+1), ierr)
  END DO

  CALL MPI_Waitall(ntasks*4, stats, MPI_STATUSES_IGNORE, ierr)
  stats = MPI_REQUEST_NULL

  put_loc = put_loc_init
  DO i=0,ntasks-1
    IF(.NOT. PartIsNeighbour(i)) CYCLE

    DO j=1,NRXFPointEx(i) % rcount
      !new_id = NRXFPointEx(i) % RecvIDs(j) <- not needed, order is unchanged
      NRXF%A(:,put_loc) = NRXFPointEx(i) % R(j*3-2 : j*3)
      put_loc = put_loc + 1
    END DO
  END DO

  IF(DebugMode) PRINT *,myid,' NRXF pstart, np putloc ppost: ',NRXF%pstrt, NRXF%NP,put_loc

  !Pass UT & UTM
  DO i=0,ntasks-1
    IF(.NOT. PartIsNeighbour(i)) CYCLE

    PointEx(i) % S = 0.0
    PointEx(i) % R = 0.0
    DO j=1,PointEx(i) % scount
      id = PointEx(i) % SendIDs(j)
      PointEx(i) % S(j*12 - 11) = UT%M(6*id - 5)
      PointEx(i) % S(j*12 - 10) = UT%M(6*id - 4)
      PointEx(i) % S(j*12 - 9) = UT%M(6*id - 3)
      PointEx(i) % S(j*12 - 8) = UT%M(6*id - 2)
      PointEx(i) % S(j*12 - 7) = UT%M(6*id - 1)
      PointEx(i) % S(j*12 - 6) = UT%M(6*id - 0)
      PointEx(i) % S(j*12 - 5) = UTM%M(6*id - 5)
      PointEx(i) % S(j*12 - 4) = UTM%M(6*id - 4)
      PointEx(i) % S(j*12 - 3) = UTM%M(6*id - 3)
      PointEx(i) % S(j*12 - 2) = UTM%M(6*id - 2)
      PointEx(i) % S(j*12 - 1) = UTM%M(6*id - 1)
      PointEx(i) % S(j*12 - 0) = UTM%M(6*id - 0)
    END DO

    CALL MPI_ISend(PointEx(i) % S, PointEx(i) % scount*12, MPI_DOUBLE_PRECISION, &
         i, 203, MPI_COMM_WORLD, stats(i*2),ierr)
    CALL MPI_IRecv(PointEx(i) % R, PointEx(i) % rcount*12, MPI_DOUBLE_PRECISION, &
         i, 203, MPI_COMM_WORLD, stats(i*2+1),ierr)
    
  END DO

  CALL MPI_Waitall(ntasks*4, stats, MPI_STATUSES_IGNORE, ierr)
  stats = MPI_REQUEST_NULL

  IF(PrintTimes) THEN
    CALL CPU_TIME(T2)
    PRINT *,myid,' ExchangeProxPoints Part 1 time: ',T2-T1,' secs'
    CALL CPU_TIME(T1)
  END IF
  !Store UT & UTM
  ! PointEx(i) % rcount = InvPartInfo(i) % Pcount, and the values of
  ! InvPartInfo(i) % ProxIDs are all present in PointEx(i) % RecvIDs
  DO i=0,ntasks-1
    IF(.NOT. PartIsNeighbour(i)) CYCLE

    CALL sort_int2(InvPartInfo(i) % ProxIDs, InvPartInfo(i) % ProxLocs, InvPartInfo(i) % Pcount)

    !Deal with loss of particles
    !e.g. if we lost particles 7 & 25:
    ! pcount=7,ProxIDs:  5  7  9  11 17  20  25
    ! rcount=5,RecvIDs:  5  9  11 17 20
    IF(InvPartInfo(i) % PCount > PointEx(i) % rcount) THEN
      IF(DebugMode) PRINT *,myid,' debug, detecting loss of nearby points from part ',i,&
           ' pcount,rcount: ',InvPartInfo(i) % PCount, PointEx(i) % rcount
      rmcount = 0
      DO j=1,InvPartInfo(i) % Pcount
        loss = .FALSE.
        IF(j-rmcount > PointEx(i) % rcount) THEN !this catches loss of e.g. '25'
          loss = .TRUE.
        ELSE IF( InvPartInfo(i) % ProxIDs(j) /= PointEx(i) % RecvIDs(j-rmcount)) THEN !this catches '7'
          loss = .TRUE.
        END IF
        IF(loss) THEN
          rmcount = rmcount + 1
          loc = InvPartInfo(i) % ProxLocs(j)

          IF(DebugMode) PRINT *,myid,' debug lost: ',InvPartInfo(i) % ProxIDs(j),&
               ' from ',i,' loc: ',loc,NRXF%A(:,loc),NRXF%PartInfo(:,loc)

          NRXF%A(:,loc) = 0.0
          NRXF%PartInfo(:,loc) = -1
          !NRXF%NP = NRXF%NP - 1
          UT%A(6*loc-5: 6*loc) = 0.0
          UTM%A(6*loc-5: 6*loc) = 0.0
          InvPartInfo(i) % ProxIDs(j) = -1
        END IF
      END DO

      InvPartInfo(i) % PCount = InvPartInfo(i) % PCount - rmcount
      !Tidy up ProxLocs
      InvPartInfo(i) % ProxLocs(1:InvPartInfo(i) % PCount) = &
           PACK(InvPartInfo(i) % ProxLocs,InvPartInfo(i) % ProxIDs /= -1)
      InvPartInfo(i) % ProxLocs(InvPartInfo(i) % PCount+1:UBOUND(InvPartInfo(i) % ProxLocs,1)) = -1
      !And ProxIDs
      InvPartInfo(i) % ProxIDs(1:InvPartInfo(i) % PCount) = &
           PACK(InvPartInfo(i) % ProxIDs,InvPartInfo(i) % ProxIDs /= -1)
      InvPartInfo(i) % ProxIDs(InvPartInfo(i) % PCount+1:UBOUND(InvPartInfo(i) % ProxIDs,1)) = -1

    ELSE IF(InvPartInfo(i) % PCount < PointEx(i) % rcount) THEN
      CALL FatalError("ExchangeProxPoints - got more points than we ought")
    END IF

    DO j=1,PointEx(i) % rcount
      IF(PointEx(i) % RecvIDs(j) /= InvPartInfo(i) % ProxIDs(j)) THEN
        PRINT *,myid,' debug, j,recvid,proxid: ',j, PointEx(i) % RecvIDs(j), &
             InvPartInfo(i) % ProxIDs(j)
        PRINT *,myid,' debug, j,all recvid ,proxid: ',j, PointEx(i) % RecvIDs(:), &
             InvPartInfo(i) % ProxIDs(:)
        PRINT *,myid,' debug, j,SIZE recvid ,proxid: ',j, SIZE(PointEx(i) % RecvIDs), &
             SIZE(InvPartInfo(i) % ProxIDs), InvPartInfo(i) % Pcount, PointEx(i) % rcount
        CALL FatalError("ExchangeProxPoints: Incorrect sorting assumption")
      END IF
      put_loc = InvPartInfo(i) % ProxLocs(j)

      UT%A(6*put_loc - 5) = PointEx(i) % R(j*12 - 11)
      UT%A(6*put_loc - 4) = PointEx(i) % R(j*12 - 10)
      UT%A(6*put_loc - 3) = PointEx(i) % R(j*12 - 9)
      UT%A(6*put_loc - 2) = PointEx(i) % R(j*12 - 8)
      UT%A(6*put_loc - 1) = PointEx(i) % R(j*12 - 7)
      UT%A(6*put_loc - 0) = PointEx(i) % R(j*12 - 6)

      UTM%A(6*put_loc - 5) = PointEx(i) % R(j*12 - 5)
      UTM%A(6*put_loc - 4) = PointEx(i) % R(j*12 - 4)
      UTM%A(6*put_loc - 3) = PointEx(i) % R(j*12 - 3)
      UTM%A(6*put_loc - 2) = PointEx(i) % R(j*12 - 2)
      UTM%A(6*put_loc - 1) = PointEx(i) % R(j*12 - 1)
      UTM%A(6*put_loc - 0) = PointEx(i) % R(j*12 - 0)

      IF(DebugMode) THEN
        IF(j==3) THEN
          PRINT *,myid,' 3rd from ',i,' is ',PointEx(i) % RecvIDs(j),' put_loc: ',&
               put_loc,' UT%A: ',UT%A(6*put_loc-5 : 6*put_loc),' NRXF: ',NRXF%A(:,put_loc)
        END IF
      END IF

    END DO
  END DO

  IF(PrintTimes) THEN
    CALL CPU_TIME(T2)
    PRINT *,myid,' ExchangeProxPoints Part 2 time: ',T2-T1,' secs'

    CALL CPU_TIME(tend)
    PRINT *,myid,'Exchange Prox Points took: ',tend -  tstrt,' secs'
  END IF

  IF(DebugMode) PRINT *,myid,' successfully completed ExchangeProxPoints.'
  CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

END SUBROUTINE ExchangeProxPoints

!Pass EFS values for shared beams between partitions
!This is only required because the beam EFS values are randomly
!generated. The larger partition ID sends to the lower.
!Note that it would be sufficient (in current implementation)
!to send an MPI_LOGICAL 'broken' instead of DP EFS value, but 
!this way is flexible to future change in EFS
SUBROUTINE ExchangeEFS(NANS, NANPart, NRXF, InvPartInfo, EFS)


  INTEGER, ALLOCATABLE :: NANS(:,:), NANPart(:)
  REAL(KIND=dp), ALLOCATABLE :: EFS(:)
  TYPE(NRXF_t) :: NRXF
  !-----------------------------
  INTEGER :: i,j,k,neigh,counter,N1,N2,NTOT,ierr,recvtot
  INTEGER, ALLOCATABLE :: stats(:)
  REAL(KIND=dp) :: T1,T2
  LOGICAL :: Send
  TYPE(PointEx_t), ALLOCATABLE :: PointEx(:)
  TYPE(InvPartInfo_t) :: InvPartInfo(0:)

  ALLOCATE(stats(0:ntasks*2-1), PointEx(0:ntasks-1))
  stats = MPI_REQUEST_NULL

  IF(SIZE(EFS) /= SIZE(NANPart) .OR. SIZE(EFS) /= SIZE(NANS,2)) &
       CALL FatalError("Size mismatch in ExchangeEFS")
  NTOT = SIZE(NANPart)

  DO i=0,ntasks-1
    IF(InvPartInfo(i) % ccount == 0) CYCLE
    neigh = InvPartInfo(i) % NID ! == i...
    Send = (neigh < myid) 

    counter = 0
    DO j=1,NTOT
      IF(NANPart(j) /= neigh) CYCLE
      counter = counter + 1
    END DO

    IF(DebugMode) PRINT *,myid,' sharing ',counter,' EFS with ',neigh

    PointEx(i) % scount = counter
    PointEx(i) % rcount = counter

    IF(Send) THEN
      ALLOCATE(PointEx(i) % S(counter), &
           PointEx(i) % SendIDs(counter*2))

      counter = 0
      DO j=1,NTOT
        IF(NANPart(j) /= neigh) CYCLE
        counter = counter + 1
        PointEx(i) % S(counter) = EFS(j)
        PointEx(i) % SendIDs(counter*2-1) = NRXF%PartInfo(2,NANS(1,j))
        PointEx(i) % SendIDs(counter*2)   = NRXF%PartInfo(2,NANS(2,j))
      END DO
    ELSE
      ALLOCATE(PointEx(i) % R(counter), &
           PointEx(i) % RecvIDs(counter * 2))
    END IF

    IF(Send) THEN
      CALL MPI_ISend(PointEx(i) % S,counter,MPI_DOUBLE_PRECISION,neigh,&
           190,MPI_COMM_WORLD,stats(i*2),ierr)

      CALL MPI_ISend(PointEx(i) % SendIDs,counter*2,MPI_INTEGER,neigh,&
           191,MPI_COMM_WORLD,stats(i*2+1),ierr)

    ELSE
      CALL MPI_IRecv(PointEx(i) % R,counter,MPI_DOUBLE_PRECISION,neigh,&
           190,MPI_COMM_WORLD,stats(i*2),ierr)

      CALL MPI_IRecv(PointEx(i) % RecvIDs,counter*2,MPI_INTEGER,neigh,&
           191,MPI_COMM_WORLD,stats(i*2+1),ierr)

    END IF
  END DO

  IF(DebugMode) PRINT *,myid,' finished EFS send'
  CALL MPI_WaitAll(ntasks*2, stats, MPI_STATUSES_IGNORE, ierr)
  IF(DebugMode) PRINT *,myid,' finished EFS recv'

  CALL CPU_Time(T1)

  !First replace otherpart's N1 with our N1 node ID
  recvtot = 0
  DO i=0,ntasks-1
    IF(i<=myid) CYCLE
    recvtot = recvtot + PointEx(i) % rcount
    DO j=1,PointEx(i) % rcount
      N1 = PointEx(i) % RecvIDs(j*2)
      DO k=1,InvPartInfo(i) % CCount
        IF(InvPartInfo(i) % ConnIDs(k) == N1) THEN
          N1 = InvPartInfo(i) % ConnLocs(k)
          EXIT
        END IF
      END DO
      PointEx(i) % RecvIDs(j*2) = N1
    END DO
  END DO

  counter = 0
  DO k=1,NTOT
    IF(NANPart(k) <= myid) CYCLE
    neigh = NANPart(k)

    DO j=1,PointEx(neigh) % rcount
      N2 = PointEx(neigh) % RecvIDs(j*2-1) !our ID
      N1 = PointEx(neigh) % RecvIDs(j*2)   !otherpart ID
      IF(NANS(2,k) /= N2) CYCLE
      IF(NANS(1,k) /= N1) CYCLE

      EFS(k) = PointEx(neigh) % R(j)
      counter = counter + 1
      EXIT
    END DO
  END DO

  IF(counter /= recvtot) THEN
    PRINT *,myid,' counter recvtot: ',counter, recvtot
    PRINT *,myid,' OH NO! '; STOP
  END IF

  CALL CPU_Time(T2)

  IF(PrintTimes) PRINT *,myid,' EFS replacement took: ',T2-T1,' secs'
  CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

END SUBROUTINE ExchangeEFS

!Use octree search to find nodes which are nearby
!Direct contact is identified by circ
!TODO - take care of IsLost particles
SUBROUTINE FindNearbyParticles(NRXF, UT, NN, BBox,search_dist,ND,NDL)

  USE Octree

  TYPE(NRXF_t) :: NRXF
  TYPE(UT_t) :: UT
  INTEGER :: NN,ND
  INTEGER, ALLOCATABLE :: NDL(:,:)
  REAL(KIND=dp) :: search_dist, LNN, BBox(6)
  !-----------------------
  !octree stuff
  type(point_type), allocatable :: points(:)
  INTEGER i,j,cnt,npoints, num_ngb, totsize
  INTEGER, ALLOCATABLE :: seed(:), ngb_ids(:), point_loc(:)
  REAL(KIND=dp) :: x(3), dx(3),oct_bbox(2,3),eps,T1,T2
  REAL(KIND=dp) :: tstrt,tend

  REAL(KIND=dp) :: X1,X2,Y1,Y2,Z1,Z2
  INTEGER :: id
  CALL CPU_Time(tstrt)

  eps = 1.0
!  search_dist = 1.87 * SCL
  npoints = COUNT(NRXF%PartInfo(1,:) /= -1)
  totsize = SIZE(NRXF%PartInfo,2)

  IF(DebugMode) PRINT *, myid, 'debug nrxf count: ',nn,npoints,SIZE(NRXF%PartInfo,2)

  IF(.NOT. ALLOCATED(NDL)) ALLOCATE(NDL(2,npoints*12)) !guess size
  ALLOCATE(points(npoints),point_loc(npoints)) 

  point_loc = 0
  ND = 0
  NDL = 0

  !buffer by search_dist to ensure all points belonging to this partition
  !are contained, + eps to avoid floating point issues.
  DO i=1,3
    oct_bbox(1,i) = BBox(i*2-1) - search_dist - eps 
    oct_bbox(2,i) = BBox(i*2)   + search_dist + eps
  END DO

  CALL Octree_init(max_depth=10,max_num_point=6,bbox=oct_bbox)

  cnt = 0
  DO i=1,totsize
    IF(NRXF%PartInfo(1,i) == -1) CYCLE
    cnt = cnt+1
    point_loc(cnt) = i
    Points(cnt) % id = cnt
    Points(cnt) % x(1) = NRXF%A(1,i) + UT%A(6*I-5)
    Points(cnt) % x(2) = NRXF%A(2,i) + UT%A(6*I-4)
    Points(cnt) % x(3) = NRXF%A(3,i) + UT%A(6*I-3)
  END DO

  IF(PrintTimes) CALL CPU_TIME(T1)
  CALL Octree_build(Points)

  IF(PrintTimes) THEN
    CALL CPU_TIME(T2)
    PRINT *,myid,'FindNearbyParticles: Time to build octree: ',T2-T1
    CALL CPU_TIME(T1)
  END IF

  ALLOCATE(ngb_ids(1000))

  DO i=1,npoints

    IF(NRXF%PartInfo(1,point_loc(i)) /= myid) CYCLE

    num_ngb = 0
    ngb_ids = 0
    CALL Octree_search(points(i) % x, search_dist, num_ngb, ngb_ids)

    !Debug info
    IF(num_ngb > 50) THEN
      id = point_loc(i)
      X1 = NRXF%A(1,id) + UT%A(id*6 - 5)
      Y1 = NRXF%A(2,id) + UT%A(id*6 - 4)
      Z1 = NRXF%A(3,id) + UT%A(id*6 - 3)
      PRINT *,myid,'Node ',point_loc(i),' has ',num_ngb,' neighbours: ',X1,Y1,Z1
      DO j=1,num_ngb
        id = point_loc(ngb_ids(j))
        X2 = NRXF%A(1,id) + UT%A(id*6 - 5)
        Y2 = NRXF%A(2,id) + UT%A(id*6 - 4)
        Z2 = NRXF%A(3,id) + UT%A(id*6 - 3)
        PRINT *,myid,' Node ',i,' neigh ',j,id,'PartInfo: ',NRXF%PartInfo(:,id),' is :',X2,Y2,Z2
      END DO
    END IF

    DO j=1,num_ngb
      !Only save each pair once
      !This works w/ other partition points too because they are 
      !guaranteed to be higher up in the array than our own (>NN)
      IF(point_loc(ngb_ids(j)) <= point_loc(i)) CYCLE

      ND = ND + 1
      IF(ND > SIZE(NDL,2)) CALL ExpandIntArray(NDL)
      NDL(1,ND) = point_loc(i)
      NDL(2,ND) = point_loc(ngb_ids(j))
    END DO
  END DO

  IF(PrintTimes) THEN
    CALL CPU_TIME(T2)
    PRINT *,myid,'FindNearbyParticles: Time to search octree: ',T2-T1
  END IF

  CALL Octree_final()

  IF(PrintTimes) THEN
    CALL CPU_TIME(tend)
    PRINT *,myid,'Finding nearby particles took: ',tend - tstrt,' secs'
  END IF

END SUBROUTINE FindNearbyParticles


!Use octree search to find nodes which are in contact
!TODO - take care of IsLost particles
SUBROUTINE FindBeams(xo, ip, SCL, CN, nbeams, searchdist_in)
  
  USE Octree

  REAL(KIND=dp), ALLOCATABLE :: xo(:,:)
  REAL(KIND=dp) :: SCL
  REAL(KIND=dp), OPTIONAL :: searchdist_in
  INTEGER :: ip,nbeams,max_neigh
  TYPE(Conn_t), ALLOCATABLE :: CN(:)
  !-----------------------
  !octree stuff
  type(point_type), allocatable :: points(:)
  INTEGER i,j,k, num_ngb,counter
  integer, allocatable :: seed(:), ngb_ids(:)
  REAL(KIND=dp) :: x(3), dx(3),oct_bbox(2,3),eps
  REAL(KIND=dp) :: dist, max_dist, min_dist, searchdist

  IF(PRESENT(searchdist_in)) THEN
    !If a search distance is given, we can't assume there'll only
    !be 12 neighbours - go for 30
    searchdist = searchdist_in
    max_neigh = 30
  ELSE 
    !Default usage - we are looking for particles which share a beam
    !FCC lattice = 12 neighbours
    searchdist = SCL * (1.6**0.5)
    max_neigh = 12
  END IF

  ALLOCATE(points(ip),&
       CN(ip))

  DO i=1,ip
    ALLOCATE(CN(i) % Conn(max_neigh))
    CN(i) % NCN = 0
    CN(i) % Conn = 0
    CN(i) % ID = i
  END DO

  eps = 1.0

  nbeams = 0

  DO i=1,3
    oct_bbox(1,i) = MINVAL(xo(i,1:ip)) - eps !buffer bbox to ensure all points contained
    oct_bbox(2,i) = MAXVAL(xo(i,1:ip)) + eps
  END DO

  !TODO - determine optimal depth - SCL vs bbox length?
  !  mdepth 20, m points 6 was determined for a particular case
  !   it doubled the speed w.r.t. 20,3 (i.e. max_num_points is important)
  !   we *expect* 12 neighbours for every node... 
  CALL Octree_init(max_depth=10,max_num_point=6,bbox=oct_bbox)

  DO i=1,ip
    Points(i) % id = i
    Points(i) % x(1) = xo(1,i)
    Points(i) % x(2) = xo(2,i)
    Points(i) % x(3) = xo(3,i)
  END DO

  CALL Octree_build(Points)

  ALLOCATE(ngb_ids(100))

  DO i=1,ip
    num_ngb = 0
    ngb_ids = 0
    CALL Octree_search(points(i) % x, searchdist, num_ngb, ngb_ids)
    counter = 0

    DO j=1,num_ngb
      IF(ngb_ids(j) == i) CYCLE
      counter = counter+1
      IF(counter > max_neigh) THEN 
        CALL Warn("FindBeams found too many neighbours, ignoring...")
        EXIT
      END IF

      CN(i) % NCN = CN(i) % NCN + 1
      CN(i) % Conn(counter) = ngb_ids(j)
      ! PRINT *,i, j, ' ngb_id ',ngb_ids(j)
      IF(ngb_ids(j) > i) nbeams = nbeams+1
    END DO
    IF(DebugMode .AND. myid==0 .AND. MOD(i,1000)==0)&
         PRINT *,myid,' node ',i,' nobeams: ',num_ngb
  END DO

  IF(DebugMode) PRINT *,myid,' sum(ncn) ',SUM(CN%NCN),' nbeams*2 ',nbeams*2
CALL Octree_final()

END SUBROUTINE FindBeams

!Use octree search to find nearest N particles to each particle
!TODO - take care of IsLost particles
SUBROUTINE FindNNearest(xo, ip, nfind, SCL, CN)
  
  USE Octree

  REAL(KIND=dp), ALLOCATABLE :: xo(:,:)
  REAL(KIND=dp) :: SCL
  INTEGER :: ip,nfind
  TYPE(Conn_t), ALLOCATABLE :: CN(:)
  !-----------------------
  !octree stuff
  type(point_type), allocatable :: points(:)
  INTEGER i,j,k, num_ngb,counter
  INTEGER, ALLOCATABLE :: seed(:), ngb_ids(:)
  REAL(KIND=dp), ALLOCATABLE :: ngb_dists(:)
  REAL(KIND=dp) :: x(3), dx(3),oct_bbox(2,3),eps
  REAL(KIND=dp) :: dist, max_dist, min_dist, searchdist

  IF(nfind > ip) CALL FatalError("Requested too many nearest particles for particle count!")

  eps = 1.0

  ALLOCATE(points(ip),&
       CN(ip))

  DO i=1,ip
    ALLOCATE(CN(i) % Conn(nfind))
    CN(i) % NCN = 0
    CN(i) % Conn = 0
    CN(i) % ID = i
  END DO

  !Get the octree bounding box
  DO i=1,3
    oct_bbox(1,i) = MINVAL(xo(i,1:ip)) - eps !buffer bbox to ensure all points contained
    oct_bbox(2,i) = MAXVAL(xo(i,1:ip)) + eps
  END DO

  CALL Octree_init(max_depth=10,max_num_point=6,bbox=oct_bbox)

  DO i=1,ip
    Points(i) % id = i
    Points(i) % x(1) = xo(1,i)
    Points(i) % x(2) = xo(2,i)
    Points(i) % x(3) = xo(3,i)
  END DO

  CALL Octree_build(Points)

  ALLOCATE(ngb_ids(MAX(10000,10*nfind)), ngb_dists(MAX(10000,10*nfind)))

  !Need to iterate over progressively larger search distances, only
  !searching for those which don't already have sufficient neighbours
  !We can begin with those which are connected by beams

  !Start with beams
  searchdist = SCL * (1.6**0.5)

  DO WHILE(.TRUE.) 
    DO i=1,ip

      !Already got sufficient neighbours for this node
      IF(CN(i) % NCN == nfind) CYCLE

      num_ngb = 0
      ngb_ids = 0
      ngb_dists = 0.0
      CALL Octree_search(points(i) % x, searchdist, num_ngb, ngb_ids, ngb_dists=ngb_dists)

      !Didn't find enough this time (including self) - need to cycle anyway 
      !so don't bother putting neighbours
      IF(num_ngb <= nfind) CYCLE

      CALL sort_real2(ngb_dists, ngb_ids, num_ngb)

      counter = 0
      DO j=1,num_ngb
        IF(ngb_ids(j) == i) CYCLE !should be the first sorted
        counter = counter+1

        CN(i) % Conn(counter) = ngb_ids(j)
        IF(counter == nfind) EXIT
      END DO
      CN(i) % NCN = counter

    END DO
    
    !quit if we're done
    IF(ALL(CN%NCN == nfind)) EXIT

    !else double search dist and try again
    searchdist = searchdist * 1.1
  END DO


CALL Octree_final()

END SUBROUTINE FindNNearest

!Use octree search to find connections between separate clusters of particles
!TODO - take care of IsLost particles
SUBROUTINE FindClusterConns(xo, CN, Cluster, nfind, SCL)
  
  USE Octree

  REAL(KIND=dp) :: xo(:,:)
  REAL(KIND=dp) :: SCL
  TYPE(Conn_t), ALLOCATABLE :: CN(:)
  INTEGER :: Cluster(:)
  INTEGER :: nfind
  !-----------------------
  !octree stuff
  type(point_type), allocatable :: points(:)
  INTEGER i,j,k,l,ip,num_ngb,counter,np,id, NClust, conn_seek,N1,N2
  INTEGER, ALLOCATABLE :: ngb_ids(:), mins(:), min_conns(:,:,:), &
       ClusterMap(:)
  REAL(KIND=dp), ALLOCATABLE :: ngb_dists(:)
  REAL(KIND=dp) :: oct_bbox(2,3),eps
  REAL(KIND=dp) :: dist, max_dist, min_dist, searchdist
  LOGICAL, ALLOCATABLE :: ClusterActive(:)
  TYPE(Conn_t), ALLOCATABLE :: ClustCN(:)

  eps = 1.0
  conn_seek = 5
  ip = SIZE(Cluster)
  IF(nfind > ip) CALL FatalError("Requested too many nearest particles for particle count!")

  ALLOCATE(points(ip))

  NClust = MAXVAL(Cluster)

  !Initialize and build the octree
  !----------------------------

  !Get the octree bounding box
  DO i=1,3
    oct_bbox(1,i) = MINVAL(xo(i,1:ip)) - eps !buffer bbox to ensure all points contained
    oct_bbox(2,i) = MAXVAL(xo(i,1:ip)) + eps
  END DO

  CALL Octree_init(max_depth=10,max_num_point=6,bbox=oct_bbox)

  DO i=1,ip
    Points(i) % id = i
    Points(i) % tag = Cluster(i)
    Points(i) % x(1) = xo(1,i)
    Points(i) % x(2) = xo(2,i)
    Points(i) % x(3) = xo(3,i)
  END DO

  CALL Octree_build(Points)

  ALLOCATE(ngb_ids(MAX(1000,10*nfind)), ngb_dists(MAX(1000,10*nfind)), ClusterActive(NClust))
  ClusterActive = .TRUE.


  !Find connections and redefine clusters iteratively until only 1 left
  !--------------------------
  DO WHILE(.TRUE.) 
    np = COUNT(Cluster /= 1)

    DO i=1,NClust
      IF(.NOT. ANY(Cluster == i)) ClusterActive(i) = .FALSE.
    END DO

    ALLOCATE(ClustCN(np))

    counter = 0
    DO i=1,ip
      !Set the tag for each point for the octree search
      Points(i) % tag = Cluster(i)

      IF(Cluster(i) == 1) CYCLE
      counter = counter + 1
      ClustCN(counter) % ID = i
      ClustCN(counter)%NCN = 0
      ALLOCATE(ClustCN(counter) % Conn(nfind),&
           ClustCN(counter) % Dists(nfind))
    END DO

    !Need to iterate over progressively larger search distances, only
    !searching for those which don't already have sufficient neighbours
    !Start with 5 * SCL
    searchdist = SCL * 5.0

    DO WHILE(.TRUE.) 
      DO i=1,np

        id = ClustCN(i) % ID

        !Already got sufficient neighbours for this node
        IF(ClustCN(i) % NCN == nfind) CYCLE

        num_ngb = 0
        ngb_ids = 0
        ngb_dists = 0.0
        CALL Octree_search(points(id) % x, searchdist, num_ngb, ngb_ids, ngb_dists=ngb_dists, &
             tag=Cluster(id), same_tag=.FALSE.)

        !Didn't find enough this time (including self) - need to cycle anyway 
        !so don't bother putting neighbours
        IF(num_ngb <= nfind) CYCLE

        CALL sort_real2(ngb_dists, ngb_ids, num_ngb)

        counter = 0
        DO j=1,num_ngb
          IF(ngb_ids(j) == id) CYCLE !should be the first sorted
          counter = counter+1

          ClustCN(i) % Conn(counter) = ngb_ids(j)
          ClustCN(i) % Dists(counter) = ngb_dists(j)
          IF(counter == nfind) EXIT
        END DO
        ClustCN(i) % NCN = counter

      END DO

      !quit if we're done
      IF(ALL(ClustCN%NCN == nfind)) EXIT

      !else double search dist and try again
      searchdist = searchdist * 2.0
      IF(DebugMode) PRINT *,myid,' octree searching dist: ',searchdist
    END DO

    !For each cluster, find the N closest connections to another cluster
    ALLOCATE(mins(conn_seek), min_conns(NClust,conn_seek,2))
    min_conns = 0

    DO i=2,NClust !don't do cluster 1 (glacier)
      IF(.NOT. ClusterActive(i)) CYCLE
      mins = HUGE(mins)
      DO j=1,np
        IF(Cluster(ClustCN(j) % ID) /= i) CYCLE
        DO k=1,ClustCN(j) % NCN
          !Cycle stored mins comparing, store if reqd
          !Need another array to hold the IDs
          DO l=1,conn_seek
            IF(ClustCN(j) % Dists(k) < mins(l)) THEN
              mins(l) = ClustCN(j) % Dists(k)
              min_conns(i,l,1) = ClustCN(j) % ID
              min_conns(i,l,2) = ClustCN(j) % Conn(k)
              EXIT
            END IF
          END DO
        END DO
      END DO
      IF(DebugMode)  PRINT *,'Cluster ',i,' min dists: ',mins

      DO j=1,conn_seek
        IF(DebugMode)  PRINT *,'Cluster ',i,' conn ',j,' : ',min_conns(i,j,:), Cluster(min_conns(i,j,:))
      END DO
    END DO

    !Add connections to CN, NCN
    DO i=2, NClust
      IF(.NOT. ClusterActive(i)) CYCLE

      DO j=1,conn_seek
        N1 = min_conns(i,j,1)
        N2 = min_conns(i,j,2)

        IF(.NOT. ANY(CN(N1) % Conn == N2)) THEN
          CN(N1) % NCN = CN(N1) % NCN + 1
          IF(CN(N1) % NCN > SIZE(CN(N1) % Conn)) CALL ExpandIntArray(CN(N1) % Conn)
          CN(N1) % Conn(CN(N1) % NCN) = N2
        END IF

        IF(.NOT. ANY(CN(N2) % Conn == N1)) THEN
          CN(N2) % NCN = CN(N2) % NCN + 1
          IF(CN(N2) % NCN > SIZE(CN(N2) % Conn)) CALL ExpandIntArray(CN(N2) % Conn)
          CN(N2) % Conn(CN(N2) % NCN) = N1
        END IF

      END DO
    END DO

    !Update clusters
    ALLOCATE(ClusterMap(NClust))
    DO i=1,NClust
      ClusterMap(i) = i
    END DO

    !Cycle smallest to largest cluster, joining
    DO i=NClust,2,-1
      IF(.NOT. ClusterActive(i)) CYCLE
      DO j=1,conn_seek
        ClusterMap(i) = MIN(ClusterMap(Cluster(min_conns(i,j,2))),ClusterMap(i))
      END DO
    END DO

    DO i=1,ip
      Cluster(i) = ClusterMap(Cluster(i))
    END DO

    DEALLOCATE(ClustCN, mins, min_conns, ClusterMap)

    IF(ALL(Cluster == 1)) EXIT
  END DO

  CALL Octree_final()

END SUBROUTINE FindClusterConns


!Checks that connectivity info between particles
!goes both ways (i.e. a proper dual-graph)
SUBROUTINE EnsureDualGraph(CN)
  TYPE(Conn_t), ALLOCATABLE :: CN(:)
  !----------------------------------
  INTEGER :: i,j,n,NN,neigh

  NN = SIZE(CN)
  n = MAXVAL(CN%NCN)

  DO i=1,NN
    DO j=1,CN(i)%NCN
      neigh = CN(i) % Conn(j)
      IF(.NOT. ANY(CN(neigh) % Conn == i)) THEN
        CN(neigh) % NCN = CN(neigh) % NCN + 1
        IF(CN(neigh) % NCN > SIZE(CN(neigh) % Conn)) CALL ExpandIntArray(CN(neigh) % Conn)
        CN(neigh) % Conn(CN(neigh) % NCN) = i
      END IF
    END DO
  END DO

END SUBROUTINE EnsureDualGraph

SUBROUTINE FindCollisions(SI,ND,NN,NRXF,UT,FRX,FRY,FRZ, &
     T,IS,WE,EFC,FXF,FXC,NDL,LNN)

  USE TypeDefs
  USE Utils

  IMPLICIT NONE
  REAL(KIND=dp) ::  X1,X2,Y1,Y2,Z1,Z2
  REAL(KIND=dp) ::  T1,T2
  REAL(KIND=dp), ALLOCATABLE :: EFC(:)
  REAL(KIND=dp) ::  SX,SY,SZ,SUM,T,WE(:),L0
  REAL(KIND=dp) ::  DDEL,DWE,OWE,ESUM,LNN
  REAL(KIND=dp) ::  LS,LS2,DEL,SCL,ViscForce
  INTEGER ierr,FXC,ND
  INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),comm
  INTEGER, ALLOCATABLE :: FXF(:,:),NDL(:,:)
  REAL(KIND=dp) ::  RC,RCX,RCY,RCZ,FRX(NN),FRY(NN),FRZ(NN)
  INTEGER NTOT,I,N1,N2,IS,NN
  TYPE(UT_t) :: UT
  TYPE(NRXF_t) :: NRXF
  LOGICAL :: own(2)
  TYPE(SimInfo_t) :: SI

  IF(PrintTimes) CALL CPU_TIME(T1)

  SCL = SI%SCL
  ViscForce = SI%ViscForce

  IF(.NOT. ALLOCATED(FXF)) ALLOCATE(FXF(2,NN*12))
  FXF = 0 !particle proximity info
  
  FRX = 0.0
  FRY=0.0
  FRZ=0.0
  WE=0.0
  FXC = 0
  FXF = 0

  DO I=1,ND
    N1=NDL(1,I)
    N2=NDL(2,I)

    !Particle may no longer be nearby
    IF(NRXF%PartInfo(1,N1) == -1 .OR. NRXF%PartInfo(1,N2) == -1) CYCLE

    own(1) = N1 <= NN
    own(2) = N2 <= NN

    IF(.NOT. ANY(own)) CYCLE !don't own either particle

    X1=NRXF%A(1,N1)+UT%A(6*N1-5)
    Y1=NRXF%A(2,N1)+UT%A(6*N1-4)
    Z1=NRXF%A(3,N1)+UT%A(6*N1-3)
    X2=NRXF%A(1,N2)+UT%A(6*N2-5)
    Y2=NRXF%A(2,N2)+UT%A(6*N2-4)
    Z2=NRXF%A(3,N2)+UT%A(6*N2-3)
    IF (ABS(X1-X2).LE.LNN.AND.ABS(Y1-Y2).LE.LNN.AND.ABS(Z1-Z2).LE.LNN) THEN
      RC=SQRT((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
      RCX=(X1-X2)/RC
      RCY=(Y1-Y2)/RC
      RCZ=(Z1-Z2)/RC

      !TODO - when this was parallel, FRX,Y,Z were only saved for *our* nodes (not other parts)
      !       does this matter??

      IF(FXC+1 > SIZE(FXF,2)) CALL ExpandIntArray(FXF)

      IF (RC.LT.LNN) THEN

        IF(own(1)) THEN !One of our own particles
          FRX(N1)=FRX(N1)+EFC(N1)*(LNN-RC)**1.5*RCX
          FRY(N1)=FRY(N1)+EFC(N1)*(LNN-RC)**1.5*RCY
          FRZ(N1)=FRZ(N1)+EFC(N1)*(LNN-RC)**1.5*RCZ
        END IF
        IF(own(2)) THEN !One of our own particles
          FRX(N2)=FRX(N2)-EFC(N2)*(LNN-RC)**1.5*RCX
          FRY(N2)=FRY(N2)-EFC(N2)*(LNN-RC)**1.5*RCY
          FRZ(N2)=FRZ(N2)-EFC(N2)*(LNN-RC)**1.5*RCZ
        END IF

        IF(own(2)) THEN
          WE(N2)=WE(N2)+0.4*EFC(N2)*(LNN-RC)**2.5
        ELSE
          WE(N1)=WE(N1)+0.4*EFC(N1)*(LNN-RC)**2.5
        END IF

        FXC=FXC+1
        FXF(1,FXC)=N1 
        FXF(2,FXC)=N2 
      ENDIF

      !if almost touching, add forces but don't register interaction
      !NOTE - TODO - is this an attractive force??
      IF (RC.GT.LNN.AND.RC.LT.LNN+SI%ViscDist*SCL) THEN
        IF(own(1)) THEN
          FRX(N1)=FRX(N1)+SCL**2.0*ViscForce*(LNN-RC)*RCX
          FRY(N1)=FRY(N1)+SCL**2.0*ViscForce*(LNN-RC)*RCY
          FRZ(N1)=FRZ(N1)+SCL**2.0*ViscForce*(LNN-RC)*RCZ
        END IF
        IF(own(2)) THEN
          FRX(N2)=FRX(N2)-SCL**2.0*ViscForce*(LNN-RC)*RCX
          FRY(N2)=FRY(N2)-SCL**2.0*ViscForce*(LNN-RC)*RCY
          FRZ(N2)=FRZ(N2)-SCL**2.0*ViscForce*(LNN-RC)*RCZ
        END IF
        IF(own(2)) THEN
          WE(N2)=WE(N2)+SCL**2.0*0.5*ViscForce*(LNN-RC)**2.0
        ELSE
          WE(N1)=WE(N1)+SCL**2.0*0.5*ViscForce*(LNN-RC)**2.0
        END IF
      ENDIF

    ENDIF
  ENDDO

  IF(PrintTimes) THEN
    CALL CPU_TIME(T2)
    PRINT *,myid,' Finding collisions took: ',T2-T1,' secs'
  END IF

  RETURN

END SUBROUTINE FindCollisions


FUNCTION PInBBox(i,NRXF, UT, BBox, Buff) RESULT(InBB)
  INTEGER :: i
  TYPE(NRXF_t) :: NRXF
  TYPE(UT_t) :: UT
  REAL(KIND=dp) :: BBox(6)
  REAL(KIND=dp), OPTIONAL :: Buff
  LOGICAL :: InBB
  !------------------------
  REAL(KIND=dp) :: X,Y,Z, Buffer

  IF(PRESENT(Buff)) THEN
    Buffer = Buff
  ELSE
    Buffer = 200.0
  END IF

  X = NRXF%A(1,i) + UT%A(6*i - 5)
  Y = NRXF%A(2,i) + UT%A(6*i - 4)
  Z = NRXF%A(3,i) + UT%A(6*i - 3)
  
  InBB = (X > (BBox(1)-Buffer) .AND. X < (BBox(2)+Buffer) .AND. &
       Y > (BBox(3)-Buffer) .AND. Y < (BBox(4)+Buffer) .AND. &
       Z > (BBox(5)-Buffer) .AND. Z < (BBox(6)+Buffer))
END FUNCTION PInBBox

SUBROUTINE ConnectRegions(xo, CN, SCL) 
  TYPE(Conn_t), ALLOCATABLE :: CN(:)
  REAL(KIND=dp) :: xo(:,:)
  REAL(KIND=dp) :: SCL

  !--------------------------
  INTEGER, ALLOCATABLE :: Cluster(:)
  INTEGER :: NN

  !Identify disconnected regions
  Cluster = GetClusters(CN)

  IF(ALL(Cluster == 1)) RETURN

  !Connect clusters with new bonds
  CALL FindClusterConns(xo, CN, Cluster, 5, SCL)

END SUBROUTINE ConnectRegions

!Function takes neighbourhood info in CN/NCN format
!and finds clusters of connected particles.
FUNCTION GetClusters(CN) RESULT(Cluster)
  TYPE(Conn_t), ALLOCATABLE :: CN(:)
  !--------------------------
  INTEGER, ALLOCATABLE :: Cluster(:),WorkArr(:),ClusterSizes(:),&
       ClusterIDs(:),ClusterMap(:)
  INTEGER :: i,NTOT,N1,N2,NClust, Starticle,NN

  NN = SIZE(CN)

  ALLOCATE(Cluster(NN))
  Cluster = 0
  Cluster(1) = 1
  NClust = 1
  Starticle = 1

  DO WHILE(.TRUE.)

    !Recursively mark all connected particles
    CALL MarkNeighbours(Cluster, Starticle, CN)

    !If any unmarked, start a new cluster
    IF(ANY(Cluster == 0)) THEN
      NClust = NClust + 1
      !Find a particle to mark w/ new cluster
      Starticle = MINLOC(Cluster,1)
      Cluster(Starticle) = NClust
    ELSE
      !All marked, exit the loop
      EXIT
    END IF

  END DO

  !Reorder clusters by size (largest = 1 - i.e. the largely intact glacier)
  ALLOCATE(ClusterSizes(NClust),ClusterIDs(NClust),ClusterMap(NClust))
  DO i=1,NClust
    ClusterSizes(i) = COUNT(Cluster==i)
    ClusterIDs(i) = i
  END DO

  CALL sort_int2(ClusterSizes,ClusterIDs,NClust)
  
  DO i=1,NClust
    ClusterMap(ClusterIDs(NClust-i+1)) = i
  END DO

  DO i=1,NN
    Cluster(i) = ClusterMap(Cluster(i))
  END DO

  CONTAINS

  RECURSIVE SUBROUTINE MarkNeighbours(Cluster, Starticle, CN)
    TYPE(Conn_t), ALLOCATABLE :: CN(:)
    INTEGER :: Cluster(:), Starticle
    !-----------------------
    INTEGER :: i

    DO i = 1,CN(Starticle) % NCN
      IF(Cluster(CN(Starticle) % Conn(i)) == NClust) CYCLE !already found

      IF(Cluster(CN(Starticle) % Conn(i)) /= 0) THEN
        CALL FatalError("GetClusters: Incomplete Dual Graph!")
      END IF

      Cluster(CN(Starticle) % Conn(i)) = NClust
      CALL MarkNeighbours(Cluster, CN(Starticle) % Conn(i), CN)
    END DO

  END SUBROUTINE MarkNeighbours

END FUNCTION GetClusters

END MODULE Lattice
