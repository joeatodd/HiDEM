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
SUBROUTINE FIBG3(NN,NTOT,NANS,NRXF,NANPart,particles_G,NCN,CN,CNPart,InvPartInfo,neighcount,l,&
     wrkdir,geomfile,SCL,grid,melta,wl,UC,StrictDomain)

  IMPLICIT NONE
  INCLUDE 'na90.dat'

!Real*8,ALLOCATABLE :: surf(:,:),bed(:,:)
Real*8 surf(-100:2000,-100:2000),bed(-100:2000,-100:2000),melt(-100:2000,-100:2000)
Real*8 :: x,y,s1,b1,b2,u1,grid,m1,melta,wl,UC,z1
Real*8 :: box,b,SCL
INTEGER :: l,NN,i,j
INTEGER :: N1,N2,xk,yk,neighcount,NTOT
INTEGER, ALLOCATABLE :: NCN(:),CN(:,:),CNPart(:,:), particles_G(:),NANS(:,:),NANPart(:)
CHARACTER(LEN=256) :: wrkdir,geomfile
LOGICAL :: StrictDomain
!TYPE(NTOT_t) :: NTOT
TYPE(NRXF2_t) :: NRXF
TYPE(InvPartInfo_t), ALLOCATABLE :: InvPartInfo(:)
!Open(300,file='mass.dat',STATUS='OLD')

 12    FORMAT(2I8,' ',2F14.7)
 13    FORMAT(4F14.7)
 14    FORMAT(3F14.7)

surf=-100.0
bed=1000.0

!TODO pass through file once to check extent, then allocate bed, surf, melt
! then reread
OPEN(400,file=TRIM(geomfile),STATUS='OLD')
READ(400,*) N2
DO I=1,N2
  READ(400,*) x,y,s1,b1,b2,z1
  xk=INT(x/grid)
  yk=INT(y/grid)
  IF (xk.GE.-100.AND.yk.GE.-100) bed(xk,yk)=b1
  IF (xk.GE.-100.AND.yk.GE.-100) surf(xk,yk)=s1
  IF (xk.GE.-100.AND.yk.GE.-100) melt(xk,yk)=melta*0.0
ENDDO
CLOSE(400)

!box is never actualy used...
!b is used, which is box/l, so L never actually enters into this, so it's only
!used for the number of vertical layers
box=2.0d0**(2.0d0/3.0d0)*REAL(l,8) ! box size equal to fcc ground state
CALL Initializefcc(NN,NTOT,NANS,NRXF,NANPart,particles_G, NCN, CN, CNPart, InvPartInfo,&
     neighcount,box,l,wrkdir,SCL,surf,bed,melt,grid,wl,UC,StrictDomain)

CLOSE(400)

END SUBROUTINE FIBG3

!---------------------------------------------------------------!


SUBROUTINE Initializefcc(NN,NTOT,NANS,NRXF,NANPart,particles_G, NCN, CN, CNPart, &
     InvPartInfo,neighcount,box,l,wrkdir,SCL,surf,bed,melt,grid,wl,UC,StrictDomain)

  !USE INOUT

  IMPLICIT NONE

  INCLUDE 'na90.dat'
  INCLUDE 'mpif.h'

Real*8 surf(-100:2000,-100:2000),bed(-100:2000,-100:2000),melt(-100:2000,-100:2000)
Real*8 b,x0(3,4),box,SCL
REAL*8 gridminx, gridmaxx, gridminy, gridmaxy,T1,T2
REAL*8 z,x,y,sint,bint,mint,grid,wl,lc,UC,UCV, efficiency
REAL*8 :: X1,X2,Y1,Y2,Z1,Z2,RC,rc_min,rc_max
REAL*8, ALLOCATABLE :: xo(:,:),work_arr(:,:)

INTEGER i,j,k,n,K1,k2,l,ip,NN,nb,xk,yk,nx,ny,nbeams,ierr,pown
INTEGER neighcount,NTOT
INTEGER, ALLOCATABLE :: NCN_All(:), CN_All(:,:),NCN(:),CN(:,:),CNPart(:,:),&
     particles_G(:),particles_L(:),neighparts(:), NANS(:,:),NANPart(:)
CHARACTER(LEN=256) :: wrkdir
LOGICAL :: StrictDomain
LOGICAL, ALLOCATABLE :: SharedNode(:)

!metis stuff
INTEGER :: objval,counter
INTEGER, ALLOCATABLE :: metis_options(:), particlepart(:),counters(:)
INTEGER, POINTER :: vwgt=>NULL(),vsize=>NULL(),adjwgt=>NULL(),xadj(:),adjncy(:),countptr
REAL*8, POINTER :: tpwgts=>NULL(),ubvec=>NULL()

!TYPE(NTOT_t) :: NTOT
TYPE(NRXF2_t), TARGET :: NRXF
TYPE(NRXF_t) :: NRXFold
TYPE(InvPartInfo_t), ALLOCATABLE, TARGET :: InvPartInfo(:)

11    FORMAT(2I8,' ',2F14.7)
13    FORMAT(4F14.7)

ALLOCATE(xo(3,NOMA))

b=SCL*box/REAL(l,8)  ! the size of the unit cell is the box length divided by l
x0(:,:)=b/2.0d0; x0(:,1)=0.0d0; x0(3,2)=0.0d0; x0(2,3)=0.0d0; x0(1,4)=0.0d0 

!Use bed and surf to determine the extent of the domain
!This is not perfect, just checks that there is SCL dist between surf and bed+melt

gridminx = HUGE(gridminx); gridminy = HUGE(gridminy)
gridmaxx = -HUGE(gridmaxx); gridmaxy = -HUGE(gridmaxy)
DO i=LBOUND(surf,1), UBOUND(surf,1)

  DO j=LBOUND(surf,2), UBOUND(surf,2)
    IF(.NOT. surf(i,j) > (bed(i,j) + melt(i,j))) CYCLE
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

!nx (ny) is the number of boxes in the x (y) direction
ip=0
DO i=1,nx !x step
      DO j=1,ny !y step
	DO k=-2,l !vertical layer
	  Do k1=1,4

             !These x,y coords are divided by grid just for the interp
             x=(x0(1,k1) + REAL(i-1)*b)/grid
             y=(x0(2,k1) + REAL(j-1)*b)/grid
             xk=INT((x0(1,k1) + REAL(i-1)*b)/grid)
             yk=INT((x0(2,k1) + REAL(j-1)*b)/grid)

             !just beyond edge of geom def
             IF(ANY(surf(xk:xk+1,yk:yk+1) == bed(xk:xk+1,yk:yk+1)) .AND. StrictDomain) CYCLE

             Call BIPINT(x-xk,y-yk,bed(xk,yk),bed(xk,yk+1),bed(xk+1,yk),bed(xk+1,yk+1),bint)
             Call BIPINT(x-xk,y-yk,surf(xk,yk),surf(xk,yk+1),surf(xk+1,yk),surf(xk+1,yk+1),sint)
             Call BIPINT(x-xk,y-yk,melt(xk,yk),melt(xk,yk+1),melt(xk+1,yk),melt(xk+1,yk+1),mint)
!             bint=bed(xk,yk)+(x-xk)*(bed(xk+1,yk)-bed(xk,yk))+(y-yk)*(bed(xk,yk+1)-bed(xk,yk))
!             sint=surf(xk,yk)+(x-xk)*(surf(xk+1,yk)-surf(xk,yk))+(y-yk)*(surf(xk,yk+1)-surf(xk,yk))

              ! these are the actual point coords
              z=x0(3,k1) + REAL(k-1)*b
              y=x0(2,k1) + REAL(j-1)*b
              x=x0(1,k1) + REAL(i-1)*b

!             write(1510+myid,13) 40.0*x,40.0*y,bint,sint
!             If (bed(xk,yk).ne.0.0.and.bed(xk,yk+1).ne.0.0.and.bed(xk+1,yk).ne.0.0&
              ! .AND.bed(xk+1,yk+1).NE.0.0) THEN

              !TODO - unhardcode this
              ! If (z.ge.bint+mint.and.z.le.sint.and.((sint-bint).gt.4.0*SCL.or.&
              ! (ABS(z-wl).LT.4.0*SCL.AND.bint.LT.wl))) THEN
             IF (z.GE.bint+mint .AND. z.LT.sint .AND. (sint-(bint+mint)).GT.SCL) THEN

             ! undercut shape functions
             ! lc=4420.0+1.5e-04*(x-3300.0)**2+0.42*exp((x-3700.0)/2.0e+02)
             ! UCV=lc-1500.0*exp(-(x-3500.0)**2/50000.0)
             ! If (y.lt.lc-UC.or.y.gt.lc.or.z.gt.bint+3.0*SQRT(y-(lc-UC)).or.z.ge.WL-20.0) then
             ! If (y.lt.lc-UC.or.z.gt.bint+3.0*SQRT(y-(lc-UC)).or.z.ge.WL-40.0) then
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
               CALL Warn("Doubling size of xo")
             END IF

	     xo(1,ip) = x0(1,k1) + REAL(i-1)*b	
	     xo(2,ip) = x0(2,k1) + REAL(j-1)*b
	     xo(3,ip) = x0(3,k1) + REAL(k-1)*b
             ! EndIF
             ! EndIF
             EndIF
!             EndIF
	  EndDo
	EndDo
      EndDo
End Do

IF(DebugMode) PRINT *,myid,' Done generating particles: ',ip
IF(DebugMode) PRINT *,myid,' Finding connections...'

CALL CPU_TIME(T1)
CALL FindBeams(xo, ip, SCL, NCN_All, CN_All, nbeams)
CALL CPU_TIME(T2)

PRINT *,myid,' Done finding connections: ',T2-T1,' secs'

ALLOCATE(ParticlePart(ip))
IF(myid==0) THEN

  PRINT *,'About to metis'
  ALLOCATE(metis_options(40),xadj(ip+1),adjncy(nbeams*2))
  
  !Put particle connections into CRS format
  xadj = 0
  adjncy = 0
  counter = 0
  xadj(1) = 1
  DO i=1,ip
    DO j=1,NCN_All(i)
      counter = counter + 1
      adjncy(counter) = CN_All(i,j)
    END DO
    xadj(i+1) = counter+1
  END DO

  !TODO - ensure compatibility (REAL and INT size)
  CALL METIS_SetDefaultOptions(metis_options)
  metis_options(18) = 1 !fortran array numbering

  !Use METIS to partition the particles
  CALL METIS_PartGraphKway(ip,1,xadj,adjncy,vwgt,vsize,adjwgt,ntasks,&
       tpwgts, ubvec, metis_options,objval,particlepart)

  particlepart = particlepart - 1 !mpi partitions are zero indexed

  IF(DebugMode) THEN
    PRINT *,'obvjal: ',objval
    PRINT *,'particlepart: ',MINVAL(particlepart), MAXVAL(particlepart), &
         COUNT(particlepart==0),COUNT(particlepart==ntasks-1)
    PRINT *,'--- METIS SUCCESS ---'
  END IF
  DEALLOCATE(xadj,adjncy)
END IF

IF(DebugMode) PRINT *,myid,' debug size particlepart: ',SIZE(particlepart), ip

CALL MPI_BCast(particlepart,ip,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

NN = COUNT(particlepart == myid)
ALLOCATE(particles_G(NN),NCN(NN),&
     CN(NN,12),CNpart(NN,12),&
     neighparts(ntasks),SharedNode(NN),&
     counters(ntasks),particles_L(ip))

particles_G = 0
neighparts = -1 !partitions we share a connection with
neighcount = 0
SharedNode = .FALSE.
counters = 0
particles_L = 0

!Construct all local nodenums:
!note - same results on every partition
DO i=1,ip
  counters(particlepart(i)+1) = counters(particlepart(i)+1) + 1
  particles_L(i) = counters(particlepart(i)+1)
END DO
IF(counters(myid+1) /= NN) CALL FatalError("Programming error in determining local NN")


!Construct local arrays of: 
! - node connections (CN) <= global NN!
! - number of connections (NCN)
! - other node partitions (CNpart)
! - all partitions we share connections with (neighparts)

!Allocate the structure holding the point data
CALL PointDataInit(NRXF, NN, part_expand)

counter = 0
DO i=1,ip
  IF(particlepart(i)==myid) THEN
    counter = counter + 1

    particles_G(counter) = i

    NRXF%M(:,counter) = xo(:,i) !our points
    NRXF%PartInfo(1,counter) = myid !belong to our partition
    NRXF%PartInfo(2,counter) = counter !with localID = counter
    
    NCN(counter) = NCN_All(i)
    DO j=1,NCN(counter)
      !CN_All, NCN_All
      CN(counter,j) = CN_All(i,j)
      pown = particlepart(CN_All(i,j))
      CNpart(counter,j) = pown

      !Add partition to list of neighbours if not already found
      IF((.NOT. ANY(neighparts == pown)) .AND. (myid /= pown)) THEN
        neighcount = neighcount + 1
        neighparts(neighcount) = pown
      END IF

    END DO
  END IF
END DO

IF(DebugMode) PRINT *,'DEBUG ',myid,' SUM PARTICLES_L ',SUM(particles_L),' min, max: ',&
     MINVAL(particles_L),MAXVAL(particles_L)
IF(DebugMode) PRINT *,'DEBUG ',myid,' counted ',counter,' nodes, NN: ',NN
IF(DebugMode) PRINT *,'DEBUG ',myid,' min max ncn: ',MINVAL(NCN),MAXVAL(NCN)
IF(DebugMode) PRINT *,'DEBUG ',myid,' neighparts: ',neighparts(1:neighcount)

DO i=1,NN
  IF(ANY(CNPart(i,1:NCN(i)) /= myid)) SharedNode(i) = .TRUE.
END DO

IF(DebugMode) PRINT *,myid,' shares ',COUNT(SharedNode),' of ',NN,' nodes.'

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
  END IF

  DO i=1,NN
    DO j=1,NCN(i)
      !Count each beam only once - except across boundaries, both need to count
      IF(CN(i,j) > i .OR. CNPart(i,j) /= myid) THEN
        counter = counter + 1
        IF(k==2) THEN
          NANS(1,counter) = particles_L(CN(i,j)) !Note - folows convention N1 = other part
          NANS(2,counter) = i  !NANS = their/ourNN, ourNN, otherPart (usually myid though!)
          NANpart(counter) = CNPart(i,j)
        END IF
      END IF
    END DO
  END DO
END DO

!Construct lookup arrays in NRXF for across partition beams
!Note this fills the %PartInfo but doesn't actually fill the points
!This is done by ExchangeConnPoints
CALL InvPartInfoInit(InvPartInfo, neighparts)

counter = NRXF%cstrt - 1
DO j=1,neighcount
  n = neighparts(j)
  IF(DebugMode) PRINT *,myid, 'neighbour ',j,' is ',InvPartInfo(n) % NID
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

    !TODO - expand invpartinfo % connids etc if neccessary
    InvPartInfo(n) % ConnIDs(countptr) = NANS(1,i)
    InvPartInfo(n) % ConnLocs(countptr) = counter
    
    IF(counter > SIZE(NRXF%PartInfo,2)) CALL ResizePointData(NRXF,1.5_8)
    NRXF%PartInfo(1,counter) = InvPartInfo(n) % NID
    NRXF%PartInfo(2,counter) = NANS(1,i)
  END DO
END DO
NRXF%NC = counter-NRXF%cstrt

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
! OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FS'//na(myid),STATUS='UNKNOWN')
! DO I=1,NTOT
!   WRITE(117+myid,*) NANS(1,I),NANS(2,I),NANS(3,I),&
!        NRXF%M(1,NANS(1,I)),NRXF%M(2,NANS(1,I)),&
!        NRXF%M(3,NANS(1,I)),NRXF%M(1,NANS(2,I)),&
!        NRXF%M(2,NANS(2,I)),NRXF%M(3,NANS(2,I)),EFS % M(I)
! END DO
! CLOSE (117+myid)

!TODO - keep these?
DEALLOCATE(NCN_All, CN_All)

CALL ExchangeConnPoints(NANS, NRXF, InvPartInfo)

!Write out my particle to nodfil
12    FORMAT(I8,' ',4F14.7,2I8)
OPEN(510+myid,file=TRIM(wrkdir)//'/NODFIL2'//na(myid))
DO i=1,NN
  WRITE(510+myid,12) i,NRXF%A(:,i),1.0
END DO
DO i=NRXF%cstrt, NRXF%cstrt + NRXF%NC
  WRITE(510+myid,12) i,NRXF%A(:,i),1.0,NRXF%PartInfo(:,i)
END DO
CLOSE(510+myid)

WRITE(*,*) 'NN=',NN,myid

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

!-----------------------------------------------------

Subroutine BIPINT(x,y,f11,f12,f21,f22,fint)
Implicit none
Real*8 :: x,y,f11,f12,f21,f22,fint
fint=f11*(1.0-x)*(1.0-y)+f21*x*(1.0-y)+f12*(1.0-x)*y+f22*x*y
End Subroutine

!-----------------------------------------------------

Subroutine BIPINTN(x,y,f11,f12,f21,f22,dix,diy,diz,grid)
Implicit none
Real*8 :: x,y,f11,f12,f21,f22,dix,diy,diz,grid
REAL*8 norm
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
End Subroutine

 SUBROUTINE GetBBoxes(NRXF, UT, NN, BBox, PBBox)

   IMPLICIT NONE
   INCLUDE 'mpif.h'

   INTEGER ::  NN,i
   TYPE(UT2_t) :: UT
   TYPE(NRXF2_t) :: NRXF
   !-----------------
   INTEGER :: ierr
   REAL*8 :: BBox(6),workarr(6*ntasks),X,Y,Z
   REAL*8, ALLOCATABLE :: PBBox(:,:)
   REAL*8 :: minx,maxx,miny,maxy,minz,maxz

   IF(.NOT. ALLOCATED(PBBox)) ALLOCATE(PBBox(6,ntasks))

   minx = HUGE(minx)
   miny = HUGE(miny)
   minz = HUGE(minz)
   maxx = -HUGE(maxx)
   maxy = -HUGE(maxy)
   maxz = -HUGE(maxz)
   
   DO i=1,NN
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

   DO i=1,ntasks
     PBBox(:,i) = workarr(((i-1)*6)+1:i*6)
   END DO

 END SUBROUTINE GetBBoxes

!Determine and pass point information between partitions (based on beams/connections)
SUBROUTINE ExchangeConnPoints(NANS, NRXF, InvPartInfo, UT, passNRXF)

  INCLUDE 'mpif.h'
  
  INTEGER, ALLOCATABLE :: NANS(:,:)
  TYPE(NRXF2_t) :: NRXF
  TYPE(InvPartInfo_t) :: InvPartInfo(0:)
  TYPE(UT2_t), OPTIONAL :: UT
  LOGICAL, OPTIONAL :: passNRXF
  !--------------------
  INTEGER :: i,j,id,counter,loc,neigh,getcount,sendcount,ierr
  INTEGER, ALLOCATABLE :: stats(:)
  TYPE(PointEx_t), ALLOCATABLE :: PointEx(:)
  LOGICAL :: doNRXF=.TRUE.
  
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
      ALLOCATE(PointEx(i) % S(6*sendcount),&
           PointEx(i) % R(6*getcount))


      PointEx(i) % S = 0.0
      PointEx(i) % R = 0.0

      DO j=1,sendcount
        id = PointEx(i) % SendIDs(j)
        PointEx(i) % S(j*6-5 : j*6) = UT%M(6*id-5 : 6*id)
      END DO

      CALL MPI_ISend(PointEx(i) % S, sendcount*6, MPI_DOUBLE_PRECISION,neigh,&
           201,MPI_COMM_WORLD,stats(i*2), ierr)

      CALL MPI_IRecv(PointEx(i) % R, getcount*6, MPI_DOUBLE_PRECISION,neigh,&
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
      
      ! TODO ***************
      ! IF((counter + getcount)*6 > SIZE(UT%P)) THEN
      !   PRINT *,myid,' debug, resizing UT'
      !   CALL ResizePointData(UT,1.5_8)
      ! END IF
      
      DO j=1,getcount
        loc = InvPartInfo(i) % ConnLocs(j)
        UT%A(6*loc-5 : 6*loc) = PointEx(i) % R(j*6-5 : j*6)
      END DO
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE ExchangeConnPoints

!Determine and pass point information between partitions (based on proximity)
SUBROUTINE ExchangeProxPoints(NRXF, UT, NN, SCL)

  USE INOUT
  USE UTILS

  INCLUDE 'mpif.h'

  TYPE(NRXF2_t) :: NRXF
  TYPE(UT2_t) :: UT
  INTEGER :: NN
  REAL*8 :: SCL
  !--------------------
  INTEGER :: i,j,id,cnt,neigh,neighcount,getcount,stat(MPI_STATUS_SIZE),ierr
  INTEGER, ALLOCATABLE :: WorkInt(:),stats(:)
  REAL*8, ALLOCATABLE :: PBBox(:,:)
  REAL*8 :: BBox(6)
  TYPE(PointEx_t), ALLOCATABLE :: PointEx(:)

  ALLOCATE(PointEx(ntasks), stats(ntasks*2))
  stats = MPI_REQUEST_NULL

  DO i=1,ntasks
    IF(i==myid+1) CYCLE
    ALLOCATE(PointEx(i) % SendIDs(NN/2))
    PointEx(i) % partid = i-1
    PointEx(i) % scount = 0
    PointEx(i) % rcount = 0
  END DO

  !Pass bounding boxes
  CALL GetBBoxes(NRXF, UT, NN,  BBox, PBBox)


  !Count points in each BB
  DO i=1,ntasks
    IF(i==myid+1) CYCLE
    cnt = 0
    DO j=1,NN
      IF(PInBBox(j,NRXF,UT,PBBox(:,i),SCL*2.0)) THEN
        PointEx(i) % scount = PointEx(i) % scount + 1
        cnt = PointEx(i) % scount

        IF(cnt > SIZE(PointEx(i) % SendIDs)) CALL ExpandIntArray(&
             PointEx(i) % SendIDs)

        PointEx(i) % SendIDs(cnt) = j
      END IF
    END DO

    ALLOCATE(PointEx(i) % S(3 * PointEx(i) % scount))
    CALL MPI_ISend(PointEx(i) % scount,1,MPI_INTEGER,i-1,201,MPI_COMM_WORLD,stats(i*2-1),ierr)
    CALL MPI_ISend(PointEx(i) % SendIDs(1:cnt),cnt,MPI_INTEGER,i-1,202,MPI_COMM_WORLD,stats(i*2),ierr)
  END DO


  !Receive count and IDs sent above
  DO i=1,ntasks
    IF(i==myid+1) CYCLE

    CALL MPI_Recv(PointEx(i) % rcount, 1, MPI_INTEGER, i-1,201,MPI_COMM_WORLD, stat, ierr)
    cnt = PointEx(i) % rcount

    ALLOCATE(PointEx(i) % R(3 * cnt),&
         PointEx(i) % RecvIDs(cnt))

    CALL MPI_Recv(PointEx(i) % RecvIDs(1:cnt),cnt, MPI_INTEGER, i-1,202,MPI_COMM_WORLD, stat, ierr)

  END DO

  !Wait for the previous non-blocking sends, then reset stats
  CALL MPI_Waitall(ntasks*2, stats, MPI_STATUSES_IGNORE, ierr)
  stats = MPI_REQUEST_NULL

  !Pass NRXF
  DO i=1,ntasks

    IF(i==myid+1) CYCLE
    DO j=1,PointEx(i) % scount
      id = PointEx(i) % SendIDs(j)
      PointEx(i) % S(j*3 - 2) = NRXF%M(1,id)
      PointEx(i) % S(j*3 - 1) = NRXF%M(2,id)
      PointEx(i) % S(j*3 - 0) = NRXF%M(3,id)
    END DO

    CALL MPI_ISend(PointEx(i) % S, PointEx(i) % scount*3, MPI_DOUBLE_PRECISION, &
         i-1, 203, MPI_COMM_WORLD, stats(i*2-1), ierr)
    CALL MPI_IRecv(PointEx(i) % R, PointEx(i) % rcount*3, MPI_DOUBLE_PRECISION, &
         i-1, 203, MPI_COMM_WORLD, stats(i*2), ierr)
  END DO

  CALL MPI_Waitall(ntasks*2, stats, MPI_STATUSES_IGNORE, ierr)
  stats = MPI_REQUEST_NULL


  !Store NRXF and Pass UT
  DO i=1,ntasks
    IF(i==myid+1) CYCLE

    !TODO - STORE NRXF DATA SOMEWHERE    

    PointEx(i) % S = 0.0
    PointEx(i) % R = 0.0
    DO j=1,PointEx(i) % scount
      id = PointEx(i) % SendIDs(j)
      PointEx(i) % S(j*3 - 2) = UT%M(6*id - 5)
      PointEx(i) % S(j*3 - 1) = UT%M(6*id - 4)
      PointEx(i) % S(j*3 - 0) = UT%M(6*id - 3)
    END DO

    CALL MPI_ISend(PointEx(i) % S, PointEx(i) % scount*3, MPI_DOUBLE_PRECISION, &
         i-1, 203, MPI_COMM_WORLD, stats(i*2-1),ierr)
    CALL MPI_IRecv(PointEx(i) % R, PointEx(i) % rcount*3, MPI_DOUBLE_PRECISION, &
         i-1, 203, MPI_COMM_WORLD, stats(i*2),ierr)
    
  END DO

  CALL MPI_Waitall(ntasks*2, stats, MPI_STATUSES_IGNORE, ierr)
  stats = MPI_REQUEST_NULL

  !Store UT
  DO i=1,ntasks
    IF(i==myid+1) CYCLE

    !TODO - DO SOMETHING WITH UT

  END DO

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

  INCLUDE 'mpif.h'

  INTEGER, ALLOCATABLE :: NANS(:,:), NANPart(:)
  REAL*8, ALLOCATABLE :: EFS(:)
  TYPE(NRXF2_t) :: NRXF
  !-----------------------------
  INTEGER :: i,j,k,neigh,counter,N1,N2,NTOT,ierr
  INTEGER, ALLOCATABLE :: stats(:)
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

      CALL MPI_IRecv(PointEx(i) % RecvIDs,counter*2,MPI_DOUBLE_PRECISION,neigh,&
           191,MPI_COMM_WORLD,stats(i*2+1),ierr)

    END IF
  END DO

  CALL MPI_WaitAll(ntasks*2, stats, MPI_STATUSES_IGNORE, ierr)

  DO i=0,ntasks-1
    IF(InvPartInfo(i) % ccount == 0) CYCLE
    neigh = InvPartInfo(i) % NID
    Send = (neigh < myid) 
    IF(Send) CYCLE

    counter = 0
    DO j=1,PointEx(i) % rcount
      N2 = PointEx(i) % RecvIDs(j*2-1) !our ID
      N1 = PointEx(i) % RecvIDs(j*2)   !otherpart ID

      DO k=1,InvPartInfo(i) % CCount
        IF(InvPartInfo(i) % ConnIDs(k) == N1) THEN
          N1 = InvPartInfo(i) % ConnLocs(k)
          EXIT
        END IF
      END DO

      DO k=1,NTOT
        IF(NANPart(k) /= neigh) CYCLE
        IF(NANS(1,k) /= N1) CYCLE
        IF(NANS(2,k) /= N2) CYCLE
        PRINT *,myid,' EFS ',k,' changing from ',EFS(k),' to ',PointEx(i) % R(j)
        EFS(k) = PointEx(i) % R(j)
        counter = counter + 1
      END DO
    END DO

    IF(counter /= PointEx(i) % rcount) &
         CALL FatalError("Programming error: count mismatch in EFS comms")
    IF(DebugMode) PRINT *,myid,' received ',counter,PointEx(i) % rcount,' EFS values from ',neigh
  END DO

  CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

END SUBROUTINE ExchangeEFS

!Use octree search to find nodes which are in contact
SUBROUTINE FindCollisions(NRXF, UT, NN, BBox,SCL,LNN)

  USE Octree

  TYPE(NRXF2_t) :: NRXF
  TYPE(UT2_t) :: UT
  INTEGER :: NN
  REAL*8 :: SCL, LNN, BBox(6)
  !-----------------------
  !octree stuff
  type(point_type), allocatable :: points(:)
  INTEGER i,j, npoints, num_ngb
  integer, allocatable :: seed(:), ngb_ids(:)
  REAL*8 :: x(3), dx(3),oct_bbox(2,3),eps

  REAL*8 :: dist, max_dist, min_dist

  eps = 1.0

  npoints = NN !TODO - NN => include parallel
  ALLOCATE(points(npoints)) 

  DO i=1,3
    oct_bbox(1,i) = BBox(i*2-1) - eps !buffer to ensure all points contained
    oct_bbox(2,i) = BBox(i*2)   + eps
  END DO

  CALL Octree_init(max_depth=20,max_num_point=6,bbox=oct_bbox)

  DO i=1,npoints
    Points(i) % id = i
    Points(i) % x(1) = NRXF%A(1,i) + UT%A(6*I-5)
    Points(i) % x(2) = NRXF%A(2,i) + UT%A(6*I-4)
    Points(i) % x(3) = NRXF%A(3,i) + UT%A(6*I-3)
  END DO

  CALL Octree_build(Points)

  ALLOCATE(ngb_ids(100))

  DO i=1,npoints
    num_ngb = 0
    ngb_ids = 0
    CALL Octree_search(points(i) % x, SCL*1.87, num_ngb, ngb_ids)
    max_dist = -HUGE(max_dist)
    min_dist = HUGE(max_dist)
    DO j=1,num_ngb
      IF(ngb_ids(j) == i) CYCLE
      dist = ((((NRXF%A(1,i) + UT%A(6*i-5)) - (NRXF%A(1,ngb_ids(j)) + UT%A(6*ngb_ids(j)-5)))**2.0) + &
      (((NRXF%A(2,i) + UT%A(6*i-4)) - (NRXF%A(2,ngb_ids(j)) + UT%A(6*ngb_ids(j)-4)))**2.0) + &
      (((NRXF%A(3,i) + UT%A(6*i-3)) - (NRXF%A(3,ngb_ids(j)) + UT%A(6*ngb_ids(j)-3)))**2.0)) ** 0.5_dp

      max_dist = MAX(max_dist, dist)
      min_dist = MIN(min_dist, dist)
    END DO
    IF(DebugMode) PRINT *,myid,' node ',i,' noneigh: ',num_ngb,' max/min dist: ',max_dist, min_dist
  END DO

  CALL Octree_final()

END SUBROUTINE FindCollisions


!Use octree search to find nodes which are in contact
SUBROUTINE FindBeams(xo, ip, SCL, NCN, CN, nbeams)
  
  USE Octree

  REAL*8, ALLOCATABLE :: xo(:,:)
  REAL*8 :: SCL
  INTEGER :: ip,nbeams
  INTEGER, ALLOCATABLE :: NCN(:), CN(:,:)
  !-----------------------
  !octree stuff
  type(point_type), allocatable :: points(:)
  INTEGER i,j,k, num_ngb,counter
  integer, allocatable :: seed(:), ngb_ids(:)
  REAL*8 :: x(3), dx(3),oct_bbox(2,3),eps
  REAL*8 :: dist, max_dist, min_dist, searchdist

  ALLOCATE(points(ip),&
       NCN(ip),&
       CN(ip,12)) 

  eps = 1.0

  nbeams = 0
  CN = 0
  NCN = 0

  DO i=1,3
    oct_bbox(1,i) = MINVAL(xo(i,:)) - eps !buffer bbox to ensure all points contained
    oct_bbox(2,i) = MAXVAL(xo(i,:)) + eps
  END DO


  searchdist = SCL * (1.6**0.5)

  !TODO - determine optimal depth - SCL vs bbox length?
  !  mdepth 20, m points 6 was determined for a particular case
  !   it doubled the speed w.r.t. 20,3 (i.e. max_num_points is important)
  !   we *expect* 12 neighbours for every node... 
  CALL Octree_init(max_depth=20,max_num_point=6,bbox=oct_bbox)

  DO i=1,ip
    Points(i) % id = i
    Points(i) % x(1) = xo(1,i)
    Points(i) % x(2) = xo(2,i)
    Points(i) % x(3) = xo(3,i)
  END DO
  !octree build misses points...
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
      NCN(i) = NCN(i) + 1
      CN(i,counter) = ngb_ids(j)
      ! PRINT *,i, j, ' ngb_id ',ngb_ids(j)
      IF(ngb_ids(j) > i) nbeams = nbeams+1
    END DO
    IF(DebugMode .AND. myid==0 .AND. MOD(i,1000)==0)&
         PRINT *,myid,' node ',i,' nobeams: ',num_ngb
  END DO

  PRINT *,myid,' sum(ncn) ',SUM(ncn),' nbeams*2 ',nbeams*2
  !SUM(NCN) is too large compared to nbeams
  !For some reason, testing id(j) > i doesn't produce right number of beams
CALL Octree_final()

END SUBROUTINE FindBeams

FUNCTION PInBBox(i,NRXF, UT, BBox, Buff) RESULT(InBB)
  INTEGER :: i
  TYPE(NRXF2_t) :: NRXF
  TYPE(UT2_t) :: UT
  REAL*8 :: BBox(6)
  REAL*8, OPTIONAL :: Buff
  LOGICAL :: InBB
  !------------------------
  REAL*8 :: X,Y,Z, Buffer

  IF(PRESENT(Buff)) THEN
    Buffer = Buff
  ELSE
    Buffer = 200.0
  END IF

  X = NRXF%A(1,i) + UT%A(6*i - 5)
  Y = NRXF%A(2,i) + UT%A(6*i - 4)
  Z = NRXF%A(3,i) + UT%A(6*i - 3)
  
  InBB = (X > BBox(1)-Buffer .AND. X < BBox(2)+Buffer .AND. &
       Y > BBox(3)-Buffer .AND. Y < BBox(4)+Buffer .AND. &
       Z > BBox(5)-Buffer .AND. Z < BBox(6)+Buffer)
END FUNCTION PInBBox

END MODULE Lattice
