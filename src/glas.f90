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
     wrkdir,geomfile,SCL,grid,melta,wl,UC,StrictDomain,GeomMasked)

  IMPLICIT NONE
  INCLUDE 'na90.dat'

!Real(KIND=dp),ALLOCATABLE :: surf(:,:),bed(:,:)
Real(KIND=dp) :: surf(-100:2000,-100:2000),bed(-100:2000,-100:2000),melt(-100:2000,-100:2000)
Real(KIND=dp) :: x,y,s1,b1,b2,u1,grid,m1,melta,wl,UC,z1
Real(KIND=dp) :: box,b,SCL
INTEGER :: l,NN,i,j,mask
INTEGER :: N1,N2,xk,yk,neighcount,NTOT
INTEGER, ALLOCATABLE :: NCN(:),CN(:,:),CNPart(:,:), particles_G(:),NANS(:,:),NANPart(:)
CHARACTER(LEN=256) :: wrkdir,geomfile
LOGICAL :: StrictDomain,GeomMasked
!TYPE(NTOT_t) :: NTOT
TYPE(NRXF_t) :: NRXF
TYPE(InvPartInfo_t), ALLOCATABLE :: InvPartInfo(:)
!Open(300,file='mass.dat',STATUS='OLD')

 12    FORMAT(2I8,' ',2F14.7)
 13    FORMAT(4F14.7)
 14    FORMAT(3F14.7)

surf=-100.0_dp
bed=1000.0_dp
melta = 0.0_dp

!TODO pass through file once to check extent, then allocate bed, surf, melt
! then reread
OPEN(400,file=TRIM(geomfile),STATUS='OLD')
READ(400,*) N2
DO I=1,N2
  IF(GeomMasked) THEN
    READ(400,*) x,y,s1,b1,b2,z1,mask
  ELSE
    READ(400,*) x,y,s1,b1,b2,z1
  END IF
  xk=INT(x/grid)
  yk=INT(y/grid)
  IF (xk.GE.-100.AND.yk.GE.-100 .AND. xk.LE.2000.AND.yk.LE.2000) THEN
    bed(xk,yk)=b1
    surf(xk,yk)=s1
    melt(xk,yk)=melta*0.0
  END IF
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

Real(KIND=dp) :: surf(-100:2000,-100:2000),bed(-100:2000,-100:2000),melt(-100:2000,-100:2000)
Real(KIND=dp) :: b,x0(3,4),box,SCL
REAL(KIND=dp) :: gridminx, gridmaxx, gridminy, gridmaxy,T1,T2
REAL(KIND=dp) :: z,x,y,sint,bint,mint,grid,wl,lc,UC,UCV, efficiency
REAL(KIND=dp) :: X1,X2,Y1,Y2,Z1,Z2,RC,rc_min,rc_max
REAL(KIND=dp), ALLOCATABLE :: xo(:,:),work_arr(:,:)

INTEGER i,j,k,n,K1,k2,l,ip,NN,nb,xk,yk,nx,ny,nbeams,ierr,pown
INTEGER NTOT,neighcount
INTEGER, ALLOCATABLE :: NCN_All(:), CN_All(:,:),NCN(:),CN(:,:),CNPart(:,:),&
     particles_G(:),particles_L(:),NANS(:,:),NANPart(:)
CHARACTER(LEN=256) :: wrkdir
LOGICAL :: StrictDomain
LOGICAL, ALLOCATABLE :: SharedNode(:),neighparts(:)

!metis stuff
INTEGER :: objval,counter
INTEGER, ALLOCATABLE :: metis_options(:), particlepart(:),counters(:)
INTEGER, POINTER :: vwgt=>NULL(),vsize=>NULL(),adjwgt=>NULL(),xadj(:),adjncy(:),countptr
REAL(KIND=dp), POINTER :: tpwgts=>NULL(),ubvec=>NULL()

!TYPE(NTOT_t) :: NTOT
TYPE(NRXF_t), TARGET :: NRXF
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
               IF(DebugMode) CALL Warn("Doubling size of xo")
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

IF(PrintTimes) CALL CPU_TIME(T1)
CALL FindBeams(xo, ip, SCL, NCN_All, CN_All, nbeams)
IF(PrintTimes) CALL CPU_TIME(T2)

IF(PrintTimes) PRINT *,myid,' Done finding connections: ',T2-T1,' secs'

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
     neighparts(0:ntasks-1),SharedNode(NN),&
     counters(ntasks),particles_L(ip))

particles_G = 0
neighparts = .FALSE. !partitions we share a connection with
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

      IF(pown /= myid) neighparts(pown) = .TRUE. !partition is a neighbour
    END DO
  END IF
END DO

neighcount = COUNT(neighparts)

IF(DebugMode) PRINT *,'DEBUG ',myid,' SUM PARTICLES_L ',SUM(particles_L),' min, max: ',&
     MINVAL(particles_L),MAXVAL(particles_L)
IF(DebugMode) PRINT *,'DEBUG ',myid,' counted ',counter,' nodes, NN: ',NN
IF(DebugMode) PRINT *,'DEBUG ',myid,' min max ncn: ',MINVAL(NCN),MAXVAL(NCN)

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
      IF(particles_L(CN(i,j)) > i .OR. CNPart(i,j) /= myid) THEN
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

!-----------------------------------------------------

Subroutine BIPINT(x,y,f11,f12,f21,f22,fint)
Implicit none
Real(KIND=dp) :: x,y,f11,f12,f21,f22,fint
fint=f11*(1.0-x)*(1.0-y)+f21*x*(1.0-y)+f12*(1.0-x)*y+f22*x*y
End Subroutine

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
End Subroutine

SUBROUTINE GetBBoxes(NRXF, UT, NN, NANS, NTOT, EFS, BBox, PBBox)

   IMPLICIT NONE
   INCLUDE 'mpif.h'

   INTEGER ::  NN, NTOT
   INTEGER :: NANS(:,:)
   REAL(KIND=dp), ALLOCATABLE :: EFS(:)
   TYPE(UT_t) :: UT
   TYPE(NRXF_t) :: NRXF
   !-----------------
   INTEGER :: i,j,ierr,max_gap_loc(3),idx_ulimit, idx_llimit
   REAL(KIND=dp) :: BBox(6),workarr(6*ntasks),X,Y,Z,bbox_vol,bbox_vols(ntasks)
   REAL(KIND=dp) :: PBBox(6,0:ntasks),max_gap(3),gap(3)
   REAL(KIND=dp) :: outlier_gap_prop, outlier_arr_prop
   REAL(KIND=dp) :: minx,maxx,miny,maxy,minz,maxz,midx,midy,midz
   REAL(KIND=dp) :: dist(3,NN), pos(3,NN), sigma_dist(3,NN),mean_dist(3),std_dev(3)
   INTEGER :: dist_sort(3,NN), pos_sort(3,NN)
   LOGICAL :: outlier(3,NN), poutlier(NN),check_outliers

   check_outliers = .FALSE.

   ! outlier_gap_prop = 0.4
   ! outlier_arr_prop = 0.05

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
   !than average bbox volume will check for outlier particles
   bbox_vol = (maxx - minx) * (maxy - miny) * (maxz - minz)
   CALL MPI_AllGather(bbox_vol,1,MPI_DOUBLE_PRECISION,bbox_vols,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)

   IF(DebugMode) THEN
     PRINT *,myid,' debug bbox vol: ',bbox_vol
     IF(myid==0) PRINT *,' mean bbox vol: ',(SUM(bbox_vols)/ntasks)
   END IF

   IF(bbox_vol / (SUM(bbox_vols)/ntasks) > 2.0) THEN
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
       poutlier(i) = ANY(ABS(sigma_dist(:,i)) > 3.0)
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

     !Recompute BB w/o outliers
     minx = HUGE(minx)
     miny = HUGE(miny)
     minz = HUGE(minz)
     maxx = -HUGE(maxx)
     maxy = -HUGE(maxy)
     maxz = -HUGE(maxz)
     DO i=1,NN
       IF(poutlier(i)) CYCLE
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

     IF(ANY(poutlier) .AND. DebugMode) THEN
       PRINT *,myid,' new bbox: ',minx,maxx,miny,maxy,minz,maxz
     END IF

   END IF !check outliers

   CALL MPI_AllGather(BBox,6,MPI_DOUBLE_PRECISION,workarr,6,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)

   DO i=0,ntasks-1
     PBBox(:,i) = workarr((i*6)+1:(i+1)*6)
   END DO

 END SUBROUTINE GetBBoxes

 SUBROUTINE FindNeighbours(PBBox,PartIsNeighbour)
   REAL(KIND=dp) :: PBBox(:,0:)
   LOGICAL :: PartIsNeighbour(0:ntasks-1)
   !-------------------------------
   INTEGER :: i
   REAL(KIND=dp) :: BBox(6),Buffer

   BBox = PBBox(:,myid)
   Buffer = 200.0

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

  INCLUDE 'mpif.h'
  
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

  INCLUDE 'mpif.h'

  TYPE(NRXF_t) :: NRXF
  TYPE(UT_t) :: UT, UTM
  INTEGER :: NN
  REAL(KIND=dp) :: SCL, PBBox(6,0:ntasks-1)
  LOGICAL :: PartIsNeighbour(0:ntasks-1)
  TYPE(InvPartInfo_t) :: InvPartInfo(0:)
  !--------------------
  REAL(KIND=dp) :: T1, T2, tstrt,tend
  INTEGER :: i,j,k,id,new_id,put_loc,put_loc_init,cnt,cnt2,rmcnt,neigh,neighcount,&
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

  INCLUDE 'mpif.h'

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
SUBROUTINE FindNearbyParticles(NRXF, UT, NN, BBox,SCL,LNN,ND,NDL)

  USE Octree

  TYPE(NRXF_t) :: NRXF
  TYPE(UT_t) :: UT
  INTEGER :: NN,ND
  INTEGER, ALLOCATABLE :: NDL(:,:)
  REAL(KIND=dp) :: SCL, LNN, BBox(6)
  !-----------------------
  !octree stuff
  type(point_type), allocatable :: points(:)
  INTEGER i,j,cnt,npoints, num_ngb, totsize
  INTEGER, ALLOCATABLE :: seed(:), ngb_ids(:), point_loc(:)
  REAL(KIND=dp) :: x(3), dx(3),oct_bbox(2,3),eps,T1,T2
  REAL(KIND=dp) :: dist, max_dist, min_dist,tstrt,tend

  REAL(KIND=dp) :: X1,X2,Y1,Y2,Z1,Z2
  INTEGER :: id
  CALL CPU_Time(tstrt)

  eps = 1.0

  npoints = COUNT(NRXF%PartInfo(1,:) /= -1)
  totsize = SIZE(NRXF%PartInfo,2)

  IF(DebugMode) PRINT *, myid, 'debug nrxf count: ',nn,npoints,SIZE(NRXF%PartInfo,2)

  IF(.NOT. ALLOCATED(NDL)) ALLOCATE(NDL(2,npoints*12)) !guess size
  ALLOCATE(points(npoints),point_loc(npoints)) 

  point_loc = 0
  ND = 0
  NDL = 0

  DO i=1,3
    oct_bbox(1,i) = BBox(i*2-1) - eps !buffer to ensure all points contained
    oct_bbox(2,i) = BBox(i*2)   + eps
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
    CALL Octree_search(points(i) % x, SCL*1.87, num_ngb, ngb_ids)

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
SUBROUTINE FindBeams(xo, ip, SCL, NCN, CN, nbeams)
  
  USE Octree

  REAL(KIND=dp), ALLOCATABLE :: xo(:,:)
  REAL(KIND=dp) :: SCL
  INTEGER :: ip,nbeams
  INTEGER, ALLOCATABLE :: NCN(:), CN(:,:)
  !-----------------------
  !octree stuff
  type(point_type), allocatable :: points(:)
  INTEGER i,j,k, num_ngb,counter
  integer, allocatable :: seed(:), ngb_ids(:)
  REAL(KIND=dp) :: x(3), dx(3),oct_bbox(2,3),eps
  REAL(KIND=dp) :: dist, max_dist, min_dist, searchdist

  ALLOCATE(points(ip),&
       NCN(ip),&
       CN(ip,12)) 

  eps = 1.0

  nbeams = 0
  CN = 0
  NCN = 0

  DO i=1,3
    oct_bbox(1,i) = MINVAL(xo(i,1:ip)) - eps !buffer bbox to ensure all points contained
    oct_bbox(2,i) = MAXVAL(xo(i,1:ip)) + eps
  END DO


  searchdist = SCL * (1.6**0.5)

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

  IF(DebugMode) PRINT *,myid,' sum(ncn) ',SUM(ncn),' nbeams*2 ',nbeams*2
  !SUM(NCN) is too large compared to nbeams
  !For some reason, testing id(j) > i doesn't produce right number of beams
CALL Octree_final()

END SUBROUTINE FindBeams

SUBROUTINE FindCollisions(ND,NN,NRXF,UT,FRX,FRY,FRZ, &
     T,IS,DT,WE,EFC,FXF,FXC,NDL,LNN,SCL)

  USE TypeDefs
  USE Utils

  IMPLICIT NONE
  INCLUDE 'mpif.h'
  REAL(KIND=dp) ::  X1,X2,Y1,Y2,Z1,Z2
  REAL(KIND=dp) ::  T1,T2
  REAL(KIND=dp), ALLOCATABLE :: EFC(:)
  REAL(KIND=dp) ::  SX,SY,SZ,SUM,T,WE(:),L0
  REAL(KIND=dp) ::  DDEL,DWE,OWE,DT,ESUM,LNN
  REAL(KIND=dp) ::  LS,LS2,DEL,SCL
  INTEGER ierr,FXC,ND
  INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),comm
  INTEGER, ALLOCATABLE :: FXF(:,:),NDL(:,:)
  REAL(KIND=dp) ::  RC,RCX,RCY,RCZ,FRX(NN),FRY(NN),FRZ(NN)
  INTEGER NTOT,I,N1,N2,IS,NN
  TYPE(UT_t) :: UT
  TYPE(NRXF_t) :: NRXF
  LOGICAL :: own(2)

  IF(PrintTimes) CALL CPU_TIME(T1)

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
      IF (RC.GT.LNN.AND.RC.LT.LNN+0.04*SCL) THEN
        IF(own(1)) THEN
          FRX(N1)=FRX(N1)+SCL**2.0*1.0e+04*(LNN-RC)*RCX
          FRY(N1)=FRY(N1)+SCL**2.0*1.0e+04*(LNN-RC)*RCY
          FRZ(N1)=FRZ(N1)+SCL**2.0*1.0e+04*(LNN-RC)*RCZ
        END IF
        IF(own(2)) THEN
          FRX(N2)=FRX(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCX
          FRY(N2)=FRY(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCY
          FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCZ
        END IF
        IF(own(2)) THEN
          WE(N2)=WE(N2)+SCL**2.0*0.5e+04*(LNN-RC)**2.0
        ELSE
          WE(N1)=WE(N1)+SCL**2.0*0.5e+04*(LNN-RC)**2.0
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
  
  InBB = (X > BBox(1)-Buffer .AND. X < BBox(2)+Buffer .AND. &
       Y > BBox(3)-Buffer .AND. Y < BBox(4)+Buffer .AND. &
       Z > BBox(5)-Buffer .AND. Z < BBox(6)+Buffer)
END FUNCTION PInBBox

END MODULE Lattice
