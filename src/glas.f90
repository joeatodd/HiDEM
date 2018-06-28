!defines and makes the FCC lattice - dense packing

!ip - returns number of nodes (in this partition)
!ntasks - how many cores
!myid - this partition id
SUBROUTINE FIBG3(l,ip,ND,myid,maxx,&
     maxy,maxz,minx,miny,minz,ntasks,wrkdir,geomfile,SCL,YN,XN,grid,melta,wl,UC,StrictDomain)

USE TypeDefs

Implicit none
INCLUDE 'na90.dat'
!Real*8,ALLOCATABLE :: surf(:,:),bed(:,:)
Real*8 surf(-100:2000,-100:2000),bed(-100:2000,-100:2000),melt(-100:2000,-100:2000)
Real*8 :: x,y,s1,b1,b2,u1,grid,m1,melta,wl,UC,z1
Real*8 :: box,b,maxx,maxy,maxz,minx,miny,minz,SCL
INTEGER :: l,ip,i,j,YN,XN
INTEGER :: myid,ntasks,N1,N2,xk,yk
CHARACTER(LEN=256) :: wrkdir,geomfile
LOGICAL :: StrictDomain
TYPE(NTOT_t) :: ND

!Open(300,file='mass.dat',STATUS='OLD')
OPEN(510,file=TRIM(wrkdir)//'/NODFIL2')

 12    FORMAT(2I8,' ',2F14.7)
 13    FORMAT(4F14.7)
 14    FORMAT(3F14.7)

surf=-100.0
bed=1000.0

!TODO pass through file once to check extent, then allocate bed, surf, melt
! then reread
!TODO - error message if this file doesn't exist
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

maxx=0.0
maxy=0.0
maxz=0.0
minx=1.0e+08
miny=1.0e+08
minz=1.0e+08

!box is never actualy used...
!b is used, which is box/l, so L never actually enters into this, so it's only
!used for the number of vertical layers
box=2.0d0**(2.0d0/3.0d0)*REAL(l,8) ! box size equal to fcc ground state

CALL Initializefcc(box,l,ip,myid,maxx,maxy,maxz,minx,miny,minz,ntasks,wrkdir,&
     SCL,YN,XN,surf,bed,melt,grid,wl,UC,StrictDomain)

CLOSE(510)
CLOSE(400)

END SUBROUTINE FIBG3

!---------------------------------------------------------------!


SUBROUTINE Initializefcc(box,l,ip,myid,maxx,maxy,maxz,minx,miny,minz,ntasks,wrkdir,&
     SCL,YN,XN,surf,bed,melt,grid,wl,UC,StrictDomain)

  USE INOUT
  USE iso_c_binding

Implicit None
INCLUDE 'na90.dat'
INCLUDE 'mpif.h'

Real*8 surf(-100:2000,-100:2000),bed(-100:2000,-100:2000),melt(-100:2000,-100:2000)
Real*8 b,x0(3,4),box,maxx,maxy,maxz,minx,miny,minz,SCL
REAL*8 gridminx, gridmaxx, gridminy, gridmaxy
REAL*8 z,x,y,sint,bint,mint,grid,wl,lc,UC,UCV,YN_estimate, XN_estimate, efficiency
REAL*8 :: X1,X2,Y1,Y2,Z1,Z2,RC
REAL*8, ALLOCATABLE :: xo(:,:),work_arr(:,:)

INTEGER i,j,k,k1,k2,l,ip,NN,myid,ntasks,m,nb,YN,XN,xk,yk,nx,ny,nbeams,ierr,pown
INTEGER neighcount
INTEGER, ALLOCATABLE :: NCN_All(:), CN_All(:,:), NCN(:),CN(:,:), CNPart(:,:),&
     particles_G(:),neighparts(:)
CHARACTER(LEN=256) :: wrkdir
LOGICAL :: StrictDomain
LOGICAL, ALLOCATABLE :: SharedNode(:)

!metis stuff
INTEGER :: objval,counter
INTEGER, ALLOCATABLE :: metis_options(:), particlepart(:)
INTEGER, POINTER :: vwgt=>NULL(),vsize=>NULL(),adjwgt=>NULL(),xadj(:),adjncy(:)
REAL*8, POINTER :: tpwgts=>NULL(),ubvec=>NULL()


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
m=MOD(myid,ntasks/YN)
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
             if (xo(1,ip).gt.maxx) maxx=xo(1,ip)
             if (xo(2,ip).gt.maxy) maxy=xo(2,ip)
             if (xo(3,ip).gt.maxz) maxz=xo(3,ip)
             if (xo(1,ip).lt.minx) minx=xo(1,ip)
             if (xo(2,ip).lt.miny) miny=xo(2,ip)
             if (xo(3,ip).lt.minz) minz=xo(3,ip)
             ! EndIF
             ! EndIF
             EndIF
!             EndIF
	  EndDo
	EndDo
      EndDo
End Do

 12    FORMAT(I8,' ',4F14.7)

IF(myid==0) THEN
  DO i=1,ip
    WRITE(510,12) i,xo(:,i),1.0
  END DO
  WRITE(*,*) 'NN=',ip,myid
END IF

!Find connections:

ALLOCATE(NCN_All(ip), CN_All(ip,12))
nbeams = 0
CN_All = 0
NCN_All = 0

DO I=1,ip
  DO J=I+1,ip
    X1=xo(1,I)
    Y1=xo(2,I)
    Z1=xo(3,I)
    X2=xo(1,J)
    Y2=xo(2,J)
    Z2=xo(3,J)
    RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
    IF (RC.LT.SCL*SCL*1.6) THEN
      nbeams = nbeams + 1
      NCN_All(I)=NCN_All(I)+1
      NCN_All(J)=NCN_All(J)+1
      CN_All(I,NCN_All(I)) = J
      CN_All(J,NCN_All(J)) = I
    ENDIF
  END DO
END DO

! IF(myid==0) THEN
!   OPEN(UNIT=210,FILE=TRIM(wrkdir)//'/FSGRAPH',STATUS='UNKNOWN')
!   WRITE(210,*) ip, nbeams
!   DO I=1,ip
!     WRITE(210,*) CN_All(I,1:NCN_All(I))
!   END DO

!   CLOSE(210)
! END IF

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
  CALL METIS_PartGraphKway(ip,1,xadj,adjncy,vwgt,vsize,adjwgt,ntasks,tpwgts, ubvec, metis_options,objval,particlepart)

  particlepart = particlepart - 1 !mpi partitions are zero indexed

  PRINT *,'obvjal: ',objval
  PRINT *,'particlepart: ',MINVAL(particlepart), MAXVAL(particlepart), COUNT(particlepart==0),COUNT(particlepart==ntasks-1)
  PRINT *,'--- METIS SUCCESS ---'

  DEALLOCATE(xadj,adjncy)
END IF

CALL MPI_BCast(particlepart,ip,MPI_INTEGER,0,MPI_COMM_WORLD)

NN = COUNT(particlepart == myid)
ALLOCATE(particles_G(NN),NCN(NN),&
     CN(NN,12),CNpart(NN,12),&
     neighparts(ntasks),SharedNode(NN))

particles_G = 0
neighparts = -1 !partitions we share a connection with
neighcount = 0
SharedNode = .FALSE.

!Construct local arrays of: 
! - node connections (CN) 
! - number of above (NCN)
! - other node partitions (CNpart)
! - all partitions we share connections with (neighparts)
counter = 0
DO i=1,ip
  IF(particlepart(i)==myid) THEN
    counter = counter + 1
    particles_G(counter) = i
    NCN(counter) = NCN_All(i)
    DO j=1,NCN_All(i)
      !CN_All, NCN_All
      CN(counter,j) = CN_All(i,j)
      pown = particlepart(CN_All(i,j))
      CNpart(counter,j) = pown
      IF((.NOT. ANY(neighparts == pown)) .AND. (myid /= pown)) THEN
        neighcount = neighcount + 1
        neighparts(neighcount) = pown
      END IF
    END DO
  END IF
END DO

PRINT *,'DEBUG ',myid,' counted ',counter,' nodes, NN: ',NN
PRINT *,'DEBUG ',myid,' min max ncn: ',MINVAL(NCN),MAXVAL(NCN)
PRINT *,'DEBUG ',myid,' neighparts: ',neighparts(1:neighcount)
DO i=1,NN
  IF(ANY(CNPart(i,1:NCN(i)) /= myid)) SharedNode(i) = .TRUE.
END DO

PRINT *,myid,' shares ',COUNT(SharedNode),' of ',NN,' nodes.'

!TODO - keep these?
DEALLOCATE(NCN_All, CN_All)
STOP

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

