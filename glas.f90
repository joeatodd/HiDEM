Subroutine FIBG3(l,ip,ND,NDLC,NDRC,NDFC,NDBC,NDFRC,NDBRC,NDFLC,NDBLC,myid,maxx,maxy,maxz,minx,miny,minz,ntasks,SCL,YN,grid,melta,wl,UC)

Implicit none
INCLUDE 'na90.dat'
!Real*8,ALLOCATABLE :: surf(:,:),bed(:,:)
Real*8 surf(-100:2000,-100:2000),bed(-100:2000,-100:2000),melt(-100:2000,-100:2000)
Real*8 :: x,y,s1,b1,b2,u1,grid,m1,melta,wl,UC,z1
Real*8 :: box,xo(3,1000000),b,maxx,maxy,maxz,minx,miny,minz,SCL
Integer :: l,ND,ip,i,j,NDLC,NDRC,NDFC,NDBC,YN,NDFRC,NDBRC,NDFLC,NDBLC
Integer :: myid,ntasks,N1,N2,xk,yk

!Open(300,file='mass.dat',STATUS='OLD')
Open(510+myid,file='NODFIL2'//na(myid))

 12    FORMAT(2I8,' ',2F14.7)
 13    FORMAT(4F14.7)
 14    FORMAT(3F14.7)


!ALLOCATE (surf(-100:1000,-100:1000))
!ALLOCATE (bed(-100:1000,-100:1000))

DO i=-100,2000
DO j=-100,2000
surf(i,j)=-100.0
ENDDO
ENDDO

DO i=-100,2000
DO j=-100,2000
bed(i,j)=1000.0
ENDDO
ENDDO


Open(400,file='mass3.dat',STATUS='OLD')
READ(400,*) N2
DO I=1,N2
!READ(400,*) x,y,s1,b1,b2,m1
READ(400,*) x,y,s1,b1,b2,z1
y=y-7200.0
xk=INT(x/grid)
yk=INT(y/grid)
! if (b1.ne.0.0) then
! if (xk.ge.-100.and.yk.ge.-100) bed(xk,yk)=b1+250.0*sin(0.005*y)**8*cos(0.015*y)**2
! else
if (xk.ge.-100.and.yk.ge.-100) bed(xk,yk)=b1
! endif
if (xk.ge.-100.and.yk.ge.-100) surf(xk,yk)=s1
if (xk.ge.-100.and.yk.ge.-100) melt(xk,yk)=0.0
ENDDO
close(400)

!DO I=1,N1
!READ(300,*) x,y,s1,b1
!xk=INT(x/5.0)
!yk=INT(y/5.0)
!surf(xk,yk)=s1
!ENDDO

maxx=0.0
maxy=0.0
maxz=0.0
minx=1.0e+08
miny=1.0e+08
minz=1.0e+08

box=2.0d0**(2.0d0/3.0d0)*Dfloat(l) ! box size equal to fcc ground state

Call Initializefcc(box,l,xo,ip,myid,maxx,maxy,maxz,minx,miny,minz,ntasks,SCL,YN,surf,bed,melt,grid,wl,UC)

Call DT(ip,xo,ND,NDLC,NDRC,NDFC,NDBC,NDFRC,NDBRC,NDFLC,NDBLC,myid,ntasks,SCL,YN)

!DEALLOCATE (surf)
!DEALLOCATE (bed)

CLOSE(510+myid)
CLOSE(400)

End Subroutine

!---------------------------------------------------------------!


Subroutine Initializefcc(box,l,xo,ip,myid,maxx,maxy,maxz,minx,miny,minz,ntasks,SCL,YN,surf,bed,melt,grid,wl,UC)
Implicit None
INCLUDE 'na90.dat'
!Real*8,ALLOCATABLE :: surf(:,:),bed(:,:)
Real*8 surf(-100:2000,-100:2000),bed(-100:2000,-100:2000),melt(-100:2000,-100:2000)
Real*8 xo(3,1000000),b,x0(3,4),box,maxx,maxy,maxz,minx,miny,minz,SCL
!Real*8 z,surf(-100:3000,-100:3000),bed(-100:3000,-100:3000)
Real*8 z,x,y,sint,bint,mint,grid,wl,lc,UC,UCV
Integer i,j,k,k1,k2,l,ip,myid,ntasks,m,YN,xk,yk
! Setting up the four atom positions in one box

!Open(1510+myid,file='tt'//na(myid))
11    FORMAT(2I8,' ',2F14.7)
13    FORMAT(4F14.7)

b=SCL*box/DFloat(l)  ! the size of the unit cell is the box length divided by l
x0(:,:)=b/2.0d0; x0(:,1)=0.0d0; x0(3,2)=0.0d0; x0(2,3)=0.0d0; x0(1,4)=0.0d0 

ip=0
m=MOD(myid,ntasks/YN)
Do i=1+8*m,8+8*m
      Do j=(myid/(ntasks/YN))*8+1,(myid/(ntasks/YN)+1)*8
	Do k=-2,l
	  Do k1=1,4
             x=(x0(1,k1) + Float(i-1)*b)/grid
             y=(x0(2,k1) + Float(j-1)*b)/grid
             xk=INT((x0(1,k1) + Float(i-1)*b)/grid)
             yk=INT((x0(2,k1) + Float(j-1)*b)/grid)
             Call BIPINT(x-xk,y-yk,bed(xk,yk),bed(xk,yk+1),bed(xk+1,yk),bed(xk+1,yk+1),bint)
             Call BIPINT(x-xk,y-yk,surf(xk,yk),surf(xk,yk+1),surf(xk+1,yk),surf(xk+1,yk+1),sint)
             Call BIPINT(x-xk,y-yk,melt(xk,yk),melt(xk,yk+1),melt(xk+1,yk),melt(xk+1,yk+1),mint)
!             bint=bed(xk,yk)+(x-xk)*(bed(xk+1,yk)-bed(xk,yk))+(y-yk)*(bed(xk,yk+1)-bed(xk,yk))
!             sint=surf(xk,yk)+(x-xk)*(surf(xk+1,yk)-surf(xk,yk))+(y-yk)*(surf(xk,yk+1)-surf(xk,yk))
              z=x0(3,k1) + Float(k-1)*b
              y=x0(2,k1) + Float(j-1)*b
              x=x0(1,k1) + Float(i-1)*b
!             write(1510+myid,13) 40.0*x,40.0*y,bint,sint
!             If (bed(xk,yk).ne.0.0.and.bed(xk,yk+1).ne.0.0.and.bed(xk+1,yk).ne.0.0.and.bed(xk+1,yk+1).ne.0.0) then
!             If ((z.ge.bint+mint.or.z.ge.wl).and.z.le.sint.and.(sint-bint).gt.SCL) Then
             If (z.ge.bint+mint.and.z.lt.sint) then
!             lc=4420.0+1.5e-04*(x-3300.0)**2+0.42*exp((x-3700.0)/2.0e+02)
!             UCV=lc-1500.0*exp(-(x-3500.0)**2/50000.0)
!             If (y.lt.lc-UC.or.y.gt.lc.or.z.gt.bint+3.0*SQRT(y-(lc-UC)).or.z.ge.WL-20.0) then
!             If (y.lt.lc-UC.or.z.gt.bint+3.0*SQRT(y-(lc-UC)).or.z.ge.WL-40.0) then
!             If (y.lt.UCV.or.(z.gt.bint+3.0*sqrt(y-UCV)).or.z.ge.WL-40.0) then
             ip=ip+1
	     xo(1,ip) = x0(1,k1) + Float(i-1)*b	
	     xo(2,ip) = x0(2,k1) + Float(j-1)*b
	     xo(3,ip) = x0(3,k1) + Float(k-1)*b
             if (xo(1,ip).gt.maxx) maxx=xo(1,ip)
             if (xo(2,ip).gt.maxy) maxy=xo(2,ip)
             if (xo(3,ip).gt.maxz) maxz=xo(3,ip)
             if (xo(1,ip).lt.minx) minx=xo(1,ip)
             if (xo(2,ip).lt.miny) miny=xo(2,ip)
             if (xo(3,ip).lt.minz) minz=xo(3,ip)
!             EndIF
!             EndIF
             EndIF
!             EndIF
	  EndDo
	EndDo
      EndDo
End Do

 12    FORMAT(I8,' ',4F14.7)

Do i=1,ip
	Write(510+myid,12) i,xo(:,i),1.0
	End Do 
write(*,*) 'NN=',ip,myid

!DEALLOCATE (surf)
!DEALLOCATE (bed)

End Subroutine

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

