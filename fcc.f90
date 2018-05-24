! *************************************************************************
! *  HiDEM, A Discrete Element Model for Fracture Simulation
! *  Copyright (C) 24th May 2018 - Jan Åström
! *
! *  This program is free software: you can redistribute it and/or modify
! *  it under the terms of the GNU General Public License as published by
! *  the Free Software Foundation, either version 3 of the License, or
! *  (at your option) any later version.
! *
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
! *************************************************************************

!Not used

Subroutine FIBG3(l,ip,ND,NDLC,NDRC,NDFC,NDBC,NDFRC,NDBRC,NDFLC,NDBLC,myid,maxx,maxy,maxz,minx,miny,minz,ntasks,SCL,YN)

Implicit none
INCLUDE 'na90.dat'
Real*8 :: box,xo(3,10000000),b,maxx,maxy,maxz,minx,miny,minz,SCL
Integer :: l,ND,ip,i,NDLC,NDRC,NDFC,NDBC,YN,NDFRC,NDBRC,NDFLC,NDBLC
Integer :: myid,ntasks

Open(510+myid, file='NODFIL2'//na(myid))

maxx=0.0
maxy=0.0
maxz=0.0
minx=1.0e+08
miny=1.0e+08
minz=1.0e+08

box=2.0d0**(2.0d0/3.0d0)*Dfloat(l) ! box size equal to fcc ground state

Call Initializefcc(box,l,xo,ip,myid,maxx,maxy,maxz,minx,miny,minz,ntasks,SCL,YN)

Call DT(ip,xo,ND,NDLC,NDRC,NDFC,NDBC,NDFRC,NDBRC,NDFLC,NDBLC,myid,ntasks,SCL,YN)

!write(*,*) 'NLINKS=',ND+NDLC+NDRC
!write(*,*) 
!write(*,*) 'X',maxx
!write(*,*) 'Y',maxy
!write(*,*) 'Z',maxz

CLOSE(510+myid)

End Subroutine
!End Program fcc

!---------------------------------------------------------------!
! Subroutine to set up fcc lattice
Subroutine Initializefcc(box,l,xo,ip,myid,maxx,maxy,maxz,minx,miny,minz,ntasks,SCL,YN)
Implicit None
Real*8 xo(3,10000000),b,x0(3,4),box,maxx,maxy,maxz,minx,miny,minz,SCL
Real*8 xi(3,10000000)
Integer i,j,k,k1,k2,l,ip,myid,ntasks,m,YN

! Setting up the four atom positions in one box
b=SCL*box/DFloat(l)  ! the size of the unit cell is the box length divided by l
x0(:,:)=b/2.0d0; x0(:,1)=0.0d0; x0(3,2)=0.0d0; x0(2,3)=0.0d0; x0(1,4)=0.0d0 

ip=0
m=MOD(myid,ntasks/YN)
Do i=1+10*m,10+10*m
      Do j=(myid/(ntasks/YN))*10+1,(myid/(ntasks/YN)+1)*10
	Do k=1,l
	  Do k1=1,4
                If (((x0(2,k1)+Float(j-1)*b<maxy-10.0).or.(x0(3,k1)+Float(k-1)*b>5.0)).or.myid.lt.(YN-1)*ntasks/YN) Then
		ip=ip+1
		xi(1,ip) = x0(1,k1) + Float(i-1)*b	
		xi(2,ip) = x0(2,k1) + Float(j-1)*b
		xi(3,ip) = x0(3,k1) + Float(k-1)*b
                if (xi(1,ip).gt.maxx) maxx=xi(1,ip)
                if (xi(2,ip).gt.maxy) maxy=xi(2,ip)
                if (xi(3,ip).gt.maxz) maxz=xi(3,ip)
                if (xi(1,ip).lt.minx) minx=xi(1,ip)
                if (xi(2,ip).lt.miny) miny=xi(2,ip)
                if (xi(3,ip).lt.minz) minz=xi(3,ip)
                EndIF
	  EndDo
	EndDo
      EndDo
End Do

ip=0
m=MOD(myid,ntasks/YN)
Do i=1+10*m,10+10*m
      Do j=(myid/(ntasks/YN))*10+1,(myid/(ntasks/YN)+1)*10
	Do k=1,l
	  Do k1=1,4
                If (((x0(2,k1)+Float(j-1)*b<maxy-10.0).or.(x0(3,k1)+Float(k-1)*b>5.0)).or.myid.lt.(YN-1)*ntasks/YN) Then
		ip=ip+1
		xo(1,ip) = x0(1,k1) + Float(i-1)*b	
		xo(2,ip) = x0(2,k1) + Float(j-1)*b
		xo(3,ip) = x0(3,k1) + Float(k-1)*b
                EndIF
	  EndDo
	EndDo
      EndDo
End Do

 12    FORMAT(I8,' ',4F14.7)

Do i=1,ip
	Write(510+myid,12) i,xo(:,i),1.0
	End Do 
write(*,*) 'NN=',ip,myid
End Subroutine
!-----------------------------------------------------                                                        

Subroutine BIPINT(x,y,f11,f21,f12,f22,fint)
Implicit none
Real*8 :: x,y,f11,f12,f21,f22,fint
fint=f11*(1.0-x)*(1.0-y)+f21*x*(1.0-y)+f12*(1.0-x)*y+f22*x*y
End Subroutine

!-----------------------------------------------------                                                        

Subroutine BIPINTN(x,y,f11,f21,f12,f22,dix,diy,diz)
Implicit none
Real*8 :: x,y,f11,f12,f21,f22,dix,diy,diz
REAL*8 norm
dix=(-f11*(1.0-y)+f21*(1.0-y)-f12*y+f22*y)/40.0
diy=(-f11*(1.0-x)-f21*x+f12*(1.0-x)+f22*x)/40.0
diz=1.0
norm=SQRT(dix**2.0+diy**2.0+1.0)
dix=-dix/norm
diy=-diy/norm
diz=diz/norm
!if (dix.gt.0.3) dix=0.3
!if (diy.gt.0.3) diy=0.3
!if (dix.lt.-0.3) dix=-0.3
!if (diy.lt.-0.3) diy=-0.3
End Subroutine
