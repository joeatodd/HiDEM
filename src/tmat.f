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

	SUBROUTINE TMAT(X1,Y1,Z1,X2,Y2,Z2,RY,T)

	REAL*8 T(12,12),DX,DY,DZ,DX2,DY2,DZ2
	REAL*8 X1,Y1,Z1,X2,Y2,Z2,L,LP
	INTEGER RY
c        COMMON /BT/ T


	DX=(X2-X1)
	DY=(Y2-Y1)
	DZ=(Z2-Z1)
	DX2=(X2-X1)**2
	DY2=(Y2-Y1)**2
	DZ2=(Z2-Z1)**2
	L=SQRT(DX2+DY2+DZ2)
	LP=SQRT(DX2+DY2)

c	IF (RY.EQ.1) THEN
c	DO 20 J=1,12
c	DO 10 I=1,12
c	T(I,J)=0.0
c 10	CONTINUE
c 20	CONTINUE
c	ENDIF

	T(1,1)=DX/L
	T(1,2)=DY/L
	T(1,3)=DZ/L
	T(2,1)=-DY/LP
	T(2,2)=DX/LP
	T(2,3)=0.0
	T(3,1)=-(DZ*DX)/(L*LP)
	T(3,2)=-(DZ*DY)/(L*LP)
	T(3,3)=(DX2+DY2)/(L*LP)

	T(4,4)=DX/L
	T(4,5)=DY/L
	T(4,6)=DZ/L
	T(5,4)=-DY/LP
	T(5,5)=DX/LP
	T(5,6)=0.0
	T(6,4)=-(DZ*DX)/(L*LP)
	T(6,5)=-(DZ*DY)/(L*LP)
	T(6,6)=(DX2+DY2)/(L*LP)

	T(7,7)=DX/L
	T(7,8)=DY/L
	T(7,9)=DZ/L
	T(8,7)=-DY/LP
	T(8,8)=DX/LP
	T(8,9)=0.0
	T(9,7)=-(DZ*DX)/(L*LP)
	T(9,8)=-(DZ*DY)/(L*LP)
	T(9,9)=(DX2+DY2)/(L*LP)

	T(10,10)=DX/L
	T(10,11)=DY/L
	T(10,12)=DZ/L
	T(11,10)=-DY/LP
	T(11,11)=DX/LP
	T(11,12)=0.0
	T(12,10)=-(DZ*DX)/(L*LP)
	T(12,11)=-(DZ*DY)/(L*LP)
	T(12,12)=(DX2+DY2)/(L*LP)

	RETURN
	END
