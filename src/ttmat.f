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

	SUBROUTINE TTMAT(X1,Y1,Z1,X2,Y2,Z2,RY,TT)

	REAL*8 TT(12,12),DX,DY,DZ,DX2,DY2,DZ2
	REAL*8 X1,Y1,Z1,X2,Y2,Z2,L,LP
	INTEGER RY
c        COMMON /BTT/ TT

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
c	TT(I,J)=0.0
c 10	CONTINUE
c 20	CONTINUE
c	ENDIF

	TT(1,1)=DX/L
	TT(2,1)=DY/L
	TT(3,1)=DZ/L
	TT(1,2)=-DY/LP
	TT(2,2)=DX/LP
	TT(3,2)=0.0
	TT(1,3)=-(DZ*DX)/(L*LP)
	TT(2,3)=-(DZ*DY)/(L*LP)
	TT(3,3)=(DX2+DY2)/(L*LP)

	TT(4,4)=DX/L
	TT(5,4)=DY/L
	TT(6,4)=DZ/L
	TT(4,5)=-DY/LP
	TT(5,5)=DX/LP
	TT(6,5)=0.0
	TT(4,6)=-(DZ*DX)/(L*LP)
	TT(5,6)=-(DZ*DY)/(L*LP)
	TT(6,6)=(DX2+DY2)/(L*LP)

	TT(7,7)=DX/L
	TT(8,7)=DY/L
	TT(9,7)=DZ/L
	TT(7,8)=-DY/LP
	TT(8,8)=DX/LP
	TT(9,8)=0.0
	TT(7,9)=-(DZ*DX)/(L*LP)
	TT(8,9)=-(DZ*DY)/(L*LP)
	TT(9,9)=(DX2+DY2)/(L*LP)

	TT(10,10)=DX/L
	TT(11,10)=DY/L
	TT(12,10)=DZ/L
	TT(10,11)=-DY/LP
	TT(11,11)=DX/LP
	TT(12,11)=0.0
	TT(10,12)=-(DZ*DX)/(L*LP)
	TT(11,12)=-(DZ*DY)/(L*LP)
	TT(12,12)=(DX2+DY2)/(L*LP)

	RETURN
	END
