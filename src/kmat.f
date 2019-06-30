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
!
! KMAT returns the bond stiffness matrix
	SUBROUTINE KMAT(E,G,W,L,RY,RK)

	REAL*8 RK(12,12)
	REAL*8 E,G,W,L,W4,W2,L2,L3
	INTEGER RY
c        COMMON /KM/ RK

	W4=W**4
	W2=W**2
	L2=L**2
	L3=L**3

c	IF (RY.EQ.1) THEN
c	DO 20 J=1,12
c	DO 10 I=1,12
c	RK(I,J)=0.0
c 10	CONTINUE
c 20	CONTINUE
c	ENDIF

	RK(1,1)=E*W2/L
	RK(2,2)=E*W4/(L3)
	RK(3,3)=E*W4/(L3)
	RK(4,4)=G*W4/(6*L)
	RK(5,5)=E*W4/(3*L)
	RK(6,6)=E*W4/(3*L)
	RK(7,7)=E*W2/L
	RK(8,8)=E*W4/(L3)
	RK(9,9)=E*W4/(L3)
	RK(10,10)=G*W4/(6*L)
	RK(11,11)=E*W4/(3*L)
	RK(12,12)=E*W4/(3*L)

	RK(1,7)=-E*W2/L
	RK(7,1)=RK(1,7)

	RK(2,6)=E*W4/(2*L2)
	RK(2,8)=-E*W4/(L3)
	RK(2,12)=E*W4/(2*L2)
	RK(6,2)=RK(2,6)
	RK(8,2)=RK(2,8)
	RK(12,2)=RK(2,12)

	RK(3,5)=-E*W4/(2*L2)
	RK(3,9)=-E*W4/(L3)
	RK(3,11)=-E*W4/(2*L2)
	RK(5,3)=RK(3,5)
	RK(9,3)=RK(3,9)
	RK(11,3)=RK(3,11)

	RK(4,10)=-G*W4/(6*L)
	RK(10,4)=RK(4,10)

	RK(5,9)=E*W4/(2*L2)
	RK(5,11)=E*W4/(6*L)
	RK(9,5)=RK(5,9)
	RK(11,5)=RK(5,11)

	RK(6,8)=-E*W4/(2*L2)
	RK(6,12)=E*W4/(6*L)
	RK(8,6)=RK(6,8)
	RK(12,6)=RK(6,12)

	RK(8,12)=-E*W4/(2*L2)
	RK(12,8)=RK(8,12)

	RK(9,11)=E*W4/(2*L2)
	RK(11,9)=RK(9,11)

	RETURN
	END
