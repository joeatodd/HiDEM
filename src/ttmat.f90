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

! TTMAT represents the bond orientation/rotation matrix.

SUBROUTINE TTMAT(X1,Y1,Z1,X2,Y2,Z2,RY,TT,Transp)

  USE TypeDefs

  IMPLICIT NONE

  REAL(KIND=dp) :: TT(12,12),DX,DY,DZ,DX2,DY2,DZ2
  REAL(KIND=dp) :: X1,Y1,Z1,X2,Y2,Z2,L,LP
  INTEGER :: RY
  LOGICAL :: Transp

  DX=(X2-X1)
  DY=(Y2-Y1)
  DZ=(Z2-Z1)
  DX2=(X2-X1)**2
  DY2=(Y2-Y1)**2
  DZ2=(Z2-Z1)**2
  L=SQRT(DX2+DY2+DZ2)
  LP=SQRT(DX2+DY2)
  
  ! Need to guard against rare /0 here
  ! LP == 0 implies DX == DY == 0
  ! L == 0 implies LP == DX == DY == DZ == 0

  IF(L == 0) THEN !Two particles are in the same position

    TT = 0.0

  ELSE IF(LP == 0) THEN !Two particles are in same XY position

    TT(1,1)=0.0
    TT(2,1)=0.0
    TT(3,1)=DZ/L
    TT(1,2)=0.0
    TT(2,2)=0.0
    TT(3,2)=0.0
    TT(1,3)=0.0
    TT(2,3)=0.0
    TT(3,3)=0.0

  ELSE !Normal case

    TT(1,1)=DX/L
    TT(2,1)=DY/L
    TT(3,1)=DZ/L
    TT(1,2)=-DY/LP
    TT(2,2)=DX/LP
    TT(3,2)=0.0
    TT(1,3)=-(DZ*DX)/(L*LP)
    TT(2,3)=-(DZ*DY)/(L*LP)
    TT(3,3)=(DX2+DY2)/(L*LP)

  END IF


  TT(4,4)=TT(1,1)
  TT(5,4)=TT(2,1)
  TT(6,4)=TT(3,1)
  TT(4,5)=TT(1,2)
  TT(5,5)=TT(2,2)
  TT(6,5)=TT(3,2)
  TT(4,6)=TT(1,3)
  TT(5,6)=TT(2,3)
  TT(6,6)=TT(3,3)

  TT(7,7)=TT(1,1)
  TT(8,7)=TT(2,1)
  TT(9,7)=TT(3,1)
  TT(7,8)=TT(1,2)
  TT(8,8)=TT(2,2)
  TT(9,8)=TT(3,2)
  TT(7,9)=TT(1,3)
  TT(8,9)=TT(2,3)
  TT(9,9)=TT(3,3)

  TT(10,10)=TT(1,1)
  TT(11,10)=TT(2,1)
  TT(12,10)=TT(3,1)
  TT(10,11)=TT(1,2)
  TT(11,11)=TT(2,2)
  TT(12,11)=TT(3,2)
  TT(10,12)=TT(1,3)
  TT(11,12)=TT(2,3)
  TT(12,12)=TT(3,3)

  IF(Transp) THEN
    TT = TRANSPOSE(TT)
  END IF

  RETURN

END SUBROUTINE
