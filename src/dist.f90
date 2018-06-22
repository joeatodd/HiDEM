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

	SUBROUTINE DIST(NNA,UT,ND,NCL,NDR,NRXF, &
     	NDL,NDLL,NDLR,NDLF,NDLB,NDF,NDB,myid,ntasks,SCL,PNN,YN, &
     	NDLFL,NDLFR,NDLBR,NDLBL,NDFL,NDFR,NDBR,NDBL)

        USE TypeDefs

	IMPLICIT NONE
        include 'mpif.h'
	REAL*8 DX1,DX2,DY1,DY2,DZ1,DZ2
	REAL*8 X1,X2,Y1,Y2,Z1,Z2
	REAL*8 RC,SCL,RT
	INTEGER NCL,NDR,NDF,NDB
	INTEGER NDFL,NDBL
	INTEGER NDFR,NDBR
	INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),comm
	INTEGER myid,ntasks,ierr,YN
	INTEGER NNA,ND,I,J,PNN(0:5000)
	INTEGER NDL(2,NODC),NDLR(2,NODC),NDLL(2,NODC)
	INTEGER NDLF(2,NODC),NDLB(2,NODC)
	INTEGER NDLFL(2,NODC),NDLBL(2,NODC)
	INTEGER NDLFR(2,NODC),NDLBR(2,NODC)
        TYPE(UT_t) :: UT
        TYPE(NRXF_t) :: NRXF

!	OPEN(UNIT=10,FILE='TSR',STATUS='UNKNOWN')
!	OPEN(UNIT=11,FILE='TSL',STATUS='UNKNOWN')

	ND=0
	NDR=0
	NDB=0
	NDF=0
	NDBR=0
	NDFR=0
	NDBL=0
	NDFL=0
	NCL=0
	RT=SCL*SCL*3.5

	DO I=1,NNA-1
         DX1=UT%M(6*I-5)
         DY1=UT%M(6*I-4)
         DZ1=UT%M(6*I-3)
 	 X1=NRXF%M(1,I)+DX1
	 Y1=NRXF%M(2,I)+DY1
	 Z1=NRXF%M(3,I)+DZ1
	 DO J=I+1,NNA
          DX2=UT%M(6*J-5)
          DY2=UT%M(6*J-4)
          DZ2=UT%M(6*J-3)
	  X2=NRXF%M(1,J)+DX2
	  Y2=NRXF%M(2,J)+DY2
	  Z2=NRXF%M(3,J)+DZ2
	  RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  IF (RC.LT.RT) THEN
	  ND=ND+1
	  NDL(1,ND)=I
	  NDL(2,ND)=J
 	  ENDIF
	 END DO
        END DO

      IF (MOD(myid,ntasks/YN).ne.0) THEN
	DO I=1,PNN(myid-1)
          DX1=UT%L(6*I-5)
          DY1=UT%L(6*I-4)
          DZ1=UT%L(6*I-3)
 	  X1=NRXF%L(1,I)+DX1
	  Y1=NRXF%L(2,I)+DY1
	  Z1=NRXF%L(3,I)+DZ1
	 DO J=1,NNA
          DX2=UT%M(6*J-5)
          DY2=UT%M(6*J-4)
          DZ2=UT%M(6*J-3)
	  X2=NRXF%M(1,J)+DX2
	  Y2=NRXF%M(2,J)+DY2
	  Z2=NRXF%M(3,J)+DZ2
	  RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  IF (RC.LT.RT) THEN
	  NCL=NCL+1
	  NDLL(1,NCL)=I
	  NDLL(2,NCL)=J
 	  ENDIF
	 END DO
	 END DO
	 END IF

      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
	DO I=1,PNN(myid+1)
          DX1=UT%R(6*I-5)
          DY1=UT%R(6*I-4)
          DZ1=UT%R(6*I-3)
 	  X1=NRXF%R(1,I)+DX1
	  Y1=NRXF%R(2,I)+DY1
	  Z1=NRXF%R(3,I)+DZ1
	 DO J=1,NNA
          DX2=UT%M(6*J-5)
          DY2=UT%M(6*J-4)
          DZ2=UT%M(6*J-3)
	  X2=NRXF%M(1,J)+DX2
	  Y2=NRXF%M(2,J)+DY2
	  Z2=NRXF%M(3,J)+DZ2
	  RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  IF (RC.LT.RT) THEN
	  NDR=NDR+1
	  NDLR(1,NDR)=I
	  NDLR(2,NDR)=J
 	  ENDIF
	 END DO
	 END DO
	 END IF

      IF (myid.lt.(YN-1)*ntasks/YN) THEN
	DO I=1,PNN(myid+ntasks/YN)
          DX1=UT%F(6*I-5)
          DY1=UT%F(6*I-4)
          DZ1=UT%F(6*I-3)
 	  X1=NRXF%F(1,I)+DX1
	  Y1=NRXF%F(2,I)+DY1
	  Z1=NRXF%F(3,I)+DZ1
	 DO J=1,NNA
          DX2=UT%M(6*J-5)
          DY2=UT%M(6*J-4)
          DZ2=UT%M(6*J-3)
	  X2=NRXF%M(1,J)+DX2
	  Y2=NRXF%M(2,J)+DY2
	  Z2=NRXF%M(3,J)+DZ2
	  RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  IF (RC.LT.RT) THEN
	  NDF=NDF+1
	  NDLF(1,NDF)=I
	  NDLF(2,NDF)=J
 	  ENDIF
	 END DO
	 END DO
	 END IF


      IF (myid.ge.ntasks/YN) THEN
	DO I=1,PNN(myid-ntasks/YN)
          DX1=UT%B(6*I-5)
          DY1=UT%B(6*I-4)
          DZ1=UT%B(6*I-3)
 	  X1=NRXF%B(1,I)+DX1
	  Y1=NRXF%B(2,I)+DY1
	  Z1=NRXF%B(3,I)+DZ1
	 DO J=1,NNA
          DX2=UT%M(6*J-5)
          DY2=UT%M(6*J-4)
          DZ2=UT%M(6*J-3)
	  X2=NRXF%M(1,J)+DX2
	  Y2=NRXF%M(2,J)+DY2
	  Z2=NRXF%M(3,J)+DZ2
	  RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  IF (RC.LT.RT) THEN
	  NDB=NDB+1
	  NDLB(1,NDB)=I
	  NDLB(2,NDB)=J
 	  ENDIF
	 END DO
	 END DO
	 END IF


      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
	DO I=1,PNN(myid-ntasks/YN-1)
          DX1=UT%BL(6*I-5)
          DY1=UT%BL(6*I-4)
          DZ1=UT%BL(6*I-3)
 	  X1=NRXF%BL(1,I)+DX1
	  Y1=NRXF%BL(2,I)+DY1
	  Z1=NRXF%BL(3,I)+DZ1
	 DO J=1,NNA
          DX2=UT%M(6*J-5)
          DY2=UT%M(6*J-4)
          DZ2=UT%M(6*J-3)
	  X2=NRXF%M(1,J)+DX2
	  Y2=NRXF%M(2,J)+DY2
	  Z2=NRXF%M(3,J)+DZ2
	  RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  IF (RC.LT.RT) THEN
	  NDBL=NDBL+1
	  NDLBL(1,NDBL)=I
	  NDLBL(2,NDBL)=J
 	  ENDIF
	 END DO
	 END DO
	 END IF

      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
	DO I=1,PNN(myid-ntasks/YN+1)
          DX1=UT%BR(6*I-5)
          DY1=UT%BR(6*I-4)
          DZ1=UT%BR(6*I-3)
 	  X1=NRXF%BR(1,I)+DX1
	  Y1=NRXF%BR(2,I)+DY1
	  Z1=NRXF%BR(3,I)+DZ1
	 DO J=1,NNA
          DX2=UT%M(6*J-5)
          DY2=UT%M(6*J-4)
          DZ2=UT%M(6*J-3)
	  X2=NRXF%M(1,J)+DX2
	  Y2=NRXF%M(2,J)+DY2
	  Z2=NRXF%M(3,J)+DZ2
	  RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  IF (RC.LT.RT) THEN
	  NDBR=NDBR+1
	  NDLBR(1,NDBR)=I
	  NDLBR(2,NDBR)=J
 	  ENDIF
	 END DO
	 END DO
	 END IF


      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
	DO I=1,PNN(myid+ntasks/YN-1)
          DX1=UT%FL(6*I-5)
          DY1=UT%FL(6*I-4)
          DZ1=UT%FL(6*I-3)
 	  X1=NRXF%FL(1,I)+DX1
	  Y1=NRXF%FL(2,I)+DY1
	  Z1=NRXF%FL(3,I)+DZ1
	 DO J=1,NNA
          DX2=UT%M(6*J-5)
          DY2=UT%M(6*J-4)
          DZ2=UT%M(6*J-3)
	  X2=NRXF%M(1,J)+DX2
	  Y2=NRXF%M(2,J)+DY2
	  Z2=NRXF%M(3,J)+DZ2
	  RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  IF (RC.LT.RT) THEN
	  NDFL=NDFL+1
	  NDLFL(1,NDFL)=I
	  NDLFL(2,NDFL)=J
 	  ENDIF
	 END DO
	 END DO
	 END IF

      IF (myid.lt.(YN-1)*ntasks/YN &
           .AND.MOD(myid,ntasks/YN).NE.ntasks/YN-1) THEN
	DO I=1,PNN(myid+ntasks/YN+1)
          DX1=UT%FR(6*I-5)
          DY1=UT%FR(6*I-4)
          DZ1=UT%FR(6*I-3)
 	  X1=NRXF%FR(1,I)+DX1
	  Y1=NRXF%FR(2,I)+DY1
	  Z1=NRXF%FR(3,I)+DZ1
	 DO J=1,NNA
          DX2=UT%M(6*J-5)
          DY2=UT%M(6*J-4)
          DZ2=UT%M(6*J-3)
	  X2=NRXF%M(1,J)+DX2
	  Y2=NRXF%M(2,J)+DY2
	  Z2=NRXF%M(3,J)+DZ2
	  RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  IF (RC.LT.RT) THEN
	  NDFR=NDFR+1
	  NDLFR(1,NDFR)=I
	  NDLFR(2,NDFR)=J
 	  ENDIF
	 END DO
	 END DO
	 END IF

	RETURN
	END





