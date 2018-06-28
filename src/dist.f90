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

	SUBROUTINE DIST(NNA,UT,ND,NRXF,NDL,myid,ntasks,SCL,PNN,YN)

        USE TypeDefs

	IMPLICIT NONE
        include 'mpif.h'
	REAL*8 DX1,DX2,DY1,DY2,DZ1,DZ2
	REAL*8 X1,X2,Y1,Y2,Z1,Z2
	REAL*8 RC,SCL,RT
	INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),comm
	INTEGER myid,ntasks,ierr,YN
	INTEGER NNA,I,J,PNN(0:5000)
        TYPE(UT_t) :: UT
        TYPE(NRXF_t) :: NRXF
        TYPE(NTOT_t) :: ND
        TYPE(FXF_t) :: NDL
!	OPEN(UNIT=10,FILE='TSR',STATUS='UNKNOWN')
!	OPEN(UNIT=11,FILE='TSL',STATUS='UNKNOWN')

	ND%M=0
	ND%R=0
	ND%B=0
	ND%F=0
	ND%BR=0
	ND%FR=0
	ND%BL=0
	ND%FL=0
	ND%L=0
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
	  ND%M=ND%M+1
	  NDL%M(1,ND%M)=I
	  NDL%M(2,ND%M)=J
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
	  ND%L=ND%L+1
	  NDL%L(1,ND%L)=I
	  NDL%L(2,ND%L)=J
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
	  ND%R=ND%R+1
	  NDL%R(1,ND%R)=I
	  NDL%R(2,ND%R)=J
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
	  ND%F=ND%F+1
	  NDL%F(1,ND%F)=I
	  NDL%F(2,ND%F)=J
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
	  ND%B=ND%B+1
	  NDL%B(1,ND%B)=I
	  NDL%B(2,ND%B)=J
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
	  ND%BL=ND%BL+1
	  NDL%BL(1,ND%BL)=I
	  NDL%BL(2,ND%BL)=J
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
	  ND%BR=ND%BR+1
	  NDL%BR(1,ND%BR)=I
	  NDL%BR(2,ND%BR)=J
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
	  ND%FL=ND%FL+1
	  NDL%FL(1,ND%FL)=I
	  NDL%FL(2,ND%FL)=J
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
	  ND%FR=ND%FR+1
	  NDL%FR(1,ND%FR)=I
	  NDL%FR(2,ND%FR)=J
 	  ENDIF
	 END DO
	 END DO
	 END IF

	RETURN
	END





