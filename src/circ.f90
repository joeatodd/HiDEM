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

SUBROUTINE CIRC(ND,NN,NRXF,UT,FRX,FRY,FRZ, &
     T,IS,DT,WE,EFC,FXF,FXC,NDL,myid,ntasks,LNN,YN,SCL)

        USE TypeDefs

	IMPLICIT NONE
        include 'mpif.h'
	REAL*8 X1,X2,Y1,Y2,Z1,Z2
	REAL*8 T1,T2
	REAL*8 EFC(NOMA)
	REAL*8 SX,SY,SZ,SUM,T,WE(NOMA),L0
	REAL*8 DDEL,DWE,OWE,DT,ESUM,LNN
	REAL*8 LS,LS2,DEL,SCL
	INTEGER myid,ntasks,ierr,YN
        INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),comm
	REAL*8 RC,RCX,RCY,RCZ,FRX(NOMA),FRY(NOMA),FRZ(NOMA)
	INTEGER NTOT,I,N1,N2,IS,NN
        TYPE(UT_t) :: UT
        TYPE(NRXF_t) :: NRXF
        TYPE(NTOT_t) :: FXC,ND
        TYPE(FXF_t) :: FXF,NDL

	DO I=1,NOMA
	FRX(I)=0.0
	FRY(I)=0.0
	FRZ(I)=0.0
	WE(I)=0.0
	END DO

	FXC%M=0
	DO I=1,ND%M
	N1=NDL%M(1,I)
	N2=NDL%M(2,I)
	X1=NRXF%M(1,N1)+UT%M(6*N1-5)
	Y1=NRXF%M(2,N1)+UT%M(6*N1-4)
	Z1=NRXF%M(3,N1)+UT%M(6*N1-3)
	X2=NRXF%M(1,N2)+UT%M(6*N2-5)
	Y2=NRXF%M(2,N2)+UT%M(6*N2-4)
	Z2=NRXF%M(3,N2)+UT%M(6*N2-3)
	IF (ABS(X1-X2).LE.LNN.AND.ABS(Y1-Y2).LE.LNN.AND.ABS(Z1-Z2).LE.LNN) THEN
	RC=SQRT((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	RCX=(X1-X2)/RC
	RCY=(Y1-Y2)/RC
	RCZ=(Z1-Z2)/RC

	 IF (RC.LT.LNN) THEN
  	 FRX(N1)=FRX(N1)+EFC(N1)*(LNN-RC)**1.5*RCX
	 FRY(N1)=FRY(N1)+EFC(N1)*(LNN-RC)**1.5*RCY
	 FRZ(N1)=FRZ(N1)+EFC(N1)*(LNN-RC)**1.5*RCZ
	 FRX(N2)=FRX(N2)-EFC(N2)*(LNN-RC)**1.5*RCX
	 FRY(N2)=FRY(N2)-EFC(N2)*(LNN-RC)**1.5*RCY
	 FRZ(N2)=FRZ(N2)-EFC(N2)*(LNN-RC)**1.5*RCZ
	 WE(N2)=WE(N2)+0.4*EFC(N2)*(LNN-RC)**2.5
	 FXC%M=FXC%M+1
	 FXF%M(1,FXC%M)=N1 
	 FXF%M(2,FXC%M)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.04*SCL) THEN
         FRX(N1)=FRX(N1)+SCL**2.0*1.0e+04*(LNN-RC)*RCX
         FRY(N1)=FRY(N1)+SCL**2.0*1.0e+04*(LNN-RC)*RCY
         FRZ(N1)=FRZ(N1)+SCL**2.0*1.0e+04*(LNN-RC)*RCZ
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+04*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO

!        IF (myid.ne.0.and.myid.ne.(ntasks-1)/2+1) THEN
        IF (MOD(myid,ntasks/YN).ne.0) THEN
	FXC%L=0
	DO I=1,ND%L
	N1=NDL%L(1,I)
	N2=NDL%L(2,I)
	X1=NRXF%L(1,N1)+UT%L(6*N1-5)
	Y1=NRXF%L(2,N1)+UT%L(6*N1-4)
	Z1=NRXF%L(3,N1)+UT%L(6*N1-3)
	X2=NRXF%M(1,N2)+UT%M(6*N2-5)
	Y2=NRXF%M(2,N2)+UT%M(6*N2-4)
	Z2=NRXF%M(3,N2)+UT%M(6*N2-3)
	IF (ABS(X1-X2).LE.LNN.AND.ABS(Y1-Y2).LE.LNN.AND.ABS(Z1-Z2).LE.LNN) THEN
	RC=SQRT((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	RCX=(X1-X2)/RC
	RCY=(Y1-Y2)/RC
	RCZ=(Z1-Z2)/RC

	 IF (RC.LT.LNN) THEN
	 FRX(N2)=FRX(N2)-EFC(N2)*(LNN-RC)**1.5*RCX
	 FRY(N2)=FRY(N2)-EFC(N2)*(LNN-RC)**1.5*RCY
	 FRZ(N2)=FRZ(N2)-EFC(N2)*(LNN-RC)**1.5*RCZ
	 WE(N2)=WE(N2)+0.4*EFC(N2)*(LNN-RC)**2.5
	 FXC%L=FXC%L+1
	 FXF%L(1,FXC%L)=N1 
	 FXF%L(2,FXC%L)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.04*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+04*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

!        IF (myid.ne.ntasks-1.and.myid.ne.(ntasks-1)/2) THEN
        IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
	FXC%R=0
	DO I=1,ND%R
	N1=NDL%R(1,I)
	N2=NDL%R(2,I)
	X1=NRXF%R(1,N1)+UT%R(6*N1-5)
	Y1=NRXF%R(2,N1)+UT%R(6*N1-4)
	Z1=NRXF%R(3,N1)+UT%R(6*N1-3)
	X2=NRXF%M(1,N2)+UT%M(6*N2-5)
	Y2=NRXF%M(2,N2)+UT%M(6*N2-4)
	Z2=NRXF%M(3,N2)+UT%M(6*N2-3)
	IF (ABS(X1-X2).LE.LNN.AND.ABS(Y1-Y2).LE.LNN.AND.ABS(Z1-Z2).LE.LNN) THEN
	RC=SQRT((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	RCX=(X1-X2)/RC
	RCY=(Y1-Y2)/RC
	RCZ=(Z1-Z2)/RC

	 IF (RC.LT.LNN) THEN
! 	 write(*,*) N1,N2,RC,LNN
	 FRX(N2)=FRX(N2)-EFC(N2)*(LNN-RC)**1.5*RCX
	 FRY(N2)=FRY(N2)-EFC(N2)*(LNN-RC)**1.5*RCY
	 FRZ(N2)=FRZ(N2)-EFC(N2)*(LNN-RC)**1.5*RCZ
	 WE(N2)=WE(N2)+0.4*EFC(N2)*(LNN-RC)**2.5
	 FXC%R=FXC%R+1
	 FXF%R(1,FXC%R)=N1 
	 FXF%R(2,FXC%R)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.04*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+04*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

!------------------------------------------------

        IF (myid.lt.(YN-1)*ntasks/YN) THEN
	FXC%F=0
	DO I=1,ND%F
	N1=NDL%F(1,I)
	N2=NDL%F(2,I)
	X1=NRXF%F(1,N1)+UT%F(6*N1-5)
	Y1=NRXF%F(2,N1)+UT%F(6*N1-4)
	Z1=NRXF%F(3,N1)+UT%F(6*N1-3)
	X2=NRXF%M(1,N2)+UT%M(6*N2-5)
	Y2=NRXF%M(2,N2)+UT%M(6*N2-4)
	Z2=NRXF%M(3,N2)+UT%M(6*N2-3)
	IF (ABS(X1-X2).LE.LNN.AND.ABS(Y1-Y2).LE.LNN.AND.ABS(Z1-Z2).LE.LNN) THEN
	RC=SQRT((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	RCX=(X1-X2)/RC
	RCY=(Y1-Y2)/RC
	RCZ=(Z1-Z2)/RC

	 IF (RC.LT.LNN) THEN
	 FRX(N2)=FRX(N2)-EFC(N2)*(LNN-RC)**1.5*RCX
	 FRY(N2)=FRY(N2)-EFC(N2)*(LNN-RC)**1.5*RCY
	 FRZ(N2)=FRZ(N2)-EFC(N2)*(LNN-RC)**1.5*RCZ
	 WE(N2)=WE(N2)+0.4*EFC(N2)*(LNN-RC)**2.5
	 FXC%F=FXC%F+1
	 FXF%F(1,FXC%F)=N1 
	 FXF%F(2,FXC%F)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.04*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCZ
 	 WE(N2)=WE(N2)+SCL**2.0*0.5e+04*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
	FXC%FL=0
	DO I=1,ND%FL
	N1=NDL%FL(1,I)
	N2=NDL%FL(2,I)
	X1=NRXF%FL(1,N1)+UT%FL(6*N1-5)
	Y1=NRXF%FL(2,N1)+UT%FL(6*N1-4)
	Z1=NRXF%FL(3,N1)+UT%FL(6*N1-3)
	X2=NRXF%M(1,N2)+UT%M(6*N2-5)
	Y2=NRXF%M(2,N2)+UT%M(6*N2-4)
	Z2=NRXF%M(3,N2)+UT%M(6*N2-3)
	IF (ABS(X1-X2).LE.LNN.AND.ABS(Y1-Y2).LE.LNN.AND.ABS(Z1-Z2).LE.LNN) THEN
	RC=SQRT((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	RCX=(X1-X2)/RC
	RCY=(Y1-Y2)/RC
	RCZ=(Z1-Z2)/RC

	 IF (RC.LT.LNN) THEN
	 FRX(N2)=FRX(N2)-EFC(N2)*(LNN-RC)**1.5*RCX
	 FRY(N2)=FRY(N2)-EFC(N2)*(LNN-RC)**1.5*RCY
	 FRZ(N2)=FRZ(N2)-EFC(N2)*(LNN-RC)**1.5*RCZ
	 WE(N2)=WE(N2)+0.4*EFC(N2)*(LNN-RC)**2.5
	 FXC%FL=FXC%FL+1
	 FXF%FL(1,FXC%FL)=N1 
	 FXF%FL(2,FXC%FL)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.04*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+04*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN &
             .AND.MOD(myid,ntasks/YN).NE.ntasks/YN-1) THEN
	FXC%FR=0
	DO I=1,ND%FR
	N1=NDL%FR(1,I)
	N2=NDL%FR(2,I)
	X1=NRXF%FR(1,N1)+UT%FR(6*N1-5)
	Y1=NRXF%FR(2,N1)+UT%FR(6*N1-4)
	Z1=NRXF%FR(3,N1)+UT%FR(6*N1-3)
	X2=NRXF%M(1,N2)+UT%M(6*N2-5)
	Y2=NRXF%M(2,N2)+UT%M(6*N2-4)
	Z2=NRXF%M(3,N2)+UT%M(6*N2-3)
	IF (ABS(X1-X2).LE.LNN.AND.ABS(Y1-Y2).LE.LNN.AND.ABS(Z1-Z2).LE.LNN) THEN
	RC=SQRT((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	RCX=(X1-X2)/RC
	RCY=(Y1-Y2)/RC
	RCZ=(Z1-Z2)/RC

	 IF (RC.LT.LNN) THEN
	 FRX(N2)=FRX(N2)-EFC(N2)*(LNN-RC)**1.5*RCX
	 FRY(N2)=FRY(N2)-EFC(N2)*(LNN-RC)**1.5*RCY
	 FRZ(N2)=FRZ(N2)-EFC(N2)*(LNN-RC)**1.5*RCZ
	 WE(N2)=WE(N2)+0.4*EFC(N2)*(LNN-RC)**2.5
	 FXC%FR=FXC%FR+1
	 FXF%FR(1,FXC%FR)=N1 
	 FXF%FR(2,FXC%FR)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.04*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+04*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

        IF (myid.ge.ntasks/YN) THEN
	FXC%B=0
	DO I=1,ND%B
	N1=NDL%B(1,I)
	N2=NDL%B(2,I)
	X1=NRXF%B(1,N1)+UT%B(6*N1-5)
	Y1=NRXF%B(2,N1)+UT%B(6*N1-4)
	Z1=NRXF%B(3,N1)+UT%B(6*N1-3)
	X2=NRXF%M(1,N2)+UT%M(6*N2-5)
	Y2=NRXF%M(2,N2)+UT%M(6*N2-4)
	Z2=NRXF%M(3,N2)+UT%M(6*N2-3)
	IF (ABS(X1-X2).LE.LNN.AND.ABS(Y1-Y2).LE.LNN.AND.ABS(Z1-Z2).LE.LNN) THEN
	RC=SQRT((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	RCX=(X1-X2)/RC
	RCY=(Y1-Y2)/RC
	RCZ=(Z1-Z2)/RC

	 IF (RC.LT.LNN) THEN
	 FRX(N2)=FRX(N2)-EFC(N2)*(LNN-RC)**1.5*RCX
	 FRY(N2)=FRY(N2)-EFC(N2)*(LNN-RC)**1.5*RCY
	 FRZ(N2)=FRZ(N2)-EFC(N2)*(LNN-RC)**1.5*RCZ
	 WE(N2)=WE(N2)+0.4*EFC(N2)*(LNN-RC)**2.5
	 FXC%B=FXC%B+1
	 FXF%B(1,FXC%B)=N1 
	 FXF%B(2,FXC%B)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.04*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+04*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

        IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
	FXC%BL=0
	DO I=1,ND%BL
	N1=NDL%BL(1,I)
	N2=NDL%BL(2,I)
	X1=NRXF%BL(1,N1)+UT%BL(6*N1-5)
	Y1=NRXF%BL(2,N1)+UT%BL(6*N1-4)
	Z1=NRXF%BL(3,N1)+UT%BL(6*N1-3)
	X2=NRXF%M(1,N2)+UT%M(6*N2-5)
	Y2=NRXF%M(2,N2)+UT%M(6*N2-4)
	Z2=NRXF%M(3,N2)+UT%M(6*N2-3)
	IF (ABS(X1-X2).LE.LNN.AND.ABS(Y1-Y2).LE.LNN.AND.ABS(Z1-Z2).LE.LNN) THEN
	RC=SQRT((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	RCX=(X1-X2)/RC
	RCY=(Y1-Y2)/RC
	RCZ=(Z1-Z2)/RC

	 IF (RC.LT.LNN) THEN
	 FRX(N2)=FRX(N2)-EFC(N2)*(LNN-RC)**1.5*RCX
	 FRY(N2)=FRY(N2)-EFC(N2)*(LNN-RC)**1.5*RCY
	 FRZ(N2)=FRZ(N2)-EFC(N2)*(LNN-RC)**1.5*RCZ
	 WE(N2)=WE(N2)+0.4*EFC(N2)*(LNN-RC)**2.5
	 FXC%BL=FXC%BL+1
	 FXF%BL(1,FXC%BL)=N1 
	 FXF%BL(2,FXC%BL)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.04*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+04*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

        IF (myid.ge.ntasks/YN &
             .AND.MOD(myid,ntasks/YN).NE.ntasks/YN-1) THEN
	FXC%BR=0
	DO I=1,ND%BR
	N1=NDL%BR(1,I)
	N2=NDL%BR(2,I)
	X1=NRXF%BR(1,N1)+UT%BR(6*N1-5)
	Y1=NRXF%BR(2,N1)+UT%BR(6*N1-4)
	Z1=NRXF%BR(3,N1)+UT%BR(6*N1-3)
	X2=NRXF%M(1,N2)+UT%M(6*N2-5)
	Y2=NRXF%M(2,N2)+UT%M(6*N2-4)
	Z2=NRXF%M(3,N2)+UT%M(6*N2-3)
	IF (ABS(X1-X2).LE.LNN.AND.ABS(Y1-Y2).LE.LNN.AND.ABS(Z1-Z2).LE.LNN) THEN
	RC=SQRT((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	RCX=(X1-X2)/RC
	RCY=(Y1-Y2)/RC
	RCZ=(Z1-Z2)/RC

	 IF (RC.LT.LNN) THEN
	 FRX(N2)=FRX(N2)-EFC(N2)*(LNN-RC)**1.5*RCX
	 FRY(N2)=FRY(N2)-EFC(N2)*(LNN-RC)**1.5*RCY
	 FRZ(N2)=FRZ(N2)-EFC(N2)*(LNN-RC)**1.5*RCZ
	 WE(N2)=WE(N2)+0.4*EFC(N2)*(LNN-RC)**2.5
	 FXC%BR=FXC%BR+1
	 FXF%BR(1,FXC%BR)=N1 
	 FXF%BR(2,FXC%BR)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.04*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+04*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+04*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF


	RETURN

END SUBROUTINE CIRC






