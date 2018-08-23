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
     T,IS,DT,WE,EFC,FXF,FXC,NDL,LNN,SCL)

        USE TypeDefs
        USE Utils

	IMPLICIT NONE
        include 'mpif.h'
	REAL*8 X1,X2,Y1,Y2,Z1,Z2
	REAL*8 T1,T2
	REAL*8, ALLOCATABLE :: EFC(:)
	REAL*8 SX,SY,SZ,SUM,T,WE(NOMA),L0
	REAL*8 DDEL,DWE,OWE,DT,ESUM,LNN
	REAL*8 LS,LS2,DEL,SCL
	INTEGER ierr,FXC
        INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),comm
        INTEGER, ALLOCATABLE :: FXF(:,:)
	REAL*8 RC,RCX,RCY,RCZ,FRX(NOMA),FRY(NOMA),FRZ(NOMA)
	INTEGER NTOT,I,N1,N2,IS,NN
        TYPE(UT_t) :: UT
        TYPE(NRXF_t) :: NRXF
        TYPE(NTOT_t) :: ND
        TYPE(FXF_t) :: NDL

        ALLOCATE(FXF(2,NN*12))
        FXF = 0 !particle proximity info

        FRX = 0.0
	FRY=0.0
	FRZ=0.0
	WE=0.0
        FXC = 0
        FXF = 0

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

 !TODO - when this was parallel, FRX,Y,Z were only saved for *our* nodes (not other parts)
 !       does this matter??

 !TODO - expand FXF arr if necessary
        IF(FXC+1 > SIZE(FXF,2)) CALL ExpandIntArray(FXF)
 
	 IF (RC.LT.LNN) THEN
  	 FRX(N1)=FRX(N1)+EFC(N1)*(LNN-RC)**1.5*RCX
	 FRY(N1)=FRY(N1)+EFC(N1)*(LNN-RC)**1.5*RCY
	 FRZ(N1)=FRZ(N1)+EFC(N1)*(LNN-RC)**1.5*RCZ
	 FRX(N2)=FRX(N2)-EFC(N2)*(LNN-RC)**1.5*RCX
	 FRY(N2)=FRY(N2)-EFC(N2)*(LNN-RC)**1.5*RCY
	 FRZ(N2)=FRZ(N2)-EFC(N2)*(LNN-RC)**1.5*RCZ
	 WE(N2)=WE(N2)+0.4*EFC(N2)*(LNN-RC)**2.5
	 FXC=FXC+1
	 FXF(1,FXC)=N1 
	 FXF(2,FXC)=N2 
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

	RETURN

END SUBROUTINE CIRC






