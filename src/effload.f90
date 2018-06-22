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

	SUBROUTINE EFFLOAD(S,NTOT,NN,T,DT,M,JS,DMP,DMP2,UT,UTM,R,EN,RY, &
     FXF,FXC,VDP,DPE,EFS,NANS,NRXF,MFIL,CT,&
     myid,ntasks,FXL,FXFL,FXR,FXFR,L, &
     FXCF,FXFF,FXCB,FXFB,PNN,YN, &
     FXCFL,FXFFL,FXCBL,FXFBL, &
     FXCFR,FXFFR,FXCBR,FXFBR)

        USE INOUT
        USE TypeDefs

	IMPLICIT NONE
        include 'mpif.h'
        REAL*8 MFIL(NOMA),VDP(NOMA)
	REAL*8 A(NODM),C(NODM),F(NODM),D(NODM)
	REAL*8 LNN
	REAL*8 DUT(NODM),CT(NODC),EN(NODM),R(NODM)
	INTEGER FXF(2,NODC),FXFL(2,NODC),FXFR(2,NODC)
	INTEGER FXFF(2,NODC),FXFB(2,NODC),YN
        INTEGER FXFFL(2,NODC),FXFBL(2,NODC)
        INTEGER FXFFR(2,NODC),FXFBR(2,NODC)
	REAL*8 DPE,S,E
	REAL*8 T,DT,M,JS,L,ALF,DMP
	REAL*8 G,X1,Y1,Z1,X2,Y2,Z2,TT(12,12)
	REAL*8 DX1,DY1,DZ1,DX2,DY2,DZ2,DMP2
	REAL*8 DP,DP2
	INTEGER N,NL,NB,N1,N2,X,XL,XR,PNN(0:5000)
	INTEGER I,J,NN,RY,FXC,FXL,FXR,FXCF,FXCB
	INTEGER FXCFL,FXCBL,FXCFR,FXCBR
        INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),comm
        INTEGER myid,ntasks,ierr
        TYPE(NAN_t) :: NANS
        TYPE(EF_t) :: EFS
        TYPE(NTOT_t) :: NTOT
        TYPE(UT_t) :: UT, UTM
        TYPE(NRXF_t) :: NRXF
	DO I=1,6*NN
	A(I)=0.0
	D(I)=0.0
        END DO

 	DO X=1,NTOT%M
	N1=NANS%M(1,X)
	N2=NANS%M(2,X)
	X1=NRXF%M(1,N1)
	Y1=NRXF%M(2,N1)
	Z1=NRXF%M(3,N1)
	X2=NRXF%M(1,N2)
	Y2=NRXF%M(2,N2)
	Z2=NRXF%M(3,N2)


	DX1=UT%M(6*N1-5)
	DY1=UT%M(6*N1-4)
	DZ1=UT%M(6*N1-3)
	DX2=UT%M(6*N2-5)
	DY2=UT%M(6*N2-4)
	DZ2=UT%M(6*N2-3)

        DUT(6*N1-5)=UT%M(6*N1-5)-UTM%M(6*N1-5)
        DUT(6*N1-4)=UT%M(6*N1-4)-UTM%M(6*N1-4)
        DUT(6*N1-3)=UT%M(6*N1-3)-UTM%M(6*N1-3)
        DUT(6*N1-2)=UT%M(6*N1-2)-UTM%M(6*N1-2)
        DUT(6*N1-1)=UT%M(6*N1-1)-UTM%M(6*N1-1)
        DUT(6*N1-0)=UT%M(6*N1-0)-UTM%M(6*N1-0)
                                                        
        DUT(6*N2-5)=UT%M(6*N2-5)-UTM%M(6*N2-5)
        DUT(6*N2-4)=UT%M(6*N2-4)-UTM%M(6*N2-4)
        DUT(6*N2-3)=UT%M(6*N2-3)-UTM%M(6*N2-3)
        DUT(6*N2-2)=UT%M(6*N2-2)-UTM%M(6*N2-2)
        DUT(6*N2-1)=UT%M(6*N2-1)-UTM%M(6*N2-1)
        DUT(6*N2-0)=UT%M(6*N2-0)-UTM%M(6*N2-0)

        IF (EFS % M(X).NE.0.0) THEN
        CALL AMAT(EFS % M(X),S,EFS % M(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
             X2+DX2,Y2+DY2,Z2+DZ2,L,DUT,N1,N2,X,RY,CT)
	ELSE
        CT(12*X-11)=0.0
        CT(12*X-10)=0.0
        CT(12*X-9)=0.0
        CT(12*X-8)=0.0
        CT(12*X-7)=0.0
        CT(12*X-6)=0.0
        CT(12*X-5)=0.0
        CT(12*X-4)=0.0
        CT(12*X-3)=0.0
        CT(12*X-2)=0.0
        CT(12*X-1)=0.0
        CT(12*X-0)=0.0
        ENDIF
      END DO

      dest=myid+1
      source=myid-1
      tag=133
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) &
     CALL MPI_Send(UTM%M,6*NN,MPI_DOUBLE_PRECISION, &
     dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (MOD(myid,ntasks/YN).ne.0) &
     CALL MPI_Recv(UTM%L,6*PNN(source),MPI_DOUBLE_PRECISION, &
     source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-1
      source=myid+1
      tag=135
      IF (MOD(myid,ntasks/YN).ne.0) &
     CALL MPI_Send(UTM%M,6*NN,MPI_DOUBLE_PRECISION, &
     dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) &
     CALL MPI_Recv(UTM%R,6*PNN(source),MPI_DOUBLE_PRECISION, &
     source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid+ntasks/YN
      source=myid-ntasks/YN
      tag=137
      IF (myid.lt.(YN-1)*ntasks/YN) &
     CALL MPI_Send(UTM%M,6*NN,MPI_DOUBLE_PRECISION, &
     dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.ge.ntasks/YN) &
     CALL MPI_Recv(UTM%B,6*PNN(source),MPI_DOUBLE_PRECISION, &
     source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid+ntasks/YN+1
      source=myid-ntasks/YN-1
      tag=139
      IF (myid.lt.(YN-1)*ntasks/YN &
     .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) &
     CALL MPI_Send(UTM%M,6*NN,MPI_DOUBLE_PRECISION, &
     dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) &
     CALL MPI_Recv(UTM%BL,6*PNN(source),MPI_DOUBLE_PRECISION, &
     source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid+ntasks/YN-1
      source=myid-ntasks/YN+1
      tag=141
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) &
           CALL MPI_Send(UTM%M,6*NN,MPI_DOUBLE_PRECISION, &
           dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) &
           CALL MPI_Recv(UTM%BR,6*PNN(source),MPI_DOUBLE_PRECISION, &
           source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-ntasks/YN
      source=myid+ntasks/YN
      tag=143
      IF (myid.ge.ntasks/YN) &
      CALL MPI_Send(UTM%M,6*NN,MPI_DOUBLE_PRECISION, &
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN) &
      CALL MPI_Recv(UTM%F,6*PNN(source),MPI_DOUBLE_PRECISION, &
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-ntasks/YN+1
      source=myid+ntasks/YN-1
      tag=145
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) &
      CALL MPI_Send(UTM%M,6*NN,MPI_DOUBLE_PRECISION, &
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) &
      CALL MPI_Recv(UTM%FL,6*PNN(source),MPI_DOUBLE_PRECISION, &
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-ntasks/YN-1
      source=myid+ntasks/YN+1
      tag=147
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) &
      CALL MPI_Send(UTM%M,6*NN,MPI_DOUBLE_PRECISION, &
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN &
      .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) &
      CALL MPI_Recv(UTM%FR,6*PNN(source),MPI_DOUBLE_PRECISION, &
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

!------------------------------------------------------
        IF (MOD(myid,ntasks/YN).ne.0) THEN
  	DO X=1,NTOT%L
	XL=NTOT%M+X
	N1=NANS % L(1,X)
	N2=NANS % L(2,X)
	X1=NRXF%L(1,N1)
	Y1=NRXF%L(2,N1)
	Z1=NRXF%L(3,N1)
	X2=NRXF%M(1,N2)
	Y2=NRXF%M(2,N2)
	Z2=NRXF%M(3,N2)
	DX1=UT%L(6*N1-5)
	DY1=UT%L(6*N1-4)
	DZ1=UT%L(6*N1-3)
	DX2=UT%M(6*N2-5)
	DY2=UT%M(6*N2-4)
	DZ2=UT%M(6*N2-3)
        DUT(6*N1-5)=UT%L(6*N1-5)-UTM%L(6*N1-5)
        DUT(6*N1-4)=UT%L(6*N1-4)-UTM%L(6*N1-4)
        DUT(6*N1-3)=UT%L(6*N1-3)-UTM%L(6*N1-3)
        DUT(6*N1-2)=UT%L(6*N1-2)-UTM%L(6*N1-2)
        DUT(6*N1-1)=UT%L(6*N1-1)-UTM%L(6*N1-1)
        DUT(6*N1-0)=UT%L(6*N1-0)-UTM%L(6*N1-0)
        DUT(6*N2-5)=UT%M(6*N2-5)-UTM%M(6*N2-5)
        DUT(6*N2-4)=UT%M(6*N2-4)-UTM%M(6*N2-4)
        DUT(6*N2-3)=UT%M(6*N2-3)-UTM%M(6*N2-3)
        DUT(6*N2-2)=UT%M(6*N2-2)-UTM%M(6*N2-2)
        DUT(6*N2-1)=UT%M(6*N2-1)-UTM%M(6*N2-1)
        DUT(6*N2-0)=UT%M(6*N2-0)-UTM%M(6*N2-0)
        IF (EFS % L(X).NE.0.0) THEN
        CALL AMAT(EFS % L(X),S,EFS % L(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
        X2+DX2,Y2+DY2,Z2+DZ2,L,DUT,N1,N2,XL,RY,CT)
	ELSE
        CT(12*XL-11)=0.0
        CT(12*XL-10)=0.0
        CT(12*XL-9)=0.0
        CT(12*XL-8)=0.0
        CT(12*XL-7)=0.0
        CT(12*XL-6)=0.0
        CT(12*XL-5)=0.0
        CT(12*XL-4)=0.0
        CT(12*XL-3)=0.0
        CT(12*XL-2)=0.0
        CT(12*XL-1)=0.0
        CT(12*XL-0)=0.0
        ENDIF
	ENDDO
	ENDIF

        IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
  	DO X=1,NTOT%R
	XR=NTOT%M+NTOT%L+X
	N1=NANS % R(1,X)
	N2=NANS % R(2,X)
	X1=NRXF%R(1,N1)
	Y1=NRXF%R(2,N1)
	Z1=NRXF%R(3,N1)
	X2=NRXF%M(1,N2)
	Y2=NRXF%M(2,N2)
	Z2=NRXF%M(3,N2)
	DX1=UT%R(6*N1-5)
	DY1=UT%R(6*N1-4)
	DZ1=UT%R(6*N1-3)
	DX2=UT%M(6*N2-5)
	DY2=UT%M(6*N2-4)
	DZ2=UT%M(6*N2-3)
        DUT(6*N1-5)=UT%R(6*N1-5)-UTM%R(6*N1-5)
        DUT(6*N1-4)=UT%R(6*N1-4)-UTM%R(6*N1-4)
        DUT(6*N1-3)=UT%R(6*N1-3)-UTM%R(6*N1-3)
        DUT(6*N1-2)=UT%R(6*N1-2)-UTM%R(6*N1-2)
        DUT(6*N1-1)=UT%R(6*N1-1)-UTM%R(6*N1-1)
        DUT(6*N1-0)=UT%R(6*N1-0)-UTM%R(6*N1-0)
        DUT(6*N2-5)=UT%M(6*N2-5)-UTM%M(6*N2-5)
        DUT(6*N2-4)=UT%M(6*N2-4)-UTM%M(6*N2-4)
        DUT(6*N2-3)=UT%M(6*N2-3)-UTM%M(6*N2-3)
        DUT(6*N2-2)=UT%M(6*N2-2)-UTM%M(6*N2-2)
        DUT(6*N2-1)=UT%M(6*N2-1)-UTM%M(6*N2-1)
        DUT(6*N2-0)=UT%M(6*N2-0)-UTM%M(6*N2-0)
        IF (EFS % R(X).NE.0.0) THEN
        CALL AMAT(EFS % R(X),S,EFS % R(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
        X2+DX2,Y2+DY2,Z2+DZ2,L,DUT,N1,N2,XR,RY,CT)
	ELSE
        CT(12*XR-11)=0.0
        CT(12*XR-10)=0.0
        CT(12*XR-9)=0.0
        CT(12*XR-8)=0.0
        CT(12*XR-7)=0.0
        CT(12*XR-6)=0.0
        CT(12*XR-5)=0.0
        CT(12*XR-4)=0.0
        CT(12*XR-3)=0.0
        CT(12*XR-2)=0.0
        CT(12*XR-1)=0.0
        CT(12*XR-0)=0.0
        ENDIF
	ENDDO
	ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN) THEN
        DO X=1,NTOT%F
        XR=NTOT%M+NTOT%L+NTOT%R+X
        N1=NANS % F(1,X)
        N2=NANS % F(2,X)
        X1=NRXF%F(1,N1)
        Y1=NRXF%F(2,N1)
        Z1=NRXF%F(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
        DX1=UT%F(6*N1-5)
        DY1=UT%F(6*N1-4)
        DZ1=UT%F(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
        DUT(6*N1-5)=UT%F(6*N1-5)-UTM%F(6*N1-5)
        DUT(6*N1-4)=UT%F(6*N1-4)-UTM%F(6*N1-4)
        DUT(6*N1-3)=UT%F(6*N1-3)-UTM%F(6*N1-3)
        DUT(6*N1-2)=UT%F(6*N1-2)-UTM%F(6*N1-2)
        DUT(6*N1-1)=UT%F(6*N1-1)-UTM%F(6*N1-1)
        DUT(6*N1-0)=UT%F(6*N1-0)-UTM%F(6*N1-0)
        DUT(6*N2-5)=UT%M(6*N2-5)-UTM%M(6*N2-5)
        DUT(6*N2-4)=UT%M(6*N2-4)-UTM%M(6*N2-4)
        DUT(6*N2-3)=UT%M(6*N2-3)-UTM%M(6*N2-3)
        DUT(6*N2-2)=UT%M(6*N2-2)-UTM%M(6*N2-2)
        DUT(6*N2-1)=UT%M(6*N2-1)-UTM%M(6*N2-1)
        DUT(6*N2-0)=UT%M(6*N2-0)-UTM%M(6*N2-0)
        IF (EFS % F(X).NE.0.0) THEN
        CALL AMAT(EFS % F(X),S,EFS % F(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
        X2+DX2,Y2+DY2,Z2+DZ2,L,DUT,N1,N2,XR,RY,CT)
        ELSE
        CT(12*XR-11)=0.0
        CT(12*XR-10)=0.0
        CT(12*XR-9)=0.0
        CT(12*XR-8)=0.0
        CT(12*XR-7)=0.0
        CT(12*XR-6)=0.0
        CT(12*XR-5)=0.0
        CT(12*XR-4)=0.0
        CT(12*XR-3)=0.0
        CT(12*XR-2)=0.0
        CT(12*XR-1)=0.0
        CT(12*XR-0)=0.0
        ENDIF
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,NTOT%FL
        XR=NTOT%M+NTOT%L+NTOT%R+NTOT%F+X
        N1=NANS % FL(1,X)
        N2=NANS % FL(2,X)
        X1=NRXF%FL(1,N1)
        Y1=NRXF%FL(2,N1)
        Z1=NRXF%FL(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
        DX1=UT%FL(6*N1-5)
        DY1=UT%FL(6*N1-4)
        DZ1=UT%FL(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
        DUT(6*N1-5)=UT%FL(6*N1-5)-UTM%FL(6*N1-5)
        DUT(6*N1-4)=UT%FL(6*N1-4)-UTM%FL(6*N1-4)
        DUT(6*N1-3)=UT%FL(6*N1-3)-UTM%FL(6*N1-3)
        DUT(6*N1-2)=UT%FL(6*N1-2)-UTM%FL(6*N1-2)
        DUT(6*N1-1)=UT%FL(6*N1-1)-UTM%FL(6*N1-1)
        DUT(6*N1-0)=UT%FL(6*N1-0)-UTM%FL(6*N1-0)
        DUT(6*N2-5)=UT%M(6*N2-5)-UTM%M(6*N2-5)
        DUT(6*N2-4)=UT%M(6*N2-4)-UTM%M(6*N2-4)
        DUT(6*N2-3)=UT%M(6*N2-3)-UTM%M(6*N2-3)
        DUT(6*N2-2)=UT%M(6*N2-2)-UTM%M(6*N2-2)
        DUT(6*N2-1)=UT%M(6*N2-1)-UTM%M(6*N2-1)
        DUT(6*N2-0)=UT%M(6*N2-0)-UTM%M(6*N2-0)
        IF (EFS % FL(X).NE.0.0) THEN
        CALL AMAT(EFS % FL(X),S,EFS % FL(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
        X2+DX2,Y2+DY2,Z2+DZ2,L,DUT,N1,N2,XR,RY,CT)
        ELSE
        CT(12*XR-11)=0.0
        CT(12*XR-10)=0.0
        CT(12*XR-9)=0.0
        CT(12*XR-8)=0.0
        CT(12*XR-7)=0.0
        CT(12*XR-6)=0.0
        CT(12*XR-5)=0.0
        CT(12*XR-4)=0.0
        CT(12*XR-3)=0.0
        CT(12*XR-2)=0.0
        CT(12*XR-1)=0.0
        CT(12*XR-0)=0.0
        ENDIF
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN &
      	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,NTOT%FR
        XR=NTOT%M+NTOT%L+NTOT%R+NTOT%F+NTOT%FL+X
        N1=NANS % FR(1,X)
        N2=NANS % FR(2,X)
        X1=NRXF%FR(1,N1)
        Y1=NRXF%FR(2,N1)
        Z1=NRXF%FR(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
        DX1=UT%FR(6*N1-5)
        DY1=UT%FR(6*N1-4)
        DZ1=UT%FR(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
        DUT(6*N1-5)=UT%FR(6*N1-5)-UTM%FR(6*N1-5)
        DUT(6*N1-4)=UT%FR(6*N1-4)-UTM%FR(6*N1-4)
        DUT(6*N1-3)=UT%FR(6*N1-3)-UTM%FR(6*N1-3)
        DUT(6*N1-2)=UT%FR(6*N1-2)-UTM%FR(6*N1-2)
        DUT(6*N1-1)=UT%FR(6*N1-1)-UTM%FR(6*N1-1)
        DUT(6*N1-0)=UT%FR(6*N1-0)-UTM%FR(6*N1-0)
        DUT(6*N2-5)=UT%M(6*N2-5)-UTM%M(6*N2-5)
        DUT(6*N2-4)=UT%M(6*N2-4)-UTM%M(6*N2-4)
        DUT(6*N2-3)=UT%M(6*N2-3)-UTM%M(6*N2-3)
        DUT(6*N2-2)=UT%M(6*N2-2)-UTM%M(6*N2-2)
        DUT(6*N2-1)=UT%M(6*N2-1)-UTM%M(6*N2-1)
        DUT(6*N2-0)=UT%M(6*N2-0)-UTM%M(6*N2-0)
        IF (EFS % FR(X).NE.0.0) THEN
        CALL AMAT(EFS % FR(X),S,EFS % FR(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
        X2+DX2,Y2+DY2,Z2+DZ2,L,DUT,N1,N2,XR,RY,CT)
        ELSE
        CT(12*XR-11)=0.0
        CT(12*XR-10)=0.0
        CT(12*XR-9)=0.0
        CT(12*XR-8)=0.0
        CT(12*XR-7)=0.0
        CT(12*XR-6)=0.0
        CT(12*XR-5)=0.0
        CT(12*XR-4)=0.0
        CT(12*XR-3)=0.0
        CT(12*XR-2)=0.0
        CT(12*XR-1)=0.0
        CT(12*XR-0)=0.0
        ENDIF
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN) THEN
        DO X=1,NTOT%B
        XR=NTOT%M+NTOT%L+NTOT%R+NTOT%F+NTOT%FL+NTOT%FR+X
        N1=NANS % B(1,X)
        N2=NANS % B(2,X)
        X1=NRXF%B(1,N1)
        Y1=NRXF%B(2,N1)
        Z1=NRXF%B(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
	DX1=UT%B(6*N1-5)
        DY1=UT%B(6*N1-4)
        DZ1=UT%B(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
	DUT(6*N1-5)=UT%B(6*N1-5)-UTM%B(6*N1-5)
        DUT(6*N1-4)=UT%B(6*N1-4)-UTM%B(6*N1-4)
        DUT(6*N1-3)=UT%B(6*N1-3)-UTM%B(6*N1-3)
        DUT(6*N1-2)=UT%B(6*N1-2)-UTM%B(6*N1-2)
        DUT(6*N1-1)=UT%B(6*N1-1)-UTM%B(6*N1-1)
        DUT(6*N1-0)=UT%B(6*N1-0)-UTM%B(6*N1-0)
	DUT(6*N2-5)=UT%M(6*N2-5)-UTM%M(6*N2-5)
        DUT(6*N2-4)=UT%M(6*N2-4)-UTM%M(6*N2-4)
        DUT(6*N2-3)=UT%M(6*N2-3)-UTM%M(6*N2-3)
        DUT(6*N2-2)=UT%M(6*N2-2)-UTM%M(6*N2-2)
        DUT(6*N2-1)=UT%M(6*N2-1)-UTM%M(6*N2-1)
        DUT(6*N2-0)=UT%M(6*N2-0)-UTM%M(6*N2-0)
	IF (EFS % B(X).NE.0.0) THEN
        CALL AMAT(EFS % B(X),S,EFS % B(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
        X2+DX2,Y2+DY2,Z2+DZ2,L,DUT,N1,N2,XR,RY,CT)
        ELSE
        CT(12*XR-11)=0.0
        CT(12*XR-10)=0.0
        CT(12*XR-9)=0.0
        CT(12*XR-8)=0.0
        CT(12*XR-7)=0.0
        CT(12*XR-6)=0.0
        CT(12*XR-5)=0.0
        CT(12*XR-4)=0.0
        CT(12*XR-3)=0.0
        CT(12*XR-2)=0.0
        CT(12*XR-1)=0.0
        CT(12*XR-0)=0.0
        ENDIF
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,NTOT%BL
        XR=NTOT%M+NTOT%L+NTOT%R+NTOT%F+NTOT%FL+NTOT%FR+NTOT%B+X
        N1=NANS % BL(1,X)
        N2=NANS % BL(2,X)
        X1=NRXF%BL(1,N1)
        Y1=NRXF%BL(2,N1)
        Z1=NRXF%BL(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
	DX1=UT%BL(6*N1-5)
        DY1=UT%BL(6*N1-4)
        DZ1=UT%BL(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
	DUT(6*N1-5)=UT%BL(6*N1-5)-UTM%BL(6*N1-5)
        DUT(6*N1-4)=UT%BL(6*N1-4)-UTM%BL(6*N1-4)
        DUT(6*N1-3)=UT%BL(6*N1-3)-UTM%BL(6*N1-3)
        DUT(6*N1-2)=UT%BL(6*N1-2)-UTM%BL(6*N1-2)
        DUT(6*N1-1)=UT%BL(6*N1-1)-UTM%BL(6*N1-1)
        DUT(6*N1-0)=UT%BL(6*N1-0)-UTM%BL(6*N1-0)
	DUT(6*N2-5)=UT%M(6*N2-5)-UTM%M(6*N2-5)
        DUT(6*N2-4)=UT%M(6*N2-4)-UTM%M(6*N2-4)
        DUT(6*N2-3)=UT%M(6*N2-3)-UTM%M(6*N2-3)
        DUT(6*N2-2)=UT%M(6*N2-2)-UTM%M(6*N2-2)
        DUT(6*N2-1)=UT%M(6*N2-1)-UTM%M(6*N2-1)
        DUT(6*N2-0)=UT%M(6*N2-0)-UTM%M(6*N2-0)
	IF (EFS % BL(X).NE.0.0) THEN
        CALL AMAT(EFS % BL(X),S,EFS % BL(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
        X2+DX2,Y2+DY2,Z2+DZ2,L,DUT,N1,N2,XR,RY,CT)
        ELSE
        CT(12*XR-11)=0.0
        CT(12*XR-10)=0.0
        CT(12*XR-9)=0.0
        CT(12*XR-8)=0.0
        CT(12*XR-7)=0.0
        CT(12*XR-6)=0.0
        CT(12*XR-5)=0.0
        CT(12*XR-4)=0.0
        CT(12*XR-3)=0.0
        CT(12*XR-2)=0.0
        CT(12*XR-1)=0.0
        CT(12*XR-0)=0.0
        ENDIF
	ENDDO
        ENDIF

	IF (myid.gt.(ntasks-1)/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,NTOT%BR
        XR=NTOT%M+NTOT%L+NTOT%R+NTOT%F+NTOT%FL+NTOT%FR+NTOT%B+NTOT%BL+X
        N1=NANS % BR(1,X)
        N2=NANS % BR(2,X)
        X1=NRXF%BR(1,N1)
        Y1=NRXF%BR(2,N1)
        Z1=NRXF%BR(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
	DX1=UT%BR(6*N1-5)
        DY1=UT%BR(6*N1-4)
        DZ1=UT%BR(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
	DUT(6*N1-5)=UT%BR(6*N1-5)-UTM%BR(6*N1-5)
        DUT(6*N1-4)=UT%BR(6*N1-4)-UTM%BR(6*N1-4)
        DUT(6*N1-3)=UT%BR(6*N1-3)-UTM%BR(6*N1-3)
        DUT(6*N1-2)=UT%BR(6*N1-2)-UTM%BR(6*N1-2)
        DUT(6*N1-1)=UT%BR(6*N1-1)-UTM%BR(6*N1-1)
        DUT(6*N1-0)=UT%BR(6*N1-0)-UTM%BR(6*N1-0)
	DUT(6*N2-5)=UT%M(6*N2-5)-UTM%M(6*N2-5)
        DUT(6*N2-4)=UT%M(6*N2-4)-UTM%M(6*N2-4)
        DUT(6*N2-3)=UT%M(6*N2-3)-UTM%M(6*N2-3)
        DUT(6*N2-2)=UT%M(6*N2-2)-UTM%M(6*N2-2)
        DUT(6*N2-1)=UT%M(6*N2-1)-UTM%M(6*N2-1)
        DUT(6*N2-0)=UT%M(6*N2-0)-UTM%M(6*N2-0)
	IF (EFS % BR(X).NE.0.0) THEN
        CALL AMAT(EFS % BR(X),S,EFS % BR(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
        X2+DX2,Y2+DY2,Z2+DZ2,L,DUT,N1,N2,XR,RY,CT)
        ELSE
        CT(12*XR-11)=0.0
        CT(12*XR-10)=0.0
        CT(12*XR-9)=0.0
        CT(12*XR-8)=0.0
        CT(12*XR-7)=0.0
        CT(12*XR-6)=0.0
        CT(12*XR-5)=0.0
        CT(12*XR-4)=0.0
        CT(12*XR-3)=0.0
        CT(12*XR-2)=0.0
        CT(12*XR-1)=0.0
        CT(12*XR-0)=0.0
        ENDIF
	ENDDO
        ENDIF
!----------------------------------------------

        DO X=1,NTOT%M
        N1=NANS % M(1,X)
        N2=NANS % M(2,X)
        X1=NRXF%M(1,N1)
        Y1=NRXF%M(2,N1)
        Z1=NRXF%M(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
        DX1=UT%M(6*N1-5)
        DY1=UT%M(6*N1-4)
        DZ1=UT%M(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
        CALL TTMAT(X1+DX1,Y1+DY1,Z1+DZ1,X2+DX2,Y2+DY2,Z2+DZ2,RY,TT)
        A(6*N1-5)= A(6*N1-5) &
                  +  TT(1,1)*CT(12*X-11) + TT(1,2)*CT(12*X-10) &
                  + TT(1,3)*CT(12*X-9)
        A(6*N1-4)= A(6*N1-4) &
                  +  TT(2,1)*CT(12*X-11) + TT(2,2)*CT(12*X-10) &
                  + TT(2,3)*CT(12*X-9)
        A(6*N1-3)= A(6*N1-3) &
                  +  TT(3,1)*CT(12*X-11) + TT(3,2)*CT(12*X-10) &
                  + TT(3,3)*CT(12*X-9)
        A(6*N1-2)= A(6*N1-2) &
                  + TT(4,4)*CT(12*X-8) &
                  + TT(4,5)*CT(12*X-7) + TT(4,6)*CT(12*X-6)
        A(6*N1-1)= A(6*N1-1) &
                  + TT(5,4)*CT(12*X-8) &
                  + TT(5,5)*CT(12*X-7) + TT(5,6)*CT(12*X-6)
        A(6*N1-0)= A(6*N1-0) &
       	          + TT(6,4)*CT(12*X-8) &
                  + TT(6,5)*CT(12*X-7) + TT(6,6)*CT(12*X-6)
        A(6*N2-5)= A(6*N2-5) &
                  + TT(7,7)*CT(12*X-5) + TT(7,8)*CT(12*X-4) &
                  + TT(7,9)*CT(12*X-3)
        A(6*N2-4)= A(6*N2-4) &
                  + TT(8,7)*CT(12*X-5) + TT(8,8)*CT(12*X-4) &
                  + TT(8,9)*CT(12*X-3)
        A(6*N2-3)= A(6*N2-3) &
                  + TT(9,7)*CT(12*X-5) + TT(9,8)*CT(12*X-4) &
                  + TT(9,9)*CT(12*X-3)
        A(6*N2-2)= A(6*N2-2) &
                  + TT(10,10)*CT(12*X-2) + TT(10,11)*CT(12*X-1) &
      	          + TT(10,12)*CT(12*X-0)
        A(6*N2-1)= A(6*N2-1) &
      	          + TT(11,10)*CT(12*X-2) + TT(11,11)*CT(12*X-1) &
      	          + TT(11,12)*CT(12*X-0)
        A(6*N2-0)= A(6*N2-0) &
      	          + TT(12,10)*CT(12*X-2) + TT(12,11)*CT(12*X-1) &
      	          + TT(12,12)*CT(12*X-0)
	ENDDO

        IF (MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,NTOT%L
	XL=NTOT%M+X
        N1=NANS % L(1,X)
        N2=NANS % L(2,X)
        X1=NRXF%L(1,N1)
        Y1=NRXF%L(2,N1)
        Z1=NRXF%L(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
        DX1=UT%L(6*N1-5)
        DY1=UT%L(6*N1-4)
        DZ1=UT%L(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
        CALL TTMAT(X1+DX1,Y1+DY1,Z1+DZ1,X2+DX2,Y2+DY2,Z2+DZ2,RY,TT)
        A(6*N2-5)= A(6*N2-5) &
                  + TT(7,7)*CT(12*XL-5) + TT(7,8)*CT(12*XL-4) &
                  + TT(7,9)*CT(12*XL-3)
        A(6*N2-4)= A(6*N2-4) &
                  + TT(8,7)*CT(12*XL-5) + TT(8,8)*CT(12*XL-4) &
                  + TT(8,9)*CT(12*XL-3)
        A(6*N2-3)= A(6*N2-3) &
                  + TT(9,7)*CT(12*XL-5) + TT(9,8)*CT(12*XL-4) &
                  + TT(9,9)*CT(12*XL-3)
        A(6*N2-2)= A(6*N2-2) &
      	          + TT(10,10)*CT(12*XL-2) &
                  + TT(10,11)*CT(12*XL-1) + TT(10,12)*CT(12*XL-0)
        A(6*N2-1)= A(6*N2-1) &
      	          + TT(11,10)*CT(12*XL-2) &
                  + TT(11,11)*CT(12*XL-1) + TT(11,12)*CT(12*XL-0)
        A(6*N2-0)= A(6*N2-0) &
      	          + TT(12,10)*CT(12*XL-2) &
                  + TT(12,11)*CT(12*XL-1) + TT(12,12)*CT(12*XL-0)
	ENDDO
	ENDIF

        IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,NTOT%R
 	XR=NTOT%M+NTOT%L+X
        N1=NANS % R(1,X)
        N2=NANS % R(2,X)
        X1=NRXF%R(1,N1)
        Y1=NRXF%R(2,N1)
        Z1=NRXF%R(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
        DX1=UT%R(6*N1-5)
        DY1=UT%R(6*N1-4)
        DZ1=UT%R(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
        CALL TTMAT(X1+DX1,Y1+DY1,Z1+DZ1,X2+DX2,Y2+DY2,Z2+DZ2,RY,TT)
        A(6*N2-5)= A(6*N2-5) &
                  + TT(7,7)*CT(12*XR-5) + TT(7,8)*CT(12*XR-4) &
                  + TT(7,9)*CT(12*XR-3)
        A(6*N2-4)= A(6*N2-4) &
                  + TT(8,7)*CT(12*XR-5) + TT(8,8)*CT(12*XR-4) &
                  + TT(8,9)*CT(12*XR-3)
        A(6*N2-3)= A(6*N2-3) &
                  + TT(9,7)*CT(12*XR-5) + TT(9,8)*CT(12*XR-4) &
                  + TT(9,9)*CT(12*XR-3)
        A(6*N2-2)= A(6*N2-2) &
      	          + TT(10,10)*CT(12*XR-2) &
                  + TT(10,11)*CT(12*XR-1) + TT(10,12)*CT(12*XR-0)
        A(6*N2-1)= A(6*N2-1) &
      	          + TT(11,10)*CT(12*XR-2) &
                  + TT(11,11)*CT(12*XR-1) + TT(11,12)*CT(12*XR-0)
        A(6*N2-0)= A(6*N2-0) &
      	          + TT(12,10)*CT(12*XR-2) &
                  + TT(12,11)*CT(12*XR-1) + TT(12,12)*CT(12*XR-0)
	ENDDO
	ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN) THEN
        DO X=1,NTOT%F
        XR=NTOT%M+NTOT%L+NTOT%R+X
        N1=NANS % F(1,X)
        N2=NANS % F(2,X)
        X1=NRXF%F(1,N1)
        Y1=NRXF%F(2,N1)
        Z1=NRXF%F(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
        DX1=UT%F(6*N1-5)
        DY1=UT%F(6*N1-4)
        DZ1=UT%F(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
        CALL TTMAT(X1+DX1,Y1+DY1,Z1+DZ1,X2+DX2,Y2+DY2,Z2+DZ2,RY,TT)
        A(6*N2-5)= A(6*N2-5) &
                  + TT(7,7)*CT(12*XR-5) + TT(7,8)*CT(12*XR-4) &
                  + TT(7,9)*CT(12*XR-3)
        A(6*N2-4)= A(6*N2-4) &
                  + TT(8,7)*CT(12*XR-5) + TT(8,8)*CT(12*XR-4) &
                  + TT(8,9)*CT(12*XR-3)
        A(6*N2-3)= A(6*N2-3) &
                  + TT(9,7)*CT(12*XR-5) + TT(9,8)*CT(12*XR-4) &
                  + TT(9,9)*CT(12*XR-3)
	A(6*N2-2)= A(6*N2-2) &
                  + TT(10,10)*CT(12*XR-2) &
                  + TT(10,11)*CT(12*XR-1) + TT(10,12)*CT(12*XR-0)
	A(6*N2-1)= A(6*N2-1) &
                  + TT(11,10)*CT(12*XR-2) &
                  + TT(11,11)*CT(12*XR-1) + TT(11,12)*CT(12*XR-0)
        A(6*N2-0)= A(6*N2-0) &
                  + TT(12,10)*CT(12*XR-2) &
                  + TT(12,11)*CT(12*XR-1) + TT(12,12)*CT(12*XR-0)
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,NTOT%FL
        XR=NTOT%M+NTOT%L+NTOT%R+NTOT%F+X
        N1=NANS % FL(1,X)
        N2=NANS % FL(2,X)
        X1=NRXF%FL(1,N1)
        Y1=NRXF%FL(2,N1)
        Z1=NRXF%FL(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
        DX1=UT%FL(6*N1-5)
        DY1=UT%FL(6*N1-4)
        DZ1=UT%FL(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
        CALL TTMAT(X1+DX1,Y1+DY1,Z1+DZ1,X2+DX2,Y2+DY2,Z2+DZ2,RY,TT)
        A(6*N2-5)= A(6*N2-5) &
                  + TT(7,7)*CT(12*XR-5) + TT(7,8)*CT(12*XR-4) &
                  + TT(7,9)*CT(12*XR-3)
        A(6*N2-4)= A(6*N2-4) &
                  + TT(8,7)*CT(12*XR-5) + TT(8,8)*CT(12*XR-4) &
                  + TT(8,9)*CT(12*XR-3)
        A(6*N2-3)= A(6*N2-3) &
                  + TT(9,7)*CT(12*XR-5) + TT(9,8)*CT(12*XR-4) &
                  + TT(9,9)*CT(12*XR-3)
	A(6*N2-2)= A(6*N2-2) &
                  + TT(10,10)*CT(12*XR-2) &
                  + TT(10,11)*CT(12*XR-1) + TT(10,12)*CT(12*XR-0)
	A(6*N2-1)= A(6*N2-1) &
                  + TT(11,10)*CT(12*XR-2) &
                  + TT(11,11)*CT(12*XR-1) + TT(11,12)*CT(12*XR-0)
        A(6*N2-0)= A(6*N2-0) &
                  + TT(12,10)*CT(12*XR-2) &
                  + TT(12,11)*CT(12*XR-1) + TT(12,12)*CT(12*XR-0)
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN &
      	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,NTOT%FR
        XR=NTOT%M+NTOT%L+NTOT%R+NTOT%F+NTOT%FL+X
        N1=NANS % FR(1,X)
        N2=NANS % FR(2,X)
        X1=NRXF%FR(1,N1)
        Y1=NRXF%FR(2,N1)
        Z1=NRXF%FR(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
        DX1=UT%FR(6*N1-5)
        DY1=UT%FR(6*N1-4)
        DZ1=UT%FR(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
        CALL TTMAT(X1+DX1,Y1+DY1,Z1+DZ1,X2+DX2,Y2+DY2,Z2+DZ2,RY,TT)
        A(6*N2-5)= A(6*N2-5) &
                  + TT(7,7)*CT(12*XR-5) + TT(7,8)*CT(12*XR-4) &
                  + TT(7,9)*CT(12*XR-3)
        A(6*N2-4)= A(6*N2-4) &
                  + TT(8,7)*CT(12*XR-5) + TT(8,8)*CT(12*XR-4) &
                  + TT(8,9)*CT(12*XR-3)
        A(6*N2-3)= A(6*N2-3) &
                  + TT(9,7)*CT(12*XR-5) + TT(9,8)*CT(12*XR-4) &
                  + TT(9,9)*CT(12*XR-3)
	A(6*N2-2)= A(6*N2-2) &
                  + TT(10,10)*CT(12*XR-2) &
                  + TT(10,11)*CT(12*XR-1) + TT(10,12)*CT(12*XR-0)
	A(6*N2-1)= A(6*N2-1) &
                  + TT(11,10)*CT(12*XR-2) &
                  + TT(11,11)*CT(12*XR-1) + TT(11,12)*CT(12*XR-0)
        A(6*N2-0)= A(6*N2-0) &
                  + TT(12,10)*CT(12*XR-2) &
                  + TT(12,11)*CT(12*XR-1) + TT(12,12)*CT(12*XR-0)
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN) THEN
        DO X=1,NTOT%B
        XR=NTOT%M+NTOT%L+NTOT%R+NTOT%F+NTOT%FL+NTOT%FR+X
	N1=NANS % B(1,X)
        N2=NANS % B(2,X)
	X1=NRXF%B(1,N1)
        Y1=NRXF%B(2,N1)
        Z1=NRXF%B(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
	DX1=UT%B(6*N1-5)
        DY1=UT%B(6*N1-4)
        DZ1=UT%B(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
	CALL TTMAT(X1+DX1,Y1+DY1,Z1+DZ1,X2+DX2,Y2+DY2,Z2+DZ2,RY,TT)
	A(6*N2-5)= A(6*N2-5) &
                  + TT(7,7)*CT(12*XR-5) + TT(7,8)*CT(12*XR-4) &
                  + TT(7,9)*CT(12*XR-3)
	A(6*N2-4)= A(6*N2-4) &
                  + TT(8,7)*CT(12*XR-5) + TT(8,8)*CT(12*XR-4) &
                  + TT(8,9)*CT(12*XR-3)
        A(6*N2-3)= A(6*N2-3) &
                  + TT(9,7)*CT(12*XR-5) + TT(9,8)*CT(12*XR-4) &
                  + TT(9,9)*CT(12*XR-3)
	A(6*N2-2)= A(6*N2-2) &
                  + TT(10,10)*CT(12*XR-2) &
                  + TT(10,11)*CT(12*XR-1) + TT(10,12)*CT(12*XR-0)
	A(6*N2-1)= A(6*N2-1) &
                  + TT(11,10)*CT(12*XR-2) &
                  + TT(11,11)*CT(12*XR-1) + TT(11,12)*CT(12*XR-0)
        A(6*N2-0)= A(6*N2-0) &
                  + TT(12,10)*CT(12*XR-2) &
                  + TT(12,11)*CT(12*XR-1) + TT(12,12)*CT(12*XR-0)
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,NTOT%BL
        XR=NTOT%M+NTOT%L+NTOT%R+NTOT%F+NTOT%FL+NTOT%FR+NTOT%B+X
	N1=NANS % BL(1,X)
        N2=NANS % BL(2,X)
	X1=NRXF%BL(1,N1)
        Y1=NRXF%BL(2,N1)
        Z1=NRXF%BL(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
	DX1=UT%BL(6*N1-5)
        DY1=UT%BL(6*N1-4)
        DZ1=UT%BL(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
	CALL TTMAT(X1+DX1,Y1+DY1,Z1+DZ1,X2+DX2,Y2+DY2,Z2+DZ2,RY,TT)
	A(6*N2-5)= A(6*N2-5) &
                  + TT(7,7)*CT(12*XR-5) + TT(7,8)*CT(12*XR-4) &
                  + TT(7,9)*CT(12*XR-3)
	A(6*N2-4)= A(6*N2-4) &
                  + TT(8,7)*CT(12*XR-5) + TT(8,8)*CT(12*XR-4) &
                  + TT(8,9)*CT(12*XR-3)
        A(6*N2-3)= A(6*N2-3) &
                  + TT(9,7)*CT(12*XR-5) + TT(9,8)*CT(12*XR-4) &
                  + TT(9,9)*CT(12*XR-3)
	A(6*N2-2)= A(6*N2-2) &
                  + TT(10,10)*CT(12*XR-2) &
                  + TT(10,11)*CT(12*XR-1) + TT(10,12)*CT(12*XR-0)
	A(6*N2-1)= A(6*N2-1) &
                  + TT(11,10)*CT(12*XR-2) &
                  + TT(11,11)*CT(12*XR-1) + TT(11,12)*CT(12*XR-0)
        A(6*N2-0)= A(6*N2-0) &
                  + TT(12,10)*CT(12*XR-2) &
                  + TT(12,11)*CT(12*XR-1) + TT(12,12)*CT(12*XR-0)
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,NTOT%BR
        XR=NTOT%M+NTOT%L+NTOT%R+NTOT%F+NTOT%FL+NTOT%FR+NTOT%B+NTOT%BL+X
	N1=NANS % BR(1,X)
        N2=NANS % BR(2,X)
	X1=NRXF%BR(1,N1)
        Y1=NRXF%BR(2,N1)
        Z1=NRXF%BR(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
	DX1=UT%BR(6*N1-5)
        DY1=UT%BR(6*N1-4)
        DZ1=UT%BR(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
	CALL TTMAT(X1+DX1,Y1+DY1,Z1+DZ1,X2+DX2,Y2+DY2,Z2+DZ2,RY,TT)
	A(6*N2-5)= A(6*N2-5) &
                  + TT(7,7)*CT(12*XR-5) + TT(7,8)*CT(12*XR-4) &
                  + TT(7,9)*CT(12*XR-3)
	A(6*N2-4)= A(6*N2-4) &
                  + TT(8,7)*CT(12*XR-5) + TT(8,8)*CT(12*XR-4) &
                  + TT(8,9)*CT(12*XR-3)
        A(6*N2-3)= A(6*N2-3) &
                  + TT(9,7)*CT(12*XR-5) + TT(9,8)*CT(12*XR-4) &
                  + TT(9,9)*CT(12*XR-3)
	A(6*N2-2)= A(6*N2-2) &
                  + TT(10,10)*CT(12*XR-2) &
                  + TT(10,11)*CT(12*XR-1) + TT(10,12)*CT(12*XR-0)
	A(6*N2-1)= A(6*N2-1) &
                  + TT(11,10)*CT(12*XR-2) &
                  + TT(11,11)*CT(12*XR-1) + TT(11,12)*CT(12*XR-0)
        A(6*N2-0)= A(6*N2-0) &
                  + TT(12,10)*CT(12*XR-2) &
                  + TT(12,11)*CT(12*XR-1) + TT(12,12)*CT(12*XR-0)
	ENDDO
        ENDIF

!------------------------------------------------------------------------

	DO X=1,NN
	C(6*X-5)= (MFIL(X)/DT**2)*UTM%M(6*X-5)-(2*MFIL(X)/DT**2)*UT%M(6*X-5)
	C(6*X-4)= (MFIL(X)/DT**2)*UTM%M(6*X-4)-(2*MFIL(X)/DT**2)*UT%M(6*X-4)
	C(6*X-3)= (MFIL(X)/DT**2)*UTM%M(6*X-3)-(2*MFIL(X)/DT**2)*UT%M(6*X-3)
	C(6*X-2)= ((MFIL(X)*JS/M)/DT**2)* &
      	UTM%M(6*X-2)-(2*(MFIL(X)*JS/M)/DT**2)*UT%M(6*X-2)
	C(6*X-1)= ((MFIL(X)*JS/M)/DT**2)* &
      	UTM%M(6*X-1)-(2*(MFIL(X)*JS/M)/DT**2)*UT%M(6*X-1)
	C(6*X-0)= ((MFIL(X)*JS/M)/DT**2)* &
      	UTM%M(6*X-0)-(2*(MFIL(X)*JS/M)/DT**2)*UT%M(6*X-0)
	ENDDO

!------------------------------------------------------------------------


	DO X=1,NTOT%M
	IF (EFS % M(X).NE.0.0) THEN
	N1=NANS % M(1,X)
	N2=NANS % M(2,X)
	D(6*N1-5)=D(6*N1-5)+(DMP/DT)*((UT%M(6*N1-5)-UT%M(6*N2-5)) &
      	-(UTM%M(6*N1-5)-UTM%M(6*N2-5)))
	D(6*N1-4)=D(6*N1-4)+(DMP/DT)*((UT%M(6*N1-4)-UT%M(6*N2-4)) &
      	-(UTM%M(6*N1-4)-UTM%M(6*N2-4)))
	D(6*N1-3)=D(6*N1-3)+(DMP/DT)*((UT%M(6*N1-3)-UT%M(6*N2-3)) &
      	-(UTM%M(6*N1-3)-UTM%M(6*N2-3)))
	D(6*N1-2)=D(6*N1-2)+(DMP2/DT)*((UT%M(6*N1-2)-UT%M(6*N2-2)) &
      	-(UTM%M(6*N1-2)-UTM%M(6*N2-2)))
	D(6*N1-1)=D(6*N1-1)+(DMP2/DT)*((UT%M(6*N1-1)-UT%M(6*N2-1)) &
      	-(UTM%M(6*N1-1)-UTM%M(6*N2-1)))
	D(6*N1-0)=D(6*N1-0)+(DMP2/DT)*((UT%M(6*N1-0)-UT%M(6*N2-0)) &
      	-(UTM%M(6*N1-0)-UTM%M(6*N2-0)))
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%M(6*N1-5)) &
      	-(UTM%M(6*N2-5)-UTM%M(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%M(6*N1-4)) &
      	-(UTM%M(6*N2-4)-UTM%M(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%M(6*N1-3)) &
      	-(UTM%M(6*N2-3)-UTM%M(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%M(6*N1-2)) &
      	-(UTM%M(6*N2-2)-UTM%M(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%M(6*N1-1)) &
      	-(UTM%M(6*N2-1)-UTM%M(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%M(6*N1-0)) &
      	-(UTM%M(6*N2-0)-UTM%M(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%M(6*N1-5)) &
     	-(UTM%M(6*N2-5)-UTM%M(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%M(6*N1-4)) &
     	-(UTM%M(6*N2-4)-UTM%M(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%M(6*N1-3)) &
     	-(UTM%M(6*N2-3)-UTM%M(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%M(6*N1-2)) &
     	-(UTM%M(6*N2-2)-UTM%M(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%M(6*N1-1)) &
     	-(UTM%M(6*N2-1)-UTM%M(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%M(6*N1-0)) &
     	-(UTM%M(6*N2-0)-UTM%M(6*N1-0)))**2
	ENDIF
	ENDDO

        IF (MOD(myid,ntasks/YN).ne.0) THEN
	DO X=1,NTOT%L
	IF (EFS % L(X).NE.0.0) THEN
	N1=NANS % L(1,X)
	N2=NANS % L(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%L(6*N1-5)) &
     	-(UTM%M(6*N2-5)-UTM%L(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%L(6*N1-4)) &
     	-(UTM%M(6*N2-4)-UTM%L(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%L(6*N1-3)) &
     	-(UTM%M(6*N2-3)-UTM%L(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%L(6*N1-2)) &
     	-(UTM%M(6*N2-2)-UTM%L(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%L(6*N1-1)) &
     	-(UTM%M(6*N2-1)-UTM%L(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%L(6*N1-0)) &
     	-(UTM%M(6*N2-0)-UTM%L(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%L(6*N1-5)) &
     	-(UTM%M(6*N2-5)-UTM%L(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%L(6*N1-4)) &
     	-(UTM%M(6*N2-4)-UTM%L(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%L(6*N1-3)) &
     	-(UTM%M(6*N2-3)-UTM%L(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%L(6*N1-2)) &
     	-(UTM%M(6*N2-2)-UTM%L(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%L(6*N1-1)) &
     	-(UTM%M(6*N2-1)-UTM%L(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%L(6*N1-0)) &
     	-(UTM%M(6*N2-0)-UTM%L(6*N1-0)))**2
	ENDIF
	ENDDO
	ENDIF

        IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
	DO X=1,NTOT%R
	IF (EFS % R(X).NE.0.0) THEN
	N1=NANS % R(1,X)
	N2=NANS % R(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%R(6*N1-5)) &
     	-(UTM%M(6*N2-5)-UTM%R(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%R(6*N1-4)) &
     	-(UTM%M(6*N2-4)-UTM%R(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%R(6*N1-3)) &
     	-(UTM%M(6*N2-3)-UTM%R(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%R(6*N1-2)) &
     	-(UTM%M(6*N2-2)-UTM%R(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%R(6*N1-1)) &
     	-(UTM%M(6*N2-1)-UTM%R(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%R(6*N1-0)) &
     	-(UTM%M(6*N2-0)-UTM%R(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%R(6*N1-5)) &
     	-(UTM%M(6*N2-5)-UTM%R(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%R(6*N1-4)) &
     	-(UTM%M(6*N2-4)-UTM%R(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%R(6*N1-3)) &
     	-(UTM%M(6*N2-3)-UTM%R(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%R(6*N1-2)) &
     	-(UTM%M(6*N2-2)-UTM%R(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%R(6*N1-1)) &
     	-(UTM%M(6*N2-1)-UTM%R(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%R(6*N1-0)) &
     	-(UTM%M(6*N2-0)-UTM%R(6*N1-0)))**2
	ENDIF
	ENDDO
	ENDIF

	IF (myid.lt.(YN-1)*ntasks/YN) THEN
        DO X=1,NTOT%F
        IF (EFS % F(X).NE.0.0) THEN
	N1=NANS % F(1,X)
        N2=NANS % F(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%F(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%F(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%F(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%F(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%F(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%F(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%F(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%F(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%F(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%F(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%F(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%F(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%F(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%F(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%F(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%F(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%F(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%F(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%F(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%F(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%F(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%F(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%F(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%F(6*N1-0)))**2
	ENDIF
	ENDDO
        ENDIF

	IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,NTOT%FL
        IF (EFS % FL(X).NE.0.0) THEN
	N1=NANS % FL(1,X)
        N2=NANS % FL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%FL(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%FL(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%FL(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%FL(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%FL(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%FL(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%FL(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%FL(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%FL(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%FL(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%FL(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%FL(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%FL(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%FL(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%FL(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%FL(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%FL(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%FL(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%FL(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%FL(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%FL(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%FL(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%FL(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%FL(6*N1-0)))**2
	ENDIF
	ENDDO
        ENDIF

	IF (myid.lt.(YN-1)*ntasks/YN &
     	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,NTOT%FR
        IF (EFS % FR(X).NE.0.0) THEN
	N1=NANS % FR(1,X)
        N2=NANS % FR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%FR(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%FR(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%FR(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%FR(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%FR(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%FR(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%FR(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%FR(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%FR(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%FR(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%FR(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%FR(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%FR(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%FR(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%FR(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%FR(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%FR(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%FR(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%FR(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%FR(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%FR(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%FR(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%FR(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%FR(6*N1-0)))**2
	ENDIF
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN) THEN
        DO X=1,NTOT%B
        IF (EFS % B(X).NE.0.0) THEN
	N1=NANS % B(1,X)
        N2=NANS % B(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%B(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%B(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%B(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%B(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%B(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%B(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%B(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%B(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%B(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%B(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%B(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%B(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%B(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%B(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%B(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%B(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%B(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%B(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%B(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%B(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%B(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%B(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%B(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%B(6*N1-0)))**2
	ENDIF
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,NTOT%BL
        IF (EFS % BL(X).NE.0.0) THEN
	N1=NANS % BL(1,X)
        N2=NANS % BL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%BL(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%BL(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%BL(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%BL(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%BL(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%BL(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%BL(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%BL(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%BL(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%BL(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%BL(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%BL(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%BL(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%BL(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%BL(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%BL(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%BL(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%BL(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%BL(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%BL(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%BL(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%BL(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%BL(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%BL(6*N1-0)))**2
	ENDIF
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,NTOT%BR
        IF (EFS % BR(X).NE.0.0) THEN
	N1=NANS % BR(1,X)
        N2=NANS % BR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%BR(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%BR(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%BR(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%BR(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%BR(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%BR(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%BR(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%BR(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%BR(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%BR(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%BR(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%BR(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%BR(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%BR(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%BR(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%BR(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%BR(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%BR(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%BR(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%BR(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%BR(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%BR(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%BR(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%BR(6*N1-0)))**2
	ENDIF
	ENDDO
        ENDIF

!-----------------------------------------------------------------
	DO X=1,FXC
	N1=FXF(1,X)
	N2=FXF(2,X)
	D(6*N1-5)=D(6*N1-5)+(DMP/DT)*((UT%M(6*N1-5)-UT%M(6*N2-5)) &
     	-(UTM%M(6*N1-5)-UTM%M(6*N2-5)))
	D(6*N1-4)=D(6*N1-4)+(DMP/DT)*((UT%M(6*N1-4)-UT%M(6*N2-4)) &
     	-(UTM%M(6*N1-4)-UTM%M(6*N2-4)))
	D(6*N1-3)=D(6*N1-3)+(DMP/DT)*((UT%M(6*N1-3)-UT%M(6*N2-3)) &
     	-(UTM%M(6*N1-3)-UTM%M(6*N2-3)))
	D(6*N1-2)=D(6*N1-2)+(DMP2/DT)*((UT%M(6*N1-2)-UT%M(6*N2-2)) &
     	-(UTM%M(6*N1-2)-UTM%M(6*N2-2)))
	D(6*N1-1)=D(6*N1-1)+(DMP2/DT)*((UT%M(6*N1-1)-UT%M(6*N2-1)) &
     	-(UTM%M(6*N1-1)-UTM%M(6*N2-1)))
	D(6*N1-0)=D(6*N1-0)+(DMP2/DT)*((UT%M(6*N1-0)-UT%M(6*N2-0)) &
     	-(UTM%M(6*N1-0)-UTM%M(6*N2-0)))
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%M(6*N1-5)) &
     	-(UTM%M(6*N2-5)-UTM%M(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%M(6*N1-4)) &
     	-(UTM%M(6*N2-4)-UTM%M(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%M(6*N1-3)) &
     	-(UTM%M(6*N2-3)-UTM%M(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%M(6*N1-2)) &
     	-(UTM%M(6*N2-2)-UTM%M(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%M(6*N1-1)) &
     	-(UTM%M(6*N2-1)-UTM%M(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%M(6*N1-0)) &
     	-(UTM%M(6*N2-0)-UTM%M(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%M(6*N1-5)) &
     	-(UTM%M(6*N2-5)-UTM%M(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%M(6*N1-4)) &
     	-(UTM%M(6*N2-4)-UTM%M(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%M(6*N1-3)) &
     	-(UTM%M(6*N2-3)-UTM%M(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%M(6*N1-2)) &
     	-(UTM%M(6*N2-2)-UTM%M(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%M(6*N1-1)) &
     	-(UTM%M(6*N2-1)-UTM%M(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%M(6*N1-0)) &
     	-(UTM%M(6*N2-0)-UTM%M(6*N1-0)))**2
	ENDDO

        IF (MOD(myid,ntasks/YN).ne.0) THEN
	DO X=1,FXL
	N1=FXFL(1,X)
	N2=FXFL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%L(6*N1-5)) &
      	-(UTM%M(6*N2-5)-UTM%L(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%L(6*N1-4)) &
      	-(UTM%M(6*N2-4)-UTM%L(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%L(6*N1-3)) &
      	-(UTM%M(6*N2-3)-UTM%L(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%L(6*N1-2)) &
      	-(UTM%M(6*N2-2)-UTM%L(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%L(6*N1-1)) &
      	-(UTM%M(6*N2-1)-UTM%L(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%L(6*N1-0)) &
      	-(UTM%M(6*N2-0)-UTM%L(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%L(6*N1-5)) &
      	-(UTM%M(6*N2-5)-UTM%L(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%L(6*N1-4)) &
      	-(UTM%M(6*N2-4)-UTM%L(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%L(6*N1-3)) &
      	-(UTM%M(6*N2-3)-UTM%L(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%L(6*N1-2)) &
      	-(UTM%M(6*N2-2)-UTM%L(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%L(6*N1-1)) &
      	-(UTM%M(6*N2-1)-UTM%L(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%L(6*N1-0)) &
      	-(UTM%M(6*N2-0)-UTM%L(6*N1-0)))**2
	ENDDO
	ENDIF

        IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
	DO X=1,FXR
	N1=FXFR(1,X)
	N2=FXFR(2,X)
 	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%R(6*N1-5)) &
      	-(UTM%M(6*N2-5)-UTM%R(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%R(6*N1-4)) &
      	-(UTM%M(6*N2-4)-UTM%R(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%R(6*N1-3)) &
      	-(UTM%M(6*N2-3)-UTM%R(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%R(6*N1-2)) &
      	-(UTM%M(6*N2-2)-UTM%R(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%R(6*N1-1)) &
      	-(UTM%M(6*N2-1)-UTM%R(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%R(6*N1-0)) &
      	-(UTM%M(6*N2-0)-UTM%R(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%R(6*N1-5)) &
      	-(UTM%M(6*N2-5)-UTM%R(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%R(6*N1-4)) &
      	-(UTM%M(6*N2-4)-UTM%R(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%R(6*N1-3)) &
      	-(UTM%M(6*N2-3)-UTM%R(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%R(6*N1-2)) &
      	-(UTM%M(6*N2-2)-UTM%R(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%R(6*N1-1)) &
      	-(UTM%M(6*N2-1)-UTM%R(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%R(6*N1-0)) &
      	-(UTM%M(6*N2-0)-UTM%R(6*N1-0)))**2
	ENDDO
	ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN) THEN
        DO X=1,FXCF
        N1=FXFF(1,X)
        N2=FXFF(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%F(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%F(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%F(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%F(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%F(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%F(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%F(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%F(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%F(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%F(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%F(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%F(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%F(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%F(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%F(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%F(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%F(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%F(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%F(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%F(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%F(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%F(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%F(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%F(6*N1-0)))**2
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,FXCFL
        N1=FXFFL(1,X)
        N2=FXFFL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%FL(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%FL(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%FL(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%FL(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%FL(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%FL(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%FL(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%FL(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%FL(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%FL(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%FL(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%FL(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%FL(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%FL(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%FL(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%FL(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%FL(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%FL(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%FL(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%FL(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%FL(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%FL(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%FL(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%FL(6*N1-0)))**2
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN &
      	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,FXCFR
        N1=FXFFR(1,X)
        N2=FXFFR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%FR(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%FR(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%FR(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%FR(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%FR(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%FR(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%FR(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%FR(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%FR(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%FR(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%FR(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%FR(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%FR(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%FR(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%FR(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%FR(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%FR(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%FR(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%FR(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%FR(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%FR(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%FR(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%FR(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%FR(6*N1-0)))**2
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN) THEN
        DO X=1,FXCB
	N1=FXFB(1,X)
        N2=FXFB(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%B(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%B(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%B(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%B(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%B(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%B(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%B(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%B(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%B(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%B(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%B(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%B(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%B(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%B(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%B(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%B(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%B(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%B(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%B(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%B(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%B(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%B(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%B(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%B(6*N1-0)))**2
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,FXCBL
	N1=FXFBL(1,X)
        N2=FXFBL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%BL(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%BL(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%BL(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%BL(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%BL(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%BL(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%BL(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%BL(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%BL(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%BL(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%BL(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%BL(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%BL(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%BL(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%BL(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%BL(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%BL(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%BL(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%BL(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%BL(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%BL(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%BL(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%BL(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%BL(6*N1-0)))**2
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN &
     	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,FXCBR
	N1=FXFBR(1,X)
        N2=FXFBR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%BR(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%BR(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%BR(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%BR(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%BR(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%BR(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%BR(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%BR(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%BR(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%BR(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%BR(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%BR(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%BR(6*N1-5)) &
        -(UTM%M(6*N2-5)-UTM%BR(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%BR(6*N1-4)) &
        -(UTM%M(6*N2-4)-UTM%BR(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%BR(6*N1-3)) &
        -(UTM%M(6*N2-3)-UTM%BR(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%BR(6*N1-2)) &
        -(UTM%M(6*N2-2)-UTM%BR(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%BR(6*N1-1)) &
        -(UTM%M(6*N2-1)-UTM%BR(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%BR(6*N1-0)) &
        -(UTM%M(6*N2-0)-UTM%BR(6*N1-0)))**2
	ENDDO
        ENDIF

!---------------------------------------------------------

        DO X=1,NN
        F(6*X-5)= -(VDP(X)/(DT*2))*UTM%M(6*X-5)
        F(6*X-4)= -(VDP(X)/(DT*2))*UTM%M(6*X-4)
        F(6*X-3)= -(VDP(X)/(DT*2))*UTM%M(6*X-3)
        F(6*X-2)= -(VDP(X)/(DT*2))*UTM%M(6*X-2)
        F(6*X-1)= -(VDP(X)/(DT*2))*UTM%M(6*X-1)
        F(6*X-0)= -(VDP(X)/(DT*2))*UTM%M(6*X-0)
	END DO

	DO X=1,6*NN
	R(X)=-A(X)-C(X)-D(X)-F(X)
	EN(X)=EN(X)+A(X)*(UT%M(X)-UTM%M(X))
	END DO

	RETURN
END SUBROUTINE EFFLOAD
