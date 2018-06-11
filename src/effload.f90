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
     FXF,FXC,VDP,DPE,EF,EFL,EFR,NAN,NRXF,MFIL,CT,NANL,NANR,NTOL, &
     NTOR,NRXFL,NRXFR,UTL,UTR,myid,ntasks,FXL,FXFL,FXR,FXFR,L, &
     NTOF,NTOB,NRXFF,NRXFB,UTF,UTB,EFF,EFB, &
     FXCF,FXFF,FXCB,FXFB,NANB,NANF,PNN,YN, &
     NANFL,NANFR,NANBL,NANBR, &
     NRXFFL,NRXFFR,NRXFBL,NRXFBR, &
     UTFL,UTFR,UTBL,UTBR, &
     EFFL,EFFR,EFBL,EFBR, &
     NTOFL,NTOFR,NTOBL,NTOBR, &
     FXCFL,FXFFL,FXCBL,FXFBL, &
     FXCFR,FXFFR,FXCBR,FXFBR)

        USE INOUT

	IMPLICIT NONE
        include 'mpif.h'
        include 'param.dat'
        REAL*8 MFIL(NOMA),NRXF(3,NOMA),VDP(NOMA)
        REAL*8 NRXFL(3,NOMA),NRXFR(3,NOMA)
        REAL*8 NRXFF(3,NOMA),NRXFB(3,NOMA)
        REAL*8 NRXFFL(3,NOMA),NRXFBL(3,NOMA)
        REAL*8 NRXFFR(3,NOMA),NRXFBR(3,NOMA)
        REAL*8 UTF(NODM),UTB(NODM)
        REAL*8 UTFL(NODM),UTBL(NODM)
        REAL*8 UTFR(NODM),UTBR(NODM)
	REAL*8 A(NODM),C(NODM),F(NODM),D(NODM)
	REAL*8 UT(NODM),UTM(NODM),EF(NOCON),EFL(NOCON)
	REAL*8 UTL(NODM),UTML(NODM),EFR(NOCON)
	REAL*8 UTR(NODM),UTMR(NODM),LNN,UTMB(NODM),UTMF(NODM)
	REAL*8 UTMBL(NODM),UTMFL(NODM),UTMBR(NODM),UTMFR(NODM)
	REAL*8 DUT(NODM),CT(NODC),EN(NODM),R(NODM)
	INTEGER FXF(2,NODC),FXFL(2,NODC),FXFR(2,NODC)
	INTEGER FXFF(2,NODC),FXFB(2,NODC),YN
        INTEGER FXFFL(2,NODC),FXFBL(2,NODC)
        INTEGER FXFFR(2,NODC),FXFBR(2,NODC)
	REAL*8 DPE,S,E,EFB(NOCON),EFF(NOCON)
	REAL*8 EFBL(NOCON),EFFL(NOCON)
	REAL*8 EFBR(NOCON),EFFR(NOCON)
	REAL*8 T,DT,M,JS,L,ALF,DMP
	REAL*8 G,X1,Y1,Z1,X2,Y2,Z2,TT(12,12)
	REAL*8 DX1,DY1,DZ1,DX2,DY2,DZ2,DMP2
	REAL*8 DP,DP2
	INTEGER NTOT,N,NL,NB,N1,N2,X,XL,XR,PNN(0:5000)
	INTEGER I,J,NN,RY,FXC,FXL,FXR,FXCF,FXCB,NAN(3,NOCON)
	INTEGER FXCFL,FXCBL,FXCFR,FXCBR
	INTEGER NANL(3,NOMA),NANR(3,NOMA),NTOL,NTOR
        INTEGER NANF(3,NOMA),NANB(3,NOMA),NTOF,NTOB
	INTEGER NTOFL,NTOFR,NTOBL,NTOBR
        INTEGER NANFL(3,NOMA),NANBL(3,NOMA)
        INTEGER NANFR(3,NOMA),NANBR(3,NOMA)
        INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),comm
        INTEGER myid,ntasks,ierr

	DO I=1,6*NN
	A(I)=0.0
	D(I)=0.0
        END DO

 	DO X=1,NTOT
	N1=NAN(1,X)
	N2=NAN(2,X)
	X1=NRXF(1,N1)
	Y1=NRXF(2,N1)
	Z1=NRXF(3,N1)
	X2=NRXF(1,N2)
	Y2=NRXF(2,N2)
	Z2=NRXF(3,N2)


	DX1=UT(6*N1-5)
	DY1=UT(6*N1-4)
	DZ1=UT(6*N1-3)
	DX2=UT(6*N2-5)
	DY2=UT(6*N2-4)
	DZ2=UT(6*N2-3)

        DUT(6*N1-5)=UT(6*N1-5)-UTM(6*N1-5)
        DUT(6*N1-4)=UT(6*N1-4)-UTM(6*N1-4)
        DUT(6*N1-3)=UT(6*N1-3)-UTM(6*N1-3)
        DUT(6*N1-2)=UT(6*N1-2)-UTM(6*N1-2)
        DUT(6*N1-1)=UT(6*N1-1)-UTM(6*N1-1)
        DUT(6*N1-0)=UT(6*N1-0)-UTM(6*N1-0)
                                                        
        DUT(6*N2-5)=UT(6*N2-5)-UTM(6*N2-5)
        DUT(6*N2-4)=UT(6*N2-4)-UTM(6*N2-4)
        DUT(6*N2-3)=UT(6*N2-3)-UTM(6*N2-3)
        DUT(6*N2-2)=UT(6*N2-2)-UTM(6*N2-2)
        DUT(6*N2-1)=UT(6*N2-1)-UTM(6*N2-1)
        DUT(6*N2-0)=UT(6*N2-0)-UTM(6*N2-0)

        IF (EF(X).NE.0.0) THEN
        CALL AMAT(EF(X),S,EF(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
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
     CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION, &
     dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (MOD(myid,ntasks/YN).ne.0) &
     CALL MPI_Recv(UTML,6*PNN(source),MPI_DOUBLE_PRECISION, &
     source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-1
      source=myid+1
      tag=135
      IF (MOD(myid,ntasks/YN).ne.0) &
     CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION, &
     dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) &
     CALL MPI_Recv(UTMR,6*PNN(source),MPI_DOUBLE_PRECISION, &
     source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid+ntasks/YN
      source=myid-ntasks/YN
      tag=137
      IF (myid.lt.(YN-1)*ntasks/YN) &
     CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION, &
     dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.ge.ntasks/YN) &
     CALL MPI_Recv(UTMB,6*PNN(source),MPI_DOUBLE_PRECISION, &
     source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid+ntasks/YN+1
      source=myid-ntasks/YN-1
      tag=139
      IF (myid.lt.(YN-1)*ntasks/YN &
     .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) &
     CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION, &
     dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) &
     CALL MPI_Recv(UTMBL,6*PNN(source),MPI_DOUBLE_PRECISION, &
     source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid+ntasks/YN-1
      source=myid-ntasks/YN+1
      tag=141
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) &
           CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION, &
           dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) &
           CALL MPI_Recv(UTMBR,6*PNN(source),MPI_DOUBLE_PRECISION, &
           source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-ntasks/YN
      source=myid+ntasks/YN
      tag=143
      IF (myid.ge.ntasks/YN) &
      CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION, &
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN) &
      CALL MPI_Recv(UTMF,6*PNN(source),MPI_DOUBLE_PRECISION, &
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-ntasks/YN+1
      source=myid+ntasks/YN-1
      tag=145
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) &
      CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION, &
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) &
      CALL MPI_Recv(UTMFL,6*PNN(source),MPI_DOUBLE_PRECISION, &
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-ntasks/YN-1
      source=myid+ntasks/YN+1
      tag=147
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) &
      CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION, &
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN &
      .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) &
      CALL MPI_Recv(UTMFR,6*PNN(source),MPI_DOUBLE_PRECISION, &
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

!------------------------------------------------------
        IF (MOD(myid,ntasks/YN).ne.0) THEN
  	DO X=1,NTOL
	XL=NTOT+X
	N1=NANL(1,X)
	N2=NANL(2,X)
	X1=NRXFL(1,N1)
	Y1=NRXFL(2,N1)
	Z1=NRXFL(3,N1)
	X2=NRXF(1,N2)
	Y2=NRXF(2,N2)
	Z2=NRXF(3,N2)
	DX1=UTL(6*N1-5)
	DY1=UTL(6*N1-4)
	DZ1=UTL(6*N1-3)
	DX2=UT(6*N2-5)
	DY2=UT(6*N2-4)
	DZ2=UT(6*N2-3)
        DUT(6*N1-5)=UTL(6*N1-5)-UTML(6*N1-5)
        DUT(6*N1-4)=UTL(6*N1-4)-UTML(6*N1-4)
        DUT(6*N1-3)=UTL(6*N1-3)-UTML(6*N1-3)
        DUT(6*N1-2)=UTL(6*N1-2)-UTML(6*N1-2)
        DUT(6*N1-1)=UTL(6*N1-1)-UTML(6*N1-1)
        DUT(6*N1-0)=UTL(6*N1-0)-UTML(6*N1-0)
        DUT(6*N2-5)=UT(6*N2-5)-UTM(6*N2-5)
        DUT(6*N2-4)=UT(6*N2-4)-UTM(6*N2-4)
        DUT(6*N2-3)=UT(6*N2-3)-UTM(6*N2-3)
        DUT(6*N2-2)=UT(6*N2-2)-UTM(6*N2-2)
        DUT(6*N2-1)=UT(6*N2-1)-UTM(6*N2-1)
        DUT(6*N2-0)=UT(6*N2-0)-UTM(6*N2-0)
        IF (EFL(X).NE.0.0) THEN
        CALL AMAT(EFL(X),S,EFL(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
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
  	DO X=1,NTOR
	XR=NTOT+NTOL+X
	N1=NANR(1,X)
	N2=NANR(2,X)
	X1=NRXFR(1,N1)
	Y1=NRXFR(2,N1)
	Z1=NRXFR(3,N1)
	X2=NRXF(1,N2)
	Y2=NRXF(2,N2)
	Z2=NRXF(3,N2)
	DX1=UTR(6*N1-5)
	DY1=UTR(6*N1-4)
	DZ1=UTR(6*N1-3)
	DX2=UT(6*N2-5)
	DY2=UT(6*N2-4)
	DZ2=UT(6*N2-3)
        DUT(6*N1-5)=UTR(6*N1-5)-UTMR(6*N1-5)
        DUT(6*N1-4)=UTR(6*N1-4)-UTMR(6*N1-4)
        DUT(6*N1-3)=UTR(6*N1-3)-UTMR(6*N1-3)
        DUT(6*N1-2)=UTR(6*N1-2)-UTMR(6*N1-2)
        DUT(6*N1-1)=UTR(6*N1-1)-UTMR(6*N1-1)
        DUT(6*N1-0)=UTR(6*N1-0)-UTMR(6*N1-0)
        DUT(6*N2-5)=UT(6*N2-5)-UTM(6*N2-5)
        DUT(6*N2-4)=UT(6*N2-4)-UTM(6*N2-4)
        DUT(6*N2-3)=UT(6*N2-3)-UTM(6*N2-3)
        DUT(6*N2-2)=UT(6*N2-2)-UTM(6*N2-2)
        DUT(6*N2-1)=UT(6*N2-1)-UTM(6*N2-1)
        DUT(6*N2-0)=UT(6*N2-0)-UTM(6*N2-0)
        IF (EFR(X).NE.0.0) THEN
        CALL AMAT(EFR(X),S,EFR(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
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
        DO X=1,NTOF
        XR=NTOT+NTOL+NTOR+X
        N1=NANF(1,X)
        N2=NANF(2,X)
        X1=NRXFF(1,N1)
        Y1=NRXFF(2,N1)
        Z1=NRXFF(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
        DX1=UTF(6*N1-5)
        DY1=UTF(6*N1-4)
        DZ1=UTF(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
        DUT(6*N1-5)=UTF(6*N1-5)-UTMF(6*N1-5)
        DUT(6*N1-4)=UTF(6*N1-4)-UTMF(6*N1-4)
        DUT(6*N1-3)=UTF(6*N1-3)-UTMF(6*N1-3)
        DUT(6*N1-2)=UTF(6*N1-2)-UTMF(6*N1-2)
        DUT(6*N1-1)=UTF(6*N1-1)-UTMF(6*N1-1)
        DUT(6*N1-0)=UTF(6*N1-0)-UTMF(6*N1-0)
        DUT(6*N2-5)=UT(6*N2-5)-UTM(6*N2-5)
        DUT(6*N2-4)=UT(6*N2-4)-UTM(6*N2-4)
        DUT(6*N2-3)=UT(6*N2-3)-UTM(6*N2-3)
        DUT(6*N2-2)=UT(6*N2-2)-UTM(6*N2-2)
        DUT(6*N2-1)=UT(6*N2-1)-UTM(6*N2-1)
        DUT(6*N2-0)=UT(6*N2-0)-UTM(6*N2-0)
        IF (EFF(X).NE.0.0) THEN
        CALL AMAT(EFF(X),S,EFF(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
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
        DO X=1,NTOFL
        XR=NTOT+NTOL+NTOR+NTOF+X
        N1=NANFL(1,X)
        N2=NANFL(2,X)
        X1=NRXFFL(1,N1)
        Y1=NRXFFL(2,N1)
        Z1=NRXFFL(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
        DX1=UTFL(6*N1-5)
        DY1=UTFL(6*N1-4)
        DZ1=UTFL(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
        DUT(6*N1-5)=UTFL(6*N1-5)-UTMFL(6*N1-5)
        DUT(6*N1-4)=UTFL(6*N1-4)-UTMFL(6*N1-4)
        DUT(6*N1-3)=UTFL(6*N1-3)-UTMFL(6*N1-3)
        DUT(6*N1-2)=UTFL(6*N1-2)-UTMFL(6*N1-2)
        DUT(6*N1-1)=UTFL(6*N1-1)-UTMFL(6*N1-1)
        DUT(6*N1-0)=UTFL(6*N1-0)-UTMFL(6*N1-0)
        DUT(6*N2-5)=UT(6*N2-5)-UTM(6*N2-5)
        DUT(6*N2-4)=UT(6*N2-4)-UTM(6*N2-4)
        DUT(6*N2-3)=UT(6*N2-3)-UTM(6*N2-3)
        DUT(6*N2-2)=UT(6*N2-2)-UTM(6*N2-2)
        DUT(6*N2-1)=UT(6*N2-1)-UTM(6*N2-1)
        DUT(6*N2-0)=UT(6*N2-0)-UTM(6*N2-0)
        IF (EFFL(X).NE.0.0) THEN
        CALL AMAT(EFFL(X),S,EFFL(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
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
        DO X=1,NTOFR
        XR=NTOT+NTOL+NTOR+NTOF+NTOFL+X
        N1=NANFR(1,X)
        N2=NANFR(2,X)
        X1=NRXFFR(1,N1)
        Y1=NRXFFR(2,N1)
        Z1=NRXFFR(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
        DX1=UTFR(6*N1-5)
        DY1=UTFR(6*N1-4)
        DZ1=UTFR(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
        DUT(6*N1-5)=UTFR(6*N1-5)-UTMFR(6*N1-5)
        DUT(6*N1-4)=UTFR(6*N1-4)-UTMFR(6*N1-4)
        DUT(6*N1-3)=UTFR(6*N1-3)-UTMFR(6*N1-3)
        DUT(6*N1-2)=UTFR(6*N1-2)-UTMFR(6*N1-2)
        DUT(6*N1-1)=UTFR(6*N1-1)-UTMFR(6*N1-1)
        DUT(6*N1-0)=UTFR(6*N1-0)-UTMFR(6*N1-0)
        DUT(6*N2-5)=UT(6*N2-5)-UTM(6*N2-5)
        DUT(6*N2-4)=UT(6*N2-4)-UTM(6*N2-4)
        DUT(6*N2-3)=UT(6*N2-3)-UTM(6*N2-3)
        DUT(6*N2-2)=UT(6*N2-2)-UTM(6*N2-2)
        DUT(6*N2-1)=UT(6*N2-1)-UTM(6*N2-1)
        DUT(6*N2-0)=UT(6*N2-0)-UTM(6*N2-0)
        IF (EFFR(X).NE.0.0) THEN
        CALL AMAT(EFFR(X),S,EFFR(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
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
        DO X=1,NTOB
        XR=NTOT+NTOL+NTOR+NTOF+NTOFL+NTOFR+X
        N1=NANB(1,X)
        N2=NANB(2,X)
        X1=NRXFB(1,N1)
        Y1=NRXFB(2,N1)
        Z1=NRXFB(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
	DX1=UTB(6*N1-5)
        DY1=UTB(6*N1-4)
        DZ1=UTB(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
	DUT(6*N1-5)=UTB(6*N1-5)-UTMB(6*N1-5)
        DUT(6*N1-4)=UTB(6*N1-4)-UTMB(6*N1-4)
        DUT(6*N1-3)=UTB(6*N1-3)-UTMB(6*N1-3)
        DUT(6*N1-2)=UTB(6*N1-2)-UTMB(6*N1-2)
        DUT(6*N1-1)=UTB(6*N1-1)-UTMB(6*N1-1)
        DUT(6*N1-0)=UTB(6*N1-0)-UTMB(6*N1-0)
	DUT(6*N2-5)=UT(6*N2-5)-UTM(6*N2-5)
        DUT(6*N2-4)=UT(6*N2-4)-UTM(6*N2-4)
        DUT(6*N2-3)=UT(6*N2-3)-UTM(6*N2-3)
        DUT(6*N2-2)=UT(6*N2-2)-UTM(6*N2-2)
        DUT(6*N2-1)=UT(6*N2-1)-UTM(6*N2-1)
        DUT(6*N2-0)=UT(6*N2-0)-UTM(6*N2-0)
	IF (EFB(X).NE.0.0) THEN
        CALL AMAT(EFB(X),S,EFB(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
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
        DO X=1,NTOBL
        XR=NTOT+NTOL+NTOR+NTOF+NTOFL+NTOFR+NTOB+X
        N1=NANBL(1,X)
        N2=NANBL(2,X)
        X1=NRXFBL(1,N1)
        Y1=NRXFBL(2,N1)
        Z1=NRXFBL(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
	DX1=UTBL(6*N1-5)
        DY1=UTBL(6*N1-4)
        DZ1=UTBL(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
	DUT(6*N1-5)=UTBL(6*N1-5)-UTMBL(6*N1-5)
        DUT(6*N1-4)=UTBL(6*N1-4)-UTMBL(6*N1-4)
        DUT(6*N1-3)=UTBL(6*N1-3)-UTMBL(6*N1-3)
        DUT(6*N1-2)=UTBL(6*N1-2)-UTMBL(6*N1-2)
        DUT(6*N1-1)=UTBL(6*N1-1)-UTMBL(6*N1-1)
        DUT(6*N1-0)=UTBL(6*N1-0)-UTMBL(6*N1-0)
	DUT(6*N2-5)=UT(6*N2-5)-UTM(6*N2-5)
        DUT(6*N2-4)=UT(6*N2-4)-UTM(6*N2-4)
        DUT(6*N2-3)=UT(6*N2-3)-UTM(6*N2-3)
        DUT(6*N2-2)=UT(6*N2-2)-UTM(6*N2-2)
        DUT(6*N2-1)=UT(6*N2-1)-UTM(6*N2-1)
        DUT(6*N2-0)=UT(6*N2-0)-UTM(6*N2-0)
	IF (EFBL(X).NE.0.0) THEN
        CALL AMAT(EFBL(X),S,EFBL(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
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
        DO X=1,NTOBR
        XR=NTOT+NTOL+NTOR+NTOF+NTOFL+NTOFR+NTOB+NTOBL+X
        N1=NANBR(1,X)
        N2=NANBR(2,X)
        X1=NRXFBR(1,N1)
        Y1=NRXFBR(2,N1)
        Z1=NRXFBR(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
	DX1=UTBR(6*N1-5)
        DY1=UTBR(6*N1-4)
        DZ1=UTBR(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
	DUT(6*N1-5)=UTBR(6*N1-5)-UTMBR(6*N1-5)
        DUT(6*N1-4)=UTBR(6*N1-4)-UTMBR(6*N1-4)
        DUT(6*N1-3)=UTBR(6*N1-3)-UTMBR(6*N1-3)
        DUT(6*N1-2)=UTBR(6*N1-2)-UTMBR(6*N1-2)
        DUT(6*N1-1)=UTBR(6*N1-1)-UTMBR(6*N1-1)
        DUT(6*N1-0)=UTBR(6*N1-0)-UTMBR(6*N1-0)
	DUT(6*N2-5)=UT(6*N2-5)-UTM(6*N2-5)
        DUT(6*N2-4)=UT(6*N2-4)-UTM(6*N2-4)
        DUT(6*N2-3)=UT(6*N2-3)-UTM(6*N2-3)
        DUT(6*N2-2)=UT(6*N2-2)-UTM(6*N2-2)
        DUT(6*N2-1)=UT(6*N2-1)-UTM(6*N2-1)
        DUT(6*N2-0)=UT(6*N2-0)-UTM(6*N2-0)
	IF (EFBR(X).NE.0.0) THEN
        CALL AMAT(EFBR(X),S,EFBR(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
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

        DO X=1,NTOT
        N1=NAN(1,X)
        N2=NAN(2,X)
        X1=NRXF(1,N1)
        Y1=NRXF(2,N1)
        Z1=NRXF(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
        DX1=UT(6*N1-5)
        DY1=UT(6*N1-4)
        DZ1=UT(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
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
        DO X=1,NTOL
	XL=NTOT+X
        N1=NANL(1,X)
        N2=NANL(2,X)
        X1=NRXFL(1,N1)
        Y1=NRXFL(2,N1)
        Z1=NRXFL(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
        DX1=UTL(6*N1-5)
        DY1=UTL(6*N1-4)
        DZ1=UTL(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
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
        DO X=1,NTOR
 	XR=NTOT+NTOL+X
        N1=NANR(1,X)
        N2=NANR(2,X)
        X1=NRXFR(1,N1)
        Y1=NRXFR(2,N1)
        Z1=NRXFR(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
        DX1=UTR(6*N1-5)
        DY1=UTR(6*N1-4)
        DZ1=UTR(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
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
        DO X=1,NTOF
        XR=NTOT+NTOL+NTOR+X
        N1=NANF(1,X)
        N2=NANF(2,X)
        X1=NRXFF(1,N1)
        Y1=NRXFF(2,N1)
        Z1=NRXFF(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
        DX1=UTF(6*N1-5)
        DY1=UTF(6*N1-4)
        DZ1=UTF(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
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
        DO X=1,NTOFL
        XR=NTOT+NTOL+NTOR+NTOF+X
        N1=NANFL(1,X)
        N2=NANFL(2,X)
        X1=NRXFFL(1,N1)
        Y1=NRXFFL(2,N1)
        Z1=NRXFFL(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
        DX1=UTFL(6*N1-5)
        DY1=UTFL(6*N1-4)
        DZ1=UTFL(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
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
        DO X=1,NTOFR
        XR=NTOT+NTOL+NTOR+NTOF+NTOFL+X
        N1=NANFR(1,X)
        N2=NANFR(2,X)
        X1=NRXFFR(1,N1)
        Y1=NRXFFR(2,N1)
        Z1=NRXFFR(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
        DX1=UTFR(6*N1-5)
        DY1=UTFR(6*N1-4)
        DZ1=UTFR(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
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
        DO X=1,NTOB
        XR=NTOT+NTOL+NTOR+NTOF+NTOFL+NTOFR+X
	N1=NANB(1,X)
        N2=NANB(2,X)
	X1=NRXFB(1,N1)
        Y1=NRXFB(2,N1)
        Z1=NRXFB(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
	DX1=UTB(6*N1-5)
        DY1=UTB(6*N1-4)
        DZ1=UTB(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
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
        DO X=1,NTOBL
        XR=NTOT+NTOL+NTOR+NTOF+NTOFL+NTOFR+NTOB+X
	N1=NANBL(1,X)
        N2=NANBL(2,X)
	X1=NRXFBL(1,N1)
        Y1=NRXFBL(2,N1)
        Z1=NRXFBL(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
	DX1=UTBL(6*N1-5)
        DY1=UTBL(6*N1-4)
        DZ1=UTBL(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
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
        DO X=1,NTOBR
        XR=NTOT+NTOL+NTOR+NTOF+NTOFL+NTOFR+NTOB+NTOBL+X
	N1=NANBR(1,X)
        N2=NANBR(2,X)
	X1=NRXFBR(1,N1)
        Y1=NRXFBR(2,N1)
        Z1=NRXFBR(3,N1)
        X2=NRXF(1,N2)
        Y2=NRXF(2,N2)
        Z2=NRXF(3,N2)
	DX1=UTBR(6*N1-5)
        DY1=UTBR(6*N1-4)
        DZ1=UTBR(6*N1-3)
        DX2=UT(6*N2-5)
        DY2=UT(6*N2-4)
        DZ2=UT(6*N2-3)
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
	C(6*X-5)= (MFIL(X)/DT**2)*UTM(6*X-5)-(2*MFIL(X)/DT**2)*UT(6*X-5)
	C(6*X-4)= (MFIL(X)/DT**2)*UTM(6*X-4)-(2*MFIL(X)/DT**2)*UT(6*X-4)
	C(6*X-3)= (MFIL(X)/DT**2)*UTM(6*X-3)-(2*MFIL(X)/DT**2)*UT(6*X-3)
	C(6*X-2)= ((MFIL(X)*JS/M)/DT**2)* &
      	UTM(6*X-2)-(2*(MFIL(X)*JS/M)/DT**2)*UT(6*X-2)
	C(6*X-1)= ((MFIL(X)*JS/M)/DT**2)* &
      	UTM(6*X-1)-(2*(MFIL(X)*JS/M)/DT**2)*UT(6*X-1)
	C(6*X-0)= ((MFIL(X)*JS/M)/DT**2)* &
      	UTM(6*X-0)-(2*(MFIL(X)*JS/M)/DT**2)*UT(6*X-0)
	ENDDO

!------------------------------------------------------------------------


	DO X=1,NTOT
	IF (EF(X).NE.0.0) THEN
	N1=NAN(1,X)
	N2=NAN(2,X)
	D(6*N1-5)=D(6*N1-5)+(DMP/DT)*((UT(6*N1-5)-UT(6*N2-5)) &
      	-(UTM(6*N1-5)-UTM(6*N2-5)))
	D(6*N1-4)=D(6*N1-4)+(DMP/DT)*((UT(6*N1-4)-UT(6*N2-4)) &
      	-(UTM(6*N1-4)-UTM(6*N2-4)))
	D(6*N1-3)=D(6*N1-3)+(DMP/DT)*((UT(6*N1-3)-UT(6*N2-3)) &
      	-(UTM(6*N1-3)-UTM(6*N2-3)))
	D(6*N1-2)=D(6*N1-2)+(DMP2/DT)*((UT(6*N1-2)-UT(6*N2-2)) &
      	-(UTM(6*N1-2)-UTM(6*N2-2)))
	D(6*N1-1)=D(6*N1-1)+(DMP2/DT)*((UT(6*N1-1)-UT(6*N2-1)) &
      	-(UTM(6*N1-1)-UTM(6*N2-1)))
	D(6*N1-0)=D(6*N1-0)+(DMP2/DT)*((UT(6*N1-0)-UT(6*N2-0)) &
      	-(UTM(6*N1-0)-UTM(6*N2-0)))
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UT(6*N1-5)) &
      	-(UTM(6*N2-5)-UTM(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UT(6*N1-4)) &
      	-(UTM(6*N2-4)-UTM(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UT(6*N1-3)) &
      	-(UTM(6*N2-3)-UTM(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UT(6*N1-2)) &
      	-(UTM(6*N2-2)-UTM(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UT(6*N1-1)) &
      	-(UTM(6*N2-1)-UTM(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UT(6*N1-0)) &
      	-(UTM(6*N2-0)-UTM(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UT(6*N1-5)) &
     	-(UTM(6*N2-5)-UTM(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UT(6*N1-4)) &
     	-(UTM(6*N2-4)-UTM(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UT(6*N1-3)) &
     	-(UTM(6*N2-3)-UTM(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UT(6*N1-2)) &
     	-(UTM(6*N2-2)-UTM(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UT(6*N1-1)) &
     	-(UTM(6*N2-1)-UTM(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UT(6*N1-0)) &
     	-(UTM(6*N2-0)-UTM(6*N1-0)))**2
	ENDIF
	ENDDO

        IF (MOD(myid,ntasks/YN).ne.0) THEN
	DO X=1,NTOL
	IF (EFL(X).NE.0.0) THEN
	N1=NANL(1,X)
	N2=NANL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTL(6*N1-5)) &
     	-(UTM(6*N2-5)-UTML(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTL(6*N1-4)) &
     	-(UTM(6*N2-4)-UTML(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTL(6*N1-3)) &
     	-(UTM(6*N2-3)-UTML(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTL(6*N1-2)) &
     	-(UTM(6*N2-2)-UTML(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTL(6*N1-1)) &
     	-(UTM(6*N2-1)-UTML(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTL(6*N1-0)) &
     	-(UTM(6*N2-0)-UTML(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTL(6*N1-5)) &
     	-(UTM(6*N2-5)-UTML(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTL(6*N1-4)) &
     	-(UTM(6*N2-4)-UTML(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTL(6*N1-3)) &
     	-(UTM(6*N2-3)-UTML(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTL(6*N1-2)) &
     	-(UTM(6*N2-2)-UTML(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTL(6*N1-1)) &
     	-(UTM(6*N2-1)-UTML(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTL(6*N1-0)) &
     	-(UTM(6*N2-0)-UTML(6*N1-0)))**2
	ENDIF
	ENDDO
	ENDIF

        IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
	DO X=1,NTOR
	IF (EFR(X).NE.0.0) THEN
	N1=NANR(1,X)
	N2=NANR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTR(6*N1-5)) &
     	-(UTM(6*N2-5)-UTMR(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTR(6*N1-4)) &
     	-(UTM(6*N2-4)-UTMR(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTR(6*N1-3)) &
     	-(UTM(6*N2-3)-UTMR(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTR(6*N1-2)) &
     	-(UTM(6*N2-2)-UTMR(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTR(6*N1-1)) &
     	-(UTM(6*N2-1)-UTMR(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTR(6*N1-0)) &
     	-(UTM(6*N2-0)-UTMR(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTR(6*N1-5)) &
     	-(UTM(6*N2-5)-UTMR(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTR(6*N1-4)) &
     	-(UTM(6*N2-4)-UTMR(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTR(6*N1-3)) &
     	-(UTM(6*N2-3)-UTMR(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTR(6*N1-2)) &
     	-(UTM(6*N2-2)-UTMR(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTR(6*N1-1)) &
     	-(UTM(6*N2-1)-UTMR(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTR(6*N1-0)) &
     	-(UTM(6*N2-0)-UTMR(6*N1-0)))**2
	ENDIF
	ENDDO
	ENDIF

	IF (myid.lt.(YN-1)*ntasks/YN) THEN
        DO X=1,NTOF
        IF (EFF(X).NE.0.0) THEN
	N1=NANF(1,X)
        N2=NANF(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTF(6*N1-5)) &
        -(UTM(6*N2-5)-UTMF(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTF(6*N1-4)) &
        -(UTM(6*N2-4)-UTMF(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTF(6*N1-3)) &
        -(UTM(6*N2-3)-UTMF(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTF(6*N1-2)) &
        -(UTM(6*N2-2)-UTMF(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTF(6*N1-1)) &
        -(UTM(6*N2-1)-UTMF(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTF(6*N1-0)) &
        -(UTM(6*N2-0)-UTMF(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTF(6*N1-5)) &
        -(UTM(6*N2-5)-UTMF(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTF(6*N1-4)) &
        -(UTM(6*N2-4)-UTMF(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTF(6*N1-3)) &
        -(UTM(6*N2-3)-UTMF(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTF(6*N1-2)) &
        -(UTM(6*N2-2)-UTMF(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTF(6*N1-1)) &
        -(UTM(6*N2-1)-UTMF(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTF(6*N1-0)) &
        -(UTM(6*N2-0)-UTMF(6*N1-0)))**2
	ENDIF
	ENDDO
        ENDIF

	IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,NTOFL
        IF (EFFL(X).NE.0.0) THEN
	N1=NANFL(1,X)
        N2=NANFL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTFL(6*N1-5)) &
        -(UTM(6*N2-5)-UTMFL(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTFL(6*N1-4)) &
        -(UTM(6*N2-4)-UTMFL(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTFL(6*N1-3)) &
        -(UTM(6*N2-3)-UTMFL(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTFL(6*N1-2)) &
        -(UTM(6*N2-2)-UTMFL(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTFL(6*N1-1)) &
        -(UTM(6*N2-1)-UTMFL(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTFL(6*N1-0)) &
        -(UTM(6*N2-0)-UTMFL(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTFL(6*N1-5)) &
        -(UTM(6*N2-5)-UTMFL(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTFL(6*N1-4)) &
        -(UTM(6*N2-4)-UTMFL(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTFL(6*N1-3)) &
        -(UTM(6*N2-3)-UTMFL(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTFL(6*N1-2)) &
        -(UTM(6*N2-2)-UTMFL(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTFL(6*N1-1)) &
        -(UTM(6*N2-1)-UTMFL(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTFL(6*N1-0)) &
        -(UTM(6*N2-0)-UTMFL(6*N1-0)))**2
	ENDIF
	ENDDO
        ENDIF

	IF (myid.lt.(YN-1)*ntasks/YN &
     	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,NTOFR
        IF (EFFR(X).NE.0.0) THEN
	N1=NANFR(1,X)
        N2=NANFR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTFR(6*N1-5)) &
        -(UTM(6*N2-5)-UTMFR(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTFR(6*N1-4)) &
        -(UTM(6*N2-4)-UTMFR(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTFR(6*N1-3)) &
        -(UTM(6*N2-3)-UTMFR(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTFR(6*N1-2)) &
        -(UTM(6*N2-2)-UTMFR(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTFR(6*N1-1)) &
        -(UTM(6*N2-1)-UTMFR(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTFR(6*N1-0)) &
        -(UTM(6*N2-0)-UTMFR(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTFR(6*N1-5)) &
        -(UTM(6*N2-5)-UTMFR(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTFR(6*N1-4)) &
        -(UTM(6*N2-4)-UTMFR(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTFR(6*N1-3)) &
        -(UTM(6*N2-3)-UTMFR(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTFR(6*N1-2)) &
        -(UTM(6*N2-2)-UTMFR(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTFR(6*N1-1)) &
        -(UTM(6*N2-1)-UTMFR(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTFR(6*N1-0)) &
        -(UTM(6*N2-0)-UTMFR(6*N1-0)))**2
	ENDIF
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN) THEN
        DO X=1,NTOB
        IF (EFB(X).NE.0.0) THEN
	N1=NANB(1,X)
        N2=NANB(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTB(6*N1-5)) &
        -(UTM(6*N2-5)-UTMB(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTB(6*N1-4)) &
        -(UTM(6*N2-4)-UTMB(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTB(6*N1-3)) &
        -(UTM(6*N2-3)-UTMB(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTB(6*N1-2)) &
        -(UTM(6*N2-2)-UTMB(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTB(6*N1-1)) &
        -(UTM(6*N2-1)-UTMB(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTB(6*N1-0)) &
        -(UTM(6*N2-0)-UTMB(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTB(6*N1-5)) &
        -(UTM(6*N2-5)-UTMB(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTB(6*N1-4)) &
        -(UTM(6*N2-4)-UTMB(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTB(6*N1-3)) &
        -(UTM(6*N2-3)-UTMB(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTB(6*N1-2)) &
        -(UTM(6*N2-2)-UTMB(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTB(6*N1-1)) &
        -(UTM(6*N2-1)-UTMB(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTB(6*N1-0)) &
        -(UTM(6*N2-0)-UTMB(6*N1-0)))**2
	ENDIF
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,NTOBL
        IF (EFBL(X).NE.0.0) THEN
	N1=NANBL(1,X)
        N2=NANBL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTBL(6*N1-5)) &
        -(UTM(6*N2-5)-UTMBL(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTBL(6*N1-4)) &
        -(UTM(6*N2-4)-UTMBL(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTBL(6*N1-3)) &
        -(UTM(6*N2-3)-UTMBL(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTBL(6*N1-2)) &
        -(UTM(6*N2-2)-UTMBL(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTBL(6*N1-1)) &
        -(UTM(6*N2-1)-UTMBL(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTBL(6*N1-0)) &
        -(UTM(6*N2-0)-UTMBL(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTBL(6*N1-5)) &
        -(UTM(6*N2-5)-UTMBL(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTBL(6*N1-4)) &
        -(UTM(6*N2-4)-UTMBL(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTBL(6*N1-3)) &
        -(UTM(6*N2-3)-UTMBL(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTBL(6*N1-2)) &
        -(UTM(6*N2-2)-UTMBL(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTBL(6*N1-1)) &
        -(UTM(6*N2-1)-UTMBL(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTBL(6*N1-0)) &
        -(UTM(6*N2-0)-UTMBL(6*N1-0)))**2
	ENDIF
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,NTOBR
        IF (EFBR(X).NE.0.0) THEN
	N1=NANBR(1,X)
        N2=NANBR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTBR(6*N1-5)) &
        -(UTM(6*N2-5)-UTMBR(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTBR(6*N1-4)) &
        -(UTM(6*N2-4)-UTMBR(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTBR(6*N1-3)) &
        -(UTM(6*N2-3)-UTMBR(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTBR(6*N1-2)) &
        -(UTM(6*N2-2)-UTMBR(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTBR(6*N1-1)) &
        -(UTM(6*N2-1)-UTMBR(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTBR(6*N1-0)) &
        -(UTM(6*N2-0)-UTMBR(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTBR(6*N1-5)) &
        -(UTM(6*N2-5)-UTMBR(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTBR(6*N1-4)) &
        -(UTM(6*N2-4)-UTMBR(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTBR(6*N1-3)) &
        -(UTM(6*N2-3)-UTMBR(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTBR(6*N1-2)) &
        -(UTM(6*N2-2)-UTMBR(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTBR(6*N1-1)) &
        -(UTM(6*N2-1)-UTMBR(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTBR(6*N1-0)) &
        -(UTM(6*N2-0)-UTMBR(6*N1-0)))**2
	ENDIF
	ENDDO
        ENDIF

!-----------------------------------------------------------------
	DO X=1,FXC
	N1=FXF(1,X)
	N2=FXF(2,X)
	D(6*N1-5)=D(6*N1-5)+(DMP/DT)*((UT(6*N1-5)-UT(6*N2-5)) &
     	-(UTM(6*N1-5)-UTM(6*N2-5)))
	D(6*N1-4)=D(6*N1-4)+(DMP/DT)*((UT(6*N1-4)-UT(6*N2-4)) &
     	-(UTM(6*N1-4)-UTM(6*N2-4)))
	D(6*N1-3)=D(6*N1-3)+(DMP/DT)*((UT(6*N1-3)-UT(6*N2-3)) &
     	-(UTM(6*N1-3)-UTM(6*N2-3)))
	D(6*N1-2)=D(6*N1-2)+(DMP2/DT)*((UT(6*N1-2)-UT(6*N2-2)) &
     	-(UTM(6*N1-2)-UTM(6*N2-2)))
	D(6*N1-1)=D(6*N1-1)+(DMP2/DT)*((UT(6*N1-1)-UT(6*N2-1)) &
     	-(UTM(6*N1-1)-UTM(6*N2-1)))
	D(6*N1-0)=D(6*N1-0)+(DMP2/DT)*((UT(6*N1-0)-UT(6*N2-0)) &
     	-(UTM(6*N1-0)-UTM(6*N2-0)))
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UT(6*N1-5)) &
     	-(UTM(6*N2-5)-UTM(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UT(6*N1-4)) &
     	-(UTM(6*N2-4)-UTM(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UT(6*N1-3)) &
     	-(UTM(6*N2-3)-UTM(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UT(6*N1-2)) &
     	-(UTM(6*N2-2)-UTM(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UT(6*N1-1)) &
     	-(UTM(6*N2-1)-UTM(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UT(6*N1-0)) &
     	-(UTM(6*N2-0)-UTM(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UT(6*N1-5)) &
     	-(UTM(6*N2-5)-UTM(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UT(6*N1-4)) &
     	-(UTM(6*N2-4)-UTM(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UT(6*N1-3)) &
     	-(UTM(6*N2-3)-UTM(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UT(6*N1-2)) &
     	-(UTM(6*N2-2)-UTM(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UT(6*N1-1)) &
     	-(UTM(6*N2-1)-UTM(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UT(6*N1-0)) &
     	-(UTM(6*N2-0)-UTM(6*N1-0)))**2
	ENDDO

        IF (MOD(myid,ntasks/YN).ne.0) THEN
	DO X=1,FXL
	N1=FXFL(1,X)
	N2=FXFL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTL(6*N1-5)) &
      	-(UTM(6*N2-5)-UTML(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTL(6*N1-4)) &
      	-(UTM(6*N2-4)-UTML(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTL(6*N1-3)) &
      	-(UTM(6*N2-3)-UTML(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTL(6*N1-2)) &
      	-(UTM(6*N2-2)-UTML(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTL(6*N1-1)) &
      	-(UTM(6*N2-1)-UTML(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTL(6*N1-0)) &
      	-(UTM(6*N2-0)-UTML(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTL(6*N1-5)) &
      	-(UTM(6*N2-5)-UTML(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTL(6*N1-4)) &
      	-(UTM(6*N2-4)-UTML(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTL(6*N1-3)) &
      	-(UTM(6*N2-3)-UTML(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTL(6*N1-2)) &
      	-(UTM(6*N2-2)-UTML(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTL(6*N1-1)) &
      	-(UTM(6*N2-1)-UTML(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTL(6*N1-0)) &
      	-(UTM(6*N2-0)-UTML(6*N1-0)))**2
	ENDDO
	ENDIF

        IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
	DO X=1,FXR
	N1=FXFR(1,X)
	N2=FXFR(2,X)
 	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTR(6*N1-5)) &
      	-(UTM(6*N2-5)-UTMR(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTR(6*N1-4)) &
      	-(UTM(6*N2-4)-UTMR(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTR(6*N1-3)) &
      	-(UTM(6*N2-3)-UTMR(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTR(6*N1-2)) &
      	-(UTM(6*N2-2)-UTMR(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTR(6*N1-1)) &
      	-(UTM(6*N2-1)-UTMR(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTR(6*N1-0)) &
      	-(UTM(6*N2-0)-UTMR(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTR(6*N1-5)) &
      	-(UTM(6*N2-5)-UTMR(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTR(6*N1-4)) &
      	-(UTM(6*N2-4)-UTMR(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTR(6*N1-3)) &
      	-(UTM(6*N2-3)-UTMR(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTR(6*N1-2)) &
      	-(UTM(6*N2-2)-UTMR(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTR(6*N1-1)) &
      	-(UTM(6*N2-1)-UTMR(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTR(6*N1-0)) &
      	-(UTM(6*N2-0)-UTMR(6*N1-0)))**2
	ENDDO
	ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN) THEN
        DO X=1,FXCF
        N1=FXFF(1,X)
        N2=FXFF(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTF(6*N1-5)) &
        -(UTM(6*N2-5)-UTMF(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTF(6*N1-4)) &
        -(UTM(6*N2-4)-UTMF(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTF(6*N1-3)) &
        -(UTM(6*N2-3)-UTMF(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTF(6*N1-2)) &
        -(UTM(6*N2-2)-UTMF(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTF(6*N1-1)) &
        -(UTM(6*N2-1)-UTMF(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTF(6*N1-0)) &
        -(UTM(6*N2-0)-UTMF(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTF(6*N1-5)) &
        -(UTM(6*N2-5)-UTMF(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTF(6*N1-4)) &
        -(UTM(6*N2-4)-UTMF(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTF(6*N1-3)) &
        -(UTM(6*N2-3)-UTMF(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTF(6*N1-2)) &
        -(UTM(6*N2-2)-UTMF(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTF(6*N1-1)) &
        -(UTM(6*N2-1)-UTMF(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTF(6*N1-0)) &
        -(UTM(6*N2-0)-UTMF(6*N1-0)))**2
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,FXCFL
        N1=FXFFL(1,X)
        N2=FXFFL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTFL(6*N1-5)) &
        -(UTM(6*N2-5)-UTMFL(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTFL(6*N1-4)) &
        -(UTM(6*N2-4)-UTMFL(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTFL(6*N1-3)) &
        -(UTM(6*N2-3)-UTMFL(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTFL(6*N1-2)) &
        -(UTM(6*N2-2)-UTMFL(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTFL(6*N1-1)) &
        -(UTM(6*N2-1)-UTMFL(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTFL(6*N1-0)) &
        -(UTM(6*N2-0)-UTMFL(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTFL(6*N1-5)) &
        -(UTM(6*N2-5)-UTMFL(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTFL(6*N1-4)) &
        -(UTM(6*N2-4)-UTMFL(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTFL(6*N1-3)) &
        -(UTM(6*N2-3)-UTMFL(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTFL(6*N1-2)) &
        -(UTM(6*N2-2)-UTMFL(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTFL(6*N1-1)) &
        -(UTM(6*N2-1)-UTMFL(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTFL(6*N1-0)) &
        -(UTM(6*N2-0)-UTMFL(6*N1-0)))**2
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN &
      	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,FXCFR
        N1=FXFFR(1,X)
        N2=FXFFR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTFR(6*N1-5)) &
        -(UTM(6*N2-5)-UTMFR(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTFR(6*N1-4)) &
        -(UTM(6*N2-4)-UTMFR(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTFR(6*N1-3)) &
        -(UTM(6*N2-3)-UTMFR(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTFR(6*N1-2)) &
        -(UTM(6*N2-2)-UTMFR(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTFR(6*N1-1)) &
        -(UTM(6*N2-1)-UTMFR(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTFR(6*N1-0)) &
        -(UTM(6*N2-0)-UTMFR(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTFR(6*N1-5)) &
        -(UTM(6*N2-5)-UTMFR(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTFR(6*N1-4)) &
        -(UTM(6*N2-4)-UTMFR(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTFR(6*N1-3)) &
        -(UTM(6*N2-3)-UTMFR(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTFR(6*N1-2)) &
        -(UTM(6*N2-2)-UTMFR(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTFR(6*N1-1)) &
        -(UTM(6*N2-1)-UTMFR(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTFR(6*N1-0)) &
        -(UTM(6*N2-0)-UTMFR(6*N1-0)))**2
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN) THEN
        DO X=1,FXCB
	N1=FXFB(1,X)
        N2=FXFB(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTB(6*N1-5)) &
        -(UTM(6*N2-5)-UTMB(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTB(6*N1-4)) &
        -(UTM(6*N2-4)-UTMB(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTB(6*N1-3)) &
        -(UTM(6*N2-3)-UTMB(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTB(6*N1-2)) &
        -(UTM(6*N2-2)-UTMB(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTB(6*N1-1)) &
        -(UTM(6*N2-1)-UTMB(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTB(6*N1-0)) &
        -(UTM(6*N2-0)-UTMB(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTB(6*N1-5)) &
        -(UTM(6*N2-5)-UTMB(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTB(6*N1-4)) &
        -(UTM(6*N2-4)-UTMB(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTB(6*N1-3)) &
        -(UTM(6*N2-3)-UTMB(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTB(6*N1-2)) &
        -(UTM(6*N2-2)-UTMB(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTB(6*N1-1)) &
        -(UTM(6*N2-1)-UTMB(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTB(6*N1-0)) &
        -(UTM(6*N2-0)-UTMB(6*N1-0)))**2
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,FXCBL
	N1=FXFBL(1,X)
        N2=FXFBL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTBL(6*N1-5)) &
        -(UTM(6*N2-5)-UTMBL(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTBL(6*N1-4)) &
        -(UTM(6*N2-4)-UTMBL(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTBL(6*N1-3)) &
        -(UTM(6*N2-3)-UTMBL(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTBL(6*N1-2)) &
        -(UTM(6*N2-2)-UTMBL(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTBL(6*N1-1)) &
        -(UTM(6*N2-1)-UTMBL(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTBL(6*N1-0)) &
        -(UTM(6*N2-0)-UTMBL(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTBL(6*N1-5)) &
        -(UTM(6*N2-5)-UTMBL(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTBL(6*N1-4)) &
        -(UTM(6*N2-4)-UTMBL(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTBL(6*N1-3)) &
        -(UTM(6*N2-3)-UTMBL(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTBL(6*N1-2)) &
        -(UTM(6*N2-2)-UTMBL(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTBL(6*N1-1)) &
        -(UTM(6*N2-1)-UTMBL(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTBL(6*N1-0)) &
        -(UTM(6*N2-0)-UTMBL(6*N1-0)))**2
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN &
     	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,FXCBR
	N1=FXFBR(1,X)
        N2=FXFBR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTBR(6*N1-5)) &
        -(UTM(6*N2-5)-UTMBR(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTBR(6*N1-4)) &
        -(UTM(6*N2-4)-UTMBR(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTBR(6*N1-3)) &
        -(UTM(6*N2-3)-UTMBR(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTBR(6*N1-2)) &
        -(UTM(6*N2-2)-UTMBR(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTBR(6*N1-1)) &
        -(UTM(6*N2-1)-UTMBR(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTBR(6*N1-0)) &
        -(UTM(6*N2-0)-UTMBR(6*N1-0)))
	DPE=DPE+(DMP/DT)*((UT(6*N2-5)-UTBR(6*N1-5)) &
        -(UTM(6*N2-5)-UTMBR(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-4)-UTBR(6*N1-4)) &
        -(UTM(6*N2-4)-UTMBR(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT(6*N2-3)-UTBR(6*N1-3)) &
        -(UTM(6*N2-3)-UTMBR(6*N1-3)))**2
	DPE=DPE+(DMP2/DT)*((UT(6*N2-2)-UTBR(6*N1-2)) &
        -(UTM(6*N2-2)-UTMBR(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-1)-UTBR(6*N1-1)) &
        -(UTM(6*N2-1)-UTMBR(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT(6*N2-0)-UTBR(6*N1-0)) &
        -(UTM(6*N2-0)-UTMBR(6*N1-0)))**2
	ENDDO
        ENDIF

!---------------------------------------------------------

        DO X=1,NN
        F(6*X-5)= -(VDP(X)/(DT*2))*UTM(6*X-5)
        F(6*X-4)= -(VDP(X)/(DT*2))*UTM(6*X-4)
        F(6*X-3)= -(VDP(X)/(DT*2))*UTM(6*X-3)
        F(6*X-2)= -(VDP(X)/(DT*2))*UTM(6*X-2)
        F(6*X-1)= -(VDP(X)/(DT*2))*UTM(6*X-1)
        F(6*X-0)= -(VDP(X)/(DT*2))*UTM(6*X-0)
	END DO

	DO X=1,6*NN
	R(X)=-A(X)-C(X)-D(X)-F(X)
	EN(X)=EN(X)+A(X)*(UT(X)-UTM(X))
	END DO

	RETURN
END SUBROUTINE EFFLOAD
