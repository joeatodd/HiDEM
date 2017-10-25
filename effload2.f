	SUBROUTINE EFFLOAD2(S,NTOT,NN,T,DT,M,JS,DMP,DMP2,UT,UTM,D,EN,RY,
     1	FXF,FXC,VDP,DPE,EF,EFL,EFR,NAN,NRXF,MFIL,CT,NANL,NANR,NTOL,
     1	NTOR,NRXFL,NRXFR,UTL,UTR,myid,ntasks,FXL,FXFL,FXR,FXFR,L,
     1  NTOF,NTOB,NRXFF,NRXFB,UTF,UTB,EFF,EFB,
     1  FXCF,FXFF,FXCB,FXFB,NANB,NANF,PNN,YN,
     1	NANFL,NANFR,NANBL,NANBR,
     1	NRXFFL,NRXFFR,NRXFBL,NRXFBR,
     1	UTFL,UTFR,UTBL,UTBR,
     1	EFFL,EFFR,EFBL,EFBR,
     1  NTOFL,NTOFR,NTOBL,NTOBR,
     1  FXCFL,FXFFL,FXCBL,FXFBL,
     1	FXCFR,FXFFR,FXCBR,FXFBR)
	
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
	REAL*8 D(NODM)
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
	D(I)=0.0
	END DO


      dest=myid+1
      source=myid-1
      tag=133
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)
     1CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (MOD(myid,ntasks/YN).ne.0)
     1CALL MPI_Recv(UTML,6*PNN(source),MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid-1
      source=myid+1
      tag=135
      IF (MOD(myid,ntasks/YN).ne.0)
     1CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)
     1CALL MPI_Recv(UTMR,6*PNN(source),MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid+ntasks/YN
      source=myid-ntasks/YN
      tag=137
      IF (myid.lt.(YN-1)*ntasks/YN)
     1CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.ge.ntasks/YN)
     1CALL MPI_Recv(UTMB,6*PNN(source),MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid+ntasks/YN+1
      source=myid-ntasks/YN-1
      tag=139
      IF (myid.lt.(YN-1)*ntasks/YN
     1.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)
     1CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)
     1CALL MPI_Recv(UTMBL,6*PNN(source),MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid+ntasks/YN-1
      source=myid-ntasks/YN+1
      tag=141
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)
     1CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)
     1CALL MPI_Recv(UTMBR,6*PNN(source),MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid-ntasks/YN
      source=myid+ntasks/YN
      tag=143
      IF (myid.ge.ntasks/YN)
     1CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN)
     1CALL MPI_Recv(UTMF,6*PNN(source),MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid-ntasks/YN+1
      source=myid+ntasks/YN-1
      tag=145
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)
     1CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)
     1CALL MPI_Recv(UTMFL,6*PNN(source),MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid-ntasks/YN-1
      source=myid+ntasks/YN+1
      tag=147
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)
     1CALL MPI_Send(UTM,6*NN,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN
     1.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)
     1CALL MPI_Recv(UTMFR,6*PNN(source),MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)

	DO X=1,NTOT
	IF (EF(X).NE.0.0) THEN
	N1=NAN(1,X)
	N2=NAN(2,X)
	D(6*N1-5)=D(6*N1-5)+(DMP/DT)*((UT(6*N1-5)-UT(6*N2-5))
     1	-(UTM(6*N1-5)-UTM(6*N2-5)))
	D(6*N1-4)=D(6*N1-4)+(DMP/DT)*((UT(6*N1-4)-UT(6*N2-4))
     1	-(UTM(6*N1-4)-UTM(6*N2-4)))
	D(6*N1-3)=D(6*N1-3)+(DMP/DT)*((UT(6*N1-3)-UT(6*N2-3))
     1	-(UTM(6*N1-3)-UTM(6*N2-3)))
	D(6*N1-2)=D(6*N1-2)+(DMP2/DT)*((UT(6*N1-2)-UT(6*N2-2))
     1	-(UTM(6*N1-2)-UTM(6*N2-2)))
	D(6*N1-1)=D(6*N1-1)+(DMP2/DT)*((UT(6*N1-1)-UT(6*N2-1))
     1	-(UTM(6*N1-1)-UTM(6*N2-1)))
	D(6*N1-0)=D(6*N1-0)+(DMP2/DT)*((UT(6*N1-0)-UT(6*N2-0))
     1	-(UTM(6*N1-0)-UTM(6*N2-0)))
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UT(6*N1-5))
     1	-(UTM(6*N2-5)-UTM(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UT(6*N1-4))
     1	-(UTM(6*N2-4)-UTM(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UT(6*N1-3))
     1	-(UTM(6*N2-3)-UTM(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UT(6*N1-2))
     1	-(UTM(6*N2-2)-UTM(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UT(6*N1-1))
     1	-(UTM(6*N2-1)-UTM(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UT(6*N1-0))
     1	-(UTM(6*N2-0)-UTM(6*N1-0)))
	ENDIF
	ENDDO

        IF (MOD(myid,ntasks/YN).ne.0) THEN
	DO X=1,NTOL
	IF (EFL(X).NE.0.0) THEN
	N1=NANL(1,X)
	N2=NANL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTL(6*N1-5))
     1	-(UTM(6*N2-5)-UTML(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTL(6*N1-4))
     1	-(UTM(6*N2-4)-UTML(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTL(6*N1-3))
     1	-(UTM(6*N2-3)-UTML(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTL(6*N1-2))
     1	-(UTM(6*N2-2)-UTML(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTL(6*N1-1))
     1	-(UTM(6*N2-1)-UTML(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTL(6*N1-0))
     1	-(UTM(6*N2-0)-UTML(6*N1-0)))
	ENDIF
	ENDDO
	ENDIF

        IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
	DO X=1,NTOR
	IF (EFR(X).NE.0.0) THEN
	N1=NANR(1,X)
	N2=NANR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTR(6*N1-5))
     1	-(UTM(6*N2-5)-UTMR(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTR(6*N1-4))
     1	-(UTM(6*N2-4)-UTMR(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTR(6*N1-3))
     1	-(UTM(6*N2-3)-UTMR(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTR(6*N1-2))
     1	-(UTM(6*N2-2)-UTMR(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTR(6*N1-1))
     1	-(UTM(6*N2-1)-UTMR(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTR(6*N1-0))
     1	-(UTM(6*N2-0)-UTMR(6*N1-0)))
	ENDIF
	ENDDO
	ENDIF

	IF (myid.lt.(YN-1)*ntasks/YN) THEN
        DO X=1,NTOF
        IF (EFF(X).NE.0.0) THEN
	N1=NANF(1,X)
        N2=NANF(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTF(6*N1-5))
     1  -(UTM(6*N2-5)-UTMF(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTF(6*N1-4))
     1  -(UTM(6*N2-4)-UTMF(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTF(6*N1-3))
     1  -(UTM(6*N2-3)-UTMF(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTF(6*N1-2))
     1  -(UTM(6*N2-2)-UTMF(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTF(6*N1-1))
     1  -(UTM(6*N2-1)-UTMF(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTF(6*N1-0))
     1  -(UTM(6*N2-0)-UTMF(6*N1-0)))
	ENDIF
	ENDDO
        ENDIF

	IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,NTOFL
        IF (EFFL(X).NE.0.0) THEN
	N1=NANFL(1,X)
        N2=NANFL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTFL(6*N1-5))
     1  -(UTM(6*N2-5)-UTMFL(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTFL(6*N1-4))
     1  -(UTM(6*N2-4)-UTMFL(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTFL(6*N1-3))
     1  -(UTM(6*N2-3)-UTMFL(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTFL(6*N1-2))
     1  -(UTM(6*N2-2)-UTMFL(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTFL(6*N1-1))
     1  -(UTM(6*N2-1)-UTMFL(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTFL(6*N1-0))
     1  -(UTM(6*N2-0)-UTMFL(6*N1-0)))
	ENDIF
	ENDDO
        ENDIF

	IF (myid.lt.(YN-1)*ntasks/YN
     1	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,NTOFR
        IF (EFFR(X).NE.0.0) THEN
	N1=NANFR(1,X)
        N2=NANFR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTFR(6*N1-5))
     1  -(UTM(6*N2-5)-UTMFR(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTFR(6*N1-4))
     1  -(UTM(6*N2-4)-UTMFR(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTFR(6*N1-3))
     1  -(UTM(6*N2-3)-UTMFR(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTFR(6*N1-2))
     1  -(UTM(6*N2-2)-UTMFR(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTFR(6*N1-1))
     1  -(UTM(6*N2-1)-UTMFR(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTFR(6*N1-0))
     1  -(UTM(6*N2-0)-UTMFR(6*N1-0)))
	ENDIF
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN) THEN
        DO X=1,NTOB
        IF (EFB(X).NE.0.0) THEN
	N1=NANB(1,X)
        N2=NANB(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTB(6*N1-5))
     1  -(UTM(6*N2-5)-UTMB(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTB(6*N1-4))
     1  -(UTM(6*N2-4)-UTMB(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTB(6*N1-3))
     1  -(UTM(6*N2-3)-UTMB(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTB(6*N1-2))
     1  -(UTM(6*N2-2)-UTMB(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTB(6*N1-1))
     1  -(UTM(6*N2-1)-UTMB(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTB(6*N1-0))
     1  -(UTM(6*N2-0)-UTMB(6*N1-0)))
	ENDIF
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,NTOBL
        IF (EFBL(X).NE.0.0) THEN
	N1=NANBL(1,X)
        N2=NANBL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTBL(6*N1-5))
     1  -(UTM(6*N2-5)-UTMBL(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTBL(6*N1-4))
     1  -(UTM(6*N2-4)-UTMBL(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTBL(6*N1-3))
     1  -(UTM(6*N2-3)-UTMBL(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTBL(6*N1-2))
     1  -(UTM(6*N2-2)-UTMBL(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTBL(6*N1-1))
     1  -(UTM(6*N2-1)-UTMBL(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTBL(6*N1-0))
     1  -(UTM(6*N2-0)-UTMBL(6*N1-0)))
	ENDIF
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,NTOBR
        IF (EFBR(X).NE.0.0) THEN
	N1=NANBR(1,X)
        N2=NANBR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTBR(6*N1-5))
     1  -(UTM(6*N2-5)-UTMBR(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTBR(6*N1-4))
     1  -(UTM(6*N2-4)-UTMBR(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTBR(6*N1-3))
     1  -(UTM(6*N2-3)-UTMBR(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTBR(6*N1-2))
     1  -(UTM(6*N2-2)-UTMBR(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTBR(6*N1-1))
     1  -(UTM(6*N2-1)-UTMBR(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTBR(6*N1-0))
     1  -(UTM(6*N2-0)-UTMBR(6*N1-0)))
	ENDIF
	ENDDO
        ENDIF

c-----------------------------------------------------------------
	DO X=1,FXC
	N1=FXF(1,X)
	N2=FXF(2,X)
	D(6*N1-5)=D(6*N1-5)+(DMP/DT)*((UT(6*N1-5)-UT(6*N2-5))
     1	-(UTM(6*N1-5)-UTM(6*N2-5)))
	D(6*N1-4)=D(6*N1-4)+(DMP/DT)*((UT(6*N1-4)-UT(6*N2-4))
     1	-(UTM(6*N1-4)-UTM(6*N2-4)))
	D(6*N1-3)=D(6*N1-3)+(DMP/DT)*((UT(6*N1-3)-UT(6*N2-3))
     1	-(UTM(6*N1-3)-UTM(6*N2-3)))
	D(6*N1-2)=D(6*N1-2)+(DMP2/DT)*((UT(6*N1-2)-UT(6*N2-2))
     1	-(UTM(6*N1-2)-UTM(6*N2-2)))
	D(6*N1-1)=D(6*N1-1)+(DMP2/DT)*((UT(6*N1-1)-UT(6*N2-1))
     1	-(UTM(6*N1-1)-UTM(6*N2-1)))
	D(6*N1-0)=D(6*N1-0)+(DMP2/DT)*((UT(6*N1-0)-UT(6*N2-0))
     1	-(UTM(6*N1-0)-UTM(6*N2-0)))
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UT(6*N1-5))
     1	-(UTM(6*N2-5)-UTM(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UT(6*N1-4))
     1	-(UTM(6*N2-4)-UTM(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UT(6*N1-3))
     1	-(UTM(6*N2-3)-UTM(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UT(6*N1-2))
     1	-(UTM(6*N2-2)-UTM(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UT(6*N1-1))
     1	-(UTM(6*N2-1)-UTM(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UT(6*N1-0))
     1	-(UTM(6*N2-0)-UTM(6*N1-0)))
	ENDDO

        IF (MOD(myid,ntasks/YN).ne.0) THEN
	DO X=1,FXL
	N1=FXFL(1,X)
	N2=FXFL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTL(6*N1-5))
     1	-(UTM(6*N2-5)-UTML(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTL(6*N1-4))
     1	-(UTM(6*N2-4)-UTML(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTL(6*N1-3))
     1	-(UTM(6*N2-3)-UTML(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTL(6*N1-2))
     1	-(UTM(6*N2-2)-UTML(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTL(6*N1-1))
     1	-(UTM(6*N2-1)-UTML(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTL(6*N1-0))
     1	-(UTM(6*N2-0)-UTML(6*N1-0)))
	ENDDO
	ENDIF

        IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
	DO X=1,FXR
	N1=FXFR(1,X)
	N2=FXFR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTR(6*N1-5))
     1	-(UTM(6*N2-5)-UTMR(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTR(6*N1-4))
     1	-(UTM(6*N2-4)-UTMR(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTR(6*N1-3))
     1	-(UTM(6*N2-3)-UTMR(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTR(6*N1-2))
     1	-(UTM(6*N2-2)-UTMR(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTR(6*N1-1))
     1	-(UTM(6*N2-1)-UTMR(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTR(6*N1-0))
     1	-(UTM(6*N2-0)-UTMR(6*N1-0)))
	ENDDO
	ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN) THEN
        DO X=1,FXCF
        N1=FXFF(1,X)
        N2=FXFF(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTF(6*N1-5))
     1  -(UTM(6*N2-5)-UTMF(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTF(6*N1-4))
     1  -(UTM(6*N2-4)-UTMF(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTF(6*N1-3))
     1  -(UTM(6*N2-3)-UTMF(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTF(6*N1-2))
     1  -(UTM(6*N2-2)-UTMF(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTF(6*N1-1))
     1  -(UTM(6*N2-1)-UTMF(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTF(6*N1-0))
     1  -(UTM(6*N2-0)-UTMF(6*N1-0)))
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,FXCFL
        N1=FXFFL(1,X)
        N2=FXFFL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTFL(6*N1-5))
     1  -(UTM(6*N2-5)-UTMFL(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTFL(6*N1-4))
     1  -(UTM(6*N2-4)-UTMFL(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTFL(6*N1-3))
     1  -(UTM(6*N2-3)-UTMFL(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTFL(6*N1-2))
     1  -(UTM(6*N2-2)-UTMFL(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTFL(6*N1-1))
     1  -(UTM(6*N2-1)-UTMFL(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTFL(6*N1-0))
     1  -(UTM(6*N2-0)-UTMFL(6*N1-0)))
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN
     1	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,FXCFR
        N1=FXFFR(1,X)
        N2=FXFFR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTFR(6*N1-5))
     1  -(UTM(6*N2-5)-UTMFR(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTFR(6*N1-4))
     1  -(UTM(6*N2-4)-UTMFR(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTFR(6*N1-3))
     1  -(UTM(6*N2-3)-UTMFR(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTFR(6*N1-2))
     1  -(UTM(6*N2-2)-UTMFR(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTFR(6*N1-1))
     1  -(UTM(6*N2-1)-UTMFR(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTFR(6*N1-0))
     1  -(UTM(6*N2-0)-UTMFR(6*N1-0)))
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN) THEN
        DO X=1,FXCB
	N1=FXFB(1,X)
        N2=FXFB(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTB(6*N1-5))
     1  -(UTM(6*N2-5)-UTMB(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTB(6*N1-4))
     1  -(UTM(6*N2-4)-UTMB(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTB(6*N1-3))
     1  -(UTM(6*N2-3)-UTMB(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTB(6*N1-2))
     1  -(UTM(6*N2-2)-UTMB(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTB(6*N1-1))
     1  -(UTM(6*N2-1)-UTMB(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTB(6*N1-0))
     1  -(UTM(6*N2-0)-UTMB(6*N1-0)))
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO X=1,FXCBL
	N1=FXFBL(1,X)
        N2=FXFBL(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTBL(6*N1-5))
     1  -(UTM(6*N2-5)-UTMBL(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTBL(6*N1-4))
     1  -(UTM(6*N2-4)-UTMBL(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTBL(6*N1-3))
     1  -(UTM(6*N2-3)-UTMBL(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTBL(6*N1-2))
     1  -(UTM(6*N2-2)-UTMBL(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTBL(6*N1-1))
     1  -(UTM(6*N2-1)-UTMBL(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTBL(6*N1-0))
     1  -(UTM(6*N2-0)-UTMBL(6*N1-0)))
	ENDDO
        ENDIF

	IF (myid.ge.ntasks/YN
     1	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO X=1,FXCBR
	N1=FXFBR(1,X)
        N2=FXFBR(2,X)
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT(6*N2-5)-UTBR(6*N1-5))
     1  -(UTM(6*N2-5)-UTMBR(6*N1-5)))
        D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT(6*N2-4)-UTBR(6*N1-4))
     1  -(UTM(6*N2-4)-UTMBR(6*N1-4)))
        D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT(6*N2-3)-UTBR(6*N1-3))
     1  -(UTM(6*N2-3)-UTMBR(6*N1-3)))
        D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT(6*N2-2)-UTBR(6*N1-2))
     1  -(UTM(6*N2-2)-UTMBR(6*N1-2)))
        D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT(6*N2-1)-UTBR(6*N1-1))
     1  -(UTM(6*N2-1)-UTMBR(6*N1-1)))
        D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT(6*N2-0)-UTBR(6*N1-0))
     1  -(UTM(6*N2-0)-UTMBR(6*N1-0)))
	ENDDO
        ENDIF


	RETURN
	END





