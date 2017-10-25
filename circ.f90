	SUBROUTINE CIRC(ND,NN,NRXF,UT,FRX,FRY,FRZ, &
      T,IS,DT,WE,EFC,FXF,FXC,NDR,NCL,NDL, &
      NDLL,NDLR,NRXFL,NRXFR,UTL,UTR,myid,ntasks,FXFL,FXFR,FXL,FXR, &
      FXCF,FXFF,UTF,NRXFF,NDLF,NDF, &
      FXCB,FXFB,UTB,NRXFB,NDLB,NDB,LNN,YN, &
      NDFL,NDFR,NDBL,NDBR, &
      FXCFR,FXCFL,FXCBR,FXCBL, &
      NDLFL,NDLFR,NDLBL,NDLBR, &
      UTFL,UTFR,UTBL,UTBR, &
      NRXFBL,NRXFBR,NRXFFL,NRXFFR, &
      FXFFR,FXFFL,FXFBR,FXFBL,SCL)

	IMPLICIT NONE
        include 'mpif.h'
        include 'param.dat'
	REAL*8 X1,X2,Y1,Y2,Z1,Z2
	REAL*8 T1,T2
	REAL*8 NRXF(3,NOMA),NRXFR(3,NOMA),NRXFL(3,NOMA)
	REAL*8 NRXFF(3,NOMA),NRXFB(3,NOMA),EFC(NOMA)
	REAL*8 NRXFFL(3,NOMA),NRXFBL(3,NOMA)
	REAL*8 NRXFFR(3,NOMA),NRXFBR(3,NOMA)
	REAL*8 UT(NODM),UTR(NODM),UTL(NODM),UTB(NODM),UTF(NODM)
	REAL*8 UTBL(NODM),UTFL(NODM),UTBR(NODM),UTFR(NODM)
	REAL*8 SX,SY,SZ,SUM,T,WE(NOMA),L0
	REAL*8 DDEL,DWE,OWE,DT,ESUM,LNN
	REAL*8 LS,LS2,DEL,SCL
	INTEGER myid,ntasks,ierr,YN
        INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),comm
	INTEGER FXF(2,NODC),NDL(2,NODC),NDLR(2,NODC),NDLL(2,NODC)
	INTEGER NDLF(2,NODC),NDLB(2,NODC)
	INTEGER NDLFL(2,NODC),NDLBL(2,NODC)
	INTEGER NDLFR(2,NODC),NDLBR(2,NODC)
	REAL*8 RC,RCX,RCY,RCZ,FRX(NOMA),FRY(NOMA),FRZ(NOMA)
	INTEGER NTOT,I,N1,N2,ND,IS,NN,FXC,NCL,NDR,NDF,NDB
	INTEGER NDFL,NDFR,NDBL,NDBR
	INTEGER FXL,FXR,FXFR(2,NODC),FXFL(2,NODC)
	INTEGER FXCF,FXCB,FXFB(2,NODC),FXFF(2,NODC)
	INTEGER FXCFL,FXCBL,FXFBL(2,NODC),FXFFL(2,NODC)
	INTEGER FXCFR,FXCBR,FXFBR(2,NODC),FXFFR(2,NODC)


	DO I=1,NOMA
	FRX(I)=0.0
	FRY(I)=0.0
	FRZ(I)=0.0
	WE(I)=0.0
	END DO

	FXC=0
	DO I=1,ND
	N1=NDL(1,I)
	N2=NDL(2,I)
	X1=NRXF(1,N1)+UT(6*N1-5)
	Y1=NRXF(2,N1)+UT(6*N1-4)
	Z1=NRXF(3,N1)+UT(6*N1-3)
	X2=NRXF(1,N2)+UT(6*N2-5)
	Y2=NRXF(2,N2)+UT(6*N2-4)
	Z2=NRXF(3,N2)+UT(6*N2-3)
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
	 FXC=FXC+1
	 FXF(1,FXC)=N1 
	 FXF(2,FXC)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.000*SCL) THEN
         FRX(N1)=FRX(N1)+SCL**2.0*1.0e+06*(LNN-RC)*RCX
         FRY(N1)=FRY(N1)+SCL**2.0*1.0e+06*(LNN-RC)*RCY
         FRZ(N1)=FRZ(N1)+SCL**2.0*1.0e+06*(LNN-RC)*RCZ
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCZ
 	 WE(N2)=WE(N2)+SCL**2.0*0.5e+06*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO

c        IF (myid.ne.0.and.myid.ne.(ntasks-1)/2+1) THEN
        IF (MOD(myid,ntasks/YN).ne.0) THEN
	FXL=0
	DO I=1,NCL
	N1=NDLL(1,I)
	N2=NDLL(2,I)
	X1=NRXFL(1,N1)+UTL(6*N1-5)
	Y1=NRXFL(2,N1)+UTL(6*N1-4)
	Z1=NRXFL(3,N1)+UTL(6*N1-3)
	X2=NRXF(1,N2)+UT(6*N2-5)
	Y2=NRXF(2,N2)+UT(6*N2-4)
	Z2=NRXF(3,N2)+UT(6*N2-3)
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
	 FXL=FXL+1
	 FXFL(1,FXL)=N1 
	 FXFL(2,FXL)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.000*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+06*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

c        IF (myid.ne.ntasks-1.and.myid.ne.(ntasks-1)/2) THEN
        IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
	FXR=0
	DO I=1,NDR
	N1=NDLR(1,I)
	N2=NDLR(2,I)
	X1=NRXFR(1,N1)+UTR(6*N1-5)
	Y1=NRXFR(2,N1)+UTR(6*N1-4)
	Z1=NRXFR(3,N1)+UTR(6*N1-3)
	X2=NRXF(1,N2)+UT(6*N2-5)
	Y2=NRXF(2,N2)+UT(6*N2-4)
	Z2=NRXF(3,N2)+UT(6*N2-3)
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
	 FXR=FXR+1
	 FXFR(1,FXR)=N1 
	 FXFR(2,FXR)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.000*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+06*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

c------------------------------------------------

        IF (myid.lt.(YN-1)*ntasks/YN) THEN
	FXCF=0
	DO I=1,NDF
	N1=NDLF(1,I)
	N2=NDLF(2,I)
	X1=NRXFF(1,N1)+UTF(6*N1-5)
	Y1=NRXFF(2,N1)+UTF(6*N1-4)
	Z1=NRXFF(3,N1)+UTF(6*N1-3)
	X2=NRXF(1,N2)+UT(6*N2-5)
	Y2=NRXF(2,N2)+UT(6*N2-4)
	Z2=NRXF(3,N2)+UT(6*N2-3)
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
	 FXCF=FXCF+1
	 FXFF(1,FXCF)=N1 
	 FXFF(2,FXCF)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.000*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCZ
 	 WE(N2)=WE(N2)+SCL**2.0*0.5e+06*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
	FXCFL=0
	DO I=1,NDFL
	N1=NDLFL(1,I)
	N2=NDLFL(2,I)
	X1=NRXFFL(1,N1)+UTFL(6*N1-5)
	Y1=NRXFFL(2,N1)+UTFL(6*N1-4)
	Z1=NRXFFL(3,N1)+UTFL(6*N1-3)
	X2=NRXF(1,N2)+UT(6*N2-5)
	Y2=NRXF(2,N2)+UT(6*N2-4)
	Z2=NRXF(3,N2)+UT(6*N2-3)
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
	 FXCFL=FXCFL+1
	 FXFFL(1,FXCFL)=N1 
	 FXFFL(2,FXCFL)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.000*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+06*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

        IF (myid.lt.(YN-1)*ntasks/YN
     1	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
	FXCFR=0
	DO I=1,NDFR
	N1=NDLFR(1,I)
	N2=NDLFR(2,I)
	X1=NRXFFR(1,N1)+UTFR(6*N1-5)
	Y1=NRXFFR(2,N1)+UTFR(6*N1-4)
	Z1=NRXFFR(3,N1)+UTFR(6*N1-3)
	X2=NRXF(1,N2)+UT(6*N2-5)
	Y2=NRXF(2,N2)+UT(6*N2-4)
	Z2=NRXF(3,N2)+UT(6*N2-3)
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
	 FXCFR=FXCFR+1
	 FXFFR(1,FXCFR)=N1 
	 FXFFR(2,FXCFR)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.000*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+06*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

        IF (myid.ge.ntasks/YN) THEN
	FXCB=0
	DO I=1,NDB
	N1=NDLB(1,I)
	N2=NDLB(2,I)
	X1=NRXFB(1,N1)+UTB(6*N1-5)
	Y1=NRXFB(2,N1)+UTB(6*N1-4)
	Z1=NRXFB(3,N1)+UTB(6*N1-3)
	X2=NRXF(1,N2)+UT(6*N2-5)
	Y2=NRXF(2,N2)+UT(6*N2-4)
	Z2=NRXF(3,N2)+UT(6*N2-3)
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
	 FXCB=FXCB+1
	 FXFB(1,FXCB)=N1 
	 FXFB(2,FXCB)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.000*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+06*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

        IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
	FXCBL=0
	DO I=1,NDBL
	N1=NDLBL(1,I)
	N2=NDLBL(2,I)
	X1=NRXFBL(1,N1)+UTBL(6*N1-5)
	Y1=NRXFBL(2,N1)+UTBL(6*N1-4)
	Z1=NRXFBL(3,N1)+UTBL(6*N1-3)
	X2=NRXF(1,N2)+UT(6*N2-5)
	Y2=NRXF(2,N2)+UT(6*N2-4)
	Z2=NRXF(3,N2)+UT(6*N2-3)
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
	 FXCBL=FXCBL+1
	 FXFBL(1,FXCBL)=N1 
	 FXFBL(2,FXCBL)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.000*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+06*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF

        IF (myid.ge.ntasks/YN
     1	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
	FXCBR=0
	DO I=1,NDBR
	N1=NDLBR(1,I)
	N2=NDLBR(2,I)
	X1=NRXFBR(1,N1)+UTBR(6*N1-5)
	Y1=NRXFBR(2,N1)+UTBR(6*N1-4)
	Z1=NRXFBR(3,N1)+UTBR(6*N1-3)
	X2=NRXF(1,N2)+UT(6*N2-5)
	Y2=NRXF(2,N2)+UT(6*N2-4)
	Z2=NRXF(3,N2)+UT(6*N2-3)
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
	 FXCBR=FXCBR+1
	 FXFBR(1,FXCBR)=N1 
	 FXFBR(2,FXCBR)=N2 
	 ENDIF

         IF (RC.GT.LNN.AND.RC.LT.LNN+0.000*SCL) THEN
         FRX(N2)=FRX(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCX
         FRY(N2)=FRY(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCY
         FRZ(N2)=FRZ(N2)-SCL**2.0*1.0e+06*(LNN-RC)*RCZ
	 WE(N2)=WE(N2)+SCL**2.0*0.5e+06*(LNN-RC)**2.0
         ENDIF

	ENDIF
	ENDDO
        ENDIF


	RETURN
	END






