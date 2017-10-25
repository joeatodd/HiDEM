	SUBROUTINE DT(NNA,NRXF,ND,NDLC,NDRC,NDFC,NDBC,NDFRC,NDBRC,NDFLC,NDBLC
     1	,myid,ntasks,SCL,YN)

	IMPLICIT NONE
        include 'mpif.h'
        include 'param.dat'
        include 'na.dat'
	REAL*8 NRXF(3,NOMA),NRXFL(3,NOMA),NRXFR(3,NOMA)
	REAL*8 NRXFB(3,NOMA),NRXFF(3,NOMA)
	REAL*8 NRXFBL(3,NOMA),NRXFFL(3,NOMA)
	REAL*8 NRXFBR(3,NOMA),NRXFFR(3,NOMA)
	REAL*8 X1,X2,Y1,Y2,Z1,Z2,RC,SCL,NCN(NOMA)
	INTEGER NNA,ND,I,J,PNN(0:5000)
	INTEGER NDRC,NDLC,NDFC,NDBC,YN
	INTEGER NDFRC,NDBRC,NDFLC,NDBLC
        INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),comm
        INTEGER myid,ntasks,ierr
c        CHARACTER*3 na(0:128)
c        data na/
c     1  '000','001','002','003','004','005','006','007','008',
c     1  '009','010','011','012','013','014','015','016','017',
c     1  '018','019','020','021','022','023','024','025','026',
c     1  '027','028','029','030','031','032','033','034','035',
c     1  '036','037','038','039','040','041','042','043','044',
c     1  '045','046','047','048','049','050','051','052','053',
c     1  '054','055','056','057','058','059','060','061','062',
c     1  '063','064','065','066','067','068','069','070','071',
c     1  '072','073','074','075','076','077','078','079','080',
c     1  '081','082','083','084','085','086','087','088','089',
c     1  '090','091','092','093','094','095','096','097','098',
c     1  '099','100','101','102','103','104','105','106','107',
c     1  '108','109','110','111','112','113','114','115','116',
c     1  '117','118','119','120','121','122','123','124','125',
c     1  '126','127','128'/


 11     FORMAT(2I8,' ',7F18.8)
 12     FORMAT(6F18.8)

        OPEN(UNIT=210+myid,FILE='FS'//na(myid),STATUS='UNKNOWN')
        OPEN(UNIT=310+myid,FILE='FSL'//na(myid),STATUS='UNKNOWN')
        OPEN(UNIT=410+myid,FILE='FSR'//na(myid),STATUS='UNKNOWN')

        CALL MPI_ALLGATHER(NNA,1,MPI_INTEGER,
     1  PNN,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

	ND=0
	NDLC=0
	NDRC=0

	DO I=1,NOMA
	NCN(I)=0
	ENDDO

	DO I=1,NNA-1
 	  X1=NRXF(1,I)
	  Y1=NRXF(2,I)
	  Z1=NRXF(3,I)
	 IF (Z1.GT.70.0*(3000.0-Y1)**0.2.OR.Y1.GT.3000.0) THEN
	 DO J=I+1,NNA
	  X2=NRXF(1,J)
	  Y2=NRXF(2,J)
	  Z2=NRXF(3,J)
c	  RC=SQRT((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  IF (RC.LT.SCL*SCL*1.6) THEN
	  ND=ND+1
 	  write(210+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
c 	  write(210+myid,*) X1,Y1,Z1
c 	  write(210+myid,*) X2,Y2,Z2
c 	  write(210+myid,*) 
	  NCN(I)=NCN(I)+1
	  NCN(J)=NCN(J)+1
 	  ENDIF
	 END DO
 	 ENDIF
	END DO


      dest=myid+1
      source=myid-1

        tag=32
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)
     1CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (MOD(myid,ntasks/YN).ne.0)
     1CALL MPI_Recv(NRXFL,3*NOMA,MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid-1
      source=myid+1

        tag=33
      IF (MOD(myid,ntasks/YN).ne.0)
     1CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)
     1CALL MPI_Recv(NRXFR,3*NOMA,MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)

c      IF (myid.ne.0.and.myid.ne.(ntasks-1)/2+1) THEN
      IF (MOD(myid,ntasks/YN).ne.0) THEN
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
	IF (Z2.GT.70.0*(3000.0-Y2)**0.2.OR.Y2.GT.3000.0) THEN
        DO I=1,PNN(myid-1)
          X1=NRXFL(1,I)
          Y1=NRXFL(2,I)
          Z1=NRXFL(3,I)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDLC=NDLC+1
 	  write(310+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
c 	  write(310+myid,*) X1,Y1,Z1
c 	  write(310+myid,*) X2,Y2,Z2
c 	  write(310+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
 	 ENDIF
        END DO
      END IF

c      IF (myid.ne.ntasks-1.and.myid.ne.(ntasks-1)/2) THEN
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO I=1,PNN(myid+1)
          X1=NRXFR(1,I)
          Y1=NRXFR(2,I)
          Z1=NRXFR(3,I)
 	 IF (Z1.GT.70.0*(3000.0-Y1)**0.2.OR.Y1.GT.3000.0) THEN
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDRC=NDRC+1
 	  write(410+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
c 	  write(410+myid,*) X1,Y1,Z1
c 	  write(410+myid,*) X2,Y2,Z2
c 	  write(410+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
 	 ENDIF
        END DO
      END IF

        CLOSE (210+myid)
        CLOSE (310+myid)
        CLOSE (410+myid)


        OPEN(UNIT=310+myid,FILE='FSB'//na(myid),STATUS='UNKNOWN')
        OPEN(UNIT=410+myid,FILE='FSF'//na(myid),STATUS='UNKNOWN')


	NDFC=0
	NDBC=0

      dest=myid+ntasks/YN
      source=myid-ntasks/YN
	
        tag=34
      IF (myid.lt.(YN-1)*ntasks/YN)
     1CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.ge.ntasks/YN)
     1CALL MPI_Recv(NRXFB,3*NOMA,MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid-ntasks/YN
      source=myid+ntasks/YN

        tag=35
      IF (myid.ge.ntasks/YN)
     1CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN)
     1CALL MPI_Recv(NRXFF,3*NOMA,MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)


      IF (myid.ge.ntasks/YN) THEN
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
	IF (Z2.GT.70.0*(3000.0-Y2)**0.2.OR.Y2.GT.3000.0) THEN
        DO I=1,PNN(myid-ntasks/YN)
          X1=NRXFB(1,I)
          Y1=NRXFB(2,I)
          Z1=NRXFB(3,I)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDBC=NDBC+1
 	  write(310+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
c 	  write(310+myid,*) X1,Y1,Z1
c 	  write(310+myid,*) X2,Y2,Z2
c 	  write(310+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
 	 ENDIF
        END DO
      END IF

      IF (myid.lt.(YN-1)*ntasks/YN) THEN
        DO I=1,PNN(ntasks/YN+myid)
          X1=NRXFF(1,I)
          Y1=NRXFF(2,I)
          Z1=NRXFF(3,I)
 	 IF (Z1.GT.70.0*(3000.0-Y1)**0.2.OR.Y1.GT.3000.0) THEN
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDFC=NDFC+1
 	  write(410+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
c 	  write(410+myid,*) X1,Y1,Z1
c 	  write(410+myid,*) X2,Y2,Z2
c 	  write(410+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
 	 ENDIF
        END DO
      END IF

c	write(*,*) 'NDR',NDRC,myid
c	write(*,*) 'NDL',NDLC,myid
c	write(*,*) 'ND',ND,myid

        CLOSE (310+myid)
        CLOSE (410+myid)


        OPEN(UNIT=310+myid,FILE='FSBR'//na(myid),STATUS='UNKNOWN')
        OPEN(UNIT=410+myid,FILE='FSFR'//na(myid),STATUS='UNKNOWN')


	NDFRC=0
	NDBRC=0

      dest=myid+ntasks/YN-1
      source=myid-ntasks/YN+1
	
        tag=36
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)
     1CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)
     1CALL MPI_Recv(NRXFBR,3*NOMA,MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid-ntasks/YN-1
      source=myid+ntasks/YN+1

        tag=37
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)
     1CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN
     1.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)
     1CALL MPI_Recv(NRXFFR,3*NOMA,MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)


      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
	IF (Z2.GT.70.0*(3000.0-Y2)**0.2.OR.Y2.GT.3000.0) THEN
        DO I=1,PNN(myid-ntasks/YN+1)
          X1=NRXFBR(1,I)
          Y1=NRXFBR(2,I)
          Z1=NRXFBR(3,I)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDBRC=NDBRC+1
 	  write(310+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
c 	  write(310+myid,*) X1,Y1,Z1
c 	  write(310+myid,*) X2,Y2,Z2
c 	  write(310+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
 	 ENDIF
        END DO
      END IF

      IF (myid.lt.(YN-1)*ntasks/YN
     1.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO I=1,PNN(myid+ntasks/YN+1)
          X1=NRXFFR(1,I)
          Y1=NRXFFR(2,I)
          Z1=NRXFFR(3,I)
	 IF (Z1.GT.70.0*(3000.0-Y1)**0.2.OR.Y1.GT.3000.0) THEN
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDFRC=NDFRC+1
 	  write(410+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
c 	  write(410+myid,*) X1,Y1,Z1
c 	  write(410+myid,*) X2,Y2,Z2
c 	  write(410+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
 	 ENDIF
        END DO
      END IF

c	write(*,*) 'NDR',NDRC,myid
c	write(*,*) 'NDL',NDLC,myid
c	write(*,*) 'ND',ND,myid

        CLOSE (310+myid)
        CLOSE (410+myid)

        OPEN(UNIT=310+myid,FILE='FSBL'//na(myid),STATUS='UNKNOWN')
        OPEN(UNIT=410+myid,FILE='FSFL'//na(myid),STATUS='UNKNOWN')


	NDFLC=0
	NDBLC=0

      dest=myid+ntasks/YN+1
      source=myid-ntasks/YN-1
	
        tag=38
      IF (myid.lt.(YN-1)*ntasks/YN
     1.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)
     1CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)
     1CALL MPI_Recv(NRXFBL,3*NOMA,MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid-ntasks/YN+1
      source=myid+ntasks/YN-1

        tag=39
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)
     1CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,
     1dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)
     1CALL MPI_Recv(NRXFFL,3*NOMA,MPI_DOUBLE_PRECISION,
     1source,tag,MPI_COMM_WORLD,stat,ierr)



      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
	IF (Z2.GT.70.0*(3000.0-Y2)**0.2.OR.Y2.GT.3000.0) THEN
        DO I=1,PNN(myid-ntasks/YN-1)
          X1=NRXFBL(1,I)
          Y1=NRXFBL(2,I)
          Z1=NRXFBL(3,I)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDBLC=NDBLC+1
 	  write(310+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
c 	  write(310+myid,*) X1,Y1,Z1
c 	  write(310+myid,*) X2,Y2,Z2
c 	  write(310+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
 	 ENDIF
        END DO
      END IF


      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO I=1,PNN(myid+ntasks/YN-1)
          X1=NRXFFL(1,I)
          Y1=NRXFFL(2,I)
          Z1=NRXFFL(3,I)
	 IF (Z1.GT.70.0*(3000.0-Y1)**0.2.OR.Y1.GT.3000.0) THEN
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDFLC=NDFLC+1
 	  write(410+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
c 	  write(410+myid,*) X1,Y1,Z1
c 	  write(410+myid,*) X2,Y2,Z2
c 	  write(410+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
 	 ENDIF
        END DO
      END IF

c	write(*,*) 'NDR',NDRC,myid
c	write(*,*) 'NDL',NDLC,myid
c	write(*,*) 'ND',ND,myid

        CLOSE (310+myid)
        CLOSE (410+myid)

        OPEN(UNIT=410+myid,FILE='NCN'//na(myid),STATUS='UNKNOWN')

	DO I=1,NNA
	write(410+myid,*) NRXF(3,I), NCN(I)
	ENDDO


	RETURN
	END





