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

!Generates the info on the connections between particles (within and between partitions)
!If two particles are closer than 1.6*SCL, they are connected. This should be predictable 
!based on structure, but this code is nicely flexible (i.e. no prior info about structure required)
	SUBROUTINE DT(NNA,NRXF,ND,NDLC,NDRC,NDFC,NDBC,NDFRC,NDBRC,NDFLC,NDBLC,myid,ntasks,wrkdir,SCL,YN)

        USE INOUT

	IMPLICIT NONE
        include 'mpif.h'
        include 'param.dat'
        include 'na90.dat'
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
        CHARACTER(LEN=256) :: wrkdir

 11     FORMAT(2I8,' ',7F18.8)
 12     FORMAT(6F18.8)

        OPEN(UNIT=210+myid,FILE=TRIM(wrkdir)//'/FS'//na(myid),STATUS='UNKNOWN')
        OPEN(UNIT=310+myid,FILE=TRIM(wrkdir)//'/FSL'//na(myid),STATUS='UNKNOWN')
        OPEN(UNIT=410+myid,FILE=TRIM(wrkdir)//'/FSR'//na(myid),STATUS='UNKNOWN')

        CALL MPI_ALLGATHER(NNA,1,MPI_INTEGER,PNN,1,MPI_INTEGER,MPI_COMM_ACTIVE,ierr)

	ND=0
	NDLC=0
	NDRC=0

	DO I=1,NOMA
	NCN(I)=0
	ENDDO

	DO I=1,NNA-1
	 DO J=I+1,NNA
 	  X1=NRXF(1,I)
	  Y1=NRXF(2,I)
	  Z1=NRXF(3,I)
	  X2=NRXF(1,J)
	  Y2=NRXF(2,J)
	  Z2=NRXF(3,J)
!	  RC=SQRT((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  IF (RC.LT.SCL*SCL*1.6) THEN
	  ND=ND+1
 	  write(210+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
! 	  write(210+myid,*) X1,Y1,Z1
! 	  write(210+myid,*) X2,Y2,Z2
! 	  write(210+myid,*) 
	  NCN(I)=NCN(I)+1
	  NCN(J)=NCN(J)+1
 	  ENDIF
	 END DO
	END DO


      dest=myid+1
      source=myid-1

        tag=32
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Recv(NRXFL,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-1
      source=myid+1

        tag=33
      IF (MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Recv(NRXFR,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

!      IF (myid.ne.0.and.myid.ne.(ntasks-1)/2+1) THEN
      IF (MOD(myid,ntasks/YN).ne.0) THEN
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
        DO I=1,PNN(myid-1)
          X1=NRXFL(1,I)
          Y1=NRXFL(2,I)
          Z1=NRXFL(3,I)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDLC=NDLC+1
 	  write(310+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
! 	  write(310+myid,*) X1,Y1,Z1
! 	  write(310+myid,*) X2,Y2,Z2
! 	  write(310+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
        END DO
      END IF

!      IF (myid.ne.ntasks-1.and.myid.ne.(ntasks-1)/2) THEN
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO I=1,PNN(myid+1)
          X1=NRXFR(1,I)
          Y1=NRXFR(2,I)
          Z1=NRXFR(3,I)
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDRC=NDRC+1
 	  write(410+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
! 	  write(410+myid,*) X1,Y1,Z1
! 	  write(410+myid,*) X2,Y2,Z2
! 	  write(410+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
        END DO
      END IF

        CLOSE (210+myid)
        CLOSE (310+myid)
        CLOSE (410+myid)


        OPEN(UNIT=310+myid,FILE=TRIM(wrkdir)//'/FSB'//na(myid),STATUS='UNKNOWN')
        OPEN(UNIT=410+myid,FILE=TRIM(wrkdir)//'/FSF'//na(myid),STATUS='UNKNOWN')


	NDFC=0
	NDBC=0

      dest=myid+ntasks/YN
      source=myid-ntasks/YN
	
        tag=34
      IF (myid.lt.(YN-1)*ntasks/YN)&
      CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.ge.ntasks/YN)&
      CALL MPI_Recv(NRXFB,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-ntasks/YN
      source=myid+ntasks/YN

        tag=35
      IF (myid.ge.ntasks/YN)&
      CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN)&
      CALL MPI_Recv(NRXFF,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)


      IF (myid.ge.ntasks/YN) THEN
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
        DO I=1,PNN(myid-ntasks/YN)
          X1=NRXFB(1,I)
          Y1=NRXFB(2,I)
          Z1=NRXFB(3,I)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDBC=NDBC+1
 	  write(310+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
! 	  write(310+myid,*) X1,Y1,Z1
! 	  write(310+myid,*) X2,Y2,Z2
! 	  write(310+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
        END DO
      END IF

      IF (myid.lt.(YN-1)*ntasks/YN) THEN
        DO I=1,PNN(ntasks/YN+myid)
          X1=NRXFF(1,I)
          Y1=NRXFF(2,I)
          Z1=NRXFF(3,I)
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDFC=NDFC+1
 	  write(410+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
! 	  write(410+myid,*) X1,Y1,Z1
! 	  write(410+myid,*) X2,Y2,Z2
! 	  write(410+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
        END DO
      END IF

!	write(*,*) 'NDR',NDRC,myid
!	write(*,*) 'NDL',NDLC,myid
!	write(*,*) 'ND',ND,myid

        CLOSE (310+myid)
        CLOSE (410+myid)


        OPEN(UNIT=310+myid,FILE=TRIM(wrkdir)//'/FSBR'//na(myid),STATUS='UNKNOWN')
        OPEN(UNIT=410+myid,FILE=TRIM(wrkdir)//'/FSFR'//na(myid),STATUS='UNKNOWN')


	NDFRC=0
	NDBRC=0

      dest=myid+ntasks/YN-1
      source=myid-ntasks/YN+1
	
        tag=36
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Recv(NRXFBR,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-ntasks/YN-1
      source=myid+ntasks/YN+1

        tag=37
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN&
      .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Recv(NRXFFR,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)


      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
        DO I=1,PNN(myid-ntasks/YN+1)
          X1=NRXFBR(1,I)
          Y1=NRXFBR(2,I)
          Z1=NRXFBR(3,I)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDBRC=NDBRC+1
 	  write(310+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
! 	  write(310+myid,*) X1,Y1,Z1
! 	  write(310+myid,*) X2,Y2,Z2
! 	  write(310+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
        END DO
      END IF

      IF (myid.lt.(YN-1)*ntasks/YN&
      .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO I=1,PNN(myid+ntasks/YN+1)
          X1=NRXFFR(1,I)
          Y1=NRXFFR(2,I)
          Z1=NRXFFR(3,I)
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDFRC=NDFRC+1
 	  write(410+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
! 	  write(410+myid,*) X1,Y1,Z1
! 	  write(410+myid,*) X2,Y2,Z2
! 	  write(410+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
        END DO
      END IF

!	write(*,*) 'NDR',NDRC,myid
!	write(*,*) 'NDL',NDLC,myid
!	write(*,*) 'ND',ND,myid

        CLOSE (310+myid)
        CLOSE (410+myid)

        OPEN(UNIT=310+myid,FILE=TRIM(wrkdir)//'/FSBL'//na(myid),STATUS='UNKNOWN')
        OPEN(UNIT=410+myid,FILE=TRIM(wrkdir)//'/FSFL'//na(myid),STATUS='UNKNOWN')


	NDFLC=0
	NDBLC=0

      dest=myid+ntasks/YN+1
      source=myid-ntasks/YN-1
	
        tag=38
      IF (myid.lt.(YN-1)*ntasks/YN&
      .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Recv(NRXFBL,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-ntasks/YN+1
      source=myid+ntasks/YN-1

        tag=39
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Send(NRXF,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Recv(NRXFFL,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)



      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
        DO I=1,PNN(myid-ntasks/YN-1)
          X1=NRXFBL(1,I)
          Y1=NRXFBL(2,I)
          Z1=NRXFBL(3,I)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDBLC=NDBLC+1
 	  write(310+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
! 	  write(310+myid,*) X1,Y1,Z1
! 	  write(310+myid,*) X2,Y2,Z2
! 	  write(310+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
        END DO
      END IF


      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO I=1,PNN(myid+ntasks/YN-1)
          X1=NRXFFL(1,I)
          Y1=NRXFFL(2,I)
          Z1=NRXFFL(3,I)
         DO J=1,NNA
          X2=NRXF(1,J)
          Y2=NRXF(2,J)
          Z2=NRXF(3,J)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  NDFLC=NDFLC+1
 	  write(410+myid,11) I,J,X1,Y1,Z1,X2,Y2,Z2,1.0
! 	  write(410+myid,*) X1,Y1,Z1
! 	  write(410+myid,*) X2,Y2,Z2
! 	  write(410+myid,*) 
	  NCN(J)=NCN(J)+1
          ENDIF
         END DO
        END DO
      END IF

!	write(*,*) 'NDR',NDRC,myid
!	write(*,*) 'NDL',NDLC,myid
!	write(*,*) 'ND',ND,myid

        CLOSE (310+myid)
        CLOSE (410+myid)

!        OPEN(UNIT=410+myid,FILE=TRIM(wrkdir)//'/NCN'//na(myid),STATUS='UNKNOWN')

!	DO I=1,NNA
!	write(410+myid,*) NRXF(3,I), NCN(I)
!	ENDDO


	RETURN
	END





