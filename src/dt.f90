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
	SUBROUTINE DT(NNA,xo,ND,wrkdir,SCL,YN)

        USE INOUT
        USE TypeDefs

	IMPLICIT NONE
        include 'mpif.h'
        include 'na90.dat'

	REAL*8 :: X1,X2,Y1,Y2,Z1,Z2,RC,SCL,NCN(NOMA),xo(3,NOMA)
	INTEGER :: NNA,I,J,PNN(0:5000),YN,ierr
	INTEGER :: dest,source,tag,stat(MPI_STATUS_SIZE),comm
	CHARACTER(LEN=256) :: wrkdir
        TYPE(NRXF_t) :: NRXF
        TYPE(NTOT_t) :: ND

 11     FORMAT(2I8,' ',7F18.8)
 12     FORMAT(6F18.8)

        NRXF%M = xo

        OPEN(UNIT=210+myid,FILE=TRIM(wrkdir)//'/FS'//na(myid),STATUS='UNKNOWN')
        OPEN(UNIT=310+myid,FILE=TRIM(wrkdir)//'/FSL'//na(myid),STATUS='UNKNOWN')
        OPEN(UNIT=410+myid,FILE=TRIM(wrkdir)//'/FSR'//na(myid),STATUS='UNKNOWN')

        CALL MPI_ALLGATHER(NNA,1,MPI_INTEGER,PNN,1,MPI_INTEGER,MPI_COMM_ACTIVE,ierr)

	ND%M=0
	ND%L=0
	ND%R=0

	DO I=1,NOMA
	NCN(I)=0
	ENDDO

	DO I=1,NNA-1
	 DO J=I+1,NNA
 	  X1=NRXF%M(1,I)
	  Y1=NRXF%M(2,I)
	  Z1=NRXF%M(3,I)
	  X2=NRXF%M(1,J)
	  Y2=NRXF%M(2,J)
	  Z2=NRXF%M(3,J)
!	  RC=SQRT((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
	  IF (RC.LT.SCL*SCL*1.6) THEN
	  ND%M=ND%M+1
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
      CALL MPI_Send(NRXF%M,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Recv(NRXF%L,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-1
      source=myid+1

        tag=33
      IF (MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Send(NRXF%M,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Recv(NRXF%R,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

!      IF (myid.ne.0.and.myid.ne.(ntasks-1)/2+1) THEN
      IF (MOD(myid,ntasks/YN).ne.0) THEN
         DO J=1,NNA
          X2=NRXF%M(1,J)
          Y2=NRXF%M(2,J)
          Z2=NRXF%M(3,J)
        DO I=1,PNN(myid-1)
          X1=NRXF%L(1,I)
          Y1=NRXF%L(2,I)
          Z1=NRXF%L(3,I)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  ND%L=ND%L+1
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
          X1=NRXF%R(1,I)
          Y1=NRXF%R(2,I)
          Z1=NRXF%R(3,I)
         DO J=1,NNA
          X2=NRXF%M(1,J)
          Y2=NRXF%M(2,J)
          Z2=NRXF%M(3,J)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  ND%R=ND%R+1
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


	ND%F=0
	ND%B=0

      dest=myid+ntasks/YN
      source=myid-ntasks/YN
	
        tag=34
      IF (myid.lt.(YN-1)*ntasks/YN)&
      CALL MPI_Send(NRXF%M,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.ge.ntasks/YN)&
      CALL MPI_Recv(NRXF%B,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-ntasks/YN
      source=myid+ntasks/YN

        tag=35
      IF (myid.ge.ntasks/YN)&
      CALL MPI_Send(NRXF%M,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN)&
      CALL MPI_Recv(NRXF%F,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)


      IF (myid.ge.ntasks/YN) THEN
         DO J=1,NNA
          X2=NRXF%M(1,J)
          Y2=NRXF%M(2,J)
          Z2=NRXF%M(3,J)
        DO I=1,PNN(myid-ntasks/YN)
          X1=NRXF%B(1,I)
          Y1=NRXF%B(2,I)
          Z1=NRXF%B(3,I)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  ND%B=ND%B+1
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
          X1=NRXF%F(1,I)
          Y1=NRXF%F(2,I)
          Z1=NRXF%F(3,I)
         DO J=1,NNA
          X2=NRXF%M(1,J)
          Y2=NRXF%M(2,J)
          Z2=NRXF%M(3,J)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  ND%F=ND%F+1
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


	ND%FR=0
	ND%BR=0

      dest=myid+ntasks/YN-1
      source=myid-ntasks/YN+1
	
        tag=36
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Send(NRXF%M,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Recv(NRXF%BR,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-ntasks/YN-1
      source=myid+ntasks/YN+1

        tag=37
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Send(NRXF%M,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN&
      .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Recv(NRXF%FR,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)


      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
         DO J=1,NNA
          X2=NRXF%M(1,J)
          Y2=NRXF%M(2,J)
          Z2=NRXF%M(3,J)
        DO I=1,PNN(myid-ntasks/YN+1)
          X1=NRXF%BR(1,I)
          Y1=NRXF%BR(2,I)
          Z1=NRXF%BR(3,I)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  ND%BR=ND%BR+1
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
          X1=NRXF%FR(1,I)
          Y1=NRXF%FR(2,I)
          Z1=NRXF%FR(3,I)
         DO J=1,NNA
          X2=NRXF%M(1,J)
          Y2=NRXF%M(2,J)
          Z2=NRXF%M(3,J)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  ND%FR=ND%FR+1
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
!	write(*,*) 'ND%M',ND%M,myid

        CLOSE (310+myid)
        CLOSE (410+myid)

        OPEN(UNIT=310+myid,FILE=TRIM(wrkdir)//'/FSBL'//na(myid),STATUS='UNKNOWN')
        OPEN(UNIT=410+myid,FILE=TRIM(wrkdir)//'/FSFL'//na(myid),STATUS='UNKNOWN')


	ND%FL=0
	ND%BL=0

      dest=myid+ntasks/YN+1
      source=myid-ntasks/YN-1
	
        tag=38
      IF (myid.lt.(YN-1)*ntasks/YN&
      .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Send(NRXF%M,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Recv(NRXF%BL,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)

      dest=myid-ntasks/YN+1
      source=myid+ntasks/YN-1

        tag=39
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Send(NRXF%M,3*NOMA,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_ACTIVE,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Recv(NRXF%FL,3*NOMA,MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_ACTIVE,stat,ierr)



      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
         DO J=1,NNA
          X2=NRXF%M(1,J)
          Y2=NRXF%M(2,J)
          Z2=NRXF%M(3,J)
        DO I=1,PNN(myid-ntasks/YN-1)
          X1=NRXF%BL(1,I)
          Y1=NRXF%BL(2,I)
          Z1=NRXF%BL(3,I)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  ND%BL=ND%BL+1
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
          X1=NRXF%FL(1,I)
          Y1=NRXF%FL(2,I)
          Z1=NRXF%FL(3,I)
         DO J=1,NNA
          X2=NRXF%M(1,J)
          Y2=NRXF%M(2,J)
          Z2=NRXF%M(3,J)
          RC=((X1-X2)**2.0+(Y1-Y2)**2.0+(Z1-Z2)**2.0)
          IF (RC.LT.SCL*SCL*1.6) THEN
	  ND%FL=ND%FL+1
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
!	write(*,*) 'ND%M',ND%M,myid

        CLOSE (310+myid)
        CLOSE (410+myid)

!        OPEN(UNIT=410+myid,FILE=TRIM(wrkdir)//'/NCN'//na(myid),STATUS='UNKNOWN')

!	DO I=1,NNA
!	write(410+myid,*) NRXF%M(3,I), NCN(I)
!	ENDDO


	RETURN
	END





