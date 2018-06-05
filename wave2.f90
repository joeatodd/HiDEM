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


!NOMA = 1e+05 - number of points
!NODM = 6*NOMA - degrees of freedom
!NOCON = 4e+05 - number of beams
!NODC - 12 times beams - linear system of equations describing beam interactions

 	PROGRAM WAVE

        USE INOUT

	IMPLICIT NONE
        INCLUDE 'mpif.h'
        INCLUDE 'param.dat'
        INCLUDE 'na90.dat'
	REAL*8 EMM(NOMA),BOYZ(NOMA),BOYY(NOMA),WSY(NOMA)
	REAL*8 EN(NODM),WE(NOMA),VDP(NOMA),WSX(NOMA),BED(-100:2000,-100:2000)
	REAL*8 SUF(-100:2000,-100:2000)
	REAL*8 MFIL(NOMA),NRXF(3,NOMA),EF(NOCON),EFC(NOMA),EFL(NOCON)
	REAL*8 NRXFL(3,NOMA),NRXFR(3,NOMA),EFR(NOCON)
	REAL*8 NRXFB(3,NOMA),NRXFF(3,NOMA),EFF(NOCON),EFB(NOCON)
	REAL*8 NRXFBL(3,NOMA),UTBL(NODM),EFFR(NOCON),EFBR(NOCON)
	REAL*8 NRXFBR(3,NOMA),UTBR(NODM),EFFL(NOCON),EFBL(NOCON)
	REAL*8 NRXFFL(3,NOMA),UTFL(NODM),FBED(-100:2000,-100:2000)
	REAL*8 NRXFFR(3,NOMA),UTFR(NODM),DIX,DIY,DIZ,FRIC,UC
	REAL*8 UTM(NODM),UT(NODM),UTP(NODM),R(NODM)
	REAL*8 UTR(NODM),UTL(NODM),UTW(NODM),NRXFW(3,NOMA)
	REAL*8 UTB(NODM),UTF(NODM),ENM0,ENMS0,POR,GRID
        REAL*8 RHO,RHOW,GRAV, BedIntConst
        REAL*8 CT(NODC),FRX(NOMA),FRY(NOMA),FRZ(NOMA)
	INTEGER CCN(NOMA),CCNW(NOMA),CB(NOMA),CCB(NOMA,3)
	INTEGER NAN(3,NOCON),NDL(2,NODC),NDLL(2,NODC),NDLR(2,NODC)
	INTEGER NDLF(2,NODC),NDLB(2,NODC),NANW(3,NOCON)
	INTEGER FXF(2,NODC),NANL(3,NOMA),NANR(3,NOMA)
	INTEGER NANF(3,NOMA),NANB(3,NOMA),REST,RY0
	INTEGER NANFL(3,NOMA),NANBL(3,NOMA),PTB(0:5000),PTR(0:5000)
	INTEGER NANFR(3,NOMA),NANBR(3,NOMA),NM2
	INTEGER FXFL(2,NODC),FXFR(2,NODC),FXL,FXR,FXCF,FXCB
	INTEGER FXFF(2,NODC),FXFB(2,NODC),FXCFL,FXCBL,FXCFR,FXCBR
	INTEGER FXFFL(2,NODC),FXFBL(2,NODC)
	INTEGER FXFFR(2,NODC),FXFBR(2,NODC)
	REAL*8 M,MN,JS,DT,T,X,Y,CL,CN,E,GSUM,GSUM0
	REAL*8 DMPEN,PSUM,KIN,KIN2,PRESS,MELT
	INTEGER I,N,NL,NTOT,NN,STEPS,IX,IM,MS,N1,N2,RY
	INTEGER PN,NTOL,NTOR,NTOF,NTOB,NRY,PNN(0:5000),NTOTW(0:5000)
	INTEGER NTOFR,NTOBR,NTOFL,NTOBL,XK,YK,ZK
        REAL*8 L,ALF,MLOAD,DMP,VEL,G,S1,S2,M1,B1,B2
	REAL*8 S,LOAD,DMP2,BCE,STR,I1,I2,ZB,ZS
	REAL*8 T2,T1,TS1,TT1,TT2,ENM,WEN,ERSUM
	REAL*8 TT(11),TTCUM(11)
	INTEGER NNO,NS,SI,NNT,BCCS
	INTEGER YNOD,LS,KK,STEPS0,YN,XN
        INTEGER DST,ZNOD,J,XY,ND,O
	INTEGER P,BCC,FXC,NCL,NDR,NDF,NDB
	INTEGER NDFL,NDBL,NDFR,NDBR
	INTEGER NDLFL(2,NODC),NDLBL(2,NODC)
	INTEGER NDLFR(2,NODC),NDLBR(2,NODC)!,XIND,YIND
	REAL*8 KINS,KINS2,ENMS,MGHS,DMPENS,PSUMS
	REAL*8 WENS,GSUMS,DPES,BCES
	REAL*8 XE,DEX,DEY,RDE,MML,XI,ZI,YI
	REAL*8 Z,AVE,ROUGH,FG,SCA,MGH,MGH0,DPE
	REAL*8 KX(12,12),KY(12,12),KZ(12,12),K(12,12)
	REAL*8 DDX,DDY,DDZ,DX,DY,DZ,DL,DTX,DTY,DTZ
	REAL*8 X1,Y1,Z1,X2,Y2,Z2,DXL,DYL,DZL,DDL,RLS
	REAL*8 MAXX,MAXY,MAXZ,MINY,MINX,MINZ,MAXUT
	REAL*8 V1,V2,V3,MAXV,EF0,GL,WL,SLIN,PI,SUB
	REAL*8 SSB,CSB,SQB,LNN,SCL,DAMP1,DAMP2,DRAG,fractime
	REAL RAN(NOCON)
        INTEGER dest,source,tag,stat(MPI_STATUS_SIZE)
        INTEGER rc,myid,ntasks,ierr,SEED,SEEDI,OUTINT,RESOUTINT
        INTEGER, DIMENSION(8) :: datetime
        LOGICAL :: BedZOnly
        CHARACTER(LEN=256) INFILE, geomfile, runname, wrkdir, resdir,restname

        CALL MPI_INIT(rc)
        IF (rc /= MPI_SUCCESS) THEN
        WRITE(*,*) 'MPI initialisation failed'
        STOP
        END IF
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, rc)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, rc)

IF(myid==0) THEN
  CALL DATE_AND_TIME(VALUES=datetime)
5 FORMAT("HiDEM run starting at: ",i2.2,"/",i2.2,"/",i4.4,' ',i2.2,':',i2.2,':',i2.2)
  PRINT *, '----------------------------------------------------'
  WRITE(*,5) datetime(3),datetime(2),datetime(1),datetime(5),datetime(6),datetime(7)
END IF


 6      FORMAT(F11.6,' ',I8,' ',2F11.6)
 7	FORMAT(7F13.7)
 8	FORMAT (A)
 9	FORMAT(I4,'     ',2F11.6)
 10	FORMAT(5F28.1)
 11	FORMAT(2I8,' ',7F22.9)
 12	FORMAT(4F16.6)
 13	FORMAT(6F22.12)
 14	FORMAT(7I8)
 15	FORMAT(5F20.7)
 16	FORMAT(F11.6,' ',I4)
 17	FORMAT(10I9)
 18	FORMAT(6F16.6)
 19	FORMAT(3F24.4,' ',I8)

!writing out energy output
        OPEN(UNIT=609, FILE='HIDEM_STARTINFO',STATUS='OLD')
        READ(609,*) INFILE
        CLOSE(609)
        !INFILE = 'testinp.dat'

        CALL ReadInput(INFILE, myid, runname, wrkdir, resdir, geomfile, PRESS, MELT, UC, DT, S, GRAV, &
             RHO, RHOW, EF0, LS, SUB, GL, SLIN, MLOAD, FRIC, REST, restname, POR, SEEDI, DAMP1, DAMP2, &
             DRAG, BedIntConst, BedZOnly, OUTINT, RESOUTINT, MAXUT, SCL, WL, STEPS0,GRID,fractime)

   IF(myid==0) THEN
     OPEN(UNIT=610,FILE=TRIM(wrkdir)//'/dtop00',STATUS='UNKNOWN',POSITION='APPEND')
     OPEN(UNIT=611,FILE=TRIM(wrkdir)//'/dtop01',STATUS='UNKNOWN',POSITION='APPEND')
     OPEN(UNIT=612,FILE=TRIM(wrkdir)//'/dtopr',STATUS='UNKNOWN',POSITION='APPEND')
     OPEN(UNIT=613,FILE=TRIM(wrkdir)//'/kins2',STATUS='UNKNOWN',POSITION='APPEND')
     OPEN(UNIT=614,FILE=TRIM(wrkdir)//'/lbound',STATUS='UNKNOWN')
     OPEN(UNIT=615,FILE=TRIM(wrkdir)//'/rbound',STATUS='UNKNOWN')
     OPEN(UNIT=110,FILE=TRIM(wrkdir)//'/fib00',STATUS='UNKNOWN')
     OPEN(UNIT=800+myid,FILE=TRIM(wrkdir)//'/tbed'//na(myid),STATUS='UNKNOWN')
     !	OPEN(UNIT=1700+myid,FILE='bccs'//na(myid),STATUS='UNKNOWN')
   END IF

! S = width/thickness of the beams, scaled by SCL
! S * SCL, scales the whole system up, so beam width * 60, particle size * 60
! The unit case - 1m diam particles, 1.14m long beam, S width beam

!MN - mass of particles
!JS - moment of inertia - heaviness w.r.t. rotation
!LNN - distance (scaled) between particles
!MLOAD - maximum load - bonds break when this is exceeded - tension and bending

	S=S*SCL
	MN=SCL**3.0*RHO
	JS=SCL**2.0*MN/6.0
	LNN=SCL*1.1225  
	MLOAD=SCL*MLOAD

        SEED=SEEDI+873*myid
        CALL RMARIN(SEED,0,0)
        CALL RANMAR(RAN,NOCON)

	PI=ACOS(-1.0)
        DMP=DAMP1*SCL**3.0
        DMP2=DAMP2*SCL**3.0

!mpi stuff - boundaries
	FXC=0
	FXR=0
	FXL=0
	FXCF=0
	FXCB=0
	FXCFL=0
	FXCFR=0
	FXCBL=0
	FXCBR=0

	DO I=1,NOMA
	FXF(1,I)=0
        CCN(I)=0
	EFC(I)=SCL*EF0
	END DO

	DO I=1,NODM
	EN(I)=0.0
	END DO

!accumulative energy terms - don't zero if restarting
	IF (REST.EQ.0) THEN
	BCC=0
	BCE=0.0
	DPE=0.0
	DMPEN=0.0	
	PSUM=0.0
	END IF

!more MPI stuff...
!square partitioning, (F)orward, (B)ack, (L)eft, (R)right, (FR) Forward Right, etc 
	DO J=1,NOMA
	NRXFL(1,J)=-1000.0
	NRXFL(2,J)=-1000.0
	NRXFL(3,J)=-1000.0
	NRXFR(1,J)=-1000.0
	NRXFR(2,J)=-1000.0
	NRXFR(3,J)=-1000.0
	NRXFF(1,J)=-1000.0
	NRXFF(2,J)=-1000.0
	NRXFF(3,J)=-1000.0
	NRXFB(1,J)=-1000.0
	NRXFB(2,J)=-1000.0
	NRXFB(3,J)=-1000.0
	NRXFBL(1,J)=-1000.0
	NRXFBL(2,J)=-1000.0
	NRXFBL(3,J)=-1000.0
	NRXFBR(1,J)=-1000.0
	NRXFBR(2,J)=-1000.0
	NRXFBR(3,J)=-1000.0
	NRXFFL(1,J)=-1000.0
	NRXFFL(2,J)=-1000.0
	NRXFFL(3,J)=-1000.0
	NRXFFR(1,J)=-1000.0
	NRXFFR(2,J)=-1000.0
	NRXFFR(3,J)=-1000.0
	END DO

!inclination of the domain - not really used

	SSB=SIN(SUB*PI/2.0)
	CSB=COS(SUB*PI/2.0)
	SQB=SQRT(1.0+SIN(SUB*PI/2.0)**2)
	RLS=LS

	IF (REST.EQ.0) THEN

	RY0=1

        !Go to glas.f90 to make the grid
	CALL FIBG3(LS,NN,NTOT,NTOL,NTOR,NTOF,NTOB,NTOFR,NTOBR,NTOFL,NTOBL,&
        myid,MAXX,MAXY,MAXZ,MINX,MINY,MINZ,ntasks,wrkdir,geomfile,SCL,YN,XN,GRID,MELT,WL,UC)

	write(*,17) myid,NTOL,NTOT,NTOR,NTOF,NTOB,NTOFL,NTOFR,NTOBL,NTOBR
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
        !NN - particles in each core
        CALL MPI_ALLGATHER(NN,1,MPI_INTEGER,&
        PNN,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        CALL MPI_ALLGATHER(NTOR,1,MPI_INTEGER,&
        PTR,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        CALL MPI_ALLGATHER(NTOB,1,MPI_INTEGER,&
        PTB,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        CALL MPI_ALLGATHER(NTOT,1,MPI_INTEGER,&
        NTOTW,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)


        !MATHS!
        !Iterate over: this part nodes, edge nodes
        ! CT - see notes - it's the cumulative translation of all the beams
	DO I=1,NTOT+NTOL+NTOR+NTOF+NTOB+NTOFL+NTOFR+NTOBL+NTOBR
	CT(12*I-11)=0.0
	CT(12*I-10)=0.0
	CT(12*I-9)=0.0
	CT(12*I-8)=0.0
	CT(12*I-7)=0.0
	CT(12*I-6)=0.0
	CT(12*I-5)=0.0
	CT(12*I-4)=0.0
	CT(12*I-3)=0.0
	CT(12*I-2)=0.0
	CT(12*I-1)=0.0
	CT(12*I-0)=0.0
	END DO

	T=0 !T is time
	ENM=0.0
	WEN=0.0
	ENM0=0

	DO J=1,NN
        VDP(J)=0.0     !VDP = drag coefficient
	UT(6*J-5)=0.0  ! UT = current displacement
	UT(6*J-4)=0.0
	UT(6*J-3)=0.0
	UT(6*J-2)=0.0
	UT(6*J-1)=0.0
	UT(6*J-0)=0.0
	UTM(6*J-5)=0.0 ! UTM = previous displacement
	UTM(6*J-4)=0.0
	UTM(6*J-3)=0.0
	UTM(6*J-2)=0.0
	UTM(6*J-1)=0.0
	UTM(6*J-0)=0.0
	END DO

	ELSE

        !If restarting, read from restart file instead
	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(restname)//'_REST0'//na(myid),STATUS='OLD')

	READ(117+myid,*) NN,NTOT,NTOL,NTOR,NTOF,NTOB,BCC
	READ(117+myid,*) NTOFL,NTOFR,NTOBL,NTOBR
	READ(117+myid,*) MAXX,MAXY,MAXZ,DMPEN,ENM0
	READ(117+myid,*) DPE,BCE,MGH0,GSUM0,PSUM,T,RY0
	READ(117+myid,*) XN,YN
	CLOSE(117+myid)

        CALL MPI_ALLGATHER(NN,1,MPI_INTEGER,&
        PNN,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        CALL MPI_ALLGATHER(NTOR,1,MPI_INTEGER,&
        PTR,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        CALL MPI_ALLGATHER(NTOB,1,MPI_INTEGER,&
        PTB,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
 
        CALL MPI_ALLGATHER(NTOT,1,MPI_INTEGER,&
        NTOTW,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(restname)//'_REST1'//na(myid),STATUS='OLD')
	DO I=1,NTOT+NTOL+NTOR+NTOF+NTOB+NTOFL+NTOFR+NTOBL+NTOBR
	READ(117+myid,*) CT(12*I-11),CT(12*I-10),CT(12*I-9),&
      	CT(12*I-8),CT(12*I-7),CT(12*I-6)
	READ(117+myid,*) CT(12*I-5),CT(12*I-4),CT(12*I-3),&
      	CT(12*I-2),CT(12*I-1),CT(12*I-0)
	END DO
	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(restname)//'_REST2'//na(myid),STATUS='OLD')
	DO I=1,NN
	READ(117+myid,*) UT(6*I-5),UT(6*I-4),UT(6*I-3),&
      	UT(6*I-2),UT(6*I-1),UT(6*I-0)
	READ(117+myid,*) UTM(6*I-5),UTM(6*I-4),UTM(6*I-3),&
      	UTM(6*I-2),UTM(6*I-1),UTM(6*I-0)
	END DO
	CLOSE (117+myid)

	END IF

        !Set bed -1000.0 everywhere
        DO I=-100,2000
        DO J=-100,2000
        BED(I,J)=-1000.0
        ENDDO
        ENDDO

        !Read the geometry and friction
        !into the grids BED, SUF, and FBED
        OPEN(UNIT=400,file=TRIM(geomfile),STATUS='UNKNOWN')
        READ(400,*) NM2
	DO I=1,NM2
        READ(400,*) X,Y,S1,B2,B1,Z1
!        X=X-2000.0
!        Y=Y-7000.0
        XK=INT(X/GRID)
        YK=INT(Y/GRID)
        IF (XK.GE.-100.AND.YK.GE.-100) BED(XK,YK)=B1
        IF (XK.GE.-100.AND.YK.GE.-100) SUF(XK,YK)=S1
        IF (XK.GE.-100.AND.YK.GE.-100) FBED(XK,YK)=FRIC*SCL*SCL*Z1

        !Previously had commented from here
        !------------
        IF (YK.EQ.0.AND.XK.GE.-100) THEN
        DO J=0,20
        BED(XK,-J)=B1
        SUF(XK,-J)=S1
        FBED(XK,-J)=FRIC*SCL*SCL*Z1
        ENDDO
        ENDIF

        IF (XK.EQ.0.AND.YK.GE.-100) THEN
        DO J=0,20
        BED(-J,YK)=B1
        SUF(-J,YK)=S1
        FBED(-J,YK)=FRIC*SCL*SCL*Z1
        ENDDO
        ENDIF
        !To here
        !------------

	ENDDO
	CLOSE(400)

 !Interpolation functions to smooth input data for particles

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/NODFIL2'//na(myid),STATUS='UNKNOWN')
	DO I=1,NN
	READ(117+myid,*) IX,X,Y,Z,M
	NRXF(1,IX)=X
	NRXF(2,IX)=Y
	NRXF(3,IX)=Z
        MFIL(IX)=MN
	X=NRXF(1,IX)+UT(6*IX-5)
        Y=NRXF(2,IX)+UT(6*IX-4)
        Z=NRXF(3,IX)+UT(6*IX-3)
        I1=X/GRID-INT(X/GRID)
	I2=Y/GRID-INT(Y/GRID)
        CALL BIPINT(I1,I2,BED(INT(X/GRID),INT(Y/GRID)),BED(INT(X/GRID)+1,INT(Y/GRID)),&
        BED(INT(X/GRID),INT(Y/GRID)+1),BED(INT(X/GRID)+1,INT(Y/GRID)+1),ZB)
        
        IF (ABS(ZB-Z).LT.SCL*2.0) THEN
        CALL BIPINT(I1,I2,FBED(INT(X/GRID),INT(Y/GRID)),FBED(INT(X/GRID)+1,INT(Y/GRID)),&
        FBED(INT(X/GRID),INT(Y/GRID)+1),FBED(INT(X/GRID)+1,INT(Y/GRID)+1),VDP(IX))
!        IF (VDP(IX).GT.SCL*SCL*2.0e+07) VDP(IX)=SCL*SCL*2.0e+07
        ELSE
          IF (Z.LT.WL) THEN
          VDP(IX)=SCL*SCL*DRAG
          ELSE
          VDP(IX)=SCL*SCL*DRAG
          ENDIF
        ENDIF


	END DO
	CLOSE (117+myid)

 	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FS'//na(myid),STATUS='UNKNOWN')
        DO I=1,NTOT
          READ(117+myid,*) N1,N2,X1,Y1,Z1,X2,Y2,Z2,E
        NAN(1,I)=N1
        NAN(2,I)=N2
	I1=X1/GRID-INT(X1/GRID)
	I2=Y1/GRID-INT(Y1/GRID)
        CALL BIPINT(I1,I2,SUF(INT(X1/GRID),INT(Y1/GRID)),SUF(INT(X1/GRID)+1,INT(Y1/GRID)),&
        SUF(INT(X1/GRID),INT(Y1/GRID)+1),SUF(INT(X1/GRID)+1,INT(Y1/GRID)+1),ZS)
	 IF (REST.EQ.0) THEN
          IF (RAN(I).LT.1.0-POR.AND.ABS(Z1-ZS).LT.SLIN) THEN
 	  EF(I)=EF0
	  ELSE
 	  EF(I)=0.1
	  ENDIF
	 ELSE
	 EF(I)=E
	 ENDIF
	END DO
        CLOSE (117+myid)

        IF (MOD(myid,ntasks/YN).ne.0) THEN
 	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSL'//na(myid),STATUS='UNKNOWN')
        DO I=1,NTOL
        READ(117+myid,*) N1,N2,X1,Y1,Z1,X2,Y2,Z2,E
        NANL(1,I)=N1
        NANL(2,I)=N2
	NRXFL(1,N1)=X1
	NRXFL(2,N1)=Y1
	NRXFL(3,N1)=Z1
	END DO
        CLOSE (117+myid)
	ENDIF

        IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
 	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSR'//na(myid),STATUS='UNKNOWN')
        DO I=1,NTOR
        READ(117+myid,*) N1,N2,X1,Y1,Z1,X2,Y2,Z2,E
        NANR(1,I)=N1
        NANR(2,I)=N2
	NRXFR(1,N1)=X1
	NRXFR(2,N1)=Y1
	NRXFR(3,N1)=Z1
	I1=X1/GRID-INT(X1/GRID)
	I2=Y1/GRID-INT(Y1/GRID)
        CALL BIPINT(I1,I2,SUF(INT(X1/GRID),INT(Y1/GRID)),SUF(INT(X1/GRID)+1,INT(Y1/GRID)),&
        SUF(INT(X1/GRID),INT(Y1/GRID)+1),SUF(INT(X1/GRID)+1,INT(Y1/GRID)+1),ZS)
	 IF (REST.EQ.0) THEN
          !Setting porosity in domain - pre-existing damage
          !POR - porosity
	  IF (RAN(I).LT.1.0-POR.AND.ABS(Z1-ZS).LT.SLIN) THEN
 	  EFR(I)=EF0
	  ELSE
 	  EFR(I)=0.1
          ENDIF
	 ELSE
	 EFR(I)=E
	 ENDIF
	END DO
        CLOSE (117+myid)
	ENDIF

        dest=myid+1
        source=myid-1
        tag=131
        IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
        CALL MPI_Send(EFR,NTOR,MPI_DOUBLE_PRECISION,&
        dest,tag,MPI_COMM_WORLD,ierr)
        IF (MOD(myid,ntasks/YN).ne.0)&
        CALL MPI_Recv(EFL,NTOL,MPI_DOUBLE_PRECISION,&
        source,tag,MPI_COMM_WORLD,stat,ierr)

        IF (myid.lt.(YN-1)*ntasks/YN) THEN
 	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSF'//na(myid),STATUS='UNKNOWN')
        DO I=1,NTOF
        READ(117+myid,*) N1,N2,X1,Y1,Z1,X2,Y2,Z2,E
        NANF(1,I)=N1
        NANF(2,I)=N2
	NRXFF(1,N1)=X1
	NRXFF(2,N1)=Y1
	NRXFF(3,N1)=Z1
	I1=X1/GRID-INT(X1/GRID)
	I2=Y1/GRID-INT(Y1/GRID)
        CALL BIPINT(I1,I2,SUF(INT(X1/GRID),INT(Y1/GRID)),SUF(INT(X1/GRID)+1,INT(Y1/GRID)),&
        SUF(INT(X1/GRID),INT(Y1/GRID)+1),SUF(INT(X1/GRID)+1,INT(Y1/GRID)+1),ZS)
	 IF (REST.EQ.0) THEN
	  IF (RAN(I).LT.1.0-POR.AND.ABS(Z1-ZS).LT.SLIN) THEN
 	  EFF(I)=EF0
	  ELSE
 	  EFF(I)=0.1
          ENDIF
	 ELSE
	 EFF(I)=E
	 ENDIF
	END DO
        CLOSE (117+myid)
	ENDIF

        dest=myid+ntasks/YN
        source=myid-ntasks/YN
        tag=132
        IF (myid.lt.(YN-1)*ntasks/YN)&
        CALL MPI_Send(EFF,NTOF,MPI_DOUBLE_PRECISION,&
        dest,tag,MPI_COMM_WORLD,ierr)
        IF (myid.ge.ntasks/YN)&
        CALL MPI_Recv(EFB,NTOB,MPI_DOUBLE_PRECISION,&
        source,tag,MPI_COMM_WORLD,stat,ierr)

        IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
 	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSFL'//na(myid),STATUS='UNKNOWN')
        DO I=1,NTOFL
        READ(117+myid,*) N1,N2,X1,Y1,Z1,X2,Y2,Z2,E
        NANFL(1,I)=N1
        NANFL(2,I)=N2
	NRXFFL(1,N1)=X1
	NRXFFL(2,N1)=Y1
	NRXFFL(3,N1)=Z1
	I1=X1/GRID-INT(X1/GRID)
	I2=Y1/GRID-INT(Y1/GRID)
        CALL BIPINT(I1,I2,SUF(INT(X1/GRID),INT(Y1/GRID)),SUF(INT(X1/GRID)+1,INT(Y1/GRID)),&
        SUF(INT(X1/GRID),INT(Y1/GRID)+1),SUF(INT(X1/GRID)+1,INT(Y1/GRID)+1),ZS)
	 IF (REST.EQ.0) THEN
	  IF (RAN(I).LT.1.0-POR.AND.ABS(Z1-ZS).LT.SLIN) THEN
 	  EFFL(I)=EF0
	  ELSE
 	  EFFL(I)=0.1
          ENDIF
	 ELSE
	 EFFL(I)=E
	 ENDIF
        END DO
        CLOSE (117+myid)
       ENDIF

        dest=myid+ntasks/YN-1
        source=myid-ntasks/YN+1
        tag=133
        IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
        CALL MPI_Send(EFFL,NTOFL,MPI_DOUBLE_PRECISION,&
        dest,tag,MPI_COMM_WORLD,ierr)
        IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
        CALL MPI_Recv(EFBR,NTOBR,MPI_DOUBLE_PRECISION,&
        source,tag,MPI_COMM_WORLD,stat,ierr)


        IF (myid.lt.(YN-1)*ntasks/YN&
      	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
 	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSFR'//na(myid),STATUS='UNKNOWN')
        DO I=1,NTOFR
        READ(117+myid,*) N1,N2,X1,Y1,Z1,X2,Y2,Z2,E
        NANFR(1,I)=N1
        NANFR(2,I)=N2
	NRXFFR(1,N1)=X1
	NRXFFR(2,N1)=Y1
	NRXFFR(3,N1)=Z1
	I1=X1/GRID-INT(X1/GRID)
	I2=Y1/GRID-INT(Y1/GRID)
        CALL BIPINT(I1,I2,SUF(INT(X1/GRID),INT(Y1/GRID)),SUF(INT(X1/GRID)+1,INT(Y1/GRID)),&
        SUF(INT(X1/GRID),INT(Y1/GRID)+1),SUF(INT(X1/GRID)+1,INT(Y1/GRID)+1),ZS)
	 IF (REST.EQ.0) THEN
	  IF (RAN(I).LT.1.0-POR.AND.ABS(Z1-ZS).LT.SLIN) THEN
 	  EFFR(I)=EF0
	  ELSE
 	  EFFR(I)=0.1
          ENDIF
	 ELSE
	 EFFR(I)=E
	 ENDIF
	END DO
        CLOSE (117+myid)
       ENDIF


       !Read the original positions of particles

        dest=myid+ntasks/YN+1
        source=myid-ntasks/YN-1
        tag=134
        IF (myid.lt.(YN-1)*ntasks/YN&
        .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
        CALL MPI_Send(EFFR,NTOFR,MPI_DOUBLE_PRECISION,&
        dest,tag,MPI_COMM_WORLD,ierr)
        IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
        CALL MPI_Recv(EFBL,NTOBL,MPI_DOUBLE_PRECISION,&
        source,tag,MPI_COMM_WORLD,stat,ierr)

        IF (myid.ge.ntasks/YN) THEN
 	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSB'//na(myid),STATUS='UNKNOWN')
        DO I=1,NTOB
        READ(117+myid,*) N1,N2,X1,Y1,Z1,X2,Y2,Z2,E
        NANB(1,I)=N1
        NANB(2,I)=N2
	NRXFB(1,N1)=X1
	NRXFB(2,N1)=Y1
	NRXFB(3,N1)=Z1
	END DO
        CLOSE (117+myid)
	ENDIF

        IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
 	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSBL'//na(myid),STATUS='UNKNOWN')
        DO I=1,NTOBL
        READ(117+myid,*) N1,N2,X1,Y1,Z1,X2,Y2,Z2,E
        NANBL(1,I)=N1
        NANBL(2,I)=N2
	NRXFBL(1,N1)=X1
	NRXFBL(2,N1)=Y1
	NRXFBL(3,N1)=Z1
	END DO
        CLOSE (117+myid)
	ENDIF

        IF (myid.ge.ntasks/YN&
      	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
 	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSBR'//na(myid),STATUS='UNKNOWN')
        DO I=1,NTOBR
        READ(117+myid,*) N1,N2,X1,Y1,Z1,X2,Y2,Z2,E
        NANBR(1,I)=N1
        NANBR(2,I)=N2
	NRXFBR(1,N1)=X1
	NRXFBR(2,N1)=Y1
	NRXFBR(3,N1)=Z1
	END DO
        CLOSE (117+myid)
	ENDIF

        CALL CPU_TIME(T1)
	TS1=0.0

        TTCUM=0.0
        TT=0.0
	DO 100 RY=RY0,RY0+STEPS0 !START THE TIME LOOP

        CALL CPU_TIME(TT1)
        TT(1) = TT1
        TTCUM(1) = TTCUM(1) + (TT(1) - TT(11))

        !TODO  TIME 1

      dest=myid+1
      source=myid-1

!Send positions of particles to neighbouring partitions

        tag=142
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Send(UT,6*NN,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_WORLD,ierr)
      IF (MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Recv(UTL,6*PNN(source),MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid-1
      source=myid+1

        tag=144
      IF (MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Send(UT,6*NN,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_WORLD,ierr)
      IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Recv(UTR,6*PNN(source),MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid+ntasks/YN
      source=myid-ntasks/YN

        tag=145
      IF (myid.lt.(YN-1)*ntasks/YN)&
      CALL MPI_Send(UT,6*NN,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.ge.ntasks/YN)&
      CALL MPI_Recv(UTB,6*PNN(source),MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid-ntasks/YN
      source=myid+ntasks/YN

        tag=146
      IF (myid.ge.ntasks/YN)&
      CALL MPI_Send(UT,6*NN,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN)&
      CALL MPI_Recv(UTF,6*PNN(source),MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid-ntasks/YN-1
      source=myid+ntasks/YN+1

        tag=147
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Send(UT,6*NN,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN&
      .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Recv(UTFR,6*PNN(source),MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid-ntasks/YN+1
      source=myid+ntasks/YN-1

        tag=148
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Send(UT,6*NN,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Recv(UTFL,6*PNN(source),MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid+ntasks/YN-1
      source=myid-ntasks/YN+1

        tag=149
      IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Send(UT,6*NN,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Recv(UTBR,6*PNN(source),MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_WORLD,stat,ierr)

      dest=myid+ntasks/YN+1
      source=myid-ntasks/YN-1

        tag=150
      IF (myid.lt.(YN-1)*ntasks/YN&
      .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      CALL MPI_Send(UT,6*NN,MPI_DOUBLE_PRECISION,&
      dest,tag,MPI_COMM_WORLD,ierr)
      IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
      CALL MPI_Recv(UTBL,6*PNN(source),MPI_DOUBLE_PRECISION,&
      source,tag,MPI_COMM_WORLD,stat,ierr)

      !TODO  TIME 2
      CALL CPU_TIME(TT(2))
      TTCUM(2) = TTCUM(2) + (TT(2) - TT(1))

      !Every 250 steps, reconstruct neighbourhood list
	IF (MOD(RY,250).EQ.1.OR.RY.EQ.RY0) THEN
     	CALL DIST(NN,UT,ND,NCL,NDR,NRXF,NDL,&
      	NDLL,NDLR,NRXFL,NRXFR,UTL,UTR,NRXFF,&
      	NRXFB,UTB,UTF,NDLF,NDLB,NDF,NDB,myid,ntasks,SCL,PNN,YN,&
        NRXFBL,NRXFBR,NRXFFL,NRXFFR,NDLFL,NDLFR,NDLBR,NDLBL,&
      	NDFL,NDFR,NDBR,NDBL,UTBL,UTBR,UTFL,UTFR)
	END IF

!      write(*,17) RY,myid,ND,NCL,NDR,NDF,NDB,NDBL,NDBL,NDFR

       !TODO  TIME 3
       CALL CPU_TIME(TT(3))
       TTCUM(3) = TTCUM(3) + (TT(3) - TT(2))

       !circ checks which particles are really in contact and computes the forces
	CALL CIRC(ND,NN,NRXF,UT,FRX,FRY,FRZ,&
      	T,RY,DT,WE,EFC,FXF,FXC,NDR,NCL,NDL,&
      	NDLL,NDLR,NRXFL,NRXFR,UTL,UTR,myid,ntasks,FXFL,FXFR,FXL,FXR,&
      	FXCF,FXFF,UTF,NRXFF,NDLF,NDF,&
        FXCB,FXFB,UTB,NRXFB,NDLB,NDB,LNN,YN,&
        NDFL,NDFR,NDBL,NDBR,&
      	FXCFR,FXCFL,FXCBR,FXCBL,&
      	NDLFL,NDLFR,NDLBL,NDLBR,&
      	UTFL,UTFR,UTBL,UTBR,&
      	NRXFBL,NRXFBR,NRXFFL,NRXFFR,&
        FXFFR,FXFFL,FXFBR,FXFBL,SCL)

!      write(*,17) RY,myid,FXC,FXL,FXCB,FXCF,FXCFR,FXCFL,FXCBR,FXCBL

       !TODO  TIME 4
       CALL CPU_TIME(TT(4))
       TTCUM(4) = TTCUM(4) + (TT(4) - TT(3))

       !Calculates elastic forces from beams. Stiffness matrix K
	CALL EFFLOAD(S,NTOT,NN,T,DT,MN,JS,DMP,DMP2,UT,UTM,R,EN,RY,&
      	FXF,FXC,VDP,DPE,EF,EFL,EFR,NAN,NRXF,MFIL,CT,NANL,NANR,NTOL,NTOR,&
      	NRXFL,NRXFR,UTL,UTR,myid,ntasks,FXL,FXFL,FXR,FXFR,LNN,&
        NTOF,NTOB,NRXFF,NRXFB,UTF,UTB,EFF,EFB,&
        FXCF,FXFF,FXCB,FXFB,NANB,NANF,PNN,YN,&
      	NANFL,NANFR,NANBL,NANBR,&
      	NRXFFL,NRXFFR,NRXFBL,NRXFBR,&
      	UTFL,UTFR,UTBL,UTBR,&
      	EFFL,EFFR,EFBL,EFBR,&
        NTOFL,NTOFR,NTOBL,NTOBR,&
        FXCFL,FXFFL,FXCBL,FXFBL,&
        FXCFR,FXFFR,FXCBR,FXFBR)


        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

        CALL CPU_TIME(TT2)
	TS1=TS1+(TT2-TT1)

        !TODO  TIME 5
       CALL CPU_TIME(TT(5))
       TTCUM(5) = TTCUM(5) + (TT(5) - TT(4))

	WEN=0.0
	KIN=0.0
	KIN2=0.0
	MGH=0.0
	GSUM=0.0
	ENM=0.0

	DO 27 I=1,NN
        !TODO  TIME 6
        CALL CPU_TIME(TT(6))
        IF(I>1) THEN
          TTCUM(6) = TTCUM(6) + (TT(6) - TT(9))
        END IF

	X=NRXF(1,I)+UT(6*I-5)  !<- x,y,z actual positions, NRXF = original, UT = displacement
	Y=NRXF(2,I)+UT(6*I-4)
	Z=NRXF(3,I)+UT(6*I-3)

	WSY(I)=0.0

! If you want to add a pressure from the backwall, instead of just free.
! Something like hydrostatic pressure - scaled

!	IF (Y.LT.1000.0.AND.(X.GT.8000.0.AND.X.LT.34000.0)) THEN
!          IF (Z.GT.WL+Y*SSB) THEN
!          WSY(I)=PRESS*1.0e+01*SCL**3*(MAXZ-NRXF(3,I))
!          PSUM=PSUM+WSY(I)*(UT(6*I-4)-UTM(6*I-4))
!          ELSE
!          WSY(I)=PRESS*(1.0e+01*SCL**3*(MAXZ-NRXF(3,I))-1.1e+01*SCL**3*(WL+Y*SSB-NRXF(3,I)))
!          PSUM=PSUM+WSY(I)*(UT(6*I-4)-UTM(6*I-4))
!          ENDIF
!	ENDIF

	WSX(I)=0.0

	I1=X/GRID-INT(X/GRID)
	I2=Y/GRID-INT(Y/GRID)

       !Update the drag calculation
        CALL BIPINT(I1,I2,BED(INT(X/GRID),INT(Y/GRID)),BED(INT(X/GRID)+1,INT(Y/GRID)),&
        BED(INT(X/GRID),INT(Y/GRID)+1),BED(INT(X/GRID)+1,INT(Y/GRID)+1),ZB)
        CALL BIPINT(I1,I2,SUF(INT(X/GRID),INT(Y/GRID)),SUF(INT(X/GRID)+1,INT(Y/GRID)),&
        SUF(INT(X/GRID),INT(Y/GRID)+1),SUF(INT(X/GRID)+1,INT(Y/GRID)+1),ZS)

        IF (ABS(ZB-Z).LT.SCL*2.0) THEN
        CALL BIPINT(I1,I2,FBED(INT(X/GRID),INT(Y/GRID)),FBED(INT(X/GRID)+1,INT(Y/GRID)),&
        FBED(INT(X/GRID),INT(Y/GRID)+1),FBED(INT(X/GRID)+1,INT(Y/GRID)+1),VDP(I))
!        IF (VDP(I).GT.SCL*SCL*2.0e+07) VDP(I)=SCL*SCL*2.0e+07
        ELSE
          IF (Z.LT.WL) THEN
          VDP(I)=SCL*SCL*DRAG
          ELSE
          VDP(I)=SCL*SCL*DRAG
          ENDIF
        ENDIF

        !TODO  TIME 7
       CALL CPU_TIME(TT(7))
       TTCUM(7) = TTCUM(7) + (TT(7) - TT(6))

        CALL BIPINTN(I1,I2,BED(INT(X/GRID),INT(Y/GRID)),BED(INT(X/GRID)+1,INT(Y/GRID)),&
        BED(INT(X/GRID),INT(Y/GRID)+1),BED(INT(X/GRID)+1,INT(Y/GRID)+1),DIX,DIY,DIZ,GRID)

!Bed interaction
!Bed Normal
! FR* - forces on particles
! DIX - components of bed normal
        IF (Z.LT.ZB+SCL/2.0) THEN
          IF(BedZOnly) THEN
            FRZ(I)=FRZ(I)+BedIntConst*(ZB+SCL/2.0-Z)
            GSUM=GSUM+BedIntConst*(ZB+SCL/2.0-Z)**2
          ELSE
            FRX(I)=FRX(I)+DIX*BedIntConst*(ZB+SCL/2.0-Z)
            FRY(I)=FRY(I)+DIY*BedIntConst*(ZB+SCL/2.0-Z)
            FRZ(I)=FRZ(I)+DIZ*BedIntConst*(ZB+SCL/2.0-Z)
            GSUM=GSUM+BedIntConst*(ZB+SCL/2.0-Z)**2
          END IF
	ENDIF

!        IF (Y.GT.1e+04.AND.T.GT.200.0) THEN
!      	FRY(I)=FRY(I)+5e+03*(Y-1e+04)
!	ENDIF


       !Buoyancy calculation
        !Note, this is gravity, not just buoyancy!
	IF (Z.LT.WL+Y*SSB) THEN
        BOYZ(I)=GRAV*(-1.0 + (RHOW/RHO))*MFIL(I)*CSB
	ELSE
	BOYZ(I)=-GRAV*MFIL(I)*CSB
	ENDIF

	IF (Z.LT.WL+Y*SSB) THEN
        BOYY(I)=-GRAV*(-1.0 + (RHOW/RHO))*MFIL(I)*SSB
	ELSE
	BOYY(I)=GRAV*MFIL(I)*SSB
	ENDIF

        !TODO  TIME 8
       CALL CPU_TIME(TT(8))
       TTCUM(8) = TTCUM(8) + (TT(8) - TT(7))

 !Jan's code for smoothing the buoyant forces across the waterline
        ! IF (Z.LT.WL-0.4*SCL) THEN
        !    BOYZ(I)=1.29*MFIL(I)
        ! ELSE
        !    IF (Z.GT.WL+0.1*SCL) THEN
        !       BOYZ(I)=-9.8*MFIL(I)
        !    ELSE
        !       KK=(1.29*MFIL(I)+9.8*MFIL(I))/(0.5*SCL)
        !       BOYZ(I)=1.29*MFIL(I)-KK*(Z-(WL-0.4*SCL))
        !    ENDIF
        ! ENDIF

!Kinetic energy
	KIN=KIN+0.5*MN*((UT(6*I-5)-UTM(6*I-5))/DT)**2
	KIN=KIN+0.5*MN*((UT(6*I-4)-UTM(6*I-4))/DT)**2
	KIN=KIN+0.5*MN*((UT(6*I-3)-UTM(6*I-3))/DT)**2
	KIN2=KIN2+0.5*JS*((UT(6*I-2)-UTM(6*I-2))/DT)**2
	KIN2=KIN2+0.5*JS*((UT(6*I-1)-UTM(6*I-1))/DT)**2
	KIN2=KIN2+0.5*JS*((UT(6*I-0)-UTM(6*I-0))/DT)**2

!Calculate final UTP
	UTP(6*I-5)=(R(6*I-5)+FRX(I)+WSX(I))/((MFIL(I)/DT**2)+VDP(I)/(2.*DT))
	UTP(6*I-4)=(R(6*I-4)+FRY(I)+BOYY(I)+WSY(I))/((MFIL(I)/DT**2)+VDP(I)/(2.*DT))
	UTP(6*I-3)=(R(6*I-3)+BOYZ(I)+FRZ(I))/((MFIL(I)/DT**2)+VDP(I)/(2.*DT))
	UTP(6*I-2)=R(6*I-2)/((MFIL(I)*JS/MN)/DT**2+VDP(I)/(2.*DT))
	UTP(6*I-1)=R(6*I-1)/((MFIL(I)*JS/MN)/DT**2+VDP(I)/(2.*DT))
	UTP(6*I-0)=R(6*I-0)/((MFIL(I)*JS/MN)/DT**2+VDP(I)/(2.*DT))

        !Check that particles haven't left the domain
        !and freeze them if they have!
         ! XIND = INT((NRXF(1,I) + UTP(6*I-5))/GRID)
         ! YIND = INT((NRXF(2,I) + UTP(6*I-4))/GRID)
         ! IF(XIND > 2000 .OR. XIND < -100 .OR. YIND > 2000 .OR. YIND < -100 .OR. &
         !     ABS(UTP(6*I-5)) > MAXUT .OR. ABS(UTP(6*I-4)) > MAXUT .OR. &
         !     ABS(UTP(6*I-3)) > MAXUT) THEN
         !   UTP(6*I-5) = UT(6*I-5)
         !   UTP(6*I-4) = UT(6*I-4)
         !   UTP(6*I-3) = UT(6*I-3)
         !   UTP(6*I-2) = UT(6*I-2)
         !   UTP(6*I-1) = UT(6*I-1)
         !   UTP(6*I-0) = UT(6*I-0)
         !   PRINT *, myid, " Lost a particle : ",I," at time: ",T
         ! END IF

!	IF (X.LT.4000.0.AND.Y.LT.4000.0) THEN
!	IF ((ZB.GT.WL-2.0*SCL.OR.Z-ZB.LT.5.0*SCL).AND.ABS(Z-ZB).LT.2.0*SCL) THEN
!Freeze particles near bed if friction is too high
        IF (ABS(ZB-Z).LT.SCL*2.0.AND.(VDP(I).GE.1.0e+08*SCL*SCL.OR.T.LT.20.0)) THEN
        UTP(6*I-5)=UT(6*I-5)
        UTP(6*I-4)=UT(6*I-4)
        UTP(6*I-3)=UT(6*I-3)
        UTP(6*I-2)=UT(6*I-2)
        UTP(6*I-1)=UT(6*I-1)
        UTP(6*I-0)=UT(6*I-0)
	ENDIF

!	IF (X.LT.250.0.OR.Y.LT.250.0) THEN
        !Backplane is fixed
	IF (Y.LT.250.0) THEN
        UTP(6*I-5)=UT(6*I-5)
        UTP(6*I-4)=UT(6*I-4)
	ENDIF

        !No XY displacement of some nodes? 
        IF (MOD(myid+1,ntasks/YN).eq.0) THEN
        IF (X.GT.MAXX-250.0) THEN
        UTP(6*I-5)=UT(6*I-5)
        UTP(6*I-4)=UT(6*I-4)
	ENDIF
	ENDIF


       !Compute damping
	DMPEN=DMPEN+VDP(I)*(UTP(6*I-5)-UTM(6*I-5))**2/(4*DT)
	DMPEN=DMPEN+VDP(I)*(UTP(6*I-4)-UTM(6*I-4))**2/(4*DT)
	DMPEN=DMPEN+VDP(I)*(UTP(6*I-3)-UTM(6*I-3))**2/(4*DT)
	DMPEN=DMPEN+VDP(I)*(UTP(6*I-2)-UTM(6*I-2))**2/(4*DT)
	DMPEN=DMPEN+VDP(I)*(UTP(6*I-1)-UTM(6*I-1))**2/(4*DT)
	DMPEN=DMPEN+VDP(I)*(UTP(6*I-0)-UTM(6*I-0))**2/(4*DT)
 

       !Calculate all the energy and sum over partitions

        WEN=WEN+WE(I)

        ENM=ENM+EN(6*I-5)+EN(6*I-4)+EN(6*I-3)+EN(6*I-2)+EN(6*I-1)+EN(6*I-0)

        EMM(I)=EN(6*I-5)+EN(6*I-4)+EN(6*I-3)+EN(6*I-2)+EN(6*I-1)+EN(6*I-0)

	IF (Z.LT.WL+Y*SSB) THEN
	MGH=MGH+1.06*MFIL(I)*ABS(Z-WL-Y*SSB)/SQB
	ELSE
	MGH=MGH+9.8*MFIL(I)*ABS(Z-WL-Y*SSB)/SQB
	ENDIF

        !TODO  TIME 9
       CALL CPU_TIME(TT(9))
       TTCUM(9) = TTCUM(9) + (TT(9) - TT(8))

 27	CONTINUE

        IF (RY.EQ.1) THEN
        CALL MPI_ALLREDUCE(GSUM,GSUM0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(MGH,MGH0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        END IF

!	IF (RY.EQ.REST*STEPS0+1) THEN

	IF (RY.EQ.RY0) THEN
        CALL MPI_ALLREDUCE(ENM0,ENMS0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	ENDIF

 	IF (MOD(RY,100).EQ.0) THEN
        CALL MPI_ALLREDUCE(KIN,KINS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(KIN2,KINS2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(ENM,ENMS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(MGH,MGHS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(DMPEN,DMPENS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(WEN,WENS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(GSUM,GSUMS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(DPE,DPES,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(BCE,BCES,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(BCC,BCCS,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(PSUM,PSUMS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 	IF (myid.EQ.0) WRITE(612,*) T,WENS+ENMS+ENMS0-BCES+KINS+KINS2&
      	+MGHS-MGH0,PSUMS-DPES-DMPENS-GSUMS+GSUM0

 	IF (myid.EQ.0) WRITE(610,10) T,WENS,ENMS+ENMS0,KINS,MGHS-MGH0
 	IF (myid.EQ.0) WRITE(611,10) T,DPES,DMPENS,PSUMS,GSUMS
 	IF (myid.EQ.0) WRITE(613,19) T,KINS2,BCES,BCCS
	END IF

	T=T+DT
	
	MML=0.0

        !TODO  TIME 10
       CALL CPU_TIME(TT(10))
       TTCUM(10) = TTCUM(10) + (TT(10) - TT(9))

!Check for fracture!

!internally
        DO I=1,NTOT
	N1=NAN(1,I)
	N2=NAN(2,I)
        IF (EF(I).GT.0.0) THEN
	DDX=NRXF(1,N1)+UT(6*N1-5)-NRXF(1,N2)-UT(6*N2-5)
	DDY=NRXF(2,N1)+UT(6*N1-4)-NRXF(2,N2)-UT(6*N2-4)
	DDZ=NRXF(3,N1)+UT(6*N1-3)-NRXF(3,N2)-UT(6*N2-3)
	DX=NRXF(1,N1)-NRXF(1,N2)
	DY=NRXF(2,N1)-NRXF(2,N2)
	DZ=NRXF(3,N1)-NRXF(3,N2)
	L=SQRT(DX**2+DY**2+DZ**2)
	DL=SQRT(DDX**2+DDY**2+DDZ**2)
	V1=ABS(UT(6*N1-2)-UT(6*N2-2))
	V2=ABS(UT(6*N1-1)-UT(6*N2-1))
	V3=ABS(UT(6*N1-0)-UT(6*N2-0))
	MAXV=MAX(V1,V2,V3)
	LOAD=((DL-L)+0.05*MAXV)
	IF (LOAD.GT.MML) MML=LOAD
	 IF (LOAD.GT.MLOAD.AND.T.GT.fractime) THEN
	 BCE=BCE+0.5*EF(I)*S**2/LNN*(DL-L)**2
         EF(I)=0.0
	 BCC=BCC+1
	 ENDIF
        ENDIF
 	END DO

!        IF (myid.ne.0.and.myid.ne.(ntasks-1)/2+1) THEN
!        DO I=1,NTOL
!	N1=NANL(1,I)
!	N2=NANL(2,I)
!        IF (EFL(I).GT.0.0) THEN
!	DDX=NRXFL(1,N1)+UTL(6*N1-5)-NRXF(1,N2)-UT(6*N2-5)
!	DDY=NRXFL(2,N1)+UTL(6*N1-4)-NRXF(2,N2)-UT(6*N2-4)
!	DDZ=NRXFL(3,N1)+UTL(6*N1-3)-NRXF(3,N2)-UT(6*N2-3)
!	DX=NRXFL(1,N1)-NRXF(1,N2)
!	DY=NRXFL(2,N1)-NRXF(2,N2)
!	DZ=NRXFL(3,N1)-NRXF(3,N2)
!	L=SQRT(DX**2+DY**2+DZ**2)
!	DL=SQRT(DDX**2+DDY**2+DDZ**2)
!	V1=ABS(UTL(6*N1-2)-UT(6*N2-2))
!	V2=ABS(UTL(6*N1-1)-UT(6*N2-1))
!	V3=ABS(UTL(6*N1-0)-UT(6*N2-0))
!	MAXV=MAX(V1,V2,V3)
!	LOAD=((DL-L)+0.05*MAXV)
!	IF (LOAD.GT.MML) MML=LOAD
!	 IF (LOAD.GT.MLOAD.AND.T.GT.1.0) THEN
!	 BCE=BCE+0.5*EFL(I)*S**2/LNN*(DL-L)**2
!         EFL(I)=0.0
!	 BCC=BCC+1
!	 ENDIF
!	ENDIF
! 	END DO
!	ENDIF

        IF (myid.ne.ntasks-1.and.myid.ne.(ntasks-1)/2) THEN
        DO I=1,NTOR
	N1=NANR(1,I)
	N2=NANR(2,I)
        IF (EFR(I).GT.0.0) THEN
	DDX=NRXFR(1,N1)+UTR(6*N1-5)-NRXF(1,N2)-UT(6*N2-5)
	DDY=NRXFR(2,N1)+UTR(6*N1-4)-NRXF(2,N2)-UT(6*N2-4)
	DDZ=NRXFR(3,N1)+UTR(6*N1-3)-NRXF(3,N2)-UT(6*N2-3)
	DX=NRXFR(1,N1)-NRXF(1,N2)
	DY=NRXFR(2,N1)-NRXF(2,N2)
	DZ=NRXFR(3,N1)-NRXF(3,N2)
	L=SQRT(DX**2+DY**2+DZ**2)
	DL=SQRT(DDX**2+DDY**2+DDZ**2)
	V1=ABS(UTR(6*N1-2)-UT(6*N2-2))
	V2=ABS(UTR(6*N1-1)-UT(6*N2-1))
	V3=ABS(UTR(6*N1-0)-UT(6*N2-0))
	MAXV=MAX(V1,V2,V3)
	LOAD=((DL-L)+0.05*MAXV)
	IF (LOAD.GT.MML) MML=LOAD
	 IF (LOAD.GT.MLOAD.AND.T.GT.40.0) THEN
	 BCE=BCE+0.5*EFR(I)*S**2/LNN*(DL-L)**2
         EFR(I)=0.0
	 BCC=BCC+1
	 ENDIF
	ENDIF
 	END DO
	ENDIF

!across the boundaries...
        dest=myid+1
        source=myid-1
        tag=131
        IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
        CALL MPI_Send(EFR,NTOR,MPI_DOUBLE_PRECISION,&
        dest,tag,MPI_COMM_WORLD,ierr)
        IF (MOD(myid,ntasks/YN).ne.0)&
        CALL MPI_Recv(EFL,NTOL,MPI_DOUBLE_PRECISION,&
        source,tag,MPI_COMM_WORLD,stat,ierr)

        IF (myid.le.(YN-1)*ntasks/YN) THEN
        DO I=1,NTOF
	N1=NANF(1,I)
	N2=NANF(2,I)
        IF (EFF(I).GT.0.0) THEN
	DDX=NRXFF(1,N1)+UTF(6*N1-5)-NRXF(1,N2)-UT(6*N2-5)
	DDY=NRXFF(2,N1)+UTF(6*N1-4)-NRXF(2,N2)-UT(6*N2-4)
	DDZ=NRXFF(3,N1)+UTF(6*N1-3)-NRXF(3,N2)-UT(6*N2-3)
	DX=NRXFF(1,N1)-NRXF(1,N2)
	DY=NRXFF(2,N1)-NRXF(2,N2)
	DZ=NRXFF(3,N1)-NRXF(3,N2)
	L=SQRT(DX**2+DY**2+DZ**2)
	DL=SQRT(DDX**2+DDY**2+DDZ**2)
	V1=ABS(UTF(6*N1-2)-UT(6*N2-2))
	V2=ABS(UTF(6*N1-1)-UT(6*N2-1))
	V3=ABS(UTF(6*N1-0)-UT(6*N2-0))
	MAXV=MAX(V1,V2,V3)
	LOAD=((DL-L)+0.05*MAXV)
	IF (LOAD.GT.MML) MML=LOAD
	 IF (LOAD.GT.MLOAD.AND.T.GT.40.0) THEN
	 BCE=BCE+0.5*EFF(I)*S**2/LNN*(DL-L)**2
         EFF(I)=0.0
	 BCC=BCC+1
	 ENDIF
	ENDIF
 	END DO
	ENDIF

        dest=myid+ntasks/YN
        source=myid-ntasks/YN
        tag=132
        IF (myid.lt.(YN-1)*ntasks/YN)&
        CALL MPI_Send(EFF,NTOF,MPI_DOUBLE_PRECISION,&
        dest,tag,MPI_COMM_WORLD,ierr)
        IF (myid.ge.ntasks/YN)&
        CALL MPI_Recv(EFB,NTOB,MPI_DOUBLE_PRECISION,&
        source,tag,MPI_COMM_WORLD,stat,ierr)

        IF (myid.le.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
        DO I=1,NTOFL
	N1=NANFL(1,I)
	N2=NANFL(2,I)
        IF (EFFL(I).GT.0.0) THEN
	DDX=NRXFFL(1,N1)+UTFL(6*N1-5)-NRXF(1,N2)-UT(6*N2-5)
	DDY=NRXFFL(2,N1)+UTFL(6*N1-4)-NRXF(2,N2)-UT(6*N2-4)
	DDZ=NRXFFL(3,N1)+UTFL(6*N1-3)-NRXF(3,N2)-UT(6*N2-3)
	DX=NRXFFL(1,N1)-NRXF(1,N2)
	DY=NRXFFL(2,N1)-NRXF(2,N2)
	DZ=NRXFFL(3,N1)-NRXF(3,N2)
	L=SQRT(DX**2+DY**2+DZ**2)
	DL=SQRT(DDX**2+DDY**2+DDZ**2)
	V1=ABS(UTFL(6*N1-2)-UT(6*N2-2))
	V2=ABS(UTFL(6*N1-1)-UT(6*N2-1))
	V3=ABS(UTFL(6*N1-0)-UT(6*N2-0))
	MAXV=MAX(V1,V2,V3)
	LOAD=((DL-L)+0.05*MAXV)
	IF (LOAD.GT.MML) MML=LOAD
	 IF (LOAD.GT.MLOAD.AND.T.GT.40.0) THEN
	 BCE=BCE+0.5*EFFL(I)*S**2/LNN*(DL-L)**2
         EFFL(I)=0.0
	 BCC=BCC+1
	 ENDIF
	ENDIF
 	END DO
	ENDIF

        dest=myid+ntasks/YN-1
        source=myid-ntasks/YN+1
        tag=133
        IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
        CALL MPI_Send(EFFL,NTOFL,MPI_DOUBLE_PRECISION,&
        dest,tag,MPI_COMM_WORLD,ierr)
        IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
        CALL MPI_Recv(EFBR,NTOBR,MPI_DOUBLE_PRECISION,& 
        source,tag,MPI_COMM_WORLD,stat,ierr)

        IF (myid.le.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
        DO I=1,NTOFR
	N1=NANFR(1,I)
	N2=NANFR(2,I)
        IF (EFFR(I).GT.0.0) THEN
	DDX=NRXFFR(1,N1)+UTFR(6*N1-5)-NRXF(1,N2)-UT(6*N2-5)
	DDY=NRXFFR(2,N1)+UTFR(6*N1-4)-NRXF(2,N2)-UT(6*N2-4)
	DDZ=NRXFFR(3,N1)+UTFR(6*N1-3)-NRXF(3,N2)-UT(6*N2-3)
	DX=NRXFFR(1,N1)-NRXF(1,N2)
	DY=NRXFFR(2,N1)-NRXF(2,N2)
	DZ=NRXFFR(3,N1)-NRXF(3,N2)
	L=SQRT(DX**2+DY**2+DZ**2)
	DL=SQRT(DDX**2+DDY**2+DDZ**2)
	V1=ABS(UTFR(6*N1-2)-UT(6*N2-2))
	V2=ABS(UTFR(6*N1-1)-UT(6*N2-1))
	V3=ABS(UTFR(6*N1-0)-UT(6*N2-0))
	MAXV=MAX(V1,V2,V3)
	LOAD=((DL-L)+0.05*MAXV)
	IF (LOAD.GT.MML) MML=LOAD
	 IF (LOAD.GT.MLOAD.AND.T.GT.40.0) THEN
	 BCE=BCE+0.5*EFFR(I)*S**2/LNN*(DL-L)**2
         EFFR(I)=0.0
	 BCC=BCC+1
	 ENDIF
	ENDIF
 	END DO
	ENDIF

        dest=myid+ntasks/YN+1
        source=myid-ntasks/YN-1
        tag=134
        IF (myid.lt.(YN-1)*ntasks/YN&
        .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
        CALL MPI_Send(EFFR,NTOFR,MPI_DOUBLE_PRECISION,&
        dest,tag,MPI_COMM_WORLD,ierr)
        IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
        CALL MPI_Recv(EFBL,NTOBL,MPI_DOUBLE_PRECISION,&
        source,tag,MPI_COMM_WORLD,stat,ierr)

!        IF (myid.gt.(ntasks-1)/YN) THEN
!        DO I=1,NTOB
!	N1=NANB(1,I)
!	N2=NANB(2,I)
!        IF (EFB(I).GT.0.0) THEN
!	DDX=NRXFB(1,N1)+UTB(6*N1-5)-NRXF(1,N2)-UT(6*N2-5)
!	DDY=NRXFB(2,N1)+UTB(6*N1-4)-NRXF(2,N2)-UT(6*N2-4)
!	DDZ=NRXFB(3,N1)+UTB(6*N1-3)-NRXF(3,N2)-UT(6*N2-3)
!	DX=NRXFB(1,N1)-NRXF(1,N2)
!	DY=NRXFB(2,N1)-NRXF(2,N2)
!	DZ=NRXFB(3,N1)-NRXF(3,N2)
!	L=SQRT(DX**2+DY**2+DZ**2)
!	DL=SQRT(DDX**2+DDY**2+DDZ**2)
!	V1=ABS(UTF(6*N1-2)-UT(6*N2-2))
!	V2=ABS(UTF(6*N1-1)-UT(6*N2-1))
!	V3=ABS(UTF(6*N1-0)-UT(6*N2-0))
!	MAXV=MAX(V1,V2,V3)
!	LOAD=((DL-L)+0.05*MAXV)
!	IF (LOAD.GT.MML) MML=LOAD
!	 IF (LOAD.GT.MLOAD.AND.T.GT.1.0) THEN
!	 BCE=BCE+0.5*EFB(I)*S**2/LNN*(DL-L)**2
!         EFB(I)=0.0
!	 BCC=BCC+1
!	 ENDIF
!	ENDIF
! 	END DO
!	ENDIF

!        IF (myid.gt.(ntasks-1)/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
!        DO I=1,NTOBL
!	N1=NANBL(1,I)
!	N2=NANBL(2,I)
!        IF (EFBL(I).GT.0.0) THEN
!	DDX=NRXFBL(1,N1)+UTBL(6*N1-5)-NRXF(1,N2)-UT(6*N2-5)
!	DDY=NRXFBL(2,N1)+UTBL(6*N1-4)-NRXF(2,N2)-UT(6*N2-4)
!	DDZ=NRXFBL(3,N1)+UTBL(6*N1-3)-NRXF(3,N2)-UT(6*N2-3)
!	DX=NRXFBL(1,N1)-NRXF(1,N2)
!	DY=NRXFBL(2,N1)-NRXF(2,N2)
!	DZ=NRXFBL(3,N1)-NRXF(3,N2)
!	L=SQRT(DX**2+DY**2+DZ**2)
!	DL=SQRT(DDX**2+DDY**2+DDZ**2)
!	V1=ABS(UTFL(6*N1-2)-UT(6*N2-2))
!	V2=ABS(UTFL(6*N1-1)-UT(6*N2-1))
!	V3=ABS(UTFL(6*N1-0)-UT(6*N2-0))
!	MAXV=MAX(V1,V2,V3)
!	LOAD=((DL-L)+0.05*MAXV)
!	IF (LOAD.GT.MML) MML=LOAD
!	 IF (LOAD.GT.MLOAD.AND.T.GT.1.0) THEN
!	 BCE=BCE+0.5*EFBL(I)*S**2/LNN*(DL-L)**2
!         EFBL(I)=0.0
!	 BCC=BCC+1
!	 ENDIF
!	ENDIF
! 	END DO
!	ENDIF

!        IF (myid.gt.(ntasks-1)/YN
!     1	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
!        DO I=1,NTOBR
!	N1=NANBR(1,I)
!	N2=NANBR(2,I)
!        IF (EFBR(I).GT.0.0) THEN
!	DDX=NRXFBR(1,N1)+UTBR(6*N1-5)-NRXF(1,N2)-UT(6*N2-5)
!	DDY=NRXFBR(2,N1)+UTBR(6*N1-4)-NRXF(2,N2)-UT(6*N2-4)
!	DDZ=NRXFBR(3,N1)+UTBR(6*N1-3)-NRXF(3,N2)-UT(6*N2-3)
!	DX=NRXFBR(1,N1)-NRXF(1,N2)
!	DY=NRXFBR(2,N1)-NRXF(2,N2)
!	DZ=NRXFBR(3,N1)-NRXF(3,N2)
!	L=SQRT(DX**2+DY**2+DZ**2)
!	DL=SQRT(DDX**2+DDY**2+DDZ**2)
!	V1=ABS(UTFR(6*N1-2)-UT(6*N2-2))
!	V2=ABS(UTFR(6*N1-1)-UT(6*N2-1))
!	V3=ABS(UTFR(6*N1-0)-UT(6*N2-0))
!	MAXV=MAX(V1,V2,V3)
!	LOAD=((DL-L)+0.05*MAXV)
!	IF (LOAD.GT.MML) MML=LOAD
!	 IF (LOAD.GT.MLOAD.AND.T.GT.1.0) THEN
!	 BCE=BCE+0.5*EFBR(I)*S**2/LNN*(DL-L)**2
!         EFBR(I)=0.0
!	 BCC=BCC+1
!	 ENDIF
!	ENDIF
! 	END DO
!	ENDIF

	DO I=1,6*NN
	UTM(I)=UT(I)
	UT(I)=UTP(I)
	END DO

 !output
	IF (MOD(RY,OUTINT).EQ.1) THEN

20 FORMAT(" Times: ",11F4.1)
          WRITE(*,20) TTCUM

          dest=0

          tag=151
          IF (myid.NE.0)&
          CALL MPI_Send(UT,6*NN,MPI_DOUBLE_PRECISION,&
          dest,tag,MPI_COMM_WORLD,ierr)

          tag=152
          IF (myid.NE.0)&
          CALL MPI_Send(NAN,3*NTOT,MPI_INTEGER,&
          dest,tag,MPI_COMM_WORLD,ierr)

          tag=161
          IF (myid.NE.0)&
          CALL MPI_Send(NRXF,3*NN,MPI_DOUBLE_PRECISION,&
          dest,tag,MPI_COMM_WORLD,ierr)

          IF (myid.EQ.0) THEN
  	  NRY=INT(RY/OUTINT)
          OPEN(UNIT=910,FILE=TRIM(resdir)//'/'//TRIM(runname)//'_JYR'//na(NRY)//'.csv',STATUS='UNKNOWN')
          OPEN(UNIT=920,FILE=TRIM(resdir)//'/'//TRIM(runname)//'_STR'//na(NRY)//'.csv',STATUS='UNKNOWN')

          DO KK=1,ntasks-1

          source=KK

          tag=151
          CALL MPI_Recv(UTW,6*PNN(source),MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)

          tag=152
          CALL MPI_Recv(NANW,3*NTOTW(source),MPI_INTEGER,source,tag,MPI_COMM_WORLD,stat,ierr)

          tag=161
          CALL MPI_Recv(NRXFW,3*PNN(source),MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)

          DO I=1,PNN(source)
	  X=NRXFW(1,I)+UTW(6*I-5)
	  Y=NRXFW(2,I)+UTW(6*I-4)
	  Z=NRXFW(3,I)+UTW(6*I-3)
          WRITE(910,12) X,Y,Z
          END DO

          DO I=1,NTOTW(source)

	  X=(NRXFW(1,NANW(1,I))+UTW(6*NANW(1,I)-5)+NRXFW(1,NANW(2,I))+UTW(6*NANW(2,I)-5))/2.0

	  Y=(NRXFW(2,NANW(1,I))+UTW(6*NANW(1,I)-4)+NRXFW(2,NANW(2,I))+UTW(6*NANW(2,I)-4))/2.0

	  Z=(NRXFW(3,NANW(1,I))+UTW(6*NANW(1,I)-3)+NRXFW(3,NANW(2,I))+UTW(6*NANW(2,I)-3))/2.0

	  
	N1=NANW(1,I)
	N2=NANW(2,I)
	DDX=NRXFW(1,N1)+UTW(6*N1-5)-NRXFW(1,N2)-UTW(6*N2-5)
	DDY=NRXFW(2,N1)+UTW(6*N1-4)-NRXFW(2,N2)-UTW(6*N2-4)
	DDZ=NRXFW(3,N1)+UTW(6*N1-3)-NRXFW(3,N2)-UTW(6*N2-3)
	DX=NRXFW(1,N1)-NRXFW(1,N2)
	DY=NRXFW(2,N1)-NRXFW(2,N2)
	DZ=NRXFW(3,N1)-NRXFW(3,N2)
	L=SQRT(DX**2+DY**2+DZ**2)
	DL=SQRT(DDX**2+DDY**2+DDZ**2)
	STR=(DL-L)/L


!          IF (Z.GT.MAXZ-120.0) WRITE(920,12) X,Y,Z,STR
          WRITE(920,12) X,Y,Z,STR
          END DO

          END DO

          DO I=1,NN
	  X=NRXF(1,I)+UT(6*I-5)
	  Y=NRXF(2,I)+UT(6*I-4)
	  Z=NRXF(3,I)+UT(6*I-3)
          WRITE(910,12) X,Y,Z
          END DO

          DO I=1,NTOT
	  X=(NRXF(1,NAN(1,I))+UT(6*NAN(1,I)-5)+NRXF(1,NAN(2,I))+UT(6*NAN(2,I)-5))/2.0
	  Y=(NRXF(2,NAN(1,I))+UT(6*NAN(1,I)-4)+NRXF(2,NAN(2,I))+UT(6*NAN(2,I)-4))/2.0
	  Z=(NRXF(3,NAN(1,I))+UT(6*NAN(1,I)-3)+NRXF(3,NAN(2,I))+UT(6*NAN(2,I)-3))/2.0
	  
	N1=NAN(1,I)
	N2=NAN(2,I)
	DDX=NRXF(1,N1)+UT(6*N1-5)-NRXF(1,N2)-UT(6*N2-5)
	DDY=NRXF(2,N1)+UT(6*N1-4)-NRXF(2,N2)-UT(6*N2-4)
	DDZ=NRXF(3,N1)+UT(6*N1-3)-NRXF(3,N2)-UT(6*N2-3)
	DX=NRXF(1,N1)-NRXF(1,N2)
	DY=NRXF(2,N1)-NRXF(2,N2)
	DZ=NRXF(3,N1)-NRXF(3,N2)
	L=SQRT(DX**2+DY**2+DZ**2)
	DL=SQRT(DDX**2+DDY**2+DDZ**2)
	STR=(DL-L)/L

!          IF (Z.GT.MAXZ-120.0) WRITE(920,12) X,Y,STR
          WRITE(920,12) X,Y,Z,STR
          END DO

          ENDIF

!          CALL PSNET(NTOT,NN,myid,GL,WL,SUB,ntasks)
	  CLOSE(910)
	  CLOSE(920)
	ENDIF

	IF (MOD(RY,RESOUTINT).EQ.1) THEN
            !Add restart stuff here
            CALL WriteRestart()
        END IF


!	IF (RY.EQ.RY0+STEPS0) THEN
!        CALL CRACK(NN,NTOT,CCN,EF,NAN,NRXF,UT,WL,myid,ntasks,NANR,NTOR,EFR)
!
!        dest=0
!	tag=171
!        IF (myid.NE.0) CALL MPI_Send(CCN,NN,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,ierr)
!
!        IF (myid.EQ.0) THEN
!	OPEN(UNIT=112, FILE='frag',STATUS='UNKNOWN')
!         DO KK=1,ntasks-1
!         source=KK
!         tag=171
!         CALL MPI_Recv(CCNW,NN,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
!           DO I=1,NN
!           CCN(I)=CCN(I)+CCNW(I)
!	   END DO
!	 END DO
!          DO I=1,NN
!          IF (CCN(I).NE.0) WRITE(112,*) I,CCN(I)
!	  CCN(I)=0
!	  END DO
!	CLOSE(112)
!	ENDIF
!	ENDIF

        IF (MOD(RY,1000).EQ.0) THEN
	IF (myid.EQ.0) WRITE(*,*) RY,MML/MLOAD,BCC
!	IF (BCC.NE.0) WRITE(700+myid,16) T, BCC
	BCC=0
	ENDIF

        !TODO  TIME 11
       CALL CPU_TIME(TT(11))
       TTCUM(11) = TTCUM(11) + (TT(11) - TT(10))

 100	CONTINUE !end of time loop

        OPEN(UNIT=930,FILE=TRIM(resdir)//'/RCK',STATUS='UNKNOWN')
        DO KK=0,ntasks-1
        WRITE(930,*) KK, NTOTW(KK),PNN(KK),PTR(KK),PTB(KK)
        END DO
        CLOSE(930)

        CALL CPU_TIME(T2)
	WRITE(*,*) 'TIME=',T2-T1,TS1

        CALL WriteRestart()

!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        IF (myid.EQ.0) CALL CRACK(ntasks,NTOTW,PNN,PTR,PTB,YN)

        CALL MPI_FINALIZE(rc)
  	STOP

CONTAINS

 SUBROUTINE WriteRestart()
   ! Write out the restart files
	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(runname)//'_REST0'//na(myid),STATUS='UNKNOWN')
	WRITE(117+myid,*) NN,NTOT,NTOL,NTOR,NTOF,NTOB,BCC
	write(117+myid,*) NTOFL,NTOFR,NTOBL,NTOBR
	WRITE(117+myid,*) MAXX,MAXY,MAXZ,DMPEN,ENM+ENM0
	WRITE(117+myid,*) DPE,BCE,MGH0,GSUM0,PSUM,T,RY-1
	WRITE(117+myid,*) XN,YN
	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(runname)//'_REST1'//na(myid),STATUS='UNKNOWN')
	DO I=1,NTOT+NTOL+NTOR+NTOF+NTOB+NTOFL+NTOFR+NTOBL+NTOBR
	WRITE(117+myid,*) CT(12*I-11),CT(12*I-10),CT(12*I-9),CT(12*I-8),CT(12*I-7),CT(12*I-6)
	WRITE(117+myid,*) CT(12*I-5),CT(12*I-4),CT(12*I-3),CT(12*I-2),CT(12*I-1),CT(12*I-0)
	END DO
	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(runname)//'_REST2'//na(myid),STATUS='UNKNOWN')
	DO I=1,NN
	WRITE(117+myid,*) UT(6*I-5),UT(6*I-4),UT(6*I-3),UT(6*I-2),UT(6*I-1),UT(6*I-0)
	WRITE(117+myid,*) UTM(6*I-5),UTM(6*I-4),UTM(6*I-3),UTM(6*I-2),UTM(6*I-1),UTM(6*I-0)
	END DO
	CLOSE (117+myid)

 !NTOT - no connections
 !NAN - the particle numbers (1 and 2), 2*ntot array
 !NRXF - the original locations of all the particles
 !EF - Youngs modulus
	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FS'//na(myid),STATUS='UNKNOWN')
	DO SI=1,NTOT
	WRITE(117+myid,*) NAN(1,SI),NAN(2,SI),NRXF(1,NAN(1,SI)),&
      	   NRXF(2,NAN(1,SI)),NRXF(3,NAN(1,SI)),NRXF(1,NAN(2,SI)),&
           NRXF(2,NAN(2,SI)),NRXF(3,NAN(2,SI)),EF(SI)
	END DO
	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSL'//na(myid),STATUS='UNKNOWN')
	DO SI=1,NTOL
	WRITE(117+myid,*) NANL(1,SI),NANL(2,SI),NRXFL(1,NANL(1,SI)),&
      	   NRXFL(2,NANL(1,SI)),NRXFL(3,NANL(1,SI)),NRXF(1,NANL(2,SI)),&
           NRXF(2,NANL(2,SI)),NRXF(3,NANL(2,SI)),EFL(SI)
	END DO
	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSR'//na(myid),STATUS='UNKNOWN')
	DO SI=1,NTOR
	WRITE(117+myid,*) NANR(1,SI),NANR(2,SI),NRXFR(1,NANR(1,SI)),&
      	   NRXFR(2,NANR(1,SI)),NRXFR(3,NANR(1,SI)),NRXF(1,NANR(2,SI)),&
           NRXF(2,NANR(2,SI)),NRXF(3,NANR(2,SI)),EFR(SI)
	END DO
	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSF'//na(myid),STATUS='UNKNOWN')
	DO SI=1,NTOF
	WRITE(117+myid,*) NANF(1,SI),NANF(2,SI),NRXFF(1,NANF(1,SI)),&
      	   NRXFF(2,NANF(1,SI)),NRXFF(3,NANF(1,SI)),NRXF(1,NANF(2,SI)),&
           NRXF(2,NANF(2,SI)),NRXF(3,NANF(2,SI)),EFF(SI)
	END DO
	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSFL'//na(myid),STATUS='UNKNOWN')
	DO SI=1,NTOFL
	WRITE(117+myid,*) NANFL(1,SI),NANFL(2,SI),NRXFFL(1,NANFL(1,SI)),&
      	 NRXFFL(2,NANFL(1,SI)),NRXFFL(3,NANFL(1,SI)),NRXF(1,NANFL(2,SI)),&
         NRXF(2,NANFL(2,SI)),NRXF(3,NANFL(2,SI)),EFFL(SI)
	END DO
	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSFR'//na(myid),STATUS='UNKNOWN')
	DO SI=1,NTOFR
	WRITE(117+myid,*) NANFR(1,SI),NANFR(2,SI),NRXFFR(1,NANFR(1,SI)),&
      	 NRXFFR(2,NANFR(1,SI)),NRXFFR(3,NANFR(1,SI)),NRXF(1,NANFR(2,SI)),&
         NRXF(2,NANFR(2,SI)),NRXF(3,NANFR(2,SI)),EFFR(SI)
	END DO
	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSB'//na(myid),STATUS='UNKNOWN')
	DO SI=1,NTOB
	WRITE(117+myid,*) NANB(1,SI),NANB(2,SI),NRXFB(1,NANB(1,SI)),&
      	   NRXFB(2,NANB(1,SI)),NRXFB(3,NANB(1,SI)),NRXF(1,NANB(2,SI)),&
           NRXF(2,NANB(2,SI)),NRXF(3,NANB(2,SI)),EFB(SI)
	END DO
	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSBL'//na(myid),STATUS='UNKNOWN')
	DO SI=1,NTOBL
	WRITE(117+myid,*) NANBL(1,SI),NANBL(2,SI),NRXFBL(1,NANBL(1,SI)),&
      	 NRXFBL(2,NANBL(1,SI)),NRXFBL(3,NANBL(1,SI)),NRXF(1,NANBL(2,SI)),&
         NRXF(2,NANBL(2,SI)),NRXF(3,NANBL(2,SI)),EFBL(SI)
	END DO
	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FSBR'//na(myid),STATUS='UNKNOWN')
	DO SI=1,NTOBR
	WRITE(117+myid,*) NANBR(1,SI),NANBR(2,SI),NRXFBR(1,NANBR(1,SI)),&
      	 NRXFBR(2,NANBR(1,SI)),NRXFBR(3,NANBR(1,SI)),NRXF(1,NANBR(2,SI)),&
         NRXF(2,NANBR(2,SI)),NRXF(3,NANBR(2,SI)),EFBR(SI)
	END DO
	CLOSE (117+myid)

	IF (myid.EQ.0) THEN
	OPEN(UNIT=117,FILE=TRIM(resdir)//'/crkd',STATUS='UNKNOWN')
	DO SI=1,ntasks
	WRITE(117,*) NTOTW(SI-1),PNN(SI-1),PTR(SI-1),PTB(SI-1)
	END DO
	CLOSE (117)
	ENDIF

 END SUBROUTINE WriteRestart
END PROGRAM

