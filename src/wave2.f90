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
        USE TypeDefs
        USE Lattice
        USE Utils

	IMPLICIT NONE
        INCLUDE 'mpif.h'
        INCLUDE 'na90.dat'
	REAL*8 EMM(NOMA),BOYZ(NOMA),BOYY(NOMA),WSY(NOMA)
	REAL*8 EN(NODM),WE(NOMA),VDP(NOMA),WSX(NOMA),BED(-100:2000,-100:2000)
	REAL*8 SUF(-100:2000,-100:2000)
	REAL*8, ALLOCATABLE :: MFIL(:),EFC(:),EFS(:)
	REAL*8 FBED(-100:2000,-100:2000)
	REAL*8 DIX,DIY,DIZ,FRIC,UC
	REAL*8 UTP(NODM),R(NODM)
	REAL*8 UTW(NODM),NRXFW(3,NOMA)
	REAL*8 ENM0,ENMS0,POR,GRID
        REAL*8 RHO,RHOW,GRAV, BedIntConst
        REAL*8 CT(NODC),FRX(NOMA),FRY(NOMA),FRZ(NOMA),BBox(6)
        REAL*8, ALLOCATABLE :: PBBox(:,:)
	INTEGER CCN(NOMA),CCNW(NOMA),CB(NOMA),CCB(NOMA,3)
	INTEGER NANW(3,NOCON),mask,GEOMMASK(-100:2000,-100:2000)
	INTEGER REST,RY0,NM2,noprocs, counter
	INTEGER PTB(0:5000),PTR(0:5000)
	REAL*8 M,MN,JS,DT,T,X,Y,E,GSUM,GSUM0
	REAL*8 DMPEN,PSUM,KIN,KIN2,PRESS,MELT
	INTEGER I,N,NL,NN,STEPS,IX,IM,MS,N1,N2,P1,RY
	INTEGER PN,NRY,PNN(0:5000),NTOTW(0:5000),XK,YK,ZK
        REAL*8 L,ALF,MLOAD,DMP,VEL,G,S1,S2,M1,B1,B2
	REAL*8 S,LOAD,DMP2,BCE,STR,I1,I2,ZB,ZS
	REAL*8 T2,T1,TS1,TT1,TT2,ENM,WEN,ERSUM
	REAL*8 TT(11),TTCUM(11)
	INTEGER NNO,NS,SI,NNT,BCCS
	INTEGER YNOD,LS,KK,STEPS0, YN
        INTEGER DST,ZNOD,J,XY,O
	INTEGER P,BCC,XIND,YIND,gridratio
	REAL*8 KINS,KINS2,ENMS,MGHS,DMPENS,PSUMS
	REAL*8 WENS,GSUMS,DPES,BCES
	REAL*8 XE,DEX,DEY,RDE,MML,XI,ZI,YI
	REAL*8 Z,AVE,ROUGH,FG,SCA,MGH,MGH0,DPE
	REAL*8 KX(12,12),KY(12,12),KZ(12,12),K(12,12)
	REAL*8 DDX,DDY,DDZ,DX,DY,DZ,DL,DTX,DTY,DTZ
	REAL*8 X1,Y1,Z1,X2,Y2,Z2,DXL,DYL,DZL,DDL,RLS
	REAL*8 MAXX,MAXY,MAXZ,MAXUT
	REAL*8 V1,V2,V3,MAXV,EF0,GL,WL,SLIN,PI,SUB
	REAL*8 SSB,CSB,SQB,LNN,SCL,DAMP1,DAMP2,DRAG
        REAL*8 fractime
	REAL RAN(NOCON)
        INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),maxid,neighcount
        INTEGER rc,ntasks_init,ierr,SEED,SEEDI,OUTINT,RESOUTINT,&
             NTOT,FXC
        INTEGER, ALLOCATABLE :: NCN(:),CN(:,:),CNPart(:,:), particles_G(:),&
             neighparts(:), NANS(:,:),NANPart(:),FXF(:,:)

        INTEGER, DIMENSION(8) :: datetime
        LOGICAL :: BedZOnly,FileExists,PrintTimes,StrictDomain,DoublePrec,CSVOutput
        LOGICAL :: GeomMasked,FixLat,FixBack,doShearLine,aboveShearLine
        LOGICAL, ALLOCATABLE :: LostParticle(:)
        CHARACTER(LEN=256) INFILE, geomfile, runname, wrkdir, resdir,restname,outstr

        TYPE(NAN_t) :: NANS_prev
        TYPE(NTOT_t) :: NTOT_prev,ND
!        TYPE(EF_t) :: EFS
        TYPE(NEI_t) :: NeighbourID
        TYPE(UT_t) :: UT, UTM
        TYPE(UT2_t) :: UT2, UTM2
        TYPE(NRXF_t) :: NRXF
        TYPE(NRXF2_t) :: NRXF2
        TYPE(FXF_t) :: NDL
        TYPE(InvPartInfo_t), ALLOCATABLE :: InvPartInfo(:)

        CALL MPI_INIT(rc)
        IF (rc /= MPI_SUCCESS) THEN
        WRITE(*,*) 'MPI initialisation failed'
        STOP
        END IF
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, rc)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, rc)
        ntasks_init = ntasks

        PrintTimes = .FALSE.
        DebugMode = .TRUE.

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

        CALL ReadInput(INFILE, runname, wrkdir, resdir, geomfile, PRESS, MELT, UC, DT, S, GRAV, &
             RHO, RHOW, EF0, LS, SUB, GL, SLIN, doShearLine, MLOAD, FRIC, REST, restname, POR, &
             SEEDI, DAMP1, DAMP2, DRAG, BedIntConst, BedZOnly, OUTINT, RESOUTINT, MAXUT, SCL, &
             WL, STEPS0,GRID, fractime,StrictDomain,DoublePrec,CSVOutput,GeomMasked,FixLat,FixBack)

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
	NRXF%M(:,:)=-1000.0
	NRXF%L(:,:)=-1000.0
	NRXF%R(:,:)=-1000.0
	NRXF%F(:,:)=-1000.0
	NRXF%B(:,:)=-1000.0
	NRXF%FL(:,:)=-1000.0
	NRXF%FR(:,:)=-1000.0
	NRXF%BL(:,:)=-1000.0
	NRXF%BR(:,:)=-1000.0

!Initialize UT
        UT%M(:) = 0.0
        UT%L(:) = 0.0
        UT%R(:) = 0.0
        UT%F(:) = 0.0
        UT%B(:) = 0.0
        UT%FR(:) = 0.0
        UT%FL(:) = 0.0
        UT%BR(:) = 0.0
        UT%BL(:) = 0.0

!inclination of the domain - not really used

	SSB=SIN(SUB*PI/2.0)
	CSB=COS(SUB*PI/2.0)
	SQB=SQRT(1.0+SIN(SUB*PI/2.0)**2)
	RLS=LS

	IF (REST.EQ.0) THEN

	RY0=1

!TODO - automatically translate and rotate the input data 
!     (interp required <- either at load time or within BIPINT)
!     Write out translation and rotation matrices to REST?

        !Go to glas.f90 to make the grid
	CALL FIBG3(NN,NTOT,NANS,NRXF2,NANPart,particles_G, NCN, CN, CNPart, InvPartInfo, &
        neighcount, LS, wrkdir,geomfile,SCL,GRID,MELT,WL,UC,&
        StrictDomain)
 
        IF(DebugMode) PRINT *,myid,' made it out of FIBG3 alive!'

 !TODO - allocate these to same size as NRXF
        !Allocate point data structures & pointers
        CALL PointDataInit(UT2,NN,part_expand)
        CALL PointDataInit(UTM2,NN,part_expand)

 !TODO - test output of NCN, CN, particles_g, NRXF

        CALL MPI_BARRIER(MPI_COMM_ACTIVE,ierr)

        CALL GetBBoxes(NRXF2, UT2, NN, BBox, PBBox)

        !PRINT *, myid,' neighparts: ',neighparts(1:neighcount)

        MAXX = BBox(2)
        MAXY = BBox(4)
        MAXZ = BBox(6) !don't really use these...

        IF(DebugMode) PRINT *,myid,' Got bbox: ',PBBox(:,myid+1)

	write(*,17) myid,NTOT
        CALL MPI_BARRIER(MPI_COMM_ACTIVE,ierr)

        !NN - particles in each core
        CALL MPI_ALLGATHER(NN,1,MPI_INTEGER,&
        PNN,1,MPI_INTEGER,MPI_COMM_ACTIVE,ierr)

        !TODO - links to CSV output - replace
        CALL MPI_ALLGATHER(NTOT_prev % M,1,MPI_INTEGER,&
        NTOTW,1,MPI_INTEGER,MPI_COMM_ACTIVE,ierr)

        IF(DebugMode) PRINT *,myid,'Wave checkpoint 2'

        !MATHS!
        !Iterate over: this part nodes, edge nodes

        ! CT - see notes - it's the cumulative translation of all the beams
        CT = 0.0

	T=0 !T is time
	ENM=0.0
	WEN=0.0
	ENM0=0

        VDP = 0.0     !VDP = drag coefficient
        UT%M = 0.0    ! UT%M = current displacement
        UTM%M = 0.0   ! UTM%M = previous displacement


	ELSE         !If restarting, read from restart file instead

        INQUIRE( FILE=TRIM(wrkdir)//'/'//TRIM(restname)//'_REST0'//na(myid), EXIST=FileExists ) 
        IF(.NOT. FileExists) CALL FatalError("Running with too many cores!&
             &(some restart files don't exist)")

        OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(restname)//'_REST0'//na(myid),STATUS='OLD')
        READ(117+myid,*) NN,NTOT,BCC
        !Allocate the structure holding the point data
        CALL PointDataInit(NRXF2,NN,part_expand)
        READ(117+myid,*) NRXF2%cstrt,NRXF2%NC
        READ(117+myid,*) MAXX,MAXY,MAXZ,DMPEN,ENM0
        READ(117+myid,*) DPE,BCE,MGH0,GSUM0,PSUM,T,RY0
        CLOSE(117+myid)

        !TODO -sort this nightmare
        IF(NRXF2%cstrt + NRXF2%NC > SIZE(NRXF2%A,2)) CALL ResizePointData(NRXF2,1.5_8)

        CALL MPI_ALLGATHER(NN,1,MPI_INTEGER,&
        PNN,1,MPI_INTEGER,MPI_COMM_ACTIVE,ierr)

        !TODO - links to CSV output - replace
        CALL MPI_ALLGATHER(NTOT_prev % M,1,MPI_INTEGER,&
        NTOTW,1,MPI_INTEGER,MPI_COMM_ACTIVE,ierr)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(restname)//'_REST1'//na(myid),STATUS='OLD')
	DO I=1,NTOT
	READ(117+myid,*) CT(12*I-11),CT(12*I-10),CT(12*I-9),&
      	CT(12*I-8),CT(12*I-7),CT(12*I-6)
	READ(117+myid,*) CT(12*I-5),CT(12*I-4),CT(12*I-3),&
      	CT(12*I-2),CT(12*I-1),CT(12*I-0)
	END DO
	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(restname)//'_REST2'//na(myid),STATUS='OLD')
	DO I=1,NN
	READ(117+myid,*) UT%M(6*I-5),UT%M(6*I-4),UT%M(6*I-3),&
      	UT%M(6*I-2),UT%M(6*I-1),UT%M(6*I-0)
	READ(117+myid,*) UTM%M(6*I-5),UTM%M(6*I-4),UTM%M(6*I-3),&
      	UTM%M(6*I-2),UTM%M(6*I-1),UTM%M(6*I-0)
	END DO
	CLOSE (117+myid)

        !Read particle init pos from file if restarting
 !TODO - writing of NODFIL2 has been modified to include neighbouring connected points
 ! so this is inconsistent
	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/NODFIL2'//na(myid),STATUS='UNKNOWN')
	DO I=1,NN
          READ(117+myid,*) IX,X,Y,Z,M
          NRXF2%M(1,IX)=X
          NRXF2%M(2,IX)=Y
          NRXF2%M(3,IX)=Z
        END DO
	CLOSE (117+myid)

        END IF 
        

! ============== END IF RESTART ==============


        IF(DebugMode) PRINT *,myid,' checkpoint 3'

        !initialize some data

        !Keep track of any particles which leave the domain
        ALLOCATE(LostParticle(NN),&
             EFS(NTOT))
        LostParticle = .FALSE.

        ALLOCATE(EFC(SIZE(NRXF2%A,2))) !to hold our & other nodes
        EFC = SCL*EF0

	DO I=1,NODM
	EN(I)=0.0
	END DO


        !Set bed -1000.0 everywhere
        DO I=-100,2000
        DO J=-100,2000
        BED(I,J)=-1000.0
        ENDDO
        ENDDO

        !Geometry mask: 1=glacier, 2=fjord, 0=bedrock/outside domain
        GEOMMASK = 0

        !Read the geometry and friction
        !into the grids BED, SUF, and FBED
        OPEN(UNIT=400,file=TRIM(geomfile),STATUS='UNKNOWN')
        READ(400,*) NM2
	DO I=1,NM2
        IF(GeomMasked) THEN
          READ(400,*) X,Y,S1,B2,B1,Z1,mask
        ELSE
          READ(400,*) X,Y,S1,B2,B1,Z1
        END IF
!        X=X-2000.0
!        Y=Y-7000.0
        XK=INT(X/GRID)
        YK=INT(Y/GRID)
        IF (XK.GE.-100.AND.YK.GE.-100) BED(XK,YK)=B1
        IF (XK.GE.-100.AND.YK.GE.-100) SUF(XK,YK)=S1
        IF (XK.GE.-100.AND.YK.GE.-100) FBED(XK,YK)=FRIC*SCL*SCL*Z1
        IF (XK.GE.-100.AND.YK.GE.-100.AND.GeomMasked) GEOMMASK(XK,YK)=mask

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

        IF(DebugMode) PRINT *,myid,' checkpoint 5'

        !Initialize drag/friction
        ALLOCATE(MFIL(NN))
	DO I=1,NN

        MFIL(I)=MN
	X=NRXF2%M(1,I)+UT2%M(6*I-5)
        Y=NRXF2%M(2,I)+UT2%M(6*I-4)
        Z=NRXF2%M(3,I)+UT2%M(6*I-3)
        I1=X/GRID-INT(X/GRID)
	I2=Y/GRID-INT(Y/GRID)
        CALL BIPINT(I1,I2,BED(INT(X/GRID),INT(Y/GRID)),BED(INT(X/GRID)+1,INT(Y/GRID)),&
        BED(INT(X/GRID),INT(Y/GRID)+1),BED(INT(X/GRID)+1,INT(Y/GRID)+1),ZB)
        
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
	END DO

!This code reads in connection information for each partition
!and works out whether to randomise or receive the EFS value

 IF(DebugMode) PRINT *,myid,' about to go to ExchangeConnPoints'
 CALL ExchangeConnPoints(NANS, NRXF2, InvPartInfo, UT2)
 IF(DebugMode)  PRINT *,myid,' finished exchange.'
 

 IF(REST==1) THEN
   OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FS'//na(myid),STATUS='UNKNOWN')
   DO I=1,NTOT
     READ(117+myid,*) N1,N2,P1,X1,Y1,Z1,X2,Y2,Z2,E
     NANS(1,I)=N1
     NANS(2,I)=N2
     NANPart(I)=P1
     EFS(I)=E
   END DO
   CLOSE (117+myid)
 ELSE
   DO I=1,NTOT

     aboveShearLine = .TRUE.
     IF(doShearLine) THEN
       N1 = NANS(1,I)
       N2 = NANS(2,I)
       X1 = NRXF2%A(1,N1)
       Y1 = NRXF2%A(2,N1)
       Z1 = NRXF2%A(3,N1)

       !TODO - could save time here if SLIN not specified
       I1=X1/GRID-INT(X1/GRID)
       I2=Y1/GRID-INT(Y1/GRID)
       CALL BIPINT(I1,I2,SUF(INT(X1/GRID),INT(Y1/GRID)),SUF(INT(X1/GRID)+1,INT(Y1/GRID)),&
            SUF(INT(X1/GRID),INT(Y1/GRID)+1),SUF(INT(X1/GRID)+1,INT(Y1/GRID)+1),ZS)
     
       aboveShearLine = ABS(Z1-ZS).LT.SLIN
     END IF
     
     IF (RAN(I).LT.1.0-POR.AND.aboveShearLine) THEN
       EFS(I)=EF0
     ELSE
       EFS(I)=0.1
     ENDIF
   END DO

   !Share randomly generated EFS with other parts to avoid conflict
   CALL ExchangeEFS(NANS, NANPart, NRXF2, InvPartInfo, EFS)
 END IF


!============================================================
!================= START THE TIME LOOP ======================
!============================================================

        CALL CPU_TIME(T1)
	TS1=0.0

        TTCUM=0.0
        TT=0.0


	DO 100 RY=RY0,RY0+STEPS0 

        CALL CPU_TIME(TT1)
        TT(1) = TT1
        TTCUM(1) = TTCUM(1) + (TT(1) - TT(11))

        !TODO  TIME 1

        !TODO - Use this for passing relevant points
        CALL GetBBoxes(NRXF2, UT2, NN, BBox, PBBox)

        CALL FindCollisions(NRXF2,UT2,NN,BBox,SCL,LNN)

!============================ UT TRANSMIT ===============================

      ! dest=myid+1
      ! source=myid-1

      !   tag=142
      ! IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      ! CALL MPI_Send(UT%M,6*NN,MPI_DOUBLE_PRECISION,&
      ! dest,tag,MPI_COMM_ACTIVE,ierr)
      ! IF (MOD(myid,ntasks/YN).ne.0)&
      ! CALL MPI_Recv(UT%L,6*PNN(source),MPI_DOUBLE_PRECISION,&
      ! source,tag,MPI_COMM_ACTIVE,stat,ierr)

      !....... etc for every direction - REPLACE THIS

!============================ END UT TRANSMIT ===============================

      IF(DebugMode) PRINT *,myid,' about to go to ExchangeProxPoints'
      CALL ExchangeProxPoints(NRXF2, UT2, NN, SCL)
      IF(DebugMode) PRINT *,myid,' finished  ExchangeProxPoints'

!============================== END REPLACEMENT =============================

      !TODO  TIME 2
      CALL CPU_TIME(TT(2))
      TTCUM(2) = TTCUM(2) + (TT(2) - TT(1))

      !Every 250 steps, reconstruct neighbourhood list
	IF (MOD(RY,250).EQ.1.OR.RY.EQ.RY0) THEN
     	CALL DIST(NN,UT,ND,NRXF,NDL,SCL,PNN,YN)
	END IF

!      write(*,17) RY,myid,ND,NCL,NDR,NDF,NDB,NDBL,NDBL,NDFR

       !TODO  TIME 3
       CALL CPU_TIME(TT(3))
       TTCUM(3) = TTCUM(3) + (TT(3) - TT(2))

       !circ checks which particles are really in contact and computes the forces
	CALL CIRC(ND,NN,NRXF,UT,FRX,FRY,FRZ,&
      	T,RY,DT,WE,EFC,FXF,FXC,NDL,LNN,SCL)

!      write(*,17) RY,myid,FXC%M,FXC%L,FXC%B,FXC%F,FXC%FR,FXC%FL,FXC%BR,FXC%BL

       !TODO  TIME 4
       CALL CPU_TIME(TT(4))
       TTCUM(4) = TTCUM(4) + (TT(4) - TT(3))

       !Calculates elastic forces from beams. Stiffness matrix K
	CALL EFFLOAD(S,NTOT_prev,NN,T,DT,MN,JS,DMP,DMP2,UT,UTM,R,EN,RY,&
      	FXF,FXC,VDP,DPE,EFS,NANS_prev,NRXF,MFIL,CT,&
      	LNN,PNN,YN)


        CALL MPI_BARRIER(MPI_COMM_ACTIVE,ierr)

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


!------------------ Compute bed interaction, drag, gravity, buoyancy -------------------
!------------------   and new displacement UTP -----------------------------------------


	DO I=1,NN !Cycle over particles

        !TODO  TIME 6
        CALL CPU_TIME(TT(6))
        IF(I>1) THEN
          TTCUM(6) = TTCUM(6) + (TT(6) - TT(9))
        END IF

	X=NRXF%M(1,I)+UT%M(6*I-5)  !<- x,y,z actual positions, NRXF%M = original, UT%M = displacement
	Y=NRXF%M(2,I)+UT%M(6*I-4)
	Z=NRXF%M(3,I)+UT%M(6*I-3)

	WSY(I)=0.0

! If you want to add a pressure from the backwall, instead of just free.
! Something like hydrostatic pressure - scaled

!	IF (Y.LT.1000.0.AND.(X.GT.8000.0.AND.X.LT.34000.0)) THEN
!          IF (Z.GT.WL+Y*SSB) THEN
!          WSY(I)=PRESS*1.0e+01*SCL**3*(MAXZ-NRXF%M(3,I))
!          PSUM=PSUM+WSY(I)*(UT%M(6*I-4)-UTM%M(6*I-4))
!          ELSE
!          WSY(I)=PRESS*(1.0e+01*SCL**3*(MAXZ-NRXF%M(3,I))-1.1e+01*SCL**3*(WL+Y*SSB-NRXF%M(3,I)))
!          PSUM=PSUM+WSY(I)*(UT%M(6*I-4)-UTM%M(6*I-4))
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
	KIN=KIN+0.5*MN*((UT%M(6*I-5)-UTM%M(6*I-5))/DT)**2
	KIN=KIN+0.5*MN*((UT%M(6*I-4)-UTM%M(6*I-4))/DT)**2
	KIN=KIN+0.5*MN*((UT%M(6*I-3)-UTM%M(6*I-3))/DT)**2
	KIN2=KIN2+0.5*JS*((UT%M(6*I-2)-UTM%M(6*I-2))/DT)**2
	KIN2=KIN2+0.5*JS*((UT%M(6*I-1)-UTM%M(6*I-1))/DT)**2
	KIN2=KIN2+0.5*JS*((UT%M(6*I-0)-UTM%M(6*I-0))/DT)**2

!Calculate final UTP
	UTP(6*I-5)=(R(6*I-5)+FRX(I)+WSX(I))/((MFIL(I)/DT**2)+VDP(I)/(2.*DT))
	UTP(6*I-4)=(R(6*I-4)+FRY(I)+BOYY(I)+WSY(I))/((MFIL(I)/DT**2)+VDP(I)/(2.*DT))
	UTP(6*I-3)=(R(6*I-3)+BOYZ(I)+FRZ(I))/((MFIL(I)/DT**2)+VDP(I)/(2.*DT))
	UTP(6*I-2)=R(6*I-2)/((MFIL(I)*JS/MN)/DT**2+VDP(I)/(2.*DT))
	UTP(6*I-1)=R(6*I-1)/((MFIL(I)*JS/MN)/DT**2+VDP(I)/(2.*DT))
	UTP(6*I-0)=R(6*I-0)/((MFIL(I)*JS/MN)/DT**2+VDP(I)/(2.*DT))

        !Check that particles haven't left the domain
        !and freeze them if they have!
         XIND = INT((NRXF%M(1,I) + UTP(6*I-5))/GRID)
         YIND = INT((NRXF%M(2,I) + UTP(6*I-4))/GRID)
         IF(XIND > 2000 .OR. XIND < -100 .OR. YIND > 2000 .OR. YIND < -100 .OR. &
             ABS(UTP(6*I-5)) > MAXUT .OR. ABS(UTP(6*I-4)) > MAXUT .OR. &
             ABS(UTP(6*I-3)) > MAXUT) THEN
           UTP(6*I-5) = UT%M(6*I-5)
           UTP(6*I-4) = UT%M(6*I-4)
           UTP(6*I-3) = UT%M(6*I-3)
           UTP(6*I-2) = UT%M(6*I-2)
           UTP(6*I-1) = UT%M(6*I-1)
           UTP(6*I-0) = UT%M(6*I-0)
           IF(.NOT. LostParticle(I)) THEN
             PRINT *, myid, " Lost a particle : ",I," at time: ",T
             LostParticle(i) = .TRUE.
           END IF
         END IF

!	IF (X.LT.4000.0.AND.Y.LT.4000.0) THEN
!	IF ((ZB.GT.WL-2.0*SCL.OR.Z-ZB.LT.5.0*SCL).AND.ABS(Z-ZB).LT.2.0*SCL) THEN

!TODO - standardise and unhardcode this
!Freeze particles near bed if friction is too high
        ! IF (ABS(ZB-Z).LT.SCL*2.0.AND.(VDP(I).GE.1.0e+08*SCL*SCL.OR.T.LT.20.0)) THEN
        ! UTP(6*I-5)=UT%M(6*I-5)
        ! UTP(6*I-4)=UT%M(6*I-4)
        ! UTP(6*I-3)=UT%M(6*I-3)
        ! UTP(6*I-2)=UT%M(6*I-2)
        ! UTP(6*I-1)=UT%M(6*I-1)
        ! UTP(6*I-0)=UT%M(6*I-0)
	! ENDIF

        !This code allows glacier edge/inflow boundary to be frozen in XY
        !TODO (esp for melange) - convert to contact BC/elastic rebound?
        IF(FixBack .OR. FixLat) THEN
          XIND = INT((NRXF%M(1,I) + UTP(6*I-5))/GRID)
          YIND = INT((NRXF%M(2,I) + UTP(6*I-4))/GRID)
          gridratio = INT(MAX(SCL/GRID,1.0))

          !Freeze if near back plane
          IF(FixBack) THEN
            IF(ANY(GEOMMASK(XIND,YIND-(2*gridratio):YIND+(2*gridratio))==0)) THEN
              UTP(6*I-5)=UT%M(6*I-5)
              UTP(6*I-4)=UT%M(6*I-4)
            END IF
          END IF

          !Freeze if near edge of glacier
          IF(FixLat) THEN
            IF(ANY(GEOMMASK(XIND-(2*gridratio):XIND+(2*gridratio),YIND)==0)) THEN
              UTP(6*I-5)=UT%M(6*I-5)
              UTP(6*I-4)=UT%M(6*I-4)
            END IF
          END IF
        END IF

       !Compute damping
	DMPEN=DMPEN+VDP(I)*(UTP(6*I-5)-UTM%M(6*I-5))**2/(4*DT)
	DMPEN=DMPEN+VDP(I)*(UTP(6*I-4)-UTM%M(6*I-4))**2/(4*DT)
	DMPEN=DMPEN+VDP(I)*(UTP(6*I-3)-UTM%M(6*I-3))**2/(4*DT)
	DMPEN=DMPEN+VDP(I)*(UTP(6*I-2)-UTM%M(6*I-2))**2/(4*DT)
	DMPEN=DMPEN+VDP(I)*(UTP(6*I-1)-UTM%M(6*I-1))**2/(4*DT)
	DMPEN=DMPEN+VDP(I)*(UTP(6*I-0)-UTM%M(6*I-0))**2/(4*DT)
 

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

       END DO !loop over particles


!-------------- Gather and write energy -----------------

        IF (RY.EQ.1) THEN
        CALL MPI_ALLREDUCE(GSUM,GSUM0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_ACTIVE,ierr)
        CALL MPI_ALLREDUCE(MGH,MGH0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_ACTIVE,ierr)
        END IF

!	IF (RY.EQ.REST*STEPS0+1) THEN

	IF (RY.EQ.RY0) THEN
        CALL MPI_ALLREDUCE(ENM0,ENMS0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_ACTIVE,ierr)
	ENDIF

 	IF (MOD(RY,100).EQ.0) THEN
        CALL MPI_ALLREDUCE(KIN,KINS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_ACTIVE,ierr)
        CALL MPI_ALLREDUCE(KIN2,KINS2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_ACTIVE,ierr)
        CALL MPI_ALLREDUCE(ENM,ENMS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_ACTIVE,ierr)
        CALL MPI_ALLREDUCE(MGH,MGHS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_ACTIVE,ierr)
        CALL MPI_ALLREDUCE(DMPEN,DMPENS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_ACTIVE,ierr)
        CALL MPI_ALLREDUCE(WEN,WENS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_ACTIVE,ierr)
        CALL MPI_ALLREDUCE(GSUM,GSUMS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_ACTIVE,ierr)
        CALL MPI_ALLREDUCE(DPE,DPES,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_ACTIVE,ierr)
        CALL MPI_ALLREDUCE(BCE,BCES,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_ACTIVE,ierr)
        CALL MPI_ALLREDUCE(BCC,BCCS,1,MPI_INTEGER,MPI_SUM,MPI_COMM_ACTIVE,ierr)
        CALL MPI_ALLREDUCE(PSUM,PSUMS,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_ACTIVE,ierr)
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

!--------------- Check for fracture -------------------

!internally
!TODO - what needs to be passed here? see commented below
        DO I=1,NTOT_prev%M
	N1=NANS_prev % M(1,I)
	N2=NANS_prev % M(2,I)
        IF (EFS(I).GT.0.0) THEN
	DDX=NRXF%M(1,N1)+UT%M(6*N1-5)-NRXF%M(1,N2)-UT%M(6*N2-5)
	DDY=NRXF%M(2,N1)+UT%M(6*N1-4)-NRXF%M(2,N2)-UT%M(6*N2-4)
	DDZ=NRXF%M(3,N1)+UT%M(6*N1-3)-NRXF%M(3,N2)-UT%M(6*N2-3)
	DX=NRXF%M(1,N1)-NRXF%M(1,N2)
	DY=NRXF%M(2,N1)-NRXF%M(2,N2)
	DZ=NRXF%M(3,N1)-NRXF%M(3,N2)
	L=SQRT(DX**2+DY**2+DZ**2)
	DL=SQRT(DDX**2+DDY**2+DDZ**2)
	V1=ABS(UT%M(6*N1-2)-UT%M(6*N2-2))
	V2=ABS(UT%M(6*N1-1)-UT%M(6*N2-1))
	V3=ABS(UT%M(6*N1-0)-UT%M(6*N2-0))
	MAXV=MAX(V1,V2,V3)
	LOAD=((DL-L)+0.05*MAXV)
	IF (LOAD.GT.MML) MML=LOAD
	 IF (LOAD.GT.MLOAD.AND.T.GT.fractime) THEN
	 BCE=BCE+0.5*EFS(I)*S**2/LNN*(DL-L)**2
         EFS(I)=0.0
	 BCC=BCC+1
	 ENDIF
        ENDIF
 	END DO

!        IF (myid.ne.0.and.myid.ne.(ntasks-1)/2+1) THEN
!        DO I=1,NTOT_prev % L
!	N1=NANS_prev % L(1,I)
!	N2=NANS_prev % L(2,I)
!        IF (EFS % L(I).GT.0.0) THEN
!	DDX=NRXF%L(1,N1)+UT%L(6*N1-5)-NRXF%M(1,N2)-UT%M(6*N2-5)
!	DDY=NRXF%L(2,N1)+UT%L(6*N1-4)-NRXF%M(2,N2)-UT%M(6*N2-4)
!	DDZ=NRXF%L(3,N1)+UT%L(6*N1-3)-NRXF%M(3,N2)-UT%M(6*N2-3)
!	DX=NRXF%L(1,N1)-NRXF%M(1,N2)
!	DY=NRXF%L(2,N1)-NRXF%M(2,N2)
!	DZ=NRXF%L(3,N1)-NRXF%M(3,N2)
!	L=SQRT(DX**2+DY**2+DZ**2)
!	DL=SQRT(DDX**2+DDY**2+DDZ**2)
!	V1=ABS(UT%L(6*N1-2)-UT%M(6*N2-2))
!	V2=ABS(UT%L(6*N1-1)-UT%M(6*N2-1))
!	V3=ABS(UT%L(6*N1-0)-UT%M(6*N2-0))
!	MAXV=MAX(V1,V2,V3)
!	LOAD=((DL-L)+0.05*MAXV)
!	IF (LOAD.GT.MML) MML=LOAD
!	 IF (LOAD.GT.MLOAD.AND.T.GT.1.0) THEN
!	 BCE=BCE+0.5*EFS % L(I)*S**2/LNN*(DL-L)**2
!         EFS % L(I)=0.0
!	 BCC=BCC+1
!	 ENDIF
!	ENDIF
! 	END DO
!	ENDIF

!         IF (myid.ne.ntasks-1.and.myid.ne.(ntasks-1)/2) THEN
!         DO I=1,NTOT_prev % R
! 	N1=NANS_prev % R(1,I)
! 	N2=NANS_prev % R(2,I)
!         IF (EFS % R(I).GT.0.0) THEN
! 	DDX=NRXF%R(1,N1)+UT%R(6*N1-5)-NRXF%M(1,N2)-UT%M(6*N2-5)
! 	DDY=NRXF%R(2,N1)+UT%R(6*N1-4)-NRXF%M(2,N2)-UT%M(6*N2-4)
! 	DDZ=NRXF%R(3,N1)+UT%R(6*N1-3)-NRXF%M(3,N2)-UT%M(6*N2-3)
! 	DX=NRXF%R(1,N1)-NRXF%M(1,N2)
! 	DY=NRXF%R(2,N1)-NRXF%M(2,N2)
! 	DZ=NRXF%R(3,N1)-NRXF%M(3,N2)
! 	L=SQRT(DX**2+DY**2+DZ**2)
! 	DL=SQRT(DDX**2+DDY**2+DDZ**2)
! 	V1=ABS(UT%R(6*N1-2)-UT%M(6*N2-2))
! 	V2=ABS(UT%R(6*N1-1)-UT%M(6*N2-1))
! 	V3=ABS(UT%R(6*N1-0)-UT%M(6*N2-0))
! 	MAXV=MAX(V1,V2,V3)
! 	LOAD=((DL-L)+0.05*MAXV)
! 	IF (LOAD.GT.MML) MML=LOAD
! 	 IF (LOAD.GT.MLOAD.AND.T.GT.fractime) THEN
! 	 BCE=BCE+0.5*EFS % R(I)*S**2/LNN*(DL-L)**2
!          EFS % R(I)=0.0
! 	 BCC=BCC+1
! 	 ENDIF
! 	ENDIF
!  	END DO
! 	ENDIF

! !across the boundaries...
!         dest=myid+1
!         source=myid-1
!         tag=131
!         IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
!         CALL MPI_Send(EFS % R,NTOT_prev % R,MPI_DOUBLE_PRECISION,&
!         dest,tag,MPI_COMM_ACTIVE,ierr)
!         IF (MOD(myid,ntasks/YN).ne.0)&
!         CALL MPI_Recv(EFS % L,NTOT_prev % L,MPI_DOUBLE_PRECISION,&
!         source,tag,MPI_COMM_ACTIVE,stat,ierr)

!         IF (myid.le.(YN-1)*ntasks/YN) THEN
!         DO I=1,NTOT_prev % F
! 	N1=NANS_prev % F(1,I)
! 	N2=NANS_prev % F(2,I)
!         IF (EFS % F(I).GT.0.0) THEN
! 	DDX=NRXF%F(1,N1)+UT%F(6*N1-5)-NRXF%M(1,N2)-UT%M(6*N2-5)
! 	DDY=NRXF%F(2,N1)+UT%F(6*N1-4)-NRXF%M(2,N2)-UT%M(6*N2-4)
! 	DDZ=NRXF%F(3,N1)+UT%F(6*N1-3)-NRXF%M(3,N2)-UT%M(6*N2-3)
! 	DX=NRXF%F(1,N1)-NRXF%M(1,N2)
! 	DY=NRXF%F(2,N1)-NRXF%M(2,N2)
! 	DZ=NRXF%F(3,N1)-NRXF%M(3,N2)
! 	L=SQRT(DX**2+DY**2+DZ**2)
! 	DL=SQRT(DDX**2+DDY**2+DDZ**2)
! 	V1=ABS(UT%F(6*N1-2)-UT%M(6*N2-2))
! 	V2=ABS(UT%F(6*N1-1)-UT%M(6*N2-1))
! 	V3=ABS(UT%F(6*N1-0)-UT%M(6*N2-0))
! 	MAXV=MAX(V1,V2,V3)
! 	LOAD=((DL-L)+0.05*MAXV)
! 	IF (LOAD.GT.MML) MML=LOAD
! 	 IF (LOAD.GT.MLOAD.AND.T.GT.fractime) THEN
! 	 BCE=BCE+0.5*EFS % F(I)*S**2/LNN*(DL-L)**2
!          EFS % F(I)=0.0
! 	 BCC=BCC+1
! 	 ENDIF
! 	ENDIF
!  	END DO
! 	ENDIF

!         dest=myid+ntasks/YN
!         source=myid-ntasks/YN
!         tag=132
!         IF (myid.lt.(YN-1)*ntasks/YN)&
!         CALL MPI_Send(EFS % F,NTOT_prev % F,MPI_DOUBLE_PRECISION,&
!         dest,tag,MPI_COMM_ACTIVE,ierr)
!         IF (myid.ge.ntasks/YN)&
!         CALL MPI_Recv(EFS % B,NTOT_prev % B,MPI_DOUBLE_PRECISION,&
!         source,tag,MPI_COMM_ACTIVE,stat,ierr)

!         IF (myid.le.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
!         DO I=1,NTOT_prev % FL
! 	N1=NANS_prev % FL(1,I)
! 	N2=NANS_prev % FL(2,I)
!         IF (EFS % FL(I).GT.0.0) THEN
! 	DDX=NRXF%FL(1,N1)+UT%FL(6*N1-5)-NRXF%M(1,N2)-UT%M(6*N2-5)
! 	DDY=NRXF%FL(2,N1)+UT%FL(6*N1-4)-NRXF%M(2,N2)-UT%M(6*N2-4)
! 	DDZ=NRXF%FL(3,N1)+UT%FL(6*N1-3)-NRXF%M(3,N2)-UT%M(6*N2-3)
! 	DX=NRXF%FL(1,N1)-NRXF%M(1,N2)
! 	DY=NRXF%FL(2,N1)-NRXF%M(2,N2)
! 	DZ=NRXF%FL(3,N1)-NRXF%M(3,N2)
! 	L=SQRT(DX**2+DY**2+DZ**2)
! 	DL=SQRT(DDX**2+DDY**2+DDZ**2)
! 	V1=ABS(UT%FL(6*N1-2)-UT%M(6*N2-2))
! 	V2=ABS(UT%FL(6*N1-1)-UT%M(6*N2-1))
! 	V3=ABS(UT%FL(6*N1-0)-UT%M(6*N2-0))
! 	MAXV=MAX(V1,V2,V3)
! 	LOAD=((DL-L)+0.05*MAXV)
! 	IF (LOAD.GT.MML) MML=LOAD
! 	 IF (LOAD.GT.MLOAD.AND.T.GT.fractime) THEN
! 	 BCE=BCE+0.5*EFS % FL(I)*S**2/LNN*(DL-L)**2
!          EFS % FL(I)=0.0
! 	 BCC=BCC+1
! 	 ENDIF
! 	ENDIF
!  	END DO
! 	ENDIF

!         dest=myid+ntasks/YN-1
!         source=myid-ntasks/YN+1
!         tag=133
!         IF (myid.lt.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
!         CALL MPI_Send(EFS % FL,NTOT_prev % FL,MPI_DOUBLE_PRECISION,&
!         dest,tag,MPI_COMM_ACTIVE,ierr)
!         IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
!         CALL MPI_Recv(EFS % BR,NTOT_prev % BR,MPI_DOUBLE_PRECISION,& 
!         source,tag,MPI_COMM_ACTIVE,stat,ierr)

!         IF (myid.le.(YN-1)*ntasks/YN.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
!         DO I=1,NTOT_prev % FR
! 	N1=NANS_prev % FR(1,I)
! 	N2=NANS_prev % FR(2,I)
!         IF (EFS % FR(I).GT.0.0) THEN
! 	DDX=NRXF%FR(1,N1)+UT%FR(6*N1-5)-NRXF%M(1,N2)-UT%M(6*N2-5)
! 	DDY=NRXF%FR(2,N1)+UT%FR(6*N1-4)-NRXF%M(2,N2)-UT%M(6*N2-4)
! 	DDZ=NRXF%FR(3,N1)+UT%FR(6*N1-3)-NRXF%M(3,N2)-UT%M(6*N2-3)
! 	DX=NRXF%FR(1,N1)-NRXF%M(1,N2)
! 	DY=NRXF%FR(2,N1)-NRXF%M(2,N2)
! 	DZ=NRXF%FR(3,N1)-NRXF%M(3,N2)
! 	L=SQRT(DX**2+DY**2+DZ**2)
! 	DL=SQRT(DDX**2+DDY**2+DDZ**2)
! 	V1=ABS(UT%FR(6*N1-2)-UT%M(6*N2-2))
! 	V2=ABS(UT%FR(6*N1-1)-UT%M(6*N2-1))
! 	V3=ABS(UT%FR(6*N1-0)-UT%M(6*N2-0))
! 	MAXV=MAX(V1,V2,V3)
! 	LOAD=((DL-L)+0.05*MAXV)
! 	IF (LOAD.GT.MML) MML=LOAD
! 	 IF (LOAD.GT.MLOAD.AND.T.GT.fractime) THEN
! 	 BCE=BCE+0.5*EFS % FR(I)*S**2/LNN*(DL-L)**2
!          EFS % FR(I)=0.0
! 	 BCC=BCC+1
! 	 ENDIF
! 	ENDIF
!  	END DO
! 	ENDIF

!         dest=myid+ntasks/YN+1
!         source=myid-ntasks/YN-1
!         tag=134
!         IF (myid.lt.(YN-1)*ntasks/YN&
!         .AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
!         CALL MPI_Send(EFS % FR,NTOT_prev % FR,MPI_DOUBLE_PRECISION,&
!         dest,tag,MPI_COMM_ACTIVE,ierr)
!         IF (myid.ge.ntasks/YN.AND.MOD(myid,ntasks/YN).ne.0)&
!         CALL MPI_Recv(EFS % BL,NTOT_prev % BL,MPI_DOUBLE_PRECISION,&
!         source,tag,MPI_COMM_ACTIVE,stat,ierr)

! !        IF (myid.gt.(ntasks-1)/YN) THEN
! !        DO I=1,NTOT_prev % B
! !	N1=NANS_prev % B(1,I)
! !	N2=NANS_prev % B(2,I)
! !        IF (EFS % B(I).GT.0.0) THEN
! !	DDX=NRXF%B(1,N1)+UT%B(6*N1-5)-NRXF%M(1,N2)-UT%M(6*N2-5)
! !	DDY=NRXF%B(2,N1)+UT%B(6*N1-4)-NRXF%M(2,N2)-UT%M(6*N2-4)
! !	DDZ=NRXF%B(3,N1)+UT%B(6*N1-3)-NRXF%M(3,N2)-UT%M(6*N2-3)
! !	DX=NRXF%B(1,N1)-NRXF%M(1,N2)
! !	DY=NRXF%B(2,N1)-NRXF%M(2,N2)
! !	DZ=NRXF%B(3,N1)-NRXF%M(3,N2)
! !	L=SQRT(DX**2+DY**2+DZ**2)
! !	DL=SQRT(DDX**2+DDY**2+DDZ**2)
! !	V1=ABS(UT%F(6*N1-2)-UT%M(6*N2-2))
! !	V2=ABS(UT%F(6*N1-1)-UT%M(6*N2-1))
! !	V3=ABS(UT%F(6*N1-0)-UT%M(6*N2-0))
! !	MAXV=MAX(V1,V2,V3)
! !	LOAD=((DL-L)+0.05*MAXV)
! !	IF (LOAD.GT.MML) MML=LOAD
! !	 IF (LOAD.GT.MLOAD.AND.T.GT.1.0) THEN
! !	 BCE=BCE+0.5*EFS % B(I)*S**2/LNN*(DL-L)**2
! !         EFS % B(I)=0.0
! !	 BCC=BCC+1
! !	 ENDIF
! !	ENDIF
! ! 	END DO
! !	ENDIF

! !        IF (myid.gt.(ntasks-1)/YN.AND.MOD(myid,ntasks/YN).ne.0) THEN
! !        DO I=1,NTOT_prev % BL
! !	N1=NANS_prev % BL(1,I)
! !	N2=NANS_prev % BL(2,I)
! !        IF (EFS % BL(I).GT.0.0) THEN
! !	DDX=NRXF%BL(1,N1)+UT%BL(6*N1-5)-NRXF%M(1,N2)-UT%M(6*N2-5)
! !	DDY=NRXF%BL(2,N1)+UT%BL(6*N1-4)-NRXF%M(2,N2)-UT%M(6*N2-4)
! !	DDZ=NRXF%BL(3,N1)+UT%BL(6*N1-3)-NRXF%M(3,N2)-UT%M(6*N2-3)
! !	DX=NRXF%BL(1,N1)-NRXF%M(1,N2)
! !	DY=NRXF%BL(2,N1)-NRXF%M(2,N2)
! !	DZ=NRXF%BL(3,N1)-NRXF%M(3,N2)
! !	L=SQRT(DX**2+DY**2+DZ**2)
! !	DL=SQRT(DDX**2+DDY**2+DDZ**2)
! !	V1=ABS(UT%FL(6*N1-2)-UT%M(6*N2-2))
! !	V2=ABS(UT%FL(6*N1-1)-UT%M(6*N2-1))
! !	V3=ABS(UT%FL(6*N1-0)-UT%M(6*N2-0))
! !	MAXV=MAX(V1,V2,V3)
! !	LOAD=((DL-L)+0.05*MAXV)
! !	IF (LOAD.GT.MML) MML=LOAD
! !	 IF (LOAD.GT.MLOAD.AND.T.GT.1.0) THEN
! !	 BCE=BCE+0.5*EFS % BL(I)*S**2/LNN*(DL-L)**2
! !         EFS % BL(I)=0.0
! !	 BCC=BCC+1
! !	 ENDIF
! !	ENDIF
! ! 	END DO
! !	ENDIF

! !        IF (myid.gt.(ntasks-1)/YN
! !     1	.AND.MOD(myid,ntasks/YN).ne.ntasks/YN-1) THEN
! !        DO I=1,NTOT_prev % BR
! !	N1=NANS_prev % BR(1,I)
! !	N2=NANS_prev % BR(2,I)
! !        IF (EFS % BR(I).GT.0.0) THEN
! !	DDX=NRXF%BR(1,N1)+UT%BR(6*N1-5)-NRXF%M(1,N2)-UT%M(6*N2-5)
! !	DDY=NRXF%BR(2,N1)+UT%BR(6*N1-4)-NRXF%M(2,N2)-UT%M(6*N2-4)
! !	DDZ=NRXF%BR(3,N1)+UT%BR(6*N1-3)-NRXF%M(3,N2)-UT%M(6*N2-3)
! !	DX=NRXF%BR(1,N1)-NRXF%M(1,N2)
! !	DY=NRXF%BR(2,N1)-NRXF%M(2,N2)
! !	DZ=NRXF%BR(3,N1)-NRXF%M(3,N2)
! !	L=SQRT(DX**2+DY**2+DZ**2)
! !	DL=SQRT(DDX**2+DDY**2+DDZ**2)
! !	V1=ABS(UT%FR(6*N1-2)-UT%M(6*N2-2))
! !	V2=ABS(UT%FR(6*N1-1)-UT%M(6*N2-1))
! !	V3=ABS(UT%FR(6*N1-0)-UT%M(6*N2-0))
! !	MAXV=MAX(V1,V2,V3)
! !	LOAD=((DL-L)+0.05*MAXV)
! !	IF (LOAD.GT.MML) MML=LOAD
! !	 IF (LOAD.GT.MLOAD.AND.T.GT.1.0) THEN
! !	 BCE=BCE+0.5*EFS % BR(I)*S**2/LNN*(DL-L)**2
! !         EFS % BR(I)=0.0
! !	 BCC=BCC+1
! !	 ENDIF
! !	ENDIF
! ! 	END DO
! !	ENDIF

	DO I=1,6*NN
	UTM%M(I)=UT%M(I)
	UT%M(I)=UTP(I)
	END DO

!--------------- Start of output ----------------------

	IF (MOD(RY,OUTINT).EQ.1) THEN
          IF(PrintTimes) THEN
20          FORMAT(" Times: ",11F7.1)
            WRITE(*,20) TTCUM
          END IF

          NRY=INT(RY/OUTINT)

          IF(.NOT. CSVOutput) THEN

          CALL BinaryVTKOutput(NRY,resdir,runname,PNN,NRXF,UT,&
               UTM,NeighbourID,NANS_prev,NTOT_prev,DoublePrec)
          

          CALL BinarySTROutput(NRY,resdir,runname,NRXF,UT,&
               NeighbourID,NANS_prev,NTOT_prev,DoublePrec)
          
          ELSE
          !--------------- CSV Output ----------------------
          dest=0

          tag=151
          IF (myid.NE.0)&
          CALL MPI_Send(UT%M,6*NN,MPI_DOUBLE_PRECISION,&
          dest,tag,MPI_COMM_ACTIVE,ierr)

          tag=152
          IF (myid.NE.0)&
          CALL MPI_Send(NANS_prev % M,3*NTOT_prev%M,MPI_INTEGER,&
          dest,tag,MPI_COMM_ACTIVE,ierr)

          tag=161
          IF (myid.NE.0)&
          CALL MPI_Send(NRXF%M,3*NN,MPI_DOUBLE_PRECISION,&
          dest,tag,MPI_COMM_ACTIVE,ierr)

          IF (myid.EQ.0) THEN
          OPEN(UNIT=910,FILE=TRIM(resdir)//'/'//TRIM(runname)//'_JYR'//na(NRY)//'.csv',STATUS='UNKNOWN')
          OPEN(UNIT=920,FILE=TRIM(resdir)//'/'//TRIM(runname)//'_STR'//na(NRY)//'.csv',STATUS='UNKNOWN')

          DO KK=1,ntasks-1

          source=KK

          tag=151
          CALL MPI_Recv(UTW,6*PNN(source),MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_ACTIVE,stat,ierr)

          tag=152
          CALL MPI_Recv(NANW,3*NTOTW(source),MPI_INTEGER,source,tag,MPI_COMM_ACTIVE,stat,ierr)

          tag=161
          CALL MPI_Recv(NRXFW,3*PNN(source),MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_ACTIVE,stat,ierr)

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
	  X=NRXF%M(1,I)+UT%M(6*I-5)
	  Y=NRXF%M(2,I)+UT%M(6*I-4)
	  Z=NRXF%M(3,I)+UT%M(6*I-3)
          WRITE(910,12) X,Y,Z
          END DO

          DO I=1,NTOT_prev%M
	  X=(NRXF%M(1,NANS_prev % M(1,I))+UT%M(6*NANS_prev % M(1,I)-5)+&
            NRXF%M(1,NANS_prev % M(2,I))+UT%M(6*NANS_prev % M(2,I)-5))/2.0
	  Y=(NRXF%M(2,NANS_prev % M(1,I))+UT%M(6*NANS_prev % M(1,I)-4)+&
            NRXF%M(2,NANS_prev % M(2,I))+UT%M(6*NANS_prev % M(2,I)-4))/2.0
	  Z=(NRXF%M(3,NANS_prev % M(1,I))+UT%M(6*NANS_prev % M(1,I)-3)+&
            NRXF%M(3,NANS_prev % M(2,I))+UT%M(6*NANS_prev % M(2,I)-3))/2.0
	  
	N1=NANS_prev % M(1,I)
	N2=NANS_prev % M(2,I)
	DDX=NRXF%M(1,N1)+UT%M(6*N1-5)-NRXF%M(1,N2)-UT%M(6*N2-5)
	DDY=NRXF%M(2,N1)+UT%M(6*N1-4)-NRXF%M(2,N2)-UT%M(6*N2-4)
	DDZ=NRXF%M(3,N1)+UT%M(6*N1-3)-NRXF%M(3,N2)-UT%M(6*N2-3)
	DX=NRXF%M(1,N1)-NRXF%M(1,N2)
	DY=NRXF%M(2,N1)-NRXF%M(2,N2)
	DZ=NRXF%M(3,N1)-NRXF%M(3,N2)
	L=SQRT(DX**2+DY**2+DZ**2)
	DL=SQRT(DDX**2+DDY**2+DDZ**2)
	STR=(DL-L)/L

!          IF (Z.GT.MAXZ-120.0) WRITE(920,12) X,Y,STR
          WRITE(920,12) X,Y,Z,STR
          END DO

        ENDIF !myid==0

!          CALL PSNET(NTOT_prev%M,NN,myid,GL,WL,SUB,ntasks)
	  CLOSE(910)
	  CLOSE(920)

        END IF !CSV or Binary output
      
      ENDIF !output interval

!--------------- End of output --------------------

	IF (MOD(RY,RESOUTINT).EQ.1) THEN
            !Add restart stuff here
            CALL WriteRestart()
        END IF


!	IF (RY.EQ.RY0+STEPS0) THEN
!        CALL CRACK(NN,NTOT_prev % M,CCN,EFS % M,NANS_prev % M,NRXF%M,UT%M,WL,myid,ntasks,NANS_prev % R,NTOT_prev % R,EFS % R)
!
!        dest=0
!	tag=171
!        IF (myid.NE.0) CALL MPI_Send(CCN,NN,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_ACTIVE,ierr)
!
!        IF (myid.EQ.0) THEN
!	OPEN(UNIT=112, FILE='frag',STATUS='UNKNOWN')
!         DO KK=1,ntasks-1
!         source=KK
!         tag=171
!         CALL MPI_Recv(CCNW,NN,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_ACTIVE,stat,ierr)
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
21          FORMAT(" TStep: ",I10," Simulation time: ",F10.5," secs")
          !IF (myid.EQ.0) WRITE(*,21) RY,MML/MLOAD,BCC
          IF (myid.EQ.0) WRITE(*,21) RY,T
          BCC=0
	ENDIF

        !TODO  TIME 11
       CALL CPU_TIME(TT(11))
       TTCUM(11) = TTCUM(11) + (TT(11) - TT(10))

 100	CONTINUE !end of time loop

!=========================================================
!================= END OF TIME LOOP ======================
!=========================================================

        CALL CPU_TIME(T2)

        CALL WriteRestart()

        CALL MPI_FINALIZE(rc)

        DEALLOCATE(LostParticle)

  	STOP

CONTAINS

 SUBROUTINE WriteRestart()
   ! Write out the restart files
	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(runname)//'_REST0'//na(myid),STATUS='UNKNOWN')
	WRITE(117+myid,*) NN,NTOT,BCC
        WRITE(117+myid,*) NRXF2%cstrt, NRXF2%NC
	WRITE(117+myid,*) MAXX,MAXY,MAXZ,DMPEN,ENM+ENM0
	WRITE(117+myid,*) DPE,BCE,MGH0,GSUM0,PSUM,T,RY-1

	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(runname)//'_REST1'//na(myid),STATUS='UNKNOWN')
	DO I=1,NTOT
	WRITE(117+myid,*) CT(12*I-11),CT(12*I-10),CT(12*I-9),CT(12*I-8),CT(12*I-7),CT(12*I-6)
	WRITE(117+myid,*) CT(12*I-5),CT(12*I-4),CT(12*I-3),CT(12*I-2),CT(12*I-1),CT(12*I-0)
	END DO
	CLOSE (117+myid)

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(runname)//'_REST2'//na(myid),STATUS='UNKNOWN')
	DO I=1,NN
	WRITE(117+myid,*) UT%M(6*I-5),UT%M(6*I-4),UT%M(6*I-3),UT%M(6*I-2),UT%M(6*I-1),UT%M(6*I-0)
	WRITE(117+myid,*) UTM%M(6*I-5),UTM%M(6*I-4),UTM%M(6*I-3),UTM%M(6*I-2),UTM%M(6*I-1),UTM%M(6*I-0)
	END DO
	CLOSE (117+myid)

 !NOTE - if NN or NTOT changes, need to rewrite NODFIL2 here!

 !NTOT_prev - no connections
 !NAN - the particle numbers (1 and 2), 2*ntot array
 !NRXF%M - the original locations of all the particles
 !EF - Youngs modulus

!====================== FROM HERE ==========================

 !TODO - don't need to reread NRXF - could just communicate

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FS'//na(myid),STATUS='UNKNOWN')
	DO SI=1,NTOT
	WRITE(117+myid,*) NANS(1,SI),NANS(2,SI),NANPart(SI),&
           NRXF%M(1,NANS(1,SI)),NRXF%M(2,NANS(1,SI)),&
           NRXF%M(3,NANS(1,SI)),NRXF%M(1,NANS(2,SI)),&
           NRXF%M(2,NANS(2,SI)),NRXF%M(3,NANS(2,SI)),EFS(SI)
	END DO
	CLOSE (117+myid)

!======================== TO HERE ==========================


 END SUBROUTINE WriteRestart

END PROGRAM

