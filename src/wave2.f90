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
        USE Effl
        USE Utils

	IMPLICIT NONE
        INCLUDE 'mpif.h'
        INCLUDE 'na90.dat'
	REAL(KIND=dp), ALLOCATABLE :: EMM(:),BOYZ(:),BOYY(:),WSY(:),WSX(:)
	REAL(KIND=dp), ALLOCATABLE :: MFIL(:),EFC(:),EFS(:),VDP(:)
        REAL(KIND=dp), ALLOCATABLE :: FRX(:),FRY(:),FRZ(:),WE(:),CT(:)
        REAL(KIND=dp), ALLOCATABLE :: EN(:),UTP(:),UTW(:),R(:),NRXFW(:,:)
	REAL(KIND=dp) :: BED(-100:2000,-100:2000)
	REAL(KIND=dp) :: SUF(-100:2000,-100:2000)
	REAL(KIND=dp) :: FBED(-100:2000,-100:2000)
	REAL(KIND=dp) :: DIX,DIY,DIZ,FRIC,UC
	REAL(KIND=dp) :: ENM0,ENMS0,POR,GRID,BBox(6)
        REAL(KIND=dp) :: RHO,RHOW,GRAV, BedIntConst
        REAL(KIND=dp), ALLOCATABLE :: PBBox(:,:)
	INTEGER, ALLOCATABLE :: CCN(:),CCNW(:)
	INTEGER NANW(3,NOCON),mask,GEOMMASK(-100:2000,-100:2000)
	INTEGER REST,RY0,NM2,noprocs, counter
        INTEGER, POINTER :: countptr
	REAL(KIND=dp) M,MN,JS,DT,T,X,Y,E,GSUM,GSUM0
	REAL(KIND=dp) DMPEN,PSUM,KIN,KIN2,PRESS,MELT
	INTEGER I,N,NL,NN,NSIZE,STEPS,IX,IM,MS,N1,N2,P1,RY
        INTEGER cstrt,NC,pstrt,NP
	INTEGER PN,NRY,NTOTW(0:5000),XK,YK,ZK,Part,ID
        REAL(KIND=dp) :: L,ALF,MLOAD,DMP,VEL,G,S1,S2,M1,B1,B2
	REAL(KIND=dp) :: S,LOAD,DMP2,BCE,STR,I1,I2,ZB,ZS
	REAL(KIND=dp) :: T2,T1,TS1,TT1,TT2,ENM,WEN,ERSUM
	REAL(KIND=dp) :: TT(11),TTCUM(11)
	INTEGER NNO,NS,SI,NNT,BCCS
	INTEGER YNOD,LS,KK,STEPS0
        INTEGER DST,ZNOD,J,XY,O
	INTEGER P,BCC,XIND,YIND,gridratio
	REAL(KIND=dp) :: KINS,KINS2,ENMS,MGHS,DMPENS,PSUMS
	REAL(KIND=dp) :: WENS,GSUMS,DPES,BCES
	REAL(KIND=dp) :: XE,DEX,DEY,RDE,MML,XI,ZI,YI
	REAL(KIND=dp) :: Z,AVE,ROUGH,FG,SCA,MGH,MGH0,DPE
	REAL(KIND=dp) :: KX(12,12),KY(12,12),KZ(12,12),K(12,12)
	REAL(KIND=dp) :: DDX,DDY,DDZ,DX,DY,DZ,DL,DTX,DTY,DTZ
	REAL(KIND=dp) :: X1,Y1,Z1,X2,Y2,Z2,DXL,DYL,DZL,DDL,RLS
	REAL(KIND=dp) :: MAXX,MAXY,MAXZ,MAXUT
	REAL(KIND=dp) :: V1,V2,V3,MAXV,EF0,GL,WL,SLIN,PI,SUB
	REAL(KIND=dp) :: SSB,CSB,SQB,LNN,SCL,DAMP1,DAMP2,DRAG
        REAL(KIND=dp) :: fractime
	REAL, ALLOCATABLE :: RAN(:)
        INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),maxid,neighcount
        INTEGER rc,ntasks_init,ierr,SEED,SEEDI,OUTINT,RESOUTINT,&
             NTOT,FXC,ND
        INTEGER, ALLOCATABLE :: NCN(:),CN(:,:),CNPart(:,:), particles_G(:),&
             neighparts(:), NANS(:,:),NANPart(:),FXF(:,:), NDL(:,:),PNN(:)

        INTEGER, DIMENSION(8) :: datetime
        LOGICAL :: BedZOnly,FileExists,StrictDomain,DoublePrec,CSVOutput
        LOGICAL :: GeomMasked,FixLat,FixBack,doShearLine,aboveShearLine
        LOGICAL, ALLOCATABLE :: LostParticle(:), PartIsNeighbour(:)
        CHARACTER(LEN=256) INFILE, geomfile, runname, wrkdir, resdir,restname,outstr

!        TYPE(EF_t) :: EFS
        TYPE(NEI_t) :: NeighbourID
        TYPE(UT_t) :: UT, UTM
        TYPE(NRXF_t) :: NRXF
        TYPE(InvPartInfo_t), TARGET, ALLOCATABLE :: InvPartInfo(:)

        CALL MPI_INIT(rc)
        IF (rc /= MPI_SUCCESS) THEN
        WRITE(*,*) 'MPI initialisation failed'
        STOP
        END IF

        CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, rc)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, rc)
        ntasks_init = ntasks

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

!inclination of the domain - not really used

	SSB=SIN(SUB*PI/2.0)
	CSB=COS(SUB*PI/2.0)
	SQB=SQRT(1.0+SIN(SUB*PI/2.0)**2)
	RLS=LS

        ALLOCATE(PBBox(6,0:ntasks-1),&
             PartIsNeighbour(0:ntasks-1),&
             PNN(ntasks))

        PBBox = 0.0
        PartIsNeighbour = .FALSE.
        PNN = 0

!============= Read in the geometry from DEM ==============

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


!============= Generate the lattice if new run ==============

	IF (REST.EQ.0) THEN

	RY0=1

!TODO - automatically translate and rotate the input data 
!     (interp required <- either at load time or within BIPINT)
!     Write out translation and rotation matrices to REST?

        !Go to glas.f90 to make the grid
	CALL FIBG3(NN,NTOT,NANS,NRXF,NANPart,particles_G, NCN, CN, CNPart, InvPartInfo, &
        neighcount, LS, wrkdir,geomfile,SCL,GRID,MELT,WL,UC,&
        StrictDomain,GeomMasked)
 
        IF(DebugMode) PRINT *,myid,' made it out of FIBG3 alive!'

        !Allocate point data structures & pointers
        CALL PointDataInit(UT,NRXF)
        CALL PointDataInit(UTM,NRXF)

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

        CALL GetBBoxes(NRXF, UT, NN, NANS, NTOT, EFS, BBox, PBBox)

        MAXX = BBox(2)
        MAXY = BBox(4)
        MAXZ = BBox(6) !don't really use these...

        IF(DebugMode) PRINT *,myid,' Got bbox: ',PBBox(:,myid)

	PRINT *, 'Part: ',myid,' NN: ',NN,' NTOT: ',NTOT
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !NN - particles in each core
        CALL MPI_ALLGATHER(NN,1,MPI_INTEGER,&
        PNN,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        !TODO - links to CSV output - replace
        CALL MPI_ALLGATHER(NTOT,1,MPI_INTEGER,&
        NTOTW,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        IF(DebugMode) PRINT *,myid,'Wave checkpoint 2'

        !MATHS!
        !Iterate over: this part nodes, edge nodes

	T=0 !T is time
	ENM=0.0
	WEN=0.0
	ENM0=0

        UT%M = 0.0    ! UT%M = current displacement
        UTM%M = 0.0   ! UTM%M = previous displacement

        ALLOCATE(CT(NTOT*12),RAN(NTOT))
        CT = 0.0

        SEED=SEEDI+873*myid
        CALL RMARIN(SEED,0,0)
        CALL RANMAR(RAN,NTOT)

        !This code reads in connection information for each partition
        !and works out whether to randomise or receive the EFS value
        ALLOCATE(EFS(NTOT))
        DO I=1,NTOT
          
          aboveShearLine = .TRUE.
          IF(doShearLine) THEN
            N1 = NANS(1,I)
            N2 = NANS(2,I)
            X1 = NRXF%A(1,N1)
            Y1 = NRXF%A(2,N1)
            Z1 = NRXF%A(3,N1)

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
        CALL ExchangeEFS(NANS, NANPart, NRXF, InvPartInfo, EFS)

!======================== Read Restart Data ======================

	ELSE         !If restarting, read from restart file instead


   !TODO - not sure that NRXF%PartInfo etc is correctly filled in here...

        INQUIRE( FILE=TRIM(wrkdir)//'/'//TRIM(restname)//'_REST0'//na(myid), EXIST=FileExists ) 
        IF(.NOT. FileExists) CALL FatalError("Running with too many cores!&
             &(some restart files don't exist)")

        OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(restname)//'_REST0'//na(myid),STATUS='OLD')
        READ(117+myid,*) NN,cstrt,NC,pstrt,NP,NSIZE,NTOT,BCC
        READ(117+myid,*) MAXX,MAXY,MAXZ,DMPEN,ENM0
        READ(117+myid,*) DPE,BCE,MGH0,GSUM0,PSUM,T,RY0
        CLOSE(117+myid)

        !Allocate the structures holding the point data
        CALL PointDataInit(NRXF,NN,arrsize=NSIZE)
        NRXF % cstrt = cstrt
        NRXF % NC = NC
        NRXF % pstrt = pstrt
        !NRXF % NP = NP
        CALL PointDataInit(UT,NRXF)
        CALL PointDataInit(UTM,NRXF)

        CALL MPI_ALLGATHER(NN,1,MPI_INTEGER,&
        PNN,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        !TODO - links to CSV output - replace
        CALL MPI_ALLGATHER(NTOT,1,MPI_INTEGER,&
        NTOTW,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        ALLOCATE(CT(12*NTOT))          ! CT - cumulative translation of all the beams
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
	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/NODFIL2'//na(myid),STATUS='UNKNOWN')
	DO I=1,NN
          READ(117+myid,*) IX,X,Y,Z,M
          NRXF%M(1,IX)=X
          NRXF%M(2,IX)=Y
          NRXF%M(3,IX)=Z
          NRXF%PartInfo(1,IX) = myid
          NRXF%PartInfo(2,IX) = IX
        END DO
	CLOSE (117+myid)

        !Read in the partition & id of particles we share from other partitions
	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/ONODFIL2'//na(myid),STATUS='UNKNOWN')
	DO I=1,NC!+NP
          READ(117+myid,*) IX,Part,ID
          IF(IX >= NRXF%pstrt .OR. IX < NRXF%cstrt) THEN
            PRINT *,myid,' debug IX, pstrt, cstrt: ',IX,NRXF%pstrt, NRXF%cstrt
            CALL FatalError("Restart programming error - NRXF bounds")
          END IF
          NRXF%PartInfo(1,IX) = Part
          NRXF%PartInfo(2,IX) = ID
          PartIsNeighbour(Part) = .TRUE.
        END DO
	CLOSE (117+myid)

        !Construct InvPartInfo 
        CALL InvPartInfoInit(InvPartInfo, PartIsNeighbour)

        DO I=NRXF%cstrt, NRXF%cstrt + NRXF%NC - 1
          Part = NRXF%PartInfo(1,i)
          ID = NRXF%PartInfo(2,i)
          countptr => InvPartInfo(Part) % ccount
          countptr = countptr + 1

          IF(countptr > SIZE(InvPartInfo(Part) % ConnIDs)) THEN
            CALL ExpandIntArray(InvPartInfo(Part) % ConnIDs)
            CALL ExpandIntArray(InvPartInfo(Part) % ConnLocs)
          END IF

          InvPartInfo(Part) % ConnIDs(countptr) = ID
          InvPartInfo(Part) % ConnLocs(countptr) = I
        END DO

        ALLOCATE(EFS(NTOT), NANS(2,NTOT), NANPart(NTOT))
        OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FS'//na(myid),STATUS='UNKNOWN')
        DO I=1,NTOT
          READ(117+myid,*) N1,N2,P1,X1,Y1,Z1,X2,Y2,Z2,E
          NANS(1,I)=N1
          NANS(2,I)=N2
          NANPart(I)=P1
          EFS(I)=E
        END DO
        CLOSE (117+myid)

        IF(DebugMode) PRINT *,myid,' about to go to ExchangeConnPoints'
        CALL ExchangeConnPoints(NANS, NRXF, InvPartInfo, UT, UTM)
        IF(DebugMode)  PRINT *,myid,' finished exchange.'

        END IF 

        

! ============== END IF RESTART ==============


        IF(DebugMode) PRINT *,myid,' checkpoint 3'

        !initialize some data

        !Keep track of any particles which leave the domain
        ALLOCATE(LostParticle(NN),&
             VDP(NN),&
             FRX(NN),&
             FRY(NN),&
             FRZ(NN),&
             WE(NN),&
             EN(NN*6),&
             UTP(NN*6),&
             UTW(NN*6),&
             R(NN*6),&
             EMM(NN),&
             BOYZ(NN),&
             BOYY(NN),&
             WSY(NN),&
             WSX(NN),&
             NRXFW(3,MAXVAL(PNN)))

        !ALLOCATE(CCN(NN),CCNW(NN))

        LostParticle = .FALSE.
        VDP = 0.0     !VDP = drag coefficient
        EN = 0.0
        UTP = 0.0
        UTW = 0.0
        R = 0.0


        ALLOCATE(EFC(SIZE(NRXF%A,2))) !to hold our & other nodes
        EFC = SCL*EF0

        IF(DebugMode) PRINT *,myid,' checkpoint 5'

        !Initialize drag/friction
        ALLOCATE(MFIL(NN))
	DO I=1,NN

        MFIL(I)=MN
	X=NRXF%M(1,I)+UT%M(6*I-5)
        Y=NRXF%M(2,I)+UT%M(6*I-4)
        Z=NRXF%M(3,I)+UT%M(6*I-3)
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

 


!============================================================
!================= START THE TIME LOOP ======================
!============================================================

 IF(PrintTimes) THEN
   CALL CPU_TIME(T1)
   TS1=0.0
   TTCUM=0.0
   TT=0.0
 END IF

 DO 100 RY=RY0,RY0+STEPS0 

   IF(myid==0 .AND. PrintTimes) THEN
     PRINT *,'Time step: ',RY
     CALL CPU_TIME(TT1)
     TT(1) = TT1
     TTCUM(1) = TTCUM(1) + (TT(1) - TT(11))
     !TODO  TIME 1
   END IF

        CALL ExchangeConnPoints(NANS, NRXF, InvPartInfo, UT, UTM, .FALSE.) !Don't pass NRXF...

        CALL GetBBoxes(NRXF, UT, NN, NANS, NTOT, EFS, BBox, PBBox)

        CALL FindNeighbours(PBBox, PartIsNeighbour)

        CALL ExchangeProxPoints(NRXF, UT, UTM, NN, SCL, PBBox, InvPartInfo, PartIsNeighbour)


	IF (MOD(RY,250).EQ.1 .OR. (RY.EQ.RY0)) THEN
          IF(DebugMode) PRINT *,'About to find nearby particles'

          !Clear out invalid particles - those which were previously near our partition
          !but are no longer
          CALL ClearOldParticles(NRXF,UT,UTM,InvPartInfo)
          !Identify possible collisions
          CALL FindNearbyParticles(NRXF,UT,NN,BBox,SCL,LNN,ND,NDL)
        END IF

       !circ checks which particles are really in contact and computes the forces
	CALL FindCollisions(ND,NN,NRXF,UT,FRX,FRY,FRZ,&
          T,RY,DT,WE,EFC,FXF,FXC,NDL,LNN,SCL)

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       
        !TODO - should we still adopt the approach of finding nearby particles every N steps?
        ! - how efficient is octree vs prev strategy?

!============================ Old Prox/Collision Strategy ===============================
!Went like this:
! 1) Exchange *all* UT from all neighbours
! 2) Every 250 steps, check for near (but not neccesarily contacting points)
! 3) Every step, check for collision/interaction, compute FRX,Y,Z, store IDs FXF

!New strategy:
! 1) Use BBoxes to determine possible prox swaps (list of prox neighbour parts)
! 2) Swap nearby pionts using ExchangeProxPoints
! 3) Store those NRXF/UT values somehow?   <- issue - we don't want to swap NRXF every 
!                                             time, but sometimes NEW prox points will make this necessary
! 4) Use octree search to search for collision/interaction (no need for every 250 steps malarky)
!
! Make use of InvPartInfo to confirm which points we already have/which we need
! Note - do we need to use ExchangeConnPoints separately here, or do they form part of the same strategy? 
  
      ! dest=myid+1
      ! source=myid-1

      !   tag=142
      ! IF (MOD(myid,ntasks/YN).ne.ntasks/YN-1)&
      ! CALL MPI_Send(UT%M,6*NN,MPI_DOUBLE_PRECISION,&
      ! dest,tag,MPI_COMM_WORLD,ierr)
      ! IF (MOD(myid,ntasks/YN).ne.0)&
      ! CALL MPI_Recv(UT%L,6*PNN(source),MPI_DOUBLE_PRECISION,&
      ! source,tag,MPI_COMM_WORLD,stat,ierr)

      !....... etc for every direction - REPLACE THIS


!       !TODO  TIME 2
!       CALL CPU_TIME(TT(2))
!       TTCUM(2) = TTCUM(2) + (TT(2) - TT(1))

!       !Every 250 steps, reconstruct neighbourhood list
! 	IF (MOD(RY,250).EQ.1.OR.RY.EQ.RY0) THEN

!         ALLOCATE(NDL(2,NN*12))
!         CALL CPU_TIME(T1)
!      	CALL DIST(NN,UT,ND,NRXF,NDL,SCL,PNN,YN)
!         CALL CPU_TIME(T2)
!         PRINT *,myid,'Dist took: ',T2-T1,' secs'

! 	END IF
!         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

! !      write(*,17) RY,myid,ND,NCL,NDR,NDF,NDB,NDBL,NDBL,NDFR

!        !TODO  TIME 3
!        CALL CPU_TIME(TT(3))
!        TTCUM(3) = TTCUM(3) + (TT(3) - TT(2))

!        !circ checks which particles are really in contact and computes the forces
! 	CALL CIRC(ND,NN,NRXF,UT,FRX,FRY,FRZ,&
!       	T,RY,DT,WE,EFC,FXF,FXC,NDL,LNN,SCL)

! !      write(*,17) RY,myid,FXC%M,FXC%L,FXC%B,FXC%F,FXC%FR,FXC%FL,FXC%BR,FXC%BL

!============================ END Old Prox/Collision Strategy ===============================

       !TODO  TIME 4
        IF(PrintTimes) THEN
          CALL CPU_TIME(TT(4))
          TTCUM(4) = TTCUM(4) + (TT(4) - TT(3))
        END IF

       !Calculates elastic forces from beams. Stiffness matrix K
	CALL EFFLOAD(S,NTOT,NN,T,DT,MN,JS,DMP,DMP2,UT,UTM,R,EN,RY,&
      	FXF,FXC,VDP,DPE,EFS,NANS,NRXF,MFIL,CT,LNN)


        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

        IF(PrintTimes) THEN

          CALL CPU_TIME(TT2)
          TS1=TS1+(TT2-TT1)

          !TODO  TIME 5
          CALL CPU_TIME(TT(5))
          TTCUM(5) = TTCUM(5) + (TT(5) - TT(4))
        END IF

	WEN=0.0
	KIN=0.0
	KIN2=0.0
	MGH=0.0
	GSUM=0.0
	ENM=0.0


!------------------ Compute bed interaction, drag, gravity, buoyancy -------------------
!------------------   and new displacement UTP -----------------------------------------


	DO I=1,NN !Cycle over particles

        IF(PrintTimes) THEN
          !TODO  TIME 6
          CALL CPU_TIME(TT(6))
          IF(I>1) THEN
            TTCUM(6) = TTCUM(6) + (TT(6) - TT(9))
          END IF
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

        IF(PrintTimes) THEN
          CALL CPU_TIME(TT(7))
          TTCUM(7) = TTCUM(7) + (TT(7) - TT(6))
        END IF

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

       IF(PrintTimes) THEN
         CALL CPU_TIME(TT(8))
         TTCUM(8) = TTCUM(8) + (TT(8) - TT(7))
       END IF

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

        IF(PrintTimes) THEN
          CALL CPU_TIME(TT(9))
          TTCUM(9) = TTCUM(9) + (TT(9) - TT(8))
        END IF
       END DO !loop over particles


!-------------- Gather and write energy -----------------

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

        IF(PrintTimes) THEN
          CALL CPU_TIME(TT(10))
          TTCUM(10) = TTCUM(10) + (TT(10) - TT(9))
        END IF
!--------------- Check for fracture -------------------

!internally
!TODO - what needs to be passed here? see commented below
        DO I=1,NTOT
	N1=NANS(1,I)
	N2=NANS(2,I)
        IF (EFS(I).GT.0.0) THEN
	DDX=NRXF%A(1,N1)+UT%A(6*N1-5)-NRXF%A(1,N2)-UT%A(6*N2-5)
	DDY=NRXF%A(2,N1)+UT%A(6*N1-4)-NRXF%A(2,N2)-UT%A(6*N2-4)
	DDZ=NRXF%A(3,N1)+UT%A(6*N1-3)-NRXF%A(3,N2)-UT%A(6*N2-3)
	DX=NRXF%A(1,N1)-NRXF%A(1,N2)
	DY=NRXF%A(2,N1)-NRXF%A(2,N2)
	DZ=NRXF%A(3,N1)-NRXF%A(3,N2)
	L=SQRT(DX**2+DY**2+DZ**2)
	DL=SQRT(DDX**2+DDY**2+DDZ**2)
	V1=ABS(UT%A(6*N1-2)-UT%A(6*N2-2))
	V2=ABS(UT%A(6*N1-1)-UT%A(6*N2-1))
	V3=ABS(UT%A(6*N1-0)-UT%A(6*N2-0))
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

!--------------- Start of output ----------------------

	IF (MOD(RY,OUTINT).EQ.1) THEN
          IF(PrintTimes) THEN
20          FORMAT(" Times: ",11F7.1)
            WRITE(*,20) TTCUM
          END IF

          NRY=INT(RY/OUTINT)

          IF(.NOT. CSVOutput) THEN

          CALL BinaryVTKOutput(NRY,resdir,runname,PNN,NRXF,UT,&
               UTM,NANS,NTOT,NANPart,DoublePrec)
          

          CALL BinarySTROutput(NRY,resdir,runname,NRXF,UT,&
               NANS,NTOT,NANPart,DoublePrec)
          
          ELSE

            CALL FatalError("CSV Output not currently supported - needs to&
                 &be recoded for metis partitining")
          !--------------- CSV Output ----------------------
          dest=0

          !TODO - issue here - prev just sent our own beams, UT, NRXF, but now we would need
          !to send others too. Or just scrap CSV Output?
          tag=151
          IF (myid.NE.0)&
          CALL MPI_Send(UT%M,6*NN,MPI_DOUBLE_PRECISION,&
          dest,tag,MPI_COMM_WORLD,ierr)

          tag=152
          IF (myid.NE.0)&
          CALL MPI_Send(NANS,2*NTOT,MPI_INTEGER,&
          dest,tag,MPI_COMM_WORLD,ierr)

          tag=161
          IF (myid.NE.0)&
          CALL MPI_Send(NRXF%M,3*NN,MPI_DOUBLE_PRECISION,&
          dest,tag,MPI_COMM_WORLD,ierr)

          IF (myid.EQ.0) THEN
          OPEN(UNIT=910,FILE=TRIM(resdir)//'/'//TRIM(runname)//'_JYR'//na(NRY)//'.csv',STATUS='UNKNOWN')
          OPEN(UNIT=920,FILE=TRIM(resdir)//'/'//TRIM(runname)//'_STR'//na(NRY)//'.csv',STATUS='UNKNOWN')

          DO KK=1,ntasks-1

          source=KK

          tag=151
          CALL MPI_Recv(UTW,6*PNN(source),MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)

          tag=152
          CALL MPI_Recv(NANW,2*NTOTW(source),MPI_INTEGER,source,tag,MPI_COMM_WORLD,stat,ierr)

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
	  X=NRXF%M(1,I)+UT%M(6*I-5)
	  Y=NRXF%M(2,I)+UT%M(6*I-4)
	  Z=NRXF%M(3,I)+UT%M(6*I-3)
          WRITE(910,12) X,Y,Z
          END DO

          DO I=1,NTOT
	  X=(NRXF%M(1,NANS(1,I))+UT%M(6*NANS(1,I)-5)+&
            NRXF%M(1,NANS(2,I))+UT%M(6*NANS(2,I)-5))/2.0
	  Y=(NRXF%M(2,NANS(1,I))+UT%M(6*NANS(1,I)-4)+&
            NRXF%M(2,NANS(2,I))+UT%M(6*NANS(2,I)-4))/2.0
	  Z=(NRXF%M(3,NANS(1,I))+UT%M(6*NANS(1,I)-3)+&
            NRXF%M(3,NANS(2,I))+UT%M(6*NANS(2,I)-3))/2.0
	  
	N1=NANS(1,I)
	N2=NANS(2,I)
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

      !update UT
      DO I=1,6*NN
        UTM%M(I)=UT%M(I)
        UT%M(I)=UTP(I)
      END DO


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
21          FORMAT(" TStep: ",I10," Simulation time: ",F10.5," secs")
          !IF (myid.EQ.0) WRITE(*,21) RY,MML/MLOAD,BCC
          IF (myid.EQ.0) WRITE(*,21) RY,T
          BCC=0
	ENDIF

        IF(PrintTimes) THEN
          CALL CPU_TIME(TT(11))
          TTCUM(11) = TTCUM(11) + (TT(11) - TT(10))
        END IF

 100	CONTINUE !end of time loop

!=========================================================
!================= END OF TIME LOOP ======================
!=========================================================

        IF(PrintTimes) CALL CPU_TIME(T2)

        CALL WriteRestart()

        CALL MPI_FINALIZE(rc)

        DEALLOCATE(LostParticle)

  	STOP

CONTAINS

 SUBROUTINE WriteRestart()
        INTEGER :: counter

   ! Write out the restart files
	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/'//TRIM(runname)//'_REST0'//na(myid),STATUS='UNKNOWN')
	WRITE(117+myid,*) NRXF%NN,NRXF%cstrt,NRXF%NC,NRXF%pstrt,NRXF%NP,SIZE(NRXF%A,2),NTOT,BCC
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

        OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/ONODFIL2'//na(myid),STATUS='UNKNOWN')
        counter = 0
        DO I=NRXF%cstrt, NRXF%cstrt + NRXF%NC - 1
          IF(NRXF%PartInfo(1,i) == -1) CALL FatalError("Programming error in ONODFIL2")
          WRITE(117+myid,*) i,NRXF%PartInfo(:,i)
        END DO
	CLOSE (117+myid)

!====================== FROM HERE ==========================

 !TODO - don't need to reread NRXF - could just communicate

	OPEN(UNIT=117+myid,FILE=TRIM(wrkdir)//'/FS'//na(myid),STATUS='UNKNOWN')
	DO SI=1,NTOT
	WRITE(117+myid,*) NANS(1,SI),NANS(2,SI),NANPart(SI),&
           NRXF%A(1,NANS(1,SI)),NRXF%A(2,NANS(1,SI)),&
           NRXF%A(3,NANS(1,SI)),NRXF%A(1,NANS(2,SI)),&
           NRXF%A(2,NANS(2,SI)),NRXF%A(3,NANS(2,SI)),EFS(SI)
	END DO
	CLOSE (117+myid)

!======================== TO HERE ==========================


 END SUBROUTINE WriteRestart

END PROGRAM
