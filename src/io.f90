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

MODULE INOUT

  USE TypeDefs

  IMPLICIT NONE

CONTAINS

 SUBROUTINE ReadInput(INFILE, runname, wrkdir, resdir, geomfile, PRESS, MELT, UC, DT, S, GRAV, &
      RHO, RHOW, EF0, LS, SUB, GL, SLIN, doShearLine, MLOAD, FRIC, REST, restname, POR, SEEDI, DAMP1, &
      DAMP2, DRAG, ViscDist, ViscForce, BedIntConst, BedZOnly, OUTINT, RESOUTINT, &
      MAXUT, SCL, WL, STEPS0, GRID, fractime, &
      StrictDomain, DoublePrec, CSVOutput, GeomMasked, FixLat,FixBack, gotMelange, MelRunName)

   REAL(KIND=dp) :: PRESS, MELT, UC, DT, S, EF0, SUB, GL, SLIN, MLOAD, FRIC, POR
   REAL(KIND=dp) :: DAMP1, DAMP2, DRAG,MAXUT, SCL, WL, GRID, GRAV, RHO, RHOW, BedIntConst
   REAL(KIND=dp) :: fractime,viscforce,viscdist
   INTEGER :: REST, SEEDI, OUTINT, RESOUTINT, STEPS0, LS
   INTEGER :: readstat, i,incount
   CHARACTER(256) :: INFILE, geomfile, buff,VarName,VarValue,runname,MelRunName,wrkdir,&
        resdir,restname
   LOGICAL :: BedZOnly,StrictDomain,DoublePrec,CSVOutput,FileExists,FixLat,&
        FixBack,GeomMasked,doShearLine,gotMelange
   LOGICAL :: gotWL=.FALSE., gotSteps=.FALSE., gotSCL=.FALSE., &
        gotGrid=.FALSE.,gotName=.FALSE.,gotGeom=.FALSE.,gotRestName=.FALSE.

   OPEN(UNIT=112,FILE=infile,STATUS='old')
   incount = 0

   !Set default values
   PRESS = 0.0
   MELT = 0.0
   UC = 0.0
   DT = 1.0e-4
   S = 0.7
   EF0 = 1.0e+9
   LS = 100
   SUB = 0.0
   GL = -100.0
   SLIN = 2000.0
   doShearLine = .FALSE.
   MLOAD = 0.0002
   FRIC = 1.0
   REST = 0
   POR = 0.1
   SEEDI = 11695378
   DAMP1 = 1.0E4
   DAMP2 = 1.0E4
   DRAG = 1.0E1
   viscforce=1.0E4
   viscdist = 4.0E-2
   OUTINT = 20000
   RESOUTINT = 20000
   MAXUT = 1.0E6
   GRAV = 9.81
   RHO = 900.0
   RHOW = 1030.0
   BedIntConst = 1.0E8
   BedZOnly = .TRUE.
   wrkdir = './'
   resdir = './'
   fractime = 40.0
   StrictDomain = .TRUE.
   DoublePrec = .FALSE.
   CSVOutput = .FALSE.
   FixLat = .FALSE.
   FixBack = .TRUE.
   GeomMasked = .FALSE.
   DebugMode = .FALSE.
   PrintTimes = .FALSE.
   gotMelange = .FALSE.

   DO
     READ(112,"(A)", IOSTAT=readstat) buff
     IF(readstat > 0) STOP
     IF(readstat < 0) EXIT
     incount = incount+1

     !Ignore comments and blank lines
     IF(INDEX(TRIM(buff),'!') > 0) CYCLE
     IF (LEN_TRIM(buff) == 0) CYCLE
 
     i = INDEX(buff,'=')
     IF(i==0) PRINT *,'Format error in input file on line: ',incount

     VarName = buff(1:i-1)
     VarValue = buff(i+1:)

!     PRINT *, TRIM(ToLowerCase(VarName)),' has value: ',TRIM(VarValue)

     SELECT CASE (TRIM(ToLowerCase(VarName)))
     CASE ("density")
       READ(VarValue,*) RHO
     CASE ("water density")
       READ(VarValue,*) RHOW
     CASE("gravity")
       READ(VarValue,*) GRAV
     CASE("backwall pressure")
       READ(VarValue,*) PRESS
     CASE("submarine melt")
       READ(VarValue,*) MELT
     CASE("uc")
       READ(VarValue,*) UC
     CASE("timestep")
       READ(VarValue,*) DT
     CASE("width")
       READ(VarValue,*) S
     CASE("youngs modulus")
       READ(VarValue,*) EF0
     CASE("size")
       READ(VarValue,*) LS
     CASE("domain inclination")
       READ(VarValue,*) SUB
     CASE("water line")
       READ(VarValue,*) WL
       gotWL = .TRUE.
     CASE("grounding line")
       READ(VarValue,*) GL
     CASE("shear line")
       READ(VarValue,*) SLIN
       doShearLine = .TRUE.
     CASE("no timesteps")
       READ(VarValue,*) STEPS0
       gotSteps = .TRUE.
     CASE("max load")
       READ(VarValue,*) MLOAD
     CASE("friction scale")
       READ(VarValue,*) FRIC
     CASE("restart")
       READ(VarValue,*) REST
     CASE("scale")
       READ(VarValue,*) SCL
       gotSCL = .TRUE.
     CASE("grid")
       READ(VarValue,*) GRID
       gotGrid = .TRUE.
     CASE("porosity")
       READ(VarValue,*) POR
     CASE("random seed")
       READ(VarValue,*) SEEDI
     CASE("translational damping")
       READ(VarValue,*) DAMP1
     CASE("rotational damping")
       READ(VarValue,*) DAMP2
     CASE("viscous distance")
       READ(VarValue,*) viscdist
     CASE("viscous force")
       READ(VarValue,*) viscforce
     CASE("drag coefficient")
       READ(VarValue,*) DRAG
     CASE("output interval")
       READ(VarValue,*) OUTINT
     CASE("restart output interval")
       READ(VarValue,*) RESOUTINT
     CASE("maximum displacement")
       READ(VarValue,*) MAXUT
     CASE("run name")
       READ(VarValue,*) runname
       gotName = .TRUE.
     CASE("restart from run name")
       READ(VarValue,*) restname
       gotRestName = .TRUE.
     CASE("work directory")
       READ(VarValue,*) wrkdir
     CASE("geometry file")
       READ(VarValue,*) geomfile
       gotGeom = .TRUE.
     CASE("geometry file has mask")
       READ(VarValue,*) GeomMasked
     CASE("results directory")
       READ(VarValue,*) resdir
     CASE("bed stiffness constant")
       READ(VarValue,*) BedIntConst
     CASE("bed z only")
       READ(VarValue,*) BedZOnly
     CASE("fracture after time")
       READ(VarValue,*) fractime
     CASE("strict domain interpolation")
       READ(VarValue,*) StrictDomain
     CASE("double precision output")
       READ(VarValue,*) DoublePrec
     CASE("csv output")
       READ(VarValue,*) CSVOutput
     CASE("fixed lateral margins")
       READ(VarValue,*) FixLat
     CASE("fixed inflow margin")
       READ(VarValue,*) FixBack
     CASE("debug mode")
       READ(VarValue,*) DebugMode
     CASE("print times")
       READ(VarValue,*) PrintTimes
     CASE("melange run name")
       READ(VarValue,*) MelRunName
       gotMelange = .TRUE.
     CASE DEFAULT
       PRINT *,'Unrecognised input: ',TRIM(VarName)
       STOP
     END SELECT

   END DO

   CLOSE(112)

   IF(.NOT. gotWL) CALL FatalError("Didn't get Water Line")
   IF(.NOT. gotGrid) CALL FatalError("Didn't get Grid")
   IF(.NOT. gotSCL) CALL FatalError("Didn't get Scale")
   IF(.NOT. gotSteps) CALL FatalError("Didn't get 'No Timesteps'")
   IF(.NOT. gotName) CALL FatalError("No Run Name specified!")
   IF(.NOT. gotGeom) CALL FatalError("No Geometry File specified!")
   IF(.NOT. gotRestName .AND. REST == 1) THEN
     restname = runname
   END IF
   IF(REST == 1 .AND. gotMelange) CALL FatalError("Can't restart and load melange in same run!")

   IF(FixLat .AND. .NOT. GeomMasked) THEN
     CALL FatalError("'Fixed Lateral Margin' requires a geometry file with mask")
   END IF

   !check the geometry file exists
   INQUIRE( FILE=TRIM(geomfile), EXIST=FileExists ) 
   IF(.NOT. FileExists) CALL FatalError("Geometry input file '"//TRIM(geomfile)//"' doesn't exist!")


   IF(myid==0) THEN
     PRINT *,'--------------------Input Vars----------------------'
     WRITE(*,'(A,A)') "Run Name = ",TRIM(runname)
     IF(REST == 1) WRITE(*,'(A,A)') "Restarting from Run Name = ",TRIM(restname)
     WRITE(*,'(A,A)') "Geometry File = ",TRIM(geomfile)
     WRITE(*,'(A,A)') "Work Directory = ",TRIM(wrkdir)
     WRITE(*,'(A,A)') "Results Directory = ",TRIM(resdir)
     WRITE(*,'(A,L)') "Geometry File Has Mask = ",GeomMasked
     WRITE(*,'(A,F9.2)') "Backwall Pressure = ",PRESS
     WRITE(*,'(A,F9.2)') "Submarine Melt = ",MELT
     WRITE(*,'(A,F9.2)') "UC = ",UC
     WRITE(*,'(A,ES12.5)') "Timestep = ",DT
     WRITE(*,'(A,F9.2)') "Width = ",S
     WRITE(*,'(A,F9.2)') "Gravity = ",GRAV
     WRITE(*,'(A,F7.2)') "Density = ",RHO
     WRITE(*,'(A,F7.2)') "Water Density = ",RHOW
     WRITE(*,'(A,ES12.5)') "Youngs Modulus = ",EF0
     WRITE(*,'(A,I0)') "Size = ",LS
     WRITE(*,'(A,F9.2)') "Domain Inclination = ",SUB
     WRITE(*,'(A,F7.2)') "Grounding Line = ",GL
     WRITE(*,'(A,L)') "Do Shear Line = ",doShearLine
     WRITE(*,'(A,F7.2)') "Shear Line = ",SLIN
     WRITE(*,'(A,ES12.5)') "Max Load = ",MLOAD
     WRITE(*,'(A,ES12.5)') "Friction Scale = ",FRIC
     WRITE(*,'(A,I0)') "Restart = ",REST
     WRITE(*,'(A,F9.2)') "Porosity = ",POR
     WRITE(*,'(A,I0)') "Random Seed = ",SEEDI
     WRITE(*,'(A,ES12.5)') "Translational Damping = ",DAMP1
     WRITE(*,'(A,ES12.5)') "Rotational Damping = ",DAMP2
     WRITE(*,'(A,ES12.5)') "Drag Coefficient = ",DRAG
     WRITE(*,'(A,ES12.5)') "Bed Stiffness Constant = ",BedIntConst
     WRITE(*,'(A,L)') "Bed Z Only = ",BedZOnly
     WRITE(*,'(A,I0)') "Output Interval = ",OUTINT
     WRITE(*,'(A,I0)') "Restart Output Interval = ",RESOUTINT
     WRITE(*,'(A,ES12.5)') "Maximum Displacement = ",MAXUT
     WRITE(*,'(A,F9.2)') "Scale = ",SCL
     WRITE(*,'(A,F9.2)') "Water Line = ",WL
     WRITE(*,'(A,I0)') "No Timesteps = ",STEPS0
     WRITE(*,'(A,F9.2)') "Grid = ",GRID
     WRITE(*,'(A,F9.2)') "Fracture After Time = ",fractime
     WRITE(*,'(A,L)') "Double Precision Output = ",DoublePrec
     WRITE(*,'(A,L)') "Strict Domain Interpolation = ",StrictDomain
     WRITE(*,'(A,L)') "Fixed Lateral Margins = ",FixLat
     WRITE(*,'(A,L)') "Fixed Inflow Margin = ",FixBack
     PRINT *,'----------------------------------------------------'
   END IF
END SUBROUTINE ReadInput

SUBROUTINE BinaryVTKOutput(NRY,resdir,runname,PNN,NRXF,UT,&
     UTM,NANS,NTOT,NANPart,DoublePrec)

  USE MPI
  INCLUDE 'na90.dat'

  INTEGER :: NRY,PNN(:)
  CHARACTER(LEN=256) :: resdir, runname
  INTEGER :: NTOT, NANS(2,NTOT),NANPart(NTOT)
  TYPE(NRXF_t) :: NRXF
  TYPE(UT_t) :: UT, UTM
  LOGICAL :: DoublePrec
  !----------------------------------
  INTEGER :: NN,NNTot,NBeamsTot,counter,ms_counter,VTK_Offset
  INTEGER :: i,j,GlobalNNOffset(ntasks)
  REAL(KIND=dp) :: X,Y,Z
  CHARACTER(LEN=1024) :: output_str,datatype_str
  CHARACTER :: lfeed
  REAL(KIND=dp), ALLOCATABLE :: work_real_dp(:), displacements(:)
  REAL(KIND=sp), ALLOCATABLE :: work_real_sp(:)
  INTEGER :: fh,ierr,testsum,contig_type,realsize,intsize
  INTEGER :: Nbeams,PNbeams(ntasks),ntotal,mybeamoffset,otherbeamoffset,othertask
  INTEGER(kind=MPI_Offset_kind) :: fh_mpi_offset,fh_mpi_byte_offset, fh_starts(4), fh_mystarts(4)
  INTEGER, ALLOCATABLE :: work_int(:)
  LOGICAL :: OutputBeams,OutputDisplacement,OutputPartition

  lfeed = CHAR(10) !line feed character
  OutputDisplacement = .TRUE.
  OutputPartition = .TRUE.

  !Some MPI setup - define types and sizes
  IF(DoublePrec) THEN
    CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, realsize, ierr)
    CALL MPI_Type_Contiguous(3, MPI_DOUBLE_PRECISION, contig_type, ierr)
    datatype_str = "Float64"
  ELSE
    CALL MPI_TYPE_SIZE(MPI_REAL4, realsize, ierr)
    CALL MPI_Type_Contiguous(3, MPI_REAL4, contig_type, ierr)
    datatype_str = "Float32"
  END IF
  CALL MPI_Type_Commit(contig_type, ierr)
  CALL MPI_TYPE_SIZE(MPI_INTEGER, intsize, ierr)

  OutputBeams = .FALSE.

  !---------- Particle Info -----
  NN = PNN(myid+1)
  NNtot = SUM(PNN(1:ntasks))

  !Compute point positions
  ALLOCATE(work_real_dp(3*NN))
  work_real_dp = 0.0
  DO i=1,NN
    work_real_dp((i-1)*3 + 1) = NRXF%M(1,i)+UT%M(6*I-5)
    work_real_dp((i-1)*3 + 2) = NRXF%M(2,i)+UT%M(6*I-4)
    work_real_dp((i-1)*3 + 3) = NRXF%M(3,i)+UT%M(6*I-3)
  END DO

  !------- Beam Info ----------
  IF(OutputBeams) THEN
    !For writing node connection info (beams) 
    !need the global particle(node) numbers
    GlobalNNOffset(1) = 0
    DO i=2,ntasks
      GlobalNNOffset(i) = GlobalNNOffset(i-1) + PNN(i-1)
    END DO

    !Cross-partition beam ownership goes to lower 'mytask'
    Nbeams = COUNT(NANPart >= myid)
    !PRINT *,myid,' debug, has ',NTOT%M,' own beams, ',Nbeams,' total.'

    CALL MPI_ALLGATHER(Nbeams, 1, MPI_INTEGER, PNBeams, &
         1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    NBeamsTot = SUM(PNBeams(1:ntasks))

    ! Write all beams to work array
    ALLOCATE(work_int(Nbeams*2))

    counter = 0
    DO j=1,nbeams
      IF(NANPart(j) < myid) CYCLE
      counter = counter + 1

      mybeamoffset = GlobalNNOffset(myid+1) !TODO -these are wrong! how do we get global node numbers?
      otherbeamoffset = GlobalNNOffset(NANPart(j)+1)

      work_int(counter*2 - 1) = NRXF%PartInfo(2,NANS(1,j)) + otherbeamoffset - 1 !vtk 0 indexes the cells
      work_int(counter*2) = NANS(2,j) + mybeamoffset - 1
    END DO
  ELSE
    NBeamsTot = 0
  END IF

  !---- particle displacement ----
  IF(OutputDisplacement) THEN
    ALLOCATE(displacements(3*NN))
    displacements = 0.0
    DO i=1,NN
      displacements((i-1)*3 + 1) = UT%M(6*I-5) - UTM%M(6*I-5)
      displacements((i-1)*3 + 2) = UT%M(6*I-4) - UTM%M(6*I-4)
      displacements((i-1)*3 + 3) = UT%M(6*I-3) - UTM%M(6*I-3)
    END DO
  END IF

  !Compute offsets (global and cpu specific)
  ms_counter = 1
  fh_starts(ms_counter)=0
  fh_mystarts(ms_counter)=0

  DO i=1,ntasks
    IF(i > myid) EXIT
    fh_mystarts(ms_counter) = fh_mystarts(ms_counter) + PNN(i)*3*realsize
  END DO
  IF(myid /= 0) fh_mystarts(ms_counter) = fh_mystarts(ms_counter) + intsize !root writes an extra int at the start

  ms_counter = ms_counter + 1
  fh_starts(ms_counter) = NNTot*3*realsize + intsize
  fh_mystarts(ms_counter) = fh_starts(ms_counter)

  IF(OutputBeams) THEN

    !connectivity (particles in beams)
    DO i=1,ntasks
      IF(i > myid) EXIT
      fh_mystarts(ms_counter) = fh_mystarts(ms_counter) + PNBeams(i)*2*intsize
    END DO
    IF(myid /= 0) fh_mystarts(ms_counter) = fh_mystarts(ms_counter) + intsize !root writes an extra int at the start

    ms_counter = ms_counter + 1
    fh_starts(ms_counter) = fh_starts(ms_counter-1) + (NBeamsTot*4*intsize + intsize*3) !<- extra offset because root writes offsets & types
    fh_mystarts(ms_counter) = fh_starts(ms_counter)
  END IF

  IF(OutputDisplacement) THEN
    DO i=1,ntasks
      IF(i > myid) EXIT
      fh_mystarts(ms_counter) = fh_mystarts(ms_counter) + PNN(i)*3*realsize
    END DO
    IF(myid /= 0) fh_mystarts(ms_counter) = fh_mystarts(ms_counter) + intsize !root writes an extra int at the start

    ms_counter = ms_counter + 1
    fh_starts(ms_counter) = fh_starts(ms_counter-1) + NNTot*3*realsize + intsize
    fh_mystarts(ms_counter) = fh_starts(ms_counter)
  END IF

  IF(OutputPartition) THEN
    DO i=1,ntasks
      IF(i > myid) EXIT
      fh_mystarts(ms_counter) = fh_mystarts(ms_counter) + PNN(i)*intsize
    END DO
    IF(myid /= 0) fh_mystarts(ms_counter) = fh_mystarts(ms_counter) + intsize !root writes an extra int at the start

    ms_counter = ms_counter + 1
    fh_starts(ms_counter) = fh_starts(ms_counter-1) + NNTot*intsize + intsize
    fh_mystarts(ms_counter) = fh_starts(ms_counter)
  END IF

  CALL MPI_File_Open(MPI_COMM_WORLD,TRIM(resdir)//'/'//TRIM(runname)//'_JYR'//na(NRY)//'.vtu',&
       MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

  IF(myid==0) THEN

    VTK_Offset = 0

    !TODO - test endianness

    WRITE( output_str,'(A)') '<?xml version="1.0"?>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
   
    WRITE( output_str, '(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    WRITE( output_str,'(A)') '  <UnstructuredGrid>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    WRITE( output_str,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="',NNtot,'" NumberOfCells="',NBeamsTot,'">'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    !--------- POINTS ------------------

    WRITE( output_str,'(A)') '      <Points>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    WRITE( output_str,'(A,A,A,I0,A)') '        <DataArray type="',TRIM(datatype_str),'" Name="Position" &
         &NumberOfComponents="3" format="appended" offset="',VTK_Offset,'"/>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    VTK_Offset = VTK_Offset + NNTot*3*realsize + intsize

    WRITE( output_str,'(A)') '      </Points>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    !--------- CELLS (beams)------------------

    WRITE( output_str,'(A)') '      <Cells>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    WRITE( output_str,'(A,I0,A)') '        <DataArray type="Int32" Name="connectivity" &
         &format="appended" offset="',VTK_Offset,'"/>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
    IF(OutputBeams) VTK_Offset = VTK_Offset + NBeamsTot*2*intsize + intsize

    WRITE( output_str,'(A,I0,A)') '        <DataArray type="Int32" Name="offsets" &
         &format="appended" offset="',VTK_Offset,'"/>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
    IF(OutputBeams) VTK_Offset = VTK_Offset + NBeamsTot*intsize + intsize

    WRITE( output_str,'(A,I0,A)') '        <DataArray type="Int32" Name="types" &
         &format="appended" offset="',VTK_Offset,'"/>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
    IF(OutputBeams) VTK_Offset = VTK_Offset + NBeamsTot*intsize + intsize

    WRITE( output_str,'(A)') '      </Cells>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    !--------- POINT DATA ------------------

    WRITE( output_str,'(A)') '      <PointData>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    IF(OutputDisplacement) THEN
      WRITE( output_str,'(A,A,A,I0,A)') '        <DataArray type="',TRIM(datatype_str),'" Name="Displacement" &
           &NumberOfComponents="3" format="appended" offset="',VTK_Offset,'"/>'//lfeed
      CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

      VTK_Offset = VTK_Offset + NNTot*3*realsize + intsize
    END IF
    IF(OutputPartition) THEN
      WRITE( output_str,'(A,I0,A)') '        <DataArray type="Int32" Name="Partition" &
           &NumberOfComponents="1" format="appended" offset="',VTK_Offset,'"/>'//lfeed
      CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

      VTK_Offset = VTK_Offset + NNTot*intsize + intsize
    END IF

    WRITE( output_str,'(A)') '      </PointData>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    !----------- XML footer ----------------

    WRITE( output_str,'(A)') '    </Piece>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    WRITE( output_str,'(A)') '  </UnstructuredGrid>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    WRITE( output_str,'(A)') '  <AppendedData encoding="raw">'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    CALL MPI_File_Write(fh, "_", 1, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    !Compute the length of the header...
    CALL MPI_File_Get_Position(fh, fh_mpi_offset,ierr)
    CALL MPI_File_Get_Byte_Offset(fh, fh_mpi_offset,fh_mpi_byte_offset,ierr)
  END IF

  !... and tell everyone else
  CALL MPI_BCast(fh_mpi_byte_offset,1, MPI_OFFSET, 0, MPI_COMM_WORLD, ierr)

  !Update cpu specific start points w/ header length
  fh_starts = fh_mpi_byte_offset + fh_starts
  fh_mystarts = fh_mpi_byte_offset + fh_mystarts
  ms_counter = 1

  !Write the points (using collective I/O)
  IF(DoublePrec) THEN
    CALL MPI_File_Set_View(fh, fh_mystarts(ms_counter), MPI_DOUBLE_PRECISION, contig_type, 'native', MPI_INFO_NULL, ierr)
    IF(myid==0) CALL MPI_File_Write(fh, INT(NNtot * KIND(work_real_dp) * 3), 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
    CALL MPI_File_Write_All(fh, work_real_dp, NN, contig_type, MPI_STATUS_IGNORE, ierr)
  ELSE
    ALLOCATE(work_real_sp(SIZE(work_real_dp)))

    work_real_sp = work_real_dp
    CALL MPI_File_Set_View(fh, fh_mystarts(ms_counter), MPI_REAL4, contig_type, 'native', MPI_INFO_NULL, ierr)
    IF(myid==0) CALL MPI_File_Write(fh, INT(NNtot * KIND(work_real_sp) * 3), 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
    CALL MPI_File_Write_All(fh, work_real_sp, NN, contig_type, MPI_STATUS_IGNORE, ierr)
  END IF
  ms_counter = ms_counter + 1

  IF(OutputBeams) THEN
    !Find end of file, set view, write beam node nums, offsets, types
    CALL MPI_File_Set_View(fh, fh_mystarts(ms_counter), MPI_INTEGER, MPI_INTEGER, 'native', MPI_INFO_NULL, ierr)
    !Write byte count for connectivity
    IF(myid==0) CALL MPI_File_Write(fh, INT(NBeamsTot*KIND(work_int) * 2), 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
    CALL MPI_File_Write_All(fh, work_int, NBeams*2, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)

    !Reset the MPI I/O view to default (full file, read as bytes)
    fh_mpi_offset = 0
    CALL MPI_File_Set_View(fh, fh_mpi_offset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
    !Write beam offsets & types
    IF(myid==0) THEN
      CALL MPI_File_Seek(fh, fh_mpi_offset,MPI_SEEK_END,ierr)
      CALL MPI_File_Write(fh,NBeamsTot*KIND(work_int),1,MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      CALL MPI_File_Write(fh, (/(i*2,i=1,NBeamsTot)/),NBeamsTot,MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      CALL MPI_File_Write(fh,NBeamsTot*KIND(work_int),1,MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      CALL MPI_File_Write(fh, (/(3,i=0,NBeamsTot-1)/),NBeamsTot ,MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
    END IF
    ms_counter = ms_counter + 1
  END IF

  IF(OutputDisplacement) THEN

    !Find end of file, set view, write displacements
    IF(DoublePrec) THEN
      CALL MPI_File_Set_View(fh, fh_mystarts(ms_counter), MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
      !Write byte count for connectivity
      IF(myid==0) CALL MPI_File_Write(fh, INT(NNTot*KIND(work_real_dp) * 3), 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      CALL MPI_File_Write_All(fh, displacements, NN*3, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
    ELSE

      IF(ALLOCATED(work_real_sp)) DEALLOCATE(work_real_sp)
      ALLOCATE(work_real_sp(SIZE(displacements)))
      work_real_sp = displacements

      CALL MPI_File_Set_View(fh, fh_mystarts(ms_counter), MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)
      !Write byte count for connectivity
      IF(myid==0) CALL MPI_File_Write(fh, INT(NNTot*KIND(work_real_sp) * 3), 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      CALL MPI_File_Write_All(fh, work_real_sp, NN*3, MPI_REAL4, MPI_STATUS_IGNORE, ierr)
    END IF
    ms_counter = ms_counter + 1
  END IF

  IF(OutputPartition) THEN

    IF(ALLOCATED(work_int)) DEALLOCATE(work_int)
    ALLOCATE(work_int(NN))
    work_int = myid

    CALL MPI_File_Set_View(fh, fh_mystarts(ms_counter), MPI_INTEGER, MPI_INTEGER, 'native', MPI_INFO_NULL, ierr)
    !Write byte count for connectivity
    IF(myid==0) CALL MPI_File_Write(fh, INT(NNTot*KIND(work_int)), 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
    CALL MPI_File_Write_All(fh, work_int, NN, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)

    ms_counter = ms_counter + 1
  END IF

  !---- Writing VTU Footer -----


  !Reset the MPI I/O view to default (full file, read as bytes)
  fh_mpi_offset = 0
  CALL MPI_File_Set_View(fh, fh_mpi_offset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)

  IF(myid==0) THEN
    !Write vtu footer to end of file
    CALL MPI_File_Seek(fh, fh_mpi_offset,MPI_SEEK_END,ierr)
    WRITE( output_str,'(A)') lfeed//'  </AppendedData>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
    WRITE( output_str,'(A)') '</VTKFile>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
  END IF

  CALL MPI_File_Close(fh, ierr)


  DEALLOCATE(work_real_dp)
  IF(.NOT. DoublePrec) DEALLOCATE(work_real_sp)

END SUBROUTINE BinaryVTKOutput

SUBROUTINE BinarySTROutput(NRY,resdir,runname,NRXF,UT,&
     NANS,NTOT,NANPart,DoublePrec)

  USE MPI
  INCLUDE 'na90.dat'

  INTEGER :: NRY
  CHARACTER(LEN=256) :: resdir, runname
  INTEGER :: NTOT, NANPart(NTOT), NANS(2,NTOT)
  TYPE(UT_t), TARGET :: UT
  TYPE(NRXF_t), TARGET :: NRXF
  LOGICAL :: DoublePrec
  !----------------------------------
  INTEGER :: Nbeams,PNbeams(ntasks),NBeamsTot,counter
  INTEGER :: i,j,N1,N2
  CHARACTER(LEN=1024) :: output_str
  CHARACTER :: lfeed
  REAL(KIND=dp) :: X,Y,Z,DDX,DDY,DDZ,DX,DY,DZ,DL,L,STR
  REAL(KIND=dp), ALLOCATABLE :: work_real(:)
  REAL(KIND=sp), ALLOCATABLE :: work_real_sp(:)
  INTEGER :: fh,ierr,realsize
  INTEGER :: ntotal,othertask
  INTEGER(kind=MPI_Offset_kind) :: fh_header_offset,fh_mystart

  lfeed = CHAR(10) !line feed character

  IF(DoublePrec) THEN
    CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, realsize, ierr)
  ELSE
    CALL MPI_TYPE_SIZE(MPI_REAL4, realsize, ierr)
  END IF

  !Cross-partition beam ownership goes to lower partition no
  Nbeams = COUNT(NANPart >= myid)

  IF(DebugMode) PRINT *,myid,' debug nbeams, ntot: ',nbeams, ntot

  CALL MPI_ALLGATHER(Nbeams, 1, MPI_INTEGER, PNBeams, &
       1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
  NBeamsTot = SUM(PNBeams(1:ntasks))

  fh_mystart = 0
  DO i=1,ntasks
    IF(i > myid) EXIT
    fh_mystart = fh_mystart + PNBeams(i)*4*realsize
  END DO

  ! Write all beams to work array
  ALLOCATE(work_real(Nbeams*4))

  counter = 0
  DO j=1,NTOT
    IF(NANPart(j) < myid) CYCLE
    counter = counter + 1

    N1 = NANS(1,j)
    N2 = NANS(2,j)

    X=(NRXF%A(1,N1)+UT%A(6*N1-5)+NRXF%A(1,N2)+UT%A(6*N2-5))/2.0
    Y=(NRXF%A(2,N1)+UT%A(6*N1-4)+NRXF%A(2,N2)+UT%A(6*N2-4))/2.0
    Z=(NRXF%A(3,N1)+UT%A(6*N1-3)+NRXF%A(3,N2)+UT%A(6*N2-3))/2.0

    DDX=NRXF%A(1,N1)+UT%A(6*N1-5)-NRXF%A(1,N2)-UT%A(6*N2-5)
    DDY=NRXF%A(2,N1)+UT%A(6*N1-4)-NRXF%A(2,N2)-UT%A(6*N2-4)
    DDZ=NRXF%A(3,N1)+UT%A(6*N1-3)-NRXF%A(3,N2)-UT%A(6*N2-3)

    DX=NRXF%A(1,N1)-NRXF%A(1,N2)
    DY=NRXF%A(2,N1)-NRXF%A(2,N2)
    DZ=NRXF%A(3,N1)-NRXF%A(3,N2)
    L=SQRT(DX**2+DY**2+DZ**2)
    DL=SQRT(DDX**2+DDY**2+DDZ**2)
    STR=(DL-L)/L

    work_real(counter*4 - 3) = X
    work_real(counter*4 - 2) = Y
    work_real(counter*4 - 1) = Z
    work_real(counter*4 - 0) = STR

  END DO

  IF(counter /= nbeams) CALL FatalError("BinaryStrOutput: Programming error - wrong beam count")

  CALL MPI_File_Open(MPI_COMM_WORLD,TRIM(resdir)//'/'//TRIM(runname)//'_STR'//na(NRY)//'.bin',&
       MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

  !root write out header info (count & type)
  !everyone writes the string to get its size for offset...
  IF(DoublePrec) THEN
    WRITE( output_str,'(A,I0,A)') 'Count: ',NBeamsTot,' Type: Float64'//lfeed
  ELSE
    WRITE( output_str,'(A,I0,A)') 'Count: ',NBeamsTot,' Type: Float32'//lfeed
  END IF

  fh_header_offset = LEN_TRIM(output_str) !fortran characters are 1 byte
  fh_mystart = fh_mystart + fh_header_offset

  !... but only root writes it
  IF(myid==0) THEN
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), &
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
  END IF

  !Write the points (using collective I/O)
  IF(DoublePrec) THEN
    CALL MPI_File_Set_View(fh, fh_mystart, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
         'native', MPI_INFO_NULL, ierr)
    CALL MPI_File_Write_All(fh, work_real, Nbeams*4, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
  ELSE
    CALL MPI_File_Set_View(fh, fh_mystart, MPI_REAL4, MPI_REAL4, &
         'native', MPI_INFO_NULL, ierr)
    ALLOCATE(work_real_sp(NBeams*4))
    work_real_sp = work_real
    CALL MPI_File_Write_All(fh, work_real_sp, NBeams*4, MPI_REAL4, MPI_STATUS_IGNORE, ierr)
    DEALLOCATE(work_real_sp)
  END IF

  CALL MPI_File_Close(fh, ierr)

END SUBROUTINE BinarySTROutput

 FUNCTION ToLowerCase(from) RESULT(to)
   !------------------------------------------------------------------------------
      CHARACTER(LEN=256)  :: from
      CHARACTER(LEN=256) :: to
!------------------------------------------------------------------------------
      INTEGER :: n
      INTEGER :: i,j,nlen
      INTEGER, PARAMETER :: A=ICHAR('A'),Z=ICHAR('Z'),U2L=ICHAR('a')-ICHAR('A')

      n = LEN(to)
      DO i=LEN(from),1,-1
        IF ( from(i:i) /= ' ' ) EXIT
      END DO
      IF ( n>i ) THEN
        to(i+1:n) = ' '
        n=i
      END IF

      nlen = n
      DO i=1,nlen
        j = ICHAR( from(i:i) )
        IF ( j >= A .AND. j <= Z ) THEN
          to(i:i) = CHAR(j+U2L)
        ELSE
          to(i:i) = from(i:i)
          IF ( to(i:i)=='[') n=i-1
        END IF
      END DO

    END FUNCTION ToLowerCase

    SUBROUTINE FatalError(Message)
      CHARACTER(*) Message

      PRINT *, 'Fatal Error: ',Message
      STOP

    END SUBROUTINE FatalError

    SUBROUTINE Warn(Message)
      CHARACTER(*) Message

      PRINT *, 'WARNING: ',Message

    END SUBROUTINE Warn

  END MODULE INOUT
