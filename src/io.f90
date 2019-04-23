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

 SUBROUTINE ReadInput(INFILE, SimInfo)

   TYPE(SimInfo_t) :: SimInfo
   INTEGER :: readstat, i,incount
   CHARACTER(256) :: INFILE, buff, VarName,VarValue
   LOGICAL :: FileExists, gotWL=.FALSE., gotSteps=.FALSE., gotSCL=.FALSE., &
        gotGrid=.FALSE.,gotName=.FALSE.,gotGeom=.FALSE.,gotRestName=.FALSE., &
        gotViscDist=.FALSE.,gotViscStrength=.FALSE.

   OPEN(UNIT=112,FILE=infile,STATUS='old')
   incount = 0

   !Note - these variables have default values which
   !are set in their type definition (typedefs.f90)
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

     SELECT CASE (TRIM(ToLowerCase(VarName)))
     CASE ("density")
       READ(VarValue,*) SimInfo % RHO
     CASE ("water density")
       READ(VarValue,*) SimInfo % RHOW
     CASE("gravity")
       READ(VarValue,*) SimInfo % GRAV
     CASE("backwall pressure")
       READ(VarValue,*) SimInfo % PRESS
     CASE("submarine melt")
       READ(VarValue,*) SimInfo % MELT
     CASE("uc")
       READ(VarValue,*) SimInfo % UC
     CASE("timestep")
       READ(VarValue,*) SimInfo % DT
     CASE("width")
       READ(VarValue,*) SimInfo % S
     CASE("youngs modulus")
       READ(VarValue,*) SimInfo % EF0
     CASE("size")
       READ(VarValue,*) SimInfo % LS
     CASE("domain inclination")
       READ(VarValue,*) SimInfo % SUB
     CASE("water line")
       READ(VarValue,*) SimInfo % WL
       gotWL = .TRUE.
     CASE("grounding line")
       READ(VarValue,*) SimInfo % GL
     CASE("shear line")
       READ(VarValue,*) SimInfo % SLIN
       SimInfo % doShearLine = .TRUE.
     CASE("no timesteps")
       READ(VarValue,*) SimInfo % STEPS0
       gotSteps = .TRUE.
     CASE("max load")
       READ(VarValue,*) SimInfo % MLOAD
     CASE("friction scale")
       READ(VarValue,*) SimInfo % FRIC
     CASE("restart")
       READ(VarValue,*) SimInfo % REST
     CASE("scale")
       READ(VarValue,*) SimInfo % SCL
       gotSCL = .TRUE.
     CASE("grid")
       READ(VarValue,*) SimInfo % GRID
       gotGrid = .TRUE.
     CASE("porosity")
       READ(VarValue,*) SimInfo % POR
     CASE("random seed")
       READ(VarValue,*) SimInfo % SEEDI
     CASE("translational damping")
       READ(VarValue,*) SimInfo % DAMP1
     CASE("rotational damping")
       READ(VarValue,*) SimInfo % DAMP2
     CASE("viscoelastic")
       READ(VarValue,*) SimInfo % ViscoElastic
     CASE("viscous distance")
       READ(VarValue,*) SimInfo % viscdist
       gotViscDist = .TRUE.
     CASE("viscous strength")
       READ(VarValue,*) SimInfo % ViscStrength
       gotViscStrength = .TRUE.
     CASE("drag coefficient")
       READ(VarValue,*) SimInfo % DRAG_AIR
       SimInfo % DRAG_WATER = SimInfo % DRAG_AIR
     CASE("air drag coefficient")
       READ(VarValue,*) SimInfo % DRAG_AIR
     CASE("water drag coefficient")
       READ(VarValue,*) SimInfo % DRAG_WATER
     CASE("output interval")
       READ(VarValue,*) SimInfo % OUTINT
     CASE("restart output interval")
       READ(VarValue,*) SimInfo % RESOUTINT
     CASE("maximum displacement")
       READ(VarValue,*) SimInfo % MAXUT
     CASE("run name")
       READ(VarValue,*) SimInfo % runname
       gotName = .TRUE.
     CASE("restart from run name")
       READ(VarValue,*) SimInfo % restname
       gotRestName = .TRUE.
     CASE("work directory")
       READ(VarValue,*) SimInfo % wrkdir
     CASE("geometry file")
       READ(VarValue,*) SimInfo % geomfile
       gotGeom = .TRUE.
     CASE("geometry file has mask")
       READ(VarValue,*) SimInfo % GeomMasked
     CASE("results directory")
       READ(VarValue,*) SimInfo % resdir
     CASE("bed stiffness constant")
       READ(VarValue,*) SimInfo % BedIntConst
     CASE("bed damping factor")
       READ(VarValue,*) SimInfo % BedDampFactor
     CASE("bed z only")
       READ(VarValue,*) SimInfo % BedZOnly
     CASE("fracture after time")
       READ(VarValue,*) SimInfo % fractime
     CASE("strict domain interpolation")
       READ(VarValue,*) SimInfo % StrictDomain
     CASE("double precision output")
       READ(VarValue,*) SimInfo % DoublePrec
     CASE("csv output")
       READ(VarValue,*) SimInfo % CSVOutput
     CASE("fixed lateral margins")
       READ(VarValue,*) SimInfo % FixLat
     CASE("fixed inflow margin")
       READ(VarValue,*) SimInfo % FixBack
     CASE("debug mode")
       READ(VarValue,*) DebugMode
     CASE("print times")
       READ(VarValue,*) PrintTimes
     CASE("output displacement")
       READ(VarValue,*) SimInfo % outputDispl
     CASE("output rotation")
       READ(VarValue,*) SimInfo % outputRot
     CASE("output partition")
       READ(VarValue,*) SimInfo % outputPart
     CASE("melange run name")
       READ(VarValue,*) SimInfo % MelRunName
       SimInfo % gotMelange = .TRUE.
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

   IF(SimInfo % ViscoElastic) THEN
     IF(.NOT. gotViscDist) SimInfo % ViscDist = SimInfo % MLOAD
     IF(.NOT. gotViscStrength) SimInfo % ViscStrength = SimInfo%POR * SimInfo%EF0
   END IF

   IF(.NOT. gotRestName .AND. SimInfo % REST == 1) THEN
     SimInfo % restname = SimInfo % runname
   END IF
   IF(SimInfo % REST == 1 .AND. SimInfo % gotMelange) CALL FatalError("Can't restart and load melange in same run!")

   IF(SimInfo % FixLat .AND. .NOT. SimInfo % GeomMasked) THEN
     CALL FatalError("'Fixed Lateral Margin' requires a geometry file with mask")
   END IF

   !check the geometry file exists
   INQUIRE( FILE=TRIM(SimInfo % geomfile), EXIST=FileExists ) 
   IF(.NOT. FileExists) CALL FatalError("Geometry input file '"//TRIM(SimInfo % geomfile)//"' doesn't exist!")

   IF(myid==0) THEN
     PRINT *,'--------------------Input Vars----------------------'
     WRITE(*,'(A,A)') "Run Name = ",TRIM(SimInfo % runname)
     IF(SimInfo % REST == 1) WRITE(*,'(A,A)') "Restarting from Run Name = ",TRIM(SimInfo % restname)
     WRITE(*,'(A,A)') "Geometry File = ",TRIM(SimInfo % geomfile)
     WRITE(*,'(A,A)') "Work Directory = ",TRIM(SimInfo % wrkdir)
     WRITE(*,'(A,A)') "Results Directory = ",TRIM(SimInfo % resdir)
     WRITE(*,'(A,L)') "Geometry File Has Mask = ",SimInfo % GeomMasked
     WRITE(*,'(A,F9.2)') "Backwall Pressure = ", SimInfo % PRESS
     WRITE(*,'(A,F9.2)') "Submarine Melt = ", SimInfo % MELT
     WRITE(*,'(A,F9.2)') "UC = ", SimInfo % UC
     WRITE(*,'(A,ES12.5)') "Timestep = ", SimInfo % DT
     WRITE(*,'(A,F9.2)') "Width = ", SimInfo % S
     WRITE(*,'(A,F9.2)') "Gravity = ", SimInfo % GRAV
     WRITE(*,'(A,F7.2)') "Density = ", SimInfo % RHO
     WRITE(*,'(A,F7.2)') "Water Density = ", SimInfo % RHOW
     WRITE(*,'(A,ES12.5)') "Youngs Modulus = ", SimInfo % EF0
     WRITE(*,'(A,I0)') "Size = ", SimInfo % LS
     WRITE(*,'(A,F9.2)') "Domain Inclination = ", SimInfo % SUB
     WRITE(*,'(A,F7.2)') "Grounding Line = ", SimInfo % GL
     WRITE(*,'(A,L)') "Do Shear Line = ", SimInfo % doShearLine
     WRITE(*,'(A,F7.2)') "Shear Line = ", SimInfo % SLIN
     WRITE(*,'(A,ES12.5)') "Max Load = ", SimInfo % MLOAD
     WRITE(*,'(A,ES12.5)') "Friction Scale = ", SimInfo % FRIC
     WRITE(*,'(A,I0)') "Restart = ", SimInfo % REST
     WRITE(*,'(A,F9.2)') "Porosity = ", SimInfo % POR
     WRITE(*,'(A,I0)') "Random Seed = ", SimInfo % SEEDI
     WRITE(*,'(A,ES12.5)') "Translational Damping = ", SimInfo % DAMP1
     WRITE(*,'(A,ES12.5)') "Rotational Damping = ", SimInfo % DAMP2
     WRITE(*,'(A,ES12.5)') "Air Drag Coefficient = ", SimInfo % DRAG_AIR
     WRITE(*,'(A,ES12.5)') "Water Drag Coefficient = ", SimInfo % DRAG_WATER
     WRITE(*,'(A,ES12.5)') "Bed Stiffness Constant = ", SimInfo % BedIntConst
     WRITE(*,'(A,ES12.5)') "Bed Damping Factor = ", SimInfo % BedDampFactor
     WRITE(*,'(A,L)') "Bed Z Only = ", SimInfo % BedZOnly
     WRITE(*,'(A,I0)') "Output Interval = ", SimInfo % OUTINT
     WRITE(*,'(A,I0)') "Restart Output Interval = ", SimInfo % RESOUTINT
     WRITE(*,'(A,ES12.5)') "Maximum Displacement = ", SimInfo % MAXUT
     WRITE(*,'(A,F9.2)') "Scale = ", SimInfo % SCL
     WRITE(*,'(A,F9.2)') "Water Line = ", SimInfo % WL
     WRITE(*,'(A,I0)') "No Timesteps = ", SimInfo % STEPS0
     WRITE(*,'(A,F9.2)') "Grid = ", SimInfo % GRID
     WRITE(*,'(A,F9.2)') "Fracture After Time = ", SimInfo % fractime
     WRITE(*,'(A,L)') "Viscoelastic sim = ", SimInfo % ViscoElastic
     WRITE(*,'(A,F9.2)') "Viscous Distance = ", SimInfo % ViscDist
     WRITE(*,'(A,F9.2)') "Viscous Strength = ", SimInfo % ViscStrength
     WRITE(*,'(A,L)') "Double Precision Output = ", SimInfo % DoublePrec
     WRITE(*,'(A,L)') "Strict Domain Interpolation = ", SimInfo % StrictDomain
     WRITE(*,'(A,L)') "Fixed Lateral Margins = ", SimInfo % FixLat
     WRITE(*,'(A,L)') "Fixed Inflow Margin = ", SimInfo % FixBack
     WRITE(*,'(A,L)') "Output Displacement = ", SimInfo % outputDispl
     WRITE(*,'(A,L)') "Output Rotation = ", SimInfo % outputRot
     WRITE(*,'(A,L)') "Output Partition = ", SimInfo % outputPart
     PRINT *,'----------------------------------------------------'
   END IF
END SUBROUTINE ReadInput

SUBROUTINE BinaryVTKOutput(SI,NRY,PNN,NRXF,UT,UTM,NANS,NTOT,NANPart)

  USE MPI
  INCLUDE 'na90.dat'

  INTEGER :: NRY,PNN(:)
  INTEGER :: NTOT, NANS(2,NTOT),NANPart(NTOT)
  TYPE(NRXF_t) :: NRXF
  TYPE(UT_t) :: UT, UTM
  !----------------------------------
  INTEGER :: NN,NNTot,NBeamsTot,counter,ms_counter,VTK_Offset
  INTEGER :: i,j,GlobalNNOffset(ntasks)
  REAL(KIND=dp) :: X,Y,Z
  CHARACTER(LEN=1024) :: output_str,datatype_str
  CHARACTER :: lfeed
  REAL(KIND=dp), ALLOCATABLE :: work_real_dp(:), displacements(:), rotations(:)
  INTEGER :: fh,ierr,testsum,contig_type,realsize,intsize
  INTEGER :: Nbeams,PNbeams(ntasks),ntotal,mybeamoffset,otherbeamoffset,othertask
  INTEGER(kind=MPI_Offset_kind) :: fh_mpi_offset,fh_mpi_byte_offset, fh_starts(10), fh_mystarts(10)
  INTEGER, ALLOCATABLE :: work_int(:)
  LOGICAL :: OutputBeams,OutputDisplacement,OutputRotation,OutputPartition
  TYPE(SimInfo_t) :: SI

  lfeed = CHAR(10) !line feed character

  OutputDisplacement = SI%outputDispl
  OutputRotation = SI%outputRot
  OutputPartition = SI%outputPart
  OutputBeams = .FALSE.

  !Some MPI setup - define types and sizes
  IF(SI%DoublePrec) THEN
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
    !need the global particle(node) positions (NOT GID!)
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

  !---- particle rotation ----
  IF(OutputRotation) THEN
    ALLOCATE(rotations(3*NN))
    rotations = 0.0
    DO i=1,NN
      rotations((i-1)*3 + 1) = UT%M(6*I-2)
      rotations((i-1)*3 + 2) = UT%M(6*I-1)
      rotations((i-1)*3 + 3) = UT%M(6*I-0)
    END DO
  END IF

  !Compute offsets (global and cpu specific)
  ms_counter = 1
  fh_starts(ms_counter)=0
  fh_mystarts(ms_counter)=0

  CALL SetVTKOffsets(fh_starts, fh_mystarts, ms_counter, PNN, 3, realsize, intsize)

  IF(OutputBeams) CALL SetVTKOffsets(fh_starts, fh_mystarts, ms_counter, &
       PNBeams, 2, intsize, intsize, ((NBeamsTot*4*intsize + intsize*3)))
  IF(OutputDisplacement) CALL SetVTKOffsets(fh_starts, fh_mystarts, ms_counter, &
       PNN, 3, realsize, intsize)
  IF(OutputRotation) CALL SetVTKOffsets(fh_starts, fh_mystarts, ms_counter, PNN, &
       3, realsize, intsize)
  IF(OutputPartition) CALL SetVTKOffsets(fh_starts, fh_mystarts, ms_counter, PNN, &
       1, intsize, intsize)

  CALL MPI_File_Open(MPI_COMM_WORLD,TRIM(SI%resdir)//'/'//TRIM(SI%runname)//'_JYR'//na(NRY)//'.vtu',&
       MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

  IF(myid==0) THEN

    VTK_Offset = 0

    !TODO - test endianness

    WRITE( output_str,'(A)') '<?xml version="1.0"?>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), &
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
   
    WRITE( output_str, '(A)') '<VTKFile type="UnstructuredGrid" version=&
         &"0.1" byte_order="LittleEndian">'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), &
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    WRITE( output_str,'(A)') '  <UnstructuredGrid>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), &
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    WRITE( output_str,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="',&
         NNtot,'" NumberOfCells="',NBeamsTot,'">'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    !--------- POINTS ------------------

    WRITE( output_str,'(A)') '      <Points>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    WRITE( output_str,'(A,A,A,I0,A)') '        <DataArray type="',&
         TRIM(datatype_str),'" Name="Position" &
         &NumberOfComponents="3" format="appended" offset="',VTK_Offset,'"/>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    VTK_Offset = VTK_Offset + NNTot*3*realsize + intsize

    WRITE( output_str,'(A)') '      </Points>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    !--------- CELLS (beams)------------------

    WRITE( output_str,'(A)') '      <Cells>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    WRITE( output_str,'(A,I0,A)') '        <DataArray type="Int32" Name="connectivity" &
         &format="appended" offset="',VTK_Offset,'"/>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
    IF(OutputBeams) VTK_Offset = VTK_Offset + NBeamsTot*2*intsize + intsize

    WRITE( output_str,'(A,I0,A)') '        <DataArray type="Int32" Name="offsets" &
         &format="appended" offset="',VTK_Offset,'"/>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
    IF(OutputBeams) VTK_Offset = VTK_Offset + NBeamsTot*intsize + intsize

    WRITE( output_str,'(A,I0,A)') '        <DataArray type="Int32" Name="types" &
         &format="appended" offset="',VTK_Offset,'"/>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
    IF(OutputBeams) VTK_Offset = VTK_Offset + NBeamsTot*intsize + intsize

    WRITE( output_str,'(A)') '      </Cells>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    !--------- POINT DATA ------------------

    WRITE( output_str,'(A)') '      <PointData>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    IF(OutputDisplacement) CALL WriteVTKPointHeader(fh, "Displacement", &
         datatype_str, 3, VTK_Offset, NNTot, realsize, intsize)

    IF(OutputRotation) CALL WriteVTKPointHeader(fh, "Rotation", &
         datatype_str, 3, VTK_Offset, NNTot, realsize, intsize)

    IF(OutputPartition) CALL WriteVTKPointHeader(fh, "Partition", &
         "Int32", 1, VTK_Offset, NNTot, intsize, intsize)

    WRITE( output_str,'(A)') '      </PointData>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    !----------- XML footer ----------------

    WRITE( output_str,'(A)') '    </Piece>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    WRITE( output_str,'(A)') '  </UnstructuredGrid>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

    WRITE( output_str,'(A)') '  <AppendedData encoding="raw">'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

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

  CALL WriteRealPointDataToVTK(SI, fh, work_real_dp, fh_mystarts, ms_counter, PNN, 3)

  IF(OutputBeams) THEN
    !Find end of file, set view, write beam node nums, offsets, types
    CALL MPI_File_Set_View(fh, fh_mystarts(ms_counter), MPI_INTEGER, MPI_INTEGER,&
         'native', MPI_INFO_NULL, ierr)
    !Write byte count for connectivity
    IF(myid==0) CALL MPI_File_Write(fh, INT(NBeamsTot*KIND(work_int) * 2), 1, &
         MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
    CALL MPI_File_Write_All(fh, work_int, NBeams*2, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)

    !Reset the MPI I/O view to default (full file, read as bytes)
    fh_mpi_offset = 0
    CALL MPI_File_Set_View(fh, fh_mpi_offset, MPI_BYTE, MPI_BYTE, 'native', &
         MPI_INFO_NULL, ierr)
    !Write beam offsets & types
    IF(myid==0) THEN
      CALL MPI_File_Seek(fh, fh_mpi_offset,MPI_SEEK_END,ierr)
      CALL MPI_File_Write(fh,NBeamsTot*KIND(work_int),1,MPI_INTEGER, &
           MPI_STATUS_IGNORE, ierr)
      CALL MPI_File_Write(fh, (/(i*2,i=1,NBeamsTot)/),NBeamsTot,MPI_INTEGER, &
           MPI_STATUS_IGNORE, ierr)
      CALL MPI_File_Write(fh,NBeamsTot*KIND(work_int),1,MPI_INTEGER, &
           MPI_STATUS_IGNORE, ierr)
      CALL MPI_File_Write(fh, (/(3,i=0,NBeamsTot-1)/),NBeamsTot ,MPI_INTEGER, &
           MPI_STATUS_IGNORE, ierr)
    END IF
    ms_counter = ms_counter + 1
  END IF

  IF(OutputDisplacement) CALL WriteRealPointDataToVTK(SI, fh, displacements, fh_mystarts, &
       ms_counter, PNN, 3)
  IF(OutputRotation) CALL WriteRealPointDataToVTK(SI, fh, rotations, fh_mystarts, ms_counter, &
       PNN, 3)

  IF(OutputPartition) THEN
    IF(ALLOCATED(work_int)) DEALLOCATE(work_int)
    ALLOCATE(work_int(NN))
    work_int = myid

    CALL WriteIntPointDataToVTK(fh, work_int, fh_mystarts, ms_counter, PNN, 1)
  END IF

  !---- Writing VTU Footer -----


  !Reset the MPI I/O view to default (full file, read as bytes)
  fh_mpi_offset = 0
  CALL MPI_File_Set_View(fh, fh_mpi_offset, MPI_BYTE, MPI_BYTE, 'native', &
       MPI_INFO_NULL, ierr)

  IF(myid==0) THEN
    !Write vtu footer to end of file
    CALL MPI_File_Seek(fh, fh_mpi_offset,MPI_SEEK_END,ierr)
    WRITE( output_str,'(A)') lfeed//'  </AppendedData>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER,&
         MPI_STATUS_IGNORE, ierr)
    WRITE( output_str,'(A)') '</VTKFile>'//lfeed
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), MPI_CHARACTER,&
         MPI_STATUS_IGNORE, ierr)
  END IF

  CALL MPI_File_Close(fh, ierr)


  DEALLOCATE(work_real_dp)

END SUBROUTINE BinaryVTKOutput

SUBROUTINE SetVTKOffsets(fh_starts, fh_mystarts, ms_counter, counts, DOFs, var_size, intsize, custom_offset)

  INTEGER(kind=MPI_Offset_kind) :: fh_starts(10), fh_mystarts(10)
  INTEGER :: ms_counter, counts(:), DOFs, var_size, intsize
  INTEGER, OPTIONAL :: custom_offset
  !----------------------
  INTEGER :: i,counts_tot
  LOGICAL :: offset_override

  offset_override = PRESENT(custom_offset)
  counts_tot = SUM(counts(1:ntasks))

  DO i=1,ntasks
    IF(i > myid) EXIT
    fh_mystarts(ms_counter) = fh_mystarts(ms_counter) + counts(i)*DOFs*var_size
  END DO
  !root writes an extra int at the start
  IF(myid /= 0) fh_mystarts(ms_counter) = fh_mystarts(ms_counter) + intsize 


  ms_counter = ms_counter + 1
  IF(offset_override) THEN
    fh_starts(ms_counter) = fh_starts(ms_counter-1) + custom_offset
  ELSE
    fh_starts(ms_counter) = fh_starts(ms_counter-1) + counts_tot*DOFs*var_size + intsize
  END IF
  fh_mystarts(ms_counter) = fh_starts(ms_counter)

END SUBROUTINE SetVTKOffsets

SUBROUTINE WriteVTKPointHeader(fh, varname, datatype, DOFs, VTK_Offset, counttot, var_size, intsize)
  CHARACTER(LEN=*) :: varname, datatype
  INTEGER :: DOFs, VTK_Offset, var_size, intsize, counttot, fh
  !------------------------
  CHARACTER(LEN=1024) :: output_str
  CHARACTER :: lfeed
  INTEGER :: ierr

  lfeed= CHAR(10)

  WRITE( output_str,'(A,A,A,A,A,I0,A,I0,A)') '        <DataArray type="',&
       TRIM(datatype),'" Name="',TRIM(varname),'" NumberOfComponents="',DOFs,'" &
       &format="appended" offset="',VTK_Offset,'"/>'//lfeed
  CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str),&
       MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

  VTK_Offset = VTK_Offset + counttot*DOFs*var_size + intsize

END SUBROUTINE WriteVTKPointHeader

SUBROUTINE WriteRealPointDataToVTK(SI,fh, data_arr, fh_mystarts, ms_counter, counts, DOFs)

  REAL(KIND=dp) :: data_arr(:)
  INTEGER(kind=MPI_Offset_kind) :: fh_mystarts(10)
  INTEGER :: fh, ms_counter, DOFs, counts(:)
  TYPE(SimInfo_t) :: SI
  !-------------------------------
  INTEGER :: mycount, totcount, ierr
  REAL(KIND=sp), ALLOCATABLE :: work_real_sp(:)

  mycount = counts(myid+1)
  totcount = SUM(counts(1:ntasks))

  !Find end of file, set view, write var
  IF(SI%DoublePrec) THEN

    CALL MPI_File_Set_View(fh, fh_mystarts(ms_counter), MPI_DOUBLE_PRECISION,&
         MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
    !Root writes byte count
    IF(myid==0) CALL MPI_File_Write(fh, INT(totcount * KIND(data_arr) * DOFs), 1,&
         MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
    CALL MPI_File_Write_All(fh, data_arr, mycount*DOFs, MPI_DOUBLE_PRECISION, &
         MPI_STATUS_IGNORE, ierr)


  ELSE

    ALLOCATE(work_real_sp(SIZE(data_arr)))
    work_real_sp = REAL(data_arr,sp)

    CALL MPI_File_Set_View(fh, fh_mystarts(ms_counter), MPI_REAL4, MPI_REAL4,&
         'native', MPI_INFO_NULL, ierr)
    !Write byte count for connectivity
    IF(myid==0) CALL MPI_File_Write(fh, INT(totcount * KIND(work_real_sp) * DOFs), 1,&
         MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
    CALL MPI_File_Write_All(fh, work_real_sp, mycount*DOFs, MPI_REAL4, &
         MPI_STATUS_IGNORE, ierr)

    DEALLOCATE(work_real_sp)

  END IF

  ms_counter = ms_counter + 1

END SUBROUTINE WriteRealPointDataToVTK

SUBROUTINE WriteIntPointDataToVTK(fh, data_arr, fh_mystarts, ms_counter, counts, DOFs)

  INTEGER :: data_arr(:)
  INTEGER(kind=MPI_Offset_kind) :: fh_mystarts(10)
  INTEGER :: fh, ms_counter, DOFs, counts(:)
  !-------------------------------
  INTEGER :: mycount, totcount, ierr

  mycount = counts(myid+1)
  totcount = SUM(counts(1:ntasks))
  
  CALL MPI_File_Set_View(fh, fh_mystarts(ms_counter), MPI_INTEGER,&
       MPI_INTEGER, 'native', MPI_INFO_NULL, ierr)
  !Root writes byte count
  IF(myid==0) CALL MPI_File_Write(fh, INT(totcount * KIND(data_arr) * DOFs), 1,&
       MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
  CALL MPI_File_Write_All(fh, data_arr, mycount*DOFs, MPI_INTEGER, &
       MPI_STATUS_IGNORE, ierr)

  ms_counter = ms_counter + 1

END SUBROUTINE WriteIntPointDataToVTK

SUBROUTINE BinarySTROutput(SI,NRY,NTOT,NANPart,EFS)

  USE MPI
  INCLUDE 'na90.dat'

  INTEGER :: NRY
  INTEGER :: NTOT, NANPart(NTOT)
  TYPE(SimInfo_t) :: SI
  REAL(KIND=dp) :: EFS(:)
  !----------------------------------
  INTEGER :: Nbeams,PNbeams(ntasks),NBeamsTot,counter
  INTEGER :: i,j
  CHARACTER(LEN=1024) :: output_str
  CHARACTER :: lfeed
  REAL(KIND=dp), ALLOCATABLE :: work_real(:)
  REAL(KIND=sp), ALLOCATABLE :: work_real_sp(:)
  INTEGER :: fh,ierr,realsize
  INTEGER :: ntotal,othertask,NVars
  INTEGER(kind=MPI_Offset_kind) :: fh_header_offset,fh_mystart

  lfeed = CHAR(10) !line feed character
  NVars = 1 !EFS    - STR, XYZ obtained using .vtu file
  !Cross-partition beam ownership goes to lower partition no
  Nbeams = COUNT(NANPart >= myid)

  IF(SI%DoublePrec) THEN
    CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, realsize, ierr)
  ELSE
    CALL MPI_TYPE_SIZE(MPI_REAL4, realsize, ierr)
  END IF

  IF(DebugMode) PRINT *,myid,' debug nbeams, ntot: ',nbeams, ntot

  CALL MPI_ALLGATHER(Nbeams, 1, MPI_INTEGER, PNBeams, &
       1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
  NBeamsTot = SUM(PNBeams(1:ntasks))

  ALLOCATE(work_real(Nbeams*NVars))

  !Write our EFS values to work array
  counter = 0
  DO j=1,NTOT
    IF(NANPart(j) < myid) CYCLE
    counter = counter + 1
    work_real(counter) = EFS(j)
  END DO

  IF(counter /= nbeams) CALL FatalError("BinaryStrOutput: Programming error - wrong beam count")

  CALL MPI_File_Open(MPI_COMM_WORLD,TRIM(SI%resdir)//'/'//TRIM(SI%runname)//'_STR'//na(NRY)//'.bin',&
       MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

  !root write out header info (count & type)
  !everyone writes the string to get its size for offset...
  IF(SI%DoublePrec) THEN
    WRITE( output_str,'(A,I0,A)') 'Count: ',NBeamsTot,' Type: Float64'//lfeed
  ELSE
    WRITE( output_str,'(A,I0,A)') 'Count: ',NBeamsTot,' Type: Float32'//lfeed
  END IF

  !offsets for ints & real variables
  fh_mystart = 0

  !Var offsets
  DO i=1,ntasks
    IF(i > myid) EXIT
    fh_mystart = fh_mystart + PNBeams(i)*NVars*realsize !NVars == 1
  END DO

  !Initial string
  fh_header_offset = LEN_TRIM(output_str) !fortran characters are 1 byte

  fh_mystart = fh_mystart + fh_header_offset

  !Only root writes header
  IF(myid==0) THEN
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), &
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
  END IF

  !Write EFS (using collective I/O)
  IF(SI%DoublePrec) THEN
    CALL MPI_File_Set_View(fh, fh_mystart, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
         'native', MPI_INFO_NULL, ierr)
    CALL MPI_File_Write_All(fh, work_real, Nbeams*NVars, MPI_DOUBLE_PRECISION, &
         MPI_STATUS_IGNORE, ierr)
  ELSE
    CALL MPI_File_Set_View(fh, fh_mystart, MPI_REAL4, MPI_REAL4, &
         'native', MPI_INFO_NULL, ierr)
    ALLOCATE(work_real_sp(NBeams*NVars))
    work_real_sp = work_real
    CALL MPI_File_Write_All(fh, work_real_sp, NBeams*NVars, MPI_REAL4, MPI_STATUS_IGNORE, ierr)
    DEALLOCATE(work_real_sp)
  END IF

  CALL MPI_File_Close(fh, ierr)

END SUBROUTINE BinarySTROutput


SUBROUTINE BinarySTHOutput(SI,PNN,NRXF,NANS,NTOT,NANPart)

  USE MPI
  INCLUDE 'na90.dat'

  INTEGER :: NTOT, NANPart(NTOT), NANS(2,NTOT),PNN(:)
  TYPE(NRXF_t), TARGET :: NRXF
  TYPE(SimInfo_t) :: SI
  !----------------------------------
  INTEGER :: Nbeams,PNbeams(ntasks),NBeamsTot,counter
  INTEGER :: i,j,N1,N2,GlobalNNOffset(ntasks),otherbeamoffset,mybeamoffset
  CHARACTER(LEN=1024) :: output_str
  CHARACTER :: lfeed
  INTEGER(KIND=4), ALLOCATABLE :: work_int(:)
  INTEGER :: fh,ierr
  INTEGER(kind=MPI_Offset_kind) :: fh_header_offset,fh_mystart_int

  lfeed = CHAR(10) !line feed character
  !Cross-partition beam ownership goes to lower partition no
  Nbeams = COUNT(NANPart >= myid)

  IF(DebugMode) PRINT *,myid,' debug nbeams, ntot: ',nbeams, ntot

  CALL MPI_ALLGATHER(Nbeams, 1, MPI_INTEGER, PNBeams, &
       1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
  NBeamsTot = SUM(PNBeams(1:ntasks))

  ! Write all beams to work array
  ALLOCATE(work_int(NBeams*2))

  !Compute particle array offsets for NAN (not GID!)
  GlobalNNOffset(1) = 0
  DO i=2,ntasks
    GlobalNNOffset(i) = GlobalNNOffset(i-1) + PNN(i-1)
  END DO
  mybeamoffset = GlobalNNOffset(myid+1)

  counter = 0
  DO j=1,NTOT
    IF(NANPart(j) < myid) CYCLE
    counter = counter + 1

    N1 = NANS(1,j)
    N2 = NANS(2,j)
    otherbeamoffset = GlobalNNOffset(NANPart(j)+1)

    !Construct array of beam particle IDs
    work_int(counter*2 - 1) = NRXF%PartInfo(2,N1) + otherbeamoffset
    work_int(counter*2 - 0) = N2 + mybeamoffset
  END DO

  IF(counter /= nbeams) CALL FatalError("BinaryStrOutput: Programming error - wrong beam count")

  CALL MPI_File_Open(MPI_COMM_WORLD,TRIM(SI%resdir)//'/'//TRIM(SI%runname)//'_STH.bin',&
       MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

  !root write out header info (count & type)
  !everyone writes the string to get its size for offset...
  WRITE( output_str,'(A,I0,A,ES15.8,A,ES15.8,A,ES15.8,A)') &
       'BeamCount: ',NBeamsTot,&
       ' SCL: ',SI%SCL,&
       ' EF0: ',SI%EF0,&
       ' DT: ',SI%DT,&
       ' Type: Int32'//lfeed

  !offsets for N1, N2
  fh_mystart_int = 0

  !N1, N2 offsets
  DO i=1,ntasks
    IF(i > myid) EXIT
    fh_mystart_int = fh_mystart_int + PNBeams(i)*2*4 !4 = INT32
  END DO

  !Initial string
  fh_header_offset = LEN_TRIM(output_str) !fortran characters are 1 byte

  fh_mystart_int = fh_mystart_int + fh_header_offset

  !... Only root writes the header
  IF(myid==0) THEN
    CALL MPI_File_Write(fh, TRIM(output_str), LEN_TRIM(output_str), &
         MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
  END IF

  !Write N1, N2 (GID)
  CALL MPI_File_Set_View(fh, fh_mystart_int, MPI_INTEGER, MPI_INTEGER, &
       'native', MPI_INFO_NULL, ierr)
  CALL MPI_File_Write_All(fh, work_int, Nbeams*2, MPI_INTEGER, &
       MPI_STATUS_IGNORE, ierr)

  CALL MPI_File_Close(fh, ierr)

END SUBROUTINE BinarySTHOutput


 FUNCTION ToLowerCase(from) RESULT(to)
     !------------------------------------------------------------------------------
      CHARACTER(LEN=*)  :: from
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
