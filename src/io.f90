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

  IMPLICIT NONE

  INTEGER :: MPI_COMM_ACTIVE

CONTAINS

 SUBROUTINE ReadInput(INFILE, myid, runname, wrkdir, resdir, geomfile, PRESS, MELT, UC, DT, S, GRAV, &
      RHO, RHOW, EF0, LS, SUB, GL, SLIN, MLOAD, FRIC, REST, restname, POR, SEEDI, DAMP1, &
      DAMP2, DRAG, BedIntConst, BedZOnly, OUTINT, RESOUTINT, MAXUT, SCL, WL, STEPS0, GRID, fractime)
   REAL*8 :: PRESS, MELT, UC, DT, S, EF0, SUB, GL, SLIN, MLOAD, FRIC, POR
   REAL*8 :: DAMP1, DAMP2, DRAG,MAXUT, SCL, WL, GRID, GRAV, RHO, RHOW, BedIntConst
   REAL*8 :: fractime
   INTEGER :: REST, SEEDI, OUTINT, RESOUTINT, STEPS0, LS
   INTEGER :: myid, readstat, i,incount
   CHARACTER(256) :: INFILE, geomfile, buff,VarName,VarValue,runname,wrkdir,&
        resdir,restname
   LOGICAL :: BedZOnly
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
   MLOAD = 0.0002
   FRIC = 1.0
   REST = 0
   POR = 0.1
   SEEDI = 11695378
   DAMP1 = 1.0E4
   DAMP2 = 1.0E4
   DRAG = 1.0E1
   OUTINT = 20000
   RESOUTINT = 20000
   MAXUT = 5.0E3
   GRAV = 9.81
   RHO = 900.0
   RHOW = 1030.0
   BedIntConst = 1.0E8
   BedZOnly = .TRUE.
   wrkdir = './'
   resdir = './'
   fractime = 40.0

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
     CASE("results directory")
       READ(VarValue,*) resdir
     CASE("bed stiffness constant")
       READ(VarValue,*) BedIntConst
     CASE("bed z only")
       READ(VarValue,*) BedZOnly
     CASE("fracture after time")
       READ(VarValue,*) fractime
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

   IF(myid==0) THEN
     PRINT *,'---------------Input Vars-----------------'
     WRITE(*,'(A,A)') "Run Name = ",TRIM(runname)
     IF(REST == 1) WRITE(*,'(A,A)') "Restarting from Run Name = ",TRIM(restname)
     WRITE(*,'(A,A)') "Geometry File = ",TRIM(geomfile)
     WRITE(*,'(A,A)') "Work Directory = ",TRIM(wrkdir)
     WRITE(*,'(A,A)') "Results Directory = ",TRIM(resdir)
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
     PRINT *,'---------------End Input------------------'
   END IF
END SUBROUTINE ReadInput

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

    !Create the 'MPI_COMM_ACTIVE' communicator which includes only 
    !CPUs which are actually doing something
    SUBROUTINE RedefineMPI(noprocs,myid)

      INCLUDE 'mpif.h'

      INTEGER :: myid,i,noprocs, groupworld, groupactive,ierr
      INTEGER, ALLOCATABLE :: active_parts(:)

        CALL MPI_COMM_GROUP(MPI_COMM_WORLD, groupworld, ierr)
        ALLOCATE(active_parts(noprocs))
        DO i=1,noprocs
          active_parts(i) = i-1
        END DO

        CALL MPI_Group_Incl(groupworld, noprocs, active_parts, groupactive, ierr)
        CALL MPI_Comm_Create(MPI_COMM_WORLD, groupactive, MPI_COMM_ACTIVE, ierr)

        IF(myid < noprocs) THEN
          CALL MPI_BARRIER(MPI_COMM_ACTIVE, ierr)
        ELSE
          PRINT *,'MPI process ',myid,' terminating because not required.'
          CALL MPI_Finalize()
          STOP
        END IF

      END SUBROUTINE RedefineMPI

  END MODULE INOUT
