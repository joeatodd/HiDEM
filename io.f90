MODULE INOUT

  IMPLICIT NONE

CONTAINS

 SUBROUTINE ReadInput(INFILE, PRESS, MELT, UC, DT, S, GRAV, RHO, RHOW, EF0, LS, &
      SUB, GL, SLIN, MLOAD, FRIC, REST, POR, SEEDI, DAMP1, &
      DAMP2, DRAG, OUTINT, RESOUTINT, MAXUT, SCL, WL, STEPS0, GRID)
   REAL*8 :: PRESS, MELT, UC, DT, S, EF0, SUB, GL, SLIN, MLOAD, FRIC, POR
   REAL*8 :: DAMP1, DAMP2, DRAG,MAXUT, SCL, WL, GRID, GRAV, RHO, RHOW
   INTEGER :: REST, SEEDI, OUTINT, RESOUTINT, STEPS0, LS
   INTEGER :: readstat, i,incount
   CHARACTER(256) :: INFILE, buff,VarName,VarValue
   LOGICAL :: gotWL=.FALSE., gotSteps=.FALSE., gotSCL=.FALSE., &
        gotGrid=.FALSE.

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

   DO
     READ(112,"(A)", IOSTAT=readstat) buff
     IF(readstat > 0) STOP
     IF(readstat < 0) EXIT
     incount = incount+1

     i = INDEX(TRIM(buff),'!')
     IF(i>0) CYCLE

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

      PRINT *,Message
      STOP

    END SUBROUTINE FatalError
  END MODULE INOUT
