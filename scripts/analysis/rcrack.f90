      PROGRAM STAND

      include 'na90.dat'
      REAL AFIL(40,2000),NFIL(40,2000),DX,DY,DZ,LJ,DR,DT
      REAL FDY(2,10000),R
      REAL X,Y,Z1,Z2,LIM,STR
      INTEGER NC,T,I,L,IT,J
      CHARACTER*20 INFI

 2    FORMAT (A)
 12   FORMAT(4F18.7)

      WRITE(*,*) 'INFI?'
      READ(*,2) INFI
      WRITE(*,*) 'N?'
      READ(*,*) N
      N=N/2
      WRITE(*,*) 'LIM?'
      READ(*,*) LIM


      OPEN(UNIT=220,FILE=INFI,STATUS='UNKNOWN')
      OPEN(UNIT=230,FILE='crev.dat',STATUS='UNKNOWN')
      OPEN(UNIT=240,FILE='crev2.dat',STATUS='UNKNOWN')

      DO I=1,N
      READ(220,*) X,Y,Z
      READ(220,*) STR
      R=SQRT(X**2+Y**2)
      IF (STR.GT.LIM.AND.STR.LT.3.0.AND.Z.GT.500.0-0.015*R) WRITE(230,*) X,Y,Z,STR
      IF (STR.GT.LIM.AND.STR.LT.3.0.AND.Z.GT.500.0-0.015*R) WRITE(240,*) Y+682280.0,8855640.0-X,Z,STR
      ENDDO


      STOP
      END



















