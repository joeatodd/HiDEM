      PROGRAM STAND

      REAL*8 AFIL(40,2000),NFIL(40,2000),DX,DY,DZ,LJ,DR,DT
      REAL*8 FDY(2,10000)
      REAL*8 X,Y,B1,B2,S2,OB1
      INTEGER NC,T,I,L,IT,J
      CHARACTER*8 junk

 2    FORMAT (A)
 12   FORMAT(3F18.7)

      OPEN(UNIT=220,FILE='mass3.dat',STATUS='UNKNOWN')
      OPEN(UNIT=250,FILE='grad.dat',STATUS='UNKNOWN')

      READ(220,*) NC
      OB1=0.0
      DO I=1,NC
      READ(220,*) X,Y,S1,B1,B2
      IF (ABS(B1-S1).GT.1.0) WRITE(250,12) X,Y,B1-OB1
      OB1=B1
      ENDDO


      STOP
      END



















