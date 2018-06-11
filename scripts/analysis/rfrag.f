      PROGRAM STAND

      REAL X1,Y1,Z1,X2,Y2,Z2,F,MINX,MINY,MINZ
      REAL MAXX,MAXY,MAXZ,PI,F
      REAL SZ1,SZ2,SF,X,Y,SCL,CCAL1,CCAL2
      INTEGER N,T,I,NN
      CHARACTER*20 infil

 2    FORMAT (A)
 12   FORMAT(2I8,4F18.7)

      WRITE(*,*) 'Infil?'
      READ(*,2) infil
 
      OPEN(UNIT=800,FILE=infil,STATUS='UNKNOWN')

      WRITE(*,*) 'N?'
      READ(*,*) N

      WRITE(*,*) 'SCL?'
      READ(*,*) SCL

      CCAL1=0.0
      CCAL2=0.0
      DO T=1,N-1
      READ(800,*) I,NN
      CCAL1=CCAL1+I*NN*SCL**3.0
      ENDDO
      READ(800,*) I,NN
      CCAL2=CCAL1+I*NN*SCL**3.0

      WRITE(*,*) 'Calv=',CCAL1, CCAL1/CCAL2

      STOP
      END



















