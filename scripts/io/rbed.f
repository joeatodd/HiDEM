      PROGRAM STAND

      Real surf(-100:2000,-100:2000),bed(-100:2000,-100:2000)

      REAL X1,Y1,Z1,X2,Y2,Z2,MINX,MINY,MINZ
      REAL MAXX,MAXY,MAXZ,PI,F,MFIL(5,1000000)
      REAL SZ1,SZ2,SF,X,Y,S1,S2,B,Z,SIT
      INTEGER N,T,I,L,IT,J,SCC,SC,IX1,IY1,N2,xk,yk

 2    FORMAT (A)
 12   FORMAT(6F22.3)

      DO i=-100,2000
      DO j=-100,2000
      bed(i,j)=0.0
      ENDDO
      ENDDO


      Open(400,file='mass3.dat',STATUS='OLD')
      Open(500,file='tbed.csv',STATUS='UNKNOWN')

      READ(400,*) N2
      DO I=1,N2
      READ(400,*) x,y,s1,b1,b2,z1
      xk=INT(x/20.0)
      yk=INT(y/20.0)
      bed(xk,yk)=b2
      ENDDO
      close(400)


      DO i=0,500
      DO j=0,1000
      WRITE(500,*) i*20.0-2000.0,j*20.0-7000.0,bed(i,j)
      ENDDO
      ENDDO

      STOP
      END



















