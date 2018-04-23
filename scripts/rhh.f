      PROGRAM STAND

      REAL X1,Y1,Z1,X2,Y2,Z2,F,MINX,MINY,MINZ,MAX
      REAL STR1,SSTR1
      REAL MAXX,MAXY,MAXZ,PI,ZG(0:1000,0:1000)
      REAL SZ1,SZ2,SF,X,Y,SX1,SY1,VG(0:1000,0:1000)
      INTEGER N,T,I,L,IT,J,SCC,SC,NG(0:1000,0:1000)
      CHARACTER*20 INFI1,INFI2

 2    FORMAT (A)
 12   FORMAT(3F18.7)
 13   FORMAT(4F18.7)

      WRITE(*,*) 'INFI1'
      READ(*,2) INFI1
      WRITE(*,*) 'INFI2'
      READ(*,2) INFI2
      WRITE(*,*) 'N'
      READ(*,*) N
      N=N/2
      WRITE(*,*) 'MAX'
      READ(*,*) MAX

      OPEN(UNIT=700,FILE='str.dat',STATUS='UNKNOWN')
      OPEN(UNIT=710,FILE='str2.dat',STATUS='UNKNOWN')
      OPEN(UNIT=812,FILE='NG.dat',STATUS='UNKNOWN')
      OPEN(UNIT=800,FILE=INFI1,STATUS='UNKNOWN')
      OPEN(UNIT=900,FILE=INFI2,STATUS='UNKNOWN')


      DO J=0,1000
      DO I=0,1000
      VG(I,J)=0.0
      ZG(I,J)=0.0
      NG(I,J)=0
      ENDDO
      ENDDO

      DO T=1,N
      READ(800,*) X1,Y1,Z1
      READ(800,*) STR1
      READ(900,*) SX1,SY1,SZ1
      READ(900,*) SSTR1
      VEL=ABS(STR1-SSTR1)
!      WRITE(700,12) X1,Y1,VEL
      IX=INT((X1+SX1)/80.0)
      IY=INT((Y1+SY1)/80.0)
!       IF ((Z1+SZ1)/2.0.GT.1262.0) THEN
       IF ((Z1+SZ1)/2.0.GT.ZG(IX,IY)) ZG(IX,IY)=(Z1+SZ1)/2.0
       VG(IX,IY)=VG(IX,IY)+VEL
       NG(IX,IY)=NG(IX,IY)+1

       IF (VEL.GT.MAX/2.0) THEN
       WRITE(710,12) (X1+SX1)/2.0,(Y1+SY1)/2.0,VEL
       ENDIF

!       ENDIF
      ENDDO

      DO I=0,600
      DO J=0,600
      IF (NG(I,J).GE.3) THEN
!      write(812,*) NG(I,J)
       IF (VG(I,J)/(1.0*NG(I,J)).LT.MAX) THEN
       WRITE(700,13) I*40.0,J*40.0,ZG(I,J),VG(I,J)/(1.0*NG(I,J))
       ELSE
       WRITE(700,13) I*40.0,J*40.0,ZG(I,J),MAX
       ENDIF
      ELSE
      WRITE(700,13) I*40.0,J*40.0,0.0,MAX
      ENDIF
      ENDDO
!      WRITE(700,12)
      ENDDO
      

      STOP
      END



















