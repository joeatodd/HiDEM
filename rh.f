      PROGRAM STAND

      REAL X1,Y1,Z1,X2,Y2,Z2,F,MINX,MINY,MINZ,MAX
      REAL STR1,SSTR1
      REAL MAXX,MAXY,MAXZ,PI
      REAL SZ1,SZ2,SF,X,Y,SX1,SY1,VG(-6000:6000,-6000:6000)
      INTEGER N,T,I,L,IT,J,SCC,SC,NG(-6000:6000,-6000:6000)
      CHARACTER*20 INFI1,INFI2

 2    FORMAT (A)
 12   FORMAT(3F18.7)

      WRITE(*,*) 'INFI1'
      READ(*,2) INFI1
      WRITE(*,*) 'INFI2'
      READ(*,2) INFI2
      WRITE(*,*) 'N'
      READ(*,*) N
      N=N/2
      WRITE(*,*) 'MAX'
      READ(*,*) MAX

      OPEN(UNIT=700,FILE='str.csv',STATUS='UNKNOWN')
      OPEN(UNIT=710,FILE='str2.dat',STATUS='UNKNOWN')
      OPEN(UNIT=812,FILE='NG.dat',STATUS='UNKNOWN')
      OPEN(UNIT=800,FILE=INFI1,STATUS='UNKNOWN')
      OPEN(UNIT=900,FILE=INFI2,STATUS='UNKNOWN')


      DO J=-6000,6000
      DO I=-6000,6000
      VG(I,J)=0.0
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
!      IX=INT((X1+SX1)/50.0)
!      IY=INT((Y1+SY1)/50.0)
      IX=INT((Y1+SY1)/60.0)
      IY=INT((Z1+SZ1)/60.0)
!       IF ((X1+SX1)/2.0.GT.4500.0) THEN

       VG(IX,IY)=VG(IX,IY)+VEL
       NG(IX,IY)=NG(IX,IY)+1

       VG(IX,IY+1)=VG(IX,IY+1)+VEL
       NG(IX,IY+1)=NG(IX,IY+1)+1

       VG(IX,IY-1)=VG(IX,IY-1)+VEL
       NG(IX,IY-1)=NG(IX,IY-1)+1
!--------------------------------------
       VG(IX+1,IY)=VG(IX+1,IY)+VEL
       NG(IX+1,IY)=NG(IX+1,IY)+1

       VG(IX+1,IY+1)=VG(IX+1,IY+1)+VEL
       NG(IX+1,IY+1)=NG(IX+1,IY+1)+1

       VG(IX+1,IY-1)=VG(IX+1,IY-1)+VEL
       NG(IX+1,IY-1)=NG(IX+1,IY-1)+1
!---------------------------------------
       VG(IX-1,IY)=VG(IX-1,IY)+VEL
       NG(IX-1,IY)=NG(IX-1,IY)+1

       VG(IX-1,IY+1)=VG(IX-1,IY+1)+VEL
       NG(IX-1,IY+1)=NG(IX-1,IY+1)+1

       VG(IX-1,IY-1)=VG(IX-1,IY-1)+VEL
       NG(IX-1,IY-1)=NG(IX-1,IY-1)+1

       IF (VEL.GT.MAX) THEN
       WRITE(710,12) (X1+SX1)/2.0,(Y1+SY1)/2.0,(Z1+SZ1)/2.0
       ENDIF

!       ENDIF
      ENDDO

      DO I=0,400
      DO J=0,100
      IF (NG(I,J).GE.1) THEN
!      write(812,*) NG(I,J)
       IF (VG(I,J)/(1.0*NG(I,J)).LT.MAX) THEN
       WRITE(700,12) I*30.0,J*30.0,VG(I,J)/(1.0*NG(I,J))
       ELSE
       WRITE(700,12) I*30.0,J*30.0,MAX
       ENDIF
      ELSE
      WRITE(700,12) I*30.0,J*30.0,0.5*MAX
      ENDIF
      ENDDO
      WRITE(700,12)
      ENDDO
      

      STOP
      END



















