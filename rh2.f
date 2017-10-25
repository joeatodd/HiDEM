      PROGRAM STAND

      REAL X1,Y1,Z1,X2,Y2,Z2,MINX,MINY,MINZ,MAX,TVC
      REAL MAXX,MAXY,MAXZ,PI,F,MFIL(5,1000000),AVV,AVC
      REAL SZ1,SZ2,SF,X,Y,SX1,SY1,VG(0:1000,0:1000)
      INTEGER N,T,I,L,IT,J,SCC,SC,NG(0:1000,0:1000),TVN
      CHARACTER*20 INFI1,INFI2

 2    FORMAT (A)
 12   FORMAT(3F18.7)

      WRITE(*,*) 'INFI1'
      READ(*,2) INFI1
      WRITE(*,*) 'INFI2'
      READ(*,2) INFI2
      WRITE(*,*) 'N'
      READ(*,*) N
      WRITE(*,*) 'MAX'
      READ(*,*) MAX
!      MAX=MAX*1.0e-04

      OPEN(UNIT=700,FILE='vel.dat',STATUS='UNKNOWN')
      OPEN(UNIT=800,FILE=INFI1,STATUS='UNKNOWN')
      OPEN(UNIT=900,FILE=INFI2,STATUS='UNKNOWN')


      DO J=0,1000
      DO I=0,1000
      VG(I,J)=0.0
      NG(I,J)=0
      ENDDO
      ENDDO

      AVV=0.0
      AVC=0.0
      DO T=1,N
      READ(800,*) X1,Y1,Z1
      READ(900,*) SX1,SY1,SZ1
      VEL=SQRT((Y1-SY1)**2+(Z1-SZ1)**2)
!      WRITE(700,12) X1,Y1,VEL
!      IX=INT((X1+SX1)/200.0)
!      IY=INT((Y1+SY1)/200.0)
      IX=INT((Y1+SY1)/100.0)
      IY=INT((Z1+SZ1)/100.0)
!       IF ((Z1+SZ1)/2.0.GT.1762.0) THEN
       VG(IX,IY)=VG(IX,IY)+VEL
       NG(IX,IY)=NG(IX,IY)+1
       AVV=AVV+VEL
       AVC=AVC+1.0
!       ENDIF
      ENDDO

      WRITE(*,*) 'VEL=',AVV/AVC*3600.0*2.4
      WRITE(*,*) 'DIST=',AVV/AVC
      TVC=0.0
      TVN=0
      DO I=0,200
      DO J=0,100
      IF (NG(I,J).NE.0) THEN
       IF (VG(I,J)/(1.0*NG(I,J)).LT.MAX) THEN
!       WRITE(700,12) I*100.0,J*100.0,1.0e+04*VG(I,J)/(1.0*NG(I,J))
       WRITE(700,12) I*50.0,J*50.0,VG(I,J)/(1.0*NG(I,J))
!       TVC=TVC+1.0e+04*VG(I,J)/(1.0*NG(I,J))
       TVC=TVC+VG(I,J)/(1.0*NG(I,J))
       TVN=TVN+1
       ELSE
!       WRITE(700,12) I*100.0,J*100.0,1.0e+04*MAX
       WRITE(700,12) I*50.0,J*50.0,MAX
!       TVC=TVC+1.0e+04*MAX
       TVC=TVC+MAX
       TVN=TVN+1
       ENDIF
      ELSE
      WRITE(700,12) I*50.0,J*50.0,0.5*MAX
      ENDIF
      ENDDO
      WRITE(700,12)
      ENDDO
      

      WRITE(*,*) TVC/(1.0*TVN)
      STOP
      END



















