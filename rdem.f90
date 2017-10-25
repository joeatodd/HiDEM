      PROGRAM STAND

      include 'na90.dat'
      REAL AFIL(40,2000),NFIL(40,2000),DX,DY,DZ,LJ,DR,DT
      REAL FDY(2,10000),MAX(0:10000,0:10000)
      REAL X,Y,Z1,Z,V
      INTEGER N,NC,T,I,L,IT,J,IX,IY
      CHARACTER*20 INFI

 2    FORMAT (A)
 12   FORMAT(3F18.7)

      WRITE(*,*) 'INFI?'
      READ(*,2) INFI

      WRITE(*,*) 'N?'
      READ(*,*) N

      OPEN(UNIT=220,FILE=INFI,STATUS='OLD')
      OPEN(UNIT=230,FILE='dem.csv',STATUS='UNKNOWN')

!      DO I=1,2023006
!      READ(220,*) X,Y,Z
!      PART(1,I)=X
!      PART(2,I)=Y
!      PART(3,I)=Z
!      ENDDO

      
      DO I=0,7000
      DO J=0,7000
      MAX(I,J)=0.0
      ENDDO
      ENDDO

      DO I=1,N
      READ(220,*) X,Y,Z
      IX=INT(X/41.0)
      IY=INT(Y/41.0)
      IF (Z.GT.MAX(IX,IY)) MAX(IX,IY)=Z
      ENDDO

      DO I=0,7000
      DO J=0,7000
      IF (MAX(I,J).GT.10.0) WRITE(230,12) I*41.0,J*41.0,MAX(I,J)
      ENDDO
      ENDDO

      STOP
      END















