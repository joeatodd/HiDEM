      PROGRAM STAND

      REAL AFIL(40,2000),NFIL(40,2000),DX,DY,DZ,LJ,DR,DT
      REAL FDY(2,10000)
      INTEGER N,T,I,L,IT,J
      CHARACTER*20 INFI

 2    FORMAT (A)
 12   FORMAT(4F13.7)

      OPEN(UNIT=220,FILE='bed.dat',STATUS='OLD')
      OPEN(UNIT=210,FILE='surface.dat',STATUS='OLD')
      OPEN(UNIT=600,FILE='mass2.dat',STATUS='UNKNOWN')

!      DO T=1,1115
!      READ(210,*) DY,DZ
!      FDY(1,T)=0.0
!      FDY(2,T)=DY
!      ENDDO

      DO T=1,1093
      READ(220,*) DY,DZ
      FDY(1,T)=DZ+262.0
      ENDDO


!      CLOSE(210)
!      OPEN(UNIT=210,FILE='surface.dat',STATUS='OLD')


      DO T=1,1093
      READ(210,*) DY,DZ
        IF (DY.GE.9800.0) THEN
        DO J=-20,150
        WRITE(600,12) J*10.0,DY-10000.0,DZ+262,FDY(1,T)
        ENDDO
        END IF
      END DO

      DO T=1094,1114
      READ(210,*) DY,DZ
        IF (DY.GE.10000.0) THEN
        DO J=-20,150
        WRITE(600,12) J*10.0,DY-10000.0,DZ+262,FDY(1,1093-(T-1093))
        ENDDO
        END IF
      END DO

      DO T=1,50
        DO J=-20,150
        WRITE(600,12) J*10.0,11130.0-10000.0+10.0*T,
     1                FDY(1,1093-(T+(1114-1093))),
     1                FDY(1,1093-(T+(1114-1093)))
        ENDDO
      END DO


      CLOSE (210)
      CLOSE (220)
      CLOSE (600)


!      OPEN(UNIT=600,FILE='velocity.dat',STATUS='UNKNOWN')
!      OPEN(UNIT=700,FILE='vel.dat',STATUS='UNKNOWN')

!      DO T=1,71360
!      READ(600,*) DY,DZ,VX,VY
!        IF (DY.GE.10350.0.AND.VX.GT.50.0) THEN
!           IF (VX.GT.300) THEN
!           WRITE(700,*) DY-10350,DZ+120,VX
!           ELSE
!           WRITE(700,*) DY-10350,DZ+120,300.0
!           END IF
!        END IF
!      END DO

      STOP
      END



















