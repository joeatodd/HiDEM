      PROGRAM RC

      IMPLICIT NONE
      INCLUDE 'param.dat'
      INCLUDE 'na90.dat'
      INTEGER I,N1,N2,NN,NTOT,NTOR,PTB(0:1000),PNN(0:1000),NTOTW(0:1000)
      INTEGER CCN(5000000),CLUST(5000000),CN(5000000),SPN(0:10000)
      INTEGER MAXI,MAXCN
      REAL*8 STC,MN,TT,E,STR
      REAL*8 X1,Y1,Z1,X2,Y2,Z2
      INTEGER J,TTC,ntasks,my,PTR(0:1000),M1,M2,YN,IS,SNN


      OPEN(UNIT=112,FILE='frag',STATUS='UNKNOWN')
      OPEN(UNIT=511,FILE='maxi',STATUS='UNKNOWN')
      OPEN(UNIT=512,FILE='maxi1',STATUS='UNKNOWN')
      OPEN(UNIT=513,FILE='maxi2',STATUS='UNKNOWN')
      OPEN(UNIT=514,FILE='maxi4',STATUS='UNKNOWN')
      OPEN(UNIT=515,FILE='maxi8',STATUS='UNKNOWN')
      OPEN(UNIT=516,FILE='maxi16',STATUS='UNKNOWN')
      OPEN(UNIT=517,FILE='maxi32',STATUS='UNKNOWN')
      OPEN(UNIT=518,FILE='maxi64',STATUS='UNKNOWN')
      OPEN(UNIT=311,FILE='crkd',STATUS='UNKNOWN')

 3    FORMAT(2F11.6,' ',2I6,' ',F11.6)
 4    FORMAT(3F15.6)

      WRITE(*,*) 'ntasks,YN?'
      READ(*,*) ntasks,YN

      DO my=0,ntasks-1
      READ(311,*) NTOTW(my),PNN(my),PTR(my),PTB(my)
      END DO


      IS=0
      DO my=0,ntasks-1
      SPN(my)=IS   
      DO I=1,PNN(my)
      IS=IS+1
      CLUST(IS)=IS
      CN(IS)=0
      CCN(IS)=0
      END DO
      END DO

      SNN=0
      DO my=0,ntasks-1
      WRITE(*,*) 'myid1=',my
      OPEN(UNIT=117,FILE='FS'//na(my),STATUS='UNKNOWN')
      DO I=1,NTOTW(my)
      READ(117,*) N1,N2,X1,Y1,Z1,X2,Y2,Z2,E
      N1=N1+SNN
      N2=N2+SNN
         IF (E.NE.0.0) THEN
!         IF (E.EQ.0.0.AND.STR.GT.0.002) THEN
         M1=CLUST(N1)
         M2=CLUST(N2)
         MN=MIN(CLUST(N1),CLUST(N2))
            DO J=1+SNN,PNN(my)+SNN
            IF (CLUST(J).EQ.M1.OR.CLUST(J).EQ.M2) CLUST(J)=MN
            END DO
         CLUST(N1)=MN
         CLUST(N2)=MN
         END IF
      END DO
      CLOSE(117)
      SNN=SNN+PNN(my)
      END DO

      DO my=0,ntasks-1
      WRITE(*,*) 'myid2=',my
      OPEN(UNIT=117,FILE='FSR'//na(my),STATUS='UNKNOWN')
      DO I=1,PTR(my)
      READ(117,*) N1,N2,X1,Y1,Z1,X2,Y2,Z2,E
      N1=N1+SPN(my+1)
      N2=N2+SPN(my)
         IF (E.EQ.0.0) THEN
         M1=CLUST(N1)
         M2=CLUST(N2)
         MN=MIN(CLUST(N1),CLUST(N2))
            DO J=1,IS
            IF (CLUST(J).EQ.M1.OR.CLUST(J).EQ.M2) CLUST(J)=MN
            END DO
         CLUST(N1)=MN
         CLUST(N2)=MN
         END IF
      END DO
      CLOSE(117)
      END DO

      DO my=0,ntasks-1
      WRITE(*,*) 'myid3=',my
      OPEN(UNIT=117,FILE='FSB'//na(my),STATUS='UNKNOWN')
      DO I=1,PTB(my)
      READ(117,*) N1,N2,X1,Y1,Z1,X2,Y2,Z2,E
      N1=N1+SPN(my-ntasks/YN)
      N2=N2+SPN(my)
         IF (E.EQ.0.0) THEN
         M1=CLUST(N1)
         M2=CLUST(N2)
         MN=MIN(CLUST(N1),CLUST(N2))
            DO J=1,IS
            IF (CLUST(J).EQ.M1.OR.CLUST(J).EQ.M2) CLUST(J)=MN
            END DO
         CLUST(N1)=MN
         CLUST(N2)=MN
         END IF
      END DO
      CLOSE(117)
      END DO

      DO I=1,IS
      CN(CLUST(I))=CN(CLUST(I))+1
      END DO

      MAXCN=0
      DO I=1,IS
      CCN(CN(I))=CCN(CN(I))+1
       IF (CN(I).GT.MAXCN) THEN
       MAXCN=CN(I)
       MAXI=I
       ENDIF
      END DO

      SNN=0
      DO my=0,ntasks-1
      WRITE(*,*) 'myid4=',my
      OPEN(UNIT=117,FILE='NODFIL2'//na(my),STATUS='UNKNOWN')
       DO I=1,PNN(my)
       READ(117,*) N1,X1,Y1,Z1,E
       N1=N1+SNN
       IF (CLUST(N1).EQ.MAXI) WRITE(511,4) X1,Y1,Z1

       IF (CN(CLUST(N1)).GT.0.AND.CN(CLUST(N1)).LE.1) WRITE(512,4) X1,Y1,Z1
       IF (CN(CLUST(N1)).GT.1.AND.CN(CLUST(N1)).LE.2) WRITE(513,4) X1,Y1,Z1
       IF (CN(CLUST(N1)).GT.2.AND.CN(CLUST(N1)).LE.4) WRITE(514,4) X1,Y1,Z1
       IF (CN(CLUST(N1)).GT.4.AND.CN(CLUST(N1)).LE.8) WRITE(515,4) X1,Y1,Z1
       IF (CN(CLUST(N1)).GT.8.AND.CN(CLUST(N1)).LE.16) WRITE(516,4) X1,Y1,Z1
       IF (CN(CLUST(N1)).GT.16.AND.CN(CLUST(N1)).LE.32) WRITE(517,4) X1,Y1,Z1
       IF (CN(CLUST(N1)).GT.32.AND.CN(CLUST(N1)).LT.CN(MAXI)) WRITE(518,4) X1,Y1,Z1

       ENDDO
      CLOSE(117)
      SNN=SNN+PNN(my)
      ENDDO


      DO I=1,IS
      IF (CCN(I).NE.0) WRITE(112,*) I,CCN(I) 
      END DO


      TTC=0
      STC=0.0
      DO I=1,IS
      TTC=TTC+CCN(I)
      STC=STC+I*CCN(I)
      END DO
      STC=STC/(1.0*TTC)
      WRITE(*,*) 'STC=',STC

      CLOSE(112)
      CLOSE(511)
      CLOSE(311)

      STOP
      END





