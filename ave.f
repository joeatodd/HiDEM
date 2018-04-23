      PROGRAM STAND

      REAL*8 AFIL(2000),NFIL(2000),DX,DY,J,LJ,DR,DT
      REAL*8 STR
      INTEGER IDX,N,T,I,SDX
      CHARACTER*2 na(0:64)
        data na/'00','01','02','03','04','05','06','07','08',
     1  '09','10','11','12','13','14','15','16','17',
     1  '18','19','20','21','22','23','24','25','26',
     1  '27','28','29','30','31','32','33','34','35',
     1  '36','37','38','39','40','41','42','43','44',
     1  '45','46','47','48','49','50','51','52','53',
     1  '54','55','56','57','58','59','60','61','62','63','64'/
      CHARACTER*20 INFI

 2    FORMAT (A)
c      WRITE(*,*) 'N?'
      WRITE(*,*) 'INFI?,N'
      READ(*,2) INFI
      READ(*,*) N

      SDX=0
      OPEN(UNIT=10,FILE=INFI,STATUS='OLD')
      OPEN(UNIT=700,FILE='AR',STATUS='UNKNOWN')

      DO 5 I=1,2000
      AFIL(I)=0.0
      NFIL(I)=0.0
 5    CONTINUE


      DO 20 T=1,N
      READ(10,*) DX,DR
!      DX=(DX*30.0*30.0)/1000.0
!      DX=(DX*0.03*0.03*0.03)**0.667
!      DR=DR**0.5
      IF (DX.LT.10) THEN
      IDX=DX
      AFIL(IDX)=AFIL(IDX)+DR
      GOTO 20
      ENDIF

      
!      DX=(DX*0.03**3.0)**0.77/6.0

      DO 10 I=1,200
      J=I-50
c      J=(I+1)/7.0
c      LJ=(I-1)/7.0
      RI=EXP(J/4.0)
      ORI=EXP((J-1)/4.0)
      IF (DX.GE.ORI.AND.DX.LT.RI) THEN
      AFIL(I)=AFIL(I)+DR
      NFIL(I)=NFIL(I)+1.0
      ENDIF
 10   CONTINUE
 20   CONTINUE


      DO 23 T=1,9
!      WRITE(700,*) (T*0.03**3.0)**(1.0-T/40.0)/6.0,AFIL(T)/
!     1    (((T+0.5)*0.03**3.0)**(1.0-T/40.0)
!     1    /6.0-((T-0.5)*0.03**3.0)**(1.0-T/40.0)/6.0)
      WRITE(700,*) T,AFIL(T)
 23   CONTINUE


      SM=0
      DO 30 T=1,200
      J=T-50
      IF (NFIL(T).GE.1) THEN
!      write(*,*) NFIL(T)
!      WRITE(70,*) (EXP(J/10.)+EXP((J-1)/10.))/2.,(AFIL(T)/NFIL(T))
      WRITE(700,*) EXP((T-50.0-0.5)/4.0),
     1    AFIL(T)/(EXP((J)/4.)-EXP((J-1.0)/4.))
      ENDIF
 30   CONTINUE

      STOP
      END



















