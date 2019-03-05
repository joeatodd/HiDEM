! *************************************************************************
! *  HiDEM, A Discrete Element Model for Fracture Simulation
! *  Copyright (C) 24th May 2018 - Jan Åström
! *
! *  This program is free software: you can redistribute it and/or modify
! *  it under the terms of the GNU General Public License as published by
! *  the Free Software Foundation, either version 3 of the License, or
! *  (at your option) any later version.
! *
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
! *************************************************************************
MODULE EFFL

CONTAINS

  SUBROUTINE EFFLOAD(SI,NTOT,NN,T,M,JS,UT,UTM,R,EN,RY, &
       FXF,FXC,VDP,DPE,EFS,NANS,NRXF,MFIL,CT,LNN)

        USE INOUT
        USE TypeDefs

	IMPLICIT NONE
        REAL(KIND=dp) :: MFIL(NN),VDP(NN)
	REAL(KIND=dp) :: LNN,DPE
	REAL(KIND=dp) :: CT(NTOT*12),EN(NN*6),R(NN*6)
	REAL(KIND=dp) :: T1,T2
	REAL(KIND=dp) :: T,M,JS
        INTEGER :: NTOT,NN
        TYPE(SimInfo_t) :: SI
 !-----------------------------------------------
        REAL(KIND=dp) :: DMP,DMP2,DT,S
	REAL(KIND=dp) :: DX1,DY1,DZ1,DX2,DY2,DZ2
	REAL(KIND=dp) :: X1,Y1,Z1,X2,Y2,Z2,TT(12,12),DUT(12)
        REAL(KIND=dp),ALLOCATABLE :: EFS(:),A(:),C(:),F(:),D(:)

	INTEGER N,NL,NB,N1,N2,X,XL,XR
	INTEGER I,J,RY,FXC
        INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),comm
        INTEGER ierr
        INTEGER, ALLOCATABLE :: FXF(:,:),NANS(:,:)
        TYPE(UT_t) :: UT, UTM
        TYPE(NRXF_t) :: NRXF
        LOGICAL :: FirstTime=.TRUE.

        SAVE :: FirstTime, A,C,F,D

        CALL CPU_Time(T1)

        IF(FirstTime) THEN
          FirstTime = .FALSE.
          ALLOCATE(A(NN*6),C(NN*6),F(NN*6),D(NN*6))
        END IF

        DMP = SI%DMP
        DMP2 = SI%DMP2
        DT = SI%DT
        S = SI%S

        A = 0.0
        C = 0.0
        D = 0.0
        F = 0.0
        DUT = 0.0

 	DO X=1,NTOT
        IF (EFS(X).NE.0.0) THEN

	N1=NANS(1,X)
	N2=NANS(2,X)
	X1=NRXF%A(1,N1)
	Y1=NRXF%A(2,N1)
	Z1=NRXF%A(3,N1)
	X2=NRXF%A(1,N2)
	Y2=NRXF%A(2,N2)
	Z2=NRXF%A(3,N2)


	DX1=UT%A(6*N1-5)
	DY1=UT%A(6*N1-4)
	DZ1=UT%A(6*N1-3)
	DX2=UT%A(6*N2-5)
	DY2=UT%A(6*N2-4)
	DZ2=UT%A(6*N2-3)

        DUT(6*1-5)=UT%A(6*N1-5)-UTM%A(6*N1-5)
        DUT(6*1-4)=UT%A(6*N1-4)-UTM%A(6*N1-4)
        DUT(6*1-3)=UT%A(6*N1-3)-UTM%A(6*N1-3)
        DUT(6*1-2)=UT%A(6*N1-2)-UTM%A(6*N1-2)
        DUT(6*1-1)=UT%A(6*N1-1)-UTM%A(6*N1-1)
        DUT(6*1-0)=UT%A(6*N1-0)-UTM%A(6*N1-0)
                                                        
        DUT(6*2-5)=UT%A(6*N2-5)-UTM%A(6*N2-5)
        DUT(6*2-4)=UT%A(6*N2-4)-UTM%A(6*N2-4)
        DUT(6*2-3)=UT%A(6*N2-3)-UTM%A(6*N2-3)
        DUT(6*2-2)=UT%A(6*N2-2)-UTM%A(6*N2-2)
        DUT(6*2-1)=UT%A(6*N2-1)-UTM%A(6*N2-1)
        DUT(6*2-0)=UT%A(6*N2-0)-UTM%A(6*N2-0)

        CALL AMAT(EFS(X),S,EFS(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
             X2+DX2,Y2+DY2,Z2+DZ2,LNN,DUT,X,RY,CT,NTOT)

	ELSE

        CT(12*X-11)=0.0
        CT(12*X-10)=0.0
        CT(12*X-9)=0.0
        CT(12*X-8)=0.0
        CT(12*X-7)=0.0
        CT(12*X-6)=0.0
        CT(12*X-5)=0.0
        CT(12*X-4)=0.0
        CT(12*X-3)=0.0
        CT(12*X-2)=0.0
        CT(12*X-1)=0.0
        CT(12*X-0)=0.0
        ENDIF
      END DO

        DO X=1,NTOT
        N1=NANS(1,X)
        N2=NANS(2,X)
        X1=NRXF%A(1,N1)
        Y1=NRXF%A(2,N1)
        Z1=NRXF%A(3,N1)
        X2=NRXF%A(1,N2)
        Y2=NRXF%A(2,N2)
        Z2=NRXF%A(3,N2)
        DX1=UT%A(6*N1-5)
        DY1=UT%A(6*N1-4)
        DZ1=UT%A(6*N1-3)
        DX2=UT%A(6*N2-5)
        DY2=UT%A(6*N2-4)
        DZ2=UT%A(6*N2-3)
        CALL TTMAT(X1+DX1,Y1+DY1,Z1+DZ1,X2+DX2,Y2+DY2,Z2+DZ2,RY,TT)

        IF(N1 <= NN) THEN
        A(6*N1-5)= A(6*N1-5) &
                  +  TT(1,1)*CT(12*X-11) + TT(1,2)*CT(12*X-10) &
                  + TT(1,3)*CT(12*X-9)
        A(6*N1-4)= A(6*N1-4) &
                  +  TT(2,1)*CT(12*X-11) + TT(2,2)*CT(12*X-10) &
                  + TT(2,3)*CT(12*X-9)
        A(6*N1-3)= A(6*N1-3) &
                  +  TT(3,1)*CT(12*X-11) + TT(3,2)*CT(12*X-10) &
                  + TT(3,3)*CT(12*X-9)
        A(6*N1-2)= A(6*N1-2) &
                  + TT(4,4)*CT(12*X-8) &
                  + TT(4,5)*CT(12*X-7) + TT(4,6)*CT(12*X-6)
        A(6*N1-1)= A(6*N1-1) &
                  + TT(5,4)*CT(12*X-8) &
                  + TT(5,5)*CT(12*X-7) + TT(5,6)*CT(12*X-6)
        A(6*N1-0)= A(6*N1-0) &
       	          + TT(6,4)*CT(12*X-8) &
                  + TT(6,5)*CT(12*X-7) + TT(6,6)*CT(12*X-6)
        END IF
        IF(N2 <= NN) THEN
        A(6*N2-5)= A(6*N2-5) &
                  + TT(7,7)*CT(12*X-5) + TT(7,8)*CT(12*X-4) &
                  + TT(7,9)*CT(12*X-3)
        A(6*N2-4)= A(6*N2-4) &
                  + TT(8,7)*CT(12*X-5) + TT(8,8)*CT(12*X-4) &
                  + TT(8,9)*CT(12*X-3)
        A(6*N2-3)= A(6*N2-3) &
                  + TT(9,7)*CT(12*X-5) + TT(9,8)*CT(12*X-4) &
                  + TT(9,9)*CT(12*X-3)
        A(6*N2-2)= A(6*N2-2) &
                  + TT(10,10)*CT(12*X-2) + TT(10,11)*CT(12*X-1) &
      	          + TT(10,12)*CT(12*X-0)
        A(6*N2-1)= A(6*N2-1) &
      	          + TT(11,10)*CT(12*X-2) + TT(11,11)*CT(12*X-1) &
      	          + TT(11,12)*CT(12*X-0)
        A(6*N2-0)= A(6*N2-0) &
      	          + TT(12,10)*CT(12*X-2) + TT(12,11)*CT(12*X-1) &
      	          + TT(12,12)*CT(12*X-0)
        END IF
	ENDDO


	DO X=1,NN
	C(6*X-5)= (MFIL(X)/DT**2)*UTM%M(6*X-5)-(2*MFIL(X)/DT**2)*UT%M(6*X-5)
	C(6*X-4)= (MFIL(X)/DT**2)*UTM%M(6*X-4)-(2*MFIL(X)/DT**2)*UT%M(6*X-4)
	C(6*X-3)= (MFIL(X)/DT**2)*UTM%M(6*X-3)-(2*MFIL(X)/DT**2)*UT%M(6*X-3)
	C(6*X-2)= ((MFIL(X)*JS/M)/DT**2)* &
      	UTM%M(6*X-2)-(2*(MFIL(X)*JS/M)/DT**2)*UT%M(6*X-2)
	C(6*X-1)= ((MFIL(X)*JS/M)/DT**2)* &
      	UTM%M(6*X-1)-(2*(MFIL(X)*JS/M)/DT**2)*UT%M(6*X-1)
	C(6*X-0)= ((MFIL(X)*JS/M)/DT**2)* &
      	UTM%M(6*X-0)-(2*(MFIL(X)*JS/M)/DT**2)*UT%M(6*X-0)
	ENDDO


	DO X=1,NTOT
	IF (EFS(X).NE.0.0) THEN
	N1=NANS(1,X)
	N2=NANS(2,X)
        IF(N1 <= NN) THEN
	D(6*N1-5)=D(6*N1-5)+(DMP/DT)*((UT%A(6*N1-5)-UT%A(6*N2-5)) &
      	-(UTM%A(6*N1-5)-UTM%A(6*N2-5)))
	D(6*N1-4)=D(6*N1-4)+(DMP/DT)*((UT%A(6*N1-4)-UT%A(6*N2-4)) &
      	-(UTM%A(6*N1-4)-UTM%A(6*N2-4)))
	D(6*N1-3)=D(6*N1-3)+(DMP/DT)*((UT%A(6*N1-3)-UT%A(6*N2-3)) &
      	-(UTM%A(6*N1-3)-UTM%A(6*N2-3)))
	D(6*N1-2)=D(6*N1-2)+(DMP2/DT)*((UT%A(6*N1-2)-UT%A(6*N2-2)) &
      	-(UTM%A(6*N1-2)-UTM%A(6*N2-2)))
	D(6*N1-1)=D(6*N1-1)+(DMP2/DT)*((UT%A(6*N1-1)-UT%A(6*N2-1)) &
      	-(UTM%A(6*N1-1)-UTM%A(6*N2-1)))
	D(6*N1-0)=D(6*N1-0)+(DMP2/DT)*((UT%A(6*N1-0)-UT%A(6*N2-0)) &
      	-(UTM%A(6*N1-0)-UTM%A(6*N2-0)))
        END IF
        IF(N2 <= NN) THEN
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%A(6*N2-5)-UT%A(6*N1-5)) &
      	-(UTM%A(6*N2-5)-UTM%A(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%A(6*N2-4)-UT%A(6*N1-4)) &
      	-(UTM%A(6*N2-4)-UTM%A(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%A(6*N2-3)-UT%A(6*N1-3)) &
      	-(UTM%A(6*N2-3)-UTM%A(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%A(6*N2-2)-UT%A(6*N1-2)) &
      	-(UTM%A(6*N2-2)-UTM%A(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%A(6*N2-1)-UT%A(6*N1-1)) &
      	-(UTM%A(6*N2-1)-UTM%A(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%A(6*N2-0)-UT%A(6*N1-0)) &
      	-(UTM%A(6*N2-0)-UTM%A(6*N1-0)))
        END IF

        !TODO - investigate if doing energy calcs every N timesteps 
        ! makes things more efficient - below also
        DPE=DPE+(DMP/DT)*((UT%A(6*N2-5)-UT%A(6*N1-5)) &
             -(UTM%A(6*N2-5)-UTM%A(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%A(6*N2-4)-UT%A(6*N1-4)) &
             -(UTM%A(6*N2-4)-UTM%A(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%A(6*N2-3)-UT%A(6*N1-3)) &
             -(UTM%A(6*N2-3)-UTM%A(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT%A(6*N2-2)-UT%A(6*N1-2)) &
             -(UTM%A(6*N2-2)-UTM%A(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%A(6*N2-1)-UT%A(6*N1-1)) &
             -(UTM%A(6*N2-1)-UTM%A(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%A(6*N2-0)-UT%A(6*N1-0)) &
             -(UTM%A(6*N2-0)-UTM%A(6*N1-0)))**2

	ENDIF
	ENDDO

	DO X=1,FXC
	N1=FXF(1,X)
	N2=FXF(2,X)
        IF(N1 <= NN) THEN
	D(6*N1-5)=D(6*N1-5)+(DMP/DT)*((UT%A(6*N1-5)-UT%A(6*N2-5)) &
     	-(UTM%A(6*N1-5)-UTM%A(6*N2-5)))
	D(6*N1-4)=D(6*N1-4)+(DMP/DT)*((UT%A(6*N1-4)-UT%A(6*N2-4)) &
     	-(UTM%A(6*N1-4)-UTM%A(6*N2-4)))
	D(6*N1-3)=D(6*N1-3)+(DMP/DT)*((UT%A(6*N1-3)-UT%A(6*N2-3)) &
     	-(UTM%A(6*N1-3)-UTM%A(6*N2-3)))
	D(6*N1-2)=D(6*N1-2)+(DMP2/DT)*((UT%A(6*N1-2)-UT%A(6*N2-2)) &
     	-(UTM%A(6*N1-2)-UTM%A(6*N2-2)))
	D(6*N1-1)=D(6*N1-1)+(DMP2/DT)*((UT%A(6*N1-1)-UT%A(6*N2-1)) &
     	-(UTM%A(6*N1-1)-UTM%A(6*N2-1)))
	D(6*N1-0)=D(6*N1-0)+(DMP2/DT)*((UT%A(6*N1-0)-UT%A(6*N2-0)) &
     	-(UTM%A(6*N1-0)-UTM%A(6*N2-0)))
        END IF
        IF(N2 <= NN) THEN
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%A(6*N2-5)-UT%A(6*N1-5)) &
     	-(UTM%A(6*N2-5)-UTM%A(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%A(6*N2-4)-UT%A(6*N1-4)) &
     	-(UTM%A(6*N2-4)-UTM%A(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%A(6*N2-3)-UT%A(6*N1-3)) &
     	-(UTM%A(6*N2-3)-UTM%A(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%A(6*N2-2)-UT%A(6*N1-2)) &
     	-(UTM%A(6*N2-2)-UTM%A(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%A(6*N2-1)-UT%A(6*N1-1)) &
     	-(UTM%A(6*N2-1)-UTM%A(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%A(6*N2-0)-UT%A(6*N1-0)) &
     	-(UTM%A(6*N2-0)-UTM%A(6*N1-0)))
        END IF

        DPE=DPE+(DMP/DT)*((UT%A(6*N2-5)-UT%A(6*N1-5)) &
             -(UTM%A(6*N2-5)-UTM%A(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%A(6*N2-4)-UT%A(6*N1-4)) &
             -(UTM%A(6*N2-4)-UTM%A(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%A(6*N2-3)-UT%A(6*N1-3)) &
             -(UTM%A(6*N2-3)-UTM%A(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT%A(6*N2-2)-UT%A(6*N1-2)) &
             -(UTM%A(6*N2-2)-UTM%A(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%A(6*N2-1)-UT%A(6*N1-1)) &
             -(UTM%A(6*N2-1)-UTM%A(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%A(6*N2-0)-UT%A(6*N1-0)) &
             -(UTM%A(6*N2-0)-UTM%A(6*N1-0)))**2
	ENDDO

        DO X=1,NN
        F(6*X-5)= -(VDP(X)/(DT*2))*UTM%M(6*X-5)
        F(6*X-4)= -(VDP(X)/(DT*2))*UTM%M(6*X-4)
        F(6*X-3)= -(VDP(X)/(DT*2))*UTM%M(6*X-3)
        F(6*X-2)= -(VDP(X)/(DT*2))*UTM%M(6*X-2)
        F(6*X-1)= -(VDP(X)/(DT*2))*UTM%M(6*X-1)
        F(6*X-0)= -(VDP(X)/(DT*2))*UTM%M(6*X-0)
	END DO

        !TODO - comment equation breakdown here for future ref
	R=-A -C -D -F
	EN=EN+A*(UT%M-UTM%M)

        IF(PrintTimes) THEN
          CALL CPU_TIME(T2)
          PRINT *,myid,' Effload took: ',T2-T1,' secs'
        END IF

	RETURN
END SUBROUTINE EFFLOAD

END MODULE EFFL
