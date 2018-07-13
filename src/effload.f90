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

	SUBROUTINE EFFLOAD(S,NTOT,NN,T,DT,M,JS,DMP,DMP2,UT,UTM,R,EN,RY, &
     FXF,FXC,VDP,DPE,EFS,NANS,NRXF,MFIL,CT,L,PNN,YN)

        USE INOUT
        USE TypeDefs

	IMPLICIT NONE
        include 'mpif.h'
        REAL*8 MFIL(NOMA),VDP(NN)
	REAL*8 A(NODM),C(NODM),F(NODM),D(NODM)
	REAL*8 LNN
	REAL*8 DUT(NODM),CT(NODC),EN(NODM),R(NODM)
	REAL*8 DPE,S,E
	REAL*8 T,DT,M,JS,L,ALF,DMP
	REAL*8 G,X1,Y1,Z1,X2,Y2,Z2,TT(12,12)
	REAL*8 DX1,DY1,DZ1,DX2,DY2,DZ2,DMP2
        REAL*8,ALLOCATABLE :: EFS(:)
	INTEGER N,NL,NB,N1,N2,X,XL,XR,PNN(0:5000)
	INTEGER I,J,NN,RY,YN
        INTEGER dest,source,tag,stat(MPI_STATUS_SIZE),comm
        INTEGER ierr
        TYPE(NAN_t) :: NANS
        TYPE(NTOT_t) :: NTOT,FXC
        TYPE(UT_t) :: UT, UTM
        TYPE(NRXF_t) :: NRXF
        TYPE(FXF_t) :: FXF
	DO I=1,6*NN
	A(I)=0.0
	D(I)=0.0
        END DO

        !TODO - need to pass UTM before this point (for all conn & prox nodes?)
!------------------------------------------------------
      !Note - use of DUT here suggests N1 and N2 should be mutually 
      !exclusive sets? but they aren't - NOTE - may be now (my mods)
!----------------------------------------------

 	DO X=1,NTOT%M
	N1=NANS%M(1,X)
	N2=NANS%M(2,X)
	X1=NRXF%M(1,N1)
	Y1=NRXF%M(2,N1)
	Z1=NRXF%M(3,N1)
	X2=NRXF%M(1,N2)
	Y2=NRXF%M(2,N2)
	Z2=NRXF%M(3,N2)


	DX1=UT%M(6*N1-5)
	DY1=UT%M(6*N1-4)
	DZ1=UT%M(6*N1-3)
	DX2=UT%M(6*N2-5)
	DY2=UT%M(6*N2-4)
	DZ2=UT%M(6*N2-3)

        DUT(6*N1-5)=UT%M(6*N1-5)-UTM%M(6*N1-5)
        DUT(6*N1-4)=UT%M(6*N1-4)-UTM%M(6*N1-4)
        DUT(6*N1-3)=UT%M(6*N1-3)-UTM%M(6*N1-3)
        DUT(6*N1-2)=UT%M(6*N1-2)-UTM%M(6*N1-2)
        DUT(6*N1-1)=UT%M(6*N1-1)-UTM%M(6*N1-1)
        DUT(6*N1-0)=UT%M(6*N1-0)-UTM%M(6*N1-0)
                                                        
        DUT(6*N2-5)=UT%M(6*N2-5)-UTM%M(6*N2-5)
        DUT(6*N2-4)=UT%M(6*N2-4)-UTM%M(6*N2-4)
        DUT(6*N2-3)=UT%M(6*N2-3)-UTM%M(6*N2-3)
        DUT(6*N2-2)=UT%M(6*N2-2)-UTM%M(6*N2-2)
        DUT(6*N2-1)=UT%M(6*N2-1)-UTM%M(6*N2-1)
        DUT(6*N2-0)=UT%M(6*N2-0)-UTM%M(6*N2-0)

        IF (EFS(X).NE.0.0) THEN
        CALL AMAT(EFS(X),S,EFS(X)/2.0,X1+DX1,Y1+DY1,Z1+DZ1, &
             X2+DX2,Y2+DY2,Z2+DZ2,L,DUT,N1,N2,X,RY,CT)
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

!-------------------------------------------------------
! TODO - again, A was only filled in for *our* particles

        DO X=1,NTOT%M
        N1=NANS % M(1,X)
        N2=NANS % M(2,X)
        X1=NRXF%M(1,N1)
        Y1=NRXF%M(2,N1)
        Z1=NRXF%M(3,N1)
        X2=NRXF%M(1,N2)
        Y2=NRXF%M(2,N2)
        Z2=NRXF%M(3,N2)
        DX1=UT%M(6*N1-5)
        DY1=UT%M(6*N1-4)
        DZ1=UT%M(6*N1-3)
        DX2=UT%M(6*N2-5)
        DY2=UT%M(6*N2-4)
        DZ2=UT%M(6*N2-3)
        CALL TTMAT(X1+DX1,Y1+DY1,Z1+DZ1,X2+DX2,Y2+DY2,Z2+DZ2,RY,TT)
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
	ENDDO

!------------------------------------------------------------------------

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

!------------------------------------------------------------------------
! TODO - when this was %L,%R etc, D was only saved for our nodes (not shared)


	DO X=1,NTOT%M
	IF (EFS(X).NE.0.0) THEN
	N1=NANS % M(1,X)
	N2=NANS % M(2,X)
	D(6*N1-5)=D(6*N1-5)+(DMP/DT)*((UT%M(6*N1-5)-UT%M(6*N2-5)) &
      	-(UTM%M(6*N1-5)-UTM%M(6*N2-5)))
	D(6*N1-4)=D(6*N1-4)+(DMP/DT)*((UT%M(6*N1-4)-UT%M(6*N2-4)) &
      	-(UTM%M(6*N1-4)-UTM%M(6*N2-4)))
	D(6*N1-3)=D(6*N1-3)+(DMP/DT)*((UT%M(6*N1-3)-UT%M(6*N2-3)) &
      	-(UTM%M(6*N1-3)-UTM%M(6*N2-3)))
	D(6*N1-2)=D(6*N1-2)+(DMP2/DT)*((UT%M(6*N1-2)-UT%M(6*N2-2)) &
      	-(UTM%M(6*N1-2)-UTM%M(6*N2-2)))
	D(6*N1-1)=D(6*N1-1)+(DMP2/DT)*((UT%M(6*N1-1)-UT%M(6*N2-1)) &
      	-(UTM%M(6*N1-1)-UTM%M(6*N2-1)))
	D(6*N1-0)=D(6*N1-0)+(DMP2/DT)*((UT%M(6*N1-0)-UT%M(6*N2-0)) &
      	-(UTM%M(6*N1-0)-UTM%M(6*N2-0)))
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%M(6*N1-5)) &
      	-(UTM%M(6*N2-5)-UTM%M(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%M(6*N1-4)) &
      	-(UTM%M(6*N2-4)-UTM%M(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%M(6*N1-3)) &
      	-(UTM%M(6*N2-3)-UTM%M(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%M(6*N1-2)) &
      	-(UTM%M(6*N2-2)-UTM%M(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%M(6*N1-1)) &
      	-(UTM%M(6*N2-1)-UTM%M(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%M(6*N1-0)) &
      	-(UTM%M(6*N2-0)-UTM%M(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%M(6*N1-5)) &
     	-(UTM%M(6*N2-5)-UTM%M(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%M(6*N1-4)) &
     	-(UTM%M(6*N2-4)-UTM%M(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%M(6*N1-3)) &
     	-(UTM%M(6*N2-3)-UTM%M(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%M(6*N1-2)) &
     	-(UTM%M(6*N2-2)-UTM%M(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%M(6*N1-1)) &
     	-(UTM%M(6*N2-1)-UTM%M(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%M(6*N1-0)) &
     	-(UTM%M(6*N2-0)-UTM%M(6*N1-0)))**2
	ENDIF
	ENDDO

!-----------------------------------------------------------------
! TODO - when this was %L,%R etc, D was only saved for our nodes (not shared)

	DO X=1,FXC%M
	N1=FXF%M(1,X)
	N2=FXF%M(2,X)
	D(6*N1-5)=D(6*N1-5)+(DMP/DT)*((UT%M(6*N1-5)-UT%M(6*N2-5)) &
     	-(UTM%M(6*N1-5)-UTM%M(6*N2-5)))
	D(6*N1-4)=D(6*N1-4)+(DMP/DT)*((UT%M(6*N1-4)-UT%M(6*N2-4)) &
     	-(UTM%M(6*N1-4)-UTM%M(6*N2-4)))
	D(6*N1-3)=D(6*N1-3)+(DMP/DT)*((UT%M(6*N1-3)-UT%M(6*N2-3)) &
     	-(UTM%M(6*N1-3)-UTM%M(6*N2-3)))
	D(6*N1-2)=D(6*N1-2)+(DMP2/DT)*((UT%M(6*N1-2)-UT%M(6*N2-2)) &
     	-(UTM%M(6*N1-2)-UTM%M(6*N2-2)))
	D(6*N1-1)=D(6*N1-1)+(DMP2/DT)*((UT%M(6*N1-1)-UT%M(6*N2-1)) &
     	-(UTM%M(6*N1-1)-UTM%M(6*N2-1)))
	D(6*N1-0)=D(6*N1-0)+(DMP2/DT)*((UT%M(6*N1-0)-UT%M(6*N2-0)) &
     	-(UTM%M(6*N1-0)-UTM%M(6*N2-0)))
	D(6*N2-5)=D(6*N2-5)+(DMP/DT)*((UT%M(6*N2-5)-UT%M(6*N1-5)) &
     	-(UTM%M(6*N2-5)-UTM%M(6*N1-5)))
	D(6*N2-4)=D(6*N2-4)+(DMP/DT)*((UT%M(6*N2-4)-UT%M(6*N1-4)) &
     	-(UTM%M(6*N2-4)-UTM%M(6*N1-4)))
	D(6*N2-3)=D(6*N2-3)+(DMP/DT)*((UT%M(6*N2-3)-UT%M(6*N1-3)) &
     	-(UTM%M(6*N2-3)-UTM%M(6*N1-3)))
	D(6*N2-2)=D(6*N2-2)+(DMP2/DT)*((UT%M(6*N2-2)-UT%M(6*N1-2)) &
     	-(UTM%M(6*N2-2)-UTM%M(6*N1-2)))
	D(6*N2-1)=D(6*N2-1)+(DMP2/DT)*((UT%M(6*N2-1)-UT%M(6*N1-1)) &
     	-(UTM%M(6*N2-1)-UTM%M(6*N1-1)))
	D(6*N2-0)=D(6*N2-0)+(DMP2/DT)*((UT%M(6*N2-0)-UT%M(6*N1-0)) &
     	-(UTM%M(6*N2-0)-UTM%M(6*N1-0)))
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-5)-UT%M(6*N1-5)) &
     	-(UTM%M(6*N2-5)-UTM%M(6*N1-5)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-4)-UT%M(6*N1-4)) &
     	-(UTM%M(6*N2-4)-UTM%M(6*N1-4)))**2
        DPE=DPE+(DMP/DT)*((UT%M(6*N2-3)-UT%M(6*N1-3)) &
     	-(UTM%M(6*N2-3)-UTM%M(6*N1-3)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-2)-UT%M(6*N1-2)) &
     	-(UTM%M(6*N2-2)-UTM%M(6*N1-2)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-1)-UT%M(6*N1-1)) &
     	-(UTM%M(6*N2-1)-UTM%M(6*N1-1)))**2
        DPE=DPE+(DMP2/DT)*((UT%M(6*N2-0)-UT%M(6*N1-0)) &
     	-(UTM%M(6*N2-0)-UTM%M(6*N1-0)))**2
	ENDDO
!---------------------------------------------------------

        DO X=1,NN
        F(6*X-5)= -(VDP(X)/(DT*2))*UTM%M(6*X-5)
        F(6*X-4)= -(VDP(X)/(DT*2))*UTM%M(6*X-4)
        F(6*X-3)= -(VDP(X)/(DT*2))*UTM%M(6*X-3)
        F(6*X-2)= -(VDP(X)/(DT*2))*UTM%M(6*X-2)
        F(6*X-1)= -(VDP(X)/(DT*2))*UTM%M(6*X-1)
        F(6*X-0)= -(VDP(X)/(DT*2))*UTM%M(6*X-0)
	END DO

	DO X=1,6*NN
	R(X)=-A(X)-C(X)-D(X)-F(X)
	EN(X)=EN(X)+A(X)*(UT%M(X)-UTM%M(X))
	END DO

	RETURN
END SUBROUTINE EFFLOAD
