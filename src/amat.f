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

      SUBROUTINE AMAT(E,W,G,X1,Y1,Z1,X2,Y2,Z2,L,DUT,I,RY,CT,NTOT)
 
      INCLUDE 'param.dat'
      REAL*8 RK(12,12),B(12),CT(NTOT*12)
      REAL*8 T(12,12),TT(12,12),DUT(12)
      REAL*8 F(12),X1,Y1,Z1,X2,Y2,Z2
      REAL*8 E,W,L,G,SPRK,SRK,SF
      INTEGER I,J,S,RY,myid,NTOT
c      COMMON /KM/ RK
c      COMMON /BT/ T
c      COMMON /BTT/ TT
c      COMMON /CT/ CT
      PARAMETER(ZERO = 0.0E+0)

 
c      OPEN(UNIT=20,FILE='TEST',STATUS='UNKNOWN')
 
c       write(20,*) (X1+X2)/2.0, L
c       IF (L.LT.0.01) L=0.01
c       write(*,*) 'E=',E
c      L=1.0
 
      CALL TTMAT(X1,Y1,Z1,X2,Y2,Z2,RY,T,.TRUE.)
      CALL KMAT(E,G,W,L,RY,RK)


        B(12-11)= 
     1              T(1,1)*DUT(6*1-5) + T(1,2)*DUT(6*1-4)
     1            + T(1,3)*DUT(6*1-3)
                                                            
        B(12-10)= 
     1              T(2,1)*DUT(6*1-5) + T(2,2)*DUT(6*1-4)
     1            + T(2,3)*DUT(6*1-3)
                                                     
        B(12-9)= 
     1              T(3,1)*DUT(6*1-5) + T(3,2)*DUT(6*1-4)
     1            + T(3,3)*DUT(6*1-3)
                                                      
        B(12-8)= 
     1               T(4,4)*DUT(6*1-2)
     1            + T(4,5)*DUT(6*1-1) + T(4,6)*DUT(6*1-0)
                                                             
        B(12-7)= 
     1              T(5,4)*DUT(6*1-2)
     1            + T(5,5)*DUT(6*1-1) + T(5,6)*DUT(6*1-0)
                                                           
        B(12-6)= 
     1              T(6,4)*DUT(6*1-2)
     1            + T(6,5)*DUT(6*1-1) + T(6,6)*DUT(6*1-0)
                                                   
        B(12-5)= 
     1              T(7,7)*DUT(6*2-5) + T(7,8)*DUT(6*2-4)
     1            + T(7,9)*DUT(6*2-3)
                                                                    
        B(12-4)= 
     1              T(8,7)*DUT(6*2-5) + T(8,8)*DUT(6*2-4)
     1            + T(8,9)*DUT(6*2-3)
                                                                        
        B(12-3)= 
     1              T(9,7)*DUT(6*2-5) + T(9,8)*DUT(6*2-4)
     1            + T(9,9)*DUT(6*2-3)
                                                             
        B(12-2)=  
     1              T(10,10)*DUT(6*2-2)
     1            + T(10,11)*DUT(6*2-1) + T(10,12)*DUT(6*2-0)
                                          
        B(12-1)= 
     1              T(11,10)*DUT(6*2-2)
     1            + T(11,11)*DUT(6*2-1) + T(11,12)*DUT(6*2-0)
                                                                              
        B(12-0)= 
     1              T(12,10)*DUT(6*2-2)
     1            + T(12,11)*DUT(6*2-1) + T(12,12)*DUT(6*2-0)
c----------------------------------------------------------------
        CT(12*I-11)= CT(12*I-11)   
     1            + RK(1,1)*B(12-11)+ RK(1,7)*B(12-5)
                                                            
        CT(12*I-10)= CT(12*I-10)
     1            + RK(2,2)*B(12-10)
     1            + RK(2,6)*B(12-6)
     1            + RK(2,8)*B(12-4)
     1            + RK(2,12)*B(12-0)
                                                     
        CT(12*I-9)=  CT(12*I-9)
     1            + RK(3,3)*B(12-9)
     1            + RK(3,5)*B(12-7)
     1            + RK(3,9)*B(12-3)
     1            + RK(3,11)*B(12-1)
                                                      
        CT(12*I-8)=  CT(12*I-8)
     1            + RK(4,4)*B(12-8)
     1            + RK(4,10)*B(12-2)

                                                             
        CT(12*I-7)= CT(12*I-7)
     1            + RK(5,3)*B(12-9)
     1            + RK(5,5)*B(12-7)
     1            + RK(5,9)*B(12-3)
     1            + RK(5,11)*B(12-1)
                                                           
        CT(12*I-6)=  CT(12*I-6)
     1            + RK(6,2)*B(12-10)
     1            + RK(6,6)*B(12-6)
     1            + RK(6,8)*B(12-4)
     1            + RK(6,12)*B(12-0)
                                                   
        CT(12*I-5)=  CT(12*I-5)
     1            + RK(7,1)*B(12-11)
     1            + RK(7,7)*B(12-5)
                                                                    
        CT(12*I-4)= CT(12*I-4)
     1            + RK(8,2)*B(12-10)
     1            + RK(8,6)*B(12-6)
     1            + RK(8,8)*B(12-4)
     1            + RK(8,12)*B(12-0)
                                                                        
        CT(12*I-3)=  CT(12*I-3)
     1            + RK(9,3)*B(12-9)
     1            + RK(9,5)*B(12-7)
     1            + RK(9,9)*B(12-3)
     1            + RK(9,11)*B(12-1)
                                                             
        CT(12*I-2)=  CT(12*I-2)
     1            + RK(10,4)*B(12-8)
     1            + RK(10,10)*B(12-2)

                                          
        CT(12*I-1)= CT(12*I-1)
     1            + RK(11,3)*B(12-9)
     1            + RK(11,5)*B(12-7)
     1            + RK(11,9)*B(12-3)
     1            + RK(11,11)*B(12-1)
                                                                              
        CT(12*I-0)= CT(12*I-0)
     1            + RK(12,2)*B(12-10)
     1            + RK(12,6)*B(12-6)
     1            + RK(12,8)*B(12-4)
     1            + RK(12,12)*B(12-0)
c-----------------------------------------------------------
       
        RETURN
        END
