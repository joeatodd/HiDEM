MODULE TypeDefs

  INCLUDE 'param.dat'

  !Types are defined here to save passing large numbers of arguments.
  !Type members refer to the arrays belonging to partitions:
  !M = 'mine', F,B,R,L = Forward, Back, Right, Left, etc

  TYPE NAN_t
     INTEGER :: M(3,NOCON), L(3,NOMA), R(3,NOMA), F(3,NOMA), B(3,NOMA), &
          FR(3,NOMA), FL(3,NOMA), BR(3,NOMA), BL(3,NOMA)
  END TYPE NAN_t

  TYPE NTOT_t
     INTEGER :: M, L, R, F, B, FR, FL, BR, BL
  END TYPE NTOT_t

  TYPE NRXF_t
     REAL*8 :: M(3,NOMA), L(3,NOMA), R(3,NOMA), F(3,NOMA), B(3,NOMA), &
          FR(3,NOMA), FL(3,NOMA), BR(3,NOMA), BL(3,NOMA)
  END TYPE NRXF_t

  TYPE EF_t
     REAL*8 :: M(NOCON),L(NOCON),R(NOCON),F(NOCON),B(NOCON),FR(NOCON),&
          FL(NOCON),BR(NOCON),BL(NOCON)
  END TYPE EF_t

  TYPE NEI_t
     INTEGER :: L=-1,R=-1,F=-1,B=-1,FR=-1,FL=-1,BR=-1,BL=-1
  END TYPE NEI_t

  TYPE UT_t
     REAL*8 :: M(NODM), L(NODM), R(NODM), F(NODM), B(NODM), FR(NODM),&
          FL(NODM), BR(NODM), BL(NODM)
  END TYPE UT_t


END MODULE TypeDefs
