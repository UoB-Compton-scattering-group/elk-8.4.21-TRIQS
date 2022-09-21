        !COMPILER-GENERATED INTERFACE MODULE: Fri May 20 16:46:50 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CMFGKB__genmod
          INTERFACE 
            SUBROUTINE CMFGKB(LOT,IDO,IP,L1,LID,NA,CC,CC1,IM1,IN1,CH,CH1&
     &,IM2,IN2,WA)
              INTEGER(KIND=4) :: IN2
              INTEGER(KIND=4) :: IN1
              INTEGER(KIND=4) :: LID
              INTEGER(KIND=4) :: L1
              INTEGER(KIND=4) :: IP
              INTEGER(KIND=4) :: IDO
              INTEGER(KIND=4) :: LOT
              INTEGER(KIND=4) :: NA
              REAL(KIND=8) :: CC(2,IN1,L1,IP,IDO)
              REAL(KIND=8) :: CC1(2,IN1,LID,IP)
              INTEGER(KIND=4) :: IM1
              REAL(KIND=8) :: CH(2,IN2,L1,IDO,IP)
              REAL(KIND=8) :: CH1(2,IN2,LID,IP)
              INTEGER(KIND=4) :: IM2
              REAL(KIND=8) :: WA(IDO,IP-1,2)
            END SUBROUTINE CMFGKB
          END INTERFACE 
        END MODULE CMFGKB__genmod
