        !COMPILER-GENERATED INTERFACE MODULE: Fri May 20 16:46:50 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CMF2KB__genmod
          INTERFACE 
            SUBROUTINE CMF2KB(LOT,IDO,L1,NA,CC,IM1,IN1,CH,IM2,IN2,WA)
              INTEGER(KIND=4) :: IN2
              INTEGER(KIND=4) :: IN1
              INTEGER(KIND=4) :: L1
              INTEGER(KIND=4) :: IDO
              INTEGER(KIND=4) :: LOT
              INTEGER(KIND=4) :: NA
              REAL(KIND=8) :: CC(2,IN1,L1,IDO,2)
              INTEGER(KIND=4) :: IM1
              REAL(KIND=8) :: CH(2,IN2,L1,2,IDO)
              INTEGER(KIND=4) :: IM2
              REAL(KIND=8) :: WA(IDO,1,2)
            END SUBROUTINE CMF2KB
          END INTERFACE 
        END MODULE CMF2KB__genmod
