        !COMPILER-GENERATED INTERFACE MODULE: Fri May 20 16:46:49 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CFFTMI__genmod
          INTERFACE 
            SUBROUTINE CFFTMI(N,WSAVE,LENSAV)
              INTEGER(KIND=4) :: LENSAV
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: WSAVE(LENSAV)
            END SUBROUTINE CFFTMI
          END INTERFACE 
        END MODULE CFFTMI__genmod
