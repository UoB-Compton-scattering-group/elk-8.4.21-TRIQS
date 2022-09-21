        !COMPILER-GENERATED INTERFACE MODULE: Fri May 20 16:46:49 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CFFTND__genmod
          INTERFACE 
            SUBROUTINE CFFTND(ND,N,SGN,C)
              INTEGER(KIND=4), INTENT(IN) :: ND
              INTEGER(KIND=4), INTENT(IN) :: N(ND)
              INTEGER(KIND=4), INTENT(IN) :: SGN
              COMPLEX(KIND=8), INTENT(INOUT) :: C(*)
            END SUBROUTINE CFFTND
          END INTERFACE 
        END MODULE CFFTND__genmod
