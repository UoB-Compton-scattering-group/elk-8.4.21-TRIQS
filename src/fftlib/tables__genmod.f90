        !COMPILER-GENERATED INTERFACE MODULE: Fri May 20 16:46:50 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TABLES__genmod
          INTERFACE 
            SUBROUTINE TABLES(IDO,IP,WA)
              INTEGER(KIND=4) :: IP
              INTEGER(KIND=4) :: IDO
              REAL(KIND=8) :: WA(IDO,IP-1,2)
            END SUBROUTINE TABLES
          END INTERFACE 
        END MODULE TABLES__genmod
