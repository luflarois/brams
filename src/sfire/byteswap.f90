!***********************************************************************************************************************************
!
!                                                        B Y T E S W A P
!
!  Module:       BYTESWAP
!
!  Programmer:   Dr. David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
!  Date:         July 23, 2005
!
!  Language:     Fortran-90
!
!  Description:  This module includes several subroutines to reverse the byte order of INTEGER and REAL variables.
!                This is useful when reading binary data from a file intended for use on a compute whose byte
!                order is opposite that of the computer on which the Fortran program is to be run.
!
!                   SWAP_I2    Swap bytes of a 2-byte INTEGER
!                   SWAP_I4    Swap bytes of a 4-byte INTEGER
!                   SWAP_F4    Swap bytes of a 4-byte REAL
!                   SWAP_F8    Swap bytes of an 8-byte REAL
!
!***********************************************************************************************************************************

      MODULE byteswap


      CONTAINS


!***********************************************************************************************************************************
!  SWAP_I2
!
!  Swap bytes for a two-byte INTEGER.
!  After calling this subroutine, the input integer will be replaced by the output integer.
!***********************************************************************************************************************************

      SUBROUTINE SWAP_I2 (BYTE2)

      IMPLICIT NONE

      INTEGER(KIND=2), INTENT(IN OUT) :: BYTE2

      INTEGER(KIND=1), DIMENSION(2) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (BYTE2, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 2
         BYTE_ARR(I) = BYTE_ARR_TMP(3-I)
      END DO

      BYTE2 = TRANSFER (BYTE_ARR, BYTE2)

      RETURN

      END SUBROUTINE SWAP_I2





!***********************************************************************************************************************************
!  SWAP_I4
!
!  Swap bytes for a four-byte INTEGER.
!  After calling this subroutine, the input integer will be replaced by the output integer.
!***********************************************************************************************************************************

      SUBROUTINE SWAP_I4 (BYTE4)

      IMPLICIT NONE

      INTEGER(KIND=4), INTENT(IN OUT) :: BYTE4

      INTEGER(KIND=1), DIMENSION(4) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (BYTE4, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 4
         BYTE_ARR(I) = BYTE_ARR_TMP(5-I)
      END DO

      BYTE4 = TRANSFER (BYTE_ARR, BYTE4)

      RETURN

      END SUBROUTINE SWAP_I4





!***********************************************************************************************************************************
!  SWAP_F4
!
!  Swap bytes for a four-byte REAL.
!  After calling this subroutine, the input number will be replaced by the output number.
!***********************************************************************************************************************************

      SUBROUTINE SWAP_F4 (FLOAT4)

      IMPLICIT NONE

      REAL(KIND=4), INTENT(IN OUT) :: FLOAT4

      INTEGER(KIND=1), DIMENSION(4) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (FLOAT4, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 4
         BYTE_ARR(I) = BYTE_ARR_TMP(5-I)
      END DO

      FLOAT4 = TRANSFER (BYTE_ARR, FLOAT4)

      RETURN

      END SUBROUTINE SWAP_F4





!***********************************************************************************************************************************
!  SWAP_F8
!
!  Swap bytes for an eight-byte REAL.
!  After calling this subroutine, the input number will be replaced by the output number.
!***********************************************************************************************************************************

      SUBROUTINE SWAP_F8 (FLOAT8)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN OUT) :: FLOAT8

      INTEGER(KIND=1), DIMENSION(8) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (FLOAT8, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 8
         BYTE_ARR(I) = BYTE_ARR_TMP(9-I)
      END DO

      FLOAT8 = TRANSFER (BYTE_ARR, FLOAT8)

      RETURN

      END SUBROUTINE SWAP_F8


      END MODULE byteswap

