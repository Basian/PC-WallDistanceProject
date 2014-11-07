SUBROUTINE array_to_mat (NI,NJ,NK,VAR_C,MAT_C)    !WRITES SOLUTION TO PLOT3D FILE
USE, INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE

!-------------------------------------------
!-------------------------------------------
!  THIS ROUTINE CONVERTS A CELL CENTERED
!  ARRAY TO A CELL CENTERED MATRIX
!-------------------------------------------
!-------------------------------------------


INTEGER :: NI,NJ,NK,I,J,K,L,NBLK
!REAL*8, DIMENSION(:), POINTER :: VAR_C
REAL(C_DOUBLE), INTENT(OUT) :: VAR_C((NI-1)*(NJ-1))
REAL*8, DIMENSION(NI-1,NJ-1,5) :: MAT_C


!-----------------------------------
!	ASSEMBLE MATRIX OUTPUT
!-----------------------------------

L = 1
DO I=1,NI-1
   DO J=1,NJ-1
     
      DO K=1,5
         MAT_C(I,J,K) = VAR_C(L)
      END DO
      L = L + 1
      
   END DO
END DO


END SUBROUTINE array_to_mat






























