SUBROUTINE array_to_mat (NI,NJ,NK,VAR_C,MAT_C)    !WRITES SOLUTION TO PLOT3D FILE
IMPLICIT NONE

!-------------------------------------------
!-------------------------------------------
!  THIS ROUTINE CONVERTS A CELL CENTERED
!  ARRAY TO A CELL CENTERED MATRIX
!-------------------------------------------
!-------------------------------------------


INTEGER :: NI,NJ,NK,I,J,L,NBLK
REAL :: MACH, ALPHA, REYN, TIME
REAL, DIMENSION(:) :: VAR_C
REAL, DIMENSION(NI-1,NJ-1,1) :: MAT_C


!-----------------------------------
!	ASSEMBLE MATRIX OUTPUT
!-----------------------------------

L = 1
DO I=1,NI-1
   DO J=1,NJ-1
     
      MAT_C(I,J) = VAR_C(L)
      L = L + 1
      
   END DO
END DO


END SUBROUTINE array_to_mat






























