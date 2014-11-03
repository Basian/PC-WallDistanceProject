SUBROUTINE writevar (NBLK,NI,NJ,NK,VAR_C)    !WRITES SOLUTION TO PLOT3D FILE
IMPLICIT NONE

INTEGER :: N,NI,NJ,NK,I,J,K,NBLK,JUNK
REAL*8 :: MACH, ALPHA, REYN, TIME
REAL*8, DIMENSION(NI,NJ) :: X,Y,Z
REAL*8, DIMENSION(:) :: VAR_C
REAL*8, DIMENSION(NI,NJ,1) :: MAT_V
REAL*8, DIMENSION(NI-1,NJ-1,1) :: MAT_C


!------------------------------
! GET CELL CENTERS
!------------------------------

OPEN(UNIT=7, FILE='grid.x', FORM='UNFORMATTED')
READ(7) JUNK
READ(7) JUNK,JUNK,JUNK















MACH=1234.0
ALPHA=1234.0
REYN=1234.0
TIME=1234.0

CALL array_to_mat(NI,NJ,NK,VAR_C,MAT_C)
CALL cell_to_vertex(NI,NJ,X,Y,CT,MAT_C,MAT_V)



OPEN(UNIT=7, FILE='var.x', FORM='UNFORMATTED')
WRITE(7) NBLK
WRITE(7) NI,NJ,NK
WRITE(7) MACH, ALPHA, REYN, TIME
WRITE(7) ((( MAT_V(I,J,N), I=1,NI), J=1,NJ), N=1)

END SUBROUTINE writevar





























