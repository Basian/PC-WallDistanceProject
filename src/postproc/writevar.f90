SUBROUTINE writevar(NI,NJ,VAR_C_C,ALGORITHM)    !WRITES SOLUTION TO PLOT3D FILE
USE, INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE

INTEGER :: N,NI,NJ,NK,I,J,K,L,NBLK,ALGORITHM
INTEGER :: JUNKA,JUNKB,JUNKC,JUNKD,GRD
REAL*8 :: MACH, ALPHA, REYN, TIME
REAL*8, DIMENSION(NI,NJ) :: X,Y,Z
REAL*8, DIMENSION(NI,NJ,5) :: MAT_V
REAL*8, DIMENSION(NI-1,NJ-1,5) :: MAT_C
REAL*8, DIMENSION(NI-1,NJ-1,2) :: CC
REAL*8, DIMENSION(NI-1,NJ-1,5,2) :: PT


REAL(C_DOUBLE), INTENT(OUT) :: VAR_C_C((NI-1)*(NJ-1))



!----------------------------------------------
! GET CELL CENTERS
!----------------------------------------------
OPEN(UNIT=7, FILE='grid.x', FORM='UNFORMATTED')
READ(7) NBLK
READ(7) NI,NJ,NK


!--------------------------------
!	READ GRID
!--------------------------------
DO GRD=1,NBLK
  READ(7) ((X(I,J),I=1,NI),J=1,NJ), &
   &      ((Y(I,J),I=1,NI),J=1,NJ), &
   &      ((Z(I,J),I=1,NI),J=1,NJ)
END DO


!--------------------------------
!	COMPUTE CELL CENTERS
!--------------------------------
L = 1

DO I=1,NI-1
   DO J=1,NJ-1
      
      ! STORE THE FOUR POINTS ASSOCIATED WITH EACH CELL - COUNTER-CW, (I,J,P,X/Y), 1=X,2=Y
      ! FIFTH POINT CLOSES THE VOLUME
      PT(I,J,1,1) = X(I,J)
      PT(I,J,2,1) = X(I+1,J)
      PT(I,J,3,1) = X(I+1,J+1)
      PT(I,J,4,1) = X(I,J+1)
      PT(I,J,5,1) = X(I,J)

      PT(I,J,1,2) = Y(I,J)
      PT(I,J,2,2) = Y(I+1,J)
      PT(I,J,3,2) = Y(I+1,J+1)
      PT(I,J,4,2) = Y(I,J+1)
      PT(I,J,5,2) = Y(I,J)

      CC(I,J,1) = (PT(I,J,1,1)+PT(I,J,2,1)+PT(I,J,3,1)+PT(I,J,4,1))/4.0  !X
      CC(I,J,2) = (PT(I,J,1,2)+PT(I,J,2,2)+PT(I,J,3,2)+PT(I,J,4,2))/4.0  !Y

   END DO
END DO

CALL array_to_mat(NI,NJ,NK,VAR_C_C,MAT_C)
CALL cell_to_vertex(NI,NJ,X,Y,CC,MAT_C,MAT_V)

!----------------------------------------------
! WRITE VARIABLE
!----------------------------------------------

MACH=1234.0
ALPHA=1234.0
REYN=1234.0
TIME=1234.0

IF (ALGORITHM == 0) THEN
	OPEN(UNIT=8, FILE='varBF.x', FORM='UNFORMATTED')
ELSE IF (ALGORITHM == 1) THEN
	OPEN(UNIT=8, FILE='varAB1.x', FORM='UNFORMATTED')
ELSE
	OPEN(UNIT=8, FILE='varAB2.x', FORM='UNFORMATTED')
END IF


WRITE(8) NBLK
WRITE(8) NI,NJ,NK
WRITE(8) MACH, ALPHA, REYN, TIME
WRITE(8) ((( MAT_V(I,J,N), I=1,NI), J=1,NJ), N=1,5)


CLOSE(7)
CLOSE(8)

END SUBROUTINE writevar




