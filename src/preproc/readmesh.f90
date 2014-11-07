SUBROUTINE readmesh(NI,NJ,CC_X_OUT,CC_Y_OUT,FC_X_OUT,FC_Y_OUT)  ! READS PLOT3D MESH
USE, INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE

INTEGER :: I,J,K,L,NI,NJ,NK,NBLK,GRD
REAL*8 :: R
REAL*8, DIMENSION(:,:), ALLOCATABLE :: X,Y,Z
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: CC_MAT
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: PT
REAL*8, DIMENSION(:), POINTER :: CC_X,CC_Y,FC_X,FC_Y
TYPE(C_PTR) :: CC_X_OUT,CC_Y_OUT,FC_X_OUT,FC_Y_OUT

!--------------------------------
!	READ GRID HEADER
!--------------------------------
OPEN(UNIT=7, FILE='grid.x', FORM='UNFORMATTED')
READ(7) NBLK
READ(7) NI,NJ,NK

!--------------------------------
!	ALLOCATE MEMORY 
!--------------------------------
ALLOCATE(X(NI,NJ))
ALLOCATE(Y(NI,NJ))
ALLOCATE(Z(NI,NJ))

ALLOCATE(PT(NI-1,NJ-1,5,2))
ALLOCATE(CC_MAT(NI-1,NJ-1,2))
ALLOCATE(CC_X((NI-1)*(NJ-1)))
ALLOCATE(CC_Y((NI-1)*(NJ-1)))
ALLOCATE(FC_X((NI-1)))
ALLOCATE(FC_Y((NI-1)))

!--------------------------------
!	READ GRID
!--------------------------------
DO GRD=1,NBLK
  READ(7) ((X(I,J),I=1,NI),J=1,NJ), &
   &      ((Y(I,J),I=1,NI),J=1,NJ), &
   &      ((Z(I,J),I=1,NI),J=1,NJ)
END DO


CLOSE(7)
!--------------------------------
!	COMPUTE CELL CENTERS
!--------------------------------
OPEN(UNIT=8, FILE='centers.dat')
OPEN(UNIT=9, FILE='vertices.dat')
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

      ! CALCULATE THE CENTER OF EACH CELL, STORE IN 1D ARRAY
      CC_X(L) = (PT(I,J,1,1)+PT(I,J,2,1)+PT(I,J,3,1)+PT(I,J,4,1))/4.0  !X
      CC_Y(L) = (PT(I,J,1,2)+PT(I,J,2,2)+PT(I,J,3,2)+PT(I,J,4,2))/4.0  !Y

      CC_MAT(I,J,1) = (PT(I,J,1,1)+PT(I,J,2,1)+PT(I,J,3,1)+PT(I,J,4,1))/4.0  !X
      CC_MAT(I,J,2) = (PT(I,J,1,2)+PT(I,J,2,2)+PT(I,J,3,2)+PT(I,J,4,2))/4.0  !Y

      ! COMPUTE RADIUS
      R = SQRT(CC_X(L)**2.0 + CC_Y(L)**2.0)

      ! WRITE CENTERS TO FILE
      WRITE(8,*) CC_X(L), CC_Y(L), R

      ! Update cell index
      L = L + 1

   END DO
END DO

! WRITE VERTICES TO FILE
DO I=1,NI
	DO J=1,NJ
		WRITE(9,*) X(I,J), Y(I,J)
	END DO
END DO

!---------------------------------
!	COMPUTE FACE CENTERS
!---------------------------------
OPEN(UNIT=10, FILE='faces.dat')
L = 1

! ASSUMING THE WALL LIES ON THE J=1 GRID BOUNDARY
DO I=1,NI-1
    FC_X(L) = (X(I,1) + X(I+1,1))/2.0
    FC_Y(L) = (Y(I,1) + Y(I+1,1))/2.0

    R = SQRT(FC_X(L)**2.0 + FC_Y(L)**2.0)

    WRITE(10,*) FC_X(L), FC_Y(L), R
    L = L + 1
END DO


CC_X_OUT = C_LOC(CC_X(1))
CC_Y_OUT = C_LOC(CC_Y(1))

FC_X_OUT = C_LOC(FC_X(1))
FC_Y_OUT = C_LOC(FC_Y(1))

CLOSE(8)
CLOSE(9)
CLOSE(10)

END SUBROUTINE readmesh