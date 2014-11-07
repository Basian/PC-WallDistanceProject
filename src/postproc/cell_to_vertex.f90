SUBROUTINE cell_to_vertex (NI,NJ,X,Y,CT,VAR_C,VAR_V)

!---------------------------------------------------
!---------------------------------------------------
!  INTERPOLATES VARIABLES TO VERTICES FROM SURROUNDING 
!  CELL CENTERS FOR OUTPUT TO FILE
!---------------------------------------------------
!---------------------------------------------------
IMPLICIT NONE

INTEGER :: M,N,I,J,K,L
INTEGER, INTENT(IN) :: NI,NJ
REAL*8, DIMENSION(NI,NJ), INTENT(IN) :: X,Y
REAL*8, DIMENSION(NI-1,NJ-1,2), INTENT(IN) :: CT
REAL*8, DIMENSION(NI-1,NJ-1), INTENT(IN) :: VAR_C
REAL*8, DIMENSION(NI,NJ,5), INTENT(OUT) :: VAR_V
REAL*8, DIMENSION(5) :: QSUM
REAL*8, DIMENSION(4) :: R
REAL*8 :: RSUM,TEST

DO I=1,NI
   DO J=1,NJ

      RSUM=0.0
      DO K=1,5
         QSUM(K)=0.0
      END DO

      IF (I==1 .AND. J==1) THEN
         DO K=1,5
            VAR_V(I,J,K) = VAR_C(I,J)      !Corner one
         END DO
      ELSE IF (I==NI .AND. J==1) THEN
         DO K=1,5
            VAR_V(I,J,K) = VAR_C(I-1,1)    !Corner two
         END DO
      ELSE IF (I==NI .AND. J==NJ) THEN
         DO K=1,5
            VAR_V(I,J,K) = VAR_C(I-1,J-1)  !Corner three
         END DO
      ELSE IF (I==1 .AND. J==NJ) THEN
         DO K=1,5
            VAR_V(I,J,K) = VAR_C(I,J-1)    !Corner four
         END DO
      ELSE IF (I==1 .AND. J/=1 .AND. J/=NJ) THEN   !Left boundary 
         R(1) = SQRT((X(I,J)-CT(I,J-1,1))**2.0 + (Y(I,J)-CT(I,J-1,2))**2.0)
         R(2) = SQRT((X(I,J)-CT(I,J,1))**2.0 + (Y(I,J)-CT(I,J,2))**2.0)
         RSUM = 1/R(1) + 1/R(2)

         DO L=1,5
            QSUM(L) = 0.0
         END DO

         DO K=1,5
            QSUM(K) = QSUM(K) + VAR_C(I,J-1)/R(1)
            VAR_V(I,J,K) = (QSUM(K) + VAR_C(I,J)/R(2))/RSUM
         END DO
      ELSE IF (J==1 .AND. I/=1 .AND. I/=NI) THEN   !Bottom boundary
         R(1) = SQRT((X(I,J)-CT(I-1,J,1))**2.0 + (Y(I,J)-CT(I-1,J,2))**2.0)
         R(2) = SQRT((X(I,J)-CT(I,J,1))**2.0 + (Y(I,J)-CT(I,J,2))**2.0)
         RSUM = 1/R(1) + 1/R(2)

         DO L=1,5
            QSUM(L) = 0.0
         END DO

         DO K=1,5
            QSUM(K) = QSUM(K) + VAR_C(I-1,J)/R(1)
            VAR_V(I,J,K) = (QSUM(K) + VAR_C(I,J)/R(2))/RSUM
         END DO
      ELSE IF (I==NI .AND. J/=1 .AND. J/=NJ) THEN     !Right boundary
         R(1) = SQRT((X(I,J)-CT(I-1,J-1,1))**2.0 + (Y(I,J)-CT(I-1,J-1,2))**2.0)
         R(2) = SQRT((X(I,J)-CT(I-1,J,1))**2.0 + (Y(I,J)-CT(I-1,J,2))**2.0)
         RSUM = 1/R(1) + 1/R(2)

         DO L=1,5
            QSUM(L) = 0.0
         END DO

         DO K=1,5
            QSUM(K) = QSUM(K) + VAR_C(I-1,J-1)/R(1)
            VAR_V(I,J,K) = (QSUM(K) + VAR_C(I-1,J)/R(2))/RSUM
         END DO
      ELSE IF (J==NJ .AND. I/=1 .AND. I/=NI) THEN     !Top boundary
         R(1) = SQRT((X(I,J)-CT(I-1,J-1,1))**2.0 + (Y(I,J)-CT(I-1,J-1,2))**2.0)
         R(2) = SQRT((X(I,J)-CT(I,J-1,1))**2.0 + (Y(I,J)-CT(I,J-1,2))**2.0)
         RSUM = 1/R(1) + 1/R(2)

         DO L=1,5
            QSUM(L) = 0.0
         END DO

         DO K=1,5
            QSUM(K) = QSUM(K) + VAR_C(I-1,J-1)/R(1)
            VAR_V(I,J,K) = (QSUM(K) + VAR_C(I,J-1)/R(2))/RSUM
         END DO
      ELSE                                            !Interior points
         R(1) = SQRT((X(I,J)-CT(I-1,J-1,1))**2.0 + (Y(I,J)-CT(I-1,J-1,2))**2.0)
         R(2) = SQRT((X(I,J)-CT(I,J-1,1))**2.0 + (Y(I,J)-CT(I,J-1,2))**2.0)
         R(3) = SQRT((X(I,J)-CT(I,J,1))**2.0 + (Y(I,J)-CT(I,J,2))**2.0)
         R(4) = SQRT((X(I,J)-CT(I-1,J,1))**2.0 + (Y(I,J)-CT(I-1,J,2))**2.0)
         RSUM = 1/R(1) + 1/R(2) + 1/R(3) + 1/R(4)

         DO K=1,5
            QSUM(K) = QSUM(K) + VAR_C(I-1,J-1)/R(1)
            QSUM(K) = QSUM(K) + VAR_C(I,J-1)/R(2)
            QSUM(K) = QSUM(K) + VAR_C(I,J)/R(3)
            VAR_V(I,J,K) = (QSUM(K) + VAR_C(I-1,J)/R(4))/RSUM
         END DO
      END IF

   END DO
END DO

END SUBROUTINE cell_to_vertex








































