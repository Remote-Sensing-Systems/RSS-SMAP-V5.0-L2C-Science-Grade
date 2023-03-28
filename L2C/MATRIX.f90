! ROUTINES FOR MATRIX INVERSION AND SOLVING LINEAR EQUATIONS
! THOMAS MEISSNER
! RSS
! MARCH 2005

MODULE MATRIX_ROUTINES
   IMPLICIT NONE
   PUBLIC
   SAVE


CONTAINS

   SUBROUTINE INVERT_MAT(NDIM,AMAT,  AINV,ISING)
! INVERT A(NDIM,NDIM)
! LU DECOMPOSITION
! ISING = 1: FLAG FOR SINGULAR MATRIX
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NDIM
      REAL(8), DIMENSION(NDIM,NDIM), INTENT(IN)     :: AMAT
      REAL(8), DIMENSION(NDIM,NDIM), INTENT(OUT)    :: AINV
      INTEGER(4), INTENT(OUT)                       :: ISING  !=0 o.k.    =1 singular
      REAL(8), DIMENSION(NDIM,NDIM) :: UNITY
      INTEGER(4) :: I

      UNITY = 0.0
      DO I = 1,NDIM
         UNITY(I,I) = 1.0
      ENDDO

      CALL SOLVE_LINEQ(NDIM,NDIM,AMAT,UNITY,   AINV,ISING)
      RETURN
   END  SUBROUTINE INVERT_MAT




   SUBROUTINE SOLVE_LINEQ(NDIM,MDIM,AMAT,BMAT,   XMAT,ISING)
! SOLVE LINEAR SYSTEM: A * X = B
! LU DECOMPOSITION
! ISING = 1: FLAG FOR SINGULAR SYSTEM

      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NDIM, MDIM
      REAL(8), DIMENSION(NDIM,NDIM), INTENT(IN)     :: AMAT
      REAL(8), DIMENSION(NDIM,MDIM), INTENT(IN)     :: BMAT
      REAL(8), DIMENSION(NDIM,MDIM), INTENT(OUT)    :: XMAT
      INTEGER(4), INTENT(OUT)                       :: ISING  !=0 o.k.    =1 singular


      REAL(8), DIMENSION(NDIM,NDIM):: LMAT, UMAT
      INTEGER(4), DIMENSION(NDIM)  :: IPIVOT
      REAL(8), DIMENSION(NDIM     ):: BVEC, CVEC, XVEC
      INTEGER(4) :: I,K

      CALL LU_DEC(NDIM,AMAT,  LMAT,UMAT,IPIVOT,ISING)
      IF (ISING == 1) RETURN


      DO I = 1,MDIM
         BVEC(:) = BMAT(:,I)
         DO K=1,NDIM
            BVEC(K) = BMAT(IPIVOT(K),I)
         ENDDO
         CALL LSOLVE(NDIM,LMAT,BVEC,  CVEC)
         CALL USOLVE(NDIM,UMAT,CVEC,  XVEC)
         XMAT (:,I) = XVEC(:)
      ENDDO

      RETURN

   END SUBROUTINE SOLVE_LINEQ




   SUBROUTINE LU_DEC(NDIM,AMAT,  LMAT,UMAT,IPERMUT,IERR)
! LU DECOMPOSITION W/I PIVOT
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NDIM
      REAL(8), DIMENSION(NDIM,NDIM), INTENT(IN)  :: AMAT
      REAL(8), DIMENSION(NDIM,NDIM), INTENT(OUT) :: LMAT, UMAT
      INTEGER(4), DIMENSION(NDIM), INTENT(OUT)   :: IPERMUT
      INTEGER(4), INTENT(OUT)                    :: IERR


      REAL(8), DIMENSION(NDIM,NDIM) :: ACOP
      INTEGER(4) :: I,J,K,L,IPIV
      REAL(8) :: SUM

      REAL(8) :: XMAX, YMAX, Z
      REAL(8), PARAMETER :: SMALL = 1.0E-15

      IERR = 0
      ACOP = AMAT
      DO I = 1,NDIM
         IPERMUT(I) = I
      ENDDO



      COLUMNS: DO  J = 1,NDIM

         ROWS1: DO  I = 1,J-1
            SUM = ACOP(I,J)
            DO K=1,I-1
               SUM = SUM - ACOP(I,K)*ACOP(K,J)
            ENDDO
            ACOP(I,J) = SUM
         ENDDO ROWS1

         XMAX = 0.0
         ROWS2: DO  I = J,NDIM
            SUM = ACOP(I,J)
            DO K =1,J-1
               SUM = SUM - ACOP(I,K)*ACOP(K,J)
            ENDDO
            ACOP(I,J) = SUM

            YMAX = ABS(SUM)  ! SEARCH FOR BEST (=LARGEST) PIVOT
            IF (YMAX >= XMAX) THEN
               IPIV = I
               XMAX = YMAX
            ENDIF
         ENDDO ROWS2


         IF (J /= IPIV) THEN !CURRENT ELEMENT IS NOT THE OPTIMAL PIVOT
            !THEREFORE INTERCHANGE ROWS
            DO K=1,NDIM
               Z=ACOP(IPIV,K)
               ACOP(IPIV,K) = ACOP(J,K)
               ACOP(J,K) =Z
            ENDDO
         ENDIF
         K = IPERMUT(J)
         L = IPERMUT(IPIV)
         IPERMUT(J) = L
         IPERMUT(IPIV) = K
         ! ROW PERMUTATION INDEX

         IF (J /= NDIM) THEN
            ROWS3: DO  I=J+1,NDIM
               ACOP(I,J)=ACOP(I,J)/ACOP(J,J)
            ENDDO ROWS3
         ENDIF

      ENDDO COLUMNS


      UMAT = 0.0
      LMAT = 0.0
      DO  J=1,NDIM
         IF (ABS(ACOP(J,J)) < SMALL) THEN
            IERR = 1
            RETURN
         ENDIF
         DO I=1,NDIM
            IF (I<=J) UMAT(I,J)=ACOP(I,J)
            IF (J< I) LMAT(I,J)=ACOP(I,J)
            IF (I==J) LMAT(I,J) = 1.0
         ENDDO
      ENDDO
      RETURN
   END SUBROUTINE LU_DEC



   SUBROUTINE LSOLVE(NDIM,LMAT,BVEC,  CVEC)  ! L*c = b
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NDIM
      REAL(8), DIMENSION(NDIM,NDIM), INTENT(IN)  :: LMAT
      REAL(8), DIMENSION(NDIM)     , INTENT(IN)  :: BVEC
      REAL(8), DIMENSION(NDIM)     , INTENT(OUT) :: CVEC
      INTEGER(4) :: I

      CVEC(1) = BVEC(1)/LMAT(1,1)
      DO I = 2,NDIM
         CVEC(I) = (BVEC(I) - DOT_PRODUCT(LMAT(I,1:I-1),CVEC(1:I-1)))/ LMAT(I,I)
      ENDDO
      RETURN

   END SUBROUTINE LSOLVE



   SUBROUTINE USOLVE(NDIM,UMAT,CVEC,  YVEC)  !U*y = c
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NDIM
      REAL(8), DIMENSION(NDIM,NDIM), INTENT(IN)  :: UMAT
      REAL(8), DIMENSION(NDIM)     , INTENT(IN)  :: CVEC
      REAL(8), DIMENSION(NDIM)     , INTENT(OUT) :: YVEC
      INTEGER(4) :: I


      YVEC(NDIM) = CVEC(NDIM)/UMAT(NDIM,NDIM)
      DO I = NDIM-1,1,-1
         YVEC(I) = (CVEC(I) - DOT_PRODUCT(UMAT(I,I+1:NDIM),YVEC(I+1:NDIM)))/ UMAT(I,I)
      ENDDO
      RETURN

   END SUBROUTINE USOLVE


   SUBROUTINE INDEXINVERT(NDIM,INDX,   JNDX)
      INTEGER(4), INTENT(IN)    :: NDIM
      INTEGER(4), DIMENSION(NDIM), INTENT(IN)    :: INDX ! PERMUTATION OF 1,2, ... NDIM
      INTEGER(4), DIMENSION(NDIM), INTENT(OUT)   :: JNDX ! INVERSE PERMUTATION
      INTEGER(4) :: I,K

      DO  I =1,NDIM
         K = INDX(I)
         JNDX(K)=I
      ENDDO

   END SUBROUTINE INDEXINVERT




   SUBROUTINE ROWINTERCHANGE(NDIM,MDIM,IPIVOT,XMAT)
      IMPLICIT NONE
      INTEGER(4), INTENT(IN)    :: NDIM,MDIM ! ROW COL DIM
      INTEGER(4), DIMENSION(NDIM), INTENT(IN)    :: IPIVOT
      REAL(8)   , DIMENSION(NDIM,MDIM), INTENT(INOUT) :: XMAT
      INTEGER(4) :: I
      REAL(8), DIMENSION(MDIM) :: XVEC

      DO  I = 1,NDIM
         XVEC(:) = XMAT(IPIVOT(I),:)
         XMAT(I,:) = XVEC(:)
      ENDDO
      RETURN
   END SUBROUTINE ROWINTERCHANGE


END MODULE MATRIX_ROUTINES
