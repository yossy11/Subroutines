      ! calculate inverse matrix
      SUBROUTINE calcInverse(dim,mat,invMat)
      INTEGER dim
      DOUBLE PRECISION mat(dim,dim),invMat(dim,dim),bufMat(dim,dim),buf
      bufMat(:,:) = mat(:,:)
      invMat(:,:) = 0.0D0
      DO i=1,dim
        invMat(i,i) = 1.0D0
      END DO
      DO i=1,dim
        buf = 1.0D0/mat(i,i)
        mat(i,:) = mat(i,:)*buf
        invMat(i,:) = invMat(i,:)*buf
        DO j=1,dim
          IF (i/=j) THEN
            buf = mat(j,i)
            mat(j,:) = mat(j,:) - mat(i,:)*buf
            invMat(j,:) = invMat(j,:) - invMat(i,:)*buf
          END IF
        END DO
      END DO
      mat(:,:) = bufMat(:,:)
      RETURN
      END SUBROUTINE calcInverse


      ! calculate invariants of stress tensor
      SUBROUTINE calcInvariants(stress,invariants)
      DOUBLE PRECISION stress(6),invariants(2)
      invariants(1) = (stress(1) + stress(2) + stress(3))/3.0D0
      invariants(2) = (stress(4)**2 + stress(5)**2 + stress(6)**2 - 
     & stress(1)*stress(2) - stress(2)*stress(3) - 
     & stress(3)*stress(1))/3.0D0
      RETURN
      END SUBROUTINE calcInvariants


      ! calculate principal value of matrix
      SUBROUTINE calcPrincipal(mat,principal,invariants)
      DOUBLE PRECISION mat(3,3),principal(3),invariants(3),radius,
     & theta,PI
      PI = acos(-1.0D0)
      principal(:) = 0.0D0
      invariants(:) = 0.0D0
      DO i=1,3
        invariants(1) = invariants(1) + mat(i,i)/3.0D0
      END DO
      invariants(2) = (mat(2,3)**2 + mat(3,1)**2 + mat(1,2)**2 - 
      & mat(2,2)*mat(3,3) - mat(3,3)*mat(1,1) - mat(1,1)*mat(2,2))/3.0D0
      invariants(3) = (2*mat(1,2)*mat(2,3)*mat(3,1) + 
      & mat(1,1)*mat(2,2)*mat(3,3) - mat(1,1)*mat(2,3)**2 - 
      & mat(2,2)*mat(3,1)**2 - mat(3,3)*mat(1,2)**2)/2.0D0
      radius = 2*sqrt(invariants(1)**2 + invariants(2))
      theta = acos(((2*invariants(1)**3 + 3*invariants(1)*
      & invariants(2) + 2*invariants(3))/2.0D0)/
      & sqrt((invariants(1)**2 + invariants(2))**3))
      principal(1) = radius*cos(theta/3.0D0) + invariants(1)
      principal(2) = radius*cos((theta + 4*PI)/3.0D0) + invariants(1)
      principal(3) = radius*cos((theta + 2*PI)/3.0D0) + invariants(1)
      RETURN
      END SUBROUTINE calcPrincipal


      ! enable using same random seed, this func is just a memo
      FUNCTION random()
      DOUBLE PRECISION x
      INTEGER seedsize,seed
      ALLOCATABLE seed(:)
      CALL RANDOM_SEED(size=seedsize)
      ALLOCATE(seed(seedsize))
      CALL RANDOM_SEED(get=seed)
      CALL RANDOM_NUMBER(x)
      WRITE(*,*) x
      ! CALL sleep(1)
      CALL RANDOM_SEED(put=seed)
      CALL RANDOM_NUMBER(x)
      WRITE(*,*) x
      END FUNCTION random


      ! sort list in descending order
      RECURSIVE SUBROUTINE quickSort(vec, first, last)
      IMPLICIT NONE
      INTEGER first,last,i,j
      DOUBLE PRECISION  vec(*),x,t      
      x = vec( (first+last) / 2 )
      i = first
      j = last
      DO
        DO WHILE (vec(i) > x)
          i=i+1
        END DO
        DO WHILE (x > vec(j))
          j=j-1
        END DO
        IF (i >= j) EXIT
        t = vec(i)
        vec(i) = vec(j)
        vec(j) = t
        i=i+1
        j=j-1
      END DO
      IF (first < i - 1) THEN
        CALL quicksort(vec, first, i - 1)
      END IF
      IF (j + 1 < last) THEN
        CALL quicksort(vec, j + 1, last)
      END IF
      END SUBROUTINE quickSort
