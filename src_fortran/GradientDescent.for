      PROGRAM GradientDescent
      IMPLICIT NONE
      INTEGER n,i,iteCount
      DOUBLE PRECISION TOL,strain,stress,strains,stresses,params(3),
     & errorFunc,error
      ALLOCATABLE strains(:),stresses(:)
      PARAMETER(TOL=1.0D-5)
      n = 0
      OPEN(20, FILE='./Datas/exData.csv', STATUS='old')
      READ(20, '()')
      DO
        READ(20, *, END=100)
        n = n + 1
      END DO
  100 continue
      REWIND(20)
      ALLOCATE(strains(n),stresses(n))
      strains(:) = 0.0D0
      stresses(:) = 0.0D0
      READ(20, '()')
      DO i = 1, n
        READ(20, *) strain,stress
        strains(i) = strain
        stresses(i) = stress
      END DO
      params(1) = 615.3D0
      params(2) = 7.61D-3
      params(3) = 0.363D0

      error = errorFunc(n,strains,stresses,params)
      iteCount = 0
      DO WHILE (error>TOL)
        WRITE(*,*) error
        CALL optimizeParams(n,strains,stresses,params)
        error = errorFunc(n,strains,stresses,params)
        iteCount = iteCount + 1
        IF (iteCount>=1000) THEN
          WRITE(*,*) 'too many'
          STOP
        END IF
      END DO
      END PROGRAM GradientDescent


      DOUBLE PRECISION FUNCTION func(params,eqpStrain)
      IMPLICIT NONE
      DOUBLE PRECISION params(3),eqpStrain
      func = params(1)*(params(2) + eqpStrain)**params(3)
      IF(ISNAN(func)) THEN
        WRITE(*,*) 'damedeshita'
        STOP
      END IF
      RETURN
      END FUNCTION func


      DOUBLE PRECISION FUNCTION errorFunc(n,strains,stresses,params)
      IMPLICIT NONE
      INTEGER n,i
      DOUBLE PRECISION strains(n),stresses(n),params(3),WEIGHT,func
      errorFunc = 0.0D0
      DO i=1,n
        errorFunc = 
     &   errorFunc + (func(params,strains(i))/stresses(i))**2
      END DO
      RETURN
      END FUNCTION errorFunc


      ! optimize parameters, have to code for each function
      SUBROUTINE optimizeParams(n,strains,stresses,params)
      IMPLICIT NONE
      INTEGER n,i
      DOUBLE PRECISION strains(n),stresses(n),stress,params(3),
     & dParams(3),func,multiplier
      
      DO i=1,n
        stress = func(params,strains(i))
        multiplier = 2.0D0*(stress/stresses(i))/stresses(i)
        dParams(1) = dParams(1) + 
     &   multiplier*(params(2) + strains(i))**params(3)
        dParams(2) = dParams(2) + multiplier*params(1)*params(3)*
     &   (params(2) + strains(i))**(params(3) - 1.0D0)
        dParams(3) = 
     &   dParams(3) + multiplier*stress*LOG(params(2) + strains(i))
      END DO
      params(:) = params(:) - dParams(:)
      RETURN
      END SUBROUTINE optimizeParams