      ! This is UMAT base subroutine made by yossy11
      ! yield function : 
      ! plastic potential : 
      ! hardening rule : 
      ! flow rule : non-AFR
      ! integration algorithm : CCP(convex cutting plane)

      ! Umat interface
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     & DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     & CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     & CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     & DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     & STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     & PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      ! define variables, and their dimensions
      INTEGER i,iterationNum
      DOUBLE PRECISION TOLER,YOUNG,POISSON,lame,shearMod,eStrain(6),
     & pStrain(6),eqpStrain,totalStrain(6),eqStress,eqGStress,
     & flowStress,calc_eqStress,calc_eqGStress,calc_FlowStress,dfds(6),
     & dgds(6),H,calc_H,dLambda,F,lambda

      ! define constants
      PARAMETER(TOLER=1.0D-5,YOUNG=6.9D4,POISSON=0.33D0)

      ! this subroutine can be used only on the condition of NDI=NSHR=3
      IF (NDI/=3 .or. NSHR/=3) THEN
        WRITE(7,*) 'invalid element type, please check about that'
        CALL XIT
      END IF

      ! initialize DDSDDE
      lame = YOUNG*POISSON/((1.0D0 - 2.0D0*POISSON)*(1.0D0 + POISSON))
      shearMod = YOUNG/(2.0D0*(1.0D0 + POISSON))
      DDSDDE(:,:) = 0.0D0
      DDSDDE(1:3,1:3) = lame
      DO i=1,3
        DDSDDE(i,i) = lame + 2*shearMod
        DDSDDE(i+3,i+3) = shearMod
      END DO

      ! read STATEV
      CALL ROTSIG(STATEV(1), DROT, eStrain, 2, NDI, NSHR)
      CALL ROTSIG(STATEV(NTENS+1), DROT, pStrain, 2, NDI, NSHR)
      eqpStrain = STATEV(2*NTENS+1)
      eStrain(:) = eStrain(:) + DSTRAN(:)
      totalStrain(:) = eStrain(:) + pStrain(:)

      ! calculate trial stress
      STRESS(:) = STRESS(:) + MATMUL(DDSDDE,DSTRAN)
      iterationNum = 0
      lambda = 0.0D0

      DO WHILE (iterationNum < 1000000)
        ! calculate eqStress and flowStress
        eqStress = calc_eqStress()
        flowStress = calc_FlowStress()
        F = eqStress - flowStress
        IF (ISNAN(F)) THEN
          WRITE(7,*) "F is NaN"
          CALL XIT
        END IF

        ! if not yield
        IF (F < flowStress*TOLER) THEN
          CALL updateSTATEV(NSTATV,STATEV,eStrain,pStrain,eqpStrain)
          IF (iterationNum /= 0) THEN
            CALL updateDDSDDE()
          END IF
          RETURN
        END IF

        ! if yield
        eqGStress = calc_eqGStress()
        H = calc_H()
        CALL calc_dfds()
        CALL calc_dgds()

        dLambda = F/(DOT_PRODUCT(dfds,MATMUL(DDSDDE,dgds)) + 
     &   H*eqGStress/eqStress)
        IF (dLambda<0) THEN
          WRITE(7,*) "dLambda < 0, invalid!!"
          CALL XIT
        END IF

        lambda = lambda + dLambda

        STRESS(:) = STRESS(:) - dLambda*MATMUL(DDSDDE,dgds)
        eqpStrain = eqpStrain + dLambda*eqGStress/eqStress
        pStrain(:) = pStrain(:) + dLambda*dgds(:)
        eStrain(:) = eStrain(:) - dLambda*dgds(:)
        iterationNum = iterationNum + 1
      END DO
      WRITE(7,*) "not converged!!!"
      WRITE(7,*) "F",F
      CALL XIT
      RETURN
      END SUBROUTINE UMAT


      ! calculate inverse matrix
      SUBROUTINE calc_Inverse(dim,mat,invMat)
      IMPLICIT NONE
      INTEGER dim,i,j
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
      END SUBROUTINE calc_Inverse


      ! calculate dyad of vectors
      SUBROUTINE calc_Dyad(dim,vec1,vec2,mat)
      IMPLICIT NONE
      INTEGER dim,i,j
      DOUBLE PRECISION vec1(dim),vec2(dim),mat(dim,dim)
      DO i=1,dim
        DO j=1,dim
          mat(i,j) = vec1(i)*vec2(j)
        END DO
      END DO
      RETURN
      END SUBROUTINE calc_Dyad


      ! update all state variable
      SUBROUTINE updateSTATEV(NSTATV,STATEV,eStrain,pStrain,eqpStrain)
      IMPLICIT NONE
      INTEGER NSTATV
      DOUBLE PRECISION STATEV(NSTATV),eStrain(6),pStrain(6),eqpStrain
      STATEV(1:6)=eStrain
      STATEV(7:12)=pStrain
      STATEV(13)=eqpStrain
      RETURN
      END SUBROUTINE updateSTATEV


      ! update tangent modulus
      SUBROUTINE updateDDSDDE()
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION eqStress,eqGStress,calc_eqStress,calc_eqGStress,
     & H,calc_H,dfds(6),dgds(6),ddgdds(6,6),subDDSDDE(6,6),DDSDDE(6,6),
     & invDDSDDE(6,6),A(7,7),B(6,6),invB(6,6),C(7,7),lambda,h0,
     & subVec(7),vec1(7),vec2(7),STRESS(6),eqpStrain

      eqStress = calc_eqStress()
      eqGStress = calc_eqGStress()
      H = calc_H()
      CALL calc_dfds()
      CALL calc_dgds()
      CALL calc_ddgdds()

      subDDSDDE(:,:) = DDSDDE(:,:)
      CALL calc_Inverse(6,DDSDDE,invDDSDDE)
      invB(:,:) = invDDSDDE(:,:) + lambda*ddgdds(:,:)
      CALL calc_Inverse(6,invB,B)
      A(:,:) = 0.0D0
      A(1:6,1:6) = B(:,:)
      A(7,7) = 1.0D0
      h0 = eqGStress/eqStress
      subVec(:) = 0.0D0
      subVec(1:6) = dgds(:)
      subVec(7) = -1.0D0*h0
      vec1 = MATMUL(A,subVec)
      subVec(1:6) = dfds(:)
      subVec(7) = 0.0D0
      vec2 = MATMUL(subVec,A)
      CALL calc_Dyad(7,vec1,vec2,C)
      C(:,:) = C(:,:)/(DOT_PRODUCT(dfds,MATMUL(B,dgds)) + H*h0)
      A(:,:) = A(:,:) - C(:,:)
      DDSDDE(:,:) = A(1:6,1:6)
      DO i=1,6
        DO j=1,6
          IF (ISNAN(DDSDDE(i,j))) THEN
            DDSDDE(:,:) = subDDSDDE(:,:)
            RETURN
          END if
        END DO
      END DO
      RETURN
      END SUBROUTINE updateDDSDDE


      ! calculate equivalent stress using yield function
      DOUBLE PRECISION FUNCTION calc_eqStress()
      IMPLICIT NONE
      calc_eqStress = 0.0D0
      RETURN
      END FUNCTION calc_eqStress


      ! calculate equivalent stress using plastic potential
      DOUBLE PRECISION FUNCTION calc_eqGStress()
      IMPLICIT NONE
      calc_eqGStress = 0.0D0
      RETURN
      END FUNCTION calc_eqGStress


      ! calculate flow stress
      DOUBLE PRECISION FUNCTION calc_FlowStress()
      IMPLICIT NONE
      calc_FlowStress = 0.0D0
      RETURN
      END FUNCTION calc_FlowStress


      ! calculate first order differential of flow stress
      DOUBLE PRECISION FUNCTION calc_H()
      IMPLICIT NONE
      calc_H = 0.0D0
      RETURN
      END FUNCTION calc_H


      ! calculate first order differentials of yield function
      SUBROUTINE calc_dfds()
      IMPLICIT NONE
      RETURN
      END SUBROUTINE calc_dfds


      ! calculate first order differentials of plastic potential
      SUBROUTINE calc_dgds()
      IMPLICIT NONE
      RETURN
      END SUBROUTINE calc_dgds


      ! calculate second order differentials of plastic potential
      SUBROUTINE calc_ddgdds()
      IMPLICIT NONE
      RETURN
      END SUBROUTINE calc_ddgdds


      