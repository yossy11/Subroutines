      ! This is UMAT subroutine made by yossy11
      ! yield function : Hill48
      ! plastic potential : Hill48
      ! hardening rule : Swift
      ! flow rule : non-AFR
      ! integration algorithm : CBE(classical backward euler)
      ! ToDo: fix b0(variable used in newton-raphson iteration scheme)

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
      INTEGER i,iterationNum,numSubSteps
      DOUBLE PRECISION TOLER,YOUNG,POISSON,HARDK,HARDN,HARDSTRAIN0,lame,
     & shearMod,eStrain(6),pStrain(6),eqpStrain,totalStrain(6),
     & hillFParams(4),hillGParams(4),eqStress,calc_eqStress,eqGStress,
     & calc_eqGStress,flowStress,calc_FlowStress,F,lambda,
     & trialStress(6),trialeqpStrain,dGdS(6)

      ! define constants
      PARAMETER(TOLER=1.0D-5,YOUNG=6.9D4,POISSON=0.33D0,HARDK=646.0D0,
     & HARDN=0.227D0,HARDSTRAIN0=2.5D-2,YLDM=6,numSubSteps=100)

      ! anisotropic params
      hillFParams(1) = 0.570281847241541
      hillFParams(2) = 0.3632288531248622
      hillFParams(3) = 0.6367711468751378
      hillFParams(4) = 2.5710503236965727

      hillGParams(1) = 0.25216953733566727
      hillGParams(2) = 0.8254230293025175
      hillGParams(3) = 0.17457697069748246
      hillGParams(4) = 2.2380520016508463

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
      trialStress(:) = STRESS(:)
      trialeqpStrain = eqpStrain
      iterationNum = 0
      lambda = 0.0D0

      ! calculate eqStress and flowStress
      eqStress = calc_eqStress(hillFParams,STRESS)
      eqGStress = calc_eqGStress(hillGParams,STRESS)
      flowStress = calc_FlowStress(HARDK,HARDSTRAIN0,HARDN,eqpStrain)
      F = eqStress - flowStress
      IF (ISNAN(F)) THEN
        WRITE(7,*) "F is NaN"
        CALL XIT
      END IF

      ! return mapping method
      DO WHILE (iterationNum < numSubSteps)
        IF (iterationNum==numSubSteps-1) THEN
          WRITE(7,*) "not converged"
          CALL XIT
        END IF

        dGdS(:) = 0.0D0
        CALL calc_dGdS(hillGParams,STRESS,dGdS)

        ! if not yield
        IF (F < flowStress*TOLER) THEN
          EXIT
        ! if yield
        ELSE
          STRESS(:) = trialStress(:) - lambda*MATMUL(DDSDDE,dGdS)
          eqpStrain = trialeqpStrain + lambda*eqGStress/eqStress
        END IF
        
        CALL newton_raphson(DDSDDE,STRESS,trialStress,hillFParams,
     & hillGParams,HARDK,HARDN,HARDSTRAIN0,eqpStrain,trialeqpStrain,
     & lambda,TOLER)
        eqStress = calc_eqStress(hillFParams,STRESS)
        eqGStress = calc_eqGStress(hillGParams,STRESS)
        flowStress = calc_FlowStress(HARDK,HARDSTRAIN0,HARDN,eqpStrain)
        F = eqStress - flowStress
        IF (ISNAN(F)) THEN
          WRITE(7,*) "F is NaN"
          CALL XIT
        END IF
        iterationNum = iterationNum + 1
      END DO

      CALL calc_dGdS(hillGParams,STRESS,dGdS)
      pStrain(:) = pStrain(:) + lambda*dGdS(:)
      eStrain(:) = eStrain(:) - lambda*dGdS(:)
      CALL updateSTATEV(NSTATV,STATEV,eStrain,pStrain,eqpStrain)
      IF (iterationNum/=0) THEN
        CALL updateDDSDDE(hillFParams,hillGParams,STRESS,DDSDDE,lambda,
     &   eqStress,HARDK,HARDN,HARDSTRAIN0,eqpStrain)
      END IF
      RETURN
      END SUBROUTINE UMAT


      ! Newton-Raphson iteration
      SUBROUTINE newton_raphson(DDSDDE,STRESS,trialStress,hillFParams,
     & hillGParams,HARDK,HARDN,HARDSTRAIN0,eqpStrain,trialeqpStrain,
     & lambda,TOLER)
      INCLUDE 'ABA_PARAM.INC'
      INTEGER NRIterationNum,i,j
      DOUBLE PRECISION DDSDDE(6,6),STRESS(6),trialStress(6),
     & hillFParams(4),hillGParams(4),HARDK,HARDN,HARDSTRAIN0,eqpStrain,
     & trialeqpStrain,lambda,NRTOLER,invDDSDDE(6,6),eqStress,
     & calc_eqStress,eqGStress,calc_eqGStress,flowStress,
     & calc_FlowStress,initialF,F,dfdS(6),dGdS(6),ddGddS(6,6),H,a0(6),
     & b0,A(7,7),invA(7,7),vec1(7),vec2(7),vec3(7),numerator,
     & denominator,dLambda,increment(7),TOLER
      PARAMETER(NRTOLER=5.0D-1)
      CALL calc_Inverse(6,DDSDDE,invDDSDDE)
      eqStress = calc_eqStress(hillFParams,STRESS)
      eqGStress = calc_eqGStress(hillGParams,STRESS)
      flowStress = calc_FlowStress(HARDK,HARDSTRAIN0,HARDN,eqpStrain)
      initialF = eqStress - flowStress
      F = initialF
      IF (ISNAN(F)) THEN
        WRITE(7,*) "initialF is NaN"
        CALL XIT
      END IF
      NRIterationNum = 0
      DO WHILE (F>initialF*NRTOLER)
        IF (NRIterationNum>100000) THEN
          WRITE(7,*) "NRiteration not converged"
          CALL XIT
        END IF
        ! calculate differentials
        CALL calc_dfdS(hillFParams,STRESS,dfdS)
        CALL calc_dGdS(hillGParams,STRESS,dGdS)
        CALL calc_ddGddS(hillGParams,STRESS,dGdS,ddGddS)
        H = HARDK*HARDN*((HARDSTRAIN0 + eqpStrain)**(HARDN-1.0D0))

        ! calculate dLambda
        a0(:) = MATMUL(invDDSDDE,STRESS(:)-trialStress(:)) + 
     &   lambda*dGdS(:)
        DO i=1,6
          IF (ISNAN(a0(i))) THEN
            WRITE(7,*) "invalid a0"
            CALL XIT
          END IF
        END DO
        ! b0 = trialeqpStrain - eqpStrain + lambda*eqGStress/eqStress
        b0 = trialeqpStrain - eqpStrain + lambda

        invA(:,:) = 0.0D0
        invA(1:6,1:6) = invDDSDDE(:,:) + lambda*ddGddS(:,:)
        invA(7,7) = -1.0D0
        CALL calc_Inverse(7,invA,A)
        DO i=1,7
          DO j=1,7
            IF (ISNAN(A(i,j))) THEN
              WRITE(7,*) "invalid A"
              CALL XIT
            END IF
          END DO
        END DO
        vec1(1:6) = dfdS(:)
        vec1(7) = -1.0D0*H
        vec2(1:6) = a0(:)
        vec2(7) = b0
        numerator = F - DOT_PRODUCT(vec1,MATMUL(A,vec2))
        vec3(1:6) = dGdS(:)
        ! vec3(7) = eqGStress/eqStress
        vec3(7) = 1.0D0

        denominator = DOT_PRODUCT(vec1,MATMUL(A,vec3))
        dLambda = numerator/denominator
        IF (ISNAN(dLambda).or.dLambda<0) THEN
          WRITE(7,*) "invalid dLambda",dLambda,numerator,denominator
          WRITE(7,*) "vec1",vec1
          WRITE(7,*) "vec2",vec2
          WRITE(7,*) "vec3",vec3
          WRITE(7,*) "STRESS",STRESS
          WRITE(7,*) "lambda",lambda
          WRITE(7,*) "eqpStrain",eqpStrain
          WRITE(7,*) "eqStress",eqStress
          WRITE(7,*) "eqGStress",eqGStress
          WRITE(7,*) "flowStress",flowStress
          WRITE(7,*) "F",F
          WRITE(7,*) "ddGddS",ddGddS
          WRITE(7,*) "A",A
          WRITE(7,*) "NRIterationNum",NRIterationNum
          CALL XIT
        END IF
        lambda = lambda + dLambda

        ! update STRESS and eqpStrain
        vec2(7) = 0.0D0
        increment(:) = -1.0D0*MATMUL(A,vec2) - dLambda*MATMUL(A,vec3)
        DO i=1,7
          IF (ISNAN(increment(i))) THEN
            WRITE(7,*) "invalid increment",i
            CALL XIT
          END IF
        END DO
        STRESS(:) = STRESS(:) + increment(1:6)
        eqpStrain = eqpStrain + increment(7)
        IF (increment(7)<0) THEN
          WRITE(7,*) "increment of eqpStrain is invalid",increment(7)
          WRITE(7,*) "dLambda",dLambda
          WRITE(7,*) "left",MATMUL(A,vec2)
          WRITE(7,*) "right",MATMUL(A,vec3)
          WRITE(7,*) increment
          CALL XIT
        END IF

        ! prepare for next loop
        eqStress = calc_eqStress(hillFParams,STRESS)
        eqGStress = calc_eqGStress(hillGParams,STRESS)
        flowStress = calc_FlowStress(HARDK,HARDSTRAIN0,HARDN,eqpStrain)
        F = eqStress - flowStress
        IF (ISNAN(F)) THEN
          WRITE(7,*) "F is NaN"
          CALL XIT
        END IF
        IF (F>initialF) THEN
          WRITE(7,*) "NRiteration diverged",NRIterationNum
          CALL XIT
        END IF
        ! if not yield
        IF (F < flowStress*TOLER) THEN
          EXIT
        END IF
        NRIterationNum = NRIterationNum + 1
      END DO
      RETURN
      END SUBROUTINE newton_raphson


      ! calculate equivalent stress
      DOUBLE PRECISION FUNCTION calc_eqStress(hillParams,STRESS)
      IMPLICIT NONE
      DOUBLE PRECISION hillParams(4),STRESS(6),numerator
      numerator = hillParams(1)*(STRESS(2) - STRESS(3))**2 + 
     & hillParams(2)*(STRESS(3)-STRESS(1))**2 + 
     & hillParams(3)*(STRESS(1)-STRESS(2))**2 + 
     & 2.0D0*hillParams(4)*SUM(STRESS(4:6)**2)
      calc_eqStress = SQRT(1.5D0*numerator/SUM(hillParams(1:3)))
      RETURN
      END FUNCTION calc_eqStress


      ! calculate equivalent stress with plastic potential G
      DOUBLE PRECISION FUNCTION calc_eqGStress(hillParams,STRESS)
      IMPLICIT NONE
      DOUBLE PRECISION hillParams(4),STRESS(6),numerator
      numerator = hillParams(1)*(STRESS(2) - STRESS(3))**2 + 
     & hillParams(2)*(STRESS(3)-STRESS(1))**2 + 
     & hillParams(3)*(STRESS(1)-STRESS(2))**2 + 
     & 2.0D0*hillParams(4)*SUM(STRESS(4:6)**2)
      calc_eqGStress = SQRT(1.5D0*numerator/SUM(hillParams(1:3)))
      RETURN
      END FUNCTION calc_eqGStress


      ! calculate flow stress
      DOUBLE PRECISION FUNCTION calc_FlowStress(HARDK,HARDSTRAIN0,
     & HARDN,eqpStrain)
      IMPLICIT NONE
      DOUBLE PRECISION HARDK,HARDSTRAIN0,HARDN,eqpStrain
      calc_FlowStress = HARDK*(HARDSTRAIN0 + eqpStrain)**HARDN
      RETURN
      END FUNCTION calc_FlowStress


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
      SUBROUTINE updateDDSDDE(hillFParams,hillGParams,STRESS,DDSDDE,
     & lambda,eqStress,HARDK,HARDN,HARDSTRAIN0,eqpStrain)
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION hillFParams(4),hillGParams(4),STRESS(6),
     & DDSDDE(6,6),lambda,eqStress,HARDK,HARDN,HARDSTRAIN0,eqpStrain,
     & eqGStress,calc_eqGStress,dfdS(6),dGdS(6),ddGddS(6,6),A(7,7),
     & B(6,6),invB(6,6),invDDSDDE(6,6),h0,subVec(7),vec1(7),vec2(7),
     & H,C(7,7),subDDSDDE(6,6)
      subDDSDDE(:,:) = DDSDDE(:,:)
      eqGStress = calc_eqGStress(hillGParams,STRESS)
      CALL calc_dfdS(hillFParams,STRESS,dfdS)
      CALL calc_dGdS(hillGParams,STRESS,dGdS)
      CALL calc_ddGddS(hillGParams,STRESS,dGdS,ddGddS)
      CALL calc_Inverse(6,DDSDDE,invDDSDDE)
      invB(:,:) = invDDSDDE(:,:) + lambda*ddGddS(:,:)
      CALL calc_Inverse(6,invB,B)
      A(:,:) = 0.0D0
      A(1:6,1:6) = B(:,:)
      A(7,7) = 1.0D0
      h0 = eqGStress/eqStress
      subVec(:) = 0.0D0
      subVec(1:6) = dGdS(:)
      subVec(7) = -1.0D0*h0
      vec1 = MATMUL(A,subVec)
      subVec(1:6) = dfdS(:)
      subVec(7) = 0.0D0
      vec2 = MATMUL(subVec,A)
      CALL calc_Dyad(7,vec1,vec2,C)
      H = HARDK*HARDN*((HARDSTRAIN0 + eqpStrain)**(HARDN-1.0D0))
      C(:,:) = C(:,:)/
     & (DOT_PRODUCT(dfdS,MATMUL(B,dGdS)) + H*h0)
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


      ! calculate first derivatives of plastic potential
      SUBROUTINE calc_dGdS(hillParams,STRESS,dGdS)
      IMPLICIT NONE
      DOUBLE PRECISION hillParams(4),STRESS(6),dGdS(6),eqGStress,
     & calc_eqGStress,multiplier
      eqGStress = calc_eqGStress(hillParams,STRESS)
      multiplier = 0.75/(SUM(hillParams(1:3))*eqGStress)
      dGdS(1) = -1.0D0*2*hillParams(2)*(STRESS(3)-STRESS(1))+
     & 2*hillParams(3)*(STRESS(1)-STRESS(2))
      dGdS(2) = 2*hillParams(1)*(STRESS(2)-STRESS(3))-
     & 2*hillParams(3)*(STRESS(1)-STRESS(2))
      dGdS(3) = -1.0D0*2*hillParams(1)*(STRESS(2)-STRESS(3))+
     & 2*hillParams(2)*(STRESS(3)-STRESS(1))
      dGdS(4) = 4*hillParams(4)*STRESS(4)
      dGdS(5) = 4*hillParams(4)*STRESS(5)
      dGdS(6) = 4*hillParams(4)*STRESS(6)
      dGdS(:) = dGdS(:) * multiplier
      RETURN
      END SUBROUTINE calc_dGdS


      ! calculate second derivatives of plastic potential
      SUBROUTINE calc_ddGddS(hillParams,STRESS,dGdS,ddGddS)
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION hillParams(4),STRESS(6),dGdS(6),ddGddS(6,6),
     & eqGStress,calc_eqGStress,multiplier
      eqGStress = calc_eqGStress(hillParams,STRESS)
      multiplier = 0.75/(SUM(hillParams(1:3))*eqGStress)

      ddGddS(:,:) = 0.0D0
      ddGddS(1,1) = 2*hillParams(2) + 2*hillParams(3)
      ddGddS(2,2) = 2*hillParams(1) + 2*hillParams(3)
      ddGddS(3,3) = 2*hillParams(1) + 2*hillParams(2)
      ddGddS(4,4) = 4*hillParams(4)
      ddGddS(5,5) = 4*hillParams(4)
      ddGddS(6,6) = 4*hillParams(4)
      ddGddS(1,2) = -1.0D0*2*hillParams(3)
      ddGddS(1,3) = -1.0D0*2*hillParams(2)
      ddGddS(2,1) = -1.0D0*2*hillParams(3)
      ddGddS(2,3) = -1.0D0*2*hillParams(1)
      ddGddS(3,1) = -1.0D0*2*hillParams(2)
      ddGddS(3,2) = -1.0D0*2*hillParams(1)
      DO i=1,6
        DO j=1,6
          ddGddS(i,j) = multiplier*ddGddS(i,j) - 
     &     dGdS(i)*dGdS(j)/eqGStress
        END DO
      END DO
      RETURN
      END SUBROUTINE calc_ddGddS


      ! calculate first derivatives of yield function
      SUBROUTINE calc_dfdS(hillParams,STRESS,dfdS)
      IMPLICIT NONE
      DOUBLE PRECISION hillParams(4),STRESS(6),dfdS(6),eqStress,
     & calc_eqStress,multiplier
      eqStress = calc_eqStress(hillParams,STRESS)
      multiplier = 0.75/(SUM(hillParams(1:3))*eqStress)
      dfdS(1) = -1.0D0*2*hillParams(2)*(STRESS(3)-STRESS(1))+
     & 2*hillParams(3)*(STRESS(1)-STRESS(2))
      dfdS(2) = 2*hillParams(1)*(STRESS(2)-STRESS(3))-
     & 2*hillParams(3)*(STRESS(1)-STRESS(2))
      dfdS(3) = -1.0D0*2*hillParams(1)*(STRESS(2)-STRESS(3))+
     & 2*hillParams(2)*(STRESS(3)-STRESS(1))
      dfdS(4) = 4*hillParams(4)*STRESS(4)
      dfdS(5) = 4*hillParams(4)*STRESS(5)
      dfdS(6) = 4*hillParams(4)*STRESS(6)
      dfdS(:) = dfdS(:) * multiplier
      RETURN
      END SUBROUTINE calc_dfdS


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
