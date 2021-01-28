      ! This is UMAT subroutine made by yossy11
      ! yield function : YLD2004-18p
      ! plastic potential : Hill48
      ! hardening rule : Kim-Tuan
      ! flow rule : AFR
      ! integration algorithm : CCP(convex cutting plane)
      ! using numerical differentiation for the calculation of dfdS

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
      DOUBLE PRECISION TOLER,YOUNG,POISSON,HARDK,HARDT,HARDA,HARDH,
      HARDSTRESS0,lame,
     & shearMod,eStrain(6),pStrain(6),eqpStrain,totalStrain(6),
     & yldCPrime(6,6),yldCDbPrime(6,6),hillParams(4),eqStress,
     & flowStress,calc_eqStress,calc_FlowStress,dfdS(6),dGdS(6),
     & eqGStress,calc_eqGStress,dLambda,F,H,calc_H,lambda,YLDM

      ! define constants
      PARAMETER(TOLER=1.0D-5,YOUNG=22.58D4,POISSON=0.39D0,
     & HARDK=944.03D0,HARDT=23.58D0,HARDA=0.534D0,HARDH=0.185D0,
     & HARDSTRESS0=519.73D0,YLDM=6.0D0)

      ! anisotropic params
      yldCPrime(:,:) = 0.0D0
      yldCDbPrime(:,:) = 0.0D0

      yldCPrime(1,2) = -0.978676811742272D0
      yldCPrime(1,3) = -1.0155655137865276D0
      yldCPrime(2,1) = -0.956293501571447D0
      yldCPrime(2,3) = -0.9923044145704766D0
      yldCPrime(3,1) = -1.0379668696483502D0
      yldCPrime(3,2) = -0.9699665765886128D0
      yldCPrime(4,4) = 1.0229368224938618D0
      yldCPrime(5,5) = 1.0D0
      yldCPrime(6,6) = 1.0D0

      yldCDbPrime(1,2) = -0.978658158833893D0
      yldCDbPrime(1,3) = -1.0156280173593986D0
      yldCDbPrime(2,1) = -0.9563282984646336D0
      yldCDbPrime(2,3) = -0.992335717702875D0
      yldCDbPrime(3,1) = -1.037930467713471D0
      yldCDbPrime(3,2) = -0.9699600496984562D0
      yldCDbPrime(4,4) = 1.0229238225580262D0
      yldCDbPrime(5,5) = 1.0D0
      yldCDbPrime(6,6) = 1.0D0

      hillParams(1) = 0.36784658988120367
      hillParams(2) = 0.6211180124223603
      hillParams(3) = 0.37888198757763975
      hillParams(4) = 1.0581921244648134

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
        eqStress = calc_eqStress(YLDM,yldCPrime,yldCDbPrime,STRESS)
        flowStress = calc_FlowStress(HARDK,HARDT,HARDA,HARDH,
     &   HARDSTRESS0,eqpStrain)
        F = eqStress - flowStress
        IF (ISNAN(F)) THEN
          WRITE(7,*) "F is NaN"
          CALL XIT
        END IF

        ! if not yield
        IF (F < flowStress*TOLER) THEN
          CALL updateSTATEV(NSTATV,STATEV,eStrain,pStrain,eqpStrain)
          IF (iterationNum /= 0) THEN
            CALL updateDDSDDE(hillParams,YLDM,yldCPrime,yldCDbPrime,
     &       STRESS,DDSDDE,lambda,eqStress,HARDK,HARDT,HARDA,HARDH,
     &       HARDSTRESS0,eqpStrain)
          END IF
          RETURN
        END IF

        ! if yield
        eqGStress = calc_eqGStress(hillParams,STRESS)
        H = calc_H(HARDK,HARDT,HARDA,HARDH,HARDSTRESS0,eqpStrain)
        CALL calc_dfdS(YLDM,yldCPrime,yldCDbPrime,STRESS,dfdS)
        CALL calc_dGdS(hillParams,STRESS,dGdS)

    !     dLambda = F/(DOT_PRODUCT(dfdS,MATMUL(DDSDDE,dGdS)) + 
    !  &   H*eqGStress/eqStress)
        dLambda = F/(DOT_PRODUCT(dfdS,MATMUL(DDSDDE,dfdS)) + H)
        IF (dLambda<0) THEN
          WRITE(7,*) "dLambda < 0, invalid!!"
          CALL XIT
        END IF

        lambda = lambda + dLambda

        ! STRESS(:) = STRESS(:) - dLambda*MATMUL(DDSDDE,dGdS)
        ! eqpStrain = eqpStrain + dLambda*eqGStress/eqStress
        ! pStrain(:) = pStrain(:) + dLambda*dGdS(:)
        ! eStrain(:) = eStrain(:) - dLambda*dGdS(:)
        STRESS(:) = STRESS(:) - dLambda*MATMUL(DDSDDE,dfdS)
        eqpStrain = eqpStrain + dLambda
        pStrain(:) = pStrain(:) + dLambda*dfdS(:)
        eStrain(:) = eStrain(:) - dLambda*dfdS(:)
        iterationNum = iterationNum + 1
      END DO
      WRITE(7,*) "not converged!!!"
      WRITE(7,*) "F",F
      CALL XIT
      RETURN
      END SUBROUTINE UMAT


      ! calculate equivalent stress
      DOUBLE PRECISION FUNCTION calc_eqStress(YLDM,yldCPrime,
     & yldCDbPrime,STRESS)
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION yldCPrime(6,6),yldCDbPrime(6,6),STRESS(6),
     & yldSPriPrime(3),yldSPriDbPrime(3),yldPhi,YLDM
      yldSPriPrime(:) = 0.0D0
      yldSPriDbPrime(:) = 0.0D0
      CALL calc_Principal(yldCPrime,yldCDbPrime,yldSPriPrime,
     & yldSPriDbPrime,STRESS)
      yldPhi = 0.0D0
      DO i=1,3
        DO j=1,3
          yldPhi = yldPhi + 
     &     ABS(yldSPriPrime(i) - yldSPriDbPrime(j))**YLDM
        END DO
      END DO
      calc_eqStress = (yldPhi/4.0D0)**(1.0D0/YLDM)
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
      DOUBLE PRECISION FUNCTION calc_FlowStress(HARDK,HARDT,HARDA,HARDH,
     & HARDSTRESS0,eqpStrain)
      IMPLICIT NONE
      DOUBLE PRECISION HARDK,HARDT,HARDA,HARDH,HARDSTRESS0,eqpStrain
      calc_FlowStress = HARDSTRESS0 + 
     & HARDK*(1.0D0-EXP(-1.0D0*HARDT*eqpStrain**HARDA))*eqpStrain**HARDH
      RETURN
      END FUNCTION calc_FlowStress


      ! calculate derivative of flow stress
      DOUBLE PRECISION FUNCTION calc_H(HARDK,HARDT,HARDA,HARDH,
     & HARDSTRESS0,eqpStrain)
      IMPLICIT NONE
      DOUBLE PRECISION HARDK,HARDT,HARDA,HARDH,HARDSTRESS0,eqpStrain
      calc_H = (HARDK*HARDA*HARDT*eqpStrain**(HARDA-1.0D0))*
     & EXP(-1.0D0*HARDT*eqpStrain**HARDA)
      calc_H = calc_H + 
     & HARDK*(1.0D0-EXP(-1.0D0*HARDT*eqpStrain**HARDA))*
     & HARDH*eqpStrain**(HARDH-1.0D0)
      RETURN
      END FUNCTION calc_H


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
      SUBROUTINE updateDDSDDE(hillParams,YLDM,yldCPrime,yldCDbPrime,
     & STRESS,DDSDDE,lambda,eqStress,HARDK,HARDT,HARDA,HARDH,
     & HARDSTRESS0,eqpStrain)
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION hillParams(4),yldCPrime(6,6),yldCDbPrime(6,6),
     & STRESS(6),DDSDDE(6,6),lambda,eqStress,HARDK,HARDT,HARDA,HARDH,
     & HARDSTRESS0,eqpStrain,eqGStress,calc_eqGStress,dfdS(6),dGdS(6),
     & ddGddS(6,6),A(7,7),B(6,6),invB(6,6),invDDSDDE(6,6),h0,subVec(7),
     & vec1(7),vec2(7),H,calc_H,C(7,7),subDDSDDE(6,6),YLDM,ddfddS(6,6)
      subDDSDDE(:,:) = DDSDDE(:,:)
      eqGStress = calc_eqGStress(hillParams,STRESS)
      CALL calc_dfdS(YLDM,yldCPrime,yldCDbPrime,STRESS,dfdS)
      CALL calc_ddfddS(YLDM,yldCPrime,yldCDbPrime,STRESS,dfdS,ddfddS)
      CALL calc_dGdS(hillParams,STRESS,dGdS)
      CALL calc_ddGddS(hillParams,STRESS,dGdS,ddGddS)
      CALL calc_Inverse(6,DDSDDE,invDDSDDE)
      ! invB(:,:) = invDDSDDE(:,:) + lambda*ddGddS(:,:)
      invB(:,:) = invDDSDDE(:,:) + lambda*ddfddS(:,:)
      CALL calc_Inverse(6,invB,B)
      A(:,:) = 0.0D0
      A(1:6,1:6) = B(:,:)
      A(7,7) = 1.0D0
      ! h0 = eqGStress/eqStress
      h0 = 1.0D0
      subVec(:) = 0.0D0
      ! subVec(1:6) = dGdS(:)
      subVec(1:6) = dfdS(:)
      subVec(7) = -1.0D0*h0
      vec1 = MATMUL(A,subVec)
      subVec(1:6) = dfdS(:)
      subVec(7) = 0.0D0
      vec2 = MATMUL(subVec,A)
      CALL calc_Dyad(7,vec1,vec2,C)
      H = calc_H(HARDK,HARDT,HARDA,HARDH,HARDSTRESS0,eqpStrain)
    !   C(:,:) = C(:,:)/
    !  & (DOT_PRODUCT(dfdS,MATMUL(B,dGdS)) + H*h0)
      C(:,:) = C(:,:)/
     & (DOT_PRODUCT(dfdS,MATMUL(B,dfdS)) + H*h0)
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


      ! calculate differential of plastic potential
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


      ! calculate differential of plastic potential
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


      SUBROUTINE calc_dfdS(YLDM,yldCPrime,yldCDbPrime,STRESS,dfdS)
      IMPLICIT NONE
      INTEGER i,j,k
      DOUBLE PRECISION yldCPrime(6,6),yldCDbPrime(6,6),STRESS(6),
     & dfdS(6),DELTAX,yldSPriPrime(3),yldSPriDbPrime(3),yldPhi,YLDM,
     & subStress(6),subyldPhi,dPhidS(6),eqStress,calc_eqStress,dfdPhi
      PARAMETER(DELTAX=1.0D-5)
      yldSPriPrime(:) = 0.0D0
      yldSPriDbPrime(:) = 0.0D0
      CALL calc_Principal(yldCPrime,yldCDbPrime,yldSPriPrime,
     & yldSPriDbPrime,STRESS)
      yldPhi = 0.0D0
      DO i=1,3
        DO j=1,3
          yldPhi = yldPhi + 
     &     ABS(yldSPriPrime(i) - yldSPriDbPrime(j))**YLDM
        END DO
      END DO
      DO i=1,6
        subStress(:) = STRESS(:)
        subStress(i) = subStress(i) + DELTAX
        yldSPriPrime(:) = 0.0D0
        yldSPriDbPrime(:) = 0.0D0
        CALL calc_Principal(yldCPrime,yldCDbPrime,yldSPriPrime,
     &   yldSPriDbPrime,subStress)
        subyldPhi = 0.0D0
        DO j=1,3
          DO k=1,3
            subyldPhi = subyldPhi + 
     &       ABS(yldSPriPrime(j) - yldSPriDbPrime(k))**YLDM
          END DO
        END DO
        dPhidS(i) = (subyldPhi - yldPhi)/DELTAX
      END DO
      eqStress = calc_eqStress(YLDM,yldCPrime,yldCDbPrime,STRESS)
      dfdPhi = (eqStress**(1-YLDM))/(4.0D0*YLDM)
      dfdS(:) = dfdPhi*dPhidS(:)
      RETURN
      END SUBROUTINE calc_dfdS


      SUBROUTINE calc_ddfddS(YLDM,yldCPrime,yldCDbPrime,STRESS,dfdS,
     & ddfddS)
      IMPLICIT NONE
      INTEGER YLDM,i
      DOUBLE PRECISION yldCPrime(6,6),yldCDbPrime(6,6),STRESS(6),
     & dfdS(6),ddfddS(6,6),DELTAX,subStress(6),subdfdS(6)
      PARAMETER(DELTAX=1.0D-1)
      DO i=1,6
        subStress(:) = STRESS(:)
        subStress(i) = subStress(i) + DELTAX
        CALL calc_dfdS(YLDM,yldCPrime,yldCDbPrime,STRESS,subdfdS)
        ddfddS(:,i) = (subdfdS(:) - dfdS(:))/DELTAX
      END DO  
      RETURN
      END SUBROUTINE calc_ddfddS
      

      SUBROUTINE calc_Principal(yldCPrime,yldCDbPrime,yldSPriPrime,
     &   yldSPriDbPrime,STRESS)
      INCLUDE 'ABA_PARAM.INC'
      INTEGER i
      DOUBLE PRECISION yldCPrime(6,6),yldCDbPrime(6,6),yldSPriPrime(3),
     & yldSPriDbPrime(3),STRESS(6),yldT(6,6),yldSPrime(6),yldSDbPrime(6)
      yldT(:,:) = 0.0D0
      yldT(1:3,1:3) = -1.0D0
      DO i=1,3
        yldT(i,i) = 2.0D0
        yldT(i+3,i+3) = 3.0D0
      END DO
      yldT(:,:) = yldT(:,:)/3.0D0
      yldSPrime = MATMUL(yldCPrime,MATMUL(yldT,STRESS))
      yldSDbPrime = MATMUL(yldCDbPrime,MATMUL(yldT,STRESS))
      yldSPriPrime(:) = 0.0D0
      yldSPriDbPrime(:) = 0.0D0
      CALL SPRINC(yldSPrime,yldSPriPrime,1,3,3)
      CALL SPRINC(yldSDbPrime,yldSPriDbPrime,1,3,3)
      RETURN
      END SUBROUTINE calc_Principal


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
