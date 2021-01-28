      ! This is UMAT subroutine made by yossy11
      ! yield function : Hill48
      ! plastic potential : YLD2004-18p
      ! hardening rule : Kim-Tuan
      ! flow rule : non-AFR
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
      INTEGER i,iterationNum,YLDM
      DOUBLE PRECISION TOLER,YOUNG,POISSON,HARDK,HARDT,HARDA,HARDH,
     & HARDSTRESS0,lame,shearMod,eStrain(6),pStrain(6),eqpStrain,
     & totalStrain(6),yldCPrime(6,6),yldCDbPrime(6,6),hillParams(4),
     & eqStress,flowStress,calc_eqStress,calc_FlowStress,dfdS(6),
     & dGdS(6),eqGStress,calc_eqGStress,dLambda,F,H,calc_H,lambda

      ! define constants
      PARAMETER(TOLER=1.0D-5,YOUNG=22.58D4,POISSON=0.39D0,
     & HARDK=944.03D0,HARDT=23.58D0,HARDA=0.534D0,HARDH=0.185D0,
     & HARDSTRESS0=519.73D0,YLDM=6.0D0)
 
      ! anisotropic params
      yldCPrime(:,:) = 0.0D0
      yldCDbPrime(:,:) = 0.0D0

      yldCPrime(1,2) = 0.9630755236085117D0
      yldCPrime(1,3) = 1.0404635133877267D0
      yldCPrime(2,1) = 0.9671641135689687D0
      yldCPrime(2,3) = 0.9666250667631737D0
      yldCPrime(3,1) = 1.0364152750397875D0
      yldCPrime(3,2) = 0.9707464664422254D0
      yldCPrime(4,4) = 1.0278975580349605D0
      yldCPrime(5,5) = 1.0D0
      yldCPrime(6,6) = 1.0D0

      yldCDbPrime(1,2) = 0.9630621276929356D0
      yldCDbPrime(1,3) = 1.0405349522836391D0
      yldCDbPrime(2,1) = 0.9671784592047932D0
      yldCDbPrime(2,3) = 0.9666518844943464D0
      yldCDbPrime(3,1) = 1.0363759388697624D0
      yldCDbPrime(3,2) = 0.9707396235439515D0
      yldCDbPrime(4,4) = 1.027884484472155D0
      yldCDbPrime(5,5) = 1.0D0
      yldCDbPrime(6,6) = 1.0D0

      hillParams(1) = 0.4844687812695614D0
      hillParams(2) = 0.5441303900909464D0
      hillParams(3) = 0.4558696099090537D0
      hillParams(4) = 1.5201303973492752D0

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
        eqStress = calc_eqStress(hillParams,STRESS)
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
            CALL updateDDSDDE(YLDM,yldCPrime,yldCDbPrime,hillParams,
     &       STRESS,DDSDDE,lambda,eqStress,HARDK,HARDT,HARDA,HARDH,
     &       HARDSTRESS0,eqpStrain)
          END IF
          RETURN
        END IF

        ! if yield
        eqGStress = calc_eqGStress(YLDM,yldCPrime,yldCDbPrime,STRESS)
        H = calc_H(HARDK,HARDT,HARDA,HARDH,HARDSTRESS0,eqpStrain)
        CALL calc_dfdS(hillParams,STRESS,dfdS)
        CALL calc_dGdS(YLDM,yldCPrime,yldCDbPrime,STRESS,dGdS)

        dLambda = F/(DOT_PRODUCT(dfdS,MATMUL(DDSDDE,dGdS)) + 
     &   H*eqGStress/eqStress)
        IF (dLambda<0) THEN
          WRITE(7,*) "dLambda < 0, invalid!!"
          CALL XIT
        END IF

        lambda = lambda + dLambda

        STRESS(:) = STRESS(:) - dLambda*MATMUL(DDSDDE,dGdS)
        eqpStrain = eqpStrain + dLambda*eqGStress/eqStress
        pStrain(:) = pStrain(:) + dLambda*dGdS(:)
        eStrain(:) = eStrain(:) - dLambda*dGdS(:)
        iterationNum = iterationNum + 1
      END DO
      WRITE(7,*) "not converged!!!"
      WRITE(7,*) "F",F
      CALL XIT
      RETURN
      END SUBROUTINE UMAT


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


      ! calculate equivalent stress
      DOUBLE PRECISION FUNCTION calc_eqGStress(YLDM,yldCPrime,
     & yldCDbPrime,STRESS)
      IMPLICIT NONE
      INTEGER YLDM,i,j
      DOUBLE PRECISION yldCPrime(6,6),yldCDbPrime(6,6),STRESS(6),
     & yldSPriPrime(3),yldSPriDbPrime(3),yldPhi
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
      calc_eqGStress = (yldPhi/4.0D0)**(1.0D0/YLDM)
      RETURN
      END FUNCTION calc_eqGStress


      ! calculate flow stress
      DOUBLE PRECISION FUNCTION calc_FlowStress(HARDK,HARDT,HARDA,HARDH,
     & HARDSTRESS0,eqpStrain)
      IMPLICIT NONE
      DOUBLE PRECISION HARDK,HARDT,HARDA,HARDH,HARDSTRESS0,eqpStrain
      IF (eqpStrain == 0.0D0) THEN
        calc_FlowStress = HARDSTRESS0
        RETURN
      END IF
      calc_FlowStress = HARDSTRESS0 + 
     & HARDK*(1.0D0-EXP(-1.0D0*HARDT*eqpStrain**HARDA))*eqpStrain**HARDH
      RETURN
      END FUNCTION calc_FlowStress


      ! calculate derivative of flow stress
      DOUBLE PRECISION FUNCTION calc_H(HARDK,HARDT,HARDA,HARDH,
     & HARDSTRESS0,eqpStrain)
      IMPLICIT NONE
      DOUBLE PRECISION HARDK,HARDT,HARDA,HARDH,HARDSTRESS0,eqpStrain
      IF (eqpStrain == 0.0D0) THEN
        calc_H = 0.0D0
        RETURN
      END IF
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
      SUBROUTINE updateDDSDDE(YLDM,yldCPrime,yldCDbPrime,hillParams,
     & STRESS,DDSDDE,lambda,eqStress,HARDK,HARDT,HARDA,HARDH,
     & HARDSTRESS0,eqpStrain)
      IMPLICIT NONE
      INTEGER YLDM,i,j
      DOUBLE PRECISION yldCPrime(6,6),yldCDbPrime(6,6),hillParams(4),
     & STRESS(6),DDSDDE(6,6),lambda,eqStress,HARDK,HARDT,HARDA,HARDH,
     & HARDSTRESS0,eqpStrain,eqGStress,calc_eqGStress,dfdS(6),dGdS(6),
     & ddGddS(6,6),A(7,7),B(6,6),invB(6,6),invDDSDDE(6,6),h0,subVec(7),
     & vec1(7),vec2(7),H,calc_H,C(7,7),subDDSDDE(6,6)
      subDDSDDE(:,:) = DDSDDE(:,:)
      eqGStress = calc_eqGStress(YLDM,yldCPrime,yldCDbPrime,STRESS)
      CALL calc_dfdS(hillParams,STRESS,dfdS)
      CALL calc_dGdS(YLDM,yldCPrime,yldCDbPrime,STRESS,dGdS)
      CALL calc_ddGddS(YLDM,yldCPrime,yldCDbPrime,STRESS,dGdS,ddGddS)
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
      H = calc_H(HARDK,HARDT,HARDA,HARDH,HARDSTRESS0,eqpStrain)
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


      ! calculate second derivatives of plastic potential
      SUBROUTINE calc_ddGddS(YLDM,yldCPrime,yldCDbPrime,STRESS,dGdS,
     & ddGddS)
      IMPLICIT NONE
      INTEGER YLDM,i
      DOUBLE PRECISION yldCPrime(6,6),yldCDbPrime(6,6),STRESS(6),
     & dGdS(6),ddGddS(6,6),DELTAX,subStress(6),subdGdS(6)
      PARAMETER(DELTAX=1.0D-1)
      DO i=1,6
        subStress(:) = STRESS(:)
        subStress(i) = subStress(i) + DELTAX
        CALL calc_dGdS(YLDM,yldCPrime,yldCDbPrime,STRESS,subdGdS)
        ddGddS(:,i) = (subdGdS(:) - dGdS(:))/DELTAX
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


      ! calculate first derivatives of plastic potential
      SUBROUTINE calc_dGdS(YLDM,yldCPrime,yldCDbPrime,STRESS,dGdS)
      IMPLICIT NONE
      INTEGER YLDM,i,j,k
      DOUBLE PRECISION yldCPrime(6,6),yldCDbPrime(6,6),STRESS(6),
     & dGdS(6),DELTAX,yldSPriPrime(3),yldSPriDbPrime(3),yldPhi,
     & subStress(6),subyldPhi,dPhidS(6),eqGStress,calc_eqGStress,dGdPhi
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
      eqGStress = calc_eqGStress(YLDM,yldCPrime,yldCDbPrime,STRESS)
      dGdPhi = (eqGStress**(1-YLDM))/(4.0D0*YLDM)
      dGdS(:) = dGdPhi*dPhidS(:)
      RETURN
      END SUBROUTINE calc_dGdS


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
