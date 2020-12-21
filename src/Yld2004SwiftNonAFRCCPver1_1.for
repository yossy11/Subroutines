      ! This is UMAT subroutine made by yossy11, using YLD2004-18p, Swift, non-AFR
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
      CHARACTER(5) orientations(6)
      INTEGER i,iterationNum,YLDM
      DOUBLE PRECISION TOLER,YOUNG,POISSON,HARDK,HARDN,HARDSTRAIN0,
     & lame,shearMod,invDDSDDE(6,6),eStrain(6),pStrain(6),eqpStrain,
     & totalStrain(6),yldCPrime(6,6),yldCDbPrime(6,6),hillParams(4),
     & eqStress,flowStress,calc_eqStress,calc_FlowStress,
     & dfdS(6),calc_dPhidX,dGdS(6),ddGddS(6,6),eqGStress,calc_eqGStress,
     & dLambda,invA(7,7),A(7,7),dpStrain(6),updatedSS(7),x(7),y(7),
     & dFS(6),dFeqpStrain,dF(7),ddLambda,numeratordLambda,
     & denominatordLambda,dyadMat(6,6),F,H,prevF

      ! define constants
      PARAMETER(TOLER=1.0D-2,YOUNG=6.9D4,POISSON=0.33D0,HARDK=646.0D0,
     & HARDN=0.227D0,HARDSTRAIN0=2.5D-2,YLDM=6)

      ! anisotropic params
      yldCPrime(:,:) = 0.0D0
      yldCDbPrime(:,:) = 0.0D0

      yldCPrime(1,2) = -0.0698D0
      yldCPrime(1,3) = 0.9364D0
      yldCPrime(2,1) = 0.0791D0
      yldCPrime(2,3) = 1.0030D0
      yldCPrime(3,1) = 0.5247D0
      yldCPrime(3,2) = 1.3631D0
      yldCPrime(4,4) = 0.9543D0
      yldCPrime(5,5) = 1.0690D0
      yldCPrime(6,6) = 1.0237D0

      yldCDbPrime(1,2) = 0.9811D0
      yldCDbPrime(1,3) = 0.4767D0
      yldCDbPrime(2,1) = 0.5753D0
      yldCDbPrime(2,3) = 0.8668D0
      yldCDbPrime(3,1) = 1.1450D0
      yldCDbPrime(3,2) = -0.0792D0
      yldCDbPrime(4,4) = 1.4046D0
      yldCDbPrime(5,5) = 1.1471D0
      yldCDbPrime(6,6) = 1.0516D0

      hillParams(1) = 0.25216953733566727
      hillParams(2) = 0.8254230293025175
      hillParams(3) = 0.17457697069748246
      hillParams(4) = 2.2380520016508463

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
      CALL calc_Inverse(6,DDSDDE,invDDSDDE)

      ! read STATEV
      CALL ROTSIG(STATEV(1), DROT, eStrain, 2, NDI, NSHR)
      CALL ROTSIG(STATEV(NTENS+1), DROT, pStrain, 2, NDI, NSHR)
      eqpStrain = STATEV(2*NTENS+1)
      eStrain(:) = eStrain(:) + DSTRAN(:)
      totalStrain(:) = eStrain(:) + pStrain(:)


      ! calculate trial stress
      STRESS(:) = STRESS(:) + MATMUL(DDSDDE,DSTRAN)
      iterationNum = 0


      prevF = 10000
      
      DO WHILE (iterationNum < 1000000)
        WRITE(7,*) "here", iterationNum
        ! calculate eqStress and flowStress
        eqStress = calc_eqStress(YLDM,yldCPrime,yldCDbPrime,STRESS)
        flowStress = calc_FlowStress(HARDK,HARDSTRAIN0,HARDN,eqpStrain)
        F = eqStress - flowStress
        ! WRITE(7,*) "eqStress", eqStress
        ! WRITE(7,*) "flowStress", flowStress
        WRITE(7,*) "F", F

        ! if not yield
        IF (F < flowStress*5*TOLER) THEN
          CALL updateSTATEV(NSTATV,STATEV,eStrain,pStrain,eqpStrain)
          RETURN
        END IF

        ! IF (F > prevF) THEN
        !   CALL updateSTATEV(NSTATV,STATEV,eStrain,pStrain,eqpStrain)
        !   RETURN
        ! END IF

        prevF = F

        ! if yield
        eqGStress = calc_eqGStress(hillParams,STRESS)
        H = HARDK*HARDN*((HARDSTRAIN0 + eqpStrain)**(HARDN-1.0D0))
        CALL calc_dfdS(YLDM,yldCPrime,yldCDbPrime,STRESS,dfdS)
        CALL calc_dGdS(hillParams,STRESS,dGdS)
        ! WRITE(7,*) "H", H
        ! WRITE(7,*) "denominator", DOT_PRODUCT(dfdS,MATMUL(DDSDDE,dGdS))
        dLambda = -1.0D0*F/(-1.0D0*DOT_PRODUCT(dfdS,MATMUL(DDSDDE,dGdS))
     &    + H*eqGStress/eqStress)
        ! WRITE(7,*) "eqGStress", eqGStress
        ! WRITE(7,*) "dfdS", dfdS
        ! WRITE(7,*) "dGdS", dGdS
        ! WRITE(7,*) "dLambda", dLambda
        ! dLambda = 1.0D-4

        STRESS(:) = STRESS(:) - dLambda*MATMUL(DDSDDE,dGdS)
        eqpStrain = eqpStrain + dLambda*eqGStress/eqStress
        WRITE(7,*) "eqpStrain",eqpStrain
        pStrain(:) = pStrain(:) + dLambda*dGdS(:)
        eStrain(:) = eStrain(:) - dLambda*dGdS(:)
        iterationNum = iterationNum + 1
      END DO
      WRITE(7,*) "not converged!!!"
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


      ! calculate equivalent stress
      DOUBLE PRECISION FUNCTION calc_eqStress(YLDM,yldCPrime,
     & yldCDbPrime,STRESS)
      INCLUDE 'ABA_PARAM.INC'
      INTEGER YLDM,i,j
      DOUBLE PRECISION yldCPrime(6,6),yldCDbPrime(6,6),STRESS(6),
     & yldT(6,6),yldSPrime(6),yldSDbPrime(6),yldSPriPrime(3),
     & yldSPriDbPrime(3),yldPhi,invariants_(3)
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
      DOUBLE PRECISION hillParams(4),STRESS(6)
      calc_eqGStress = SQRT((1.5D0*(hillParams(1)*(STRESS(2) - 
     & STRESS(3))**2 + hillParams(2)*(STRESS(3)-STRESS(1))**2 + 
     & hillParams(3)*(STRESS(1)-STRESS(2))**2 + 2.0D0*hillParams(4)*
     & SUM(STRESS(4:6)**2)))/(SUM(hillParams(1:3))))
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


      ! calculate differential of plastic potential
      SUBROUTINE calc_dGdS(hillParams,STRESS,dGdS)
      IMPLICIT NONE
      DOUBLE PRECISION hillParams(4),STRESS(6),dGdS(6),eqGStress,
     & calc_eqGStress,multiplier
      eqGStress = calc_eqGStress(hillParams,STRESS)
      multiplier = 0.75/(SUM(hillParams(1:3))*eqGStress)
      dGdS(1) = 2*hillParams(2)*(STRESS(3)-STRESS(1))+2*hillParams(3)*
     & (STRESS(1)-STRESS(2))
      dGdS(2) = 2*hillParams(1)*(STRESS(2)-STRESS(3))+2*hillParams(3)*
     & (STRESS(1)-STRESS(2))
      dGdS(3) = 2*hillParams(1)*(STRESS(2)-STRESS(3))+2*hillParams(2)*
     & (STRESS(3)-STRESS(1))
      dGdS(4) = 4*hillParams(4)*STRESS(4)
      dGdS(5) = 4*hillParams(4)*STRESS(5)
      dGdS(6) = 4*hillParams(4)*STRESS(6)
      dGdS(:) = dGdS(:) * multiplier
      RETURN
      END SUBROUTINE calc_dGdS


      SUBROUTINE calc_dfdS(YLDM,yldCPrime,yldCDbPrime,STRESS,dfdS)
      INTEGER YLDM,i,j
      DOUBLE PRECISION DELTAX,yldCPrime(6,6),yldCDbPrime(6,6),STRESS(6),
     & dfdS(6),eqStress,calc_eqStress,sign,yldSPriPrime(3),
     & yldSPriDbPrime(3),dfdPhi,dPhidS(6),signPrime,signDbPrime,
     & dSds(6,6),subStress(6),yldSSubPriPrime(3),yldSSubPriDbPrime(3)
      PARAMETER(DELTAX=1.0D-6)
      
      eqStress = calc_eqStress(YLDM,yldCPrime,yldCDbPrime,STRESS)
      CALL calc_Principal(yldCPrime,yldCDbPrime,yldSPriPrime,
     & yldSPriDbPrime,STRESS)

      ! dfdPhi = (eqStress/4.0D0)**((1.0D0/YLDM)-1.0D0)/(YLDM*4.0D0)
      dfdPhi = (eqStress**(1-YLDM))/(4.0D0*YLDM)

      
      ! calculate dPhidS
      dPhidS(:) = 0.0D0
      DO i=1,3
        DO j=1,3
          IF (yldSPriPrime(i)- yldSPriDbPrime(j)>=0) THEN
            signPrime = 1.0D0
          ELSE
            signPrime = -1.0D0
          END IF
          IF (yldSPriPrime(j)- yldSPriDbPrime(i)>=0) THEN
            signDbPrime = -1.0D0
          ELSE
            signDbPrime = 1.0D0
          END IF
          dPhidS(i) = dPhidS(i) + signPrime*YLDM*ABS(yldSPriPrime(i) - 
     &     yldSPriDbPrime(j))**(YLDM-1)
          dPhidS(i+3) = dPhidS(i+3) + signDbPrime*YLDM*
     &     ABS(yldSPriPrime(j)- yldSPriDbPrime(i))**(YLDM-1)
        END DO
      END DO

      ! calculate dSds
      dSds(:,:) = 0.0D0
      DO i=1,6
        subStress(:) = STRESS(:)
        subStress(i) = subStress(i)*(1.0D0+DELTAX)
        CALL calc_Principal(yldCPrime,yldCDbPrime,yldSSubPriPrime,
     &   yldSSubPriDbPrime,subStress)
        dSds(1:3,i) = (yldSSubPriPrime(:) - yldSPriPrime(:))/DELTAX
        dSds(4:6,i) = (yldSSubPriDbPrime(:) - yldSPriDbPrime(:))/DELTAX
      END DO
      dfdS(:) = dfdPhi*MATMUL(dPhidS,dSds)
      RETURN
      END SUBROUTINE calc_dfdS


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




        