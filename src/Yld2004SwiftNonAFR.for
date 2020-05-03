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
      INTEGER i
      DOUBLE PRECISION TOLER,YOUNG,POISSON,HARDK,HARDN,HARDSTRAIN0,YLDM,
     & lame,shearMod,invDDSDDE(6,6),eStrain(6),pStrain(6),eqpStrain,
     & totalStrain(6),yldCPrime(6,6),yldCDbPrime(6,6),hillParams(4),
     & eqStress,effectiveStress,calc_eqStress,calc_EffectiveStress,
     & dfdS(6),calc_dPhidX,dGdS(6),ddGddS(6,6),eqGStress,calc_eqGStress,
     & dLambda,invA(7,7),A(7,7),dpStrain(6),updatedSS(7),x(7),y(7),
     & dFS(6),dFeqpStrain,dF(7),ddLambda,numeratordLambda,
     & denominatordLambda,dyadMat(6,6)

      ! define constants
      PARAMETER(TOLER=1.0D-6,YOUNG=7.03D4,POISSON=0.3D0,HARDK=615.3D0,
     & HARDN=0.363D0,HARDSTRAIN0=7.61D-3,YLDM=6)

      ! this subroutine can be used only on the condition of NDI=NSHR=3
      IF (NDI/=3 .or. NSHR/=3) THEN
        WRITE(6,*) 'invalid element type, please check about that'
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
      eqPStrain = STATEV(2*NTENS+1)
      eStrain(:) = eStrain(:) + DSTRAN(:)
      totalStrain(:) = eStrain(:) + pStrain(:)

      ! calculate trial stress
      STRESS(:) = STRESS(:) + MATMUL(DDSDDE,DSTRAN)

      ! anisotropic params
      yldCPrime(:,:) = 0.0D0
      yldCPrime(5,5) = 1.0D0
      yldCPrime(6,6) = 1.0D0
      yldCDbPrime(:,:) = 0.0D0
      yldCDbPrime(5,5) = 1.0D0
      yldCDbPrime(6,6) = 1.0D0

      yldCPrime(1,2) = 1.0D0
      yldCPrime(1,3) = 1.0D0
      yldCPrime(2,1) = 1.0D0
      yldCPrime(2,3) = 1.0D0
      yldCPrime(3,1) = 1.0D0
      yldCPrime(3,2) = 1.0D0
      yldCPrime(4,4) = 1.0D0

      yldCDbPrime(1,2) = 1.0D0
      yldCDbPrime(1,3) = 1.0D0
      yldCDbPrime(2,1) = 1.0D0
      yldCDbPrime(2,3) = 1.0D0
      yldCDbPrime(3,1) = 1.0D0
      yldCDbPrime(3,2) = 1.0D0
      yldCDbPrime(4,4) = 1.0D0

      hillParams(1) = 0.45360550654844645
      hillParams(2) = 0.50524270631786716
      hillParams(3) = 0.49475729368213273
      hillParams(4) = 1.5532035559766002

      ! calculate eqStress and effectiveStress
      eqStress = calc_eqStress(YLDM,yldCPrime,yldCDbPrime,STRESS)
      effectiveStress = 
     & calc_EffectiveStress(HARDK,HARDSTRAIN0,HARDN,eqpStrain)

      ! if not yield
      IF ((eqStress - effectiveStress) < effectiveStress*TOLER) THEN
        CALL updateSTATEV(NSTATV,STATEV,eStrain,pStrain,eqpStrain)
        RETURN
      END IF
      
      ! if yield
      orientations = ['sigxx','sigyy','sigzz','sigxy','sigxz','sigyz']
      DO i=1,6
        dfdS(i) = 
     &   calc_dPhidX(orientations(i),YLDM,yldCPrime,yldCDbPrime,STRESS)
      END DO
      CALL calc_dGdS(hillParams,STRESS,dGdS)
      CALL calc_ddGddS(hillParams,STRESS,dGdS,ddGddS)
      eqGStress = calc_eqGStress(hillParams,STRESS)
      dLambda = DOT_PRODUCT(dfdS,MATMUL(DDSDDE,DSTRAN))/
     & (DOT_PRODUCT(dfdS,MATMUL(DDSDDE,dGdS)) + HARDK*HARDN*
     & ((HARDSTRAIN0 + eqpStrain)**(HARDN-1.0D0))*eqGStress/eqStress)
      
      ! update states to trial step of N-R iteration
      STRESS(:) = STRESS(:) - dLambda*MATMUL(DDSDDE,dGdS)
      CALL calc_dGdS(hillParams,STRESS,dGdS)
      CALL calc_ddGddS(hillParams,STRESS,dGdS,ddGddS)
      eqStress = calc_eqStress(YLDM,yldCPrime,yldCDbPrime,STRESS)
      eqGStress = calc_eqGStress(hillParams,STRESS)
      pStrain(:) = pStrain(:) + dLambda*dGdS(:)
      eqpStrain = eqpStrain + dLambda*eqGStress/eqStress
      invA(:,:) = 0.0D0
      invA(7,7) = -1.0D0
      dpStrain = 0.0D0
      updatedSS(:) = 0.0D0

      ! N-R iteration
      DO WHILE ((eqStress - effectiveStress)>=effectiveStress*TOLER)
        x(1:6) = dLambda*dGdS(:) - dpStrain
        x(7) = dLambda*eqGStress/eqStress - updatedSS(7)
        y(1:6) = dGdS(:)
        y(7) = 1.0D0
        invA(1:6,1:6) = invDDSDDE(:,:) + dLambda*ddGddS(:,:)
        CALL calc_Inverse(6,invA,A)
        DO i=1,6
          dFS(i) = calc_dPhidX(orientations(i),YLDM,
     &     yldCPrime,yldCDbPrime,STRESS)
        END DO
        dFS(:) = (dFS(:)*eqStress**(1.0D0-YLDM))/YLDM
        dFeqpStrain = 
     &   -1.0D0*HARDK*HARDN*(HARDSTRAIN0+eqpStrain)**(HARDN-1.0D0)
        dF(1:6) = dFS(:)
        dF(7) = dFeqpStrain
        ddLambda =(eqStress - effectiveStress - 
     &   DOT_PRODUCT(dF,MATMUL(A,x)))/DOT_PRODUCT(dF,MATMUL(A,y))
        updatedSS = (-1.0D0)*(MATMUL(A,x) + ddLambda*MATMUL(A,y))
        dpStrain = (-1.0D0)*MATMUL(invDDSDDE,updatedSS(1:6))
        pStrain(:) = pStrain(:) + dpStrain(:)
        eqpStrain = eqpStrain + updatedSS(7)
        dLambda = dLambda + ddLambda
        STRESS = STRESS + updatedSS(1:6)
        eStrain(:) = totalStrain(:) - pStrain(:)
        eqStress = calc_eqStress(YLDM,yldCPrime,yldCDbPrime,STRESS)
        eqGStress = calc_eqGStress(hillParams,STRESS)
        effectiveStress = 
     &   calc_EffectiveStress(HARDK,HARDSTRAIN0,HARDN,eqpStrain)
        CALL calc_dGdS(hillParams,STRESS,dGdS)
      END DO

      ! update DDSDDE
      DO i=1,6
        dfdS(i) = 
     &   calc_dPhidX(orientations(i),YLDM,yldCPrime,yldCDbPrime,STRESS)
      END DO
      numeratordLambda = DOT_PRODUCT(dfdS,MATMUL(DDSDDE,DSTRAN))
      denominatordLambda = dLambda/numeratordLambda
      CALL calc_dyad(6,MATMUL(DDSDDE,dGdS),MATMUL(dfdS,DDSDDE),dyadMat)
      DDSDDE(:,:) = DDSDDE(:,:) - dyadMat(:,:)/denominatordLambda
      CALL updateSTATEV(NSTATV,STATEV,eStrain,pStrain,eqpStrain)
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
      IMPLICIT NONE
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
      CALL calc_Principal(yldSPrime,yldSPriPrime,invariants_)
      CALL calc_Principal(yldSDbPrime,yldSPriDbPrime,invariants_)
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
     & hillParams(3)*(STRESS(1)-STRESS(3))**2 + 2.0D0*hillParams(4)*
     & SUM(STRESS(4:6)**2)))/(SUM(hillParams(1:3))))
      RETURN
      END FUNCTION calc_eqGStress


      ! calculate effective stress
      DOUBLE PRECISION FUNCTION calc_EffectiveStress(HARDK,HARDSTRAIN0,
     & HARDN,eqpStrain)
      IMPLICIT NONE
      DOUBLE PRECISION HARDK,HARDSTRAIN0,HARDN,eqpStrain
      calc_EffectiveStress = HARDK*((HARDSTRAIN0 + eqpStrain)**HARDN)
      RETURN
      END FUNCTION calc_EffectiveStress


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


      ! calculate differential dPhi/dX e.g. dPhi/dSzz
      DOUBLE PRECISION FUNCTION calc_dPhidX(orientation,YLDM,yldCPrime,
     & yldCDbPrime,bufStress)
      IMPLICIT NONE
      CHARACTER(5) orientation
      INTEGER YLDM,i
      DOUBLE PRECISION yldCPrime(6,6),yldCDbPrime(6,6),bufStress(6),
     & yldT(6,6),yldSPrime(6),yldSDbPrime(6),yldSPriPrime(3),
     & yldSPriDbPrime(3),invariantsPrime(3),invariantsDbPrime(3),
     & dSTildedXPrime(6),dSTildedXDbPrime(6),dHdSTildePrime(3,6),
     & dHdSTildeDbPrime(3,6),dHdXPrime(3),dHdXDbPrime(3),
     & dSPridHPrime(3,3),dSPridHDbPrime(3,3),dSPridXPrime(3),
     & dSPridXDbPrime(3),dPhidSPriPrime(3),dPhidSPriDbPrime(3),
     & denominator
      yldT(:,:) = 0.0D0
      yldT(1:3,1:3) = -1.0D0
      DO i=1,3
        yldT(i,i) = 2.0D0
        yldT(i+3,i+3) = 3.0D0
      END DO
      yldT(:,:) = yldT(:,:)/3.0D0
      yldSPrime = MATMUL(yldCPrime,MATMUL(yldT,bufStress))
      yldSDbPrime = MATMUL(yldCDbPrime,MATMUL(yldT,bufStress))
      CALL calc_Principal(yldSPrime,yldSPriPrime,invariantsPrime)
      CALL calc_Principal(yldSDbPrime,yldSPriDbPrime,invariantsDbPrime)
      
      dSTildedXPrime(:) = 0.0D0
      dSTildedXDbPrime(:) = 0.0D0

      SELECT CASE (orientation)
      CASE ('sigxx')
        dSTildedXPrime(1) = (yldCPrime(1,2) + yldCPrime(1,3))/3.0D0
        dSTildedXPrime(2) = (-2.0D0*yldCPrime(2,1) + yldCPrime(2,3))/3.0D0
        dSTildedXPrime(3) = (-2.0D0*yldCPrime(3,1) + yldCPrime(3,2))/3.0D0
        dSTildedXDbPrime(1) = (yldCDbPrime(1,2) + yldCDbPrime(1,3))/3.0D0
        dSTildedXDbPrime(2) = 
     &   (-2.0D0*yldCDbPrime(2,1) + yldCDbPrime(2,3))/3.0D0
        dSTildedXDbPrime(3) = 
     &   (-2.0D0*yldCDbPrime(3,1) + yldCDbPrime(3,2))/3.0D0
      CASE ('sigyy')
        dSTildedXPrime(1) = 
     &   (-2.0D0*yldCPrime(1,2) + yldCPrime(1,3))/3.0D0
        dSTildedXPrime(2) = (yldCPrime(2,1) + yldCPrime(2,3))/3.0D0
        dSTildedXPrime(3) = 
     &   (yldCPrime(3,1) - 2.0D0*yldCPrime(3,2))/3.0D0
        dSTildedXDbPrime(1) = 
     &   (-2.0D0*yldCDbPrime(1,2) + yldCDbPrime(1,3))/3.0D0
        dSTildedXDbPrime(2) = 
     &   (yldCDbPrime(2,1) + yldCDbPrime(2,3))/3.0D0
        dSTildedXDbPrime(3) = 
     &   (yldCDbPrime(3,1) - 2.0D0*yldCDbPrime(3,2))/3.0D0
      CASE ('sigzz')
        dSTildedXPrime(1) = 
     &   (yldCPrime(1,2) - 2.0D0*yldCPrime(1,3))/3.0D0
        dSTildedXPrime(2) = 
     &   (yldCPrime(2,1) - 2.0D0*yldCPrime(2,3))/3.0D0
        dSTildedXPrime(3) = (yldCPrime(3,1) + yldCPrime(3,2))/3.0D0
        dSTildedXDbPrime(1) = 
     &   (yldCDbPrime(1,2) - 2.0D0*yldCDbPrime(1,3))/3.0D0
        dSTildedXDbPrime(2) = 
     &   (yldCDbPrime(2,1) - 2.0D0*yldCDbPrime(2,3))/3.0D0
        dSTildedXDbPrime(3) = 
     &   (yldCDbPrime(3,1) + yldCDbPrime(3,2))/3.0D0
      CASE ('sigxy')
        dSTildedXPrime(4) = yldCPrime(4,4)
        dSTildedXDbPrime(4) = yldCDbPrime(4,4)
      CASE ('sigxz')
        dSTildedXPrime(5) = 1.0D0
        dSTildedXDbPrime(5) = 1.0D0
      CASE ('sigyz')
        dSTildedXPrime(6) = 1.0D0
        dSTildedXDbPrime(6) = 1.0D0
      CASE DEFAULT
        WRITE(6,*) 'Invalid orientation'
        CALL XIT
      END SELECT

      dHdSTildePrime(1,:) = 0.0D0
      dHdSTildePrime(1,1:3) = 1.0D0/3.0D0
      dHdSTildeDbPrime(1,:) = 0.0D0
      dHdSTildeDbPrime(1,1:3) = 1.0D0/3.0D0
      
      DO i=1,3
        dHdSTildePrime(2,i) = (yldSPrime(i) - SUM(yldSPrime(1:3)))/3.0D0
        dHdSTildePrime(2,i+3) = 2.0D0*yldSPrime(i+3)/3.0D0
        dHdSTildeDbPrime(2,i) = 
     &   (yldSDbPrime(i) - SUM(yldSDbPrime(1:3)))/3.0D0
        dHdSTildeDbPrime(2,i+3) = 2.0D0*yldSDbPrime(i+3)/3.0D0
      END DO
      dHdSTildePrime(3,1) = 
     & (yldSPrime(2)*yldSPrime(3) - yldSPrime(6)**2)/2.0D0
      dHdSTildePrime(3,2) = 
     & (yldSPrime(3)*yldSPrime(1) - yldSPrime(5)**2)/2.0D0
      dHdSTildePrime(3,3) = 
     & (yldSPrime(1)*yldSPrime(2) - yldSPrime(4)**2)/2.0D0
      dHdSTildePrime(3,4) = 
     & yldSPrime(5)*yldSPrime(6) - yldSPrime(3)*yldSPrime(4)
      dHdSTildePrime(3,5) = 
     & yldSPrime(4)*yldSPrime(6) - yldSPrime(2)*yldSPrime(5)
      dHdSTildePrime(3,6) = 
     & yldSPrime(4)*yldSPrime(5) - yldSPrime(1)*yldSPrime(6)

      dHdSTildeDbPrime(3,1) = 
     & (yldSDbPrime(2)*yldSDbPrime(3) - yldSDbPrime(6)**2)/2.0D0
      dHdSTildeDbPrime(3,2) = 
     & (yldSDbPrime(3)*yldSDbPrime(1) - yldSDbPrime(5)**2)/2.0D0
      dHdSTildeDbPrime(3,3) = 
     & (yldSDbPrime(1)*yldSDbPrime(2) - yldSDbPrime(4)**2)/2.0D0
      dHdSTildeDbPrime(3,4) = 
     & yldSDbPrime(5)*yldSDbPrime(6) - yldSDbPrime(3)*yldSDbPrime(4)
      dHdSTildeDbPrime(3,5) = 
     & yldSDbPrime(4)*yldSDbPrime(6) - yldSDbPrime(2)*yldSDbPrime(5)
      dHdSTildeDbPrime(3,6) = 
     & yldSDbPrime(4)*yldSDbPrime(5) - yldSDbPrime(1)*yldSDbPrime(6)

      dHdXPrime = MATMUL(dHdSTildePrime,dSTildedXPrime)
      dHdXDbPrime = MATMUL(dHdSTildeDbPrime,dSTildedXDbPrime)
      
      dSPridHPrime(:,:) = 0.0D0
      dSPridHDbPrime(:,:) = 0.0D0
      DO i=1,3
        denominator = 3.0D0*(yldSPriPrime(i)**2 - 
     &   2.0D0*invariantsPrime(1)*yldSPriPrime(i) - invariantsPrime(2))
        IF (denominator==0.0D0) THEN
          denominator = denominator + 1.0D-6
        END IF
        dSPridHPrime(i,3) = 2.0D0/denominator
        dSPridHPrime(i,2) = 
     &   dSPridHPrime(i,3)*3.0D0*yldSPriPrime(i)/2.0D0
        dSPridHPrime(i,1) = dSPridHPrime(i,2)*yldSPriPrime(i)
        denominator = 3.0D0*(yldSPriDbPrime(i)**2 - 2.0D0*
     &   invariantsDbPrime(1)*yldSPriDbPrime(i) - invariantsDbPrime(2))
        IF (denominator==0.0D0) THEN
          denominator = denominator + 1.0D-6
        END IF
        dSPridHDbPrime(i,3) = 2.0D0/denominator
        dSPridHDbPrime(i,2) = 
     &   dSPridHDbPrime(i,3)*3.0D0*yldSPriDbPrime(i)/2.0D0
        dSPridHDbPrime(i,1) = dSPridHDbPrime(i,2)*yldSPriDbPrime(i)
      END DO
      dSPridXPrime = MATMUL(dSPridHPrime,dHdXPrime)
      dSPridXDbPrime = MATMUL(dSPridHDbPrime,dHdXDbPrime)
      dPhidSPriPrime(:) = 0.0D0
      dPhidSPriDbPrime(:) = 0.0D0
      DO i=1,3
        dPhidSPriPrime(i) = YLDM*((yldSPriPrime(i) - 
     &   yldSPriDbPrime(1))*(ABS(yldSPriPrime(i) - 
     &   yldSPriDbPrime(1))**(YLDM - 2))+(yldSPriPrime(i) - 
     &   yldSPriDbPrime(2))*(ABS(yldSPriPrime(i) - 
     &   yldSPriDbPrime(2))**(YLDM - 2))+(yldSPriPrime(i) - 
     &   yldSPriDbPrime(3))*(ABS(yldSPriPrime(i) - yldSPriDbPrime(3))**
     &   (YLDM - 2)))
        dPhidSPriDbPrime(i) = YLDM*((yldSPriDbPrime(i) - 
     &   yldSPriPrime(1))*(ABS(yldSPriDbPrime(i) - 
     &   yldSPriPrime(1))**(YLDM - 2))+(yldSPriDbPrime(i) - 
     &   yldSPriPrime(2))*(ABS(yldSPriDbPrime(i) - 
     &   yldSPriPrime(2))**(YLDM - 2))+(yldSPriDbPrime(i) - 
     &   yldSPriPrime(3))*(ABS(yldSPriDbPrime(i) - yldSPriPrime(3))**
     &   (YLDM - 2)))
      END DO
      calc_dPhidX =  DOT_PRODUCT(dPhidSPriPrime,dSPridXPrime) + 
     & DOT_PRODUCT(dPhidSPriDbPrime,dSPridXDbPrime)
      RETURN
      END FUNCTION calc_dPhidX


      ! calculate differential of plastic potential
      SUBROUTINE calc_dGdS(hillParams,STRESS,dGdS)
      IMPLICIT NONE
      INTEGER i
      DOUBLE PRECISION hillParams(4),STRESS(6),dGdS(6),calc_eqGStress
      DO i=1,3
        dGdS(i) = (hillParams(MOD(i+1,3)) + 
     &   hillParams(MOD(i+2,3)))*STRESS(i) - 
     &   hillParams(MOD(i+1,3))*STRESS(MOD(i+2,3)) - 
     &   hillParams(MOD(i+2,3))*STRESS(MOD(i+1,3))
        dGdS(i+3) = 2.0D0*hillParams(4)*STRESS(i)
      END DO
      dGdS(:) = dGdS(:)*1.5D0/
     & (SUM(hillParams(1:3))*calc_eqGStress(hillParams,STRESS))
      RETURN
      END SUBROUTINE calc_dGdS


      ! calculate second differential of plastic potential
      SUBROUTINE calc_ddGddS(hillParams,STRESS,dGdS,ddGddS)
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION hillParams(4),STRESS(6),dGdS(6),ddGddS(6,6),
     & eqGStress,calc_eqGStress
      eqGStress = calc_eqGStress(hillParams,STRESS)
      DO i=1,3
        DO j=1,6
          ddGddS(i,j) = ((hillParams(MOD(i+1,3)) + 
     &     hillParams(MOD(i+2,3)))*STRESS(i) - hillParams(MOD(i+1,3))*
     &     STRESS(MOD(i+2,3)) - hillParams(MOD(i+2,3))*
     &     STRESS(MOD(i+1,3)))*(-1.0D0)*dGdS(j)/(eqGStress**2)
          ddGddS(i+3,j) = 
     &     2.0D0*hillParams(4)*STRESS(i)*(-1.0D0)*dGdS(j)/(eqGStress**2)
          ddGddS(i+3,i+3) = 
     &     ddGddS(i+3,i+3) + 2.0D0*hillParams(4)/eqGStress
          IF (j<4) THEN
            IF (i==j) THEN
              ddGddS(i,j) = ddGddS(i,j) + (hillParams(MOD(i+1,3)) + 
     &         hillParams(MOD(i+2,3)))/eqGStress
            ELSE
              ddGddS(i,j) = ddGddS(i,j) - hillParams(6-i-j)/eqGStress
            END IF
          END IF
        END DO
      END DO
      ddGddS(:,:) = ddGddS(:,:)*1.5D0/SUM(hillParams(1:3))
      RETURN
      END SUBROUTINE calc_ddGddS


      ! calculate principal value of stress tensor
      SUBROUTINE calc_Principal(stress,principal,invariants)
      IMPLICIT NONE
      DOUBLE PRECISION stress(6),principal(3),invariants(3),PI,radius,
     & cosTheta,theta
      PARAMETER(PI=ACOS(-1.0D0))
      CALL calc_Invariants(stress,invariants)
      radius = 2*SQRT(invariants(1)**2 + invariants(2))
      cosTheta = MAX(-1.0D0,MIN(1.0D0,((2*invariants(1)**3 + 
     & 3*invariants(1)*invariants(2) + 2*invariants(3))/2.0D0)/
     & SQRT((invariants(1)**2 + invariants(2))**3)))
      theta = ACOS(cosTheta)
      principal(1) = radius*COS(theta/3.0D0) + invariants(1)
      principal(2) = radius*COS((theta + 4*PI)/3.0D0) + invariants(1)
      principal(3) = radius*COS((theta + 2*PI)/3.0D0) + invariants(1)
      RETURN
      END SUBROUTINE calc_Principal


      ! calculate invariants of stress tensor
      SUBROUTINE calc_Invariants(stress,invariants)
      IMPLICIT NONE
      DOUBLE PRECISION stress(6),invariants(3)
      invariants(1) = (stress(1) + stress(2) + stress(3))/3.0D0
      invariants(2) = (stress(4)**2 + stress(5)**2 + stress(6)**2 - 
     & stress(1)*stress(2) - stress(2)*stress(3) - 
     & stress(3)*stress(1))/3.0D0
      invariants(3) = (2.0D0*stress(4)*stress(5)*stress(6) + 
     & stress(1)*stress(2)*stress(3) - stress(1)*(stress(6)**2) - 
     & stress(2)*(stress(5)**2) - stress(3)*(stress(4)**2))/2.0D0
      RETURN
      END SUBROUTINE calc_Invariants


      ! calculate dyad of vectors
      SUBROUTINE calc_dyad(dim,vec1,vec2,mat)
      IMPLICIT NONE
      INTEGER dim,i,j
      DOUBLE PRECISION vec1(dim),vec2(dim),mat(dim,dim)
      DO i=1,dim
        DO j=1,dim
          mat(i,j) = vec1(i)*vec2(j)
        END DO
      END DO
      RETURN
      END SUBROUTINE calc_dyad
