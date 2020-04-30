      ! This is UMAT subroutine made by yossy11
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

      ! define constants
      ! ToDo: separate params to props
      PARAMETER(TOLER=1.0D-6,PI=180,YOUNG=7.03D4,POISSON=0.3D0,HARDK=615.3D0,HARDN=0.363D0,HARDSTRAIN0=7.61D-3,YLDM=2)

      ! define variables, and their dimensions
      DOUBLE PRECISION lame,shearMod,eStrain(NTENS),pStrain(NTENS),
     & eqPStrain,eqStress,effectiveStress,
     & yldCPrime(NTENS,NTENS),yldCDbPrime(NTENS,NTENS),exStress(8),exLankford(8)

      exStress(:) = PROPS(1:8)
      exLankford(:) = PROPS(9:16)

      ! initialize DDSDDE
      lame = YOUNG*POISSON/((1 - 2*POISSON)*(1 + POISSON))
      shearMod = YOUNG/(2*(1 + POISSON))
      DDSDDE(:,:) = 0.0D0
      DDSDDE(1:NDI,1:NDI) = lame
      DO i=1,NDI
        DDSDDE(i,i) = lame + 2*shearMod
      END DO
      DO i=NDI+1,NTENS
        DDSDDE(i,i) = shearMod
      END DO      

      ! read STATEV
      CALL ROTSIG(STATEV(1), DROT, eStrain, 2, NDI, NSHR)
      CALL ROTSIG(STATEV(NTENS+1), DROT, pStrain, 2, NDI, NSHR)
      eqPStrain = STATEV(2*NTENS+1)
      eStrain(:) = eStrain(:) + DSTRAN(:)

      ! calculate trial stress
      STRESS(:) = STRESS(:) + MATMUL(DDSDDE,DSTRAN)

      ! calculate eqStress and effectiveStress
      yldCPrime(:,:) = 0.0D0
      yldCDbPrime(:,:) = 0.0D0
      CALL optimizeYldC(yldCPrime,yldCDbPrime)
      eqStress = calcEqStress(YLDM,NDI,NSHR,STRESS,yldCPrime,
     & yldCDbPrime)
      effectiveStress = HARDK*((HARDSTRAIN0 + eqpStrain)**HARDN)

      ! if not yield
      IF ((eqStress - effectiveStress) < TOLER) THEN
        CALL updateSTATEV(NTENS,STATEV,eStrain,pStrain,eqpStrain)
        RETURN
      END IF
      
      ! ToDo: if yield

      CALL updateSTATEV(NTENS,STATEV,eStrain,pStrain,eqpStrain)
      RETURN
      END SUBROUTINE UMAT


      ! update all state variable
      SUBROUTINE updateSTATEV(NTENS,STATEV,eStrain,pStrain,eqpStrain)
      INTEGER NTENS
      DOUBLE PRECISION STATEV(2*NTENS+1),eStrain(NTENS),pStrain(NTENS)
     & eqpStrain
      STATEV(1:NTENS)=eStrain
      STATEV(NTENS+1:2*NTENS)=pStrain
      STATEV(2*NTENS+1)=eqpStrain
      RETURN
      END SUBROUTINE updateSTATEV


      ! optimize yldCPrime and yldCDbPrime
      SUBROUTINE optimizeYldC(NTENS,yldCPrime,yldCDbPrime)
      PARAMETER(learningRate=0.1D0,tol=1.0D-6)
      INTEGER NTENS,numIteration=50000
      DOUBLE PRECISION yldCPrime(NTENS,NTENS),yldCDbPrime(NTENS,NTENS),dCFactors(14)
      DO WHILE (errorFunc(YLDM,NDI,NSHR,yldCPrime,yldCDbPrime,exStress,exLankford)>tol)
        CALl calcDCFactors(dCFactors)
        yldCPrime(1,2) = yldCPrime(1,2) - learningRate*dCFactors(1)
        yldCPrime(1,3) = yldCPrime(1,3) - learningRate*dCFactors(2)
        yldCPrime(2,1) = yldCPrime(2,1) - learningRate*dCFactors(3)
        yldCPrime(2,3) = yldCPrime(2,3) - learningRate*dCFactors(4)
        yldCPrime(3,1) = yldCPrime(3,1) - learningRate*dCFactors(5)
        yldCPrime(3,2) = yldCPrime(3,2) - learningRate*dCFactors(6)
        yldCPrime(4,4) = yldCPrime(4,4) - learningRate*dCFactors(7)
        yldCDbPrime(1,2) = yldCDbPrime(1,2) - learningRate*dCFactors(8)
        yldCDbPrime(1,3) = yldCDbPrime(1,3) - learningRate*dCFactors(9)
        yldCDbPrime(2,1) = yldCDbPrime(2,1) - learningRate*dCFactors(10)
        yldCDbPrime(2,3) = yldCDbPrime(2,3) - learningRate*dCFactors(11)
        yldCDbPrime(3,1) = yldCDbPrime(3,1) - learningRate*dCFactors(12)
        yldCDbPrime(3,2) = yldCDbPrime(3,2) - learningRate*dCFactors(13)
        yldCDbPrime(4,4) = yldCDbPrime(4,4) - learningRate*dCFactors(14)
      END DO
      RETURN
      END SUBROUTINE optimizeYldC


      ! error function for anisotropic params
      DOUBLE PRECISION FUNCTION errorFunc(YLDM,NDI,NSHR,yldCPrime,
     & yldCDbPrime,exStress,exLankford)
      INTEGER YLDM,NDI,NSHR
      DOUBLE PRECISION yldCPrime(6,6),yldCDbPrime(6,6),exStress(8),exLankford(8)
      PI = ACOS(-1)
      WEIGHTS = 1.0D0
      WEIGHTR = 0.1D0
      WEIGHTB = 0.1D-1
      DO i=1,7
        bufStress(:) = 0.0D0
        bufStress(1) = exStress(i)*(COS((i-1)*PI/12.0D0)**2)
        bufStress(2) = exStress(i)*(SIN((i-1)*PI/12.0D0)**2)
        bufStress(6) = exStress(i)*SIN((i-1)*PI/12.0D0)*
     &   COS((i-1)*PI/12.0D0)
        eqStress = calcEqStress(YLDM,NDI,NSHR,bufStress,yldCPrime,
     &   yldCDbPrime)
        lankford = -1.0D0 - 4.0D0*YLDM/(calcDPhiDSDev('zz',YLDM,NDI,
     &   NSHR,yldCPrime,yldCDbPrime,bufStress)*exStress(i)/eqStress)
        errorFunc = errorFunc + WEIGHTS*((eqStress/exStress(i) - 1)**2)
     &    + WEIGHTR*((lankford/exLankford(i) - 1)**2)
      END DO
      bufStress(:) = 0.0D0
      bufStress(3) = exStress(8)
      eqStress = calcEqStress(YLDM,NDI,NSHR,bufStress,yldCPrime,
     & yldCDbPrime)
      lankford = calcDPhiDSDev('yy',YLDM,NDI,NSHR,yldCPrime,
     & yldCDbPrime,bufStress)/calcDPhiDSDev('xx',YLDM,NDI,
     & NSHR,yldCPrime,yldCDbPrime,bufStress)
      errorFunc = errorFunc + WEIGHTB*(((eqStress/exStress(8)) - 1)**2
     &  + (lankford/exLankford(8) - 1)**2)
      RETURN
      END FUNCTION errorFunc


      ! calculate equivalent stress
      DOUBLE PRECISION FUNCTION calcEqStress(YLDM,NDI,NSHR,STRESS,yldCPrime,
     & yldCDbPrime)
      INTEGER YLDM,NDI,NSHR
      DOUBLE PRECISION STRESS(6),yldCPrime(6,6),
     & yldCDbPrime(6,6),yldT(6,6),yldSPrime(6),yldSDbPrime(6),yldSPriPrime(3),
     & yldSPriDbPrime(3),yldPhi
      yldT(:,:) = 0.0D0
      yldT(1:NDI,1:NDI) = -1.0D0
      DO i=1,NDI
        yldT(i,i) = 2.0D0
      END DO
      DO i=NDI+1,NTENS
        yldT(i,i) = 3.0D0
      END DO
      yldT(:,:) = yldT(:,:)/3.0D0
      yldSPrime = MATMUL(yldCPrime,MATMUL(yldT,STRESS))
      yldSDbPrime = MATMUL(yldCDbPrime,MATMUL(yldT,STRESS))
      CALL SPRINC(yldSPrime, yldSPriPrime, 1, NDI, NSHR)
      CALL SPRINC(yldSDbPrime, yldSPriDbPrime, 1, NDI, NSHR)
      DO i=1,3
        DO j=1,3
          yldPhi = yldPhi + ABS(yldSPriPrime(i) - yldSPriDbPrime(j))**YLDM
        END DO
      END DO
      calcEqStress = (yldPhi/4.0D0)**(1/YLDM)
      RETURN
      END FUNCTION calcEqStress


      ! calculate differential
      DOUBLE PRECISION FUNCTION calcDPhiDSDev(orientation,YLDM,NDI,NSHR,
     & yldCPrime,yldCDbPrime,bufStress)
      yldT(:,:) = 0.0D0
      yldT(1:NDI,1:NDI) = -1.0D0
      DO i=1,NDI
        yldT(i,i) = 2.0D0
      END DO
      DO i=NDI+1,NTENS
        yldT(i,i) = 3.0D0
      END DO
      yldT(:,:) = yldT(:,:)/3.0D0
      yldSPrime = MATMUL(yldCPrime,MATMUL(yldT,bufStress))
      yldSDbPrime = MATMUL(yldCDbPrime,MATMUL(yldT,bufStress))
      CALL SPRINC(yldSPrime, yldSPriPrime, 1, NDI, NSHR)
      CALL SPRINC(yldSDbPrime, yldSPriDbPrime, 1, NDI, NSHR)
      CALL calcInvariants(yldSPrime,invariantsPrime)
      CALL calcInvariants(yldSDbPrime,invariantsDbPrime)
      dHDSTildePrime(1,:) = 1.0D0/3.0D0
      dHDSTildeDbPrime(1,:) = 1.0D0/3.0D0
      DO i=1,3
        dHDSTildePrime(2,i) = (yldSPrime(i) - SUM(yldSPrime(1:3)))/3.0D0
        dHDSTildePrime(3,i) = (PRODUCT(yldSPrime(1:3))/yldSPrime(i) - yldSPrime(7-i)**2)/2.0D0
        dHDSTildeDbPrime(2,i) = (yldSDbPrime(i) - SUM(yldSDbPrime(1:3)))/3.0D0
        dHDSTildeDbPrime(3,i) = (PRODUCT(yldSDbPrime(1:3))/yldSDbPrime(i) - yldSDbPrime(7-i)**2)/2.0D0
      END DO

      dSTildeDSDevPrime(:) = 0.0D0
      dSTildeDSDevDbPrime(:) = 0.0D0
      SELECT CASE (orientation)
      CASE ('xx')
        dSTildeDSDevPrime(1) = 0.0D0
        dSTildeDSDevPrime(2) = -1.0D0*yldCPrime(2,1)
        dSTildeDSDevPrime(3) = -1.0D0*yldCPrime(3,1)
        dSTildeDSDevDbPrime(1) = 0.0D0
        dSTildeDSDevDbPrime(2) = -1.0D0*yldCDbPrime(2,1)
        dSTildeDSDevDbPrime(3) = -1.0D0*yldCDbPrime(3,1)

      CASE ('yy')
        dSTildeDSDevPrime(1) = -1.0D0*yldCPrime(1,2)
        dSTildeDSDevPrime(2) = 0.0D0
        dSTildeDSDevPrime(3) = -1.0D0*yldCPrime(3,2)
        dSTildeDSDevDbPrime(1) = -1.0D0*yldCDbPrime(1,2)
        dSTildeDSDevDbPrime(2) = 0.0D0
        dSTildeDSDevDbPrime(3) = -1.0D0*yldCDbPrime(3,2)

      CASE ('zz')
        dSTildeDSDevPrime(1) = -1.0D0*yldCPrime(1,3)
        dSTildeDSDevPrime(2) = -1.0D0*yldCPrime(2,3)
        dSTildeDSDevPrime(3) = 0.0D0
        dSTildeDSDevDbPrime(1) = -1.0D0*yldCDbPrime(1,3)
        dSTildeDSDevDbPrime(2) = -1.0D0*yldCDbPrime(2,3)
        dSTildeDSDevDbPrime(3) = 0.0D0

      CASE ('c12')
        dSTildeDSDevPrime(1) = -1.0D0*yldSPrime(2)

      CASE ('c13')
        dSTildeDSDevPrime(1) = -1.0D0*yldSPrime(3)

      CASE ('c21')
        dSTildeDSDevPrime(2) = -1.0D0*yldSPrime(1)

      CASE ('c23')
        dSTildeDSDevPrime(2) = -1.0D0*yldSPrime(3)

      CASE ('c31')
        dSTildeDSDevPrime(3) = -1.0D0*yldSPrime(1)

      CASE ('c32')
        dSTildeDSDevPrime(3) = -1.0D0*yldSPrime(2)

      CASE ('c12Db')
        dSTildeDSDevDbPrime(1) = -1.0D0*yldSDbPrime(2)

      CASE ('c13Db')
        dSTildeDSDevDbPrime(1) = -1.0D0*yldSDbPrime(3)

      CASE ('c21Db')
        dSTildeDSDevDbPrime(2) = -1.0D0*yldSDbPrime(1)

      CASE ('c23Db')
        dSTildeDSDevDbPrime(2) = -1.0D0*yldSDbPrime(3)

      CASE ('c31Db')
        dSTildeDSDevDbPrime(3) = -1.0D0*yldSDbPrime(1)

      CASE ('c32Db')
        dSTildeDSDevDbPrime(3) = -1.0D0*yldSDbPrime(2)

      CASE ('c44')
      CASE ('c44Db')

      CASE DEFAULT
        WRITE(6,*) 'Invalid orientation'
        CALL XIT
      END SELECT

      dHDSDevPrime = MATMUL(dHDSTildePrime,dSTildeDSDevPrime)
      dHDSDevDbPrime = MATMUL(dHDSTildeDbPrime,dSTildeDSDevDbPrime)
      IF (orientation=='c44') THEN
        dHDSDevPrime(0) = 0.0D0
        dHDSDevPrime(1) = 2.0D0*(yldSPrime(4)**2)/3.0D0
        dHDSDevPrime(2) = (yldSPrime(5)*yldSPrime(6) - yldSPrime(3)*yldSPrime(4))*yldSPrime(4)
      ELSE IF (orientation=='c44Db') THEN
        dHDSDevDbPrime(0) = 0.0D0
        dHDSDevDbPrime(1) = 2.0D0*(yldSDbPrime(4)**2)/3.0D0
        dHDSDevDbPrime(2) = (yldSDbPrime(5)*yldSDbPrime(6) - yldSDbPrime(3)*yldSDbPrime(4))*yldSDbPrime(4)
      END IF
      dSPriDHPrime(:,:) = 0.0D0
      dSPriDHDbPrime(:,:) = 0.0D0
      DO i=1,3
        dSPriDHPrime(i,3) = 2.0D0/(3.0D0*(yldSPriPrime(i)**2 - 
     &   2.0D0*invariantsPrime(1)*yldSPriPrime(i) - invariantsPrime(2)))
        dSPriDHPrime(i,2) = dSPriDHPrime(i,3)*3.0D0*
     &   yldSPriPrime(i)/2.0D0
        dSPriDHPrime(i,1) = dSPriDHPrime(i,2)*yldSPriPrime(i)
        dSPriDHDbPrime(i,3) = 2.0D0/(3.0D0*(yldSPriDbPrime(i)**2 - 
     &   2.0D0*invariantsDbPrime(1)*yldSPriDbPrime(i) - invariantsDbPrime(2)))
        dSPriDHDbPrime(i,2) = dSPriDHDbPrime(i,3)*3.0D0*
     &   yldSPriDbPrime(i)/2.0D0
        dSPriDHDbPrime(i,1) = dSPriDHDbPrime(i,2)*yldSPriDbPrime(i)
      END DO
      dSPriDSDevPrime = MATMUL(dSPriDHPrime,dHDSDevPrime)
      dSPriDSDevDbPrime = MATMUL(dSPriDHDbPrime,dHDSDevDbPrime)

      dPhiDSPriPrime(:) = 0.0D0
      dPhiDSPriDbPrime(:) = 0.0D0
      DO i=1,3
        dPhiDSPriPrime(i) = YLDM*((yldSPriPrime(i) - 
     &   yldSPriDbPrime(1))*(ABS(yldSPriPrime(i) - 
     &   yldSPriDbPrime(1))**(YLDM - 2))+(yldSPriPrime(i) - 
     &   yldSPriDbPrime(2))*(ABS(yldSPriPrime(i) - 
     &   yldSPriDbPrime(2))**(YLDM - 2))+(yldSPriPrime(i) - 
     &   yldSPriDbPrime(3))*(ABS(yldSPriPrime(i) - yldSPriDbPrime(3))**
     &   (YLDM - 2)))
        dPhiDSPriDbPrime(i) = YLDM*((yldSPriDbPrime(i) - 
     &   yldSPriPrime(1))*(ABS(yldSPriDbPrime(i) - 
     &   yldSPriPrime(1))**(YLDM - 2))+(yldSPriDbPrime(i) - 
     &   yldSPriPrime(2))*(ABS(yldSPriDbPrime(i) - 
     &   yldSPriPrime(2))**(YLDM - 2))+(yldSPriDbPrime(i) - 
     &   yldSPriPrime(3))*(ABS(yldSPriDbPrime(i) - yldSPriPrime(3))**
     &   (YLDM - 2)))
      END DO
      calcDPhiDSDev =  DOT_PRODUCT(dPhiDSPriPrime,dSPriDSDevPrime) + DOT_PRODUCT(dPhiDSPriDbPrime,dSPriDSDevDbPrime)
      RETURN
      END FUNCTION calcDPhiDSDev


      SUBROUTINE calcDCFactors(dCFactors)  
      orientations = ['c12','c13','c21','c23','c31','c32','c44','c12Db','c13Db','c21Db','c23Db','c31Db','c32Db','c44Db']
      DO i=1,6
        DO j=1,7
          bufStress(:) = 0.0D0
          bufStress(1) = exStress(j)*(COS((j-1)*PI/12.0D0)**2)
          bufStress(2) = exStress(j)*(SIN((j-1)*PI/12.0D0)**2)
          bufStress(6) = exStress(j)*SIN((j-1)*PI/12.0D0)*COS((j-1)*PI/12.0D0)
          eqStress = calcEqStress(YLDM,NDI,NSHR,bufStress,yldCPrime,yldCDbPrime)
          lankford = -1.0D0 - 4.0D0*YLDM/(calcDPhiDSDev('zz',YLDM,NDI,NSHR,yldCPrime,yldCDbPrime,bufStress)*exStress(i)/eqStress)
          dFDPhi = eqStress/(YLDM*(eqStress**YLDM))
          dSDCPrime = dFDPhi * calcDPhiDSDev(orientations(i),YLDM,NDI,NSHR,yldCPrime,yldCDbPrime,bufStress)
          dSDCDbPrime = dFDPhi * calcDPhiDSDev(orientations(i+7),YLDM,NDI,NSHR,yldCPrime,yldCDbPrime,bufStress)
          dRDCPrime = -4.0D0*YLDM*dSDCPrime/(exStress(j)*calcDPhiDSDev('zz',YLDM,NDI,NSHR,yldCPrime,yldCDbPrime,bufStress))
          dRDCDbPrime = -4.0D0*YLDM*dSDCDbPrime/(exStress(j)*calcDPhiDSDev('zz',YLDM,NDI,NSHR,yldCPrime,yldCDbPrime,bufStress))
          dCFactors(i) = dCFactors(i) + WEIGHTS*2*dSDCPrime*((eqStress/exStress(j)) - 1)/exStress(j) + WEIGHTR*2*dRDCPrime*((lankford/exLankford(j)) - 1)/exLankford(j)
          dCFactors(i+7) = dCFactors(i+7) + WEIGHTS*2*dSDCDbPrime*((eqStress/exStress(j)) - 1)/exStress(j) + WEIGHTR*2*dRDCDbPrime*((lankford/exLankford(j)) - 1)/exLankford(j)
        END DO
      END DO
      bufStress(:) = 0.0D0
      bufStress(3) = exStress(8)
      eqStress = calcEqStress(YLDM,NDI,NSHR,bufStress,yldCPrime,yldCDbPrime)
      lankford = calcDPhiDSDev('yy',YLDM,NDI,NSHR,yldCPrime,
     & yldCDbPrime,bufStress)/calcDPhiDSDev('xx',YLDM,NDI,
     & NSHR,yldCPrime,yldCDbPrime,bufStress)
      dFDPhi = eqStress/(YLDM*(eqStress**YLDM))
      dSDCPrime = dFDPhi * calcDPhiDSDev(orientations(7),YLDM,NDI,NSHR,yldCPrime,yldCDbPrime,bufStress)
      dSDCDbPrime = dFDPhi * calcDPhiDSDev(orientations(14),YLDM,NDI,NSHR,yldCPrime,yldCDbPrime,bufStress)
      dCFactors(7) = dCFactors(7) + WEIGHTB*2*dSDCPrime*((eqStress/exStress(8)) - 1)/exStress(8)
      dCFactors(14) = dCFactors(14) + WEIGHTB*2*dSDCDbPrime*((eqStress/exStress(8)) - 1)/exStress(8)
      RETURN
      END SUBROUTINE calcDCFactors


      ! calculate invariants of stress tensor
      SUBROUTINE calcInvariants(stress,invariants)
      DOUBLE PRECISION stress(6),invariants(2)
      invariants(1) = (stress(1) + stress(2) + stress(3))/3.0D0
      invariants(2) = (stress(4)**2 + stress(5)**2 + stress(6)**2 - 
     & stress(1)*stress(2) - stress(2)*stress(3) - 
     & stress(3)*stress(1))/3.0D0
      RETURN
      END SUBROUTINE calcInvariants
