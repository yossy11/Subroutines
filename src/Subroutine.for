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
      PARAMETER(TOLER=1.0D-6,PI=180,YOUNG=7.03D4,POISSON=0.3D0,HARDK=615.3D0,HARDN=0.363D0,HARDSTRAIN0=7.61D-3,m=2)

      ! define variables, and their dimensions
      DOUBLE PRECISION lame,shearMod,eStrain,pStrain,eqPStrain,eqStress,
     & effectiveStress,yldT,yldCPrime,yldCDbPrime,yldPhi

      DIMENSION eStrain(NTENS),pStrain(NTENS)

      ! initialize DDSDDE
      lame = YOUNG*POISSON/((1 - 2*POISSON)*(1 + POISSON))
      shearMod = YOUNG/(2*(1 + POISSON))
      DDSDDE(:,:) = 0
      DDSDDE(1:NDI,1:NDI) = lame
      DO i=1,NDI
        DDSDDE(i,i) = lame + 2*shearMod
      ENDDO
      DO i=NDI+1,NTENS
        DDSDDE(i,i) = shearMod
      ENDDO

      ! read STATEV
      CALL ROTSIG(STATEV(1), DROT, eStrain, 2, NDI, NSHR)
      CALL ROTSIG(STATEV(NTENS+1), DROT, pStrain, 2, NDI, NSHR)
      eqPStrain = STATEV(2*NTENS+1)
      eStrain(:) = eStrain(:) + DSTRAN(:)

      ! calculate trial stress
      DO i=1,NTENS
        DO j=1,NTENS
          STRESS(i) = STRESS(i) + DDSDDE(i,j)*DSTRAN(j)
        ENDDO
      ENDDO

      ! calculate eqStress and effectiveStress
      yldT(:,:) = 0.0D0
      yldT(1:NDI,1:NDI) = -1.0D0
      DO i=1,NDI
        yldT(i,i) = 2.0D0
      ENDDO
      DO i=NDI+1,NTENS
        yldT(i,i) = 3.0D0
      ENDDO
      yldT(:,:) = yldT(:,:)/3.0D0
      yldCPrime(:,:) = 0.0D0
      yldCDbPrime(:,:) = 0.0D0

      ! ToDo: optimize C(minimize error function)

      yldSPrime = MATMUL(yldCPrime,MATMUL(yldT,STRESS))
      yldSDbPrime = MATMUL(yldCDbPrime,MATMUL(yldT,STRESS))
      ! ToDo: check the result of SPRINC
      CALL SPRINC(yldSPrime, yldSPriPrime, 1, NDI, NSHR)
      CALL SPRINC(yldSDbPrime, yldSPriDbPrime, 1, NDI, NSHR)

      DO i=1,3
        DO j=1,3
          yldPhi = yldPhi + ABS(yldSPriPrime(i) - yldSPriDbPrime(j))**m
        ENDDO
      ENDDO
      
      eqStress = (yldPhi/4)**(1/m)

      effectiveStress = HARDK*((HARDSTRAIN0 + eqpStrain)**HARDN)
      IF ((eqStress - effectiveStress) < TOLER) THEN
        CALL updateSTATEV(NTENS,STATEV,eStrain,pStrain,eqpStrain)
        RETURN
      ENDIF
      
      ! ToDo: if yield

      CALL updateSTATEV(NTENS,STATEV,eStrain,pStrain,eqpStrain)
      RETURN
      END

      ! update all state variable
      SUBROUTINE updateSTATEV(NTENS,STATEV,eStrain,pStrain,eqpStrain)
      DOUBLE PRECISION NTENS,STATEV,eStrain,pStrain,eqpStrain
      DIMENSION STATEV(4*NTENS+1),eStrain(NTENS),pStrain(NTENS)
      STATEV(1:NTENS)=eStrain
      STATEV(NTENS+1:2*NTENS)=pStrain
      STATEV(2*NTENS+1)=eqpStrain
      RETURN
      END SUBROUTINE updateSTATEV
