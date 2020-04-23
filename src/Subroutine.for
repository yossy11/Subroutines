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
      PARAMETER(TOLER=1.0D-6,PI=180,YOUNG=70300D0,POISSON=0.3)

      ! define variables, and their dimensions
      DOUBLE PRECISION lame,shearMod,eStrain,pStrain,backStress1,
     & backStress2,backStressTotal,eqStress,eqpStrain

      DIMENSION eStrain(NTENS),pStrain(NTENS),backStress1(NTENS),
     & backStress2(NTENS),backStressTotal(NTENS)

      ! initialize DDSDDE
      lame = YOUNG * POISSON / ((1 - 2 * POISSON) * (1 + POISSON))
      shearMod = YOUNG / (2 * (1 + POISSON))
      DDSDDE = 0
      DDSDDE(1:NDI,1:NDI) = lame
      DO i=1,NDI
        DDSDDE(i,i) = lame + 2 * shearMod
      ENDDO
      DO i=NDI+1,NTENS
        DDSDDE(i,i) = shearMod
      ENDDO

      ! read STATEV
      CALL ROTSIG(STATEV(1), DROT, eStrain, 2, NDI, NSHR)
      CALL ROTSIG(STATEV(NTENS+1), DROT, pStrain, 2, NDI, NSHR)
      CALL ROTSIG(STATEV(2*NTENS+1), DROT, backStress1, 1, NDI, NSHR)
      CALL ROTSIG(STATEV(3*NTENS+1), DROT, backStress2, 1, NDI, NSHR)
      backStressTotal = backStress1 + backStress2
      eqpStrain = STATEV(4*NTENS+1)
      eStrain = eStrain + DSTRAN

      !calculate trial stress
      DO i=1,NTENS
        DO j=1,NTENS
          STRESS(i) = STRESS(i) + DDSDDE(i,j) * DSTRAN(j)
        ENDDO
      ENDDO

      ! ToDo: calculate eqStress and yieldStress

      IF ((eqStress - yieldStress) .LT. TOLER) THEN
        call updateSTATEV(NTENS,STATEV,eStrain,pStrain,backStress1,
     &   backStress2,eqpStrain)
        RETURN
      ENDIF

      ! ToDo: if yield

      call updateSTATEV(NTENS,STATEV,eStrain,pStrain,backStress1,
     & backStress2,eqpStrain)
      RETURN
      END
      

      SUBROUTINE updateSTATEV(NTENS,STATEV,eStrain,pStrain,
     & backStress1,backStress2,eqpStrain)
      DOUBLE PRECISION NTENS,STATEV,eStrain,pStrain,backStress1,
     & backStress2,eqpStrain
      DIMENSION STATEV(4*NTENS+1),eStrain(NTENS),pStrain(NTENS),
     & backStress1(NTENS),backStress2(NTENS)
      STATEV(1:NTENS)=eStrain
      STATEV(NTENS+1:2*NTENS)=pStrain
      STATEV(2*NTENS+1:3*NTENS)=backStress1
      STATEV(3*NTENS+1:4*NTENS)=backStress2
      STATEV(4*NTENS+1)=eqpStrain
      RETURN
      END SUBROUTINE updateSTATEV