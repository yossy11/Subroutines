      ! This is sample UMAT subroutine made by yossy11, on 20200419
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

      PARAMETER(TOLER=1.0D-6,PI=180,YOUNG=70300D0,POISSON=0.3D0)
      DOUBLE PRECISION val,lame,shearMod

      ! initialize DDSDDE
      lame = YOUNG * POISSON / ((1 - 2 * POISSON) * (1 + POISSON))
      shearMod = YOUNG / (2 * (1 + POISSON))
      DDSDDE(:,:) = 0.0D0
      DDSDDE(1:NDI,1:NDI) = lame
      DO i=1,NDI
        DDSDDE(i,i) = lame + 2*shearMod
      END DO
      DO i=NDI+1,NTENS
        DDSDDE(i,i) = shearMod
      END DO
      
      ! update STRESS
      STRESS(:) = STRESS(:) + MATMUL(DDSDDE,DSTRAN)

      val = 5.0D0
      val = testFunc(.FALSE.,val)
      WRITE(6,*),"Hello",val
      RETURN
      END
      
      DOUBLE PRECISION FUNCTION testFunc(isAH,value,opValue)
        DOUBLE PRECISION value, opValue
        OPTIONAL opValue
        LOGICAL isAH
        IF (isAH) THEN
          IF (.NOT. PRESENT(opValue)) THEN
            WRITE(6,*),"opValue is required in testFunc"
            RETURN
          ENDIF
          testFunc = value * opValue
        ELSE
          testFunc = value ** 2.0D0
        ENDIF
        RETURN
      END FUNCTION testFunc

      SUBROUTINE testSub(isAH,value,opValue)
        DOUBLE PRECISION value, opValue
        OPTIONAL opValue
        LOGICAL isAH
        IF (isAH) THEN
          IF (.NOT. PRESENT(opValue)) THEN
            WRITE(6,*),"opValue is required in testSub"
            RETURN
          ENDIF
          value = value * opValue
        ELSE
          value = value ** 2.0D0
        ENDIF
        RETURN
      END SUBROUTINE testSub