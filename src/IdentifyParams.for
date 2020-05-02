      ! identify anisotropic params of yld2004-18p
      PROGRAM IdentifyParams
      IMPLICIT NONE
      INTEGER YLDM,i
      DOUBLE PRECISION exStress(8),exLankford(8),yldCPrime(6,6),
     & yldCDbPrime(6,6)
      PARAMETER(YLDM=6)
      exStress(1) = 115.796080569375
      exStress(2) = 116.39694377773
      exStress(3) = 112.821037305888
      exStress(4) = 114.862943596084
      exStress(5) = 115.045917778501
      exStress(6) = 114.736215791944
      exStress(7) = 118.906776179804
      exStress(8) = 118.254840921828

      exLankford(1) = 0.69935
      exLankford(2) = 0.70275
      exLankford(3) = 0.7424
      exLankford(4) = 0.77555
      exLankford(5) = 0.78975
      exLankford(6) = 0.7594
      exLankford(7) = 0.7748
      exLankford(8) = 1.112

      yldCPrime(:,:) = 0.0D0
      yldCPrime(1:3,1:3) = 0.75D0
      yldCPrime(5,5) = 1.0D0
      yldCPrime(6,6) = 1.0D0
      yldCDbPrime(:,:) = 0.0D0
      yldCDbPrime(1:3,1:3) = 1.0D0
      yldCDbPrime(5,5) = 1.0D0
      yldCDbPrime(6,6) = 1.0D0
      DO i=1,3
        yldCPrime(i,i) = 0.0D0
        yldCDbPrime(i,i) = 0.0D0
      END DO
      CALL optimize_yldC(YLDM,exStress,exLankford,yldCPrime,yldCDbPrime)
      OPEN(20,FILE='anisotropicParams.csv')
      WRITE(20,*) 'yldCPrime'
      DO i=1,6
        WRITE(20,*) yldCPrime(i,:)
      END DO
      WRITE(20,*) 'yldCDbPrime'
      DO i=1,6
        WRITE(20,*) yldCDbPrime(i,:)
      END DO
      END PROGRAM IdentifyParams


      ! optimize yldCPrime and yldCDbPrime
      SUBROUTINE optimize_yldC(YLDM,exStress,exLankford,yldCPrime,
     & yldCDbPrime)
      IMPLICIT NONE
      INTEGER YLDM,iterationCount
      DOUBLE PRECISION exStress(8),exLankford(8),yldCPrime(6,6),
     & yldCDbPrime(6,6),LEARNINGRATE,TOL,WEIGHTS,WEIGHTR,WEIGHTB,
     & dCFactors(14),error,errorFunc
      PARAMETER(LEARNINGRATE=0.1D0,TOL=1.0D-6,WEIGHTS=1.0D0,
     & WEIGHTR=1.0D-1,WEIGHTB=1.0D-2)
      dCFactors(:) = 1.0D0
      error = errorFunc(YLDM,WEIGHTS,WEIGHTR,WEIGHTB,
     & exStress,exLankford,yldCPrime,yldCDbPrime)
      DO WHILE (error>TOL)
        CALL calc_dCFactors(YLDM,WEIGHTS,WEIGHTR,WEIGHTB,exStress,
     &   exLankford,yldCPrime,yldCDbPrime,dCFactors)
        yldCPrime(1,2) = yldCPrime(1,2) - LEARNINGRATE*dCFactors(1)
        yldCPrime(1,3) = yldCPrime(1,3) - LEARNINGRATE*dCFactors(2)
        yldCPrime(2,1) = yldCPrime(2,1) - LEARNINGRATE*dCFactors(3)
        yldCPrime(2,3) = yldCPrime(2,3) - LEARNINGRATE*dCFactors(4)
        yldCPrime(3,1) = yldCPrime(3,1) - LEARNINGRATE*dCFactors(5)
        yldCPrime(3,2) = yldCPrime(3,2) - LEARNINGRATE*dCFactors(6)
        yldCPrime(4,4) = yldCPrime(4,4) - LEARNINGRATE*dCFactors(7)
        yldCDbPrime(1,2) = yldCDbPrime(1,2) - LEARNINGRATE*dCFactors(8)
        yldCDbPrime(1,3) = yldCDbPrime(1,3) - LEARNINGRATE*dCFactors(9)
        yldCDbPrime(2,1) = yldCDbPrime(2,1) - LEARNINGRATE*dCFactors(10)
        yldCDbPrime(2,3) = yldCDbPrime(2,3) - LEARNINGRATE*dCFactors(11)
        yldCDbPrime(3,1) = yldCDbPrime(3,1) - LEARNINGRATE*dCFactors(12)
        yldCDbPrime(3,2) = yldCDbPrime(3,2) - LEARNINGRATE*dCFactors(13)
        yldCDbPrime(4,4) = yldCDbPrime(4,4) - LEARNINGRATE*dCFactors(14)
        error = errorFunc(YLDM,WEIGHTS,WEIGHTR,WEIGHTB,
     &   exStress,exLankford,yldCPrime,yldCDbPrime)
        WRITE(*,*) 'error: ',error
        iterationCount = iterationCount + 1
        IF (iterationCount >= 10000) THEN
          WRITE(*,*) 'too many'
          EXIT
        END IF
      END DO
      RETURN
      END SUBROUTINE optimize_yldC


      ! error function for anisotropic params
      DOUBLE PRECISION FUNCTION errorFunc(YLDM,WEIGHTS,WEIGHTR,WEIGHTB,
     & exStress,exLankford,yldCPrime,yldCDbPrime)
      IMPLICIT NONE
      INTEGER YLDM,i
      DOUBLE PRECISION WEIGHTS,WEIGHTR,WEIGHTB,exStress(8),
     & exLankford(8),yldCPrime(6,6),yldCDbPrime(6,6),PI,bufStress(6),
     & eqStress,calc_eqStress,lankford,calc_dPhidX
      PARAMETER(PI=ACOS(-1.0D0))
      errorFunc = 0.0D0
      DO i=1,7
        bufStress(:) = 0.0D0
        bufStress(1) = exStress(i)*COS((i-1)*PI/12.0D0)**2
        bufStress(2) = exStress(i)*SIN((i-1)*PI/12.0D0)**2
        bufStress(6) = exStress(i)*SIN((i-1)*PI/12.0D0)*
     &   COS((i-1)*PI/12.0D0)
        eqStress = calc_eqStress(YLDM,yldCPrime,yldCDbPrime,bufStress)
        IF (eqStress==0.0D0.or.ISNAN(eqStress)) THEN
          WRITE(*,*) 'error occurred, invalid eqStress: ',eqStress
        END IF
        lankford = -1.0D0 - 4.0D0*YLDM/(calc_dPhidX('szz  ',YLDM,
     &   yldCPrime,yldCDbPrime,bufStress)*exStress(i)/eqStress)
        errorFunc = errorFunc + WEIGHTS*(eqStress/exStress(i) - 1)**2
     &    + WEIGHTR*(lankford/exLankford(i) - 1)**2
      END DO
      bufStress(:) = 0.0D0
      bufStress(3) = exStress(8)
      eqStress = 1.0D0
      eqStress = calc_eqStress(YLDM,yldCPrime,yldCDbPrime,bufStress)
      lankford = 
     & calc_dPhidX('syy  ',YLDM,yldCPrime,yldCDbPrime,bufStress)
     & /calc_dPhidX('sxx  ',YLDM,yldCPrime,yldCDbPrime,bufStress)
      errorFunc = errorFunc + WEIGHTB*(((eqStress/exStress(8)) - 1)**2
     &  + (lankford/exLankford(8) - 1)**2)
      RETURN
      END FUNCTION errorFunc


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
      CASE ('sxx  ')
        dSTildedXPrime(1) = 0.0D0
        dSTildedXPrime(2) = -1.0D0*yldCPrime(2,1)
        dSTildedXPrime(3) = -1.0D0*yldCPrime(3,1)
        dSTildedXDbPrime(1) = 0.0D0
        dSTildedXDbPrime(2) = -1.0D0*yldCDbPrime(2,1)
        dSTildedXDbPrime(3) = -1.0D0*yldCDbPrime(3,1)

      CASE ('syy  ')
        dSTildedXPrime(1) = -1.0D0*yldCPrime(1,2)
        dSTildedXPrime(2) = 0.0D0
        dSTildedXPrime(3) = -1.0D0*yldCPrime(3,2)
        dSTildedXDbPrime(1) = -1.0D0*yldCDbPrime(1,2)
        dSTildedXDbPrime(2) = 0.0D0
        dSTildedXDbPrime(3) = -1.0D0*yldCDbPrime(3,2)

      CASE ('szz  ')
        dSTildedXPrime(1) = -1.0D0*yldCPrime(1,3)
        dSTildedXPrime(2) = -1.0D0*yldCPrime(2,3)
        dSTildedXPrime(3) = 0.0D0
        dSTildedXDbPrime(1) = -1.0D0*yldCDbPrime(1,3)
        dSTildedXDbPrime(2) = -1.0D0*yldCDbPrime(2,3)
        dSTildedXDbPrime(3) = 0.0D0

      CASE ('c12  ')
        dSTildedXPrime(1) = -1.0D0*yldSPrime(2)

      CASE ('c13  ')
        dSTildedXPrime(1) = -1.0D0*yldSPrime(3)

      CASE ('c21  ')
        dSTildedXPrime(2) = -1.0D0*yldSPrime(1)

      CASE ('c23  ')
        dSTildedXPrime(2) = -1.0D0*yldSPrime(3)

      CASE ('c31  ')
        dSTildedXPrime(3) = -1.0D0*yldSPrime(1)

      CASE ('c32  ')
        dSTildedXPrime(3) = -1.0D0*yldSPrime(2)

      CASE ('c44  ')
        dSTildedXPrime(4) = yldSPrime(4)

      CASE ('c12Db')
        dSTildedXDbPrime(1) = -1.0D0*yldSDbPrime(2)

      CASE ('c13Db')
        dSTildedXDbPrime(1) = -1.0D0*yldSDbPrime(3)

      CASE ('c21Db')
        dSTildedXDbPrime(2) = -1.0D0*yldSDbPrime(1)

      CASE ('c23Db')
        dSTildedXDbPrime(2) = -1.0D0*yldSDbPrime(3)

      CASE ('c31Db')
        dSTildedXDbPrime(3) = -1.0D0*yldSDbPrime(1)

      CASE ('c32Db')
        dSTildedXDbPrime(3) = -1.0D0*yldSDbPrime(2)

      CASE ('c44Db')
        dSTildedXDbPrime(4) = yldSDbPrime(4)

      CASE DEFAULT
        WRITE(*,*) 'Invalid orientation'
        STOP
        ! CALL XIT
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


      SUBROUTINE calc_dCFactors(YLDM,WEIGHTS,WEIGHTR,WEIGHTB,exStress,
     & exLankford,yldCPrime,yldCDbPrime,dCFactors)
      IMPLICIT NONE
      CHARACTER(5) orientations(14)
      INTEGER YLDM,i,j
      DOUBLE PRECISION WEIGHTS,WEIGHTR,WEIGHTB,exStress(8),
     & exLankford(8),yldCPrime(6,6),yldCDbPrime(6,6),dCFactors(14),PI,
     & bufStress(6),eqStress,lankford,calc_eqStress,bufdPhidX,
     & calc_dPhidX,dFdPhi,dSdCPrime,dSdCDbPrime,dRdCPrime,dRdCDbPrime
      PARAMETER(PI=ACOS(-1.0D0))
      orientations = ['c12  ','c13  ','c21  ','c23  ','c31  ','c32  ',
     & 'c44  ','c12Db','c13Db','c21Db','c23Db','c31Db','c32Db','c44Db']
      DO i=1,6
        DO j=1,7
          bufStress(:) = 0.0D0
          bufStress(1) = exStress(j)*(COS((j-1)*PI/12.0D0)**2)
          bufStress(2) = exStress(j)*(SIN((j-1)*PI/12.0D0)**2)
          bufStress(6) = 
     &     exStress(j)*SIN((j-1)*PI/12.0D0)*COS((j-1)*PI/12.0D0)
          eqStress = calc_eqStress(YLDM,yldCPrime,yldCDbPrime,bufStress)
          bufdPhidX = 
     &     calc_dPhidX('szz  ',YLDM,yldCPrime,yldCDbPrime,bufStress)
          lankford = 
     &     -1.0D0 - 4.0D0*YLDM/(bufdPhidX*exStress(i)/eqStress)
          dFdPhi = eqStress/(YLDM*(eqStress**YLDM))
          dSdCPrime = dFdPhi*calc_dPhidX(orientations(i),YLDM,
     &     yldCPrime,yldCDbPrime,bufStress)
          dSdCDbPrime = dFdPhi*calc_dPhidX(orientations(i+7),YLDM,
     &     yldCPrime,yldCDbPrime,bufStress)
          dRdCPrime = -4.0D0*YLDM*dSdCPrime/(exStress(j)*bufdPhidX)
          dRdCDbPrime = -4.0D0*YLDM*dSdCDbPrime/(exStress(j)*bufdPhidX)
          dCFactors(i) = dCFactors(i) + WEIGHTS*2*dSdCPrime*
     &     ((eqStress/exStress(j)) - 1)/exStress(j) + WEIGHTR*2*
     &     dRdCPrime*((lankford/exLankford(j)) - 1)/exLankford(j)
          dCFactors(i+7) = dCFactors(i+7) + WEIGHTS*2*dSdCDbPrime*
     &     ((eqStress/exStress(j)) - 1)/exStress(j) + WEIGHTR*2*
     &     dRdCDbPrime*((lankford/exLankford(j)) - 1)/exLankford(j)
        END DO
      END DO
      bufStress(:) = 0.0D0
      bufStress(3) = exStress(8)
      eqStress = calc_eqStress(YLDM,yldCPrime,yldCDbPrime,bufStress)
      lankford = calc_dPhidX('syy  ',YLDM,yldCPrime,yldCDbPrime,
     & bufStress)/calc_dPhidX('sxx  ',YLDM,yldCPrime,yldCDbPrime,
     & bufStress)
      dFdPhi = eqStress/(YLDM*(eqStress**YLDM))
      dSdCPrime = dFdPhi*calc_dPhidX(orientations(7),
     & YLDM,yldCPrime,yldCDbPrime,bufStress)
      dSdCDbPrime = dFdPhi*calc_dPhidX(orientations(14),
     & YLDM,yldCPrime,yldCDbPrime,bufStress)
      dCFactors(7) = dCFactors(7) + 
     & WEIGHTB*2*dSdCPrime*((eqStress/exStress(8)) - 1)/exStress(8)
      dCFactors(14) = dCFactors(14) + 
     & WEIGHTB*2*dSdCDbPrime*((eqStress/exStress(8)) - 1)/exStress(8)
      RETURN
      END SUBROUTINE calc_dCFactors


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
  
