      ! This is UMAT base subroutine made by yossy11
      ! yield function : Yld2004-18P
      ! plastic potential : Yld2004-18P
      ! hardening rule : Swift
      ! flow rule : non-AFR
      ! integration algorithm : CCP(convex cutting plane)

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
      DOUBLE PRECISION TOLER,YOUNG,POISSON,lame,shearMod,eStrain(6),
     & pStrain(6),eqpStrain,totalStrain(6),eqStress,eqGStress,
     & flowStress,calc_eqStress,calc_eqGStress,calc_FlowStress,dfds(6),
     & dgds(6),H,calc_H,dLambda,F,lambda

      DOUBLE PRECISION HARDK,HARDN,HARDSTRAIN0,YLDM,f_yldCPrime(6,6),
     & f_yldCDbPrime(6,6),g_yldCPrime(6,6),g_yldCDbPrime(6,6)

      ! define constants
      PARAMETER(TOLER=1.0D-5,YOUNG=6.9D4,POISSON=0.33D0)

      ! define hardening constants
      PARAMETER(HARDK=646.0D0,HARDN=0.227D0,HARDSTRAIN0=2.5D-2,
     & YLDM=8.0D0)

      ! anisotropic params
      f_yldCPrime(:,:) = 0.0D0
      f_yldCDbPrime(:,:) = 0.0D0

      f_yldCPrime(1,2) = -0.056940277665336686D0
      f_yldCPrime(1,3) = -0.9859971641591349D0
      f_yldCPrime(2,1) = -0.1986277232332533D0
      f_yldCPrime(2,3) = -0.9313609018430112D0
      f_yldCPrime(3,1) = -0.5871124291038275D0
      f_yldCPrime(3,2) = -1.330942716974344D0
      f_yldCPrime(4,4) = 0.944649833438939D0
      f_yldCPrime(5,5) = 1.090740680806561D0
      f_yldCPrime(6,6) = 1.0925155792832137D0

      f_yldCDbPrime(1,2) = -1.0763323731674885D0
      f_yldCDbPrime(1,3) = -0.4323653254398786D0
      f_yldCDbPrime(2,1) = -0.613688607514406D0
      f_yldCDbPrime(2,3) = -0.9659074930180515D0
      f_yldCDbPrime(3,1) = -1.1519059774170437D0
      f_yldCDbPrime(3,2) = 0.08362991282554624D0
      f_yldCDbPrime(4,4) = 1.4567012369610057D0
      f_yldCDbPrime(5,5) = 1.1752482654026506D0
      f_yldCDbPrime(6,6) = 1.116415832572442D0

      g_yldCPrime(:,:) = 0.0D0
      g_yldCDbPrime(:,:) = 0.0D0

      g_yldCPrime(1,2) = 0.0681890550605269D0
      g_yldCPrime(1,3) = -0.8728581783781829D0
      g_yldCPrime(2,1) = -0.07883817293470376D0
      g_yldCPrime(2,3) = -0.9800868482532625D0
      g_yldCPrime(3,1) = -0.6396161030250869D0
      g_yldCPrime(3,2) = -1.416739220885185D0
      g_yldCPrime(4,4) = 1.0767566538602826D0
      g_yldCPrime(5,5) = 1.069D0
      g_yldCPrime(6,6) = 1.0237D0

      g_yldCDbPrime(1,2) = -0.936165230037105D0
      g_yldCDbPrime(1,3) = -0.7827698684119997D0
      g_yldCDbPrime(2,1) = -0.6837450696952991D0
      g_yldCDbPrime(2,3) = -0.7782197960468684D0
      g_yldCDbPrime(3,1) = -1.1210834461178005D0
      g_yldCDbPrime(3,2) = 0.0930311520854942D0
      g_yldCDbPrime(4,4) = 1.1682692085440296D0
      g_yldCDbPrime(5,5) = 1.1471D0
      g_yldCDbPrime(6,6) = 1.0516D0

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
        eqStress = calc_eqStress(YLDM,f_yldCPrime,f_yldCDbPrime,STRESS)
        flowStress = calc_FlowStress(HARDK,HARDSTRAIN0,HARDN,eqpStrain)
        F = eqStress - flowStress
        IF (ISNAN(F)) THEN
          WRITE(7,*) "F is NaN"
          CALL XIT
        END IF

        ! if not yield
        IF (F < flowStress*TOLER) THEN
          CALL updateSTATEV(NSTATV,STATEV,eStrain,pStrain,eqpStrain)
          IF (iterationNum /= 0) THEN
            CALL updateDDSDDE(YLDM,f_yldCPrime,f_yldCDbPrime,
     &       g_yldCPrime,g_yldCDbPrime,HARDK,HARDSTRAIN0,HARDN,STRESS,
     &       eqpStrain,DDSDDE)
          END IF
          RETURN
        END IF

        ! if yield
        eqGStress = calc_eqGStress(YLDM,g_yldCPrime,g_yldCDbPrime,STRESS)
        H = calc_H(HARDK,HARDSTRAIN0,HARDN,eqpStrain)
        CALL calc_dfds(YLDM,f_yldCPrime,f_yldCDbPrime,STRESS,dfds)
        CALL calc_dgds(YLDM,g_yldCPrime,g_yldCDbPrime,STRESS,dgds)

        dLambda = F/(DOT_PRODUCT(dfds,MATMUL(DDSDDE,dgds)) + 
     &   H*eqGStress/eqStress)
        IF (dLambda<0) THEN
          WRITE(7,*) "dLambda < 0, invalid!!"
          CALL XIT
        END IF

        lambda = lambda + dLambda

        STRESS(:) = STRESS(:) - dLambda*MATMUL(DDSDDE,dgds)
        eqpStrain = eqpStrain + dLambda*eqGStress/eqStress
        pStrain(:) = pStrain(:) + dLambda*dgds(:)
        eStrain(:) = eStrain(:) - dLambda*dgds(:)
        iterationNum = iterationNum + 1
      END DO
      WRITE(7,*) "not converged!!!"
      WRITE(7,*) "F",F
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
      SUBROUTINE updateDDSDDE(YLDM,f_yldCPrime,f_yldCDbPrime,
     & g_yldCPrime,g_yldCDbPrime,HARDK,HARDSTRAIN0,HARDN,STRESS,
     & eqpStrain,DDSDDE)
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION eqStress,eqGStress,calc_eqStress,calc_eqGStress,
     & H,calc_H,dfds(6),dgds(6),ddgdds(6,6),subDDSDDE(6,6),DDSDDE(6,6),
     & invDDSDDE(6,6),A(7,7),B(6,6),invB(6,6),C(7,7),lambda,h0,
     & subVec(7),vec1(7),vec2(7),STRESS(6),eqpStrain
        
      DOUBLE PRECISION YLDM,f_yldCPrime(6,6),f_yldCDbPrime(6,6),
     & g_yldCPrime(6,6),g_yldCDbPrime(6,6),HARDK,HARDSTRAIN0,HARDN

      eqStress = calc_eqStress(YLDM,f_yldCPrime,f_yldCDbPrime,STRESS)
      eqGStress = calc_eqGStress(YLDM,g_yldCPrime,g_yldCDbPrime,STRESS)
      H = calc_H(HARDK,HARDSTRAIN0,HARDN,eqpStrain)
      CALL calc_dfds(YLDM,f_yldCPrime,f_yldCDbPrime,STRESS,dfds)
      CALL calc_dgds(YLDM,g_yldCPrime,g_yldCDbPrime,STRESS,dgds)
      CALL calc_ddgdds(YLDM,g_yldCPrime,g_yldCDbPrime,STRESS,dgds,ddgdds)

      subDDSDDE(:,:) = DDSDDE(:,:)
      CALL calc_Inverse(6,DDSDDE,invDDSDDE)
      invB(:,:) = invDDSDDE(:,:) + lambda*ddgdds(:,:)
      CALL calc_Inverse(6,invB,B)
      A(:,:) = 0.0D0
      A(1:6,1:6) = B(:,:)
      A(7,7) = 1.0D0
      h0 = eqGStress/eqStress
      subVec(:) = 0.0D0
      subVec(1:6) = dgds(:)
      subVec(7) = -1.0D0*h0
      vec1 = MATMUL(A,subVec)
      subVec(1:6) = dfds(:)
      subVec(7) = 0.0D0
      vec2 = MATMUL(subVec,A)
      CALL calc_Dyad(7,vec1,vec2,C)
      C(:,:) = C(:,:)/(DOT_PRODUCT(dfds,MATMUL(B,dgds)) + H*h0)
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


      ! calculate equivalent stress using yield function
      DOUBLE PRECISION FUNCTION calc_eqStress(YLDM,f_yldCPrime,
     & f_yldCDbPrime,STRESS)
      IMPLICIT NONE
      DOUBLE PRECISION YLDM,f_yldCPrime(6,6),f_yldCDbPrime(6,6),
     & STRESS(6),calc_yld2004_eqStress
      calc_eqStress = calc_yld2004_eqStress(YLDM,f_yldCPrime,
     & f_yldCDbPrime,STRESS)
      RETURN
      END FUNCTION calc_eqStress


      ! calculate equivalent stress using plastic potential
      DOUBLE PRECISION FUNCTION calc_eqGStress(YLDM,g_yldCPrime,
     & g_yldCDbPrime,STRESS)
      IMPLICIT NONE
      DOUBLE PRECISION YLDM,g_yldCPrime(6,6),g_yldCDbPrime(6,6),
     & STRESS(6),calc_yld2004_eqStress
      calc_eqGStress = calc_yld2004_eqStress(YLDM,g_yldCPrime,
     & g_yldCDbPrime,STRESS)
      RETURN
      END FUNCTION calc_eqGStress


      ! calculate flow stress
      DOUBLE PRECISION FUNCTION calc_FlowStress(HARDK,HARDSTRAIN0,HARDN,
     & eqpStrain)
      IMPLICIT NONE
      DOUBLE PRECISION Swift,HARDK,HARDSTRAIN0,HARDN,eqpStrain
      calc_FlowStress = Swift(HARDK,HARDSTRAIN0,HARDN,eqpStrain)
      RETURN
      END FUNCTION calc_FlowStress


      ! calculate first order differential of flow stress
      DOUBLE PRECISION FUNCTION calc_H(HARDK,HARDSTRAIN0,HARDN,eqpStrain)
      IMPLICIT NONE
      DOUBLE PRECISION calc_Swift_differential,HARDK,HARDSTRAIN0,HARDN,
     & eqpStrain
      calc_H = calc_Swift_differential(HARDK,HARDSTRAIN0,HARDN,eqpStrain)
      RETURN
      END FUNCTION calc_H


      ! calculate first order differentials of yield function
      SUBROUTINE calc_dfds(YLDM,f_yldCPrime,f_yldCDbPrime,STRESS,dfds)
      IMPLICIT NONE
      DOUBLE PRECISION dfds(6)
      DOUBLE PRECISION YLDM,f_yldCPrime(6,6),f_yldCDbPrime(6,6),
     & STRESS(6)
      CALL calc_yld2004_1st_differential(YLDM,f_yldCPrime,f_yldCDbPrime,
     & STRESS,dfds)
      RETURN
      END SUBROUTINE calc_dfds


      ! calculate first order differentials of plastic potential
      SUBROUTINE calc_dgds(YLDM,g_yldCPrime,g_yldCDbPrime,STRESS,dgds)
      IMPLICIT NONE
      DOUBLE PRECISION dgds(6)
      DOUBLE PRECISION YLDM,g_yldCPrime(6,6),g_yldCDbPrime(6,6),
     & STRESS(6)
      CALL calc_yld2004_1st_differential(YLDM,g_yldCPrime,g_yldCDbPrime,
     & STRESS,dgds)
      RETURN
      END SUBROUTINE calc_dgds


      ! calculate second order differentials of plastic potential
      SUBROUTINE calc_ddgdds(YLDM,g_yldCPrime,g_yldCDbPrime,STRESS,dgds,
     & ddgdds)
      IMPLICIT NONE
      DOUBLE PRECISION ddgdds(6,6)
      DOUBLE PRECISION YLDM,g_yldCPrime(6,6),g_yldCDbPrime(6,6),
     & STRESS(6),dgds(6)
      CALL calc_yld2004_2nd_differential(YLDM,g_yldCPrime,g_yldCDbPrime,
     & STRESS,dgds,ddgdds)
      RETURN
      END SUBROUTINE calc_ddgdds


      ! calculate principal value of linealy transformed matrix
      SUBROUTINE calc_yld2004_Principal(yldCPrime,yldCDbPrime,
     & yldSPriPrime,yldSPriDbPrime,STRESS)
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
      END SUBROUTINE calc_yld2004_Principal


      ! calculate equivalent stress using Yld2004-18P
      DOUBLE PRECISION FUNCTION calc_yld2004_eqStress(YLDM,yldCPrime,
     & yldCDbPrime,STRESS)
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION YLDM,yldCPrime(6,6),yldCDbPrime(6,6),STRESS(6),
     & yldSPriPrime(3),yldSPriDbPrime(3),yldPhi
      yldSPriPrime(:) = 0.0D0
      yldSPriDbPrime(:) = 0.0D0
      CALL calc_yld2004_Principal(yldCPrime,yldCDbPrime,yldSPriPrime,
     & yldSPriDbPrime,STRESS)
      yldPhi = 0.0D0
      DO i=1,3
        DO j=1,3
          yldPhi = yldPhi + 
     &     ABS(yldSPriPrime(i) - yldSPriDbPrime(j))**YLDM
        END DO
      END DO
      calc_yld2004_eqStress = (yldPhi/4.0D0)**(1.0D0/YLDM)
      RETURN
      END FUNCTION calc_yld2004_eqStress


      ! calculate first order differential of Yld2004-18P
      SUBROUTINE calc_yld2004_1st_differential(YLDM,yldCPrime,
     & yldCDbPrime,STRESS,firstDiff)
      IMPLICIT NONE
      INTEGER i,j,k
      DOUBLE PRECISION YLDM,yldCPrime(6,6),yldCDbPrime(6,6),STRESS(6),
     & firstDiff(6),DELTAX,yldSPriPrime(3),yldSPriDbPrime(3),yldPhi,
     & subStress(6),subyldPhi,dPhidS(6),eqStress,calc_yld2004_eqStress,
     & dfdPhi
      PARAMETER(DELTAX=1.0D-5)
      yldSPriPrime(:) = 0.0D0
      yldSPriDbPrime(:) = 0.0D0
      CALL calc_yld2004_Principal(yldCPrime,yldCDbPrime,yldSPriPrime,
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
        CALL calc_yld2004_Principal(yldCPrime,yldCDbPrime,yldSPriPrime,
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
      eqStress = calc_yld2004_eqStress(YLDM,yldCPrime,yldCDbPrime,STRESS)
      dfdPhi = (eqStress**(1-YLDM))/(4.0D0*YLDM)
      firstDiff(:) = dfdPhi*dPhidS(:)
      RETURN
      END SUBROUTINE calc_yld2004_1st_differential


      ! calculate second order differential of Yld2004-18P
      SUBROUTINE calc_yld2004_2nd_differential(YLDM,yldCPrime,
     & yldCDbPrime,STRESS,firstDiff,secondDiff)
      IMPLICIT NONE
      INTEGER i
      DOUBLE PRECISION YLDM,yldCPrime(6,6),yldCDbPrime(6,6),STRESS(6),
     & firstDiff(6),secondDiff(6,6),DELTAX,subStress(6),subfirstDiff(6)
      PARAMETER(DELTAX=1.0D-1)
      DO i=1,6
        subStress(:) = STRESS(:)
        subStress(i) = subStress(i) + DELTAX
        CALL calc_yld2004_1st_differential(YLDM,yldCPrime,yldCDbPrime,
     &   subStress,subfirstDiff)
        secondDiff(:,i) = (subfirstDiff(:) - firstDiff(:))/DELTAX
      END DO  
      RETURN
      END SUBROUTINE calc_yld2004_2nd_differential


      ! calculate flow stress
      DOUBLE PRECISION FUNCTION Swift(HARDK,HARDSTRAIN0,
     & HARDN,eqpStrain)
      IMPLICIT NONE
      DOUBLE PRECISION HARDK,HARDSTRAIN0,HARDN,eqpStrain
      Swift = HARDK*(HARDSTRAIN0 + eqpStrain)**HARDN
      RETURN
      END FUNCTION Swift


      ! calculate first order differential of Kim-Tuan flow stress
      DOUBLE PRECISION FUNCTION calc_Swift_differential(HARDK,
     & HARDSTRAIN0,HARDN,eqpStrain)
      IMPLICIT NONE
      DOUBLE PRECISION HARDK,HARDSTRAIN0,HARDN,eqpStrain
      calc_Swift_differential = HARDK*HARDN*
     & ((HARDSTRAIN0 + eqpStrain)**(HARDN-1.0D0))
      RETURN
      END FUNCTION calc_Swift_differential
