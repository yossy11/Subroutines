      !----------------------------------------
      ! YLD2004-18P
      !----------------------------------------
      
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
     &   STRESS,subfirstDiff)
        secondDiff(:,i) = (subfirstDiff(:) - firstDiff(:))/DELTAX
      END DO  
      RETURN
      END SUBROUTINE calc_yld2004_2nd_differential


      !----------------------------------------
      ! Hill48
      !----------------------------------------

      ! calculate equivalent stress using Hill48
      DOUBLE PRECISION FUNCTION calc_hill48_eqStress(hillParams,STRESS)
      IMPLICIT NONE
      DOUBLE PRECISION hillParams(4),STRESS(6),numerator
      numerator = hillParams(1)*(STRESS(2) - STRESS(3))**2 + 
     & hillParams(2)*(STRESS(3)-STRESS(1))**2 + 
     & hillParams(3)*(STRESS(1)-STRESS(2))**2 + 
     & 2.0D0*hillParams(4)*SUM(STRESS(4:6)**2)
      calc_eqStress = SQRT(1.5D0*numerator/SUM(hillParams(1:3)))
      RETURN
      END FUNCTION calc_hill48_eqStress


      ! calculate differential of plastic potential
      SUBROUTINE calc_hill48_1st_differential(hillParams,STRESS,
     & firstDiff)
      IMPLICIT NONE
      DOUBLE PRECISION hillParams(4),STRESS(6),firstDiff(6),eqStress,
     & calc_hill48_eqStress,multiplier
      eqStress = calc_hill48_eqStress(hillParams,STRESS)
      multiplier = 0.75/(SUM(hillParams(1:3))*eqStress)
      firstDiff(1) = -1.0D0*2*hillParams(2)*(STRESS(3)-STRESS(1))+
     & 2*hillParams(3)*(STRESS(1)-STRESS(2))
      firstDiff(2) = 2*hillParams(1)*(STRESS(2)-STRESS(3))-
     & 2*hillParams(3)*(STRESS(1)-STRESS(2))
      firstDiff(3) = -1.0D0*2*hillParams(1)*(STRESS(2)-STRESS(3))+
     & 2*hillParams(2)*(STRESS(3)-STRESS(1))
      firstDiff(4) = 4*hillParams(4)*STRESS(4)
      firstDiff(5) = 4*hillParams(4)*STRESS(5)
      firstDiff(6) = 4*hillParams(4)*STRESS(6)
      firstDiff(:) = firstDiff(:) * multiplier
      RETURN
      END SUBROUTINE calc_hill48_1st_differential


      ! calculate differential of plastic potential
      SUBROUTINE calc_hill48_2nd_differential(hillParams,STRESS,
     & firstDiff,secondDiff)
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION hillParams(4),STRESS(6),firstDiff(6),
     & secondDiff(6,6),eqStress,calc_hill48_eqStress,multiplier
      eqStress = calc_hill48_eqStress(hillParams,STRESS)
      multiplier = 0.75/(SUM(hillParams(1:3))*eqStress)

      secondDiff(:,:) = 0.0D0
      secondDiff(1,1) = 2*hillParams(2) + 2*hillParams(3)
      secondDiff(2,2) = 2*hillParams(1) + 2*hillParams(3)
      secondDiff(3,3) = 2*hillParams(1) + 2*hillParams(2)
      secondDiff(4,4) = 4*hillParams(4)
      secondDiff(5,5) = 4*hillParams(4)
      secondDiff(6,6) = 4*hillParams(4)
      secondDiff(1,2) = -1.0D0*2*hillParams(3)
      secondDiff(1,3) = -1.0D0*2*hillParams(2)
      secondDiff(2,1) = -1.0D0*2*hillParams(3)
      secondDiff(2,3) = -1.0D0*2*hillParams(1)
      secondDiff(3,1) = -1.0D0*2*hillParams(2)
      secondDiff(3,2) = -1.0D0*2*hillParams(1)
      DO i=1,6
        DO j=1,6
          secondDiff(i,j) = multiplier*secondDiff(i,j) - 
     &     firstDiff(i)*firstDiff(j)/eqStress
        END DO
      END DO
      RETURN
      END SUBROUTINE calc_hill48_2nd_differential