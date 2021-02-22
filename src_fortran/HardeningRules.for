      !----------------------------------------
      ! Kim-Tuan
      !----------------------------------------
      
      ! calculate flow stress using Kim-Tuan rule
      DOUBLE PRECISION FUNCTION KimTuan(HARDK,HARDT,HARDA,HARDH,
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
      END FUNCTION KimTuan


      ! calculate first order differential of Kim-Tuan flow stress
      DOUBLE PRECISION FUNCTION calc_KimTuan_differential(HARDK,HARDT,
     & HARDA,HARDH,HARDSTRESS0,eqpStrain)
      IMPLICIT NONE
      DOUBLE PRECISION HARDK,HARDT,HARDA,HARDH,HARDSTRESS0,eqpStrain
      IF (eqpStrain == 0.0D0) THEN
        calc_KimTuan_differential = 0.0D0
        RETURN
      END IF
      calc_KimTuan_differential = 
     & (HARDK*HARDA*HARDT*eqpStrain**(HARDA-1.0D0))*
     & EXP(-1.0D0*HARDT*eqpStrain**HARDA)
      calc_KimTuan_differential = calc_KimTuan_differential + 
     & HARDK*(1.0D0-EXP(-1.0D0*HARDT*eqpStrain**HARDA))*
     & HARDH*eqpStrain**(HARDH-1.0D0)
      RETURN
      END FUNCTION calc_KimTuan_differential


      !----------------------------------------
      ! Swift
      !----------------------------------------

      ! calculate flow stress
      DOUBLE PRECISION FUNCTION Swift(HARDK,HARDSTRAIN0,
     & HARDN,eqpStrain)
      IMPLICIT NONE
      DOUBLE PRECISION HARDK,HARDSTRAIN0,HARDN,eqpStrain
      calc_FlowStress = HARDK*(HARDSTRAIN0 + eqpStrain)**HARDN
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