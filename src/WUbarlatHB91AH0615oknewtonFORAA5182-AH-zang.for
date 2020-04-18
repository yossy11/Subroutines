C***********************************************************************
C***UMAT FOR ABAQUS/STANDARD,LARGE DEFORMATION FORMULATION        ******
C***                                                              ******
C***WRITTEN BY WU BOXUN     THE UNIVERSITY OF TOKYO               ******
C***NON-ASSOCIATED FLOW RULES   &YLD91  HILL 1948                 ******
C*** ALGORITHM   BACKWARD EULER       RETURN MAPPING              ******
C***YIELD FUNCTION F     POTENTIAL FUNCTION G                     ******
C***********************************************************************
C
*******USER SUBROUTINE 
C
        SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,
        DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,
        NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,
        DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
       INCLUDE 'ABA_PARAM.INC'
C
       CHARACTER*80 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C   DEFINE THE VARIABLES
       PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=5.D-1,
     1 ENUMAX=0.4999D0,TOLER=1.0D-6,FOUR=4.D0,SIX=6.D0,NINE=9.D0,
     2 NEWTON=80,ONEHALF=1.5D0,PI=180)	 

       DOUBLE PRECISION EMOD,ENU,EG,EG2,EG3,EBULK3,ELAME,YIELD0,YIELDD
C   DEFINE THE INTERMEDIATE VARIABLES
C   HH-HARDENING RATE  HW-FACTOR
       DOUBLE PRECISION SIGMA_X,SIGMA_Y,SIGMA_Z,TAO_XY,TAO_YZ,TAO_ZX,ALPHAF,ALPHAG,NBLFM1,NBLFM2,DNBL,NBLFM,NBLFZ,CNBL,JCOBFM,BBB,GAA
C	
       DOUBLE PRECISION EQPSTRESS
       DOUBLE PRECISION FA0,FB0,FC0,FA90,FB90,FC90,FA45,FB45,FC45,
     1 FAB,FBB,FCB,SIGMA0,SIGMA90,SIGMA45,SIGMAB,HARD0,HARD90,HARD45,
     2 HARDB,TH1,TH2,TH3,TH4
       DOUBLE PRECISION SIGMA0K,SIGMA90K,SIGMA45K,SIGMABK,HARD0K,
     1 HARD90K,HARD45K,HARDBK,TH1K,TH2K,TH3K,TH4K,FAK,FBK,FCK,FHK,
     2 FGK,FFK        
C   VARIABLES IN F FUNCTION
       DOUBLE PRECISION FSIGMA_SX,FSIGMA_SY,FSIGMA_SZ,
	1 FTAO_SYZ,FTAO_SZX,FTAO_SXY,FEQSTRESS,FSUM,
	2 FSA,FSB,FSC,FSD,FSE,FSF
C   VARIABLES IN G FUNCTION
       DOUBLE PRECISION GSIGMA_SX,GSIGMA_SY,GSIGMA_SZ,
	1 GTAO_SYZ,GTAO_SZX,GTAO_SXY,GEQSTRESS,GSUM,GH,GF,GG,GL,GM,GN,
	2 GSA,GSB,GSC,GSD,GSE,GSF
C	DEFINE STATE VARIABLES
       DOUBLE PRECISION EQPLAS
C     DEFINE UPDATA VARIABLES	   
       DOUBLE PRECISION DEQPLAS		
C     NEW VARIABLES (UN-FINISHED)
       DIMENSION EELAS(6),EPLAS(6),PSTRES(6),DDFDDS(6),DDGDDS(6),
     1 JCOBFZ(6,6),DDBDDS(6,6),DSTRES(6),DDGDS(6),DDSTIFF(6,6),
	2 DDEDDB(6,6),Q(6,6),INVQ(6,6),R(6,6),NBLFMT(6), 
     3 R0(6),NBLQR(6),ATRT(6),RB(6),IS(6),JS(6),PHI(20),
     4 TSTRESS(6),KSTRESS(6),KDSTRES(6),MSTRESS(6),NF(6),NG(6),
     5 NDDGDS(6)
       DOUBLE PRECISION EELAS,EPLAS,PSTRES,DDFDDS,DDGDDS,
     1 JCOBFZ,DDBDDS,DSTRES,DDGDS,DDSTIFF,
	2 DDEDDB,Q,INVQ,R,NBLFMT,
     3 R0,NBLQR,ATRT,RB,IS,JS,TSTRESS,PHI,KSTRESS,KDSTRES,
     4 MSTRESS,NG,NF,FB,FC,FK,HARDK,KDNBL,NDNBL,IDNBL,TDNBL
C***FLD****FLD**************FLD************************
      DOUBLE PRECISION SMAX,SMIN,FLDF,ETOTAL,FLDW
      DIMENSION ETOTAL(6)
       
       INTEGER NVALUE,L
C***
C    YLD91 F FUNCTION 
       DIMENSION FPHI_S(3),FS_I(3,3),FI_L(3,6),FL_SIGMA(6,6),
     1 FM36X66(3,6),FM33X36(3,6),FPHI_SIGMA(6)
       DIMENSION FS_LAMDA(3),FI_LAMDA(3),FL_LAMDA(6),FL_P(6,6),
     2  FP_LAMDA(6)
       DOUBLE PRECISION FSIGMA_LAMDA,FPHI_LAMDA,FS_LAMDA,
     1 FI_LAMDA,FL_LAMDA,FL_P,FP_LAMDA
       DIMENSION GPHI_S(3),GS_I(3,3),GI_L(3,6),GL_SIGMA(6,6),
     1 GM36X66(3,6),GM33X36(3,6),GPHI_SIGMA(6)
       DIMENSION GPHI_SS(3,3),GS1_II(3,3),GS2_II(3,3),
     1 GS3_II(3,3),GI2_LL(6,6),GI3_LL(6,6)
       DIMENSION GTEMP1(6,6),GTEMP2(6,6),GTEMP3(6,6) 
       DOUBLE PRECISION FPHI_S,FS_I,FI_L,FL_SIGMA,FM36X66,FM33X36,
     1 GPHI_S,GS_I,GI_L,GL_SIGMA,GM36X66,GM33X36,GPHI_SS,GS1_II,GS2_II,
     2 GS3_II,GI2_LL,GI3_LL,GTEMP1,GTEMP2,GTEMP3,FPHI_SIGMA,GPHI_SIGMA
       DOUBLE PRECISION SIGMA_SA,SIGMA_SB,SIGMA_SC,SIGMA_SF,SIGMA_SG,
     1 SIGMA_SH 
       DOUBLE PRECISION FYA,FYB,FYC,FYG,FYH,FYF,GYA,GYB,GYC GYG GYH GYF,
     1 FYAK,FYBK,FYCK,FYGK,FYHK,FYFK,
     1 FYM,GYM,FI1,FI2,FI3,GI1,GI2,GI3,FS1,FS2,FS3,GS1,GS2,GS3,
     2 FL1,FL2,FL3,FL4,FL5,FL6,GL1,GL2,GL3,GL4,GL5,GL6,
     3 FRHO,FTHETA,FRHO3,GRHO,GTHETA,GRHO3
       DOUBLE PRECISION AAH,BBH,CCH,HHH,GGH,FFH
       DOUBLE PRECISION AAHK,BBHK,CCHK,HHHK,GGHK,FFHK
       DIMENSION x(4),w(4)
       DOUBLE PRECISION x,w
C***
C*************BACKSTRESS********************
C
      DOUBLE PRECISION CI,CII,C,V,EQPLASOLD,BI,BII,BIII,DHIDS,F,
     1 CNBLFM,CNBLFZ,MAFZ,MAFM,RA,RBI,RBII
C      
C      
      DIMENSION AA(6),BBI(6),BBII(6),AI(6),AII(6),A(6),DAI(6),
     1 DAII(6),DA(6),A11(6,6),A13(6,6),A14(6,6),A31(6,6),A31s(6,6),
     2 A33(6,6),A33s(6,6),A34(6,6),AIOLD(6),AIIOLD(6),HI(6),
     3 HII(6),A41(6,6),A43(6,6),A44(6,6),A44s(6,6), 
     4 INVDDSTIFF(6,6),BI(6,6),INVXX(6,6),XX(6,6),LEFT(6),RIGHT(6),
     5 A11AA(6),A13BBI(6),A14BBII(6),A31AA(6),A33BBI(6),
     6 A34BBII(6),A41AA(6),A43BBI(6),A44BBII(6),AFZ(6),
     7 A11DDGDDS(6),A13HI(6),A14HII(6),A31DDGDDS(6),A33HI(6),
     8 A34HII(6),A41DDGDDS(6),A43HI(6),A44HII(6),AFM(6),
     9 PBSTRES(6),EPLASOLD(6),KBSTRES(6),TBSTRES(6)
C      
C      
      DOUBLE PRECISION AA,BBI,BBII,AI,AII,A,DAI,
     1 DAII,DA,A11,A13,A14,A31,
     2 A33,A33s,A34,AIOLD,AIIOLD,HI,
     3 HII,A41,A43,A44,A44s,
     4 INVDDSTIFF,INVXX,XX,LEFT,RIGHT,
     5 A11AA,A13BBI,A14BBII,A31AA,A33BBI,
     6 A34BBII,A41AA,A43BBI,A44BBII,AFZ,
     7 A11DDGDDS,A13HI,A14HII,A31DDGDDS,A33HI,
     8 A34HII,A41DDGDDS,A43HI,A44HII,AFM,PBSTRES,KBSTRES,
     9 EPLASOLD,XDNBL      	   
c
       
C***********************************************************
C   
      EMOD=PROPS(1)
	ENU=MIN(PROPS(2),ENUMAX)
      FYM=PROPS(3)
C      
       FA0=PROPS(4)
	FB0=PROPS(5)
      FC0=PROPS(6)
c      FD0=PROPS(7)
      
      FA45=PROPS(8)
      FB45=PROPS(9)
	FC45=PROPS(10)
      FD45=PROPS(11)

      FA90=PROPS(12)
      FB90=PROPS(13)
      FC90=PROPS(14)
      FD90=PROPS(15)

      FAB=PROPS(16)
      FBB=PROPS(17)
      FCB=PROPS(18)
      FDB=PROPS(19)
C
	GH=PROPS(20)
	GF=PROPS(21)
	GG=PROPS(22)
	GL=PROPS(23)
	GM=PROPS(23)
	GN=PROPS(23)
C      

      CI=PROPS(24)
      CII=PROPS(25)
      V=PROPS(26)
C
      C=CII/CI
C	  
      EG=EMOD/(TWO*(ONE+ENU))
C   3K	  
	EBULK3=EMOD/(ONE-TWO*ENU)  
C   2G 3G
	EG2=TWO*EG
	EG3=THREE*EG
      ELAME=(EBULK3-EG2)/THREE	
C     INITIALIZE THE MATRIX(NOT FINISHED)
      DO K1=1,NTENS
         DO K2=1,NTENS
            DDSDDE(K1,K2)=ZERO
            JCOBFZ(K1,K2)=ZERO	
            DDBDDS(K1,K2)=ZERO 			
			DDSTIFF(K1,K2)=ZERO
			Q(K1,K2)=ZERO
			INVQ(K1,K2)=ZERO
		    R(K1,K2)=ZERO
         ENDDO
	  PSTRES(K1)=ZERO 
	  PBSTRES(K1)=ZERO
	  DSTRES(K1)=ZERO
C	  
        AI(K1)=ZERO
        AII(K1)=ZERO
        A(K1)=ZERO
        DAI(K1)=ZERO
        DAII(K1)=ZERO
        DA(K1)=ZERO
C		
        TSTRESS(K1)=ZERO
        TBSTRES(K1)=ZERO
        KSTRESS(K1)=ZERO
        KBSTRES(K1)=ZERO
        KDSTRES(K1)=ZERO
	  DDGDS(K1)=ZERO
        NBLFMT(K1)=ZERO
        NBLQR(K1)=ZERO
        R0(K1)=ZERO
        ATRT(K1)=ZERO
        RB(K1)=ZERO
      ENDDO
C   ELASTIC MATRIX
      DO K1=1,NDI
        DO K2=1,NDI
        DDSDDE(K1,K2)=ELAME
        ENDDO
        DDSDDE(K1,K1)=EG2+ELAME
      ENDDO
      DO K1=NDI+1,NTENS
        DDSDDE(K1,K1)=EG
      ENDDO
C   SAVE THE ELASTIC MATRIX TO DDSTIFF
      DO K1=1,NTENS
         DO K2=1,NTENS
		 DDSTIFF(K1,K2)=DDSDDE(K1,K2)
	   ENDDO
      ENDDO	 
C	READ THE EQUIVALENT PLASTIC STRAIN  
      CALL ROTSIG(STATEV(1),DROT,EELAS,2,NDI,NSHR)
      CALL ROTSIG(STATEV(NTENS+1),DROT,EPLAS,2,NDI,NSHR)
C
      DO K1=1,NTENS
      EELAS(K1)=STATEV(K1)
      EPLAS(K1)=STATEV(K1+NTENS)
C      
      ETOTAL(K1)=STATEV(K1+40)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         AI(K1)=STATEV(K1+30)
         AII(K1)=STATEV(K1+36)
         A(K1)=AI(K1)+AII(K1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	  
      ENDDO
C     READ EQPLAS    INITIAL VALUE IS 0
      EQPLAS=STATEV(13) 	  
C      WRITE(7,*),'STRESS(1)',STRESS(1)
C
C*****************************************************
C  CALCULATE THE TRIAL STRESS--PSTRES  // ON B POINT
C*****************************************************
       DO K1=1,NTENS
	   DO K2=1,NTENS
		 DSTRES(K1)=DSTRES(K1)+DDSDDE(K1,K2)*DSTRAN(K2)
         ENDDO
         EELAS(K1)=EELAS(K1)+DSTRAN(K1)
      ENDDO 
C      
C      WRITE(6,*),'EELAS',EELAS(1)
C      
      DO K1=1,NTENS
          PSTRES(K1)=STRESS(K1)+DSTRES(K1)
		  PBSTRES(K1)=PSTRES(K1)-A(K1)
      ENDDO
C   READ STRESS VARIABLES
      SIGMA_X=PBSTRES(1)
	SIGMA_Y=PBSTRES(2)
      SIGMA_Z=PBSTRES(3)
      TAO_XY=PBSTRES(4)
	TAO_ZX=PBSTRES(5)
      TAO_YZ=PBSTRES(6)

	GSUM=ONE/(GF+GG+GH)      
C      
C     CALCULATE THE WORK-HARDENING CURVE
      
       
      SIGMA0=FA0+FB0*(1-exp(-FC0*EQPLAS))-CI*(1-exp(-V*EQPLAS))/V
      SIGMA90=FA90-FB90*EXP(-FC90*EQPLAS)+FD90*EQPLAS      
      SIGMA45=FA45-FB45*EXP(-FC45*EQPLAS)+FD45*EQPLAS      
      SIGMAB=FAB-FBB*EXP(-FCB*EQPLAS)+FDB*EQPLAS
C      
C     
C      WRITE(6,*),SIGMA0,SIGMA45,SIGMA90,SIGMAB,'SIGMA'
      IF(EQPLAS.EQ.0.0)THEN
      HARD0=EMOD
      HARD90=EMOD
      HARD45=EMOD
      HARDB=EMOD
      ELSE

      HARD0=FB0*FC0*exp(-FC0*EQPLAS)-CI*exp(-V*EQPLAS)          
      HARD45=FC45*FB45*EXP(-FC45*EQPLAS)+FD45
      HARD90=FC90*FB90*EXP(-FC90*EQPLAS)+FD90
      HARDB=FCB*FBB*EXP(-FCB*EQPLAS)+FDB  
      ENDIF
C      
C	DEFINE THE DERIVATION STRESS IN F FUNCTION
      SIGMA_SA=SIGMA_Y-SIGMA_Z
      SIGMA_SB=SIGMA_Z-SIGMA_X
      SIGMA_SC=SIGMA_X-SIGMA_Y
      SIGMA_SH=TAO_XY
      SIGMA_SG=TAO_ZX
      SIGMA_SF=TAO_YZ
c*************************************************************
      CALL PARACAL(SIGMA0,SIGMA90,SIGMAB,SIGMA45,FYM,x)
      FYA=x(1)
      FYB=x(2)
      FYC=x(3)
      FYH=x(4)
      FYG=x(4)
      FYF=x(4)
C      write(7,*),'a,b,c,h,M',FYA,FYB,FYC,FYH,FYM
      CALL PARACALAH(FYM,w,SIGMA0,SIGMA90,SIGMAB,SIGMA45,
     1 FYA,FYB,FYC,FYH,HARD0,HARD90,HARDB,HARD45)
c      write(6,*),PSTRES
C      WRITE(6,*),'FYM,x,w',FYM,x,w
c       write(6,*),x
C     
      AAH=w(1)
      BBH=w(2)
      CCH=w(3)
      HHH=w(4)
      GGH=w(4)
      FFH=w(4)
      DO K1=1,NTENS
          FP_LAMDA(K1)=ZERO
      ENDDO
      FP_LAMDA(1)=AAH
      FP_LAMDA(2)=BBH
      FP_LAMDA(3)=CCH
      FP_LAMDA(4)=HHH
      FP_LAMDA(5)=GGH
      FP_LAMDA(6)=FFH
c     
C      DO K1=1,NTENS
C          WRITE(6,*),'FP_LAMDA(K)',FP_LAMDA(K1)
C      ENDDO
C      
      FL1=(FYC*SIGMA_SC-FYB*SIGMA_SB)/THREE
      FL2=(FYA*SIGMA_SA-FYC*SIGMA_SC)/THREE
      FL3=(FYB*SIGMA_SB-FYA*SIGMA_SA)/THREE
      FL4=FYH*SIGMA_SH
      FL5=FYG*SIGMA_SG
      FL6=FYF*SIGMA_SF
C      WRITE(6,*),'L1-L6',FL1,FL2,FL3,FL4
C     
      FI1=ZERO
      FI2=(FL4**2+FL5**2+FL6**2)/3+(FL1**2+FL2**2+FL3**2)/6
      FI3=(FL1*FL2*FL3+2*FL4*FL5*FL6-FL3*FL4**2-FL2*FL5**2-FL1*FL6**2)/2
C     
C      WRITE(6,*),'I2',FI2,'I3',FI3 
      FRHO=FI2**(ONEHALF)
      FTHETA=ACOSD(FI3/(FI2**ONEHALF))
      FRHO3=FI2**(-0.5)
C      FRHO3=FRHO**(-1/3)
C      WRITE(6,*),'theta',FTHETA,'RHO',FRHO
C
      FS1=TWO*SQRT(FI2)*COSD(FTHETA/THREE)
      FS2=TWO*SQRT(FI2)*COSD((2*PI-FTHETA)/THREE)
      FS3=TWO*SQRT(FI2)*COSD((2*PI+FTHETA)/THREE)
C
C ******     
      FPHI=ABS(FS1-FS2)**FYM+ABS(FS2-FS3)**FYM+ABS(FS3-FS1)**FYM      
      FEQSTRESS=(FPHI/2)**(1/FYM)
C ******	  
C      WRITE(6,*),'FEQSTRESS',FEQSTRESS
C     DEFINE THE DERIVATION STRESS IN G FUNCTION
	 GSIGMA_SX=(GH*(SIGMA_X-SIGMA_Y)+GG*(SIGMA_X-SIGMA_Z))*GSUM
	 GSIGMA_SY=(GF*(SIGMA_Y-SIGMA_Z)+GH*(SIGMA_Y-SIGMA_X))*GSUM
	 GSIGMA_SZ=(GG*(SIGMA_Z-SIGMA_X)+GF*(SIGMA_Z-SIGMA_Y))*GSUM
	 GTAO_SYZ=GL*TAO_YZ*GSUM
	 GTAO_SZX=GM*TAO_ZX*GSUM
	 GTAO_SXY=GN*TAO_XY*GSUM
C     OTHER STRESS IN G FUNCTION 	 
	 GEQSTRESS=GF*(SIGMA_Y-SIGMA_Z)**TWO+GG*(SIGMA_Z-SIGMA_X)**TWO
     1 +GH*(SIGMA_X-SIGMA_Y)**TWO+TWO*GL*TAO_YZ**TWO+
     2 TWO*GM*TAO_ZX**TWO+TWO*GN*TAO_XY**TWO
	 GEQSTRESS=DSQRT(THREE*GSUM*GEQSTRESS/TWO)       
C
C
C     CALCULATE THE YIELD STRESS AT B POINT
 
C  
C
       YIELDD=SIGMA0
c      ENDIF
C  
	FB=FEQSTRESS-YIELDD
C      WRITE(6,*),'FS1',FS1,'FS2',FS2,'FS3',FS3
C      WRITE(6,*),'FEQSTRESS1',FEQSTRESS,'GEQSTRESS1',GEQSTRESS
      
C,'FB',FB,'YIELDD',YIELDD,
C            WRITE(6,*),'EQPLAS',EQPLAS    
C	IF FB>0 PLASTIC  FB<0 ELASTIC
*********************************************      
      IF(FB.GT.TOLER)THEN
*********************************************       
C      WRITE(6,*),'FB--FEQSTRESS-yieldd',FB,FEQSTRESS,YIELDD
C	FIRST ORDER F-SIGMA 
C
C	
C          WRITE(6,*),'FB',FB 
       DO K1=1,20
         PHI(K1)=ZERO
        ENDDO
       PHI(1)=(1.0D-1)
       PHI(2)=(1.0D-2)
       PHI(3)=(1.0D-3)
       PHI(4)=(1.0D-4)
       PHI(5)=(1.0D-5)
C       
      DO K1=1,NTENS
         TSTRESS(K1)=PSTRES(K1)
         TBSTRES(K1)=PBSTRES(K1)
      END DO
C**************************
      IF ((FI2.EQ.ZERO).and.(FI3.EQ.ZERO)) THEN
      DO K1=1,NTENS
	   DDFDDS(K1)=ZERO
      ENDDO
C         
         FSIGMA_LAMDA=ZERO
C
      GOTO 1234
      ELSE
C**************************      
c      WRITE(6,*),'FI2',FI2,'FI3',FI3
c      WRITE(6,*),'FRHO',FRHO,'FTHETA',FTHETA
c      WRITE(6,*),'FEQSTRESS',FEQSTRESS,'GEQSTRESS',GEQSTRESS    
      DO K1=1,NDI
	  FPHI_S(K1)=ZERO
      ENDDO
        FPHI_S(1)=((FS1-FS2)*(ABS(FS1-FS2)**(FYM-2))-
     1 (FS3-FS1)*(ABS(FS3-FS1)**(FYM-2)))*FYM
        FPHI_S(2)=((FS2-FS3)*(ABS(FS2-FS3)**(FYM-2))-
     1 (FS1-FS2)*(ABS(FS1-FS2)**(FYM-2)))*FYM
        FPHI_S(3)=((FS3-FS1)*(ABS(FS3-FS1)**(FYM-2))-
     1 (FS2-FS3)*(ABS(FS2-FS3)**(FYM-2)))*FYM
C        write(6,*),'FS1',FS1,'FS2',FS2,'FS3',FS3
C
      DO K1=1,NDI
	     DO K2=1,NDI
	     FS_I(K1,K2)=ZERO
	     ENDDO
	  ENDDO
        FS_I(1,1)=0
        FS_I(1,2)=FRHO3*sind(2*FTHETA/3)/sind(FTHETA)
        FS_I(1,3)=2*FRHO3*FRHO3*sind(FTHETA/3)/(3*sind(FTHETA))
        FS_I(2,1)=0
        FS_I(2,2)=FRHO3*sind(2*(FTHETA+PI)/3)/sind(FTHETA)
        FS_I(2,3)=2*FRHO3*FRHO3*sind((FTHETA-2*PI)/3)/(3*sind(FTHETA))
        FS_I(3,1)=0
        FS_I(3,2)=FRHO3*sind(2*(FTHETA-PI)/3)/sind(FTHETA)
        FS_I(3,3)=2*FRHO3*FRHO3*sind((FTHETA+2*PI)/3)/(3*sind(FTHETA))
C     
      DO K1=1,NDI
	     DO K2=1,NTENS
	     FI_L(K1,K2)=ZERO
	     ENDDO
	  ENDDO 
        FI_L(2,1)=FL1/3
	  FI_L(2,2)=FL2/3
	  FI_L(2,3)=FL3/3
	  FI_L(2,4)=2*FL4/3
	  FI_L(2,5)=2*FL5/3
	  FI_L(2,6)=2*FL6/3
	  FI_L(3,1)=(FL2*FL3-FL6*FL6)/2
	  FI_L(3,2)=(FL1*FL3-FL5*FL5)/2
	  FI_L(3,3)=(FL1*FL2-FL4*FL4)/2
	  FI_L(3,4)=FL5*FL6-FL3*FL4
	  FI_L(3,5)=FL4*FL6-FL2*FL5
	  FI_L(3,6)=FL4*FL5-FL1*FL6
C	  
      DO K1=1,NTENS
	     DO K2=1,NTENS
	     FL_SIGMA(K1,K2)=ZERO
	     ENDDO
	  ENDDO 
	  FL_SIGMA(1,1)=(FYB+FYC)/3
	  FL_SIGMA(1,2)=-FYC/3
	  FL_SIGMA(1,3)=-FYB/3
	  FL_SIGMA(2,1)=-FYC/3
	  FL_SIGMA(2,2)=(FYA+FYC)/3
	  FL_SIGMA(2,3)=-FYA/3
	  FL_SIGMA(3,1)=-FYB/3
	  FL_SIGMA(3,2)=-FYA/3
	  FL_SIGMA(3,3)=(FYA+FYB)/3
	  FL_SIGMA(4,4)=FYH
	  FL_SIGMA(5,5)=FYG
	  FL_SIGMA(6,6)=FYF
C       
C	  
	  DO 177 K1=1,NDI
	  DO 177 K2=1,NTENS
	     FM36X66(K1,K2)=ZERO
	     DO 178 K3=1,NTENS
		 FM36X66(K1,K2)=FM36X66(K1,K2)+FI_L(K1,K3)*FL_SIGMA(K3,K2)
178        CONTINUE
C        WRITE(6,*),'FM36X66',FM36X66(K1,K2),K1,K2   
177     CONTINUE   
C        
C       
C      
	  DO 777 K1=1,NDI
	  DO 777 K2=1,NTENS
	     FM33X36(K1,K2)=ZERO
	     DO 888 K3=1,3
		 FM33X36(K1,K2)=FM33X36(K1,K2)+FS_I(K1,K3)*FM36X66(K3,K2)
888        CONTINUE
777     CONTINUE
C        
C     
       DO K1=1,NTENS
       FPHI_SIGMA(K1)=ZERO
       ENDDO	 
       DO 555 K1=1,NTENS
        DO 666 K2=1,NDI
        FPHI_SIGMA(K1)=FPHI_SIGMA(K1)+FPHI_S(K2)*FM33X36(K2,K1)
C
666     CONTINUE
555   CONTINUE     		
C      
      DO K1=1,NTENS
C          WRITE(6,*),'FPHI_SIGMA(K1)',FPHI_SIGMA(K1)
      ENDDO
C
C
C     1  'FPHI_SIGMA(2)',FPHI_SIGMA(2)
      DO K1=1,NTENS
	 DDFDDS(K1)=ZERO
      ENDDO
	 DO K1=1,NTENS
C	 **********************************
	 DDFDDS(K1)=FPHI_SIGMA(K1)/(2*FYM*(FEQSTRESS**(FYM-1)))
C            
C       WRITE(6,*),'DDFDDS(K)',DDFDDS(K1)
       ENDDO
C************************************************************
C*******CALCULATE THE ANISOTROPIC HARDENING PARAMETER*******
       DO K1=1,NTENS
           DO K2=1,NTENS
           FL_P(K1,K2)=ZERO
           ENDDO
       ENDDO
       FL_P(1,2)=(SIGMA_X-SIGMA_Z)/3
       FL_P(1,3)=(SIGMA_X-SIGMA_Y)/3
       FL_P(2,1)=(SIGMA_Y-SIGMA_Z)/3
       FL_P(2,3)=(SIGMA_Y-SIGMA_X)/3
       FL_P(3,1)=(SIGMA_Z-SIGMA_Y)/3
       FL_P(3,2)=(SIGMA_Z-SIGMA_X)/3
       FL_P(4,4)=SIGMA_XY
       FL_P(5,5)=SIGMA_ZX
       FL_P(6,6)=SIGMA_YZ
C      
      DO K1=1,NTENS
       FL_LAMDA(K1)=ZERO
       ENDDO	 
       DO 149 K1=1,NTENS
        DO 249 K2=1,NTENS
        FL_LAMDA(K1)=FL_LAMDA(K1)+FP_LAMDA(K2)*FL_P(K2,K1)
C
249     CONTINUE
149   CONTINUE 
C      
       DO K1=1,NDI
       FI_LAMDA(K1)=ZERO
       ENDDO	 
       DO 500 K1=1,NDI
        DO 501 K2=1,NTENS
        FI_LAMDA(K1)=FI_LAMDA(K1)+FL_LAMDA(K2)*FI_L(K2,K1)
C
501     CONTINUE
500   CONTINUE 
C
      DO K1=1,NDI
       FS_LAMDA(K1)=ZERO
       ENDDO	 
       DO 510 K1=1,NDI
        DO 511 K2=1,NDI
        FS_LAMDA(K1)=FS_LAMDA(K1)+FI_LAMDA(K2)*FS_I(K2,K1)
C
511     CONTINUE
510   CONTINUE 
C     
      FPHI_LAMDA=ZERO
      DO K1=1,NDI
      FPHI_LAMDA=FPHI_LAMDA+FS_LAMDA(K1)* FPHI_S(K1)
      ENDDO
C    
      FSIGMA_LAMDA=FPHI_LAMDA/(2*FYM*(FEQSTRESS**(FYM-1)))
C**********************************************************************       
C      
C***** IF FI2 FI3 EQ ZERO
      ENDIF
1234  CONTINUE

      
C	FIRST ORDER G-SIGMA
      ALPHAG=THREE/(TWO*GEQSTRESS)
      DDGDDS(1)=ALPHAG*GSIGMA_SX
	DDGDDS(2)=ALPHAG*GSIGMA_SY
	DDGDDS(3)=ALPHAG*GSIGMA_SZ
      DDGDDS(4)=ALPHAG*TWO*GTAO_SXY
	DDGDDS(5)=ALPHAG*TWO*GTAO_SZX
	DDGDDS(6)=ALPHAG*TWO*GTAO_SYZ
C
       DO K1=1,NTENS
c       WRITE(6,*),'DDFDDS(K1)',DDFDDS(K1)
       ENDDO
C**************************************************
C    
      DO K1=1,NTENS
C         WRITE(6,*),'DDFDDS(K1)',DDFDDS(K1),'DDGDDS(K1)',DDGDDS(K1)
      ENDDO 
C      WRITE(6,*),'GSIGMA_SX',GSIGMA_SX,'FL1',FL1
C***************************************************      
C
C     CALCULATE [D][DDGDDS]=DDGDS
      DO K1=1,NTENS
	     DO K2=1,NTENS
		   DDGDS(K1)=DDGDS(K1)+DDSDDE(K1,K2)*DDGDDS(K2)
           ENDDO
      ENDDO
C      WRITE(7,*),'DDGDS',DDGDS(1)
C***
      DO K1=1,NTENS   
          HI(K1)=CI*(PSTRES(K1)-AI(K1)-AII(K1))/GEQSTRESS-V*AI(K1)
          HII(K1)=CII*(PSTRES(K1)-AI(K1)-AII(K1))/GEQSTRESS
      ENDDO
C*******      *******         *******              *********
C      
C     CALCULATE NBLFM2--A'                 *****QUESTION

       BBB=1
C      
       BBB=GEQSTRESS/FEQSTRESS
	 NBLFM3=HARD0*BBB
       NBLFM2=FSIGMA_LAMDA*BBB
c       NBLFM2=zero
C     CALCULATE NBLFM1 F-D-G
      DO K1=1,NTENS
         NBLFM1=NBLFM1+DDFDDS(K1)*DDGDS(K1)
      END DO
C     CALCULATE DNBL (D LAMEDA)
C      DNBL=FB/(NBLFM1-NBLFM2+NBLFM3)
      DNBL=ZERO
C      
      IDNBL=DNBL
C    
      DO KK=1,1
C  K CHANGE  T UNCHANGE
      KDNBL=EQPLAS
      TDNBL=EQPLAS
C      WRITE(6,*),'********',KDNBL,TDNBL
C
      IF(KK.GT.1) THEN
      DO K1=1,NTENS
          DDFDDS(K1)=NF(K1)
          DDGDDS(K1)=NG(K1)
          DDGDS(K1)=NDDGDS(K1)
      ENDDO
      DNBL=NDNBL
      ELSE
      DNBL=IDNBL    
      CONTINUE
      ENDIF
C      
C	PULL THE STRESS TO THE YIELD SURFACE*****DDGDS IS ON B POINT
      DO K1=1,NTENS
         KSTRESS(K1)=TSTRESS(K1)-DNBL*DDGDS(K1)
         KBSTRES(K1)=TBSTRES(K1)-DNBL*DDGDS(K1)
         EPLASOLD(K1)=EPLAS(K1)
          AIOLD(K1)=AI(K1)
          AIIOLD(K1)=AII(K1)
      END DO
C      WRITE(6,*),'FEQ',FEQSTRESS,'GEQ',GEQSTRESS  
C
C*********************************************************************
C     DO  NEWTON ITERATION FOR THE SATISFIED STRESS
C*********************************************************************	  
      DO KNEWRTON=1,NEWTON
	   SIGMA_X=KSTRESS(1)-AI(1)-AII(1)
	   SIGMA_Y=KSTRESS(2)-AI(2)-AII(2)
         SIGMA_Z=KSTRESS(3)-AI(3)-AII(3)
         TAO_XY=KSTRESS(4)-AI(4)-AII(4)
         TAO_ZX=KSTRESS(5)-AI(5)-AII(5)	   
         TAO_YZ=KSTRESS(6)-AI(6)-AII(6)
C         WRITE(7,*),'STRESS@C XYZ',SIGMA_X, SIGMA_Y, SIGMA_Z
C     
      SIGMA0K=FD0*FA0*(TDNBL+FB0)**FC0+
     1 (1-FD0)*(FE0+FF0*(1-EXP(-FG0*TDNBL)))
C      
      SIGMA90K=FD90*FA90*(TDNBL+FB90)**FC90+
     1 (1-FD90)*(FE90+FF90*(1-EXP(-FG90*TDNBL)))
C      
      SIGMA45K=FD45*FA45*(TDNBL+FB45)**FC45+
     1 (1-FD45)*(FE45+FF45*(1-EXP(-FG45*TDNBL)))
C      

      SIGMA0K=FA0+FB0*(1-exp(-FC0*TDNBL))-CI*(1-exp(-V*TDNBL))/V   
      SIGMA90K=FA90-FB90*EXP(-FC90*TDNBL)+FD90*TDNBL     
      SIGMA45K=FA45-FB45*EXP(-FC45*TDNBL)+FD45*TDNBL      
      SIGMABK=FAB-FBB*EXP(-FCB*TDNBL)+FDB*TDNBL
C      
      IF(TDNBL.EQ.0.0)THEN
      HARD0K=EMOD
      HARD90K=EMOD
      HARD45K=EMOD
      HARDBK=EMOD
      ELSE
      HARD0K=FB0*FC0*exp(-FC0*TDNBL)-CI*exp(-V*TDNBL)
      HARD45K=FC45*FB45*EXP(-FC45*TDNBL)+FD45
      HARD90K=FC90*FB90*EXP(-FC90*TDNBL)+FD90
      HARDBK=FCB*FBB*EXP(-FCB*TDNBL)+FDB       
      ENDIF
      
      
C      
      CALL PARACAL(SIGMA0K,SIGMA90K,SIGMABK,SIGMA45K,FYM,x)
C     
      FYAK=x(1)
      FYBK=x(2)
      FYCK=x(3)
      FYHK=x(4)
      FYGK=x(4)
      FYFK=x(4)
c      write(6,*),FYA,FYB,FYC,FYH,'FYAK,FYBK',FYAK,FYBK,FYCK,FYHK
      CALL PARACALAH(FYM,w,SIGMA0K,SIGMA90K,SIGMABK,SIGMA45K,
     1 FYAK,FYBK,FYCK,FYHK,HARD0K,HARD90K,HARDBK,HARD45K)
c       write(6,*),x
C
      AAHK=w(1)
      BBHK=w(2)
      CCHK=w(3)
      HHHK=w(4)
      GGHK=w(4)
      FFHK=w(4)
      DO K1=1,NTENS
          FP_LAMDA(K1)=ZERO
      ENDDO
      FP_LAMDA(1)=AAHK
      FP_LAMDA(2)=BBHK
      FP_LAMDA(3)=CCHK
      FP_LAMDA(4)=HHHK
      FP_LAMDA(5)=GGHK
      FP_LAMDA(6)=FFHK
C      
C******         ********        ********      
      YIELDK=FA0+FB0*(1-exp(-FC0*KDNBL))-CI*(1-exp(-V*KDNBL))/V  
C      
C      
C******      
      SIGMA_SA=SIGMA_Y-SIGMA_Z
      SIGMA_SB=SIGMA_Z-SIGMA_X
      SIGMA_SC=SIGMA_X-SIGMA_Y
      SIGMA_SF=TAO_YZ
      SIGMA_SG=TAO_ZX
      SIGMA_SH=TAO_XY
C     
      FL1=(FYCK*SIGMA_SC-FYBK*SIGMA_SB)/THREE
      FL2=(FYAK*SIGMA_SA-FYCK*SIGMA_SC)/THREE
      FL3=(FYBK*SIGMA_SB-FYAK*SIGMA_SA)/THREE
      FL4=FYHK*SIGMA_SH
      FL5=FYGK*SIGMA_SG
      FL6=FYFK*SIGMA_SF
C     
      FI1=ZERO
      FI2=(FL4**2+FL5**2+FL6**2)/3+(FL1**2+FL2**2+FL3**2)/6
      FI3=(FL1*FL2*FL3+2*FL4*FL5*FL6-FL3*FL4**2-FL2*FL5**2-FL1*FL6**2)/2
C     
      FRHO=FI2**(ONEHALF)
      FTHETA=ACOSD(FI3/(FI2**ONEHALF))
      FRHO3=FI2**(-0.5)
C
      FS1=TWO*SQRT(FI2)*COSD(FTHETA/THREE)
      FS2=TWO*SQRT(FI2)*COSD((2*PI-FTHETA)/THREE)
      FS3=TWO*SQRT(FI2)*COSD((2*PI+FTHETA)/THREE)
C
C ******     
      FPHI=ABS(FS1-FS2)**FYM+ABS(FS2-FS3)**FYM+ABS(FS3-FS1)**FYM      
      FEQSTRESS=(FPHI/2)**(1/FYM)
c
C**********
C      
C     CALCULATE G FUNCTION ON POINT C
      GSIGMA_SX=(GH*(SIGMA_X-SIGMA_Y)+GG*(SIGMA_X-SIGMA_Z))*GSUM
	GSIGMA_SY=(GF*(SIGMA_Y-SIGMA_Z)+GH*(SIGMA_Y-SIGMA_X))*GSUM
	GSIGMA_SZ=(GG*(SIGMA_Z-SIGMA_X)+GF*(SIGMA_Z-SIGMA_Y))*GSUM
	GTAO_SYZ=GL*TAO_YZ*GSUM
	GTAO_SZX=GM*TAO_ZX*GSUM
	GTAO_SXY=GN*TAO_XY*GSUM
	GEQSTRESS=GF*(SIGMA_Y-SIGMA_Z)**TWO+GG*(SIGMA_Z-SIGMA_X)**TWO
     1 +GH*(SIGMA_X-SIGMA_Y)**TWO+TWO*GL*TAO_YZ**TWO+
     2 TWO*GM*TAO_ZX**TWO+TWO*GN*TAO_XY**TWO
	GEQSTRESS=DSQRT(THREE*GSUM*GEQSTRESS/TWO)
C
C      WRITE(6,*),'FEQ',FEQSTRESS,'GEQ',GEQSTRESS
C*******************************************************
      IF ((FI2.EQ.ZERO).and.(FI3.EQ.ZERO)) THEN
      DO K1=1,NTENS
	   DDFDDS(K1)=ZERO
      ENDDO
      FSIGMA_LAMDA=ZERO
      GOTO 2234
C
C********************************************************      
      ELSE
C      FIRST ORDER F-SIGMA
      DO K1=1,NDI
	  FPHI_S(K1)=ZERO
      ENDDO
        FPHI_S(1)=((FS1-FS2)*(ABS(FS1-FS2)**(FYM-2))-
     1 (FS3-FS1)*(ABS(FS3-FS1)**(FYM-2)))*FYM
        FPHI_S(2)=((FS2-FS3)*(ABS(FS2-FS3)**(FYM-2))-
     1 (FS1-FS2)*(ABS(FS1-FS2)**(FYM-2)))*FYM
        FPHI_S(3)=((FS3-FS1)*(ABS(FS3-FS1)**(FYM-2))-
     1 (FS2-FS3)*(ABS(FS2-FS3)**(FYM-2)))*FYM
C
      DO K1=1,NDI
	     DO K2=1,NDI
	     FS_I(K1,K2)=ZERO
	     ENDDO
      ENDDO
      FS_I(1,1)=ZERO
      FS_I(1,2)=FRHO3*SIND(2*FTHETA/THREE)/SIND(FTHETA)
      FS_I(1,3)=2*FRHO3*FRHO3*SIND(FTHETA/THREE)/(3*SIND(FTHETA))
      FS_I(2,1)=ZERO
      FS_I(2,2)=FRHO3*SIND(2*(FTHETA+PI)/THREE)/SIND(FTHETA)
      FS_I(2,3)=2*FRHO3*FRHO3*SIND((FTHETA-2*PI)/THREE)/(3*SIND(FTHETA))
      FS_I(3,1)=ZERO
      FS_I(3,2)=FRHO3*SIND(2*(FTHETA-PI)/THREE)/SIND(FTHETA)
      FS_I(3,3)=2*FRHO3*FRHO3*SIND((FTHETA+2*PI)/THREE)/(3*SIND(FTHETA))
C
C     
      DO K1=1,NDI
	     DO K2=1,NTENS
	     FI_L(K1,K2)=ZERO
	     ENDDO
      ENDDO 
      FI_L(2,1)=FL1/3
      FI_L(2,2)=FL2/3
      FI_L(2,3)=FL3/3
      FI_L(2,4)=2*FL4/3
      FI_L(2,5)=2*FL5/3
      FI_L(2,6)=2*FL6/3
      FI_L(3,1)=(FL2*FL3-FL6*FL6)/2
      FI_L(3,2)=(FL1*FL3-FL5*FL5)/2
      FI_L(3,3)=(FL1*FL2-FL4*FL4)/2
      FI_L(3,4)=FL5*FL6-FL3*FL4
      FI_L(3,5)=FL4*FL6-FL2*FL5
      FI_L(3,6)=FL4*FL5-FL1*FL6
C	  
      DO K1=1,NTENS
	     DO K2=1,NTENS
	     FL_SIGMA(K1,K2)=ZERO
	     ENDDO
      ENDDO 
      FL_SIGMA(1,1)=(FYBK+FYCK)/3
      FL_SIGMA(1,2)=-FYCK/3
      FL_SIGMA(1,3)=-FYBK/3
      FL_SIGMA(2,1)=-FYCK/3
      FL_SIGMA(2,2)=(FYAK+FYCK)/3
      FL_SIGMA(2,3)=-FYAK/3
      FL_SIGMA(3,1)=-FYBK/3
      FL_SIGMA(3,2)=-FYAK/3
      FL_SIGMA(3,3)=(FYAK+FYBK)/3
      FL_SIGMA(4,4)=FYHK
      FL_SIGMA(5,5)=FYGK
      FL_SIGMA(6,6)=FYFK
C	  
      DO 509 K1=1,NDI
      DO 509 K2=1,NTENS
	     FM36X66(K1,K2)=ZERO
	     DO 508 K3=1,NTENS
		 FM36X66(K1,K2)=FM36X66(K1,K2)+FI_L(K1,K3)*FL_SIGMA(K3,K2)
508     CONTINUE
509   CONTINUE   
C      
	  DO 517 K1=1,NDI
	  DO 517 K2=1,NTENS
	     FM33X36(K1,K2)=ZERO
	     DO 516 K3=1,NDI
		 FM33X36(K1,K2)=FM33X36(K1,K2)+FS_I(K1,K3)*FM36X66(K3,K2)
516     CONTINUE
517   CONTINUE
C     
       DO K1=1,NTENS
       FPHI_S(K1)=ZERO
       ENDDO	 
       DO 526 K1=1,NTENS
        DO 525 K2=1,NDI
       FPHI_SIGMA(K1)=FPHI_SIGMA(K1)+FPHI_S(K1)*FM33X36(K2,K1)
525     CONTINUE
526    CONTINUE     		
C
      DO K1=1,NTENS
	 DDFDDS(K1)=ZERO
      ENDDO
	 DO K1=1,NTENS
C	 **********************************
	 DDFDDS(K1)=FPHI_SIGMA(K1)/(2*FYM*(FEQSTRESS**(FYM-1)))
C	 ********************************** 
       ENDDO
C
C  ********************************************************    
C*******CALCULATE THE ANISOTROPIC HARDENING PARAMETER*******
       DO K1=1,NTENS
           DO K2=1,NTENS
           FL_P(K1,K2)=ZERO
           ENDDO
       ENDDO
       FL_P(1,2)=(SIGMA_X-SIGMA_Z)/3
       FL_P(1,3)=(SIGMA_X-SIGMA_Y)/3
       FL_P(2,1)=(SIGMA_Y-SIGMA_Z)/3
       FL_P(2,3)=(SIGMA_Y-SIGMA_X)/3
       FL_P(3,1)=(SIGMA_Z-SIGMA_Y)/3
       FL_P(3,2)=(SIGMA_Z-SIGMA_X)/3
       FL_P(4,4)=SIGMA_XY
       FL_P(5,5)=SIGMA_ZX
       FL_P(6,6)=SIGMA_YZ
C      
      DO K1=1,NTENS
       FL_LAMDA(K1)=ZERO
       ENDDO	 
       DO 833 K1=1,NTENS
        DO 834 K2=1,NTENS
        FL_LAMDA(K1)=FL_LAMDA(K1)+FP_LAMDA(K2)*FL_P(K2,K1)
C
834     CONTINUE
833   CONTINUE 
C      
       DO K1=1,NDI
       FI_LAMDA(K1)=ZERO
       ENDDO	 
       DO 5001 K1=1,NDI
        DO 5011 K2=1,NTENS
        FI_LAMDA(K1)=FI_LAMDA(K1)+FL_LAMDA(K2)*FI_L(K2,K1)
C
5011     CONTINUE
5001   CONTINUE 
C
      DO K1=1,NDI
       FS_LAMDA(K1)=ZERO
       ENDDO	 
       DO 854 K1=1,NDI
        DO 855 K2=1,NDI
        FS_LAMDA(K1)=FS_LAMDA(K1)+FI_LAMDA(K2)*FS_I(K2,K1)
C
855     CONTINUE
854   CONTINUE 
C     
      FPHI_LAMDA=ZERO
      DO K1=1,NDI
      FPHI_LAMDA=FPHI_LAMDA+FS_LAMDA(K1)* FPHI_S(K1)
      ENDDO
C    
      FSIGMA_LAMDA=FPHI_LAMDA/(2*FYM*(FEQSTRESS**(FYM-1)))  
C      IF ((FI2.EQ.ZERO).and.(FI3.EQ.ZERO)) 2234
C2234  CONTINUE 
      ENDIF
2234  CONTINUE       
C	FIRST ORDER G-SIGMA
      ALPHAG=THREE/(TWO*GEQSTRESS)
      DDGDDS(1)=ALPHAG*GSIGMA_SX
	DDGDDS(2)=ALPHAG*GSIGMA_SY
	DDGDDS(3)=ALPHAG*GSIGMA_SZ
      DDGDDS(4)=ALPHAG*TWO*GTAO_SXY
	DDGDDS(5)=ALPHAG*TWO*GTAO_SZX
	DDGDDS(6)=ALPHAG*TWO*GTAO_SYZ	 
C
      DO K1=1,NTENS
C         WRITE(6,*),'DDFDDS(K1)',DDFDDS(K1),'DDGDDS(K1)',DDGDDS(K1)
      ENDDO
C**************************************************
C***************************************************  
C     CALCULATE [D][DDGDDS]=DDGDS ON POINT C. INITIALIZE 0
      DO K1=1,NTENS  
	   DDGDS(K1)=ZERO
      ENDDO	
      DO K1=1,NTENS
	   DO K2=1,NTENS
		 DDGDS(K1)=DDGDS(K1)+DDSTIFF(K1,K2)*DDGDDS(K2)
         ENDDO
      ENDDO	 
C******CALCULATE YIELD SURFACE AND PARAMETERS***********************
       BBB=GF*((GG*DDGDDS(2)-GH*DDGDDS(3))*GAA)**TWO+
     1 GG*((GH*DDGDDS(3)-GF*DDGDDS(1))*GAA)**TWO+
     2 GH*((GF*DDGDDS(1)-GG*DDGDDS(2))*GAA)**TWO
C
       BBB=BBB+DDGDDS(4)*DDGDDS(4)*HALF/GN+
     1 DDGDDS(5)*DDGDDS(5)*HALF/GM+
     2 DDGDDS(6)*DDGDDS(6)*HALF/GL
C
       BBB=DSQRT(TWO*(GF+GG+GH)*BBB/THREE)
c
       BBB=1
       BBB=FEQSTRESS/GEQSTRESS
C       WRITE(7,*),'NBLFM1&2',NBLFM1,NBLFM2
C       WRITE(7,*),'DNBL',DNBL
C
      DO K1=1,NTENS   
          HI(K1)=CI*(KSTRESS(K1)-AI(K1)-AII(K1))/GEQSTRESS-V*AI(K1)
          HII(K1)=CII*(KSTRESS(K1)-AI(K1)-AII(K1))/GEQSTRESS
      ENDDO
C*** 
C************************************
C     
      V4=abs(2*dsqrt(((FYAK-FYBK)/12)**2+(FYHK/2)**2))**FYM+
     1 abs((FYAK+FYBK)/4-dsqrt(((FYAK-FYBK)/12)**2+(FYHK/2)**2))**FYM+
     2 abs((FYAK+FYBK)/4+dsqrt(((FYAK-FYBK)/12)**2+(FYHK/2)**2))**FYM-
     3 2*(SIGMA0/SIGMA45)**FYM
      STATEV(28)=V4
      STATEV(30)=BBB
      STATEV(29)=FK
C      WRITE(6,*),'KNEWTON',KNEWRTON,'FB,FC',FB,FC,'eqplas',EQPLAS
C      WRITE(6,*),'YIELDK',YIELDK
C************JUDGEMENT WHETHER IT SATISFIED AND END THE LOOP	   
C      IF(FC.LT.(TOLER+PHI(KK))*YIELDK) GOTO 100
C      IF(FC.LT.(TOLER)*YIELDK) GOTO 100
       
c
C***************
C      CAL HARD
c     
      NBLFM2=FSIGMA_LAMDA*BBB
c            
	NBLFM3=HARD0K*BBB
C***********	  
C     FROM HERE CALCULATE Q Q-1 AND OTHER TEMP     
C     SENCOND ORDER G-SIGMA ri
       DDBDDS(1,1)=ONEHALF*((GH+GG)*GSUM/GEQSTRESS-ONEHALF*GSIGMA_SX*
     1 GSIGMA_SX/(GEQSTRESS**3))
       DDBDDS(1,2)=ONEHALF*((-GH)*GSUM/GEQSTRESS-ONEHALF*GSIGMA_SX*
	1 GSIGMA_SY/(GEQSTRESS**3))
	 DDBDDS(1,3)=ONEHALF*((-GG)*GSUM/GEQSTRESS-ONEHALF*GSIGMA_SX*
	1 GSIGMA_SZ/(GEQSTRESS**3))
	 DDBDDS(1,4)=ONEHALF*(-THREE*GSIGMA_SX*GTAO_SXY/(GEQSTRESS**3))
	 DDBDDS(1,5)=ONEHALF*(-THREE*GSIGMA_SX*GTAO_SZX/(GEQSTRESS**3)) 
	 DDBDDS(1,6)=ONEHALF*(-THREE*GSIGMA_SX*GTAO_SYZ/(GEQSTRESS**3)) 
C      
       DDBDDS(2,1)=ONEHALF*((-GH)*GSUM/GEQSTRESS-ONEHALF*GSIGMA_SY*GSIGMA_SX/(GEQSTRESS**3))
       DDBDDS(2,2)=ONEHALF*((GF+GH)*GSUM/GEQSTRESS-ONEHALF*GSIGMA_SY*GSIGMA_SY/(GEQSTRESS**3))
	 DDBDDS(2,3)=ONEHALF*((-GF)*GSUM/GEQSTRESS-ONEHALF*GSIGMA_SY*
	1 GSIGMA_SZ/(GEQSTRESS**3))
	 DDBDDS(2,4)=ONEHALF*(-THREE*GSIGMA_SY*GTAO_SXY/(GEQSTRESS**3))
	 DDBDDS(2,5)=ONEHALF*(-THREE*GSIGMA_SY*GTAO_SZX/(GEQSTRESS**3)) 
	 DDBDDS(2,6)=ONEHALF*(-THREE*GSIGMA_SY*GTAO_SYZ/(GEQSTRESS**3)) 
C
       DDBDDS(3,1)=ONEHALF*((-GG)*GSUM/GEQSTRESS-ONEHALF*GSIGMA_SZ*
	1 GSIGMA_SX/(GEQSTRESS**3))
       DDBDDS(3,2)=ONEHALF*((-GF)*GSUM/GEQSTRESS-ONEHALF*GSIGMA_SZ*
     1 GSIGMA_SY/(GEQSTRESS**3))
	 DDBDDS(3,3)=ONEHALF*((GG+GF)*GSUM/GEQSTRESS-ONEHALF*GSIGMA_SZ*
	1 GSIGMA_SZ/(GEQSTRESS**3))
	 DDBDDS(3,4)=ONEHALF*(-THREE*GSIGMA_SZ*GTAO_SXY/(GEQSTRESS**3))
	 DDBDDS(3,5)=ONEHALF*(-THREE*GSIGMA_SZ*GTAO_SZX/(GEQSTRESS**3)) 
	 DDBDDS(3,6)=ONEHALF*(-THREE*GSIGMA_SZ*GTAO_SYZ/(GEQSTRESS**3))
C***
       DDBDDS(4,1)=ONEHALF*(-THREE*GTAO_SXY*GSIGMA_SX/(GEQSTRESS**3))
       DDBDDS(4,2)=ONEHALF*(-THREE*GTAO_SXY*GSIGMA_SY/(GEQSTRESS**3))
	 DDBDDS(4,3)=ONEHALF*(-THREE*GTAO_SXY*GSIGMA_SZ/(GEQSTRESS**3))
	 DDBDDS(4,4)=THREE*(GN*GSUM/GEQSTRESS-
	1 THREE*GTAO_SXY*GTAO_SXY/(GEQSTRESS**3))
	 DDBDDS(4,5)=THREE*(-THREE*GTAO_SXY*GTAO_SZX/(GEQSTRESS**3))
	 DDBDDS(4,6)=THREE*(-THREE*GTAO_SXY*GTAO_SYZ/(GEQSTRESS**3))	   
C***
       DDBDDS(5,1)=ONEHALF*(-THREE*GTAO_SZX*GSIGMA_SX/(GEQSTRESS**3))
       DDBDDS(5,2)=ONEHALF*(-THREE*GTAO_SZX*GSIGMA_SY/(GEQSTRESS**3))
	 DDBDDS(5,3)=ONEHALF*(-THREE*GTAO_SZX*GSIGMA_SZ/(GEQSTRESS**3))
	 DDBDDS(5,4)=THREE*(-THREE*GTAO_SZX*GTAO_SXY/(GEQSTRESS**3))
	 DDBDDS(5,5)=THREE*(GM*GSUM/GEQSTRESS-
	1 THREE*GTAO_SZX*GTAO_SZX/(GEQSTRESS**3))
	 DDBDDS(5,6)=THREE*(-THREE*GTAO_SZX*GTAO_SYZ/(GEQSTRESS**3))
C***
       DDBDDS(6,1)=ONEHALF*(-THREE*GTAO_SYZ*GSIGMA_SX/(GEQSTRESS**3))
       DDBDDS(6,2)=ONEHALF*(-THREE*GTAO_SYZ*GSIGMA_SY/(GEQSTRESS**3))
	 DDBDDS(6,3)=ONEHALF*(-THREE*GTAO_SYZ*GSIGMA_SZ/(GEQSTRESS**3))
	 DDBDDS(6,4)=THREE*(-THREE*GTAO_SYZ*GTAO_SXY/(GEQSTRESS**3))
	 DDBDDS(6,5)=THREE*(-THREE*GTAO_SYZ*GTAO_SZX/(GEQSTRESS**3))
	 DDBDDS(6,6)=THREE*(GL*GSUM/GEQSTRESS-
	1 THREE*GTAO_SYZ*GTAO_SYZ/(GEQSTRESS**3))
C
C	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
****************          RESIDUAL           **************** 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C----------------------------------------------------------------------
      DO K1=1,NTENS   
          AA(K1)=-EPLAS(K1)+EPLASOLD(K1)+DNBL*DDGDDS(K1)
          BBI(K1)=-AI(K1)+AIOLD(K1)+DNBL*HI(K1)
          BBII(K1)=-AII(K1)+AIIOLD(K1)+DNBL*HII(K1)
      ENDDO    
      AAA=0
      BBBI=0
      BBBII=0
      DO K1=1,NTENS
          AAA=AAA+AA(K1)*AA(K1)
          BBBI=BBBI+BBI(K1)*BBI(K1)
          BBBII=BBBII+BBII(K1)*BBII(K1)
      ENDDO
      AAA=sqrt(AAA/6)
      BBBI=sqrt(BBBI/6)
      BBBII=sqrt(BBBII/6)
C***	  
      	FK=FEQSTRESS-YIELDK
C***		
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
****************       RESIDUAL  END         **************** 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C**************************************  BI**************
      DO K1=1,NTENS
          DO K2=1,NTENS
             BI(K1,K2)=DNBL*DDBDDS(K1,K2)
          ENDDO
      ENDDO 	
C**************************************  BI**************
C
C****************************DHIDS
c 
      DHIDS=ZERO
c ****************************************************     
C
      DO K1=1,NTENS
		 DHIDS=DHIDS+DDGDDS(K1)*(KSTRESS(K1)-AI(K1)-AII(K1))
      ENDDO	
c      
      DHIDS=CI*(GEQSTRESS-DHIDS)/(GEQSTRESS**2)
c
C****************************DHIDS
c
      BII=DNBL*DHIDS*BBB
C**************************************  BII**************
c
C**************************************  BIII**************
      BIII=DNBL*V+1
C**************************************  BIII**************
c
c
C**************************************  INVXX**************    
c      
C****************************INVDDSTIFF      
      DO K1=1,NTENS
         DO K2=1,NTENS
		 INVDDSTIFF(K1,K2)=DDSTIFF(K1,K2)
	   ENDDO
      ENDDO	
c      
      CALL BRINV(INVDDSTIFF,6,L,IS,JS)
C****************************INVDDSTIFF 
c      
C****************************XX      
      DO K1=1,NTENS
          DO K2=1,NTENS
       XX(K1,K2)=BIII*BI(K1,K2)+(BII+BIII+C*BII*BIII)*INVDDSTIFF(K1,K2)
          ENDDO
      ENDDO 
C****************************XX
      DO K1=1,NTENS
         DO K2=1,NTENS
		 INVXX(K1,K2)=XX(K1,K2)
	   ENDDO
      ENDDO	
c      
      CALL BRINV(INVXX,6,L,IS,JS)
C**************************************  INVXX**************
c
c
C**************************************  A11**************
      DO K1=1,NTENS
          DO K2=1,NTENS            
             A11(K1,K2)=ZERO
             A13(K1,K2)=ZERO
             A14(K1,K2)=ZERO
             A31(K1,K2)=ZERO
             A31s(K1,K2)=ZERO
             A33(K1,K2)=ZERO
             A34(K1,K2)=ZERO
             A41(K1,K2)=ZERO
             A43(K1,K2)=ZERO
             A44s(K1,K2)=ZERO
             A44(K1,K2)=ZERO
          ENDDO
      ENDDO 
CC    TAKE ZERO     
      DO K1=1,NTENS
          DO K2=1,NTENS
             A11(K1,K2)=(BIII+BII*(1+C*BIII))*INVXX(K1,K2)
          ENDDO
      ENDDO 	
C**************************************  A11**************
c
c
C**************************************  A13**************
c      
c      
      DO K1=1,NTENS
          DO K2=1,NTENS
              DO K3=1,NTENS
              A13(K1,K2)=A13(K1,K2)-INVXX(K1,K3)*BI(K3,K2)
              ENDDO
          ENDDO
      ENDDO 	
C**************************************  A13**************
c
c
C**************************************  A14************** 
      DO K1=1,NTENS
          DO K2=1,NTENS
              A14(K1,K2)=A13(K1,K2)*BIII
          ENDDO
      ENDDO 	
C**************************************  A14**************
c
c
C**************************************  A31************** 
      DO K1=1,NTENS
          DO K2=1,NTENS
              A31(K1,K2)=INVXX(K1,K2)*BII
          ENDDO
      ENDDO 	
C**************************************  A31**************
c
c
C**************************************  A33************** 
C**************************** A33s
      DO K1=1,NTENS
          DO K2=1,NTENS
              A33s(K1,K2)=BI(K1,K2)+(1+C*BII)*INVDDSTIFF(K1,K2)
          ENDDO
      ENDDO 
C**************************** A33s          
c      
      DO K1=1,NTENS
          DO K2=1,NTENS
              DO K3=1,NTENS
              A33(K1,K2)=A33(K1,K2)-INVXX(K1,K3)*A33s(K3,K2)
              ENDDO
          ENDDO
      ENDDO 	
C**************************************  A33**************
c
c
C**************************************  A34************** 
c 
      DO K1=1,NTENS
          DO K2=1,NTENS
              DO K3=1,NTENS
      A34(K1,K2)=A34(K1,K2)+BII*INVXX(K1,K3)*INVDDSTIFF(K3,K2)
              ENDDO
          ENDDO
      ENDDO 	
C**************************************  A34**************
c
c
C**************************************  A41************** 
      DO K1=1,NTENS
          DO K2=1,NTENS
              A41(K1,K2)=INVXX(K1,K2)*C*BII*BIII
          ENDDO
      ENDDO 	
C**************************************  A41**************
c
c
C**************************************  A43************** 	   
c
      DO K1=1,NTENS
          DO K2=1,NTENS
              DO K3=1,NTENS
          A43(K1,K2)=A43(K1,K2)+C*BII*INVXX(K1,K3)*INVDDSTIFF(K3,K2)
              ENDDO
          ENDDO
      ENDDO 	
C**************************************  A43**************
c
C**************************************  A44************** 
C**************************** A44s
      DO K1=1,NTENS
          DO K2=1,NTENS
              A44s(K1,K2)=BIII*BI(K1,K2)+(BII+BIII)*INVDDSTIFF(K1,K2)
          ENDDO
      ENDDO 
C**************************** A44s
c
      DO K1=1,NTENS
          DO K2=1,NTENS
              DO K3=1,NTENS
              A44(K1,K2)=A44(K1,K2)-INVXX(K1,K3)*A44s(K3,K2)
              ENDDO
          ENDDO
      ENDDO 	
C**************************************  A44**************
     
C**************************************  CNBLFZ      
      DO K1=1,NTENS
       A11AA(K1)=0
       A13BBI(K1)=0
       A14BBII(K1)=0
c       
       A31AA(K1)=0
       A33BBI(K1)=0
       A34BBII(K1)=0
c       
       A41AA(K1)=0
       A43BBI(K1)=0
       A44BBII(K1)=0
      ENDDO
c      
      DO K1=1,NTENS
	     DO K2=1,NTENS
		    A11AA(K1)= A11AA(K1)+A11(K1,K2)*AA(K2)
              A13BBI(K1)= A13BBI(K1)+A13(K1,K2)*BBI(K2)
              A14BBII(K1)= A14BBII(K1)+A14(K1,K2)*BBII(K2)
c              
              A31AA(K1)= A31AA(K1)+A31(K1,K2)*AA(K2)
              A33BBI(K1)= A33BBI(K1)+A33(K1,K2)*BBI(K2)
              A34BBII(K1)= A34BBII(K1)+A34(K1,K2)*BBII(K2)
c              
              A41AA(K1)= A41AA(K1)+A41(K1,K2)*AA(K2)
              A43BBI(K1)= A43BBI(K1)+A43(K1,K2)*BBI(K2)
              A44BBII(K1)= A44BBII(K1)+A44(K1,K2)*BBII(K2)
           ENDDO
      ENDDO
c      
      DO K1=1,NTENS
       AFZ(K1)=-A11AA(K1)-A13BBI(K1)-A14BBII(K1)+
     1 A31AA(K1)+A33BBI(K1)+A34BBII(K1)+
     2 A41AA(K1)+A43BBI(K1)+A44BBII(K1)
      ENDDO
c      
      DO K1=1,NTENS
          MAFZ=ZERO
      ENDDO
c      
      DO K1=1,NTENS
          MAFZ=MAFZ+DDFDDS(K1)*AFZ(K1)
      ENDDO
C**************************************
C**************************************      
C      (FK=FEQSTRESS-YIELDD)
C**************************************
C**************************************
      CNBLFZ=FK+MAFZ
c
C**************************************  CNBLFZ
c
C**************************************  CNBLFM     
      DO K1=1,NTENS
       A11DDGDDS(K1)=ZERO
       A13HI(K1)=ZERO
       A14HII(K1)=ZERO
c       
       A31DDGDDS(K1)=ZERO
       A33HI(K1)=ZERO
       A34HII(K1)=ZERO
c       
       A41DDGDDS(K1)=ZERO
       A43HI(K1)=ZERO
       A44HII(K1)=ZERO
      ENDDO
c      
      DO K1=1,NTENS
	     DO K2=1,NTENS
		    A11DDGDDS(K1)=A11DDGDDS(K1)+A11(K1,K2)*DDGDDS(K2)
              A13HI(K1)=A13HI(K1)+A13(K1,K2)*HI(K2)
              A14HII(K1)=A14HII(K1)+A14(K1,K2)*HII(K2)
c              
              A31DDGDDS(K1)=A31DDGDDS(K1)+A31(K1,K2)*DDGDDS(K2)
              A33HI(K1)=A33HI(K1)+A33(K1,K2)*HI(K2)
              A34HII(K1)=A34HII(K1)+A34(K1,K2)*HII(K2)
c              
              A41DDGDDS(K1)=A41DDGDDS(K1)+A41(K1,K2)*DDGDDS(K2)
              A43HI(K1)=A43HI(K1)+A43(K1,K2)*HI(K2)
              A44HII(K1)=A44HII(K1)+A44(K1,K2)*HII(K2)
           ENDDO
      ENDDO
c
      DO K1=1,NTENS
       AFM(K1)=A11DDGDDS(K1)+A13HI(K1)+A14HII(K1)-
     1 A31DDGDDS(K1)-A33HI(K1)-A34HII(K1)-
     2 A41DDGDDS(K1)-A43HI(K1)-A44HII(K1)
      ENDDO
c      
      DO K1=1,NTENS
          MAFM=ZERO
      ENDDO
c      
      DO K1=1,NTENS
          MAFM=MAFM+DDFDDS(K1)*AFM(K1)
      ENDDO
c      
      CNBLFM=MAFM-NBLFM2+NBLFM3
C**************************************  CNBLFM
      CNBL=CNBLFZ/CNBLFM
C      write(7,*),'cnbl',CNBL
CCCCC--------------------------------------------------CCCCCC
****************           UPDATE            **************** 
CCCCC--------------------------------------------------CCCCCC
c      
C*****************************UPDATE LAMEDA*********
       DNBL=DNBL+CNBL
C------------------------------------------------
        KDNBL=TDNBL+DNBL*BBB
C---------------------------------  EQPLAS
C**************************************  BACKSTRESS/OLD*******
      DO K1=1,NTENS
          DAI(K1)=ZERO
          DAII(K1)=ZERO
      ENDDO
C      
      DO K1=1,NTENS
          DAI(K1)=-(A31AA(K1)+A33BBI(K1)+A34BBII(K1))-
     1 (A31DDGDDS(K1)+A33HI(K1)+A34HII(K1))*CNBL
      ENDDO
c      
      DO K1=1,NTENS
          DAII(K1)=-(A41AA(K1)+A43BBI(K1)+A44BBII(K1))-
     1 (A41DDGDDS(K1)+A43HI(K1)+A44HII(K1))*CNBL
      ENDDO
c      
      DO K1=1,NTENS
          AI(K1)=AI(K1)+DAI(K1)
          AII(K1)=AII(K1)+DAII(K1)
      ENDDO
C**************************************  BACKSTRESS*******
c
c
C**************************************  STRES*******
      DO K1=1,NTENS
          DSTRES(K1)=-(A11AA(K1)+A13BBI(K1)+A14BBII(K1))-
     1 (A11DDGDDS(K1)+A13HI(K1)+A14HII(K1))*CNBL
      ENDDO
C      write(7,*),'DSTRES1',DSTRES(1)
c      
      DO K1=1,NTENS   
          KSTRESS(K1)=KSTRESS(K1)+DSTRES(K1)  
      ENDDO
C**************************************  STRES*******
c      
c      
C**************************************  EPLAS/OLD*******
      DO K1=1,NTENS
          DO K2=1,NTENS
	     EPLAS(K1)=EPLAS(K1)-INVDDSTIFF(K1,K2)*DSTRES(K2)
          ENDDO
      ENDDO
C**************************************  EPLAS*******
c      
c      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
****************        UPDATE END           **************** 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  

c----------------------------------------------------------
C
C******************************************************
       IF(FK.LT.(TOLER)*YIELDK) GOTO 100
c        IF(FK.LT.(TOLER+PHI(KK))*YIELDK) GOTO 100
C       
C************************************************************       
C     WITH THE DO KNEWTONS    
C       WRITE(7,*),'**********************k',KNEWTON
      ENDDO
C    
C       WRITE(6,*),KK
C    AFTER NEWTON INTERATION  IF NOT CONVERGE 
       WRITE(7,2)NEWTON,FC
2     FORMAT(//,'***WARNING-PLASTICITY ALGORITHM DID NOT
     1 CONVERGE AFTER',I3,'ITERRATIONS',2X,'FC EQUALS',E11.4)
C ****
C     OUT THE LOOP
100   CONTINUE
C***SAVE DDFDDS DDGDDS
       DO K1=1,NTENS
          NF(K1)=DDFDDS(K1)
          NG(K1)=DDGDDS(K1)
          NDDGDS(K1)=DDGDS(K1)
      ENDDO
C
      NDNBL=DNBL
C    LOOP KK
      ENDDO
C    
C*****      UPDATE STRESS   ********
      DO K1=1,NTENS
          STRESS(K1)=KSTRESS(K1)
      ENDDO
      EQPLAS=KDNBL
C      
C
C*********UPDATE EELAS,EPLAS
C
      DO 53 K1=1,NTENS
C	     EPLAS(K1)=EPLAS(K1)+DNBL*DDGDDS(K1)
		 EELAS(K1)=EELAS(K1)-DNBL*DDGDDS(K1)
53    CONTINUE
C      WRITE(6,*),'EQPLAS',EQPLAS
C       
             
C***
C*********************************************************
C    CALCULATE THE STIFFNESS TANGENT MATRIX       ********
C*********************************************************
C
      DO K1=1,NTENS
          LEFT(K1)=ZERO
      ENDDO
C      
      DO K1=1,NTENS
          LEFT(K1)=A11DDGDDS(K1)+A13HI(K1)+A14HII(K1)
      ENDDO
C      
      DO K1=1,NTENS
          RIGHT(K1)=ZERO
      ENDDO
C      
      DO K1=1,NTENS
          DO K2=1,NTENS
       RIGHT(K1)=RIGHT(K1)+DDFDDS(K2)*(A11(K1,K2)-A13(K1,K2)
     1 -A14(K1,K2))
          ENDDO
      ENDDO
C      
      DO K1=1,NTENS
         DO K2=1,NTENS
         JCOBFZ(K1,K2)=LEFT(K1)*RIGHT(K2)
         ENDDO
      ENDDO
C      
      JCOBFM=CNBLFM
C      
      DO K1=1,NTENS
         DO K2=1,NTENS
         DDSDDE(K1,K2)=A11(K1,K2)-JCOBFZ(K1,K2)/JCOBFM
         ENDDO
      ENDDO
C
C------------------------------------------------
C     IF(FB.LT.YIELDD) THEN
C  ***************************************      
      ELSE 
      DO K1=1,NTENS
	     STRESS(K1)=PSTRES(K1)
      ENDDO
C    
      ENDIF
C 
C*********OUTPUT***************************
      DO K1=1,NTENS
      STATEV(K1)=EELAS(K1)
      STATEV(K1+NTENS)=EPLAS(K1)
      ENDDO
      STATEV(13)=EQPLAS
      STATEV(14)=FEQSTRESS
      STATEV(15)=GEQSTRESS
      DO K1=1,NTENS
          STATEV(K1+15)=STRESS(K1)
      ENDDO    
C****************************************************   
      DO K1=1,NTENS
          STATEV(K1+30)=AI(K1)
          STATEV(K1+36)=AII(K1)
      ENDDO     
C***
      
C***  CALCULATE FLD FACTOR
 
       DO K1=1,NTENS
         ETOTAL(K1)=EELAS(K1)+EPLAS(K1) 
         STATEV(K1+40)=ETOTAL(K1)
       ENDDO  
       
       SMAX=0.5*(ETOTAL(1)+ETOTAL(2))+SQRT((0.5*(ETOTAL(1)-ETOTAL(2)))
     1 **2+(ETOTAL(4)*0.5)**2) 
       SMIN=0.5*(ETOTAL(1)+ETOTAL(2))-SQRT((0.5*(ETOTAL(1)-ETOTAL(2)))
     1 **2+(ETOTAL(4)*0.5)**2)
       IF ((SMIN.GE.-0.136).and.(SMIN.LE.0.275)) THEN
       FLDF=0.23-0.04*COS(15.03*SMIN)-0.027*SIN(15.03*SMIN)
       ELSE
           IF(SMIN.LT.-0.136) THEN
           FLDF=-2*SMIN
           ELSE
                IF(SMIN.GT.0.275) THEN
                FLDF=SMIN 
                ENDIF
           ENDIF     
       ENDIF 
       FLDW=SMAX/FLDF
       
      STATEV(50)=FLDW
      STATEV(51)=SMAX
      STATEV(52)=SMIN
      STATEV(53)=FLDF
c******fld************************************************      
      
      RETURN
      END
C     
C
      SUBROUTINE BRINV(A,N,L,IS,JS)
      DIMENSION A(N,N),IS(N),JS(N)
      DOUBLE PRECISION A,T,D
      L=1
      DO 3100 K=1,N
          D=0.0
          DO 310 I=K,N
          DO 310 J=K,N
           IF(ABS(A(I,J)).GT.D)THEN
               D=ABS(A(I,J))
               IS(K)=I
               JS(K)=J
           ENDIF
310   CONTINUE
      IF((D+1.0).EQ.1.0)THEN
          L=0
          WRITE(*,320)
          RETURN
      ENDIF
320   FORMAT(1X,'ERR**NOT INV')
      DO 330 J=1,N
          T=A(K,J)
          A(K,J)=A(IS(K),J)
          A(IS(K),J)=T
330   CONTINUE
      DO 340 I=1,N
      T=A(I,K)
      A(I,K)=A(I,JS(K))
      A(I,JS(K))=T
340   CONTINUE
      A(K,K)=1/A(K,K)
      DO 350 J=1,N
          IF(J.NE.K)THEN
              A(K,J)=A(K,J)*A(K,K)
          ENDIF
350   CONTINUE
      DO 370 I=1,N
          IF(I.NE.K)THEN
              DO 360 J=1,N
                IF(J.NE.K)THEN
                   A(I,J)=A(I,J)-A(I,K)*A(K,J)
                ENDIF
360   CONTINUE
      ENDIF        
370   CONTINUE
      DO 380 I=1,N
          IF(I.NE.K)THEN
              A(I,K)=-A(I,K)*A(K,K)
          ENDIF
380   CONTINUE
3100  CONTINUE      
      DO 3130 K=N,1,-1
          DO 3110 J=1,N
           T=A(K,J)
           A(K,J)=A(JS(K),J)
           A(JS(K),J)=T
3110      CONTINUE
      DO 3120 I=1,N
          T=A(I,K)
          A(I,K)=A(I,IS(K))
          A(I,IS(K))=T
3120  CONTINUE
3130  CONTINUE
      RETURN
      END SUBROUTINE BRINV
C**************************************************************
      MODULE m_gauss
      CONTAINS
      SUBROUTINE lineq(A,b,x,N)
      integer::i,k,N
      integer::id_max 
      DIMENSION A(N,N),b(N),x(N)
      DIMENSION Aup(N,N),bup(N)
      DOUBLE PRECISION A,b,x,Aup,bup
C
      DIMENSION Ab(N,N+1)
      DIMENSION vtemp1(N+1),vtemp2(N+1)
      DOUBLE PRECISION Ab,vtemp1,vtemp2
      Ab(1:N,1:N)=A
      Ab(:,N+1)=b
      DO k=1,N-1
      elmax=dabs(Ab(k,k))
      id_max=k 
C
	   DO i=k+1,n
         IF (dabs(Ab(i,k))>elmax) THEN
         elmax=Ab(i,k)
         id_max=i
         END IF          
         END DO
C
      vtemp1=Ab(k,:)
      vtemp2=Ab(id_max,:)
      Ab(k,:)=vtemp2
      Ab(id_max,:)=vtemp1   
!#########################################################
          DO i=k+1,N
          temp=Ab(i,k)/Ab(k,k)
          Ab(i,:)=Ab(i,:)-temp*Ab(k,:)
         END DO
      END DO
      Aup(:,:)=Ab(1:N,1:N)
      bup(:)=Ab(:,N+1)
C
      call uptri(Aup,bup,x,n)
      END SUBROUTINE lineq

       SUBROUTINE uptri(A,b,x,N)
C       implicit real*8(a-z)
       integer::i,j,N
       DIMENSION A(N,N),b(N),x(N)
       DOUBLE PRECISION A,b,x
       x(N)=b(N)/A(N,N)
C
       DO i=n-1,1,-1 
       x(i)=b(i)
       DO j=i+1,N
       x(i)=x(i)-a(i,j)*x(j)
       END DO
       x(i)=x(i)/A(i,i)
       END DO

       END SUBROUTINE uptri

       END MODULE m_gauss
C
       MODULE  m_newton

       CONTAINS
       SUBROUTINE  solve(x,s00,s90,s45,sb,mm)  
       USE m_gauss
C
C       implicit real*8(a-z)
       integer::I,itmax=50000
       integer::N=4,j
       DIMENSION x(4),f(4),dx(4),df(4,4)
       DOUBLE PRECISION x,f,dx,df,s00,s45,s90,sb,mm
C  
        x=(/1d0,1d0,1d0,1d0/)
        tol=1d-7

        DO i=1,itmax
C
        call  func(f,x,mm,s00,s90,sb,s45)
        call  jac(df,x,mm)
        call  lineq(df,-f,dx,N)
!  x=x+dx
        DO j=1,N
        x(j)=x(j)+dx(j)
        enddo     
C
        dx4=dsqrt(dx(1)**2+dx(2)**2+dx(3)**2+dx(4)**2)
	  if (dx4<tol) exit
c	  if (dx4<tol) exit 
!    write(11,*),'dx4',dx4
!------
        END DO
        END SUBROUTINE solve

        SUBROUTINE func(f,x,mm,s00,s90,sb,s45)

C         implicit real*8(a-z)
         DIMENSION x(4),f(4)
         DOUBLE PRECISION x,f,mm,s00,s90,sb,s45
         DO i=1,4
          f(i)=0
           ENDDO
C
       f(1)=(abs(x(2)-x(3)))**mm+(abs(x(2)*2+x(3)))**mm+
     1 (abs(x(2)+2*x(3)))**mm-2*(3*s00/s00)**mm
       f(2)=(abs(x(1)-x(3)))**mm+(abs(x(1)+2*x(3)))**mm+
     2  (abs(2*x(1)+x(3)))**mm-2*(3*s00/s90)**mm
       f(3)=(abs(x(2)+2*x(1)))**mm+(abs(x(1)+2*x(2)))**mm+
     3  (abs(x(1)-x(2)))**mm-2*(3*s00/sb)**mm
       f(4)=(abs(2*sqrt(((x(1)-x(2))/12)**2+(x(4)/2)**2)))**mm+
     4  (abs((x(1)+x(2))/4-sqrt(((x(1)-x(2))/12)**2+(x(4)/2)**2)))**mm+
     5  (abs((x(1)+x(2))/4+sqrt(((x(1)-x(2))/12)**2+(x(4)/2)**2)))**mm
     6  -2*(s00/s45)**mm
      END SUBROUTINE func

       SUBROUTINE  jac(df,x,mm)
       DIMENSION x(4),df(4,4)
       DOUBLE PRECISION x,df,mm
       DO i=1,4
           DO j=1,4
        df(i,j)=0
           ENDDO
       ENDDO    
C
       df(1,1)=0
       df(1,2)=mm*(abs(x(2)-x(3)))**(mm-2)*(x(2)-x(3))+
     1 mm*(abs(x(2)+2*x(3)))**(mm-2)*(x(2)+2*x(3))+
     2  2*mm*(abs(2*x(2)+x(3)))**(mm-2)*(2*x(2)+x(3))
       df(1,3)=2*mm*(abs(x(2)+2*x(3)))**(mm-2)*(x(2)+2*x(3))
     1  -(mm*abs(x(2)-x(3)))**(mm-2)*(x(2)-x(3))+
     2  mm*(abs(2*x(2)+x(3)))**(mm-2)*(2*x(2)+x(3))
       df(1,4)=0

       df(2,1)=mm*(abs(x(1)-x(3)))**(mm-2)*(x(1)-x(3))+
     1  mm*(abs(x(1)+2*x(3)))**(mm-2)*(x(1)+2*x(3))
     2  +2*mm*(abs(2*x(1)+x(3)))**(mm-2)*(2*x(1)+x(3))
       df(2,2)=0
       df(2,3)=2*mm*(abs(x(1)+2*x(3)))**(mm-2)*(x(1)+2*x(3))-
     1 mm*(abs(x(1)-x(3)))**(mm-2)*(x(1)-x(3))+
     2 mm*(abs(2*x(1)+x(3)))**(mm-2)*(2*x(1)+x(3))
       df(2,4)=0
 
       df(3,1)=mm*(abs(x(1)-x(2)))**(mm-2)*(x(1)-x(2))+
     1 mm*(abs(x(1)+2*x(2)))**(mm-2)*(x(1)+2*x(2))
     2 +2*mm*(abs(2*x(1)+x(2)))**(mm-2)*(2*x(1)+x(2))
       df(3,2)=2*mm*(abs(x(1)+2*x(2)))**(mm-2)*(x(1)+2*x(2))-
     1 mm*(abs(x(1)-x(2)))**(mm-2)*(x(1)-x(2))+
     2 mm*(abs(2*x(1)+x(2)))**(mm-2)*(2*x(1)+x(2))
       df(3,3)=0
       df(3,4)=0


      df(4,3)=0

      df(4,1)=mm*(abs(x(1)/4+x(2)/4-(x(4)**2/4+(x(1)-x(2))**2/144)**
     1 (1/2)))**(mm-2)*(x(1)/4-(x(1)-x(2))/(144*(x(4)**2/4+(x(1)-x(2))
     2 **2/144)**(1/2)))*(x(1)/4+x(2)/4-(x(4)**2/4+(x(1)-x(2))**2/144)
     3 **(1/2))+mm*(abs(x(1)/4+x(2)/4+(x(4)**2/4+(x(1)-x(2))**2/144)**
     4 (1/2)))**(mm-2)*(x(1)/4+(x(1)-x(2))/(144*(x(4)**2/4+(x(1)-x(2))
     5 **2/144)**(1/2)))*(x(1)/4+x(2)/4+(x(4)**2/4+(x(1)-x(2))**2/144)
     6 **(1/2))+(2**mm*mm*((x(4)**2/4+(x(1)-x(2))**2/144)**(1/2))**
     7 (mm-1)*(x(1)-x(2)))/(144*(x(4)**2/4+(x(1)-x(2))**2/144)**(1/2))
      df(4,2)=mm*(abs(x(1)/4+x(2)/4-(x(4)**2/4+(x(1)-x(2))**2/144)**
     1 (1/2)))**(mm-2)*(x(2)/4+(x(1)-x(2))/(144*(x(4)**2/4+(x(1)-x(2))
     2 **2/144)**(1/2)))*(x(1)/4+x(2)/4-(x(4)**2/4+(x(1)-x(2))**2/144)
     3 **(1/2))+mm*(abs(x(1)/4+x(2)/4+(x(4)**2/4+(x(1)-x(2))**2/144)
     4 **(1/2)))**(mm-2)*(x(2)/4-(x(1)-x(2))/(144*(x(4)**2/4+(x(1)-x(2))
     5 **2/144)**(1/2)))*(x(1)/4+x(2)/4+(x(4)**2/4+(x(1)-x(2))**2/144)
     6 **(1/2))-(2**mm*mm*((x(4)**2/4+(x(1)-x(2))**2/144)**(1/2))
     7 **(mm-1)*(x(1)-x(2)))/(144*(x(4)**2/4+(x(1)-x(2))**2/144)**(1/2))
      df(4,4)=mm*(abs(x(1)/4+x(2)/4+(x(4)**2/4+(x(1)-x(2))**2/144)
     1 **(1/2)))**(mm-2)*(x(2)/4+x(4)/(4*(x(4)**2/4+(x(1)-x(2))**2/144)
     2 **(1/2)))*(x(1)/4+x(2)/4+(x(4)**2/4+(x(1)-x(2))**2/144)**(1/2))
     3 +(2**mm*x(4)*mm*((x(4)**2/4+(x(1)-x(2))**2/144)**(1/2))**
     4 (mm-1))/(4*(x(4)**2/4+(x(1)-x(2))**2/144)**(1/2))-(x(4)*
     5 mm*(abs(x(1)/4+x(2)/4-(x(4)**2/4+(x(1)-x(2))**2/144)**(1/2)))**
     6 (mm-2)*(x(1)/4+x(2)/4-(x(4)**2/4+(x(1)-x(2))**2/144)**(1/2)))/
     7 (4*(x(4)**2/4+(x(1)-x(2))**2/144)**(1/2))
! write(11,*),'df in jac',df

      END SUBROUTINE jac

      END MODULE m_newton
     
      SUBROUTINE PARACAL(SIGMA0,SIGMA90,SIGMAB,SIGMA45,FYM,x)
      USE m_newton
      DIMENSION x(4)
      DOUBLE PRECISION SIGMA0,SIGMA90,SIGMAB,SIGMA45,FYM   
      DOUBLE PRECISION x,s00,s45,s90,sb,mm
      s00=SIGMA0
      s90=SIGMA90
      s45=SIGMA45
      sb=SIGMAB
      mm=FYM
      call solve(x,s00,s90,s45,sb,mm)
      RETURN
      END SUBROUTINE PARACAL
C****************************************************************  
C****************************************************************
      module m_gaussAH
      contains
      subroutine lineqAH(A,b,w,N)
      integer::i,k,N
      integer::id_max  
      real*8::A(N,N),b(N),w(N)
      real*8::Aup(N,N),bup(N)
      real*8::Ab(N,N+1)
      real*8::vtemp1(N+1),vtemp2(N+1)
      Ab(1:N,1:N)=A
      Ab(:,N+1)=b
      do k=1,N-1
      elmax=dabs(Ab(k,k))
      id_max=k 
C
	do i=k+1,n
      if (dabs(Ab(i,k))>elmax) then
         elmax=Ab(i,k)
         id_max=i
      end if          
      end do
C
      vtemp1=Ab(k,:)
      vtemp2=Ab(id_max,:)
      Ab(k,:)=vtemp2
      Ab(id_max,:)=vtemp1   
      do i=k+1,N
      temp=Ab(i,k)/Ab(k,k)
       Ab(i,:)=Ab(i,:)-temp*Ab(k,:)
       end do
      end do
      Aup(:,:)=Ab(1:N,1:N)
      bup(:)=Ab(:,N+1)

      call uptriAH(Aup,bup,w,n)
      end subroutine lineqAH

      subroutine uptriAH(A,b,w,N)
      integer::i,j,N
       real*8::A(N,N),b(N),w(N)
      w(N)=b(N)/A(N,N)
      do i=n-1,1,-1 
       w(i)=b(i)
       do j=i+1,N
       w(i)=w(i)-a(i,j)*w(j)
       end do
       w(i)=w(i)/A(i,i)
       end do

       end subroutine uptriAH

       end module m_gaussAH
c      
       MODULE  m_newtonAH

       CONTAINS
       SUBROUTINE  solveAH(w,mm,s00,s90,sb,s45,aa,bb,cc,hh,hard0,
     1 hard90,hardb,hard45)  
       USE m_gaussAH
C
C       implicit real*8(a-z)
       integer::I,itmax=50000
       integer::N=4,j
       DIMENSION w(4),f(4),dw(4),df(4,4)
       DOUBLE PRECISION w,f,dw,df,s00,s45,s90,sb,mm
       DOUBLE PRECISION aa,bb,cc,hh,hard0,hard90,hardb,hard45
C  
        w=(/1d0,1d0,1d0,1d0/)
        tol=1d-8

        DO i=1,itmax
C
        call  funcAH(f,w,mm,s00,s90,sb,s45,aa,bb,cc,hh,hard0,
     1 hard90,hardb,hard45)
        call  jacAH(df,w,mm,aa,bb,cc,hh)
        call  lineqAH(df,-f,dw,N)
!  x=x+dx
        DO j=1,N
        w(j)=w(j)+dw(j)
        enddo     
C
        dw4=dsqrt(dw(1)**2+dw(2)**2+dw(3)**2+dw(4)**2)
	  if (dw4<tol) exit 
!    write(11,*),'dx4',dx4
!------
        END DO
        END SUBROUTINE solveAH

        SUBROUTINE funcAH(f,w,mm,s00,s90,sb,s45,aa,bb,cc,hh,hard0,
     1 hard90,hardb,hard45)
      DIMENSION w(4),f(4)
      DOUBLE PRECISION w,f,mm,s00,s90,sb,s45,aa,bb,cc,hh,hard0,
     1 hard90,hardb,hard45,R1,R2,R3,R4,TT 
      DO i=1,4
      f(i)=0
      ENDDO
      R1=2*(3**mm)*(s00/s00)**(mm-1)*(hard0*s00-hard0*s00)/(s00**2)
      R2=2*(3**mm)*(s00/s90)**(mm-1)*(hard0*s90-hard90*s00)/(s90**2)
      R3=2*(3**mm)*(s00/sb)**(mm-1)*(hard0*sb-hardb*s00)/(sb**2)
      R4=2*(s00/s45)**(mm-1)*(hard0*s45-hard45*s00)/(s45**2)
      TT=dsqrt(((aa-bb)/12)**2+(hh/2)**2)
c      write(*,*),mm,s00,aa,hard0,R1,R4,TT 

       f(1)=(bb-cc)*(abs(bb-cc))**(mm-2)*(w(2)-w(3))+(2*bb+cc)*
     1 (abs(2*bb+cc))**(mm-2)*(2*w(2)+w(3))+(bb+2*cc)*
     2  (abs(bb+2*cc))**(mm-2)*(w(2)+2*w(3))-R1
       f(2)=(cc-aa)*(abs(cc-aa))**(mm-2)*(w(3)-w(1))+(2*aa+cc)*
     1  (abs(2*aa+cc))**(mm-2)*(2*w(1)+w(3))+(aa+2*cc)*
     2  (abs(aa+2*cc))**(mm-2)*(w(1)+2*w(3))-R2
       f(3)=(aa-bb)*(abs(aa-bb))**(mm-2)*(w(1)-w(2))+(2*bb+aa)*
     1  (abs(2*bb+aa))**(mm-2)*(2*w(2)+w(1))+(bb+2*aa)*
     2  (abs(bb+2*aa))**(mm-2)*(w(2)+2*w(1))-R3
       f(4)=(2**mm)*TT**(mm-1)*((aa-bb)*(w(1)-w(2))/72+hh*w(4)/2)/(2*TT)
     1  +((aa+bb)/4-TT)*(abs((aa+bb)/4-TT))**(mm-2)*((w(1)+w(2))/4-
     2  ((aa-bb)*(w(1)-w(2))/72+hh*w(4)/2)/(2*TT))+((aa+bb)/4+TT)*
     3  (abs((aa+bb)/4+TT))**(mm-2)*((w(1)+w(2))/4+
     4  ((aa-bb)*(w(1)-w(2))/72+hh*w(4)/2)/(2*TT))-R4
!write(11,*),aa,bb,cc,hh,R1,R2,R3,R4,'********************',w
      END SUBROUTINE funcAH

      SUBROUTINE  jacAH(df,w,mm,aa,bb,cc,hh)
      DIMENSION w(4),df(4,4)
       DOUBLE PRECISION w,df,mm,aa,bb,cc,hh,TT
       TT=dsqrt(((aa-bb)/12)**2+(hh/2)**2)
       DO i=1,4
           DO j=1,4
        df(i,j)=0
           ENDDO
       ENDDO 
C      
       df(1,1)=0
       df(1,2)=(bb-cc)*(abs(bb-cc))**(mm-2)+(2*bb+cc)*
     1  (abs(2*bb+cc))**(mm-2)*2+(bb+2*cc)*(abs(bb+2*cc))**(mm-2)
       df(1,3)=-(bb-cc)*(abs(bb-cc))**(mm-2)+(2*bb+cc)*
     2  (abs(2*bb+cc))**(mm-2)+(bb+2*cc)*(abs(bb+2*cc))**(mm-2)*2
       df(1,4)=0
   
       df(2,1)=-(cc-aa)*(abs(cc-aa))**(mm-2)+(2*aa+cc)*
     1  (abs(2*aa+cc))**(mm-2)*2+(aa+2*cc)*(abs(aa+2*cc))**(mm-2)
       df(2,2)=0
       df(2,3)=(cc-aa)*(abs(cc-aa))**(mm-2)+(2*aa+cc)*
     2  (abs(2*aa+cc))**(mm-2)+(aa+2*cc)*(abs(aa+2*cc))**(mm-2)*2
       df(2,4)=0
       
       df(3,1)=(aa-bb)*(abs(aa-bb))**(mm-2)+(2*bb+aa)*
     1  (abs(2*bb+aa))**(mm-2)+(bb+2*aa)*(abs(bb+2*aa))**(mm-2)*2
       df(3,2)=-(aa-bb)*(abs(aa-bb))**(mm-2)+(2*bb+aa)*
     2  (abs(2*bb+aa))**(mm-2)*2+(bb+2*aa)*(abs(bb+2*aa))**(mm-2)
       df(3,3)=0
       df(3,4)=0
       
       df(4,1)=(2**mm)*TT**(mm-1)*((aa-bb)/72)/(2*TT)+((aa+bb)/4-TT)*
     1  (abs((aa+bb)/4-TT))**(mm-2)*(1/4-((aa-bb)/72)/(2*TT))+
     2  ((aa+bb)/4+TT)*(abs((aa+bb)/4+TT))**(mm-2)*
     3  (1/4+((aa-bb)/72)/(2*TT))
       df(4,2)=(2**mm)*TT**(mm-1)*(-(aa-bb)/72)/(2*TT)+
     1  ((aa+bb)/4-TT)*(abs((aa+bb)/4-TT))**(mm-2)*
     2  (1/4+((aa-bb)/72)/(2*TT))+mm*((aa+bb)/4+TT)*
     3  (abs((aa+bb)/4+TT))**(mm-2)*(1/4-((aa-bb)/72)/(2*TT))
       df(4,3)=0
       df(4,4)=(2**mm)*TT**(mm-1)*(hh/2)/(2*TT)+((aa+bb)/4-TT)*
     1  (abs((aa+bb)/4-TT))**(mm-2)*(-hh/2)/(2*TT)+((aa+bb)/4+TT)*
     2  (abs((aa+bb)/4+TT))**(mm-2)*(hh/2)/(2*TT)
!write(11,*),df,'df in jacAH'
      END SUBROUTINE jacAH  
      END MODULE m_newtonAH
C
C
      SUBROUTINE PARACALAH(FYM,w,SIGMA0,SIGMA90,SIGMAB,SIGMA45,
     1 FYA,FYB,FYC,FYH,HARD0,HARD90,HARDB,HARD45)
      USE m_newtonAH
      DIMENSION w(4)
      DOUBLE PRECISION SIGMA0,SIGMA90,SIGMAB,SIGMA45,FYM
      DOUBLE PRECISION FYA,FYB,FYC,FYH
      DOUBLE PRECISION HARD0,HARD90,HARDB,HARD45 
      DOUBLE PRECISION w,s00,s45,s90,sb,mm
      DOUBLE PRECISION aa,bb,cc,hh
      DOUBLE PRECISION PPP,QQQ,R4,XXX,TT
      s00=SIGMA0
      s90=SIGMA90
      s45=SIGMA45
      sb=SIGMAB
      hard0=HARD0
      hard90=HARD90
      hardb=HARDB
      hard45=HARD45
      aa=FYA
      bb=FYB
      cc=FYC
      hh=FYH
      mm=FYM
C      WRITE(*,*),x,mm,s00,s90,sb,s45,aa,bb,cc,hh,hard0,
C     1 hard90,hardb,hard45
       call solveAH(w,mm,s00,s90,sb,s45,aa,bb,cc,hh,hard0,
     1 hard90,hardb,hard45)  

      RETURN
      END SUBROUTINE PARACALAH
C
     