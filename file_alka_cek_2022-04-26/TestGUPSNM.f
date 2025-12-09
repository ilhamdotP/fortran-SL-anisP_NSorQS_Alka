      PROGRAM SNM        
c     ****************************************************
c     verified 14 August 2022 check "GUP"
c     
c     ****************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      INTEGER I
      CHARACTER(LEN=40) CI,FILENAME
      PI  = 3.14159265358979D0
      HC  = 197.327D0

c       BETA = 3.0D-7
c       BETA = 2.5D-7
        BETA = 1.0D-7

C------------------------------------------------------------------------------
C nuclear matter Saturation density !!!

c    BSP      
       KFO =1.31D0*HC

C------------------------------------------------------------------------------
       RHO = 2.D0*KFO*KFO*KFO/(3.D0*PI*PI)
     &     -4.D0*BETA*KFO**(5.D0)/(5.D0*PI*PI)
 
c      OPEN(unit=1,status='unknown',file='BE_SNM_BSP30.dat')
c      OPEN(unit=2,status='unknown',file='EOS_PNM_FSU2HX.dat')
c      OPEN(unit=3,status='unknown',file='SOS_PNM_FSU2HX.dat')

c 100  FORMAT (I3)
c        WRITE(CI,100)NINT(BETA*1D8)
c        FILENAME="basic_BSR12_SNM_GUP"//TRIM(ADJUSTL(CI))//".dat"
cc        FILENAME="basic_BSR12_PNM_GUP"//TRIM(ADJUSTL(CI))//".dat"
cc        WRITE(*,*)"File output:"
cc        WRITE(*,*)FILENAME
c        OPEN(unit=1,status='unknown',file=TRIM(FILENAME))
      
       open(unit=1,status='unknown',file='data_SNM.dat')
c       open(unit=1,status='unknown',file='data_PNM.dat')


        DO 10 I=83,1,-1        
        RHBO=I*0.1D0 
c         Write(*,*)I

c         RHBO=2.5D0
c          DO 10 I=97,103        
c         RHBO=I * 0.01D0
             
         RB=RHBO*RHO    

           
c         KF  = (1.5D0*PI*PI*RB)**(1.D0/3.D0)      
c     &      +(1.5D0*PI*PI*RB)*2.D0*BETA/5.D0

c         
C     CALL ENERGY DENSITY AND BINDING ENERGY
C     

        
         CALL FED(RHBO,ED,EB,EL,BETA)

C     
C     CALL PRESSURE
C     
  
         PRESS= RHO*RHBO*RHBO*func1(RHBO,BETA)/(HC*HC*HC)

c     CALL saturation  properties
         
         KNCM= 9.0D0*RHBO*RHBO*func2(RHBO,BETA)
         KJN0= 27.0D0*RHBO*RHBO*RHBO*func3(RHBO,BETA)
    
         SYME= ASIM(RHBO,BETA)
         LSY= 3.0D0*RHBO*lsym(RHBO,BETA)
         KSY= 9.0D0*RHBO*RHBO*ksym(RHBO,BETA)
         
         KASY=KSY-6.0D0*LSY
         KSAT2=KASY-KJN0/KNCM*LSY

C     RESULTS
C     
C     ****************************************************
C     RB is Baryon density,  EB is binding energy, ED is energy density
C     PRESS is pressure
C     KNCM is incompresibility
c     RB in fm^-3, KF in fm^-1, EB in MeV, ED  and PRESS in MEV/fm^3
C     LSY and KSY Slope and curvature Symmetry energy CALCULATION

C     ****************************************************

c       WRITE(1,*)RHBO,(RB/(HC*HC*HC)),EB,PRESS
       WRITE(1,*)RHBO,(RB/(HC*HC*HC)),EB,KNCM,PRESS,SYME,LSY,KJN0,KSY

       WRITE(*,*)RHBO,(RB/(HC*HC*HC)),EB,PRESS

      
 10   CONTINUE

         STOP
         END
c-----------------------------------------------------------------
      INCLUDE "parset(untuk-pnm-snm).f"
c------------------------------------------------------------------      
c-----------------------------------------------------------------------------
     
      FUNCTION ksym(xa,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL lsym
      h  = 1.D-4
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h    
      ksym = (lsym(x4,BETA)-8.D0*lsym(x2,BETA)
     . +8.D0*lsym(x1,BETA)-lsym(x3,BETA))/(12.D0*h)
      RETURN
      END


      FUNCTION lsym(xa,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL ASIM
      h  = 1.D-4
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h    
      lsym = (ASIM(x4,BETA)-8.D0*ASIM(x2,BETA)
     . +8.D0*ASIM(x1,BETA)-ASIM(x3,BETA))/(12.D0*h)
      RETURN
      END
C------------------------------------------------------------------------------
C      
C     SYMMETRY ENERGY CALCULATIONS
C
C      
      Function ASIM(RHBO,BETA)     
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      HC  = 197.327D0
C------------------------------------------------------------------------------
C nuclear matter Saturation density !!!

c     BSP B
       KFO =1.303D0*HC
            
C-----------------------------------------------------------------------------
      RHO = 2D0*KFO*KFO*KFO/(3D0*PI*PI)-4D0*BETA*KFO**5/(5D0*PI*PI)  
      R0 = RHBO*RHO  
      KF  = (1.5D0*PI*PI*R0)**(1.D0/3.D0)
     .     +(1.5D0*PI*PI*R0)*2.D0*BETA/5.D0  
              
      KFP= (1.5D0*PI*PI*R0)**(1.D0/3.D0)
     .    +(1.5D0*PI*PI*R0)*2.D0*BETA/5.D0      
      KFN= (1.5D0*PI*PI*R0)**(1.D0/3.D0)   
     .    +(1.5D0*PI*PI*R0)*2.D0*BETA/5.D0

       CALL SET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)



      CALL FRG2(KFP,KFN,SIG,DEL,V0,B0,MN,MP,BETA)
      MBE=(MN+MP)/2.D0

      RHS=KF*DSQRT((KF*KF)+(MBE*MBE))
      V = DLOG((KF+DSQRT((KF*KF)+(MBE*MBE)))/MBE)
      RHS = RHS - MBE*MBE*V
      RHS = RHS * MBE/(PI*PI)
      
      DRHS1 = DSQRT((KF*KF)+(MBE*MBE))
      DRHS = (3.D0*KF*MBE**(5.D0)
     &     +KF**(3.D0)*MBE**(3.D0)-10.D0*KF**(5.D0)*MBE)/DRHS1
      DRHS = DRHS-3.D0*MBE**(5.D0)*V
      DRHS = DRHS/(24.D0*PI*PI)
      
      RHS = RHS+2*BETA*DRHS
      
      MRS = (MR*MR+2.D0*F2*SIG+G4*V0*V0+G3*SIG*SIG)
      EFS = (KF*KF+MBE*MBE)**(1D0/2D0) !ini perlu ditransformasi ga?
      ASIM = KF*KF/(6.D0*EFS)+ GR*GR*KF*KF*KF/(12.D0*PI*PI*MRS)       
      MBR = EFS*EFS*(1+3.D0*GD*GD/(MD*MD)*(RHS/MBE-R0/EFS)) 
      ASIM = ASIM-0.5D0*GD*GD/(MD*MD)*MBE*MBE*R0/MBR

c      WRITE(*,*)(KF/HC),SIG,V0,DEL,MBE
      RETURN
      END  


C------------------------------------------------------------------------------

C     JANOL    
      FUNCTION func3(xa,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL func2
     
      
      h=1.D-4
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h    
      func3=(func2(x4,BETA)-8.D0*func2(x2,BETA)
     . +8.D0*func2(x1,BETA)-func2(x3,BETA))
      func3=func3/(12.D0*h)
      RETURN
      END


C------------------------------------------------------------------------------

C     INCOMPRESSIBILITY CALCULATION    
      FUNCTION func2(xa,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL func1
        
      h=1.D-4
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h    
      func2=(func1(x4,BETA)-8.D0*func1(x2,BETA)
     . +8.D0*func1(x1,BETA)-func1(x3,BETA))
      func2=func2/(12.D0*h)
      RETURN
      END
C------------------------------------------------------------------------------
C     PRESSURE CALCULATION    
      FUNCTION func1(xa,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL func
      h  = 1.D-4
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h    
      func1 = (func(x4,BETA)-8.D0*func(x2,BETA)
     . +8.D0*func(x1,BETA)-func(x3,BETA))/(12.D0*h)
      RETURN
      END
      

      
      FUNCTION func(RHBO,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      CALL FED(RHBO,ED,EB,EL,BETA)
      func=EB
      RETURN
      END

C------------------------------------------------------------------------------
C      

      
C     ENERGY CALCULATIONS      
      SUBROUTINE FED(RHBO,ED,EB,EL,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      HC  = 197.327D0

      print*,RHBO
C------------------------------------------------------------------------------
C nuclear matter Saturation density !!!

c     BSP B
       KFO =1.303D0*HC

      
C-----------------------------------------------------------------------------
    
      RHO = 2.D0*KFO*KFO*KFO/(3.D0*PI*PI)
      DRHO = -4.D0*KFO*KFO*KFO*KFO*KFO/(5.D0*PI*PI)
      RHO = RHO + BETA*DRHO   
      
      RB  = RHBO*RHO  
      CALL SET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)



C     ****************************************************
C      CHEK NUCLEAR MATTER
C     **************************************************** 
c---------------------------------------------------
c  aktifkan untuk kasus Symetric Nuclear Matter(SNM)
      RP = 0.5*RB
      RN = 0.5*RB
c---------------------------------------------------
c  aktifkan untuk kasus Pure Neutron Matter(PNM)
c       RP = 0.0*RB
c       RN = 1.0*RB
c---------------------------------------------------
      RE =0.0*RB
      RM=0.0*RB

      KFE=(3.D0*PI*PI*RE)**(1.D0/3.D0)
     &    +(3.D0*PI*PI*RE)*2.D0*BETA/5.D0      
      KFM=(3.D0*PI*PI*RM)**(1.D0/3.D0)
     &    +(3.D0*PI*PI*RM)*2.D0*BETA/5.D0 
      KFP=(3.D0*PI*PI*RP)**(1.D0/3.D0)
     &    +(3.D0*PI*PI*RP)*2.D0*BETA/5.D0       
      KFN=(3.D0*PI*PI*RN)**(1.D0/3.D0)
     &    +(3.D0*PI*PI*RN)*2.D0*BETA/5.D0 
       EE=0.0D0
       EM=0.0D0
c     **************************************************** 

      

      CALL FRG2(KFP,KFN,SIG,DEL,V0,B0,MN,MP,BETA)

c      WRITE(*,*)RHBO,SIG,DEL,V0  
       RP = KFP*KFP*KFP/(3.D0*PI*PI)-BETA*2.D0*KFP**(5.D0)/(5.D0*PI*PI) 
       RN = KFN*KFN*KFN/(3.D0*PI*PI)-BETA*2.D0*KFN**(5.D0)/(5.D0*PI*PI) 
       
   

      U = 0.5D0*MS*MS*SIG*SIG
      U = U + 0.5D0*MD*MD*DEL*DEL
      U = U - 0.5D0*MV*MV*V0*V0
      U = U - 0.5D0*MR*MR*B0*B0
      U = U + 0.333D0*B2*SIG*SIG*SIG
      U = U + 0.25D0*B3*SIG*SIG*SIG*SIG
      U = U - 0.25D0*C1*V0*V0*V0*V0
      U = U - D2*SIG*V0*V0  
      U = U - F2*SIG*B0*B0
      U = U - 0.5D0*D3*SIG*SIG*V0*V0

      U = U - 0.5D0*G3*SIG*SIG*B0*B0
      U = U - 0.5D0*G4*V0*V0*B0*B0
      U = U - 0.25D0*G5*B0*B0*B0*B0

      EP = KFP*DSQRT((KFP*KFP)+(MP*MP))
      EP = EP*(2.D0*KFP*KFP+MP*MP)
      TEMP = DLOG((KFP+DSQRT((KFP*KFP)+(MP*MP)))/MP)
      EP = EP - MP*MP*MP*MP*TEMP
      EP = EP/(8.D0*PI*PI)
      
      DEP = DSQRT((KFP*KFP)+(MP*MP))
      DEP = DEP*(3.D0*KFP*MP**(4.D0)
     &    -2.D0*KFP**3.D0*MP**(2.D0)-56.D0*KFP**(5.D0))
      DEP = DEP-3.D0*MP**(6.D0)*TEMP
      DEP = DEP/(144.D0*PI*PI)
      
      EP = EP+BETA*DEP
      

      EN = KFN * DSQRT((KFN*KFN)+(MN*MN))
      EN = EN * (2.D0*KFN*KFN+MN*MN)
      TEMN = DLOG((KFN+DSQRT((KFN*KFN)+(MN*MN)))/MN)
      EN = EN - MN*MN*MN*MN*TEMN
      EN = EN/(8.D0*PI*PI)
      
      DEN = DSQRT((KFN*KFN)+(MN*MN))
      DEN = DEN*(3.D0*KFN*MN**(4.D0)
     &      -2.D0*KFN**(3.D0)*MN**(2.D0)-56.D0*KFN**(5.D0))
      DEN = DEN-3.D0*MN**(6.D0)*TEMN
      DEN = DEN/(144.D0*PI*PI)
      
      EN = EN+BETA*DEN
      
     
      EE = KFE * DSQRT((KFE*KFE)+(ME*ME))
      EE = EE * (2.D0*KFE*KFE+ME*ME)
      TEME = DLOG((KFE+DSQRT((KFE*KFE)+(ME*ME)))/ME)
      EE = EE - ME*ME*ME*ME*TEME
      EE = EE /(8.D0*PI*PI)
      
      DEE = DSQRT((KFE*KFE)+(ME*ME))
      DEE = DEE*(3.D0*KFE*ME**(4.D0)
     &    -2.D0*KFE**(3.D0)*ME**(2.D0)-56.D0*KFE**(5.D0))
      DEE = DEE-3.D0*ME**(6.D0)*TEME
      DEE = DEE/(144.D0*PI*PI)
      
      EE = EE+BETA*DEE
            

      EM = KFM * DSQRT((KFM*KFM)+(MO*MO))
      EM = EM * (2.D0*KFM*KFM+MO*MO)
      TEMM = DLOG((KFM+DSQRT((KFM*KFM)+(MO*MO)))/MO)
      EM = EM - MO*MO*MO*MO*TEMM
      EM = EM /(8.D0*PI*PI)
      
      DEM = DSQRT((KFM*KFM)+(MO*MO))
      DEM = DEM*(3.D0*KFM*MO**(4.D0)
     &     -2.D0*KFM**(3.D0)*MO**(2.D0)-56.D0*KFM**(5.D0))
      DEM = DEM-3.D0*MO**(6.D0)*TEMM
      DEM = DEM/(144.D0*PI*PI)
      
      EM = EM+BETA*DEM



      E = EP+EN+EE+EM+GV*V0*(RP+RN)+0.5D0*GR*B0*(RP-RN)


      ED = (E + U )/(HC*HC*HC)
      
      EB = (E + U )/RB

     
      EB = EB - MB

      EL= (EE+EM)/(RB*HC*HC*HC)
c      WRITE(*,*)U,SIG,V0
      RETURN
      END

 
c------------------------------------------------------------------------------ 



      SUBROUTINE FRG2(KFP,KFN,SIG,DEL,V0,B0,MN,MP,BETA)
C
C   meson fields calculation
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      MUE = 1.D+9      
      CALL SET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)  


       RP = KFP*KFP*KFP/(3.D0*PI*PI)-BETA*2.D0*KFP**(5.D0)/(5.D0*PI*PI) 
       RN = KFN*KFN*KFN/(3.D0*PI*PI)-BETA*2.D0*KFN**(5.D0)/(5.D0*PI*PI)   

c      PRESIG = -10.0D0
c      PREV0  = 10.0D0
c      PREDEL = 0.0D0

      OPEN(unit=4,form='formatted',status='old',
     $     file='mesfac')
      READ(4,*)PRESIG,PREV0,PREDEL 
      CLOSE(4)
  
 20   SIG = PRESIG-(1.D-2/MUE)*FS(PRESIG,PREV0,PREDEL,MB,KFP,
     &     KFN,GS,GR,GD,MS,MR,B2,B3,D2,D3,F2,G3,G4,BETA)
      V0 = PREV0-(1.D0/MUE)*FS1(PREV0,PRESIG,MV,MR,GR,GV,KFP,
     &     KFN,C1,D2,D3,F2,G3,G4,BETA)     
      DEL = PREDEL-(1.D0/MUE)* FS4(PREDEL,PRESIG,MD,GS,GD,KFP
     . ,KFN,MB,BETA)

      IF(((SIG-PRESIG)*(SIG-PRESIG)).LT.1.D-20 .AND.
     &     ((DEL-PREDEL)*(DEL-PREDEL)).LT.1.D-20  .AND.  
     &     ((V0-PREV0)*(V0-PREV0)).LT.1.D-20) GOTO 8
      PRESIG = SIG
      PREV0 = V0    
      PREDEL = DEL
c      write (*,*) SIG,V0,DEL
      GOTO 20
 8    CONTINUE 
      B0 = -0.5D0*(RN-RP)*GR/(MR*MR+2.D0*F2*SIG+G4*V0*V0+G3*SIG*SIG)
      MN = MB+GS*SIG-GD*DEL
      MP = MB+GS*SIG+GD*DEL 

       OPEN(unit=6,form='formatted',status='replace',file='mesfac')
       WRITE(6,*)SIG,V0,DEL
       CLOSE(6)
      RETURN
      END

c------------------------------------------------------------------------------

      FUNCTION FS(SIG,V0,DEL,MB,KFP,KFN,GS,GR,GD,MS,MR,B2,B3,D2,
     &            D3,F2,G3,G4,BETA)
C
C  sigma meson calculation
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI = 3.14159265358979D0
      RP = 1.D0*KFP*KFP*KFP/(3.D0*PI*PI)
      RN = 1.D0*KFN*KFN*KFN/(3.D0*PI*PI)
      B0 = -0.5D0*(RN-RP)*GR/(MR*MR+2.D0*F2*SIG+G4*V0*V0+G3*SIG*SIG)  
      MN= MB+GS*SIG-GD*DEL
      MP = MB +GS*SIG+GD*DEL
      FS = MS*MS*SIG
      FS = FS + B2*SIG*SIG
      FS = FS + B3*SIG*SIG*SIG
      FS = FS - D2*V0*V0
      FS = FS - D3*SIG*V0*V0
      FS = FS - F2*B0*B0
      FS = FS - G3*SIG*B0*B0

      VN = KFN*DSQRT(((KFN*KFN)+(MN*MN)))
      VIN = DLOG((KFN+DSQRT((KFN*KFN)+(MN*MN)))/MN)
      VN = VN - MN*MN*VIN
      VN = VN*MN
      VN = VN/(2.D0*PI*PI)
      
      DVN = DSQRT((KFN*KFN)+(MN*MN))
      DVN = (3.D0*KFN*MN**(5.D0)+KFN**(3.D0)*MN**(3.D0)
     &    -10.D0*KFN**(5.D0)*MN)/DVN
      DVN = DVN-3.D0*MN**(5.D0)*VIN
      DVN = DVN/(24.D0*PI*PI)
      
      VN = VN+BETA*DVN
      
      

      VP = KFP*DSQRT(((KFP*KFP)+(MP*MP)))
      VIP = DLOG((KFP+DSQRT((KFP*KFP)+(MP*MP)))/MP)
      VP = VP - MP*MP*VIP
      VP = VP*MP
      VP = VP/(2.D0*PI*PI)
      
      DVP = DSQRT(((KFP*KFP)+(MP*MP)))
      DVP = (3.D0*KFP*MP**(5.D0)+KFP**(3.D0)*MP**(3.D0)
     &    -10.D0*KFP**(5.D0)*MP)/DVN
      DVP = DVP-3.D0*MP**(5.D0)*VIP
      DVP = DVP/(24.D0*PI*PI)
      
      VP = VP+BETA*DVP

      V = (VP+VN)*GS
      FS= FS + V
c      write(*,*)(RP+RN),(VP+VN)
      RETURN
      END
      
      
c------------------------------------------------------------------------------      
      
      FUNCTION FS1(V0,SIG,MV,MR,GR,GV,KFP,KFN,C1,D2,D3,F2,G3,G4,BETA)
C
C omega meson calculation
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      RP = KFP*KFP*KFP/(3.D0*PI*PI)-BETA*2.D0*KFP**(5.D0)/(5D0*PI*PI) 
      RN = KFN*KFN*KFN/(3.D0*PI*PI)-BETA*2.D0*KFN**(5.D0)/(5D0*PI*PI)   
      B0 = -0.5D0*(RN-RP)*GR/(MR*MR+2.D0*F2*SIG+G4*V0*V0+G3*SIG*SIG)


      FS1 = MV*MV*V0
      FS1 = FS1 - GV*(RP+RN)
      FS1 = FS1 + 2.D0*D2*SIG*V0
      FS1 = FS1 + D3*SIG*SIG*V0
      FS1 = FS1 + C1*V0*V0*V0
      FS1 = FS1 + G4*V0*B0*B0
      RETURN
      END


       FUNCTION FS3(B0,V0,SIG,MV,MR,GR,GV,KFP,KFN,C1,D2,D3,F2,G3,
     %             G4,G5,BETA)
C
C rho meson calculation
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      RP = KFP*KFP*KFP/(3.D0*PI*PI)-BETA*2.D0*KFP**(5.D0)/(5.D0*PI*PI) 
      RN = KFN*KFN*KFN/(3.D0*PI*PI)-BETA*2.D0*KFN**(5.D0)/(5.D0*PI*PI)   
 
      FS3 = MR*MR*B0
      FS3 = FS3 - 0.5D0*GR*(RP-RN)
      FS3 = FS3 + 2.D0*F2*SIG*B0
      FS3 = FS3 + F3*SIG*SIG*B0
      FS3 = FS3 + G5*B0*B0*B0
      FS3 = FS3 + G4*V0*V0*B0
      RETURN
      END
c-----------------------------------------------------------------------------

      FUNCTION FS4(DEL,SIG,MD,GS,GD,KFP,KFN,MB,BETA)
C
C  Delta meson calculation
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      MN = MB+GS*SIG-GD*DEL
      MP = MB+GS*SIG+GD*DEL
 
      VN = KFN*DSQRT(((KFN*KFN)+(MN*MN)))
      VIN = DLOG((KFN+DSQRT((KFN*KFN)+(MN*MN)))/MN)
      VN = VN - MN*MN*VIN
      VN = VN*MN
      VN = VN/(2.D0*PI*PI)
      
      DVN = DSQRT((KFN*KFN)+(MN*MN))
      DVN = (3.D0*KFN*MN**(5.D0)+KFN**(3.D0)*MN**(3.D0)
     &       -10.D0*KFN**(5.D0)*MN)/DVN
      DVN = DVN-3.D0*MN**(5.D0)*VIN
      DVN = DVN/(24.D0*PI*PI)
      
      VN = VN+BETA*DVN
      

      VP = KFP*DSQRT(((KFP*KFP)+(MP*MP)))
      VIP = DLOG((KFP+DSQRT((KFP*KFP)+(MP*MP)))/MP)
      VP = VP - MP*MP*VIP
      VP = VP*MP
      VP = VP/(2.D0*PI*PI)
      
      DVP = DSQRT(((KFP*KFP)+(MP*MP)))
      DVP = (3.D0*KFP*MP**(5.D0)+KFP**(3.D0)*MP**(3.D0)
     &     -10.D0*KFP**(5.D0)*MP)/DVN
      DVP = DVP-3.D0*MP**(5.D0)*VIP
      DVP = DVP/(24.D0*PI*PI)
      
      VP = VP+BETA*DVP
 
      V = (VP-VN)*GD
      FS4 = MD*MD*DEL+V 
      RETURN
      END
c-----------------------------------------------------------------------------
 
