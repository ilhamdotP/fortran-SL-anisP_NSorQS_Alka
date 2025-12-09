      PROGRAM NSM        
c     ****************************************************
c     verified 5 March 2022 check "new analitic pressure"
c     ****************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      INTEGER I
      CHARACTER(LEN = 1000) filename, frmt
      PI  = 3.14159265358979D0
      HC  = 197.327D0
      frmt = "(A15,ES6.0,A4)"
      
      
      
      BETA = 2.D-7
      

C------------------------------------------------------------------------------
C nuclear matter Saturation density !!!

c     FSU 
       KFO =1.306D0*HC

C------------------------------------------------------------------------------
c      RHO = 2D0*KFO*KFO*KFO/(3D0*PI*PI)
      
      RHO = 2D0*KFO*KFO*KFO/(3D0*PI*PI)-4D0*BETA*KFO**5/(5D0*PI*PI)
 
     
      write(filename,frmt)"EOS_BSPB_GUP", BETA, ".dat"
c      write(*,*)filename
      OPEN(unit=1,file=trim(filename),status='unknown')
c      OPEN(unit=2,status='unknown',file='SOS_NSM_FSU2HZ2.dat')




c     DO 10 I=1,85
      DO 10 I=1,850
        RHBO=I * 0.01D0  
         write(*,*)I    
         RB=RHBO*RHO    


C     CALL ENERGY DENSITY AND BINDING ENERGY
C     

         CALL FED(RHBO,ED,EB,EL,PT,PAVG,PRAD,PTAN,BETA)
 
c     CALL other properties
c         CALL FERMIM(RB,KFP,KFN,KFE,KFM,YE,YN,YP,YM)
c         CALL FRG2(KFP,KFN,SIG,DEL,V0,B0,MN,MP)
c         CALL SPOSN(RHBO,DPDE,CL,SOS,ED,PRESSCR,BETA)
c         DEDEN=DE(RHBO)
c         DPRESS=DP(RHBO)
C     
C     CALL PRESSURE
C     
         PRESS= PT
c         ABC=func1(RHBO)

c         PRESS= RHO*RHBO*RHBO*func1(RHBO,BETA)/(HC*HC*HC)
c         KNCM= 9.0D0*RHBO*RHBO*func2(RHBO,BETA)
c         KJN0= 27.0D0*RHBO*RHBO*RHBO*func3(RHBO,BETA)
C    
c         SYME= ASIM(RHBO,BETA)
c         LSY= 3.0D0*RHBO*lsym(RHBO,BETA)
c         KSY= 9.0D0*RHBO*RHBO*ksym(RHBO,BETA)
         
c         KASY=KSY-6.0D0*LSY
c         KSAT2=KASY-KJN0/KNCM*LSY



C     RESULTS
C     
C     ****************************************************
C     RB is Baryon density,  EB is binding energy, ED is energy density
C     PRESS is pressure
C     KNCM is incompresibility
c     RB in fm^-3, KF in fm^-1, EB in MeV, ED  and PRESS in MEV/fm^3

C     ****************************************************
       WRITE(1,*)RHBO,(RB/(HC*HC*HC)),ED,PRESS,PRAD,PTAN
c       WRITE(2,*)RHBO,SOS,CL,PRESSCR,ED
c       WRITE(*,*)RHBO,CL,SOS,ED,PRESSCR,PRESS
c       WRITE(*,*)RHBO,PT,PAVG,PRAD,PTAN
       WRITE(*,*)RHBO,PRAD,PTAN
      
 10   CONTINUE

         STOP
         END
c-----------------------------------------------------------------
c      BSP B in parset(untuk-pnm-snm).f
c      BKA22 in parset(untuk-nsm).f
      INCLUDE "parset(untuk-pnm-snm).f" 
c      INCLUDE "parset(untuk-nsm).f" 
c------------------------------------------------------------------      
      FUNCTION gunc1(xa)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      CALL FED(xa,ED,EB,EL,PT,PAVG,PRAD,PTAN,BETA)
      gunc1=ED
      RETURN
      END
      FUNCTION gunc2(xa)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)     
       CALL FED(xa,ED,EB,EL,PT,PAVG,PRAD,PTAN,BETA)
      gunc2=PT
      RETURN
      END

      FUNCTION DE(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL gunc1
      h  = 1.D-4
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h    
      DE = (gunc1(x4)-8.D0*gunc1(x2)+8.D0*gunc1(x1)-gunc1(x3))/(12.D0*h)
      RETURN
      END

      FUNCTION DP(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL gunc2
      h  = 1.D-4
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h    
      DP = (gunc2(x4)-8.D0*gunc2(x2)+8.D0*gunc2(x1)-gunc2(x3))/(12.D0*h)
      RETURN
      END
      
      SUBROUTINE SPOSN(RHBO,DPDE,CL,SOS,ED,PRESSCR,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL DE,DP,gunc2
      CL=1.D0/DSQRT(3.D0)
      DPRESS=DP(RHBO)
      DEDEN=DE(RHBO)
      DPDE=DPRESS/DEDEN
      
      IF (DPDE .GE. 0.0D0) THEN 
         SOS=DSQRT(DPDE)
      ELSE
         SOS=0.0D0
      ENDIF
c  Note: be carefull here token SOS=0 means unstable (imaginer speed of sound)!!!         
      IF (DSQRT((SOS-CL)*(SOS-CL)) .LT. 2.0D-3) THEN
         CALL FED(RHBO,ED,EB,EL,PT,PAVG,PRAD,PTAN,BETA)
         EDX=ED
         PRX=gunc2(RHBO)
      END IF
      
      IF (SOS .GT. CL ) THEN
         PRESSCR=CL*(ED-EDX)+PRX
      ELSE   
         PRESSCR=gunc2(RHBO)
      END IF


      RETURN
      END
c------------------------------------------------------------------------------


C------------------------------------------------------------------------------
C      

C------------------------------------------------------------------------------
C     PRESSURE CALCULATION    
      FUNCTION func1(xa,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL func
      h  = 1.D-8
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
      CALL FED(RHBO,ED,EB,EL,PT,PAVG,PRAD,PTAN,BETA)
      func=EB
      RETURN
      END

C------------------------------------------------------------------------------
C

      
C     ENERGY CALCULATIONS      
      SUBROUTINE FED(RHBO,ED,EB,EL,PT,PAVG,PRAD,PTAN,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      HC  = 197.327D0


C------------------------------------------------------------------------------
C nuclear matter Saturation density !!!

c     FSU 
       KFO =1.306D0*HC

C-----------------------------------------------------------------------------
    
      RHO = 2.D0*KFO*KFO*KFO/(3.D0*PI*PI) 
      DRHO = -4.D0*KFO*KFO*KFO*KFO*KFO/(5.D0*PI*PI)
      RHO = RHO + BETA*DRHO   
      
      RB  = RHBO*RHO  
      CALL SET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)
c      write(*,*) "call parset in FED"
c      PRINT*, MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
c     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD
c      STOP


C     ****************************************************
C      CHEK NUCLEAR MATTER
C     **************************************************** 
c      RP = 0.0*RB
c      RN = 1.0*RB
c      RE =0.0*RB
c      RM=0.0*RB

c      KFE=(3.D0*PI*PI*RE)**(1.D0/3.D0)
c      KFM=(3.D0*PI*PI*RM)**(1.D0/3.D0)
c      KFP=(3.D0*PI*PI*RP)**(1.D0/3.D0)
c      KFN=(3.D0*PI*PI*RN)**(1.D0/3.D0)
c       EE=0.0D0
c       EM=0.0D0
c     **************************************************** 

      CALL FERMIM(RB,KFP,KFN,KFE,KFM,YE,YN,YP,YM,BETA)

      CALL FRG2(KFP,KFN,SIG,DEL,V0,B0,MN,MP,BETA)

c      RP = KFP*KFP*KFP/(3.D0*PI*PI)
c      RN = KFN*KFN*KFN/(3.D0*PI*PI)   
      RP = KFP*KFP*KFP/(3.D0*PI*PI)-BETA*2D0*KFP**5/(5D0*PI*PI) 
      RN = KFN*KFN*KFN/(3.D0*PI*PI)-BETA*2D0*KFN**5/(5D0*PI*PI)     
      

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
c---------------------------------------------------------------------------
c     should be modified      
      U = U - 0.5D0*G3*SIG*SIG*B0*B0
      U = U - 0.5D0*G4*V0*V0*B0*B0
      U = U - 0.25D0*G5*B0*B0*B0*B0
c-------------------------------------------------------------------------
      EP = KFP*DSQRT((KFP*KFP)+(MP*MP))
      EP = EP*(2.D0*KFP*KFP+MP*MP)
      TEMP = DLOG((KFP+DSQRT((KFP*KFP)+(MP*MP)))/MP)
      EP = EP - MP*MP*MP*MP*TEMP
      EP = EP/(8.D0*PI*PI)

      DEP = DSQRT((KFP*KFP)+(MP*MP))
      DEP = DEP*(3.D0*KFP*MP**4-2*KFP**3*MP**2-56.D0*KFP**5)
      DEP = DEP-3*MP**6*TEMP
      DEP = DEP/(144.D0*PI*PI)
      
      EP = EP+BETA*DEP
      
      EN = KFN * DSQRT((KFN*KFN)+(MN*MN))
      EN = EN * (2.D0*KFN*KFN+MN*MN)
      TEMN = DLOG((KFN+DSQRT((KFN*KFN)+(MN*MN)))/MN)
      EN = EN - MN*MN*MN*MN*TEMN
      EN = EN/(8.D0*PI*PI)
      
      DEN = DSQRT((KFN*KFN)+(MN*MN))
      DEN = DEN*(3.D0*KFN*MN**4-2*KFN**3*MN**2-56.D0*KFN**5)
      DEN = DEN-3*MN**6*TEMN
      DEN = DEN/(144.D0*PI*PI)
      
      EN = EN+BETA*DEN


      EE = KFE * DSQRT((KFE*KFE)+(ME*ME))
      EE = EE * (2.D0*KFE*KFE+ME*ME)
      TEME = DLOG((KFE+DSQRT((KFE*KFE)+(ME*ME)))/ME)
      EE = EE - ME*ME*ME*ME*TEME
      EE = EE /(8.D0*PI*PI)
      
      DEE = DSQRT((KFE*KFE)+(ME*ME))
      DEE = DEE*(3.D0*KFE*ME**4-2*KFE**3*ME**2-56.D0*KFE**5)
      DEE = DEE-3*ME**6*TEME
      DEE = DEE/(144.D0*PI*PI)
      
      EE = EE+BETA*DEE

      EM = KFM * DSQRT((KFM*KFM)+(MO*MO))
      EM = EM * (2.D0*KFM*KFM+MO*MO)
      TEMM = DLOG((KFM+DSQRT((KFM*KFM)+(MO*MO)))/MO)
      EM = EM - MO*MO*MO*MO*TEMM
      EM = EM /(8.D0*PI*PI)
      
      DEM = DSQRT((KFM*KFM)+(MO*MO))
      DEM = DEM*(3.D0*KFM*MO**4-2*KFM**3*MO**2-56.D0*KFM**5)
      DEM = DEM-3*MO**6*TEMM
      DEM = DEM/(144.D0*PI*PI)
      
      EM = EM+BETA*DEM



      E = EP+EN+EE+EM+GV*V0*(RP+RN)+0.5D0*GR*B0*(RP-RN)

      ED = (E + U )/(HC*HC*HC)
      
      EB = (E + U )/RB

     
      EB = EB - MB

      EL= (EE+EM)/(RB*HC*HC*HC)

c------------------------------------------------------
c     PRESSURE Analitic
c-----------------------------------------------------      

      PP = KFP*DSQRT((KFP*KFP)+(MP*MP))
      PP = PP*(2.D0*KFP*KFP-3.0D0*MP*MP)
      TEMP2 = DLOG((KFP+DSQRT((KFP*KFP)+(MP*MP)))/MP)
      PP = PP + 3.0D0*MP*MP*MP*MP*TEMP2
      PP = PP/(8.D0*PI*PI)   
      
      DPP = (15.*KFP*MP**6+5.*KFP**3*MP**4-2.*KFP**5*MP**2
     .     -40.*KFP**7)/DSQRT((KFP*KFP)+(MP*MP))
      DPP = DPP-15*MP**6*TEMP2
      DPP = DPP/(144.D0*PI*PI)

      PP = PP+BETA*DPP

      DMUP = - KFP**4/(3*DSQRT((KFP*KFP)+(MP*MP)))
      MUP  = DSQRT((KFP*KFP)+(MP*MP)) + GV*V0 + 0.5D0*GR*B0
      RHOP = KFP*KFP*KFP/(3.D0*PI*PI)
      DRHOP= -2D0*KFP**5/(5D0*PI*PI)
      SBRP = -2./3.*(DMUP*RHOP+MUP*DRHOP-(DEP+DPP/3))*BETA

      PN = KFN*DSQRT((KFN*KFN)+(MN*MN))
      PN = PN*(2.D0*KFN*KFN-3.0D0*MN*MN)
      TEMN2 = DLOG((KFN+DSQRT((KFN*KFN)+(MN*MN)))/MN)
      PN = PN + 3.0D0*MN*MN*MN*MN*TEMN2
      PN = PN/(8.D0*PI*PI)
      
      DPN = (15.*KFN*MN**6+5.*KFN**3*MN**4-2.*KFN**5*MN**2
     .     -40.*KFN**7)/DSQRT((KFN*KFN)+(MN*MN))
      DPN = DPN-15*MN**6*TEMN2
      DPN = DPN/(144.D0*PI*PI)

      PN = PN+BETA*DPN
      
      DMUN = - KFN**4/(3*DSQRT((KFN*KFN)+(MN*MN)))
      MUN  = DSQRT((KFN*KFN)+(MN*MN)) + GV*V0 - 0.5D0*GR*B0
      RHON = KFN*KFN*KFN/(3.D0*PI*PI)
      DRHON= -2D0*KFN**5/(5D0*PI*PI)
      SBRN = -2./3.*(DMUN*RHON+MUN*DRHON-(DEN+DPN/3))*BETA    
          
      
      PE = KFE * DSQRT((KFE*KFE)+(ME*ME))
      PE = PE * (2.D0*KFE*KFE-3.0D0*ME*ME)
      TEME2 = DLOG((KFE+DSQRT((KFE*KFE)+(ME*ME)))/ME)
      PE = PE + 3.0D0*ME*ME*ME*ME*TEME2
      PE = PE /(8.D0*PI*PI)

      DPE = (15.*KFE*ME**6+5.*KFE**3*ME**4-2.*KFE**5*ME**2
     .     -40.*KFE**7)/DSQRT((KFE*KFE)+(ME*ME))
      DPE = DPE-15*ME**6*TEME2
      DPE = DPE/(144.D0*PI*PI)

      PE = PE+BETA*DPE

      DMUE = - KFE**4/(3*DSQRT((KFE*KFE)+(ME*ME)))
      MUE  = DSQRT((KFE*KFE)+(ME*ME))
      RHOE = KFE*KFE*KFE/(3.D0*PI*PI)
      DRHOE= -2D0*KFE**5/(5D0*PI*PI)
      SBRE = -2./3.*(DMUE*RHOE+MUE*DRHOE-(DEE+DPE/3))*BETA

      PM = KFM * DSQRT((KFM*KFM)+(MO*MO))
      PM = PM * (2.D0*KFM*KFM-3.0D0*MO*MO)
      TEMM2 = DLOG((KFM+DSQRT((KFM*KFM)+(MO*MO)))/MO)
      PM = PM + 3.0D0*MO*MO*MO*MO*TEMM2
      PM = PM /(8.D0*PI*PI)
      
      DPM = (15.*KFM*MO**6+5.*KFM**3*MO**4-2.*KFM**5*MO**2
     .     -40.*KFM**7)/DSQRT((KFM*KFM)+(MO*MO))
      DPM = DPM-15*MO**6*TEMM2
      DPM = DPM/(144.D0*PI*PI)

      PM = PM+BETA*DPM      

      DMUM = - KFM**4/(3*DSQRT((KFM*KFM)+(MO*MO)))
      MUM  = DSQRT((KFM*KFM)+(MO*MO))
      RHOM = KFM*KFM*KFM/(3.D0*PI*PI)
      DRHOM= -2D0*KFM**5/(5D0*PI*PI)
      SBRM = -2./3.*(DMUM*RHOM+MUM*DRHOM-(DEM+DPM/3))*BETA


      P=(PP+PN+PE+PM)/3.0D0
      SBR = (SBRP+SBRN+SBRE+SBRM)
      
      PRADS = P + 2.*SBR/3.
      PTANS = PRADS - SBR
      
C      write(*,*)P,SBR!P,SBRN,SBRE,SBRM,"P,SBR"
      
      PAVGS = (PRADS+2*PTANS)/3. 
      
      PT = (P - U )/(HC*HC*HC)
      PRAD = (PRADS - U )/(HC*HC*HC)
      PTAN = (PTANS - U )/(HC*HC*HC)
      PAVG = (PAVGS - U )/(HC*HC*HC)      

      RETURN
      END

 
c------------------------------------------------------------------------------ 

      SUBROUTINE  FERMIM(RB,KFPF,KFNF,KFEF,KFMF,YE,YN,YP,YM,BETA)
C
C   Fermi momenta calculations
C

      IMPLICIT DOUBLE PRECISION (A-H,K-Z)       

      PI  = 3.14159265358979D0
      MUE = 100.D0      
      CALL SET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)  
C      write(*,*) "in FERMIM"     
       
       YN0=0.8D0
       YP0=0.2D0
 
 50    RBN = YN0*RB
         RBP = YP0*RB
         KFP =(3.D0*PI*PI*RBP)**(1D0/3D0)
     .    +(3.D0*PI*PI*RBP)*2.D0*BETA/5.D0
         KFN =(3.D0*PI*PI*RBN)**(1D0/3D0)
     .    +(3.D0*PI*PI*RBN)*2.D0*BETA/5.D0 
     
     
         CALL FRG2(KFP,KFN,SIG,DEL,V0,B0,MN,MP,BETA)   
C         PRINT*, "FINISH CALLING FRG2 IN FERMIM"
        
         PREYP = 0.0D0
 20      YP = PREYP-(1.D0/MUE)* 
     .      FYP(PREYP,MN,MP,RB,ME,MO,SIG,DEL,V0,B0,BETA)
C         print*, YP, PREYP
         IF (  YP .ne. YP ) THEN
C            PRINT*, 'YP is NaN!'
C            print*, PREYP,(YP-PREYP)**2
            GOTO 91
         ENDIF
         IF(((YP-PREYP)*(YP-PREYP)).LT.1.D-24) GOTO 9
         PREYP = YP
         GOTO 20
         
 91      YP=PREYP   
C         STOP 
         
 9       CONTINUE
C         print*, "in FERMIM:",YP
             
C        
C     Fraction of each contituent
C
         A = (3.D0*PI*PI*(1-YP))**(2D0/3D0)
         B = (3.D0*PI*PI*YP)**(2D0/3D0)
         CN = (MN*MN)/RB**(2D0/3D0)
         CP = (MP*MP)/RB**(2D0/3D0)
         D = (ME*ME)/RB**(2D0/3D0)
         F = (MO*MO)/RB**(2D0/3D0)        
         G = 1.D0/(3.D0*PI*PI)
         MZ= (MR*MR+2.D0*F2*SIG+G4*V0*V0+G3*SIG*SIG+G5*B0*B0)
         R = 0.5D0*GR*GR*RB**(2.D0/3.D0)*(1.D0-2.D0*YP)/MZ
C------------------------------------------------------------------
C GT = GUP Term
         GT = (B*B/DSQRT(B+CP)-A*A/DSQRT(B+CN))
     .        *2./15.*BETA*RB**(2./3.)
C------------------------------------------------------------------     
         FZ = (DSQRT(A+CN)-DSQRT(B+CP)+R+GT)**2
         YN = (1.D0-YP)
         YE = G*(FZ-D)**(3D0/2D0)
     .        *(1.+(FZ-D)*4.*BETA/5.)**(3D0/2D0)
         IF ((FZ-F) .GT. 0D0) THEN
         YM = G*(FZ-F)**(3D0/2D0)
     .        *(1.+(FZ-F)*4.*BETA/5.)**(3D0/2D0)
         ELSE 
            YM = 0D0
         ENDIF  
    
        
        IF (  YP .ne. YP ) THEN
            PRINT*, 'YP is NaN! again... T__T'
            print*,YP,(YP-YP0)**2
c            STOP
            GOTO 41
         ENDIF
       IF(((YP-YP0)*(YP-YP0)).LT.1.D-6 ) GOTO 40
        YP0=YP
        YN0=YN
       GOTO 50
 41    YP=YP0
 40    CONTINUE
 
       KFEF= (3.D0*YE*RB*PI*PI)**(1D0/3D0)+2.*BETA*(3.*YE*RB*PI*PI)/5.
       KFPF= (3.D0*YP*RB*PI*PI)**(1D0/3D0)+2.*BETA*(3.*YP*RB*PI*PI)/5.
       KFNF= (3.D0*YN*RB*PI*PI)**(1D0/3D0)+2.*BETA*(3.*YN*RB*PI*PI)/5.
       KFMF= (3.D0*YM*RB*PI*PI)**(1D0/3D0)+2.*BETA*(3.*YM*RB*PI*PI)/5.
c      WRITE(*,*)YP,YN,YE,YM,B0
      RETURN
      END
c-------------------------------------------------------------------

      SUBROUTINE FRG2(KFP,KFN,SIG,DEL,V0,B0,MN,MP,BETA)
C
C   meson fields calculation
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      MUE = 1.D+9      
      CALL SET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)  
C      write(*,*) "in FRG2"
       

       RP = KFP*KFP*KFP/(3.D0*PI*PI)-BETA*2D0*KFP**5/(5D0*PI*PI) 
       RN = KFN*KFN*KFN/(3.D0*PI*PI)-BETA*2D0*KFN**5/(5D0*PI*PI)      

C     initial guess
       
c      PREDEL = (RP-RN)*GD/(MD*MD)
c      PRESIG = (RN+RP)*GS/(MS*MS)
c      PREV0  = (RN+RP)*GV/(MV*MV)
      
c     MR2EF=MR*MR+2.D0*F2*PRESIG+G4*PREV0*PREV0+G3*PRESIG*PRESIG
c      MR2EF=MR*MR
c      PREB0  = 0.5D0*(RP-RN)*GR/MR2EF

      PREDEL = 0.0D0
      PRESIG = 0.0D0
      PREV0  = 0.0D0
      PREB0  = 0.0D0
      
 20   SIG = PRESIG-(1.D0/MUE)*FS(PRESIG,PREV0,PREDEL,PREB0,MB,KFP,
     &     KFN,GS,GR,GD,MS,MR,B2,B3,D2,D3,F2,G3,G4,BETA)
      V0 = PREV0-(1.D0/MUE)*FS1(PREV0,PRESIG,PREB0,MV,MR,GR,GV,KFP,
     &     KFN,C1,D2,D3,F2,G3,G4,BETA)     
      B0 = PREB0-(1.D0/MUE)*FS3(PREB0,PREV0,PRESIG,MV,MR,GR,GV,
     &     KFP,KFN,F2,G3,G4,G5,BETA)
      DEL = PREDEL-(1.D0/MUE)*FS4(PREDEL,PRESIG,MD,GS,GD,KFP,KFN,MB,
     &     BETA)
C      write(*,*) "looping",(SIG-PRESIG)**2,(DEL-PREDEL)**2,
C     &      (B0-PREB0)**2,(V0-PREV0)**2
        IF (  SIG .ne. SIG ) THEN
C            PRINT*, 'SiG is NaN!'
C            print*, (SIG-PRESIG)**2,(DEL-PREDEL)**2,
C     &      (B0-PREB0)**2,(V0-PREV0)**2
            GOTO 81
        ELSEIF (  DEL .ne. DEL ) THEN
            PRINT*, 'DEL is NaN!'
            print*, (SIG-PRESIG)**2,(DEL-PREDEL)**2,
     &      (B0-PREB0)**2,(V0-PREV0)**2
            GOTO 81
        ELSEIF (  B0 .ne. B0 ) THEN
            PRINT*, 'B0 is NaN!'
            print*, (SIG-PRESIG)**2,(DEL-PREDEL)**2,
     &      (B0-PREB0)**2,(V0-PREV0)**2
            GOTO 81
        ELSEIF (  V0 .ne. V0 ) THEN
            PRINT*, 'V0 is NaN!'
            print*, (SIG-PRESIG)**2,(DEL-PREDEL)**2,
     &      (B0-PREB0)**2,(V0-PREV0)**2
            GOTO 81
        ENDIF
      IF(((SIG-PRESIG)*(SIG-PRESIG)).LT.1.D-28 .AND.
     &   ((DEL-PREDEL)*(DEL-PREDEL)).LT.1.D-28 .AND.
     &   ((B0-PREB0)*(B0-PREB0)).LT.1.D-28     .AND. 
     &   ((V0-PREV0)*(V0-PREV0)).LT.1.D-28   ) GOTO 8

      PRESIG = SIG
      PREV0 = V0
      PREB0= B0
      PREDEL = DEL      
      GOTO 20
      
 81   SIG = PRESIG
      V0 = PREV0
      B0= PREB0
      DEL = PREDEL 
      
C      print*, "in FRG2:",SIG,DEL,V0,B0
      
      
 8    CONTINUE 
               
      MN = MB+GS*SIG-GD*DEL
      MP = MB+GS*SIG+GD*DEL 

      RETURN
      END

c------------------------------------------------------------------------------

      FUNCTION FS(SIG,V0,DEL,B0,MB,KFP,KFN,GS,GR,GD,MS,MR,B2,B3,D2,
     &            D3,F2,G3,G4,BETA)
C
C  sigma meson calculation
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI = 3.14159265358979D0
      RP = 1.D0*KFP*KFP*KFP/(3.D0*PI*PI)-2D0*BETA*KFP**5/(5D0*PI*PI)
      RN = 1.D0*KFN*KFN*KFN/(3.D0*PI*PI)-2D0*BETA*KFN**5/(5D0*PI*PI)
    
      MN = MB + GS*SIG-GD*DEL
      MP = MB + GS*SIG+GD*DEL
      FS = MS*MS*SIG
      FS = FS + B2*SIG*SIG
      FS = FS + B3*SIG*SIG*SIG
      FS = FS - D2*V0*V0
      FS = FS - D3*SIG*V0*V0
c-------------------------------------------------------------------------
c  should be modifed      
      FS = FS - F2*B0*B0
      FS = FS - G3*SIG*B0*B0
c------------------------------------------------------------------------
      VN = KFN*DSQRT(((KFN*KFN)+(MN*MN)))
      VIN = DLOG((KFN+DSQRT((KFN*KFN)+(MN*MN)))/MN)
      VN = VN - MN*MN*VIN
      VN = VN*MN
      VN = VN/(2.D0*PI*PI)
      
      DVN = DSQRT((KFN*KFN)+(MN*MN))
      DVN = (3.D0*KFN*MN**5+KFN**3*MN**3-10.D0*KFN**5*MN)/DVN
      DVN = DVN-3.D0*MN**5*VIN
      DVN = DVN/(24.D0*PI*PI)
      
      VN = VN+BETA*DVN
      
      
      VP = KFP*DSQRT(((KFP*KFP)+(MP*MP)))
      VIP = DLOG((KFP+DSQRT((KFP*KFP)+(MP*MP)))/MP)
      VP = VP - MP*MP*VIP
      VP = VP*MP
      VP = VP/(2.D0*PI*PI)
      
      DVP = DSQRT(((KFP*KFP)+(MP*MP)))
      DVP = (3.D0*KFP*MP**5+KFP**3*MP**3-10.D0*KFP**5*MP)/DVN
      DVP = DVP-3.D0*MP**5*VIP
      DVP = DVP/(24.D0*PI*PI)
      
      VP = VP+BETA*DVP

      V = (VP+VN)*GS
      FS= FS + V

      RETURN
      END
      
      
c------------------------------------------------------------------------------      
      
      FUNCTION FS1(V0,SIG,B0,MV,MR,GR,GV,KFP,KFN,C1,D2,D3,F2,G3,G4,
     &             BETA)
C
C omega meson calculation
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      RP = 1.D0*KFP*KFP*KFP/(3.D0*PI*PI)-2D0*BETA*KFP**5/(5D0*PI*PI)
      RN = 1.D0*KFN*KFN*KFN/(3.D0*PI*PI)-2D0*BETA*KFN**5/(5D0*PI*PI)

      FS1 = MV*MV*V0
      FS1 = FS1 - GV*(RP+RN)
      FS1 = FS1 + 2.D0*D2*SIG*V0
      FS1 = FS1 + D3*SIG*SIG*V0
c--------------------------------------------------------------------
c should be modifed      
      FS1 = FS1 + C1*V0*V0*V0
      FS1 = FS1 + G4*V0*B0*B0
c-------------------------------------------------------------------      
      RETURN
      END


       FUNCTION FS3(B0,V0,SIG,MV,MR,GR,GV,KFP,KFN,F2,G3,G4,G5,BETA)
C
C rho meson calculation
C
c--------------------------------------------------------------------
c     should be modifed
       
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      RP = KFP*KFP*KFP/(3.D0*PI*PI)-BETA*2D0*KFP**5/(5D0*PI*PI) 
      RN = KFN*KFN*KFN/(3.D0*PI*PI)-BETA*2D0*KFN**5/(5D0*PI*PI) 
 
      FS3 = MR*MR*B0
      FS3 = FS3 - 0.5D0*GR*(RP-RN)
      FS3 = FS3 + 2.D0*F2*SIG*B0
      FS3 = FS3 + G3*SIG*SIG*B0
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
      DVN = (3.D0*KFN*MN**5+KFN**3*MN**3-10.D0*KFN**5*MN)/DVN
      DVN = DVN-3.D0*MN**5*VIN
      DVN = DVN/(24.D0*PI*PI)
      
      VN = VN+BETA*DVN
      

      VP = KFP*DSQRT(((KFP*KFP)+(MP*MP)))
      VIP = DLOG((KFP+DSQRT((KFP*KFP)+(MP*MP)))/MP)
      VP = VP - MP*MP*VIP
      VP = VP*MP
      VP = VP/(2.D0*PI*PI)
      
      DVP = DSQRT(((KFP*KFP)+(MP*MP)))
      DVP = (3.D0*KFP*MP**5+KFP**3*MP**3-10.D0*KFP**5*MP)/DVN
      DVP = DVP-3.D0*MP**5*VIP
      DVP = DVP/(24.D0*PI*PI)
      
      VP = VP+BETA*DVP
 
      V = (VP-VN)*GD
      FS4 = MD*MD*DEL+V 
      RETURN
      END
c-----------------------------------------------------------------------------
 
      FUNCTION FYP(YP,MN,MP,RB,ME,MO,SIG,DEL,V0,B0,BETA)
C
C   Function of fraction
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      CALL SET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)
 
      A = (3.D0*PI*PI*(1.D0-YP))**(2D0/3D0)
      B = (3.D0*PI*PI*YP)**(2D0/3D0)
      CN = (MN*MN)/RB**(2D0/3D0)
      CP = (MP*MP)/RB**(2D0/3D0)
      D = (ME*ME)/RB**(2D0/3D0)
      F = (MO*MO)/RB**(2D0/3D0)    
      G = 1.D0/(3.D0*PI*PI)
      MZ= (MR*MR+2.D0*F2*SIG+G4*V0*V0+G3*SIG*SIG+G5*B0*B0)
      R = 0.5D0*GR*GR*RB**(2.D0/3.D0)*(1.D0-2.D0*YP)/MZ   
C------------------------------------------------------------------
C GT = GUP Term

      GT = (B*B/DSQRT(B+CP)-A*A/DSQRT(B+CN))*2./15.*BETA*RB**(2./3.)
C------------------------------------------------------------------           
      H = (DSQRT(A+CN)-DSQRT(B+CP)+R+GT)**2
      IF ((H-F).GE.0.D0) THEN 
      GF = 1.+4.*BETA*(H-D)/15.
         FYP = -G**(2./3.)*(((H-D)*(1.+4.*BETA*(H-D)/15.))**(3D0/2D0)
     &         +((H-F)*1.+4./15.*BETA*(H-F))**(3D0/2D0))**(2D0/3D0)
     &         +YP**(2D0/3D0)
      ELSE     
         FYP =-G**(2./3.)*((H-D)*(1.+4.*BETA*(H-D)/15.))**(3D0/2D0)
     &         +YP**(2D0/3D0)
      ENDIF
      RETURN
      END

c-----------------------------------------------------------------------------
