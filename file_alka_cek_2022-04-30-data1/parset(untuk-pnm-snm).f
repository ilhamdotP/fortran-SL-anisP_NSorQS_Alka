c
c     BKA22 and BSR12
c
C--------------------------------------------------------------
C    BKA 22

      SUBROUTINE ASET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
C     Initialisation
      PI  = 3.14159265358979D0
      HC  = 197.327D0
      ME  = 0.511D0
      MO  = 105.6D0

      GS = 0.8462D0*4.D0*PI
      GV = 1.1089D0*4.D0*PI
c----------------------------------------------------------
      GR = 1.0302D0*4.D0*PI
      GD = 0.0D0
      ETR = -3.9294D0
c---------------------------------------------------------- 
      K2  = -1.5500D0
      K3  =  2.13451D0
      ET2 = -0.1555D0
      ET3 = 0.0697D0

      ETR3 =  0.0D0
      ETR4 = 0.0D0

      EP3 = 5.8253D0
      EPR3= 0.0D0

      MB  = 939.0D0
  
 
      MS  =  0.5302D0*MB
      MV  =  782.D0
      MR  =  770.D0
      MD  =  980.D0

   
  
      B2 = GS*MS*MS*K2/(2.D0*MB)
      B3 = GS*GS*MS*MS*K3/(6.D0*MB*MB)
      C1 = GV*GV*EP3/(6.D0)
      D2 = GS*MV*MV*ET2/(2.D0*MB)
      D3 = GS*GS*MV*MV*ET3/(2.D0*MB*MB)  
      F2 = ETR*GS*MR*MR/(2.D0*MB)
      G3 = GS*GS*MR*MR*ETR3/(2.D0*MB*MB)
      G4 = GV*GV*MR*MR*ETR4/(2.D0*MB*MB)
      G5 = GR*GR*EPR3/(6.D0)     


c--------------------------------------------------------------
      RETURN
      END 
c---------------------------------------------------------------
c  BSR12


      SUBROUTINE DSET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
C     Initialisation
      PI  = 3.14159265358979D0
      HC  = 197.327D0
      ME  = 0.511D0
      MO  = 105.6D0

      GS = 0.8450D0*4.D0*PI
      GV = 1.1051D0*4.D0*PI
c----------------------------------------------------------
      GR = 0.8725D0*4.D0*PI
      GD = 0.0D0
      ETR = -0.8923D0
c---------------------------------------------------------- 
      K2  = -1.4858D0
      K3  =  1.9261D0
      ET2 = -0.1469D0
      ET3 = 0.0029D0

      ETR3 =  4.5617D0
      ETR4 = 1.036D0

      EP3 = 5.7855D0
      EPR3= 0.0D0

      MB  = 939.0D0
  
 
      MS  =  499.1246D0
      MV  =  782.5D0
      MR  =  770.D0
      MD  =  980.D0

   
  
      B2 = GS*MS*MS*K2/(2.D0*MB)
      B3 = GS*GS*MS*MS*K3/(6.D0*MB*MB)
      C1 = GV*GV*EP3/(6.D0)
      D2 = GS*MV*MV*ET2/(2.D0*MB)
      D3 = GS*GS*MV*MV*ET3/(2.D0*MB*MB)  
      F2 = ETR*GS*MR*MR/(2.D0*MB)
      G3 = GS*GS*MR*MR*ETR3/(2.D0*MB*MB)
      G4 = GV*GV*MR*MR*ETR4/(2.D0*MB*MB)
      G5 = GR*GR*EPR3/(6.D0)     


c--------------------------------------------------------------
      RETURN
      END 
c---------------------------------------------------------------
C
C Parameter BSP B
C  KFO =1.303
C
      SUBROUTINE SET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
C     Initialisation
      PI  = 3.14159265358979D0
      HC  = 197.327D0
      ME  = 0.511D0
      MO  = 105.6D0
      MB  = 939.20D0
    
      K2  = -1.0681D0
      K3  = 14.9857D0
      ET2 = -0.0872D0
      ET3 = 3.1265D0

      EP3 = 0.0D0
      EPR3=0.0D0

      MS = 0.5384D0*MB
      MV  =0.8333D0*MB
      MR  =0.8200D0*MB
      MD  = 980D0

      ETR = 0.0D0   
      ETR3 = 0.0D0
      ETR4 =  53.7642D0
 

      GS = 0.8764D0*4.D0*PI
      GV = 1.1481D0*4.D0*PI
      GR = 1.0508D0*4.D0*PI
 
      GD = 0.0D0
     
     

      B2 = GS*MS*MS*K2/(2.D0*MB)
      B3 = GS*GS*MS*MS*K3/(6.D0*MB*MB)
      C1 = GV*GV*EP3/(6.D0)
      D2 = GS*MV*MV*ET2/(2.D0*MB)
      D3 = GS*GS*MV*MV*ET3/(2.D0*MB*MB)  
      F2 = ETR*GS*MR*MR/(2.D0*MB)
      G3 = GS*GS*MR*MR*ETR3/(2.D0*MB*MB)
      G4 = GV*GV*MR*MR*ETR4/(2.D0*MB*MB)
      G5 = GR*GR*EPR3/(6.D0)     


c--------------------------------------------------------------
      RETURN
      END 
c---------------------------------------------------------------
