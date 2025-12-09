      PROGRAM TOVSOLVER                   
C     *********************************************************
C     
C     TEST1 linear EOS+corr initial boundary 28 Mei 2020
C                  
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IL, IN    
      DIMENSION YA(10), EK(4,10), Y(10)


c      OPEN (unit=2,STATUS='unknown',FILE='CatatanB145SC.dat')
       OPEN (unit=3,STATUS='unknown',FILE='radmass-TOViso.dat')
       OPEN (unit=1,STATUS='unknown',FILE='profil-TOViso.dat')

C     IM = NUMBER OF EQUATIONS Y(I)=Pressure, Y(2)=NS Mass and Y(3)=E density

      HC  = 197.327D0
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
      IM=3
      IN=IM-1
c---------------------------------------------------------------------
c     XL=1.0D-4
C      XL=1.0D-2
C      ALPHA=5.0D-3
C      XC=XL*ALPHA      
      XC=1.0D-3

c---------------------------------------------------------------------      
       DO 10 IL=2,600,2
       FIXEDIL=300
       !IL=FIXEDIL
       PC=1.D0*IL
 
       YA(10)=XL
c       PC=800.D0
       EDC=FED(PC)
       FIXP=1.D0-2.D0*XL*XL*PI*GS*(3.D0*PC+EDC)/3.D0
c     FIXM=1.D0+XL*XL/(XC*XC)-XL*XL*XL/(XC*XC*XC)*DATANH(XC/XL)
c     we used series expansion above because in fortran77 no function arctanh.
c     In fortran 90/95 the function is exist, but because the code is written
c     with 77 version, if we compile with gfortran, the code is not always
c     stable                   !!!       
       !FIXM=2.D0/3.D0-(XC/XL)**2/5.D0-(XC/XL)**4/7.D0-(XC/XL)**6/9.D0
       FIXM=1.D0
       IF (ABS(ALPHA)<1)  FIXM=1.D0 + (XC/XL)**(-2)  
     &     - LOG((1 + (XC/XL))/(1 - (XC/XL)))/(2.D0*(XC/XL)**3)

c       FIXP=1.0D0
c       FIXM=1.0D0
       
       PCC=PC !-2.D0*PI*GS*XC*XC*(PC+EDC)*(3.D0*PC+EDC)*FIXP/3.D0
       MCC=4.D0*PI*XC*XC*XC*EDC/3.D0 !*FIXM

       
       
      Y(1)=PCC

      Y(2)=MCC
c-----------------------------------------------------------------------
c     Y(1)=PC 
c      Y(2)=1.0D-8
      
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      
      P0=Y(1)
 
      Y(3)=FED(P0)
    

     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D0
      NS=8
      !XL=30.0D3

      H=PU/NS
c     XP should be larger than XC in order to avoid the unphysical behavior
c      near center!!    
      XP=1.0D0
      HH=H/(2.0D0)

      IF (XP.LT.XC) THEN
            WRITE(*,*) "XP=",XP," is NOT larger than XC=",XC
            STOP
      END IF

C     LINE NUMBER INITIALIZATION

      LI=0
 

 28   LI=LI+1

C     XB=OLD radius, XP=NEW radius, XM=MIDPOINT radius

      DO N=1,NS
         XB=XP
         XP=XP+H
         XM=XB+HH

C    COMPUTE K1, L1

         J=1
         DO I=1,IM
            YA(I)=Y(I)
         END DO
         XA=XB
   
         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K2, L2

         J=2
         
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
            P0=YA(1)
            
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      


         YA(3)=FED(P0)

         XA=XM

         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K3, L3

         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
     

          P0=YA(1)
          
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

  
            YA(3)=FED(P0)          

         XA=XM


         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K4, L4

         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO

          P0=YA(1)
         
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

 

            YA(3)=FED(P0)
     

         XA=XP

         CALL FUNCT(EK,J,YA,XA,H)

C    4-TH ORDER RUNGE KUTTA SCHEME

         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
          P0=Y(1)
         
          
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

  


          Y(3)=FED(P0)
          
       END DO
       
       XA=XP
       PRESS=Y(1)
       MASST=Y(2)
       EDEN=Y(3)
       XX=(2.D0*GS*MASST*MSS/XA)
     &   *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MASST*MSS))
     &   /(1.D0-2.D0*GS*MASST*MSS/XA)
       TPA=SQRT(1.D0+(XL/XA)**2*XX)

      !WRITE(*,*)IL,LI,(XP/1.D3),Y(1),Y(2),Y(3)
      IF (IL.EQ.FIXEDIL) WRITE(1,*)IL,LI,(XP/1.D3),Y(1),Y(2),TPA
     

       PS=Y(1)
       PMIN=1.0D-8
c       PMIN=2.0D-5
      
      IF (MCC.LT.0.0D0) THEN
            WRITE(*,*) "MCC negative"
            WRITE(*,*) Y(2)
            GOTO 10
      ELSE IF (ABS(PS).GT.10*ABS(PCC)) THEN
            WRITE(*,*) "ABS(P) =",ABS(PS),"> 10 ABS(PCC) =",10*ABS(PCC)
            GOTO 10
      ELSE IF (LI.GE.1000000) THEN
            WRITE(*,*) "Max iteration = 1000000 reached"
            WRITE(*,*)"p=",Y(1),"m=",Y(2),"rho=",Y(3)
            GOTO 10
      END IF

      IF (PS .GT. PMIN  ) GOTO 28
    


c      WRITE(2,*)IL,(XP/1.D3), (Y(I),I=1,IM)
        WRITE(3,*)IL,(EDC/1.D3),(XP/1.D3),MASST,
     &             2.0D0*GS*Y(2)*MSS/XP     
c       write(*,*)X
        WRITE(*,*)IL,(EDC/1.D3),(XP/1.D3),MASST,
     &             2.0D0*GS*Y(2)*MSS/XP
 
   

  10   CONTINUE
      
      
      END

      SUBROUTINE FUNCT(EK,J,YA,XA,H)
C     *********************************************************
C     DEFINES TOV EQUATIONS 
C
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(4,10), YA(10)
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
c--------------------------------------------------------------
c     XL=1.0D-4
      XL=YA(10)
      EDEN=YA(3)
      PRESS=YA(1)
      MASST=YA(2)

       XX=(2.D0*GS*MASST*MSS/XA)
     &   *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MASST*MSS))
     &   /(1.D0-2.D0*GS*MASST*MSS/XA)
c-------------------------------------------------------------      
c       TPA=0.0D0
c       TMA=1.0D0
c      OPEN (unit=7,STATUS='unknown',FILE='cekTMAXL10mm.dat')
c      write(7,*)TPA,FF,GG,TMA
c-------------------------------------------------------------      
  
      EK(J,1)=-EDEN*(1.D0+PRESS/EDEN)*XX/(2.D0*XA)*H
      
      EK(J,2)=4.D0*PI*XA*XA*EDEN*H/MSS 

      RETURN
      END

 
c---------------------------------------------------------------------       
c     GUP EOS BETA=0.D0
c----------------------------------------------------------------
  
      
      FUNCTION SIG(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      
      SIG = 0.D0
   
      RETURN
      END
      
      FUNCTION FED(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      
       IF ( P0 .GT. 50.D0 ) THEN
	 
     	FED= 252.7351203007992D0 + 2.0517873998109133D0*P0
   
      ELSE IF ( P0 .GT.  0.5658D0 .AND. P0 .LE. 50.D0 ) THEN
 
       FED= 40.72822519163616D0 + 172.88564172842007D0*P0 - 
     -  129.75665741187734D0*P0**2 + 
     -  57.21198616046201D0*P0**3 - 15.6150708008523D0*P0**4 + 
     -  2.8354313027185887D0*P0**5 - 0.3602847077751246D0*P0**6 + 
     -  0.03316049472454524D0*P0**7 - 0.0022626453096747964D0*P0**8 + 
     -  0.00011613173543073938D0*P0**9 - 4.5162772145274745D-6*P0**10 + 
     -  1.3312210942573273D-7*P0**11 - 2.952609024402367D-9*P0**12 + 
     -  4.845080539678656D-11*P0**13 - 5.701316398676422D-13*P0**14 + 
     -  4.5477089999551786D-15*P0**15 - 2.2015407469797255D-17*P0**16 + 
     -  4.881799449675595D-20*P0**17
        
       ELSE IF (P0 .GT. 4.99313436D-4 .AND. P0 .LE. 0.5658D0) THEN
   

       FED= 0.2138614869659046D0 + 768.9378464621362D0*P0 - 
     -  6179.085226412627D0*P0**2 + 
     -  31192.024180740616D0*P0**3 - 85373.61029759394D0*P0**4 + 
     -  117008.53995402799D0*P0**5 - 62844.457041598056D0*P0**6

       ELSE
     
       FED=  0.00020663104786823767D0 + 985.7550962048155D0*P0 - 
     -  6.452649548410661D6*P0**2 + 4.045493650683346D10*P0**3 - 
     -  1.2422017897554195D14*P0**4 + 1.765493095354902D17*P0**5 - 
     -  9.1529417012139D19*P0**6
 
      
       END IF
   
      RETURN
      END
c-----------------------------------------------------------------------      
