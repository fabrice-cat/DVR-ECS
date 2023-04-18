module coulombwave

IMPLICIT NONE;

CONTAINS



! *** CALCUL EN DOUBLE PRECISION DES FONCTIONS COULOMBIENNES REGULIERES 
!     ET DES FONCTIONS DE RICATTI-BESSEL D'ARGUMENTS REELS.             
!                                                                       
!     -1000. < ETA < -.001  OU  .001 < ETA <1000.  OU  ETA=0.           
!     0. < RHO < 20000.                                                 
!                                                                       
! *** D'APRES A.R.BARNETT ET AL.                                        
!             COULOMB WAVE FUNCTIONS FOR ALL REAL ETA AND RHO.          
!             COMPUTER PHYSICS COMMUNICATIONS 8 (1974) 377-395.         
!                                                                       
! *** POUR LES CONVENTIONS MATHEMATIQUES ADOPTEES, VOIR AUSSI           
!     ABRAMOVITCH AND STEGUN, P538.                                     
!                                                                       
!                                                                       
! *** LISTE DES ARGUMENTS DE RCWFN.                                     
!     *****************************                                     
!                                                                       
! *** RHO,ETA                     REAL*8                                
! *** MINL0,MAXL                  ENTIER. ON POSE MINL=IABS(MINL0).     
! *** FC,FCP,GC,GCP               TABLEAUX REAL*8                       
!                                 DIM MAXL+1 AU MOINS SI MINL0>=0.      
!                                 DIM MAXL+1-MINL AU MOINS SINON.       
!     LA FONCTION REGULIERE FL, ET SA DERIVEE FPL,                      
!     LA FONCTION IRREGULIERE GL ET SA DERIVEE GPL,                     
!     SONT CALCULEES POUR L=MINL,...MAXL.                               
!     LE RESULTAT EST RENVOYE DANS LES TABLEAUX FC,FCP,GC,GCP.          
!       I) MINL0 >= 0: STOCKAGE NORMAL.                                 
!             FC (L+1) = FL                                             
!             FCP(L+1) = FPL                                            
!             GC (L+1) = GL                                             
!             GCP(L+1) = GPL,  L=MINL,...MAXL.                          
!      II) MINL0 < 0: STOCKAGE CONDENSE.                                
!             FC (L+1-MINL) = FL                                        
!             FCP(L+1-MINL) = FPL                                       
!             GC (L+1-MINL) = GL                                        
!             GCP(L+1-MINL) = GPL,  L=MINL,...MAXL.                     
! *** ACCUR                       REAL*8                                
!     PRECISION DEMANDEE.                                               
! *** KOD                         ENTIER                                
!     CODE D'ERREUR.                                                    
!     SI UNE ERREUR SE PRODUIT, FC,FCP,GC,GCP SONT RETOURNES NULS,      
!     ET KOD INDIQUE LA NATURE DE L'ERREUR.                             
!                                                                       
! *** KOD=0                                                             
!     CALCUL SANS INCIDENTS.                                            
! *** KOD=1                                                             
!     ETA TROP PETIT EN VALEUR ABSOLUE.                                 
!     NON CONVERGENCE DE LA FRACTION CONTINUE DONNANT (GP+I*FP)/(G+I*F).
! *** KOD=2                                                             
!     RHO TROP GRAND.                                                   
!     NON CONVERGENCE DE LA FRACTION CONTINUE DONNANT FP/F.             
! *** KOD=3                                                             
!     REGION CLASSIQUEMENT INTERDITE EN L=MAXL.                         
!     (MAXL-MINL) TROP GRAND. RISQUE D'OVERFLOWS.                       
! *** KOD=4                                                             
!     REGION CLASSIQUEMENT INTERDITE EN L=MINL.                         
!     CALCUL POSSIBLE POUR L=MINL.                                      
!     CALCUL IMPOSSIBLE POUR L=MAXL.                                    
!     (MAXL-MINL) TROP GRAND. RISQUE D'OVERFLOWS.                       
! *** KOD=5                                                             
!     REGION CLASSIQUEMENT INTERDITE EN L=MINL.                         
!     CALCUL IMPOSSIBLE POUR L=MINL.                                    
!     RHO TROP PETIT OU MINL TROP GRAND.                                
!     RISQUE D'OVERFLOWS DANS L'INTEGRATION DE G,GP A PARTIR            
!     DU POINT RHO=TURN.                                                
! *** KOD=6                                                             
!     REGION CLASSIQUEMENT INTERDITE EN L=MINL.                         
!     CALCUL IMPOSSIBLE POUR L=MINL.                                    
!     TROP DE PAS D'INTEGRATION.                                        
! *** KOD=7                                                             
!     REGION CLASSIQUEMENT INTERDITE EN L=MINL.                         
!     CALCUL IMPOSSIBLE POUR L=MINL.                                    
!     SINGULARITE DE G,GP LORS DE L'INTEGRATION.                        
!     ARGUMENTS INCORRECTS EN ENTREE DE RCWFN.                          
!                                                                       
!                                                                       
!                                                                       
      SUBROUTINE RCWFN(RHO,ETA,MINL0,MAXL,FC,FCP,GC,GCP,ACCUR,KOD)      
      
      real*8,intent(in) :: RHO,ETA,ACCUR;
      integer,intent(in) :: MINL0,MAXL;
 
      integer :: KOD,MINL,IBORNE,KTR,KTRP,LMAX,IPIC,L,NSTEP,LMIN1,LP;

      real*8 :: K,OVRFLW,UNDFLW,ACC,R,ETA2,TURN,XLL1,PL,PMX,FP,DK,DEL,F,P,Q,AR,AI,BR,BI,WI,DR,DI,DP,DQ;
      real*8 :: H,T,G,GP,FPF,GF,GPF,W,TF,TFP,ETAR,RHO2,D,X,Y,DYDX,DELTA,DELABS;
                                   

      real*8,DIMENSION(MAXL+1) :: FC,FCP,GC,GCP;

                                                    
!      DIMENSION FC(MAXL+1),FCP(MAXL+1),GC(MAXL+1),GCP(MAXL+1)
!**      DATA OVRFLW/1.0D+60/,UNDFLW/1.D-36/                               
      DATA OVRFLW/1.0D+300/;
      DATA UNDFLW/1.D-300/;                               
!**      FLOAT(I)=DFLOAT(I)                                                
!**      EXP(X)=DEXP(X)                                                    
!**      ALOG(X)=DLOG(X)                                                   
!**      ABS(X)=DABS(X)                                                    
!**      IFIX(X)=IDINT(X)                                                  
!**      SQRT(X)=DSQRT(X)                                                  
! *** COULOMB WAVEFUNCTIONS CALCULATED AT R = RHO BY THE                
! *** CONTINUED-FRACTION METHOD OF STEED   MINL,MAXL ARE ACTUAL L-VALUES
! *** SEE BARNETT FENG STEED AND GOLDFARB COMPUTER PHYSICS COMMUN 1974  
      KOD=0                                                             
      MINL=IABS(MINL0)                                                  
      IBORNE=0                                                          
      IF(MINL0.LT.0) IBORNE=MINL                                        
      ACC  = ACCUR                                                      
      IF(ACC.GT.(1.0d-6)) ACC=1.0d-6                                      
      IF(ACC.LT.(1.0d-14)) ACC=1.0d-14                                    
      R    = RHO                                                        
      KTR  = 1                                                          
      LMAX = MAXL                                                       
      LMIN1= MINL + 1                                                   
      XLL1 = DFLOAT(MINL*LMIN1)                                          
      ETA2 = ETA*ETA                                                    
      TURN = ETA + DSQRT(ETA2 + XLL1)                                    
      IF(R.LT.TURN) KTR=-1                                              
      KTRP = KTR                                                        
      GO TO 2                                                           
1     R    = TURN                                                       
      TF   = F                                                          
      TFP  = FP                                                         
      LMAX = MINL                                                       
      KTRP = 1                                                          
2     ETAR = ETA*R                                                      
      RHO2 =   R*R                                                      
      PL   = DFLOAT(LMAX + 1)                                            
      PMX  = PL + 0.5d0                                                   
! *** CONTINUED FRACTION FOR FP(MAXL)/F(MAXL)  XL IS F  XLPRIME IS FP **
      FP  = ETA/PL + PL/R                                               
      DK  = ETAR*2.0d0                                                    
      DEL = 0.0d0;                                                         
      D   = 0.0d0;                                                         
      F   = 1.0d0;                                                         
      K   = (PL*PL - PL + ETAR)*(2.0d0*PL - 1.0d0)                          
      IF(DABS(PL*PL+PL+ETAR).GT.UNDFLW) GO TO 3                         
1000  R   = R + ACC                                                     
      GO TO 2                                                           
3     H   = (PL*PL + ETA2)*(1.0d0 - PL*PL)*RHO2                           
      K   = K + DK + PL*PL*6.0d0                                          
      D   =  1.0d0/(D*H + K)                                              
      DEL =  DEL*(D*K - 1.0d0)                                            
      IF(PL.LT.PMX) DEL = -R*(PL*PL + ETA2)*(PL + 1.0d0)*D/PL             
      PL  = PL + 1.0d0                                                    
      FP  = FP + DEL                                                    
      IF(D.LT.0.0) F = -F                                               
      IF(PL.GT.20000.) GOTO 102                                         
      IF(DABS(FP).LT.UNDFLW) GO TO 1000                                 
      IF(ABS(DEL/FP).GE.ACC) GO TO 3                                    
      FP  = F*FP                                                        
      IF( LMAX.EQ.MINL) GO TO 5                                         
      IPIC=LMAX+1-IBORNE                                                
      FC (IPIC) = F                                                     
      FCP(IPIC) = FP                                                    
! *** DOWNWARD RECURSION TO MINL FOR F AND FP, ARRAYS GC,GCP ARE STORAGE
      L  = LMAX                                                         
      DO 4 LP  = LMIN1,LMAX                                             
      PL = DFLOAT(L)                                                     
      IPIC=L-IBORNE                                                     
      GC (IPIC+1) = ETA/PL + PL/R                                       
      GCP(IPIC+1) = DSQRT(ETA2 + PL*PL)/PL                               
      FC (IPIC)   = (GC(IPIC+1)*FC(IPIC+1) + FCP(IPIC+1))/GCP(IPIC+1)   
      FCP(IPIC)   =  GC(IPIC+1)*FC(IPIC)   - GCP(IPIC+1)*FC(IPIC+1)     
      IF(DABS(FCP(IPIC)).GT.OVRFLW) GOTO 103                             
4     L  = L - 1                                                        
      IPIC=LMIN1-IBORNE                                                 
      F  = FC (IPIC)                                                    
      FP = FCP(IPIC)                                                    
5     IF(KTRP.EQ.-1) GO TO 1                                            
! *** REPEAT FOR R = TURN IF RHO LT TURN                                
! *** NOW OBTAIN P + I.Q FOR MINL FROM CONTINUED FRACTION (32)          
! *** REAL ARITHMETIC TO FACILITATE CONVERSION TO IBM USING REAL*8      
      P  = 0.0d0                                                          
      Q  = R - ETA                                                      
      PL = 0.0d0                                                          
      AR = -(ETA2 + XLL1)                                               
      AI =   ETA                                                        
      BR = 2.0d0*Q                                                        
      BI = 2.0d0                                                          
      WI = 2.0d0*ETA                                                      
      DR =   BR/(BR*BR + BI*BI)                                         
      DI =  -BI/(BR*BR + BI*BI)                                         
      DP = -(AR*DI + AI*DR)                                             
      DQ =  (AR*DR - AI*DI)                                             
6     P  =  P + DP                                                      
      Q  =  Q + DQ                                                      
      PL = PL + 2.0d0                                                     
      AR = AR + PL                                                      
      AI = AI + WI                                                      
      BI = BI + 2.0d0                                                     
      D  = AR*DR - AI*DI + BR                                           
      DI = AI*DR + AR*DI + BI                                           
      T  = 1.0d0/(D*D + DI*DI)                                            
      DR =  T*D                                                         
      DI = -T*DI                                                        
      H  = BR*DR - BI*DI - 1.0d0                                          
      K  = BI*DR + BR*DI                                                
      T  = DP*H  - DQ*K                                                 
      DQ = DP*K  + DQ*H                                                 
      DP = T                                                            
      IF(PL.GT.46000.) GOTO 101                                         
      IF(DABS(DP)+DABS(DQ).GE.(DABS(P)+DABS(Q))*ACC) GO TO 6                
      P  = P/R                                                          
      Q  = Q/R                                                          
! *** SOLVE FOR FP,G,GP AND NORMALISE F  AT L=MINL                      
      G  = (FP - P*F)/Q                                                 
      GP = P*G - Q*F                                                    
      FPF=FP/F                                                          
      GF =G/F                                                           
      GPF=GP/F                                                          
      W=1.0/(DABS(F)*DSQRT(FPF*GF-GPF))                                   
      G  = W*G                                                          
      GP = W*GP                                                         
      IF(KTR.EQ.1) GO TO 8                                              
      F  = TF                                                           
      FP = TFP                                                          
      LMAX = MAXL                                                       
! *** INTEGRATION OF G(MINL) AND GP(MINL) A R DECROISSANT DEPUIS TURN   
! *** PAR UNE METHODE D'EXTRAPOLATION A PAS VARIABLE.                   
! *** CF. DIPLOMARBEIT VON HANS-GUNTER HUSSELS,                         
! ***     AIMABLEMENT COMMUNIQUE PAR LE PROF. DR. R BURLIRSCH.          
      X=TURN                                                            
      Y=G                                                               
      DYDX=GP                                                           
      H=.1                                                              
      NSTEP=0                                                           
   21 DELTA=RHO-X                                                       
      DELABS=DABS(DELTA)                                                
      IF(DELABS/RHO.LT.1.D-14) GOTO 22                                  
      H=DSIGN(DMIN1(DABS(H),DELABS),DELTA)                              
      CALL DIFSYL(ETA+ETA,XLL1,ACC,H,X,Y,DYDX)                          
      IF(DABS(DYDX).GT.OVRFLW) GOTO 105                                  
      IF(H.EQ.0.0) GOTO 107                                             
      NSTEP=NSTEP+1                                                     
      IF(NSTEP.GT.20000) GOTO 106                                       
      GOTO 21                                                           
   22 G=Y                                                               
      GP=DYDX                                                           
      IF(DSQRT(DABS(FP))*DSQRT(DABS(GP)).GT.DSQRT(OVRFLW)) GOTO 104          
      W  = 1.0/(FP*G - F*GP)                                            
! *** UPWARD RECURSION FROM GC(MINL) AND GCP(MINL),STORED VALUES ARE R,S
! *** RENORMALISE FC,FCP FOR EACH L-VALUE                               
8     IPIC=LMIN1-IBORNE                                                 
      GC (IPIC) = G                                                     
      GCP(IPIC) = GP                                                    
      IF(LMAX.EQ.MINL) GO TO 10                                         
      DO  9  L = LMIN1,LMAX                                             
      IPIC=L-IBORNE                                                     
      T        = GC(IPIC+1)                                             
      GC (IPIC+1) = (GC(IPIC)*GC (IPIC+1) - GCP(IPIC))/GCP(IPIC+1)      
      GCP(IPIC+1) =  GC(IPIC)*GCP(IPIC+1) - GC(IPIC+1)*T                
      IF(ABS(GCP(IPIC+1)).GT.OVRFLW) GOTO 103                           
      FC (IPIC+1) = W*FC (IPIC+1)                                       
9     FCP(IPIC+1) = W*FCP(IPIC+1)                                       
      IPIC=LMIN1-IBORNE                                                 
      FC (IPIC) = FC (IPIC)*W                                           
      FCP(IPIC) = FCP(IPIC)*W                                           
      RETURN                                                            
10    IPIC=LMIN1-IBORNE                                                 
      FC (IPIC) = W*F                                                   
      FCP(IPIC) = W*FP                                                  
      RETURN                                                            
11    W  = 0.0                                                          
      G  = 0.0                                                          
      GP = 0.0                                                          
      GO TO 8                                                           
  101 KOD=1                                                             
      GOTO 11                                                           
  102 KOD=2                                                             
      GOTO 11                                                           
  103 KOD=3                                                             
      GOTO 11                                                           
  104 KOD=4                                                             
      GOTO 11                                                           
  105 KOD=5                                                             
      GOTO 11                                                           
  106 KOD=6                                                             
      GOTO 11                                                           
  107 KOD=7                                                             
      GOTO 11                                                           
      END SUBROUTINE RCWFN  


                                                               
      SUBROUTINE DIFSYL(ETA2,XLL1,EPS,H,X,Y,DYDX)                       

      real*8 :: ETA2,XLL1,EPS,H,X,Y,DYDX;
 
      integer :: J,JTI,L,M,NI,I,N,JR,JS,K;

      real*8 :: COU,B,G,C,U,V,DZ2,XN,DY2,TA,B1,W;
                                   
      REAL*4 :: FA,FV,ETA,FY,FS,FLOAT,SNGL;

      real*8,DIMENSION(2) :: YR,YS,S; 
      real*8,DIMENSION(4) :: EP; 
      real*8,DIMENSION(7) :: D; 
      real*8,DIMENSION(2,7) :: DT;
               
      LOGICAL KONV,BO,KL,GR;
                                           
      DATA EP/0.4E-1,0.16E-2,0.64E-4,0.256E-5/                          
      DATA N/1/
                                                         
      COU(X)=ETA2/X+XLL1/(X*X)-1.D0                                     
      JTI=0                                                             
      FY=1.0d0;                                                             
      ETA=SNGL(DABS(EPS))                                               
      IF(ETA.LT.(1.0D-11)) ETA=1.0D-11                                      
      DZ2=COU(X)*Y                                                      
   10 XN=X+H                                                            
      BO=.FALSE.                                                        
      S(1)=DABS(Y)                                                      
      S(2)=DABS(DYDX)                                                   
      M=1                                                               
      JR=2                                                              
      JS=3                                                              
      DO 260 J=1,10                                                     
      IF(.NOT.BO) GOTO 200                                              
      D(2)=1.777777777777778D0                                          
      D(4)=7.111111111111111D0                                          
      D(6)=2.844444444444444D1                                          
      GOTO 201                                                          
  200 D(2)=2.25D0                                                       
      D(4)=9.D0                                                         
      D(6)=3.6D1                                                        
  201 IF(J.LE.7) GOTO 202                                               
      L=7                                                               
      D(7)=6.4D1                                                        
      GOTO 203                                                          
  202 L=J                                                               
      D(L)=M*M                                                          
  203 KONV=L.GT.3                                                       
      B=H/DFLOAT(M)                                                     
      G=B*0.5D00                                                        
      YS(2)=DYDX+G*DZ2                                                  
      YS(1)=Y+B*YS(2)                                                   
      M=M-1                                                             
      IF(IABS(M).LT.1) GOTO 221                                         
      DO 220 K=1,M                                                      
      DY2=COU(X+DFLOAT(K)*B)*YS(1)                                      
      I=1                                                               
      U=DABS(YS(I))                                                     
      IF(U.GT.S(I)) S(I)=U                                              
      NI=2                                                              
      V=YS(NI)+B*DY2                                                    
      YS(I)=YS(I)+B*V                                                   
      YS(NI)=V                                                          
      V=DABS(V)                                                         
      IF(V.GT.S(NI)) S(NI)=V                                            
  220 CONTINUE                                                          
  221 DY2=COU(XN)*YS(1)                                                 
      YS(2)=YS(2)+G*DY2                                                 
      KL=L.LT.2                                                         
      GR=L.GT.5                                                         
      FS=0.                                                             
      M=N+N                                                             
      DO 233 I=1,M                                                      
      V=DT(I,1)                                                         
      C=YS(I)                                                           
      U=DABS(C)                                                         
      IF(U.GT.S(I)) S(I)=U                                              
      DT(I,1)=C                                                         
      TA=C                                                              
      IF(KL) GOTO 233                                                   
      DO 231 K=2,L                                                      
      B1=D(K)*V                                                         
      B=B1-C                                                            
      W=C-V                                                             
      U=V                                                               
      IF(B.EQ.0.D0) GOTO 230                                            
      B=W/B                                                             
      U=C*B                                                             
      C=B1*B                                                            
  230 V=DT(I,K)                                                         
      DT(I,K)=U                                                         
  231 TA=U+TA                                                           
      IF(.NOT.KONV) GOTO 232                                            
      IF(DABS(YR(I)-TA).GT.S(I)*ETA) KONV=.FALSE.                       
  232 IF(GR.OR.S(I).EQ.0.D0) GOTO 233                                   
      FV=DABS(W)/S(I)                                                   
      IF(FS.LT.FV) FS=FV                                                
  233 YR(I)=TA                                                          
      IF(FS.EQ.0.D0) GOTO 250                                           
      FA=FY                                                             
      K=L-1                                                             
      FY=(EP(K)/FS)**(1./FLOAT(L+K))                                    
      IF(L.EQ.2) GOTO 240                                               
      IF(FY.LT.0.7*FA) GOTO 250                                         
  240 IF(FY.GT.0.7) GOTO 250                                            
      H=H*FY                                                            
      JTI=JTI+1                                                         
      IF(JTI.GT.5) GOTO 30                                              
      GOTO 10                                                           
  250 IF(KONV) GOTO 20                                                  
      D(3)=4.D0                                                         
      D(5)=1.6D1                                                        
      BO=.NOT.BO                                                        
      M=JR                                                              
      JR=JS                                                             
  260 JS=M+M                                                            
      H=H*0.5D0                                                         
      GOTO 10                                                           
   20 X=XN                                                              
      H=H*FY                                                            
      Y=YR(1)                                                           
      DYDX=YR(2)                                                        
      RETURN                                                            
   30 H=0.D00                                                           
      RETURN                                                            
      END SUBROUTINE DIFSYL             
      





      complex*16 function cgamma(z)

!  *********************************************************************
!  *                                                                   *
!  *    evaluation of the gamma function for a complex argument.       *
!  *                                                                   *
!  *    algorithm 404 collected algorithms from acm.                   *
!  *    algorithm appeared in comm. acm, vol. 14, no. 01, p. 048.      *
!  *                                                                   *
!  ********************************************************************* 

      implicit none

!  ****** VARIABLE DECLARATION ******

      integer        m, j, i

      real*8         ppi, x, y, tol, xdist

      real*8         c(12)

      complex*16     z, zm, t, tt, sum, term, den, pi, a

      logical        reflek

!  ****** EVALUATES GAMMA FUNCTION ******

!  ---- initialization

      ppi=4.0d0*atan(1.0d0)
      pi = ppi
      x = dreal(z)
      y = dimag(z)

!  ---- tol = limit of precision of computer system in single precisi

      tol = 1.0d-15
      reflek = .true.

!  ---- determine whether z is too close to a pole
!  ---- check whether too close to origin

      if(x.ge.tol) go to 20

!  ---- find the nearest pole and compute distance to it

      xdist = x-int(x-.5d0)
      zm = dcmplx(xdist,y)
      if(abs(zm).ge.tol) go to 10

!  ---- if z is too close to a pole, cgamma is set equal to 1.e15

      cgamma = (1.d15,0.d0)
      return

!  ---- for real(z) negative employ the reflection formula
!  ---- gamma(z) = pi/(sin(pi*z)*gamma(1-z))
!  ---- and compute gamma(1-z).  note reflek is a tag to indicate tha
!  ---- this relation must be used later.

10    if(x.ge.0.0d0) go to 20
      reflek = .false.
      z = dcmplx(1.0d0,0.0d0)-z
      x = 1.0d0-x
      y = -y

!  ---- if z is not too close to a pole, make real(z)>10 and arg(z)<p

20    m = 0
40    if(x.ge.10.d0) go to 50
      x = x + 1.0d0
      m = m + 1
      go to 40
50    if(abs(y).lt.x) go to 60
      x = x + 1.0d0
      m = m + 1
      go to 50
60    t = dcmplx(x,y)
      tt = t*t
      den = t

!  ---- coefficients in stirling*s approximation for ln(gamma(t))

      c(1) = 1.0d0/12.0d0
      c(2) = -1.0d0/360.0d0
      c(3) = 1.0d0/1260.0d0
      c(4) = -1.0d0/1680.0d0
      c(5) = 1.0d0/1188.0d0
      c(6) = -691.0d0/360360.0d0
      c(7) = 1.0d0/156.0d0
      c(8) = -3617.0d0/122400.0d0
      c(9) = 43867.0d0/244188.0d0
      c(10) = -174611.0d0/125400.0d0
      c(11) = 77683.0d0/5796.0d0

      sum = (t-0.5d0)*log(t)-t+dcmplx(0.5d0*log(2.0d0*ppi),0.0d0)
      j = 1
70    term = c(j)/den

!  ---- test real and imaginary parts of ln(gamma(z)) separately for
!  ---- convergence.  if z is real skip imaginary part of check.

      if(abs(dreal(term)/dreal(sum)).ge.tol) go to 80
      if(y.eq.0.0d0) go to 100
      if(abs(dimag(term)/dimag(sum)).lt.tol) go to 100
80    sum = sum + term
      j = j + 1
      den = den*tt

!  ---- test for nonconvergence

      if(j.eq.12) go to 90
      go to 70

!  ---- stirling*s series did not converge.  print error message and
!  ---- procede.

90    print*,'Stirling series did not converged --> code stopped'
      stop

!  ---- recursion relation used to obtain ln(gamma(z))
!  ---- ln(gamma(z)) = ln(gamma(z+m)/(z*(z+1)*...*(z+m-1)))
!  ---- = ln(gamma(z+m)-ln(z)-ln(z+1)-...-ln(z+m

100   if(m.eq.0) go to 120
      do 110 i = 1,m
      a = dcmplx(i*1.d0-1.d0,0.0d0)
110   sum = sum-log(z+a)

!  ---- check to see if reflection formula should be used

120   if(reflek) go to 130
      sum = log(pi/sin(pi*z))-sum
      z = (1.0d0,0.0d0) -z
130   cgamma = exp(sum)
      return
      end function cgamma      
      
END MODULE coulombwave
