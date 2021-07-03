c-------------------------------------------------------------------- 
c                    Set Initial Conditions for N-body simulations
c                    Uses Jeans equation for spherically symmetric
c                    configuration with isotropic velocities
c  
c-------------------------------------------------------------------- 
      INCLUDE '../nbody1.h'
      REAL*8    Ekin,Epot,Etot,yc(3)

                          ! open files
      Open(1,file='../DATA/ic.dat',status='unknown')

      n=np
                          ! set parameters
      do i=1,3         ! offset for the center of snapshot
         yc(i)   = 0.
      EndDo 


      amasstotal=amasstotal0

      Call Initial(aMassTotal)

 20   format(f8.3,i8,6g12.5)
 30   format(3x,'Time',T13,'Step',T19,'Ekin',T30,'Epot',
     &          T42,'Etot',T53,'Ekin/Etot')
 40   format(3x,'Time',T13,'Step',T19,'TotalM',
     &          T28,'MCentral',
     &          T37,'Rmax',T44,'Raverage',
     &          15(' Mass   '))
      Stop
      End
c-------------------------------------------------------------------- 
c              make initial conditions
      SUBROUTINE Initial(aMassTotal)
c              Use Jeans equation to set equilibrium spherical
c              configuration with isotropic velocities.
c                     Coords: 1-3 = x,y,z
c                             4-6 = vx,vy,vz
c                             10  = mass
c-------------------------------------------------------------------- 
      INCLUDE '../nbody1.h'
      Real*8  sMass,sx,sy,sz,s1,s2,s3,dMass,dSigma,dfSigma

      Nseed = 121071            ! Initial seed for random numbers
      sMass =0.
      sx    =0.
      sy    =0.
      sz    =0.
      s1    =0.
      s2    =0.
      s3    =0.
      ntry  =0                   ! counter of the total number of trials

      x = Rmininit               ! find Max of weighting function
      rMax = 0.                  ! Need it for rejection method
      rMin = 0.
      do while(x.lt.Rmaxinit)
         xx = weight(x)
         rMax =max(rMax,xx)
          x =x*1.05
c          x1 = dMass(dble(x))
c          x2 = dSigma(x)
c          x3 = dfSigma(dble(x))
      EndDo 


         Do i=1,Np               ! make coordinates of particles
 10         radius = (RANDD(Nseed))**2*Rmaxinit
            reject = (weight(radius)-rMin)/rMax
            reject = min(max(reject,0.),1.)
            ntry   = ntry + 1
c                  write (*,'(6g12.4)')  radius,weight(radius),amassp(radius)
            If(reject.le.RANDD(Nseed))Goto10
                                                ! take this particle
            Coords(10,i) =amassp(radius)
            sMass = sMass +  Coords(10,i)
            phi    = 2.*pi*RANDD(Nseed)         ! create angles
            cs     = RANDD(Nseed)
            ss     = sign(asin(cs), RANDD(Nseed)-0.5)
            cs     = cos(ss)
            Coords(1,i) = radius*cs*cos(phi)    ! assign coordinates
            Coords(2,i) = radius*cs*sin(phi)
            Coords(3,i) = radius*sin(ss)

            sx = sx +Coords(1,i)**2*Coords(10,i)
            sy = sy +Coords(2,i)**2*Coords(10,i)
            sz = sz +Coords(3,i)**2*Coords(10,i)
            s1 = s1 +Coords(1,i)*Coords(10,i)
            s2 = s2 +Coords(2,i)*Coords(10,i)
            s3 = s3 +Coords(3,i)*Coords(10,i)
         EndDo 

         aMtot = dMass(dble(Rmaxinit))         ! total analytical mass
 
         Do i=1,Np               ! make velocities of particles
            Coords(10,i) =Coords(10,i)/sMass  ! renormalize masses to unity
            radius = sqrt(Coords(1,i)**2+Coords(2,i)**2+Coords(3,i)**2)
            sig_r =sqrt(dSigma(radius)/aMtot) ! normalized rms radial V

            do k=4,6
               Coords(k,i) = sig_r*GAUSS(Nseed)
            EndDo 
         EndDo 

         do i=1,np
          write(1,101)(coords(k,i),k=1,6),0.,0.,0.,coords(10,i)
         enddo

         sx = sqrt(sx/Smass)
         sy = sqrt(sy/Smass)
         sz = sqrt(sz/Smass)
         s1 = s1/Smass
         s2 = s2/Smass
         s3 = s3/Smass
         write (*,'(2(a,i8),a,g12.4)')' Done Init: Np=',Np,
     &                     ' Ntrials=',ntry,' Mass=',sMass
         write (*,'(a,3g12.4)')' Center =',s1,s2,s3
         write (*,'(a,3g12.4)')' Inertia=',sx,sy,sz
         aMassTotal = sMass

 101   format(10f18.10)

      Return
      End
c-------------------------------------------------------------------- 
C		
      FUNCTION density(x)
C                             density profile for NFW
c-------------------------------------------------------------------- 
      INCLUDE '../nbody1.h'
         density =1./(x/rs)/(1.+(x/rs))**2
      Return
      End
c-------------------------------------------------------------------- 
C		
      FUNCTION weight(x)
C                  weights for rejection of particles in Initial
c-------------------------------------------------------------------- 
      INCLUDE '../nbody1.h'
         weight =density(x)/amassp(x) *sqrt(x**5)
      Return
      End
c-------------------------------------------------------------------- 
C		
      FUNCTION amassp(x)
C                       masses of individual particles
c-------------------------------------------------------------------- 
      INCLUDE '../nbody1.h'
         amassp =(1.+(x/rs)**2)*1.e-5 
c         amassp =1.e-5 
      Return
      End
c-------------------------------------------------------------------- 
C		
      REAL*8 FUNCTION dSigma(x)
C                       radial velocity dispersion
c-------------------------------------------------------------------- 
      INCLUDE '../nbody1.h'
      REAL*8   dfSigma
      EXTERNAL dfSigma

      If(x.ge.Rmaxinit*0.999)Then
         dSigma =0.
         return
      Else
         CALL DGAUS8(dfSigma,dble(x),dble(Rmaxinit),1.d-5,dSigma,IERR)
         dSigma =dSigma /density(x)
      EndIf 
      Return
      End
c-------------------------------------------------------------------- 
C		
      REAL*8 FUNCTION dfSigma(x)
C                       integrant for dSigma
c-------------------------------------------------------------------- 
      INCLUDE '../nbody1.h'
      REAL*8 x,dMass
         dfSigma = dMass(x)/x**2*density(real(x))
      Return
      End
c-------------------------------------------------------------------- 
C		
      REAL*8 FUNCTION dMass(x)
C                       mass inside given radius
c-------------------------------------------------------------------- 
      INCLUDE '../nbody1.h'
      REAL*8 x,dfMass,ERR
      EXTERNAL dfMass

         CALL DGAUS8a(dfMass,1.d-10,x,1.d-5,dMass,IERR)
         dMass =pi4*dMass
      Return
      End
c-------------------------------------------------------------------- 
C		
      REAL*8 FUNCTION dfMass(x)
C                       integrant for dMass
c-------------------------------------------------------------------- 
      INCLUDE '../nbody1.h'
      REAL*8 x
         dfMass = x**2*density(real(x))
      Return
      End
C--------------------------------------
C		normal random numbers
      FUNCTION GAUSS(M)
C                    Uses RANDd() + 5 random numbers with correction
C                    Old code used in many simulations
C   N=1e8  sigma= 0.99556     mean= 0.10542E-03
c  sigma Frac(>sig) Frac(<-sig)   True       n(>)    n(<)
c  1.00 0.1589     0.1589     0.1587     15894145 15889594
c  1.50 0.6704E-01 0.6702E-01 0.6681E-01  6703736  6702430
c  2.00 0.2258E-01 0.2256E-01 0.2275E-01  2258038  2256030
c  2.50 0.5841E-02 0.5843E-02 0.6210E-02   584148   584264
c  3.00 0.1031E-02 0.1028E-02 0.1350E-02   103141   102783
c  4.00 0.5900E-06 0.6100E-06 0.3170E-04       59       61
C
C--------------------------------------
      X=0.
      DO  I=1,5
         X=X+RANDd(M)
      EndDo
      X2   =1.5491933*(X-2.5)
      GAUSS=X2*(1.-0.01*(3.-X2*X2))
      RETURN
      END
C--------------------------------------
C		normal random numbers
      FUNCTION GAUSS2(M)
C                        Uses ranndd for homogenous rand numbers
C                          sum of 12 random numbers + corrections
C                          It is  slow: very (10 times than GAUSS3)
C                          Quality is very good: upto 4 sigma
c  sigma Frac(>sig) Frac(<-sig)   True       n(>)    n(<)
c  1.00 0.1586     0.1586     0.1587     15857116 15857711
c  1.50 0.6678E-01 0.6680E-01 0.6681E-01  6678110  6679725
c  2.00 0.2277E-01 0.2275E-01 0.2275E-01  2277430  2274868
c  2.50 0.6197E-02 0.6189E-02 0.6210E-02   619732   618895
c  3.00 0.1327E-02 0.1319E-02 0.1350E-02   132713   131929
c  4.00 0.2664E-04 0.2727E-04 0.3170E-04     2664     2727
C                   
c-------------------------------------------------------------------- 
      X=0.
      DO  I=1,12
         X=X+Randd(M)
      EndDo
      X2 = X-6.0
      GAUSS2   =X2 *(1.-0.0045*(3.-X2*X2))
      RETURN
      END

C------------------------------------------------
C                                              random number generator
      FUNCTION RANDd(M)
C------------------------------------------------
      DATA LC,AM,KI,K1,K2,K3,K4,L1,L2,L3,L4
     +	/453815927,2147483648.,2147483647,536870912,131072,256,
     +	 16777216,4,16384,8388608,128/
      ML=M/K1*K1
      M1=(M-ML)*L1
      ML=M/K2*K2
      M2=(M-ML)*L2
      ML=M/K3*K3
      M3=(M-ML)*L3
      ML=M/K4*K4
      M4=(M-ML)*L4
      M5=KI-M
      IF(M1.GE.M5)M1=M1-KI-1
      ML=M+M1
      M5=KI-ML
      IF(M2.GE.M5)M2=M2-KI-1
      ML=ML+M2
      M5=KI-ML
      IF(M3.GE.M5)M3=M3-KI-1
      ML=ML+M3
      M5=KI-ML
      IF(M4.GE.M5)M4=M4-KI-1
      ML=ML+M4
      M5=KI-ML
      IF(LC.GE.M5)ML=ML-KI-1
      M=ML+LC
      RANDd=M/AM
      RETURN
      END
c-------------------------------------------------------------------- 
      SUBROUTINE DGAUS8 (FUN, A, B, ERR, ANS, IERR)
C***BEGIN PROLOGUE  DGAUS8
C***PURPOSE  Integrate a real function of one variable over a finite
C            interval using an adaptive 8-point Legendre-Gauss
C            algorithm.  Intended primarily for high accuracy
C            integration or integration of smooth functions.
C***LIBRARY   SLATEC
C***CATEGORY  H2A1A1
C***TYPE      DOUBLE PRECISION (GAUS8-S, DGAUS8-D)
C***KEYWORDS  ADAPTIVE QUADRATURE, AUTOMATIC INTEGRATOR,
C             GAUSS QUADRATURE, NUMERICAL INTEGRATION
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract  *** a DOUBLE PRECISION routine ***
C        DGAUS8 integrates real functions of one variable over finite
C        intervals using an adaptive 8-point Legendre-Gauss algorithm.
C        DGAUS8 is intended primarily for high accuracy integration
C        or integration of smooth functions.
C
C        The maximum number of significant digits obtainable in ANS
C        is the smaller of 18 and the number of digits carried in
C        double precision arithmetic.
C
C     Description of Arguments
C
C        Input--* FUN, A, B, ERR are DOUBLE PRECISION *
C        FUN - name of external function to be integrated.  This name
C              must be in an EXTERNAL statement in the calling program.
C              FUN must be a DOUBLE PRECISION function of one DOUBLE
C              PRECISION argument.  The value of the argument to FUN
C              is the variable of integration which ranges from A to B.
C        A   - lower limit of integration
C        B   - upper limit of integration (may be less than A)
C        ERR - is a requested pseudorelative error tolerance.  Normally
C              pick a value of ABS(ERR) so that DTOL .LT. ABS(ERR) .LE.
C              1.0D-3 where DTOL is the larger of 1.0D-18 and the
C              double precision unit roundoff D1MACH(4).  ANS will
C              normally have no more error than ABS(ERR) times the
C              integral of the absolute value of FUN(X).  Usually,
C              smaller values of ERR yield more accuracy and require
C              more function evaluations.
C
C              A negative value for ERR causes an estimate of the
C              absolute error in ANS to be returned in ERR.  Note that
C              ERR must be a variable (not a constant) in this case.
C              Note also that the user must reset the value of ERR
C              before making any more calls that use the variable ERR.
C
C        Output--* ERR,ANS are double precision *
C        ERR - will be an estimate of the absolute error in ANS if the
C              input value of ERR was negative.  (ERR is unchanged if
C              the input value of ERR was non-negative.)  The estimated
C              error is solely for information to the user and should
C              not be used as a correction to the computed integral.
C        ANS - computed value of integral
C        IERR- a status code
C            --Normal codes
C               1 ANS most likely meets requested error tolerance,
C                 or A=B.
C              -1 A and B are too nearly equal to allow normal
C                 integration.  ANS is set to zero.
C            --Abnormal code
C               2 ANS probably does not meet requested error tolerance.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, I1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   810223  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  DGAUS8
      INTEGER IERR, K, KML, KMX, L, LMN, LMX, LR, MXL, NBITS,
     1 NIB, NLMN, NLMX
      INTEGER I1MACH
      DOUBLE PRECISION A,AA,AE,ANIB,ANS,AREA,B,C,CE,EE,EF,
     1 EPS, ERR, EST, GL, GLR, GR, HH, SQ2, TOL, VL, VR, W1, W2, W3,
     2 W4, X1, X2, X3, X4, X, H,
     3 D1MACH5, D1MACH4
      DOUBLE PRECISION D1MACH, G8, FUN
      DIMENSION AA(60), HH(60), LR(60), VL(60), GR(60)
      SAVE X1, X2, X3, X4, W1, W2, W3, W4, SQ2,
     1 NLMN, KMX, KML
      DATA X1, X2, X3, X4/
     1     1.83434642495649805D-01,     5.25532409916328986D-01,
     2     7.96666477413626740D-01,     9.60289856497536232D-01/
      DATA W1, W2, W3, W4/
     1     3.62683783378361983D-01,     3.13706645877887287D-01,
     2     2.22381034453374471D-01,     1.01228536290376259D-01/
      DATA SQ2/1.41421356D0/
      DATA NLMN/1/,KMX/5000/,KML/6/
      G8(X,H)=H*((W1*(FUN(X-X1*H) + FUN(X+X1*H))
     1           +W2*(FUN(X-X2*H) + FUN(X+X2*H)))
     2          +(W3*(FUN(X-X3*H) + FUN(X+X3*H))
     3           +W4*(FUN(X-X4*H) + FUN(X+X4*H))))
C***FIRST EXECUTABLE STATEMENT  DGAUS8
C
C     Initialize
C
      D1MACH4 = 1.d-17 
c      K = I1MACH(14)
c      ANIB = D1MACH(5)*K/0.30102000D0
c      NBITS = ANIB
c      NLMX = MIN(60,(NBITS*5)/8)
      NBITS = 16 
      NLMX = MIN(60,(NBITS*5)/8)
      ANS = 0.0D0
      IERR = 1
      CE = 0.0D0
      IF (A .EQ. B) GO TO 140
      LMX = NLMX
      LMN = NLMN
      IF (B .EQ. 0.0D0) GO TO 10
      IF (SIGN(1.0D0,B)*A .LE. 0.0D0) GO TO 10
      C = ABS(1.0D0-A/B)
      IF (C .GT. 0.1D0) GO TO 10
      IF (C .LE. 0.0D0) GO TO 140
      ANIB = 0.5D0 - LOG(C)/0.69314718D0
      NIB = ANIB
      LMX = MIN(NLMX,NBITS-NIB-7)
      IF (LMX .LT. 1) GO TO 130
      LMN = MIN(LMN,LMX)
   10 TOL = MAX(ABS(ERR),2.0D0**(5-NBITS))/2.0D0

      IF (ERR .EQ. 0.0D0) TOL = SQRT(D1MACH4)
      EPS = TOL
      HH(1) = (B-A)/4.0D0
      AA(1) = A
      LR(1) = 1
      L = 1
      EST = G8(AA(L)+2.0D0*HH(L),2.0D0*HH(L))

      K = 8
      AREA = ABS(EST)
      EF = 0.5D0
      MXL = 0
C
C     Compute refined estimates, estimate the error, etc.
C
   20 GL = G8(AA(L)+HH(L),HH(L))
      GR(L) = G8(AA(L)+3.0D0*HH(L),HH(L))
      K = K + 16
      AREA = AREA + (ABS(GL)+ABS(GR(L))-ABS(EST))
C     IF (L .LT .LMN) GO TO 11
      GLR = GL + GR(L)
      EE = ABS(EST-GLR)*EF
      AE = MAX(EPS*AREA,TOL*ABS(GLR))
      IF (EE-AE) 40, 40, 50
   30 MXL = 1
   40 CE = CE + (EST-GLR)
      IF (LR(L)) 60, 60, 80
C
C     Consider the left half of this level
C
   50 IF (K .GT. KMX) LMX = KML
      IF (L .GE. LMX) GO TO 30
      L = L + 1
      EPS = EPS*0.5D0
      EF = EF/SQ2
      HH(L) = HH(L-1)*0.5D0
      LR(L) = -1
      AA(L) = AA(L-1)
      EST = GL
      GO TO 20
C
C     Proceed to right half at this level
C
   60 VL(L) = GLR
   70 EST = GR(L-1)
      LR(L) = 1
      AA(L) = AA(L) + 4.0D0*HH(L)
      GO TO 20
C
C     Return one level
C
   80 VR = GLR
   90 IF (L .LE. 1) GO TO 120
      L = L - 1
      EPS = EPS*2.0D0
      EF = EF*SQ2
      IF (LR(L)) 100, 100, 110
  100 VL(L) = VL(L+1) + VR
      GO TO 70
  110 VR = VL(L+1) + VR
      GO TO 90
C
C     Exit
C
  120 ANS = VR
      IF ((MXL.EQ.0) .OR. (ABS(CE).LE.2.0D0*TOL*AREA)) GO TO 140
      IERR = 2
c      write (*,*)  'DGAUS8: ',
c     +   'ANS is probably insufficiently accurate.'
      GO TO 140
  130 IERR = -1
      write (*,*)   'DGAUS8: ',
     +   'A and B are too nearly equal to allow normal integration.'
  140 IF (ERR .LT. 0.0D0) ERR = CE
      RETURN
      END
c-------------------------------------------------------------------- 
      SUBROUTINE DGAUS8a (FUN, A, B, ERR, ANS, IERR)
      INTEGER IERR, K, KML, KMX, L, LMN, LMX, LR, MXL, NBITS,
     1 NIB, NLMN, NLMX
      INTEGER I1MACH
      DOUBLE PRECISION A,AA,AE,ANIB,ANS,AREA,B,C,CE,EE,EF,
     1 EPS, ERR, EST, GL, GLR, GR, HH, SQ2, TOL, VL, VR, W1, W2, W3,
     2 W4, X1, X2, X3, X4, X, H,
     3 D1MACH5, D1MACH4
      DOUBLE PRECISION D1MACH, G8, FUN
      DIMENSION AA(60), HH(60), LR(60), VL(60), GR(60)
      SAVE X1, X2, X3, X4, W1, W2, W3, W4, SQ2,
     1 NLMN, KMX, KML
      DATA X1, X2, X3, X4/
     1     1.83434642495649805D-01,     5.25532409916328986D-01,
     2     7.96666477413626740D-01,     9.60289856497536232D-01/
      DATA W1, W2, W3, W4/
     1     3.62683783378361983D-01,     3.13706645877887287D-01,
     2     2.22381034453374471D-01,     1.01228536290376259D-01/
      DATA SQ2/1.41421356D0/
      DATA NLMN/1/,KMX/5000/,KML/6/
      G8(X,H)=H*((W1*(FUN(X-X1*H) + FUN(X+X1*H))
     1           +W2*(FUN(X-X2*H) + FUN(X+X2*H)))
     2          +(W3*(FUN(X-X3*H) + FUN(X+X3*H))
     3           +W4*(FUN(X-X4*H) + FUN(X+X4*H))))
C***FIRST EXECUTABLE STATEMENT  DGAUS8
C
C     Initialize
C
      D1MACH4 = 1.d-17 
      NBITS = 16 
      NLMX = MIN(60,(NBITS*5)/8)
      ANS = 0.0D0
      IERR = 1
      CE = 0.0D0
      IF (A .EQ. B) GO TO 140
      LMX = NLMX
      LMN = NLMN
      IF (B .EQ. 0.0D0) GO TO 10
      IF (SIGN(1.0D0,B)*A .LE. 0.0D0) GO TO 10
      C = ABS(1.0D0-A/B)
      IF (C .GT. 0.1D0) GO TO 10
      IF (C .LE. 0.0D0) GO TO 140
      ANIB = 0.5D0 - LOG(C)/0.69314718D0
      NIB = ANIB
      LMX = MIN(NLMX,NBITS-NIB-7)
      IF (LMX .LT. 1) GO TO 130
      LMN = MIN(LMN,LMX)
   10 TOL = MAX(ABS(ERR),2.0D0**(5-NBITS))/2.0D0

      IF (ERR .EQ. 0.0D0) TOL = SQRT(D1MACH4)
      EPS = TOL
      HH(1) = (B-A)/4.0D0
      AA(1) = A
      LR(1) = 1
      L = 1
      EST = G8(AA(L)+2.0D0*HH(L),2.0D0*HH(L))

      K = 8
      AREA = ABS(EST)
      EF = 0.5D0
      MXL = 0
C
C     Compute refined estimates, estimate the error, etc.
C
   20 GL = G8(AA(L)+HH(L),HH(L))
      GR(L) = G8(AA(L)+3.0D0*HH(L),HH(L))
      K = K + 16
      AREA = AREA + (ABS(GL)+ABS(GR(L))-ABS(EST))
C     IF (L .LT .LMN) GO TO 11
      GLR = GL + GR(L)
      EE = ABS(EST-GLR)*EF
      AE = MAX(EPS*AREA,TOL*ABS(GLR))
      IF (EE-AE) 40, 40, 50
   30 MXL = 1
   40 CE = CE + (EST-GLR)
      IF (LR(L)) 60, 60, 80
C
C     Consider the left half of this level
C
   50 IF (K .GT. KMX) LMX = KML
      IF (L .GE. LMX) GO TO 30
      L = L + 1
      EPS = EPS*0.5D0
      EF = EF/SQ2
      HH(L) = HH(L-1)*0.5D0
      LR(L) = -1
      AA(L) = AA(L-1)
      EST = GL
      GO TO 20
C
C     Proceed to right half at this level
C
   60 VL(L) = GLR
   70 EST = GR(L-1)
      LR(L) = 1
      AA(L) = AA(L) + 4.0D0*HH(L)
      GO TO 20
C
C     Return one level
C
   80 VR = GLR
   90 IF (L .LE. 1) GO TO 120
      L = L - 1
      EPS = EPS*2.0D0
      EF = EF*SQ2
      IF (LR(L)) 100, 100, 110
  100 VL(L) = VL(L+1) + VR
      GO TO 70
  110 VR = VL(L+1) + VR
      GO TO 90
C
C     Exit
C
  120 ANS = VR
      IF ((MXL.EQ.0) .OR. (ABS(CE).LE.2.0D0*TOL*AREA)) GO TO 140
      IERR = 2
c      write (*,*)  'DGAUS8: ',
c     +   'ANS is probably insufficiently accurate.'
      GO TO 140
  130 IERR = -1
      write (*,*)   'DGAUS8: ',
     +   'A and B are too nearly equal to allow normal integration.'
  140 IF (ERR .LT. 0.0D0) ERR = CE
      RETURN
      END
