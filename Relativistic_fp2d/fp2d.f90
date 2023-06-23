!C    FULLY RELATIVISTIC COLLISIONAL OPERATOR BRAAMS-KARNEY
!C    N1 NUMBER OF POINTS IN THE MOMENTUM P
!C    N2 NUMBER OF POINTS IN COS(THETA)
!C    TRAPO=1. TRAPPING INCLUDED IN THE CALCULATION
!C    TRAPO=0. TRAPPING EXCLUDED IN THE CALCULATION
!C    MODE=5 USE FIVE POINTS DIAGONAL SOLVER
!C    MODE=9 USE NINE POINTS DIAGONAL SOLVER
!C    APARAM PARAMETERS FOR CONVERGENCE (USUALLY 60)
!C    NITS NUMBER OF ITERATIONS NEEDED
!C    PMAX MAXIMUM MOMENTUM
!C    TRUNC=1 INCLUDES THE COLLISION TERM FOR MOMENTUM CONSERVATION
!C    PMIN MINIMUM VALUE OF THE MOMENTUM ,USUALLY ZERO.
!C    TE ELECTRON temperature in keV
!C    DT TIME STEP
!C    MINIMUM VALUE OF THE VELOCITY FOR THE LOWER HYBRID SPECTRUM
!C    MAXIMUM VALUE OF THE VELOCITY FOR THE LOWER HYBRID SPECTRUM
!C    DDIF QUASILINEAR DIFFUSION OPERATOR FOR THE LOWER HYBRID WAVE
!C    OMC RATIO OF THE CYCLOTRON FREQUENCY OVER THE WAVE FREQUENCY
!C    DPERP COEFFICIENT FOR THE ELECTRON CYCLOTRON DIFFUSION COEFFICIENT
!C    DNRES RESONANCE REFRACTIVE INDEX FOR THE ELECTRON CYCLOTRON WAVE
!C    DNPAR WIDTH OF THE ELECTRON CYCLOTRON RESONANCE
!C    T(I,J) DISTRIBUTION FUNCTION TO BE CALCULATED
!C    F(I) MAXWELLIAN DISTRIBUTION
!C    DCY(I,J) QUASILINEAR DIFFUSION COEFFICIENT FOR THE ELECTRON
!C    CYCLOTRON WAVE
!C    EPS ELECTRIC FIELD NORMALIZED TO THE DREICER FIELD
!C    pdens  electron density in 10**13/cm**3
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N1,N2,N1M,NITS,IFAIL,MODE
      PARAMETER(N1=1001,N2=1000,N1M=N1,ntime=110000)
      REAL*8 A(N1M,N2),B(N1M,N2),C(N1M,N2),D(N1M,N2),E(N1M,N2),Q(N1M,N2)
      REAL*8 R(N1M,N2),T(N1M,N2),WRKSP1(N1M,N2),WRKSP2(N1M,N2),XP(N1),YMU(N2)
      REAL*8 DIFF(N1,N2),T1(N1),F(N1),AP1(N1),AP2(N1),AP3(N1),AB1(N1),AB2(N1)
      REAL*8 AB3(N1),AB4(N1),AB5(N1),AB6(N1),AB7(N1),AB8(N1),AB9(N1),AB10(N1)
      REAL*8 AB11(N1),AB12(N1),AB13(N1),AB14(N1),SIG(N1),T2(N1),DFDP(N1,N2)
      REAL*8 DFDM(N1,N2),DE(N1),DF(N1),DFS(N1,N2),ER(N1,N2),DFW(N1,N2),BI(N1)
      REAL*8 APH(N1),FPH(N1),FMH(N1),AMH(N1),EJ0(N2),EJ2(N2),DIF(N1,N2),TRP(N2)
      REAL*8 TRPP(N2),ctime(ntime),curtim(ntime),pabsti(ntime)
      REAL*8 AM(N1M,N2),AP(N1M,N2),EM(N1M,N2),EP(N1M,N2),DCY(N1,N2)
      REAL*8 WRKSP3(N1M,N2),WRKSP4(N1M,N2),ALPHAMAX
!C    FULLY RELATIVISTIC COLLISIONAL OPERATOR BRAAMS-KARNEY
      REAL*8 S11ABF
      EXTERNAL S11ABF
!C     EXTERNAL D03UAF
      pi=3.14159265D0
!C    TRAPO=1. INCLUDE TRAPPING
	TRAPO=1.
!C     WRITE(*,*) 'MODE 5 or 9 ?'
!C     READ(6,*) MODE
	MODE=9
!C     WRITE(*,*) 'APARAM ?'
!C     READ(6,*) APARAM
	APARAM=1000. 
!C     WRITE(*,*) 'NITS ?'
!C     READ(6,*) NITS
	NITS=10000
!C     WRITE(*,*) 'PMAX ?'
!C     READ(6,*) PMAX
	PMAX=100.
      write(*,*),'pmax=',pmax
!C     WRITE(*,*) 'TRUNC E-E COLLISION (0 OR 1) ?'
!C     READ(6,*) TRUNC
	TRUNC=1.
      PMIN=0.
      DP=(PMAX-PMIN)/FLOAT(N1-1)
      DMU=2./FLOAT(N2-2-1)
!C     WRITE(*,*) 'TE ?'
!C     READ(6,*) TE
      
       pause
      TE=0.5
       write(*,*) 'te=',te
      BTH=0.04424*DSQRT(TE)
	TI=0.05
      TRAT=TE/TI
!C     WRITE(*,*) 'DT ?'
!C     READ(6,*) DTT
	DTT=0.05
      DT=DTT
!C      WRITE(*,*) 'E ?'
!C     READ(6,*) EEPS
	EEPS=0.01
      elfield=eeps
      COLLF=0.0
      Z=2.0
!C    RMIN=RADIUS OF THE FLUX SURFACE
      RMIN=5.
!C    RMAJ=MAJOR RADIUS OF THE TOKAMAK
      RMAJ=43.
      NPRINT=1
      MPRINT=0
      CC=DSQRT(1840.D0*TE/TI)
      CCC=DSQRT(1840.D0)
!C    PARAMETERS FOR THE LOWER HYBRID WAVE CURRENT DRIVE PROBLEM
!C    DIFFUSION MATRIX IN DIF(I,J)
!C     WRITE(*,*) 'VMIN ?'
!C     READ(6,*) V1
	V1=3.5
!C     V2=-3.5
!C     WRITE(*,*) 'VMAX ?'
!C     READ(6,*) V2
	V2=13.
!C     V1=-12.
      V1i=V1
!C     WRITE(*,*) 'DQL ?'
!C     READ(6,*) DDIF
	DDIF=0.0d0
      pdens=1.d0
!C    PARAMETERS FOR ELECTRON CYCLOTRON HEATING
      HARMON=1.
      COSTHETA=1.
      THE=DACOS(COSTHETA)
      DTHE=THE*180./pi
      WRITE(6,*) DTHE
      DPERP=0.000
!C    DELTA IS THE SHAFRANOV SHIFT
      DELTA=0.0
      OMC=1.07
      DNRES=0.2
      DNPAR=0.1
!C    PARAMETERS FOR THE TRAPPING EFFECT
      RRATIO=RMIN/RMAJ
      TRAP=DSQRT(2.*RRATIO/(1.+RRATIO))
      THE=DACOS(TRAP)
      DTHE=THE*180./pi
      WRITE(6,*) DTHE,THE
      TRAP2=TRAP*TRAP
	WRITE(6,*) 'TRAP=',TRAP,' TRAP2=',TRAP2
      RATIOB=(1.+DELTA/RMAJ+RRATIO)/(1.+DELTA/RMAJ+RRATIO*COSTHETA)
      OMC=OMC/(1.+DELTA/RMAJ+RRATIO*COSTHETA)
      IF(TRAPO.EQ.0.0D0) RATIOB=1.
	DTRA=DSQRT(1.D0/RATIOB)
	DTRA=DASIN(DTRA)
      WRITE(6,*) RATIOB,DTRA
      WRITE(6,*) OMC
!C     pause
      N1M1=N1-1
      N2M1=N2-1
      DO 1 I=1,N1
  1   XP(I)=DP*FLOAT(I-1)+PMIN
      DO 2 J=2,N2M1
  2   YMU(J)=-1.D0+DMU*FLOAT(J-2)
      YMU(1)=YMU(3)
      YMU(N2)=YMU(N2-2)
      CONST=2./DSQRT(pi)
      DO 46 J=2,N2M1
      ABSM=DABS(YMU(J))
      YMU2=YMU(J)*YMU(J)
      IF(ABSM.LE.TRAP) THEN
      ARG=YMU2/TRAP2
      ARG1=1.-ARG
      ELSE
      ARG=TRAP2/YMU2
      ARG1=1.-ARG
      ENDIF
      EK=1.38629436112+0.09666344259*ARG1+0.03590092383*(ARG1**2) &
       + 0.03742563713*(ARG1**3)+0.01451196212*(ARG1**4)+(0.5+ &
      .12498593597*ARG1+.06880248576*(ARG1**2)+.03328355346* &
      (ARG1**3)+0.00441787012*(ARG1**4))*DLOG(1.D0/ARG1)
      EE=1.+0.44325141463*ARG1+0.0626060122*(ARG1**2)+ &
      .04757383546*(ARG1**3)+.01736506451*(ARG1**4)+(.24998368310 &
      *ARG1+.09200180037*(ARG1**2)+ &
      .04069697526*(ARG1**3)+.00526449639*(ARG1**4))*DLOG(1.D0/ARG1)
      IF(ABSM.LE.TRAP) THEN
      EJ0(J)=EK*DABS(YMU(J)/TRAP)
      EJ2(J)=(EK-EE)*DABS(YMU(J)/TRAP)
      ELSE
      EJ0(J)=EK
      EJ2(J)=(EK-EE)*YMU2/TRAP2
      ENDIF
      EJJ=EJ0(J)-0.5*TRAP2*EJ2(J)
	TRPP(J)=EJJ*2.0/pi
!C	WRITE(6,*) J,TRPP(J)
      EJ2(J)=-TRAP2*(EJ2(J)-.5*(2.*(TRAP2+YMU2)*EJ2(J)-YMU2*EJ0(J))/3.)/EJJ
      EJ0(J)=EJJ
!C     IF(YMU2.GT.TRAP2) THEN
      IF(ABSM.LE.TRAP) THEN
      EJ0(J)=0.
      ELSE
      EJ0(J)=(pi/2.)*(1.+RRATIO)/(EJ0(J)*DSQRT(1.D0-RRATIO*RRATIO))
      ENDIF
!C    REMOVE TRAPPING
      IF(TRAPO.EQ.0.0D0) EJ0(J)=1.
      IF(TRAPO.EQ.0.0D0) EJ2(J)=0.
!CF
      IF(TRAPO.EQ.0.0D0) TRPP(J)=1.
  46  CONTINUE
      EJ0(1)=EJ0(3)
      EJ0(N2)=EJ0(N2-2)
      EJ2(1)=EJ2(3)
      EJ2(N2)=EJ2(N2-2)
	TRPP(1)=TRPP(3)
	TRPP(N2)=TRPP(N2-2)
      DO 3 I=1,N1
      DO 3 J=2,N2M1
      ABSM=DABS(YMU(J))
      YMUSQ=YMU(J)*YMU(J)
      GAMA=DSQRT(1.D0+XP(I)*XP(I)*BTH*BTH)
      T(I,J)=DEXP(-XP(I)*XP(I)/(1.+GAMA))/(2.*pi)**1.5
      F(I)=DEXP(-XP(I)*XP(I)/(1.+GAMA))/(2.*pi)**1.5
      if(F(i).lt.1.D-150) then
          F(i)=1.D-150
      end if
      CNSTT=1.D0-RATIOB*(1.D0-YMU(J)*YMU(J))
      IF(CNSTT.GT.0.0.AND.CNSTT.LE.1.) THEN
      CNSTT=DSQRT(CNSTT)
      ELSE
      CNSTT=0.
      ENDIF
      IF(YMU(J).LT.0.0) CNSTT=-CNSTT
      PPAR=XP(I)*CNSTT
      PPERP=XP(I)*XP(I)-PPAR*PPAR
      IF(PPAR.EQ.0.0) THEN
      RESN=10.
      DCY(I,J)=0.
      ELSE
      RESN=(GAMA-OMC)/(PPAR*BTH)
      DCY(I,J)=PPERP*DPERP*DEXP(-(RESN-DNRES)**2/DNPAR**2)
      IF(CNSTT.NE.0.) CNSTA=YMU(J)/(CNSTT*TRPP(J))
      IF(TRAPO.GT.0.) DCY(I,J)=DCY(I,J)*CNSTA
      ENDIF
      VPAR=PPAR/GAMA
      DIF(I,J)=0.
      IF(VPAR.GE.V1.AND.VPAR.LE.V2) DIF(I,J)=DDIF
  3   CONTINUE
      DO 44 I=1,N1
      T(I,1)=T(I,3)
      T(I,N2)=T(I,N2-2)
	DCY(I,1)=DCY(I,3)
	DCY(I,N2)=DCY(I,N2-2)
      Q(I,N2)=0.
      Q(I,1)=0.
      A(I,N2)=0.
      B(I,N2)=0.
      C(I,N2)=0.
      D(I,N2)=0.0
      E(I,N2)=0.
      A(I,1)=0.
      B(I,1)=0.
      C(I,1)=0.
      D(I,1)=0.
!C9 diagonales
  	  AM(I,1)=0.
  	  AP(I,1)=0.
  	  EM(I,1)=0.
  	  EP(I,1)=0.
  44  E(I,1)=0.
      DO 45 J=1,N2
      Q(N1,J)=0.
      Q(1,J)=0.
      A(N1,J)=0.
      B(N1,J)=0.
      C(N1,J)=0.
      D(N1,J)=0.
      E(N1,J)=0.
!C9 diagonales
  	  AM(N1,J)=0.
  	  AP(N1,J)=0.
  	  EM(N1,J)=0.
  	  EP(N1,J)=0.
      A(1,J)=0.0
      B(1,J)=0.
      C(1,J)=0.
      D(1,J)=0.
!C9 diagonales
      AM(1,J)=0.
  	  AP(1,J)=0.
  	  EM(1,J)=0.
  	  EP(1,J)=0.
  45  E(1,J)=0.
      A(1,1)=0.0
      B(1,1)=0.0
      C(1,1)=0.0
      D(1,1)=0.0
      E(1,1)=0.0
!C9 diagonales
      AM(1,1)=0.
  	  AP(1,1)=0.
  	  EM(1,1)=0.
  	  EP(1,1)=0.
      IIFAIL=0.
      DO 75 I=1,N1
      XSIG=XP(I)*BTH
      SIG(I)=S11ABF(XSIG,IIFAIL)
  75  CONTINUE
      AP1(1)=0.
      AP2(1)=0.
      AP3(N1)=0.
!C    DO 54 I=2,N1
!C    GAMA=DSQRT(1.+XP(I)*XP(I)*BTH*BTH)
!C    V=XP(I)/GAMA
!C    AP1(I)=AP1(I-1)+XP(I)*F(I)*DP*V*4.*pi
!C    AP2(I)=AP2(I-1)+XP(I)*F(I)*DP*V*4.*pi*(1.-
!C    1GAMA*SIG(I)/(BTH*XP(I)))
!C 54  CONTINUE
      AP1(2)=4.*pi*(F(1)+F(2))*(DP**3)/8.
      AP1(3)=4.*pi*(F(1)*XP(1)**2/3.+4.*F(2)*XP(2)**2/3. + F(3)*XP(3)**2/3.)*DP
      AP1(4)=4.*pi*(3.*F(1)*XP(1)**2/8.+9.*F(2)*XP(2)**2/8. + 9.*F(3)*XP(3)**2/8.)*DP
      AP2(2)=4.*pi*BTH*BTH*(F(1)+F(2))*(DP**5)/(12.*5.)
      AP2(3)=4.*pi*BTH*BTH*(F(1)*XP(1)**4/3.+4.*F(2)*XP(2)**4/3.+ F(3)*XP(3)**4/3.)*DP/6.
      AP2(4)=4.*pi*BTH*BTH*(3.*F(1)*XP(1)**4/8.+9.*F(2) * XP(2)**4/8.+9.*F(3)*XP(3)**4/8.)*DP/6.
      AP1(1)=0.
      AP2(1)=0.
      AP3(N1)=0.
      DO 54 I=3,N1
      GAMA=DSQRT(1.+XP(I)*XP(I)*BTH*BTH)
      GAMA1=DSQRT(1.+XP(I-1)*XP(I-1)*BTH*BTH)
      V=XP(I)/GAMA
      V1=XP(I-1)/GAMA1
      AP1(I)=AP1(I-1)+0.5*XP(I)*F(I)*DP*V*4.*pi+0.5*XP(I-1)*F(I-1)*DP*V1*4.*pi
      AP2(I)=AP2(I-1)+0.5*XP(I)*F(I)*DP*V*4.*pi*(1.-GAMA*SIG(I)/(BTH*XP(I))) &
      +0.5*XP(I-1)*F(I-1)*DP*V1*4.*pi*(1.-GAMA1*SIG(I-1)/(BTH*XP(I-1)))
 54   CONTINUE
      DO 58 I=1,N1M1
      K=N1-I
      AP3(K)=AP3(K+1)+0.5*XP(K)*F(K)*DP*4.*pi + 0.5*XP(K+1)*F(K+1)*DP*4.*pi
  58  CONTINUE
      AB1(1)=0.
      AB2(1)=0.
      AB3(1)=0.
      AB4(1)=0.
      AB5(1)=0.
      AB6(1)=0.
      AB7(1)=0.
      AB1(2)=0.
      AB1(2)=4.*pi*(F(1)+F(2))*(DP**3)/8.
      AB2(2)=0.
      AB2(2)=4.*pi*(F(1)+F(2))*(DP**5)/(32.)
      AB3(2)=0.
      AB4(1)=0.
      AB5(1)=0.
      AB4(2)=0.
      AB5(2)=0.
      DO 59 I=3,N1
      GAMA=DSQRT(1.+XP(I)*XP(I)*BTH*BTH)
      V=XP(I)/GAMA
      GAMA1=DSQRT(1.+XP(I-1)*XP(I-1)*BTH*BTH)
      V1=XP(I-1)/GAMA1
      AB1(I)=AB1(I-1)+0.5*DP*F(I)*XP(I)*XP(I)*4.*pi &
     +0.5*DP*F(I-1)*XP(I-1)*XP(I-1)*4.*pi
      AB2(I)=AB2(I-1)+0.5*DP*F(I)*XP(I)*XP(I)*XP(I)*XP(I)*4.*pi &
     +0.5*DP*F(I-1)*XP(I-1)*XP(I-1)*XP(I-1)*XP(I-1)*4.*pi
      AB3(I)=AB3(I-1)+.5*DP*F(I)*XP(I)*XP(I)*4.*pi*(-3.*GAMA+SIG(I) &
      *(3./(BTH*XP(I))+2.*BTH*XP(I)))/GAMA &
     +0.5*DP*F(I-1)*XP(I-1)*XP(I-1)*4.*pi*(-3.*GAMA1+SIG(I-1) &
     *(3./(BTH*XP(I-1))+2.*BTH*XP(I-1)))/GAMA1
      AB4(I)=AB4(I-1)+0.5*DP*F(I)*XP(I)*XP(I)*4.*pi*(GAMA-SIG(I) &
     /(BTH*XP(I))-2.*GAMA*BTH*BTH*XP(I)*XP(I)/3.)/GAMA &
      +0.5*DP*F(I-1)*XP(I-1)*XP(I-1)*4.*pi*(GAMA1-SIG(I-1) &
     /(BTH*XP(I-1))-2.*GAMA1*BTH*BTH*XP(I-1)*XP(I-1)/3.)/GAMA1
      AB5(I)=AB5(I-1) &
     +.5*DP*F(I)*4.*pi*XP(I)*XP(I)*(1.-SIG(I)/(GAMA*BTH*XP(I))) &
     +.5*DP*F(I-1)*4.*pi*XP(I-1)*XP(I-1)*(1.-SIG(I-1)/(GAMA1*BTH*XP(I-1)))
  59  CONTINUE
      AB8(N1)=0.
      AB9(N1)=0.
      AB10(N1)=0.
      AB8(1)=0.
      AB9(1)=0.
      AB10(1)=0.
      N1M2=N1M1-1
      DO 66 I=1,N1M2
      K=N1-I
      GAMA=DSQRT(1.+XP(K)*XP(K)*BTH*BTH)
      V=XP(K)/GAMA
      GAMA1=DSQRT(1.+XP(K+1)*XP(K+1)*BTH*BTH)
      V1=XP(K+1)/GAMA1
      AB8(K)=AB8(K+1)+0.5*DP*XP(K)*XP(K)*F(K)*4.*pi/V &
      +0.5*DP*XP(K+1)*XP(K+1)*F(K+1)*4.*pi/V1
      AB9(K)=AB9(K+1)+0.5*DP*XP(K)*XP(K)*F(K)*4.*pi/(V*GAMA*GAMA) &
      +.5*DP*XP(K+1)*XP(K+1)*F(K+1)*4.*pi/(V1*GAMA1*GAMA1)
      AB10(K)=AB10(K+1)+0.5*DP*F(K)*XP(K)*XP(K)*4.*pi*V &
      +0.5*DP*F(K+1)*XP(K+1)*XP(K+1)*4.*pi*V1
  66  CONTINUE
      AB8(1)=AB8(2)
      AB9(1)=AB9(2)
      AB10(1)=AB10(2)
      IIFAIL=0
      DO 400 I=2,N1M1
      PPH=XP(I)+DP/2.
      PPH2=PPH*PPH
      PMH=XP(I)-DP/2.
      PMH2=PMH*PMH
      GAMA=DSQRT(1.D0+XP(I)*XP(I)*BTH*BTH)
      V=XP(I)/GAMA
      GAMAPH=DSQRT(1.D0+PPH*PPH*BTH*BTH)
      VPH=PPH/GAMAPH
      GAMAMH=DSQRT(1.D0+PMH*PMH*BTH*BTH)
      VMH=PMH/GAMAMH
      XSIG=PPH*BTH
      SIGPH=S11ABF(XSIG,IIFAIL)
      XSIG=PMH*BTH
      SIGMH=S11ABF(XSIG,IIFAIL)
      S=V/DSQRT(2.D0)
      SPH=VPH/DSQRT(2.D0)
      SMH=VMH/DSQRT(2.D0)
      EXPPH=DEXP(-PPH2/(1.+GAMAPH))/(2.*pi)**1.5
      FPH(I)=(AP1(I)+PPH*EXPPH*VPH*0.5*DP*4.*pi)/(VPH**2)+(AP2(I) &
      +0.5*DP*EXPPH*PPH*VPH*4.*pi*(1.-SIGPH*GAMAPH/(BTH*PPH)))/PPH2 &
      +(AP3(I)-PPH*EXPPH*(DP/2.)*4.*pi)*(1.-SIGPH/(BTH*PPH*GAMAPH))/VPH
      APH(I)=FPH(I)/VPH+Z/(VPH**3*CC*CC)
      IF(I.EQ.2) APH(2)=4.*pi*F(2)/3.+Z/(VPH**3*CC*CC)
      FPH(I)=FPH(I)+Z*TE/(TI*VPH*VPH*CC*CC)
      EXPMH=DEXP(-PMH2/(1.+GAMAMH))/(2.*pi)**1.5
      FMH(I)=(AP1(I)-PMH*EXPMH*VMH*0.5*DP*4.*pi)/(VMH**2)+(AP2(I) &
      -0.5*DP*EXPMH*PMH*VMH*4.*pi*(1.-SIGMH*GAMAMH/(BTH*PMH)))/PMH2 &
      +(AP3(I)+PMH*EXPMH*(DP/2.)*4.*pi)*(1.-SIGMH/(BTH*PMH*GAMAMH))/VMH
      AMH(I)=FMH(I)/VMH+Z/(VMH**3*CC*CC)
      IF(I.EQ.2) AMH(2)=4.*pi*F(2)/3.+Z/(VPH**3*CC*CC)
      FMH(I)=FMH(I)+Z*TE/(TI*VMH*VMH*CC*CC)
      FU2=1.-1./(2.*S*S)
      FU2C=1.-1./(2.*S*S*CC*CC)
      FU=(S*DEXP(-S*S)+Z*S*DEXP(-S*S*CC*CC)/CC)*CONST
      BI(I)=AB1(I)/(2.*V)-AB2(I)/(6.*XP(I)*XP(I)*V) &
      +AB3(I)/(BTH*BTH*XP(I)*XP(I)*V*GAMA*GAMA*8.)
      BI(I)=BI(I)-AB4(I) &
      /(V*4.*BTH*BTH*XP(I)*XP(I))-AB5(I)/(V*4.*GAMA*GAMA)+AB8(I)/2.
      BI(I)=BI(I)+AB9(I)*(-GAMA*GAMA/6.)
      BI(I)=BI(I)+AB9(I)*(-3.*GAMA+SIG(I)* &
      (3./(BTH*XP(I))+2.*BTH*XP(I)))/(8.*GAMA*BTH*BTH*XP(I)*XP(I))
      BI(I)=BI(I)-AB9(I)*GAMA*(GAMA-SIG(I)/(BTH*XP(I)) &
     -2.*GAMA*BTH*BTH*XP(I)*XP(I)/3.)/(4.*BTH*BTH*XP(I)*XP(I))
      BI(I)=BI(I)-AB10(I)*(GAMA-SIG(I)/(BTH*XP(I)))/(4.*GAMA*XP(I) &
     *XP(I))+Z*(1.-1./(CC*CC*V*V))/(2.*V)
 400  CONTINUE
      IIFAIL=0
      DO 701 J=1,N2
      DIF(1,J)=0.
 701  CONTINUE
      DO 702 I=2,N1
      DIF(I,1)=DIF(I,3)
      DIF(I,N2)=DIF(I,N2-2)
 702  CONTINUE
      IT=1
      itpr=1
!C    IF DIFF(I,J) IS TIME DEPENDENT, STARTS 30 CONTINUE HERE
!C  30 CONTINUE
      DO 707 J=1,N2
      DO 707 I=1,N1
      DIFF(I,J)=DIF(I,J)
 707  CONTINUE
      DO 4 I=2,N1M1
      PPH=XP(I)+DP/2.
      PPH2=PPH*PPH
      PMH=XP(I)-DP/2.
      PMH2=PMH*PMH
      GAMA=DSQRT(1.D0+XP(I)*XP(I)*BTH*BTH)
      GAMAB=GAMA*RATIOB
      V=XP(I)/GAMA
      GAMAPH=DSQRT(1.D0+PPH*PPH*BTH*BTH)
	  GAMAPHB=GAMAPH*RATIOB
      VPH=PPH/GAMAPH
      GAMAMH=DSQRT(1.D0+PMH*PMH*BTH*BTH)
	  GAMAMHB=GAMAMH*RATIOB
      VMH=PMH/GAMAMH
      XSIG=PPH*BTH
      SIGPH=S11ABF(XSIG,IIFAIL)
      XSIG=PMH*BTH
      SIGMH=S11ABF(XSIG,IIFAIL)
      S=V/DSQRT(2.D0)
      SPH=VPH/DSQRT(2.D0)
      SMH=VMH/DSQRT(2.D0)
      DO 4 J=2,N2M1
      CNSTT=1.D0-RATIOB*(1.D0-YMU(J)*YMU(J))
      IF(CNSTT.GT.0.0.AND.CNSTT.LE.1.) THEN
      CNSTT=DSQRT(CNSTT)
      ELSE
      CNSTT=0.
      ENDIF
      DABMU=CNSTT
      IF(YMU(J).LT.0.0) CNSTT=-CNSTT
      YPH=YMU(J)+DMU/2.
      IF(J.EQ.N2M1) YPH=YMU(N2M1)-DMU/2.
      DYPH=(1.D0-RATIOB*(1.D0-YPH*YPH))
      IF(DYPH.GT.0.0.AND.DYPH.LE.1.) THEN
        DYPH=DSQRT(DYPH)
      ELSE
        DYPH=0.
      ENDIF
      YMH=YMU(J)-DMU/2.
      IF(J.EQ.2) YMH=YMU(2)+DMU/2.
      DYMH=(1.D0-RATIOB*(1.D0-YMH*YMH))
      IF(DYMH.GT.0.0.AND.DYMH.LE.1.) THEN
        DYMH=DSQRT(DYMH)
      ELSE
        DYMH=0.
      ENDIF
      FBI=1.+EJ2(J)
      EPS=EEPS*EJ0(J)
      FBI1=1.+0.5*(EJ2(J)+EJ2(J-1))
      EPS1=EEPS*(EJ0(J)+EJ0(J-1))*0.5
      FBI2=1.+0.5*(EJ2(J)+EJ2(J+1))
      EPS2=EEPS*(EJ0(J)+EJ0(J+1))*0.5
      IF(DYMH.NE.0.) THEN
        WW=-EPS1*DMU*XP(I)/(BI(I)*FBI1+.5*(DIFF(I,J)+DIFF(I,J-1))*(1.- &
        YMH*YMH)+0.5*(DCY(I,J)+DCY(I,J-1))*(YMH*GAMAB-(GAMAB-OMC)/YMH)**2 &
        !C    2/(GAMA*DABS(YMH)*XP(I)))
        /(GAMAB*DYMH*XP(I)))
      ELSE
        WW=-EPS1*DMU*XP(I)/(BI(I)*FBI1+.5*(DIFF(I,J)+DIFF(I,J-1))*(1.-YMH*YMH))
      ENDIF
      IF(WW.NE.0.D0) THEN
        DELJ1=1./WW-1./(DEXP(WW)-1.)
      ELSE
        DELJ1=0.5
      END IF
      IF(DYPH.NE.0.) THEN
        WW=-EPS2*DMU*XP(I)/(BI(I)*FBI2+.5*(DIFF(I,J)+DIFF(I,J+1))*(1.- &
        YPH*YPH)+0.5*(DCY(I,J)+DCY(I,J+1))*(YPH*GAMAB-(GAMAB-OMC)/YPH)**2 &
        !C    2/(GAMA*DABS(YPH)*XP(I)))
        /(GAMAB*DYPH*XP(I)))
      ELSE
        WW=-EPS2*DMU*XP(I)/(BI(I)*FBI2+.5*(DIFF(I,J)+DIFF(I,J+1))* &
        (1.-YPH*YPH))
      ENDIF
      IF(WW.NE.0.D0) THEN
        DELJ2=1./WW-1./(DEXP(WW)-1.)
      ELSE
        DELJ2=0.5
      END IF
      BF=FPH(I)-EPS*YMU(J)
      IF(DABMU.NE.0) THEN
        CF=APH(I)+0.5*(DIFF(I,J)+DIFF(I+1,J))*YMU(J)*YMU(J) &
        +0.5*(DCY(I,J)+DCY(I+1,J))*GAMAPHB*(1.-YMU(J)*YMU(J))/(PPH*DABMU)
      ELSE
        CF=APH(I)+0.5*(DIFF(I,J)+DIFF(I+1,J))*YMU(J)*YMU(J)
      ENDIF
      WW=DP*BF/CF
      DELP=1./WW-1./(DEXP(WW)-1.)
      BF=FMH(I)-EPS*YMU(J)
      IF(DABMU.NE.0.) THEN
        CF=AMH(I)+0.5*(DIFF(I,J)+DIFF(I-1,J))*YMU(J)*YMU(J) &
        +0.5*(DCY(I,J)+DCY(I-1,J))*GAMAMHB*(1.-YMU(J)*YMU(J))/(PMH*DABMU)
      ELSE
        CF=AMH(I)+0.5*(DIFF(I,J)+DIFF(I-1,J))*YMU(J)*YMU(J)
      ENDIF
      WW=DP*BF/CF
      DELM=1./WW-1./(DEXP(WW)-1.)
      CMH=0.5*(DIFF(I,J)+DIFF(I,J-1))*(1.-YMH*YMH)/DMU
      CPH=0.5*(DIFF(I,J+1)+DIFF(I,J))*(1.-YPH*YPH)/DMU
      CDIFF=XP(I)*(DIFF(I,J+1)*YMU(J+1)*(1.-YMU(J+1)*YMU(J+1))*TRPP(J+1) &
      -DIFF(I,J-1)*YMU(J-1)*(1.-YMU(J-1)*YMU(J-1))*TRPP(J-1)) &
      /(4.*DP*DMU*TRPP(J))
      CNTE=(DIFF(I,J)+XP(I)*0.5*(DIFF(I+1,J)-DIFF(I-1,J))/DP)/(2.*DMU)
      CNSTT=1.D0-RATIOB*(1.D0-YMU(J-1)*YMU(J-1))
      IF(CNSTT.GT.0.0.AND.CNSTT.LE.1.) THEN
          CNSTT=DSQRT(CNSTT)
      ELSE
          CNSTT=0.
      ENDIF
      DABM1=CNSTT
      !C    DABM1=DABS(YMU(J-1))
      CNSTT=1.D0-RATIOB*(1.D0-YMU(J+1)*YMU(J+1))
      IF(CNSTT.GT.0.0.AND.CNSTT.LE.1.) THEN
          CNSTT=DSQRT(CNSTT)
      ELSE
          CNSTT=0.
      ENDIF
      DABP1=CNSTT
      DCYPH=(DCY(I+1,J)+DCY(I,J))
      DCYMH=(DCY(I,J)+DCY(I-1,J))
      IF(DABMU.NE.0.) THEN
          CNTCY1=-(1.-YMU(J)*YMU(J))*(-YMU(J)*PPH*0.5*(DIFF(I+1,J)+ &
          DIFF(I,J))+.5*DCYPH*(OMC+GAMAPHB*(YMU(J)*YMU(J)-1.))/(YMU(J)* &
          DABMU))/(4.*DP*DMU)-(1.-YMU(J)*YMU(J))*(YMU(J)*PMH*0.5*(DIFF(I,J)+ &
          DIFF(I-1,J))-0.5*DCYMH*(OMC+GAMAMHB*(YMU(J)*YMU(J)-1.))/(YMU(J)* &
          DABMU))/(4.*DP*DMU)
      ELSE
          CNTCY1=-(1.-YMU(J)*YMU(J))*(-YMU(J)*PPH*0.5*(DIFF(I+1,J)+ &
          DIFF(I,J)))-(1.-YMU(J)*YMU(J))*YMU(J)*PMH*0.5*(DIFF(I,J)+ &
          DIFF(I-1,J))
      ENDIF
!CF    EC TERMS NOT CROSSED
	  DLAMM1=0.5*(TRPP(J)+TRPP(J-1))/TRPP(J)
	  DLAMP1=0.5*(TRPP(J)+TRPP(J+1))/TRPP(J)
      IF(DYMH.NE.0.) THEN
          CNTY2=(1.-YMH*YMH)*DLAMM1*0.5*(DCY(I,J-1)+DCY(I,J))* &
          (YMH*GAMAB-(GAMAB-OMC)/YMH)**2/(GAMAB*XP(I)*DMU*DMU*DYMH)
      ELSE
          CNTY2=0.
      ENDIF
      A(I,J)=-(1.-YMH*YMH)*DLAMM1*(FBI1*BI(I)/DMU+EPS1*XP(I)*DELJ1+CMH)/DMU + CNTCY1-CNTY2
      !CM     1+YMU(J)*(1.-YMU(J)*YMU(J))*CNTE+CNTCY1-CNTY2
      IF(DYPH.NE.0..AND.DYMH.NE.0.)THEN
          CY1=(1.-YPH*YPH)*DLAMP1*(-XP(I)*YPH*0.5*(DIFF(I,J+1)+DIFF(I,J)) &
          +.5*(DCY(I,J+1)+DCY(I,J))*(OMC+GAMAB*(YPH*YPH-1.))/(YPH*DYPH)) /(4.*DP*DMU) &
          +(1.-YMH*YMH)*DLAMM1*(XP(I)*YMH*0.5*(DIFF(I,J)+DIFF(I,J-1))-0.5* &
          (DCY(I,J)+DCY(I,J-1))*(OMC+GAMAB*(YMH*YMH-1.))/(YMH*DYMH))/(4.*DP*DMU)
      ELSEIF (DYPH.EQ.0..AND.DYMH.NE.0.) THEN
           CY1=(1.-YPH*YPH)*DLAMP1*(-XP(I)*YPH*0.5*(DIFF(I,J+1)+DIFF(I,J))) &
           /(4.*DP*DMU)+(1.-YMH*YMH)*DLAMM1*(XP(I)*YMH*0.5*(DIFF(I,J)+ &
           DIFF(I,J-1))-0.5*(DCY(I,J)+DCY(I,J-1))*(OMC+GAMAB*(YMH*YMH-1.)) &
           /(YMH*DYMH))/(4.*DP*DMU)
      ELSEIF(DYPH.NE.0..AND.DYMH.EQ.0.)THEN
           CY1=(1.-YPH*YPH)*DLAMP1*(-XP(I)*YPH*0.5*(DIFF(I,J+1)+DIFF(I,J)) &
           +.5*(DCY(I,J+1)+DCY(I,J))*(OMC+GAMAB*(YPH*YPH-1.))/(YPH*DYPH)) &
           /(4.*DP*DMU) &
           +(1.-YMH*YMH)*DLAMM1*(XP(I)*YMH*0.5*(DIFF(I,J)+DIFF(I,J-1))) &
           /(4.*DP*DMU)
      ELSE
           CY1=(1.-YPH*YPH)*DLAMP1*(-XP(I)*YPH*0.5*(DIFF(I,J+1)+DIFF(I,J))) &
           /(4.*DP*DMU)+(1.-YMH*YMH)*DLAMM1*XP(I)*YMH*0.5*(DIFF(I,J)+ &
           DIFF(I,J-1))/(4.*DP*DMU)
      ENDIF
      IF(DABMU.NE.0.) THEN
          CY2=PMH*.5*DCYMH*GAMAMHB*(1.-YMU(J)*YMU(J))/(DABMU*DP*DP)
      ELSE
          CY2=0.
	  ENDIF
      B(I,J)=-PMH2*YMU(J)*EPS*DELM/DP-(DIFF(I,J)+DIFF(I-1,J))*YMU(J)*YMU(J) &
      *PMH2/(2.*DP*DP)-(AMH(I)/DP-FMH(I)*DELM)*PMH2/DP-CY1-CY2
      !  2*PMH2/(2.*DP*DP)-(AMH(I)/DP-FMH(I)*DELM)*PMH2/DP+CDIFF-CY1-CY2
      CNTCC=+((DIFF(I+1,J)+DIFF(I,J))*PPH2+(DIFF(I,J) + DIFF(I-1,J))*PMH2)*YMU(J)*YMU(J)/(2.*DP*DP)
      CNTC=CNTCC-(1.-YPH*YPH)*DLAMP1*(-FBI2*BI(I)/DMU-EPS2*XP(I)*DELJ2-CPH)/DMU &
      +(1.-YMH*YMH)*DLAMM1*(FBI1*BI(I)/DMU-EPS1*XP(I)*(1.-DELJ1)+CMH)/DMU
      CYY=0.
      IF(DABMU.NE.0.) THEN
          if (DYMH.NE.0.AND.DYPH.NE.0.) then
                   CYY=PPH*.5*DCYPH*GAMAPHB*(1.-YMU(J)*YMU(J))/(DABMU*DP*DP) &
                   +PMH*.5*DCYMH*GAMAMHB*(1.-YMU(J)*YMU(J))/(DABMU*DP*DP) &
                   +(1.-YMH*YMH)*DLAMM1*0.5*(DCY(I,J-1)+DCY(I,J))* &
                  (YMH*GAMAB-(GAMAB-OMC)/YMH)**2/(GAMAB*XP(I)*DMU*DMU*DYMH) &
                  +(1.-YPH*YPH)*DLAMP1*0.5*(DCY(I,J+1)+DCY(I,J))* &
                  (YPH*GAMAB-(GAMAB-OMC)/YPH)**2/(GAMAB*XP(I)*DMU*DMU*DYPH)
		  endif
      ELSE
        CYY=0.
      ENDIF
      C(I,J)=EPS*YMU(J)*(PPH2*DELP-PMH2*(1.-DELM))/DP &
      +CNTC-PPH2*(-APH(I)/DP+FPH(I)*DELP)/DP &
      +PMH2*(AMH(I)/DP+FMH(I)*(1.-DELM))/DP+CYY+XP(I)*XP(I)/DT
      IF(DABMU.NE.0.) THEN
        CY22=PPH*.5*DCYPH*GAMAPHB*(1.-YMU(J)*YMU(J))/(DABMU*DP*DP)
      ELSE
        CY22=0.
      ENDIF
      D(I,J)=PPH2*YMU(J)*(EPS*(1.-DELP)-(DIFF(I+1,J)+DIFF(I,J))*YMU(J) &
      /(2.*DP))/DP-PPH2*(APH(I)/DP+FPH(I)*(1.-DELP))/DP+CY1-CY22
      IF(DYPH.NE.0.) THEN
          CNTY22=(1.-YPH*YPH)*DLAMP1*0.5*(DCY(I,J+1)+DCY(I,J))* &
          (YPH*GAMAB-(GAMAB-OMC)/YPH)**2/(GAMAB*XP(I)*DMU*DMU*DYPH)
      ELSE
      CNTY22=0.
      ENDIF
      E(I,J)=-(1.-YPH*YPH)*DLAMP1*(FBI2*BI(I)/DMU-EPS2*XP(I)*(1.-DELJ2)+CPH)/DMU-CNTCY1-CNTY22
      !C    1/DMU-YMU(J)*(1.-YMU(J)*YMU(J))*CNTE-CNTCY1-CNTY22
      IF(DABMU.NE.0..AND.DYMH.NE.0.) THEN
        AM(I,J)=-(1.-YMU(J)*YMU(J))*(YMU(J)*PMH*0.5*(DIFF(I-1,J)+ &
        DIFF(I,J))-0.5*DCYMH*(OMC+GAMAMHB*(YMU(J)*YMU(J)-1.))/(YMU(J)* &
        DABMU))/(4.*DP*DMU)-DLAMM1*(1.-YMH*YMH)*(XP(I)*YMH*0.5*(DIFF(I,J)+ &
        DIFF(I,J-1))-0.5*(DCY(I,J)+DCY(I,J-1))*(OMC+GAMAB*(YMH*YMH-1.))/ &
        (YMH*DYMH))/(4.*DP*DMU)
        AP(I,J)=-(1.-YMU(J)*YMU(J))*(-YMU(J)*PPH*0.5*(DIFF(I+1,J)+ &
        DIFF(I,J))+0.5*DCYPH*(OMC+GAMAPHB*(YMU(J)*YMU(J)-1.))/(YMU(J)* &
        DABMU))/(4.*DP*DMU)+DLAMM1*(1.-YMH*YMH)*(XP(I)*YMH*0.5*(DIFF(I,J)+ &
        DIFF(I,J-1))-0.5*(DCY(I,J)+DCY(I,J-1))*(OMC+GAMAB*(YMH*YMH-1.))/ &
        (YMH*DYMH))/(4.*DP*DMU)
      ELSEIF (DABMU.EQ.0..AND.DYMH.NE.0.) THEN
        AM(I,J)=-(1.-YMU(J)*YMU(J))*YMU(J)*PMH*0.5*(DIFF(I-1,J)+DIFF(I,J)) &
        /(4.*DP*DMU)-DLAMM1*(1.-YMH*YMH)*(XP(I)*YMH*0.5*(DIFF(I,J)+ &
        DIFF(I,J-1))-0.5*(DCY(I,J)+DCY(I,J-1))*(OMC+GAMAB*(YMH*YMH-1.))/(YMH*DYMH))/(4.*DP*DMU)
        AP(I,J)=-(1.-YMU(J)*YMU(J))*(-YMU(J)*PPH*0.5*(DIFF(I+1,J)+DIFF(I,J))) &
        /(4.*DP*DMU)+DLAMM1*(1.-YMH*YMH)*(XP(I)*YMH*0.5*(DIFF(I,J)+ &
        DIFF(I,J-1))-0.5*(DCY(I,J)+DCY(I,J-1))*(OMC+GAMAB*(YMH*YMH-1.))/ &
        (YMH*DYMH))/(4.*DP*DMU)
       ELSEIF (DABMU.NE.0..AND.DYMH.EQ.0.) THEN
          AM(I,J)=-(1.-YMU(J)*YMU(J))*(YMU(J)*PMH*0.5*(DIFF(I-1,J)+ &
                  DIFF(I,J))-0.5*DCYMH*(OMC+GAMAMHB*(YMU(J)*YMU(J)-1.))/(YMU(J)* &
                  DABMU))/(4.*DP*DMU)-DLAMM1*(1.-YMH*YMH)*XP(I)*YMH*0.5*(DIFF(I,J)+ &
                  DIFF(I,J-1))/(4.*DP*DMU)
          AP(I,J)=-(1.-YMU(J)*YMU(J))*(-YMU(J)*PPH*0.5*(DIFF(I+1,J)+ &
                  DIFF(I,J))+0.5*DCYPH*(OMC+GAMAPHB*(YMU(J)*YMU(J)-1.))/(YMU(J)* &
                  DABMU))/(4.*DP*DMU)+DLAMM1*(1.-YMH*YMH)*XP(I)*YMH*0.5*(DIFF(I,J)+ &
                  DIFF(I,J-1))/(4.*DP*DMU)
      ELSE
            AM(I,J)=-(1.-YMU(J)*YMU(J))*YMU(J)*PMH*0.5*(DIFF(I-1,J)+ &
                    DIFF(I,J)) &
                    /(4.*DP*DMU)-DLAMM1*(1.-YMH*YMH)*XP(I)*YMH*0.5*(DIFF(I,J)+ &
                     DIFF(I,J-1))/(4.*DP*DMU)
            AP(I,J)=-(1.-YMU(J)*YMU(J))*(-YMU(J)*PPH*0.5*(DIFF(I+1,J) + DIFF(I,J))) &
                    /(4.*DP*DMU)+DLAMM1*(1.-YMH*YMH)*XP(I)*YMH*0.5*(DIFF(I,J) + DIFF(I,J-1))/(4.*DP*DMU)
        ENDIF
     IF(DABMU.NE.0..AND.DYPH.NE.0.) THEN
          EM(I,J)=(1.-YMU(J)*YMU(J))*(YMU(J)*PMH*0.5*(DIFF(I-1,J)+ &
         DIFF(I,J))-0.5*DCYMH*(OMC+GAMAMHB*(YMU(J)*YMU(J)-1.))/(YMU(J)* &
         DABMU))/(4.*DP*DMU)-DLAMP1*(1.-YPH*YPH)*(-XP(I)*YPH*0.5* &
         (DIFF(I,J)+DIFF(I,J+1))+0.5*(DCY(I,J)+DCY(I,J+1))*(OMC+GAMAB* &
         (YPH*YPH-1.))/(YPH*DYPH))/(4.*DP*DMU)
	  
          EP(I,J)=(1.-YMU(J)*YMU(J))*(-YMU(J)*PPH*0.5*(DIFF(I+1,J)+ &
          DIFF(I,J))+0.5*DCYPH*(OMC+GAMAPHB*(YMU(J)*YMU(J)-1.))/(YMU(J)* &
          DABMU))/(4.*DP*DMU)+DLAMP1*(1.-YPH*YPH)*(-XP(I)*YPH*0.5* &
          (DIFF(I,J)+DIFF(I,J+1))+0.5*(DCY(I,J)+DCY(I,J+1))*(OMC+GAMAB* &
          (YPH*YPH-1.))/(YPH*DYPH))/(4.*DP*DMU)
      ELSEIF(DABMU.EQ.0..AND.DYPH.NE.0.)THEN
          EM(I,J)=(1.-YMU(J)*YMU(J))*(YMU(J)*PMH*0.5*(DIFF(I-1,J)+ &
          DIFF(I,J))) &
          /(4.*DP*DMU)-DLAMP1*(1.-YPH*YPH)*(-XP(I)*YPH*0.5*(DIFF(I,J)+ &
          DIFF(I,J+1))+0.5*(DCY(I,J)+DCY(I,J+1))*(OMC+GAMAB*(YPH*YPH-1.))/ &
          (YPH*DYPH))/(4.*DP*DMU)
          EP(I,J)=(1.-YMU(J)*YMU(J))*(-YMU(J)*PPH*0.5*(DIFF(I+1,J)+ DIFF(I,J))) &
          /(4.*DP*DMU)+DLAMP1*(1.-YPH*YPH)*(-XP(I)*YPH*0.5*(DIFF(I,J)+ &
          DIFF(I,J+1))+0.5*(DCY(I,J)+DCY(I,J+1))*(OMC+GAMAB*(YPH*YPH-1.))/ &
          (YPH*DYPH))/(4.*DP*DMU)
      ELSEIF (DABMU.NE.0..AND.DYPH.EQ.0.) THEN
        EM(I,J)=(1.-YMU(J)*YMU(J))*(YMU(J)*PMH*0.5*(DIFF(I-1,J)+ &
            DIFF(I,J))-0.5*DCYMH*(OMC+GAMAMHB*(YMU(J)*YMU(J)-1.))/(YMU(J)* &
            DABMU))/(4.*DP*DMU)-DLAMP1*(1.-YPH*YPH)*(-XP(I)*YPH*0.5* &
            (DIFF(I,J)+DIFF(I,J+1)))/(4.*DP*DMU) 
          EP(I,J)=(1.-YMU(J)*YMU(J))*(-YMU(J)*PPH*0.5*(DIFF(I+1,J)+ &
          DIFF(I,J))+0.5*DCYPH*(OMC+GAMAPHB*(YMU(J)*YMU(J)-1.))/(YMU(J)* &
          DABMU))/(4.*DP*DMU)+DLAMP1*(1.-YPH*YPH)*(-XP(I)*YPH*0.5* &
          (DIFF(I,J)+DIFF(I,J+1)))/(4.*DP*DMU)
      ELSE
          EM(I,J)=(1.-YMU(J)*YMU(J))*(YMU(J)*PMH*0.5*(DIFF(I-1,J)+DIFF(I,J))) &
          /(4.*DP*DMU)-DLAMP1*(1.-YPH*YPH)*(-XP(I)*YPH*0.5*(DIFF(I,J)+ DIFF(I,J+1)))/(4.*DP*DMU)
          EP(I,J)=(1.-YMU(J)*YMU(J))*(-YMU(J)*PPH*0.5*(DIFF(I+1,J)+DIFF(I,J))) &
          /(4.*DP*DMU)+DLAMP1*(1.-YPH*YPH)*(-XP(I)*YPH*0.5*(DIFF(I,J)+ &
          DIFF(I,J+1)))/(4.*DP*DMU)
     ENDIF
     !Cmode 5 diagonales
     IF(MODE.EQ.5) THEN
     !C           compute RF dependent (and time-independent if no FW) part
     !C            of matrix of coefficient Q /2/A.13
                Q(I,J)=T(I,J)*XP(I)*XP(I)/DT &
                  -AM(I,J)*T(I-1,J-1)-AP(I,J)*T(I+1,J-1) &
                   -EP(I,J)*T(I+1,J+1)-EM(I,J)*T(i-1,J+1)
     ENDIF
     !Cmode 9 diagonales

     IF (MODE.EQ.9) THEN
     !C        compute matrix of coefficient Q
                Q(I,J)=T(I,J)*XP(I)*XP(I)/DT
     ENDIF
     !Cwhatever the choice, 5 or 9 points, the coefficients A,B,C,D,E,Q and
     !Cparts thereof have been computed

  4   CONTINUE
      DO 444 I=1,N1
      Q(I,1)=Q(I,3)
	  Q(I,N2)=Q(I,N2-2)
      A(I,1)=0.
      E(I,N2)=0.
      B(I,1)=B(I,3)
      C(I,1)=C(I,3)
      C(I,N2)=C(I,N2-2)
      D(I,1)=D(I,3)
      D(I,N2)=D(I,N2-2)
      A(I,N2)=A(I,N2-2)
      E(I,1)=E(I,3)
      !C9 diagonales
      AM(I,1) = AM(I,3)
      AP(I,1) = AP(I,3)
      EM(I,1) = EM(I,3)
      EP(I,1) = EP(I,3)
      AM(I,N2) = AM(I,N2-2)
      AP(I,N2) = AP(I,N2-2)
      EM(I,N2) = EM(I,N2-2)
      EP(I,N2) = EP(I,N2-2)
 444  CONTINUE
      DO 445 J=1,N2
!Cla ligne suivante a ete ajoute pour la condition aux limites au centre
      C(1,J)=1.
      D(1,J)=-F(1)/F(2)
 445  CONTINUE
!C    IT=1
!C    FOR LHCD ALONE (NO FAST WAVE), IF DIFF(I,J) AND EPS ARE CONSTANT
!C    START LOOP HERE
 30   CONTINUE
      DINT=0.
      DINT2=0.
      DO 151 I=1,N1
      DINT=0.
      DO 150 J=2,N2M1
      DINT=DINT+T(I,J)*YMU(J)
      DINT2=DINT2+T(I,J)*(3.*YMU(J)*YMU(J)-1.)
 150  CONTINUE
      T2(I)=DINT2*DMU*5./4.
      T1(I)=DINT*DMU*3./2.
 151  CONTINUE
      AB1(1)=0.
      AB2(1)=0.
      AB3(1)=0.
      AB4(1)=0.
      AB5(1)=0.
      AB6(1)=0.
      AB7(1)=0.
      AB8(1)=0.
      AB9(1)=0.
      AB10(1)=0.
      AB1(2)=0.
      AB2(2)=0.
      AB3(2)=0.
      AB4(2)=0.
      AB5(2)=0.
      AB6(2)=0.
      AB7(2)=0.
      AB8(2)=0.
      AB9(2)=0.
      AB10(2)=0.
      DO 152 I=3,N1
      XP2=XP(I)*XP(I)
      GAMA=DSQRT(1.+XP2*BTH*BTH)
      V=XP(I)/GAMA
      GAMA1=DSQRT(1.+XP(I-1)*XP(I-1)*BTH*BTH)
      V1=XP(I-1)/GAMA1
      AB1(I)=AB1(I-1)+0.5*DP*T1(I)*(XP(I)**3)*4.*pi+0.5*DP*T1(I-1)*(XP(I-1)**3)*4.*pi
	  
      AB2(I)=AB2(I-1)+0.5*DP*T1(I)*(XP(I)**5)*4.*pi/GAMA+0.5*DP*T1(I-1)*(XP(I-1)**5)*4.*pi/GAMA1
	  
      AB3(I)=AB3(I-1)+0.5*DP*T1(I)*XP(I)*4.*pi*(GAMA-SIG(I)/(BTH*XP(I)))/GAMA &
      +0.5*DP*T1(I-1)*XP(I-1)*4.*pi*(GAMA1-SIG(I-1)/(BTH*XP(I-1)))/GAMA1
      AB4(I)=AB4(I-1)+0.5*DP*T1(I)*(XP(I)**3)*4.*pi/GAMA+0.5*DP*T1(I-1)*(XP(I-1)**3)*4.*pi/GAMA1
	  
      AB5(I)=AB5(I-1)+0.5*DP*T1(I)*(XP(I)**3)*4.*pi*(GAMA-SIG(I)/(BTH*XP(I))-2.*GAMA*BTH*BTH*XP(I)*XP(I)/3.) &
      /(GAMA*BTH*BTH*XP(I)*XP(I))+0.5*DP*T1(I-1)*(XP(I-1)**3)*4.*pi*(GAMA1-SIG(I-1) &
      /(BTH*XP(I-1))-2.*GAMA1*BTH*BTH*XP(I-1)*XP(I-1)/3.)/(GAMA1*BTH*BTH*XP(I-1)*XP(I-1))
	  
      AB6(I)=AB6(I-1)+0.5*DP*T1(I)*(XP(I)**3)*4.*pi &
      *(-3.*GAMA*SIG(I)/(BTH*BTH*XP(I)*XP(I))+3./(BTH*XP(I)) &
      +BTH*XP(I))/(6.*GAMA*BTH*XP(I)) &
      +0.5*DP*T1(I-1)*(XP(I-1)**3)*4.*pi &
      *(-3.*GAMA1*SIG(I-1)/(BTH*BTH*XP(I-1)*XP(I-1))+3./(BTH*XP(I-1))+BTH*XP(I-1))/(6.*GAMA1*BTH*XP(I-1))
	  
      AB8(I)=AB8(I-1)+.5*DP*T1(I)*4.*pi*(XP(I)**3)*(-3.*GAMA+SIG(I) &
      *(3./(BTH*XP(I))+2.*BTH*XP(I)))/(2.*GAMA*BTH*BTH*XP(I)*XP(I)) &
      +0.5*DP*T1(I-1)*4.*pi*(XP(I-1)**3)*(-3.*GAMA1+SIG(I-1)*(3./  &
      (BTH*XP(I-1))+2.*BTH*XP(I-1)))/(2.*GAMA1*BTH*BTH*XP(I-1)*XP(I-1))
	  
      AB7(I)=AB7(I-1)+.5*DP*T1(I)*4.*pi*(XP(I)**3)*(-3.*GAMA*SIG(I) &
      /(BTH*BTH*XP(I)*XP(I))+3./(BTH*XP(I))+BTH*XP(I)-2.*(BTH*XP(I))**3/5.)/(2.*GAMA*BTH*BTH*BTH*XP(I)) &
      +.5*DP*T1(I-1)*4.*pi*(XP(I-1)**3)*(-3.*GAMA1*SIG(I-1) &
      /(BTH*BTH*XP(I-1)*XP(I-1))+3./(BTH*XP(I-1))+BTH*XP(I-1)-2.*(BTH*XP(I-1))**3/5.)/(2.*GAMA1*BTH*BTH*BTH*XP(I-1))
	  
      AB9(I)=AB9(I-1)+0.5*DP*T1(I)*XP(I)*4.*pi*(GAMA*SIG(I)/(BTH*XP(I))-1.)/GAMA &
      +0.5*DP*T1(I-1)*XP(I-1)*4.*pi*(GAMA1*SIG(I-1)/(BTH*XP(I-1))-1.)/GAMA1
	  
      AB10(I)=AB10(I-1)+0.5*DP*T1(I)*(XP(I)**3)*4.*pi*(GAMA*SIG(I) &
      *(15./(BTH*XP(I))**2+6.)-15./(BTH*XP(I))-11.*BTH*XP(I))/ &
      (12.*GAMA*BTH*BTH*BTH*XP(I))+0.5*DP*T1(I-1)*(XP(I-1)**3)*4.*pi*(GAMA1*SIG(I-1) &
      *(15./(BTH*XP(I-1))**2+6.)-15./(BTH*XP(I-1))-11.*BTH*XP(I-1))/ &
      (12.*GAMA1*BTH*BTH*BTH*XP(I-1))
 152  CONTINUE
      AB11(N1)=0.
      AB12(N1)=0.
      AB13(N1)=0.
      AB14(N1)=0.
      DO 155 I=1,N1M1
      K=N1-I
      XP2=XP(K)*XP(K)
      GAMA=DSQRT(1.+XP2*BTH*BTH)
      GAMA1=DSQRT(1.+XP(K+1)*XP(K+1)*BTH*BTH)
      AB11(K)=AB11(K+1)+0.5*DP*T1(K)*4.*pi/GAMA+0.5*DP*T1(K+1)*4.*pi/GAMA1
	  
      AB12(K)=AB12(K+1)+0.5*DP*T1(K)*4.*pi+0.5*DP*T1(K+1)*4.*pi
	  
      AB13(K)=AB13(K+1)+0.5*DP*T1(K)*(XP(K)**2)*4.*pi+0.5*DP*T1(K+1)*(XP(K+1)**2)*4.*pi
	  
      AB14(K)=AB14(K+1)+0.5*DP*T1(K)*(XP(K)**2)*4.*pi/GAMA+.5*DP*T1(K+1)*(XP(K+1)**2)*4.*pi/GAMA1
 155  CONTINUE
      DO 502 I=2,N1M1
      DO 502 J=2,N2M1
!Cmode 5 diagonales
            IF(MODE.EQ.5) THEN
                Q(I,J)=T(I,J)*XP(I)*XP(I)/DT &
                   -AM(I,J)*T(I-1,J-1)-AP(I,J)*T(I+1,J-1) &
                    -EP(I,J)*T(I+1,J+1)-EM(I,J)*T(i-1,J+1)
            ENDIF
!Cmode 9 diagonales
            IF (MODE.EQ.9) THEN
                Q(I,J)=T(I,J)*XP(I)*XP(I)/DT
            ENDIF
!C

 502  CONTINUE
      DO 503 I=2,N1M1
      GAMA=DSQRT(1.D0+XP(I)*XP(I)*BTH*BTH)
      V=XP(I)/GAMA
	  CNT=-2.*AB1(I)*GAMA/3.+AB2(I)*GAMA/5.+AB3(I)+AB4(I)/3. &
      +(-AB5(I)+AB7(I)+AB8(I)-AB10(I))*GAMA
      CNT=CNT+AB6(I)*(GAMA*XP(I)*XP(I)-5.)+AB9(I)*XP(I)*XP(I)
      CNT=CNT/(GAMA*XP(I)*XP(I))
      CNT1=+AB11(I)*(1./3.+(GAMA-SIG(I)/(BTH*XP(I)))/(XP(I)*XP(I)) &
      -5.*(-3.*GAMA*SIG(I) &
      /(BTH*BTH*XP(I)*XP(I))+3./(BTH*XP(I))+BTH*XP(I))/(6.*BTH*XP(I))) 
      CNT1=CNT1+AB12(I)*(-2.*GAMA/3.+XP(I)*XP(I)/5.-(GAMA-SIG(I) &
      /(BTH*XP(I))-2.*GAMA*BTH*BTH*XP(I)*XP(I)/3.)/(BTH*BTH*XP(I) &
      *XP(I)))
      CNT1=CNT1+AB12(I)*(-3.*GAMA*SIG(I)/(BTH*BTH*XP(I)*XP(I))+3. &
      /(BTH*XP(I))+BTH*XP(I)-2.*(BTH*XP(I))**3/5.)/(2.*BTH*BTH*BTH*XP(I))
      CNT1=CNT1+AB12(I)*(-3.*GAMA+SIG(I)*(3./(BTH*XP(I))+2.*BTH*XP(I)))/(2.*BTH*BTH*XP(I)*XP(I))
      CNT1=CNT1-AB12(I)*(GAMA*SIG(I)*(15./(BTH*BTH*XP(I)*XP(I))+6.)-15./(BTH*XP(I))-11.*BTH*XP(I))/(12.*BTH*XP(I)*BTH*BTH)
      CNT1=CNT1+AB13(I)*(-3.*GAMA*SIG(I)/(BTH*BTH*XP(I)*XP(I))+3./(BTH*XP(I))+BTH*XP(I))/(6.*BTH*XP(I))
      CNT1=CNT1+AB14(I)*(GAMA*SIG(I)/(BTH*XP(I))-1.)/(XP(I)*XP(I))
      CNT1=CNT1*XP(I)/GAMA
      CNT=(CNT*F(I)+CNT1*F(I))+4.*pi*F(I)*T1(I)/GAMA
      CNT=CNT*XP(I)*XP(I)
      DO 503 J=2,N2M1
      Q(I,J)=Q(I,J)+EJ0(J)*CNT*YMU(J)*TRUNC
 503  CONTINUE
      DO 513 I=1,N1
      Q(I,1)=Q(I,3)
      Q(I,N2)=Q(I,N2-2)
!C9 diagonales
      AM(I,1) = AM(I,3)
      AP(I,1) = AP(I,3)
      EM(I,1) = EM(I,3)
      EP(I,1) = EP(I,3)
      AM(I,N2) = AM(I,N2-2)
      AP(I,N2) = AP(I,N2-2)
      EM(I,N2) = EM(I,N2-2)
      EP(I,N2) = EP(I,N2-2)
 513  CONTINUE
      IF(TRAPO.EQ.0.0D0) GO TO 1701
      N2H=N2/2
      DO I=1,N1
      DO J=1,N2H
          ABSM=DABS(YMU(J))
          IF(ABSM.LT.TRAP) A(I,J)=A(I,N2+1-J)
          IF(ABSM.LT.TRAP) B(I,J)=B(I,N2+1-J)
          IF(ABSM.LT.TRAP) C(I,J)=C(I,N2+1-J)
          IF(ABSM.LT.TRAP) D(I,J)=D(I,N2+1-J)
          IF(ABSM.LT.TRAP) E(I,J)=E(I,N2+1-J)
          IF(ABSM.LT.TRAP) Q(I,J)=Q(I,N2+1-J)
          IF(ABSM.LT.TRAP) AM(I,J)=AM(I,N2+1-J)
          IF(ABSM.LT.TRAP) AP(I,J)=AP(I,N2+1-J)
          IF(ABSM.LT.TRAP) EM(I,J)=EM(I,N2+1-J)
          IF(ABSM.LT.TRAP) EP(I,J)=EP(I,N2+1-J)
      ENDDO
      ENDDO
 1701 CONTINUE
      DO 80 J=2,N2M1
      DO 60 I=2,N1-1
      IF (MODE.EQ.5) THEN
          !Cmode 5 diagonales
          IF(C(I,J).NE.0.D0) THEN
              R(I,J)=Q(I,J)-A(I,J)*T(I,J-1)-B(I,J)*T(I-1,J)-C(I,J)*T(I,J)- &
              D(I,J)*T(I+1,J)-E(I,J)*T(I,J+1)-COLLF*XP(I)*XP(I)*T(I,J)
          ELSE
            R(I,J)=Q(I,J)-T(I,J)*XP(I)*XP(I)/DT
          END IF
      END IF
      IF (MODE.EQ.9) THEN
          !Cmode 9 diagonales
          IF(C(I,J).NE.0.D0) THEN
              R(I,J)=Q(I,J)-A(I,J)*T(I,J-1)-B(I,J)*T(I-1,J)-C(I,J)*T(I,J)- &
              D(I,J)*T(I+1,J)-E(I,J)*T(I,J+1) &
              -COLLF*XP(I)*XP(I)*T(I,J) &
              -AM(I,J)*T(I-1,J-1)-AP(I,J)*T(I+1,J-1)-EM(I,J)*T(I-1,J+1) &
              -EP(I,J)*T(I+1,J+1)
          ELSE
              R(I,J)=Q(I,J)-T(I,J)*XP(I)*XP(I)/DT
          END IF
      END IF
  60  CONTINUE
  80  CONTINUE
      DO 81 I=1,N1M1
          IF(C(I,N2).NE.0.0D0) THEN
              R(I,1)=R(I,3)
              R(I,N2)=R(I,N2-2)
          ELSE
              R(I,N2)=Q(I,N2)-T(I,N2)*XP(I)*XP(I)/DT
          END IF
  81  CONTINUE
  55  FORMAT(2I10,2F16.6)
      DO 500 J=1,N2
          Q(N1,J)=T(N1,J)*XP(N1)*XP(N1)/DT
          Q(1,J)=0.
          R(1,J)=0.
          R(N1,J)=0.
 500  CONTINUE
      IFAIL=0
!C    CALL D03UAF(N1,N2,N1M,A,B,C,D,E,APARAM,IT,R,WRKSP1,WRKSP2,IFAIL)
	  IF (MODE.EQ.5) THEN
	      ALPHAMAX = 1.- 2.*APARAM/(N1M1*N1M1+N2M1*N2M1)
          call sip5d(N1,N2,A,B,C,D,E,R,IT,ALPHAMAX,WRKSP1,WRKSP2)
      END IF
      IF (MODE.EQ.9) THEN
	      ALPHAMAX = 1.- 2.*APARAM/(N1M1*N1M1+N2M1*N2M1)
          call sip9d(N1,N2,AM,A,AP,B,C,D,EM,E,EP,R,IT,ALPHAMAX,WRKSP1,WRKSP2,WRKSP3,WRKSP4)
      ENDIF
      MPRINT=MPRINT+1
      IT=IT+1
!C    IF YOU WANT DT TO REMAIN CONSTANT AT EACH ITERATION, COMMENT THE
!C    FOLLOWING CARD
!C    DT=DTT+(IT-1)*DTT/NITS
      DO 6 J=1,N2
      DO 6 I=1,N1
        T(I,J)=T(I,J)+R(I,J)
   6  CONTINUE
      DO 501 I=1,N1
        T(I,1)=T(I,3)
 501    T(I,N2)=T(I,N2-2)
      IF(TRAPO.EQ.0.0D0) GO TO 1700
      N2H=N2/2
      DO I=1,N1
      DO J=1,N2H
      ABSM=DABS(YMU(J))
      CNST=1.D0-RATIOB*(1.D0-YMU(J)*YMU(J))
      IF(CNST.GT.0.0.AND.CNST.LE.1.) THEN
      CNST=DSQRT(CNST)
      ELSE
      CNST=0.
      ENDIF
      IF(ABSM.LT.TRAP) T(I,N2+1-J)=(T(I,N2+1-J)+T(I,J))/2.
      IF(ABSM.LT.TRAP) T(I,J)=T(I,N2+1-J)
!C     IF(CNST.EQ.0.0) T(I,J)=T(I,N2+1-J)
      ENDDO
      ENDDO
	GO TO 1700
      DO J=2,N2H
      ABSM=DABS(YMU(J))
      ABSM1=DABS(YMU(J+1))
      CNST=1.D0-RATIOB*(1.D0-YMU(J)*YMU(J))
      IF(CNST.GT.0.0.AND.CNST.LE.1.) THEN
      CNST=DSQRT(CNST)
      ELSE
      CNST=0.
      ENDIF
      CNST1=1.D0-RATIOB*(1.D0-YMU(J+1)*YMU(J+1))
      IF(CNST1.GT.0.0.AND.CNST1.LE.1.) THEN
      CNST1=DSQRT(CNST1)
      ELSE
      CNST1=0.
      ENDIF
      DO I=1,N1M1
      IF(ABSM1.LE.TRAP.AND.ABSM.GT.TRAP) T(I,J+1)=(T(I,J+1)+T(I,J))/2.
      IF(ABSM1.LE.TRAP.AND.ABSM.GT.TRAP) T(I,J)=(T(I,J+1)+T(I,J))/2.
      IF(ABSM1.LE.TRAP.AND.ABSM.GT.TRAP) T(I,J-1)=(T(I,J-1)+T(I,J))/2.
      ENDDO
      ENDDO
      
!C    do i=1,n1
!C         do j=1,n2
!C             if(T(i,j).gt.T(1,1)) then
!C                 T(i,j)=1.d-200
!C                 WRITE(6,*) 'i=',i,' j=',j
!C             endif
!C         enddo
!C     enddo
 1700 CONTINUE
      IF(IT.EQ.NITS) GO TO 101
      IF(MPRINT.LT.NPRINT) GO TO 100
      MPRINT=0.
  51  FORMAT(5X,4F16.8,I10)
  101 DENS=0.0
      CUR=0.
      CURR=0.
      CURRR=0.
      DO 5 J=2,N2M1
      DO 5 I=2,N1M1
      ABSM=DABS(YMU(J))
      XPH=XP(I)-DP/2.
      GAMAH=DSQRT(1.D0+XPH*XPH*BTH*BTH)
      GAMA=DSQRT(1.D0+XP(I)*XP(I)*BTH*BTH)
      CNSTT=1.D0-RATIOB*(1.D0-YMU(J)*YMU(J))
      IF(CNSTT.GT.0.0.AND.CNSTT.LE.1.) THEN
      CNSTT=DSQRT(CNSTT)
      ELSE
      CNSTT=0.
      ENDIF
!C    IF(YMU(J).LT.0.0.AND.ABSM.GT.TRAP) CNSTT=-CNSTT
      IF(YMU(J).LT.0.0) CNSTT=-CNSTT
      PPAR=XP(I)*CNSTT
      VPAR=PPAR/GAMA
	IF(TRAPO.EQ.0.) THEN
      CUR=CUR+T(I,J)*PPAR*XP(I)*XP(I)/GAMA
	ELSE
!C    CUR=CUR+2.*T(I,J)*PPAR*XP(I)*XP(I)/GAMA
!C    1+2.*(T(I,J)+T(I-1,J))*PPARH*XPH*XPH/GAMAH
!C!C     IF(CNSTT.LT.0.0) CURRR=CURRR+T(I,J)*PPAR*XP(I)*XP(I)/GAMA
	IF(ABSM.GT.TRAP.AND.YMU(J).LT.0.) CURRR=CURRR+T(I,J)*PPAR*XP(I)*XP(I)/GAMA
!C!C     IF(CNSTT.GT.0.0) CUR=CUR+T(I,J)*PPAR*XP(I)*XP(I)/GAMA
	IF(ABSM.GT.TRAP.AND.YMU(J).GT.0.) CUR=CUR+T(I,J)*PPAR*XP(I)*XP(I)/GAMA
!C!C     IF(CNSTT.EQ.0.0) CURR=CURR+T(I,J)*PPAR*XP(I)*XP(I)/GAMA
	IF(ABSM.LE.TRAP) CURR=CURR+T(I,J)*PPAR*XP(I)*XP(I)/GAMA
	ENDIF
      DENS=DENS+2.*T(I,J)*XP(I)*XP(I)+2.*(T(I,J)+T(I-1,J))*XPH*XPH
   5  CONTINUE
!C    CALCULATE DERIVATIVE OF THE DISTRIBUTION FUNCTION USING CUBIC SPLINE
      DO 506 J=1,N2
      DE(1)=0.
      DF(1)=0.
      DO 507 I=2,N1M1
      DFI=3.*(T(I+1,J)-T(I-1,J))/DP
      DE(I)=-1./(DE(I-1)+4.)
 507  DF(I)=(DF(I-1)-DFI)*DE(I)
      DFDP(N1,J)=0.0
      DO 508 I=1,N1M1
      K=N1-I
 508  DFDP(K,J)=DE(K)*DFDP(K+1,J)+DF(K)
 506  CONTINUE
      DO 526 I=1,N1
      DE(1)=0.
      DF(1)=0.0
      DO 527 J=2,N2M1
      DFI=3.*(T(I,J+1)-T(I,J-1))/DMU
      DE(J)=-1./(DE(J-1)+4.)
 527  DF(J)=(DF(J-1)-DFI)*DE(J)
      DFDM(I,N2)=0.0
      DO 528 J=1,N2M1
      K=N2-J
 528  DFDM(I,K)=DE(K)*DFDM(I,K+1)+DF(K)
 526  CONTINUE
      PABS1=0.
	PABSR=0.
	PABSRR=0.
	PABSRRR=0.
      DO 505 J=2,N2M1
      DO 505 I=2,N1M1
      ABSM=DABS(YMU(J))
      GAMA=DSQRT(1.D0+XP(I)*XP(I)*BTH*BTH)
      GAMAB=GAMA*RATIOB
      CNSTT=1.D0-RATIOB*(1.D0-YMU(J)*YMU(J))
      IF(CNSTT.GT.0.0.AND.CNSTT.LE.1.) THEN
      CNSTT=DSQRT(CNSTT)
      ELSE
      CNSTT=0.
      ENDIF
      IF(YMU(J).LT.0.0) CNSTT=-CNSTT
      PPAR=XP(I)*CNSTT
      PPERP=XP(I)*XP(I)-PPAR*PPAR
      VPAR=PPAR/GAMA
      IF(PPAR.EQ.0.) GO TO 3000
	IF(TRAPO.EQ.0.) THEN
      PABSR=PABSR+(XP(I)**3)*DIFF(I,J)*YMU(J)*(YMU(J)*DFDP(I,J)+ &
      (1.-YMU(J)*YMU(J))*DFDM(I,J)/XP(I))/GAMA &
      +(XP(I)**3)*DCY(I,J)*(1.-YMU(J)*YMU(J))*(GAMAB*DFDP(I,J)-(OMC- &
      GAMAB*(1.-YMU(J)*YMU(J)))*DFDM(I,J)/(XP(I)*YMU(J)))/(GAMA* &
      DABS(PPAR))
    ELSE
!C!C     IF(CNSTT.LT.0.0)
	IF(ABSM.GE.TRAP.AND.YMU(J).LT.0.) &
      PABSRRR=PABSRRR+(XP(I)**3)*DIFF(I,J)*YMU(J)*(YMU(J)*DFDP(I,J)+ &
      (1.-YMU(J)*YMU(J))*DFDM(I,J)/XP(I))*TRPP(J)/GAMA+TRPP(J)* &
      (XP(I)**3)*DCY(I,J)*(1.-YMU(J)*YMU(J))*(GAMAB*DFDP(I,J)-(OMC-GAMAB &
      *(1.-YMU(J)*YMU(J)))*DFDM(I,J)/(XP(I)*YMU(J)))/(GAMA*DABS(PPAR))
!C!C     IF(CNSTT.GT.0.0)
	IF(ABSM.GE.TRAP.AND.YMU(J).GT.0.) &
      PABSR=PABSR+(XP(I)**3)*DIFF(I,J)*YMU(J)*(YMU(J)*DFDP(I,J)+ &
      (1.-YMU(J)*YMU(J))*DFDM(I,J)/XP(I))*TRPP(J)/GAMA+TRPP(J)* &
      (XP(I)**3)*DCY(I,J)*(1.-YMU(J)*YMU(J))*(GAMAB*DFDP(I,J)-(OMC-GAMAB &
      *(1.-YMU(J)*YMU(J)))*DFDM(I,J)/(XP(I)*YMU(J)))/(GAMA*DABS(PPAR)) 
!C!C  IF(CNSTT.EQ.0.0)
	IF(ABSM.LT.TRAP) &
      PABSRR=PABSRR+(XP(I)**3)*DIFF(I,J)*YMU(J)*(YMU(J)*DFDP(I,J)+ &
      (1.-YMU(J)*YMU(J))*DFDM(I,J)/XP(I))*TRPP(J)/GAMA+TRPP(J)* &
      (XP(I)**3)*DCY(I,J)*(1.-YMU(J)*YMU(J))*(GAMAB*DFDP(I,J)-(OMC-GAMAB &
      *(1.-YMU(J)*YMU(J)))*DFDM(I,J)/(XP(I)*YMU(J)))/(GAMA*DABS(PPAR))
	ENDIF
 3000 CONTINUE
 505  CONTINUE
      if( IT.eq.NITS) then
!C    THE FOLLOWING LOOP CALCULATE AND WRITE THE RUNAWAY RATE AT DIFFERENT
!C    RADII RAD
      DO 517 I=100,N1,25
      EPS=EEPS*EJ0(J)
      RUN1=0.
      RUN=0.
      RAD=XP(I)
      GAMA=DSQRT(1.D0+XP(I)*XP(I)*BTH*BTH)
      GAMAPH=GAMA
      V=XP(I)/GAMA
      VPH=V
      PPH=XP(I)
      XSIG=PPH*BTH
      SIGPH=S11ABF(XSIG,IIFAIL)
      PPH2=PPH*PPH
      SPH=V/DSQRT(2.D0)
      EXPPH=DEXP(-PPH2/(1.+GAMAPH))/(2.*pi)**1.5
      FPHI=AP1(I)/(VPH**2)+AP2(I)/PPH2+AP3(I)*(1.-SIGPH/(BTH*PPH*GAMAPH))/VPH
      APHI=FPHI/VPH+Z/(VPH**3*CC*CC)
      FPHI=FPHI+Z*TE/(TI*VPH*VPH*CC*CC)
      DO 516 J=2,N2M1
      RUN=RUN+(EPS*YMU(J)*T(I,J) - APHI*(T(I+1,J)-T(I-1,J))/(2.*DP)-FPHI*T(I,J))*DMU
      RUN1=RUN1+(EPS*YMU(J)*T(I,J) - APHI*DFDP(I,J)-FPHI*T(I,J))*DMU
  516 CONTINUE
      RUN=-RUN*2.*pi*RAD*RAD
      RUN1=-RUN1*2.*pi*RAD*RAD
      WRITE(6,518) RAD,gama
      WRITE(6,518) RUN,RUN1
      rdrd1=rad
      gmgm1=gama
  517 CONTINUE
      endif
 518  FORMAT(5X,F15.6,2x,G15.6,2x,G15.6 )
      DENS=DENS*2.*pi*DP*DMU/6.
!C     CUR=CUR*2.*pi*DP*DMU/6.
      CUR=CUR*2.*pi*DP*DMU
      CURR=CURR*2.*pi*DP*DMU
      CURRR=CURRR*2.*pi*DP*DMU
      TCURR=CUR+CURR+CURRR
      PABSR=-PABSR*2.*pi*DP*DMU
      PABSRR=-PABSRR*2.*pi*DP*DMU
      PABSRRR=-PABSRRR*2.*pi*DP*DMU
      PABS1=PABSR+PABSRR+PABSRRR
      IF(PABS1.EQ.0.) GO TO 5181
      RJP1=TCURR/PABS1
 5181 CONTINUE
!C    TO RENORMALIZE THE SOLUTION TO KEEP DENSITY STRICTLY 1
      DO 300 J=1,N2
      DO 300 I=1,N1
  300 T(I,J)=T(I,J)/DENS
!C    RENORMALIZED VALUES
      CUR=CUR/DENS 
      CURR=CURR/DENS
      CURRR=CURRR/DENS
      TCURR=TCURR/DENS
      PABSR=PABSR/DENS
      PABSRR=PABSRR/DENS
      PABSRRR=PABSRRR/DENS
      PABS1=PABS1/DENS
      DENS=DENS/DENS
  100 CONTINUE
!C     WRITE(6,*) IT,DENS,TCURR,PABS1,RJP1
!C     WRITE(6,*) CUR,CURR,CURRR,eeps
  50  FORMAT(I4,2X,4F15.6)
      if(it.eq.1) ctime(it)=dt
      if(it.gt.1) ctime(it)=ctime(it-1)+dt
      curtim(it)=tcurr
      pabsti(it)=pabs1
      IF(IT/100.ne.itpr) GO TO 928
      WRITE(6,*) IT,DENS,TCURR,PABS1,RJP1
      WRITE(6,*) CUR,CURR,CURRR,eeps
      itpr=itpr+1
 928  continue
      if (IT.eq.NITS/5) then
      OPEN(UNIT=11,FILE='distr02.dat')
      write(11,975) N1
      write(11,975) N2
      write(11,975) nits
      write(11,974) (XP(i),i=1,n1)
!C           write(11,*)  '  '
      write(11,974) (YMU(j),j=2,n2-1)
!C           write(11,*)  '  '
      write(11,974) ((T(I,J),I=1,N1),J=2,N2M1)
!C           write(11,*)  '  '
      write(11,974) (ctime(i),i=1,nits)
!C     write(11,*)  '  '
      write(11,974) (curtim(i),i=1,nits)
!C           write(11,*)  '  '
      write(11,974) (pabsti(i),i=1,nits)
      
      CLOSE(11)
      endif      
      if (IT.eq.2*NITS/5) then
      OPEN(UNIT=11,FILE='distr04.dat')
      write(11,975) N1
      write(11,975) N2
      write(11,975) nits
      write(11,974) (XP(i),i=1,n1)
!C           write(11,*)  '  '
      write(11,974) (YMU(j),j=2,n2-1)
!C           write(11,*)  '  '
      write(11,974) ((T(I,J),I=1,N1),J=2,N2M1)
!C           write(11,*)  '  '
      write(11,974) (ctime(i),i=1,nits)
!C     write(11,*)  '  '
      write(11,974) (curtim(i),i=1,nits)
!C           write(11,*)  '  '
      write(11,974) (pabsti(i),i=1,nits)
      
      CLOSE(11)
      endif      
      if (IT.eq.3*NITS/5) then
      OPEN(UNIT=11,FILE='distr06.dat')
      write(11,975) N1
      write(11,975) N2
      write(11,975) nits
      write(11,974) (XP(i),i=1,n1)
!C           write(11,*)  '  '
      write(11,974) (YMU(j),j=2,n2-1)
!C           write(11,*)  '  '
      write(11,974) ((T(I,J),I=1,N1),J=2,N2M1)
!C           write(11,*)  '  '
      write(11,974) (ctime(i),i=1,nits)
!C     write(11,*)  '  '
      write(11,974) (curtim(i),i=1,nits)
!C           write(11,*)  '  '
      write(11,974) (pabsti(i),i=1,nits)
      
      CLOSE(11)
      endif      
      if (IT.eq.4*NITS/5) then
      OPEN(UNIT=11,FILE='distr08.dat')
      write(11,975) N1
      write(11,975) N2
      write(11,975) nits
      write(11,974) (XP(i),i=1,n1)
!C           write(11,*)  '  '
      write(11,974) (YMU(j),j=2,n2-1)
!C           write(11,*)  '  '
      write(11,974) ((T(I,J),I=1,N1),J=2,N2M1)
!C           write(11,*)  '  '
      write(11,974) (ctime(i),i=1,nits)
!C     write(11,*)  '  '
      write(11,974) (curtim(i),i=1,nits)
!C           write(11,*)  '  '
      write(11,974) (pabsti(i),i=1,nits)
      
      CLOSE(11)
      endif      
      
      
      IF(IT.LT.NITS) GO TO 30

  56  FORMAT(5X,6F14.8)
!C    WRITE THE DISTRIBUTION FUNCTION IN A FILE
!C     OPEN(UNIT=11,FILE='distr.dat',FORM='UNFORMATTED')
!Cc	  OPEN(UNIT=11,FILE='distr.dat',STATUS='UNKNOWN')
!C     Write(11) N1,N2,nits
!C     write(11) XP,YMU
!C     WRITE(11) ((T(I,J),I=1,N1),J=2,N2M1)
!C     WRITE(11) ctime,curtim,pabsti
!C     CLOSE(11)
!C    WRITE THE DISTRIBUTION FUNCTION IN A FILE
      do i=1,n1
      do j=1,n2
      if (abs(T(i,j).le. 1.d-100)) T(i,j)=1.d-100
      enddo
      enddo
  975 format(I10)
  974 Format(6(e20.13,1X))
      OPEN(UNIT=11,FILE='distr.dat')
      write(11,975) N1
      write(11,975) N2
      write(11,975) nits
      write(11,974) (XP(i),i=1,n1)
!C           write(11,*)  '  '
      write(11,974) (YMU(j),j=2,n2-1)
!C           write(11,*)  '  '
      write(11,974) ((T(I,J),I=1,N1),J=2,N2M1)
!C           write(11,*)  '  '
      write(11,974) (ctime(i),i=1,nits)
!C     write(11,*)  '  '
      write(11,974) (curtim(i),i=1,nits)
!C           write(11,*)  '  '
      write(11,974) (pabsti(i),i=1,nits)
      
      CLOSE(11)



      OPEN(UNIT=11,FILE='res.txt',STATUS='UNKNOWN')
!C	  OPEN(UNIT=11,FILE='distr.dat',STATUS='UNKNOWN')
      Write(11,*) 'N1=',N1,' N2=',n2
      Write(11,*) 'nits=',nits 
      Write(11,*) 'PMAX=',pmax,' PMIN=',pmin
      Write(11,*) 'TRAPO=',TRAPO,' elfield=',elfield
      write(11,*) 'RMIN=',rmin,' RMAJ=',RMAJ
      write(11,*) 'Te=',Te,'Ti=',ti
      write(11,*) 'pdens=',pdens
      WRITE(11,*) 'Ddif=',ddif
      write(11,*) 'V1=',v1i,' V2=',v2
      WRITE(11,*) 'TCURR=',tcurr,' pabs1=',pabs1
      write(11,*) 'cur=',cur,' curr=',curr,' currr=',currr
      write(11,*) 'pabsr=',pabsr,' pabsrr=',pabsrr,' pabsrrr=',pabsrrr
      write(11,*) 'pmx=',rdrd1,'gam=',gmgm1
      write(11,*) 'run=',run,' run1=',run1
      CLOSE(11)

 5000 CONTINUE
!C      pause
      STOP
      END



      SUBROUTINE sip5d(m,n,a1,a2,cc,a3,a4,r,kt,alphamax,w1,w2)

   	  integer m,n,kt,iacc
   	  integer k,itc,itc0,itc1
   	  integer i,ii,j,jj,ij,itj,ijb,ibj,ib,it,jb,n1,m1
      double precision a1(m,n),a2(m,n),a3(m,n),a4(m,n),cc(m,n),r(m,n)
	  double precision alpha,alphamax,alp(9),w1(m,n),w2(m,n)
      double precision eijb,eibj,fijb,fibj,rijb,ribj,sb,sc,sd
      double precision sbeijb,scfibj
!c
      data alp(1),alp(2),alp(3),alp(4),alp(5),alp(6),alp(7),alp(8) &
      ,alp(9)/1.0D0,0.625D0,0.25D0,0.875D0,0.5D0,0.125D0,0.75D0, &
      0.375D0,0.0D0/
!c

	  n1 = n + 1
	  m1 = m + 1

!c
	  itc = mod(kt-1,2) + 1
	  if (itc.eq.0)	itc = 2

	  itc0 = 2 - itc
	  itc1 = itc - 1

!CAcceleration parameter

	  iacc = mod(kt-1,18)
	  if (iacc.lt.0) iacc = iacc + 18

	  iacc = iacc/2 + 1
	  alpha = 1.0 - (1.0 - alphamax)**alp(iacc)

!CApproximate factorisation and inversion of the lower triangular matrix

	  do 10 jj=1,n

		j = itc0*jj + itc1*(n1-jj)
		jb = j - itc0 + itc1

		do 10 i=1,m

			ib = i - 1

			if (cc(i,j).eq.0.0) goto 100

			if ((jb.eq.0).or.(jb.eq.n1)) then
				eijb = 0.0D0
				fijb = 0.0D0
				rijb = 0.0D0

			else
				eijb = w1(i,jb)
				fijb = w2(i,jb)
				rijb = r(i,jb)

			endif

			if (ib.eq.0) then
				eibj = 0.0D0
				fibj = 0.0D0
				ribj = 0.0D0

			else
				eibj = w1(ib,j)
				fibj = w2(ib,j)
				ribj = r(ib,j)

			endif

			sb = (itc0*a1(i,j) + itc1*a4(i,j))/(1 + alpha*eijb)
			sc = a2(i,j)/(1 + alpha*fibj)

			sbeijb = sb*eijb
			scfibj = sc*fibj

!c			sd = 1.0D0 - sb*fijb - sc*eibj + alpha*(sbeijb + scfibj)
	        sd = cc(i,j) - sb*fijb - sc*eibj + alpha*(sbeijb + scfibj)

			w1(i,j) = (a3(i,j) - alpha*sbeijb)/sd
			w2(i,j) = (itc0*a4(i,j) + itc1*a1(i,j) - alpha*scfibj)/sd



!CInversion of the lower triangular matrix

			r(i,j) = (r(i,j) - sb*rijb - sc*ribj)/sd

			goto 10

  100       w1(i,j) = 0.0
            w2(i,j) = 0.0

   10 continue

!CInversion of the upper triangular matrix

	  do 20 jj=1,n

		j = itc0*(n1-jj) + itc1*jj
		jb = j + itc0 - itc1

		do 20 ii=1,m

			i = m1 - ii
			it = i + 1

			if ((jb.ne.0).and.(jb.ne.n1))	then
			    r(i,j) = r(i,j) - w2(i,j)*r(i,jb)
			endif

			if (it.ne.m1) then
				r(i,j) = r(i,j) - w1(i,j)*r(it,j)
		    endif


   20 continue

      return
      end


      SUBROUTINE sip9d(m,n,a1,a2,a3,a4,cc,a5,a6,a7,a8,r,kt,alphamax,w1,w2,w3,w4)
	  integer n,m,kt,iacc
	  integer itc,itc00,itc01,itc10,itc11
	  integer i,ii,j,jj,ij
	  integer it,itj,ijb,ibj,ib,jb
	  integer ibjb,itjb
	  integer z1,z2,z3,z4
	  double precision a1(m,n),a2(m,n),a3(m,n),a4(m,n),a5(m,n),a6(m,n)
	  double precision a7(m,n),a8(m,n),cc(m,n),r(m,n),w1(m,n)
	  double precision alpha,alphamax,alp(9),w2(m,n),w3(m,n),w4(m,n)
	  double precision fibjb,fijb,fitjb,gibjb,gijb,gitjb,hibjb
	  double precision hijb,hitjb,iibjb,iijb,iitjb,ribjb,rijb,ritjb

	  double precision fibj,gibj,hibj,iibj,ribj
	  double precision sa,sb,sc,sd,se
!c
      data alp(1),alp(2),alp(3),alp(4),alp(5),alp(6),alp(7),alp(8) &
      ,alp(9)/1.0D0,0.625D0,0.25D0,0.875D0,0.5D0,0.125D0,0.75D0, &
      0.375D0,0.0D0/
!c

	  n1 = n + 1
	  m1 = m + 1

!C    Alternate direction technique

	  itc = mod(kt-1,4) + 1
	  if (itc.eq.0) itc = 4

	  if (itc.eq.1) then
		itc00 = 1
		itc01 = 0
		itc10 = 1
		itc11 = 0
	  endif

	  if (itc.eq.2) then
		itc00 = 0
		itc01 = 1
		itc10 = 1
		itc11 = 0
	  endif

	  if (itc.eq.3) then
		itc00 = 1
		itc01 = 0
		itc10 = 0
		itc11 = 1
	  endif

	  if (itc.eq.4) then
		itc00 = 0
		itc01 = 1
		itc10 = 0
		itc11 = 1
	  endif

	  z1 = itc00*itc10
	  z2 = itc01*itc10
	  z3 = itc00*itc11
	  z4 = itc01*itc11

!CAcceleration parameter

	  iacc = mod(kt-1,36)
	  if (iacc.lt.0) iacc = iacc +36

	  iacc = iacc/4 + 1
	  alpha = 1.0 - (1.0 - alphamax)**alp(iacc)

!CApproximate factorisation and inversion of the lower triangular matrix

	  do 10 jj=1,n

		j = itc00*jj + itc01*(n1-jj)
		jb = j - itc00 + itc01

		do 10 ii=1,m

			i = itc10*ii + itc11*(m1-ii)
			ib = i - itc10 + itc11
			it = i + itc10 - itc11

			if (cc(i,j).eq.0.0) goto 100

			if ((jb.eq.0).or.(jb.eq.n1)) then
				fibjb = 0.0
				fijb = 0.0
				fitjb = 0.0
				gibjb = 0.0
				gijb = 0.0

				gitjb = 0.0
				hibjb = 0.0
				hijb = 0.0

				hitjb = 0.0
				iibjb = 0.0
				iijb = 0.0

				iitjb = 0.0
				ribjb = 0.0
				rijb = 0.0

				ritjb = 0.0
			else
				if ((ib.eq.0).or.(ib.eq.m1)) then
					fibjb = 0.0
					gibjb = 0.0
					hibjb = 0.0
					iibjb = 0.0
					ribjb = 0.0
				else
					fibjb = w1(ib,jb)
					gibjb = w2(ib,jb)
					hibjb = w3(ib,jb)
					iibjb = w4(ib,jb)
					ribjb = r(ib,jb)
				endif

				if ((it.eq.0).or.(it.eq.m1)) then
					fitjb = 0.0
					gitjb = 0.0
					hitjb = 0.0
					iitjb = 0.0
					ritjb = 0.0
				else
					fitjb = w1(it,jb)
					gitjb = w2(it,jb)
					hitjb = w3(it,jb)
					iitjb = w4(it,jb)
					ritjb = r(it,jb)
				endif

				fijb = w1(i,jb)
				gijb = w2(i,jb)
				hijb = w3(i,jb)
				iijb = w4(i,jb)
				rijb = r(i,jb)
			endif


			if ((ib.eq.0).or.(ib.eq.m1)) then
				fibj = 0
				gibj = 0
				hibj = 0
				iibj = 0
				ribj = 0
			else
				fibj = w1(ib,j)
				gibj = w2(ib,j)
				hibj = w3(ib,j)
				iibj = w4(ib,j)
				ribj = r(ib,j)
			endif

			sa = a1(i,j)*z1 + a6(i,j)*z2 + a3(i,j)*z3 + a8(i,j)*z4
			sb = (a2(i,j)*z1 + a7(i,j)*z2 + a2(i,j)*z3 + a7(i,j)*z4) - sa*fibjb
			sc = ((a3(i,j)*z1 + a8(i,j)*z2 + a1(i,j)*z3+a6(i,j)*z4) - sb*fijb)	/(1+alpha*fitjb)
			sd = ((a4(i,j)*z1 + a4(i,j)*z2 + a5(i,j)*z3 +a5(i,j)*z4) - sb*gijb - sa*hibjb - 2*alpha*sa*gibjb)/(1 + alpha*gibj)

			se = cc(i,j) - sd*fibj - sc*gitjb - sb*hijb - sa*iibjb + alpha*(sc*fitjb + sa*gibjb + sc*iitjb + sd*gibj)

			w1(i,j) = ((a5(i,j)*z1 + a5(i,j)*z2 + a4(i,j)*z3 + a4(i,j)*z4) - sb*iijb - sc*hitjb - alpha*sc*(2*iitjb + fitjb))/se
			w2(i,j) = ((a6(i,j)*z1 + a1(i,j)*z2 + a8(i,j)*z3 + a3(i,j)*z4) - sd*hibj - alpha*sd*gibj)/se
			w3(i,j) = ((a7(i,j)*z1 + a2(i,j)*z2 + a7(i,j)*z3 + a2(i,j)*z4) - sd*iibj)/se
			w4(i,j) = (a8(i,j)*z1 + a3(i,j)*z2 + a6(i,j)*z3 + a1(i,j)*z4)/se


!c Inversion of the lower triangular matrix

			r(i,j) = (r(i,j) - sa*ribjb - sb*rijb - sc*ritjb - sd*ribj)/se

			goto 10

  100       w1(i,j) = 0.0
            w2(i,j) = 0.0
			w3(i,j) = 0.0
            w4(i,j) = 0.0

   10 continue


!CInversion of the upper triangular matrix */

      do 20 jj=1,n

        j = itc00*(n1-jj) + itc01*jj
        jb = j + itc00 - itc01

        do 20 ii=1,m

            i = itc10*(m1-ii) + itc11*ii
            ib = i + itc10 - itc11
            it = i - itc10 + itc11

            if ((jb.ne.0).and.(jb.ne.n1)) then

                if ((ib.ne.0).and.(ib.ne.m1)) then
                    r(i,j) = r(i,j) - w4(i,j)*r(ib,jb)
                endif

                if ((it.ne.0).and.(it.ne.m1)) then
                    r(i,j) = r(i,j) - w2(i,j)*r(it,jb)
                endif

                r(i,j) = r(i,j) - w3(i,j)*r(i,jb)
            endif


            if ((ib.ne.0).and.(ib.ne.m1)) then
                r(i,j) = r(i,j) - w1(i,j)*r(ib,j)
            endif


   20 continue

      return
      end


      DOUBLE PRECISION FUNCTION S11ABF(X,IFAIL)
!C    MARK 5A REVISED - NAG COPYRIGHT 1976
!C    MARK 11.5(F77) REVISED. (SEPT 1985.)
!C    ARCSINH(X)
!C
!C    .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
!C    .. Local Scalars ..
      DOUBLE PRECISION                 LN2, T, XHI, Y
!C    .. Intrinsic Functions ..
      INTRINSIC                        ABS, LOG, SIGN, SQRT
!C    .. Data statements ..
!C    PRECISION DEPENDENT CONSTANTS
!C08   DATA XHI,LN2/1.0D+5,6.9314718D-1/
!C12   DATA XHI,LN2/1.0D+7,6.93147180560D-1/
!C14   DATA XHI,LN2/1.0D+8,6.9314718055995D-1/
      DATA XHI,LN2/1.0D+9,6.931471805599453D-1/
!C18   DATA XHI,LN2/1.0D+10,6.93147180559945309D-1/
!C    .. Executable Statements ..
!C
!C    NO FAILURE EXITS
      IFAIL = 0
      T = ABS(X)
!C    TEST LARGE RANGE
      IF (T.GT.XHI) GO TO 20
      T = X*X
!C
!C    TEST FOR MIDDLE RANGE
      IF (T.GT.1.0D0) GO TO 40
!C    EXPANSION ARGUMENT
      T = 2.0D0*T - 1.0D0
!C
!C     * EXPANSION (0012) *
!C
!C    EXPANSION (0012) EVALUATED AS Y(T)  --PRECISION 08E
!C08   Y = ((((((((+2.5373657D-6)*T-8.8982681D-6)*T+2.6962688D-5)
!C08  *    *T-1.0390000D-4)*T+4.2316523D-4)*T-1.8337913D-3)
!C08  *    *T+9.0042679D-3)*T-5.7366614D-2)*T + 9.3122986D-1
!C
!C    EXPANSION (0012) EVALUATED AS Y(T)  --PRECISION 12E
!C12   Y = (((((((((((((-6.08021584775D-9)*T+1.98638973110D-8)
!C12  *    *T-4.57508091809D-8)*T+1.58893891259D-7)
!C12  *    *T-5.83039413992D-7)*T+2.05819238929D-6)
!C12  *    *T-7.40199001696D-6)*T+2.74058634726D-5)
!C12  *    *T-1.05071718894D-4)*T+4.23002689824D-4)
!C12  *    *T-1.83345864783D-3)*T+9.00428855313D-3)
!C12  *    *T-5.73666392630D-2)*T + 9.31229859453D-1
!C
!C    EXPANSION (0012) EVALUATED AS Y(T)  --PRECISION 14E
!C14   Y = (((((((((((((((-5.8308026585469D-10)*T+1.8762742119651D-9)
!C14  *    *T-3.8936648507940D-9)*T+1.3296937569163D-8)
!C14  *    *T-4.9030635676319D-8)*T+1.6792346090382D-7)
!C14  *    *T-5.8053399097437D-7)*T+2.0520358645332D-6)
!C14  *    *T-7.4030149627429D-6)*T+2.7408018256297D-5)
!C14  *    *T-1.0507150365535D-4)*T+4.2300233069384D-4)
!C14  *    *T-1.8334586677591D-3)*T+9.0042885755717D-3)
!C14  *    *T-5.7366639262479D-2)*T + 9.3122985945271D-1
!C
!C    EXPANSION (0012) EVALUATED AS Y(T)  --PRECISION 16E
      Y = (((((((((((((((+1.810792296549804D-11 &
         *T-5.731943029121004D-11)*T+1.008344962167889D-10) &
         *T-3.394726871170490D-10)*T+1.299779213740398D-9) &
         *T-4.319978113584910D-9)*T+1.432753532351304D-8) &
         *T-4.863477336087045D-8)*T+1.670117348345774D-7) &
         *T-5.807433412373489D-7)*T+2.052474396638805D-6) &
         *T-7.402952157663977D-6)*T+2.740790473603819D-5) &
         *T-1.050715136470630D-4)*T+4.230023450529706D-4) &
         *T-1.833458667045431D-3)*T+9.004288574881897D-3
      Y = (Y*T-5.736663926249348D-2)*T + 9.312298594527122D-1
!C
!C    EXPANSION (0012) EVALUATED AS Y(T)  --PRECISION 18E
!C18   Y = (((((((((((((((+1.83000936086660710D-12
!C18  *    *T-5.74532833323253760D-12)*T+8.95787616116500070D-12)
!C18  *    *T-3.00291207083554898D-11)*T+1.20278345675996627D-10)
!C18  *    *T-3.94053306282758078D-10)*T+1.27690409672956527D-9)
!C18  *    *T-4.26028056137241518D-9)*T+1.43437981020129308D-8)
!C18  *    *T-4.86735767698085752D-8)*T+1.67004579212037421D-7)
!C18  *    *T-5.80728097040980303D-7)*T+2.05247631332341414D-6)
!C18  *    *T-7.40295567555544635D-6)*T+2.74079044411636393D-5)
!C18  *    *T-1.05071513207326607D-4)*T+4.23002345076007640D-4
!C18   Y = (((Y*T-1.83345866707041640D-3)*T+9.00428857488119850D-3)
!C18  *    *T-5.73666392624930609D-2)*T + 9.31229859452712177D-1
!C
!C
      S11ABF = X*Y
      RETURN
!C
!C    LARGE X
   20 S11ABF = SIGN(DLOG(T)+LN2,X)
      RETURN
!C
!C    MIDDLE X
   40 S11ABF = SIGN(DLOG(SQRT(T+1.0D0)+ABS(X)),X)
      RETURN
!C
      END



