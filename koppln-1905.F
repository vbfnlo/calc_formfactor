C========= subroutine KOPPLN  ==============================================
C
C  Adapted from RKWW by R.Kleiss
C
C  H --> gamma gamma decay added by D. Rainwater,
C        from HHG, Mar 97
c        included D5 operators by TFigy, June 04
C  H --> b bbar decay modified by D. Rainwater,
C        from HHG, Mar 98
C
c  H --> V*V* -> all added by T.Figy, Oct 2003
c
C  Calculate the fermion-fermion-boson and three-boson couplings
C  in the standard model, without colour factors.
C  also calculate mass**2 and mass*width of the bosons.
C  V   = FFB vector coupling
C  A   = FFB axial-vector coupling
C  CLR = FFB left- and righthanded couplings
C  B   = BBB coupling
C  The fermion indices are:  1: neutrino of electron (muon,tau)
C                            2: electron (muon,tau)
C                            3: up quark (charm,top)
C                            4: down quark (strange,bottom)
C  The boson indices are     1: photon
C                            2: Z0 boson
C                            3: W+ boson
C                            4: W- boson
C                            5: gluon
C                            6: Higgs
C  The helicity indices are  -1: lefthanded
C                            +1: righthanded
C  No Kobayashi-Maskawa mixing is implemented here.
C  The W+- are identified as outgoing from the BBB vertex.
C  XM2 = boson mass**2
C  XMG = boson mass*width
C  Also the branching ratios of W and Z decay are calculted in the SM
C
C  The fundamental parameters:
C  ALFA = The Q.E.D. fine structure constant
C  SIN2W = sin**2 of the weak mixing angle
C  ALFAS = the QCD fine structure constant
C  XMZ = mass of the Z0 boson in GeV
C  XMW = mass of the W+- bosons in GeV
C  XMH = mass of the Higgs boson in GeV
C  XMT = mass of the top quark in GeV
C
C  if XMW and SIN2W are <=0, they will be calculated from XMZ and GFermi
C
C  printout is steered by INFO:
C   <=0 :   no printout
C   <=1 :   printout of couplings only
C    =2 :   printout of input parameters and branching ratios
C   else:   printout of couplings and W/Z branching ratios
C
C===========================================================================
C
      SUBROUTINE KOPPLN(INFO)
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
      double complex  f_tau_t, f_tau_b, f_tau_W, F_t, F_b, F_W, dum
      double precision  beta0, beta1, gamma0, gamma1
      double precision  mllb
      double precision  zero
      parameter ( zero = 0.0d0 )
      integer  info, i, j, k, iflag
      logical ldebug, lanom
      parameter (lanom=.false.)
c
c INPUT
c
c      double precision  alfas, xmt, alfa, xmz, xmw, sin2w, xmh, gfermi
c
c      common /bkopin/   alfas, xmt, alfa, xmz, xmw, sin2w, xmh
c
c      data alfas  /0.1185d0/
c      data xmt    /175.d0/
c      data alfa   /-1.0d0/
c      data xmz    /91.188d0/
c      data xmw    /-1.d0/
c      data sin2w  /0.2315d0/
c      data xmh    /120.d0/
c      data gfermi /1.16639d-5/
      include 'global.inc'
      include 'koppln.inc'
      include 'an_couplings_4v.inc'
c
c OUTPUT
c
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A, XM2s(6),XMGs(6),
     & BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     & BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,BHGAM
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      COMMON /BKOPOUshort/ XM2s,XMGs
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM

c     ahvv(i,v1,v2) i = 1,2,3 are coefficients of eq 1 of 
c     Physics Letters B 591,297
      double complex ahvv(3,4,4)
      common/tensorhvv/ ahvv
      real*8 treefac, loopfac 
      common/lhcoup/ treefac, loopfac
c
c old version of output format
c
      double precision OCLR, OXM2, OXMG, OB(5,5,5),OV(4,5),OA(4,5)
      COMMON /BKOPPL/ OCLR(4,5,-1:1),OXM2(5),OXMG(5)
c     &                OB, OV, OA
C
C local variable
C
      double precision e,g,s,c,g2,z,w,q,
     & gwh, gzh, decw, gwne, gwud, gwtb, decz, gznn, gzee, gzuu, gzdd,
     & gztt, xgz, heavy, fh, y1,y2,x, xmwsq,xgw,xlw,xlz,dhzz,dhwwst,
     & ghww,ghzz,ghtt, xmb,xmc,xmtau, betac,betab,betat, qcdc,qcdb,qcdt,
     & ghcc,ghbb,ghtau,ghgg,ghgam,tau_t,eta_p,eta_m,tau_b,tau_w,xgh
      double precision xmMU,betaMU,ghmu,gh4f
      double precision BHGAMZ,GHGAMZ
      common/BRANCH2/BHGAMZ
      real*8 xmax,xmin
      double complex II_1,II_2,A_sum,A_w,A_t,A_b
      double precision cos2w,tan2w,cotw,lamb_b,lamb_t,lamb_w
c
      logical luseSM
      parameter(luseSM = .false.) ! .true. implies only SM
c
c PHASE SPACE REDUCTION FACTOR FOR HEAVY FERMIONIC DECAYS
c
      HEAVY(X,Y1,Y2) = ( 1.d0 - .5d0*(Y1**2 + Y2**2) - 
     &           0.5d0*(Y1**2 - Y2**2)**2 +
     &           3.d0*Y1*Y2*((X**2 - 1.d0)/(X**2 + 1.d0)) )*
     &           SQRT((1.d0 - Y1**2 - Y2**2)**2 - 4.d0*Y1**2*Y2**2)
c
c mass dependence for H--->W*W,Z*Z (From Gunion et al., Higgs Hunters Guide)
c
c      FH(X) = - (1.d0 - X**2)*(23.5d0*X**2 - 6.5d0 + 1.d0/X**2) 
c     &        - 3.d0*(1.d0 - 6.d0*X**2 + 4.d0*X**4)*LOG(X)
c     &        + 3.d0*(1.d0 - 8.d0*X**2 + 20.d0*X**4)/
c     &                                   SQRT(4.d0*X**2 - 1.d0)
c     &        * ACOS( (3.d0*X**2 - 1.d0)/(2.d0*X**3) )
c      PI = 4.d0*ATAN(1.d0)
      IFLAG = 0
c
c for 1/alpha_QED < 136 the Weak angle is calculated from the input
c values of alpha(assumed at scale m_W**2), m_Z and G_Fermi
c
      ldebug = .false.   ! alfa.eq.0.122431
      IF ( ALFA.GT.1.d0/136.d0 .and. sin2w.le.0 ) THEN
        XMWSQ = XMZ**2/2.d0 + SQRT(XMZ**4/4.d0 - 
     &              PI*ALFA*XMZ**2/SQRT(2.d0)/GFERMI) 
        SIN2W = 1.d0 - XMWSQ/XMZ**2
      ENDIF
c
c variable input: sin**2(theta_W) or W mass
c
      IF ( SIN2W.LE.zero .AND. XMW.LE.zero ) THEN
        SIN2W = 0.2312d0
        XMW = XMZ*SQRT(1.d0 - SIN2W)
      ELSEIF ( SIN2W.LE.zero ) THEN
        SIN2W = 1.d0 - (XMW/XMZ)**2
      ELSEIF ( XMW.LE.zero ) THEN
        XMW = XMZ*SQRT(1.d0 - SIN2W)
      ENDIF

      G = SQRT(4*PI*ALFAS)
      S = SQRT(SIN2W)
      C = SQRT(1.d0 - SIN2W)
c     initialize ahvv to zero
      A_sum=(0.0d0,0.0d0)
      dum =(0.0d0,0.0d0)
c      call fillanomhcoupl(1,A_sum,dum,loopfac)
      do i = 1,4
         do j = 1,4
            ahvv(1,i,j) = (0.0d0,0.0d0)
            ahvv(2,i,j) = (0.0d0,0.0d0)
            ahvv(3,i,j) = (0.0d0,0.0d0)
         enddo
      enddo
c     fill ahvv coefficients h-w-w and h-z-z
      if(lanom) call fillanomhcoupl(2,A_sum,dum,loopfac)
     
c
c G2 = E/S IS CALCULATED FROM GFERMI AND THE Z MASS
c
      G2 = SQRT(8.d0*GFERMI/SQRT(2.d0))*XMZ*C
      IF ( ALFA.GT.1.d0/136.d0 ) THEN
        E = SQRT(4.d0*PI*ALFA)
        G2 = E/S
      ELSEIF ( ALFA.GT.1.d0/138.d0 ) THEN
        E = SQRT(4.d0*PI*ALFA)
        IFLAG = 1
      ELSE
        E = G2*S
        ALFA = E**2/(4.d0*PI)
      ENDIF

      Z = G2/4.d0/C
      W = G2/SQRT(8.d0)
      Q = G2*C
c
c THE FFB VECTOR COUPLING
c
      V(1,1) = zero
      V(2,1) = -E
      V(3,1) = 2.d0*E/3.d0
      V(4,1) = -E/3.d0
      V(1,2) =  Z
      V(2,2) = -Z*(1.d0 - 4.d0*SIN2W)
      V(3,2) =  Z*(1.d0 - 8.d0/3.d0*SIN2W)
      V(4,2) = -Z*(1.d0 - 4.d0/3.d0*SIN2W)
      DO 1 K = 1,4
         DO 1 J = 3,4
            V(K,J) = W
 1    continue
      V(1,5) = zero
      V(2,5) = zero
      V(3,5) = G
      V(4,5) = G
c
c THE FFB AXIAL-VECTOR COUPLING
c
      DO 2 K = 1,4
         DO 2 J = 1,5,4
            A(K,J) = zero
 2    continue
      A(1,2) = -Z
      A(2,2) =  Z
      A(3,2) = -Z
      A(4,2) =  Z
      do K = 1,4
         do J = 3,4
            A(K,J) = -W
         enddo
      enddo
c
c THE FFB LEFT- AND RIGHTHANDED COUPLINGS
c
      DO 4 K = 1,4
      DO 4 J = 1,5
      CLR(K,J,0) = zero
      CLR(K,J,-1)= V(K,J)-A(K,J)
      CLR(K,J,1) = V(K,J)+A(K,J)
    4 continue
      if (ldebug) then
         print*,' electron L coupling '
         print*, (clr(2,i,-1),i=1,5)
         print*,' electron R coupling '
         print*, (clr(2,i,+1),i=1,5)
      endif
c
c THE BBB COUPLINGS
c
      DO 5 I = 1,6
         DO 5 J = 1,6
            DO 5 K = 1,6
               B(I,J,K) = zero
 5    continue
      B(3,4,1) = -E
      B(4,1,3) = -E
      B(1,3,4) = -E
      B(1,4,3) =  E
      B(4,3,1) =  E
      B(3,1,4) =  E
      B(3,4,2) = -Q
      B(4,2,3) = -Q
      B(2,3,4) = -Q
      B(2,4,3) =  Q
      B(4,3,2) =  Q
      B(3,2,4) =  Q
      B(5,5,5) =  G
c
c Higgs couplings to W and Z
c
      GWH = G2
      GZH = GWH/C**2
      B(2,2,6) = GZH
      B(2,6,2) = GZH
      B(6,2,2) = GZH
      B(3,4,6) = GWH
      B(4,3,6) = GWH
      B(3,6,4) = GWH
      B(4,6,3) = GWH
      B(6,3,4) = GWH
      B(6,4,3) = GWH
      
c      treefac = 1d0
      treefac = 1d0/sqrt(1+2d0*fphid_0/(sqrt(2d0)*GFERMI))
      ahvv(1,2,2) = treefac*B(6,2,2)*xmw *(-1.0d0,0.0d0) !a1 are sm hzz and hww couplings
      do i=3,4
         j = 7-i
         ahvv(1,i,j) = treefac*B(6,i,j)*xmw *(-1.0d0,0.0d0)
      enddo
c
c Calculate the total decay widths and branching ratios
c of the heavy bosons
c
c W first
c
      DECW = XMW/12.d0/PI
      GWNE = DECW*( V(1,3)**2 + A(1,3)**2 )
      GWUD = DECW*( V(3,3)**2 + A(3,3)**2 )*3.d0*(1.d0 + ALFAS/PI)
      GWTB = zero
      IF(XMT.GT.XMW) GOTO 31
      GWTB = GWUD*HEAVY( V(3,3)/A(3,3), XMT/XMW, zero )
   31 XGW  = 3.d0*GWNE + 2.d0*GWUD + GWTB
      BWNE = GWNE/XGW
      BWUD = GWUD/XGW
      BWTB = GWTB/XGW
c
c Z next
c
      DECZ = XMZ/12.d0/PI
      GZNN = DECZ*( V(1,2)**2 + A(1,2)**2 )
      GZEE = DECZ*( V(2,2)**2 + A(2,2)**2 )
      GZUU = DECZ*( V(3,2)**2 + A(3,2)**2 )*3.d0*(1.d0 + ALFAS/PI)
      GZDD = DECZ*( V(4,2)**2 + A(4,2)**2 )*3.d0*(1.d0 + ALFAS/PI)
      GZTT = zero
      IF(XMT.GT.XMZ/2.d0) GOTO 32
      GZTT = GZUU*HEAVY( V(3,2)/A(3,2) , XMT/XMZ , XMT/XMZ )
   32 XGZ  = 3.d0*GZNN + 3.d0*GZEE + 2.d0*GZUU + 3.d0*GZDD + GZTT
      BZNN = GZNN/XGZ
      BZEE = GZEE/XGZ
      BZUU = GZUU/XGZ
      BZDD = GZDD/XGZ
      BZTT = GZTT/XGZ
c
c     mass and width info 
c
      XM2(1) = zero
      XM2(2) = XMZ**2
      XM2(3) = XMW**2
      XM2(4) = XM2(3)
      XM2(5) = zero
      xm2(6) = xmh**2
c
      XMG(1) = zero
      XMG(2) = XMZ*XGZ
      XMG(3) = XMW*XGW
      XMG(4) = XMG(3)
      XMG(5) = zero
c     cannot compute xmg(6) at this point
c
c dummy xm2 and xmg for commbon block BKOPOUshort 
c
      do i = 1,6
         xm2s(i) = xm2(i)
      end do
      do i = 1,5
         xmgs(i) = xmg(i)
      end do     
c
c Next the Higgs decays to WW, ZZ, bb, cc, tau tau and gamma gamma
c
      XLZ = (XMZ/XMH)**2
      XLW = (XMW/XMH)**2
c      IF ( ABS(XLZ-0.25d0).LE.0.01d0 ) THEN
c        WRITE(6,21) 'Warning: H-->ZZ width not reliable near threshold'
c      ELSEIF ( ABS(XLW-0.25d0).LE.0.01d0 ) THEN
c        WRITE(6,21) 'Warning: H-->WW width not reliable near threshold'
c      ENDIF
c      DHZZ = GFERMI/(16.d0*PI*SQRT(2.d0)) * XMH**3
c      DHWWST = XMW**4 * GFERMI**2 /(16.d0*PI**3)
c      IF ( XLZ.LT.0.25d0 ) THEN
c        GHZZ = DHZZ*SQRT(1.d0 - 4.d0*XLZ)*
c     &              (12.d0*XLZ**2 - 4.d0*XLZ + 1.d0)
c      ELSEIF ( XLZ.LT.1.d0 ) THEN
c        GHZZ = DHWWST/4.d0 * (7.d0-40.d0/3.d0*SIN2W + 
c     &                          160.d0/9.d0*SIN2W**2)/C**4
c     &         * XMH*FH(XMZ/XMH)
c      ELSE
c        GHZZ = zero
c      ENDIF
c      print*,'ghzz=',ghzz
c
c     Higgs decays to WW and ZZ
c     subroutines below use  BKOPOUshort common block
c     First h -> Z*Z*
      call initbosid(2)
      call calZ(0.0d0,xmin)      
      call calZ(xmh**2,xmax)
c      xmin = 0.0d0              !no bw
c      xmax = xmh**2
      if (ldebug) print*,' xmin,xmax = ',xmin,xmax

      call quad2d(xmin,xmax,ghzz)
      if(ghzz.lt.0.0d0) ghzz = 0.0d0
!      print*,'For h -> Z*Z*: ghzz = ',ghzz
c      print*,'xmin,xmax=',xmin,xmax
c
c     h -> W*W*
      call initbosid(4)
      call calZ(0.0d0,xmin)      
      call calZ(xmh**2,xmax)
c      xmin = 0.0d0
c      xmax = xmh**2
      if (ldebug) print*,' xmin,xmax = ',xmin,xmax
         
      call quad2d(xmin,xmax,ghww)
      if(ghww.lt.0.0d0) ghww = 0.0d0
!      print*,'ghww = ',ghww 
c      print*,'xmin,xmax=',xmin,xmax
c
c
c      IF ( XLW.LT.0.25d0 ) THEN
c        GHWW = 2.d0*DHZZ*SQRT(1.d0 - 4.d0*XLW)*
c     &                 (12.d0*XLW**2 - 4.d0*XLW + 1.d0)
c      ELSEIF ( XLZ.LT.1.d0 ) THEN
c        GHWW = 3.d0*DHWWST*XMH*FH(XMW/XMH)
c      ELSE
c        GHWW = zero
c      ENDIF
c
c Now the Higgs partial decay widths 
c
      IF ( XMH.LT.60.d0 ) THEN
         XMH = 60.d0
      ENDIF

      XMB = 4.5d0                  ! decay into bbar first
      BETAB = SQRT(1.d0 - 4.d0*(XMB/XMH)**2)
      beta0 = 11.d0 - 10.d0/3.d0
      beta1 = 102.d0 - 190.d0/3.d0
      gamma0 = -8.d0
      gamma1 = -404.d0/3.d0 + 200.d0/9.d0
      mllb = 1.d0 + ( gamma0*beta1*log(2.d0*log(xmh/.2d0)) + 
     &               (gamma0*beta1 - gamma1*beta0) )
     &                / (2.d0*beta0**3 * 2.d0*log(xmh/.2d0) )
      QCDB = ( LOG(XMH/0.2d0) / LOG(2.d0*XMB/0.2d0) )**(gamma0/beta0) *
     &       ( 1.d0 + 34.d0/3.d0/beta0/log(xmh/.2d0) ) 
C                    new C factor (Higgs Hunter's Guide p. 20)
C     &       ( 1.d0 + 3.d0*ALFAS/PI )     !old C factor
      GHBB = 3.d0*GFERMI*XMB**2/(4.d0*PI*SQRT(2.d0)) * BETAB**3 * 
     &                      XMH * QCDB * mllb**2

      IF ( XMH.GT.2.d0*XMT ) THEN
        BETAT = SQRT(1.d0 - 4.d0*(XMT/XMH)**2)
        beta0 = 7.d0  !11.d0 - 12.d0/3.d0
        QCDT = ( LOG(2.d0*XMT/0.2d0) / LOG(XMH/0.2d0) )**(8.d0/beta0)*
     &       ( 1.d0 + 3.d0*ALFAS/PI )
        GHTT = 3.d0*GFERMI*XMT**2/(4.d0*PI*SQRT(2.d0)) * BETAT**3 * 
     &              XMH * QCDT
      ELSE
        GHTT = zero
      ENDIF

      XMC = 1.37d0                 ! now H ---> c cbar 
      BETAC = SQRT(1.d0 - 4.d0*(XMC/XMH)**2)
      QCDC = ( LOG(2.d0*XMC/0.2d0) / LOG(XMH/0.2d0) )**(24.d0/(33.d0 
     &                                                   - 2.d0*5.d0)) *
     &       ( 1.d0 + 3.d0*ALFAS/PI )
      GHCC = 3.d0*GFERMI*XMC**2/(4.d0*PI*SQRT(2.d0)) * BETAC**3 * 
     &                                XMH * QCDC

      XMtau = 1.777d0              !  and  H ----> tau+ tau-
      BETAt = SQRT(1.d0 - 4.d0*(XMtau/XMH)**2)
      GHtau = GFERMI*XMtau**2/(4.d0*PI*SQRT(2.d0)) * BETAt**3 * XMH 
c      
      XMmu = 0.106d0              !  and  H ----> mu+ mu-
      BETAmu = SQRT(1.d0 - 4.d0*(XMmu/XMH)**2)
      GHmu = GFERMI*XMmu**2/(4.d0*PI*SQRT(2.d0)) * BETAmu**3 * XMH 
c
c for  H --> gg  I use the simplified formula in the heavy top approx.
c
      GHGG  = 1.5d-10 * XMH**3 * (ALFAS/0.15d0)**2  
c
c for  H --> gamma gamma  I use p. 25 of Higgs Hunter's Guide,
c using t,b and W loops only
c Physics Letters B 318,1 eq 12 with cpodd d5 operator
c      GHGAM = (alfa*G2/xmw)**2 * (xmh/pi)**3 / 1024.d0
      ghgam = alfa*g2/(4.0d0*pi*xmw)
      tau_t = 4.d0*(xmt/xmh)**2
      if ( tau_t.lt.1.d0 ) then
         eta_p = 1.d0 + sqrt(1.d0 - tau_t)
         eta_m = 1.d0 - sqrt(1.d0 - tau_t)
         f_tau_t = -0.25d0 * ( log(eta_p/eta_m) + pi*(0d0,-1.d0) )**2
      else
         f_tau_t = ( asin(sqrt(1.d0/tau_t)) )**2
      end if
      F_t = -2.d0*tau_t* (1.d0 + (1.d0 - tau_t)*f_tau_t)

      tau_b = 4.d0*(xmb/xmh)**2
      if ( tau_b.lt.1.d0 ) then
         eta_p = 1.d0 + sqrt(1.d0 - tau_b)
         eta_m = 1.d0 - sqrt(1.d0 - tau_b)
         f_tau_b = -0.25d0 * ( log(eta_p/eta_m) + pi*(0d0,-1.d0) )**2
      else
         f_tau_b = ( asin(sqrt(1.d0/tau_b)) )**2
      end if
      F_b = -2.d0*tau_b* (1.d0 + (1.d0 - tau_b)*f_tau_b)

      tau_W = 4.d0*(xmw/xmh)**2
      if ( tau_W.lt.1.d0 ) then
         eta_p = 1.d0 + sqrt(1.d0 - tau_W)
         eta_m = 1.d0 - sqrt(1.d0 - tau_W)
         f_tau_W = -0.25d0 * ( log(eta_p/eta_m) + pi*(0d0,-1.d0) )**2
      else
         f_tau_W = ( asin(sqrt(1.d0/tau_W)) )**2
      end if
      F_W = 2.d0 + 3.d0*tau_W*(1.d0 + (2.d0 - tau_W)*f_tau_W)

      dum = 3.d0*(4.d0*F_t + F_b)/9.d0 + F_W

      dum = ghgam*dum
      ahvv(2,1,1) = dum 
c
      if(lanom) call fillanomhcoupl(2,A_sum,dum,loopfac)
c
      if(luseSM) then
         ghgam = xmh**3/(64.0d0*pi)*(dble(dum)**2+dimag(dum)**2)
      else
         ghgam = xmh**3/(64.0d0*pi) * 
     1             ( dble(ahvv(2,1,1))**2+dimag(ahvv(2,1,1))**2+
     2               dble(ahvv(3,1,1))**2+dimag(ahvv(3,1,1))**2 )
      endif
c      GHGAM = GHGAM * ( dble(dum)**2 + dimag(dum)**2 )
c
c for H -> gamma Z (on-shell)
c cp even and odd d5 operators
c
c     Using Higgs Hunters guide
c     t,b ,W loops are include in SM calculation
c     compute At,Ab,Aw
      cos2w = 1.0d0 - sin2w
      tan2w = sin2w/cos2w
      cotw = sqrt(1.0d0/tan2w)
c
c     b quark loop
      tau_b = 4.d0*(xmb/xmh)**2
      lamb_b = 4.0d0*(xmb/xmz)**2
      call compints(tau_b,lamb_b,II_1,II_2)
c   
      A_b = -3.0d0 *2.0d0*V(4,1)*V(4,2)*(2.0d0*c/(g2*E))*(II_1 - II_2)
      A_b = A_b/(s*c) 
c     top quark loop
      tau_t = 4.d0*(xmt/xmh)**2
      lamb_t = 4.0d0*(xmt/xmz)**2
      call compints(tau_t,lamb_t,II_1,II_2)
c   
      A_t = -3.0d0 *2.0d0*V(3,1)*V(3,2)*(2.0d0*c/(g2*E))*(II_1 - II_2)
      A_t = A_t/(s*c) 
c     W loop
      tau_w = 4.0d0*(xmw/xmh)**2
      lamb_w = 4.0d0*(xmw/xmz)**2
      call compints(tau_w,lamb_w,II_1,II_2)
      A_w = 4.0d0*(3.0d0-tan2w)*II_2
      A_w = A_w + ((1.0d0+2.0d0/tau_w)*tan2w
     $     -(5.0d0+2.0d0/tau_w))*II_1
      A_w = -cotw*A_w
      A_sum = A_w + A_t + A_b
      A_sum = A_sum * alfa*g2/(4.0d0*pi*xmw)
      ahvv(2,1,2) = -A_sum 
      ahvv(2,2,1) = -A_sum
c
c      print*,'A_sum=',A_sum
      if(lanom) call fillanomhcoupl(2,A_sum,dum,loopfac)

c     
      if(xmh.gt.xmz) then
         ghgamz = (dble(ahvv(2,1,2))**2+dimag(ahvv(2,1,2))**2 +
     1       dble(ahvv(3,1,2))**2+dimag(ahvv(3,1,2))**2) *
     2       (xmh**2-xmz**2)**3
         ghgamz = ghgamz/(32.0d0*pi*xmh**3) !check this
      else
         ghgamz = 0.0d0
      endif
!      print*,'For h-> gamma Z:ghgamz=',ghgamz

      if(luseSM) then
         ghgamz = 0.0d0
         goto 50
      endif

      if(lanom) then
         call flavor_sum(gh4f)
         print*,'For h-> (gamma/Z)*(gamma/Z)*: gh4f=',gh4f
c     replace ghzz with full calculation
         ghzz = gh4f
c     h-> gam ffbar
         if(xmh.gt.xmz) then
            call width_hgamff(ghgamz)
         else
            ghgamz = 0.0d0
         endif
         
         print*,'For h-> (Z/gamma)* gamma : ghgamz=',ghgamz
      endif
   
c
c
c H total decay width
c
 50   XGH = GHWW + GHZZ + GHBB + GHTT + GHCC + GHTAU + GHMU 
     $      + GHGG + GHGAM + GHGAMZ
c
c H branching ratios
c
      BHWW  = GHWW/XGH
      BHZZ  = GHZZ/XGH
      BHGG  = GHGG/XGH
      BHTT  = GHTT/XGH
      BHBB  = GHBB/XGH
      BHTAU = GHTAU/XGH
      BHMU  = GHMU/XGH
      BHCC  = GHCC/XGH
      BHGAM = GHGAM/XGH
      BHGAMZ = GHGAMZ/XGH
c      print*,'BHGAMZ=',bhgamz
c
c width info for higgs
c
      XMG(6) = XMH*XGH 
c
c dummy xmg for commbon block BKOPOUshort for higgs
c
      xmgs(6) = xmg(6)      
c
c print out the coupling constants
c
    9 FORMAT(/,
     &' ****************************************************'/
     &' * WARNING: chosen values of alpha_QED and W and Z  *'/
     &' * couplings may upset gauge theory cancellations!  *'/
     &' ****************************************************'//)
   10 FORMAT(20(1H-),a,20(1H-))
      IF(INFO.LE.0) THEN
        IF (IFLAG.EQ.1) write(6,9)
        RETURN
      ENDIF
      write(6,*) ''
      write(6,10) ' START OF KOPPEL INFO '
      IF (IFLAG.EQ.1) write(6,9)
      write(6,11) 1./ALFA,GFERMI,SIN2W,ALFAS,XMZ,XMW,XMH,XMT
   11 FORMAT(/,
     & ' I CALCULATE THE COUPLINGS IN THE STANDARD MODEL.'/,
     & ' THE INPUT PARAMETERS ARE :'/,
     & ' Q.E.D. FINE STRUCTURE CONSTANT  :  1/',F7.3/,
     & ' G_FERMI                         :',1X,F14.12,' GeV^-2'/,
     & ' SIN**2 OF THE WEAK MIXING ANGLE :',F8.5/,
     & ' Q.C.D. FINE STRUCTURE CONSTANT  :',F8.5/,
     & ' MASS OF THE Z0 BOSON            :',F8.3,' GEV'/,
     & ' MASS OF THE W+/- BOSONS         :',F8.3,' GEV'/,
     & ' MASS OF THE HIGGS BOSON         :',F8.3,' GEV'/,
     & ' MASS OF THE TOP QUARK           :',F8.3,' GEV'/,
     & ' (ALL OTHER FERMIONS HAVE ZERO MASS)' )
!      write(6,12) C,S,Z,W,Q
!   12 FORMAT(/' C,S,Z,W,Q =',5F10.6)
!      write(6,*)

      if (info.ne.2) then
         write(6,13) ((V(J,K),K=1,5),J=1,4)
         write(6,14) ((A(J,K),K=1,5),J=1,4)
         write(6,15) ((CLR(J,K,-1),K=1,5),J=1,4)
         write(6,16) ((CLR(J,K,1),K=1,5),J=1,4)
         write(6,17) (((B(I,J,K),K=1,6),J=1,6),I=1,6)
      endif
   13 FORMAT(' FFB VECTOR COUPLING ( ROW=F, COL= B )'/,4(5F10.6/))
   14 FORMAT(' FFB AXIAL-VECTOR COUPLING ( ROW=F, COL= B )'/,4(5F10.6/))
   15 FORMAT(' FFB LEFTHANDED COUPLING ( ROW=F, COL= B )'/,4(5F10.6/))
   16 FORMAT(' FFB RIGHTHANDED COUPLING ( ROW=F, COL= B )'/,4(5F10.6/))
   17 FORMAT(' BBB NON-ABELIAN COUPLINGS'/,6( 6(6F10.6/)/) )

      IF (INFO.LE.1) RETURN
!      WRITE(6,18) XGW,BWNE*100.d0,BWUD*100.d0,BWTB*100.d0,
!     &            XGZ,BZNN*100.d0,BZEE*100.d0,BZUU*100.d0,
!     &                BZDD*100.d0,BZTT*100.d0,
!     &            XGH,BHWW*100.d0,BHZZ*100.d0,BHGG*100.d0,BHTT*100.d0,
!     &                BHBB*100.d0,BHCC*100.d0,BHMU*100.d0,BHTAU*100.d0,
!     &            BHGAM*100.d0,BHGAMZ*100.0d0
!   18 FORMAT('THE WIDTHS OF THE HEAVY BOSONS, AND THE WAY',
!     & ' THEY ARE BUILT UP:'/,
!     & ' TOTAL WIDTH OF W+ OR W-       =',F7.3,' GEV'/,
!     & ' BR FOR LEPTON DECAY           =',F7.3,' %'/,
!     & ' BR FOR MASSLESS QUARK DECAY   =',F7.3,' %'/,
!     & ' BR FOR TOP-BOTTOM DECAY       =',F7.3,' %'/,
!     & ' TOTAL WIDTH OF Z0             =',F7.3,' GEV'/,
!     & ' BR FOR NEUTRINO DECAY         =',F7.3,' %'/,
!     & ' BR FOR LEPTON DECAY           =',F7.3,' %'/,
!     & ' BR FOR LIGHT UP DECAY         =',F7.3,' %'/,
!     & ' BR FOR DOWN DECAY             =',F7.3,' %'/,
!     & ' BR FOR TOP DECAY              =',F7.3,' %'/,
!     & ' TOTAL WIDTH OF HIGGS          =',F9.5,' GEV'/,
!     & ' BR FOR WW DECAY               =',F7.3,' %'/,
!     & ' BR FOR ZZ DECAY               =',F7.3,' %'/,
!     & ' BR FOR GLUE DECAY             =',F7.3,' %'/,
!     & ' BR FOR TOP DECAY              =',F7.3,' %'/,
!     & ' BR FOR B DECAY                =',F7.3,' %'/,
!     & ' BR FOR CHARM DECAY            =',F7.3,' %'/,
!     & ' BR FOR MU DECAY               =',F7.3,' %'/,
!     $ ' BR FOR TAU DECAY              =',F7.3,' %'/,
!     & ' BR FOR PHOTON DECAY           =',F7.3,' %'/
!     & ' BR FOR Photon Z decay         =',F7.3,' %')
!      write(6,19) (XM2(J),XMG(J),J=1,6)
!   19 FORMAT(/,' MASS**2 AND MASS*WIDTH VALUES'/,6(2F10.2/))
      write(6,*) ''
      write(6,20)
      write(6,*) ''
   20 FORMAT(20(1H-),' END OF KOPPEL INFO   ',20(1H-))
   21 FORMAT(//' ',63('*'),/2(' *',61(' '),'*'/),' *       ',A,5(' '),
     &       '*'/,2(' *',61(' '),'*'/),' ',63('*') )
c
c assign old output format variables 
c
      DO 41 I = 1,5
        OXM2(I) = XM2(I)
        OXMG(I) = XMG(I)
        DO 42 J = 1,4
          OV(J,I) = V(J,I)
          OA(J,I) = A(J,I)
          DO 42 K = -1,1
            OCLR(J,I,K) = CLR(J,I,K)
 42   CONTINUE
      DO 41 J = 1,5
         DO 41 K = 1,5
            OB(I,J,K) = B(I,J,K)
 41   CONTINUE
      if (ldebug) then
         print ('(a,5f10.3)'),' oxm2 = ',oxm2
         print ('(a,5f10.3)'),'  xm2 = ',(xm2(i),i=1,5)
         print ('(a,3f10.3)'),' ratio = ',(oxm2(i)/xm2(i),i=2,4)
         print ('(a,5f10.3)'),' oxmg = ',oxmg
         print ('(a,5f10.3)'),'  xmg = ',(xmg(i),i=1,5)
         print ('(a,3f10.3)'),' ratio = ',(oxmg(i)/xmg(i),i=2,4)
         DO I = 1,5
            do j = 1,4
               DO K = -1,1,2
                  print('(3i4,2f15.6)'),i,j,k,OCLR(J,I,K),CLR(J,I,K)
                  if (abs(CLR(J,I,K)).gt.1d-10) then
                     print*,' ratio = ',OCLR(J,I,K)/CLR(J,I,K)
                  endif
               enddo
            enddo
         enddo
      endif
c
c done
c
      if (ldebug) stop
      RETURN
      END
      
c                         
c     The subroutine below calculate the partial width for H-> W*W* semi-
c     analytically.
      
   

      function fun(xx,yy)
      implicit none

      real*8 lamb1,lamb2,xx,yy,fun,
     $     fun5de,fun5do,funsm5de,funsm
c     
      double precision g5hvv(2,1:4,1:4),lambda5
      common/hcoupl/g5hvv     ,lambda5
      real*8 xmh, gfermil
      parameter(gfermil = 1.16639d-5)
c      double precision  alfas, xmt, alfa, xmz, xmw, sin2w, xmh, gfermi
c      common /bkopin/   alfas, xmt, alfa, xmz, xmw, sin2w, xmh
      real*8 XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      real*8 PI,xmv,xmw
      parameter( PI = 3.141592653589793d0)
      real*8 lambda,constants
      real*8 Q2(2),prop1,prop2,DV(2:4,2)
      integer id,bos,bos1,bos2,i,j
      common /partid/ id
c      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
c      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      double complex ahvv(3,4,4)
      common/tensorhvv/ ahvv
      logical luseSM
      parameter(luseSM = .false.)
c     
      bos = id
      if(id.eq.4) then
         bos1 = 3
         bos2 = 4
      else
         bos1 = 2
         bos2 = 2
      endif
      xmv = sqrt(xm2s(bos1))
      xmh = sqrt(xm2s(6))
      xmw = sqrt(xm2s(4))
c      print*,'mh=',xmh
      call calQ2(xx,Q2(1))
c      Q2(1) = xx 
      lamb1 = Q2(1)/xmh**2
c      
      call calQ2(yy,Q2(2))
c      Q2(2) = yy
      lamb2 = Q2(2)/xmh**2
c
      call lambda1(1.0d0,lamb1,lamb2,lambda) !call 1 time for all
      funsm = lambda * (1.0d0 + lamb1**2 + lamb2**2 + 
     1       10.0d0*lamb1*lamb2 -  2.0d0 * lamb1 -  2.0d0 * lamb2)
c
      if(luseSM) then
         constants = gfermil * xmh**3/(sqrt(2.0d0)*(PI**3) * 8.0d0)
      else            
         constants = dble(ahvv(1,bos1,bos2))**2 + 
     1              dimag(ahvv(1,bos1,bos2))**2
         constants = constants * xmh**3/((PI**3) * 64.0d0 * xmv**4)
      endif

      funsm = funsm * constants
      
c     cp odd contribution
      if(luseSM) then
         fun = funsm
         return
      endif

      fun5do = lambda * (1.0d0 + lamb1**2 + lamb2**2 - 
     1     2.0d0*lamb1*lamb2 -  2.0d0 * lamb1 -  2.0d0 * lamb2)
     2     *lamb1*lamb2
c
      constants = g5hvv(2,bos1,bos2)**2 * 
     $     xmh**7/( 8.0d0 * xmv**4 * (PI**3) * (lambda5**2))
      fun5do = fun5do * constants
c     cp even contribution 

      fun5de = lambda * (1.0d0 + lamb1**2 + lamb2**2 + 
     1     4.0d0*lamb1*lamb2 -  2.0d0 * lamb1 -  2.0d0 * lamb2)
     2     *lamb1*lamb2
      constants = g5hvv(1,bos1,bos2)**2 * 
     $     xmh**7/( 8.0d0 * xmv**4 * (PI**3) * (lambda5**2))
      fun5de = fun5de*constants
c     interefence of SM and CP even 

      funsm5de = lambda * (1.0d0 - lamb1 - lamb2) 
     2     *lamb1*lamb2
c
c id = 1 for photon, 2 for Z , 3 for W , 4 for W
c
    
      constants = -g5hvv(1,bos1,bos2)*dble(ahvv(1,bos1,bos2))* 
     $     3.0d0*xmh**5/( 8.0d0 * xmv**4 * (PI**3) * lambda5)
      funsm5de = funsm5de * constants
      fun = funsm + funsm5de + fun5do + fun5de !sum 
c      fun = fun * xmgs(bos1)*xmgs(bos2) !for no BW map
c      do i = 2,4
c         do j=1,2 ! i is bos and j is 1 or 2 
c            Dv(i,j) = 1.0d0/((Q2(j) - xm2s(i))**2 + xmgs(i)**2)
c         enddo
c      enddo
c      
c      fun = fun * Dv(bos1,1) * Dv(bos2,2)

      return 
      end
c
      subroutine calQ2(x,Q2)
      implicit none
      real*8 XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      real*8 x,Q2
      integer id                !id = 3,4 W's 2 is Z
      common /partid/ id
c      for W bosons
      Q2 = xm2s(id) + xmgs(id) * dtan(x)
      return
      end
c
      subroutine calZ(QQ2,Z)
      implicit none
      real*8 XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      real*8 QQ2,Z
      integer id
      common /partid/ id

c      print*,'mv=',sqrt(xm2s(id))
      Z = datan((QQ2 - xm2s(id))/xmgs(id))
c      
      return
      end
c
c     
c
      subroutine lambda1(x,y,z,lambda)
      implicit none
      real*8 x,y,z,lambda2, lambda
c
      lambda2 = x**2 + y**2 + z**2 - 2.0d0*x*y - 2.0d0*y*z - 2.0d0*x*z
      if (lambda2.gt.0d0) then
         lambda = sqrt(lambda2)
      else
         lambda = 0d0
      endif
      return
      end

      subroutine quad2d(x1,x2,ss)
      implicit none
      real*8 eps
      parameter(eps = -1.0d-5)
      real*8 ss,x1,x2,h,gaus2
      external h
      external gaus2
c
      ss = gaus2(h,x1,x2,eps)   !partial width for h-> V*V* -> all
      return
      end
c
      function ff(yy)
      implicit none
      real*8 ff,yy,fun,x,y
      common/xy/ x,y

      y = yy
      
      ff = fun(x,y)
      return
      end
c
      function h(xx)
      implicit none
      real*8 eps
      parameter(eps = 1.0d-7)
      real*8 h,xx,y1,y2,x,y,ff
      common /xy/ x,y
      real*8 ss,gaus
      external ff
      external gaus
      real*8 Q2,q2max
      double precision  xmh
      real*8 XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      integer id,boson
      common /partid/ id
c     
      xmh = sqrt(xm2s(6))
      x = xx
      call calQ2(x,Q2)
c      Q2 = x
      q2max = (xmh - sqrt(Q2))**2
      call calZ(q2max,y2) 
c      y2 = q2max                !no BW map
      call calZ(0.0d0,y1)
c      y1 = 0.0d0

      ss = gaus(ff,y1,y2,eps)
      if(id.eq.2) then          !symmetry factor of h->ZZ
         h =ss/2.0d0
      else
         h = ss
      endif

      return
      end

      subroutine initbosid(boson)
      integer id,boson
      common /partid/ id

      id = boson
      end
cccc
c
c
* THIS IS ITERATIVE INTEGRATION PROCEDURE                               
* ORIGINATES  PROBABLY FROM CERN LIBRARY                                
* IT SUBDIVIDES INEGRATION RANGE UNTIL REQUIRED PRECISION IS REACHED    
* PRECISION IS A DIFFERENCE FROM 8 AND 16 POINT GAUSS ITEGR. RESULT     
* EEPS POSITIVE TREATED AS ABSOLUTE PRECISION                           
* EEPS NEGATIVE TREATED AS RELATIVE PRECISION                           
      DOUBLE PRECISION FUNCTION GAUS(F,A,B,EEPS)                        
*     *************************                                         
      IMPLICIT REAL*8(A-H,O-Z)                                          
      IMPLICIT INTEGER(I-N)                                          
      EXTERNAL F                                                        
      DIMENSION W(12),X(12)                                             
      DATA CONST /1.0D-19/                                              
      DATA W                                                            
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,   
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,   
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,   
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/   
      DATA X                                                            
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,   
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,   
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,   
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/   
      EPS=dABS(EEPS)                                                    
      DELTA=CONST*dABS(A-B)                                             
      GAUS=0.d0                                                         
      AA=A                                                              
    5 Y=B-AA                                                            
      IF(dABS(Y) .LE. DELTA) RETURN                                     
    2 BB=AA+Y                                                           
      C1=0.5d0*(AA+BB)                                                  
      C2=C1-AA                                                          
      S8=0.d0                                                           
      S16=0.d0                                                          
      DO 1 I=1,4                                                        
      U=X(I)*C2                                                         
    1 S8=S8+W(I)*(F(C1+U)+F(C1-U))                                      
      DO 3 I=5,12                                                       
      U=X(I)*C2                                                         
    3 S16=S16+W(I)*(F(C1+U)+F(C1-U))                                    
      S8=S8*C2                                                          
      S16=S16*C2                                                        
      IF(EEPS.LT.0D0) THEN                                              
        IF(dABS(S16-S8) .GT. EPS*dABS(S16)) GO TO 4                     
      ELSE                                                              
        IF(dABS(S16-S8) .GT. EPS) GO TO 4                               
      ENDIF                                                             
      GAUS=GAUS+S16                                                     
      AA=BB                                                             
      GO TO 5                                                           
    4 Y=0.5d0*Y                                                         
      IF(dABS(Y) .GT. DELTA) GOTO 2                                     
      PRINT 7                                                           
      GAUS=0.d0                                                         
      RETURN                                                            
    7 FORMAT(1X,36HGAUS  ... TOO HIGH ACCURACY REQUIRED)                
      END                                                               
                                                                        
      DOUBLE PRECISION FUNCTION GAUS2(F,A,B,EEPS)                       
*     *************************                                         
      IMPLICIT REAL*8(A-H,O-Z)                                          
      IMPLICIT INTEGER(I-N)                                          
      EXTERNAL F                                                        
      DIMENSION W(12),X(12)                                             
      DATA CONST /1.0D-19/                                              
      DATA W                                                            
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,   
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,   
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,   
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/   
      DATA X                                                            
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,   
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,   
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,   
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/   
      EPS=dABS(EEPS)                                                    
      DELTA=CONST*dABS(A-B)                                             
      GAUS2=0.D0                                                        
      AA=A                                                              
    5 Y=B-AA                                                            
      IF(dABS(Y) .LE. DELTA) RETURN                                     
    2 BB=AA+Y                                                           
      C1=0.5d0*(AA+BB)                                                  
      C2=C1-AA                                                          
      S8=0.d0                                                           
      S16=0.d0                                                          
      DO 1 I=1,4                                                        
      U=X(I)*C2                                                         
    1 S8=S8+W(I)*(F(C1+U)+F(C1-U))                                      
      DO 3 I=5,12                                                       
      U=X(I)*C2                                                         
    3 S16=S16+W(I)*(F(C1+U)+F(C1-U))                                    
      S8=S8*C2                                                          
      S16=S16*C2                                                        
      IF(EEPS.LT.0D0) THEN                                              
        IF(dABS(S16-S8) .GT. EPS*dABS(S16)) GO TO 4                     
      ELSE                                                              
        IF(dABS(S16-S8) .GT. EPS) GO TO 4                               
      ENDIF                                                             
      GAUS2=GAUS2+S16                                                   
      AA=BB                                                             
      GO TO 5                                                           
    4 Y=0.5d0*Y                                                         
      IF(dABS(Y) .GT. DELTA) GOTO 2                                     
      PRINT 7                                                           
      GAUS2=0.D0                                                        
      RETURN                                                            
    7 FORMAT(1X,36HGAUS2 ... TOO HIGH ACCURACY REQUIRED)                
      END                                                               

      
C This subroutine computes the parametric integrals in the 
c Higgs Hunters Guide on pg 29
      subroutine compints(a,b,I1,I2)
      implicit none
      real*8 a,b
      double complex I1,I2
      double complex f,g
      external f
      external g

      I1 = a*b/(2.0d0*(a-b)) + a**2 * b**2/(2.0*(a-b)**2) * (f(a)-f(b))
     $     + a**2 * b/((a-b)**2) *(g(a) - g(b))
      I2 = -a*b/(2.0d0*(a -b))* (f(a) - f(b))
      return
      end
c
      double complex function g(aa)
      implicit none
      real*8 aa,nup,num
      real*8 pi
      parameter(pi =3.141592653589793d0) 

      if(aa.ge.1.0d0) then
         g = sqrt(aa -1.0d0)*ASIN(1.0d0/sqrt(aa))
      else
         nup = 1.0d0 + sqrt(1.0d0 - aa)
         num = 1.0d0 - sqrt(1.0d0 - aa)
         g = 0.5d0* sqrt(1.0d0 - aa)*(log(nup/num)*(1.d0,0.d0)
     $        + pi*(0.d0,-1.d0))
      endif
      return
      end
c
      double complex function f(aa)
      implicit none
      real*8 aa,nup,num
      real*8 pi
      parameter(pi =3.141592653589793d0) 
      
      if(aa.ge.1.0d0) then
         f = (ASIN(sqrt(1.0d0/aa)))**2
      
      else
         nup = 1.0d0 + sqrt(1.0d0 - aa)
         num = 1.0d0 - sqrt(1.0d0 - aa)
         f = -0.25d0*(log(nup/num)
     $        + pi*(0.d0,-1.d0))**2
      endif
      return
      end



c
c	all couplings, widths, and masses are computed from 
c	values stored in common/ bkopout/
c
c	edited by Barbara Jager on March 10, 2005

	subroutine ctrans
	
c***********************************************************************
c 	This subroutine sets up the couplings, masses and widths 
c	needed by HELAS as they are given in koppln.f
c***********************************************************************

	implicit none      
 
c common blocks which need to be filled for HELAS (masses, widths, couplings):	
	
	include 'coupl.inc'
			
      double precision g2,ee2,sw,cw,vv,sc2,pi		
      parameter ( pi = 3.141592653589793d0 )
      logical loutput
      parameter (loutput = .false.)
			
c common block with input values (only needed for top mass xmt):

      double precision  alphas, xmt, alfa, xmz, xmw,s2w, xmh
      common /bkopin/   alphas, xmt, alfa, xmz, xmw, s2w, xmh
	
c common block which has been filled in koppln.f and contains 
c FFB, BBB, and Higgs couplings (but not in form needed by HELAS):
	
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
           

c fill basic parameters:

	G      =  B(5,5,5)
	GG(1)  = -G
	GG(2)  = -G
	alfas  =  G**2/(pi*4.d0)  ! strong coupling at fixed ren. scale MZ

	ee     = -CLR(2,1,-1)
	ee2    =  ee**2
	alpha  =  ee2/(pi*4.d0)
	
	g2     = sqrt(2.d0)*CLR(2,3,-1)
	cw     = B(2,4,3)/g2
	sin2w  = 1.d0 - cw**2
	sw     = sqrt(sin2w)
	
	
c fill masses:

	hmass = sqrt(xm2(6))
	wmass = sqrt(xm2(3))
	zmass = sqrt(xm2(2))
	amass = 0.d0
	tmass = xmt
c	tmass = 175.d0
c	bmass = 0.d0
	bmass = 4.5d0	! value given in qqhqq.f
c	bmass = 4.7d0	! value from switchmom.f here 	 
	cmass = 0.d0
	lmass = 0.d0
	
	
c  fill widths:			

	hwidth = xmg(6)/hmass 
	wwidth = xmg(3)/wmass
	zwidth = xmg(2)/zmass
	awidth = 0.d0	
	lwidth =  0d0
	twidth = 1.6d0	! no top emerges in WW-prod. -> value here 'arbitrary'
	
c  auxiliary variables:	
	
	vv  = 2.d0*wmass*sw/ee
	sc2 = sin2w*(1.d0-sin2w)
	
c  fill gauge couplings: 
c
c	fermion-fermion-vector couplings: 

	gal(1) = -dcmplx(clr(2,1,-1))	! photon
	gal(2) = -dcmplx(clr(2,1, 1))
	gau(1) = -dcmplx(clr(3,1,-1))	
	gau(2) = -dcmplx(clr(3,1, 1))	
	gad(1) = -dcmplx(clr(4,1,-1))	
	gad(2) = -dcmplx(clr(4,1, 1))	
	
	gwf(1) = -dcmplx(clr(1,3,-1))	! W+- boson (same for all fermions)
	gwf(2) = -dcmplx(clr(1,3, 1))	! zero
	
	gzn(1) = -dcmplx(clr(1,2,-1))	! Z boson
	gzn(2) = -dcmplx(clr(1,2, 1))	! zero	
	gzl(1) = -dcmplx(clr(2,2,-1))
	gzl(2) = -dcmplx(clr(2,2, 1))	
	gzu(1) = -dcmplx(clr(3,2,-1))
	gzu(2) = -dcmplx(clr(3,2, 1))	
	gzd(1) = -dcmplx(clr(4,2,-1))
	gzd(2) = -dcmplx(clr(4,2, 1))
	
c  	vector boson couplings:

	gw   = ee/sw
	gwwa = -b(3,4,1)
	gwwz = -b(3,4,2)
	
c	gauge - higgs boson couplings:

	gwwh  = dcmplx(b(3,4,6)*wmass)	
	gzzh  = dcmplx(b(2,2,6)*wmass)	
	ghhh  = dcmplx(-hmass**2/vv*3d0)
	gwwhh = dcmplx( ee2/sin2w*0.5d0 )
        gzzhh = dcmplx( ee2/sc2*0.5d0)
        ghhhh = ghhh/vv
	
	
c-----------------------------------------------

c for reference:
c	info on couplings returned by ctrans:

      if (loutput) then
	write(6,*) 'couplings returned by ctrans:'
	write(6,*) 'alpha,ee, sin2w, alfas,g,gg:'
	write(6,*)  alpha,ee, sin2w, alfas,g,gg
	write(6,*) 
	write(6,*) 'hmass, wmass, zmass, amass:'
	write(6,*)  hmass, wmass, zmass, amass
	write(6,*) 
	write(6,*) 'tmass, bmass, lmass, cmass:'
	write(6,*)  tmass, bmass, lmass, cmass
	write(6,*) 
	write(6,*) 'hwidth, wwidth, zwidth:' 
	write(6,*)  hwidth, wwidth, zwidth
	write(6,*) 
	write(6,*) 'twidth, lwidth, awidth:'
	write(6,*)  twidth, lwidth, awidth
	write(6,*) 
	write(6,*) 'gal   , gad   , gau   , gwf:'
	write(6,*)  gal   , gad   , gau   , gwf
	write(6,*) 
	write(6,*) 'gzn   , gzl   , gzd   , gzu:'
	write(6,*)  gzn   , gzl   , gzd   , gzu
	write(6,*) 
	write(6,*) 'gw, gwwa, gwwz:'
	write(6,*)  gw, gwwa, gwwz
	write(6,*) 
	write(6,*) 'gwwh, gzzh, gwwhh, gzzhh, ghhh, ghhhh:'
	write(6,*)  gwwh, gzzh, gwwhh, gzzhh, ghhh, ghhhh
	write(6,*) 
      endif
	
c---------------------------------------------------------------------------	
	

	return
	end
	
	
	
	
	
	
	
