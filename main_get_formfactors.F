
      program get_formfactors

c computes form factors for anomalous quartic gauge couplings from vector boson scattering

      implicit none
      include "coupl.inc"
      include "an_couplings_4v.inc"

      real*8 pi
      parameter(pi =3.141592653589793d0)     
      integer k

      real*8 real_pwave_ww, real_pwave_wpwp, real_pwave_wz, real_pwave_zz, real_pwave_wwzz, real_pwave_zzww
      external real_pwave_ww, real_pwave_wpwp, real_pwave_wz, real_pwave_zz, real_pwave_wwzz, real_pwave_zzww
      real*8 real_pwave_aa, real_pwave_zzaa, real_pwave_wwaa, real_pwave_wwza, real_pwave_zaaa
      external real_pwave_aa, real_pwave_zzaa, real_pwave_wwaa, real_pwave_wwza, real_pwave_zaaa
      real*8 real_pwave_aaww, real_pwave_zaww, real_pwave_zaza, real_pwave_aaza, real_pwave_aazz, real_pwave_zazz, real_pwave_zzza
      external real_pwave_aaww, real_pwave_zaww, real_pwave_zaza, real_pwave_aaza, real_pwave_aazz,real_pwave_zazz,real_pwave_zzza
      real*8 real_pwave_wa, real_pwave_wzwa, real_pwave_wawz
      external real_pwave_wa, real_pwave_wzwa, real_pwave_wawz

      real*8 im_pwave_ww, im_pwave_wpwp, im_pwave_wz, im_pwave_zz, im_pwave_wwzz, im_pwave_zzww
      external im_pwave_ww, im_pwave_wpwp, im_pwave_wz, im_pwave_zz, im_pwave_wwzz, im_pwave_zzww
      real*8 im_pwave_aa, im_pwave_zzaa, im_pwave_wwaa, im_pwave_wwza, im_pwave_zaaa
      external im_pwave_aa, im_pwave_zzaa, im_pwave_wwaa, im_pwave_wwza, im_pwave_zaaa
      real*8 im_pwave_aaww, im_pwave_zaww, im_pwave_zaza, im_pwave_aaza, im_pwave_aazz, im_pwave_zazz, im_pwave_zzza
      external im_pwave_aaww, im_pwave_zaww, im_pwave_zaza, im_pwave_aaza, im_pwave_aazz,im_pwave_zazz,im_pwave_zzza
      real*8 im_pwave_wa, im_pwave_wzwa, im_pwave_wawz
      external im_pwave_wa, im_pwave_wzwa, im_pwave_wawz

      real*8 ffscale_ww, ffscale_wpwp, ffscale_wz, ffscale_zz, ffscale_wwzz, ffscale_zzww
      real*8 ffscale_aa, ffscale_zzaa, ffscale_wwaa, ffscale_wwza, ffscale_zaaa, ffscale_zzza, ffscale_zaza
      real*8 ffscale_wa, ffscale_wzwa
      
      real*8 ffscale_ww_diag, ffscale_wpwp_diag, ffscale_wz_diag, ffscale_zz_diag
      real*8 ffscale_aa_diag, ffscale_za_diag, ffscale_wa_diag

      real*8 ffscale_diag_q2, ffscale_diag_q1, ffscale_diag_q0
      real*8 ffscale_diag_q1_a, ffscale_diag_q0_a

      real*8 unitviol_ww, unitviol_wpwp, unitviol_wz, unitviol_zz, unitviol_wwzz, unitviol_zzww
      real*8 unitviol_aa, unitviol_zzaa, unitviol_wwaa, unitviol_wwza, unitviol_zaaa, unitviol_zzza, unitviol_zaza
      real*8 unitviol_wa, unitviol_wzwa

      real*8 unitviol_ww_diag, unitviol_wpwp_diag, unitviol_wz_diag, unitviol_zz_diag
      real*8 unitviol_aa_diag, unitviol_za_diag, unitviol_wa_diag
      

      real*8 unitviol_diag_q2, unitviol_diag_q1, unitviol_diag_q0
      real*8 unitviol_diag_q1_a, unitviol_diag_q0_a

      real*8 searchff_maxpwave_ind, searchff_diag_ind, searchff_diag_q0, searchff_diag_q1_a, searchff_diag_q0_a
      external searchff_maxpwave_ind, searchff_diag_ind, searchff_diag_q0, searchff_diag_q1_a, searchff_diag_q0_a

      ! transfer center-of-masss energy, partial wave and helicity information
      ! to amplitude subroutine
      real*8 sc
      common /cms/ sc
      integer pw
      INTEGER NHEL(4)
      common /partwave/ pw, nhel

      pw=0        ! which partial wave should be computed? (0,1,2)

      if (pw.gt.0) then
            print*, ""
            print*, ""
            print*, "Warning:"
            print*, "Not all helicity combinations are tested for higher partial waves!"
            print*, ""
            print*, ""
      endif

      call readinput 



      ! setting center-of-mass energy to max. value
      sc=sqrts**2



      print*, ""
      print*, ""
      print*, ""
      print*, ""
      print*, "---------------------------------------------------------"
      print*, ""
      print*, "Starting calculation..."
      print*, ""
      print*, "---------------------------------------------------------"
      print*, ""


      open(11,file="ff_output.dat")

      write(11,*) "! Maximal center-of-mass energy in TeV:"
      write(11,"(A,F9.1,A)") " sqrt_s =", sqrts/1000d0
      write(11,*) ""
      write(11,*) "! The following anomalous couplings to / Lambda^4 have been set != 0 :"
      if (fs0_0.ne.0d0) write(11,"(A,E12.3,A)") " f_S0 =", fs0_0, " GeV^-4"
      if (fs1_0.ne.0d0) write(11,"(A,E12.3,A)") " f_S1 =", fs1_0, " GeV^-4"
      if (fM0_0.ne.0d0) write(11,"(A,E12.3,A)") " f_M0 =", fm0_0, " GeV^-4"
      if (fM1_0.ne.0d0) write(11,"(A,E12.3,A)") " f_M1 =", fm1_0, " GeV^-4"
      if (fM2_0.ne.0d0) write(11,"(A,E12.3,A)") " f_M2 =", fm2_0, " GeV^-4"
      if (fM3_0.ne.0d0) write(11,"(A,E12.3,A)") " f_M3 =", fm3_0, " GeV^-4"
      if (fM4_0.ne.0d0) write(11,"(A,E12.3,A)") " f_M4 =", fm4_0, " GeV^-4"
      if (fM5_0.ne.0d0) write(11,"(A,E12.3,A)") " f_M5 =", fm5_0, " GeV^-4"
      if (fM6_0.ne.0d0) write(11,"(A,E12.3,A)") " f_M6 =", fm6_0, " GeV^-4"
      if (fM7_0.ne.0d0) write(11,"(A,E12.3,A)") " f_M7 =", fm7_0, " GeV^-4"
      if (fT0_0.ne.0d0) write(11,"(A,E12.3,A)") " f_T0 =", ft0_0, " GeV^-4"
      if (fT1_0.ne.0d0) write(11,"(A,E12.3,A)") " f_T1 =", ft1_0, " GeV^-4"
      if (fT2_0.ne.0d0) write(11,"(A,E12.3,A)") " f_T2 =", ft2_0, " GeV^-4"
      if (fT5_0.ne.0d0) write(11,"(A,E12.3,A)") " f_T5 =", ft5_0, " GeV^-4"
      if (fT6_0.ne.0d0) write(11,"(A,E12.3,A)") " f_T6 =", ft6_0, " GeV^-4"
      if (fT7_0.ne.0d0) write(11,"(A,E12.3,A)") " f_T7 =", ft7_0, " GeV^-4"
      if (fT8_0.ne.0d0) write(11,"(A,E12.3,A)") " f_T8 =", ft8_0, " GeV^-4"
      if (fT9_0.ne.0d0) write(11,"(A,E12.3,A)") " f_T9 =", ft9_0, " GeV^-4"
      if (fb_0.ne.0d0) write(11,"(A,E12.3,A)") " f_B =", fb_0, " GeV^-2"
      if (fw_0.ne.0d0) write(11,"(A,E12.3,A)") " f_W =", fw_0, " GeV^-2"
      if (fwww_0.ne.0d0) write(11,"(A,E12.3,A)") " f_WWW =", fwww_0, " GeV^-2"
      if (fww_0.ne.0d0) write(11,"(A,E12.3,A)") " f_WW =", fWw_0, " GeV^-2"
      if (fphid_0.ne.0d0) write(11,"(A,E12.3,A)") " f_PHID =", fphid_0, " GeV^-2"
      if (fbtilde_0.ne.0d0) write(11,"(A,E12.3,A)") " f_BT =", fbtilde_0, " GeV^-2"
      if (fwwwtilde_0.ne.0d0) write(11,"(A,E12.3,A)") " f_WWWT =", fwwwtilde_0, " GeV^-2"
      if (fwwtilde_0.ne.0d0) write(11,"(A,E12.3,A)") " f_WWT =", fWwtilde_0, " GeV^-2"
      if (fdwtilde_0.ne.0d0) write(11,"(A,E12.3,A)") " f_DWT =", fdwtilde_0, " GeV^-2"
      write(11,*) ""

      ! WW -> WW
      ! due to singularity in t-channel photon exchange integration range 0.1° < theta < 180°
      ffscale_ww = searchff_maxpwave_ind(real_pwave_ww, "W+W- -> W+W-", 0.1d0/180d0*pi, unitviol_ww)
      write(11,*) "! Form factor from on-shell W+W- -> W+W- scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_WWWW =", ffscale_ww

      ! W+W+ -> W+W+
      ! due to singularity in t-channel photon exchange integration range 0.1° < theta < 180°
      ffscale_wpwp = searchff_maxpwave_ind(real_pwave_wpwp, "W+W+ -> W+W+", 0.1d0/180d0*pi, unitviol_wpwp)
      write(11,*) "! Form factor from on-shell W+W+ -> W+W+ scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_WPWP =", ffscale_wpwp

      ! WW -> ZZ
      ffscale_wwzz = searchff_maxpwave_ind(real_pwave_wwzz, "W+W- -> ZZ  ", 0d0/180d0*pi, unitviol_wwzz)
      write(11,*) "! Form factor from on-shell W+W- -> ZZ scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_WWZZ =", ffscale_wwzz

      ! WZ -> WZ
      ffscale_wz = searchff_maxpwave_ind(real_pwave_wz, "W+Z  -> W+Z ", 0d0/180d0*pi, unitviol_wz)
      write(11,*) "! Form factor from on-shell W+Z -> W+Z scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_WZWZ =", ffscale_wz

      ! ZZ -> ZZ
      ffscale_zz = searchff_maxpwave_ind(real_pwave_zz, "ZZ   -> ZZ  ", 0d0/180d0*pi, unitviol_zz)
      write(11,*) "! Form factor from on-shell ZZ -> ZZ scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_ZZZZ =", ffscale_zz


      ! Channels with external photons:

      ! WW -> AA
      ffscale_wwaa = searchff_maxpwave_ind(real_pwave_wwaa, "W+W- -> AA  ", 0d0/180d0*pi, unitviol_wwaa)
      write(11,*) "! Form factor from on-shell W+W- -> AA scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_WWAA =", ffscale_wwaa

      ! WW -> ZA
      ffscale_wwza = searchff_maxpwave_ind(real_pwave_wwza, "W+W- -> ZA  ", 0d0/180d0*pi, unitviol_wwza)
      write(11,*) "! Form factor from on-shell W+W- -> ZA scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_WWZA =", ffscale_wwza

      ! WA -> WA
      ffscale_wa = searchff_maxpwave_ind(real_pwave_wa, "W+A  -> W+A ", 0d0/180d0*pi, unitviol_wa)
      write(11,*) "! Form factor from on-shell W+A -> W+A scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_WAWA =", ffscale_wa

      ! WZ -> WA
      ffscale_wzwa = searchff_maxpwave_ind(real_pwave_wzwa, "W+Z  -> W+A ", 0d0/180d0*pi, unitviol_wzwa)
      write(11,*) "! Form factor from on-shell W+Z -> W+A scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_WZWA =", ffscale_wzwa

      ! AA -> AA
      ffscale_aa = searchff_maxpwave_ind(real_pwave_aa, "AA   -> AA  ", 0d0/180d0*pi, unitviol_aa)
      write(11,*) "! Form factor from on-shell AA -> AA scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_AAAA =", ffscale_aa

      ! ZA -> AA
      ffscale_zaaa = searchff_maxpwave_ind(real_pwave_zaaa, "ZA   -> AA  ", 0d0/180d0*pi, unitviol_zaaa)
      write(11,*) "! Form factor from on-shell ZA -> AA scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_ZAAA =", ffscale_zaaa

      ! ZA -> ZA
      ffscale_zaza = searchff_maxpwave_ind(real_pwave_zaza, "ZA   -> ZA  ", 0d0/180d0*pi, unitviol_zaza)
      write(11,*) "! Form factor from on-shell ZA -> ZA scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_ZAZA =", ffscale_zaza

      ! ZZ -> AA
      ffscale_zzaa = searchff_maxpwave_ind(real_pwave_zzaa, "ZZ   -> AA  ", 0d0/180d0*pi, unitviol_zzaa)
      write(11,*) "! Form factor from on-shell ZZ -> AA scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_ZZAA =", ffscale_zzaa

      ! ZZ -> ZA
      ffscale_zzza = searchff_maxpwave_ind(real_pwave_zzza, "ZZ   -> ZA  ", 0d0/180d0*pi, unitviol_zzza)
      write(11,*) "! Form factor from on-shell ZZ -> ZA scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_ZZZA =", ffscale_zzza



      write(11,*) ""



      ! due to singularity in t-channel photon exchange integration range 0.1° < theta < 180°
      ffscale_ww_diag = searchff_diag_ind(real_pwave_ww, im_pwave_ww, "W+W- -> W+W-", 0.1d0/180d0*pi, unitviol_ww_diag)
      write(11,*) "! Form factor from on-shell W+W- -> W+W- scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_WWWW_diag =", ffscale_ww_diag

      ! due to singularity in t-channel photon exchange integration range 0.1° < theta < 180°
      ffscale_wpwp_diag = searchff_diag_ind(real_pwave_wpwp, im_pwave_wpwp, "W+W+ -> W+W+", 0.1d0/180d0*pi, unitviol_wpwp_diag)
      write(11,*) "! Form factor from on-shell W+W+ -> W+W+ scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_WPWP_diag =", ffscale_wpwp_diag

      ffscale_zz_diag = searchff_diag_ind(real_pwave_zz, im_pwave_zz, "ZZ   -> ZZ  ", 0d0/180d0*pi, unitviol_zz_diag)
      write(11,*) "! Form factor from on-shell ZZ -> ZZ scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_ZZZZ_diag =", ffscale_zz_diag

      ffscale_wz_diag = searchff_diag_ind(real_pwave_wz, im_pwave_wz, "W+Z  -> W+Z ", 0d0/180d0*pi, unitviol_wz_diag)
      write(11,*) "! Form factor from on-shell W+Z -> W+Z scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_WZWZ_diag =", ffscale_wz_diag


      ! Channels with external photons:

      ffscale_aa_diag = searchff_diag_ind(real_pwave_aa, im_pwave_aa, "AA   -> AA  ", 0d0/180d0*pi, unitviol_aa_diag)
      write(11,*) "! Form factor from on-shell AA -> AA scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_AAAA_diag =", ffscale_aa_diag

      ffscale_za_diag = searchff_diag_ind(real_pwave_zaza, im_pwave_zaza, "ZA   -> ZA  ", 0d0/180d0*pi, unitviol_za_diag)
      write(11,*) "! Form factor from on-shell ZA -> ZA scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_ZAZA_diag =", ffscale_za_diag

      ffscale_wa_diag = searchff_diag_ind(real_pwave_wa, im_pwave_wa, "W+A  -> W+A ", 0d0/180d0*pi, unitviol_wa_diag)
      write(11,*) "! Form factor from on-shell W+A -> W+A scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_WAWA_diag =", ffscale_wa_diag



      write(11,*) ""



      ! due to singularity in t-channel photon exchange integration range 0.1° < theta < 180°
      ! ffscale_diag_q2 = searchff_diag_ind(real_pwave_wpwp, im_pwave_wpwp, "VV->VV, Q=2 ", 0.1d0/180d0*pi, unitviol_diag_q2)
      ffscale_diag_q2  = ffscale_wpwp_diag
      unitviol_diag_q2 = unitviol_wpwp_diag
      write(11,*) "! Form factor from Q = 2 on-shell VV scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_VVVV_Q_2 =", ffscale_diag_q2

      ! ffscale_diag_q1 = searchff_diag_ind(real_pwave_wz, im_pwave_wz, "VV->VV, Q=1 ", 0d0/180d0*pi, unitviol_diag_q1)
      ffscale_diag_q1  = ffscale_wz_diag
      unitviol_diag_q1 = unitviol_wz_diag
      write(11,*) "! Form factor from Q = 1 on-shell VV scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_VVVV_Q_1 =", ffscale_diag_q1

      ffscale_diag_q0 = searchff_diag_Q0(unitviol_diag_q0)
      write(11,*) "! Form factor from Q = 0 on-shell VV scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_VVVV_Q_0 =", ffscale_diag_q0


      ! Channels with external photons:

      ffscale_diag_q1_a = searchff_diag_Q1_a(unitviol_diag_q1_a)
      write(11,*) "! Form factor from Q = 1 on-shell VV scattering in GeV (with photons), diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_VVVV_Q_1_A =", ffscale_diag_q1_a

      ffscale_diag_q0_a = searchff_diag_Q0_a(unitviol_diag_q0_a)
      write(11,*) "! Form factor from Q = 0 on-shell VV scattering in GeV (with photons), diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_VVVV_Q_0_A =", ffscale_diag_q0_a



      write(11,*) ""

      write(11,"(A,F9.1,A)") "(FF = 0 => no FF needed for this channel)"

      close(11)

      print*, ""
      print*, ""
      print*, ""
      print*, "****************************************************************************************"
      print*, "****************************************************************************************"
      print*, ""
      print*, "VBFNLO tool for form factor calculation version 1.3.0"
      print*, ""
      print*, "Summary:"
      print*, ""
      print*, "Tested tree-level unitarity behaviour for zeroth partial wave in VV scattering."
      print*, "Unitarity criteria: The absolute value of the real part of the zeroth partial wave"
      print*, "has to be < 0.5 ." 
      print*, "Information on the different calculation modes can be found in the README."
      print*, ""
      print*, "The form factor scales determined assume the following form factor shape:"
      print*, ""
      print*, "             F = 1 / ( 1 + s / FFscale^2 )^FFexp"
      print*, ""
      write(*,"(A,F4.1)") " with FFexp =", ffexponent
      print*, ""
      write(*,"(A,F9.1,A)") " Tested center-of-mass energy up to: ", sqrts/1000d0, " TeV"
      print*, ""
      print*, "The following anomalous couplings to / Lambda^4 have been set != 0 :"
      print*, ""
      if (fs0_0.ne.0d0) write(*,"(A,E12.3,A)") " f_S0 =", fs0_0, " GeV^-4"
      if (fs1_0.ne.0d0) write(*,"(A,E12.3,A)") " f_S1 =", fs1_0, " GeV^-4"
      if (fM0_0.ne.0d0) write(*,"(A,E12.3,A)") " f_M0 =", fm0_0, " GeV^-4"
      if (fM1_0.ne.0d0) write(*,"(A,E12.3,A)") " f_M1 =", fm1_0, " GeV^-4"
      if (fM2_0.ne.0d0) write(*,"(A,E12.3,A)") " f_M2 =", fm2_0, " GeV^-4"
      if (fM3_0.ne.0d0) write(*,"(A,E12.3,A)") " f_M3 =", fm3_0, " GeV^-4"
      if (fM4_0.ne.0d0) write(*,"(A,E12.3,A)") " f_M4 =", fm4_0, " GeV^-4"
      if (fM5_0.ne.0d0) write(*,"(A,E12.3,A)") " f_M5 =", fm5_0, " GeV^-4"
      if (fM6_0.ne.0d0) write(*,"(A,E12.3,A)") " f_M6 =", fm6_0, " GeV^-4"
      if (fM7_0.ne.0d0) write(*,"(A,E12.3,A)") " f_M7 =", fm7_0, " GeV^-4"
      if (fT0_0.ne.0d0) write(*,"(A,E12.3,A)") " f_T0 =", ft0_0, " GeV^-4"
      if (fT1_0.ne.0d0) write(*,"(A,E12.3,A)") " f_T1 =", ft1_0, " GeV^-4"
      if (fT2_0.ne.0d0) write(*,"(A,E12.3,A)") " f_T2 =", ft2_0, " GeV^-4"
      if (fT5_0.ne.0d0) write(*,"(A,E12.3,A)") " f_T5 =", ft5_0, " GeV^-4"
      if (fT6_0.ne.0d0) write(*,"(A,E12.3,A)") " f_T6 =", ft6_0, " GeV^-4"
      if (fT7_0.ne.0d0) write(*,"(A,E12.3,A)") " f_T7 =", ft7_0, " GeV^-4"
      if (fT8_0.ne.0d0) write(*,"(A,E12.3,A)") " f_T8 =", ft8_0, " GeV^-4"
      if (fT9_0.ne.0d0) write(*,"(A,E12.3,A)") " f_T9 =", ft9_0, " GeV^-4"
      if (fb_0.ne.0d0) write(*,"(A,E12.3,A)") " f_B =", fb_0, " GeV^-2"
      if (fw_0.ne.0d0) write(*,"(A,E12.3,A)") " f_W =", fw_0, " GeV^-2"
      if (fwww_0.ne.0d0) write(*,"(A,E12.3,A)") " f_WWW =", fwww_0, " GeV^-2"
      if (fww_0.ne.0d0) write(*,"(A,E12.3,A)") " f_WW =", fWw_0, " GeV^-2"
      if (fphid_0.ne.0d0) write(*,"(A,E12.3,A)") " f_PHID =", fphid_0, " GeV^-2"
      if (fbtilde_0.ne.0d0) write(*,"(A,E12.3,A)") " f_BT =", fbtilde_0, " GeV^-2"
      if (fwwwtilde_0.ne.0d0) write(*,"(A,E12.3,A)") " f_WWWT =", fwwwtilde_0, " GeV^-2"
      if (fwwtilde_0.ne.0d0) write(*,"(A,E12.3,A)") " f_WWT =", fWwtilde_0, " GeV^-2"
      if (fdwtilde_0.ne.0d0) write(*,"(A,E12.3,A)") " f_DWT =", fdwtilde_0, " GeV^-2"

      print*, ""
      print*, ""
      print*, ""
      print*, "Results for each channel, taking only the helicity combination with the largest "
      print*, "contribution to the zeroth partial wave into account:"
      print*, ""

      if (ffscale_ww.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WWWW =", ffscale_ww, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_ww/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+W- -> W+W- scattering found."
      endif
      if (ffscale_wpwp.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WPWP =", ffscale_wpwp, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_wpwp/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+W+ -> W+W+ scattering found."
      endif
      if (ffscale_wwzz.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WWZZ =", ffscale_wwzz, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_wwzz/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+W- -> ZZ scattering found."
      endif
      if (ffscale_wz.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WZWZ =", ffscale_wz, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_wz/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+Z -> W+Z scattering found."
      endif
      if (ffscale_zz.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_ZZZZ =", ffscale_zz, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_zz/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in ZZ -> ZZ scattering found."
      endif

      print*, ""
      print*, "Channels with external photons:"
      if (ffscale_wwaa.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WWAA =", ffscale_wwaa, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_wwaa/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+W- -> AA scattering found."
      endif
      if (ffscale_wwza.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WWZA =", ffscale_wwza, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_wwza/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+W- -> ZA scattering found."
      endif
      if (ffscale_wa.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WAWA =", ffscale_wa, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_wa/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+A -> W+A scattering found."
      endif
      if (ffscale_wzwa.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WZWA =", ffscale_wzwa, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_wzwa/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+Z -> W+A scattering found."
      endif
      if (ffscale_aa.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_AAAA =", ffscale_aa, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_aa/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in AA -> AA scattering found."
      endif
      if (ffscale_zaaa.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_ZAAA =", ffscale_zaaa, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_zaaa/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in ZA -> AA scattering found."
      endif
      if (ffscale_zaza.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_ZAZA =", ffscale_zaza, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_zaza/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in ZA -> ZA scattering found."
      endif
      if (ffscale_zzaa.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_ZZAA =", ffscale_zzaa, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_zzaa/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in ZZ -> AA scattering found."
      endif
      if (ffscale_zzza.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_ZZZA =", ffscale_zzza, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_zzza/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in ZZ -> ZA scattering found."
      endif

      print*, ""
      print*, ""
      print*, ""
      print*, "Results for each channel, taking contributions from all helicity combinations to"
      print*, "the zeroth partial wave into account by diagonalizing the T-matrix:"
      print*, ""

      if (ffscale_ww_diag.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WWWW_diag =", ffscale_ww_diag, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_ww_diag/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+W- -> W+W- scattering found."
      endif
      if (ffscale_wpwp_diag.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WPWP_diag =", ffscale_wpwp_diag, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_wpwp_diag/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+W+ -> W+W+ scattering found."
      endif
      if (ffscale_zz_diag.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_ZZZZ_diag =", ffscale_zz_diag, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_zz_diag/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in ZZ -> ZZ scattering found."
      endif
      if (ffscale_wz_diag.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WZWZ_diag =", ffscale_wz_diag, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_wz_diag/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+Z -> W+Z scattering found."
      endif

      print*, ""
      print*, "Channels with external photons:"
      if (ffscale_aa_diag.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_AAAA_diag =", ffscale_aa_diag, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_aa_diag/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in AA -> AA scattering found."
      endif
      if (ffscale_za_diag.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_ZAZA_diag =", ffscale_za_diag, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_za_diag/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in ZA -> ZA scattering found."
      endif
      if (ffscale_wa_diag.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WAWA_diag =", ffscale_wa_diag, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_wa_diag/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+A -> W+A scattering found."
      endif

      print*, ""
      print*, ""
      print*, ""
      print*, "Results taking contributions from all channels involving W & Z bosons with all helicity "
      print*, "combinations to the zeroth partial wave into account by diagonalizing the T-matrix:"
      print*, ""
      print*, "Channels with electric charge Q = 2: W+W+ -> W+W+"
      print*, "Channels with electric charge Q = 1: W+Z  -> W+Z"
      print*, "Channels with electric charge Q = 0: W+W- -> W+W-, W+W- -> ZZ, ZZ -> W+W-, ZZ -> ZZ"
      print*, ""

      if (ffscale_diag_q2.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_VVVV_Q_2 =", ffscale_diag_q2, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_diag_q2/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in Q = 2 channels found."
      endif
      if (ffscale_diag_q1.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_VVVV_Q_1 =", ffscale_diag_q1, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_diag_q1/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in Q = 1 channels found."
      endif
      if (ffscale_diag_q0.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_VVVV_Q_0 =", ffscale_diag_q0, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_diag_q0/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in Q = 0 channels found."
      endif

      print*, ""
      print*, "Including external photons:"
      print*, "Channels with electric charge Q = 1: {W+Z, W+A -> {W+Z, W+A}"
      print*, "Channels with electric charge Q = 0: {W+W-, ZZ, ZA, AA} -> {W+W-, ZZ, ZA, AA}"
      print*, ""

      if (ffscale_diag_q1_a.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_VVVV_Q_1_A =", ffscale_diag_q1_a, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_diag_q1_a/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in Q = 1 channels including photons found."
      endif
      if (ffscale_diag_q0_a.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_VVVV_Q_0_A =", ffscale_diag_q0_a, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_diag_q0_a/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in Q = 0 channels including photons found."
      endif

      print*, ""
      print*, ""
      print*, "****************************************************************************************"
      print*, "****************************************************************************************"
      print*, ""
      print*, ""

!      if (fm2_0.gt.0d0 .or. fm3_0.gt.0d0) then
!        print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!        print*, ""
!        print*, "Remark:"
!        print*, "The operators M2 and M3 do have only small contributions to the zeroth partial wave"
!        print*, "for WZ->WZ scattering. Consider using the form factor scale from WW->ZZ scattering!"
!        print*, ""
!        print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!        print*, ""
!        print*, ""
!      endif

      end 


