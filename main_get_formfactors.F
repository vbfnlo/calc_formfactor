
      program get_formfactors

c computes form factors for anomalous quartic gauge couplings from vector boson scattering

      implicit none
      include "coupl.inc"
      include "an_couplings_4v.inc"

      real*8 pi
      parameter(pi =3.141592653589793d0)     
      integer k

      real*8 real_pwave_ww, real_pwave_wz, real_pwave_zz, real_pwave_wwzz, real_pwave_zzww
      external real_pwave_ww, real_pwave_wz, real_pwave_zz, real_pwave_wwzz, real_pwave_zzww
      real*8 im_pwave_ww, im_pwave_wz, im_pwave_zz, im_pwave_wwzz, im_pwave_zzww
      external im_pwave_ww, im_pwave_wz, im_pwave_zz, im_pwave_wwzz, im_pwave_zzww
      real*8 ffscale_ww, ffscale_wz, ffscale_zz, ffscale_wwzz, ffscale_zzww
      real*8 ffscale_ww_diag, ffscale_wz_diag, ffscale_zz_diag, ffscale_wwzz_diag, ffscale_zzww_diag
      real*8 ffscale_diag_q1, ffscale_diag_q0
      real*8 unitviol_ww, unitviol_wz, unitviol_zz, unitviol_wwzz, unitviol_zzww
      real*8 unitviol_ww_diag, unitviol_wz_diag, unitviol_zz_diag, unitviol_wwzz_diag, unitviol_zzww_diag
      real*8 unitviol_diag_q1, unitviol_diag_q0

      real*8 searchff_maxpwave_ind, searchff_diag_ind, searchff_diag_q0
      external searchff_maxpwave_ind, searchff_diag_ind, searchff_diag_q0

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

      ! due to singularity in t-channel photon exchange integration range 0.1° < theta < 180°
      ffscale_ww = searchff_maxpwave_ind(real_pwave_ww, "W+W- -> W+W-", 0.1d0/180d0*pi, unitviol_ww)
      write(11,*) "! Form factor from on-shell W+W- -> W+W- scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_WWWW =", ffscale_ww

      ! ZZ -> WW is the same
      ffscale_wwzz = searchff_maxpwave_ind(real_pwave_wwzz, "W+W- -> ZZ  ", 0d0/180d0*pi, unitviol_wwzz)
      write(11,*) "! Form factor from on-shell W+W- -> ZZ scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_WWZZ =", ffscale_wwzz

      ffscale_zz = searchff_maxpwave_ind(real_pwave_zz, "ZZ   -> ZZ  ", 0d0/180d0*pi, unitviol_zz)
      write(11,*) "! Form factor from on-shell ZZ -> ZZ scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_ZZZZ =", ffscale_zz

      ffscale_wz = searchff_maxpwave_ind(real_pwave_wz, "W+Z  -> W+Z ", 0d0/180d0*pi, unitviol_wz)
      write(11,*) "! Form factor from on-shell W+Z -> W+Z scattering in GeV:"
      write(11,"(A,F9.0,A)") " FFscale_WZWZ =", ffscale_wz

      write(11,*) ""

      ffscale_ww_diag = searchff_diag_ind(real_pwave_ww, im_pwave_ww, "W+W- -> W+W-", 0.1d0/180d0*pi, unitviol_ww_diag)
      write(11,*) "! Form factor from on-shell W+W- -> W+W- scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_WWWW_diag =", ffscale_ww_diag

      ffscale_zz_diag = searchff_diag_ind(real_pwave_zz, im_pwave_zz, "ZZ   -> ZZ  ", 0d0/180d0*pi, unitviol_zz_diag)
      write(11,*) "! Form factor from on-shell ZZ -> ZZ scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_ZZZZ_diag =", ffscale_zz_diag

      ffscale_wz_diag = searchff_diag_ind(real_pwave_wz, im_pwave_wz, "W+Z  -> W+Z ", 0d0/180d0*pi, unitviol_wz_diag)
      write(11,*) "! Form factor from on-shell W+Z -> W+Z scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_WZWZ_diag =", ffscale_wz_diag

      write(11,*) ""

      ffscale_diag_q1 = searchff_diag_ind(real_pwave_wz, im_pwave_wz, "VV->VV, Q=1 ", 0d0/180d0*pi, unitviol_diag_q1)
      write(11,*) "! Form factor from Q = 1 on-shell VV scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_VVVV_Q_1 =", ffscale_diag_q1

      ffscale_diag_q0 = searchff_diag_Q0(unitviol_diag_q0)
      write(11,*) "! Form factor from Q = 0 on-shell VV scattering in GeV, diagonalizing the T-matrix:"
      write(11,"(A,F9.0,A)") " FFscale_VVVV_Q_0 =", ffscale_diag_q0

      write(11,*) ""

      write(11,"(A,F9.1,A)") "(FF = 0 => no FF needed for this channel)"

      close(11)

      print*, ""
      print*, ""
      print*, ""
      print*, "****************************************************************************************"
      print*, "****************************************************************************************"
      print*, ""
      print*, "VBFNLO tool for form factor calculation version 1.0.0"
      print*, ""
      print*, "Summary:"
      print*, ""
      print*, "Tested tree-level unitarity behaviour for zeroth partial wave in VV scattering."
      print*, "Unitarity criteria: The absolute value of the real part of the zeroth partial wave"
      print*, "has to be < 0.5 ." 
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
      if (fs0_0.gt.0d0) write(*,"(A,E12.3,A)") " f_S0 =", fs0_0, " GeV^-4"
      if (fs1_0.gt.0d0) write(*,"(A,E12.3,A)") " f_S1 =", fs1_0, " GeV^-4"
      if (fM0_0.gt.0d0) write(*,"(A,E12.3,A)") " f_M0 =", fm0_0, " GeV^-4"
      if (fM1_0.gt.0d0) write(*,"(A,E12.3,A)") " f_M1 =", fm1_0, " GeV^-4"
      if (fM2_0.gt.0d0) write(*,"(A,E12.3,A)") " f_M2 =", fm2_0, " GeV^-4"
      if (fM3_0.gt.0d0) write(*,"(A,E12.3,A)") " f_M3 =", fm3_0, " GeV^-4"
      if (fM4_0.gt.0d0) write(*,"(A,E12.3,A)") " f_M4 =", fm4_0, " GeV^-4"
      if (fM5_0.gt.0d0) write(*,"(A,E12.3,A)") " f_M5 =", fm5_0, " GeV^-4"
      if (fM6_0.gt.0d0) write(*,"(A,E12.3,A)") " f_M6 =", fm6_0, " GeV^-4"
      if (fM7_0.gt.0d0) write(*,"(A,E12.3,A)") " f_M7 =", fm7_0, " GeV^-4"
      if (fT0_0.gt.0d0) write(*,"(A,E12.3,A)") " f_T0 =", ft0_0, " GeV^-4"
      if (fT1_0.gt.0d0) write(*,"(A,E12.3,A)") " f_T1 =", ft1_0, " GeV^-4"
      if (fT2_0.gt.0d0) write(*,"(A,E12.3,A)") " f_T2 =", ft2_0, " GeV^-4"
      if (fT5_0.gt.0d0) write(*,"(A,E12.3,A)") " f_T5 =", ft5_0, " GeV^-4"
      if (fT6_0.gt.0d0) write(*,"(A,E12.3,A)") " f_T6 =", ft6_0, " GeV^-4"
      if (fT7_0.gt.0d0) write(*,"(A,E12.3,A)") " f_T7 =", ft7_0, " GeV^-4"
      if (fT8_0.gt.0d0) write(*,"(A,E12.3,A)") " f_T8 =", ft8_0, " GeV^-4"
      if (fT9_0.gt.0d0) write(*,"(A,E12.3,A)") " f_T9 =", ft9_0, " GeV^-4"
      if (fb_0.gt.0d0) write(*,"(A,E12.3,A)") " f_B =", fb_0, " GeV^-2"
      if (fw_0.gt.0d0) write(*,"(A,E12.3,A)") " f_W =", fw_0, " GeV^-2"
      if (fwww_0.gt.0d0) write(*,"(A,E12.3,A)") " f_WWW =", fwww_0, " GeV^-2"
      if (fww_0.gt.0d0) write(*,"(A,E12.3,A)") " f_WW =", fWw_0, " GeV^-2"
      if (fbtilde_0.gt.0d0) write(*,"(A,E12.3,A)") " f_BT =", fbtilde_0, " GeV^-2"
      if (fwwwtilde_0.gt.0d0) write(*,"(A,E12.3,A)") " f_WWWT =", fwwwtilde_0, " GeV^-2"
      if (fwwtilde_0.gt.0d0) write(*,"(A,E12.3,A)") " f_WWT =", fWwtilde_0, " GeV^-2"
      if (fdwtilde_0.gt.0d0) write(*,"(A,E12.3,A)") " f_DWT =", fdwtilde_0, " GeV^-2"
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
      if (ffscale_wwzz.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WWZZ =", ffscale_wwzz, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_wwzz/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+W- -> ZZ scattering found."
      endif
      if (ffscale_zz.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_ZZZZ =", ffscale_zz, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_zz/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in ZZ -> ZZ scattering found."
      endif
      if (ffscale_wz.gt.0d0) then
        write(*,"(A,F9.0,A,A,F5.1,A)") " FFscale_WZWZ =", ffscale_wz, " GeV",
     &      "    ( without FF: |Re(pwave_0)| > 0.5 at ", unitviol_wz/1000d0, " TeV )"
      else
        print*, "No tree-level unitarity violation in W+Z -> W+Z scattering found."
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
      print*, ""
      print*, ""
      print*, "Results taking contributions from all channels involving W & Z bosons with all helicity "
      print*, "combinations to the zeroth partial wave into account by diagonalizing the T-matrix:"
      print*, ""
      print*, "Channels with electric charge Q = 1: W+Z  -> W+Z"
      print*, "Channels with electric charge Q = 0: W+W- -> W+W-, W+W- -> ZZ, ZZ -> W+W-, ZZ -> ZZ"
      print*, ""
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
      print*, ""
      print*, "****************************************************************************************"
      print*, "****************************************************************************************"
      print*, ""
      print*, ""

      if (fm2_0.gt.0d0 .or. fm3_0.gt.0d0) then
        print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print*, ""
        print*, "Remark:"
        print*, "The operators M2 and M3 do have only small contributions to the zeroth partial wave"
        print*, "for WZ->WZ scattering. Consider using the form factor scale from WW->ZZ scattering!"
        print*, ""
        print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print*, ""
        print*, ""
      endif

      end 


