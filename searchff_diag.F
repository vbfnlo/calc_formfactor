      real*8 function searchff_diag_Q1_A(unitviol)
      implicit none

      real*8 real_pwave_wz, real_pwave_wa, real_pwave_wawz, real_pwave_wzwa
      external real_pwave_wz, real_pwave_wa, real_pwave_wawz, real_pwave_wzwa
      real*8 im_pwave_wz, im_pwave_wa, im_pwave_wawz, im_pwave_wzwa
      external im_pwave_wz, im_pwave_wa, im_pwave_wawz, im_pwave_wzwa

      include "coupl.inc"
      include "global.inc"
      include "an_couplings_4v.inc"

      real*8 real_pwave_value, max_pwave_value, real_pwave_error, minangle
      real*8 im_pwave_value, im_pwave_error, unitviol

      integer mat_i, mat_j, i, j, ier
      complex*16 mat(0:4,0:4)
      double precision mat_acc

      real*8 ff_lowerlimit, ff_upperlimit, ffmassscale_old
      real*8 energy

      ! transfer center-of-mass energy, partial wave and helicity information
      ! to amplitude subroutine
      real*8 sc
      common /cms/ sc
      integer pw
      INTEGER NHEL(4)
      common /partwave/ pw, nhel

      select case (outputdetails)
      case (0)
            print*, ""
            print*, "Checking tree-level unitarity violation with on-shell VV->VV scattering"
            print*, "including all Q=1 channels involving W / Z / photon using the largest"
            print*, "helicity combination of the zeroth partial wave..."
            print*, ""
      case (1)
            print*, ""
            print*, ""
            print*, "*********************************************************************"
            print*, "Checking tree-level unitarity violation with on-shell VV scattering,"
            print*, "including channels W+Z -> W+Z, W+Z -> W+A, W+A -> W+Z and W+A -> W+A:"
            print*, "*********************************************************************"
            print*, ""
            print*, "Unitarity criteria: The largest real part of the eigenvalues of the"
            print*, "T-matrix for the zeroth partial wave has to be < 0.5 ."
      case default
            print*, ""
            print*, ""
            print*, ""
            print*, ""
            print*, ""
            print*, "*********************************************************************"
            print*, ""
            print*, "Checking tree-level unitarity violation with on-shell VV scattering,"
            print*, "including channels W+Z -> W+Z, W+Z -> W+A, W+A -> W+Z and W+A -> W+A:"
            print*, ""
            print*, "*********************************************************************"
            print*, ""
            print*, "Testing tree-level unitarity behaviour for zeroth partial wave in"
            print*, "Q=1 VV scattering, diagonalizing the T-matrix."
            print*, "Unitarity criteria: The largest real part of the eigenvalues for"
            print*, "the zeroth partial wave has to be < 0.5 ." 
            print*, ""
      end select

      ! first check if form factor is neccessary:
      formfact = .false.
      energy = 100d0
      max_pwave_value = 0d0
      do while (energy .lt. sqrts .and. abs(max_pwave_value).le.0.5d0)
        max_pwave_value = 0d0
        energy = energy + 100d0
        sc = energy**2
        ! pwave0: only h1=h2 & h3=h4 allowed
        do mat_i=0,4
          do mat_j=0,4
             if (mat_i.eq.4) then
               nhel(1) = 1
             else
               nhel(1) = mod(mat_i,3)-1
             endif
             if (mat_j.eq.4) then
               nhel(3) = 1
             else
               nhel(3) = mod(mat_j,3)-1
             endif
             nhel(2) = nhel(1)
             nhel(4) = nhel(3)
             if (mat_i.lt.3 .and. mat_j.lt.3) then
               call gaussint(real_pwave_wz,0d0/180d0*pi, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_wz,0d0/180d0*pi, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.3 .and. mat_j.lt.3) then
               call gaussint(real_pwave_wawz,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_wawz,0d0, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.lt.3 .and. mat_j.ge.3) then
               call gaussint(real_pwave_wzwa,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_wzwa,0d0, pi, im_pwave_value,im_pwave_error)
             else
               call gaussint(real_pwave_wa,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_wa,0d0, pi, im_pwave_value,im_pwave_error)
             endif
             mat(mat_i,mat_j) = real_pwave_value + (0.,1.)*im_pwave_value
          enddo
        enddo
!        print*, ""
!        do i=0,4
!          print*, ((mat(i,j)), j=0,4)
!        enddo
!        print*, ""
        mat_acc = 1d-6
        call eigenvalues(mat,5,mat_acc,10000,ier)
        do i=0,4
          if (abs(dreal(mat(i,i))).gt.abs(max_pwave_value)) max_pwave_value = abs(dreal(mat(i,i)))
        enddo
      enddo

      if (outputdetails.ge.1) print*, ""
      if (abs(max_pwave_value).gt.0.5) then
        if (outputdetails.ge.1) write(*,"(A,A,F9.0,A)") " Tree-level unitarity violation without",
     &      " form factor at ", energy, " GeV"
        unitviol = energy
      else
        if (outputdetails.ge.1) write(*,"(A,A,F9.0,A)") " No tree-level unitarity violation without",
     &      " form factor up to ", energy, " GeV"
        unitviol = 0d0
      endif
      if (outputdetails.ge.2) then
        print*, ""
        print*, ""
      endif

      if (abs(max_pwave_value).gt.0.5) then
        if (outputdetails.ge.1) then
          print*, ""
          write(*,"(A,F9.0,A)") " Searching optimal Lambda_FF for energies up to  ", sqrts," GeV ..."
          print*, ""
        endif
        formfact = .true.
        ffmassscale   = -1d0
        ff_upperlimit = 2000d0
        ff_lowerlimit = 0d0
        ffmassscale_old   = 0d0
        do while ( abs(abs(ffmassscale/ffmassscale_old)-1d0).gt.0.001 .or. abs(max_pwave_value).gt.0.5 )
          if (ff_upperlimit.le.10d0) then
            print*, ""
            print*, "Error:"
            print*, "Exponent of form factor probably too weak, would result in Lambda_FF < 10 GeV."
            print*, "Exiting!"
            print*, ""
            stop
          endif
          if (ffmassscale.lt.0d0) then    !init
            ffmassscale   = ff_upperlimit
          else
            ffmassscale_old = ffmassscale
            if (abs(max_pwave_value).gt.0.5) then
              ff_upperlimit = ffmassscale
              ffmassscale   = (ffmassscale + ff_lowerlimit) / 2d0
            else
              ff_lowerlimit = ffmassscale
              if (ff_upperlimit.le.ffmassscale) then
                ff_upperlimit = 2d0*ffmassscale
                ffmassscale   = ff_upperlimit
              else
                ffmassscale   = (ff_upperlimit + ffmassscale) / 2d0
              endif
            endif
          endif
          max_pwave_value = 0d0
          energy = 100d0
          do while (energy .lt. sqrts .and. abs(max_pwave_value).le.0.5d0)
            max_pwave_value = 0d0
            energy = energy + 100d0
            sc = energy**2
            ! pwave0: only h1=h2 & h3=h4 allowed
            do mat_i=0,4
                do mat_j=0,4
                   if (mat_i.eq.4) then
                     nhel(1) = 1
                   else
                     nhel(1) = mod(mat_i,3)-1
                   endif
                   if (mat_j.eq.4) then
                     nhel(3) = 1
                   else
                     nhel(3) = mod(mat_j,3)-1
                   endif
                   nhel(2) = nhel(1)
                   nhel(4) = nhel(3)
                   if (mat_i.lt.3 .and. mat_j.lt.3) then
                     call gaussint(real_pwave_wz,0d0/180d0*pi, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_wz,0d0/180d0*pi, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.ge.3 .and. mat_j.lt.3) then
                     call gaussint(real_pwave_wawz,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_wawz,0d0, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.lt.3 .and. mat_j.ge.3) then
                     call gaussint(real_pwave_wzwa,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_wzwa,0d0, pi, im_pwave_value,im_pwave_error)
                   else
                     call gaussint(real_pwave_wa,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_wa,0d0, pi, im_pwave_value,im_pwave_error)
                   endif
                   mat(mat_i,mat_j) = real_pwave_value + (0.,1.)*im_pwave_value
                enddo
            enddo
            mat_acc = 1d-6
            call eigenvalues(mat,5,mat_acc,10000,ier)
            do i=0,4
              if (abs(dreal(mat(i,i))).gt.abs(max_pwave_value)) max_pwave_value = abs(dreal(mat(i,i)))
            enddo
          enddo
          if (ldebugout) write(*,"(A,F11.2,A,F9.4,A,F9.0,A)") " Lambda_FF: ", ffmassscale, "   ---   Re(pwave_0): ", 
     &      abs(max_pwave_value), "   ---   at ", energy, " GeV"
        enddo
        if (outputdetails.eq.2) then
              print*, ""
              print*, ""
              print*, ""
              print*, "************************************************************************************"
              print*, "* Result for VV scattering, combining Q=1 channels and diagonalizing the T-matrix: *"
              print*, "************************************************************************************"
              print*, ""
              write(*,"(A,F9.0,A)") " Lambda_FF =", ffmassscale, " GeV"
              print*, ""
              write(*,"(A)") " prevents tree-level unitarity violation for Q=1 on-shell VV scattering"
              write(*,"(A,F9.0,A)") " up to sqrt(s) =", sqrts," GeV."
              print*, ""
        endif
        if (outputdetails.eq.1) write(*,"(A,F9.0,A)") " Result: Lambda_FF =", ffmassscale, " GeV"
        searchff_diag_Q1_A = ffmassscale
      else
        if (outputdetails.eq.2) then
              print*, ""
              print*, "No form factor is needed for this set of parameters!"
              print*, "Tree-level unitarity is preserved for Q=1 on-shell VV scattering"
              write(*,"(A,F9.0,A)") " up to sqrt(s) =", sqrts," GeV."
              print*, ""
        endif
        if (outputdetails.eq.1) then
            print*, "No form factor is needed for this set of parameters, tree level unitarity is"
            print*, "preserved for the zeroth partial wave!"
        endif
        searchff_diag_Q1_A = 0d0
      endif

      return
      end 



c**************************************************************************************************************
c**************************************************************************************************************
c**************************************************************************************************************



      real*8 function searchff_diag_Q0(unitviol)
      implicit none

      real*8 real_pwave_ww, real_pwave_zz, real_pwave_wwzz, real_pwave_zzww
      external real_pwave_ww, real_pwave_zz, real_pwave_wwzz, real_pwave_zzww
      real*8 im_pwave_ww, im_pwave_zz, im_pwave_wwzz, im_pwave_zzww
      external im_pwave_ww, im_pwave_zz, im_pwave_wwzz, im_pwave_zzww

      include "coupl.inc"
      include "global.inc"
      include "an_couplings_4v.inc"

      real*8 real_pwave_value, max_pwave_value, real_pwave_error, minangle
      real*8 im_pwave_value, im_pwave_error, unitviol

      integer mat_i, mat_j, i, j, ier
      complex*16 mat(0:5,0:5)
      double precision mat_acc

      real*8 ff_lowerlimit, ff_upperlimit, ffmassscale_old
      real*8 energy

      ! transfer center-of-mass energy, partial wave and helicity information
      ! to amplitude subroutine
      real*8 sc
      common /cms/ sc
      integer pw
      INTEGER NHEL(4)
      common /partwave/ pw, nhel

      select case (outputdetails)
      case (0)
            print*, ""
            print*, "Checking tree-level unitarity violation with on-shell VV->VV scattering"
            print*, "including all Q=0 channels involving W and Z bosons using the largest"
            print*, "helicity combination of the zeroth partial wave..."
            print*, ""
      case (1)
            print*, ""
            print*, ""
            print*, "*********************************************************************"
            print*, "Checking tree-level unitarity violation with on-shell VV scattering,"
            print*, "including channels W+W- -> W+W-, W+W- -> ZZ, ZZ -> W+W- and ZZ -> ZZ:"
            print*, "*********************************************************************"
            print*, ""
            print*, "Unitarity criteria: The largest real part of the eigenvalues of the"
            print*, "T-matrix for the zeroth partial wave has to be < 0.5 ."
      case default
            print*, ""
            print*, ""
            print*, ""
            print*, ""
            print*, ""
            print*, "*********************************************************************"
            print*, ""
            print*, "Checking tree-level unitarity violation with on-shell VV scattering,"
            print*, "including channels W+W- -> W+W-, W+W- -> ZZ, ZZ -> W+W- and ZZ -> ZZ:"
            print*, ""
            print*, "*********************************************************************"
            print*, ""
            print*, "Testing tree-level unitarity behaviour for zeroth partial wave in"
            print*, "Q=0 VV scattering, diagonalizing the T-matrix."
            print*, "Unitarity criteria: The largest real part of the eigenvalues for"
            print*, "the zeroth partial wave has to be < 0.5 ." 
            print*, ""
      end select

      ! first check if form factor is neccessary:
      formfact = .false.
      energy = 100d0
      max_pwave_value = 0d0
      do while (energy .lt. sqrts .and. abs(max_pwave_value).le.0.5d0)
        max_pwave_value = 0d0
        energy = energy + 100d0
        sc = energy**2
        ! pwave0: only h1=h2 & h3=h4 allowed
        do mat_i=0,5
          do mat_j=0,5
             nhel(1) = mod(mat_i,3)-1
             nhel(2) = mod(mat_i,3)-1
             nhel(3) = mod(mat_j,3)-1
             nhel(4) = mod(mat_j,3)-1
             if (mat_i.lt.3 .and. mat_j.lt.3) then
               call gaussint(real_pwave_ww,0.1d0/180d0*pi, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_ww,0.1d0/180d0*pi, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.3 .and. mat_j.lt.3) then
               call gaussint(real_pwave_zzww,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_zzww,0d0, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.lt.3 .and. mat_j.ge.3) then
               call gaussint(real_pwave_wwzz,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_wwzz,0d0, pi, im_pwave_value,im_pwave_error)
             else
               call gaussint(real_pwave_zz,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_zz,0d0, pi, im_pwave_value,im_pwave_error)
             endif
             mat(mat_i,mat_j) = real_pwave_value + (0.,1.)*im_pwave_value
          enddo
        enddo
!        print*, ""
!        do i=0,5
!          print*, ((mat(i,j)), j=0,5)
!        enddo
!        print*, ""
        mat_acc = 1d-6
        call eigenvalues(mat,6,mat_acc,10000,ier)
        do i=0,5
          if (abs(dreal(mat(i,i))).gt.abs(max_pwave_value)) max_pwave_value = abs(dreal(mat(i,i)))
        enddo
      enddo

      if (outputdetails.ge.1) print*, ""
      if (abs(max_pwave_value).gt.0.5) then
        if (outputdetails.ge.1) write(*,"(A,A,F9.0,A)") " Tree-level unitarity violation without",
     &      " form factor at ", energy, " GeV"
        unitviol = energy
      else
        if (outputdetails.ge.1) write(*,"(A,A,F9.0,A)") " No tree-level unitarity violation without",
     &      " form factor up to ", energy, " GeV"
        unitviol = 0d0
      endif
      if (outputdetails.ge.2) then
        print*, ""
        print*, ""
      endif

      if (abs(max_pwave_value).gt.0.5) then
        if (outputdetails.ge.1) then
          print*, ""
          write(*,"(A,F9.0,A)") " Searching optimal Lambda_FF for energies up to  ", sqrts," GeV ..."
          print*, ""
        endif
        formfact = .true.
        ffmassscale   = -1d0
        ff_upperlimit = 2000d0
        ff_lowerlimit = 0d0
        ffmassscale_old   = 0d0
        do while ( abs(abs(ffmassscale/ffmassscale_old)-1d0).gt.0.001 .or. abs(max_pwave_value).gt.0.5 )
          if (ff_upperlimit.le.10d0) then
            print*, ""
            print*, "Error:"
            print*, "Exponent of form factor probably too weak, would result in Lambda_FF < 10 GeV."
            print*, "Exiting!"
            print*, ""
            stop
          endif
          if (ffmassscale.lt.0d0) then    !init
            ffmassscale   = ff_upperlimit
          else
            ffmassscale_old = ffmassscale
            if (abs(max_pwave_value).gt.0.5) then
              ff_upperlimit = ffmassscale
              ffmassscale   = (ffmassscale + ff_lowerlimit) / 2d0
            else
              ff_lowerlimit = ffmassscale
              if (ff_upperlimit.le.ffmassscale) then
                ff_upperlimit = 2d0*ffmassscale
                ffmassscale   = ff_upperlimit
              else
                ffmassscale   = (ff_upperlimit + ffmassscale) / 2d0
              endif
            endif
          endif
          max_pwave_value = 0d0
          energy = 100d0
          do while (energy .lt. sqrts .and. abs(max_pwave_value).le.0.5d0)
            max_pwave_value = 0d0
            energy = energy + 100d0
            sc = energy**2
            ! pwave0: only h1=h2 & h3=h4 allowed
            do mat_i=0,5
              do mat_j=0,5
                 nhel(1) = mod(mat_i,3)-1
                 nhel(2) = mod(mat_i,3)-1
                 nhel(3) = mod(mat_j,3)-1
                 nhel(4) = mod(mat_j,3)-1
                 if (mat_i.lt.3 .and. mat_j.lt.3) then
                   call gaussint(real_pwave_ww,0.1d0/180d0*pi, pi, real_pwave_value,real_pwave_error)
                   call gaussint(im_pwave_ww,0.1d0/180d0*pi, pi, im_pwave_value,im_pwave_error)
                 elseif (mat_i.ge.3 .and. mat_j.lt.3) then
                   call gaussint(real_pwave_zzww,0d0, pi, real_pwave_value,real_pwave_error)
                   call gaussint(im_pwave_zzww,0d0, pi, im_pwave_value,im_pwave_error)
                 elseif (mat_i.lt.3 .and. mat_j.ge.3) then
                   call gaussint(real_pwave_wwzz,0d0, pi, real_pwave_value,real_pwave_error)
                   call gaussint(im_pwave_wwzz,0d0, pi, im_pwave_value,im_pwave_error)
                 else
                   call gaussint(real_pwave_zz,0d0, pi, real_pwave_value,real_pwave_error)
                   call gaussint(im_pwave_zz,0d0, pi, im_pwave_value,im_pwave_error)
                 endif
                 mat(mat_i,mat_j) = real_pwave_value + (0.,1.)*im_pwave_value
              enddo
            enddo
            mat_acc = 1d-6
            call eigenvalues(mat,6,mat_acc,10000,ier)
            do i=0,5
              if (abs(dreal(mat(i,i))).gt.abs(max_pwave_value)) max_pwave_value = abs(dreal(mat(i,i)))
            enddo
          enddo
          if (ldebugout) write(*,"(A,F11.2,A,F9.4,A,F9.0,A)") " Lambda_FF: ", ffmassscale, "   ---   Re(pwave_0): ", 
     &      abs(max_pwave_value), "   ---   at ", energy, " GeV"
        enddo
        if (outputdetails.eq.2) then
              print*, ""
              print*, ""
              print*, ""
              print*, "************************************************************************************"
              print*, "* Result for VV scattering, combining Q=0 channels and diagonalizing the T-matrix: *"
              print*, "************************************************************************************"
              print*, ""
              write(*,"(A,F9.0,A)") " Lambda_FF =", ffmassscale, " GeV"
              print*, ""
              write(*,"(A)") " prevents tree-level unitarity violation for Q=0 on-shell VV scattering"
              write(*,"(A,F9.0,A)") " up to sqrt(s) =", sqrts," GeV."
              print*, ""
        endif
        if (outputdetails.eq.1) write(*,"(A,F9.0,A)") " Result: Lambda_FF =", ffmassscale, " GeV"
        searchff_diag_Q0 = ffmassscale
      else
        if (outputdetails.eq.2) then
              print*, ""
              print*, "No form factor is needed for this set of parameters!"
              print*, "Tree-level unitarity is preserved for Q=0 on-shell VV scattering"
              write(*,"(A,F9.0,A)") " up to sqrt(s) =", sqrts," GeV."
              print*, ""
        endif
        if (outputdetails.eq.1) then
            print*, "No form factor is needed for this set of parameters, tree level unitarity is"
            print*, "preserved for the zeroth partial wave!"
        endif
        searchff_diag_Q0 = 0d0
      endif

      return
      end 
c**************************************************************************************************************
c**************************************************************************************************************
c**************************************************************************************************************



      real*8 function searchff_diag_Q0_a(unitviol)
      implicit none

      real*8 real_pwave_ww, real_pwave_zz, real_pwave_wwzz, real_pwave_zzww
      external real_pwave_ww, real_pwave_zz, real_pwave_wwzz, real_pwave_zzww
      real*8 real_pwave_aa, real_pwave_zzaa, real_pwave_wwaa, real_pwave_wwza, real_pwave_zaaa
      external real_pwave_aa, real_pwave_zzaa, real_pwave_wwaa, real_pwave_wwza, real_pwave_zaaa
      real*8 real_pwave_aaww, real_pwave_zaww, real_pwave_zaza, real_pwave_aaza, real_pwave_aazz, real_pwave_zazz, real_pwave_zzza
      external real_pwave_aaww, real_pwave_zaww, real_pwave_zaza, real_pwave_aaza, real_pwave_aazz,real_pwave_zazz,real_pwave_zzza

      real*8 im_pwave_ww, im_pwave_zz, im_pwave_wwzz, im_pwave_zzww
      external im_pwave_ww, im_pwave_zz, im_pwave_wwzz, im_pwave_zzww
      real*8 im_pwave_aa, im_pwave_zzaa, im_pwave_wwaa, im_pwave_wwza, im_pwave_zaaa
      external im_pwave_aa, im_pwave_zzaa, im_pwave_wwaa, im_pwave_wwza, im_pwave_zaaa
      real*8 im_pwave_aaww, im_pwave_zaww, im_pwave_zaza, im_pwave_aaza, im_pwave_aazz, im_pwave_zazz, im_pwave_zzza
      external im_pwave_aaww, im_pwave_zaww, im_pwave_zaza, im_pwave_aaza, im_pwave_aazz,im_pwave_zazz,im_pwave_zzza

      include "coupl.inc"
      include "global.inc"
      include "an_couplings_4v.inc"

      real*8 real_pwave_value, max_pwave_value, real_pwave_error, minangle
      real*8 im_pwave_value, im_pwave_error, unitviol

      integer mat_i, mat_j, i, j, ier
      complex*16 mat(0:9,0:9)
      double precision mat_acc

      real*8 ff_lowerlimit, ff_upperlimit, ffmassscale_old
      real*8 energy

      ! transfer center-of-mass energy, partial wave and helicity information
      ! to amplitude subroutine
      real*8 sc
      common /cms/ sc
      integer pw
      INTEGER NHEL(4)
      common /partwave/ pw, nhel

      select case (outputdetails)
      case (0)
            print*, ""
            print*, "Checking tree-level unitarity violation with on-shell VV->VV scattering"
            print*, "including all Q=0 channels involving W / Z / photon using the largest"
            print*, "helicity combination of the zeroth partial wave..."
            print*, ""
      case (1)
            print*, ""
            print*, ""
            print*, "*********************************************************************"
            print*, "Checking tree-level unitarity violation with on-shell VV scattering,"
            print*, "including channels {W+W-, ZZ, ZA, AA} -> {W+W-, ZZ, ZA, AA}:"
            print*, "*********************************************************************"
            print*, ""
            print*, "Unitarity criteria: The largest real part of the eigenvalues of the"
            print*, "T-matrix for the zeroth partial wave has to be < 0.5 ."
      case default
            print*, ""
            print*, ""
            print*, ""
            print*, ""
            print*, ""
            print*, "*********************************************************************"
            print*, ""
            print*, "Checking tree-level unitarity violation with on-shell VV scattering,"
            print*, "including channels {W+W-, ZZ, ZA, AA} -> {W+W-, ZZ, ZA, AA}:"
            print*, ""
            print*, "*********************************************************************"
            print*, ""
            print*, "Testing tree-level unitarity behaviour for zeroth partial wave in"
            print*, "Q=0 VV scattering, diagonalizing the T-matrix."
            print*, "Unitarity criteria: The largest real part of the eigenvalues for"
            print*, "the zeroth partial wave has to be < 0.5 ." 
            print*, ""
      end select

      ! first check if form factor is neccessary:
      formfact = .false.
      energy = 100d0
      max_pwave_value = 0d0
      do while (energy .lt. sqrts .and. abs(max_pwave_value).le.0.5d0)
        max_pwave_value = 0d0
        energy = energy + 100d0
        sc = energy**2
        ! pwave0: only h1=h2 & h3=h4 allowed
        do mat_i=0,9
          do mat_j=0,9
             if (mat_i.le.5) then
               nhel(1) = mod(mat_i,3)-1
             elseif (mod(mat_i,2).eq.1) then
               nhel(1) = 1
             else
               nhel(1) = -1
             endif
             if (mat_j.le.5) then
               nhel(3) = mod(mat_j,3)-1
             elseif (mod(mat_j,2).eq.1) then
               nhel(3) = 1
             else
               nhel(3) = -1
             endif
             nhel(2) = nhel(1)
             nhel(4) = nhel(3)

             if (mat_i.lt.3 .and. mat_j.lt.3) then
               call gaussint(real_pwave_ww,0.1d0/180d0*pi, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_ww,0.1d0/180d0*pi, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.3 .and. mat_i.lt.6 .and. mat_j.lt.3) then
               call gaussint(real_pwave_zzww,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_zzww,0d0, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.6 .and. mat_i.lt.8 .and. mat_j.lt.3) then
               call gaussint(real_pwave_zaww,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_zaww,0d0, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.8 .and. mat_j.lt.3) then
               call gaussint(real_pwave_aaww,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_aaww,0d0, pi, im_pwave_value,im_pwave_error)

             elseif (mat_i.lt.3 .and. mat_j.ge.3 .and. mat_j.lt.6) then
               call gaussint(real_pwave_wwzz,0.1d0/180d0*pi, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_wwzz,0.1d0/180d0*pi, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.3 .and. mat_i.lt.6 .and. mat_j.ge.3 .and. mat_j.lt.6) then
               call gaussint(real_pwave_zz,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_zz,0d0, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.6 .and. mat_i.lt.8 .and. mat_j.ge.3 .and. mat_j.lt.6) then
               call gaussint(real_pwave_zazz,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_zazz,0d0, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.8 .and. mat_j.ge.3 .and. mat_j.lt.6) then
               call gaussint(real_pwave_aazz,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_aazz,0d0, pi, im_pwave_value,im_pwave_error)

             elseif (mat_i.lt.3 .and. mat_j.ge.6 .and. mat_j.lt.8) then
               call gaussint(real_pwave_wwza,0.1d0/180d0*pi, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_wwza,0.1d0/180d0*pi, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.3 .and. mat_i.lt.6 .and. mat_j.ge.6 .and. mat_j.lt.8) then
               call gaussint(real_pwave_zzza,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_zzza,0d0, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.6 .and. mat_i.lt.8 .and. mat_j.ge.6 .and. mat_j.lt.8) then
               call gaussint(real_pwave_zaza,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_zaza,0d0, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.8 .and. mat_j.ge.6 .and. mat_j.lt.8) then
               call gaussint(real_pwave_aaza,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_aaza,0d0, pi, im_pwave_value,im_pwave_error)

             elseif (mat_i.lt.3 .and. mat_j.ge.8) then
               call gaussint(real_pwave_wwaa,0.1d0/180d0*pi, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_wwaa,0.1d0/180d0*pi, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.3 .and. mat_i.lt.6 .and. mat_j.ge.8) then
               call gaussint(real_pwave_zzaa,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_zzaa,0d0, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.6 .and. mat_i.lt.8 .and. mat_j.ge.8) then
               call gaussint(real_pwave_zaaa,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_zaaa,0d0, pi, im_pwave_value,im_pwave_error)
             elseif (mat_i.ge.8 .and. mat_j.ge.8) then
               call gaussint(real_pwave_aa,0d0, pi, real_pwave_value,real_pwave_error)
               call gaussint(im_pwave_aa,0d0, pi, im_pwave_value,im_pwave_error)
             endif

             mat(mat_i,mat_j) = real_pwave_value + (0.,1.)*im_pwave_value
          enddo
        enddo
!        print*, ""
!        do i=0,9
!          print*, ((mat(i,j)), j=0,5)
!        enddo
!        print*, ""
        mat_acc = 1d-6
        call eigenvalues(mat,10,mat_acc,10000,ier)
        do i=0,9
          if (abs(dreal(mat(i,i))).gt.abs(max_pwave_value)) max_pwave_value = abs(dreal(mat(i,i)))
        enddo
      enddo

      if (outputdetails.ge.1) print*, ""
      if (abs(max_pwave_value).gt.0.5) then
        if (outputdetails.ge.1) write(*,"(A,A,F9.0,A)") " Tree-level unitarity violation without",
     &      " form factor at ", energy, " GeV"
        unitviol = energy
      else
        if (outputdetails.ge.1) write(*,"(A,A,F9.0,A)") " No tree-level unitarity violation without",
     &      " form factor up to ", energy, " GeV"
        unitviol = 0d0
      endif
      if (outputdetails.ge.2) then
        print*, ""
        print*, ""
      endif

      if (abs(max_pwave_value).gt.0.5) then
        if (outputdetails.ge.1) then
          print*, ""
          write(*,"(A,F9.0,A)") " Searching optimal Lambda_FF for energies up to  ", sqrts," GeV ..."
          print*, ""
        endif
        formfact = .true.
        ffmassscale   = -1d0
        ff_upperlimit = 2000d0
        ff_lowerlimit = 0d0
        ffmassscale_old   = 0d0
        do while ( abs(abs(ffmassscale/ffmassscale_old)-1d0).gt.0.001 .or. abs(max_pwave_value).gt.0.5 )
          if (ff_upperlimit.le.10d0) then
            print*, ""
            print*, "Error:"
            print*, "Exponent of form factor probably too weak, would result in Lambda_FF < 10 GeV."
            print*, "Exiting!"
            print*, ""
            stop
          endif
          if (ffmassscale.lt.0d0) then    !init
            ffmassscale   = ff_upperlimit
          else
            ffmassscale_old = ffmassscale
            if (abs(max_pwave_value).gt.0.5) then
              ff_upperlimit = ffmassscale
              ffmassscale   = (ffmassscale + ff_lowerlimit) / 2d0
            else
              ff_lowerlimit = ffmassscale
              if (ff_upperlimit.le.ffmassscale) then
                ff_upperlimit = 2d0*ffmassscale
                ffmassscale   = ff_upperlimit
              else
                ffmassscale   = (ff_upperlimit + ffmassscale) / 2d0
              endif
            endif
          endif
          max_pwave_value = 0d0
          energy = 100d0
          do while (energy .lt. sqrts .and. abs(max_pwave_value).le.0.5d0)
            max_pwave_value = 0d0
            energy = energy + 100d0
            sc = energy**2
            ! pwave0: only h1=h2 & h3=h4 allowed
            do mat_i=0,9
                do mat_j=0,9
                   if (mat_i.le.5) then
                     nhel(1) = mod(mat_i,3)-1
                   elseif (mod(mat_i,2).eq.1) then
                     nhel(1) = 1
                   else
                     nhel(1) = -1
                   endif
                   if (mat_j.le.5) then
                     nhel(3) = mod(mat_j,3)-1
                   elseif (mod(mat_j,2).eq.1) then
                     nhel(3) = 1
                   else
                     nhel(3) = -1
                   endif
                   nhel(2) = nhel(1)
                   nhel(4) = nhel(3)

                   if (mat_i.lt.3 .and. mat_j.lt.3) then
                     call gaussint(real_pwave_ww,0.1d0/180d0*pi, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_ww,0.1d0/180d0*pi, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.ge.3 .and. mat_i.lt.6 .and. mat_j.lt.3) then
                     call gaussint(real_pwave_zzww,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_zzww,0d0, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.ge.6 .and. mat_i.lt.8 .and. mat_j.lt.3) then
                     call gaussint(real_pwave_zaww,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_zaww,0d0, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.ge.8 .and. mat_j.lt.3) then
                     call gaussint(real_pwave_aaww,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_aaww,0d0, pi, im_pwave_value,im_pwave_error)

                   elseif (mat_i.lt.3 .and. mat_j.ge.3 .and. mat_j.lt.6) then
                     call gaussint(real_pwave_wwzz,0.1d0/180d0*pi, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_wwzz,0.1d0/180d0*pi, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.ge.3 .and. mat_i.lt.6 .and. mat_j.ge.3 .and. mat_j.lt.6) then
                     call gaussint(real_pwave_zz,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_zz,0d0, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.ge.6 .and. mat_i.lt.8 .and. mat_j.ge.3 .and. mat_j.lt.6) then
                     call gaussint(real_pwave_zazz,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_zazz,0d0, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.ge.8 .and. mat_j.ge.3 .and. mat_j.lt.6) then
                     call gaussint(real_pwave_aazz,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_aazz,0d0, pi, im_pwave_value,im_pwave_error)

                   elseif (mat_i.lt.3 .and. mat_j.ge.6 .and. mat_j.lt.8) then
                     call gaussint(real_pwave_wwza,0.1d0/180d0*pi, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_wwza,0.1d0/180d0*pi, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.ge.3 .and. mat_i.lt.6 .and. mat_j.ge.6 .and. mat_j.lt.8) then
                     call gaussint(real_pwave_zzza,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_zzza,0d0, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.ge.6 .and. mat_i.lt.8 .and. mat_j.ge.6 .and. mat_j.lt.8) then
                     call gaussint(real_pwave_zaza,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_zaza,0d0, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.ge.8 .and. mat_j.ge.6 .and. mat_j.lt.8) then
                     call gaussint(real_pwave_aaza,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_aaza,0d0, pi, im_pwave_value,im_pwave_error)

                   elseif (mat_i.lt.3 .and. mat_j.ge.8) then
                     call gaussint(real_pwave_wwaa,0.1d0/180d0*pi, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_wwaa,0.1d0/180d0*pi, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.ge.3 .and. mat_i.lt.6 .and. mat_j.ge.8) then
                     call gaussint(real_pwave_zzaa,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_zzaa,0d0, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.ge.6 .and. mat_i.lt.8 .and. mat_j.ge.8) then
                     call gaussint(real_pwave_zaaa,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_zaaa,0d0, pi, im_pwave_value,im_pwave_error)
                   elseif (mat_i.ge.8 .and. mat_j.ge.8) then
                     call gaussint(real_pwave_aa,0d0, pi, real_pwave_value,real_pwave_error)
                     call gaussint(im_pwave_aa,0d0, pi, im_pwave_value,im_pwave_error)
                   endif

                   mat(mat_i,mat_j) = real_pwave_value + (0.,1.)*im_pwave_value
                enddo
            enddo
            mat_acc = 1d-6
            call eigenvalues(mat,10,mat_acc,10000,ier)
            do i=0,9
              if (abs(dreal(mat(i,i))).gt.abs(max_pwave_value)) max_pwave_value = abs(dreal(mat(i,i)))
            enddo
          enddo
          if (ldebugout) write(*,"(A,F11.2,A,F9.4,A,F9.0,A)") " Lambda_FF: ", ffmassscale, "   ---   Re(pwave_0): ", 
     &      abs(max_pwave_value), "   ---   at ", energy, " GeV"
        enddo
        if (outputdetails.eq.2) then
              print*, ""
              print*, ""
              print*, ""
              print*, "************************************************************************************"
              print*, "* Result for VV scattering, combining Q=0 channels and diagonalizing the T-matrix: *"
              print*, "************************************************************************************"
              print*, ""
              write(*,"(A,F9.0,A)") " Lambda_FF =", ffmassscale, " GeV"
              print*, ""
              write(*,"(A)") " prevents tree-level unitarity violation for Q=0 on-shell VV scattering"
              write(*,"(A,F9.0,A)") " up to sqrt(s) =", sqrts," GeV."
              print*, ""
        endif
        if (outputdetails.eq.1) write(*,"(A,F9.0,A)") " Result: Lambda_FF =", ffmassscale, " GeV"
        searchff_diag_Q0_A = ffmassscale
      else
        if (outputdetails.eq.2) then
              print*, ""
              print*, "No form factor is needed for this set of parameters!"
              print*, "Tree-level unitarity is preserved for Q=0 on-shell VV scattering"
              write(*,"(A,F9.0,A)") " up to sqrt(s) =", sqrts," GeV."
              print*, ""
        endif
        if (outputdetails.eq.1) then
            print*, "No form factor is needed for this set of parameters, tree level unitarity is"
            print*, "preserved for the zeroth partial wave!"
        endif
        searchff_diag_Q0_A = 0d0
      endif

      return
      end 

