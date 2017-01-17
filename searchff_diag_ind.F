
      real*8 function searchff_diag_ind(pwave_func_real, pwave_func_im, channel, minangle, unitviol)
      implicit none

      real*8 pwave_func_real, pwave_func_im
      external pwave_func_real, pwave_func_im

      character*12 channel

      include "coupl.inc"
      include "global.inc"
      include "an_couplings_4v.inc"

      real*8 real_pwave_value, max_pwave_value, real_pwave_error, minangle
      real*8 im_pwave_value, im_pwave_error, unitviol

      integer mat_i, mat_j, i, j, ier
      complex*16 mat(0:2,0:2)
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
            print*, "Checking tree-level unitarity violation with on-shell ", channel," scattering"
            print*, "diagonalizing the T-matrix of the zeroth partial wave..."
            print*, ""
      case (1)
            print*, ""
            print*, ""
            print*, "******************************************************************************"
            print*, "Checking tree-level unitarity violation with on-shell ", channel," scattering:"
            print*, "******************************************************************************"
            print*, ""
            print*, "Unitarity criteria: The largest real part of the eigenvalues of the"
            print*, "T-matrix for the zeroth partial wave has to be < 0.5 ."
      case default
            print*, ""
            print*, ""
            print*, ""
            print*, ""
            print*, ""
            print*, "******************************************************************************"
            print*, ""
            print*, "Checking tree-level unitarity violation with on-shell ", channel," scattering:"
            print*, ""
            print*, "******************************************************************************"
            print*, ""
            print*, "Testing tree-level unitarity behaviour for zeroth partial wave in"
            print*, channel, " scattering, diagonalizing the T-matrix."
            print*, "Unitarity criteria: The largest real part of the eigenvalues for"
            print*, "the zeroth partial wave has to be < 0.5 ." 
            print*, ""
      end select

      ! first check if form factor is neccessary:
      formfact = .false.
      energy = 100d0
      do while (energy .lt. sqrts .and. abs(max_pwave_value).le.0.5d0)
        max_pwave_value = 0d0
        energy = energy + 100d0
        sc = energy**2
        ! pwave0: only h1=h2 & h3=h4 allowed
        do mat_i=0,2
          do mat_j=0,2
             nhel(1) = mat_i-1
             nhel(2) = mat_i-1
             nhel(3) = mat_j-1
             nhel(4) = mat_j-1
             call gaussint(pwave_func_real,minangle, pi, real_pwave_value,real_pwave_error)
             call gaussint(pwave_func_im,minangle, pi, im_pwave_value,im_pwave_error)
             mat(mat_i,mat_j) = real_pwave_value + (0.,1.)*im_pwave_value
          enddo
        enddo
!        print*, ""
!        do i=0,2
!          print*, (mat(i,j), j=0,2)
!        enddo
!        print*, ""
        mat_acc = 1d-6
        call eigenvalues(mat,3,mat_acc,10000,ier)
        max_pwave_value = max(abs(dreal(mat(0,0))),abs(dreal(mat(1,1))),abs(dreal(mat(2,2))))
!        write(*,*) ier, mat_acc
!        write(*,*) (mat(i,i),i=0,2)
!        print*, ""
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
          do while (energy .lt. sqrts .and. max_pwave_value.le.0.5d0)
            max_pwave_value = 0d0
            energy = energy + 100d0
            sc = energy**2
            ! pwave0: only h1=h2 & h3=h4 allowed
            do mat_i=0,2
              do mat_j=0,2
                 nhel(1) = mat_i-1
                 nhel(2) = mat_i-1
                 nhel(3) = mat_j-1
                 nhel(4) = mat_j-1
                 call gaussint(pwave_func_real,minangle, pi, real_pwave_value,real_pwave_error)
                 call gaussint(pwave_func_im,minangle, pi, im_pwave_value,im_pwave_error)
                 mat(mat_i,mat_j) = real_pwave_value + (0.,1.)*im_pwave_value
              enddo
            enddo
            mat_acc = 1d-6
            call eigenvalues(mat,3,mat_acc,10000,ier)
            max_pwave_value = max(abs(dreal(mat(0,0))),abs(dreal(mat(1,1))),abs(dreal(mat(2,2))))
          enddo
          if (ldebugout) write(*,"(A,F11.2,A,F9.4,A,F9.0,A)") " Lambda_FF: ", ffmassscale, "   ---   Re(pwave_0): ", 
     &      abs(max_pwave_value), "   ---   at ", energy, " GeV"
        enddo
        if (outputdetails.eq.2) then
              print*, ""
              print*, ""
              print*, ""
              print*, "********************************************************************"
              print*, "* Result for ", channel, " scattering, diagonalizing the T-matrix: *"
              print*, "********************************************************************"
              print*, ""
              write(*,"(A,F9.0,A)") " Lambda_FF =", ffmassscale, " GeV"
              print*, ""
              write(*,"(A,A,A)") " prevents tree-level unitarity violation for on-shell ", channel," scattering"
              write(*,"(A,F9.0,A)") " up to sqrt(s) =", sqrts," GeV."
              print*, ""
        endif
        if (outputdetails.eq.1) write(*,"(A,F9.0,A)") " Result: Lambda_FF =", ffmassscale, " GeV"
        searchff_diag_ind = ffmassscale
      else
        if (outputdetails.eq.2) then
              print*, ""
              print*, "No form factor is needed for this set of parameters!"
              write(*,"(A,A,A)") " Tree-level unitarity is preserved for on-shell ", channel," scattering"
              write(*,"(A,F9.0,A)") " up to sqrt(s) =", sqrts," GeV."
              print*, ""
        endif        
        if (outputdetails.eq.1) then
            print*, "No form factor is needed for this set of parameters, tree level unitarity is"
            print*, "preserved for the zeroth partial wave!"
        endif
        searchff_diag_ind = 0d0
      endif

      return
      end 

