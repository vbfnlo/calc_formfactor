
      real*8 function searchff_maxpwave_ind(pwave_func_real, channel, minangle, unitviol)
      implicit none

      real*8 pwave_func_real
      external pwave_func_real

      character*12 channel

      include "coupl.inc"
      include "global.inc"
      include "an_couplings_4v.inc"

      real*8 real_pwave_value, max_pwave_value, real_pwave_error, minangle
      integer helcount
      logical photon_init, photon_final

      real*8 ff_lowerlimit, ff_upperlimit, ffmassscale_old
      real*8 energy, unitviol

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
            print*, "using the largest helicity combination of the zeroth partial wave..."
            print*, ""
      case (1)
            print*, ""
            print*, ""
            print*, "******************************************************************************"
            print*, "Checking tree-level unitarity violation with on-shell ", channel," scattering:"
            print*, "******************************************************************************"
            print*, ""
            print*, "Unitarity criteria: The real part of the contribution of the largest"
            print*, "helicity combination to the zeroth partial wave has to be < 0.5 ."
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
            print*, channel, " scattering individually for each helicity combination."
            print*, "Unitarity criteria: For each helicity combination the"
            print*, "absolute value of the real part of the zeroth partial wave"
            print*, "has to be < 0.5 ." 
            print*, ""
      end select
      
      if (channel.eq."AA   -> AA  " .or. channel.eq."ZA   -> AA  " .or.
     &    channel.eq."AA   -> ZA  " .or. channel.eq."ZA   -> ZA  " .or. channel.eq."W+A  -> W+A ") then
         photon_init  = .true.
         photon_final = .true.
         if (ldebugout) then
            print*, ""
            print*, "Photons in initial and final state: No longitudinal polarizations!"
            print*, ""
         endif
      elseif (channel.eq."ZZ   -> ZA  " .or. channel.eq."ZZ   -> AA  " .or. 
     &        channel.eq."W+W- -> AA  " .or. channel.eq."W+W- -> ZA  " .or. channel.eq."W+Z  -> W+A ") then
         photon_init  = .false.
         photon_final = .true.
         if (ldebugout) then
            print*, ""
            print*, "Photons in initial state: No longitudinal polarizations there!"
            print*, ""
         endif
      elseif (channel.eq."ZA   -> ZZ  " .or. channel.eq."AA   -> ZZ  " .or. 
     &        channel.eq."AA   -> W+W-" .or. channel.eq."ZA   -> W+W-" .or. channel.eq."W+A  -> W+Z ") then
         photon_init  = .true.
         photon_final = .false.
         if (ldebugout) then
            print*, ""
            print*, "Photons in final state: No longitudinal polarizations there!"
            print*, ""
         endif
      else
         photon_init  = .false.
         photon_final = .false.
      endif

      ! first check if form factor is neccessary:
      formfact = .false.
      max_pwave_value = 0d0
      energy = 100d0
      do while (energy .lt. sqrts .and. abs(max_pwave_value).le.0.5d0)
        max_pwave_value = 0d0
        energy = energy + 100d0
        sc = energy**2
        ! partial wave 0: only h1=h2 & h3=h4 allowed
        ! with photon in initial or final state: no longitudinal polarization
        do helcount = 0,8
           nhel(1) = mod(helcount/3,3)-1
           nhel(2) = mod(helcount/3,3)-1
           nhel(3) = mod(helcount,3)-1
           nhel(4) = mod(helcount,3)-1
           if (.not.((photon_init.and.(nhel(1).eq.0)) .or. (photon_final.and.(nhel(3).eq.0)))) then
              call gaussint(pwave_func_real,minangle, pi, real_pwave_value,real_pwave_error)
              if (abs(real_pwave_value) .gt. abs(max_pwave_value)) max_pwave_value = real_pwave_value
           endif
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
            ! partial wave 0: only h1=h2 & h3=h4 allowed
            ! with photon in initial or final state: no longitudinal polarization
            do helcount=0,8
               nhel(1) = mod(helcount/3,3)-1
               nhel(2) = mod(helcount/3,3)-1
               nhel(3) = mod(helcount,3)-1
               nhel(4) = mod(helcount,3)-1
               if (.not.((photon_init.and.(nhel(1).eq.0)) .or. (photon_final.and.(nhel(3).eq.0)))) then
                  call gaussint(pwave_func_real,minangle, pi, real_pwave_value,real_pwave_error)
                  if (abs(real_pwave_value) .gt. abs(max_pwave_value)) max_pwave_value = real_pwave_value
               endif
            enddo
          enddo
          if (ldebugout) write(*,"(A,F11.2,A,F9.4,A,F9.0,A)") " Lambda_FF: ", ffmassscale, "   ---   Re(pwave_0): ", 
     &      abs(max_pwave_value), "   ---   at ", energy, " GeV"
        enddo
        if (outputdetails.eq.2) then
              print*, ""
              print*, ""
              print*, ""
              print*, "***************************************"
              print*, "* Result for ", channel, " scattering: *"
              print*, "***************************************"
              print*, ""
              write(*,"(A,F9.0,A)") " Lambda_FF =", ffmassscale, " GeV"
              print*, ""
              write(*,"(A,A,A)") " prevents tree-level unitarity violation for on-shell ", channel," scattering"
              write(*,"(A,F9.0,A)") " up to sqrt(s) =", sqrts," GeV."
              print*, ""
        endif
        if (outputdetails.eq.1) write(*,"(A,F9.0,A)") " Result: Lambda_FF =", ffmassscale, " GeV"
        searchff_maxpwave_ind = ffmassscale
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
        searchff_maxpwave_ind = 0d0
      endif

      return
      end 


