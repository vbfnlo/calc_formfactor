ccc
c $Id: gauss.F,v 1.3 2004/04/16 13:14:21 mrauch Exp $
ccc
c Iterative Gaussian integration procedure
ccc
c Author: K.S. Koelbig, Cernlib program D103
c         adapted for HadCalc by Michael Rauch
ccc

        subroutine gaussint(integrand,lower,upper,integral,error)
        implicit none
        double precision integrand, lower, upper, integral
        external integrand

        integer nregions, neval, fail
        double precision error
 
        double precision clower, cupper, cdiff
        common /cgauss/ clower,cupper, cdiff

        clower = lower
        cupper = upper
        cdiff = cupper-clower

        call gaussad(1,integrand,1d-4,1d-10,2,
     &                     0,1000000,
     &                     nregions,neval,fail,integral,error)
        integral = integral * cdiff
        end

        subroutine gaussad(ndim,integrand,epsrel,epsabs,flags,
     &                     mineval,maxeval,
     &                     nregions,neval,fail,integral,error)
        implicit none

        integer ndim,flags,mineval,maxeval,nregions,neval,fail
        double precision epsrel,epsabs,integral,error
        double precision integrand
        external integrand

        integer i,j
        double precision f1,f2,s8,s16,diff,
     &                   lower,upper,c1,c2,u,
     &                   delta,zero,one,W,X
        logical accuracy

        double precision clower, cupper, cdiff
        common /cgauss/ clower,cupper, cdiff


        dimension W(12),X(12)
        data delta /1.d-12/
        data zero /0.d0/
        data one /1.d0/
        data W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
        data X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/

c initialize return values
        integral=0.d0
        error=0.d0
        nregions=0
        neval=0
        fail=0

!         if (mod(flags,4) .gt. 1) then
!           write(*,*) "Gaussad input parameters:"
!           write(*,*) "epsrel ", epsrel
!           write(*,*) "epsabs ", epsabs
!           write(*,*) "flags ", flags
!           write(*,*) "mineval ", mineval
!           write(*,*) "maxeval ", maxeval
!         endif

c This routine can only compute one-dimensional integrals
        if (ndim .ne. 1) then
          fail=-1
          return
        endif


        lower=zero

        do while( (one-lower) .gt. delta) 
          diff=one-lower
          accuracy=.FALSE.
          nregions=nregions+1

          do while(accuracy .eqv. .FALSE.)

            if ((dabs(diff) .le. delta) .or. (neval .gt. maxeval)) then
              fail=1
              return
            endif

            upper=lower+diff
            c1=0.5d0*(lower+upper)
            c2=c1-lower
            s8=0.d0
            s16=0.d0

            do j=1,4
              u=X(j)*c2
              f1 = integrand((c1+u)*cdiff+clower)
              f2 = integrand((c1-u)*cdiff+clower)
              s8=s8 + W(j)*( f1 + f2 )
            enddo

            do j=5,12
              u=X(j)*c2
              f1 = integrand((c1+u)*cdiff+clower)
              f2 = integrand((c1-u)*cdiff+clower)
              s16=s16 + W(j)*( f1 + f2 )
            enddo

            s8=s8*c2
            s16=s16*c2


c Check whether we have reached the desired accuracy
            accuracy=.TRUE.
            if ((epsabs .lt. dabs(s16-s8)) .and. 
     &          (epsrel*dabs(s16) .lt. dabs(s16-s8) ) ) then
              accuracy=.FALSE.
            endif

c some statistics
            neval=neval+24
            if (mod(flags,4) .ge. 3) then
              write(*,*) "Region: ", lower, " - ", upper
              write(*,101) "  ", s16, 
     &                     " +- ", dabs(s16-s8)
101              format( A,E16.8E3,A,E16.8E3)
            endif

            diff=0.5d0*diff

          enddo

          if (mod(flags,4) .ge. 1) then
!             write(*,102) "Integrating region ",lower," - ",upper," done"
102            format(A,E16.8E3,A,E16.8E3,A)
          endif

          integral=integral+s16
          error=error+dabs(s16-s8)
          lower=upper

        enddo

        if (neval .lt. mineval) then
          fail=1
        endif

        end

