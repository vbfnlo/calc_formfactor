!* call d-functions from amplitude routine:

!      double precision dfunction
!      external dfunction

!      pwave_wwzz_anom = 1.d0/(32.d0*pi) * amp_wwzz * sin(theta) * DFUNCTION(pw, nhel(1)-nhel(2), nhel(3)-nhel(4), theta)

!* helicity selection in main routine has to be adapted / extended


      REAL*8 FUNCTION DFUNCTION(dfktj,dfktlambda,dfktmu,theta) 

c   calculates d functions which are needed to calculate the coefficients
c   of the partial wave analysis

      implicit none
       
      double precision theta, dfktlambda, dfktmu, dfktj

      if (dfktj.eq.0) then
        DFUNCTION = 1.0d0

      elseif (dfktj.eq.1) then
        if ((dfktlambda.eq.0).and.(dfktmu.eq.0)) then
          DFUNCTION = cos(theta)
        elseif ((dfktlambda.eq.1).and.(dfktmu.eq.0)) then
          DFUNCTION = -sin(theta)/(sqrt(2.0d0))
        elseif ((dfktlambda.eq.0).and.(dfktmu.eq.-1)) then
          DFUNCTION = -sin(theta)/(sqrt(2.0d0))
        elseif ((dfktlambda.eq.0).and.(dfktmu.eq.1)) then
          DFUNCTION = sin(theta)/(sqrt(2.0d0))
        elseif ((dfktlambda.eq.-1).and.(dfktmu.eq.0)) then
          DFUNCTION = sin(theta)/(sqrt(2.0d0))
        elseif ((dfktlambda.eq.1).and.(dfktmu.eq.1)) then
          DFUNCTION = (1.0d0+cos(theta))/(2.0d0)
        elseif ((dfktlambda.eq.1).and.(dfktmu.eq.-1)) then
          DFUNCTION = (1.0d0-cos(theta))/(2.0d0)
        elseif ((dfktlambda.eq.-1).and.(dfktmu.eq.1)) then
          DFUNCTION = (1.0d0-cos(theta))/(2.0d0)
        elseif ((dfktlambda.eq.-1).and.(dfktmu.eq.-1)) then
          DFUNCTION = (1.0d0+cos(theta))/(2.0d0)

       else
       write(*,*) "d function not available for this helicity combination"
        endif

      elseif (dfktj.eq.2) then

        if ((dfktlambda.eq.2).and.(dfktmu.eq.2)) then
          DFUNCTION = ((1.0d0+cos(theta))/(2.0d0))**2
        elseif ((dfktlambda.eq.2).and.(dfktmu.eq.1)) then
          DFUNCTION = -((1.0d0+cos(theta))/(2.0d0))*sin(theta)
        elseif ((dfktlambda.eq.2).and.(dfktmu.eq.0)) then
          DFUNCTION = sqrt(6.0d0)/4.0d0*(sin(theta))**2
        elseif ((dfktlambda.eq.2).and.(dfktmu.eq.-1)) then
          DFUNCTION = -((1.0d0-cos(theta))/(2.0d0))*sin(theta)
        elseif ((dfktlambda.eq.2).and.(dfktmu.eq.-2)) then
          DFUNCTION = ((1.0d0-cos(theta))/(2.0d0))**2

        elseif ((dfktlambda.eq.1).and.(dfktmu.eq.2)) then
          DFUNCTION = ((1.0d0+cos(theta))/(2.0d0))*sin(theta)
        elseif ((dfktlambda.eq.1).and.(dfktmu.eq.1)) then
          DFUNCTION = ((1.0d0+cos(theta))/(2.0d0))*
     &                (2.0d0*cos(theta)-1.0d0)
        elseif ((dfktlambda.eq.1).and.(dfktmu.eq.0)) then
          DFUNCTION = -sqrt(3.0d0/2.0d0)*sin(theta)*cos(theta)
        elseif ((dfktlambda.eq.1).and.(dfktmu.eq.-1)) then
          DFUNCTION = ((1.0d0-cos(theta))/(2.0d0))*
     &                (2.0d0*cos(theta)+1.0d0)
        elseif ((dfktlambda.eq.1).and.(dfktmu.eq.-2)) then
          DFUNCTION = -((1.0d0-cos(theta))/(2.0d0))*sin(theta)

        elseif ((dfktlambda.eq.0).and.(dfktmu.eq.2)) then
          DFUNCTION = sqrt(6.0d0)/4.0d0*(sin(theta))**2
        elseif ((dfktlambda.eq.0).and.(dfktmu.eq.1)) then
          DFUNCTION = sqrt(3.0d0/2.0d0)*sin(theta)*cos(theta)
        elseif ((dfktlambda.eq.0).and.(dfktmu.eq.0)) then
          DFUNCTION = (3.0d0/2.0d0*((cos(theta))**2)-1.0d0/2.0d0)
        elseif ((dfktlambda.eq.0).and.(dfktmu.eq.-1)) then
          DFUNCTION = -sqrt(3.0d0/2.0d0)*sin(theta)*cos(theta)
        elseif ((dfktlambda.eq.0).and.(dfktmu.eq.-2)) then
          DFUNCTION = sqrt(6.0d0)/4.0d0*(sin(theta))**2

        elseif ((dfktlambda.eq.-1).and.(dfktmu.eq.2)) then
          DFUNCTION = ((1.0d0-cos(theta))/(2.0d0))*sin(theta)
        elseif ((dfktlambda.eq.-1).and.(dfktmu.eq.1)) then
          DFUNCTION = ((1.0d0-cos(theta))/(2.0d0))*
     &                (2.0d0*cos(theta)+1.0d0)
        elseif ((dfktlambda.eq.-1).and.(dfktmu.eq.0)) then
          DFUNCTION = sqrt(3.0d0/2.0d0)*sin(theta)*cos(theta)
        elseif ((dfktlambda.eq.-1).and.(dfktmu.eq.-1)) then
          DFUNCTION = ((1.0d0+cos(theta))/(2.0d0))*
     &                (2.0d0*cos(theta)-1.0d0)
        elseif ((dfktlambda.eq.-1).and.(dfktmu.eq.-2)) then
          DFUNCTION = -((1.0d0+cos(theta))/(2.0d0))*sin(theta)

        elseif ((dfktlambda.eq.-2).and.(dfktmu.eq.2)) then
          DFUNCTION = ((1.0d0-cos(theta))/(2.0d0))**2
        elseif ((dfktlambda.eq.-2).and.(dfktmu.eq.1)) then
          DFUNCTION = ((1.0d0-cos(theta))/(2.0d0))*sin(theta)
        elseif ((dfktlambda.eq.-2).and.(dfktmu.eq.0)) then
          DFUNCTION = sqrt(6.0d0)/4.0d0*(sin(theta))**2
        elseif ((dfktlambda.eq.-2).and.(dfktmu.eq.-1)) then
          DFUNCTION = ((1.0d0+cos(theta))/(2.0d0))*sin(theta)
        elseif ((dfktlambda.eq.-2).and.(dfktmu.eq.-2)) then
          DFUNCTION = ((1.0d0+cos(theta))/(2.0d0))**2

       else
       write(*,*) "d function not available for this helicity combination"

       endif

       else
       write(*,*) "d function not available for j > 2"

      endif


      end
