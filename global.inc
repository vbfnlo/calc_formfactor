

c     Print additional output useful for debugging?
      logical ldebugout
      parameter (ldebugout=.false.)

c
c     Print short output or output with more details during calculation?
c
c     outputdetails = 0: (almost) no output during calculation
c     outputdetails = 1: short output during calculation
c     outputdetails = 2: very detailed output during calculation
c
      integer outputdetails
      common /output/outputdetails

      double precision pi
      parameter (pi = 3.141592653589793d0)

