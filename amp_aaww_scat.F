
      double complex function pwave_aaww_anom(theta)
      implicit none
      include "coupl.inc"
      
      double precision theta
      double precision s1
      complex*16 amp_aaww
      real*8 p(0:3,4)
      real*8 pi
      parameter(pi =3.141592653589793d0)  

      REAL*8 T
      INTEGER JC(4)

      ! get helicities and partial wave from main program
      INTEGER NHEL(4), IHEL
      integer pw
      common /partwave/ pw, nhel

      ! get center-of-mass energy squared from main program
      real*8 sc
      common /cms/ sc

      complex*16 aaww_anom
      external aaww_anom

      s1=sc

c   AA -> WW scattering

c defining all in- and outgoing momenta of the bosons      
c the scattering is assumed to be in the y-z-plane         
         

       p(0,1)=sqrt(s1)/2d0
       p(1,1)=0d0
       p(2,1)=0d0
       p(3,1)=Sqrt(s1/4d0)
	   
       p(0,2)=sqrt(s1)/2d0
       p(1,2)=0d0
       p(2,2)=0d0
       p(3,2)=-p(3,1)
       
       p(0,3)=sqrt(s1)/2d0
       p(1,3)=0d0
       p(2,3)=Sqrt(s1/4d0 - wmass**2)*sin(theta)
       p(3,3)=Sqrt(s1/4d0 - wmass**2)*cos(theta)
       
       p(0,4)=sqrt(s1)/2
       p(1,4)=0d0
       p(2,4)=-p(2,3)
       p(3,4)=-p(3,3)   

      DO IHEL=1,4
         JC(IHEL) = +1
      ENDDO

      amp_aaww = aaww_anom(P,NHEL,JC)

      if (pw.eq.2) then
            pwave_aaww_anom = 1.d0/(32.d0*pi) * amp_aaww * sin(theta) * 1./2.*(3.*cos(theta)**2-1)   ! pwave 2
      elseif (pw.eq.1) then
            pwave_aaww_anom = 1.d0/(32.d0*pi) * amp_aaww * sin(theta) * cos(theta)                   ! pwave 1
      else
            pwave_aaww_anom = 1.d0/(32.d0*pi) * amp_aaww * sin(theta)                                ! pwave 0
      endif

      ! symmetry factor for identical particles from wave function normalization:
      pwave_aaww_anom = pwave_aaww_anom / sqrt(2d0)

      end


      double precision function real_pwave_aaww(theta)
            double precision theta
            double complex amp, pwave_aaww_anom
            amp = pwave_aaww_anom(theta)
            real_pwave_aaww = DREAL(amp)
      end

      double precision function im_pwave_aaww(theta)
            double precision theta
            double complex amp, pwave_aaww_anom
            amp = pwave_aaww_anom(theta)
            im_pwave_aaww = DREAL(amp * DCMPLX(0.d0,-1.d0))
      end



