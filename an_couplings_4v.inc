      real*8 fs0, fs1, ft5, ft0, ft1, ft6, fm0, fm2, fm6
      real*8 fm3, fm4, fm1, fm5, fm7, ft2, ft7, ft8, ft9
      real*8 fwww, fw, fb, fww, fbb, fbtilde, fwwtilde, fbbtilde, fwwwtilde, fdwtilde
      real*8 fphid
      real*8 fphi1, fbw, fdw, fwtilde, fbwtilde
      real*8 fs0_0, fs1_0, ft5_0, ft0_0, ft1_0, ft6_0, fm0_0, fm2_0, fm6_0
      real*8 fm3_0, fm4_0, fm1_0, fm5_0, fm7_0, ft2_0, ft7_0, ft8_0, ft9_0
      real*8 fwww_0, fw_0, fb_0, fww_0, fbb_0, fbtilde_0, fwwtilde_0, fbbtilde_0, fwwwtilde_0, fdwtilde_0
      real*8 fphid_0
      real*8 fphi1_0, fbw_0, fdw_0, fwtilde_0, fbwtilde_0
      real*8 ffvalue, ffexponent, ffmassscale
      logical formfact
  
c   All anomalous coupling constants are meant to be already divided by lambda^4, e.g.: fw:= fw/lambda^4
      
      common /anomal_v4/ fs0, fs1, ft0, ft1, ft5, ft6, fm0, fm2, fm6
      common /anomal_v4/ fm3, fm4, fm1, fm5, fm7, ft2, ft7, ft8, ft9
      common /anomal_v4/ fwww, fw, fb, fww, fbb, fbtilde, fwwtilde, fbbtilde, fwwwtilde, fdwtilde
      common /anomal_v4/ fphid
      common /anomal_v4/ fs0_0, fs1_0, ft0_0, ft1_0, ft5_0, ft6_0, fm0_0, fm2_0, fm6_0
      common /anomal_v4/ fm3_0, fm4_0, fm1_0, fm5_0, fm7_0, ft2_0, ft7_0, ft8_0, ft9_0
      common /anomal_v4/ fwww_0, fw_0, fb_0, fww_0, fbb_0, fbtilde_0, fwwtilde_0, fbbtilde_0, fwwwtilde_0, fdwtilde_0
      common /anomal_v4/ fphid_0
      common /anomal_v4/ ffvalue, ffexponent, ffmassscale
      common /anomal_v4/ formfact
      

      real*8 sqrts
      common /maxenergy/ sqrts
