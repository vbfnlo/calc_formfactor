C----------------------  subroutine readinput ----------------------------
C
      subroutine readinput
      implicit none
c
c this subroutine reads input parameters from the input data file 'input.dat'
c and stores it in common blocks; additioal input parameters can be changed
c in the file koppln.inc
      include "global.inc"
      include "coupl.inc"
      include "an_couplings_4v.inc"

      print*, ""
      print*, "Reading in anomalous couplings parameter:"
      print*, ""

      call loadfile("input.dat",.true.)
      call read_real("SQRT_S",sqrts,50000d0)
      call read_real("FFEXP",ffexponent,2d0)
      call read_real("FS0",fs0,0.0d0)
      call read_real("FS1",fs1,0.0d0)
      call read_real("FM0",fm0,0.0d0)
      call read_real("FM1",fm1,0.0d0)
      call read_real("FM2",fm2,0.0d0)
      call read_real("FM3",fm3,0.0d0)
      call read_real("FM4",fm4,0.0d0)
      call read_real("FM5",fm5,0.0d0)
      call read_real("FM6",fm6,0.0d0)
      call read_real("FM7",fm7,0.0d0)
      call read_real("FT0",ft0,0.0d0)
      call read_real("FT1",ft1,0.0d0)
      call read_real("FT2",ft2,0.0d0)
      call read_real("FT5",ft5,0.0d0)
      call read_real("FT6",ft6,0.0d0)
      call read_real("FT7",ft7,0.0d0)
      call read_real("FT8",ft8,0.0d0)
      call read_real("FT9",ft9,0.0d0)
      call read_real("FB", fb, 0.0d0)
      call read_real("FW", fw, 0.0d0)
      call read_real("FWWW", fwww, 0.0d0)
      call read_real("FWW", fWw, 0.0d0)
      call read_real("FBT", fbtilde, 0.0d0)
      call read_real("FWWWT", fwwwtilde, 0.0d0)
      call read_real("FWWT", fWwtilde, 0.0d0)
      call read_real("FDWT", fdwtilde, 0.0d0)
      call read_real("FPHID", fphid, 0.0d0)

      call read_int("OUTPUT_DETAILS", outputdetails, 0)
      call closefile

      fs0_0=fs0
      fs1_0=fs1
      fm0_0=fm0
      fm1_0=fm1
      fm2_0=fm2
      fm3_0=fm3
      fm4_0=fm4
      fm5_0=fm5
      fm6_0=fm6
      fm7_0=fm7
      ft0_0=ft0
      ft1_0=ft1
      ft2_0=ft2
      ft5_0=ft5
      ft6_0=ft6
      ft7_0=ft7
      ft8_0=ft8
      ft9_0=ft9
      fw_0=fw
      fb_0=fb
      fwww_0=fwww
      fww_0=fww
      fbtilde_0=fbtilde
      fwwtilde_0=fwwtilde
      fwwwtilde_0=fwwwtilde
      fdwtilde_0=fdwtilde
      fphid_0=fphid

c
c set Standard Model couplings. SM input is in koppln.inc
c

      CALL KOPPLN(2)
      call ctrans

      print*, ""
      print*, "Setting widths of H, Z, W to zero to prevent violation of"
      print*, "gauge invariance at large energies due to MadGraph-like"
      print*, "complex-mass-scheme."
      print*, ""
      HWIDTH = 0d0
      ZWIDTH = 0d0
      WWIDTH = 0d0



      end
