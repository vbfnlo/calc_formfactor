c======================================================
c
c  Define common block containing KOPPLN inputs.
c
c  Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c  Dave Rainwater,    <rain@pheno.physics.wisc.edu>
c  Initial version:  1997 Mar 3
c
c======================================================
c
      double precision  alfas, xmt, alfa, xmz, xmw, sin2w, xmh, gfermi

      common /bkopin/   alfas, xmt, alfa, xmz, xmw, sin2w, xmh

      data xmh    /126.d0/                ! Higgs boson mass; do not set >= 2 * m_W, as widths are 
                                          ! set to zero! => on-shell VV -> H -> HH would cause
                                          ! problems!

      data xmz    /91.1876d0/             ! Z boson mass
      data xmw    /80.398d0/              ! W boson mass
      data gfermi /1.16639d-5/
      data alfa   /-1d0/                  ! choose alfa and sin2w negative for G_mu scheme
      data sin2w  /-1d0/                  ! => m_Z, m_W and G_F are used as EW input

      data alfas  /0.1184d0/              ! not important for tree-level VV scattering
      data xmt    /172.4d0/               ! not important for tree-level VV scattering


!     Issue a "make" after changing anything in here!

