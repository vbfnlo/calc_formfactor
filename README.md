# VBFNLO utility to calculate form factors

* Version: 1.4.0
* Release date: 18 January 2017
* License: GPL v2 (see file COPYING)

---

* Contributions from:
  C. Englert, B. Feigl, J. Frank, M. Rauch, O. Schlimpert, F. Schissler
* https://www.itp.kit.edu/vbfnlo
* vbfnlo@itp.kit.edu

--- 

This program belongs to the program package VBFNLO and can calculate
input parameters needed for anomalous gauge boson coupling studies
with VBFNLO. 
As especially the pure operators for anomalous quartic gauge boson 
couplings lead to a violation of tree-level unitarity within the energy 
range of the LHC, special care has to be taken to avoid this unphysical
behaviour. Within VBFNLO we have opted for the use of a dipole form 
factor and this tool can calculate the maximal form factor scale 
Lambda_FF which is allowed for a given input of coupling parameters, 
assuming the form factor shape 

            FF = 1 / ( 1 + s / Lambda_FF^2 )^FFexp .

The form factor is determined by calculating on-shell VV scattering and 
computing the zeroth partial wave of the amplitude. As unitarity criterion 
the absolute value of the real part of the zeroth partial wave has to be 
below 0.5 [Barger, Cheung, Han, Phillips; Phys.Rev. D42 (1990) 3052].
Each channel in VV -> VV scattering (with V = W / Z / gamma) is checked 
individually. Additionally, channels with the same electrical charge Q of 
the VV system are combined [Gounaris, Layssac, Renard; 
Phys.Lett. B332 (1994) 146 (their definition of the partial wave expansion 
differs from ours by a factor of 2)].

For more details please refer to
* diploma thesis of B. Feigl (in German only)
* diploma thesis of O. Schlimpert (in German only)
at https://www.itp.kit.edu/publications/diploma .

The full list of operators with anomalous gauge boson couplings which is 
implemented in VBFNLO and probed here can be found in the manual of VBFNLO 
version 2.7 (and later) at 
https://www.itp.kit.edu/vbfnlo .

Further reading on tree-level unitarity with anomalous couplings:
* Barger, Cheung, Han, Phillips; Phys.Rev. D42 (1990) 3052
* Baur, Zeppenfeld; Phys.Lett. B201 (1988) 383
* Baur, Zeppenfeld; Nucl.Phys. B308 (1988) 127
* Gounaris, Layssac, Renard; Phys.Lett. B332 (1994) 146; arXiv:hep-ph/9311370
* Gounaris, Layssac, Paschalis, Renard; Z.Phys. C66 (1995) 619; arXiv:hep-ph/9409260



## Usage

* run `make` in the source directory
* the binary and input file are located in subfolder "run"
* modify `run/input.dat`
* run `./get-formfactor` from the folder "run"

The following compilers have been tested and compile this tool properly:
gfortran (default), g77, ifort.
If you want to switch the compiler, uncomment the appropriate lines in the 
makefile.

By running `./batchrun` all operators can be checked sequentially with a 
given coupling strength. Issuing the command without arguments displays 
a short manual.



## Input

In the file `run/input.dat` the following parameters can be set:
* Maximum center-of-mass energy that should be considered
* Desired exponent of the form factor.
  You need at least
  - FFexp = 1 for the dimension-6 operators and
  - FFexp = 2 for the dimension-8 operators.
* Coupling constants of the anomalous couplings



## Output

The program calculates the energy at which tree-level unitarity would
be violated without a form factor and the form factor scale that
ensures tree-level unitarity up to the given energy.
The form factor scales are also written out to the file `run/ff_output.dat`.

The form factor scales are determined in three ways:

* Considering only the helicity combination with the largest contribution
  to the zeroth partial wave, separately for each VV -> VV channel, with
  V including W-bosons, Z-bosons and photons.
  This is the simplest approach. Only taking the largest helicity combination
  may be slightly too optimistic. This option is useful for checking the relative
  importance of different channels.

* Considering all helicity combinations with contributions to the zeroth
  partial wave by diagonalizing the T-matrix, separately for each VV -> VV channel,
  where V includes W-bosons, Z-bosons and photons.
  This option gives form factor scales that are more restrictive than option 1.
  The values are valid for processes which involve only one channel.

* Considering all helicity combinations with contributions to the zeroth
  partial wave by diagonalizing the T-matrix for all VV -> VV channels 
  with V being a W-boson, Z-boson or photon.
  Using this approach gives the most stringent values for form factor and 
  unitarity violation scale. They also hold for processes with more than one
  quartic vertex contributing.


## Changelog

Version 1.4.0:
* Add dimension 8 operator O_S,2

Version 1.3.0:
* Bugfix in printout of coupling constants: negative values are now printed
  as well.

Version 1.3 beta:
* Include operator O_phiD .

Version 1.2.1:
* Bugfix for form factor output option 1:
  When using FFexp larger than minimal value, the form factor scale for some
  operators was too high. Output of option 2 and 3 was not affected.
* README: Added some guidelines how to use the results from the different
  calculation options.

Version 1.2.0 beta:
* Include all photon channels for all three calculation options.
* Include W+W+ -> W+W+ channel.

Version 1.1.0 beta:
* Some channels with external photons have been included for the first
  form factor calculation option.

Version 1.0.0:
* First public release.

Version 1.00 beta:
* Cleaned up output to screen.
* Now ifort and g77 work as well.
* Beta to first public release.

Version 0.30:
* Scan energy range instead of only checking at maximal energy.
* Set maximal energy to 14 TeV.
* Stop program when resulting form factor scale would be very small
  (FF exponent too small).
* Small shell script for batch runs added.

Version 0.20:
* First version to include different channels simultaneously.
* Set boson widths to zero to avoid gauge invariance violation at high energies.

Version 0.12:
* Prepared code for diagonalizing T-matrix
* Split main routine
* Include WW->ZZ channel

Version 0.11:
* New integration routine
* Lower integration bound for WW 0.1°

Version 0.10:
* Initial version in svn

