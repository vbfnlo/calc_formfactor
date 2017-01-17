
# to change the compiler, comment out the two lines appropriate for your favourite compiler:

FC     = gfortran
FFLAGS = -g -fno-automatic -ffixed-line-length-none -O2
# gfortran debug flags:
# FFLAGS = -Wall -Wsurprising -Wno-unused -finit-integer=-1000000 -fbounds-check -Wtabs -g -fno-automatic -ffixed-line-length-none -ggdb3 -O0

# FC     = g77
# FFLAGS = -g -fno-automatic -ffixed-line-length-none -O2

# FC     = ifort
# FFLAGS = -save -extend_source



OBJECTSD = ./main_get_formfactors.o \
           ./amp_ww_scat.o\
           ./wwww_anom.o \
           ./amp_wz_scat.o\
           ./wzwz_anom.o \
           ./amp_wwzz_scat.o\
           ./wwzz_anom.o \
           ./amp_zzww_scat.o\
           ./zzww_anom.o \
           ./amp_zz_scat.o\
           ./zzzz_anom.o \
           ./anomal3.o \
           ./anomal4.o \
           ./dotproducts.o \
           ./gauss.o \
           ./vxxxxx.o \
	     ./readinput.o \
           ./readinput2.o \
           ./koppln-1905.o \
           ./eigenvalues.o \
           ./searchff_diag.o \
           ./searchff_diag_ind.o \
           ./searchff_maxpwave_ind.o

TARGETS = get-formfactor


get-formfactor:  $(OBJECTSD)
	$(FC) -o  run/get-formfactor $(FFLAGS) $(OBJECTSD)

#
# dependencies
#
main_get_formfactors.o: coupl.inc an_couplings_4v.inc main_get_formfactors.F
readinput.o: coupl.inc an_couplings_4v.inc readinput.F global.inc
readinput2.o: readinput.inc readinput2.F
koppln-1905.o: coupl.inc koppln.inc global.inc koppln-1905.F
anomal4.o: anomal4.F coupl.inc an_couplings_4v.inc
anomal3.o: anomal3.F coupl.inc an_couplings_4v.inc
amp_ww_scat.o: amp_ww_scat.F coupl.inc
wwww_anom.o: wwww_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_wz_scat.o: amp_wz_scat.F coupl.inc
wzwz_anom.o: wzwz_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_wwzz_scat.o: amp_wwzz_scat.F coupl.inc
wwzz_anom.o: wwzz_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_zzww_scat.o: amp_zzww_scat.F coupl.inc
zzww_anom.o: zzww_anom.F coupl.inc genps.inc an_couplings_4v.inc
searchff_diag.o: searchff_diag.F coupl.inc an_couplings_4v.inc global.inc
searchff_diag_ind.o: searchff_diag_ind.F coupl.inc an_couplings_4v.inc global.inc
searchff_maxpwave_ind.o: searchff_maxpwave_ind.F coupl.inc an_couplings_4v.inc global.inc


.F.o:
	$(FC) -c -o  $*.o  $(FFLAGS) $*.F

# ----------------------------------------------------------------------------

clean:
	rm -f *.o core *~ *.x run/get-formfactor


