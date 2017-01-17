
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
           ./amp_wpwp_scat.o\
           ./wpwpwpwp_anom.o \
           ./amp_wz_scat.o\
           ./wzwz_anom.o \
           ./amp_wzwa_scat.o\
           ./wzwa_anom.o \
           ./amp_wawz_scat.o\
           ./wawz_anom.o \
           ./amp_wa_scat.o\
           ./wawa_anom.o \
           ./amp_wwzz_scat.o\
           ./wwzz_anom.o \
           ./amp_zzww_scat.o\
           ./zzww_anom.o \
           ./amp_zz_scat.o\
           ./zzzz_anom.o \
           ./amp_wwaa_scat.o\
           ./wwaa_anom.o \
           ./amp_aaww_scat.o\
           ./aaww_anom.o \
           ./amp_wwza_scat.o\
           ./wwza_anom.o \
           ./amp_zaww_scat.o\
           ./zaww_anom.o \
           ./amp_zzaa_scat.o\
           ./zzaa_anom.o \
           ./amp_aa_scat.o\
           ./aaaa_anom.o \
           ./amp_zaaa_scat.o\
           ./zaaa_anom.o \
           ./amp_zazz_scat.o\
           ./zazz_anom.o \
           ./amp_aazz_scat.o\
           ./aazz_anom.o \
           ./amp_zzza_scat.o\
           ./zzza_anom.o \
           ./amp_zaza_scat.o\
           ./zaza_anom.o \
           ./amp_aaza_scat.o\
           ./aaza_anom.o \
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
amp_wpwp_scat.o: amp_wpwp_scat.F coupl.inc
wpwpwpwp_anom.o: wpwpwpwp_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_wz_scat.o: amp_wz_scat.F coupl.inc
wzwz_anom.o: wzwz_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_wzwa_scat.o: amp_wzwa_scat.F coupl.inc
wzwa_anom.o: wzwa_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_wawz_scat.o: amp_wawz_scat.F coupl.inc
wawz_anom.o: wawz_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_wa_scat.o: amp_wa_scat.F coupl.inc
wawa_anom.o: wawa_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_wwzz_scat.o: amp_wwzz_scat.F coupl.inc
wwzz_anom.o: wwzz_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_zzww_scat.o: amp_zzww_scat.F coupl.inc
zzww_anom.o: zzww_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_zz_scat.o: amp_zz_scat.F coupl.inc
zzzz_anom.o: zzzz_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_wwaa_scat.o: amp_wwaa_scat.F coupl.inc
wwaa_anom.o: wwaa_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_aaww_scat.o: amp_aaww_scat.F coupl.inc
aaww_anom.o: aaww_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_wwza_scat.o: amp_wwza_scat.F coupl.inc
wwza_anom.o: wwza_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_zaww_scat.o: amp_zaww_scat.F coupl.inc
zaww_anom.o: zaww_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_zzaa_scat.o: amp_zzaa_scat.F coupl.inc
zzaa_anom.o: zzaa_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_zaaa_scat.o: amp_zaaa_scat.F coupl.inc
zaaa_anom.o: zaaa_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_zzza_scat.o: amp_zzza_scat.F coupl.inc
zzza_anom.o: zzza_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_zazz_scat.o: amp_zazz_scat.F coupl.inc
zazz_anom.o: zazz_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_aazz_scat.o: amp_aazz_scat.F coupl.inc
aazz_anom.o: aazz_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_zaza_scat.o: amp_zaza_scat.F coupl.inc
zaza_anom.o: zaza_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_aaza_scat.o: amp_aaza_scat.F coupl.inc
aaza_anom.o: aaza_anom.F coupl.inc genps.inc an_couplings_4v.inc
amp_aa_scat.o: amp_aa_scat.F coupl.inc
aaaa_anom.o: aaaa_anom.F coupl.inc genps.inc an_couplings_4v.inc
searchff_diag.o: searchff_diag.F coupl.inc an_couplings_4v.inc global.inc
searchff_diag_ind.o: searchff_diag_ind.F coupl.inc an_couplings_4v.inc global.inc
searchff_maxpwave_ind.o: searchff_maxpwave_ind.F coupl.inc an_couplings_4v.inc global.inc


.F.o:
	$(FC) -c -o  $*.o  $(FFLAGS) $*.F

# ----------------------------------------------------------------------------

clean:
	rm -f *.o core *~ *.x run/get-formfactor


