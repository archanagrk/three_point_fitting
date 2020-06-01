#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = -I@top_srcdir@/lib \
	 @CXXFLAGS@ @SEMBLE_CXXFLAGS@ @ADAT_CXXFLAGS@ @FITTING_LIB_CXXFLAGS@ @ITPP_CXXFLAGS@

AM_LDFLAGS  = -L@top_builddir@/lib \
  	 @LDFLAGS@ @ADAT_LDFLAGS@ @FITTING_LIB_LDFLAGS@ @SEMBLE_LDFLAGS@

LDADD = -lfitting @SEMBLE_LIBS@ @ADAT_LIBS@ @ITPP_LIBS@ @LIBS@ @FITTING_LIBS@

#
# The programs to build
# 
## NOTE: if you want the programs to be build by default when 'make' is run
##   add them to bin_PROGRAMS.  But, if you don't want them to be installed
##   in @prefix@/bin by 'make install' then add them to nodist_PROGRAMS.
##   Another option is to add it to check_PROGRAMS, which means it will only
##   be built by 'make check' and is not installed.


bin_PROGRAMS = fit_zt ensemble_fit_three_pt_fns make_elab fit_formfac_qsq

check_PROGRAMS = test_three_pt_fit

##############################################################

ensemble_fit_three_pt_fns_SOURCES   = exe/ensemble_fit_three_pt_fns.cc lib/three_point_timeslice_corr_functions.cc lib/three_point_timeslice_fitting.cc io/edb_reader.cc io/key_struct.cc\
													lib/three_point_timeslice_corr_functions_factory.cc lib/three_pt_fit_quality_factory.cc \
													lib/three_pt_fit_quality.cc lib/rotations.cc\
													io/xml_reader.cc io/write_data.cc

fit_zt_SOURCES		= test/fit_zt.cc test/z_timeslice_fitting.cc test/z_timeslice_functions.cc \
												test/z_timeslice_functions_factory.cc

fit_formfac_qsq_SOURCES    = exe/fit_formfac_qsq.cc lib/formfac_qsq_fitting.cc lib/formfac_qsq_functions.cc \
										 lib/formfac_qsq_functions_factory.cc io/read_formfac_qsq_data.cc

test_three_pt_fit_SOURCES	= test/test_three_pt_fit.cc lib/three_point_timeslice_corr_functions.cc lib/three_point_timeslice_fitting.cc \
                           lib/three_point_timeslice_corr_functions_factory.cc lib/three_pt_fit_quality_factory.cc \
                           lib/three_pt_fit_quality.cc

make_elab_SOURCES = utility/make_elab.cc 
##############################################################


#
# The following rule will always try to rebuild libadat.a when you 
# try to compile anything else
#
####DEPENDENCIES = build_lib ${lib_LIBRARIES}
DEPENDENCIES = ${lib_LIBRARIES}
${bin_PROGRAMS}: ${DEPENDENCIES}
${EXTRA_PROGRAMS}: ${DEPENDENCIES}

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_lib
build_lib:
	cd @top_builddir@/lib ; $(MAKE)
