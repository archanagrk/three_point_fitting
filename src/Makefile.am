#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = -I@top_srcdir@/src \
	 @CXXFLAGS@ @ADAT_CXXFLAGS@ @FITTING_LIB_CXXFLAGS@ @ITPP_CXXFLAGS@

AM_LDFLAGS  = -L@top_builddir@/src \
  	 @LDFLAGS@ @ADAT_LDFLAGS@ @FITTING_LIB_LDFLAGS@

LDADD = -lfitting @ADAT_LIBS@ @ITPP_LIBS@ @LIBS@ @FITTING_LIBS@

#
# The programs to build
# 
## NOTE: if you want the programs to be build by default when 'make' is run
##   add them to bin_PROGRAMS.  But, if you don't want them to be installed
##   in @prefix@/bin by 'make install' then add them to nodist_PROGRAMS.
##   Another option is to add it to check_PROGRAMS, which means it will only
##   be built by 'make check' and is not installed.


bin_PROGRAMS = fit_zt fit_const


##############################################################

fit_const_SOURCES   = fit_const.cc three_point_timeslice_corr_functions.cc three_point_timeslice_fitting.cc \
													three_point_timeslice_corr_functions_factory.cc three_pt_fit_quality_factory.cc \
													three_pt_fit_quality.cc

fit_zt_SOURCES		= fit_zt.cc z_timeslice_fitting.cc z_timeslice_functions.cc \
												z_timeslice_functions_factory.cc

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
