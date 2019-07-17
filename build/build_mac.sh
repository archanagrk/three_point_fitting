#!/bin/bash

../configure --prefix="${HOME}/LQCDSoftware/three_pt_analysis/three_pt_fit_install" \
             --with-adat="${HOME}/LQCDSoftware/Kinematic_factor/adat_devel_install" \
              --with-itpp="${HOME}/LQCDSoftware/Scattering/itpp_install" \
              --with-fitting="${HOME}/LQCDSoftware/three_pt_analysis/fitting_lib_install" \
              --with-semble="${HOME}/LQCDSoftware/three_pt_analysis/semble_install" \
             #CXXFLAGS="-Wall -g -O3 -mtune=native -fopenmp -fbounds-check" CFLAGS="-Wall -g -O3 -mtune=native -fopenmp -fbounds-check" \
             CXXFLAGS="-Wall -O3 -mtune=native -fopenmp -fbounds-check" CFLAGS="-Wall -O3 -mtune=native -fopenmp -fbounds-check"\
             CXX=g++-9 CC=gcc-9
             #--with-gcc-8="/usr/local/Cellar/gcc/8.2.0/"
