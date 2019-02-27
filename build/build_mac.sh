#!/bin/bash

../configure --prefix="${HOME}/three_pt_analysis/three_pt_fit_install" \
             --with-adat="${HOME}/Scattering/adat_install" \
              --with-itpp="${HOME}/Scattering/itpp_install" \
              --with-fitting="${HOME}/3pt_analysis/fitting_lib_install" \
             CXXFLAGS="-O3 -mtune=native -fopenmp -fbounds-check" CFLAGS="-O3 -mtune=native -fopenmp -fbounds-check" \
             CXX=g++-8 CC=gcc-8
             #--with-gcc-8="/usr/local/Cellar/gcc/8.2.0/"
