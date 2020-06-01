#ifndef __WRITE_DATA__
#define __WRITE_DATA__

#include "lib/three_point_timeslice_fitting.h"
#include "edb_reader.h"


//************************************************************************
void print_output(string& path, string& name, prefactor& pf, fit_three_point_output& output_rl, fit_three_point_output& output_im,
                string& Zfile, map<string, vector< pair< pair<ENSEM::EnsemReal, ENSEM::EnsemReal> , ENSEM::EnsemReal >>>& fqsq,
                map<string, vector<ENSEM::EnsemReal> >& favg, map<string, vector<ENSEM::EnsemReal> >& fmean);


//************************************************************************
void print_extra_data(string& xmlini, string& Zfile, map<string, vector< pair< pair<ENSEM::EnsemReal, ENSEM::EnsemReal> , ENSEM::EnsemReal >>>& fqsq,
                map<string, vector<ENSEM::EnsemReal> >& favg);

#endif