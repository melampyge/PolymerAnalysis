
/* analysis on static structure factor */

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "omp.h"

#include "data_structures.hpp"
#include "basic.hpp"
#include "read_write.hpp"

#define pi M_PI

/////////////////////////////////////////////////////////////////////////////////

class AnalyseStaticStruct {
public:
  
  Simulation sim;
  Polymers polymers;
  std::tuple<std::vector<double>, std::vector<double> > results;
  
  AnalyseStaticStruct(const char *datafilename, char *forc);
  ~AnalyseStaticStruct();
  void perform_analysis ();
  void write_analysis_results (const char *outfilepath);

  std::tuple<std::vector<double>, std::vector<double> > calc_static_struct_lin_k (const double * const *x,
									    const double * const *y);
  std::tuple<std::vector<double>, std::vector<double> > calc_static_struct (const double * const *x,
									    const double * const *y);
};

/////////////////////////////////////////////////////////////////////////////////
