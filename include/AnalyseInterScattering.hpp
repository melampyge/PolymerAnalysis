
/* analysis on intermediate scattering function */

/////////////////////////////////////////////////////////////////////////////////

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
#define SUBTRACT_COM 		// whether to subtract center of mass from positions

/////////////////////////////////////////////////////////////////////////////////

class AnalyseInterScattering {
public:
  
  Simulation sim;
  Polymers polymers;
  std::tuple<std::vector<double>, std::vector<double> > results;
  
  AnalyseInterScattering(const char *datafilename, char *forc);
  ~AnalyseInterScattering();
  void perform_analysis ();
  void write_analysis_results (const char *outfilepath);

  std::tuple<std::vector<double>, std::vector<double> > calc_inter_scattering (const double * const *x,
									       const double * const *y);
					  
};

/////////////////////////////////////////////////////////////////////////////////
