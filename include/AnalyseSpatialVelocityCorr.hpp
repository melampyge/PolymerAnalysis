
/* analysis on spatial velocity correlation */

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
#include <tuple>

#include "omp.h"

#include "data_structures.hpp"
#include "basic.hpp"
#include "read_write.hpp"

#define pi M_PI
#define SUBTRACT_COM 		// whether to subtract center of mass from positions

/////////////////////////////////////////////////////////////////////////////////

class AnalyseSpatialVelocityCorr {
public:
  
  Simulation sim;
  Polymers polymers;
  std::vector<double> results;
  
  AnalyseSpatialVelocityCorr(const char *datafilename, char *forc);
  ~AnalyseSpatialVelocityCorr();
  void perform_analysis ();
  void write_analysis_results (const char *outfilepath);

  std::tuple<std::vector<double>, std::vector<double> > calc_velocity (const double * const *x,
								       const double * const *y,
								       int delta, int nvels);
  std::vector<double> calc_sp_vel_corr (const double * const *x, 
					const double * const *y,
					const std::vector<double> &vx,
					const std::vector<double> &vy,
					int delta, int ndata, int nvels);
							    
};

/////////////////////////////////////////////////////////////////////////////////
