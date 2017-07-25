
/* analysis on the density histogram */

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <stack>
#include <algorithm>

#include "data_structures.hpp"
#include "basic.hpp"
#include "read_write.hpp"

#define pi M_PI

/////////////////////////////////////////////////////////////////////////////////

class AnalyseDensityHistogram {
public:
  
  Simulation sim;
  Beads beads;
  std::vector<double> results;
  
  AnalyseDensityHistogram(const char *datafilename, char *forc);
  ~AnalyseDensityHistogram();
  void perform_analysis ();
  void write_analysis_results (const char *outfilepath);

  std::vector<double> calc_densities (const double * const *x, 
				 const double * const *y);
};

/////////////////////////////////////////////////////////////////////////////////
