
/* analysis on pair correlation function */

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

#include "data_structures.hpp"
#include "basic.hpp"
#include "read_write.hpp"

#define pi M_PI

/////////////////////////////////////////////////////////////////////////////////

class AnalysePairCorr {
public:
  
  Simulation sim;
  Polymers polymers;
  std::vector<double> results;
  
  AnalysePairCorr(const char *datafilename, char *forc);
  ~AnalysePairCorr();
  void perform_analysis ();
  void write_analysis_results (const char *outfilepath);

  std::vector<double> calc_pair_corr (const double * const *x,
				      const double * const *y);
};

/////////////////////////////////////////////////////////////////////////////////
