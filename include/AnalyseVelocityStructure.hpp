
/* analysis on velocity structure factors */

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

class AnalyseVelocityStructure {
public:
  
  Simulation sim;
  Polymers polymers;
  std::vector<double> results;
  
  AnalyseVelocityStructure(const char *datafilename, char *forc);
  ~AnalyseVelocityStructure();
  void perform_analysis ();
  void write_analysis_results (const char *outfilepath);
							    
};

/////////////////////////////////////////////////////////////////////////////////
