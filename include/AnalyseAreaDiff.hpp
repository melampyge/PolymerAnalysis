
/* analysis on the averaged area difference
   per cell and per simulation */

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <map>
#include <set>
#include <stack>
#include <algorithm>
#include <tuple>

#include "data_structures.hpp"
#include "basic.hpp"
#include "read_write.hpp"

#define pi M_PI

/////////////////////////////////////////////////////////////////////////////////

class AnalyseAreaDiff {
  public:
    
    Simulation sim;
    Beads beads;
    double results;
    
    AnalyseAreaDiff(const char *datafilename, char *forc);
    ~AnalyseAreaDiff();
    void perform_analysis ();
    void write_analysis_results (const char *outfilepath);

    std::vector<double> calc_initial_areas (const std::vector<int> & nbpc, 
					    int ncells);
    double calc_avg_area_diff (const double * const *x, 
			       const double * const *y);
					    
};

/////////////////////////////////////////////////////////////////////////////////
