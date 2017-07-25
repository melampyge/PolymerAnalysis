
/* analysis on lattice order */

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
#include "cartesian.hpp"
#include "kdtree.hpp"

/////////////////////////////////////////////////////////////////////////////////

class AnalyseLatticeOrder {
  public:
    
    Simulation sim;
    Polymers polymers;
    double results;
    
    AnalyseLatticeOrder(const char *datafilename, char *forc);
    ~AnalyseLatticeOrder();
    void perform_analysis ();
    void write_analysis_results (const char *outfilepath); 

    void populate_extended_pos_array (std::vector<Point<2> > &points, int npoints, int ncells, 
				      int step, const double * const *x, 
				      const double * const *y, double lx, double ly);
    double calc_glob_order (const double * const *x, 
			    const double * const *y, 
			    int nn);    
};

/////////////////////////////////////////////////////////////////////////////////
