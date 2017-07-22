
/* analysis on neighbours
 such as number of neighbours */

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <map>
#include <set>
#include <stack>
#include <algorithm>

#include "data_structures.hpp"
#include "basic.hpp"
#include "read_write.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
class AnalyseNeighbours : public Analyser {
public:
  
  T sim;
  Beads *beads;
  
  AnalyseNeighbours(T sim_, Beads *beads_);
  ~AnalyseNeighbours();
  double perform_analysis();
  void write_analysis_results(double data, string outfilepath);

  void AnalyseNeighbours::build_linked_cell_list(const double * const *x,
                                                 const double * const *y,
                                                 std::vector<int> & heads,
                                                 std::vector<int> & llist,
                                                 double wbin, int nboxes, int nsize,
                                                 int step, int natoms, double l);
 
  double AnalyseNeighbours::calc_num_neighbours (const double * const *x,
                                                 const double * const *y);
  
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////
