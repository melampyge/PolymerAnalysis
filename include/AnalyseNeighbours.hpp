
/* analysis on neighbours
 such as number of neighbours */

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

/////////////////////////////////////////////////////////////////////////////////

class AnalyseNeighbours {
public:
  
  Simulation sim;
  Beads beads;
  std::tuple<double, double> results;
  
  AnalyseNeighbours(const char *datafilename, char *forc);
  ~AnalyseNeighbours();
  std::tuple<double, double> perform_analysis ();
  void write_analysis_results (const char *outfilepath,
                               const char *outfilepath_2);
  void build_linked_cell_list(const double * const *x,
                             const double * const *y,
                             std::vector<int> & heads,
                             std::vector<int> & llist,
                             const double wbin, const int nboxes,
                             const int nsize, const int step,
                             const int natoms, const double l);
  std::set<int> get_neighs_of_a_cell (const int start, const int end,
                                const double * const *x,
                                const double * const *y,
                                const int step, const double rcut,
                                const int nboxes,
                                const double lx, const double ly,
                                const std::vector<int> & cid,
                                const std::vector<int> & heads,
                                const std::vector<int> & llist);
  std::map<int, std::set<int> > build_neigh_list (const double * const *x,
                                       const double * const *y,
                                       const std::vector<int> & cid,
                                       const std::vector<int> & nbpp,
                                       const std::vector<int> & heads,
                                       const std::vector<int> & llist,
                                       const int npols,
                                       const double rcut,
                                       const int nsteps,
                                       const int step,
                                       const int nboxes,
                                       const double lx,
                                       const double ly);
  std::tuple<double, double> calc_num_neighbours (const double * const *x,
                                             const double * const *y);
};

/////////////////////////////////////////////////////////////////////////////////
