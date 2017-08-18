
/* analysis on vorticity field */

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
#include <tuple>

#include "data_structures.hpp"
#include "basic.hpp"
#include "read_write.hpp"

#define pi M_PI

/////////////////////////////////////////////////////////////////////////////////

class AnalyseVorticity {
public:
  
  Simulation sim;
  Polymers polymers;
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, int, int, double> results;
  
  AnalyseVorticity(const char *datafilename, char *forc);
  ~AnalyseVorticity();
  void perform_analysis ();
  void write_analysis_results (const char *outfilepath, const char *outfilepath_2);

  std::tuple<std::vector<double>, std::vector<double> > calc_velocity (const double * const *x,
								       const double * const *y,
								       int delta, int nvels);
							  
  void add_velocity_to_bins (int ix_bin, int iy_bin, int step,
			    std::vector<double> &vx_bin, std::vector<double> &vy_bin,
			    std::vector<int> &ncells_per_bin_per_step,
			    double vxc, double vyc, 
			    int nbins, int ncells);
			    
  std::tuple<std::vector<double>, std::vector<double>, std::vector<std::vector<int> > > calc_velocity_per_bin (
			      const std::vector<double> &vx, 
			      const std::vector<double> &vy, 
			      const double * const *x,
			      const double * const *y,
			      double wbin, int nbins, double woverlap,
			      double lx, double ly, 
			      int ncells, int nsteps);
			      
  std::tuple<std::vector<double>, double> calc_vorticity (const std::vector<std::vector<int> > &ncells_per_bin,
							  const std::vector<double> &vx_bin, 
							  const std::vector<double> &vy_bin,
							  double wbin, int nbins, 
							  double lx, double ly, 
							  int ncells, int nsteps);
};

/////////////////////////////////////////////////////////////////////////////////

