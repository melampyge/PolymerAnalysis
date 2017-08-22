
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
#include "gen_velocity_grid.hpp"

#define pi M_PI
typedef  std::tuple<std::vector<std::vector<std::vector<double> > >, 
    std::vector<std::vector<std::vector<double> > >, 
    std::vector<std::vector<std::vector<double> > >,
    int, int, double, double,
    std::vector<double>, std::vector<double> > VorticityContainer;
 
/////////////////////////////////////////////////////////////////////////////////

class AnalyseVorticity {
public:
  
  Simulation sim;
  Polymers polymers;
  VorticityContainer results;
  
  AnalyseVorticity(const char *datafilename, char *forc);
  ~AnalyseVorticity();
  void perform_analysis ();
  void write_analysis_results (const char *outfilepath, const char *outfilepath_2,
      const char *outfilepath_3, const char *outfilepath_4);

  std::tuple<std::vector<std::vector<std::vector<double> > >, double, double,
    std::vector<double>, std::vector<double> > calc_vorticity (
		const std::vector<std::vector<std::vector<double> > > &vx_bin, 
		const std::vector<std::vector<std::vector<double> > > &vy_bin,
		double wbin, int nbins, double lx, double ly, int nsteps); 
};

/////////////////////////////////////////////////////////////////////////////////

