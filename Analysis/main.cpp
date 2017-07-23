
/* main module to start a particular analysis */

/////////////////////////////////////////////////////////////////////////////////

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
#include <tuple>
#include <stack>
#include <algorithm>

#include "read_write.hpp"
#include "data_structures.hpp"
#include "basic.hpp"
#include "analyser.hpp"
#include "AnalyseNeighbours.hpp"

#define pi M_PI

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {
  /* argv[1] -> datafilename         (path of the data file)
     argv[2] -> beadsorpols          (make a choice between reading beads or polymers)
     argv[3] -> filsorcells          (make a choice of sim type between filaments or cells)
     argv[4...] -> outfilename(s)... (path of the analysis results file)
   
   */
  
  // get command line arguments
  
  string datafilename = argv[1];
  string beadsorpols = argv[2];
  string filsorcells = argv[3];
  string outfilename = argv[4];
  string outfilename_2 = argv[5];
  
  // load the data based on the chosen options
  
  if (filsorcells == "cells") {
    SimulationCells sim(datafilename);
  }
  else if (filsorcells == "filaments") {
    SimulationFilaments sim(datafilename);
  }
  else {
    throw runtime_error("filsorcells should be filaments or cells!");
  }
  
  if (beadsorpols == "beads") {
    Beads data(datafilename, sim);
  }
  else if (beadsorpols == "polymers") {
    Polymers data(datafilename, sim);
  }
  else {
    throw runtime_error("beadsorpols should be beads or polymers!");
  }
  
  // print information
  
  cout << "\nData is loaded successfully for the following file: \n" <<
    datafilename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "npols = " << sim.npols << endl;
  cout << "nbeads = " << sim.nbeads << endl;
  
  // perform the analysis
  
  AnalyseNeighbours analysis(sim, data&);
  analysis.perform_analysis();
  
  // write the computed results
  
  analysis.write_analysis_results(outfilename, outfilename_2);
  
  return 0;
}  

/////////////////////////////////////////////////////////////////////////////////

