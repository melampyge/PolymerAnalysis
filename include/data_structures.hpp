
/* data structures used in the analysis */

/////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <stdexcept>

#include "read_write.hpp"
#include "basic.hpp"

/////////////////////////////////////////////////////////////////////////////////

class Simulation {
  /* base class for creating containers of general simulation data */
public:
  
  Simulation() {};
  Simulation(std::string filename, std::string filsorcells);
  ~Simulation() {};
  
  
  std::string simtype = "";
  
  // general simulation parameters
  
  int nsteps, nbeads, npols;
  double bond_length, dt, sigma, lx, ly;
  std::vector<int> nbpp;
  double kT, gamma_0;
  
  // read write functions
  
  void read_filaments_data(std::string filename);
  void read_cells_data(std::string filename);
  
  // model specific parameters
  
  double density, kappa, fp, eps, fm, areak;
};

/////////////////////////////////////////////////////////////////////////////////

class Beads {
  /* container for bead data */
  public:
  
    Beads() {};
    Beads(std::string filename, Simulation sim);
    ~Beads();
    double **x;
    double **y;
    int nsteps;
    int nbeads;
  
    std::vector<int> get_pol_ids(int npols, const std::vector<int> &nbpp);
  
};

/////////////////////////////////////////////////////////////////////////////////

class Polymers {
  /* container for polymer data */
  public:
  
    Polymers() {};
    Polymers(std::string filename, Simulation sim, std::string forc);
    ~Polymers();
    double **x;
    double **y;
    int nsteps;
    int npols;
  
};

/////////////////////////////////////////////////////////////////////////////////
