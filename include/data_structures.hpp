
/* data structures used in the analysis */

/////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <vector>
#include <cstring>

#include "read_write.hpp"
#include "basic.hpp"

/////////////////////////////////////////////////////////////////////////////////

class Simulation {
  /* base class for creating containers of general simulation data */
  public:
    
    Simulation() {};
    Simulation(const char *filename, char *filsorcells);
    ~Simulation() {};
    
    // type of the simulation --filaments or cells--
    
    char *simtype = "";
    
    // general simulation parameters
    
    int nsteps, nbeads, npols;
    double bond_length, dt, sigma, lx, ly;
    std::vector<int> nbpp;
    double kT, gamma_0;
    
    // read write functions
    
    void read_filaments_data(const char *filename);
    void read_cells_data(const char *filename);
    
    // model specific parameters
    
    double density, kappa, fp, eps, fm, areak;
};

/////////////////////////////////////////////////////////////////////////////////

class Beads {
  /* container for bead data */
  public:
  
    Beads() {};
    Beads(const char *filename, Simulation sim);
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
    Polymers(const char *filename, Simulation sim, const char *forc);
    ~Polymers();
    double **x;
    double **y;
    int nsteps;
    int npols;
  
};

/////////////////////////////////////////////////////////////////////////////////
