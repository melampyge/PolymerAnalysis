
/* data structures used in the analysis */

/////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <iostream>
#include <vector>

#include "read_write.hpp"
#include "basic.hpp"

/////////////////////////////////////////////////////////////////////////////////

class Simulation {
  /* abstract class for creating containers of general simulation data */
  public:
    
    Simulation() {};
    ~Simulation() {};
    int nsteps, nbeads, npols;
    double bond_length, dt, sigma, lx, ly;
    std::vector<int> nbpp;
    double kT = 1.;
    double gamma_0 = 1.;
    virtual void read_sim_info() const = 0;
  
};

/////////////////////////////////////////////////////////////////////////////////

class SimulationFilaments : public Simulation {
  /* container for general simulation data of filamentous polymers */
  public:
    
    SimulationFilaments(std::string filename);
    ~SimulationFilaments() {};
    double density, kappa, fp;
    void read_sim_data(std::string filename);
  
};

SimulationFilaments::SimulationFilaments(std::string filename) {
  
  // read in general simulation data

  nsteps = nbeads = npols = 0;
  bond_length = dt = sigma = lx = ly = 0.;
  density = kappa = fp = 0.;

  read_sim_data(filename);
  
  return;
}

void SimulationFilaments::read_sim_data(std::string filename) {
  /* read in the general simualtion data */
  
  const char *fl1 = filename.c_str(); 	// type conversion needed for hdf5
 
  // open the file pointer
  
  hid_t file = H5Fopen(fl1, H5F_ACC_RDONLY, H5P_DEFAULT);

  // create buffers to get single data --note that buffer needs to be an array of size 1 for single entries--
  
  int i_buffer[1];
  i_buffer[0] = 0;
  double d_buffer[1];
  d_buffer[0] = 0.;
  
  // read in the simulation parameters

  nsteps = read_integer_data(file, "/sim/nsteps", i_buffer);
  nbeads = read_integer_data(file, "/sim/nbeads", i_buffer);
  npols = read_integer_data(file, "/sim/npols", i_buffer);  
  int nbpp_buffer[npols];
  for (int i = 0; i < npols; i++) nbpp_buffer[i] = 0;
  read_integer_array(filename, "/sim/nbpp", nbpp_buffer);  
  for (int i = 0; i < npols; i++) nbpp.push_back(nbpp_buffer[i]);
  }
  // read in the model parameters
  
  bond_length = read_double_data(file, "/params/bl", d_buffer);
  sigma = read_double_data(file, "/params/sigma", d_buffer);
  dt = read_double_data(file, "/params/dt", d_buffer);
    
  // read in the box info

  lx = read_double_data(file, "/params/lx", d_buffer);
  ly = read_double_data(file, "/params/ly", d_buffer);  
  
  // read in model specific parameters
  
  density = read_double_data(file, "/params/density", d_buffer);
  kappa = read_double_data(file, "/params/kappa", d_buffer);
  fp = read_double_data(file, "/params/fp", d_buffer);
  
  H5Fclose(file);
  
  return;  
}

/////////////////////////////////////////////////////////////////////////////////

class SimulationCells : public Simulation {
  /* container for general simulation data for cells */
  public:
    
    SimulationCells(std::string filename);
    ~SimulationCells() {};
    double density, eps, fm, areak, kappa;
    void read_sim_data(std::string filename);
  
};

SimulationCells::SimulationCells(std::string filename) {
  
  // read in general simulation data

  nsteps = nbeads = npols = 0;
  bond_length = dt = sigma = lx = ly = 0.;
  density = eps = fm = areak = kappa = 0.;

  read_sim_data(filename);
  
  return;
}

void SimulationCells::read_sim_data(std::string filename) {
  /* read in the general simualtion data */
  
  const char *fl1 = filename.c_str(); 	// type conversion needed for hdf5
 
  // open the file pointer
  
  hid_t file = H5Fopen(fl1, H5F_ACC_RDONLY, H5P_DEFAULT);

  // create buffers to get single data --note that buffer needs to be an array of size 1 for single entries--
  
  int i_buffer[1];
  i_buffer[0] = 0;
  double d_buffer[1];
  d_buffer[0] = 0.;
  
  // read in the simulation parameters

  nsteps = read_integer_data(file, "/info/nsteps", i_buffer);
  nbeads = read_integer_data(file, "/info/nbeads", i_buffer);
  npols = read_integer_data(file, "/info/ncells", i_buffer);  
  int nsamp = read_integer_data(file, "/info/nsamp", i_buffer);
  int nbpp_buffer[npols];
  for (int i = 0; i < npols; i++) nbpp_buffer[i] = 0;
  read_integer_array(filename, "/cells/nbpc", nbpp_buffer);  
  for (int i = 0; i < npols; i++) nbpp.push_back(nbpp_buffer[i]);
  }
  
  // read in the model parameters
  
  bond_length = read_double_data(file, "/params/bl", d_buffer);
  sigma = read_double_data(file, "/params/sigma", d_buffer);  
  dt = read_double_data(file, "/info/dt", d_buffer);
  dt *= nsamp;

  // read in the box info

  lx = read_double_data(file, "/info/box/x", d_buffer);
  ly = read_double_data(file, "/info/box/y", d_buffer);

  // read in the model specific parameters
  
  eps = read_double_data(file, "/param/eps", d_buffer);
  density = read_double_data(file, "/param/rho", d_buffer);
  fm = read_double_data(file, "/param/fp", d_buffer);
  areak = read_double_data(file, "/param/areak", d_buffer);
  kappa = read_double_data(file, "/param/kappa", d_buffer);
  
  H5Fclose(file);
  
  return; 
}

/////////////////////////////////////////////////////////////////////////////////

template<class T>
class Beads {
  /* container for bead data */
  public:
    
    Beads(std::string filename, T sim);
    ~Beads();
    double **x;
    double **y;
  
    std::vector<int> get_pol_ids();
  
};

template<class T>
Beads::Beads(std::string filename, T sim) {
  
  allocate_2d_array(x, sim.nsteps, sim.nbeads);
  allocate_2d_array(y, sim.nsteps, sim.nbeads);
  read_all_pos_data(filename, x, y, sim.nsteps, sim.nbeads, "/beads/xu");
 
  return;
}

Beads::~Beads() {
  
  deallocate_2d_array(x, sim.nsteps);
  deallocate_2d_array(y, sim.nsteps);
  
  return;
}

std::vector<int> Beads::get_pol_ids () {
  
  std::vector<int> pid(sim.nbeads, 0);
  int k = 0;
  for (int n = 0; n < sim.npols; n++) {
    std::fill(pid.begin()+k, pid.begin()+k+nbpp[n], n);
    k += nbpp[n];
  }
  
  return pid;
}

/////////////////////////////////////////////////////////////////////////////////

template<class T>
class Polymers {
  /* container for polymer data */
  public:
    
    Polymers(std::string filename, T sim, std::string sim_type);
    ~Polymers();
    double **x;
    double **y;
    T sim;
  
};

template<class T>
Polymers::Polymers(std::string filename, T sim, std::string sim_type) {
  
  T sim = sim;
  allocate_2d_array(x, sim.nsteps, sim.nbeads);
  allocate_2d_array(y, sim.nsteps, sim.nbeads);
  std::string data_path = "/pols/comu";
  if (sim_type == "cells") 
    data_path = "/cells/comu";
  read_all_pos_data(filename, x, y, sim.nsteps, sim.npols, data_path);
 
}

Polymers::~Polymers() {
  
  deallocate_2d_array(x, sim.nsteps);
  deallocate_2d_array(y, sim.nsteps);
  
}

/////////////////////////////////////////////////////////////////////////////////
