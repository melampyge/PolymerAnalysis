
/* data structures used in the analysis */

/////////////////////////////////////////////////////////////////////////////////

#include "data_structures.hpp"

/////////////////////////////////////////////////////////////////////////////////

Simulation::Simulation(const char *filename, char *filsorcells) {
  
  nsteps = nbeads = npols = 0;
  bond_length = dt = sigma = lx = ly = 0.;
  kT = gamma_0 = 1.;
  density = kappa = fp = eps = fm = areak = 0.;
    
  if (strcmp(filsorcells, "filaments") == 0) {
    read_filaments_data(filename);
    simtype = filsorcells;
  }
  else if (strcmp(filsorcells, "cells") == 0) {
    read_cells_data(filename);
    simtype = filsorcells;
  }

  return;
}

void Simulation::read_filaments_data(const char *filename) {
  /* read in the general filaments simualtion data */
    
  // open the file pointer
  
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  // read in the simulation parameters
  
  nsteps = read_single_int_data(file, "/sim/nsteps");
  nbeads = read_single_int_data(file, "/sim/nbeads");
  npols = read_single_int_data(file, "/sim/npols");
  int nbpp_buffer[npols];
  for (int i = 0; i < npols; i++) nbpp_buffer[i] = 0;
  read_int_array_data(filename, "/sim/nbpp", nbpp_buffer);
  for (int i = 0; i < npols; i++) nbpp.push_back(nbpp_buffer[i]);
  
  // read in the model parameters
  
  bond_length = read_single_double_data(file, "/params/bl");
  sigma = read_single_double_data(file, "/params/sigma");
  dt = read_single_double_data(file, "/params/dt");
  
  // read in the box info
  
  lx = read_single_double_data(file, "/params/lx");
  ly = read_single_double_data(file, "/params/ly");
  
  // read in model specific parameters
  
  density = read_single_double_data(file, "/params/density");
  kappa = read_single_double_data(file, "/params/kappa");
  fp = read_single_double_data(file, "/params/fp");
  
  H5Fclose(file);
  
  return;
}

void Simulation::read_cells_data(const char *filename) {
  /* read in the general cells simualtion data */
    
  // open the file pointer
  
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  // read in the simulation parameters
  
  nsteps = read_single_int_data(file, "/info/nsteps");
  nbeads = read_single_int_data(file, "/info/nbeads");
  npols = read_single_int_data(file, "/info/ncells");
  int nsamp = read_single_int_data(file, "/info/nsamp");
  
  int nbpp_buffer[npols];
  for (int i = 0; i < npols; i++) nbpp_buffer[i] = 0;
  read_int_array_data(filename, "/cells/nbpc", nbpp_buffer);
  for (int i = 0; i < npols; i++) nbpp.push_back(nbpp_buffer[i]);
  
  // read in the model parameters
  
  bond_length = read_single_double_data(file, "/param/bl");
  sigma = read_single_double_data(file, "/param/sigma");
  dt = read_single_double_data(file, "/info/dt");
  dt *= nsamp;
  
  // read in the box info
  
  lx = read_single_double_data(file, "/info/box/x");
  ly = read_single_double_data(file, "/info/box/y");
  
  // read in the model specific parameters
  
  eps = read_single_double_data(file, "/param/eps");
  density = read_single_double_data(file, "/param/rho");
  fm = read_single_double_data(file, "/param/fp");
  areak = read_single_double_data(file, "/param/areak");  
  kappa = read_single_double_data(file, "/param/kappa");  
  
  H5Fclose(file);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

Beads::Beads(const char *filename, Simulation sim) {
  
  allocate_2d_array(x, sim.nsteps, sim.nbeads);
  allocate_2d_array(y, sim.nsteps, sim.nbeads);
  read_all_pos_data(filename, x, y, sim.nsteps, sim.nbeads, "/beads/xu");
  nsteps = sim.nsteps;
  nbeads = sim.nbeads;
  
  return;
}

Beads::~Beads() {
  
  deallocate_2d_array(x, nsteps);
  deallocate_2d_array(y, nsteps);
  
  return;
}

std::vector<int> Beads::get_pol_ids (int npols, const std::vector<int> &nbpp) {
  
  std::vector<int> pid(nbeads, 0);
  int k = 0;
  for (int n = 0; n < npols; n++) {
    std::fill(pid.begin()+k, pid.begin()+k+nbpp[n], n);
    k += nbpp[n];
  }
  
  return pid;
}

/////////////////////////////////////////////////////////////////////////////////

Polymers::Polymers(const char *filename, Simulation sim, const char *forc) {
  
  char *data_path = "/cells/comu";
  if (strcmp(forc, "filaments") == 0)
    data_path = "/pols/comu";
  allocate_2d_array(x, sim.nsteps, sim.npols);
  allocate_2d_array(y, sim.nsteps, sim.npols);
  read_all_pos_data(filename, x, y, sim.nsteps, sim.npols, data_path);
  nsteps = sim.nsteps;
  npols = sim.npols;
  
  return;
}

Polymers::~Polymers() {
  
  deallocate_2d_array(x, nsteps);
  deallocate_2d_array(y, nsteps);
  
}

/////////////////////////////////////////////////////////////////////////////////
