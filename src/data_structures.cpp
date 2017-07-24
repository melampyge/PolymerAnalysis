
/* data structures used in the analysis */

/////////////////////////////////////////////////////////////////////////////////

#include "data_structures.hpp"

/////////////////////////////////////////////////////////////////////////////////

Simulation::Simulation(std::string filename, std::string filsorcells) {
  
  nsteps = nbeads = npols = 0;
  bond_length = dt = sigma = lx = ly = 0.;
  kT = gamma_0 = 1.;
  density = kappa = fp = eps = fm = areak = 0.;
  
  if (filsorcells == "filaments") {
    read_filaments_data(filename);
    simtype = filsorcells;
  }
  else if (filsorcells == "cells") {
    read_cells_data(filename);
    simtype = filsorcells;
  }
  else {
    throw std::runtime_error("filsorcells should be filaments or cells!");
  }

  return;
}

void Simulation::read_filaments_data(std::string filename) {
  /* read in the general filaments simualtion data */
  
  const char *fl1 = filename.c_str(); 	// type conversion needed for hdf5
  
  // open the file pointer
  
  hid_t file = H5Fopen(fl1, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  // create buffers to get single data
  // --note that buffer needs to be an array of size 1 for single entries--
  
  int i_buffer[1];
  i_buffer[0] = 0;
  double d_buffer[1];
  d_buffer[0] = 0.;
  
  // read in the simulation parameters
  
  nsteps = read_single_data(file, "/sim/nsteps", i_buffer);
  nbeads = read_single_data(file, "/sim/nbeads", i_buffer);
  npols = read_single_data(file, "/sim/npols", i_buffer);
  int nbpp_buffer[npols];
  for (int i = 0; i < npols; i++) nbpp_buffer[i] = 0;
  read_array_data(filename, "/sim/nbpp", nbpp_buffer);
  for (int i = 0; i < npols; i++) nbpp.push_back(nbpp_buffer[i]);
  
  // read in the model parameters
  
  bond_length = read_single_data(file, "/params/bl", d_buffer);
  sigma = read_single_data(file, "/params/sigma", d_buffer);
  dt = read_single_data(file, "/params/dt", d_buffer);
  
  // read in the box info
  
  lx = read_single_data(file, "/params/lx", d_buffer);
  ly = read_single_data(file, "/params/ly", d_buffer);
  
  // read in model specific parameters
  
  density = read_single_data(file, "/params/density", d_buffer);
  kappa = read_single_data(file, "/params/kappa", d_buffer);
  fp = read_single_data(file, "/params/fp", d_buffer);
  
  H5Fclose(file);
  
  return;
}

void Simulation::read_cells_data(std::string filename) {
  /* read in the general cells simualtion data */
  
  const char *fl1 = filename.c_str(); 	// type conversion needed for hdf5
  
  // open the file pointer
  
  hid_t file = H5Fopen(fl1, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  // create buffers to get single data --note that buffer needs to be an array of size 1 for single entries--
  
  int i_buffer[1];
  i_buffer[0] = 0;
  double d_buffer[1];
  d_buffer[0] = 0.;
  
  // read in the simulation parameters
  
  nsteps = read_single_data(file, "/info/nsteps", i_buffer);
  nbeads = read_single_data(file, "/info/nbeads", i_buffer);
  npols = read_single_data(file, "/info/ncells", i_buffer);
  int nsamp = read_single_data(file, "/info/nsamp", i_buffer);
  int nbpp_buffer[npols];
  for (int i = 0; i < npols; i++) nbpp_buffer[i] = 0;
  read_array_data(filename, "/cells/nbpc", nbpp_buffer);
  for (int i = 0; i < npols; i++) nbpp.push_back(nbpp_buffer[i]);
  
  // read in the model parameters
  
  bond_length = read_single_data(file, "/params/bl", d_buffer);
  sigma = read_single_data(file, "/params/sigma", d_buffer);
  dt = read_single_data(file, "/info/dt", d_buffer);
  dt *= nsamp;
  
  // read in the box info
  
  lx = read_single_data(file, "/info/box/x", d_buffer);
  ly = read_single_data(file, "/info/box/y", d_buffer);
  
  // read in the model specific parameters
  
  eps = read_single_data(file, "/param/eps", d_buffer);
  density = read_single_data(file, "/param/rho", d_buffer);
  fm = read_single_data(file, "/param/fp", d_buffer);
  areak = read_single_data(file, "/param/areak", d_buffer);
  kappa = read_single_data(file, "/param/kappa", d_buffer);
  
  H5Fclose(file);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

Beads::Beads(std::string filename, Simulation sim) {
  
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

Polymers::Polymers(std::string filename, Simulation sim, std::string forc) {
  
  std::string data_path = "/cells/comu";
  if (forc == "filaments")
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
