
#pragma once

#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include "hdf5.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
T read_single_data (hid_t file, std::string path_in_file, T *buffer) {
  /* wrapper to read singular data from hdf5 file 
  --note that buffer needs to be an array of size 1 for single entries-- */
   
  const char *fl = path_in_file.c_str(); 
  hid_t dataset = H5Dopen(file, fl, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);

  return buffer[0];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
void read_array_data (std::string filename, std::string path_in_file, T *buffer) {
  /* wrapper to read array data from hdf5 file 
  --note that buffer needs to be the array size-- */

  const char *fl1 = filename.c_str();
  const char *fl2 = path_in_file.c_str();

  // open the file pointer
  
  hid_t file = H5Fopen(fl1, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  // read the array
  
  hid_t dataset = H5Dopen(file, fl2, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);
  H5Fclose(file);

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

/* read the position data at a single timestep */  
void read_single_pos_data(int step, hid_t dataset, hid_t dataspace, double **x, double **y, int natoms);
  
/* read position data in hdf5 format all at once */
void read_all_pos_data(std::string filename, double **x, double **y, int nsteps, int natoms, std::string datapath);

/* read the polarity data at a single timestep */
void read_single_pol_data(int step, hid_t dataset, hid_t dataspace, double **pol, int natoms);

/* read polarity data in hdf5 format all at once */
void read_all_pol_data(std::string filename, double **pol, int nsteps, int natoms, std::string datapath);
  
/* write analysis data with single parameter to the outfile */
void write_single_analysis_data(double data, std::string outpath);
  
/* write the 1d analysis data to the outfile */
void write_1d_analysis_data(double *x, int ndata, std::string outpath);
  
/* write the 1d vector analysis data to the outfile */
void write_1d_vec_analysis_data(const std::vector<double> &x, int ndata, std::string outpath);
  
/* write the 2d analysis data to the outfile */
void write_2d_analysis_data(double *x, double *y, double ndata, std::string outpath);
    
/* write multid analysis data to the outfile */
void write_multid_analysis_data(const std::vector<double> &w, const std::vector<double> &vx, const std::vector<double> &vy, const int ndata, const int nsteps, std::string outpath);
