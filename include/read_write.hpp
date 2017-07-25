
/* i/o operations used in the analysis */

/////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <fstream>
#include <vector>
#include <cmath>
#include "hdf5.h"

/////////////////////////////////////////////////////////////////////////////////

/* wrapper to read singular integer data from hdf5 file 
--note that buffer needs to be an array of size 1 for single entries-- */
int read_single_int_data (hid_t file, const char *path_in_file);
  
/////////////////////////////////////////////////////////////////////////////////

/* wrapper to read singular double data from hdf5 file 
--note that buffer needs to be an array of size 1 for single entries-- */
double read_single_double_data (hid_t file, const char *path_in_file);

/////////////////////////////////////////////////////////////////////////////////

/* wrapper to read integer array data from hdf5 file 
--note that buffer needs to be the array size-- */
void read_int_array_data (const char *filename, const char *path_in_file, int *buffer);

/////////////////////////////////////////////////////////////////////////////////

/* wrapper to read double array data from hdf5 file 
--note that buffer needs to be the array size-- */
void read_double_array_data (const char *filename, const char *path_in_file, double *buffer);

/////////////////////////////////////////////////////////////////////////////////

/* read the position data at a single timestep */
void read_single_pos_data (int step, hid_t dataset, hid_t dataspace,
                           double **x, double **y, int natoms);

/////////////////////////////////////////////////////////////////////////////////

/* read position data in hdf5 format all at once */
void read_all_pos_data (const char *filename, double **x, double **y, int nsteps,
                        int natoms, const char *datapath);

/////////////////////////////////////////////////////////////////////////////////

/* write analysis data with single parameter to the outfile */
void write_single_analysis_data (double data, const char *outpath);

/////////////////////////////////////////////////////////////////////////////////

/* write the 1d analysis data to the outfile */
void write_1d_analysis_data (const std::vector<double> &x,
                             const char *outpath);

/////////////////////////////////////////////////////////////////////////////////

/* write the 2d analysis data to the outfile */
void write_2d_analysis_data (const std::vector<double> &x,
			     const std::vector<double> &y,
                             const char *outpath);

/////////////////////////////////////////////////////////////////////////////////

/* write multid analysis data to the outfile */
void write_multid_analysis_data (const std::vector<double> &w,
                                 const std::vector<double> &vx,
                                 const std::vector<double> &vy,
                                 const int ndata, const int nsteps,
                                 const char *outpath);

/////////////////////////////////////////////////////////////////////////////////
