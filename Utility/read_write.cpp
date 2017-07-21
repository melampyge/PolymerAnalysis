
/* helper function in data reading and writing */

// COMPILATION AND RUN COMMANDS:
// g++ -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/read_write.cpp -lhdf5 -o read_write
// -

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "read_write.hpp"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_single_pos_data (int step, hid_t dataset, hid_t dataspace, double **x, double **y, int natoms) {
  /* read the position data at a single timestep */
 
  // read the data in the x direction
  
  // define the hyperslab in the dataset
  /* we are gonna reduce the data in the dataset from 3 dimensions to two 2 dimensional arrays */
    
  hsize_t offset[3];
  offset[0] = step; offset[1] = 0; offset[2] = 0;
  
  hsize_t count[3];
  count[0] = 1; count[1] = 1; count[2] = natoms;
  
  // select a 2D hyperslab from the original 3D dataset
  
  herr_t status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

  // define memory dataspace
  
  hsize_t dimsm[3];	// dimensions and sizes in each dimension
  dimsm[0] = 1; dimsm[1] = 1; dimsm[2] = natoms;
  hid_t memspace = H5Screate_simple(3, dimsm, NULL);
  
  // define memory hyperslab
  
  hsize_t offset_out[3];
  offset_out[0] = 0; offset_out[1] = 0; offset_out[2] = 0;
  
  hsize_t count_out[3];
  count_out[0] = 1; count_out[1] = 1; count_out[2] = natoms;
  
  status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
  
  // read data from hyperslab in the file into the hyperslab in memory 
  
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, x[step]); 
      
  H5Sclose(memspace);
  
  // read the data in the y direction
  
  offset[0] = step; offset[1] = 1; offset[2] = 0;
  count[0] = 1; count[1] = 1; count[2] = natoms;
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);   
  memspace = H5Screate_simple(3, dimsm, NULL);
  status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, y[step]);   
  H5Sclose(memspace);
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_all_pos_data (string filename, double **x, double **y, int nsteps, int natoms, string datapath) {
  /* read position data in hdf5 format all at once */
 
  const char *fl1 = filename.c_str();
  const char *fl2 = datapath.c_str();

  // open the file pointer
  
  hid_t file = H5Fopen(fl1, H5F_ACC_RDONLY, H5P_DEFAULT);

  // get the dataset in the file
  /* DATASET (H5D) is the raw data (either singular or in arrays) */
  
  hid_t dataset = H5Dopen(file, fl2, H5P_DEFAULT);
  
  // get dataspace of the selected dataset
  /* DATASPACE (H5S) describes the number of dimensions and the size of the dataset in those dimension */
  
  hid_t dataspace = H5Dget_space(dataset);
 
  /* READ THE DATA
  load the data into the memory step by step 
  register the data at each step to the array
  */
  
  for (int step = 0; step < nsteps; step++) {
    read_single_pos_data(step, dataset, dataspace, x, y, natoms);
  }

  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Fclose(file);
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_single_analysis_data (double data, string outpath) {
  /* write analysis data with single parameter to the outfile */
 
  const char *outfl = outpath.c_str();
  ofstream fl(outfl);
  fl << data << endl;
  fl.close();
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_1d_analysis_data (double *x, int ndata, string outpath) {
  /* write the 1d analysis data to the outfile */
  
  const char *outfl = outpath.c_str();
  ofstream fl(outfl);
  for (int j = 0; j < ndata; j++) {
    fl << j << "\t\t" << x[j] << endl;
  }
  fl.close();
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_1d_vec_analysis_data (const vector<double> &x, int ndata, string outpath) {
  /* write the 1d vector analysis data to the outfile */
  
  const char *outfl = outpath.c_str();
  ofstream fl(outfl);
  for (int j = 0; j < ndata; j++) {
    fl << j << "\t\t" << x[j] << endl;
  }
  fl.close();
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_2d_analysis_data (double *x, double *y, double ndata, string outpath) {
  /* write the 2d analysis data to the outfile */
  
  const char *outfl = outpath.c_str();
  ofstream fl(outfl);
  for (int j = 0; j < ndata; j++) {
    fl << x[j] << "\t\t" << y[j] << endl;
  }
  fl.close();
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_multid_analysis_data (const vector<double> &w, 
			    const vector<double> &vx, const vector<double> &vy, 
			    const int ndata, const int nsteps, string outpath) {
  /* write multid analysis data to the outfile */

  const char *outfl = outpath.c_str();
  ofstream fl(outfl);
  for (int step = 0; step < nsteps; step++) {
    for (int i = 0; i < ndata; i++) {
      for (int j = 0; j < ndata; j++) {
	int idx = step*ndata*ndata + ndata*i + j;
	fl << step << "\t" << i << "\t" << j << "\t" << w[idx] << "\t" << vx[idx] << "\t" << vy[idx] << "\t" << sqrt(vx[idx]*vx[idx]+vy[idx]*vy[idx]) << endl;
      }
    }
  }	
    
  fl.close();
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
