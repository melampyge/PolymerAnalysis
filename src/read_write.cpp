
/* i/o operations used in the analysis */

/////////////////////////////////////////////////////////////////////////////////

#include "read_write.hpp"

/////////////////////////////////////////////////////////////////////////////////

int read_single_int_data (hid_t file, const char *path_in_file) {
  /* wrapper to read singular integer data from hdf5 file 
  --note that buffer needs to be an array of size 1 for single entries-- */
   
  int buffer[1];
  buffer[0] = 0;
  hid_t dataset = H5Dopen(file, path_in_file, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL,
                          H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);

  return buffer[0];
}

/////////////////////////////////////////////////////////////////////////////////

double read_single_double_data (hid_t file, const char *path_in_file) {
  /* wrapper to read singular double data from hdf5 file 
  --note that buffer needs to be an array of size 1 for single entries-- */
   
  double buffer[1];
  buffer[0] = 0;
  hid_t dataset = H5Dopen(file, path_in_file, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                          H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);

  return buffer[0];
}

/////////////////////////////////////////////////////////////////////////////////

void read_int_array_data (const char *filename, const char *path_in_file, int *buffer) {
  /* wrapper to read integer array data from hdf5 file 
  --note that buffer needs to be the array size-- */

  // open the file pointer
  
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  // read the array
  
  hid_t dataset = H5Dopen(file, path_in_file, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL,
                          H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);
  H5Fclose(file);

  return;
}

/////////////////////////////////////////////////////////////////////////////////

void read_double_array_data (const char *filename, const char *path_in_file, double *buffer) {
  /* wrapper to read double array data from hdf5 file 
  --note that buffer needs to be the array size-- */

  // open the file pointer
  
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  // read the array
  
  hid_t dataset = H5Dopen(file, path_in_file, H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                          H5S_ALL, H5P_DEFAULT, buffer);

  H5Dclose(dataset);
  H5Fclose(file);

  return;
}

/////////////////////////////////////////////////////////////////////////////////

void read_single_pos_data (int step, hid_t dataset, hid_t dataspace,
                           double **x, double **y, int natoms) {
  /* read the position data at a single timestep */
  
  // read the data in the x direction
  
  // define the hyperslab in the dataset
  /* we are gonna reduce the data in the dataset 
   from 3 dimensions to two 2 dimensional arrays */
  
  hsize_t offset[3];
  offset[0] = step; offset[1] = 0; offset[2] = 0;
  
  hsize_t count[3];
  count[0] = 1; count[1] = 1; count[2] = natoms;
  
  // select a 2D hyperslab from the original 3D dataset
  
  herr_t status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                      offset, NULL, count, NULL);
  
  // define memory dataspace
  
  hsize_t dimsm[3];	// dimensions and sizes in each dimension
  dimsm[0] = 1; dimsm[1] = 1; dimsm[2] = natoms;
  hid_t memspace = H5Screate_simple(3, dimsm, NULL);
  
  // define memory hyperslab
  
  hsize_t offset_out[3];
  offset_out[0] = 0; offset_out[1] = 0; offset_out[2] = 0;
  
  hsize_t count_out[3];
  count_out[0] = 1; count_out[1] = 1; count_out[2] = natoms;
  
  status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                               offset_out, NULL, count_out, NULL);
  
  // read data from hyperslab in the file into the hyperslab in memory
  
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace,
                   dataspace, H5P_DEFAULT, x[step]);
  
  H5Sclose(memspace);
  
  // read the data in the y direction
  
  offset[0] = step; offset[1] = 1; offset[2] = 0;
  count[0] = 1; count[1] = 1; count[2] = natoms;
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                               offset, NULL, count, NULL);
  memspace = H5Screate_simple(3, dimsm, NULL);
  status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                               offset_out, NULL, count_out, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace,
                   dataspace, H5P_DEFAULT, y[step]);
  H5Sclose(memspace);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

void read_all_pos_data (const char *filename, double **x, double **y, int nsteps,
                        int natoms, const char *datapath) {
  /* read position data in hdf5 format all at once */

  // open the file pointer
  
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  // get the dataset in the file
  /* DATASET (H5D) is the raw data (either singular or in arrays) */
  
  hid_t dataset = H5Dopen(file, datapath, H5P_DEFAULT);
  
  // get dataspace of the selected dataset
  /* DATASPACE (H5S) describes the number of dimensions 
   and the size of the dataset in those dimension */
  
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

/////////////////////////////////////////////////////////////////////////////////

void write_single_analysis_data (double data, const char *outpath) {
  /* write analysis data with single parameter to the outfile */
  
  std::ofstream fl(outpath);
  fl << data << std::endl;
  fl.close();
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

void write_1d_analysis_data (const std::vector<double> &x,
                             const char *outpath) {
  /* write the 1d analysis data to the outfile */
  
  std::ofstream fl(outpath);
  for (int j = 0; j < x.size(); j++) {
    fl << j << "\t\t" << x[j] << std::endl;
  }
  fl.close();
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

void write_2d_analysis_data (const std::vector<double> &x,
			     const std::vector<double> &y, 
			     const char *outpath) {
  /* write the 2d analysis data to the outfile */
  
  std::ofstream fl(outpath);
  for (int j = 0; j < x.size(); j++) {
    fl << x[j] << "\t\t" << y[j] << std::endl;
  }
  fl.close();
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

void write_8d_analysis_data (const std::vector<double> &v1,
                             const std::vector<double> &v2,
                             const std::vector<double> &v3,
                             const std::vector<double> &v4,
                             const std::vector<double> &v5,
                             const std::vector<double> &v6,
                             const std::vector<double> &v7,
                             const std::vector<double> &v8,
                             const char *outpath) {
  /* write 8d analysis data to a single outfile */

  std::ofstream fl(outpath);
  for (unsigned j = 0; j < v1.size(); j++) {
    fl << v1[j] << "\t" << v2[j] << "\t" << v3[j] << "\t" << v4[j] <<
      "\t" << v5[j] << "\t" << v6[j] << "\t" << v7[j] << "\t" <<
      v8[j] << std::endl;
  }
  fl.close();

}

/////////////////////////////////////////////////////////////////////////////////

void write_multid_analysis_data (const std::vector<double> &w,
                                 const std::vector<double> &vx,
                                 const std::vector<double> &vy,
                                 const int ndata, const int nsteps,
                                 const char *outpath) {
  /* write multid analysis data to the outfile */
  
  std::ofstream fl(outpath);
  for (int step = 0; step < nsteps; step++) {
    for (int i = 0; i < ndata; i++) {
      for (int j = 0; j < ndata; j++) {
        int idx = step*ndata*ndata + ndata*i + j;
        fl << step << "\t" << i << "\t" << j << "\t" << w[idx] << "\t"
        << vx[idx] << "\t" << vy[idx] << "\t"
        << sqrt(vx[idx]*vx[idx]+vy[idx]*vy[idx]) << std::endl;
      }
    }
  }
  
  fl.close();
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

