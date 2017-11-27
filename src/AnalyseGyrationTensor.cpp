
/* analysis on gyration tensor */

/* 
 *
 * The idea is to reduce the per bead data 
 * into per polymer data. 
 * Analysis can be conducted later in Python.
 *
 */
 
/////////////////////////////////////////////////////////////////////////////////

#include "AnalyseGyrationTensor.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

AnalyseGyrationTensor::AnalyseGyrationTensor (const char *datafilename, char *forc) : 
    sim(datafilename, forc), beads(datafilename, sim) {
 
  // print information
  
  cout << "\nData is loaded successfully for the following file: \n" <<
    datafilename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "npols = " << sim.npols << endl;
  cout << "nbeads = " << sim.nbeads << endl;
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

AnalyseGyrationTensor::~AnalyseGyrationTensor () { }

/////////////////////////////////////////////////////////////////////////////////

void AnalyseGyrationTensor::perform_analysis () {
 
  results = calc_gyration_tensor(beads.x, beads.y);
  
  return;
}


/////////////////////////////////////////////////////////////////////////////////

void AnalyseGyrationTensor::write_analysis_results (const char *outfilepath) {

  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  vector<vector<double> > xcm, ycm, r11, r12, r22;
  vector<double> time;
  tie(time, xcm, ycm, r11, r12, r22) = results;
  write_gyration_tensor_data(time, xcm, ycm, 
      r11, r12, r22, outfilepath);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

gyrationTensorType AnalyseGyrationTensor::calc_gyration_tensor (const double * const *x,
						const double * const *y) {
  /* calculate the gyration tensor */
 
  // calculate the center of mass of pols

  vector<vector<double> > xcm(sim.nsteps, vector<double>(sim.npols, 0.));
  vector<vector<double> > ycm(sim.nsteps, vector<double>(sim.npols, 0.));
  vector<double> time(sim.nsteps, 0.);

  cout << "Calculating center of mass" << endl;
  for (int step = 0; step < sim.nsteps; step++) {
    int bead_id = 0;
    time[step] = sim.dt*step;

    for (int n = 0; n < sim.npols; n++) {
      for (int j = 0; j < sim.nbpp[n]; j++) {
        
        xcm[step][n] += x[step][bead_id];
        ycm[step][n] += y[step][bead_id];
        bead_id++;
      }   // beads per polymer

      xcm[step][n] /= sim.nbpp[n];
      ycm[step][n] /= sim.nbpp[n];
    }     // polymers 
  }       // time steps
  cout << "Finished center of mass calculation" << endl;

  // calculate the gyration tensor

  vector<vector<double> > r11(sim.nsteps, vector<double>(sim.npols, 0.));
  vector<vector<double> > r12(sim.nsteps, vector<double>(sim.npols, 0.));
  vector<vector<double> > r22(sim.nsteps, vector<double>(sim.npols, 0.));
  
  cout << "Calculating the gyration tensor" << endl;
  for (int step = 0; step < sim.nsteps; step++) {
    int bead_id = 0;

    for (int n = 0; n < sim.npols; n++) {
      for (int j = 0; j < sim.nbpp[n]; j++) {
        
        double dx = x[step][bead_id] - xcm[step][n];
        double dy = y[step][bead_id] - ycm[step][n];
        r11[step][n] += dx*dx;
        r12[step][n] += dx*dy;
        r22[step][n] += dy*dy;
        bead_id++;
      }    // beads per polymer

      r11[step][n] /= sim.nbpp[n];
      r12[step][n] /= sim.nbpp[n];
      r22[step][n] /= sim.nbpp[n];
    }      // polymers
  }        // time steps
  cout << "Finished the gyration tensor calculation" << endl;

  return make_tuple(time, xcm, ycm, r11, r12, r22);
}

/////////////////////////////////////////////////////////////////////////////////

