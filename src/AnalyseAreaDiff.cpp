
/* analysis on the averaged area difference
   per cell and per simulation */

/////////////////////////////////////////////////////////////////////////////////

#include "AnalyseAreaDiff.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

AnalyseAreaDiff::AnalyseAreaDiff (const char *datafilename, char *forc) : 
    sim(datafilename, forc), beads(datafilename, sim) {
 
  // print information
  
  cout << "\nData is loaded successfully for the following file: \n" <<
    datafilename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "npols = " << sim.npols << endl;
  cout << "nbeads = " << sim.nbeads << endl;
  
  results = 0.;
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

AnalyseAreaDiff::~AnalyseAreaDiff () { }

/////////////////////////////////////////////////////////////////////////////////

void AnalyseAreaDiff::perform_analysis () {
  
  results = calc_avg_area_diff(beads.x, beads.y);
  
  return;
}


/////////////////////////////////////////////////////////////////////////////////

void AnalyseAreaDiff::write_analysis_results (const char *outfilepath) {

  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  write_single_analysis_data(results, outfilepath);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

vector<double> AnalyseAreaDiff::calc_initial_areas (const vector<int> & nbpc, 
						    int ncells) {
  /* calculate the initial target areas of cells */
  
  vector<double> A0(ncells, 0.);
  for (int j = 0; j < ncells; j++) {
    double rad = nbpc[j]*0.5/(2.*pi);
    A0[j] = 0.9*pi*rad*rad;
  }
  
  return A0;
}

/////////////////////////////////////////////////////////////////////////////////

double AnalyseAreaDiff::calc_avg_area_diff (const double * const *x, 
					    const double * const *y) {
  /* calculate the area covered by cells per frame 
   and return the averaged area over all the frames */

  vector<double> init_areas = calc_initial_areas(sim.nbpp, sim.npols);
  
  double adiff = 0.;
  
  for (int step = 0; step < sim.nsteps; step++) {
    
    cout << "steps / nsteps : " << step << " / " << sim.nsteps << endl;
    
    double adiff_per_cell = 0.;
    int k = 0;
    
    for (int n = 0; n < sim.npols; n++) {
      double area = 0.;
      
      for (int j = 0; j < sim.nbpp[n]-1; j++) {
        area += x[step][k]*y[step][k+1] - x[step][k+1]*y[step][k];
        k++;
      }		// beads of the cell
      
      area += x[step][k]*y[step][k-sim.nbpp[n]+1] - x[step][k-sim.nbpp[n]+1]*y[step][k];
      k++;
      area /= 2.;
      adiff_per_cell += area - init_areas[n];
      
    }	 	// cells
    
    adiff += adiff_per_cell/sim.npols;
    
  }     	// timesteps
  
  adiff /= sim.nsteps;
  
  return adiff;
}

/////////////////////////////////////////////////////////////////////////////////

