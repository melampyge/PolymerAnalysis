
/* analysis on pair correlation function */

/////////////////////////////////////////////////////////////////////////////////

#include "AnalysePairCorr.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

AnalysePairCorr::AnalysePairCorr (const char *datafilename, char *forc) : 
    sim(datafilename, forc), polymers(datafilename, sim, forc) {
 
  // print information
  
  cout << "\nData is loaded successfully for the following file: \n" <<
    datafilename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "npols = " << sim.npols << endl;
  cout << "nbeads = " << sim.nbeads << endl;
  
  fill(results.begin(), results.end(), 0.);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

AnalysePairCorr::~AnalysePairCorr () { }

/////////////////////////////////////////////////////////////////////////////////

void AnalysePairCorr::perform_analysis () {
  
  results = calc_pair_corr(polymers.x, polymers.y);
  
  return;
}


/////////////////////////////////////////////////////////////////////////////////

void AnalysePairCorr::write_analysis_results (const char *outfilepath) {

  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  write_1d_analysis_data(results, outfilepath);
  
  return;
  
}

/////////////////////////////////////////////////////////////////////////////////

vector<double> AnalysePairCorr::calc_pair_corr (const double * const *x,
						const double * const *y) {
  /* calculate pair correlation function */

  // set the maximum distance possible between two polymers
  /* A long standing problem here has to do with periodic boundary conditions 
  * Is the maximum distance allowed box size over 2 per direction,
  * how about the diagonal distance? 
  */ 
  
  int max_distance = static_cast<int>(sim.lx/2.+2);
  vector<int> cnt_distances(max_distance, 0);
  
  for (int step = 0; step < sim.nsteps; step++) {
    
    cout << "step / nsteps: " << step << " / " << sim.nsteps << endl;
    
    for (int i = 0; i < sim.npols-1; i++) {
      
      for (int j = i+1; j < sim.npols; j++) {
	
	// calculate the distance between the particles
	
	double dx = x[step][j] - x[step][i];
	dx = get_min_img_dist(dx, sim.lx);
	double dy = y[step][j] - y[step][i];
	dy = get_min_img_dist(dy, sim.ly);
	double dr = sqrt(dx*dx + dy*dy);
	
	// increase the histogram depending on the distance-based bin
	
	int bin = inearbyint(dr);
	if (bin < max_distance) {
	  cnt_distances[bin] += 2;
	}
	
      } // particle j
      
    } //particle i
    
  } // timesteps
  
  // normalize the data
  
  vector<double> gr(max_distance, 0.);
  double area = sim.lx*sim.ly;
  double number_density = sim.npols/area;
  double jfactor = 1./(2*pi*sim.npols*number_density);
  double delr = 1.;
  for (int j = 1; j < max_distance; j++) {
    double factor = jfactor/(delr*j);
    cnt_distances[j] /= sim.nsteps;
    gr[j] = factor*cnt_distances[j];
  }

  return gr;
}

/////////////////////////////////////////////////////////////////////////////////
