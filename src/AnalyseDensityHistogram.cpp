
/* analysis on the density histogram */

/////////////////////////////////////////////////////////////////////////////////

#include "AnalyseDensityHistogram.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

AnalyseDensityHistogram::AnalyseDensityHistogram (const char *datafilename, char *forc) : 
    sim(datafilename, forc), beads(datafilename, sim) {
 
  // print information
  
  cout << "\nData is loaded successfully for the following file: \n" <<
    datafilename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "npols = " << sim.npols << endl;
  cout << "nbeads = " << sim.nbeads << endl;
  
  get_img_pos(beads.x, beads.y, 
	      sim.nsteps, sim.nbeads, sim.lx, sim.ly);    
  fill(results.begin(), results.end(), 0.);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

AnalyseDensityHistogram::~AnalyseDensityHistogram () { }

/////////////////////////////////////////////////////////////////////////////////

void AnalyseDensityHistogram::perform_analysis () {
  
  results = calc_densities(beads.x, beads.y);
  
  return;
}


/////////////////////////////////////////////////////////////////////////////////

void AnalyseDensityHistogram::write_analysis_results (const char *outfilepath) {

  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  write_1d_analysis_data(results, outfilepath);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

vector<double> AnalyseDensityHistogram::calc_densities (const double * const *x, 
							const double * const *y) {
  /* calculate densities per bin per frame */
  
  // set bin properties
  
  double avg_cell_diameter = 40.*sim.bond_length/(2.*pi);
  double bin_length = 2.5*avg_cell_diameter;
  int nbins = static_cast<int>(floor(sim.lx/bin_length));
  cout << "number of bins : " << nbins << endl;
  int nsize = nbins*nbins;
  bin_length = sim.lx/nbins;
  vector<double> densities(nsize, 0.);
  
  // calculate the density per bin per frame
  
  for (int step = 0; step < sim.nsteps; step++) {
    cout << "step / nsteps : " << step << " / " << sim.nsteps << endl;
    
    for (int j = 0; j < sim.nbeads; j++) {
      int xbin = get_bin_number(x[step][j], bin_length, nbins);
      int ybin = get_bin_number(y[step][j], bin_length, nbins);
      int bin = xbin*nbins + ybin;
      densities[bin]++;
    }
  }
  
  // normalization
  
  for (int j = 0; j < nsize; j++) 
    densities[j] /= (sim.nsteps*bin_length*bin_length);
  
  return densities;
}

/////////////////////////////////////////////////////////////////////////////////

