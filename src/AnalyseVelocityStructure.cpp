
/* analysis on velocity structure factors */

/////////////////////////////////////////////////////////////////////////////////

#include "AnalyseVelocityStructure.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

AnalyseVelocityStructure::AnalyseVelocityStructure (const char *datafilename, char *forc) : 
    sim(datafilename, forc), polymers(datafilename, sim, forc) {
 
  // print information
  
  cout << "\nData is loaded successfully for the following file: \n" <<
    datafilename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "npols = " << sim.npols << endl;
  cout << "nbeads = " << sim.nbeads << endl;
  
  get_img_pos(polymers.x, polymers.y, 
	      sim.nsteps, sim.npols, sim.lx, sim.ly);  
  fill(results.begin(), results.end(), 0.);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

AnalyseVelocityStructure::~AnalyseVelocityStructure () { }

/////////////////////////////////////////////////////////////////////////////////

void AnalyseVelocityStructure::perform_analysis () {


  
  return;
}


/////////////////////////////////////////////////////////////////////////////////

void AnalyseVelocityStructure::write_analysis_results (const char *outfilepath) {

  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  write_1d_analysis_data(results, outfilepath);
  
  return;
  
}

/////////////////////////////////////////////////////////////////////////////////

void AnalyseVelocityStructure::calc_vel_structure(const double * const *x,
						  const double * const *y) {
  
  int delta = 8;                     		   // number of data points between two steps
					       	   // to calculate velocity
  double delt = delta*sim.dt;
  int nvels = sim.nsteps-delta;       		   // number of data points in the velocity array
  int longest_dist = static_cast<int>(sim.lx+2);   // longest distance allowed by the sim. box
  
  for (int step = 0; step < nvels; step++) {
    for (int j = 0; j < sim.npols; j++) {
      double xj = x[step][j];
      double yj = y[step][j];
      double vxj = get_min_img_dist(x[step+delta][j]-xj)/delt;
      double vyj = get_min_img_dist(y[step+delta][j]-yj)/delt;
      
      for (int n = 0; sim.npols; n++) {
	double xn = x[step][n];
	double yn = y[step][n];
	double vxn = get_min_img_dist(x[step+delta][n]-xn)/delt;
	double vyn = get_min_img_dist(y[step+delta][n]-yn)/delt;
	
	double dvx = vxj-vxn;
	double dvy = vyj-vyn;
      
	double dx = xj-xn;
	double dy = yj-yn; 
	double dr = sqrt(dx*dx + dy*dy);	
	
	double dxpar = dx/dr;
	double dypar = dy/dr; 
	
	double par = dvx*dxpar + dvy*dypar;
	double perp = dvx*dxperp + dvy*dyperp;
	
      }    // second polymer
    }      // first polymer
  } 	   // timestep
  
  
}

/////////////////////////////////////////////////////////////////////////////////
