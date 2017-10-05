
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
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

AnalyseVelocityStructure::~AnalyseVelocityStructure () { }

/////////////////////////////////////////////////////////////////////////////////

void AnalyseVelocityStructure::perform_analysis () {

  results = calc_vel_structure(polymers.x, polymers.y);
  
  return;
}


/////////////////////////////////////////////////////////////////////////////////

void AnalyseVelocityStructure::write_analysis_results (const char *outfilepath) {

  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  vector<double> v1, v2, v3, v4, v5, v6, v7, v8;
  tie(v1, v2, v3, v4, v5, v6, v7, v8) = results;
  write_8d_analysis_data(v1, v2, v3, v4, v5, v6, v7, v8, outfilepath);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

tuple8 AnalyseVelocityStructure::calc_vel_structure(const double * const *x,
		                                      				  const double * const *y) {
  
  int delta = 4;                     		   // number of data points between two steps
					       	                         // to calculate velocity
  double delt = delta*sim.dt;
  int nvels = sim.nsteps-delta;       		 // number of data points in the velocity array
  int longest_dist = static_cast<int>(sim.lx+2);  // longest distance allowed by the sim. box

  double avg_vsq = 0.;

  vector<double> cvpar(longest_dist, 0);
  vector<double> cvperp(longest_dist, 0);
  vector<double> cvpar_second(longest_dist, 0);
  vector<double> cvperp_second(longest_dist, 0);
  vector<double> cvpar_third(longest_dist, 0);
  vector<double> cvperp_third(longest_dist, 0);
  vector<double> cvpar_fourth(longest_dist, 0);
  vector<double> cvperp_fourth(longest_dist, 0);

  for (int step = 0; step < nvels; step++) {

    cout << step << " / " << nvels << endl;

    vector<int> cnn_per_step(longest_dist, 0);
    vector<double> cvpar_per_step(longest_dist, 0);
    vector<double> cvperp_per_step(longest_dist, 0);
    vector<double> cvpar_second_per_step(longest_dist, 0);
    vector<double> cvperp_second_per_step(longest_dist, 0);
    vector<double> cvpar_third_per_step(longest_dist, 0);
    vector<double> cvperp_third_per_step(longest_dist, 0);
    vector<double> cvpar_fourth_per_step(longest_dist, 0);
    vector<double> cvperp_fourth_per_step(longest_dist, 0);

    double avg_vsq_per_step = 0.;

    for (int j = 0; j < sim.npols-1; j++) {

      // position and velocity of first polymer

      double xj = x[step][j];
      double yj = y[step][j];
      double vxj = get_min_img_dist(x[step+delta][j]-xj, sim.lx)/delt;
      double vyj = get_min_img_dist(y[step+delta][j]-yj, sim.ly)/delt;
      avg_vsq_per_step += vxj*vxj + vyj*vyj; 

      for (int n = j+1; n < sim.npols; n++) {
        
        // position and velocity of second polymer

	      double xn = x[step][n];
	      double yn = y[step][n];
	      double vxn = get_min_img_dist(x[step+delta][n]-xn, sim.lx)/delt;
	      double vyn = get_min_img_dist(y[step+delta][n]-yn, sim.ly)/delt;
	      
        // parallel and perpendicular velocity increments
         
	      double dvx = vxn-vxj;
	      double dvy = vyn-vyj;
            
	      double dx = get_min_img_dist(xn-xj, sim.lx);
	      double dy = get_min_img_dist(yn-yj, sim.ly); 
	      double dr = sqrt(dx*dx + dy*dy);	
	      
	      double dxpar = dx/dr;
	      double dypar = dy/dr; 
	      double par = dvx*dxpar + dvy*dypar;
	      
	      double dxperp = -dy/dr;
	      double dyperp = dx/dr; 	
	      double perp = dvx*dxperp + dvy*dyperp;
      
        // moments of the velocity increments
        
        int rbin = inearbyint(dr);
        cnn_per_step[rbin] += 2;
        cvpar_per_step[rbin] += 2*par;
        cvperp_per_step[rbin] += 2*perp;
        cvpar_second_per_step[rbin] += 2*par*par;
        cvperp_second_per_step[rbin] += 2*perp*perp;
        cvpar_third_per_step[rbin] += 2*par*par*par;
        cvperp_third_per_step[rbin] += 2*perp*perp*perp;
        cvpar_fourth_per_step[rbin] += 2*par*par*par*par;
        cvperp_fourth_per_step[rbin] += 2*perp*perp*perp*perp;
      }   // second polymer
    }     // first polymer

    double vxj_last = get_min_img_dist(
        x[step+delta][sim.npols-1]-x[step][sim.npols-1], sim.lx)/delt;
    double vyj_last = get_min_img_dist(
        y[step+delta][sim.npols-1]-y[step][sim.npols-1], sim.ly)/delt;
    avg_vsq_per_step += vxj_last*vxj_last + vyj_last*vyj_last;
    avg_vsq += avg_vsq_per_step/sim.npols;
    
    // average the statistical moments of the velocity increments
    
    for (int k = 0; k < longest_dist; k++) {
      if (cnn_per_step[k] != 0) {
        cvpar[k] = cvpar_per_step[k]/cnn_per_step[k];
        cvperp[k] = cvperp_per_step[k]/cnn_per_step[k]; 
        cvpar_second[k] = cvpar_second_per_step[k]/cnn_per_step[k]; 
        cvperp_second[k] = cvperp_second_per_step[k]/cnn_per_step[k]; 
        cvpar_third[k] = cvpar_third_per_step[k]/cnn_per_step[k]; 
        cvperp_third[k] = cvperp_third_per_step[k]/cnn_per_step[k]; 
        cvpar_fourth[k] = cvpar_fourth_per_step[k]/cnn_per_step[k]; 
        cvperp_fourth[k] = cvperp_fourth_per_step[k]/cnn_per_step[k]; 
      }
    }

  } 	    // timestep
  
  // normalize the statistical moments
  
  double avg_v = pow(avg_vsq/nvels, 0.5);
  for (int k = 0; k < longest_dist; k++) {
    cvpar[k] /= (nvels*avg_v);
    cvperp[k] /= (nvels*avg_v);
    cvpar_second[k] /= (nvels*avg_v*avg_v);
    cvperp_second[k] /= (nvels*avg_v*avg_v);
    cvpar_third[k] /= (nvels*avg_v*avg_v*avg_v);
    cvperp_third[k] /= (nvels*avg_v*avg_v*avg_v);
    cvpar_fourth[k] /= (nvels*avg_v*avg_v*avg_v*avg_v);
    cvperp_fourth[k] /= (nvels*avg_v*avg_v*avg_v*avg_v); 
  }
 
  return make_tuple(cvpar, cvperp, cvpar_second, cvperp_second, 
      cvpar_third, cvperp_third, cvpar_fourth, cvperp_fourth);
}

/////////////////////////////////////////////////////////////////////////////////

