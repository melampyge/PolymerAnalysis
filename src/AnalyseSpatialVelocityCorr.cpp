
/* analysis on spatial velocity correlation */

/////////////////////////////////////////////////////////////////////////////////

#include "AnalyseSpatialVelocityCorr.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

AnalyseSpatialVelocityCorr::AnalyseSpatialVelocityCorr (const char *datafilename, char *forc) : 
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

AnalyseSpatialVelocityCorr::~AnalyseSpatialVelocityCorr () { }

/////////////////////////////////////////////////////////////////////////////////

void AnalyseSpatialVelocityCorr::perform_analysis () {

  int delta = 4;                     		   // number of data points between two steps
					       	   // to calculate velocity
  int nvels = sim.nsteps-delta;       		   // number of data points in the velocity array
  int longest_dist = static_cast<int>(sim.lx+2);   // longest distance allowed by the sim. box
  vector<double> vx;
  vector<double> vy;
  tie(vx, vy) = calc_velocity(polymers.x, polymers.y, delta, nvels);
  results = calc_sp_vel_corr(polymers.x, polymers.y, vx, vy, delta, longest_dist, nvels);
  
  return;
}


/////////////////////////////////////////////////////////////////////////////////

void AnalyseSpatialVelocityCorr::write_analysis_results (const char *outfilepath) {

  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  write_1d_analysis_data(results, outfilepath);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

tuple<vector<double>, vector<double> > AnalyseSpatialVelocityCorr::calc_velocity (const double * const *x,
										  const double * const *y,
										  int delta, int nvels) {
  /* calculate velocities with dt as delta */

  // set variables related to the analysis
  
  cout << "Calculating velocities" << endl;

  vector<double> vx(nvels*sim.npols, 0.);
  vector<double> vy(nvels*sim.npols, 0.);
  
  double deltaDt = delta*sim.dt;
  for (int i = 0; i < nvels; i++) {
    
    #ifdef SUBTRACT_COM
    long double comvx = 0.;
    long double comvy = 0.;
    #endif
    
    for (int j = 0; j < sim.npols; j++) {
      
      // note that UNWRAPPED COORDS ARE ASSUMED!
      
      double dx = x[i+delta][j] - x[i][j];
      vx[i*sim.npols+j] = dx/deltaDt;
      
      double dy = y[i+delta][j] - y[i][j];
      vy[i*sim.npols+j] = dy/deltaDt;
      
      #ifdef SUBTRACT_COM
      comvx += vx[i*sim.npols+j];
      comvy += vy[i*sim.npols+j];
      #endif
      
    }   // cell loop
    
    #ifdef SUBTRACT_COM
    comvx /= sim.npols;
    comvy /= sim.npols;    
    
    for (int j = 0; j < sim.npols; j++) {
      vx[i*sim.npols+j] -= comvx;
      vy[i*sim.npols+j] -= comvy;
    } 	// cells loop
    #endif
    
  }     // velocity step loop
  
  return make_tuple(vx, vy);
}

/////////////////////////////////////////////////////////////////////////////////

vector<double> AnalyseSpatialVelocityCorr::calc_sp_vel_corr (const double * const *x, 
							    const double * const *y,
							    const vector<double> &vx,
							    const vector<double> &vy,
							    int delta, int ndata, int nvels) {
  /* calculate spatial velocity correlation
  with the following definition (Wysocki, et. al.):
  Cvv(dr) = <v_i(r)*v_j(r+dr)>/(<v_i(r)^2/N>*(del_ri*del_rj))
  */
 
  cout << "Calculating velocity correlation" << endl;

  // note that ndata here refers to the longest distance (longest_dist)
  // note that steps here refer to velocity data points
  // velocity time unit and position time unit are different
  
  vector<double> cvv(ndata, 0.);
  
  for (int step = 0; step < nvels; step++) {
    
    cout << step << "\t" << nvels << endl;

    vector<double> cnorm_per_step(ndata, 0.);
    vector<double> cvv_per_step(ndata, 0.);
    vector<int> cnn_per_step(ndata, 0.);
    
    for (int j1 = 0; j1 < sim.npols-1; j1++) {
      for (int j2 = j1+1; j2 < sim.npols; j2++) {
        
        // calculate the min. img. distance between the cells
        
        double dx = x[step][j2] - x[step][j1];
        dx = get_min_img_dist(dx, sim.lx);
        double dy = y[step][j2] - y[step][j1];
        dy = get_min_img_dist(dy, sim.ly);
        double dr = sqrt(dx*dx + dy*dy);
        int ibin = inearbyint(dr);
        
        cvv_per_step[ibin] += 2*(vx[step*sim.npols+j1]*vx[step*sim.npols+j2] + 
            vy[step*sim.npols+j1]*vy[step*sim.npols+j2]);
	      cnorm_per_step[ibin] += vx[step*sim.npols+j1]*vx[step*sim.npols+j1] + 
          vy[step*sim.npols+j1]*vy[step*sim.npols+j1] + 
          vx[step*sim.npols+j2]*vx[step*sim.npols+j2] + 
          vy[step*sim.npols+j2]*vy[step*sim.npols+j2];	
        cnn_per_step[ibin] += 1;

      } // inner cells loop
      
    }  // outer cells loop
    
    
    for (int i = 0; i < ndata; i++) {
      if (cnn_per_step[i] != 0 && cnorm_per_step[i] != 0.) {
	      cvv[i] += cvv_per_step[i]/cnorm_per_step[i];
      }
    }
    
  }  // timestep loop
  
  // normalize
  
  for (int j = 0; j < ndata; j++) {
    cvv[j] /= nvels;
  }
  
  return cvv;
}

/////////////////////////////////////////////////////////////////////////////////

