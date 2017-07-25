
/* analysis on lattice order */

/////////////////////////////////////////////////////////////////////////////////

#include "AnalyseLatticeOrder.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

AnalyseLatticeOrder::AnalyseLatticeOrder (const char *datafilename, char *forc) : 
    sim(datafilename, forc), polymers(datafilename, sim, forc) {
 
  // print information
  
  cout << "\nData is loaded successfully for the following file: \n" <<
    datafilename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "npols = " << sim.npols << endl;
  cout << "nbeads = " << sim.nbeads << endl;
  
  get_img_pos(polymers.x, polymers.y, 
	      sim.nsteps, sim.npols, sim.lx, sim.ly);  
	      
  results = 0.;
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

AnalyseLatticeOrder::~AnalyseLatticeOrder () { }

/////////////////////////////////////////////////////////////////////////////////

void AnalyseLatticeOrder::perform_analysis () {
  
  int nn = 6;			// number of nearest neighbor points to investigate
  results = calc_glob_order(polymers.x, polymers.y, nn);  
  
  return;
}


/////////////////////////////////////////////////////////////////////////////////

void AnalyseLatticeOrder::write_analysis_results (const char *outfilepath) {

  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  write_single_analysis_data(results, outfilepath);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

void AnalyseLatticeOrder::populate_extended_pos_array (vector<Point<2> > &points, int npoints, int ncells, 
						      int step, const double * const *x, 
						      const double * const *y, double lx, double ly) {
  /* populate an extended array with all the image points */
    
  for (int i = 0; i < ncells; i++) {
    points[i].x[0] = x[step][i];
    points[i].x[1] = y[step][i];
    
    points[i+ncells].x[0] = x[step][i] - lx;
    points[i+ncells].x[1] = y[step][i] + ly;
    
    points[i+2*ncells].x[0] = x[step][i];
    points[i+2*ncells].x[1] = y[step][i] + ly;
    
    points[i+3*ncells].x[0] = x[step][i] + lx;
    points[i+3*ncells].x[1] = y[step][i] + ly; 
    
    points[i+4*ncells].x[0] = x[step][i] + lx;
    points[i+4*ncells].x[1] = y[step][i];
    
    points[i+5*ncells].x[0] = x[step][i] + lx;
    points[i+5*ncells].x[1] = y[step][i] - ly;  
    
    points[i+6*ncells].x[0] = x[step][i];
    points[i+6*ncells].x[1] = y[step][i] - ly; 
    
    points[i+7*ncells].x[0] = x[step][i] - lx;
    points[i+7*ncells].x[1] = y[step][i] - ly; 
    
    points[i+8*ncells].x[0] = x[step][i] - lx;
    points[i+8*ncells].x[1] = y[step][i];     
  }
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

double AnalyseLatticeOrder::calc_glob_order (const double * const *x, 
					     const double * const *y, 
					     int nn) {
  /* calculate global lattice order for nn */

  // set size of the extended array with all the periodic images
  
  const int npoints = 9*sim.npols;	
  double order_param = 0.; 
  
  for (int step = 0; step < sim.nsteps; step++) {
    
    cout << "step / nsteps: " << step << " / " << sim.nsteps << endl;
  
    // populate a vector of points in an extended way with all the image particles around the central box per time step

    vector<Point<2> > points(npoints);
    populate_extended_pos_array(points, npoints, sim.npols, step, 
				x, y, sim.lx, sim.ly);
    
    // build the kdtree 
    
    KDtree<2> kdtree(points);
    
    // calculate the global lattice order
    
    long double cost = 0.;
    long double sint = 0.;
    
    for (int j = 0; j < sim.npols; j++) {
      
      // find the nn number of nearest neighbors
      
      int neighbor_idx[nn];
      double neighbor_pos[nn];
      for (int ns = 0; ns < nn; ns++) {
	neighbor_idx[ns] = 0;  
	neighbor_pos[ns] = 0.;
      }
      kdtree.nnearest(j, neighbor_idx, neighbor_pos, nn);
      
      for (int k = 0; k < nn; k++) {
	
	double dx = dist_xy(points[neighbor_idx[k]], points[j], 0);
	double dy = dist_xy(points[neighbor_idx[k]], points[j], 1);

	double theta = nn*atan2(dy, dx);
	
	cost += cos(nn*theta);
	sint += sin(nn*theta);
	
      } // neighbors loop
      
    }  // cells loop
    
    cost /= (nn*sim.npols);
    sint /= (nn*sim.npols);
    
    order_param += cost*cost + sint*sint;
  
  }  // timestep loop

  order_param /= sim.nsteps;
  
  return order_param;
}

/////////////////////////////////////////////////////////////////////////////////
