
/* analysis on vorticity field */

/////////////////////////////////////////////////////////////////////////////////

#include "AnalyseVorticity.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

AnalyseVorticity::AnalyseVorticity (const char *datafilename, char *forc) : 
    sim(datafilename, forc), polymers(datafilename, sim, forc) {
 
  // print information
  
  cout << "\nData is loaded successfully for the following file: \n" <<
    datafilename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "npols = " << sim.npols << endl;
  cout << "nbeads = " << sim.nbeads << endl;
  
  vector<double> wbin = {0.};
  vector<double> vxbin = {0.};
  vector<double> vybin = {0.};  
  int nbins = 0;
  int nvels = 0;
  double enstrophy = 0.;
  results = make_tuple(wbin, vxbin, vybin, nbins, nvels, enstrophy);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

AnalyseVorticity::~AnalyseVorticity () { }

/////////////////////////////////////////////////////////////////////////////////

void AnalyseVorticity::perform_analysis () {

  int delta = 4;                       // number of data points between two steps
				       // to calculate velocity
  int nvels = sim.nsteps-delta;        // number of data points in the velocity arrays  
  vector<double> vx;
  vector<double> vy;
  tie(vx, vy) = calc_velocity(polymers.x, polymers.y, delta, nvels);
   
  get_img_pos(polymers.x, polymers.y, 
	      sim.nsteps, sim.npols, sim.lx, sim.ly);  
	        
  double wbin = 16.; 			// bin width
  double woverlap = wbin*1./2.;		// bin overlap distance (currently %50)
  int nbins = int(sim.lx/wbin + 0.5);
  wbin = sim.lx/nbins;
  vector<double> vx_bin;
  vector<double> vy_bin;
  vector<vector<int> > ncells_per_bin;
  tie(vx_bin, vy_bin, ncells_per_bin) = calc_velocity_per_bin(vx, vy, polymers.x, polymers.y, 
						   	     wbin, nbins, woverlap,
							       sim.lx, sim.ly, sim.npols, nvels);
  vector<double> w_bin;
  double enstrophy;
  tie(w_bin, enstrophy) = calc_vorticity(ncells_per_bin, 
			   vx_bin, vy_bin, wbin, nbins, 
			   sim.lx, sim.ly, sim.npols, nvels);  
  results = make_tuple(w_bin, vx_bin, vy_bin, nbins, nvels, enstrophy);
  
  return;
}


/////////////////////////////////////////////////////////////////////////////////

void AnalyseVorticity::write_analysis_results (const char *outfilepath, 
					       const char *outfilepath_2) {

  vector<double> w_bin;
  vector<double> vx_bin;
  vector<double> vy_bin;
  int nbins;
  int nvels;
  double enstrophy;
  tie(w_bin, vx_bin, vy_bin, nbins, nvels, enstrophy) = results;  
  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath_2 << endl;    
  write_multid_analysis_data(w_bin, vx_bin, vy_bin, nbins, nvels, outfilepath);
  write_single_analysis_data(enstrophy, outfilepath_2);
  
  return;
  
}

/////////////////////////////////////////////////////////////////////////////////

tuple<vector<double>, vector<double> > AnalyseVorticity::calc_velocity (const double * const *x,
									const double * const *y,
									int delta, int nvels) {
  /* calculate velocities with dt as delta */

  // set variables related to the analysis
  
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
    
    for (int j = 0; j < ncells; j++) {
      vx[i*sim.npols+j] -= comvx;
      vy[i*sim.npols+j] -= comvy;
    } 	// cells loop
    #endif
    
  }     // velocity step loop
  
  return make_tuple(vx, vy);
}

/////////////////////////////////////////////////////////////////////////////////

void AnalyseVorticity::add_velocity_to_bins (int ix_bin, int iy_bin, int step,
			   vector<double> &vx_bin, vector<double> &vy_bin,
			   vector<int> &ncells_per_bin_per_step,
			   double vxc, double vyc, 
			   int nbins, int ncells) {

  int idx = step*nbins*nbins + ix_bin*nbins + iy_bin;

  vx_bin[idx] += vxc;
  vy_bin[idx] += vyc;
  ncells_per_bin_per_step[ix_bin*nbins+iy_bin]++;
  
  return; 
}

/////////////////////////////////////////////////////////////////////////////////

tuple<vector<double>, vector<double>, vector<vector<int> > > AnalyseVorticity::calc_velocity_per_bin (const vector<double> &vx, 
			    const vector<double> &vy, 
			    const double * const *x,
			    const double * const *y,
			    double wbin, int nbins, double woverlap,
			    double lx, double ly, 
			    int ncells, int nsteps) {
  /* calculate velocities for each bin for all timesteps */

  vector<double> vx_bin(nsteps*nbins*nbins, 0.0);
  vector<double> vy_bin(nsteps*nbins*nbins, 0.0);
  vector<vector<int> > ncells_per_bin(nsteps, vector<int>(nbins*nbins, 0.0));
  
  for (int step = 0; step < nsteps; step++) {
    
    vector<int> ncells_per_bin_per_step(nbins*nbins, 0.0);  

    // calculate the velocities 
    
    for (int j = 0; j < ncells; j++) {
      
      double vx_cell = vx[step*ncells+j];
      double vy_cell = vy[step*ncells+j];
      
      // add to the central bin
      
      int ix_bin = get_bin_number(x[step][j], wbin, nbins);
      int iy_bin = get_bin_number(y[step][j], wbin, nbins);

      add_velocity_to_bins(ix_bin, iy_bin, step, vx_bin, vy_bin,
			   ncells_per_bin_per_step, vx_cell, vy_cell, 
			   nbins, ncells);

      // add to the neighboring bins because of overlaps
      
      int xsub = fmod(x[step][j], wbin);
      int ysub = fmod(y[step][j], wbin);
      
      int ixneigh, iyneigh;
      
      if (xsub > woverlap) {
	ixneigh = (ix_bin+1) % nbins;
      }
      else {
	ixneigh = (ix_bin-1) % nbins;
      }
      if (ysub > woverlap) {
	iyneigh = (iy_bin+1) % nbins;
      }
      else {
	iyneigh = (iy_bin-1) % nbins;
      }
      
      if (ixneigh < 0) 
	ixneigh = nbins-1;
      if (iyneigh < 0)
	iyneigh = nbins-1;
            
      add_velocity_to_bins(ixneigh, iy_bin, step, vx_bin, vy_bin,
			   ncells_per_bin_per_step, vx_cell, vy_cell, nbins, ncells);      
      add_velocity_to_bins(ix_bin, iyneigh, step, vx_bin, vy_bin,
			   ncells_per_bin_per_step, vx_cell, vy_cell, nbins, ncells);      
      add_velocity_to_bins(ixneigh, iyneigh, step, vx_bin, vy_bin,
			   ncells_per_bin_per_step, vx_cell, vy_cell, nbins, ncells);      
            
    }	  // cells
            
    // normalize the velocities per bin
    
    for (int i = 0; i < nbins; i++) {
      for (int j = 0; j < nbins; j++) {
	int idx = step*nbins*nbins + i*nbins + j;
	if (ncells_per_bin_per_step[i*nbins+j] > 0) {
	  vx_bin[idx] /= ncells_per_bin_per_step[i*nbins+j];
	  vy_bin[idx] /= ncells_per_bin_per_step[i*nbins+j];	
	}
      }	 // ybins
    }	// xbins
    
    ncells_per_bin[step] = ncells_per_bin_per_step;
    
  }	  // timesteps
  
  return make_tuple(vx_bin, vy_bin, ncells_per_bin);
}

/////////////////////////////////////////////////////////////////////////////////

tuple<vector<double>, double> AnalyseVorticity::calc_vorticity (const vector<vector<int> > &ncells_per_bin,
								const vector<double> &vx_bin, 
								const vector<double> &vy_bin,
								double wbin, int nbins, 
								double lx, double ly, 
								int ncells, int nsteps) {
  /* calculate the vorticity for all timesteps with finite difference differentiation */
  
  double enstrophy = 0.0;
  vector<double> w_bin(nsteps*nbins*nbins, 0.0);
  
  for (int step = 0; step < nsteps; step++) {
    
    double ens_per_step = 0.0;
        
    for (int i = 0; i < nbins; i++) {
      for (int j = 0; j< nbins; j++) {

	if (ncells_per_bin[step][i*nbins+j] > 0) {
	  double rho = wbin/2./ncells_per_bin[step][i*nbins+j];
	  
	  int xfi = (i+1) % nbins;
	  int xbi = (i-1) % nbins;
	  if (xbi < 0) 
	    xbi = nbins-1;
	  double wx = vy_bin[step*nbins*nbins+xfi*nbins+j] - vy_bin[step*nbins*nbins+xbi*nbins+j];
	  
	  int yfi = (j+1) % nbins;
	  int ybi = (j-1) % nbins;
	  if (ybi < 0) 
	    ybi = nbins-1;
	  double wy = vx_bin[step*nbins*nbins+i*nbins+yfi] - vx_bin[step*nbins*nbins+i*nbins+ybi];
	  
	  w_bin[step*nbins*nbins+i*nbins+j] = (wx-wy)*rho;
	  
	  ens_per_step += (wx-wy)*(wx-wy)*rho*rho;
	}
	
      }		// ybins
    }		// xbins
    
    enstrophy += ens_per_step/nbins/nbins;
    
  }		// timesteps
  
  enstrophy /= (2.*nsteps);
  
  return make_tuple(w_bin, enstrophy);
}

/////////////////////////////////////////////////////////////////////////////////

