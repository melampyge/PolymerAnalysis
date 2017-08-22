
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
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

AnalyseVorticity::~AnalyseVorticity () { }

/////////////////////////////////////////////////////////////////////////////////

void AnalyseVorticity::perform_analysis () {

  int delta = 4;                       // number of data points between two steps
				                               // to calculate velocity
  int nvels = sim.nsteps-delta;        // number of data points in the velocity arrays  
   
  double wbin = 12.; 			                          // bin width
  int nbins = static_cast<int>(sim.lx/wbin + 0.5);  // number of bins 
  
  vector<vector<vector<double> > > vx_bin;
  vector<vector<vector<double> > > vy_bin;
  tie(vx_bin, vy_bin) = get_velocity_grid(
      polymers.x, polymers.y, nvels, delta,
      sim.dt, sim.lx, sim.ly, sim.npols, 
      wbin, nbins);
  
  vector<vector<vector<double> > > w_bin;
  double energy;
  double enstrophy;
  vector<double> energy_per_step;
  vector<double> enstrophy_per_step;
  tie(w_bin, energy, enstrophy, energy_per_step, 
      enstrophy_per_step) = calc_vorticity(
			   vx_bin, vy_bin, wbin, nbins, 
			   sim.lx, sim.ly, nvels);  
  results = make_tuple(w_bin, vx_bin, vy_bin, nbins, nvels, 
      energy, enstrophy, energy_per_step, enstrophy_per_step);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

void AnalyseVorticity::write_analysis_results (const char *outfilepath, 
					       const char *outfilepath_2, const char *outfilepath_3, 
                 const char *outfilepath_4) {

  vector<vector<vector<double> > > w_bin;
  vector<vector<vector<double> > > vx_bin;
  vector<vector<vector<double> > > vy_bin;
  int nbins;
  int nvels;
  double energy;
  double enstrophy;
  vector<double> energy_per_step;
  vector<double> enstrophy_per_step;
  tie(w_bin, vx_bin, vy_bin, nbins, nvels, 
      energy, enstrophy, energy_per_step, enstrophy_per_step) = results;  
  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath_2 << endl;    
  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath_3 << endl;    
  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath_4 << endl;    
  write_vorticity_analysis_data(w_bin, vx_bin, vy_bin, nbins, nvels, outfilepath);
  write_single_analysis_data(energy, outfilepath_2);
  write_single_analysis_data(enstrophy, outfilepath_3);
  write_2d_analysis_data(energy_per_step, enstrophy_per_step, outfilepath_4);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

tuple<vector<vector<vector<double> > >, double, double, vector<double>, vector<double> > AnalyseVorticity::calc_vorticity (
		const vector<vector<vector<double> > > &vx_bin, 
		const vector<vector<vector<double> > > &vy_bin,
		double wbin, int nbins, double lx, double ly, int nsteps) {
  /* calculate the vorticity for all timesteps with finite difference differentiation */

  cout << "Calculating the vorticity field" << endl;

  vector<double> energy_per_step(nsteps, 0.);
  vector<double> enstrophy_per_step(nsteps, 0.);
  double energy = 0.0;
  double enstrophy = 0.0;
  vector<vector<vector<double> > > w_bin(nsteps, vector<vector<double> >(nbins, vector<double>(nbins, 0.)));
  
  for (int step = 0; step < nsteps; step++) {

    cout << step << " / " << nsteps << endl;  

    double eng_per_step = 0.0;
    double ens_per_step = 0.0;
        
    for (int i = 0; i < nbins; i++) {
      for (int j = 0; j < nbins; j++) {

	      int xfi = (i+1) % nbins;
	      int xbi = (i-1) % nbins;
	      if (xbi < 0) 
	        xbi = nbins-1;
	      double wx = vy_bin[step][xfi][j] - vy_bin[step][xbi][j];
	  
	      int yfi = (j+1) % nbins;
	      int ybi = (j-1) % nbins;
	      if (ybi < 0) 
	        ybi = nbins-1;
	      double wy = vx_bin[step][i][yfi] - vx_bin[step][i][ybi];
	      
        w_bin[step][i][j] = wx-wy;
	      ens_per_step += (wx-wy)*(wx-wy);

        eng_per_step += vx_bin[step][i][j]*vx_bin[step][i][j] + 
          vy_bin[step][i][j]*vy_bin[step][i][j];
	
      }		// ybins
    }		  // xbins
    
    energy_per_step[step] = eng_per_step/(2*nbins*nbins);
    enstrophy_per_step[step] = ens_per_step/(2*nbins*nbins);
    energy += eng_per_step/(2*nbins*nbins);
    enstrophy += ens_per_step/(2*nbins*nbins);

  }		    // timesteps
  
  enstrophy /= nsteps;
  energy /= nsteps;
  for (int j = 0; j < nsteps; j++) {
    energy_per_step[j] /= nsteps;
    enstrophy_per_step[j] /= nsteps;
  }

  return make_tuple(w_bin, energy, enstrophy, energy_per_step, enstrophy_per_step);
}

/////////////////////////////////////////////////////////////////////////////////

