
/* analysis on intermediate scattering function */

/////////////////////////////////////////////////////////////////////////////////

#include "AnalyseInterScattering.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

AnalyseInterScattering::AnalyseInterScattering (const char *datafilename, char *forc) : 
    sim(datafilename, forc), polymers(datafilename, sim, forc) {
 
  // print information
  
  cout << "\nData is loaded successfully for the following file: \n" <<
    datafilename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "npols = " << sim.npols << endl;
  cout << "nbeads = " << sim.nbeads << endl;
   
  vector<double> result_1 = {0.};
  vector<double> result_2 = {0.};
  results = make_tuple(result_1, result_2);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

AnalyseInterScattering::~AnalyseInterScattering () { }

/////////////////////////////////////////////////////////////////////////////////

void AnalyseInterScattering::perform_analysis () {

  results = calc_inter_scattering(polymers.x, polymers.y);
  
  return;
}


/////////////////////////////////////////////////////////////////////////////////

void AnalyseInterScattering::write_analysis_results (const char *outfilepath) {

  vector<double> result_1;
  vector<double> result_2;
  tie(result_1, result_2) = results;  
  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  write_2d_analysis_data(result_1, result_2, outfilepath);
  
  return;
  
}

/////////////////////////////////////////////////////////////////////////////////

tuple<vector<double>, vector<double> > AnalyseInterScattering::calc_inter_scattering (const double * const *x,
										      const double * const *y) {
  /* calculate the self part/incoherent part of the intermediate scattering function */
   
  // populate the circular orientation of the wavevector

  int ndelay = sim.nsteps/2;
  vector<double> delay(ndelay, 0.);
  vector<double> Fs(ndelay, 0.);
  
  double kvector = 2*pi/(7.0*sim.sigma); 	
  int nks = 12;
  double kxs[nks], kys[nks];
  for (int j = 0; j < nks; j++) { 
    kxs[j] = cos(2*pi*j/nks)*kvector;
    kys[j] = sin(2*pi*j/nks)*kvector;
  }

  // allocate the arrays
  
  vector<double> avg_over_qvector(ndelay, 0.);
  
  // calculate the intermediate scattering function
  
  delay[0] = 0.;
  Fs[0] = 1.;
  
  for (int d = 1; d < ndelay; d++) {
        
      // sum over all wavevectors with a circular discretization
      
      for (int j = 0; j < nks; j++) {
	
	double term_to_avg = 0.;
	
	// sum over different time origins
	
	for (int step = 0; step < sim.nsteps-ndelay; step++) {
	  
	  double cost = 0.;

	  omp_set_num_threads(4);

	  #ifdef SUBTRACT_COM
	  // calculate center of mass drift
	  
	  double dcomx = 0.;
	  double dcomy = 0.;
	  #pragma omp parallel for reduction(+:dcomx,dcomy)
	  for (int p = 0; p < sim.npols; p++) {
	    double dx = x[d+step][p] - x[step][p];
	    double dy = y[d+step][p] - y[step][p];
	    dcomx += dx;
	    dcomy += dy;
	  }
	  
	  dcomx /= sim.npols;
	  dcomy /= sim.npols;
	  #endif
	  
	  // sum over particles
	  
	  #pragma omp parallel for reduction(+:cost)	  
	  for (int p = 0; p < sim.npols; p++) {
	    
	    double dx = x[d+step][p] - x[step][p];
	    double dy = y[d+step][p] - y[step][p];
	    #ifdef SUBTRACT_COM
	    dx -= dcomx;
	    dy -= dcomy;
	    #endif
	    double dotp = kxs[j]*dx + kys[j]*dy;
	    cost += cos(dotp);
	    
	  } // particles
	  
	  term_to_avg += cost;
	
	} // time origins
	
	term_to_avg /= (sim.npols*(sim.nsteps-ndelay));
	avg_over_qvector[d] += term_to_avg;
	
      } // wavectors
   
  } // delay
  
  // perform the normalizations
  
  for (int d = 1; d < ndelay; d++) {
    delay[d] = d*sim.dt; 
    Fs[d] = avg_over_qvector[d]/nks;
    
  }

  return make_tuple(delay, Fs);
}

/////////////////////////////////////////////////////////////////////////////////
