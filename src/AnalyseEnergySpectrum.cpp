
/* analysis on kinetic energy spectrum */

/////////////////////////////////////////////////////////////////////////////////

#include "AnalyseEnergySpectrum.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

AnalyseEnergySpectrum::AnalyseEnergySpectrum (const char *datafilename, char *forc) : 
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

AnalyseEnergySpectrum::~AnalyseEnergySpectrum () { }

/////////////////////////////////////////////////////////////////////////////////

void AnalyseEnergySpectrum::perform_analysis () {

  int delta = 8;                     		   // number of data points between two steps
					       	                         // to calculate velocity
  int nvels = sim.nsteps-delta;       	   // number of data points in the velocity array
  int longest_dist = static_cast<int>(sim.lx+2);   // longest distance allowed by the sim. box
  vector<double> vx;
  vector<double> vy;
  tie(vx, vy) = calc_velocity(polymers.x, polymers.y, delta, nvels);
  vector<double> cvv = calc_sp_vel_corr(polymers.x, polymers.y, vx, vy, delta, longest_dist, nvels);
  results = calc_energy_spectrum(cvv, nvels, sim.npols);  
  //results = calc_energy_spectrum_direct(delta, nvels, sim.npols);  

  return;
}


/////////////////////////////////////////////////////////////////////////////////

void AnalyseEnergySpectrum::write_analysis_results (const char *outfilepath) {

  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  vector<double> result_1, result_2;
  tie(result_1, result_2) = results;
  write_2d_analysis_data(result_1, result_2, outfilepath);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

tuple<vector<double>, vector<double> > AnalyseEnergySpectrum::calc_velocity (const double * const *x,
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
      
    }   // polymers loop
    
    #ifdef SUBTRACT_COM
    comvx /= sim.npols;
    comvy /= sim.npols;    
    
    for (int j = 0; j < sim.npols; j++) {
      vx[i*sim.npols+j] -= comvx;
      vy[i*sim.npols+j] -= comvy;
    } 	// polymers loop
    #endif
    
  }     // velocity step loop
  
  return make_tuple(vx, vy);
}

/////////////////////////////////////////////////////////////////////////////////

vector<double> AnalyseEnergySpectrum::calc_sp_vel_corr (const double * const *x, 
                                      							    const double * const *y,
							                                          const vector<double> &vx,
							                                          const vector<double> &vy,
							                                          int delta, int ndata, int nsteps) {
  /* calculate spatial velocity correlation PER STEP with 2D bins */
  
  cout << "Calculating velocity correlations" << endl;
  
  // note that ndata here refers to the longest distance (longest_dist)
  // note that steps here refer to velocity data points
  
  vector<double> cvv(nsteps*ndata*ndata, 0);

  for (int step = 0; step < nsteps; step++) {
  
    vector<double> cvv_per_step(ndata*ndata, 0);
    vector<int> cnn_per_step(ndata*ndata, 0);
    
    for (int j1 = 0; j1 < sim.npols-1; j1++) {
      for (int j2 = j1+1; j2 < sim.npols; j2++) {
        
        // calculate the min. img. distance between the cells
        
        double dx = x[step][j2] - x[step][j1];
        dx = get_min_img_dist(dx, sim.lx);
        double dy = y[step][j2] - y[step][j1];
        dy = get_min_img_dist(dy, sim.ly);
        int ixbin = abs(inearbyint(dx));
        int iybin = abs(inearbyint(dy));
        
        if (ixbin < ndata && ixbin >= 0 && iybin < ndata && iybin >= 0) {
          cvv_per_step[ixbin*ndata+iybin] += vx[step*sim.npols+j1]*vx[step*sim.npols+j2] 
            + vy[step*sim.npols+j1]*vy[step*sim.npols+j2];
          cnn_per_step[ixbin*ndata+iybin] += 1;
        }
        else {
          cout << ixbin << "\t" << iybin << endl;
        }

      } // inner polymers loop
      
    }  // outer polymers loop
    
    for (int i = 0; i < ndata; i++) {
      for (int j = 0; j < ndata; j++) {
        if (cnn_per_step[i*ndata+j] != 0) {
          cvv[step*ndata*ndata+i*ndata+j] = cvv_per_step[i*ndata+j]/cnn_per_step[i*ndata+j];
        }
      }
    }
    
  }  // timestep loop
  
  return cvv;
}

/////////////////////////////////////////////////////////////////////////////////

tuple<vector<double>, vector<double> > AnalyseEnergySpectrum::calc_energy_spectrum(const vector<double> &cvv,
                                                                                   int nsteps, int natoms) {
  /* calculate the energy spectrum */
			       
  cout << "Calculating the energy spectrum" << endl;
  
  // note that ndata here refers to the longest distance (longest_dist)
  // note that steps here refer to velocity data points
  
  // set preliminary data
  
  double delk = 2*pi/sim.lx;
  int ndata = static_cast<int>(sim.lx);
  double lamda_min = 4.; 
  double lamda_max = sim.lx;
  double kmax = 2*pi/lamda_min;
  double kmin = 2*pi/lamda_max;
  int Nmax = 100;
   
  kmax = log(kmax);
  kmin = log(kmin);  
  double wbin = (kmax-kmin)/Nmax;

  vector<double> ek(Nmax, 0.);
  vector<double> kvec(Nmax);

  int nks = 12;
  vector<double> kxs(nks);
  vector<double> kys(nks);
  for (int j = 0; j < nks; j++) {
    kxs[j] = cos(2*pi*j/nks);
    kys[j] = sin(2*pi*j/nks);
  }

  cout << "delta k = " << delk << "\n" << "wbin = " << wbin << "\n" <<
    "kmin = " << exp(kmin) << "\n" << "kmax = " << exp(kmax) << "\n" << 
    "Nmax = " << Nmax << "\n" << endl;
  
  // populate log-spaced absolute value of the k vectors
  
  for (int j = 0; j < Nmax; j++) kvec[j] = kmin + j*wbin;
  for (int j = 0; j < Nmax; j++) kvec[j] = exp(kvec[j]);
  
  for (int step = 0; step < nsteps; step++) {
    
    for (int k = 0; k < Nmax; k++) {
      double cost = 0.0;
      double sint = 0.0;

      for (int n = 0; n < nks; n++) {
        
        double kxval = kvec[k]*kxs[n];
        double kyval = kvec[k]*kys[n];
        omp_set_num_threads(4);
        #pragma omp parallel for reduction(+:cost,sint)
        
        for (int i = 0; i < ndata; i++) {
          for (int j = 0; j < ndata; j++) {
	
	          double dotp = kxval*i + kyval*j;
            double cv_val = cvv[step*ndata*ndata+i*ndata+j];
	          cost += cos(dotp)*cv_val;    
	          sint += sin(dotp)*cv_val;

          } // y distance loop
	
        } // x distance loop

      } // circular k loop

      ek[k] += sqrt(cost*cost + sint*sint)/nks;
	
    } // absolute value k loop  
  
  } // timesteps
  
  // perform the normalization
  
  for (int j = 0; j < Nmax; j++)  { 
    ek[j] = j*ek[j]/(nsteps*2*pi); 
  }  

  return make_tuple(kvec, ek);
}

/////////////////////////////////////////////////////////////////////////////////

double * AnalyseEnergySpectrum::calc_velocity_per_step (const double * const *x,
										  const double * const *y,
										  int delta, int step) {
  /* calculate velocities with dt as delta */

  // set variables related to the analysis
  
  double deltaDt = delta*sim.dt;
  vector<double> vx(sim.npols, 0.);
  vector<double> vy(sim.npols, 0.);

  #ifdef SUBTRACT_COM
  long double comvx = 0.;
  long double comvy = 0.;
  #endif
  
  for (int j = 0; j < sim.npols; j++) {
    
    // note that UNWRAPPED COORDS ARE ASSUMED!
    
    double dx = x[step+delta][j] - x[step][j];
    vx[j] = dx/deltaDt;
    
    double dy = y[step+delta][j] - y[step][j];
    vy[j] = dy/deltaDt;
    
    #ifdef SUBTRACT_COM
    comvx += vx[j];
    comvy += vy[j];
    #endif
    
  }   // polymers loop
  
  #ifdef SUBTRACT_COM
  comvx /= sim.npols;
  comvy /= sim.npols;    
  
  for (int j = 0; j < sim.npols; j++) {
    vx[j] -= comvx;
    vy[j] -= comvy;
  } 	// polymers loop
  #endif
  
  double *in = (double *)malloc(sizeof(double)*sim.npols*sim.npols);
  for (int i = 0; i < sim.npols; i++) { 
    in[i] = vx[i];
    in[sim.npols+i] = vy[i];
  }

  return in;
}

/////////////////////////////////////////////////////////////////////////////////

tuple<vector<double>, vector<double> > AnalyseEnergySpectrum::calc_energy_spectrum_direct(int delta, int nsteps, int natoms) {
  /* calculate the energy spectrum by taking the autocorrelation of Fourier transformed velocities */
			       
  cout << "Calculating the energy spectrum directly in Fourier space" << endl;
  
  // note that ndata here refers to the longest distance (longest_dist)
  // note that steps here refer to velocity data points
  
  int nx = 2;
  int ny = (natoms%2 == 0 ? natoms/2 : natoms/2+1);
  int kmax = inearbyint(sqrt(2*pi*2*pi*ny*ny/(natoms*natoms)));
  cout << nx << "\t" << ny << "\t" << kmax << endl;

  vector<double> kvec(kmax, 0);
  vector<double> ek(kmax, 0);

  for (int i = 1; i < ny; i++) {
    double kx = 2*pi*i/sim.npols;

    for (int j = ny+1; j < 2*ny; j++) {
      double ky = 2*pi*(j-ny)/sim.npols;
      double k = sqrt(kx*kx + ky*ky);
      
      int kbin = inearbyint(k);
      kvec[kbin] = kbin;
    }
  }

  for (int step = 0; step < nsteps; step++) {
    double *in = calc_velocity_per_step(polymers.x, polymers.y, delta, step); 
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny);
    fftw_plan p = fftw_plan_dft_r2c_2d(nx, ny, in, out, FFTW_FORWARD);
   
    fftw_execute(p);
    
    vector<int> nkcnt(kmax, 0);
    vector<double> ek_per_step(kmax, 0);
    for (int i = 1; i < ny; i++) {
      double kx = 2*pi*i/sim.npols;
      double ukx2 =  out[i][0]*out[i][0] + out[i][1]*out[i][1];

      for (int j = ny+1; j < 2*ny; j++) {
        double ky = 2*pi*(j-ny)/sim.npols;
        double uky2 = out[j][0]*out[j][0] + out[j][1]*out[j][1];
        
        double k = sqrt(kx*kx + ky*ky);
        int kbin = inearbyint(k);

        nkcnt[kbin]++;
        ek_per_step[kbin] += ukx2 + uky2;
      }
    }

    for (int j = 0; j < kmax; j++) 
      ek[j] = ek_per_step[j]/nkcnt[j];

    fftw_destroy_plan(p);
    free(in); fftw_free(out);
  }

  for (int j = 0; j < kmax; j++)
    ek[j] = j*ek[j]/(2*nsteps);

  return make_tuple(kvec, ek);
}

/////////////////////////////////////////////////////////////////////////////////

