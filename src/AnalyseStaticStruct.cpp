
/* analysis on static structure factor */

/////////////////////////////////////////////////////////////////////////////////

#include "AnalyseStaticStruct.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

AnalyseStaticStruct::AnalyseStaticStruct (const char *datafilename, char *forc) : 
    sim(datafilename, forc), polymers(datafilename, sim, forc) {
 
  // print information
  
  cout << "\nData is loaded successfully for the following file: \n" <<
    datafilename << endl;
  cout << "nsteps = " << sim.nsteps << endl;
  cout << "npols = " << sim.npols << endl;
  cout << "nbeads = " << sim.nbeads << endl;
  
  get_img_pos(polymers.x, polymers.y, 
	      sim.nsteps, sim.npols, sim.lx, sim.ly);  
  vector<double> result_1 = {0.};
  vector<double> result_2 = {0.};
  results = make_tuple(result_1, result_2);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////

AnalyseStaticStruct::~AnalyseStaticStruct () { }

/////////////////////////////////////////////////////////////////////////////////

void AnalyseStaticStruct::perform_analysis () {

  results = calc_static_struct(polymers.x, polymers.y);
  
  return;
}


/////////////////////////////////////////////////////////////////////////////////

void AnalyseStaticStruct::write_analysis_results (const char *outfilepath) {

  vector<double> result_1;
  vector<double> result_2;
  tie(result_1, result_2) = results;  
  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  write_2d_analysis_data(result_1, result_2, outfilepath);
  
  return;
  
}

/////////////////////////////////////////////////////////////////////////////////

tuple<vector<double>, vector<double> > AnalyseStaticStruct::calc_static_struct (const double * const *x,
										const double * const *y) {
  /* calculate static structure factor per timestep with a running average */

  // allocate the arrays and set preliminary information

  double delk = 2.*pi/sim.lx;
  int njump = 10;
  int ndata = static_cast<int>(sim.lx);
  vector<double> kxvec(ndata, 0.);
  vector<double> kyvec(ndata, 0.);
  for (int j = 0; j < ndata; j++) {
    kxvec[j] = delk*j;  
    kyvec[j] = delk*j; 
  }
  
  double mink = 0.;
  double maxk = ceil(sqrt(kxvec[ndata-1]*kxvec[ndata-1] + kyvec[ndata-1]*kyvec[ndata-1]));
  int Nmax = static_cast<int>(maxk);
  vector<double> Sk(Nmax, 0.);
  vector<double> kvec(Nmax, 0.);
  for (int j = 0; j < Nmax; j++) {
    kvec[j] = j;
  }
  
  vector<int> kcount(Nmax, 0);
  
  for (int step = 0; step < sim.nsteps; step += njump) {
    
    cout << "step / nsteps: " << step << " / " << sim.nsteps << endl;
    
    for (int kx = 0; kx < ndata; kx++) {
	double kxt = kxvec[kx]; 

	for (int ky = 0; ky < ndata; ky++) {
	  double costerm = 0.;
	  double sinterm = 0.;
	  double kyt = kyvec[ky];

	  omp_set_num_threads(4);
	  #pragma omp parallel for reduction(+:costerm,sinterm)
	  for (int j = 0; j < sim.npols; j++) {
	    double dotp = kxt*x[step][j] + kyt*y[step][j];
	    costerm += cos(dotp);
	    sinterm += sin(dotp);
    
	  } // particle loop

	  int knorm = static_cast<int>(sqrt(kxt*kxt + kyt*kyt));
	  kcount[knorm]++;
	  Sk[knorm] += costerm*costerm + sinterm*sinterm;

	} // ky loop
		
    } // kx loop 
  
  } // timesteps
  
  // perform the normalization
  
  double totalData = sim.nsteps/njump;
  for (int j = 0; j < Nmax; j++)  { 
    Sk[j] /= (sim.npols*kcount[j]); 
  }  
  
  return make_tuple(kvec, Sk);
}

/////////////////////////////////////////////////////////////////////////////////

