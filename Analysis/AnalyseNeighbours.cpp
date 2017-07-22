
/* analysis on neighbours
 such as number of neighbours */

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "AnalyseNeighbours.hpp"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

AnalyseNeighbours::AnalyseNeighbours (T sim_, Beads *beads_) {
  
  T sim = sim_;
  Beads *beads = beads_;
  
}

AnalyseNeighbours::~AnalyseNeighbours () { }

double AnalyseNeighbours::perform_analysis () {
  
  double avg_num_neigh = calc_num_neighbours(beads->x, beads->y);
  
  return avg_num_neigh;
}

void AnalyseNeighbours::write_analysis_results (double data, string outfilepath) {

  cout << "Writing the analysis results to the following file: \n" <<
    outfilepath << endl;
  write_single_analysis_data(data, outfilepath);
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void AnalyseNeighbours::build_linked_cell_list(const double * const *x,
                                               const double * const *y,
                                               vector<int> & heads,
                                               vector<int> & llist,
                                               double wbin, int nboxes, int nsize,
                                               int step, int natoms, double l) {
  /* build a linked cell/bucket list */
  
  for (int j = 0; j < natoms; j++) {
    int xbin = get_bin_number(get_single_img_pos(x[step][j], l), wbin, nboxes);
    int ybin = get_bin_number(get_single_img_pos(y[step][j], l), wbin, nboxes);
    int bin = xbin*nboxes + ybin;
    
    llist[j] = heads[bin];
    heads[bin] = j;
  }
  
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

double AnalyseNeighbours::calc_num_neighbours (const double * const *x,
                                               const double * const *y) {
  /* calculate densities per bin per frame */
  
  double avg_num_neigh = 0.;
  
  // set linked cell/box list parameters
  
  double rcut = 2.*sim.sigma;
  int nboxes = int(floor(sim.lx/rcut));
  int nsize = nboxes*nboxes;
  rcut = sim.lx/nboxes;
  vector<int> heads(nsize, -1);
  vector<int> llist(sim.nbeads*100, -1);
  
  vector<int> num_neighs_per_cell(sim.ncells, 0);
  
  for (int step = 0; step < sim.nsteps; step++) {
    cout << "step / nsteps : " << step << " / " << sim.nsteps << endl;
    
    // build the linked cell/box list
    
    build_linked_cell_list(x, y, heads, llist,
                           rcut, nboxes, nsize,
                           step, sim.nbeads, sim.lx);
    
    // count the number of neighbours of each cell
    
    vector<int> num_neighs_per_cell_per_time(sim.ncells, 0);
    int k = 0;
    for (int n = 0; n < sim.ncells; n++) {
      for (int j = 0; j < sim.nbpc[j]; j++) {
        int xbin = get_bin_number(get_single_img_pos(x[step][k], sim.lx), rcut, nboxes);
        int ybin = get_bin_number(get_single_img_pos(y[step][k], sim.ly), rcut, nboxes);
        //int bin = xbin*nboxes + ybin;
        
        // loop over the neighbouring bins
        
        for (int ix = -1; ix < 2; ix++) {
          int xneighbin = (ix+xbin)%nboxes;
          
          for (int iy = -1; iy < 2; iy++) {
            int yneighbin = (iy+ybin)%nboxes;
            int neighbin = xneighbin*nboxes + yneighbin;
            
            // restore the head bead
            
            
            
          }       // y neighbour bin loop
        }         // x neighbor bin loop
        k++;
      }           // nbpc loop
    }             // cell loop
    
  }               // timestep loop
  
  return avg_num_neigh;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

