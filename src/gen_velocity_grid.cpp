
#include "gen_velocity_grid.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

void add_velocity_to_the_gridpnt(
    vector<vector<vector<double> > > &vxb,
    vector<vector<vector<double> > > &vyb, vector<vector<int> > &nb, 
    double vxp, double vyp, int i, int j, int k) {
  /* add velocity to the corresponding grid point */

  vxb[i][j][k] += vxp;
  vyb[i][j][k] += vyp;
  nb[j][k] += 1;

  return;
}

/////////////////////////////////////////////////////////////////////////////////

tuple<int, int> get_neighbour_bin_indices(int binidx, int nbins) {
  /* get the neighbouring bin indices */

  int forward_bin = (binidx+1) % nbins;
  int backward_bin = (binidx-1) % nbins;

  if (backward_bin < 0) backward_bin = nbins-1;
 
  return make_tuple(forward_bin, backward_bin);
}

/////////////////////////////////////////////////////////////////////////////////

tuple<int, int> determine_neighbour_bins(double xsub, double ysub,
    int nbins, double lower_bound, double upper_bound) {
  /* determine the neighboring bins to create overlaps between the bins */

  int ix_neigh_bin, iy_neigh_bin; 
  ix_neigh_bin = iy_neigh_bin = 0;

  if (xsub > upper_bound) {
    ix_neigh_bin = 1;  
  }
  else if (xsub < lower_bound) { 
    ix_neigh_bin = 2; 
  }
  
  if (ysub > upper_bound) {
    iy_neigh_bin = 1;  
  }
  else if (ysub < lower_bound) { 
    iy_neigh_bin = 2; 
  }
  
  return make_tuple(ix_neigh_bin, iy_neigh_bin);
}

/////////////////////////////////////////////////////////////////////////////////

tuple<vector<vector<vector<double> > >, vector<vector<vector<double> > > > get_velocity_grid (
    const double * const *x,
    const double * const *y, int nvels, int delta, 
    double dt, double lx, double ly, int npols,
    double wbin, int nbins) {
  /* calculate velocities of polymers on overlapping bins for all timesteps */

  cout << "Generating the velocity grid" << endl;

  double deltaDt = dt*delta;
  double woverlap = wbin*0.75;
  double wupper = woverlap;
  double wlower = wbin-woverlap;
  wbin = lx/nbins;
  vector<vector<vector<double> > > vx_bin(nvels, vector<vector<double> >(nbins, vector<double>(nbins, 0.)));
  vector<vector<vector<double> > > vy_bin(nvels, vector<vector<double> >(nbins, vector<double>(nbins, 0.)));
  
  for (int step = 0; step < nvels; step++) {
    
    cout << step << " / " << nvels << endl;

    vector<vector<int> > npols_bin(nbins, vector<int>(nbins, 0));
  
    for (int j = 0; j < npols; j++) {
      
      // calculate the velocity of the polymers 
      
      double vx_pol = get_min_img_dist(x[step+delta][j]-x[step][j], lx)/deltaDt; 
      double vy_pol = get_min_img_dist(y[step+delta][j]-y[step][j], ly)/deltaDt; 
      
      // get the central bin

      double x_pol_in_box = get_single_img_pos(x[step][j], lx);
      int ix_bin = get_bin_number(x_pol_in_box, wbin, nbins);
      double y_pol_in_box = get_single_img_pos(y[step][j], ly);
      int iy_bin =  get_bin_number(y_pol_in_box, wbin, nbins);
      
      // add the velocity to the central bin
      
      add_velocity_to_the_gridpnt(vx_bin, vy_bin, npols_bin, vx_pol, vy_pol, step, ix_bin, iy_bin);

      // get the neighbouring bin 
      // based on directions (up-down, left-right, middle) : 1-2-0

      int xsub = fmod(x_pol_in_box, wbin);
      int ysub = fmod(y_pol_in_box, wbin);

      int ix_neigh_bin, iy_neigh_bin;
      tie(ix_neigh_bin, iy_neigh_bin) = determine_neighbour_bins(xsub, ysub, 
          nbins, wlower, wupper); 
      
      int ix_forward_bin, ix_backward_bin; 
      tie(ix_forward_bin, ix_backward_bin) = get_neighbour_bin_indices(ix_bin, nbins);
      int iy_forward_bin, iy_backward_bin;
      tie(iy_forward_bin, iy_backward_bin) = get_neighbour_bin_indices(iy_bin, nbins);
      
      // add the velocity to the neighbouring/overlapping bins
      
      if (ix_neigh_bin == 1) {  // on left
          add_velocity_to_the_gridpnt(vx_bin, vy_bin, npols_bin, vx_pol, vy_pol, 
              step, ix_forward_bin, iy_bin);  
        if (iy_neigh_bin == 1) { // on left and up
          add_velocity_to_the_gridpnt(vx_bin, vy_bin, npols_bin, vx_pol, vy_pol, 
              step, ix_forward_bin, iy_forward_bin);  
        }
        else if (iy_neigh_bin == 2) { // on left and down
          add_velocity_to_the_gridpnt(vx_bin, vy_bin, npols_bin, vx_pol, vy_pol, 
              step, ix_forward_bin, iy_backward_bin);  
        }
      } // left is over

      if (ix_neigh_bin == 2) {  // on right
          add_velocity_to_the_gridpnt(vx_bin, vy_bin, npols_bin, vx_pol, vy_pol, 
              step, ix_backward_bin, iy_bin);  
        if (iy_neigh_bin == 1) { // on right and up
          add_velocity_to_the_gridpnt(vx_bin, vy_bin, npols_bin, vx_pol, vy_pol, 
              step, ix_backward_bin, iy_forward_bin);  
        }
        else if (iy_neigh_bin == 2) { // on right and down
          add_velocity_to_the_gridpnt(vx_bin, vy_bin, npols_bin, vx_pol, vy_pol, 
              step, ix_backward_bin, iy_backward_bin);  
        }
      } // right is over

    }   // polymers

    // normalize within this timestep

    for (int i = 0; i < nbins; i++) {
      for (int j = 0; j < nbins; j++) {
        if (npols_bin[i][j] != 0) {
          vx_bin[step][i][j] /= npols_bin[i][j];
          vy_bin[step][i][j] /= npols_bin[i][j];
        }
      }
    }
        
  }        // timestep 

  return make_tuple(vx_bin, vy_bin);
}

/////////////////////////////////////////////////////////////////////////////////

