
/* main module to start a particular analysis */

/////////////////////////////////////////////////////////////////////////////////

#include "analyser.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {
  /* argv[1] -> datafilename         (path of the data file)
     argv[2] -> beadsorpols          (make a choice between reading beads or polymers)
     argv[3] -> filsorcells          (make a choice of sim type between filaments or cells)
     argv[4...] -> outfilename(s)... (path of the analysis results file)
   
   */
  
  // get command line arguments
  
  string datafilename = argv[1];
  string beadsorpols = argv[2];
  string filsorcells = argv[3];
  string outfilename = argv[4];
  string outfilename_2 = argv[5];
  
  // perform the analysis
  
  AnalyseNeighbours analysis(datafilename, filsorcells);
  analysis.perform_analysis();
  
  // write the computed results
  
  analysis.write_analysis_results(outfilename, outfilename_2);
  
  return 0;
}  

/////////////////////////////////////////////////////////////////////////////////

