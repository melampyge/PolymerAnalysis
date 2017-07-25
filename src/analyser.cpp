
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
  
  char *datafilename = argv[1];
  char *beadsorpols = argv[2];
  char *filsorcells = argv[3];
  char *outfilename = argv[4];
  //char *outfilename_2 = argv[5];
  
  cout << "Data file: " << datafilename << endl;
  cout << "Analysis to be conducted on: " << beadsorpols << endl;
  cout << "Simulation type: " << filsorcells << endl;
  cout << "Analysis results output file path: " << outfilename << endl;
  //cout << "Analysis results output file path: " << outfilename_2 << endl;
  
  // perform the analysis
  
//   AnalyseNeighbours analysis(datafilename, filsorcells);
/*  AnalysePairCorr analysis(datafilename, filsorcells);*/
/*  AnalyseStaticStruct analysis(datafilename, filsorcells); */
/*  AnalyseLatticeOrder analysis(datafilename, filsorcells);*/
//   AnalyseAreaDiff analysis(datafilename, filsorcells);
  AnalyseDensityHistogram analysis(datafilename, filsorcells);

  analysis.perform_analysis();
  
  
  // write the computed results
  
  //analysis.write_analysis_results(outfilename, outfilename_2);
  analysis.write_analysis_results(outfilename);
  
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////

