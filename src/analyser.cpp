
/* main module to start a particular analysis */

/////////////////////////////////////////////////////////////////////////////////

#include "analyser.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {
  /* argv[1] -> datafilename         (path of the data file)
     argv[2] -> analysisname 	     (name of the analysis type)
     argv[3] -> beadsorpols          (make a choice between reading beads or polymers)
     argv[4] -> filsorcells          (make a choice of sim type between filaments or cells)
     argv[5...] -> outfilename(s)... (path of the analysis results file)
   */
  
  // get command line arguments
  
  char *datafilename = argv[1];
  char *analysisname = argv[2];
  char *beadsorpols = argv[3];
  char *filsorcells = argv[4];
  char *outfilename = argv[5];
  
  cout << "Data file: " << datafilename << endl;
  cout << "Analysis type: " << analysisname << endl;
  cout << "Analysis to be conducted on: " << beadsorpols << endl;
  cout << "Simulation type: " << filsorcells << endl;
  cout << "Analysis results output file path: " << outfilename << endl;
  
  // choose analysis type
  
  if (strcmp(analysisname, "Num_neighbour") == 0) {   
    char *outfilename_2 = argv[6];  
    cout << "Analysis results output file path: " << outfilename_2 << endl;
    AnalyseNeighbours analysis(datafilename, filsorcells);
    analysis.perform_analysis();    
    analysis.write_analysis_results(outfilename, outfilename_2);
  }
  else if (strcmp(analysisname, "Vorticity") == 0) {   
    char *outfilename_2 = argv[6]; 
    char *outfilename_3 = argv[7];
    char *outfilename_4 = argv[8];
    cout << "Analysis results output file path: " << outfilename_2 << endl;
    cout << "Analysis results output file path: " << outfilename_3 << endl;
    cout << "Analysis results output file path: " << outfilename_4 << endl;
    AnalyseVorticity analysis(datafilename, filsorcells);
    analysis.perform_analysis();    
    analysis.write_analysis_results(outfilename, outfilename_2, 
        outfilename_3, outfilename_4);
  }
  else if (strcmp(analysisname, "Static_struct") == 0) { 	
    AnalyseStaticStruct analysis(datafilename, filsorcells);
    analysis.perform_analysis();  
    analysis.write_analysis_results(outfilename);    
  }
  else if (strcmp(analysisname, "Pair_corr") == 0) {
    AnalysePairCorr analysis(datafilename, filsorcells);
    analysis.perform_analysis();  
    analysis.write_analysis_results(outfilename);     
  }
  else if (strcmp(analysisname, "Lattice_order") == 0) {
    AnalyseLatticeOrder analysis(datafilename, filsorcells);
    analysis.perform_analysis();  
    analysis.write_analysis_results(outfilename);      
  }  
  else if (strcmp(analysisname, "Area_diff") == 0) {
    AnalyseAreaDiff analysis(datafilename, filsorcells);
    analysis.perform_analysis();  
    analysis.write_analysis_results(outfilename);      
  }  
  else if (strcmp(analysisname, "Density_histo") == 0) {
    AnalyseDensityHistogram analysis(datafilename, filsorcells);
    analysis.perform_analysis();  
    analysis.write_analysis_results(outfilename);      
  }   
  else if (strcmp(analysisname, "Inter_scattering_subt") == 0) {
    AnalyseInterScattering analysis(datafilename, filsorcells);
    analysis.perform_analysis();
    analysis.write_analysis_results(outfilename);
  }
  else if (strcmp(analysisname, "Sp_velocity_corr_subt") == 0) {
    AnalyseSpatialVelocityCorr analysis(datafilename, filsorcells);
    analysis.perform_analysis();
    analysis.write_analysis_results(outfilename);
  }
  else if (strcmp(analysisname, "Velocity_structure") == 0) {
    AnalyseVelocityStructure analysis(datafilename, filsorcells);
    analysis.perform_analysis();
    analysis.write_analysis_results(outfilename);
  }
  else if (strcmp(analysisname, "Energy_spectrum") == 0) {
    AnalyseEnergySpectrum analysis(datafilename, filsorcells);
    analysis.perform_analysis();
    analysis.write_analysis_results(outfilename);
  }

  return 0;
}

/////////////////////////////////////////////////////////////////////////////////

