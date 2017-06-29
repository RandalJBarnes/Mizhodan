//=============================================================================
// write_results.cpp
//
//    Read in the observation data from the user-specified file.
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "engine.h"
#include "write_results.h"

//-----------------------------------------------------------------------------
void write_results( const std::string& resultsfilename, std::vector<ResultRecord> results ) {
   // Open the results file.
   std::ofstream resultsfile( resultsfilename );
   if ( resultsfile.fail() ) {
      std::stringstream message;
      message << "Could not open <" << resultsfilename << "> for output.";
      throw InvalidResultsFile(message.str());
   }

   // Write out the header line to the results file.
   resultsfile << "ID,X,Y,Zhat,Kstd" << std::endl;

   // Write out the results.
   resultsfile << std::setprecision(std::numeric_limits<long double>::digits10 + 1);

   for ( unsigned n = 0; n < results.size(); ++n ) {
      resultsfile << results[n].id << ',';
      resultsfile << results[n].x  << ',';
      resultsfile << results[n].y  << ',';
      resultsfile << results[n].zhat << ',';
      resultsfile << results[n].kstd;
      resultsfile << std::endl;
   }
   resultsfile.close();
}
