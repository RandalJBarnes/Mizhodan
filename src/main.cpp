//=============================================================================
// main.cpp
//
//    The driver for the command line version of the user interface.
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#include <cstring>
#include <ctime>
#include <iostream>

#include "engine.h"
#include "now.h"
#include "numerical_constants.h"
#include "read_obs.h"
#include "read_targets.h"
#include "version.h"
#include "write_results.h"


//-----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
   // Check the command line.
   switch (argc) {
      case 1: {
         Usage();
         return 0;
      }
      case 2: {
         if ( strcmp(argv[1], "--help") == 0 )
            Help();
         else if ( strcmp(argv[1], "--version") == 0 )
            Version();
         else
            Usage();
         return 0;
      }
      case 7: {
         Banner( std::cout );
         break;
      }
      default: {
         Usage();
         return 1;
      }
   }

   // Get and check the semi-variogram nugget effect.
   double nugget = atof( argv[1] );
   if ( nugget <= EPS ) {
      std::cerr << "ERROR: nugget = " << argv[1] << " is not valid;  0 < nugget." << std::endl;
      std::cerr << std::endl;
      Usage();
      return 2;
   }

   // Get and check the semi-variogram sill.
   double sill = atof( argv[2] );
   if ( sill <= EPS ) {
      std::cerr << "ERROR: sill = " << argv[2] << " is not valid;  0 < sill." << std::endl;
      std::cerr << std::endl;
      Usage();
      return 2;
   }

   // Get and check the semi-variogram range.
   double range = atof( argv[3] );
   if ( range <= EPS ) {
      std::cerr << "ERROR: range = " << argv[3] << " is not valid;  0 < range." << std::endl;
      std::cerr << std::endl;
      Usage();
      return 2;
   }

   // Read in the observation data from the specified file.
   std::vector<ObsRecord> obs;

   try {
      obs = read_obs( argv[4] );
      std::cout << obs.size() << " data records read from <" << argv[4] << ">." << std::endl;
   }
   catch (InvalidObsFile& e) {
      std::cerr << e.what() << std::endl;
      return 3;
   }
   catch (InvalidObsRecord& e) {
      std::cerr << e.what() << std::endl;
      return 3;
   }

   // Read in the target data from the specified input data file.
   std::vector<TargetRecord> targets;

   try {
      targets = read_targets( argv[5] );
      std::cout << targets.size() << " target locations read from <" << argv[5] << ">." << std::endl;
   }
   catch (InvalidTargetsFile& e) {
      std::cerr << e.what() << std::endl;
      return 3;
   }
   catch (InvalidTargetRecord& e) {
      std::cerr << e.what() << std::endl;
      return 3;
   }

   // Execute all of the computations.
   std::vector<ResultRecord> results;
   try {
       results = Engine(nugget, sill, range, obs, targets);
   }
   catch (NoTargetsSpecified& e) {
      std::cerr << e.what() << std::endl;
      return 4;
   }
   catch (CholeskyDecompositionFailed& e) {
      std::cerr << e.what() << std::endl;
      return 4;
   }
   catch (...) {
      std::cerr << "The Mizhodan Engine failed for an unknown reason." << std::endl;
      throw;
   }

   // Write out the results to the specified output data file.
   try {
      write_results( argv[6], results );
      std::cout << "Results file <" << argv[6] << "> created. " << std::endl;
   }
   catch (InvalidResultsFile& e) {
      std::cerr << e.what() << std::endl;
      return 5;
   }

   // Successful termination.
   double elapsed = static_cast<double>(clock())/CLOCKS_PER_SEC;
   std::cout << "elapsed time: " << std::fixed << elapsed << " seconds." << std::endl;
   std::cout << std::endl;

   // Terminate execution.
	return 0;
}
