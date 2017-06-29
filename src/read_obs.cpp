//=============================================================================
// read_obs.cpp
//
//    Read in the observation data from the user-specified file.
//
// notes:
// o  This function uses Ben Strasser's "fast-cpp-csv-parser" to read in the
//    .csv input file. See
//
//       https://github.com/ben-strasser/fast-cpp-csv-parser
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
#include <sstream>

#include "../include/csv.h"
#include "read_obs.h"

//-----------------------------------------------------------------------------
std::vector<ObsRecord> read_obs( const std::string& obsfilename ) {
   std::vector<ObsRecord> obs;

   try {
      io::CSVReader<4,
         io::trim_chars<' ', '\t'>,
         io::no_quote_escape<','>,
         io::throw_on_overflow,
         io::single_and_empty_line_comment<'!','#'>> in(obsfilename);

      std::string id;
      double X, Y, Z;

      while (in.read_row(id,X,Y,Z)){
         ObsRecord s = { id, X, Y, Z };
         obs.push_back(s);
      }
   }
   catch (io::error::can_not_open_file& e) {
      std::stringstream message;
      message << "Could not open <" << obsfilename << "> for input.";
      throw InvalidObsFile(message.str());
   }
   catch (...) {
      std::stringstream message;
      message << "Reading the observation data failed on line " << obs.size()+1 << " of file " << obsfilename << ".";
      throw InvalidObsRecord(message.str());
   }

   return obs;
}
