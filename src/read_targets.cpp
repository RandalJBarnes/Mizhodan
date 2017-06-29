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
#include "read_targets.h"

//-----------------------------------------------------------------------------
std::vector<TargetRecord> read_targets( const std::string& targetsfilename ) {
   std::vector<TargetRecord> targets;

   try {
      io::CSVReader<3,
         io::trim_chars<' ', '\t'>,
         io::no_quote_escape<','>,
         io::throw_on_overflow,
         io::single_and_empty_line_comment<'!','#'>> in(targetsfilename);

      std::string id;
      double X, Y;

      while (in.read_row(id,X,Y)){
         TargetRecord s = { id, X, Y };
         targets.push_back(s);
      }
   }
   catch (io::error::can_not_open_file& e) {
      std::stringstream message;
      message << "Could not open <" << targetsfilename << "> for input.";
      throw InvalidTargetsFile(message.str());
   }
   catch (...) {
      std::stringstream message;
      message << "Reading the target data failed on line " << targets.size()+1 << " of file " << targetsfilename << ".";
      throw InvalidTargetRecord(message.str());
   }

   return targets;
}
