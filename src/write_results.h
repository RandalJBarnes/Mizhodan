//=============================================================================
// write_results.h
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#ifndef WRITE_RESULTS_H
#define WRITE_RESULTS_H

#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

//-----------------------------------------------------------------------------
class InvalidResultsFile : public std::runtime_error {
   public :
      InvalidResultsFile( const std::string& message ) : std::runtime_error(message) {
      }
};

//-----------------------------------------------------------------------------
void write_results( const std::string& outfilename, std::vector<ResultRecord> results );


//=============================================================================
#endif  // WRITE_RESULTS_H
