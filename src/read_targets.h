//=============================================================================
// read_targets.h
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#ifndef read_targets_H
#define read_targets_H

#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

//-----------------------------------------------------------------------------
class InvalidTargetsFile : public std::runtime_error {
   public :
      InvalidTargetsFile( const std::string& message ) : std::runtime_error(message) {
      }
};

class InvalidTargetRecord : public std::runtime_error {
   public :
      InvalidTargetRecord( const std::string& message ) : std::runtime_error(message) {
      }
};

//-----------------------------------------------------------------------------
struct TargetRecord{
   std::string id;
   double x;
   double y;
};

std::vector<TargetRecord> read_targets( const std::string& targetsfilename );


//=============================================================================
#endif  // read_targets_H
