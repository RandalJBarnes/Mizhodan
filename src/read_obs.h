//=============================================================================
// read_obs.h
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#ifndef read_obs_H
#define read_obs_H

#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

//-----------------------------------------------------------------------------
class InvalidObsFile : public std::runtime_error {
   public :
      InvalidObsFile( const std::string& message ) : std::runtime_error(message) {
      }
};

class InvalidObsRecord : public std::runtime_error {
   public :
      InvalidObsRecord( const std::string& message ) : std::runtime_error(message) {
      }
};

//-----------------------------------------------------------------------------
struct ObsRecord{
   std::string id;
   double x;
   double y;
   double z;
};

std::vector<ObsRecord> read_obs( const std::string& obsfilename );


//=============================================================================
#endif  // read_obs_H
