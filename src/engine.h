//=============================================================================
// engine.h
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#ifndef ENGINE_H
#define ENGINE_H

#include <stdexcept>
#include <vector>

#include "read_obs.h"
#include "read_targets.h"

//-----------------------------------------------------------------------------
class NoTargetsSpecified : public std::runtime_error {
   public :
      NoTargetsSpecified( const std::string& message ) : std::runtime_error(message) {
      }
};

class TooFewObservations : public std::runtime_error {
   public :
      TooFewObservations( const std::string& message ) : std::runtime_error(message) {
      }
};

class TooManyObservations : public std::runtime_error {
   public :
      TooManyObservations( const std::string& message ) : std::runtime_error(message) {
      }
};

class CholeskyDecompositionFailed : public std::runtime_error {
   public :
      CholeskyDecompositionFailed( const std::string& message ) : std::runtime_error(message) {
      }
};


//-----------------------------------------------------------------------------
struct ResultRecord {
   std::string id;
   double x;
   double y;
   double zhat;
   double kstd;
};

//-----------------------------------------------------------------------------
std::vector<ResultRecord> Engine(
   double nugget,
   double sill,
   double range,
   std::vector<ObsRecord> obs,
   std::vector<TargetRecord> targets
);


//=============================================================================
#endif  // ENGINE_H
