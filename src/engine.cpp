//=============================================================================
// engine.cpp
//
//    Compute the boomerang statistic for each measured location using the
//    user-specified semi-variogram model.
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#include <cassert>
#include <iomanip>
#include <math.h>
#include <numeric>
#include <sstream>

#include "engine.h"
#include "matrix.h"
#include "linear_systems.h"

//=============================================================================
//
//=============================================================================
std::vector<ResultRecord> Engine(
   double nugget,
   double sill,
   double range,
   std::vector<ObsRecord> obs,
   std::vector<TargetRecord> targets )
{
   // Manifest constants.
   const int MINIMUM_COUNT = 10;
   const int MAXIMUM_COUNT = 500;

   const int M = targets.size();
   if (M < 1) {
      throw NoTargetsSpecified("No targets were specified.");
   }

   const int N = obs.size();

   if (N < MINIMUM_COUNT) {
      std::stringstream message;
      message << "There must be at least " << MINIMUM_COUNT << " observations.";
      throw TooFewObservations(message.str());
   }

   if (N > MAXIMUM_COUNT) {
      std::stringstream message;
      message << "There must be no more than " << MAXIMUM_COUNT << " observations.";
      throw TooManyObservations(message.str());
   }

   // Create the matrix of observed values.
   Matrix Z(N, 1);
   for (int n = 0; n < N; ++n)
      Z(n,0) = obs[n].z;

   // Create the covariance matrix for all of the observations.
   Matrix C(N, N, sill);
   for (int i = 0; i < N-1; ++i) {
      for (int j = i+1; j < N; ++j) {
         double h = hypot( obs[i].x-obs[j].x, obs[i].y-obs[j].y );
         C(i,j) = (sill-nugget)*exp(-3.0*h/range);
         C(j,i) = C(i,j);
      }
   }

   // Solve the Ordinary Kriging system.
   Matrix L;
   if (!CholeskyDecomposition(C,L)) {
      throw CholeskyDecompositionFailed("Cholesky decomposition of the Kriging system failed.");
   }

   // Precompute the v matrix.
   Matrix ones(N, 1, 1.0);
   Matrix v;
   CholeskySolve(L,ones,v);
   double sumv = Sum(v);

   // Pass through the set of observations one at a time.
   std::vector<ResultRecord> results(M);

   for (int m = 0; m < M; ++m) {
      // Setup the Ordinary Kriging right-hand-side.
      Matrix b(N,1);
      for (int n = 0; n < N; ++n) {
         double h = hypot(targets[m].x - obs[n].x, targets[m].y - obs[n].y);
         b(n,0) = (sill - nugget) * exp(-3.0 * h / range);
      }

      // Solve the Ordinary Kriging system.
      Matrix u;
      CholeskySolve(L,b,u);

      double lambda = (Sum(u) - 1) / sumv;

      Matrix lv;
      Multiply_aM(lambda, v, lv);

      Matrix w;
      Subtract_MM(u, lv, w);

      double zhat = DotProduct(w, Z);
      double kstd = sqrt( sill - DotProduct(b, w) - lambda );

      results[m].id   = targets[m].id;
      results[m].x    = targets[m].x;
      results[m].y    = targets[m].y;
      results[m].zhat = zhat;
      results[m].kstd = kstd;
   }
   return results;
}
