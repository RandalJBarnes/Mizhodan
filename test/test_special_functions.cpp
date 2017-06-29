//=============================================================================
// test_special_functions.cpp
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
#include <cmath>

#include "test_special_functions.h"
#include "unit_test.h"
#include "..\src\special_functions.h"

//-----------------------------------------------------------------------------
// Hide all of the testing details inside an unnamed namespace. This allows me
// to create many small unit tests with polluting the global namespace.
//-----------------------------------------------------------------------------
namespace{
   const double TOLERANCE = 1e-6;

   //--------------------------------------------------------------------------
   // TestBeta
   //--------------------------------------------------------------------------
   bool TestBeta()
   {
      const double a[] = { 0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
      const double b[] = {10.5, 8.5, 6.5, 4.5, 2.5, 0.5};
      const int N = sizeof(a)/sizeof(double);

      // These test values were computed using Matlab's beta.
      const double y[] = {
         0.553539364153514,
         0.0342748832293199,
         0.00949150612504242,
         0.00766990393942822,
         0.021475731030399,
         0.773126317094364 };

      bool flag = true;

      for (int i = 0; i < N; ++i)
         flag &= isClose(Beta(a[i],b[i]), y[i], TOLERANCE);

      return flag;
   }

   //--------------------------------------------------------------------------
   // TestIncompleteBeta
   //--------------------------------------------------------------------------
   bool TestIncompleteBeta()
   {
      const double x[] = { 0.0, 0.1, 0.2, 0.3, 0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0};
      const double a[] = { 0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5};
      const double b[] = {10.5, 9.5, 8.5, 7.5, 6.5,  5.5,  4.5,  3.5,  2.5,  1.5,  0.5};
      const int N = sizeof(x)/sizeof(double);

      // These test values were computed using Matlab's betainc.
      const double y[] = {
         0.0,
         0.172547939103671,
         0.126378352888963,
         0.104744581678078,
         0.0964455758818249,
         0.0976021400100797,
         0.107992218900607,
         0.13064528349167,
         0.174376209213594,
         0.266917195063557,
         1.0 };

      bool flag = true;

      for (int i = 0; i < N; ++i)
         flag &= isClose(IncompleteBeta(x[i],a[i],b[i]), y[i], TOLERANCE);

      double yy = IncompleteBeta(0.954356616956718,16.5,2.5);
      flag &= isClose(yy, 0.9, TOLERANCE);

      return flag;
   }

   //--------------------------------------------------------------------------
   // TestIncompleteBetaInv
   //--------------------------------------------------------------------------
   bool TestIncompleteBetaInv()
   {
      const double p[] = { 0.1, 0.2, 0.3, 0.4, 0.5,  0.6,  0.7,  0.8,  0.9,  0.99, 0.999 };
      const double a[] = { 0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5 };
      const double b[] = {10.5, 9.5, 8.5, 7.5, 6.5,  5.5,  4.5,  3.5,  2.5,  1.5,  0.5 };
      const int N = sizeof(p)/sizeof(double);

      // These test values were computed using Matlab's betaincinv.
      const double x[] = {
         0.000769755304685,
         0.108267393149978,
         0.270734224119940,
         0.428332524440065,
         0.569705886332540,
         0.692808782635794,
         0.797922728147059,
         0.885438907619142,
         0.954356616956718,
         0.996942039627465,
         0.999999961217840 };

      bool flag = true;

      for (int i = 0; i < N; ++i)
         flag &= isClose(IncompleteBetaInv(p[i],a[i],b[i]), x[i], TOLERANCE);

      return flag;
   }

   //--------------------------------------------------------------------------
   // TestGamma
   //--------------------------------------------------------------------------
   bool TestGamma()
   {
      const double x[] = {-5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5 };
      const int N = sizeof(x)/sizeof(double);

      // These test values were computed using Matlab's gamma.
      const double y[] = {
         0.0109126547819099,
         -0.0600196013005042,
         0.270088205852269,
         -0.945308720482942,
         2.36327180120735,
         -3.54490770181103,
         1.77245385090552,
         1.0,
         0.886226925452758,
         1.0,
         1.32934038817914,
         2.0,
         3.32335097044784,
         6.0,
         11.6317283965675,
         24.0,
         52.3427777845535 };

      bool flag = true;

      for (int i = 0; i < N; ++i)
         flag &= isClose(Gamma(x[i]), y[i], TOLERANCE);

      return flag;
   }

   //--------------------------------------------------------------------------
   // TestIncompleteGamma
   //--------------------------------------------------------------------------
   bool TestIncompleteGamma()
   {
      const double x[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
      const double a[] = {10, 9, 8, 7, 6, 5, 4, 3, 2,  1};
      const int N = sizeof(x)/sizeof(double);

      // These test values were computed using Matlab's gammainc.
      const double y[] = {
         0.000000111425478,
         0.000237447328261,
         0.011904503856357,
         0.110673978402574,
         0.384039345166937,
         0.714943499683369,
         0.918234583755278,
         0.986246032255997,
         0.998765901959133,
         0.999954600070238 };

      bool flag = true;

      for (int i = 0; i < N; ++i)
         flag &= isClose(IncompleteGamma(x[i],a[i]), y[i], TOLERANCE);

      return flag;
   }

   //--------------------------------------------------------------------------
   // TestIncompleteGammaInv
   //--------------------------------------------------------------------------
   bool TestIncompleteGammaInv()
   {
      const double p[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99};
      const double a[] = {8.5, 7.5, 6.5, 5.5, 4.5, 3.5, 2.5, 1.5, 0.5, 0.01};
      const int N = sizeof(p)/sizeof(double);

      // These test values were computed using Matlab's gammaincinv.
      const double x[] = {
         5.042593167309668,
         5.153479503312642,
         4.962841207473447,
         4.618642711920757,
         4.171416346126477,
         3.641603816420171,
         3.032214992077453,
         2.320813838043723,
         1.352771727047708,
         0.265052550251590 };

      bool flag = true;

      for (int i = 0; i < N; ++i)
         flag &= isClose(IncompleteGammaInv(p[i],a[i]), x[i], TOLERANCE);

      return flag;
   }


   //--------------------------------------------------------------------------
   // TestGaussianCDF
   //--------------------------------------------------------------------------
   bool TestGaussianCDF()
   {
      const double x[] = {-4, -3, -2, -1, 0, 1, 2, 3, 4};
      const int N = sizeof(x)/sizeof(double);

      // These test values we computed using MATLAB's normcdf.
      const double y[] = {
         3.167124183312e-005,
         0.0013498980316301,
         0.0227501319481792,
         0.158655253931457,
         0.5,
         0.841344746068543,
         0.977249868051821,
         0.99865010196837,
         0.999968328758167 };

      bool flag = true;

      for (int i = 0; i < N; ++i)
         flag &= isClose(GaussianCDF(x[i]), y[i], TOLERANCE);

      return flag;
   }

   //--------------------------------------------------------------------------
   // TestGaussianCDFInv
   //--------------------------------------------------------------------------
   bool TestGaussianCDFInv()
   {
      const double p[] = {0.0001, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9};
      const int N = sizeof(p)/sizeof(double);

      // The test values we computed using MATLAB's normcdfinv.
      const double z[] = {
         -3.71901648545568,
         -3.09023230616781,
         -2.32634787404084,
         -1.2815515655446,
         -0.841621233572914,
         -0.524400512708041,
         -0.2533471031358,
         0.2533471031358,
         0.524400512708041,
         0.841621233572914,
         1.2815515655446 };

      bool flag = true;

      for (int i = 0; i < N; ++i)
         flag &= isClose(GaussianCDFInv(p[i]), z[i], TOLERANCE);

      return flag;
   }
}

//-----------------------------------------------------------------------------
// test_SpecialFunctions
//-----------------------------------------------------------------------------
std::pair<int,int> test_SpecialFunctions()
{
   int nsucc = 0;
   int nfail = 0;

   TALLY( TestBeta() );
   TALLY( TestIncompleteBeta() );
   TALLY( TestIncompleteBetaInv() );
   TALLY( TestGamma() );
   TALLY( TestIncompleteGamma() );
   TALLY( TestIncompleteGammaInv() );
   TALLY( TestGaussianCDF() );
   TALLY( TestGaussianCDFInv() );

   return std::make_pair( nsucc, nfail );
}
