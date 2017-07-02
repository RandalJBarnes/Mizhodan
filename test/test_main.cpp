//=============================================================================
// test_main.cpp
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#include <iostream>

#include "test_engine.h"
#include "test_linear_systems.h"
#include "test_matrix.h"
#include "test_special_functions.h"

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
int main()
{
   int nsucc = 0;
   int nfail = 0;

   std::pair<int,int> counts;

   counts = test_Engine();
   nsucc += counts.first;
   nfail += counts.second;

   counts = test_LinearSystems();
   nsucc += counts.first;
   nfail += counts.second;

   counts = test_Matrix();
   nsucc += counts.first;
   nfail += counts.second;

   counts = test_SpecialFunctions();
   nsucc += counts.first;
   nfail += counts.second;

   if (nfail > 0)
      std::cerr << "MIZHODAN TESTS: nsucc = " << nsucc << '\t' << "nfail = " << nfail << std::endl;
   else
      std::cerr << "MIZHODAN TESTS: All " << nsucc << " tests passed." << std::endl;
}
