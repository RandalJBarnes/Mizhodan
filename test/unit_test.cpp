//=============================================================================
// unit_test.cpp
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

//-----------------------------------------------------------------------------
bool isClose( double x, double y, double tol )
{
   return( abs(x-y) < tol );
}

//-----------------------------------------------------------------------------
bool Check( bool test, int line, const char* file )
{
   if (!test)
      std::cerr << "FAILED TEST: " << "\t LINE: " << line << "\t FILE: " << file << std::endl;

   return test;
}
