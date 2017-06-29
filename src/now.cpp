//=============================================================================
// now.cpp
//
//    Return the current time and date as a character string.
//
//    This is really silly that C++ does not offer such a function built-in
//    to the standard library.
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#include <ctime>
#include <string>

#include "now.h"

//-----------------------------------------------------------------------------
// Now
//
//    Return the current date and time as a string.
//-----------------------------------------------------------------------------
std::string Now()
{
   time_t rawtime;
   struct tm *timeinfo;

   time( &rawtime );
   timeinfo = localtime( &rawtime );

   std::string str( asctime(timeinfo) );
   str = str.substr(0, str.find('\n'));

   return str;
}
