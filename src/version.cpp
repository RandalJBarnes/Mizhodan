//=============================================================================
// version.cpp
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

#include "version.h"

namespace
{
   const char* VERSION {"29 June 2017 [Beta]"};
}

//-----------------------------------------------------------------------------
void Banner( std::ostream& ost)
{
   ost <<
      "-------------------------------------------- \n"
      "Mizhodan (" << VERSION << ") \n"
      "\n"
      "Randal Barnes, University of Minnesota \n"
      "William Olsen,  Dakota County, Minnesota \n"
      "-------------------------------------------- \n"
   << std::endl;
}

//-----------------------------------------------------------------------------
void Help()
{
   Version();

   std::cout <<
      "   A basic two-dimensional Ordinary Kriging interpolator using an \n"
      "   exponential variogram model and all of the data. \n"
   << std::endl;

   Usage();

   std::cout <<
      "Arguments: \n"
      "   <nugget>        The <nugget>, or 'nugget effect', is the discontinuity \n"
      "                   in the variogram at a lag of 0. The <nugget> quantifies \n"
      "                   the variance of the sampling and measurement errors and \n"
      "                   the hyper-local spatial variation. The units of the \n"
      "                   <nugget> are the units of the observed values squared. \n"
      "                   The <nugget> must be a strictly positive value. \n"
      "\n"
      "   <sill>          The <sill> is the value (height) at which variogram \n"
      "                   levels out. With the exponential model used in Mizhodan, \n"
      "                   the variogram approaches the <sill> asymptotically. \n"
      "                   In the common geostatistical framework (i.e. a second \n"
      "                   order stationary model), the sill equals the variance \n"
      "                   of the underlying population. The units of the <sill> \n"
      "                   are the units of the observed values squared. The <sill> \n"
      "                   must be a strictly positive value. \n"
      "\n"
      "   <range>         The <range> is the separation distance at which we \n"
      "                   model two observations as essentially uncorrelated. \n"
      "                   With the exponential model used in Mizhodan, this is \n"
      "                   the separation distance at which the variogram reaches \n"
      "                   95% of the <sill>. The units of the <range> are the \n"
      "                   units of the observations locations. The <range> must \n"
      "                   be a strictly positive value. \n"
      "\n"
      "   <obs file>      The <obs file> is the name of the file (including any \n"
      "                   necessary path information and the .csv file extension) \n"
      "                   containing the observation data. \n"
      "\n"
      "   <targets file>  The <target file> is the name of the file (including any \n"
      "                   necessary path information and the .csv file extension) \n"
      "                   containing the targets data. \n"
      "\n"

      "   <results file>  The <results filename> is the name of the file (including \n"
      "                   any necessary path information and the file extension) \n"
      "                   where Mizhodan will write all of the program results. \n"
      "                   If the specified <output file> already exists, it will \n"
      "                   be overwritten. \n"
   << std::endl;

   std::cout <<
      "Example: \n"
      "   Mizhodan 3 25 3500 obs.csv target.csv results.csv \n"
   << std::endl;

   std::cout <<
      "Observation File: \n"
      "   All of the observation head data are supplied by this .csv file. \n"
      "\n"
      "   The observation file contains no header line. The observation file may \n"
      "   include blank lines, which are ignored. The observation may include \n"
      "   comment lines, which are identified by an octothorpe (#) in the first \n"
      "   column of the line. \n"
      "\n"
      "   The observation file contains one line for each head observation. Each \n"
      "   line in the observation file has four fields. \n"
      "\n"
      "   <ID>            The observation identification string. The ID string can \n"
      "                   contain numbers, letters, underscores, and internal spaces. \n"
      "                   The ID may not contain commas. \n"
      "\n"
      "   <x>             The x-coordinate [L] of observation location. \n"
      "\n"
      "   <y>             The y-coordinate [L] of observation location. \n"
      "\n"
      "   <z>             The observation value. at location (x,y). \n"
      "\n"
      "   Each of the four fields must separated by a single comma. Spaces and tabs \n"
      "   at the start and end of fields are trimmed. \n"
   << std::endl;

   std::cout <<
      "Targets File: \n"
      "   All of the estimation targets are identified by this .csv file. \n"
      "\n"
      "   The target file contains no header line. The target file may include blank \n"
      "   lines, which are ignored. The observation may include comment lines, which \n"
      "   are identified by an octothorpe (#) in the first column of the line. \n"
      "\n"
      "   The targets file contains one line for each target location. Each \n"
      "   line in the target file has three fields. \n"
      "\n"
      "   <ID>            The target identification string. The ID string can \n"
      "                   contain numbers, letters, underscores, and internal spaces. \n"
      "                   The ID may not contain commas. \n"
      "\n"
      "   <x>             The x-coordinate [L] of target location. \n"
      "\n"
      "   <y>             The y-coordinate [L] of target location. \n"
      "\n"
      "   Each of the three fields must separated by a single comma. Spaces and tabs \n"
      "   at the start and end of fields are trimmed. \n"
   << std::endl;

   std::cout <<
      "Results File: \n"
      "   All of the program results are written to the results .csv file. \n"
      "\n"
      "   The results file contains one header line with five comma separated \n"
      "   text fields containing the field titles. \n"
      "\n"
      "   The rest of the output file comprises one line for each target. \n"
      "   Each line has five fields. \n"
      "\n"
      "   <ID>            The target identification string. \n"
      "\n"
      "   <x>             The x-coordinate for the location of the target. \n"
      "\n"
      "   <y>             The y-coordinate for the location of the target. \n"
      "\n"
      "   <Zhat>          The interpolated value at the target location. \n"
      "\n"
      "   <Kstd>          The 'standard error' of the interpolated value at the \n"
      "                   target location. The <Kstd> is the square root \n"
      "                   Ordinary Kriging variance. \n"
   << std::endl;

   std::cout <<
      "Notes: \n"
      "   o  An exponential variogram model is used. \n"
      "         gamma(h) = <nugget> + (<sill>-<nugget>)*(1 - exp(-3h/<range>)) \n"
      "\n"
      "   o  The project name 'Mizhodan' is the Ojibwe word for the inanimate \n"
      "      transitive verb 'hit it (in shooting)'. See [http://ojibwe.lib.umn.edu]. \n"
   << std::endl;

   std::cout <<
      "Authors: \n"
      "   Randal Barnes, University of Minnesota \n"
      "   William Olsen, Dakota County, Minnesota \n"
   << std::endl;
}

//-----------------------------------------------------------------------------
void Usage()
{
   std::cout <<
      "Usage: \n"
      "   Mizhodan <nugget> <sill> <range> <obs file> <targets file> <results file> \n"
      "   Mizhodan --help \n"
      "   Mizhodan --version \n"
   << std::endl;
}

//-----------------------------------------------------------------------------
void Version()
{
   std::cout << std::endl;
   std::cout << "Mizhodan (" << VERSION << ")" << std::endl;
   std::cout << std::endl;
}
