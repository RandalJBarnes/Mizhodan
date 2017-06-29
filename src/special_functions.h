//=============================================================================
// special_functions.h
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#ifndef SPEICAL_FUNCTIONS_H
#define SPECIAL_FUNCTIONS_H

#include <iostream>

double Beta( double a, double b );
double IncompleteBeta( double x, double a, double b );
double IncompleteBetaInv( double p, double a, double b );

double Gamma( double x );
double IncompleteGamma( double x, double a );
double IncompleteGammaInv( double p, double a );

double GaussianCDF( double x );
double GaussianCDFInv( double p );

//=============================================================================
#endif  // SPECIAL_FUNCTIONS_H
