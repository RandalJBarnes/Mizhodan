//=============================================================================
// numerical_constants.h
//
//    A small collection of numerical constants that have proven useful.
//
// notes:
// o  All constants are given with 34 significant digits in anticipation of
//    future higher precision floating point representations.  Only the first
//    15, or so, digits are significant with <double>.
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#ifndef NUMERICAL_CONSTANTS_H
#define NUMERICAL_CONSTANTS_H

#include <limits>

namespace{
   const double ONE_PI               = 3.141592653589793238462643383279503;
   const double TWO_PI               = 6.283185307179586476925286766559006;
   const double FOUR_PI              = 12.56637061435917295385057353311801;

   const double HALF_PI              = 1.570796326794896619231321691639751;
   const double QUARTER_PI           = 0.7853981633974483096156608458198757;

   const double ONE_OVER_PI          = 0.3183098861837906715377675267450287;
   const double ONE_OVER_TWO_PI      = 0.1591549430918953357688837633725144;
   const double ONE_OVER_FOUR_PI     = 0.07957747154594766788444188168625718;

   const double SQRT_PI              = 1.772453850905516027298167483341145;
   const double SQRT_TWO_PI          = 2.506628274631000502415765284811045;

   const double ONE_OVER_SQRT_PI     = 0.5641895835477562869480794515607726;
   const double TWO_OVER_SQRT_PI     = 1.128379167095512573896158903121545;
   const double ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399460599343819;

   const double RAD_TO_DEG           = 57.29577951308232087679815481410517;
   const double DEG_TO_RAD           = 0.01745329251994329576923690768488613;

   const double LN_TEN               = 2.302585092994045684017991454684364;
   const double EULER_E              = 2.718281828459045235360287471352662;
   const double EULER_GAMMA          = 0.5772156649015328606065120900824024;

   const double GOLDEN_RATIO         = 1.618033988749894848204586834365638;

   const double SQRT_TWO             = 1.414213562373095048801688724209698;
   const double LN_TWO               = 0.6931471805599453094172321214581766;

   const double EPS = std::numeric_limits<double>::epsilon();
   const double INF = std::numeric_limits<double>::max();
}

//=============================================================================
#endif  // NUMERICAL_CONSTANTS_H
