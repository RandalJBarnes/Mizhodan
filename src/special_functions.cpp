//=============================================================================
// special_functions.cpp
//
//    A small collection of special functions.
//
// references:
// o  S. Zhang and J. Jin, 1996, Computation of Special Functions, John Wiley
//    and Sons, ISBN 0-471-11963-6.
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

#include "numerical_constants.h"
#include "special_functions.h"

//-----------------------------------------------------------------------------
// Beta
//
//    \Beta(a,b) = \int_0^1 x^{a-1} (1-x)^{b-1} dx
//               = \frac{\Gamma(a) \Gamma(a)}{\Gamma(a+b)}
//
// notes:
// o  a > 0 and b > 0.
//-----------------------------------------------------------------------------
double Beta( double a, double b )
{
   assert( a > 0 && b > 0 );

   return Gamma(a)*Gamma(b)/Gamma(a+b);
}

//-----------------------------------------------------------------------------
// IncompleteBeta
//
//    I_x(a,b) = 1/\Beta(a,b) \int_0^x t^{a-1} (1-t)^{b-1} dt
//
// notes:
// o  This evaluation is based upon the continued fraction representation.
//    The formulae are taken from Section 3.5 of Zhang and Jin 1996.
//
// o  See also Section 11.3 (pages 288-298) of Temme 1996.
//
// references:
// o  Nico M. Temme, 1996, Special Functions an Introduction to the Classical
//    Functions of Mathematical Physics, John Wiley and Sons,
//    ISBN 0-471-11313-1.
//
// o  S. Zhang and J. Jin, 1996, Computation of Special Functions, John Wiley
//    and Sons, ISBN 0-471-11963-6.
//-----------------------------------------------------------------------------
double IncompleteBeta( double x, double a, double b )
{
   assert( a > 0 && b > 0 );
   assert( x >= 0 && x <= 1);

   // Check the end points.
   if (abs(x-1.0) < EPS) return 1.0;
   if (abs(x) < EPS) return 0.0;

   // Compute using the continued fraction.
   const int M = 20;    // 2*M+1 terms in the continued fraction.
   double T = 0;

   if (x < a/(a+b)) {                                       // Zhang (3.5.7)
      for (int n = 2*M+1; n >= 1; --n) {
         double d;

         if (n%2 == 1) {
            double m = (n-1)/2;
            d = -(a+m)*(a+b+m)/(a+2*m)/(a+2*m+1) * x;
         }
         else {
            double m = n/2;
            d = m*(b-m)/(a+2*m-1)/(a+2*m) * x;
         }
         T = d/(1+T);
      }
      T = 1/(1+T);

      return pow(x,a) * pow(1-x,b) / (a*Beta(a,b)) * T;
   }
   else {                                                   // Zhang (3.5.9)
      for (int n = 2*M+1; n >= 1; --n) {
         double d;

         if (n%2 == 1) {
            double m = (n-1)/2;
            d = -(b+m)*(a+b+m)/(b+2*m)/(b+2*m+1) * (1-x);
         }
         else {
            double m = n/2;
            d = m*(a-m)/(b+2*m-1)/(b+2*m) * (1-x);
         }
         T = d/(1+T);
      }
      T = 1/(1+T);

      return 1 - pow(x,a) * pow(1-x,b) / (b*Beta(a,b)) * T;
   }
}

//-----------------------------------------------------------------------------
// IncompleteBetaInv
//
//    Returns x such that I_x(a,b) = p for an argument 0 <= p <= 1;
//
// notes:
// o  Start with a few iterations of a simple bisection scheme, then switch
//    to Halley's method.
//-----------------------------------------------------------------------------
double IncompleteBetaInv( double p, double a, double b )
{
   assert( a > 0 && b > 0 );
   assert( p >= 0 && p <= 1 );

   // Check the end points.
   if (abs(p-1.0) <= EPS) return 1.0;
   if (abs(p) <= EPS) return 0.0;

   // Start with a bisection scheme.
   double xL = 0.0;
   double xR = 1.0;
   double x  = 0.5;
   for (int j = 0; j < 12; ++j) {
      x = (xL+xR)/2;
      if ( IncompleteBeta(x,a,b) > p )
         xR = x;
      else
         xL = x;
   }

   // Halley's iterations.
   double ba = Beta(a,b);

   for (int j = 0; j < 12; ++j) {
      double f   = IncompleteBeta(x,a,b) - p;         // error
      double df  = pow(x,a-1) * pow(1-x,b-1) / ba;    // dI/dx
      double ddf = df * ( (a-1)/x - (b-1)/(1-x) );    // d^2P/dx^2

      double deltax = f/(df - f*ddf/(2*df));          // Halley's iteration.
      double xnew = x - deltax;

      if (xnew >= 1.0)
         x = (1+x)/2;
      else if (xnew <= 0.0)
         x = x/2;
      else
         x = xnew;

      if (fabs(deltax) < EPS*x ) break;
   }

   return x;
}

//-----------------------------------------------------------------------------
// Gamma
//
//    \Gamma(x) = \int_0^{\infty} t^{x-1} e^{-t} dt
//
// notes:
// o  This function is based on the FORTRAN subroutine presented in Section
//    3.1.5 (page 49-50) of Zhang and Jin 1996.
//
// references:
// o  S. Zhang and J. Jin, 1996, Computation of Special Functions, John Wiley
//    and Sons, ISBN 0-471-11963-6.
//-----------------------------------------------------------------------------
double Gamma( double x )
{
   double ga;

   const double g[] = {
      1.0,
      0.5772156649015329,
     -0.6558780715202538,
     -0.420026350340952e-1,
      0.1665386113822915,
     -0.421977345555443e-1,
     -0.9621971527877e-2,
      0.7218943246663e-2,
     -0.11651675918591e-2,
     -0.2152416741149e-3,
      0.1280502823882e-3,
     -0.201348547807e-4,
     -0.12504934821e-5,
      0.1133027232e-5,
     -0.2056338417e-6,
      0.6116095e-8,
      0.50020075e-8,
     -0.11812746e-8,
      0.1043427e-9,
      0.77823e-11,
     -0.36968e-11,
      0.51e-12,
     -0.206e-13,
     -0.54e-14,
      0.14e-14 };

   // Protect against overflow.
   if (x > 171.0) return INF;

   // Handle the special case of an integer argument.
   if (abs(x-floor(x)) <= EPS )
   {
      // When x == n > 0, use (3.1.5).
      if (x > 0.0)
      {
         ga = 1.0;
         for (int k = 2; k < x; ++k) {
            ga *= k;
         }
      }
      else
      {
         ga = INF;
      }
   }
   else
   {
      double r,z;

      // When |X| > 1, use (3.1.9).
      if (fabs(x) > 1.0) {
         z = fabs(x);
         int m = static_cast<int>(z);
         r = 1.0;

         for (int k = 1; k <= m; ++k) {
            r *= (z-k);
         }
         z -= m;
      }
      else
      {
         z = x;
      }

      // Calculate 1/Gamma(x) using (3.1.15)
      double gr = g[24];

      for (int k = 23; k >= 0; --k) {
         gr = gr*z + g[k];
      }
      ga = 1.0/(gr*z);

      // When |x| > 1, use (3.1.9)
      if (fabs(x) > 1.0) {
         ga *= r;
         if (x < 0.0) {
             ga = -ONE_PI/(x*ga*sin(ONE_PI*x));
         }
      }
   }

   return ga;
}

//-----------------------------------------------------------------------------
// IncompleteGamma
//
//    P(a,x) = 1/\Gamma(a) \int_0^x t^{a-1} e^{-t} dt
//
// notes:
// o  This function is motivated by the FORTRAN subroutine presented in
//    Section 3.4 (page 63) of Zhang and Jin 1996.  See also Section 11.2
//    of Temme 1996.
//
// references:
// o  Nico M. Temme, 1996, Special Functions an Introduction to the Classical
//    Functions of Mathematical Physics, John Wiley and Sons,
//    ISBN 0-471-11313-1.
//
// o  S. Zhang and J. Jin, 1996, Computation of Special Functions, John Wiley
//    and Sons, ISBN 0-471-11963-6.
//-----------------------------------------------------------------------------
double IncompleteGamma( double x, double a )
{
   assert( a > 0 && a < 170 );
   assert( x >= 0 );

   if (abs(x) <= EPS) {                           // special case.
      return 0.0;
   }
   else if (x <= 1+a)                        // Zhang (3.4.4)
   {
      double s = 1/a;
      double r = s;

      for (int k = 1; k <= 60; ++k) {
         r *= x/(a+k);
         s += r;
         if (fabs(r/s) < 1e-15) break;
      }
      return s*exp(a*log(x)-x) / Gamma(a);
   }
   else                                      // Zhang (3.4.11)
   {
      double T = 0;
      for (int k = 60; k >= 1; --k)
         T = (k-a)/(1+ k/(x+T));
      return 1 - exp(a*log(x)-x)/(x+T) / Gamma(a);
   }
}

//-----------------------------------------------------------------------------
// IncompleteGammaInv
//
//    Returns x such that P(a,x) = p for an argument 0 <= p < 1;
//
// notes:
// o  This routine is based upon an initial guess followed by multiple
//    iterations of Halley's method.
//
// o  The initial guess is given by Press et al. (2007), as
//
//       x = a t^3
//
//    where
//
//       t = 1 - d - GaussianInv(p) * \sqrt(d)
//       d = 1/(9a)
//
// references:
// o  W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, 2007,
//    Numerical Recipes, 3rd edition, Cambridge Unversity Press, 1256 pp.,
//    ISBN 978-0521880688.
//-----------------------------------------------------------------------------
double IncompleteGammaInv( double p, double a )
{
   assert( p >= 0 && p <= 1 );

   double x;

   // End point cases.
   if (abs(p) <= EPS) return 0.0;
   if (abs(p-1.0) <= EPS) return INF;

   // Initial guess from Press et al. (2007).
   if (a > 1) {
      double d = 1/(9*a);
      double t = 1 - d - GaussianCDFInv(p) * sqrt(d);
      x = a*t*t*t;
   }
   else {
      double t = 1.0 - (0.253 + 0.12*a)*a;
      if (p < t)
         x = pow(p/t,1/a);
      else
         x = 1 - log(1 - (p-t)/(1-t));
   }

   // Halley's iterations.
   double ga = Gamma(a);

   for (int j = 0; j < 12; ++j) {
      if (x <= 0.0) return 0.0;                    // x is too small to compute accurately.

      double f   = IncompleteGamma(x,a) - p;       // error
      double df  = pow(x,a-1) * exp(-x) / ga;      // dP/dx
      double ddf = df * ( (a-1)/x - 1 );           // d^2P/dx^2

      double deltax = f/(df - f*ddf/(2*df));       // Halley's iteration.
      x -= deltax;

      if (fabs(deltax) < EPS*x ) break;
   }

   return x;
}

//-----------------------------------------------------------------------------
// GaussianCDF
//
//    Returns the value of the Standard Normal cumulative distribution
//    function at the argument.
//
// notes:
// o  This routine is relatively slow, but the associated error is less than
//    1e-15 for all x.
//
// references:
// o  Marsaglia, George, 2004, Evaluating the Normal Distribution, Journal of
//    Statistical Software, v. 11, n. 4, June 2004.  Available on-line at
//    http://www.jstatsoft.org/v11/a05/paper.
//-----------------------------------------------------------------------------
double GaussianCDF(double x)
{
   if (x < -8.0)
      return 0.0;
   else if (x > 8.0)
      return 1.0;
   else
   {
      long double s=x, t=0, b=x, q=x*x, i=1;
      while (abs(s-t) > EPS)
         s = (t=s) + (b *= q/(i+=2));
      return static_cast<double>( 0.5 + s*exp(-0.5*q - 0.91893853320467274178L) );
   }
}

//-----------------------------------------------------------------------------
// GaussianCDFInv
//
// notes:
// o  Based on (26.2.23) on page 933 of Abramowitz and Stegun 1972. This
//    yields an error < 4.5e4.
//
// o  Applying Halley's "one-point third-order" iterative update yields
//    almost full machine precision.
//
// o  The idea of applying Halley's formula is based on a C program written by
//    Jeremy Lea (Research Engineer, Infrastructure Engineering, Transportek,
//    CSIR, South Africa), which I found online at:
//    http://home.online.no/~pjacklam/notes/invnorm/#C.
//
// references:
// o  Milton Abromowitz and Irene Stegun, 1972, Handobook of Mathematical
//    Functions, Dover Publications, ISBN 978-0-486-61272-0.
//
// o  For more on Halley's formula (the same guy with the comet), see
//    http://en.wikipedia.org/wiki/Simple_rational_approximation.
//-----------------------------------------------------------------------------
double GaussianCDFInv(double p)
{
   assert(p>0 && p<1);

   const double c[] = { 2.515517, 0.802853, 0.010328 };
   const double d[] = { 1.432788, 0.189269, 0.001308 };

   double q = (p<0.5 ? p : 1-p);

   double t = sqrt(-2*log(q));
   double num = c[0] + (c[1] + c[2]*t)*t;
   double den = 1 + (d[0] + (d[1] + d[2]*t)*t)*t;
   double u = -t + num/den;

   t = GaussianCDF(u) - q;             // error
   t = t * SQRT_TWO_PI*exp(u*u/2);     // f(u)/df(u)
   u = u - t/(1 + u*t/2);              // Halley's update formula

   return (p<0.5 ? u : -u);
}
