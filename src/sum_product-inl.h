//=============================================================================
// sum_product-inl.h
//
//    A simple implementation of a core linear algebra computational component.
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#ifndef SUM_PRODUCT_H
#define SUM_PRODUCT_H

//-----------------------------------------------------------------------------
// This routine computes a dot product between two vectors.
//
// Arguments:
//
//    n     total number of elements in each vector.
//    x     pointer to the first element of the first vector.
//    y     pointer to the first element of the second vector.
//-----------------------------------------------------------------------------
inline double SumProduct( int n, const double* x, const double* y )
{
   double Sum = 0.0;

   for (int i=0; i<n; ++i)
      Sum += (*x++) * (*y++);

   return Sum;
}

//-----------------------------------------------------------------------------
// This routine computes a dot product between two vectors.  Both vectors
// allow for a non-unit stride.
//
// Arguments:
//
//    n     total number of elements in each vector.
//    x     pointer to the first element of the first vector.
//    dx    stride between subsequent elements in the first vector.
//    y     pointer to the first element of the second vector.
//    dy    stride between subsequent elements in the second vector.
//-----------------------------------------------------------------------------
inline double SumProduct( int n, const double* x, int dx, const double* y, int dy )
{
   double Sum = 0.0;

   for (int i=0; i<n; ++i)
   {
      Sum += (*x) * (*y);
      x += dx;
      y += dy;
   }

   return Sum;
}

//-----------------------------------------------------------------------------
// This routine computes a dot product between two vectors.  The y vector
// allows for a non-unit stride.
//
// Arguments:
//
//    n     total number of elements in each vector.
//    x     pointer to the first element of the first vector.
//    y     pointer to the first element of the second vector.
//    dy    stride between subsequent elements in the second vector.
//-----------------------------------------------------------------------------
inline double SumProduct( int n, const double* x, const double* y, int dy )
{
   double Sum = 0.0;

   for (int i=0; i<n; ++i)
   {
      Sum += (*x++) * (*y);
      y += dy;
   }

   return Sum;
}

//-----------------------------------------------------------------------------
// This routine computes a dot product between two vectors.  The x vector
// allows for a non-unit stride.
//
// Arguments:
//
//    n     total number of elements in each vector.
//    x     pointer to the first element of the first vector.
//    dx    stride between subsequent elements in the first vector.
//    y     pointer to the first element of the second vector.
//-----------------------------------------------------------------------------
inline double SumProduct( int n, const double* x, int dx, const double* y )
{
   double Sum = 0.0;

   for (int i=0; i<n; ++i)
   {
      Sum += (*x) * (*y++);
      x += dx;
   }

   return Sum;
}

//-----------------------------------------------------------------------------
// This routine computes a dot product between a vector and itself.
//
// Arguments:
//
//    n     total number of elements in each vector.
//    x     pointer to the first element of the vector.
//-----------------------------------------------------------------------------
inline double SumProduct( int n, const double* x )
{
   double Sum = 0.0;

   for (int i=0; i<n; ++i)
   {
      Sum += (*x) * (*x);
      ++x;
   }

   return Sum;
}

//-----------------------------------------------------------------------------
// This routine computes a dot product between a vector and itself.  The
// vector allows for a non-unit stride.
//
// Arguments:
//
//    n     total number of elements in each vector.
//    x     pointer to the first element of the vector.
//    dx    stride between subsequent elements in the vector.
//-----------------------------------------------------------------------------
inline double SumProduct( int n, const double* x, int dx )
{
   double Sum = 0.0;

   for (int i=0; i<n; ++i)
   {
      Sum += (*x) * (*x);
      x += dx;
   }

   return Sum;
}

//=============================================================================
#endif  // SUM_PRODUCT_H
