//=============================================================================
// matrix.h
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    2 July 2017
//=============================================================================
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>

//=============================================================================
// Matrix
//=============================================================================
class Matrix
{
public:
   // Life cycle
   Matrix();                                          // null constructor
   Matrix( const Matrix& A );                         // copy constructor
   Matrix( const std::vector<double>v );              // constructor w/ std:vector

   Matrix( int nrows, int ncols );                    // dimensioned constructor
   Matrix( int nrows, int ncols, double a );          // constructor w/ scalar fill
   Matrix( int nrows, int ncols, const double* a );   // constructor w/ array fill
   Matrix( const std::string& str );

   ~Matrix();                                         // destructor
   void Resize( int nrows, int ncols );               // destructive resize.

   // Operators
   Matrix& operator=( const Matrix& A );              // assignment operator
   Matrix& operator=( double a );                     // scalar assignment

   double& operator()( int row, int col );            // mutable access
   double  operator()( int row, int col ) const;      // const access

   // Inquiry.
   int nRows() const;                                 // return the row size
   int nCols() const;                                 // return the column size

   // Access to the raw storage.
   const double* Base() const;                        // r/o access
   const double* Base( int row, int col ) const;      // r/o access

   double* Base();                                    // r/w access
   double* Base( int row, int col );                  // r/w access with offset

   // STL-like iterators.
   const double* begin() const;                       // r/o access
   const double* end() const;                         // r/o access

   double* begin();                                   // r/w access
   double* end();                                     // r/w access

private:
   int     m_nRows;                                   // allocated # of rows
   int     m_nCols;                                   // allocated # of columns
   double* m_Data;                                    // allocated memory
};


//=============================================================================
// IO Stream
//=============================================================================
std::ostream& operator << ( std::ostream& ostr, const Matrix& A );


//=============================================================================
// Matrix sums, measures and norms.
//=============================================================================
void ColumnSum( const Matrix& A, Matrix& x );
void RowSum( const Matrix& A, Matrix& x );

int Length( const Matrix& A );         // max( nRows, nCols )
double Trace( const Matrix& A );       // sum of the diagonal elements

double Sum( const Matrix& A );         // sum of all of the elements
double SumAbs( const Matrix& A );      // sum of the abs of all of the elements

double MaxAbs( const Matrix& A );      // maximum absolute value
double L1Norm( const Matrix& A );      // max column sum of abs
double LInfNorm( const Matrix& A );    // max row sum of abs
double FNorm( const Matrix& A );       // sqrt of sum of squares


//=============================================================================
// Unary Matrix operations.
//=============================================================================
void Transpose( const Matrix& A, Matrix& C );                        // C = A'
void Negative(  const Matrix& A, Matrix& C );                        // C = -A
void Identity( Matrix& A, int n );                                   // A = I(n)

//=============================================================================
// Slice Matrix operations.
//=============================================================================
void Slice( const Matrix& A, const std::vector<int>& row_flag, const std::vector<int>& col_flag, Matrix& C );
void SliceRows( const Matrix& A, const std::vector<int>& row_flag, Matrix& C );

//=============================================================================
// scalar/Matrix arithmetic routines.
//=============================================================================
void Add_aM( double a, const Matrix& A, Matrix& C );                 // C = a+A
void Subtract_aM( double a, const Matrix& A, Matrix& C );            // C = a-A
void Multiply_aM( double a, const Matrix& A, Matrix& C );            // C = a*A

//=============================================================================
// Matrix/Matrix addition and subtraction.
//=============================================================================
void Add_MM     ( const Matrix& A, const Matrix& B, Matrix& C );     // C = A + B
void Subtract_MM( const Matrix& A, const Matrix& B, Matrix& C );     // C = A - B

//=============================================================================
// Matrix/Matrix multiplication
//=============================================================================
void Multiply_MM  ( const Matrix& A, const Matrix& B, Matrix& C );   // C = AB
void Multiply_MtM ( const Matrix& A, const Matrix& B, Matrix& C );   // C = A'B
void Multiply_MMt ( const Matrix& A, const Matrix& B, Matrix& C );   // C = AB'
void Multiply_MtMt( const Matrix& A, const Matrix& B, Matrix& C );   // C = A'B'

//=============================================================================
// Dot Products
//=============================================================================
double DotProduct( const Matrix& A, const Matrix& B );

//=============================================================================
// Quadratic forms
//=============================================================================
double QuadraticForm_MtMM( const Matrix& a, const Matrix& B, const Matrix& c );  // a'Bc
double QuadraticForm_MMM ( const Matrix& a, const Matrix& B, const Matrix& c );  // aBc

//=============================================================================
// Matrix comparison
//=============================================================================
bool isSquare(const Matrix& A);
bool isCongruent( const Matrix& A, const Matrix& B );
bool isClose( const Matrix& A, const Matrix& B, double tol );

//=============================================================================
// isVector
//=============================================================================
bool isRow( const Matrix& A );
bool isCol( const Matrix& A );
bool isVector( const Matrix& A );

//=============================================================================
#endif  // MATRIX_H
