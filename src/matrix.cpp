//=============================================================================
// matrix.cpp
//
//    A minimal matrix class with some basic operations and arithmetic.  This
//    matrix class is explicitly based on <double>. This matrix class is NOT
//    a template.
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    2 July 2017
//=============================================================================
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

#include "matrix.h"
#include "sum_product-inl.h"

//=============================================================================
// Matrix
//=============================================================================

//-----------------------------------------------------------------------------
// Null constructor.
//-----------------------------------------------------------------------------
Matrix::Matrix()
:  m_nRows( 0 ),
   m_nCols( 0 ),
   m_Data( nullptr ) 
{
}

//-----------------------------------------------------------------------------
// Copy constructor.
//-----------------------------------------------------------------------------
Matrix::Matrix( const Matrix& A )
:  m_nRows( 0 ),
   m_nCols( 0 ),
   m_Data( nullptr )
{
   if ( A.nRows() > 0 && A.nCols() > 0 ) {
      m_nRows = A.nRows();
      m_nCols = A.nCols();
      m_Data  = new double[ m_nRows*m_nCols ];
      memcpy( m_Data, A.Base(), sizeof(double)*m_nRows*m_nCols );
   }
}

//-----------------------------------------------------------------------------
// constructor from an std:vector
//-----------------------------------------------------------------------------
Matrix::Matrix( const std::vector<double>v )
:  m_nRows( 0 ),
   m_nCols( 0 ),
   m_Data( nullptr )
{
   if ( v.size() > 0 ) {
      m_nRows = v.size();
      m_nCols = 1;
      m_Data  = new double[ m_nRows ];

      for (int k = 0; k < m_nRows; ++k)
         m_Data[k] = v[k];
   }
}

//-----------------------------------------------------------------------------
// Dimensioned constructor, with zero fill.
//-----------------------------------------------------------------------------
Matrix::Matrix( int nrows, int ncols )
:  m_nRows( 0 ),
   m_nCols( 0 ),
   m_Data( nullptr )
{
   assert( nrows >= 0 && ncols >= 0 );

   m_nRows = nrows;
   m_nCols = ncols;
   m_Data  = new double[ m_nRows*m_nCols ];
   memset( m_Data, 0, sizeof(double)*m_nRows*m_nCols );
}

//-----------------------------------------------------------------------------
// Constructor with scalar fill.
//-----------------------------------------------------------------------------
Matrix::Matrix( int nrows, int ncols, double a )
:  m_nRows( 0 ),
   m_nCols( 0 ),
   m_Data( nullptr )
{
   assert( nrows >= 0 && ncols >= 0 );

   m_nRows = nrows;
   m_nCols = ncols;
   m_Data  = new double[ m_nRows*m_nCols ];

   for (int i = 0; i < nrows; ++i)
      for (int j = 0; j < ncols; ++j)
         m_Data[i*ncols + j] = a;
}

//-----------------------------------------------------------------------------
// Constructor with array fill.
//-----------------------------------------------------------------------------
Matrix::Matrix( int nrows, int ncols, const double* data )
:  m_nRows( 0 ),
   m_nCols( 0 ),
   m_Data( nullptr )
{
   assert( nrows >= 0 && ncols >= 0 );

   m_nRows = nrows;
   m_nCols = ncols;
   m_Data  = new double[ m_nRows*m_nCols ];
   memcpy( m_Data, data, sizeof(double)*m_nRows*m_nCols );
}

//-----------------------------------------------------------------------------
// Constructor with string-based initialization.
//
//    Columns are separated by comma, rows are separated by semicolons.
//    Missing values are filled with zeros.  For example,
//
//       Matrix A("1,2,3;4,5,6")
//       Matrix B("1,,;,,6");
//       Matrix C("1,2,3;4,5,6;");
//
//    would create
//
//       A = [ 1  2  3 ]      B = [ 1  0  0 ]      C = [ 1  2  3 ]
//           [ 4  5  6 ]          [ 0  0  6 ]          [ 4  5  6 ]
//                                                     [ 0  0  0 ]
//
//    Any token that cannot be interpreted as a valid double is set to zero.
//-----------------------------------------------------------------------------
Matrix::Matrix( const std::string& str )
:  m_nRows( 0 ),
   m_nCols( 0 ),
   m_Data( nullptr )
{
   assert( str.find_first_not_of("-0123456789eE.,; \t") == std::string::npos );

   // Parse the string, storing the values in a vector of vectors.
   std::vector< std::vector< double > > rows;

   std::string::size_type beg_line = str.find_first_not_of(" \t", 0);
   while (beg_line != std::string::npos)
   {
      std::string::size_type end_line = str.find_first_of(";", beg_line);
      std::string line = str.substr(beg_line, end_line-beg_line);

      std::vector< double > row;

      std::string::size_type beg_token = line.find_first_not_of(" \t", 0);
      while (beg_token != std::string::npos)
      {
         std::string::size_type end_token = line.find_first_of(",", beg_token);
         std::string token = line.substr(beg_token, end_token-beg_token);

         std::istringstream iss( token );
         double value = 0.0;
         if ((iss >> value).fail()) value = 0.0;
         row.push_back( value );

         if (end_token == std::string::npos) break;
         beg_token = str.find_first_not_of(" \t", end_token+1);
      }
      rows.push_back( row );

      if (end_line == std::string::npos) break;
      beg_line = str.find_first_not_of(" \t", end_line+1);
   }

   // Construct the Matrix from the parsed values.
   m_nRows = rows.size();

   m_nCols = 0;
   for (std::vector< std::vector< double > >::const_iterator i = rows.begin(); i != rows.end(); ++i) {
      if ( static_cast<int>(i->size()) > m_nCols) m_nCols = i->size();
   }

   m_Data  = new double[ m_nRows*m_nCols ];
   memset( m_Data, 0, sizeof(double)*m_nRows*m_nCols );

   for (std::vector<std::vector<double>>::const_iterator i = rows.begin(); i != rows.end(); ++i)
      for (std::vector<double>::const_iterator j = i->begin(); j != i->end(); ++j) {
         int k = (i - rows.begin())*m_nCols + (j - i->begin());
         m_Data[k] = *j;
      }
}

//-----------------------------------------------------------------------------
// Destructor.
//-----------------------------------------------------------------------------
Matrix::~Matrix()
{
   delete [] m_Data;

   m_nRows = 0;
   m_nCols = 0;
   m_Data  = nullptr;
}

//-----------------------------------------------------------------------------
// Destructive resize.
//
//    The resized Matrix is filled with zeros.
//-----------------------------------------------------------------------------
void Matrix::Resize( int nrows, int ncols )
{
   // Check the arguments.
   assert( nrows >= 0 && ncols >= 0 );

   // Reallocate memory if necessary.
   if (m_nRows != nrows || m_nCols != ncols) {
      delete [] m_Data;

      if ( nrows > 0 && ncols > 0 ) {
         m_nRows = nrows;
         m_nCols = ncols;
         m_Data  = new double[ m_nRows*m_nCols ];
      }
      else {
         m_nRows = 0;
         m_nCols = 0;
         m_Data  = nullptr;
      }
   }

   memset( m_Data, 0, sizeof(double)*m_nRows*m_nCols );
}

//-----------------------------------------------------------------------------
// Assignment operator.
//-----------------------------------------------------------------------------
Matrix& Matrix::operator=( const Matrix& A )
{
   // Check for self-assignment.
   if ( this == &A ) return *this;

   // Commensurate memory allocation.
   Resize( A.nRows(), A.nCols() );

   // Copy the data.
   if ( m_nRows*m_nCols > 0 )
      memcpy( m_Data, A.Base(), sizeof(double)*m_nRows*m_nCols );

   return *this;
}

//-----------------------------------------------------------------------------
// Scalar assignment operator.
//-----------------------------------------------------------------------------
Matrix& Matrix::operator=( double a )
{
   for (int k = 0; k < m_nRows*m_nCols; ++k)
      m_Data[k] = a;

   return *this;
}

//-----------------------------------------------------------------------------
// Non-constant element access operator (put).
//-----------------------------------------------------------------------------
double& Matrix::operator()( int row, int col )
{
   assert( row >= 0 && row < m_nRows );
   assert( col >= 0 && col < m_nCols );

   return m_Data[ row*m_nCols + col ];
}

//-----------------------------------------------------------------------------
// Constant element access operator (get).
//-----------------------------------------------------------------------------
double Matrix::operator()( int row, int col ) const
{
   assert( row >= 0 && row < m_nRows );
   assert( col >= 0 && col < m_nCols );

   return m_Data[ row*m_nCols + col ];
}

//-----------------------------------------------------------------------------
// Number of rows.
//-----------------------------------------------------------------------------
int Matrix::nRows() const
{
   return m_nRows;
}

//-----------------------------------------------------------------------------
// Number of columns.
//-----------------------------------------------------------------------------
int Matrix::nCols() const
{
   return m_nCols;
}

//-----------------------------------------------------------------------------
// Read only access to raw storage.
//-----------------------------------------------------------------------------
const double* Matrix::Base() const
{
   return m_Data;
}

//-----------------------------------------------------------------------------
// Read only access to raw storage with an offset.
//-----------------------------------------------------------------------------
const double* Matrix::Base( int row, int col ) const
{
   assert( row >= 0 && row < m_nRows );
   assert( col >= 0 && col < m_nCols );

   return m_Data + row*m_nCols + col;
}

//-----------------------------------------------------------------------------
// Read/Write access to raw storage.
//-----------------------------------------------------------------------------
double* Matrix::Base()
{
   return m_Data;
}

//-----------------------------------------------------------------------------
// Read/Write access to raw storage with an offset.
//-----------------------------------------------------------------------------
double* Matrix::Base( int row, int col )
{
   assert( row >= 0 && row < m_nRows );
   assert( col >= 0 && col < m_nCols );

   return m_Data + row*m_nCols + col;
}

//-----------------------------------------------------------------------------
// Read only STL-conforming begin() iterators.
//-----------------------------------------------------------------------------
const double* Matrix::begin() const
{
   return m_Data;
}

const double* Matrix::end() const
{
   return m_Data + m_nRows*m_nCols;
}

//-----------------------------------------------------------------------------
// Read/write STL-conforming end() iterator.
//-----------------------------------------------------------------------------
double* Matrix::begin()
{
   return m_Data;
}

double* Matrix::end()
{
   return m_Data + m_nRows*m_nCols;
}


//=============================================================================
// I/O routines.
//=============================================================================

//-----------------------------------------------------------------------------
// Output operator.
//-----------------------------------------------------------------------------
std::ostream& operator <<( std::ostream& ostr, const Matrix& A )
{
   // Output the Matrix.
   for (int i = 0; i < A.nRows(); ++i) {
      for (int j = 0; j < A.nCols(); ++j) {
         ostr << std::fixed << std::setw(12) << std::setprecision(3) << A(i,j);
      }
      ostr << std::endl;
   }

   return ostr;
}


//=============================================================================
// Matrix measures and norms.
//=============================================================================

//-----------------------------------------------------------------------------
// Row Matrix of column sums.
//-----------------------------------------------------------------------------
void ColumnSum( const Matrix& A, Matrix& x )
{
   x.Resize( 1, A.nCols() );
   x = 0.0;

   for (int i = 0; i < A.nRows(); ++i)
      for (int j = 0; j < A.nCols(); ++j)
         x(0,j) += A(i,j);
}

//-----------------------------------------------------------------------------
// Column Matrix of row sums.
//-----------------------------------------------------------------------------
void RowSum( const Matrix& A, Matrix& x )
{
   x.Resize( A.nRows(), 1 );
   x = 0.0;

   for (int j = 0; j < A.nCols(); ++j)
      for (int i = 0; i < A.nRows(); ++i)
         x(i,0) += A(i,j);
}

//-----------------------------------------------------------------------------
// Matrix Length = max( nRows, nCols )
//-----------------------------------------------------------------------------
int Length( const Matrix& A )
{
   return std::max( A.nRows(), A.nCols() );
}

//-----------------------------------------------------------------------------
// Matrix Trace.
//
//    The trace of a square Matrix is the sum of the diagonal elements.
//    The trace is not defined for a non-square Matrix.
//
// See Golub and Van Loan, 1996, p. 310.
//-----------------------------------------------------------------------------
double Trace( const Matrix& A )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );
   assert( A.nRows() == A.nCols() );

   // Compute the Matrix trace: i.e. sum of the diagonal elements.
   double Sum  = 0.0;
   for (int i = 0; i < A.nRows(); ++i)
      Sum += A(i,i);

   return Sum;
}

//-----------------------------------------------------------------------------
// Matrix Sum.
//
//    Sum all of the elements in the matrix.
//-----------------------------------------------------------------------------
double Sum( const Matrix& A )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );

   return std::accumulate( A.begin(), A.end(), 0.0 );
}

//-----------------------------------------------------------------------------
// Matrix SumAbs.
//
//    Sum the absolute values all of the elements in the matrix.
//-----------------------------------------------------------------------------
double SumAbs( const Matrix& A )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );

   return std::accumulate( A.begin(), A.end(), 0.0, [](double a, double b){return a + abs(b);} );
}

//-----------------------------------------------------------------------------
// Matrix MaxAbs.
//
//    The MaxAbs of a Matrix is the maximum of the absolute values of the
//    elements.
//-----------------------------------------------------------------------------
double MaxAbs( const Matrix& A )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );

   return std::accumulate( A.begin(), A.end(), 0.0, [](double a, double b){return std::max(a, abs(b));} );
}

//-----------------------------------------------------------------------------
// Matrix L1Norm.
//
//    The L1 norm of a Matrix is the maximum column sum of the absolute values
//    of the elements.
//
// See Golub and Van Loan, 1996, p. 56, (2.3.9).
//-----------------------------------------------------------------------------
double L1Norm( const Matrix& A )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );

   // Compute the L1 norm.
   double MaxColSum = 0;

   for (int j = 0; j < A.nCols(); ++j) {
      double Sum = 0;
      for (int i = 0; i < A.nRows(); ++i)
         Sum += fabs( A(i,j) );
      if (Sum > MaxColSum) MaxColSum = Sum;
   }

   return MaxColSum;
}

//-----------------------------------------------------------------------------
// Matrix LInfNorm.
//
//    The L-infinity norm of a Matrix is the maximum row sum of the absolute
//    values of the elements.
//
// See Golub and Van Loan, 1996, p. 56, (2.3.10).
//-----------------------------------------------------------------------------
double LInfNorm( const Matrix& A )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );

   // Compute the LInf norm.
   double MaxRowSum = 0;

   for (int i = 0; i < A.nRows(); ++i) {
      double Sum = 0;
      for (int j = 0; j < A.nCols(); ++j)
         Sum += fabs( A(i,j) );
      if (Sum > MaxRowSum) MaxRowSum = Sum;
   }

   return MaxRowSum;
}

//-----------------------------------------------------------------------------
// Matrix FNorm.
//
//    The Frobenius norm of a Matrix is the square root of the sum of the
//    squares of the elements.
//
// See Golub and Van Loan, 1996, p. 55, (2.3.1).
//-----------------------------------------------------------------------------
double FNorm( const Matrix& A )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );

   return sqrt( std::accumulate(A.begin(), A.end(), 0.0, [](double a, double b){return a + b*b;}) );
}


//=============================================================================
// Unary Matrix operations.
//=============================================================================

//-----------------------------------------------------------------------------
// Matrix transpose : C = A'
//-----------------------------------------------------------------------------
void Transpose( const Matrix& A, Matrix& C )
{
   // Check the arguments.
   assert( A.nCols() > 0 && A.nRows() > 0 );

   // Commensurate memory allocation.
   Matrix At( A.nCols(), A.nRows() );

   // Set the transpose.
   for (int i = 0; i < A.nRows(); ++i)
      for (int j = 0; j < A.nCols(); ++j)
         At(j,i) = A(i,j);

   C = At;
}

//-----------------------------------------------------------------------------
// Matrix negation : C = -A
//-----------------------------------------------------------------------------
void Negative( const Matrix& A, Matrix& C )
{
   // Check the arguments.
   assert( A.nCols() > 0 && A.nRows() > 0 );

   C.Resize( A.nRows(), A.nCols() );
   std::transform( A.begin(), A.end(), C.begin(), [](double a){return -(a);});
}

//-----------------------------------------------------------------------------
// Reset a Matrix to an n by n identity Matrix.
//-----------------------------------------------------------------------------
void Identity( Matrix& A, int n )
{
   // Make A a square n x n Matrix
   A.Resize( n, n );

   // Set up the identity Matrix.
   for (int i = 0; i < n; ++i)
      A(i,i) = 1;
}

//-----------------------------------------------------------------------------
// Slice Matrix operations.
//-----------------------------------------------------------------------------
void Slice( const Matrix& A, const std::vector<int>& row_flag, const std::vector<int>& col_flag, Matrix& C )
{
   assert( int(row_flag.size()) == A.nRows() );
   assert( int(col_flag.size()) == A.nCols() );

   int nRows = std::count_if( row_flag.begin(), row_flag.end(), [](int i){return i != 0;} );
   int nCols = std::count_if( col_flag.begin(), col_flag.end(), [](int i){return i != 0;} );
   C.Resize( nRows, nCols );

   if( nRows*nCols > 0 ) {
      int row = 0;
      for (int i = 0; i < A.nRows(); ++i) {
         if (row_flag[i] != 0) {
            int col = 0;
            for (int j = 0; j < A.nCols(); ++j) {
               if (col_flag[j] != 0) {
                  C(row, col) = A(i,j);
                  ++col;
               }
            }
            ++row;
         }
      }
   }
}

//-----------------------------------------------------------------------------
// SliceRows Matrix operations.
//-----------------------------------------------------------------------------
void SliceRows( const Matrix& A, const std::vector<int>& row_flag, Matrix& C ) {
   assert( int(row_flag.size()) == A.nRows() );

   int nRows = std::count_if( row_flag.begin(), row_flag.end(), [](int i) {
      return i != 0;
   } );
   int nCols = A.nCols();
   C.Resize( nRows, nCols );

   if( nRows*nCols > 0 ) {
      int row = 0;
      for (int i = 0; i < A.nRows(); ++i) {
         if (row_flag[i] != 0) {
            for (int j = 0; j < A.nCols(); ++j) {
               C(row, j) = A(i,j);
            }
            ++row;
         }
      }
   }
}


//=============================================================================
// scalar/Matrix arithmetic routines.
//=============================================================================

//-----------------------------------------------------------------------------
// scalar/Matrix addition:  C = a+A (term-by-term)
//-----------------------------------------------------------------------------
void Add_aM( double a, const Matrix& A, Matrix& C )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );

   // Do the update:  C = aA
   C = A;
   double* p = C.Base();

   for (int i = 0; i < C.nRows(); ++i)
      for (int j = 0; j < C.nCols(); ++j) {
         (*p) = a+(*p);
         ++p;
      }
}

//-----------------------------------------------------------------------------
// scalar/Matrix subtraction:  C = a-A (term-by-term)
//-----------------------------------------------------------------------------
void Subtract_aM( double a, const Matrix& A, Matrix& C )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );

   // Do the update:  C = aA
   C = A;
   double* p = C.Base();

   for (int i = 0; i < C.nRows(); ++i)
      for (int j = 0; j < C.nCols(); ++j) {
         (*p) = a-(*p);
         ++p;
      }
}

//-----------------------------------------------------------------------------
// scalar/Matrix multiplication:  C = a*A (term-by-term)
//-----------------------------------------------------------------------------
void Multiply_aM( double a, const Matrix& A, Matrix& C )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );

   // Do the update:  C = aA
   C = A;
   double* p = C.Base();

   for (int i = 0; i < C.nRows(); ++i)
      for (int j = 0; j < C.nCols(); ++j) {
         (*p) = a*(*p);
         ++p;
      }
}


//=============================================================================
// Matrix arithmetic routines.
//=============================================================================

//-----------------------------------------------------------------------------
// Matrix addition:  C = A + B
//-----------------------------------------------------------------------------
void Add_MM( const Matrix& A, const Matrix& B, Matrix& C )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );
   assert( B.nRows() > 0 && B.nCols() > 0 );
   assert( A.nRows() == B.nRows() && A.nCols() == B.nCols() );

   // Commensurate memory allocation.
   C.Resize( A.nRows(), A.nCols() );

   // Compute the Matrix addition:  C = A + B
   const double* p = A.Base();
   const double* q = B.Base();
   double*       r = C.Base();

   for (int i = 0; i < A.nRows(); ++i)
      for (int j = 0; j < A.nCols(); ++j)
         (*r++) = (*p++) + (*q++);
}

//-----------------------------------------------------------------------------
// Matrix subtraction:  C = A - B
//-----------------------------------------------------------------------------
void Subtract_MM( const Matrix& A, const Matrix& B, Matrix& C )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );
   assert( B.nRows() > 0 && B.nCols() > 0 );
   assert( A.nRows() == B.nRows() && A.nCols() == B.nCols() );

   // Commensurate memory allocation.
   C.Resize( A.nRows(), A.nCols() );

   // Compute the Matrix subtraction:  C = A - B
   const double* p = A.Base();
   const double* q = B.Base();
   double*       r = C.Base();

   for (int i = 0; i < A.nRows(); ++i)
      for (int j = 0; j < A.nCols(); ++j)
         (*r++) = (*p++) - (*q++);
}


//=============================================================================
// Matrix/Matrix multiplication routines.
//=============================================================================

//-----------------------------------------------------------------------------
// Matrix = Matrix/Matrix multiply:  C = AB
//-----------------------------------------------------------------------------
void Multiply_MM( const Matrix& A, const Matrix& B, Matrix& C )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );
   assert( B.nRows() > 0 && B.nCols() > 0 );
   assert( A.nCols() == B.nRows() );

   // Commensurate memory allocation.
   Matrix AB( A.nRows(), B.nCols() );

   // Compute the Matrix product.
   for (int i = 0; i < A.nRows(); ++i)
      for (int j = 0; j < B.nCols(); ++j)
         AB(i,j) = SumProduct( A.nCols(), A.Base(i,0), B.Base(0,j), B.nCols() );

   C = AB;
}

//-----------------------------------------------------------------------------
// Matrix = Matrix/Matrix multiply:  C = A'B
//-----------------------------------------------------------------------------
void Multiply_MtM( const Matrix& A, const Matrix& B, Matrix& C )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );
   assert( B.nRows() > 0 && B.nCols() > 0 );
   assert( A.nRows() == B.nRows() );

   // Commensurate memory allocation.
   Matrix AtB( A.nCols(), B.nCols() );

   // Compute the Matrix product.
   for (int i = 0; i < A.nCols(); ++i)
      for (int j = 0; j < B.nCols(); ++j)
         AtB(i,j) = SumProduct( A.nRows(), A.Base(0,i), A.nCols(), B.Base(0,j), B.nCols() );

   C = AtB;
}

//-----------------------------------------------------------------------------
// Matrix = Matrix/Matrix multiply:  C = AB'
//-----------------------------------------------------------------------------
void Multiply_MMt( const Matrix& A, const Matrix& B, Matrix& C )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );
   assert( B.nRows() > 0 && B.nCols() > 0 );
   assert( A.nCols() == B.nCols() );

   // Commensurate memory allocation.
   Matrix ABt( A.nRows(), B.nRows() );

   // Compute the Matrix product.
   for (int i = 0; i < A.nRows(); ++i)
      for (int j = 0; j < B.nRows(); ++j)
         ABt(i,j) = SumProduct( A.nCols(), A.Base(i,0), B.Base(j,0) );

   C = ABt;
}

//-----------------------------------------------------------------------------
// Matrix = Matrix/Matrix multiply:  C = A'B'
//-----------------------------------------------------------------------------
void Multiply_MtMt( const Matrix& A, const Matrix& B, Matrix& C )
{
   // Check the arguments.
   assert( A.nRows() > 0 && A.nCols() > 0 );
   assert( B.nRows() > 0 && B.nCols() > 0 );
   assert( A.nRows() == B.nCols() );

   // Commensurate memory allocation.
   Matrix AtBt( A.nCols(), B.nRows() );

   // Compute the Matrix product.
   for (int i = 0; i < A.nCols(); ++i)
      for (int j=0; j < B.nRows(); ++j)
         AtBt(i,j) = SumProduct( A.nRows(), A.Base(0,i), A.nCols(), B.Base(j,0) );

   C = AtBt;
}

//-----------------------------------------------------------------------------
// Dot product = A'B
//-----------------------------------------------------------------------------
double DotProduct( const Matrix& A, const Matrix& B )
{
   // Check the arguments.
   assert( isVector(A) && isVector(B) );
   assert( Length(A) == Length(B) );

   return SumProduct( Length(A), A.Base(), B.Base() );
}

//-----------------------------------------------------------------------------
// Quadratic form = a' B c
//-----------------------------------------------------------------------------
double QuadraticForm_MtMM( const Matrix& a, const Matrix& B, const Matrix& c )
{
   // Check the arguments.
   assert( a.nRows() > 0 && a.nCols() == 1 );
   assert( a.nRows() == B.nRows() );
   assert( c.nRows() > 0 && c.nCols() == 1 );
   assert( B.nCols() == c.nRows() );

   Matrix Bc;
   Multiply_MM(B,c,Bc);

   Matrix atBc;
   Multiply_MtM(a,Bc,atBc);

   return atBc(0,0);
}

//-----------------------------------------------------------------------------
// Quadratic form = a B c
//-----------------------------------------------------------------------------
double QuadraticForm_MMM( const Matrix& a, const Matrix& B, const Matrix& c )
{
   // Check the arguments.
   assert( a.nRows() == 1 && a.nCols() > 0 );
   assert( a.nCols() == B.nRows() );
   assert( c.nRows() > 0 && c.nCols() == 1 );
   assert( B.nCols() == c.nRows() );

   Matrix Bc;
   Multiply_MM(B,c,Bc);

   Matrix aBc;
   Multiply_MM(a,Bc,aBc);

   return aBc(0,0);
}

//=============================================================================
// Matrix comparison
//=============================================================================

//-----------------------------------------------------------------------------
bool isSquare(const Matrix& A) {
   if (A.nRows() > 0 && A.nRows() == A.nCols())
      return true;
   else
      return false;
}

//-----------------------------------------------------------------------------
bool isCongruent( const Matrix& A, const Matrix& B )
{
   // Compare the sizes
   if (A.nRows() == B.nRows() && A.nCols() == B.nCols())
      return true;
   else
      return false;
}

//-----------------------------------------------------------------------------
bool isClose( const Matrix& A, const Matrix& B, double tol )
{
   // Compare the sizes first.
   if (A.nRows() != B.nRows() || A.nCols() != B.nCols())
      return false;

   // Compare the contents.
   const double* p = A.Base();
   const double* q = B.Base();
   int n = A.nRows() * A.nCols();

   for (int k = 0; k<n; k++)
      if (abs((*p++) - (*q++)) > tol) return false;

   return true;
}

//=============================================================================
// isRow/Col/Vector
//=============================================================================

//-----------------------------------------------------------------------------
bool isRow( const Matrix& A )
{
   return A.nRows() == 1 && A.nCols() > 0;
}

//-----------------------------------------------------------------------------
bool isCol( const Matrix& A )
{
   return A.nCols() == 1 && A.nRows() > 0;
}

//-----------------------------------------------------------------------------
bool isVector( const Matrix& A )
{
   return isRow(A) || isCol(A);
}
