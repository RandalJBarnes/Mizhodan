//=============================================================================
// linear_systems.cpp
//
//    A minimal set of decomposition and solution routines for systems of
//    linear equations.
//
// references:
// o  Golub, G.H., and Van Loan, C.F., 1996, MATRIX COMPUTATIONS, 3rd Edition,
//    Johns Hopkins University Press, Baltimore, Maryland, 694 pp.
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#include "linear_systems.h"

#include <cassert>
#include <cmath>

#include "sum_product-inl.h"

namespace{
   double MIN_DIVISOR = 1e-12;
}

//=============================================================================
// CholeskyDecomposition
//
//    Compute the Cholesky decomposition of the symmetric positive definite
//    Matrix "A".
//
// Arguments:
//
//    A     on entrance, a symmetric positive definite Matrix.
//
//    L     on exit, the lower triangular Matrix L where A = LL'.
//
// Return:
//
//    true  if the decomposition was completed successfully;
//    false if not.
//
// Notes:
//
// o  This routine is based upon Golub and Van Loan, 1996, Algorithm 4.2-1,
//    page 144.
//
// o  Only the lower triangular portion of A is accessed, so only the lower
//    triangular portion needs to be filled.
//
// o  The routine CholeskySolve is this routine's complementary pair.
//
// References:
//
// o  Golub, G.H., and Van Loan, C.F., 1996, MATRIX COMPUTATIONS, 3rd Edition,
//    Johns Hopkins University Press, Baltimore, Maryland, 694 pp.
//=============================================================================
bool CholeskyDecomposition( const Matrix& A, Matrix& L )
{
   // Validate the arguments.
   assert(isSquare(A));

   // Define local constants.
   const int N = A.nRows();

   // Carry out the Cholesky decomposition on Matrix "A".
   L = A;
   for (int j = 0; j < N; ++j) {
      if (j > 0) {
         for (int k = j; k < N; ++k)
            L(k,j) -= SumProduct(j, L.Base(j,0), L.Base(k,0));
      }

      if (L(j,j) < MIN_DIVISOR) return false;
      L(j,j) = sqrt(L(j,j));

      for (int k = j+1; k < N; ++k) {
         L(k,j) /= L(j,j);
         L(j,k) = 0.0;
      }
   }
   return true;
}

//=============================================================================
// CholeskySolve
//
//    This routine solves the system of linear equations given by "LL' x = b",
//    using the Cholesky factorizion of matrix "LL' = A" and forward
//    elimination followed by back substitution.
//
// Arguments:
//
//    L     the Cholesky decomposition of a symmetric positive definite
//          matrix A = LL'.
//    b     on entrance, the right hand side of the system of equations;
//          on exit, the solution.
//
// Notes:
//
// o  This routine is based upon Golub and Van Loan, 1983, Algorithm 4.1-1,
//    page 53.
//
// o  The Cholesky decomposition MUST be successfully carried out before
//    calling this routine.
//
// o  This routine works in the following manner
//
//       L y = b  then  L'x = y
//
//    however, in both sub-systems the soultion overwrites the right hand
//    side vector b.
//
// References:
//
//    Golub, G.H., and Van Loan, C.F., 1983, MATRIX COMPUTATIONS, Johns
//    Hopkins University Press, Baltimore, Maryland, 476 pp.
//
//=============================================================================
void CholeskySolve( const Matrix& L, const Matrix& b, Matrix& x )
{
   // Validate the arguments.
   assert( L.nRows() == L.nCols() );
   assert( b.nRows() == L.nRows() );

   // Define local constants.
   const int N = L.nRows();

   // Solve L y = b using forward elimination.
   x = b;

   double Sum;
   for (int i = 0; i < N; i++) {
      Sum = x(i,0);
      for (int j = 0; j < i; ++j)
         Sum -= L(i,j) * x(j,0);

      x(i,0) = Sum / L(i,i);
   }

   // Solve L' x = y using back substitution.
   // See Golub and Van Loan, 1983, Algorithm 4.1-2, page 53.
   for (int i = N-1; i >= 0; --i) {
      Sum = x(i,0);
      for (int j=i+1; j<N; ++j)
         Sum -= L(j,i) * x(j,0);

      x(i,0) = Sum / L(i,i);
   }
}

//=============================================================================
// CholeskyInverse
//
//    Return the inverse of a real, symmetric, positive definite Matrix A
//    whose Cholesky decomposition is given by L.
//
// Arguments:
//    L     on entrance, the Cholesky decomposition of a real, symmetric,
//          positive definite Matrix.
//
//    Ainv  on exit, the inverse of A.
//
// Notes:
// o  The computation of the inverse is based upon the standard Cholesky
//    decompostion.
//
// References:
// o  Stewart, G., 1998, "Matrix Algorithms - Volume I: Basic Decompositions",
//    SIAM, Philadelphia, 458pp., ISBN 0-89871-414-1.
//=============================================================================
void CholeskyInverse( const Matrix& L, Matrix& Ainv ) 
{
   assert( L.nRows() > 0 );
   assert( L.nRows() == L.nCols() );
   const int N = L.nRows();

   Matrix LL(L);

   // Invert L in place; remember that L is lower triangular.
   for (int k = 0; k < N; ++k) {
      LL(k,k) = 1.0/LL(k,k);

      for (int i = 0; i < k; ++i)
         LL(k,i) = -LL(k,k) * SumProduct( k-i, LL.Base(i,i), N, LL.Base(k,i) );
   }

   // A = L L' --> Ainv = (L')~ L~ = (L~)' L~
   Multiply_MtM( LL, LL, Ainv );
}

//=============================================================================
// RSPDInv
//
//    Return the inverse of a real, symmetric, positive definite Matrix A.
//
// Arguments:
//    A     on entrance, a real, symmetric, positive definite Matrix.
//
//    Ainv  on exit, the inverse of A.
//
// Notes:
// o  The computation of the inverse is based upon the standard Cholesky
//    decompostion.
//
// o  Only the lower triangular portion of A is accessed, so only the lower
//    triangular portion needs to be filled.
//
// o  The inversion in place for an lower triangular Matrix is based upon
//    the pseudo-code given in Stewart (1998, p. 179).
//
// o  The matrices A and Ainv may be the same space in memory.
//
// References:
// o  Stewart, G., 1998, "Matrix Algorithms - Volume I: Basic Decompositions",
//    SIAM, Philadelphia, 458pp., ISBN 0-89871-414-1.
//=============================================================================
bool RSPDInv( const Matrix& A, Matrix& Ainv )
{
   assert(isSquare(A));
   const int N = A.nRows();

   // Compute the Cholesky decomposition of "A", putting the result in "L".
   Matrix L;
   CholeskyDecomposition(A,L);

   // Invert L in place; remember that L is lower triangular.
   for (int k = 0; k < N; ++k) {
      L(k,k) = 1.0/L(k,k);

      for (int i = 0; i < k; ++i)
         L(k,i) = -L(k,k) * SumProduct( k-i, L.Base(i,i), N, L.Base(k,i) );
   }

   // A = L L' --> Ainv = (L')~ L~ = (L~)' L~
   Multiply_MtM( L, L, Ainv );

   return true;
}


//=============================================================================
// LeastSquaresSolve
//
// Purpose:
//    This routine computes the least-squares solution to the overdetermined
//    system of linear equations in a Matrix-based form:
//
//       A X = B
//
//    using a modified Gram-Schmidt orthogonalization algorithm.
//
// Arguments:
//    A        (m x n) coefficient Matrix.
//    B        (m x p) right-hand-side column Matrix.
//    X        (n x p) solution column Matrix.
//
// Return:
//    true     if the solution is computed, and false otherwise.
//
// Notes:
// o  Matrix A must have at least as many rows as columns (i.e. m >= n), and
//    it must have full column rank:  i.e. rank(A) = n.
//
// o  If Matrix A is rank deficient, rank(A) < n, then false is returned.
//
// o  This routine uses a modified Gram-Schmidt algorithm.  See Golub and
//    Van Loan (1996, Algorithm 5.2.5).  The Matrix A is rewritten using a
//    QR factorization:
//
//       A = Q R
//
//    where R (n x n) is upper triangular, and Q (m x n) with orthogonal
//    columns
//
//       Q'Q = I
//
//    Thus, the solution to the least squares problem can be given by
//    computing a simple back-substitution solution to
//
//       RX = Q'B
//
// o  However, following the recommendation of Golub and Van Loan (1996),
//    Section 5.3.5, the error properties of the solution can be improved
//    significantly by computing the factorization of the augmented Matrix
//
//       [A,B]  =  [Q,S] [R,Z]
//                       [0,P]
//
//              =  [QR, QZ+SP]
//
//    where the augmented [Q,S] has orthogonal columns:
//
//       [Q,S]' [Q,S] = [Q'Q   Q'S] = [I,0]
//                      [S'Q   S'S]   [0,I]
//
//    Thus,
//
//       B  =  QZ + SP
//
//    and
//
//       Q'B  =  Q'QZ + Q'SP
//            =  IZ + 0P
//            =  Z
//
//    We solve for X using back-substitute on Z:
//
//       R X  =  Z
//
// References:
// o  Golub, G. H., and C. F. Van Loan, 1996, MATRIX COMPUTATIONS (3rd
//    Edition), Johns Hopkins University Press, Baltimore Maryland,
//    ISBN 0-8018-5414-8.
//=============================================================================
bool LeastSquaresSolve( const Matrix& A, const Matrix& B, Matrix& X )
{
   assert(A.nRows() == B.nRows());

   // Setup the necessary dimension constants.
   const int M = A.nRows();
   const int N = A.nCols();
   const int P = B.nCols();

   // Allocate the space for the solution.
   X.Resize(N,P);

   // By design, this algorithm operates on A and B in place.  That is, it is
   // a destructive routine.  To eliminate side-effects, we work work on copies
   // of A and B.  This is a bit slower, but much safer.
   Matrix AA( A );
   Matrix BB( B );

   // Allocate the necessary local memory for the upper-triangular Matrix R.
   // The augmenting column "z" will be stored in "X".
   Matrix R(N,N);

   // Carry out the modified Gram-Schmidt orthogonalization:  i.e. Golub and
   // Van Loan (1996), Algorithm 5.2.5. applied to the augmented coefficient
   // Matrix.
   for (int k = 0; k < N; ++k) {
      double Sum = SumProduct(M, AA.Base(0,k), AA.nCols());
      if (Sum < MIN_DIVISOR) return false;

      R(k,k) = sqrt(Sum);
      for (int i = 0; i < M; ++i)
         AA(i,k) /= R(k,k);

      for (int j = k+1; j < N; ++j) {
         R(k,j) = SumProduct(M, AA.Base(0,k), AA.nCols(), AA.Base(0,j), AA.nCols());
         for (int i = 0; i < M; ++i)
            AA(i,j) -= AA(i,k)*R(k,j);
      }

      for (int p = 0; p < P; ++p) {
         X(k,p) = SumProduct(M, AA.Base(0,k), AA.nCols(), BB.Base(0,p), BB.nCols());
         for (int i = 0; i < M; ++i)
            BB(i,p) -= AA(i,k)*X(k,p);
      }
   }

   // Compute X = R~Z using back-substitution: e.g. Golub and Van Loan (1996)
   // Algorithm 3.1.2. Recall that the augmenting Matrix "z" is stored in "X".
   if( abs(R(N-1,N-1)) < MIN_DIVISOR ) return false;

   for (int p = 0; p < P; ++p)
      X(N-1,p) /= R(N-1,N-1);

   for (int i = N-2; i >= 0; --i) {
      if (abs(R(i,i)) < MIN_DIVISOR ) return false;

      for (int p = 0; p < P; ++p)
         X(i,p) = (X(i,p) - SumProduct(N-i-1, R.Base(i,i+1), X.Base(i+1,p), X.nCols())) / R(i,i);
   }
   return true;
}


//=============================================================================
// AffineTransformation
//
// Purpose:
//    Compute the affine transformation of each row of A, that is
//
//       D(i,:) = A(i,:)*B + C
//
//    This is a row-by-row operation.  This transformation preserves the
//    dimension of the row vectors.
//
// Arguments:
//    A     (MxN) Matrix containing the rows to be transformed.
//    B     (NxN) rotation and scale Matrix.
//    C     (1xN) shift Matrix.
//    D     (MxN) Matrix contraining the transofrmed rows (on exit).
//=============================================================================
void AffineTransformation( const Matrix& A, const Matrix& B, const Matrix& C, Matrix& D )
{
   assert( A.nCols() == B.nRows() );
   assert( B.nRows() == B.nCols() );
   assert( C.nRows() == 1 );
   assert( C.nCols() == B.nCols() );

   Matrix DD;
   Multiply_MM(A,B,DD);

   for (int i = 0; i < DD.nRows(); ++i) {
      for (int j = 0; j < DD.nCols(); ++j) {
         DD(i,j) += C(0,j);
      }
   }
   D = DD;
}
