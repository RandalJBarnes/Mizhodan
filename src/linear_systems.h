//=============================================================================
// linear_systems.h
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#ifndef LINEAR_SYSTEMS_H
#define LINEAR_SYSTEMS_H

#include "matrix.h"


//=============================================================================
//
//=============================================================================
bool CholeskyDecomposition( const Matrix& A, Matrix& L );
void CholeskySolve( const Matrix& L, const Matrix& b, Matrix& x );
void CholeskyInverse( const Matrix& L, Matrix& Ainv );

bool RSPDInv( const Matrix& A, Matrix& Ainv );
bool LeastSquaresSolve( const Matrix& A, const Matrix& B, Matrix& X );

void AffineTransformation( const Matrix& A, const Matrix& B, const Matrix& C, Matrix& D );


//=============================================================================
#endif  // LINEAR_SYSTEMS_H
