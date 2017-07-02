//=============================================================================
// test_linear_systems.cpp
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    2 July 2017
//=============================================================================
#include <utility>

#include "test_linear_systems.h"
#include "unit_test.h"
#include "..\src\linear_systems.h"

//-----------------------------------------------------------------------------
// Hide all of the testing details inside an unnamed namespace. This allows me
// to create many small unit tests with polluting the global namespace.
//-----------------------------------------------------------------------------
namespace{
   const double TOLERANCE = 1e-9;

   //--------------------------------------------------------------------------
   // TestCholeskyDecomposition
   //--------------------------------------------------------------------------
   bool TestCholeskyDecomposition()
   {
      Matrix A("4,6,4,4; 6,10,9,7; 4,9,17,11; 4,7,11,18");
      Matrix L;
      CholeskyDecomposition(A,L);
      Matrix B("2,0,0,0; 3,1,0,0; 2,3,2,0; 2,1,2,3");

      return CHECK( isClose(L, B, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestCholeskySolve
   //--------------------------------------------------------------------------
   bool TestCholeskySolve()
   {
      Matrix A("4,6,4,4; 6,10,9,7; 4,9,17,11; 4,7,11,18");
      Matrix L;
      CholeskyDecomposition(A,L);
      Matrix B("44; 81; 117; 123");
      Matrix X;
      CholeskySolve(L,B,X);
      Matrix Z("1;2;3;4");

      return CHECK( isClose(X, Z, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestCholeskyInverse
   //--------------------------------------------------------------------------
   bool TestCholeskyInverse()
   {
      Matrix A("4,6,4,4; 6,10,9,7; 4,9,17,11; 4,7,11,18");
      Matrix B("945,-690,174,-48; -690,532,-140,32; 174,-140,52,-16; -48,32,-16,16");

      Matrix L;
      CholeskyDecomposition(A, L);

      Matrix Ainv;
      CholeskyInverse(L, Ainv);

      Matrix C;
      Multiply_aM(1.0/144.0, B, C);

      return CHECK( isClose(Ainv, C, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestRSPDInv
   //--------------------------------------------------------------------------
   bool TestRSPDInv()
   {
      Matrix A("4,6,4,4; 6,10,9,7; 4,9,17,11; 4,7,11,18");
      Matrix B("945,-690,174,-48; -690,532,-140,32; 174,-140,52,-16; -48,32,-16,16");
      Matrix Ainv;
      Multiply_aM(1.0/144.0, B, Ainv);
      RSPDInv(A,B);

      return CHECK( isClose(Ainv, B, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestLeastSquaresSolve
   //--------------------------------------------------------------------------
   bool TestLeastSquaresSolve()
   {
      Matrix A("5,2,8,1; 4,6,5,5; 7,1,1,3; 2,6,1,1; 4,6,7,4; 8,6,4,2; 5,8,7,1; 7,8,2,2; 6,7,5,2; 5,5,6,2");
      Matrix B("1,7,1; 6,7,2; 3,3,2; 5,2,5; 6,5,5; 4,6,1; 5,4,8; 4,2,6; 1,8,6; 4,1,1");
      Matrix X;
      LeastSquaresSolve(A,B,X);
      Matrix C("-0.122286918422277,0.266063484829536,-0.0575443373772838; 0.464217553042304,-0.0279214573318259,0.846505417553293; -0.00883317831785533,0.470311201138176,-0.027798955351842; 0.836316520297104,0.470195843209534,-0.259472798611811");

      return CHECK( isClose(X, C, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestAffineTransformation
   //--------------------------------------------------------------------------
   bool TestAffineTransformation()
   {
      Matrix A("7,8,6; 6,3,7; 6,1,6; 2,1,4; 1,8,8; 8,2,6; 5,5,6; 6,6,2");
      Matrix B("7,2,4; 5,1,2; 5,7,7");
      Matrix C("6,2,8");
      Matrix D;
      AffineTransformation(A,B,C,D);
      Matrix DD("125,66,94; 98,66,87; 83,57,76; 45,35,46; 93,68,84; 102,62,86; 96,59,80; 88,34,58");

      return CHECK( isClose(D, DD, TOLERANCE) );
   }

}

//-----------------------------------------------------------------------------
// test_LinearSystems
//-----------------------------------------------------------------------------
std::pair<int,int> test_LinearSystems()
{
   int nsucc = 0;
   int nfail = 0;

   TALLY( TestCholeskyDecomposition() );
   TALLY( TestCholeskySolve() );
   TALLY( TestCholeskyInverse() );
   TALLY( TestRSPDInv() );
   TALLY( TestLeastSquaresSolve() );
   TALLY( TestAffineTransformation() );

   return std::make_pair( nsucc, nfail );
}
