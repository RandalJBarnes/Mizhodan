//=============================================================================
// test_matrix.cpp
//
//    Test each of the class methods.  These tests are not exhaustive, but
//    every class method is exercised and the results compared to a known
//    correct result.
//
//    Test each of the functions associated with the Matrix class.  These
//    tests are not exhaustive, but every function is exercised and the
//    results compared to a known correct result.
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    2 July 2017
//=============================================================================
#include <iomanip>
#include <utility>

#include "test_matrix.h"
#include "unit_test.h"
#include "..\src\matrix.h"

//-----------------------------------------------------------------------------
// Hide all of the testing details inside an unnamed namespace. This allows me
// to create many small unit tests with polluting the global namespace.
//-----------------------------------------------------------------------------
namespace{
   const double TOLERANCE = 1e-9;

   //--------------------------------------------------------------------------
   // TestMatrixNullConstructor
   //--------------------------------------------------------------------------
   bool TestMatrixNullConstructor()
   {
      Matrix A;

      bool flag = true;

      flag &= CHECK( A.nRows() == 0 );
      flag &= CHECK( A.nCols() == 0 );

      return flag;
   }

   //--------------------------------------------------------------------------
   // TestMatrixCopyConstructor
   //--------------------------------------------------------------------------
   bool TestMatrixCopyConstructor()
   {
      Matrix A("1,2,3;4,5,6");
      Matrix B( A );

      return CHECK( isClose(A, B, TOLERANCE) );
   }


   //--------------------------------------------------------------------------
   // TestMatrixCopyConstructor
   //--------------------------------------------------------------------------
   bool TestMatrixConstructorFromVector()
   {
      Matrix A("1;2;3;4;5;6");
      std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
      Matrix B( v );

      return CHECK( isClose(A, B, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixDimensionedConstructor
   //--------------------------------------------------------------------------
   bool TestMatrixDimensionedConstructor()
   {
      Matrix A(2,3);
      Matrix B("0,0,0;0,0,0");

      return CHECK( isClose(A, B, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixConstructorWithScalarFill
   //--------------------------------------------------------------------------
   bool TestMatrixConstructorWithScalarFill()
   {
      Matrix A(2,3, 1.2);
      Matrix B("1.2,1.2,1.2;1.2,1.2,1.2");

      return CHECK( isClose(A, B, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixConstructorWithArrayFill
   //--------------------------------------------------------------------------
   bool TestMatrixConstructorWithArrayFill()
   {
      double A_data[] = {1.0,2.0,3.0,4.0,5.0,6.0};
      Matrix A(2,3, A_data);
      Matrix B("1,2,3;4,5,6");

      return CHECK( isClose(A, B, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixConstructorWithStringFill
   //--------------------------------------------------------------------------
   bool TestMatrixConstructorWithStringFill()
   {
      Matrix A("1,,;4,5,");
      Matrix B("1,0,0;4,5,0");

      return CHECK( isClose(A, B, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixDestructiveResize
   //--------------------------------------------------------------------------
   bool TestMatrixDestructiveResize()
   {
      Matrix A("1,2,3;4,5,6");
      A.Resize(2,2);
      Matrix B("0,0;0,0");

      return CHECK( isClose(A, B, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixAssignmentOperator
   //--------------------------------------------------------------------------
   bool TestMatrixAssignmentOperator()
   {
      Matrix A("1,2,3;4,5,6");
      Matrix B("0,1,1,0");
      B = A;

      return CHECK( isClose(A, B, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixScalarAssignment
   //--------------------------------------------------------------------------
   bool TestMatrixScalarAssignment()
   {
      Matrix A("1,2,3;4,5,6");
      A = 0.0;
      Matrix B("0,0,0;0,0,0");

      return CHECK( isClose(A, B, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixAccess
   //--------------------------------------------------------------------------
   bool TestMatrixAccess()
   {
      Matrix A(2,3);
      Matrix B("1,2,3;4,5,6");

      for (int i = 0; i < 2; ++i)
         for (int j = 0; j < 3; ++j)
            A(i,j) = B(i,j);

      return CHECK( isClose(A, B, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixRowAndColumnSize
   //--------------------------------------------------------------------------
   bool TestMatrixRowAndColumnSize()
   {
      Matrix A("1,2,3;4,5,6");

      bool flag = true;

      flag &= CHECK( A.nRows() == 2 );
      flag &= CHECK( A.nCols() == 3 );

      return flag;
   }

   //--------------------------------------------------------------------------
   // TestMatrixAccessToRawStorage
   //--------------------------------------------------------------------------
   bool TestMatrixAccessToRawStorage()
   {
      Matrix A("1,2,3;4,5,6");
      Matrix B(2,3);

      const double* pA = A.Base();
      double* pB = B.Base();

      for (int i = 0; i < 6; ++i)
         *pB++ = *pA++;

      return CHECK( isClose(A, B, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixAccessToRawStorageWithOffset
   //--------------------------------------------------------------------------
   bool TestMatrixAccessToRawStorageWithOffset()
   {
      Matrix A("1,2,3;4,5,6");
      Matrix B(2,3);

      for (int i = 0; i < 2; ++i)
         for (int j = 0; j < 3; ++j)
            *B.Base(i,j) = *A.Base(i,j);

      return CHECK( isClose(A, B, TOLERANCE) );
   }


   //--------------------------------------------------------------------------
   // TestMatrixColumnSum
   //--------------------------------------------------------------------------
   bool TestMatrixColumnSum()
   {
      Matrix A("1,2,3;4,5,6;7,8,9");
      Matrix x;
      ColumnSum( A, x );
      Matrix col_sum("12,15,18");

      return CHECK( isClose(x, col_sum, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixRowSum
   //--------------------------------------------------------------------------
   bool TestMatrixRowSum()
   {
      Matrix A("1,2,3;4,5,6;7,8,9");
      Matrix x;
      RowSum( A, x );
      Matrix row_sum("6;15;24");

      return CHECK( isClose(x, row_sum, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixLength
   //--------------------------------------------------------------------------
   bool TestMatrixLength()
   {
      Matrix A("1,2,3,4;5,6,7,8");
      Matrix B("1,2;3,4;5,6");
      Matrix C(0, 0);

      bool flag = true;
      flag &= CHECK( Length(A) == 4 );
      flag &= CHECK( Length(B) == 3 );
      flag &= CHECK( Length(C) == 0 );
      return flag;
   }

   //--------------------------------------------------------------------------
   // TestMatrixTrace
   //--------------------------------------------------------------------------
   bool TestMatrixTrace()
   {
      Matrix A("1,2,3;4,5,6;7,8,9");

      return CHECK( isClose(Trace(A), 15.0, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixSum
   //--------------------------------------------------------------------------
   bool TestMatrixSum()
   {
      Matrix A("1,2,3;4,5,6;7,8,9");

      return CHECK( isClose(Sum(A), 45.0, TOLERANCE) );
   }

//--------------------------------------------------------------------------
   // TestMatrixSumAbs
   //--------------------------------------------------------------------------
   bool TestMatrixSumAbs()
   {
      Matrix A("-1,2,-3;4,-5,6;-7,-8,9");

      return CHECK( isClose(SumAbs(A), 45.0, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixMaxAbs
   //--------------------------------------------------------------------------
   bool TestMatrixMaxAbs()
   {
      Matrix A("-1,2,-3;4,-5,6;-7,8,-9");

      return CHECK( isClose(MaxAbs(A), 9.0, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixL1Norm
   //--------------------------------------------------------------------------
   bool TestMatrixL1Norm()
   {
      Matrix A("1,2,3;4,5,6;7,8,9");

      return CHECK( isClose(L1Norm(A), 18.0, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixLInfNorm
   //--------------------------------------------------------------------------
   bool TestMatrixLInfNorm()
   {
      Matrix A("1,2,3;4,5,6;7,8,9");

      return CHECK( isClose(LInfNorm(A), 24.0, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixFNorm
   //--------------------------------------------------------------------------
   bool TestMatrixFNorm()
   {
      Matrix A("1,2,3;4,5,6;7,8,9");

      return CHECK( isClose(FNorm(A), 16.8819430161341, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixTranspose
   //--------------------------------------------------------------------------
   bool TestMatrixTranspose()
   {
      Matrix A("1,2,3;4,5,6;7,8,9");
      Matrix C;
      Transpose( A, C );
      Matrix At("1,4,7; 2,5,8; 3,6,9");

      return CHECK( isClose(C, At, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixNegative
   //--------------------------------------------------------------------------
   bool TestMatrixNegative()
   {
      Matrix A("1,2,3;4,5,6;7,8,9");
      Matrix C;
      Negative( A, C );
      Matrix B("-1,-2,-3;-4,-5,-6;-7,-8,-9");

      return CHECK( isClose(B, C, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixIdentity
   //--------------------------------------------------------------------------
   bool TestMatrixIdentity()
   {
      Matrix A(3,2,1.0);
      Identity( A, 4 );
      Matrix I("1,0,0,0; 0,1,0,0; 0,0,1,0; 0,0,0,1");

      return CHECK( isClose(A, I, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixSlice
   //--------------------------------------------------------------------------
   bool TestMatrixSlice()
   {
      Matrix A("1,2,3,4;5,6,7,8;9,10,11,12");
      Matrix B;

      std::vector<int> col_flag = { 1, 0, 1, 0 };
      std::vector<int> row_flag = { 1, 0, 1 };

      Slice( A, row_flag, col_flag, B );
      Matrix C("1,3;9,11");

      return CHECK( isClose(B, C, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixSliceRow
   //--------------------------------------------------------------------------
   bool TestMatrixSliceRow()
   {
      Matrix A("1,2,3,4;5,6,7,8;9,10,11,12");
      Matrix B;

      std::vector<int> row_flag = { 1, 0, 1 };

      SliceRows( A, row_flag, B );
      Matrix C("1,2,3,4;9,10,11,12");

      return CHECK( isClose(B, C, TOLERANCE) );
   }
   
   //--------------------------------------------------------------------------
   // TestMatrixAdd_aM
   //--------------------------------------------------------------------------
   bool TestMatrixAdd_aM()
   {
      Matrix A("1,2,3;4,5,6");
      Matrix B;
      Add_aM(2,A,B);
      Matrix C("3,4,5;6,7,8");

      return CHECK( isClose(B, C, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixSubtract_aM
   //--------------------------------------------------------------------------
   bool TestMatrixSubtract_aM()
   {
      Matrix A("1,2,3;4,5,6");
      Matrix B;
      Subtract_aM(2,A,B);
      Matrix C("1,0,-1;-2,-3,-4");

      return CHECK( isClose(B, C, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixMultiply_aM
   //--------------------------------------------------------------------------
   bool TestMatrixMultiply_aM()
   {
      Matrix A("1,2,3;4,5,6");
      Matrix B;
      Multiply_aM(2,A,B);
      Matrix Ax2("2,4,6;8,10,12");

      return CHECK( isClose(B, Ax2, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixAdd_MM
   //--------------------------------------------------------------------------
   bool TestMatrixAdd_MM()
   {
      Matrix A("1,2,3;4,5,6");
      Matrix B("1,0,1;0,0,1");
      Matrix C;
      Add_MM(A,B,C);
      Matrix ApB("2,2,4;4,5,7");

      return CHECK( isClose(C, ApB, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixSubtract_MM
   //--------------------------------------------------------------------------
   bool TestMatrixSubtract_MM()
   {
      Matrix A("1,2,3;4,5,6");
      Matrix B("1,0,1;0,0,1");
      Matrix C;
      Subtract_MM(A,B,C);
      Matrix AmB("0,2,2;4,5,5");

      return CHECK( isClose(C, AmB, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixMultiply_MM
   //--------------------------------------------------------------------------
   bool TestMatrixMultiply_MM()
   {
      Matrix A("1,2,3;4,5,6");
      Matrix B("1,2;3,4;5,6");
      Matrix C;
      Multiply_MM(A,B,C);
      Matrix AxB("22,28; 49,64");

      return CHECK( isClose(C, AxB, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixMultiply_MtM
   //--------------------------------------------------------------------------
   bool TestMatrixMultiply_MtM()
   {
      Matrix A("1,4;2,5;3,6");
      Matrix B("1,2;3,4;5,6");
      Matrix C;
      Multiply_MtM(A,B,C);
      Matrix AtxB("22,28; 49,64");

      return CHECK( isClose(C, AtxB, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixMultiply_MMt
   //--------------------------------------------------------------------------
   bool TestMatrixMultiply_MMt()
   {
      Matrix A("1,2,3;4,5,6");
      Matrix B("1,3,5;2,4,6");
      Matrix C;
      Multiply_MMt(A,B,C);
      Matrix AxBt("22,28; 49,64");

      return CHECK( isClose(C, AxBt, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixMultiply_MtMt
   //--------------------------------------------------------------------------
   bool TestMatrixMultiply_MtMt()
   {
      Matrix A("1,4;2,5;3,6");
      Matrix B("1,3,5;2,4,6");
      Matrix C;
      Multiply_MtMt(A,B,C);
      Matrix AtxBt("22,28; 49,64");

      return CHECK( isClose(C, AtxBt, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestDotProduct
   //--------------------------------------------------------------------------
   bool TestDotProduct()
   {
      Matrix A("1,2,3,4");
      Matrix B("1;2;3;4");
      Matrix C("4,3,2,1");
      Matrix D("4;3;2;1");

      bool flag = true;
      flag &= CHECK( isClose(DotProduct(A, A), 30, TOLERANCE) );
      flag &= CHECK( isClose(DotProduct(A, B), 30, TOLERANCE) );
      flag &= CHECK( isClose(DotProduct(A, C), 20, TOLERANCE) );
      flag &= CHECK( isClose(DotProduct(A, D), 20, TOLERANCE) );
      return flag;
   }

   //--------------------------------------------------------------------------
   // TestMatrixQuadraticForm_MtMM
   //--------------------------------------------------------------------------
   bool TestMatrixQuadraticForm_MtMM()
   {
      Matrix a("1;2;3");
      Matrix B("1,2,3;4,5,6;7,8,9");
      Matrix c("4;5;6");

      double q = QuadraticForm_MtMM(a,B,c);

      return CHECK( isClose(q, 552.0, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixQuadraticForm_MMM
   //--------------------------------------------------------------------------
   bool TestMatrixQuadraticForm_MMM()
   {
      Matrix a("1,2,3");
      Matrix B("1,2,3;4,5,6;7,8,9");
      Matrix c("4;5;6");

      double q = QuadraticForm_MMM(a,B,c);

      return CHECK( isClose(q, 552.0, TOLERANCE) );
   }

   //--------------------------------------------------------------------------
   // TestMatrixisSquare
   //--------------------------------------------------------------------------
   bool TestMatrixisSquare()
   {
      Matrix A("1,2,3;4,5,6;7,8,9");
      Matrix B("1,2,3;4,5,6");
      Matrix C;

      bool flag = true;
      flag &= CHECK( isSquare(A) );
      flag &= CHECK( !isSquare(B) );
      flag &= CHECK( !isSquare(C) );
      return flag;
   }


   //--------------------------------------------------------------------------
   // TestMatrixisCongruent
   //--------------------------------------------------------------------------
   bool TestMatrixisCongruent()
   {
      Matrix A("1,4;2,5;3,6");
      Matrix B("2,5;3,8;1,4");
      Matrix C("1,3,5;2,4,6");

      bool flag = true;
      flag &= CHECK( isCongruent(A, A) );
      flag &= CHECK( isCongruent(A, B) );
      flag &= CHECK( !isCongruent(A, C) );
      return flag;
   }

   //--------------------------------------------------------------------------
   // TestMatrixisClose
   //--------------------------------------------------------------------------
   bool TestMatrixisClose()
   {
      Matrix A("1,4;2,5;3,6");
      Matrix B("1,4;2,5;3,8");
      Matrix C("1,3,5;2,4,6");

      bool flag = true;

      flag &= CHECK( isClose(A, A, TOLERANCE) );
      flag &= CHECK( !isClose(A, B, 1.0) );
      flag &= CHECK( !isClose(A, C, TOLERANCE) );

      return flag;
   }

   //--------------------------------------------------------------------------
   // TestisRow
   //--------------------------------------------------------------------------
   bool TestMatrixisRow()
   {
      Matrix A("1,2,3,4");
      Matrix B("1;2;3;4");
      Matrix C("1,2;3,4");

      bool flag = true;

      flag &= CHECK( isRow(A) );
      flag &= CHECK( !isRow(B) );
      flag &= CHECK( !isRow(C) );
      return flag;
   }

   //--------------------------------------------------------------------------
   // TestisCol
   //--------------------------------------------------------------------------
   bool TestMatrixisCol()
   {
      Matrix A("1,2,3,4");
      Matrix B("1;2;3;4");
      Matrix C("1,2;3,4");

      bool flag = true;

      flag &= CHECK( !isCol(A) );
      flag &= CHECK( isCol(B) );
      flag &= CHECK( !isCol(C) );
      return flag;
   }

   //--------------------------------------------------------------------------
   // TestisVector
   //--------------------------------------------------------------------------
   bool TestMatrixisVector()
   {
      Matrix A("1,2,3,4");
      Matrix B("1;2;3;4");
      Matrix C("1,2;3,4");

      bool flag = true;

      flag &= CHECK( isVector(A) );
      flag &= CHECK( isVector(B) );
      flag &= CHECK( !isVector(C) );
      return flag;
   }
}

//-----------------------------------------------------------------------------
// test_Matrix
//-----------------------------------------------------------------------------
std::pair<int,int> test_Matrix()
{
   int nsucc = 0;
   int nfail = 0;

   TALLY( TestMatrixNullConstructor() );
   TALLY( TestMatrixCopyConstructor() );
   TALLY( TestMatrixConstructorFromVector() );
   TALLY( TestMatrixDimensionedConstructor() );
   TALLY( TestMatrixConstructorWithScalarFill() );
   TALLY( TestMatrixConstructorWithArrayFill() );
   TALLY( TestMatrixConstructorWithStringFill() );
   TALLY( TestMatrixDestructiveResize() );
   TALLY( TestMatrixAssignmentOperator() );
   TALLY( TestMatrixScalarAssignment() );
   TALLY( TestMatrixAccess() );
   TALLY( TestMatrixRowAndColumnSize() );
   TALLY( TestMatrixAccessToRawStorage() );
   TALLY( TestMatrixAccessToRawStorageWithOffset() );
   TALLY( TestMatrixColumnSum() );
   TALLY( TestMatrixRowSum() );
   TALLY( TestMatrixLength() );
   TALLY( TestMatrixTrace() );
   TALLY( TestMatrixSum() );
   TALLY( TestMatrixSumAbs() );
   TALLY( TestMatrixMaxAbs() );
   TALLY( TestMatrixL1Norm() );
   TALLY( TestMatrixLInfNorm() );
   TALLY( TestMatrixFNorm() );
   TALLY( TestMatrixTranspose() );
   TALLY( TestMatrixNegative() );
   TALLY( TestMatrixIdentity() );
   TALLY( TestMatrixSlice() );
   TALLY( TestMatrixSliceRow() );
   TALLY( TestMatrixAdd_aM() );
   TALLY( TestMatrixSubtract_aM() );
   TALLY( TestMatrixMultiply_aM() );
   TALLY( TestMatrixAdd_MM() );
   TALLY( TestMatrixSubtract_MM() );
   TALLY( TestMatrixMultiply_MM() );
   TALLY( TestMatrixMultiply_MtM() );
   TALLY( TestMatrixMultiply_MMt() );
   TALLY( TestMatrixMultiply_MtMt() );
   TALLY( TestDotProduct() );
   TALLY( TestMatrixQuadraticForm_MtMM() );
   TALLY( TestMatrixQuadraticForm_MMM() );
   TALLY( TestMatrixisSquare() );
   TALLY( TestMatrixisCongruent() );
   TALLY( TestMatrixisClose() );
   TALLY( TestMatrixisRow() );
   TALLY( TestMatrixisCol() );
   TALLY( TestMatrixisVector() );

   return std::make_pair( nsucc, nfail );
}
