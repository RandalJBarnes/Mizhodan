//=============================================================================
// unit_test.h
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil, Environmental, and Geo- Engineering
//    University of Minnesota
//
// version:
//    29 June 2017
//=============================================================================
#ifndef UNIT_TEST_H
#define UNIT_TEST_H


//=============================================================================
bool isClose( double x, double y, double tol );
bool Check( bool test, int line, const char* file );

#define CHECK(X) Check( (X), __LINE__, __FILE__ )
#define TALLY(X) ( (X) ? ++nsucc : ++nfail );

//=============================================================================
#endif  // UNIT_TEST_H
